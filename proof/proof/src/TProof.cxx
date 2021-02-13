// @(#)root/proof:$Id: a2a50e759072c37ccbc65ecbcce735a76de86e95 $
// Author: Fons Rademakers   13/02/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
/**
  \defgroup proof PROOF

  Classes defining the Parallel ROOT Facility, PROOF, a framework for parallel analysis of ROOT TTrees.

*/

/**
  \defgroup proofkernel PROOF kernel Libraries
  \ingroup proof

  The PROOF kernel libraries (libProof, libProofPlayer, libProofDraw) contain the classes defining
  the kernel of the PROOF facility, i.e. the protocol and the utilities to steer data processing
  and handling of results.

*/

/** \class TProof
\ingroup proofkernel

This class controls a Parallel ROOT Facility, PROOF, cluster.
It fires the worker servers, it keeps track of how many workers are
running, it keeps track of the workers running status, it broadcasts
messages to all workers, it collects results, etc.

*/

#include <stdlib.h>
#include <fcntl.h>
#include <errno.h>
#ifdef WIN32
#   include <io.h>
#   include <sys/stat.h>
#   include <sys/types.h>
#   include "snprintf.h"
#else
#   include <unistd.h>
#endif
#include <vector>

#include "RConfigure.h"
#include "Riostream.h"
#include "Getline.h"
#include "TBrowser.h"
#include "TChain.h"
#include "TCondor.h"
#include "TDSet.h"
#include "TError.h"
#include "TEnv.h"
#include "TEntryList.h"
#include "TEventList.h"
#include "TFile.h"
#include "TFileInfo.h"
#include "TFTP.h"
#include "THashList.h"
#include "TInterpreter.h"
#include "TKey.h"
#include "TMap.h"
#include "TMath.h"
#include "TMessage.h"
#include "TMonitor.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TProof.h"
#include "TProofNodeInfo.h"
#include "TProofOutputFile.h"
#include "TVirtualProofPlayer.h"
#include "TVirtualPacketizer.h"
#include "TProofServ.h"
#include "TPluginManager.h"
#include "TQueryResult.h"
#include "TRandom.h"
#include "TRegexp.h"
#include "TROOT.h"
#include "TSlave.h"
#include "TSocket.h"
#include "TSortedList.h"
#include "TSystem.h"
#include "TTree.h"
#include "TUrl.h"
#include "TFileCollection.h"
#include "TDataSetManager.h"
#include "TDataSetManagerFile.h"
#include "TMacro.h"
#include "TSelector.h"
#include "TPRegexp.h"
#include "TPackMgr.h"

#include <mutex>

TProof *gProof = 0;

// Rotating indicator
char TProofMergePrg::fgCr[4] = {'-', '\\', '|', '/'};

TList   *TProof::fgProofEnvList = 0;          // List of env vars for proofserv
TPluginHandler *TProof::fgLogViewer = 0;      // Log viewer handler

ClassImp(TProof);

//----- PROOF Interrupt signal handler -----------------------------------------
////////////////////////////////////////////////////////////////////////////////
/// TProof interrupt handler.

Bool_t TProofInterruptHandler::Notify()
{
   if (!fProof->IsTty() || fProof->GetRemoteProtocol() < 22) {

      // Cannot ask the user : abort any remote processing
      fProof->StopProcess(kTRUE);

   } else {
      // Real stop or request to switch to asynchronous?
      const char *a = 0;
      if (fProof->GetRemoteProtocol() < 22) {
         a = Getline("\nSwitch to asynchronous mode not supported remotely:"
                     "\nEnter S/s to stop, Q/q to quit, any other key to continue: ");
      } else {
         a = Getline("\nEnter A/a to switch asynchronous, S/s to stop, Q/q to quit,"
                     " any other key to continue: ");
      }
      if (a[0] == 'Q' || a[0] == 'S' || a[0] == 'q' || a[0] == 's') {

         Info("Notify","Processing interrupt signal ... %c", a[0]);

         // Stop or abort any remote processing
         Bool_t abort = (a[0] == 'Q' || a[0] == 'q') ? kTRUE : kFALSE;
         fProof->StopProcess(abort);

      } else if ((a[0] == 'A' || a[0] == 'a') && fProof->GetRemoteProtocol() >= 22) {
         // Stop any remote processing
         fProof->GoAsynchronous();
      }
   }

   return kTRUE;
}

//----- Input handler for messages from TProofServ -----------------------------
////////////////////////////////////////////////////////////////////////////////
/// Constructor

TProofInputHandler::TProofInputHandler(TProof *p, TSocket *s)
                   : TFileHandler(s->GetDescriptor(),1),
                     fSocket(s), fProof(p)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Handle input

Bool_t TProofInputHandler::Notify()
{
   fProof->CollectInputFrom(fSocket);
   return kTRUE;
}


//------------------------------------------------------------------------------

ClassImp(TSlaveInfo);

////////////////////////////////////////////////////////////////////////////////
/// Used to sort slaveinfos by ordinal.

Int_t TSlaveInfo::Compare(const TObject *obj) const
{
   if (!obj) return 1;

   const TSlaveInfo *si = dynamic_cast<const TSlaveInfo*>(obj);

   if (!si) return fOrdinal.CompareTo(obj->GetName());

   const char *myord = GetOrdinal();
   const char *otherord = si->GetOrdinal();
   while (myord && otherord) {
      Int_t myval = atoi(myord);
      Int_t otherval = atoi(otherord);
      if (myval < otherval) return 1;
      if (myval > otherval) return -1;
      myord = strchr(myord, '.');
      if (myord) myord++;
      otherord = strchr(otherord, '.');
      if (otherord) otherord++;
   }
   if (myord) return -1;
   if (otherord) return 1;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Used to compare slaveinfos by ordinal.

Bool_t TSlaveInfo::IsEqual(const TObject* obj) const
{
   if (!obj) return kFALSE;
   const TSlaveInfo *si = dynamic_cast<const TSlaveInfo*>(obj);
   if (!si) return kFALSE;
   return (strcmp(GetOrdinal(), si->GetOrdinal()) == 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Print slave info. If opt = "active" print only the active
/// slaves, if opt="notactive" print only the not active slaves,
/// if opt = "bad" print only the bad slaves, else
/// print all slaves.

void TSlaveInfo::Print(Option_t *opt) const
{
   TString stat = fStatus == kActive ? "active" :
                  fStatus == kBad ? "bad" :
                  "not active";

   Bool_t newfmt = kFALSE;
   TString oo(opt);
   if (oo.Contains("N")) {
      newfmt = kTRUE;
      oo.ReplaceAll("N","");
   }
   if (oo == "active" && fStatus != kActive) return;
   if (oo == "notactive" && fStatus != kNotActive) return;
   if (oo == "bad" && fStatus != kBad) return;

   if (newfmt) {
      TString msd, si, datadir;
      if (!(fMsd.IsNull())) msd.Form("| msd: %s ", fMsd.Data());
      if (!(fDataDir.IsNull())) datadir.Form("| datadir: %s ", fDataDir.Data());
      if (fSysInfo.fCpus > 0) {
         si.Form("| %s, %d cores, %d MB ram", fHostName.Data(),
               fSysInfo.fCpus, fSysInfo.fPhysRam);
      } else {
         si.Form("| %s", fHostName.Data());
      }
      Printf("Worker: %9s %s %s%s| %s", fOrdinal.Data(), si.Data(), msd.Data(), datadir.Data(), stat.Data());

   } else {
      TString msd  = fMsd.IsNull() ? "<null>" : fMsd.Data();

      std::cout << "Slave: "          << fOrdinal
         << "  hostname: "     << fHostName
         << "  msd: "          << msd
         << "  perf index: "   << fPerfIndex
         << "  "               << stat
         << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Setter for fSysInfo

void TSlaveInfo::SetSysInfo(SysInfo_t si)
{
   fSysInfo.fOS       = si.fOS;          // OS
   fSysInfo.fModel    = si.fModel;       // computer model
   fSysInfo.fCpuType  = si.fCpuType;     // type of cpu
   fSysInfo.fCpus     = si.fCpus;        // number of cpus
   fSysInfo.fCpuSpeed = si.fCpuSpeed;    // cpu speed in MHz
   fSysInfo.fBusSpeed = si.fBusSpeed;    // bus speed in MHz
   fSysInfo.fL2Cache  = si.fL2Cache;     // level 2 cache size in KB
   fSysInfo.fPhysRam  = si.fPhysRam;     // Physical RAM
}

ClassImp(TProof);

//------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TMergerInfo::~TMergerInfo()
{
   // Just delete the list, the objects are owned by other list
   if (fWorkers) {
      fWorkers->SetOwner(kFALSE);
      SafeDelete(fWorkers);
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Increase number of already merged workers by 1

void TMergerInfo::SetMergedWorker()
{
   if (AreAllWorkersMerged())
      Error("SetMergedWorker", "all workers have been already merged before!");
   else
      fMergedWorkers++;
}

////////////////////////////////////////////////////////////////////////////////
/// Add new worker to the list of workers to be merged by this merger

void TMergerInfo::AddWorker(TSlave *sl)
{
   if (!fWorkers)
      fWorkers = new TList();
   if (fWorkersToMerge == fWorkers->GetSize()) {
      Error("AddWorker", "all workers have been already assigned to this merger");
      return;
   }
   fWorkers->Add(sl);
}

////////////////////////////////////////////////////////////////////////////////
/// Return if merger has already merged all workers, i.e. if it has finished its merging job

Bool_t TMergerInfo::AreAllWorkersMerged()
{
   return (fWorkersToMerge == fMergedWorkers);
}

////////////////////////////////////////////////////////////////////////////////
/// Return if the determined number of workers has been already assigned to this merger

Bool_t TMergerInfo::AreAllWorkersAssigned()
{
      if (!fWorkers)
         return kFALSE;

      return (fWorkers->GetSize() == fWorkersToMerge);
}

////////////////////////////////////////////////////////////////////////////////
/// This a private API function.
/// It checks whether the connection string contains a PoD cluster protocol.
/// If it does, then the connection string will be changed to reflect
/// a real PROOF connection string for a PROOF cluster managed by PoD.
/// PoD: http://pod.gsi.de .
/// Return -1 if the PoD request failed; return 0 otherwise.

static Int_t PoDCheckUrl(TString *_cluster)
{
   if ( !_cluster )
      return 0;

   // trim spaces from both sides of the string
   *_cluster = _cluster->Strip( TString::kBoth );
   // PoD protocol string
   const TString pod_prot("pod");

   // URL test
   // TODO: The URL test is to support remote PoD servers (not managed by pod-remote)
   TUrl url( _cluster->Data() );
   if( pod_prot.CompareTo(url.GetProtocol(), TString::kIgnoreCase) )
      return 0;

   // PoD cluster is used
   // call pod-info in a batch mode (-b).
   // pod-info will find either a local PoD cluster or
   // a remote one, manged by pod-remote.
   *_cluster = gSystem->GetFromPipe("pod-info -c -b");
   if( 0 == _cluster->Length() ) {
      Error("PoDCheckUrl", "PoD server is not running");
      return -1;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a PROOF environment. Starting PROOF involves either connecting
/// to a master server, which in turn will start a set of slave servers, or
/// directly starting as master server (if master = ""). Masterurl is of
/// the form: [proof[s]://]host[:port]. Conffile is the name of the config
/// file describing the remote PROOF cluster (this argument alows you to
/// describe different cluster configurations).
/// The default is proof.conf. Confdir is the directory where the config
/// file and other PROOF related files are (like motd and noproof files).
/// Loglevel is the log level (default = 1). User specified custom config
/// files will be first looked for in $HOME/.conffile.

TProof::TProof(const char *masterurl, const char *conffile, const char *confdir,
               Int_t loglevel, const char *alias, TProofMgr *mgr)
       : fUrl(masterurl)
{
   // Default initializations
   InitMembers();

   // This may be needed during init
   fManager = mgr;

   // Default server type
   fServType = TProofMgr::kXProofd;

   // Default query mode
   fQueryMode = kSync;

   // Parse the main URL, adjusting the missing fields and setting the relevant
   // bits
   ResetBit(TProof::kIsClient);
   ResetBit(TProof::kIsMaster);

   // Protocol and Host
   if (!masterurl || strlen(masterurl) <= 0) {
      fUrl.SetProtocol("proof");
      fUrl.SetHost("__master__");
   } else if (!(strstr(masterurl, "://"))) {
      fUrl.SetProtocol("proof");
   }
   // Port
   if (fUrl.GetPort() == TUrl(" ").GetPort())
      fUrl.SetPort(TUrl("proof:// ").GetPort());

   // Make sure to store the FQDN, so to get a solid reference for subsequent checks
   if (!strcmp(fUrl.GetHost(), "__master__"))
      fMaster = fUrl.GetHost();
   else if (!strlen(fUrl.GetHost()))
      fMaster = gSystem->GetHostByName(gSystem->HostName()).GetHostName();
   else
      fMaster = gSystem->GetHostByName(fUrl.GetHost()).GetHostName();

   // Server type
   if (strlen(fUrl.GetOptions()) > 0) {
      TString opts(fUrl.GetOptions());
      if (!(strncmp(fUrl.GetOptions(),"std",3))) {
         fServType = TProofMgr::kProofd;
         opts.Remove(0,3);
         fUrl.SetOptions(opts.Data());
      } else if (!(strncmp(fUrl.GetOptions(),"lite",4))) {
         fServType = TProofMgr::kProofLite;
         opts.Remove(0,4);
         fUrl.SetOptions(opts.Data());
      }
   }

   // Instance type
   fMasterServ = kFALSE;
   SetBit(TProof::kIsClient);
   ResetBit(TProof::kIsMaster);
   if (fMaster == "__master__") {
      fMasterServ = kTRUE;
      ResetBit(TProof::kIsClient);
      SetBit(TProof::kIsMaster);
   } else if (fMaster == "prooflite") {
      // Client and master are merged
      fMasterServ = kTRUE;
      SetBit(TProof::kIsMaster);
   }
   // Flag that we are a client
   if (TestBit(TProof::kIsClient))
      if (!gSystem->Getenv("ROOTPROOFCLIENT")) gSystem->Setenv("ROOTPROOFCLIENT","");

   Init(masterurl, conffile, confdir, loglevel, alias);

   // If the user was not set, get it from the master
   if (strlen(fUrl.GetUser()) <= 0) {
      TString usr, emsg;
      if (Exec("gProofServ->GetUser()", "0", kTRUE) == 0) {
         TObjString *os = fMacroLog.GetLineWith("const char");
         if (os) {
            Ssiz_t fst =  os->GetString().First('\"');
            Ssiz_t lst =  os->GetString().Last('\"');
            usr = os->GetString()(fst+1, lst-fst-1);
         } else {
            emsg = "could not find 'const char *' string in macro log";
         }
      } else {
         emsg = "could not retrieve user info";
      }
      if (!emsg.IsNull()) {
         // Get user logon name
         UserGroup_t *pw = gSystem->GetUserInfo();
         if (pw) {
            usr = pw->fUser;
            delete pw;
         }
         Warning("TProof", "%s: using local default %s", emsg.Data(), usr.Data());
      }
      // Set the user name in the main URL
      fUrl.SetUser(usr.Data());
   }

   // If called by a manager, make sure it stays in last position
   // for cleaning
   if (mgr) {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfSockets()->Remove(mgr);
      gROOT->GetListOfSockets()->Add(mgr);
   }

   // Old-style server type: we add this to the list and set the global pointer
   if (IsProofd() || TestBit(TProof::kIsMaster))
      if (!gROOT->GetListOfProofs()->FindObject(this))
         gROOT->GetListOfProofs()->Add(this);

   // Still needed by the packetizers: needs to be changed
   gProof = this;
}

////////////////////////////////////////////////////////////////////////////////
/// Protected constructor to be used by classes deriving from TProof
/// (they have to call Init themselves and override StartSlaves
/// appropriately).
///
/// This constructor simply closes any previous gProof and sets gProof
/// to this instance.

TProof::TProof() : fUrl(""), fServType(TProofMgr::kXProofd)
{
   // Default initializations
   InitMembers();

   if (!gROOT->GetListOfProofs()->FindObject(this))
      gROOT->GetListOfProofs()->Add(this);

   gProof = this;
}

////////////////////////////////////////////////////////////////////////////////
/// Default initializations

void TProof::InitMembers()
{
   fValid = kFALSE;
   fTty = kFALSE;
   fRecvMessages = 0;
   fSlaveInfo = 0;
   fMasterServ = kFALSE;
   fSendGroupView = kFALSE;
   fIsPollingWorkers = kFALSE;
   fLastPollWorkers_s = -1;
   fActiveSlaves = 0;
   fInactiveSlaves = 0;
   fUniqueSlaves = 0;
   fAllUniqueSlaves = 0;
   fNonUniqueMasters = 0;
   fActiveMonitor = 0;
   fUniqueMonitor = 0;
   fAllUniqueMonitor = 0;
   fCurrentMonitor = 0;
   fBytesRead = 0;
   fRealTime = 0;
   fCpuTime = 0;
   fIntHandler = 0;
   fProgressDialog = 0;
   fProgressDialogStarted = kFALSE;
   SetBit(kUseProgressDialog);
   fPlayer = 0;
   fFeedback = 0;
   fChains = 0;
   fDSet = 0;
   fNotIdle = 0;
   fSync = kTRUE;
   fRunStatus = kRunning;
   fIsWaiting = kFALSE;
   fRedirLog = kFALSE;
   fLogFileW = 0;
   fLogFileR = 0;
   fLogToWindowOnly = kFALSE;
   fSaveLogToMacro = kFALSE;
   fMacroLog.SetName("ProofLogMacro");

   fWaitingSlaves = 0;
   fQueries = 0;
   fOtherQueries = 0;
   fDrawQueries = 0;
   fMaxDrawQueries = 1;
   fSeqNum = 0;

   fSessionID = -1;
   fEndMaster = kFALSE;

   fPackMgr = 0;
   fEnabledPackagesOnCluster = 0;

   fInputData = 0;

   fPrintProgress = 0;

   fLoadedMacros = 0;

   fProtocol = -1;
   fSlaves = 0;
   fTerminatedSlaveInfos = 0;
   fBadSlaves = 0;
   fAllMonitor = 0;
   fDataReady = kFALSE;
   fBytesReady = 0;
   fTotalBytes = 0;
   fAvailablePackages = 0;
   fEnabledPackages = 0;
   fRunningDSets = 0;

   fCollectTimeout = -1;

   fManager = 0;
   fQueryMode = kSync;
   fDynamicStartup = kFALSE;

   fMergersSet = kFALSE;
   fMergersByHost = kFALSE;
   fMergers = 0;
   fMergersCount = -1;
   fLastAssignedMerger = 0;
   fWorkersToMerge = 0;
   fFinalizationRunning = kFALSE;

   fPerfTree = "";

   fWrksOutputReady = 0;

   fSelector = 0;

   fPrepTime = 0.;

   // Check if the user defined a list of environment variables to send over:
   // include them into the dedicated list
   if (gSystem->Getenv("PROOF_ENVVARS")) {
      TString envs(gSystem->Getenv("PROOF_ENVVARS")), env, envsfound;
      Int_t from = 0;
      while (envs.Tokenize(env, from, ",")) {
         if (!env.IsNull()) {
            if (!gSystem->Getenv(env)) {
               Warning("Init", "request for sending over undefined environemnt variable '%s' - ignoring", env.Data());
            } else {
               if (!envsfound.IsNull()) envsfound += ",";
               envsfound += env;
               TProof::DelEnvVar(env);
               TProof::AddEnvVar(env, gSystem->Getenv(env));
            }
         }
      }
      if (envsfound.IsNull()) {
         Warning("Init", "none of the requested env variables were found: '%s'", envs.Data());
      } else {
         Info("Init", "the following environment variables have been added to the list to be sent to the nodes: '%s'", envsfound.Data());
      }
   }

   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Clean up PROOF environment.

TProof::~TProof()
{
   if (fChains) {
      while (TChain *chain = dynamic_cast<TChain*> (fChains->First()) ) {
         // remove "chain" from list
         chain->SetProof(0);
         RemoveChain(chain);
      }
   }

   // remove links to packages enabled on the client
   if (TestBit(TProof::kIsClient)) {
      // iterate over all packages
      TList *epl = fPackMgr->GetListOfEnabled();
      if (epl) {
         TIter nxp(epl);
         while (TObjString *pck = (TObjString *)(nxp())) {
            FileStat_t stat;
            if (gSystem->GetPathInfo(pck->String(), stat) == 0) {
               // check if symlink, if so unlink
               // NOTE: GetPathInfo() returns 1 in case of symlink that does not point to
               // existing file or to a directory, but if fIsLink is true the symlink exists
               if (stat.fIsLink)
                  gSystem->Unlink(pck->String());
            }
         }
         epl->Delete();
         delete epl;
      }
   }

   Close();
   SafeDelete(fIntHandler);
   SafeDelete(fSlaves);
   SafeDelete(fActiveSlaves);
   SafeDelete(fInactiveSlaves);
   SafeDelete(fUniqueSlaves);
   SafeDelete(fAllUniqueSlaves);
   SafeDelete(fNonUniqueMasters);
   SafeDelete(fTerminatedSlaveInfos);
   SafeDelete(fBadSlaves);
   SafeDelete(fAllMonitor);
   SafeDelete(fActiveMonitor);
   SafeDelete(fUniqueMonitor);
   SafeDelete(fAllUniqueMonitor);
   SafeDelete(fSlaveInfo);
   SafeDelete(fChains);
   SafeDelete(fPlayer);
   SafeDelete(fFeedback);
   SafeDelete(fWaitingSlaves);
   SafeDelete(fAvailablePackages);
   SafeDelete(fEnabledPackages);
   SafeDelete(fLoadedMacros);
   SafeDelete(fPackMgr);
   SafeDelete(fRecvMessages);
   SafeDelete(fInputData);
   SafeDelete(fRunningDSets);
   SafeDelete(fEnabledPackagesOnCluster);
   if (fWrksOutputReady) {
      fWrksOutputReady->SetOwner(kFALSE);
      delete fWrksOutputReady;
   }
   if (fQueries) {
      fQueries->Delete();
      delete fQueries;
   }

   // remove file with redirected logs
   if (TestBit(TProof::kIsClient)) {
      if (fLogFileR)
         fclose(fLogFileR);
      if (fLogFileW)
         fclose(fLogFileW);
      if (fLogFileName.Length() > 0)
         gSystem->Unlink(fLogFileName);
   }

   // Remove for the global list
   gROOT->GetListOfProofs()->Remove(this);
   // ... and from the manager list
   if (fManager && fManager->IsValid())
      fManager->DiscardSession(this);

   if (gProof && gProof == this) {
      // Set previous as default
      TIter pvp(gROOT->GetListOfProofs(), kIterBackward);
      while ((gProof = (TProof *)pvp())) {
         if (gProof->InheritsFrom(TProof::Class()))
            break;
      }
   }

   // For those interested in our destruction ...
   Emit("~TProof()");
   Emit("CloseWindow()");
}

////////////////////////////////////////////////////////////////////////////////
/// Start the PROOF environment. Starting PROOF involves either connecting
/// to a master server, which in turn will start a set of slave servers, or
/// directly starting as master server (if master = ""). For a description
/// of the arguments see the TProof ctor. Returns the number of started
/// master or slave servers, returns 0 in case of error, in which case
/// fValid remains false.

Int_t TProof::Init(const char *, const char *conffile,
                   const char *confdir, Int_t loglevel, const char *alias)
{
   R__ASSERT(gSystem);

   fValid = kFALSE;

   // Connected to terminal?
   fTty = (isatty(0) == 0 || isatty(1) == 0) ? kFALSE : kTRUE;

   // If in attach mode, options is filled with additional info
   Bool_t attach = kFALSE;
   if (strlen(fUrl.GetOptions()) > 0) {
      attach = kTRUE;
      // A flag from the GUI
      TString opts = fUrl.GetOptions();
      if (opts.Contains("GUI")) {
         SetBit(TProof::kUsingSessionGui);
         opts.Remove(opts.Index("GUI"));
         fUrl.SetOptions(opts);
      }
   }

   if (TestBit(TProof::kIsMaster)) {
      // Fill default conf file and conf dir
      if (!conffile || !conffile[0])
         fConfFile = kPROOF_ConfFile;
      if (!confdir  || !confdir[0])
         fConfDir  = kPROOF_ConfDir;
      // The group; the client receives it in the kPROOF_SESSIONTAG message
      if (gProofServ) fGroup = gProofServ->GetGroup();
   } else {
      fConfDir     = confdir;
      fConfFile    = conffile;
   }

   // Analysise the conffile field
   if (fConfFile.Contains("workers=0")) fConfFile.ReplaceAll("workers=0", "masteronly");
   ParseConfigField(fConfFile);

   fWorkDir        = gSystem->WorkingDirectory();
   fLogLevel       = loglevel;
   fProtocol       = kPROOF_Protocol;
   fSendGroupView  = kTRUE;
   fImage          = fMasterServ ? "" : "<local>";
   fIntHandler     = 0;
   fStatus         = 0;
   fRecvMessages   = new TList;
   fRecvMessages->SetOwner(kTRUE);
   fSlaveInfo      = 0;
   fChains         = new TList;
   fAvailablePackages = 0;
   fEnabledPackages = 0;
   fRunningDSets   = 0;
   fEndMaster      = TestBit(TProof::kIsMaster) ? kTRUE : kFALSE;
   fInputData      = 0;
   ResetBit(TProof::kNewInputData);
   fPrintProgress  = 0;

   // Timeout for some collect actions
   fCollectTimeout = gEnv->GetValue("Proof.CollectTimeout", -1);

   // Should the workers be started dynamically; default: no
   fDynamicStartup = gEnv->GetValue("Proof.DynamicStartup", kFALSE);

   // Default entry point for the data pool is the master
   if (TestBit(TProof::kIsClient))
      fDataPoolUrl.Form("root://%s", fMaster.Data());
   else
      fDataPoolUrl = "";

   fProgressDialog        = 0;
   fProgressDialogStarted = kFALSE;

   // Default alias is the master name
   TString      al = (alias) ? alias : fMaster.Data();
   SetAlias(al);

   // Client logging of messages from the master and slaves
   fRedirLog = kFALSE;
   if (TestBit(TProof::kIsClient)) {
      fLogFileName.Form("%s/ProofLog_%d", gSystem->TempDirectory(), gSystem->GetPid());
      if ((fLogFileW = fopen(fLogFileName, "w")) == 0)
         Error("Init", "could not create temporary logfile");
      if ((fLogFileR = fopen(fLogFileName, "r")) == 0)
         Error("Init", "could not open temp logfile for reading");
   }
   fLogToWindowOnly = kFALSE;

   // Status of cluster
   fNotIdle = 0;
   // Query type
   fSync = (attach) ? kFALSE : kTRUE;
   // Not enqueued
   fIsWaiting = kFALSE;

   // Counters
   fBytesRead = 0;
   fRealTime = 0;
   fCpuTime = 0;

   // List of queries
   fQueries = 0;
   fOtherQueries = 0;
   fDrawQueries = 0;
   fMaxDrawQueries = 1;
   fSeqNum = 0;

   // Remote ID of the session
   fSessionID = -1;

   // Part of active query
   fWaitingSlaves = 0;

   // Make remote PROOF player
   fPlayer = 0;
   MakePlayer();

   fFeedback = new TList;
   fFeedback->SetOwner();
   fFeedback->SetName("FeedbackList");
   AddInput(fFeedback);

   // sort slaves by descending performance index
   fSlaves           = new TSortedList(kSortDescending);
   fActiveSlaves     = new TList;
   fInactiveSlaves   = new TList;
   fUniqueSlaves     = new TList;
   fAllUniqueSlaves  = new TList;
   fNonUniqueMasters = new TList;
   fBadSlaves        = new TList;
   fAllMonitor       = new TMonitor;
   fActiveMonitor    = new TMonitor;
   fUniqueMonitor    = new TMonitor;
   fAllUniqueMonitor = new TMonitor;
   fCurrentMonitor   = 0;

   fTerminatedSlaveInfos = new TList;
   fTerminatedSlaveInfos->SetOwner(kTRUE);

   fLoadedMacros            = 0;
   fPackMgr                 = 0;

   // Enable optimized sending of streamer infos to use embedded backward/forward
   // compatibility support between different ROOT versions and different versions of
   // users classes
   Bool_t enableSchemaEvolution = gEnv->GetValue("Proof.SchemaEvolution",1);
   if (enableSchemaEvolution) {
      TMessage::EnableSchemaEvolutionForAll();
   } else {
      Info("TProof", "automatic schema evolution in TMessage explicitly disabled");
   }

   if (IsMaster()) {
      // to make UploadPackage() method work on the master as well.
      if (gProofServ) fPackMgr = gProofServ->GetPackMgr();
   } else {

      TString sandbox;
      if (GetSandbox(sandbox, kTRUE) != 0) {
         Error("Init", "failure asserting sandbox directory %s", sandbox.Data());
         return 0;
      }

      // Package Dir
      TString packdir = gEnv->GetValue("Proof.PackageDir", "");
      if (packdir.IsNull())
         packdir.Form("%s/%s", sandbox.Data(), kPROOF_PackDir);
      if (AssertPath(packdir, kTRUE) != 0) {
         Error("Init", "failure asserting directory %s", packdir.Data());
         return 0;
      }
      fPackMgr = new TPackMgr(packdir);
      if (gDebug > 0)
         Info("Init", "package directory set to %s", packdir.Data());
   }

   if (!IsMaster()) {
      // List of directories where to look for global packages
      TString globpack = gEnv->GetValue("Proof.GlobalPackageDirs","");
      TProofServ::ResolveKeywords(globpack);
      Int_t nglb = TPackMgr::RegisterGlobalPath(globpack);
      if (gDebug > 0)
         Info("Init", " %d global package directories registered", nglb);
   }

   // Master may want dynamic startup
   if (fDynamicStartup) {
      if (!IsMaster()) {
         // If on client - start the master
         if (!StartSlaves(attach))
            return 0;
      }
   } else {

      // Master Only mode (for operations requiring only the master, e.g. dataset browsing,
      // result retrieving, ...)
      Bool_t masterOnly = gEnv->GetValue("Proof.MasterOnly", kFALSE);
      if (!IsMaster() || !masterOnly) {
         // Start slaves (the old, static, per-session way)
         if (!StartSlaves(attach))
            return 0;
         // Client: Is Master in dynamic startup mode?
         if (!IsMaster()) {
            Int_t dyn = 0;
            GetRC("Proof.DynamicStartup", dyn);
            if (dyn != 0) fDynamicStartup = kTRUE;
         }
      }
   }
   // we are now properly initialized
   fValid = kTRUE;

   // De-activate monitor (will be activated in Collect)
   fAllMonitor->DeActivateAll();

   // By default go into parallel mode
   Int_t nwrk = GetRemoteProtocol() > 35 ? -1 : 9999;
   TNamed *n = 0;
   if (TProof::GetEnvVars() &&
      (n = (TNamed *) TProof::GetEnvVars()->FindObject("PROOF_NWORKERS"))) {
      TString s(n->GetTitle());
      if (s.IsDigit()) nwrk = s.Atoi();
   }
   GoParallel(nwrk, attach);

   // Send relevant initial state to slaves
   if (!attach)
      SendInitialState();
   else if (!IsIdle())
      // redirect log
      fRedirLog = kTRUE;

   // Done at this point, the alias will be communicated to the coordinator, if any
   if (TestBit(TProof::kIsClient))
      SetAlias(al);

   SetActive(kFALSE);

   if (IsValid()) {

      // Activate input handler
      ActivateAsyncInput();

      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfSockets()->Add(this);
   }

   AskParallel();

   return fActiveSlaves->GetSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the sandbox path from ' Proof.Sandbox' or the alternative var 'rc'.
/// Use the existing setting or the default if nothing is found.
/// If 'assert' is kTRUE, make also sure that the path exists.
/// Return 0 on success, -1 on failure

Int_t TProof::GetSandbox(TString &sb, Bool_t assert, const char *rc)
{
   // Get it from 'rc', if defined
   if (rc && strlen(rc)) sb = gEnv->GetValue(rc, sb);
   // Or use the default 'rc'
   if (sb.IsNull()) sb = gEnv->GetValue("Proof.Sandbox", "");
   // If nothing found , use the default
   if (sb.IsNull()) sb.Form("~/%s", kPROOF_WorkDir);
   // Expand special settings
   if (sb == ".") {
      sb = gSystem->pwd();
   } else if (sb == "..") {
      sb = gSystem->GetDirName(gSystem->pwd());
   }
   gSystem->ExpandPathName(sb);

   // Assert the path, if required
   if (assert && AssertPath(sb, kTRUE) != 0) return -1;
   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// The config file field may contain special instructions which need to be
/// parsed at the beginning, e.g. for debug runs with valgrind.
/// Several options can be given separated by a ','

void TProof::ParseConfigField(const char *config)
{
   TString sconf(config), opt;
   Ssiz_t from = 0;
   Bool_t cpuPin = kFALSE;

   // Analysise the field
   const char *cq = (IsLite()) ? "\"" : "";
   while (sconf.Tokenize(opt, from, ",")) {
      if (opt.IsNull()) continue;

      if (opt.BeginsWith("valgrind")) {
         // Any existing valgrind setting? User can give full settings, which we fully respect,
         // or pass additional options for valgrind by prefixing 'valgrind_opts:'. For example,
         //    TProof::AddEnvVar("PROOF_MASTER_WRAPPERCMD", "valgrind_opts:--time-stamp --leak-check=full"
         // will add option "--time-stamp --leak-check=full" to our default options
         TString mst, top, sub, wrk, all;
         TList *envs = fgProofEnvList;
         TNamed *n = 0;
         if (envs) {
            if ((n = (TNamed *) envs->FindObject("PROOF_WRAPPERCMD")))
               all = n->GetTitle();
            if ((n = (TNamed *) envs->FindObject("PROOF_MASTER_WRAPPERCMD")))
               mst = n->GetTitle();
            if ((n = (TNamed *) envs->FindObject("PROOF_TOPMASTER_WRAPPERCMD")))
               top = n->GetTitle();
            if ((n = (TNamed *) envs->FindObject("PROOF_SUBMASTER_WRAPPERCMD")))
               sub = n->GetTitle();
            if ((n = (TNamed *) envs->FindObject("PROOF_SLAVE_WRAPPERCMD")))
               wrk = n->GetTitle();
         }
         if (all != "" && mst == "") mst = all;
         if (all != "" && top == "") top = all;
         if (all != "" && sub == "") sub = all;
         if (all != "" && wrk == "") wrk = all;
         if (all != "" && all.BeginsWith("valgrind_opts:")) {
            // The field is used to add an option Reset the setting
            Info("ParseConfigField","valgrind run: resetting 'PROOF_WRAPPERCMD':"
                                    " must be set again for next run , if any");
            TProof::DelEnvVar("PROOF_WRAPPERCMD");
         }
         TString var, cmd;
         cmd.Form("%svalgrind -v --suppressions=<rootsys>/etc/valgrind-root.supp", cq);
         TString mstlab("NO"), wrklab("NO");
         Bool_t doMaster = (opt == "valgrind" || (opt.Contains("master") &&
                           !opt.Contains("topmaster") && !opt.Contains("submaster")))
                         ? kTRUE : kFALSE;
         if (doMaster) {
            if (!IsLite()) {
               // Check if we have to add a var
               if (mst == "" || mst.BeginsWith("valgrind_opts:")) {
                  mst.ReplaceAll("valgrind_opts:","");
                  var.Form("%s --log-file=<logfilemst>.valgrind.log %s", cmd.Data(), mst.Data());
                  TProof::AddEnvVar("PROOF_MASTER_WRAPPERCMD", var);
                  mstlab = "YES";
               } else if (mst != "") {
                  mstlab = "YES";
               }
            } else {
               if (opt.Contains("master")) {
                  Warning("ParseConfigField",
                        "master valgrinding does not make sense for PROOF-Lite: ignoring");
                  opt.ReplaceAll("master", "");
                  if (!opt.Contains("workers")) return;
               }
               if (opt == "valgrind" || opt == "valgrind=") opt = "valgrind=workers";
            }
         }
         if (opt.Contains("topmaster")) {
            // Check if we have to add a var
            if (top == "" || top.BeginsWith("valgrind_opts:")) {
               top.ReplaceAll("valgrind_opts:","");
               var.Form("%s --log-file=<logfilemst>.valgrind.log %s", cmd.Data(), top.Data());
               TProof::AddEnvVar("PROOF_TOPMASTER_WRAPPERCMD", var);
               mstlab = "YES";
            } else if (top != "") {
               mstlab = "YES";
            }
         }
         if (opt.Contains("submaster")) {
            // Check if we have to add a var
            if (sub == "" || sub.BeginsWith("valgrind_opts:")) {
               sub.ReplaceAll("valgrind_opts:","");
               var.Form("%s --log-file=<logfilemst>.valgrind.log %s", cmd.Data(), sub.Data());
               TProof::AddEnvVar("PROOF_SUBMASTER_WRAPPERCMD", var);
               mstlab = "YES";
            } else if (sub != "") {
               mstlab = "YES";
            }
         }
         if (opt.Contains("=workers") || opt.Contains("+workers")) {
            // Check if we have to add a var
            if (wrk == "" || wrk.BeginsWith("valgrind_opts:")) {
               wrk.ReplaceAll("valgrind_opts:","");
               var.Form("%s --log-file=<logfilewrk>.__valgrind__.log %s%s", cmd.Data(), wrk.Data(), cq);
               TProof::AddEnvVar("PROOF_SLAVE_WRAPPERCMD", var);
               TString nwrks("2");
               Int_t inw = opt.Index('#');
               if (inw != kNPOS) {
                  nwrks = opt(inw+1, opt.Length());
                  if (!nwrks.IsDigit()) nwrks = "2";
               }
               // Set the relevant variables
               if (!IsLite()) {
                  TProof::AddEnvVar("PROOF_NWORKERS", nwrks);
               } else {
                  gEnv->SetValue("ProofLite.Workers", nwrks.Atoi());
               }
               wrklab = nwrks;
               // Register the additional worker log in the session file
               // (for the master this is done automatically)
               TProof::AddEnvVar("PROOF_ADDITIONALLOG", "__valgrind__.log*");
            } else if (wrk != "") {
               wrklab = "ALL";
            }
         }
         // Increase the relevant timeouts
         if (!IsLite()) {
            TProof::AddEnvVar("PROOF_INTWAIT", "5000");
            gEnv->SetValue("Proof.SocketActivityTimeout", 6000);
         } else {
            gEnv->SetValue("ProofLite.StartupTimeOut", 5000);
         }
         // Warn for slowness
         Printf(" ");
         if (!IsLite()) {
            Printf(" ---> Starting a debug run with valgrind (master:%s, workers:%s)", mstlab.Data(), wrklab.Data());
         } else {
            Printf(" ---> Starting a debug run with valgrind (workers:%s)", wrklab.Data());
         }
         Printf(" ---> Please be patient: startup may be VERY slow ...");
         Printf(" ---> Logs will be available as special tags in the log window (from the progress dialog or TProof::LogViewer()) ");
         Printf(" ---> (Reminder: this debug run makes sense only if you are running a debug version of ROOT)");
         Printf(" ");

      } else if (opt.BeginsWith("igprof-pp")) {

         // IgProf profiling on master and worker. PROOF does not set the
         // environment for you: proper environment variables (like PATH and
         // LD_LIBRARY_PATH) should be set externally

         Printf("*** Requested IgProf performance profiling ***");
         TString addLogExt = "__igprof.pp__.log";
         TString addLogFmt = "igprof -pk -pp -t proofserv.exe -o %s.%s";
         TString tmp;

         if (IsLite()) {
            addLogFmt.Append("\"");
            addLogFmt.Prepend("\"");
         }

         tmp.Form(addLogFmt.Data(), "<logfilemst>", addLogExt.Data());
         TProof::AddEnvVar("PROOF_MASTER_WRAPPERCMD",  tmp.Data());

         tmp.Form(addLogFmt.Data(), "<logfilewrk>", addLogExt.Data());
         TProof::AddEnvVar("PROOF_SLAVE_WRAPPERCMD", tmp.Data() );

         TProof::AddEnvVar("PROOF_ADDITIONALLOG", addLogExt.Data());

      } else if (opt.BeginsWith("cpupin=")) {
         // Enable CPU pinning. Takes as argument the list of processor IDs
         // that will be used in order. Processor IDs are numbered from 0,
         // use likwid to see how they are organized. A possible parameter
         // format would be:
         //
         //   cpupin=3+4+0+9+10+22+7
         //
         // Only the specified processor IDs will be used in a round-robin
         // fashion, dealing with the fact that you can request more workers
         // than the number of processor IDs you have specified.
         //
         // To use all available processors in their order:
         //
         //   cpupin=*

         opt.Remove(0, 7);

         // Remove any char which is neither a number nor a plus '+'
         for (Ssiz_t i=0; i<opt.Length(); i++) {
            Char_t c = opt[i];
            if ((c != '+') && ((c < '0') || (c > '9')))
               opt[i] = '_';
         }
         opt.ReplaceAll("_", "");
         TProof::AddEnvVar("PROOF_SLAVE_CPUPIN_ORDER", opt);
         cpuPin = kTRUE;
      } else if (opt.BeginsWith("workers=")) {

         // Request for a given number of workers (within the max) or worker
         // startup combination:
         //      workers=5         start max 5 workers (or less, if less are assigned)
         //      workers=2x        start max 2 workers per node (or less, if less are assigned)
         opt.ReplaceAll("workers=","");
         TProof::AddEnvVar("PROOF_NWORKERS", opt);
      }
   }

   // In case of PROOF-Lite, enable CPU pinning when requested (Linux only)
   #ifdef R__LINUX
   if (IsLite() && cpuPin) {
      Printf("*** Requested CPU pinning ***");
      const TList *ev = GetEnvVars();
      const char *pinCmd = "taskset -c <cpupin>";
      TString val;
      TNamed *p;
      if (ev && (p = dynamic_cast<TNamed *>(ev->FindObject("PROOF_SLAVE_WRAPPERCMD")))) {
         val = p->GetTitle();
         val.Insert(val.Length()-1, " ");
         val.Insert(val.Length()-1, pinCmd);
      }
      else {
         val.Form("\"%s\"", pinCmd);
      }
      TProof::AddEnvVar("PROOF_SLAVE_WRAPPERCMD", val.Data());
   }
   #endif
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that 'path' exists; if 'writable' is kTRUE, make also sure
/// that the path is writable

Int_t TProof::AssertPath(const char *inpath, Bool_t writable)
{
   if (!inpath || strlen(inpath) <= 0) {
      Error("AssertPath", "undefined input path");
      return -1;
   }

   TString path(inpath);
   gSystem->ExpandPathName(path);

   if (gSystem->AccessPathName(path, kFileExists)) {
      if (gSystem->mkdir(path, kTRUE) != 0) {
         Error("AssertPath", "could not create path %s", path.Data());
         return -1;
      }
   }
   // It must be writable
   if (gSystem->AccessPathName(path, kWritePermission) && writable) {
      if (gSystem->Chmod(path, 0666) != 0) {
         Error("AssertPath", "could not make path %s writable", path.Data());
         return -1;
      }
   }

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set manager and schedule its destruction after this for clean
/// operations.

void TProof::SetManager(TProofMgr *mgr)
{
   fManager = mgr;

   if (mgr) {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfSockets()->Remove(mgr);
      gROOT->GetListOfSockets()->Add(mgr);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Works on the master node only.
/// It starts workers on the machines in workerList and sets the paths,
/// packages and macros as on the master.
/// It is a subbstitute for StartSlaves(...)
/// The code is mostly the master part of StartSlaves,
/// with the parallel startup removed.

Int_t TProof::AddWorkers(TList *workerList)
{
   if (!IsMaster()) {
      Error("AddWorkers", "AddWorkers can only be called on the master!");
      return -1;
   }

   if (!workerList || !(workerList->GetSize())) {
      Error("AddWorkers", "empty list of workers!");
      return -2;
   }

   // Code taken from master part of StartSlaves with the parllel part removed

   fImage = gProofServ->GetImage();
   if (fImage.IsNull())
      fImage.Form("%s:%s", TUrl(gSystem->HostName()).GetHostFQDN(), gProofServ->GetWorkDir());

   // Get all workers
   UInt_t nSlaves = workerList->GetSize();
   UInt_t nSlavesDone = 0;
   Int_t ord = 0;

   // Loop over all new workers and start them (if we had already workers it means we are
   // increasing parallelism or that is not the first time we are called)
   Bool_t goMoreParallel = (fSlaves->GetEntries() > 0) ? kTRUE : kFALSE;

   // A list of TSlave objects for workers that are being added
   TList *addedWorkers = new TList();
   if (!addedWorkers) {
      // This is needed to silence Coverity ...
      Error("AddWorkers", "cannot create new list for the workers to be added");
      return -2;
   }
   addedWorkers->SetOwner(kFALSE);
   TListIter next(workerList);
   TObject *to;
   TProofNodeInfo *worker;
   TSlaveInfo *dummysi = new TSlaveInfo();
   while ((to = next())) {
      // Get the next worker from the list
      worker = (TProofNodeInfo *)to;

      // Read back worker node info
      const Char_t *image = worker->GetImage().Data();
      const Char_t *workdir = worker->GetWorkDir().Data();
      Int_t perfidx = worker->GetPerfIndex();
      Int_t sport = worker->GetPort();
      if (sport == -1)
         sport = fUrl.GetPort();

      // Create worker server
      TString fullord;
      if (worker->GetOrdinal().Length() > 0) {
         fullord.Form("%s.%s", gProofServ->GetOrdinal(), worker->GetOrdinal().Data());
      } else {
         fullord.Form("%s.%d", gProofServ->GetOrdinal(), ord);
      }

      // Remove worker from the list of workers terminated gracefully
      dummysi->SetOrdinal(fullord);
      TSlaveInfo *rmsi = (TSlaveInfo *)fTerminatedSlaveInfos->Remove(dummysi);
      SafeDelete(rmsi);

      // Create worker server
      TString wn(worker->GetNodeName());
      if (wn == "localhost" || wn.BeginsWith("localhost.")) wn = gSystem->HostName();
      TUrl u(TString::Format("%s:%d", wn.Data(), sport));
      // Add group info in the password firdl, if any
      if (strlen(gProofServ->GetGroup()) > 0) {
         // Set also the user, otherwise the password is not exported
         if (strlen(u.GetUser()) <= 0)
            u.SetUser(gProofServ->GetUser());
         u.SetPasswd(gProofServ->GetGroup());
      }
      TSlave *slave = 0;
      if (worker->IsWorker()) {
         slave = CreateSlave(u.GetUrl(), fullord, perfidx, image, workdir);
      } else {
         slave = CreateSubmaster(u.GetUrl(), fullord,
                                 image, worker->GetMsd(), worker->GetNWrks());
      }

      // Add to global list (we will add to the monitor list after
      // finalizing the server startup)
      Bool_t slaveOk = kTRUE;
      fSlaves->Add(slave);
      if (slave->IsValid()) {
         addedWorkers->Add(slave);
      } else {
         slaveOk = kFALSE;
         fBadSlaves->Add(slave);
         Warning("AddWorkers", "worker '%s' is invalid", slave->GetOrdinal());
      }

      PDB(kGlobal,3)
         Info("AddWorkers", "worker on host %s created"
              " and added to list (ord: %s)", worker->GetName(), slave->GetOrdinal());

      // Notify opening of connection
      nSlavesDone++;
      TMessage m(kPROOF_SERVERSTARTED);
      m << TString("Opening connections to workers") << nSlaves
        << nSlavesDone << slaveOk;
      gProofServ->GetSocket()->Send(m);

      ord++;
   } //end of the worker loop
   SafeDelete(dummysi);

   // Cleanup
   SafeDelete(workerList);

   nSlavesDone = 0;

   // Here we finalize the server startup: in this way the bulk
   // of remote operations are almost parallelized
   TIter nxsl(addedWorkers);
   TSlave *sl = 0;
   while ((sl = (TSlave *) nxsl())) {

      // Finalize setup of the server
      if (sl->IsValid())
          sl->SetupServ(TSlave::kSlave, 0);

      // Monitor good slaves
      Bool_t slaveOk = kTRUE;
      if (sl->IsValid()) {
         fAllMonitor->Add(sl->GetSocket());
      PDB(kGlobal,3)
         Info("AddWorkers", "worker on host %s finalized"
              " and added to list", sl->GetOrdinal());
      } else {
         slaveOk = kFALSE;
         fBadSlaves->Add(sl);
      }

      // Notify end of startup operations
      nSlavesDone++;
      TMessage m(kPROOF_SERVERSTARTED);
      m << TString("Setting up worker servers") << nSlaves
        << nSlavesDone << slaveOk;
      gProofServ->GetSocket()->Send(m);
   }

   // Now set new state on the added workers (on all workers for simplicity)
   // use fEnabledPackages, fLoadedMacros,
   // gSystem->GetDynamicPath() and gSystem->GetIncludePath()
   // no need to load packages that are only loaded and not enabled (dyn mode)
   Int_t nwrk = GetRemoteProtocol() > 35 ? -1 : 9999;
   TNamed *n = 0;
   if (TProof::GetEnvVars() &&
      (n = (TNamed *) TProof::GetEnvVars()->FindObject("PROOF_NWORKERS"))) {
      TString s(n->GetTitle());
      if (s.IsDigit()) nwrk = s.Atoi();
   }

   if (fDynamicStartup && goMoreParallel) {

      PDB(kGlobal, 3)
         Info("AddWorkers", "will invoke GoMoreParallel()");
      Int_t nw = GoMoreParallel(nwrk);
      PDB(kGlobal, 3)
         Info("AddWorkers", "GoMoreParallel()=%d", nw);

   }
   else {
      // Not in Dynamic Workers mode
      PDB(kGlobal, 3)
         Info("AddWorkers", "will invoke GoParallel()");
      GoParallel(nwrk, kFALSE, 0);
   }

   // Set worker processing environment
   SetupWorkersEnv(addedWorkers, goMoreParallel);

   // Update list of current workers
   PDB(kGlobal, 3)
      Info("AddWorkers", "will invoke SaveWorkerInfo()");
   SaveWorkerInfo();

   // Inform the client that the number of workers has changed
   if (fDynamicStartup && gProofServ) {
      PDB(kGlobal, 3)
         Info("AddWorkers", "will invoke SendParallel()");
      gProofServ->SendParallel(kTRUE);

      if (goMoreParallel && fPlayer) {
         // In case we are adding workers dynamically to an existing process, we
         // should invoke a special player's Process() to set only added workers
         // to the proper state
         PDB(kGlobal, 3)
            Info("AddWorkers", "will send the PROCESS message to selected workers");
         fPlayer->JoinProcess(addedWorkers);
         // Update merger counters (new workers are not yet active)
         fMergePrg.SetNWrks(fActiveSlaves->GetSize() + addedWorkers->GetSize());
      }
   }

   // Cleanup
   delete addedWorkers;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set up packages, loaded macros, include and lib paths ...

void TProof::SetupWorkersEnv(TList *addedWorkers, Bool_t increasingWorkers)
{
   TList *server_packs = gProofServ ? gProofServ->GetEnabledPackages() : nullptr;
   // Packages
   TList *packs = server_packs ? server_packs : GetEnabledPackages();
   if (packs && packs->GetSize() > 0) {
      TIter nxp(packs);
      TPair *pck = 0;
      while ((pck = (TPair *) nxp())) {
         // Upload and Enable methods are intelligent and avoid
         // re-uploading or re-enabling of a package to a node that has it.
         if (fDynamicStartup && increasingWorkers) {
            // Upload only on added workers
            PDB(kGlobal, 3)
               Info("SetupWorkersEnv", "will invoke UploadPackage() and EnablePackage() on added workers");
            if (UploadPackage(pck->GetName(), kUntar, addedWorkers) >= 0)
               EnablePackage(pck->GetName(), (TList *) pck->Value(), kTRUE, addedWorkers);
         } else {
            PDB(kGlobal, 3)
               Info("SetupWorkersEnv", "will invoke UploadPackage() and EnablePackage() on all workers");
            if (UploadPackage(pck->GetName()) >= 0)
               EnablePackage(pck->GetName(), (TList *) pck->Value(), kTRUE);
         }
      }
   }

   if (server_packs) {
      server_packs->Delete();
      delete server_packs;
   }

   // Loaded macros
   if (fLoadedMacros) {
      TIter nxp(fLoadedMacros);
      TObjString *os = 0;
      while ((os = (TObjString *) nxp())) {
         PDB(kGlobal, 3) {
            Info("SetupWorkersEnv", "will invoke Load() on selected workers");
            Printf("Loading a macro : %s", os->GetName());
         }
         Load(os->GetName(), kTRUE, kTRUE, addedWorkers);
      }
   }

   // Dynamic path
   TString dyn = gSystem->GetDynamicPath();
   dyn.ReplaceAll(":", " ");
   dyn.ReplaceAll("\"", " ");
   PDB(kGlobal, 3)
      Info("SetupWorkersEnv", "will invoke AddDynamicPath() on selected workers");
   AddDynamicPath(dyn, kFALSE, addedWorkers, kFALSE); // Do not Collect

   // Include path
   TString inc = gSystem->GetIncludePath();
   inc.ReplaceAll("-I", " ");
   inc.ReplaceAll("\"", " ");
   PDB(kGlobal, 3)
      Info("SetupWorkersEnv", "will invoke AddIncludePath() on selected workers");
   AddIncludePath(inc, kFALSE, addedWorkers, kFALSE);  // Do not Collect

   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Used for shuting down the workres after a query is finished.
/// Sends each of the workers from the workerList, a kPROOF_STOP message.
/// If the workerList == 0, shutdown all the workers.

Int_t TProof::RemoveWorkers(TList *workerList)
{
   if (!IsMaster()) {
      Error("RemoveWorkers", "RemoveWorkers can only be called on the master!");
      return -1;
   }

   fFileMap.clear(); // This could be avoided if CopyFromCache was used in SendFile

   if (!workerList) {
      // shutdown all the workers
      TIter nxsl(fSlaves);
      TSlave *sl = 0;
      while ((sl = (TSlave *) nxsl())) {
         // Shut down the worker assumig that it is not processing
         TerminateWorker(sl);
      }

   } else {
      if (!(workerList->GetSize())) {
         Error("RemoveWorkers", "The list of workers should not be empty!");
         return -2;
      }

      // Loop over all the workers and stop them
      TListIter next(workerList);
      TObject *to;
      TProofNodeInfo *worker;
      while ((to = next())) {
         TSlave *sl = 0;
         if (!strcmp(to->ClassName(), "TProofNodeInfo")) {
            // Get the next worker from the list
            worker = (TProofNodeInfo *)to;
            TIter nxsl(fSlaves);
            while ((sl = (TSlave *) nxsl())) {
               // Shut down the worker assumig that it is not processing
               if (sl->GetName() == worker->GetNodeName())
                  break;
            }
         } else if (to->InheritsFrom(TSlave::Class())) {
            sl = (TSlave *) to;
         } else {
            Warning("RemoveWorkers","unknown object type: %s - it should be"
                  " TProofNodeInfo or inheriting from TSlave", to->ClassName());
         }
         // Shut down the worker assumig that it is not processing
         if (sl) {
            if (gDebug > 0)
               Info("RemoveWorkers","terminating worker %s", sl->GetOrdinal());
            TerminateWorker(sl);
         }
      }
   }

   // Update also the master counter
   if (gProofServ && fSlaves->GetSize() <= 0) gProofServ->ReleaseWorker("master");

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Start up PROOF slaves.

Bool_t TProof::StartSlaves(Bool_t attach)
{
   // If this is a master server, find the config file and start slave
   // servers as specified in the config file
   if (TestBit(TProof::kIsMaster)) {

      Int_t pc = 0;
      TList *workerList = new TList;
      // Get list of workers
      if (gProofServ->GetWorkers(workerList, pc) == TProofServ::kQueryStop) {
         TString emsg("no resource currently available for this session: please retry later");
         if (gDebug > 0) Info("StartSlaves", "%s", emsg.Data());
         gProofServ->SendAsynMessage(emsg.Data());
         return kFALSE;
      }
      // Setup the workers
      if (AddWorkers(workerList) < 0)
         return kFALSE;

   } else {

      // create master server
      Printf("Starting master: opening connection ...");
      TSlave *slave = CreateSubmaster(fUrl.GetUrl(), "0", "master", 0);

      if (slave->IsValid()) {

         // Notify
         fprintf(stderr,"Starting master:"
                        " connection open: setting up server ...             \r");
         StartupMessage("Connection to master opened", kTRUE, 1, 1);

         if (!attach) {

            // Set worker interrupt handler
            slave->SetInterruptHandler(kTRUE);

            // Finalize setup of the server
            slave->SetupServ(TSlave::kMaster, fConfFile);

            if (slave->IsValid()) {

               // Notify
               Printf("Starting master: OK                                     ");
               StartupMessage("Master started", kTRUE, 1, 1);

               // check protocol compatibility
               // protocol 1 is not supported anymore
               if (fProtocol == 1) {
                  Error("StartSlaves",
                        "client and remote protocols not compatible (%d and %d)",
                        kPROOF_Protocol, fProtocol);
                  slave->Close("S");
                  delete slave;
                  return kFALSE;
               }

               fSlaves->Add(slave);
               fAllMonitor->Add(slave->GetSocket());

               // Unset worker interrupt handler
               slave->SetInterruptHandler(kFALSE);

               // Set interrupt PROOF handler from now on
               fIntHandler = new TProofInterruptHandler(this);

               // Give-up after 5 minutes
               Int_t rc = Collect(slave, 300);
               Int_t slStatus = slave->GetStatus();
               if (slStatus == -99 || slStatus == -98 || rc == 0) {
                  fSlaves->Remove(slave);
                  fAllMonitor->Remove(slave->GetSocket());
                  if (slStatus == -99)
                     Error("StartSlaves", "no resources available or problems setting up workers (check logs)");
                  else if (slStatus == -98)
                     Error("StartSlaves", "could not setup output redirection on master");
                  else
                     Error("StartSlaves", "setting up master");
                  slave->Close("S");
                  delete slave;
                  return 0;
               }

               if (!slave->IsValid()) {
                  fSlaves->Remove(slave);
                  fAllMonitor->Remove(slave->GetSocket());
                  slave->Close("S");
                  delete slave;
                  Error("StartSlaves",
                        "failed to setup connection with PROOF master server");
                  return kFALSE;
               }

               if (!gROOT->IsBatch() && TestBit(kUseProgressDialog)) {
                  if ((fProgressDialog =
                     gROOT->GetPluginManager()->FindHandler("TProofProgressDialog")))
                     if (fProgressDialog->LoadPlugin() == -1)
                        fProgressDialog = 0;
               }
            } else {
               // Notify
               Printf("Starting master: failure");
            }
         } else {

            // Notify
            Printf("Starting master: OK                                     ");
            StartupMessage("Master attached", kTRUE, 1, 1);

            if (!gROOT->IsBatch() && TestBit(kUseProgressDialog)) {
               if ((fProgressDialog =
                  gROOT->GetPluginManager()->FindHandler("TProofProgressDialog")))
                  if (fProgressDialog->LoadPlugin() == -1)
                     fProgressDialog = 0;
            }

            fSlaves->Add(slave);
            fIntHandler = new TProofInterruptHandler(this);
         }

      } else {
         delete slave;
         // Notify only if verbosity is on: most likely the failure has already been notified
         if (gDebug > 0)
            Error("StartSlaves", "failed to create (or connect to) the PROOF master server");
         return kFALSE;
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Close all open slave servers.
/// Client can decide to shutdown the remote session by passing option is 'S'
/// or 's'. Default for clients is detach, if supported. Masters always
/// shutdown the remote counterpart.

void TProof::Close(Option_t *opt)
{
   {  std::lock_guard<std::recursive_mutex> lock(fCloseMutex);

      fValid = kFALSE;
      if (fSlaves) {
         if (fIntHandler)
            fIntHandler->Remove();

         TIter nxs(fSlaves);
         TSlave *sl = 0;
         while ((sl = (TSlave *)nxs()))
            sl->Close(opt);

         fActiveSlaves->Clear("nodelete");
         fUniqueSlaves->Clear("nodelete");
         fAllUniqueSlaves->Clear("nodelete");
         fNonUniqueMasters->Clear("nodelete");
         fBadSlaves->Clear("nodelete");
         fInactiveSlaves->Clear("nodelete");
         fSlaves->Delete();
      }
   }

   {  R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfSockets()->Remove(this);

      if (fChains) {
         while (TChain *chain = dynamic_cast<TChain*> (fChains->First()) ) {
            // remove "chain" from list
            chain->SetProof(0);
            RemoveChain(chain);
         }
      }

      if (IsProofd()) {

         gROOT->GetListOfProofs()->Remove(this);
         if (gProof && gProof == this) {
            // Set previous proofd-related as default
            TIter pvp(gROOT->GetListOfProofs(), kIterBackward);
            while ((gProof = (TProof *)pvp())) {
               if (gProof->IsProofd())
                  break;
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new TSlave of type TSlave::kSlave.
/// Note: creation of TSlave is private with TProof as a friend.
/// Derived classes must use this function to create slaves.

TSlave *TProof::CreateSlave(const char *url, const char *ord,
                            Int_t perf, const char *image, const char *workdir)
{
   TSlave* sl = TSlave::Create(url, ord, perf, image,
                               this, TSlave::kSlave, workdir, 0);

   if (sl->IsValid()) {
      sl->SetInputHandler(new TProofInputHandler(this, sl->GetSocket()));
      // must set fParallel to 1 for slaves since they do not
      // report their fParallel with a LOG_DONE message
      sl->fParallel = 1;
   }

   return sl;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a new TSlave of type TSlave::kMaster.
/// Note: creation of TSlave is private with TProof as a friend.
/// Derived classes must use this function to create slaves.

TSlave *TProof::CreateSubmaster(const char *url, const char *ord,
                                const char *image, const char *msd, Int_t nwk)
{
   TSlave *sl = TSlave::Create(url, ord, 100, image, this,
                               TSlave::kMaster, 0, msd, nwk);

   if (sl->IsValid()) {
      sl->SetInputHandler(new TProofInputHandler(this, sl->GetSocket()));
   }

   return sl;
}

////////////////////////////////////////////////////////////////////////////////
/// Find slave that has TSocket s. Returns 0 in case slave is not found.

TSlave *TProof::FindSlave(TSocket *s) const
{
   TSlave *sl;
   TIter   next(fSlaves);

   while ((sl = (TSlave *)next())) {
      if (sl->IsValid() && sl->GetSocket() == s)
         return sl;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add to the fUniqueSlave list the active slaves that have a unique
/// (user) file system image. This information is used to transfer files
/// only once to nodes that share a file system (an image). Submasters
/// which are not in fUniqueSlaves are put in the fNonUniqueMasters
/// list. That list is used to trigger the transferring of files to
/// the submaster's unique slaves without the need to transfer the file
/// to the submaster.

void TProof::FindUniqueSlaves()
{
   fUniqueSlaves->Clear();
   fUniqueMonitor->RemoveAll();
   fAllUniqueSlaves->Clear();
   fAllUniqueMonitor->RemoveAll();
   fNonUniqueMasters->Clear();

   TIter next(fActiveSlaves);

   while (TSlave *sl = dynamic_cast<TSlave*>(next())) {
      if (fImage == sl->fImage) {
         if (sl->GetSlaveType() == TSlave::kMaster) {
            fNonUniqueMasters->Add(sl);
            fAllUniqueSlaves->Add(sl);
            fAllUniqueMonitor->Add(sl->GetSocket());
         }
         continue;
      }

      TIter next2(fUniqueSlaves);
      TSlave *replace_slave = 0;
      Bool_t add = kTRUE;
      while (TSlave *sl2 = dynamic_cast<TSlave*>(next2())) {
         if (sl->fImage == sl2->fImage) {
            add = kFALSE;
            if (sl->GetSlaveType() == TSlave::kMaster) {
               if (sl2->GetSlaveType() == TSlave::kSlave) {
                  // give preference to master
                  replace_slave = sl2;
                  add = kTRUE;
               } else if (sl2->GetSlaveType() == TSlave::kMaster) {
                  fNonUniqueMasters->Add(sl);
                  fAllUniqueSlaves->Add(sl);
                  fAllUniqueMonitor->Add(sl->GetSocket());
               } else {
                  Error("FindUniqueSlaves", "TSlave is neither Master nor Slave");
                  R__ASSERT(0);
               }
            }
            break;
         }
      }

      if (add) {
         fUniqueSlaves->Add(sl);
         fAllUniqueSlaves->Add(sl);
         fUniqueMonitor->Add(sl->GetSocket());
         fAllUniqueMonitor->Add(sl->GetSocket());
         if (replace_slave) {
            fUniqueSlaves->Remove(replace_slave);
            fAllUniqueSlaves->Remove(replace_slave);
            fUniqueMonitor->Remove(replace_slave->GetSocket());
            fAllUniqueMonitor->Remove(replace_slave->GetSocket());
         }
      }
   }

   // will be actiavted in Collect()
   fUniqueMonitor->DeActivateAll();
   fAllUniqueMonitor->DeActivateAll();
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of slaves as described in the config file.

Int_t TProof::GetNumberOfSlaves() const
{
   return fSlaves->GetSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of active slaves, i.e. slaves that are valid and in
/// the current computing group.

Int_t TProof::GetNumberOfActiveSlaves() const
{
   return fActiveSlaves->GetSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of inactive slaves, i.e. slaves that are valid but not in
/// the current computing group.

Int_t TProof::GetNumberOfInactiveSlaves() const
{
   return fInactiveSlaves->GetSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of unique slaves, i.e. active slaves that have each a
/// unique different user files system.

Int_t TProof::GetNumberOfUniqueSlaves() const
{
   return fUniqueSlaves->GetSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of bad slaves. This are slaves that we in the config
/// file, but refused to startup or that died during the PROOF session.

Int_t TProof::GetNumberOfBadSlaves() const
{
   return fBadSlaves->GetSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Ask the for the statistics of the slaves.

void TProof::AskStatistics()
{
   if (!IsValid()) return;

   Broadcast(kPROOF_GETSTATS, kActive);
   Collect(kActive, fCollectTimeout);
}

////////////////////////////////////////////////////////////////////////////////
/// Get statistics about CPU time, real time and bytes read.
/// If verbose, print the resuls (always available via GetCpuTime(), GetRealTime()
/// and GetBytesRead()

void TProof::GetStatistics(Bool_t verbose)
{
   if (fProtocol > 27) {
      // This returns the correct result
      AskStatistics();
   } else {
      // AskStatistics is buggy: parse the output of Print()
      RedirectHandle_t rh;
      gSystem->RedirectOutput(fLogFileName, "a", &rh);
      Print();
      gSystem->RedirectOutput(0, 0, &rh);
      TMacro *mp = GetLastLog();
      if (mp) {
         // Look for global directories
         TIter nxl(mp->GetListOfLines());
         TObjString *os = 0;
         while ((os = (TObjString *) nxl())) {
            TString s(os->GetName());
            if (s.Contains("Total MB's processed:")) {
               s.ReplaceAll("Total MB's processed:", "");
               if (s.IsFloat()) fBytesRead = (Long64_t) s.Atof() * (1024*1024);
            } else if (s.Contains("Total real time used (s):")) {
               s.ReplaceAll("Total real time used (s):", "");
               if (s.IsFloat()) fRealTime = s.Atof();
            } else if (s.Contains("Total CPU time used (s):")) {
               s.ReplaceAll("Total CPU time used (s):", "");
               if (s.IsFloat()) fCpuTime = s.Atof();
            }
         }
         delete mp;
      }
   }

   if (verbose) {
      Printf(" Real/CPU time (s): %.3f / %.3f; workers: %d; processed: %.2f MBs",
             GetRealTime(), GetCpuTime(), GetParallel(), float(GetBytesRead())/(1024*1024));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Ask the for the number of parallel slaves.

void TProof::AskParallel()
{
   if (!IsValid()) return;

   Broadcast(kPROOF_GETPARALLEL, kActive);
   Collect(kActive, fCollectTimeout);
}

////////////////////////////////////////////////////////////////////////////////
/// Ask the master for the list of queries.

TList *TProof::GetListOfQueries(Option_t *opt)
{
   if (!IsValid() || TestBit(TProof::kIsMaster)) return (TList *)0;

   Bool_t all = ((strchr(opt,'A') || strchr(opt,'a'))) ? kTRUE : kFALSE;
   TMessage m(kPROOF_QUERYLIST);
   m << all;
   Broadcast(m, kActive);
   Collect(kActive, fCollectTimeout);

   // This should have been filled by now
   return fQueries;
}

////////////////////////////////////////////////////////////////////////////////
/// Number of queries processed by this session

Int_t TProof::GetNumberOfQueries()
{
   if (fQueries)
      return fQueries->GetSize() - fOtherQueries;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set max number of draw queries whose results are saved

void TProof::SetMaxDrawQueries(Int_t max)
{
   if (max > 0) {
      if (fPlayer)
         fPlayer->SetMaxDrawQueries(max);
      fMaxDrawQueries = max;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get max number of queries whose full results are kept in the
/// remote sandbox

void TProof::GetMaxQueries()
{
   TMessage m(kPROOF_MAXQUERIES);
   m << kFALSE;
   Broadcast(m, kActive);
   Collect(kActive, fCollectTimeout);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to the list of query results in the player

TList *TProof::GetQueryResults()
{
   return (fPlayer ? fPlayer->GetListOfResults() : (TList *)0);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to the full TQueryResult instance owned by the player
/// and referenced by 'ref'. If ref = 0 or "", return the last query result.

TQueryResult *TProof::GetQueryResult(const char *ref)
{
   return (fPlayer ? fPlayer->GetQueryResult(ref) : (TQueryResult *)0);
}

////////////////////////////////////////////////////////////////////////////////
/// Ask the master for the list of queries.
/// Options:
///           "A"     show information about all the queries known to the
///                   server, i.e. even those processed by other sessions
///           "L"     show only information about queries locally available
///                   i.e. already retrieved. If "L" is specified, "A" is
///                   ignored.
///           "F"     show all details available about queries
///           "H"     print help menu
/// Default ""

void TProof::ShowQueries(Option_t *opt)
{
   Bool_t help = ((strchr(opt,'H') || strchr(opt,'h'))) ? kTRUE : kFALSE;
   if (help) {

      // Help

      Printf("+++");
      Printf("+++ Options: \"A\" show all queries known to server");
      Printf("+++          \"L\" show retrieved queries");
      Printf("+++          \"F\" full listing of query info");
      Printf("+++          \"H\" print this menu");
      Printf("+++");
      Printf("+++ (case insensitive)");
      Printf("+++");
      Printf("+++ Use Retrieve(<#>) to retrieve the full"
             " query results from the master");
      Printf("+++     e.g. Retrieve(8)");

      Printf("+++");

      return;
   }

   if (!IsValid()) return;

   Bool_t local = ((strchr(opt,'L') || strchr(opt,'l'))) ? kTRUE : kFALSE;

   TObject *pq = 0;
   if (!local) {
      GetListOfQueries(opt);

      if (!fQueries) return;

      TIter nxq(fQueries);

      // Queries processed by other sessions
      if (fOtherQueries > 0) {
         Printf("+++");
         Printf("+++ Queries processed during other sessions: %d", fOtherQueries);
         Int_t nq = 0;
         while (nq++ < fOtherQueries && (pq = nxq()))
            pq->Print(opt);
      }

      // Queries processed by this session
      Printf("+++");
      Printf("+++ Queries processed during this session: selector: %d, draw: %d",
              GetNumberOfQueries(), fDrawQueries);
      while ((pq = nxq()))
         pq->Print(opt);

   } else {

      // Queries processed by this session
      Printf("+++");
      Printf("+++ Queries processed during this session: selector: %d, draw: %d",
              GetNumberOfQueries(), fDrawQueries);

      // Queries available locally
      TList *listlocal = fPlayer ? fPlayer->GetListOfResults() : (TList *)0;
      if (listlocal) {
         Printf("+++");
         Printf("+++ Queries available locally: %d", listlocal->GetSize());
         TIter nxlq(listlocal);
         while ((pq = nxlq()))
            pq->Print(opt);
      }
   }
   Printf("+++");
}

////////////////////////////////////////////////////////////////////////////////
/// See if the data is ready to be analyzed.

Bool_t TProof::IsDataReady(Long64_t &totalbytes, Long64_t &bytesready)
{
   if (!IsValid()) return kFALSE;

   TList submasters;
   TIter nextSlave(GetListOfActiveSlaves());
   while (TSlave *sl = dynamic_cast<TSlave*>(nextSlave())) {
      if (sl->GetSlaveType() == TSlave::kMaster) {
         submasters.Add(sl);
      }
   }

   fDataReady = kTRUE; //see if any submasters set it to false
   fBytesReady = 0;
   fTotalBytes = 0;
   //loop over submasters and see if data is ready
   if (submasters.GetSize() > 0) {
      Broadcast(kPROOF_DATA_READY, &submasters);
      Collect(&submasters);
   }

   bytesready = fBytesReady;
   totalbytes = fTotalBytes;

   EmitVA("IsDataReady(Long64_t,Long64_t)", 2, totalbytes, bytesready);

   PDB(kGlobal,2)
      Info("IsDataReady", "%lld / %lld (%s)",
           bytesready, totalbytes, fDataReady?"READY":"NOT READY");

   return fDataReady;
}

////////////////////////////////////////////////////////////////////////////////
/// Send interrupt to master or slave servers.

void TProof::Interrupt(EUrgent type, ESlaves list)
{
   if (!IsValid()) return;

   TList *slaves = 0;
   if (list == kAll)       slaves = fSlaves;
   if (list == kActive)    slaves = fActiveSlaves;
   if (list == kUnique)    slaves = fUniqueSlaves;
   if (list == kAllUnique) slaves = fAllUniqueSlaves;

   if (slaves->GetSize() == 0) return;

   TSlave *sl;
   TIter   next(slaves);

   while ((sl = (TSlave *)next())) {
      if (sl->IsValid()) {

         // Ask slave to progate the interrupt request
         sl->Interrupt((Int_t)type);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns number of slaves active in parallel mode. Returns 0 in case
/// there are no active slaves. Returns -1 in case of error.

Int_t TProof::GetParallel() const
{
   if (!IsValid()) return -1;

   // iterate over active slaves and return total number of slaves
   TIter nextSlave(GetListOfActiveSlaves());
   Int_t nparallel = 0;
   while (TSlave* sl = dynamic_cast<TSlave*>(nextSlave()))
      if (sl->GetParallel() >= 0)
         nparallel += sl->GetParallel();

   return nparallel;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns list of TSlaveInfo's. In case of error return 0.

TList *TProof::GetListOfSlaveInfos()
{
   if (!IsValid()) return 0;

   if (fSlaveInfo == 0) {
      fSlaveInfo = new TSortedList(kSortDescending);
      fSlaveInfo->SetOwner();
   } else {
      fSlaveInfo->Delete();
   }

   TList masters;
   TIter next(GetListOfSlaves());
   TSlave *slave;

   while ((slave = (TSlave *) next()) != 0) {
      if (slave->GetSlaveType() == TSlave::kSlave) {
         const char *name = IsLite() ? gSystem->HostName() : slave->GetName();
         TSlaveInfo *slaveinfo = new TSlaveInfo(slave->GetOrdinal(),
                                                name,
                                                slave->GetPerfIdx());
         fSlaveInfo->Add(slaveinfo);

         TIter nextactive(GetListOfActiveSlaves());
         TSlave *activeslave;
         while ((activeslave = (TSlave *) nextactive())) {
            if (TString(slaveinfo->GetOrdinal()) == activeslave->GetOrdinal()) {
               slaveinfo->SetStatus(TSlaveInfo::kActive);
               break;
            }
         }

         TIter nextbad(GetListOfBadSlaves());
         TSlave *badslave;
         while ((badslave = (TSlave *) nextbad())) {
            if (TString(slaveinfo->GetOrdinal()) == badslave->GetOrdinal()) {
               slaveinfo->SetStatus(TSlaveInfo::kBad);
               break;
            }
         }
         // Get system info if supported
         if (slave->IsValid()) {
            if (slave->GetSocket()->Send(kPROOF_GETSLAVEINFO) == -1)
               MarkBad(slave, "could not send kPROOF_GETSLAVEINFO message");
            else
               masters.Add(slave);
         }

      } else if (slave->GetSlaveType() == TSlave::kMaster) {
         if (slave->IsValid()) {
            if (slave->GetSocket()->Send(kPROOF_GETSLAVEINFO) == -1)
               MarkBad(slave, "could not send kPROOF_GETSLAVEINFO message");
            else
               masters.Add(slave);
         }
      } else {
         Error("GetSlaveInfo", "TSlave is neither Master nor Slave");
         R__ASSERT(0);
      }
   }
   if (masters.GetSize() > 0) Collect(&masters);

   return fSlaveInfo;
}

////////////////////////////////////////////////////////////////////////////////
/// Activate slave server list.

void TProof::Activate(TList *slaves)
{
   TMonitor *mon = fAllMonitor;
   mon->DeActivateAll();

   slaves = !slaves ? fActiveSlaves : slaves;

   TIter next(slaves);
   TSlave *sl;
   while ((sl = (TSlave*) next())) {
      if (sl->IsValid())
         mon->Activate(sl->GetSocket());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Activate (on == TRUE) or deactivate (on == FALSE) all sockets
/// monitored by 'mon'.

void TProof::SetMonitor(TMonitor *mon, Bool_t on)
{
   TMonitor *m = (mon) ? mon : fCurrentMonitor;
   if (m) {
      if (on)
         m->ActivateAll();
      else
         m->DeActivateAll();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast the group priority to all workers in the specified list. Returns
/// the number of workers the message was successfully sent to.
/// Returns -1 in case of error.

Int_t TProof::BroadcastGroupPriority(const char *grp, Int_t priority, TList *workers)
{
   if (!IsValid()) return -1;

   if (workers->GetSize() == 0) return 0;

   int   nsent = 0;
   TIter next(workers);

   TSlave *wrk;
   while ((wrk = (TSlave *)next())) {
      if (wrk->IsValid()) {
         if (wrk->SendGroupPriority(grp, priority) == -1)
            MarkBad(wrk, "could not send group priority");
         else
            nsent++;
      }
   }

   return nsent;
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast the group priority to all workers in the specified list. Returns
/// the number of workers the message was successfully sent to.
/// Returns -1 in case of error.

Int_t TProof::BroadcastGroupPriority(const char *grp, Int_t priority, ESlaves list)
{
   TList *workers = 0;
   if (list == kAll)       workers = fSlaves;
   if (list == kActive)    workers = fActiveSlaves;
   if (list == kUnique)    workers = fUniqueSlaves;
   if (list == kAllUnique) workers = fAllUniqueSlaves;

   return BroadcastGroupPriority(grp, priority, workers);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the merge progress notificator

void TProof::ResetMergePrg()
{
   fMergePrg.Reset(fActiveSlaves->GetSize());
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast a message to all slaves in the specified list. Returns
/// the number of slaves the message was successfully sent to.
/// Returns -1 in case of error.

Int_t TProof::Broadcast(const TMessage &mess, TList *slaves)
{
   if (!IsValid()) return -1;

   if (!slaves || slaves->GetSize() == 0) return 0;

   int   nsent = 0;
   TIter next(slaves);

   TSlave *sl;
   while ((sl = (TSlave *)next())) {
      if (sl->IsValid()) {
         if (sl->GetSocket()->Send(mess) == -1)
            MarkBad(sl, "could not broadcast request");
         else
            nsent++;
      }
   }

   return nsent;
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast a message to all slaves in the specified list (either
/// all slaves or only the active slaves). Returns the number of slaves
/// the message was successfully sent to. Returns -1 in case of error.

Int_t TProof::Broadcast(const TMessage &mess, ESlaves list)
{
   TList *slaves = 0;
   if (list == kAll)       slaves = fSlaves;
   if (list == kActive)    slaves = fActiveSlaves;
   if (list == kUnique)    slaves = fUniqueSlaves;
   if (list == kAllUnique) slaves = fAllUniqueSlaves;

   return Broadcast(mess, slaves);
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast a character string buffer to all slaves in the specified
/// list. Use kind to set the TMessage what field. Returns the number of
/// slaves the message was sent to. Returns -1 in case of error.

Int_t TProof::Broadcast(const char *str, Int_t kind, TList *slaves)
{
   TMessage mess(kind);
   if (str) mess.WriteString(str);
   return Broadcast(mess, slaves);
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast a character string buffer to all slaves in the specified
/// list (either all slaves or only the active slaves). Use kind to
/// set the TMessage what field. Returns the number of slaves the message
/// was sent to. Returns -1 in case of error.

Int_t TProof::Broadcast(const char *str, Int_t kind, ESlaves list)
{
   TMessage mess(kind);
   if (str) mess.WriteString(str);
   return Broadcast(mess, list);
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast an object to all slaves in the specified list. Use kind to
/// set the TMEssage what field. Returns the number of slaves the message
/// was sent to. Returns -1 in case of error.

Int_t TProof::BroadcastObject(const TObject *obj, Int_t kind, TList *slaves)
{
   TMessage mess(kind);
   mess.WriteObject(obj);
   return Broadcast(mess, slaves);
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast an object to all slaves in the specified list. Use kind to
/// set the TMEssage what field. Returns the number of slaves the message
/// was sent to. Returns -1 in case of error.

Int_t TProof::BroadcastObject(const TObject *obj, Int_t kind, ESlaves list)
{
   TMessage mess(kind);
   mess.WriteObject(obj);
   return Broadcast(mess, list);
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast a raw buffer of specified length to all slaves in the
/// specified list. Returns the number of slaves the buffer was sent to.
/// Returns -1 in case of error.

Int_t TProof::BroadcastRaw(const void *buffer, Int_t length, TList *slaves)
{
   if (!IsValid()) return -1;

   if (slaves->GetSize() == 0) return 0;

   int   nsent = 0;
   TIter next(slaves);

   TSlave *sl;
   while ((sl = (TSlave *)next())) {
      if (sl->IsValid()) {
         if (sl->GetSocket()->SendRaw(buffer, length) == -1)
            MarkBad(sl, "could not send broadcast-raw request");
         else
            nsent++;
      }
   }

   return nsent;
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast a raw buffer of specified length to all slaves in the
/// specified list. Returns the number of slaves the buffer was sent to.
/// Returns -1 in case of error.

Int_t TProof::BroadcastRaw(const void *buffer, Int_t length, ESlaves list)
{
   TList *slaves = 0;
   if (list == kAll)       slaves = fSlaves;
   if (list == kActive)    slaves = fActiveSlaves;
   if (list == kUnique)    slaves = fUniqueSlaves;
   if (list == kAllUnique) slaves = fAllUniqueSlaves;

   return BroadcastRaw(buffer, length, slaves);
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast file to all workers in the specified list. Returns the number of workers
/// the buffer was sent to.
/// Returns -1 in case of error.

Int_t TProof::BroadcastFile(const char *file, Int_t opt, const char *rfile, TList *wrks)
{
   if (!IsValid()) return -1;

   if (wrks->GetSize() == 0) return 0;

   int   nsent = 0;
   TIter next(wrks);

   TSlave *wrk;
   while ((wrk = (TSlave *)next())) {
      if (wrk->IsValid()) {
         if (SendFile(file, opt, rfile, wrk) < 0)
            Error("BroadcastFile",
                  "problems sending file to worker %s (%s)",
                  wrk->GetOrdinal(), wrk->GetName());
         else
            nsent++;
      }
   }

   return nsent;
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast file to all workers in the specified list. Returns the number of workers
/// the buffer was sent to.
/// Returns -1 in case of error.

Int_t TProof::BroadcastFile(const char *file, Int_t opt, const char *rfile, ESlaves list)
{
   TList *wrks = 0;
   if (list == kAll)       wrks = fSlaves;
   if (list == kActive)    wrks = fActiveSlaves;
   if (list == kUnique)    wrks = fUniqueSlaves;
   if (list == kAllUnique) wrks = fAllUniqueSlaves;

   return BroadcastFile(file, opt, rfile, wrks);
}

////////////////////////////////////////////////////////////////////////////////
/// Release the used monitor to be used, making sure to delete newly created
/// monitors.

void TProof::ReleaseMonitor(TMonitor *mon)
{
   if (mon && (mon != fAllMonitor) && (mon != fActiveMonitor)
           && (mon != fUniqueMonitor) && (mon != fAllUniqueMonitor)) {
      delete mon;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Collect responses from slave sl. Returns the number of slaves that
/// responded (=1).
/// If timeout >= 0, wait at most timeout seconds (timeout = -1 by default,
/// which means wait forever).
/// If defined (>= 0) endtype is the message that stops this collection.

Int_t TProof::Collect(const TSlave *sl, Long_t timeout, Int_t endtype, Bool_t deactonfail)
{
   Int_t rc = 0;

   TMonitor *mon = 0;
   if (!sl->IsValid()) return 0;

   if (fCurrentMonitor == fAllMonitor) {
      mon = new TMonitor;
   } else {
      mon = fAllMonitor;
      mon->DeActivateAll();
   }
   mon->Activate(sl->GetSocket());

   rc = Collect(mon, timeout, endtype, deactonfail);
   ReleaseMonitor(mon);
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Collect responses from the slave servers. Returns the number of slaves
/// that responded.
/// If timeout >= 0, wait at most timeout seconds (timeout = -1 by default,
/// which means wait forever).
/// If defined (>= 0) endtype is the message that stops this collection.

Int_t TProof::Collect(TList *slaves, Long_t timeout, Int_t endtype, Bool_t deactonfail)
{
   Int_t rc = 0;

   TMonitor *mon = 0;

   if (fCurrentMonitor == fAllMonitor) {
      mon = new TMonitor;
   } else {
      mon = fAllMonitor;
      mon->DeActivateAll();
   }
   TIter next(slaves);
   TSlave *sl;
   while ((sl = (TSlave*) next())) {
      if (sl->IsValid())
         mon->Activate(sl->GetSocket());
   }

   rc = Collect(mon, timeout, endtype, deactonfail);
   ReleaseMonitor(mon);
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Collect responses from the slave servers. Returns the number of slaves
/// that responded.
/// If timeout >= 0, wait at most timeout seconds (timeout = -1 by default,
/// which means wait forever).
/// If defined (>= 0) endtype is the message that stops this collection.

Int_t TProof::Collect(ESlaves list, Long_t timeout, Int_t endtype, Bool_t deactonfail)
{
   Int_t rc = 0;
   TMonitor *mon = 0;

   if (list == kAll)       mon = fAllMonitor;
   if (list == kActive)    mon = fActiveMonitor;
   if (list == kUnique)    mon = fUniqueMonitor;
   if (list == kAllUnique) mon = fAllUniqueMonitor;
   if (fCurrentMonitor == mon) {
      // Get a copy
      mon = new TMonitor(*mon);
   }
   mon->ActivateAll();

   rc = Collect(mon, timeout, endtype, deactonfail);
   ReleaseMonitor(mon);
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Collect responses from the slave servers. Returns the number of messages
/// received. Can be 0 if there are no active slaves.
/// If timeout >= 0, wait at most timeout seconds (timeout = -1 by default,
/// which means wait forever).
/// If defined (>= 0) endtype is the message that stops this collection.
/// Collect also stops its execution from time to time to check for new
/// workers in Dynamic Startup mode.

Int_t TProof::Collect(TMonitor *mon, Long_t timeout, Int_t endtype, Bool_t deactonfail)
{
   Int_t collectId = gRandom->Integer(9999);

   PDB(kCollect, 3)
      Info("Collect", ">>>>>> Entering collect responses #%04d", collectId);

   // Reset the status flag and clear the messages in the list, if any
   fStatus = 0;
   fRecvMessages->Clear();

   Long_t actto = (Long_t)(gEnv->GetValue("Proof.SocketActivityTimeout", -1) * 1000);

   if (!mon->GetActive(actto)) return 0;

   DeActivateAsyncInput();

   // Used by external code to know what we are monitoring
   TMonitor *savedMonitor = 0;
   if (fCurrentMonitor) {
      savedMonitor = fCurrentMonitor;
      fCurrentMonitor = mon;
   } else {
      fCurrentMonitor = mon;
      fBytesRead = 0;
      fRealTime  = 0.0;
      fCpuTime   = 0.0;
   }

   // We want messages on the main window during synchronous collection,
   // but we save the present status to restore it at the end
   Bool_t saveRedirLog = fRedirLog;
   if (!IsIdle() && !IsSync())
      fRedirLog = kFALSE;

   int cnt = 0, rc = 0;

   // Timeout counter
   Long_t nto = timeout;
   PDB(kCollect, 2)
      Info("Collect","#%04d: active: %d", collectId, mon->GetActive());

   // On clients, handle Ctrl-C during collection
   if (fIntHandler)
      fIntHandler->Add();

   // Sockets w/o activity during the last 'sto' millisecs are deactivated
   Int_t nact = 0;
   Long_t sto = -1;
   Int_t nsto = 60;
   Int_t pollint = gEnv->GetValue("Proof.DynamicStartupPollInt", (Int_t) kPROOF_DynWrkPollInt_s);
   mon->ResetInterrupt();
   while ((nact = mon->GetActive(sto)) && (nto < 0 || nto > 0)) {

      // Dump last waiting sockets, if in debug mode
      PDB(kCollect, 2) {
         if (nact < 4) {
            TList *al = mon->GetListOfActives();
            if (al && al->GetSize() > 0) {
               Info("Collect"," %d node(s) still active:", al->GetSize());
               TIter nxs(al);
               TSocket *xs = 0;
               while ((xs = (TSocket *)nxs())) {
                  TSlave *wrk = FindSlave(xs);
                  if (wrk)
                     Info("Collect","   %s (%s)", wrk->GetName(), wrk->GetOrdinal());
                  else
                     Info("Collect","   %p: %s:%d", xs, xs->GetInetAddress().GetHostName(),
                                                        xs->GetInetAddress().GetPort());
               }
            }
            delete al;
         }
      }

      // Preemptive poll for new workers on the master only in Dynamic Mode and only
      // during processing (TODO: should work on Top Master only)
      if (TestBit(TProof::kIsMaster) && !IsIdle() && fDynamicStartup && !fIsPollingWorkers &&
         ((fLastPollWorkers_s == -1) || (time(0)-fLastPollWorkers_s >= pollint))) {
         fIsPollingWorkers = kTRUE;
         if (PollForNewWorkers() > 0) DeActivateAsyncInput();
         fLastPollWorkers_s = time(0);
         fIsPollingWorkers = kFALSE;
         PDB(kCollect, 1)
            Info("Collect","#%04d: now active: %d", collectId, mon->GetActive());
      }

      // Wait for a ready socket
      PDB(kCollect, 3)
         Info("Collect", "Will invoke Select() #%04d", collectId);
      TSocket *s = mon->Select(1000);

      if (s && s != (TSocket *)(-1)) {
         // Get and analyse the info it did receive
         rc = CollectInputFrom(s, endtype, deactonfail);
         if (rc  == 1 || (rc == 2 && !savedMonitor)) {
            // Deactivate it if we are done with it
            mon->DeActivate(s);
            TList *al = mon->GetListOfActives();
            PDB(kCollect, 2)
               Info("Collect","#%04d: deactivating %p (active: %d, %p)", collectId,
                              s, mon->GetActive(),
                              al->First());
            delete al;
         } else if (rc == 2) {
            // This end message was for the saved monitor
            // Deactivate it if we are done with it
            if (savedMonitor) {
               savedMonitor->DeActivate(s);
               TList *al = mon->GetListOfActives();
               PDB(kCollect, 2)
                  Info("Collect","save monitor: deactivating %p (active: %d, %p)",
                                 s, savedMonitor->GetActive(),
                                 al->First());
               delete al;
            }
         }

         // Update counter (if no error occured)
         if (rc >= 0)
            cnt++;
      } else {
         // If not timed-out, exit if not stopped or not aborted
         // (player exits status is finished in such a case); otherwise,
         // we still need to collect the partial output info
         if (!s)
            if (fPlayer && (fPlayer->GetExitStatus() == TVirtualProofPlayer::kFinished))
               mon->DeActivateAll();
         // Decrease the timeout counter if requested
         if (s == (TSocket *)(-1) && nto > 0)
            nto--;
      }

      // Check if there are workers with ready output to be sent and ask the first to send it
      if (IsMaster() && fWrksOutputReady && fWrksOutputReady->GetSize() > 0) {
         // Maximum number of concurrent sendings
         Int_t mxws = gEnv->GetValue("Proof.ControlSendOutput", 1);
         if (TProof::GetParameter(fPlayer->GetInputList(), "PROOF_ControlSendOutput", mxws) != 0)
            mxws = gEnv->GetValue("Proof.ControlSendOutput", 1);
         TIter nxwr(fWrksOutputReady);
         TSlave *wrk = 0;
         while (mxws && (wrk = (TSlave *) nxwr())) {
            if (!wrk->TestBit(TSlave::kOutputRequested)) {
               // Ask worker for output
               TMessage sendoutput(kPROOF_SENDOUTPUT);
               PDB(kCollect, 2)
                  Info("Collect", "worker %s was asked to send its output to master",
                                                wrk->GetOrdinal());
               if (wrk->GetSocket()->Send(sendoutput) != 1) {
                  wrk->SetBit(TSlave::kOutputRequested);
                  mxws--;
               }
            } else {
               // Count
               mxws--;
            }
         }
      }

      // Check if we need to check the socket activity (we do it every 10 cycles ~ 10 sec)
      sto = -1;
      if (--nsto <= 0) {
         sto = (Long_t) actto;
         nsto = 60;
      }

   } // end loop over active monitors

   // If timed-out, deactivate the remaining sockets
   if (nto == 0) {
      TList *al = mon->GetListOfActives();
      if (al && al->GetSize() > 0) {
         // Notify the name of those which did timeout
         Info("Collect"," %d node(s) went in timeout:", al->GetSize());
         TIter nxs(al);
         TSocket *xs = 0;
         while ((xs = (TSocket *)nxs())) {
            TSlave *wrk = FindSlave(xs);
            if (wrk)
               Info("Collect","   %s", wrk->GetName());
            else
               Info("Collect","   %p: %s:%d", xs, xs->GetInetAddress().GetHostName(),
                                                  xs->GetInetAddress().GetPort());
         }
      }
      delete al;
      mon->DeActivateAll();
   }

   // Deactivate Ctrl-C special handler
   if (fIntHandler)
      fIntHandler->Remove();

   // make sure group view is up to date
   SendGroupView();

   // Restore redirection setting
   fRedirLog = saveRedirLog;

   // Restore the monitor
   fCurrentMonitor = savedMonitor;

   ActivateAsyncInput();

   PDB(kCollect, 3)
      Info("Collect", "<<<<<< Exiting collect responses #%04d", collectId);

   return cnt;
}

////////////////////////////////////////////////////////////////////////////////
/// Asks the PROOF Serv for new workers in Dynamic Startup mode and activates
/// them. Returns the number of new workers found, or <0 on errors.

Int_t TProof::PollForNewWorkers()
{
   // Requests for worker updates
   Int_t dummy = 0;
   TList *reqWorkers = new TList();
   reqWorkers->SetOwner(kFALSE);

   if (!TestBit(TProof::kIsMaster)) {
      Error("PollForNewWorkers", "Can't invoke: not on a master -- should not happen!");
      return -1;
   }
   if (!gProofServ) {
      Error("PollForNewWorkers", "No ProofServ available -- should not happen!");
      return -1;
   }

   gProofServ->GetWorkers(reqWorkers, dummy, kTRUE);  // last 2 are dummy

   // List of new workers only (TProofNodeInfo)
   TList *newWorkers = new TList();
   newWorkers->SetOwner(kTRUE);

   TIter next(reqWorkers);
   TProofNodeInfo *ni;
   TString fullOrd;
   while (( ni = dynamic_cast<TProofNodeInfo *>(next()) )) {

      // Form the full ordinal
      fullOrd.Form("%s.%s", gProofServ->GetOrdinal(), ni->GetOrdinal().Data());

      TIter nextInner(fSlaves);
      TSlave *sl;
      Bool_t found = kFALSE;
      while (( sl = dynamic_cast<TSlave *>(nextInner()) )) {
         if ( strcmp(sl->GetOrdinal(), fullOrd.Data()) == 0 ) {
            found = kTRUE;
            break;
         }
      }

      if (found) delete ni;
      else {
         newWorkers->Add(ni);
         PDB(kGlobal, 1)
            Info("PollForNewWorkers", "New worker found: %s:%s",
               ni->GetNodeName().Data(), fullOrd.Data());
      }
   }

   delete reqWorkers;  // not owner

   Int_t nNewWorkers = newWorkers->GetEntries();

   // Add the new workers
   if (nNewWorkers > 0) {
      PDB(kGlobal, 1)
         Info("PollForNewWorkers", "Requesting to add %d new worker(s)", newWorkers->GetEntries());
      Int_t rv = AddWorkers(newWorkers);
      if (rv < 0) {
         Error("PollForNewWorkers", "Call to AddWorkers() failed (got %d < 0)", rv);
         return -1;
      }
      // Don't delete newWorkers: AddWorkers() will do that
   }
   else {
      PDB(kGlobal, 2)
         Info("PollForNewWorkers", "No new worker found");
      delete newWorkers;
   }

   return nNewWorkers;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove links to objects in list 'ol' from gDirectory

void TProof::CleanGDirectory(TList *ol)
{
   if (ol) {
      TIter nxo(ol);
      TObject *o = 0;
      while ((o = nxo()))
         gDirectory->RecursiveRemove(o);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Collect and analyze available input from socket s.
/// Returns 0 on success, -1 if any failure occurs.

Int_t TProof::CollectInputFrom(TSocket *s, Int_t endtype, Bool_t deactonfail)
{
   TMessage *mess;

   Int_t recvrc = 0;
   if ((recvrc = s->Recv(mess)) < 0) {
      PDB(kCollect,2)
         Info("CollectInputFrom","%p: got %d from Recv()", s, recvrc);
      Bool_t bad = kTRUE;
      if (recvrc == -5) {
         // Broken connection: try reconnection
         if (fCurrentMonitor) fCurrentMonitor->Remove(s);
         if (s->Reconnect() == 0) {
            if (fCurrentMonitor) fCurrentMonitor->Add(s);
            bad = kFALSE;
         }
      }
      if (bad)
         MarkBad(s, "problems receiving a message in TProof::CollectInputFrom(...)");
      // Ignore this wake up
      return -1;
   }
   if (!mess) {
      // we get here in case the remote server died
      MarkBad(s, "undefined message in TProof::CollectInputFrom(...)");
      return -1;
   }
   Int_t rc = 0;

   Int_t what = mess->What();
   TSlave *sl = FindSlave(s);
   rc = HandleInputMessage(sl, mess, deactonfail);
   if (rc == 1 && (endtype >= 0) && (what != endtype))
      // This message was for the base monitor in recursive case
      rc = 2;

   // We are done successfully
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Analyze the received message.
/// Returns 0 on success (1 if this the last message from this socket), -1 if
/// any failure occurs.

Int_t TProof::HandleInputMessage(TSlave *sl, TMessage *mess, Bool_t deactonfail)
{
   char str[512];
   TObject *obj;
   Int_t rc = 0;

   if (!mess || !sl) {
      Warning("HandleInputMessage", "given an empty message or undefined worker");
      return -1;
   }
   Bool_t delete_mess = kTRUE;
   TSocket *s = sl->GetSocket();
   if (!s) {
      Warning("HandleInputMessage", "worker socket is undefined");
      return -1;
   }

   // The message type
   Int_t what = mess->What();

   PDB(kCollect,3)
      Info("HandleInputMessage", "got type %d from '%s'", what, sl->GetOrdinal());

   switch (what) {

      case kMESS_OK:
         // Add the message to the list
         fRecvMessages->Add(mess);
         delete_mess = kFALSE;
         break;

      case kMESS_OBJECT:
         if (fPlayer) fPlayer->HandleRecvHisto(mess);
         break;

      case kPROOF_FATAL:
         {  TString msg;
            if ((mess->BufferSize() > mess->Length()))
               (*mess) >> msg;
            if (msg.IsNull()) {
               MarkBad(s, "received kPROOF_FATAL");
            } else {
               MarkBad(s, msg);
            }
         }
         if (fProgressDialogStarted) {
            // Finalize the progress dialog
            Emit("StopProcess(Bool_t)", kTRUE);
         }
         break;

      case kPROOF_STOP:
         // Stop collection from this worker
         Info("HandleInputMessage", "received kPROOF_STOP from %s: disabling any further collection this worker",
              sl->GetOrdinal());
         rc = 1;
         break;

      case kPROOF_GETTREEHEADER:
         // Add the message to the list
         fRecvMessages->Add(mess);
         delete_mess = kFALSE;
         rc = 1;
         break;

      case kPROOF_TOUCH:
         // send a request for touching the remote admin file
         {
            sl->Touch();
         }
         break;

      case kPROOF_GETOBJECT:
         // send slave object it asks for
         mess->ReadString(str, sizeof(str));
         obj = gDirectory->Get(str);
         if (obj)
            s->SendObject(obj);
         else
            s->Send(kMESS_NOTOK);
         break;

      case kPROOF_GETPACKET:
         {
            PDB(kGlobal,2)
               Info("HandleInputMessage","%s: kPROOF_GETPACKET", sl->GetOrdinal());
            TDSetElement *elem = 0;
            elem = fPlayer ? fPlayer->GetNextPacket(sl, mess) : 0;

            if (elem != (TDSetElement*) -1) {
               TMessage answ(kPROOF_GETPACKET);
               answ << elem;
               s->Send(answ);

               while (fWaitingSlaves != 0 && fWaitingSlaves->GetSize()) {
                  TPair *p = (TPair*) fWaitingSlaves->First();
                  s = (TSocket*) p->Key();
                  TMessage *m = (TMessage*) p->Value();

                  elem = fPlayer ? fPlayer->GetNextPacket(sl, m) : 0;
                  if (elem != (TDSetElement*) -1) {
                     TMessage a(kPROOF_GETPACKET);
                     a << elem;
                     s->Send(a);
                     // remove has to happen via Links because TPair does not have
                     // a Compare() function and therefore RemoveFirst() and
                     // Remove(TObject*) do not work
                     fWaitingSlaves->Remove(fWaitingSlaves->FirstLink());
                     delete p;
                     delete m;
                  } else {
                     break;
                  }
               }
            } else {
               if (fWaitingSlaves == 0) fWaitingSlaves = new TList;
               fWaitingSlaves->Add(new TPair(s, mess));
               delete_mess = kFALSE;
            }
         }
         break;

      case kPROOF_LOGFILE:
         {
            Int_t size;
            (*mess) >> size;
            PDB(kGlobal,2)
               Info("HandleInputMessage","%s: kPROOF_LOGFILE: size: %d", sl->GetOrdinal(), size);
            RecvLogFile(s, size);
         }
         break;

      case kPROOF_LOGDONE:
         (*mess) >> sl->fStatus >> sl->fParallel;
         PDB(kCollect,2)
            Info("HandleInputMessage","%s: kPROOF_LOGDONE: status %d  parallel %d",
                 sl->GetOrdinal(), sl->fStatus, sl->fParallel);
         if (sl->fStatus != 0) {
            // Return last nonzero status
            fStatus = sl->fStatus;
            // Deactivate the worker, if required
            if (deactonfail) DeactivateWorker(sl->fOrdinal);
         }
         // Remove from the workers-ready list
         if (fWrksOutputReady && fWrksOutputReady->FindObject(sl)) {
            sl->ResetBit(TSlave::kOutputRequested);
            fWrksOutputReady->Remove(sl);
         }
         rc = 1;
         break;

      case kPROOF_GETSTATS:
         {
            (*mess) >> sl->fBytesRead >> sl->fRealTime >> sl->fCpuTime
                  >> sl->fWorkDir >> sl->fProofWorkDir;
            PDB(kCollect,2)
               Info("HandleInputMessage", "kPROOF_GETSTATS: %s", sl->fWorkDir.Data());
            TString img;
            if ((mess->BufferSize() > mess->Length()))
               (*mess) >> img;
            // Set image
            if (img.IsNull()) {
               if (sl->fImage.IsNull())
                  sl->fImage.Form("%s:%s", TUrl(sl->fName).GetHostFQDN(),
                                           sl->fProofWorkDir.Data());
            } else {
               sl->fImage = img;
            }
            PDB(kGlobal,2)
               Info("HandleInputMessage",
                        "kPROOF_GETSTATS:%s image: %s", sl->GetOrdinal(), sl->GetImage());

            fBytesRead += sl->fBytesRead;
            fRealTime  += sl->fRealTime;
            fCpuTime   += sl->fCpuTime;
            rc = 1;
         }
         break;

      case kPROOF_GETPARALLEL:
         {
            Bool_t async = kFALSE;
            (*mess) >> sl->fParallel;
            if ((mess->BufferSize() > mess->Length()))
               (*mess) >> async;
            rc = (async) ? 0 : 1;
         }
         break;

      case kPROOF_CHECKFILE:
         {  // New servers (>= 5.22) send the status
            if ((mess->BufferSize() > mess->Length())) {
               (*mess) >> fCheckFileStatus;
            } else {
               // Form old servers this meant success (failure was signaled with the
               // dangerous kPROOF_FATAL)
               fCheckFileStatus = 1;
            }
            rc = 1;
         }
         break;

      case kPROOF_SENDFILE:
         {  // New server: signals ending of sendfile operation
            rc = 1;
         }
         break;

      case kPROOF_PACKAGE_LIST:
         {
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_PACKAGE_LIST: enter");
            Int_t type = 0;
            (*mess) >> type;
            switch (type) {
            case TProof::kListEnabledPackages:
               SafeDelete(fEnabledPackages);
               fEnabledPackages = (TList *) mess->ReadObject(TList::Class());
               if (fEnabledPackages) {
                  fEnabledPackages->SetOwner();
               } else {
                  Error("HandleInputMessage",
                        "kPROOF_PACKAGE_LIST: kListEnabledPackages: TList not found in message!");
               }
               break;
            case TProof::kListPackages:
               SafeDelete(fAvailablePackages);
               fAvailablePackages = (TList *) mess->ReadObject(TList::Class());
               if (fAvailablePackages) {
                  fAvailablePackages->SetOwner();
               } else {
                  Error("HandleInputMessage",
                        "kPROOF_PACKAGE_LIST: kListPackages: TList not found in message!");
               }
               break;
            default:
               Error("HandleInputMessage", "kPROOF_PACKAGE_LIST: unknown type: %d", type);
            }
         }
         break;

      case kPROOF_SENDOUTPUT:
         {
            // We start measuring the merging time
            fPlayer->SetMerging();

            // Worker is ready to send output: make sure the relevant bit is reset
            sl->ResetBit(TSlave::kOutputRequested);
            PDB(kGlobal,2)
               Info("HandleInputMessage","kPROOF_SENDOUTPUT: enter (%s)", sl->GetOrdinal());
            // Create the list if not yet done
            if (!fWrksOutputReady) {
               fWrksOutputReady = new TList;
               fWrksOutputReady->SetOwner(kFALSE);
            }
            fWrksOutputReady->Add(sl);
         }
         break;

      case kPROOF_OUTPUTOBJECT:
         {
            // We start measuring the merging time
            fPlayer->SetMerging();

            PDB(kGlobal,2)
               Info("HandleInputMessage","kPROOF_OUTPUTOBJECT: enter");
            Int_t type = 0;
            const char *prefix = gProofServ ? gProofServ->GetPrefix() : "Lite-0";
            if (!TestBit(TProof::kIsClient) && !fMergersSet && !fFinalizationRunning) {
               Info("HandleInputMessage", "finalization on %s started ...", prefix);
               fFinalizationRunning = kTRUE;
            }

            while ((mess->BufferSize() > mess->Length())) {
               (*mess) >> type;
               // If a query result header, add it to the player list
               if (fPlayer) {
                  if (type == 0) {
                     // Retrieve query result instance (output list not filled)
                     TQueryResult *pq =
                        (TQueryResult *) mess->ReadObject(TQueryResult::Class());
                     if (pq) {
                        // Add query to the result list in TProofPlayer
                        fPlayer->AddQueryResult(pq);
                        fPlayer->SetCurrentQuery(pq);
                        // And clear the output list, as we start merging a new set of results
                        if (fPlayer->GetOutputList())
                           fPlayer->GetOutputList()->Clear();
                        // Add the unique query tag as TNamed object to the input list
                        // so that it is available in TSelectors for monitoring
                        TString qid = TString::Format("%s:%s",pq->GetTitle(),pq->GetName());
                        if (fPlayer->GetInputList()->FindObject("PROOF_QueryTag"))
                           fPlayer->GetInputList()->Remove(fPlayer->GetInputList()->FindObject("PROOF_QueryTag"));
                        fPlayer->AddInput(new TNamed("PROOF_QueryTag", qid.Data()));
                     } else {
                        Warning("HandleInputMessage","kPROOF_OUTPUTOBJECT: query result missing");
                     }
                  } else if (type > 0) {
                     // Read object
                     TObject *o = mess->ReadObject(TObject::Class());
                     // Increment counter on the client side
                     fMergePrg.IncreaseIdx();
                     TString msg;
                     Bool_t changed = kFALSE;
                     msg.Form("%s: merging output objects ... %s", prefix, fMergePrg.Export(changed));
                     if (gProofServ) {
                        gProofServ->SendAsynMessage(msg.Data(), kFALSE);
                     } else if (IsTty() || changed) {
                        fprintf(stderr, "%s\r", msg.Data());
                     }
                     // Add or merge it
                     if ((fPlayer->AddOutputObject(o) == 1)) {
                        // Remove the object if it has been merged
                        SafeDelete(o);
                     }
                     if (type > 1) {
                        // Update the merger progress info
                        fMergePrg.DecreaseNWrks();
                        if (TestBit(TProof::kIsClient) && !IsLite()) {
                           // In PROOFLite this has to be done once only in TProofLite::Process
                           TQueryResult *pq = fPlayer->GetCurrentQuery();
                           if (pq) {
                              pq->SetOutputList(fPlayer->GetOutputList(), kFALSE);
                              // Add input objects (do not override remote settings, if any)
                              TObject *xo = 0;
                              TIter nxin(fPlayer->GetInputList());
                              // Servers prior to 5.28/00 do not create the input list in the TQueryResult
                              if (!pq->GetInputList()) pq->SetInputList(new TList());
                              while ((xo = nxin()))
                                 if (!pq->GetInputList()->FindObject(xo->GetName()))
                                    pq->AddInput(xo->Clone());
                              // If the last object, notify the GUI that the result arrived
                              QueryResultReady(TString::Format("%s:%s", pq->GetTitle(), pq->GetName()));
                           }
                           // Processing is over
                           UpdateDialog();
                        }
                     }
                  }
               } else {
                  Warning("HandleInputMessage", "kPROOF_OUTPUTOBJECT: player undefined!");
               }
            }
         }
         break;

      case kPROOF_OUTPUTLIST:
         {
            // We start measuring the merging time

            PDB(kGlobal,2)
               Info("HandleInputMessage","%s: kPROOF_OUTPUTLIST: enter", sl->GetOrdinal());
            TList *out = 0;
            if (fPlayer) {
               fPlayer->SetMerging();
               if (TestBit(TProof::kIsMaster) || fProtocol < 7) {
                  out = (TList *) mess->ReadObject(TList::Class());
               } else {
                  TQueryResult *pq =
                     (TQueryResult *) mess->ReadObject(TQueryResult::Class());
                  if (pq) {
                     // Add query to the result list in TProofPlayer
                     fPlayer->AddQueryResult(pq);
                     fPlayer->SetCurrentQuery(pq);
                     // To avoid accidental cleanups from anywhere else
                     // remove objects from gDirectory and clone the list
                     out = pq->GetOutputList();
                     CleanGDirectory(out);
                     out = (TList *) out->Clone();
                     // Notify the GUI that the result arrived
                     QueryResultReady(TString::Format("%s:%s", pq->GetTitle(), pq->GetName()));
                  } else {
                     PDB(kGlobal,2)
                        Info("HandleInputMessage",
                             "%s: kPROOF_OUTPUTLIST: query result missing", sl->GetOrdinal());
                  }
               }
               if (out) {
                  out->SetOwner();
                  fPlayer->AddOutput(out); // Incorporate the list
                  SafeDelete(out);
               } else {
                  PDB(kGlobal,2)
                     Info("HandleInputMessage",
                          "%s: kPROOF_OUTPUTLIST: outputlist is empty", sl->GetOrdinal());
               }
            } else {
               Warning("HandleInputMessage",
                       "%s: kPROOF_OUTPUTLIST: player undefined!", sl->GetOrdinal());
            }
            // On clients at this point processing is over
            if (TestBit(TProof::kIsClient) && !IsLite())
               UpdateDialog();
         }
         break;

      case kPROOF_QUERYLIST:
         {
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_QUERYLIST: enter");
            (*mess) >> fOtherQueries >> fDrawQueries;
            if (fQueries) {
               fQueries->Delete();
               delete fQueries;
               fQueries = 0;
            }
            fQueries = (TList *) mess->ReadObject(TList::Class());
         }
         break;

      case kPROOF_RETRIEVE:
         {
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_RETRIEVE: enter");
            TQueryResult *pq =
               (TQueryResult *) mess->ReadObject(TQueryResult::Class());
            if (pq && fPlayer) {
               fPlayer->AddQueryResult(pq);
               // Notify the GUI that the result arrived
               QueryResultReady(TString::Format("%s:%s", pq->GetTitle(), pq->GetName()));
            } else {
               PDB(kGlobal,2)
                  Info("HandleInputMessage",
                       "kPROOF_RETRIEVE: query result missing or player undefined");
            }
         }
         break;

      case kPROOF_MAXQUERIES:
         {
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_MAXQUERIES: enter");
            Int_t max = 0;

            (*mess) >> max;
            Printf("Number of queries fully kept remotely: %d", max);
         }
         break;

      case kPROOF_SERVERSTARTED:
         {
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_SERVERSTARTED: enter");

            UInt_t tot = 0, done = 0;
            TString action;
            Bool_t st = kTRUE;

            (*mess) >> action >> tot >> done >> st;

            if (TestBit(TProof::kIsClient)) {
               if (tot) {
                  TString type = (action.Contains("submas")) ? "submasters"
                                                             : "workers";
                  Int_t frac = (Int_t) (done*100.)/tot;
                  char msg[512] = {0};
                  if (frac >= 100) {
                     snprintf(msg, 512, "%s: OK (%d %s)                 \n",
                             action.Data(),tot, type.Data());
                  } else {
                     snprintf(msg, 512, "%s: %d out of %d (%d %%)\r",
                             action.Data(), done, tot, frac);
                  }
                  if (fSync)
                     fprintf(stderr,"%s", msg);
                  else
                     NotifyLogMsg(msg, 0);
               }
               // Notify GUIs
               StartupMessage(action.Data(), st, (Int_t)done, (Int_t)tot);
            } else {

               // Just send the message one level up
               TMessage m(kPROOF_SERVERSTARTED);
               m << action << tot << done << st;
               gProofServ->GetSocket()->Send(m);
            }
         }
         break;

      case kPROOF_DATASET_STATUS:
         {
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_DATASET_STATUS: enter");

            UInt_t tot = 0, done = 0;
            TString action;
            Bool_t st = kTRUE;

            (*mess) >> action >> tot >> done >> st;

            if (TestBit(TProof::kIsClient)) {
               if (tot) {
                  TString type = "files";
                  Int_t frac = (Int_t) (done*100.)/tot;
                  char msg[512] = {0};
                  if (frac >= 100) {
                     snprintf(msg, 512, "%s: OK (%d %s)                 \n",
                             action.Data(),tot, type.Data());
                  } else {
                     snprintf(msg, 512, "%s: %d out of %d (%d %%)\r",
                             action.Data(), done, tot, frac);
                  }
                  if (fSync)
                     fprintf(stderr,"%s", msg);
                  else
                     NotifyLogMsg(msg, 0);
               }
               // Notify GUIs
               DataSetStatus(action.Data(), st, (Int_t)done, (Int_t)tot);
            } else {

               // Just send the message one level up
               TMessage m(kPROOF_DATASET_STATUS);
               m << action << tot << done << st;
               gProofServ->GetSocket()->Send(m);
            }
         }
         break;

      case kPROOF_STARTPROCESS:
         {
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_STARTPROCESS: enter");

            // For Proof-Lite this variable is the number of workers and is set
            // by the player
            if (!IsLite()) {
               fNotIdle = 1;
               fIsWaiting = kFALSE;
            }

            // Redirect the output, if needed
            fRedirLog = (fSync) ? fRedirLog : kTRUE;

            // The signal is used on masters by XrdProofdProtocol to catch
            // the start of processing; on clients it allows to update the
            // progress dialog
            if (!TestBit(TProof::kIsMaster)) {

               // This is the end of preparation
               fQuerySTW.Stop();
               fPrepTime = fQuerySTW.RealTime();
               PDB(kGlobal,2) Info("HandleInputMessage","Preparation time: %f s", fPrepTime);

               TString selec;
               Int_t dsz = -1;
               Long64_t first = -1, nent = -1;
               (*mess) >> selec >> dsz >> first >> nent;
               // Start or reset the progress dialog
               if (!gROOT->IsBatch()) {
                  if (fProgressDialog &&
                      !TestBit(kUsingSessionGui) && TestBit(kUseProgressDialog)) {
                     if (!fProgressDialogStarted) {
                        fProgressDialog->ExecPlugin(5, this,
                                                   selec.Data(), dsz, first, nent);
                        fProgressDialogStarted = kTRUE;
                     } else {
                        ResetProgressDialog(selec, dsz, first, nent);
                     }
                  }
                  ResetBit(kUsingSessionGui);
               }
            }
         }
         break;

      case kPROOF_ENDINIT:
         {
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_ENDINIT: enter");

            if (TestBit(TProof::kIsMaster)) {
               if (fPlayer)
                  fPlayer->SetInitTime();
            }
         }
         break;

      case kPROOF_SETIDLE:
         {
            PDB(kGlobal,2)
               Info("HandleInputMessage","kPROOF_SETIDLE from '%s': enter (%d)", sl->GetOrdinal(), fNotIdle);

            // The session is idle
            if (IsLite()) {
               if (fNotIdle > 0) {
                  fNotIdle--;
                  PDB(kGlobal,2)
                     Info("HandleInputMessage", "%s: got kPROOF_SETIDLE", sl->GetOrdinal());
               } else {
                  Warning("HandleInputMessage",
                          "%s: got kPROOF_SETIDLE but no running workers ! protocol error?",
                          sl->GetOrdinal());
               }
            } else {
               fNotIdle = 0;
               // Check if the query has been enqueued
               if ((mess->BufferSize() > mess->Length()))
                  (*mess) >> fIsWaiting;
            }
         }
         break;

      case kPROOF_QUERYSUBMITTED:
         {
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_QUERYSUBMITTED: enter");

            // We have received the sequential number
            (*mess) >> fSeqNum;
            Bool_t sync = fSync;
            if ((mess->BufferSize() > mess->Length()))
               (*mess) >> sync;
            if (sync !=  fSync && fSync) {
               // The server required to switch to asynchronous mode
               Activate();
               fSync = kFALSE;
            }
            DisableGoAsyn();
            // Check if the query has been enqueued
            fIsWaiting = kTRUE;
            // For Proof-Lite this variable is the number of workers and is set by the player
            if (!IsLite())
               fNotIdle = 1;

            rc = 1;
         }
         break;

      case kPROOF_SESSIONTAG:
         {
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_SESSIONTAG: enter");

            // We have received the unique tag and save it as name of this object
            TString stag;
            (*mess) >> stag;
            SetName(stag);
            // In the TSlave object
            sl->SetSessionTag(stag);
            // Server may have also sent the group
            if ((mess->BufferSize() > mess->Length()))
               (*mess) >> fGroup;
            // Server may have also sent the user
            if ((mess->BufferSize() > mess->Length())) {
               TString usr;
               (*mess) >> usr;
               if (!usr.IsNull()) fUrl.SetUser(usr.Data());
            }
         }
         break;

      case kPROOF_FEEDBACK:
         {
            PDB(kGlobal,2)
               Info("HandleInputMessage","kPROOF_FEEDBACK: enter");
            TList *out = (TList *) mess->ReadObject(TList::Class());
            out->SetOwner();
            if (fPlayer)
               fPlayer->StoreFeedback(sl, out); // Adopts the list
            else
               // Not yet ready: stop collect asap
               rc = 1;
         }
         break;

      case kPROOF_AUTOBIN:
         {
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_AUTOBIN: enter");

            TString name;
            Double_t xmin, xmax, ymin, ymax, zmin, zmax;

            (*mess) >> name >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;

            if (fPlayer) fPlayer->UpdateAutoBin(name,xmin,xmax,ymin,ymax,zmin,zmax);

            TMessage answ(kPROOF_AUTOBIN);

            answ << name << xmin << xmax << ymin << ymax << zmin << zmax;

            s->Send(answ);
         }
         break;

      case kPROOF_PROGRESS:
         {
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_PROGRESS: enter");

            if (GetRemoteProtocol() > 25) {
               // New format
               TProofProgressInfo *pi = 0;
               (*mess) >> pi;
               fPlayer->Progress(sl,pi);
            } else if (GetRemoteProtocol() > 11) {
               Long64_t total, processed, bytesread;
               Float_t initTime, procTime, evtrti, mbrti;
               (*mess) >> total >> processed >> bytesread
                     >> initTime >> procTime
                     >> evtrti >> mbrti;
               if (fPlayer)
                  fPlayer->Progress(sl, total, processed, bytesread,
                                    initTime, procTime, evtrti, mbrti);

            } else {
               // Old format
               Long64_t total, processed;
               (*mess) >> total >> processed;
               if (fPlayer)
                  fPlayer->Progress(sl, total, processed);
            }
         }
         break;

      case kPROOF_STOPPROCESS:
         {
            // This message is sent from a worker that finished processing.
            // We determine whether it was asked to finish by the
            // packetizer or stopped during processing a packet
            // (by TProof::RemoveWorkers() or by an external signal).
            // In the later case call packetizer->MarkBad.
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_STOPPROCESS: enter");

            Long64_t events = 0;
            Bool_t abort = kFALSE;
            TProofProgressStatus *status = 0;

            if ((mess->BufferSize() > mess->Length()) && (fProtocol > 18)) {
               (*mess) >> status >> abort;
            } else if ((mess->BufferSize() > mess->Length()) && (fProtocol > 8)) {
               (*mess) >> events >> abort;
            } else {
               (*mess) >> events;
            }
            if (fPlayer) {
               if (fProtocol > 18) {
                  TList *listOfMissingFiles = 0;
                  if (!(listOfMissingFiles = (TList *)GetOutput("MissingFiles"))) {
                     listOfMissingFiles = new TList();
                     listOfMissingFiles->SetName("MissingFiles");
                     if (fPlayer)
                        fPlayer->AddOutputObject(listOfMissingFiles);
                  }
                  if (fPlayer->GetPacketizer()) {
                     Int_t ret =
                        fPlayer->GetPacketizer()->AddProcessed(sl, status, 0, &listOfMissingFiles);
                     if (ret > 0)
                        fPlayer->GetPacketizer()->MarkBad(sl, status, &listOfMissingFiles);
                     // This object is now owned by the packetizer
                     status = 0;
                  }
                  if (status) fPlayer->AddEventsProcessed(status->GetEntries());
               } else {
                  fPlayer->AddEventsProcessed(events);
               }
            }
            SafeDelete(status);
            if (!TestBit(TProof::kIsMaster))
               Emit("StopProcess(Bool_t)", abort);
            break;
         }

      case kPROOF_SUBMERGER:
         {
            PDB(kGlobal,2) Info("HandleInputMessage", "kPROOF_SUBMERGER: enter");
            HandleSubmerger(mess, sl);
         }
         break;

      case kPROOF_GETSLAVEINFO:
         {
            PDB(kGlobal,2) Info("HandleInputMessage", "kPROOF_GETSLAVEINFO: enter");

            Bool_t active = (GetListOfActiveSlaves()->FindObject(sl) != 0);
            Bool_t bad = (GetListOfBadSlaves()->FindObject(sl) != 0);
            TList* tmpinfo = 0;
            (*mess) >> tmpinfo;
            if (tmpinfo == 0) {
               Error("HandleInputMessage", "kPROOF_GETSLAVEINFO: no list received!");
            } else {
               tmpinfo->SetOwner(kFALSE);
               Int_t nentries = tmpinfo->GetSize();
               for (Int_t i=0; i<nentries; i++) {
                  TSlaveInfo* slinfo =
                     dynamic_cast<TSlaveInfo*>(tmpinfo->At(i));
                  if (slinfo) {
                     // If PROOF-Lite
                     if (IsLite()) slinfo->fHostName = gSystem->HostName();
                     // Check if we have already a instance for this worker
                     TIter nxw(fSlaveInfo);
                     TSlaveInfo *ourwi = 0;
                     while ((ourwi = (TSlaveInfo *)nxw())) {
                        if (!strcmp(ourwi->GetOrdinal(), slinfo->GetOrdinal())) {
                           ourwi->SetSysInfo(slinfo->GetSysInfo());
                           ourwi->fHostName = slinfo->GetName();
                           if (slinfo->GetDataDir() && (strlen(slinfo->GetDataDir()) > 0))
                              ourwi->fDataDir = slinfo->GetDataDir();
                           break;
                        }
                     }
                     if (!ourwi) {
                        fSlaveInfo->Add(slinfo);
                     } else {
                        slinfo = ourwi;
                     }
                     if (slinfo->fStatus != TSlaveInfo::kBad) {
                        if (!active) slinfo->SetStatus(TSlaveInfo::kNotActive);
                        if (bad) slinfo->SetStatus(TSlaveInfo::kBad);
                     }
                     if (sl->GetMsd() && (strlen(sl->GetMsd()) > 0))
                        slinfo->fMsd = sl->GetMsd();
                  }
               }
               delete tmpinfo;
               rc = 1;
            }
         }
         break;

      case kPROOF_VALIDATE_DSET:
         {
            PDB(kGlobal,2)
               Info("HandleInputMessage", "kPROOF_VALIDATE_DSET: enter");
            TDSet* dset = 0;
            (*mess) >> dset;
            if (!fDSet)
               Error("HandleInputMessage", "kPROOF_VALIDATE_DSET: fDSet not set");
            else
               fDSet->Validate(dset);
            delete dset;
         }
         break;

      case kPROOF_DATA_READY:
         {
            PDB(kGlobal,2) Info("HandleInputMessage", "kPROOF_DATA_READY: enter");
            Bool_t dataready = kFALSE;
            Long64_t totalbytes, bytesready;
            (*mess) >> dataready >> totalbytes >> bytesready;
            fTotalBytes += totalbytes;
            fBytesReady += bytesready;
            if (dataready == kFALSE) fDataReady = dataready;
         }
         break;

      case kPROOF_PING:
         // do nothing (ping is already acknowledged)
         break;

      case kPROOF_MESSAGE:
         {
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_MESSAGE: enter");

            // We have received the unique tag and save it as name of this object
            TString msg;
            (*mess) >> msg;
            Bool_t lfeed = kTRUE;
            if ((mess->BufferSize() > mess->Length()))
               (*mess) >> lfeed;

            if (TestBit(TProof::kIsClient)) {

               if (fSync) {
                  // Notify locally
                  fprintf(stderr,"%s%c", msg.Data(), (lfeed ? '\n' : '\r'));
               } else {
                  // Notify locally taking care of redirection, windows logs, ...
                  NotifyLogMsg(msg, (lfeed ? "\n" : "\r"));
               }
            } else {

               // The message is logged for debugging purposes.
               fprintf(stderr,"%s%c", msg.Data(), (lfeed ? '\n' : '\r'));
               if (gProofServ) {
                  // We hide it during normal operations
                  gProofServ->FlushLogFile();

                  // And send the message one level up
                  gProofServ->SendAsynMessage(msg, lfeed);
               }
            }
         }
         break;

      case kPROOF_VERSARCHCOMP:
         {
            TString vac;
            (*mess) >> vac;
            PDB(kGlobal,2) Info("HandleInputMessage","kPROOF_VERSARCHCOMP: %s", vac.Data());
            Int_t from = 0;
            TString vers, archcomp;
            if (vac.Tokenize(vers, from, "|"))
               vac.Tokenize(archcomp, from, "|");
            sl->SetArchCompiler(archcomp);
            vers.ReplaceAll(":","|");
            sl->SetROOTVersion(vers);
         }
         break;

      default:
         {
            Error("HandleInputMessage", "unknown command received from '%s' (what = %d)",
                                        sl->GetOrdinal(), what);
         }
         break;
   }

   // Cleanup
   if (delete_mess)
      delete mess;

   // We are done successfully
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Process a message of type kPROOF_SUBMERGER

void TProof::HandleSubmerger(TMessage *mess, TSlave *sl)
{
   // Message sub-type
   Int_t type = 0;
   (*mess) >> type;
   TSocket *s = sl->GetSocket();

   switch (type) {
      case kOutputSent:
         {
            if (IsEndMaster()) {
               Int_t  merger_id = -1;
               (*mess) >> merger_id;

               PDB(kSubmerger, 2)
                  Info("HandleSubmerger", "kOutputSent: Worker %s:%d:%s had sent its output to merger #%d",
                                          sl->GetName(), sl->GetPort(), sl->GetOrdinal(), merger_id);

               if (!fMergers || fMergers->GetSize() <= merger_id) {
                  Error("HandleSubmerger", "kOutputSize: #%d not in list ", merger_id);
                  break;
               }
               TMergerInfo * mi = (TMergerInfo *) fMergers->At(merger_id);
               mi->SetMergedWorker();
               if (mi->AreAllWorkersMerged()) {
                  mi->Deactivate();
                  if (GetActiveMergersCount() == 0) {
                     fMergers->Clear();
                     delete fMergers;
                     fMergersSet = kFALSE;
                     fMergersCount = -1;
                     fLastAssignedMerger = 0;
                     PDB(kSubmerger, 2) Info("HandleSubmerger", "all mergers removed ... ");
                  }
               }
            } else {
               PDB(kSubmerger, 2) Error("HandleSubmerger","kOutputSent: received not on endmaster!");
            }
         }
         break;

      case kMergerDown:
         {
            Int_t  merger_id = -1;
            (*mess) >> merger_id;

            PDB(kSubmerger, 2) Info("HandleSubmerger", "kMergerDown: #%d ", merger_id);

            if (!fMergers || fMergers->GetSize() <= merger_id) {
               Error("HandleSubmerger", "kMergerDown: #%d not in list ", merger_id);
               break;
            }

            TMergerInfo * mi = (TMergerInfo *) fMergers->At(merger_id);
            if (!mi->IsActive()) {
               break;
            } else {
               mi->Deactivate();
            }

            // Stop the invalid merger in the case it is still listening
            TMessage stop(kPROOF_SUBMERGER);
            stop << Int_t(kStopMerging);
            stop <<  0;
            s->Send(stop);

            // Ask for results from merger (only original results from this node as worker are returned)
            AskForOutput(mi->GetMerger());

            // Ask for results from all workers assigned to this merger
            TIter nxo(mi->GetWorkers());
            TObject * o = 0;
            while ((o = nxo())) {
               AskForOutput((TSlave *)o);
            }
            PDB(kSubmerger, 2) Info("HandleSubmerger", "kMergerDown:%d: exit", merger_id);
         }
         break;

      case kOutputSize:
         {
            if (IsEndMaster()) {
               PDB(kSubmerger, 2)
                  Info("HandleSubmerger", "worker %s reported as finished ", sl->GetOrdinal());

               const char *prefix = gProofServ ? gProofServ->GetPrefix() : "Lite-0";
               if (!fFinalizationRunning) {
                  Info("HandleSubmerger", "finalization on %s started ...", prefix);
                  fFinalizationRunning = kTRUE;
               }

               Int_t  output_size = 0;
               Int_t  merging_port = 0;
               (*mess) >> output_size >> merging_port;

               PDB(kSubmerger, 2) Info("HandleSubmerger",
                                       "kOutputSize: Worker %s:%d:%s reports %d output objects (+ available port %d)",
                                       sl->GetName(), sl->GetPort(), sl->GetOrdinal(), output_size, merging_port);
               TString msg;
               if (!fMergersSet) {

                  Int_t activeWorkers = fCurrentMonitor ? fCurrentMonitor->GetActive() : GetNumberOfActiveSlaves();

                  // First pass - setting number of mergers according to user or dynamically
                  fMergersCount = -1; // No mergers used if not set by user
                  TParameter<Int_t> *mc = dynamic_cast<TParameter<Int_t> *>(GetParameter("PROOF_UseMergers"));
                  if (mc) fMergersCount = mc->GetVal(); // Value set by user
                  TParameter<Int_t> *mh = dynamic_cast<TParameter<Int_t> *>(GetParameter("PROOF_MergersByHost"));
                  if (mh) fMergersByHost = (mh->GetVal() != 0) ? kTRUE : kFALSE; // Assign submergers by hostname

                  // Mergers count specified by user but not valid
                  if (fMergersCount < 0 || (fMergersCount > (activeWorkers/2) )) {
                     msg.Form("%s: Invalid request: cannot start %d mergers for %d workers",
                              prefix, fMergersCount, activeWorkers);
                     if (gProofServ)
                        gProofServ->SendAsynMessage(msg);
                     else
                        Printf("%s",msg.Data());
                     fMergersCount = 0;
                  }
                  // Mergers count will be set dynamically
                  if ((fMergersCount == 0) && (!fMergersByHost)) {
                     if (activeWorkers > 1) {
                        fMergersCount = TMath::Nint(TMath::Sqrt(activeWorkers));
                        if (activeWorkers / fMergersCount < 2)
                           fMergersCount = (Int_t) TMath::Sqrt(activeWorkers);
                     }
                     if (fMergersCount > 1)
                        msg.Form("%s: Number of mergers set dynamically to %d (for %d workers)",
                                 prefix, fMergersCount, activeWorkers);
                     else {
                        msg.Form("%s: No mergers will be used for %d workers",
                                 prefix, activeWorkers);
                        fMergersCount = -1;
                     }
                     if (gProofServ)
                        gProofServ->SendAsynMessage(msg);
                     else
                        Printf("%s",msg.Data());
                  } else if (fMergersByHost) {
                     // We force mergers at host level to minimize network traffic
                     if (activeWorkers > 1) {
                        fMergersCount = 0;
                        THashList hosts;
                        TIter nxwk(fSlaves);
                        TObject *wrk = 0;
                        while ((wrk = nxwk())) {
                           if (!hosts.FindObject(wrk->GetName())) {
                              hosts.Add(new TObjString(wrk->GetName()));
                              fMergersCount++;
                           }
                        }
                     }
                     if (fMergersCount > 1)
                        msg.Form("%s: Number of mergers set to %d (for %d workers), one for each slave host",
                                 prefix, fMergersCount, activeWorkers);
                     else {
                        msg.Form("%s: No mergers will be used for %d workers",
                                 prefix, activeWorkers);
                        fMergersCount = -1;
                     }
                     if (gProofServ)
                        gProofServ->SendAsynMessage(msg);
                     else
                        Printf("%s",msg.Data());
                  } else {
                     msg.Form("%s: Number of mergers set by user to %d (for %d workers)",
                              prefix, fMergersCount, activeWorkers);
                     if (gProofServ)
                        gProofServ->SendAsynMessage(msg);
                     else
                        Printf("%s",msg.Data());
                  }

                  // We started merging; we call it here because fMergersCount is still the original number
                  // and can be saved internally
                  fPlayer->SetMerging(kTRUE);

                  // Update merger counters (new workers are not yet active)
                  fMergePrg.SetNWrks(fMergersCount);

                  if (fMergersCount > 0) {

                     fMergers = new TList();
                     fLastAssignedMerger = 0;
                     // Total number of workers, which will not act as mergers ('pure workers')
                     fWorkersToMerge = (activeWorkers - fMergersCount);
                     // Establish the first merger
                     if (!CreateMerger(sl, merging_port)) {
                        // Cannot establish first merger
                        AskForOutput(sl);
                        fWorkersToMerge--;
                        fMergersCount--;
                     }
                     if (IsLite()) fMergePrg.SetNWrks(fMergersCount);
                  } else {
                     AskForOutput(sl);
                  }
                  fMergersSet = kTRUE;
               } else {
                  // Multiple pass
                  if (fMergersCount == -1) {
                     // No mergers. Workers send their outputs directly to master
                     AskForOutput(sl);
                  } else {
                     if ((fRedirectNext > 0 ) && (!fMergersByHost)) {
                        RedirectWorker(s, sl, output_size);
                        fRedirectNext--;
                     } else {
                        Bool_t newMerger = kTRUE;
                        if (fMergersByHost) {
                           TIter nxmg(fMergers);
                           TMergerInfo *mgi = 0;
                           while ((mgi = (TMergerInfo *) nxmg())) {
                              if (!strcmp(sl->GetName(), mgi->GetMerger()->GetName())) {
                                 newMerger = kFALSE;
                                 break;
                              }
                           }
                        }
                        if ((fMergersCount > fMergers->GetSize()) && newMerger) {
                           // Still not enough mergers established
                           if (!CreateMerger(sl, merging_port)) {
                              // Cannot establish a merger
                              AskForOutput(sl);
                              fWorkersToMerge--;
                              fMergersCount--;
                           }
                        } else
                           RedirectWorker(s, sl, output_size);
                     }
                  }
               }
            } else {
               Error("HandleSubMerger","kOutputSize received not on endmaster!");
            }
         }
      break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Redirect output of worker sl to some merger

void TProof::RedirectWorker(TSocket *s, TSlave * sl, Int_t output_size)
{
   Int_t merger_id = -1;

   if (fMergersByHost) {
      for (Int_t i = 0; i < fMergers->GetSize(); i++) {
         TMergerInfo *mgi = (TMergerInfo *)fMergers->At(i);
         if (!strcmp(sl->GetName(), mgi->GetMerger()->GetName())) {
            merger_id = i;
            break;
         }
      }
   } else {
      merger_id  = FindNextFreeMerger();
   }

   if (merger_id == -1) {
      // No free merger (probably it had crashed before)
      AskForOutput(sl);
   } else {
      TMessage sendoutput(kPROOF_SUBMERGER);
      sendoutput << Int_t(kSendOutput);
      PDB(kSubmerger, 2)
         Info("RedirectWorker", "redirecting worker %s to merger %d", sl->GetOrdinal(), merger_id);

       PDB(kSubmerger, 2) Info("RedirectWorker", "redirecting output to merger #%d",  merger_id);
       if (!fMergers || fMergers->GetSize() <= merger_id) {
          Error("RedirectWorker", "#%d not in list ", merger_id);
          return;
       }
       TMergerInfo * mi = (TMergerInfo *) fMergers->At(merger_id);

       TString hname = (IsLite()) ? "localhost" : mi->GetMerger()->GetName();
       sendoutput <<  merger_id;
       sendoutput << hname;
       sendoutput << mi->GetPort();
       s->Send(sendoutput);
       mi->AddMergedObjects(output_size);
       mi->AddWorker(sl);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return a merger, which is both active and still accepts some workers to be
/// assigned to it. It works on the 'round-robin' basis.

Int_t TProof::FindNextFreeMerger()
{
   while (fLastAssignedMerger < fMergers->GetSize() &&
         (!((TMergerInfo*)fMergers->At(fLastAssignedMerger))->IsActive() ||
           ((TMergerInfo*)fMergers->At(fLastAssignedMerger))->AreAllWorkersAssigned())) {
      fLastAssignedMerger++;
   }

   if (fLastAssignedMerger == fMergers->GetSize()) {
      fLastAssignedMerger = 0;
   } else {
      return fLastAssignedMerger++;
   }

   while (fLastAssignedMerger < fMergers->GetSize() &&
         (!((TMergerInfo*)fMergers->At(fLastAssignedMerger))->IsActive() ||
           ((TMergerInfo*)fMergers->At(fLastAssignedMerger))->AreAllWorkersAssigned())) {
      fLastAssignedMerger++;
   }

   if (fLastAssignedMerger == fMergers->GetSize()) {
      return -1;
   } else {
      return fLastAssignedMerger++;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Master asks for output from worker sl

void TProof::AskForOutput(TSlave *sl)
{
   TMessage sendoutput(kPROOF_SUBMERGER);
   sendoutput << Int_t(kSendOutput);

   PDB(kSubmerger, 2) Info("AskForOutput",
                           "worker %s was asked to send its output to master",
                            sl->GetOrdinal());

   sendoutput << -1;
   sendoutput << TString("master");
   sendoutput << -1;
   sl->GetSocket()->Send(sendoutput);
   if (IsLite()) fMergePrg.IncreaseNWrks();
}

////////////////////////////////////////////////////////////////////////////////
/// Final update of the progress dialog

void TProof::UpdateDialog()
{
   if (!fPlayer) return;

   // Handle abort ...
   if (fPlayer->GetExitStatus() == TVirtualProofPlayer::kAborted) {
      if (fSync)
         Info("UpdateDialog",
              "processing was aborted - %lld events processed",
              fPlayer->GetEventsProcessed());

      if (GetRemoteProtocol() > 11) {
         // New format
         Progress(-1, fPlayer->GetEventsProcessed(), -1, -1., -1., -1., -1.);
      } else {
         Progress(-1, fPlayer->GetEventsProcessed());
      }
      Emit("StopProcess(Bool_t)", kTRUE);
   }

   // Handle stop ...
   if (fPlayer->GetExitStatus() == TVirtualProofPlayer::kStopped) {
      if (fSync)
         Info("UpdateDialog",
              "processing was stopped - %lld events processed",
              fPlayer->GetEventsProcessed());

      if (GetRemoteProtocol() > 25) {
         // New format
         Progress(-1, fPlayer->GetEventsProcessed(), -1, -1., -1., -1., -1., -1, -1, -1.);
      } else if (GetRemoteProtocol() > 11) {
         Progress(-1, fPlayer->GetEventsProcessed(), -1, -1., -1., -1., -1.);
      } else {
         Progress(-1, fPlayer->GetEventsProcessed());
      }
      Emit("StopProcess(Bool_t)", kFALSE);
   }

   // Final update of the dialog box
   if (GetRemoteProtocol() > 25) {
      // New format
      EmitVA("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)",
              10, (Long64_t)(-1), (Long64_t)(-1), (Long64_t)(-1),(Float_t)(-1.),(Float_t)(-1.),
                  (Float_t)(-1.),(Float_t)(-1.),(Int_t)(-1),(Int_t)(-1),(Float_t)(-1.));
   } else if (GetRemoteProtocol() > 11) {
      // New format
      EmitVA("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)",
               7, (Long64_t)(-1), (Long64_t)(-1), (Long64_t)(-1),
                  (Float_t)(-1.),(Float_t)(-1.),(Float_t)(-1.),(Float_t)(-1.));
   } else {
      EmitVA("Progress(Long64_t,Long64_t)", 2, (Long64_t)(-1), (Long64_t)(-1));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Activate the a-sync input handler.

void TProof::ActivateAsyncInput()
{
   TIter next(fSlaves);
   TSlave *sl;

   while ((sl = (TSlave*) next()))
      if (sl->GetInputHandler())
         sl->GetInputHandler()->Add();
}

////////////////////////////////////////////////////////////////////////////////
/// De-activate a-sync input handler.

void TProof::DeActivateAsyncInput()
{
   TIter next(fSlaves);
   TSlave *sl;

   while ((sl = (TSlave*) next()))
      if (sl->GetInputHandler())
         sl->GetInputHandler()->Remove();
}

////////////////////////////////////////////////////////////////////////////////
/// Get the active mergers count

Int_t TProof::GetActiveMergersCount()
{
   if (!fMergers) return 0;

   Int_t active_mergers = 0;

   TIter mergers(fMergers);
   TMergerInfo *mi = 0;
   while ((mi = (TMergerInfo *)mergers())) {
      if (mi->IsActive()) active_mergers++;
   }

   return active_mergers;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new merger

Bool_t TProof::CreateMerger(TSlave *sl, Int_t port)
{
   PDB(kSubmerger, 2)
      Info("CreateMerger", "worker %s will be merger ", sl->GetOrdinal());

   PDB(kSubmerger, 2) Info("CreateMerger","Begin");

   if (port <= 0) {
      PDB(kSubmerger,2)
         Info("CreateMerger", "cannot create merger on port %d - exit", port);
      return kFALSE;
   }

   Int_t workers = -1;
   if (!fMergersByHost) {
      Int_t mergersToCreate = fMergersCount - fMergers->GetSize();
      // Number of pure workers, which are not simply divisible by mergers
      Int_t rest = fWorkersToMerge % mergersToCreate;
      // We add one more worker for each of the first 'rest' mergers being established
      if (rest > 0 && fMergers->GetSize() < rest) {
         rest = 1;
      } else {
         rest = 0;
      }
      workers = (fWorkersToMerge / mergersToCreate) + rest;
   } else {
      Int_t workersOnHost = 0;
      for (Int_t i = 0; i < fActiveSlaves->GetSize(); i++) {
         if(!strcmp(sl->GetName(), fActiveSlaves->At(i)->GetName())) workersOnHost++;
      }
      workers = workersOnHost - 1;
   }

   TString msg;
   msg.Form("worker %s on host %s will be merger for %d additional workers", sl->GetOrdinal(), sl->GetName(), workers);

   if (gProofServ) {
      gProofServ->SendAsynMessage(msg);
   } else {
      Printf("%s",msg.Data());
   }
   TMergerInfo * merger = new TMergerInfo(sl, port, workers);

   TMessage bemerger(kPROOF_SUBMERGER);
   bemerger << Int_t(kBeMerger);
   bemerger <<  fMergers->GetSize();
   bemerger <<  workers;
   sl->GetSocket()->Send(bemerger);

   PDB(kSubmerger,2) Info("CreateMerger",
                          "merger #%d (port: %d) for %d workers started",
                          fMergers->GetSize(), port, workers);

   fMergers->Add(merger);
   fWorkersToMerge = fWorkersToMerge - workers;

   fRedirectNext = workers / 2;

   PDB(kSubmerger, 2) Info("CreateMerger", "exit");
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a bad slave server to the bad slave list and remove it from
/// the active list and from the two monitor objects. Assume that the work
/// done by this worker was lost and ask packerizer to reassign it.

void TProof::MarkBad(TSlave *wrk, const char *reason)
{
   std::lock_guard<std::recursive_mutex> lock(fCloseMutex);

   // We may have been invalidated in the meanwhile: nothing to do in such a case
   if (!IsValid()) return;

   if (!wrk) {
      Error("MarkBad", "worker instance undefined: protocol error? ");
      return;
   }

   // Local URL
   static TString thisurl;
   if (thisurl.IsNull()) {
      if (IsMaster()) {
         Int_t port = gEnv->GetValue("ProofServ.XpdPort",-1);
         thisurl = TUrl(gSystem->HostName()).GetHostFQDN();
         if (port > 0) thisurl += TString::Format(":%d", port);
      } else {
         thisurl.Form("%s@%s:%d", fUrl.GetUser(), fUrl.GetHost(), fUrl.GetPort());
      }
   }

   if (!reason || (strcmp(reason, kPROOF_TerminateWorker) && strcmp(reason, kPROOF_WorkerIdleTO))) {
      // Message for notification
      const char *mastertype = (gProofServ && gProofServ->IsTopMaster()) ? "top master" : "master";
      TString src = IsMaster() ? Form("%s at %s", mastertype, thisurl.Data()) : "local session";
      TString msg;
      msg.Form("\n +++ Message from %s : marking %s:%d (%s) as bad\n +++ Reason: %s",
               src.Data(), wrk->GetName(), wrk->GetPort(), wrk->GetOrdinal(),
               (reason && strlen(reason)) ? reason : "unknown");
      Info("MarkBad", "%s", msg.Data());
      // Notify one level up, if the case
      // Add some hint for diagnostics
      if (gProofServ) {
         msg += TString::Format("\n\n +++ Most likely your code crashed on worker %s at %s:%d.\n",
                     wrk->GetOrdinal(), wrk->GetName(), wrk->GetPort());
      } else {
         msg += TString::Format("\n\n +++ Most likely your code crashed\n");
      }
      msg += TString::Format(" +++ Please check the session logs for error messages either using\n");
      msg += TString::Format(" +++ the 'Show logs' button or executing\n");
      msg += TString::Format(" +++\n");
      if (gProofServ) {
         msg += TString::Format(" +++ root [] TProof::Mgr(\"%s\")->GetSessionLogs()->"
                                "Display(\"%s\",0)\n\n", thisurl.Data(), wrk->GetOrdinal());
         gProofServ->SendAsynMessage(msg, kTRUE);
      } else {
         msg += TString::Format(" +++ root [] TProof::Mgr(\"%s\")->GetSessionLogs()->"
                                "Display(\"*\")\n\n", thisurl.Data());
         Printf("%s", msg.Data());
      }
   } else if (reason) {
      if (gDebug > 0 && strcmp(reason, kPROOF_WorkerIdleTO)) {
         Info("MarkBad", "worker %s at %s:%d asked to terminate",
                         wrk->GetOrdinal(), wrk->GetName(), wrk->GetPort());
      }
   }

   if (IsMaster() && reason) {
      if (strcmp(reason, kPROOF_TerminateWorker)) {
         // if the reason was not a planned termination
         TList *listOfMissingFiles = 0;
         if (!(listOfMissingFiles = (TList *)GetOutput("MissingFiles"))) {
            listOfMissingFiles = new TList();
            listOfMissingFiles->SetName("MissingFiles");
            if (fPlayer)
               fPlayer->AddOutputObject(listOfMissingFiles);
         }
         // If a query is being processed, assume that the work done by
         // the worker was lost and needs to be reassigned.
         TVirtualPacketizer *packetizer = fPlayer ? fPlayer->GetPacketizer() : 0;
         if (packetizer) {
            // the worker was lost so do resubmit the packets
            packetizer->MarkBad(wrk, 0, &listOfMissingFiles);
         }
      } else {
         // Tell the coordinator that we are gone
         if (gProofServ) {
            TString ord(wrk->GetOrdinal());
            Int_t id = ord.Last('.');
            if (id != kNPOS) ord.Remove(0, id+1);
            gProofServ->ReleaseWorker(ord.Data());
         }
      }
   } else if (TestBit(TProof::kIsClient) && reason && !strcmp(reason, kPROOF_WorkerIdleTO)) {
      // We are invalid after this
      fValid = kFALSE;
   }

   fActiveSlaves->Remove(wrk);
   FindUniqueSlaves();

   fAllMonitor->Remove(wrk->GetSocket());
   fActiveMonitor->Remove(wrk->GetSocket());

   fSendGroupView = kTRUE;

   if (IsMaster()) {
      if (reason && !strcmp(reason, kPROOF_TerminateWorker)) {
         // if the reason was a planned termination then delete the worker and
         // remove it from all the lists
         fSlaves->Remove(wrk);
         fBadSlaves->Remove(wrk);
         fActiveSlaves->Remove(wrk);
         fInactiveSlaves->Remove(wrk);
         fUniqueSlaves->Remove(wrk);
         fAllUniqueSlaves->Remove(wrk);
         fNonUniqueMasters->Remove(wrk);

         // we add it to the list of terminated slave infos instead, so that it
         // stays available in the .workers persistent file
         TSlaveInfo *si = new TSlaveInfo(
            wrk->GetOrdinal(),
            Form("%s@%s:%d", wrk->GetUser(), wrk->GetName(), wrk->GetPort()),
            0, "", wrk->GetWorkDir());
         if (!fTerminatedSlaveInfos->Contains(si)) fTerminatedSlaveInfos->Add(si);
         else delete si;

         delete wrk;
      } else {
         fBadSlaves->Add(wrk);
         fActiveSlaves->Remove(wrk);
         fUniqueSlaves->Remove(wrk);
         fAllUniqueSlaves->Remove(wrk);
         fNonUniqueMasters->Remove(wrk);
         if (fCurrentMonitor) fCurrentMonitor->DeActivate(wrk->GetSocket());
         wrk->Close();
         // Update the mergers count, if needed
         if (fMergersSet) {
            Int_t mergersCount = -1;
            TParameter<Int_t> *mc = dynamic_cast<TParameter<Int_t> *>(GetParameter("PROOF_UseMergers"));
            if (mc) mergersCount = mc->GetVal(); // Value set by user
            // Mergers count is set dynamically: recalculate it
            if (mergersCount == 0) {
               Int_t activeWorkers = fCurrentMonitor ? fCurrentMonitor->GetActive() : GetNumberOfActiveSlaves();
               if (activeWorkers > 1) {
                  fMergersCount = TMath::Nint(TMath::Sqrt(activeWorkers));
                  if (activeWorkers / fMergersCount < 2)
                     fMergersCount = (Int_t) TMath::Sqrt(activeWorkers);
               }
            }
         }
      }

      // Update session workers files
      SaveWorkerInfo();
   } else {
      // On clients the proof session should be removed from the lists
      // and deleted, since it is not valid anymore
      fSlaves->Remove(wrk);
      if (fManager)
         fManager->DiscardSession(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add slave with socket s to the bad slave list and remove if from
/// the active list and from the two monitor objects.

void TProof::MarkBad(TSocket *s, const char *reason)
{
   std::lock_guard<std::recursive_mutex> lock(fCloseMutex);

   // We may have been invalidated in the meanwhile: nothing to do in such a case
   if (!IsValid()) return;

   TSlave *wrk = FindSlave(s);
   MarkBad(wrk, reason);
}

////////////////////////////////////////////////////////////////////////////////
/// Ask an active worker 'wrk' to terminate, i.e. to shutdown

void TProof::TerminateWorker(TSlave *wrk)
{
   if (!wrk) {
      Warning("TerminateWorker", "worker instance undefined: protocol error? ");
      return;
   }

   // Send stop message
   if (wrk->GetSocket() && wrk->GetSocket()->IsValid()) {
      TMessage mess(kPROOF_STOP);
      wrk->GetSocket()->Send(mess);
   } else {
      if (gDebug > 0)
         Info("TerminateWorker", "connection to worker is already down: cannot"
                                 " send termination message");
   }

   // This is a bad worker from now on
   MarkBad(wrk, kPROOF_TerminateWorker);
}

////////////////////////////////////////////////////////////////////////////////
/// Ask an active worker 'ord' to terminate, i.e. to shutdown

void TProof::TerminateWorker(const char *ord)
{
   if (ord && strlen(ord) > 0) {
      Bool_t all = (ord[0] == '*') ? kTRUE : kFALSE;
      if (IsMaster()) {
         TIter nxw(fSlaves);
         TSlave *wrk = 0;
         while ((wrk = (TSlave *)nxw())) {
            if (all || !strcmp(wrk->GetOrdinal(), ord)) {
               TerminateWorker(wrk);
               if (!all) break;
            }
         }
      } else {
         TMessage mess(kPROOF_STOP);
         mess << TString(ord);
         Broadcast(mess);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Ping PROOF. Returns 1 if master server responded.

Int_t TProof::Ping()
{
   return Ping(kActive);
}

////////////////////////////////////////////////////////////////////////////////
/// Ping PROOF slaves. Returns the number of slaves that responded.

Int_t TProof::Ping(ESlaves list)
{
   TList *slaves = 0;
   if (list == kAll)       slaves = fSlaves;
   if (list == kActive)    slaves = fActiveSlaves;
   if (list == kUnique)    slaves = fUniqueSlaves;
   if (list == kAllUnique) slaves = fAllUniqueSlaves;

   if (slaves->GetSize() == 0) return 0;

   int   nsent = 0;
   TIter next(slaves);

   TSlave *sl;
   while ((sl = (TSlave *)next())) {
      if (sl->IsValid()) {
         if (sl->Ping() == -1) {
            MarkBad(sl, "ping unsuccessful");
         } else {
            nsent++;
         }
      }
   }

   return nsent;
}

////////////////////////////////////////////////////////////////////////////////
/// Ping PROOF slaves. Returns the number of slaves that responded.

void TProof::Touch()
{
   TList *slaves = fSlaves;

   if (slaves->GetSize() == 0) return;

   TIter next(slaves);

   TSlave *sl;
   while ((sl = (TSlave *)next())) {
      if (sl->IsValid()) {
         sl->Touch();
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Print status of PROOF cluster.

void TProof::Print(Option_t *option) const
{
   TString secCont;

   if (TestBit(TProof::kIsClient)) {
      Printf("Connected to:             %s (%s)", GetMaster(),
                                             IsValid() ? "valid" : "invalid");
      Printf("Port number:              %d", GetPort());
      Printf("User:                     %s", GetUser());
      Printf("ROOT version|rev:         %s|%s", gROOT->GetVersion(), gROOT->GetGitCommit());
      Printf("Architecture-Compiler:    %s-%s", gSystem->GetBuildArch(),
                                                gSystem->GetBuildCompilerVersion());
      TSlave *sl = (TSlave *)fActiveSlaves->First();
      if (sl) {
         TString sc;
         if (sl->GetSocket()->GetSecContext())
            Printf("Security context:         %s",
                                      sl->GetSocket()->GetSecContext()->AsString(sc));
         Printf("Proofd protocol version:  %d", sl->GetSocket()->GetRemoteProtocol());
      } else {
         Printf("Security context:         Error - No connection");
         Printf("Proofd protocol version:  Error - No connection");
      }
      Printf("Client protocol version:  %d", GetClientProtocol());
      Printf("Remote protocol version:  %d", GetRemoteProtocol());
      Printf("Log level:                %d", GetLogLevel());
      Printf("Session unique tag:       %s", IsValid() ? GetSessionTag() : "");
      Printf("Default data pool:        %s", IsValid() ? GetDataPoolUrl() : "");
      if (IsValid())
         const_cast<TProof*>(this)->SendPrint(option);
   } else {
      const_cast<TProof*>(this)->AskStatistics();
      if (IsParallel())
         Printf("*** Master server %s (parallel mode, %d workers):",
                gProofServ->GetOrdinal(), GetParallel());
      else
         Printf("*** Master server %s (sequential mode):",
                gProofServ->GetOrdinal());

      Printf("Master host name:           %s", gSystem->HostName());
      Printf("Port number:                %d", GetPort());
      if (strlen(gProofServ->GetGroup()) > 0) {
         Printf("User/Group:                 %s/%s", GetUser(), gProofServ->GetGroup());
      } else {
         Printf("User:                       %s", GetUser());
      }
      TString ver;
      ver.Form("%s|%s", gROOT->GetVersion(), gROOT->GetGitCommit());
      if (gSystem->Getenv("ROOTVERSIONTAG"))
         ver.Form("%s|%s", gROOT->GetVersion(), gSystem->Getenv("ROOTVERSIONTAG"));
      Printf("ROOT version|rev|tag:       %s", ver.Data());
      Printf("Architecture-Compiler:      %s-%s", gSystem->GetBuildArch(),
                                                  gSystem->GetBuildCompilerVersion());
      Printf("Protocol version:           %d", GetClientProtocol());
      Printf("Image name:                 %s", GetImage());
      Printf("Working directory:          %s", gSystem->WorkingDirectory());
      Printf("Config directory:           %s", GetConfDir());
      Printf("Config file:                %s", GetConfFile());
      Printf("Log level:                  %d", GetLogLevel());
      Printf("Number of workers:          %d", GetNumberOfSlaves());
      Printf("Number of active workers:   %d", GetNumberOfActiveSlaves());
      Printf("Number of unique workers:   %d", GetNumberOfUniqueSlaves());
      Printf("Number of inactive workers: %d", GetNumberOfInactiveSlaves());
      Printf("Number of bad workers:      %d", GetNumberOfBadSlaves());
      Printf("Total MB's processed:       %.2f", float(GetBytesRead())/(1024*1024));
      Printf("Total real time used (s):   %.3f", GetRealTime());
      Printf("Total CPU time used (s):    %.3f", GetCpuTime());
      if (TString(option).Contains("a", TString::kIgnoreCase) && GetNumberOfSlaves()) {
         Printf("List of workers:");
         TList masters;
         TIter nextslave(fSlaves);
         while (TSlave* sl = dynamic_cast<TSlave*>(nextslave())) {
            if (!sl->IsValid()) continue;

            if (sl->GetSlaveType() == TSlave::kSlave) {
               sl->Print(option);
            } else if (sl->GetSlaveType() == TSlave::kMaster) {
               TMessage mess(kPROOF_PRINT);
               mess.WriteString(option);
               if (sl->GetSocket()->Send(mess) == -1)
                  const_cast<TProof*>(this)->MarkBad(sl, "could not send kPROOF_PRINT request");
               else
                  masters.Add(sl);
            } else {
               Error("Print", "TSlave is neither Master nor Worker");
               R__ASSERT(0);
            }
         }
         const_cast<TProof*>(this)->Collect(&masters, fCollectTimeout);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Extract from opt information about output handling settings.
/// The understood keywords are:
///     of=<file>, outfile=<file>         output file location
///     ds=<dsname>, dataset=<dsname>     dataset name ('of' and 'ds' are
///                                       mutually exclusive,execution stops
///                                       if both are found)
///     sft[=<opt>], savetofile[=<opt>]   control saving to file
///
/// For 'mvf', the <opt> integer has the following meaning:
///     <opt> = <how>*10 + <force>
///             <force> = 0      save to file if memory threshold is reached
///                              (the memory threshold is set by the cluster
///                              admin); in case an output file is defined, the
///                              files are merged at the end;
///                       1      save results to file.
///             <how> =   0      save at the end of the query
///                       1      save results after each packet (to reduce the
///                              loss in case of crash).
///
/// Setting 'ds' automatically sets 'mvf=1'; it is still possible to set 'mvf=11'
/// to save results after each packet.
///
/// The separator from the next option is either a ' ' or a ';'
///
/// All recognized settings are removed from the input string opt.
/// If action == 0, set up the output file accordingly, if action == 1 clean related
/// output file settings.
/// If the final target file is local then 'target' is set to the final local path
/// when action == 0 and used to retrieve the file with TFile::Cp when action == 1.
///
/// Output file settings are in the form
///
///       <previous_option>of=name <next_option>
///       <previous_option>outfile=name,...;<next_option>
///
/// The separator from the next option is either a ' ' or a ';'
/// Called interanally by TProof::Process.
///
/// Returns 0 on success, -1 on error.

Int_t TProof::HandleOutputOptions(TString &opt, TString &target, Int_t action)
{
   TString outfile, dsname, stfopt;
   if (action == 0) {
      TString tagf, tagd, tags, oo;
      Ssiz_t from = 0, iof = kNPOS, iod = kNPOS, ios = kNPOS;
      while (opt.Tokenize(oo, from, "[; ]")) {
         if (oo.BeginsWith("of=")) {
            tagf = "of=";
            iof = opt.Index(tagf);
         } else if (oo.BeginsWith("outfile=")) {
            tagf = "outfile=";
            iof = opt.Index(tagf);
         } else if (oo.BeginsWith("ds")) {
            tagd = "ds";
            iod = opt.Index(tagd);
         } else if (oo.BeginsWith("dataset")) {
            tagd = "dataset";
            iod = opt.Index(tagd);
         } else if (oo.BeginsWith("stf")) {
            tags = "stf";
            ios = opt.Index(tags);
         } else if (oo.BeginsWith("savetofile")) {
            tags = "savetofile";
            ios = opt.Index(tags);
         }
      }
      // Check consistency
      if (iof != kNPOS && iod != kNPOS) {
         Error("HandleOutputOptions", "options 'of'/'outfile' and 'ds'/'dataset' are incompatible!");
         return -1;
      }

      // Check output file first
      if (iof != kNPOS) {
         from = iof + tagf.Length();
         if (!opt.Tokenize(outfile, from, "[; ]") || outfile.IsNull()) {
            Error("HandleOutputOptions", "could not extract output file settings string! (%s)", opt.Data());
            return -1;
         }
         // For removal from original options string
         tagf += outfile;
      }
      // Check dataset
      if (iod != kNPOS) {
         from = iod + tagd.Length();
         if (!opt.Tokenize(dsname, from, "[; ]"))
            if (gDebug > 0) Info("HandleOutputOptions", "no dataset name found: use default");
         // For removal from original options string
         tagd += dsname;
         // The name may be empty or beginning with a '='
         if (dsname.BeginsWith("=")) dsname.Replace(0, 1, "");
         if (dsname.Contains("|V")) {
            target = "ds|V";
            dsname.ReplaceAll("|V", "");
         }
         if (dsname.IsNull()) dsname = "dataset_<qtag>";
      }
      // Check stf
      if (ios != kNPOS) {
         from = ios + tags.Length();
         if (!opt.Tokenize(stfopt, from, "[; ]"))
            if (gDebug > 0) Info("HandleOutputOptions", "save-to-file not found: use default");
         // For removal from original options string
         tags += stfopt;
         // It must be digit
         if (!stfopt.IsNull()) {
            if (stfopt.BeginsWith("=")) stfopt.Replace(0,1,"");
            if (!stfopt.IsNull()) {
               if (!stfopt.IsDigit()) {
                  Error("HandleOutputOptions", "save-to-file option must be a digit! (%s)", stfopt.Data());
                  return -1;
               }
            } else {
               // Default
               stfopt = "1";
            }
         } else {
            // Default
            stfopt = "1";
         }
      }
      // Remove from original options string
      opt.ReplaceAll(tagf, "");
      opt.ReplaceAll(tagd, "");
      opt.ReplaceAll(tags, "");
   }

   // Parse now
   if (action == 0) {
      // Output file
      if (!outfile.IsNull()) {
         if (!outfile.BeginsWith("master:")) {
            if (gSystem->AccessPathName(gSystem->GetDirName(outfile.Data()), kWritePermission)) {
               Warning("HandleOutputOptions",
                     "directory '%s' for the output file does not exists or is not writable:"
                     " saving to master", gSystem->GetDirName(outfile.Data()).Data());
               outfile.Form("master:%s", gSystem->BaseName(outfile.Data()));
            } else {
               if (!IsLite()) {
                  // The target file is local, so we need to retrieve it
                  target = outfile;
                  if (!stfopt.IsNull()) {
                     outfile.Form("master:%s", gSystem->BaseName(target.Data()));
                  } else {
                     outfile = "";
                  }
               }
            }
         }
         if (outfile.BeginsWith("master:")) {
            outfile.ReplaceAll("master:", "");
            if (outfile.IsNull() || !gSystem->IsAbsoluteFileName(outfile)) {
               // Get the master data dir
               TString ddir, emsg;
               if (!IsLite()) {
                  if (Exec("gProofServ->GetDataDir()", "0", kTRUE) == 0) {
                     TObjString *os = fMacroLog.GetLineWith("const char");
                     if (os) {
                        Ssiz_t fst =  os->GetString().First('\"');
                        Ssiz_t lst =  os->GetString().Last('\"');
                        ddir = os->GetString()(fst+1, lst-fst-1);
                     } else {
                        emsg = "could not find 'const char *' string in macro log! cannot continue";
                     }
                  } else {
                     emsg = "could not retrieve master data directory info! cannot continue";
                  }
                  if (!emsg.IsNull()) {
                     Error("HandleOutputOptions", "%s", emsg.Data());
                     return -1;
                  }
               }
               if (!ddir.IsNull()) ddir += "/";
               if (outfile.IsNull()) {
                  outfile.Form("%s<file>", ddir.Data());
               } else {
                  outfile.Insert(0, TString::Format("%s", ddir.Data()));
               }
            }
         }
         // Set the parameter
         if (!outfile.IsNull()) {
            if (!outfile.BeginsWith("of:")) outfile.Insert(0, "of:");
            SetParameter("PROOF_DefaultOutputOption", outfile.Data());
         }
      }
      // Dataset creation
      if (!dsname.IsNull()) {
         dsname.Insert(0, "ds:");
         // Set the parameter
         SetParameter("PROOF_DefaultOutputOption", dsname.Data());
         // Check the Save-To-File option
         if (!stfopt.IsNull()) {
            Int_t ostf = (Int_t) stfopt.Atoi();
            if (ostf%10 <= 0) {
               Warning("HandleOutputOptions", "Dataset required bu Save-To-File disabled: enabling!");
               stfopt.Form("%d", ostf+1);
            }
         } else {
            // Minimal setting
            stfopt = "1";
         }
      }
      // Save-To-File options
      if (!stfopt.IsNull()) {
         // Set the parameter
         SetParameter("PROOF_SavePartialResults", (Int_t) stfopt.Atoi());
      }
   } else {
      // Retrieve the file, if required
      if (GetOutputList()) {
         if (target == "ds|V") {
            // Find the dataset
            dsname = "";
            TIter nxo(GetOutputList());
            TObject *o = 0;
            while ((o = nxo())) {
               if (o->InheritsFrom(TFileCollection::Class())) {
                  VerifyDataSet(o->GetName());
                  dsname = o->GetName();
                  break;
               }
            }
            if (!dsname.IsNull()) {
               TFileCollection *fc = GetDataSet(dsname);
               if (fc) {
                  fc->Print();
               } else {
                  Warning("HandleOutputOptions", "could not retrieve TFileCollection for dataset '%s'", dsname.Data());
               }
            } else {
               Warning("HandleOutputOptions", "dataset not found!");
            }
         } else {
            Bool_t targetcopied = kFALSE;
            TProofOutputFile *pf = 0;
            if (!target.IsNull())
               pf = (TProofOutputFile *) GetOutputList()->FindObject(gSystem->BaseName(target.Data()));
            if (pf) {
               // Copy the file
               if (strcmp(TUrl(pf->GetOutputFileName(), kTRUE).GetUrl(),
                          TUrl(target, kTRUE).GetUrl())) {
                  if (TFile::Cp(pf->GetOutputFileName(), target)) {
                     Printf(" Output successfully copied to %s", target.Data());
                     targetcopied = kTRUE;
                  } else {
                     Warning("HandleOutputOptions", "problems copying output to %s", target.Data());
                  }
               }
            }
            TFile *fout = 0;
            TObject *o = 0;
            TIter nxo(GetOutputList());
            Bool_t swapcopied = kFALSE;
            while ((o = nxo())) {
               TProofOutputFile *pof = dynamic_cast<TProofOutputFile *>(o);
               if (pof) {
                  if (pof->TestBit(TProofOutputFile::kSwapFile) && !target.IsNull()) {
                     if (pof == pf && targetcopied) continue;
                     // Copy the file
                     if (strcmp(TUrl(pf->GetOutputFileName(), kTRUE).GetUrl(),
                                TUrl(target, kTRUE).GetUrl())) {
                        if (TFile::Cp(pof->GetOutputFileName(), target)) {
                           Printf(" Output successfully copied to %s", target.Data());
                           swapcopied = kTRUE;
                        } else {
                           Warning("HandleOutputOptions", "problems copying output to %s", target.Data());
                        }
                     }
                  } else if (pof->IsRetrieve()) {
                     // Retrieve this file to the local path indicated in the title
                     if (strcmp(TUrl(pf->GetOutputFileName(), kTRUE).GetUrl(),
                                TUrl(pof->GetTitle(), kTRUE).GetUrl())) {
                        if (TFile::Cp(pof->GetOutputFileName(), pof->GetTitle())) {
                           Printf(" Output successfully copied to %s", pof->GetTitle());
                        } else {
                           Warning("HandleOutputOptions",
                                 "problems copying %s to %s", pof->GetOutputFileName(), pof->GetTitle());
                        }
                     }
                  }
               }
            }
            if (!target.IsNull() && !swapcopied) {
               if (!fout && !pf) {
                  fout = TFile::Open(target, "RECREATE");
                  if (!fout || (fout && fout->IsZombie())) {
                     SafeDelete(fout);
                     Warning("HandleOutputOptions", "problems opening output file %s", target.Data());
                  }
               }
               if (fout) {
                  nxo.Reset();
                  while ((o = nxo())) {
                     TProofOutputFile *pof = dynamic_cast<TProofOutputFile *>(o);
                     if (!pof) {
                        // Write the object to the open output file
                        o->Write();
                     }
                  }
               }
            }
            // Clean-up
            if (fout) {
               fout->Close();
               SafeDelete(fout);
               Printf(" Output saved to %s", target.Data());
            }
         }
      }
      // Remove the parameter
      DeleteParameters("PROOF_DefaultOutputOption");
      // Remove the parameter
      DeleteParameters("PROOF_SavePartialResults");
   }
   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Extract from opt in optfb information about wanted feedback settings.
/// Feedback are removed from the input string opt.
/// If action == 0, set up feedback accordingly, if action == 1 clean related
/// feedback settings (using info in optfb, if available, or reparsing opt).
///
/// Feedback requirements are in the form
///
///       <previous_option>fb=name1,name2,name3,... <next_option>
///       <previous_option>feedback=name1,name2,name3,...;<next_option>
///
/// The special name 'stats' triggers feedback about events and packets.
/// The separator from the next option is either a ' ' or a ';'.
/// Called interanally by TProof::Process.

void TProof::SetFeedback(TString &opt, TString &optfb, Int_t action)
{
   Ssiz_t from = 0;
   if (action == 0 || (action == 1 && optfb.IsNull())) {
      TString tag("fb=");
      Ssiz_t ifb = opt.Index(tag);
      if (ifb == kNPOS) {
         tag = "feedback=";
         ifb = opt.Index(tag);
      }
      if (ifb == kNPOS) return;
      from = ifb + tag.Length();

      if (!opt.Tokenize(optfb, from, "[; ]") || optfb.IsNull()) {
         Warning("SetFeedback", "could not extract feedback string! Ignoring ...");
         return;
      }
      // Remove from original options string
      tag += optfb;
      opt.ReplaceAll(tag, "");
   }

   // Parse now
   TString nm, startdraw, stopdraw;
   from = 0;
   while (optfb.Tokenize(nm, from, ",")) {
      // Special name first
      if (nm == "stats") {
         if (action == 0) {
            startdraw.Form("gDirectory->Add(new TStatsFeedback((TProof *)%p))", this);
            gROOT->ProcessLine(startdraw.Data());
            SetParameter("PROOF_StatsHist", "");
            AddFeedback("PROOF_EventsHist");
            AddFeedback("PROOF_PacketsHist");
            AddFeedback("PROOF_ProcPcktHist");
         } else {
            stopdraw.Form("TObject *o = gDirectory->FindObject(\"%s\"); "
                          " if (o && strcmp(o->ClassName(), \"TStatsFeedback\")) "
                          " { gDirectory->Remove(o); delete o; }", GetSessionTag());
            gROOT->ProcessLine(stopdraw.Data());
            DeleteParameters("PROOF_StatsHist");
            RemoveFeedback("PROOF_EventsHist");
            RemoveFeedback("PROOF_PacketsHist");
            RemoveFeedback("PROOF_ProcPcktHist");
         }
      } else {
         if (action == 0) {
            // Enable or
            AddFeedback(nm);
            startdraw.Form("gDirectory->Add(new TDrawFeedback((TProof *)%p))", this);
            gROOT->ProcessLine(startdraw.Data());
         } else {
            // ... or disable
            RemoveFeedback(nm);
            stopdraw.Form("TObject *o = gDirectory->FindObject(\"%s\"); "
                          " if (o && strcmp(o->ClassName(), \"TDrawFeedback\")) "
                          " { gDirectory->Remove(o); delete o; }", GetSessionTag());
            gROOT->ProcessLine(stopdraw.Data());
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Process a data set (TDSet) using the specified selector (.C) file or
/// Tselector object
/// Entry- or event-lists should be set in the data set object using
/// TDSet::SetEntryList.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.

Long64_t TProof::Process(TDSet *dset, const char *selector, Option_t *option,
                         Long64_t nentries, Long64_t first)
{
   if (!IsValid() || !fPlayer) return -1;

   // Set PROOF to running state
   SetRunStatus(TProof::kRunning);

   TString opt(option), optfb, outfile;
   // Enable feedback, if required
   if (opt.Contains("fb=") || opt.Contains("feedback=")) SetFeedback(opt, optfb, 0);
   // Define output file, either from 'opt' or the default one
   if (HandleOutputOptions(opt, outfile, 0) != 0) return -1;

   // Resolve query mode
   fSync = (GetQueryMode(opt) == kSync);

   if (fSync && (!IsIdle() || IsWaiting())) {
      // Already queued or processing queries: switch to asynchronous mode
      Info("Process", "session is in waiting or processing status: switch to asynchronous mode");
      fSync = kFALSE;
      opt.ReplaceAll("SYNC","");
      opt += "ASYN";
   }

   // Cleanup old temporary datasets
   if ((IsIdle() && !IsWaiting()) && fRunningDSets && fRunningDSets->GetSize() > 0) {
      fRunningDSets->SetOwner(kTRUE);
      fRunningDSets->Delete();
   }

   // deactivate the default application interrupt handler
   // ctrl-c's will be forwarded to PROOF to stop the processing
   TSignalHandler *sh = 0;
   if (fSync) {
      if (gApplication)
         sh = gSystem->RemoveSignalHandler(gApplication->GetSignalHandler());
   }

   // Make sure we get a fresh result
   fOutputList.Clear();

   // Make sure that the workers ready list is empty
   if (fWrksOutputReady) {
      fWrksOutputReady->SetOwner(kFALSE);
      fWrksOutputReady->Clear();
   }

   // Make sure the selector path is in the macro path
   TProof::AssertMacroPath(selector);

   // Reset time measurements
   fQuerySTW.Reset();

   Long64_t rv = -1;
   if (selector && strlen(selector)) {
      rv = fPlayer->Process(dset, selector, opt.Data(), nentries, first);
   } else if (fSelector) {
      rv = fPlayer->Process(dset, fSelector, opt.Data(), nentries, first);
   } else {
      Error("Process", "neither a selecrot file nor a selector object have"
                       " been specified: cannot process!");
   }

   // This is the end of merging
   fQuerySTW.Stop();
   Float_t rt = fQuerySTW.RealTime();
   // Update the query content
   TQueryResult *qr = GetQueryResult();
   if (qr) {
      qr->SetTermTime(rt);
      qr->SetPrepTime(fPrepTime);
   }

   // Disable feedback, if required
   if (!optfb.IsNull()) SetFeedback(opt, optfb, 1);
   // Finalise output file settings (opt is ignored in here)
   if (HandleOutputOptions(opt, outfile, 1) != 0) return -1;

   // Retrieve status from the output list
   if (rv >= 0) {
      TParameter<Long64_t> *sst =
        (TParameter<Long64_t> *) fOutputList.FindObject("PROOF_SelectorStatus");
      if (sst) rv = sst->GetVal();
   }

   if (fSync) {
      // reactivate the default application interrupt handler
      if (sh)
         gSystem->AddSignalHandler(sh);
      // Save the performance info, if required
      if (!fPerfTree.IsNull()) {
         if (SavePerfTree() != 0) Error("Process", "saving performance info ...");
         // Must be re-enabled each time
         SetPerfTree(0);
      }
   }

   return rv;
}

////////////////////////////////////////////////////////////////////////////////
/// Process a data set (TFileCollection) using the specified selector (.C) file
/// or TSelector object.
/// The default tree is analyzed (i.e. the first one found). To specify another
/// tree, the default tree can be changed using TFileCollection::SetDefaultMetaData .
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.

Long64_t TProof::Process(TFileCollection *fc, const char *selector,
                         Option_t *option, Long64_t nentries, Long64_t first)
{
   if (!IsValid() || !fPlayer) return -1;

   if (fProtocol < 17) {
      Info("Process", "server version < 5.18/00:"
                      " processing of TFileCollection not supported");
      return -1;
   }

   // We include the TFileCollection to the input list and we create a
   // fake TDSet with infor about it
   TDSet *dset = new TDSet(TString::Format("TFileCollection:%s", fc->GetName()), 0, 0, "");
   fPlayer->AddInput(fc);


   Long64_t retval = -1;
   if (selector && strlen(selector)) {
      retval = Process(dset, selector, option, nentries, first);
   } else if (fSelector) {
      retval = Process(dset, fSelector, option, nentries, first);
   } else {
      Error("Process", "neither a selecrot file nor a selector object have"
                       " been specified: cannot process!");
   }
   fPlayer->GetInputList()->Remove(fc); // To avoid problems in future

   // Cleanup
   if (IsLite() && !fSync) {
      if (!fRunningDSets) fRunningDSets = new TList;
      fRunningDSets->Add(dset);
   } else {
      delete dset;
   }

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Process a dataset which is stored on the master with name 'dsetname'.
/// The syntax for dsetname is name[#[dir/]objname], e.g.
///   "mydset"       analysis of the first tree in the top dir of the dataset
///                  named "mydset"
///   "mydset#T"     analysis tree "T" in the top dir of the dataset
///                  named "mydset"
///   "mydset#adir/T" analysis tree "T" in the dir "adir" of the dataset
///                  named "mydset"
///   "mydset#adir/" analysis of the first tree in the dir "adir" of the
///                  dataset named "mydset"
/// The component 'name' in its more general form contains also the group and
/// user name following "/<group>/<user>/<dsname>". Each of these components
/// can contain one or more wildcards '*', in which case all the datasets matching
/// the expression are added together as a global dataset (wildcard support has
/// been added in version 5.27/02).
/// The last argument 'elist' specifies an entry- or event-list to be used as
/// event selection.
/// It is also possible (starting w/ version 5.27/02) to run on multiple datasets
/// at once in a more flexible way that the one provided by wildcarding. There
/// are three possibilities:
///    1) specifying the dataset names separated by the OR operator '|', e.g.
///          dsetname = "<dset1>|<dset2>|<dset3>|..."
///       in this case the datasets are a seen as a global unique dataset
///    2) specifying the dataset names separated by a ',' or a ' ', e.g.
///          dsetname = "<dset1>,<dset2> <dset3>,..."
///       in this case the datasets are processed one after the other and the
///       selector is notified when switching dataset via a bit in the current
///       processed element.
///    3) giving the path of a textfile where the dataset names are specified
///       on one or multiple lines; the lines found are joined as in 1), unless
///       the filepath is followed by a ',' (i.e. p->Process("datasets.txt,",...)
///       with the dataset names listed in 'datasets.txt') in which case they are
///       treated as in 2); the file is open in raw mode with TFile::Open and
///       therefore it cane be remote, e.g. on a Web server.
/// Each <dsetj> has the format specified above for the single dataset processing,
/// included wildcarding (the name of the tree and subdirectory must be same for
/// all the datasets).
/// In the case of multiple datasets, 'elist' is treated a global entry list.
/// It is possible to specify per-dataset entry lists using the syntax
///   "mydset[#adir/[T]]?enl=entrylist"
/// or
///   "mydset[#adir/[T]]<<entrylist"
/// Here 'entrylist' is a tag identifying, in the order :
///   i. a named entry-list in the input list or in the input data list
///  ii. a named entry-list in memory (in gDirectory)
/// iii. the path of a file containing the entry-list to be used
/// In the case ii) and iii) the entry-list object(s) is(are) added to the input
/// data list.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.

Long64_t TProof::Process(const char *dsetname, const char *selector,
                         Option_t *option, Long64_t nentries,
                         Long64_t first, TObject *elist)
{
   if (fProtocol < 13) {
      Info("Process", "processing 'by name' not supported by the server");
      return -1;
   }

   TString dsname, fname(dsetname);
   // If the 'dsetname' corresponds to an existing and readable file we will try to
   // interpretate its content as names of datasets to be processed. One line can contain
   // more datasets, separated by ',' or '|'. By default the dataset lines will be added
   // (i.e. joined as in option '|'); if the file name ends with ',' the dataset lines are
   // joined with ','.
   const char *separator = (fname.EndsWith(",")) ? "," : "|";
   if (!strcmp(separator, ",") || fname.EndsWith("|")) fname.Remove(fname.Length()-1, 1);
   if (!(gSystem->AccessPathName(fname, kReadPermission))) {
      TUrl uf(fname, kTRUE);
      uf.SetOptions(TString::Format("%sfiletype=raw", uf.GetOptions()));
      TFile *f = TFile::Open(uf.GetUrl());
      if (f && !(f->IsZombie())) {
         const Int_t blen = 8192;
         char buf[blen];
         Long64_t rest = f->GetSize();
         while (rest > 0) {
            Long64_t len = (rest > blen - 1) ? blen - 1 : rest;
            if (f->ReadBuffer(buf, len)) {
               Error("Process", "problems reading from file '%s'", fname.Data());
               dsname = "";
               break;
            }
            buf[len] = '\0';
            dsname += buf;
            rest -= len;
         }
         f->Close();
         SafeDelete(f);
         // We fail if a failure occured
         if (rest > 0) return -1;
      } else {
         Error("Process", "could not open file '%s'", fname.Data());
         return -1;
      }
   }
   if (dsname.IsNull()) {
      dsname = dsetname;
   } else {
      // Remove trailing '\n'
      if (dsname.EndsWith("\n")) dsname.Remove(dsname.Length()-1, 1);
      // Replace all '\n' with the proper separator
      dsname.ReplaceAll("\n", separator);
      if (gDebug > 0) {
         Info("Process", "processing multi-dataset read from file '%s':", fname.Data());
         Info("Process", "  '%s'", dsname.Data());
      }
   }

   TString names(dsname), name, enl, newname;
   // If multi-dataset check if server supports it
   if (fProtocol < 28 && names.Index(TRegexp("[, |]")) != kNPOS) {
      Info("Process", "multi-dataset processing not supported by the server");
      return -1;
   }

   TEntryList *el = 0;
   TString dsobj, dsdir;
   Int_t from = 0;
   while (names.Tokenize(name, from, "[, |]")) {

      newname = name;
      // Extract the specific entry-list, if any
      enl = "";
      Int_t ienl = name.Index("?enl=");
      if (ienl == kNPOS) {
         ienl = name.Index("<<");
         if (ienl != kNPOS) {
            newname.Remove(ienl);
            ienl += strlen("<<");
         }
      } else {
         newname.Remove(ienl);
         ienl += strlen("?enl=");
      }

      // Check the name syntax first
      TString obj, dir("/");
      Int_t idxc = newname.Index("#");
      if (idxc != kNPOS) {
         Int_t idxs = newname.Index("/", 1, idxc, TString::kExact);
         if (idxs != kNPOS) {
            obj = newname(idxs+1, newname.Length());
            dir = newname(idxc+1, newname.Length());
            dir.Remove(dir.Index("/") + 1);
            newname.Remove(idxc);
         } else {
            obj = newname(idxc+1, newname.Length());
            newname.Remove(idxc);
         }
      } else if (newname.Index(":") != kNPOS && newname.Index("://") == kNPOS) {
         // protection against using ':' instead of '#'
         Error("Process", "bad name syntax (%s): please use"
                          " a '#' after the dataset name", name.Data());
         dsname.ReplaceAll(name, "");
         continue;
      }
      if (dsobj.IsNull() && dsdir.IsNull()) {
         // The first one specifies obj and dir
         dsobj = obj;
         dsdir = dir;
      } else if (obj != dsobj || dir != dsdir) {
         // Inconsistent specification: not supported
         Warning("Process", "'obj' or 'dir' specification not consistent w/ the first given: ignore");
      }
      // Process the entry-list name, if any
      if (ienl != kNPOS) {
         // Get entrylist name or path
         enl = name(ienl, name.Length());
         el = 0;
         TObject *oel = 0;
         // If not in the input list ...
         TList *inpl = GetInputList();
         if (inpl && (oel = inpl->FindObject(enl))) el = dynamic_cast<TEntryList *>(oel);
         // ... check the heap
         if (!el && gDirectory && (oel = gDirectory->FindObject(enl))) {
            if ((el = dynamic_cast<TEntryList *>(oel))) {
               // Add to the input list (input data not available on master where
               // this info will be processed)
               if (fProtocol >= 28)
                  if (!(inpl->FindObject(el->GetName()))) AddInput(el);
            }
         }
         // If not in the heap, check a file, if any
         if (!el) {
            if (!gSystem->AccessPathName(enl)) {
               TFile *f = TFile::Open(enl);
               if (f && !(f->IsZombie()) && f->GetListOfKeys()) {
                  TIter nxk(f->GetListOfKeys());
                  TKey *k = 0;
                  while ((k = (TKey *) nxk())) {
                     if (!strcmp(k->GetClassName(), "TEntryList")) {
                        if (!el) {
                           if ((el = dynamic_cast<TEntryList *>(f->Get(k->GetName())))) {
                              // Add to the input list (input data not available on master where
                              // this info will be processed)
                              if (fProtocol >= 28) {
                                 if (!(inpl->FindObject(el->GetName()))) {
                                    el = (TEntryList *) el->Clone();
                                    AddInput(el);
                                 }
                              } else {
                                 el = (TEntryList *) el->Clone();
                              }
                           }
                        } else if (strcmp(el->GetName(), k->GetName())) {
                           Warning("Process", "multiple entry lists found in file '%s': the first one is taken;\n"
                                              "if this is not what you want, load first the content in memory"
                                              "and select it by name  ", enl.Data());
                        }
                     }
                  }
               } else {
                  Warning("Process","file '%s' cannot be open or is empty - ignoring", enl.Data());
               }
            }
         }
         // Transmit the information
         if (fProtocol >= 28) {
            newname += "?enl=";
            if (el) {
               // An entry list object is avalaible in the input list: add its name
               newname += el->GetName();
            } else {
               // The entry list object was not found: send the name, the future entry list manager will
               // find it on the server side
               newname += enl;
            }
         }
      }
      // Adjust the name for this dataset
      dsname.ReplaceAll(name, newname);
   }

   // Create the dataset object
   TDSet *dset = new TDSet(dsname, dsobj, dsdir);
   // Set entry list
   if (el && fProtocol < 28) {
      dset->SetEntryList(el);
   } else {
      dset->SetEntryList(elist);
   }
   // Run
   Long64_t retval = -1;
   if (selector && strlen(selector)) {
      retval = Process(dset, selector, option, nentries, first);
   } else if (fSelector) {
      retval = Process(dset, fSelector, option, nentries, first);
   } else {
      Error("Process", "neither a selector file nor a selector object have"
                       " been specified: cannot process!");
   }
   // Cleanup
   if (IsLite() && !fSync) {
      if (!fRunningDSets) fRunningDSets = new TList;
      fRunningDSets->Add(dset);
   } else {
      delete dset;
   }

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Generic (non-data based) selector processing: the Process() method of the
/// specified selector (.C) or TSelector object is called 'n' times.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.

Long64_t TProof::Process(const char *selector, Long64_t n, Option_t *option)
{
   if (!IsValid()) return -1;

   if (fProtocol < 16) {
      Info("Process", "server version < 5.17/04: generic processing not supported");
      return -1;
   }

   // Fake data set
   TDSet *dset = new TDSet;
   dset->SetBit(TDSet::kEmpty);

   Long64_t retval = -1;
   if (selector && strlen(selector)) {
      retval = Process(dset, selector, option, n);
   } else if (fSelector) {
      retval = Process(dset, fSelector, option, n);
   } else {
      Error("Process", "neither a selector file nor a selector object have"
                       " been specified: cannot process!");
   }

   // Cleanup
   if (IsLite() && !fSync) {
      if (!fRunningDSets) fRunningDSets = new TList;
      fRunningDSets->Add(dset);
   } else {
      delete dset;
   }
   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Process a data set (TDSet) using the specified selector object.
/// Entry- or event-lists should be set in the data set object using
/// TDSet::SetEntryList.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.

Long64_t TProof::Process(TDSet *dset, TSelector *selector, Option_t *option,
                         Long64_t nentries, Long64_t first)
{
   if (fProtocol < 34) {
      Error("Process", "server version < 5.33/02:"
                       "processing by object not supported");
      return -1;
   }
   if (!selector) {
      Error("Process", "selector object undefined!");
      return -1;
   }
   fSelector = selector;
   Long64_t rc = Process(dset, (const char*)0, option, nentries, first);
   fSelector = 0;
   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Process a data set (TFileCollection) using the specified selector object
/// The default tree is analyzed (i.e. the first one found). To specify another
/// tree, the default tree can be changed using TFileCollection::SetDefaultMetaData .
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.

Long64_t TProof::Process(TFileCollection *fc, TSelector *selector,
                         Option_t *option, Long64_t nentries, Long64_t first)
{
   if (fProtocol < 34) {
      Error("Process", "server version < 5.33/02:"
                       "processing by object not supported");
      return -1;
   }
   if (!selector) {
      Error("Process", "selector object undefined!");
      return -1;
   }
   fSelector = selector;
   Long64_t rc = Process(fc, (const char*)0, option, nentries, first);
   fSelector = 0;
   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Process with name of dataset and TSelector object

Long64_t TProof::Process(const char *dsetname, TSelector *selector,
                         Option_t *option, Long64_t nentries,
                         Long64_t first, TObject *elist)
{
   if (fProtocol < 34) {
      Error("Process", "server version < 5.33/02:"
                       "processing by object not supported");
      return -1;
   }
   if (!selector) {
      Error("Process", "selector object undefined!");
      return -1;
   }
   fSelector = selector;
   Long64_t rc = Process(dsetname, (const char*)0, option, nentries, first, elist);
   fSelector = 0;
   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Generic (non-data based) selector processing: the Process() method of the
/// specified selector is called 'n' times.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.

Long64_t TProof::Process(TSelector *selector, Long64_t n, Option_t *option)
{
   if (fProtocol < 34) {
      Error("Process", "server version < 5.33/02:"
                       "processing by object not supported");
      return -1;
   }
   if (!selector) {
      Error("Process", "selector object undefined!");
      return -1;
   }
   fSelector = selector;
   Long64_t rc = Process((const char*)0, n, option);
   fSelector = 0;
   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Get reference for the qry-th query in fQueries (as
/// displayed by ShowQueries).

Int_t TProof::GetQueryReference(Int_t qry, TString &ref)
{
   ref = "";
   if (qry > 0) {
      if (!fQueries)
         GetListOfQueries();
      if (fQueries) {
         TIter nxq(fQueries);
         TQueryResult *qr = 0;
         while ((qr = (TQueryResult *) nxq()))
            if (qr->GetSeqNum() == qry) {
               ref.Form("%s:%s", qr->GetTitle(), qr->GetName());
               return 0;
            }
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Finalize the qry-th query in fQueries.
/// If force, force retrieval if the query is found in the local list
/// but has already been finalized (default kFALSE).
/// If query < 0, finalize current query.
/// Return 0 on success, -1 on error

Long64_t TProof::Finalize(Int_t qry, Bool_t force)
{
   if (fPlayer) {
      if (qry > 0) {
         TString ref;
         if (GetQueryReference(qry, ref) == 0) {
            return Finalize(ref, force);
         } else {
            Info("Finalize", "query #%d not found", qry);
         }
      } else {
         // The last query
         return Finalize("", force);
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Finalize query with reference ref.
/// If force, force retrieval if the query is found in the local list
/// but has already been finalized (default kFALSE).
/// If ref = 0, finalize current query.
/// Return 0 on success, -1 on error

Long64_t TProof::Finalize(const char *ref, Bool_t force)
{
   if (fPlayer) {
      // Get the pointer to the query
      TQueryResult *qr = (ref && strlen(ref) > 0) ? fPlayer->GetQueryResult(ref)
                                                  : GetQueryResult();
      Bool_t retrieve = kFALSE;
      TString xref(ref);
      if (!qr) {
         if (!xref.IsNull()) {
            retrieve =  kTRUE;
         }
      } else {
         if (qr->IsFinalized()) {
            if (force) {
               retrieve = kTRUE;
            } else {
               Info("Finalize","query already finalized:"
                     " use Finalize(<qry>,kTRUE) to force new retrieval");
               qr = 0;
            }
         } else {
            retrieve = kTRUE;
            xref.Form("%s:%s", qr->GetTitle(), qr->GetName());
         }
      }
      if (retrieve) {
         Retrieve(xref.Data());
         qr = fPlayer->GetQueryResult(xref.Data());
      }
      if (qr)
         return fPlayer->Finalize(qr);
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Send retrieve request for the qry-th query in fQueries.
/// If path is defined save it to path.

Int_t TProof::Retrieve(Int_t qry, const char *path)
{
   if (qry > 0) {
      TString ref;
      if (GetQueryReference(qry, ref) == 0)
         return Retrieve(ref, path);
      else
         Info("Retrieve", "query #%d not found", qry);
   } else {
      Info("Retrieve","positive argument required - do nothing");
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Send retrieve request for the query specified by ref.
/// If path is defined save it to path.
/// Generic method working for all queries known by the server.

Int_t TProof::Retrieve(const char *ref, const char *path)
{
   if (ref) {
      TMessage m(kPROOF_RETRIEVE);
      m << TString(ref);
      Broadcast(m, kActive);
      Collect(kActive, fCollectTimeout);

      // Archive it locally, if required
      if (path) {

         // Get pointer to query
         TQueryResult *qr = fPlayer ? fPlayer->GetQueryResult(ref) : 0;

         if (qr) {

            TFile *farc = TFile::Open(path,"UPDATE");
            if (!farc || (farc && !(farc->IsOpen()))) {
               Info("Retrieve", "archive file cannot be open (%s)", path);
               return 0;
            }
            farc->cd();

            // Update query status
            qr->SetArchived(path);

            // Write to file
            qr->Write();

            farc->Close();
            SafeDelete(farc);

         } else {
            Info("Retrieve", "query not found after retrieve");
            return -1;
         }
      }

      return 0;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Send remove request for the qry-th query in fQueries.

Int_t TProof::Remove(Int_t qry, Bool_t all)
{
   if (qry > 0) {
      TString ref;
      if (GetQueryReference(qry, ref) == 0)
         return Remove(ref, all);
      else
         Info("Remove", "query #%d not found", qry);
   } else {
      Info("Remove","positive argument required - do nothing");
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Send remove request for the query specified by ref.
/// If all = TRUE remove also local copies of the query, if any.
/// Generic method working for all queries known by the server.
/// This method can be also used to reset the list of queries
/// waiting to be processed: for that purpose use ref == "cleanupqueue".

Int_t TProof::Remove(const char *ref, Bool_t all)
{
   if (all) {
      // Remove also local copies, if any
      if (fPlayer)
         fPlayer->RemoveQueryResult(ref);
   }

   if (IsLite()) return 0;

   if (ref) {
      TMessage m(kPROOF_REMOVE);
      m << TString(ref);
      Broadcast(m, kActive);
      Collect(kActive, fCollectTimeout);
      return 0;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Send archive request for the qry-th query in fQueries.

Int_t TProof::Archive(Int_t qry, const char *path)
{
   if (qry > 0) {
      TString ref;
      if (GetQueryReference(qry, ref) == 0)
         return Archive(ref, path);
      else
         Info("Archive", "query #%d not found", qry);
   } else {
      Info("Archive","positive argument required - do nothing");
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Send archive request for the query specified by ref.
/// Generic method working for all queries known by the server.
/// If ref == "Default", path is understood as a default path for
/// archiving.

Int_t TProof::Archive(const char *ref, const char *path)
{
   if (ref) {
      TMessage m(kPROOF_ARCHIVE);
      m << TString(ref) << TString(path);
      Broadcast(m, kActive);
      Collect(kActive, fCollectTimeout);
      return 0;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Send cleanup request for the session specified by tag.

Int_t TProof::CleanupSession(const char *sessiontag)
{
   if (sessiontag) {
      TMessage m(kPROOF_CLEANUPSESSION);
      m << TString(sessiontag);
      Broadcast(m, kActive);
      Collect(kActive, fCollectTimeout);
      return 0;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Change query running mode to the one specified by 'mode'.

void TProof::SetQueryMode(EQueryMode mode)
{
   fQueryMode = mode;

   if (gDebug > 0)
      Info("SetQueryMode","query mode is set to: %s", fQueryMode == kSync ?
           "Sync" : "Async");
}

////////////////////////////////////////////////////////////////////////////////
/// Find out the query mode based on the current setting and 'mode'.

TProof::EQueryMode TProof::GetQueryMode(Option_t *mode) const
{
   EQueryMode qmode = fQueryMode;

   if (mode && (strlen(mode) > 0)) {
      TString m(mode);
      m.ToUpper();
      if (m.Contains("ASYN")) {
         qmode = kAsync;
      } else if (m.Contains("SYNC")) {
         qmode = kSync;
      }
   }

   if (gDebug > 0)
      Info("GetQueryMode","query mode is set to: %s", qmode == kSync ?
           "Sync" : "Async");

   return qmode;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute the specified drawing action on a data set (TDSet).
/// Event- or Entry-lists should be set in the data set object using
/// TDSet::SetEntryList.
/// Returns -1 in case of error or number of selected events otherwise.

Long64_t TProof::DrawSelect(TDSet *dset, const char *varexp,
                            const char *selection, Option_t *option,
                            Long64_t nentries, Long64_t first)
{
   if (!IsValid() || !fPlayer) return -1;

   // Make sure that asynchronous processing is not active
   if (!IsIdle()) {
      Info("DrawSelect","not idle, asynchronous Draw not supported");
      return -1;
   }
   TString opt(option);
   Int_t idx = opt.Index("ASYN", 0, TString::kIgnoreCase);
   if (idx != kNPOS)
      opt.Replace(idx,4,"");

   return fPlayer->DrawSelect(dset, varexp, selection, opt, nentries, first);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute the specified drawing action on a data set which is stored on the
/// master with name 'dsetname'.
/// The syntax for dsetname is name[#[dir/]objname], e.g.
///   "mydset"       analysis of the first tree in the top dir of the dataset
///                  named "mydset"
///   "mydset#T"     analysis tree "T" in the top dir of the dataset
///                  named "mydset"
///   "mydset#adir/T" analysis tree "T" in the dir "adir" of the dataset
///                  named "mydset"
///   "mydset#adir/" analysis of the first tree in the dir "adir" of the
///                  dataset named "mydset"
/// The last argument 'enl' specifies an entry- or event-list to be used as
/// event selection.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.

Long64_t TProof::DrawSelect(const char *dsetname, const char *varexp,
                            const char *selection, Option_t *option,
                            Long64_t nentries, Long64_t first, TObject *enl)
{
   if (fProtocol < 13) {
      Info("Process", "processing 'by name' not supported by the server");
      return -1;
   }

   TString name(dsetname);
   TString obj;
   TString dir = "/";
   Int_t idxc = name.Index("#");
   if (idxc != kNPOS) {
      Int_t idxs = name.Index("/", 1, idxc, TString::kExact);
      if (idxs != kNPOS) {
         obj = name(idxs+1, name.Length());
         dir = name(idxc+1, name.Length());
         dir.Remove(dir.Index("/") + 1);
         name.Remove(idxc);
      } else {
         obj = name(idxc+1, name.Length());
         name.Remove(idxc);
      }
   } else if (name.Index(":") != kNPOS && name.Index("://") == kNPOS) {
      // protection against using ':' instead of '#'
      Error("DrawSelect", "bad name syntax (%s): please use"
                       " a '#' after the dataset name", dsetname);
      return -1;
   }

   TDSet *dset = new TDSet(name, obj, dir);
   // Set entry-list, if required
   dset->SetEntryList(enl);
   Long64_t retval = DrawSelect(dset, varexp, selection, option, nentries, first);
   delete dset;
   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Send STOPPROCESS message to master and workers.

void TProof::StopProcess(Bool_t abort, Int_t timeout)
{
   PDB(kGlobal,2)
      Info("StopProcess","enter %d", abort);

   if (!IsValid())
      return;

   // Flag that we have been stopped
   ERunStatus rst = abort ? TProof::kAborted : TProof::kStopped;
   SetRunStatus(rst);

   if (fPlayer)
      fPlayer->StopProcess(abort, timeout);

   // Stop any blocking 'Collect' request; on masters we do this only if
   // aborting; when stopping, we still need to receive the results
   if (TestBit(TProof::kIsClient) || abort)
      InterruptCurrentMonitor();

   if (fSlaves->GetSize() == 0)
      return;

   // Notify the remote counterpart
   TSlave *sl;
   TIter   next(fSlaves);
   while ((sl = (TSlave *)next()))
      if (sl->IsValid())
         // Ask slave to progate the stop/abort request
         sl->StopProcess(abort, timeout);
}

////////////////////////////////////////////////////////////////////////////////
/// Signal to disable related switches

void TProof::DisableGoAsyn()
{
   Emit("DisableGoAsyn()");
}

////////////////////////////////////////////////////////////////////////////////
/// Send GOASYNC message to the master.

void TProof::GoAsynchronous()
{
   if (!IsValid()) return;

   if (GetRemoteProtocol() < 22) {
      Info("GoAsynchronous", "functionality not supported by the server - ignoring");
      return;
   }

   if (fSync && !IsIdle()) {
      TMessage m(kPROOF_GOASYNC);
      Broadcast(m);
   } else {
      Info("GoAsynchronous", "either idle or already in asynchronous mode - ignoring");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Receive the log file of the slave with socket s.

void TProof::RecvLogFile(TSocket *s, Int_t size)
{
   const Int_t kMAXBUF = 16384;  //32768  //16384  //65536;
   char buf[kMAXBUF];

   // If macro saving is enabled prepare macro
   if (fSaveLogToMacro && fMacroLog.GetListOfLines()) {
      fMacroLog.GetListOfLines()->SetOwner(kTRUE);
      fMacroLog.GetListOfLines()->Clear();
   }

   // Append messages to active logging unit
   Int_t fdout = -1;
   if (!fLogToWindowOnly) {
      fdout = (fRedirLog) ? fileno(fLogFileW) : fileno(stdout);
      if (fdout < 0) {
         Warning("RecvLogFile", "file descriptor for outputs undefined (%d):"
                 " will not log msgs", fdout);
         return;
      }
      lseek(fdout, (off_t) 0, SEEK_END);
   }

   Int_t  left, rec, r;
   Long_t filesize = 0;

   while (filesize < size) {
      left = Int_t(size - filesize);
      if (left >= kMAXBUF)
         left = kMAXBUF-1;
      rec = s->RecvRaw(&buf, left);
      filesize = (rec > 0) ? (filesize + rec) : filesize;
      if (!fLogToWindowOnly && !fSaveLogToMacro) {
         if (rec > 0) {

            char *p = buf;
            r = rec;
            while (r) {
               Int_t w;

               w = write(fdout, p, r);

               if (w < 0) {
                  SysError("RecvLogFile", "error writing to unit: %d", fdout);
                  break;
               }
               r -= w;
               p += w;
            }
         } else if (rec < 0) {
            Error("RecvLogFile", "error during receiving log file");
            break;
         }
      }
      if (rec > 0) {
         buf[rec] = 0;
         EmitVA("LogMessage(const char*,Bool_t)", 2, buf, kFALSE);
         // If macro saving is enabled add to TMacro
         if (fSaveLogToMacro) fMacroLog.AddLine(buf);
      }
   }

   // If idle restore logs to main session window
   if (fRedirLog && IsIdle() && !TestBit(TProof::kIsMaster))
      fRedirLog = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Notify locally 'msg' to the appropriate units (file, stdout, window)
/// If defined, 'sfx' is added after 'msg' (typically a line-feed);

void TProof::NotifyLogMsg(const char *msg, const char *sfx)
{
   // Must have somenthing to notify
   Int_t len = 0;
   if (!msg || (len = strlen(msg)) <= 0)
      return;

   // Get suffix length if any
   Int_t lsfx = (sfx) ? strlen(sfx) : 0;

   // Append messages to active logging unit
   Int_t fdout = -1;
   if (!fLogToWindowOnly) {
      fdout = (fRedirLog) ? fileno(fLogFileW) : fileno(stdout);
      if (fdout < 0) {
         Warning("NotifyLogMsg", "file descriptor for outputs undefined (%d):"
                 " will not notify msgs", fdout);
         return;
      }
      lseek(fdout, (off_t) 0, SEEK_END);
   }

   if (!fLogToWindowOnly) {
      // Write to output unit (stdout or a log file)
      if (len > 0) {
         char *p = (char *)msg;
         Int_t r = len;
         while (r) {
            Int_t w = write(fdout, p, r);
            if (w < 0) {
               SysError("NotifyLogMsg", "error writing to unit: %d", fdout);
               break;
            }
            r -= w;
            p += w;
         }
         // Add a suffix, if requested
         if (lsfx > 0)
            if (write(fdout, sfx, lsfx) != lsfx)
               SysError("NotifyLogMsg", "error writing to unit: %d", fdout);
      }
   }
   if (len > 0) {
      // Publish the message to the separate window (if the latter is missing
      // the message will just get lost)
      EmitVA("LogMessage(const char*,Bool_t)", 2, msg, kFALSE);
   }

   // If idle restore logs to main session window
   if (fRedirLog && IsIdle())
      fRedirLog = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Log a message into the appropriate window by emitting a signal.

void TProof::LogMessage(const char *msg, Bool_t all)
{
   PDB(kGlobal,1)
      Info("LogMessage","Enter ... %s, 'all: %s", msg ? msg : "",
           all ? "true" : "false");

   if (gROOT->IsBatch()) {
      PDB(kGlobal,1) Info("LogMessage","GUI not started - use TProof::ShowLog()");
      return;
   }

   if (msg)
      EmitVA("LogMessage(const char*,Bool_t)", 2, msg, all);

   // Re-position at the beginning of the file, if requested.
   // This is used by the dialog when it re-opens the log window to
   // provide all the session messages
   if (all)
      lseek(fileno(fLogFileR), (off_t) 0, SEEK_SET);

   const Int_t kMAXBUF = 32768;
   char buf[kMAXBUF];
   Int_t len;
   do {
      while ((len = read(fileno(fLogFileR), buf, kMAXBUF-1)) < 0 &&
             TSystem::GetErrno() == EINTR)
         TSystem::ResetErrno();

      if (len < 0) {
         Error("LogMessage", "error reading log file");
         break;
      }

      if (len > 0) {
         buf[len] = 0;
         EmitVA("LogMessage(const char*,Bool_t)", 2, buf, kFALSE);
      }

   } while (len > 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Send to all active slaves servers the current slave group size
/// and their unique id. Returns number of active slaves.
/// Returns -1 in case of error.

Int_t TProof::SendGroupView()
{
   if (!IsValid()) return -1;
   if (TestBit(TProof::kIsClient)) return 0;
   if (!fSendGroupView) return 0;
   fSendGroupView = kFALSE;

   TIter   next(fActiveSlaves);
   TSlave *sl;

   int  bad = 0, cnt = 0, size = GetNumberOfActiveSlaves();
   char str[32];

   while ((sl = (TSlave *)next())) {
      snprintf(str, 32, "%d %d", cnt, size);
      if (sl->GetSocket()->Send(str, kPROOF_GROUPVIEW) == -1) {
         MarkBad(sl, "could not send kPROOF_GROUPVIEW message");
         bad++;
      } else
         cnt++;
   }

   // Send the group view again in case there was a change in the
   // group size due to a bad slave

   if (bad) SendGroupView();

   return GetNumberOfActiveSlaves();
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to extract the filename (if any) form a CINT command.
/// Returns kTRUE and the filename in 'fn'; returns kFALSE if not found or not
/// appliable.

Bool_t TProof::GetFileInCmd(const char *cmd, TString &fn)
{
   TString s = cmd;
   s = s.Strip(TString::kBoth);

   if (s.Length() > 0 &&
      (s.BeginsWith(".L") || s.BeginsWith(".x") || s.BeginsWith(".X"))) {
      TString file = s(2, s.Length());
      TString acm, arg, io;
      fn = gSystem->SplitAclicMode(file, acm, arg, io);
      if (!fn.IsNull())
         return kTRUE;
   }

   // Not found
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Send command to be executed on the PROOF master and/or slaves.
/// If plusMaster is kTRUE then exeucte on slaves and master too.
/// Command can be any legal command line command. Commands like
/// ".x file.C" or ".L file.C" will cause the file file.C to be send
/// to the PROOF cluster. Returns -1 in case of error, >=0 in case of
/// succes.

Int_t TProof::Exec(const char *cmd, Bool_t plusMaster)
{
   return Exec(cmd, kActive, plusMaster);
}

////////////////////////////////////////////////////////////////////////////////
/// Send command to be executed on the PROOF master and/or slaves.
/// Command can be any legal command line command. Commands like
/// ".x file.C" or ".L file.C" will cause the file file.C to be send
/// to the PROOF cluster. Returns -1 in case of error, >=0 in case of
/// succes.

Int_t TProof::Exec(const char *cmd, ESlaves list, Bool_t plusMaster)
{
   if (!IsValid()) return -1;

   TString s = cmd;
   s = s.Strip(TString::kBoth);

   if (!s.Length()) return 0;

   // check for macro file and make sure the file is available on all slaves
   TString filename;
   if (TProof::GetFileInCmd(s.Data(), filename)) {
      char *fn = gSystem->Which(TROOT::GetMacroPath(), filename, kReadPermission);
      if (fn) {
         if (GetNumberOfUniqueSlaves() > 0) {
            if (SendFile(fn, kAscii | kForward | kCpBin) < 0) {
               Error("Exec", "file %s could not be transfered", fn);
               delete [] fn;
               return -1;
            }
         } else {
            TString scmd = s(0,3) + fn;
            Int_t n = SendCommand(scmd, list);
            delete [] fn;
            return n;
         }
      } else {
         Error("Exec", "macro %s not found", filename.Data());
         return -1;
      }
      delete [] fn;
   }

   if (plusMaster) {
      if (IsLite()) {
         gROOT->ProcessLine(cmd);
      } else {
         DeactivateWorker("*");
         Int_t res = SendCommand(cmd, list);
         ActivateWorker("restore");
         if (res < 0)
            return res;
      }
   }
   return SendCommand(cmd, list);
}

////////////////////////////////////////////////////////////////////////////////
/// Send command to be executed on node of ordinal 'ord' (use "0" for master).
/// Command can be any legal command line command. Commands like
/// ".x file.C" or ".L file.C" will cause the file file.C to be send
/// to the PROOF cluster.
/// If logtomacro is TRUE the text result of the action is saved in the fMacroLog
/// TMacro, accessible via TMacro::GetMacroLog();
/// Returns -1 in case of error, >=0 in case of succes.

Int_t TProof::Exec(const char *cmd, const char *ord, Bool_t logtomacro)
{
   if (!IsValid()) return -1;

   TString s = cmd;
   s = s.Strip(TString::kBoth);

   if (!s.Length()) return 0;

   Int_t res = 0;
   if (IsLite()) {
      gROOT->ProcessLine(cmd);
   } else {
      Bool_t oldRedirLog = fRedirLog;
      fRedirLog = kTRUE;
      // Deactivate all workers
      DeactivateWorker("*");
      fRedirLog = kFALSE;
      // Reactivate the target ones, if needed
      if (strcmp(ord, "master") && strcmp(ord, "0")) ActivateWorker(ord);
      // Honour log-to-macro-saving settings
      Bool_t oldSaveLog = fSaveLogToMacro;
      fSaveLogToMacro = logtomacro;
      res = SendCommand(cmd, kActive);
      fSaveLogToMacro = oldSaveLog;
      fRedirLog = kTRUE;
      ActivateWorker("restore");
      fRedirLog = oldRedirLog;
   }
   // Done
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Send command to be executed on the PROOF master and/or slaves.
/// Command can be any legal command line command, however commands
/// like ".x file.C" or ".L file.C" will not cause the file.C to be
/// transfered to the PROOF cluster. In that case use TProof::Exec().
/// Returns the status send by the remote server as part of the
/// kPROOF_LOGDONE message. Typically this is the return code of the
/// command on the remote side. Returns -1 in case of error.

Int_t TProof::SendCommand(const char *cmd, ESlaves list)
{
   if (!IsValid()) return -1;

   Broadcast(cmd, kMESS_CINT, list);
   Collect(list);

   return fStatus;
}

////////////////////////////////////////////////////////////////////////////////
/// Get value of environment variable 'env' on node 'ord'

TString TProof::Getenv(const char *env, const char *ord)
{
   // The command to be executed
   TString cmd = TString::Format("gSystem->Getenv(\"%s\")", env);
   if (Exec(cmd.Data(), ord, kTRUE) != 0) return TString("");
   // Get the line
   TObjString *os = fMacroLog.GetLineWith("const char");
   if (os) {
      TString info;
      Ssiz_t from = 0;
      os->GetString().Tokenize(info, from, "\"");
      os->GetString().Tokenize(info, from, "\"");
      if (gDebug > 0) Printf("%s: '%s'", env, info.Data());
      return info;
   }
   return TString("");
}

////////////////////////////////////////////////////////////////////////////////
/// Get into 'env' the value of integer RC env variable 'rcenv' on node 'ord'

Int_t TProof::GetRC(const char *rcenv, Int_t &env, const char *ord)
{
   // The command to be executed
   TString cmd = TString::Format("if (gEnv->Lookup(\"%s\")) { gEnv->GetValue(\"%s\",\"\"); }", rcenv, rcenv);
   // Exectute the command saving the logs to macro
   if (Exec(cmd.Data(), ord, kTRUE) != 0) return -1;
   // Get the line
   TObjString *os = fMacroLog.GetLineWith("const char");
   Int_t rc = -1;
   if (os) {
      Ssiz_t fst =  os->GetString().First('\"');
      Ssiz_t lst =  os->GetString().Last('\"');
      TString info = os->GetString()(fst+1, lst-fst-1);
      if (info.IsDigit()) {
         env = info.Atoi();
         rc = 0;
         if (gDebug > 0)
            Printf("%s: %d", rcenv, env);
      }
   }
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Get into 'env' the value of double RC env variable 'rcenv' on node 'ord'

Int_t TProof::GetRC(const char *rcenv, Double_t &env, const char *ord)
{
   // The command to be executed
   TString cmd = TString::Format("if (gEnv->Lookup(\"%s\")) { gEnv->GetValue(\"%s\",\"\"); }", rcenv, rcenv);
   // Exectute the command saving the logs to macro
   if (Exec(cmd.Data(), ord, kTRUE) != 0) return -1;
   // Get the line
   TObjString *os = fMacroLog.GetLineWith("const char");
   Int_t rc = -1;
   if (os) {
      Ssiz_t fst =  os->GetString().First('\"');
      Ssiz_t lst =  os->GetString().Last('\"');
      TString info = os->GetString()(fst+1, lst-fst-1);
      if (info.IsFloat()) {
         env = info.Atof();
         rc = 0;
         if (gDebug > 0)
            Printf("%s: %f", rcenv, env);
      }
   }
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Get into 'env' the value of string RC env variable 'rcenv' on node 'ord'

Int_t TProof::GetRC(const char *rcenv, TString &env, const char *ord)
{
   // The command to be executed
   TString cmd = TString::Format("if (gEnv->Lookup(\"%s\")) { gEnv->GetValue(\"%s\",\"\"); }", rcenv, rcenv);
   // Exectute the command saving the logs to macro
   if (Exec(cmd.Data(), ord, kTRUE) != 0) return -1;
   // Get the line
   TObjString *os = fMacroLog.GetLineWith("const char");
   Int_t rc = -1;
   if (os) {
      Ssiz_t fst =  os->GetString().First('\"');
      Ssiz_t lst =  os->GetString().Last('\"');
      env = os->GetString()(fst+1, lst-fst-1);
      rc = 0;
      if (gDebug > 0)
         Printf("%s: %s", rcenv, env.Data());
   }
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Transfer the current state of the master to the active slave servers.
/// The current state includes: the current working directory, etc.
/// Returns the number of active slaves. Returns -1 in case of error.

Int_t TProof::SendCurrentState(TList *list)
{
   if (!IsValid()) return -1;

   // Go to the new directory, reset the interpreter environment and
   // tell slave to delete all objects from its new current directory.
   Broadcast(gDirectory->GetPath(), kPROOF_RESET, list);

   return GetParallel();
}

////////////////////////////////////////////////////////////////////////////////
/// Transfer the current state of the master to the active slave servers.
/// The current state includes: the current working directory, etc.
/// Returns the number of active slaves. Returns -1 in case of error.

Int_t TProof::SendCurrentState(ESlaves list)
{
   if (!IsValid()) return -1;

   // Go to the new directory, reset the interpreter environment and
   // tell slave to delete all objects from its new current directory.
   Broadcast(gDirectory->GetPath(), kPROOF_RESET, list);

   return GetParallel();
}

////////////////////////////////////////////////////////////////////////////////
/// Transfer the initial (i.e. current) state of the master to all
/// slave servers. Currently the initial state includes: log level.
/// Returns the number of active slaves. Returns -1 in case of error.

Int_t TProof::SendInitialState()
{
   if (!IsValid()) return -1;

   SetLogLevel(fLogLevel, gProofDebugMask);

   return GetNumberOfActiveSlaves();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if a file needs to be send to the slave. Use the following
/// algorithm:
///   - check if file appears in file map
///     - if yes, get file's modtime and check against time in map,
///       if modtime not same get md5 and compare against md5 in map,
///       if not same return kTRUE.
///     - if no, get file's md5 and modtime and store in file map, ask
///       slave if file exists with specific md5, if yes return kFALSE,
///       if no return kTRUE.
/// The options 'cpopt' define if to copy things from cache to sandbox and what.
/// To retrieve from the cache the binaries associated with the file TProof::kCpBin
/// must be set in cpopt; the default is copy everything.
/// Returns kTRUE in case file needs to be send, returns kFALSE in case
/// file is already on remote node.

Bool_t TProof::CheckFile(const char *file, TSlave *slave, Long_t modtime, Int_t cpopt)
{
   Bool_t sendto = kFALSE;

   // create worker based filename
   TString sn = slave->GetName();
   sn += ":";
   sn += slave->GetOrdinal();
   sn += ":";
   sn += gSystem->BaseName(file);

   // check if file is in map
   FileMap_t::const_iterator it;
   if ((it = fFileMap.find(sn)) != fFileMap.end()) {
      // file in map
      MD5Mod_t md = (*it).second;
      if (md.fModtime != modtime) {
         TMD5 *md5 = TMD5::FileChecksum(file);
         if (md5) {
            if ((*md5) != md.fMD5) {
               sendto       = kTRUE;
               md.fMD5      = *md5;
               md.fModtime  = modtime;
               fFileMap[sn] = md;
               // When on the master, the master and/or slaves may share
               // their file systems and cache. Therefore always make a
               // check for the file. If the file already exists with the
               // expected md5 the kPROOF_CHECKFILE command will cause the
               // file to be copied from cache to slave sandbox.
               if (TestBit(TProof::kIsMaster)) {
                  sendto = kFALSE;
                  TMessage mess(kPROOF_CHECKFILE);
                  mess << TString(gSystem->BaseName(file)) << md.fMD5 << cpopt;
                  slave->GetSocket()->Send(mess);

                  fCheckFileStatus = 0;
                  Collect(slave, fCollectTimeout, kPROOF_CHECKFILE);
                  sendto = (fCheckFileStatus == 0) ? kTRUE : kFALSE;
               }
            }
            delete md5;
         } else {
            Error("CheckFile", "could not calculate local MD5 check sum - dont send");
            return kFALSE;
         }
      }
   } else {
      // file not in map
      TMD5 *md5 = TMD5::FileChecksum(file);
      MD5Mod_t md;
      if (md5) {
         md.fMD5      = *md5;
         md.fModtime  = modtime;
         fFileMap[sn] = md;
         delete md5;
      } else {
         Error("CheckFile", "could not calculate local MD5 check sum - dont send");
         return kFALSE;
      }
      TMessage mess(kPROOF_CHECKFILE);
      mess << TString(gSystem->BaseName(file)) << md.fMD5 << cpopt;
      slave->GetSocket()->Send(mess);

      fCheckFileStatus = 0;
      Collect(slave, fCollectTimeout, kPROOF_CHECKFILE);
      sendto = (fCheckFileStatus == 0) ? kTRUE : kFALSE;
   }

   return sendto;
}

////////////////////////////////////////////////////////////////////////////////
/// Send a file to master or slave servers. Returns number of slaves
/// the file was sent to, maybe 0 in case master and slaves have the same
/// file system image, -1 in case of error.
/// If defined, send to worker 'wrk' only.
/// If defined, the full path of the remote path will be rfile.
/// If rfile = "cache" the file is copied to the remote cache instead of the sandbox
/// (to copy to the cache on a different name use rfile = "cache:newname").
/// The mask 'opt' is an or of ESendFileOpt:
///
///       kAscii  (0x0)      if set true ascii file transfer is used
///       kBinary (0x1)      if set true binary file transfer is used
///       kForce  (0x2)      if not set an attempt is done to find out
///                          whether the file really needs to be downloaded
///                          (a valid copy may already exist in the cache
///                          from a previous run); the bit is set by
///                          UploadPackage, since the check is done elsewhere.
///       kForward (0x4)     if set, ask server to forward the file to slave
///                          or submaster (meaningless for slave servers).
///       kCpBin   (0x8)     Retrieve from the cache the binaries associated
///                          with the file
///       kCp      (0x10)    Retrieve the files from the cache
///

Int_t TProof::SendFile(const char *file, Int_t opt, const char *rfile, TSlave *wrk)
{
   if (!IsValid()) return -1;

   // Use the active slaves list ...
   TList *slaves = (rfile && !strcmp(rfile, "cache")) ? fUniqueSlaves : fActiveSlaves;
   // ... or the specified slave, if any
   if (wrk) {
      slaves = new TList();
      slaves->Add(wrk);
   }

   if (slaves->GetSize() == 0) return 0;

#ifndef R__WIN32
   Int_t fd = open(file, O_RDONLY);
#else
   Int_t fd = open(file, O_RDONLY | O_BINARY);
#endif
   if (fd < 0) {
      SysError("SendFile", "cannot open file %s", file);
      return -1;
   }

   // Get info about the file
   Long64_t size = -1;
   Long_t id, flags, modtime = 0;
   if (gSystem->GetPathInfo(file, &id, &size, &flags, &modtime) == 1) {
      Error("SendFile", "cannot stat file %s", file);
      close(fd);
      return -1;
   }
   if (size == 0) {
      Error("SendFile", "empty file %s", file);
      close(fd);
      return -1;
   }

   // Decode options
   Bool_t bin   = (opt & kBinary)  ? kTRUE : kFALSE;
   Bool_t force = (opt & kForce)   ? kTRUE : kFALSE;
   Bool_t fw    = (opt & kForward) ? kTRUE : kFALSE;

   // Copy options
   Int_t cpopt  = 0;
   if ((opt & kCp)) cpopt |= kCp;
   if ((opt & kCpBin)) cpopt |= (kCp | kCpBin);

   const Int_t kMAXBUF = 32768;  //16384  //65536;
   char buf[kMAXBUF];
   Int_t nsl = 0;

   TIter next(slaves);
   TSlave *sl;
   TString fnam(rfile);
   if (fnam == "cache") {
      fnam += TString::Format(":%s", gSystem->BaseName(file));
   } else if (fnam.IsNull()) {
      fnam = gSystem->BaseName(file);
   }
   // List on which we will collect the results
   fStatus = 0;
   while ((sl = (TSlave *)next())) {
      if (!sl->IsValid())
         continue;

      Bool_t sendto = force ? kTRUE : CheckFile(file, sl, modtime, cpopt);
      // Don't send the kPROOF_SENDFILE command to real slaves when sendto
      // is false. Masters might still need to send the file to newly added
      // slaves.
      PDB(kPackage,2) {
         const char *snd = (sl->fSlaveType == TSlave::kSlave && sendto) ? "" : "not";
         Info("SendFile", "%s sending file %s to: %s:%s (%d)", snd,
                    file, sl->GetName(), sl->GetOrdinal(), sendto);
      }
      if (sl->fSlaveType == TSlave::kSlave && !sendto)
         continue;
      // The value of 'size' is used as flag remotely, so we need to
      // reset it to 0 if we are not going to send the file
      Long64_t siz = sendto ? size : 0;
      snprintf(buf, kMAXBUF, "%s %d %lld %d", fnam.Data(), bin, siz, fw);
      if (sl->GetSocket()->Send(buf, kPROOF_SENDFILE) == -1) {
         MarkBad(sl, "could not send kPROOF_SENDFILE request");
         continue;
      }

      if (sendto) {

         lseek(fd, 0, SEEK_SET);

         Int_t len;
         do {
            while ((len = read(fd, buf, kMAXBUF)) < 0 && TSystem::GetErrno() == EINTR)
               TSystem::ResetErrno();

            if (len < 0) {
               SysError("SendFile", "error reading from file %s", file);
               Interrupt(kSoftInterrupt, kActive);
               close(fd);
               return -1;
            }

            if (len > 0 && sl->GetSocket()->SendRaw(buf, len) == -1) {
               SysError("SendFile", "error writing to slave %s:%s (now offline)",
                        sl->GetName(), sl->GetOrdinal());
               MarkBad(sl, "sendraw failure");
               sl = 0;
               break;
            }

         } while (len > 0);

         nsl++;
      }
      // Wait for the operation to be done
      if (sl)
         Collect(sl, fCollectTimeout, kPROOF_SENDFILE);
   }

   close(fd);

   // Cleanup temporary list, if any
   if (slaves != fActiveSlaves && slaves != fUniqueSlaves)
      SafeDelete(slaves);

   // We return failure is at least one unique worker failed
   return (fStatus != 0) ? -1 : nsl;
}

////////////////////////////////////////////////////////////////////////////////
/// Sends an object to master and workers and expect them to send back a
/// message with the output of its TObject::Print(). Returns -1 on error, the
/// number of workers that received the objects on success.

Int_t TProof::Echo(const TObject *obj)
{
   if (!IsValid() || !obj) return -1;
   TMessage mess(kPROOF_ECHO);
   mess.WriteObject(obj);
   return Broadcast(mess);
}

////////////////////////////////////////////////////////////////////////////////
/// Sends a string to master and workers and expect them to echo it back to
/// the client via a message. It is a special case of the generic Echo()
/// that works with TObjects. Returns -1 on error, the number of workers that
/// received the message on success.

Int_t TProof::Echo(const char *str)
{
   TObjString *os = new TObjString(str);
   Int_t rv = Echo(os);
   delete os;
   return rv;
}

////////////////////////////////////////////////////////////////////////////////
/// Send object to master or slave servers. Returns number of slaves object
/// was sent to, -1 in case of error.

Int_t TProof::SendObject(const TObject *obj, ESlaves list)
{
   if (!IsValid() || !obj) return -1;

   TMessage mess(kMESS_OBJECT);

   mess.WriteObject(obj);
   return Broadcast(mess, list);
}

////////////////////////////////////////////////////////////////////////////////
/// Send print command to master server. Returns number of slaves message
/// was sent to. Returns -1 in case of error.

Int_t TProof::SendPrint(Option_t *option)
{
   if (!IsValid()) return -1;

   Broadcast(option, kPROOF_PRINT, kActive);
   return Collect(kActive, fCollectTimeout);
}

////////////////////////////////////////////////////////////////////////////////
/// Set server logging level.

void TProof::SetLogLevel(Int_t level, UInt_t mask)
{
   char str[32];
   fLogLevel        = level;
   gProofDebugLevel = level;
   gProofDebugMask  = (TProofDebug::EProofDebugMask) mask;
   snprintf(str, 32, "%d %u", level, mask);
   Broadcast(str, kPROOF_LOGLEVEL, kAll);
}

////////////////////////////////////////////////////////////////////////////////
/// Switch ON/OFF the real-time logging facility. When this option is
/// ON, log messages from processing are sent back as they come, instead of
/// being sent back at the end in one go. This may help debugging or monitoring
/// in some cases, but, depending on the amount of log, it may have significant
/// consequencies on the load over the network, so it must be used with care.

void TProof::SetRealTimeLog(Bool_t on)
{
   if (IsValid()) {
      TMessage mess(kPROOF_REALTIMELOG);
      mess << on;
      Broadcast(mess);
   } else {
      Warning("SetRealTimeLog","session is invalid - do nothing");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Tell PROOF how many slaves to use in parallel. If random is TRUE a random
/// selection is done (if nodes is less than the available nodes).
/// Returns the number of parallel slaves. Returns -1 in case of error.

Int_t TProof::SetParallelSilent(Int_t nodes, Bool_t random)
{
   if (!IsValid()) return -1;

   if (TestBit(TProof::kIsMaster)) {
      if (!fDynamicStartup) GoParallel(nodes, kFALSE, random);
      return SendCurrentState();
   } else {
      if (nodes < 0) {
         PDB(kGlobal,1) Info("SetParallelSilent", "request all nodes");
      } else {
         PDB(kGlobal,1) Info("SetParallelSilent", "request %d node%s", nodes,
                                                  nodes == 1 ? "" : "s");
      }
      TMessage mess(kPROOF_PARALLEL);
      mess << nodes << random;
      Broadcast(mess);
      Collect(kActive, fCollectTimeout);
      Int_t n = GetParallel();
      PDB(kGlobal,1) Info("SetParallelSilent", "got %d node%s", n, n == 1 ? "" : "s");
      return n;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Tell PROOF how many slaves to use in parallel. Returns the number of
/// parallel slaves. Returns -1 in case of error.

Int_t TProof::SetParallel(Int_t nodes, Bool_t random)
{
   // If delayed startup reset settings, if required
   if (fDynamicStartup && nodes < 0) {
      if (gSystem->Getenv("PROOF_NWORKERS")) gSystem->Unsetenv("PROOF_NWORKERS");
   }

   Int_t n = SetParallelSilent(nodes, random);
   if (TestBit(TProof::kIsClient)) {
      if (n < 1) {
         Printf("PROOF set to sequential mode");
      } else {
         TString subfix = (n == 1) ? "" : "s";
         if (random)
            subfix += ", randomly selected";
         Printf("PROOF set to parallel mode (%d worker%s)", n, subfix.Data());
      }
   } else if (fDynamicStartup && nodes >= 0) {
      if (gSystem->Getenv("PROOF_NWORKERS")) gSystem->Unsetenv("PROOF_NWORKERS");
      gSystem->Setenv("PROOF_NWORKERS", TString::Format("%d", nodes));
   }
   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Add nWorkersToAdd workers to current list of workers. This function is
/// works on the master only, and only when an analysis is ongoing. A message
/// is sent back to the client when we go "more" parallel.
/// Returns -1 on error, number of total (not added!) workers on success.

Int_t TProof::GoMoreParallel(Int_t nWorkersToAdd)
{
   if (!IsValid() || !IsMaster() || IsIdle()) {
      Error("GoMoreParallel", "can't invoke here -- should not happen!");
      return -1;
   }
   if (!gProofServ && !IsLite()) {
      Error("GoMoreParallel", "no ProofServ available nor Lite -- should not happen!");
      return -1;
   }

   TSlave *sl = 0x0;
   TIter next( fSlaves );
   Int_t nAddedWorkers = 0;

   while (((nAddedWorkers < nWorkersToAdd) || (nWorkersToAdd == -1)) &&
          (( sl = dynamic_cast<TSlave *>( next() ) ))) {

      // If worker is of an invalid type, break everything: it should not happen!
      if ((sl->GetSlaveType() != TSlave::kSlave) &&
          (sl->GetSlaveType() != TSlave::kMaster)) {
         Error("GoMoreParallel", "TSlave is neither a Master nor a Slave: %s:%s",
            sl->GetName(), sl->GetOrdinal());
         R__ASSERT(0);
      }

      // Skip current worker if it is not a good candidate
      if ((!sl->IsValid()) || (fBadSlaves->FindObject(sl)) ||
          (strcmp("IGNORE", sl->GetImage()) == 0)) {
         PDB(kGlobal, 2)
            Info("GoMoreParallel", "Worker %s:%s won't be considered",
               sl->GetName(), sl->GetOrdinal());
         continue;
      }

      // Worker is good but it is already active: skip it
      if (fActiveSlaves->FindObject(sl)) {
         Info("GoMoreParallel", "Worker %s:%s is already active: skipping",
            sl->GetName(), sl->GetOrdinal());
         continue;
      }

      //
      // From here on: worker is a good candidate
      //

      if (sl->GetSlaveType() == TSlave::kSlave) {
         sl->SetStatus(TSlave::kActive);
         fActiveSlaves->Add(sl);
         fInactiveSlaves->Remove(sl);
         fActiveMonitor->Add(sl->GetSocket());
         nAddedWorkers++;
         PDB(kGlobal, 2)
            Info("GoMoreParallel", "Worker %s:%s marked as active!",
               sl->GetName(), sl->GetOrdinal());
      }
      else {
         // Can't add masters dynamically: this should not happen!
         Error("GoMoreParallel", "Dynamic addition of master is not supported");
         R__ASSERT(0);
      }

   } // end loop over all slaves

   // Get slave status (will set the slaves fWorkDir correctly)
   PDB(kGlobal, 3)
      Info("GoMoreParallel", "Will invoke AskStatistics() -- implies a Collect()");
   AskStatistics();

   // Find active slaves with unique image
   PDB(kGlobal, 3)
      Info("GoMoreParallel", "Will invoke FindUniqueSlaves()");
   FindUniqueSlaves();

   // Send new group-view to slaves
   PDB(kGlobal, 3)
      Info("GoMoreParallel", "Will invoke SendGroupView()");
   SendGroupView();

   PDB(kGlobal, 3)
      Info("GoMoreParallel", "Will invoke GetParallel()");
   Int_t nTotalWorkers = GetParallel();

   // Notify the client that we've got more workers, and print info on
   // Master's log as well
   TString s;
   s.Form("PROOF just went more parallel (%d additional worker%s, %d worker%s total)",
      nAddedWorkers, (nAddedWorkers == 1) ? "" : "s",
      nTotalWorkers, (nTotalWorkers == 1) ? "" : "s");
   if (gProofServ) gProofServ->SendAsynMessage(s);
   Info("GoMoreParallel", "%s", s.Data());

   return nTotalWorkers;
}

////////////////////////////////////////////////////////////////////////////////
/// Go in parallel mode with at most "nodes" slaves. Since the fSlaves
/// list is sorted by slave performace the active list will contain first
/// the most performant nodes. Returns the number of active slaves.
/// If random is TRUE, and nodes is less than the number of available workers,
/// a random selection is done.
/// Returns -1 in case of error.

Int_t TProof::GoParallel(Int_t nodes, Bool_t attach, Bool_t random)
{
   if (!IsValid()) return -1;

   fActiveSlaves->Clear();
   fActiveMonitor->RemoveAll();

   // Prepare the list of candidates first.
   // Algorithm depends on random option.
   TSlave *sl = 0;
   TList *wlst = new TList;
   TIter nxt(fSlaves);
   fInactiveSlaves->Clear();
   while ((sl = (TSlave *)nxt())) {
      if (sl->IsValid() && !fBadSlaves->FindObject(sl)) {
         if (strcmp("IGNORE", sl->GetImage()) == 0) continue;
         if ((sl->GetSlaveType() != TSlave::kSlave) &&
             (sl->GetSlaveType() != TSlave::kMaster)) {
            Error("GoParallel", "TSlave is neither Master nor Slave");
            R__ASSERT(0);
         }
         // Good candidate
         wlst->Add(sl);
         // Set it inactive
         fInactiveSlaves->Add(sl);
         sl->SetStatus(TSlave::kInactive);
      }
   }
   Int_t nwrks = (nodes < 0 || nodes > wlst->GetSize()) ? wlst->GetSize() : nodes;
   int cnt = 0;
   fEndMaster = TestBit(TProof::kIsMaster) ? kTRUE : kFALSE;
   while (cnt < nwrks) {
      // Random choice, if requested
      if (random) {
         Int_t iwrk = (Int_t) (gRandom->Rndm() * wlst->GetSize());
         sl = (TSlave *) wlst->At(iwrk);
      } else {
         // The first available
         sl = (TSlave *) wlst->First();
      }
      if (!sl) {
         Error("GoParallel", "attaching to candidate!");
         break;
      }
      // Remove from the list
      wlst->Remove(sl);

      Int_t slavenodes = 0;
      if (sl->GetSlaveType() == TSlave::kSlave) {
         sl->SetStatus(TSlave::kActive);
         fActiveSlaves->Add(sl);
         fInactiveSlaves->Remove(sl);
         fActiveMonitor->Add(sl->GetSocket());
         slavenodes = 1;
      } else if (sl->GetSlaveType() == TSlave::kMaster) {
         fEndMaster = kFALSE;
         TMessage mess(kPROOF_PARALLEL);
         if (!attach) {
            Int_t nn = (nodes < 0) ? -1 : nodes-cnt;
            mess << nn;
         } else {
            // To get the number of slaves
            mess.SetWhat(kPROOF_LOGFILE);
            mess << -1 << -1;
         }
         if (sl->GetSocket()->Send(mess) == -1) {
            MarkBad(sl, "could not send kPROOF_PARALLEL or kPROOF_LOGFILE request");
            slavenodes = 0;
         } else {
            Collect(sl, fCollectTimeout);
            if (sl->IsValid()) {
               sl->SetStatus(TSlave::kActive);
               fActiveSlaves->Add(sl);
               fInactiveSlaves->Remove(sl);
               fActiveMonitor->Add(sl->GetSocket());
               if (sl->GetParallel() > 0) {
                  slavenodes = sl->GetParallel();
               } else {
                  // Sequential mode: the master acts as a worker
                  slavenodes = 1;
               }
            } else {
               MarkBad(sl, "collect failed after kPROOF_PARALLEL or kPROOF_LOGFILE request");
               slavenodes = 0;
            }
         }
      }
      // 'slavenodes' may be different than 1 in multimaster setups
      cnt += slavenodes;
   }

   // Cleanup list
   wlst->SetOwner(0);
   SafeDelete(wlst);

   // Get slave status (will set the slaves fWorkDir correctly)
   AskStatistics();

   // Find active slaves with unique image
   FindUniqueSlaves();

   // Send new group-view to slaves
   if (!attach)
      SendGroupView();

   Int_t n = GetParallel();

   if (TestBit(TProof::kIsClient)) {
      if (n < 1)
         printf("PROOF set to sequential mode\n");
      else
         printf("PROOF set to parallel mode (%d worker%s)\n",
                n, n == 1 ? "" : "s");
   }

   PDB(kGlobal,1) Info("GoParallel", "got %d node%s", n, n == 1 ? "" : "s");
   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// List contents of the data directory in the sandbox.
/// This is the place where files produced by the client queries are kept

void TProof::ShowData()
{
   if (!IsValid() || !fManager) return;

   // This is run via the manager
   fManager->Find("~/data", "-type f", "all");
}

////////////////////////////////////////////////////////////////////////////////
/// Remove files for the data directory.
/// The option 'what' can take the values:
///     kPurge                 remove all files and directories under '~/data'
///     kUnregistered          remove only files not in registered datasets (default)
///     kDataset               remove files belonging to dataset 'dsname'
/// User is prompt for confirmation, unless kForceClear is ORed with the option

void TProof::ClearData(UInt_t what, const char *dsname)
{
   if (!IsValid() || !fManager) return;

   // Check whether we need to prompt
   TString prompt, a("Y");
   Bool_t force = (what & kForceClear) ? kTRUE : kFALSE;
   Bool_t doask = (!force && IsTty()) ? kTRUE : kFALSE;

   // If all just send the request
   if ((what & TProof::kPurge)) {
      // Prompt, if requested
      if (doask && !Prompt("Do you really want to remove all data files")) return;
      if (fManager->Rm("~/data/*", "-rf", "all") < 0)
         Warning("ClearData", "problems purging data directory");
      return;
   } else if ((what & TProof::kDataset)) {
      // We must have got a name
      if (!dsname || strlen(dsname) <= 0) {
         Error("ClearData", "dataset name mandatory when removing a full dataset");
         return;
      }
      // Check if the dataset is registered
      if (!ExistsDataSet(dsname)) {
         Error("ClearData", "dataset '%s' does not exists", dsname);
         return;
      }
      // Get the file content
      TFileCollection *fc = GetDataSet(dsname);
      if (!fc) {
         Error("ClearData", "could not retrieve info about dataset '%s'", dsname);
         return;
      }
      // Prompt, if requested
      TString pmpt = TString::Format("Do you really want to remove all data files"
                                     " of dataset '%s'", dsname);
      if (doask && !Prompt(pmpt.Data())) return;

      // Loop through the files
      Bool_t rmds = kTRUE;
      TIter nxf(fc->GetList());
      TFileInfo *fi = 0;
      Int_t rfiles = 0, nfiles = fc->GetList()->GetSize();
      while ((fi = (TFileInfo *) nxf())) {
         // Fill the host info
         TString host, file;
         // Take info from the current url
         if (!(fi->GetFirstUrl())) {
            Error("ClearData", "GetFirstUrl() returns NULL for '%s' - skipping",
                               fi->GetName());
            continue;
         }
         TUrl uf(*(fi->GetFirstUrl()));
         file = uf.GetFile();
         host = uf.GetHost();
         // Now search for any "file:" url
         Int_t nurl = fi->GetNUrls();
         fi->ResetUrl();
         TUrl *up = 0;
         while (nurl-- && fi->NextUrl()) {
            up = fi->GetCurrentUrl();
            if (!strcmp(up->GetProtocol(), "file")) {
               TString opt(up->GetOptions());
               if (opt.BeginsWith("node=")) {
                  host=opt;
                  host.ReplaceAll("node=","");
                  file = up->GetFile();
                  break;
               }
            }
         }
         // Issue a remove request now
         if (fManager->Rm(file.Data(), "-f", host.Data()) != 0) {
            Error("ClearData", "problems removing '%s'", file.Data());
            // Some files not removed: keep the meta info about this dataset
            rmds = kFALSE;
         }
         rfiles++;
         ClearDataProgress(rfiles, nfiles);
      }
      fprintf(stderr, "\n");
      if (rmds) {
         // All files were removed successfully: remove also the dataset meta info
         RemoveDataSet(dsname);
      }
   } else if (what & TProof::kUnregistered) {

      // Get the existing files
      TString outtmp("ProofClearData_");
      FILE *ftmp = gSystem->TempFileName(outtmp);
      if (!ftmp) {
         Error("ClearData", "cannot create temp file for logs");
         return;
      }
      fclose(ftmp);
      RedirectHandle_t h;
      gSystem->RedirectOutput(outtmp.Data(), "w", &h);
      ShowData();
      gSystem->RedirectOutput(0, 0, &h);
      // Parse the output file now
      std::ifstream in;
      in.open(outtmp.Data());
      if (!in.is_open()) {
         Error("ClearData", "could not open temp file for logs: %s", outtmp.Data());
         gSystem->Unlink(outtmp);
         return;
      }
      // Go through
      Int_t nfiles = 0;
      TMap *afmap = new TMap;
      TString line, host, file;
      Int_t from = 0;
      while (in.good()) {
         line.ReadLine(in);
         if (line.IsNull()) continue;
         while (line.EndsWith("\n")) { line.Strip(TString::kTrailing, '\n'); }
         from = 0;
         host = "";
         if (!line.Tokenize(host, from, "| ")) continue;
         file = "";
         if (!line.Tokenize(file, from, "| ")) continue;
         if (!host.IsNull() && !file.IsNull()) {
            TList *fl = (TList *) afmap->GetValue(host.Data());
            if (!fl) {
               fl = new TList();
               fl->SetName(host);
               afmap->Add(new TObjString(host), fl);
            }
            fl->Add(new TObjString(file));
            nfiles++;
            PDB(kDataset,2)
               Info("ClearData", "added info for: h:%s, f:%s", host.Data(), file.Data());
         } else {
            Warning("ClearData", "found incomplete line: '%s'", line.Data());
         }
      }
      // Close and remove the file
      in.close();
      gSystem->Unlink(outtmp);

      // Get registered data files
      TString sel = TString::Format("/%s/%s/", GetGroup(), GetUser());
      TMap *fcmap = GetDataSets(sel);
      if (!fcmap || (fcmap && fcmap->GetSize() <= 0)) {
         PDB(kDataset,1)
         Warning("ClearData", "no dataset beloning to '%s'", sel.Data());
         SafeDelete(fcmap);
      }

      // Go thorugh and prepare the lists per node
      TString opt;
      TObjString *os = 0;
      if (fcmap) {
         TIter nxfc(fcmap);
         while ((os = (TObjString *) nxfc())) {
            TFileCollection *fc = 0;
            if ((fc = (TFileCollection *) fcmap->GetValue(os))) {
               TFileInfo *fi = 0;
               TIter nxfi(fc->GetList());
               while ((fi = (TFileInfo *) nxfi())) {
                  // Get special "file:" url
                  fi->ResetUrl();
                  Int_t nurl = fi->GetNUrls();
                  TUrl *up = 0;
                  while (nurl-- && fi->NextUrl()) {
                     up = fi->GetCurrentUrl();
                     if (!strcmp(up->GetProtocol(), "file")) {
                        opt = up->GetOptions();
                        if (opt.BeginsWith("node=")) {
                           host=opt;
                           host.ReplaceAll("node=","");
                           file = up->GetFile();
                           PDB(kDataset,2)
                              Info("ClearData", "found: host: %s, file: %s", host.Data(), file.Data());
                           // Remove this from the full list, if there
                           TList *fl = (TList *) afmap->GetValue(host.Data());
                           if (fl) {
                              TObjString *fn = (TObjString *) fl->FindObject(file.Data());
                              if (fn) {
                                 fl->Remove(fn);
                                 SafeDelete(fn);
                                 nfiles--;
                              } else {
                                 Warning("ClearData",
                                         "registered file '%s' not found in the full list!",
                                         file.Data());
                              }
                           }
                           break;
                        }
                     }
                  }
               }
            }
         }
         // Clean up the the received map
         if (fcmap) fcmap->SetOwner(kTRUE);
         SafeDelete(fcmap);
      }
      // List of the files to be removed
      Info("ClearData", "%d unregistered files to be removed:", nfiles);
      afmap->Print();
      // Prompt, if requested
      TString pmpt = TString::Format("Do you really want to remove all %d"
                                     " unregistered data files", nfiles);
      if (doask && !Prompt(pmpt.Data())) return;

      // Remove one by one; we may implement a bloc remove in the future
      Int_t rfiles = 0;
      TIter nxls(afmap);
      while ((os = (TObjString *) nxls())) {
         TList *fl = 0;
         if ((fl = (TList *) afmap->GetValue(os))) {
            TIter nxf(fl);
            TObjString *fn = 0;
            while ((fn = (TObjString *) nxf())) {
               // Issue a remove request now
               if (fManager->Rm(fn->GetName(), "-f", os->GetName()) != 0) {
                  Error("ClearData", "problems removing '%s' on host '%s'",
                                     fn->GetName(), os->GetName());
               }
               rfiles++;
               ClearDataProgress(rfiles, nfiles);
            }
         }
      }
      fprintf(stderr, "\n");
      // Final cleanup
      afmap->SetOwner(kTRUE);
      SafeDelete(afmap);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Prompt the question 'p' requiring an answer y,Y,n,N
/// Return kTRUE is the answer was y or Y, kFALSE in all other cases.

Bool_t TProof::Prompt(const char *p)
{
   TString pp(p);
   if (!pp.Contains("?")) pp += "?";
   if (!pp.Contains("[y/N]")) pp += " [y/N]";
   TString a = Getline(pp.Data());
   if (a != "\n" && a[0] != 'y' &&  a[0] != 'Y' &&  a[0] != 'n' &&  a[0] != 'N') {
      Printf("Please answer y, Y, n or N");
      // Unclear answer: assume negative
      return kFALSE;
   } else if (a == "\n" || a[0] == 'n' ||  a[0] == 'N') {
      // Explicitly Negative answer
      return kFALSE;
   }
   // Explicitly Positive answer
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Progress bar for clear data

void TProof::ClearDataProgress(Int_t r, Int_t t)
{
   fprintf(stderr, "[TProof::ClearData] Total %5d files\t|", t);
   for (Int_t l = 0; l < 20; l++) {
      if (r > 0 && t > 0) {
         if (l < 20*r/t)
            fprintf(stderr, "=");
         else if (l == 20*r/t)
            fprintf(stderr, ">");
         else if (l > 20*r/t)
            fprintf(stderr, ".");
      } else
         fprintf(stderr, "=");
   }
   fprintf(stderr, "| %.02f %%      \r", 100.0*(t ? (r/t) : 1));
}

////////////////////////////////////////////////////////////////////////////////
/// List contents of file cache. If all is true show all caches also on
/// slaves. If everything is ok all caches are to be the same.

void TProof::ShowCache(Bool_t all)
{
   if (!IsValid()) return;

   TMessage mess(kPROOF_CACHE);
   mess << Int_t(kShowCache) << all;
   Broadcast(mess, kUnique);

   if (all) {
      TMessage mess2(kPROOF_CACHE);
      mess2 << Int_t(kShowSubCache) << all;
      Broadcast(mess2, fNonUniqueMasters);

      Collect(kAllUnique, fCollectTimeout);
   } else {
      Collect(kUnique, fCollectTimeout);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove file from all file caches. If file is 0 or "" or "*", remove all
/// the files

void TProof::ClearCache(const char *file)
{
   if (!IsValid()) return;

   TMessage mess(kPROOF_CACHE);
   mess << Int_t(kClearCache) << TString(file);
   Broadcast(mess, kUnique);

   TMessage mess2(kPROOF_CACHE);
   mess2 << Int_t(kClearSubCache) << TString(file);
   Broadcast(mess2, fNonUniqueMasters);

   Collect(kAllUnique);

   // clear file map so files get send again to remote nodes
   fFileMap.clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Exec system command 'cmd'. If fdout > -1, append the output to fdout.

void TProof::SystemCmd(const char *cmd, Int_t fdout)
{
   if (fdout < 0) {
      // Exec directly the command
      gSystem->Exec(cmd);
   } else {
      // Exec via a pipe
      FILE *fin = gSystem->OpenPipe(cmd, "r");
      if (fin) {
         // Now we go
         char line[2048];
         while (fgets(line, 2048, fin)) {
            Int_t r = strlen(line);
            if (r > 0) {
               if (write(fdout, line, r) < 0) {
                  ::Warning("TProof::SystemCmd",
                            "errno %d writing to file descriptor %d",
                            TSystem::GetErrno(), fdout);
               }
            } else {
               // Done
               break;
            }
         }
         gSystem->ClosePipe(fin);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// List contents of package directory. If all is true show all package
/// directories also on slaves. If everything is ok all package directories
/// should be the same. If redir is kTRUE the result is redirected to the log
/// file (option available for internal actions).

void TProof::ShowPackages(Bool_t all, Bool_t redirlog)
{
   if (!IsValid()) return;

   Bool_t oldredir = fRedirLog;
   if (redirlog) fRedirLog = kTRUE;

   // Active logging unit
   FILE *fout = (fRedirLog) ? fLogFileW : stdout;
   if (!fout) {
      Warning("ShowPackages", "file descriptor for outputs undefined (%p):"
              " will not log msgs", fout);
      return;
   }
   lseek(fileno(fout), (off_t) 0, SEEK_END);

   if (TestBit(TProof::kIsClient)) {
      fPackMgr->Show();
   }

   // Nothing more to do if we are a Lite-session
   if (IsLite()) {
      fRedirLog = oldredir;
      return;
   }

   TMessage mess(kPROOF_CACHE);
   mess << Int_t(kShowPackages) << all;
   Broadcast(mess, kUnique);

   if (all) {
      TMessage mess2(kPROOF_CACHE);
      mess2 << Int_t(kShowSubPackages) << all;
      Broadcast(mess2, fNonUniqueMasters);

      Collect(kAllUnique, fCollectTimeout);
   } else {
      Collect(kUnique, fCollectTimeout);
   }
   // Restore logging option
   fRedirLog = oldredir;
}

////////////////////////////////////////////////////////////////////////////////
/// List which packages are enabled. If all is true show enabled packages
/// for all active slaves. If everything is ok all active slaves should
/// have the same packages enabled.

void TProof::ShowEnabledPackages(Bool_t all)
{
   if (!IsValid()) return;

   if (TestBit(TProof::kIsClient)) {
      fPackMgr->ShowEnabled(TString::Format("*** Enabled packages on client on %s\n",
                            gSystem->HostName()));
   }

   // Nothing more to do if we are a Lite-session
   if (IsLite()) return;

   TMessage mess(kPROOF_CACHE);
   mess << Int_t(kShowEnabledPackages) << all;
   Broadcast(mess);
   Collect(kActive, fCollectTimeout);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all packages.
/// Returns 0 in case of success and -1 in case of error.

Int_t TProof::ClearPackages()
{
   if (!IsValid()) return -1;

   if (UnloadPackages() == -1)
      return -1;

   if (DisablePackages() == -1)
      return -1;

   return fStatus;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a specific package.
/// Returns 0 in case of success and -1 in case of error.

Int_t TProof::ClearPackage(const char *package)
{
   if (!IsValid()) return -1;

   if (!package || !package[0]) {
      Error("ClearPackage", "need to specify a package name");
      return -1;
   }

   // if name, erroneously, is a par pathname strip off .par and path
   TString pac = package;
   if (pac.EndsWith(".par"))
      pac.Remove(pac.Length()-4);
   pac = gSystem->BaseName(pac);

   if (UnloadPackage(pac) == -1)
      return -1;

   if (DisablePackage(pac) == -1)
      return -1;

   return fStatus;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a specific package.
/// Returns 0 in case of success and -1 in case of error.

Int_t TProof::DisablePackage(const char *pack)
{
   if (!IsValid()) return -1;

   if (!pack || strlen(pack) <= 0) {
      Error("DisablePackage", "need to specify a package name");
      return -1;
   }

   // if name, erroneously, is a par pathname strip off .par and path
   TString pac = pack;
   if (pac.EndsWith(".par"))
      pac.Remove(pac.Length()-4);
   pac = gSystem->BaseName(pac);

   if (fPackMgr->Remove(pack) < 0)
      Warning("DisablePackage", "problem removing locally package '%s'", pack);

   // Nothing more to do if we are a Lite-session
   if (IsLite()) return 0;

   Int_t st = -1;
   Bool_t done = kFALSE;
   if (fManager) {
      // Try to do it via XROOTD (new way)
      TString path;
      path.Form("~/packages/%s", pack);
      if (fManager->Rm(path, "-rf", "all") != -1) {
         path.Append(".par");
         if (fManager->Rm(path, "-f", "all") != -1) {
            done = kTRUE;
            st = 0;
         }
      }
   }
   if (!done) {
      // Try via TProofServ (old way)
      TMessage mess(kPROOF_CACHE);
      mess << Int_t(kDisablePackage) << pac;
      Broadcast(mess, kUnique);

      TMessage mess2(kPROOF_CACHE);
      mess2 << Int_t(kDisableSubPackage) << pac;
      Broadcast(mess2, fNonUniqueMasters);

      Collect(kAllUnique);
      st = fStatus;
   }

   // Done
   return st;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all packages.
/// Returns 0 in case of success and -1 in case of error.

Int_t TProof::DisablePackages()
{
   if (!IsValid()) return -1;

   // remove all packages on client
   if (fPackMgr->Remove(nullptr) < 0)
      Warning("DisablePackages", "problem removing packages locally");

   // Nothing more to do if we are a Lite-session
   if (IsLite()) return 0;

   Int_t st = -1;
   Bool_t done = kFALSE;
   if (fManager) {
      // Try to do it via XROOTD (new way)
      if (fManager->Rm("~/packages/*", "-rf", "all") != -1) {
         done = kTRUE;
         st = 0;
      }
   }
   if (!done) {

      TMessage mess(kPROOF_CACHE);
      mess << Int_t(kDisablePackages);
      Broadcast(mess, kUnique);

      TMessage mess2(kPROOF_CACHE);
      mess2 << Int_t(kDisableSubPackages);
      Broadcast(mess2, fNonUniqueMasters);

      Collect(kAllUnique);
      st = fStatus;
   }

   // Done
   return st;
}

////////////////////////////////////////////////////////////////////////////////
/// Build specified package. Executes the PROOF-INF/BUILD.sh
/// script if it exists on all unique nodes. If opt is kBuildOnSlavesNoWait
/// then submit build command to slaves, but don't wait
/// for results. If opt is kCollectBuildResults then collect result
/// from slaves. To be used on the master.
/// If opt = kBuildAll (default) then submit and wait for results
/// (to be used on the client).
/// Returns 0 in case of success and -1 in case of error.

Int_t TProof::BuildPackage(const char *package,
                           EBuildPackageOpt opt, Int_t chkveropt, TList *workers)
{
   if (!IsValid()) return -1;

   if (!package || !package[0]) {
      Error("BuildPackage", "need to specify a package name");
      return -1;
   }

   // if name, erroneously, is a par pathname strip off .par and path
   TString pac = package;
   if (pac.EndsWith(".par"))
      pac.Remove(pac.Length()-4);
   pac = gSystem->BaseName(pac);

   Bool_t buildOnClient = kTRUE;
   if (opt == kDontBuildOnClient) {
      buildOnClient = kFALSE;
      opt = kBuildAll;
   }
   // Prepare the local package
   TString pdir;
   Int_t st = 0;

   if (opt <= kBuildAll && (!IsLite() || !buildOnClient)) {
      if (workers) {
         TMessage mess(kPROOF_CACHE);
         mess << Int_t(kBuildPackage) << pac << chkveropt;
         Broadcast(mess, workers);

      } else {
         TMessage mess(kPROOF_CACHE);
         mess << Int_t(kBuildPackage) << pac << chkveropt;
         Broadcast(mess, kUnique);

         TMessage mess2(kPROOF_CACHE);
         mess2 << Int_t(kBuildSubPackage) << pac << chkveropt;
         Broadcast(mess2, fNonUniqueMasters);
      }
   }

   if (opt >= kBuildAll) {
      // by first forwarding the build commands to the master and slaves
      // and only then building locally we build in parallel
      if (buildOnClient) {
         st = fPackMgr->Build(pac, chkveropt);
      }


      fStatus = 0;
      if (!IsLite() || !buildOnClient) {

        // On the master, workers that fail are deactivated
        // Bool_t deactivateOnFailure = (IsMaster()) ? kTRUE : kFALSE;
         if (workers) {
//            Collect(workers, -1, -1, deactivateOnFailure);
            Collect(workers);
         } else {
            Collect(kAllUnique);
         }
      }

      if (fStatus < 0 || st < 0)
         return -1;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Load specified package. Executes the PROOF-INF/SETUP.C script
/// on all active nodes. If notOnClient = true, don't load package
/// on the client. The default is to load the package also on the client.
/// The argument 'loadopts' specify a list of objects to be passed to the SETUP.
/// The objects in the list must be streamable; the SETUP macro will be executed
/// like this: SETUP.C(loadopts).
/// Returns 0 in case of success and -1 in case of error.

Int_t TProof::LoadPackage(const char *package, Bool_t notOnClient,
   TList *loadopts, TList *workers)
{
   if (!IsValid()) return -1;

   if (!package || !package[0]) {
      Error("LoadPackage", "need to specify a package name");
      return -1;
   }

   // if name, erroneously, is a par pathname strip off .par and path
   TString pac = package;
   if (pac.EndsWith(".par"))
      pac.Remove(pac.Length()-4);
   pac = gSystem->BaseName(pac);

   if (!notOnClient && TestBit(TProof::kIsClient))
      if (fPackMgr->Load(package, loadopts) == -1) return -1;

   TMessage mess(kPROOF_CACHE);
   mess << Int_t(kLoadPackage) << pac;
   if (loadopts) mess << loadopts;

   // On the master, workers that fail are deactivated
   Bool_t deactivateOnFailure = (IsMaster()) ? kTRUE : kFALSE;

   Bool_t doCollect = (fDynamicStartup && !IsIdle()) ? kFALSE : kTRUE;

   if (workers) {
      PDB(kPackage, 3)
         Info("LoadPackage", "Sending load message to selected workers only");
      Broadcast(mess, workers);
      if (doCollect) Collect(workers, -1, -1, deactivateOnFailure);
   } else {
      Broadcast(mess);
      Collect(kActive, -1, -1, deactivateOnFailure);
   }

   return fStatus;
}

////////////////////////////////////////////////////////////////////////////////
/// Unload specified package.
/// Returns 0 in case of success and -1 in case of error.

Int_t TProof::UnloadPackage(const char *package)
{
   if (!IsValid()) return -1;

   if (!package || !package[0]) {
      Error("UnloadPackage", "need to specify a package name");
      return -1;
   }

   // if name, erroneously, is a par pathname strip off .par and path
   TString pac = package;
   if (pac.EndsWith(".par"))
      pac.Remove(pac.Length()-4);
   pac = gSystem->BaseName(pac);

   if (fPackMgr->Unload(package) < 0)
      Warning("UnloadPackage", "unable to remove symlink to %s", package);

   // Nothing more to do if we are a Lite-session
   if (IsLite()) return 0;

   TMessage mess(kPROOF_CACHE);
   mess << Int_t(kUnloadPackage) << pac;
   Broadcast(mess);
   Collect();

   return fStatus;
}

////////////////////////////////////////////////////////////////////////////////
/// Unload all packages.
/// Returns 0 in case of success and -1 in case of error.

Int_t TProof::UnloadPackages()
{
   if (!IsValid()) return -1;

   if (TestBit(TProof::kIsClient)) {
      if (fPackMgr->Unload(0) < 0) return -1;
   }

   // Nothing more to do if we are a Lite-session
   if (IsLite()) return 0;

   TMessage mess(kPROOF_CACHE);
   mess << Int_t(kUnloadPackages);
   Broadcast(mess);
   Collect();

   return fStatus;
}

////////////////////////////////////////////////////////////////////////////////
/// Enable specified package. Executes the PROOF-INF/BUILD.sh
/// script if it exists followed by the PROOF-INF/SETUP.C script.
/// In case notOnClient = true, don't enable the package on the client.
/// The default is to enable packages also on the client.
/// If specified, enables packages only on the specified workers.
/// Returns 0 in case of success and -1 in case of error.
/// Provided for backward compatibility.

Int_t TProof::EnablePackage(const char *package, Bool_t notOnClient,
   TList *workers)
{
   return EnablePackage(package, (TList *)0, notOnClient, workers);
}

////////////////////////////////////////////////////////////////////////////////
/// Enable specified package. Executes the PROOF-INF/BUILD.sh
/// script if it exists followed by the PROOF-INF/SETUP.C script.
/// In case notOnClient = true, don't enable the package on the client.
/// The default is to enable packages also on the client.
/// It is is possible to specify options for the loading step via 'loadopts';
/// the string will be passed passed as argument to SETUP.
/// Special option 'chkv=<o>' (or 'checkversion=<o>') can be used to control
/// plugin version checking during building: possible choices are:
///     off         no check; failure may occur at loading
///     on          check ROOT version [default]
///     svn         check ROOT version and Git commit SHA1.
/// (Use ';', ' ' or '|' to separate 'chkv=<o>' from the rest.)
/// If specified, enables packages only on the specified workers.
/// Returns 0 in case of success and -1 in case of error.

Int_t TProof::EnablePackage(const char *package, const char *loadopts,
                            Bool_t notOnClient, TList *workers)
{
   TList *optls = 0;
   if (loadopts && strlen(loadopts)) {
      if (fProtocol > 28) {
         TObjString *os = new TObjString(loadopts);
         // Filter out 'checkversion=off|on|svn' or 'chkv=...'
         os->String().ReplaceAll("checkversion=", "chkv=");
         Ssiz_t fcv = kNPOS, lcv = kNPOS;
         if ((fcv = os->String().Index("chkv=")) !=  kNPOS) {
            TRegexp re("[; |]");
            if ((lcv = os->String().Index(re, fcv)) == kNPOS) {
               lcv = os->String().Length();
            }
            TString ocv = os->String()(fcv, lcv - fcv);
            Int_t cvopt = -1;
            if (ocv.EndsWith("=off") || ocv.EndsWith("=0"))
               cvopt = (Int_t) TPackMgr::kDontCheck;
            else if (ocv.EndsWith("=on") || ocv.EndsWith("=1"))
               cvopt = (Int_t) TPackMgr::kCheckROOT;
            else
               Warning("EnablePackage", "'checkversion' option unknown from argument: '%s' - ignored", ocv.Data());
            if (cvopt > -1) {
               if (gDebug > 0)
                  Info("EnablePackage", "setting check version option from argument: %d", cvopt);
               optls = new TList;
               optls->Add(new TParameter<Int_t>("PROOF_Package_CheckVersion", (Int_t) cvopt));
               // Remove the special option from; we leave a separator if there were two (one before and one after)
               if (lcv != kNPOS && fcv == 0) ocv += os->String()[lcv];
               if (fcv > 0 && os->String().Index(re, fcv - 1) == fcv - 1) os->String().Remove(fcv - 1, 1);
               os->String().ReplaceAll(ocv.Data(), "");
            }
         }
         if (!os->String().IsNull()) {
            if (!optls) optls = new TList;
            optls->Add(new TObjString(os->String().Data()));
         }
         if (optls) optls->SetOwner(kTRUE);
      } else {
         // Notify
         Warning("EnablePackage", "remote server does not support options: ignoring the option string");
      }
   }
   // Run
   Int_t rc = EnablePackage(package, optls, notOnClient, workers);
   // Clean up
   SafeDelete(optls);
   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Enable specified package. Executes the PROOF-INF/BUILD.sh
/// script if it exists followed by the PROOF-INF/SETUP.C script.
/// In case notOnClient = true, don't enable the package on the client.
/// The default is to enable packages also on the client.
/// It is is possible to specify a list of objects to be passed to the SETUP
/// functions via 'loadopts'; the objects must be streamable.
/// Returns 0 in case of success and -1 in case of error.

Int_t TProof::EnablePackage(const char *package, TList *loadopts,
                            Bool_t notOnClient, TList *workers)
{
   if (!IsValid()) return -1;

   if (!package || !package[0]) {
      Error("EnablePackage", "need to specify a package name");
      return -1;
   }

   // if name, erroneously, is a par pathname strip off .par and path
   TString pac = package;
   if (pac.EndsWith(".par"))
      pac.Remove(pac.Length()-4);
   pac = gSystem->BaseName(pac);

   EBuildPackageOpt opt = kBuildAll;
   if (notOnClient)
      opt = kDontBuildOnClient;

   // Get check version option; user settings have priority
   Int_t chkveropt = TPackMgr::kCheckROOT;
   TString ocv = gEnv->GetValue("Proof.Package.CheckVersion", "");
   if (!ocv.IsNull()) {
      if (ocv == "off" || ocv == "0")
         chkveropt = (Int_t) TPackMgr::kDontCheck;
      else if (ocv == "on" || ocv == "1")
         chkveropt = (Int_t) TPackMgr::kCheckROOT;
      else
         Warning("EnablePackage", "'checkversion' option unknown from rootrc: '%s' - ignored", ocv.Data());
   }
   if (loadopts) {
      TParameter<Int_t> *pcv = (TParameter<Int_t> *) loadopts->FindObject("PROOF_Package_CheckVersion");
      if (pcv) {
         chkveropt = pcv->GetVal();
         loadopts->Remove(pcv);
         delete pcv;
      }
   }
  if (gDebug > 0)
      Info("EnablePackage", "using check version option: %d", chkveropt);

   if (BuildPackage(pac, opt, chkveropt, workers) == -1)
      return -1;

   TList *optls = (loadopts && loadopts->GetSize() > 0) ? loadopts : 0;
   if (optls && fProtocol <= 28) {
      Warning("EnablePackage", "remote server does not support options: ignoring the option list");
      optls = 0;
   }

   if (LoadPackage(pac, notOnClient, optls, workers) == -1)
      return -1;

   // Record the information for later usage (simulation of dynamic start on PROOF-Lite)
   if (!fEnabledPackagesOnCluster) {
      fEnabledPackagesOnCluster = new TList;
      fEnabledPackagesOnCluster->SetOwner();
   }
   if (!fEnabledPackagesOnCluster->FindObject(pac)) {
      TPair *pck = (optls && optls->GetSize() > 0) ? new TPair(new TObjString(pac), optls->Clone())
                                                   : new TPair(new TObjString(pac), 0);
      fEnabledPackagesOnCluster->Add(pck);
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Download a PROOF archive (PAR file) from the master package repository.
/// The PAR file is downloaded in the current directory or in the directory
/// specified by 'dstdir'. If a package with the same name already exists
/// at destination, a check on the MD5 sum is done and the user warned or
/// prompted for action, depending is the file is equal or different.
/// Returns 0 in case of success and -1 in case of error.

Int_t TProof::DownloadPackage(const char *pack, const char *dstdir)
{
   if (!fManager || !(fManager->IsValid())) {
      Error("DownloadPackage", "the manager is undefined!");
      return -1;
   }

   // Create the default source and destination paths
   TString parname(gSystem->BaseName(pack)), src, dst;
   if (!parname.EndsWith(".par")) parname += ".par";
   src.Form("packages/%s", parname.Data());
   if (!dstdir || strlen(dstdir) <= 0) {
      dst.Form("./%s", parname.Data());
   } else {
      // Check the destination directory
      FileStat_t st;
      if (gSystem->GetPathInfo(dstdir, st) != 0) {
         // Directory does not exit: create it
         if (gSystem->mkdir(dstdir, kTRUE) != 0) {
            Error("DownloadPackage",
                  "could not create the destination directory '%s' (errno: %d)",
                   dstdir, TSystem::GetErrno());
            return -1;
         }
      } else if (!R_ISDIR(st.fMode) && !R_ISLNK(st.fMode)) {
         Error("DownloadPackage",
               "destination path '%s' exist but is not a directory!", dstdir);
         return -1;
      }
      dst.Form("%s/%s", dstdir, parname.Data());
   }

   // Make sure the source file exists
   FileStat_t stsrc;
   RedirectHandle_t rh;
   if (gSystem->RedirectOutput(fLogFileName, "a", &rh) != 0)
      Warning("DownloadPackage", "problems redirecting output to '%s'", fLogFileName.Data());
   Int_t rc = fManager->Stat(src, stsrc);
   if (gSystem->RedirectOutput(0, 0, &rh) != 0)
      Warning("DownloadPackage", "problems restoring output");
   if (rc != 0) {
      // Check if there is another possible source
      ShowPackages(kFALSE, kTRUE);
      TMacro *mp = GetLastLog();
      if (mp) {
         // Look for global directories
         Bool_t isGlobal = kFALSE;
         TIter nxl(mp->GetListOfLines());
         TObjString *os = 0;
         TString globaldir;
         while ((os = (TObjString *) nxl())) {
            TString s(os->GetName());
            if (s.Contains("*** Global Package cache")) {
               // Get the directory
               s.Remove(0, s.Last(':') + 1);
               s.Remove(s.Last(' '));
               globaldir = s;
               isGlobal = kTRUE;
            } else if (s.Contains("*** Package cache")) {
               isGlobal = kFALSE;
               globaldir = "";
            }
            // Check for the package
            if (isGlobal && s.Contains(parname)) {
               src.Form("%s/%s", globaldir.Data(), parname.Data());
               break;
            }
         }
         // Cleanup
         delete mp;
      }
   }

   // Do it via the manager
   if (fManager->GetFile(src, dst, "silent") != 0) {
      Error("DownloadPackage", "problems downloading '%s' (src:%s, dst:%s)",
                                pack, src.Data(), dst.Data());
      return -1;
   } else {
      Info("DownloadPackage", "'%s' cross-checked against master repository (local path: %s)",
                                pack, dst.Data());
   }
   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Upload a PROOF archive (PAR file). A PAR file is a compressed
/// tar file with one special additional directory, PROOF-INF
/// (blatantly copied from Java's jar format). It must have the extension
/// .par. A PAR file can be directly a binary or a source with a build
/// procedure. In the PROOF-INF directory there can be a build script:
/// BUILD.sh to be called to build the package, in case of a binary PAR
/// file don't specify a build script or make it a no-op. Then there is
/// SETUP.C which sets the right environment variables to use the package,
/// like LD_LIBRARY_PATH, etc.
/// The 'opt' allows to specify whether the .PAR should be just unpacked
/// in the existing dir (opt = kUntar, default) or a remove of the existing
/// directory should be executed (opt = kRemoveOld), so triggering a full
/// re-build. The option if effective only for PROOF protocol > 8 .
/// The lab 'dirlab' (e.g. 'G0') indicates that the package is to uploaded to
/// an alternative global directory for global usage. This may require special
/// privileges.
/// If download is kTRUE and the package is not found locally, then it is downloaded
/// from the master repository.
/// Returns 0 in case of success and -1 in case of error.

Int_t TProof::UploadPackage(const char *pack, EUploadPackageOpt opt,
   TList *workers)
{
   if (!IsValid()) return -1;

   // Remote PAR ?
   TFile::EFileType ft = TFile::GetType(pack);
   Bool_t remotepar = (ft == TFile::kWeb || ft == TFile::kNet) ? kTRUE : kFALSE;

   TString par(pack), base, name;
   if (par.EndsWith(".par")) {
      base = gSystem->BaseName(par);
      name = base(0, base.Length() - strlen(".par"));
   } else {
      name = gSystem->BaseName(par);
      base.Form("%s.par", name.Data());
      par += ".par";
   }

   // Default location is the local working dir; then the package dir
   gSystem->ExpandPathName(par);
   if (gSystem->AccessPathName(par, kReadPermission)) {
      Int_t xrc = -1;
      if (!remotepar) xrc = TPackMgr::FindParPath(fPackMgr, name, par);
      if (xrc == 0) {
         // Package is in the global dirs
         if (gDebug > 0)
            Info("UploadPackage", "global package found (%s): no upload needed",
                                  par.Data());
         return 0;
      } else if (xrc < 0) {
         Error("UploadPackage", "PAR file '%s' not found", par.Data());
         return -1;
      }
   }

   // Strategy:
   // On the client:
   // get md5 of package and check if it is different
   // from the one stored in the local package directory. If it is lock
   // the package directory and copy the package, unlock the directory.
   // On the masters:
   // get md5 of package and check if it is different from the
   // one stored on the remote node. If it is different lock the remote
   // package directory and use TFTP or SendFile to ftp the package to the
   // remote node, unlock the directory.


   if (TestBit(TProof::kIsClient)) {
      Bool_t rmold = (opt == TProof::kRemoveOld) ? kTRUE : kFALSE;
      if (fPackMgr->Install(par, rmold) < 0) {
         Error("UploadPackage", "installing '%s' failed", gSystem->BaseName(par));
         return -1;
      }
   }

   // Nothing more to do if we are a Lite-session
   if (IsLite()) return 0;

   TMD5 *md5 = fPackMgr->ReadMD5(name);

   TString smsg;
   if (remotepar && GetRemoteProtocol() > 36) {
      smsg.Form("+%s", par.Data());
   } else {
      smsg.Form("+%s", base.Data());
   }

   TMessage mess(kPROOF_CHECKFILE);
   mess << smsg << (*md5);
   TMessage mess2(kPROOF_CHECKFILE);
   smsg.Replace(0, 1, "-");
   mess2 << smsg << (*md5);
   TMessage mess3(kPROOF_CHECKFILE);
   smsg.Replace(0, 1, "=");
   mess3 << smsg << (*md5);

   delete md5;

   if (fProtocol > 8) {
      // Send also the option
      mess << (UInt_t) opt;
      mess2 << (UInt_t) opt;
      mess3 << (UInt_t) opt;
   }

   // Loop over all slaves with unique fs image, or to a selected
   // list of workers, if specified
   if (!workers)
      workers = fUniqueSlaves;
   TIter next(workers);
   TSlave *sl = 0;
   while ((sl = (TSlave *) next())) {
      if (!sl->IsValid())
         continue;

      sl->GetSocket()->Send(mess);

      fCheckFileStatus = 0;
      Collect(sl, fCollectTimeout, kPROOF_CHECKFILE);
      if (fCheckFileStatus == 0) {

         if (fProtocol > 5) {
            // remote directory is locked, upload file over the open channel
            smsg.Form("%s/%s/%s", sl->GetProofWorkDir(), kPROOF_PackDir, base.Data());
            if (SendFile(par, (kBinary | kForce | kCpBin | kForward), smsg.Data(), sl) < 0) {
               Error("UploadPackage", "%s: problems uploading file %s",
                                      sl->GetOrdinal(), par.Data());
               return -1;
            }
         } else {
            // old servers receive it via TFTP
            TFTP ftp(TString("root://")+sl->GetName(), 1);
            if (!ftp.IsZombie()) {
               smsg.Form("%s/%s", sl->GetProofWorkDir(), kPROOF_PackDir);
               ftp.cd(smsg.Data());
               ftp.put(par, base.Data());
            }
         }

         // install package and unlock dir
         sl->GetSocket()->Send(mess2);
         fCheckFileStatus = 0;
         Collect(sl, fCollectTimeout, kPROOF_CHECKFILE);
         if (fCheckFileStatus == 0) {
            Error("UploadPackage", "%s: unpacking of package %s failed",
                                   sl->GetOrdinal(), base.Data());
            return -1;
         }
      }
   }

   // loop over all other master nodes
   TIter nextmaster(fNonUniqueMasters);
   TSlave *ma;
   while ((ma = (TSlave *) nextmaster())) {
      if (!ma->IsValid())
         continue;

      ma->GetSocket()->Send(mess3);

      fCheckFileStatus = 0;
      Collect(ma, fCollectTimeout, kPROOF_CHECKFILE);
      if (fCheckFileStatus == 0) {
         // error -> package should have been found
         Error("UploadPackage", "package %s did not exist on submaster %s",
               base.Data(), ma->GetOrdinal());
         return -1;
      }
   }

   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Make sure that the directory path contained by macro is in the macro path

void TProof::AssertMacroPath(const char *macro)
{
   static TString macrop(gROOT->GetMacroPath());
   if (macro && strlen(macro) > 0) {
      TString dirn = gSystem->GetDirName(macro);
      if (!macrop.Contains(dirn)) {
         macrop += TString::Format("%s:", dirn.Data());
         gROOT->SetMacroPath(macrop);
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Load the specified macro on master, workers and, if notOnClient is
/// kFALSE, on the client. The macro file is uploaded if new or updated.
/// Additional files to be uploaded (or updated, if needed) can be specified
/// after a comma, e.g. "mymacro.C+,thisheader.h,thatheader.h".
/// If existing in the same directory, a header basename(macro).h or .hh, is also
/// uploaded.
/// The default is to load the macro also on the client; notOnClient can be used
/// to avoid loading on the client.
/// On masters, if uniqueWorkers is kTRUE, the macro is loaded on unique workers
/// only, and collection is not done; if uniqueWorkers is kFALSE, collection
/// from the previous request is done, and broadcasting + collection from the
/// other workers is done.
/// The wrks arg can be used on the master to limit the set of workers.
/// Returns 0 in case of success and -1 in case of error.

Int_t TProof::Load(const char *macro, Bool_t notOnClient, Bool_t uniqueWorkers,
                   TList *wrks)
{
   if (!IsValid()) return -1;

   if (!macro || !macro[0]) {
      Error("Load", "need to specify a macro name");
      return -1;
   }

   // Make sure the path is in the macro path
   TProof::AssertMacroPath(macro);

   if (TestBit(TProof::kIsClient) && !wrks) {

      // Extract the file implementation name first
      TString addsname, implname = macro;
      Ssiz_t icom = implname.Index(",");
      if (icom != kNPOS) {
         addsname = implname(icom + 1, implname.Length());
         implname.Remove(icom);
      }
      TString basemacro = gSystem->BaseName(implname), mainmacro(implname);
      TString bmsg(basemacro), acmode, args, io;
      implname = gSystem->SplitAclicMode(implname, acmode, args, io);

      // Macro names must have a standard format
      Int_t dot = implname.Last('.');
      if (dot == kNPOS) {
         Info("Load", "macro '%s' does not contain a '.': do nothing", macro);
         return -1;
      }

      // Is there any associated header file
      Bool_t hasHeader = kTRUE;
      TString headname = implname;
      headname.Remove(dot);
      headname += ".h";
      if (gSystem->AccessPathName(headname, kReadPermission)) {
         TString h = headname;
         headname.Remove(dot);
         headname += ".hh";
         if (gSystem->AccessPathName(headname, kReadPermission)) {
            hasHeader = kFALSE;
            if (gDebug > 0)
               Info("Load", "no associated header file found: tried: %s %s",
                            h.Data(), headname.Data());
         }
      }

      // Is there any additional file ?
      TString addincs;
      TList addfiles;
      if (!addsname.IsNull()) {
         TString fn;
         Int_t from = 0;
         while (addsname.Tokenize(fn, from, ",")) {
            if (gSystem->AccessPathName(fn, kReadPermission)) {
               Error("Load", "additional file '%s' not found", fn.Data());
               return -1;
            }
            // Create the additional include statement
            if (!notOnClient) {
               TString dirn = gSystem->GetDirName(fn);
               if (addincs.IsNull()) {
                  addincs.Form("-I%s", dirn.Data());
               } else if (!addincs.Contains(dirn)) {
                  addincs += TString::Format(" -I%s", dirn.Data());
               }
            }
            // Remember these files ...
            addfiles.Add(new TObjString(fn));
         }
      }

      // Send files now; the md5 check is run here; see SendFile for more
      // details.
      if (SendFile(implname, kAscii | kForward , "cache") == -1) {
         Error("Load", "problems sending implementation file %s", implname.Data());
         return -1;
      }
      if (hasHeader)
         if (SendFile(headname, kAscii | kForward , "cache") == -1) {
            Error("Load", "problems sending header file %s", headname.Data());
            return -1;
         }
      // Additional files
      if (addfiles.GetSize() > 0) {
         TIter nxfn(&addfiles);
         TObjString *os = 0;
         while ((os = (TObjString *) nxfn())) {
            // These files need to be available everywhere, cache and sandbox
            if (SendFile(os->GetName(), kAscii | kForward, "cache") == -1) {
               Error("Load", "problems sending additional file %s", os->GetName());
               return -1;
            }
            // Add the base names to the message broadcasted
            bmsg += TString::Format(",%s", gSystem->BaseName(os->GetName()));
         }
         addfiles.SetOwner(kTRUE);
      }

      // The files are now on the workers: now we send the loading request
      TMessage mess(kPROOF_CACHE);
      if (GetRemoteProtocol() < 34) {
         mess << Int_t(kLoadMacro) << basemacro;
         // This may be needed
         AddIncludePath("../../cache");
      } else {
         mess << Int_t(kLoadMacro) << bmsg;
      }
      Broadcast(mess, kActive);

      // Load locally, if required
      if (!notOnClient) {
         // Mofify the include path
         TString oldincs = gSystem->GetIncludePath();
         if (!addincs.IsNull()) gSystem->AddIncludePath(addincs);

         // By first forwarding the load command to the master and workers
         // and only then loading locally we load/build in parallel
         gROOT->ProcessLine(TString::Format(".L %s", mainmacro.Data()));

         // Restore include path
         if (!addincs.IsNull()) gSystem->SetIncludePath(oldincs);

         // Update the macro path
         TString mp(TROOT::GetMacroPath());
         TString np = gSystem->GetDirName(macro);
         if (!np.IsNull()) {
            np += ":";
            if (!mp.BeginsWith(np) && !mp.Contains(":"+np)) {
               Int_t ip = (mp.BeginsWith(".:")) ? 2 : 0;
               mp.Insert(ip, np);
               TROOT::SetMacroPath(mp);
               if (gDebug > 0)
                  Info("Load", "macro path set to '%s'", TROOT::GetMacroPath());
            }
         }
      }

      // Wait for master and workers to be done
      Collect(kActive);

      if (IsLite()) {
         PDB(kGlobal, 1) Info("Load", "adding loaded macro: %s", macro);
         if (!fLoadedMacros) {
            fLoadedMacros = new TList();
            fLoadedMacros->SetOwner();
         }
         // if wrks is specified the macro should already be loaded on the master.
         fLoadedMacros->Add(new TObjString(macro));
      }

   } else {
      // On master

      // The files are now on the workers: now we send the loading request first
      // to the unique workers, so that the eventual compilation occurs only once.
      TString basemacro = gSystem->BaseName(macro);
      TMessage mess(kPROOF_CACHE);

      if (uniqueWorkers) {
         mess << Int_t(kLoadMacro) << basemacro;
         if (wrks) {
            Broadcast(mess, wrks);
            Collect(wrks);
         } else {
            Broadcast(mess, kUnique);
         }
      } else {
         // Wait for the result of the previous sending
         Collect(kUnique);

         // We then send a tuned loading request to the other workers
         TList others;
         TSlave *wrk = 0;
         TIter nxw(fActiveSlaves);
         while ((wrk = (TSlave *)nxw())) {
            if (!fUniqueSlaves->FindObject(wrk)) {
               others.Add(wrk);
            }
         }

         // Do not force compilation, if it was requested
         Int_t ld = basemacro.Last('.');
         if (ld != kNPOS) {
            Int_t lpp = basemacro.Index("++", ld);
            if (lpp != kNPOS) basemacro.Replace(lpp, 2, "+");
         }
         mess << Int_t(kLoadMacro) << basemacro;
         Broadcast(mess, &others);
         Collect(&others);
      }

      PDB(kGlobal, 1) Info("Load", "adding loaded macro: %s", macro);
      if (!fLoadedMacros) {
         fLoadedMacros = new TList();
         fLoadedMacros->SetOwner();
      }
      // if wrks is specified the macro should already be loaded on the master.
      if (!wrks)
         fLoadedMacros->Add(new TObjString(macro));
   }

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add 'libpath' to the lib path search.
/// Multiple paths can be specified at once separating them with a comma or
/// a blank.
/// Return 0 on success, -1 otherwise

Int_t TProof::AddDynamicPath(const char *libpath, Bool_t onClient, TList *wrks,
   Bool_t doCollect)
{
   if ((!libpath || !libpath[0])) {
      if (gDebug > 0)
         Info("AddDynamicPath", "list is empty - nothing to do");
      return 0;
   }

   // Do it also on clients, if required
   if (onClient)
      HandleLibIncPath("lib", kTRUE, libpath);

   TMessage m(kPROOF_LIB_INC_PATH);
   m << TString("lib") << (Bool_t)kTRUE;

   // Add paths
   if (libpath && strlen(libpath)) {
      m << TString(libpath);
   } else {
      m << TString("-");
   }

   // Tell the server to send back or not
   m << (Int_t)doCollect;

   // Forward the request
   if (wrks) {
      Broadcast(m, wrks);
      if (doCollect)
         Collect(wrks, fCollectTimeout);
   } else {
      Broadcast(m);
      Collect(kActive, fCollectTimeout);
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add 'incpath' to the inc path search.
/// Multiple paths can be specified at once separating them with a comma or
/// a blank.
/// Return 0 on success, -1 otherwise

Int_t TProof::AddIncludePath(const char *incpath, Bool_t onClient, TList *wrks,
   Bool_t doCollect)
{
   if ((!incpath || !incpath[0])) {
      if (gDebug > 0)
         Info("AddIncludePath", "list is empty - nothing to do");
      return 0;
   }

   // Do it also on clients, if required
   if (onClient)
      HandleLibIncPath("inc", kTRUE, incpath);

   TMessage m(kPROOF_LIB_INC_PATH);
   m << TString("inc") << (Bool_t)kTRUE;

   // Add paths
   if (incpath && strlen(incpath)) {
      m << TString(incpath);
   } else {
      m << TString("-");
   }

   // Tell the server to send back or not
   m << (Int_t)doCollect;

   // Forward the request
   if (wrks) {
      Broadcast(m, wrks);
      if (doCollect)
         Collect(wrks, fCollectTimeout);
   } else {
      Broadcast(m);
      Collect(kActive, fCollectTimeout);
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove 'libpath' from the lib path search.
/// Multiple paths can be specified at once separating them with a comma or
/// a blank.
/// Return 0 on success, -1 otherwise

Int_t TProof::RemoveDynamicPath(const char *libpath, Bool_t onClient)
{
   if ((!libpath || !libpath[0])) {
      if (gDebug > 0)
         Info("RemoveDynamicPath", "list is empty - nothing to do");
      return 0;
   }

   // Do it also on clients, if required
   if (onClient)
      HandleLibIncPath("lib", kFALSE, libpath);

   TMessage m(kPROOF_LIB_INC_PATH);
   m << TString("lib") <<(Bool_t)kFALSE;

   // Add paths
   if (libpath && strlen(libpath))
      m << TString(libpath);
   else
      m << TString("-");

   // Forward the request
   Broadcast(m);
   Collect(kActive, fCollectTimeout);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove 'incpath' from the inc path search.
/// Multiple paths can be specified at once separating them with a comma or
/// a blank.
/// Return 0 on success, -1 otherwise

Int_t TProof::RemoveIncludePath(const char *incpath, Bool_t onClient)
{
   if ((!incpath || !incpath[0])) {
      if (gDebug > 0)
         Info("RemoveIncludePath", "list is empty - nothing to do");
      return 0;
   }

   // Do it also on clients, if required
   if (onClient)
      HandleLibIncPath("in", kFALSE, incpath);

   TMessage m(kPROOF_LIB_INC_PATH);
   m << TString("inc") << (Bool_t)kFALSE;

   // Add paths
   if (incpath && strlen(incpath))
      m << TString(incpath);
   else
      m << TString("-");

   // Forward the request
   Broadcast(m);
   Collect(kActive, fCollectTimeout);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle lib, inc search paths modification request

void TProof::HandleLibIncPath(const char *what, Bool_t add, const char *dirs)
{
   TString type(what);
   TString path(dirs);

   // Check type of action
   if ((type != "lib") && (type != "inc")) {
      Error("HandleLibIncPath","unknown action type: %s - protocol error?", type.Data());
      return;
   }

   // Separators can be either commas or blanks
   path.ReplaceAll(","," ");

   // Decompose lists
   TObjArray *op = 0;
   if (path.Length() > 0 && path != "-") {
      if (!(op = path.Tokenize(" "))) {
         Warning("HandleLibIncPath","decomposing path %s", path.Data());
         return;
      }
   }

   if (add) {

      if (type == "lib") {

         // Add libs
         TIter nxl(op, kIterBackward);
         TObjString *lib = 0;
         while ((lib = (TObjString *) nxl())) {
            // Expand path
            TString xlib = lib->GetName();
            gSystem->ExpandPathName(xlib);
            // Add to the dynamic lib search path if it exists and can be read
            if (!gSystem->AccessPathName(xlib, kReadPermission)) {
               TString newlibpath = gSystem->GetDynamicPath();
               // In the first position after the working dir
               Int_t pos = 0;
               if (newlibpath.BeginsWith(".:"))
                  pos = 2;
               if (newlibpath.Index(xlib) == kNPOS) {
                  newlibpath.Insert(pos,TString::Format("%s:", xlib.Data()));
                  gSystem->SetDynamicPath(newlibpath);
               }
            } else {
               if (gDebug > 0)
                  Info("HandleLibIncPath",
                       "libpath %s does not exist or cannot be read - not added", xlib.Data());
            }
         }

      } else {

         // Add incs
         TIter nxi(op);
         TObjString *inc = 0;
         while ((inc = (TObjString *) nxi())) {
            // Expand path
            TString xinc = inc->GetName();
            gSystem->ExpandPathName(xinc);
            // Add to the dynamic lib search path if it exists and can be read
            if (!gSystem->AccessPathName(xinc, kReadPermission)) {
               TString curincpath = gSystem->GetIncludePath();
               if (curincpath.Index(xinc) == kNPOS)
                  gSystem->AddIncludePath(TString::Format("-I%s", xinc.Data()));
            } else
               if (gDebug > 0)
                   Info("HandleLibIncPath",
                        "incpath %s does not exist or cannot be read - not added", xinc.Data());
         }
      }


   } else {

      if (type == "lib") {

         // Remove libs
         TIter nxl(op);
         TObjString *lib = 0;
         while ((lib = (TObjString *) nxl())) {
            // Expand path
            TString xlib = lib->GetName();
            gSystem->ExpandPathName(xlib);
            // Remove from the dynamic lib search path
            TString newlibpath = gSystem->GetDynamicPath();
            newlibpath.ReplaceAll(TString::Format("%s:", xlib.Data()),"");
            gSystem->SetDynamicPath(newlibpath);
         }

      } else {

         // Remove incs
         TIter nxi(op);
         TObjString *inc = 0;
         while ((inc = (TObjString *) nxi())) {
            TString newincpath = gSystem->GetIncludePath();
            newincpath.ReplaceAll(TString::Format("-I%s", inc->GetName()),"");
            // Remove the interpreter path (added anyhow internally)
            newincpath.ReplaceAll(gInterpreter->GetIncludePath(),"");
            gSystem->SetIncludePath(newincpath);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get from the master the list of names of the packages available.

TList *TProof::GetListOfPackages()
{
   if (!IsValid())
      return (TList *)0;

   TMessage mess(kPROOF_CACHE);
   mess << Int_t(kListPackages);
   Broadcast(mess);
   Collect(kActive, fCollectTimeout);

   return fAvailablePackages;
}

////////////////////////////////////////////////////////////////////////////////
/// Get from the master the list of names of the packages enabled.

TList *TProof::GetListOfEnabledPackages()
{
   if (!IsValid())
      return (TList *)0;

   TMessage mess(kPROOF_CACHE);
   mess << Int_t(kListEnabledPackages);
   Broadcast(mess);
   Collect(kActive, fCollectTimeout);

   return fEnabledPackages;
}

////////////////////////////////////////////////////////////////////////////////
/// Print a progress bar on stderr. Used in batch mode.

void TProof::PrintProgress(Long64_t total, Long64_t processed,
                           Float_t procTime, Long64_t bytesread)
{
   if (fPrintProgress) {
      Bool_t redirlog = fRedirLog;
      fRedirLog = kFALSE;
      // Call the external function
      (*fPrintProgress)(total, processed, procTime, bytesread);
      fRedirLog = redirlog;
      return;
   }

   fprintf(stderr, "[TProof::Progress] Total %lld events\t|", total);

   for (int l = 0; l < 20; l++) {
      if (total > 0) {
         if (l < 20*processed/total)
            fprintf(stderr, "=");
         else if (l == 20*processed/total)
            fprintf(stderr, ">");
         else if (l > 20*processed/total)
            fprintf(stderr, ".");
      } else
         fprintf(stderr, "=");
   }
   Float_t evtrti = (procTime > 0. && processed > 0) ? processed / procTime : -1.;
   Float_t mbsrti = (procTime > 0. && bytesread > 0) ? bytesread / procTime : -1.;
   TString sunit("B/s");
   if (evtrti > 0.) {
      Float_t remainingTime = (total >= processed) ? (total - processed) / evtrti : -1;
      if (mbsrti > 0.) {
         const Float_t toK = 1024., toM = 1048576., toG = 1073741824.;
         if (mbsrti >= toG) {
            mbsrti /= toG;
            sunit = "GB/s";
         } else if (mbsrti >= toM) {
            mbsrti /= toM;
            sunit = "MB/s";
         } else if (mbsrti >= toK) {
            mbsrti /= toK;
            sunit = "kB/s";
         }
         fprintf(stderr, "| %.02f %% [%.1f evts/s, %.1f %s, time left: %.1f s]\r",
                (total ? ((100.0*processed)/total) : 100.0), evtrti, mbsrti, sunit.Data(), remainingTime);
      } else {
         fprintf(stderr, "| %.02f %% [%.1f evts/s, time left: %.1f s]\r",
                (total ? ((100.0*processed)/total) : 100.0), evtrti, remainingTime);
      }
   } else {
      fprintf(stderr, "| %.02f %%\r",
              (total ? ((100.0*processed)/total) : 100.0));
   }
   if (processed >= total) {
      fprintf(stderr, "\n Query processing time: %.1f s\n", procTime);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get query progress information. Connect a slot to this signal
/// to track progress.

void TProof::Progress(Long64_t total, Long64_t processed)
{
   if (fPrintProgress) {
      // Call the external function
      return (*fPrintProgress)(total, processed, -1., -1);
   }

   PDB(kGlobal,1)
      Info("Progress","%2f (%lld/%lld)", 100.*processed/total, processed, total);

   if (gROOT->IsBatch()) {
      // Simple progress bar
      if (total > 0)
         PrintProgress(total, processed);
   } else {
      EmitVA("Progress(Long64_t,Long64_t)", 2, total, processed);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get query progress information. Connect a slot to this signal
/// to track progress.

void TProof::Progress(Long64_t total, Long64_t processed, Long64_t bytesread,
                      Float_t initTime, Float_t procTime,
                      Float_t evtrti, Float_t mbrti)
{
   PDB(kGlobal,1)
      Info("Progress","%lld %lld %lld %f %f %f %f", total, processed, bytesread,
                                initTime, procTime, evtrti, mbrti);

   if (gROOT->IsBatch()) {
      // Simple progress bar
      if (total > 0)
         PrintProgress(total, processed, procTime, bytesread);
   } else {
      EmitVA("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)",
             7, total, processed, bytesread, initTime, procTime, evtrti, mbrti);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get query progress information. Connect a slot to this signal
/// to track progress.

void TProof::Progress(Long64_t total, Long64_t processed, Long64_t bytesread,
                      Float_t initTime, Float_t procTime,
                      Float_t evtrti, Float_t mbrti, Int_t actw, Int_t tses, Float_t eses)
{
   PDB(kGlobal,1)
      Info("Progress","%lld %lld %lld %f %f %f %f %d %f", total, processed, bytesread,
                                initTime, procTime, evtrti, mbrti, actw, eses);

   if (gROOT->IsBatch()) {
      // Simple progress bar
      if (total > 0)
         PrintProgress(total, processed, procTime, bytesread);
   } else {
      EmitVA("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)",
             10, total, processed, bytesread, initTime, procTime, evtrti, mbrti, actw, tses, eses);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get list of feedback objects. Connect a slot to this signal
/// to monitor the feedback object.

void TProof::Feedback(TList *objs)
{
   PDB(kGlobal,1)
      Info("Feedback","%d objects", objs->GetSize());
   PDB(kFeedback,1) {
      Info("Feedback","%d objects", objs->GetSize());
      objs->ls();
   }

   Emit("Feedback(TList *objs)", (Long_t) objs);
}

////////////////////////////////////////////////////////////////////////////////
/// Close progress dialog.

void TProof::CloseProgressDialog()
{
   PDB(kGlobal,1)
      Info("CloseProgressDialog",
           "called: have progress dialog: %d", fProgressDialogStarted);

   // Nothing to do if not there
   if (!fProgressDialogStarted)
      return;

   Emit("CloseProgressDialog()");
}

////////////////////////////////////////////////////////////////////////////////
/// Reset progress dialog.

void TProof::ResetProgressDialog(const char *sel, Int_t sz, Long64_t fst,
                                 Long64_t ent)
{
   PDB(kGlobal,1)
      Info("ResetProgressDialog","(%s,%d,%lld,%lld)", sel, sz, fst, ent);

   EmitVA("ResetProgressDialog(const char*,Int_t,Long64_t,Long64_t)",
          4, sel, sz, fst, ent);
}

////////////////////////////////////////////////////////////////////////////////
/// Send startup message.

void TProof::StartupMessage(const char *msg, Bool_t st, Int_t done, Int_t total)
{
   PDB(kGlobal,1)
      Info("StartupMessage","(%s,%d,%d,%d)", msg, st, done, total);

   EmitVA("StartupMessage(const char*,Bool_t,Int_t,Int_t)",
          4, msg, st, done, total);
}

////////////////////////////////////////////////////////////////////////////////
/// Send dataset preparation status.

void TProof::DataSetStatus(const char *msg, Bool_t st, Int_t done, Int_t total)
{
   PDB(kGlobal,1)
      Info("DataSetStatus","(%s,%d,%d,%d)", msg, st, done, total);

   EmitVA("DataSetStatus(const char*,Bool_t,Int_t,Int_t)",
          4, msg, st, done, total);
}

////////////////////////////////////////////////////////////////////////////////
/// Send or notify data set status

void TProof::SendDataSetStatus(const char *action, UInt_t done,
                               UInt_t tot, Bool_t st)
{
   if (IsLite()) {
      if (tot) {
         TString type = "files";
         Int_t frac = (Int_t) (done*100.)/tot;
         char msg[512] = {0};
         if (frac >= 100) {
            snprintf(msg, 512, "%s: OK (%d %s)                 \n",
                     action,tot, type.Data());
         } else {
            snprintf(msg, 512, "%s: %d out of %d (%d %%)\r",
                     action, done, tot, frac);
         }
         if (fSync)
            fprintf(stderr,"%s", msg);
         else
            NotifyLogMsg(msg, 0);
      }
      return;
   }

   if (TestBit(TProof::kIsMaster)) {
      TMessage mess(kPROOF_DATASET_STATUS);
      mess << TString(action) << tot << done << st;
      gProofServ->GetSocket()->Send(mess);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Notify availability of a query result.

void TProof::QueryResultReady(const char *ref)
{
   PDB(kGlobal,1)
      Info("QueryResultReady","ref: %s", ref);

   Emit("QueryResultReady(const char*)",ref);
}

////////////////////////////////////////////////////////////////////////////////
/// Validate a TDSet.

void TProof::ValidateDSet(TDSet *dset)
{
   if (dset->ElementsValid()) return;

   TList nodes;
   nodes.SetOwner();

   TList slholder;
   slholder.SetOwner();
   TList elemholder;
   elemholder.SetOwner();

   // build nodelist with slaves and elements
   TIter nextSlave(GetListOfActiveSlaves());
   while (TSlave *sl = dynamic_cast<TSlave*>(nextSlave())) {
      TList *sllist = 0;
      TPair *p = dynamic_cast<TPair*>(nodes.FindObject(sl->GetName()));
      if (!p) {
         sllist = new TList;
         sllist->SetName(sl->GetName());
         slholder.Add(sllist);
         TList *elemlist = new TList;
         elemlist->SetName(TString(sl->GetName())+"_elem");
         elemholder.Add(elemlist);
         nodes.Add(new TPair(sllist, elemlist));
      } else {
         sllist = dynamic_cast<TList*>(p->Key());
      }
      if (sllist) sllist->Add(sl);
   }

   // add local elements to nodes
   TList nonLocal; // list of nonlocal elements
   // make two iterations - first add local elements - then distribute nonlocals
   for (Int_t i = 0; i < 2; i++) {
      Bool_t local = i>0?kFALSE:kTRUE;
      TIter nextElem(local ? dset->GetListOfElements() : &nonLocal);
      while (TDSetElement *elem = dynamic_cast<TDSetElement*>(nextElem())) {
         if (elem->GetValid()) continue;
         TPair *p = dynamic_cast<TPair*>(local?nodes.FindObject(TUrl(elem->GetFileName()).GetHost()):nodes.At(0));
         if (p) {
            TList *eli = dynamic_cast<TList*>(p->Value());
            TList *sli = dynamic_cast<TList*>(p->Key());
            if (eli && sli) {
               eli->Add(elem);

               // order list by elements/slave
               TPair *p2 = p;
               Bool_t stop = kFALSE;
               while (!stop) {
                  TPair *p3 = dynamic_cast<TPair*>(nodes.After(p2->Key()));
                  if (p3) {
                     TList *p3v = dynamic_cast<TList*>(p3->Value());
                     TList *p3k = dynamic_cast<TList*>(p3->Key());
                     if (p3v && p3k) {
                        Int_t nelem = p3v->GetSize();
                        Int_t nsl = p3k->GetSize();
                        if (nelem*sli->GetSize() < eli->GetSize()*nsl) p2 = p3;
                        else stop = kTRUE;
                     }
                  } else {
                     stop = kTRUE;
                  }
               }

               if (p2!=p) {
                  nodes.Remove(p->Key());
                  nodes.AddAfter(p2->Key(), p);
               }
            } else {
               Warning("ValidateDSet", "invalid values from TPair! Protocol error?");
               continue;
            }

         } else {
            if (local) {
               nonLocal.Add(elem);
            } else {
               Warning("ValidateDSet", "no node to allocate TDSetElement to - ignoring");
            }
         }
      }
   }

   // send to slaves
   TList usedslaves;
   TIter nextNode(&nodes);
   SetDSet(dset); // set dset to be validated in Collect()
   while (TPair *node = dynamic_cast<TPair*>(nextNode())) {
      TList *slaves = dynamic_cast<TList*>(node->Key());
      TList *setelements = dynamic_cast<TList*>(node->Value());
      if (!slaves || !setelements) continue;
      // distribute elements over the slaves
      Int_t nslaves = slaves->GetSize();
      Int_t nelements = setelements->GetSize();
      for (Int_t i=0; i<nslaves; i++) {

         TDSet copyset(dset->GetType(), dset->GetObjName(),
                       dset->GetDirectory());
         for (Int_t j = (i*nelements)/nslaves;
                    j < ((i+1)*nelements)/nslaves;
                    j++) {
            TDSetElement *elem =
               dynamic_cast<TDSetElement*>(setelements->At(j));
            if (elem) {
               copyset.Add(elem->GetFileName(), elem->GetObjName(),
                           elem->GetDirectory(), elem->GetFirst(),
                           elem->GetNum(), elem->GetMsd());
            }
         }

         if (copyset.GetListOfElements()->GetSize()>0) {
            TMessage mesg(kPROOF_VALIDATE_DSET);
            mesg << &copyset;

            TSlave *sl = dynamic_cast<TSlave*>(slaves->At(i));
            if (sl) {
               PDB(kGlobal,1) Info("ValidateDSet",
                                 "Sending TDSet with %d elements to slave %s"
                                 " to be validated",
                                 copyset.GetListOfElements()->GetSize(),
                                 sl->GetOrdinal());
               sl->GetSocket()->Send(mesg);
               usedslaves.Add(sl);
            }
         }
      }
   }

   PDB(kGlobal,1)
      Info("ValidateDSet","Calling Collect");
   Collect(&usedslaves);
   SetDSet(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Add data objects that might be needed during the processing of
/// the selector (see Process()). This object can be very large, so they
/// are distributed in an optimized way using a dedicated file.
/// If push is TRUE the input data are sent over even if no apparent change
/// occured to the list.

void TProof::AddInputData(TObject *obj, Bool_t push)
{
   if (obj) {
      if (!fInputData) fInputData = new TList;
      if (!fInputData->FindObject(obj)) {
         fInputData->Add(obj);
         SetBit(TProof::kNewInputData);
      }
   }
   if (push) SetBit(TProof::kNewInputData);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove obj form the input data list; if obj is null (default), clear the
/// input data info.

void TProof::ClearInputData(TObject *obj)
{
   if (!obj) {
      if (fInputData) {
         fInputData->SetOwner(kTRUE);
         SafeDelete(fInputData);
      }
      ResetBit(TProof::kNewInputData);

      // Also remove any info about input data in the input list
      TObject *o = 0;
      TList *in = GetInputList();
      while ((o = GetInputList()->FindObject("PROOF_InputDataFile")))
         in->Remove(o);
      while ((o = GetInputList()->FindObject("PROOF_InputData")))
         in->Remove(o);

      // ... and reset the file
      fInputDataFile = "";
      gSystem->Unlink(kPROOF_InputDataFile);

   } else if (fInputData) {
      Int_t sz = fInputData->GetSize();
      while (fInputData->FindObject(obj))
         fInputData->Remove(obj);
      // Flag for update, if anything changed
      if (sz != fInputData->GetSize())
         SetBit(TProof::kNewInputData);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove obj 'name' form the input data list;

void TProof::ClearInputData(const char *name)
{
   TObject *obj = (fInputData && name) ? fInputData->FindObject(name) : 0;
   if (obj) ClearInputData(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the file to be used to optimally distribute the input data objects.
/// If the file exists the object in the file are added to those in the
/// fInputData list. If the file path is null, a default file will be created
/// at the moment of sending the processing request with the content of
/// the fInputData list. See also SendInputDataFile.

void TProof::SetInputDataFile(const char *datafile)
{
   if (datafile && strlen(datafile) > 0) {
      if (fInputDataFile != datafile && strcmp(datafile, kPROOF_InputDataFile))
         SetBit(TProof::kNewInputData);
      fInputDataFile = datafile;
   } else {
      if (!fInputDataFile.IsNull())
         SetBit(TProof::kNewInputData);
      fInputDataFile = "";
   }
   // Make sure that the chosen file is readable
   if (fInputDataFile != kPROOF_InputDataFile && !fInputDataFile.IsNull() &&
      gSystem->AccessPathName(fInputDataFile, kReadPermission)) {
      fInputDataFile = "";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Send the input data objects to the master; the objects are taken from the
/// dedicated list and / or the specified file.
/// If the fInputData is empty the specified file is sent over.
/// If there is no specified file, a file named "inputdata.root" is created locally
/// with the content of fInputData and sent over to the master.
/// If both fInputData and the specified file are not empty, a copy of the file
/// is made locally and augmented with the content of fInputData.

void TProof::SendInputDataFile()
{
   // Prepare the file
   TString dataFile;
   PrepareInputDataFile(dataFile);

   // Send it, if not empty
   if (dataFile.Length() > 0) {

      Info("SendInputDataFile", "broadcasting %s", dataFile.Data());
      BroadcastFile(dataFile.Data(), kBinary, "cache", kActive);

      // Set the name in the input list
      TString t = TString::Format("cache:%s", gSystem->BaseName(dataFile));
      AddInput(new TNamed("PROOF_InputDataFile", t.Data()));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Prepare the file with the input data objects to be sent the master; the
/// objects are taken from the dedicated list and / or the specified file.
/// If the fInputData is empty the specified file is sent over.
/// If there is no specified file, a file named "inputdata.root" is created locally
/// with the content of fInputData and sent over to the master.
/// If both fInputData and the specified file are not empty, a copy of the file
/// is made locally and augmented with the content of fInputData.

void TProof::PrepareInputDataFile(TString &dataFile)
{
   // Save info about new data for usage in this call;
   Bool_t newdata = TestBit(TProof::kNewInputData) ? kTRUE : kFALSE;
   // Next time we need some change
   ResetBit(TProof::kNewInputData);

   // Check the list
   Bool_t list_ok = (fInputData && fInputData->GetSize() > 0) ? kTRUE : kFALSE;
   // Check the file
   Bool_t file_ok = kFALSE;
   if (fInputDataFile != kPROOF_InputDataFile && !fInputDataFile.IsNull() &&
      !gSystem->AccessPathName(fInputDataFile, kReadPermission)) {
      // It must contain something
      TFile *f = TFile::Open(fInputDataFile);
      if (f && f->GetListOfKeys() && f->GetListOfKeys()->GetSize() > 0)
         file_ok = kTRUE;
   }

   // Remove any info about input data in the input list
   TObject *o = 0;
   TList *in = GetInputList();
   while ((o = GetInputList()->FindObject("PROOF_InputDataFile")))
      in->Remove(o);
   while ((o = GetInputList()->FindObject("PROOF_InputData")))
      in->Remove(o);

   // We must have something to send
   dataFile = "";
   if (!list_ok && !file_ok) return;

   // Three cases:
   if (file_ok && !list_ok) {
      // Just send the file
      dataFile = fInputDataFile;
   } else if (!file_ok && list_ok) {
      fInputDataFile = kPROOF_InputDataFile;
      // Nothing to do, if no new data
      if (!newdata && !gSystem->AccessPathName(fInputDataFile)) return;
      // Create the file first
      TFile *f = TFile::Open(fInputDataFile, "RECREATE");
      if (f) {
         f->cd();
         TIter next(fInputData);
         TObject *obj;
         while ((obj = next())) {
            obj->Write(0, TObject::kSingleKey, 0);
         }
         f->Close();
         SafeDelete(f);
      } else {
         Error("PrepareInputDataFile", "could not (re-)create %s", fInputDataFile.Data());
         return;
      }
      dataFile = fInputDataFile;
   } else if (file_ok && list_ok) {
      dataFile = kPROOF_InputDataFile;
      // Create the file if not existing or there are new data
      if (newdata || gSystem->AccessPathName(dataFile)) {
         // Cleanup previous file if obsolete
         if (!gSystem->AccessPathName(dataFile))
            gSystem->Unlink(dataFile);
         if (dataFile != fInputDataFile) {
            // Make a local copy first
            if (gSystem->CopyFile(fInputDataFile, dataFile, kTRUE) != 0) {
               Error("PrepareInputDataFile", "could not make local copy of %s", fInputDataFile.Data());
               return;
            }
         }
         // Add the input data list
         TFile *f = TFile::Open(dataFile, "UPDATE");
         if (f) {
            f->cd();
            TIter next(fInputData);
            TObject *obj = 0;
            while ((obj = next())) {
               obj->Write(0, TObject::kSingleKey, 0);
            }
            f->Close();
            SafeDelete(f);
         } else {
            Error("PrepareInputDataFile", "could not open %s for updating", dataFile.Data());
            return;
         }
      }
   }

   //  Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Add objects that might be needed during the processing of
/// the selector (see Process()).

void TProof::AddInput(TObject *obj)
{
   if (fPlayer) fPlayer->AddInput(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Clear input object list.

void TProof::ClearInput()
{
   if (fPlayer) fPlayer->ClearInput();

   // the system feedback list is always in the input list
   AddInput(fFeedback);
}

////////////////////////////////////////////////////////////////////////////////
/// Get input list.

TList *TProof::GetInputList()
{
   return (fPlayer ? fPlayer->GetInputList() : (TList *)0);
}

////////////////////////////////////////////////////////////////////////////////
/// Get specified object that has been produced during the processing
/// (see Process()).

TObject *TProof::GetOutput(const char *name)
{

   if (TestBit(TProof::kIsMaster))
      // Can be called by MarkBad on the master before the player is initialized
      return (fPlayer) ? fPlayer->GetOutput(name) : (TObject *)0;

   // This checks also associated output files
   return (GetOutputList()) ? GetOutputList()->FindObject(name) : (TObject *)0;
}

////////////////////////////////////////////////////////////////////////////////
/// Find object 'name' in list 'out' or in the files specified in there

TObject *TProof::GetOutput(const char *name, TList *out)
{
   TObject *o = 0;
   if (!name || (name && strlen(name) <= 0) ||
       !out || (out && out->GetSize() <= 0)) return o;
   if ((o = out->FindObject(name))) return o;

   // For the time being we always check for all the files; this may require
   // some caching
   TProofOutputFile *pf = 0;
   TIter nxo(out);
   while ((o = nxo())) {
      if ((pf = dynamic_cast<TProofOutputFile *> (o))) {
         TFile *f = 0;
         if (!(f = (TFile *) gROOT->GetListOfFiles()->FindObject(pf->GetOutputFileName()))) {
            TString fn = TString::Format("%s/%s", pf->GetDir(), pf->GetFileName());
            f = TFile::Open(fn.Data());
            if (!f || (f && f->IsZombie())) {
               ::Warning("TProof::GetOutput", "problems opening file %s", fn.Data());
            }
         }
         if (f && (o = f->Get(name))) return o;
      }
   }

   // Done, unsuccessfully
   return o;
}

////////////////////////////////////////////////////////////////////////////////
/// Get list with all object created during processing (see Process()).

TList *TProof::GetOutputList()
{
   if (fOutputList.GetSize() > 0) return &fOutputList;
   if (fPlayer) {
      fOutputList.AttachList(fPlayer->GetOutputList());
      return &fOutputList;
   }
   return (TList *)0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set input list parameter. If the parameter is already
/// set it will be set to the new value.

void TProof::SetParameter(const char *par, const char *value)
{
   if (!fPlayer) {
      Warning("SetParameter", "player undefined! Ignoring");
      return;
   }

   TList *il = fPlayer->GetInputList();
   TObject *item = il->FindObject(par);
   if (item) {
      il->Remove(item);
      delete item;
   }
   il->Add(new TNamed(par, value));
}

////////////////////////////////////////////////////////////////////////////////
/// Set an input list parameter.

void TProof::SetParameter(const char *par, Int_t value)
{
   if (!fPlayer) {
      Warning("SetParameter", "player undefined! Ignoring");
      return;
   }

   TList *il = fPlayer->GetInputList();
   TObject *item = il->FindObject(par);
   if (item) {
      il->Remove(item);
      delete item;
   }
   il->Add(new TParameter<Int_t>(par, value));
}

////////////////////////////////////////////////////////////////////////////////
/// Set an input list parameter.

void TProof::SetParameter(const char *par, Long_t value)
{
   if (!fPlayer) {
      Warning("SetParameter", "player undefined! Ignoring");
      return;
   }

   TList *il = fPlayer->GetInputList();
   TObject *item = il->FindObject(par);
   if (item) {
      il->Remove(item);
      delete item;
   }
   il->Add(new TParameter<Long_t>(par, value));
}

////////////////////////////////////////////////////////////////////////////////
/// Set an input list parameter.

void TProof::SetParameter(const char *par, Long64_t value)
{
   if (!fPlayer) {
      Warning("SetParameter", "player undefined! Ignoring");
      return;
   }

   TList *il = fPlayer->GetInputList();
   TObject *item = il->FindObject(par);
   if (item) {
      il->Remove(item);
      delete item;
   }
   il->Add(new TParameter<Long64_t>(par, value));
}

////////////////////////////////////////////////////////////////////////////////
/// Set an input list parameter.

void TProof::SetParameter(const char *par, Double_t value)
{
   if (!fPlayer) {
      Warning("SetParameter", "player undefined! Ignoring");
      return;
   }

   TList *il = fPlayer->GetInputList();
   TObject *item = il->FindObject(par);
   if (item) {
      il->Remove(item);
      delete item;
   }
   il->Add(new TParameter<Double_t>(par, value));
}

////////////////////////////////////////////////////////////////////////////////
/// Get specified parameter. A parameter set via SetParameter() is either
/// a TParameter or a TNamed or 0 in case par is not defined.

TObject *TProof::GetParameter(const char *par) const
{
   if (!fPlayer) {
      Warning("GetParameter", "player undefined! Ignoring");
      return (TObject *)0;
   }

   TList *il = fPlayer->GetInputList();
   return il->FindObject(par);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete the input list parameters specified by a wildcard (e.g. PROOF_*)
/// or exact name (e.g. PROOF_MaxSlavesPerNode).

void TProof::DeleteParameters(const char *wildcard)
{
   if (!fPlayer) return;

   if (!wildcard) wildcard = "";
   TRegexp re(wildcard, kTRUE);
   Int_t nch = strlen(wildcard);

   TList *il = fPlayer->GetInputList();
   if (il) {
      TObject *p = 0;
      TIter next(il);
      while ((p = next())) {
         TString s = p->GetName();
         if (nch && s != wildcard && s.Index(re) == kNPOS) continue;
         il->Remove(p);
         delete p;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Show the input list parameters specified by the wildcard.
/// Default is the special PROOF control parameters (PROOF_*).

void TProof::ShowParameters(const char *wildcard) const
{
   if (!fPlayer) return;

   if (!wildcard) wildcard = "";
   TRegexp re(wildcard, kTRUE);
   Int_t nch = strlen(wildcard);

   TList *il = fPlayer->GetInputList();
   TObject *p;
   TIter next(il);
   while ((p = next())) {
      TString s = p->GetName();
      if (nch && s != wildcard && s.Index(re) == kNPOS) continue;
      if (p->IsA() == TNamed::Class()) {
         Printf("%s\t\t\t%s", s.Data(), p->GetTitle());
      } else if (p->IsA() == TParameter<Long_t>::Class()) {
         Printf("%s\t\t\t%ld", s.Data(), dynamic_cast<TParameter<Long_t>*>(p)->GetVal());
      } else if (p->IsA() == TParameter<Long64_t>::Class()) {
         Printf("%s\t\t\t%lld", s.Data(), dynamic_cast<TParameter<Long64_t>*>(p)->GetVal());
      } else if (p->IsA() == TParameter<Double_t>::Class()) {
         Printf("%s\t\t\t%f", s.Data(), dynamic_cast<TParameter<Double_t>*>(p)->GetVal());
      } else {
         Printf("%s\t\t\t%s", s.Data(), p->GetTitle());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add object to feedback list.

void TProof::AddFeedback(const char *name)
{
   PDB(kFeedback, 3)
      Info("AddFeedback", "Adding object \"%s\" to feedback", name);
   if (fFeedback->FindObject(name) == 0)
      fFeedback->Add(new TObjString(name));
}

////////////////////////////////////////////////////////////////////////////////
/// Remove object from feedback list.

void TProof::RemoveFeedback(const char *name)
{
   TObject *obj = fFeedback->FindObject(name);
   if (obj != 0) {
      fFeedback->Remove(obj);
      delete obj;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clear feedback list.

void TProof::ClearFeedback()
{
   fFeedback->Delete();
}

////////////////////////////////////////////////////////////////////////////////
/// Show items in feedback list.

void TProof::ShowFeedback() const
{
   if (fFeedback->GetSize() == 0) {
      Info("","no feedback requested");
      return;
   }

   fFeedback->Print();
}

////////////////////////////////////////////////////////////////////////////////
/// Return feedback list.

TList *TProof::GetFeedbackList() const
{
   return fFeedback;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a tree header (a tree with nonexisting files) object for
/// the DataSet.

TTree *TProof::GetTreeHeader(TDSet *dset)
{
   TList *l = GetListOfActiveSlaves();
   TSlave *sl = (TSlave*) l->First();
   if (sl == 0) {
      Error("GetTreeHeader", "No connection");
      return 0;
   }

   TSocket *soc = sl->GetSocket();
   TMessage msg(kPROOF_GETTREEHEADER);

   msg << dset;

   soc->Send(msg);

   TMessage *reply;
   Int_t d = -1;
   if (fProtocol >= 20) {
      Collect(sl, fCollectTimeout, kPROOF_GETTREEHEADER);
      reply = (TMessage *) fRecvMessages->First();
   } else {
      d = soc->Recv(reply);
   }
   if (!reply) {
      Error("GetTreeHeader", "Error getting a replay from the master.Result %d", (int) d);
      return 0;
   }

   TString s1;
   TTree *t = 0;
   (*reply) >> s1;
   if (s1 == "Success")
      (*reply) >> t;

   PDB(kGlobal, 1) {
      if (t) {
         Info("GetTreeHeader", "%s, message size: %d, entries: %d",
                               s1.Data(), reply->BufferSize(), (int) t->GetMaxEntryLoop());
      } else {
         Info("GetTreeHeader", "tree header retrieval failed");
      }
   }
   delete reply;

   return t;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw feedback creation proxy. When accessed via TProof avoids
/// link dependency on libProofPlayer.

TDrawFeedback *TProof::CreateDrawFeedback()
{
   return (fPlayer ? fPlayer->CreateDrawFeedback(this) : (TDrawFeedback *)0);
}

////////////////////////////////////////////////////////////////////////////////
/// Set draw feedback option.

void TProof::SetDrawFeedbackOption(TDrawFeedback *f, Option_t *opt)
{
   if (fPlayer) fPlayer->SetDrawFeedbackOption(f, opt);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete draw feedback object.

void TProof::DeleteDrawFeedback(TDrawFeedback *f)
{
   if (fPlayer) fPlayer->DeleteDrawFeedback(f);
}

////////////////////////////////////////////////////////////////////////////////
///   FIXME: to be written

TList *TProof::GetOutputNames()
{
   return 0;
/*
   TMessage msg(kPROOF_GETOUTPUTLIST);
   TList* slaves = fActiveSlaves;
   Broadcast(msg, slaves);
   TMonitor mon;
   TList* outputList = new TList();

   TIter    si(slaves);
   TSlave   *slave;
   while ((slave = (TSlave*)si.Next()) != 0) {
      PDB(kGlobal,4) Info("GetOutputNames","Socket added to monitor: %p (%s)",
          slave->GetSocket(), slave->GetName());
      mon.Add(slave->GetSocket());
   }
   mon.ActivateAll();
   ((TProof*)gProof)->DeActivateAsyncInput();
   ((TProof*)gProof)->fCurrentMonitor = &mon;

   while (mon.GetActive() != 0) {
      TSocket *sock = mon.Select();
      if (!sock) {
         Error("GetOutputList","TMonitor::.Select failed!");
         break;
      }
      mon.DeActivate(sock);
      TMessage *reply;
      if (sock->Recv(reply) <= 0) {
         MarkBad(slave, "receive failed after kPROOF_GETOUTPUTLIST request");
//         Error("GetOutputList","Recv failed! for slave-%d (%s)",
//               slave->GetOrdinal(), slave->GetName());
         continue;
      }
      if (reply->What() != kPROOF_GETOUTPUTNAMES ) {
//         Error("GetOutputList","unexpected message %d from slawe-%d (%s)",  reply->What(),
//               slave->GetOrdinal(), slave->GetName());
         MarkBad(slave, "wrong reply to kPROOF_GETOUTPUTLIST request");
         continue;
      }
      TList* l;

      (*reply) >> l;
      TIter next(l);
      TNamed *n;
      while ( (n = dynamic_cast<TNamed*> (next())) ) {
         if (!outputList->FindObject(n->GetName()))
            outputList->Add(n);
      }
      delete reply;
   }
   ((TProof*)gProof)->fCurrentMonitor = 0;

   return outputList;
*/
}

////////////////////////////////////////////////////////////////////////////////
/// Build the PROOF's structure in the browser.

void TProof::Browse(TBrowser *b)
{
   b->Add(fActiveSlaves, fActiveSlaves->Class(), "fActiveSlaves");
   b->Add(&fMaster, fMaster.Class(), "fMaster");
   b->Add(fFeedback, fFeedback->Class(), "fFeedback");
   b->Add(fChains, fChains->Class(), "fChains");

   if (fPlayer) {
      b->Add(fPlayer->GetInputList(), fPlayer->GetInputList()->Class(), "InputList");
      if (fPlayer->GetOutputList())
         b->Add(fPlayer->GetOutputList(), fPlayer->GetOutputList()->Class(), "OutputList");
      if (fPlayer->GetListOfResults())
         b->Add(fPlayer->GetListOfResults(),
               fPlayer->GetListOfResults()->Class(), "ListOfResults");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set a new PROOF player.

void TProof::SetPlayer(TVirtualProofPlayer *player)
{
   if (fPlayer)
      delete fPlayer;
   fPlayer = player;
};

////////////////////////////////////////////////////////////////////////////////
/// Construct a TProofPlayer object. The player string specifies which
/// player should be created: remote, slave, sm (supermaster) or base.
/// Default is remote. Socket is needed in case a slave player is created.

TVirtualProofPlayer *TProof::MakePlayer(const char *player, TSocket *s)
{
   if (!player)
      player = "remote";

   SetPlayer(TVirtualProofPlayer::Create(player, this, s));
   return GetPlayer();
}

////////////////////////////////////////////////////////////////////////////////
/// Add chain to data set

void TProof::AddChain(TChain *chain)
{
   fChains->Add(chain);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove chain from data set

void TProof::RemoveChain(TChain *chain)
{
   fChains->Remove(chain);
}

////////////////////////////////////////////////////////////////////////////////
/// Ask for remote logs in the range [start, end]. If start == -1 all the
/// messages not yet received are sent back.

void TProof::GetLog(Int_t start, Int_t end)
{
   if (!IsValid() || TestBit(TProof::kIsMaster)) return;

   TMessage msg(kPROOF_LOGFILE);

   msg << start << end;

   Broadcast(msg, kActive);
   Collect(kActive, fCollectTimeout);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill a TMacro with the log lines since the last reading (fLogFileR)
/// Return (TMacro *)0 if no line was logged.
/// The returned TMacro must be deleted by the caller.

TMacro *TProof::GetLastLog()
{
   TMacro *maclog = 0;

   // Save present offset
   off_t nowlog = lseek(fileno(fLogFileR), (off_t) 0, SEEK_CUR);
   if (nowlog < 0) {
      SysError("GetLastLog",
               "problem lseeking log file to current position (errno: %d)", TSystem::GetErrno());
      return maclog;
   }

   // Get extremes
   off_t startlog = nowlog;
   off_t endlog = lseek(fileno(fLogFileR), (off_t) 0, SEEK_END);
   if (endlog < 0) {
      SysError("GetLastLog",
               "problem lseeking log file to end position (errno: %d)", TSystem::GetErrno());
      return maclog;
   }

   // Perhaps nothing to log
   UInt_t tolog = (UInt_t)(endlog - startlog);
   if (tolog <= 0) return maclog;

   // Set starting point
   if (lseek(fileno(fLogFileR), startlog, SEEK_SET) < 0) {
      SysError("GetLastLog",
               "problem lseeking log file to start position (errno: %d)", TSystem::GetErrno());
      return maclog;
   }

   // Create the output object
   maclog = new TMacro;

   // Now we go
   char line[2048];
   Int_t wanted = (tolog > sizeof(line)) ? sizeof(line) : tolog;
   while (fgets(line, wanted, fLogFileR)) {
      Int_t r = strlen(line);
      if (r > 0) {
         if (line[r-1] == '\n') line[r-1] = '\0';
         maclog->AddLine(line);
      } else {
         // Done
         break;
      }
      tolog -= r;
      wanted = (tolog > sizeof(line)) ? sizeof(line) : tolog;
   }

   // Restore original pointer
   if (lseek(fileno(fLogFileR), nowlog, SEEK_SET) < 0) {
      Warning("GetLastLog",
              "problem lseeking log file to original position (errno: %d)", TSystem::GetErrno());
   }

   // Done
   return maclog;
}

////////////////////////////////////////////////////////////////////////////////
/// Display log of query pq into the log window frame

void TProof::PutLog(TQueryResult *pq)
{
   if (!pq) return;

   TList *lines = pq->GetLogFile()->GetListOfLines();
   if (lines) {
      TIter nxl(lines);
      TObjString *l = 0;
      while ((l = (TObjString *)nxl()))
         EmitVA("LogMessage(const char*,Bool_t)", 2, l->GetName(), kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Display on screen the content of the temporary log file for query
/// in reference

void TProof::ShowLog(const char *queryref)
{
   // Make sure we have all info (GetListOfQueries retrieves the
   // head info only)
   Retrieve(queryref);

   if (fPlayer) {
      if (queryref) {
         if (fPlayer->GetListOfResults()) {
            TIter nxq(fPlayer->GetListOfResults());
            TQueryResult *qr = 0;
            while ((qr = (TQueryResult *) nxq()))
               if (strstr(queryref, qr->GetTitle()) &&
                   strstr(queryref, qr->GetName()))
                  break;
            if (qr) {
               PutLog(qr);
               return;
            }

         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Display on screen the content of the temporary log file.
/// If qry == -2 show messages from the last (current) query.
/// If qry == -1 all the messages not yet displayed are shown (default).
/// If qry == 0, all the messages in the file are shown.
/// If qry  > 0, only the messages related to query 'qry' are shown.
/// For qry != -1 the original file offset is restored at the end

void TProof::ShowLog(Int_t qry)
{
   // Save present offset
   off_t nowlog = lseek(fileno(fLogFileR), (off_t) 0, SEEK_CUR);
   if (nowlog < 0) {
      SysError("ShowLog", "problem lseeking log file (errno: %d)", TSystem::GetErrno());
      return;
   }

   // Get extremes
   off_t startlog = nowlog;
   off_t endlog = lseek(fileno(fLogFileR), (off_t) 0, SEEK_END);
   if (endlog < 0) {
      SysError("ShowLog", "problem lseeking log file (errno: %d)", TSystem::GetErrno());
      return;
   }

   lseek(fileno(fLogFileR), nowlog, SEEK_SET);
   if (qry == 0) {
      startlog = 0;
      lseek(fileno(fLogFileR), (off_t) 0, SEEK_SET);
   } else if (qry != -1) {

      TQueryResult *pq = 0;
      if (qry == -2) {
         // Pickup the last one
         pq = (GetQueryResults()) ? ((TQueryResult *)(GetQueryResults()->Last())) : 0;
         if (!pq) {
            GetListOfQueries();
            if (fQueries)
               pq = (TQueryResult *)(fQueries->Last());
         }
      } else if (qry > 0) {
         TList *queries = GetQueryResults();
         if (queries) {
            TIter nxq(queries);
            while ((pq = (TQueryResult *)nxq()))
               if (qry == pq->GetSeqNum())
                  break;
         }
         if (!pq) {
            queries = GetListOfQueries();
            TIter nxq(queries);
            while ((pq = (TQueryResult *)nxq()))
               if (qry == pq->GetSeqNum())
                  break;
         }
      }
      if (pq) {
         PutLog(pq);
         return;
      } else {
         if (gDebug > 0)
            Info("ShowLog","query %d not found in list", qry);
         qry = -1;
      }
   }

   // Number of bytes to log
   UInt_t tolog = (UInt_t)(endlog - startlog);

   // Perhaps nothing
   if (tolog <= 0) {
      // Set starting point
      lseek(fileno(fLogFileR), startlog, SEEK_SET);
   }

   // Now we go
   Int_t np = 0;
   char line[2048];
   Int_t wanted = (tolog > sizeof(line)) ? sizeof(line) : tolog;
   while (fgets(line, wanted, fLogFileR)) {

      Int_t r = strlen(line);
      if (!SendingLogToWindow()) {
         if (line[r-1] != '\n') line[r-1] = '\n';
         if (r > 0) {
            char *p = line;
            while (r) {
               Int_t w = write(fileno(stdout), p, r);
               if (w < 0) {
                  SysError("ShowLog", "error writing to stdout");
                  break;
               }
               r -= w;
               p += w;
            }
         }
         tolog -= strlen(line);
         np++;

         // Ask if more is wanted
         if (!(np%10)) {
            const char *opt = Getline("More (y/n)? [y]");
            if (opt[0] == 'n')
               break;
         }

         // We may be over
         if (tolog <= 0)
            break;

         // Update wanted bytes
         wanted = (tolog > sizeof(line)) ? sizeof(line) : tolog;
      } else {
         // Log to window
         if (line[r-1] == '\n') line[r-1] = 0;
         LogMessage(line, kFALSE);
      }
   }
   if (!SendingLogToWindow()) {
      // Avoid screwing up the prompt
      if (write(fileno(stdout), "\n", 1) != 1)
         SysError("ShowLog", "error writing to stdout");
   }

   // Restore original pointer
   if (qry > -1)
      lseek(fileno(fLogFileR), nowlog, SEEK_SET);
}

////////////////////////////////////////////////////////////////////////////////
/// Set session with 'id' the default one. If 'id' is not found in the list,
/// the current session is set as default

void TProof::cd(Int_t id)
{
   if (GetManager()) {
      TProofDesc *d = GetManager()->GetProofDesc(id);
      if (d) {
         if (d->GetProof()) {
            gProof = d->GetProof();
            return;
         }
      }

      // Id not found or undefined: set as default this session
      gProof = this;
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Detach this instance to its proofserv.
/// If opt is 'S' or 's' the remote server is shutdown

void TProof::Detach(Option_t *opt)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) return;

   // Get worker and socket instances
   TSlave *sl = (TSlave *) fActiveSlaves->First();
   TSocket *s = 0;
   if (!sl || !(sl->IsValid()) || !(s = sl->GetSocket())) {
      Error("Detach","corrupted worker instance: wrk:%p, sock:%p", sl, s);
      return;
   }

   Bool_t shutdown = (strchr(opt,'s') || strchr(opt,'S')) ? kTRUE : kFALSE;

   // If processing, try to stop processing first
   if (shutdown && !IsIdle()) {
      // Remove pending requests
      Remove("cleanupqueue");
      // Do not wait for ever, but al least 20 seconds
      Long_t timeout = gEnv->GetValue("Proof.ShutdownTimeout", 60);
      timeout = (timeout > 20) ? timeout : 20;
      // Send stop signal
      StopProcess(kFALSE, (Long_t) (timeout / 2));
      // Receive results
      Collect(kActive, timeout);
   }

   // Avoid spurious messages: deactivate new inputs ...
   DeActivateAsyncInput();

   // ... and discard existing ones
   sl->FlushSocket();

   // Close session (we always close the connection)
   Close(opt);

   // Close the progress dialog, if any
   if (fProgressDialogStarted)
      CloseProgressDialog();

   // Update info in the table of our manager, if any
   if (GetManager() && GetManager()->QuerySessions("L")) {
      TIter nxd(GetManager()->QuerySessions("L"));
      TProofDesc *d = 0;
      while ((d = (TProofDesc *)nxd())) {
         if (d->GetProof() == this) {
            d->SetProof(0);
            GetManager()->QuerySessions("L")->Remove(d);
            break;
         }
      }
   }

   // Invalidate this instance
   fValid = kFALSE;

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Set an alias for this session. If reconnection is supported, the alias
/// will be communicated to the remote coordinator so that it can be recovered
/// when reconnecting

void TProof::SetAlias(const char *alias)
{
   // Set it locally
   TNamed::SetTitle(alias);
   if (TestBit(TProof::kIsMaster))
      // Set the name at the same value
      TNamed::SetName(alias);

   // Nothing to do if not in contact with coordinator
   if (!IsValid()) return;

   if (!IsProofd() && TestBit(TProof::kIsClient)) {
      TSlave *sl = (TSlave *) fActiveSlaves->First();
      if (sl)
         sl->SetAlias(alias);
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// *** This function is deprecated and will disappear in future versions ***
/// *** It is just a wrapper around TFile::Cp.
/// *** Please use TProofMgr::UploadFiles.
///
/// Upload a set of files and save the list of files by name dataSetName.
/// The 'files' argument is a list of TFileInfo objects describing the files
/// as first url.
/// The mask 'opt' is a combination of EUploadOpt:
///   kAppend             (0x1)   if set true files will be appended to
///                               the dataset existing by given name
///   kOverwriteDataSet   (0x2)   if dataset with given name exited it
///                               would be overwritten
///   kNoOverwriteDataSet (0x4)   do not overwirte if the dataset exists
///   kOverwriteAllFiles  (0x8)   overwrite all files that may exist
///   kOverwriteNoFiles   (0x10)  overwrite none
///   kAskUser            (0x0)   ask user before overwriteng dataset/files
/// The default value is kAskUser.
/// The user will be asked to confirm overwriting dataset or files unless
/// specified opt provides the answer!
/// If kOverwriteNoFiles is set, then a pointer to TList must be passed as
/// skippedFiles argument. The function will add to this list TFileInfo
/// objects describing all files that existed on the cluster and were
/// not uploaded.
///
/// Communication Summary
/// Client                             Master
///    |------------>DataSetName----------->|
///    |<-------kMESS_OK/kMESS_NOTOK<-------| (Name OK/file exist)
/// (*)|-------> call RegisterDataSet ------->|
/// (*) - optional

Int_t TProof::UploadDataSet(const char *, TList *, const char *, Int_t, TList *)
{
   Printf(" *** WARNING: this function is obsolete: it has been replaced by TProofMgr::UploadFiles ***");

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// *** This function is deprecated and will disappear in future versions ***
/// *** It is just a wrapper around TFile::Cp.
/// *** Please use TProofMgr::UploadFiles.
///
/// Upload a set of files and save the list of files by name dataSetName.
/// The mask 'opt' is a combination of EUploadOpt:
///   kAppend             (0x1)   if set true files will be appended to
///                               the dataset existing by given name
///   kOverwriteDataSet   (0x2)   if dataset with given name exited it
///                               would be overwritten
///   kNoOverwriteDataSet (0x4)   do not overwirte if the dataset exists
///   kOverwriteAllFiles  (0x8)   overwrite all files that may exist
///   kOverwriteNoFiles   (0x10)  overwrite none
///   kAskUser            (0x0)   ask user before overwriteng dataset/files
/// The default value is kAskUser.
/// The user will be asked to confirm overwriting dataset or files unless
/// specified opt provides the answer!
/// If kOverwriteNoFiles is set, then a pointer to TList must be passed as
/// skippedFiles argument. The function will add to this list TFileInfo
/// objects describing all files that existed on the cluster and were
/// not uploaded.
///

Int_t TProof::UploadDataSet(const char *, const char *, const char *, Int_t, TList *)
{
   Printf(" *** WARNING: this function is obsolete: it has been replaced by TProofMgr::UploadFiles ***");

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// *** This function is deprecated and will disappear in future versions ***
/// *** It is just a wrapper around TFile::Cp.
/// *** Please use TProofMgr::UploadFiles.
///
/// Upload files listed in "file" to PROOF cluster.
/// Where file = name of file containing list of files and
/// dataset = dataset name and opt is a combination of EUploadOpt bits.
/// Each file description (line) can include wildcards.
/// Check TFileInfo compatibility

Int_t TProof::UploadDataSetFromFile(const char *, const char *, const char *, Int_t, TList *)
{
   Printf(" *** WARNING: this function is obsolete: it has been replaced by TProofMgr::UploadFiles ***");

    // Done
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Register the 'dataSet' on the cluster under the current
/// user, group and the given 'dataSetName'.
/// If a dataset with the same name already exists the action fails unless 'opts'
/// contains 'O', in which case the old dataset is overwritten, or contains 'U',
/// in which case 'newDataSet' is added to the existing dataset (duplications are
/// ignored, if any).
/// If 'opts' contains 'V' the dataset files are also verified (if the dataset manager
/// is configured to allow so). By default the dataset is not verified.
/// If 'opts' contains 'T' the in the dataset object (status bits, meta,...)
/// is trusted, i.e. not reset (if the dataset manager is configured to allow so).
/// If 'opts' contains 'S' validation would be run serially (meaningful only if
/// validation is required).
/// Returns kTRUE on success.

Bool_t TProof::RegisterDataSet(const char *dataSetName,
                               TFileCollection *dataSet, const char *optStr)
{
   // Check TFileInfo compatibility
   if (fProtocol < 17) {
      Info("RegisterDataSet",
           "functionality not available: the server does not have dataset support");
      return kFALSE;
   }

   if (!dataSetName || strlen(dataSetName) <= 0) {
      Info("RegisterDataSet", "specifying a dataset name is mandatory");
      return kFALSE;
   }

   Bool_t parallelverify = kFALSE;
   TString sopt(optStr);
   if (sopt.Contains("V") && fProtocol >= 34 && !sopt.Contains("S")) {
      // We do verification in parallel later on; just register for now
      parallelverify = kTRUE;
      sopt.ReplaceAll("V", "");
   }
   // This would screw up things remotely, make sure is not there
   sopt.ReplaceAll("S", "");

   TMessage mess(kPROOF_DATASETS);
   mess << Int_t(kRegisterDataSet);
   mess << TString(dataSetName);
   mess << sopt;
   mess.WriteObject(dataSet);
   Broadcast(mess);

   Bool_t result = kTRUE;
   Collect();
   if (fStatus != 0) {
      Error("RegisterDataSet", "dataset was not saved");
      result = kFALSE;
      return result;
   }

   // If old server or not verifying in parallel we are done
   if (!parallelverify) return result;

   // If we are here it means that we will verify in parallel
   sopt += "V";
   if (VerifyDataSet(dataSetName, sopt) < 0){
      Error("RegisterDataSet", "problems verifying dataset '%s'", dataSetName);
      return kFALSE;
   }

   // We are done
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set/Change the name of the default tree. The tree name may contain
/// subdir specification in the form "subdir/name".
/// Returns 0 on success, -1 otherwise.

Int_t TProof::SetDataSetTreeName(const char *dataset, const char *treename)
{
   // Check TFileInfo compatibility
   if (fProtocol < 23) {
      Info("SetDataSetTreeName", "functionality not supported by the server");
      return -1;
   }

   if (!dataset || strlen(dataset) <= 0) {
      Info("SetDataSetTreeName", "specifying a dataset name is mandatory");
      return -1;
   }

   if (!treename || strlen(treename) <= 0) {
      Info("SetDataSetTreeName", "specifying a tree name is mandatory");
      return -1;
   }

   TUri uri(dataset);
   TString fragment(treename);
   if (!fragment.BeginsWith("/")) fragment.Insert(0, "/");
   uri.SetFragment(fragment);

   TMessage mess(kPROOF_DATASETS);
   mess << Int_t(kSetDefaultTreeName);
   mess << uri.GetUri();
   Broadcast(mess);

   Collect();
   if (fStatus != 0) {
      Error("SetDataSetTreeName", "some error occured: default tree name not changed");
      return -1;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Lists all datasets that match given uri.
/// The 'optStr' can contain a comma-separated list of servers for which the
/// information is wanted. If ':lite:' (case insensitive) is specified in 'optStr'
/// only the global information in the TFileCollection is retrieved; useful to only
/// get the list of available datasets.

TMap *TProof::GetDataSets(const char *uri, const char *optStr)
{
   if (fProtocol < 15) {
      Info("GetDataSets",
           "functionality not available: the server does not have dataset support");
      return 0;
   }
   if (fProtocol < 31 && strstr(optStr, ":lite:"))
      Warning("GetDataSets", "'lite' option not supported by the server");

   TMessage mess(kPROOF_DATASETS);
   mess << Int_t(kGetDataSets);
   mess << TString(uri ? uri : "");
   mess << TString(optStr ? optStr : "");
   Broadcast(mess);
   Collect(kActive, fCollectTimeout);

   TMap *dataSetMap = 0;
   if (fStatus != 0) {
      Error("GetDataSets", "error receiving datasets information");
   } else {
      // Look in the list
      TMessage *retMess = (TMessage *) fRecvMessages->First();
      if (retMess && retMess->What() == kMESS_OK) {
         if (!(dataSetMap = (TMap *)(retMess->ReadObject(TMap::Class()))))
            Error("GetDataSets", "error receiving datasets");
      } else
         Error("GetDataSets", "message not found or wrong type (%p)", retMess);
   }

   return dataSetMap;
}

////////////////////////////////////////////////////////////////////////////////
/// Shows datasets in locations that match the uri.
/// By default shows the user's datasets and global ones

void TProof::ShowDataSets(const char *uri, const char* optStr)
{
   if (fProtocol < 15) {
      Info("ShowDataSets",
           "functionality not available: the server does not have dataset support");
      return;
   }

   TMessage mess(kPROOF_DATASETS);
   mess << Int_t(kShowDataSets);
   mess << TString(uri ? uri : "");
   mess << TString(optStr ? optStr : "");
   Broadcast(mess);

   Collect(kActive, fCollectTimeout);
   if (fStatus != 0)
      Error("ShowDataSets", "error receiving datasets information");
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if 'dataset' exists, kFALSE otherwise

Bool_t TProof::ExistsDataSet(const char *dataset)
{
   if (fProtocol < 15) {
      Info("ExistsDataSet", "functionality not available: the server has an"
                            " incompatible version of TFileInfo");
      return kFALSE;
   }

   if (!dataset || strlen(dataset) <= 0) {
      Error("ExistsDataSet", "dataset name missing");
      return kFALSE;
   }

   TMessage msg(kPROOF_DATASETS);
   msg << Int_t(kCheckDataSetName) << TString(dataset);
   Broadcast(msg);
   Collect(kActive, fCollectTimeout);
   if (fStatus == -1) {
      // The dataset exists
      return kTRUE;
   }
   // The dataset does not exists
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the content of the dataset cache, if any (matching 'dataset', if defined).

void TProof::ClearDataSetCache(const char *dataset)
{
   if (fProtocol < 28) {
      Info("ClearDataSetCache", "functionality not available on server");
      return;
   }

   TMessage msg(kPROOF_DATASETS);
   msg << Int_t(kCache) << TString(dataset) << TString("clear");
   Broadcast(msg);
   Collect(kActive, fCollectTimeout);
   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Display the content of the dataset cache, if any (matching 'dataset', if defined).

void TProof::ShowDataSetCache(const char *dataset)
{
   if (fProtocol < 28) {
      Info("ShowDataSetCache", "functionality not available on server");
      return;
   }

   TMessage msg(kPROOF_DATASETS);
   msg << Int_t(kCache) << TString(dataset) << TString("show");
   Broadcast(msg);
   Collect(kActive, fCollectTimeout);
   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a list of TFileInfo objects describing the files of the specified
/// dataset.
/// To get the short version (containing only the global meta information)
/// specify optStr = "S:" or optStr = "short:".
/// To get the sub-dataset of files located on a given server(s) specify
/// the list of servers (comma-separated) in the 'optStr' field.

TFileCollection *TProof::GetDataSet(const char *uri, const char *optStr)
{
   if (fProtocol < 15) {
      Info("GetDataSet", "functionality not available: the server has an"
                         " incompatible version of TFileInfo");
      return 0;
   }

   if (!uri || strlen(uri) <= 0) {
      Info("GetDataSet", "specifying a dataset name is mandatory");
      return 0;
   }

   TMessage nameMess(kPROOF_DATASETS);
   nameMess << Int_t(kGetDataSet);
   nameMess << TString(uri);
   nameMess << TString(optStr ? optStr: "");
   if (Broadcast(nameMess) < 0)
      Error("GetDataSet", "sending request failed");

   Collect(kActive, fCollectTimeout);
   TFileCollection *fileList = 0;
   if (fStatus != 0) {
      Error("GetDataSet", "error receiving datasets information");
   } else {
      // Look in the list
      TMessage *retMess = (TMessage *) fRecvMessages->First();
      if (retMess && retMess->What() == kMESS_OK) {
         if (!(fileList = (TFileCollection*)(retMess->ReadObject(TFileCollection::Class()))))
            Error("GetDataSet", "error reading list of files");
      } else
         Error("GetDataSet", "message not found or wrong type (%p)", retMess);
   }

   return fileList;
}

////////////////////////////////////////////////////////////////////////////////
/// display meta-info for given dataset usi

void TProof::ShowDataSet(const char *uri, const char* opt)
{
   TFileCollection *fileList = 0;
   if ((fileList = GetDataSet(uri))) {
      fileList->Print(opt);
      delete fileList;
   } else
      Warning("ShowDataSet","no such dataset: %s", uri);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the specified dataset from the PROOF cluster.
/// Files are not deleted.

Int_t TProof::RemoveDataSet(const char *uri, const char* optStr)
{
   TMessage nameMess(kPROOF_DATASETS);
   nameMess << Int_t(kRemoveDataSet);
   nameMess << TString(uri?uri:"");
   nameMess << TString(optStr?optStr:"");
   if (Broadcast(nameMess) < 0)
      Error("RemoveDataSet", "sending request failed");
   Collect(kActive, fCollectTimeout);

   if (fStatus != 0)
      return -1;
   else
      return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Find datasets, returns in a TList all found datasets.

TList* TProof::FindDataSets(const char* /*searchString*/, const char* /*optStr*/)
{
   Error ("FindDataSets", "not yet implemented");
   return (TList *) 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Allows users to request staging of a particular dataset. Requests are
/// saved in a special dataset repository and must be honored by the endpoint.

Bool_t TProof::RequestStagingDataSet(const char *dataset)
{
   if (fProtocol < 35) {
      Error("RequestStagingDataSet",
         "functionality not supported by the server");
      return kFALSE;
   }

   TMessage mess(kPROOF_DATASETS);
   mess << Int_t(kRequestStaging);
   mess << TString(dataset);
   Broadcast(mess);

   Collect();
   if (fStatus != 0) {
      Error("RequestStagingDataSet", "staging request was unsuccessful");
      return kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Cancels a dataset staging request. Returns kTRUE on success, kFALSE on
/// failure. Dataset not found equals to a failure.

Bool_t TProof::CancelStagingDataSet(const char *dataset)
{
   if (fProtocol < 36) {
      Error("CancelStagingDataSet",
         "functionality not supported by the server");
      return kFALSE;
   }

   TMessage mess(kPROOF_DATASETS);
   mess << Int_t(kCancelStaging);
   mess << TString(dataset);
   Broadcast(mess);

   Collect();
   if (fStatus != 0) {
      Error("CancelStagingDataSet", "cancel staging request was unsuccessful");
      return kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Obtains a TFileCollection showing the staging status of the specified
/// dataset. A valid dataset manager and dataset staging requests repository
/// must be present on the endpoint.

TFileCollection *TProof::GetStagingStatusDataSet(const char *dataset)
{
   if (fProtocol < 35) {
      Error("GetStagingStatusDataSet",
         "functionality not supported by the server");
      return NULL;
   }

   TMessage nameMess(kPROOF_DATASETS);
   nameMess << Int_t(kStagingStatus);
   nameMess << TString(dataset);
   if (Broadcast(nameMess) < 0) {
      Error("GetStagingStatusDataSet", "sending request failed");
      return NULL;
   }

   Collect(kActive, fCollectTimeout);
   TFileCollection *fc = NULL;

   if (fStatus < 0) {
      Error("GetStagingStatusDataSet", "problem processing the request");
   }
   else if (fStatus == 0) {
      TMessage *retMess = (TMessage *)fRecvMessages->First();
      if (retMess && (retMess->What() == kMESS_OK)) {
         fc = (TFileCollection *)(
           retMess->ReadObject(TFileCollection::Class()) );
         if (!fc)
            Error("GetStagingStatusDataSet", "error reading list of files");
      }
      else {
         Error("GetStagingStatusDataSet",
            "response message not found or wrong type (%p)", retMess);
      }
   }
   //else {}

   return fc;
}

////////////////////////////////////////////////////////////////////////////////
/// Like GetStagingStatusDataSet, but displays results immediately.

void TProof::ShowStagingStatusDataSet(const char *dataset, const char *opt)
{
   TFileCollection *fc = GetStagingStatusDataSet(dataset);
   if (fc) {
      fc->Print(opt);
      delete fc;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Verify if all files in the specified dataset are available.
/// Print a list and return the number of missing files.
/// Returns -1 in case of error.

Int_t TProof::VerifyDataSet(const char *uri, const char *optStr)
{
   if (fProtocol < 15) {
      Info("VerifyDataSet", "functionality not available: the server has an"
                            " incompatible version of TFileInfo");
      return -1;
   }

   // Sanity check
   if (!uri || (uri && strlen(uri) <= 0)) {
      Error("VerifyDataSet", "dataset name is is mandatory");
      return -1;
   }

   Int_t nmissingfiles = 0;

   TString sopt(optStr);
   if (fProtocol < 34 || sopt.Contains("S")) {
      sopt.ReplaceAll("S", "");
      Info("VerifyDataSet", "Master-only verification");
      TMessage nameMess(kPROOF_DATASETS);
      nameMess << Int_t(kVerifyDataSet);
      nameMess << TString(uri);
      nameMess << sopt;
      Broadcast(nameMess);

      Collect(kActive, fCollectTimeout);

      if (fStatus < 0) {
         Info("VerifyDataSet", "no such dataset %s", uri);
         return  -1;
      } else
         nmissingfiles = fStatus;
      return nmissingfiles;
   }

   // Request for parallel verification: can only be done if we have workers
   if (!IsParallel() && !fDynamicStartup) {
      Error("VerifyDataSet", "PROOF is in sequential mode (no workers): cannot do parallel verification.");
      Error("VerifyDataSet", "Either start PROOF with some workers or force sequential adding 'S' as option.");
      return -1;
   }

   // Do parallel verification
   return VerifyDataSetParallel(uri, optStr);
}

////////////////////////////////////////////////////////////////////////////////
/// Internal function for parallel dataset verification used TProof::VerifyDataSet and
/// TProofLite::VerifyDataSet

Int_t TProof::VerifyDataSetParallel(const char *uri, const char *optStr)
{
   Int_t nmissingfiles = 0;

   // Let PROOF master prepare node-files map
   SetParameter("PROOF_FilesToProcess", Form("dataset:%s", uri));

   // Use TPacketizerFile
   TString oldpack;
   if (TProof::GetParameter(GetInputList(), "PROOF_Packetizer", oldpack) != 0) oldpack = "";
   SetParameter("PROOF_Packetizer", "TPacketizerFile");

   // Add dataset name
   SetParameter("PROOF_VerifyDataSet", uri);
   // Add options
   SetParameter("PROOF_VerifyDataSetOption", optStr);
   SetParameter("PROOF_SavePartialResults", (Int_t)0);
   Int_t oldifiip = -1;
   if (TProof::GetParameter(GetInputList(), "PROOF_IncludeFileInfoInPacket", oldifiip) != 0) oldifiip = -1;
   SetParameter("PROOF_IncludeFileInfoInPacket", (Int_t)1);

   // TO DO : figure out mss and stageoption
   const char* mss="";
   SetParameter("PROOF_MSS", mss);
   const char* stageoption="";
   SetParameter("PROOF_StageOption", stageoption);

   // Process verification in parallel
   Process("TSelVerifyDataSet", (Long64_t) 1);

   // Restore packetizer
   if (!oldpack.IsNull())
      SetParameter("PROOF_Packetizer", oldpack);
   else
      DeleteParameters("PROOF_Packetizer");

   // Delete or restore parameters
   DeleteParameters("PROOF_FilesToProcess");
   DeleteParameters("PROOF_VerifyDataSet");
   DeleteParameters("PROOF_VerifyDataSetOption");
   DeleteParameters("PROOF_MSS");
   DeleteParameters("PROOF_StageOption");
   if (oldifiip > -1) {
      SetParameter("PROOF_IncludeFileInfoInPacket", oldifiip);
   } else {
      DeleteParameters("PROOF_IncludeFileInfoInPacket");
   }
   DeleteParameters("PROOF_SavePartialResults");

   // Merge outputs
   Int_t nopened = 0;
   Int_t ntouched = 0;
   Bool_t changed_ds = kFALSE;

   TIter nxtout(GetOutputList());
   TObject* obj;
   TList *lfiindout = new TList;
   while ((obj = nxtout())) {
      TList *l = dynamic_cast<TList *>(obj);
      if (l && TString(l->GetName()).BeginsWith("PROOF_ListFileInfos_")) {
         TIter nxt(l);
         TFileInfo *fiindout = 0;
         while ((fiindout = (TFileInfo*) nxt())) {
            lfiindout->Add(fiindout);
         }
      }
      // Add up number of disppeared files
      TParameter<Int_t>* pdisappeared = dynamic_cast<TParameter<Int_t>*>(obj);
      if ( pdisappeared && TString(pdisappeared->GetName()).BeginsWith("PROOF_NoFilesDisppeared_")) {
         nmissingfiles += pdisappeared->GetVal();
      }
      TParameter<Int_t>* pnopened = dynamic_cast<TParameter<Int_t>*>(obj);
      if (pnopened && TString(pnopened->GetName()).BeginsWith("PROOF_NoFilesOpened_")) {
         nopened += pnopened->GetVal();
      }
      TParameter<Int_t>* pntouched = dynamic_cast<TParameter<Int_t>*>(obj);
      if (pntouched && TString(pntouched->GetName()).BeginsWith("PROOF_NoFilesTouched_")) {
         ntouched += pntouched->GetVal();
      }
      TParameter<Bool_t>* pchanged_ds = dynamic_cast<TParameter<Bool_t>*>(obj);
      if (pchanged_ds && TString(pchanged_ds->GetName()).BeginsWith("PROOF_DataSetChanged_")) {
         if (pchanged_ds->GetVal() == kTRUE) changed_ds = kTRUE;
      }
   }

   Info("VerifyDataSetParallel", "%s: changed? %d (# files opened = %d, # files touched = %d,"
                                 " # missing files = %d)",
                                 uri, changed_ds, nopened, ntouched, nmissingfiles);
   // Done
   return nmissingfiles;
}

////////////////////////////////////////////////////////////////////////////////
/// returns a map of the quotas of all groups

TMap *TProof::GetDataSetQuota(const char* optStr)
{
   if (IsLite()) {
      Info("UploadDataSet", "Lite-session: functionality not implemented");
      return (TMap *)0;
   }

   TMessage mess(kPROOF_DATASETS);
   mess << Int_t(kGetQuota);
   mess << TString(optStr?optStr:"");
   Broadcast(mess);

   Collect(kActive, fCollectTimeout);
   TMap *groupQuotaMap = 0;
   if (fStatus < 0) {
      Info("GetDataSetQuota", "could not receive quota");
   } else {
      // Look in the list
      TMessage *retMess = (TMessage *) fRecvMessages->First();
      if (retMess && retMess->What() == kMESS_OK) {
         if (!(groupQuotaMap = (TMap*)(retMess->ReadObject(TMap::Class()))))
            Error("GetDataSetQuota", "error getting quotas");
      } else
         Error("GetDataSetQuota", "message not found or wrong type (%p)", retMess);
   }

   return groupQuotaMap;
}

////////////////////////////////////////////////////////////////////////////////
/// shows the quota and usage of all groups
/// if opt contains "U" shows also distribution of usage on user-level

void TProof::ShowDataSetQuota(Option_t* opt)
{
   if (fProtocol < 15) {
      Info("ShowDataSetQuota",
           "functionality not available: the server does not have dataset support");
      return;
   }

   if (IsLite()) {
      Info("UploadDataSet", "Lite-session: functionality not implemented");
      return;
   }

   TMessage mess(kPROOF_DATASETS);
   mess << Int_t(kShowQuota);
   mess << TString(opt?opt:"");
   Broadcast(mess);

   Collect();
   if (fStatus != 0)
      Error("ShowDataSetQuota", "error receiving quota information");
}

////////////////////////////////////////////////////////////////////////////////
/// If in active in a monitor set ready state

void TProof::InterruptCurrentMonitor()
{
   if (fCurrentMonitor)
      fCurrentMonitor->Interrupt();
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that the worker identified by the ordinal number 'ord' is
/// in the active list. The request will be forwarded to the master
/// in direct contact with the worker. If needed, this master will move
/// the worker from the inactive to the active list and rebuild the list
/// of unique workers.
/// Use ord = "*" to activate all inactive workers.
/// The string 'ord' can also be a comma-separated list of ordinal numbers the
/// status of which will be modified at once.
/// Return <0 if something went wrong (-2 if at least one worker was not found)
/// or the number of workers with status change (on master; 0 on client).

Int_t TProof::ActivateWorker(const char *ord, Bool_t save)
{
   return ModifyWorkerLists(ord, kTRUE, save);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the worker identified by the ordinal number 'ord' from the
/// the active list. The request will be forwarded to the master
/// in direct contact with the worker. If needed, this master will move
/// the worker from the active to the inactive list and rebuild the list
/// of unique workers.
/// Use ord = "*" to deactivate all active workers.
/// The string 'ord' can also be a comma-separated list of ordinal numbers the
/// status of which will be modified at once.
/// Return <0 if something went wrong (-2 if at least one worker was not found)
/// or the number of workers with status change (on master; 0 on client).

Int_t TProof::DeactivateWorker(const char *ord, Bool_t save)
{
   return ModifyWorkerLists(ord, kFALSE, save);
}

////////////////////////////////////////////////////////////////////////////////
/// Modify the worker active/inactive list by making the worker identified by
/// the ordinal number 'ord' active (add == TRUE) or inactive (add == FALSE).
/// The string 'ord' can also be a comma-separated list of ordinal numbers the
/// status of which will be modified at once.
/// If needed, the request will be forwarded to the master in direct contact
/// with the worker. The end-master will move the worker from one list to the
/// other active and rebuild the list of unique active workers.
/// Use ord = "*" to deactivate all active workers.
/// If save is TRUE the current active list is saved before any modification is
/// done; re-running with ord = "restore" restores the saved list
/// Return <0 if something went wrong (-2 if at least one worker was not found)
/// or the number of workers with status change (on master; 0 on client).

Int_t TProof::ModifyWorkerLists(const char *ord, Bool_t add, Bool_t save)
{
   // Make sure the input make sense
   if (!ord || strlen(ord) <= 0) {
      Info("ModifyWorkerLists",
           "an ordinal number - e.g. \"0.4\" or \"*\" for all - is required as input");
      return -1;
   }
   if (gDebug > 0)
      Info("ModifyWorkerLists", "ord: '%s' (add: %d, save: %d)", ord, add, save);

   Int_t nwc = 0;
   Bool_t restoring = !strcmp(ord, "restore") ? kTRUE : kFALSE;
   if (IsEndMaster()) {
      if (restoring) {
         // We are asked to restore the previous settings
         nwc = RestoreActiveList();
      } else {
         if (save) SaveActiveList();
      }
   }

   Bool_t allord = strcmp(ord, "*") ? kFALSE : kTRUE;

   // Check if this is for us
   if (TestBit(TProof::kIsMaster) && gProofServ) {
      if (!allord &&
         strncmp(ord, gProofServ->GetOrdinal(), strlen(gProofServ->GetOrdinal())))
         return 0;
   }

   Bool_t fw = kTRUE;    // Whether to forward one step down
   Bool_t rs = kFALSE;   // Whether to rescan for unique workers

   // Appropriate list pointing
   TList *in = (add) ? fInactiveSlaves : fActiveSlaves;
   TList *out = (add) ? fActiveSlaves : fInactiveSlaves;

   if (IsEndMaster() && !restoring) {
      // Create the hash list of ordinal numbers
      THashList *ords = 0;
      if (!allord) {
         ords = new THashList();
         const char *masterord = (gProofServ) ? gProofServ->GetOrdinal() : "0";
         TString oo(ord), o;
         Int_t from = 0;
         while(oo.Tokenize(o, from, ","))
            if (o.BeginsWith(masterord)) ords->Add(new TObjString(o));
      }
      // We do not need to send forward
      fw = kFALSE;
      // Look for the worker in the initial list
      TObject *os = 0;
      TSlave *wrk = 0;
      if (in->GetSize() > 0) {
         TIter nxw(in);
         while ((wrk = (TSlave *) nxw())) {
            os = 0;
            if (allord || (ords && (os = ords->FindObject(wrk->GetOrdinal())))) {
               // Add it to the final list
               if (!out->FindObject(wrk)) {
                  out->Add(wrk);
                  if (add)
                     fActiveMonitor->Add(wrk->GetSocket());
               }
               // Remove it from the initial list
               in->Remove(wrk);
               if (!add) {
                  fActiveMonitor->Remove(wrk->GetSocket());
                  wrk->SetStatus(TSlave::kInactive);
               } else
                  wrk->SetStatus(TSlave::kActive);
               // Count
               nwc++;
               // Nothing to forward (ord is unique)
               fw = kFALSE;
               // Rescan for unique workers (active list modified)
               rs = kTRUE;
               // We may be done, if not option 'all'
               if (!allord && ords) {
                  if (os) ords->Remove(os);
                  if (ords->GetSize() == 0) break;
                  SafeDelete(os);
               }
            }
         }
      }
      // If some worker not found, notify it if at the end
      if (!fw && ords && ords->GetSize() > 0) {
         TString oo;
         TIter nxo(ords);
         while ((os = nxo())) {
            TIter nxw(out);
            while ((wrk = (TSlave *) nxw()))
               if (!strcmp(os->GetName(), wrk->GetOrdinal())) break;
            if (!wrk) {
               if (!oo.IsNull()) oo += ",";
               oo += os->GetName();
            }
         }
         if (!oo.IsNull()) {
            Warning("ModifyWorkerLists", "worker(s) '%s' not found!", oo.Data());
            nwc = -2;
         }
      }
      // Cleanup hash list
      if (ords) {
         ords->Delete();
         SafeDelete(ords);
      }
   }

   // Rescan for unique workers
   if (rs)
      FindUniqueSlaves();

   // Forward the request one step down, if needed
   Int_t action = (add) ? (Int_t) kActivateWorker : (Int_t) kDeactivateWorker;
   if (fw) {
      if (fProtocol > 32) {
         TMessage mess(kPROOF_WORKERLISTS);
         mess << action << TString(ord);
         Broadcast(mess);
         Collect(kActive, fCollectTimeout);
         if (fStatus != 0) {
            nwc = (fStatus < nwc) ? fStatus : nwc;
            if (fStatus == -2) {
               if (gDebug > 0)
                  Warning("ModifyWorkerLists", "request not completely full filled");
            } else {
               Error("ModifyWorkerLists", "request failed");
            }
         }
      } else {
         TString oo(ord), o;
         if (oo.Contains(","))
            Warning("ModifyWorkerLists", "block request not supported by server: splitting into pieces ...");
         Int_t from = 0;
         while(oo.Tokenize(o, from, ",")) {
            TMessage mess(kPROOF_WORKERLISTS);
            mess << action << o;
            Broadcast(mess);
            Collect(kActive, fCollectTimeout);
         }
      }
   }
   // Done
   return nwc;
}

////////////////////////////////////////////////////////////////////////////////
/// Save current list of active workers

void TProof::SaveActiveList()
{
   if (!fActiveSlavesSaved.IsNull()) fActiveSlavesSaved = "";
   if (fInactiveSlaves->GetSize() == 0) {
      fActiveSlavesSaved = "*";
   } else {
      TIter nxw(fActiveSlaves);
      TSlave *wk = 0;
      while ((wk = (TSlave *)nxw())) { fActiveSlavesSaved += TString::Format("%s,", wk->GetOrdinal()); }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Restore saved list of active workers

Int_t TProof::RestoreActiveList()
{
   // Clear the current active list
   DeactivateWorker("*", kFALSE);
   // Restore the previous active list
   if (!fActiveSlavesSaved.IsNull())
      return ActivateWorker(fActiveSlavesSaved, kFALSE);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Start a PROOF session on a specific cluster. If cluster is 0 (the
/// default) then the PROOF Session Viewer GUI pops up and 0 is returned.
/// If cluster is "lite://" we start a PROOF-lite session.
/// If cluster is "" (empty string) then we connect to the cluster specified
/// by 'Proof.LocalDefault', defaulting to "lite://".
/// If cluster is "pod://" (case insensitive), then we connect to a PROOF cluster
/// managed by PROOF on Demand (PoD, http://pod.gsi.de ).
/// Via conffile a specific PROOF config file in the confir directory can be specified.
/// Use loglevel to set the default loging level for debugging.
/// The appropriate instance of TProofMgr is created, if not
/// yet existing. The instantiated TProof object is returned.
/// Use TProof::cd() to switch between PROOF sessions.
/// For more info on PROOF see the TProof ctor.

TProof *TProof::Open(const char *cluster, const char *conffile,
                                          const char *confdir, Int_t loglevel)
{
   const char *pn = "TProof::Open";

   // Make sure libProof and dependents are loaded and TProof can be created,
   // dependents are loaded via the information in the [system].rootmap file
   if (!cluster) {

      TPluginManager *pm = gROOT->GetPluginManager();
      if (!pm) {
         ::Error(pn, "plugin manager not found");
         return 0;
      }

      if (gROOT->IsBatch()) {
         ::Error(pn, "we are in batch mode, cannot show PROOF Session Viewer");
         return 0;
      }
      // start PROOF Session Viewer
      TPluginHandler *sv = pm->FindHandler("TSessionViewer", "");
      if (!sv) {
         ::Error(pn, "no plugin found for TSessionViewer");
         return 0;
      }
      if (sv->LoadPlugin() == -1) {
         ::Error(pn, "plugin for TSessionViewer could not be loaded");
         return 0;
      }
      sv->ExecPlugin(0);
      return 0;

   } else {

      TString clst(cluster);

      // Check for PoD cluster
      if (PoDCheckUrl( &clst ) < 0) return 0;

      if (clst.BeginsWith("workers=")) clst.Insert(0, "lite:///?");
      if (clst.BeginsWith("tunnel=")) clst.Insert(0, "/?");

      // Parse input URL
      TUrl u(clst);

      // *** GG, 060711: this does not seem to work any more (at XrdClient level)
      // *** to be investigated (it is not really needed; static tunnels work).
      // Dynamic tunnel:
      // Parse any tunning info ("<cluster>/?tunnel=[<tunnel_host>:]tunnel_port)
      TString opts(u.GetOptions());
      if (!opts.IsNull()) {
         Int_t it = opts.Index("tunnel=");
         if (it != kNPOS) {
            TString sport = opts(it + strlen("tunnel="), opts.Length());
            TString host("127.0.0.1");
            Int_t port = -1;
            Int_t ic = sport.Index(":");
            if (ic != kNPOS) {
               // Isolate the host
               host = sport(0, ic);
               sport.Remove(0, ic + 1);
            }
            if (!sport.IsDigit()) {
               // Remove the non digit part
               TRegexp re("[^0-9]");
               Int_t ind = sport.Index(re);
               if (ind != kNPOS)
                  sport.Remove(ind);
            }
            // Set the port
            if (sport.IsDigit())
               port = sport.Atoi();
            if (port > 0) {
               // Set the relevant variables
               ::Info("TProof::Open","using tunnel at %s:%d", host.Data(), port);
               gEnv->SetValue("XNet.SOCKS4Host", host);
               gEnv->SetValue("XNet.SOCKS4Port", port);
            } else {
               // Warn parsing problems
               ::Warning("TProof::Open",
                         "problems parsing tunnelling info from options: %s", opts.Data());
            }
         }
      }

      // Find out if we are required to attach to a specific session
      Int_t locid = -1;
      Bool_t create = kFALSE;
      if (opts.Length() > 0) {
         if (opts.BeginsWith("N",TString::kIgnoreCase)) {
            create = kTRUE;
            opts.Remove(0,1);
            u.SetOptions(opts);
         } else if (opts.IsDigit()) {
            locid = opts.Atoi();
         }
      }

      // Attach-to or create the appropriate manager
      TProofMgr *mgr = TProofMgr::Create(u.GetUrl());

      TProof *proof = 0;
      if (mgr && mgr->IsValid()) {

         // If XProofd we always attempt an attach first (unless
         // explicitly not requested).
         Bool_t attach = (create || mgr->IsProofd() || mgr->IsLite()) ? kFALSE : kTRUE;
         if (attach) {
            TProofDesc *d = 0;
            if (locid < 0)
               // Get the list of sessions
               d = (TProofDesc *) mgr->QuerySessions("")->First();
            else
               d = (TProofDesc *) mgr->GetProofDesc(locid);
            if (d) {
               proof = (TProof*) mgr->AttachSession(d);
               if (!proof || !proof->IsValid()) {
                  if (locid)
                     ::Error(pn, "new session could not be attached");
                  SafeDelete(proof);
               }
            }
         }

         // start the PROOF session
         if (!proof) {
            proof = (TProof*) mgr->CreateSession(conffile, confdir, loglevel);
            if (!proof || !proof->IsValid()) {
               ::Error(pn, "new session could not be created");
               SafeDelete(proof);
            }
         }
      }
      return proof;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get instance of the effective manager for 'url'
/// Return 0 on failure.

TProofMgr *TProof::Mgr(const char *url)
{
   if (!url)
      return (TProofMgr *)0;

   // Attach or create the relevant instance
   return TProofMgr::Create(url);
}

////////////////////////////////////////////////////////////////////////////////
/// Wrapper around TProofMgr::Reset(...).

void TProof::Reset(const char *url, Bool_t hard)
{
   if (url) {
      TProofMgr *mgr = TProof::Mgr(url);
      if (mgr && mgr->IsValid())
         mgr->Reset(hard);
      else
         ::Error("TProof::Reset",
                 "unable to initialize a valid manager instance");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get environemnt variables.

const TList *TProof::GetEnvVars()
{
   return fgProofEnvList;
}

////////////////////////////////////////////////////////////////////////////////
/// Add an variable to the list of environment variables passed to proofserv
/// on the master and slaves

void TProof::AddEnvVar(const char *name, const char *value)
{
   if (gDebug > 0) ::Info("TProof::AddEnvVar","%s=%s", name, value);

   if (fgProofEnvList == 0) {
      // initialize the list if needed
      fgProofEnvList = new TList;
      fgProofEnvList->SetOwner();
   } else {
      // replace old entries with the same name
      TObject *o = fgProofEnvList->FindObject(name);
      if (o != 0) {
         fgProofEnvList->Remove(o);
      }
   }
   fgProofEnvList->Add(new TNamed(name, value));
}

////////////////////////////////////////////////////////////////////////////////
/// Remove an variable from the list of environment variables passed to proofserv
/// on the master and slaves

void TProof::DelEnvVar(const char *name)
{
   if (fgProofEnvList == 0) return;

   TObject *o = fgProofEnvList->FindObject(name);
   if (o != 0) {
      fgProofEnvList->Remove(o);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the list of environment variables passed to proofserv
/// on the master and slaves

void TProof::ResetEnvVars()
{
   if (fgProofEnvList == 0) return;

   SafeDelete(fgProofEnvList);
}

////////////////////////////////////////////////////////////////////////////////
/// Save information about the worker set in the file .workers in the working
/// dir. Called each time there is a change in the worker setup, e.g. by
/// TProof::MarkBad().

void TProof::SaveWorkerInfo()
{
   // We must be masters
   if (TestBit(TProof::kIsClient))
      return;

   // We must have a server defined
   if (!gProofServ) {
      Error("SaveWorkerInfo","gProofServ undefined");
      return;
   }

   // The relevant lists must be defined
   if (!fSlaves && !fBadSlaves) {
      Warning("SaveWorkerInfo","all relevant worker lists is undefined");
      return;
   }

   // Create or truncate the file first
   TString fnwrk = gSystem->GetDirName(gProofServ->GetSessionDir())+"/.workers";
   FILE *fwrk = fopen(fnwrk.Data(),"w");
   if (!fwrk) {
      Error("SaveWorkerInfo",
            "cannot open %s for writing (errno: %d)", fnwrk.Data(), errno);
      return;
   }

   // Do we need to register an additional line for another log?
   TString addlogext;
   TString addLogTag;
   if (gSystem->Getenv("PROOF_ADDITIONALLOG")) {
      addlogext = gSystem->Getenv("PROOF_ADDITIONALLOG");
      TPMERegexp reLogTag("^__(.*)__\\.log");  // $
      if (reLogTag.Match(addlogext) == 2) {
         addLogTag = reLogTag[1];
      }
      else {
         addLogTag = "+++";
      }
      if (gDebug > 0)
         Info("SaveWorkerInfo", "request for additional line with ext: '%s'",  addlogext.Data());
   }

   // Used to eliminate datetime and PID from workdir to obtain log file name
   TPMERegexp re("(.*?)-[0-9]+-[0-9]+$");

   // Loop over the list of workers (active is any worker not flagged as bad)
   TIter nxa(fSlaves);
   TSlave *wrk = 0;
   TString logfile;
   while ((wrk = (TSlave *) nxa())) {
      Int_t status = (fBadSlaves && fBadSlaves->FindObject(wrk)) ? 0 : 1;
      logfile = wrk->GetWorkDir();
      if (re.Match(logfile) == 2) logfile = re[1];
      else continue;  // invalid (should not happen)
      // Write out record for this worker
      fprintf(fwrk,"%s@%s:%d %d %s %s.log\n",
                   wrk->GetUser(), wrk->GetName(), wrk->GetPort(), status,
                   wrk->GetOrdinal(), logfile.Data());
      // Additional line, if required
      if (addlogext.Length() > 0) {
         fprintf(fwrk,"%s@%s:%d %d %s(%s) %s.%s\n",
                     wrk->GetUser(), wrk->GetName(), wrk->GetPort(), status,
                     wrk->GetOrdinal(), addLogTag.Data(), logfile.Data(), addlogext.Data());
      }

   }

   // Loop also over the list of bad workers (if they failed to startup they are not in
   // the overall list
   TIter nxb(fBadSlaves);
   while ((wrk = (TSlave *) nxb())) {
      logfile = wrk->GetWorkDir();
      if (re.Match(logfile) == 2) logfile = re[1];
      else continue;  // invalid (should not happen)
      if (!fSlaves->FindObject(wrk)) {
         // Write out record for this worker
         fprintf(fwrk,"%s@%s:%d 0 %s %s.log\n",
                     wrk->GetUser(), wrk->GetName(), wrk->GetPort(),
                     wrk->GetOrdinal(), logfile.Data());
      }
   }

   // Eventually loop over the list of gracefully terminated workers: we'll get
   // logfiles from those workers as well. They'll be shown with a special
   // status of "2"
   TIter nxt(fTerminatedSlaveInfos);
   TSlaveInfo *sli;
   while (( sli = (TSlaveInfo *)nxt() )) {
      logfile = sli->GetDataDir();
      if (re.Match(logfile) == 2) logfile = re[1];
      else continue;  // invalid (should not happen)
      fprintf(fwrk, "%s 2 %s %s.log\n",
              sli->GetName(), sli->GetOrdinal(), logfile.Data());
      // Additional line, if required
      if (addlogext.Length() > 0) {
         fprintf(fwrk, "%s 2 %s(%s) %s.%s\n",
                 sli->GetName(), sli->GetOrdinal(), addLogTag.Data(),
                 logfile.Data(), addlogext.Data());
      }
   }

   // Close file
   fclose(fwrk);

   // We are done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the value from the specified parameter from the specified collection.
/// Returns -1 in case of error (i.e. list is 0, parameter does not exist
/// or value type does not match), 0 otherwise.

Int_t TProof::GetParameter(TCollection *c, const char *par, TString &value)
{
   TObject *obj = c ? c->FindObject(par) : (TObject *)0;
   if (obj) {
      TNamed *p = dynamic_cast<TNamed*>(obj);
      if (p) {
         value = p->GetTitle();
         return 0;
      }
   }
   return -1;

}

////////////////////////////////////////////////////////////////////////////////
/// Get the value from the specified parameter from the specified collection.
/// Returns -1 in case of error (i.e. list is 0, parameter does not exist
/// or value type does not match), 0 otherwise.

Int_t TProof::GetParameter(TCollection *c, const char *par, Int_t &value)
{
   TObject *obj = c ? c->FindObject(par) : (TObject *)0;
   if (obj) {
      TParameter<Int_t> *p = dynamic_cast<TParameter<Int_t>*>(obj);
      if (p) {
         value = p->GetVal();
         return 0;
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the value from the specified parameter from the specified collection.
/// Returns -1 in case of error (i.e. list is 0, parameter does not exist
/// or value type does not match), 0 otherwise.

Int_t TProof::GetParameter(TCollection *c, const char *par, Long_t &value)
{
   TObject *obj = c ? c->FindObject(par) : (TObject *)0;
   if (obj) {
      TParameter<Long_t> *p = dynamic_cast<TParameter<Long_t>*>(obj);
      if (p) {
         value = p->GetVal();
         return 0;
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the value from the specified parameter from the specified collection.
/// Returns -1 in case of error (i.e. list is 0, parameter does not exist
/// or value type does not match), 0 otherwise.

Int_t TProof::GetParameter(TCollection *c, const char *par, Long64_t &value)
{
   TObject *obj = c ? c->FindObject(par) : (TObject *)0;
   if (obj) {
      TParameter<Long64_t> *p = dynamic_cast<TParameter<Long64_t>*>(obj);
      if (p) {
         value = p->GetVal();
         return 0;
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the value from the specified parameter from the specified collection.
/// Returns -1 in case of error (i.e. list is 0, parameter does not exist
/// or value type does not match), 0 otherwise.

Int_t TProof::GetParameter(TCollection *c, const char *par, Double_t &value)
{
   TObject *obj = c ? c->FindObject(par) : (TObject *)0;
   if (obj) {
      TParameter<Double_t> *p = dynamic_cast<TParameter<Double_t>*>(obj);
      if (p) {
         value = p->GetVal();
         return 0;
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that dataset is in the form to be processed. This may mean
/// retrieving the relevant info from the dataset manager or from the
/// attached input list.
/// Returns 0 on success, -1 on error

Int_t TProof::AssertDataSet(TDSet *dset, TList *input,
                            TDataSetManager *mgr, TString &emsg)
{
   emsg = "";

   // We must have something to process
   if (!dset || !input || !mgr) {
      emsg.Form("invalid inputs (%p, %p, %p)", dset, input, mgr);
      return -1;
   }

   TList *datasets = new TList;
   TFileCollection *dataset = 0;
   TString lookupopt;
   TString dsname(dset->GetName());

   // First extract the "entry list" part on the global name, if any
   TString dsns(dsname), enlname;
   Ssiz_t eli = dsns.Index("?enl=");
   if (eli != kNPOS) {
      enlname = dsns(eli + strlen("?enl="), dsns.Length());
      dsns.Remove(eli, dsns.Length()-eli);
   }

   // The dataset maybe in the form of a TFileCollection in the input list
   if (dsname.BeginsWith("TFileCollection:")) {
      // Isolate the real name
      dsname.ReplaceAll("TFileCollection:", "");
      // Get the object
      dataset = (TFileCollection *) input->FindObject(dsname);
      if (!dataset) {
         emsg.Form("TFileCollection %s not found in input list", dset->GetName());
         return -1;
      }
      // Remove from everywhere
      input->RecursiveRemove(dataset);
      // Add it to the local list
      datasets->Add(new TPair(dataset, new TObjString(enlname.Data())));
      // Make sure we lookup everything (unless the client or the administrator
      // required something else)
      if (TProof::GetParameter(input, "PROOF_LookupOpt", lookupopt) != 0) {
         lookupopt = gEnv->GetValue("Proof.LookupOpt", "all");
         input->Add(new TNamed("PROOF_LookupOpt", lookupopt.Data()));
      }
   }

   // This is the name we parse for additional specifications, such directory
   // and object name; for multiple datasets we assume that the directory and
   // and object name are the same for all datasets
   TString dsnparse;
   // The received message included an empty dataset, with only the name
   // defined: assume that a dataset, stored on the PROOF master by that
   // name, should be processed.
   if (!dataset) {

      TFileCollection *fc = nullptr;

      // Check if the entry list and dataset name are valid. If they have spaces,
      // commas, or pipes, they are not considered as valid and we revert to the
      // "multiple datasets" case
      TRegexp rg("[, |]");
      Bool_t validEnl = (enlname.Index(rg) == kNPOS) ? kTRUE : kFALSE;
      Bool_t validSdsn = (dsns.Index(rg) == kNPOS) ? kTRUE : kFALSE;

      if (validEnl && validSdsn && (( fc = mgr->GetDataSet(dsns) ))) {

         //
         // String corresponds to ONE dataset only
         //

         TIter nxfi(fc->GetList());
         TFileInfo *fi;
         while (( fi = (TFileInfo *)nxfi() ))
            fi->SetTitle(dsns.Data());
         dataset = fc;
         dsnparse = dsns;  // without entry list

         // Adds the entry list (or empty string if not specified)
         datasets->Add( new TPair(dataset, new TObjString( enlname.Data() )) );

      } else {

         //
         // String does NOT correspond to one dataset: check if many datasets
         // were specified instead
         //

         dsns = dsname.Data();
         TString dsn1;
         Int_t from1 = 0;
         while (dsns.Tokenize(dsn1, from1, "[, ]")) {
            TString dsn2;
            Int_t from2 = 0;
            while (dsn1.Tokenize(dsn2, from2, "|")) {
               enlname = "";
               Int_t ienl = dsn2.Index("?enl=");
               if (ienl != kNPOS) {
                  enlname = dsn2(ienl + 5, dsn2.Length());
                  dsn2.Remove(ienl);
               }
               if ((fc = mgr->GetDataSet(dsn2.Data()))) {
                  // Save dataset name in TFileInfo's title to use it in TDset
                  TIter nxfi(fc->GetList());
                  TFileInfo *fi;
                  while ((fi = (TFileInfo *) nxfi())) { fi->SetTitle(dsn2.Data()); }
                  dsnparse = dsn2;
                  if (!dataset) {
                     // This is our dataset
                     dataset = fc;
                  } else {
                     // Add it to the dataset
                     dataset->Add(fc);
                     SafeDelete(fc);
                  }
               }
            }
            // The dataset name(s) in the first element
            if (dataset) {
               if (dataset->GetList()->First())
                  ((TFileInfo *)(dataset->GetList()->First()))->SetTitle(dsn1.Data());
               // Add it to the local list
               datasets->Add(new TPair(dataset, new TObjString(enlname.Data())));
            }
            // Reset the pointer
            dataset = 0;
         }

      }

      //
      // At this point the dataset(s) to be processed, if any, are found in the
      // "datasets" variable
      //

      if (!datasets || datasets->GetSize() <= 0) {
         emsg.Form("no dataset(s) found on the master corresponding to: %s", dsname.Data());
         return -1;
      } else {
         // Make 'dataset' to point to the first one in the list
         if (!(dataset = (TFileCollection *) ((TPair *)(datasets->First()))->Key())) {
            emsg.Form("dataset pointer is null: corruption? - aborting");
            return -1;
         }
      }
      // Apply the lookup option requested by the client or the administartor
      // (by default we trust the information in the dataset)
      if (TProof::GetParameter(input, "PROOF_LookupOpt", lookupopt) != 0) {
         lookupopt = gEnv->GetValue("Proof.LookupOpt", "stagedOnly");
         input->Add(new TNamed("PROOF_LookupOpt", lookupopt.Data()));
      }
   } else {
      // We were given a named, single, TFileCollection
      dsnparse = dsname;
   }

   // Logic for the subdir/obj names: try first to see if the dataset name contains
   // some info; if not check the settings in the TDSet object itself; if still empty
   // check the default tree name / path in the TFileCollection object; if still empty
   // use the default as the flow will determine
   TString dsTree;
   // Get the [subdir/]tree, if any
   mgr->ParseUri(dsnparse.Data(), 0, 0, 0, &dsTree);
   if (dsTree.IsNull()) {
      // Use what we have in the original dataset; we need this to locate the
      // meta data information
      dsTree += dset->GetDirectory();
      dsTree += dset->GetObjName();
   }
   if (!dsTree.IsNull() && dsTree != "/") {
      TString tree(dsTree);
      Int_t idx = tree.Index("/");
      if (idx != kNPOS) {
         TString dir = tree(0, idx+1);
         tree.Remove(0, idx);
         dset->SetDirectory(dir);
      }
      dset->SetObjName(tree);
   } else {
      // Use the default obj name from the TFileCollection
      dsTree = dataset->GetDefaultTreeName();
   }

   // Pass dataset server mapping instructions, if any
   TList *srvmapsref = TDataSetManager::GetDataSetSrvMaps();
   TList *srvmapslist = srvmapsref;
   TString srvmaps;
   if (TProof::GetParameter(input, "PROOF_DataSetSrvMaps", srvmaps) == 0) {
      srvmapslist = TDataSetManager::ParseDataSetSrvMaps(srvmaps);
      if (gProofServ) {
         TString msg;
         if (srvmapsref && !srvmapslist) {
            msg.Form("+++ Info: dataset server mapping(s) DISABLED by user");
         } else if (srvmapsref && srvmapslist && srvmapslist != srvmapsref) {
            msg.Form("+++ Info: dataset server mapping(s) modified by user");
         } else if (!srvmapsref && srvmapslist) {
            msg.Form("+++ Info: dataset server mapping(s) added by user");
         }
         gProofServ->SendAsynMessage(msg.Data());
      }
   }

   // Flag multi-datasets
   if (datasets->GetSize() > 1) dset->SetBit(TDSet::kMultiDSet);
   // Loop over the list of datasets
   TList *listOfMissingFiles = new TList;
   TEntryList *entrylist = 0;
   TPair *pair = 0;
   TIter nxds(datasets);
   while ((pair = (TPair *) nxds())) {
      // File Collection
      dataset = (TFileCollection *) pair->Key();
      // Entry list, if any
      TEntryList *enl = 0;
      TObjString *os = (TObjString *) pair->Value();
      if (strlen(os->GetName())) {
         if (!(enl = dynamic_cast<TEntryList *>(input->FindObject(os->GetName())))) {
            if (gProofServ)
               gProofServ->SendAsynMessage(TString::Format("+++ Warning:"
                                           " entry list %s not found", os->GetName()));
         }
         if (enl && (!(enl->GetLists()) || enl->GetLists()->GetSize() <= 0)) {
            if (gProofServ)
               gProofServ->SendAsynMessage(TString::Format("+++ Warning:"
                                           " no sub-lists in entry-list!"));
         }
      }
      TList *missingFiles = new TList;
      TSeqCollection* files = dataset->GetList();
      if (gDebug > 0) files->Print();
      Bool_t availableOnly = (lookupopt != "all") ? kTRUE : kFALSE;
      if (dset->TestBit(TDSet::kMultiDSet)) {
         TDSet *ds = new TDSet(dataset->GetName(), dset->GetObjName(), dset->GetDirectory());
         ds->SetSrvMaps(srvmapslist);
         if (!ds->Add(files, dsTree, availableOnly, missingFiles)) {
            emsg.Form("error integrating dataset %s", dataset->GetName());
            continue;
         }
         // Add the TDSet object to the multi-dataset
         dset->Add(ds);
         // Add entry list if any
         if (enl) ds->SetEntryList(enl);
      } else {
         dset->SetSrvMaps(srvmapslist);
         if (!dset->Add(files, dsTree, availableOnly, missingFiles)) {
            emsg.Form("error integrating dataset %s", dataset->GetName());
            continue;
         }
         if (enl) entrylist = enl;
      }
      if (missingFiles) {
         // The missing files objects have to be removed from the dataset
         // before delete.
         TIter next(missingFiles);
         TObject *file;
         while ((file = next())) {
            dataset->GetList()->Remove(file);
            listOfMissingFiles->Add(file);
         }
         missingFiles->SetOwner(kFALSE);
         missingFiles->Clear();
      }
      SafeDelete(missingFiles);
   }
   // Cleanup; we need to do this because pairs do no delete their content
   nxds.Reset();
   while ((pair = (TPair *) nxds())) {
      if (pair->Key()) delete pair->Key();
      if (pair->Value()) delete pair->Value();
   }
   datasets->SetOwner(kTRUE);
   SafeDelete(datasets);

   // Cleanup the server mapping list, if created by the user
   if (srvmapslist && srvmapslist != srvmapsref) {
      srvmapslist->SetOwner(kTRUE);
      SafeDelete(srvmapslist);
   }

   // Set the global entrylist, if required
   if (entrylist) dset->SetEntryList(entrylist);

   // Make sure it will be sent back merged with other similar lists created
   // during processing; this list will be transferred by the player to the
   // output list, once the latter has been created (see TProofPlayerRemote::Process)
   if (listOfMissingFiles && listOfMissingFiles->GetSize() > 0) {
      listOfMissingFiles->SetName("MissingFiles");
      input->Add(listOfMissingFiles);
   }

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Save input data file from 'cachedir' into the sandbox or create a the file
/// with input data objects

Int_t TProof::SaveInputData(TQueryResult *qr, const char *cachedir, TString &emsg)
{
   TList *input = 0;

   // We must have got something to process
   if (!qr || !(input = qr->GetInputList()) ||
       !cachedir || strlen(cachedir) <= 0) return 0;

   // There must be some input data or input data file
   TNamed *data = (TNamed *) input->FindObject("PROOF_InputDataFile");
   TList *inputdata = (TList *) input->FindObject("PROOF_InputData");
   if (!data && !inputdata) return 0;
   // Default dstination filename
   if (!data)
      input->Add((data = new TNamed("PROOF_InputDataFile", kPROOF_InputDataFile)));

   TString dstname(data->GetTitle()), srcname;
   Bool_t fromcache = kFALSE;
   if (dstname.BeginsWith("cache:")) {
      fromcache = kTRUE;
      dstname.ReplaceAll("cache:", "");
      srcname.Form("%s/%s", cachedir, dstname.Data());
      if (gSystem->AccessPathName(srcname)) {
         emsg.Form("input data file not found in cache (%s)", srcname.Data());
         return -1;
      }
   }

   // If from cache, just move the cache file
   if (fromcache) {
      if (gSystem->CopyFile(srcname, dstname, kTRUE) != 0) {
         emsg.Form("problems copying %s to %s", srcname.Data(), dstname.Data());
         return -1;
      }
   } else {
      // Create the file
      if (inputdata && inputdata->GetSize() > 0) {
         TFile *f = TFile::Open(dstname.Data(), "RECREATE");
         if (f) {
            f->cd();
            inputdata->Write();
            f->Close();
            delete f;
         } else {
            emsg.Form("could not create %s", dstname.Data());
            return -1;
         }
      } else {
         emsg.Form("no input data!");
         return -1;
      }
   }
   ::Info("TProof::SaveInputData", "input data saved to %s", dstname.Data());

   // Save the file name and clean up the data list
   data->SetTitle(dstname);
   if (inputdata) {
      input->Remove(inputdata);
      inputdata->SetOwner();
      delete inputdata;
   }

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Send the input data file to the workers

Int_t TProof::SendInputData(TQueryResult *qr, TProof *p, TString &emsg)
{
   TList *input = 0;

   // We must have got something to process
   if (!qr || !(input = qr->GetInputList())) return 0;

   // There must be some input data or input data file
   TNamed *inputdata = (TNamed *) input->FindObject("PROOF_InputDataFile");
   if (!inputdata) return 0;

   TString fname(inputdata->GetTitle());
   if (gSystem->AccessPathName(fname)) {
      emsg.Form("input data file not found in sandbox (%s)", fname.Data());
      return -1;
   }

   // PROOF session must available
   if (!p || !p->IsValid()) {
      emsg.Form("TProof object undefined or invalid: protocol error!");
      return -1;
   }

   // Send to unique workers and submasters
   p->BroadcastFile(fname, TProof::kBinary, "cache");

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the input data from the file defined in the input list

Int_t TProof::GetInputData(TList *input, const char *cachedir, TString &emsg)
{
   // We must have got something to process
   if (!input || !cachedir || strlen(cachedir) <= 0) return 0;

   // There must be some input data or input data file
   TNamed *inputdata = (TNamed *) input->FindObject("PROOF_InputDataFile");
   if (!inputdata) return 0;

   TString fname;
   fname.Form("%s/%s", cachedir, inputdata->GetTitle());
   if (gSystem->AccessPathName(fname)) {
      emsg.Form("input data file not found in cache (%s)", fname.Data());
      return -1;
   }

   // List of added objects (for proper cleaning ...)
   TList *added = new TList;
   added->SetName("PROOF_InputObjsFromFile");
   // Read the input data into the input list
   TFile *f = TFile::Open(fname.Data());
   if (f) {
      TList *keys = (TList *) f->GetListOfKeys();
      if (!keys) {
         emsg.Form("could not get list of object keys from file");
         return -1;
      }
      TIter nxk(keys);
      TKey *k = 0;
      while ((k = (TKey *)nxk())) {
         TObject *o = f->Get(k->GetName());
         if (o) {
            input->Add(o);
            added->Add(o);
         }
      }
      // Add the file as last one
      if (added->GetSize() > 0) {
         added->Add(f);
         input->Add(added);
      } else {
         // Cleanup the file now
         f->Close();
         delete f;
      }
   } else {
      emsg.Form("could not open %s", fname.Data());
      return -1;
   }

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Start the log viewer window usign the plugin manager

void TProof::LogViewer(const char *url, Int_t idx)
{
   if (!gROOT->IsBatch()) {
      // Get the handler, if not yet done
      if (!fgLogViewer) {
         if ((fgLogViewer =
            gROOT->GetPluginManager()->FindHandler("TProofProgressLog"))) {
            if (fgLogViewer->LoadPlugin() == -1) {
               fgLogViewer = 0;
               ::Error("TProof::LogViewer", "cannot load the relevant plug-in");
               return;
            }
         }
      }
      if (fgLogViewer) {
         // Execute the plug-in
         TString u = (url && strlen(url) <= 0) ? "lite" : url;
         fgLogViewer->ExecPlugin(2, u.Data(), idx);
      }
   } else {
      if (url && strlen(url) > 0) {
         ::Info("TProof::LogViewer",
                "batch mode: use TProofLog *pl = TProof::Mgr(\"%s\")->GetSessionLogs(%d)", url, idx);
      } else if (url && strlen(url) <= 0) {
         ::Info("TProof::LogViewer",
                "batch mode: use TProofLog *pl = TProof::Mgr(\"lite\")->GetSessionLogs(%d)", idx);
      } else {
         ::Info("TProof::LogViewer",
                "batch mode: use TProofLog *pl = TProof::Mgr(\"<master>\")->GetSessionLogs(%d)", idx);
      }
   }
   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/Disable the graphic progress dialog.
/// By default the dialog is enabled

void TProof::SetProgressDialog(Bool_t on)
{
   if (on)
      SetBit(kUseProgressDialog);
   else
      ResetBit(kUseProgressDialog);
}

////////////////////////////////////////////////////////////////////////////////
/// Show information about missing files during query described by 'qr' or the
/// last query if qr is null (default).
/// A short summary is printed in the end.

void TProof::ShowMissingFiles(TQueryResult *qr)
{
   TQueryResult *xqr = (qr) ? qr : GetQueryResult();
   if (!xqr) {
      Warning("ShowMissingFiles", "no (last) query found: do nothing");
      return;
   }

   // Get the list, if any
   TList *missing = (xqr->GetOutputList()) ? (TList *) xqr->GetOutputList()->FindObject("MissingFiles") : 0;
   if (!missing) {
      Info("ShowMissingFiles", "no files missing in query %s:%s", xqr->GetTitle(), xqr->GetName());
      return;
   }

   Int_t nmf = 0, ncf = 0;
   Long64_t msz = 0, mszzip = 0, mev = 0;
   // Scan the list
   TFileInfo *fi = 0;
   TIter nxf(missing);
   while ((fi = (TFileInfo *) nxf())) {
      char status = 'M';
      if (fi->TestBit(TFileInfo::kCorrupted)) {
         ncf++;
         status = 'C';
      } else {
         nmf++;
      }
      TFileInfoMeta *im = fi->GetMetaData();
      if (im) {
         if (im->GetTotBytes() > 0) msz += im->GetTotBytes();
         if (im->GetZipBytes() > 0) mszzip += im->GetZipBytes();
         mev += im->GetEntries();
         Printf(" %d. (%c) %s %s %lld", ncf+nmf, status, fi->GetCurrentUrl()->GetUrl(), im->GetName(), im->GetEntries());
      } else {
         Printf(" %d. (%c) %s '' -1", ncf+nmf, status, fi->GetCurrentUrl()->GetUrl());
      }
   }

   // Final notification
   if (msz <= 0) msz = -1;
   if (mszzip <= 0) mszzip = -1;
   Double_t xf = (Double_t)mev / (mev + xqr->GetEntries()) ;
   if (msz > 0. || mszzip > 0.) {
      Printf(" +++ %d file(s) missing, %d corrupted, i.e. %lld unprocessed events -->"
             " about %.2f%% of the total (%lld bytes, %lld zipped)",
             nmf, ncf, mev, xf * 100., msz, mszzip);
   } else {
      Printf(" +++ %d file(s) missing, %d corrupted, i.e. %lld unprocessed events -->"
             " about %.2f%% of the total", nmf, ncf, mev, xf * 100.);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get a TFileCollection with the files missing in the query described by 'qr'
/// or the last query if qr is null (default).
/// Return a null pointer if none were found, for whatever reason.
/// The caller is responsible for the returned object.

TFileCollection *TProof::GetMissingFiles(TQueryResult *qr)
{
   TFileCollection *fc = 0;

   TQueryResult *xqr = (qr) ? qr : GetQueryResult();
   if (!xqr) {
      Warning("GetMissingFiles", "no (last) query found: do nothing");
      return fc;
   }

   // Get the list, if any
   TList *missing = (xqr->GetOutputList()) ? (TList *) xqr->GetOutputList()->FindObject("MissingFiles") : 0;
   if (!missing) {
      if (gDebug > 0)
         Info("ShowMissingFiles", "no files missing in query %s:%s", xqr->GetTitle(), xqr->GetName());
      return fc;
   }

   // Create collection: name is <dsname>.m<j>, where 'j' is the first giving a non existing name
   TString fcname("unknown");
   TDSet *ds = (TDSet *) xqr->GetInputObject("TDSet");
   if (ds) {
      fcname.Form("%s.m0", ds->GetName());
      Int_t j = 1;
      while (gDirectory->FindObject(fcname) && j < 1000)
         fcname.Form("%s.m%d", ds->GetName(), j++);
   }
   fc = new TFileCollection(fcname, "Missing Files");
   if (ds) fc->SetDefaultTreeName(ds->GetObjName());
   // Scan the list
   TFileInfo *fi = 0;
   TIter nxf(missing);
   while ((fi = (TFileInfo *) nxf())) {
      fc->Add((TFileInfo *) fi->Clone());
   }
   fc->Update();
   // Done
   return fc;
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/Disable saving of the performance tree

void TProof::SetPerfTree(const char *pf, Bool_t withWrks)
{
   if (pf && strlen(pf) > 0) {
      fPerfTree = pf;
      SetParameter("PROOF_StatsHist", "");
      SetParameter("PROOF_StatsTrace", "");
      if (withWrks) SetParameter("PROOF_SlaveStatsTrace", "");
      Info("SetPerfTree", "saving of the performance tree enabled (%s)", fPerfTree.Data());
   } else {
      fPerfTree = "";
      DeleteParameters("PROOF_StatsHist");
      DeleteParameters("PROOF_StatsTrace");
      DeleteParameters("PROOF_SlaveStatsTrace");
      Info("SetPerfTree", "saving of the performance tree disabled");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save performance information from TPerfStats to file 'pf'.
/// If 'ref' is defined, do it for query 'ref'.
/// Return 0 on sucecss, -1 in case of any error

Int_t TProof::SavePerfTree(const char *pf, const char *ref)
{
   if (!IsValid()) {
      Error("SafePerfTree", "this TProof instance is invalid!");
      return -1;
   }

   TList *outls = GetOutputList();
   TString sref;
   if (ref && strlen(ref) > 0) {
      if (!fPlayer) {
         Error("SafePerfTree", "requested to use query '%s' but player instance undefined!", ref);
         return -1;
      }
      TQueryResult *qr = fPlayer->GetQueryResult(ref);
      if (!qr) {
         Error("SafePerfTree", "TQueryResult instance for query '%s' could not be retrieved", ref);
         return -1;
      }
      outls = qr->GetOutputList();
      sref.Form(" for requested query '%s'", ref);
   }
   if (!outls || (outls && outls->GetSize() <= 0)) {
      Error("SafePerfTree", "outputlist%s undefined or empty", sref.Data());
      return -1;
   }

   TString fn = fPerfTree;
   if (pf && strlen(pf)) fn = pf;
   if (fn.IsNull()) fn = "perftree.root";

   TFile f(fn, "RECREATE");
   if (f.IsZombie()) {
      Error("SavePerfTree", "could not open file '%s' for writing", fn.Data());
   } else {
      f.cd();
      TIter nxo(outls);
      TObject* obj = 0;
      while ((obj = nxo())) {
         TString objname(obj->GetName());
         if (objname.BeginsWith("PROOF_")) {
            // Must list the objects since other PROOF_ objects exist
            // besides timing objects
            if (objname == "PROOF_PerfStats" ||
                objname == "PROOF_PacketsHist" ||
                objname == "PROOF_EventsHist" ||
                objname == "PROOF_NodeHist" ||
                objname == "PROOF_LatencyHist" ||
                objname == "PROOF_ProcTimeHist" ||
                objname == "PROOF_CpuTimeHist")
               obj->Write();
         }
      }
      f.Close();
   }
   Info("SavePerfTree", "performance information%s saved in %s ...", sref.Data(), fn.Data());

   // Done
   return 0;
}
