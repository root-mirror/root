// @(#)root/proofplayer:$Id$
// Author: Maarten Ballintijn    18/03/02

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TPacketizer
\ingroup proofkernel

This class generates packets to be processed on PROOF worker servers.
A packet is an event range (begin entry and number of entries) or
object range (first object and number of objects) in a TTree
(entries) or a directory (objects) in a file.
Packets are generated taking into account the performance of the
remote machine, the time it took to process a previous packet on
the remote machine, the locality of the database files, etc.

*/

#include "TPacketizer.h"

#include "Riostream.h"
#include "TDSet.h"
#include "TEnv.h"
#include "TError.h"
#include "TEventList.h"
#include "TEntryList.h"
#include "TMap.h"
#include "TMessage.h"
#include "TMonitor.h"
#include "TNtupleD.h"
#include "TObject.h"
#include "TParameter.h"
#include "TPerfStats.h"
#include "TProofDebug.h"
#include "TProof.h"
#include "TProofPlayer.h"
#include "TProofServ.h"
#include "TSlave.h"
#include "TSocket.h"
#include "TTimer.h"
#include "TUrl.h"
#include "TClass.h"
#include "TMath.h"

//
// The following three utility classes manage the state of the
// work to be performed and the slaves involved in the process.
// A list of TFileNode(s) describes the hosts with files, each
// has a list of TFileStat(s) keeping the state for each TDSet
// element (file).
//
// The list of TSlaveStat(s) keep track of the work (being) done
// by each slave
//


//------------------------------------------------------------------------------

class TPacketizer::TFileStat : public TObject {

private:
   Bool_t         fIsDone;       // is this element processed
   TFileNode     *fNode;         // my FileNode
   TDSetElement  *fElement;      // location of the file and its range
   Long64_t       fNextEntry;    // cursor in the range, -1 when done

public:
   TFileStat(TFileNode *node, TDSetElement *elem);

   Bool_t         IsDone() const {return fIsDone;}
   void           SetDone() {fIsDone = kTRUE;}
   TFileNode     *GetNode() const {return fNode;}
   TDSetElement  *GetElement() const {return fElement;}
   Long64_t       GetNextEntry() const {return fNextEntry;}
   void           MoveNextEntry(Long64_t step) {fNextEntry += step;}
};


TPacketizer::TFileStat::TFileStat(TFileNode *node, TDSetElement *elem)
   : fIsDone(kFALSE), fNode(node), fElement(elem), fNextEntry(elem->GetFirst())
{
}


//------------------------------------------------------------------------------

class TPacketizer::TFileNode : public TObject {

private:
   TString        fNodeName;        // FQDN of the node
   TList         *fFiles;           // TDSetElements (files) stored on this node
   TObject       *fUnAllocFileNext; // cursor in fFiles
   TList         *fActFiles;        // files with work remaining
   TObject       *fActFileNext;     // cursor in fActFiles
   Int_t          fMySlaveCnt;      // number of slaves running on this node
   Int_t          fSlaveCnt;        // number of external slaves processing files on this node

public:
   TFileNode(const char *name);
   ~TFileNode() { delete fFiles; delete fActFiles; }

   void        IncMySlaveCnt() { fMySlaveCnt++; }
   void        IncSlaveCnt(const char *slave) { if (fNodeName != slave) fSlaveCnt++; }
   void        DecSlaveCnt(const char *slave) { if (fNodeName != slave) fSlaveCnt--; R__ASSERT(fSlaveCnt >= 0); }
   Int_t       GetSlaveCnt() const {return fMySlaveCnt + fSlaveCnt;}
   Int_t       GetNumberOfActiveFiles() const { return fActFiles->GetSize(); }
   Bool_t      IsSortable() const { return kTRUE; }

   const char *GetName() const { return fNodeName.Data(); }

   void Add(TDSetElement *elem)
   {
      TFileStat *f = new TFileStat(this,elem);
      fFiles->Add(f);
      if (fUnAllocFileNext == 0) fUnAllocFileNext = fFiles->First();
   }

   TFileStat *GetNextUnAlloc()
   {
      TObject *next = fUnAllocFileNext;

      if (next != 0) {
         // make file active
         fActFiles->Add(next);
         if (fActFileNext == 0) fActFileNext = fActFiles->First();

         // move cursor
         fUnAllocFileNext = fFiles->After(fUnAllocFileNext);
      }

      return (TFileStat *) next;
   }

   TFileStat *GetNextActive()
   {
      TObject *next = fActFileNext;

      if (fActFileNext != 0) {
         fActFileNext = fActFiles->After(fActFileNext);
         if (fActFileNext == 0) fActFileNext = fActFiles->First();
      }

      return (TFileStat *) next;
   }

   void RemoveActive(TFileStat *file)
   {
      if (fActFileNext == file) fActFileNext = fActFiles->After(file);
      fActFiles->Remove(file);
      if (fActFileNext == 0) fActFileNext = fActFiles->First();
   }

   Int_t Compare(const TObject *other) const
   {
      // Must return -1 if this is smaller than obj, 0 if objects are equal
      // and 1 if this is larger than obj.
      const TFileNode *obj = dynamic_cast<const TFileNode*>(other);
      if (!obj) {
         Error("Compare", "input is not a TPacketizer::TFileNode object");
         return 0;
      }

      Int_t myVal = GetSlaveCnt();
      Int_t otherVal = obj->GetSlaveCnt();
      if (myVal < otherVal) {
         return -1;
      } else if (myVal > otherVal) {
         return 1;
      } else {
         return 0;
      }
   }

   void Print(Option_t *) const
   {
      std::cout << "OBJ: " << IsA()->GetName() << "\t" << fNodeName
           << "\tMySlaveCount " << fMySlaveCnt
           << "\tSlaveCount " << fSlaveCnt << std::endl;
   }

   void Reset()
   {
      fUnAllocFileNext = fFiles->First();
      fActFiles->Clear();
      fActFileNext = 0;
      fSlaveCnt = 0;
      fMySlaveCnt = 0;
   }
};


TPacketizer::TFileNode::TFileNode(const char *name)
   : fNodeName(name), fFiles(new TList), fUnAllocFileNext(0),fActFiles(new TList),
     fActFileNext(0), fMySlaveCnt(0), fSlaveCnt(0)
{
   // Constructor

   fFiles->SetOwner();
   fActFiles->SetOwner(kFALSE);
}


//------------------------------------------------------------------------------

class TPacketizer::TSlaveStat : public TVirtualPacketizer::TVirtualSlaveStat {

friend class TPacketizer;

private:
   TFileNode     *fFileNode;     // corresponding node or 0
   TFileStat     *fCurFile;      // file currently being processed
   TDSetElement  *fCurElem;      // TDSetElement currently being processed
   TProofProgressStatus *AddProcessed(TProofProgressStatus *st);
public:
   TSlaveStat(TSlave *slave);
   ~TSlaveStat();

   TFileNode  *GetFileNode() const { return fFileNode; }

   void        SetFileNode(TFileNode *node) { fFileNode = node; }
};


TPacketizer::TSlaveStat::TSlaveStat(TSlave *slave)
   : fFileNode(0), fCurFile(0), fCurElem(0)
{
   fSlave = slave;
   fStatus = new TProofProgressStatus();
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup

TPacketizer::TSlaveStat::~TSlaveStat()
{
   SafeDelete(fStatus);
}

TProofProgressStatus *TPacketizer::TSlaveStat::AddProcessed(TProofProgressStatus *st)
{
   // Update the status info to the 'st'.
   // return the difference (*st - *fStatus)

   if (st) {
      // The entriesis not correct in 'st'
      Long64_t lastEntries = st->GetEntries() - fStatus->GetEntries();
      // The last proc time should not be added
      fStatus->SetLastProcTime(0.);
      // Get the diff
      TProofProgressStatus *diff = new TProofProgressStatus(*st - *fStatus);
      *fStatus += *diff;
      // Set the correct value
      fStatus->SetLastEntries(lastEntries);
      return diff;
   } else {
      Error("AddProcessed", "status arg undefined");
      return 0;
   }
}

//------------------------------------------------------------------------------

ClassImp(TPacketizer);

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TPacketizer::TPacketizer(TDSet *dset, TList *slaves, Long64_t first,
                         Long64_t num, TList *input, TProofProgressStatus *st)
            : TVirtualPacketizer(input, st)
{
   PDB(kPacketizer,1) Info("TPacketizer", "Enter (first %lld, num %lld)", first, num);

   // Init pointer members
   fPackets = 0;
   fUnAllocated = 0;
   fActive = 0;
   fFileNodes = 0;
   fMaxPerfIdx = 1;
   fMaxSlaveCnt = 0;
   fHeuristicPSiz = kFALSE;
   fDefMaxWrkNode = kTRUE;

   if (!fProgressStatus) {
      Error("TPacketizer", "No progress status");
      return;
   }

   Long_t maxSlaveCnt = 0;
   if (TProof::GetParameter(input, "PROOF_MaxSlavesPerNode", maxSlaveCnt) == 0) {
      if (maxSlaveCnt < 0) {
         Warning("TPacketizer", "PROOF_MaxSlavesPerNode must be positive");
         maxSlaveCnt = 0;
      }
      if (maxSlaveCnt > 0) fDefMaxWrkNode = kFALSE;
   } else {
      // Try also with Int_t (recently supported in TProof::SetParameter)
      Int_t mxslcnt = -1;
      if (TProof::GetParameter(input, "PROOF_MaxSlavesPerNode", mxslcnt) == 0) {
         if (mxslcnt < 0) {
            Warning("TPacketizer", "PROOF_MaxSlavesPerNode must be positive");
            mxslcnt = 0;
         }
         maxSlaveCnt = (Long_t) mxslcnt;
         if (maxSlaveCnt > 0) fDefMaxWrkNode = kFALSE;
      }
   }
   if (!maxSlaveCnt) {
      maxSlaveCnt = gEnv->GetValue("Packetizer.MaxWorkersPerNode", slaves->GetSize());
      if (maxSlaveCnt != slaves->GetSize()) fDefMaxWrkNode = kFALSE;
   }
   if (maxSlaveCnt > 0) {
      fMaxSlaveCnt = maxSlaveCnt;
      PDB(kPacketizer,1)
         Info("TPacketizer", "setting max number of workers per node to %ld", fMaxSlaveCnt);
   }

   fPackets = new TList;
   fPackets->SetOwner();

   fFileNodes = new TList;
   fFileNodes->SetOwner();
   fUnAllocated = new TList;
   fUnAllocated->SetOwner(kFALSE);
   fActive = new TList;
   fActive->SetOwner(kFALSE);


   fValid = kTRUE;

   // Resolve end-point urls to optmize distribution
   // dset->Lookup(); // moved to TProofPlayerRemote::Process

   // Split into per host entries
   dset->Reset();
   TDSetElement *e;
   while ((e = (TDSetElement*)dset->Next())) {
      if (e->GetValid()) continue;

      TUrl url = e->GetFileName();

      // Map non URL filenames to dummy host
      TString host;
      if ( !url.IsValid() ||
          (strncmp(url.GetProtocol(),"root", 4) &&
           strncmp(url.GetProtocol(),"file", 4)) ) {
         host = "no-host";
      } else if ( url.IsValid() && !strncmp(url.GetProtocol(),"file", 4)) {
         host = "localhost";
         url.SetProtocol("root");
      } else {
         host = url.GetHost();
      }
      // Get full name for local hosts
      if (host.Contains("localhost") || host == "127.0.0.1") {
         url.SetHost(gSystem->HostName());
         host = url.GetHostFQDN();
      }

      TFileNode *node = (TFileNode*) fFileNodes->FindObject( host );

      if (node == 0) {
         node = new TFileNode(host);
         fFileNodes->Add(node);
      }

      node->Add( e );
   }

   fSlaveStats = new TMap;
   fSlaveStats->SetOwner(kFALSE);

   // Record initial available workers
   Int_t nwrks = AddWorkers(slaves);
   Info("TPacketizer", "Initial number of workers: %d", nwrks);

   // Setup file & filenode structure
   Reset();
   // Optimize the number of files to be open when running on subsample
   Int_t validateMode = 0;
   Int_t gprc = TProof::GetParameter(input, "PROOF_ValidateByFile", validateMode);
   Bool_t byfile = (gprc == 0 && validateMode > 0 && num > -1) ? kTRUE : kFALSE;
   if (num > -1)
      PDB(kPacketizer,2)
         Info("TPacketizer",
              "processing subset of entries: validating by file? %s", byfile ? "yes": "no");
   ValidateFiles(dset, slaves, num, byfile);

   if (!fValid) return;

   // apply global range (first,num) to dset and rebuild structure
   // ommitting TDSet elements that are not needed

   Int_t files = 0;
   fTotalEntries = 0;
   fUnAllocated->Clear();  // avoid dangling pointers
   fActive->Clear();
   fFileNodes->Clear();    // then delete all objects
   PDB(kPacketizer,2) Info("TPacketizer", "processing range: first %lld, num %lld", first, num);

   dset->Reset();
   Long64_t cur = 0;
   while (( e = (TDSetElement*)dset->Next())) {

      // Skip invalid or missing file; It will be moved
      // from the dset to the 'MissingFiles' list in the player.
      if (!e->GetValid()) continue;

      // The dataset name, if any
      if (fDataSet.IsNull() && e->GetDataSet() && strlen(e->GetDataSet()))
         fDataSet = e->GetDataSet();

      TUrl url = e->GetFileName();
      Long64_t eFirst = e->GetFirst();
      Long64_t eNum = e->GetNum();
      PDB(kPacketizer,2)
         Info("TPacketizer", " --> '%s'", e->GetFileName());
      PDB(kPacketizer,2)
         Info("TPacketizer", " --> first %lld, num %lld (cur %lld)", eFirst, eNum, cur);

      if (!e->GetEntryList()){
         // this element is before the start of the global range, skip it
         if (cur + eNum < first) {
            cur += eNum;
            PDB(kPacketizer,2)
               Info("TPacketizer", " --> skip element cur %lld", cur);
            continue;
         }

         // this element is after the end of the global range, skip it
         if (num != -1 && (first+num <= cur)) {
            cur += eNum;
            PDB(kPacketizer,2)
               Info("TPacketizer", " --> drop element cur %lld", cur);
            continue; // break ??
         }

         Bool_t inRange = kFALSE;
         if (cur <= first || (num != -1 && (first+num <= cur+eNum))) {

            if (cur <= first) {
               // If this element contains the start of the global range
               // adjust its start and number of entries
               e->SetFirst( eFirst + (first - cur) );
               e->SetNum( e->GetNum() - (first - cur) );
               PDB(kPacketizer,2)
                  Info("TPacketizer", " --> adjust start %lld and end %lld",
                       eFirst + (first - cur), first + num - cur);
               inRange = kTRUE;
            }
            if (num != -1 && (first+num <= cur+eNum)) {
               // If this element contains the end of the global range
               // adjust its number of entries
               e->SetNum( first + num - e->GetFirst() - cur );
               PDB(kPacketizer,2)
                  Info("TPacketizer", " --> adjust end %lld", first + num - cur);
               inRange = kTRUE;
            }

         } else {
            // Increment the counter ...
            PDB(kPacketizer,2)
               Info("TPacketizer", " --> increment 'cur' by %lld", eNum);
            cur += eNum;
         }
         // Re-adjust eNum and cur, if needed
         if (inRange) {
            cur += eNum;
            eNum = e->GetNum();
         }

      } else {
         TEntryList *enl = dynamic_cast<TEntryList *>(e->GetEntryList());
         if (enl) {
            eNum = enl->GetN();
         } else {
            TEventList *evl = dynamic_cast<TEventList *>(e->GetEntryList());
            eNum = evl ? evl->GetN() : eNum;
         }
         if (!eNum)
            continue;
      }
      PDB(kPacketizer,2)
         Info("TPacketizer", " --> next cur %lld", cur);

      // Map non URL filenames to dummy host
      TString host;
      if ( !url.IsValid() ||
          (strncmp(url.GetProtocol(),"root", 4) &&
           strncmp(url.GetProtocol(),"file", 4)) ) {
         host = "no-host";
      } else if ( url.IsValid() && !strncmp(url.GetProtocol(),"file", 4)) {
         host = "localhost";
         url.SetProtocol("root");
      } else {
         host = url.GetHostFQDN();
      }
      // Get full name for local hosts
      if (host.Contains("localhost") || host == "127.0.0.1") {
         url.SetHost(gSystem->HostName());
         host = url.GetHostFQDN();
      }

      TFileNode *node = (TFileNode*) fFileNodes->FindObject( host );

      if ( node == 0 ) {
         node = new TFileNode( host );
         fFileNodes->Add( node );
      }

      ++files;
      fTotalEntries += eNum;
      node->Add(e);
      PDB(kPacketizer,2) e->Print("a");
   }

   PDB(kPacketizer,1)
      Info("TPacketizer", "processing %lld entries in %d files on %d hosts",
                         fTotalEntries, files, fFileNodes->GetSize());

   // Set the total number for monitoring
   if (gPerfStats)
      gPerfStats->SetNumEvents(fTotalEntries);

   Reset();

   if (fFileNodes->GetSize() == 0) {
      Info("TPacketizer", "no valid or non-empty file found: setting invalid");
      // No valid files: set invalid and return
      fValid = kFALSE;
      return;
   }

   // Below we provide a possibility to change the way packet size is
   // calculated or define the packet size directly.
   // fPacketAsAFraction can be interpreted as follows:
   // assuming all slaves have equal processing rate,
   // packet size is (#events processed by 1 slave) / fPacketSizeAsAFraction.
   // It substitutes 20 in the old formula to calculate the fPacketSize:
   // fPacketSize = fTotalEntries / (20 * nslaves)
   Long_t packetAsAFraction = 20;
   if (TProof::GetParameter(input, "PROOF_PacketAsAFraction", packetAsAFraction) == 0)
      Info("Process", "using alternate fraction of query time as a packet Size: %ld",
           packetAsAFraction);
   fPacketAsAFraction = (Int_t)packetAsAFraction;

   fPacketSize = 1;
   if (TProof::GetParameter(input, "PROOF_PacketSize", fPacketSize) == 0) {
      Info("Process","using alternate packet size: %lld", fPacketSize);
   } else {
      // Heuristic for starting packet size
      fHeuristicPSiz = kTRUE;
      Int_t nslaves = fSlaveStats->GetSize();
      if (nslaves > 0) {
         fPacketSize = fTotalEntries / (fPacketAsAFraction * nslaves);
         if (fPacketSize < 1) fPacketSize = 1;
      } else {
         fPacketSize = 1;
      }
   }

   PDB(kPacketizer,1) Info("TPacketizer", "Base Packetsize = %lld", fPacketSize);

   if (!fValid)
      SafeDelete(fProgress);

   PDB(kPacketizer,1) Info("TPacketizer", "Return");
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TPacketizer::~TPacketizer()
{
   if (fSlaveStats) {
      fSlaveStats->DeleteValues();
   }

   SafeDelete(fPackets);
   SafeDelete(fSlaveStats);
   SafeDelete(fUnAllocated);
   SafeDelete(fActive);
   SafeDelete(fFileNodes);
}

////////////////////////////////////////////////////////////////////////////////
/// Adds new workers. Returns the number of workers added, or -1 on failure.

Int_t TPacketizer::AddWorkers(TList *workers)
{
   if (!workers) {
      Error("AddWorkers", "Null list of new workers!");
      return -1;
   }

   Int_t curNumOfWrks = fSlaveStats->GetEntries();

   TSlave *sl;
   TIter next(workers);
   while (( sl = dynamic_cast<TSlave*>(next()) ))
      if (!fSlaveStats->FindObject(sl)) {
         fSlaveStats->Add(sl, new TSlaveStat(sl));
         fMaxPerfIdx = sl->GetPerfIdx() > fMaxPerfIdx ? sl->GetPerfIdx() : fMaxPerfIdx;
      }

   // If heuristic (and new workers) set the packet size
   Int_t nwrks = fSlaveStats->GetSize();
   if (fHeuristicPSiz && nwrks > curNumOfWrks) {
      if (nwrks > 0) {
         fPacketSize = fTotalEntries / (fPacketAsAFraction * nwrks);
         if (fPacketSize < 1) fPacketSize = 1;
      } else {
         fPacketSize = 1;
      }
   }

   // Update the max number that can access one file node if the default is used
   if (fDefMaxWrkNode && nwrks > fMaxSlaveCnt) fMaxSlaveCnt = nwrks;

   // Done
   return nwrks;
}

////////////////////////////////////////////////////////////////////////////////
/// Get next unallocated file.

TPacketizer::TFileStat *TPacketizer::GetNextUnAlloc(TFileNode *node)
{
   TFileStat *file = 0;

   if (node != 0) {
      file = node->GetNextUnAlloc();
      if (file == 0) RemoveUnAllocNode(node);
   } else {
      while (file == 0 && ((node = NextUnAllocNode()) != 0)) {
         file = node->GetNextUnAlloc();
         if (file == 0) RemoveUnAllocNode(node);
      }
   }

   if (file != 0) {
      // if needed make node active
      if (fActive->FindObject(node) == 0) {
         fActive->Add(node);
      }
   }

   return file;
}

////////////////////////////////////////////////////////////////////////////////
/// Get next unallocated node.

TPacketizer::TFileNode *TPacketizer::NextUnAllocNode()
{
   fUnAllocated->Sort();
   PDB(kPacketizer,2) {
      std::cout << "TPacketizer::NextUnAllocNode()" << std::endl;
      fUnAllocated->Print();
   }

   TFileNode *fn = (TFileNode*) fUnAllocated->First();
   if (fn != 0 && fMaxSlaveCnt > 0 && fn->GetSlaveCnt() >= fMaxSlaveCnt) {
      PDB(kPacketizer,1) Info("NextUnAllocNode", "reached workers per node limit (%ld)",
                              fMaxSlaveCnt);
      fn = 0;
   }

   return fn;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove unallocated node.

void TPacketizer::RemoveUnAllocNode(TFileNode * node)
{
   fUnAllocated->Remove(node);
}

////////////////////////////////////////////////////////////////////////////////
/// Get next active file.

TPacketizer::TFileStat *TPacketizer::GetNextActive()
{
   TFileNode *node;
   TFileStat *file = 0;

   while (file == 0 && ((node = NextActiveNode()) != 0)) {
         file = node->GetNextActive();
         if (file == 0) RemoveActiveNode(node);
   }

   return file;
}

////////////////////////////////////////////////////////////////////////////////
/// Get next active node.

TPacketizer::TFileNode *TPacketizer::NextActiveNode()
{
   fActive->Sort();
   PDB(kPacketizer,2) {
      Printf("TPacketizer::NextActiveNode : ----------------------");
      fActive->Print();
   }

   TFileNode *fn = (TFileNode*) fActive->First();
   if (fn != 0 && fMaxSlaveCnt > 0 && fn->GetSlaveCnt() >= fMaxSlaveCnt) {
      PDB(kPacketizer,1)
         Info("NextActiveNode", "reached workers per node limit (%ld)", fMaxSlaveCnt);
      fn = 0;
   }

   return fn;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove file from the list of actives.

void TPacketizer::RemoveActive(TFileStat *file)
{
   TFileNode *node = file->GetNode();

   node->RemoveActive(file);
   if (node->GetNumberOfActiveFiles() == 0) RemoveActiveNode(node);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove node from the list of actives.

void TPacketizer::RemoveActiveNode(TFileNode *node)
{
   fActive->Remove(node);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the internal datastructure for packet distribution.

void TPacketizer::Reset()
{
   fUnAllocated->Clear();
   fUnAllocated->AddAll(fFileNodes);

   fActive->Clear();

   TIter files(fFileNodes);
   TFileNode *fn;
   while ((fn = (TFileNode*) files.Next()) != 0) {
      fn->Reset();
   }

   TIter slaves(fSlaveStats);
   TObject *key;
   while ((key = slaves.Next()) != 0) {
      TSlaveStat *slstat = (TSlaveStat*) fSlaveStats->GetValue(key);
      if (slstat) {
         fn = (TFileNode*) fFileNodes->FindObject(slstat->GetName());
         if (fn != 0 ) {
            slstat->SetFileNode(fn);
            fn->IncMySlaveCnt();
         }
         slstat->fCurFile = 0;
      } else {
         Warning("Reset", "TSlaveStat associated to key '%s' is NULL", key->GetName());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check existence of file/dir/tree an get number of entries.
/// Assumes the files have been setup.

void TPacketizer::ValidateFiles(TDSet *dset, TList *slaves, Long64_t maxent, Bool_t byfile)
{
   TMap     slaves_by_sock;
   TMonitor mon;
   TList    workers;


   // Setup the communication infrastructure

   workers.AddAll(slaves);
   TIter    si(slaves);
   TSlave  *slm = 0;
   while ((slm = (TSlave*)si.Next()) != 0) {
      PDB(kPacketizer,3)
         Info("ValidateFiles","socket added to monitor: %p (%s)",
              slm->GetSocket(), slm->GetName());
      mon.Add(slm->GetSocket());
      slaves_by_sock.Add(slm->GetSocket(), slm);
      PDB(kPacketizer,1)
         Info("ValidateFiles",
              "mon: %p, wrk: %p, sck: %p", &mon, slm, slm->GetSocket());
   }

   mon.DeActivateAll();

   ((TProof*)gProof)->DeActivateAsyncInput();

   // Some monitoring systems (TXSocketHandler) need to know this
   ((TProof*)gProof)->fCurrentMonitor = &mon;

   // Preparing for client notification
   TString msg("Validating files");
   UInt_t n = 0;
   UInt_t tot = dset->GetListOfElements()->GetSize();
   Bool_t st = kTRUE;

   Long64_t totent = 0, nopenf = 0;
   while (kTRUE) {

      // send work
      while( TSlave *s = (TSlave*)workers.First() ) {

         workers.Remove(s);

         // find a file

         TSlaveStat *slstat = (TSlaveStat*)fSlaveStats->GetValue(s);
         if (!slstat) {
            Error("ValidateFiles", "TSlaveStat associated to slave '%s' is NULL", s->GetName());
            continue;
         }
         TFileNode *node = 0;
         TFileStat *file = 0;

         // try its own node first
         if ( (node = slstat->GetFileNode()) != 0 ) {
            file = GetNextUnAlloc(node);
            if ( file == 0 ) {
               slstat->SetFileNode(0);
            }
         }

         // look for a file on any other node if necessary
         if (file == 0) {
            file = GetNextUnAlloc();
         }

         if ( file != 0 ) {
            // files are done right away
            RemoveActive(file);

            slstat->fCurFile = file;
            TDSetElement *elem = file->GetElement();
            Long64_t entries = elem->GetEntries(kTRUE, kFALSE);
            if (entries < 0 || strlen(elem->GetTitle()) <= 0) {
               // This is decremented when we get the reply
               file->GetNode()->IncSlaveCnt(slstat->GetName());
               TMessage m(kPROOF_GETENTRIES);
               m << dset->IsTree()
               << TString(elem->GetFileName())
               << TString(elem->GetDirectory())
               << TString(elem->GetObjName());

               s->GetSocket()->Send( m );
               mon.Activate(s->GetSocket());
               PDB(kPacketizer,2)
                  Info("ValidateFiles",
                       "sent to worker-%s (%s) via %p GETENTRIES on %s %s %s %s",
                       s->GetOrdinal(), s->GetName(), s->GetSocket(),
                       dset->IsTree() ? "tree" : "objects", elem->GetFileName(),
                       elem->GetDirectory(), elem->GetObjName());
            } else {
               // Fill the info
               elem->SetTDSetOffset(entries);
               if (entries > 0) {
                  // Most likely valid
                  elem->SetValid();
                  if (!elem->GetEntryList()) {
                     if (elem->GetFirst() > entries) {
                        Error("ValidateFiles",
                              "first (%lld) higher then number of entries (%lld) in %s",
                              elem->GetFirst(), entries, elem->GetFileName());
                        // disable element
                        slstat->fCurFile->SetDone();
                        elem->Invalidate();
                        dset->SetBit(TDSet::kSomeInvalid);
                     }
                     if (elem->GetNum() == -1) {
                        elem->SetNum(entries - elem->GetFirst());
                     } else if (elem->GetFirst() + elem->GetNum() > entries) {
                        Warning("ValidateFiles", "num (%lld) + first (%lld) larger then number of"
                                 " keys/entries (%lld) in %s", elem->GetNum(), elem->GetFirst(),
                                 entries, elem->GetFileName());
                        elem->SetNum(entries - elem->GetFirst());
                     }
                     PDB(kPacketizer,2)
                        Info("ValidateFiles",
                             "found elem '%s' with %lld entries", elem->GetFileName(), entries);
                  }
               }
               // Notify the client
               n++;
               gProof->SendDataSetStatus(msg, n, tot, st);

               // This worker is ready for the next validation
               workers.Add(s);
            }
         }
      }

      // Check if there is anything to wait for
      if (mon.GetActive() == 0) {
         if (byfile && maxent > 0 && totent > 0) {
            // How many files do we still need ?
            Long64_t nrestf = (maxent - totent) * nopenf / totent ;
            if (nrestf <= 0 && maxent > totent) nrestf = 1;
            if (nrestf > 0) {
               PDB(kPacketizer,3)
                  Info("ValidateFiles", "{%lld, %lld, %lld): needs to validate %lld more files",
                                         maxent, totent, nopenf, nrestf);
               si.Reset();
               while ((slm = (TSlave *) si.Next()) && nrestf--) {
                  workers.Add(slm);
               }
               continue;
            } else {
               PDB(kPacketizer,3)
                  Info("ValidateFiles", "no need to validate more files");
               break;
            }
         } else {
            break;
         }
      }

      PDB(kPacketizer,3) {
         Info("ValidateFiles", "waiting for %d workers:", mon.GetActive());
         TList *act = mon.GetListOfActives();
         TIter next(act);
         TSocket *s = 0;
         while ((s = (TSocket*) next())) {
            Info("ValidateFiles", "found sck: %p", s);
            TSlave *sl = (TSlave *) slaves_by_sock.GetValue(s);
            if (sl)
               Info("ValidateFiles", "   worker-%s (%s)", sl->GetOrdinal(), sl->GetName());
         }
         delete act;
      }

      TSocket *sock = mon.Select();
      // If we have been interrupted break
      if (!sock) {
         Error("ValidateFiles", "selection has been interrupted - STOP");
         mon.DeActivateAll();
         fValid = kFALSE;
         break;
      }
      mon.DeActivate(sock);

      PDB(kPacketizer,3) Info("ValidateFiles", "select returned: %p", sock);

      TSlave *slave = (TSlave *) slaves_by_sock.GetValue( sock );
      if (!sock->IsValid()) {
         // A socket got invalid during validation
         Error("ValidateFiles", "worker-%s (%s) got invalid - STOP",
               slave->GetOrdinal(), slave->GetName());
         ((TProof*)gProof)->MarkBad(slave);
         fValid = kFALSE;
         break;
      }

      TMessage *reply;

      if ( sock->Recv(reply) <= 0 ) {
         // Help! lost a slave?
         ((TProof*)gProof)->MarkBad(slave);
         fValid = kFALSE;
         Error("ValidateFiles", "Recv failed! for worker-%s (%s)",
               slave->GetOrdinal(), slave->GetName());
         continue;
      }

      if (reply->What() != kPROOF_GETENTRIES) {
         // Not what we want: handover processing to the central machinery
         Int_t what = reply->What();
         ((TProof*)gProof)->HandleInputMessage(slave, reply);
         if (what == kPROOF_FATAL) {
             Error("ValidateFiles", "kPROOF_FATAL from worker-%s (%s)",
                                    slave->GetOrdinal(), slave->GetName());
             fValid = kFALSE;
         } else {
            // Reactivate the socket
            mon.Activate(sock);
         }
         // Get next message
         continue;
      }

      TSlaveStat *slavestat = (TSlaveStat*) fSlaveStats->GetValue( slave );
      TDSetElement *e = slavestat->fCurFile->GetElement();
      slavestat->fCurFile->GetNode()->DecSlaveCnt(slavestat->GetName());
      Long64_t entries;

      (*reply) >> entries;

      // Extract object name, if there
      if ((reply->BufferSize() > reply->Length())) {
         TString objname;
         (*reply) >> objname;
         e->SetTitle(objname);
      }

      e->SetTDSetOffset(entries);
      if ( entries > 0 ) {

         // This dataset element is most likely valid
         e->SetValid();

         //if (!e->GetEventList()) {
         if (!e->GetEntryList()){
            if ( e->GetFirst() > entries ) {
               Error("ValidateFiles", "first (%lld) higher then number of entries (%lld) in %s",
                                      e->GetFirst(), entries, e->GetFileName());

               // Invalidate the element
               slavestat->fCurFile->SetDone();
               e->Invalidate();
               dset->SetBit(TDSet::kSomeInvalid);
            }

            if ( e->GetNum() == -1 ) {
               e->SetNum( entries - e->GetFirst() );
            } else if ( e->GetFirst() + e->GetNum() > entries ) {
               Error("ValidateFiles",
                     "num (%lld) + first (%lld) larger then number of keys/entries (%lld) in %s",
                     e->GetNum(), e->GetFirst(), entries, e->GetFileName());
               e->SetNum(entries - e->GetFirst());
            }
         }

         // Count
         totent += entries;
         nopenf++;

         // Notify the client
         n++;
         gProof->SendDataSetStatus(msg, n, tot, st);

      } else {

         Error("ValidateFiles", "cannot get entries for %s (", e->GetFileName() );
         //
         // Need to fix this with a user option to allow incomplete file sets (rdm)
         //
         //fValid = kFALSE; // all element must be readable!
         if (gProofServ) {
            TMessage m(kPROOF_MESSAGE);
            m << TString(Form("Cannot get entries for file: %s - skipping", e->GetFileName()));
            gProofServ->GetSocket()->Send(m);
         }

         // Invalidate the element
         e->Invalidate();
         dset->SetBit(TDSet::kSomeInvalid);
      }
      PDB(kPacketizer,3) Info("ValidateFiles", " %lld events validated", totent);

      // Ready for the next job, unless we have enough files
      if (maxent < 0 || ((totent < maxent) && !byfile))
         workers.Add(slave);
   }

   // report std. output from slaves??

   ((TProof*)gProof)->ActivateAsyncInput();

   // This needs to be reset
   ((TProof*)gProof)->fCurrentMonitor = 0;

   // No reason to continue if invalid
   if (!fValid)
      return;


   // compute the offset for each file element
   Long64_t offset = 0;
   Long64_t newOffset = 0;
   TIter next(dset->GetListOfElements());
   TDSetElement *el;
   while ( (el = dynamic_cast<TDSetElement*> (next())) ) {
      newOffset = offset + el->GetTDSetOffset();
      el->SetTDSetOffset(offset);
      offset = newOffset;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get entries processed by the specified slave.

Long64_t TPacketizer::GetEntriesProcessed(TSlave *slave) const
{
   if ( fSlaveStats == 0 ) return 0;

   TSlaveStat *slstat = (TSlaveStat*) fSlaveStats->GetValue( slave );

   if ( slstat == 0 ) return 0;

   return slstat->GetEntriesProcessed();
}

////////////////////////////////////////////////////////////////////////////////
/// Get Estimation of the current rate; just summing the current rates of
/// the active workers

Float_t TPacketizer::GetCurrentRate(Bool_t &all)
{
   all = kTRUE;
   // Loop over the workers
   Float_t currate = 0.;
   if (fSlaveStats && fSlaveStats->GetSize() > 0) {
      TIter nxw(fSlaveStats);
      TObject *key;
      while ((key = nxw()) != 0) {
         TSlaveStat *slstat = (TSlaveStat *) fSlaveStats->GetValue(key);
         if (slstat && slstat->GetProgressStatus() && slstat->GetEntriesProcessed() > 0) {
            // Sum-up the current rates
            currate += slstat->GetProgressStatus()->GetCurrentRate();
         } else {
            all = kFALSE;
         }
      }
   }
   // Done
   return currate;
}

////////////////////////////////////////////////////////////////////////////////
/// Get next packet

TDSetElement *TPacketizer::GetNextPacket(TSlave *sl, TMessage *r)
{
   if ( !fValid ) {
      return 0;
   }

   // Find worker

   TSlaveStat *slstat = (TSlaveStat*) fSlaveStats->GetValue( sl );

   R__ASSERT( slstat != 0 );

   PDB(kPacketizer,1)
      Info("GetNextPacket","worker-%s (%s)", sl->GetOrdinal(), sl->GetName());
   // update stats & free old element

   Bool_t firstPacket = kFALSE;
   if ( slstat->fCurElem != 0 ) {
      Double_t latency = 0., proctime = 0., proccpu = 0.;
      Long64_t bytesRead = -1;
      Long64_t totalEntries = -1;
      Long64_t totev = 0;
      Long64_t numev = slstat->fCurElem->GetNum();

      fPackets->Add(slstat->fCurElem);

      if (sl->GetProtocol() > 18) {
         TProofProgressStatus *status = 0;
         (*r) >> latency;
         (*r) >> status;

         // Calculate the progress made in the last packet
         TProofProgressStatus *progress = 0;
         if (status) {
            // upadte the worker status
            numev = status->GetEntries() - slstat->GetEntriesProcessed();
            progress = slstat->AddProcessed(status);
            if (progress) {
               // (*fProgressStatus) += *progress;
               proctime = progress->GetProcTime();
               proccpu  = progress->GetCPUTime();
               totev  = status->GetEntries(); // for backward compatibility
               bytesRead  = progress->GetBytesRead();
               delete progress;
            }
            delete status;
         } else
             Error("GetNextPacket", "no status came in the kPROOF_GETPACKET message");
      } else {

         (*r) >> latency >> proctime >> proccpu;

         // only read new info if available
         if (r->BufferSize() > r->Length()) (*r) >> bytesRead;
         if (r->BufferSize() > r->Length()) (*r) >> totalEntries;
         if (r->BufferSize() > r->Length()) (*r) >> totev;

         numev = totev - slstat->GetEntriesProcessed();
         if (numev > 0)  slstat->GetProgressStatus()->IncEntries(numev);
         if (bytesRead > 0) slstat->GetProgressStatus()->IncBytesRead(bytesRead);
         if (numev > 0 || bytesRead > 0) slstat->GetProgressStatus()->SetLastUpdate();
      }

      if (fProgressStatus) {
         if (numev > 0)  fProgressStatus->IncEntries(numev);
         if (bytesRead > 0)  fProgressStatus->IncBytesRead(bytesRead);
         if (numev > 0 || bytesRead > 0) fProgressStatus->SetLastUpdate();
      }
      PDB(kPacketizer,2)
         Info("GetNextPacket","worker-%s (%s): %lld %7.3lf %7.3lf %7.3lf %lld",
                              sl->GetOrdinal(), sl->GetName(),
                              numev, latency, proctime, proccpu, bytesRead);

      if (gPerfStats)
         gPerfStats->PacketEvent(sl->GetOrdinal(), sl->GetName(), slstat->fCurElem->GetFileName(),
                                 numev, latency, proctime, proccpu, bytesRead);

      slstat->fCurElem = 0;
      if (fProgressStatus && fProgressStatus->GetEntries() == fTotalEntries) {
         HandleTimer(0);   // Send last timer message
         delete fProgress; fProgress = 0;
      }
   } else {
      firstPacket = kTRUE;
   }

   if ( fStop ) {
      HandleTimer(0);
      return 0;
   }

   // get a file if needed

   TFileStat *file = slstat->fCurFile;

   if ( file != 0 && file->IsDone() ) {
      file->GetNode()->DecSlaveCnt(slstat->GetName());
      if (gPerfStats)
         gPerfStats->FileEvent(sl->GetOrdinal(), sl->GetName(), file->GetNode()->GetName(),
                               file->GetElement()->GetFileName(), kFALSE);
      file = 0;
   }
   // Reset the current file field
   slstat->fCurFile = file;

   if (!file) {

      // Try its own node first
      if (slstat->GetFileNode() != 0) {
         file = GetNextUnAlloc(slstat->GetFileNode());
         if (!file) {
            slstat->SetFileNode(0);
         }
      }

      // try to find an unused filenode first
      if (!file) {
         file = GetNextUnAlloc();
      }

      // then look at the active filenodes
      if (!file) {
         file = GetNextActive();
      }

      if (!file) return 0;

      slstat->fCurFile = file;
      file->GetNode()->IncSlaveCnt(slstat->GetName());
      if (gPerfStats)
         gPerfStats->FileEvent(sl->GetOrdinal(), sl->GetName(),
                               file->GetNode()->GetName(),
                               file->GetElement()->GetFileName(), kTRUE);
   }

   // get a packet

   TDSetElement *base = file->GetElement();
   Long64_t num = Long64_t(fPacketSize*(Float_t)slstat->fSlave->GetPerfIdx()/fMaxPerfIdx);
   if (num < 1) num = 1;

   Long64_t first = file->GetNextEntry();
   Long64_t last = base->GetFirst() + base->GetNum();

   if ( first + num >= last ) {
      num = last - first;
      file->SetDone(); // done

      // delete file from active list (unalloc list is single pass, no delete needed)
      RemoveActive(file);

   } else {
      file->MoveNextEntry(num);
   }


   slstat->fCurElem = CreateNewPacket(base, first, num);
   if (base->GetEntryList())
      slstat->fCurElem->SetEntryList(base->GetEntryList(), first, num);

   // Flag the first packet of a new run (dataset)
   if (firstPacket)
      slstat->fCurElem->SetBit(TDSetElement::kNewRun);
   else
      slstat->fCurElem->ResetBit(TDSetElement::kNewRun);

   PDB(kPacketizer,2)
      Info("GetNextPacket","%s: %s %lld %lld", sl->GetOrdinal(), base->GetFileName(), first, num);

   return slstat->fCurElem;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of workers still processing

Int_t TPacketizer::GetActiveWorkers()
{
   Int_t actw = 0;
   TIter nxw(fSlaveStats);
   TObject *key;
   while ((key = nxw())) {
      TSlaveStat *wrkstat = (TSlaveStat *) fSlaveStats->GetValue(key);
      if (wrkstat && wrkstat->fCurFile) actw++;
   }
   // Done
   return actw;
}
