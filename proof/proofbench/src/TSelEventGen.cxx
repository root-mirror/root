// @(#)root/proof:$Id$
// Author: Sangsu Ryu 22/06/2010

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TSelEventGen
\ingroup proofbench

Selector for event file generation.
List of files to be generated for each node is provided by client.
And list of files generated is sent back.
Existing files are reused if not forced to be regenerated.

*/

#define TSelEventGen_cxx

#include "TSelEventGen.h"
#include "TParameter.h"
#include "TProofNodeInfo.h"
#include "TProofBenchTypes.h"
#include "TProof.h"
#include "TMap.h"
#include "TDSet.h"
#include "TFileInfo.h"
#include "TFile.h"
#include "TSortedList.h"
#include "TRandom.h"
#include "Event.h"
#include "TProofServ.h"
#include "TMacro.h"

ClassImp(TSelEventGen);

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TSelEventGen::TSelEventGen()
             : fBaseDir(""), fNEvents(100000), fNTracks(100), fNTracksMax(-1),
               fRegenerate(kFALSE), fTotalGen(0), fFilesGenerated(0),
               fGenerateFun(0), fChain(0)
{
   if (gProofServ){
      fBaseDir=gProofServ->GetDataDir();
      // Two directories up
      fBaseDir.Remove(fBaseDir.Last('/'));
      fBaseDir.Remove(fBaseDir.Last('/'));
   }
   else{
      fBaseDir="";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// The Begin() function is called at the start of the query.
/// When running with PROOF Begin() is only called on the client.
/// The tree argument is deprecated (on PROOF 0 is passed).

void TSelEventGen::Begin(TTree *)
{
   TString option = GetOption();
   // Determine the test type
   TMap *filemap = dynamic_cast<TMap *>
                     (fInput->FindObject("PROOF_FilesToProcess"));
   if (filemap) {
      //Info("Begin", "dumping the file map:");
      //filemap->Print();
   } else {
      if (fInput->FindObject("PROOF_FilesToProcess")) {
         Error("Begin", "object 'PROOF_FilesToProcess' found but not a map"
              " (%s)", fInput->FindObject("PROOF_FilesToProcess")->ClassName());
      } else {
         Error("Begin", "object 'PROOF_FilesToProcess' not found");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// The SlaveBegin() function is called after the Begin() function.
/// When running with PROOF SlaveBegin() is called on each slave server.
/// The tree argument is deprecated (on PROOF 0 is passed).

void TSelEventGen::SlaveBegin(TTree *tree)
{
   Init(tree);

   TString option = GetOption();

   //get parameters

   Bool_t found_basedir=kFALSE;
   Bool_t found_nevents=kFALSE;
   Bool_t found_ntracks=kFALSE;
   Bool_t found_ntrkmax=kFALSE;
   Bool_t found_regenerate=kFALSE;

   TIter nxt(fInput);
   TString sinput;
   TObject *obj;

   while ((obj = nxt())){

      sinput=obj->GetName();
      //Info("SlaveBegin", "Input list: %s", sinput.Data());

      if (sinput.Contains("PROOF_BenchmarkBaseDir")){
         TNamed* a = dynamic_cast<TNamed*>(obj);
         if (a){
            TString bdir = a->GetTitle();
            if (!bdir.IsNull()){
               TUrl u(bdir, kTRUE);
               Bool_t isLocal = !strcmp(u.GetProtocol(), "file") ? kTRUE : kFALSE;
               if (isLocal && !gSystem->IsAbsoluteFileName(u.GetFile()))
                  u.SetFile(TString::Format("%s/%s", fBaseDir.Data(), u.GetFile()));
               if (isLocal) {
                  if ((gSystem->AccessPathName(u.GetFile()) &&
                     gSystem->mkdir(u.GetFile(), kTRUE) == 0) ||
                     !gSystem->AccessPathName(u.GetFile(), kWritePermission)) {
                     // Directory is writable
                     fBaseDir = u.GetFile();
                     Info("SlaveBegin", "using base directory \"%s\"", fBaseDir.Data());
                  } else{
                     // Directory is not writable or not available, use default directory
                     Warning("SlaveBegin", "\"%s\" directory is not writable or not existing,"
                           " using default directory: %s",
                           bdir.Data(), fBaseDir.Data());
                  }
               } else {
                  // We assume the user knows what it does
                  fBaseDir = bdir;
                  Info("SlaveBegin", "using non local base directory \"%s\"", fBaseDir.Data());
               }
            } else{
               Info("SlaveBegin", "using default directory: %s",
                                   fBaseDir.Data());
            }
            found_basedir=kTRUE;
         }
         else{
            Error("SlaveBegin", "PROOF_BenchmarkBaseDir not type TNamed");
         }
         continue;
      }
      if (sinput.Contains("PROOF_BenchmarkNEvents")){
         TParameter<Long64_t>* a=dynamic_cast<TParameter<Long64_t>*>(obj);
         if (a){
            fNEvents= a->GetVal();
            found_nevents=kTRUE;
         }
         else{
            Error("SlaveBegin", "PROOF_BenchmarkEvents not type TParameter"
                                "<Long64_t>*");
         }
         continue;
      }
      if (sinput.Contains("PROOF_BenchmarkNTracks")){
         TParameter<Int_t>* a=dynamic_cast<TParameter<Int_t>*>(obj);
         if (a){
            fNTracks=a->GetVal();
            found_ntracks=kTRUE;
         }
         else{
            Error("SlaveBegin", "PROOF_BenchmarkNTracks not type TParameter"
                                "<Int_t>*");
         }
         continue;
      }
      if (sinput.Contains("PROOF_BenchmarkNTracksMax")){
         TParameter<Int_t>* a=dynamic_cast<TParameter<Int_t>*>(obj);
         if (a){
            fNTracksMax=a->GetVal();
            found_ntrkmax=kTRUE;
         }
         else{
            Error("SlaveBegin", "PROOF_BenchmarkNTracksMax not type TParameter"
                                "<Int_t>*");
         }
         continue;
      }
      if (sinput.Contains("PROOF_BenchmarkRegenerate")){
         TParameter<Int_t>* a=dynamic_cast<TParameter<Int_t>*>(obj);
         if (a){
            fRegenerate=a->GetVal();
            found_regenerate=kTRUE;
         }
         else{
            Error("SlaveBegin", "PROOF_BenchmarkRegenerate not type TParameter"
                                "<Int_t>*");
         }
         continue;
      }
      if (sinput.Contains("PROOF_GenerateFun")){
         TNamed *a = dynamic_cast<TNamed*>(obj);
         if (!(fGenerateFun = dynamic_cast<TMacro *>(fInput->FindObject(a->GetTitle())))) {
            Error("SlaveBegin", "PROOF_GenerateFun requires the TMacro object in the input list");
         }
         continue;
      }
   }

   if (!found_basedir){
      Warning("SlaveBegin", "PROOF_BenchmarkBaseDir not found; using default:"
                            " %s", fBaseDir.Data());
   }
   if (!found_nevents){
      Warning("SlaveBegin", "PROOF_BenchmarkNEvents not found; using default:"
                            " %lld", fNEvents);
   }
   if (!found_ntracks){
      Warning("SlaveBegin", "PROOF_BenchmarkNTracks not found; using default:"
                            " %d", fNTracks);
   }
   if (!found_ntrkmax){
      Warning("SlaveBegin", "PROOF_BenchmarkNTracksMax not found; using default:"
                            " %d", fNTracksMax);
   } else if (fNTracksMax <= fNTracks) {
      Warning("SlaveBegin", "PROOF_BenchmarkNTracksMax must be larger then"
                            " fNTracks=%d ; ignoring", fNTracks);
      fNTracksMax = -1;
      found_ntrkmax = kFALSE;
   }
   if (!found_regenerate){
      Warning("SlaveBegin", "PROOF_BenchmarkRegenerate not found; using"
                            " default: %d", fRegenerate);
   }

   fFilesGenerated = new TList();

   TString hostname(TUrl(gSystem->HostName()).GetHostFQDN());
   TString thisordinal = gProofServ ? gProofServ->GetOrdinal() : "n.d";
   TString sfilegenerated =
      TString::Format("PROOF_FilesGenerated_%s_%s", hostname.Data(), thisordinal.Data());
   fFilesGenerated->SetName(sfilegenerated);
   fFilesGenerated->SetOwner(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
///Generate files for IO-bound run
///Input parameters
///   filename: The name of the file to be generated
///   sizenevents: Either the number of events to generate when
///                filetype==kPBFileBenchmark
///                or the size of the file to generate when
///                filetype==kPBFileCleanup
///Returns
///   Either Number of entries in the file when
///   filetype==kPBFileBenchmark
///   or bytes written when filetype==kPBFileCleanup
///return 0 in case error

Long64_t TSelEventGen::GenerateFiles(const char *filename, Long64_t sizenevents)
{
   Long64_t nentries=0;
   TDirectory* savedir = gDirectory;
   //printf("current dir=%s\n", gDirectory->GetPath());

   TFile *f = TFile::Open(filename, "RECREATE");

   savedir->cd();

   if (!f || f->IsZombie()) return 0;

   Event *event=new Event();
   Event *ep = event;
   TTree* eventtree= new TTree("EventTree", "Event Tree");
   eventtree->SetDirectory(f);

   const Int_t buffersize=32000;
   eventtree->Branch("event", "Event", &ep, buffersize, 1);
   eventtree->AutoSave();

   Long64_t i=0;
   Long64_t size_generated=0;

//   f->SetCompressionLevel(0); //no compression
   Int_t ntrks = fNTracks;

   Info("GenerateFiles", "Generating %s", filename);
   while (sizenevents--){
      //event->Build(i++,fNTracksBench,0);
      if (fNTracksMax > fNTracks) {
         // Required to smear the number of tracks between [min,max]
         ntrks = fNTracks + gRandom->Integer(fNTracksMax - fNTracks);
      }
      event->Build(i++, ntrks, 0);
      size_generated+=eventtree->Fill();
   }
   nentries=eventtree->GetEntries();
   Info("GenerateFiles", "%s generated with %lld entries", filename, nentries);
   savedir = gDirectory;

   f = eventtree->GetCurrentFile();
   f->cd();
   eventtree->Write();
   eventtree->SetDirectory(0);

   f->Close();
   delete f;
   f = 0;
   eventtree->Delete();
   event->Delete();
   savedir->cd();

   return nentries;
}

////////////////////////////////////////////////////////////////////////////////
/// The Process() function is called for each entry in the tree (or possibly
/// keyed object in the case of PROOF) to be processed. The entry argument
/// specifies which entry in the currently loaded tree is to be processed.
/// It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
/// to read either all or the required parts of the data. When processing
/// keyed objects with PROOF, the object is already loaded and is available
/// via the fObject pointer.
///
/// This function should contain the "body" of the analysis. It can contain
/// simple or elaborate selection criteria, run algorithms on the data
/// of the event and typically fill histograms.

Bool_t TSelEventGen::Process(Long64_t entry)
{
   // WARNING when a selector is used with a TChain, you must use
   //  the pointer to the current TTree to call GetEntry(entry).
   //  The entry is always the local entry number in the current tree.
   //  Assuming that fChain is the pointer to the TChain being processed,
   //  use fChain->GetTree()->GetEntry(entry).

   TDSetElement *fCurrent = 0;
   TPair *elemPair = 0;
   if (fInput && (elemPair = dynamic_cast<TPair *>
                      (fInput->FindObject("PROOF_CurrentElement")))) {
      if ((fCurrent = dynamic_cast<TDSetElement *>(elemPair->Value()))) {
         Info("Process", "entry %lld: file: '%s'", entry, fCurrent->GetName());
      } else {
         Error("Process", "entry %lld: no file specified!", entry);
         return kFALSE;
      }
   }

   // Generate
   TString filename(fCurrent->GetName());
   if (!fBaseDir.IsNull()) {
      if (fBaseDir.Contains("<fn>")) {
         filename = fBaseDir;
         filename.ReplaceAll("<fn>", fCurrent->GetName());
      } else {
         filename.Form("%s/%s", fBaseDir.Data(), fCurrent->GetName());
      }
   }
   TString fndset(filename);

   // Set the Url for remote access
   TString seed = TString::Format("%s/%s", gSystem->HostName(), filename.Data()), dsrv;
   TUrl basedirurl(filename, kTRUE);
   if (!strcmp(basedirurl.GetProtocol(), "file")) {
      TProofServ::GetLocalServer(dsrv);
      TProofServ::FilterLocalroot(fndset, dsrv);
   }

   //generate files
   Long64_t neventstogenerate = fNEvents;

   Long64_t entries_file=0;
   Long64_t filesize=0;
   Bool_t filefound=kFALSE;
   FileStat_t filestat;
   TUUID uuid;
   if (!fRegenerate && !gSystem->GetPathInfo(filename, filestat)) { //stat'ed
      TFile *f = TFile::Open(filename);
      if (f && !f->IsZombie()){
         TTree* t = (TTree *) f->Get("EventTree");
         if (t) {
            entries_file = t->GetEntries();
            if (entries_file == neventstogenerate) {
               // File size seems to be correct, skip generation
               Info("Process", "bench file (%s, entries=%lld) exists:"
                               " skipping generation.", filename.Data(), entries_file);
               filesize = f->GetSize();
               uuid = f->GetUUID();
               filefound = kTRUE;
            }
         }
         f->Close();
      }
      SafeDelete(f);
   }

   // Make sure there is enough space left of the device, if local
   TString bdir(fBaseDir);
   bdir.ReplaceAll("<fn>", "");
   if (!gSystem->AccessPathName(bdir)) {
      Long_t devid, devbsz, devbtot, devbfree;
      gSystem->GetFsInfo(bdir, &devid, &devbsz, &devbtot, &devbfree);
      // Must be more than 10% of space and at least 1 GB
      Long_t szneed = 1024 * 1024 * 1024, tomb = 1024 * 1024;
      if (devbfree * devbsz < szneed || devbfree < 0.1 * devbtot) {
         Error("Process", "not enough free space on device (%ld MB < {%ld, %ld} MB):"
                          " skipping generation of: %s",
                          (devbfree * devbsz) / tomb,
                          szneed / tomb, (Long_t) (0.1 * devbtot * devbsz / tomb),
                          filename.Data());
         fStatus = TSelector::kAbortFile;
      }
   }

   if (!filefound) {  // Generate
      gRandom->SetSeed(static_cast<UInt_t>(TMath::Hash(seed)));
      if (fGenerateFun) {
         TString fargs = TString::Format("\"%s\",%lld", filename.Data(), neventstogenerate);
         entries_file = (Long64_t) fGenerateFun->Exec(fargs);
      } else {
         entries_file = GenerateFiles(filename, neventstogenerate);
      }

      TFile *f = TFile::Open(filename);
      if (f && !f->IsZombie()) {
         filesize = f->GetSize();
         uuid = f->GetUUID();
         f->Close();
      } else {
         Error("Process", "can not open generated file: %s", filename.Data());
         fStatus = TSelector::kAbortFile;
         return kFALSE;
      }

      SafeDelete(f);
   }

   // Add meta data to the file info
   TFileInfoMeta* fimeta = new TFileInfoMeta("/EventTree", "TTree", entries_file);
   TMD5* md5 = 0;
   if (!strcmp(TUrl(filename,kTRUE).GetProtocol(), "file"))
      md5 = TMD5::FileChecksum(filename);
   TString md5s = (md5) ? md5->AsString() : "";
   TFileInfo *fi = new TFileInfo(TString::Format("%s%s", dsrv.Data(), fndset.Data()),
                                 filesize, uuid.AsString(), md5s.Data(), fimeta);
   SafeDelete(md5);

   // Mark it as staged
   fi->SetBit(TFileInfo::kStaged);

   // Add the fileinfo to the list
   if (fFilesGenerated) fFilesGenerated->Add(fi);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// The SlaveTerminate() function is called after all entries or objects
/// have been processed. When running with PROOF SlaveTerminate() is called
/// on each slave server

void TSelEventGen::SlaveTerminate()
{
   if (fFilesGenerated && fFilesGenerated->GetSize() > 0) {
      fOutput->Add(fFilesGenerated);
      Info("SlaveTerminate",
              "list '%s' of files generated by this worker added to the output list",
              fFilesGenerated->GetName());
   } else {
      if (!fFilesGenerated) {
         Warning("SlaveTerminate", "no list of generated files defined!");
      } else {
         Warning("SlaveTerminate", "list of generated files is empty!");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// The Terminate() function is the last function to be called during
/// a query. It always runs on the client, it can be used to present
/// the results graphically or save the results to file.

void TSelEventGen::Terminate()
{
}

////////////////////////////////////////////////////////////////////////////////

void TSelEventGen::Print(Option_t *) const
{
   Printf("fNEvents=%lld", fNEvents);
   Printf("fBaseDir=%s", fBaseDir.Data());
   Printf("fNTracks=%d", fNTracks);
   Printf("fRegenerate=%d", fRegenerate);
}

