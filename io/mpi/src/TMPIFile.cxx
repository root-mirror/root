// @(#)root/io:$Id$
// Author: Amit Bashyal, August 2018

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TMPIFile.h"
#include "TFileCacheWrite.h"
#include "TKey.h"
#include "THashTable.h"
#include "TMath.h"

ClassImp(TMPIFile);

/** \class TMPIFile
A TFile implimentation that uses MPI to enable multiple parallel MPI
processes to write to a single file.

### Example usage

Begin_Macro (source)
{
   Int_t n_collectors = 2
   TMPIFile *newfile = new TMPIFile("mpi_output.root", "RECREATE", n_collectors);

   if (newfile->IsCollector())
      // collector rank, listens for incoming data
      newfile->RunCollector();
   else{
      // worker rank

      // generate data objects

      // syncronize data across ranks
      newfile->Sync();
   }

   newfile->Close();
}
End_Macro

See TMPIFile class for the list of functions
*/

const Int_t MIN_FILE_NUM = 2;

////////////////////////////////////////////////////////////////////////////////
/// TMPIFile constructor
///
/// See TMemFile for constructor explanation of the syntax.
///
/// \param[split] is the number of collectors to use

TMPIFile::TMPIFile(const char *name, char *buffer, Long64_t size, Option_t *option, Int_t split, const char *ftitle,
                   Int_t compress)
   : TMemFile(name, buffer, size, option, ftitle, compress), fSplitLevel(split), fMPIColor(0), fMPIRequest(0),
     fSendBuf(0)
{
   // check that split is set to reasonable value
   CheckSplitLevel();
   // split the MPI COMM_WORLD into sub-communicators
   // one for each split level
   SplitMPIComm();
}

////////////////////////////////////////////////////////////////////////////////
/// TMPIFile constructor
///
/// See TMemFile for constructor explanation of the syntax.
///
/// \param[split] is the number of collectors to use

TMPIFile::TMPIFile(const char *name, Option_t *option, Int_t split, const char *ftitle, Int_t compress)
   : TMemFile(name, option, ftitle, compress), fSplitLevel(split), fMPIColor(0), fMPIRequest(0), fSendBuf(0)
{
   // check that split is set to reasonable value
   CheckSplitLevel();
   // split the MPI COMM_WORLD into sub-communicators
   // one for each split level
   SplitMPIComm();
}

////////////////////////////////////////////////////////////////////////////////
/// TMPIFile destructor
///
/// This ensures the TMPIFile::Close function is called and any MPI
/// communicators are freed.

TMPIFile::~TMPIFile()
{
   // Sub communicators should be freed
   Int_t finalized = 0;
   MPI_Finalized(&finalized);
   if (!finalized && (fSplitLevel > 1)) {
      MPI_Comm_free(&fSubComm);
   }
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// As worker ranks exit, they send the collector empty messages.
/// This counter keeps track of the number of empty messages the collector
/// has received. Thereby the collector knows when all workers have exited
/// and it can exit

void TMPIFile::UpdateEndProcess()
{
   fEndProcess++;
}

////////////////////////////////////////////////////////////////////////////////
/// This is the core of the Collector rank which listens for incoming
/// messages from Worker ranks. The Collector

void TMPIFile::RunCollector(Bool_t cache)
{
   // update the user set filename with the current process ID and Rank ID
   // this ensures collectors do not overwrite one anothers files
   this->SetOutputName();
   Info("RunCollector", "writing to filename: %s", fMPIFilename.Data());
   THashTable mergers;

   Int_t client_Id = 0;
   std::vector<char> buffer(0);

   // loop until all other ranks in the subcommunicator have exited
   while (fEndProcess != fMPILocalSize - 1) {
      // Info("RunCollector","process counter %i",fEndProcess);
      // check if message has been received
      MPI_Status status;
      // this call blocks until a message is received
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, fSubComm, &status);

      // get bytes received
      Int_t number_bytes;
      MPI_Get_count(&status, MPI_CHAR, &number_bytes);
      buffer.resize(number_bytes);
      char *buf = buffer.data();

      Int_t source = status.MPI_SOURCE;
      Int_t tag = status.MPI_TAG;

      // retrieve the message
      MPI_Recv(buf, number_bytes, MPI_CHAR, source, tag, fSubComm, MPI_STATUS_IGNORE);

      // empty message signifies a Worker exited
      if (number_bytes == 0) {
         this->UpdateEndProcess();
      } else {
         // create a TMemFile from the buffer
         TMemFile *transient = new TMemFile(fMPIFilename, buf, number_bytes, "UPDATE");
         if (transient->IsZombie()) {
            Error("RunCollector", "Failed to create TMemFile from buffer");
         }
         // match compression settings of this TMPIFile object
         transient->SetCompressionSettings(this->GetCompressionSettings());

         // retrieve existing output file object
         ParallelFileMerger *info = (ParallelFileMerger *)mergers.FindObject(fMPIFilename);
         // if exiting file does not exist, create a new one
         if (!info) {
            info = new ParallelFileMerger(fMPIFilename, this->GetCompressionSettings(), cache);
            // add file to hash table
            mergers.Add(info);
         }

         // first merge needs extra care
         if (info->NeedInitialMerge(transient)) {
            info->InitialMerge(transient);
         }

         // merge the data
         info->RegisterClient(client_Id, transient);
         info->Merge();
         transient = 0;

         client_Id++;
      }
      buffer.resize(0);
   }

   if (fEndProcess == fMPILocalSize - 1) {
      mergers.Delete();
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor for ParallelFileMerger class

TMPIFile::ParallelFileMerger::ParallelFileMerger(const char *filename, Int_t compression_settings, Bool_t writeCache)
   : fFilename(filename), fClientsContact(0), fMerger(kFALSE, kTRUE)
{
   fMerger.SetPrintLevel(0);
   if (!fMerger.OutputFile(filename, "RECREATE")) {
      Error("ParallelFileMerger", "Cannot recreate the output file");
   }
   fMerger.GetOutputFile()->SetCompressionSettings(compression_settings);
   if (writeCache) {
      new TFileCacheWrite(fMerger.GetOutputFile(), 32 * 1024 * 1024);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Deconstructor for ParallelFileMerger class
///
/// taken from root/tutorials/net/parallelMergeServer.C

TMPIFile::ParallelFileMerger::~ParallelFileMerger()
{
   for (auto& client : fClients)
      delete client.GetFile();
}

void TMPIFile::ParallelFileMerger::DeleteObject(TDirectory *dir, Bool_t withReset)
{
   if (dir == 0)
      return;

   TIter nextkey(dir->GetListOfKeys());
   TKey *key;
   while ((key = (TKey *)nextkey())) {
      TClass *cl = TClass::GetClass(key->GetClassName());
      if (cl->InheritsFrom(TDirectory::Class())) {
         TDirectory *subdir = (TDirectory *)dir->GetList()->FindObject(key->GetName());
         if (!subdir) {
            subdir = (TDirectory *)key->ReadObj();
         }
         DeleteObject(subdir, withReset);
      } else {
         Bool_t todelete = kFALSE;
         if (withReset) {
            todelete = (0 != cl->GetResetAfterMerge());
         } else {
            todelete = (0 == cl->GetResetAfterMerge());
         }
         if (todelete) {
            key->Delete();
            dir->GetListOfKeys()->Remove(key);
            delete key;
         }
      }
   }
}

Bool_t TMPIFile::ParallelFileMerger::NeedInitialMerge(TDirectory *dir)
{
   if (dir == 0)
      return kFALSE;
   TIter nextkey(dir->GetListOfKeys());
   TKey *key;
   while ((key = (TKey *)nextkey())) {
      TClass *cl = TClass::GetClass(key->GetClassName());
      if (cl->InheritsFrom(TDirectory::Class())) {
         TDirectory *subdir = (TDirectory *)dir->GetList()->FindObject(key->GetName());
         if (!subdir) {
            subdir = (TDirectory *)key->ReadObj();
         }
         if (NeedInitialMerge(subdir)) {
            return kTRUE;
         }
      } else {
         if (0 != cl->GetResetAfterMerge()) {
            return kTRUE;
         }
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Initial merge of the input to copy the resetable object (TTree) into the output
/// and remove them from the input file.
///
/// taken from root/tutorials/net/parallelMergeServer.C

Bool_t TMPIFile::ParallelFileMerger::InitialMerge(TFile *input)
{
   // Initial merge of the input to copy the resetable object (TTree) into the
   // output and remove them from the input file.
   fMerger.AddFile(input);
   Bool_t result =
      fMerger.PartialMerge(TFileMerger::kIncremental | TFileMerger::kResetable | TFileMerger::kKeepCompression);
   DeleteObject(input, kTRUE);
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Merge the current inputs into the output file.
///
/// taken from root/tutorials/net/parallelMergeServer.C

Bool_t TMPIFile::ParallelFileMerger::Merge()
{
   // Remove object that can *not* be incrementally merge and will *not* be reset by the client code.
   DeleteObject(fMerger.GetOutputFile(), kFALSE);
   for (UInt_t f = 0; f < fClients.size(); ++f) {
      fMerger.AddFile(fClients[f].GetFile());
   }
   Bool_t result = fMerger.PartialMerge(TFileMerger::kAllIncremental | TFileMerger::kKeepCompression);

   // Remove any 'resetable' object (like TTree) from the input file so that they
   // will not be re-merged.  Keep only the object that always need to be
   // re-merged (Histograms).
   for (UInt_t f = 0; f < fClients.size(); ++f) {
      if (fClients[f].GetFile()) {
         DeleteObject(fClients[f].GetFile(), kTRUE);
      } else {
         // We back up the file (probably due to memory constraint)
         TFile *file = TFile::Open(fClients[f].GetLocalName(), "UPDATE");
         if (file->IsZombie()) {
            Error("Merge", "output file unavailable");
         }
         // Remove object that can be incrementally merge and will be reset by the client code.
         DeleteObject(file, kTRUE);
         file->Write();
         delete file;
      }
   }
   fLastMerge = TTimeStamp();
   fNClientsContact = 0;
   fClientsContact.Clear();

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Register that a client has sent a file.
///
/// taken from root/tutorials/net/parallelMergeServer.C

void TMPIFile::ParallelFileMerger::RegisterClient(UInt_t clientID, TFile *file)
{
   ++fNClientsContact;
   fClientsContact.SetBitNumber(clientID);
   TMPIClientInfo ntcl(std::string(fFilename).c_str(), clientID);
   if (fClients.size() < clientID + 1) {
      fClients.push_back(ntcl);
   }
   fClients[clientID].SetFile(file);
}

////////////////////////////////////////////////////////////////////////////////
/// Return true, if enough client have reported
///
/// In the case of TMPIFile this happens everytime a client/worker sends
/// the buffer (tested).
///
/// taken from root/tutorials/net/parallelMergeServer.C

Bool_t TMPIFile::ParallelFileMerger::NeedMerge(Float_t clientThreshold)
{

   if (fClients.size() == 0) {
      return kFALSE;
   }

   // Calculate average and rms of the time between the last 2 contacts.
   Double_t sum = 0.;
   Double_t sum2 = 0.;
   for (UInt_t c = 0; c < fClients.size(); ++c) {
      sum += fClients[c].GetTimeSincePrevContact();
      sum2 += fClients[c].GetTimeSincePrevContact() * fClients[c].GetTimeSincePrevContact();
   }
   Double_t avg = sum / fClients.size();
   Double_t sigma = sum2 ? TMath::Sqrt(sum2 / fClients.size() - avg * avg) : 0;
   Double_t target = avg + 2 * sigma;
   TTimeStamp now;
   if ((now.AsDouble() - fLastMerge.AsDouble()) > target) {
      return kTRUE;
   }
   Float_t cut = clientThreshold * fClients.size();
   return fClientsContact.CountBits() > cut || fNClientsContact > 2 * cut;
}



////////////////////////////////////////////////////////////////////////////////
/// return True if this is the Collector rank, otherwise False

Bool_t TMPIFile::IsCollector()
{
   return !fMPILocalRank;
}

////////////////////////////////////////////////////////////////////////////////
/// Called by the Workers only: Copies the current content in memory and
/// sends it asynchronously to the Collector for merging and writing to disk.
///

void TMPIFile::CreateBufferAndSend()
{
   if (this->IsCollector()) {
      Error("CreateBufferAndSend", " should not be called by a collector");
      return;
   }
   this->Write();
   Int_t count = this->GetEND();
   fSendBuf = new char[count];
   this->CopyTo(fSendBuf, count);
   MPI_Isend(fSendBuf, count, MPI_CHAR, 0, fMPIColor, fSubComm, &fMPIRequest);
}

////////////////////////////////////////////////////////////////////////////////
/// For Workers: Creates an empty buffer and sends it to the Collector.
/// This indicates the completion of the worker.

void TMPIFile::CreateEmptyBufferAndSend()
{
   if (this->IsCollector()) {
      return;
   }

   if (!IsReceived()) {
      MPI_Wait(&fMPIRequest, MPI_STATUS_IGNORE);
   }
   delete[] fSendBuf; // empty the buffer once received by master
   fSendBuf = nullptr;
   MPI_Send(fSendBuf, 0, MPI_CHAR, 0, fMPIColor, fSubComm);
}

////////////////////////////////////////////////////////////////////////////////
/// Called by the Workers only: Called periodically by workers and triggers
/// the sending of data to the Collector for writing.

void TMPIFile::Sync()
{
   // check if the previous send request is accepted by master.
   if (!IsReceived()) {
      MPI_Wait(&fMPIRequest, MPI_STATUS_IGNORE);
   }
   delete[] fSendBuf; // empty the buffer once received by master
   fSendBuf = nullptr;
   CreateBufferAndSend();
   this->ResetAfterMerge((TFileMergeInfo *)0);
}

////////////////////////////////////////////////////////////////////////////////
/// Closes the file. For Worker ranks, this function will signal to the
/// Collector that the Worker has exited. It also closes the inherited TMemFile.

void TMPIFile::Close(Option_t *option)
{
   if (IsOpen()) {
      // sends empty buffer
      CreateEmptyBufferAndSend();
      // call parent close function
      TMemFile::Close(option);
      
      // check to see that MPI has not already been finalized
      Int_t finalized = 0;
      MPI_Finalized(&finalized);
      if (!finalized) {
         MPI_Finalize();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called by the Collector only: edits the input filename from the user
/// to append the rank ID of the Collector so that each collector has a
/// unique filename.

void TMPIFile::SetOutputName()
{
   std::string _filename = this->GetName();

   ULong_t found = _filename.rfind(".root");
   if (found != std::string::npos) {
      _filename.resize(found);
   }
   fMPIFilename = _filename;
   fMPIFilename += "_";
   fMPIFilename += fMPIColor;
   fMPIFilename += ".root";
}

////////////////////////////////////////////////////////////////////////////////
/// Checks that the split level is more than one.
/// There must be at least one Worker and one Collector rank.

void TMPIFile::CheckSplitLevel()
{
   if (fSplitLevel < 1) {
      Error("CheckSplitLevel", "At least one collector is required instead of %d", fSplitLevel);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called by all ranks to create the sub communicators (if more than one
/// rank).
///

void TMPIFile::SplitMPIComm()
{
   // Initialize MPI if it is not already initialized...
   Int_t flag;
   MPI_Initialized(&flag);
   if (!flag) {
      MPI_Init(NULL, NULL);
   }
   // get global size and current global rank
   MPI_Comm_size(MPI_COMM_WORLD, &fMPIGlobalSize);
   MPI_Comm_rank(MPI_COMM_WORLD, &fMPIGlobalRank);

   if (MIN_FILE_NUM * fSplitLevel > fMPIGlobalSize) {
      Error("TMPIFile",
            "Number of Output File is larger than number of Processors Allocated."
            " Number of processors should be two times larger than outpts. For %d outputs at least %d "
            "should be allocated instead of %d",
            fSplitLevel, MIN_FILE_NUM * fSplitLevel, fMPIGlobalSize);
   }

   // using one collector
   if (fSplitLevel == 1) {
      fSubComm = MPI_COMM_WORLD;
   }
   // using more than one collector
   else {
      // number of ranks per sub-communicator
      Int_t comm_size = fMPIGlobalSize / fSplitLevel;
      if (fMPIGlobalSize % fSplitLevel != 0) {
         comm_size++;
      }
      fMPIColor = fMPIGlobalRank / comm_size;
      // split the COMM_WORLD communicator by color
      MPI_Comm_split(MPI_COMM_WORLD, fMPIColor, fMPIGlobalRank, &fSubComm);
   }
   // get the sub-communicator size and rank
   MPI_Comm_size(fSubComm, &fMPILocalSize);
   MPI_Comm_rank(fSubComm, &fMPILocalRank);
}

////////////////////////////////////////////////////////////////////////////////
/// Checks the member MPI_REQEUST object to see if a message has been received.

Bool_t TMPIFile::IsReceived()
{
   if (!fMPIRequest) {
      return kTRUE;
   }
   Int_t flag = 0;
   MPI_Test(&fMPIRequest, &flag, MPI_STATUS_IGNORE);
   return flag;
}
