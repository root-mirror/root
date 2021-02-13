// @(#)root/tree:$Id$
// Author: Anna Kreshuk 17/03/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TEntryListFromFile
\ingroup tree

Manages entry lists from different files, when they are not loaded
in memory at the same time.

This entry list should only be used when processing a TChain (see
TChain::SetEntryList() function). File naming convention:
- by default, filename_elist.root is used, where filename is the
  name of the chain element.
- xxx$xxx.root - $ sign is replaced by the name of the chain element
If the list name is not specified (by passing filename_elist.root/listname to
the TChain::SetEntryList() function, the first object of class TEntryList
in the file is taken.
It is assumed that there are as many lists, as there are chain elements,
and they are in the same order.

If one of the list files can't be opened, or there is an error reading a list
from the file, this list is skipped and the entry loop continues on the next
list.
*/

#include "TEntryListFromFile.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TKey.h"
#include "TError.h"
#include "TTree.h"

ClassImp(TEntryListFromFile);

TEntryListFromFile::TEntryListFromFile(): TEntryList(),
   fListFileName(""), fListName(""), fNFiles(0), fListOffset(0), fFile(0), fFileNames(0)
{
   // default constructor.

}

////////////////////////////////////////////////////////////////////////////////
/// File naming convention:
/// - by default, filename_elist.root is used, where filename is the
///   name of the chain element
/// - xxx$xxx.root - $ sign is replaced by the name of the chain element
///
/// The TObjArray of chain elements is set by the TEntryListFromFile::SetFileNames()
/// function.
///
/// If the list name is not specified, the first object of class TEntryList
/// in the file is taken.
///
/// nfiles is the total number of files to process

TEntryListFromFile::TEntryListFromFile(const char *filename, const char *listname, Int_t nfiles) : TEntryList(),
   fListFileName(filename), fListName(listname), fNFiles(nfiles), fListOffset(0), fFile(0), fFileNames(0)
{
   fListOffset = new Long64_t[fNFiles+1];
   fListOffset[0]=0;
   for (Int_t i=1; i<fNFiles+1; i++){
      fListOffset[i]=TTree::kMaxEntries;
   }
   fN = TTree::kMaxEntries;
}

////////////////////////////////////////////////////////////////////////////////
/// d-tor

TEntryListFromFile::~TEntryListFromFile()
{
   delete [] fListOffset;
   fListOffset = 0;
   delete fFile;
   fFile = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns entry \#index
/// See also Next() for a faster alternative

Long64_t TEntryListFromFile::GetEntry(Int_t index)
{
   if (index<0) return -1;

   if (index > fListOffset[fNFiles] && fListOffset[fNFiles]!=TTree::kMaxEntries){
      Error("GetEntry", "Index value is too large\n");
      return -1;
   }

   if (index==fLastIndexQueried+1)
      return Next();

   Int_t itree =0;
   while (!fCurrent && itree<fNFiles){
      LoadList(itree);
      itree++;
   }
   if (itree == fNFiles){
      Error("GetEntry", "All lists are empty\n");
      return -1;
   }

   if (index < fListOffset[fTreeNumber]) {
      //this entry is in one of previously opened lists
      itree=0;
      for (itree=0; itree<fTreeNumber; itree++){
         if (index >= fListOffset[itree] && fListOffset[itree]!=fListOffset[itree+1])
            break;
      }
      LoadList(itree);
   }
   else if (index >= fListOffset[fTreeNumber+1]){
      //this entry is in one of following lists
      itree = fTreeNumber;
      while (itree < fNFiles){
         itree++;
         if (fListOffset[itree+1]==TTree::kMaxEntries){
            //this list hasn't been loaded yet
            LoadList(itree);
         }
         if (index < fListOffset[itree+1]){
            //the entry is in this list
            break;
         }
      }
      if (fTreeNumber == fNFiles){
         Error("GetEntry", "Entry number is too big\n");
         return -1;
      }
      if (fTreeNumber!=itree)
         LoadList(itree);
   }
   //now the entry is in the currently opened list
   Long64_t localentry = index - fListOffset[fTreeNumber];
   Long64_t retentry = fCurrent->GetEntry(localentry);
   fLastIndexQueried = index;
   fLastIndexReturned = retentry;
   return retentry;

}

////////////////////////////////////////////////////////////////////////////////
/// Return the entry corresponding to the index parameter and the
/// number of the tree, where this entry is

Long64_t TEntryListFromFile::GetEntryAndTree(Int_t index, Int_t &treenum)
{
   Long64_t result = GetEntry(index);
   treenum = fTreeNumber;
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the total number of entries in the list.
/// If some lists have not been loaded, loads them.

Long64_t TEntryListFromFile::GetEntries()
{
   if (fN==TTree::kMaxEntries){
      for (Int_t i=0; i<fNFiles; i++){
         if (fListOffset[i+1]==TTree::kMaxEntries){
            LoadList(i);
         }
      }
   }
   fN = fListOffset[fNFiles];
   fLastIndexQueried = -3;
   return fN;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the next entry in the list.
/// Faster than GetEntry()

Long64_t TEntryListFromFile::Next()
{
   Int_t itree =0;
   while (!fCurrent && itree<fNFiles){
      LoadList(itree);
      itree++;
   }
   if (itree == fNFiles){
      Error("Next", "All lists are empty\n");
      return -1;
   }

   Long64_t retentry = fCurrent->Next();
   if (retentry<0){
      if (fLastIndexQueried == fListOffset[fTreeNumber+1]-1){
         //requested entry is in the next list
         if (fTreeNumber == fNFiles -1){
            // Error("Next", "No more entries, last list\n");
            return -1;
         }
         do{
            //load the next non-empty list. fTreeNumber is changed by LoadList()
            fTreeNumber++;
            LoadList(fTreeNumber);
         } while (fListOffset[fTreeNumber+1]==fListOffset[fTreeNumber] && fTreeNumber<fNFiles-1);
         if (fTreeNumber == fNFiles -1 && fListOffset[fTreeNumber+1]==fListOffset[fTreeNumber]){
            //no more lists
            return -1;
         }
         retentry = fCurrent->Next();
      } else {
         Error("Next", "Something wrong with reading the current list, even though thefile #%d and the list exist\n", fTreeNumber);
         return -1;
      }

   }

   fLastIndexQueried++;
   fLastIndexReturned = retentry;
   return retentry;

}

////////////////////////////////////////////////////////////////////////////////
/// Loads the list \#listnumber
/// This is the only function that can modify fCurrent and fFile data members

Int_t TEntryListFromFile::LoadList(Int_t listnumber)
{
   //first close the current list
   if (fCurrent){
      if (fFile) {
         delete fFile;
         fFile = 0;
         fCurrent = 0;
      }
   }

   R__ASSERT(fFileNames);

   //find the right name
   //get the name of the corresponding chain element (with the treenumber=listnumber)
   TNamed *nametitle = (TNamed*)fFileNames->At(listnumber);
   TString filename_short = nametitle->GetTitle();
   if (filename_short.Contains(".root")){
      filename_short.Remove(filename_short.Length()-5, 5);
   }
   if (!strcmp(fListFileName.Data(), "")){
      //no name supplied, use the one of the chain file
      filename_short.Append("_elist.root");
      //printf("filename: %s\n", filename_short.Data());
      fFile = TFile::Open(filename_short.Data());
   } else {
      TString filename = fListFileName;
      filename.ReplaceAll("$", filename_short);
      //printf("filename: %s\n", filename.Data());
      fFile = TFile::Open(filename.Data());
   }

   if (!fFile || fFile->IsZombie()){
      if (fFile) {
         delete fFile;
         fFile = 0;
      }
      fCurrent = 0;
      fListOffset[listnumber+1] = fListOffset[listnumber];
      return -1;
   }

   if (!strcmp(fListName.Data(), "")){
      TKey *key;
      TIter nextkey(fFile->GetListOfKeys());
      while ((key=(TKey*)nextkey())){
         if (strcmp("TEntryList", key->GetClassName())==0){
            //found an object of class TEntryList
            fCurrent = (TEntryList*)key->ReadObj();
         }
      }
   } else {
      fCurrent = (TEntryList*)fFile->Get(fListName.Data());
   }

   if (!fCurrent){
      Error("LoadList", "List %s not found in the file\n", fListName.Data());
      fCurrent = 0;
      fListOffset[listnumber+1]=fListOffset[listnumber];
      return -1;
   }
   fTreeNumber = listnumber;
   Long64_t nentries = fCurrent->GetN();
   if (fListOffset[fTreeNumber+1] != (fListOffset[fTreeNumber] + nentries)) {
      fListOffset[fTreeNumber+1] = fListOffset[fTreeNumber] + nentries;
      fN = fListOffset[fNFiles];
   }

   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Print info about this list

void TEntryListFromFile::Print(const Option_t* option) const
{
   printf("total number of files: %d\n", fNFiles);
   TFile *f;
   TEntryList *el=0;
   if (fFileNames==0) {
      Error("Print","fFileNames was not set properly.");
   } else {
      for (Int_t listnumber=0; listnumber<fNFiles; listnumber++){
         TNamed *nametitle = (TNamed*)fFileNames->At(listnumber);
         TString filename_short = nametitle->GetTitle();
         if (filename_short.Contains(".root")){
            filename_short.Remove(filename_short.Length()-5, 5);
         }
         if (!strcmp(fListFileName.Data(), "")){
            //no name supplied, use the one of the chain file
            filename_short.Append("_elist.root");
            //printf("filename: %s\n", filename_short.Data());
            f = TFile::Open(filename_short.Data());
         } else {
            TString filename = fListFileName;
            filename.ReplaceAll("$", filename_short);
            //printf("filename: %s\n", filename.Data());
            f = TFile::Open(filename.Data());
         }
         if (f && !f->IsZombie()){
            if (!strcmp(fListName.Data(), "")){
               TKey *key;
               TIter nextkey(f->GetListOfKeys());
               while ((key=(TKey*)nextkey())){
                  if (strcmp("TEntryList", key->GetClassName())==0){
                     //found an object of class TEntryList
                     el = (TEntryList*)key->ReadObj();
                  }
               }
            } else {
               el = (TEntryList*)f->Get(fListName.Data());
            }
            if (el)
               el->Print(option);
         }
      }
   }

}
