// @(#)root/tree:$Id$
// Author: Rene Brun   03/02/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TChain
\ingroup tree

A chain is a collection of files containing TTree objects.
When the chain is created, the first parameter is the default name
for the Tree to be processed later on.

Enter a new element in the chain via the TChain::Add function.
Once a chain is defined, one can use the normal TTree functions
to Draw,Scan,etc.

Use TChain::SetBranchStatus to activate one or more branches for all
the trees in the chain.
*/

#include "TChain.h"

#include <iostream>
#include <cfloat>

#include "TBranch.h"
#include "TBrowser.h"
#include "TBuffer.h"
#include "TChainElement.h"
#include "TClass.h"
#include "TColor.h"
#include "TCut.h"
#include "TError.h"
#include "TFile.h"
#include "TFileInfo.h"
#include "TFriendElement.h"
#include "TLeaf.h"
#include "TList.h"
#include "TObjString.h"
#include "TPluginManager.h"
#include "TROOT.h"
#include "TRegexp.h"
#include "TSelector.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeCache.h"
#include "TUrl.h"
#include "TVirtualIndex.h"
#include "TEventList.h"
#include "TEntryList.h"
#include "TEntryListFromFile.h"
#include "TFileStager.h"
#include "TFilePrefetch.h"
#include "TVirtualMutex.h"
#include "TVirtualPerfStats.h"
#include "strlcpy.h"
#include "snprintf.h"

ClassImp(TChain);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TChain::TChain()
: TTree()
, fTreeOffsetLen(100)
, fNtrees(0)
, fTreeNumber(-1)
, fTreeOffset(0)
, fCanDeleteRefs(kFALSE)
, fTree(0)
, fFile(0)
, fFiles(0)
, fStatus(0)
, fProofChain(0)
{
   fTreeOffset = new Long64_t[fTreeOffsetLen];
   fFiles = new TObjArray(fTreeOffsetLen);
   fStatus = new TList();
   fTreeOffset[0]  = 0;
   if (gDirectory) gDirectory->Remove(this);
   gROOT->GetListOfSpecials()->Add(this);
   fFile = 0;
   fDirectory = 0;

   // Reset PROOF-related bits
   ResetBit(kProofUptodate);
   ResetBit(kProofLite);

   // Add to the global list
   gROOT->GetListOfDataSets()->Add(this);

   // Make sure we are informed if the TFile is deleted.
   R__LOCKGUARD(gROOTMutex);
   gROOT->GetListOfCleanups()->Add(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a chain.
///
/// A TChain is a collection of TFile objects.
/// the first parameter "name" is the name of the TTree object
/// in the files added with Add.
/// Use TChain::Add to add a new element to this chain.
///
/// In case the Tree is in a subdirectory, do, eg:
/// ~~~ {.cpp}
///     TChain ch("subdir/treename");
/// ~~~
/// Example:
///  Suppose we have 3 files f1.root, f2.root and f3.root. Each file
///  contains a TTree object named "T".
/// ~~~ {.cpp}
///     TChain ch("T");  creates a chain to process a Tree called "T"
///     ch.Add("f1.root");
///     ch.Add("f2.root");
///     ch.Add("f3.root");
///     ch.Draw("x");
/// ~~~
/// The Draw function above will process the variable "x" in Tree "T"
/// reading sequentially the 3 files in the chain ch.
///
/// The TChain data structure:
///
/// Each TChainElement has a name equal to the tree name of this TChain
/// and a title equal to the file name. So, to loop over the
/// TFiles that have been added to this chain:
/// ~~~ {.cpp}
///     TObjArray *fileElements=chain->GetListOfFiles();
///     TIter next(fileElements);
///     TChainElement *chEl=0;
///     while (( chEl=(TChainElement*)next() )) {
///        TFile f(chEl->GetTitle());
///        ... do something with f ...
///     }
/// ~~~

TChain::TChain(const char* name, const char* title)
:TTree(name, title, /*splitlevel*/ 99, nullptr)
, fTreeOffsetLen(100)
, fNtrees(0)
, fTreeNumber(-1)
, fTreeOffset(0)
, fCanDeleteRefs(kFALSE)
, fTree(0)
, fFile(0)
, fFiles(0)
, fStatus(0)
, fProofChain(0)
{
   //
   //*-*

   fTreeOffset = new Long64_t[fTreeOffsetLen];
   fFiles = new TObjArray(fTreeOffsetLen);
   fStatus = new TList();
   fTreeOffset[0]  = 0;
   fFile = 0;

   // Reset PROOF-related bits
   ResetBit(kProofUptodate);
   ResetBit(kProofLite);

   R__LOCKGUARD(gROOTMutex);

   // Add to the global lists
   gROOT->GetListOfSpecials()->Add(this);
   gROOT->GetListOfDataSets()->Add(this);

   // Make sure we are informed if the TFile is deleted.
   gROOT->GetListOfCleanups()->Add(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TChain::~TChain()
{
   bool rootAlive = gROOT && !gROOT->TestBit(TObject::kInvalidObject);

   if (rootAlive) {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfCleanups()->Remove(this);
   }

   SafeDelete(fProofChain);
   fStatus->Delete();
   delete fStatus;
   fStatus = 0;
   fFiles->Delete();
   delete fFiles;
   fFiles = 0;

   //first delete cache if exists
   auto tc = fFile && fTree ? fTree->GetReadCache(fFile) : nullptr;
   if (tc) {
      delete tc;
      fFile->SetCacheRead(0, fTree);
   }

   delete fFile;
   fFile = 0;
   // Note: We do *not* own the tree.
   fTree = 0;
   delete[] fTreeOffset;
   fTreeOffset = 0;

   // Remove from the global lists
   if (rootAlive) {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfSpecials()->Remove(this);
      gROOT->GetListOfDataSets()->Remove(this);
   }

   // This is the same as fFile, don't delete it a second time.
   fDirectory = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add all files referenced by the passed chain to this chain.
/// The function returns the total number of files connected.

Int_t TChain::Add(TChain* chain)
{
   if (!chain) return 0;

   // Check for enough space in fTreeOffset.
   if ((fNtrees + chain->GetNtrees()) >= fTreeOffsetLen) {
      fTreeOffsetLen += 2 * chain->GetNtrees();
      Long64_t* trees = new Long64_t[fTreeOffsetLen];
      for (Int_t i = 0; i <= fNtrees; i++) {
         trees[i] = fTreeOffset[i];
      }
      delete[] fTreeOffset;
      fTreeOffset = trees;
   }
   chain->GetEntries(); //to force the computation of nentries
   TIter next(chain->GetListOfFiles());
   Int_t nf = 0;
   TChainElement* element = 0;
   while ((element = (TChainElement*) next())) {
      Long64_t nentries = element->GetEntries();
      if (fTreeOffset[fNtrees] == TTree::kMaxEntries) {
         fTreeOffset[fNtrees+1] = TTree::kMaxEntries;
      } else {
         fTreeOffset[fNtrees+1] = fTreeOffset[fNtrees] + nentries;
      }
      fNtrees++;
      fEntries += nentries;
      TChainElement* newelement = new TChainElement(element->GetName(), element->GetTitle());
      newelement->SetPacketSize(element->GetPacketSize());
      newelement->SetNumberEntries(nentries);
      fFiles->Add(newelement);
      nf++;
   }
   if (fProofChain)
      // This updates the proxy chain when we will really use PROOF
      ResetBit(kProofUptodate);

   return nf;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new file to this chain.
///
/// Argument name may have either of two set of formats. The first:
/// ~~~ {.cpp}
///     [//machine]/path/file_name[?query[#tree_name]]
///  or [//machine]/path/file_name.root[.oext][/tree_name]
/// ~~~
/// If tree_name is missing the chain name will be assumed. Tagging the
/// tree_name with a slash [/tree_name] is only supported for backward
/// compatibility; it requires the file name to contain the string '.root'
/// and its use is deprecated.
/// Wildcard treatment is triggered by the any of the special characters []*?
/// which may be used in the file name, eg. specifying "xxx*.root" adds
/// all files starting with xxx in the current file system directory.
///
/// Alternatively name may have the format of a url, eg.
/// ~~~ {.cpp}
///         root://machine/path/file_name[?query[#tree_name]]
///     or  root://machine/path/file_name
///     or  root://machine/path/file_name.root[.oext]/tree_name
///     or  root://machine/path/file_name.root[.oext]/tree_name?query
/// ~~~
/// where "query" is to be interpreted by the remote server. Wildcards may be
/// supported in urls, depending on the protocol plugin and the remote server.
/// http or https urls can contain a query identifier without tree_name, but
/// generally urls can not be written with them because of ambiguity with the
/// wildcard character. (Also see the documentation for TChain::AddFile,
/// which does not support wildcards but allows the url name to contain query).
/// Again, tagging the tree_name with a slash [/tree_name] is only supported
/// for backward compatibility; it requires the file name ot contain the string
/// '.root' and its use is deprecated.
///
/// NB. To add all the files of a TChain to a chain, use Add(TChain *chain).
///
/// A. if nentries <= 0, the file is connected and the tree header read
///    in memory to get the number of entries.
///
/// B. if (nentries > 0, the file is not connected, nentries is assumed to be
///    the number of entries in the file. In this case, no check is made that
///    the file exists and the Tree existing in the file. This second mode
///    is interesting in case the number of entries in the file is already stored
///    in a run data base for example.
///
/// C. if (nentries == TTree::kMaxEntries) (default), the file is not connected.
///    the number of entries in each file will be read only when the file
///    will need to be connected to read an entry.
///    This option is the default and very efficient if one process
///    the chain sequentially. Note that in case TChain::GetEntry(entry)
///    is called and entry refers to an entry in the 3rd file, for example,
///    this forces the Tree headers in the first and second file
///    to be read to find the number of entries in these files.
///    Note that if one calls TChain::GetEntriesFast() after having created
///    a chain with this default, GetEntriesFast will return TTree::kMaxEntries!
///    TChain::GetEntries will force of the Tree headers in the chain to be
///    read to read the number of entries in each Tree.
///
/// D. The TChain data structure
///    Each TChainElement has a name equal to the tree name of this TChain
///    and a title equal to the file name. So, to loop over the
///    TFiles that have been added to this chain:
/// ~~~ {.cpp}
///        TObjArray *fileElements=chain->GetListOfFiles();
///        TIter next(fileElements);
///        TChainElement *chEl=0;
///        while (( chEl=(TChainElement*)next() )) {
///           TFile f(chEl->GetTitle());
///           ... do something with f ...
///        }
/// ~~~
/// Return value:
///
/// - If nentries>0 (including the default of TTree::kMaxEntries) and no
///   wildcarding is used, ALWAYS returns 1 without regard to whether
///   the file exists or contains the correct tree.
///
/// - If wildcarding is used, regardless of the value of nentries,
///   returns the number of files matching the name without regard to
///   whether they contain the correct tree.
///
/// - If nentries<=0 and wildcarding is not used, return 1 if the file
///  exists and contains the correct tree and 0 otherwise.

Int_t TChain::Add(const char* name, Long64_t nentries /* = TTree::kMaxEntries */)
{
   TString basename, treename, query, suffix;
   ParseTreeFilename(name, basename, treename, query, suffix, kTRUE);

   // case with one single file
   if (!basename.MaybeWildcard()) {
      return AddFile(name, nentries);
   }

   // wildcarding used in name
   Int_t nf = 0;

   Int_t slashpos = basename.Last('/');
   TString directory;
   if (slashpos>=0) {
      directory = basename(0,slashpos); // Copy the directory name
      basename.Remove(0,slashpos+1);      // and remove it from basename
   } else {
      directory = gSystem->UnixPathName(gSystem->WorkingDirectory());
   }

   const char *file;
   const char *epath = gSystem->ExpandPathName(directory.Data());
   void *dir = gSystem->OpenDirectory(epath);
   delete [] epath;
   if (dir) {
      //create a TList to store the file names (not yet sorted)
      TList l;
      TRegexp re(basename,kTRUE);
      while ((file = gSystem->GetDirEntry(dir))) {
         if (!strcmp(file,".") || !strcmp(file,"..")) continue;
         TString s = file;
         if ( (basename!=file) && s.Index(re) == kNPOS) continue;
         l.Add(new TObjString(file));
      }
      gSystem->FreeDirectory(dir);
      //sort the files in alphanumeric order
      l.Sort();
      TIter next(&l);
      TObjString *obj;
      while ((obj = (TObjString*)next())) {
         file = obj->GetName();
         nf += AddFile(TString::Format("%s/%s%s",directory.Data(),file,suffix.Data()),nentries);
      }
      l.Delete();
   }
   if (fProofChain)
      // This updates the proxy chain when we will really use PROOF
      ResetBit(kProofUptodate);

   return nf;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new file to this chain.
///
/// Filename formats are similar to TChain::Add. Wildcards are not
/// applied. urls may also contain query and fragment identifiers
/// where the tree name can be specified in the url fragment.
///
/// eg.
/// ~~~ {.cpp}
///     root://machine/path/file_name[?query[#tree_name]]
///     root://machine/path/file_name.root[.oext]/tree_name[?query]
/// ~~~
/// If tree_name is given as a part of the file name it is used to
/// as the name of the tree to load from the file. Otherwise if tname
/// argument is specified the chain will load the tree named tname from
/// the file, otherwise the original treename specified in the TChain
/// constructor will be used.
/// Tagging the tree_name with a slash [/tree_name] is only supported for
/// backward compatibility; it requires the file name ot contain the string
/// '.root' and its use is deprecated.
///
/// A. If nentries <= 0, the file is opened and the tree header read
///    into memory to get the number of entries.
///
/// B. If nentries > 0, the file is not opened, and nentries is assumed
///    to be the number of entries in the file. In this case, no check
///    is made that the file exists nor that the tree exists in the file.
///    This second mode is interesting in case the number of entries in
///    the file is already stored in a run database for example.
///
/// C. If nentries == TTree::kMaxEntries (default), the file is not opened.
///    The number of entries in each file will be read only when the file
///    is opened to read an entry.  This option is the default and very
///    efficient if one processes the chain sequentially.  Note that in
///    case GetEntry(entry) is called and entry refers to an entry in the
///    third file, for example, this forces the tree headers in the first
///    and second file to be read to find the number of entries in those
///    files.  Note that if one calls GetEntriesFast() after having created
///    a chain with this default, GetEntriesFast() will return TTree::kMaxEntries!
///    Using the GetEntries() function instead will force all of the tree
///    headers in the chain to be read to read the number of entries in
///    each tree.
///
/// D. The TChain data structure
///    Each TChainElement has a name equal to the tree name of this TChain
///    and a title equal to the file name. So, to loop over the
///    TFiles that have been added to this chain:
/// ~~~ {.cpp}
///         TObjArray *fileElements=chain->GetListOfFiles();
///         TIter next(fileElements);
///         TChainElement *chEl=0;
///         while (( chEl=(TChainElement*)next() )) {
///            TFile f(chEl->GetTitle());
///            ... do something with f ...
///         }
/// ~~~
/// The function returns 1 if the file is successfully connected, 0 otherwise.

Int_t TChain::AddFile(const char* name, Long64_t nentries /* = TTree::kMaxEntries */, const char* tname /* = "" */)
{
   if(name==0 || name[0]=='\0') {
      Error("AddFile", "No file name; no files connected");
      return 0;
   }

   const char *treename = GetName();
   if (tname && strlen(tname) > 0) treename = tname;

   TString basename, tn, query, suffix;
   ParseTreeFilename(name, basename, tn, query, suffix, kFALSE);

   if (!tn.IsNull()) {
      treename = tn.Data();
   }

   Int_t nch = basename.Length() + query.Length();
   char *filename = new char[nch+1];
   strlcpy(filename,basename.Data(),nch+1);
   strlcat(filename,query.Data(),nch+1);

   //Check enough space in fTreeOffset
   if (fNtrees+1 >= fTreeOffsetLen) {
      fTreeOffsetLen *= 2;
      Long64_t *trees = new Long64_t[fTreeOffsetLen];
      for (Int_t i=0;i<=fNtrees;i++) trees[i] = fTreeOffset[i];
      delete [] fTreeOffset;
      fTreeOffset = trees;
   }

   // Open the file to get the number of entries.
   Int_t pksize = 0;
   if (nentries <= 0) {
      TFile* file;
      {
         TDirectory::TContext ctxt;
         file = TFile::Open(filename);
      }
      if (!file || file->IsZombie()) {
         delete file;
         file = 0;
         delete[] filename;
         filename = 0;
         return 0;
      }

      // Check that tree with the right name exists in the file.
      // Note: We are not the owner of obj, the file is!
      TObject* obj = file->Get(treename);
      if (!obj || !obj->InheritsFrom(TTree::Class())) {
         Error("AddFile", "cannot find tree with name %s in file %s", treename, filename);
         delete file;
         file = 0;
         delete[] filename;
         filename = 0;
         return 0;
      }
      TTree* tree = (TTree*) obj;
      nentries = tree->GetEntries();
      pksize = tree->GetPacketSize();
      // Note: This deletes the tree we fetched.
      delete file;
      file = 0;
   }

   if (nentries > 0) {
      if (nentries != TTree::kMaxEntries) {
         fTreeOffset[fNtrees+1] = fTreeOffset[fNtrees] + nentries;
         fEntries += nentries;
      } else {
         fTreeOffset[fNtrees+1] = TTree::kMaxEntries;
         fEntries = TTree::kMaxEntries;
      }
      fNtrees++;

      TChainElement* element = new TChainElement(treename, filename);
      element->SetPacketSize(pksize);
      element->SetNumberEntries(nentries);
      fFiles->Add(element);
   } else {
      Warning("AddFile", "Adding tree with no entries from file: %s", filename);
   }

   delete [] filename;
   if (fProofChain)
      // This updates the proxy chain when we will really use PROOF
      ResetBit(kProofUptodate);

   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Add all files referenced in the list to the chain. The object type in the
/// list must be either TFileInfo or TObjString or TUrl .
/// The function return 1 if successful, 0 otherwise.

Int_t TChain::AddFileInfoList(TCollection* filelist, Long64_t nfiles /* = TTree::kMaxEntries */)
{
   if (!filelist)
      return 0;
   TIter next(filelist);

   TObject *o = 0;
   Long64_t cnt=0;
   while ((o = next())) {
      // Get the url
      TString cn = o->ClassName();
      const char *url = 0;
      if (cn == "TFileInfo") {
         TFileInfo *fi = (TFileInfo *)o;
         url = (fi->GetCurrentUrl()) ? fi->GetCurrentUrl()->GetUrl() : 0;
         if (!url) {
            Warning("AddFileInfoList", "found TFileInfo with empty Url - ignoring");
            continue;
         }
      } else if (cn == "TUrl") {
         url = ((TUrl*)o)->GetUrl();
      } else if (cn == "TObjString") {
         url = ((TObjString*)o)->GetName();
      }
      if (!url) {
         Warning("AddFileInfoList", "object is of type %s : expecting TFileInfo, TUrl"
                                 " or TObjString - ignoring", o->ClassName());
         continue;
      }
      // Good entry
      cnt++;
      AddFile(url);
      if (cnt >= nfiles)
         break;
   }
   if (fProofChain) {
      // This updates the proxy chain when we will really use PROOF
      ResetBit(kProofUptodate);
   }

   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a TFriendElement to the list of friends of this chain.
///
/// A TChain has a list of friends similar to a tree (see TTree::AddFriend).
/// You can add a friend to a chain with the TChain::AddFriend method, and you
/// can retrieve the list of friends with TChain::GetListOfFriends.
/// This example has four chains each has 20 ROOT trees from 20 ROOT files.
/// ~~~ {.cpp}
///     TChain ch("t"); // a chain with 20 trees from 20 files
///     TChain ch1("t1");
///     TChain ch2("t2");
///     TChain ch3("t3");
/// ~~~
/// Now we can add the friends to the first chain.
/// ~~~ {.cpp}
///     ch.AddFriend("t1")
///     ch.AddFriend("t2")
///     ch.AddFriend("t3")
/// ~~~
/// \image html tchain_friend.png
///
///
/// The parameter is the name of friend chain (the name of a chain is always
/// the name of the tree from which it was created).
/// The original chain has access to all variable in its friends.
/// We can use the TChain::Draw method as if the values in the friends were
/// in the original chain.
/// To specify the chain to use in the Draw method, use the syntax:
/// ~~~ {.cpp}
///     <chainname>.<branchname>.<varname>
/// ~~~
/// If the variable name is enough to uniquely identify the variable, you can
/// leave out the chain and/or branch name.
/// For example, this generates a 3-d scatter plot of variable "var" in the
/// TChain ch versus variable v1 in TChain t1 versus variable v2 in TChain t2.
/// ~~~ {.cpp}
///     ch.Draw("var:t1.v1:t2.v2");
/// ~~~
/// When a TChain::Draw is executed, an automatic call to TTree::AddFriend
/// connects the trees in the chain. When a chain is deleted, its friend
/// elements are also deleted.
///
/// The number of entries in the friend must be equal or greater to the number
/// of entries of the original chain. If the friend has fewer entries a warning
/// is given and the resulting histogram will have missing entries.
/// For additional information see TTree::AddFriend.

TFriendElement* TChain::AddFriend(const char* chain, const char* dummy /* = "" */)
{
   if (!fFriends) {
      fFriends = new TList();
   }
   TFriendElement* fe = new TFriendElement(this, chain, dummy);

   R__ASSERT(fe); // There used to be a "if (fe)" test ... Keep this assert until we are sure that fe is never null

   fFriends->Add(fe);

   if (fProofChain)
      // This updates the proxy chain when we will really use PROOF
      ResetBit(kProofUptodate);

   // We need to invalidate the loading of the current tree because its list
   // of real friends is now obsolete.  It is repairable only from LoadTree.
   InvalidateCurrentTree();

   TTree* tree = fe->GetTree();
   if (!tree) {
      Warning("AddFriend", "Unknown TChain %s", chain);
   }
   return fe;
}

////////////////////////////////////////////////////////////////////////////////
/// Add the whole chain or tree as a friend of this chain.

TFriendElement* TChain::AddFriend(const char* chain, TFile* dummy)
{
   if (!fFriends) fFriends = new TList();
   TFriendElement *fe = new TFriendElement(this,chain,dummy);

   R__ASSERT(fe); // There used to be a "if (fe)" test ... Keep this assert until we are sure that fe is never null

   fFriends->Add(fe);

   if (fProofChain)
      // This updates the proxy chain when we will really use PROOF
      ResetBit(kProofUptodate);

   // We need to invalidate the loading of the current tree because its list
   // of real friend is now obsolete.  It is repairable only from LoadTree
   InvalidateCurrentTree();

   TTree *t = fe->GetTree();
   if (!t) {
      Warning("AddFriend","Unknown TChain %s",chain);
   }
   return fe;
}

////////////////////////////////////////////////////////////////////////////////
/// Add the whole chain or tree as a friend of this chain.

TFriendElement* TChain::AddFriend(TTree* chain, const char* alias, Bool_t /* warn = kFALSE */)
{
   if (!chain) return nullptr;
   if (!fFriends) fFriends = new TList();
   TFriendElement *fe = new TFriendElement(this,chain,alias);
   R__ASSERT(fe);

   fFriends->Add(fe);

   if (fProofChain)
      // This updates the proxy chain when we will really use PROOF
      ResetBit(kProofUptodate);

   // We need to invalidate the loading of the current tree because its list
   // of real friend is now obsolete.  It is repairable only from LoadTree
   InvalidateCurrentTree();

   TTree *t = fe->GetTree();
   if (!t) {
      Warning("AddFriend","Unknown TChain %s",chain->GetName());
   }
   chain->RegisterExternalFriend(fe);
   return fe;
}

////////////////////////////////////////////////////////////////////////////////
/// Browse the contents of the chain.

void TChain::Browse(TBrowser* b)
{
   TTree::Browse(b);
}

////////////////////////////////////////////////////////////////////////////////
/// When closing a file during the chain processing, the file
/// may be closed with option "R" if flag is set to kTRUE.
/// by default flag is kTRUE.
/// When closing a file with option "R", all TProcessIDs referenced by this
/// file are deleted.
/// Calling TFile::Close("R") might be necessary in case one reads a long list
/// of files having TRef, writing some of the referenced objects or TRef
/// to a new file. If the TRef or referenced objects of the file being closed
/// will not be referenced again, it is possible to minimize the size
/// of the TProcessID data structures in memory by forcing a delete of
/// the unused TProcessID.

void TChain::CanDeleteRefs(Bool_t flag /* = kTRUE */)
{
   fCanDeleteRefs = flag;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the packet descriptor string.

void TChain::CreatePackets()
{
   TIter next(fFiles);
   TChainElement* element = 0;
   while ((element = (TChainElement*) next())) {
      element->CreatePackets();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Override the TTree::DirectoryAutoAdd behavior:
/// we never auto add.

void TChain::DirectoryAutoAdd(TDirectory * /* dir */)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Draw expression varexp for selected entries.
/// Returns -1 in case of error or number of selected events in case of success.
///
/// This function accepts TCut objects as arguments.
/// Useful to use the string operator +, example:
/// ~~~{.cpp}
///    ntuple.Draw("x",cut1+cut2+cut3);
/// ~~~
///

Long64_t TChain::Draw(const char* varexp, const TCut& selection,
                      Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   if (fProofChain) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         SetProof(kTRUE, kTRUE);
      fProofChain->SetEventList(fEventList);
      fProofChain->SetEntryList(fEntryList);
      return fProofChain->Draw(varexp, selection, option, nentries, firstentry);
   }

   return TChain::Draw(varexp, selection.GetTitle(), option, nentries, firstentry);
}

////////////////////////////////////////////////////////////////////////////////
/// Process all entries in this chain and draw histogram corresponding to
/// expression varexp.
/// Returns -1 in case of error or number of selected events in case of success.

Long64_t TChain::Draw(const char* varexp, const char* selection,
                      Option_t* option,Long64_t nentries, Long64_t firstentry)
{
   if (fProofChain) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         SetProof(kTRUE, kTRUE);
      fProofChain->SetEventList(fEventList);
      fProofChain->SetEntryList(fEntryList);
      return fProofChain->Draw(varexp, selection, option, nentries, firstentry);
   }
   GetPlayer();
   if (LoadTree(firstentry) < 0) return 0;
   return TTree::Draw(varexp,selection,option,nentries,firstentry);
}

////////////////////////////////////////////////////////////////////////////////
/// See TTree::GetReadEntry().

TBranch* TChain::FindBranch(const char* branchname)
{
   if (fProofChain && !(fProofChain->TestBit(kProofLite))) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         SetProof(kTRUE, kTRUE);
      return fProofChain->FindBranch(branchname);
   }
   if (fTree) {
      return fTree->FindBranch(branchname);
   }
   LoadTree(0);
   if (fTree) {
      return fTree->FindBranch(branchname);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// See TTree::GetReadEntry().

TLeaf* TChain::FindLeaf(const char* searchname)
{
   if (fProofChain && !(fProofChain->TestBit(kProofLite))) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         SetProof(kTRUE, kTRUE);
      return fProofChain->FindLeaf(searchname);
   }
   if (fTree) {
      return fTree->FindLeaf(searchname);
   }
   LoadTree(0);
   if (fTree) {
      return fTree->FindLeaf(searchname);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the expanded value of the alias.  Search in the friends if any.

const char* TChain::GetAlias(const char* aliasName) const
{
   const char* alias = TTree::GetAlias(aliasName);
   if (alias) {
      return alias;
   }
   if (fTree) {
      return fTree->GetAlias(aliasName);
   }
   const_cast<TChain*>(this)->LoadTree(0);
   if (fTree) {
      return fTree->GetAlias(aliasName);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to the branch name in the current tree.

TBranch* TChain::GetBranch(const char* name)
{
   if (fProofChain && !(fProofChain->TestBit(kProofLite))) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         SetProof(kTRUE, kTRUE);
      return fProofChain->GetBranch(name);
   }
   if (fTree) {
      return fTree->GetBranch(name);
   }
   LoadTree(0);
   if (fTree) {
      return fTree->GetBranch(name);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// See TTree::GetReadEntry().

Bool_t TChain::GetBranchStatus(const char* branchname) const
{
   if (fProofChain && !(fProofChain->TestBit(kProofLite))) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         Warning("GetBranchStatus", "PROOF proxy not up-to-date:"
                                    " run TChain::SetProof(kTRUE, kTRUE) first");
      return fProofChain->GetBranchStatus(branchname);
   }
   return TTree::GetBranchStatus(branchname);
}

////////////////////////////////////////////////////////////////////////////////
/// Return an iterator over the cluster of baskets starting at firstentry.
///
/// This iterator is not yet supported for TChain object.

TTree::TClusterIterator TChain::GetClusterIterator(Long64_t /* firstentry */)
{
   Fatal("GetClusterIterator","Not support for TChain object");
   return TTree::GetClusterIterator(-1);
}

////////////////////////////////////////////////////////////////////////////////
/// Return absolute entry number in the chain.
/// The input parameter entry is the entry number in
/// the current tree of this chain.

Long64_t TChain::GetChainEntryNumber(Long64_t entry) const
{
   return entry + fTreeOffset[fTreeNumber];
}

////////////////////////////////////////////////////////////////////////////////
/// Return the total number of entries in the chain.
/// In case the number of entries in each tree is not yet known,
/// the offset table is computed.

Long64_t TChain::GetEntries() const
{
   if (fProofChain && !(fProofChain->TestBit(kProofLite))) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         Warning("GetEntries", "PROOF proxy not up-to-date:"
                               " run TChain::SetProof(kTRUE, kTRUE) first");
      return fProofChain->GetEntries();
   }
   if (fEntries == TTree::kMaxEntries) {
      const_cast<TChain*>(this)->LoadTree(TTree::kMaxEntries-1);
   }
   return fEntries;
}

////////////////////////////////////////////////////////////////////////////////
/// Get entry from the file to memory.
///
/// - getall = 0 : get only active branches
/// - getall = 1 : get all branches
///
/// Return the total number of bytes read,
/// 0 bytes read indicates a failure.

Int_t TChain::GetEntry(Long64_t entry, Int_t getall)
{
   Long64_t treeReadEntry = LoadTree(entry);
   if (treeReadEntry < 0) {
      return 0;
   }
   if (!fTree) {
      return 0;
   }
   return fTree->GetEntry(treeReadEntry, getall);
}

////////////////////////////////////////////////////////////////////////////////
/// Return entry number corresponding to entry.
///
/// if no TEntryList set returns entry
/// else returns entry \#entry from this entry list and
/// also computes the global entry number (loads all tree headers)

Long64_t TChain::GetEntryNumber(Long64_t entry) const
{

   if (fEntryList){
      Int_t treenum = 0;
      Long64_t localentry = fEntryList->GetEntryAndTree(entry, treenum);
      //find the global entry number
      //same const_cast as in the GetEntries() function
      if (localentry<0) return -1;
      if (treenum != fTreeNumber){
         if (fTreeOffset[treenum]==TTree::kMaxEntries){
            for (Int_t i=0; i<=treenum; i++){
               if (fTreeOffset[i]==TTree::kMaxEntries)
                  (const_cast<TChain*>(this))->LoadTree(fTreeOffset[i-1]);
            }
         }
         //(const_cast<TChain*>(this))->LoadTree(fTreeOffset[treenum]);
      }
      Long64_t globalentry = fTreeOffset[treenum] + localentry;
      return globalentry;
   }
   return entry;
}

////////////////////////////////////////////////////////////////////////////////
/// Return entry corresponding to major and minor number.
///
/// The function returns the total number of bytes read.
/// If the Tree has friend trees, the corresponding entry with
/// the index values (major,minor) is read. Note that the master Tree
/// and its friend may have different entry serial numbers corresponding
/// to (major,minor).

Int_t TChain::GetEntryWithIndex(Int_t major, Int_t minor)
{
   Long64_t serial = GetEntryNumberWithIndex(major, minor);
   if (serial < 0) return -1;
   return GetEntry(serial);
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the current file.
/// If no file is connected, the first file is automatically loaded.

TFile* TChain::GetFile() const
{
   if (fFile) {
      return fFile;
   }
   // Force opening the first file in the chain.
   const_cast<TChain*>(this)->LoadTree(0);
   return fFile;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the leaf name in the current tree.

TLeaf* TChain::GetLeaf(const char* branchname, const char *leafname)
{
   if (fProofChain && !(fProofChain->TestBit(kProofLite))) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         SetProof(kTRUE, kTRUE);
      return fProofChain->GetLeaf(branchname, leafname);
   }
   if (fTree) {
      return fTree->GetLeaf(branchname, leafname);
   }
   LoadTree(0);
   if (fTree) {
      return fTree->GetLeaf(branchname, leafname);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the leaf name in the current tree.

TLeaf* TChain::GetLeaf(const char* name)
{
   if (fProofChain && !(fProofChain->TestBit(kProofLite))) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         SetProof(kTRUE, kTRUE);
      return fProofChain->GetLeaf(name);
   }
   if (fTree) {
      return fTree->GetLeaf(name);
   }
   LoadTree(0);
   if (fTree) {
      return fTree->GetLeaf(name);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the list of branches of the current tree.
///
/// Warning: If there is no current TTree yet, this routine will open the
/// first in the chain.
///
/// Returns 0 on failure.

TObjArray* TChain::GetListOfBranches()
{
   if (fProofChain && !(fProofChain->TestBit(kProofLite))) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         SetProof(kTRUE, kTRUE);
      return fProofChain->GetListOfBranches();
   }
   if (fTree) {
      return fTree->GetListOfBranches();
   }
   LoadTree(0);
   if (fTree) {
      return fTree->GetListOfBranches();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the list of leaves of the current tree.
///
/// Warning: May set the current tree!

TObjArray* TChain::GetListOfLeaves()
{
   if (fProofChain && !(fProofChain->TestBit(kProofLite))) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         SetProof(kTRUE, kTRUE);
      return fProofChain->GetListOfLeaves();
   }
   if (fTree) {
      return fTree->GetListOfLeaves();
   }
   LoadTree(0);
   if (fTree) {
      return fTree->GetListOfLeaves();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return maximum of column with name columname.

Double_t TChain::GetMaximum(const char* columname)
{
   Double_t theMax = -DBL_MAX;
   for (Int_t file = 0; file < fNtrees; file++) {
      Long64_t first = fTreeOffset[file];
      LoadTree(first);
      Double_t curmax = fTree->GetMaximum(columname);
      if (curmax > theMax) {
         theMax = curmax;
      }
   }
   return theMax;
}

////////////////////////////////////////////////////////////////////////////////
/// Return minimum of column with name columname.

Double_t TChain::GetMinimum(const char* columname)
{
   Double_t theMin = DBL_MAX;
   for (Int_t file = 0; file < fNtrees; file++) {
      Long64_t first = fTreeOffset[file];
      LoadTree(first);
      Double_t curmin = fTree->GetMinimum(columname);
      if (curmin < theMin) {
         theMin = curmin;
      }
   }
   return theMin;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of branches of the current tree.
///
/// Warning: May set the current tree!

Int_t TChain::GetNbranches()
{
   if (fTree) {
      return fTree->GetNbranches();
   }
   LoadTree(0);
   if (fTree) {
      return fTree->GetNbranches();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// See TTree::GetReadEntry().

Long64_t TChain::GetReadEntry() const
{
   if (fProofChain && !(fProofChain->TestBit(kProofLite))) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         Warning("GetBranchStatus", "PROOF proxy not up-to-date:"
                                    " run TChain::SetProof(kTRUE, kTRUE) first");
      return fProofChain->GetReadEntry();
   }
   return TTree::GetReadEntry();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the chain weight.
///
/// By default the weight is the weight of the current tree.
/// However, if the weight has been set in TChain::SetWeight()
/// with the option "global", then that weight will be returned.
///
/// Warning: May set the current tree!

Double_t TChain::GetWeight() const
{
   if (TestBit(kGlobalWeight)) {
      return fWeight;
   } else {
      if (fTree) {
         return fTree->GetWeight();
      }
      const_cast<TChain*>(this)->LoadTree(0);
      if (fTree) {
         return fTree->GetWeight();
      }
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the TTree to be reloaded as soon as possible.  In particular this
/// is needed when adding a Friend.
///
/// If the tree has clones, copy them into the chain
/// clone list so we can change their branch addresses
/// when necessary.
///
/// This is to support the syntax:
/// ~~~ {.cpp}
///     TTree* clone = chain->GetTree()->CloneTree(0);
/// ~~~

void TChain::InvalidateCurrentTree()
{
   if (fTree && fTree->GetListOfClones()) {
      for (TObjLink* lnk = fTree->GetListOfClones()->FirstLink(); lnk; lnk = lnk->Next()) {
         TTree* clone = (TTree*) lnk->GetObject();
         AddClone(clone);
      }
   }
   fTreeNumber = -1;
   fTree = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Dummy function.
/// It could be implemented and load all baskets of all trees in the chain.
/// For the time being use TChain::Merge and TTree::LoadBasket
/// on the resulting tree.

Int_t TChain::LoadBaskets(Long64_t /*maxmemory*/)
{
   Error("LoadBaskets", "Function not yet implemented for TChain.");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Find the tree which contains entry, and set it as the current tree.
///
/// Returns the entry number in that tree.
///
/// The input argument entry is the entry serial number in the whole chain.
///
/// In case of error, LoadTree returns a negative number:
///   * -1: The chain is empty.
///   * -2: The requested entry number is less than zero or too large for the chain.
///       or too large for the large TTree.
///   * -3: The file corresponding to the entry could not be correctly open
///   * -4: The TChainElement corresponding to the entry is missing or
///       the TTree is missing from the file.
///   * -5: Internal error, please report the circumstance when this happen
///       as a ROOT issue.
///   * -6: An error occurred within the notify callback.
///
/// Note: This is the only routine which sets the value of fTree to
///       a non-zero pointer.

Long64_t TChain::LoadTree(Long64_t entry)
{
   // We already have been visited while recursively looking
   // through the friends tree, let's return.
   if (kLoadTree & fFriendLockStatus) {
      return 0;
   }

   if (!fNtrees) {
      // -- The chain is empty.
      return -1;
   }

   if ((entry < 0) || ((entry > 0) && (entry >= fEntries && entry!=(TTree::kMaxEntries-1) ))) {
      // -- Invalid entry number.
      if (fTree) fTree->LoadTree(-1);
      fReadEntry = -1;
      return -2;
   }

   // Find out which tree in the chain contains the passed entry.
   Int_t treenum = fTreeNumber;
   if ((fTreeNumber == -1) || (entry < fTreeOffset[fTreeNumber]) || (entry >= fTreeOffset[fTreeNumber+1]) || (entry==TTree::kMaxEntries-1)) {
      // -- Entry is *not* in the chain's current tree.
      // Do a linear search of the tree offset array.
      // FIXME: We could be smarter by starting at the
      //        current tree number and going forwards,
      //        then wrapping around at the end.
      for (treenum = 0; treenum < fNtrees; treenum++) {
         if (entry < fTreeOffset[treenum+1]) {
            break;
         }
      }
   }

   // Calculate the entry number relative to the found tree.
   Long64_t treeReadEntry = entry - fTreeOffset[treenum];
   fReadEntry = entry;

   // If entry belongs to the current tree return entry.
   if (fTree && treenum == fTreeNumber) {
      // First set the entry the tree on its owns friends
      // (the friends of the chain will be updated in the
      // next loop).
      fTree->LoadTree(treeReadEntry);
      if (fFriends) {
         // The current tree has not changed but some of its friends might.
         //
         TIter next(fFriends);
         TFriendLock lock(this, kLoadTree);
         TFriendElement* fe = 0;
         while ((fe = (TFriendElement*) next())) {
            TTree* at = fe->GetTree();
            // If the tree is a
            // direct friend of the chain, it should be scanned
            // used the chain entry number and NOT the tree entry
            // number (treeReadEntry) hence we do:
            at->LoadTreeFriend(entry, this);
         }
         Bool_t needUpdate = kFALSE;
         if (fTree->GetListOfFriends()) {
            for(auto fetree : ROOT::Detail::TRangeStaticCast<TFriendElement>(*fTree->GetListOfFriends())) {
               if (fetree->IsUpdated()) {
                  needUpdate = kTRUE;
                  fetree->ResetUpdated();
               }
            }
         }
         if (needUpdate) {
            // Update the branch/leaf addresses and
            // thelist of leaves in all TTreeFormula of the TTreePlayer (if any).

            // Set the branch statuses for the newly opened file.
            TChainElement *frelement;
            TIter fnext(fStatus);
            while ((frelement = (TChainElement*) fnext())) {
               Int_t status = frelement->GetStatus();
               fTree->SetBranchStatus(frelement->GetName(), status);
            }

            // Set the branch addresses for the newly opened file.
            fnext.Reset();
            while ((frelement = (TChainElement*) fnext())) {
               void* addr = frelement->GetBaddress();
               if (addr) {
                  TBranch* br = fTree->GetBranch(frelement->GetName());
                  TBranch** pp = frelement->GetBranchPtr();
                  if (pp) {
                     // FIXME: What if br is zero here?
                     *pp = br;
                  }
                  if (br) {
                     if (!frelement->GetCheckedType()) {
                        Int_t res = CheckBranchAddressType(br, TClass::GetClass(frelement->GetBaddressClassName()),
                                                         (EDataType) frelement->GetBaddressType(), frelement->GetBaddressIsPtr());
                        if ((res & kNeedEnableDecomposedObj) && !br->GetMakeClass()) {
                           br->SetMakeClass(kTRUE);
                        }
                        frelement->SetDecomposedObj(br->GetMakeClass());
                        frelement->SetCheckedType(kTRUE);
                     }
                     // FIXME: We may have to tell the branch it should
                     //        not be an owner of the object pointed at.
                     br->SetAddress(addr);
                     if (TestBit(kAutoDelete)) {
                        br->SetAutoDelete(kTRUE);
                     }
                  }
               }
            }
            if (fPlayer) {
               fPlayer->UpdateFormulaLeaves();
            }
            // Notify user if requested.
            if (fNotify) {
               if(!fNotify->Notify()) return -6;
            }
         }
      }
      return treeReadEntry;
   }

   if (fExternalFriends) {
      for(auto external_fe : ROOT::Detail::TRangeStaticCast<TFriendElement>(*fExternalFriends)) {
         external_fe->MarkUpdated();
      }
   }

   // Delete the current tree and open the new tree.
   TTreeCache* tpf = 0;
   // Delete file unless the file owns this chain!
   // FIXME: The "unless" case here causes us to leak memory.
   if (fFile) {
      if (!fDirectory->GetList()->FindObject(this)) {
         if (fTree) {
            // (fFile != 0 && fTree == 0) can happen when
            // InvalidateCurrentTree is called (for example from
            // AddFriend).  Having fTree === 0 is necessary in that
            // case because in some cases GetTree is used as a check
            // to see if a TTree is already loaded.
            // However, this prevent using the following to reuse
            // the TTreeCache object.
            tpf = fTree->GetReadCache(fFile);
            if (tpf) {
               tpf->ResetCache();
            }

            fFile->SetCacheRead(0, fTree);
            // If the tree has clones, copy them into the chain
            // clone list so we can change their branch addresses
            // when necessary.
            //
            // This is to support the syntax:
            //
            //      TTree* clone = chain->GetTree()->CloneTree(0);
            //
            // We need to call the invalidate exactly here, since
            // we no longer need the value of fTree and it is
            // about to be deleted.
            InvalidateCurrentTree();
         }

         if (fCanDeleteRefs) {
            fFile->Close("R");
         }
         delete fFile;
         fFile = 0;
      } else {
         // If the tree has clones, copy them into the chain
         // clone list so we can change their branch addresses
         // when necessary.
         //
         // This is to support the syntax:
         //
         //      TTree* clone = chain->GetTree()->CloneTree(0);
         //
         if (fTree) InvalidateCurrentTree();
      }
   }

   TChainElement* element = (TChainElement*) fFiles->At(treenum);
   if (!element) {
      if (treeReadEntry) {
         return -4;
      }
      // Last attempt, just in case all trees in the chain have 0 entries.
      element = (TChainElement*) fFiles->At(0);
      if (!element) {
         return -4;
      }
   }

   // FIXME: We leak memory here, we've just lost the open file
   //        if we did not delete it above.
   {
      TDirectory::TContext ctxt;
      fFile = TFile::Open(element->GetTitle());
      if (fFile) fFile->SetBit(kMustCleanup);
   }

   // ----- Begin of modifications by MvL
   Int_t returnCode = 0;
   if (!fFile || fFile->IsZombie()) {
      if (fFile) {
         delete fFile;
         fFile = 0;
      }
      // Note: We do *not* own fTree.
      fTree = 0;
      returnCode = -3;
   } else {
      if (fPerfStats)
         fPerfStats->SetFile(fFile);

      // Note: We do *not* own fTree after this, the file does!
      fTree = dynamic_cast<TTree*>(fFile->Get(element->GetName()));
      if (!fTree) {
         // Now that we do not check during the addition, we need to check here!
         Error("LoadTree", "Cannot find tree with name %s in file %s", element->GetName(), element->GetTitle());
         delete fFile;
         fFile = 0;
         // We do not return yet so that 'fEntries' can be updated with the
         // sum of the entries of all the other trees.
         returnCode = -4;
      }
   }

   fTreeNumber = treenum;
   // FIXME: We own fFile, we must be careful giving away a pointer to it!
   // FIXME: We may set fDirectory to zero here!
   fDirectory = fFile;

   // Reuse cache from previous file (if any).
   if (tpf) {
      if (fFile) {
         // FIXME: fTree may be zero here.
         tpf->UpdateBranches(fTree);
         tpf->ResetCache();
         fFile->SetCacheRead(tpf, fTree);
      } else {
         // FIXME: One of the file in the chain is missing
         // we have no place to hold the pointer to the
         // TTreeCache.
         delete tpf;
         tpf = 0;
      }
   } else {
      if (fCacheUserSet) {
         this->SetCacheSize(fCacheSize);
      }
   }

   // Check if fTreeOffset has really been set.
   Long64_t nentries = 0;
   if (fTree) {
      nentries = fTree->GetEntries();
   }

   if (fTreeOffset[fTreeNumber+1] != (fTreeOffset[fTreeNumber] + nentries)) {
      fTreeOffset[fTreeNumber+1] = fTreeOffset[fTreeNumber] + nentries;
      fEntries = fTreeOffset[fNtrees];
      element->SetNumberEntries(nentries);
      // Below we must test >= in case the tree has no entries.
      if (entry >= fTreeOffset[fTreeNumber+1]) {
         if ((fTreeNumber < (fNtrees - 1)) && (entry < fTreeOffset[fTreeNumber+2])) {
            // The request entry is not in the tree 'fTreeNumber' we will need
            // to look further.

            // Before moving on, let's record the result.
            element->SetLoadResult(returnCode);

            // Before trying to read the file file/tree, notify the user
            // that we have switched trees if requested; the user might need
            // to properly account for the number of files/trees even if they
            // have no entries.
            if (fNotify) {
               if(!fNotify->Notify()) return -6;
            }

            // Load the next TTree.
            return LoadTree(entry);
         } else {
            treeReadEntry = fReadEntry = -2;
         }
      }
   }


   if (!fTree) {
      // The Error message already issued.  However if we reach here
      // we need to make sure that we do not use fTree.
      //
      // Force a reload of the tree next time.
      fTreeNumber = -1;

      element->SetLoadResult(returnCode);
      return returnCode;
   }
   // ----- End of modifications by MvL

   // Copy the chain's clone list into the new tree's
   // clone list so that branch addresses stay synchronized.
   if (fClones) {
      for (TObjLink* lnk = fClones->FirstLink(); lnk; lnk = lnk->Next()) {
         TTree* clone = (TTree*) lnk->GetObject();
         ((TChain*) fTree)->TTree::AddClone(clone);
      }
   }

   // Since some of the friends of this chain might simple trees
   // (i.e., not really chains at all), we need to execute this
   // before calling LoadTree(entry) on the friends (so that
   // they use the correct read entry number).

   // Change the new current tree to the new entry.
   Long64_t loadResult = fTree->LoadTree(treeReadEntry);
   if (loadResult == treeReadEntry) {
      element->SetLoadResult(0);
   } else {
      // This is likely to be an internal error, if treeReadEntry was not in range
      // (or intentionally -2 for TChain::GetEntries) then something happened
      // that is very odd/surprising.
      element->SetLoadResult(-5);
   }


   // Change the chain friends to the new entry.
   if (fFriends) {
      // An alternative would move this code to each of the function
      // calling LoadTree (and to overload a few more).
      TIter next(fFriends);
      TFriendLock lock(this, kLoadTree);
      TFriendElement* fe = 0;
      while ((fe = (TFriendElement*) next())) {
         TTree* t = fe->GetTree();
         if (!t) continue;
         if (t->GetTreeIndex()) {
            t->GetTreeIndex()->UpdateFormulaLeaves(0);
         }
         if (t->GetTree() && t->GetTree()->GetTreeIndex()) {
            t->GetTree()->GetTreeIndex()->UpdateFormulaLeaves(GetTree());
         }
         t->LoadTreeFriend(entry, this);
         TTree* friend_t = t->GetTree();
         if (friend_t) {
            auto localfe = fTree->AddFriend(t, fe->GetName());
            localfe->SetBit(TFriendElement::kFromChain);
         }
      }
   }

   fTree->SetMakeClass(fMakeClass);
   fTree->SetMaxVirtualSize(fMaxVirtualSize);

   SetChainOffset(fTreeOffset[fTreeNumber]);

   // Set the branch statuses for the newly opened file.
   TIter next(fStatus);
   while ((element = (TChainElement*) next())) {
      Int_t status = element->GetStatus();
      fTree->SetBranchStatus(element->GetName(), status);
   }

   // Set the branch addresses for the newly opened file.
   next.Reset();
   while ((element = (TChainElement*) next())) {
      void* addr = element->GetBaddress();
      if (addr) {
         TBranch* br = fTree->GetBranch(element->GetName());
         TBranch** pp = element->GetBranchPtr();
         if (pp) {
            // FIXME: What if br is zero here?
            *pp = br;
         }
         if (br) {
            if (!element->GetCheckedType()) {
               Int_t res = CheckBranchAddressType(br, TClass::GetClass(element->GetBaddressClassName()),
                                                  (EDataType) element->GetBaddressType(), element->GetBaddressIsPtr());
               if ((res & kNeedEnableDecomposedObj) && !br->GetMakeClass()) {
                  br->SetMakeClass(kTRUE);
               }
               element->SetDecomposedObj(br->GetMakeClass());
               element->SetCheckedType(kTRUE);
            }
            // FIXME: We may have to tell the branch it should
            //        not be an owner of the object pointed at.
            br->SetAddress(addr);
            if (TestBit(kAutoDelete)) {
               br->SetAutoDelete(kTRUE);
            }
         }
      }
   }

   // Update the addresses of the chain's cloned trees, if any.
   if (fClones) {
      for (TObjLink* lnk = fClones->FirstLink(); lnk; lnk = lnk->Next()) {
         TTree* clone = (TTree*) lnk->GetObject();
         CopyAddresses(clone);
      }
   }

   // Update list of leaves in all TTreeFormula's of the TTreePlayer (if any).
   if (fPlayer) {
      fPlayer->UpdateFormulaLeaves();
   }

   // Notify user we have switched trees if requested.
   if (fNotify) {
      if(!fNotify->Notify()) return -6;
   }

   // Return the new local entry number.
   return treeReadEntry;
}

////////////////////////////////////////////////////////////////////////////////
/// Check / locate the files in the chain.
/// By default only the files not yet looked up are checked.
/// Use force = kTRUE to check / re-check every file.

void TChain::Lookup(Bool_t force)
{
   TIter next(fFiles);
   TChainElement* element = 0;
   Int_t nelements = fFiles->GetEntries();
   printf("\n");
   printf("TChain::Lookup - Looking up %d files .... \n", nelements);
   Int_t nlook = 0;
   TFileStager *stg = 0;
   while ((element = (TChainElement*) next())) {
      // Do not do it more than needed
      if (element->HasBeenLookedUp() && !force) continue;
      // Count
      nlook++;
      // Get the Url
      TUrl elemurl(element->GetTitle(), kTRUE);
      // Save current options and anchor
      TString anchor = elemurl.GetAnchor();
      TString options = elemurl.GetOptions();
      // Reset options and anchor
      elemurl.SetOptions("");
      elemurl.SetAnchor("");
      // Locate the file
      TString eurl(elemurl.GetUrl());
      if (!stg || !stg->Matches(eurl)) {
         SafeDelete(stg);
         {
            TDirectory::TContext ctxt;
            stg = TFileStager::Open(eurl);
         }
         if (!stg) {
            Error("Lookup", "TFileStager instance cannot be instantiated");
            break;
         }
      }
      Int_t n1 = (nelements > 100) ? (Int_t) nelements / 100 : 1;
      if (stg->Locate(eurl.Data(), eurl) == 0) {
         if (nlook > 0 && !(nlook % n1)) {
            printf("Lookup | %3d %% finished\r", 100 * nlook / nelements);
            fflush(stdout);
         }
         // Get the effective end-point Url
         elemurl.SetUrl(eurl);
         // Restore original options and anchor, if any
         elemurl.SetOptions(options);
         elemurl.SetAnchor(anchor);
         // Save it into the element
         element->SetTitle(elemurl.GetUrl());
         // Remember
         element->SetLookedUp();
      } else {
         // Failure: remove
         fFiles->Remove(element);
         if (gSystem->AccessPathName(eurl))
            Error("Lookup", "file %s does not exist\n", eurl.Data());
         else
            Error("Lookup", "file %s cannot be read\n", eurl.Data());
      }
   }
   if (nelements > 0)
      printf("Lookup | %3d %% finished\n", 100 * nlook / nelements);
   else
      printf("\n");
   fflush(stdout);
   SafeDelete(stg);
}

////////////////////////////////////////////////////////////////////////////////
/// Loop on nentries of this chain starting at firstentry.  (NOT IMPLEMENTED)

void TChain::Loop(Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   Error("Loop", "Function not yet implemented");

   if (option || nentries || firstentry) { }  // keep warnings away

#if 0
   if (LoadTree(firstentry) < 0) return;

   if (firstentry < 0) firstentry = 0;
   Long64_t lastentry = firstentry + nentries -1;
   if (lastentry > fEntries-1) {
      lastentry = fEntries -1;
   }

   GetPlayer();
   GetSelector();
   fSelector->Start(option);

   Long64_t entry = firstentry;
   Int_t tree,e0,en;
   for (tree=0;tree<fNtrees;tree++) {
      e0 = fTreeOffset[tree];
      en = fTreeOffset[tree+1] - 1;
      if (en > lastentry) en = lastentry;
      if (entry > en) continue;

      LoadTree(entry);
      fSelector->BeginFile();

      while (entry <= en) {
         fSelector->Execute(fTree, entry - e0);
         entry++;
      }
      fSelector->EndFile();
   }

   fSelector->Finish(option);
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// List the chain.

void TChain::ls(Option_t* option) const
{
   TObject::ls(option);
   TIter next(fFiles);
   TChainElement* file = 0;
   TROOT::IncreaseDirLevel();
   while ((file = (TChainElement*)next())) {
      file->ls(option);
   }
   TROOT::DecreaseDirLevel();
}

////////////////////////////////////////////////////////////////////////////////
/// Merge all the entries in the chain into a new tree in a new file.
///
/// See important note in the following function Merge().
///
/// If the chain is expecting the input tree inside a directory,
/// this directory is NOT created by this routine.
///
/// So in a case where we have:
/// ~~~ {.cpp}
///     TChain ch("mydir/mytree");
///     ch.Merge("newfile.root");
/// ~~~
/// The resulting file will have not subdirectory. To recreate
/// the directory structure do:
/// ~~~ {.cpp}
///     TFile* file = TFile::Open("newfile.root", "RECREATE");
///     file->mkdir("mydir")->cd();
///     ch.Merge(file);
/// ~~~

Long64_t TChain::Merge(const char* name, Option_t* option)
{
   TFile *file = TFile::Open(name, "recreate", "chain files", 1);
   return Merge(file, 0, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Merge all chains in the collection.  (NOT IMPLEMENTED)

Long64_t TChain::Merge(TCollection* /* list */, Option_t* /* option */ )
{
   Error("Merge", "not implemented");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Merge all chains in the collection.  (NOT IMPLEMENTED)

Long64_t TChain::Merge(TCollection* /* list */, TFileMergeInfo *)
{
   Error("Merge", "not implemented");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Merge all the entries in the chain into a new tree in the current file.
///
/// Note: The "file" parameter is *not* the file where the new
///       tree will be inserted.  The new tree is inserted into
///       gDirectory, which is usually the most recently opened
///       file, or the directory most recently cd()'d to.
///
/// If option = "C" is given, the compression level for all branches
/// in the new Tree is set to the file compression level.  By default,
/// the compression level of all branches is the original compression
/// level in the old trees.
///
/// If basketsize > 1000, the basket size for all branches of the
/// new tree will be set to basketsize.
///
/// Example using the file generated in $ROOTSYS/test/Event
/// merge two copies of Event.root
/// ~~~ {.cpp}
///     gSystem.Load("libEvent");
///     TChain ch("T");
///     ch.Add("Event1.root");
///     ch.Add("Event2.root");
///     ch.Merge("all.root");
/// ~~~
/// If the chain is expecting the input tree inside a directory,
/// this directory is NOT created by this routine.
///
/// So if you do:
/// ~~~ {.cpp}
///     TChain ch("mydir/mytree");
///     ch.Merge("newfile.root");
/// ~~~
/// The resulting file will not have subdirectories.  In order to
/// preserve the directory structure do the following instead:
/// ~~~ {.cpp}
///     TFile* file = TFile::Open("newfile.root", "RECREATE");
///     file->mkdir("mydir")->cd();
///     ch.Merge(file);
/// ~~~
/// If 'option' contains the word 'fast' the merge will be done without
/// unzipping or unstreaming the baskets (i.e., a direct copy of the raw
/// bytes on disk).
///
/// When 'fast' is specified, 'option' can also contains a
/// sorting order for the baskets in the output file.
///
/// There is currently 3 supported sorting order:
/// ~~~ {.cpp}
///     SortBasketsByOffset (the default)
///     SortBasketsByBranch
///     SortBasketsByEntry
/// ~~~
/// When using SortBasketsByOffset the baskets are written in
/// the output file in the same order as in the original file
/// (i.e. the basket are sorted on their offset in the original
/// file; Usually this also means that the baskets are sorted
/// on the index/number of the _last_ entry they contain)
///
/// When using SortBasketsByBranch all the baskets of each
/// individual branches are stored contiguously.  This tends to
/// optimize reading speed when reading a small number (1->5) of
/// branches, since all their baskets will be clustered together
/// instead of being spread across the file.  However it might
/// decrease the performance when reading more branches (or the full
/// entry).
///
/// When using SortBasketsByEntry the baskets with the lowest
/// starting entry are written first.  (i.e. the baskets are
/// sorted on the index/number of the first entry they contain).
/// This means that on the file the baskets will be in the order
/// in which they will be needed when reading the whole tree
/// sequentially.
///
/// ## IMPORTANT Note 1: AUTOMATIC FILE OVERFLOW
///
/// When merging many files, it may happen that the resulting file
/// reaches a size > TTree::fgMaxTreeSize (default = 100 GBytes).
/// In this case the current file is automatically closed and a new
/// file started.  If the name of the merged file was "merged.root",
/// the subsequent files will be named "merged_1.root", "merged_2.root",
/// etc.  fgMaxTreeSize may be modified via the static function
/// TTree::SetMaxTreeSize.
/// When in fast mode, the check and switch is only done in between each
/// input file.
///
/// ## IMPORTANT Note 2: The output file is automatically closed and deleted.
///
/// This is required because in general the automatic file overflow described
/// above may happen during the merge.
/// If only the current file is produced (the file passed as first argument),
/// one can instruct Merge to not close and delete the file by specifying
/// the option "keep".
///
/// The function returns the total number of files produced.
/// To check that all files have been merged use something like:
/// ~~~ {.cpp}
///     if (newchain->GetEntries()!=oldchain->GetEntries()) {
///        ... not all the file have been copied ...
///     }
/// ~~~

Long64_t TChain::Merge(TFile* file, Int_t basketsize, Option_t* option)
{
   // We must have been passed a file, we will use it
   // later to reset the compression level of the branches.
   if (!file) {
      // FIXME: We need an error message here.
      return 0;
   }

   // Options
   Bool_t fastClone = kFALSE;
   TString opt = option;
   opt.ToLower();
   if (opt.Contains("fast")) {
      fastClone = kTRUE;
   }

   // The chain tree must have a list of branches
   // because we may try to change their basket
   // size later.
   TObjArray* lbranches = GetListOfBranches();
   if (!lbranches) {
      // FIXME: We need an error message here.
      return 0;
   }

   // The chain must have a current tree because
   // that is the one we will clone.
   if (!fTree) {
      // -- LoadTree() has not yet been called, no current tree.
      // FIXME: We need an error message here.
      return 0;
   }

   // Copy the chain's current tree without
   // copying any entries, we will do that later.
   TTree* newTree = CloneTree(0);
   if (!newTree) {
      // FIXME: We need an error message here.
      return 0;
   }

   // Strip out the (potential) directory name.
   // FIXME: The merged chain may or may not have the
   //        same name as the original chain.  This is
   //        bad because the chain name determines the
   //        names of the trees in the chain by default.
   newTree->SetName(gSystem->BaseName(GetName()));

   // FIXME: Why do we do this?
   newTree->SetAutoSave(2000000000);

   // Circularity is incompatible with merging, it may
   // force us to throw away entries, which is not what
   // we are supposed to do.
   newTree->SetCircular(0);

   // Reset the compression level of the branches.
   if (opt.Contains("c")) {
      TBranch* branch = 0;
      TIter nextb(newTree->GetListOfBranches());
      while ((branch = (TBranch*) nextb())) {
         branch->SetCompressionSettings(file->GetCompressionSettings());
      }
   }

   // Reset the basket size of the branches.
   if (basketsize > 1000) {
      TBranch* branch = 0;
      TIter nextb(newTree->GetListOfBranches());
      while ((branch = (TBranch*) nextb())) {
         branch->SetBasketSize(basketsize);
      }
   }

   // Copy the entries.
   if (fastClone) {
      if ( newTree->CopyEntries( this, -1, option ) < 0 ) {
         // There was a problem!
         Error("Merge", "TTree has not been cloned\n");
      }
   } else {
      newTree->CopyEntries( this, -1, option );
   }

   // Write the new tree header.
   newTree->Write();

   // Get our return value.
   Int_t nfiles = newTree->GetFileNumber() + 1;

   // Close and delete the current file of the new tree.
   if (!opt.Contains("keep")) {
      // Delete the currentFile and the TTree object.
      delete newTree->GetCurrentFile();
   }
   return nfiles;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the tree url or filename and other information from the name
///
/// A treename and a url's query section is split off from name. The
/// splitting depends on whether the resulting filename is to be
/// subsequently treated for wildcards or not, since the question mark is
/// both the url query identifier and a wildcard. Wildcard matching is not
/// done in this method itself.
/// ~~~ {.cpp}
///     [xxx://host]/a/path/file_name[?query[#treename]]
/// ~~~
///
/// The following way to specify the treename is still supported with the
/// constrain that the file name contains the sub-string '.root'.
/// This is now deprecated and will be removed in future versions.
/// ~~~ {.cpp}
///     [xxx://host]/a/path/file.root[.oext][/treename]
///     [xxx://host]/a/path/file.root[.oext][/treename][?query]
/// ~~~
///
/// Note that in a case like this
/// ~~~ {.cpp}
///     [xxx://host]/a/path/file#treename
/// ~~~
/// i.e. anchor but no options (query), the filename will be the full path, as
/// the anchor may be the internal file name of an archive. Use '?#treename' to
/// pass the treename if the query field is empty.
///
/// \param[in] name        is the original name
/// \param[in] wildcards   indicates if the resulting filename will be treated for
///                        wildcards. For backwards compatibility, with most protocols
///                        this flag suppresses the search for the url fragment
///                        identifier and limits the query identifier search to cases
///                        where the tree name is given as a trailing slash-separated
///                        string at the end of the file name.
/// \param[out] filename   the url or filename to be opened or matched
/// \param[out] treename   the treename, which may be found in a url fragment section
///                        as a trailing part of the name (deprecated).
///                        If not found this will be empty.
/// \param[out] query      is the url query section, including the leading question
///                        mark. If not found or the query section is only followed by
///                        a fragment this will be empty.
/// \param[out] suffix     the portion of name which was removed to from filename.

void TChain::ParseTreeFilename(const char *name, TString &filename, TString &treename, TString &query, TString &suffix,
                               Bool_t) const
{
   Ssiz_t pIdx = kNPOS;
   filename.Clear();
   treename.Clear();
   query.Clear();
   suffix.Clear();

   // General case
   TUrl url(name, kTRUE);
   filename = (strcmp(url.GetProtocol(), "file")) ? url.GetUrl() : url.GetFileAndOptions();

   TString fn = url.GetFile();
   // Extract query, if any
   if (url.GetOptions() && (strlen(url.GetOptions()) > 0))
      query.Form("?%s", url.GetOptions());
   // The treename can be passed as anchor
   if (url.GetAnchor() && (strlen(url.GetAnchor()) > 0)) {
      // Support "?#tree_name" and "?query#tree_name"
      // "#tree_name" (no '?' is for tar archives)
      if (!query.IsNull() || strstr(name, "?#")) {
         treename = url.GetAnchor();
      } else {
         // The anchor is part of the file name
         fn = url.GetFileAndOptions();
      }
   }
   // Suffix
   suffix = url.GetFileAndOptions();
   // Get options from suffix by removing the file name
   suffix.Replace(suffix.Index(fn), fn.Length(), "");
   // Remove the options suffix from the original file name
   filename.Replace(filename.Index(suffix), suffix.Length(), "");

   // Special case: [...]file.root/treename
   static const char *dotr = ".root";
   static Ssiz_t dotrl = strlen(dotr);
   // Find the last one
   Ssiz_t js = filename.Index(dotr);
   while (js != kNPOS) {
      pIdx = js;
      js = filename.Index(dotr, js + 1);
   }
   if (pIdx != kNPOS) {
      static const char *slash = "/";
      static Ssiz_t slashl = strlen(slash);
      // Find the last one
      Ssiz_t ppIdx = filename.Index(slash, pIdx + dotrl);
      if (ppIdx != kNPOS) {
         // Good treename with the old receipe
         treename = filename(ppIdx + slashl, filename.Length());
         filename.Remove(ppIdx + slashl - 1);
         suffix.Insert(0, TString::Format("/%s", treename.Data()));
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print the header information of each tree in the chain.
/// See TTree::Print for a list of options.

void TChain::Print(Option_t *option) const
{
   TIter next(fFiles);
   TChainElement *element;
   while ((element = (TChainElement*)next())) {
      Printf("******************************************************************************");
      Printf("*Chain   :%-10s: %-54s *", GetName(), element->GetTitle());
      Printf("******************************************************************************");
      TFile *file = TFile::Open(element->GetTitle());
      if (file && !file->IsZombie()) {
         TTree *tree = (TTree*)file->Get(element->GetName());
         if (tree) tree->Print(option);
      }
      delete file;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Process all entries in this chain, calling functions in filename.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.
/// See TTree::Process.

Long64_t TChain::Process(const char *filename, Option_t *option, Long64_t nentries, Long64_t firstentry)
{
   if (fProofChain) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         SetProof(kTRUE, kTRUE);
      fProofChain->SetEventList(fEventList);
      fProofChain->SetEntryList(fEntryList);
      return fProofChain->Process(filename, option, nentries, firstentry);
   }

   if (LoadTree(firstentry) < 0) {
      return 0;
   }
   return TTree::Process(filename, option, nentries, firstentry);
}

////////////////////////////////////////////////////////////////////////////////
/// Process this chain executing the code in selector.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.

Long64_t TChain::Process(TSelector* selector, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   if (fProofChain) {
      // Make sure the element list is uptodate
      if (!TestBit(kProofUptodate))
         SetProof(kTRUE, kTRUE);
      fProofChain->SetEventList(fEventList);
      fProofChain->SetEntryList(fEntryList);
      return fProofChain->Process(selector, option, nentries, firstentry);
   }

   return TTree::Process(selector, option, nentries, firstentry);
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that obj (which is being deleted or will soon be) is no
/// longer referenced by this TTree.

void TChain::RecursiveRemove(TObject *obj)
{
   if (fFile == obj) {
      fFile = 0;
      fDirectory = 0;
      fTree = 0;
   }
   if (fDirectory == obj) {
      fDirectory = 0;
      fTree = 0;
   }
   if (fTree == obj) {
      fTree = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a friend from the list of friends.

void TChain::RemoveFriend(TTree* oldFriend)
{
   // We already have been visited while recursively looking
   // through the friends tree, let return

   if (!fFriends) {
      return;
   }

   TTree::RemoveFriend(oldFriend);

   if (fProofChain)
      // This updates the proxy chain when we will really use PROOF
      ResetBit(kProofUptodate);

   // We need to invalidate the loading of the current tree because its list
   // of real friends is now obsolete.  It is repairable only from LoadTree.
   InvalidateCurrentTree();
}

////////////////////////////////////////////////////////////////////////////////
/// Resets the state of this chain.

void TChain::Reset(Option_t*)
{
   delete fFile;
   fFile = 0;
   fNtrees         = 0;
   fTreeNumber     = -1;
   fTree           = 0;
   fFile           = 0;
   fFiles->Delete();
   fStatus->Delete();
   fTreeOffset[0]  = 0;
   TChainElement* element = new TChainElement("*", "");
   fStatus->Add(element);
   fDirectory = 0;

   TTree::Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Resets the state of this chain after a merge (keep the customization but
/// forget the data).

void TChain::ResetAfterMerge(TFileMergeInfo *info)
{
   fNtrees         = 0;
   fTreeNumber     = -1;
   fTree           = 0;
   fFile           = 0;
   fFiles->Delete();
   fTreeOffset[0]  = 0;

   TTree::ResetAfterMerge(info);
}

////////////////////////////////////////////////////////////////////////////////
/// Save TChain as a C++ statements on output stream out.
/// With the option "friend" save the description of all the
/// TChain's friend trees or chains as well.

void TChain::SavePrimitive(std::ostream &out, Option_t *option)
{
   static Int_t chCounter = 0;

   TString chName = gInterpreter->MapCppName(GetName());
   if (chName.IsNull())
      chName = "_chain";
   ++chCounter;
   chName += chCounter;

   TString opt = option;
   opt.ToLower();

   out << "   TChain *" << chName.Data() << " = new TChain(\"" << GetName() << "\");" << std::endl;

   if (opt.Contains("friend")) {
      opt.ReplaceAll("friend", "");
      for (TObject *frel : *fFriends) {
         TTree *frtree = ((TFriendElement *)frel)->GetTree();
         if (dynamic_cast<TChain *>(frtree)) {
            if (strcmp(frtree->GetName(), GetName()) != 0)
               --chCounter; // make friends get the same chain counter
            frtree->SavePrimitive(out, opt.Data());
            out << "   " << chName.Data() << "->AddFriend(\"" << frtree->GetName() << "\");" << std::endl;
         } else { // ordinary friend TTree
            TDirectory *file = frtree->GetDirectory();
            if (file && dynamic_cast<TFile *>(file))
               out << "   " << chName.Data() << "->AddFriend(\"" << frtree->GetName() << "\", \"" << file->GetName()
                   << "\");" << std::endl;
         }
      }
   }
   out << std::endl;

   for (TObject *el : *fFiles) {
      TChainElement *chel = (TChainElement *)el;
      // Save tree file if it is really loaded to the chain
      if (chel->GetLoadResult() == 0 && chel->GetEntries() != 0) {
         if (chel->GetEntries() == TTree::kMaxEntries) // tree number of entries is not yet known
            out << "   " << chName.Data() << "->AddFile(\"" << chel->GetTitle() << "\");" << std::endl;
         else
            out << "   " << chName.Data() << "->AddFile(\"" << chel->GetTitle() << "\"," << chel->GetEntries() << ");"
                << std::endl;
      }
   }
   out << std::endl;

   if (GetMarkerColor() != 1) {
      if (GetMarkerColor() > 228) {
         TColor::SaveColor(out, GetMarkerColor());
         out << "   " << chName.Data() << "->SetMarkerColor(ci);" << std::endl;
      } else
         out << "   " << chName.Data() << "->SetMarkerColor(" << GetMarkerColor() << ");" << std::endl;
   }
   if (GetMarkerStyle() != 1) {
      out << "   " << chName.Data() << "->SetMarkerStyle(" << GetMarkerStyle() << ");" << std::endl;
   }
   if (GetMarkerSize() != 1) {
      out << "   " << chName.Data() << "->SetMarkerSize(" << GetMarkerSize() << ");" << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Loop on tree and print entries passing selection.
/// - If varexp is 0 (or "") then print only first 8 columns.
/// - If varexp = "*" print all columns.
/// - Otherwise a columns selection can be made using "var1:var2:var3".
/// See TTreePlayer::Scan for more information.

Long64_t TChain::Scan(const char* varexp, const char* selection, Option_t* option, Long64_t nentries, Long64_t firstentry)
{
   if (LoadTree(firstentry) < 0) {
      return 0;
   }
   return TTree::Scan(varexp, selection, option, nentries, firstentry);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the global branch kAutoDelete bit.
///
/// When LoadTree loads a new Tree, the branches for which
/// the address is set will have the option AutoDelete set
/// For more details on AutoDelete, see TBranch::SetAutoDelete.

void TChain::SetAutoDelete(Bool_t autodelete)
{
   if (autodelete) {
      SetBit(kAutoDelete, 1);
   } else {
      SetBit(kAutoDelete, 0);
   }
}

Int_t TChain::SetCacheSize(Long64_t cacheSize)
{
   // Set the cache size of the underlying TTree,
   // See TTree::SetCacheSize.
   // Returns  0 cache state ok (exists or not, as appropriate)
   //         -1 on error

   Int_t res = 0;

   // remember user has requested this cache setting
   fCacheUserSet = kTRUE;

   if (fTree) {
      res = fTree->SetCacheSize(cacheSize);
   } else {
      // If we don't have a TTree yet only record the cache size wanted
      res = 0;
   }
   fCacheSize = cacheSize; // Record requested size.
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the addresses of the branch.

void TChain::ResetBranchAddress(TBranch *branch)
{
   TChainElement* element = (TChainElement*) fStatus->FindObject(branch->GetName());
   if (element) {
      element->SetBaddress(0);
   }
   if (fTree) {
      fTree->ResetBranchAddress(branch);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the addresses of the branches.

void TChain::ResetBranchAddresses()
{
   TIter next(fStatus);
   TChainElement* element = 0;
   while ((element = (TChainElement*) next())) {
      element->SetBaddress(0);
   }
   if (fTree) {
      fTree->ResetBranchAddresses();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set branch address.
///
/// \param[in] bname    is the name of a branch.
/// \param[in] add      is the address of the branch.
/// \param[in] ptr
///
/// Note: See the comments in TBranchElement::SetAddress() for a more
/// detailed discussion of the meaning of the add parameter.
///
/// IMPORTANT REMARK:
///
/// In case TChain::SetBranchStatus is called, it must be called
/// BEFORE calling this function.
///
/// See TTree::CheckBranchAddressType for the semantic of the return value.

Int_t TChain::SetBranchAddress(const char *bname, void* add, TBranch** ptr)
{
   Int_t res = kNoCheck;

   // Check if bname is already in the status list.
   // If not, create a TChainElement object and set its address.
   TChainElement* element = (TChainElement*) fStatus->FindObject(bname);
   if (!element) {
      element = new TChainElement(bname, "");
      fStatus->Add(element);
   }
   element->SetBaddress(add);
   element->SetBranchPtr(ptr);
   // Also set address in current tree.
   // FIXME: What about the chain clones?
   if (fTreeNumber >= 0) {
      TBranch* branch = fTree->GetBranch(bname);
      if (ptr) {
         *ptr = branch;
      }
      if (branch) {
         res = CheckBranchAddressType(branch, TClass::GetClass(element->GetBaddressClassName()), (EDataType) element->GetBaddressType(), element->GetBaddressIsPtr());
         if ((res & kNeedEnableDecomposedObj) && !branch->GetMakeClass()) {
            branch->SetMakeClass(kTRUE);
         }
         element->SetDecomposedObj(branch->GetMakeClass());
         element->SetCheckedType(kTRUE);
         if (fClones) {
            void* oldAdd = branch->GetAddress();
            for (TObjLink* lnk = fClones->FirstLink(); lnk; lnk = lnk->Next()) {
               TTree* clone = (TTree*) lnk->GetObject();
               TBranch* cloneBr = clone->GetBranch(bname);
               if (cloneBr && (cloneBr->GetAddress() == oldAdd)) {
                  // the clone's branch is still pointing to us
                  cloneBr->SetAddress(add);
                  if ((res & kNeedEnableDecomposedObj) && !cloneBr->GetMakeClass()) {
                     cloneBr->SetMakeClass(kTRUE);
                  }
               }
            }
         }

         branch->SetAddress(add);
      } else {
         Error("SetBranchAddress", "unknown branch -> %s", bname);
         return kMissingBranch;
      }
   } else {
      if (ptr) {
         *ptr = 0;
      }
   }
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if bname is already in the status list, and if not, create a TChainElement object and set its address.
/// See TTree::CheckBranchAddressType for the semantic of the return value.
///
/// Note: See the comments in TBranchElement::SetAddress() for a more
/// detailed discussion of the meaning of the add parameter.

Int_t TChain::SetBranchAddress(const char* bname, void* add, TClass* realClass, EDataType datatype, Bool_t isptr)
{
   return SetBranchAddress(bname, add, 0, realClass, datatype, isptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if bname is already in the status list, and if not, create a TChainElement object and set its address.
/// See TTree::CheckBranchAddressType for the semantic of the return value.
///
/// Note: See the comments in TBranchElement::SetAddress() for a more
/// detailed discussion of the meaning of the add parameter.

Int_t TChain::SetBranchAddress(const char* bname, void* add, TBranch** ptr, TClass* realClass, EDataType datatype, Bool_t isptr)
{
   TChainElement* element = (TChainElement*) fStatus->FindObject(bname);
   if (!element) {
      element = new TChainElement(bname, "");
      fStatus->Add(element);
   }
   if (realClass) {
      element->SetBaddressClassName(realClass->GetName());
   }
   element->SetBaddressType((UInt_t) datatype);
   element->SetBaddressIsPtr(isptr);
   element->SetBranchPtr(ptr);
   return SetBranchAddress(bname, add, ptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Set branch status to Process or DoNotProcess
///
/// \param[in] bname     is the name of a branch. if bname="*", apply to all branches.
/// \param[in] status    = 1  branch will be processed,
///                      = 0  branch will not be processed
/// \param[out] found
///
///  See IMPORTANT REMARKS in TTree::SetBranchStatus and TChain::SetBranchAddress
///
///  If found is not 0, the number of branch(es) found matching the regular
///  expression is returned in *found AND the error message 'unknown branch'
///  is suppressed.

void TChain::SetBranchStatus(const char* bname, Bool_t status, UInt_t* found)
{
   // FIXME: We never explicitly set found to zero!

   // Check if bname is already in the status list,
   // if not create a TChainElement object and set its status.
   TChainElement* element = (TChainElement*) fStatus->FindObject(bname);
   if (element) {
      fStatus->Remove(element);
   } else {
      element = new TChainElement(bname, "");
   }
   fStatus->Add(element);
   element->SetStatus(status);
   // Also set status in current tree.
   if (fTreeNumber >= 0) {
      fTree->SetBranchStatus(bname, status, found);
   } else if (found) {
      *found = 1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove reference to this chain from current directory and add
/// reference to new directory dir. dir can be 0 in which case the chain
/// does not belong to any directory.

void TChain::SetDirectory(TDirectory* dir)
{
   if (fDirectory == dir) return;
   if (fDirectory) fDirectory->Remove(this);
   fDirectory = dir;
   if (fDirectory) {
      fDirectory->Append(this);
      fFile = fDirectory->GetFile();
   } else {
      fFile = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the input entry list (processing the entries of the chain will then be
/// limited to the entries in the list).
/// This function finds correspondence between the sub-lists of the TEntryList
/// and the trees of the TChain.
/// By default (opt=""), both the file names of the chain elements and
/// the file names of the TEntryList sublists are expanded to full path name.
/// If opt = "ne", the file names are taken as they are and not expanded

void TChain::SetEntryList(TEntryList *elist, Option_t *opt)
{
   if (fEntryList){
      //check, if the chain is the owner of the previous entry list
      //(it happens, if the previous entry list was created from a user-defined
      //TEventList in SetEventList() function)
      if (fEntryList->TestBit(kCanDelete)) {
         TEntryList *tmp = fEntryList;
         fEntryList = 0; // Avoid problem with RecursiveRemove.
         delete tmp;
      } else {
         fEntryList = 0;
      }
   }
   if (!elist){
      fEntryList = 0;
      fEventList = 0;
      return;
   }
   if (!elist->TestBit(kCanDelete)){
      //this is a direct call to SetEntryList, not via SetEventList
      fEventList = 0;
   }
   if (elist->GetN() == 0){
      fEntryList = elist;
      return;
   }
   if (fProofChain){
      //for processing on proof, event list and entry list can't be
      //set at the same time.
      fEventList = 0;
      fEntryList = elist;
      return;
   }

   Int_t ne = fFiles->GetEntries();
   Int_t listfound=0;
   TString treename, filename;

   TEntryList *templist = 0;
   for (Int_t ie = 0; ie<ne; ie++){
      auto chainElement = (TChainElement*)fFiles->UncheckedAt(ie);
      treename = chainElement->GetName();
      filename = chainElement->GetTitle();
      templist = elist->GetEntryList(treename, filename, opt);
      if (templist) {
         listfound++;
         templist->SetTreeNumber(ie);
      }
   }

   if (listfound == 0){
      Error("SetEntryList", "No list found for the trees in this chain");
      fEntryList = 0;
      return;
   }
   fEntryList = elist;
   TList *elists = elist->GetLists();
   Bool_t shift = kFALSE;
   TIter next(elists);

   //check, if there are sub-lists in the entry list, that don't
   //correspond to any trees in the chain
   while((templist = (TEntryList*)next())){
      if (templist->GetTreeNumber() < 0){
         shift = kTRUE;
         break;
      }
   }
   fEntryList->SetShift(shift);

}

////////////////////////////////////////////////////////////////////////////////
/// Set the input entry list (processing the entries of the chain will then be
/// limited to the entries in the list). This function creates a special kind
/// of entry list (TEntryListFromFile object) that loads lists, corresponding
/// to the chain elements, one by one, so that only one list is in memory at a time.
///
/// If there is an error opening one of the files, this file is skipped and the
/// next file is loaded
///
/// File naming convention:
///
/// - by default, filename_elist.root is used, where filename is the
///   name of the chain element
/// - xxx$xxx.root - $ sign is replaced by the name of the chain element
///
/// If the list name is not specified (by passing filename_elist.root/listname to
/// the TChain::SetEntryList() function, the first object of class TEntryList
/// in the file is taken.
///
/// It is assumed, that there are as many list files, as there are elements in
/// the chain and they are in the same order

void TChain::SetEntryListFile(const char *filename, Option_t * /*opt*/)
{

   if (fEntryList){
      //check, if the chain is the owner of the previous entry list
      //(it happens, if the previous entry list was created from a user-defined
      //TEventList in SetEventList() function)
      if (fEntryList->TestBit(kCanDelete)) {
         TEntryList *tmp = fEntryList;
         fEntryList = 0; // Avoid problem with RecursiveRemove.
         delete tmp;
      } else {
         fEntryList = 0;
      }
   }

   fEventList = 0;

   TString basename(filename);

   Int_t dotslashpos = basename.Index(".root/");
   TString behind_dot_root = "";
   if (dotslashpos>=0) {
      // Copy the list name specification
      behind_dot_root = basename(dotslashpos+6,basename.Length()-dotslashpos+6);
      // and remove it from basename
      basename.Remove(dotslashpos+5);
   }
   fEntryList = new TEntryListFromFile(basename.Data(), behind_dot_root.Data(), fNtrees);
   fEntryList->SetBit(kCanDelete, kTRUE);
   fEntryList->SetDirectory(0);
   ((TEntryListFromFile*)fEntryList)->SetFileNames(fFiles);
}

////////////////////////////////////////////////////////////////////////////////
/// This function transfroms the given TEventList into a TEntryList
///
/// NOTE, that this function loads all tree headers, because the entry numbers
/// in the TEventList are global and have to be recomputed, taking into account
/// the number of entries in each tree.
///
/// The new TEntryList is owned by the TChain and gets deleted when the chain
/// is deleted. This TEntryList is returned by GetEntryList() function, and after
/// GetEntryList() function is called, the TEntryList is not owned by the chain
/// any more and will not be deleted with it.

void TChain::SetEventList(TEventList *evlist)
{
   fEventList = evlist;
   if (fEntryList) {
      if (fEntryList->TestBit(kCanDelete)) {
         TEntryList *tmp = fEntryList;
         fEntryList = 0; // Avoid problem with RecursiveRemove.
         delete tmp;
      } else {
         fEntryList = 0;
      }
   }

   if (!evlist) {
      fEntryList = 0;
      fEventList = 0;
      return;
   }

   if(fProofChain) {
      //on proof, fEventList and fEntryList shouldn't be set at the same time
      if (fEntryList){
         //check, if the chain is the owner of the previous entry list
         //(it happens, if the previous entry list was created from a user-defined
         //TEventList in SetEventList() function)
         if (fEntryList->TestBit(kCanDelete)){
            TEntryList *tmp = fEntryList;
            fEntryList = 0; // Avoid problem with RecursiveRemove.
            delete tmp;
         } else {
            fEntryList = 0;
         }
      }
      return;
   }

   char enlistname[100];
   snprintf(enlistname,100, "%s_%s", evlist->GetName(), "entrylist");
   TEntryList *enlist = new TEntryList(enlistname, evlist->GetTitle());
   enlist->SetDirectory(0);

   Int_t nsel = evlist->GetN();
   Long64_t globalentry, localentry;
   const char *treename;
   const char *filename;
   if (fTreeOffset[fNtrees-1]==TTree::kMaxEntries){
      //Load all the tree headers if the tree offsets are not known
      //It is assumed here, that loading the last tree will load all
      //previous ones
      printf("loading trees\n");
      (const_cast<TChain*>(this))->LoadTree(evlist->GetEntry(evlist->GetN()-1));
   }
   for (Int_t i=0; i<nsel; i++){
      globalentry = evlist->GetEntry(i);
      //add some protection from globalentry<0 here
      Int_t treenum = 0;
      while (globalentry>=fTreeOffset[treenum])
         treenum++;
      treenum--;
      localentry = globalentry - fTreeOffset[treenum];
      // printf("globalentry=%lld, treeoffset=%lld, localentry=%lld\n", globalentry, fTreeOffset[treenum], localentry);
      treename = ((TNamed*)fFiles->At(treenum))->GetName();
      filename = ((TNamed*)fFiles->At(treenum))->GetTitle();
      //printf("entering for tree %s %s\n", treename, filename);
      enlist->SetTree(treename, filename);
      enlist->Enter(localentry);
   }
   enlist->SetBit(kCanDelete, kTRUE);
   enlist->SetReapplyCut(evlist->GetReapplyCut());
   SetEntryList(enlist);
}

////////////////////////////////////////////////////////////////////////////////
/// Change the name of this TChain.

void TChain::SetName(const char* name)
{
   {
      // Should this be extended to include the call to TTree::SetName?
      R__WRITE_LOCKGUARD(ROOT::gCoreMutex); // Take the lock once rather than 3 times.
      gROOT->GetListOfCleanups()->Remove(this);
      gROOT->GetListOfSpecials()->Remove(this);
      gROOT->GetListOfDataSets()->Remove(this);
   }
   TTree::SetName(name);
   {
      // Should this be extended to include the call to TTree::SetName?
      R__WRITE_LOCKGUARD(ROOT::gCoreMutex); // Take the lock once rather than 3 times.
      gROOT->GetListOfCleanups()->Add(this);
      gROOT->GetListOfSpecials()->Add(this);
      gROOT->GetListOfDataSets()->Add(this);
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Set number of entries per packet for parallel root.

void TChain::SetPacketSize(Int_t size)
{
   fPacketSize = size;
   TIter next(fFiles);
   TChainElement *element;
   while ((element = (TChainElement*)next())) {
      element->SetPacketSize(size);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/Disable PROOF processing on the current default Proof (gProof).
///
/// "Draw" and "Processed" commands will be handled by PROOF.
/// The refresh and gettreeheader are meaningful only if on == kTRUE.
/// If refresh is kTRUE the underlying fProofChain (chain proxy) is always
/// rebuilt (even if already existing).
/// If gettreeheader is kTRUE the header of the tree will be read from the
/// PROOF cluster: this is only needed for browsing and should be used with
/// care because it may take a long time to execute.

void TChain::SetProof(Bool_t on, Bool_t refresh, Bool_t gettreeheader)
{
   if (!on) {
      // Disable
      SafeDelete(fProofChain);
      // Reset related bit
      ResetBit(kProofUptodate);
   } else {
      if (fProofChain && !refresh &&
         (!gettreeheader || (gettreeheader && fProofChain->GetTree()))) {
         return;
      }
      SafeDelete(fProofChain);
      ResetBit(kProofUptodate);

      // Make instance of TChainProof via the plugin manager
      TPluginHandler *h;
      if ((h = gROOT->GetPluginManager()->FindHandler("TChain", "proof"))) {
         if (h->LoadPlugin() == -1)
            return;
         if (!(fProofChain = reinterpret_cast<TChain *>(h->ExecPlugin(2, this, gettreeheader))))
            Error("SetProof", "creation of TProofChain failed");
         // Set related bits
         SetBit(kProofUptodate);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set chain weight.
///
/// The weight is used by TTree::Draw to automatically weight each
/// selected entry in the resulting histogram.
/// For example the equivalent of
/// ~~~ {.cpp}
///     chain.Draw("x","w")
/// ~~~
/// is
/// ~~~ {.cpp}
///     chain.SetWeight(w,"global");
///     chain.Draw("x");
/// ~~~
/// By default the weight used will be the weight
/// of each Tree in the TChain. However, one can force the individual
/// weights to be ignored by specifying the option "global".
/// In this case, the TChain global weight will be used for all Trees.

void TChain::SetWeight(Double_t w, Option_t* option)
{
   fWeight = w;
   TString opt = option;
   opt.ToLower();
   ResetBit(kGlobalWeight);
   if (opt.Contains("global")) {
      SetBit(kGlobalWeight);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.

void TChain::Streamer(TBuffer& b)
{
   if (b.IsReading()) {
      // Remove using the 'old' name.
      {
         R__LOCKGUARD(gROOTMutex);
         gROOT->GetListOfCleanups()->Remove(this);
      }

      UInt_t R__s, R__c;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         b.ReadClassBuffer(TChain::Class(), this, R__v, R__s, R__c);
      } else {
         //====process old versions before automatic schema evolution
         TTree::Streamer(b);
         b >> fTreeOffsetLen;
         b >> fNtrees;
         fFiles->Streamer(b);
         if (R__v > 1) {
            fStatus->Streamer(b);
            fTreeOffset = new Long64_t[fTreeOffsetLen];
            b.ReadFastArray(fTreeOffset,fTreeOffsetLen);
         }
         b.CheckByteCount(R__s, R__c, TChain::IsA());
         //====end of old versions
      }
      // Re-add using the new name.
      {
         R__LOCKGUARD(gROOTMutex);
         gROOT->GetListOfCleanups()->Add(this);
      }

   } else {
      b.WriteClassBuffer(TChain::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Dummy function kept for back compatibility.
/// The cache is now activated automatically when processing TTrees/TChain.

void TChain::UseCache(Int_t /* maxCacheSize */, Int_t /* pageSize */)
{
}
