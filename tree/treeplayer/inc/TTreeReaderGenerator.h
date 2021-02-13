// @(#)root/treeplayer:$Id$
// Author: Akos Hajdu 22/06/2015

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers and al.        *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreeReaderGenerator
#define ROOT_TTreeReaderGenerator

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TTreeReaderGenerator                                                 //
//                                                                      //
// Generate a Selector using the TTreeReader interface                  //
// (TTreeReaderValue, TTreeReaderArray) to access the data in the tree. //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TTreeGeneratorBase.h"

#include "TNamed.h"
#include <vector>

class TBranch;
class TBranchElement;
class TLeaf;

namespace ROOT {
namespace Internal {

   // 0 for the general case, 1 when this a split clases inside a TClonesArray,
   // 2 when this is a split classes inside an STL container.
   enum ELocation { kOut=0, kClones, kSTL };

   class TTreeReaderDescriptor : public TObject {
   public:
      enum ReaderType { kValue, kArray };
      ReaderType fType;    // Type of the reader: Value or Array
      TString fDataType;   // Data type of reader
      TString fName;       // Reader name
      TString fBranchName; // Branch corresponding to the reader

      TTreeReaderDescriptor(ReaderType type, TString dataType, TString name, TString branchName) :
         fType(type),
         fDataType(dataType),
         fName(name),
         fBranchName(branchName) { }
   };

   class TBranchDescriptor : public TNamed {
   public:
      ELocation             fIsClones;       // Type of container
      TString               fContainerName;  // Name of the container
      TString               fBranchName;     // Name of the branch
      TString               fSubBranchPrefix;// Prefix (e.g. if the branch name is "A." the prefix is "A"
      TVirtualStreamerInfo *fInfo;           // Streamer info
      TBranchDescriptor    *fParent;         // Descriptor of the parent branch (NULL for topmost)

      TBranchDescriptor(const char *type, TVirtualStreamerInfo *info,
                        const char *branchname, const char *subBranchPrefix, ELocation isclones,
                        const TString &containerName, TBranchDescriptor *parent = 0) :
         TNamed(type,type),
         fIsClones(isclones),
         fContainerName(containerName),
         fBranchName(branchname),
         fSubBranchPrefix(subBranchPrefix),
         fInfo(info),
         fParent(parent)
         {
            if (fSubBranchPrefix.Length() && fSubBranchPrefix[fSubBranchPrefix.Length() - 1] == '.') {
               fSubBranchPrefix.Remove(fSubBranchPrefix.Length()-1);
            }
         }

      Bool_t IsClones() const { return fIsClones == kClones; }

      Bool_t IsSTL() const { return fIsClones == kSTL; }
   };

   class TTreeReaderGenerator : public TTreeGeneratorBase
   {
      TString               fClassname;         // Class name of the selector
      TList                 fListOfReaders;     // List of readers
      Bool_t                fIncludeAllLeaves;  // Should all leaves be included
      Bool_t                fIncludeAllTopmost; // Should all topmost branches be included
      std::vector<TString>  fIncludeLeaves;     // Branches whose leaves should be included
      std::vector<TString>  fIncludeStruct;     // Branches whom should be included

      void   AddReader(TTreeReaderDescriptor::ReaderType type, TString dataType, TString name,
                       TString branchName, TBranchDescriptor *parent = 0, Bool_t isLeaf = kTRUE);
      UInt_t AnalyzeBranches(TBranchDescriptor *desc, TBranchElement *branch, TVirtualStreamerInfo *info);
      UInt_t AnalyzeBranches(TBranchDescriptor *desc, TIter &branches, TVirtualStreamerInfo *info);
      UInt_t AnalyzeOldBranch(TBranch *branch);
      UInt_t AnalyzeOldLeaf(TLeaf *leaf, Int_t nleaves);
      Bool_t BranchNeedsReader(TString branchName, TBranchDescriptor *parent, Bool_t isLeaf);

      void   ParseOptions();
      void   AnalyzeTree(TTree *tree);
      void   WriteSelector();

   public:
      TTreeReaderGenerator(TTree* tree, const char *classname, Option_t *option);
   };
}
}

#endif
