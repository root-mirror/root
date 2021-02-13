// @(#)root/treeplayer:$Id$
// Author: Akos Hajdu 22/06/2015

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers and al.        *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TTreeReaderGenerator.h"
#include <algorithm>
#include <stdio.h>
#include <fstream>

#include "TBranchElement.h"
#include "TChain.h"
#include "TClass.h"
#include "TClassEdit.h"
#include "TClonesArray.h"
#include "TError.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TLeafC.h"
#include "TLeafObject.h"
#include "TROOT.h"
#include "TStreamerElement.h"
#include "TStreamerInfo.h"
#include "TTree.h"
#include "TObjString.h"
#include "TVirtualCollectionProxy.h"
#include "TVirtualStreamerInfo.h"

namespace ROOT {
namespace Internal {

   ////////////////////////////////////////////////////////////////////////////////
   /// Constructor. Analyzes the tree and writes selector.

   TTreeReaderGenerator::TTreeReaderGenerator(TTree* tree, const char *classname, Option_t *option) :
      TTreeGeneratorBase(tree, option), fClassname(classname),
      fIncludeAllLeaves(kFALSE), fIncludeAllTopmost(kFALSE)
   {
      ParseOptions();
      AnalyzeTree(fTree);
      WriteSelector();
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Add a reader to the generated code.

   void TTreeReaderGenerator::AddReader(TTreeReaderDescriptor::ReaderType type, TString dataType, TString name,
                                        TString branchName, TBranchDescriptor *parent, Bool_t isLeaf)
   {
      if(BranchNeedsReader(branchName, parent, isLeaf)) {
         // Ignore unknown types
         if (dataType.EqualTo("")) {
            Warning("TTreeReaderGenerator::AddReader", "Ignored branch %s because type is unsupported.", branchName.Data());
            return;
         }
         // Loop through existing readers to check duplicate branch names
         TIter next(&fListOfReaders);
         TTreeReaderDescriptor *descriptor;
         while ( ( descriptor = (TTreeReaderDescriptor*)next() ) ) {
            if (descriptor->fBranchName.EqualTo(branchName)) {
               Warning("AddReader", "Ignored branch %s because a branch with the same name already exists. "
                                    "TTreeReader requires an unique name for the branches. You may need to "
                                    "put a dot at the end of the name of top-level branches.", branchName.Data());
               return;
            }
         }
         name.ReplaceAll('.', '_'); // Replace dots with underscore
         // Remove array dimensions from name
         while (name.Index('[') >= 0 && name.Index(']') >= 0 && name.Index(']') > name.Index('[')) {
            name.Remove(name.Index('['), name.Index(']') - name.Index('[') + 1);
         }
         fListOfReaders.Add( new TTreeReaderDescriptor(type, dataType, name, branchName) );
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Analyse sub-branches of 'branch' recursively and extract readers.

   UInt_t TTreeReaderGenerator::AnalyzeBranches(TBranchDescriptor *desc, TBranchElement *branch, TVirtualStreamerInfo *info)
   {
      if (info==0) info = branch->GetInfo();

      TIter branches(branch->GetListOfBranches());

      return AnalyzeBranches(desc, branches, info);
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Analyse sub-branches 'branches' recursively and extract readers.

   UInt_t TTreeReaderGenerator::AnalyzeBranches(TBranchDescriptor *desc, TIter &branches, TVirtualStreamerInfo *info)
   {
      UInt_t lookedAt = 0;     // Number of sub-branches analyzed
      ELocation outer_isclones = kOut; // Is the parent branch a container
      TString containerName;   // Container name
      TString subBranchPrefix; // Prefix of sub-branch (if the elements and sub-branches do not match).
      Bool_t skipped = false;  // Should the branch be skipped

      // Check for containers (TClonesArray or STL)
      {
         TIter peek = branches;
         TBranchElement *branch = (TBranchElement*)peek();
         if (desc && desc->IsClones()) {
            outer_isclones = kClones;
            containerName = "TClonesArray";
         } else if (desc && desc->IsSTL()) {
            outer_isclones = kSTL;
            containerName = desc->fContainerName;
         } else if (!desc && branch && branch->GetBranchCount() == branch->GetMother()) {
            if ( ((TBranchElement*)(branch->GetMother()))->GetType()==3)  {
               outer_isclones = kClones;
               containerName = "TClonesArray";
            } else {
               outer_isclones = kSTL;
               containerName = branch->GetMother()->GetClassName();
            }
         }
         // FIXME: this is wrong because 'branch' is the first sub-branch and even though
         // it can be a collection, it does not mean that the other sub-branches are
         // collections as well.
         /* else if (branch->GetType() == 3) {
            outer_isclones = kClones;
            containerName = "TClonesArray";
         } else if (branch->GetType() == 4) {
            outer_isclones = kSTL;
            containerName = branch->GetMother()->GetSubBranch(branch)->GetClassName();
         }*/
         if (desc) {
            subBranchPrefix = desc->fSubBranchPrefix;
         } else {
            TBranchElement *mom = (TBranchElement*)branch->GetMother();
            subBranchPrefix = mom->GetName();
            if (subBranchPrefix[subBranchPrefix.Length()-1]=='.') {
               subBranchPrefix.Remove(subBranchPrefix.Length()-1);
            } else if (mom->GetType()!=3 && mom->GetType() != 4) {
               subBranchPrefix = "";
            }
         }
      }

      // Loop through elements (i.e., sub-branches). The problem is that the elements
      // and sub-branches do not always match. For example suppose that class A contains
      // members x and y, and there is a class B inheriting from A and containing an extra
      // member z. It is possible that B has two elements: A and z but three sub-branches
      // x, y and z. Therefore, the branch iterator is treated differently.
      TIter elements( info->GetElements() );
      for( TStreamerElement *element = (TStreamerElement*)elements();
           element;
           element = (TStreamerElement*)elements() )
      {
         Bool_t isBase = false;     // Does the element correspond to a base class
         Bool_t usedBranch = kTRUE; // Does the branch correspond to the element (i.e., they match)
         Bool_t isLeaf = true;    // Is the branch a leaf (i.e. no sub-branches)
         TIter peek = branches;     // Iterator for sub-branches
         // Always start with the first available sub-branch and if it does not match the element,
         // try the next ones
         TBranchElement *branch = (TBranchElement*)peek();
         // There is a problem if there are more elements than branches
         if (branch==0) {
            if (desc) {
               Error("AnalyzeBranches","Ran out of branches when looking in branch %s, class %s",
                     desc->fBranchName.Data(), info->GetName());
            } else {
               Error("AnalyzeBranches","Ran out of branches when looking in class %s, element %s",
                     info->GetName(), element->GetName());
            }
            return lookedAt;
         }

         if (info->GetClass()->GetCollectionProxy() && strcmp(element->GetName(),"This")==0) {
            continue; // Skip the artifical streamer element.
         }

         if (element->GetType() == -1) {
            continue; // This is an ignored TObject base class.
         }

         // Get branch name
         TString branchname = branch->GetName();
         TString branchEndName;
         {
            TLeaf *leaf = (TLeaf*)branch->GetListOfLeaves()->At(0);
            if (leaf && outer_isclones == kOut
                && !(branch->GetType() == 3 || branch->GetType() == 4)) branchEndName = leaf->GetName();
            else branchEndName = branch->GetName();
            Int_t pos;
            pos = branchEndName.Index(".");
            if (pos!=-1) {
               if (subBranchPrefix.Length() && branchEndName.BeginsWith(subBranchPrefix)) {
                  branchEndName.Remove(0, subBranchPrefix.Length() + 1);
               }
            }
         }

         TString dataType; // Data type of reader
         TTreeReaderDescriptor::ReaderType readerType = TTreeReaderDescriptor::ReaderType::kValue;
         Bool_t ispointer = false;
         ELocation isclones = outer_isclones; // Is the actual sub-branch a collection (inherit from parent branch)
         // Get data type
         switch(element->GetType()) {
            // Built-in types
            case TVirtualStreamerInfo::kBool:    { dataType = "Bool_t";         readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kChar:    { dataType = "Char_t";         readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kShort:   { dataType = "Short_t";        readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kInt:     { dataType = "Int_t";          readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kLong:    { dataType = "Long_t";         readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kLong64:  { dataType = "Long64_t";       readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kFloat:   { dataType = "Float_t";        readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kFloat16: { dataType = "Float16_t";      readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kDouble:  { dataType = "Double_t";       readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kDouble32:{ dataType = "Double32_t";     readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kUChar:   { dataType = "UChar_t";        readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kUShort:  { dataType = "unsigned short"; readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kUInt:    { dataType = "unsigned int";   readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kULong:   { dataType = "ULong_t";        readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kULong64: { dataType = "ULong64_t";      readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            case TVirtualStreamerInfo::kBits:    { dataType = "unsigned int";   readerType = TTreeReaderDescriptor::ReaderType::kValue; break; }
            // Character arrays
            case TVirtualStreamerInfo::kCharStar: { dataType = "Char_t"; readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            // Array of built-in types [8]
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kBool:    { dataType = "Bool_t";         readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kChar:    { dataType = "Char_t";         readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kShort:   { dataType = "Short_t";        readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kInt:     { dataType = "Int_t";          readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kLong:    { dataType = "Long_t";         readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kLong64:  { dataType = "Long64_t";       readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kFloat:   { dataType = "Float_t";        readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kFloat16: { dataType = "Float16_t";      readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kDouble:  { dataType = "Double_t";       readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kDouble32:{ dataType = "Double32_t";     readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kUChar:   { dataType = "UChar_t";        readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kUShort:  { dataType = "unsigned short"; readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kUInt:    { dataType = "unsigned int";   readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kULong:   { dataType = "ULong_t";        readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kULong64: { dataType = "ULong64_t";      readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kBits:    { dataType = "unsigned int";   readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            // Array of built-in types [n]
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kBool:    { dataType = "Bool_t";         readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kChar:    { dataType = "Char_t";         readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kShort:   { dataType = "Short_t";        readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kInt:     { dataType = "Int_t";          readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kLong:    { dataType = "Long_t";         readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kLong64:  { dataType = "Long64_t";       readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kFloat:   { dataType = "Float_t";        readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kFloat16: { dataType = "Float16_t";      readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kDouble:  { dataType = "Double_t";       readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kDouble32:{ dataType = "Double32_t";     readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kUChar:   { dataType = "UChar_t";        readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kUShort:  { dataType = "unsigned short"; readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kUInt:    { dataType = "unsigned int";   readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kULong:   { dataType = "ULong_t";        readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kULong64: { dataType = "ULong64_t";      readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            case TVirtualStreamerInfo::kOffsetP + TVirtualStreamerInfo::kBits:    { dataType = "unsigned int";   readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            // array counter [n]
            case TVirtualStreamerInfo::kCounter: { dataType = "Int_t"; readerType = TTreeReaderDescriptor::ReaderType::kArray; break; }
            // other stuff (containers, classes, ...)
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kObjectp:
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kObjectP:
            case TVirtualStreamerInfo::kObjectp:
            case TVirtualStreamerInfo::kObjectP:
            case TVirtualStreamerInfo::kAnyp:
            case TVirtualStreamerInfo::kAnyP:
            case TVirtualStreamerInfo::kSTL + TVirtualStreamerInfo::kObjectp:
            case TVirtualStreamerInfo::kSTL + TVirtualStreamerInfo::kObjectP:
            // set as pointers and fall through to the next switches
               ispointer = true;
            case TVirtualStreamerInfo::kOffsetL + TVirtualStreamerInfo::kObject:
               // This means an array of objects, but then fall through
               readerType = TTreeReaderDescriptor::ReaderType::kArray;
            case TVirtualStreamerInfo::kObject:
            case TVirtualStreamerInfo::kTString:
            case TVirtualStreamerInfo::kTNamed:
            case TVirtualStreamerInfo::kTObject:
            case TVirtualStreamerInfo::kAny:
            case TVirtualStreamerInfo::kBase:
            case TVirtualStreamerInfo::kSTL: {
               TClass *cl = element->GetClassPointer();
               R__ASSERT(cl);
               dataType = cl->GetName();
               // Check for containers
               if (cl == TClonesArray::Class()) { // TClonesArray
                  isclones = kClones;
                  containerName = "TClonesArray";
                  if (outer_isclones != kOut) { // If the parent is already a collection
                     isclones = outer_isclones;
                     dataType = "TClonesArray";
                  } else {
                     readerType = TTreeReaderDescriptor::ReaderType::kArray;
                     dataType = GetContainedClassName(branch, element, ispointer);
                  }
               } else if (cl->GetCollectionProxy()) { // STL collection
                  isclones = kSTL;
                  containerName = cl->GetName();
                  if (outer_isclones != kOut || containerName.EqualTo("vector<bool>")) {
                     // If the parent is already a collection we can only add this collection as a whole
                     // Also TTreeReaderArray does currently not support vectors of bool so
                     // that need to be added as a whole.
                     // Also getting the inner type of vector would return "unsigned" so the full
                     // name has to be compared
                     isclones = outer_isclones;
                     dataType = cl->GetName();
                  } else {
                     readerType = TTreeReaderDescriptor::ReaderType::kArray;
                     TClass *valueClass = cl->GetCollectionProxy()->GetValueClass();
                     if (valueClass) { // Get class inside container
                        dataType = valueClass->GetName();
                     }
                     else { // Get built-in type inside container
                        TDataType *valueClassBuiltIn = TDataType::GetDataType(cl->GetCollectionProxy()->GetType());
                        if (valueClassBuiltIn) dataType = valueClassBuiltIn->GetName();
                        else Error("AnalyzeBranches", "Could not get type from collection %s in branch %s", cl->GetName(), branch->GetName());
                     }
                  }
               }

               // Analyze the actual element
               // Possible cases:
               // - Base class
               //   - Element and branch matches
               //     - Non-split: do nothing
               //     - Split: recurse with sub-sub-branches
               //   - Element and branch does not match: recurse with same branches
               // - Not base class
               //   - Element and branch matches
               //     - Non-split: do nothing
               //     - Split: recurse with sub-sub-branches
               //   - Element and branch does not match: recurse with same branches
               TBranch *parent = branch->GetMother()->GetSubBranch(branch);
               TVirtualStreamerInfo *objInfo = 0;
               if (branch->GetListOfBranches()->GetEntries()) {
                  objInfo = ((TBranchElement*)branch->GetListOfBranches()->At(0))->GetInfo();
               } else {
                  objInfo = branch->GetInfo();
               }
               if (element->IsBase()) { // Base class
                  isBase = true;
                  if (cl == TObject::Class() && info->GetClass()->CanIgnoreTObjectStreamer())
                  {
                     continue; // Ignore TObject
                  }

                  TBranchDescriptor *bdesc = 0;

                  if (branchEndName == element->GetName()) { // The element and the branch matches
                     if (branch->GetListOfBranches()->GetEntries() == 0) { // The branch contains a non-split base class
                        // FIXME: nothing to do in such cases, because readers cannot access
                        // non-split members and a reader for the whole branch will be added
                     } else { // The branch contains a split base class
                        Int_t pos = branchname.Last('.');
                        if (pos != -1) {
                           branchname.Remove(pos);
                        }
                        TString local_prefix = desc ? desc->fSubBranchPrefix : TString(parent->GetName());
                        bdesc = new TBranchDescriptor(cl->GetName(), objInfo, branchname.Data(), local_prefix.Data(),
                                                      isclones, containerName, desc);
                        // Recurse: analyze sub-branches of the sub-branch
                        lookedAt += AnalyzeBranches(bdesc, branch, objInfo);
                        isLeaf = false;

                     }
                  } else { // The element and the branch does not match, we need to loop over the next branches
                     Int_t pos = branchname.Last('.');
                     if (pos != -1) {
                        branchname.Remove(pos);
                     }
                     TString local_prefix = desc ? desc->fSubBranchPrefix : TString(parent->GetName());
                     objInfo = GetBaseClass(element);
                     if (objInfo == 0) {
                        continue; // There is no data in this base class
                     }
                     cl = objInfo->GetClass();
                     bdesc = new TBranchDescriptor(cl->GetName(), objInfo, branchname.Data(), local_prefix.Data(),
                                                    isclones, containerName, desc);
                     usedBranch = kFALSE;
                     // Recurse: analyze the sub-elements with the same branches
                     lookedAt += AnalyzeBranches(bdesc, branches, objInfo);
                  }
                  delete bdesc;
               } else { // Not base class
                  TBranchDescriptor *bdesc = 0;
                  if (branchEndName == element->GetName()) { // The element and the branch matches
                     if (branch->GetListOfBranches()->GetEntries() == 0) { // The branch contains a non-split class
                        // FIXME: nothing to do in such cases, because readers cannot access
                        // non-split members and a reader for the whole branch will be added
                     } else { // The branch contains a split class
                        if (isclones != kOut) {
                           // We have to guess the version number!
                           cl = TClass::GetClass(dataType);
                           objInfo = GetStreamerInfo(branch, branch->GetListOfBranches(), cl);
                        }
                        bdesc = new TBranchDescriptor(cl->GetName(), objInfo, branch->GetName(), branch->GetName(),
                                                      isclones, containerName, desc);
                        // Recurse: analyze sub-branches of the sub-branch
                        lookedAt += AnalyzeBranches(bdesc, branch, objInfo);
                        isLeaf = false;
                     }
                  } else { // The element and the branch does not match, we need to loop over the next branches
                     TString local_prefix = desc ? desc->fSubBranchPrefix : TString(parent->GetName());
                     if (local_prefix.Length()) local_prefix += ".";
                     local_prefix += element->GetName();
                     objInfo = branch->GetInfo();
                     Int_t pos = branchname.Last('.');
                     if (pos != -1) {
                        branchname.Remove(pos);
                     }
                     if (isclones != kOut) {
                        // We have to guess the version number!
                        cl = TClass::GetClass(dataType);
                        objInfo = GetStreamerInfo(branch, branches, cl);
                     }
                     bdesc = new TBranchDescriptor(cl->GetName(), objInfo, branchname.Data(), local_prefix.Data(),
                                                   isclones, containerName, desc);
                     usedBranch = kFALSE;
                     skipped = kTRUE;
                     // Recurse: analyze the sub-elements with the same branches
                     lookedAt += AnalyzeBranches(bdesc, branches, objInfo);
                  }
                  delete bdesc;
               }

               break;
            }
            default:
               Error("AnalyzeBranch", "Unsupported type for %s (%d).", branch->GetName(), element->GetType());
         }

         if (!isBase && !skipped) { // Add reader for the whole branch
            if (outer_isclones != kOut && readerType == TTreeReaderDescriptor::ReaderType::kArray) {
               Error("AnalyzeBranch", "Arrays inside collections are not supported yet (branch: %s).", branch->GetName());
            } else {
               if (outer_isclones != kOut || isclones != kOut) {
                  readerType = TTreeReaderDescriptor::ReaderType::kArray;
               }
               AddReader(readerType, dataType, branch->GetName(), branch->GetName(), desc, isLeaf);
            }
         }

         // If the branch was used, jump to the next
         if (usedBranch) {
            branches.Next();
            ++lookedAt;
         }
      }

      return lookedAt;
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Analyze branch and add the variables found. The number of analyzed
   /// sub-branches is returned.

   UInt_t TTreeReaderGenerator::AnalyzeOldBranch(TBranch *branch)
   {
      UInt_t extraLookedAt = 0;
      TString prefix;

      TString branchName = branch->GetName();

      TObjArray *leaves = branch->GetListOfLeaves();
      Int_t nleaves = leaves ? leaves->GetEntriesFast() : 0;

      // Loop through leaves and analyze them
      for(int l=0;l<nleaves;l++) {
         TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(l);
         extraLookedAt += AnalyzeOldLeaf(leaf, nleaves);
      }

      return extraLookedAt;
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Analyze the leaf and add the variables found.

   UInt_t TTreeReaderGenerator::AnalyzeOldLeaf(TLeaf *leaf, Int_t nleaves)
   {
      if (leaf->IsA()==TLeafObject::Class()) {
         Error("AnalyzeOldLeaf","TLeafObject not supported yet");
         return 0;
      }

      TString leafTypeName = leaf->GetTypeName();
      Int_t pos = leafTypeName.Last('_');
      //if (pos != -1) leafTypeName.Remove(pos); // FIXME: this is not required since it makes Float_t -> Float

      // Analyze dimensions
      UInt_t dim = 0;
      std::vector<Int_t> maxDim;

      TString dimensions;
      TString temp = leaf->GetName();
      pos = temp.Index("[");
      if (pos != -1) {
         if (pos) temp.Remove(0, pos);
         dimensions.Append(temp);
      }
      temp = leaf->GetTitle();
      pos = temp.Index("[");
      if (pos != -1) {
         if (pos) temp.Remove(0, pos);
         dimensions.Append(temp);
      }

      Int_t dimlen = dimensions.Length();

      if (dimlen) {
         const char *current = dimensions.Data();

         Int_t index;
         Int_t scanindex ;
         while (current) {
            current++;
            if (current[0] == ']') {
               maxDim.push_back(-1); // maxDim[dim] = -1; // Loop over all elements;
            } else {
               scanindex = sscanf(current,"%d",&index);
               if (scanindex) {
                  maxDim.push_back(index); // maxDim[dim] = index;
               } else {
                  maxDim.push_back(-2); // maxDim[dim] = -2; // Index is calculated via a variable.
               }
            }
            dim ++;
            current = (char*)strstr( current, "[" );
         }
      }

      if (dim == 0 && leaf->IsA() == TLeafC::Class()) {
         dim = 1; // For C style strings
      }

      TTreeReaderDescriptor::ReaderType type = TTreeReaderDescriptor::ReaderType::kValue;
      TString dataType;
      switch (dim) {
         case 0: {
            type = TTreeReaderDescriptor::ReaderType::kValue;
            dataType = leafTypeName;
            break;
         }
         case 1: {
            type = TTreeReaderDescriptor::ReaderType::kArray;
            dataType = leafTypeName;
            break;
         }
         default: {
            // TODO: transform this
            /*type = "TArrayProxy<";
            for(Int_t ind = dim - 2; ind > 0; --ind) {
               type += "TMultiArrayType<";
            }
            type += "TArrayType<";
            type += leaf->GetTypeName();
            type += ",";
            type += maxDim[dim-1];
            type += "> ";
            for(Int_t ind = dim - 2; ind > 0; --ind) {
               type += ",";
               type += maxDim[ind];
               type += "> ";
            }
            type += ">";*/
            break;
         }
      }

      // If there are multiple leaves (leaflist) the name of the branch is
      // <branch_name>.<leaf_name>
      // (Otherwise the brach name does not change)
      TString branchName = leaf->GetBranch()->GetName();
      if (nleaves > 1) {
         branchName.Form("%s.%s", leaf->GetBranch()->GetName(), leaf->GetName());
      }

      AddReader(type, dataType, leaf->GetName(), branchName, 0, kTRUE);

      return 0;
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Check whether a branch should have a corresponding reader added, depending
   /// on the options provided by the user.

   Bool_t TTreeReaderGenerator::BranchNeedsReader(TString branchName, TBranchDescriptor *parent, Bool_t isLeaf)
   {
      if (isLeaf) { // Branch is a leaf
         // Include if all leaves should be included or it is contained in any of the lists.
         if (fIncludeAllLeaves) return kTRUE;
         if (std::find(fIncludeLeaves.begin(), fIncludeLeaves.end(), branchName) != fIncludeLeaves.end()) return kTRUE;
         if (std::find(fIncludeStruct.begin(), fIncludeStruct.end(), branchName) != fIncludeStruct.end()) return kTRUE;
         if (!parent) { // Branch is topmost (top-level leaf)
            if (fIncludeAllTopmost) return kTRUE;
         } else {       // Branch is not topmost
            while (parent) { // Check if any parent is in the list of "include as leaves"
               if (std::find(fIncludeLeaves.begin(), fIncludeLeaves.end(), parent->fBranchName) != fIncludeLeaves.end()) {
                  return kTRUE;
               }
               parent = parent->fParent;
            }
         }
      } else {      // Branch is not a leaf (has sub-branches)
         if (std::find(fIncludeStruct.begin(), fIncludeStruct.end(), branchName) != fIncludeStruct.end()) return kTRUE;
         if (!parent) { // Branch is topmost
            if (fIncludeAllTopmost) return kTRUE;
         }
      }
      return false;
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Parse the user options.

   void TTreeReaderGenerator::ParseOptions() {
      if (fOptionStr.EqualTo("")) { // Empty -> include all leaves
         fIncludeAllLeaves = kTRUE;
      } else if (fOptionStr.EqualTo("@")) { // "@" -> include all topmost
         fIncludeAllTopmost = kTRUE;
      } else { // Otherwise split at ";" to get names
         TObjArray *tokens = fOptionStr.Tokenize(TString(";"));
         for (Int_t i = 0; i < tokens->GetEntries(); ++i) {
            TString token = ((TObjString*)tokens->At(i))->GetString();
            if ( token.Length() == 0 || (token.Length() == 1 && token[0] == '@') ) {
               Warning("ParseOptions", "Ignored empty branch name in option string.");
            } else if (token[0] == '@') { // "@X" means include X as a whole
               token = TString(token.Data()+1);
               fIncludeStruct.push_back(token);
            } else {                      // "X"  means include leaves of X
               fIncludeLeaves.push_back(token);
            }
            if (!fTree->GetBranch(token)) { // Display a warning for non-existing branch names
               Warning("ParseOptions", "Tree %s does not contain a branch named %s.", fTree->GetName(), token.Data());
            }
         }
         delete tokens;
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Analyze tree and extract readers.

   void TTreeReaderGenerator::AnalyzeTree(TTree *tree)
   {
      TIter next(tree->GetListOfBranches());
      TBranch *branch;

      // Loop through branches
      while ( (branch = (TBranch*)next()) ) {
         TVirtualStreamerInfo *info = 0;
         // Get the name and the class of the branch
         const char *branchName = branch->GetName();
         const char *branchClassName = branch->GetClassName();
         TClass *cl = TClass::GetClass(branchClassName);

         // Add headers for user classes
         if (branchClassName && strlen(branchClassName)) {
            AddHeader(cl);
         }

         TString type = "unknown";   // Type of branch
         ELocation isclones = kOut;  // Type of container
         TString containerName = ""; // Name of container
         TBranchDescriptor *desc = 0;
         // Check whether the branch is a container
         if (cl) {
            // Check if it is a TClonesArray
            if (cl == TClonesArray::Class()) {
               isclones = kClones;
               containerName = "TClonesArray";
               if (branch->IsA()==TBranchElement::Class()) {
                  // Get the class inside the TClonesArray
                  const char *cname = ((TBranchElement*)branch)->GetClonesName();
                  TClass *ncl = TClass::GetClass(cname);
                  if (ncl) {
                     cl = ncl;
                     info = GetStreamerInfo(branch, branch->GetListOfBranches(), cl);
                  } else {
                     Error("AnalyzeTree",
                           "Introspection of TClonesArray in older file not implemented yet.");
                  }
               } else {
                  TClonesArray **ptr = (TClonesArray**)branch->GetAddress();
                  TClonesArray *clones = 0;
                  if (ptr==0) {
                     clones = new TClonesArray;
                     branch->SetAddress(&clones);
                     ptr = &clones;
                  }
                  branch->GetEntry(0);
                  TClass *ncl = *ptr ? (*ptr)->GetClass() : 0;
                  if (ncl) {
                     cl = ncl;
                  } else {
                     Error("AnalyzeTree",
                           "Introspection of TClonesArray for %s failed.",branch->GetName());
                  }
               }
            // Check if it is an STL container
            } else if (cl->GetCollectionProxy()) {
               isclones = kSTL;
               containerName = cl->GetName();
               // Get the type inside container
               if (cl->GetCollectionProxy()->GetValueClass()) { // Class inside container
                  cl = cl->GetCollectionProxy()->GetValueClass();
               } else { // RAW type (or missing class) inside container
                  // TODO: CheckForMissingClass?
                  // TTreeReaderArray does currently not support vectors of bool so that need to
                  // be added as a TTreeReaderValue<vector<bool>>. Also getting the inner type of
                  // vector would return "unsigned" so the full name has to be compared.
                  if (containerName.EqualTo("vector<bool>")) {
                     AddReader(TTreeReaderDescriptor::ReaderType::kValue,
                            containerName,
                            branch->GetName(), branch->GetName(), 0, kTRUE);
                  } else { // Otherwise we can generate a TTreeReaderArray with the inner type
                     AddReader(TTreeReaderDescriptor::ReaderType::kArray,
                            TDataType::GetDataType(cl->GetCollectionProxy()->GetType())->GetName(),
                            branch->GetName(), branch->GetName(), 0, kTRUE);
                  }
                  continue; // Nothing else to with this branch in these cases
               }
            }

            // Check class inside container and create a descriptor or add a reader
            if (cl) {
               if (cl->TestBit(TClass::kIsEmulation) || branchName[strlen(branchName)-1] == '.' || branch->GetSplitLevel()) {
                  TBranchElement *be = dynamic_cast<TBranchElement*>(branch);
                  TVirtualStreamerInfo *beinfo = (be && isclones == kOut)
                     ? be->GetInfo() : cl->GetStreamerInfo(); // the 2nd hand need to be fixed
                  // Create descriptor
                  desc = new TBranchDescriptor(cl->GetName(), beinfo, branchName, branchName, isclones, containerName);
                  info = beinfo;
               } else {
                  // Add a reader for non-split classes
                  AddReader(isclones == kOut ?
                              TTreeReaderDescriptor::ReaderType::kValue
                            : TTreeReaderDescriptor::ReaderType::kArray,
                            cl->GetName(), branchName, branchName, 0, kTRUE);
                  // TODO: can't we just put a continue here?
               }
            }
         }

         // Analyze sub-branches (if exist) and add readers
         if (branch->GetListOfBranches()->GetEntries() == 0) { // Branch is non-splitted
            if (cl) { // Non-split object
               if (desc) { // If there is a descriptor add reader (otherwise
                           // it was already added).
                  AddReader(isclones == kOut ?
                              TTreeReaderDescriptor::ReaderType::kValue
                            : TTreeReaderDescriptor::ReaderType::kArray,
                            desc->GetName(), desc->fBranchName, desc->fBranchName, 0, kTRUE);
               }
            } else { // Top-level RAW type
               AnalyzeOldBranch(branch); // Analyze branch and extract readers
            }
         } else { // Branch is splitted
            TIter subnext( branch->GetListOfBranches() );
            if (desc) {
               // Analyze sub-branches and extract readers
               TBranchElement *branchElem = dynamic_cast<TBranchElement*>(branch);
               if (branchElem) {
                  AnalyzeBranches(desc, branchElem, info);
               } else {
                  Error("AnalyzeTree", "Cannot analyze branch %s because it is not a TBranchElement.", branchName);
               }
               // Also add a reader for the whole branch
               AddReader(isclones == kOut ?
                              TTreeReaderDescriptor::ReaderType::kValue
                            : TTreeReaderDescriptor::ReaderType::kArray,
                            desc->GetName(), desc->fBranchName, desc->fBranchName, 0, kFALSE);
            }
         }
         delete desc;
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Generate code for selector class.

   void TTreeReaderGenerator::WriteSelector()
   {
      // If no name is given, set to default (name of the tree)
      if (!fClassname) fClassname = fTree->GetName();

      TString treefile;
      if (fTree->GetDirectory() && fTree->GetDirectory()->GetFile()) {
         treefile = fTree->GetDirectory()->GetFile()->GetName();
      } else {
         treefile = "Memory Directory";
      }
      // In the case of a chain, the GetDirectory information usually does
      // pertain to the Chain itself but to the currently loaded tree.
      // So we can not rely on it.
      Bool_t ischain = fTree->InheritsFrom(TChain::Class());
      Bool_t isHbook = fTree->InheritsFrom("THbookTree");
      if (isHbook)
         treefile = fTree->GetTitle();

      //======================Generate classname.h=====================
      TString thead;
      thead.Form("%s.h", fClassname.Data());
      std::ofstream ofs (thead, std::ofstream::out);
      if (!ofs) {
         Error("WriteSelector","cannot open output file %s", thead.Data());
         return;
      }
      // Print header
      TDatime td;
      ofs <<
R"CODE(//////////////////////////////////////////////////////////
// This class has been automatically generated on
// )CODE" << td.AsString() << R"CODE( by ROOT version )CODE" << gROOT->GetVersion() << std::endl;
      if (!ischain) {
         ofs << "// from TTree " << fTree->GetName() << "/" << fTree->GetTitle() << std::endl
             << "// found on file: " << treefile << std::endl;
      } else {
         ofs << "// from TChain " << fTree->GetName() << "/" << fTree->GetTitle() << std::endl;
      }
      ofs <<
R"CODE(//////////////////////////////////////////////////////////

#ifndef )CODE" << fClassname << R"CODE(_h
#define )CODE" << fClassname << R"CODE(_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
)CODE";
      if (isHbook) ofs << "#include <THbookFile.h>" << std::endl;
      ofs <<
R"CODE(#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
)CODE";

      TIter next(&fListOfHeaders);
      TObject *header;
      while ( (header = next()) ) {
         ofs << header->GetTitle() << std::endl;
      }
      ofs << std::endl << std::endl;

      // Generate class declaration with TTreeReaderValues and Arrays
      ofs <<
R"CODE(class )CODE" << fClassname << R"CODE( : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
)CODE";
      next = &fListOfReaders;
      TTreeReaderDescriptor *descriptor;
      while ( ( descriptor = (TTreeReaderDescriptor*)next() ) ) {
         ofs << "   TTreeReader" << (descriptor->fType == TTreeReaderDescriptor::ReaderType::kValue ? "Value" : "Array")
                                 << "<" << descriptor->fDataType
                                 << "> " << descriptor->fName
                                 << " = {fReader, \"" << descriptor->fBranchName << "\"};" << std::endl;
      }
      // Generate class member functions prototypes
      ofs <<
R"CODE(

   )CODE" << fClassname << R"CODE((TTree * /*tree*/ =0) { }
   virtual ~)CODE" << fClassname << R"CODE(() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef()CODE" << fClassname << R"CODE(,0);

};

#endif

#ifdef )CODE" << fClassname << R"CODE(_cxx
void )CODE" << fClassname << R"CODE(::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t )CODE" << fClassname << R"CODE(::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef )CODE" << fClassname << R"CODE(_cxx
)CODE";
      ofs.close();

      //======================Generate classname.C=====================
      TString tcimp;
      tcimp.Form("%s.C", fClassname.Data());
      std::ofstream ofsc (tcimp, std::ofstream::out);
      if (!ofsc) {
         Error("WriteSelector","cannot open output file %s", tcimp.Data());
         return;
      }

      ofsc <<
R"CODE(#define )CODE" << fClassname << R"CODE(_cxx
// The class definition in )CODE" << fClassname << R"CODE(.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process(")CODE" << fClassname << R"CODE(.C")
// root> T->Process(")CODE" << fClassname << R"CODE(.C","some options")
// root> T->Process(")CODE" << fClassname << R"CODE(.C+")
//


#include ")CODE" << thead << R"CODE("
#include <TH2.h>
#include <TStyle.h>

void )CODE" << fClassname << R"CODE(::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void )CODE" << fClassname << R"CODE(::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t )CODE" << fClassname << R"CODE(::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);

   return kTRUE;
}

void )CODE" << fClassname << R"CODE(::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void )CODE" << fClassname << R"CODE(::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

})CODE";
      ofsc.close();
   }

} // namespace Internal
} // namespace ROOT
