// @(#)root/treeplayer:$Id$
// Author: Akos Hajdu 13/08/2015

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers and al.        *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TTreeGeneratorBase.h"

#include "TBranchElement.h"
#include "TClass.h"
#include "TClassEdit.h"
#include "TClonesArray.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TStreamerElement.h"
#include "TStreamerInfo.h"
#include "TTree.h"
#include "TError.h"
#include "TVirtualCollectionProxy.h"
#include "TVirtualStreamerInfo.h"

namespace ROOT {
namespace Internal {

   ////////////////////////////////////////////////////////////////////////////////
   /// Constructor.

   TTreeGeneratorBase::TTreeGeneratorBase(TTree *tree, const char *option) : fTree(tree), fOptionStr(option) { }

   ////////////////////////////////////////////////////////////////////////////////
   /// Add a header inclusion request. If the header is already included it will
   /// not be included again.

   void TTreeGeneratorBase::AddHeader(TClass *cl)
   {
      if (cl==0) return;

      // Check if already included
      TObject *obj = fListOfHeaders.FindObject(cl->GetName());
      if (obj) return;

      TString directive;

      // Extract inner class from collection
      if (cl->GetCollectionProxy() && cl->GetCollectionProxy()->GetValueClass()) {
         AddHeader( cl->GetCollectionProxy()->GetValueClass() );
      }

      // Construct directive
      Int_t stlType;
      if (0 == strcmp(cl->GetName(), "string")) { // Check if it is a string
         directive = "#include <string>\n";
      } else if (cl->GetCollectionProxy() && (stlType = cl->GetCollectionType())) { // Check if it is an STL container
         const char *what = "";
         switch(stlType)  {
            case  ROOT::kSTLvector:            what = "vector"; break;
            case  ROOT::kSTLlist:              what = "list"; break;
            case  ROOT::kSTLforwardlist:       what = "forward_list"; break;
            case -ROOT::kSTLdeque:             // same as positive
            case  ROOT::kSTLdeque:             what = "deque"; break;
            case -ROOT::kSTLmap:               // same as positive
            case  ROOT::kSTLmap:               what = "map"; break;
            case -ROOT::kSTLmultimap:          // same as positive
            case  ROOT::kSTLmultimap:          what = "map"; break;
            case -ROOT::kSTLset:               // same as positive
            case  ROOT::kSTLset:               what = "set"; break;
            case -ROOT::kSTLmultiset:          // same as positive
            case  ROOT::kSTLmultiset:          what = "set"; break;
            case -ROOT::kSTLunorderedset:      // same as positive
            case  ROOT::kSTLunorderedset:      what = "unordered_set"; break;
            case -ROOT::kSTLunorderedmultiset: // same as positive
            case  ROOT::kSTLunorderedmultiset: what = "unordered_multiset"; break;
            case -ROOT::kSTLunorderedmap:      // same as positive
            case  ROOT::kSTLunorderedmap:      what = "unordered_map"; break;
            case -ROOT::kSTLunorderedmultimap: // same as positive
            case  ROOT::kSTLunorderedmultimap: what = "unordered_multimap"; break;
         }
         if (what[0]) {
            directive = "#include <";
            directive.Append(what);
            directive.Append(">\n");
         }
      } else if (cl->GetDeclFileName() && strlen(cl->GetDeclFileName()) ) { // Custom file
         const char *filename = cl->GetDeclFileName();

         if (!filename) return;

#ifdef R__WIN32
         TString inclPath("include;prec_stl"); // GetHtml()->GetIncludePath());
#else
         TString inclPath("include:prec_stl"); // GetHtml()->GetIncludePath());
#endif
         Ssiz_t posDelim = 0;
         TString inclDir;
         TString sIncl(filename);
#ifdef R__WIN32
         const char* pdelim = ";";
         static const char ddelim = '\\';
#else
         const char* pdelim = ":";
         static const char ddelim = '/';
#endif
         while (inclPath.Tokenize(inclDir, posDelim, pdelim))
         {
            if (sIncl.BeginsWith(inclDir)) {
               filename += inclDir.Length();
               if (filename[0] == ddelim || filename[0] == '/') {
                  ++filename;
               }
               break;
            }
         }
         directive = Form("#include \"%s\"\n",filename);
      } else if (TClassEdit::IsStdPair(cl->GetName())) {
         TClassEdit::TSplitType split(cl->GetName());
         // 4 elements expected: "pair", "first type name", "second type name", "trailing stars"
         // However legacy code had a test for 3, we will leave it here until
         // a test is developed (or found :) ) that exercise these lines of code.
         if (split.fElements.size() == 3 || split.fElements.size() == 4) {
            for (int arg = 1; arg < 3; ++arg) {
               TClass* clArg = TClass::GetClass(split.fElements[arg].c_str());
               if (clArg) AddHeader(clArg);
            }
         }
      }
      // Add directive (if it is not added already)
      if (directive.Length()) {
         TIter i( &fListOfHeaders );
         for(TNamed *n = (TNamed *)i(); n; n = (TNamed*)i()) {
            if (directive == n->GetTitle()) {
               return;
            }
         }
         fListOfHeaders.Add(new TNamed(cl->GetName(), directive.Data()));
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Add a header inclusion request. If the header is already included it will
   /// not be included again.

   void TTreeGeneratorBase::AddHeader(const char *classname)
   {
      AddHeader(TClass::GetClass(classname));
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Get name of class inside a container.

   TString TTreeGeneratorBase::GetContainedClassName(TBranchElement *branch, TStreamerElement *element, Bool_t ispointer)
   {
      TString cname = branch->GetClonesName();
      if (cname.Length()==0) {
         // We may have any unsplit clones array
         Long64_t i = branch->GetTree()->GetReadEntry();
         if (i<0) i = 0;
         branch->GetEntry(i);
         char *obj = branch->GetObject();

         TBranchElement *parent = (TBranchElement*)branch->GetMother()->GetSubBranch(branch);
         const char *pclname = parent->GetClassName();

         TClass *clparent = TClass::GetClass(pclname);
         // TClass *clm = TClass::GetClass(GetClassName());
         Int_t lOffset = 0; // offset in the local streamerInfo.
         if (clparent) {
            const char *ename = 0;
            if (element) {
               ename = element->GetName();
               lOffset = clparent->GetStreamerInfo()->GetOffset(ename);
            } else {
               lOffset = 0;
            }
         }
         else Error("AnalyzeBranch", "Missing parent for %s.", branch->GetName());

         TClonesArray *arr;
         if (ispointer) {
            arr = (TClonesArray*)*(void**)(obj+lOffset);
         } else {
            arr = (TClonesArray*)(obj+lOffset);
         }
         cname = arr->GetClass()->GetName();
      }
      if (cname.Length()==0) {
         Error("AnalyzeBranch",
            "Introspection of TClonesArray in older file not implemented yet.");
      }
      return cname;
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Check if element is a base class and if yes, return the base class.

   TVirtualStreamerInfo *TTreeGeneratorBase::GetBaseClass(TStreamerElement *element)
   {
      TStreamerBase *base = dynamic_cast<TStreamerBase*>(element);
      if (base) {
         TVirtualStreamerInfo *info = base->GetBaseStreamerInfo();
         if (info) return info;
      }
      return 0;
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Return the correct TStreamerInfo of class 'cl' in the list of branches
   /// (current) [Assuming these branches correspond to a flattened version of
   /// the class.]

   TVirtualStreamerInfo *TTreeGeneratorBase::GetStreamerInfo(TBranch *branch, TIter current, TClass *cl)
   {
      TVirtualStreamerInfo *objInfo = 0;
      TBranchElement *b = 0;
      TString cname = cl->GetName();

      while( ( b = (TBranchElement*)current() ) ) {
         if ( cname == b->GetInfo()->GetName() ) {
            objInfo = b->GetInfo();
            break;
         }
      }
      if (objInfo == 0 && branch->GetTree()->GetDirectory()->GetFile()) {
         const TList *infolist = branch->GetTree()->GetDirectory()->GetFile()->GetStreamerInfoCache();
         if (infolist) {
            TVirtualStreamerInfo *i = (TVirtualStreamerInfo *)infolist->FindObject(cname);
            if (i) {
               // NOTE: Is this correct for Foreigh classes?
               objInfo = (TVirtualStreamerInfo *)cl->GetStreamerInfo(i->GetClassVersion());
            }
         }
      }
      if (objInfo == 0) {
         // We still haven't found it ... this is likely to be an STL collection .. anyway, use the current StreamerInfo.
         objInfo = cl->GetStreamerInfo();
      }
      return objInfo;
   }

} // namespace Internal
} // namespace ROOT
