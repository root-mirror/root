// @(#)root/treeplayer:$Id$
// Author: Rene Brun   19/01/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TTreeFormula.h"
#include "TList.h"
#include "TTree.h"
#include "TBuffer.h"
#include "TBranch.h"
#include "TBranchObject.h"
#include "TBranchElement.h"
#include "TClonesArray.h"
#include "TLeafB.h"
#include "TLeafC.h"
#include "TLeafElement.h"
#include "TLeafObject.h"
#include "TMethodCall.h"
#include "TCutG.h"
#include "TRandom.h"
#include "TInterpreter.h"
#include "TDataType.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"
#include "TArrayI.h"
#include "TAxis.h"
#include "TError.h"
#include "TVirtualCollectionProxy.h"
#include "TString.h"
#include "TMath.h"

#include "TVirtualRefProxy.h"
#include "TTreeFormulaManager.h"
#include "TFormLeafInfo.h"
#include "TMethod.h"
#include "TFormLeafInfoReference.h"
#include "strlcpy.h"
#include "snprintf.h"
#include "TEntryList.h"

#include <cctype>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <typeinfo>
#include <algorithm>

const Int_t kMaxLen     = 1024;

/** \class TTreeFormula
Used to pass a selection expression to the Tree drawing routine. See TTree::Draw

A TreeFormula can contain any arithmetic expression including
standard operators and mathematical functions separated by operators.
Examples of valid expression:
~~~{.cpp}
     "x<y && sqrt(z)>3.2"
~~~
TTreeFormula now relies on a variety of TFormLeafInfo classes to handle the
reading of the information. Here is the list of theses classes:
  - TFormLeafInfo
  - TFormLeafInfoDirect
  - TFormLeafInfoNumerical
  - TFormLeafInfoClones
  - TFormLeafInfoCollection
  - TFormLeafInfoPointer
  - TFormLeafInfoMethod
  - TFormLeafInfoMultiVarDim
  - TFormLeafInfoMultiVarDimDirect
  - TFormLeafInfoCast

The following method are available from the TFormLeafInfo interface:

 -  AddOffset(Int_t offset, TStreamerElement* element)
 -  GetCounterValue(TLeaf* leaf) : return the size of the array pointed to.
 -  GetObjectAddress(TLeafElement* leaf) : Returns the the location of the object pointed to.
 -  GetMultiplicity() : Returns info on the variability of the number of elements
 -  GetNdata(TLeaf* leaf) : Returns the number of elements
 -  GetNdata() : Used by GetNdata(TLeaf* leaf)
 -  GetValue(TLeaf *leaf, Int_t instance = 0) : Return the value
 -  GetValuePointer(TLeaf *leaf, Int_t instance = 0) : Returns the address of the value
 -  GetLocalValuePointer(TLeaf *leaf, Int_t instance = 0) : Returns the address of the value of 'this' LeafInfo
 -  IsString()
 -  ReadValue(char *where, Int_t instance = 0) : Internal function to interpret the location 'where'
 -  Update() : react to the possible loading of a shared library.
*/

ClassImp(TTreeFormula);

////////////////////////////////////////////////////////////////////////////////

inline static void R__LoadBranch(TBranch* br, Long64_t entry, Bool_t quickLoad)
{
   if (!quickLoad || (br->GetReadEntry() != entry)) {
      br->GetEntry(entry);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// \class TDimensionInfo
/// A small helper class to help in keeping track of the array
/// dimensions encountered in the analysis of the expression.

class TDimensionInfo : public TObject {
public:
   Int_t fCode;  // Location of the leaf in TTreeFormula::fCode
   Int_t fOper;  // Location of the Helper using the leaf in TTreeFormula::fOper
   Int_t fSize;
   TFormLeafInfoMultiVarDim* fMultiDim;
   TDimensionInfo(Int_t code, Int_t oper, Int_t size, TFormLeafInfoMultiVarDim* multiDim)
      : fCode(code), fOper(oper), fSize(size), fMultiDim(multiDim) {};
   ~TDimensionInfo() {};
};

////////////////////////////////////////////////////////////////////////////////

TTreeFormula::TTreeFormula(): ROOT::v5::TFormula(), fQuickLoad(kFALSE), fNeedLoading(kTRUE),
   fDidBooleanOptimization(kFALSE), fDimensionSetup(0)

{
   // Tree Formula default constructor

   fTree         = 0;
   fLookupType   = 0;
   fNindex       = 0;
   fNcodes       = 0;
   fAxis         = 0;
   fHasCast      = 0;
   fManager      = 0;
   fMultiplicity = 0;
   fConstLD      = 0;

   Int_t j,k;
   for (j=0; j<kMAXCODES; j++) {
      fNdimensions[j] = 0;
      fCodes[j] = 0;
      fNdata[j] = 1;
      fHasMultipleVarDim[j] = kFALSE;
      for (k = 0; k<kMAXFORMDIM; k++) {
         fIndexes[j][k] = -1;
         fCumulSizes[j][k] = 1;
         fVarIndexes[j][k] = 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Normal TTree Formula Constructor

TTreeFormula::TTreeFormula(const char *name,const char *expression, TTree *tree)
   :ROOT::v5::TFormula(), fTree(tree), fQuickLoad(kFALSE), fNeedLoading(kTRUE),
    fDidBooleanOptimization(kFALSE), fDimensionSetup(0)
{
   Init(name,expression);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor used during the expansion of an alias

TTreeFormula::TTreeFormula(const char *name,const char *expression, TTree *tree,
                           const std::vector<std::string>& aliases)
   :ROOT::v5::TFormula(), fTree(tree), fQuickLoad(kFALSE), fNeedLoading(kTRUE),
    fDidBooleanOptimization(kFALSE), fDimensionSetup(0), fAliasesUsed(aliases)
{
   Init(name,expression);
}

////////////////////////////////////////////////////////////////////////////////
/// Initialization called from the constructors.

void TTreeFormula::Init(const char*name, const char* expression)
{
   TDirectory *const savedir=gDirectory;

   fNindex       = kMAXFOUND;
   fLookupType   = new Int_t[fNindex];
   fNcodes       = 0;
   fMultiplicity = 0;
   fAxis         = 0;
   fHasCast      = 0;
   fConstLD      = 0;
   Int_t i,j,k;
   fManager      = new TTreeFormulaManager;
   fManager->Add(this);

   for (j=0; j<kMAXCODES; j++) {
      fNdimensions[j] = 0;
      fLookupType[j] = kDirect;
      fCodes[j] = 0;
      fNdata[j] = 1;
      fHasMultipleVarDim[j] = kFALSE;
      for (k = 0; k<kMAXFORMDIM; k++) {
         fIndexes[j][k] = -1;
         fCumulSizes[j][k] = 1;
         fVarIndexes[j][k] = 0;
      }
   }

   fDimensionSetup = new TList;

   if (Compile(expression)) {
      fTree = 0; fNdim = 0;
      if(savedir) savedir->cd();
      return;
   }

   if (fNcodes >= kMAXFOUND) {
      Warning("TTreeFormula","Too many items in expression:%s",expression);
      fNcodes = kMAXFOUND;
   }
   SetName(name);

   for (i=0;i<fNoper;i++) {

      if (GetAction(i)==kDefinedString) {
         Int_t string_code = GetActionParam(i);
         TLeaf *leafc = (TLeaf*)fLeaves.UncheckedAt(string_code);
         if (!leafc) continue;

         // We have a string used as a string

         // This dormant portion of code would be used if (when?) we allow the histogramming
         // of the integral content (as opposed to the string content) of strings
         // held in a variable size container delimited by a null (as opposed to
         // a fixed size container or variable size container whose size is controlled
         // by a variable).  In GetNdata, we will then use strlen to grab the current length.
         //fCumulSizes[i][fNdimensions[i]-1] = 1;
         //fUsedSizes[fNdimensions[i]-1] = -TMath::Abs(fUsedSizes[fNdimensions[i]-1]);
         //fUsedSizes[0] = - TMath::Abs( fUsedSizes[0]);

         if (fNoper == 1) {
            // If the string is by itself, then it can safely be histogrammed as
            // in a string based axis.  To histogram the number inside the string
            // just make it part of a useless expression (for example: mystring+0)
            SetBit(kIsCharacter);
         }
         continue;
      }
      if (GetAction(i)==kJump && GetActionParam(i)==(fNoper-1)) {
         // We have cond ? string1 : string2
         if (IsString(fNoper-1)) SetBit(kIsCharacter);
      }
   }
   if (fNoper == 1 && GetAction(0)==kStringConst) {
      SetBit(kIsCharacter);
   }
   if (fNoper==1 && GetAction(0)==kAliasString) {
      TTreeFormula *subform = static_cast<TTreeFormula*>(fAliases.UncheckedAt(0));
      R__ASSERT(subform);
      if (subform->IsString()) SetBit(kIsCharacter);
   } else if (fNoper==2 && GetAction(0)==kAlternateString) {
      TTreeFormula *subform = static_cast<TTreeFormula*>(fAliases.UncheckedAt(0));
      R__ASSERT(subform);
      if (subform->IsString()) SetBit(kIsCharacter);
   }

   fManager->Sync();

   // Let's verify the indexes and dies if we need to.
   Int_t k0,k1;
   for(k0 = 0; k0 < fNcodes; k0++) {
      for(k1 = 0; k1 < fNdimensions[k0]; k1++ ) {
         // fprintf(stderr,"Saw %d dim %d and index %d\n",k1, fFixedSizes[k0][k1], fIndexes[k0][k1]);
         if ( fIndexes[k0][k1]>=0 && fFixedSizes[k0][k1]>=0
              && fIndexes[k0][k1]>=fFixedSizes[k0][k1]) {
            Error("TTreeFormula",
                  "Index %d for dimension #%d in %s is too high (max is %d)",
                  fIndexes[k0][k1],k1+1, expression,fFixedSizes[k0][k1]-1);
            fTree = 0; fNdim = 0;
            if(savedir) savedir->cd();
            return;
         }
      }
   }

   // Create a list of uniques branches to load.
   for(k=0; k<fNcodes ; k++) {
      TLeaf *leaf = k <= fLeaves.GetLast() ? (TLeaf*)fLeaves.UncheckedAt(k) : 0;
      TBranch *branch = 0;
      if (leaf) {
         branch = leaf->GetBranch();
         if (fBranches.FindObject(branch)) branch = 0;
      }
      fBranches.AddAtAndExpand(branch,k);
   }

   if (IsInteger(kFALSE)) SetBit(kIsInteger);

   if (TestBit(TTreeFormula::kNeedEntries)) {
      // Call TTree::GetEntries() to insure that it is already calculated.
      // This will need to be done anyway at the first iteration and insure
      // that it will not mess up the branch reading (because TTree::GetEntries
      // opens all the file in the chain and 'stays' on the last file.

      Long64_t readentry = fTree->GetReadEntry();
      Int_t treenumber = fTree->GetTreeNumber();
      fTree->GetEntries();
      if (treenumber != fTree->GetTreeNumber()) {
         if (readentry >= 0) {
            fTree->LoadTree(readentry);
         }
         UpdateFormulaLeaves();
      } else {
         if (readentry >= 0) {
            fTree->LoadTree(readentry);
         }
      }

   }

   if(savedir) savedir->cd();
}

////////////////////////////////////////////////////////////////////////////////
/// Tree Formula default destructor.

TTreeFormula::~TTreeFormula()
{
   if (fManager) {
      fManager->Remove(this);
      if (fManager->fFormulas.GetLast()<0) {
         delete fManager;
         fManager = 0;
      }
   }
   // Objects in fExternalCuts are not owned and should not be deleted
   // fExternalCuts.Clear();
   fLeafNames.Delete();
   fDataMembers.Delete();
   fMethods.Delete();
   fAliases.Delete();
   if (fLookupType) delete [] fLookupType;
   for (int j=0; j<fNcodes; j++) {
      for (int k = 0; k<fNdimensions[j]; k++) {
         if (fVarIndexes[j][k]) delete fVarIndexes[j][k];
         fVarIndexes[j][k] = 0;
      }
   }
   if (fDimensionSetup) {
      fDimensionSetup->Delete();
      delete fDimensionSetup;
   }
   delete[] fConstLD;
}

////////////////////////////////////////////////////////////////////////////////
/// This method is used internally to decode the dimensions of the variables.

void TTreeFormula::DefineDimensions(Int_t code, Int_t size,
                                    TFormLeafInfoMultiVarDim * info,
                                    Int_t& virt_dim) {
   if (info) {
      fManager->EnableMultiVarDims();
      //if (fIndexes[code][info->fDim]<0) { // removed because the index might be out of bounds!
         info->fVirtDim = virt_dim;
         fManager->AddVarDims(virt_dim); // if (!fVarDims[virt_dim]) fVarDims[virt_dim] = new TArrayI;
      //}
   }

   Int_t vsize = 0;
   bool scalarindex = false;

   if (fIndexes[code][fNdimensions[code]]==-2) {
      TTreeFormula *indexvar = fVarIndexes[code][fNdimensions[code]];
      // ASSERT(indexvar!=0);
      Int_t index_multiplicity = indexvar->GetMultiplicity();
      switch (index_multiplicity) {
         case  0:
            scalarindex = true;
            vsize = 1;
            break;
         case -1:
         case  2:
            vsize = indexvar->GetNdata();
            break;
         case  1:
            vsize = -1;
            break;
      };
   } else vsize = size;

   fCumulSizes[code][fNdimensions[code]] = size;

   if ( !scalarindex && fIndexes[code][fNdimensions[code]] < 0 ) {
      fManager->UpdateUsedSize(virt_dim, vsize);
   }

   fNdimensions[code] ++;

}

////////////////////////////////////////////////////////////////////////////////
/// This method is used internally to decode the dimensions of the variables.

Int_t TTreeFormula::RegisterDimensions(const char *info, Int_t code)
{
   // We assume that there are NO white spaces in the info string
   const char * current;
   Int_t size, scanindex, vardim;

   current = info;
   vardim = 0;
   // the next value could be before the string but
   // that's okay because the next operation is ++
   // (this is to avoid (?) a if statement at the end of the
   // loop)
   if (current[0] != '[') current--;
   while (current) {
      current++;
      scanindex = sscanf(current,"%d",&size);
      // if scanindex is 0 then we have a name index thus a variable
      // array (or TClonesArray!).

      if (scanindex==0) size = -1;

      vardim += RegisterDimensions(code, size);

      if (fNdimensions[code] >= kMAXFORMDIM) {
         // NOTE: test that fNdimensions[code] is NOT too big!!

         break;
      }
      current = (char*)strstr( current, "[" );
   }
   return vardim;
}


////////////////////////////////////////////////////////////////////////////////
/// This method stores the dimension information for later usage.

Int_t TTreeFormula::RegisterDimensions(Int_t code, Int_t size, TFormLeafInfoMultiVarDim * multidim) {
   TDimensionInfo * info = new TDimensionInfo(code,fNoper,size,multidim);
   fDimensionSetup->Add(info);
   fCumulSizes[code][fNdimensions[code]] = size;
   fNdimensions[code] ++;
   return (size==-1) ? 1 : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// This method is used internally to decode the dimensions of the variables.

Int_t TTreeFormula::RegisterDimensions(Int_t code, TFormLeafInfo *leafinfo,
                                       TFormLeafInfo * /* maininfo */,
                                       Bool_t useCollectionObject) {
   Int_t ndim, size, current, vardim;
   vardim = 0;

   const TStreamerElement * elem = leafinfo->fElement;
   TClass* c = elem ? elem->GetClassPointer() : 0;

   TFormLeafInfoMultiVarDim * multi = dynamic_cast<TFormLeafInfoMultiVarDim * >(leafinfo);
   if (multi) {
      // We have a second variable dimensions
      fManager->EnableMultiVarDims();
      multi->fDim = fNdimensions[code];
      return RegisterDimensions(code, -1, multi);
   }
   if (elem->IsA() == TStreamerBasicPointer::Class()) {

      if (elem->GetArrayDim()>0) {

         ndim = elem->GetArrayDim();
         size = elem->GetMaxIndex(0);
         vardim += RegisterDimensions(code, -1);
      } else {
         ndim = 1;
         size = -1;
      }

      TStreamerBasicPointer *array = (TStreamerBasicPointer*)elem;
      TClass *cl = leafinfo->fClass;
      Int_t offset;
      TStreamerElement* counter = ((TStreamerInfo*)cl->GetStreamerInfo())->GetStreamerElement(array->GetCountName(),offset);
#if 1
      leafinfo->fCounter = new TFormLeafInfo(cl,offset,counter);
#else /* Code is not ready yet see revision 14078 */
      if (maininfo==0 || maininfo==leafinfo || 1) {
         leafinfo->fCounter = new TFormLeafInfo(cl,offset,counter);
      } else {
         leafinfo->fCounter = maininfo->DeepCopy();
         TFormLeafInfo *currentinfo = leafinfo->fCounter;
         while(currentinfo->fNext && currentinfo->fNext->fNext) currentinfo=currentinfo->fNext;
         delete currentinfo->fNext;
         currentinfo->fNext = new TFormLeafInfo(cl,offset,counter);
      }
#endif
   } else if (!useCollectionObject && elem->GetClassPointer() == TClonesArray::Class() ) {

      ndim = 1;
      size = -1;

      TClass * clonesClass = TClonesArray::Class();
      Int_t c_offset;
      TStreamerElement *counter = ((TStreamerInfo*)clonesClass->GetStreamerInfo())->GetStreamerElement("fLast",c_offset);
      leafinfo->fCounter = new TFormLeafInfo(clonesClass,c_offset,counter);

   } else if (!useCollectionObject && elem->GetClassPointer() && elem->GetClassPointer()->GetCollectionProxy() ) {

      if ( typeid(*leafinfo) == typeid(TFormLeafInfoCollection) ) {
         ndim = 1;
         size = -1;
      } else {
         R__ASSERT( fHasMultipleVarDim[code] );
         ndim = 1;
         size = 1;
      }

   } else if ( c && c->GetReferenceProxy() && c->GetReferenceProxy()->HasCounter() ) {
      ndim = 1;
      size = -1;
   } else if (elem->GetArrayDim()>0) {

      ndim = elem->GetArrayDim();
      size = elem->GetMaxIndex(0);

   } else if ( elem->GetNewType()== TStreamerInfo::kCharStar) {

      // When we implement being able to read the length from
      // strlen, we will have:
      // ndim = 1;
      // size = -1;
      // until then we more or so die:
      ndim = 1;
      size = 1; //NOTE: changed from 0

   } else return 0;

   current = 0;
   do {
      vardim += RegisterDimensions(code, size);

      if (fNdimensions[code] >= kMAXFORMDIM) {
         // NOTE: test that fNdimensions[code] is NOT too big!!

         break;
      }
      current++;
      size = elem->GetMaxIndex(current);
   } while (current<ndim);

   return vardim;
}

////////////////////////////////////////////////////////////////////////////////
/// This method is used internally to decode the dimensions of the variables.

Int_t TTreeFormula::RegisterDimensions(Int_t code, TBranchElement *branch) {
   TBranchElement * leafcount2 = branch->GetBranchCount2();
   if (leafcount2) {
      // With have a second variable dimensions
      TBranchElement *leafcount = dynamic_cast<TBranchElement*>(branch->GetBranchCount());

      R__ASSERT(leafcount); // The function should only be called on a functional TBranchElement object

      fManager->EnableMultiVarDims();
      TFormLeafInfoMultiVarDim * info = new TFormLeafInfoMultiVarDimDirect();
      fDataMembers.AddAtAndExpand(info, code);
      fHasMultipleVarDim[code] = kTRUE;

      info->fCounter = new TFormLeafInfoDirect(leafcount);
      info->fCounter2 = new TFormLeafInfoDirect(leafcount2);
      info->fDim = fNdimensions[code];
      //if (fIndexes[code][info->fDim]<0) {
      //  info->fVirtDim = virt_dim;
      //  if (!fVarDims[virt_dim]) fVarDims[virt_dim] = new TArrayI;
      //}
      return RegisterDimensions(code, -1, info);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// This method is used internally to decode the dimensions of the variables.

Int_t TTreeFormula::RegisterDimensions(Int_t code, TLeaf *leaf) {
   Int_t numberOfVarDim = 0;

   // Let see if we can understand the structure of this branch.
   // Usually we have: leafname[fixed_array] leaftitle[var_array]\type
   // (with fixed_array that can be a multi-dimension array.
   const char *tname = leaf->GetTitle();
   char *leaf_dim = (char*)strstr( tname, "[" );

   const char *bname = leaf->GetBranch()->GetName();
   char *branch_dim = (char*)strstr(bname,"[");
   if (branch_dim) branch_dim++; // skip the '['

   Bool_t isString  = kFALSE;
   if  (leaf->IsA() == TLeafElement::Class()) {
      Int_t type =((TBranchElement*)leaf->GetBranch())->GetStreamerType();
      isString =    (type == TStreamerInfo::kOffsetL+TStreamerInfo::kChar)
                 || (type == TStreamerInfo::kCharStar);
   } else {
      isString = (leaf->IsA() == TLeafC::Class());
   }
   if (leaf_dim) {
      leaf_dim++; // skip the '['
      if (!branch_dim || strncmp(branch_dim,leaf_dim,strlen(branch_dim))) {
         // then both are NOT the same so do the leaf title first:
         numberOfVarDim += RegisterDimensions( leaf_dim, code);
      } else if (branch_dim && strncmp(branch_dim,leaf_dim,strlen(branch_dim))==0
                 && strlen(leaf_dim)>strlen(branch_dim)
                 && (leaf_dim+strlen(branch_dim))[0]=='[') {
         // we have extra info in the leaf title
         numberOfVarDim += RegisterDimensions( leaf_dim+strlen(branch_dim)+1, code);
      }
   }
   if (branch_dim) {
      // then both are NOT same so do the branch name next:
      if (isString) {
         numberOfVarDim += RegisterDimensions( code, 1);
      } else {
         numberOfVarDim += RegisterDimensions( branch_dim, code);
      }
   }
   if (leaf->IsA() == TLeafElement::Class()) {
      TBranchElement* branch = (TBranchElement*) leaf->GetBranch();
      if (branch->GetBranchCount2()) {

         if (!branch->GetBranchCount()) {
            Warning("DefinedVariable",
                    "Noticed an incorrect in-memory TBranchElement object (%s).\nIt has a BranchCount2 but no BranchCount!\nThe result might be incorrect!",
                    branch->GetName());
            return numberOfVarDim;
         }

         // Switch from old direct style to using a TLeafInfo
         if (fLookupType[code] == kDataMember)
            Warning("DefinedVariable",
                    "Already in kDataMember mode when handling multiple variable dimensions");
         fLookupType[code] = kDataMember;

         // Feed the information into the Dimensions system
         numberOfVarDim += RegisterDimensions( code, branch);

      }
   }
   return numberOfVarDim;
}

////////////////////////////////////////////////////////////////////////////////
/// This method check for treat the case where expression contains $Atl and load up
/// both fAliases and fExpr.
/// We return:
/// -  -1 in case of failure
/// -  0 in case we did not find $Alt
/// -  the action number in case of success.

Int_t TTreeFormula::DefineAlternate(const char *expression)
{
   static const char *altfunc = "Alt$(";
   static const char *minfunc = "MinIf$(";
   static const char *maxfunc = "MaxIf$(";
   Int_t action = 0;
   Int_t start = 0;

   if (   strncmp(expression,altfunc,strlen(altfunc))==0
       && expression[strlen(expression)-1]==')' ) {
      action = kAlternate;
      start = strlen(altfunc);
   }
   if (   strncmp(expression,maxfunc,strlen(maxfunc))==0
       && expression[strlen(expression)-1]==')' ) {
      action = kMaxIf;
      start = strlen(maxfunc);
   }
   if (   strncmp(expression,minfunc,strlen(minfunc))==0
       && expression[strlen(expression)-1]==')' ) {
      action = kMinIf;
      start = strlen(minfunc);
   }

   if (action) {
      TString full = expression;
      TString part1;
      TString part2;
      int paran = 0;
      int instr = 0;
      int brack = 0;
      for(unsigned int i=start;i<strlen(expression);++i) {
         switch (expression[i]) {
            case '(': paran++; break;
            case ')': paran--; break;
            case '"': instr = instr ? 0 : 1; break;
            case '[': brack++; break;
            case ']': brack--; break;
         };
         if (expression[i]==',' && paran==0 && instr==0 && brack==0) {
            part1 = full( start, i-start );
            part2 = full( i+1, full.Length() -1 - (i+1) );
            break; // out of the for loop
         }
      }
      if (part1.Length() && part2.Length()) {
         TTreeFormula *primary = new TTreeFormula("primary",part1,fTree);
         TTreeFormula *alternate = new TTreeFormula("alternate",part2,fTree);

         short isstring = 0;

         if (action == kAlternate) {
            if (alternate->GetManager()->GetMultiplicity() != 0 ) {
               Error("DefinedVariable","The 2nd arguments in %s can not be an array (%s,%d)!",
                     expression,alternate->GetTitle(),
                     alternate->GetManager()->GetMultiplicity());
               return -1;
            }

            // Should check whether we have strings.
            if (primary->IsString()) {
               if (!alternate->IsString()) {
                  Error("DefinedVariable",
                        "The 2nd arguments in %s has to return the same type as the 1st argument (string)!",
                        expression);
                  return -1;
               }
               isstring = 1;
            } else if (alternate->IsString()) {
               Error("DefinedVariable",
                     "The 2nd arguments in %s has to return the same type as the 1st argument (numerical type)!",
                     expression);
               return -1;
            }
         } else {
            primary->GetManager()->Add( alternate );
            primary->GetManager()->Sync();
            if (primary->IsString() || alternate->IsString()) {
               if (!alternate->IsString()) {
                  Error("DefinedVariable",
                        "The arguments of %s can not be strings!",
                        expression);
                  return -1;
               }
            }
         }

         fAliases.AddAtAndExpand(primary,fNoper);
         fExpr[fNoper] = "";
         SetAction(fNoper, (Int_t)action + isstring, 0 );
         ++fNoper;

         fAliases.AddAtAndExpand(alternate,fNoper);
         return (Int_t)kAlias + isstring;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Decompose 'expression' as pointing to something inside the leaf
/// Returns:
/// -  -2  Error: some information is missing (message already printed)
/// -  -1  Error: Syntax is incorrect (message already printed)
/// -  0
/// -  >0 the value returns is the action code.

Int_t TTreeFormula::ParseWithLeaf(TLeaf* leaf, const char* subExpression, Bool_t final, UInt_t paran_level, TObjArray& castqueue, Bool_t useLeafCollectionObject, const char* fullExpression)
{
   Int_t action = 0;

   Int_t numberOfVarDim = 0;
   char *current;

   char  scratch[kMaxLen]; scratch[0] = '\0';
   char     work[kMaxLen];    work[0] = '\0';

   const char *right = subExpression;
   TString name = fullExpression;

   TBranch *branch = leaf ? leaf->GetBranch() : 0;
   Long64_t readentry = fTree->GetTree()->GetReadEntry();
   if (readentry < 0) readentry=0;

   Bool_t useLeafReferenceObject = false;
   Int_t code = fNcodes-1;

   // Make a check to prevent problem with some corrupted files (missing TStreamerInfo).
   if (leaf && leaf->IsA()==TLeafElement::Class()) {
      TBranchElement *br = 0;
      if( branch->IsA() ==  TBranchElement::Class() )
      {
         br = ((TBranchElement*)branch);

         if ( br->GetInfo() == 0 ) {
            Error("DefinedVariable","Missing StreamerInfo for %s.  We will be unable to read!",
                  name.Data());
            return -2;
         }
      }

      TBranch *bmom = branch->GetMother();
      if( bmom->IsA() == TBranchElement::Class() )
      {
         TBranchElement *mom = (TBranchElement*)br->GetMother();
         if (mom!=br) {
            if (mom->GetInfo()==0) {
               Error("DefinedVariable","Missing StreamerInfo for %s."
                     "  We will be unable to read!",
                     mom->GetName());
               return -2;
            }
            if ((mom->GetType()) < -1 && !mom->GetAddress()) {
               Error("DefinedVariable", "Address not set when the type of the branch is negative for for %s.  We will be unable to read!", mom->GetName());
               return -2;
            }
         }
      }
   }

   // We need to record the location in the list of leaves because
   // the tree might actually be a chain and in that case the leaf will
   // change from tree to tree!.

   // Let's reconstruct the name of the leaf, including the possible friend alias
   TTree *realtree = fTree->GetTree();
   const char* alias = 0;
   if (leaf) {
      if (realtree) alias = realtree->GetFriendAlias(leaf->GetBranch()->GetTree());
      if (!alias && realtree!=fTree) {
         // Let's try on the chain
         alias = fTree->GetFriendAlias(leaf->GetBranch()->GetTree());
      }
   }
   if (alias) snprintf(scratch,kMaxLen-1,"%s.%s",alias,leaf->GetName());
   else if (leaf) strlcpy(scratch,leaf->GetName(),kMaxLen);

   TTree *tleaf = realtree;
   if (leaf) {
      tleaf = leaf->GetBranch()->GetTree();
      fCodes[code] = tleaf->GetListOfLeaves()->IndexOf(leaf);
      const char *mother_name = leaf->GetBranch()->GetMother()->GetName();
      TString br_extended_name; // Could do ( strlen(mother_name)+strlen( leaf->GetBranch()->GetName() ) + 2 )
      if (leaf->GetBranch()!=leaf->GetBranch()->GetMother()) {
         if (mother_name[strlen(mother_name)-1]!='.') {
            br_extended_name = mother_name;
            br_extended_name.Append('.');
         }
      }
      br_extended_name.Append( leaf->GetBranch()->GetName() );
      Ssiz_t dim = br_extended_name.First('[');
      if (dim >= 0) br_extended_name.Remove(dim);

      TNamed *named = new TNamed(scratch,br_extended_name.Data());
      fLeafNames.AddAtAndExpand(named,code);
      fLeaves.AddAtAndExpand(leaf,code);
   }

   // If the leaf belongs to a friend tree which has an index, we might
   // be in the case where some entry do not exist.
   if (tleaf != realtree && tleaf->GetTreeIndex()) {
      // reset the multiplicity
      if (fMultiplicity >= 0) fMultiplicity = 1;
   }

   // Analyze the content of 'right'

   // Try to find out the class (if any) of the object in the leaf.
   TClass * cl = 0;
   TFormLeafInfo *maininfo = 0;
   TFormLeafInfo *previnfo = 0;
   Bool_t unwindCollection = kFALSE;
   const static TClassRef stdStringClass = TClass::GetClass("string");

   if (leaf==0) {
      TNamed *names = (TNamed*)fLeafNames.UncheckedAt(code);
      fLeafNames.AddAt(0,code);
      TTree *what = (TTree*)fLeaves.UncheckedAt(code);
      fLeaves.AddAt(0,code);

      cl = what ? what->IsA() : TTree::Class();
      maininfo = new TFormLeafInfoTTree(fTree,names->GetName(),what);
      previnfo = maininfo;

      delete names;
   } else if (leaf->InheritsFrom(TLeafObject::Class()) ) {
      TBranchObject *bobj = (TBranchObject*)leaf->GetBranch();
      cl = TClass::GetClass(bobj->GetClassName());
   } else if (leaf->InheritsFrom(TLeafElement::Class())) {
      TBranchElement *branchEl = (TBranchElement *)leaf->GetBranch();
      branchEl->SetupAddresses();
      TStreamerInfo *info = branchEl->GetInfo();
      TStreamerElement *element = 0;
      Int_t type = branchEl->GetStreamerType();
      switch(type) {
         case TStreamerInfo::kBase:
         case TStreamerInfo::kObject:
         case TStreamerInfo::kTString:
         case TStreamerInfo::kTNamed:
         case TStreamerInfo::kTObject:
         case TStreamerInfo::kAny:
         case TStreamerInfo::kAnyP:
         case TStreamerInfo::kAnyp:
         case TStreamerInfo::kSTL:
         case TStreamerInfo::kSTLp:
         case TStreamerInfo::kObjectp:
         case TStreamerInfo::kObjectP: {
            element = info->GetElement(branchEl->GetID());
            if (element) cl = element->GetClassPointer();
         }
         break;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kSTL:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kSTLp:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kAny:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kAnyp:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kAnyP:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kObjectp:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kObjectP:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kObject:  {
            element = info->GetElement(branchEl->GetID());
            if (element){
               cl = element->GetClassPointer();
            }
         }
         break;
         case -1: {
            cl = info->GetClass();
         }
         break;
      }

      // If we got a class object, we need to verify whether it is on a
      // split TClonesArray sub branch.
      if (cl && branchEl->GetBranchCount()) {
         if (branchEl->GetType()==31) {
            // This is inside a TClonesArray.

            if (!element) {
               Warning("DefinedVariable",
                       "Missing TStreamerElement in object in TClonesArray section");
               return -2;
            }
            TFormLeafInfo* clonesinfo = new TFormLeafInfoClones(cl, 0, element, kTRUE);

            // The following code was commented out because in THIS case
            // the dimension are actually handled by parsing the title and name of the leaf
            // and branch (see a little further)
            // The dimension needs to be handled!
            // numberOfVarDim += RegisterDimensions(code,clonesinfo);

            maininfo = clonesinfo;

            // We skip some cases because we can assume we have an object.
            Int_t offset=0;
            info->GetStreamerElement(element->GetName(),offset);
            if (type == TStreamerInfo::kObjectp ||
                  type == TStreamerInfo::kObjectP ||
                  type == TStreamerInfo::kOffsetL + TStreamerInfo::kObjectp ||
                  type == TStreamerInfo::kOffsetL + TStreamerInfo::kObjectP ||
                  type == TStreamerInfo::kSTLp ||
                  type == TStreamerInfo::kAnyp ||
                  type == TStreamerInfo::kAnyP ||
                  type == TStreamerInfo::kOffsetL + TStreamerInfo::kSTLp ||
                  type == TStreamerInfo::kOffsetL + TStreamerInfo::kAnyp ||
                  type == TStreamerInfo::kOffsetL + TStreamerInfo::kAnyP) {
               previnfo = new TFormLeafInfoPointer(cl,offset+branchEl->GetOffset(),element);
            } else {
               previnfo = new TFormLeafInfo(cl,offset+branchEl->GetOffset(),element);
            }
            maininfo->fNext = previnfo;
            unwindCollection = kTRUE;

         } else if (branchEl->GetType()==41) {

            // This is inside a Collection

            if (!element) {
               Warning("DefinedVariable","Missing TStreamerElement in object in Collection section");
               return -2;
            }
            // First we need to recover the collection.
            TBranchElement *count = branchEl->GetBranchCount();
            TFormLeafInfo* collectioninfo;
            if ( count->GetID() >= 0 ) {
               TStreamerElement *collectionElement =
                  count->GetInfo()->GetElement(count->GetID());
               TClass *collectionCl = collectionElement->GetClassPointer();

               collectioninfo =
                  new TFormLeafInfoCollection(collectionCl, 0, collectionElement, kTRUE);
            } else {
               TClass *collectionCl = TClass::GetClass(count->GetClassName());
               collectioninfo =
                  new TFormLeafInfoCollection(collectionCl, 0, collectionCl, kTRUE);
            }

            // The following code was commented out because in THIS case
            // the dimension are actually handled by parsing the title and name of the leaf
            // and branch (see a little further)
            // The dimension needs to be handled!
            // numberOfVarDim += RegisterDimensions(code,clonesinfo);

            maininfo = collectioninfo;

            // We skip some cases because we can assume we have an object.
            Int_t offset=0;
            info->GetStreamerElement(element->GetName(),offset);
            if (type == TStreamerInfo::kObjectp ||
                  type == TStreamerInfo::kObjectP ||
                  type == TStreamerInfo::kOffsetL + TStreamerInfo::kObjectp ||
                  type == TStreamerInfo::kOffsetL + TStreamerInfo::kObjectP ||
                  type == TStreamerInfo::kSTLp ||
                  type == TStreamerInfo::kAnyp ||
                  type == TStreamerInfo::kAnyP ||
                  type == TStreamerInfo::kOffsetL + TStreamerInfo::kSTLp ||
                  type == TStreamerInfo::kOffsetL + TStreamerInfo::kAnyp ||
                  type == TStreamerInfo::kOffsetL + TStreamerInfo::kAnyP) {
               previnfo = new TFormLeafInfoPointer(cl,offset+branchEl->GetOffset(),element);
            } else {
               previnfo = new TFormLeafInfo(cl,offset+branchEl->GetOffset(),element);
            }
            maininfo->fNext = previnfo;
            unwindCollection = kTRUE;
         }
      } else if ( branchEl->GetType()==3) {
         TFormLeafInfo* clonesinfo;
         if (useLeafCollectionObject) {
            clonesinfo = new TFormLeafInfoCollectionObject(cl);
         } else {
            clonesinfo = new TFormLeafInfoClones(cl, 0, kTRUE);
            // The dimension needs to be handled!
            numberOfVarDim += RegisterDimensions(code,clonesinfo,maininfo,useLeafCollectionObject);

         }
         maininfo = clonesinfo;
         previnfo = maininfo;

      } else if (!useLeafCollectionObject && branchEl->GetType()==4) {

         TFormLeafInfo* collectioninfo;
         if (useLeafCollectionObject) {
            collectioninfo = new TFormLeafInfoCollectionObject(cl);
         } else {
            collectioninfo = new TFormLeafInfoCollection(cl, 0, cl, kTRUE);
            // The dimension needs to be handled!
            numberOfVarDim += RegisterDimensions(code,collectioninfo,maininfo,useLeafCollectionObject);
         }

         maininfo = collectioninfo;
         previnfo = maininfo;

      } else if (branchEl->GetStreamerType()==-1 && cl && cl->GetCollectionProxy()) {

         if (useLeafCollectionObject) {

            TFormLeafInfo *collectioninfo = new TFormLeafInfoCollectionObject(cl);
            maininfo = collectioninfo;
            previnfo = collectioninfo;

         } else {
            TFormLeafInfo *collectioninfo = new TFormLeafInfoCollection(cl, 0, cl, kTRUE);
            // The dimension needs to be handled!
            numberOfVarDim += RegisterDimensions(code,collectioninfo,maininfo,kFALSE);

            maininfo = collectioninfo;
            previnfo = collectioninfo;

            if (cl->GetCollectionProxy()->GetValueClass()!=0 &&
                cl->GetCollectionProxy()->GetValueClass()->GetCollectionProxy()!=0) {

               TFormLeafInfo *multi = new TFormLeafInfoMultiVarDimCollection(cl,0,
                     cl->GetCollectionProxy()->GetValueClass(),collectioninfo);

               fHasMultipleVarDim[code] = kTRUE;
               numberOfVarDim += RegisterDimensions(code,multi,maininfo,kFALSE);
               previnfo->fNext = multi;
               cl = cl->GetCollectionProxy()->GetValueClass();
               multi->fNext =  new TFormLeafInfoCollection(cl, 0, cl, false);
               previnfo = multi->fNext;

            }
            if (cl->GetCollectionProxy()->GetValueClass()==0 &&
                cl->GetCollectionProxy()->GetType()>0) {

               previnfo->fNext =
                        new TFormLeafInfoNumerical(cl->GetCollectionProxy());
               previnfo = previnfo->fNext;
            } else {
               // nothing to do
            }
         }

      } else if (strlen(right)==0 && cl && element && final) {

         TClass *elemCl = element->GetClassPointer();
         if (!useLeafCollectionObject
             && elemCl && elemCl->GetCollectionProxy()
             && elemCl->GetCollectionProxy()->GetValueClass()
             && elemCl->GetCollectionProxy()->GetValueClass()->GetCollectionProxy()) {

            TFormLeafInfo *collectioninfo =
               new TFormLeafInfoCollection(cl, 0, elemCl);

            // The dimension needs to be handled!
            numberOfVarDim += RegisterDimensions(code,collectioninfo,maininfo,kFALSE);

            maininfo = collectioninfo;
            previnfo = collectioninfo;

            TFormLeafInfo *multi =
               new TFormLeafInfoMultiVarDimCollection(elemCl, 0,
                                                      elemCl->GetCollectionProxy()->GetValueClass(),
                                                      collectioninfo);

            fHasMultipleVarDim[code] = kTRUE;
            numberOfVarDim += RegisterDimensions(code,multi,maininfo,kFALSE);
            previnfo->fNext = multi;
            cl = elemCl->GetCollectionProxy()->GetValueClass();
            multi->fNext =  new TFormLeafInfoCollection(cl, 0, cl, false);
            previnfo = multi->fNext;

            if (cl->GetCollectionProxy()->GetValueClass()==0 &&
                cl->GetCollectionProxy()->GetType()>0) {

               previnfo->fNext =
                  new TFormLeafInfoNumerical(cl->GetCollectionProxy());
               previnfo = previnfo->fNext;
            }

         } else if (!useLeafCollectionObject
               && elemCl && elemCl->GetCollectionProxy()
               && elemCl->GetCollectionProxy()->GetValueClass()==0
               && elemCl->GetCollectionProxy()->GetType()>0) {

            // At this point we have an element which is inside a class (which is not
            // a collection) and this element of a collection of numerical type.
            // (Note: it is not possible to have more than one variable dimension
            //  unless we were supporting variable size C-style array of collection).

            TFormLeafInfo* collectioninfo =
               new TFormLeafInfoCollection(cl, 0, elemCl);

            // The dimension needs to be handled!
            numberOfVarDim += RegisterDimensions(code,collectioninfo,maininfo,kFALSE);

            collectioninfo->fNext =
               new TFormLeafInfoNumerical(elemCl->GetCollectionProxy());

            maininfo = collectioninfo;
            previnfo = maininfo->fNext;

         }  else if (!useLeafCollectionObject
                     && elemCl && elemCl->GetCollectionProxy()) {
            if (elemCl->GetCollectionProxy()->GetValueClass()==TString::Class()) {
               right = "Data()";
            } else if (elemCl->GetCollectionProxy()->GetValueClass()==stdStringClass) {
               right = "c_str()";
            }

         } else if (!element->IsaPointer()) {

            maininfo = new TFormLeafInfoDirect(branchEl);
            previnfo = maininfo;

         }
      }
      else if ( cl && cl->GetReferenceProxy() )  {
         if ( useLeafCollectionObject || fullExpression[0] == '@' || fullExpression[strlen(scratch)] == '@' ) {
            useLeafReferenceObject = true;
         }
         else  {
            if ( !maininfo )  {
               maininfo = previnfo = new TFormLeafInfoReference(cl, element, 0);
               numberOfVarDim += RegisterDimensions(code,maininfo,maininfo,kFALSE);
            }
            TVirtualRefProxy *refproxy = cl->GetReferenceProxy();
            for(Long64_t i=0; i<leaf->GetBranch()->GetEntries()-readentry; ++i)  {
               R__LoadBranch(leaf->GetBranch(), readentry+i, fQuickLoad);
               void *refobj = maininfo->GetValuePointer(leaf,0);
               if (refobj) {
                  cl = refproxy->GetValueClass(refobj);
               }
               if ( cl ) break;
            }
            if ( !cl )  {
               Error("DefinedVariable","Failed to access class type of reference target (%s)",element->GetName());
               return -1;
            }
         }
      }
   } else {
      // Regular/old TLeaf, there should not be anything afterward ...
      if (subExpression && subExpression[0]) {
        Error("ParseWithLeaf", "Found a numerical leaf but the name has trailing characters: \"%s\"", subExpression);
        return -1;
      }
   }

   // Treat the dimension information in the leaf name, title and 2nd branch count
   if (leaf) numberOfVarDim += RegisterDimensions(code,leaf);

   if (cl) {
      if (unwindCollection) {
         // So far we should get here only if we encounter a split collection of a class that contains
         // directly a collection.
         R__ASSERT(numberOfVarDim==1 && maininfo);

         if (!useLeafCollectionObject && cl && cl->GetCollectionProxy()) {
            TFormLeafInfo *multi =
               new TFormLeafInfoMultiVarDimCollection(cl, 0, cl, maininfo);
            fHasMultipleVarDim[code] = kTRUE;
            numberOfVarDim += RegisterDimensions(code,multi,maininfo,kFALSE);
            previnfo->fNext = multi;

            multi->fNext =  new TFormLeafInfoCollection(cl, 0, cl, false);
            previnfo = multi->fNext;

            if (cl->GetCollectionProxy()->GetValueClass()==0 &&
                cl->GetCollectionProxy()->GetType()>0) {

               previnfo->fNext =
                  new TFormLeafInfoNumerical(cl->GetCollectionProxy());
               previnfo = previnfo->fNext;
            }
         } else if (!useLeafCollectionObject && cl == TClonesArray::Class()) {

            TFormLeafInfo *multi =
               new TFormLeafInfoMultiVarDimClones(cl, 0, cl, maininfo);
            fHasMultipleVarDim[code] = kTRUE;
            numberOfVarDim += RegisterDimensions(code,multi,maininfo,kFALSE);
            previnfo->fNext = multi;

            multi->fNext =  new TFormLeafInfoClones(cl, 0, false);
            previnfo = multi->fNext;
         }
      }
      Int_t offset=0;
      if (cl == TString::Class() && strcmp(right,"fData")==0) {
         // For backward compatibility replace TString::fData which no longer exist
         // by a call to TString::Data()
         right = "Data()";
      }
      Int_t nchname = strlen(right);
      TFormLeafInfo *leafinfo = 0;
      TStreamerElement* element = 0;

      // Let see if the leaf was attempted to be casted.
      // Since there would have been something like
      // ((cast_class*)leafname)->....  we need to use
      // paran_level+1
      // Also we disable this functionality in case of TClonesArray
      // because it is not yet allowed to have 'inheritance' (or virtuality)
      // in play in a TClonesArray.
      {
         TClass * casted = (TClass*) castqueue.At(paran_level+1);
         if (casted && cl != TClonesArray::Class()) {
            if ( ! casted->InheritsFrom(cl) ) {
               Error("DefinedVariable","%s does not inherit from %s.  Casting not possible!",
                     casted->GetName(),cl->GetName());
               return -2;
            }
            leafinfo = new TFormLeafInfoCast(cl,casted);
            fHasCast = kTRUE;
            if (maininfo==0) {
               maininfo = leafinfo;
            }
            if (previnfo==0) {
               previnfo = leafinfo;
            } else {
               previnfo->fNext = leafinfo;
               previnfo = leafinfo;
            }
            leafinfo = 0;

            cl = casted;
            castqueue.AddAt(0,paran_level);
         }
      }
      Int_t i;
      Bool_t prevUseCollectionObject = useLeafCollectionObject;
      Bool_t useCollectionObject = useLeafCollectionObject;
      Bool_t useReferenceObject = useLeafReferenceObject;
      Bool_t prevUseReferenceObject = useLeafReferenceObject;
      for (i=0, current = &(work[0]); i<=nchname;i++ ) {
         // We will treated the terminator as a token.
         if (right[i] == '(') {
            // Right now we do not allow nested parenthesis
            do {
               *current++ = right[i++];
            } while(right[i]!=')' && right[i]);
            *current++ = right[i];
            *current='\0';
            char *params = strchr(work,'(');
            if (params) {
               *params = 0; params++;
            } else params = (char *) ")";
            if (cl==0) {
               Error("DefinedVariable","Can not call '%s' with a class",work);
               return -1;
            }
            if (!cl->HasDataMemberInfo() && !cl->GetCollectionProxy()) {
               Error("DefinedVariable","Class probably unavailable:%s",cl->GetName());
               return -2;
            }
            if (!useCollectionObject && cl == TClonesArray::Class()) {
               // We are not interested in the ClonesArray object but only
               // in its contents.
               // We need to retrieve the class of its content.

               TBranch *clbranch = leaf->GetBranch();
               R__LoadBranch(clbranch,readentry,fQuickLoad);
               TClonesArray * clones;
               if (previnfo) clones = (TClonesArray*)previnfo->GetLocalValuePointer(leaf,0);
               else {
                  Bool_t top = (clbranch==((TBranchElement*)clbranch)->GetMother()
                                 || !leaf->IsOnTerminalBranch());
                  TClass *mother_cl;
                  if (leaf->IsA()==TLeafObject::Class()) {
                     // in this case mother_cl is not really used
                     mother_cl = cl;
                  } else {
                     mother_cl = ((TBranchElement*)clbranch)->GetInfo()->GetClass();
                  }
                  TFormLeafInfo* clonesinfo = new TFormLeafInfoClones(mother_cl, 0, top);

                  // The dimension needs to be handled!
                  numberOfVarDim += RegisterDimensions(code,clonesinfo,maininfo,kFALSE);

                  previnfo = clonesinfo;
                  maininfo = clonesinfo;

                  clones = (TClonesArray*)clonesinfo->GetLocalValuePointer(leaf,0);
               }
               TClass * inside_cl = clones->GetClass();
               cl = inside_cl;

            }
            else if (!useCollectionObject && cl && cl->GetCollectionProxy() ) {

               // We are NEVER (for now!) interested in the ClonesArray object but only
               // in its contents.
               // We need to retrieve the class of its content.

               if (previnfo==0) {

                  Bool_t top = (branch==((TBranchElement*)branch)->GetMother()
                                 || !leaf->IsOnTerminalBranch());

                  TClass *mother_cl;
                  if (leaf->IsA()==TLeafObject::Class()) {
                     // in this case mother_cl is not really used
                     mother_cl = cl;
                  } else {
                     mother_cl = ((TBranchElement*)branch)->GetInfo()->GetClass();
                  }

                  TFormLeafInfo* collectioninfo =
                     new TFormLeafInfoCollection(mother_cl, 0,cl,top);
                  // The dimension needs to be handled!
                  numberOfVarDim += RegisterDimensions(code,collectioninfo,maininfo,kFALSE);

                  previnfo = collectioninfo;
                  maininfo = collectioninfo;

               }

               TClass * inside_cl = cl->GetCollectionProxy()->GetValueClass();
               if (inside_cl) cl = inside_cl;
               else if (cl->GetCollectionProxy()->GetType()>0) {
                  Warning("DefinedVariable","Can not call method on content of %s in %s\n",
                           cl->GetName(),name.Data());
                  return -2;
               }
            }
            TMethodCall *method  = 0;
            if (cl==0) {
               Error("DefinedVariable",
                     "Could not discover the TClass corresponding to (%s)!",
                     right);
               return -2;
            } else if (cl==TClonesArray::Class() && strcmp(work,"size")==0) {
               method = new TMethodCall(cl, "GetEntriesFast", "");
            } else if (cl->GetCollectionProxy() && strcmp(work,"size")==0) {
               if (maininfo==0) {
                  TFormLeafInfo* collectioninfo=0;
                  if (useLeafCollectionObject) {

                     Bool_t top = (branch==((TBranchElement*)branch)->GetMother()
                                 || !leaf->IsOnTerminalBranch());
                     collectioninfo = new TFormLeafInfoCollectionObject(cl,top);
                  }
                  maininfo=previnfo=collectioninfo;
               }
               leafinfo = new TFormLeafInfoCollectionSize(cl);
               cl = 0;
            } else {
               if (!cl->HasDataMemberInfo()) {
                  Error("DefinedVariable",
                        "Can not call method %s on class without dictionary (%s)!",
                        right,cl->GetName());
                  return -2;
               }
               method = new TMethodCall(cl, work, params);
            }
            if (method) {
               if (!method->GetMethod()) {
                  Error("DefinedVariable","Unknown method:%s in %s",right,cl->GetName());
                  return -1;
               }
               switch(method->ReturnType()) {
                  case TMethodCall::kLong:
                     leafinfo = new TFormLeafInfoMethod(cl,method);
                     cl = 0;
                     break;
                  case TMethodCall::kDouble:
                     leafinfo = new TFormLeafInfoMethod(cl,method);
                     cl = 0;
                     break;
                  case TMethodCall::kString:
                     leafinfo = new TFormLeafInfoMethod(cl,method);
                     // 1 will be replaced by -1 when we know how to use strlen
                     numberOfVarDim += RegisterDimensions(code,1); //NOTE: changed from 0
                     cl = 0;
                     break;
                  case TMethodCall::kOther:
                     {
                        leafinfo = new TFormLeafInfoMethod(cl,method);
                        cl = TFormLeafInfoMethod::ReturnTClass(method);
                     }
                     break;
                  default:
                     Error("DefineVariable","Method %s from %s has an impossible return type %d",
                           work,cl->GetName(), (Int_t)method->ReturnType());
                     return -2;
               }
            }
            if (maininfo==0) {
               maininfo = leafinfo;
            }
            if (previnfo==0) {
               previnfo = leafinfo;
            } else {
               previnfo->fNext = leafinfo;
               previnfo = leafinfo;
            }
            leafinfo = 0;
            current = &(work[0]);
            *current = 0;
            prevUseCollectionObject = kFALSE;
            prevUseReferenceObject = kFALSE;
            useCollectionObject = kFALSE;

            if (cl && cl->GetCollectionProxy()) {
               if (numberOfVarDim>1) {
                  Warning("DefinedVariable","TTreeFormula support only 2 level of variables size collections.  Assuming '@' notation for the collection %s.",
                     cl->GetName());
                  leafinfo = new TFormLeafInfo(cl,0,0);
                  useCollectionObject = kTRUE;
               } else if (numberOfVarDim==0) {
                  R__ASSERT(maininfo);
                  leafinfo = new TFormLeafInfoCollection(cl,0,cl);
                  numberOfVarDim += RegisterDimensions(code,leafinfo,maininfo,kFALSE);
               } else if (numberOfVarDim==1) {
                  R__ASSERT(maininfo);
                  leafinfo =
                     new TFormLeafInfoMultiVarDimCollection(cl,0,
                     (TStreamerElement*)0,maininfo);
                  previnfo->fNext = leafinfo;
                  previnfo = leafinfo;
                  leafinfo = new TFormLeafInfoCollection(cl,0,cl);

                  fHasMultipleVarDim[code] = kTRUE;
                  numberOfVarDim += RegisterDimensions(code,leafinfo,maininfo,kFALSE);
               }
               previnfo->fNext = leafinfo;
               previnfo = leafinfo;
               leafinfo = 0;
            }
            continue;
         } else if (right[i] == ')') {
            // We should have the end of a cast operator.  Let's introduce a TFormLeafCast
            // in the chain.
            TClass * casted = (TClass*) ((int(--paran_level)>=0) ? castqueue.At(paran_level) : 0);
            if (casted) {
               leafinfo = new TFormLeafInfoCast(cl,casted);
               fHasCast = kTRUE;

               if (maininfo==0) {
                  maininfo = leafinfo;
               }
               if (previnfo==0) {
                  previnfo = leafinfo;
               } else {
                  previnfo->fNext = leafinfo;
                  previnfo = leafinfo;
               }
               leafinfo = 0;
               current = &(work[0]);
               *current = 0;

               cl = casted;
               continue;

            }
         } else if (i > 0 && (right[i] == '.' || right[i] == '[' || right[i] == '\0') ) {
            // A delimiter happened let's see if what we have seen
            // so far does point to a data member.
            Bool_t needClass = kTRUE;
            *current = '\0';

            // skip it all if there is nothing to look at
            if (strlen(work)==0) continue;

            prevUseCollectionObject = useCollectionObject;
            prevUseReferenceObject = useReferenceObject;
            if (work[0]=='@') {
               useReferenceObject = kTRUE;
               useCollectionObject = kTRUE;
               Int_t l = 0;
               for(l=0;work[l+1]!=0;++l) work[l] = work[l+1];
               work[l] = '\0';
            } else if (work[strlen(work)-1]=='@') {
               useReferenceObject = kTRUE;
               useCollectionObject = kTRUE;
               work[strlen(work)-1] = '\0';
            } else {
               useReferenceObject = kFALSE;
               useCollectionObject = kFALSE;
            }

            Bool_t mustderef = kFALSE;
            if ( !prevUseReferenceObject && cl && cl->GetReferenceProxy() )  {
               R__LoadBranch(leaf->GetBranch(), readentry, fQuickLoad);
               if ( !maininfo )  {
                  maininfo = previnfo = new TFormLeafInfoReference(cl, element, offset);
                  if ( cl->GetReferenceProxy()->HasCounter() )  {
                     numberOfVarDim += RegisterDimensions(code,-1);
                  }
                  prevUseReferenceObject = kFALSE;
               } else {
                  previnfo->fNext = new TFormLeafInfoReference(cl, element, offset);
                  previnfo = previnfo->fNext;
               }
               TVirtualRefProxy *refproxy = cl->GetReferenceProxy();
               cl = 0;
               for(Long64_t entry=0; entry<leaf->GetBranch()->GetEntries()-readentry; ++entry)  {
                  R__LoadBranch(leaf->GetBranch(), readentry+i, fQuickLoad);
                  void *refobj = maininfo->GetValuePointer(leaf,0);
                  if (refobj) {
                     cl = refproxy->GetValueClass(refobj);
                  }
                  if ( cl ) break;
               }
               needClass = kFALSE;
               mustderef = kTRUE;
            }
            else if (!prevUseCollectionObject && cl == TClonesArray::Class()) {
               // We are not interested in the ClonesArray object but only
               // in its contents.
               // We need to retrieve the class of its content.

               TBranch *clbranch = leaf->GetBranch();
               R__LoadBranch(clbranch,readentry,fQuickLoad);
               TClonesArray * clones;
               if (maininfo) {
                  clones = (TClonesArray*)maininfo->GetValuePointer(leaf,0);
               } else {
                  // we have a unsplit TClonesArray leaves
                  // or we did not yet match any of the sub-branches!

                  TClass *mother_cl;
                  if (leaf->IsA()==TLeafObject::Class()) {
                     // in this case mother_cl is not really used
                     mother_cl = cl;
                  } else {
                     mother_cl = ((TBranchElement*)clbranch)->GetInfo()->GetClass();
                  }

                  TFormLeafInfo* clonesinfo = new TFormLeafInfoClones(mother_cl, 0);
                  // The dimension needs to be handled!
                  numberOfVarDim += RegisterDimensions(code,clonesinfo,maininfo,kFALSE);

                  mustderef = kTRUE;
                  previnfo = clonesinfo;
                  maininfo = clonesinfo;

                  if (clbranch->GetListOfBranches()->GetLast()>=0) {
                     if (clbranch->IsA() != TBranchElement::Class()) {
                        Error("DefinedVariable","Unimplemented usage of ClonesArray");
                        return -2;
                     }
                     //clbranch = ((TBranchElement*)clbranch)->GetMother();
                     clones = (TClonesArray*)((TBranchElement*)clbranch)->GetObject();
                  } else
                     clones = (TClonesArray*)clonesinfo->GetLocalValuePointer(leaf,0);
               }
               // NOTE clones can be zero!
               if (clones==0) {
                  Warning("DefinedVariable",
                          "TClonesArray object was not retrievable for %s!",
                          name.Data());
                  return -1;
               }
               TClass * inside_cl = clones->GetClass();
#if 1
               cl = inside_cl;
#else
/* Maybe we should make those test lead to warning messages */
               if (1 || inside_cl) cl = inside_cl;
               // if inside_cl is nul ... we have a problem of inconsistency :(
               if (0 && strlen(work)==0) {
                  // However in this case we have NO content :(
                  // so let get the number of objects
                  //strcpy(work,"fLast");
               }
#endif
            } else if (!prevUseCollectionObject && cl && cl->GetCollectionProxy() ) {

               // We are NEVER interested in the Collection object but only
               // in its contents.
               // We need to retrieve the class of its content.

               TBranch *clbranch = leaf->GetBranch();
               R__LoadBranch(clbranch,readentry,fQuickLoad);

               if (maininfo==0) {

                  // we have a unsplit Collection leaf
                  // or we did not yet match any of the sub-branches!

                  TClass *mother_cl;
                  if (leaf->IsA()==TLeafObject::Class()) {
                     // in this case mother_cl is not really used
                     mother_cl = cl;
                  } else {
                     mother_cl = ((TBranchElement*)clbranch)->GetInfo()->GetClass();
                  }

                  TFormLeafInfo* collectioninfo =
                     new TFormLeafInfoCollection(mother_cl, 0, cl);
                  // The dimension needs to be handled!
                  numberOfVarDim += RegisterDimensions(code,collectioninfo,maininfo,kFALSE);

                  mustderef = kTRUE;
                  previnfo = collectioninfo;
                  maininfo = collectioninfo;

               } //else if (clbranch->GetStreamerType()==0) {

               //}

               TClass * inside_cl = cl->GetCollectionProxy()->GetValueClass();

               if (!inside_cl) {
                  Error("DefinedVariable","Could you not find the inner class for %s with coll type = %d",
                        cl->GetName(),cl->GetCollectionProxy()->GetType());
               }
               if (!inside_cl && cl->GetCollectionProxy()->GetType() > 0) {
                  Warning("DefinedVariable","No data member in content of %s in %s\n",
                           cl->GetName(),name.Data());
               }
               cl = inside_cl;
               // if inside_cl is nul ... we have a problem of inconsistency.
            }

            if (!cl) {
               if (leaf) leaf->GetBranch()->Print();
               Warning("DefinedVariable","Missing class for %s!",name.Data());
            } else {
               element = ((TStreamerInfo*)cl->GetStreamerInfo())->GetStreamerElement(work,offset);
            }

            if (!element && !prevUseCollectionObject) {
               // We allow for looking for a data member inside a class inside
               // a TClonesArray without mentioning the TClonesArrays variable name
               TIter next( cl->GetStreamerInfo()->GetElements() );
               TStreamerElement * curelem;
               while ((curelem = (TStreamerElement*)next())) {
                  if (curelem->GetClassPointer() ==  TClonesArray::Class()) {
                     Int_t clones_offset = 0;
                     ((TStreamerInfo*)cl->GetStreamerInfo())->GetStreamerElement(curelem->GetName(),clones_offset);
                     TFormLeafInfo* clonesinfo =
                        new TFormLeafInfo(cl, clones_offset, curelem);
                     TClonesArray * clones;
                     R__LoadBranch(leaf->GetBranch(),readentry,fQuickLoad);

                     if (previnfo) {
                        previnfo->fNext = clonesinfo;
                        clones = (TClonesArray*)maininfo->GetValuePointer(leaf,0);
                        previnfo->fNext = 0;
                     } else {
                        clones = (TClonesArray*)clonesinfo->GetLocalValuePointer(leaf,0);
                     }

                     TClass *sub_cl = clones->GetClass();
                     if (sub_cl) element = ((TStreamerInfo*)sub_cl->GetStreamerInfo())->GetStreamerElement(work,offset);
                     delete clonesinfo;

                     if (element) {
                        leafinfo = new TFormLeafInfoClones(cl,clones_offset,curelem);
                        numberOfVarDim += RegisterDimensions(code,leafinfo,maininfo,kFALSE);
                        if (maininfo==0) maininfo = leafinfo;
                        if (previnfo==0) previnfo = leafinfo;
                        else {
                           previnfo->fNext = leafinfo;
                           previnfo = leafinfo;
                        }
                        leafinfo = 0;
                        cl = sub_cl;
                        break;
                     }
                  } else if (curelem->GetClassPointer() && curelem->GetClassPointer()->GetCollectionProxy()) {

                     Int_t coll_offset = 0;
                     ((TStreamerInfo*)cl->GetStreamerInfo())->GetStreamerElement(curelem->GetName(),coll_offset);

                     TClass *sub_cl =
                        curelem->GetClassPointer()->GetCollectionProxy()->GetValueClass();
                     if (sub_cl) {
                        element = ((TStreamerInfo*)sub_cl->GetStreamerInfo())->GetStreamerElement(work,offset);
                     }
                     if (element) {
                        if (numberOfVarDim>1) {
                           Warning("DefinedVariable","TTreeFormula support only 2 level of variables size collections.  Assuming '@' notation for the collection %s.",
                                   curelem->GetName());
                           leafinfo = new TFormLeafInfo(cl,coll_offset,curelem);
                           useCollectionObject = kTRUE;
                        } else if (numberOfVarDim==1) {
                           R__ASSERT(maininfo);
                           leafinfo =
                              new TFormLeafInfoMultiVarDimCollection(cl,coll_offset,
                                                                     curelem,maininfo);
                           fHasMultipleVarDim[code] = kTRUE;
                           leafinfo->fNext = new TFormLeafInfoCollection(cl,coll_offset,curelem);
                           numberOfVarDim += RegisterDimensions(code,leafinfo,maininfo,kFALSE);
                        } else {
                           leafinfo = new TFormLeafInfoCollection(cl,coll_offset,curelem);
                           numberOfVarDim += RegisterDimensions(code,leafinfo,maininfo,kFALSE);
                        }
                        if (maininfo==0) maininfo = leafinfo;
                        if (previnfo==0) previnfo = leafinfo;
                        else {
                           previnfo->fNext = leafinfo;
                           previnfo = leafinfo;
                        }
                        if (leafinfo->fNext) {
                           previnfo = leafinfo->fNext;
                        }
                        leafinfo = 0;
                        cl = sub_cl;
                        break;
                     }
                  }
               }

            }

            if (element) {
               Int_t type = element->GetNewType();
               if (type<60 && type!=0) {
                  // This is a basic type ...
                  if (numberOfVarDim>=1 && type>40) {
                     // We have a variable array within a variable array!
                     leafinfo = new TFormLeafInfoMultiVarDim(cl,offset,element,maininfo);
                     fHasMultipleVarDim[code] = kTRUE;
                  } else {
                     if (leafinfo && type<=40 ) {
                        leafinfo->AddOffset(offset,element);
                     } else {
                        leafinfo = new TFormLeafInfo(cl,offset,element);
                     }
                  }
               } else {
                  Bool_t object = kFALSE;
                  Bool_t pointer = kFALSE;
                  Bool_t objarr = kFALSE;
                  switch(type) {
                     case TStreamerInfo::kObjectp:
                     case TStreamerInfo::kObjectP:
                     case TStreamerInfo::kSTLp:
                     case TStreamerInfo::kAnyp:
                     case TStreamerInfo::kAnyP:
                     case TStreamerInfo::kOffsetL + TStreamerInfo::kObjectP:
                     case TStreamerInfo::kOffsetL + TStreamerInfo::kSTLp:
                     case TStreamerInfo::kOffsetL + TStreamerInfo::kAnyp:
                     case TStreamerInfo::kOffsetL + TStreamerInfo::kObjectp:
                     case TStreamerInfo::kOffsetL + TStreamerInfo::kAnyP:
                        pointer = kTRUE;
                        break;
                     case TStreamerInfo::kBase:
                     case TStreamerInfo::kAny :
                     case TStreamerInfo::kSTL:
                     case TStreamerInfo::kObject:
                     case TStreamerInfo::kTString:
                     case TStreamerInfo::kTNamed:
                     case TStreamerInfo::kTObject:
                        object = kTRUE;
                        break;
                     case TStreamerInfo::kOffsetL + TStreamerInfo::kAny:
                     case TStreamerInfo::kOffsetL + TStreamerInfo::kSTL:
                     case TStreamerInfo::kOffsetL + TStreamerInfo::kObject:
                        objarr = kTRUE;
                        break;
                     case TStreamerInfo::kStreamer:
                     case TStreamerInfo::kStreamLoop:
                        // Unsupported case.
                        Error("DefinedVariable",
                              "%s is a datamember of %s BUT is not yet of a supported type (%d)",
                              right,cl ? cl->GetName() : "unknown class",type);
                        return -2;
                     default:
                        // Unknown and Unsupported case.
                        Error("DefinedVariable",
                              "%s is a datamember of %s BUT is not of a unknown type (%d)",
                              right,cl ? cl->GetName() : "unknown class",type);
                        return -2;
                  }

                  if (object && !useCollectionObject &&
                      ( element->GetClassPointer() ==  TClonesArray::Class()
                        || element->GetClassPointer()->GetCollectionProxy() ) )
                  {
                     object = kFALSE;
                  }
                  if (object && leafinfo) {
                     leafinfo->AddOffset(offset,element);
                  } else if (objarr) {
                     // This is an embedded array of objects. We can not increase the offset.
                     leafinfo = new TFormLeafInfo(cl,offset,element);
                     mustderef = kTRUE;
                  } else {

                     if (!useCollectionObject && element->GetClassPointer() ==  TClonesArray::Class()) {

                        leafinfo = new TFormLeafInfoClones(cl,offset,element);
                        mustderef = kTRUE;

                     } else if (!useCollectionObject &&  element->GetClassPointer()
                                && element->GetClassPointer()->GetCollectionProxy()) {

                        mustderef = kTRUE;
                        if (numberOfVarDim>1) {
                           Warning("DefinedVariable","TTreeFormula support only 2 level of variables size collections.  Assuming '@' notation for the collection %s.",
                                   element->GetName());
                           leafinfo = new TFormLeafInfo(cl,offset,element);
                           useCollectionObject = kTRUE;
                        } else if (numberOfVarDim==1) {
                           R__ASSERT(maininfo);
                           leafinfo =
                              new TFormLeafInfoMultiVarDimCollection(cl,offset,element,maininfo);

                           fHasMultipleVarDim[code] = kTRUE;
                           //numberOfVarDim += RegisterDimensions(code,leafinfo);
                           //cl = cl->GetCollectionProxy()->GetValueClass();

                           //if (maininfo==0) maininfo = leafinfo;
                           //if (previnfo==0) previnfo = leafinfo;
                           //else {
                           //   previnfo->fNext = leafinfo;
                           //   previnfo = leafinfo;
                           //}
                           leafinfo->fNext =  new TFormLeafInfoCollection(cl, offset, element);
                           if (element->GetClassPointer()->GetCollectionProxy()->GetValueClass()==0) {
                              TFormLeafInfo *info = new TFormLeafInfoNumerical(
                                 element->GetClassPointer()->GetCollectionProxy());
                              if (leafinfo->fNext) leafinfo->fNext->fNext = info;
                              else leafinfo->fNext = info;
                           }
                        } else {
                           leafinfo =  new TFormLeafInfoCollection(cl, offset, element);

                           TClass *elemCl = element->GetClassPointer();
                           TClass *valueCl = elemCl->GetCollectionProxy()->GetValueClass();
                           if (!maininfo) maininfo = leafinfo;

                           if (valueCl!=0 && valueCl->GetCollectionProxy()!=0) {

                              numberOfVarDim += RegisterDimensions(code,leafinfo,maininfo,kFALSE);
                              if (previnfo==0) previnfo = leafinfo;
                              else {
                                 previnfo->fNext = leafinfo;
                                 previnfo = leafinfo;
                              }
                              leafinfo = new TFormLeafInfoMultiVarDimCollection(elemCl,0,
                                 elemCl->GetCollectionProxy()->GetValueClass(),maininfo);
                              //numberOfVarDim += RegisterDimensions(code,previnfo->fNext);
                              fHasMultipleVarDim[code] = kTRUE;
                              //previnfo = previnfo->fNext;
                              leafinfo->fNext =  new TFormLeafInfoCollection(elemCl,0,
                                 valueCl);
                              elemCl = valueCl;
                           }
                           if (elemCl->GetCollectionProxy() &&
                               elemCl->GetCollectionProxy()->GetValueClass()==0) {
                              TFormLeafInfo *info = new TFormLeafInfoNumerical(elemCl->GetCollectionProxy());
                              if (leafinfo->fNext) leafinfo->fNext->fNext = info;
                              else leafinfo->fNext = info;
                           }
                        }
                     } else if ( (object || pointer) && !useReferenceObject && element->GetClassPointer()->GetReferenceProxy() ) {
                        TClass* c = element->GetClassPointer();
                        R__LoadBranch(leaf->GetBranch(),readentry,fQuickLoad);
                        if ( object )  {
                           leafinfo = new TFormLeafInfoReference(c, element, offset);
                        }
                        else  {
                           leafinfo = new TFormLeafInfoPointer(cl,offset,element);
                           leafinfo->fNext = new TFormLeafInfoReference(c, element, 0);
                        }
                        //if ( c->GetReferenceProxy()->HasCounter() )  {
                        //   numberOfVarDim += RegisterDimensions(code,-1);
                        //}
                        prevUseReferenceObject = kFALSE;
                        needClass = kFALSE;
                        mustderef = kTRUE;
                     } else if (pointer) {
                        // this is a pointer to be followed.
                        leafinfo = new TFormLeafInfoPointer(cl,offset,element);
                        mustderef = kTRUE;
                     } else {
                        // this is an embedded object.
                        R__ASSERT(object);
                        leafinfo = new TFormLeafInfo(cl,offset,element);
                     }
                  }
               }
            } else {
               if (cl) Error("DefinedVariable","%s is not a datamember of %s",work,cl->GetName());
               // no else, we warned earlier that the class was missing.
               return -1;
            }

            numberOfVarDim += RegisterDimensions(code,leafinfo,maininfo,useCollectionObject); // Note or useCollectionObject||prevUseColectionObject
            if (maininfo==0) {
               maininfo = leafinfo;
            }
            if (previnfo==0) {
               previnfo = leafinfo;
            } else if (previnfo!=leafinfo) {
               previnfo->fNext = leafinfo;
               previnfo = leafinfo;
            }
            while (previnfo->fNext) previnfo = previnfo->fNext;

            if ( right[i] != '\0' )  {
               if ( !needClass && mustderef )   {
                  maininfo->SetBranch(leaf->GetBranch());
                  char *ptr = (char*)maininfo->GetValuePointer(leaf,0);
                  TFormLeafInfoReference* refInfo = 0;
                  if ( !maininfo->IsReference() )  {
                     for( TFormLeafInfo* inf = maininfo->fNext; inf; inf = inf->fNext )  {
                        if ( inf->IsReference() )  {
                           refInfo = (TFormLeafInfoReference*)inf;
                        }
                     }
                  }
                  else {
                     refInfo = (TFormLeafInfoReference*)maininfo;
                  }
                  if ( refInfo )  {
                     cl = refInfo->GetValueClass(ptr);
                     if ( !cl )  {
                        Error("DefinedVariable","Failed to access class type of reference target (%s)",element->GetName());
                        return -1;
                     }
                     element = ((TStreamerInfo*)cl->GetStreamerInfo())->GetStreamerElement(work,offset);
                  }
                  else  {
                     Error("DefinedVariable","Failed to access class type of reference target (%s)",element->GetName());
                     return -1;
                  }
               }
               else if ( needClass )  {
                  cl = element->GetClassPointer();
               }
            }
            if (mustderef) leafinfo = 0;
            current = &(work[0]);
            *current = 0;
            R__ASSERT(right[i] != '[');  // We are supposed to have removed all dimensions already!

            if (cl == TString::Class() && strcmp(right+i+1,"fData") == 0) {
               // For backward compatibility replace TString::fData which no longer exist
               // by a call to TString::Data()
               right = ".Data()";
               i = 0;
               nchname = strlen(right);
            }

         } else
            *current++ = right[i];
      }
      if (maininfo) {
         fDataMembers.AddAtAndExpand(maininfo,code);
         if (leaf) fLookupType[code] = kDataMember;
         else fLookupType[code] = kTreeMember;
      }
   }

   if (strlen(work)!=0) {
      // We have something left to analyze.  Let's make this an error case!
      return -1;
   }

   TClass *objClass = EvalClass(code);
   if (objClass && !useLeafCollectionObject && objClass->GetCollectionProxy() && objClass->GetCollectionProxy()->GetValueClass()) {
      TFormLeafInfo *last = 0;
      if ( SwitchToFormLeafInfo(code) ) {

         last = (TFormLeafInfo*)fDataMembers.At(code);

         if (!last) return action;
         while (last->fNext) { last = last->fNext; }

      }
      if (last && last->GetClass() != objClass) {
         TClass *mother_cl;
         if (leaf->IsA()==TLeafObject::Class()) {
            // in this case mother_cl is not really used
            mother_cl = cl;
         } else {
            mother_cl = ((TBranchElement*)branch)->GetInfo()->GetClass();
         }

         TFormLeafInfo* collectioninfo = new TFormLeafInfoCollection(mother_cl, 0, objClass, kFALSE);
         // The dimension needs to be handled!
         numberOfVarDim += RegisterDimensions(code,collectioninfo,maininfo,kFALSE);
         last->fNext = collectioninfo;
      }
      numberOfVarDim += RegisterDimensions(code,1); //NOTE: changed from 0
      objClass = objClass->GetCollectionProxy()->GetValueClass();
   }
   if (IsLeafString(code) || objClass == TString::Class() || objClass == stdStringClass) {

      TFormLeafInfo *last = 0;
      if ( SwitchToFormLeafInfo(code) ) {

         last = (TFormLeafInfo*)fDataMembers.At(code);

         if (!last) return action;
         while (last->fNext) { last = last->fNext; }

      }
      const char *funcname = 0;
      if (objClass == TString::Class()) {
         funcname = "Data";
         //tobetested: numberOfVarDim += RegisterDimensions(code,1,0); // Register the dim of the implied char*
      } else if (objClass == stdStringClass) {
         funcname = "c_str";
         //tobetested: numberOfVarDim += RegisterDimensions(code,1,0); // Register the dim of the implied char*
      }
      if (funcname) {
         TMethodCall *method = new TMethodCall(objClass, funcname, "");
         if (last) {
            last->fNext = new TFormLeafInfoMethod(objClass,method);
         } else {
            fDataMembers.AddAtAndExpand(new TFormLeafInfoMethod(objClass,method),code);
            if (leaf) fLookupType[code] = kDataMember;
            else fLookupType[code] = kTreeMember;
         }
      }
      return kDefinedString;
   }

   if (objClass) {
      TMethodCall *method = new TMethodCall(objClass, "AsDouble", "");
      if (method->IsValid()
          && (method->ReturnType() == TMethodCall::kLong || method->ReturnType() == TMethodCall::kDouble)) {

         TFormLeafInfo *last = 0;
         if (SwitchToFormLeafInfo(code)) {
            last = (TFormLeafInfo*)fDataMembers.At(code);
            // Improbable case
            if (!last) {
               delete method;
               return action;
            }
            while (last->fNext) { last = last->fNext; }
         }
         if (last) {
            last->fNext = new TFormLeafInfoMethod(objClass,method);
         } else {
            fDataMembers.AddAtAndExpand(new TFormLeafInfoMethod(objClass,method),code);
            if (leaf) fLookupType[code] = kDataMember;
            else fLookupType[code] = kTreeMember;
         }

         return kDefinedVariable;
      }
      delete method;
      method = new TMethodCall(objClass, "AsString", "");
      if (method->IsValid()
          && method->ReturnType() == TMethodCall::kString) {

         TFormLeafInfo *last = 0;
         if (SwitchToFormLeafInfo(code)) {
            last = (TFormLeafInfo*)fDataMembers.At(code);
            // Improbable case
            if (!last) {
               delete method;
               return action;
            }
            while (last->fNext) { last = last->fNext; }
         }
         if (last) {
            last->fNext = new TFormLeafInfoMethod(objClass,method);
         } else {
            fDataMembers.AddAtAndExpand(new TFormLeafInfoMethod(objClass,method),code);
            if (leaf) fLookupType[code] = kDataMember;
            else fLookupType[code] = kTreeMember;
         }

         //tobetested: numberOfVarDim += RegisterDimensions(code,1,0); // Register the dim of the implied char*
         return kDefinedString;
      }
      if (method->IsValid()
          && method->ReturnType() == TMethodCall::kOther) {

         TClass *rcl = TFormLeafInfoMethod::ReturnTClass(method);
         if ((rcl == TString::Class() || rcl == stdStringClass) ) {

            TFormLeafInfo *last = 0;
            if (SwitchToFormLeafInfo(code)) {
               last = (TFormLeafInfo*)fDataMembers.At(code);
               // Improbable case
               if (!last) {
                  delete method;
                  return action;
               }
               while (last->fNext) { last = last->fNext; }
            }
            if (last) {
               last->fNext = new TFormLeafInfoMethod(objClass,method);
               last = last->fNext;
            } else {
               last = new TFormLeafInfoMethod(objClass,method);
               fDataMembers.AddAtAndExpand(last,code);
               if (leaf) fLookupType[code] = kDataMember;
               else fLookupType[code] = kTreeMember;
            }

            objClass = rcl;

            const char *funcname = 0;
            if (objClass == TString::Class()) {
               funcname = "Data";
            } else if (objClass == stdStringClass) {
               funcname = "c_str";
            }
            if (funcname) {
               method = new TMethodCall(objClass, funcname, "");
               last->fNext = new TFormLeafInfoMethod(objClass,method);
            }
            return kDefinedString;
         }
      }
      delete method;
   }

   return action;
}

////////////////////////////////////////////////////////////////////////////////
/// Look for the leaf corresponding to the start of expression.
/// It returns the corresponding leaf if any.
/// It also modify the following arguments:
///
/// -  leftover: contain from expression that was not used to determine the leaf
/// -  final:
///    * paran_level: number of un-matched open parenthesis
///    * cast_queue: list of cast to be done
///    * aliases: list of aliases used
/// -  Return <0 in case of failure
///
/// -  Return 0 if a leaf has been found
/// -  Return 2 if info about the TTree itself has been requested.

Int_t TTreeFormula::FindLeafForExpression(const char* expression, TLeaf*& leaf, TString& leftover, Bool_t& final, UInt_t& paran_level, TObjArray& castqueue, std::vector<std::string>& aliasUsed, Bool_t& useLeafCollectionObject, const char* fullExpression)
{
   // Later on we will need to read one entry, let's make sure
   // it is a real entry.
   if (fTree->GetTree()==0) {
      fTree->LoadTree(0);
      if (fTree->GetTree()==0) return -1;
   }
   Long64_t readentry = fTree->GetTree()->GetReadEntry();
   if (readentry < 0) readentry=0;
   const char *cname = expression;
   char    first[kMaxLen];    first[0] = '\0';
   char   second[kMaxLen*2]; second[0] = '\0';
   char    right[kMaxLen*2];  right[0] = '\0';
   char     work[kMaxLen];     work[0] = '\0';
   char     left[kMaxLen];     left[0] = '\0';
   char  scratch[kMaxLen*5];
   char scratch2[kMaxLen*5];
   std::string currentname;
   Int_t previousdot = 0;
   char *current;
   TLeaf *tmp_leaf=0;
   TBranch *branch=0, *tmp_branch=0;
   Int_t nchname = strlen(cname);
   Int_t i;
   Bool_t foundAtSign = kFALSE;
   Bool_t startWithParan = kFALSE;

   for (i=0, current = &(work[0]); i<=nchname && !final;i++ ) {
      // We will treated the terminator as a token.
      *current++ = cname[i];

      if (cname[i] == '(') {
         ++paran_level;

         if (current==work+1) {
            // If the expression starts with a parenthesis, we are likely
            // to have a cast operator inside.
            startWithParan = kTRUE;
            current--;
         }
         continue;
         //i++;
         //while( cname[i]!=')' && cname[i] ) {
         //   *current++ = cname[i++];
         //}
         //*current++ = cname[i];
         ////*current = 0;
         //continue;
      }
      if (cname[i] == ')') {
         if (paran_level==0) {
            Error("DefinedVariable","Unmatched parenthesis in %s",fullExpression);
            return -1;
         }
         paran_level--;

         if (startWithParan) {
            startWithParan = kFALSE; // the next match wont be against the starting parenthesis.

            // Let's see if work is a classname and thus we have a cast.
            *(--current) = 0;
            TString cast_name = gInterpreter->TypeName(work);
            TClass *cast_cl = TClass::GetClass(cast_name);
            if (cast_cl) {
               // We must have a cast
               castqueue.AddAtAndExpand(cast_cl,paran_level);
               current = &(work[0]);
               *current = 0;
               //            Warning("DefinedVariable","Found cast to %s",cast_fullExpression);
               continue;
            } else if (gROOT->GetType(cast_name)) {
               // We reset work
               current = &(work[0]);
               *current = 0;
               Warning("DefinedVariable",
                       "Casting to primary types like \"%s\" is not supported yet",cast_name.Data());
               continue;
            }
            *(current++)=')';
         }

         *current='\0';
         char *params = strchr(work,'(');
         if (params) {
            *params = 0; params++;

            if (branch && !leaf) {
               // We have a branch but not a leaf.  We are likely to have found
               // the top of split branch.
               if (BranchHasMethod(0, branch, work, params, readentry)) {
                  //fprintf(stderr, "Does have a method %s for %s.\n", work, branch->GetName());
               }
            }

            // What we have so far might be a member function of one of the
            // leaves that are not split (for example "GetNtrack" for the Event class).
            TIter next(fTree->GetIteratorOnAllLeaves());
            TLeaf* leafcur = 0;
            while (!leaf && (leafcur = (TLeaf*) next())) {
               TBranch* br = leafcur->GetBranch();
               Bool_t yes = BranchHasMethod(leafcur, br, work, params, readentry);
               if (yes) {
                  leaf = leafcur;
                  //fprintf(stderr, "Does have a method %s for %s found in leafcur %s.\n", work, leafcur->GetBranch()->GetName(), leafcur->GetName());
               }
            }
            if (!leaf) {
               // Check for an alias.
               if (strlen(left) && left[strlen(left)-1]=='.') left[strlen(left)-1]=0;
               const char *aliasValue = fTree->GetAlias(left);
               if (aliasValue && strcspn(aliasValue,"+*/-%&!=<>|")==strlen(aliasValue)) {
                  // First check whether we are using this alias recursively (this would
                  // lead to an infinite recursion.
                  if (find(aliasUsed.begin(),
                     aliasUsed.end(),
                     left) != aliasUsed.end()) {
                        Error("DefinedVariable",
                           "The substitution of the branch alias \"%s\" by \"%s\" in \"%s\" failed\n"\
                           "\tbecause \"%s\" is used [recursively] in its own definition!",
                           left,aliasValue,fullExpression,left);
                        return -3;
                  }
                  aliasUsed.push_back(left);
                  TString newExpression = aliasValue;
                  newExpression += (cname+strlen(left));
                  Int_t res = FindLeafForExpression(newExpression, leaf, leftover, final, paran_level,
                     castqueue, aliasUsed, useLeafCollectionObject, fullExpression);
                  if (res<0) {
                     Error("DefinedVariable",
                        "The substitution of the alias \"%s\" by \"%s\" failed.",left,aliasValue);
                     return -3;
                  }
                  return res;
               }

               // This is actually not really any error, we probably received something
               // like "abs(some_val)", let ROOT::v5::TFormula decompose it first.
               return -1;
            }
            //         if (!leaf->InheritsFrom(TLeafObject::Class()) ) {
            // If the leaf that we found so far is not a TLeafObject then there is
            // nothing we would be able to do.
            //   Error("DefinedVariable","Need a TLeafObject to call a function!");
            // return -1;
            //}
            // We need to recover the info not used.
            strlcpy(right,work,2*kMaxLen);
            strncat(right,"(",2*kMaxLen-1-strlen(right));
            strncat(right,params,2*kMaxLen-1-strlen(right));
            final = kTRUE;

            // Record in 'i' what we consumed
            i += strlen(params);

            // we reset work
            current = &(work[0]);
            *current = 0;
            break;
         }
      }
      if (cname[i] == '.' || cname[i] == '\0' || cname[i] == ')') {
         // A delimiter happened let's see if what we have seen
         // so far does point to a leaf.
         *current = '\0';

         Int_t len = strlen(work);
         if (work[0]=='@') {
            foundAtSign = kTRUE;
            Int_t l = 0;
            for(l=0;work[l+1]!=0;++l) work[l] = work[l+1];
            work[l] = '\0';
            --current;
         } else if (len>=2 && work[len-2]=='@') {
            foundAtSign = kTRUE;
            work[len-2] = cname[i];
            work[len-1] = '\0';
            --current;
         } else {
            foundAtSign = kFALSE;
         }

         if (left[0]==0) strlcpy(left,work,kMaxLen);
         if (!leaf && !branch) {
            // So far, we have not found a matching leaf or branch.
            strlcpy(first,work,kMaxLen);

            std::string treename(first);
            if (treename.size() && treename[treename.size()-1]=='.') {
               treename.erase(treename.size()-1);
            }
            if (treename== "This" /* || treename == fTree->GetName() */ ) {
               // Request info about the TTree object itself,
               TNamed *named = new TNamed(fTree->GetName(),fTree->GetName());
               fLeafNames.AddAtAndExpand(named,fNcodes);
               fLeaves.AddAtAndExpand(fTree,fNcodes);
               if (cname[i]) leftover = &(cname[i+1]);
               return 2;
            }
            // The following would allow to access the friend by name
            // however, it would also prevent the access of the leaves
            // within the friend.  We could use the '@' notation here
            // however this would not be aesthetically pleasing :(
            // What we need to do, is add the ability to look ahead to
            // the next 'token' to decide whether we to access the tree
            // or its leaf.
            //} else {
            //   TTree *tfriend = fTree->GetFriend(treename.c_str());
            //   TTree *realtree = fTree->GetTree();
            //   if (!tfriend && realtree != fTree){
            //      // If it is a chain and we did not find a friend,
            //      // let's try with the internal tree.
            //      tfriend = realtree->GetFriend(treename.c_str());
            //   }
            //   if (tfriend) {
            //      TNamed *named = new TNamed(treename.c_str(),tfriend->GetName());
            //      fLeafNames.AddAtAndExpand(named,fNcodes);
            //      fLeaves.AddAtAndExpand(tfriend,fNcodes);
            //      if (cname[i]) leftover = &(cname[i+1]);
            //      return 2;
            //   }
            //}

            branch = fTree->FindBranch(first);
            leaf = fTree->FindLeaf(first);

            // Now look with the delimiter removed (we looked with it first
            // because a dot is allowed at the end of some branches).
            if (cname[i]) first[strlen(first)-1]='\0';
            if (!branch) branch = fTree->FindBranch(first);
            if (!leaf) leaf = fTree->FindLeaf(first);
            TClass* cl = 0;
            if ( branch && branch->InheritsFrom(TBranchElement::Class()) ) {
               int offset=0;
               TBranchElement* bElt = (TBranchElement*)branch;
               TStreamerInfo* info  = bElt->GetInfo();
               TStreamerElement* element = info ? info->GetStreamerElement(first,offset) : 0;
               if (element) cl = element->GetClassPointer();
               if ( cl && !cl->GetReferenceProxy() ) cl = 0;
            }
            if ( cl )  {  // We have a reference class here....
               final = kTRUE;
               useLeafCollectionObject = foundAtSign;
               // we reset work
               current = &(work[0]);
               *current = 0;
            }
            else if (branch && (foundAtSign || cname[i] != 0)  ) {
               // Since we found a branch and there is more information in the name,
               // we do NOT look at the 'IsOnTerminalBranch' status of the leaf
               // we found ... yet!

               if (leaf==0) {
                  // Note we do not know (yet?) what (if anything) to do
                  // for a TBranchObject branch.
                  if (branch->InheritsFrom(TBranchElement::Class()) ) {
                     Int_t type = ((TBranchElement*)branch)->GetType();
                     if ( type == 3 || type ==4) {
                        // We have a Collection branch.
                        leaf = (TLeaf*)branch->GetListOfLeaves()->At(0);
                        if (foundAtSign) {
                           useLeafCollectionObject = foundAtSign;
                           foundAtSign = kFALSE;
                           current = &(work[0]);
                           *current = 0;
                           ++i;
                           break;
                        }
                     }
                  }
               }

               // we reset work
               useLeafCollectionObject = foundAtSign;
               foundAtSign = kFALSE;
               current = &(work[0]);
               *current = 0;
            } else if (leaf || branch) {
               if (leaf && branch) {
                  // We found both a leaf and branch matching the request name
                  // let's see which one is the proper one to use! (On annoying case
                  // is that where the same name is repeated ( varname.varname )

                  // We always give priority to the branch
                  // leaf = 0;
               }
               if (leaf && leaf->IsOnTerminalBranch()) {
                  // This is a non-object leaf, it should NOT be specified more except for
                  // dimensions.
                  final = kTRUE;
               }
               // we reset work
               current = &(work[0]);
               *current = 0;
            } else {
               // What we have so far might be a data member of one of the
               // leaves that are not split (for example "fNtrack" for the Event class.
               TLeaf *leafcur = GetLeafWithDatamember(first,work,readentry);
               if (leafcur) {
                  leaf = leafcur;
                  branch = leaf->GetBranch();
                  if (leaf->IsOnTerminalBranch()) {
                     final = kTRUE;
                     strlcpy(right,first,kMaxLen);
                     //We need to put the delimiter back!
                     if (foundAtSign) strncat(right,"@",2*kMaxLen-1-strlen(right));
                     if (cname[i]=='.') strncat(right,".",2*kMaxLen-1-strlen(right));

                     // We reset work
                     current = &(work[0]);
                     *current = 0;
                  };
               } else if (cname[i] == '.') {
                  // If we have a branch that match a name preceded by a dot
                  // then we assume we are trying to drill down the branch
                  // Let look if one of the top level branch has a branch with the name
                  // we are looking for.
                  TBranch *branchcur;
                  TIter next( fTree->GetListOfBranches() );
                  while(!branch && (branchcur=(TBranch*)next()) ) {
                     branch = branchcur->FindBranch(first);
                  }
                  if (branch) {
                     // We reset work
                     current = &(work[0]);
                     *current = 0;
                  }
               }
            }
         } else {  // correspond to if (leaf || branch)
            if (final) {
               Error("DefinedVariable", "Unexpected control flow!");
               return -1;
            }

            // No dot is allowed in subbranches and leaves, so
            // we always remove it in the present case.
            if (cname[i]) work[strlen(work)-1] = '\0';
            snprintf(scratch,sizeof(scratch),"%s.%s",first,work);
            snprintf(scratch2,sizeof(scratch2),"%s.%s.%s",first,second,work);

            if (previousdot) {
               currentname = &(work[previousdot+1]);
            }

            // First look for the current 'word' in the list of
            // leaf of the
            if (branch) {
               tmp_leaf = branch->FindLeaf(work);
               if (!tmp_leaf)  tmp_leaf = branch->FindLeaf(scratch);
               if (!tmp_leaf)  tmp_leaf = branch->FindLeaf(scratch2);
               if (!tmp_leaf)  tmp_leaf = branch->FindLeaf(currentname.c_str());
            }
            if (tmp_leaf && tmp_leaf->IsOnTerminalBranch() ) {
               // This is a non-object leaf, it should NOT be specified more except for
               // dimensions.
               final = kTRUE;
            }

            if (branch) {
               tmp_branch = branch->FindBranch(work);
               if (!tmp_branch) tmp_branch = branch->FindBranch(scratch);
               if (!tmp_branch) tmp_branch = branch->FindBranch(scratch2);
               if (!tmp_branch) tmp_branch = branch->FindBranch(currentname.c_str());
            }
            if (tmp_branch) {
               branch=tmp_branch;

               // NOTE: Should we look for a leaf within here?
               if (!final) {
                  tmp_leaf = branch->FindLeaf(work);
                  if (!tmp_leaf)  tmp_leaf = branch->FindLeaf(scratch);
                  if (!tmp_leaf)  tmp_leaf = branch->FindLeaf(scratch2);
                  if (!tmp_leaf)  tmp_leaf = branch->FindLeaf(currentname.c_str());
                  if (tmp_leaf && tmp_leaf->IsOnTerminalBranch() ) {
                     // This is a non-object leaf, it should NOT be specified
                     // more except for dimensions.
                     final = kTRUE;
                     leaf = tmp_leaf;
                  }
               }
            }
            if (tmp_leaf) {
               // Something was found.
               if (second[0]) strncat(second,".",2*kMaxLen-1-strlen(second));
               strncat(second,work,2*kMaxLen-1-strlen(second));
               leaf = tmp_leaf;
               useLeafCollectionObject = foundAtSign;
               foundAtSign = kFALSE;

               // we reset work
               current = &(work[0]);
               *current = 0;
            } else {
               //We need to put the delimiter back!
               if (strlen(work)) {
                  if (foundAtSign) {
                     Int_t where = strlen(work);
                     work[where] = '@';
                     work[where+1] = cname[i];
                     ++current;
                     previousdot = where+1;
                  } else {
                     previousdot = strlen(work);
                     work[strlen(work)] = cname[i];
                  }
               } else --current;
            }
         }
      }
   }

   // Copy the left over for later use.
   if (strlen(work)) {
      strncat(right,work,2*kMaxLen-1-strlen(right));
   }

   if (i<nchname) {
      if (strlen(right) && right[strlen(right)-1]!='.' && cname[i]!='.') {
         // In some cases we remove a little to fast the period, we add
         // it back if we need.  It is assumed that 'right' and the rest of
         // the name was cut by a delimiter, so this should be safe.
         strncat(right,".",2*kMaxLen-1-strlen(right));
      }
      strncat(right,&cname[i],2*kMaxLen-1-strlen(right));
   }

   if (!final && branch) {
      if (!leaf) {
         leaf = (TLeaf*)branch->GetListOfLeaves()->UncheckedAt(0);
         if (!leaf) return -1;
      }
      final = leaf->IsOnTerminalBranch();
   }

   if (leaf && leaf->InheritsFrom(TLeafObject::Class()) ) {
      if (strlen(right)==0) strlcpy(right,work,2*kMaxLen);
   }

   if (leaf==0 && left[0]!=0) {
      if (left[strlen(left)-1]=='.') left[strlen(left)-1]=0;

      // Check for an alias.
      const char *aliasValue = fTree->GetAlias(left);
      if (aliasValue && strcspn(aliasValue,"()[]+*/-%&!=<>|")==strlen(aliasValue)) {
         // First check whether we are using this alias recursively (this would
         // lead to an infinite recursion).
         if (find(aliasUsed.begin(),
                  aliasUsed.end(),
                  left) != aliasUsed.end()) {
            Error("DefinedVariable",
                  "The substitution of the branch alias \"%s\" by \"%s\" in \"%s\" failed\n"\
                  "\tbecause \"%s\" is used [recursively] in its own definition!",
                  left,aliasValue,fullExpression,left);
            return -3;
         }
         aliasUsed.push_back(left);
         TString newExpression = aliasValue;
         newExpression += (cname+strlen(left));
         Int_t res = FindLeafForExpression(newExpression, leaf, leftover, final, paran_level,
                                           castqueue, aliasUsed, useLeafCollectionObject, fullExpression);
         if (res<0) {
            Error("DefinedVariable",
                  "The substitution of the alias \"%s\" by \"%s\" failed.",left,aliasValue);
            return -3;
         }
         return res;
      }
   }
   leftover = right;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if name is in the list of Tree/Branch leaves.
///
/// This member function redefines the function in ROOT::v5::TFormula
/// If a leaf has a name corresponding to the argument name, then
/// returns a new code.
///
/// A TTreeFormula may contain more than one variable.
/// For each variable referenced, the pointers to the corresponding
/// branch and leaf is stored in the object arrays fBranches and fLeaves.
///
/// name can be :
///  - Leaf_Name (simple variable or data member of a ClonesArray)
///  - Branch_Name.Leaf_Name
///  - Branch_Name.Method_Name
///  - Leaf_Name[index]
///  - Branch_Name.Leaf_Name[index]
///  - Branch_Name.Leaf_Name[index1]
///  - Branch_Name.Leaf_Name[][index2]
///  - Branch_Name.Leaf_Name[index1][index2]
///
/// New additions:
///  - Branch_Name.Leaf_Name[OtherLeaf_Name]
///  - Branch_Name.Datamember_Name
///  - '.' can be replaced by '->'
///
///   and
///  - Branch_Name[index1].Leaf_Name[index2]
///  - Leaf_name[index].Action().OtherAction(param)
///  - Leaf_name[index].Action()[val].OtherAction(param)
///
/// The expected returns values are
/// -  -2 :  the name has been recognized but won't be usable
/// -  -1 :  the name has not been recognized
/// -  >=0 :  the name has been recognized, return the internal code for this name.

Int_t TTreeFormula::DefinedVariable(TString &name, Int_t &action)
{

   action = kDefinedVariable;
   if (!fTree) return -1;

   fNpar = 0;
   if (name.Length() > kMaxLen) return -1;
   Int_t i,k;

   if (name == "Entry$") {
      Int_t code = fNcodes++;
      fCodes[code] = 0;
      fLookupType[code] = kIndexOfEntry;
      return code;
   }
   if (name == "LocalEntry$") {
      Int_t code = fNcodes++;
      fCodes[code] = 0;
      fLookupType[code] = kIndexOfLocalEntry;
      return code;
   }
   if (name == "Entries$") {
      Int_t code = fNcodes++;
      fCodes[code] = 0;
      fLookupType[code] = kEntries;
      SetBit(kNeedEntries);
      fManager->SetBit(kNeedEntries);
      return code;
   }
   if (name == "LocalEntries$") {
      Int_t code = fNcodes++;
      fCodes[code] = 0;
      fLookupType[code] = kLocalEntries;
      SetBit(kNeedEntries); // FIXME: necessary?
      fManager->SetBit(kNeedEntries); // FIXME: necessary?
      return code;
   }
   if (name == "Iteration$") {
      Int_t code = fNcodes++;
      fCodes[code] = 0;
      fLookupType[code] = kIteration;
      return code;
   }
   if (name == "Length$") {
      Int_t code = fNcodes++;
      fCodes[code] = 0;
      fLookupType[code] = kLength;
      return code;
   }
   static const char *lenfunc = "Length$(";
   if (strncmp(name.Data(),"Length$(",strlen(lenfunc))==0
       && name[name.Length()-1]==')') {

      TString subform = name.Data()+strlen(lenfunc);
      subform.Remove( subform.Length() - 1 );
      TTreeFormula *lengthForm = new TTreeFormula("lengthForm",subform,fTree);
      fAliases.AddAtAndExpand(lengthForm,fNoper);
      Int_t code = fNcodes++;
      fCodes[code] = 0;
      fLookupType[code] = kLengthFunc;
      return code;
   }
   static const char *minfunc = "Min$(";
   if (strncmp(name.Data(),"Min$(",strlen(minfunc))==0
       && name[name.Length()-1]==')') {

      TString subform = name.Data()+strlen(minfunc);
      subform.Remove( subform.Length() - 1 );
      TTreeFormula *minForm = new TTreeFormula("minForm",subform,fTree);
      fAliases.AddAtAndExpand(minForm,fNoper);
      Int_t code = fNcodes++;
      fCodes[code] = 0;
      fLookupType[code] = kMin;
      return code;
   }
   static const char *maxfunc = "Max$(";
   if (strncmp(name.Data(),"Max$(",strlen(maxfunc))==0
       && name[name.Length()-1]==')') {

      TString subform = name.Data()+strlen(maxfunc);
      subform.Remove( subform.Length() - 1 );
      TTreeFormula *maxForm = new TTreeFormula("maxForm",subform,fTree);
      fAliases.AddAtAndExpand(maxForm,fNoper);
      Int_t code = fNcodes++;
      fCodes[code] = 0;
      fLookupType[code] = kMax;
      return code;
   }
   static const char *sumfunc = "Sum$(";
   if (strncmp(name.Data(),"Sum$(",strlen(sumfunc))==0
       && name[name.Length()-1]==')') {

      TString subform = name.Data()+strlen(sumfunc);
      subform.Remove( subform.Length() - 1 );
      TTreeFormula *sumForm = new TTreeFormula("sumForm",subform,fTree);
      fAliases.AddAtAndExpand(sumForm,fNoper);
      Int_t code = fNcodes++;
      fCodes[code] = 0;
      fLookupType[code] = kSum;
      return code;
   }



   // Check for $Alt(expression1,expression2)
   Int_t res = DefineAlternate(name.Data());
   if (res!=0) {
      // There was either a syntax error or we found $Alt
      if (res<0) return res;
      action = res;
      return 0;
   }

   // Find the top level leaf and deal with dimensions

   char    cname[kMaxLen];  strlcpy(cname,name.Data(),kMaxLen);
   char     dims[kMaxLen];   dims[0] = '\0';

   Bool_t final = kFALSE;

   UInt_t paran_level = 0;
   TObjArray castqueue;

   // First, it is easier to remove all dimensions information from 'cname'
   Int_t cnamelen = strlen(cname);
   for(i=0,k=0; i<cnamelen; ++i, ++k) {
      if (cname[i] == '[') {
         int bracket = i;
         int bracket_level = 1;
         int j;
         for (j=++i; j<cnamelen && (bracket_level>0 || cname[j]=='['); j++, i++) {
            if (cname[j]=='[') bracket_level++;
            else if (cname[j]==']') bracket_level--;
         }
         if (bracket_level != 0) {
            //Error("DefinedVariable","Bracket unbalanced");
            return -1;
         }
         strncat(dims,&cname[bracket],j-bracket);
         //k += j-bracket;
      }
      if (i!=k) cname[k] = cname[i];
   }
   cname[k]='\0';

   Bool_t useLeafCollectionObject = kFALSE;
   TString leftover;
   TLeaf *leaf = 0;
   {
      std::vector<std::string> aliasSofar = fAliasesUsed;
      res = FindLeafForExpression(cname, leaf, leftover, final, paran_level, castqueue, aliasSofar, useLeafCollectionObject, name);
   }
   if (res<0) return res;

   if (!leaf && res!=2) {
      // Check for an alias.
      const char *aliasValue = fTree->GetAlias(cname);
      if (aliasValue) {
         // First check whether we are using this alias recursively (this would
         // lead to an infinite recursion.
         if (find(fAliasesUsed.begin(),
                  fAliasesUsed.end(),
                  cname) != fAliasesUsed.end()) {
            Error("DefinedVariable",
                  "The substitution of the alias \"%s\" by \"%s\" failed\n"\
                  "\tbecause \"%s\" is recursively used in its own definition!",
                  cname,aliasValue,cname);
            return -3;
         }


         if (strcspn(aliasValue,"()+*/-%&!=<>|")!=strlen(aliasValue)) {
            // If the alias contains an operator, we need to use a nested formula
            // (since DefinedVariable must only add one entry to the operation's list).

            // Need to check the aliases used so far
            std::vector<std::string> aliasSofar = fAliasesUsed;
            aliasSofar.push_back( cname );

            TString subValue( aliasValue );
            if (dims[0]) {
               subValue += dims;
            }

            TTreeFormula *subform = new TTreeFormula(cname,subValue,fTree,aliasSofar); // Need to pass the aliases used so far.

            if (subform->GetNdim()==0) {
               delete subform;
               Error("DefinedVariable",
                     "The substitution of the alias \"%s\" by \"%s\" failed.",cname,aliasValue);
               return -3;
            }

            fManager->Add(subform);
            fAliases.AddAtAndExpand(subform,fNoper);

            if (subform->IsString()) {
               action = kAliasString;
               return 0;
            } else {
               action = kAlias;
               return 0;
            }
         } else { /* assumes strcspn(aliasValue,"[]")!=strlen(aliasValue) */
            TString thisAlias( aliasValue );
            thisAlias += dims;
            Int_t aliasRes = DefinedVariable(thisAlias,action);
            if (aliasRes<0) {
               // We failed but DefinedVariable has not printed why yet.
               // and because we want those to be printed _before_ the notice
               // of the failure of the substitution, we need to print them here.
               if (aliasRes==-1) {
                  Error("Compile", " Bad numerical expression : \"%s\"",thisAlias.Data());
               } else if (aliasRes==-2) {
                  Error("Compile", " Part of the Variable \"%s\" exists but some of it is not accessible or useable",thisAlias.Data());

               }
               Error("DefinedVariable",
                     "The substitution of the alias \"%s\" by \"%s\" failed.",cname,aliasValue);
               return -3;
            }
            return aliasRes;
         }
      }
   }


   if (leaf || res==2) {

      if (leaf && leaf->GetBranch() && leaf->GetBranch()->TestBit(kDoNotProcess)) {
         Error("DefinedVariable","the branch \"%s\" has to be enabled to be used",leaf->GetBranch()->GetName());
         return -2;
      }

      Int_t code = fNcodes++;

      // If needed will now parse the indexes specified for
      // arrays.
      if (dims[0]) {
         char *current = &( dims[0] );
         Int_t dim = 0;
         TString varindex;
         Int_t index;
         Int_t scanindex ;
         while (current) {
            current++;
            if (current[0] == ']') {
               fIndexes[code][dim] = -1; // Loop over all elements;
            } else {
               scanindex = sscanf(current,"%d",&index);
               if (scanindex) {
                  fIndexes[code][dim] = index;
               } else {
                  fIndexes[code][dim] = -2; // Index is calculated via a variable.
                  varindex = current;
                  char *end = (char*)(varindex.Data());
                  for(char bracket_level = 0;*end!=0;end++) {
                     if (*end=='[') bracket_level++;
                     if (bracket_level==0 && *end==']') break;
                     if (*end==']') bracket_level--;
                  }
                  *end = '\0';
                  fVarIndexes[code][dim] = new TTreeFormula("index_var",
                                                            varindex,
                                                            fTree);
                  if (fVarIndexes[code][dim]->GetNdim() == 0) {
                     // Parsing failed for the index, let's stop here ....
                     return -1;
                  }
                  current += strlen(varindex)+1; // move to the end of the index array
               }
            }
            dim ++;
            if (dim >= kMAXFORMDIM) {
               // NOTE: test that dim this is NOT too big!!
               break;
            }
            current = (char*)strstr( current, "[" );
         }
      }

      // Now that we have cleaned-up the expression, let's compare it to the content
      // of the leaf!

      res = ParseWithLeaf(leaf,leftover,final,paran_level,castqueue,useLeafCollectionObject,name);
      if (res<0) return res;
      if (res>0) action = res;
      return code;
   }

//*-*- May be a graphical cut ?
   TCutG *gcut = (TCutG*)gROOT->GetListOfSpecials()->FindObject(name.Data());
   if (gcut) {
      if (gcut->GetObjectX()) {
         if(!gcut->GetObjectX()->InheritsFrom(TTreeFormula::Class()))
            gcut->SetObjectX(nullptr);
      }
      if (gcut->GetObjectY()) {
         if(!gcut->GetObjectY()->InheritsFrom(TTreeFormula::Class()))
            gcut->SetObjectY(nullptr);
      }

      Int_t code = fNcodes;

      if (strlen(gcut->GetVarX()) && strlen(gcut->GetVarY()) ) {

         TTreeFormula *fx = new TTreeFormula("f_x",gcut->GetVarX(),fTree);
         gcut->SetObjectX(fx);

         TTreeFormula *fy = new TTreeFormula("f_y",gcut->GetVarY(),fTree);
         gcut->SetObjectY(fy);

         fCodes[code] = -2;

      } else if (strlen(gcut->GetVarX())) {

         // Let's build the equivalent formula:
         // min(gcut->X) <= VarX <= max(gcut->Y)
         Double_t min = 0;
         Double_t max = 0;
         Int_t n = gcut->GetN();
         Double_t *x = gcut->GetX();
         min = max = x[0];
         for(Int_t i2 = 1; i2<n; i2++) {
            if (x[i2] < min) min = x[i2];
            if (x[i2] > max) max = x[i2];
         }
         TString formula = "(";
         formula += min;
         formula += "<=";
         formula += gcut->GetVarX();
         formula += " && ";
         formula += gcut->GetVarX();
         formula += "<=";
         formula += max;
         formula += ")";

         TTreeFormula *fx = new TTreeFormula("f_x",formula.Data(),fTree);
         gcut->SetObjectX(fx);

         fCodes[code] = -1;

      } else {

         Error("DefinedVariable","Found a TCutG without leaf information (%s)",
               gcut->GetName());
         return -1;

      }

      fExternalCuts.AddAtAndExpand(gcut,code);
      fNcodes++;
      fLookupType[code] = -1;
      return code;
   }

   //may be an entrylist
   TEntryList *elist = dynamic_cast<TEntryList*> (gDirectory->Get(name.Data()));
   if (elist) {
      Int_t code = fNcodes;
      fCodes[code] = 0;
      fExternalCuts.AddAtAndExpand(elist, code);
      fNcodes++;
      fLookupType[code] = kEntryList;
      return code;

   }

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the leaf (if any) which contains an object containing
/// a data member which has the name provided in the arguments.

TLeaf* TTreeFormula::GetLeafWithDatamember(const char* topchoice, const char* nextchoice, Long64_t readentry) const
{
   TClass * cl = 0;
   TIter nextleaf (fTree->GetIteratorOnAllLeaves());
   TFormLeafInfo* clonesinfo = 0;
   TLeaf *leafcur;
   while ((leafcur = (TLeaf*)nextleaf())) {
      // The following code is used somewhere else, we need to factor it out.

      // Here since we are interested in data member, we want to consider only
      // 'terminal' branch and leaf.
      cl = 0;
      if (leafcur->InheritsFrom(TLeafObject::Class()) &&
          leafcur->GetBranch()->GetListOfBranches()->Last()==0) {
         TLeafObject *lobj = (TLeafObject*)leafcur;
         cl = lobj->GetClass();
      } else if (leafcur->InheritsFrom(TLeafElement::Class()) && leafcur->IsOnTerminalBranch()) {
         TLeafElement * lElem = (TLeafElement*) leafcur;
         if (lElem->IsOnTerminalBranch()) {
            TBranchElement *branchEl = (TBranchElement *)leafcur->GetBranch();
            Int_t type = branchEl->GetStreamerType();
            if (type==-1) {
               cl = branchEl->GetInfo() ? branchEl->GetInfo()->GetClass() : 0;
            } else if (type>60 || type==0) {
               // Case of an object data member.  Here we allow for the
               // variable name to be omitted.  Eg, for Event.root with split
               // level 1 or above  Draw("GetXaxis") is the same as Draw("fH.GetXaxis()")
               TStreamerElement* element = branchEl->GetInfo()->GetElement(branchEl->GetID());
               if (element) cl = element->GetClassPointer();
               else cl = 0;
            }
         }

      }
      if (clonesinfo) { delete clonesinfo; clonesinfo = 0; }
      if (cl ==  TClonesArray::Class()) {
         // We have a unsplit TClonesArray leaves
         // In this case we assume that cl is the class in which the TClonesArray
         // belongs.
         R__LoadBranch(leafcur->GetBranch(),readentry,fQuickLoad);
         TClonesArray * clones;

         TBranch *branch = leafcur->GetBranch();
         if  (   branch->IsA()==TBranchElement::Class()
                 && ((TBranchElement*)branch)->GetType()==31) {

            // We have an unsplit TClonesArray as part of a split TClonesArray!

            // Let's not dig any further.  If the user really wants a data member
            // inside the nested TClonesArray, it has to specify it explicitly.

            continue;

         } else {
            Bool_t toplevel = (branch == branch->GetMother());
            clonesinfo = new TFormLeafInfoClones(cl, 0, toplevel);
            clones = (TClonesArray*)clonesinfo->GetLocalValuePointer(leafcur,0);
         }
         if (clones) cl = clones->GetClass();
      } else if (cl && cl->GetCollectionProxy()) {

         // We have a unsplit Collection leaves
         // In this case we assume that cl is the class in which the TClonesArray
         // belongs.

         TBranch *branch = leafcur->GetBranch();
         if  (   branch->IsA()==TBranchElement::Class()
                 && ((TBranchElement*)branch)->GetType()==41) {

            // We have an unsplit Collection as part of a split Collection!

            // Let's not dig any further.  If the user really wants a data member
            // inside the nested Collection, it has to specify it explicitly.

            continue;

         } else {
            clonesinfo = new TFormLeafInfoCollection(cl, 0);
         }
         cl = cl->GetCollectionProxy()->GetValueClass();
      }
      if (cl) {
         // Now that we have the class, let's check if the topchoice is of its datamember
         // or if the nextchoice is a datamember of one of its datamember.
         Int_t offset;
         TStreamerInfo* info =  (TStreamerInfo*)cl->GetStreamerInfo();
         TStreamerElement* element = info?info->GetStreamerElement(topchoice,offset):0;
         if (!element) {
            TIter nextel( cl->GetStreamerInfo()->GetElements() );
            TStreamerElement * curelem;
            while ((curelem = (TStreamerElement*)nextel())) {

               if (curelem->GetClassPointer() ==  TClonesArray::Class()) {
                  // In case of a TClonesArray we need to load the data and read the
                  // clonesArray object before being able to look into the class inside.
                  // We need to do that because we are never interested in the TClonesArray
                  // itself but only in the object inside.
                  TBranch *branch = leafcur->GetBranch();
                  TFormLeafInfo *leafinfo = 0;
                  if (clonesinfo) {
                     leafinfo = clonesinfo;
                  } else if (branch->IsA()==TBranchElement::Class()
                             && ((TBranchElement*)branch)->GetType()==31) {
                     // Case of a sub branch of a TClonesArray
                     TBranchElement *branchEl = (TBranchElement*)branch;
                     TStreamerInfo *bel_info = branchEl->GetInfo();
                     TClass * mother_cl = ((TBranchElement*)branch)->GetInfo()->GetClass();
                     TStreamerElement *bel_element =
                        bel_info->GetElement(branchEl->GetID());
                     leafinfo = new TFormLeafInfoClones(mother_cl, 0, bel_element, kTRUE);
                  }

                  Int_t clones_offset = 0;
                  ((TStreamerInfo*)cl->GetStreamerInfo())->GetStreamerElement(curelem->GetName(),clones_offset);
                  TFormLeafInfo* sub_clonesinfo = new TFormLeafInfo(cl, clones_offset, curelem);
                  if (leafinfo)
                     if (leafinfo->fNext) leafinfo->fNext->fNext = sub_clonesinfo;
                     else leafinfo->fNext = sub_clonesinfo;
                  else leafinfo = sub_clonesinfo;

                  R__LoadBranch(branch,readentry,fQuickLoad);

                  TClonesArray * clones = (TClonesArray*)leafinfo->GetValuePointer(leafcur,0);

                  delete leafinfo; clonesinfo = 0;
                  // If TClonesArray object does not exist we have no information, so let go
                  // on.  This is a weakish test since the TClonesArray object might exist in
                  // the next entry ... In other word, we ONLY rely on the information available
                  // in entry #0.
                  if (!clones) continue;
                  TClass *sub_cl = clones->GetClass();

                  // Now that we finally have the inside class, let's query it.
                  element = ((TStreamerInfo*)sub_cl->GetStreamerInfo())->GetStreamerElement(nextchoice,offset);
                  if (element) break;
               } // if clones array
               else if (curelem->GetClassPointer() && curelem->GetClassPointer()->GetCollectionProxy()) {

                  TClass *sub_cl = curelem->GetClassPointer()->GetCollectionProxy()->GetValueClass();

                  while(sub_cl && sub_cl->GetCollectionProxy())
                     sub_cl = sub_cl->GetCollectionProxy()->GetValueClass();

                  // Now that we finally have the inside class, let's query it.
                  if (sub_cl) element = ((TStreamerInfo*)sub_cl->GetStreamerInfo())->GetStreamerElement(nextchoice,offset);
                  if (element) break;

               }
            } // loop on elements
         }
         if (element) break;
         else cl = 0;
      }
   }
   delete clonesinfo;
   if (cl) {
      return leafcur;
   } else {
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the leaf (if any) of the tree with contains an object of a class
/// having a method which has the name provided in the argument.

Bool_t TTreeFormula::BranchHasMethod(TLeaf* leafcur, TBranch* branch, const char* method, const char* params, Long64_t readentry) const
{
   TClass *cl = 0;
   TLeafObject* lobj = 0;

   // Since the user does not want this branch to be loaded anyway, we just
   // skip it.  This prevents us from warning the user that the method might
   // be on a disabled branch.  However, and more usefully, this allows the
   // user to avoid error messages from branches that cannot be currently
   // read without warnings/errors.

   if (branch->TestBit(kDoNotProcess)) {
      return kFALSE;
   }

   // FIXME: The following code is used somewhere else, we need to factor it out.
   if (branch->InheritsFrom(TBranchObject::Class())) {
      lobj = (TLeafObject*) branch->GetListOfLeaves()->At(0);
      cl = lobj->GetClass();
   } else if (branch->InheritsFrom(TBranchElement::Class())) {
      TBranchElement* branchEl = (TBranchElement*) branch;
      Int_t type = branchEl->GetStreamerType();
      if (type == -1) {
         cl = branchEl->GetInfo() ? branchEl->GetInfo()->GetClass() : 0;
      } else if (type > 60) {
         // Case of an object data member.  Here we allow for the
         // variable name to be omitted.  Eg, for Event.root with split
         // level 1 or above  Draw("GetXaxis") is the same as Draw("fH.GetXaxis()")
         TStreamerElement* element = branchEl->GetInfo()->GetElement(branchEl->GetID());
         if (element) {
            cl = element->GetClassPointer();
         } else {
            cl = 0;
         }
         if ((cl == TClonesArray::Class()) && (branchEl->GetType() == 31)) {
            // we have a TClonesArray inside a split TClonesArray,
            // Let's not dig any further.  If the user really wants a data member
            // inside the nested TClonesArray, it has to specify it explicitly.
            cl = 0;
         }
         // NOTE do we need code for Collection here?
      }
   }

   if (cl == TClonesArray::Class()) {
      // We might be try to call a method of the top class inside a
      // TClonesArray.
      // Since the leaf was not terminal, we might have a split or
      // unsplit and/or top leaf/branch.
      TClonesArray* clones = 0;
      R__LoadBranch(branch, readentry, fQuickLoad);
      if (branch->InheritsFrom(TBranchObject::Class())) {
         clones = (TClonesArray*) lobj->GetObject();
      } else if (branch->InheritsFrom(TBranchElement::Class())) {
         // We do not know exactly where the leaf of the TClonesArray is
         // in the hierarchy but we still need to get the correct class
         // holder.
         TBranchElement* bc = (TBranchElement*) branch;
         if (bc == bc->GetMother()) {
            // Top level branch
            //clones = *((TClonesArray**) bc->GetAddress());
            clones = (TClonesArray*) bc->GetObject();
         } else if (!leafcur || !leafcur->IsOnTerminalBranch()) {
            TStreamerElement* element = bc->GetInfo()->GetElement(bc->GetID());
            if (element->IsaPointer()) {
               clones = *((TClonesArray**) bc->GetAddress());
               //clones = *((TClonesArray**) bc->GetObject());
            } else {
               //clones = (TClonesArray*) bc->GetAddress();
               clones = (TClonesArray*) bc->GetObject();
            }
         }
         if (!clones) {
            R__LoadBranch(bc, readentry, fQuickLoad);
            TClass* mother_cl;
            mother_cl = bc->GetInfo()->GetClass();
            TFormLeafInfo* clonesinfo = new TFormLeafInfoClones(mother_cl, 0);
            // if (!leafcur) { leafcur = (TLeaf*) branch->GetListOfLeaves()->At(0); }
            clones = (TClonesArray*) clonesinfo->GetLocalValuePointer(leafcur, 0);
            // cl = clones->GetClass();
            delete clonesinfo;
         }
      } else {
         Error("BranchHasMethod","A TClonesArray was stored in a branch type no yet support (i.e. neither TBranchObject nor TBranchElement): %s",branch->IsA()->GetName());
         return kFALSE;
      }
      cl = clones ? clones->GetClass() : 0;
   } else if (cl && cl->GetCollectionProxy()) {
      cl = cl->GetCollectionProxy()->GetValueClass();
   }

   if (cl) {
      if (cl->GetClassInfo()) {
         if (cl->GetMethodAllAny(method)) {
            // Let's try to see if the function we found belongs to the current
            // class.  Note that this implementation currently can not work if
            // one the argument is another leaf or data member of the object.
            // (Anyway we do NOT support this case).
            TMethodCall methodcall(cl, method, params);
            if (methodcall.GetMethod()) {
               // We have a method that works.
               // We will use it.
               return kTRUE;
            }
         }
      }
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Now let calculate what physical instance we really need.
/// Some redundant code is used to speed up the cases where
/// they are no dimensions.
///
/// We know that instance is less that fCumulUsedSize[0] so
/// we can skip the modulo when virt_dim is 0.

Int_t TTreeFormula::GetRealInstance(Int_t instance, Int_t codeindex) {
      Int_t real_instance = 0;
      Int_t virt_dim;

      Bool_t check = kFALSE;
      if (codeindex<0) {
         codeindex = 0;
         check = kTRUE;
      }

      TFormLeafInfo * info = 0;
      Int_t max_dim = fNdimensions[codeindex];
      if ( max_dim ) {
         virt_dim = 0;
         max_dim--;

         if (!fManager->fMultiVarDim) {
            if (fIndexes[codeindex][0]>=0) {
               real_instance = fIndexes[codeindex][0] * fCumulSizes[codeindex][1];
            } else {
               Int_t local_index;
               local_index = ( instance / fManager->fCumulUsedSizes[virt_dim+1]);
               if (fIndexes[codeindex][0]==-2) {
                  // NOTE: Should we check that this is a valid index?
                  if (check) {
                     Int_t index_real_instance = fVarIndexes[codeindex][0]->GetRealInstance(local_index,-1);
                     if (index_real_instance >= fVarIndexes[codeindex][0]->fNdata[0]) {
                        // out of bounds
                        return fNdata[0]+1;
                     }
                  }
                  if (fDidBooleanOptimization && local_index!=0) {
                     // Force the loading of the index.
                     fVarIndexes[codeindex][0]->LoadBranches();
                  }
                  local_index = (Int_t)fVarIndexes[codeindex][0]->EvalInstance(local_index);
                  if (local_index<0) {
                     Error("EvalInstance","Index %s is out of bound (%d) in formula %s",
                           fVarIndexes[codeindex][0]->GetTitle(),
                           local_index,
                           GetTitle());
                     return fNdata[0]+1;
                  }
               }
               real_instance = local_index * fCumulSizes[codeindex][1];
               virt_dim ++;
            }
         } else {
            // NOTE: We assume that ONLY the first dimension of a leaf can have a variable
            // size AND contain the index for the size of yet another sub-dimension.
            // I.e. a variable size array inside a variable size array can only have its
            // size vary with the VERY FIRST physical dimension of the leaf.
            // Thus once the index of the first dimension is found, all other dimensions
            // are fixed!

            // NOTE: We could unroll some of this loops to avoid a few tests.
            if (fHasMultipleVarDim[codeindex]) {
               info = (TFormLeafInfo *)(fDataMembers.At(codeindex));
               // if (info && info->GetVarDim()==-1) info = 0;
            }
            Int_t local_index;

            switch (fIndexes[codeindex][0]) {
            case -2:
               if (fDidBooleanOptimization && instance!=0) {
                  // Force the loading of the index.
                  fVarIndexes[codeindex][0]->LoadBranches();
               }
               local_index = (Int_t)fVarIndexes[codeindex][0]->EvalInstance(instance);
               if (local_index<0) {
                  Error("EvalInstance","Index %s is out of bound (%d) in formula %s",
                        fVarIndexes[codeindex][0]->GetTitle(),
                        local_index,
                        GetTitle());
                  local_index = 0;
               }
               break;
            case -1: {
                  if (instance <= fRealInstanceCache.fInstanceCache) {
                     fRealInstanceCache.fLocalIndexCache = 0;
                     fRealInstanceCache.fVirtAccumCache = 0;
                  }
                  fRealInstanceCache.fInstanceCache = instance;
                  local_index = fRealInstanceCache.fLocalIndexCache;
                  Int_t virt_accum = fRealInstanceCache.fVirtAccumCache;

                  Int_t maxloop = fManager->fCumulUsedVarDims->GetSize();
                  if (maxloop == 0) {
                     local_index--;
                     instance = fNdata[0]+1; // out of bounds.
                     if (check) return fNdata[0]+1;
                  } else {
                     do {
                        virt_accum += fManager->fCumulUsedVarDims->GetArray()[local_index];
                        local_index++;
                     } while( instance >= virt_accum && local_index<maxloop);
                     local_index--;
                     // update the cache
                     fRealInstanceCache.fVirtAccumCache = virt_accum - fManager->fCumulUsedVarDims->GetArray()[local_index];
                     fRealInstanceCache.fLocalIndexCache = local_index;

                     if (local_index==(maxloop-1) && (instance >= virt_accum)) {
                        instance = fNdata[0]+1; // out of bounds.
                        if (check) return fNdata[0]+1;
                     } else {
                        if (fManager->fCumulUsedVarDims->At(local_index)) {
                           instance -= (virt_accum - fManager->fCumulUsedVarDims->At(local_index));
                        } else {
                           instance = fNdata[0]+1; // out of bounds.
                           if (check) return fNdata[0]+1;
                        }
                     }
                  }
                  virt_dim ++;
               }
               break;
            default:
               local_index = fIndexes[codeindex][0];
            }

            // Inform the (appropriate) MultiVarLeafInfo that the clones array index is
            // local_index.

            if (fManager->fVarDims[kMAXFORMDIM]) {
               fManager->fCumulUsedSizes[kMAXFORMDIM] = fManager->fVarDims[kMAXFORMDIM]->At(local_index);
            } else {
               fManager->fCumulUsedSizes[kMAXFORMDIM] = fManager->fUsedSizes[kMAXFORMDIM];
            }
            for(Int_t d = kMAXFORMDIM-1; d>0; d--) {
               if (fManager->fVarDims[d]) {
                  fManager->fCumulUsedSizes[d] = fManager->fCumulUsedSizes[d+1] * fManager->fVarDims[d]->At(local_index);
               } else {
                  fManager->fCumulUsedSizes[d] = fManager->fCumulUsedSizes[d+1] * fManager->fUsedSizes[d];
               }
            }
            if (info) {
               // When we have multiple variable dimensions, the LeafInfo only expect
               // the instance after the primary index has been set.
               info->SetPrimaryIndex(local_index);
               real_instance = 0;

               // Let's update fCumulSizes for the rest of the code.
               Int_t vdim = info->GetVarDim();
               Int_t isize = info->GetSize(local_index);
               if (fIndexes[codeindex][vdim]>=0) {
                  info->SetSecondaryIndex(fIndexes[codeindex][vdim]);
               }
               if  (isize!=1 && fIndexes[codeindex][vdim]>isize) {
                  // We are out of bounds!
                  return fNdata[0]+1;
               }
               fCumulSizes[codeindex][vdim] =  isize*fCumulSizes[codeindex][vdim+1];
               for(Int_t k=vdim -1; k>0; --k) {
                  fCumulSizes[codeindex][k] = fCumulSizes[codeindex][k+1]*fFixedSizes[codeindex][k];
               }
            } else {
               real_instance = local_index * fCumulSizes[codeindex][1];
            }
         }
         if (max_dim>0) {
            for (Int_t dim = 1; dim < max_dim; dim++) {
               if (fIndexes[codeindex][dim]>=0) {
                  real_instance += fIndexes[codeindex][dim] * fCumulSizes[codeindex][dim+1];
               } else {
                  Int_t local_index;
                  if (virt_dim && fManager->fCumulUsedSizes[virt_dim]>1) {
                     local_index = ( ( instance % fManager->fCumulUsedSizes[virt_dim] )
                                     / fManager->fCumulUsedSizes[virt_dim+1]);
                  } else {
                     local_index = ( instance / fManager->fCumulUsedSizes[virt_dim+1]);
                  }
                  if (fIndexes[codeindex][dim]==-2) {
                     // NOTE: Should we check that this is a valid index?
                     if (fDidBooleanOptimization && local_index!=0) {
                        // Force the loading of the index.
                        fVarIndexes[codeindex][dim]->LoadBranches();
                     }
                     local_index = (Int_t)fVarIndexes[codeindex][dim]->EvalInstance(local_index);
                     if (local_index<0 ||
                         local_index>=(fCumulSizes[codeindex][dim]/fCumulSizes[codeindex][dim+1])) {
                        Error("EvalInstance","Index %s is out of bound (%d/%d) in formula %s",
                              fVarIndexes[codeindex][dim]->GetTitle(),
                              local_index,
                              (fCumulSizes[codeindex][dim]/fCumulSizes[codeindex][dim+1]),
                              GetTitle());
                        local_index = (fCumulSizes[codeindex][dim]/fCumulSizes[codeindex][dim+1])-1;
                     }
                  }
                  real_instance += local_index * fCumulSizes[codeindex][dim+1];
                  virt_dim ++;
               }
            }
            if (fIndexes[codeindex][max_dim]>=0) {
               if (!info) real_instance += fIndexes[codeindex][max_dim];
            } else {
               Int_t local_index;
               if (virt_dim && fManager->fCumulUsedSizes[virt_dim]>1) {
                  local_index = instance % fManager->fCumulUsedSizes[virt_dim];
               } else {
                  local_index = instance;
               }
               if (info && local_index>=fCumulSizes[codeindex][max_dim]) {
                  // We are out of bounds! [Multiple var dims, See same message a few line above]
                  return fNdata[0]+1;
               }
               if (fIndexes[codeindex][max_dim]==-2) {
                  if (fDidBooleanOptimization && local_index!=0) {
                     // Force the loading of the index.
                     fVarIndexes[codeindex][max_dim]->LoadBranches();
                  }
                  local_index = (Int_t)fVarIndexes[codeindex][max_dim]->EvalInstance(local_index);
                  if (local_index<0 ||
                         local_index>=fCumulSizes[codeindex][max_dim]) {
                     Error("EvalInstance","Index %s is of out bound (%d/%d) in formula %s",
                           fVarIndexes[codeindex][max_dim]->GetTitle(),
                           local_index,
                           fCumulSizes[codeindex][max_dim],
                           GetTitle());
                     local_index = fCumulSizes[codeindex][max_dim]-1;
                  }
               }
               real_instance += local_index;
            }
         } // if (max_dim-1>0)
      } // if (max_dim)

      return real_instance;
}

////////////////////////////////////////////////////////////////////////////////
///  Evaluate the class of this treeformula.
///
///  If the 'value' of this formula is a simple pointer to an object,
///  this function returns the TClass corresponding to its type.

TClass* TTreeFormula::EvalClass() const
{
   if (fNoper != 1 || fNcodes <=0 ) return 0;

   return EvalClass(0);
}

////////////////////////////////////////////////////////////////////////////////
///  Evaluate the class of the operation oper.
///
///  If the 'value' in the requested operation is a simple pointer to an object,
///  this function returns the TClass corresponding to its type.

TClass* TTreeFormula::EvalClass(Int_t oper) const
{
   TLeaf *leaf = (TLeaf*)fLeaves.UncheckedAt(oper);
   switch(fLookupType[oper]) {
      case kDirect: {
         if (leaf->IsA()==TLeafObject::Class()) {
            return ((TLeafObject*)leaf)->GetClass();
         } else if ( leaf->IsA()==TLeafElement::Class()) {
            TBranchElement * branch = (TBranchElement*)((TLeafElement*)leaf)->GetBranch();
            TStreamerInfo * info = branch->GetInfo();
            Int_t id = branch->GetID();
            if (id>=0) {
               if (info==0 || !info->IsCompiled()) {
                  // we probably do not have a way to know the class of the object.
                  return 0;
               }
               TStreamerElement* elem = (TStreamerElement*)info->GetElement(id);
               if (elem==0) {
                  // we probably do not have a way to know the class of the object.
                  return 0;
               } else {
                  return elem->GetClass();
               }
            } else return TClass::GetClass( branch->GetClassName() );
         } else {
            return 0;
         }
      }
      case kMethod: return 0; // kMethod is deprecated so let's no waste time implementing this.
      case kTreeMember:
      case kDataMember: {
         TObject *obj = fDataMembers.UncheckedAt(oper);
         if (!obj) return 0;
         return ((TFormLeafInfo*)obj)->GetClass();
      }

      default: return 0;
   }


}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate this treeformula.
///
///  Return the address of the object pointed to by the formula.
///  Return 0 if the formula is not a single object
///  The object type can be retrieved using by call EvalClass();

void* TTreeFormula::EvalObject(int instance)
{
   if (fNoper != 1 || fNcodes <=0 ) return 0;


   switch (fLookupType[0]) {
      case kIndexOfEntry:
      case kIndexOfLocalEntry:
      case kEntries:
      case kLocalEntries:
      case kLength:
      case kLengthFunc:
      case kIteration:
      case kEntryList:
         return 0;
   }

   TLeaf *leaf = (TLeaf*)fLeaves.UncheckedAt(0);

   Int_t real_instance = GetRealInstance(instance,0);

   if (instance==0 || fNeedLoading) {
      fNeedLoading = kFALSE;
      R__LoadBranch(leaf->GetBranch(),
                    leaf->GetBranch()->GetTree()->GetReadEntry(),
                    fQuickLoad);
   }
   else if (real_instance>=fNdata[0]) return 0;
   if (fAxis) {
      return 0;
   }
   switch(fLookupType[0]) {
      case kDirect: {
         if (real_instance) {
            Warning("EvalObject","Not yet implement for kDirect and arrays (for %s).\nPlease contact the developers",GetName());
         }
         return leaf->GetValuePointer();
      }
      case kMethod: return GetValuePointerFromMethod(0,leaf);
      case kTreeMember:
      case kDataMember: return ((TFormLeafInfo*)fDataMembers.UncheckedAt(0))->GetValuePointer(leaf,real_instance);
      default: return 0;
   }


}


////////////////////////////////////////////////////////////////////////////////
/// Eval the instance as a string.

const char* TTreeFormula::EvalStringInstance(Int_t instance)
{
   const Int_t kMAXSTRINGFOUND = 10;
   const char *stringStack[kMAXSTRINGFOUND];

   if (fNoper==1 && fNcodes>0 && IsString()) {
      TLeaf *leaf = (TLeaf*)fLeaves.UncheckedAt(0);

      Int_t real_instance = GetRealInstance(instance,0);

      if (instance==0 || fNeedLoading) {
         fNeedLoading = kFALSE;
         TBranch *branch = leaf->GetBranch();
         R__LoadBranch(branch,branch->GetTree()->GetReadEntry(),fQuickLoad);
      } else if (real_instance>=fNdata[0]) {
         return 0;
      }

      if (fLookupType[0]==kDirect) {
         return (char*)leaf->GetValuePointer();
      } else {
         return  (char*)GetLeafInfo(0)->GetValuePointer(leaf,real_instance);
      }
   }

   EvalInstance(instance,stringStack);

   return stringStack[0];
}

#define TT_EVAL_INIT                                                                            \
   TLeaf *leaf = (TLeaf*)fLeaves.UncheckedAt(0);                                                \
                                                                                                \
   const Int_t real_instance = GetRealInstance(instance,0);                                     \
                                                                                                \
   if (instance==0) fNeedLoading = kTRUE;                                                       \
   if (real_instance>=fNdata[0]) return 0;                                                       \
                                                                                                \
   /* Since the only operation in this formula is reading this branch,                          \
      we are guaranteed that this function is first called with instance==0 and                 \
      hence we are guaranteed that the branch is always properly read */                        \
                                                                                                \
   if (fNeedLoading) {                                                                          \
      fNeedLoading = kFALSE;                                                                    \
      TBranch *br = leaf->GetBranch();                                                          \
      Long64_t tentry = br->GetTree()->GetReadEntry();                                          \
      R__LoadBranch(br,tentry,fQuickLoad);                                                      \
   }                                                                                            \
                                                                                                \
   if (fAxis) {                                                                                 \
      char * label;                                                                             \
      /* This portion is a duplicate (for speed reason) of the code                             \
         located  in the main for loop at "a tree string" (and in EvalStringInstance) */        \
      if (fLookupType[0]==kDirect) {                                                            \
         label = (char*)leaf->GetValuePointer();                                                \
      } else {                                                                                  \
         label = (char*)GetLeafInfo(0)->GetValuePointer(leaf,instance);                         \
      }                                                                                         \
      Int_t bin = fAxis->FindBin(label);                                                        \
      return bin-0.5;                                                                           \
   }

#define TREE_EVAL_INIT                                                                          \
   const Int_t real_instance = GetRealInstance(instance,0);                                     \
                                                                                                \
   if (real_instance>=fNdata[0]) return 0;                                                      \
                                                                                                \
   if (fAxis) {                                                                                 \
      char * label;                                                                             \
      /* This portion is a duplicate (for speed reason) of the code                             \
         located  in the main for loop at "a tree string" (and in EvalStringInstance) */        \
      label = (char*)GetLeafInfo(0)->GetValuePointer((TLeaf*)0x0,instance);                     \
      Int_t bin = fAxis->FindBin(label);                                                        \
      return bin-0.5;                                                                           \
   }

#define TT_EVAL_INIT_LOOP                                                                       \
   TLeaf *leaf = (TLeaf*)fLeaves.UncheckedAt(code);                                             \
                                                                                                \
   /* Now let calculate what physical instance we really need.  */                              \
   const Int_t real_instance = GetRealInstance(instance,code);                                  \
                                                                                                \
   if (willLoad) {                                                                              \
      TBranch *branch = (TBranch*)fBranches.UncheckedAt(code);                                  \
      if (branch) {                                                                             \
         Long64_t treeEntry = branch->GetTree()->GetReadEntry();                                \
         R__LoadBranch(branch,treeEntry,fQuickLoad);                                            \
      } else if (fDidBooleanOptimization) {                                                     \
         branch = leaf->GetBranch();                                                            \
         Long64_t treeEntry = branch->GetTree()->GetReadEntry();                                \
         if (branch->GetReadEntry() != treeEntry) branch->GetEntry( treeEntry );                \
      }                                                                                         \
   } else {                                                                                     \
      /* In the cases where we are behind (i.e. right of) a potential boolean optimization      \
         this tree variable reading may have not been executed with instance==0 which would     \
         result in the branch being potentially not read in. */                                 \
      if (fDidBooleanOptimization) {                                                            \
         TBranch *br = leaf->GetBranch();                                                       \
         Long64_t treeEntry = br->GetTree()->GetReadEntry();                                    \
         if (br->GetReadEntry() != treeEntry) br->GetEntry( treeEntry );                        \
      }                                                                                         \
   }                                                                                            \
   if (real_instance>=fNdata[code]) return 0;

#define TREE_EVAL_INIT_LOOP                                                                     \
   /* Now let calculate what physical instance we really need.  */                              \
   const Int_t real_instance = GetRealInstance(instance,code);                                  \
                                                                                                \
   if (real_instance>=fNdata[code]) return 0;


template<typename T> T Summing(TTreeFormula *sum) {
   Int_t len = sum->GetNdata();
   T res = 0;
   for (int i=0; i<len; ++i) res += sum->EvalInstance<T>(i);
   return res;
}

template<typename T> T  FindMin(TTreeFormula *arr) {
   Int_t len = arr->GetNdata();
   T res = 0;
   if (len) {
      res = arr->EvalInstance<T>(0);
      for (int i=1; i<len; ++i) {
         T val = arr->EvalInstance<T>(i);
         if (val < res) {
            res = val;
         }
      }
   }
   return res;
}

template<typename T> T FindMax(TTreeFormula *arr) {
   Int_t len = arr->GetNdata();
   T res = 0;
   if (len) {
      res = arr->EvalInstance<T>(0);
      for (int i=1; i<len; ++i) {
         T val = arr->EvalInstance(i);
         if (val > res) {
            res = val;
         }
      }
   }
   return res;
}

template<typename T> T FindMin(TTreeFormula *arr, TTreeFormula *condition) {
   Int_t len = arr->GetNdata();
   T res = 0;
   if (len) {
      int i = 0;
      T condval;
      do {
         condval = condition->EvalInstance<T>(i);
         ++i;
      } while (!condval && i<len);
      if (!condval && i==len) {
         return 0;
      }
      if (i!=1) {
         // Insure the loading of the branch.
         arr->EvalInstance<T>(0);
      }
      // Now we know that i>0 && i<len and cond==true
      res = arr->EvalInstance<T>(i-1);
      for (; i<len; ++i) {
         condval = condition->EvalInstance<T>(i);
         if (condval) {
            T val = arr->EvalInstance<T>(i);
            if (val < res) {
               res = val;
            }
         }
      }
   }
   return res;
}

template<typename T> T FindMax(TTreeFormula *arr, TTreeFormula *condition) {
   Int_t len = arr->GetNdata();
   T res = 0;
   if (len) {
      int i = 0;
      T condval;
      do {
         condval = condition->EvalInstance<T>(i);
         ++i;
      } while (!condval && i<len);
      if (!condval && i==len) {
         return 0;
      }
      if (i!=1) {
         // Insure the loading of the branch.
         arr->EvalInstance<T>(0);
      }
      // Now we know that i>0 && i<len and cond==true
      res = arr->EvalInstance<T>(i-1);
      for (; i<len; ++i) {
         condval = condition->EvalInstance<T>(i);
         if (condval) {
            T val = arr->EvalInstance<T>(i);
            if (val > res) {
               res = val;
            }
         }
      }
   }
   return res;
}

namespace {

template <typename T> T fmod_local(T x, T y) { return fmod(x,y); }
template <> Long64_t fmod_local(Long64_t x, Long64_t y) { return fmod((LongDouble_t)x,(LongDouble_t)y); }

template<typename T> inline void SetMethodParam(TMethodCall *method, T p) { method->SetParam(p); }
template<> void SetMethodParam(TMethodCall *method, LongDouble_t p) { method->SetParam((Double_t)p); }

}

template<typename T> inline T TTreeFormula::GetConstant(Int_t k) { return fConst[k]; }
template<> inline LongDouble_t TTreeFormula::GetConstant(Int_t k) {
   if( !fConstLD ) {
      // create LD version of the constants list by rescanning all literals used in the expression
      fConstLD = new LongDouble_t[fNconst];
      for (Int_t op=0; op<fNoper ; ++op) {
         const Int_t oper = GetOper()[op];
         if( (oper >> kTFOperShift) == kConstant ) {
            int i = (oper & kTFOperMask);
            if( !strncmp(fExpr[op], "0x", 2) || !strncmp(fExpr[op], "0X", 2) ) {
               ULong64_t val;
               sscanf( fExpr[op], "%llx", &val );
               fConstLD[i] = (LongDouble_t)val;
            } else {
               sscanf( fExpr[op], "%Lg", &fConstLD[i] );
            }
         }
      }
   }
   return fConstLD[k];
}
template<> inline Long64_t TTreeFormula::GetConstant(Int_t k) { return (Long64_t)GetConstant<LongDouble_t>(k); }

////////////////////////////////////////////////////////////////////////////////
/// Evaluate this treeformula.

template<typename T>
T TTreeFormula::EvalInstance(Int_t instance, const char *stringStackArg[])
{
// Note that the redundancy and structure in this code is tailored to improve
// efficiencies.
   if (TestBit(kMissingLeaf)) return 0;
   if (fNoper == 1 && fNcodes > 0) {

      switch (fLookupType[0]) {
         case kDirect:     {
            TT_EVAL_INIT;
            return leaf->GetTypedValue<T>(real_instance);
         }
         case kMethod:     {
            TT_EVAL_INIT;
            ((TFormLeafInfo*)fDataMembers.UncheckedAt(0))->SetBranch(leaf->GetBranch());
            return GetValueFromMethod(0,leaf);
         }
         case kDataMember: {
            TT_EVAL_INIT;
            ((TFormLeafInfo*)fDataMembers.UncheckedAt(0))->SetBranch(leaf->GetBranch());
            return ((TFormLeafInfo*)fDataMembers.UncheckedAt(0))->GetTypedValue<T>(leaf,real_instance);
         }
         case kTreeMember: {
            TREE_EVAL_INIT;
            return ((TFormLeafInfo*)fDataMembers.UncheckedAt(0))->GetTypedValue<T>((TLeaf*)0x0,real_instance);
         }
         case kIndexOfEntry: return (T)fTree->GetReadEntry();
         case kIndexOfLocalEntry: return (T)fTree->GetTree()->GetReadEntry();
         case kEntries:      return (T)fTree->GetEntries();
         case kLocalEntries: return (T)fTree->GetTree()->GetEntries();
         case kLength:       return fManager->fNdata;
         case kLengthFunc:   return ((TTreeFormula*)fAliases.UncheckedAt(0))->GetNdata();
         case kIteration:    return instance;
         case kSum:          return Summing<T>((TTreeFormula*)fAliases.UncheckedAt(0));
         case kMin:          return FindMin<T>((TTreeFormula*)fAliases.UncheckedAt(0));
         case kMax:          return FindMax<T>((TTreeFormula*)fAliases.UncheckedAt(0));
         case kEntryList: {
            TEntryList *elist = (TEntryList*)fExternalCuts.At(0);
            return elist->Contains(fTree->GetTree()->GetReadEntry());
         }
         case -1: break;
      }
      switch (fCodes[0]) {
         case -2: {
            TCutG *gcut = (TCutG*)fExternalCuts.At(0);
            TTreeFormula *fx = (TTreeFormula *)gcut->GetObjectX();
            TTreeFormula *fy = (TTreeFormula *)gcut->GetObjectY();
            if (fDidBooleanOptimization) {
               fx->ResetLoading();
               fy->ResetLoading();
            }
            T xcut = fx->EvalInstance<T>(instance);
            T ycut = fy->EvalInstance<T>(instance);
            return gcut->IsInside(xcut,ycut);
         }
         case -1: {
            TCutG *gcut = (TCutG*)fExternalCuts.At(0);
            TTreeFormula *fx = (TTreeFormula *)gcut->GetObjectX();
            if (fDidBooleanOptimization) {
               fx->ResetLoading();
            }
            return fx->EvalInstance<T>(instance);
         }
         default: return 0;
      }
   }

   T tab[kMAXFOUND];
   const Int_t kMAXSTRINGFOUND = 10;
   const char *stringStackLocal[kMAXSTRINGFOUND];
   const char **stringStack = stringStackArg?stringStackArg:stringStackLocal;

   const Bool_t willLoad = (instance==0 || fNeedLoading); fNeedLoading = kFALSE;
   if (willLoad) fDidBooleanOptimization = kFALSE;

   Int_t pos  = 0;
   Int_t pos2 = 0;
   for (Int_t i=0; i<fNoper ; ++i) {

      const Int_t oper = GetOper()[i];
      const Int_t newaction = oper >> kTFOperShift;

      if (newaction<kDefinedVariable) {
         // ROOT::v5::TFormula operands.

         // one of the most used cases
         if (newaction==kConstant) { pos++; tab[pos-1] = GetConstant<T>(oper & kTFOperMask); continue; }

         switch(newaction) {

            case kEnd        : return tab[0];
            case kAdd        : pos--; tab[pos-1] += tab[pos]; continue;
            case kSubstract  : pos--; tab[pos-1] -= tab[pos]; continue;
            case kMultiply   : pos--; tab[pos-1] *= tab[pos]; continue;
            case kDivide     : pos--; if (tab[pos] == 0) tab[pos-1] = 0; //  division by 0
                                      else               tab[pos-1] /= tab[pos];
                                      continue;
            case kModulo     : {pos--;
                                Long64_t int1((Long64_t)tab[pos-1]);
                                Long64_t int2((Long64_t)tab[pos]);
                                tab[pos-1] = T(int1 % int2);
                                continue;}

            case kcos  : tab[pos-1] = TMath::Cos(tab[pos-1]); continue;
            case ksin  : tab[pos-1] = TMath::Sin(tab[pos-1]); continue;
            case ktan  : if (TMath::Cos(tab[pos-1]) == 0) {tab[pos-1] = 0;} // { tangente indeterminee }
                         else tab[pos-1] = TMath::Tan(tab[pos-1]);
                         continue;
            case kacos : if (TMath::Abs(tab[pos-1]) > 1) {tab[pos-1] = 0;} //  indetermination
                         else tab[pos-1] = TMath::ACos(tab[pos-1]);
                         continue;
            case kasin : if (TMath::Abs(tab[pos-1]) > 1) {tab[pos-1] = 0;} //  indetermination
                         else tab[pos-1] = TMath::ASin(tab[pos-1]);
                         continue;
            case katan : tab[pos-1] = TMath::ATan(tab[pos-1]); continue;
            case kcosh : tab[pos-1] = TMath::CosH(tab[pos-1]); continue;
            case ksinh : tab[pos-1] = TMath::SinH(tab[pos-1]); continue;
            case ktanh : if (TMath::CosH(tab[pos-1]) == 0) {tab[pos-1] = 0;} // { tangente indeterminee }
                         else tab[pos-1] = TMath::TanH(tab[pos-1]);
                         continue;
            case kacosh: if (tab[pos-1] < 1) {tab[pos-1] = 0;} //  indetermination
                         else tab[pos-1] = TMath::ACosH(tab[pos-1]);
                         continue;
            case kasinh: tab[pos-1] = TMath::ASinH(tab[pos-1]); continue;
            case katanh: if (TMath::Abs(tab[pos-1]) > 1) {tab[pos-1] = 0;} // indetermination
                         else tab[pos-1] = TMath::ATanH(tab[pos-1]);
                         continue;
            case katan2: pos--; tab[pos-1] = TMath::ATan2(tab[pos-1],tab[pos]); continue;

            case kfmod : pos--; tab[pos-1] = fmod_local(tab[pos-1],tab[pos]); continue;
            case kpow  : pos--; tab[pos-1] = TMath::Power(tab[pos-1],tab[pos]); continue;
            case ksq   : tab[pos-1] = tab[pos-1]*tab[pos-1]; continue;
            case ksqrt : tab[pos-1] = TMath::Sqrt(TMath::Abs(tab[pos-1])); continue;

            case kstrstr : pos2 -= 2; pos++;if (strstr(stringStack[pos2],stringStack[pos2+1])) tab[pos-1]=1;
                                            else tab[pos-1]=0;
                           continue;

            case kmin : pos--; tab[pos-1] = std::min(tab[pos-1],tab[pos]); continue;
            case kmax : pos--; tab[pos-1] = std::max(tab[pos-1],tab[pos]); continue;

            case klog  : if (tab[pos-1] > 0) tab[pos-1] = TMath::Log(tab[pos-1]);
                         else {tab[pos-1] = 0;} //{indetermination }
                          continue;
            case kexp  : { Double_t dexp = tab[pos-1];
                           if (dexp < -700) {tab[pos-1] = 0; continue;}
                           if (dexp >  700) {tab[pos-1] = TMath::Exp(700); continue;}
                           tab[pos-1] = TMath::Exp(dexp); continue;
                         }
            case klog10: if (tab[pos-1] > 0) tab[pos-1] = TMath::Log10(tab[pos-1]);
                         else {tab[pos-1] = 0;} //{indetermination }
                         continue;

            case kpi   : pos++; tab[pos-1] = TMath::ACos(-1); continue;

            case kabs  : tab[pos-1] = TMath::Abs(tab[pos-1]); continue;
            case ksign : if (tab[pos-1] < 0) tab[pos-1] = -1; else tab[pos-1] = 1;
                         continue;
            case kint  : tab[pos-1] = T(Long64_t(tab[pos-1])); continue;
            case kSignInv: tab[pos-1] = -1 * tab[pos-1]; continue;
            case krndm : pos++; tab[pos-1] = gRandom->Rndm(); continue;

            case kAnd  : pos--; if (tab[pos-1]!=0 && tab[pos]!=0) tab[pos-1]=1;
                                else tab[pos-1]=0;
                         continue;
            case kOr   : pos--; if (tab[pos-1]!=0 || tab[pos]!=0) tab[pos-1]=1;
                                else tab[pos-1]=0;
                         continue;

            case kEqual      : pos--; tab[pos-1] = (tab[pos-1] == tab[pos]) ? 1 : 0; continue;
            case kNotEqual   : pos--; tab[pos-1] = (tab[pos-1] != tab[pos]) ? 1 : 0; continue;
            case kLess       : pos--; tab[pos-1] = (tab[pos-1] <  tab[pos]) ? 1 : 0; continue;
            case kGreater    : pos--; tab[pos-1] = (tab[pos-1] >  tab[pos]) ? 1 : 0; continue;
            case kLessThan   : pos--; tab[pos-1] = (tab[pos-1] <= tab[pos]) ? 1 : 0; continue;
            case kGreaterThan: pos--; tab[pos-1] = (tab[pos-1] >= tab[pos]) ? 1 : 0; continue;
            case kNot        :        tab[pos-1] = (tab[pos-1] !=        0) ? 0 : 1; continue;

            case kStringEqual : pos2 -= 2; pos++; if (!strcmp(stringStack[pos2+1],stringStack[pos2])) tab[pos-1]=1;
                                                  else tab[pos-1]=0;
                                continue;
            case kStringNotEqual: pos2 -= 2; pos++;if (strcmp(stringStack[pos2+1],stringStack[pos2])) tab[pos-1]=1;
                                                   else tab[pos-1]=0;
                                  continue;

            case kBitAnd    : pos--; tab[pos-1]= ((ULong64_t) tab[pos-1]) & ((ULong64_t) tab[pos]); continue;
            case kBitOr     : pos--; tab[pos-1]= ((ULong64_t) tab[pos-1]) | ((ULong64_t) tab[pos]); continue;
            case kLeftShift : pos--; tab[pos-1]= ((ULong64_t) tab[pos-1]) <<((ULong64_t) tab[pos]); continue;
            case kRightShift: pos--; tab[pos-1]= ((ULong64_t) tab[pos-1]) >>((ULong64_t) tab[pos]); continue;

            case kJump   : i = (oper & kTFOperMask); continue;
            case kJumpIf : {
               pos--;
               if (!tab[pos]) {
                  i = (oper & kTFOperMask);
                  // If we skip the left (true) side of the if statement we may,
                  // skip some of the branch loading (since we remove duplicate branch
                  // request (in TTreeFormula constructor) and so we need to force the
                  // loading here.
                  if (willLoad) fDidBooleanOptimization = kTRUE;
               }
               continue;
            }

            case kStringConst: {
               // String
               pos2++; stringStack[pos2-1] = (char*)fExpr[i].Data();
               if (fAxis) {
                  // See TT_EVAL_INIT
                  Int_t bin = fAxis->FindBin(stringStack[pos2-1]);
                  return bin;
               }
               continue;
            }

            case kBoolOptimize: {
               // boolean operation optimizer

               int param = (oper & kTFOperMask);
               Bool_t skip = kFALSE;
               int op = param % 10; // 1 is && , 2 is ||

               if (op == 1 && (!tab[pos-1]) ) {
                  // &&: skip the right part if the left part is already false

                  skip = kTRUE;

                  // Preserve the existing behavior (i.e. the result of a&&b is
                  // either 0 or 1)
                  tab[pos-1] = 0;

               } else if (op == 2 && tab[pos-1] ) {
                  // ||: skip the right part if the left part is already true

                  skip = kTRUE;

                  // Preserve the existing behavior (i.e. the result of a||b is
                  // either 0 or 1)
                  tab[pos-1] = 1;
               }

               if (skip) {
                  int toskip = param / 10;
                  i += toskip;
                  if (willLoad) fDidBooleanOptimization = kTRUE;
               }
               continue;
            }

            case kFunctionCall: {
               // an external function call

               int param = (oper & kTFOperMask);
               int fno   = param / 1000;
               int nargs = param % 1000;

               // Retrieve the function
               TMethodCall *method = (TMethodCall*)fFunctions.At(fno);

               // Set the arguments
               method->ResetParam();
               if (nargs) {
                  UInt_t argloc = pos-nargs;
                  for(Int_t j=0;j<nargs;j++,argloc++,pos--) {
                     SetMethodParam(method, tab[argloc]);
                  }
               }
               pos++;
               Double_t ret = 0;
               method->Execute(ret);
               tab[pos-1] = ret; // check for the correct conversion!

               continue;
            }

//         case kParameter:    { pos++; tab[pos-1] = fParams[(oper & kTFOperMask)]; continue; }
         }

      } else {
         // TTreeFormula operands.

         // a tree variable (the most used case).

         if (newaction == kDefinedVariable) {

            const Int_t code = (oper & kTFOperMask);
            const Int_t lookupType = fLookupType[code];
            switch (lookupType) {
               case kIndexOfEntry: tab[pos++] = (T)fTree->GetReadEntry(); continue;
               case kIndexOfLocalEntry: tab[pos++] = (T)fTree->GetTree()->GetReadEntry(); continue;
               case kEntries:      tab[pos++] = (T)fTree->GetEntries(); continue;
               case kLocalEntries: tab[pos++] = (T)fTree->GetTree()->GetEntries(); continue;
               case kLength:       tab[pos++] = fManager->fNdata; continue;
               case kLengthFunc:   tab[pos++] = ((TTreeFormula*)fAliases.UncheckedAt(i))->GetNdata(); continue;
               case kIteration:    tab[pos++] = instance; continue;
               case kSum:          tab[pos++] = Summing<T>((TTreeFormula*)fAliases.UncheckedAt(i)); continue;
               case kMin:          tab[pos++] = FindMin<T>((TTreeFormula*)fAliases.UncheckedAt(i)); continue;
               case kMax:          tab[pos++] = FindMax<T>((TTreeFormula*)fAliases.UncheckedAt(i)); continue;

               case kDirect:     { TT_EVAL_INIT_LOOP; tab[pos++] = leaf->GetTypedValue<T>(real_instance); continue; }
               case kMethod:     { TT_EVAL_INIT_LOOP; tab[pos++] = GetValueFromMethod(code,leaf); continue; }
               case kDataMember: { TT_EVAL_INIT_LOOP; tab[pos++] = ((TFormLeafInfo*)fDataMembers.UncheckedAt(code))->
                                          GetTypedValue<T>(leaf,real_instance); continue; }
               case kTreeMember: { TREE_EVAL_INIT_LOOP; tab[pos++] = ((TFormLeafInfo*)fDataMembers.UncheckedAt(code))->
                                          GetTypedValue<T>((TLeaf*)0x0,real_instance); continue; }
               case kEntryList: { TEntryList *elist = (TEntryList*)fExternalCuts.At(code);
                  tab[pos++] = elist->Contains(fTree->GetReadEntry());
                  continue;}
               case -1: break;
               default: tab[pos++] = 0; continue;
            }
            switch (fCodes[code]) {
               case -2: {
                  TCutG *gcut = (TCutG*)fExternalCuts.At(code);
                  TTreeFormula *fx = (TTreeFormula *)gcut->GetObjectX();
                  TTreeFormula *fy = (TTreeFormula *)gcut->GetObjectY();
                  if (fDidBooleanOptimization) {
                     fx->ResetLoading();
                     fy->ResetLoading();
                  }
                  T xcut = fx->EvalInstance<T>(instance);
                  T ycut = fy->EvalInstance<T>(instance);
                  tab[pos++] = gcut->IsInside(xcut,ycut);
                  continue;
               }
               case -1: {
                  TCutG *gcut = (TCutG*)fExternalCuts.At(code);
                  TTreeFormula *fx = (TTreeFormula *)gcut->GetObjectX();
                  if (fDidBooleanOptimization) {
                     fx->ResetLoading();
                  }
                  tab[pos++] = fx->EvalInstance<T>(instance);
                  continue;
               }
               default: {
                  tab[pos++] = 0;
                  continue;
               }
            }
         }
         switch(newaction) {

            // a TTree Variable Alias (i.e. a sub-TTreeFormula)
            case kAlias: {
               int aliasN = i;
               TTreeFormula *subform = static_cast<TTreeFormula*>(fAliases.UncheckedAt(aliasN));
               R__ASSERT(subform);

               subform->fDidBooleanOptimization = fDidBooleanOptimization;
               T param = subform->EvalInstance<T>(instance);

               tab[pos] = param; pos++;
               continue;
            }
            // a TTree Variable Alias String (i.e. a sub-TTreeFormula)
            case kAliasString: {
               int aliasN = i;
               TTreeFormula *subform = static_cast<TTreeFormula*>(fAliases.UncheckedAt(aliasN));
               R__ASSERT(subform);

               pos2++;
               subform->fDidBooleanOptimization = fDidBooleanOptimization;
               stringStack[pos2-1] = subform->EvalStringInstance(instance);
               continue;
            }
            case kMinIf: {
               int alternateN = i;
               TTreeFormula *primary = static_cast<TTreeFormula*>(fAliases.UncheckedAt(alternateN));
               TTreeFormula *condition = static_cast<TTreeFormula*>(fAliases.UncheckedAt(alternateN+1));
               T param = FindMin<T>(primary,condition);
               ++i; // skip the place holder for the condition
               tab[pos] = param; pos++;
               continue;
            }
            case kMaxIf: {
               int alternateN = i;
               TTreeFormula *primary = static_cast<TTreeFormula*>(fAliases.UncheckedAt(alternateN));
               TTreeFormula *condition = static_cast<TTreeFormula*>(fAliases.UncheckedAt(alternateN+1));
               T param = FindMax<T>(primary,condition);
               ++i; // skip the place holder for the condition
               tab[pos] = param; pos++;
               continue;
            }

            // a TTree Variable Alternate (i.e. a sub-TTreeFormula)
            case kAlternate: {
               int alternateN = i;
               TTreeFormula *primary = static_cast<TTreeFormula*>(fAliases.UncheckedAt(alternateN));

               // First check whether we are in range for the primary formula
               if (instance < primary->GetNdata()) {

                  T param = primary->EvalInstance<T>(instance);

                  ++i; // skip the alternate value.

                  tab[pos] = param; pos++;
               } else {
                  // The primary is not in range, we will calculate the alternate value
                  // via the next operation (which will be a intentional).

                  // kAlias no operations
               }
               continue;
            }
            case kAlternateString: {
               int alternateN = i;
               TTreeFormula *primary = static_cast<TTreeFormula*>(fAliases.UncheckedAt(alternateN));

               // First check whether we are in range for the primary formula
               if (instance < primary->GetNdata()) {

                  pos2++;
                  stringStack[pos2-1] = primary->EvalStringInstance(instance);

                  ++i; // skip the alternate value.

               } else {
                  // The primary is not in range, we will calculate the alternate value
                  // via the next operation (which will be a kAlias).

                  // intentional no operations
               }
               continue;
            }

            // a tree string
            case kDefinedString: {
               Int_t string_code = (oper & kTFOperMask);
               TLeaf *leafc = (TLeaf*)fLeaves.UncheckedAt(string_code);

               // Now let calculate what physical instance we really need.
               const Int_t real_instance = GetRealInstance(instance,string_code);

               if (instance==0 || fNeedLoading) {
                  fNeedLoading = kFALSE;
                  TBranch *branch = leafc->GetBranch();
                  Long64_t readentry = branch->GetTree()->GetReadEntry();
                  R__LoadBranch(branch,readentry,fQuickLoad);
               } else {
                  // In the cases where we are behind (i.e. right of) a potential boolean optimization
                  // this tree variable reading may have not been executed with instance==0 which would
                  // result in the branch being potentially not read in.
                  if (fDidBooleanOptimization) {
                     TBranch *br = leafc->GetBranch();
                     Long64_t treeEntry = br->GetTree()->GetReadEntry();
                     R__LoadBranch(br,treeEntry,kTRUE);
                  }
                  if (real_instance>=fNdata[string_code]) return 0;
               }
               pos2++;
               if (fLookupType[string_code]==kDirect) {
                  stringStack[pos2-1] = (char*)leafc->GetValuePointer();
               } else {
                  stringStack[pos2-1] = (char*)GetLeafInfo(string_code)->GetValuePointer(leafc,real_instance);
               }
               continue;
            }

         }
      }

      R__ASSERT(i<fNoper);
   }

   //std::cout << __PRETTY_FUNCTION__ << "  returning " << tab[0] << std::endl;
   return tab[0];
}

// Template instantiations
template double TTreeFormula::EvalInstance<double> (int, char const**);
template long double TTreeFormula::EvalInstance<long double> (int, char const**);
template long long TTreeFormula::EvalInstance<long long> (int, char const**);

////////////////////////////////////////////////////////////////////////////////
/// Return DataMember corresponding to code.
///
///  function called by TLeafObject::GetValue
///  with the value of fLookupType computed in TTreeFormula::DefinedVariable

TFormLeafInfo *TTreeFormula::GetLeafInfo(Int_t code) const
{
   return (TFormLeafInfo *)fDataMembers.UncheckedAt(code);

}

////////////////////////////////////////////////////////////////////////////////
/// Return leaf corresponding to serial number n.

TLeaf *TTreeFormula::GetLeaf(Int_t n) const
{
   return (TLeaf*)fLeaves.UncheckedAt(n);
}

////////////////////////////////////////////////////////////////////////////////
/// Return methodcall corresponding to code.
///
/// function called by TLeafObject::GetValue
/// with the value of fLookupType computed in TTreeFormula::DefinedVariable

TMethodCall *TTreeFormula::GetMethodCall(Int_t code) const
{
   return (TMethodCall *)fMethods.UncheckedAt(code);

}

////////////////////////////////////////////////////////////////////////////////
/// Return number of available instances in the formula.

Int_t TTreeFormula::GetNdata()
{
  return fManager->GetNdata();
}

////////////////////////////////////////////////////////////////////////////////
/// Return result of a leafobject method.

Double_t TTreeFormula::GetValueFromMethod(Int_t i, TLeaf* leaf) const
{
   TMethodCall* m = GetMethodCall(i);

   if (!m) {
      return 0.0;
   }

   void* thisobj = 0;
   if (leaf->InheritsFrom(TLeafObject::Class())) {
      thisobj = ((TLeafObject*) leaf)->GetObject();
   } else {
      TBranchElement* branch = (TBranchElement*) ((TLeafElement*) leaf)->GetBranch();
      Int_t id = branch->GetID();
      // FIXME: This is wrong for a top-level branch.
      Int_t offset = 0;
      if (id > -1) {
         TStreamerInfo* info = branch->GetInfo();
         if (info) {
            offset = info->GetElementOffset(id);
         } else {
            Warning("GetValueFromMethod", "No streamer info for branch %s.", branch->GetName());
         }
      }
      if (id < 0) {
         char* address = branch->GetObject();
         thisobj = address;
      } else {
         //char* address = branch->GetAddress();
         char* address = branch->GetObject();
         if (address) {
            thisobj = *((char**) (address + offset));
         } else {
            // FIXME: If the address is not set, the object won't be either!
            thisobj = branch->GetObject();
         }
      }
   }

   TMethodCall::EReturnType r = m->ReturnType();

   if (r == TMethodCall::kLong) {
      Longptr_t l = 0;
      m->Execute(thisobj, l);
      return (Double_t) l;
   }

   if (r == TMethodCall::kDouble) {
      Double_t d = 0.0;
      m->Execute(thisobj, d);
      return d;
   }

   m->Execute(thisobj);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return result of a leafobject method.

void* TTreeFormula::GetValuePointerFromMethod(Int_t i, TLeaf* leaf) const
{
   TMethodCall* m = GetMethodCall(i);

   if (!m) {
      return 0;
   }

   void* thisobj;
   if (leaf->InheritsFrom(TLeafObject::Class())) {
      thisobj = ((TLeafObject*) leaf)->GetObject();
   } else {
      TBranchElement* branch = (TBranchElement*) ((TLeafElement*) leaf)->GetBranch();
      Int_t id = branch->GetID();
      Int_t offset = 0;
      if (id > -1) {
         TStreamerInfo* info = branch->GetInfo();
         if (info) {
            offset = info->GetElementOffset(id);
         } else {
            Warning("GetValuePointerFromMethod", "No streamer info for branch %s.", branch->GetName());
         }
      }
      if (id < 0) {
         char* address = branch->GetObject();
         thisobj = address;
      } else {
         //char* address = branch->GetAddress();
         char* address = branch->GetObject();
         if (address) {
            thisobj = *((char**) (address + offset));
         } else {
            // FIXME: If the address is not set, the object won't be either!
            thisobj = branch->GetObject();
         }
      }
   }

   TMethodCall::EReturnType r = m->ReturnType();

   if (r == TMethodCall::kLong) {
      Longptr_t l = 0;
      m->Execute(thisobj, l);
      return 0;
   }

   if (r == TMethodCall::kDouble) {
      Double_t d = 0.0;
      m->Execute(thisobj, d);
      return 0;
   }

   if (r == TMethodCall::kOther) {
      char* c = 0;
      m->Execute(thisobj, &c);
      return c;
   }

   m->Execute(thisobj);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return TRUE if the formula corresponds to one single Tree leaf
/// and this leaf is short, int or unsigned short, int
/// When a leaf is of type integer or string, the generated histogram is forced
/// to have an integer bin width

Bool_t TTreeFormula::IsInteger(Bool_t fast) const
{
   if (fast) {
      if (TestBit(kIsInteger)) return kTRUE;
      else                     return kFALSE;
   }

   if (fNoper==2 && GetAction(0)==kAlternate) {
      TTreeFormula *subform = static_cast<TTreeFormula*>(fAliases.UncheckedAt(0));
      R__ASSERT(subform);
      return subform->IsInteger(kFALSE);
   }

   if (GetAction(0)==kMinIf || GetAction(0)==kMaxIf) {
      return kFALSE;
   }

   if (fNoper > 1) return kFALSE;

   if (GetAction(0)==kAlias) {
      TTreeFormula *subform = static_cast<TTreeFormula*>(fAliases.UncheckedAt(0));
      R__ASSERT(subform);
      return subform->IsInteger(kFALSE);
   }

   if (fLeaves.GetEntries() != 1) {
      switch (fLookupType[0]) {
         case kIndexOfEntry:
         case kIndexOfLocalEntry:
         case kEntries:
         case kLocalEntries:
         case kLength:
         case kLengthFunc:
         case kIteration:
           return kTRUE;
         case kSum:
         case kMin:
         case kMax:
         case kEntryList:
         default:
           return kFALSE;
      }
   }

   if (EvalClass()==TBits::Class()) return kTRUE;

   if (IsLeafInteger(0) || IsLeafString(0)) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return TRUE if the leaf corresponding to code is short, int or unsigned
/// short, int When a leaf is of type integer, the generated histogram is
/// forced to have an integer bin width

Bool_t TTreeFormula::IsLeafInteger(Int_t code) const
{
   TLeaf *leaf = (TLeaf*)fLeaves.At(code);
   if (!leaf) {
      switch (fLookupType[code]) {
         case kIndexOfEntry:
         case kIndexOfLocalEntry:
         case kEntries:
         case kLocalEntries:
         case kLength:
         case kLengthFunc:
         case kIteration:
            return kTRUE;
         case kSum:
         case kMin:
         case kMax:
         case kEntryList:
         default:
            return kFALSE;
      }
   }
   if (fAxis) return kTRUE;
   TFormLeafInfo * info;
   switch (fLookupType[code]) {
      case kMethod:
      case kTreeMember:
      case kDataMember:
         info = GetLeafInfo(code);
         return info->IsInteger();
      case kDirect:
         break;
   }
   if (!strcmp(leaf->GetTypeName(),"Int_t"))    return kTRUE;
   if (!strcmp(leaf->GetTypeName(),"Short_t"))  return kTRUE;
   if (!strcmp(leaf->GetTypeName(),"UInt_t"))   return kTRUE;
   if (!strcmp(leaf->GetTypeName(),"UShort_t")) return kTRUE;
   if (!strcmp(leaf->GetTypeName(),"Bool_t"))   return kTRUE;
   if (!strcmp(leaf->GetTypeName(),"Char_t"))   return kTRUE;
   if (!strcmp(leaf->GetTypeName(),"UChar_t"))  return kTRUE;
   if (!strcmp(leaf->GetTypeName(),"Long64_t"))  return kTRUE;
   if (!strcmp(leaf->GetTypeName(),"ULong64_t"))  return kTRUE;
   if (!strcmp(leaf->GetTypeName(),"string"))   return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return TRUE if the formula is a string

Bool_t TTreeFormula::IsString() const
{
   // See TTreeFormula::Init for the setting of kIsCharacter.
   return TestBit(kIsCharacter);
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the expression at the index 'oper' is to be treated as
/// as string.

Bool_t TTreeFormula::IsString(Int_t oper) const
{
   if (ROOT::v5::TFormula::IsString(oper)) return kTRUE;
   if (GetAction(oper)==kDefinedString) return kTRUE;
   if (GetAction(oper)==kAliasString) return kTRUE;
   if (GetAction(oper)==kAlternateString) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return TRUE if the leaf or data member corresponding to code is a string

Bool_t  TTreeFormula::IsLeafString(Int_t code) const
{
   TLeaf *leaf = (TLeaf*)fLeaves.At(code);
   TFormLeafInfo * info;
   if (fLookupType[code]==kTreeMember) {
      info = GetLeafInfo(code);
      return info->IsString();
   }

   switch(fLookupType[code]) {
      case kDirect:
         if ( !leaf->IsUnsigned() && (leaf->InheritsFrom(TLeafC::Class()) || leaf->InheritsFrom(TLeafB::Class()) ) ) {
            // Need to find out if it is an 'array' or a pointer.
            if (leaf->GetLenStatic() > 1) return kTRUE;

            // Now we need to differentiate between a variable length array and
            // a TClonesArray.
            if (leaf->GetLeafCount()) {
               const char* indexname = leaf->GetLeafCount()->GetName();
               if (indexname[strlen(indexname)-1] == '_' ) {
                  // This in a clones array
                  return kFALSE;
               } else {
                  // this is a variable length char array
                  return kTRUE;
               }
            }
            return kFALSE;
         } else if (leaf->InheritsFrom(TLeafElement::Class())) {
            TBranchElement * br = (TBranchElement*)leaf->GetBranch();
            Int_t bid = br->GetID();
            if (bid < 0) return kFALSE;
            if (br->GetInfo()==0 || !br->GetInfo()->IsCompiled()) {
               // Case where the file is corrupted is some ways.
               // We can not get to the actual type of the data
               // let's assume it is NOT a string.
               return kFALSE;
            }
            TStreamerElement * elem = (TStreamerElement*) br->GetInfo()->GetElement(bid);
            if (!elem) {
               // Case where the file is corrupted is some ways.
               // We can not get to the actual type of the data
               // let's assume it is NOT a string.
               return kFALSE;
            }
            if (elem->GetNewType()== TStreamerInfo::kOffsetL +kChar_t) {
               // Check whether a specific element of the string is specified!
               if (fIndexes[code][fNdimensions[code]-1] != -1) return kFALSE;
               return kTRUE;
            }
            if ( elem->GetNewType()== TStreamerInfo::kCharStar) {
               // Check whether a specific element of the string is specified!
               if (fNdimensions[code] && fIndexes[code][fNdimensions[code]-1] != -1) return kFALSE;
               return kTRUE;
            }
            return kFALSE;
         } else {
            return kFALSE;
         }
      case kMethod:
         //TMethodCall *m = GetMethodCall(code);
         //TMethodCall::EReturnType r = m->ReturnType();
         return kFALSE;
      case kDataMember:
         info = GetLeafInfo(code);
         return info->IsString();
      default:
         return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return value of variable as a string
///
/// -  mode = -2 : Print line with ***
/// -  mode = -1 : Print column names
/// -  mode = 0  : Print column values

char *TTreeFormula::PrintValue(Int_t mode) const
{
   return PrintValue(mode,0);
}

////////////////////////////////////////////////////////////////////////////////
/// Return value of variable as a string
///
/// -  mode = -2 : Print line with ***
/// -  mode = -1 : Print column names
/// -  mode = 0  : Print column values
///
///  decform contains the requested format (with the same convention as printf).

char *TTreeFormula::PrintValue(Int_t mode, Int_t instance, const char *decform) const
{
   const int kMAXLENGTH = 1024;
   static char value[kMAXLENGTH];

   if (mode == -2) {
      for (int i = 0; i < kMAXLENGTH-1; i++)
         value[i] = '*';
      value[kMAXLENGTH-1] = 0;
   } else if (mode == -1) {
      snprintf(value, kMAXLENGTH-1, "%s", GetTitle());
   } else if (mode == 0) {
      if ( (fNstring && fNval==0 && fNoper==1) || IsString() )
      {
         const char * val = 0;
         if (GetAction(0)==kStringConst) {
            val = fExpr[0].Data();
         } else if (instance<fNdata[0]) {
            if (fNoper == 1) {
               if (fLookupType[0]==kTreeMember) {
                  val = (char*)GetLeafInfo(0)->GetValuePointer((TLeaf*)0x0,instance);
               } else {
                  TLeaf *leaf = (TLeaf*)fLeaves.UncheckedAt(0);
                  TBranch *branch = leaf->GetBranch();
                  Long64_t readentry = branch->GetTree()->GetReadEntry();
                  R__LoadBranch(branch,readentry,fQuickLoad);
                  if (fLookupType[0]==kDirect && fNoper==1) {
                     val = (const char*)leaf->GetValuePointer();
                  } else {
                     val = ((TTreeFormula*)this)->EvalStringInstance(instance);
                  }
               }
            } else {
               val = ((TTreeFormula*)this)->EvalStringInstance(instance);
            }
         }
         if (val) {
            strlcpy(value, val, kMAXLENGTH);
         } else {
            value[0] = '\0';
         }
         value[kMAXLENGTH-1] = 0;
      } else {
         //NOTE: This is terrible form ... but is forced upon us by the fact that we can not
         //use the mutable keyword AND we should keep PrintValue const.
         Int_t real_instance = ((TTreeFormula*)this)->GetRealInstance(instance,-1);
         if (real_instance<fNdata[0]) {
            Ssiz_t len = strlen(decform);
            Char_t outputSizeLevel = 1;
            char *expo = 0;
            if (len>2) {
               switch (decform[len-2]) {
                  case 'l':
                  case 'L': {
                     outputSizeLevel = 2;
                     if (len>3 && tolower(decform[len-3])=='l') {
                        outputSizeLevel = 3;
                     }
                     break;
                  }
                  case 'h': outputSizeLevel = 0; break;
               }
            }
            switch(decform[len-1]) {
               case 'c':
               case 'd':
               case 'i':
               {
                  switch (outputSizeLevel) {
                     case 0:  snprintf(value,kMAXLENGTH,Form("%%%s",decform),(Short_t)((TTreeFormula*)this)->EvalInstance(instance)); break;
                     case 2:  snprintf(value,kMAXLENGTH,Form("%%%s",decform),(Long_t)((TTreeFormula*)this)->EvalInstance(instance)); break;
                     case 3:  snprintf(value,kMAXLENGTH,Form("%%%s",decform),(Long64_t)((TTreeFormula*)this)->EvalInstance<LongDouble_t>(instance)); break;
                     case 1:
                     default: snprintf(value,kMAXLENGTH,Form("%%%s",decform),(Int_t)((TTreeFormula*)this)->EvalInstance(instance)); break;
                  }
                  break;
               }
               case 'o':
               case 'x':
               case 'X':
               case 'u':
               {
                  switch (outputSizeLevel) {
                     case 0:  snprintf(value,kMAXLENGTH,Form("%%%s",decform),(UShort_t)((TTreeFormula*)this)->EvalInstance(instance)); break;
                     case 2:  snprintf(value,kMAXLENGTH,Form("%%%s",decform),(ULong_t)((TTreeFormula*)this)->EvalInstance(instance)); break;
                     case 3:  snprintf(value,kMAXLENGTH,Form("%%%s",decform),(ULong64_t)((TTreeFormula*)this)->EvalInstance<LongDouble_t>(instance)); break;
                     case 1:
                     default: snprintf(value,kMAXLENGTH,Form("%%%s",decform),(UInt_t)((TTreeFormula*)this)->EvalInstance(instance)); break;
                  }
                  break;
               }
               case 'f':
               case 'e':
               case 'E':
               case 'g':
               case 'G':
               {
                  switch (outputSizeLevel) {
                     case 2:  snprintf(value,kMAXLENGTH,Form("%%%s",decform),(LongDouble_t)((TTreeFormula*)this)->EvalInstance<LongDouble_t>(instance)); break;
                     case 1:
                     default: snprintf(value,kMAXLENGTH,Form("%%%s",decform),((TTreeFormula*)this)->EvalInstance(instance)); break;
                  }
                  expo = strchr(value,'e');
                  break;
               }
               default:
                  snprintf(value,kMAXLENGTH,Form("%%%sg",decform),((TTreeFormula*)this)->EvalInstance(instance));
                  expo = strchr(value,'e');
            }
            if (expo) {
               // If there is an exponent we may be longer than planned.
               // so let's trim off the excess precision!
               UInt_t declen = atoi(decform);
               if (strlen(value)>declen) {
                  UInt_t off = strlen(value)-declen;
                  char *start = expo - off;
                  UInt_t vlen = strlen(expo);
                  for(UInt_t z=0;z<=vlen;++z) {
                     start[z] = expo[z];
                  }
                  //strcpy(expo-off,expo);
               }
            }
         } else {
            if (isalpha(decform[strlen(decform)-1])) {
               TString short_decform(decform);
               short_decform.Remove(short_decform.Length()-1);
               snprintf(value,kMAXLENGTH,Form(" %%%sc",short_decform.Data()),' ');
            } else {
               snprintf(value,kMAXLENGTH,Form(" %%%sc",decform),' ');
            }

         }
      }
   }
   return &value[0];
}

////////////////////////////////////////////////////////////////////////////////
/// Tell the formula that we are going to request a new entry.

void TTreeFormula::ResetLoading()
{
   fNeedLoading = kTRUE;
   fDidBooleanOptimization = kFALSE;

   for(Int_t i=0; i<fNcodes; ++i) {
      UInt_t max_dim = fNdimensions[i];
      for(UInt_t dim=0; dim<max_dim ;++dim) {
         if (fVarIndexes[i][dim]) {
            fVarIndexes[i][dim]->ResetLoading();
         }
      }
   }
   Int_t n = fAliases.GetLast();
   if ( fNoper < n ) {
      n = fNoper;
   }
   for(Int_t k=0; k <= n; ++k) {
      TTreeFormula *f = static_cast<TTreeFormula*>(fAliases.UncheckedAt(k));
      if (f) {
         f->ResetLoading();
      }
   }
   for (int i=0; i<fExternalCuts.GetSize(); i++) {
      auto c = dynamic_cast<TCutG*>(fExternalCuts.At(i));
      if (c) {
         ((TTreeFormula *)(c->GetObjectX()))->ResetLoading();
         ((TTreeFormula *)(c->GetObjectY()))->ResetLoading();
      }
   }
   fRealInstanceCache.fInstanceCache = 0;
   fRealInstanceCache.fLocalIndexCache = 0;
   fRealInstanceCache.fVirtAccumCache = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the axis (in particular get the type).

void TTreeFormula::SetAxis(TAxis *axis)
{
   if (!axis) {fAxis = 0; return;}
   if (IsString()) {
      fAxis = axis;
      if (fNoper==1 && GetAction(0)==kAliasString){
         TTreeFormula *subform = static_cast<TTreeFormula*>(fAliases.UncheckedAt(0));
         R__ASSERT(subform);
         subform->SetAxis(axis);
      } else if (fNoper==2 && GetAction(0)==kAlternateString){
         TTreeFormula *subform = static_cast<TTreeFormula*>(fAliases.UncheckedAt(0));
         R__ASSERT(subform);
         subform->SetAxis(axis);
      }
      // Since the bin are corresponding to 'string', we currently must also set
      // the axis to align the bins exactly on integer boundaries.
      axis->SetBit(TAxis::kIsInteger);
   } else if (IsInteger()) {
      axis->SetBit(TAxis::kIsInteger);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TTreeFormula.

void TTreeFormula::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TTreeFormula::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      ROOT::v5::TFormula::Streamer(R__b);
      R__b >> fTree;
      R__b >> fNcodes;
      R__b.ReadFastArray(fCodes, fNcodes);
      R__b >> fMultiplicity;
      Int_t instance;
      R__b >> instance; //data member removed
      R__b >> fNindex;
      if (fNindex) {
         fLookupType = new Int_t[fNindex];
         R__b.ReadFastArray(fLookupType, fNindex);
      }
      fMethods.Streamer(R__b);
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TTreeFormula::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Try to 'demote' a string into an array bytes.  If this is not possible,
/// return false.

Bool_t TTreeFormula::StringToNumber(Int_t oper)
{
   Int_t code = GetActionParam(oper);
   if (GetAction(oper)==kDefinedString && fLookupType[code]==kDirect) {
      if (oper>0 && GetAction(oper-1)==kJump) {
         // We are the second hand of a ternary operator, let's not do the fixing.
         return kFALSE;
      }
      TLeaf *leaf = (TLeaf*)fLeaves.At(code);
      if (leaf &&  (leaf->InheritsFrom(TLeafC::Class()) || leaf->InheritsFrom(TLeafB::Class()) ) ) {
         SetAction(oper, kDefinedVariable, code );
         fNval++;
         fNstring--;
         return kTRUE;
      }
   }
   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// This function is called TTreePlayer::UpdateFormulaLeaves, itself
/// called by TChain::LoadTree when a new Tree is loaded.
/// Because Trees in a TChain may have a different list of leaves, one
/// must update the leaves numbers in the TTreeFormula used by the TreePlayer.
///
/// A safer alternative would be to recompile the whole thing .... However
/// currently compile HAS TO be called from the constructor!

void TTreeFormula::UpdateFormulaLeaves()
{
   Int_t nleaves = fLeafNames.GetEntriesFast();
   ResetBit( kMissingLeaf );
   for (Int_t i=0;i<nleaves;i++) {
      if (!fTree) break;
      if (!fLeafNames[i]) continue;

      TLeaf *leaf = fTree->GetLeaf(fLeafNames[i]->GetTitle(),fLeafNames[i]->GetName());
      fLeaves[i] = leaf;
      if (fBranches[i] && leaf) {
         fBranches[i] = leaf->GetBranch();
         // Since sometimes we might no read all the branches for all the entries, we
         // might sometimes only read the branch count and thus reset the collection
         // but might not read the data branches, to insure that a subsequent read
         // from TTreeFormula will properly load the data branches even if fQuickLoad is true,
         // we reset the entry of all branches in the TTree.
         ((TBranch*)fBranches[i])->ResetReadEntry();
      }
      if (leaf==0) SetBit( kMissingLeaf );
   }
   for (Int_t j=0; j<kMAXCODES; j++) {
      for (Int_t k = 0; k<kMAXFORMDIM; k++) {
         if (fVarIndexes[j][k]) {
            fVarIndexes[j][k]->UpdateFormulaLeaves();
         }
      }
      if (fLookupType[j]==kDataMember || fLookupType[j]==kTreeMember) GetLeafInfo(j)->Update();
      if (j<fNval && fCodes[j]<0) {
         TCutG *gcut = (TCutG*)fExternalCuts.At(j);
         if (gcut) {
            TTreeFormula *fx = (TTreeFormula *)gcut->GetObjectX();
            TTreeFormula *fy = (TTreeFormula *)gcut->GetObjectY();
            if (fx) fx->UpdateFormulaLeaves();
            if (fy) fy->UpdateFormulaLeaves();
         }
      }
   }
   for(Int_t k=0;k<fNoper;k++) {
      const Int_t oper = GetOper()[k];
      switch(oper >> kTFOperShift) {
         case kAlias:
         case kAliasString:
         case kAlternate:
         case kAlternateString:
         case kMinIf:
         case kMaxIf:
         {
            TTreeFormula *subform = static_cast<TTreeFormula*>(fAliases.UncheckedAt(k));
            R__ASSERT(subform);
            subform->UpdateFormulaLeaves();
            break;
         }
         case kDefinedVariable:
         {
            Int_t code = GetActionParam(k);
            if (fCodes[code]==0) switch(fLookupType[code]) {
               case kLengthFunc:
               case kSum:
               case kMin:
               case kMax:
               {
                  TTreeFormula *subform = static_cast<TTreeFormula*>(fAliases.UncheckedAt(k));
                  R__ASSERT(subform);
                  subform->UpdateFormulaLeaves();
                  break;
               }
               default:
                  break;
            }
         }
         default:
            break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Populate the TTreeFormulaManager with the dimension information.

void TTreeFormula::ResetDimensions() {
   Int_t i,k;

   // Now that we saw all the expressions and variables AND that
   // we know whether arrays of chars are treated as string or
   // not, we can properly setup the dimensions.
   TIter next(fDimensionSetup);
   Int_t last_code = -1;
   Int_t virt_dim = 0;
   for(TDimensionInfo * info; (info = (TDimensionInfo*)next()); ) {
      if (last_code!=info->fCode) {
         // We know that the list is ordered by code number then by
         // dimension.  Thus a different code means that we need to
         // restart at the lowest dimensions.
         virt_dim = 0;
         last_code = info->fCode;
         fNdimensions[last_code] = 0;
      }

      if (GetAction(info->fOper)==kDefinedString) {

         // We have a string used as a string (and not an array of number)
         // We need to determine which is the last dimension and skip it.
         TDimensionInfo *nextinfo = (TDimensionInfo*)next();
         while(nextinfo && nextinfo->fCode==info->fCode) {
            DefineDimensions(info->fCode,info->fSize, info->fMultiDim, virt_dim);
            nextinfo = (TDimensionInfo*)next();
         }
         if (!nextinfo) break;

         info = nextinfo;
         virt_dim = 0;
         last_code = info->fCode;
         fNdimensions[last_code] = 0;

         info->fSize = 1; // Maybe this should actually do nothing!
      }


      DefineDimensions(info->fCode,info->fSize, info->fMultiDim, virt_dim);
   }

   fMultiplicity = 0;
   for(i=0;i<fNoper;i++) {
      Int_t action = GetAction(i);

      if (action==kMinIf || action==kMaxIf) {
         // Skip/Ignore the 2nd args
         ++i;
         continue;
      }
      if (action==kAlias || action==kAliasString) {
         TTreeFormula *subform = static_cast<TTreeFormula*>(fAliases.UncheckedAt(i));
         R__ASSERT(subform);
         switch(subform->GetMultiplicity()) {
            case 0: break;
            case 1: fMultiplicity = 1; break;
            case 2: if (fMultiplicity!=1) fMultiplicity = 2; break;
         }
         fManager->Add(subform);
         // since we are addint to this manager 'subform->ResetDimensions();'
         // will be called a little latter
         continue;
      }
      if (action==kDefinedString) {
      //if (fOper[i] >= 105000 && fOper[i]<110000) {
         // We have a string used as a string

         // This dormant portion of code would be used if (when?) we allow the histogramming
         // of the integral content (as opposed to the string content) of strings
         // held in a variable size container delimited by a null (as opposed to
         // a fixed size container or variable size container whose size is controlled
         // by a variable).  In GetNdata, we will then use strlen to grab the current length.
         //fCumulSizes[i][fNdimensions[i]-1] = 1;
         //fUsedSizes[fNdimensions[i]-1] = -TMath::Abs(fUsedSizes[fNdimensions[i]-1]);
         //fUsedSizes[0] = - TMath::Abs( fUsedSizes[0]);

         //continue;
      }
   }

   for (i=0;i<fNcodes;i++) {
      if (fCodes[i] < 0) {
         TCutG *gcut = (TCutG*)fExternalCuts.At(i);
         if (!gcut) continue;
         TTreeFormula *fx = (TTreeFormula *)gcut->GetObjectX();
         TTreeFormula *fy = (TTreeFormula *)gcut->GetObjectY();

         if (fx) {
            switch(fx->GetMultiplicity()) {
               case 0: break;
               case 1: fMultiplicity = 1; break;
               case 2: if (fMultiplicity!=1) fMultiplicity = 2; break;
            }
            fManager->Add(fx);
         }
         if (fy) {
            switch(fy->GetMultiplicity()) {
               case 0: break;
               case 1: fMultiplicity = 1; break;
               case 2: if (fMultiplicity!=1) fMultiplicity = 2; break;
            }
            fManager->Add(fy);
         }

         continue;
      }

      if (fLookupType[i]==kIteration) {
         fMultiplicity = 1;
         continue;
      }

      TLeaf *leaf = i <= fLeaves.GetLast() ? (TLeaf*)fLeaves.UncheckedAt(i) : 0;
      if (!leaf) continue;

      // Reminder of the meaning of fMultiplicity:
      //  -1: Only one or 0 element per entry but contains variable length
      //      -array! (Only used for TTreeFormulaManager)
      //   0: Only one element per entry, no variable length array
      //   1: loop over the elements of a variable length array
      //   2: loop over elements of fixed length array (nData is the same for all entry)

      if (leaf->GetLeafCount()) {
         // We assume only one possible variable length dimension (the left most)
         fMultiplicity = 1;
      } else if (fLookupType[i]==kDataMember) {
         TFormLeafInfo * leafinfo = GetLeafInfo(i);
         TStreamerElement * elem = leafinfo->fElement;
         if (fMultiplicity!=1) {
            if (leafinfo->HasCounter() ) fMultiplicity = 1;
            else if (elem && elem->GetArrayDim()>0) fMultiplicity = 2;
            else if (leaf->GetLenStatic()>1) fMultiplicity = 2;
         }
      } else {
         if (leaf->GetLenStatic()>1 && fMultiplicity!=1) fMultiplicity = 2;
      }
      if (fMultiplicity!=1) {
         // If the leaf belongs to a friend tree which has an index, we might
         // be in the case where some entry do not exist.

         TTree *realtree = fTree ? fTree->GetTree() : 0;
         TTree *tleaf = leaf->GetBranch()->GetTree();
         if (tleaf && tleaf != realtree && tleaf->GetTreeIndex()) {
            // Reset the multiplicity if we have a friend tree with an index.
            fMultiplicity = 1;
         }
      }

      Int_t virt_dim2 = 0;
      for (k = 0; k < fNdimensions[i]; k++) {
         // At this point fCumulSizes[i][k] actually contain the physical
         // dimension of the k-th dimensions.
         if ( (fCumulSizes[i][k]>=0) && (fIndexes[i][k] >= fCumulSizes[i][k]) ) {
            // unreachable element requested:
            fManager->CancelDimension(virt_dim2); // fCumulUsedSizes[virt_dim2] = 0;
         }
         if ( fIndexes[i][k] < 0 ) virt_dim2++;
         fFixedSizes[i][k] = fCumulSizes[i][k];
      }

      // Add up the cumulative size
      for (k = fNdimensions[i]; (k > 0); k--) {
         // NOTE: When support for inside variable dimension is added this
         // will become inaccurate (since one of the value in the middle of the chain
         // is unknown until GetNdata is called.
         fCumulSizes[i][k-1] *= TMath::Abs(fCumulSizes[i][k]);
      }
      // NOTE: We assume that the inside variable dimensions are dictated by the
      // first index.
      if (fCumulSizes[i][0]>0) fNdata[i] = fCumulSizes[i][0];

      //for (k = 0; k<kMAXFORMDIM; k++) {
      //   if (fVarIndexes[i][k]) fManager->Add(fVarIndexes[i][k]);
      //}

   }
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that all the branches have been loaded properly.

void TTreeFormula::LoadBranches()
{
   Int_t i;
   for (i=0; i<fNoper ; ++i) {
      TLeaf *leaf = (TLeaf*)fLeaves.UncheckedAt(i);
      if (leaf==0) continue;

      TBranch *br = leaf->GetBranch();
      Long64_t treeEntry = br->GetTree()->GetReadEntry();
      R__LoadBranch(br,treeEntry,kTRUE);

      TTreeFormula *alias = (TTreeFormula*)fAliases.UncheckedAt(i);
      if (alias) alias->LoadBranches();

      Int_t max_dim = fNdimensions[i];
      for (Int_t dim = 0; dim < max_dim; ++dim) {
         if (fVarIndexes[i][dim]) fVarIndexes[i][dim]->LoadBranches();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the actual dimension for the current entry.

Bool_t TTreeFormula::LoadCurrentDim() {
   Int_t size;
   Bool_t outofbounds = kFALSE;

   for (Int_t i=0;i<fNcodes;i++) {
      if (fCodes[i] < 0) continue;

      // NOTE: Currently only the leafcount can indicates a dimension that
      // is physically variable.  So only the left-most dimension is variable.
      // When an API is introduced to be able to determine a variable inside dimensions
      // one would need to add a way to recalculate the values of fCumulSizes for this
      // leaf.  This would probably require the addition of a new data member
      // fSizes[kMAXCODES][kMAXFORMDIM];
      // Also note that EvalInstance expect all the values (but the very first one)
      // of fCumulSizes to be positive.  So indicating that a physical dimension is
      // variable (expected for the first one) can NOT be done via negative values of
      // fCumulSizes.

      TLeaf *leaf = (TLeaf*)fLeaves.UncheckedAt(i);
      if (!leaf) {
         switch(fLookupType[i]) {
            case kDirect:
            case kMethod:
            case kTreeMember:
            case kDataMember:
               fNdata[i] = 0;
               outofbounds = kTRUE;
         }
         continue;
      }

      TTree *realtree = fTree->GetTree();
      TTree *tleaf = leaf->GetBranch()->GetTree();
      if (tleaf && tleaf != realtree && tleaf->GetTreeIndex()) {
         if (tleaf->GetReadEntry() < 0) {
            fNdata[i] = 0;
            outofbounds = kTRUE;
            continue;
         } else {
            fNdata[i] = fCumulSizes[i][0];
         }
      }
      Bool_t hasBranchCount2 = kFALSE;
      if (leaf->GetLeafCount()) {
         TLeaf* leafcount = leaf->GetLeafCount();
         TBranch *branchcount = leafcount->GetBranch();
         TFormLeafInfo * info = 0;
         if (leaf->IsA() == TLeafElement::Class()) {
            //if branchcount address not yet set, GetEntry will set the address
            // read branchcount value
            Long64_t readentry = leaf->GetBranch()->GetTree()->GetReadEntry();
            if (readentry < 0) readentry=0;
            if (!branchcount->GetAddress()) {
               R__LoadBranch(branchcount, readentry, fQuickLoad);
            } else {
               // Since we do not read the full branch let's reset the read entry number
               // so that a subsequent read from TTreeFormula will properly load the full
               // object even if fQuickLoad is true.
               branchcount->TBranch::GetEntry(readentry);
               branchcount->ResetReadEntry();
            }

            size = ((TBranchElement*)branchcount)->GetNdata();
            // Reading the size as above is correct only when the branchcount
            // is of streamer type kCounter which require the underlying data
            // member to be signed integral type.

            TBranchElement* branch = (TBranchElement*) leaf->GetBranch();
            if (branch->GetAddress() == 0) {
               // Humm there is no space reserve to write the data,
               // the data member is likely 'removed' from the class
               // layout, so rather than crashing by accessing
               // random memory, make it clear we can't read it.
               size = 0;
            }

            // NOTE: could be sped up
            if (fHasMultipleVarDim[i]) {// info && info->GetVarDim()>=0) {
               info = (TFormLeafInfo* )fDataMembers.At(i);
               if (branch->GetBranchCount2()) R__LoadBranch(branch->GetBranchCount2(),readentry,fQuickLoad);
               else R__LoadBranch(branch,readentry,fQuickLoad);

               // Here we need to add the code to take in consideration the
               // double variable length
               // We fill up the array of sizes in the TLeafInfo:
               info->LoadSizes(branch);
               hasBranchCount2 = kTRUE;
               if (info->GetVirtVarDim()>=0) info->UpdateSizes(fManager->fVarDims[info->GetVirtVarDim()]);

               // Refresh the fCumulSizes[i] to have '1' for the
               // double variable dimensions
               Int_t vdim = info->GetVarDim();
               fCumulSizes[i][vdim] =  fCumulSizes[i][vdim+1];
               for(Int_t k=vdim -1; k>=0; k--) {
                  fCumulSizes[i][k] = fCumulSizes[i][k+1]*fFixedSizes[i][k];
               }
               // Update fCumulUsedSizes
               // UpdateMultiVarSizes(vdim,info,i)
               //Int_t fixed = fCumulSizes[i][vdim+1];
               //for(Int_t k=vdim - 1; k>=0; k++) {
               //   Int_t fixed *= fFixedSizes[i][k];
               //   for(Int_t l=0;l<size; l++) {
               //     fCumulSizes[i][k] += info->GetSize(l) * fixed;
               //}
            }
         } else {
            Long64_t readentry = leaf->GetBranch()->GetTree()->GetReadEntry();
            if (readentry < 0) readentry=0;
            R__LoadBranch(branchcount,readentry,fQuickLoad);
            size = leaf->GetLen() / leaf->GetLenStatic();
         }
         if (hasBranchCount2) {
            // We assume that fCumulSizes[i][1] contains the product of the fixed sizes
            fNdata[i] = fCumulSizes[i][1] * ((TFormLeafInfo *)fDataMembers.At(i))->GetSumOfSizes();
         } else {
            fNdata[i] = size * fCumulSizes[i][1];
         }
         if (fIndexes[i][0]==-1) {
            // Case where the index is not specified AND the 1st dimension has a variable
            // size.
            if (fManager->fUsedSizes[0]==1 || (size<fManager->fUsedSizes[0]) ) fManager->fUsedSizes[0] = size;
            if (info && fIndexes[i][info->GetVarDim()]>=0) {
               for(Int_t j=0; j<size; j++) {
                  if (fIndexes[i][info->GetVarDim()] >= info->GetSize(j)) {
                     info->SetSize(j,0);
                     if (size>fManager->fCumulUsedVarDims->GetSize()) fManager->fCumulUsedVarDims->Set(size);
                     fManager->fCumulUsedVarDims->AddAt(-1,j);
                  } else if (fIndexes[i][info->GetVarDim()]>=0) {
                     // There is an index and it is not too large
                     info->SetSize(j,1);
                     if (size>fManager->fCumulUsedVarDims->GetSize()) fManager->fCumulUsedVarDims->Set(size);
                     fManager->fCumulUsedVarDims->AddAt(1,j);
                  }
               }
            }
         } else if (fIndexes[i][0] >= size) {
            // unreachable element requested:
            fManager->fUsedSizes[0] = 0;
            fNdata[i] = 0;
            outofbounds = kTRUE;
         } else if (hasBranchCount2) {
            TFormLeafInfo *info2;
            info2 = (TFormLeafInfo *)fDataMembers.At(i);
            if (fIndexes[i][0]<0
                || fIndexes[i][info2->GetVarDim()] >= info2->GetSize(fIndexes[i][0])) {
               // unreachable element requested:
               fManager->fUsedSizes[0] = 0;
               fNdata[i] = 0;
               outofbounds = kTRUE;
            }
         }
      } else if (fLookupType[i]==kDataMember) {
         TFormLeafInfo *leafinfo = (TFormLeafInfo*)fDataMembers.UncheckedAt(i);
         if (leafinfo->HasCounter()) {
            TBranch *branch = leaf->GetBranch();
            Long64_t readentry = branch->GetTree()->GetReadEntry();
            if (readentry < 0) readentry=0;
            R__LoadBranch(branch,readentry,fQuickLoad);
            size = (Int_t) leafinfo->GetCounterValue(leaf);
            if (fIndexes[i][0]==-1) {
               // Case where the index is not specified AND the 1st dimension has a variable
               // size.
               if (fManager->fUsedSizes[0]==1 || (size<fManager->fUsedSizes[0]) ) {
                  fManager->fUsedSizes[0] = size;
               }
            } else if (fIndexes[i][0] >= size) {
               // unreachable element requested:
               fManager->fUsedSizes[0] = 0;
               fNdata[i] = 0;
               outofbounds = kTRUE;
            } else {
               fNdata[i] = size*fCumulSizes[i][1];
            }
            Int_t vdim = leafinfo->GetVarDim();
            if (vdim>=0) {
               // Here we need to add the code to take in consideration the
               // double variable length
               // We fill up the array of sizes in the TLeafInfo:
               // here we can assume that branch is a TBranch element because the other style does NOT support this type
               // of complexity.
               leafinfo->LoadSizes(branch);
               hasBranchCount2 = kTRUE;
               if (fIndexes[i][0]==-1&&fIndexes[i][vdim] >= 0) {
                  for(int z=0; z<size; ++z) {
                     if (fIndexes[i][vdim] >= leafinfo->GetSize(z)) {
                        leafinfo->SetSize(z,0);
                        // --fManager->fUsedSizes[0];
                     } else if (fIndexes[i][vdim] >= 0 ) {
                        leafinfo->SetSize(z,1);
                     }
                  }
               }
               leafinfo->UpdateSizes(fManager->fVarDims[vdim]);

               // Refresh the fCumulSizes[i] to have '1' for the
               // double variable dimensions
               fCumulSizes[i][vdim] =  fCumulSizes[i][vdim+1];
               for(Int_t k=vdim -1; k>=0; k--) {
                  fCumulSizes[i][k] = fCumulSizes[i][k+1]*fFixedSizes[i][k];
               }
               fNdata[i] = fCumulSizes[i][1] * leafinfo->GetSumOfSizes();
            } else {
               fNdata[i] = size * fCumulSizes[i][1];
            }
         } else if (leafinfo->GetMultiplicity()==-1) {
            TBranch *branch = leaf->GetBranch();
            Long64_t readentry = branch->GetTree()->GetReadEntry();
            if (readentry < 0) readentry=0;
            R__LoadBranch(branch,readentry,fQuickLoad);
            if (leafinfo->GetNdata(leaf)==0) {
               outofbounds = kTRUE;
            }
         }
      }
      // However we allow several dimensions that virtually vary via the size of their
      // index variables.  So we have code to recalculate fCumulUsedSizes.
      TFormLeafInfo * info = 0;
      if (fLookupType[i]!=kDirect) {
         info = (TFormLeafInfo *)fDataMembers.At(i);
      }
      for(Int_t k=0, virt_dim=0; k < fNdimensions[i]; k++) {
         if (fIndexes[i][k]<0) {
            if (info && fIndexes[i][k]==-2 && fVarIndexes[i][k]->GetManager()->GetMultiplicity()==0) {
               // Index and thus local size provided by a "index variable of size 1"
               Int_t index = fVarIndexes[i][k]->EvalInstance(0);
               Int_t index_size = info->GetSize(index);
               if (fManager->fUsedSizes[virt_dim]==1 || (index_size!=1 && index_size<fManager->fUsedSizes[virt_dim]) )
                  fManager->fUsedSizes[virt_dim] = index_size;
            } else if (fIndexes[i][k]==-2 && fManager->fVirtUsedSizes[virt_dim]<0) {

               // if fVirtUsedSize[virt_dim] is positive then VarIndexes[i][k]->GetNdata()
               // is always the same and has already been factored in fUsedSize[virt_dim]
               Int_t index_size = fVarIndexes[i][k]->GetNdata();
               if (index_size==1) {
                  // We could either have a variable size array which is currently of size one
                  // or a single element that might or not might not be present (and is currently present!)
                  if (fVarIndexes[i][k]->GetManager()->GetMultiplicity()==1) {
                     if (index_size<fManager->fUsedSizes[virt_dim]) fManager->fUsedSizes[virt_dim] = index_size;
                  }

               } else if (fManager->fUsedSizes[virt_dim]==-fManager->fVirtUsedSizes[virt_dim] ||
                          index_size<fManager->fUsedSizes[virt_dim]) {
                  fManager->fUsedSizes[virt_dim] = index_size;
               }

            } else if (hasBranchCount2 && info && k==info->GetVarDim()) {
               // NOTE: We assume the indexing of variable sizes on the first index!
               if (fIndexes[i][0]>=0) {
                  Int_t index_size = info->GetSize(fIndexes[i][0]);
                  if (fManager->fUsedSizes[virt_dim]==1 || (index_size!=1 && index_size<fManager->fUsedSizes[virt_dim]) )
                     fManager->fUsedSizes[virt_dim] = index_size;
               }
            }
            virt_dim++;
         } else if (hasBranchCount2 && info && k==info->GetVarDim()) {

            // nothing to do, at some point I thought this might be useful:
            // if (fIndexes[i][k]>=0) {
            //    index = info->GetSize(fIndexes[i][k]);
            //    if (fManager->fUsedSizes[virt_dim]==1 || (index!=1 && index<fManager->fUsedSizes[virt_dim]) )
            //      fManager->fUsedSizes[virt_dim] = index;
            //    virt_dim++;
            // }

         }
      }
   }
   return ! outofbounds;



}

void TTreeFormula::Convert(UInt_t oldversion)
{
   // Convert the fOper of a TTTreeFormula version fromVersion to the current in memory version

   enum { kOldAlias           = /*ROOT::v5::TFormula::kVariable*/ 100000+10000+1,
          kOldAliasString     = kOldAlias+1,
          kOldAlternate       = kOldAlias+2,
          kOldAlternateString = kOldAliasString+2
   };

   for (int k=0; k<fNoper; k++) {
      // First hide from ROOT::v5::TFormula convertion

      Int_t action = GetOper()[k];

      switch (action) {

         case kOldAlias:            GetOper()[k] = -kOldAlias; break;
         case kOldAliasString:      GetOper()[k] = -kOldAliasString; break;
         case kOldAlternate:        GetOper()[k] = -kOldAlternate; break;
         case kOldAlternateString:  GetOper()[k] = -kOldAlternateString; break;
      }
   }

   ROOT::v5::TFormula::Convert(oldversion);

   for (int i=0,offset=0; i<fNoper; i++) {
      Int_t action = GetOper()[i+offset];

      switch (action) {
         case -kOldAlias:            SetAction(i, kAlias, 0); break;
         case -kOldAliasString:      SetAction(i, kAliasString, 0); break;
         case -kOldAlternate:        SetAction(i, kAlternate, 0); break;
         case -kOldAlternateString:  SetAction(i, kAlternateString, 0); break;
      }
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Convert the underlying lookup method from the direct technique
/// (dereferencing the address held by the branch) to the method using
/// TFormLeafInfo.  This is in particular useful in the case where we
/// need to append an additional TFormLeafInfo (for example to call a
/// method).
/// Return false if the switch was unsuccessful (basically in the
/// case of an old style split tree).

Bool_t TTreeFormula::SwitchToFormLeafInfo(Int_t code)
{
   TFormLeafInfo *last = 0;
   TLeaf *leaf = (TLeaf*)fLeaves.At(code);
   if (!leaf) return kFALSE;

   if (fLookupType[code]==kDirect) {
      if (leaf->InheritsFrom(TLeafElement::Class())) {
         TBranchElement * br = (TBranchElement*)leaf->GetBranch();
         if (br->GetType()==31) {
            // sub branch of a TClonesArray
            TStreamerInfo *info = br->GetInfo();
            TClass* cl = info->GetClass();
            TStreamerElement *element = (TStreamerElement *)info->GetElement(br->GetID());
            TFormLeafInfo* clonesinfo = new TFormLeafInfoClones(cl, 0, element, kTRUE);
            Int_t offset;
            info->GetStreamerElement(element->GetName(),offset);
            clonesinfo->fNext = new TFormLeafInfo(cl,offset+br->GetOffset(),element);
            last = clonesinfo->fNext;
            fDataMembers.AddAtAndExpand(clonesinfo,code);
            fLookupType[code]=kDataMember;

         } else if (br->GetType()==41) {
            // sub branch of a Collection

            TBranchElement *count = br->GetBranchCount();
            TFormLeafInfo* collectioninfo;
            if ( count->GetID() >= 0 ) {
               TStreamerElement *collectionElement =
                  (TStreamerElement *)count->GetInfo()->GetElement(count->GetID());
               TClass *collectionCl = collectionElement->GetClassPointer();

               collectioninfo =
                  new TFormLeafInfoCollection(collectionCl, 0, collectionElement, kTRUE);
            } else {
               TClass *collectionCl = TClass::GetClass(count->GetClassName());
               collectioninfo =
                  new TFormLeafInfoCollection(collectionCl, 0, collectionCl, kTRUE);
            }

            TStreamerInfo *info = br->GetInfo();
            TClass* cl = info->GetClass();
            TStreamerElement *element = (TStreamerElement *)info->GetElement(br->GetID());
            Int_t offset;
            info->GetStreamerElement(element->GetName(),offset);
            collectioninfo->fNext = new TFormLeafInfo(cl,offset+br->GetOffset(),element);
            last = collectioninfo->fNext;
            fDataMembers.AddAtAndExpand(collectioninfo,code);
            fLookupType[code]=kDataMember;

         } else if (br->GetID()<0) {
            return kFALSE;
         } else {
            last = new TFormLeafInfoDirect(br);
            fDataMembers.AddAtAndExpand(last,code);
            fLookupType[code]=kDataMember;
         }
      } else {
         //last = new TFormLeafInfoDirect(br);
         //fDataMembers.AddAtAndExpand(last,code);
         //fLookupType[code]=kDataMember;
         return kFALSE;
      }
   }
   return kTRUE;
}
