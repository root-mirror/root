// @(#)root/io:$Id$
// Author: Rene Brun   12/10/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**
\class TStreamerInfo TStreamerInfo.cxx
\ingroup IO

Describes a persistent version of a class.

A ROOT file contains the list of TStreamerInfo objects for all the
class versions written to this file.
When reading a file, all the TStreamerInfo objects are read back in
memory and registered to the TClass list of TStreamerInfo.
One can see the list and contents of the TStreamerInfo on a file
with, e.g.,
~~~{.cpp}
   TFile f("myfile.root");
   f.ShowStreamerInfo();
~~~
A TStreamerInfo is a list of TStreamerElement objects (one per data
member or base class).
When streaming an object, the system (TClass) loops on all the
TStreamerElement objects and calls the appropriate function for each
element type.
*/

#include "TStreamerInfo.h"
#include "TFile.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TStreamerElement.h"
#include "TClass.h"
#include "TClassEdit.h"
#include "TClassTable.h"
#include "TDataMember.h"
#include "TDataType.h"
#include "TRealData.h"
#include "TBaseClass.h"
#include "TBuffer.h"
#include "TArrayC.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TArrayD.h"
#include "TArrayS.h"
#include "TArrayL.h"
#include "TError.h"
#include "TEnum.h"
#include "TRef.h"
#include "TProcessID.h"
#include "TSystem.h"
#include "TObjString.h"
#include "snprintf.h"

#include "TStreamer.h"
#include "TContainerConverters.h"
#include "TCollectionProxyFactory.h"
#include "TVirtualCollectionProxy.h"
#include "TInterpreter.h"

#include "TMemberInspector.h"

#include "TMakeProject.h"

#include "TSchemaRuleSet.h"
#include "TSchemaRule.h"

#include "TVirtualMutex.h"

#include "TStreamerInfoActions.h"

#include <memory>
#include <array>

std::atomic<Int_t> TStreamerInfo::fgCount{0};

const Int_t kMaxLen = 1024;

ClassImp(TStreamerInfo);

static void R__TObjArray_InsertAt(TObjArray *arr, TObject *obj, Int_t at)
{
   // Slide by one.
   Int_t last = arr->GetLast();
   arr->AddAtAndExpand(arr->At(last),last+1);
   for(Int_t ind = last-1; ind >= at; --ind) {
      arr->AddAt( arr->At(ind), ind+1);
   };
   arr->AddAt( obj, at);
}

static void R__TObjArray_InsertAt(TObjArray *arr, std::vector<TStreamerArtificial*> &objs, Int_t at)
{
   // Slide by enough.
   Int_t offset = objs.size();
   Int_t last = arr->GetLast();
   arr->AddAtAndExpand(arr->At(last),last+offset);
   for(Int_t ind = last-1; ind >= at; --ind) {
      arr->AddAt( arr->At(ind), ind+offset);
   };
   for(size_t ins = 0; ins < objs.size(); ++ins) {
      arr->AddAt(objs[ins], at+ins);
   }
}

static void R__TObjArray_InsertAfter(TObjArray *arr, TObject *newobj, TObject *oldobj)
{
   // Slide by one.
   Int_t last = arr->GetLast();
   Int_t at = 0;
   while (at<last && arr->At(at) != oldobj) {
      ++at;
   }
   ++at; // we found the object, insert after it
   R__TObjArray_InsertAt(arr, newobj, at);
}

static void R__TObjArray_InsertBefore(TObjArray *arr, TObject *newobj, TObject *oldobj)
{
   // Slide by one.
   Int_t last = arr->GetLast();
   Int_t at = 0;
   while (at<last && arr->At(at) != oldobj) {
      ++at;
   }
   R__TObjArray_InsertAt(arr, newobj, at);
}

enum class EUniquePtrOffset : char
   {
      kNA = 0,
      kZero = 1,
      kNonZero = 2
   };

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TStreamerInfo::TStreamerInfo()
{
   fNumber   = -2;
   fClass    = 0;
   fElements = 0;
   fComp     = 0;
   fCompFull = 0;
   fCompOpt  = 0;
   fCheckSum = 0;
   fNdata    = 0;
   fNfulldata= 0;
   fNslots   = 0;
   fSize     = 0;
   fClassVersion = 0;
   fOnFileClassVersion = 0;
   fOldVersion = Class()->GetClassVersion();
   fNVirtualInfoLoc = 0;
   fVirtualInfoLoc = 0;

   fReadObjectWise = 0;
   fReadMemberWise = 0;
   fReadMemberWiseVecPtr = 0;
   fReadText = 0;
   fWriteObjectWise = 0;
   fWriteMemberWise = 0;
   fWriteMemberWiseVecPtr = 0;
   fWriteText = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TStreamerInfo object.

TStreamerInfo::TStreamerInfo(TClass *cl)
: TVirtualStreamerInfo(cl)
{
   fgCount++;
   fNumber   = fgCount;
   fClass    = cl;
   fElements = new TObjArray();
   fComp     = 0;
   fCompFull = 0;
   fCompOpt  = 0;
   fCheckSum = 0;
   fNdata    = 0;
   fNfulldata= 0;
   fNslots   = 0;
   fSize     = 0;
   fClassVersion = fClass->GetClassVersion();
   fOnFileClassVersion = 0;
   fOldVersion = Class()->GetClassVersion();
   fNVirtualInfoLoc = 0;
   fVirtualInfoLoc = 0;

   fReadObjectWise = 0;
   fReadMemberWise = 0;
   fReadMemberWiseVecPtr = 0;
   fReadText = 0;
   fWriteObjectWise = 0;
   fWriteMemberWise = 0;
   fWriteMemberWiseVecPtr = 0;
   fWriteText = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// TStreamerInfo dtor.

TStreamerInfo::~TStreamerInfo()
{
   // Note: If a StreamerInfo is loaded from a file and is the same information
   // as an existing one, it is assigned the same "unique id" and we need to double
   // check before removing it from the global list.
   if (fNumber >=0 && gROOT->GetListOfStreamerInfo()->At(fNumber) == this)
      gROOT->GetListOfStreamerInfo()->RemoveAt(fNumber);

   delete [] fComp;     fComp     = 0;
   delete [] fCompFull; fCompFull = 0;
   delete [] fCompOpt;  fCompOpt  = 0;
   delete [] fVirtualInfoLoc; fVirtualInfoLoc =0;

   delete fReadObjectWise;
   delete fReadMemberWise;
   delete fReadMemberWiseVecPtr;
   delete fReadText;
   delete fWriteObjectWise;
   delete fWriteMemberWise;
   delete fWriteMemberWiseVecPtr;
   delete fWriteText;

   if (!fElements) return;
   fElements->Delete();
   delete fElements; fElements=0;
}

////////////////////////////////////////////////////////////////////////////////
/// Makes sure kBuildOldUsed set once Build or BuildOld finishes.
/// Makes sure kBuildRunning reset once Build finishes.

namespace {
   struct TPreventRecursiveBuildGuard {
      TPreventRecursiveBuildGuard(TStreamerInfo* info): fInfo(info) {
         fInfo->SetBit(TStreamerInfo::kBuildRunning);
         fInfo->SetBit(TStreamerInfo::kBuildOldUsed);
      }
      ~TPreventRecursiveBuildGuard() {
         fInfo->ResetBit(TStreamerInfo::kBuildOldUsed);
         fInfo->ResetBit(TStreamerInfo::kBuildRunning);
      }
      TStreamerInfo* fInfo;
   };

}

////////////////////////////////////////////////////////////////////////////////
/// Build the I/O data structure for the current class version.
///
/// A list of TStreamerElement derived classes is built by scanning
/// one by one the list of data members of the analyzed class.
void TStreamerInfo::Build(Bool_t isTransient)
{
   // Did another thread already do the work?
   if (fIsCompiled) return;

   R__LOCKGUARD(gInterpreterMutex);

   // Did another thread already do the work while we were waiting ..
   if (fIsCompiled) return;

   // Has Build already been run?
   if (fIsBuilt) return;

   // Are we recursing on ourself?
   if (TestBit(TStreamerInfo::kBuildRunning)) return;

   // This is used to avoid unwanted recursive call to Build or BuildOld.
   TPreventRecursiveBuildGuard buildGuard(this);

   if (fClass->GetCollectionProxy()) {
      TVirtualCollectionProxy *proxy = fClass->GetCollectionProxy();
      TString title;
      if (proxy->GetValueClass()) {
         title.Form("<%s%s> Used to call the proper TStreamerInfo case",proxy->GetValueClass()->GetName(),proxy->HasPointers() ? "*" : "");
      } else {
         title .Form("<%s%s> Used to call the proper TStreamerInfo case",TDataType::GetTypeName(proxy->GetType()),proxy->HasPointers() ? "*" : "");
      }
      TStreamerElement* element = new TStreamerSTL("This", title.Data(), 0, fClass->GetName(), *proxy, 0);
      fElements->Add(element);
      Compile();
      fCheckSum = fClass->GetCheckSum();
      fIsBuilt = kTRUE;
      return;
   }

   TStreamerElement::Class()->IgnoreTObjectStreamer();

   fClass->BuildRealData(nullptr, isTransient);

   fCheckSum = fClass->GetCheckSum();

   Bool_t needAllocClass = kFALSE;
   Bool_t wasCompiled = fComp != 0;
   ROOT::TSchemaRuleSet::TMatches rules;
   if (fClass->GetSchemaRules()) {
       rules = fClass->GetSchemaRules()->FindRules(fClass->GetName(), fClassVersion);
   }

   //
   // Iterate over base classes.
   //

   // ROOT-9808: Here we skip the investigations of the base classes in case
   // this is a pair, otherwise, on some STL implementations, it can happen that
   // pair has mother classes which are an internal implementation detail and
   // would result in bogus messages printed on screen.
   if (!TClassEdit::IsStdPair(fClass->GetName())) {
      const bool isCollection = fClass->GetCollectionProxy();
      const bool isString = !strcmp(fClass->GetName(), "string");
      TBaseClass* base = 0;
      TIter nextb(fClass->GetListOfBases());
      while ((base = (TBaseClass*)nextb())) {
         TStreamerElement* element = 0;
         Int_t offset = base->GetDelta();
         if (offset == kMissing) {
            continue;
         }
         if (offset == kNeedObjectForVirtualBaseClass) {
            if (!isTransient)
               Error("Build()", "Cannot stream virtual base %s of class %s",
                     base->GetName(), fClass->GetName());
            continue;
         }
         const char* bname  = base->GetName();
         const char* btitle = base->GetTitle();
         // this case appears with STL collections as base class.
         if (!strcmp(bname, "string")) {
            element = new TStreamerSTLstring(bname, btitle, offset, bname, kFALSE);
         } else if (base->IsSTLContainer()) {
            TVirtualCollectionProxy *proxy = base->GetClassPointer()->GetCollectionProxy();
            if (proxy) element = new TStreamerSTL(bname, btitle, offset, bname, *proxy, kFALSE);
            else       element = new TStreamerSTL(bname, btitle, offset, bname, 0, kFALSE);
            if (fClass->IsLoaded() && ((TStreamerSTL*)element)->GetSTLtype() != ROOT::kSTLvector) {
               if (!element->GetClassPointer()->IsLoaded()) {
                  if (!isTransient)
                     Error("Build","The class \"%s\" is compiled and its base class \"%s\" is a collection and we do not have a dictionary for it, we will not be able to read or write this base class.",GetName(),bname);
                  delete element;
                  continue;
               }
            }
         } else {
            element = new TStreamerBase(bname, btitle, offset, isTransient);
            TClass* clm = element->GetClassPointer();
            if (!clm) {
               // We have no information about the class yet, except that since it
               // is a base class, we know it is a class.  So let's create it (in v5
               // it would have been created as a side effect of the dictionary of
               // for the derived class having a forward declaration of the base class).
               clm = new TClass(bname,1,TClass::kForwardDeclared, true /*silent*/);
               Warning("Build", "%s: base class %s has no streamer or dictionary it will not be saved", GetName(), clm->GetName());
               element->Init(0);
            } else {
               // Now part of the TStreamerBase constructor.
               // clm->GetStreamerInfo();
               if ((clm == TObject::Class()) && fClass->CanIgnoreTObjectStreamer()) {
                  // -- An ignored TObject base class.
                  // Note: The TClass kIgnoreTObjectStreamer == BIT(15), but
                  // the TStreamerInfo kIgnoreTobjectStreamer == BIT(13) which
                  // is confusing.
                  SetBit(kIgnoreTObjectStreamer);
                  // Flag the element to be ignored by setting its type to -1.
                  // This flag will be used later by Compile() to prevent this
                  // element from being inserted into the compiled info.
                  element->SetType(-1);
               }
               if (!isTransient && !clm->IsLoaded() && !(isCollection || isString)) {
                  // Don't complain about the base classes of collections nor of
                  // std::string.
                  Warning("Build", "%s: base class %s has no streamer or dictionary it will not be saved", GetName(), clm->GetName());
               }
            }
         }
         if (element) {
            fElements->Add(element);
         }
      } // end of base class loop
   }

   //
   // Iterate over data members.
   //

   Int_t dsize;
   TDataMember* dm = 0;
   std::string typeNameBuf;
   std::string trueTypeNameBuf;
   TIter nextd(fClass->GetListOfDataMembers());
   while ((dm = (TDataMember*) nextd())) {
      if (fClass->GetClassVersion() == 0) {
         continue;
      }
      if (!dm->IsPersistent()) {
         continue;
      }
      if (dm->Property() & kIsUnionMember) {
         continue;
      }
      TMemberStreamer* streamer = 0;
      Int_t offset = GetDataMemberOffset(dm, streamer);
      if (offset == kMissing) {
         continue;
      }
      TStreamerElement* element = 0;
      dsize = 0;

      // Save some useful variables
      const char* dmName  = dm->GetName();
      const char* dmTitle = dm->GetTitle();
      const char* dmType  = dm->GetTypeName();
      const char* dmFull  = dm->GetTrueTypeName(); // Used to be GetFullTypeName ...
      Bool_t dmIsPtr = dm->IsaPointer();
      TDataType* dt(nullptr);
      Int_t ndim = dm->GetArrayDim();
      std::array<Int_t, 5> maxIndices; // 5 is the maximum supported in TStreamerElement::SetMaxIndex
      Bool_t isStdArray(kFALSE);

      // Let's treat the unique_ptr case
      bool nameChanged;
      trueTypeNameBuf = typeNameBuf = TClassEdit::GetNameForIO(dmFull, TClassEdit::EModType::kNone, &nameChanged);
      if (nameChanged) {
         if (TClassEdit::IsUniquePtr(dmFull)) {
            dmIsPtr = true;
         }
         while(typeNameBuf.back() == '*') typeNameBuf.pop_back();
         dmFull = trueTypeNameBuf.c_str();
         dmType = typeNameBuf.c_str();
      }
      if ((isStdArray = TClassEdit::IsStdArray(dmType))){ // We tackle the std array case
         TClassEdit::GetStdArrayProperties(dmType,
                                           typeNameBuf,
                                           maxIndices,
                                           ndim);
         trueTypeNameBuf = typeNameBuf;
         while(typeNameBuf.back() == '*') typeNameBuf.pop_back();
         dmFull = dmType = typeNameBuf.c_str();
         dt = gROOT->GetType(dmType);
      }

      TDataMember* dmCounter = 0;
      if (dmIsPtr) {
         //
         // look for a pointer data member with a counter
         // in the comment string, like so:
         //
         //      int n;
         //      double* MyArray; //[n]
         //
         const char* lbracket = TVirtualStreamerInfo::GetElementCounterStart(dmTitle);
         const char* rbracket = ::strchr(dmTitle, ']');
         if (lbracket && rbracket) {
            const char* counterName = dm->GetArrayIndex();
            TRealData* rdCounter = (TRealData*) fClass->GetListOfRealData()->FindObject(counterName);
            if (!rdCounter || rdCounter->TestBit(TRealData::kTransient)) {
               if (!isTransient)
                  Error("Build", "%s, discarding: %s %s, illegal %s\n", GetName(), dmFull, dmName, dmTitle);
               continue;
            }
            dmCounter = rdCounter->GetDataMember();
            TDataType* dtCounter = dmCounter->GetDataType();
            Bool_t isInteger = dtCounter && ((dtCounter->GetType() == 3) || (dtCounter->GetType() == 13));
            if (!dtCounter || !isInteger) {
               if (!isTransient)
                  Error("Build", "%s, discarding: %s %s, illegal [%s] (must be Int_t)\n", GetName(), dmFull, dmName, counterName);
               continue;
            }
            TStreamerBasicType* bt = TStreamerInfo::GetElementCounter(counterName, dmCounter->GetClass());
            if (!bt) {
               if (dmCounter->GetClass()->Property() & kIsAbstract) {
                  continue;
               }
               if (!isTransient)
                  Error("Build", "%s, discarding: %s %s, illegal [%s] must be placed before \n", GetName(), dmFull, dmName, counterName);
               continue;
            }
         }
      }
      if (!dt && !isStdArray) dt = dm->GetDataType();
      if (dt) {
         // found a basic type
         Int_t dtype = dt->GetType();
         dsize = dt->Size();
         if (!dmCounter && (strstr(dmFull, "char*") || strstr(dmFull, "Char_t*"))) {
            dtype = kCharStar;
            dsize = sizeof(char*);
         }
         if (dtype == kOther_t || dtype == kNoType_t) {
            if (!isTransient)
               Error("Build", "%s, unknown type: %s %s", GetName(), dmFull, dmName);
            continue;
         } else if (dmIsPtr && (dtype != kCharStar)) {
            if (dmCounter) {
               // data member is pointer to an array of basic types
               element = new TStreamerBasicPointer(dmName, dmTitle, offset, dtype, dm->GetArrayIndex(), dmCounter->GetClass()->GetName(), dmCounter->GetClass()->GetClassVersion(), dmFull);
            } else {
               if ((fName == "TString") || (fName == "TClass")) {
                  continue;
               }
               if (!isTransient)
                  Error("Build", "%s, discarding: %s %s, no [dimension]\n", GetName(), dmFull, dmName);
               continue;
            }
         } else {
            // data member is a basic type
            if ((fClass == TObject::Class()) && !strcmp(dmName, "fBits")) {
               //printf("found fBits, changing dtype from %d to 15\n", dtype);
               dtype = kBits;
            }
            // Here we treat data members such as int, float, double[4]
            element = new TStreamerBasicType(dmName, dmTitle, offset, dtype, dmFull);
         }
      } else {
         // try STL container or string
         static const char* full_string_name = "basic_string<char,char_traits<char>,allocator<char> >";
         if (!strcmp(dmType, "string") || !strcmp(dmType, "std::string") || !strcmp(dmType, full_string_name)) {
            element = new TStreamerSTLstring(dmName, dmTitle, offset, dmFull, dmIsPtr);
         } else if (dm->IsSTLContainer()) {
            TVirtualCollectionProxy *proxy = TClass::GetClass(dmType /* the underlying type */)->GetCollectionProxy();
            if (proxy) element = new TStreamerSTL(dmName, dmTitle, offset, dmFull, *proxy, dmIsPtr);
            else element = new TStreamerSTL(dmName, dmTitle, offset, dmFull, dmFull, dmIsPtr);
            bool hasCustomAlloc = proxy ? proxy->GetProperties() & TVirtualCollectionProxy::kCustomAlloc : kFALSE;
            if (((TStreamerSTL*)element)->GetSTLtype() != ROOT::kSTLvector || hasCustomAlloc) {
               auto printErrorMsg = [&](const char* category)
                  {
                     if (!isTransient)
                        Error("Build","The class \"%s\" is %s and for its data member \"%s\" we do not have a dictionary for the collection \"%s\". Because of this, we will not be able to read or write this data member.",GetName(), category, dmName, dmType);
                  };
               if (fClass->IsLoaded()) {
                  if (!element->GetClassPointer()->IsLoaded()) {
                     printErrorMsg("compiled");
                     delete element;
                     continue;
                  }
               } else if (fClass->GetState() == TClass::kInterpreted) {
                  if (element->GetClassPointer()->GetCollectionProxy()->GetProperties() & TVirtualCollectionProxy::kIsEmulated) {
                     printErrorMsg("interpreted");
                     delete element;
                     continue;
                  }
               }
            }
         } else {
            TClass* clm = TClass::GetClass(dmType);
            if (!clm) {
               if (!isTransient)
                  Error("Build", "%s, unknown type: %s %s\n", GetName(), dmFull, dmName);
               continue;
            }
            if (isStdArray) {
               // We do not want to rebuild the streamerinfo of an std::array<T,N> asking the dm->GetUnitSize(), but rather of T only.

               dsize = clm->Size();
            }
            if (dmIsPtr) {
               // a pointer to a class
               if (dmCounter) {
                  element = new TStreamerLoop(dmName, dmTitle, offset, dm->GetArrayIndex(), dmCounter->GetClass()->GetName(), dmCounter->GetClass()->GetClassVersion(), dmFull);
               } else {
                  if (clm->IsTObject()) {
                     element = new TStreamerObjectPointer(dmName, dmTitle, offset, dmFull);
                  } else {
                     element = new TStreamerObjectAnyPointer(dmName, dmTitle, offset, dmFull);
                     if (!isTransient && !streamer && !clm->GetStreamer() && !clm->IsLoaded() && !clm->fIsSyntheticPair) {
                        Error("Build", "%s: %s has no streamer or dictionary, data member %s will not be saved",
                              GetName(), dmFull, dmName);
                     }
                  }
               }
            } else if (clm->IsTObject()) {
               element = new TStreamerObject(dmName, dmTitle, offset, dmFull);
            } else if ((clm == TString::Class()) && !dmIsPtr) {
               element = new TStreamerString(dmName, dmTitle, offset);
            } else {
               element = new TStreamerObjectAny(dmName, dmTitle, offset, dmFull);
               if (!isTransient && !streamer && !clm->GetStreamer() && !clm->IsLoaded() && !clm->fIsSyntheticPair) {
                  Warning("Build", "%s: %s has no streamer or dictionary, data member \"%s\" will not be saved",
                          GetName(), dmFull, dmName);
               }
            }
         }
      }
      if (!element) {
         // If we didn't make an element, there is nothing to do.
         continue;
      }
      if (!dsize) {
         dsize = dm->GetUnitSize();
      }
      for (Int_t i = 0; i < ndim; ++i) {
         auto maxIndex = 0;
         if (isStdArray) maxIndex = maxIndices[i];
         else maxIndex = dm->GetMaxIndex(i);
         element->SetMaxIndex(i, maxIndex);
      }
      element->SetArrayDim(ndim);
      // If the datamember was a int[4] this is 4, if double[3][2] 3*2=6
      Int_t narr = element->GetArrayLength();
      if (!narr) {
         narr = 1;
      }
      element->SetSize(dsize*narr);
      element->SetStreamer(streamer);
      if (!streamer) {
         Int_t k = element->GetType();
         if (k == kStreamer) {
            //if ((k == kSTL) || (k == kSTL + kOffsetL) || (k == kStreamer) || (k == kStreamLoop))
            element->SetType(-1);
         }
      }

      if ( !wasCompiled && (rules && rules.HasRuleWithSource( element->GetName(), kTRUE )) ) {
         needAllocClass = kTRUE;

         // If this is optimized to re-use TStreamerElement(s) in case of variable renaming,
         // then we must revisit the code in TBranchElement::InitInfo that recalculate the
         // fID (i.e. the index of the TStreamerElement to be used for streaming).

         TStreamerElement *cached = element;
         // Now that we are caching the unconverted element, we do not assign it to the real type even if we could have!
         if (element->GetNewType()>0 /* intentionally not including base class for now */
             && rules && !rules.HasRuleWithTarget( element->GetName(), kTRUE ) )
         {
            TStreamerElement *copy = (TStreamerElement*)element->Clone();
            fElements->Add(copy);
            copy->SetBit(TStreamerElement::kRepeat);
            cached = copy;

            // Warning("BuildOld","%s::%s is not set from the version %d of %s (You must add a rule for it)\n",GetName(), element->GetName(), GetClassVersion(), GetName() );
         } else {
            // If the element is just cached and not repeat, we need to inject an element
            // to insure the writing.
            TStreamerElement *writecopy = (TStreamerElement*)element->Clone();
            fElements->Add(element);
            writecopy->SetBit(TStreamerElement::kWrite);
            writecopy->SetNewType( writecopy->GetType() );
            writecopy->SetOffset( element->GetOffset() );
            // Put the write element after the read element (that does caching).
            element = writecopy;
         }
         cached->SetBit(TStreamerElement::kCache);
         cached->SetNewType( cached->GetType() );
      }

      fElements->Add(element);
   } // end of member loop

   // Now add artificial TStreamerElement (i.e. rules that creates new members or set transient members).
   InsertArtificialElements(rules);

   if (needAllocClass) {
      TStreamerInfo *infoalloc  = (TStreamerInfo *)Clone(TString::Format("%s@@%d",GetName(),GetClassVersion()));
      if (!infoalloc) {
         if (!isTransient)
            Error("Build","Could you create a TStreamerInfo for %s\n",TString::Format("%s@@%d",GetName(),GetClassVersion()).Data());
      } else {
         // Tell clone we should rerun BuildOld
         infoalloc->SetBit(kBuildOldUsed,false);
         // Temporarily mark it as built to avoid the BuildCheck from removing
         // Technically we only need to do this for the 'current' StreamerInfo
         fIsBuilt = kTRUE;
         infoalloc->BuildCheck();
         infoalloc->BuildOld();
         fIsBuilt = kFALSE;
         TClass *allocClass = infoalloc->GetClass();

         {
            TIter next(fElements);
            TStreamerElement* element;
            while ((element = (TStreamerElement*) next())) {
               if (element->TestBit(TStreamerElement::kRepeat) && element->IsaPointer()) {
                  TStreamerElement *other = (TStreamerElement*) infoalloc->GetElements()->FindObject(element->GetName());
                  if (other) {
                     other->SetBit(TStreamerElement::kDoNotDelete);
                  }
               }
            }
            infoalloc->GetElements()->Compress();
         }
         {
            TIter next(fElements);
            TStreamerElement* element;
            while ((element = (TStreamerElement*) next())) {
            if (element->TestBit(TStreamerElement::kCache)) {
               element->SetOffset(infoalloc->GetOffset(element->GetName()));
            }
            }
         }

         TStreamerElement *el = new TStreamerArtificial("@@alloc","", 0, TStreamerInfo::kCacheNew, allocClass->GetName());
         R__TObjArray_InsertAt( fElements, el, 0 );

         el = new TStreamerArtificial("@@dealloc","", 0, TStreamerInfo::kCacheDelete, allocClass->GetName());
         fElements->Add( el );
      }
   }

   //
   // Make a more compact version.
   //
   Compile();
   fIsBuilt = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if built and consistent with the class dictionary.
/// This method is called by TFile::ReadStreamerInfo.

void TStreamerInfo::BuildCheck(TFile *file /* = 0 */, Bool_t load /* = kTRUE */)
{
   R__LOCKGUARD(gInterpreterMutex);

   fClass = TClass::GetClass(GetName(), load);
   if (!fClass) {
      // fClassVersion should have been a Version_t and/or Version_t
      // should have been an Int_t.  Changing the on-file format
      // of the StreamerInfo is 'hard' (for forward compatibility), so
      // leave it as is for now.
      fClass = new TClass(GetName(), (Version_t)fClassVersion);
      fClass->SetBit(TClass::kIsEmulation);

      // Case of a custom collection (the user provided a CollectionProxy
      // for a class that is not an STL collection).
      if (GetElements()->GetEntries() == 1) {
         TObject *element = GetElements()->UncheckedAt(0);
         Bool_t isstl = element && strcmp("This",element->GetName())==0;
         if (isstl) {
            if (element->GetTitle()[0] == '<') {
               // We know the content.
               TString content = element->GetTitle();
               Int_t level = 1;
               for(Int_t c = 1; c < content.Length(); ++c) {
                  if (content[c] == '<') ++level;
                  else if (content[c] == '>') --level;
                  if (level == 0) {
                     content.Remove(c+1);
                     break;
                  }
               }
               content.Prepend("vector");
               TClass *clequiv = TClass::GetClass(content);
               TVirtualCollectionProxy *proxy = clequiv->GetCollectionProxy();
               if (gDebug > 1)
                  Info("BuildCheck",
                       "Update the collection proxy of the class \"%s\" \n"
                       "\tto be similar to \"%s\".",
                    GetName(),content.Data());
               fClass->CopyCollectionProxy( *proxy );
            } else {
               Warning("BuildCheck", "\n\
   The class %s had a collection proxy when written but it is not an STL\n \
   collection and we did not record the type of the content of the collection.\n \
   We will claim the content is a bool (i.e. no data will be read).",
                       GetName());
            }
         }
      }

  } else {
      if (fClass->GetCollectionType() > ROOT::kNotSTL) {
         if (TClassEdit::IsSTLCont(fClass->GetName())) {
            // We have a collection that is indeed an STL collection,
            // we know we don't need its streamerInfo.
            SetBit(kCanDelete);
            return;
         }
      }
      if (fClass->fIsSyntheticPair) {
         // The format never change, no need to import the old StreamerInof
         // (which anyway would have issue when being setup due to the lack
         // of TDataMember in the TClass.
         SetBit(kCanDelete);
         return;
      }

      if (0 == strcmp("string",fClass->GetName())) {
         // We know we do not need any offset check for a string
         SetBit(kCanDelete);
         return;
      }

      const TObjArray *array = fClass->GetStreamerInfos();
      TStreamerInfo* info = 0;

      if (fClass->TestBit(TClass::kIsEmulation) && array->GetEntries()==0) {
         // We have an emulated class that has no TStreamerInfo, this
         // means it was created to insert a (default) rule.  Consequently
         // the error message about the missing dictionary was not printed.
         // For consistency, let's print it now!

         ::Warning("TClass::TClass", "no dictionary for class %s is available", GetName());
      }

      // Case of a custom collection (the user provided a CollectionProxy
      // for a class that is not an STL collection).
      if (GetElements()->GetEntries() == 1) {
         TObject *element = GetElements()->UncheckedAt(0);
         Bool_t isstl = element && strcmp("This",element->GetName())==0;
         if (isstl && !fClass->GetCollectionProxy()) {
            if (element->GetTitle()[0] == '<') {
               // We know the content.
               TString content = element->GetTitle();
               Int_t level = 1;
               for(Int_t c = 1; c < content.Length(); ++c) {
                  if (content[c] == '<') ++level;
                  else if (content[c] == '>') --level;
                  if (level == 0) {
                     content.Remove(c+1);
                     break;
                  }
               }
               content.Prepend("vector");
               TClass *clequiv = TClass::GetClass(content);
               TVirtualCollectionProxy *proxy = clequiv->GetCollectionProxy();
               if (gDebug > 1)
                  Info("BuildCheck",
                       "Update the collection proxy of the class \"%s\" \n"
                       "\tto be similar to \"%s\".",
                    GetName(),content.Data());
               fClass->CopyCollectionProxy( *proxy );
            } else {
               Warning("BuildCheck", "\n\
   The class %s had a collection proxy when written but it is not an STL\n \
   collection and we did not record the type of the content of the collection.\n \
   We will claim the content is a bool (i.e. no data will be read).",
                       GetName());
            }
            SetBit(kCanDelete);
            return;
         }
      }

      // If the user has not specified a class version (this _used to_
      // always be the case when the class is Foreign) or if the user
      // has specified a version to be explicitly 1. [We can not
      // distinguish the two cases using the information in the "on
      // file" StreamerInfo.]

      Bool_t searchOnChecksum = kFALSE;
      if (fClass->IsLoaded() && fClass->GetClassVersion() >= 2) {
         // We know for sure that the user specified the version.

         if (fOnFileClassVersion >= 2) {
            // The class version was specified when the object was
            // written

            searchOnChecksum = kFALSE;

         } else {
            // The class version was not specified when the object was
            // written OR it was specified to be 1.

            searchOnChecksum = kTRUE;
         }
      } else if (fClass->IsLoaded() && !fClass->IsForeign()) {
         // We are in the case where the class has a Streamer function.
         // and fClass->GetClassVersion is 1, we still assume that the
         // Class Version is specified (to be one).

         searchOnChecksum = kFALSE;

      } else if (fClass->IsLoaded() /* implied: && fClass->IsForeign() */ ) {
         // We are in the case of a Foreign class with no specified
         // class version.

         searchOnChecksum = kTRUE;

      }
      else {
         // We are in the case of an 'emulated' class.

         if (fOnFileClassVersion >= 2) {
            // The class version was specified when the object was
            // written

            searchOnChecksum = kFALSE;

         } else {
            // The class version was not specified when the object was
            // written OR it was specified to be 1.

            searchOnChecksum = kTRUE;

            TStreamerInfo* v1 = (TStreamerInfo*) array->At(1);
            if (v1) {
               if (fCheckSum != v1->GetCheckSum()) {
                  fClassVersion = array->GetLast() + 1;
               }
            }
         }
      }

      if (!searchOnChecksum) {
         if (fClassVersion < (array->GetEntriesFast() - 1)) {
            info = (TStreamerInfo*) array->At(fClassVersion);
         }
      } else {
         Int_t ninfos = array->GetEntriesFast() - 1;
         for (Int_t i = -1; i < ninfos; ++i) {
            info = (TStreamerInfo*) array->UncheckedAt(i);
            if (!info) {
               continue;
            }
            if (fCheckSum == info->GetCheckSum() && (info->GetOnFileClassVersion()==1 || info->GetOnFileClassVersion()==0)) {
               // We must match on the same checksum, an existing TStreamerInfo
               // for one of the 'unversioned' class layout (i.e. version was 1).
               fClassVersion = i;
               break;
            }
            info = 0;
         }
         if (info==0) {
            // Find an empty slot.
            ninfos = array->GetEntriesFast() - 1;
            Int_t slot = 1; // Start of Class version 1.
            while ((slot < ninfos) && (array->UncheckedAt(slot) != 0)) {
               ++slot;
            }
            fClassVersion = slot;
         }
      }

      // NOTE: Should we check if the already existing info is the same as
      // the current one? Yes
      // In case a class (eg Event.h) has a TClonesArray of Tracks, it could be
      // that the old info does not have the class name (Track) in the data
      // member title. Set old title to new title
      if (info) {
         // We found an existing TStreamerInfo for our ClassVersion
         Bool_t match = kTRUE;
         Bool_t done = kFALSE;
         Bool_t oldIsNonVersioned = kFALSE;
         if (fClassVersion!=0 && !fClass->TestBit(TClass::kWarned) && (fClassVersion == info->GetClassVersion()) && (fCheckSum != info->GetCheckSum())) {
            // The TStreamerInfo's checksum is different from the checksum for the compile class.

            match = kFALSE;
            oldIsNonVersioned = info->fOnFileClassVersion==1 && info->fClassVersion != 1;

            if (fClass->IsLoaded() && (fClassVersion == fClass->GetClassVersion()) && fClass->HasDataMemberInfo()) {
               // In the case where the read-in TStreamerInfo does not
               // match in the 'current' in memory TStreamerInfo for
               // a non foreign class (we can not get here if this is
               // a foreign class so we do not need to test it),
               // we need to add this one more test since the CINT behaviour
               // with enums changed over time, so verify the checksum ignoring
               // members of type enum. We also used to not count the //[xyz] comment
               // in the checksum, so test for that too.
               if (  (fCheckSum == fClass->GetCheckSum() || fClass->MatchLegacyCheckSum(fCheckSum) )
                     &&(info->GetCheckSum() == fClass->GetCheckSum() || fClass->MatchLegacyCheckSum(info->GetCheckSum()))
                     )
                  {
                     match = kTRUE;
                  }
               if (fOldVersion <= 2) {
                  // Names of STL base classes was modified in vers==3. Allocators removed
                  // (We could be more specific (see test for the same case below)
                  match = kTRUE;
               }
               if (!match && CompareContent(0,info,kFALSE,kFALSE,file)) {
                  match = kTRUE;
               }
#ifdef TEST_FOR_BACKWARD_COMPATIBILITY_ABSTRACT_CLASSES
               if (!match && file->GetVersion() < 51800 && fClass && (fClass->Property() & kIsAbstract)
                   && fClass->GetListOfDataMembers()->GetEntries() != 0)
               {
                  // In some instances of old files (v5.17 and less), some StreamerInfo for
                  // an abstract class where not written correctly, and add no
                  // data member listed.  If in addition one of the data member
                  // was declared using a typedef _and_ the current class definition
                  // uses a different typedef, we are unable to recalculate the
                  // checksum as it was, because the information is missing from
                  // the StreamerInfo, and for the same reason CompareContent can
                  // not know whether this is okay or not ...
                  //
                  // Since this is such an unlikely scenario, let's complain
                  // about it anyway (The class layout *may* have changed, we
                  // don't know).

                  // if (this has only base classes) {
                  //    match = kTRUE;
                  // }
               }
#endif
            } else {
               // The on-file TStreamerInfo's checksum differs from the checksum of a TStreamerInfo on another file.

               match = kFALSE;
               oldIsNonVersioned = info->fOnFileClassVersion==1 && info->fClassVersion != 1;

               // In the case where the read-in TStreamerInfo does not
               // match in the 'current' in memory TStreamerInfo for
               // a non foreign class (we can not get here if this is
               // a foreign class so we do not need to test it),
               // we need to add this one more test since the CINT behaviour
               // with enums changed over time, so verify the checksum ignoring
               // members of type enum. We also used to not count the //[xyz] comment
               // in the checksum, so test for that too.
               if (fCheckSum == info->GetCheckSum(TClass::kCurrentCheckSum)
                   || info->MatchLegacyCheckSum(fCheckSum)
                   || GetCheckSum(TClass::kCurrentCheckSum) == info->fCheckSum
                   || MatchLegacyCheckSum(info->GetCheckSum())
                   || GetCheckSum(TClass::kCurrentCheckSum) == info->GetCheckSum(TClass::kCurrentCheckSum))
                  {
                     match = kTRUE;
                  }
               if (fOldVersion <= 2) {
                  // Names of STL base classes was modified in vers==3. Allocators removed
                  // (We could be more specific (see test for the same case below)
                  match = kTRUE;
               }
               if (!match && CompareContent(0,info,kFALSE,kFALSE,file)) {
                  match = kTRUE;
               }
            }
         }
         if (info->IsBuilt()) {
            SetBit(kCanDelete);
            fNumber = info->GetNumber();
            Int_t nel = fElements->GetEntriesFast();
            TObjArray* elems = info->GetElements();
            TStreamerElement* e1 = 0;
            TStreamerElement* e2 = 0;
            for (Int_t i = 0; i < nel; ++i) {
               e1 = (TStreamerElement*) fElements->UncheckedAt(i);
               e2 = (TStreamerElement*) elems->At(i);
               if (!e1 || !e2) {
                  continue;
               }
               if (strlen(e1->GetTitle()) != strlen(e2->GetTitle())) {
                  e2->SetTitle(e1->GetTitle());
               }
            }

            done = kTRUE;
         } else {
            fClass->RemoveStreamerInfo(fClassVersion);
            info = 0;
         }
         TString origin;
         if (!match && !fClass->TestBit(TClass::kWarned)) {
            if (oldIsNonVersioned) {
               if (file) {
                  Warning("BuildCheck", "\n\
   The class %s transitioned from not having a specified class version\n\
   to having a specified class version (the current class version is %d).\n\
   However too many different non-versioned layouts of the class have been\n\
   loaded so far.  This prevent the proper reading of objects written with\n\
   the class layout version %d, in particular from the file:\n\
   %s.\n\
   To work around this issue, load fewer 'old' files in the same ROOT session.",
                          GetName(),fClass->GetClassVersion(),fClassVersion,file->GetName());
               } else {
                  Warning("BuildCheck", "\n\
   The class %s transitioned from not having a specified class version\n\
   to having a specified class version (the current class version is %d).\n\
   However too many different non-versioned layouts of the class have been\n\
   loaded so far.  This prevent the proper reading of objects written with\n\
   the class layout version %d.\n\
   To work around this issue, load fewer 'old' files in the same ROOT session.",
                          GetName(),fClass->GetClassVersion(),fClassVersion);
               }
            } else {
               if (file) {
                  if (done) {
                     Warning("BuildCheck", "\n\
   The StreamerInfo for version %d of class %s read from the file %s\n\
   has a different checksum than the previously loaded StreamerInfo.\n\
   Reading objects of type %s from the file %s \n\
   (and potentially other files) might not work correctly.\n\
   Most likely the version number of the class was not properly\n\
   updated [See ClassDef(%s,%d)].",
                             fClassVersion, GetName(), file->GetName(), GetName(), file->GetName(), GetName(), fClassVersion);
                  } else {
                     Warning("BuildCheck", "\n\
   The StreamerInfo from %s does not match existing one (%s:%d)\n\
   The existing one has not been used yet and will be discarded.\n\
   Reading the file %s will work properly, however writing object of\n\
   type %s will not work properly.  Most likely the version number\n\
   of the class was not properly updated [See ClassDef(%s,%d)].",
                             file->GetName(), GetName(), fClassVersion,file->GetName(),GetName(), GetName(), fClassVersion);
                  }
               } else {
                  if (done) {
                     Warning("BuildCheck", "\n\
   The StreamerInfo for version %d of class %s\n\
   has a different checksum than the previously loaded StreamerInfo.\n\
   Reading objects of type %s\n\
   (and potentially other files) might not work correctly.\n\
   Most likely the version number of the class was not properly\n\
   updated [See ClassDef(%s,%d)].",
                             fClassVersion, GetName(), GetName(), GetName(), fClassVersion);
                  } else {
                     Warning("BuildCheck", "\n\
   The StreamerInfo from %s does not match existing one (%s:%d)\n\
   The existing one has not been used yet and will be discarded.\n\
   Reading should work properly, however writing object of\n\
   type %s will not work properly.  Most likely the version number\n\
   of the class was not properly updated [See ClassDef(%s,%d)].",
                             file->GetName(), GetName(), fClassVersion, GetName(), GetName(), fClassVersion);
                  }
               }
            }
            CompareContent(0,info,kTRUE,kTRUE,file);
            fClass->SetBit(TClass::kWarned);
         }
         if (done) {
            return;
         }
      }
      // The slot was free, however it might still be reserved for the current
      // loaded version of the class
      if (fClass->IsLoaded()
          && fClass->HasDataMemberInfo()
          && (fClassVersion != 0) // We don't care about transient classes
          && (fClassVersion == fClass->GetClassVersion())
          && (fCheckSum != fClass->GetCheckSum())) {

         // If the old TStreamerInfo matches the in-memory one when we either
         //   - ignore the members of type enum
         // or
         //   - ignore the comments annotation (//[xyz])
         // we can accept the old TStreamerInfo.

         if (!fClass->MatchLegacyCheckSum(fCheckSum)) {

            Bool_t warn = !fClass->TestBit(TClass::kWarned);
            if (warn) {
               warn = !CompareContent(fClass,0,kFALSE,kFALSE,file);
            }
#ifdef TEST_FOR_BACKWARD_COMPATIBILITY_ABSTRACT_CLASSES
            if (warn && file->GetVersion() < 51800 && fClass && (fClass->Property() & kIsAbstract)
                && fClass->GetListOfDataMembers()->GetEntries() != 0)
            {
               // In some instances of old files (v5.17 and less), some StreamerInfo for
               // an abstract class where not written correctly, and add no
               // data member listed.  If in addition one of the data member
               // was declared using a typedef _and_ the current class definition
               // uses a different typedef, we are unable to recalculate the
               // checksum as it was, because the information is missing from
               // the StreamerInfo, and for the same reason CompareContent can
               // not know whether this is okay or not ...
               //
               // Since this is such an unlikely scenario, let's complain
               // about it anyway (The class layout *may* have changed, we
               // don't know).

               // if (this has only base classes) {
               //    warn = kFALSE;
               // }
            }
#endif // TEST_FOR_BACKWARD_COMPATIBILITY
            if (warn && (fOldVersion <= 2)) {
               // Names of STL base classes was modified in vers==3. Allocators removed
               //
               TIter nextBC(fClass->GetListOfBases());
               TBaseClass* bc = 0;
               while ((bc = (TBaseClass*) nextBC())) {
                  if (bc->GetClassPointer()->GetCollectionType()) {
                     warn = kFALSE;
                  }
               }
            }
            if (warn) {
               if (file) {
                  Warning("BuildCheck", "\n\
   The StreamerInfo of class %s read from file %s\n\
   has the same version (=%d) as the active class but a different checksum.\n\
   You should update the version to ClassDef(%s,%d).\n\
   Do not try to write objects with the current class definition,\n\
   the files will not be readable.\n", GetName(), file->GetName(), fClassVersion, GetName(), fClassVersion + 1);
               } else {
                  Warning("BuildCheck", "\n\
   The StreamerInfo of class %s \n\
   has the same version (=%d) as the active class but a different checksum.\n\
   You should update the version to ClassDef(%s,%d).\n\
   Do not try to write objects with the current class definition,\n\
   the files will not be readable.\n", GetName(), fClassVersion, GetName(), fClassVersion + 1);
               }
               CompareContent(fClass,0,kTRUE,kTRUE,file);
               fClass->SetBit(TClass::kWarned);
            }
         } else {
            if (!fClass->IsVersioned()) {
               Fatal("BuildCheck", "\n\
   The StreamerInfo of unversioned class %s \n\
   has the same version (=%d) as the active class but an old checksum.\n\
   This should not happen. An assert will follow.\n", GetName(), fClassVersion);
            }
         }
      }
      if (!fClass->IsLoaded() &&  this->fOnFileClassVersion>1)
      {
         ROOT::ResetClassVersion(fClass,(const char*)-1, this->fClassVersion);
      }
   }
   // FIXME: This code can never execute because Build() calls
   // TStreamerElement::Class()->IgnoreTObjectStreamer()
   // so our bits are never saved to the file.
   if (TestBit(kIgnoreTObjectStreamer)) {
      fClass->IgnoreTObjectStreamer();
   }
   if ((fClassVersion < -1) || (fClassVersion > 65000)) {
      printf("ERROR reading TStreamerInfo: %s fClassVersion=%d\n", GetName(), fClassVersion);
      SetBit(kCanDelete);
      fNumber = -1;
      return;
   }

   if (!fClass->TestBit(TClass::kWarned)
       && fClass->GetState() >= TClass::kInterpreted
       && GetCheckSum() != fClass->GetCheckSum()
       && fClassVersion == fClass->GetClassVersion()) {
      // We got here, thus we are a perfect alias for the current streamerInfo,
      // but we might had odd v5 style name spelling, so let's prefer the
      // current one.
      auto maininfo = fClass->GetStreamerInfo();
      if (maininfo) {
         fNumber = maininfo->GetNumber(); // For ReadStreamerInfo to record the expected slot.
      }
      SetBit(kCanDelete);
      return;
   }

   fClass->RegisterStreamerInfo(this);
   ++fgCount;
   fNumber = fgCount;

   // Since we just read this streamerInfo from file, it has already been built.
   fIsBuilt = kTRUE;

   //add to the global list of StreamerInfo
   TObjArray* infos = (TObjArray*) gROOT->GetListOfStreamerInfo();
   infos->AddAtAndExpand(this, fNumber);
}

////////////////////////////////////////////////////////////////////////////////
/// Create an Emulation TStreamerInfo object.

void TStreamerInfo::BuildEmulated(TFile *file)
{
   R__LOCKGUARD(gInterpreterMutex);

   TString duName;
   R__ASSERT(file);
   Int_t fv = file->GetVersion()%100000;
   R__ASSERT(fv < 30000);
   fClassVersion = -1;
   fCheckSum = 2001;
   TObjArray *elements = GetElements();
   Int_t ndata = elements ? elements->GetEntries() : 0;
   for (Int_t i=0;i < ndata;i++) {
      TStreamerElement *element = (TStreamerElement*)elements->UncheckedAt(i);
      if (!element) break;
      int ty = element->GetType();
      if (ty < kChar || ty >kULong+kOffsetL)    continue;
      if (ty == kLong)                         element->SetType(kInt);
      if (ty == kULong)                         element->SetType(kUInt);
      if (ty == kLong + kOffsetL)                element->SetType(kInt + kOffsetL);
      if (ty == kULong + kOffsetL)                element->SetType(kUInt + kOffsetL);
      if (ty <= kULong)                         continue;
      duName = element->GetName();
      duName.Append("QWERTY");
      TStreamerBasicType *bt = new TStreamerBasicType(duName, "", 0, kInt,"Int_t");
      {for (int j=ndata-1;j>=i;j--) {elements->AddAtAndExpand(elements->At(j),j+1);}}
      elements->AddAt(bt,i);
      ndata++;
      i++;
   }
   BuildOld();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if we can build this for foreign class - do we have some rules
/// to do that.

Bool_t TStreamerInfo::BuildFor( const TClass *in_memory_cl )
{
   R__LOCKGUARD(gInterpreterMutex);

   if( !in_memory_cl || !in_memory_cl->GetSchemaRules() ) {
      return kFALSE;
   }

   auto rules = in_memory_cl->GetSchemaRules()->FindRules( GetName(), fOnFileClassVersion, fCheckSum );

   if( rules.empty() && !in_memory_cl->GetCollectionType() ) {
      Warning( "BuildFor", "The build of %s streamer info for %s has been requested, but no matching conversion rules were specified", GetName(), in_memory_cl->GetName() );
      return kFALSE;
   }

   fClass = const_cast<TClass*>(in_memory_cl);

   return kTRUE;
}


namespace {
////////////////////////////////////////////////////////////////////////////////
/// Helper function for BuildOld
   Bool_t ClassWasMovedToNamespace(TClass *oldClass, TClass *newClass)
   {
      // Returns true if oldClass is the same as newClass but newClass is in a
      // namespace (and oldClass was not in a namespace).

      if (oldClass == 0 || newClass == 0) return kFALSE;

      UInt_t newlen = strlen(newClass->GetName());
      UInt_t oldlen = strlen(oldClass->GetName());

      const char *oldname = oldClass->GetName();
      for (UInt_t i = oldlen, done = false, nest = 0; (i>0) && !done ; --i) {
         switch (oldClass->GetName()[i-1]) {
            case '>' : ++nest; break;
            case '<' : if (nest==0) return kFALSE; // the name is not well formed, give up.
                       --nest; break;
            case ':' : if (nest == 0) oldname= &(oldClass->GetName()[i]); done = kTRUE; break;
         }
      }
      oldlen = strlen(oldname);
      if (!(strlen(newClass->GetName()) > strlen(oldClass->GetName()))) {
         return kFALSE;
      }

      const char* newEnd = & (newClass->GetName()[newlen-oldlen]);

      if (0 != strcmp(newEnd, oldname)) {
         return kFALSE;
      }

      Int_t oldv = oldClass->GetStreamerInfo()->GetClassVersion();

      if (newClass->GetStreamerInfos() && oldv < newClass->GetStreamerInfos()->GetSize() && newClass->GetStreamerInfos()->At(oldv) && strcmp(newClass->GetStreamerInfos()->At(oldv)->GetName(), oldClass->GetName()) != 0) {
         // The new class has already a TStreamerInfo for the the same version as
         // the old class and this was not the result of an import.  So we do not
         // have a match
         return kFALSE;
      }
      return kTRUE;
   }

////////////////////////////////////////////////////////////////////////////////
/// Import the streamerInfo from oldClass to newClass.
///
/// In case of conflict, returns the version number of the StreamerInfo
/// with the conflict.
/// Return 0 in case of success
   Int_t ImportStreamerInfo(TClass *oldClass, TClass *newClass) {

      TIter next(oldClass->GetStreamerInfos());
      TStreamerInfo *info;
      while ((info = (TStreamerInfo*)next())) {
         info = (TStreamerInfo*)info->Clone();
         if (!info) {
            Error("ImportStreamerInfo","Unable to clone the StreamerInfo for %s.",(*next)->GetName());
         } else {
            info->SetClass(newClass);
            Int_t oldv = info->GetClassVersion();
            if (oldv > newClass->GetStreamerInfos()->GetSize() || newClass->GetStreamerInfos()->At(oldv) == 0) {
               // All is good.
               newClass->RegisterStreamerInfo(info);
            } else {
               // We verify that we are consistent and that
               //   newcl->GetStreamerInfos()->UncheckedAt(info->GetClassVersion)
               // is already the same as info.
               if (strcmp(newClass->GetStreamerInfos()->At(oldv)->GetName(),
                          oldClass->GetName()) != 0) {
                  // The existing StreamerInfo does not already come from OldClass.
                  // This is a real problem!
                  return oldv;
               }
            }
         }
      }
      return 0;
   }

   Bool_t ContainerMatchTClonesArray(TClass *newClass)
   {
      // Return true if newClass is a likely valid conversion from
      // a TClonesArray

      return newClass->GetCollectionProxy()
             && newClass->GetCollectionProxy()->GetValueClass()
             && !newClass->GetCollectionProxy()->HasPointers();
   }

   Bool_t CollectionMatch(const TClass *oldClass, const TClass* newClass)
   {
      // Return true if oldClass and newClass points to 2 compatible collection.
      // i.e. they contains the exact same type.

      TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();
      TVirtualCollectionProxy *newProxy = newClass->GetCollectionProxy();

      TClass *oldContent = oldProxy->GetValueClass();
      TClass *newContent = newProxy->GetValueClass();

      Bool_t contentMatch = kFALSE;
      if (oldContent) {
         if (oldContent == newContent) {
            contentMatch = kTRUE;
         } else if (newContent) {
            TString oldFlatContent( TMakeProject::UpdateAssociativeToVector(oldContent->GetName()) );
            TString newFlatContent( TMakeProject::UpdateAssociativeToVector(newContent->GetName()) );
            if (oldFlatContent == newFlatContent) {
               contentMatch = kTRUE;
            }
         } else {
            contentMatch = kFALSE;
         }
      } else {
         contentMatch = (newContent==0);
      }

      if (contentMatch) {
         if ((oldContent==0 && oldProxy->GetType() == newProxy->GetType())
             ||(oldContent && oldProxy->HasPointers() == newProxy->HasPointers())) {
            // We have compatibles collections (they have the same content)!
            return kTRUE;
         }
      }
      return kFALSE;
   }

   Bool_t CollectionMatchFloat16(const TClass *oldClass, const TClass* newClass)
   {
      // Return true if oldClass and newClass points to 2 compatible collection.
      // i.e. they contains the exact same type.

      TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();
      TVirtualCollectionProxy *newProxy = newClass->GetCollectionProxy();

      if (oldProxy->GetValueClass() == 0 && newProxy->GetValueClass() == 0
          && (oldProxy->GetType() == kFloat_t || oldProxy->GetType() == kFloat16_t)
          && (newProxy->GetType() == kFloat_t || newProxy->GetType() == kFloat16_t )) {
            // We have compatibles collections (they have the same content)!
         return (oldClass->GetCollectionType() == newClass->GetCollectionType());
      }
      return kFALSE;
   }

   Bool_t CollectionMatchDouble32(const TClass *oldClass, const TClass* newClass)
   {
      // Return true if oldClass and newClass points to 2 compatible collection.
      // i.e. they contains the exact same type.

      TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();
      TVirtualCollectionProxy *newProxy = newClass->GetCollectionProxy();

      if (oldProxy->GetValueClass() == 0 && newProxy->GetValueClass() == 0
          && (oldProxy->GetType() == kDouble_t || oldProxy->GetType() == kDouble32_t)
          && (newProxy->GetType() == kDouble_t || newProxy->GetType() == kDouble32_t )) {
            // We have compatibles collections (they have the same content)!
         return (oldClass->GetCollectionType() == newClass->GetCollectionType());
      }
      return kFALSE;
   }

   Bool_t CollectionMatchLong64(const TClass *oldClass, const TClass* newClass)
   {
      // Return true if oldClass and newClass points to 2 compatible collection.
      // i.e. they contains the exact same type.

      TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();
      TVirtualCollectionProxy *newProxy = newClass->GetCollectionProxy();

      if (oldProxy->GetValueClass() == 0 && newProxy->GetValueClass() == 0
          && (oldProxy->GetType() == kLong_t || oldProxy->GetType() == kLong64_t)
          && (newProxy->GetType() == kLong_t || newProxy->GetType() == kLong64_t )) {
         // We have compatibles collections (they have the same content)!
         return (oldClass->GetCollectionType() == newClass->GetCollectionType());
      }
      return kFALSE;
   }

   Bool_t CollectionMatchULong64(const TClass *oldClass, const TClass* newClass)
   {
      // Return true if oldClass and newClass points to 2 compatible collection.
      // i.e. they contains the exact same type.

      TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();
      TVirtualCollectionProxy *newProxy = newClass->GetCollectionProxy();

      if (oldProxy->GetValueClass() == 0 && newProxy->GetValueClass() == 0
          && (oldProxy->GetType() == kULong_t || oldProxy->GetType() == kULong64_t)
          && (newProxy->GetType() == kULong_t || newProxy->GetType() == kULong64_t )) {
         // We have compatibles collections (they have the same content)!
         return (oldClass->GetCollectionType() == newClass->GetCollectionType());
      }
      return kFALSE;
   }

   TClass *FindAlternate(TClass *context, const std::string &i_name, std::string& newName)
   {
      // Return a class whose has the name as oldClass and can be found
      // within the scope of the class 'context'.

      // First strip any 'const ' prefix or trailing '*'.
      std::string name(i_name);
      newName.clear();
      if (name.compare(0,6,"const ")==0) {
         newName = "const ";
         name.erase(0,6);
      }
      std::string suffix;
      UInt_t nstars = 0;
      while(name[name.length()-nstars-1]=='*') {
         ++nstars;
         suffix.append("*");
      }
      if (nstars) {
         name.erase(name.length()-nstars,nstars);
      }

      std::string alternate(context->GetName());
      alternate.append("::");
      alternate.append(name);

      TClass *altcl = TClass::GetClass(alternate.c_str(),/*load=*/ false,true);
      if (altcl) {
         newName.append(altcl->GetName());
         newName.append(suffix);
         return altcl;
      }

      size_t ctxt_cursor = strlen(context->GetName());
      for (size_t level = 0; ctxt_cursor != 0; --ctxt_cursor) {
         switch (context->GetName()[ctxt_cursor]) {
            case '<': --level; break;
            case '>': ++level; break;
            case ':': if (level == 0) {
               // we encountered a scope not within a template
               // parameter.
               alternate.clear();
               alternate.append(context->GetName(),ctxt_cursor+1);
               alternate.append(name);
               altcl = TClass::GetClass(alternate.c_str(),/*load=*/ false,true);
               if (altcl) {
                  newName.append(altcl->GetName());
                  newName.append(suffix);
                  return altcl;
               }
            }
         }
      }
      newName.clear();
      return 0;
   }

   TClass *FixCollectionV5(TClass *context, TClass *oldClass, TClass *newClass)
   {
      assert(oldClass->GetCollectionProxy() && newClass->GetCollectionProxy());

      TVirtualCollectionProxy *old = oldClass->GetCollectionProxy();
      TVirtualCollectionProxy *current = newClass->GetCollectionProxy();
      Int_t stlkind = old->GetCollectionType();

      if (stlkind == ROOT::kSTLmap || stlkind == ROOT::kSTLmultimap) {

         if (current->GetValueClass() == nullptr) {
            // This should really never happen (the content of map should always
            // be a pair and thus have a TClass ... so let's just give up ...
            // It actually happens in the case where one of the member is an
            // enum that is part of dictionary payload that is not yet
            // autoloaded.
            return nullptr;
         }
         TVirtualStreamerInfo *info = current->GetValueClass()->GetStreamerInfo();
         if (info->GetElements()->GetEntries() != 2) {
            return oldClass;
         }
         TStreamerElement *f = (TStreamerElement*) info->GetElements()->At(0);
         TStreamerElement *s = (TStreamerElement*) info->GetElements()->At(1);

         info = old->GetValueClass()->GetStreamerInfo();
         assert(info->GetElements()->GetEntries() == 2);
         TStreamerElement *of = (TStreamerElement*) info->GetElements()->At(0);
         TStreamerElement *os = (TStreamerElement*) info->GetElements()->At(1);

         TClass *firstNewCl  = f ? f->GetClass() : 0;
         TClass *secondNewCl = s ? s->GetClass() : 0;

         TClass *firstOldCl  = of ? of->GetClass() : 0;
         TClass *secondOldCl = os ? os->GetClass() : 0;

         if ((firstNewCl && !firstOldCl) || (secondNewCl && !secondOldCl))
         {
            std::vector<std::string> inside;
            int nestedLoc;
            TClassEdit::GetSplit( oldClass->GetName(), inside, nestedLoc, TClassEdit::kLong64 );

            TClass *firstAltCl = firstOldCl;
            TClass *secondAltCl = secondOldCl;
            std::string firstNewName;
            std::string secondNewName;
            if (firstNewCl && !firstOldCl) {
               firstAltCl = FindAlternate(context, inside[1], firstNewName);
            } else if (firstAltCl) {
               firstNewName = firstAltCl->GetName();
            } else {
               firstNewName = inside[1];
            }
            if (secondNewCl && !secondOldCl) {
               secondAltCl = FindAlternate(context, inside[2], secondNewName);
            } else if (secondAltCl) {
               secondNewName = secondAltCl->GetName();
            } else {
               secondNewName  = inside[2];
            }
            if ((firstNewCl && firstAltCl != firstOldCl) ||
                (secondNewCl && secondAltCl != secondOldCl) ) {

               // Need to produce new name.
               std::string alternate = inside[0];
               alternate.append("<");
               alternate.append(firstNewName);
               alternate.append(",");
               alternate.append(secondNewName);
               // We are intentionally dropping any further arguments,
               // they would be using the wrong typename and would also be
               // somewhat superflous since this is for the old layout.
               if (alternate[alternate.length()-1]=='>') {
                  alternate.append(" ");
               }
               alternate.append(">");
               return TClass::GetClass(alternate.c_str(),true,true);
            }
         }

      } else if (current->GetValueClass() && !old->GetValueClass()
          && old->GetType() == kInt_t) {

         // The old CollectionProxy claims it contains int (or enums) while
         // the new one claims to contain a class.  It is likely that we have
         // in the collection name a class (typedef) name that is missing its
         // scope.  Let's try to check.

         std::vector<std::string> inside;
         int nestedLoc;
         TClassEdit::GetSplit( oldClass->GetName(), inside, nestedLoc, TClassEdit::kLong64 );

         // Now let's if we can find this missing type.
         std::string newName;
         TClass *altcl = FindAlternate(context, inside[1], newName);

         if (altcl) {
            std::string alternate = inside[0];
            alternate.append("<");
            alternate.append(newName);
            // We are intentionally dropping any further arguments,
            // they would be using the wrong typename and would also be
            // somewhat superflous since this is for the old layout.
            if (alternate[alternate.length()-1]=='>') {
               alternate.append(" ");
            }
            alternate.append(">");
            return TClass::GetClass(alternate.c_str(),true,true);
         }
      }
      return 0;
   }

   // Makes sure kBuildOldUsed set once BuildOld finishes
   struct TBuildOldGuard {
      TBuildOldGuard(TStreamerInfo* info): fInfo(info) {
         fInfo->SetBit(TStreamerInfo::kBuildRunning);
      }
      ~TBuildOldGuard() {
         fInfo->ResetBit(TStreamerInfo::kBuildRunning);
         fInfo->SetBit(TStreamerInfo::kBuildOldUsed);
      }
      TStreamerInfo* fInfo;
   };
}

////////////////////////////////////////////////////////////////////////////////
/// rebuild the TStreamerInfo structure

void TStreamerInfo::BuildOld()
{
   R__LOCKGUARD(gInterpreterMutex);

   if ( TestBit(kBuildOldUsed) ) return;

   // Are we recursing on ourself?
   if (TestBit(TStreamerInfo::kBuildRunning)) return;

   // This is used to avoid unwanted recursive call to Build and make sure
   // that we record the execution of BuildOld.
   TBuildOldGuard buildOldGuard(this);

   if (gDebug > 0) {
      printf("\n====>Rebuilding TStreamerInfo for class: %s, version: %d\n", GetName(), fClassVersion);
   }

   Bool_t wasCompiled = IsCompiled();

   if (fClass->GetClassVersion() == fClassVersion) {
      if (!fClass->HasInterpreterInfo() || fClass->GetCollectionType() || TClassEdit::IsSTLBitset(GetName()))
      {
         // Handle emulated classes and STL containers specially.
         // in this case BuildRealData would call BuildOld for this same
         // TStreamerInfo to be able to build the real data on it.
      } else {
         fClass->BuildRealData();
      }
   }
   else {
      // This is to support the following case
      //  Shared library: Event v2
      //  calling cl->GetStreamerInfo(1)->BuildOld();  (or equivalent)
      //  which calls cl->BuildReadData()
      //  which set fRealData to some value
      //  then call Event()
      //  which call cl->GetStreamerInfo()
      //  which call cl->BuildRealData();
      //  which returns immediately (upon seeing fRealData!=0)
      //  then the main StreamerInfo build using the partial content of fRealData
      //  then BuildRealData returns
      //  then GetStreamerInfo() returns
      //  then Event() returns
      //  then fRealData is finished being populated
      //  then this function continue,
      //  then it uses the main streamerInfo
      //  .... which is incomplete.
      //
      //  Instead we force the creation of the main streamerInfo object
      //  before the creation of fRealData.
      fClass->GetStreamerInfo();
   }

   TIter next(fElements);
   TStreamerElement* element;
   Int_t offset = 0;
   TMemberStreamer* streamer = 0;

   constexpr size_t kSizeOfPtr = sizeof(void*);

   int nBaze = 0;

   if ((fElements->GetEntries() == 1) && !strcmp(fElements->At(0)->GetName(), "This")) {
      if (fClass->GetCollectionProxy())  {
         element = (TStreamerElement*)next();
         element->SetNewType( element->GetType() );
         element->SetNewClass( fClass );
      } else if (((TStreamerElement*)fElements->At(0))->GetType() == TStreamerInfo::kSTL &&
                 strcmp( ((TStreamerElement*)fElements->At(0))->GetTypeName(),GetName()) != 0) {
         // We have a collection that was proxied but does not have a collection proxy,
         // let's put one in place just for fun ... humm however we have no clue what is the value
         // type ....

         // For now wild guess ....

      }
   }

   TClass *allocClass = 0;
   TStreamerInfo *infoalloc = 0;

   //---------------------------------------------------------------------------
   // Get schema rules for this class
   /////////////////////////////////////////////////////////////////////////////

   ROOT::TSchemaRuleSet::TMatches rules;
   const ROOT::TSchemaRuleSet* ruleSet = fClass->GetSchemaRules();

   if (ruleSet) rules = ruleSet->FindRules( GetName(), fOnFileClassVersion, fCheckSum );

   Bool_t shouldHaveInfoLoc = fClass->TestBit(TClass::kIsEmulation) && !TClassEdit::IsStdClass(fClass->GetName());
   Int_t virtualInfoLocAlloc = 0;
   fNVirtualInfoLoc = 0;
   delete [] fVirtualInfoLoc;
   fVirtualInfoLoc = 0;

   while ((element = (TStreamerElement*) next())) {
      if (element->IsA()==TStreamerArtificial::Class()
          || element->TestBit(TStreamerElement::kCache) )
      {
         // Prevent BuildOld from modifying existing ArtificialElement (We need to review when and why BuildOld
         // needs to be re-run; it might be needed if the 'current' class change (for example from being an onfile
         // version to being a version loaded from a shared library) and we thus may have to remove the artifical
         // element at the beginning of BuildOld)

         continue;
      };

      element->SetNewType(element->GetType());
      if (element->IsBase()) {
         //---------------------------------------------------------------------
         // Dealing with nonSTL bases
         ///////////////////////////////////////////////////////////////////////

         if (element->IsA() == TStreamerBase::Class()) {
            TStreamerBase* base = (TStreamerBase*) element;
#if defined(PROPER_IMPLEMEMANTION_OF_BASE_CLASS_RENAMING)
            TClassRef baseclass =  fClass->GetBaseClass( base->GetName() );
#else
            // Currently the base class renaming does not work, so we use the old
            // version of the code which essentially disable the next if(!baseclass ..
            // statement.
            // During the TStreamerElement's Init an emulated TClass might be replaced
            // by one from the dictionary, we use a TClassRef to be informed of the change.
            TClassRef baseclass =  base->GetClassPointer();
#endif

            //------------------------------------------------------------------
            // We do not have this base class - check if we're renaming
            ////////////////////////////////////////////////////////////////////

            if( !baseclass && !fClass->TestBit( TClass::kIsEmulation ) ) {
               const ROOT::TSchemaRule* rule = (rules ? rules.GetRuleWithSource( base->GetName() ) : 0);

               //---------------------------------------------------------------
               // No renaming, sorry
               /////////////////////////////////////////////////////////////////

               if( !rule ) {
                  Error("BuildOld", "Could not find base class: %s for %s and could not find any matching rename rule\n", base->GetName(), GetName());
                  continue;
               }

               //----------------------------------------------------------------
               // Find a new target class
               /////////////////////////////////////////////////////////////////

               const TObjArray* targets = rule->GetTarget();
               if( !targets ) {
                  Error("BuildOld", "Could not find base class: %s for %s, renaming rule was found but is malformed\n", base->GetName(), GetName());
               }
               TString newBaseClass = ((TObjString*)targets->At(0))->GetString();
               baseclass = TClass::GetClass( newBaseClass );
               base->SetNewBaseClass( baseclass );
            }
            //-------------------------------------------------------------------
            // No base class in emulated mode
            ////////////////////////////////////////////////////////////////////

            else if( !baseclass ) {
               baseclass = base->GetClassPointer();
               if (!baseclass) {
                  Warning("BuildOld", "Missing base class: %s skipped", base->GetName());
                  // FIXME: Why is the version number 1 here? Answer: because we don't know any better at this point
                  baseclass = new TClass(element->GetName(), 1, 0, 0, -1, -1);
                  element->Update(0, baseclass);
               }
            }
            baseclass->BuildRealData();

            // Calculate the offset using the 'real' base class name (as opposed to the
            // '@@emulated' in the case of the emulation of an abstract base class.
            Int_t baseOffset = fClass->GetBaseClassOffset(baseclass);

            // Deal with potential schema evolution (renaming) of the base class.
            if (baseOffset < 0) {

               // See if this base element can be converted into one of
               // the existing base class.
               TList* listOfBases = fClass->GetListOfBases();
               if (listOfBases) {
                  TBaseClass* bc = 0;
                  TIter nextBC(fClass->GetListOfBases());
                  while ((bc = (TBaseClass*) nextBC())) {
                     TClass *in_memory_bcl = bc->GetClassPointer();
                     if (in_memory_bcl && in_memory_bcl->GetSchemaRules()) {
                        auto baserule = in_memory_bcl->GetSchemaRules()->FindRules( base->GetName(), base->GetBaseVersion(), base->GetBaseCheckSum() );
                        if (!baserule.empty()) {
                           base->SetNewBaseClass(in_memory_bcl);
                           baseOffset = bc->GetDelta();

                        }
                     }
                  }
               }
            }
            // We need to initialize the element now, as we need the
            // correct StreamerInfo next.
            element->Init(this);

            // Force the StreamerInfo "Compilation" of the base classes first. This is necessary in
            // case the base class contains a member used as an array dimension in the derived classes.
            TStreamerInfo* infobase;
            if (fClass->TestBit(TClass::kIsEmulation) && (baseclass->Property() & kIsAbstract)) {
               Int_t version = base->GetBaseVersion();
               if (version >= 0 || base->GetBaseCheckSum() == 0) {
                  infobase = (TStreamerInfo*)baseclass->GetStreamerInfoAbstractEmulated(version);
               } else {
                  infobase = (TStreamerInfo*)baseclass->FindStreamerInfoAbstractEmulated(base->GetBaseCheckSum());
               }
               if (infobase) baseclass = infobase->GetClass();
            }
            else {
               infobase = (TStreamerInfo*)base->GetBaseStreamerInfo();
            }

            if (infobase && infobase->fComp == 0) {
               infobase->BuildOld();
            }

            if (infobase && shouldHaveInfoLoc && baseclass->TestBit(TClass::kIsEmulation) ) {
               if ( (fNVirtualInfoLoc + infobase->fNVirtualInfoLoc) > virtualInfoLocAlloc ) {
                  ULong_t *store = fVirtualInfoLoc;
                  virtualInfoLocAlloc = 16 * ( (fNVirtualInfoLoc + infobase->fNVirtualInfoLoc) / 16 + 1);
                  fVirtualInfoLoc = new ULong_t[virtualInfoLocAlloc];
                  if (store) {
                     memcpy(fVirtualInfoLoc, store, sizeof(ULong_t)*fNVirtualInfoLoc);
                     delete [] store;
                  }
               }
               for (int nloc = 0; nloc < infobase->fNVirtualInfoLoc; ++nloc) {
                  fVirtualInfoLoc[ fNVirtualInfoLoc + nloc ] = baseOffset + infobase->fVirtualInfoLoc[nloc];
               }
               fNVirtualInfoLoc += infobase->fNVirtualInfoLoc;
            }


            {
               if (baseOffset < 0) {
                  element->SetNewType(-1);
               }
            }
            element->SetOffset(baseOffset);
            offset += baseclass->Size();

            continue;
         } else {
            // Not a base elem but still base, string or STL as a base
            nBaze++;
            TList* listOfBases = fClass->GetListOfBases();
            Int_t baseOffset = -1;
            Int_t asize = 0;
            if (listOfBases) {
               // Do a search for the classname and some of its alternatives spelling.

               TBaseClass* bc = 0;
               TIter nextBC(fClass->GetListOfBases());
               while ((bc = (TBaseClass*) nextBC())) {
                  if (strchr(bc->GetName(), '<') || !strcmp(bc->GetName(),"string")) {
                     TString bcName(TClassEdit::ShortType(bc->GetName(), TClassEdit::kDropStlDefault).c_str());
                     TString elName(TClassEdit::ShortType(element->GetTypeName(), TClassEdit::kDropStlDefault).c_str());
                     if (bcName == elName) {
                        break;
                     }
                  }
               }

               if (!bc) {
                  // Error("BuildOld", "Could not find STL base class: %s for %s\n", element->GetName(), GetName());
                  offset = kMissing;
                  element->SetOffset(kMissing);
                  element->SetNewType(-1);
                  continue;
               } else if (bc->GetClassPointer()->GetCollectionProxy()
                          && !bc->GetClassPointer()->IsLoaded()
                          && bc->GetClassPointer()->GetCollectionProxy()->GetCollectionType() != ROOT::kSTLvector) {
                  Error("BuildOld","The class \"%s\" is compiled and its base class \"%s\" is a collection and we do not have a dictionary for it, we will not be able to read or write this base class.",GetName(),bc->GetName());
                  offset = kMissing;
                  element->SetOffset(kMissing);
                  element->SetNewType(-1);
                  continue;
               }
               baseOffset = bc->GetDelta();
               asize = bc->GetClassPointer()->Size();

            } else if (fClass->TestBit( TClass::kIsEmulation )) {
               // Do a search for the classname and some of its alternatives spelling.

               TStreamerInfo* newInfo = (TStreamerInfo*) fClass->GetStreamerInfos()->At(fClass->GetClassVersion());
               if (newInfo == this) {
                  baseOffset = offset;
                  asize = element->GetSize();
               } else if (newInfo) {
                  TIter newElems( newInfo->GetElements() );
                  TStreamerElement *newElement;
                  while( (newElement = (TStreamerElement*)newElems()) ) {
                     const char *newElName = newElement->GetName();
                     if (newElement->IsBase() && (strchr(newElName,'<') || !strcmp(newElName,"string")) ) {
                        TString bcName(TClassEdit::ShortType(newElName, TClassEdit::kDropStlDefault).c_str());
                        TString elName(TClassEdit::ShortType(element->GetTypeName(), TClassEdit::kDropStlDefault).c_str());
                        if (bcName == elName) {
                           break;
                        }
                     }
                  }
                  if (!newElement) {
                     Error("BuildOld", "Could not find STL base class: %s for %s\n", element->GetName(), GetName());
                     continue;
                  }
                  baseOffset = newElement->GetOffset();
                  asize = newElement->GetSize();
               }
            }
            if (baseOffset == -1) {
               TClass* cb = element->GetClassPointer();
               if (!cb) {
                  element->SetNewType(-1);
                  continue;
               }
               asize = cb->Size();
               baseOffset = fClass->GetBaseClassOffset(cb);
            }

            //  we know how to read but do we know where to read?
            if (baseOffset < 0) {
               element->SetNewType(-1);
               continue;
            }
            element->SetOffset(baseOffset);
            offset += asize;
            element->Init(this);
            continue;
         } // if element is of type TStreamerBase or not.
      } // if (element->IsBase())

      // If we get here, this means that we looked at all the base classes.
      if (shouldHaveInfoLoc && fNVirtualInfoLoc==0) {
         fNVirtualInfoLoc = 1;
         fVirtualInfoLoc = new ULong_t[1]; // To allow for a single delete statement.
         fVirtualInfoLoc[0] = offset;
         offset += sizeof(TStreamerInfo*);
      }

      TDataMember* dm = 0;

      std::string typeNameBuf;
      const char* dmType = nullptr;
      Bool_t dmIsPtr = false;
      TDataType* dt(nullptr);
      Int_t ndim = 0 ; //dm->GetArrayDim();
      std::array<Int_t, 5> maxIndices; // 5 is the maximum supported in TStreamerElement::SetMaxIndex
      Bool_t isStdArray(kFALSE);

      // First set the offset and sizes.
      if (fClass->GetState() <= TClass::kEmulated) {
         // Note the initilization in this case are
         // delayed until __after__ the schema evolution
         // section, just in case the info has changed.

         // We are in the emulated case
         streamer = 0;
         element->Init(this);
      } else {
         // The class is known to Cling (and thus is not emulated)
         // and we need to use the real offsets.
         // However we may not have a 'proper' TClass for it
         // (in which case IsLoaded will be false and GetImplFileLine will be -1)

         // First look for the data member in the current class
         dm = (TDataMember*) fClass->GetListOfDataMembers()->FindObject(element->GetName());
         if (dm && dm->IsPersistent()) {
            fClass->BuildRealData();
            streamer = 0;
            offset = GetDataMemberOffset(dm, streamer);
            element->SetOffset(offset);
            element->Init(this);

            // Treat unique pointers and std arrays
            dmType = dm->GetTypeName();
            dmIsPtr = dm->IsaPointer();
            Bool_t nameChanged;
            typeNameBuf = TClassEdit::GetNameForIO(dmType, TClassEdit::EModType::kNone, &nameChanged);
            if (nameChanged) {
               dmIsPtr = TClassEdit::IsUniquePtr(dmType);
               dmType = typeNameBuf.c_str();
            }
            if ((isStdArray = TClassEdit::IsStdArray(dmType))){ // We tackle the std array case
               TClassEdit::GetStdArrayProperties(dmType,
                                                 typeNameBuf,
                                                 maxIndices,
                                                 ndim);
               dmType = typeNameBuf.c_str();
               dt = gROOT->GetType(dmType);
            }

            // We have a loaded class, let's make sure that if we have a collection
            // it is also loaded.
            TString dmClassName = TClassEdit::ShortType(dmType,TClassEdit::kDropStlDefault).c_str();
            dmClassName = dmClassName.Strip(TString::kTrailing, '*');
            if (dmClassName.Index("const ")==0) dmClassName.Remove(0,6);
            TClass *elemDm = ! (dt || dm->IsBasic()) ? TClass::GetClass(dmClassName.Data()) : 0;
            if (elemDm && elemDm->GetCollectionProxy()
                && !elemDm->IsLoaded()
                && elemDm->GetCollectionProxy()->GetCollectionType() != ROOT::kSTLvector) {
               Error("BuildOld","The class \"%s\" is compiled and for its data member \"%s\", we do not have a dictionary for the collection \"%s\", we will not be able to read or write this data member.",GetName(),dm->GetName(),elemDm->GetName());
               offset = kMissing;
               element->SetOffset(kMissing);
               element->SetNewType(-1);
            }
            element->SetStreamer(streamer);
            int narr = element->GetArrayLength();
            if (!narr) {
               narr = 1;
            }
            int dsize = dm->GetUnitSize();
            element->SetSize(dsize*narr);
         } else {
            // We did not find it, let's look for it in the base classes via TRealData
            TRealData* rd = fClass->GetRealData(element->GetName());
            if (rd && rd->GetDataMember()) {
               element->SetOffset(rd->GetThisOffset());
               element->Init(this);
               dm = rd->GetDataMember();
               dmType = dm->GetTypeName();
               dmIsPtr = dm->IsaPointer();
               int narr = element->GetArrayLength();
               if (!narr) {
                  narr = 1;
               }
               int dsize = dm->GetUnitSize();
               element->SetSize(dsize*narr);
            }
         }
      } // Class corresponding to StreamerInfo is emulated or not.

      // Now let's deal with Schema evolution
      Int_t newType = -1;
      TClassRef newClass;

      if (dm && dm->IsPersistent()) {
         auto theType = isStdArray ? dt : dm->GetDataType();
         if (theType) {
            Bool_t isArray = isStdArray || element->GetArrayLength() >= 1;
            Bool_t hasCount = element->HasCounter();
            // data member is a basic type
            if ((fClass == TObject::Class()) && !strcmp(dm->GetName(), "fBits")) {
               //printf("found fBits, changing dtype from %d to 15\n", dtype);
               newType = kBits;
            } else {
               // All the values of EDataType have the same semantic in EReadWrite
               newType = (EReadWrite)theType->GetType();
            }
            if ((newType == ::kChar_t) && dmIsPtr && !isArray && !hasCount) {
               newType = ::kCharStar;
            } else if (dmIsPtr) {
               newType += kOffsetP;
            } else if (isArray) {
               newType += kOffsetL;
            }
         }
         if (newType == -1) {
            newClass = TClass::GetClass(dmType);
         }
      } else {
         // Either the class is not loaded or the data member is gone
         if (!fClass->IsLoaded()) {
            TStreamerInfo* newInfo = (TStreamerInfo*) fClass->GetStreamerInfos()->At(fClass->GetClassVersion());
            if (newInfo && (newInfo != this)) {
               TStreamerElement* newElems = (TStreamerElement*) newInfo->GetElements()->FindObject(element->GetName());
               newClass = newElems ?  newElems->GetClassPointer() : 0;
               if (newClass == 0) {
                  newType = newElems ? newElems->GetType() : -1;
                  if (!(newType < kObject)) {
                     // sanity check.
                     newType = -1;
                  }
               }
            } else {
               newClass = element->GetClassPointer();
               if (newClass.GetClass() == 0) {
                  newType = element->GetType();
                  if (!(newType < kObject)) {
                     // sanity check.
                     newType = -1;
                  }
               }
            }
         }
      }

      if (newType > 0) {
         // Case of a numerical type
         if (element->GetType() >= TStreamerInfo::kObject) {
            // Old type was not a numerical type.
            element->SetNewType(-2);
         } else if (element->GetType() != newType) {
            element->SetNewType(newType);
            if (gDebug > 0) {
               // coverity[mixed_enums] - All the values of EDataType have the same semantic in EReadWrite
               Info("BuildOld", "element: %s %s::%s has new type: %s/%d", element->GetTypeName(), GetName(), element->GetName(), dm ? dm->GetFullTypeName() : TDataType::GetTypeName((EDataType)newType), newType);
            }
         }
      } else if (newClass.GetClass()) {
         // Sometime BuildOld is called again.
         // In that case we might already have fix up the streamer element.
         // So we need to go back to the original information!
         newClass.Reset();
         TClass* oldClass = TClass::GetClass(TClassEdit::ShortType(element->GetTypeName(), TClassEdit::kDropTrailStar).c_str());
         if (oldClass == newClass.GetClass()) {
            // Nothing to do, also in the unique_ptr case :)
         } else if (ClassWasMovedToNamespace(oldClass, newClass.GetClass())) {
            Int_t oldv;
            if (0 != (oldv = ImportStreamerInfo(oldClass, newClass.GetClass()))) {
                Warning("BuildOld", "Can not properly load the TStreamerInfo from %s into %s due to a conflict for the class version %d", oldClass->GetName(), newClass->GetName(), oldv);
            } else {
               element->SetTypeName(newClass->GetName());
               if (gDebug > 0) {
                  Warning("BuildOld", "element: %s::%s %s has new type %s", GetName(), element->GetTypeName(), element->GetName(), newClass->GetName());
               }
            }
         } else if (oldClass == TClonesArray::Class()) {
            if (ContainerMatchTClonesArray(newClass.GetClass())) {
               Int_t elemType = element->GetType();
               Bool_t isPrealloc = (elemType == kObjectp) || (elemType == kAnyp) || (elemType == (kObjectp + kOffsetL)) || (elemType == (kAnyp + kOffsetL));
               element->Update(oldClass, newClass.GetClass());
               TVirtualCollectionProxy *cp = newClass->GetCollectionProxy();
               TConvertClonesArrayToProxy *ms = new TConvertClonesArrayToProxy(cp, element->IsaPointer(), isPrealloc);
               element->SetStreamer(ms);

               // When the type is kObject, the TObject::Streamer is used instead
               // of the TStreamerElement's streamer.  So let force the usage
               // of our streamer
               if (element->GetType() == kObject) {
                  element->SetNewType(kAny);
                  element->SetType(kAny);
               }
               if (gDebug > 0) {
                  Warning("BuildOld","element: %s::%s %s has new type %s", GetName(), element->GetTypeName(), element->GetName(), newClass->GetName());
               }
            } else {
               element->SetNewType(-2);
            }
         } else if (oldClass && oldClass->GetCollectionProxy() && newClass->GetCollectionProxy()) {
            {
               TClass *oldFixedClass = FixCollectionV5(GetClass(),oldClass,newClass);
               if (oldFixedClass && oldFixedClass != oldClass) {
                  element->Update(oldClass,oldFixedClass);
                  oldClass = oldFixedClass;
               }
            }
            if (CollectionMatch(oldClass, newClass)) {
               Int_t oldkind = oldClass->GetCollectionType();
               Int_t newkind = newClass->GetCollectionType();

               if ( (oldkind==ROOT::kSTLmap || oldkind==ROOT::kSTLmultimap) &&
                    (newkind!=ROOT::kSTLmap && newkind!=ROOT::kSTLmultimap) ) {

                  Int_t elemType = element->GetType();
                  Bool_t isPrealloc = (elemType == kObjectp) || (elemType == kAnyp) || (elemType == (kObjectp + kOffsetL)) || (elemType == (kAnyp + kOffsetL));

                  TClassStreamer *streamer2 = newClass->GetStreamer();
                  if (streamer2) {
                     TConvertMapToProxy *ms = new TConvertMapToProxy(streamer2, element->IsaPointer(), isPrealloc);
                     if (ms && ms->IsValid()) {
                        element->SetStreamer(ms);
                        switch( element->GetType() ) {
                           //case TStreamerInfo::kSTLvarp:           // Variable size array of STL containers.
                        case TStreamerInfo::kSTLp:                // Pointer to container with no virtual table (stl) and no comment
                        case TStreamerInfo::kSTLp + TStreamerInfo::kOffsetL:     // array of pointers to container with no virtual table (stl) and no comment
                           element->SetNewType(-2);
                           break;
                        case TStreamerInfo::kSTL:             // container with no virtual table (stl) and no comment
                        case TStreamerInfo::kSTL + TStreamerInfo::kOffsetL:  // array of containers with no virtual table (stl) and no comment
                           break;
                        }
                     } else {
                        delete ms;
                     }
                  }
                  element->Update(oldClass, newClass.GetClass());

               } else if ( (newkind==ROOT::kSTLmap || newkind==ROOT::kSTLmultimap) &&
                           (oldkind!=ROOT::kSTLmap && oldkind!=ROOT::kSTLmultimap) ) {
                  element->SetNewType(-2);
               } else {
                  element->Update(oldClass, newClass.GetClass());
               }
               // Is this needed ? : element->SetSTLtype(newelement->GetSTLtype());
               if (gDebug > 0) {
                  Warning("BuildOld","element: %s::%s %s has new type %s", GetName(), element->GetTypeName(), element->GetName(), newClass->GetName());
               }
            } else if (CollectionMatchFloat16(oldClass,newClass)) {
               // Actually nothing to do, since both are the same collection of double in memory.
            } else if (CollectionMatchDouble32(oldClass,newClass)) {
               // Actually nothing to do, since both are the same collection of double in memory.
            } else if (CollectionMatchLong64(oldClass,newClass)) {
               // Not much to do since both are the same collection of 8 bits entities on file.
               element->Update(oldClass, newClass.GetClass());
            } else if (CollectionMatchULong64(oldClass,newClass)) {
               // Not much to do since both are the same collection of 8 bits unsigned entities on file
               element->Update(oldClass, newClass.GetClass());
            } else if (newClass->GetSchemaRules()->HasRuleWithSourceClass( oldClass->GetName() )) {
               //------------------------------------------------------------------------
               // We can convert one type to another (at least for some of the versions).
               /////////////////////////////////////////////////////////////////

               element->SetNewClass( newClass );
            } else {
               element->SetNewType(-2);
            }

         } else if(oldClass &&
                   newClass.GetClass() &&
                   newClass->GetSchemaRules() &&
                   newClass->GetSchemaRules()->HasRuleWithSourceClass( oldClass->GetName() ) ) {
            //------------------------------------------------------------------------
            // We can convert one type to another (at least for some of the versions).
            ////////////////////////////////////////////////////////////////////

            element->SetNewClass( newClass );
         } else {
            element->SetNewType(-2);
         }
         // Humm we still need to make sure we have the same 'type' (pointer, embedded object, array, etc..)
         Bool_t cannotConvert = kFALSE;
         if (element->GetNewType() != -2) {
            if (dm) {
               if (dmIsPtr) {
                  if (strncmp(dm->GetTitle(),"->",2)==0) {
                     // We are fine, nothing to do.
                     if (newClass->IsTObject()) {
                        newType = kObjectp;
                     } else if (newClass->GetCollectionProxy()) {
                        newType = kSTLp;
                     } else {
                        newType = kAnyp;
                     }
                  } else {
                     if (TClass::GetClass(dm->GetTypeName())->IsTObject()) {
                        newType = kObjectP;
                     } else if (newClass->GetCollectionProxy()) {
                        newType = kSTLp;
                     } else {
                        newType = kAnyP;
                     }
                  }
               } else {
                  if (newClass->GetCollectionProxy()) {
                     newType = kSTL;
                  } else if (newClass == TString::Class()) {
                     newType = kTString;
                  } else if (newClass == TObject::Class()) {
                     newType = kTObject;
                  } else if (newClass == TNamed::Class()) {
                     newType = kTNamed;
                  } else if (newClass->IsTObject()) {
                     newType = kObject;
                  } else {
                     newType = kAny;
                  }
               }
               if ((!dmIsPtr || newType==kSTLp) && (isStdArray ? ndim : dm->GetArrayDim()) > 0) {
                  newType += kOffsetL;
               }
            } else if (!fClass->IsLoaded()) {
               TStreamerInfo* newInfo = (TStreamerInfo*) fClass->GetStreamerInfos()->At(fClass->GetClassVersion());
               if (newInfo && (newInfo != this)) {
                  TStreamerElement* newElems = (TStreamerElement*) newInfo->GetElements()->FindObject(element->GetName());
                  if (newElems) {
                     newType = newElems->GetType();
                  }
               } else {
                  newType = element->GetType();
               }
            }
            if (element->GetType() == kSTL
                || ((element->GetType() == kObject || element->GetType() == kAny || element->GetType() == kObjectp || element->GetType() == kAnyp)
                    && oldClass == TClonesArray::Class()))
            {
               cannotConvert = (newType != kSTL && newType != kObject && newType != kAny && newType != kSTLp && newType != kObjectp && newType != kAnyp);

            } else if (element->GetType() == kSTLp  || ((element->GetType() == kObjectP || element->GetType() == kAnyP) && oldClass == TClonesArray::Class()) )
            {
               cannotConvert = (newType != kSTL && newType != kObject && newType != kAny && newType != kSTLp && newType != kObjectP && newType != kAnyP);

            } else if (element->GetType() == kSTL + kOffsetL
                || ((element->GetType() == kObject + kOffsetL|| element->GetType() == kAny + kOffsetL|| element->GetType() == kObjectp+ kOffsetL || element->GetType() == kAnyp+ kOffsetL)
                    && oldClass == TClonesArray::Class()))
            {
               cannotConvert = (newType != kSTL + kOffsetL && newType != kObject+ kOffsetL && newType != kAny+ kOffsetL && newType != kSTLp+ kOffsetL && newType != kObjectp+ kOffsetL && newType != kAnyp+ kOffsetL);

            } else if (element->GetType() == kSTLp + kOffsetL || ((element->GetType() == kObjectP+ kOffsetL || element->GetType() == kAnyP+ kOffsetL) && oldClass == TClonesArray::Class()) )
            {
               cannotConvert = (newType != kSTL+ kOffsetL && newType != kObject+ kOffsetL && newType != kAny+ kOffsetL && newType != kSTLp + kOffsetL&& newType != kObjectP+ kOffsetL && newType != kAnyP+ kOffsetL);

            } else if ((element->GetType() == kObjectp || element->GetType() == kAnyp
                 || element->GetType() == kObject || element->GetType() == kAny
                 || element->GetType() == kTObject || element->GetType() == kTNamed || element->GetType() == kTString )) {
               // We had Type* ... ; //-> or Type ...;
               // this is completely compatible with the same and with a embedded object.
               if (newType != -1) {
                  if (newType == kObjectp || newType == kAnyp
                      || newType == kObject || newType == kAny
                      || newType == kTObject || newType == kTNamed || newType == kTString) {
                     // We are fine, no transformation to make
                     element->SetNewType(newType);
                  } else {
                     // We do not support this yet.
                     cannotConvert = kTRUE;
                  }
               } else {
                  // We have no clue
                  printf("%s We have no clue\n", dm->GetName());
                  cannotConvert = kTRUE;
               }
            } else if (element->GetType() == kObjectP || element->GetType() == kAnyP) {
               if (newType != -1) {
                  if (newType == kObjectP || newType == kAnyP ) {
                     // nothing to do}
                  } else {
                     cannotConvert = kTRUE;
                  }
               } else {
                  // We have no clue
                  cannotConvert = kTRUE;
               }
            }
         }
         if (cannotConvert) {
            element->SetNewType(-2);
            if (gDebug > 0) {
               // coverity[mixed_enums] - All the values of EDataType have the same semantic in EReadWrite
               Info("BuildOld", "element: %s %s::%s has new type: %s/%d", element->GetTypeName(), GetName(), element->GetName(), dm ? dm->GetFullTypeName() : TDataType::GetTypeName((EDataType)newType), newType);
            }
         }
      } else {
         element->SetNewType(-1);
         offset = kMissing;
         element->SetOffset(kMissing);
      }

      if (offset != kMissing && fClass->GetState() <= TClass::kEmulated && !fClass->fIsSyntheticPair) {
         // Note the initialization in this case are
         // delayed until __after__ the schema evolution
         // section, just in case the info has changed.

         // The class is NOT known to Cling, i.e. is emulated,
         // and we need to use the calculated offset.

         Int_t asize;
         if (element->GetType() == TStreamerInfo::kSTL &&
             strcmp(element->GetName(),"This") == 0 &&
             strcmp(element->GetTypeName(),GetName()) == 0 &&
             !fClass->GetCollectionProxy()) {
            // Humm .. we are missing the collection Proxy
            // for a proxied (custom) collection ... avoid
            // an infinite recursion and take a wild guess
            asize = sizeof(std::vector<int>);
         } else {
            // Regular case
            asize = element->GetSize();
         }
         // align the non-basic data types (required on alpha and IRIX!!)
         if ((offset % kSizeOfPtr) != 0) {
            offset = offset - (offset % kSizeOfPtr) + kSizeOfPtr;
         }
         element->SetOffset(offset);
         offset += asize;
      }

      if (!wasCompiled && rules) {
         if (rules.HasRuleWithSource( element->GetName(), kTRUE ) ) {

            if (allocClass == 0) {
               infoalloc  = (TStreamerInfo *)Clone(TString::Format("%s@@%d",GetName(),GetOnFileClassVersion()));
               if (!infoalloc) {
                  Error("BuildOld","Unable to create the StreamerInfo for %s.",TString::Format("%s@@%d",GetName(),GetOnFileClassVersion()).Data());
               } else {
                  infoalloc->SetBit(kBuildOldUsed,false);
                  infoalloc->BuildCheck();
                  infoalloc->BuildOld();
                  allocClass = infoalloc->GetClass();
               }
            }

            // Now that we are caching the unconverted element, we do not assign it to the real type even if we could have!
            if (element->GetNewType()>0 /* intentionally not including base class for now */
                && !rules.HasRuleWithTarget( element->GetName(), kTRUE ) ) {

               TStreamerElement *copy = (TStreamerElement*)element->Clone();
               R__TObjArray_InsertBefore( fElements, copy, element );
               next(); // move the cursor passed the insert object.
               copy->SetBit(TStreamerElement::kRepeat);
               element = copy;

               // Warning("BuildOld","%s::%s is not set from the version %d of %s (You must add a rule for it)\n",GetName(), element->GetName(), GetClassVersion(), GetName() );
            } else {
               // If the element is just cached and not repeat, we need to inject an element
               // to insure the writing.
               TStreamerElement *writecopy = (TStreamerElement*)element->Clone();
               R__TObjArray_InsertAfter( fElements, writecopy, element );
               next(); // move the cursor passed the insert object.
               writecopy->SetBit(TStreamerElement::kWrite);
               writecopy->SetNewType( writecopy->GetType() );
               writecopy->SetOffset(element->GetOffset());
            }
            element->SetBit(TStreamerElement::kCache);
            element->SetNewType( element->GetType() );
            element->SetOffset(infoalloc ? infoalloc->GetOffset(element->GetName()) : 0);
         } else if (rules.HasRuleWithTarget( element->GetName(), kTRUE ) ) {
            // The data member exist in the onfile StreamerInfo and there is a rule
            // that has the same member 'only' has a target ... so this means we are
            // asked to ignore the input data ...
            if (element->GetType() == kCounter) {
               // If the element is a counter, we will need its value to read
               // other data member, so let's do so (by not disabling it) even
               // if the value will be over-written by a rule.
            } else {
               element->SetOffset(kMissing);
            }
         }
      } else if (rules && rules.HasRuleWithTarget( element->GetName(), kTRUE ) ) {
         // The data member exist in the onfile StreamerInfo and there is a rule
         // that has the same member 'only' has a target ... so this means we are
         // asked to ignore the input data ...
         if (element->GetType() == kCounter) {
            // If the element is a counter, we will need its value to read
            // other data member, so let's do so (by not disabling it) even
            // if the value will be over-written by a rule.
         } else {
            element->SetOffset(kMissing);
         }
      }

      if (element->GetNewType() == -2) {
         Warning("BuildOld", "Cannot convert %s::%s from type: %s to type: %s, skip element", GetName(), element->GetName(), element->GetTypeName(), newClass ? newClass->GetName() : (dm ? dm->GetFullTypeName() : "unknown") );
      }
   }

   // If we get here, this means that there no data member after the last base class
   // (or no base class at all).
   if (shouldHaveInfoLoc && fNVirtualInfoLoc==0) {
      fNVirtualInfoLoc = 1;
      fVirtualInfoLoc = new ULong_t[1]; // To allow for a single delete statement.
      fVirtualInfoLoc[0] = offset;
      offset += sizeof(TStreamerInfo*);
   }

   // change order , move "bazes" to the end. Workaround old bug
   if ((fOldVersion <= 2) && nBaze) {
      SetBit(kRecovered);
      TObjArray& arr = *fElements;
      TObjArray tai(nBaze);
      int narr = arr.GetLast() + 1;
      int iel;
      int jel = 0;
      int kel = 0;
      for (iel = 0; iel < narr; ++iel) {
         element = (TStreamerElement*) arr[iel];
         if (element->IsBase() && (element->IsA() != TStreamerBase::Class())) {
            tai[kel++] = element;
         } else {
            arr[jel++] = element;
         }
      }
      for (kel = 0; jel < narr;) {
         arr[jel++] = tai[kel++];
      }
   }

   // Now add artificial TStreamerElement (i.e. rules that creates new members or set transient members).
   if (!wasCompiled) InsertArtificialElements(rules);

   if (!wasCompiled && allocClass) {

      TStreamerElement *el = new TStreamerArtificial("@@alloc","", 0, TStreamerInfo::kCacheNew, allocClass->GetName());
      R__TObjArray_InsertAt( fElements, el, 0 );

      el = new TStreamerArtificial("@@dealloc","", 0, TStreamerInfo::kCacheDelete, allocClass->GetName());
      fElements->Add( el );
   }

   Compile();
}

////////////////////////////////////////////////////////////////////////////////
/// If opt contains 'built', reset this StreamerInfo as if Build or BuildOld
/// was never called on it (useful to force their re-running).

void TStreamerInfo::Clear(Option_t *option)
{
   TString opt = option;
   opt.ToLower();

   if (opt.Contains("build")) {
      R__LOCKGUARD(gInterpreterMutex);

      delete [] fComp;     fComp    = 0;
      delete [] fCompFull; fCompFull= 0;
      delete [] fCompOpt;  fCompOpt = 0;
      fNdata = 0;
      fNfulldata = 0;
      fNslots= 0;
      fSize = 0;
      ResetIsCompiled();
      ResetBit(kBuildOldUsed);

      if (fReadObjectWise) fReadObjectWise->fActions.clear();
      if (fReadMemberWise) fReadMemberWise->fActions.clear();
      if (fReadMemberWiseVecPtr) fReadMemberWiseVecPtr->fActions.clear();
      if (fReadText) fReadText->fActions.clear();
      if (fWriteObjectWise) fWriteObjectWise->fActions.clear();
      if (fWriteMemberWise) fWriteMemberWise->fActions.clear();
      if (fWriteMemberWiseVecPtr) fWriteMemberWiseVecPtr->fActions.clear();
      if (fWriteText) fWriteText->fActions.clear();
   }
}

namespace {
   // TMemberInfo
   // Local helper class to be able to compare data member represented by
   // 2 distinct TStreamerInfos
   class TMemberInfo {
   public:
      TClass  *fParent;
      TString fName;
      TString fClassName;
      TString fComment;
      Int_t   fDataType;

      TMemberInfo(TClass *parent) : fParent(parent) {};

      void SetDataType(Int_t datatype) {
         fDataType = datatype;
      }

      void SetName(const char *name) {
         fName = name;
      }
      void SetClassName(const char *name) {
         fClassName = TClassEdit::ResolveTypedef(TClassEdit::ShortType( name, TClassEdit::kDropStlDefault | TClassEdit::kDropStd ).c_str(),kTRUE);
      }
      void SetComment(const char *title) {
         const char *left = strstr(title,"[");
         if (left) {
            const char *right = strstr(left,"]");
            if (right) {
               ++left;
               fComment.Append(left,right-left);
            }
         }
      }
      void Clear() {
         fName.Clear();
         fClassName.Clear();
         fComment.Clear();
      }
      /* Hide this not yet used implementation to suppress warnings message
       from icc 11
      Bool_t operator==(const TMemberInfo &other) {
         return fName==other.fName
            && fClassName == other.fClassName
            && fComment == other.fComment;
      }
       */
      Bool_t operator!=(const TMemberInfo &other) {
         if (fName != other.fName) return kTRUE;
         if (fDataType < TStreamerInfo::kObject) {
            // For simple type, let compare the data type
            if (fDataType != other.fDataType) {
               if ( (fDataType == 4 && other.fDataType == 16)
                    || (fDataType == 16 && other.fDataType == 4) ) {
                  // long and 'long long' have the same  file format
               } else if ( (fDataType == 14 && other.fDataType == 17)
                           || (fDataType == 17 && other.fDataType == 14) ) {
                  // unsigned long and 'unsigned long long' have the same  file format
               } else if ( (fDataType == 3 && other.fDataType == 6)
                          ||(fDataType == 6 && other.fDataType == 3) ){
                  // Int_t and kCounter.  As the switch from Int_t (3) to
                  // kCounter (6) might be triggered by a derived class using
                  // the field as an array size, the class itself has no
                  // control on what the field type really use.
               } else {
                  return kTRUE;
               }
            }
         } else if (fClassName != other.fClassName) {
            if ( (fClassName == "long" && (other.fClassName == "long long" || other.fClassName == "Long64_t"))
                  || ( (fClassName == "long long" || fClassName == "Long64_t") && other.fClassName == "long") ) {
               // This is okay both have the same on file format.
            } else if ( (fClassName == "unsigned long" && (other.fClassName == "unsigned long long" || other.fClassName == "ULong64_t"))
                       || ( (fClassName == "unsigned long long" || fClassName == "ULong64_t") && other.fClassName == "unsigned long") ) {
               // This is okay both have the same on file format.
            } else if (TClassEdit::IsSTLCont(fClassName)) {
               TString name = TClassEdit::ShortType( fClassName, TClassEdit::kDropStlDefault );
               TString othername = TClassEdit::ShortType( other.fClassName, TClassEdit::kDropStlDefault );
               if (name != othername) {
                  TClass *cl = TClass::GetClass(name);
                  TClass *otherCl = TClass::GetClass(othername);
                  if (!CollectionMatch(cl,otherCl)) {
                     TClass *oldFixedClass = FixCollectionV5(fParent,cl,otherCl);
                     if (!oldFixedClass || !CollectionMatch(oldFixedClass,otherCl)) {
                        return kTRUE;
                     }
                  }
               }
            } else {
               return kTRUE;
            }
         }
         return fComment != other.fComment;
      }
   };
}

////////////////////////////////////////////////////////////////////////////////
/// Emulated a call ShowMembers() on the obj of this class type, passing insp and parent.

void TStreamerInfo::CallShowMembers(const void* obj, TMemberInspector &insp, Bool_t isTransient) const
{
   TIter next(fElements);
   TStreamerElement* element = (TStreamerElement*) next();

   TString elementName;

   for (; element; element = (TStreamerElement*) next()) {

      // Skip elements which have not been allocated memory.
      if (element->GetOffset() == kMissing) {
         continue;
      }

      char* eaddr = ((char*)obj) + element->GetOffset();

      if (element->IsBase()) {
         // Nothing to do this round.
      } else if (element->IsaPointer()) {
         elementName.Form("*%s",element->GetFullName());
         insp.Inspect(fClass, insp.GetParent(), elementName.Data(), eaddr, isTransient);
      } else {
         insp.Inspect(fClass, insp.GetParent(), element->GetFullName(), eaddr, isTransient);
         Int_t etype = element->GetType();
         switch(etype) {
            case kObject:
            case kAny:
            case kTObject:
            case kTString:
            case kTNamed:
            case kSTL:
            {
               TClass *ecl = element->GetClassPointer();
               if (ecl && (fClass!=ecl /* This happens 'artificially for stl container see the use of "This" */)) {
                  insp.InspectMember(ecl, eaddr, TString(element->GetName()) + ".", isTransient);
               }
               break;
            }
         } // switch(etype)
      } // if IsaPointer()
   } // Loop over elements

   // And now do the base classes
   next.Reset();
   element = (TStreamerElement*) next();
   for (; element; element = (TStreamerElement*) next()) {
      if (element->IsBase()) {
         // Skip elements which have not been allocated memory.
         if (element->GetOffset() == kMissing) {
            continue;
         }

         char* eaddr = ((char*)obj) + element->GetOffset();

         TClass *ecl = element->GetClassPointer();
         if (ecl) {
            ecl->CallShowMembers(eaddr, insp, isTransient);
         }
      } // If is a abse
   } // Loop over elements
}

////////////////////////////////////////////////////////////////////////////////
/// Make a clone of an object using the Streamer facility.
/// If newname is specified, this will be the name of the new object.

TObject *TStreamerInfo::Clone(const char *newname) const
{
   TStreamerInfo *newinfo = (TStreamerInfo*)TNamed::Clone(newname);
   if (newname && newname[0] && fName != newname) {
      TObjArray *newelems = newinfo->GetElements();
      Int_t ndata = newelems->GetEntries();
      for(Int_t i = 0; i < ndata; ++i) {
         TObject *element = newelems->UncheckedAt(i);
         if (element->IsA() == TStreamerLoop::Class()) {
            TStreamerLoop *eloop = (TStreamerLoop*)element;
            if (fName == eloop->GetCountClass()) {
               eloop->SetCountClass(newname);
               eloop->Init();
            }
         } else if (element->IsA() == TStreamerBasicPointer::Class()) {
            TStreamerBasicPointer *eptr = (TStreamerBasicPointer*)element;
            if (fName == eptr->GetCountClass()) {
               eptr->SetCountClass(newname);
               eptr->Init();
            }
         }
      }
   }
   ++fgCount;
   newinfo->fNumber = fgCount;
   return newinfo;
}

////////////////////////////////////////////////////////////////////////////////
/// Return True if the current StreamerInfo in cl or info is equivalent to this TStreamerInfo.
///
/// In this context 'Equivalent' means the same number of persistent data member which the same actual C++ type and
/// the same name.
/// If 'warn' is true, Warning message are printed to explicit the differences.
/// If 'complete' is false, stop at the first error, otherwise continue until all members have been checked.

Bool_t TStreamerInfo::CompareContent(TClass *cl, TVirtualStreamerInfo *info, Bool_t warn, Bool_t complete, TFile *file)
{
   Bool_t result = kTRUE;
   R__ASSERT( (cl==0 || info==0) && (cl!=0 || info!=0) /* must compare to only one thing! */);

   TString name;
   TString type;
   TStreamerElement *el;
   TStreamerElement *infoel = 0;

   TIter next(GetElements());
   TIter infonext((TList*)0);
   TIter basenext((TList*)0);
   TIter membernext((TList*)0);
   if (info) {
      infonext = info->GetElements();
   }
   if (cl) {
      TList *tlb = cl->GetListOfBases();
      if (tlb) {   // Loop over bases
         basenext = tlb;
      }
      tlb = cl->GetListOfDataMembers();
      if (tlb) {
         membernext = tlb;
      }
   }

   // First let's compare base classes
   Bool_t done = kFALSE;
   TString localClass;
   TString otherClass;
   while(!done) {
      localClass.Clear();
      otherClass.Clear();
      el = (TStreamerElement*)next();
      if (el && el->IsBase()) {
         localClass = el->GetName();
      } else {
         el = 0;
      }
      if (cl) {
         TBaseClass *tbc = (TBaseClass*)basenext();
         if (tbc) {
            otherClass = tbc->GetName();
         } else if (el==0) {
            done = kTRUE;
            break;
         }
      } else {
         infoel = (TStreamerElement*)infonext();
         if (infoel && infoel->IsBase()) {
            otherClass = infoel->GetName();
         } else if (el==0) {
            done = kTRUE;
            break;
         }
      }
      if (TClassEdit::IsSTLCont(localClass)) {
         localClass = TClassEdit::ShortType( localClass, TClassEdit::kDropStlDefault );
         otherClass = TClassEdit::ShortType( otherClass, TClassEdit::kDropStlDefault );
      }
      // Need to normalize the name
      if (localClass != otherClass) {
         if (warn) {
            if (el==0) {
               Warning("CompareContent",
                       "The in-memory layout version %d for class '%s' has a base class (%s) that the on-file layout version %d does not have.",
                       GetClassVersion(), GetName(), otherClass.Data(), GetClassVersion());
            } else if (otherClass.Length()==0) {
               Warning("CompareContent",
                       "The on-file layout version %d for class '%s'  has a base class (%s) that the in-memory layout version %d does not have",
                       GetClassVersion(), GetName(), localClass.Data(), GetClassVersion());
            } else {
               Warning("CompareContent",
                       "One base class of the on-file layout version %d and of the in memory layout version %d for '%s' is different: '%s' vs '%s'",
                       GetClassVersion(), GetClassVersion(), GetName(), localClass.Data(), otherClass.Data());
            }
         }
         if (!complete) return kFALSE;
         result = result && kFALSE;
      }
      if (cl) {
         TStreamerBase *localBase = dynamic_cast<TStreamerBase*>(el);
         if (!localBase) continue;
         // We already have localBaseClass == otherBaseClass
         TClass *otherBaseClass = localBase->GetClassPointer();
         if (!otherBaseClass) continue;
         if (otherBaseClass->IsVersioned() && localBase->GetBaseVersion() != otherBaseClass->GetClassVersion()) {
            TString msg;
            msg.Form("   The StreamerInfo of class %s read from %s%s\n"
                     "   has the same version (=%d) as the active class but a different checksum.\n"
                     "   You should update the version to ClassDef(%s,%d).\n"
                     "   The objects on this file might not be readable because:\n"
                     "   The in-memory layout version %d for class '%s' has a base class (%s) with version %d but the on-file layout version %d recorded the version number %d for this base class (%s).",
                     GetName(), file ? "file " : "", file ? file->GetName() : "", fClassVersion, GetName(), fClassVersion + 1,
                     GetClassVersion(), GetName(), otherClass.Data(), otherBaseClass->GetClassVersion(),
                     GetClassVersion(), localBase->GetBaseVersion(), localClass.Data());
            TStreamerBase *otherBase = (TStreamerBase*)cl->GetStreamerInfo()->GetElements()->FindObject(otherClass);
            otherBase->SetErrorMessage(msg);

         } else if (!otherBaseClass->IsVersioned() && localBase->GetBaseCheckSum() != otherBaseClass->GetCheckSum()) {
            TVirtualStreamerInfo *localBaseInfo = otherBaseClass->FindStreamerInfo(localBase->GetBaseCheckSum());
            if (!localBaseInfo) {
               // We are likely in the situation where the base class comes after the derived
               // class in the TFile's list of StreamerInfo, so it has not yet been loaded,
               // let's see if it is there.
               const TList *list = file->GetStreamerInfoCache();
               localBaseInfo = list ? (TStreamerInfo*)list->FindObject(localBase->GetName()) : 0;
            }
            if (!localBaseInfo) {
               TString msg;
               msg.Form("   The StreamerInfo of the base class %s (of class %s) read from %s%s\n"
                        "   refers to a checksum (%x) that can not be found neither in memory nor in the file.\n",
                        otherBaseClass->GetName(), localClass.Data(),
                        file ? "file " : "", file ? file->GetName() : "",
                        localBase->GetBaseCheckSum()
                        );
               TStreamerBase *otherBase = (TStreamerBase*)cl->GetStreamerInfo()->GetElements()->FindObject(otherClass);
               otherBase->SetErrorMessage(msg);
               continue;
            }
            if (localBaseInfo->CompareContent(otherBaseClass,0,kFALSE,kFALSE,file) ) {
               // They are equivalent, no problem.
               continue;
            }
            TString msg;
            msg.Form("   The StreamerInfo of class %s read from %s%s\n"
                     "   has the same version (=%d) as the active class but a different checksum.\n"
                     "   You should update the version to ClassDef(%s,%d).\n"
                     "   The objects on this file might not be readable because:\n"
                     "   The in-memory layout version %d for class '%s' has a base class (%s) with checksum %x but the on-file layout version %d recorded the checksum value %x for this base class (%s).",
                     GetName(), file ? "file " : "", file ? file->GetName() : "", fClassVersion, GetName(), fClassVersion + 1,
                     GetClassVersion(), GetName(), otherClass.Data(), otherBaseClass->GetCheckSum(),
                     GetClassVersion(), localBase->GetBaseCheckSum(), localClass.Data());
            TStreamerBase *otherBase = (TStreamerBase*)cl->GetStreamerInfo()->GetElements()->FindObject(otherClass);
            otherBase->SetErrorMessage(msg);
         }
      } else {
         TStreamerBase *localBase = dynamic_cast<TStreamerBase*>(el);
         TStreamerBase *otherBase = dynamic_cast<TStreamerBase*>(infoel);
         if (!localBase || !otherBase) continue;

         // We already have localBaseClass == otherBaseClass
         TClass *otherBaseClass = localBase->GetClassPointer();
         if (otherBaseClass->IsVersioned() && localBase->GetBaseVersion() != otherBase->GetBaseVersion()) {
            TString msg;
            msg.Form("   The StreamerInfo of class %s read from %s%s\n"
                     "   has the same version (=%d) as the active class but a different checksum.\n"
                     "   You should update the version to ClassDef(%s,%d).\n"
                     "   The objects on this file might not be readable because:\n"
                     "   The in-memory layout version %d for class '%s' has a base class (%s) with version %d but the on-file layout version %d recorded the version number %d for this base class (%s).",
                     GetName(), file ? "file " : "", file ? file->GetName() : "", fClassVersion, GetName(), fClassVersion + 1,
                     GetClassVersion(), GetName(), otherClass.Data(), otherBase->GetBaseVersion(),
                     GetClassVersion(), localBase->GetBaseVersion(), localClass.Data());
            otherBase->SetErrorMessage(msg);

         } else if (!otherBaseClass->IsVersioned() && localBase->GetBaseCheckSum() != otherBase->GetBaseCheckSum())
         {
            TVirtualStreamerInfo *localBaseInfo = otherBaseClass->FindStreamerInfo(localBase->GetBaseCheckSum());
            TVirtualStreamerInfo *otherBaseInfo = otherBaseClass->FindStreamerInfo(otherBase->GetBaseCheckSum());
            if (localBaseInfo == otherBaseInfo ||
                localBaseInfo->CompareContent(0,otherBaseInfo,kFALSE,kFALSE,file) ) {
               // They are equivalent, no problem.
               continue;
            }
            TString msg;
            msg.Form("   The StreamerInfo of class %s read from %s%s\n"
                     "   has the same version (=%d) as the active class but a different checksum.\n"
                     "   You should update the version to ClassDef(%s,%d).\n"
                     "   The objects on this file might not be readable because:\n"
                     "   The in-memory layout version %d for class '%s' has a base class (%s) with checksum %x but the on-file layout version %d recorded the checksum value %x for this base class (%s).",
                     GetName(), file ? "file " : "", file ? file->GetName() : "", fClassVersion, GetName(), fClassVersion + 1,
                     GetClassVersion(), GetName(), otherClass.Data(), otherBase->GetBaseCheckSum(),
                     GetClassVersion(), localBase->GetBaseCheckSum(), localClass.Data());
            otherBase->SetErrorMessage(msg);
         }
      }
   }
   if (!result && !complete) {
      return result;
   }
   // Next the datamembers
   done = kFALSE;
   next.Reset();
   infonext.Reset();

   TMemberInfo local(GetClass());
   TMemberInfo other(cl ? cl : info->GetClass());
   UInt_t idx = 0;
   while(!done) {
      local.Clear();
      other.Clear();
      el = (TStreamerElement*)next();
      while (el && (el->IsBase() || el->IsA() == TStreamerArtificial::Class())) {
         el = (TStreamerElement*)next();
         ++idx;
      }
      if (el) {
         local.SetName( el->GetName() );
         local.SetClassName( el->GetTypeName() );
         local.SetComment( el->GetTitle() );
         local.SetDataType( el->GetType() );
      }
      if (cl) {
         TDataMember *tdm = (TDataMember*)membernext();
         while(tdm && ( (!tdm->IsPersistent()) || (tdm->Property()&kIsStatic) || (el && local.fName != tdm->GetName()) )) {
            tdm = (TDataMember*)membernext();
         }
         if (tdm) {
            other.SetName( tdm->GetName() );
            other.SetClassName( tdm->GetTrueTypeName() );
            other.SetComment( tdm->GetTitle() );
            if (tdm->GetDataType()) {
               // Need to update the type for arrays.
               if (tdm->IsaPointer()) {
                  if (tdm->GetDataType()->GetType() == TVirtualStreamerInfo::kChar && !tdm->GetArrayDim() && tdm->GetArrayIndex()[0]==0) {
                     other.SetDataType( TVirtualStreamerInfo::kCharStar );
                  } else {
                     other.SetDataType( tdm->GetDataType()->GetType() + TVirtualStreamerInfo::kOffsetP);
                  }
               } else {
                  if (tdm->GetArrayDim()) {
                     other.SetDataType( tdm->GetDataType()->GetType() + TVirtualStreamerInfo::kOffsetL);
                  } else {
                     other.SetDataType( tdm->GetDataType()->GetType() );
                  }
               }
            }
         } else if (el==0) {
            done = kTRUE;
            break;
         }
      } else {
         infoel = (TStreamerElement*)infonext();
         while (infoel && (infoel->IsBase() || infoel->IsA() == TStreamerArtificial::Class())) {
            infoel = (TStreamerElement*)infonext();
         }
         if (infoel) {
            other.SetName( infoel->GetName() );
            other.SetClassName( infoel->GetTypeName() );
            other.SetComment( infoel->GetTitle() );
            other.SetDataType( infoel->GetType() );
         } else if (el==0) {
            done = kTRUE;
            break;
         }
      }
      if (local!=other) {
         if (warn) {
            if (!el) {
               Warning("CompareContent","The following data member of\nthe in-memory layout version %d of class '%s' is missing from \nthe on-file layout version %d:\n"
                       "   %s %s; //%s"
                       ,GetClassVersion(), GetName(), GetClassVersion()
                       ,other.fClassName.Data(),other.fName.Data(),other.fComment.Data());

            } else if (other.fName.Length()==0) {
               Warning("CompareContent","The following data member of\nthe in-memory layout version %d of class '%s' is missing from \nthe on-file layout version %d:\n"
                       "   %s %s; //%s"
                       ,GetClassVersion(), GetName(), GetClassVersion()
                       ,local.fClassName.Data(),local.fName.Data(),local.fComment.Data());
            } else {
               Warning("CompareContent","The following data member of\nthe on-file layout version %d of class '%s' differs from \nthe in-memory layout version %d:\n"
                       "   %s %s; //%s\n"
                       "vs\n"
                       "   %s %s; //%s"
                       ,GetClassVersion(), GetName(), GetClassVersion()
                       ,local.fClassName.Data(),local.fName.Data(),local.fComment.Data()
                       ,other.fClassName.Data(),other.fName.Data(),other.fComment.Data());
            }
         }
         result = result && kFALSE;
         if (!complete) return result;
      }
      ++idx;
   }
   return result;
}


////////////////////////////////////////////////////////////////////////////////
/// Compute total size of all persistent elements of the class

void TStreamerInfo::ComputeSize()
{
   if (this == fClass->GetCurrentStreamerInfo()) {
      if (fClass->GetState() >= TClass::kInterpreted || fClass->fIsSyntheticPair) {
         fSize = fClass->GetClassSize();
         return;
      }
   }

   TStreamerElement *element = (TStreamerElement*)fElements->Last();
   //faster and more precise to use last element offset +size
   //on 64 bit machines, offset may be forced to be a multiple of 8 bytes
   fSize = element ? element->GetOffset() + element->GetSize() : 0;
   if (fNVirtualInfoLoc > 0 && (fVirtualInfoLoc[0]+sizeof(TStreamerInfo*)) >= (ULong_t)fSize) {
      fSize = fVirtualInfoLoc[0] + sizeof(TStreamerInfo*);
   }

   // On some platform and in some case of layout non-basic data types needs
   // to be aligned.  So let's be on the safe side and align on the size of
   // the pointers.  (Question: is that the right thing on x32 ABI ?)
   constexpr size_t kSizeOfPtr = sizeof(void*);
   if ((fSize % kSizeOfPtr) != 0) {
      fSize = fSize - (fSize % kSizeOfPtr) + kSizeOfPtr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively mark streamer infos for writing to a file.
///
/// Will force this TStreamerInfo to the file and also
/// all the dependencies.
/// If argument force > 0 the loop on class dependencies is forced.
/// This function is called when streaming a class that contains
/// a null pointer. In this case, the TStreamerInfo for the class
/// with the null pointer must be written to the file and also all
/// the TStreamerInfo of all the classes referenced by the class.
/// We must be given a file to write to.

void TStreamerInfo::ForceWriteInfo(TFile* file, Bool_t force)
{
   if (!file || fNumber < 0) {
      return;
   }
   // Get the given file's list of streamer infos marked for writing.
   TArrayC* cindex = file->GetClassIndex();
   //the test below testing fArray[fNumber]>1 is to avoid a recursivity
   //problem in some cases like:
   //        class aProblemChild: public TNamed {
   //        aProblemChild *canBeNull;
   //        };
   if ( // -- Done if already marked, and we are not forcing, or forcing is blocked.
      (cindex->fArray[fNumber] && !force) || // Streamer info is already marked, and not forcing, or
      (cindex->fArray[fNumber] > 1) // == 2 means ignore forcing to prevent infinite recursion.
   ) {
      return;
   }
   // We do not want to write streamer info to the file
   // for std::string.
   static TClassRef string_classref("string");
   if (fClass == string_classref) { // We are std::string.
      return;
   }
   // We do not want to write streamer info to the file
   // for STL containers.
   if (fClass==0) {
      // Build or BuildCheck has not been called yet.
      // Let's use another means of checking.
      if (fElements && fElements->GetEntries()==1 && strcmp("This",fElements->UncheckedAt(0)->GetName())==0) {
         // We are an STL collection.
         return;
      }
   } else if (fClass->GetCollectionProxy()) { // We are an STL collection.
      return;
   }
   // Mark ourselves for output, and block
   // forcing to prevent infinite recursion.
   cindex->fArray[fNumber] = 2;
   // Signal the file that the marked streamer info list has changed.
   cindex->fArray[0] = 1;
   // Recursively mark the streamer infos for
   // all of our elements.
   TIter next(fElements);
   TStreamerElement* element = (TStreamerElement*) next();
   for (; element; element = (TStreamerElement*) next()) {
      if (element->IsTransient()) continue;
      TClass* cl = element->GetClassPointer();
      if (cl) {
         TVirtualStreamerInfo* si = 0;
         if (cl->Property() & kIsAbstract) {
            // If the class of the element is abstract, register the
            // TStreamerInfo only if it has already been built.
            // Otherwise call cl->GetStreamerInfo() would generate an
            // incorrect StreamerInfo.
            si = cl->GetCurrentStreamerInfo();
         } else {
            si = cl->GetStreamerInfo();
         }
         if (si) {
            si->ForceWriteInfo(file, force);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Assuming that obj points to (the part of) an object that is of the
/// type described by this streamerInfo, return the actual type of the
/// object (i.e. the type described by this streamerInfo is a base class
/// of the actual type of the object.
/// This routine should only be called if the class described by this
/// StreamerInfo is 'emulated'.

TClass *TStreamerInfo::GetActualClass(const void *obj) const
{
   R__ASSERT(!fClass->IsLoaded());

   if (fNVirtualInfoLoc != 0) {
      TStreamerInfo *allocator = *(TStreamerInfo**)( (const char*)obj + fVirtualInfoLoc[0] );
      if (allocator) return allocator->GetClass();
   }
   return (TClass*)fClass;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the checksum passed as argument is one of the checksum
/// value produced by the older checksum calculation algorithm.

Bool_t TStreamerInfo::MatchLegacyCheckSum(UInt_t checksum) const
{
   for(UInt_t i = 1; i < TClass::kLatestCheckSum; ++i) {
      if ( checksum == GetCheckSum( (TClass::ECheckSum) i) ) return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Recalculate the checksum of this TStreamerInfo based on its code.
///
/// The class ckecksum is used by the automatic schema evolution algorithm
/// to uniquely identify a class version.
/// The check sum is built from the names/types of base classes and
/// data members.
/// The valid range of code is determined by ECheckSum.
///   - kNoEnum:  data members of type enum are not counted in the checksum
///   - kNoRange: return the checksum of data members and base classes, not including the ranges and array size found in comments.
///   - kWithTypeDef: use the sugared type name in the calculation.
///
/// This is needed for backward compatibility.
/// ### WARNING
/// This function must be kept in sync with TClass::GetCheckSum.
/// They are both used to handle backward compatibility and should both return the same values.
/// TStreamerInfo uses the information in TStreamerElement while TClass uses the information
/// from TClass::GetListOfBases and TClass::GetListOfDataMembers.
/// Original algorithm from Victor Perevovchikov (perev@bnl.gov).

UInt_t TStreamerInfo::GetCheckSum(TClass::ECheckSum code) const
{
   // kCurrentCheckSum (0) should be kept for backward compatibility, to be
   // able to use the inequality checks, we need to set the code to the largest
   // value.
   if (code == TClass::kCurrentCheckSum) code = TClass::kLatestCheckSum;

   UInt_t id = 0;

   int il;
   TString name = GetName();
   TString type;
   il = name.Length();
   for (int i=0; i<il; i++) id = id*3+name[i];

   TIter next(GetElements());
   TStreamerElement *el;
   // Here we skip he base classes in case this is a pair or STL collection,
   // otherwise, on some STL implementations, it can happen that pair has
   // base classes which are an internal implementation detail.
   if (!fClass->GetCollectionProxy() && !TClassEdit::IsStdPair(fClass->GetName())) {
      while ( (el=(TStreamerElement*)next())) { // loop over bases
         if (el->IsBase()) {
            name = el->GetName();
            il = name.Length();
            for (int i=0; i<il; i++) id = id*3+name[i];
            if (code > TClass::kNoBaseCheckSum && el->IsA() == TStreamerBase::Class()) {
               TStreamerBase *base = (TStreamerBase*)el;
               id = id*3 + base->GetBaseCheckSum();
            }
         }
      } /* End of Base Loop */
   }

   next.Reset();
   while ( (el=(TStreamerElement*)next()) ) {
      if (el->IsBase()) continue;

      // humm can we tell if a TStreamerElement is an enum?
      // Maybe something like:
      Bool_t isenum = kFALSE;
      if ( el->GetType()==3 && gROOT->GetType(el->GetTypeName())==0) {
         // If the type is not an enum but a typedef to int then
         // el->GetTypeName() should be return 'int'
         isenum = kTRUE;
      }
      if ( (code > TClass::kNoEnum) && isenum) id = id*3 + 1;

      name = el->GetName();  il = name.Length();

      int i;
      for (i=0; i<il; i++) id = id*3+name[i];

      if (code == TClass::kReflex || code == TClass::kReflexNoComment) {
         // With TClass::kReflexV5 we do not want the Long64 in the name
         // nor any typedef.
         type = TClassEdit::ResolveTypedef(el->GetTypeName(),kTRUE);

      } else if (code <= TClass::kWithTypeDef) {
         // humm ... In the streamerInfo we only have the desugared/normalized
         // names, so we are unable to calculate the name with typedefs ...
         // except for the case of the ROOT typedef (Int_t, etc.) which are
         // kept by TClassEdit::ResolveTypedef(typeName) but not by TCling's
         // normalization ...
         //
         type = el->GetTypeName();
      } else {
         type = TClassEdit::GetLong64_Name(TClassEdit::ResolveTypedef(el->GetTypeName(),kTRUE));
      }
      if (TClassEdit::IsSTLCont(type)) {
         type = TClassEdit::ShortType( type, TClassEdit::kDropStlDefault | TClassEdit::kLong64 );
      }
      if (code == TClass::kReflex || code == TClass::kReflexNoComment) {
         type.ReplaceAll("ULong64_t","unsigned long long");
         type.ReplaceAll("Long64_t","long long");
         type.ReplaceAll("signed char","char");
         type.ReplaceAll("<signed char","<char");
         type.ReplaceAll(",signed char",",char");
         if (type=="signed char") type = "char";
      }

      il = type.Length();
      for (i=0; i<il; i++) id = id*3+type[i];

      int dim = el->GetArrayDim();
      if (dim) {
         for (i=0;i<dim;i++) id = id*3+el->GetMaxIndex(i);
      }


      if (code > TClass::kNoRange) {
         const char *left;
         if (code > TClass::kNoRangeCheck)
            left = TVirtualStreamerInfo::GetElementCounterStart(el->GetTitle());
         else
            left = strstr(el->GetTitle(),"[");
         if (left) {
            const char *right = strstr(left,"]");
            if (right) {
               ++left;
               while (left != right) {
                  id = id*3 + *left;
                  ++left;
               }
            }
         }
      }
   }
   return id;
}

////////////////////////////////////////////////////////////////////////////////

static void R__WriteConstructorBody(FILE *file, TIter &next)
{
   TStreamerElement *element = 0;
   next.Reset();
   while ((element = (TStreamerElement*)next())) {
      if (element->GetType() == TVirtualStreamerInfo::kObjectp || element->GetType() == TVirtualStreamerInfo::kObjectP ||
          element->GetType() == TVirtualStreamerInfo::kAnyp || element->GetType() == TVirtualStreamerInfo::kAnyP ||
          element->GetType() == TVirtualStreamerInfo::kCharStar || element->GetType() == TVirtualStreamerInfo::kSTLp ||
          element->GetType() == TVirtualStreamerInfo::kStreamLoop) {
         if(element->GetArrayLength() <= 1) {
            fprintf(file,"   %s = 0;\n",element->GetName());
         } else {
            fprintf(file,"   memset(%s,0,%d);\n",element->GetName(),element->GetSize());
         }
      }
      if (TVirtualStreamerInfo::kOffsetP <= element->GetType() && element->GetType() < TVirtualStreamerInfo::kObject ) {
         fprintf(file,"   %s = 0;\n",element->GetName());
      }
   }
}

static constexpr int str_length(const char* str)
{
    return *str ? 1 + str_length(str + 1) : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the element is auto_ptr or unique_ptr

static bool R__IsUniquePtr(TStreamerElement *element) {

   constexpr auto auto_ptr_len = str_length("auto_ptr<");
   constexpr auto unique_ptr_len = str_length("unique_ptr<");

   const char *name = element->GetTypeNameBasic();

   return ((strncmp(name, "auto_ptr<", auto_ptr_len) == 0)
           || (strncmp(name, "unique_ptr<", unique_ptr_len) == 0));
}

////////////////////////////////////////////////////////////////////////////////
/// Write down the body of the 'move' constructor.

static void R__WriteMoveConstructorBody(FILE *file, const TString &protoname, TIter &next)
{
   TStreamerElement *element = 0;
   next.Reset();
   Bool_t atstart = kTRUE;
   while ((element = (TStreamerElement*)next())) {
      if (element->IsBase()) {
         if (atstart) { fprintf(file,"   : "); atstart = kFALSE; }
         else fprintf(file,"   , ");
         fprintf(file, "%s(const_cast<%s &>( rhs ))\n", element->GetName(),protoname.Data());
      } else {
         if (element->GetArrayLength() <= 1) {
            if (atstart) { fprintf(file,"   : "); atstart = kFALSE; }
            else fprintf(file,"   , ");
            if (R__IsUniquePtr(element)) {
               fprintf(file, "%s(const_cast<%s &>( rhs ).%s.release() )\n",element->GetName(),protoname.Data(),element->GetName());
            } else {
               fprintf(file, "%s(const_cast<%s &>( rhs ).%s)\n",element->GetName(),protoname.Data(),element->GetName());
            }
         }
      }
   }
   fprintf(file,"{\n");
   fprintf(file,"   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).\n");
   fprintf(file,"   // Use at your own risk!\n");
   fprintf(file,"   (void)rhs; // avoid warning about unused parameter\n");
   next.Reset();
   Bool_t defMod = kFALSE;
   while ((element = (TStreamerElement*)next())) {
      if (element->GetType() == TVirtualStreamerInfo::kObjectp || element->GetType() == TVirtualStreamerInfo::kObjectP||
          element->GetType() == TVirtualStreamerInfo::kAnyp || element->GetType() == TVirtualStreamerInfo::kAnyP
          || element->GetType() == TVirtualStreamerInfo::kAnyPnoVT)
      {
         if (!defMod) { fprintf(file,"   %s &modrhs = const_cast<%s &>( rhs );\n",protoname.Data(),protoname.Data()); defMod = kTRUE; };
         const char *ename = element->GetName();
         const char *colon2 = strstr(ename,"::");
         if (colon2) ename = colon2+2;
         if(element->GetArrayLength() <= 1) {
            fprintf(file,"   modrhs.%s = 0;\n",ename);
         } else {
            fprintf(file,"   memset(modrhs.%s,0,%d);\n",ename,element->GetSize());
         }
      } else {
         const char *ename = element->GetName();
         if (element->GetType() == kCharStar) {
            if (!defMod) {
               fprintf(file,"   %s &modrhs = const_cast<%s &>( rhs );\n",protoname.Data(),protoname.Data()); defMod = kTRUE;
            };
            fprintf(file,"   modrhs.%s = 0;\n",ename);
         } else if (TVirtualStreamerInfo::kOffsetP <= element->GetType() && element->GetType() < TVirtualStreamerInfo::kObject ) {
            if (!defMod) {
               fprintf(file,"   %s &modrhs = const_cast<%s &>( rhs );\n",protoname.Data(),protoname.Data()); defMod = kTRUE;
            };
            fprintf(file,"   modrhs.%s = 0;\n",ename);
         } else if (element->GetArrayLength() > 1) {
            // FIXME: Need to add support for variable length array.
            if (element->GetArrayDim() == 1) {
               fprintf(file,"   for (Int_t i=0;i<%d;i++) %s[i] = rhs.%s[i];\n",element->GetArrayLength(),ename,ename);
            } else if (element->GetArrayDim() >= 2) {
               fprintf(file,"   for (Int_t i=0;i<%d;i++) (&(%s",element->GetArrayLength(),ename);
               for (Int_t d = 0; d < element->GetArrayDim(); ++d) {
                  fprintf(file,"[0]");
               }
               fprintf(file,"))[i] = (&(rhs.%s",ename);
               for (Int_t d = 0; d < element->GetArrayDim(); ++d) {
                  fprintf(file,"[0]");
               }
               fprintf(file,"))[i];\n");
            }
         } else if (element->GetType() == TVirtualStreamerInfo::kSTLp) {
            if (!defMod) { fprintf(file,"   %s &modrhs = const_cast<%s &>( rhs );\n",protoname.Data(),protoname.Data()); defMod = kTRUE; };
            fprintf(file,"   modrhs.%s = 0;\n",ename);
         } else if (element->GetType() == TVirtualStreamerInfo::kSTL) {
            if (!defMod) {
               fprintf(file,"   %s &modrhs = const_cast<%s &>( rhs );\n",protoname.Data(),protoname.Data()); defMod = kTRUE;
            }
            TClass *cle = element->GetClassPointer();
            TVirtualCollectionProxy *proxy = cle ? element->GetClassPointer()->GetCollectionProxy() : 0;
            std::string method_name = "clear";
            if (!element->TestBit(TStreamerElement::kDoNotDelete) && proxy && (((TStreamerSTL*)element)->GetSTLtype() == ROOT::kSTLbitset)) {
                method_name = "reset";
            }
            if (element->IsBase()) {
               fprintf(file,"   modrhs.%s();\n", method_name.c_str());
            } else {
               fprintf(file,"   modrhs.%s.%s();\n",ename, method_name.c_str());
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

static void R__WriteDestructorBody(FILE *file, TIter &next)
{
   TStreamerElement *element = 0;
   next.Reset();
   while ((element = (TStreamerElement*)next())) {
      if (element->GetType() == TVirtualStreamerInfo::kObjectp || element->GetType() == TVirtualStreamerInfo::kObjectP||
          element->GetType() == TVirtualStreamerInfo::kAnyp || element->GetType() == TVirtualStreamerInfo::kAnyP
          || element->GetType() == TVirtualStreamerInfo::kAnyPnoVT)
      {
         const char *ename = element->GetName();
         const char *colon2 = strstr(ename,"::");
         if (colon2) ename = colon2+2;
         if (element->TestBit(TStreamerElement::kDoNotDelete)) {
            if(element->GetArrayLength() <= 1) {
               fprintf(file,"   %s = 0;\n",ename);
            } else {
               fprintf(file,"   memset(%s,0,%d);\n",ename,element->GetSize());
            }
         } else {
            if(element->GetArrayLength() <= 1) {
               fprintf(file,"   delete %s;   %s = 0;\n",ename,ename);
            } else {
               fprintf(file,"   for (Int_t i=0;i<%d;i++) delete %s[i];   memset(%s,0,%d);\n",element->GetArrayLength(),ename,ename,element->GetSize());
            }
         }
      }
      if (element->GetType() == TVirtualStreamerInfo::kCharStar) {
         const char *ename = element->GetName();
         if (element->TestBit(TStreamerElement::kDoNotDelete)) {
            fprintf(file,"   %s = 0;\n",ename);
         } else {
            fprintf(file,"   delete [] %s;   %s = 0;\n",ename,ename);
         }
      }
      if (TVirtualStreamerInfo::kOffsetP <= element->GetType() && element->GetType() < TVirtualStreamerInfo::kObject ) {
         const char *ename = element->GetName();
         if (element->TestBit(TStreamerElement::kDoNotDelete)) {
            fprintf(file,"   %s = 0;\n",ename);
         } else if (element->HasCounter()) {
            fprintf(file,"   delete %s;   %s = 0;\n",ename,ename);
         } else {
            fprintf(file,"   delete [] %s;   %s = 0;\n",ename,ename);
         }
      }
      if (element->GetType() == TVirtualStreamerInfo::kSTL || element->GetType() == TVirtualStreamerInfo::kSTLp) {
         const char *ename = element->GetName();
         const char *prefix = "";
         if ( element->GetType() == TVirtualStreamerInfo::kSTLp ) {
            prefix = "*";
         } else if ( element->IsBase() ) {
            ename = "this";
         }
         TClass *cle = element->GetClassPointer();
         TVirtualCollectionProxy *proxy = cle ? element->GetClassPointer()->GetCollectionProxy() : 0;
         if (!element->TestBit(TStreamerElement::kDoNotDelete) && proxy) {
            Int_t stltype = ((TStreamerSTL*)element)->GetSTLtype();

            if (proxy->HasPointers()) {
               fprintf(file,"   std::for_each( (%s %s).rbegin(), (%s %s).rend(), DeleteObjectFunctor() );\n",prefix,ename,prefix,ename);
               //fprintf(file,"      %s::iterator iter;\n");
               //fprintf(file,"      %s::iterator begin = (%s %s).begin();\n");
               //fprintf(file,"      %s::iterator end (%s %s).end();\n");
               //fprintf(file,"      for( iter = begin; iter != end; ++iter) { delete *iter; }\n");
            } else {
               if (stltype == ROOT::kSTLmap || stltype == ROOT::kSTLmultimap) {
                  TString enamebasic = TMakeProject::UpdateAssociativeToVector(element->GetTypeNameBasic());
                  std::vector<std::string> inside;
                  int nestedLoc;
                  TClassEdit::GetSplit(enamebasic, inside, nestedLoc, TClassEdit::kLong64);
                  if ((!inside[1].empty() && inside[1][inside[1].size()-1]=='*')
                     || (!inside[2].empty() && inside[2][inside[2].size()-1]=='*')) {
                     fprintf(file,"   std::for_each( (%s %s).rbegin(), (%s %s).rend(), DeleteObjectFunctor() );\n",prefix,ename,prefix,ename);
                  }
               }
               }
         }
         if ( prefix[0] ) {
            fprintf(file,"   delete %s;   %s = 0;\n",ename,ename);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write the Declaration of class.

void TStreamerInfo::GenerateDeclaration(FILE *fp, FILE *sfp, const TList *subClasses, Bool_t top)
{
   if (fClassVersion == -3) {
      return;
   }

   Bool_t needGenericTemplate = fElements==0 || fElements->GetEntries() == 0;
   Bool_t isTemplate = kFALSE;
   const char *clname = GetName();
   TString template_protoname;
   if (strchr(clname, ':')) {
      // We might have a namespace in front of the classname.
      Int_t len = strlen(clname);
      const char *name = clname;
      UInt_t nest = 0;
      UInt_t pr_pos = 0;
      for (Int_t cur = 0; cur < len; ++cur) {
         switch (clname[cur]) {
            case '<':
               ++nest;
               pr_pos = cur;
               isTemplate = kTRUE;
               break;
            case '>':
               if (nest == 0) { cur = len; continue; } // the name is not well formed, give up.
               --nest;
               break;
            case ':': {
               if (nest == 0 && clname[cur+1] == ':') {
                  // We have a scope
                  isTemplate = kFALSE;
                  name = clname + cur + 2;
               }
               break;
            }
         }
      }
      if (isTemplate) {
         template_protoname.Append(clname,pr_pos);
      }
      clname = name;
   } else {
      const char *where = strstr(clname, "<");
      isTemplate = where != 0;
      if (isTemplate) {
         template_protoname.Append(clname,where-clname);
      }
   }

   if (needGenericTemplate && isTemplate) {
      TString templateName(TMakeProject::GetHeaderName("template "+template_protoname,0));
      fprintf(fp, "#ifndef %s_h\n", templateName.Data());
      fprintf(fp, "#define %s_h\n", templateName.Data());
   }

   TString protoname;
   UInt_t numberOfNamespaces = TMakeProject::GenerateClassPrefix(fp, GetName(), top, protoname, 0, kFALSE, needGenericTemplate);

   // Generate class statement with base classes.
   TStreamerElement *element;
   TIter next(fElements);
   Int_t nbase = 0;
   while ((element = (TStreamerElement*)next())) {
      if (!element->IsBase()) continue;
      nbase++;
      const char *ename = element->GetName();
      if (nbase == 1) fprintf(fp," : public %s",ename);
      else            fprintf(fp," , public %s",ename);
   }
   fprintf(fp," {\n");

   // Generate forward declaration nested classes.
   if (subClasses && subClasses->GetEntries()) {
      Bool_t needheader = true;

      TIter subnext(subClasses);
      TStreamerInfo *subinfo;
      Int_t len = strlen(GetName());
      while ((subinfo = (TStreamerInfo*)subnext())) {
         if (strncmp(GetName(),subinfo->GetName(),len)==0 && (subinfo->GetName()[len]==':') ) {
            if (subinfo->GetName()[len+1]==':' && strstr(subinfo->GetName()+len+2,":")==0) {
               if (needheader) {
                  fprintf(fp,"\npublic:\n");
                  fprintf(fp,"// Nested classes forward declaration.\n");
                  needheader = false;
               }
               TString sub_protoname;
               UInt_t sub_numberOfClasses = 0;
               UInt_t sub_numberOfNamespaces;
               if (subinfo->GetClassVersion() == -3) {
                  sub_numberOfNamespaces = TMakeProject::GenerateClassPrefix(fp, subinfo->GetName() + len+2, kFALSE, sub_protoname, &sub_numberOfClasses, 3);
               } else {
                  sub_numberOfNamespaces = TMakeProject::GenerateClassPrefix(fp, subinfo->GetName() + len+2, kFALSE, sub_protoname, &sub_numberOfClasses, kFALSE);
                  fprintf(fp, ";\n");
               }

               for (UInt_t i = 0;i < sub_numberOfClasses;++i) {
                  fprintf(fp, "}; // end of class.\n");
               }
               if (sub_numberOfNamespaces > 0) {
                  Error("GenerateDeclaration","Nested classes %s thought to be inside a namespace inside the class %s",subinfo->GetName(),GetName());
               }
            }
         }
      }
   }

   fprintf(fp,"\npublic:\n");
   fprintf(fp,"// Nested classes declaration.\n");

   // Generate nested classes.
   if (subClasses && subClasses->GetEntries()) {
      TIter subnext(subClasses,kIterBackward);
      TStreamerInfo *subinfo;
      Int_t len = strlen(GetName());
      while ((subinfo = (TStreamerInfo*)subnext())) {
         if (strncmp(GetName(),subinfo->GetName(),len)==0 && (subinfo->GetName()[len]==':')) {
            if (subinfo->GetName()[len+1]==':' && strstr(subinfo->GetName()+len+2,":")==0) {
               subinfo->GenerateDeclaration(fp, sfp, subClasses, kFALSE);
            }
         }
      }
   }

   fprintf(fp,"\npublic:\n");
   fprintf(fp,"// Data Members.\n");

   {
      // Generate data members.
      TString name(128);
      Int_t ltype = 12;
      Int_t ldata = 10;
      Int_t lt,ld,is;
      TString line;
      line.Resize(kMaxLen);
      next.Reset();
      while ((element = (TStreamerElement*)next())) {

         if (element->IsBase()) continue;
         const char *ename = element->GetName();

         name = ename;
         for (Int_t i=0;i < element->GetArrayDim();i++) {
            name += TString::Format("[%d]",element->GetMaxIndex(i));
         }
         name += ";";
         ld = name.Length();

         TString enamebasic = element->GetTypeNameBasic();
         if (element->IsA() == TStreamerSTL::Class()) {
            // If we have a map, multimap, set or multiset,
            // and the key is a class, we need to replace the
            // container by a vector since we don't have the
            // comparator function.
            Int_t stltype = ((TStreamerSTL*)element)->GetSTLtype();
            switch (stltype) {
               case ROOT::kSTLmap:
               case ROOT::kSTLmultimap:
               case ROOT::kSTLset:
               case ROOT::kSTLmultiset:
               case ROOT::kSTLunorderedset:
               case ROOT::kSTLunorderedmultiset:
               {
                  enamebasic = TMakeProject::UpdateAssociativeToVector(enamebasic);
               }
               default:
                  // nothing to do.
                  break;
            }
         } else if (strncmp(enamebasic.Data(), "auto_ptr<", strlen("auto_ptr<")) == 0) {
            enamebasic = TMakeProject::UpdateAssociativeToVector(enamebasic);
         }

         lt = enamebasic.Length();

         line = "   ";
         line += enamebasic;
         if (lt>=ltype) ltype = lt+1;

         for (is = 3+lt; is < (3+ltype); ++is) line += ' ';

         line += name;
         if (element->IsaPointer() && !strchr(line,'*')) line[2+ltype] = '*';

         if (ld>=ldata) ldata = ld+1;
         for (is = 3+ltype+ld; is < (3+ltype+ldata); ++is) line += ' ';

         line += "   //";
         line += element->GetTitle();
         fprintf(fp,"%s\n",line.Data());
      }
   }
   if (needGenericTemplate && isTemplate) {
      // Generate default functions, ClassDef and trailer.
      fprintf(fp,"\n   %s() {\n",protoname.Data());
      R__WriteConstructorBody(fp,next);
      fprintf(fp,"   }\n");
      fprintf(fp,"   %s(%s && ) = default;\n",protoname.Data(),protoname.Data());
      fprintf(fp,"   %s(const %s & rhs )\n",protoname.Data(),protoname.Data());
      R__WriteMoveConstructorBody(fp,protoname,next);
      fprintf(fp,"   }\n");
      fprintf(fp,"   virtual ~%s() {\n",protoname.Data());
      R__WriteDestructorBody(fp,next);
      fprintf(fp,"   }\n\n");

   } else {
      // Generate default functions, ClassDef and trailer.
      fprintf(fp,"\n   %s();\n",protoname.Data());
      fprintf(fp,"   %s(%s && ) = default;\n",protoname.Data(),protoname.Data());
      fprintf(fp,"   %s(const %s & );\n",protoname.Data(),protoname.Data());
      fprintf(fp,"   virtual ~%s();\n\n",protoname.Data());

      // Add the implementations to the source.cxx file.
      TString guard( TMakeProject::GetHeaderName( GetName(), 0, kTRUE ) );
      fprintf(sfp,"#ifndef %s_cxx\n",guard.Data());
      fprintf(sfp,"#define %s_cxx\n",guard.Data());
      fprintf(sfp,"%s::%s() {\n",GetName(),protoname.Data());
      R__WriteConstructorBody(sfp,next);
      fprintf(sfp,"}\n");

      fprintf(sfp,"%s::%s(const %s & rhs)\n",GetName(),protoname.Data(),protoname.Data());
      R__WriteMoveConstructorBody(sfp,protoname,next);
      fprintf(sfp,"}\n");

      fprintf(sfp,"%s::~%s() {\n",GetName(),protoname.Data());
      R__WriteDestructorBody(sfp,next);
      fprintf(sfp,"}\n");
      fprintf(sfp,"#endif // %s_cxx\n\n",guard.Data());
   }

   TClass *cl = gROOT->GetClass(GetName());
   if (fClassVersion > 1 || (cl && cl->IsTObject()) ) {
      // add 1 to class version in case we didn't manage reproduce the class layout to 100%.
      if (fClassVersion == 0) {
         // If the class was declared 'transient', keep it that way.
         fprintf(fp,"   ClassDef(%s,%d); // Generated by MakeProject.\n",protoname.Data(),0);
      } else {
         fprintf(fp,"   ClassDef(%s,%d); // Generated by MakeProject.\n",protoname.Data(),fClassVersion + 1);
      }
   }
   fprintf(fp,"};\n");

   for(UInt_t i=0;i<numberOfNamespaces;++i) {
      fprintf(fp,"} // namespace\n");
   }

   if (needGenericTemplate && isTemplate) {
      fprintf(fp,"#endif // generic template declaration\n");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add to the header file, the \#include need for this class.

UInt_t TStreamerInfo::GenerateIncludes(FILE *fp, char *inclist, const TList *extrainfos)
{
   if (inclist[0]==0) {
      // Always have this include for ClassDef.
      TMakeProject::AddInclude( fp, "Rtypes.h", kFALSE, inclist);
   }
   UInt_t ninc = 0;

   const char *clname = GetName();
   if (strchr(clname,'<')) {
      // This is a template, we need to check the template parameter.
      ninc += TMakeProject::GenerateIncludeForTemplate(fp, clname, inclist, kFALSE, extrainfos);
   }

   TString name(1024);
   Int_t ltype = 10;
   Int_t ldata = 10;
   Int_t lt;
   Int_t ld;
   TIter next(fElements);
   TStreamerElement *element;
   Bool_t incRiostream = kFALSE;
   while ((element = (TStreamerElement*)next())) {
      //if (element->IsA() == TStreamerBase::Class()) continue;
      const char *ename = element->GetName();
      const char *colon2 = strstr(ename,"::");
      if (colon2) ename = colon2+2;
      name = ename;
      for (Int_t i=0;i < element->GetArrayDim();i++) {
         name += TString::Format("[%d]",element->GetMaxIndex(i));
      }
      ld = name.Length();
      lt = strlen(element->GetTypeName());
      if (ltype < lt) ltype = lt;
      if (ldata < ld) ldata = ld;

      //must include Riostream.h in case of an STL container
      if (!incRiostream && element->InheritsFrom(TStreamerSTL::Class())) {
         incRiostream = kTRUE;
         TMakeProject::AddInclude( fp, "Riostream.h", kFALSE, inclist);
      }

      //get include file name if any
      const char *include = element->GetInclude();
      if (!include[0]) continue;

      Bool_t greater = (include[0]=='<');
      include++;

      if (strncmp(include,"include/",8)==0) {
         include += 8;
      }
      if (strncmp(include,"include\\",9)==0) {
         include += 9;
      }
      if (TClassEdit::IsStdPair(element->GetTypeName())) {
         TMakeProject::AddInclude( fp, "utility", kTRUE, inclist);
      } else if (strncmp(element->GetTypeName(),"auto_ptr<",strlen("auto_ptr<"))==0) {
         TMakeProject::AddInclude( fp, "memory", kTRUE, inclist);
      } else {
         TString incName( include, strlen(include)-1 );
         incName = TMakeProject::GetHeaderName(incName,extrainfos);
         TMakeProject::AddInclude( fp, incName.Data(), greater, inclist);
      }

      if (strchr(element->GetTypeName(),'<')) {
         // This is a template, we need to check the template parameter.
         ninc += TMakeProject::GenerateIncludeForTemplate(fp, element->GetTypeName(), inclist, kFALSE, extrainfos);
      }
   }
   return ninc;
}

////////////////////////////////////////////////////////////////////////////////
/// Generate header file for the class described by this TStreamerInfo
/// the function is called by TFile::MakeProject for each class in the file

Int_t TStreamerInfo::GenerateHeaderFile(const char *dirname, const TList *subClasses, const TList *extrainfos)
{
   // if (fClassVersion == -4) return 0;
   if ((fClass && fClass->GetCollectionType()) || TClassEdit::IsSTLCont(GetName())) return 0;
   if (TClassEdit::IsStdPair(GetName())) return 0;
   if (strncmp(GetName(),"auto_ptr<",strlen("auto_ptr<"))==0) return 0;

   TClass *cl = TClass::GetClass(GetName());
   if (cl) {
      if (cl->HasInterpreterInfo()) return 0; // skip known classes
   }
   Bool_t isTemplate = kFALSE;
   if (strchr(GetName(),':')) {
      UInt_t len = strlen(GetName());
      UInt_t nest = 0;
      UInt_t scope = 0;
      for(UInt_t i=len; i>0; --i) {
         switch(GetName()[i]) {
            case '>': ++nest; if (scope==0) { isTemplate = kTRUE; } break;
            case '<': --nest; break;
            case ':':
               if (nest==0 && GetName()[i-1]==':') {
                  // We have a scope
                  TString nsname(GetName(), i-1);
                  cl = gROOT->GetClass(nsname);
                  if (cl && (cl->Size()!=0 || (cl->Size()==0 && !cl->HasInterpreterInfo() /*empty 'base' class on file*/))) {
                     // This class is actually nested.
                     return 0;
                  } else if (cl == 0 && extrainfos != 0) {
                     TStreamerInfo *clinfo = (TStreamerInfo*)extrainfos->FindObject(nsname);
                     if (clinfo && clinfo->GetClassVersion() == -5) {
                        // This class is actually nested.
                        return 0;
                     }
                  }
                  ++scope;
               }
               break;
         }
      }
   }
   Bool_t needGenericTemplate = isTemplate && (fElements==0 || fElements->GetEntries()==0);

   if (gDebug) printf("generating code for class %s\n",GetName());

   // Open the file

   TString headername( TMakeProject::GetHeaderName( GetName(), extrainfos ) );
   TString filename;
   filename.Form("%s/%s.h",dirname,headername.Data());

   FILE *fp = fopen(filename.Data(),"w");
   if (!fp) {
      Error("MakeProject","Cannot open output file:%s\n",filename.Data());
      return 0;
   }

   filename.Form("%s/%sProjectHeaders.h",dirname,gSystem->BaseName(dirname));
   FILE *allfp = fopen(filename.Data(),"a");
   if (!allfp) {
      Error("MakeProject","Cannot open output file:%s\n",filename.Data());
      fclose(fp);
      return 0;
   }
   fprintf(allfp,"#include \"%s.h\"\n", headername.Data());
   fclose(allfp);

   char *inclist = new char[50000];
   inclist[0] = 0;

   // Generate class header.
   TDatime td;
   fprintf(fp,"//////////////////////////////////////////////////////////\n");
   fprintf(fp,"//   This class has been generated by TFile::MakeProject\n");
   fprintf(fp,"//     (%s by ROOT version %s)\n",td.AsString(),gROOT->GetVersion());
   fprintf(fp,"//      from the StreamerInfo in file %s\n",gDirectory->GetFile()->GetName());
   fprintf(fp,"//////////////////////////////////////////////////////////\n");
   fprintf(fp,"\n");
   fprintf(fp,"\n");
   fprintf(fp,"#ifndef %s_h\n",headername.Data());
   fprintf(fp,"#define %s_h\n",headername.Data());
   TMakeProject::GenerateForwardDeclaration(fp, GetName(), inclist, kFALSE, needGenericTemplate, extrainfos);
   fprintf(fp,"\n");

   UInt_t ninc = 0;
   ninc += GenerateIncludes(fp, inclist, extrainfos);
   if (subClasses) {
      TIter subnext(subClasses);
      TStreamerInfo *subinfo;
      while ((subinfo = (TStreamerInfo*)subnext())) {
         ninc = subinfo->GenerateIncludes(fp, inclist, extrainfos);
      }
   }
   fprintf(fp,"\n");

   TString sourcename; sourcename.Form( "%s/%sProjectSource.cxx", dirname, gSystem->BaseName(dirname) );
   FILE *sfp = fopen( sourcename.Data(), "a" );
   if (sfp) {
      GenerateDeclaration(fp, sfp, subClasses);
   } else {
      Error("GenerateHeaderFile","Could not open %s for appending",sourcename.Data());
   }
   TMakeProject::GeneratePostDeclaration(fp, this, inclist);

   fprintf(fp,"#endif\n");

   delete [] inclist;
   fclose(fp);
   if (sfp) fclose(sfp);
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute data member offset.
/// Return pointer to the Streamer function if one exists

Int_t TStreamerInfo::GetDataMemberOffset(TDataMember *dm, TMemberStreamer *&streamer) const
{
   TIter nextr(fClass->GetListOfRealData());
   char dmbracket[256];
   snprintf(dmbracket,255,"%s[",dm->GetName());
   Int_t offset = kMissing;
   if (!fClass->IsLoaded()) {
      // If the 'class' is not loaded, we do not have a TClass bootstrap and thus
      // the 'RealData' might not have enough information because of the lack
      // of proper ShowMember implementation.
      if (! (dm->Property() & kIsStatic) ) {
         // Give an offset only to non-static members.
         offset = dm->GetOffset();
      }
   }
   TRealData *rdm;
   while ((rdm = (TRealData*)nextr())) {
      char *rdmc = (char*)rdm->GetName();
      //next statement required in case a class and one of its parent class
      //have data members with the same name
      if (dm->IsaPointer() && rdmc[0] == '*') rdmc++;

      if (rdm->GetDataMember() != dm) continue;
      if (strcmp(rdmc,dm->GetName()) == 0) {
         offset   = rdm->GetThisOffset();
         streamer = rdm->GetStreamer();
         break;
      }
      if (strcmp(rdm->GetName(),dm->GetName()) == 0) {
         if (rdm->IsObject()) {
            offset = rdm->GetThisOffset();
            streamer = rdm->GetStreamer();
            break;
         }
      }
      if (strstr(rdm->GetName(),dmbracket)) {
         offset   = rdm->GetThisOffset();
         streamer = rdm->GetStreamer();
         break;
      }
   }
   return offset;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the offset of the data member as indicated by this StreamerInfo.

Int_t TStreamerInfo::GetOffset(const char *elementName) const
{
   if (elementName == 0) return 0;

   Int_t offset = 0;
   TStreamerElement *elem = (TStreamerElement*)fElements->FindObject(elementName);
   if (elem) offset = elem->GetOffset();

   return offset;
}

////////////////////////////////////////////////////////////////////////////////
/// Return total size of all persistent elements of the class (with offsets).

Int_t TStreamerInfo::GetSize() const
{
   return fSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Return total size of all persistent elements of the class
/// use GetSize if you want to get the real size in memory.

Int_t TStreamerInfo::GetSizeElements() const
{
   TIter next(fElements);
   TStreamerElement *element;
   Int_t asize = 0;
   while ((element = (TStreamerElement*)next())) {
      asize += element->GetSize();
   }
   return asize;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the StreamerElement of "datamember" inside our
/// class or any of its base classes.
///
/// The offset information
/// contained in the StreamerElement is related to its immediately
/// containing class, so we return in 'offset' the offset inside
/// our class.

TStreamerElement* TStreamerInfo::GetStreamerElement(const char* datamember, Int_t& offset) const
{
   if (!fElements) {
      return 0;
   }

   // Look first at the data members and base classes
   // of our class.
   TStreamerElement* element = (TStreamerElement*) fElements->FindObject(datamember);
   if (element) {
      offset = element->GetOffset();
      return element;
   }

   // Not found, so now try the data members and base classes
   // of the base classes of our class.
   if (fClass->HasDataMemberInfo()) {
      // Our class has a dictionary loaded, use it to search the base classes.
      TStreamerElement* base_element = 0;
      TBaseClass* base = 0;
      TClass* base_cl = 0;
      Int_t base_offset = 0;
      Int_t local_offset = 0;
      TIter nextb(fClass->GetListOfBases());
      // Iterate on list of base classes.
      while ((base = (TBaseClass*) nextb())) {
         base_cl = TClass::GetClass(base->GetName());
         base_element = (TStreamerElement*) fElements->FindObject(base->GetName());
         if (!base_cl || !base_element) {
            continue;
         }
         base_offset = base_element->GetOffset();
         element = ((TStreamerInfo*)base_cl->GetStreamerInfo())->GetStreamerElement(datamember, local_offset);
         if (element) {
            offset = base_offset + local_offset;
            return element;
         }
      }
   } else {
      // Our class's dictionary is not loaded. Search through the base class streamer elements.
      TIter next(fElements);
      TStreamerElement* curelem = 0;
      while ((curelem = (TStreamerElement*) next())) {
         if (curelem->InheritsFrom(TStreamerBase::Class())) {
            TClass* baseClass = curelem->GetClassPointer();
            if (!baseClass) {
               continue;
            }
            Int_t base_offset = curelem->GetOffset();
            Int_t local_offset = 0;
            TStreamerInfo *baseInfo;
            if (baseClass->Property() & kIsAbstract) {
               baseInfo = (TStreamerInfo*)baseClass->GetStreamerInfoAbstractEmulated();
            } else {
               baseInfo = (TStreamerInfo*)baseClass->GetStreamerInfo();
            }
            if (baseInfo) element = baseInfo->GetStreamerElement(datamember, local_offset);
            if (element) {
               offset = base_offset + local_offset;
               return element;
            }
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// <b>Obsolete</b>: this routine is obsolete and should not longer be used.
///
/// TStreamerInfo  holds two types of data structures
///   - TObjArray* fElements; containing the list of all TStreamerElement
///       objects for this class version.
///   - ULong_t*  fElem;  containing the preprocessed information
///       by TStreamerInfo::Compile In case consecutive data members
///       are of the same type, the Compile function declares the consecutive
///       elements as one single element in fElems.
///
/// Example with the class TAttLine:
/// ~~~{.cpp}
/// TClass::GetClass("TAttLine")->GetStreamerInfo()->ls(); produces;
///  StreamerInfo for class: TAttLine, version=1
///   short        fLineColor      offset=  4 type= 2 line color
///   short        fLineStyle      offset=  6 type= 2 line style
///   short        fLineWidth      offset=  8 type= 2 line width
///   i= 0, fLineColor      type= 22, offset=  4, len=3, method=0
/// ~~~
///  For I/O implementations (eg. XML) , one has to know the original name
///  of the data member. This function can be used to return a pointer
///  to the original TStreamerElement object corresponding to the j-th
///  element of a compressed array in fElems.
///  Parameters description:
///    - i: the serial number in array fElem
///    - j: the element number in the array of consecutive types
///  In the above example the class TAttLine has 3 consecutive data members
///  of the same type "short". Compile makes one single array of 3 elements.
///  To access the TStreamerElement for the second element
///  of this array, one can call:
/// ~~~{.cpp}
/// auto el = GetStreamerElementReal(0,1);
/// auto membername = el->GetName();
/// ~~~
///  This function is typically called from TBuffer, TXmlBuffer.

TStreamerElement* TStreamerInfo::GetStreamerElementReal(Int_t i, Int_t j) const
{
   ::Obsolete("TStreamerInfo::GetStreamerElementReal", "v5-34-20", "v6-00-02");

   if (i < 0 || i >= fNdata) return 0;
   if (j < 0) return 0;
   if (!fElements) return 0;
   TStreamerElement *se = (TStreamerElement*)fCompOpt[i]->fElem;
   if (!se) return 0;
   Int_t nelems = fElements->GetEntriesFast();
   for (Int_t ise=0;ise < nelems;ise++) {
      if (se != (TStreamerElement*)fElements->UncheckedAt(ise)) continue;
      if (ise+j >= nelems) return 0;
      return (TStreamerElement*)fElements->UncheckedAt(ise+j);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the value from inside a collection.

template <typename T>
T TStreamerInfo::GetTypedValueAux(Int_t type, void *ladd, Int_t k, Int_t len)
{
   if (type>=kConv && type<kSTL) {
      type -= kConv;
   }
   switch (type) {
      // basic types
      case kBool:              {Bool_t *val   = (Bool_t*)ladd;   return T(*val);}
      case kChar:              {Char_t *val   = (Char_t*)ladd;   return T(*val);}
      case kShort:             {Short_t *val  = (Short_t*)ladd;  return T(*val);}
      case kInt:               {Int_t *val    = (Int_t*)ladd;    return T(*val);}
      case kLong:              {Long_t *val   = (Long_t*)ladd;   return T(*val);}
      case kLong64:            {Long64_t *val = (Long64_t*)ladd; return T(*val);}
      case kFloat:             {Float_t *val  = (Float_t*)ladd;  return T(*val);}
      case kFloat16:           {Float_t *val  = (Float_t*)ladd;  return T(*val);}
      case kDouble:            {Double_t *val = (Double_t*)ladd; return T(*val);}
      case kDouble32:          {Double_t *val = (Double_t*)ladd; return T(*val);}
      case kUChar:             {UChar_t *val  = (UChar_t*)ladd;  return T(*val);}
      case kUShort:            {UShort_t *val = (UShort_t*)ladd; return T(*val);}
      case kUInt:              {UInt_t *val   = (UInt_t*)ladd;   return T(*val);}
      case kULong:             {ULong_t *val  = (ULong_t*)ladd;  return T(*val);}
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
      case kULong64:           {Long64_t *val = (Long64_t*)ladd;  return T(*val);}
#else
      case kULong64:           {ULong64_t *val= (ULong64_t*)ladd; return T(*val);}
#endif
      case kBits:              {UInt_t *val   = (UInt_t*)ladd;   return T(*val);}

         // array of basic types  array[8]
      case kOffsetL + kBool:    {Bool_t *val   = (Bool_t*)ladd;   return T(val[k]);}
      case kOffsetL + kChar:    {Char_t *val   = (Char_t*)ladd;   return T(val[k]);}
      case kOffsetL + kShort:   {Short_t *val  = (Short_t*)ladd;  return T(val[k]);}
      case kOffsetL + kInt:     {Int_t *val    = (Int_t*)ladd;    return T(val[k]);}
      case kOffsetL + kLong:    {Long_t *val   = (Long_t*)ladd;   return T(val[k]);}
      case kOffsetL + kLong64:  {Long64_t *val = (Long64_t*)ladd; return T(val[k]);}
      case kOffsetL + kFloat:   {Float_t *val  = (Float_t*)ladd;  return T(val[k]);}
      case kOffsetL + kFloat16: {Float_t *val  = (Float_t*)ladd;  return T(val[k]);}
      case kOffsetL + kDouble:  {Double_t *val = (Double_t*)ladd; return T(val[k]);}
      case kOffsetL + kDouble32:{Double_t *val = (Double_t*)ladd; return T(val[k]);}
      case kOffsetL + kUChar:   {UChar_t *val  = (UChar_t*)ladd;  return T(val[k]);}
      case kOffsetL + kUShort:  {UShort_t *val = (UShort_t*)ladd; return T(val[k]);}
      case kOffsetL + kUInt:    {UInt_t *val   = (UInt_t*)ladd;   return T(val[k]);}
      case kOffsetL + kULong:   {ULong_t *val  = (ULong_t*)ladd;  return T(val[k]);}
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
      case kOffsetL + kULong64: {Long64_t *val = (Long64_t*)ladd;  return T(val[k]);}
#else
      case kOffsetL + kULong64:{ULong64_t *val= (ULong64_t*)ladd; return T(val[k]);}
#endif

#define READ_ARRAY(TYPE_t)                               \
         {                                               \
            Int_t sub_instance, index;                   \
            Int_t instance = k;                          \
            if (len) {                                   \
               index = instance / len;                   \
               sub_instance = instance % len;            \
            } else {                                     \
               index = instance;                         \
               sub_instance = 0;                         \
            }                                            \
            TYPE_t **val =(TYPE_t**)(ladd);              \
            return T((val[sub_instance])[index]); \
         }

         // pointer to an array of basic types  array[n]
      case kOffsetP + kBool_t:    READ_ARRAY(Bool_t)
      case kOffsetP + kChar_t:    READ_ARRAY(Char_t)
      case kOffsetP + kShort_t:   READ_ARRAY(Short_t)
      case kOffsetP + kInt_t:     READ_ARRAY(Int_t)
      case kOffsetP + kLong_t:    READ_ARRAY(Long_t)
      case kOffsetP + kLong64_t:  READ_ARRAY(Long64_t)
      case kOffsetP + kFloat16_t:
      case kOffsetP + kFloat_t:   READ_ARRAY(Float_t)
      case kOffsetP + kDouble32_t:
      case kOffsetP + kDouble_t:  READ_ARRAY(Double_t)
      case kOffsetP + kUChar_t:   READ_ARRAY(UChar_t)
      case kOffsetP + kUShort_t:  READ_ARRAY(UShort_t)
      case kOffsetP + kUInt_t:    READ_ARRAY(UInt_t)
      case kOffsetP + kULong_t:   READ_ARRAY(ULong_t)
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
      case kOffsetP + kULong64_t: READ_ARRAY(Long64_t)
#else
      case kOffsetP + kULong64_t: READ_ARRAY(ULong64_t)
#endif

          // array counter //[n]
      case kCounter:           {Int_t *val    = (Int_t*)ladd;    return T(*val);}
   }
   return 0;
}



template Double_t TStreamerInfo::GetTypedValue(char *pointer, Int_t i, Int_t j, Int_t len) const;
template Long64_t TStreamerInfo::GetTypedValue(char *pointer, Int_t i, Int_t j, Int_t len) const;
template LongDouble_t TStreamerInfo::GetTypedValue(char *pointer, Int_t i, Int_t j, Int_t len) const;

////////////////////////////////////////////////////////////////////////////////
/// Return value of element i in object at pointer.
/// The function may be called in two ways:
///   - method1 len < 0: i is assumed to be the TStreamerElement number i in StreamerInfo
///   - method2 len >= 0: i is the type, address of variable is directly pointer.

template <typename T>
T TStreamerInfo::GetTypedValue(char *pointer, Int_t i, Int_t j, Int_t len) const
{
   char *ladd;
   Int_t atype;
   if (len >= 0) {
      ladd  = pointer;
      atype = i;
   } else {
      if (i < 0) return 0;
      ladd  = pointer + fCompFull[i]->fOffset;
      atype = fCompFull[i]->fNewType;
      len = fCompFull[i]->fElem->GetArrayLength();
      if (atype == kSTL) {
         TClass *newClass = fCompFull[i]->fElem->GetNewClass();
         if (newClass == 0) {
            newClass = fCompFull[i]->fElem->GetClassPointer();
         }
         TClass *innerClass = newClass->GetCollectionProxy()->GetValueClass();
         if (innerClass) {
            return 0; // We don't know which member of the class we would want.
         } else {
            TVirtualCollectionProxy *proxy = newClass->GetCollectionProxy();
            // EDataType is a subset of TStreamerInfo::EReadWrite
            atype = (TStreamerInfo::EReadWrite)proxy->GetType();
            TVirtualCollectionProxy::TPushPop pop(proxy,ladd);
            Int_t nc = proxy->Size();
            if (j >= nc) return 0;
            char *element_ptr = (char*)proxy->At(j);
            return GetTypedValueAux<T>(atype,element_ptr,0,1);
         }
      }
   }
   return GetTypedValueAux<T>(atype,ladd,j,len);
}

////////////////////////////////////////////////////////////////////////////////

template Double_t TStreamerInfo::GetTypedValueClones<Double_t>(TClonesArray *clones, Int_t i, Int_t j, int k, Int_t eoffset) const;
template Long64_t TStreamerInfo::GetTypedValueClones(TClonesArray *clones, Int_t i, Int_t j, int k, Int_t eoffset) const;
template LongDouble_t TStreamerInfo::GetTypedValueClones(TClonesArray *clones, Int_t i, Int_t j, int k, Int_t eoffset) const;

template <typename T>
T TStreamerInfo::GetTypedValueClones(TClonesArray *clones, Int_t i, Int_t j, int k, Int_t eoffset) const
{
   //  return value of element i in object number j in a TClonesArray and eventually
   // element k in a sub-array.

   Int_t nc = clones->GetEntriesFast();
   if (j >= nc) return 0;

   char *pointer = (char*)clones->UncheckedAt(j);
   char *ladd    = pointer + eoffset + fCompFull[i]->fOffset;
   return GetTypedValueAux<T>(fCompFull[i]->fType,ladd,k,((TStreamerElement*)fCompFull[i]->fElem)->GetArrayLength());
}

template Double_t TStreamerInfo::GetTypedValueSTL(TVirtualCollectionProxy *cont, Int_t i, Int_t j, int k, Int_t eoffset) const;
template Long64_t TStreamerInfo::GetTypedValueSTL(TVirtualCollectionProxy *cont, Int_t i, Int_t j, int k, Int_t eoffset) const;
template LongDouble_t TStreamerInfo::GetTypedValueSTL(TVirtualCollectionProxy *cont, Int_t i, Int_t j, int k, Int_t eoffset) const;

////////////////////////////////////////////////////////////////////////////////
/// Return value of element i in object number j in a TClonesArray and eventually
/// element k in a sub-array.

template <typename T>
T TStreamerInfo::GetTypedValueSTL(TVirtualCollectionProxy *cont, Int_t i, Int_t j, int k, Int_t eoffset) const
{
   Int_t nc = cont->Size();
   if (j >= nc) return 0;

   char *pointer = (char*)cont->At(j);
   char *ladd    = pointer + eoffset + fCompFull[i]->fOffset;
   return GetTypedValueAux<T>(fCompFull[i]->fType,ladd,k,((TStreamerElement*)fCompFull[i]->fElem)->GetArrayLength());
}

template Double_t TStreamerInfo::GetTypedValueSTLP(TVirtualCollectionProxy *cont, Int_t i, Int_t j, int k, Int_t eoffset) const;
template Long64_t TStreamerInfo::GetTypedValueSTLP(TVirtualCollectionProxy *cont, Int_t i, Int_t j, int k, Int_t eoffset) const;
template LongDouble_t TStreamerInfo::GetTypedValueSTLP(TVirtualCollectionProxy *cont, Int_t i, Int_t j, int k, Int_t eoffset) const;

////////////////////////////////////////////////////////////////////////////////
/// Return value of element i in object number j in a TClonesArray and eventually
/// element k in a sub-array.

template <typename T>
T TStreamerInfo::GetTypedValueSTLP(TVirtualCollectionProxy *cont, Int_t i, Int_t j, int k, Int_t eoffset) const
{
   Int_t nc = cont->Size();

   if (j >= nc) return 0;

   char **ptr = (char**)cont->At(j);
   char *pointer = *ptr;

   char *ladd    = pointer + eoffset + fCompFull[i]->fOffset;
   return GetTypedValueAux<T>(fCompFull[i]->fType,ladd,k,((TStreamerElement*)fCompFull[i]->fElem)->GetArrayLength());
}

////////////////////////////////////////////////////////////////////////////////
/// Insert new members as expressed in the array of TSchemaRule(s).

void TStreamerInfo::InsertArtificialElements(std::vector<const ROOT::TSchemaRule*> &rules)
{
   if (rules.empty()) return;

   TIter next(fElements);
   UInt_t count = 0;

   for(auto rule : rules) {
      if( rule->IsRenameRule() || rule->IsAliasRule() )
         continue;
      next.Reset();
      TStreamerElement *element;
      while ((element = (TStreamerElement*) next())) {
         if ( rule->HasTarget( element->GetName() ) ) {

            // Check whether this is an 'attribute' rule.
            if ( rule->GetAttributes()[0] != 0 ) {
               TString attr( rule->GetAttributes() );
               attr.ToLower();
               if (attr.Contains("owner")) {
                  if (attr.Contains("notowner")) {
                     element->SetBit(TStreamerElement::kDoNotDelete);
                  } else {
                     element->ResetBit(TStreamerElement::kDoNotDelete);
                  }
               }

            }
            break;
         }
      }

      // NOTE: Before adding the rule we should check that the source do
      // existing in this StreamerInfo.
      const TObjArray *sources = rule->GetSource();
      TIter input(sources);
      TObject *src;
      while((src = input())) {
         if ( !GetElements()->FindObject(src->GetName()) ) {
            // Missing source.
#if 0 // Don't warn about not activating the rule.  If don't warn the user can
      // have more flexibility in specifying when the rule applies and relying
      // on both the version number *and* the presence of the source members.
      // Activating this warning would for example mean that we need to carefully
      // tweak $ROOTSYS/etc/class.rules.
            TString ruleStr;
            rule->AsString(ruleStr);
            Warning("InsertArtificialElements","For class %s in StreamerInfo %d is missing the source data member %s when trying to apply the rule:\n   %s",
                   GetName(),GetClassVersion(),src->GetName(),ruleStr.Data());
            rule = 0;
#endif
            break;
         }
      }

      if (!rule) continue;

      TStreamerArtificial *newel;
      typedef std::vector<TStreamerArtificial*> vec_t;
      vec_t toAdd;

      if (rule->GetTarget()==0) {
         TString newName;
         newName.Form("%s_rule%d",fClass->GetName(),count);
         newel = new TStreamerArtificial(newName,"",
                                         fClass->GetDataMemberOffset(newName),
                                         TStreamerInfo::kArtificial,
                                         "void");
         newel->SetBit(TStreamerElement::kWholeObject);
         newel->SetReadFunc( rule->GetReadFunctionPointer() );
         newel->SetReadRawFunc( rule->GetReadRawFunctionPointer() );
         toAdd.push_back(newel);
      } else {
         toAdd.reserve(rule->GetTarget()->GetEntries());
         TObjString * objstr = (TObjString*)(rule->GetTarget()->At(0));
         if (objstr) {
            TString newName = objstr->String();
            TString realDataName;
            if ( TDataMember* dm = fClass->GetDataMember( newName ) ) {
               TRealData::GetName(realDataName,dm);
               newel = new TStreamerArtificial(realDataName,"",
                                               fClass->GetDataMemberOffset(newName),
                                               TStreamerInfo::kArtificial,
                                               fClass->GetDataMember( newName )->GetTypeName());
               newel->SetReadFunc( rule->GetReadFunctionPointer() );
               newel->SetReadRawFunc( rule->GetReadRawFunctionPointer() );
               toAdd.push_back(newel);
            } else {
               // This would be a completely new member (so it would need to be cached)
               // TOBEDONE
            }
            for(Int_t other = 1; other < rule->GetTarget()->GetEntries(); ++other) {
               objstr = (TObjString*)(rule->GetTarget()->At(other));
               if (objstr) {
                  newName = objstr->String();
                  if ( TDataMember* dm = fClass->GetDataMember( newName ) ) {
                     TRealData::GetName(realDataName,dm);
                     newel = new TStreamerArtificial(realDataName,"",
                                                     fClass->GetDataMemberOffset(newName),
                                                     TStreamerInfo::kArtificial,
                                                     fClass->GetDataMember( newName )->GetTypeName());
                     toAdd.push_back(newel);
                  }
               }
            }
         } // For each target of the rule
      }
      // Now find we with need to add them
      TIter s_iter(rule->GetSource());
      Int_t loc = -1;
      while( TObjString *s = (TObjString*)s_iter() ) {
         for(Int_t i = fElements->GetLast(); i >= 0 && (i+1) >= loc; --i) {
            if (s->String() == fElements->UncheckedAt(i)->GetName()) {
               if (loc == -1 || (i+1)>loc) {
                  loc = i+1;
               }
            }
         }
      }
      if (loc == -1) {
         // Verify if the last one is not 'skipped'.
         for(Int_t i = fElements->GetLast(); i >= 0 && (i+1) >= loc; --i) {
            if ( ((TStreamerElement*)fElements->UncheckedAt(i))->GetNewType() != -2 ) {
               break;
            }
            loc = i;
         }
      }
      if (loc == -1) {
         for(vec_t::iterator iter = toAdd.begin(); iter != toAdd.end(); ++iter) {
            fElements->Add(*iter);
         }
      } else {
         R__TObjArray_InsertAt(fElements, toAdd, loc);
      }
   } // None of the target of the rule are on file.
}

////////////////////////////////////////////////////////////////////////////////
///  List the TStreamerElement list and also the precomputed tables
///  if option contains the string "incOrig", also prints the original
///  (non-optimized elements in the list of compiled elements.

void TStreamerInfo::ls(Option_t *option) const
{
   if (fClass && (fName != fClass->GetName())) {
      if (fClass->IsVersioned()) {
         Printf("\nStreamerInfo for conversion to %s from: %s, version=%d, checksum=0x%x",fClass->GetName(),GetName(),fClassVersion,GetCheckSum());
      } else {
         Printf("\nStreamerInfo for conversion to %s from: %s, checksum=0x%x",fClass->GetName(),GetName(),GetCheckSum());
      }
   } else {
      if (!fClass || fClass->IsVersioned()) {
         Printf("\nStreamerInfo for class: %s, version=%d, checksum=0x%x",GetName(),fClassVersion,GetCheckSum());
      } else {
         Printf("\nStreamerInfo for class: %s, checksum=0x%x",GetName(),GetCheckSum());
      }
   }

   if (fElements) {
      TIter    next(fElements);
      TObject *obj;
      while ((obj = next()))
         obj->ls(option);
   }
   if (strstr(option,"full") != 0) {
      for (Int_t i=0; i < fNfulldata; ++i) {
         TStreamerElement *element = (TStreamerElement*)fCompFull[i]->fElem;
         TString sequenceType;
         element->GetSequenceType(sequenceType);
         // by definition of the loop (i+1) <= fNdata
         if (sequenceType.Length()) {
            sequenceType.Prepend(" [");
            sequenceType += "]";
         }
         Printf("   i=%2d, %-15s type=%3d, offset=%3d, len=%d, method=%ld%s",
                i,element->GetName(),fCompFull[i]->fType,fCompFull[i]->fOffset,fCompFull[i]->fLength,fCompFull[i]->fMethod,
                sequenceType.Data());
      }

   } else {
      Bool_t wantOrig = strstr(option,"incOrig") != 0;
      Bool_t optimized = kFALSE;
      for (Int_t i=0,j=0;i < fNdata;++i,++j) {
         TStreamerElement *element = (TStreamerElement*)fCompOpt[i]->fElem;
         TString sequenceType;
         element->GetSequenceType(sequenceType);
         // by definition of the loop (i+1) <= fNdata
         optimized = TStreamerInfo::kOffsetL < fCompOpt[i]->fType && fCompOpt[i]->fType < TStreamerInfo::kOffsetP && fCompOpt[i]->fLength > fCompOpt[i]->fElem->GetArrayLength();
         if (optimized) {
            // This was optimized.
            if (sequenceType.Length() != 0) {
               sequenceType += ',';
            }
            sequenceType += "optimized";
         }
         if (sequenceType.Length()) {
            sequenceType.Prepend(" [");
            sequenceType += "]";
         }
         Printf("   i=%2d, %-15s type=%3d, offset=%3d, len=%d, method=%ld%s",
                i,element->GetName(),fCompOpt[i]->fType,fCompOpt[i]->fOffset,fCompOpt[i]->fLength,fCompOpt[i]->fMethod,
                sequenceType.Data());
         if (optimized && wantOrig) {
            Bool_t done;
            do {
               element = (TStreamerElement*)fCompFull[j]->fElem;
               element->GetSequenceType(sequenceType);
               if (sequenceType.Length()) {
                  sequenceType.Prepend(" [");
                  sequenceType += "]";
               }
               Printf("      j=%2d, %-15s type=%3d, offset=%3d, len=%d, method=%ld%s",
                      j,element->GetName(),fCompFull[j]->fType,fCompFull[j]->fOffset,fCompFull[j]->fLength,fCompFull[j]->fMethod,
                      sequenceType.Data());
               ++j;
               done = j >= fNfulldata || ( (i+1 < fNdata) && fCompOpt[i+1]->fElem == fCompFull[j+1]->fElem );
            } while (!done);

         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// An emulated object is created at address obj, if obj is null we
/// allocate memory for the object.

void* TStreamerInfo::New(void *obj)
{
   //???FIX ME: What about varying length array elements?

   char* p = (char*) obj;

   TIter next(fElements);

   if (!p) {
      // Allocate and initialize the memory block.
      p = new char[fSize];
      memset(p, 0, fSize);
   }

   next.Reset();
   TStreamerElement* element = (TStreamerElement*) next();

   for (; element; element = (TStreamerElement*) next()) {

      // Skip elements which have not been allocated memory.
      if (element->GetOffset() == kMissing) {
         continue;
      }

      // Skip elements for which we do not have any class
      // information.  FIXME: Document how this could happen.
      TClass* cle = element->GetClassPointer();
      if (!cle) {
         continue;
      }

      char* eaddr = p + element->GetOffset();
      Int_t etype = element->GetType();

      //cle->GetStreamerInfo(); //necessary in case "->" is not specified

      switch (etype) {

         case kAnyP:
         case kObjectP:
         case kSTLp:
         {
            // Initialize array of pointers with null pointers.
            char** r = (char**) eaddr;
            Int_t len = element->GetArrayLength();
            for (Int_t i = 0; i < len; ++i) {
               r[i] = 0;
            }
         }
         break;

         case kObjectp:
         case kAnyp:
         {
            // If the option "->" is given in the data member comment field
            // it is assumed that the object exists before reading data in,
            // so we create an object.
            if (cle != TClonesArray::Class()) {
               void** r = (void**) eaddr;
               *r = cle->New();
            } else {
               // In the case of a TClonesArray, the class name of
               // the contained objects must be specified in the
               // data member comment in this format:
               //    TClonesArray* myVar; //->(className)
               const char* title = element->GetTitle();
               const char* bracket1 = strrchr(title, '(');
               const char* bracket2 = strrchr(title, ')');
               if (bracket1 && bracket2 && (bracket2 != (bracket1 + 1))) {
                  Int_t len = bracket2 - (bracket1 + 1);
                  char* clonesClass = new char[len+1];
                  clonesClass[0] = '\0';
                  strncat(clonesClass, bracket1 + 1, len);
                  void** r = (void**) eaddr;
                  *r = (void*) new TClonesArray(clonesClass);
                  delete[] clonesClass;
               } else {
                  //Warning("New", "No class name found for TClonesArray initializer in data member comment (expected \"//->(className)\"");
                  void** r = (void**) eaddr;
                  *r = (void*) new TClonesArray();
               }
            }
         }
         break;

         case kBase:
         {
            if (cle->Property() & kIsAbstract) {
               TVirtualStreamerInfo *einfo = cle->GetStreamerInfoAbstractEmulated();
               if (einfo) einfo->New(eaddr);
            } else {
               cle->New(eaddr);
            }
            break;
         }
         case kObject:
         case kAny:
         case kTObject:
         case kTString:
         case kTNamed:
         {
            cle->New(eaddr);
         }
         break;

         case kSTL:
         {
            if (strcmp(element->GetName(),"This")==0 &&
                !cle->GetCollectionProxy()) {
               // missing information, avoid infinite loop
               // by doing nothing ....
            } else {
               cle->New(eaddr);
            }
         }
         break;

         case kObject + kOffsetL:
         case kAny + kOffsetL:
         case kTObject + kOffsetL:
         case kTString + kOffsetL:
         case kTNamed + kOffsetL:
         case kSTL + kOffsetL:
         {
            Int_t size = cle->Size();
            char* r = eaddr;
            Int_t len = element->GetArrayLength();
            for (Int_t i = 0; i < len; ++i, r += size) {
               cle->New(r);
            }
         }
         break;

      } // switch etype
   } // for TIter next(fElements)

   for(int nbase = 0; nbase < fNVirtualInfoLoc; ++nbase) {
      *(TStreamerInfo**)(p + fVirtualInfoLoc[nbase]) = this;
   }
   return p;
}

////////////////////////////////////////////////////////////////////////////////
/// An array of emulated objects is created at address ary, if ary is null,
/// we allocate memory for the array.

void* TStreamerInfo::NewArray(Long_t nElements, void *ary)
{
   if (fClass == 0) {
      Error("NewArray", "TClass pointer is null!");
      return 0;
   }

   Int_t size = fClass->Size();

   char* p = (char*) ary;

   if (!p) {
      Long_t len = nElements * size + sizeof(Long_t)*2;
      p = new char[len];
      memset(p, 0, len);
   }

   // Store the array cookie
   Long_t* r = (Long_t*) p;
   r[0] = size;
   r[1] = nElements;
   char* dataBegin = (char*) &r[2];

   // Do a placement new for each element.
   p = dataBegin;
   for (Long_t cnt = 0; cnt < nElements; ++cnt) {
      New(p);
      p += size;
   } // for nElements

   return dataBegin;
}


#define DeleteBasicPointer(addr,element,name)                           \
   {                                                                    \
      name **f = (name**)(addr);                                        \
      int n = element->GetArrayLength() ? element->GetArrayLength() : 1;\
      for(int j=0;j<n;j++) {                                            \
         delete [] f[j];                                                \
         f[j] = 0;                                                      \
      }                                                                 \
   }

////////////////////////////////////////////////////////////////////////////////
/// Internal part of the destructor.
/// Destruct each of the datamembers in the same order
/// as the implicit destructor would.

void TStreamerInfo::DestructorImpl(void* obj, Bool_t dtorOnly)
{
   R__ASSERT(obj != 0);

   char *p = (char*)obj;

   Int_t nelements = fElements->GetEntriesFast();
   //for (; ele; ele = (TStreamerElement*) next())
   for (Int_t elenum = nelements - 1; elenum >= 0; --elenum) {
      TStreamerElement* ele = (TStreamerElement*) fElements->UncheckedAt(elenum);
      if (ele->GetOffset() == kMissing) continue;
      char* eaddr = p + ele->GetOffset();


      Int_t etype = ele->GetType();

      switch(etype) {
         case TStreamerInfo::kOffsetP + TStreamerInfo::kBool:   DeleteBasicPointer(eaddr,ele,Bool_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kChar:   DeleteBasicPointer(eaddr,ele,Char_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kShort:  DeleteBasicPointer(eaddr,ele,Short_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kInt:    DeleteBasicPointer(eaddr,ele,Int_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kLong:   DeleteBasicPointer(eaddr,ele,Long_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kLong64: DeleteBasicPointer(eaddr,ele,Long64_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat16:
         case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat:  DeleteBasicPointer(eaddr,ele,Float_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble32:
         case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble: DeleteBasicPointer(eaddr,ele,Double_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUChar:  DeleteBasicPointer(eaddr,ele,UChar_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUShort: DeleteBasicPointer(eaddr,ele,UShort_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUInt:   DeleteBasicPointer(eaddr,ele,UInt_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kULong:  DeleteBasicPointer(eaddr,ele,ULong_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kULong64:DeleteBasicPointer(eaddr,ele,ULong64_t);  continue;
         case TStreamerInfo::kCharStar:                         DeleteBasicPointer(eaddr,ele,Char_t);  continue;
      }



      TClass* cle = ele->GetClassPointer();
      if (!cle) continue;


      if (etype == kObjectp || etype == kAnyp) {
         // Destroy an array of pre-allocated objects.
         Int_t len = ele->GetArrayLength();
         if (!len) {
            len = 1;
         }
         void** r = (void**) eaddr;
         for (Int_t j = len - 1; j >= 0; --j) {
            if (r[j]) {
               cle->Destructor(r[j]);
               r[j] = 0;
            }
         }
      }

      if ((etype == kObjectP || etype == kAnyP || etype == kSTLp) && !ele->TestBit(TStreamerElement::kDoNotDelete)) {
         // Destroy an array of pointers to not-pre-allocated objects.
         Int_t len = ele->GetArrayLength();
         if (!len) {
            len = 1;
         }
         void** r = (void**) eaddr;
         for (Int_t j = len - 1; j >= 0; --j) {
            if (r[j]) {
               cle->Destructor(r[j]);
               r[j] = 0;
            }
         }
      }

      if (etype == kBase) {
         if (cle->Property() & kIsAbstract) {
            TVirtualStreamerInfo *einfo = cle->GetStreamerInfoAbstractEmulated();
            if (einfo) einfo->Destructor(eaddr, kTRUE);
         } else {
            cle->Destructor(eaddr, kTRUE);
         }
     }

      if (etype == kObject || etype == kAny ||
          etype == kTObject || etype == kTString || etype == kTNamed) {
         // A data member is destroyed, but not deleted.
         cle->Destructor(eaddr, kTRUE);
     }

      if (etype == kSTL) {
         // A data member is destroyed, but not deleted.
         TVirtualCollectionProxy *pr = cle->GetCollectionProxy();
         if (!pr) {
            if (strcmp(ele->GetName(),"This")==0) {
               // missing information, avoid infinite loop
               // by doing nothing ....
            } else {
               cle->Destructor(eaddr, kTRUE);
            }
         } else {
            if (ele->TestBit(TStreamerElement::kDoNotDelete)) {
               TVirtualCollectionProxy::TPushPop env(cle->GetCollectionProxy(), eaddr); // used for both this 'clear' and the 'clear' inside destructor.
               cle->GetCollectionProxy()->Clear(); // empty the collection without deleting the pointer
               pr->Destructor(eaddr, kTRUE);
            } else {
               pr->Destructor(eaddr, kTRUE);
            }
         }
      }

      if (etype == kObject  + kOffsetL || etype == kAny     + kOffsetL ||
          etype == kTObject + kOffsetL || etype == kTString + kOffsetL ||
          etype == kTNamed  + kOffsetL || etype == kSTL     + kOffsetL) {
         // For a data member which is an array of objects, we
         // destroy the objects, but do not delete them.
         Int_t len = ele->GetArrayLength();
         Int_t size = cle->Size();
         char* r = eaddr + (size * (len - 1));
         for (Int_t j = len - 1; j >= 0; --j, r -= size) {
            cle->Destructor(r, kTRUE);
         }
      }
   } // iter over elements

   if (!dtorOnly) {
      delete[] p;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Emulated destructor for this class.
///
/// An emulated object is destroyed at address p.
/// Destruct each of the datamembers in the same order
/// as the implicit destructor would.

void TStreamerInfo::Destructor(void* obj, Bool_t dtorOnly)
{
   // Do nothing if passed a null pointer.
   if (obj == 0) return;

   char* p = (char*) obj;

   if (!dtorOnly && fNVirtualInfoLoc) {
      // !dtorOnly is used to filter out the case where this is called for
      // a base class or embedded object of the outer most class.
      TStreamerInfo *allocator = *(TStreamerInfo**)(p + fVirtualInfoLoc[0]);
      if (allocator != this) {

         Int_t baseoffset = allocator->GetClass()->GetBaseClassOffset(GetClass());

         p -= baseoffset;
         allocator->DestructorImpl(p, kFALSE);
         return;
      }
   }
   DestructorImpl(p, dtorOnly);
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy an array of emulated objects, with optional delete.

void TStreamerInfo::DeleteArray(void* ary, Bool_t dtorOnly)
{
   // Do nothing if passed a null pointer.
   if (ary == 0) return;

   //???FIX ME: What about varying length arrays?

   Long_t* r = (Long_t*) ary;
   Long_t arrayLen = r[-1];
   Long_t size = r[-2];
   char* memBegin = (char*) &r[-2];

   char* p = ((char*) ary) + ((arrayLen - 1) * size);
   for (Long_t cnt = 0; cnt < arrayLen; ++cnt, p -= size) {
      // Destroy each element, but do not delete it.
      Destructor(p, kTRUE);
   } // for arrayItemSize

   if (!dtorOnly) {
      delete[] memBegin;
   }
}

////////////////////////////////////////////////////////////////////////////////
///  print value of element i in object at pointer
///  The function may be called in two ways:
///    -method1  len < 0
///           i is assumed to be the TStreamerElement number i in StreamerInfo
///    -method2  len >= 0
///           i is the type
///           address of variable is directly pointer.
///           len is the number of elements to be printed starting at pointer.

void TStreamerInfo::PrintValue(const char *name, char *pointer, Int_t i, Int_t len, Int_t lenmax) const
{
   char *ladd;
   Int_t atype,aleng;
   printf(" %-15s = ",name);

   TStreamerElement * aElement  = 0;
   Int_t *count  = 0;
   if (len >= 0) {
      ladd  = pointer;
      atype = i;
      aleng = len;
   } else        {
      if (i < 0) {
         if (pointer==0) {
            printf("NULL\n");
         } else {
            const static TClassRef stringClass("string");
            if (fClass == stringClass) {
               std::string *st = (std::string*)(pointer);
               printf("%s\n",st->c_str());
            } else if (fClass == TString::Class()) {
               TString *st = (TString*)(pointer);
               printf("%s\n",st->Data());
            } else {
               printf("(%s*)0x%lx\n",GetName(),(ULong_t)pointer);
            }
         }
         return;
      }
      ladd  = pointer + fCompFull[i]->fOffset;
      atype = fCompFull[i]->fNewType;
      aleng = fCompFull[i]->fLength;
      aElement  = (TStreamerElement*)fCompFull[i]->fElem;
      count = (Int_t*)(pointer+fCompFull[i]->fMethod);
   }
   if (aleng > lenmax) aleng = lenmax;

   PrintValueAux(ladd,atype,aElement,aleng,count);
   printf("\n");
}

////////////////////////////////////////////////////////////////////////////////
/// Print value of element i in a TClonesArray.

void TStreamerInfo::PrintValueClones(const char *name, TClonesArray *clones, Int_t i, Int_t eoffset, Int_t lenmax) const
{
   if (!clones) {printf(" %-15s = \n",name); return;}
   printf(" %-15s = ",name);
   Int_t nc = clones->GetEntriesFast();
   if (nc > lenmax) nc = lenmax;

   Int_t offset = eoffset + fCompFull[i]->fOffset;
   TStreamerElement *aElement  = (TStreamerElement*)fCompFull[i]->fElem;
   int aleng = fCompFull[i]->fLength;
   if (aleng > lenmax) aleng = lenmax;

   for (Int_t k=0;k < nc;k++) {
      char *pointer = (char*)clones->UncheckedAt(k);
      char *ladd = pointer+offset;
      Int_t *count = (Int_t*)(pointer+fCompFull[i]->fMethod);
      PrintValueAux(ladd,fCompFull[i]->fNewType,aElement, aleng, count);
      if (k < nc-1) printf(", ");
   }
   printf("\n");
}

////////////////////////////////////////////////////////////////////////////////
/// Print value of element i in a TClonesArray.

void TStreamerInfo::PrintValueSTL(const char *name, TVirtualCollectionProxy *cont, Int_t i, Int_t eoffset, Int_t lenmax) const
{
   if (!cont) {printf(" %-15s = \n",name); return;}
   printf(" %-15s = ",name);
   Int_t nc = cont->Size();
   if (nc > lenmax) nc = lenmax;

   Int_t offset = eoffset + fCompFull[i]->fOffset;
   TStreamerElement *aElement  = (TStreamerElement*)fCompFull[i]->fElem;
   int aleng = fCompFull[i]->fLength;
   if (aleng > lenmax) aleng = lenmax;

   for (Int_t k=0;k < nc;k++) {
      char *pointer = (char*)cont->At(k);
      char *ladd = pointer+offset;
      Int_t *count = (Int_t*)(pointer+fCompFull[i]->fMethod);
      PrintValueAux(ladd,fCompFull[i]->fNewType,aElement, aleng, count);
      if (k < nc-1) printf(", ");
   }
   printf("\n");
}

////////////////////////////////////////////////////////////////////////////////
/// Replace the TClass this streamerInfo is pointing to (belongs to)

void TStreamerInfo::SetClass(TClass *newcl)
{
   if (newcl) {
      // This is mostly (but not only) for the artificial "This" streamerElement for an stl collection.
      Update(fClass, newcl);
   }
   fClass = newcl;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStreamerInfo.

void TStreamerInfo::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      fOldVersion = R__v;
      if (R__v > 1) {
         //R__b.ReadClassBuffer(TStreamerInfo::Class(), this, R__v, R__s, R__c);
         R__b.ClassBegin(TStreamerInfo::Class(), R__v);
         R__b.ClassMember("TNamed");
         TNamed::Streamer(R__b);
         fName = TClassEdit::GetLong64_Name( fName.Data() ).c_str();
         R__b.ClassMember("fCheckSum","UInt_t");
         R__b >> fCheckSum;
         R__b.ClassMember("fClassVersion","Int_t");
         R__b >> fClassVersion;
         fOnFileClassVersion = fClassVersion;
         R__b.ClassMember("fElements","TObjArray*");
         R__b >> fElements;
         R__b.ClassEnd(TStreamerInfo::Class());
         R__b.SetBufferOffset(R__s+R__c+sizeof(UInt_t));
         ResetBit(kIsCompiled);
         ResetBit(kBuildOldUsed);
         ResetBit(kBuildRunning);

         if (R__b.GetParent() && R__b.GetVersionOwner() < 50000)
         {
            // In some older files, the type of the TStreamerElement was not
            // as we (now) expect.
            Int_t nobjects = fElements->GetEntriesFast();
            TClass *basic = TStreamerBasicType::Class();
            for (Int_t i = 0; i < nobjects; i++) {
               TStreamerElement *el = (TStreamerElement*)fElements->UncheckedAt(i);
               TStreamerElement *rel = 0;
               if ( el->IsA() == basic ) {
                  switch (el->GetType()) {
                     default: break; /* nothing */
                     case TStreamerInfo::kObject: /*61*/
                        rel = new TStreamerObject(el->GetName(),el->GetTitle(),el->GetOffset(),el->GetTypeName());
                        break;
                     case TStreamerInfo::kAny: /*62*/
                        rel = new TStreamerObjectAny(el->GetName(),el->GetTitle(),el->GetOffset(),el->GetTypeName());
                        break;
                     case TStreamerInfo::kObjectp: /* 63 */
                        rel = new TStreamerObjectPointer(el->GetName(),el->GetTitle(),el->GetOffset(),el->GetTypeName());
                        break;
                     case TStreamerInfo::kObjectP: /* 64 */
                        rel = new TStreamerObjectPointer(el->GetName(),el->GetTitle(),el->GetOffset(),el->GetTypeName());
                        break;
                     case TStreamerInfo::kTString: /* 65 */
                        rel = new TStreamerObject(el->GetName(),el->GetTitle(),el->GetOffset(),el->GetTypeName());
                        break;
                  }
                  if (rel) {
                     (*fElements)[i] = rel;
                     delete el;
                  }
               }
            }
         }
         return;
      }
      //====process old versions before automatic schema evolution
      TNamed::Streamer(R__b);
      fName = TClassEdit::GetLong64_Name( fName.Data() ).c_str();
      R__b >> fCheckSum;
      R__b >> fClassVersion;
      fOnFileClassVersion = fClassVersion;
      R__b >> fElements;
      R__b.CheckByteCount(R__s, R__c, TStreamerInfo::IsA());
   } else {
      R__c = R__b.WriteVersion(TStreamerInfo::IsA(), kTRUE);
      R__b.ClassBegin(TStreamerInfo::Class());
      R__b.ClassMember("TNamed");
      TNamed::Streamer(R__b);
      R__b.ClassMember("fCheckSum","UInt_t");
      R__b << fCheckSum;
      R__b.ClassMember("fClassVersion","Int_t");
      R__b << ((fClassVersion > 0) ? fClassVersion : -fClassVersion);

      //------------------------------------------------------------------------
      // Stream only non-artificial streamer elements
      //////////////////////////////////////////////////////////////////////////

      R__b.ClassMember("fElements","TObjArray*");
      {
         TObjArray elements(fElements->GetEntriesFast());
         TStreamerElement *el;
         Int_t nobjects = fElements->GetEntriesFast();
         for (Int_t i = 0; i < nobjects; i++) {
            el = (TStreamerElement *)fElements->UncheckedAt(i);
            if (el != 0 && (el->IsA() == TStreamerArtificial::Class() || el->TestBit(TStreamerElement::kRepeat))) {
               // skip
            } else if (el != 0 && (el->TestBit(TStreamerElement::kCache) && !el->TestBit(TStreamerElement::kWrite))) {
               // skip
            } else if (el != 0) {
               elements.AddLast(el);
            }
         }
         R__b.WriteObjectAny(&elements, TObjArray::Class(), kFALSE);
      }
      R__b.ClassEnd(TStreamerInfo::Class());
      R__b.SetByteCount(R__c, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Mark the classindex of the current file as using this TStreamerInfo.
/// This function is deprecated and its functionality is now done by
/// the overloads of TBuffer::TagStreamerInfo.

void TStreamerInfo::TagFile(TFile *file)
{
   if (file) {
      // If the value of the atomic is kFALSE (equal to expected), change its value
      // to kTRUE and return true. Leave it as it is otherwise and return false.
      static std::atomic<Bool_t> onlyonce(kFALSE);
      Bool_t expected = kFALSE;
      if (onlyonce.compare_exchange_strong(expected,kTRUE)) {
         Warning("TagFile","This function is deprecated, use TBuffer::TagStreamerInfo instead");
      }
      TArrayC *cindex = file->GetClassIndex();
      Int_t nindex = cindex->GetSize();
      if (fNumber < 0 || fNumber >= nindex) {
         Error("TagFile","StreamerInfo: %s number: %d out of range[0,%d] in file: %s",
               GetName(),fNumber,nindex,file->GetName());
         return;
      }
      if (cindex->fArray[fNumber] == 0) {
         cindex->fArray[0]       = 1;
         cindex->fArray[fNumber] = 1;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

#ifdef DOLOOP
#undef DOLOOP
#endif
#define DOLOOP for (k = 0, pointer = arr[0]; k < narr; pointer = arr[++k])

namespace {
   static void PrintCR(int j,Int_t aleng, UInt_t ltype)
   {
      if (j == aleng-1) printf("\n");
      else {
         printf(", ");
         if (j%ltype == ltype-1) printf("\n                    ");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///  print value of element  in object at pointer, type atype, leng aleng or *count
///  The function may be called in two ways:
///    -method1  len < 0
///           i is assumed to be the TStreamerElement number i in StreamerInfo
///    -method2  len >= 0
///           i is the type
///           address of variable is directly pointer.
///           len is the number of elements to be printed starting at pointer.

void TStreamerInfo::PrintValueAux(char *ladd, Int_t atype, TStreamerElement *aElement, Int_t aleng, Int_t *count)
{
   int j;

   //assert(!((kOffsetP + kChar) <= atype && atype <= (kOffsetP + kBool) && count == 0));
   switch (atype) {
      // basic types
      case kBool:              {Bool_t    *val = (Bool_t*   )ladd; printf("%d" ,*val);  break;}
      case kChar:              {Char_t    *val = (Char_t*   )ladd; printf("%d" ,*val);  break;}
      case kShort:             {Short_t   *val = (Short_t*  )ladd; printf("%d" ,*val);  break;}
      case kInt:               {Int_t     *val = (Int_t*    )ladd; printf("%d" ,*val);  break;}
      case kLong:              {Long_t    *val = (Long_t*   )ladd; printf("%ld",*val);  break;}
      case kLong64:            {Long64_t  *val = (Long64_t* )ladd; printf("%lld",*val);  break;}
      case kFloat:             {Float_t   *val = (Float_t*  )ladd; printf("%f" ,*val);  break;}
      case kFloat16:           {Float_t   *val = (Float_t*  )ladd; printf("%f" ,*val);  break;}
      case kDouble:            {Double_t  *val = (Double_t* )ladd; printf("%g" ,*val);  break;}
      case kDouble32:          {Double_t  *val = (Double_t* )ladd; printf("%g" ,*val);  break;}
      case kUChar:             {UChar_t   *val = (UChar_t*  )ladd; printf("%u" ,*val);  break;}
      case kUShort:            {UShort_t  *val = (UShort_t* )ladd; printf("%u" ,*val);  break;}
      case kUInt:              {UInt_t    *val = (UInt_t*   )ladd; printf("%u" ,*val);  break;}
      case kULong:             {ULong_t   *val = (ULong_t*  )ladd; printf("%lu",*val);  break;}
      case kULong64:           {ULong64_t *val = (ULong64_t*)ladd; printf("%llu",*val);  break;}
      case kBits:              {UInt_t    *val = (UInt_t*   )ladd; printf("%d" ,*val);  break;}

         // array of basic types  array[8]
      case kOffsetL + kBool:    {Bool_t    *val = (Bool_t*   )ladd; for(j=0;j<aleng;j++) { printf("%c " ,(char)val[j]); PrintCR(j,aleng,20); } break;}
      case kOffsetL + kChar:    {Char_t    *val = (Char_t*   )ladd; for(j=0;j<aleng;j++) { printf("%c " ,val[j]); PrintCR(j,aleng,20); } break;}
      case kOffsetL + kShort:   {Short_t   *val = (Short_t*  )ladd; for(j=0;j<aleng;j++) { printf("%d " ,val[j]); PrintCR(j,aleng,10); } break;}
      case kOffsetL + kInt:     {Int_t     *val = (Int_t*    )ladd; for(j=0;j<aleng;j++) { printf("%d " ,val[j]); PrintCR(j,aleng,10); } break;}
      case kOffsetL + kLong:    {Long_t    *val = (Long_t*   )ladd; for(j=0;j<aleng;j++) { printf("%ld ",val[j]); PrintCR(j,aleng, 5); } break;}
      case kOffsetL + kLong64:  {Long64_t  *val = (Long64_t* )ladd; for(j=0;j<aleng;j++) { printf("%lld ",val[j]);PrintCR(j,aleng, 5); } break;}
      case kOffsetL + kFloat:   {Float_t   *val = (Float_t*  )ladd; for(j=0;j<aleng;j++) { printf("%f " ,val[j]); PrintCR(j,aleng, 5); } break;}
      case kOffsetL + kFloat16: {Float_t   *val = (Float_t*  )ladd; for(j=0;j<aleng;j++) { printf("%f " ,val[j]); PrintCR(j,aleng, 5); } break;}
      case kOffsetL + kDouble:  {Double_t  *val = (Double_t* )ladd; for(j=0;j<aleng;j++) { printf("%g " ,val[j]); PrintCR(j,aleng, 5); } break;}
      case kOffsetL + kDouble32:{Double_t  *val = (Double_t* )ladd; for(j=0;j<aleng;j++) { printf("%g " ,val[j]); PrintCR(j,aleng, 5); } break;}
      case kOffsetL + kUChar:   {UChar_t   *val = (UChar_t*  )ladd; for(j=0;j<aleng;j++) { printf("%u " ,val[j]); PrintCR(j,aleng,20); } break;}
      case kOffsetL + kUShort:  {UShort_t  *val = (UShort_t* )ladd; for(j=0;j<aleng;j++) { printf("%u " ,val[j]); PrintCR(j,aleng,10); } break;}
      case kOffsetL + kUInt:    {UInt_t    *val = (UInt_t*   )ladd; for(j=0;j<aleng;j++) { printf("%u " ,val[j]); PrintCR(j,aleng, 5); } break;}
      case kOffsetL + kULong:   {ULong_t   *val = (ULong_t*  )ladd; for(j=0;j<aleng;j++) { printf("%lu ",val[j]); PrintCR(j,aleng, 5); } break;}
      case kOffsetL + kULong64: {ULong64_t *val = (ULong64_t*)ladd; for(j=0;j<aleng;j++) { printf("%llu ",val[j]);PrintCR(j,aleng, 5); } break;}
      case kOffsetL + kBits:    {UInt_t    *val = (UInt_t*   )ladd; for(j=0;j<aleng;j++) { printf("%d " ,val[j]); PrintCR(j,aleng, 5); } break;}

         // pointer to an array of basic types  array[n]
      case kOffsetP + kBool:    {Bool_t   **val = (Bool_t**  )ladd; for(j=0;j<*count;j++) { printf("%d " ,(*val)[j]);  PrintCR(j,aleng,20); } break;}
      case kOffsetP + kChar:    {Char_t   **val = (Char_t**  )ladd; for(j=0;j<*count;j++) { printf("%d " ,(*val)[j]);  PrintCR(j,aleng,20); } break;}
      case kOffsetP + kShort:   {Short_t  **val = (Short_t** )ladd; for(j=0;j<*count;j++) { printf("%d " ,(*val)[j]);  PrintCR(j,aleng,10); } break;}
      case kOffsetP + kInt:     {Int_t    **val = (Int_t**   )ladd; for(j=0;j<*count;j++) { printf("%d " ,(*val)[j]);  PrintCR(j,aleng,10); } break;}
      case kOffsetP + kLong:    {Long_t   **val = (Long_t**  )ladd; for(j=0;j<*count;j++) { printf("%ld ",(*val)[j]);  PrintCR(j,aleng, 5); } break;}
      case kOffsetP + kLong64:  {Long64_t **val = (Long64_t**)ladd; for(j=0;j<*count;j++) { printf("%lld ",(*val)[j]); PrintCR(j,aleng, 5); } break;}
      case kOffsetP + kFloat:   {Float_t  **val = (Float_t** )ladd; for(j=0;j<*count;j++) { printf("%f " ,(*val)[j]);  PrintCR(j,aleng, 5); } break;}
      case kOffsetP + kFloat16: {Float_t  **val = (Float_t** )ladd; for(j=0;j<*count;j++) { printf("%f " ,(*val)[j]);  PrintCR(j,aleng, 5); } break;}
      case kOffsetP + kDouble:  {Double_t **val = (Double_t**)ladd; for(j=0;j<*count;j++) { printf("%g " ,(*val)[j]);  PrintCR(j,aleng, 5); } break;}
      case kOffsetP + kDouble32:{Double_t **val = (Double_t**)ladd; for(j=0;j<*count;j++) { printf("%g " ,(*val)[j]);  PrintCR(j,aleng, 5); } break;}
      case kOffsetP + kUChar:   {UChar_t  **val = (UChar_t** )ladd; for(j=0;j<*count;j++) { printf("%u " ,(*val)[j]);  PrintCR(j,aleng,20); } break;}
      case kOffsetP + kUShort:  {UShort_t **val = (UShort_t**)ladd; for(j=0;j<*count;j++) { printf("%u " ,(*val)[j]);  PrintCR(j,aleng,10); } break;}
      case kOffsetP + kUInt:    {UInt_t   **val = (UInt_t**  )ladd; for(j=0;j<*count;j++) { printf("%u " ,(*val)[j]);  PrintCR(j,aleng, 5); } break;}
      case kOffsetP + kULong:   {ULong_t  **val = (ULong_t** )ladd; for(j=0;j<*count;j++) { printf("%lu ",(*val)[j]);  PrintCR(j,aleng, 5); } break;}
      case kOffsetP + kULong64: {ULong64_t**val = (ULong64_t**)ladd; for(j=0;j<*count;j++){ printf("%llu ",(*val)[j]); PrintCR(j,aleng, 5); } break;}

         // array counter //[n]
      case kCounter:            {Int_t *val    = (Int_t*)ladd;    printf("%d",*val);  break;}
         // char *
      case kCharStar:{
         char **val = (char**)ladd;
         if (*val) printf("%s",*val);
         break;
      }
      // Class *  derived from TObject with comment field  //->
      case kObjectp: {
         TObject **obj = (TObject**)(ladd);
         TStreamerObjectPointer *el = (TStreamerObjectPointer*)aElement;
         printf("(%s*)%lx",el ? el->GetClass()->GetName() : "unknown_type",(Long_t)(*obj));
         break;
      }

      // Class*   derived from TObject
      case kObjectP: {
         TObject **obj = (TObject**)(ladd);
         TStreamerObjectPointer *el = (TStreamerObjectPointer*)aElement;
         printf("(%s*)%lx",el ? el->GetClass()->GetName() : "unknown_type",(Long_t)(*obj));
         break;
      }

      // Class    derived from TObject
      case kObject:  {
         TObject *obj = (TObject*)(ladd);
         printf("%s",obj->GetName());
         break;
      }

      // Special case for TString, TObject, TNamed
      case kTString: {
         TString *st = (TString*)(ladd);
         printf("%s",st->Data());
         break;
      }
      case kTObject: {
         TObject *obj = (TObject*)(ladd);
         printf("%s",obj->GetName());
         break;
      }
      case kTNamed:  {
         TNamed *named = (TNamed*) (ladd);
         printf("%s/%s",named->GetName(),named->GetTitle());
         break;
      }

      // Class *  not derived from TObject with comment field  //->
      case kAnyp:    {
         TObject **obj = (TObject**)(ladd);
         TStreamerObjectAnyPointer *el = (TStreamerObjectAnyPointer*)aElement;
         printf("(%s*)0x%lx",el ? el->GetClass()->GetName() : "unknown_type",(Long_t)(*obj));
         break;
      }

      // Class*   not derived from TObject
      case kAnyP:    {
         TObject **obj = (TObject**)(ladd);
         TStreamerObjectAnyPointer *el = (TStreamerObjectAnyPointer*)aElement;
         printf("(%s*)0x%lx",el ? el->GetClass()->GetName() : "unknown_type",(Long_t)(*obj));
         break;
      }
      // Any Class not derived from TObject
      case kOffsetL + kObjectp:
      case kOffsetL + kObjectP:
      case kAny:     {
         printf("printing kAny case (%d)",atype);
//         if (aElement) {
//            TMemberStreamer *pstreamer = aElement->GetStreamer();
//            if (pstreamer == 0) {
//               //printf("ERROR, Streamer is null\n");
//               //aElement->ls();
//               break;
//            }
//            //(*pstreamer)(b,ladd,0);
//         }
         break;
      }
      // Base Class
      case kBase:    {
         printf("printing kBase case (%d)",atype);
         //aElement->ReadBuffer(b,pointer);
         break;
      }

      case kOffsetL + kObject:
      case kOffsetL + kTString:
      case kOffsetL + kTObject:
      case kOffsetL + kTNamed:
      case kStreamer: {
         printf("printing kStreamer case (%d)",atype);
//         TMemberStreamer *pstreamer = aElement->GetStreamer();
//         if (pstreamer == 0) {
//            //printf("ERROR, Streamer is null\n");
//            //aElement->ls();
//            break;
//         }
//         //UInt_t start,count;
//         //b.ReadVersion(&start, &count);
//         //(*pstreamer)(b,ladd,0);
//         //b.CheckByteCount(start,count,IsA());
         break;
      }

      case kStreamLoop: {
         printf("printing kStreamLoop case (%d)",atype);
//         TMemberStreamer *pstreamer = aElement->GetStreamer();
//         if (pstreamer == 0) {
//            //printf("ERROR, Streamer is null\n");
//            //aElement->ls();
//            break;
//         }
         //Int_t *counter = (Int_t*)(count);
         //UInt_t start,count;
         ///b.ReadVersion(&start, &count);
         //(*pstreamer)(b,ladd,*counter);
         //b.CheckByteCount(start,count,IsA());
         break;
      }
      case kSTL: {
         if (aElement) {
            static TClassRef stringClass("string");
            if (ladd && aElement->GetClass() == stringClass) {
               std::string *st = (std::string*)(ladd);
               printf("%s",st->c_str());
            } else {
               printf("(%s*)0x%lx",aElement->GetClass()->GetName(),(Long_t)(ladd));
            }
         } else {
            printf("(unknown_type*)0x%lx",(Long_t)(ladd));
         }
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///function called by the TClass constructor when replacing an emulated class
///by the real class

void TStreamerInfo::Update(const TClass *oldcl, TClass *newcl)
{
   TStreamerElement *element;
   TIter nextElement(GetElements());
   while ((element = (TStreamerElement*)nextElement())) {
      element->Update(oldcl,newcl);
   }
   for (Int_t i=0;i < fNslots;i++) {
      fComp[i].Update(oldcl,newcl);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update the TClass pointer cached in this object.

void TStreamerInfo::TCompInfo::Update(const TClass *oldcl, TClass *newcl)
{
   if (fType != -1) {
      if (fClass == oldcl || strcmp(fClassName, newcl->GetName()) == 0)
         fClass = newcl;
      else if (fClass == 0 && TClassTable::GetDict(fClassName))
         fClass = TClass::GetClass(fClassName);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Generate emulated collection proxy for a given class.

TVirtualCollectionProxy*
TStreamerInfo::GenEmulatedProxy(const char* class_name, Bool_t silent)
{
   return TCollectionProxyFactory::GenEmulatedProxy(class_name, silent);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate emulated class streamer for a given collection class.

TClassStreamer*
TStreamerInfo::GenEmulatedClassStreamer(const char* class_name, Bool_t silent)
{
   return TCollectionProxyFactory::GenEmulatedClassStreamer(class_name, silent);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate proxy from static functions.

TVirtualCollectionProxy*
TStreamerInfo::GenExplicitProxy( const ::ROOT::TCollectionProxyInfo &info, TClass *cl )
{
   return TCollectionProxyFactory::GenExplicitProxy(info, cl);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate class streamer from static functions.

TClassStreamer*
TStreamerInfo::GenExplicitClassStreamer( const ::ROOT::TCollectionProxyInfo &info, TClass *cl )
{
   return TCollectionProxyFactory::GenExplicitClassStreamer(info, cl);
}

//
// Utility functions
//
static TStreamerElement* R__CreateEmulatedElement(const char *dmName, const std::string &dmFull, Int_t offset, bool silent)
{
   // Create a TStreamerElement for the type 'dmFull' and whose data member name is 'dmName'.

   TString s1( TClassEdit::ShortType(dmFull.c_str(),0) );
   TString dmType( TClassEdit::ShortType(dmFull.c_str(),1) );
   Bool_t dmIsPtr = (s1 != dmType);
   const char *dmTitle = "Emulation";

   TDataType *dt = gROOT->GetType(dmType);
   if (dt && dt->GetType() > 0 ) {  // found a basic type
      Int_t dsize,dtype;
      dtype = dt->GetType();
      dsize = dt->Size();
      if (dmIsPtr && dtype != kCharStar) {
         if (!silent)
            Error("Pair Emulation Building","%s is not yet supported in pair emulation",
                  dmFull.c_str());
         return 0;
      } else {
         TStreamerElement *el = new TStreamerBasicType(dmName,dmTitle,offset,dtype,dmFull.c_str());
         el->SetSize(dsize);
         return el;
      }
   } else {

      static const char *full_string_name = "basic_string<char,char_traits<char>,allocator<char> >";
      if (strcmp(dmType,"string") == 0 || strcmp(dmType,"std::string") == 0 || strcmp(dmType,full_string_name)==0 ) {
         return new TStreamerSTLstring(dmName,dmTitle,offset,dmFull.c_str(),dmIsPtr);
      }
      if (TClassEdit::IsSTLCont(dmType)) {
         return new TStreamerSTL(dmName,dmTitle,offset,dmFull.c_str(),dmFull.c_str(),dmIsPtr);
      }
      TClass *clm = TClass::GetClass(dmType);
      if (!clm) {
         auto enumdesc = TEnum::GetEnum(dmType, TEnum::kNone);
         if (enumdesc) {
            auto dtype = enumdesc->GetUnderlyingType();
            auto el = new TStreamerBasicType(dmName, dmTitle, offset, dtype, dmFull.c_str());
            auto datatype = TDataType::GetDataType(dtype);
            if (datatype)
               el->SetSize(datatype->Size());
            else
               el->SetSize(sizeof(int)); // Default size of enums.
            return el;
         }
         return nullptr;
      }
      if (clm->GetState() <= TClass::kForwardDeclared)
         return nullptr;
      // a pointer to a class
      if ( dmIsPtr ) {
         if (clm->IsTObject()) {
            return new TStreamerObjectPointer(dmName,dmTitle,offset,dmFull.c_str());
         } else {
            return new TStreamerObjectAnyPointer(dmName,dmTitle,offset,dmFull.c_str());
         }
      }
      // a class
      if (clm->IsTObject()) {
         return new TStreamerObject(dmName,dmTitle,offset,dmFull.c_str());
      } else if(clm == TString::Class() && !dmIsPtr) {
         return new TStreamerString(dmName,dmTitle,offset);
      } else {
         return new TStreamerObjectAny(dmName,dmTitle,offset,dmFull.c_str());
      }
   }
}

// \brief Generate the TClass and TStreamerInfo for the requested pair.
// This creates a TVirtualStreamerInfo for the pair and trigger the BuildCheck/Old to
// provoke the creation of the corresponding TClass.  This relies on the dictionary for
// std::pair<const int, int> to already exist (or the interpreter information being available)
// as it is used as a template.
TVirtualStreamerInfo *TStreamerInfo::GenerateInfoForPair(const std::string &firstname, const std::string &secondname, bool silent, size_t hint_pair_offset, size_t hint_pair_size)
{
   // Generate a TStreamerInfo for a std::pair<fname,sname>
   // This TStreamerInfo is then used as if it was read from a file to generate
   // and emulated TClass.

   if (hint_pair_offset && hint_pair_offset >= hint_pair_size) {
      const char *msg = "problematic";
      if (hint_pair_offset == hint_pair_size)
         msg = "the same";
      else if (hint_pair_offset > hint_pair_size) {
         if (hint_pair_size == 0)
            msg = "not specified";
         else
            msg = "smaller";
      }
      Error("GenerateInfoForPair",
            "Called with inconsistent offset and size. For \"std::pair<%s,%s>\" the requested offset is %ld but the size is %s (%ld)",
            firstname.c_str(), secondname.c_str(), (long)hint_pair_offset, msg, (long)hint_pair_offset);
      return nullptr;
   }
   TStreamerInfo *i = (TStreamerInfo*)TClass::GetClass("pair<const int,int>")->GetStreamerInfo()->Clone();
   std::string pname = "pair<" + firstname + "," + secondname;
   pname += (pname[pname.length()-1]=='>') ? " >" : ">";
   i->SetName(pname.c_str());
   i->SetClass(nullptr);
   i->GetElements()->Delete();
   TStreamerElement *fel = R__CreateEmulatedElement("first", firstname, 0, silent);
   Int_t size = 0;
   if (fel) {
      i->GetElements()->Add( fel );

      size = fel->GetSize();
      Int_t sp = sizeof(void *);
      //align the non-basic data types (required on alpha and IRIX!!)
      if (size%sp != 0) size = size - size%sp + sp;
   } else {
      delete i;
      return 0;
   }
   if (hint_pair_offset)
      size = hint_pair_offset;
   TStreamerElement *second = R__CreateEmulatedElement("second", secondname, size, silent);
   if (second) {
      i->GetElements()->Add( second );
   } else {
      delete i;
      return 0;
   }
   Int_t oldlevel = gErrorIgnoreLevel;
   // Hide the warning about the missing pair dictionary.
   gErrorIgnoreLevel = kError;
   i->BuildCheck(nullptr, kFALSE); // Skipping the loading part (it would leads to infinite recursion on this very routine)
   gErrorIgnoreLevel = oldlevel;
   // In the state emulated, BuildOld would recalculate the offset and undo the offset update.
   // Note: we should consider adding a new state just for this (the hints indicates that we are mapping a compiled class but
   // then we would have to investigate all use of the state with <= and >= condition to make sure they are still appropriate).
   if (hint_pair_size) {
      i->GetClass()->SetClassSize(hint_pair_size);
      i->GetClass()->fIsSyntheticPair = kTRUE;
   }

   i->BuildOld();

   if (hint_pair_size)
      i->GetClass()->SetClassSize(hint_pair_size);
   return i;
}

TVirtualStreamerInfo *TStreamerInfo::GenerateInfoForPair(const std::string &pairclassname, bool silent, size_t hint_pair_offset, size_t hint_pair_size)
{
   const static int pairlen = strlen("pair<");
   if (pairclassname.compare(0, pairlen, "pair<") != 0) {
      if (!silent)
         Error("GenerateInfoForPair", "The class name passed is not a pair: %s", pairclassname.c_str());
      return nullptr;
   }

   std::vector<std::string> inside;
   int nested = 0;
   int num = TClassEdit::GetSplit(pairclassname.c_str(), inside, nested);
   if (num != 4) {
      if (!silent)
         Error("GenerateInfoForPair", "Could not find the pair arguments in %s", pairclassname.c_str());
      return nullptr;
   }

   return GenerateInfoForPair(inside[1], inside[2], silent, hint_pair_offset, hint_pair_size);
}
