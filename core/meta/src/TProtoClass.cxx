// @(#)root/meta:$
// Author: Axel Naumann 2014-05-02

/*************************************************************************
 * Copyright (C) 1995-2014, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TProtoClass
Persistent version of a TClass.
*/

#include "TProtoClass.h"

#include "TBaseClass.h"
#include "TClass.h"
#include "TClassEdit.h"
#include "TDataMember.h"
#include "TEnum.h"
#include "TInterpreter.h"
#include "TList.h"
#include "TListOfDataMembers.h"
#include "TListOfEnums.h"
#include "TRealData.h"
#include "TError.h"
#include "TVirtualCollectionProxy.h"

#include <cassert>
#include <unordered_map>

////////////////////////////////////////////////////////////////////////////////
/// Initialize a TProtoClass from a TClass.

TProtoClass::TProtoClass(TClass* cl):
   TNamed(*cl), fBase(cl->GetListOfBases()),
   fEnums(cl->GetListOfEnums()), fSizeof(cl->Size()), fCheckSum(cl->fCheckSum),
   fCanSplit(cl->fCanSplit), fStreamerType(cl->fStreamerType), fProperty(cl->fProperty),
   fClassProperty(cl->fClassProperty)
{
   if (cl->Property() & kIsNamespace){
      //fData=new TListOfDataMembers();
      fEnums=nullptr;
      //fPRealData=nullptr;
      fOffsetStreamer=0;
      return;
   }
   TListOfEnums *enums = dynamic_cast<TListOfEnums*>(fEnums);
   if (enums && !enums->fIsLoaded) {
      // Make sure all the enum information is loaded
      enums->Load();
   }
   // initialize list of data members (fData)
   TList * dataMembers = cl->GetListOfDataMembers();
   if (dataMembers && dataMembers->GetSize() > 0) {
      fData.reserve(dataMembers->GetSize() );
      for (auto * obj : *dataMembers) {
         TDataMember * dm = dynamic_cast<TDataMember*>(obj);
         fData.push_back(dm);
      }
   }

   fPRealData.reserve(100);
   class DepClassDedup {
      std::vector<TString> &fDepClasses;
      std::unordered_map<std::string, int> fDepClassIdx;
   public:
      DepClassDedup(std::vector<TString> &depClasses): fDepClasses(depClasses)
      {
         R__ASSERT(fDepClasses.empty() && "Expected fDepClasses to be empty before fililng it!");
      }

      ~DepClassDedup()
      {
         if (fDepClasses.size() != fDepClassIdx.size())
            ::Error("TProtoClass::DepClassDedup::~DepClassDedup",
                    "Mismatching size of fDepClasses and index map! Please report.");
      }

      int GetIdx(const char *name) {
         auto itins = fDepClassIdx.insert({name, fDepClasses.size()});
         if (itins.second) {
            fDepClasses.emplace_back(name);
         }
         return itins.first->second;
      }
   } depClassDedup(fDepClasses);

   if (!cl->GetCollectionProxy()) {
      // Build the list of RealData before we access it:
      cl->BuildRealData(0, true /*isTransient*/);
      // The data members are ordered as follows:
      // - this class's data members,
      // - foreach base: base class's data members.
      for (auto realDataObj: *cl->GetListOfRealData()) {
         TRealData *rd = (TRealData*)realDataObj;
         if (!rd->GetDataMember())
            continue;
         TProtoRealData protoRealData(rd);

         if (TClass* clRD = rd->GetDataMember()->GetClass())
            protoRealData.fClassIndex = depClassDedup.GetIdx(clRD->GetName());

         protoRealData.SetFlag(TProtoRealData::kIsTransient, rd->TestBit(TRealData::kTransient));

         fPRealData.emplace_back(protoRealData);
      }

      // if (gDebug > 2) {
         // for (const auto &data : fPRealData) {
            // const auto classType = dataPtr->IsA();
            // const auto dataName = data.fName;
            // const auto dataClass = data.fClass;
            // Info("TProtoClass","Data is a protorealdata: %s - class %s - transient %d", dataName.Data(),dataClass.Data(),data.fIsTransient);
            //if (!dataClass.IsNull()
            // if (classType == TProtoRealData::Class())
            //    Info("TProtoClass","Data is a protorealdata: %s", dataPtrName);
            // if (classType == TObjString::Class())
            //    Info("TProtoClass","Data is a objectstring: %s", dataPtrName);
            // if (dataPtr->TestBit(TRealData::kTransient))
            //    Info("TProtoClass","And is transient");
         // }
      // }
   } else if (cl->GetCollectionProxy()->GetProperties() & TVirtualCollectionProxy::kIsEmulated) {
      // The collection proxy is emulated has the wrong size.
      if (cl->HasInterpreterInfo())
         fSizeof = gCling->ClassInfo_Size(cl->GetClassInfo());
      else
         fSizeof = -1;
   }

   cl->CalculateStreamerOffset();
   fOffsetStreamer = cl->fOffsetStreamer;
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TProtoClass::~TProtoClass()
{
   Delete();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete the containers that are usually owned by their TClass.
/// if (fPRealData) fPRealData->Delete(opt);
/// delete fPRealData; fPRealData = 0;

void TProtoClass::Delete(Option_t* opt /*= ""*/) {
   if (fBase) fBase->Delete(opt);
   delete fBase; fBase = 0;

   for (auto dm: fData) {
      delete dm;
   }

   if (fEnums) fEnums->Delete(opt);
   delete fEnums; fEnums = 0;

   if (gErrorIgnoreLevel==-2) printf("Delete the protoClass %s \n",GetName());
}

////////////////////////////////////////////////////////////////////////////////
/// Move data from this TProtoClass into cl.
/// Returns 'false' if nothing was done.  This can happen in the case where
/// there is more than one dictionary for the same entity.  Note having
/// duplicate dictionary is acceptable for namespace or STL collections.

Bool_t TProtoClass::FillTClass(TClass* cl) {
   if (cl->fRealData || cl->fBase.load() || cl->fData.load() || cl->fEnums.load() || cl->fCanSplit >= 0 ||
       cl->fProperty != (-1)) {

      if (cl->fState == TClass::kHasTClassInit)
         // The class has dictionary, has gone through some initialization and is now being requested
         // to be filled by a TProtoClass.
         // This can be due to:
         //   (a) A duplicate dictionary for a class (with or without a rootpcm associated with)
         //   (b) The TClass was created before the registration of the rootpcm ** and ** it was
         //       attempted to be used before this registration

         // This is technically an error
         // but we either already warned that there is a 2nd dictionary for the class (in TClassTable::Add)
         // or this is the same (but now emptied) TProtoClass instance as before.
         // We return false, since we are doing no actual change to the TClass instance and thus
         // if a caller was hoping for 'improvement' in the state of the TClass instance, it did not
         // happen.
         return kFALSE;

      if (cl->GetCollectionType() != ROOT::kNotSTL) {
         // We are in the case of collection, duplicate dictionary are allowed
         // (and even somewhat expected since they can be auto asked for).
         // They do not always have a TProtoClass for them.  In particular
         // the one that are pre-generated in the ROOT build (in what was once
         // called the cintdlls) do not have a pcms, neither does vector<string>
         // which is part of libCore proper.
         if (gDebug > 0)
            Info("FillTClass", "Returning w/o doing anything. %s is a STL collection.",cl->GetName());
         return kFALSE;
      }
      if (cl->fProperty != -1 && (cl->fProperty & kIsNamespace)) {
         if (gDebug > 0)
            Info("FillTClass", "Returning w/o doing anything. %s is a namespace.",cl->GetName());
         return kFALSE;
      }
      Error("FillTClass", "TClass %s already initialized!", cl->GetName());
      return kFALSE;
   }
   if (cl->fHasRootPcmInfo) {
      Fatal("FillTClass", "Filling TClass %s a second time but none of the info is in the TClass instance ... ", cl->GetName());
   }
   if (gDebug > 1) Info("FillTClass","Loading TProtoClass for %s - %s",cl->GetName(),GetName());

   if (fPRealData.size() > 0) {

      // A first loop to retrieve the mother classes before starting to
      // fill this TClass instance. This is done in order to avoid recursions
      // for example in presence of daughter and mother class present in two
      // dictionaries compiled in two different libraries which are not linked
      // one with each other.
      for (auto element: fPRealData) {
         //if (element->IsA() == TObjString::Class()) {
         if (element.IsAClass() ) {
            if (gDebug > 1) Info("","Treating beforehand mother class %s",GetClassName(element.fClassIndex));
            TInterpreter::SuspendAutoParsing autoParseRaii(gInterpreter);

            TClass::GetClass(GetClassName(element.fClassIndex));
         }
      }
   }


   //this->Dump();

   // Copy only the TClass bits.
   // not bit 13 and below and not bit 24 and above, just Bits 14 - 23
   UInt_t newbits = TestBits(0x00ffc000);
   cl->ResetBit(0x00ffc000);
   cl->SetBit(newbits);

   cl->fName  = this->fName;
   cl->fTitle = this->fTitle;
   cl->fBase = fBase;

   // fill list of data members in TClass
   //if (cl->fData) { cl->fData->Delete(); delete cl->fData;  }
   cl->fData = new TListOfDataMembers(fData);
   // for (auto * dataMember : fData) {
   //    //printf("add data member for class %s - member %s \n",GetName(), dataMember->GetName() );
   //    cl->fData->Add(dataMember);
   // }
   // // set loaded bit to true to avoid re-loading the data members
   // cl->fData->SetIsLoaded();*

   //cl->fData = (TListOfDataMembers*)fData;

   // The TDataMember were passed along.
   fData.clear();

   // We need to fill enums one by one to initialise the internal map which is
   // transient
   {
      auto temp = new TListOfEnums();
      if (fEnums) {
         for (TObject* enumAsTObj : *fEnums){
            temp->Add((TEnum*) enumAsTObj);
         }
         // We did not transfer the container itself, let remove it from memory without deleting its content.
         fEnums->Clear();
         delete fEnums;
         fEnums = nullptr;
      }
      cl->fEnums = temp;
   }

   if (cl->fSizeof != -1 && cl->fSizeof != fSizeof) {
      Error("FillTClass",
            "For %s the sizeof provided by GenerateInitInstance (%d) is different from the one provided by TProtoClass (%d)",
            cl->GetName(), cl->fSizeof, fSizeof);
   } else
      cl->fSizeof = fSizeof;
   cl->fCheckSum = fCheckSum;
   cl->fCanSplit = fCanSplit;
   cl->fProperty = fProperty;
   cl->fClassProperty = fClassProperty;
   cl->fStreamerType = fStreamerType;

   // Update pointers to TClass
   if (cl->fBase.load()) {
      for (auto base: *cl->fBase) {
         ((TBaseClass*)base)->SetClass(cl);
      }
   }
   if (cl->fData) {
      for (auto dm: *cl->fData) {
         ((TDataMember*)dm)->SetClass(cl);
      }
      ((TListOfDataMembers*)cl->fData)->SetClass(cl);
   }
   if (cl->fEnums.load()) {
      for (auto en: *cl->fEnums) {
         ((TEnum*)en)->SetClass(cl);
      }
      ((TListOfEnums*)cl->fEnums)->SetClass(cl);
   }


   TClass* currentRDClass = cl;
   TRealData * prevRealData = 0;
   int prevLevel = 0;
   bool first = true;
   if (fPRealData.size()  > 0) {
      size_t element_next_idx = 0;
      for (auto element: fPRealData) {
         ++element_next_idx;
         //if (element->IsA() == TObjString::Class()) {
         if (element.IsAClass() ) {
            // We now check for the TClass entry, w/o loading. Indeed we did that above.
            // If the class is not found, it means that really it was not selected and we
            // replace it with an empty placeholder with the status of kForwardDeclared.
            // Interactivity will be of course possible but if IO is attempted, a warning
            // will be issued.
            TInterpreter::SuspendAutoParsing autoParseRaii(gInterpreter);

            const char *classname = GetClassName(element.fClassIndex);

            // Disable autoparsing which might be triggered by the use of ResolvedTypedef
            // and the fallback new TClass() below.
            currentRDClass = TClass::GetClass(classname, false /* Load */ );
            //printf("element is a class - name %s  - index %d  %s \n ",currentRDClass->GetName(), element.fClassIndex, GetClassName(element.fClassIndex) );
            if (!currentRDClass && !element.TestFlag(TProtoRealData::kIsTransient)) {

               if (TClassEdit::IsStdPair(classname) && element.fDMIndex == 0 && fPRealData.size() > element_next_idx) {
                  size_t hint_offset = fPRealData[element_next_idx].fOffset - element.fOffset;
                  size_t hint_size = 0;
                  // Now find the size.
                  size_t end = element_next_idx + 1;
                  while (end < fPRealData.size() && fPRealData[end].fLevel > element.fLevel)
                     ++end;
                  if (end < fPRealData.size()) {
                     hint_size = fPRealData[end].fOffset - element.fOffset;
                  } else {
                     hint_size = fSizeof - element.fOffset;
                  }
                  currentRDClass = TClass::GetClass(classname, true, false, hint_offset, hint_size);
               }
               if (!currentRDClass) {
                  if (gDebug > 1)
                     Info("FillTClass()",
                          "Cannot find TClass for %s; Creating an empty one in the kForwardDeclared state.", classname);
                  currentRDClass = new TClass(classname, 1, TClass::kForwardDeclared, true /*silent*/);
               }
            }
         }
         //else {
         if (!currentRDClass) continue;
         //TProtoRealData* prd = (TProtoRealData*)element;
         // pass a previous real data only if depth

         if (TRealData* rd = element.CreateRealData(currentRDClass, cl,prevRealData, prevLevel)) {
            if (first) {
               //LM: need to do here because somehow fRealData is destroyed when calling TClass::GetListOfDataMembers()
               if (cl->fRealData) {
                  Info("FillTClass","Real data for class %s is not empty - make a new one",cl->GetName() );
                  delete cl->fRealData;
               }
               cl->fRealData = new TList(); // FIXME: this should really become a THashList!
               first = false;
            }

            cl->fRealData->AddLast(rd);
            prevRealData = rd;
            prevLevel = element.fLevel;

         }
         //}
      }
   }
   else {
      if (cl->fRealData) {
         Info("FillTClas","Real data for class %s is not empty - make a new one. Class has no Proto-realdata",cl->GetName() );
         delete cl->fRealData;
      }
      cl->fRealData = new TList(); // FIXME: this should really become a THashList!
   }

   cl->SetStreamerImpl();

   // set to zero in order not to delete when protoclass is deleted
   fBase = 0;
   //fData = 0;
   fEnums = 0;

   fPRealData.clear();
   fPRealData.shrink_to_fit();  // to reset the underlying allocate space

   // if (fPRealData) fPRealData->Delete();
   // delete fPRealData;
   // fPRealData = 0;

   cl->fHasRootPcmInfo = kTRUE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////

TProtoClass::TProtoRealData::TProtoRealData(const TRealData* rd):
   //TNamed(rd->GetDataMember()->GetName(), rd->GetName()),
   //TNamed(),
   //fName(rd->GetDataMember()->GetName()),
   //fTitle(rd->GetName()),
   fOffset(rd->GetThisOffset()),
   fDMIndex(-1),
   fLevel(0),
   fClassIndex(-1),
   fStatusFlag(0)
{
   TDataMember * dm = rd->GetDataMember();
   assert(rd->GetDataMember());
   TClass * cl = dm->GetClass();
   assert(cl != NULL);
   fDMIndex = DataMemberIndex(cl,dm->GetName());
   //printf("Index of data member %s for class %s is %d \n",dm->GetName(), cl->GetName() , fDMIndex);
   TString fullDataMemberName = rd->GetName(); // full data member name (e.g. fXaxis.fNbins)
   fLevel = fullDataMemberName.CountChar('.');

   if (fullDataMemberName.Contains("*") ) SetFlag(kIsPointer);

   // Initialize this from a TRealData object.
   SetFlag(kIsObject, rd->IsObject());
   SetFlag(kIsTransient, rd->TestBit(TRealData::kTransient) );
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor to pin vtable.
///if (gErrorIgnoreLevel==-2) printf("destroy real data %s - ",GetName());

TProtoClass::TProtoRealData::~TProtoRealData()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TRealData from this, with its data member coming from dmClass.
/// find data member from protoclass

TRealData* TProtoClass::TProtoRealData::CreateRealData(TClass* dmClass,
                                                       TClass* parent, TRealData *prevData, int prevLevel) const
{

   //TDataMember* dm = (TDataMember*)dmClass->GetListOfDataMembers()->FindObject(fName);
   TDataMember* dm = TProtoClass::FindDataMember(dmClass, fDMIndex);

   if (!dm && dmClass->GetState()!=TClass::kForwardDeclared && !dmClass->fIsSyntheticPair) {
      ::Error("CreateRealData",
            "Cannot find data member # %d of class %s for parent %s!", fDMIndex, dmClass->GetName(),
            parent->GetName());
      return nullptr;
   }

   // here I need to re-construct the realdata full name (e.g. fAxis.fNbins)

   TString realMemberName;
   // keep an empty name if data member is not found
   if (dm) realMemberName = dm->GetName();
   else if (dmClass->fIsSyntheticPair) {
      realMemberName = (fDMIndex == 0) ? "first" : "second";
   }
   if (TestFlag(kIsPointer) )
      realMemberName = TString("*")+realMemberName;
   else if (dm){
      if (dm->GetArrayDim() > 0) {
         // in case of array (like fMatrix[2][2] we need to add max index )
         // this only in case of it os not a pointer
         for (int idim = 0; idim < dm->GetArrayDim(); ++idim)
            realMemberName += TString::Format("[%d]",dm->GetMaxIndex(idim) );
      } else if (TClassEdit::IsStdArray(dm->GetTypeName())) {
         std::string typeNameBuf;
         Int_t ndim = dm->GetArrayDim();
         std::array<Int_t, 5> maxIndices; // 5 is the maximum supported in TStreamerElement::SetMaxIndex
         TClassEdit::GetStdArrayProperties(dm->GetTypeName(),
                                           typeNameBuf,
                                           maxIndices,
                                           ndim);
         for (Int_t idim = 0; idim < ndim; ++idim) {
            realMemberName += TString::Format("[%d]",maxIndices[idim] );
         }
      }
   }

   if (prevData && fLevel > 0 ) {
      if (fLevel-prevLevel == 1) // I am going down 1 level
         realMemberName = TString::Format("%s.%s",prevData->GetName(), realMemberName.Data() );
      else if (fLevel <= prevLevel) { // I am at the same level
                                      // need to strip out prev name
         std::string prevName = prevData->GetName();
         // we strip the prev data member name from the full name
         std::string parentName;
         for (int i = 0; i < prevLevel-fLevel+1; ++i) {
            parentName = prevName.substr(0, prevName.find_last_of(".") );
            prevName = parentName;
         }

         // now we need to add the current name
         realMemberName =  TString::Format("%s.%s",parentName.c_str(), realMemberName.Data() );
      }
   }

   //printf("adding new realdata for class %s : %s - %s   %d    %d   \n",dmClass->GetName(), realMemberName.Data(), dm->GetName(),fLevel, fDMIndex  );

   TRealData* rd = new TRealData(realMemberName, fOffset, dm);
   if (TestFlag(kIsTransient)) {
      rd->SetBit(TRealData::kTransient);
   }
   rd->SetIsObject(TestFlag(kIsObject) );
   return rd;
}

////////////////////////////////////////////////////////////////////////////////

Int_t TProtoClass::DataMemberIndex(TClass * cl, const char * name)
{
   TList * dmList = cl->GetListOfDataMembers();

   // we cannot use IndexOf because order is guaranteed only for non-static data member
   Int_t index = 0;
   for ( auto * obj : *dmList) {
      TDataMember * dm = (TDataMember *) obj;
      if (!dm ) continue;
      if (dm->Property() & kIsStatic) continue;
      if ( TString(dm->GetName()) == TString(name) )
         return index;
      index++;
   }
   ::Error("TProtoClass::DataMemberIndex","data member %s is not found in class %s",name, cl->GetName());
   dmList->ls();
   return -1;
}
////////////////////////////////////////////////////////////////////////////////

TDataMember * TProtoClass::FindDataMember(TClass * cl, Int_t index)
{
   TList * dmList = cl->GetListOfDataMembers(false);

   // we cannot use IndexOf because order is guaranteed only for non-static data member
   Int_t i = 0;
   for ( auto * obj : *dmList) {
      TDataMember * dm = (TDataMember *) obj;
      if (!dm ) continue;
      if (dm->Property() & kIsStatic) continue;
      if (i == index)
         return dm;
      i++;
   }
   if (cl->GetState()!=TClass::kForwardDeclared && !cl->fIsSyntheticPair)
      ::Error("TProtoClass::FindDataMember","data member with index %d is not found in class %s",index,cl->GetName());
   return nullptr;
}
