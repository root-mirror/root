// @(#)root/treeplayer:$Id$
// Author: Philippe Canal 01/06/2004

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers and al.        *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TFormLeafInfo
This class is a small helper class to implement reading a data member
on an object stored in a TTree.

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

#include "TFormLeafInfo.h"

#include "TROOT.h"
#include "TArrayI.h"
#include "TClonesArray.h"
#include "TError.h"
#include "TInterpreter.h"
#include "TLeafObject.h"
#include "TMethod.h"
#include "TMethodCall.h"
#include "TTree.h"
#include "TVirtualCollectionProxy.h"
#include "TClassEdit.h"


#define INSTANTIATE_READVAL(CLASS) \
   template Double_t CLASS::ReadValueImpl<Double_t>(char*, Int_t);  \
   template Long64_t CLASS::ReadValueImpl<Long64_t>(char*, Int_t);  \
   template LongDouble_t CLASS::ReadValueImpl<LongDouble_t>(char*, Int_t)  // no semicolon


#define INSTANTIATE_GETVAL(CLASS) \
   template Double_t CLASS::GetValueImpl<Double_t>(TLeaf*, Int_t); \
   template Long64_t CLASS::GetValueImpl<Long64_t>(TLeaf*, Int_t); \
   template LongDouble_t CLASS::GetValueImpl<LongDouble_t>(TLeaf*, Int_t)  // no semicolon

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfo::TFormLeafInfo(TClass* classptr, Long_t offset,
                             TStreamerElement* element) :
     fClass(classptr),fOffset(offset),fElement(element),
     fCounter(0), fNext(0),fMultiplicity(0)
{
   if (fClass) fClassName = fClass->GetName();
   if (fElement) {
      fElementName = fElement->GetName();
   }
}

////////////////////////////////////////////////////////////////////////////////
///Constructor.

TFormLeafInfo::TFormLeafInfo(const TFormLeafInfo& orig) : TObject(orig),fClass(orig.fClass),fOffset(orig.fOffset),fElement(orig.fElement),fCounter(0),fNext(0),fClassName(orig.fClassName),fElementName(orig.fElementName),fMultiplicity(orig.fMultiplicity)
{
   // Deep copy the pointers.
   if (orig.fCounter) fCounter = orig.fCounter->DeepCopy();
   if (orig.fNext) fNext = orig.fNext->DeepCopy();
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe assignment operator.

TFormLeafInfo &TFormLeafInfo::operator=(const TFormLeafInfo &other)
{
   TFormLeafInfo tmp(other);
   Swap(tmp);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////

void TFormLeafInfo::Swap(TFormLeafInfo& other)
{
   std::swap(fClass,other.fClass);
   std::swap(fOffset,other.fOffset);
   std::swap(fElement,other.fElement);
   std::swap(fCounter,other.fCounter);
   std::swap(fNext,other.fNext);
   TString tmp(fClassName);
   fClassName = other.fClassName;
   other.fClassName = tmp;

   tmp = fElementName;
   fElementName = other.fElementName;
   other.fElementName = tmp;
   std::swap(fMultiplicity,other.fMultiplicity);
}

////////////////////////////////////////////////////////////////////////////////
/// Make a complete copy of this FormLeafInfo and all its content.

TFormLeafInfo* TFormLeafInfo::DeepCopy() const
{
   return new TFormLeafInfo(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete this object and all its content

TFormLeafInfo::~TFormLeafInfo()
{
   delete fCounter;
   delete fNext;
}

////////////////////////////////////////////////////////////////////////////////
/// Increase the offset of this element.  This intended to be the offset
/// from the start of the object to which the data member belongs.

void TFormLeafInfo::AddOffset(Int_t offset, TStreamerElement* element)
{
   fOffset += offset;
   fElement = element;
   if (fElement ) {
      //         fElementClassOwnerName = cl->GetName();
      fElementName.Append(".").Append(element->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current length of the array.

Int_t TFormLeafInfo::GetArrayLength()
{
   Int_t len = 1;
   if (fNext) len = fNext->GetArrayLength();
   if (fElement) {
      Int_t elen = fElement->GetArrayLength();
      if (elen || fElement->IsA() == TStreamerBasicPointer::Class() )
         len *= fElement->GetArrayLength();
   }
   return len;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the class of the underlying data.

TClass* TFormLeafInfo::GetClass() const
{
   if (fNext) return fNext->GetClass();
   if (fElement) return fElement->GetClassPointer();
   return fClass;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the the location of the object pointed to.
/// Modify instance if the object is part of an array.

char* TFormLeafInfo::GetObjectAddress(TLeafElement* leaf, Int_t& instance)
{
   TBranchElement* branch = (TBranchElement*) leaf->GetBranch();
   Int_t id = branch->GetID();
   if (id < 0) {
      // Branch is a top-level branch.
      if (branch->GetTree()->GetMakeClass()) {
         // Branch belongs to a MakeClass tree.
         return branch->GetAddress();
      } else {
         return branch->GetObject();
      }
   }
   TStreamerInfo* info = branch->GetInfo();
   Int_t offset = 0;
   if (id > -1) {
      // Branch is *not* a top-level branch.
      offset = info->TStreamerInfo::GetElementOffset(id);
   }
   char* address = 0;
   // Branch is *not* a top-level branch.
   if (branch->GetTree()->GetMakeClass()) {
      // Branch belongs to a MakeClass tree.
      address = (char*) branch->GetAddress();
   } else {
      address = (char*) branch->GetObject();
   }
   char* thisobj = 0;
   if (!address) {
      // FIXME: This makes no sense, if the branch address is not set, then object will not be set either.
      thisobj = branch->GetObject();
   } else {
      Int_t type = -1;
      if (id > -1) {
         // Note this is somewhat slow
         type = info->TStreamerInfo::GetElement(id)->GetNewType();
      }
      switch (type) {
         case TStreamerInfo::kOffsetL + TStreamerInfo::kObjectp:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kObjectP:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kSTLp:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kAnyp:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kAnyP:
            Error("GetValuePointer", "Type (%d) not yet supported\n", type);
            break;

         case TStreamerInfo::kOffsetL + TStreamerInfo::kObject:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kAny:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kSTL:
         {
            // An array of objects.
            Int_t index;
            Int_t sub_instance;
            Int_t len = GetArrayLength();
            if (len) {
               index = instance / len;
               sub_instance = instance % len;
            } else {
               index = instance;
               sub_instance = 0;
            }
            thisobj = address + offset + (index * fClass->Size());
            instance = sub_instance;
            break;
         }

         case TStreamerInfo::kBase:
         case TStreamerInfo::kObject:
         case TStreamerInfo::kTString:
         case TStreamerInfo::kTNamed:
         case TStreamerInfo::kTObject:
         case TStreamerInfo::kAny:
         case TStreamerInfo::kSTL:
            // A single object.
            thisobj = address + offset;
            break;

         case TStreamerInfo::kBool:
         case TStreamerInfo::kChar:
         case TStreamerInfo::kUChar:
         case TStreamerInfo::kShort:
         case TStreamerInfo::kUShort:
         case TStreamerInfo::kInt:
         case TStreamerInfo::kUInt:
         case TStreamerInfo::kLong:
         case TStreamerInfo::kULong:
         case TStreamerInfo::kLong64:
         case TStreamerInfo::kULong64:
         case TStreamerInfo::kFloat:
         case TStreamerInfo::kFloat16:
         case TStreamerInfo::kDouble:
         case TStreamerInfo::kDouble32:
         case TStreamerInfo::kLegacyChar:
         case TStreamerInfo::kCounter:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kBool:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kChar:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUChar:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kShort:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUShort:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kInt:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUInt:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kLong:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kULong:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kLong64:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kULong64:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat16:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble32:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kLegacyChar:
            // A simple type, or an array of a simple type.
            thisobj = address + offset;
            break;

         default:
            // Everything else is a pointer to something.
            thisobj = *((char**) (address + offset));
            break;
         }
   }
   return thisobj;
}

////////////////////////////////////////////////////////////////////////////////
/// Reminder of the meaning of fMultiplicity:
///  - -1: Only one or 0 element per entry but contains variable length array!
///  -  0: Only one element per entry, no variable length array
///  -  1: loop over the elements of a variable length array
///  -  2: loop over elements of fixed length array (nData is the same for all entry)

Int_t TFormLeafInfo::GetMultiplicity()
{
   // Currently only TFormLeafInfoCast uses this field.
   return fMultiplicity;
}

Int_t TFormLeafInfo::GetNdata(TLeaf* leaf)
{
   // Get the number of element in the entry.

   GetCounterValue(leaf);
   GetValue(leaf);
   return GetNdata();
}

////////////////////////////////////////////////////////////////////////////////
/// Get the number of element in the entry.

Int_t TFormLeafInfo::GetNdata()
{
   if (fNext) return fNext->GetNdata();
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if any of underlying data has a array size counter

Bool_t TFormLeafInfo::HasCounter() const
{
   Bool_t result = kFALSE;
   if (fNext) result = fNext->HasCounter();
   return fCounter!=0 || result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the underlying data is a string

Bool_t TFormLeafInfo::IsString() const
{
   if (fNext) return fNext->IsString();
   if (!fElement) return kFALSE;

   switch (fElement->GetNewType()) {
      // basic types
      case kChar_t:
         // This is new in ROOT 3.02/05
         return kFALSE;
      case TStreamerInfo::kOffsetL + kChar_t:
         // This is new in ROOT 3.02/05
         return kTRUE;
      case TStreamerInfo::kCharStar:
         return kTRUE;
      default:
         return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the underlying data is an integral value

Bool_t TFormLeafInfo::IsInteger() const
{
   if (fNext) return fNext->IsInteger();
   if (!fElement) return kFALSE;

   Int_t atype = fElement->GetNewType();
   if (TStreamerInfo::kOffsetL < atype &&
       atype < TStreamerInfo::kOffsetP ) {
      atype -= TStreamerInfo::kOffsetL;
   } else if (TStreamerInfo::kOffsetP < atype &&
              atype < TStreamerInfo::kObject) {
      atype -= TStreamerInfo::kOffsetP;
   }

   switch (atype) {
      // basic types
      case TStreamerInfo::kLegacyChar:
      case TStreamerInfo::kBool:
      case TStreamerInfo::kChar:
      case TStreamerInfo::kUChar:
      case TStreamerInfo::kShort:
      case TStreamerInfo::kUShort:
      case TStreamerInfo::kInt:
      case TStreamerInfo::kUInt:
      case TStreamerInfo::kLong:
      case TStreamerInfo::kULong:
      case TStreamerInfo::kLong64:
      case TStreamerInfo::kULong64:
         return kTRUE;
      case TStreamerInfo::kCharStar:
         return kTRUE; // For consistency with the leaf list method and proper axis setting
      case TStreamerInfo::kFloat:
      case TStreamerInfo::kFloat16:
      case TStreamerInfo::kDouble:
      case TStreamerInfo::kDouble32:
         return kFALSE;
      default:
         return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Method for multiple variable dimensions.

Int_t TFormLeafInfo::GetPrimaryIndex()
{
   if (fNext) return fNext->GetPrimaryIndex();
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the index of the dimension which varies
/// for each elements of an enclosing array (typically a TClonesArray)

Int_t TFormLeafInfo::GetVarDim()
{
   if (fNext) return fNext->GetVarDim();
   else return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the virtual index (for this expression) of the dimension which varies
/// for each elements of an enclosing array (typically a TClonesArray)

Int_t TFormLeafInfo::GetVirtVarDim()
{
   if (fNext) return fNext->GetVirtVarDim();
   else return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// For the current entry, and the value 'index' for the main array,
/// return the size of the secondary variable dimension of the 'array'.

Int_t TFormLeafInfo::GetSize(Int_t index)
{
   if (fNext) return fNext->GetSize(index);
   else return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Total all the elements that are available for the current entry
/// for the secondary variable dimension.

Int_t TFormLeafInfo::GetSumOfSizes()
{
   if (fNext) return fNext->GetSumOfSizes();
   else return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Load the current array sizes

void TFormLeafInfo::LoadSizes(TBranch* branch)
{
   if (fNext) fNext->LoadSizes(branch);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the primary index value

void TFormLeafInfo::SetPrimaryIndex(Int_t index)
{
   if (fNext) fNext->SetPrimaryIndex(index);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the primary index value

void TFormLeafInfo::SetSecondaryIndex(Int_t index)
{
   if (fNext) fNext->SetSecondaryIndex(index);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the current size of the arrays

void TFormLeafInfo::SetSize(Int_t index, Int_t val)
{
   if (fNext) fNext->SetSize(index, val);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the current sizes of the arrays

void TFormLeafInfo::UpdateSizes(TArrayI *garr)
{
   if (fNext) fNext->UpdateSizes(garr);
}

////////////////////////////////////////////////////////////////////////////////
/// We reloading all cached information in case the underlying class
/// information has changed (for example when changing from the 'emulated'
/// class to the real class.

Bool_t TFormLeafInfo::Update()
{
   if (fClass) {
      TClass * new_class = TClass::GetClass(fClassName);
      if (new_class==fClass) {
         if (fNext) fNext->Update();
         if (fCounter) fCounter->Update();
         return kFALSE;
      }
      fClass = new_class;
   }
   if (fElement && fClass) {
      TClass *cl = fClass;
      // We have to drill down the element name within the class.
      Int_t offset,i;
      TStreamerElement* element;
      char * current;
      Int_t nchname = fElementName.Length();
      char * work = new char[nchname+2];
      for (i=0, current = &(work[0]), fOffset=0; i<nchname+1;i++ ) {
         if (i==nchname || fElementName[i]=='.') {
            // A delimiter happened let's see if what we have seen
            // so far does point to a data member.
            *current = '\0';
            element = ((TStreamerInfo*)cl->GetStreamerInfo())->GetStreamerElement(work,offset);
            if (element) {
               Int_t type = element->GetNewType();
               if (type<60) {
                  fOffset += offset;
               } else if (type == TStreamerInfo::kBase ||
                          type == TStreamerInfo::kAny ||
                          type == TStreamerInfo::kObject ||
                          type == TStreamerInfo::kTString  ||
                          type == TStreamerInfo::kTNamed  ||
                          type == TStreamerInfo::kTObject ||
                          type == TStreamerInfo::kObjectp ||
                          type == TStreamerInfo::kObjectP ||
                          type == TStreamerInfo::kOffsetL + TStreamerInfo::kObjectp ||
                          type == TStreamerInfo::kOffsetL + TStreamerInfo::kObjectP ||
                          type == TStreamerInfo::kAnyp ||
                          type == TStreamerInfo::kAnyP ||
                          type == TStreamerInfo::kOffsetL + TStreamerInfo::kAnyp ||
                          type == TStreamerInfo::kOffsetL + TStreamerInfo::kAnyP ||
                          type == TStreamerInfo::kOffsetL + TStreamerInfo::kSTLp ||
                          type == TStreamerInfo::kSTL ||
                          type == TStreamerInfo::kSTLp ) {
                  fOffset += offset;
                  cl = element->GetClassPointer();
               }
               fElement = element;
               current = &(work[0]);
            }
         } else {
            if (i<nchname) *current++ = fElementName[i];
         }
      }
      delete [] work;
   }
   if (fNext) fNext->Update();
   if (fCounter) fCounter->Update();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///  Return the size of the underlying array for the current entry in the TTree.

Int_t TFormLeafInfo::GetCounterValue(TLeaf* leaf) {
   if (!fCounter) {
      if (fNext && fNext->HasCounter()) {
         char *where = (char*)GetLocalValuePointer(leaf,0);
         return fNext->ReadCounterValue(where);
      } else return 1;
   }
   return (Int_t)fCounter->GetValue(leaf);
}

////////////////////////////////////////////////////////////////////////////////
///  Return the size of the underlying array for the current entry in the TTree.

Int_t TFormLeafInfo::ReadCounterValue(char* where)
{
   if (!fCounter) {
      if (fNext) {
         char *next = (char*)GetLocalValuePointer(where,0);
         return fNext->ReadCounterValue(next);
      } else return 1;
   }
   return (Int_t)fCounter->ReadValue(where,0);
}

////////////////////////////////////////////////////////////////////////////////
/// returns the address of the value pointed to by the
/// TFormLeafInfo.

void* TFormLeafInfo::GetLocalValuePointer(TLeaf *leaf, Int_t instance)
{
   char *thisobj = 0;
   if (leaf->InheritsFrom(TLeafObject::Class()) ) {
      thisobj = (char*)((TLeafObject*)leaf)->GetObject();
   } else {
      thisobj = GetObjectAddress((TLeafElement*)leaf, instance); // instance might be modified
   }
   if (!thisobj) return 0;
   return GetLocalValuePointer(thisobj, instance);
}

////////////////////////////////////////////////////////////////////////////////
/// returns the address of the value pointed to by the
/// serie of TFormLeafInfo.

void* TFormLeafInfo::GetValuePointer(TLeaf *leaf, Int_t instance)
{
   char *thisobj = (char*)GetLocalValuePointer(leaf,instance);
   if (fNext) return fNext->GetValuePointer(thisobj,instance);
   else return thisobj;
}

////////////////////////////////////////////////////////////////////////////////
/// returns the address of the value pointed to by the
/// TFormLeafInfo.

void* TFormLeafInfo::GetValuePointer(char *thisobj, Int_t instance)
{
   char *where = (char*)GetLocalValuePointer(thisobj,instance);
   if (fNext) return fNext->GetValuePointer(where,instance);
   else return where;
}

////////////////////////////////////////////////////////////////////////////////
/// returns the address of the value pointed to by the
/// TFormLeafInfo.

void* TFormLeafInfo::GetLocalValuePointer(char *thisobj, Int_t instance)
{
   if (fElement==0 || thisobj==0) return thisobj;

   switch (fElement->GetNewType()) {
      // basic types
      case TStreamerInfo::kBool:
      case TStreamerInfo::kChar:
      case TStreamerInfo::kUChar:
      case TStreamerInfo::kShort:
      case TStreamerInfo::kUShort:
      case TStreamerInfo::kInt:
      case TStreamerInfo::kUInt:
      case TStreamerInfo::kLong:
      case TStreamerInfo::kULong:
      case TStreamerInfo::kLong64:
      case TStreamerInfo::kULong64:
      case TStreamerInfo::kFloat:
      case TStreamerInfo::kFloat16:
      case TStreamerInfo::kDouble:
      case TStreamerInfo::kDouble32:
      case TStreamerInfo::kLegacyChar:
      case TStreamerInfo::kCounter:
         return (Int_t*)(thisobj+fOffset);

         // array of basic types  array[8]
      case TStreamerInfo::kOffsetL + TStreamerInfo::kBool:
         {Bool_t *val   = (Bool_t*)(thisobj+fOffset);      return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kChar:
         {Char_t *val   = (Char_t*)(thisobj+fOffset);      return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kShort:
         {Short_t *val   = (Short_t*)(thisobj+fOffset);    return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kInt:
         {Int_t *val     = (Int_t*)(thisobj+fOffset);      return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kLong:
         {Long_t *val    = (Long_t*)(thisobj+fOffset);     return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kLong64:
         {Long64_t *val  = (Long64_t*)(thisobj+fOffset);   return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat:
         {Float_t *val   = (Float_t*)(thisobj+fOffset);    return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat16:
         {Float_t *val   = (Float_t*)(thisobj+fOffset);    return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble:
         {Double_t *val  = (Double_t*)(thisobj+fOffset);   return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble32:
         {Double_t *val  = (Double_t*)(thisobj+fOffset);   return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kUChar:
         {UChar_t *val   = (UChar_t*)(thisobj+fOffset);    return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kUShort:
         {UShort_t *val  = (UShort_t*)(thisobj+fOffset);   return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kUInt:
         {UInt_t *val    = (UInt_t*)(thisobj+fOffset);     return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kULong:
         {ULong_t *val   = (ULong_t*)(thisobj+fOffset);    return &(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kULong64:
         {ULong64_t *val  = (ULong64_t*)(thisobj+fOffset); return &(val[instance]);}

#define GET_ARRAY(TYPE_t)                                         \
         {                                                        \
            Int_t len, sub_instance, index;                       \
            if (fNext) len = fNext->GetArrayLength();             \
            else len = 1;                                         \
            if (len) {                                            \
               index = instance / len;                            \
               sub_instance = instance % len;                     \
            } else {                                              \
               index = instance;                                  \
               sub_instance = 0;                                  \
            }                                                     \
            TYPE_t **val     = (TYPE_t**)(thisobj+fOffset);       \
            return &((val[sub_instance])[index]);                 \
         }

         // pointer to an array of basic types  array[n]
      case TStreamerInfo::kOffsetP + TStreamerInfo::kBool:    GET_ARRAY(Bool_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kChar:    GET_ARRAY(Char_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kShort:   GET_ARRAY(Short_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kInt:     GET_ARRAY(Int_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kLong:    GET_ARRAY(Long_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kLong64:  GET_ARRAY(Long64_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat16:
      case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat:   GET_ARRAY(Float_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble32:
      case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble:  GET_ARRAY(Double_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kUChar:   GET_ARRAY(UChar_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kUShort:  GET_ARRAY(UShort_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kUInt:    GET_ARRAY(UInt_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kULong:   GET_ARRAY(ULong_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kULong64: GET_ARRAY(ULong64_t)

      case TStreamerInfo::kCharStar:
         {char **stringp = (char**)(thisobj+fOffset); return *stringp;}

      case TStreamerInfo::kObjectp:
      case TStreamerInfo::kObjectP:
      case TStreamerInfo::kAnyp:
      case TStreamerInfo::kAnyP:
      case TStreamerInfo::kSTLp:
         {TObject **obj = (TObject**)(thisobj+fOffset);   return *obj; }

      case TStreamerInfo::kObject:
      case TStreamerInfo::kTString:
      case TStreamerInfo::kTNamed:
      case TStreamerInfo::kTObject:
      case TStreamerInfo::kAny:
      case TStreamerInfo::kBase:
      case TStreamerInfo::kSTL:
         {TObject *obj = (TObject*)(thisobj+fOffset);   return obj; }

      case TStreamerInfo::kOffsetL + TStreamerInfo::kTObject:
      case TStreamerInfo::kOffsetL + TStreamerInfo::kSTL:
      case TStreamerInfo::kOffsetL + TStreamerInfo::kAny: {
         char *loc = thisobj+fOffset;

         Int_t len, index;
         //Int_t sub_instance;

         if (fNext) len = fNext->GetArrayLength();
         else len = 1;
         if (len) {
            index = instance / len;
            // sub_instance = instance % len;
         } else {
            index = instance;
            // sub_instance = 0;
         }

         loc += index*fElement->GetClassPointer()->Size();

         TObject *obj = (TObject*)(loc);
         return obj;
      }

      case TStreamerInfo::kOffsetL + TStreamerInfo::kObjectp:
      case TStreamerInfo::kOffsetL + TStreamerInfo::kObjectP:
      case TStreamerInfo::kOffsetL + TStreamerInfo::kAnyp:
      case TStreamerInfo::kOffsetL + TStreamerInfo::kAnyP:
      case TStreamerInfo::kOffsetL + TStreamerInfo::kSTLp:
         {TObject *obj = (TObject*)(thisobj+fOffset);   return obj; }

      case kOther_t:
      default:        return 0;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Return result of a leafobject method.

template <typename T>
T TFormLeafInfo::GetValueImpl(TLeaf *leaf, Int_t instance)
{
   char *thisobj = 0;
   if (leaf->InheritsFrom(TLeafObject::Class()) ) {
      thisobj = (char*)((TLeafObject*)leaf)->GetObject();
   } else {
      thisobj = GetObjectAddress((TLeafElement*)leaf, instance); // instance might be modified
   }
   if (thisobj==0) return 0;
   return ReadTypedValue<T>(thisobj,instance);
}

INSTANTIATE_GETVAL(TFormLeafInfo);
INSTANTIATE_READVAL(TFormLeafInfo);

////////////////////////////////////////////////////////////////////////////////
/// Read the value at the given memory location

template <typename T>
T TFormLeafInfo::ReadValueImpl(char *thisobj, Int_t instance)
{
   if ( !thisobj )  {
      Error("ReadValue","Invalid data address: result will be wrong");
      return 0.0;  // Or should throw exception/print error ?
   }
   if (fNext) {
      char *nextobj = thisobj+fOffset;
      Int_t sub_instance = instance;
      Int_t type = fElement->GetNewType();
      if (type==TStreamerInfo::kOffsetL + TStreamerInfo::kObject ||
          type==TStreamerInfo::kOffsetL + TStreamerInfo::kSTL ||
          type==TStreamerInfo::kOffsetL + TStreamerInfo::kAny) {
         Int_t index;
         Int_t len = fNext->GetArrayLength();
         if (len) {
            index = instance / len;
            sub_instance = instance % len;
         } else {
            index = instance;
            sub_instance = 0;
         }
         nextobj += index*fElement->GetClassPointer()->Size();
      }
      return fNext->ReadTypedValue<T>(nextobj,sub_instance);
   }
   //   return fInfo->ReadValue(thisobj+fOffset,fElement->GetNewType(),instance,1);
   switch (fElement->GetNewType()) {
         // basic types
      case TStreamerInfo::kBool:       return (T)(*(Bool_t*)(thisobj+fOffset));
      case TStreamerInfo::kChar:       return (T)(*(Char_t*)(thisobj+fOffset));
      case TStreamerInfo::kUChar:      return (T)(*(UChar_t*)(thisobj+fOffset));
      case TStreamerInfo::kShort:      return (T)(*(Short_t*)(thisobj+fOffset));
      case TStreamerInfo::kUShort:     return (T)(*(UShort_t*)(thisobj+fOffset));
      case TStreamerInfo::kInt:        return (T)(*(Int_t*)(thisobj+fOffset));
      case TStreamerInfo::kUInt:       return (T)(*(UInt_t*)(thisobj+fOffset));
      case TStreamerInfo::kLong:       return (T)(*(Long_t*)(thisobj+fOffset));
      case TStreamerInfo::kULong:      return (T)(*(ULong_t*)(thisobj+fOffset));
      case TStreamerInfo::kLong64:     return (T)(*(Long64_t*)(thisobj+fOffset));
      case TStreamerInfo::kULong64:    return (T)(*(Long64_t*)(thisobj+fOffset)); //cannot cast to ULong64_t with VC++6
      case TStreamerInfo::kFloat:      return (T)(*(Float_t*)(thisobj+fOffset));
      case TStreamerInfo::kFloat16:    return (T)(*(Float_t*)(thisobj+fOffset));
      case TStreamerInfo::kDouble:     return (T)(*(Double_t*)(thisobj+fOffset));
      case TStreamerInfo::kDouble32:   return (T)(*(Double_t*)(thisobj+fOffset));
      case TStreamerInfo::kLegacyChar: return (T)(*(char*)(thisobj+fOffset));
      case TStreamerInfo::kCounter:    return (T)(*(Int_t*)(thisobj+fOffset));

         // array of basic types  array[8]
      case TStreamerInfo::kOffsetL + TStreamerInfo::kBool:
         {Bool_t *val    = (Bool_t*)(thisobj+fOffset);    return T(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kChar:
         {Char_t *val    = (Char_t*)(thisobj+fOffset);    return T(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kShort:
         {Short_t *val   = (Short_t*)(thisobj+fOffset);   return T(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kInt:
         {Int_t *val     = (Int_t*)(thisobj+fOffset);     return T(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kLong:
         {Long_t *val    = (Long_t*)(thisobj+fOffset);    return T(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kLong64:
         {Long64_t *val  = (Long64_t*)(thisobj+fOffset);  return T(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat16:
         {Float_t *val   = (Float_t*)(thisobj+fOffset);   return T(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat:
         {Float_t *val   = (Float_t*)(thisobj+fOffset);   return T(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble:
         {Double_t *val  = (Double_t*)(thisobj+fOffset);  return T(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble32:
         {Double_t *val  = (Double_t*)(thisobj+fOffset);  return T(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kUChar:
         {UChar_t *val   = (UChar_t*)(thisobj+fOffset);   return T(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kUShort:
         {UShort_t *val  = (UShort_t*)(thisobj+fOffset);  return T(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kUInt:
         {UInt_t *val    = (UInt_t*)(thisobj+fOffset);    return T(val[instance]);}
      case TStreamerInfo::kOffsetL + TStreamerInfo::kULong:
         {ULong_t *val   = (ULong_t*)(thisobj+fOffset);   return T(val[instance]);}
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
      case TStreamerInfo::kOffsetL + TStreamerInfo::kULong64:
         {Long64_t *val = (Long64_t*)(thisobj+fOffset);   return T(val[instance]);}
#else
      case TStreamerInfo::kOffsetL + kULong64_t:
         {ULong64_t *val = (ULong64_t*)(thisobj+fOffset); return T(val[instance]);}
#endif

#define READ_ARRAY(TYPE_t)                               \
         {                                               \
            Int_t len, sub_instance, index;              \
            len = GetArrayLength();                      \
            if (len) {                                   \
               index = instance / len;                   \
               sub_instance = instance % len;            \
            } else {                                     \
               index = instance;                         \
               sub_instance = 0;                         \
            }                                            \
            TYPE_t **val =(TYPE_t**)(thisobj+fOffset);   \
            return T((val[sub_instance])[index]); \
         }

         // pointer to an array of basic types  array[n]
      case TStreamerInfo::kOffsetP + TStreamerInfo::kBool:    READ_ARRAY(Bool_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kChar:    READ_ARRAY(Char_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kShort:   READ_ARRAY(Short_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kInt:     READ_ARRAY(Int_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kLong:    READ_ARRAY(Long_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kLong64:  READ_ARRAY(Long64_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat16:
      case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat:   READ_ARRAY(Float_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble32:
      case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble:  READ_ARRAY(Double_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kUChar:   READ_ARRAY(UChar_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kUShort:  READ_ARRAY(UShort_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kUInt:    READ_ARRAY(UInt_t)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kULong:   READ_ARRAY(ULong_t)
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
      case TStreamerInfo::kOffsetP + TStreamerInfo::kULong64: READ_ARRAY(Long64_t)
#else
      case TStreamerInfo::kOffsetP + TStreamerInfo::kULong64: READ_ARRAY(ULong64_t)
#endif

      case kOther_t:
      default:        return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoDirect
/// A small helper class to implement reading a data
/// member on an object stored in a TTree.

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoDirect::TFormLeafInfoDirect(TBranchElement * from) :
   TFormLeafInfo(from->GetInfo()->GetClass(),0,
                 from->GetInfo()->GetElement(from->GetID()))
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this object and its content.

TFormLeafInfo* TFormLeafInfoDirect::DeepCopy() const
{
   return new TFormLeafInfoDirect(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Read the value at the given memory location

Double_t TFormLeafInfoDirect::ReadValue(char * /*where*/, Int_t /*instance*/)
{
   Error("ReadValue","Should not be used in a TFormLeafInfoDirect");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the underlying value.

template <typename T>
T TFormLeafInfoDirect::GetValueImpl(TLeaf *leaf, Int_t instance)
{
   return leaf->GetTypedValue<T>(instance);
}

INSTANTIATE_GETVAL(TFormLeafInfoDirect);

////////////////////////////////////////////////////////////////////////////////
/// Return the address of the underlying value.

void* TFormLeafInfoDirect::GetLocalValuePointer(TLeaf *leaf, Int_t instance)
{
   if (leaf->IsA() != TLeafElement::Class()) {
      return leaf->GetValuePointer();
   } else {
      return GetObjectAddress((TLeafElement*)leaf, instance); // instance might be modified
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Note this should probably never be executed.

void* TFormLeafInfoDirect::GetLocalValuePointer(char *thisobj, Int_t instance)
{
   return TFormLeafInfo::GetLocalValuePointer(thisobj,instance);
}

////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoNumerical
/// A small helper class to implement reading a numerical value inside a collection

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoNumerical::TFormLeafInfoNumerical(EDataType kind) :
   TFormLeafInfo(0,0,0),
   fKind(kind), fIsBool(kFALSE)
{
   fElement = new TStreamerElement("data","in collection", 0, fKind, "");
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a TFormLeafInfo for the numerical type contained in the collection.

TFormLeafInfoNumerical::TFormLeafInfoNumerical(TVirtualCollectionProxy *collection) :
   TFormLeafInfo(0,0,0),
   fKind(kNoType_t), fIsBool(kFALSE)
{
   if (collection) {
      fKind = (EDataType)collection->GetType();
      if (fKind == TStreamerInfo::kOffsetL + TStreamerInfo::kChar) {
         // Could be a bool
         if (strcmp( collection->GetCollectionClass()->GetName(), "vector<bool>") == 0
             || strncmp( collection->GetCollectionClass()->GetName(), "bitset<", strlen("bitset<") ) ==0 ) {
            fIsBool = kTRUE;
            fKind = (EDataType)18;
         }
      }
   }
   fElement = new TStreamerElement("data","in collection", 0, fKind, "");
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoNumerical::TFormLeafInfoNumerical(const TFormLeafInfoNumerical& orig) :
   TFormLeafInfo(orig),
   fKind(orig.fKind), fIsBool(kFALSE)
{
   fElement = new TStreamerElement("data","in collection", 0, fKind, "");
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe swap.

void TFormLeafInfoNumerical::Swap(TFormLeafInfoNumerical& other)
{
   TFormLeafInfo::Swap(other);
   std::swap(fKind,other.fKind);
   std::swap(fIsBool,other.fIsBool);
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe assignment operator.

TFormLeafInfoNumerical &TFormLeafInfoNumerical::operator=(const TFormLeafInfoNumerical &other)
{
   TFormLeafInfoNumerical tmp(other);
   Swap(tmp);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the object and all its content.

TFormLeafInfo* TFormLeafInfoNumerical::DeepCopy() const
{
   return new TFormLeafInfoNumerical(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TFormLeafInfoNumerical::~TFormLeafInfoNumerical()
{
   delete fElement;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the underlying data is a string

Bool_t TFormLeafInfoNumerical::IsString() const
{
   if (fIsBool) return kFALSE;
   return TFormLeafInfo::IsString();
}

////////////////////////////////////////////////////////////////////////////////
/// We reloading all cached information in case the underlying class
/// information has changed (for example when changing from the 'emulated'
/// class to the real class.

Bool_t TFormLeafInfoNumerical::Update()
{
   //R__ASSERT(fNext==0);

   if (fCounter) return fCounter->Update();
   return kFALSE;
}

namespace {
   TStreamerElement *R__GetFakeClonesElem() {
      static TStreamerElement gFakeClonesElem("begin","fake",0,
                                              TStreamerInfo::kAny,
                                              "TClonesArray");
      return &gFakeClonesElem;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoClones
/// A small helper class to implement reading a data member
/// on a TClonesArray object stored in a TTree.

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoClones::TFormLeafInfoClones(TClass* classptr, Long_t offset) :
   TFormLeafInfo(classptr,offset,R__GetFakeClonesElem()),fTop(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoClones::TFormLeafInfoClones(TClass* classptr, Long_t offset,
                                         Bool_t top) :
   TFormLeafInfo(classptr,offset,R__GetFakeClonesElem()),fTop(top)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoClones::TFormLeafInfoClones(TClass* classptr, Long_t offset,
                                         TStreamerElement* element,
                                         Bool_t top) :
   TFormLeafInfo(classptr,offset,element),fTop(top)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Deep Copy constructor.

TFormLeafInfoClones::TFormLeafInfoClones(const TFormLeafInfoClones &orig) :
   TFormLeafInfo(orig), fTop(orig.fTop)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe swap.

void TFormLeafInfoClones::Swap(TFormLeafInfoClones &other)
{
   TFormLeafInfo::Swap(other);
   std::swap(fTop,other.fTop);
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe assignement operator

TFormLeafInfoClones &TFormLeafInfoClones::operator=(const TFormLeafInfoClones &orig)
{
   TFormLeafInfoClones tmp(orig);
   Swap(tmp);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current size of the the TClonesArray

Int_t TFormLeafInfoClones::GetCounterValue(TLeaf* leaf)
{
   if (!fCounter) {
      TClass *clonesClass = TClonesArray::Class();
      Int_t c_offset = 0;
      TStreamerElement *counter = ((TStreamerInfo*)clonesClass->GetStreamerInfo())->GetStreamerElement("fLast",c_offset);
      fCounter = new TFormLeafInfo(clonesClass,c_offset,counter);
   }
   return (Int_t)fCounter->ReadValue((char*)GetLocalValuePointer(leaf)) + 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current size of the the TClonesArray

Int_t TFormLeafInfoClones::ReadCounterValue(char* where)
{
   if (!fCounter) {
      TClass *clonesClass = TClonesArray::Class();
      Int_t c_offset = 0;
      TStreamerElement *counter = ((TStreamerInfo*)clonesClass->GetStreamerInfo())->GetStreamerElement("fLast",c_offset);
      fCounter = new TFormLeafInfo(clonesClass,c_offset,counter);
   }
   return (Int_t)fCounter->ReadValue(where) + 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the value of the underlying data member inside the
/// clones array.

template <typename T>
T TFormLeafInfoClones::ReadValueImpl(char *where, Int_t instance)
{
   if (fNext==0) return 0;
   Int_t len,index,sub_instance;
   len = fNext->GetArrayLength();
   if (len) {
      index = instance / len;
      sub_instance = instance % len;
   } else {
      index = instance;
      sub_instance = 0;
   }
   TClonesArray * clones = (TClonesArray*)where;
   if (!clones) return 0;
   // Note we take advantage of having only one physically variable
   // dimension:
   char * obj = (char*)clones->UncheckedAt(index);
   return fNext->ReadTypedValue<T>(obj,sub_instance);
}

INSTANTIATE_GETVAL(TFormLeafInfoClones);
INSTANTIATE_READVAL(TFormLeafInfoClones);

////////////////////////////////////////////////////////////////////////////////
/// Return the pointer to the clonesArray

void* TFormLeafInfoClones::GetLocalValuePointer(TLeaf *leaf, Int_t /*instance*/)
{
   TClonesArray * clones;
   if (fTop) {
      if (leaf->InheritsFrom(TLeafObject::Class()) ) {
         clones = (TClonesArray*)((TLeafObject*)leaf)->GetObject();
      } else {
         clones = (TClonesArray*)((TBranchElement*)leaf->GetBranch())->GetObject();
      }
   } else {
      clones = (TClonesArray*)TFormLeafInfo::GetLocalValuePointer(leaf);
   }
   return clones;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the address of the underlying current value

void* TFormLeafInfoClones::GetLocalValuePointer(char *where, Int_t instance)
{
   return TFormLeafInfo::GetLocalValuePointer(where,instance);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the value of the underlying data member inside the
/// clones array.

template <typename T>
T TFormLeafInfoClones::GetValueImpl(TLeaf *leaf, Int_t instance)
{
   if (fNext==0) return 0;
   Int_t len,index,sub_instance;
   len = (fNext->fElement==0)? 0 : fNext->GetArrayLength();
   Int_t primary = fNext->GetPrimaryIndex();
   if (len) {
      index = instance / len;
      sub_instance = instance % len;
   } else if (primary>=0) {
      index = primary;
      sub_instance = instance;
   } else {
      index = instance;
      sub_instance = 0;
   }
   TClonesArray *clones = (TClonesArray*)GetLocalValuePointer(leaf);
   if (clones==0) return 0;

   // Note we take advantage of having only one physically variable
   // dimension:
   char * obj = (char*)clones->UncheckedAt(index);
   return fNext->ReadTypedValue<T>(obj,sub_instance);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the pointer to the clonesArray

void * TFormLeafInfoClones::GetValuePointer(TLeaf *leaf, Int_t instance)
{
   TClonesArray * clones = (TClonesArray*)GetLocalValuePointer(leaf);
   if (fNext && clones) {
      // Same as in TFormLeafInfoClones::GetValue
      Int_t len,index,sub_instance;
      len = (fNext->fElement==0)? 0 : fNext->GetArrayLength();
      if (len) {
         index = instance / len;
         sub_instance = instance % len;
      } else {
         index = instance;
         sub_instance = 0;
      }
      return fNext->GetValuePointer((char*)clones->UncheckedAt(index),
                                    sub_instance);
   }
   return clones;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the pointer to the clonesArray

void * TFormLeafInfoClones::GetValuePointer(char *where, Int_t instance)
{
   TClonesArray * clones = (TClonesArray*) where;
   if (fNext) {
      // Same as in TFormLeafInfoClones::GetValue
      Int_t len,index,sub_instance;
      len = (fNext->fElement==0)? 0 : fNext->GetArrayLength();
      if (len) {
         index = instance / len;
         sub_instance = instance % len;
      } else {
         index = instance;
         sub_instance = 0;
      }
      return fNext->GetValuePointer((char*)clones->UncheckedAt(index),
                                    sub_instance);
   }
   return clones;
}

////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoCollectionObject
/// A small helper class to implement reading a data member
/// on a TClonesArray object stored in a TTree.

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoCollectionObject::TFormLeafInfoCollectionObject(TClass* classptr, Bool_t top) :
   TFormLeafInfo(classptr,0,R__GetFakeClonesElem()),fTop(top)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoCollectionObject::TFormLeafInfoCollectionObject(const TFormLeafInfoCollectionObject &orig) :
   TFormLeafInfo(orig),fTop(orig.fTop)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe swap.

void TFormLeafInfoCollectionObject::Swap(TFormLeafInfoCollectionObject &other)
{
   TFormLeafInfo::Swap(other);
   std::swap(fTop,other.fTop);
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe assignement operator

TFormLeafInfoCollectionObject &TFormLeafInfoCollectionObject::operator=(const TFormLeafInfoCollectionObject &orig)
{
   TFormLeafInfoCollectionObject tmp(orig);
   Swap(tmp);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current size of the the TClonesArray

Int_t TFormLeafInfoCollectionObject::GetCounterValue(TLeaf* /* leaf */)
{
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the value of the underlying data member inside the
/// clones array.

Double_t TFormLeafInfoCollectionObject::ReadValue(char * /* where */, Int_t /* instance */)
{
   R__ASSERT(0);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the pointer to the clonesArray

void* TFormLeafInfoCollectionObject::GetLocalValuePointer(TLeaf *leaf, Int_t /*instance*/)
{
   void* collection;
   if (fTop) {
      if (leaf->InheritsFrom(TLeafObject::Class()) ) {
         collection = ((TLeafObject*)leaf)->GetObject();
      } else {
         collection = ((TBranchElement*)leaf->GetBranch())->GetObject();
      }
   } else {
      collection = TFormLeafInfo::GetLocalValuePointer(leaf);
   }
   return collection;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the address of the underlying current value

void* TFormLeafInfoCollectionObject::GetLocalValuePointer(char *where, Int_t instance)
{
   return TFormLeafInfo::GetLocalValuePointer(where,instance);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the value of the underlying data member inside the
/// clones array.

template <typename T>
T TFormLeafInfoCollectionObject::GetValueImpl(TLeaf *leaf, Int_t instance)
{
   char * obj = (char*)GetLocalValuePointer(leaf);

   if (fNext==0) return 0;
   return fNext->ReadTypedValue<T>(obj,instance);
}

INSTANTIATE_GETVAL(TFormLeafInfoCollectionObject);

////////////////////////////////////////////////////////////////////////////////
/// Return the pointer to the clonesArray

void * TFormLeafInfoCollectionObject::GetValuePointer(TLeaf *leaf, Int_t instance)
{
   void *collection = GetLocalValuePointer(leaf);
   if (fNext) {
      return fNext->GetValuePointer((char*)collection,instance);
   }
   return collection;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the pointer to the clonesArray

void * TFormLeafInfoCollectionObject::GetValuePointer(char *where, Int_t instance)
{
   if (fNext) {
      return fNext->GetValuePointer(where,instance);
   }
   return where;
}

////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoCollection
/// A small helper class to implement reading a data
/// member on a generic collection object stored in a TTree.

////////////////////////////////////////////////////////////////////////////////
/// Cosntructor.

TFormLeafInfoCollection::TFormLeafInfoCollection(TClass* classptr,
                                                 Long_t offset,
                                                 TStreamerElement* element,
                                                 Bool_t top) :
   TFormLeafInfo(classptr,offset,element),
   fTop(top),
   fCollClass( 0),
   fCollProxy( 0),
   fLocalElement( 0)
{
   if (element) {
      fCollClass = element->GetClass();
   } else if (classptr) {
      fCollClass = classptr;
   }
   if (fCollClass
       && fCollClass!=TClonesArray::Class()
       && fCollClass->GetCollectionProxy()) {

      fCollProxy = fCollClass->GetCollectionProxy()->Generate();
      fCollClassName = fCollClass->GetName();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoCollection::TFormLeafInfoCollection(TClass* motherclassptr,
                                                 Long_t offset,
                                                 TClass* elementclassptr,
                                                 Bool_t top) :
   TFormLeafInfo(motherclassptr,offset,
                 new TStreamerElement("collection","in class",
                                      0,
                                      TStreamerInfo::kAny,
                                      elementclassptr
                                      ? elementclassptr->GetName()
                                      : ( motherclassptr
                                          ? motherclassptr->GetName()
                                          : "Unknwon")
                                      ) ),
   fTop(top),
   fCollClass( 0),
   fCollProxy( 0) ,
   fLocalElement( fElement )
{
   if (elementclassptr) {
      fCollClass = elementclassptr;
   } else if (motherclassptr) {
      fCollClass = motherclassptr;
   }
   if (fCollClass
       && fCollClass!=TClonesArray::Class()
       && fCollClass->GetCollectionProxy())
      {
         fCollProxy = fCollClass->GetCollectionProxy()->Generate();
         fCollClassName = fCollClass->GetName();
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoCollection::TFormLeafInfoCollection() :
   TFormLeafInfo(),
   fTop(kFALSE),
   fCollClass( 0),
   fCollProxy( 0),
   fLocalElement( 0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoCollection::TFormLeafInfoCollection(const TFormLeafInfoCollection& orig) :
   TFormLeafInfo(orig),
   fTop( orig.fTop),
   fCollClass( orig.fCollClass ),
   fCollClassName( orig.fCollClassName ),
   fCollProxy( orig.fCollProxy ? orig.fCollProxy->Generate() : 0 ),
   fLocalElement( 0 ) // humm why not initialize it?
{
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe swap.

void TFormLeafInfoCollection::Swap(TFormLeafInfoCollection &other)
{
   TFormLeafInfo::Swap(other);
   std::swap(fTop,other.fTop);
   std::swap(fClass,other.fClass);
   std::swap(fCollClassName,other.fCollClassName);
   std::swap(fCollProxy,other.fCollProxy);
   std::swap(fLocalElement,other.fLocalElement);
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe assignment operator.

TFormLeafInfoCollection &TFormLeafInfoCollection::operator=(const TFormLeafInfoCollection &other)
{
   TFormLeafInfoCollection tmp(other);
   Swap(tmp);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TFormLeafInfoCollection::~TFormLeafInfoCollection()
{
   delete fCollProxy;
   delete fLocalElement;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy of the object and its content.

TFormLeafInfo* TFormLeafInfoCollection::DeepCopy() const
{
   return new TFormLeafInfoCollection(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// We reloading all cached information in case the underlying class
/// information has changed (for example when changing from the 'emulated'
/// class to the real class.

Bool_t TFormLeafInfoCollection::Update()
{
   Bool_t changed = kFALSE;
   TClass * new_class = TClass::GetClass(fCollClassName);
   if (new_class!=fCollClass) {
      delete fCollProxy; fCollProxy = 0;
      fCollClass = new_class;
      if (fCollClass && fCollClass->GetCollectionProxy()) {
         fCollProxy = fCollClass->GetCollectionProxy()->Generate();
      }
      changed = kTRUE;
   }
   return changed || TFormLeafInfo::Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the underlying data has a array size counter

Bool_t TFormLeafInfoCollection::HasCounter() const
{
   return fCounter!=0 || fCollProxy!=0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current size of the the TClonesArray

Int_t TFormLeafInfoCollection::GetCounterValue(TLeaf* leaf)
{
   void *ptr = GetLocalValuePointer(leaf);

   if (fCounter) { return (Int_t)fCounter->ReadValue((char*)ptr); }

   R__ASSERT(fCollProxy);
   if (ptr==0) return 0;
   TVirtualCollectionProxy::TPushPop helper(fCollProxy, ptr);
   return (Int_t)fCollProxy->Size();
}

////////////////////////////////////////////////////////////////////////////////
///  Return the size of the underlying array for the current entry in the TTree.

Int_t TFormLeafInfoCollection::ReadCounterValue(char* where)
{
   if (fCounter) { return (Int_t)fCounter->ReadValue(where); }
   R__ASSERT(fCollProxy);
   if (where==0) return 0;
   void *ptr = GetLocalValuePointer(where,0);
   TVirtualCollectionProxy::TPushPop helper(fCollProxy, ptr);
   return (Int_t)fCollProxy->Size();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current size of the the TClonesArray

Int_t TFormLeafInfoCollection::GetCounterValue(TLeaf* leaf, Int_t instance)
{
   void *ptr = GetLocalValuePointer(leaf,instance);
   if (fCounter) {
      return (Int_t)fCounter->ReadValue((char*)ptr);
   }
   R__ASSERT(fCollProxy);
   if (ptr==0) return 0;
   TVirtualCollectionProxy::TPushPop helper(fCollProxy, ptr);
   return (Int_t)fCollProxy->Size();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the value of the underlying data member inside the
/// clones array.

template <typename T>
T TFormLeafInfoCollection::ReadValueImpl(char *where, Int_t instance)
{
   if (fNext==0) return 0;
   UInt_t len,index,sub_instance;
   len = (fNext->fElement==0)? 0 : fNext->GetArrayLength();
   Int_t primary = fNext->GetPrimaryIndex();
   if (len) {
      index = instance / len;
      sub_instance = instance % len;
   } else if (primary>=0) {
      index = primary;
      sub_instance = instance;
   } else {
      index = instance;
      sub_instance = 0;
   }

   R__ASSERT(fCollProxy);
   void *ptr = GetLocalValuePointer(where,instance);
   TVirtualCollectionProxy::TPushPop helper(fCollProxy, ptr);

   // Note we take advantage of having only one physically variable
   // dimension:

   char * obj = (char*)fCollProxy->At(index);
   if (fCollProxy->HasPointers()) obj = *(char**)obj;
   return fNext->ReadTypedValue<T>(obj,sub_instance);
}

INSTANTIATE_GETVAL(TFormLeafInfoCollection);
INSTANTIATE_READVAL(TFormLeafInfoCollection);

////////////////////////////////////////////////////////////////////////////////
/// Return the pointer to the clonesArray

void* TFormLeafInfoCollection::GetLocalValuePointer(TLeaf *leaf, Int_t /*instance*/)
{
   void *collection;
   if (fTop) {
      if (leaf->InheritsFrom(TLeafObject::Class()) ) {
         collection = ((TLeafObject*)leaf)->GetObject();
      } else {
         collection = ((TBranchElement*)leaf->GetBranch())->GetObject();
      }
   } else {
      collection = TFormLeafInfo::GetLocalValuePointer(leaf);
   }
   return collection;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the address of the local value

void* TFormLeafInfoCollection::GetLocalValuePointer(char *where, Int_t instance)
{
   return TFormLeafInfo::GetLocalValuePointer(where,instance);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the value of the underlying data member inside the
/// clones array.

template <typename T>
T TFormLeafInfoCollection::GetValueImpl(TLeaf *leaf, Int_t instance)
{
   if (fNext==0) return 0;
   Int_t len,index,sub_instance;
   len = (fNext->fElement==0)? 0 : fNext->GetArrayLength();
   Int_t primary = fNext->GetPrimaryIndex();
   if (len) {
      index = instance / len;
      sub_instance = instance % len;
   } else if (primary>=0) {
      index = primary;
      sub_instance = instance;
   } else {
      index = instance;
      sub_instance = 0;
   }

   R__ASSERT(fCollProxy);
   void *coll = GetLocalValuePointer(leaf);
   TVirtualCollectionProxy::TPushPop helper(fCollProxy,coll);

   // Note we take advantage of having only one physically variable
   // dimension:
   char * obj = (char*)fCollProxy->At(index);
   if (obj==0) return 0;
   if (fCollProxy->HasPointers()) obj = *(char**)obj;
   if (obj==0) return 0;
   return fNext->ReadTypedValue<T>(obj,sub_instance);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the pointer to the clonesArray

void * TFormLeafInfoCollection::GetValuePointer(TLeaf *leaf, Int_t instance)
{
   R__ASSERT(fCollProxy);

   void *collection = GetLocalValuePointer(leaf);

   if (fNext) {
      // Same as in TFormLeafInfoClones::GetValue
      Int_t len,index,sub_instance;
      if (fNext->fElement &&
         (fNext->fNext || !fNext->IsString()) ) {
         len = fNext->GetArrayLength();
      } else {
         len = 0;
      }
      if (len) {
         index = instance / len;
         sub_instance = instance % len;
      } else {
         index = instance;
         sub_instance = 0;
      }
      TVirtualCollectionProxy::TPushPop helper(fCollProxy,collection);
      char * obj = (char*)fCollProxy->At(index);
      if (fCollProxy->HasPointers()) obj = *(char**)obj;
      return fNext->GetValuePointer(obj,sub_instance);
   }
   return collection;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the pointer to the clonesArray

void * TFormLeafInfoCollection::GetValuePointer(char *where, Int_t instance)
{
   R__ASSERT(fCollProxy);

   void *collection = where;

   if (fNext) {
      // Same as in TFormLeafInfoClones::GetValue
      Int_t len,index,sub_instance;
      len = (fNext->fElement==0)? 0 : fNext->GetArrayLength();
      if (len) {
         index = instance / len;
         sub_instance = instance % len;
      } else {
         index = instance;
         sub_instance = 0;
      }
      TVirtualCollectionProxy::TPushPop helper(fCollProxy,collection);
      char * obj = (char*)fCollProxy->At(index);
      if (fCollProxy->HasPointers()) obj = *(char**)obj;
      return fNext->GetValuePointer(obj,sub_instance);
   }
   return collection;
}

////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoCollectionSize
/// Used to return the size of a collection

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoCollectionSize::TFormLeafInfoCollectionSize(TClass* classptr) :
   TFormLeafInfo(), fCollClass(classptr), fCollProxy(0)
{
   if (fCollClass
       && fCollClass!=TClonesArray::Class()
       && fCollClass->GetCollectionProxy()) {

      fCollProxy = fCollClass->GetCollectionProxy()->Generate();
      fCollClassName = fCollClass->GetName();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoCollectionSize::TFormLeafInfoCollectionSize(
   TClass* classptr,Long_t offset,TStreamerElement* element) :
   TFormLeafInfo(classptr,offset,element), fCollClass(element->GetClassPointer()), fCollProxy(0)
{
   if (fCollClass
       && fCollClass!=TClonesArray::Class()
       && fCollClass->GetCollectionProxy()) {

      fCollProxy = fCollClass->GetCollectionProxy()->Generate();
      fCollClassName = fCollClass->GetName();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoCollectionSize::TFormLeafInfoCollectionSize() :
   TFormLeafInfo(), fCollClass(0), fCollProxy(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoCollectionSize::TFormLeafInfoCollectionSize(
   const TFormLeafInfoCollectionSize& orig) :  TFormLeafInfo(),
      fCollClass(orig.fCollClass),
      fCollClassName(orig.fCollClassName),
      fCollProxy(orig.fCollProxy?orig.fCollProxy->Generate():0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe assignment operator.

TFormLeafInfoCollectionSize &TFormLeafInfoCollectionSize::operator=(const TFormLeafInfoCollectionSize &other)
{
   TFormLeafInfoCollectionSize tmp(other);
   Swap(tmp);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe swap.

void TFormLeafInfoCollectionSize::Swap(TFormLeafInfoCollectionSize &other)
{
   TFormLeafInfo::Swap(other);
   std::swap(fCollClass,other.fCollClass);
   std::swap(fCollClassName,other.fCollClassName);
   std::swap(fCollProxy,other.fCollProxy);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TFormLeafInfoCollectionSize::~TFormLeafInfoCollectionSize()
{
   delete fCollProxy;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the object and all of its content.

TFormLeafInfo* TFormLeafInfoCollectionSize::DeepCopy() const
{
   return new TFormLeafInfoCollectionSize(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// We reloading all cached information in case the underlying class
/// information has changed (for example when changing from the 'emulated'
/// class to the real class.

Bool_t TFormLeafInfoCollectionSize::Update()
{
   Bool_t changed = kFALSE;
   TClass *new_class = TClass::GetClass(fCollClassName);
   if (new_class!=fCollClass) {
      delete fCollProxy; fCollProxy = 0;
      fCollClass = new_class;
      if (fCollClass && fCollClass->GetCollectionProxy()) {
         fCollProxy = fCollClass->GetCollectionProxy()->Generate();
      }
      changed = kTRUE;
   }
   return changed;
}

////////////////////////////////////////////////////////////////////////////////
/// Not implemented.

void *TFormLeafInfoCollectionSize::GetValuePointer(TLeaf * /* leaf */, Int_t  /* instance */)
{
   Error("GetValuePointer","This should never be called");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Not implemented.

void *TFormLeafInfoCollectionSize::GetValuePointer(char  * /* from */, Int_t  /* instance */)
{
   Error("GetValuePointer","This should never be called");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Not implemented.

void *TFormLeafInfoCollectionSize::GetLocalValuePointer(TLeaf * /* leaf */, Int_t  /* instance */)
{
   Error("GetLocalValuePointer","This should never be called");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Not implemented.

void *TFormLeafInfoCollectionSize::GetLocalValuePointer( char * /* from */, Int_t  /* instance */)
{
   Error("GetLocalValuePointer","This should never be called");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the value of the underlying pointer data member

Double_t  TFormLeafInfoCollectionSize::ReadValue(char *where, Int_t /* instance */)
{
   R__ASSERT(fCollProxy);
   if (where==0) return 0;
   void *ptr = fElement ? TFormLeafInfo::GetLocalValuePointer(where) : where;
   TVirtualCollectionProxy::TPushPop helper(fCollProxy, ptr);
   return (Int_t)fCollProxy->Size();
}

////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoPointer
/// A small helper class to implement reading a data
/// member by following a pointer inside a branch of TTree.

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoPointer::TFormLeafInfoPointer(TClass* classptr,
                                           Long_t offset,
                                           TStreamerElement* element) :
   TFormLeafInfo(classptr,offset,element)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the object and all of its contnet.

TFormLeafInfo* TFormLeafInfoPointer::DeepCopy() const
{
   return new TFormLeafInfoPointer(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the value of the underlying pointer data member

template <typename T>
T  TFormLeafInfoPointer::ReadValueImpl(char *where, Int_t instance)
{
   if (!fNext) return 0;
   char * whereoffset = where+fOffset;
   switch (fElement->GetNewType()) {
      // basic types
      case TStreamerInfo::kObjectp:
      case TStreamerInfo::kObjectP:
      case TStreamerInfo::kAnyp:
      case TStreamerInfo::kAnyP:
      case TStreamerInfo::kSTLp:
      {TObject **obj = (TObject**)(whereoffset);
      return obj && *obj ? fNext->ReadTypedValue<T>((char*)*obj,instance) : 0; }

      case TStreamerInfo::kObject:
      case TStreamerInfo::kTString:
      case TStreamerInfo::kTNamed:
      case TStreamerInfo::kTObject:
      case TStreamerInfo::kAny:
      case TStreamerInfo::kBase:
      case TStreamerInfo::kSTL:
         {
            TObject *obj = (TObject*)(whereoffset);
            return fNext->ReadTypedValue<T>((char*)obj,instance);
         }

      case TStreamerInfo::kOffsetL + TStreamerInfo::kTObject:
      case TStreamerInfo::kOffsetL + TStreamerInfo::kSTL:
      case TStreamerInfo::kOffsetL + TStreamerInfo::kAny:
         {
            Int_t len, index, sub_instance;

            if (fNext) len = fNext->GetArrayLength();
            else len = 1;
            if (len) {
               index = instance / len;
               sub_instance = instance % len;
            } else {
               index = instance;
               sub_instance = 0;
            }

            whereoffset += index*fElement->GetClassPointer()->Size();

            TObject *obj = (TObject*)(whereoffset);
            return fNext->ReadTypedValue<T>((char*)obj,sub_instance);
         }

      case TStreamerInfo::kOffsetL + TStreamerInfo::kObjectp:
      case TStreamerInfo::kOffsetL + TStreamerInfo::kObjectP:
      case TStreamerInfo::kOffsetL + TStreamerInfo::kAnyp:
      case TStreamerInfo::kOffsetL + TStreamerInfo::kAnyP:
      case TStreamerInfo::kOffsetL + TStreamerInfo::kSTLp:
         {
            TObject *obj = (TObject*)(whereoffset);
            return fNext->ReadTypedValue<T>((char*)obj,instance);
         }

      case kOther_t:
      default:        return 0;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Return the value of the underlying pointer data member

template <typename T>
T  TFormLeafInfoPointer::GetValueImpl(TLeaf *leaf, Int_t instance)
{
   if (!fNext) return 0;
   char * where = (char*)GetLocalValuePointer(leaf,instance);
   if (where==0) return 0;
   return fNext->ReadTypedValue<T>(where,instance);
}

INSTANTIATE_GETVAL(TFormLeafInfoPointer);
INSTANTIATE_READVAL(TFormLeafInfoPointer);

////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoMethod
/// Asmall helper class to implement executing a method
/// of an object stored in a TTree

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoMethod::TFormLeafInfoMethod( TClass* classptr,
                                          TMethodCall *method) :
   TFormLeafInfo(classptr,0,0),fMethod(method),
   fResult(0), fCopyFormat(),fDeleteFormat(),fValuePointer(0),fIsByValue(kFALSE)
{
   if (method) {
      fMethodName = method->GetMethodName();
      fParams = method->GetParams();
      TMethodCall::EReturnType r = fMethod->ReturnType();
      if (r == TMethodCall::kOther) {
         const char* rtype = fMethod->GetMethod()->GetReturnTypeName();
         Long_t rprop = fMethod->GetMethod()->Property();
         if (rtype[strlen(rtype)-1]!='*' &&
             rtype[strlen(rtype)-1]!='&' &&
             !(rprop & (kIsPointer|kIsReference)) ) {
            fCopyFormat = "new ";
            fCopyFormat += rtype;
            fCopyFormat += "(*(";
            fCopyFormat += rtype;
            fCopyFormat += "*)0x%lx)";

            fDeleteFormat  = "delete (";
            fDeleteFormat += rtype;
            fDeleteFormat += "*)0x%lx";

            fIsByValue = kTRUE;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoMethod::TFormLeafInfoMethod(const TFormLeafInfoMethod& orig)
   : TFormLeafInfo(orig)
{
   fMethodName = orig.fMethodName;
   fParams = orig.fParams ;
   fResult = orig.fResult;
   if (orig.fMethod) {
      fMethod = new TMethodCall();
      fMethod->Init(orig.fMethod->GetMethod());
   } else {
      fMethod = 0;
   }
   fCopyFormat = orig.fCopyFormat;
   fDeleteFormat = orig.fDeleteFormat;
   fValuePointer = 0;
   fIsByValue = orig.fIsByValue;
}


////////////////////////////////////////////////////////////////////////////////
/// Exception safe swap.

void TFormLeafInfoMethod::Swap(TFormLeafInfoMethod &other)
{
   TFormLeafInfo::Swap(other);
   std::swap(fMethod,other.fMethod);
   std::swap(fMethodName,other.fMethodName);
   std::swap(fParams,other.fParams);
   std::swap(fResult,other.fResult);
   std::swap(fCopyFormat,other.fCopyFormat);
   std::swap(fDeleteFormat,other.fDeleteFormat);
   std::swap(fValuePointer,other.fValuePointer);
   std::swap(fIsByValue,other.fIsByValue);
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe assignment operator.

TFormLeafInfoMethod &TFormLeafInfoMethod::operator=(const TFormLeafInfoMethod &other)
{
   TFormLeafInfoMethod tmp(other);
   Swap(tmp);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TFormLeafInfoMethod::~TFormLeafInfoMethod()
{
   if (fValuePointer) {
      gInterpreter->Calc(Form(fDeleteFormat.Data(),fValuePointer));
   }
   delete fMethod;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the object and all its content.

TFormLeafInfo* TFormLeafInfoMethod::DeepCopy() const
{
   return new TFormLeafInfoMethod(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the TClass corresponding to the return type of the function
/// if it is an object type or if the return type is a reference (&) then
/// return the TClass corresponding to the referencee.

TClass *TFormLeafInfoMethod::ReturnTClass(TMethodCall *mc)
{
   if (!mc || !mc->GetMethod())
      return nullptr;

   std::string return_type;

   if (0 == strcmp(mc->GetMethod()->GetReturnTypeName(), "void"))
      return nullptr;

   R__WRITE_LOCKGUARD(ROOT::gCoreMutex);

   {
      TInterpreter::SuspendAutoLoadingRAII autoloadOff(gInterpreter);
      TClassEdit::GetNormalizedName(return_type, mc->GetMethod()->GetReturnTypeName());
   }
   // Beyhond this point we no longer 'need' the lock.
   // How TClass::GetClass will take at least the read lock to search
   // So keeping it just a little longer is likely to be faster
   // than releasing and retaking it.

   return_type = gInterpreter->TypeName(return_type.c_str());

   if (return_type == "void")
      return nullptr;

   return TClass::GetClass(return_type.c_str());
}

////////////////////////////////////////////////////////////////////////////////
/// Return the type of the underlying return value

TClass* TFormLeafInfoMethod::GetClass() const
{
   if (fNext) return fNext->GetClass();
   TMethodCall::EReturnType r = fMethod->ReturnType();
   if (r!=TMethodCall::kOther) return 0;

   return ReturnTClass(fMethod);
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the return value is integral.

Bool_t TFormLeafInfoMethod::IsInteger() const
{
   TMethodCall::EReturnType r = fMethod->ReturnType();
   if (r == TMethodCall::kLong) {
      return kTRUE;
   } else return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the return value is a string.

Bool_t TFormLeafInfoMethod::IsString() const
{
   if (fNext) return fNext->IsString();

   TMethodCall::EReturnType r = fMethod->ReturnType();
   return (r==TMethodCall::kString);
}

////////////////////////////////////////////////////////////////////////////////
/// We reloading all cached information in case the underlying class
/// information has changed (for example when changing from the 'emulated'
/// class to the real class.

Bool_t TFormLeafInfoMethod::Update()
{
   if (!TFormLeafInfo::Update()) return kFALSE;
   delete fMethod;
   fMethod = new TMethodCall(fClass, fMethodName, fParams);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// This is implemented here because some compiler want ALL the
/// signature of an overloaded function to be re-implemented.

void *TFormLeafInfoMethod::GetLocalValuePointer( TLeaf *from,
                                                 Int_t instance)
{
   return TFormLeafInfo::GetLocalValuePointer( from, instance);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the address of the lcoal underlying value.

void *TFormLeafInfoMethod::GetLocalValuePointer(char *from,
                                                Int_t /*instance*/)
{
   void *thisobj = from;
   if (!thisobj) return 0;

   TMethodCall::EReturnType r = fMethod->ReturnType();
   fResult = 0;

   if (r == TMethodCall::kLong) {
      Long_t l = 0;
      fMethod->Execute(thisobj, l);
      fResult = (Double_t) l;
      // Get rid of temporary return object.
      gInterpreter->ClearStack();
      return &fResult;

   } else if (r == TMethodCall::kDouble) {
      Double_t d = 0;
      fMethod->Execute(thisobj, d);
      fResult = (Double_t) d;
      // Get rid of temporary return object.
      gInterpreter->ClearStack();
      return &fResult;

   } else if (r == TMethodCall::kString) {
      char *returntext = 0;
      fMethod->Execute(thisobj,&returntext);
      gInterpreter->ClearStack();
      return returntext;

   } else if (r == TMethodCall::kOther) {
      char * char_result = 0;
      if (fIsByValue) {
         if (fValuePointer) {
            gROOT->ProcessLine(Form(fDeleteFormat.Data(),fValuePointer));
            fValuePointer = 0;
         }
      }
      fMethod->Execute(thisobj, &char_result);
      if (fIsByValue) {
         fValuePointer = (char*)gInterpreter->Calc(Form(fCopyFormat.Data(),char_result));
         char_result = (char*)fValuePointer;
      }
      gInterpreter->ClearStack();
      return char_result;

   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute the method on the given address

template <typename T>
T TFormLeafInfoMethod::ReadValueImpl(char *where, Int_t instance)
{
   void *thisobj = where;
   if (!thisobj) return 0;

   TMethodCall::EReturnType r = fMethod->ReturnType();
   T result = 0;

   if (r == TMethodCall::kLong) {
      Long_t l = 0;
      fMethod->Execute(thisobj, l);
      result = (T) l;

   } else if (r == TMethodCall::kDouble) {
      Double_t d = 0;
      fMethod->Execute(thisobj, d);
      result = (T) d;

   } else if (r == TMethodCall::kString) {
      char *returntext = 0;
      fMethod->Execute(thisobj,&returntext);
      result = T((Long_t) returntext);

   } else if (fNext) {
      char * char_result = 0;
      fMethod->Execute(thisobj, &char_result);
      result = fNext->ReadTypedValue<T>(char_result,instance);

   } else fMethod->Execute(thisobj);

   // Get rid of temporary return object.
   gInterpreter->ClearStack();
   return result;
}

INSTANTIATE_READVAL(TFormLeafInfoMethod);

////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoMultiVarDim
/// A helper class to implement reading a
/// data member on a variable size array inside a TClonesArray object stored in
/// a TTree.  This is the version used when the data member is inside a
/// non-split object.

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoMultiVarDim::TFormLeafInfoMultiVarDim( TClass* classptr,
                                                    Long_t offset,
                                                    TStreamerElement* element,
                                                    TFormLeafInfo* parent) :
   TFormLeafInfo(classptr,offset,element),fNsize(0),fCounter2(0),fSumOfSizes(0),
   fDim(0),fVirtDim(-1),fPrimaryIndex(-1),fSecondaryIndex(-1)
{
   if (element && element->InheritsFrom(TStreamerBasicPointer::Class())) {
      TStreamerBasicPointer * elem = (TStreamerBasicPointer*)element;

      Int_t counterOffset = 0;
      TStreamerElement* counter = ((TStreamerInfo*)classptr->GetStreamerInfo())->GetStreamerElement(elem->GetCountName(),counterOffset);
      if (!parent) return;
      fCounter2 = parent->DeepCopy();
      TFormLeafInfo ** next = &(fCounter2->fNext);
      while(*next != 0) next = &( (*next)->fNext);
      *next = new TFormLeafInfo(classptr,counterOffset,counter);

   } else Error("Constructor","Called without a proper TStreamerElement");
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoMultiVarDim::TFormLeafInfoMultiVarDim() :
   TFormLeafInfo(0,0,0),fNsize(0),fCounter2(0),fSumOfSizes(0),
   fDim(0),fVirtDim(-1),fPrimaryIndex(-1),fSecondaryIndex(-1)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoMultiVarDim::TFormLeafInfoMultiVarDim(const TFormLeafInfoMultiVarDim& orig) : TFormLeafInfo(orig)
{
   fNsize = orig.fNsize;
   orig.fSizes.Copy(fSizes);
   fCounter2 = orig.fCounter2?orig.fCounter2->DeepCopy():0;
   fSumOfSizes = orig.fSumOfSizes;
   fDim = orig.fDim;
   fVirtDim = orig.fVirtDim;
   fPrimaryIndex = orig.fPrimaryIndex;
   fSecondaryIndex = orig.fSecondaryIndex;
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe swap.

void TFormLeafInfoMultiVarDim::Swap(TFormLeafInfoMultiVarDim &other)
{
   TFormLeafInfo::Swap(other);
   std::swap(fNsize,other.fNsize);
   std::swap(fSizes,other.fSizes);
   std::swap(fSumOfSizes,other.fSumOfSizes);
   std::swap(fDim,other.fDim);
   std::swap(fVirtDim,other.fVirtDim);
   std::swap(fPrimaryIndex,other.fPrimaryIndex);
   std::swap(fSecondaryIndex,other.fSecondaryIndex);
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe assignment operator.

TFormLeafInfoMultiVarDim &TFormLeafInfoMultiVarDim::operator=(const TFormLeafInfoMultiVarDim &other)
{
   TFormLeafInfoMultiVarDim tmp(other);
   Swap(tmp);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the object and all its content.

TFormLeafInfo* TFormLeafInfoMultiVarDim::DeepCopy() const
{
   return new TFormLeafInfoMultiVarDim(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TFormLeafInfoMultiVarDim:: ~TFormLeafInfoMultiVarDim()
{
   delete fCounter2;
}

/* The proper indexing and unwinding of index is done by prior leafinfo in the chain. */
//virtual Double_t  TFormLeafInfoMultiVarDim::ReadValue(char *where, Int_t instance = 0) {
//   return TFormLeafInfo::ReadValue(where,instance);
//}

////////////////////////////////////////////////////////////////////////////////
/// Load the current array sizes.

void TFormLeafInfoMultiVarDim::LoadSizes(TBranch* branch)
{
   if (fElement) {
      TLeaf *leaf = (TLeaf*)branch->GetListOfLeaves()->At(0);
      if (fCounter) fNsize = (Int_t)fCounter->GetValue(leaf);
      else fNsize = fCounter2->GetCounterValue(leaf);
      if (fNsize > fSizes.GetSize()) fSizes.Set(fNsize);
      fSumOfSizes = 0;
      for (Int_t i=0; i<fNsize; i++) {
         Int_t size = (Int_t)fCounter2->GetValue(leaf,i);
         fSumOfSizes += size;
         fSizes.AddAt( size, i );
      }
      return;
   }
   if (!fCounter2 || !fCounter) return;
   TBranchElement *br = dynamic_cast<TBranchElement*>(branch);
   R__ASSERT(br);
   fNsize = br->GetBranchCount()->GetNdata();
   if (fNsize > fSizes.GetSize()) fSizes.Set(fNsize);
   fSumOfSizes = 0;
   for (Int_t i=0; i<fNsize; i++) {
      Int_t size = (Int_t)fCounter2->GetValue((TLeaf*)br->GetBranchCount2()->GetListOfLeaves()->At(0),i);
      fSumOfSizes += size;
      fSizes.AddAt( size, i );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the index vlaue of the primary index.

Int_t TFormLeafInfoMultiVarDim::GetPrimaryIndex()
{
   return fPrimaryIndex;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the size of the requested sub-array.

Int_t TFormLeafInfoMultiVarDim::GetSize(Int_t index)
{
   if (index >= fSizes.GetSize()) {
      return -1;
   } else {
      return fSizes.At(index);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the current value of the primary index.

void TFormLeafInfoMultiVarDim::SetPrimaryIndex(Int_t index)
{
   fPrimaryIndex = index;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the current value of the primary index.

void TFormLeafInfoMultiVarDim::SetSecondaryIndex(Int_t index)
{
   fSecondaryIndex = index;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the sizes of the sub-array.

void TFormLeafInfoMultiVarDim::SetSize(Int_t index, Int_t val)
{
   fSumOfSizes += (val - fSizes.At(index));
   fSizes.AddAt(val,index);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the total size.

Int_t TFormLeafInfoMultiVarDim::GetSumOfSizes()
{
   return fSumOfSizes;
}

////////////////////////////////////////////////////////////////////////////////

Double_t TFormLeafInfoMultiVarDim::GetValue(TLeaf * /*leaf*/,
                                            Int_t /*instance*/)
{
   /* The proper indexing and unwinding of index need to be done by prior leafinfo in the chain. */
   Error("GetValue","This should never be called");
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the index of the dimension which varies
/// for each elements of an enclosing array (typically a TClonesArray)

Int_t TFormLeafInfoMultiVarDim::GetVarDim()
{
   return fDim;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the virtual index (for this expression) of the dimension which varies
/// for each elements of an enclosing array (typically a TClonesArray)

Int_t TFormLeafInfoMultiVarDim::GetVirtVarDim()
{
   return fVirtDim;
}

////////////////////////////////////////////////////////////////////////////////
/// We reloading all cached information in case the underlying class
/// information has changed (for example when changing from the 'emulated'
/// class to the real class.

Bool_t TFormLeafInfoMultiVarDim::Update()
{
   Bool_t res = TFormLeafInfo::Update();
   if (fCounter2) fCounter2->Update();
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Update the sizes of the arrays.

void TFormLeafInfoMultiVarDim::UpdateSizes(TArrayI *garr)
{
   if (!garr) return;
   if (garr->GetSize()<fNsize) garr->Set(fNsize);
   for (Int_t i=0; i<fNsize; i++) {
      Int_t local = fSizes.At(i);
      Int_t global = garr->At(i);
      if (global==0 || local<global) global = local;
      garr->AddAt(global,i);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoMultiVarDimDirect
/// A small helper class to implement reading
/// a data member on a variable size array inside a TClonesArray object stored
/// in a TTree.  This is the version used for split access

////////////////////////////////////////////////////////////////////////////////
/// Copy the object and all its content.

TFormLeafInfo* TFormLeafInfoMultiVarDimDirect::DeepCopy() const
{
   return new TFormLeafInfoMultiVarDimDirect(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the undersying value.

template <typename T>
T TFormLeafInfoMultiVarDimDirect::GetValueImpl(TLeaf *leaf, Int_t instance)
{
   return ((TLeafElement*)leaf)->GetTypedValueSubArray<T>(fPrimaryIndex,instance);
}

INSTANTIATE_GETVAL(TFormLeafInfoMultiVarDimDirect);

////////////////////////////////////////////////////////////////////////////////
/// Not implemented.

Double_t TFormLeafInfoMultiVarDimDirect::ReadValue(char * /*where*/, Int_t /*instance*/)
{
   Error("ReadValue","This should never be called");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoMultiVarDimCollection
/// A small helper class to implement reading
/// a data member on a variable size array inside a TClonesArray object stored
/// in a TTree.  This is the version used for split access

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoMultiVarDimCollection::TFormLeafInfoMultiVarDimCollection(
   TClass* motherclassptr,
   Long_t offset,
   TClass* elementclassptr,
   TFormLeafInfo *parent) :
   TFormLeafInfoMultiVarDim(motherclassptr,offset,
                 new TStreamerElement("collection","in class",
                                      0,
                                      TStreamerInfo::kAny,
                                      elementclassptr
                                      ? elementclassptr->GetName()
                                      : ( motherclassptr
                                          ? motherclassptr->GetName()
                                          : "Unknwon")
                                          )
                                          )
{
   R__ASSERT(parent);
   fCounter = parent->DeepCopy();
   fCounter2 = parent->DeepCopy();
   TFormLeafInfo ** next = &(fCounter2->fNext);
   while(*next != 0) next = &( (*next)->fNext);
   *next = new TFormLeafInfoCollectionSize(elementclassptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoMultiVarDimCollection::TFormLeafInfoMultiVarDimCollection(
   TClass* motherclassptr,
   Long_t offset,
   TStreamerElement* element,
   TFormLeafInfo *parent) :
   TFormLeafInfoMultiVarDim(motherclassptr,offset,element)
{
   R__ASSERT(parent && element);
   fCounter = parent->DeepCopy();
   fCounter2 = parent->DeepCopy();
   TFormLeafInfo ** next = &(fCounter2->fNext);
   while(*next != 0) next = &( (*next)->fNext);
   *next = new TFormLeafInfoCollectionSize(motherclassptr,offset,element);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the object and all its content.

TFormLeafInfo* TFormLeafInfoMultiVarDimCollection::DeepCopy() const
{
   return new TFormLeafInfoMultiVarDimCollection(*this);
}

////////////////////////////////////////////////////////////////////////////////

Double_t TFormLeafInfoMultiVarDimCollection::GetValue(TLeaf * /* leaf */,
                                                      Int_t /* instance */)
{
   /* The proper indexing and unwinding of index need to be done by prior leafinfo in the chain. */
   Error("GetValue","This should never be called");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Load the current array sizes.

void TFormLeafInfoMultiVarDimCollection::LoadSizes(TBranch* branch)
{
   R__ASSERT(fCounter2);

   TLeaf *leaf = (TLeaf*)branch->GetListOfLeaves()->At(0);
   fNsize = (Int_t)fCounter->GetCounterValue(leaf);

   if (fNsize > fSizes.GetSize()) fSizes.Set(fNsize);
   fSumOfSizes = 0;
   for (Int_t i=0; i<fNsize; i++) {
      Int_t size = (Int_t)fCounter2->GetValue(leaf,i);
      fSumOfSizes += size;
      fSizes.AddAt( size, i );
   }
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the value of the underlying data.

template <typename T>
T TFormLeafInfoMultiVarDimCollection::ReadValueImpl(char *where, Int_t instance)
{
   if (fSecondaryIndex>=0) {
      UInt_t len = fNext->GetArrayLength();
      if (len) {
         instance = fSecondaryIndex*len;
      } else {
         instance = fSecondaryIndex;
      }
   }
   return fNext->ReadTypedValue<T>(where,instance);
}

INSTANTIATE_READVAL(TFormLeafInfoMultiVarDimCollection);


////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoMultiVarDimClones
/// A small helper class to implement reading
/// a data member on a variable size array inside a TClonesArray object stored
/// in a TTree.  This is the version used for split access

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoMultiVarDimClones::TFormLeafInfoMultiVarDimClones(
   TClass* motherclassptr,
   Long_t offset,
   TClass* elementclassptr,
   TFormLeafInfo *parent) :
   TFormLeafInfoMultiVarDim(motherclassptr,offset,
                 new TStreamerElement("clones","in class",
                                      0,
                                      TStreamerInfo::kAny,
                                      elementclassptr
                                      ? elementclassptr->GetName()
                                      : ( motherclassptr
                                          ? motherclassptr->GetName()
                                          : "Unknwon")
                                          )
                                          )
{
   R__ASSERT(parent);
   fCounter = parent->DeepCopy();
   fCounter2 = parent->DeepCopy();
   TFormLeafInfo ** next = &(fCounter2->fNext);
   while(*next != 0) next = &( (*next)->fNext);
   *next = new TFormLeafInfoClones(elementclassptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoMultiVarDimClones::TFormLeafInfoMultiVarDimClones(
   TClass* motherclassptr,
   Long_t offset,
   TStreamerElement* element,
   TFormLeafInfo *parent) :
   TFormLeafInfoMultiVarDim(motherclassptr,offset,element)
{
   R__ASSERT(parent && element);
   fCounter = parent->DeepCopy();
   fCounter2 = parent->DeepCopy();
   TFormLeafInfo ** next = &(fCounter2->fNext);
   while(*next != 0) next = &( (*next)->fNext);
   *next = new TFormLeafInfoClones(motherclassptr,offset,element);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the object and all its data.

TFormLeafInfo* TFormLeafInfoMultiVarDimClones::DeepCopy() const
{
   return new TFormLeafInfoMultiVarDimClones(*this);
}

////////////////////////////////////////////////////////////////////////////////

Double_t TFormLeafInfoMultiVarDimClones::GetValue(TLeaf * /* leaf */,
                                                      Int_t /* instance */)
{
   /* The proper indexing and unwinding of index need to be done by prior leafinfo in the chain. */
   Error("GetValue","This should never be called");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Load the current array sizes.

void TFormLeafInfoMultiVarDimClones::LoadSizes(TBranch* branch)
{
   R__ASSERT(fCounter2);

   TLeaf *leaf = (TLeaf*)branch->GetListOfLeaves()->At(0);
   fNsize = (Int_t)fCounter->GetCounterValue(leaf);

   if (fNsize > fSizes.GetSize()) fSizes.Set(fNsize);
   fSumOfSizes = 0;
   for (Int_t i=0; i<fNsize; i++) {
      TClonesArray *clones = (TClonesArray*)fCounter2->GetValuePointer(leaf,i);
      if (clones) {
         Int_t size = clones->GetEntries();
         fSumOfSizes += size;
         fSizes.AddAt( size, i );
      }
   }
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the value of the underlying data.

template <typename T>
T TFormLeafInfoMultiVarDimClones::ReadValueImpl(char *where, Int_t instance)
{
   if (fSecondaryIndex>=0) {
      UInt_t len = fNext->GetArrayLength();
      if (len) {
         instance = fSecondaryIndex*len;
      } else {
         instance = fSecondaryIndex;
      }
   }
   return fNext->ReadTypedValue<T>(where,instance);
}

INSTANTIATE_READVAL(TFormLeafInfoMultiVarDimClones);

////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoCast
/// A small helper class to implement casting an object to
/// a different type (equivalent to dynamic_cast)

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoCast::TFormLeafInfoCast(TClass* classptr, TClass* casted) :
   TFormLeafInfo(classptr),fCasted(casted),fGoodCast(kTRUE)
{
   if (casted) { fCastedName = casted->GetName(); }
   fMultiplicity = -1;
   fIsTObject = fClass->IsTObject() && fCasted->IsLoaded();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFormLeafInfoCast::TFormLeafInfoCast(const TFormLeafInfoCast& orig) :
   TFormLeafInfo(orig)
{
   fCasted = orig.fCasted;
   fCastedName = orig.fCastedName;
   fGoodCast = orig.fGoodCast;
   fIsTObject = orig.fIsTObject;
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe swap.

void TFormLeafInfoCast::Swap(TFormLeafInfoCast &other)
{
   TFormLeafInfo::Swap(other);
   std::swap(fCasted,other.fCasted);
   std::swap(fCastedName,other.fCastedName);
   std::swap(fGoodCast,other.fGoodCast);
   std::swap(fIsTObject,other.fIsTObject);
}

////////////////////////////////////////////////////////////////////////////////
/// Exception safe assignment operator.

TFormLeafInfoCast &TFormLeafInfoCast::operator=(const TFormLeafInfoCast &other)
{
   TFormLeafInfoCast tmp(other);
   Swap(tmp);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the object and all its content.

TFormLeafInfo* TFormLeafInfoCast::DeepCopy() const
{
   return new TFormLeafInfoCast(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TFormLeafInfoCast::~TFormLeafInfoCast()
{
}

// Currently only implemented in TFormLeafInfoCast
Int_t TFormLeafInfoCast::GetNdata()
{
   // Get the number of element in the entry.

   if (!fGoodCast) return 0;
   if (fNext) return fNext->GetNdata();
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Read the value at the given memory location

template <typename T>
T TFormLeafInfoCast::ReadValueImpl(char *where, Int_t instance)
{
   if (!fNext) return 0;

   // First check that the real class inherits from the
   // casted class
   // First assume TObject ...
   if ( fIsTObject && !((TObject*)where)->InheritsFrom(fCasted) ) {
      fGoodCast = kFALSE;
      return 0;
   } else {
      // We know we have a TBranchElement and we need to find out the
      // real class name.
   }
   fGoodCast = kTRUE;
   return fNext->ReadTypedValue<T>(where,instance);
}

INSTANTIATE_READVAL(TFormLeafInfoCast);


////////////////////////////////////////////////////////////////////////////////
/// We reloading all cached information in case the underlying class
/// information has changed (for example when changing from the 'emulated'
/// class to the real class.

Bool_t TFormLeafInfoCast::Update()
{
   if (fCasted) {
      TClass * new_class = TClass::GetClass(fCastedName);
      if (new_class!=fCasted) {
         fCasted = new_class;
      }
   }
   return TFormLeafInfo::Update();
}

////////////////////////////////////////////////////////////////////////////////
/// \class TFormLeafInfoTTree
/// A small helper class to implement reading
/// from the containing TTree object itself.

TFormLeafInfoTTree::TFormLeafInfoTTree(TTree *tree, const char *alias, TTree *current) :
TFormLeafInfo( TTree::Class(), 0, 0 ), fTree(tree),fCurrent(current),fAlias(alias)
{
   // Constructor.

   if (fCurrent==0) fCurrent = fTree->GetFriend(alias);
}

TFormLeafInfoTTree::TFormLeafInfoTTree(const TFormLeafInfoTTree& orig) :
   TFormLeafInfo(orig)
{
   // Copy Constructor.
   fTree    = orig.fTree;
   fAlias   = orig.fAlias;
   fCurrent = orig.fCurrent;
}

TFormLeafInfo* TFormLeafInfoTTree::DeepCopy() const
{
   // Copy the object and all its content.

   return new TFormLeafInfoTTree(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// returns the address of the value pointed to by the
/// TFormLeafInfo.

void* TFormLeafInfoTTree::GetLocalValuePointer(TLeaf *, Int_t instance)
{
   return GetLocalValuePointer((char*)fCurrent,instance);
}

////////////////////////////////////////////////////////////////////////////////
/// Return result of a leafobject method.

template <typename T>
T TFormLeafInfoTTree::GetValueImpl(TLeaf *, Int_t instance)
{
   return ReadTypedValue<T>((char*)fCurrent,instance);
}

////////////////////////////////////////////////////////////////////////////////
/// Return result of a leafobject method.

template <typename T>
T TFormLeafInfoTTree::ReadValueImpl(char *thisobj, Int_t instance)
{
   if (fElement) return TFormLeafInfo::ReadTypedValue<T>(thisobj,instance);
   else if (fNext) return fNext->ReadTypedValue<T>(thisobj,instance);
   else return 0;
}

INSTANTIATE_GETVAL(TFormLeafInfoTTree);
INSTANTIATE_READVAL(TFormLeafInfoTTree);

////////////////////////////////////////////////////////////////////////////////
/// Update after a change of file in a chain

Bool_t TFormLeafInfoTTree::Update()
{
   if (fAlias.Length() && fAlias != fTree->GetName()) {
      fCurrent = fTree->GetFriend(fAlias.Data());
   }
   return fCurrent && TFormLeafInfo::Update();
}
