// @(#)root/meta:$Id$
// Author: Rene Brun   12/10/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TROOT.h"
#include "TStreamerElement.h"
#include "TVirtualStreamerInfo.h"
#include "TBuffer.h"
#include "TClass.h"
#include "TClassEdit.h"
#include "TClassStreamer.h"
#include "TClassTable.h"
#include "TBaseClass.h"
#include "TDataMember.h"
#include "TDataType.h"
#include "TRealData.h"
#include "ThreadLocalStorage.h"
#include "TList.h"
#include "TRef.h"
#include "TInterpreter.h"
#include "TError.h"
#include "TObjArray.h"
#include "TVirtualMutex.h"
#include "TVirtualCollectionProxy.h"
#include "strlcpy.h"
#include "snprintf.h"

#include <string>

using namespace std;

const Int_t kMaxLen = 1024;

static TString &IncludeNameBuffer() {
   TTHREAD_TLS_DECL_ARG(TString,includeName,kMaxLen);
   return includeName;
}

static TString ExtractClassName(const TString &type_name)
{
   TString className = type_name.Strip(TString::kTrailing, '*');
   if (className.Index("const ")==0) className.Remove(0,6);
   return className;
}
////////////////////////////////////////////////////////////////////////////////
/// Helper function to initialize the 'index/counter' value of
/// the Pointer streamerElements.  If directive is a StreamerInfo and it correspond to the
/// same class a 'countClass' the streamerInfo is used instead of the current StreamerInfo of the TClass
/// for 'countClass'.

static TStreamerBasicType *InitCounter(const char *countClass, const char *countName, TVirtualStreamerInfo *directive)
{
   TStreamerBasicType *counter = 0;

   TClass *cl = TClass::GetClass(countClass);

   if (directive) {

      if (directive->GetClass() == cl) {
         // The info we have been passed is indeed describing the counter holder, just look there.

         TStreamerElement *element = (TStreamerElement *)directive->GetElements()->FindObject(countName);
         if (!element) return 0;
         if (element->IsA() != TStreamerBasicType::Class()) return 0;
         counter = (TStreamerBasicType*)element;

      } else {
         if (directive->GetClass()->GetListOfRealData()) {
            TRealData* rdCounter = (TRealData*) directive->GetClass()->GetListOfRealData()->FindObject(countName);
            if (!rdCounter) return 0;
            TDataMember *dmCounter = rdCounter->GetDataMember();
            cl = dmCounter->GetClass();
         } else {
            TStreamerElement *element = (TStreamerElement *)directive->GetElements()->FindObject(countName);
            if (!element) return 0;
            if (element->IsA() != TStreamerBasicType::Class()) return 0;
            cl = directive->GetClass();
         }
         if (cl==0) return 0;
         counter = TVirtualStreamerInfo::GetElementCounter(countName,cl);
      }
   } else {

      if (cl==0) return 0;
      counter = TVirtualStreamerInfo::GetElementCounter(countName,cl);
   }

   //at this point the counter may be declared to be skipped
   if (counter) {
      if (counter->GetType() < TVirtualStreamerInfo::kCounter) counter->SetType(TVirtualStreamerInfo::kCounter);
   }
   return counter;
}

////////////////////////////////////////////////////////////////////////////////
/// Parse comments to search for a range specifier of the style:
///  [xmin,xmax] or [xmin,xmax,nbits]
///  [0,1]
///  [-10,100];
///  [-pi,pi], [-pi/2,pi/4],[-2pi,2*pi]
///  [-10,100,16]
///  [0,0,8]
/// if nbits is not specified, or nbits <2 or nbits>32 it is set to 32
/// if (xmin==0 and xmax==0 and nbits <=16) the double word will be converted
/// to a float and its mantissa truncated to nbits significative bits.
///
///  see comments in TBufferFile::WriteDouble32.

static void GetRange(const char *comments, Double_t &xmin, Double_t &xmax, Double_t &factor)
{
   const Double_t kPi =3.14159265358979323846 ;
   factor = xmin = xmax = 0;
   if (!comments) return;
   const char *left = strstr(comments,"[");
   if (!left) return;
   const char *right = strstr(left,"]");
   if (!right) return;
   const char *comma = strstr(left,",");
   if (!comma || comma > right) {
      //may be first bracket was a dimension specifier
      left = strstr(right,"[");
      if (!left) return;
      right = strstr(left,"]");
      if (!right) return;
      comma = strstr(left,",");
      if (!comma || comma >right) return;
   }
   //search if nbits is specified
   const char *comma2 = 0;
   if (comma) comma2 = strstr(comma+1,",");
   if (comma2 > right) comma2 = 0;
   Int_t nbits = 32;
   if (comma2) {
      TString sbits(comma2+1,right-comma2-1);
      sscanf(sbits.Data(),"%d",&nbits);
      if (nbits < 2 || nbits > 32) {
         ::Error("GetRange","Illegal specification for the number of bits; %d. reset to 32.",nbits);
         nbits = 32;
      }
      right = comma2;
   }
   TString range(left+1,right-left-1);
   TString sxmin(left+1,comma-left-1);
   sxmin.ToLower();
   sxmin.ReplaceAll(" ","");
   if (sxmin.Contains("pi")) {
      if      (sxmin.Contains("2pi"))   xmin = 2*kPi;
      else if (sxmin.Contains("2*pi"))  xmin = 2*kPi;
      else if (sxmin.Contains("twopi")) xmin = 2*kPi;
      else if (sxmin.Contains("pi/2"))  xmin = kPi/2;
      else if (sxmin.Contains("pi/4"))  xmin = kPi/4;
      else if (sxmin.Contains("pi"))    xmin = kPi;
      if (sxmin.Contains("-"))          xmin = -xmin;
   } else {
      sscanf(sxmin.Data(),"%lg",&xmin);
   }
   TString sxmax(comma+1,right-comma-1);
   sxmax.ToLower();
   sxmax.ReplaceAll(" ","");
   if (sxmax.Contains("pi")) {
      if      (sxmax.Contains("2pi"))   xmax = 2*kPi;
      else if (sxmax.Contains("2*pi"))  xmax = 2*kPi;
      else if (sxmax.Contains("twopi")) xmax = 2*kPi;
      else if (sxmax.Contains("pi/2"))  xmax = kPi/2;
      else if (sxmax.Contains("pi/4"))  xmax = kPi/4;
      else if (sxmax.Contains("pi"))    xmax = kPi;
      if (sxmax.Contains("-"))          xmax = -xmax;
   } else {
      sscanf(sxmax.Data(),"%lg",&xmax);
   }
   UInt_t bigint;
   if (nbits < 32)  bigint = 1<<nbits;
   else             bigint = 0xffffffff;
   if (xmin < xmax) factor = bigint/(xmax-xmin);
   if (xmin >= xmax && nbits <15) xmin = nbits+0.1;
}

ClassImp(TStreamerElement);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TStreamerElement::TStreamerElement()
{
   fType        = 0;
   fSize        = 0;
   fNewType     = 0;
   fArrayDim    = 0;
   fArrayLength = 0;
   fStreamer    = 0;
   fOffset      = 0;
   fClassObject = (TClass*)(-1);
   fNewClass    = 0;
   fTObjectOffset = 0;
   fFactor      = 0;
   fXmin        = 0;
   fXmax        = 0;
   for (Int_t i=0;i<5;i++) fMaxIndex[i] = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TStreamerElement object.

TStreamerElement::TStreamerElement(const char *name, const char *title, Int_t offset, Int_t dtype, const char *typeName)
        : TNamed(name,title)
{
   fOffset      = offset;
   fType        = dtype;
   fSize        = 0;
   fNewType     = fType;
   fArrayDim    = 0;
   fArrayLength = 0;
   if (typeName && !strcmp(typeName, "BASE")) {
      // TStreamerBase case; fTypeName should stay "BASE".
      fTypeName = typeName;
   } else {
      //must protect call into the interpreter
      R__LOCKGUARD(gInterpreterMutex);
      fTypeName    = TClassEdit::ResolveTypedef(typeName);
   }
   fStreamer    = 0;
   fClassObject = (TClass*)(-1);
   fNewClass    = 0;
   fTObjectOffset = 0;
   fFactor      = 0;
   fXmin        = 0;
   fXmax        = 0;
   for (Int_t i=0;i<5;i++) fMaxIndex[i] = 0;
   if (fTypeName == "Float16_t" || fTypeName == "Float16_t*") {
      GetRange(title,fXmin,fXmax,fFactor);
      if (fFactor > 0 || fXmin > 0) SetBit(kHasRange);
   }
   if (fTypeName == "Double32_t" || fTypeName == "Double32_t*") {
      GetRange(title,fXmin,fXmax,fFactor);
      if (fFactor > 0 || fXmin > 0) SetBit(kHasRange);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TStreamerElement dtor.

TStreamerElement::~TStreamerElement()
{
}


////////////////////////////////////////////////////////////////////////////////
/// Returns true if the element cannot be split, false otherwise.
/// An element cannot be split if the corresponding class member has
/// the special characters "||" as the first characters in the
/// comment field.

Bool_t TStreamerElement::CannotSplit() const
{
   if (GetTitle()[0] != 0 && strspn(GetTitle(),"||") == 2) return kTRUE;
   TClass *cl = GetClassPointer();
   if (!cl) return kFALSE;  //basic type

   static TClassRef clonesArray("TClonesArray");
   if (IsaPointer() && cl != clonesArray && !cl->GetCollectionProxy()) return kTRUE;

   switch(fType) {
      case TVirtualStreamerInfo::kAny    +TVirtualStreamerInfo::kOffsetL:
      case TVirtualStreamerInfo::kObject +TVirtualStreamerInfo::kOffsetL:
      case TVirtualStreamerInfo::kTObject+TVirtualStreamerInfo::kOffsetL:
      case TVirtualStreamerInfo::kTString+TVirtualStreamerInfo::kOffsetL:
      case TVirtualStreamerInfo::kTNamed +TVirtualStreamerInfo::kOffsetL:
         return kTRUE;
   }

   if ( !cl->CanSplit() ) return kTRUE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a pointer to the TClass of this element.

TClass *TStreamerElement::GetClassPointer() const
{
   if (fClassObject!=(TClass*)(-1)) return fClassObject;

   TString className(ExtractClassName(fTypeName));
   bool quiet = (fType == TVirtualStreamerInfo::kArtificial);
   ((TStreamerElement*)this)->fClassObject = TClass::GetClass(className, kTRUE, quiet);
   return fClassObject;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the TExec id for the EXEC instruction in the comment field
/// of a TRef data member.

Int_t TStreamerElement::GetExecID() const
{
   //check if element is a TRef or TRefArray
   if (strncmp(fTypeName.Data(),"TRef",4) != 0) return 0;

   //if the UniqueID of this element has already been set, we assume
   //that it contains the exec id of a TRef object.
   if (GetUniqueID()) return GetUniqueID();

   //check if an Exec is specified in the comment field
   char *action = (char*)strstr(GetTitle(),"EXEC:");
   if (!action) return 0;
   Int_t nch = strlen(action)+1;
   char *caction = new char[nch];
   strlcpy(caction,action+5,nch);
   char *blank = (char*)strchr(caction,' ');
   if (blank) *blank = 0;
   //we have found the Exec name in the comment
   //we register this Exec to the list of Execs.
   Int_t index = TRef::AddExec(caction);
   delete [] caction;
   //we save the Exec index as the uniqueid of this STreamerElement
   const_cast<TStreamerElement*>(this)->SetUniqueID(index+1);
   return index+1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return element name including dimensions, if any
/// Note that this function stores the name into a static array.
/// You should copy the result.

const char *TStreamerElement::GetFullName() const
{
   TTHREAD_TLS_DECL_ARG(TString,name,kMaxLen);
   char cdim[20];
   name = GetName();
   for (Int_t i=0;i<fArrayDim;i++) {
      snprintf(cdim,19,"[%d]",fMaxIndex[i]);
      name += cdim;
   }
   return name;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill type with the string representation of sequence
/// information including 'cached','repeat','write' or
/// 'nodelete'.

void TStreamerElement::GetSequenceType(TString &sequenceType) const
{
   sequenceType.Clear();
   Bool_t first = kTRUE;
   if (TestBit(TStreamerElement::kWholeObject)) {
      if (!first) sequenceType += ",";
      first = kFALSE;
      sequenceType += "wholeObject";
   }
   if (TestBit(TStreamerElement::kCache)) {
      first = kFALSE;
      sequenceType += "cached";
   }
   if (TestBit(TStreamerElement::kRepeat)) {
      if (!first) sequenceType += ",";
      first = kFALSE;
      sequenceType += "repeat";
   }
   if (TestBit(TStreamerElement::kDoNotDelete)) {
      if (!first) sequenceType += ",";
      first = kFALSE;
      sequenceType += "nodelete";
   }
   if (TestBit(TStreamerElement::kWrite)) {
      if (!first) sequenceType += ",";
      first = kFALSE;
      sequenceType += "write";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns size of this element in bytes.

Int_t TStreamerElement::GetSize() const
{
   return fSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the local streamer object.

TMemberStreamer *TStreamerElement::GetStreamer() const
{
   return fStreamer;
}

////////////////////////////////////////////////////////////////////////////////
/// Return type name of this element
/// in case the type name is not a standard basic type, return
/// the basic type name known to CINT.

const char *TStreamerElement::GetTypeNameBasic() const
{
   TDataType *dt = gROOT->GetType(fTypeName.Data());
   if (fType < 1 || fType > 55) return fTypeName.Data();
   if (dt && dt->GetType() > 0) return fTypeName.Data();
   Int_t dtype = fType%20;
   return TDataType::GetTypeName((EDataType)dtype);
}

////////////////////////////////////////////////////////////////////////////////
/// Initliaze the element.

void TStreamerElement::Init(TVirtualStreamerInfo *)
{
   fClassObject = GetClassPointer();
   if (fClassObject && fClassObject->IsTObject()) {
      fTObjectOffset = fClassObject->GetBaseClassOffset(TObject::Class());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// The early 3.00/00 and 3.01/01 versions used to store
/// dm->GetTypeName instead of dm->GetFullTypename
/// if this case is detected, the element type name is modified.

Bool_t TStreamerElement::IsOldFormat(const char *newTypeName)
{
   //if (!IsaPointer()) return kFALSE;
   if (!strstr(newTypeName,fTypeName.Data())) return kFALSE;
   //if (!strstr(fTypeName.Data(),newTypeName)) return kFALSE;
   fTypeName = newTypeName;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if the element represent a base class.

Bool_t TStreamerElement::IsBase() const
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if the element represent an entity that is not written
/// to the disk (transient members, cache allocator/deallocator, etc.)

Bool_t TStreamerElement::IsTransient() const
{
   if (fType == TVirtualStreamerInfo::kArtificial) {
      // if (((const TStreamerArtificial*)this)->GetWriteFunc() == 0)
         return kTRUE;
   }
   if (fType == TVirtualStreamerInfo::kCacheNew) return kTRUE;
   if (fType == TVirtualStreamerInfo::kCacheDelete) return kTRUE;
   if (fType == TVirtualStreamerInfo::kCache) return kTRUE;
   if (fType == TVirtualStreamerInfo::kMissing) return kTRUE;
   if (TVirtualStreamerInfo::kSkip <= fType && fType < TVirtualStreamerInfo::kConv) return kTRUE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Print the content of the element.

void TStreamerElement::ls(Option_t *) const
{
   TString temp(GetTypeName());
   if (IsaPointer() && !fTypeName.Contains("*")) temp += "*";

   TString sequenceType;
   GetSequenceType(sequenceType);
   if (sequenceType.Length()) {
      sequenceType.Prepend(" (");
      sequenceType += ") ";
   }
   printf("  %-14s %-15s offset=%3d type=%2d %s%-20s\n",
          temp.Data(),GetFullName(),fOffset,fType,sequenceType.Data(),
          GetTitle());
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of array dimensions.

void TStreamerElement::SetArrayDim(Int_t dim)
{
   fArrayDim = dim;
   if (dim) fType += TVirtualStreamerInfo::kOffsetL;
   fNewType = fType;
}

////////////////////////////////////////////////////////////////////////////////
///set maximum index for array with dimension dim

void TStreamerElement::SetMaxIndex(Int_t dim, Int_t max)
{
   if (dim < 0 || dim > 4) return;
   fMaxIndex[dim] = max;
   if (fArrayLength == 0)  fArrayLength  = max;
   else                    fArrayLength *= max;
}

////////////////////////////////////////////////////////////////////////////////
///set pointer to Streamer function for this element

void TStreamerElement::SetStreamer(TMemberStreamer *streamer)
{
   fStreamer = streamer;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStreamerElement.

void TStreamerElement::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      //NOTE that when reading, one cannot use Class()->ReadBuffer
      // TBuffer::Class methods used for reading streamerinfos from SQL database
      // Any changes of class structure should be reflected by them starting from version 4

      R__b.ClassBegin(TStreamerElement::Class(), R__v);
      R__b.ClassMember("TNamed");
      TNamed::Streamer(R__b);
      R__b.ClassMember("fType","Int_t");
      R__b >> fType;
      R__b.ClassMember("fSize","Int_t");
      R__b >> fSize;
      R__b.ClassMember("fArrayLength","Int_t");
      R__b >> fArrayLength;
      R__b.ClassMember("fArrayDim","Int_t");
      R__b >> fArrayDim;
      R__b.ClassMember("fMaxIndex","Int_t", 5);
      if (R__v == 1) R__b.ReadStaticArray(fMaxIndex);
      else           R__b.ReadFastArray(fMaxIndex,5);
      R__b.ClassMember("fTypeName","TString");
      fTypeName.Streamer(R__b);
      if (fType==11&&(fTypeName=="Bool_t"||fTypeName=="bool")) fType = 18;
      if (R__v > 1) {
         SetUniqueID(0);
         //check if element is a TRef or TRefArray
         GetExecID();
      }
      if (R__v <= 2 && this->IsA()==TStreamerBasicType::Class()) {
         // In TStreamerElement v2, fSize was holding the size of
         // the underlying data type.  In later version it contains
         // the full length of the data member.
         TDataType *type = gROOT->GetType(GetTypeName());
         if (type && fArrayLength) fSize = fArrayLength * type->Size();
      }
      if (R__v == 3) {
         R__b >> fXmin;
         R__b >> fXmax;
         R__b >> fFactor;
         if (fFactor > 0) SetBit(kHasRange);
      }
      if (R__v > 3) {
         if (TestBit(kHasRange)) GetRange(GetTitle(),fXmin,fXmax,fFactor);
      }
      //R__b.CheckByteCount(R__s, R__c, TStreamerElement::IsA());
      R__b.ClassEnd(TStreamerElement::Class());
      R__b.SetBufferOffset(R__s+R__c+sizeof(UInt_t));

      ResetBit(TStreamerElement::kCache);
      ResetBit(TStreamerElement::kWrite);
   } else {
      R__b.WriteClassBuffer(TStreamerElement::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
///function called by the TClass constructor when replacing an emulated class
///by the real class

void TStreamerElement::Update(const TClass *oldClass, TClass *newClass)
{
   if (fClassObject == oldClass) {
      fClassObject = newClass;
      if (fClassObject && fClassObject->IsTObject()) {
         fTObjectOffset = fClassObject->GetBaseClassOffset(TObject::Class());
      }
   } else if (fClassObject == nullptr) {
      // Well since some emulated class is replaced by a real class, we can
      // assume a new library has been loaded.  If this is the case, we should
      // check whether the class now exist (this would be the case for example
      // for reading STL containers).

      TString classname(ExtractClassName(fTypeName));

      if (classname == newClass->GetName()) {
         fClassObject = newClass;
         if (fClassObject && fClassObject->IsTObject()) {
            fTObjectOffset = fClassObject->GetBaseClassOffset(TObject::Class());
         }
      } else if (TClassTable::GetDict(classname)) {
         fClassObject = (TClass*)-1;
         GetClassPointer(); //force fClassObject
         if (fClassObject && fClassObject->IsTObject()) {
            fTObjectOffset = fClassObject->GetBaseClassOffset(TObject::Class());
         }
      }
   }
}

//______________________________________________________________________________

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStreamerBase implement the streamer of the base class               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

ClassImp(TStreamerBase);

////////////////////////////////////////////////////////////////////////////////

TStreamerBase::TStreamerBase() :
   // Abuse TStreamerElement data member that is not used by TStreamerBase
   fBaseCheckSum( *( (UInt_t*)&(fMaxIndex[1]) ) ),
   fStreamerFunc(0), fConvStreamerFunc(0), fStreamerInfo(0)
{
   // Default ctor.

   fBaseClass = (TClass*)(-1);
   fBaseVersion = 0;
   fNewBaseClass = 0;
}

////////////////////////////////////////////////////////////////////////////////

TStreamerBase::TStreamerBase(const char *name, const char *title, Int_t offset, Bool_t isTransient)
   : TStreamerElement(name,title,offset,TVirtualStreamerInfo::kBase,"BASE"),
     // Abuse TStreamerElement data member that is not used by TStreamerBase
     fBaseCheckSum( *( (UInt_t*)&(fMaxIndex[1]) ) ),
     fStreamerFunc(0), fConvStreamerFunc(0), fStreamerInfo(0)

{
   // Create a TStreamerBase object.

   if (strcmp(name,"TObject") == 0) fType = TVirtualStreamerInfo::kTObject;
   if (strcmp(name,"TNamed")  == 0) fType = TVirtualStreamerInfo::kTNamed;
   fNewType = fType;
   fBaseClass = TClass::GetClass(GetName());
   if (fBaseClass) {
      if (fBaseClass->IsVersioned()) {
         fBaseVersion = fBaseClass->GetClassVersion();
      } else {
         fBaseVersion = -1;
      }
      fBaseCheckSum = fBaseClass->GetCheckSum();
   } else {
      fBaseVersion = 0;
   }
   fNewBaseClass = 0;
   Init(isTransient);
}

////////////////////////////////////////////////////////////////////////////////
/// TStreamerBase dtor

TStreamerBase::~TStreamerBase()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a pointer to the TClass of this element.

TClass *TStreamerBase::GetClassPointer() const
{
   if (fBaseClass!=(TClass*)(-1)) return fBaseClass;
   ((TStreamerBase*)this)->fBaseClass = TClass::GetClass(GetName());
   return fBaseClass;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns size of baseclass in bytes.

Int_t TStreamerBase::GetSize() const
{
   TClass *cl = GetClassPointer();
   if (cl) return cl->Size();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Setup the element.

void TStreamerBase::Init(TVirtualStreamerInfo *)
{
   Init(kFALSE);
}

void TStreamerBase::Init(Bool_t isTransient)
{
   fBaseClass = TClass::GetClass(GetName());
   if (!fBaseClass) return;

   InitStreaming(isTransient);
}

////////////////////////////////////////////////////////////////////////////////
/// Setup the fStreamerFunc and fStreamerinfo

void TStreamerBase::InitStreaming(Bool_t isTransient)
{
   if (fNewBaseClass) {
      fStreamerFunc = fNewBaseClass->GetStreamerFunc();
      fConvStreamerFunc = fNewBaseClass->GetConvStreamerFunc();
      if (fBaseVersion > 0 || fBaseCheckSum == 0) {
         fStreamerInfo = fNewBaseClass->GetConversionStreamerInfo(fBaseClass,fBaseVersion);
      } else {
         fStreamerInfo = fNewBaseClass->FindConversionStreamerInfo(fBaseClass,fBaseCheckSum);
      }
   } else if (fBaseClass && fBaseClass != (TClass*)-1) {
      fStreamerFunc = fBaseClass->GetStreamerFunc();
      fConvStreamerFunc = fBaseClass->GetConvStreamerFunc();
      if (fBaseVersion >= 0 || fBaseCheckSum == 0) {
         fStreamerInfo = fBaseClass->GetStreamerInfo(fBaseVersion, isTransient);
      } else {
         fStreamerInfo = fBaseClass->FindStreamerInfo(fBaseCheckSum, isTransient);
      }
   } else {
      fStreamerFunc = 0;
      fConvStreamerFunc = 0;
      fStreamerInfo = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if the element represent a base class.

Bool_t TStreamerBase::IsBase() const
{
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the proper include for this element.

const char *TStreamerBase::GetInclude() const
{
   if (GetClassPointer() && fBaseClass->HasInterpreterInfo()) {
      IncludeNameBuffer().Form("\"%s\"",fBaseClass->GetDeclFileName());
   } else {
      std::string shortname( TClassEdit::ShortType( GetName(), 1 ) );
      IncludeNameBuffer().Form("\"%s.h\"",shortname.c_str());
   }
   return IncludeNameBuffer();
}

////////////////////////////////////////////////////////////////////////////////
/// Print the content of the element.

void TStreamerBase::ls(Option_t *) const
{
   TString sequenceType;
   GetSequenceType(sequenceType);
   if (sequenceType.Length()) {
      sequenceType.Prepend(" (");
      sequenceType += ") ";
   }
   printf("  %-14s %-15s offset=%3d type=%2d %s%-20s\n",GetFullName(),GetTypeName(),fOffset,fType,sequenceType.Data(),GetTitle());
}

////////////////////////////////////////////////////////////////////////////////
/// Read the content of the buffer.

Int_t TStreamerBase::ReadBuffer (TBuffer &b, char *pointer)
{
   if (fConvStreamerFunc) {
      // We have a custom Streamer member function, we must use it.
      fConvStreamerFunc(b,pointer+fOffset,fNewBaseClass ? fBaseClass : nullptr);
   } else if (fStreamerFunc) {
      // We have a custom Streamer member function, we must use it.
      fStreamerFunc(b,pointer+fOffset);
   } else {
      // We don't have a custom Streamer member function. That still doesn't mean
      // that there is no streamer - it could be an external one:
      // If the old base class has an adopted streamer we take that
      // one instead of the new base class:
      if( fNewBaseClass ) {
         TClassStreamer* extstrm = fNewBaseClass->GetStreamer();
         if (extstrm) {
            // The new base class has an adopted streamer:
            extstrm->SetOnFileClass(fBaseClass);
            (*extstrm)(b, pointer);
         } else {
            b.ReadClassBuffer( fNewBaseClass, pointer+fOffset, fBaseClass );
         }
      } else {
         TClassStreamer* extstrm = fBaseClass->GetStreamer();
         if (extstrm) {
            // The class has an adopted streamer:
            (*extstrm)(b, pointer);
         } else {
            b.ReadClassBuffer( fBaseClass, pointer+fOffset );
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStreamerBase.

void TStreamerBase::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);

      R__b.ClassBegin(TStreamerBase::Class(), R__v);

      R__b.ClassMember("TStreamerElement");
      TStreamerElement::Streamer(R__b);
      // If the class owning the TStreamerElement and the base class are not
      // loaded, on the file their streamer info might be in the following
      // order (derived class,base class) and hence the base class is not
      // yet emulated.
      fBaseClass = (TClass*)-1;
      fNewBaseClass = 0;
      // Eventually we need a v3 that stores directly fBaseCheckSum (and
      // a version of TStreamerElement should not stored fMaxIndex)
      if (R__v > 2) {
         R__b.ClassMember("fBaseVersion","Int_t");
         R__b >> fBaseVersion;
      } else {
         // could have been: fBaseVersion = GetClassPointer()->GetClassVersion();
         fBaseClass = TClass::GetClass(GetName());
         fBaseVersion = fBaseClass->GetClassVersion();
      }
      R__b.ClassEnd(TStreamerBase::Class());
      R__b.SetBufferOffset(R__s+R__c+sizeof(UInt_t));
   } else {
      R__b.WriteClassBuffer(TStreamerBase::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
///Function called by the TClass constructor when replacing an emulated class
///by the real class.

void TStreamerBase::Update(const TClass *oldClass, TClass *newClass)
{
   TStreamerElement::Update(oldClass, newClass);

   if (fBaseClass == oldClass) {
      fBaseClass = newClass;
      InitStreaming(kFALSE);
   } else if (fBaseClass == nullptr) {
      if (fName == newClass->GetName()) {
         fBaseClass = newClass;
         InitStreaming(kFALSE);
      } else if (TClassTable::GetDict(fName)) {
         fBaseClass = TClass::GetClass(fName);
         InitStreaming(kFALSE);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write the base class into the buffer.

Int_t TStreamerBase::WriteBuffer (TBuffer &b, char *pointer)
{
   if (fStreamerFunc) {
      // We have a custom Streamer member function, we must use it.
      fStreamerFunc(b,pointer+fOffset);
   } else {
      // We don't have a custom Streamer member function. That still doesn't mean
      // that there is no streamer - it could be an external one:
      // If the old base class has an adopted streamer we take that
      // one instead of the new base class:
      if (fNewBaseClass) {
         TClassStreamer* extstrm = fNewBaseClass->GetStreamer();
         if (extstrm) {
            // The new base class has an adopted streamer:
            extstrm->SetOnFileClass(fBaseClass);
            (*extstrm)(b, pointer);
            return 0;
         } else {
            fNewBaseClass->WriteBuffer(b,pointer+fOffset);
            return 0;
         }
      } else {
         TClassStreamer* extstrm = fBaseClass->GetStreamer();
         if (extstrm) {
            (*extstrm)(b, pointer);
            return 0;
         } else {
            fBaseClass->WriteBuffer(b,pointer+fOffset);
            return 0;
         }
      }
   }
   return 0;
}

//______________________________________________________________________________

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStreamerBasicPointer implements the streamering of pointer to       //
// fundamental types.                                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

ClassImp(TStreamerBasicPointer);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TStreamerBasicPointer::TStreamerBasicPointer() : fCountVersion(0),fCountName(),fCountClass(),fCounter(0)
{
   fCounter = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TStreamerBasicPointer object.

TStreamerBasicPointer::TStreamerBasicPointer(const char *name, const char *title, Int_t offset, Int_t dtype, const char *countName, const char *countClass, Int_t countVersion, const char *typeName)
   : TStreamerElement(name,title,offset,dtype,typeName)
{
   fType += TVirtualStreamerInfo::kOffsetP;
   fCountName    = countName;
   fCountClass   = countClass;
   fCountVersion = countVersion;  //currently unused
   Init();
//   printf("BasicPointer Init:%s, countName=%s, countClass=%s, countVersion=%d, fCounter=%x\n",
//      name,countName,countClass,countVersion,fCounter);
}

////////////////////////////////////////////////////////////////////////////////
/// TStreamerBasicPointer dtor.

TStreamerBasicPointer::~TStreamerBasicPointer()
{
}

////////////////////////////////////////////////////////////////////////////////
/// return offset of counter

ULong_t TStreamerBasicPointer::GetMethod() const
{
   if (!fCounter) ((TStreamerBasicPointer*)this)->Init();
   if (!fCounter) return 0;
   // FIXME: does not suport multiple inheritance for counter in base class.
   // This is wrong in case counter is not in the same class or one of
   // the left most (non virtual) base classes.  For the other we would
   // really need to use the object coming from the list of real data.
   // (and even that need analysis for virtual base class).
   return (ULong_t)fCounter->GetOffset();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns size of basicpointer in bytes.

Int_t TStreamerBasicPointer::GetSize() const
{
   if (fArrayLength) return fArrayLength*sizeof(void *);
   return sizeof(void *);
}

////////////////////////////////////////////////////////////////////////////////
/// Setup the element.
/// If directive is a StreamerInfo and it correspond to the
/// same class a 'countClass' the streamerInfo is used instead of the current StreamerInfo of the TClass
/// for 'countClass'.

void TStreamerBasicPointer::Init(TVirtualStreamerInfo *directive)
{
   fCounter = InitCounter( fCountClass, fCountName, directive );
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of array dimensions.

void TStreamerBasicPointer::SetArrayDim(Int_t dim)
{
   fArrayDim = dim;
   //if (dim) fType += TVirtualStreamerInfo::kOffsetL;
   fNewType = fType;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStreamerBasicPointer.

void TStreamerBasicPointer::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TStreamerBasicPointer::Class(), this, R__v, R__s, R__c);
         //Init();
         //fCounter = InitCounter( fCountClass, fCountName );
         return;
      }
      //====process old versions before automatic schema evolution
      TStreamerElement::Streamer(R__b);
      R__b >> fCountVersion;
      fCountName.Streamer(R__b);
      fCountClass.Streamer(R__b);
      R__b.SetBufferOffset(R__s+R__c+sizeof(UInt_t));
   } else {
      R__b.WriteClassBuffer(TStreamerBasicPointer::Class(),this);
   }
}


//______________________________________________________________________________

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStreamerLoop implement streaming of a few construct that require    //
// looping over the data member and are not convered by other case      //
// (most deprecated).                                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

ClassImp(TStreamerLoop);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TStreamerLoop::TStreamerLoop() : fCountVersion(0),fCountName(),fCountClass(),fCounter(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TStreamerLoop object.

TStreamerLoop::TStreamerLoop(const char *name, const char *title, Int_t offset, const char *countName, const char *countClass, Int_t countVersion, const char *typeName)
        : TStreamerElement(name,title,offset,TVirtualStreamerInfo::kStreamLoop,typeName)
{
   fCountName    = countName;
   fCountClass   = countClass;
   fCountVersion = countVersion;  //currently unused
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// TStreamerLoop dtor.

TStreamerLoop::~TStreamerLoop()
{
}

////////////////////////////////////////////////////////////////////////////////
/// return address of counter

ULong_t TStreamerLoop::GetMethod() const
{
   //if (!fCounter) {
   //   Init();
   //   if (!fCounter) return 0;
   //}
   if (!fCounter) return 0;
   return (ULong_t)fCounter->GetOffset();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns size of counter in bytes.

Int_t TStreamerLoop::GetSize() const
{
   if (fArrayLength) return fArrayLength*sizeof(void*);
   return sizeof(void*);
}

////////////////////////////////////////////////////////////////////////////////
/// Setup the element.
/// If directive is a StreamerInfo and it correspond to the
/// same class a 'countClass' the streamerInfo is used instead of the current StreamerInfo of the TClass
/// for 'countClass'.

void TStreamerLoop::Init(TVirtualStreamerInfo *directive)
{
   fCounter = InitCounter( fCountClass, fCountName, directive );
}

////////////////////////////////////////////////////////////////////////////////
/// Return the proper include for this element.

const char *TStreamerLoop::GetInclude() const
{
   IncludeNameBuffer().Form("<%s>","TString.h"); //to be generalized
   return IncludeNameBuffer();
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStreamerLoop.

void TStreamerLoop::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TStreamerLoop::Class(), this, R__v, R__s, R__c);
         //Init();
         return;
      }
      //====process old versions before automatic schema evolution
      TStreamerElement::Streamer(R__b);
      R__b >> fCountVersion;
      fCountName.Streamer(R__b);
      fCountClass.Streamer(R__b);
      R__b.SetBufferOffset(R__s+R__c+sizeof(UInt_t));
   } else {
      R__b.WriteClassBuffer(TStreamerLoop::Class(),this);
   }
}


//______________________________________________________________________________

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStreamerBasicType implement streaming of fundamental types (int,    //
// float, etc.).                                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

ClassImp(TStreamerBasicType);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TStreamerBasicType::TStreamerBasicType() : fCounter(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TStreamerBasicType object.

TStreamerBasicType::TStreamerBasicType(const char *name, const char *title, Int_t offset, Int_t dtype, const char *typeName)
        : TStreamerElement(name,title,offset,dtype,typeName),fCounter(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// TStreamerBasicType dtor.

TStreamerBasicType::~TStreamerBasicType()
{
}

////////////////////////////////////////////////////////////////////////////////
/// return address of counter

ULong_t TStreamerBasicType::GetMethod() const
{
   if (fType ==  TVirtualStreamerInfo::kCounter ||
       fType == (TVirtualStreamerInfo::kCounter+TVirtualStreamerInfo::kSkip)) return (ULong_t)&fCounter;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns size of this element in bytes.

Int_t TStreamerBasicType::GetSize() const
{
   return fSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStreamerBasicType.

void TStreamerBasicType::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TStreamerBasicType::Class(), this, R__v, R__s, R__c);
      } else {
         //====process old versions before automatic schema evolution
         TStreamerElement::Streamer(R__b);
         R__b.CheckByteCount(R__s, R__c, TStreamerBasicType::IsA());
      }
      Int_t type = fType;
      if (TVirtualStreamerInfo::kOffsetL < type && type < TVirtualStreamerInfo::kOffsetP) {
         type -= TVirtualStreamerInfo::kOffsetL;
      }
      switch(type) {
         // basic types
         case TVirtualStreamerInfo::kBool:     fSize = sizeof(Bool_t);    break;
         case TVirtualStreamerInfo::kShort:    fSize = sizeof(Short_t);   break;
         case TVirtualStreamerInfo::kInt:      fSize = sizeof(Int_t);     break;
         case TVirtualStreamerInfo::kLong:     fSize = sizeof(Long_t);    break;
         case TVirtualStreamerInfo::kLong64:   fSize = sizeof(Long64_t);  break;
         case TVirtualStreamerInfo::kFloat:    fSize = sizeof(Float_t);   break;
         case TVirtualStreamerInfo::kFloat16:  fSize = sizeof(Float_t);   break;
         case TVirtualStreamerInfo::kDouble:   fSize = sizeof(Double_t);  break;
         case TVirtualStreamerInfo::kDouble32: fSize = sizeof(Double_t);  break;
         case TVirtualStreamerInfo::kUChar:    fSize = sizeof(UChar_t);   break;
         case TVirtualStreamerInfo::kUShort:   fSize = sizeof(UShort_t);  break;
         case TVirtualStreamerInfo::kUInt:     fSize = sizeof(UInt_t);    break;
         case TVirtualStreamerInfo::kULong:    fSize = sizeof(ULong_t);   break;
         case TVirtualStreamerInfo::kULong64:  fSize = sizeof(ULong64_t); break;
         case TVirtualStreamerInfo::kBits:     fSize = sizeof(UInt_t);    break;
         case TVirtualStreamerInfo::kCounter:  fSize = sizeof(Int_t);     break;
         case TVirtualStreamerInfo::kChar:     fSize = sizeof(Char_t);    break;
         case TVirtualStreamerInfo::kCharStar: fSize = sizeof(Char_t*);   break;
         default:          return; // If we don't change the size let's not remultiply it.
      }
      if (fArrayLength) fSize *= GetArrayLength();
   } else {
      R__b.WriteClassBuffer(TStreamerBasicType::Class(),this);
   }
}



//______________________________________________________________________________

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStreamerObject implements streaming of embedded objects whose type  //
// inherits from TObject.                                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

ClassImp(TStreamerObject);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TStreamerObject::TStreamerObject()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TStreamerObject object.

TStreamerObject::TStreamerObject(const char *name, const char *title, Int_t offset, const char *typeName)
        : TStreamerElement(name,title,offset,0,typeName)
{
   fType = TVirtualStreamerInfo::kObject;
   if (strcmp(typeName,"TObject") == 0) fType = TVirtualStreamerInfo::kTObject;
   if (strcmp(typeName,"TNamed")  == 0) fType = TVirtualStreamerInfo::kTNamed;
   fNewType = fType;
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// TStreamerObject dtor.

TStreamerObject::~TStreamerObject()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Setup the element.

void TStreamerObject::Init(TVirtualStreamerInfo *)
{
   fClassObject = GetClassPointer();
   if (fClassObject && fClassObject->IsTObject()) {
      fTObjectOffset = fClassObject->GetBaseClassOffset(TObject::Class());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the proper include for this element.

const char *TStreamerObject::GetInclude() const
{
   TClass *cl = GetClassPointer();
   if (cl && cl->HasInterpreterInfo()) {
      IncludeNameBuffer().Form("\"%s\"",cl->GetDeclFileName());
   } else {
      std::string shortname( TClassEdit::ShortType( GetTypeName(), 1 ) );
      IncludeNameBuffer().Form("\"%s.h\"",shortname.c_str());
   }
   return IncludeNameBuffer();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns size of object class in bytes.

Int_t TStreamerObject::GetSize() const
{
   TClass *cl = GetClassPointer();
   Int_t classSize = 8;
   if (cl) classSize = cl->Size();
   if (fArrayLength) return fArrayLength*classSize;
   return classSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStreamerObject.

void TStreamerObject::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TStreamerObject::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TStreamerElement::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TStreamerObject::IsA());
   } else {
      R__b.WriteClassBuffer(TStreamerObject::Class(),this);
   }
}


//______________________________________________________________________________

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStreamerObjectAny implement streaming of embedded object not        //
// inheriting from TObject.                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

ClassImp(TStreamerObjectAny);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TStreamerObjectAny::TStreamerObjectAny()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TStreamerObjectAny object.

TStreamerObjectAny::TStreamerObjectAny(const char *name, const char *title, Int_t offset, const char *typeName)
        : TStreamerElement(name,title,offset,TVirtualStreamerInfo::kAny,typeName)
{
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// TStreamerObjectAny dtor.

TStreamerObjectAny::~TStreamerObjectAny()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Setup the element.

void TStreamerObjectAny::Init(TVirtualStreamerInfo *)
{
   fClassObject = GetClassPointer();
   if (fClassObject && fClassObject->IsTObject()) {
      fTObjectOffset = fClassObject->GetBaseClassOffset(TObject::Class());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the proper include for this element.

const char *TStreamerObjectAny::GetInclude() const
{
   TClass *cl = GetClassPointer();
   if (cl && cl->HasInterpreterInfo()) {
      IncludeNameBuffer().Form("\"%s\"",cl->GetDeclFileName());
   } else {
      std::string shortname( TClassEdit::ShortType( GetTypeName(), 1 ) );
      IncludeNameBuffer().Form("\"%s.h\"",shortname.c_str());
   }
   return IncludeNameBuffer();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns size of anyclass in bytes.

Int_t TStreamerObjectAny::GetSize() const
{
   TClass *cl = GetClassPointer();
   Int_t classSize = 8;
   if (cl) classSize = cl->Size();
   if (fArrayLength) return fArrayLength*classSize;
   return classSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStreamerObjectAny.

void TStreamerObjectAny::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TStreamerObjectAny::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TStreamerElement::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TStreamerObjectAny::IsA());
   } else {
      R__b.WriteClassBuffer(TStreamerObjectAny::Class(),this);
   }
}



//______________________________________________________________________________

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStreamerObjectPointer implements streaming of pointer to object     //
// inheriting from TObject.                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

ClassImp(TStreamerObjectPointer);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TStreamerObjectPointer::TStreamerObjectPointer()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TStreamerObjectPointer object.

TStreamerObjectPointer::TStreamerObjectPointer(const char *name, const char *title,
                                               Int_t offset, const char *typeName)
   : TStreamerElement(name,title,offset,TVirtualStreamerInfo::kObjectP,typeName)
{
   if (strncmp(title,"->",2) == 0) fType = TVirtualStreamerInfo::kObjectp;
   fNewType = fType;
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// TStreamerObjectPointer dtor.

TStreamerObjectPointer::~TStreamerObjectPointer()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Setup the element.

void TStreamerObjectPointer::Init(TVirtualStreamerInfo *)
{
   fClassObject = GetClassPointer();
   if (fClassObject && fClassObject->IsTObject()) {
      fTObjectOffset = fClassObject->GetBaseClassOffset(TObject::Class());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the proper include for this element.

const char *TStreamerObjectPointer::GetInclude() const
{
   TClass *cl = GetClassPointer();
   if (cl && cl->HasInterpreterInfo()) {
      IncludeNameBuffer().Form("\"%s\"",cl->GetDeclFileName());
   } else {
      std::string shortname( TClassEdit::ShortType( GetTypeName(), 1 ) );
      IncludeNameBuffer().Form("\"%s.h\"",shortname.c_str());
   }

   return IncludeNameBuffer();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns size of objectpointer in bytes.

Int_t TStreamerObjectPointer::GetSize() const
{
   if (fArrayLength) return fArrayLength*sizeof(void *);
   return sizeof(void *);
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of array dimensions.

void TStreamerObjectPointer::SetArrayDim(Int_t dim)
{
   fArrayDim = dim;
   //if (dim) fType += TVirtualStreamerInfo::kOffsetL;
   fNewType = fType;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStreamerObjectPointer.

void TStreamerObjectPointer::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TStreamerObjectPointer::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TStreamerElement::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TStreamerObjectPointer::IsA());
   } else {
      R__b.WriteClassBuffer(TStreamerObjectPointer::Class(),this);
   }
}


//______________________________________________________________________________

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStreamerObjectPointerAny implements streaming of pointer to object  //
// not inheriting from TObject.                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

ClassImp(TStreamerObjectAnyPointer);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TStreamerObjectAnyPointer::TStreamerObjectAnyPointer()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TStreamerObjectAnyPointer object.

TStreamerObjectAnyPointer::TStreamerObjectAnyPointer(const char *name, const char *title,
                                                     Int_t offset, const char *typeName)
   : TStreamerElement(name,title,offset,TVirtualStreamerInfo::kAnyP,typeName)
{
   if (strncmp(title,"->",2) == 0) fType = TVirtualStreamerInfo::kAnyp;
   fNewType = fType;
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// TStreamerObjectAnyPointer dtor.

TStreamerObjectAnyPointer::~TStreamerObjectAnyPointer()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Setup the element.

void TStreamerObjectAnyPointer::Init(TVirtualStreamerInfo *)
{
   fClassObject = GetClassPointer();
   if (fClassObject && fClassObject->IsTObject()) {
      fTObjectOffset = fClassObject->GetBaseClassOffset(TObject::Class());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the proper include for this element.

const char *TStreamerObjectAnyPointer::GetInclude() const
{
   TClass *cl = GetClassPointer();
   if (cl && cl->HasInterpreterInfo()) {
      IncludeNameBuffer().Form("\"%s\"",cl->GetDeclFileName());
   } else {
      std::string shortname( TClassEdit::ShortType( GetTypeName(), 1 ) );
      IncludeNameBuffer().Form("\"%s.h\"",shortname.c_str());
   }

   return IncludeNameBuffer();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns size of objectpointer in bytes.

Int_t TStreamerObjectAnyPointer::GetSize() const
{
   if (fArrayLength) return fArrayLength*sizeof(void *);
   return sizeof(void *);
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of array dimensions.

void TStreamerObjectAnyPointer::SetArrayDim(Int_t dim)
{
   fArrayDim = dim;
   //if (dim) fType += TVirtualStreamerInfo::kOffsetL;
   fNewType = fType;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStreamerObjectAnyPointer.

void TStreamerObjectAnyPointer::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TStreamerObjectAnyPointer::Class(), this);
   } else {
      R__b.WriteClassBuffer(TStreamerObjectAnyPointer::Class(),this);
   }
}


//______________________________________________________________________________

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TSreamerString implements streaming of TString.                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

ClassImp(TStreamerString);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TStreamerString::TStreamerString()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TStreamerString object.

TStreamerString::TStreamerString(const char *name, const char *title, Int_t offset)
        : TStreamerElement(name,title,offset,TVirtualStreamerInfo::kTString,"TString")
{
}

////////////////////////////////////////////////////////////////////////////////
/// TStreamerString dtor.

TStreamerString::~TStreamerString()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Return the proper include for this element.

const char *TStreamerString::GetInclude() const
{
   IncludeNameBuffer().Form("<%s>","TString.h");
   return IncludeNameBuffer();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns size of anyclass in bytes.

Int_t TStreamerString::GetSize() const
{
   if (fArrayLength) return fArrayLength*sizeof(TString);
   return sizeof(TString);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStreamerString.

void TStreamerString::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TStreamerString::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TStreamerElement::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TStreamerString::IsA());
   } else {
      R__b.WriteClassBuffer(TStreamerString::Class(),this);
   }
}

//______________________________________________________________________________

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStreamerSTL implements streamer of STL container.                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

ClassImp(TStreamerSTL);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TStreamerSTL::TStreamerSTL() : fSTLtype(0),fCtype(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TStreamerSTL object.

TStreamerSTL::TStreamerSTL(const char *name, const char *title, Int_t offset,
                           const char *typeName, const TVirtualCollectionProxy &proxy, Bool_t dmPointer)
        : TStreamerElement(name,title,offset,ROOT::kSTLany,typeName)
{
   fTypeName = TClassEdit::ShortType(fTypeName,TClassEdit::kDropStlDefault).c_str();

  if (name==typeName /* intentional pointer comparison */
      || strcmp(name,typeName)==0) {
      // We have a base class.
      fName = fTypeName;
   }
   fSTLtype = proxy.GetCollectionType();
   fCtype   = 0;

   if (dmPointer) fSTLtype += TVirtualStreamerInfo::kOffsetP;

   if (fSTLtype == ROOT::kSTLbitset) {
      // Nothing to check
   } else if (proxy.GetValueClass()) {
      if (proxy.HasPointers()) fCtype = TVirtualStreamerInfo::kObjectp;
      else                     fCtype = TVirtualStreamerInfo::kObject;
   } else {
      fCtype = proxy.GetType();
      if (proxy.HasPointers()) fCtype += TVirtualStreamerInfo::kOffsetP;
   }
   if (TStreamerSTL::IsaPointer()) fType = TVirtualStreamerInfo::kSTLp;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TStreamerSTL object.

TStreamerSTL::TStreamerSTL(const char *name, const char *title, Int_t offset,
                           const char *typeName, const char *trueType, Bool_t dmPointer)
        : TStreamerElement(name,title,offset,ROOT::kSTLany,typeName)
{
   const char *t = trueType;
   if (!t || !*t) t = typeName;

   fTypeName = TClassEdit::ShortType(fTypeName,TClassEdit::kDropStlDefault).c_str();

   if (name==typeName /* intentional pointer comparison */
       || strcmp(name,typeName)==0) {
      // We have a base class.
      fName = fTypeName;
   }

   Int_t nch = strlen(t);
   char *s = new char[nch+1];
   strlcpy(s,t,nch+1);
   char *sopen  = strchr(s,'<');
   if (sopen == 0) {
      Fatal("TStreamerSTL","For %s, the type name (%s) is seemingly not a template (template argument not found)", name, s);
      return;
   }
   *sopen  = 0; sopen++;
   // We are looking for the first arguments of the STL container, because
   // this arguments can be a templates we need to count the < and >
   char* current=sopen;
   for(int count = 0; *current!='\0'; current++) {
      if (*current=='<') count++;
      if (*current=='>') {
         if (count==0) break;
         count--;
      }
      if (*current==',' && count==0) break;
   }
   char *sclose = current; *sclose = 0; sclose--;
   char *sconst = strstr(sopen,"const ");
   char *sbracket = strstr(sopen,"<");
   if (sconst && (sbracket==0 || sconst < sbracket)) {
      // the string "const" may be part of the classname!
      char *pconst = sconst-1;
      if (*pconst == ' ' || *pconst == '<' || *pconst == '*' || *pconst == '\0') sopen = sconst + 5;
   }
   fSTLtype = TClassEdit::STLKind(s);
   fCtype   = 0;
   if (fSTLtype == ROOT::kNotSTL) { delete [] s; return;}
   if (dmPointer) fSTLtype += TVirtualStreamerInfo::kOffsetP;

   // find STL contained type
   while (*sopen==' ') sopen++;
   Bool_t isPointer = kFALSE;
   // Find stars outside of any template definitions in the
   // first template argument.
   char *star = strrchr(sopen,'>');
   if (star) star = strchr(star,'*');
   else star = strchr(sopen,'*');
   if (star) {
      isPointer = kTRUE;
      *star = 0;
      sclose = star - 1;
   }
   while (*sclose == ' ') {*sclose = 0; sclose--;}


   TDataType *dt = (TDataType*)gROOT->GetListOfTypes()->FindObject(sopen);
   if (fSTLtype == ROOT::kSTLbitset) {
      // Nothing to check
   } else if (dt) {
      fCtype = dt->GetType();
      if (isPointer) fCtype += TVirtualStreamerInfo::kOffsetP;
   } else {
     // this could also be a nested enums ... which should work ... be let's see.
      TClass *cl = TClass::GetClass(sopen);
      if (cl) {
         if (isPointer) fCtype = TVirtualStreamerInfo::kObjectp;
         else           fCtype = TVirtualStreamerInfo::kObject;
      } else {
         if (gCling->ClassInfo_IsEnum(sopen)) {
            if (isPointer) fCtype += TVirtualStreamerInfo::kOffsetP;
         } else {
            if(strcmp(sopen,"string")) {
               // This case can happens when 'this' is a TStreamerElement for
               // a STL container containing something for which we do not have
               // a TVirtualStreamerInfo (This happens in particular is the collection
               // objects themselves are always empty) and we do not have the
               // dictionary/shared library for the container.
               if (GetClassPointer() && GetClassPointer()->IsLoaded()) {
                  Warning("TStreamerSTL","For %s we could not find any information about the type %s %d %s",fTypeName.Data(),sopen,fSTLtype,s);
               }
            }
         }
      }
   }
   delete [] s;

   if (TStreamerSTL::IsaPointer()) fType = TVirtualStreamerInfo::kSTLp;
}

////////////////////////////////////////////////////////////////////////////////
/// TStreamerSTL dtor.

TStreamerSTL::~TStreamerSTL()
{
}

////////////////////////////////////////////////////////////////////////////////
/// We can not split STL's which are inside a variable size array.
/// At least for now.

Bool_t TStreamerSTL::CannotSplit() const
{
   if (IsaPointer()) {
      if (GetTitle()[0]=='[') return kTRUE;  // can not split variable size array
      return kTRUE;
   }

   if (GetArrayDim()>=1 && GetArrayLength()>1) return kTRUE;

   if (TStreamerElement::CannotSplit()) return kTRUE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the data member is a pointer.

Bool_t TStreamerSTL::IsaPointer() const
{
   const char *type_name = GetTypeName();
   if ( type_name[strlen(type_name)-1]=='*' ) return kTRUE;
   else return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if the element represent a base class.

Bool_t TStreamerSTL::IsBase() const
{
   TString ts(GetName());

   if (strcmp(ts.Data(),GetTypeName())==0) return kTRUE;
   if (strcmp(ts.Data(),GetTypeNameBasic())==0) return kTRUE;
   return kFALSE;
}
////////////////////////////////////////////////////////////////////////////////
/// Returns size of STL container in bytes.

Int_t TStreamerSTL::GetSize() const
{
   // Since the STL collection might or might not be emulated and that the
   // sizeof the object depends on this, let's just always retrieve the
   // current size!
   TClass *cl = GetClassPointer();
   UInt_t size = 0;
   if (cl==0) {
      if (!TestBit(kWarned)) {
         Error("GetSize","Could not find the TClass for %s.\n"
               "This is likely to have been a typedef, if possible please declare it in CINT to work around the issue\n",fTypeName.Data());
         const_cast<TStreamerSTL*>(this)->SetBit(kWarned);
      }
   } else {
      size = cl->Size();
   }

   if (fArrayLength) return fArrayLength*size;
   return size;
}

////////////////////////////////////////////////////////////////////////////////
/// Print the content of the element.

void TStreamerSTL::ls(Option_t *) const
{
   TString name(kMaxLen);
   TString cdim;
   name = GetName();
   for (Int_t i=0;i<fArrayDim;i++) {
      cdim.Form("[%d]",fMaxIndex[i]);
      name += cdim;
   }
   TString sequenceType;
   GetSequenceType(sequenceType);
   if (sequenceType.Length()) {
      sequenceType.Prepend(" (");
      sequenceType += ") ";
   }
   printf("  %-14s %-15s offset=%3d type=%2d %s,stl=%d, ctype=%d, %-20s\n",
          GetTypeName(),name.Data(),fOffset,fType,sequenceType.Data(),
          fSTLtype,fCtype,GetTitle());
}

////////////////////////////////////////////////////////////////////////////////
/// Return the proper include for this element.

const char *TStreamerSTL::GetInclude() const
{
   if      (fSTLtype == ROOT::kSTLvector)            IncludeNameBuffer().Form("<%s>","vector");
   else if (fSTLtype == ROOT::kSTLlist)              IncludeNameBuffer().Form("<%s>","list");
   else if (fSTLtype == ROOT::kSTLforwardlist)       IncludeNameBuffer().Form("<%s>","forward_list");
   else if (fSTLtype == ROOT::kSTLdeque)             IncludeNameBuffer().Form("<%s>","deque");
   else if (fSTLtype == ROOT::kSTLmap)               IncludeNameBuffer().Form("<%s>","map");
   else if (fSTLtype == ROOT::kSTLmultimap)          IncludeNameBuffer().Form("<%s>","map");
   else if (fSTLtype == ROOT::kSTLset)               IncludeNameBuffer().Form("<%s>","set");
   else if (fSTLtype == ROOT::kSTLmultiset)          IncludeNameBuffer().Form("<%s>","set");
   else if (fSTLtype == ROOT::kSTLunorderedset)      IncludeNameBuffer().Form("<%s>","unordered_set");
   else if (fSTLtype == ROOT::kSTLunorderedmultiset) IncludeNameBuffer().Form("<%s>","unordered_set");
   else if (fSTLtype == ROOT::kSTLunorderedmap)      IncludeNameBuffer().Form("<%s>","unordered_map");
   else if (fSTLtype == ROOT::kSTLunorderedmultimap) IncludeNameBuffer().Form("<%s>","unordered_map");
   else if (fSTLtype == ROOT::kSTLbitset)            IncludeNameBuffer().Form("<%s>","bitset");
   return IncludeNameBuffer();
}

////////////////////////////////////////////////////////////////////////////////
/// Set pointer to Streamer function for this element
/// NOTE: we do not take ownership

void TStreamerSTL::SetStreamer(TMemberStreamer  *streamer)
{
   fStreamer = streamer;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStreamerSTL.

void TStreamerSTL::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TStreamerSTL::Class(), this, R__v, R__s, R__c);
      } else {
         //====process old versions before automatic schema evolution
         TStreamerElement::Streamer(R__b);
         R__b >> fSTLtype;
         R__b >> fCtype;
         R__b.CheckByteCount(R__s, R__c, TStreamerSTL::IsA());
      }
      if (fSTLtype == ROOT::kSTLmultimap || fSTLtype == ROOT::kSTLset) {
         // For a long time those where inverted in TStreamerElement
         // compared to the other definitions.  When we moved to version '4',
         // this got standardized, but we now need to fix it.

         if (fTypeName.BeginsWith("std::set") || fTypeName.BeginsWith("set")) {
            fSTLtype = ROOT::kSTLset;
         } else if (fTypeName.BeginsWith("std::multimap") || fTypeName.BeginsWith("multimap")) {
            fSTLtype = ROOT::kSTLmultimap;
         }
      }

      if (IsaPointer()) fType = TVirtualStreamerInfo::kSTLp;
      else fType = TVirtualStreamerInfo::kSTL;
      if (GetArrayLength() > 0) {
         fType += TVirtualStreamerInfo::kOffsetL;
      }
      if (R__b.GetParent()) { // Avoid resetting during a cloning.
         if (fCtype==TVirtualStreamerInfo::kObjectp || fCtype==TVirtualStreamerInfo::kAnyp || fCtype==TVirtualStreamerInfo::kObjectP || fCtype==TVirtualStreamerInfo::kAnyP) {
            SetBit(kDoNotDelete); // For backward compatibility
         } else if ( fSTLtype == ROOT::kSTLmap || fSTLtype == ROOT::kSTLmultimap) {
            // Here we would like to set the bit only if one of the element of the pair is a pointer,
            // however we have no easy to determine this short of parsing the class name.
            SetBit(kDoNotDelete); // For backward compatibility
         }
      }
      return;
   } else {
      // To enable forward compatibility we actually save with the old value
      TStreamerSTL tmp;
      // Hand coded copy constructor since the 'normal' one are intentionally
      // deleted.
      tmp.fName = fName;
      tmp.fTitle = fTitle;
      tmp.fType = TVirtualStreamerInfo::kStreamer;
      tmp.fSize = fSize;
      tmp.fArrayLength = fArrayLength;
      for(int i = 0; i < 5; ++i)
         tmp.fMaxIndex[i] = fMaxIndex[i];
      tmp.fTypeName = fTypeName;
      tmp.fSTLtype = fSTLtype;
      tmp.fCtype = fCtype;
      R__b.WriteClassBuffer(TStreamerSTL::Class(), &tmp);
   }
}

//______________________________________________________________________________

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStreamerSTLstring implements streaming std::string.                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

ClassImp(TStreamerSTLstring);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TStreamerSTLstring::TStreamerSTLstring()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TStreamerSTLstring object.

TStreamerSTLstring::TStreamerSTLstring(const char *name, const char *title, Int_t offset,
                                       const char *typeName, Bool_t dmPointer)
        : TStreamerSTL()
{
   SetName(name);
   SetTitle(title);

   if (dmPointer) {
      fType = TVirtualStreamerInfo::kSTLp;
   } else {
      fType = TVirtualStreamerInfo::kSTL;
   }

   fNewType = fType;
   fOffset  = offset;
   fSTLtype = ROOT::kSTLstring;
   fCtype   = ROOT::kSTLstring;
   fTypeName= typeName;

}

////////////////////////////////////////////////////////////////////////////////
/// TStreamerSTLstring dtor.

TStreamerSTLstring::~TStreamerSTLstring()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Return the proper include for this element.

const char *TStreamerSTLstring::GetInclude() const
{
   IncludeNameBuffer() = "<string>";
   return IncludeNameBuffer();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns size of anyclass in bytes.

Int_t TStreamerSTLstring::GetSize() const
{
   if (fArrayLength) return fArrayLength*sizeof(string);
   return sizeof(string);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TStreamerSTLstring.

void TStreamerSTLstring::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TStreamerSTLstring::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TStreamerSTL::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TStreamerSTLstring::IsA());
   } else {
      R__b.WriteClassBuffer(TStreamerSTLstring::Class(),this);
   }
}

//______________________________________________________________________________

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TStreamerArtificial implements StreamerElement injected by a TSchemaRule. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

ClassImp(TStreamerSTLstring);

void TStreamerArtificial::Streamer(TBuffer& /* R__b */)
{
   // Avoid streaming the synthetic/artificial streamer elements.

   // Intentionally, nothing to do at all.
   return;
}

ROOT::TSchemaRule::ReadFuncPtr_t     TStreamerArtificial::GetReadFunc()
{
   // Return the read function if any.

   return fReadFunc;
}

ROOT::TSchemaRule::ReadRawFuncPtr_t  TStreamerArtificial::GetReadRawFunc()
{
   // Return the raw read function if any.

   return fReadRawFunc;
}
