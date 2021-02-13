// @(#)root/meta:$Id$
// Author: Fons Rademakers   04/02/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TDataMember

All ROOT classes may have RTTI (run time type identification) support
added. The data is stored in so called DICTIONARY (look at TDictionary).
Information about a class is stored in TClass.
This information may be obtained via the cling api - see class TCling.
TClass has a list of TDataMember objects providing information about all
data members of described class.

\image html base_classinfo.png

TDataMember provides information about name of data member, its type,
and comment field string. It also tries to find the TMethodCall objects
responsible for getting/setting a value of it, and gives you pointers
to these methods. This gives you a unique possibility to access
protected and private (!) data members if only methods for doing that
are defined.

These methods could either be specified in a comment field, or found
out automatically by ROOT: here's an example:
suppose you have a class definition:
~~~ {.cpp}
        class MyClass{
            private:
                Float_t fX1;
                    ...
            public:
                void    SetX1(Float_t x) {fX1 = x;};
                Float_t GetX1()          {return fX1;};
                    ...
        }
~~~
Look at the data member name and method names: a data member name has
a prefix letter (f) and has a base name X1 . The methods for getting and
setting this value have names which consist of string Get/Set and the
same base name. This convention of naming data fields and methods which
access them allows TDataMember find this methods by itself completely
automatically. To make this description complete, one should know,
that names that are automatically recognized may be also:
for data fields: either fXXX or fIsXXX; and for getter function
GetXXX() or IsXXX() [where XXX is base name].

As an example of using it let's analyse a few lines which get and set
a fEditable field in TCanvas:
~~~ {.cpp}
    TCanvas     *c  = new TCanvas("c");   // create a canvas
    TClass      *cl = c->IsA();            // get its class description object.

    TDataMember *dm = cl->GetDataMember("fEditable"); //This is our data member

    TMethodCall *getter = dm->GetterMethod(c); //find a method that gets value!
    Long_t l;   // declare a storage for this value;

    getter->Execute(c,"",l);  // Get this Value !!!! It will appear in l !!!


    TMethodCall *setter = dm->SetterMethod(c);
    setter->Execute(c,"0",);   // Set Value 0 !!!
~~~

This trick is widely used in ROOT TContextMenu and dialogs for obtaining
current values and put them as initial values in dialog fields.

If you don't want to follow the convention of naming used by ROOT
you still could benefit from Getter/Setter method support: the solution
is to instruct ROOT what the names of these routines are.
The way to do it is putting this information in a comment string to a data
field in your class declaration:

~~~ {.cpp}
    class MyClass{
        Int_t mydata;  //  *OPTIONS={GetMethod="Get";SetMethod="Set"}
         ...
        Int_t Get() const { return mydata;};
        void  Set(Int_t i) {mydata=i;};
        }
~~~

However, this getting/setting functions are not the only feature of
this class. The next point is providing lists of possible settings
for the concerned data member. The idea is to have a list of possible
options for this data member, with strings identifying them. This
is used in dialogs with parameters to set - for details see
TMethodArg, TRootContextMenu, TContextMenu. This list not only specifies
the allowed value, but also provides strings naming the options.
Options are managed via TList of TOptionListItem objects. This list
is also  created automatically: if a data type is an enum type,
the list will have items describing every enum value, and named
according to enum name. If type is Bool_t, two options "On" and "Off"
with values 0 and 1 are created. For other types you need to instruct
ROOT about possible options. The way to do it is the same as in case of
specifying getter/setter method: a comment string to a data field in
Your header file with class definition.
The most general format of this string is:
~~~ {.cpp}
*OPTIONS={GetMethod="getter";SetMethod="setter";Items=(it1="title1",it2="title2", ... ) }
~~~

While parsing this string ROOT firstly looks for command-tokens:
GetMethod, SetMethod, Items; They must be preceded by string
*OPTIONS= , enclosed by {} and separated by semicolons ";".
All command token should have a form TOKEN=VALUE.
All tokens are optional.
The names of getter and setter method must be enclosed by double-quote
marks (") .
Specifications of Items is slightly more complicated: you need to
put token ITEMS= and then enclose all options in curly brackets "()".
You separate options by comas ",".
Each option item may have one of the following forms:
~~~ {.cpp}
         IntegerValue  = "Text Label"

         EnumValue     = "Text Label"

        "TextValue" = Text Label"

~~~

One can specify values as Integers or Enums - when data field is an
Integer, Float or Enum type; as texts - for char (more precisely:
Option_t).

As mentioned above - this information are mainly used by contextmenu,
but also in Dump() and Inspect() methods and by the THtml class.
*/

#include "TDataMember.h"

#include "Strlen.h"
#include "strtok.h"
#include "strlcpy.h"
#include "TBuffer.h"
#include "TClass.h"
#include "TClassEdit.h"
#include "TDataType.h"
#include "TEnum.h"
#include "TEnumConstant.h"
#include "TGlobal.h"
#include "TInterpreter.h"
#include "TIterator.h"
#include "TList.h"
#include "TListOfDataMembers.h"
#include "TMethod.h"
#include "TMethodCall.h"
#include "TRealData.h"
#include "TROOT.h"
#include "TVirtualMutex.h"

#include <cassert>
#include <cctype>
#include <cstdlib>


ClassImp(TDataMember);

////////////////////////////////////////////////////////////////////////////////
/// Default TDataMember ctor. TDataMembers are constructed in TClass
/// via a call to TCling::CreateListOfDataMembers(). It parses the comment
/// string, initializes optionlist and getter/setter methods.

TDataMember::TDataMember(DataMemberInfo_t *info, TClass *cl) : TDictionary()
{
   fInfo        = info;
   fClass       = cl;
   fDataType    = 0;
   fOptions     = 0;
   fValueSetter = 0;
   fValueGetter = 0;
   fOffset      = -1;
   fProperty    = -1;
   fSTLCont     = -1;
   fArrayDim    = -1;
   fArrayMaxIndex=0;
   if (!fInfo && !fClass) return; // default ctor is called

   Init(false);
}

////////////////////////////////////////////////////////////////////////////////
/// Routines called by the constructor and Update to reset the member's
/// information.
/// afterReading is set when initializing after reading through Streamer().

void TDataMember::Init(bool afterReading)
{
   if (!afterReading) {
      // Initialize from fInfo
      if (!fInfo || !gInterpreter->DataMemberInfo_IsValid(fInfo))
         return;
      // Also sets names.
      Property();
   }
   const char *t = GetTitle();
   if (t && t[0] != '!')
      SetBit(kObjIsPersistent);
   fDataType = 0;
   if (IsBasic() || IsEnum()) {
      if (IsBasic()) {
         const char *name = GetFullTypeName();
         if (strcmp(name, "unsigned char") != 0 &&
             strncmp(name, "unsigned short", sizeof ("unsigned short")) != 0 &&
             strcmp(name, "unsigned int") != 0 &&
             strncmp(name, "unsigned long", sizeof ("unsigned long")) != 0)
            // strncmp() also covers "unsigned long long"
            name = GetTypeName();
         fDataType = gROOT->GetType(name);

         if (fDataType==0) {
            // humm we did not find it ... maybe it's a typedef that has not been loaded yet.
            // (this can happen if the executable does not have a TApplication object).
            fDataType = gROOT->GetType(name,kTRUE);
         }
      } else {
         fDataType = gROOT->GetType("Int_t", kTRUE); // In rare instance we are called before Int_t has been added to the list of types in TROOT, the kTRUE insures it is there.
      }
      //         if (!fDataType)
      //            Error("TDataMember", "basic data type %s not found in list of basic types",
      //                  GetTypeName());
   }


   if (afterReading) {
      // Options are streamed; can't build TMethodCall for getters and setters
      // because we deserialize a TDataMember when we do not have interpreter
      // data. Thus do an early return.
      return;
   }

   if (strstr(GetTitle(), "*OPTION={")) {
      // Delay setting fOptions until it's used: the enum constants might
      // not have been added as members yet.

   // if option string does not exist but it's an Enum - parse it!!!!
   } else if (IsEnum()) {
      fOptions = new TList();
      if (TEnum* enumDict = TEnum::GetEnum(GetTypeName()) ){
         TIter iEnumConst(enumDict->GetConstants());
         while (TEnumConstant* enumConst = (TEnumConstant*)iEnumConst()) {
            TOptionListItem *it
               = new TOptionListItem(this, enumConst->GetValue(),0,0,
                                     enumConst->GetName(),enumConst->GetName());
            fOptions->Add(it);
         }
      }

   // and the case od Bool_t : we add items "ON" and "Off"
   } else if (!strncmp(GetFullTypeName(),"Bool_t",6)){

      fOptions = new TList();
      TOptionListItem *it = new TOptionListItem(this,1,0,0,"ON",0);
      fOptions->Add(it);
      it = new TOptionListItem(this,0,0,0,"Off",0);
      fOptions->Add(it);

   } else fOptions = 0;

}

////////////////////////////////////////////////////////////////////////////////
/// copy constructor

TDataMember::TDataMember(const TDataMember& dm) :
  TDictionary(dm),
  fInfo(gCling->DataMemberInfo_FactoryCopy(dm.fInfo)),
  fClass(dm.fClass),
  fDataType(dm.fDataType),
  fOffset(dm.fOffset),
  fSTLCont(dm.fSTLCont),
  fProperty(dm.fProperty),
  fArrayDim(dm.fArrayDim),
  fArrayMaxIndex( dm.fArrayDim ? new Int_t[dm.fArrayDim] : 0),
  fArrayIndex(dm.fArrayIndex),
  fTypeName(dm.fTypeName),
  fFullTypeName(dm.fFullTypeName),
  fTrueTypeName(dm.fTrueTypeName),
  fValueGetter(0),
  fValueSetter(0),
  fOptions(dm.fOptions ? (TList*)dm.fOptions->Clone() : 0)
{
   for(Int_t d = 0; d < fArrayDim; ++d)
      fArrayMaxIndex[d] = dm.fArrayMaxIndex[d];
}

////////////////////////////////////////////////////////////////////////////////
/// assignment operator

TDataMember& TDataMember::operator=(const TDataMember& dm)
{
   if(this!=&dm) {
      gCling->DataMemberInfo_Delete(fInfo);
      delete fValueSetter;
      delete fValueGetter;
      if (fOptions) {
         fOptions->Delete();
         delete fOptions;
         fOptions = 0;
      }

      TDictionary::operator=(dm);
      fInfo= gCling->DataMemberInfo_FactoryCopy(dm.fInfo);
      fClass=dm.fClass;
      fDataType=dm.fDataType;
      fOffset=dm.fOffset;
      fSTLCont=dm.fSTLCont;
      fProperty=dm.fProperty;
      fArrayDim = dm.fArrayDim;
      fArrayMaxIndex = dm.fArrayDim ? new Int_t[dm.fArrayDim] : 0;
      for(Int_t d = 0; d < fArrayDim; ++d)
         fArrayMaxIndex[d] = dm.fArrayMaxIndex[d];
      fArrayIndex = dm.fArrayIndex;
      fTypeName=dm.fTypeName;
      fFullTypeName=dm.fFullTypeName;
      fTrueTypeName=dm.fTrueTypeName;
      fOptions = dm.fOptions ? (TList*)dm.fOptions->Clone() : 0;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// TDataMember dtor deletes adopted CINT DataMemberInfo object.

TDataMember::~TDataMember()
{
   delete [] fArrayMaxIndex;
   gCling->DataMemberInfo_Delete(fInfo);
   delete fValueSetter;
   delete fValueGetter;
   if (fOptions) {
      fOptions->Delete();
      delete fOptions;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of array dimensions.

Int_t TDataMember::GetArrayDim() const
{
   if (fArrayDim<0 && fInfo) {
      R__LOCKGUARD(gInterpreterMutex);
      TDataMember *dm = const_cast<TDataMember*>(this);
      dm->fArrayDim = gCling->DataMemberInfo_ArrayDim(fInfo);
      // fArrayMaxIndex should be zero
      if (dm->fArrayDim) {
         dm->fArrayMaxIndex = new Int_t[fArrayDim];
         for(Int_t dim = 0; dim < dm->fArrayDim; ++dim) {
            dm->fArrayMaxIndex[dim] = gCling->DataMemberInfo_MaxIndex(fInfo,dim);
         }
      }
   }
   return fArrayDim;
}

////////////////////////////////////////////////////////////////////////////////
/// If the data member is pointer and has a valid array size in its comments
/// GetArrayIndex returns a string pointing to it;
/// otherwise it returns an empty string.

const char *TDataMember::GetArrayIndex() const
{
   if (!IsaPointer()) return "";
   if (fArrayIndex.Length()==0 && fInfo) {
      R__LOCKGUARD(gInterpreterMutex);
      TDataMember *dm = const_cast<TDataMember*>(this);
      const char* val = gCling->DataMemberInfo_ValidArrayIndex(fInfo);
      if (val) dm->fArrayIndex = val;
      else dm->fArrayIndex.Append((Char_t)0); // Make length non-zero but string still empty.
   }
   return fArrayIndex;
}

////////////////////////////////////////////////////////////////////////////////

TDictionary::DeclId_t TDataMember::GetDeclId() const
{
   if (fInfo) return gInterpreter->GetDeclId(fInfo);
   else return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return maximum index for array dimension "dim".

Int_t TDataMember::GetMaxIndex(Int_t dim) const
{
   if (fArrayDim<0 && fInfo) {
      return gCling->DataMemberInfo_MaxIndex(fInfo,dim);
   } else {
      if (dim < 0 || dim >= fArrayDim) return -1;
      return fArrayMaxIndex[dim];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get type of data member, e,g.: "class TDirectory*" -> "TDirectory".

const char *TDataMember::GetTypeName() const
{
   if (fProperty==(-1)) Property();
   return fTypeName.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Get full type description of data member, e,g.: "class TDirectory*".

const char *TDataMember::GetFullTypeName() const
{
   if (fProperty==(-1)) Property();

   return fFullTypeName.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Get full type description of data member, e,g.: "class TDirectory*".

const char *TDataMember::GetTrueTypeName() const
{
   return fTrueTypeName.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Get offset from "this".

Long_t TDataMember::GetOffset() const
{
   if (fOffset>=0) return fOffset;

   R__LOCKGUARD(gInterpreterMutex);
   //case of an interpreted or emulated class
   if (fClass->GetDeclFileLine() < 0) {
      ((TDataMember*)this)->fOffset = gCling->DataMemberInfo_Offset(fInfo);
      return fOffset;
   }
   //case of a compiled class
   //Note that the offset cannot be computed in case of an abstract class
   //for which the list of real data has not yet been computed via
   //a real daughter class.
   TString dmbracket;
   dmbracket.Form("%s[",GetName());
   fClass->BuildRealData();
   TIter next(fClass->GetListOfRealData());
   TRealData *rdm;
   Int_t offset = 0;
   while ((rdm = (TRealData*)next())) {
      char *rdmc = (char*)rdm->GetName();
      //next statement required in case a class and one of its parent class
      //have data members with the same name
      if (this->IsaPointer() && rdmc[0] == '*') rdmc++;

      if (rdm->GetDataMember() != this) continue;
      if (strcmp(rdmc,GetName()) == 0) {
         offset = rdm->GetThisOffset();
         break;
      }
      if (strcmp(rdm->GetName(),GetName()) == 0) {
         if (rdm->IsObject()) {
            offset = rdm->GetThisOffset();
            break;
         }
      }
      if (strstr(rdm->GetName(),dmbracket.Data())) {
         offset = rdm->GetThisOffset();
         break;
      }
   }
   ((TDataMember*)this)->fOffset = offset;
   return fOffset;
}

////////////////////////////////////////////////////////////////////////////////
/// Get offset from "this" using the information in CINT only.

Long_t TDataMember::GetOffsetCint() const
{
   if (fOffset>=0) return fOffset;

   R__LOCKGUARD(gInterpreterMutex);
   TDataMember *dm = const_cast<TDataMember*>(this);

   if (dm->IsValid()) return gCling->DataMemberInfo_Offset(dm->fInfo);
   else return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the sizeof the underlying type of the data member
/// (i.e. if the member is an array sizeof(member)/length)

Int_t TDataMember::GetUnitSize() const
{
   if (IsaPointer()) return sizeof(void*);
   if (IsEnum()    ) return sizeof(Int_t);
   if (IsBasic()   ) return GetDataType()->Size();

   TClass *cl = TClass::GetClass(GetTypeName());
   if (!cl) cl = TClass::GetClass(GetTrueTypeName());
   if ( cl) return cl->Size();

   Warning("GetUnitSize","Can not determine sizeof(%s)",GetTypeName());
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if data member is a basic type, e.g. char, int, long...

Bool_t TDataMember::IsBasic() const
{
   if (fProperty == -1) Property();
   return (fProperty & kIsFundamental) ? kTRUE : kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if data member is an enum.

Bool_t TDataMember::IsEnum() const
{
   if (fProperty == -1) Property();
   return (fProperty & kIsEnum) ? kTRUE : kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if data member is a pointer.

Bool_t TDataMember::IsaPointer() const
{
   if (fProperty == -1) Property();
   return (fProperty & kIsPointer) ? kTRUE : kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// The return type is defined in TDictionary (kVector, kList, etc.)

int TDataMember::IsSTLContainer()
{
   if (fSTLCont != -1) return fSTLCont;
   R__LOCKGUARD(gInterpreterMutex);
   fSTLCont = TClassEdit::UnderlyingIsSTLCont(GetTrueTypeName());
   return fSTLCont;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if this data member object is pointing to a currently
/// loaded data member.  If a function is unloaded after the TDataMember
/// is created, the TDataMember will be set to be invalid.

Bool_t TDataMember::IsValid()
{
   if (fOffset >= 0) return kTRUE;

   // Register the transaction when checking the validity of the object.
   if (!fInfo && UpdateInterpreterStateMarker()) {
      DeclId_t newId = gInterpreter->GetDataMember(fClass->GetClassInfo(), fName);
      if (newId) {
         DataMemberInfo_t *info
            = gInterpreter->DataMemberInfo_Factory(newId, fClass->GetClassInfo());
         Update(info);
         // We need to make sure that the list of data member is properly
         // informed and updated.
         TListOfDataMembers *lst = dynamic_cast<TListOfDataMembers*>(fClass->GetListOfDataMembers());
         lst->Update(this);
      }
      return newId != 0;
   }
   return fInfo != 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get property description word. For meaning of bits see EProperty.

Long_t TDataMember::Property() const
{
   if (fProperty!=(-1)) return fProperty;

   R__LOCKGUARD(gInterpreterMutex);
   TDataMember *t = (TDataMember*)this;

   if (!fInfo || !gCling->DataMemberInfo_IsValid(fInfo)) return 0;
   int prop  = gCling->DataMemberInfo_Property(fInfo);
   int propt = gCling->DataMemberInfo_TypeProperty(fInfo);
   t->fProperty = prop|propt;

   t->fFullTypeName = TClassEdit::GetLong64_Name(gCling->DataMemberInfo_TypeName(fInfo));
   t->fTrueTypeName = TClassEdit::GetLong64_Name(gCling->DataMemberInfo_TypeTrueName(fInfo));
   t->fTypeName     = TClassEdit::GetLong64_Name(gCling->TypeName(fTrueTypeName));

   t->fName  = gCling->DataMemberInfo_Name(fInfo);
   t->fTitle = gCling->DataMemberInfo_Title(fInfo);

   return fProperty;
}


////////////////////////////////////////////////////////////////////////////////
/// Build TOptionListItems from the member comment `*OPTION={`

void TDataMember::ExtractOptionsFromComment()
{
   if (fOptions)
      return;

   const char *optTitle = strstr(GetTitle(), "*OPTION={");
   if (!optTitle)
      return;

   // If option string exist in comment - we'll parse it and create
   // list of options

   // Option-list string has a form:
   // *OPTION={GetMethod="GetXXX";SetMethod="SetXXX";
   //          Items=(0="NULL ITEM","one"="First Item",kRed="Red Item")}
   //
   // As one can see it is possible to specify value as either numerical
   // value , string  or enum.
   // One can also specify implicitly names of Getter/Setter methods.

   char cmt[2048];
   char opt[2048];
   const char *ptr1    = 0;
   char *ptr2    = 0;
   char *ptr3    = 0;
   char *tok     = 0;
   Int_t cnt     = 0;
   Int_t token_cnt;
   Int_t i;

   strlcpy(cmt,GetTitle(),2048);

   char *opt_ptr = strstr(cmt, "*OPTION={");

   // If we found it - parsing...

   //let's cut the part lying between {}
   char *rest;
   ptr1 = R__STRTOK_R(opt_ptr, "{}", &rest); // starts tokenizing:extracts "*OPTION={"
   if (ptr1 == 0) {
      Fatal("TDataMember","Internal error, found \"*OPTION={\" but not \"{}\" in %s.",GetTitle());
      return;
   }
   ptr1 = R__STRTOK_R(nullptr, "{}", &rest); // And now we have what we need in ptr1!!!
   if (ptr1 == 0) {
      Fatal("TDataMember","Internal error, found \"*OPTION={\" but not \"{}\" in %s.",GetTitle());
      return;
   }

   //and save it:
   strlcpy(opt,ptr1,2048);

   // Let's extract sub-tokens extracted by ';' sign.
   // We'll put'em in an array for convenience;
   // You have to do it in this manner because you cannot use nested tokenizing

   char *tokens[256];           // a storage for these sub-tokens.
   token_cnt = 0;
   cnt       = 0;

   do {                          //tokenizing loop
      ptr1 = R__STRTOK_R((char *)(cnt++ ? nullptr : opt), ";", &rest);
      if (ptr1){
         Int_t nch = strlen(ptr1)+1;
         tok=new char[nch];
         strlcpy(tok,ptr1,nch);
         tokens[token_cnt]=tok;
         token_cnt++;
      }
   } while (ptr1);

   // OK! Now let's check whether we have Get/Set methods encode in any string
   for (i=0;i<token_cnt;i++) {
      if (strstr(tokens[i],"GetMethod")) {
         ptr1 = R__STRTOK_R(tokens[i], "\"", &rest); // tokenizing-strip text "GetMethod"
         if (ptr1 == 0) {
            Fatal("TDataMember","Internal error, found \"GetMethod\" but not \"\\\"\" in %s.",GetTitle());
            return;
         }
         ptr1 = R__STRTOK_R(nullptr, "\"", &rest); // tokenizing - name is in ptr1!
         if (ptr1 == 0) {
            Fatal("TDataMember","Internal error, found \"GetMethod\" but not \"\\\"\" in %s.",GetTitle());
            return;
         }

         if (GetClass()->GetMethod(ptr1,"")) // check whether such method exists
            // FIXME: wrong in case called derives via multiple inheritance from this class
            fValueGetter = new TMethodCall(GetClass(),ptr1,"");

         continue; //next item!
      }

      if (strstr(tokens[i],"SetMethod")) {
         ptr1 = R__STRTOK_R(tokens[i], "\"", &rest);
         if (ptr1 == 0) {
            Fatal("TDataMember","Internal error, found \"SetMethod\" but not \"\\\"\" in %s.",GetTitle());
            return;
         }
         ptr1 = R__STRTOK_R(nullptr, "\"", &rest); // name of Setter in ptr1
         if (ptr1 == 0) {
            Fatal("TDataMember","Internal error, found \"SetMethod\" but not \"\\\"\" in %s.",GetTitle());
            return;
         }
         if (GetClass()->GetMethod(ptr1,"1"))
            // FIXME: wrong in case called derives via multiple inheritance from this class
            fValueSetter = new TMethodCall(GetClass(),ptr1,"1");
      }
   }

   //Now let's parse option strings...

   Int_t  opt_cnt    = 0;
   std::unique_ptr<TList> optionlist{new TList()};       //storage for options strings

   for (i=0;i<token_cnt;i++) {
      if (strstr(tokens[i],"Items")) {
         ptr1 = R__STRTOK_R(tokens[i], "()", &rest);
         if (ptr1 == 0) {
            Fatal("TDataMember","Internal error, found \"Items\" but not \"()\" in %s.",GetTitle());
            return;
         }
         ptr1 = R__STRTOK_R(nullptr, "()", &rest);
         if (ptr1 == 0) {
            Fatal("TDataMember","Internal error, found \"Items\" but not \"()\" in %s.",GetTitle());
            return;
         }

         char opts[2048];  //and save it!
         strlcpy(opts,ptr1,2048);

         //now parse it...
         //firstly we just store strings like: xxx="Label Name"
         //We'll store it in TOptionListItem objects, because they're derived
         //from TObject and thus can be stored in TList.
         //It's not elegant but works.
         do {
            ptr1 = R__STRTOK_R(opt_cnt++ ? nullptr : opts, ",", &rest); // options extraction
            if (ptr1) {
               TOptionListItem *it = new TOptionListItem(this,1,0,0,ptr1,"");
               optionlist->Add(it);
            }
         } while(ptr1);

      }
   }

   //having all options extracted and put into list, we finally can parse
   //them to create a list of options...

   fOptions = new TList();                //create the list

   TIter next(optionlist.get());                //we'll iterate through all
                                          //strings containing options
   TOptionListItem *it  = 0;
   TOptionListItem *it1 = 0;
   while ((it=(TOptionListItem*)next())) {

      ptr1 = it->fOptName;  // We will change the value of OptName ... but it is fine since we delete the object at the end of the loop.
      Bool_t islabel = (ptr1[0]=='\"');   // value is label or numerical?
      ptr2 = R__STRTOK_R((char *)ptr1, "=\"", &rest); // extract LeftHandeSide
      ptr3 = R__STRTOK_R(nullptr, "=\"", &rest);            // extract RightHandedSize

      if (islabel) {
         it1=new TOptionListItem(this,-9999,0,0,ptr3,ptr2);
         fOptions->Add(it1);
      }  else {

         char *strtolResult;
         Long_t l = std::strtol(ptr1, &strtolResult, 10);
         bool isnumber = (strtolResult != ptr1);

         if (!isnumber) {
            TGlobal *enumval = gROOT->GetGlobal(ptr1, kTRUE);
            if (enumval) {
               Int_t *value = (Int_t *)(enumval->GetAddress());
               // We'll try to find global enum existing in ROOT...
               l = (Long_t)(*value);
            } else if (IsEnum()) {
               TObject *obj = fClass->GetListOfDataMembers(false)->FindObject(ptr1);
               if (obj)
                  l = ((TEnumConstant *)obj)->GetValue();
               else
                  l = gInterpreter->Calc(Form("%s;", ptr1));
            } else {
               Fatal("TDataMember", "Internal error, couldn't recognize enum/global value %s.", ptr1);
            }
         }

         it1 = new TOptionListItem(this,l,0,0,ptr3,ptr1);
         fOptions->Add(it1);
      }

      optionlist->Remove(it);         //delete this option string from list
      delete it;                      // and dispose of it.

   }

   // Garbage collection

   //And dispose tokens string...
   for (i=0;i<token_cnt;i++)
      if(tokens[i])
         delete [] tokens[i];
}


////////////////////////////////////////////////////////////////////////////////
/// Returns list of options - list of TOptionListItems

TList *TDataMember::GetOptions()
{
   if (!fOptions)
      ExtractOptionsFromComment();
   return fOptions;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a TMethodCall method responsible for getting the value
/// of data member. The cl argument specifies the class of the object
/// which will be used to call this method (in case of multiple
/// inheritance TMethodCall needs to know this to calculate the proper
/// offset).

TMethodCall *TDataMember::GetterMethod(TClass *cl)
{
   if (!fValueGetter || cl) {

      R__LOCKGUARD(gInterpreterMutex);

      if (!cl) cl = fClass;

      if (fValueGetter) {
         TString methodname = fValueGetter->GetMethodName();
         delete fValueGetter;
         fValueGetter = new TMethodCall(cl, methodname.Data(), "");

      } else {
         // try to guess Getter function:
         // we strip the fist character of name of data field ('f') and then
         // try to find the name of Getter by applying "Get", "Is" or "Has"
         // as a prefix

         const char *dataname = GetName();

         TString gettername;
         gettername.Form( "Get%s", dataname+1);
         if (GetClass()->GetMethod(gettername, ""))
            return fValueGetter = new TMethodCall(cl, gettername, "");
         gettername.Form( "Is%s", dataname+1);
         if (GetClass()->GetMethod(gettername, ""))
            return fValueGetter = new TMethodCall(cl, gettername, "");
         gettername.Form( "Has%s", dataname+1);
         if (GetClass()->GetMethod(gettername, ""))
            return fValueGetter = new TMethodCall(cl, gettername, "");
      }
   }

   return fValueGetter;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a TMethodCall method responsible for setting the value
/// of data member. The cl argument specifies the class of the object
/// which will be used to call this method (in case of multiple
/// inheritance TMethodCall needs to know this to calculate the proper
/// offset).

TMethodCall *TDataMember::SetterMethod(TClass *cl)
{
   if (!fValueSetter || cl) {

      R__LOCKGUARD(gInterpreterMutex);

      if (!cl) cl = fClass;

      if (fValueSetter) {

         TString methodname = fValueSetter->GetMethodName();
         TString params = fValueSetter->GetParams();
         delete fValueSetter;
         fValueSetter = new TMethodCall(cl, methodname.Data(), params.Data());

      } else {

         // try to guess Setter function:
         // we strip the fist character of name of data field ('f') and then
         // try to find the name of Setter by applying "Set" as a prefix

         const char *dataname = GetName();

         TString settername;
         settername.Form( "Set%s", dataname+1);
         if (strstr(settername, "Is")) settername.Form( "Set%s", dataname+3);
         if (GetClass()->GetMethod(settername, "1"))
            fValueSetter = new TMethodCall(cl, settername, "1");
         if (!fValueSetter)
            if (GetClass()->GetMethod(settername, "true"))
               fValueSetter = new TMethodCall(cl, settername, "true");
      }
   }

   return fValueSetter;
}

////////////////////////////////////////////////////////////////////////////////
/// Update the TFunction to reflect the new info.
///
/// This can be used to implement unloading (info == 0) and then reloading
/// (info being the 'new' decl address).

Bool_t TDataMember::Update(DataMemberInfo_t *info)
{
   R__LOCKGUARD(gInterpreterMutex);

   if (fInfo) gCling->DataMemberInfo_Delete(fInfo);
   SafeDelete(fValueSetter);
   SafeDelete(fValueGetter);
   if (fOptions) {
      fOptions->Delete();
      SafeDelete(fOptions);
   }

   if (info == 0) {
      fOffset      = -1;
      fProperty    = -1;
      fSTLCont     = -1;
      fArrayDim    = -1;
      delete [] fArrayMaxIndex;
      fArrayMaxIndex=0;
      fArrayIndex.Clear();

      fInfo = 0;
      return kTRUE;
   } else {
      fInfo = info;
      Init(false);
      return kTRUE;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of TDataMember. Forces calculation of all cached
/// (and persistent) values.

void TDataMember::Streamer(TBuffer& b) {
   if (b.IsReading()) {
      b.ReadClassBuffer(Class(), this);
      Init(true /*reading*/);
   } else {
      // Writing.
      if (fProperty & kIsStatic) {
         // We have a static member and in this case fOffset contains the
         // actual address in memory of the data, it will be different everytime,
         // let's not record it.
         fOffset = -1;
      } else {
         GetOffset();
      }
      IsSTLContainer();
      GetArrayDim();
      GetArrayIndex();
      Property(); // also calculates fTypeName and friends
      b.WriteClassBuffer(Class(), this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TOptionListItem::TOptionListItem(TDataMember *d, Long_t val, Long_t valmask,
                 Long_t tglmask,const char *name, const char *label)
{
   fDataMember    = d;
   fValue         = val;
   fValueMaskBit  = valmask;
   fToggleMaskBit = tglmask;
   if (name) {
      fOptName = name;
   }

   if (label) {
      fOptLabel = label;
   }
}
