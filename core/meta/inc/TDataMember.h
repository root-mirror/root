// @(#)root/meta:$Id$
// Author: Fons Rademakers   04/02/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TDataMember
#define ROOT_TDataMember


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TDataMember                                                          //
//                                                                      //
// Dictionary interface for a class data member.                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TDictionary.h"

class TList;
class TClass;
class TDataType;
class TMethodCall;

class TDataMember : public TDictionary {

private:
   enum EStatusBits {
      kObjIsPersistent = BIT(2)
   };

   DataMemberInfo_t   *fInfo;         //!pointer to CINT data member info
   TClass             *fClass;        //!pointer to the class
   TDataType          *fDataType;     //!pointer to data basic type descriptor

   Long_t              fOffset;       //offset
   Int_t               fSTLCont;      //STL type
   Long_t              fProperty;     //Property
   Int_t               fArrayDim;     //Number of array dimensions
   Int_t              *fArrayMaxIndex;//[fArrayDim] Maximum index for each dimension
   TString             fArrayIndex;   //String representation of the index variable name

   TString             fTypeName;     //data member type, e,g.: "class TDirectory*" -> "TDirectory".
   TString             fFullTypeName; //full type description of data member, e,g.: "class TDirectory*".
   TString             fTrueTypeName; //full type description with no typedef

   // The following fields allows to access all (even private) datamembers and
   // provide a possibility of having options with names and strings.
   // These options are defined in a comment to a field!
   TMethodCall        *fValueGetter;  //!method that returns a value;
   TMethodCall        *fValueSetter;  //!method which sets value;
   TList              *fOptions;      //list of possible values 0=no restrictions

   void Init(bool afterReading);

   void ExtractOptionsFromComment();

protected:
   TDataMember(const TDataMember&);
   TDataMember& operator=(const TDataMember&);

public:

   TDataMember(DataMemberInfo_t *info = 0, TClass *cl = 0);
   virtual       ~TDataMember();
   Int_t          GetArrayDim() const;
   DeclId_t       GetDeclId() const;
   Int_t          GetMaxIndex(Int_t dim) const;
   TClass        *GetClass() const { return fClass; }
   TDataType     *GetDataType() const { return fDataType; } //only for basic type
   Long_t         GetOffset() const;
   Long_t         GetOffsetCint() const;
   const char    *GetTypeName() const;
   const char    *GetFullTypeName() const;
   const char    *GetTrueTypeName() const;
   const char    *GetArrayIndex() const;
   Int_t          GetUnitSize() const;
   TList         *GetOptions();
   TMethodCall   *SetterMethod(TClass *cl);
   TMethodCall   *GetterMethod(TClass *cl = 0);

   Bool_t         IsBasic() const;
   Bool_t         IsEnum() const;
   Bool_t         IsaPointer() const;
   Bool_t         IsPersistent() const { return TestBit(kObjIsPersistent); }
   Int_t          IsSTLContainer();
   Bool_t         IsValid();
   Long_t         Property() const;
   void           SetClass(TClass* cl) { fClass = cl; }
   virtual bool   Update(DataMemberInfo_t *info);

   ClassDef(TDataMember,2)  //Dictionary for a class data member
};


// This class implements one option in options list. All Data members are public
// for convenience reasons.

class TOptionListItem : public TObject {

public:
   TDataMember     *fDataMember;     //!Data member to which this option belongs
   Long_t           fValue;          //Numerical value assigned to option
   Long_t           fValueMaskBit;   //Not used yet: bitmask used when option is a toggle group
   Long_t           fToggleMaskBit;  //Not used yet: bitmask used when toggling value
   TString          fOptName;        //Text assigned to option which appears in option menu
   TString          fOptLabel;       //Text (or enum) value assigned to option.
   TOptionListItem():
      fDataMember(0), fValue(0), fValueMaskBit(0), fToggleMaskBit(0)
   {}
   TOptionListItem(TDataMember *m,Long_t val, Long_t valmask, Long_t tglmask,
                   const char *name, const char *label);

   ClassDef(TOptionListItem,2); //Element in the list of options.
};

#endif
