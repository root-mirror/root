// @(#)root/meta:$Id$
// Author: Axel Naumann 2014-04-28

/*************************************************************************
 * Copyright (C) 1995-2014, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef ROOT_TProtoClass
#define ROOT_TProtoClass

#include "TNamed.h"

#include <vector>

class TClass;
class TList;
class TRealData;
#ifdef R__LESS_INCLUDES
class TDataMember;
#else
#include "TDataMember.h"
#endif

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TProtoClass                                                          //
//                                                                      //
// Stores enough information to create a TClass from a dictionary       //
// without interpreter information.                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TProtoClass: public TNamed {
public:
   struct TProtoRealData  {
      Long_t fOffset;     // data member offset
      Int_t  fDMIndex;    // index of data member in vector of data members
      Int_t  fLevel;      // member level (0 : belong to this class, 1 is a data member of a data member object, etc...)
      Int_t  fClassIndex; // index of class belonging to in list of dep classes
      char   fStatusFlag; // status of the real data member (if bit 0 set is an object, if bit 1 set is transient if bit 2 set is a pointer)

      enum EStatusFlags {
         kIsObject    = BIT(0),    // member is object
         kIsTransient = BIT(1),    // data member is transient
         kIsPointer   = BIT(2),    // data member is a pointer
         kBitMask     = 0x000000ff
      };

   public:
      bool IsAClass() const { return fClassIndex >= 0; }
      TProtoRealData() : fOffset(0), fDMIndex(-1), fLevel(0), fClassIndex(-1), fStatusFlag(0) {}
      TProtoRealData(const TRealData *rd);
      virtual ~TProtoRealData();
      TRealData *CreateRealData(TClass *currentClass, TClass *parent, TRealData * parentData, int prevLevel) const;

      Bool_t TestFlag(UInt_t f) const { return (Bool_t) ((fStatusFlag & f) != 0); }
      void SetFlag(UInt_t f, Bool_t on = kTRUE) {
         if (on)
            fStatusFlag |= f & kBitMask;
         else
            fStatusFlag &= ~(f & kBitMask);
      }

      ClassDef(TProtoRealData, 3);//Persistent version of TRealData
   };

private:
   TList   *fBase;     // List of base classes
   TList   *fEnums;    // List of enums in this scope
   std::vector<TProtoRealData> fPRealData;  // List of TProtoRealData
   std::vector<TDataMember *>  fData;       // collection of data members
   std::vector<TString>        fDepClasses; // list of dependent classes
   Int_t    fSizeof;         // Size of the class
   UInt_t   fCheckSum;       //checksum of data members and base classes
   Int_t    fCanSplit;       // Whether this class can be split
   Int_t    fStreamerType;   // Which streaming method to use
   Long_t   fProperty;       // Class properties, see EProperties
   Long_t   fClassProperty;  // Class C++ properties, see EClassProperties
   Long_t   fOffsetStreamer; // Offset to streamer function

   TProtoClass(const TProtoClass &) = delete;
   TProtoClass &operator=(const TProtoClass &) = delete;

   const char * GetClassName(Int_t index) const { return (index >= 0) ? fDepClasses[index].Data() : 0; }

   // compute index of data member in the list
   static Int_t DataMemberIndex(TClass * cl, const char * name);
   // find data member  given an index
   static TDataMember * FindDataMember(TClass * cl,  Int_t index);

public:
   TProtoClass():
      fBase(0), fEnums(0), fSizeof(0), fCheckSum(0), fCanSplit(0),
      fStreamerType(0), fProperty(0), fClassProperty(0),
      fOffsetStreamer(0) {
   }

   TProtoClass(TProtoClass *pc);
   TProtoClass(TClass *cl);
   virtual ~TProtoClass();

   Bool_t FillTClass(TClass *pcl);
   const TList *GetListOfEnums() {
      return fEnums;
   };
   void Delete(Option_t *opt = "");

   int GetSize() { return fSizeof; }
   TList * GetBaseList() { return fBase; }
   //TList * GetDataList() { return fData; }
   TList * GetEnumList() { return fEnums; }
   std::vector<TProtoRealData> & GetPRDList() { return fPRealData; }
   std::vector<TDataMember *> & GetData() { return fData; }
   std::vector<TString> & GetDepClasses() { return fDepClasses; }


   ClassDef(TProtoClass, 2); //Persistent TClass
};

#endif
