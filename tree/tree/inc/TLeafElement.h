// @(#)root/tree:$Id$
// Author: Rene Brun   14/01/2001

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TLeafElement
#define ROOT_TLeafElement


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TLeafElement                                                          //
//                                                                      //
// A TLeaf for a general object derived from TObject.                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <atomic>

#include "TDataType.h"
#include "TLeaf.h"
#include "TBranchElement.h"

class TMethodCall;

class TLeafElement : public TLeaf {

protected:
   char               *fAbsAddress;   ///<! Absolute leaf Address
   Int_t               fID;           ///<  element serial number in fInfo
   Int_t               fType;         ///<  leaf type
   mutable std::atomic<DeserializeType> fDeserializeTypeCache{ DeserializeType::kInvalid }; ///<! Cache of the type of deserialization.
   mutable std::atomic<EDataType> fDataTypeCache{EDataType::kOther_t}; ///<! Cache of the EDataType of deserialization.

private:
   virtual Int_t       GetOffsetHeaderSize() const {return 1;}

public:
   TLeafElement();
   TLeafElement(TBranch *parent, const char *name, Int_t id, Int_t type);
   virtual ~TLeafElement();

   virtual Bool_t   CanGenerateOffsetArray() { return fLeafCount && fLenType; }
   virtual Int_t   *GenerateOffsetArrayBase(Int_t /*base*/, Int_t /*events*/) { return nullptr; }
   virtual DeserializeType GetDeserializeType() const;

           Int_t    GetID() const { return fID; }
   virtual TString  GetFullName() const;
   virtual Int_t    GetLen() const {return ((TBranchElement*)fBranch)->GetNdata()*fLen;}
   TMethodCall     *GetMethodCall(const char *name);
   virtual Int_t    GetMaximum() const {return ((TBranchElement*)fBranch)->GetMaximum();}
   virtual Int_t    GetNdata() const {return ((TBranchElement*)fBranch)->GetNdata()*fLen;}
   virtual const char *GetTypeName() const {return ((TBranchElement*)fBranch)->GetTypeName();}

   virtual Double_t     GetValue(Int_t i=0) const { return ((TBranchElement*)fBranch)->GetValue(i, fLen, kFALSE);}
   virtual Long64_t     GetValueLong64(Int_t i = 0) const { return ((TBranchElement*)fBranch)->GetTypedValue<Long64_t>(i, fLen, kFALSE); }
   virtual LongDouble_t GetValueLongDouble(Int_t i = 0) const { return ((TBranchElement*)fBranch)->GetTypedValue<LongDouble_t>(i, fLen, kFALSE); }
   template<typename T> T GetTypedValueSubArray(Int_t i=0, Int_t j=0) const {return ((TBranchElement*)fBranch)->GetTypedValue<T>(i, j, kTRUE);}

   virtual bool     ReadBasketFast(TBuffer&, Long64_t);
 
   virtual void    *GetValuePointer() const { return ((TBranchElement*)fBranch)->GetValuePointer(); }
   virtual Bool_t   IncludeRange(TLeaf *);
   virtual Bool_t   IsOnTerminalBranch() const;
   virtual void     PrintValue(Int_t i=0) const {((TBranchElement*)fBranch)->PrintValue(i);}
   virtual void     SetLeafCount(TLeaf *leaf) {fLeafCount = leaf;}

   ClassDef(TLeafElement,1);  //A TLeaf for a general object derived from TObject.
};

#endif
