// @(#)root/tree:$Id$
// Author: Rene Brun   12/01/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TLeaf
#define ROOT_TLeaf


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TLeaf                                                                //
//                                                                      //
// A TTree object is a list of TBranch.                                 //
// A TBranch object is a list of TLeaf.  In most cases, the TBranch     //
// will have one TLeaf.                                                 //
// A TLeaf describes the branch data types and holds the data.          //
//                                                                      //
// A few notes about the data held by the leaf.  It can contain:        //
//   1 a single object or primitive (e.g., one float),                  //
//   2 a fixed-number of objects (e.g., each entry has two floats).     //
//     The number of elements per entry is saved in `fLen`.             //
//   3 a dynamic number of primitives.  The number of objects in each   //
//     entry is saved in the `fLeafCount` branch.                       //
//                                                                      //
// Note options (2) and (3) can combined - if fLeafCount says an entry  //
// has 3 elements and fLen is 2, then there will be 6 objects in that   //
// entry.                                                               //
//                                                                      //
// Additionally, `fNdata` is transient and generated on read to         //
// determine the necessary size of a buffer to hold event data;         //
// depending on the call-site, it may be sized larger than the number   //
// of elements                                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TNamed.h"

#include <vector>

#ifdef R__LESS_INCLUDES
class TBranch;
#else
#include "TBranch.h"
#endif

class TClonesArray;
class TBrowser;

class TLeaf : public TNamed {

private:

   virtual Int_t GetOffsetHeaderSize() const {return 0;}

protected:

   using Counts_t = std::vector<Int_t>;
   struct LeafCountValues {
      Counts_t fValues;
      Long64_t fStartEntry{-1}; ///<! entry number of corresponding to element 0 of the vector.
   };

   Int_t            fNdata;           ///<! Number of elements in fAddress data buffer.
   Int_t            fLen;             ///<  Number of fixed length elements in the leaf's data.
   Int_t            fLenType;         ///<  Number of bytes for this data type
   Int_t            fOffset;          ///<  Offset in ClonesArray object (if one)
   Bool_t           fIsRange;         ///<  (=kTRUE if leaf has a range, kFALSE otherwise).  This is equivalent to being a 'leafcount'.  For a TLeafElement the range information is actually store in the TBranchElement.
   Bool_t           fIsUnsigned;      ///<  (=kTRUE if unsigned, kFALSE otherwise)
   TLeaf           *fLeafCount;       ///<  Pointer to Leaf count if variable length (we do not own the counter)
   TBranch         *fBranch;          ///<! Pointer to supporting branch (we do not own the branch)
   LeafCountValues *fLeafCountValues; ///<! Cache of collection/array sizes

   TLeaf(const TLeaf&);
   TLeaf& operator=(const TLeaf&);

  template <typename T> struct GetValueHelper {
    static T Exec(const TLeaf *leaf, Int_t i = 0) { return leaf->GetValue(i); }
  };

  Int_t *GenerateOffsetArrayBase(Int_t base, Int_t events) const; // For leaves containing fixed-size objects (no
                                                                  // polymorphism!), this will generate an appropriate
                                                                  // offset array.


public:
   enum EStatusBits {
      kIndirectAddress = BIT(11), ///< Data member is a pointer to an array of basic types.
      kNewValue = BIT(12)         ///< Set if we own the value buffer and so must delete it ourselves.
   };

   enum class DeserializeType {
      kInvalid = 0,      // Invalid deserialization information.
      kDestructive,      // Deserialization of this Leaf requires a separate output buffer.
      kInPlace,          // Deserialization can be done directly in the input buffer.
      kZeroCopy,         // In-memory and on-disk representation of this object are identical.
   };

   TLeaf();
   TLeaf(TBranch *parent, const char *name, const char *type);
   virtual ~TLeaf();

   virtual void     Browse(TBrowser *b);
   virtual Bool_t   CanGenerateOffsetArray() {return fLeafCount;} // overload and return true if this leaf can generate its own offset array.
   virtual void     Export(TClonesArray *, Int_t) {}
   virtual void     FillBasket(TBuffer &b);
   virtual Int_t   *GenerateOffsetArray(Int_t base, Int_t events) { return GenerateOffsetArrayBase(base, events); }
   TBranch         *GetBranch() const { return fBranch; }
   virtual DeserializeType GetDeserializeType() const { return DeserializeType::kDestructive; }
   virtual TString  GetFullName() const;
   ///  If this leaf stores a variable-sized array or a multi-dimensional array whose last dimension has variable size,
   ///  return a pointer to the TLeaf that stores such size. Return a nullptr otherwise.
   virtual TLeaf   *GetLeafCount() const { return fLeafCount; }
   virtual TLeaf   *GetLeafCounter(Int_t &countval) const;

   virtual const Counts_t *GetLeafCountValues(Long64_t start, Long64_t len);

   virtual Int_t    GetLen() const;
   /// Return the fixed length of this leaf.
   /// If the leaf stores a fixed-length array, this is the size of the array.
   /// If the leaf stores a non-array or a variable-sized array, this method returns 1.
   /// If the leaf stores an array with 2 or more dimensions, this method returns the total number of elements in the
   /// dimensions with static length: for example for float[3][2][] it would return 6.
   virtual Int_t    GetLenStatic() const { return fLen; }
   virtual Int_t    GetLenType() const { return fLenType; }
   virtual Int_t    GetMaximum() const { return 0; }
   virtual Int_t    GetMinimum() const { return 0; }
   virtual Int_t    GetNdata() const { return fNdata; }
   virtual Int_t    GetOffset() const { return fOffset; }
   virtual void    *GetValuePointer() const { return 0; }
   virtual const char *GetTypeName() const { return ""; }

   virtual Double_t GetValue(Int_t i = 0) const;
   virtual Long64_t GetValueLong64(Int_t i = 0) const { return GetValue(i); }         // overload only when it matters.
   virtual LongDouble_t GetValueLongDouble(Int_t i = 0) const { return GetValue(i); } // overload only when it matters.
   template <typename T> T GetTypedValue(Int_t i = 0) const { return GetValueHelper<T>::Exec(this, i); }

   virtual Bool_t   IncludeRange(TLeaf *) { return kFALSE; } // overload to copy/set fMinimum and fMaximum to include/be wide than those of the parameter
   virtual void     Import(TClonesArray *, Int_t) {}
   virtual Bool_t   IsOnTerminalBranch() const { return kTRUE; }
   virtual Bool_t   IsRange() const { return fIsRange; }
   virtual Bool_t   IsUnsigned() const { return fIsUnsigned; }
   virtual void     PrintValue(Int_t i = 0) const;
   virtual void     ReadBasket(TBuffer &) {}
   virtual void     ReadBasketExport(TBuffer &, TClonesArray *, Int_t) {}
   virtual bool     ReadBasketFast(TBuffer&, Long64_t) { return false; }  // Read contents of leaf into a user-provided buffer.
   virtual bool     ReadBasketSerialized(TBuffer&, Long64_t) { return false; }  // Read contents of leaf into a user-provided buffer
   virtual void     ReadValue(std::istream & /*s*/, Char_t /*delim*/ = ' ') {
      Error("ReadValue", "Not implemented!");
   }
           Int_t    ResetAddress(void *add, Bool_t destructor = kFALSE);
   virtual void     SetAddress(void *add = 0);
   virtual void     SetBranch(TBranch *branch) { fBranch = branch; }
   virtual void     SetLeafCount(TLeaf *leaf);
   virtual void     SetLen(Int_t len = 1) { fLen = len; }
   virtual void     SetOffset(Int_t offset = 0) { fOffset = offset; }
   virtual void     SetRange(Bool_t range = kTRUE) { fIsRange = range; }
   virtual void     SetUnsigned() { fIsUnsigned = kTRUE; }

   ClassDef(TLeaf, 2); // Leaf: description of a Branch data type
};


template <> struct TLeaf::GetValueHelper<Long64_t> {
   static Long64_t Exec(const TLeaf *leaf, Int_t i = 0) { return leaf->GetValueLong64(i); }
};
template <> struct TLeaf::GetValueHelper<ULong64_t> {
   static ULong64_t Exec(const TLeaf *leaf, Int_t i = 0) { return (ULong64_t)leaf->GetValueLong64(i); }
};
template <> struct TLeaf::GetValueHelper<LongDouble_t> {
   static LongDouble_t Exec(const TLeaf *leaf, Int_t i = 0) { return leaf->GetValueLongDouble(i); }
};


inline Double_t TLeaf::GetValue(Int_t /*i = 0*/) const { return 0.0; }
inline void     TLeaf::PrintValue(Int_t /* i = 0*/) const {}
inline void     TLeaf::SetAddress(void* /* add = 0 */) {}

#endif
