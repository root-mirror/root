// @(#)root/treeplayer:$Id$
// Author: Philippe Canal 01/06/2004

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers and al.        *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFormLeafInfo
#define ROOT_TFormLeafInfo

#include "TObject.h"

#include "TLeafElement.h"

#include "TArrayI.h"
#include "TDataType.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"


// declare the extra versions of GetValue() plus templated implementation
#define DECLARE_GETVAL \
   virtual Double_t  GetValue(TLeaf *leaf, Int_t instance = 0)               \
       { return GetValueImpl<Double_t>(leaf, instance); }                    \
   virtual Long64_t  GetValueLong64(TLeaf *leaf, Int_t instance = 0)         \
       { return GetValueImpl<Long64_t>(leaf, instance); }                    \
   virtual LongDouble_t  GetValueLongDouble(TLeaf *leaf, Int_t instance = 0) \
       { return GetValueImpl<LongDouble_t>(leaf, instance); }                \
   template<typename T> T  GetValueImpl(TLeaf *leaf, Int_t instance = 0)   // no semicolon


// declare the extra versions of ReadValue() plus templated implementation
#define DECLARE_READVAL \
   virtual Double_t ReadValue(char *where, Int_t instance = 0)               \
       { return ReadValueImpl<Double_t>(where, instance); }                  \
   virtual Long64_t ReadValueLong64(char *where, Int_t instance = 0)         \
       { return ReadValueImpl<Long64_t>(where, instance); }                  \
   virtual LongDouble_t ReadValueLongDouble(char *where, Int_t instance = 0) \
       { return ReadValueImpl<LongDouble_t>(where, instance); }              \
   template<typename T> T  ReadValueImpl(char *where, Int_t instance = 0)  // no semicolon




class TFormLeafInfo : public TObject {
public:
   // Constructors
   TFormLeafInfo(TClass* classptr = 0, Longptr_t offset = 0,
                 TStreamerElement* element = 0);
   TFormLeafInfo(const TFormLeafInfo& orig);
   virtual TFormLeafInfo* DeepCopy() const;
   virtual ~TFormLeafInfo();

   void Swap(TFormLeafInfo &other);
   TFormLeafInfo &operator=(const TFormLeafInfo &orig);

   // Data Members
   TClass           *fClass;   ///<! This is the class of the data pointed to
   //TStreamerInfo  *fInfo;    ///<! == fClass->GetStreamerInfo()
   Longptr_t         fOffset;  ///<! Offset of the data pointed inside the class fClass
   TStreamerElement *fElement; ///<! Descriptor of the data pointed to.
         //Warning, the offset in fElement is NOT correct because it does not take into
         //account base classes and nested objects (which fOffset does).
   TFormLeafInfo    *fCounter;
   TFormLeafInfo    *fNext;    ///< follow this to grab the inside information
   TString fClassName;
   TString fElementName;

protected:
   Int_t fMultiplicity;
public:

   virtual void AddOffset(Int_t offset, TStreamerElement* element);

   virtual Int_t GetArrayLength();
   virtual TClass*   GetClass() const;
   virtual Int_t     GetCounterValue(TLeaf* leaf);
   virtual Int_t     ReadCounterValue(char *where);

   char* GetObjectAddress(TLeafElement* leaf, Int_t &instance);

   Int_t GetMultiplicity();

   // Currently only implemented in TFormLeafInfoCast
   Int_t GetNdata(TLeaf* leaf);
   virtual Int_t GetNdata();

   virtual void     *GetValuePointer(TLeaf *leaf, Int_t instance = 0);
   virtual void     *GetValuePointer(char  *from, Int_t instance = 0);
   virtual void     *GetLocalValuePointer(TLeaf *leaf, Int_t instance = 0);
   virtual void     *GetLocalValuePointer( char *from, Int_t instance = 0);

   virtual Bool_t    HasCounter() const;
   virtual Bool_t    IsString() const;

   virtual Bool_t    IsInteger() const;
   virtual Bool_t    IsReference() const  {  return kFALSE; }

   // Method for multiple variable dimensions.
   virtual Int_t GetPrimaryIndex();
   virtual Int_t GetVarDim();
   virtual Int_t GetVirtVarDim();
   virtual Int_t GetSize(Int_t index);
   virtual Int_t GetSumOfSizes();
   virtual void  LoadSizes(TBranch* branch);
   virtual void  SetPrimaryIndex(Int_t index);
   virtual void  SetSecondaryIndex(Int_t index);
   virtual void  SetSize(Int_t index, Int_t val);
   virtual void  SetBranch(TBranch* br)  { if ( fNext ) fNext->SetBranch(br); }
   virtual void  UpdateSizes(TArrayI *garr);

   virtual Bool_t    Update();

   DECLARE_GETVAL;
   DECLARE_READVAL;

   template <typename T> struct ReadValueHelper {
      static T Exec(TFormLeafInfo *leaf, char *where, Int_t instance) {
         return leaf->ReadValue(where, instance);
      }
   };
   template <typename T > T ReadTypedValue(char *where, Int_t instance = 0) {
      return ReadValueHelper<T>::Exec(this, where, instance);
   }

   template <typename T> struct GetValueHelper {
      static T Exec(TFormLeafInfo *linfo, TLeaf *leaf, Int_t instance) {
         return linfo->GetValue(leaf, instance);
      }
   };
   template <typename T > T GetTypedValue(TLeaf *leaf, Int_t instance = 0) {
      return GetValueHelper<T>::Exec(this, leaf, instance);
   }
};


template <> struct TFormLeafInfo::ReadValueHelper<Long64_t> {
   static Long64_t Exec(TFormLeafInfo *leaf, char *where, Int_t instance) { return leaf->ReadValueLong64(where, instance); }
};
template <> struct TFormLeafInfo::ReadValueHelper<ULong64_t> {
   static ULong64_t Exec(TFormLeafInfo *leaf, char *where, Int_t instance) { return (ULong64_t)leaf->ReadValueLong64(where, instance); }
};
template <> struct TFormLeafInfo::ReadValueHelper<LongDouble_t> {
   static LongDouble_t Exec(TFormLeafInfo *leaf, char *where, Int_t instance) { return leaf->ReadValueLongDouble(where, instance); }
};

template <> struct TFormLeafInfo::GetValueHelper<Long64_t> {
   static Long64_t Exec(TFormLeafInfo *linfo, TLeaf *leaf, Int_t instance) { return linfo->GetValueLong64(leaf, instance); }
};
template <> struct TFormLeafInfo::GetValueHelper<ULong64_t> {
   static ULong64_t Exec(TFormLeafInfo *linfo, TLeaf *leaf, Int_t instance) { return (ULong64_t)linfo->GetValueLong64(leaf, instance); }
};
template <> struct TFormLeafInfo::GetValueHelper<LongDouble_t> {
   static LongDouble_t Exec(TFormLeafInfo *linfo, TLeaf *leaf, Int_t instance) { return linfo->GetValueLongDouble(leaf, instance); }
};

// TFormLeafInfoDirect is a small helper class to implement reading a data
// member on an object stored in a TTree.

class TFormLeafInfoDirect : public TFormLeafInfo {
public:
   TFormLeafInfoDirect(TBranchElement * from);
   // The implicit default constructor's implementation is correct.

   virtual TFormLeafInfo* DeepCopy() const;

   DECLARE_GETVAL;
   virtual void     *GetLocalValuePointer(TLeaf *leaf, Int_t instance = 0);
   virtual void     *GetLocalValuePointer(char *thisobj, Int_t instance = 0);

   virtual Double_t  ReadValue(char * /*where*/, Int_t /*instance*/= 0);
   virtual Long64_t  ReadValueLong64(char *where, Int_t i= 0) { return ReadValue(where, i); }
   virtual LongDouble_t  ReadValueLongDouble(char *where, Int_t i= 0) { return ReadValue(where, i); }

};

// TFormLeafInfoNumerical is a small helper class to implement reading a
// numerical value inside a collection

class TFormLeafInfoNumerical : public TFormLeafInfo {
   EDataType fKind;
   Bool_t    fIsBool;
public:
   TFormLeafInfoNumerical(TVirtualCollectionProxy *holder_of);
   TFormLeafInfoNumerical(EDataType kind);
   TFormLeafInfoNumerical(const TFormLeafInfoNumerical& orig);

   virtual TFormLeafInfo* DeepCopy() const;
   void Swap(TFormLeafInfoNumerical &other);
   TFormLeafInfoNumerical &operator=(const TFormLeafInfoNumerical &orig);

   virtual ~TFormLeafInfoNumerical();

   virtual Bool_t    IsString() const;
   virtual Bool_t    Update();
};

// TFormLeafInfoCollectionObject
// This class is used when we are interested by the collection it self and
// it is split.

class TFormLeafInfoCollectionObject : public TFormLeafInfo {
   Bool_t fTop;  //If true, it indicates that the branch itself contains
public:
   TFormLeafInfoCollectionObject(TClass* classptr = 0, Bool_t fTop = kTRUE);
   TFormLeafInfoCollectionObject(const TFormLeafInfoCollectionObject &orig);

   void Swap(TFormLeafInfoCollectionObject &other);
   TFormLeafInfoCollectionObject &operator=(const TFormLeafInfoCollectionObject &orig);

   virtual TFormLeafInfo* DeepCopy() const {
      return new TFormLeafInfoCollectionObject(*this);
   }

   DECLARE_GETVAL;
   virtual Int_t     GetCounterValue(TLeaf* leaf);
   virtual Double_t  ReadValue(char *where, Int_t instance = 0);
   virtual Long64_t  ReadValueLong64(char *where, Int_t i= 0) { return ReadValue(where, i); }
   virtual LongDouble_t  ReadValueLongDouble(char *where, Int_t i= 0) { return ReadValue(where, i); }
   virtual void     *GetValuePointer(TLeaf *leaf, Int_t instance = 0);
   virtual void     *GetValuePointer(char  *thisobj, Int_t instance = 0);
   virtual void     *GetLocalValuePointer(TLeaf *leaf, Int_t instance = 0);
   virtual void     *GetLocalValuePointer(char  *thisobj, Int_t instance = 0);
};

// TFormLeafInfoClones is a small helper class to implement reading a data
// member on a TClonesArray object stored in a TTree.

class TFormLeafInfoClones : public TFormLeafInfo {
   Bool_t fTop;  //If true, it indicates that the branch itself contains
public:
   //either the clonesArrays or something inside the clonesArray
   TFormLeafInfoClones(TClass* classptr = 0, Longptr_t offset = 0);
   TFormLeafInfoClones(TClass* classptr, Longptr_t offset, Bool_t top);
   TFormLeafInfoClones(TClass* classptr, Longptr_t offset, TStreamerElement* element,
                       Bool_t top = kFALSE);
   TFormLeafInfoClones(const TFormLeafInfoClones &orig);

   void Swap(TFormLeafInfoClones &other);
   TFormLeafInfoClones &operator=(const TFormLeafInfoClones &orig);

   virtual TFormLeafInfo* DeepCopy() const {
      return new TFormLeafInfoClones(*this);
   }

   DECLARE_GETVAL;
   DECLARE_READVAL;
   virtual Int_t     GetCounterValue(TLeaf* leaf);
   virtual Int_t     ReadCounterValue(char *where);
   virtual void     *GetValuePointer(TLeaf *leaf, Int_t instance = 0);
   virtual void     *GetValuePointer(char  *thisobj, Int_t instance = 0);
   virtual void     *GetLocalValuePointer(TLeaf *leaf, Int_t instance = 0);
   virtual void     *GetLocalValuePointer(char  *thisobj, Int_t instance = 0);
};

// TFormLeafInfoCollection is a small helper class to implement reading a data member
// on a generic collection object stored in a TTree.

class TFormLeafInfoCollection : public TFormLeafInfo {
   Bool_t fTop;  //If true, it indicates that the branch itself contains
   //either the clonesArrays or something inside the clonesArray
   TClass                  *fCollClass;
   TString                  fCollClassName;
   TVirtualCollectionProxy *fCollProxy;
   TStreamerElement        *fLocalElement;
public:

   TFormLeafInfoCollection(TClass* classptr,
                           Longptr_t offset,
                           TStreamerElement* element,
                           Bool_t top = kFALSE);

   TFormLeafInfoCollection(TClass* motherclassptr,
                           Longptr_t offset = 0,
                           TClass* elementclassptr = 0,
                           Bool_t top = kFALSE);

   TFormLeafInfoCollection();
   TFormLeafInfoCollection(const TFormLeafInfoCollection& orig);

   ~TFormLeafInfoCollection();

   void Swap(TFormLeafInfoCollection &other);
   TFormLeafInfoCollection &operator=(const TFormLeafInfoCollection &orig);

   virtual TFormLeafInfo* DeepCopy() const;

   virtual Bool_t    Update();

   DECLARE_GETVAL;
   DECLARE_READVAL;
   virtual Int_t     GetCounterValue(TLeaf* leaf);
   virtual Int_t     ReadCounterValue(char* where);
   virtual Int_t     GetCounterValue(TLeaf* leaf, Int_t instance);
   virtual Bool_t    HasCounter() const;
   virtual void     *GetValuePointer(TLeaf *leaf, Int_t instance = 0);
   virtual void     *GetValuePointer(char  *thisobj, Int_t instance = 0);
   virtual void     *GetLocalValuePointer(TLeaf *leaf, Int_t instance = 0);
   virtual void     *GetLocalValuePointer(char  *thisobj, Int_t instance = 0);
};

// TFormLeafInfoCollectionSize is used to return the size of a collection

class TFormLeafInfoCollectionSize : public TFormLeafInfo {
   TClass                  *fCollClass;
   TString                  fCollClassName;
   TVirtualCollectionProxy *fCollProxy;
public:
   TFormLeafInfoCollectionSize(TClass*);
   TFormLeafInfoCollectionSize(TClass* classptr,Longptr_t offset,TStreamerElement* element);
   TFormLeafInfoCollectionSize();
   TFormLeafInfoCollectionSize(const TFormLeafInfoCollectionSize& orig);

   ~TFormLeafInfoCollectionSize();

   void Swap(TFormLeafInfoCollectionSize &other);
   TFormLeafInfoCollectionSize &operator=(const TFormLeafInfoCollectionSize &orig);

   virtual TFormLeafInfo* DeepCopy() const;

   virtual Bool_t    Update();

   virtual void     *GetValuePointer(TLeaf *leaf, Int_t instance = 0);
   virtual void     *GetValuePointer(char  *from, Int_t instance = 0);
   virtual void     *GetLocalValuePointer(TLeaf *leaf, Int_t instance = 0);
   virtual void     *GetLocalValuePointer( char *from, Int_t instance = 0);
   virtual Double_t  ReadValue(char *where, Int_t instance = 0);
   virtual Long64_t  ReadValueLong64(char *where, Int_t i= 0) { return ReadValue(where, i); }
   virtual LongDouble_t  ReadValueLongDouble(char *where, Int_t i= 0) { return ReadValue(where, i); }
};

// TFormLeafInfoPointer is a small helper class to implement reading a data
// member by following a pointer inside a branch of TTree.

class TFormLeafInfoPointer : public TFormLeafInfo {
public:
   TFormLeafInfoPointer(TClass* classptr = 0, Longptr_t offset = 0,
                        TStreamerElement* element = 0);
   // The default copy constructor is the right implementation.

   virtual TFormLeafInfo* DeepCopy() const;

   DECLARE_GETVAL;
   DECLARE_READVAL;
};

// TFormLeafInfoMethod is a small helper class to implement executing a method
// of an object stored in a TTree

class TFormLeafInfoMethod : public TFormLeafInfo {

   TMethodCall *fMethod;
   TString      fMethodName;
   TString      fParams;
   Double_t     fResult;
   TString      fCopyFormat;
   TString      fDeleteFormat;
   void        *fValuePointer;
   Bool_t       fIsByValue;

public:
   static TClass *ReturnTClass(TMethodCall *mc);

   TFormLeafInfoMethod(TClass* classptr = 0, TMethodCall *method = 0);
   TFormLeafInfoMethod(const TFormLeafInfoMethod& orig);
   ~TFormLeafInfoMethod();

   void Swap(TFormLeafInfoMethod &other);
   TFormLeafInfoMethod &operator=(const TFormLeafInfoMethod &orig);

   virtual TFormLeafInfo* DeepCopy() const;

   DECLARE_READVAL;
   virtual TClass*  GetClass() const;
   virtual void    *GetLocalValuePointer( TLeaf *from, Int_t instance = 0);
   virtual void    *GetLocalValuePointer(char *from, Int_t instance = 0);
   virtual Bool_t   IsInteger() const;
   virtual Bool_t   IsString() const;
   virtual Bool_t   Update();
};

// TFormLeafInfoMultiVarDim is a small helper class to implement reading a
// data member on a variable size array inside a TClonesArray object stored in
// a TTree.  This is the version used when the data member is inside a
// non-split object.

class TFormLeafInfoMultiVarDim : public TFormLeafInfo {
public:
   Int_t fNsize;
   TArrayI fSizes;           // Array of sizes of the variable dimension
   TFormLeafInfo *fCounter2; // Information on how to read the secondary dimensions
   Int_t fSumOfSizes;        // Sum of the content of fSizes
   Int_t fDim;               // physical number of the dimension that is variable
   Int_t fVirtDim;           // number of the virtual dimension to which this object correspond.
   Int_t fPrimaryIndex;      // Index of the dimensions that is indexing the second dimension's size
   Int_t fSecondaryIndex;    // Index of the second dimension

protected:
   TFormLeafInfoMultiVarDim(TClass* classptr, Longptr_t offset,
                            TStreamerElement* element) : TFormLeafInfo(classptr,offset,element),fNsize(0),fSizes(),fCounter2(0),fSumOfSizes(0),fDim(0),fVirtDim(0),fPrimaryIndex(-1),fSecondaryIndex(-1) {}

public:
   TFormLeafInfoMultiVarDim(TClass* classptr, Longptr_t offset,
                            TStreamerElement* element, TFormLeafInfo* parent);
   TFormLeafInfoMultiVarDim();
   TFormLeafInfoMultiVarDim(const TFormLeafInfoMultiVarDim& orig);
   ~TFormLeafInfoMultiVarDim();

   void Swap(TFormLeafInfoMultiVarDim &other);
   TFormLeafInfoMultiVarDim &operator=(const TFormLeafInfoMultiVarDim &orig);

   virtual TFormLeafInfo* DeepCopy() const;

   /* The proper indexing and unwinding of index is done by prior leafinfo in the chain. */
   //virtual Double_t  ReadValue(char *where, Int_t instance = 0) {
   //   return TFormLeafInfo::ReadValue(where,instance);
   //}

   virtual void     LoadSizes(TBranch* branch);
   virtual Int_t    GetPrimaryIndex();
   virtual void     SetPrimaryIndex(Int_t index);
   virtual void     SetSecondaryIndex(Int_t index);
   virtual void     SetSize(Int_t index, Int_t val);
   virtual Int_t    GetSize(Int_t index);
   virtual Int_t    GetSumOfSizes();
   virtual Double_t GetValue(TLeaf * /*leaf*/, Int_t /*instance*/ = 0);
   virtual Long64_t  GetValueLong64(TLeaf *leaf, Int_t i= 0) { return GetValue(leaf, i); }
   virtual LongDouble_t  GetValueLongDouble(TLeaf *leaf, Int_t i= 0) { return GetValue(leaf, i); }
   virtual Int_t    GetVarDim();
   virtual Int_t    GetVirtVarDim();
   virtual Bool_t   Update();
   virtual void     UpdateSizes(TArrayI *garr);
};

// TFormLeafInfoMultiVarDimDirect is a small helper class to implement reading
// a data member on a variable size array inside a TClonesArray object stored
// in a TTree.  This is the version used for split access

class TFormLeafInfoMultiVarDimDirect : public TFormLeafInfoMultiVarDim {
public:
   // The default constructor are the correct implementation.

   virtual TFormLeafInfo* DeepCopy() const;

   DECLARE_GETVAL;
   virtual Double_t  ReadValue(char * /*where*/, Int_t /*instance*/ = 0);
   virtual Long64_t  ReadValueLong64(char *where, Int_t i= 0) { return ReadValue(where, i); }
   virtual LongDouble_t  ReadValueLongDouble(char *where, Int_t i= 0) { return ReadValue(where, i); }
};

// TFormLeafInfoMultiVarDimCollection is a small helper class to implement reading
// a data member which is a collection inside a TClonesArray or collection object
// stored in a TTree.  This is the version used for split access
//
class TFormLeafInfoMultiVarDimCollection : public TFormLeafInfoMultiVarDim {
public:
   TFormLeafInfoMultiVarDimCollection(TClass* motherclassptr, Longptr_t offset,
      TClass* elementclassptr, TFormLeafInfo *parent);
   TFormLeafInfoMultiVarDimCollection(TClass* classptr, Longptr_t offset,
      TStreamerElement* element, TFormLeafInfo* parent);
   // The default copy constructor is the right implementation.

   virtual TFormLeafInfo* DeepCopy() const;

   virtual Int_t GetArrayLength() { return 0; }
   virtual void      LoadSizes(TBranch* branch);
   virtual Double_t  GetValue(TLeaf *leaf, Int_t instance = 0);
   virtual Long64_t  GetValueLong64(TLeaf *leaf, Int_t i= 0) { return GetValue(leaf, i); }
   virtual LongDouble_t  GetValueLongDouble(TLeaf *leaf, Int_t i= 0) { return GetValue(leaf, i); }
   DECLARE_READVAL;
};

// TFormLeafInfoMultiVarDimClones is a small helper class to implement reading
// a data member which is a TClonesArray inside a TClonesArray or collection object
// stored in a TTree.  This is the version used for split access
//
class TFormLeafInfoMultiVarDimClones : public TFormLeafInfoMultiVarDim {
public:
   TFormLeafInfoMultiVarDimClones(TClass* motherclassptr, Longptr_t offset,
      TClass* elementclassptr, TFormLeafInfo *parent);
   TFormLeafInfoMultiVarDimClones(TClass* classptr, Longptr_t offset,
      TStreamerElement* element, TFormLeafInfo* parent);
   // The default copy constructor is the right implementation.

   virtual TFormLeafInfo* DeepCopy() const;

   virtual Int_t GetArrayLength() { return 0; }
   virtual void      LoadSizes(TBranch* branch);
   virtual Double_t  GetValue(TLeaf *leaf, Int_t instance = 0);
   virtual Long64_t  GetValueLong64(TLeaf *leaf, Int_t i= 0) { return GetValue(leaf, i); }
   virtual LongDouble_t  GetValueLongDouble(TLeaf *leaf, Int_t i= 0) { return GetValue(leaf, i); }
   DECLARE_READVAL;
};

// TFormLeafInfoCast is a small helper class to implement casting an object to
// a different type (equivalent to dynamic_cast)

class TFormLeafInfoCast : public TFormLeafInfo {
public:
   TClass *fCasted;     //! Pointer to the class we are trying to case to
   TString fCastedName; //! Name of the class we are casting to.
   Bool_t  fGoodCast;   //! Marked by ReadValue.
   Bool_t  fIsTObject;  //! Indicated whether the fClass inherits from TObject.

   TFormLeafInfoCast(TClass* classptr = 0, TClass* casted = 0);
   TFormLeafInfoCast(const TFormLeafInfoCast& orig);
   virtual ~TFormLeafInfoCast();

   void Swap(TFormLeafInfoCast &other);
   TFormLeafInfoCast &operator=(const TFormLeafInfoCast &orig);

   virtual TFormLeafInfo* DeepCopy() const;

   DECLARE_READVAL;
   // Currently only implemented in TFormLeafInfoCast
   virtual Int_t GetNdata();
   virtual Bool_t    Update();
};

// TFormLeafInfoTTree is a small helper class to implement reading
// from the containing TTree object itself.

class TFormLeafInfoTTree : public TFormLeafInfo {
   TTree   *fTree;
   TTree   *fCurrent;
   TString  fAlias;

public:
   TFormLeafInfoTTree(TTree *tree = 0, const char *alias = 0, TTree *current = 0);
   TFormLeafInfoTTree(const TFormLeafInfoTTree& orig);

   void Swap(TFormLeafInfoTTree &other);
   TFormLeafInfoTTree &operator=(const TFormLeafInfoTTree &orig);

   virtual TFormLeafInfo* DeepCopy() const;

   using TFormLeafInfo::GetLocalValuePointer;
   using TFormLeafInfo::GetValue;

   DECLARE_GETVAL;
   DECLARE_READVAL;
   virtual void     *GetLocalValuePointer(TLeaf *leaf, Int_t instance = 0);
   virtual Bool_t    Update();
};


#endif /* ROOT_TFormLeafInfo */

