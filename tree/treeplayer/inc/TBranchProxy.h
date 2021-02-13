// @(#)root/treeplayer:$Id$
// Author: Philippe Canal 01/06/2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers and al.        *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBranchProxy
#define ROOT_TBranchProxy

#include "TBranchProxyDirector.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TError.h"
#include "TVirtualCollectionProxy.h"
#include "TNotifyLink.h"

#include <algorithm>
#include <string>
#include <iostream>

class TBranch;
class TStreamerElement;

// Note we could protect the arrays more by introducing a class TArrayWrapper<class T> which somehow knows
// its internal dimensions and check for them ...
// template <class T> TArrayWrapper {
// public:
//    TArrayWrapper(void *where, int dim1);
//    const T operator[](int i) {
//       if (i>=dim1) return 0;
//       return where[i];
//    };
// };
// 2D array would actually be a wrapper of a wrapper i.e. has a method TArrayWrapper<T> operator[](int i);

namespace ROOT {
namespace Internal {
   //_______________________________________________
   // String builder to be used in the constructors.
   class TBranchProxyHelper {
   public:
      TString fName;
      TBranchProxyHelper(const char *left,const char *right = 0) :
         fName() {
         if (left) {
            fName = left;
            if (left[0]&&right && fName[fName.Length()-1]!='.') fName += ".";
         }
         if (right) {
            fName += right;
         }
      }
      operator const char*() { return fName.Data(); };
   };

   class TTreeReaderValueBase;
} // namespace Internal


namespace Detail {
   class TBranchProxy {
   protected:
      Internal::TBranchProxyDirector *fDirector; // contain pointer to TTree and entry to be read

      Bool_t        fInitialized : 1;
      const Bool_t  fIsMember : 1;    // true if we proxy an unsplit data member
      Bool_t        fIsClone : 1;     // true if we proxy the inside of a TClonesArray
      Bool_t        fIsaPointer : 1;  // true if we proxy a data member of pointer type
      Bool_t        fHasLeafCount : 1;// true if we proxy a variable size leaf of a leaflist

      const TString fBranchName;  // name of the branch to read
      TBranchProxy *fParent;      // Proxy to a parent object

      const TString fDataMember;  // name of the (eventual) data member being proxied


      TString           fClassName;     // class name of the object pointed to by the branch
      TClass           *fClass;         // class name of the object pointed to by the branch
      TStreamerElement *fElement;
      Int_t             fMemberOffset;
      Int_t             fOffset;        // Offset inside the object
      Int_t             fArrayLength;   // Number of element if the data is an array

      TBranch *fBranch;       // branch to read
      union {
         TBranch *fBranchCount;  // eventual auxiliary branch (for example holding the size)
         TLeaf   *fLeafCount;    // eventual auxiliary leaf (for example holding the size)
      };

      TNotifyLink<TBranchProxy> fNotify; // Callback object used by the TChain to update this proxy

      Long64_t fRead;     // Last entry read

      void    *fWhere;    // memory location of the data
      TVirtualCollectionProxy *fCollection; // Handle to the collection containing the data chunk.

   public:
      virtual void Print();

      TBranchProxy();
      TBranchProxy(Internal::TBranchProxyDirector* boss, const char* top, const char* name = 0);
      TBranchProxy(Internal::TBranchProxyDirector* boss, const char *top, const char *name, const char *membername);
      TBranchProxy(Internal::TBranchProxyDirector* boss, TBranchProxy *parent, const char* membername, const char* top = 0, const char* name = 0);
      TBranchProxy(Internal::TBranchProxyDirector* boss, TBranch* branch, const char* membername);
      TBranchProxy(Internal::TBranchProxyDirector* boss, const char* branchname, TBranch* branch, const char* membername);
      virtual ~TBranchProxy();

      TBranchProxy* GetProxy() { return this; }
      const char* GetBranchName() const { return fBranchName; }

      void Reset();

      Bool_t Notify() {
         fRead = -1;
         return Setup();
      }

      Bool_t Setup();

      Bool_t IsInitialized() {
         return fInitialized;
         // return fLastTree && fCurrentTreeNumber == fDirector->GetTree()->GetTreeNumber() && fLastTree == fDirector->GetTree();
      }

      Bool_t IsaPointer() const {
         return fIsaPointer;
      }

      Bool_t Read() {
         if (R__unlikely(fDirector==0)) return false;

         auto treeEntry = fDirector->GetReadEntry();
         if (treeEntry != fRead) {
            if (!IsInitialized()) {
               if (!Setup()) {
                  ::Error("TBranchProxy::Read","%s",Form("Unable to initialize %s\n",fBranchName.Data()));
                  return kFALSE;
               }
            }
            Bool_t result = kTRUE;
            if (fParent) {
               result = fParent->Read();
            } else {
               if (fBranchCount) {
                  result &= (-1 != fBranchCount->GetEntry(treeEntry));
               }
               result &= (-1 != fBranch->GetEntry(treeEntry));
            }
            fRead = treeEntry;
            if (R__unlikely(fCollection)) {
               fCollection->PopProxy(); // works even if no proxy env object was set.
               if (IsaPointer()) {
                  fCollection->PushProxy( *(void**)fWhere );
               } else {
                  fCollection->PushProxy( fWhere );
               }
            }
            return result;
         } else {
            return IsInitialized();
         }
      }

private:
      friend class ROOT::Internal::TTreeReaderValueBase;

      enum class EReadType {
         kDefault,
         kNoDirector,
         kReadParentNoCollection,
         kReadParentCollectionNoPointer,
         kReadParentCollectionPointer,
         kReadNoParentNoBranchCountCollectionPointer,
         kReadNoParentNoBranchCountCollectionNoPointer,
         kReadNoParentNoBranchCountNoCollection,
         kReadNoParentBranchCountCollectionPointer,
         kReadNoParentBranchCountCollectionNoPointer,
         kReadNoParentBranchCountNoCollection
      };

      EReadType GetReadType() {
         if (fParent) {
            if (!fCollection) {
               return EReadType::kReadParentNoCollection;
            } else {
               if (IsaPointer()) {
                  return EReadType::kReadParentCollectionPointer;
               } else {
                  return EReadType::kReadParentCollectionNoPointer;
               }
            }
         } else {
            if (fBranchCount) {
               if (fCollection) {
                  if (IsaPointer()) {
                     return EReadType::kReadNoParentBranchCountCollectionPointer;
                  } else {
                     return EReadType::kReadNoParentBranchCountCollectionNoPointer;
                  }
               } else {
                  return EReadType::kReadNoParentBranchCountNoCollection;
               }

            } else {
               if (fCollection) {
                  if (IsaPointer()) {
                     return EReadType::kReadNoParentNoBranchCountCollectionPointer;
                  } else {
                     return EReadType::kReadNoParentNoBranchCountCollectionNoPointer;
                  }
               } else {
                  return EReadType::kReadNoParentNoBranchCountNoCollection;
               }
            }
         }
         return EReadType::kDefault;
      }

      Bool_t ReadNoDirector() {
         return false;
      }

      Bool_t ReadParentNoCollection() {
         auto treeEntry = fDirector->GetReadEntry();
         if (treeEntry != fRead) {
            const Bool_t result = fParent->Read();
            fRead = treeEntry;
            return result;
         } else {
            return IsInitialized();
         }
      }

      Bool_t ReadParentCollectionNoPointer() {
         auto treeEntry = fDirector->GetReadEntry();
         if (treeEntry != fRead) {
            const Bool_t result = fParent->Read();
            fRead = treeEntry;
            fCollection->PopProxy(); // works even if no proxy env object was set.
            fCollection->PushProxy( fWhere );
            return result;
         } else {
            return IsInitialized();
         }
      }

      Bool_t ReadParentCollectionPointer() {
         auto treeEntry = fDirector->GetReadEntry();
         if (treeEntry != fRead) {
            const Bool_t result = fParent->Read();
            fRead = treeEntry;
            fCollection->PopProxy(); // works even if no proxy env object was set.
            fCollection->PushProxy( *(void**)fWhere );
            return result;
         } else {
            return IsInitialized();
         }
      }

      Bool_t ReadNoParentNoBranchCountCollectionPointer() {
         auto treeEntry = fDirector->GetReadEntry();
         if (treeEntry != fRead) {
            Bool_t result = (-1 != fBranch->GetEntry(treeEntry));
            fRead = treeEntry;
            fCollection->PopProxy(); // works even if no proxy env object was set.
            fCollection->PushProxy( *(void**)fWhere );
            return result;
         } else {
            return IsInitialized();
         }
      }

      Bool_t ReadNoParentNoBranchCountCollectionNoPointer() {
         auto treeEntry = fDirector->GetReadEntry();
         if (treeEntry != fRead) {
            Bool_t result = (-1 != fBranch->GetEntry(treeEntry));
            fRead = treeEntry;
            fCollection->PopProxy(); // works even if no proxy env object was set.
            fCollection->PushProxy( fWhere );
            return result;
         } else {
            return IsInitialized();
         }
      }

      Bool_t ReadNoParentNoBranchCountNoCollection() {
         auto treeEntry = fDirector->GetReadEntry();
         if (treeEntry != fRead) {
            Bool_t result = (-1 != fBranch->GetEntry(treeEntry));
            fRead = treeEntry;
            return result;
         } else {
            return IsInitialized();
         }
      }

      Bool_t ReadNoParentBranchCountCollectionPointer() {
         auto treeEntry = fDirector->GetReadEntry();
         if (treeEntry != fRead) {
            Bool_t result = (-1 != fBranchCount->GetEntry(treeEntry));
            result &= (-1 != fBranch->GetEntry(treeEntry));
            fRead = treeEntry;
            fCollection->PopProxy(); // works even if no proxy env object was set.
            fCollection->PushProxy( *(void**)fWhere );
            return result;
         } else {
            return IsInitialized();
         }
      }

      Bool_t ReadNoParentBranchCountCollectionNoPointer() {
         auto treeEntry = fDirector->GetReadEntry();
         if (treeEntry != fRead) {
            Bool_t result = (-1 != fBranchCount->GetEntry(treeEntry));
            result &= (-1 != fBranch->GetEntry(treeEntry));
            fRead = treeEntry;
            fCollection->PopProxy(); // works even if no proxy env object was set.
            fCollection->PushProxy( fWhere );
            return result;
         } else {
            return IsInitialized();
         }
      }

      Bool_t ReadNoParentBranchCountNoCollection() {
         auto treeEntry = fDirector->GetReadEntry();
         if (treeEntry != fRead) {
            Bool_t result = (-1 != fBranchCount->GetEntry(treeEntry));
            result &= (-1 != fBranch->GetEntry(treeEntry));
            fRead = treeEntry;
            return result;
         } else {
            return IsInitialized();
         }
      }

public:

      Bool_t ReadEntries() {
         if (R__unlikely(fDirector==0)) return false;

         auto treeEntry = fDirector->GetReadEntry();
         if (treeEntry != fRead) {
            if (!IsInitialized()) {
               if (!Setup()) {
                  ::Error("TBranchProxy::ReadEntries","%s",Form("Unable to initialize %s\n",fBranchName.Data()));
                  return false;
               }
            }
            if (fParent) fParent->ReadEntries();
            else {
               if (fBranchCount) {
                  fBranchCount->TBranch::GetEntry(treeEntry);
               }
               fBranch->TBranch::GetEntry(treeEntry);
            }
            // NO - we only read the entries, not the contained objects!
            // fRead = treeEntry;
         }
         return IsInitialized();
      }

      virtual Int_t GetEntries() {
         if (!ReadEntries()) return 0;
         if (!fHasLeafCount) {
            return *(Int_t*)fLeafCount->GetValuePointer();
         } else {
            return 1;
         }
      }

      virtual Int_t GetArrayLength() {
         return fArrayLength;
      }

      TClass *GetClass() {
         if (fDirector==0) return 0;

         if (fDirector->GetReadEntry() != fRead) {
            if (!IsInitialized()) {
               if (!Setup()) {
                  return 0;
               }
            }
         }
         return fClass;
      }

      void* GetWhere() const { return fWhere; } // intentionally non-virtual

      /// Return the address of the element number i. Returns `nullptr` for non-collections. It assumed that Setip() has
      /// been called.
      virtual void *GetAddressOfElement(UInt_t /*i*/) {
         return nullptr;
      }

      TVirtualCollectionProxy *GetCollection() { return fCollection; }

      // protected:
      virtual  void *GetStart(UInt_t /*i*/=0) {
         // return the address of the start of the object being proxied. Assumes
         // that Setup() has been called.

         if (fParent) {
            fWhere = ((unsigned char*)fParent->GetStart()) + fMemberOffset;
         }
         if (IsaPointer()) {
            if (fWhere) return *(void**)fWhere;
            else return 0;
         } else {
            return fWhere;
         }
      }

      void *GetClaStart(UInt_t i=0) {
         // return the address of the start of the object being proxied. Assumes
         // that Setup() has been called.  Assumes the object containing this data
         // member is held in TClonesArray.

         char *location;

         if (fIsClone) {

            TClonesArray *tca;
            tca = (TClonesArray*)GetStart();

            if (!tca || tca->GetLast()<(Int_t)i) return 0;

            location = (char*)tca->At(i);

            return location;

         } else if (fParent) {

            //tcaloc = ((unsigned char*)fParent->GetStart());
            location = (char*)fParent->GetClaStart(i);

         } else {

            void *tcaloc;
            tcaloc = fWhere;
            TClonesArray *tca;
            tca = (TClonesArray*)tcaloc;

            if (tca->GetLast()<(Int_t)i) return 0;

            location = (char*)tca->At(i);
         }

         if (location) location += fOffset;
         else return 0;

         if (IsaPointer()) {
            return *(void**)(location);
         } else {
            return location;
         }

      }

      void *GetStlStart(UInt_t i=0) {
         // return the address of the start of the object being proxied. Assumes
         // that Setup() has been called.  Assumes the object containing this data
         // member is held in STL Collection.

         char *location=0;

         if (fCollection) {

            if (fCollection->Size()<i) return 0;

            location = (char*)fCollection->At(i);

            // return location;

         } else if (fParent) {

            //tcaloc = ((unsigned char*)fParent->GetStart());
            location = (char*)fParent->GetStlStart(i);

         } else {

            R__ASSERT(0);
            //void *tcaloc;
            //tcaloc = fWhere;
            //TClonesArray *tca;
            //tca = (TClonesArray*)tcaloc;

            //if (tca->GetLast()<i) return 0;

            //location = (char*)tca->At(i);
         }

         if (location) location += fOffset;
         else return 0;

         if (IsaPointer()) {
            return *(void**)(location);
         } else {
            return location;
         }

      }

      Int_t GetOffset() { return fOffset; }
   };
} // namespace Detail

namespace Internal {

   //____________________________________________________________________________________________
   // Concrete Implementation of the branch proxy around the data members which are array of char
   class TArrayCharProxy : public Detail::TBranchProxy {
   public:
      void Print() override {
         TBranchProxy::Print();
         std::cout << "fWhere " << fWhere << std::endl;
         if (fWhere) std::cout << "value? " << *(unsigned char*)GetStart() << std::endl;
      }

      using TBranchProxy::TBranchProxy;
      TArrayCharProxy() = default; // work around bug in GCC < 7
      ~TArrayCharProxy() override = default;

      void *GetAddressOfElement(UInt_t i) final {
         if (!Read()) return nullptr;
         unsigned char* str = (unsigned char*)GetStart();
         return str + i;
      }

      unsigned char At(UInt_t i) {
         static unsigned char default_val = {};
         if (unsigned char* elAddr = (unsigned char*)GetAddressOfElement(i)) {
            // should add out-of bound test
            return *elAddr;
         }
         return default_val;
      }

      unsigned char operator [](Int_t i) {
         return At(i);
      }

      unsigned char operator [](UInt_t i) {
         return At(i);
      }

      operator const char*() {
         if (!Read()) return "";
         return (const char*)GetStart();
      }

      const char* Data() {
         if (!Read()) return "";
         return (const char*)GetStart();
      }

      const char* c_str() {
         if (!Read()) return "";
         return (const char*)GetStart();
      }

      operator std::string() {
         if (!Read()) return "";
         return std::string((const char*)GetStart());
      }

   };

   //_______________________________________________________
   // Base class for the proxy around object in TClonesArray.
   class TClaProxy : public Detail::TBranchProxy {
   public:
      void Print() override {
         TBranchProxy::Print();
         std::cout << "fWhere " << fWhere << std::endl;
         if (fWhere) {
            if (IsaPointer()) {
               std::cout << "location " << *(TClonesArray**)fWhere << std::endl;
            } else {
               std::cout << "location " << fWhere << std::endl;
            }
         }
      }

      using TBranchProxy::TBranchProxy;
      TClaProxy() = default; // work around bug in GCC < 7
      ~TClaProxy() override = default;

      const TClonesArray* GetPtr() {
         if (!Read()) return 0;
         return (TClonesArray*)GetStart();
      }

      Int_t GetEntries() override {
         if (!ReadEntries()) return 0;
         TClonesArray *arr = (TClonesArray*)GetStart();
         if (arr) return arr->GetEntries();
         return 0;
      }

      void *GetAddressOfElement(UInt_t i) final {
         if (!Read()) return nullptr;
         if (fWhere==0) return nullptr;
         return GetClaStart(i);
      }

      const TClonesArray* operator->() { return GetPtr(); }

   };

   //_______________________________________________
   // Base class for the proxy around STL containers.
   class TStlProxy : public Detail::TBranchProxy {
   public:
      void Print() override {
         TBranchProxy::Print();
         std::cout << "fWhere " << fWhere << std::endl;
         if (fWhere) {
            if (IsaPointer()) {
               std::cout << "location " << *(TClonesArray**)fWhere << std::endl;
            } else {
               std::cout << "location " << fWhere << std::endl;
            }
         }
      }

      using TBranchProxy::TBranchProxy;
      TStlProxy() = default; // work around bug in GCC < 7
      ~TStlProxy() override = default;

      TVirtualCollectionProxy* GetPtr() {
         if (!Read()) return 0;
         return GetCollection();
      }

      Int_t GetEntries() override {
         if (!ReadEntries()) return 0;
         return GetPtr()->Size();
      }

      void *GetAddressOfElement(UInt_t i) final {
         if (!Read()) return nullptr;
         if (fWhere==0) return nullptr;
         return GetStlStart(i);
      }

      const TVirtualCollectionProxy* operator->() { return GetPtr(); }

   };

   //______________________________________
   // Template of the proxy around objects.
   template <class T>
   class TImpProxy : public Detail::TBranchProxy {
   public:
      void Print() override {
         TBranchProxy::Print();
         std::cout << "fWhere " << fWhere << std::endl;
         if (fWhere) std::cout << "value? " << *(T*)GetStart() << std::endl;
      }

      using TBranchProxy::TBranchProxy;
      TImpProxy() = default; // work around bug in GCC < 7
      ~TImpProxy() override = default;

      operator T() {
         if (!Read()) return 0;
         return *(T*)GetStart();
      }

      // For now explicitly disable copying into the value (i.e. the proxy is read-only).
      TImpProxy(T) = delete;
      TImpProxy &operator=(T) = delete;

   };

   //____________________________________________
   // Helper template to be able to determine and
   // use array dimentsions.
   template <class T, int d = 0> struct TArrayType {
      typedef T type_t;
      typedef T array_t[d];
      static constexpr int gSize = d;
   };
   //____________________________________________
   // Helper class for proxy around multi dimension array
   template <class T> struct TArrayType<T,0> {
      typedef T type_t;
      typedef T array_t;
      static constexpr int gSize = 0;
   };
   //____________________________________________
   // Helper class for proxy around multi dimension array
   template <class T, int d> struct TMultiArrayType {
      typedef typename T::type_t type_t;
      typedef typename T::array_t array_t[d];
      static constexpr int gSize = d;
   };

   //____________________________________________
   // Template for concrete implementation of proxy around array of T
   template <class T>
   class TArrayProxy : public Detail::TBranchProxy {
   public:
      using TBranchProxy::TBranchProxy;
      TArrayProxy() = default; // work around bug in GCC < 7
      ~TArrayProxy() override = default;

      typedef typename T::array_t array_t;
      typedef typename T::type_t type_t;

      void Print() override {
         TBranchProxy::Print();
         std::cout << "fWhere " << GetWhere() << std::endl;
         if (GetWhere()) std::cout << "value? " << *(type_t*)GetWhere() << std::endl;
      }

      Int_t GetEntries() override {
         return T::gSize;
      }

      void *GetAddressOfElement(UInt_t i) final {
         if (!Read()) return nullptr;
         if (array_t *arr = (array_t*)((type_t*)(GetStart())))
            return &arr[i];
         return nullptr;
      }

      const array_t &At(UInt_t i) {
         static array_t default_val;
         // should add out-of bound test
         if (array_t *arr = (array_t*)GetAddressOfElement(i))
            return *arr;
         return default_val;
      }

      const array_t &operator [](Int_t i) { return At(i); }
      const array_t &operator [](UInt_t i) { return At(i); }
   };

   //_____________________________________________________________________________________
   // Template of the Concrete Implementation of the branch proxy around TClonesArray of T
   template <class T>
   class TClaImpProxy : public TClaProxy {
   public:

      // void Print() override {
      //    TClaProxy::Print();
      // }

      using TClaProxy::TClaProxy;
      TClaImpProxy() = default; // work around bug in GCC < 7
      ~TClaImpProxy() override = default;

      const T& At(UInt_t i) {
         static T default_val;
         if (void* addr = GetAddressOfElement(i))
            return *(T*)addr;
         return default_val;
      }

      const T& operator [](Int_t i) { return At(i); }
      const T& operator [](UInt_t i) { return At(i); }

      // For now explicitly disable copying into the value (i.e. the proxy is read-only).
      TClaImpProxy(T) = delete;
      TClaImpProxy &operator=(T) = delete;

   };

   //_________________________________________________________________________________________
   // Template of the Concrete Implementation of the branch proxy around an stl container of T
   template <class T>
   class TStlImpProxy : public TStlProxy {
   public:

      // void Print() override {
      //    TBranchProxy::Print();
      // }

      using TStlProxy::TStlProxy;
      TStlImpProxy() = default; // work around bug in GCC < 7
      ~TStlImpProxy() override = default;

      const T& At(UInt_t i) {
         static T default_val;
         if (void* addr = GetAddressOfElement(i))
            return *(T*)addr;
         return default_val;
      }

      const T& operator [](Int_t i) { return At(i); }
      const T& operator [](UInt_t i) { return At(i); }

      // For now explicitly disable copying into the value (i.e. the proxy is read-only).
      TStlImpProxy(T) = delete;
      TStlImpProxy &operator=(T) = delete;

   };

   //_________________________________________________________________________________________________
   // Template of the Concrete Implementation of the branch proxy around an TClonesArray of array of T
   template <class T>
   class TClaArrayProxy : public TClaProxy {
   public:
      typedef typename T::array_t array_t;
      typedef typename T::type_t type_t;

      // void Print() override {
      //    TClaProxy::Print();
      // }

      using TClaProxy::TClaProxy;
      TClaArrayProxy() = default; // work around bug in GCC < 7
      ~TClaArrayProxy() override = default;

      /* const */  array_t *At(UInt_t i) {
         static array_t default_val;
         if (array_t* ptr = (array_t*)GetAddressOfElement(i))
            return ptr; // no de-ref!

         return &default_val;
      }

      /* const */ array_t *operator [](Int_t i) { return At(i); }
      /* const */ array_t *operator [](UInt_t i) { return At(i); }
   };


   //__________________________________________________________________________________________________
   // Template of the Concrete Implementation of the branch proxy around an stl container of array of T
   template <class T>
   class TStlArrayProxy : public TStlProxy {
   public:
      typedef typename T::array_t array_t;
      typedef typename T::type_t type_t;

      // void Print() override {
      //    TBranchProxy::Print();
      // }

      using TStlProxy::TStlProxy;
      TStlArrayProxy() = default; // work around bug in GCC < 7
      ~TStlArrayProxy() override = default;

      /* const */  array_t *At(UInt_t i) {
         static array_t default_val;
         if (array_t* ptr = (array_t*)GetAddressOfElement(i))
            return ptr; // no de-ref!
         return &default_val;
      }

      /* const */ array_t *operator [](Int_t i) { return At(i); }
      /* const */ array_t *operator [](UInt_t i) { return At(i); }
   };

   //TImpProxy<TObject> d;
   typedef TImpProxy<Double_t>   TDoubleProxy;   // Concrete Implementation of the branch proxy around the data members which are double
   typedef TImpProxy<Double32_t> TDouble32Proxy; // Concrete Implementation of the branch proxy around the data members which are double32
   typedef TImpProxy<Float_t>    TFloatProxy;    // Concrete Implementation of the branch proxy around the data members which are float
   typedef TImpProxy<Float16_t>  TFloat16Proxy;  // Concrete Implementation of the branch proxy around the data members which are float16
   typedef TImpProxy<UInt_t>     TUIntProxy;     // Concrete Implementation of the branch proxy around the data members which are unsigned int
   typedef TImpProxy<ULong_t>    TULongProxy;    // Concrete Implementation of the branch proxy around the data members which are unsigned long
   typedef TImpProxy<ULong64_t>  TULong64Proxy;  // Concrete Implementation of the branch proxy around the data members which are unsigned long long
   typedef TImpProxy<UShort_t>   TUShortProxy;   // Concrete Implementation of the branch proxy around the data members which are unsigned short
   typedef TImpProxy<UChar_t>    TUCharProxy;    // Concrete Implementation of the branch proxy around the data members which are unsigned char
   typedef TImpProxy<Int_t>      TIntProxy;      // Concrete Implementation of the branch proxy around the data members which are int
   typedef TImpProxy<Long_t>     TLongProxy;     // Concrete Implementation of the branch proxy around the data members which are long
   typedef TImpProxy<Long64_t>   TLong64Proxy;   // Concrete Implementation of the branch proxy around the data members which are long long
   typedef TImpProxy<Short_t>    TShortProxy;    // Concrete Implementation of the branch proxy around the data members which are short
   typedef TImpProxy<Char_t>     TCharProxy;     // Concrete Implementation of the branch proxy around the data members which are char
   typedef TImpProxy<Bool_t>     TBoolProxy;     // Concrete Implementation of the branch proxy around the data members which are bool

   typedef TArrayProxy<TArrayType<Double_t> >   TArrayDoubleProxy;   // Concrete Implementation of the branch proxy around the data members which are array of double
   typedef TArrayProxy<TArrayType<Double32_t> > TArrayDouble32Proxy; // Concrete Implementation of the branch proxy around the data members which are array of double32
   typedef TArrayProxy<TArrayType<Float_t> >    TArrayFloatProxy;    // Concrete Implementation of the branch proxy around the data members which are array of float
   typedef TArrayProxy<TArrayType<Float16_t> >  TArrayFloat16Proxy;  // Concrete Implementation of the branch proxy around the data members which are array of float16
   typedef TArrayProxy<TArrayType<UInt_t> >     TArrayUIntProxy;     // Concrete Implementation of the branch proxy around the data members which are array of unsigned int
   typedef TArrayProxy<TArrayType<ULong_t> >    TArrayULongProxy;    // Concrete Implementation of the branch proxy around the data members which are array of unsigned long
   typedef TArrayProxy<TArrayType<ULong64_t> >  TArrayULong64Proxy;  // Concrete Implementation of the branch proxy around the data members which are array of unsigned long long
   typedef TArrayProxy<TArrayType<UShort_t> >   TArrayUShortProxy;   // Concrete Implementation of the branch proxy around the data members which are array of unsigned short
   typedef TArrayProxy<TArrayType<UChar_t> >    TArrayUCharProxy;    // Concrete Implementation of the branch proxy around the data members which are array of unsigned char
   typedef TArrayProxy<TArrayType<Int_t> >      TArrayIntProxy;      // Concrete Implementation of the branch proxy around the data members which are array of int
   typedef TArrayProxy<TArrayType<Long_t> >     TArrayLongProxy;     // Concrete Implementation of the branch proxy around the data members which are array of long
   typedef TArrayProxy<TArrayType<Long64_t> >   TArrayLong64Proxy;   // Concrete Implementation of the branch proxy around the data members which are array of long long
   typedef TArrayProxy<TArrayType<UShort_t> >   TArrayShortProxy;    // Concrete Implementation of the branch proxy around the data members which are array of short
   //specialized ! typedef TArrayProxy<TArrayType<Char_t> >  TArrayCharProxy; // Concrete Implementation of the branch proxy around the data members which are array of char
   typedef TArrayProxy<TArrayType<Bool_t> >     TArrayBoolProxy;     // Concrete Implementation of the branch proxy around the data members which are array of bool

   typedef TClaImpProxy<Double_t>   TClaDoubleProxy;   // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are double
   typedef TClaImpProxy<Double32_t> TClaDouble32Proxy; // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are double32
   typedef TClaImpProxy<Float_t>    TClaFloatProxy;    // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are float
   typedef TClaImpProxy<Float16_t>  TClaFloat16Proxy;  // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are float16
   typedef TClaImpProxy<UInt_t>     TClaUIntProxy;     // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are unsigned int
   typedef TClaImpProxy<ULong_t>    TClaULongProxy;    // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are unsigned long
   typedef TClaImpProxy<ULong64_t>  TClaULong64Proxy;  // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are unsigned long long
   typedef TClaImpProxy<UShort_t>   TClaUShortProxy;   // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are unsigned short
   typedef TClaImpProxy<UChar_t>    TClaUCharProxy;    // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are unsigned char
   typedef TClaImpProxy<Int_t>      TClaIntProxy;      // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are int
   typedef TClaImpProxy<Long_t>     TClaLongProxy;     // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are long
   typedef TClaImpProxy<Long64_t>   TClaLong64Proxy;   // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are long long
   typedef TClaImpProxy<Short_t>    TClaShortProxy;    // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are short
   typedef TClaImpProxy<Char_t>     TClaCharProxy;     // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are char
   typedef TClaImpProxy<Bool_t>     TClaBoolProxy;     // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are bool

   typedef TClaArrayProxy<TArrayType<Double_t> >    TClaArrayDoubleProxy;   // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of double
   typedef TClaArrayProxy<TArrayType<Double32_t> >  TClaArrayDouble32Proxy; // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of double32
   typedef TClaArrayProxy<TArrayType<Float_t> >     TClaArrayFloatProxy;    // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of float
   typedef TClaArrayProxy<TArrayType<Float16_t> >   TClaArrayFloat16Proxy;  // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of float16
   typedef TClaArrayProxy<TArrayType<UInt_t> >      TClaArrayUIntProxy;     // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of unsigned int
   typedef TClaArrayProxy<TArrayType<ULong_t> >     TClaArrayULongProxy;    // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of unsigned long
   typedef TClaArrayProxy<TArrayType<ULong64_t> >   TClaArrayULong64Proxy;  // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of unsigned long long
   typedef TClaArrayProxy<TArrayType<UShort_t> >    TClaArrayUShortProxy;   // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of unsigned short
   typedef TClaArrayProxy<TArrayType<UChar_t> >     TClaArrayUCharProxy;    // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of nsigned char
   typedef TClaArrayProxy<TArrayType<Int_t> >       TClaArrayIntProxy;      // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of int
   typedef TClaArrayProxy<TArrayType<Long_t> >      TClaArrayLongProxy;     // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of long
   typedef TClaArrayProxy<TArrayType<Long64_t> >    TClaArrayLong64Proxy;   // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of long long
   typedef TClaArrayProxy<TArrayType<UShort_t> >    TClaArrayShortProxy;    // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of short
   typedef TClaArrayProxy<TArrayType<Char_t> >      TClaArrayCharProxy;     // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of char
   typedef TClaArrayProxy<TArrayType<Bool_t> >      TClaArrayBoolProxy;     // Concrete Implementation of the branch proxy around the data members of object in TClonesArray which are array of bool
   //specialized ! typedef TClaArrayProxy<TArrayType<Char_t> >  TClaArrayCharProxy;

   typedef TStlImpProxy<Double_t>   TStlDoubleProxy;   // Concrete Implementation of the branch proxy around an stl container of double
   typedef TStlImpProxy<Double32_t> TStlDouble32Proxy; // Concrete Implementation of the branch proxy around an stl container of double32
   typedef TStlImpProxy<Float_t>    TStlFloatProxy;    // Concrete Implementation of the branch proxy around an stl container of float
   typedef TStlImpProxy<Float16_t>  TStlFloat16Proxy;  // Concrete Implementation of the branch proxy around an stl container of float16_t
   typedef TStlImpProxy<UInt_t>     TStlUIntProxy;     // Concrete Implementation of the branch proxy around an stl container of unsigned int
   typedef TStlImpProxy<ULong_t>    TStlULongProxy;    // Concrete Implementation of the branch proxy around an stl container of unsigned long
   typedef TStlImpProxy<ULong64_t>  TStlULong64Proxy;  // Concrete Implementation of the branch proxy around an stl container of unsigned long long
   typedef TStlImpProxy<UShort_t>   TStlUShortProxy;   // Concrete Implementation of the branch proxy around an stl container of unsigned short
   typedef TStlImpProxy<UChar_t>    TStlUCharProxy;    // Concrete Implementation of the branch proxy around an stl container of unsigned char
   typedef TStlImpProxy<Int_t>      TStlIntProxy;      // Concrete Implementation of the branch proxy around an stl container of int
   typedef TStlImpProxy<Long_t>     TStlLongProxy;     // Concrete Implementation of the branch proxy around an stl container of long
   typedef TStlImpProxy<Long64_t>   TStlLong64Proxy;   // Concrete Implementation of the branch proxy around an stl container of long long
   typedef TStlImpProxy<Short_t>    TStlShortProxy;    // Concrete Implementation of the branch proxy around an stl container of short
   typedef TStlImpProxy<Char_t>     TStlCharProxy;     // Concrete Implementation of the branch proxy around an stl container of char
   typedef TStlImpProxy<Bool_t>     TStlBoolProxy;     // Concrete Implementation of the branch proxy around an stl container of bool

   typedef TStlArrayProxy<TArrayType<Double_t> >    TStlArrayDoubleProxy;   // Concrete Implementation of the branch proxy around an stl container of double
   typedef TStlArrayProxy<TArrayType<Double32_t> >  TStlArrayDouble32Proxy; // Concrete Implementation of the branch proxy around an stl container of double32
   typedef TStlArrayProxy<TArrayType<Float_t> >     TStlArrayFloatProxy;    // Concrete Implementation of the branch proxy around an stl container of float
   typedef TStlArrayProxy<TArrayType<Float16_t> >   TStlArrayFloat16Proxy;  // Concrete Implementation of the branch proxy around an stl container of float16_t
   typedef TStlArrayProxy<TArrayType<UInt_t> >      TStlArrayUIntProxy;     // Concrete Implementation of the branch proxy around an stl container of unsigned int
   typedef TStlArrayProxy<TArrayType<ULong_t> >     TStlArrayULongProxy;    // Concrete Implementation of the branch proxy around an stl container of usigned long
   typedef TStlArrayProxy<TArrayType<ULong64_t> >   TStlArrayULong64Proxy;  // Concrete Implementation of the branch proxy around an stl contained of unsigned long long
   typedef TStlArrayProxy<TArrayType<UShort_t> >    TStlArrayUShortProxy;   // Concrete Implementation of the branch proxy around an stl container of unisgned short
   typedef TStlArrayProxy<TArrayType<UChar_t> >     TStlArrayUCharProxy;    // Concrete Implementation of the branch proxy around an stl container of unsingned char
   typedef TStlArrayProxy<TArrayType<Int_t> >       TStlArrayIntProxy;      // Concrete Implementation of the branch proxy around an stl container of int
   typedef TStlArrayProxy<TArrayType<Long_t> >      TStlArrayLongProxy;     // Concrete Implementation of the branch proxy around an stl container of long
   typedef TStlArrayProxy<TArrayType<Long64_t> >    TStlArrayLong64Proxy;   // Concrete Implementation of the branch proxy around an stl container of long long
   typedef TStlArrayProxy<TArrayType<UShort_t> >    TStlArrayShortProxy;    // Concrete Implementation of the branch proxy around an stl container of UShort_t
   typedef TStlArrayProxy<TArrayType<Char_t> >      TStlArrayCharProxy;     // Concrete Implementation of the branch proxy around an stl container of char
   typedef TStlArrayProxy<TArrayType<Bool_t> >      TStlArrayBoolProxy;     // Concrete Implementation of the branch proxy around an stl container of bool

} // namespace Internal

// Reasonably backward compatible.
using Detail::TBranchProxy;

} // namespace ROOT

#endif

