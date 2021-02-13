// @(#)root/tree:$Id$
// Author: Axel Naumann, 2010-08-02

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreeReaderValue
#define ROOT_TTreeReaderValue

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TTreeReaderValue                                                       //
//                                                                        //
// A simple interface for reading data from trees or chains.              //
//                                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TString.h"
#include "TDictionary.h"
#include "TBranchProxy.h"

#include <type_traits>
#include <vector>
#include <string>

class TBranch;
class TBranchElement;
class TLeaf;
class TTreeReader;

namespace ROOT {
namespace Internal {

/** \class TTreeReaderValueBase
Base class of TTreeReaderValue.
*/

   class TTreeReaderValueBase {
   public:

      // Status flags, 0 is good
      enum ESetupStatus {
         kSetupNotSetup = -7, /// No initialization has happened yet.
         kSetupTreeDestructed = -8, /// The TTreeReader has been destructed / not set.
         kSetupMakeClassModeMismatch = -9, // readers disagree on whether TTree::SetMakeBranch() should be on
         kSetupMissingCounterBranch = -6, /// The array cannot find its counter branch: Array[CounterBranch]
         kSetupMissingBranch = -5, /// The specified branch cannot be found.
         kSetupInternalError = -4, /// Some other error - hopefully the error message helps.
         kSetupMissingDictionary = -3, /// To read this branch, we need a dictionary.
         kSetupMismatch = -2, /// Mismatch of branch type and reader template type.
         kSetupNotACollection = -1, /// The branch class type is not a collection.
         kSetupMatch = 0, /// This branch has been set up, branch data type and reader template type match, reading should succeed.
         kSetupMatchBranch = 7, /// This branch has been set up, branch data type and reader template type match, reading should succeed.
         //kSetupMatchConversion = 1, /// This branch has been set up, the branch data type can be converted to the reader template type, reading should succeed.
         //kSetupMatchConversionCollection = 2, /// This branch has been set up, the data type of the branch's collection elements can be converted to the reader template type, reading should succeed.
         //kSetupMakeClass = 3, /// This branch has been set up, enabling MakeClass mode for it, reading should succeed.
         // kSetupVoidPtr = 4,
         kSetupNoCheck = 5,
         kSetupMatchLeaf = 6 /// This branch (or TLeaf, really) has been set up, reading should succeed.
      };
      enum EReadStatus {
         kReadSuccess = 0, // data read okay
         kReadNothingYet, // data now yet accessed
         kReadError // problem reading data
      };

      EReadStatus ProxyRead() { return (this->*fProxyReadFunc)(); }

      EReadStatus ProxyReadDefaultImpl();

      typedef Bool_t (ROOT::Detail::TBranchProxy::*BranchProxyRead_t)();
      template <BranchProxyRead_t Func>
      ROOT::Internal::TTreeReaderValueBase::EReadStatus ProxyReadTemplate();

      /// Return true if the branch was setup \em and \em read correctly.
      /// Use GetSetupStatus() to only check the setup status.
      Bool_t IsValid() const { return fProxy && 0 == (int)fSetupStatus && 0 == (int)fReadStatus; }
      /// Return this TTreeReaderValue's setup status.
      /// Use this method to check e.g. whether the TTreeReaderValue is correctly setup and ready for reading.
      ESetupStatus GetSetupStatus() const { return fSetupStatus; }
      virtual EReadStatus GetReadStatus() const { return fReadStatus; }

      /// If we are reading a leaf, return the corresponding TLeaf.
      TLeaf* GetLeaf() { return fLeaf; }

      void* GetAddress();

      const char* GetBranchName() const { return fBranchName; }

      virtual ~TTreeReaderValueBase();

   protected:
      TTreeReaderValueBase(TTreeReader* reader, const char* branchname, TDictionary* dict);
      TTreeReaderValueBase(const TTreeReaderValueBase&);
      TTreeReaderValueBase& operator=(const TTreeReaderValueBase&);

      void RegisterWithTreeReader();
      void NotifyNewTree(TTree* newTree);

      TBranch* SearchBranchWithCompositeName(TLeaf *&myleaf, TDictionary *&branchActualType, std::string &err);
      virtual void CreateProxy();
      static const char* GetBranchDataType(TBranch* branch,
                                           TDictionary* &dict,
                                           TDictionary const *curDict);

      virtual const char* GetDerivedTypeName() const = 0;

      Detail::TBranchProxy* GetProxy() const { return fProxy; }

      void MarkTreeReaderUnavailable() { fTreeReader = 0; fSetupStatus = kSetupTreeDestructed; }

      /// Stringify the template argument.
      static std::string GetElementTypeName(const std::type_info& ti);

      int          fHaveLeaf : 1; // Whether the data is in a leaf
      int          fHaveStaticClassOffsets : 1; // Whether !fStaticClassOffsets.empty()
      EReadStatus  fReadStatus : 2; // read status of this data access
      ESetupStatus fSetupStatus = kSetupNotSetup; // setup status of this data access
      TString      fBranchName; // name of the branch to read data from.
      TString      fLeafName;
      TTreeReader* fTreeReader; // tree reader we belong to
      TDictionary* fDict; // type that the branch should contain
      Detail::TBranchProxy* fProxy = nullptr; // proxy for this branch, owned by TTreeReader
      TLeaf*       fLeaf = nullptr;
      std::vector<Long64_t> fStaticClassOffsets;
      typedef EReadStatus (TTreeReaderValueBase::*Read_t)();
      Read_t fProxyReadFunc = &TTreeReaderValueBase::ProxyReadDefaultImpl;      ///<! Pointer to the Read implementation to use.

      // FIXME: re-introduce once we have ClassDefInline!
      //ClassDef(TTreeReaderValueBase, 0);//Base class for accessors to data via TTreeReader

      friend class ::TTreeReader;
   };

} // namespace Internal
} // namespace ROOT


template <typename T>
class R__CLING_PTRCHECK(off) TTreeReaderValue final: public ROOT::Internal::TTreeReaderValueBase {
// R__CLING_PTRCHECK is disabled because pointer / types are checked by CreateProxy().

public:
   using NonConstT_t = typename std::remove_const<T>::type;
   TTreeReaderValue() = delete;
   TTreeReaderValue(TTreeReader& tr, const char* branchname):
      TTreeReaderValueBase(&tr, branchname,
                           TDictionary::GetDictionary(typeid(NonConstT_t))) {}

   /// Return a pointer to the value of the current entry.
   /// Return a nullptr and print an error if no entry has been loaded yet.
   /// The returned address is guaranteed to stay constant while a given TTree is being read from a given file,
   /// unless the branch addresses are manipulated directly (e.g. through TTree::SetBranchAddress()).
   /// The address might also change when the underlying TTree/TFile is switched, e.g. when a TChain switches files.
   T *Get()
   {
      if (!fProxy) {
         Error("TTreeReaderValue::Get()", "Value reader not properly initialized, did you call "
                                          "TTreeReader::Set(Next)Entry() or TTreeReader::Next()?");
         return nullptr;
      }
      void *address = GetAddress(); // Needed to figure out if it's a pointer
      return fProxy->IsaPointer() ? *(T **)address : (T *)address;
   }

   /// Return a pointer to the value of the current entry.
   /// Equivalent to Get().
   T* operator->() { return Get(); }

   /// Return a reference to the value of the current entry.
   /// Equivalent to dereferencing the pointer returned by Get(). Behavior is undefined if no entry has been loaded yet.
   /// Most likely a crash will occur.
   T& operator*() { return *Get(); }

protected:
   // FIXME: use IsA() instead once we have ClassDefTInline
   /// Get the template argument as a string.
   virtual const char* GetDerivedTypeName() const {
      static const std::string sElementTypeName = GetElementTypeName(typeid(T));
      return sElementTypeName.data();
   }

   // FIXME: re-introduce once we have ClassDefTInline!
   //ClassDefT(TTreeReaderValue, 0);//Accessor to data via TTreeReader
};

#endif // ROOT_TTreeReaderValue
