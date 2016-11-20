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
// TTreeReaderValue                                                    //
//                                                                        //
// A simple interface for reading data from trees or chains.              //
//                                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TString
#include "TString.h"
#endif
#ifndef ROOT_TDictionary
#include "TDictionary.h"
#endif
#ifndef ROOT_TBranchProxy
#include "TBranchProxy.h"
#endif

#include <type_traits>

class TBranch;
class TBranchElement;
class TLeaf;
class TTreeReader;

namespace ROOT {
namespace Internal {

   class TTreeReaderValueBase {
   public:

      // Status flags, 0 is good
      enum ESetupStatus {
         kSetupNotSetup = -7, /// No initialization has happened yet.
         kSetupTreeDestructed = -8, /// The TTreeReader has been destructed / not set.
         kSetupMakeClassModeMismatch = -7, // readers disagree on whether TTree::SetMakeBranch() should be on
         kSetupMissingCounterBranch = -6, /// The array cannot find its counter branch: Array[CounterBranch]
         kSetupMissingBranch = -5, /// The specified branch cannot be found.
         kSetupInternalError = -4, /// Some other error - hopefully the error message helps.
         kSetupMissingDictionary = -3, /// To read this branch, we need a dictionary.
         kSetupMismatch = -2, /// Mismatch of branch type and reader template type.
         kSetupNotACollection = -1, /// The branch class type is not a collection.
         kSetupMatch = 0, /// This branch has been set up, branch data type and reader template type match, reading should succeed.
         kSetupMatchBranch = 0, /// This branch has been set up, branch data type and reader template type match, reading should succeed.
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

      EReadStatus ProxyRead();

      Bool_t IsValid() const { return fProxy && 0 == (int)fSetupStatus && 0 == (int)fReadStatus; }
      ESetupStatus GetSetupStatus() const { return fSetupStatus; }
      virtual EReadStatus GetReadStatus() const { return fReadStatus; }

      TLeaf* GetLeaf();

      void* GetAddress();

      const char* GetBranchName() const { return fBranchName; }

   protected:
      TTreeReaderValueBase(TTreeReader* reader = 0, const char* branchname = 0, TDictionary* dict = 0);
      TTreeReaderValueBase(const TTreeReaderValueBase&);
      TTreeReaderValueBase& operator=(const TTreeReaderValueBase&);

      virtual ~TTreeReaderValueBase();

      virtual void CreateProxy();
      const char* GetBranchDataType(TBranch* branch,
                                    TDictionary* &dict) const;

      virtual const char* GetDerivedTypeName() const = 0;

      Detail::TBranchProxy* GetProxy() const { return fProxy; }

      void MarkTreeReaderUnavailable() { fTreeReader = 0; fSetupStatus = kSetupTreeDestructed; }

      /// Stringify the template argument.
      static std::string GetElementTypeName(const std::type_info& ti);

      TString      fBranchName; // name of the branch to read data from.
      TString      fLeafName;
      TTreeReader* fTreeReader; // tree reader we belong to
      TDictionary* fDict; // type that the branch should contain
      Detail::TBranchProxy* fProxy; // proxy for this branch, owned by TTreeReader
      TLeaf*       fLeaf;
      Int_t        fLastTreeNumber; // Tree index (in a TChain) that the TLeaf* belongs to.
      ESetupStatus fSetupStatus; // setup status of this data access
      EReadStatus  fReadStatus; // read status of this data access
      std::vector<Long64_t> fStaticClassOffsets;

      // FIXME: re-introduce once we have ClassDefInline!
      //ClassDef(TTreeReaderValueBase, 0);//Base class for accessors to data via TTreeReader

      friend class ::TTreeReader;
   };

} // namespace Internal
} // namespace ROOT


template <typename T>
class TTreeReaderValue: public ROOT::Internal::TTreeReaderValueBase {
public:
   using NonConstT_t = typename std::remove_const<T>::type;
   TTreeReaderValue() {}
   TTreeReaderValue(TTreeReader& tr, const char* branchname):
      TTreeReaderValueBase(&tr, branchname,
                           TDictionary::GetDictionary(typeid(NonConstT_t))) {}

   T* Get() {
      if (!fProxy){
         Error("Get()", "Value reader not properly initialized, did you remember to call TTreeReader.Set(Next)Entry()?");
         return 0;
      }
      void *address = GetAddress(); // Needed to figure out if it's a pointer
      return fProxy->IsaPointer() ? *(T**)address : (T*)address; }
   T* operator->() { return Get(); }
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
