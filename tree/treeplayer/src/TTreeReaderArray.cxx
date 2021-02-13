// @(#)root/treeplayer:$Id$
// Author: Axel Naumann, 2011-09-28

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers and al.        *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TTreeReaderArray.h"

#include "TBranchClones.h"
#include "TBranchElement.h"
#include "TBranchRef.h"
#include "TBranchSTL.h"
#include "TBranchObject.h"
#include "TBranchProxyDirector.h"
#include "TClassEdit.h"
#include "TFriendElement.h"
#include "TFriendProxy.h"
#include "TLeaf.h"
#include "TList.h"
#include "TROOT.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"
#include "TTreeReader.h"
#include "TGenCollectionProxy.h"
#include "TRegexp.h"
#include "ROOT/RMakeUnique.hxx"

#include <memory>

// pin vtable
ROOT::Internal::TVirtualCollectionReader::~TVirtualCollectionReader() {}

namespace {
   using namespace ROOT::Internal;

   // Reader interface for clones arrays
   class TClonesReader: public TVirtualCollectionReader {
   public:
      ~TClonesReader() {}
      TClonesArray* GetCA(ROOT::Detail::TBranchProxy* proxy) {
         if (!proxy->Read()){
            fReadStatus = TTreeReaderValueBase::kReadError;
            Error("TClonesReader::GetCA()", "Read error in TBranchProxy.");
            return 0;
         }
         fReadStatus = TTreeReaderValueBase::kReadSuccess;
         return (TClonesArray*) proxy->GetWhere();
      }
      virtual size_t GetSize(ROOT::Detail::TBranchProxy* proxy) {
         TClonesArray *myClonesArray = GetCA(proxy);
         if (myClonesArray){
            return myClonesArray->GetEntries();
         }
         else return 0;
      }
      virtual void* At(ROOT::Detail::TBranchProxy* proxy, size_t idx) {
         TClonesArray *myClonesArray = GetCA(proxy);
         if (myClonesArray){
            return myClonesArray->UncheckedAt(idx);
         }
         else return 0;
      }
   };

   // Reader interface for STL
   class TSTLReader final: public TVirtualCollectionReader {
   public:
      ~TSTLReader() {}
      TVirtualCollectionProxy* GetCP(ROOT::Detail::TBranchProxy* proxy) {
         if (!proxy->Read()) {
            fReadStatus = TTreeReaderValueBase::kReadError;
            Error("TSTLReader::GetCP()", "Read error in TBranchProxy.");
            return 0;
         }
         if (!proxy->GetWhere()) {
            Error("TSTLReader::GetCP()", "Logic error, proxy object not set in TBranchProxy.");
            return 0;
         }
         fReadStatus = TTreeReaderValueBase::kReadSuccess;
         return (TVirtualCollectionProxy*) proxy->GetCollection();
      }

      virtual size_t GetSize(ROOT::Detail::TBranchProxy* proxy) {
         TVirtualCollectionProxy *myCollectionProxy = GetCP(proxy);
         if (!myCollectionProxy) return 0;
         return myCollectionProxy->Size();
      }

      virtual void* At(ROOT::Detail::TBranchProxy* proxy, size_t idx) {
         TVirtualCollectionProxy *myCollectionProxy = GetCP(proxy);
         if (!myCollectionProxy) return 0;
         if (myCollectionProxy->HasPointers()){
            return *(void**)myCollectionProxy->At(idx);
         }
         else {
            return myCollectionProxy->At(idx);
         }
      }
   };

   class TCollectionLessSTLReader final: public TVirtualCollectionReader {
   private:
      TVirtualCollectionProxy *fLocalCollection;
   public:
      TCollectionLessSTLReader(TVirtualCollectionProxy *proxy) : fLocalCollection(proxy) {}

      TVirtualCollectionProxy* GetCP(ROOT::Detail::TBranchProxy* proxy) {
         if (!proxy->Read()) {
            fReadStatus = TTreeReaderValueBase::kReadError;
            Error("TCollectionLessSTLReader::GetCP()", "Read error in TBranchProxy.");
            return 0;
         }
         if (!proxy->GetWhere()) {
            Error("TCollectionLessSTLReader::GetCP()", "Logic error, proxy object not set in TBranchProxy.");
            return 0;
         }
         fReadStatus = TTreeReaderValueBase::kReadSuccess;
         return fLocalCollection;
      }

      virtual size_t GetSize(ROOT::Detail::TBranchProxy* proxy) {
         TVirtualCollectionProxy *myCollectionProxy = GetCP(proxy);
         if (!myCollectionProxy) return 0;
         /// In the case of std::vector<bool> `PushProxy` also creates a temporary bool variable the address of which
         /// is returned from these calls.
         myCollectionProxy->PopProxy();
         myCollectionProxy->PushProxy(proxy->GetWhere());
         return myCollectionProxy->Size();
      }

      virtual void* At(ROOT::Detail::TBranchProxy* proxy, size_t idx) {
         TVirtualCollectionProxy *myCollectionProxy = GetCP(proxy);
         if (!myCollectionProxy) return 0;
         // Here we do not use a RAII but we empty the proxy to then fill it.
         // This is done because we are returning a pointer and we need to keep
         // alive the memory it points to.
         myCollectionProxy->PopProxy();
         myCollectionProxy->PushProxy(proxy->GetWhere());
         if (myCollectionProxy->HasPointers()){
            return *(void**)myCollectionProxy->At(idx);
         } else {
            return myCollectionProxy->At(idx);
         }
      }
   };


   // Reader interface for leaf list
   // SEE TTreeProxyGenerator.cxx:1319: '//We have a top level raw type'
   class TObjectArrayReader: public TVirtualCollectionReader {
   private:
      Int_t fBasicTypeSize;
   public:
      TObjectArrayReader() : fBasicTypeSize(-1) { }
      ~TObjectArrayReader() {}
      TVirtualCollectionProxy* GetCP(ROOT::Detail::TBranchProxy* proxy) {
         if (!proxy->Read()){
            fReadStatus = TTreeReaderValueBase::kReadError;
            Error("TObjectArrayReader::GetCP()", "Read error in TBranchProxy.");
            return 0;
         }
         fReadStatus = TTreeReaderValueBase::kReadSuccess;
         return (TVirtualCollectionProxy*) proxy->GetCollection();
      }
      virtual size_t GetSize(ROOT::Detail::TBranchProxy* proxy) {
         TVirtualCollectionProxy *myCollectionProxy = GetCP(proxy);
         if (!myCollectionProxy) return 0;
         return myCollectionProxy->Size();
      }
      virtual void* At(ROOT::Detail::TBranchProxy* proxy, size_t idx) {
         if (!proxy->Read()) return 0;

         Int_t objectSize;
         void *array = (void*)proxy->GetStart();

         if (fBasicTypeSize == -1){
            TClass *myClass = proxy->GetClass();
            if (!myClass){
               Error("TObjectArrayReader::At()", "Cannot get class info from branch proxy.");
               return 0;
            }
            objectSize = myClass->GetClassSize();
         }
         else {
            objectSize = fBasicTypeSize;
         }
         return (void*)((Byte_t*)array + (objectSize * idx));
      }

      void SetBasicTypeSize(Int_t size){
         fBasicTypeSize = size;
      }
   };

   template <class BASE>
   class TUIntOrIntReader: public BASE {
   private:
      // The index can be of type int or unsigned int.
      std::unique_ptr<TTreeReaderValueBase> fSizeReader;
      bool fIsUnsigned = false;

   protected:
      template <class T>
      TTreeReaderValue<T>& GetSizeReader() {
         return *static_cast<TTreeReaderValue<T>*>(fSizeReader.get());
      }

   public:
      template <class... ARGS>
      TUIntOrIntReader(TTreeReader *treeReader, const char *leafName,
                       ARGS&&... args):
         BASE(std::forward<ARGS>(args)...)
      {
         std::string foundLeafName = leafName;
         TLeaf* sizeLeaf = treeReader->GetTree()->FindLeaf(foundLeafName.c_str());

         if (!sizeLeaf) {
            // leafName might be "top.currentParent.N". But "N" might really be "top.N"!
            // Strip parents until we find the leaf.
            std::string leafNameNoParent = leafName;
            std::string parent;
            auto posLastDot = leafNameNoParent.rfind('.');
            if (posLastDot != leafNameNoParent.npos) {
               parent = leafNameNoParent.substr(0, posLastDot);
               leafNameNoParent.erase(0, posLastDot + 1);
            }

            do {
               if (!sizeLeaf && !parent.empty()) {
                  auto posLastDotParent = parent.rfind('.');
                  if (posLastDotParent != parent.npos)
                     parent = parent.substr(0, posLastDot);
                  else
                     parent.clear();
               }

               foundLeafName = parent;
               if (!parent.empty())
                  foundLeafName += ".";
               foundLeafName += leafNameNoParent;
               sizeLeaf = treeReader->GetTree()->FindLeaf(foundLeafName.c_str());
            } while (!sizeLeaf && !parent.empty());
         }

         if (!sizeLeaf) {
            Error("TUIntOrIntReader", "Cannot find leaf count for %s or any parent branch!", leafName);
            return;
         }

         fIsUnsigned = sizeLeaf->IsUnsigned();
         if (fIsUnsigned) {
            fSizeReader.reset(new TTreeReaderValue<UInt_t>(*treeReader, foundLeafName.c_str()));
         } else {
            fSizeReader.reset(new TTreeReaderValue<Int_t>(*treeReader, foundLeafName.c_str()));
         }
      }

      size_t GetSize(ROOT::Detail::TBranchProxy* /*proxy*/) override {
         if (fIsUnsigned)
            return *GetSizeReader<UInt_t>();
         return *GetSizeReader<Int_t>();
      }
   };

   class TArrayParameterSizeReader: public TUIntOrIntReader<TObjectArrayReader> {
   public:
      TArrayParameterSizeReader(TTreeReader *treeReader, const char *branchName):
         TUIntOrIntReader<TObjectArrayReader>(treeReader, branchName) {}
   };

   // Reader interface for fixed size arrays
   class TArrayFixedSizeReader : public TObjectArrayReader {
   private:
      Int_t fSize;

   public:
      TArrayFixedSizeReader(Int_t sizeArg) : fSize(sizeArg) {}

      virtual size_t GetSize(ROOT::Detail::TBranchProxy* /*proxy*/) { return fSize; }
   };

   class TBasicTypeArrayReader final: public TVirtualCollectionReader {
   public:
      ~TBasicTypeArrayReader() {}

      TVirtualCollectionProxy* GetCP (ROOT::Detail::TBranchProxy *proxy) {
         if (!proxy->Read()){
            fReadStatus = TTreeReaderValueBase::kReadError;
            Error("TBasicTypeArrayReader::GetCP()", "Read error in TBranchProxy.");
            return 0;
         }
         fReadStatus = TTreeReaderValueBase::kReadSuccess;
         return (TVirtualCollectionProxy*) proxy->GetCollection();
      }

      virtual size_t GetSize(ROOT::Detail::TBranchProxy* proxy){
         TVirtualCollectionProxy *myCollectionProxy = GetCP(proxy);
         if (!myCollectionProxy) return 0;
         return myCollectionProxy->Size();
      }

      virtual void* At(ROOT::Detail::TBranchProxy* proxy, size_t idx){
         TVirtualCollectionProxy *myCollectionProxy = GetCP(proxy);
         if (!myCollectionProxy) return 0;
         return (Byte_t*)myCollectionProxy->At(idx) + proxy->GetOffset();
      }
   };

   class TBasicTypeClonesReader final: public TClonesReader {
   private:
      Int_t fOffset;
   public:
      TBasicTypeClonesReader(Int_t offsetArg) : fOffset(offsetArg) {}

      virtual void* At(ROOT::Detail::TBranchProxy* proxy, size_t idx){
         TClonesArray *myClonesArray = GetCA(proxy);
         if (!myClonesArray) return 0;
         return (Byte_t*)myClonesArray->At(idx) + fOffset;
      }
   };

   class TLeafReader : public TVirtualCollectionReader {
   private:
      TTreeReaderValueBase *fValueReader;
      Int_t fElementSize;
   public:
      TLeafReader(TTreeReaderValueBase *valueReaderArg) : fValueReader(valueReaderArg), fElementSize(-1) {}

      virtual size_t GetSize(ROOT::Detail::TBranchProxy* /*proxy*/){
         TLeaf *myLeaf = fValueReader->GetLeaf();
         return myLeaf ? myLeaf->GetLen() : 0; // Error will be printed by GetLeaf
      }

      virtual void* At(ROOT::Detail::TBranchProxy* /*proxy*/, size_t idx){
         ProxyRead();
         void *address = fValueReader->GetAddress();
         if (fElementSize == -1){
            TLeaf *myLeaf = fValueReader->GetLeaf();
            if (!myLeaf) return 0; // Error will be printed by GetLeaf
            fElementSize = myLeaf->GetLenType();
         }
         return (Byte_t*)address + (fElementSize * idx);
      }

   protected:
      void ProxyRead(){
         fValueReader->ProxyRead();
      }
   };

   class TLeafParameterSizeReader: public TUIntOrIntReader<TLeafReader> {
   public:
      TLeafParameterSizeReader(TTreeReader *treeReader, const char *leafName,
                               TTreeReaderValueBase *valueReaderArg) :
         TUIntOrIntReader<TLeafReader>(treeReader, leafName, valueReaderArg) {}

      size_t GetSize(ROOT::Detail::TBranchProxy* proxy) override {
         ProxyRead();
         return TUIntOrIntReader<TLeafReader>::GetSize(proxy);
      }
   };
}



ClassImp(TTreeReaderArrayBase);

////////////////////////////////////////////////////////////////////////////////
/// Create the proxy object for our branch.

void ROOT::Internal::TTreeReaderArrayBase::CreateProxy()
{
   if (fProxy) {
      return;
   }

   fSetupStatus = kSetupInternalError; // Fallback; set to something concrete below.
   if (!fTreeReader) {
      Error("TTreeReaderArrayBase::CreateProxy()", "TTreeReader object not set / available for branch %s!",
            fBranchName.Data());
      fSetupStatus = kSetupTreeDestructed;
      return;
   }
   if (!fDict) {
      TBranch* br = fTreeReader->GetTree()->GetBranch(fBranchName);
      const char* brDataType = "{UNDETERMINED}";
      if (br) {
         TDictionary* dictUnused = 0;
         brDataType = GetBranchDataType(br, dictUnused, fDict);
      }
      Error("TTreeReaderArrayBase::CreateProxy()", "The template argument type T of %s accessing branch %s (which contains data of type %s) is not known to ROOT. You will need to create a dictionary for it.",
            GetDerivedTypeName(), fBranchName.Data(), brDataType);
      fSetupStatus = kSetupMissingDictionary;
      return;
   }

   // Access a branch's collection content (not the collection itself)
   // through a proxy.
   // Search for the branchname, determine what it contains, and wire the
   // TBranchProxy representing it to us so we can access its data.

   TDictionary* branchActualType = 0;
   TBranch* branch = nullptr;
   TLeaf *myLeaf = nullptr;
   if (!GetBranchAndLeaf(branch, myLeaf, branchActualType))
      return;

   if (!fDict) {
      Error("TTreeReaderArrayBase::CreateProxy()",
            "No dictionary for branch %s.", fBranchName.Data());
      return;
   }

   TNamedBranchProxy* namedProxy = fTreeReader->FindProxy(fBranchName);
   if (namedProxy) {
      if (namedProxy->GetContentDict() == fDict) {
         fSetupStatus = kSetupMatch;
         fProxy = namedProxy->GetProxy();
         SetImpl(branch, myLeaf);
         return;
      }

      // Update named proxy's dictionary
      if (!namedProxy->GetContentDict()) {
         namedProxy->SetContentDict(fDict);
         fProxy = namedProxy->GetProxy();
         if (fProxy)
            fSetupStatus = kSetupMatch;
      } else {
         Error("TTreeReaderArrayBase::CreateProxy()",
               "Type ambiguity (want %s, have %s) for branch %s.",
               fDict->GetName(), namedProxy->GetContentDict()->GetName(), fBranchName.Data());
      }
   }
   else {
      TString membername;

      bool isTopLevel = branch->GetMother() == branch;
      if (!isTopLevel) {
         membername = strrchr(branch->GetName(), '.');
         if (membername.IsNull()) {
            membername = branch->GetName();
         }
      }
      auto director = fTreeReader->fDirector;
      // Determine if the branch is actually in a Friend TTree and if so which.
      if (branch->GetTree() != fTreeReader->GetTree()->GetTree()) {
         // It is in a friend, let's find the 'index' in the list of friend ...
         int index = -1;
         int current = 0;
         for(auto fe : TRangeDynCast<TFriendElement>( fTreeReader->GetTree()->GetTree()->GetListOfFriends())) {
            if (branch->GetTree() == fe->GetTree()) {
               index = current;
            }
            ++current;
         }
         if (index == -1) {
            Error("TTreeReaderArrayBase::CreateProxy()", "The branch %s is contained in a Friend TTree that is not directly attached to the main.\n"
                  "This is not yet supported by TTreeReader.",
                  fBranchName.Data());
            return;
         }
         TFriendProxy *feproxy = nullptr;
         if ((size_t)index < fTreeReader->fFriendProxies.size()) {
            feproxy = fTreeReader->fFriendProxies.at(index);
         }
         if (!feproxy) {
            feproxy = new ROOT::Internal::TFriendProxy(director, fTreeReader->GetTree(), index);
            fTreeReader->fFriendProxies.resize(index+1);
            fTreeReader->fFriendProxies.at(index) = feproxy;
         }
         director = feproxy->GetDirector();
      }
      namedProxy = new TNamedBranchProxy(director, branch, fBranchName, membername);
      fTreeReader->AddProxy(namedProxy);
      fProxy = namedProxy->GetProxy();
      if (fProxy)
         fSetupStatus = kSetupMatch;
      else
         fSetupStatus = kSetupMismatch;
   }

   if (!myLeaf){
      TString branchActualTypeName;
      const char* nonCollTypeName = GetBranchContentDataType(branch, branchActualTypeName, branchActualType);
      if (nonCollTypeName) {
         Error("TTreeReaderArrayBase::CreateContentProxy()", "The branch %s contains data of type %s, which should be accessed through a TTreeReaderValue< %s >.",
               fBranchName.Data(), nonCollTypeName, nonCollTypeName);
         if (fSetupStatus == kSetupInternalError)
            fSetupStatus = kSetupNotACollection;
         fProxy = 0;
         return;
      }
      if (!branchActualType) {
         if (branchActualTypeName.IsNull()) {
            Error("TTreeReaderArrayBase::CreateContentProxy()", "Cannot determine the type contained in the collection of branch %s. That's weird - please report!",
                  fBranchName.Data());
         } else {
            Error("TTreeReaderArrayBase::CreateContentProxy()", "The branch %s contains data of type %s, which does not have a dictionary.",
                  fBranchName.Data(), branchActualTypeName.Data());
            if (fSetupStatus == kSetupInternalError)
               fSetupStatus = kSetupMissingDictionary;
         }
         fProxy = 0;
         return;
      }

      if (fDict != branchActualType) {
         Error("TTreeReaderArrayBase::CreateContentProxy()", "The branch %s contains data of type %s. It cannot be accessed by a TTreeReaderArray<%s>",
               fBranchName.Data(), branchActualType->GetName(), fDict->GetName());
         if (fSetupStatus == kSetupInternalError || fSetupStatus >= 0)
            fSetupStatus = kSetupMismatch;

         // Update named proxy's dictionary
         if (!namedProxy->GetContentDict()) {
            namedProxy->SetContentDict(fDict);
         }

         // fProxy = 0;
         // return;
      }
   }

   SetImpl(branch, myLeaf);
}

////////////////////////////////////////////////////////////////////////////////
/// Determine the branch / leaf and its type; reset fProxy / fSetupStatus on error.

bool ROOT::Internal::TTreeReaderArrayBase::GetBranchAndLeaf(TBranch* &branch, TLeaf* &myLeaf,
                                                            TDictionary* &branchActualType) {
   myLeaf = nullptr;
   branch = fTreeReader->GetTree()->GetBranch(fBranchName);
   if (branch)
      return true;

   if (!fBranchName.Contains(".")) {
      Error("TTreeReaderArrayBase::GetBranchAndLeaf()", "The tree does not have a branch called %s. You could check with TTree::Print() for available branches.", fBranchName.Data());
      fSetupStatus = kSetupMissingBranch;
      fProxy = 0;
      return false;
   }

   TRegexp leafNameExpression ("\\.[a-zA-Z0-9_]+$");
   TString leafName (fBranchName(leafNameExpression));
   TString branchName = fBranchName(0, fBranchName.Length() - leafName.Length());
   branch = fTreeReader->GetTree()->GetBranch(branchName);
   if (!branch){
      Error("TTreeReaderArrayBase::GetBranchAndLeaf()", "The tree does not have a branch called %s. You could check with TTree::Print() for available branches.", fBranchName.Data());
      fSetupStatus = kSetupMissingBranch;
      fProxy = 0;
      return false;
   }

   myLeaf = branch->GetLeaf(TString(leafName(1, leafName.Length())));
   if (!myLeaf){
      Error("TTreeReaderArrayBase::GetBranchAndLeaf()", "The tree does not have a branch, nor a sub-branch called %s. You could check with TTree::Print() for available branches.", fBranchName.Data());
      fSetupStatus = kSetupMissingBranch;
      fProxy = 0;
      return false;
   }

   TDictionary *tempDict = TDictionary::GetDictionary(myLeaf->GetTypeName());
   if (!tempDict){
      Error("TTreeReaderArrayBase::GetBranchAndLeaf()", "Failed to get the dictionary for %s.", myLeaf->GetTypeName());
      fSetupStatus = kSetupMissingDictionary;
      fProxy = 0;
      return false;
   }

   if (tempDict->IsA() == TDataType::Class() && TDictionary::GetDictionary(((TDataType*)tempDict)->GetTypeName()) == fDict){
      //fLeafOffset = myLeaf->GetOffset() / 4;
      branchActualType = fDict;
      fLeaf = myLeaf;
      fBranchName = branchName;
      fLeafName = leafName(1, leafName.Length());
      fHaveLeaf = (fLeafName.Length() > 0);
      fSetupStatus = kSetupMatchLeaf;
   }
   else {
      Error("TTreeReaderArrayBase::GetBranchAndLeaf()", "Leaf of type %s cannot be read by TTreeReaderValue<%s>.", myLeaf->GetTypeName(), fDict->GetName());
      fProxy = 0;
      fSetupStatus = kSetupMismatch;
      return false;
   }
   return true;
}




////////////////////////////////////////////////////////////////////////////////
/// Create the TVirtualCollectionReader object for our branch.

void ROOT::Internal::TTreeReaderArrayBase::SetImpl(TBranch* branch, TLeaf* myLeaf)
{
   if (fImpl)
      return;

   // Access a branch's collection content (not the collection itself)
   // through a proxy.
   // Search for the branchname, determine what it contains, and wire the
   // TBranchProxy representing it to us so we can access its data.
   // A proxy for branch must not have been created before (i.e. check
   // fProxies before calling this function!)

   if (myLeaf){
      if (!myLeaf->GetLeafCount()){
         fImpl = std::make_unique<TLeafReader>(this);
      }
      else {
         TString leafFullName = myLeaf->GetBranch()->GetName();
         leafFullName += ".";
         leafFullName += myLeaf->GetLeafCount()->GetName();
         fImpl = std::make_unique<TLeafParameterSizeReader>(fTreeReader, leafFullName.Data(), this);
      }
      fSetupStatus = kSetupMatchLeaf;
   }
   else if (branch->IsA() == TBranchElement::Class()) {
      TBranchElement* branchElement = ((TBranchElement*)branch);

      TStreamerInfo *streamerInfo = branchElement->GetInfo();
      Int_t id = branchElement->GetID();

      if (id >= 0){ // Not root node?
         // Int_t offset = streamerInfo->GetOffsets()[id];
         TStreamerElement *element = (TStreamerElement*)streamerInfo->GetElements()->At(id);
         // Bool_t isPointer = element->IsaPointer();
         // TClass *classPointer = element->GetClassPointer();

         if (fSetupStatus == kSetupInternalError)
            fSetupStatus = kSetupMatch;
         if (element->IsA() == TStreamerSTL::Class()){
            fImpl = std::make_unique<TSTLReader>();
         }
         else if (element->IsA() == TStreamerObject::Class()){
            //fImpl = new TObjectArrayReader(); // BArray[12]

            if (element->GetClass() == TClonesArray::Class()){
               fImpl = std::make_unique<TClonesReader>();
            }
            else if (branchElement->GetType() == TBranchElement::kSTLMemberNode){
               fImpl = std::make_unique<TBasicTypeArrayReader>();
            }
            else if (branchElement->GetType() == TBranchElement::kClonesMemberNode){
               // TBasicTypeClonesReader should work for object
               fImpl = std::make_unique<TBasicTypeClonesReader>(element->GetOffset());
            }
            else {
               fImpl = std::make_unique<TArrayFixedSizeReader>(element->GetArrayLength());
            }
         }
         else if (element->IsA() == TStreamerLoop::Class()) {
            fImpl = std::make_unique<TArrayParameterSizeReader>(fTreeReader, branchElement->GetBranchCount()->GetName());
         }
         else if (element->IsA() == TStreamerBasicType::Class()){
            if (branchElement->GetType() == TBranchElement::kSTLMemberNode){
               fImpl = std::make_unique<TBasicTypeArrayReader>();
            }
            else if (branchElement->GetType() == TBranchElement::kClonesMemberNode){
               fImpl = std::make_unique<TBasicTypeClonesReader>(element->GetOffset());
            }
            else {
               fImpl = std::make_unique<TArrayFixedSizeReader>(element->GetArrayLength());
               ((TObjectArrayReader*)fImpl.get())->SetBasicTypeSize(((TDataType*)fDict)->Size());
            }
         }
         else if (element->IsA() == TStreamerBasicPointer::Class()) {
            fImpl = std::make_unique<TArrayParameterSizeReader>(fTreeReader, branchElement->GetBranchCount()->GetName());
            ((TArrayParameterSizeReader*)fImpl.get())->SetBasicTypeSize(((TDataType*)fDict)->Size());
         }
         else if (element->IsA() == TStreamerBase::Class()){
            fImpl = std::make_unique<TClonesReader>();
         } else {
            Error("TTreeReaderArrayBase::SetImpl()",
                  "Cannot read branch %s: unhandled streamer element type %s",
                  fBranchName.Data(), element->IsA()->GetName());
            fSetupStatus = kSetupInternalError;
         }
      }
      else { // We are at root node?
         if (branchElement->GetClass()->GetCollectionProxy()){
            fImpl = std::make_unique<TCollectionLessSTLReader>(branchElement->GetClass()->GetCollectionProxy());
         }
      }
   } else if (branch->IsA() == TBranch::Class()) {
      auto topLeaf = branch->GetLeaf(branch->GetName());
      if (!topLeaf) {
         Error("TTreeReaderArrayBase::SetImpl", "Failed to get the top leaf from the branch");
         fSetupStatus = kSetupMissingBranch;
         return;
      }
      // We could have used GetLeafCounter, but it does not work well with Double32_t and Float16_t: ROOT-10149
      auto sizeLeaf = topLeaf->GetLeafCount();
      if (fSetupStatus == kSetupInternalError)
         fSetupStatus = kSetupMatch;
      if (!sizeLeaf) {
         fImpl = std::make_unique<TArrayFixedSizeReader>(topLeaf->GetLenStatic());
      }
      else {
         fImpl = std::make_unique<TArrayParameterSizeReader>(fTreeReader, sizeLeaf->GetName());
      }
      ((TObjectArrayReader*)fImpl.get())->SetBasicTypeSize(((TDataType*)fDict)->Size());
   } else if (branch->IsA() == TBranchClones::Class()) {
      Error("TTreeReaderArrayBase::SetImpl", "Support for branches of type TBranchClones not implemented");
      fSetupStatus = kSetupInternalError;
   } else if (branch->IsA() == TBranchObject::Class()) {
      Error("TTreeReaderArrayBase::SetImpl", "Support for branches of type TBranchObject not implemented");
      fSetupStatus = kSetupInternalError;
   } else if (branch->IsA() == TBranchSTL::Class()) {
      Error("TTreeReaderArrayBase::SetImpl", "Support for branches of type TBranchSTL not implemented");
      fImpl = std::make_unique<TSTLReader>();
      fSetupStatus = kSetupInternalError;
   } else if (branch->IsA() == TBranchRef::Class()) {
      Error("TTreeReaderArrayBase::SetImpl", "Support for branches of type TBranchRef not implemented");
      fSetupStatus = kSetupInternalError;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Access a branch's collection content (not the collection itself)
/// through a proxy.
/// Retrieve the type of data contained in the collection stored by branch;
/// put its dictionary into dict, If there is no dictionary, put its type
/// name into contentTypeName.
/// The contentTypeName is set to NULL if the branch does not
/// contain a collection; in that case, the type of the branch is returned.
/// In all other cases, NULL is returned.

const char* ROOT::Internal::TTreeReaderArrayBase::GetBranchContentDataType(TBranch* branch,
                                                                 TString& contentTypeName,
                                                                 TDictionary* &dict)
{
   dict = nullptr;
   contentTypeName = "";
   if (branch->IsA() == TBranchElement::Class()) {
      TBranchElement* brElement = (TBranchElement*)branch;
      if (brElement->GetType() == 4
          || brElement->GetType() == 3) {
         TVirtualCollectionProxy* collProxy = brElement->GetCollectionProxy();
         if (collProxy) {
            TClass *myClass = collProxy->GetValueClass();
            if (!myClass){
               Error("TTreeReaderArrayBase::GetBranchContentDataType()", "Could not get value class.");
               return 0;
            }
            dict = TDictionary::GetDictionary(myClass->GetName());
            if (!dict) dict = TDataType::GetDataType(collProxy->GetType());
         }
         if (!dict) {
            // We don't know the dictionary, thus we need the content's type name.
            // Determine it.
            if (brElement->GetType() == 3) {
               contentTypeName = brElement->GetClonesName();
               dict = TDictionary::GetDictionary(brElement->GetClonesName());
               return 0;
            }
            // STL:
            TClassEdit::TSplitType splitType(brElement->GetClassName());
            int isSTLCont = splitType.IsSTLCont();
            if (!isSTLCont) {
               Error("TTreeReaderArrayBase::GetBranchContentDataType()", "Cannot determine STL collection type of %s stored in branch %s", brElement->GetClassName(), branch->GetName());
               return brElement->GetClassName();
            }
            bool isMap = isSTLCont == ROOT::kSTLmap
               || isSTLCont == ROOT::kSTLmultimap;
            if (isMap) contentTypeName = "std::pair< ";
            contentTypeName += splitType.fElements[1];
            if (isMap) {
               contentTypeName += splitType.fElements[2];
               contentTypeName += " >";
            }
            return 0;
         }
         return 0;
      } else if (brElement->GetType() == 31
                 || brElement->GetType() == 41) {
         // it's a member, extract from GetClass()'s streamer info
         TClass* clData = 0;
         EDataType dtData = kOther_t;
         int ExpectedTypeRet = brElement->GetExpectedType(clData, dtData);
         if (ExpectedTypeRet == 0) {
            dict = clData;
            if (!dict) {
               dict = TDataType::GetDataType(dtData);
            }
            if (!dict) {
               Error("TTreeReaderArrayBase::GetBranchContentDataType()", "The branch %s contains a data type %d for which the dictionary cannot be retrieved.",
                     branch->GetName(), (int)dtData);
               contentTypeName = TDataType::GetTypeName(dtData);
               return 0;
            }
            return 0;
         } else if (ExpectedTypeRet == 1) {
            int brID = brElement->GetID();
            if (brID == -1) {
               // top
               Error("TTreeReaderArrayBase::GetBranchContentDataType()", "The branch %s contains data of type %s for which the dictionary does not exist. It's needed.",
                     branch->GetName(), brElement->GetClassName());
               contentTypeName = brElement->GetClassName();
               return 0;
            }
            // Either the data type name doesn't have an EDataType entry
            // or the streamer info doesn't have a TClass* attached.
            TStreamerElement* element =
               (TStreamerElement*) brElement->GetInfo()->GetElement(brID);
            contentTypeName = element->GetTypeName();
            return 0;
         }
         /* else (ExpectedTypeRet == 2)*/
         // The streamer info entry cannot be found.
         // TBranchElement::GetExpectedType() has already complained.
         return "{CANNOT DETERMINE TBranchElement DATA TYPE}";
      }
      else if (brElement->GetType() == TBranchElement::kLeafNode){
         TStreamerInfo *streamerInfo = brElement->GetInfo();
         Int_t id = brElement->GetID();

         if (id >= 0){
            TStreamerElement *element = (TStreamerElement*)streamerInfo->GetElements()->At(id);

            if (element->IsA() == TStreamerSTL::Class()){
               TClass *myClass = brElement->GetCurrentClass();
               if (!myClass){
                  Error("TTreeReaderArrayBase::GetBranchDataType()", "Could not get class from branch element.");
                  return 0;
               }
               TVirtualCollectionProxy *myCollectionProxy = myClass->GetCollectionProxy();
               if (!myCollectionProxy){
                  Error("TTreeReaderArrayBase::GetBranchDataType()", "Could not get collection proxy from STL class");
                  return 0;
               }
               // Try getting the contained class
               dict = myCollectionProxy->GetValueClass();
               // If it fails, try to get the contained type as a primitive type
               if (!dict) dict = TDataType::GetDataType(myCollectionProxy->GetType());
               if (!dict){
                  Error("TTreeReaderArrayBase::GetBranchDataType()", "Could not get valueClass from collectionProxy.");
                  return 0;
               }
               contentTypeName = dict->GetName();
               return 0;
            }
            else if (element->IsA() == TStreamerObject::Class() && !strcmp(element->GetTypeName(), "TClonesArray")){
               if (!fProxy->Setup() || !fProxy->Read()){
                  Error("TTreeReaderArrayBase::GetBranchContentDataType()", "Failed to get type from proxy, unable to check type");
                  contentTypeName = "UNKNOWN";
                  dict = 0;
                  return contentTypeName;
               }
               TClonesArray *myArray = (TClonesArray*)fProxy->GetWhere();
               dict = myArray->GetClass();
               contentTypeName = dict->GetName();
               return 0;
            }
            else {
               dict = brElement->GetCurrentClass();
               if (!dict) {
                  TDictionary *myDataType = TDictionary::GetDictionary(brElement->GetTypeName());
                  dict = TDataType::GetDataType((EDataType)((TDataType*)myDataType)->GetType());
               }
               contentTypeName = brElement->GetTypeName();
               return 0;
            }
         }
         if (brElement->GetCurrentClass() == TClonesArray::Class()){
            contentTypeName = "TClonesArray";
            Warning("TTreeReaderArrayBase::GetBranchContentDataType()", "Not able to check type correctness, ignoring check");
            dict = fDict;
            fSetupStatus = kSetupNoCheck;
         }
         else if (!dict && (branch->GetSplitLevel() == 0 || brElement->GetClass()->GetCollectionProxy())){
            // Try getting the contained class
            dict = brElement->GetClass()->GetCollectionProxy()->GetValueClass();
            // If it fails, try to get the contained type as a primitive type
            if (!dict) dict = TDataType::GetDataType(brElement->GetClass()->GetCollectionProxy()->GetType());
            if (dict) contentTypeName = dict->GetName();
            return 0;
         }
         else if (!dict){
            dict = brElement->GetClass();
            contentTypeName = dict->GetName();
            return 0;
         }

         return 0;
      }
      return 0;
   } else if (branch->IsA() == TBranch::Class()
              || branch->IsA() == TBranchObject::Class()
              || branch->IsA() == TBranchSTL::Class()) {
      const char* dataTypeName = branch->GetClassName();
      if ((!dataTypeName || !dataTypeName[0])
          && branch->IsA() == TBranch::Class()) {
         auto myLeaf = branch->GetLeaf(branch->GetName());
         if (myLeaf){
            auto myDataType = TDictionary::GetDictionary(myLeaf->GetTypeName());
            if (myDataType && myDataType->IsA() == TDataType::Class()){
               auto typeEnumConstant = EDataType(((TDataType*)myDataType)->GetType());
               // We need to consider Double32_t and Float16_t as dounle and float respectively
               // since this is the type the user uses to instantiate the TTreeReaderArray template.
               if (typeEnumConstant == kDouble32_t) typeEnumConstant = kDouble_t;
               else if (typeEnumConstant == kFloat16_t) typeEnumConstant = kFloat_t;
               dict = TDataType::GetDataType(typeEnumConstant);
               contentTypeName = myLeaf->GetTypeName();
               return 0;
            }
         }

         // leaflist. Can't represent.
         Error("TTreeReaderArrayBase::GetBranchContentDataType()", "The branch %s was created using a leaf list and cannot be represented as a C++ type. Please access one of its siblings using a TTreeReaderArray:", branch->GetName());
         TIter iLeaves(branch->GetListOfLeaves());
         TLeaf* leaf = 0;
         while ((leaf = (TLeaf*) iLeaves())) {
            Error("TTreeReaderArrayBase::GetBranchContentDataType()", "   %s.%s", branch->GetName(), leaf->GetName());
         }
         return 0;
      }
      if (dataTypeName) dict = TDictionary::GetDictionary(dataTypeName);
      if (branch->IsA() == TBranchSTL::Class()){
         Warning("TTreeReaderArrayBase::GetBranchContentDataType()", "Not able to check type correctness, ignoring check");
         dict = fDict;
         fSetupStatus = kSetupNoCheck;
         return 0;
      }
      return dataTypeName;
   } else if (branch->IsA() == TBranchClones::Class()) {
      dict = TClonesArray::Class();
      return "TClonesArray";
   } else if (branch->IsA() == TBranchRef::Class()) {
      // Can't represent.
      Error("TTreeReaderArrayBase::GetBranchContentDataType()", "The branch %s is a TBranchRef and cannot be represented as a C++ type.", branch->GetName());
      return 0;
   } else {
      Error("TTreeReaderArrayBase::GetBranchContentDataType()", "The branch %s is of type %s - something that is not handled yet.", branch->GetName(), branch->IsA()->GetName());
      return 0;
   }

   return 0;
}
