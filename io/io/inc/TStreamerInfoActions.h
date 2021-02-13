// @(#)root/io:$Id$
// Author: Philippe Canal 05/2010

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TStreamerInfoActions
#define ROOT_TStreamerInfoActions

#include <vector>
#include <ROOT/RMakeUnique.hxx>

#include "TStreamerInfo.h"
#include "TVirtualArray.h"

/**
\class TStreamerInfoActions::TConfiguration
\ingroup IO
*/

namespace TStreamerInfoActions {

   /// Base class of the Configurations.
   class TConfiguration {
   protected:
   public:
      typedef TStreamerInfo::TCompInfo_t TCompInfo_t;
      TVirtualStreamerInfo *fInfo;    ///< TStreamerInfo form which the action is derived
      UInt_t                fElemId;  ///< Identifier of the TStreamerElement
      TCompInfo_t          *fCompInfo;///< Access to compiled information (for legacy code)
      Int_t                 fOffset;  ///< Offset within the object
      UInt_t                fLength;  ///< Number of element in a fixed length array.
   public:
      TConfiguration(TVirtualStreamerInfo *info, UInt_t id, TCompInfo_t *compinfo, Int_t offset) : fInfo(info), fElemId(id), fCompInfo(compinfo), fOffset(offset),fLength(1) {};
      TConfiguration(TVirtualStreamerInfo *info, UInt_t id, TCompInfo_t *compinfo, Int_t offset, UInt_t length) : fInfo(info), fElemId(id), fCompInfo(compinfo), fOffset(offset),fLength(length) {};
      virtual ~TConfiguration() {};

      virtual void AddToOffset(Int_t delta);
      virtual void SetMissing();

      virtual TConfiguration *Copy() { return new TConfiguration(*this); }

      virtual void Print() const;
      virtual void PrintDebug(TBuffer &buffer, void *object) const;
   };

   /// Base class of the Configurations for the member wise looping routines.
   class TLoopConfiguration {
   public:
      TVirtualCollectionProxy *fProxy = nullptr;
   public:
      TLoopConfiguration() = default;
      TLoopConfiguration(TVirtualCollectionProxy *proxy) : fProxy(proxy) {}

      // virtual void PrintDebug(TBuffer &buffer, void *object) const;
      virtual ~TLoopConfiguration() {};
      virtual void Print() const;
      virtual void *GetFirstAddress(void *start, const void *end) const = 0;
      virtual TLoopConfiguration* Copy() const = 0; // { return new TLoopConfiguration(*this); }
      virtual TVirtualCollectionProxy* GetCollectionProxy() const { return fProxy; }
   };

   typedef TVirtualCollectionProxy::Next_t Next_t;

   typedef Int_t (*TStreamerInfoAction_t)(TBuffer &buf, void *obj, const TConfiguration *conf);
   typedef Int_t (*TStreamerInfoVecPtrLoopAction_t)(TBuffer &buf, void *iter, const void *end, const TConfiguration *conf);
   typedef Int_t (*TStreamerInfoLoopAction_t)(TBuffer &buf, void *iter, const void *end, const TLoopConfiguration *loopconf, const TConfiguration *conf);

   class TConfiguredAction : public TObject {
   public:
      union {
         TStreamerInfoAction_t           fAction;
         TStreamerInfoVecPtrLoopAction_t fVecPtrLoopAction;
         TStreamerInfoLoopAction_t       fLoopAction;
      };
      TConfiguration              *fConfiguration;
   private:
      // assignment operator must be the default because the 'copy' constructor is actually a move constructor and must be used.
   public:
      TConfiguredAction() : fAction(0), fConfiguration(0) {}
      TConfiguredAction(const TConfiguredAction &rval) : TObject(rval), fAction(rval.fAction), fConfiguration(rval.fConfiguration)
      {
         // WARNING: Technically this is a move constructor ...
         const_cast<TConfiguredAction&>(rval).fConfiguration = 0;
      }
      TConfiguredAction &operator=(const TConfiguredAction &rval)
      {
         // WARNING: Technically this is a move assignment!.

         TConfiguredAction tmp(rval); // this does a move.
         TObject::operator=(tmp);     // we are missing TObject::Swap
         std::swap(fAction,tmp.fAction);
         std::swap(fConfiguration,tmp.fConfiguration);
         return *this;
      };

      TConfiguredAction(TStreamerInfoAction_t action, TConfiguration *conf) : fAction(action), fConfiguration(conf)
      {
         // Usual constructor.
      }
      TConfiguredAction(TStreamerInfoVecPtrLoopAction_t action, TConfiguration *conf) : fVecPtrLoopAction(action), fConfiguration(conf)
      {
         // Usual constructor.
      }
      TConfiguredAction(TStreamerInfoLoopAction_t action, TConfiguration *conf) : fLoopAction(action), fConfiguration(conf)
      {
         // Usual constructor.
      }
      ~TConfiguredAction() {
         // Usual destructor.
         // Idea: the configuration ownership might be moved to a single list so that
         // we can shared them between the optimized and non-optimized list of actions.
         delete fConfiguration;
      }
      void PrintDebug(TBuffer &buffer, void *object) const;

      inline Int_t operator()(TBuffer &buffer, void *object) const {
         return fAction(buffer, object, fConfiguration);
      }

      inline Int_t operator()(TBuffer &buffer, void *start_collection, const void *end_collection) const {
         return fVecPtrLoopAction(buffer, start_collection, end_collection, fConfiguration);
      }

      inline Int_t operator()(TBuffer &buffer, void *start_collection, const void *end_collection, const TLoopConfiguration *loopconf) const {
         return fLoopAction(buffer, start_collection, end_collection, loopconf, fConfiguration);
      }

      ClassDef(TConfiguredAction,0); // A configured action
   };

   struct TIDNode;
   using TIDs = std::vector<TIDNode>;

   // Hold information about unfolded/extracted StreamerElement for
   // a sub-object
   struct TNestedIDs {
      TNestedIDs() = default;
      TNestedIDs(TStreamerInfo *info, Int_t offset) : fInfo(info), fOffset(offset) {}
      ~TNestedIDs() {
         if (fOwnOnfileObject)
            delete fOnfileObject;
      }
      TStreamerInfo *fInfo = nullptr; ///< Not owned.
      TVirtualArray *fOnfileObject = nullptr;
      Bool_t         fOwnOnfileObject = kFALSE;
      Int_t          fOffset;
      TIDs           fIDs;
   };

   // A 'node' in the list of StreamerElement ID, either
   // the index of the element in the current streamerInfo
   // or a set of unfolded/extracted StreamerElement for a sub-object.
   struct TIDNode {
      TIDNode() = default;
      TIDNode(Int_t id) : fElemID(id), fElement(nullptr), fInfo(nullptr) {}
      TIDNode(TStreamerInfo *info, Int_t offset) : fElemID(-1), fElement(nullptr), fInfo(nullptr)  {
         fNestedIDs = std::make_unique<TNestedIDs>(info, offset);
      }
      Int_t fElemID = -1;
      TStreamerElement *fElement = nullptr;
      TStreamerInfo *fInfo = nullptr;
      std::unique_ptr<TNestedIDs> fNestedIDs;
   };

   typedef std::vector<TConfiguredAction> ActionContainer_t;
   class TActionSequence : public TObject {
      TActionSequence() {};
   public:
      struct SequencePtr;
      using SequenceGetter_t = SequencePtr(*)(TStreamerInfo *info, TVirtualCollectionProxy *collectionProxy, TClass *originalClass);

      TActionSequence(TVirtualStreamerInfo *info, UInt_t maxdata) : fStreamerInfo(info), fLoopConfig(0) { fActions.reserve(maxdata); };
      ~TActionSequence() {
         delete fLoopConfig;
      }

      template <typename action_t>
      void AddAction( action_t action, TConfiguration *conf ) {
         fActions.push_back( TConfiguredAction(action, conf) );
      }
      void AddAction(const TConfiguredAction &action ) {
         fActions.push_back( action );
      }

      TVirtualStreamerInfo *fStreamerInfo; ///< StreamerInfo used to derive these actions.
      TLoopConfiguration   *fLoopConfig;   ///< If this is a bundle of memberwise streaming action, this configures the looping
      ActionContainer_t     fActions;

      void AddToOffset(Int_t delta);
      void SetMissing();

      TActionSequence *CreateCopy();
      static TActionSequence *CreateReadMemberWiseActions(TVirtualStreamerInfo *info, TVirtualCollectionProxy &proxy);
      static TActionSequence *CreateWriteMemberWiseActions(TVirtualStreamerInfo *info, TVirtualCollectionProxy &proxy);
      TActionSequence *CreateSubSequence(const std::vector<Int_t> &element_ids, size_t offset);

      TActionSequence *CreateSubSequence(const TIDs &element_ids, size_t offset, SequenceGetter_t create);
      void AddToSubSequence(TActionSequence *sequence, const TIDs &element_ids, Int_t offset, SequenceGetter_t create);

      void Print(Option_t * = "") const;

      // Maybe owner unique_ptr
      struct SequencePtr {
         TStreamerInfoActions::TActionSequence *fSequence = nullptr;
         Bool_t fOwner = kFALSE;

         SequencePtr() = default;

         SequencePtr(SequencePtr &&from) : fSequence(from.fSequence), fOwner(from.fOwner) {
            from.fOwner = false;
         }

         SequencePtr(TStreamerInfoActions::TActionSequence *sequence,  Bool_t owner) : fSequence(sequence), fOwner(owner) {}

         ~SequencePtr() {
            if (fOwner) delete fSequence;
         }

         // Accessor to the pointee.
         TStreamerInfoActions::TActionSequence &operator*() const {
            return *fSequence;
         }

         // Accessor to the pointee
         TStreamerInfoActions::TActionSequence *operator->() const {
            return fSequence;
         }

         // Return true is the pointee is not nullptr.
         operator bool() {
            return fSequence != nullptr;
         }
      };

      // SequenceGetter_t implementations

      static SequencePtr ReadMemberWiseActionsCollectionGetter(TStreamerInfo *info, TVirtualCollectionProxy * /* collectionProxy */, TClass * /* originalClass */) {
         auto seq = info->GetReadMemberWiseActions(kTRUE);
         return {seq, kFALSE};
      }
      static SequencePtr ConversionReadMemberWiseActionsViaProxyGetter(TStreamerInfo *info, TVirtualCollectionProxy *collectionProxy, TClass *originalClass) {
         auto seq = collectionProxy->GetConversionReadMemberWiseActions(originalClass, info->GetClassVersion());
         return {seq, kFALSE};
      }
      static SequencePtr ReadMemberWiseActionsViaProxyGetter(TStreamerInfo *info, TVirtualCollectionProxy *collectionProxy, TClass * /* originalClass */) {
         auto seq = collectionProxy->GetReadMemberWiseActions(info->GetClassVersion());
         return {seq, kFALSE};
      }
      static SequencePtr ReadMemberWiseActionsCollectionCreator(TStreamerInfo *info, TVirtualCollectionProxy *collectionProxy, TClass * /* originalClass */) {
         auto seq = TStreamerInfoActions::TActionSequence::CreateReadMemberWiseActions(info,*collectionProxy);
         return {seq, kTRUE};
      }
      // Creator5() = Creator1;
      static SequencePtr ReadMemberWiseActionsGetter(TStreamerInfo *info, TVirtualCollectionProxy * /* collectionProxy */, TClass * /* originalClass */) {
         auto seq = info->GetReadMemberWiseActions(kFALSE);
         return {seq, kFALSE};
      }

      static SequencePtr WriteMemberWiseActionsCollectionGetter(TStreamerInfo *info, TVirtualCollectionProxy * /* collectionProxy */, TClass * /* originalClass */) {
         auto seq = info->GetWriteMemberWiseActions(kTRUE);
         return {seq, kFALSE};
      }
      static SequencePtr WriteMemberWiseActionsViaProxyGetter(TStreamerInfo *, TVirtualCollectionProxy *collectionProxy, TClass * /* originalClass */) {
         auto seq = collectionProxy->GetWriteMemberWiseActions();
         return {seq, kFALSE};
      }
      static SequencePtr WriteMemberWiseActionsCollectionCreator(TStreamerInfo *info, TVirtualCollectionProxy *collectionProxy, TClass * /* originalClass */) {
         auto seq = TStreamerInfoActions::TActionSequence::CreateWriteMemberWiseActions(info,*collectionProxy);
         return {seq, kTRUE};
      }
      // Creator5() = Creator1;
      static SequencePtr WriteMemberWiseActionsGetter(TStreamerInfo *info, TVirtualCollectionProxy * /* collectionProxy */, TClass * /* originalClass */) {
         auto seq = info->GetWriteMemberWiseActions(kFALSE);
         return {seq, kFALSE};
      }
      ClassDef(TActionSequence,0);
   };

}

#endif // ROOT_TStreamerInfoActions


