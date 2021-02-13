// @(#)root/meta:$Id$
// Author: Rene Brun   07/01/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TCheckHashRecursiveRemoveConsistency
#define ROOT_TCheckHashRecursiveRemoveConsistency

#include "TBaseClass.h"
#include "TClass.h"
#include "TError.h"
#include "TMethod.h"

#include <list>
#include <utility>
#include <iostream>
#include <mutex>

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCheckHashRecursiveRemoveConsistency                                 //
//                                                                      //
// Utility class to discover whether a class that overload              //
// TObject::Hash also (as required) calls RecursiveRemove in its        //
// destructor.                                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

namespace ROOT {
namespace Internal {

class TCheckHashRecursiveRemoveConsistency : public TObject {
public:
   struct Value {
      ULong_t  fRecordedHash;
      TObject *fObjectPtr;
   };
   using Value_t = Value; // std::pair<ULong_t, TObject*>;

   std::list<Value> fCont;
   std::mutex       fMutex;

public:
   // Default constructor.  Adds object to the list of
   // cleanups.
   TCheckHashRecursiveRemoveConsistency()
   {
      SetBit(kMustCleanup);
      gROOT->GetListOfCleanups()->Add(this);
   }

   // Destructor.  This class does not overload
   // Hash so it can rely on the base class to call
   // RecursiveRemove (and hence remove this from the list
   // of cleanups).
   ~TCheckHashRecursiveRemoveConsistency()
   {
      // ... unless the mechanism is disabled in which case
      // we need to do it explicitly.
      if (!gROOT->MustClean())
         gROOT->GetListOfCleanups()->Remove(this);
   }

   void Add(TObject *obj)
   {
      obj->SetBit(kMustCleanup);
      auto hashValue = obj->Hash(); // This might/will take the ROOT lock.

      std::unique_lock<std::mutex> lock(fMutex);
      fCont.push_back(Value_t{hashValue, obj});
   }

   void RecursiveRemove(TObject *obj)
   {
      // Since we use std::list, a remove (from another thread)
      // would invalidate out iterator and taking the write lock
      // 'only' inside the loop would suspend this thread and lead
      // another reader or write go on; consequently we would need
      // to re-find the object we are wanting to remove.
      std::unique_lock<std::mutex> lock(fMutex);

      // std::cout << "Recursive Remove called for: " << obj << '\n';
      for (auto p = fCont.begin(); p != fCont.end(); ++p) {
         if (p->fObjectPtr == obj) {
            // std::cout << " Found object with hash = " << p->fRecordedHash << '\n';
            // std::cout << " Current hash = " << obj->Hash() << '\n';
            if (p->fRecordedHash == obj->Hash())
               fCont.erase(p);
            // else
            // std::cout << " Error: the recorded hash and the one returned by Hash are distinct.\n";
            break;
         }
      }
   }

   void SlowRemove(TObject *obj)
   {
      std::unique_lock<std::mutex> lock(fMutex);

      for (auto p = fCont.begin(); p != fCont.end(); ++p) {
         if (p->fObjectPtr == obj) {
            fCont.erase(p);
            break;
         }
      }
   }

   enum EResult {
      kInconsistent,
      kInconclusive,
      kConsistent
   };

   EResult CheckRecursiveRemove(TClass &classRef)
   {
      if (!classRef.HasDefaultConstructor() || classRef.Property() & kIsAbstract)
         return kInconclusive; // okay that's probably a false negative ...

      auto size = fCont.size();
      TObject *obj = (TObject *)classRef.DynamicCast(TObject::Class(), classRef.New(TClass::kDummyNew));
      if (!obj || (!gROOT->MustClean() && obj->TestBit(kIsReferenced) && obj->GetUniqueID() != 0)) {
         // Clean up is disable and the object is such that we wont be able to 'mark' it
         // as needing a clean up anyway, so we can not actually test it.
         return kInconclusive;
      }
      ROOT::Internal::SetRequireCleanup(*obj);
      Add(obj);
      delete obj;

      if (fCont.size() != size) {
         // std::cerr << "Error: old= " << size << " new=" << fCont.size() << '\n';
         // std::cerr << "Error " << classRef.GetName() <<
         //   " or one of its base classes override TObject::Hash but does not call TROOT::CallRecursiveRemoveIfNeeded
         //   in its destructor.\n";
         SlowRemove(obj);
         return kInconsistent;
      } else {
         return kConsistent;
      }
   }

   TClass *FindMissingRecursiveRemove(TClass &classRef)
   {

      if (classRef.HasLocalHashMember() && CheckRecursiveRemove(classRef) != kConsistent) {
         return &classRef;
      }

      for (auto base : ROOT::Detail::TRangeStaticCast<TBaseClass>(classRef.GetListOfBases())) {
         TClass *baseCl = base->GetClassPointer();
         TClass *res = FindMissingRecursiveRemove(*baseCl);
         if (res)
            return res;
      }
      return nullptr;
   }

   bool VerifyRecursiveRemove(const char *classname)
   {
      TClass *classPtr = TClass::GetClass(classname);
      if (classPtr)
         return VerifyRecursiveRemove(*classPtr);
      else
         return true;
   }

   EResult HasConsistentHashMember(TClass &classRef)
   {
      // Use except if the class is non-default/abstract and HasLocalHashMember.
      if (classRef.fRuntimeProperties) {
         // We already did this testing for this class.
         return classRef.HasConsistentHashMember() ? kConsistent : kInconsistent;
      }

      if (classRef.HasLocalHashMember())
         return CheckRecursiveRemove(classRef);

      EResult baseResult = kConsistent;
      for (auto base : ROOT::Detail::TRangeStaticCast<TBaseClass>(classRef.GetListOfBases())) {
         TClass *baseCl = base->GetClassPointer();

         if (baseCl->HasLocalHashMember() &&
            (!baseCl->HasDefaultConstructor() || baseCl->Property() & kIsAbstract))
         {
            // We won't be able to check the base class, we need to (try) to check
            // this class even-though it does not have a local HashMember.
            return CheckRecursiveRemove(classRef);
         }
         auto baseConsistency = HasConsistentHashMember(*baseCl);
         if (baseConsistency == kInconsistent) {
            baseResult = kInconsistent;
         } else if (baseConsistency == kInconclusive) {
            return CheckRecursiveRemove(classRef);
         }
      }
      return baseResult;
   }

   bool VerifyRecursiveRemove(TClass &classRef)
   {
      // If the class does not inherit from TObject, the setup is always 'correct'
      // (or more exactly does not matter).
      if (!classRef.IsTObject())
         return true;

      if (classRef.HasLocalHashMember() &&
          (!classRef.HasDefaultConstructor() || classRef.Property() & kIsAbstract))
         // We won't be able to check, so assume the worst but don't issue any
         // error message.
         return false;

      if (HasConsistentHashMember(classRef) != kConsistent) {
         TClass *failing = FindMissingRecursiveRemove(classRef);

         // Because ClassDefInline does not yet support class template on all platforms,
         // we have no ClassDef and thus can not get a good message from TObject::Error.
         constexpr const char *funcName = "ROOT::Internal::TCheckHashRecursiveRemoveConsistency::CheckRecursiveRemove";
         if (failing) {
            ::Error(funcName,
                    "The class %s overrides TObject::Hash but does not call TROOT::RecursiveRemove in its destructor (seen while checking %s).",
                    failing->GetName(),classRef.GetName());
         } else {
            ::Error(funcName, "The class %s "
                              "or one of its base classes override TObject::Hash but does not call "
                              "TROOT::CallRecursiveRemoveIfNeeded in its destructor.\n",
                    classRef.GetName());
         }
         return false;
      }
      return true;
   }

   static bool Check(TClass &classRef)
   {
      TCheckHashRecursiveRemoveConsistency checker;
      return checker.VerifyRecursiveRemove(classRef);
   }

   ClassDefInline(TCheckHashRecursiveRemoveConsistency, 0);
};

} // namespace Internal
} // namespace ROOT

#endif // ROOT__TCheckHashRecursiveRemoveConsistency
