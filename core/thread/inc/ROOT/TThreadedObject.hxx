// @(#)root/thread:$Id$
// Author: Danilo Piparo, CERN  11/2/2016

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TThreadedObject
#define ROOT_TThreadedObject

#include "ROOT/TSpinMutex.hxx"
#include "TDirectory.h"
#include "TError.h"
#include "TList.h"
#include "TROOT.h"


#include <algorithm>
#include <exception>
#include <deque>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

class TH1;

namespace ROOT {

   /**
    * \class ROOT::TNumSlots
    * \brief Defines the number of threads in some of ROOT's interfaces.
    */
   struct TNumSlots {
      unsigned int fVal; // number of slots
      friend bool operator==(TNumSlots lhs, TNumSlots rhs) { return lhs.fVal == rhs.fVal; }
      friend bool operator!=(TNumSlots lhs, TNumSlots rhs) { return lhs.fVal != rhs.fVal; }
   };

   namespace Internal {

      namespace TThreadedObjectUtils {


         template<typename T, bool ISHISTO = std::is_base_of<TH1,T>::value>
         struct Detacher{
            static T* Detach(T* obj) {
               return obj;
            }
         };

         template<typename T>
         struct Detacher<T, true>{
            static T* Detach(T* obj) {
               obj->SetDirectory(nullptr);
               obj->ResetBit(kMustCleanup);
               return obj;
            }
         };

         /// Return a copy of the object or a "Clone" if the copy constructor is not implemented.
         template<class T, bool isCopyConstructible = std::is_copy_constructible<T>::value>
         struct Cloner {
            static T *Clone(const T *obj, TDirectory* d = nullptr) {
               T* clone;
               if (d){
                  TDirectory::TContext ctxt(d);
                  clone = new T(*obj);
               } else {
                  clone = new T(*obj);
               }
               return Detacher<T>::Detach(clone);
            }
         };

         template<class T>
         struct Cloner<T, false> {
            static T *Clone(const T *obj, TDirectory* d = nullptr) {
               T* clone;
               if (d){
                  TDirectory::TContext ctxt(d);
                  clone = (T*)obj->Clone();
               } else {
                  clone = (T*)obj->Clone();
               }
               return clone;
            }
         };

         template <class T, bool ISHISTO = std::is_base_of<TH1, T>::value>
         struct DirCreator {
            static TDirectory *Create()
            {
               static unsigned dirCounter = 0;
               const std::string dirName = "__TThreaded_dir_" + std::to_string(dirCounter++) + "_";
               return gROOT->mkdir(dirName.c_str());
            }
         };

         template <class T>
         struct DirCreator<T, true> {
            static TDirectory *Create() { return nullptr; }
         };

      } // End of namespace TThreadedObjectUtils
   } // End of namespace Internal

   namespace TThreadedObjectUtils {

      template<class T>
      using MergeFunctionType = std::function<void(std::shared_ptr<T>, std::vector<std::shared_ptr<T>>&)>;

      /// Merge TObjects
      template<class T>
      void MergeTObjects(std::shared_ptr<T> target, std::vector<std::shared_ptr<T>> &objs)
      {
         if (!target) return;
         TList objTList;
         // Cannot do better than this
         for (auto obj : objs) {
            if (obj && obj != target) objTList.Add(obj.get());
         }
         target->Merge(&objTList);
      }
   } // end of namespace TThreadedObjectUtils

   /**
    * \class ROOT::TThreadedObject
    * \brief A wrapper to make object instances thread private, lazily.
    * \tparam T Class of the object to be made thread private (e.g. TH1F)
    * \ingroup Multicore
    *
    * A wrapper which makes objects thread private. The methods of the underlying
    * object can be invoked via the arrow operator. The object is created in
    * a specific thread lazily, i.e. upon invocation of one of its methods.
    * The correct object pointer from within a particular thread can be accessed
    * with the overloaded arrow operator or with the Get method.
    * In case an elaborate thread management is in place, e.g. in presence of
    * stream of operations or "processing slots", it is also possible to
    * manually select the correct object pointer explicitly.
    */
   template<class T>
   class TThreadedObject {
   public:
      /// The initial number of empty processing slots that a TThreadedObject is constructed with by default.
      /// Deprecated: TThreadedObject grows as more slots are required.
      static constexpr const TNumSlots fgMaxSlots{64};

      TThreadedObject(const TThreadedObject&) = delete;

      /// Construct the TThreadedObject with initSlots empty slots and the "model" of the thread private objects.
      /// \param initSlots Set the initial number of slots of the TThreadedObject.
      /// \tparam ARGS Arguments of the constructor of T
      ///
      /// This form of the constructor is useful to manually pre-set the content of a given number of slots
      /// when used in combination with TThreadedObject::SetAtSlot().
      template <class... ARGS>
      TThreadedObject(TNumSlots initSlots, ARGS &&... args) : fIsMerged(false)
      {
         const auto nSlots = initSlots.fVal;
         fObjPointers.resize(nSlots);

         // create at least one directory (we need it for fModel), plus others as needed by the size of fObjPointers
         fDirectories.emplace_back(Internal::TThreadedObjectUtils::DirCreator<T>::Create());
         for (auto i = 1u; i < nSlots; ++i)
            fDirectories.emplace_back(Internal::TThreadedObjectUtils::DirCreator<T>::Create());

         TDirectory::TContext ctxt(fDirectories[0]);
         fModel.reset(Internal::TThreadedObjectUtils::Detacher<T>::Detach(new T(std::forward<ARGS>(args)...)));
      }

      /// Construct the TThreadedObject and the "model" of the thread private objects.
      /// \tparam ARGS Arguments of the constructor of T
      template<class ...ARGS>
      TThreadedObject(ARGS&&... args) : TThreadedObject(fgMaxSlots, args...) { }

      /// Return the number of currently available slot.
      ///
      /// The method is safe to call concurrently to other TThreadedObject methods.
      /// Note that slots could be available but contain no data (i.e. a nullptr) if
      /// they have not been used yet.
      unsigned GetNSlots() const
      {
         std::lock_guard<ROOT::TSpinMutex> lg(fSpinMutex);
         return fObjPointers.size();
      }

      /// Access a particular processing slot.
      ///
      /// This method is thread-safe as long as concurrent calls request different slots (i.e. pass a different
      /// argument) and no thread accesses slot `i` via the arrow operator, so mixing usage of GetAtSlot
      /// with usage of the arrow operator can be dangerous.
      std::shared_ptr<T> GetAtSlot(unsigned i)
      {
         std::size_t nAvailableSlots;
         {
            // fObjPointers can grow due to a concurrent operation on this TThreadedObject, need to lock
            std::lock_guard<ROOT::TSpinMutex> lg(fSpinMutex);
            nAvailableSlots = fObjPointers.size();
         }

         if (i >= nAvailableSlots) {
            Warning("TThreadedObject::GetAtSlot", "This slot does not exist.");
            return nullptr;
         }

         auto &objPointer = fObjPointers[i];
         if (!objPointer)
            objPointer.reset(Internal::TThreadedObjectUtils::Cloner<T>::Clone(fModel.get(), fDirectories[i]));
         return objPointer;
      }

      /// Set the value of a particular slot.
      ///
      /// This method is thread-safe as long as concurrent calls access different slots (i.e. pass a different
      /// argument) and no thread accesses slot `i` via the arrow operator, so mixing usage of SetAtSlot
      /// with usage of the arrow operator can be dangerous.
      void SetAtSlot(unsigned i, std::shared_ptr<T> v)
      {
         std::size_t nAvailableSlots;
         {
            // fObjPointers can grow due to a concurrent operation on this TThreadedObject, need to lock
            std::lock_guard<ROOT::TSpinMutex> lg(fSpinMutex);
            nAvailableSlots = fObjPointers.size();
         }

         if (i >= nAvailableSlots) {
            Warning("TThreadedObject::SetAtSlot", "This slot does not exist, doing nothing.");
            return;
         }

         fObjPointers[i] = v;
      }

      /// Access a particular slot which corresponds to a single thread.
      /// This is in general faster than the GetAtSlot method but it is
      /// responsibility of the caller to make sure that the slot exists
      /// and to check that the contained object is initialized (and not
      /// a nullptr).
      std::shared_ptr<T> GetAtSlotUnchecked(unsigned i) const
      {
         return fObjPointers[i];
      }

      /// Access a particular slot which corresponds to a single thread.
      /// This overload is faster than the GetAtSlotUnchecked method but
      /// the caller is responsible to make sure that the slot exists, to
      /// check that the contained object is initialized and that the returned
      /// pointer will not outlive the TThreadedObject that returned it, which
      /// maintains ownership of the actual object.
      T* GetAtSlotRaw(unsigned i) const
      {
         return fObjPointers[i].get();
      }

      /// Access the pointer corresponding to the current slot. This method is
      /// not adequate for being called inside tight loops as it implies a
      /// lookup in a mapping between the threadIDs and the slot indices.
      /// A good practice consists in copying the pointer onto the stack and
      /// proceed with the loop as shown in this work item (psudo-code) which
      /// will be sent to different threads:
      /// ~~~{.cpp}
      /// auto workItem = [](){
      ///    auto objPtr = tthreadedObject.Get();
      ///    for (auto i : ROOT::TSeqI(1000)) {
      ///       // tthreadedObject->FastMethod(i); // don't do this! Inefficient!
      ///       objPtr->FastMethod(i);
      ///    }
      /// }
      /// ~~~
      std::shared_ptr<T> Get()
      {
         return GetAtSlot(GetThisSlotNumber());
      }

      /// Access the wrapped object and allow to call its methods.
      T *operator->()
      {
         return Get().get();
      }

      /// Merge all the thread private objects. Can be called once: it does not
      /// create any new object but destroys the present bookkeping collapsing
      /// all objects into the one at slot 0.
      std::shared_ptr<T> Merge(TThreadedObjectUtils::MergeFunctionType<T> mergeFunction = TThreadedObjectUtils::MergeTObjects<T>)
      {
         // We do not return if we already merged.
         if (fIsMerged) {
            Warning("TThreadedObject::Merge", "This object was already merged. Returning the previous result.");
            return fObjPointers[0];
         }
         // need to convert to std::vector because historically mergeFunction requires a vector
         auto vecOfObjPtrs = std::vector<std::shared_ptr<T>>(fObjPointers.begin(), fObjPointers.end());
         mergeFunction(fObjPointers[0], vecOfObjPtrs);
         fIsMerged = true;
         return fObjPointers[0];
      }

      /// Merge all the thread private objects. Can be called many times. It
      /// does create a new instance of class T to represent the "Sum" object.
      /// This method is not thread safe: correct or acceptable behaviours
      /// depend on the nature of T and of the merging function.
      std::unique_ptr<T> SnapshotMerge(TThreadedObjectUtils::MergeFunctionType<T> mergeFunction = TThreadedObjectUtils::MergeTObjects<T>)
      {
         if (fIsMerged) {
            Warning("TThreadedObject::SnapshotMerge", "This object was already merged. Returning the previous result.");
            return std::unique_ptr<T>(Internal::TThreadedObjectUtils::Cloner<T>::Clone(fObjPointers[0].get()));
         }
         auto targetPtr = Internal::TThreadedObjectUtils::Cloner<T>::Clone(fModel.get());
         std::shared_ptr<T> targetPtrShared(targetPtr, [](T *) {});
         // need to convert to std::vector because historically mergeFunction requires a vector
         auto vecOfObjPtrs = std::vector<std::shared_ptr<T>>(fObjPointers.begin(), fObjPointers.end());
         mergeFunction(targetPtrShared, vecOfObjPtrs);
         return std::unique_ptr<T>(targetPtr);
      }

   private:
      std::unique_ptr<T> fModel;                         ///< Use to store a "model" of the object
      // std::deque's guarantee that references to the elements are not invalidated when appending new slots
      std::deque<std::shared_ptr<T>> fObjPointers;       ///< An object pointer per slot
      // If the object is a histogram, we also create dummy directories that the histogram associates with
      // so we do not pollute gDirectory
      std::deque<TDirectory*> fDirectories;              ///< A TDirectory per slot
      std::map<std::thread::id, unsigned> fThrIDSlotMap; ///< A mapping between the thread IDs and the slots
      mutable ROOT::TSpinMutex fSpinMutex;               ///< Protects concurrent access to fThrIDSlotMap, fObjPointers
      bool fIsMerged : 1;                                ///< Remember if the objects have been merged already

      /// Get the slot number for this threadID, make a slot if needed
      unsigned GetThisSlotNumber()
      {
         const auto thisThreadID = std::this_thread::get_id();
         std::lock_guard<ROOT::TSpinMutex> lg(fSpinMutex);
         const auto thisSlotNumIt = fThrIDSlotMap.find(thisThreadID);
         if (thisSlotNumIt != fThrIDSlotMap.end())
            return thisSlotNumIt->second;
         const auto newIndex = fThrIDSlotMap.size();
         fThrIDSlotMap[thisThreadID] = newIndex;
         R__ASSERT(newIndex <= fObjPointers.size() && "This should never happen, we should create new slots as needed");
         if (newIndex == fObjPointers.size()) {
            fDirectories.emplace_back(Internal::TThreadedObjectUtils::DirCreator<T>::Create());
            fObjPointers.emplace_back(nullptr);
         }
         return newIndex;
      }
   };

   template<class T>
   constexpr const TNumSlots TThreadedObject<T>::fgMaxSlots;

} // End ROOT namespace

////////////////////////////////////////////////////////////////////////////////
/// Print a TThreadedObject at the prompt:

namespace cling {
   template<class T>
   std::string printValue(ROOT::TThreadedObject<T> *val)
   {
      auto model = ((std::unique_ptr<T>*)(val))->get();
      std::ostringstream ret;
      ret << "A wrapper to make object instances thread private, lazily. "
          << "The model which is replicated is " << printValue(model);
      return ret.str();
   }
}


#endif
