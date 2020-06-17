/// \file RClusterPool.cxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2020-03-11
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RClusterPool.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RPageStorage.hxx>

#include <TError.h>

#include <algorithm>
#include <chrono>
#include <future>
#include <iostream>
#include <memory>
#include <mutex>
#include <set>
#include <utility>

ROOT::Experimental::Detail::RClusterPool::RClusterPool(RPageSource *pageSource, unsigned int size)
   : fPageSource(pageSource)
   , fPool(size, nullptr)
   , fThreadIo(&RClusterPool::ExecLoadClusters, this)
{
   R__ASSERT(size > 0);
   fWindowPre = 0;
   fWindowPost = size;
   // Large pools maintain a small look-back window together with the large look-ahead window
   while ((1u << fWindowPre) < (fWindowPost - (fWindowPre + 1))) {
      fWindowPre++;
      fWindowPost--;
   }

}

ROOT::Experimental::Detail::RClusterPool::~RClusterPool()
{
   {
      // Controlled shutdown of the I/O thread
      std::unique_lock<std::mutex> lock(fLockWorkQueue);
      fWorkQueue.emplace(RWorkItem());
      fCvHasWork.notify_one();
   }
   fThreadIo.join();
}

void ROOT::Experimental::Detail::RClusterPool::ExecLoadClusters()
{
   while (true) {
      std::vector<RWorkItem> workItems;
      {
         std::unique_lock<std::mutex> lock(fLockWorkQueue);
         fCvHasWork.wait(lock, [&]{ return !fWorkQueue.empty(); });
         while (!fWorkQueue.empty()) {
            workItems.emplace_back(std::move(fWorkQueue.front()));
            fWorkQueue.pop();
         }
      }

      for (auto &item : workItems) {
         if (item.fClusterId == kInvalidDescriptorId)
            return;

         // TODO(jblomer): the page source needs to be capable of loading multiple clusters in one go
         auto cluster = fPageSource->LoadCluster(item.fClusterId);

         // Meanwhile, the user might have requested clusters outside the look-ahead window, so that we don't
         // need the cluster anymore, in which case we simply discard it before moving it to the pool
         bool discard = false;
         {
            std::lock_guard<std::mutex> lockGuardInFlightClusters(fLockInFlightClusters);
            for (auto &inFlight : fInFlightClusters) {
               if (inFlight.fClusterId != item.fClusterId)
                  continue;
               discard = inFlight.fIsExpired;
               break;
            }
         }
         if (discard)
            cluster.reset();

         item.fPromise.set_value(std::move(cluster));
      }
   } // while (true)
}

std::shared_ptr<ROOT::Experimental::Detail::RCluster>
ROOT::Experimental::Detail::RClusterPool::FindInPool(DescriptorId_t clusterId)
{
   for (const auto &cptr : fPool) {
      if (cptr && (cptr->GetId() == clusterId))
         return cptr;
   }
   return nullptr;
}

size_t ROOT::Experimental::Detail::RClusterPool::FindFreeSlot() {
   auto N = fPool.size();
   for (unsigned i = 0; i < N; ++i) {
      if (!fPool[i])
         return i;
   }

   R__ASSERT(false);
   return N;
}

std::shared_ptr<ROOT::Experimental::Detail::RCluster>
ROOT::Experimental::Detail::RClusterPool::GetCluster(ROOT::Experimental::DescriptorId_t clusterId)
{
   const auto &desc = fPageSource->GetDescriptor();

   // Determine previous cluster ids that we keep if they happen to be in the pool
   std::set<DescriptorId_t> cidKeep;
   auto prev = clusterId;
   for (unsigned int i = 0; i < fWindowPre; ++i) {
      prev = desc.FindPrevClusterId(prev);
      if (prev == kInvalidDescriptorId)
         break;
      cidKeep.insert(prev);
   }

   // Determine following cluster ids that we want to make available
   std::set<DescriptorId_t> cidProvide{clusterId};
   auto next = clusterId;
   for (unsigned int i = 0; i < fWindowPost - 1; ++i) {
      next = desc.FindNextClusterId(next);
      if (next == kInvalidDescriptorId)
         break;
      cidProvide.insert(next);
   }

   // Clear the cache from clusters not the in the look-ahead or the look-back window
   for (auto &cptr : fPool) {
      if (!cptr)
         continue;
      if (cidProvide.count(cptr->GetId()) > 0)
         continue;
      if (cidKeep.count(cptr->GetId()) > 0)
         continue;
      cptr.reset();
   }

   // Move clusters that meanwhile arrived into cache pool
   {
      std::lock_guard<std::mutex> lockGuardInFlightClusters(fLockInFlightClusters);
      for (auto itr = fInFlightClusters.begin(); itr != fInFlightClusters.end(); ) {
         R__ASSERT(itr->fFuture.valid());
         itr->fIsExpired = (cidProvide.count(itr->fClusterId) == 0) && (cidKeep.count(itr->fClusterId) == 0);

         if (itr->fFuture.wait_for(std::chrono::seconds(0)) != std::future_status::ready) {
            cidProvide.erase(itr->fClusterId);
            ++itr;
            continue;
         }

         auto cptr = itr->fFuture.get();
         // If cptr is nullptr, the cluster expired and was released by the I/O thread
         if (cptr) {
            if (itr->fIsExpired) {
               cptr.reset();
            } else {
               cidProvide.erase(cptr->GetId());
               auto idxFreeSlot = FindFreeSlot();
               fPool[idxFreeSlot] = std::move(cptr);
            }
         }
         itr = fInFlightClusters.erase(itr);
      }

      // Deterime clusters which get triggered for background loading
      for (auto &cptr : fPool) {
         if (cptr)
            cidProvide.erase(cptr->GetId());
      }

      // Update the work queue and the in-flight cluster list with new requests
      // TODO(jblomer): we should ensure that clusterId is given first to the I/O thread.  That is usually the
      // case but it's not ensured by the code
      std::unique_lock<std::mutex> lockWorkQueue(fLockWorkQueue);
      for (auto id : cidProvide) {
         RWorkItem workItem;
         workItem.fClusterId = id;

         RInFlightCluster inFlightCluster;
         inFlightCluster.fClusterId = id;
         inFlightCluster.fFuture = workItem.fPromise.get_future();
         fInFlightClusters.emplace_back(std::move(inFlightCluster));

         fWorkQueue.emplace(std::move(workItem));
         if (fWorkQueue.size() == 1)
            fCvHasWork.notify_one();
      }
      if (fWorkQueue.size() > 0)
         fCvHasWork.notify_one();
   } // work queue and in-flight clusters lock guards

   // Fast exit: the cluster happens to be already present in the cache pool
   auto result = FindInPool(clusterId);
   if (result)
      return result;

   // Otherwise it must have been triggered for loading by now, so block and wait
   auto slot = FindFreeSlot();
   auto itr = fInFlightClusters.begin();
   for (; itr != fInFlightClusters.end(); ++itr) {
      if (itr->fClusterId != clusterId)
         continue;
      fPool[slot] = itr->fFuture.get();
      break;
   }
   R__ASSERT(itr != fInFlightClusters.end());

   std::lock_guard<std::mutex> lockGuardInFlightClusters(fLockInFlightClusters);
   fInFlightClusters.erase(itr);

   return fPool[slot];
}
