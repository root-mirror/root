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
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RPageStorage.hxx>

#include <TError.h>

#include <algorithm>
#include <chrono>
#include <future>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <utility>

bool ROOT::Experimental::Detail::RClusterPool::RInFlightCluster::operator <(const RInFlightCluster &other) const
{
   if (fClusterId == other.fClusterId) {
      if (fColumns.size() == other.fColumns.size()) {
         for (auto itr1 = fColumns.begin(), itr2 = other.fColumns.begin(); itr1 != fColumns.end(); ++itr1, ++itr2) {
            if (*itr1 == *itr2)
               continue;
            return *itr1 < *itr2;
         }
         // *this == other
         return false;
      }
      return fColumns.size() < other.fColumns.size();
   }
   return fClusterId < other.fClusterId;
}

ROOT::Experimental::Detail::RClusterPool::RClusterPool(RPageSource &pageSource, unsigned int size)
   : fPageSource(pageSource)
   , fPool(size)
   , fThreadIo(&RClusterPool::ExecReadClusters, this)
   , fThreadUnzip(&RClusterPool::ExecUnzipClusters, this)
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
      fReadQueue.emplace(RReadItem());
      fCvHasReadWork.notify_one();
   }
   fThreadIo.join();

   {
      // Controlled shutdown of the unzip thread
      std::unique_lock<std::mutex> lock(fLockUnzipQueue);
      fUnzipQueue.emplace(RUnzipItem());
      fCvHasUnzipWork.notify_one();
   }
   fThreadUnzip.join();
}

void ROOT::Experimental::Detail::RClusterPool::ExecUnzipClusters()
{
   while (true) {
      std::vector<RUnzipItem> unzipItems;
      {
         std::unique_lock<std::mutex> lock(fLockUnzipQueue);
         fCvHasUnzipWork.wait(lock, [&]{ return !fUnzipQueue.empty(); });
         while (!fUnzipQueue.empty()) {
            unzipItems.emplace_back(std::move(fUnzipQueue.front()));
            fUnzipQueue.pop();
         }
      }

      for (auto &item : unzipItems) {
         if (!item.fCluster)
            return;

         fPageSource.UnzipCluster(item.fCluster.get());

         // Afterwards the GetCluster() method in the main thread can pick-up the cluster
         item.fPromise.set_value(std::move(item.fCluster));
      }
   } // while (true)
}

void ROOT::Experimental::Detail::RClusterPool::ExecReadClusters()
{
   while (true) {
      std::vector<RReadItem> readItems;
      {
         std::unique_lock<std::mutex> lock(fLockWorkQueue);
         fCvHasReadWork.wait(lock, [&]{ return !fReadQueue.empty(); });
         while (!fReadQueue.empty()) {
            readItems.emplace_back(std::move(fReadQueue.front()));
            fReadQueue.pop();
         }
      }

      // At this point, but not before, we can create an RPageSink with the same metadata held by fPageSource
      // Before this point, fPageSource still doesn't have all the header metadata so it's pointless.
      if (!fPageSink){

         // HARDCODED FILENAME FOR THE CACHED RNTUPLE
         std::string_view cachedntuplepath = "cachedntuple.root";

         fPageSink = RPageSink::Create(fPageSource.GetNTupleName(), cachedntuplepath);
         auto modelptr = fPageSource.GetDescriptor().GenerateModel()->Clone();
         fPageSink->Create(*modelptr);
      }

      // Track the number of entries seen so far
      // Needed in the call to fPageSink->CommitCluster that will be done in CacheCluster
      ClusterSize_t entriessofar{0};
      for (auto &item : readItems) {
         if (item.fClusterId == kInvalidDescriptorId){
            // Need to commit the cached dataset if present
            fPageSink->CommitDataset();
            return;
         }

         // TODO(jblomer): the page source needs to be capable of loading multiple clusters in one go
         auto cluster = fPageSource.LoadCluster(item.fClusterId, item.fColumns);

         // Cache the current cluster
         auto CacheCluster = [&](DescriptorId_t clusterId){
            const auto &clusterDesc = fPageSource.GetDescriptor().GetClusterDescriptor(clusterId);
            const auto clusterentries = clusterDesc.GetNEntries();

            // Traverse columns in cluster
            for (auto columnId: clusterDesc.GetColumnIds()){

               const auto &pageRange = clusterDesc.GetPageRange(columnId);
               std::uint32_t firstElementInPage = 0;


               RPageStorage::RSealedPage sealedPage;
               // Traverse pages in column
               for (const auto &pi : pageRange.fPageInfos) {

                  auto buffer = std::make_unique<unsigned char []>(pi.fLocator.fBytesOnStorage);
                  sealedPage.fBuffer = buffer.get();
                  fPageSource.LoadSealedPage(columnId, RClusterIndex(clusterId, firstElementInPage), sealedPage);
                  firstElementInPage += pi.fNElements;

                  // Commit page
                  fPageSink->CommitSealedPage(columnId, sealedPage);

               }

            }
            entriessofar += clusterentries;
            fPageSink->CommitCluster(entriessofar);

         };
         CacheCluster(cluster->GetId());

         // Meanwhile, the user might have requested clusters outside the look-ahead window, so that we don't
         // need the cluster anymore, in which case we simply discard it right away, before moving it to the pool
         bool discard = false;
         {
            std::unique_lock<std::mutex> lock(fLockWorkQueue);
            for (auto &inFlight : fInFlightClusters) {
               if (inFlight.fClusterId != item.fClusterId)
                  continue;
               discard = inFlight.fIsExpired;
               break;
            }
         }
         if (discard) {
            cluster.reset();
            item.fPromise.set_value(std::move(cluster));
         } else {
            // Hand-over the loaded cluster pages to the unzip thread
            std::unique_lock<std::mutex> lock(fLockUnzipQueue);
            fUnzipQueue.emplace(RUnzipItem{std::move(cluster), std::move(item.fPromise)});
            fCvHasUnzipWork.notify_one();
         }
      }
   } // while (true)
}

ROOT::Experimental::Detail::RCluster *
ROOT::Experimental::Detail::RClusterPool::FindInPool(DescriptorId_t clusterId) const
{
   for (const auto &cptr : fPool) {
      if (cptr && (cptr->GetId() == clusterId))
         return cptr.get();
   }
   return nullptr;
}

size_t ROOT::Experimental::Detail::RClusterPool::FindFreeSlot() const
{
   auto N = fPool.size();
   for (unsigned i = 0; i < N; ++i) {
      if (!fPool[i])
         return i;
   }

   R__ASSERT(false);
   return N;
}


namespace {

/// Helper class for the (cluster, column list) pairs that should be loaded in the background
class RProvides {
   using DescriptorId_t = ROOT::Experimental::DescriptorId_t;
   using ColumnSet_t = ROOT::Experimental::Detail::RPageSource::ColumnSet_t;

private:
   std::map<DescriptorId_t, ColumnSet_t> fMap;

public:
   void Insert(DescriptorId_t clusterId, const ColumnSet_t &columns)
   {
      fMap.emplace(clusterId, columns);
   }

   bool Contains(DescriptorId_t clusterId) {
      return fMap.count(clusterId) > 0;
   }

   void Erase(DescriptorId_t clusterId, const ColumnSet_t &columns)
   {
      auto itr = fMap.find(clusterId);
      if (itr == fMap.end())
         return;
      ColumnSet_t d;
      std::copy_if(itr->second.begin(), itr->second.end(), std::inserter(d, d.end()),
         [&columns] (DescriptorId_t needle) { return columns.count(needle) == 0; });
      if (d.empty()) {
         fMap.erase(itr);
      } else {
         itr->second = d;
      }
   }

   decltype(fMap)::iterator begin() { return fMap.begin(); }
   decltype(fMap)::iterator end() { return fMap.end(); }
};

} // anonymous namespace

ROOT::Experimental::Detail::RCluster *
ROOT::Experimental::Detail::RClusterPool::GetCluster(
   DescriptorId_t clusterId, const RPageSource::ColumnSet_t &columns)
{
   const auto &desc = fPageSource.GetDescriptor();

   // Determine previous cluster ids that we keep if they happen to be in the pool
   std::set<DescriptorId_t> keep;
   auto prev = clusterId;
   for (unsigned int i = 0; i < fWindowPre; ++i) {
      prev = desc.FindPrevClusterId(prev);
      if (prev == kInvalidDescriptorId)
         break;
      keep.insert(prev);
   }

   // Determine following cluster ids and the column ids that we want to make available
   RProvides provide;
   provide.Insert(clusterId, columns);
   auto next = clusterId;
   // TODO(jblomer): instead of a fixed-sized window, eventually we should determine the window size based on
   // a user-defined memory limit.  The size of the preloaded data can be determined at the beginning of
   // GetCluster from the descriptor and the current contents of fPool.
   for (unsigned int i = 1; i < fWindowPost; ++i) {
      next = desc.FindNextClusterId(next);
      if (next == kInvalidDescriptorId)
         break;
      provide.Insert(next, columns);
   }

   // Clear the cache from clusters not the in the look-ahead or the look-back window
   for (auto &cptr : fPool) {
      if (!cptr)
         continue;
      if (provide.Contains(cptr->GetId()) > 0)
         continue;
      if (keep.count(cptr->GetId()) > 0)
         continue;
      cptr.reset();
   }

   // Move clusters that meanwhile arrived into cache pool
   {
      // This lock is held during iteration over several data structures: the collection of in-flight clusters,
      // the current pool of cached clusters, and the set of cluster ids to be preloaded.
      // All three collections are expected to be small (certainly < 100, more likely < 10).  All operations
      // are non-blocking and moving around small items (pointers, ids, etc).  Thus the overall locking time should
      // still be reasonably small and the lock is rarely taken (usually once per cluster).
      std::lock_guard<std::mutex> lockGuard(fLockWorkQueue);

      for (auto itr = fInFlightClusters.begin(); itr != fInFlightClusters.end(); ) {
         R__ASSERT(itr->fFuture.valid());
         itr->fIsExpired = !provide.Contains(itr->fClusterId) && (keep.count(itr->fClusterId) == 0);

         if (itr->fFuture.wait_for(std::chrono::seconds(0)) != std::future_status::ready) {
            // Remove the set of columns that are already scheduled for being loaded
            provide.Erase(itr->fClusterId, itr->fColumns);
            ++itr;
            continue;
         }

         auto cptr = itr->fFuture.get();
         // If cptr is nullptr, the cluster expired previously and was released by the I/O thread
         if (!cptr || itr->fIsExpired) {
            cptr.reset();
            itr = fInFlightClusters.erase(itr);
            continue;
         }

         // We either put a fresh cluster into a free slot or we merge the cluster with an existing one
         auto existingCluster = FindInPool(cptr->GetId());
         if (existingCluster) {
            existingCluster->Adopt(std::move(*cptr));
         } else {
            auto idxFreeSlot = FindFreeSlot();
            fPool[idxFreeSlot] = std::move(cptr);
         }
         itr = fInFlightClusters.erase(itr);
      }

      // Determine clusters which get triggered for background loading
      for (auto &cptr : fPool) {
         if (!cptr)
            continue;
         provide.Erase(cptr->GetId(), cptr->GetAvailColumns());
      }

      // Update the work queue and the in-flight cluster list with new requests. We already hold the work queue
      // mutex
      // TODO(jblomer): we should ensure that clusterId is given first to the I/O thread.  That is usually the
      // case but it's not ensured by the code
      for (const auto &kv : provide) {
         R__ASSERT(!kv.second.empty());

         RReadItem readItem;
         readItem.fClusterId = kv.first;
         readItem.fColumns = kv.second;

         RInFlightCluster inFlightCluster;
         inFlightCluster.fClusterId = kv.first;
         inFlightCluster.fColumns = kv.second;
         inFlightCluster.fFuture = readItem.fPromise.get_future();
         fInFlightClusters.emplace_back(std::move(inFlightCluster));

         fReadQueue.emplace(std::move(readItem));
      }
      if (fReadQueue.size() > 0)
         fCvHasReadWork.notify_one();
   } // work queue lock guard

   return WaitFor(clusterId, columns);
}


ROOT::Experimental::Detail::RCluster *
ROOT::Experimental::Detail::RClusterPool::WaitFor(
   DescriptorId_t clusterId, const RPageSource::ColumnSet_t &columns)
{
   while (true) {
      // Fast exit: the cluster happens to be already present in the cache pool
      auto result = FindInPool(clusterId);
      if (result) {
         bool hasMissingColumn = false;
         for (auto cid : columns) {
            if (result->ContainsColumn(cid))
               continue;

            hasMissingColumn = true;
            break;
         }
         if (!hasMissingColumn)
            return result;
      }

      // Otherwise the missing data must have been triggered for loading by now, so block and wait
      decltype(fInFlightClusters)::iterator itr;
      {
         std::lock_guard<std::mutex> lockGuardInFlightClusters(fLockWorkQueue);
         itr = fInFlightClusters.begin();
         for (; itr != fInFlightClusters.end(); ++itr) {
            if (itr->fClusterId == clusterId)
               break;
         }
         R__ASSERT(itr != fInFlightClusters.end());
         // Note that the fInFlightClusters is accessed concurrently only by the I/O thread.  The I/O thread
         // never changes the structure of the in-flight clusters array (it does not add, remove, or swap elements).
         // Therefore, it is safe to access the element pointed to by itr here even after fLockWorkQueue
         // is released.  We need to release the lock before potentially blocking on the cluster future.
      }

      auto cptr = itr->fFuture.get();
      if (result) {
         result->Adopt(std::move(*cptr));
      } else {
         auto idxFreeSlot = FindFreeSlot();
         fPool[idxFreeSlot] = std::move(cptr);
      }

      std::lock_guard<std::mutex> lockGuardInFlightClusters(fLockWorkQueue);
      fInFlightClusters.erase(itr);
   }
}
