/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RIoUring
#define ROOT_RIoUring

#include <cstdint>
#include <cstring>
#include <stdexcept>

#include <liburing.h>
#include <liburing/io_uring.h>

#include "TError.h"

namespace ROOT {
namespace Internal {

class RIoUring {
private:
   struct io_uring fRing;
   std::uint32_t fDepth = 0;

   static bool CheckIsAvailable() {
      try {
         RIoUring(1);
         return true;
      }
      catch (const std::runtime_error& err) {
         Warning("RIoUring", "io_uring is not available\n%s", err.what());
      }
      return false;
   }

public:
   // Create an io_uring instance. The ring selects an appropriate queue depth. which can be queried
   // afterwards using GetQueueDepth(). The depth is typically 1024 or lower. Throws an exception if
   // ring setup fails.
   RIoUring() {
      std::uint32_t queueDepth = 1024;
      int ret;
      while (true) {
         ret = io_uring_queue_init(queueDepth, &fRing, 0 /* no flags */);
         if (ret == 0) {
            fDepth = queueDepth;
            break; // ring setup succeeded
         }
         if (ret != -ENOMEM) {
            throw std::runtime_error("Error initializing io_uring: " + std::string(std::strerror(-ret)));
         }
         // try again with a smaller queue for ENOMEM
         queueDepth /= 2;
         if (queueDepth == 0) {
            throw std::runtime_error("Fatal Error: failed to allocate memory for the smallest possible "
               "io_uring instance. 'memlock' memory has been exhausted for this user");
         }
      }
   }

   // Create a io_uring instance that can hold at least `entriesHint` submission entries. The actual
   // queue depth is rounded up to the next power of 2. Throws an exception if ring setup fails.
   explicit RIoUring(std::uint32_t entriesHint) {
      struct io_uring_params params = {}; /* zero initialize param struct, no flags */
      int ret = io_uring_queue_init_params(entriesHint, &fRing, &params);
      if (ret != 0) {
         throw std::runtime_error("Error initializing io_uring: " + std::string(std::strerror(-ret)));
      }
      fDepth = params.sq_entries;
   }

   RIoUring(const RIoUring&) = delete;
   RIoUring& operator=(const RIoUring&) = delete;

   ~RIoUring() {
      // todo(max) try submitting any pending events before exiting
      io_uring_queue_exit(&fRing);
   }

   /// Check if io_uring is available on this system.
   static bool IsAvailable() {
      static const bool available = RIoUring::CheckIsAvailable();
      return available;
   }

   std::uint32_t GetQueueDepth() {
      return fDepth;
   }

   /// Access the raw io_uring instance.
   struct io_uring *GetRawRing() {
      return &fRing;
   }

   /// Basic read event composed of IO data and a target file descriptor.
   struct RReadEvent {
      /// The destination for reading
      void *fBuffer = nullptr;
      /// The file offset
      std::uint64_t fOffset = 0;
      /// The number of desired bytes
      std::size_t fSize = 0;
      /// The number of actually read bytes, set by the RIoUring instance
      std::size_t fOutBytes = 0;
      /// The file descriptor
      int fFileDes = -1;
   };

   /// Submit a number of read events and wait for completion. Events are submitted in batches if
   /// the number of events is larger than the submission queue depth.
   void SubmitReadsAndWait(RReadEvent* readEvents, unsigned int nReads) {
      unsigned int batch = 0;
      unsigned int batchSize = fDepth;
      unsigned int readPos = 0;

      while (readPos < nReads) {
         if (readPos + batchSize > nReads) {
            batchSize = nReads - readPos;
         }
         // prep reads
         struct io_uring_sqe *sqe;
         for (std::size_t i = readPos; i < readPos + batchSize; ++i) {
            sqe = io_uring_get_sqe(&fRing);
            if (!sqe) {
               throw std::runtime_error("batch " + std::to_string(batch) + ": "
                  + "get SQE failed for read request '" + std::to_string(i)
                  + "', error: " + std::string(strerror(errno)));
            }
            if (readEvents[i].fFileDes == -1) {
               throw std::runtime_error("batch " + std::to_string(batch) + ": "
                  + "bad fd (-1) for read request '" + std::to_string(i) + "'");
            }
            if (readEvents[i].fBuffer == nullptr) {
               throw std::runtime_error("batch " + std::to_string(batch) + ": "
                  + "null read buffer for read request '" + std::to_string(i) + "'");
            }
            io_uring_prep_read(sqe,
               readEvents[i].fFileDes,
               readEvents[i].fBuffer,
               readEvents[i].fSize,
               readEvents[i].fOffset
            );
            sqe->flags |= IOSQE_ASYNC; // maximize read event throughput
            sqe->user_data = i;
         }

         // todo(max) check for any difference between submit vs. submit and wait for large nReq
         int submitted = io_uring_submit_and_wait(&fRing, batchSize);
         if (submitted <= 0) {
            throw std::runtime_error("batch " + std::to_string(batch) + ": "
               "ring submit failed, error: " + std::string(strerror(errno)));
         }
         if (submitted != static_cast<int>(batchSize)) {
            throw std::runtime_error("ring submitted " + std::to_string(submitted) +
               " events but requested " + std::to_string(batchSize));
         }
         // reap reads
         struct io_uring_cqe *cqe;
         int ret;
         for (int i = 0; i < submitted; ++i) {
            ret = io_uring_wait_cqe(&fRing, &cqe);
            if (ret < 0) {
               throw std::runtime_error("wait cqe failed, error: " + std::string(std::strerror(-ret)));
            }
            auto index = reinterpret_cast<std::size_t>(io_uring_cqe_get_data(cqe));
            if (index >= nReads) {
               throw std::runtime_error("bad cqe user data: " + std::to_string(index));
            }
            if (cqe->res < 0) {
               throw std::runtime_error("batch " + std::to_string(batch) + ": "
                  + "read failed for ReadEvent[" + std::to_string(index) + "], "
                  "error: " + std::string(std::strerror(-cqe->res)));
            }
            readEvents[index].fOutBytes = static_cast<std::size_t>(cqe->res);
            io_uring_cqe_seen(&fRing, cqe);
         }
         readPos += batchSize;
         batch += 1;
      }
      return;
   }
};

} // namespace Internal
} // namespace ROOT

#endif
