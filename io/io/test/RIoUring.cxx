#include "io_test.hxx"

#include "ROOT/RIoUring.hxx"
#include "ROOT/RRawFileUnix.hxx"

using RIoUring = ROOT::Internal::RIoUring;
using RIOVec = RRawFile::RIOVec;
using RRawFileUnix = ROOT::Internal::RRawFileUnix;

namespace {

std::vector<RIOVec> make_iovecs(int n, unsigned int fileSize) {
   std::vector<RIOVec> iovecs;
   for (int i = 0; i < n; ++i) {
      RIOVec io;
      io.fBuffer = malloc(4096 * 9);
      io.fOffset = std::rand() % fileSize;
      io.fSize = 4096 * 8;
      iovecs.push_back(io);
   }
   return iovecs;
}

} // anonymous namespace

TEST(RRawFileUnix, ReadV)
{
   auto file = "test_uring_readv";
   auto filesize = 2 << 20;
   FileRaii fileGuard(file, std::string(filesize, 'a')); // ~2MB
   auto f = RRawFileUnix::Create(file);

   auto nReq = 2000; // demo submission batching, uring size is usually around 1024

   auto iovecs = make_iovecs(nReq, filesize);
   f->ReadV(iovecs.data(), nReq);

   for (auto iovec: iovecs) {
      for (std::size_t i = 0; i < iovec.fOutBytes; ++i) {
         EXPECT_EQ('a', ((unsigned char*)iovec.fBuffer)[i]);
      }
      free(iovec.fBuffer);
   }
}

TEST(RawUring, NopRoundTrip)
{
   struct io_uring ring;
   int ret = io_uring_queue_init(
      4 /* queue depth */,
      &ring,
      0 /* no setup flags */
   );

   // ring setup succeeded
   ASSERT_EQ(ret, 0);

   // can make sqes
   struct io_uring_sqe *sqe;
   sqe = io_uring_get_sqe(&ring);
   ASSERT_NE(sqe, (io_uring_sqe*) NULL);

   // can submit sqes to the ring
   io_uring_prep_nop(sqe);
   unsigned long udata = 42;
   io_uring_sqe_set_data(sqe, (void*) udata);
   ret = io_uring_submit(&ring);
   // submitted a single sqe
   ASSERT_EQ(ret, 1);

   // successful wait for cqe
   struct io_uring_cqe *cqe;
   ret = io_uring_wait_cqe(&ring, &cqe);
   ASSERT_EQ(ret, 0);
   // cqe userdata matches sqe userdata
   ASSERT_EQ((unsigned long) io_uring_cqe_get_data(cqe), udata);

   io_uring_cqe_seen(&ring, cqe);
   io_uring_queue_exit(&ring);
}

TEST(RawUring, FileRegistration)
{
   auto file = "test_uring_readv";
   auto filesize = 2 << 20;
   FileRaii fileGuard(file, std::string(filesize, 'a')); // ~2MB
   RRawFileUnix f(file, RRawFile::ROptions());
   // files are opened lazily, force file open via GetSize
   auto size = f.GetSize();

   unsigned int nReads = 100;
   auto iovecs = make_iovecs(nReads, size);

   RIoUring ring(nReads);
   EXPECT_EQ(ring.GetQueueDepth(), 128); // queue depth rounds up to next power of 2

   auto r = ring.GetRawRing();
   auto fd = f.GetFd();
   io_uring_register_files(r, &fd, 1);
   auto fixed_file = 0; // 1st entry in fixed file array => offset 0
   {
      struct io_uring_sqe *sqe;
      for (std::size_t i = 0; i < nReads; ++i) {
         sqe = io_uring_get_sqe(r);
         io_uring_prep_read(sqe,
            fixed_file,
            iovecs[i].fBuffer,
            iovecs[i].fSize,
            iovecs[i].fOffset
         );
         sqe->flags |= IOSQE_FIXED_FILE;
         sqe->user_data = i;
      }
      int submitted = io_uring_submit_and_wait(r, nReads);
      if (submitted <= 0) {
         throw std::runtime_error("ring submit failed, error: " + std::string(strerror(errno)));
      }
      // reap reads
      struct io_uring_cqe *cqe;
      int ret;
      for (int i = 0; i < submitted; ++i) {
         ret = io_uring_wait_cqe(r, &cqe);
         if (ret < 0) {
            throw std::runtime_error("wait cqe failed, error: " + std::string(std::strerror(-ret)));
         }
         auto index = reinterpret_cast<std::size_t>(io_uring_cqe_get_data(cqe));
         if (index >= nReads) {
            throw std::runtime_error("bad cqe user data: " + std::to_string(index));
         }
         if (cqe->res < 0) {
            throw std::runtime_error("read failed for ReadEvent[" + std::to_string(index) + "], "
               "error: " + std::string(std::strerror(-cqe->res)));
         }
         iovecs[index].fOutBytes = static_cast<std::size_t>(cqe->res);
         io_uring_cqe_seen(r, cqe);
      }
   }
   for (auto iovec: iovecs) {
      for (std::size_t i = 0; i < iovec.fOutBytes; ++i) {
         EXPECT_EQ('a', ((unsigned char*)iovec.fBuffer)[i]);
      }
      free(iovec.fBuffer);
   }
}
