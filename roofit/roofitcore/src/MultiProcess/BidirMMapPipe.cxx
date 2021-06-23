/** @file BidirMMapPipe.cxx
 *
 * implementation of BidirMMapPipe, a class which forks off a child process
 * and serves as communications channel between parent and child
 *
 * @author Manuel Schiller <manuel.schiller@nikhef.nl>
 * @author Patrick Bos <p.bos@esciencecenter.nl>
 * @author Inti Pelupessy <i.pelupessy@esciencecenter.nl>
 * @date 2013-2018
 */
#ifndef _WIN32
#include <map>
#include <cerrno>
#include <limits>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <exception>

#include <poll.h>
#include <fcntl.h>
#include <signal.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/socket.h>

#include <MultiProcess/BidirMMapPipe.h>
#include <sstream>

#define BEGIN_NAMESPACE_ROOFIT namespace RooFit {
#define END_NAMESPACE_ROOFIT }

BEGIN_NAMESPACE_ROOFIT

/// namespace for implementation details of BidirMMapPipe
namespace BidirMMapPipe_impl {
  // static function:
  int BidirMMapPipeException::dostrerror_r(int err, char* buf, std::size_t sz,
                                           int (*f)(int, char*, std::size_t)) {
    return f(err, buf, sz);
  }
  
  const char* BidirMMapPipeException::what() const noexcept {
    return m_buf;
  }
  
    BidirMMapPipeException::BidirMMapPipeException(const char* msg, int err)
    {
        std::size_t msgsz = std::strlen(msg);
        if (msgsz) {
            msgsz = std::min(msgsz, std::size_t(s_sz));
            std::copy(msg, msg + msgsz, m_buf);
            if (msgsz < s_sz) { m_buf[msgsz] = ':'; ++msgsz; }
            if (msgsz < s_sz) { m_buf[msgsz] = ' '; ++msgsz; }
        }
        if (msgsz < s_sz) {
            // UGLY: GNU and POSIX cannot agree on prototype and behaviour, so
            // have to sort it out with overloads
            dostrerror_r(err, &m_buf[msgsz], s_sz - msgsz, ::strerror_r);
        }
        m_buf[s_sz - 1] = 0; // enforce zero-termination
    }

    int BidirMMapPipeException::dostrerror_r(int err, char* buf,
            std::size_t sz, char* (*f)(int, char*, std::size_t))
    {
        buf[0] = 0;
        char *tmp = f(err, buf, sz);
        if (tmp && tmp != buf) {
            std::strncpy(buf, tmp, sz);
            buf[sz - 1] = 0;
            if (std::strlen(tmp) > sz - 1) return ERANGE;
        }
        return 0;
    }

    /** @brief class representing the header structure in an mmapped page
     *
     * @author Manuel Schiller <manuel.schiller@nikhef.nl>
     * @date 2013-07-07
     *
     * contains a field to put pages into a linked list, a field for the size
     * of the data being transmitted, and a field for the position until which
     * the data has been read
     */
    class Page
    {
        private:
            // use as small a data type as possible to maximise payload area
            // of pages
            short m_next;               ///< next page in list (in pagesizes)
            unsigned short m_size;      ///< size of payload (in bytes)
            unsigned short m_pos;       ///< index of next byte in payload area
            /// copy construction forbidden
            Page(const Page&) {}
            /// assigment forbidden
            Page& operator=(const Page&) = delete;
        public:
            /// constructor
            Page() : m_next(0), m_size(0), m_pos(0)
            {
                // check that short is big enough - must be done at runtime
                // because the page size is not known until runtime
                assert(std::numeric_limits<unsigned short>::max() >=
                        PageChunk::pagesize());
            }
            /// set pointer to next page
            void setNext(const Page* p);
            /// return pointer to next page
            Page* next() const;
            /// return reference to size field
            unsigned short& size() { return m_size; }
            /// return size (of payload data)
            unsigned size() const { return m_size; }
            /// return reference to position field
            unsigned short& pos() { return m_pos; }
            /// return position
            unsigned pos() const { return m_pos; }
            /// return pointer to first byte in payload data area of page
            inline unsigned char* begin() const
            { return reinterpret_cast<unsigned char*>(const_cast<Page*>(this))
                + sizeof(Page); }
            /// return pointer to first byte in payload data area of page
            inline unsigned char* end() const
            { return reinterpret_cast<unsigned char*>(const_cast<Page*>(this))
                + PageChunk::pagesize(); }
            /// return the capacity of the page
            static unsigned capacity()
            { return PageChunk::pagesize() - sizeof(Page); }
            /// true if page empty
            bool empty() const { return !m_size; }
            /// true if page partially filled
            bool filled() const { return !empty(); }
            /// free space left (to be written to)
            unsigned free() const { return capacity() - m_size; }
            /// bytes remaining to be read
            unsigned remaining() const { return m_size - m_pos; }
            /// true if page completely full
            bool full() const { return !free(); }
    };

    void Page::setNext(const Page* p)
    {
        if (!p) {
            m_next = 0;
        } else {
            const char* p1 = reinterpret_cast<char*>(this);
            const char* p2 = reinterpret_cast<const char*>(p);
            std::ptrdiff_t tmp = p2 - p1;
            // difference must be divisible by page size
            assert(!(tmp % PageChunk::pagesize()));
            tmp /= static_cast<std::ptrdiff_t>(PageChunk::pagesize());
            m_next = tmp;
            // no truncation when saving in a short
            assert(m_next == tmp);
            // final check: next() must return p
            assert(next() == p);
        }
    }

    Page* Page::next() const
    {
        if (!m_next) return 0;
        char* ptmp = reinterpret_cast<char*>(const_cast<Page*>(this));
        ptmp += std::ptrdiff_t(m_next) * PageChunk::pagesize();
        return reinterpret_cast<Page*>(ptmp);
    }

    /** @brief class representing a page pool
     *
     * @author Manuel Schiller <manuel.schiller@nikhef.nl>
     * @date 2013-07-24
     *
     * pool of mmapped pages (on systems which support it, on all others, the
     * functionality is emulated with dynamically allocated memory)
     *
     * in most operating systems there is a limit to how many mappings any one
     * process is allowed to request; for this reason, we mmap a relatively
     * large amount up front, and then carve off little pieces as we need them
     *
     * Moreover, some systems have too large a physical page size in their MMU
     * for the code to handle (we want offsets and lengths to fit into 16
     * bits), so we carve such big physical pages into smaller logical Pages
     * if needed. The largest logical page size is currently 16 KiB.
     */
    class PagePool {
        private:
            /// convenience typedef
            typedef BidirMMapPipeException Exception;

            enum {
                minsz = 7, ///< minimum chunk size (just below 1 << minsz bytes)
                maxsz = 20, ///< maximum chunk size (just below 1 << maxsz bytes)
                szincr = 1 ///< size class increment (sz = 1 << (minsz + k * szincr))
            };
            /// a chunk of memory in the pool
            typedef BidirMMapPipe_impl::PageChunk Chunk;
            /// list of chunks
            typedef std::list<Chunk*> ChunkList;

            friend class BidirMMapPipe_impl::PageChunk;
        public:
            /// convenience typedef
            typedef PageChunk::MMapVariety MMapVariety;
            /// constructor
            PagePool(unsigned nPagesPerGroup);
            /// destructor
            ~PagePool();
            /// pop a free element out of the pool
            Pages pop();

            /// return (logical) page size of the system
            static unsigned pagesize() { return PageChunk::pagesize(); }
            /// return variety of mmap supported on the system
            static MMapVariety mmapVariety()
            { return PageChunk::mmapVariety(); }

            /// return number of pages per group (ie. as returned by pop())
            unsigned nPagesPerGroup() const { return m_nPgPerGrp; }

            /// zap the pool (unmap all but Pages p)
            void zap(Pages& p);
            void clear_freelist();

        private:
            /// list of chunks used by the pool
            ChunkList m_chunks;
            /// list of chunks used by the pool which are not full
            ChunkList m_freelist;
            /// chunk size map (histogram of chunk sizes)
            unsigned m_szmap[(maxsz - minsz) / szincr];
            /// current chunk size
            int m_cursz;
            /// page group size
            unsigned m_nPgPerGrp;

            /// adjust _cursz to current largest block
            void updateCurSz(int sz, int incr);
            /// find size of next chunk to allocate (in a hopefully smart way)
            int nextChunkSz() const;
            /// release a chunk
            void putOnFreeList(Chunk* chunk);
            /// release a chunk
            void release(Chunk* chunk);
    };

    Pages::Pages(PageChunk* parent, Page* pages, unsigned npg) :
        m_pimpl(new impl)
    {
        assert(npg < 256);
        m_pimpl->m_parent = parent;
        m_pimpl->m_pages = pages;
        m_pimpl->m_refcnt = 1;
        m_pimpl->m_npages = npg;
        /// initialise pages
        for (unsigned i = 0; i < m_pimpl->m_npages; ++i) new(page(i)) Page();
    }

    unsigned PageChunk::s_physpgsz = PageChunk::getPageSize();
    unsigned PageChunk::s_pagesize = std::min(PageChunk::s_physpgsz, 16384u);
    PageChunk::MMapVariety PageChunk::s_mmapworks = PageChunk::Unknown;

    Pages::~Pages()
    {
        if (m_pimpl && !--(m_pimpl->m_refcnt)) {
            if (m_pimpl->m_parent) m_pimpl->m_parent->push(*this);
            delete m_pimpl;
        }
    }

    Pages::Pages(const Pages& other) :
        m_pimpl(other.m_pimpl)
    { ++(m_pimpl->m_refcnt); }

    Pages& Pages::operator=(const Pages& other)
    {
        if (&other == this) return *this;
        if (!--(m_pimpl->m_refcnt)) {
            if (m_pimpl->m_parent) m_pimpl->m_parent->push(*this);
            delete m_pimpl;
        }
        m_pimpl = other.m_pimpl;
        ++(m_pimpl->m_refcnt);
        return *this;
    }

    unsigned Pages::pagesize() { return PageChunk::pagesize(); }

    Page* Pages::page(unsigned pgno) const
    {
        assert(pgno < m_pimpl->m_npages);
        unsigned char* pptr =
            reinterpret_cast<unsigned char*>(m_pimpl->m_pages);
        pptr += pgno * pagesize();
        return reinterpret_cast<Page*>(pptr);
    }

    unsigned Pages::pageno(Page* p) const
    {
        const unsigned char* pptr =
            reinterpret_cast<const unsigned char*>(p);
        const unsigned char* bptr =
            reinterpret_cast<const unsigned char*>(m_pimpl->m_pages);
        assert(0 == ((pptr - bptr) % pagesize()));
        const unsigned nr = (pptr - bptr) / pagesize();
        assert(nr < m_pimpl->m_npages);
        return nr;
    }

    unsigned PageChunk::getPageSize()
    {
        // find out page size of system
        long pgsz = sysconf(_SC_PAGESIZE);
        if (-1 == pgsz) throw Exception("sysconf", errno);
        if (pgsz > 512 && pgsz > long(sizeof(Page)))
            return pgsz;

        // in case of failure or implausible value, use a safe default: 4k
        // page size, and do not try to mmap
        s_mmapworks = Copy;
        return 1 << 12;
    }

    PageChunk::PageChunk(PagePool* parent,
            unsigned length, unsigned nPgPerGroup) :
        m_begin(dommap(length)),
        m_end(reinterpret_cast<void*>(
                    reinterpret_cast<unsigned char*>(m_begin) + length)),
        m_parent(parent), m_nPgPerGrp(nPgPerGroup), m_nUsedGrp(0)
    {
        // ok, push groups of pages onto freelist here
        unsigned char* p = reinterpret_cast<unsigned char*>(m_begin);
        unsigned char* pend = reinterpret_cast<unsigned char*>(m_end);
        while (p < pend) {
            m_freelist.push_back(reinterpret_cast<void*>(p));
            p += nPgPerGroup * PagePool::pagesize();
        }
    }

    PageChunk::~PageChunk()
    {
        if (m_parent) assert(empty());
        if (m_begin) domunmap(m_begin, len());
    }

    bool PageChunk::contains(const Pages& p) const
    { return p.m_pimpl->m_parent == this; }

    Pages PageChunk::pop()
    {
        assert(!m_freelist.empty());
        void* p = m_freelist.front();
        m_freelist.pop_front();
        ++m_nUsedGrp;
        return Pages(this, reinterpret_cast<Page*>(p), m_nPgPerGrp);
    }

    void PageChunk::push(const Pages& p)
    {
        assert(contains(p));
        bool wasempty = m_freelist.empty();
        m_freelist.push_front(reinterpret_cast<void*>(p[0u]));
        --m_nUsedGrp;
        if (m_parent) {
            // notify parent if we need to be put on the free list again
            if (wasempty) m_parent->putOnFreeList(this);
            // notify parent if we're empty
            if (empty()) return m_parent->release(this);
        }
    }

    void* PageChunk::dommap(unsigned len)
    {
        assert(len && 0 == (len % s_physpgsz));
        // ok, the idea here is to try the different methods of mmapping, and
        // choose the first one that works. we have four flavours:
        // 1 - anonymous mmap (best)
        // 2 - mmap of /dev/zero (about as good as anonymous mmap, but a tiny
        //     bit more tedious to set up, since you need to open/close a
        //     device file)
        // 3 - mmap of a temporary file (very tedious to set up - need to
        //     create a temporary file, delete it, make the underlying storage
        //     large enough, then mmap the fd and close it)
        // 4 - if all those fail, we malloc the buffers, and copy the data
        //     through the OS (then we're no better than normal pipes)
        static bool msgprinted = false;
        if (Anonymous == s_mmapworks || Unknown == s_mmapworks) {
#if defined(MAP_ANONYMOUS)
#undef MYANONFLAG
#define MYANONFLAG MAP_ANONYMOUS
#elif defined(MAP_ANON)
#undef MYANONFLAG
#define MYANONFLAG MAP_ANON
#else
#undef MYANONFLAG
#endif
#ifdef MYANONFLAG
            void* retVal = ::mmap(0, len, PROT_READ | PROT_WRITE,
                    MYANONFLAG | MAP_SHARED, -1, 0);
            if (MAP_FAILED == retVal) {
                if (Anonymous == s_mmapworks) throw Exception("mmap", errno);
            } else {
                assert(Unknown == s_mmapworks || Anonymous == s_mmapworks);
                s_mmapworks = Anonymous;
                if (BidirMMapPipe::debugflag() && !msgprinted) {
                    std::cerr << "   INFO: In " << __func__ << " (" <<
                        __FILE__ << ", line " << __LINE__ <<
                        "): anonymous mmapping works, excellent!" <<
                        std::endl;
                    msgprinted = true;
                }
                return retVal;
            }
#endif
#undef MYANONFLAG
        }
        if (DevZero == s_mmapworks || Unknown == s_mmapworks) {
            // ok, no anonymous mappings supported directly, so try to map
            // /dev/zero which has much the same effect on many systems
            int fd = ::open("/dev/zero", O_RDWR);
            if (-1 == fd)
                throw Exception("open /dev/zero", errno);
            void* retVal = ::mmap(0, len,
                    PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            if (MAP_FAILED == retVal) {
                int errsv = errno;
                ::close(fd);
                if (DevZero == s_mmapworks) throw Exception("mmap", errsv);
            } else {
                assert(Unknown == s_mmapworks || DevZero == s_mmapworks);
                s_mmapworks = DevZero;
            }
            if (-1 == ::close(fd))
                throw Exception("close /dev/zero", errno);
            if (BidirMMapPipe::debugflag() && !msgprinted) {
                std::cerr << "   INFO: In " << __func__ << " (" << __FILE__ <<
                    ", line " << __LINE__ << "): mmapping /dev/zero works, "
                    "very good!" << std::endl;
                msgprinted = true;
            }
            return retVal;
        }
        if (FileBacked == s_mmapworks || Unknown == s_mmapworks) {
            char name[] = "/tmp/BidirMMapPipe-XXXXXX";
            int fd;
            // open temp file
            if (-1 == (fd = ::mkstemp(name))) throw Exception("mkstemp", errno);
            // remove it, but keep fd open
            if (-1 == ::unlink(name)) {
                int errsv = errno;
                ::close(fd);
                throw Exception("unlink", errsv);
            }
            // make it the right size: lseek
            if (-1 == ::lseek(fd, len - 1, SEEK_SET)) {
                int errsv = errno;
                ::close(fd);
                throw Exception("lseek", errsv);
            }
            // make it the right size: write a byte
            if (1 != ::write(fd, name, 1)) {
                int errsv = errno;
                ::close(fd);
                throw Exception("write", errsv);
            }
            // do mmap
            void* retVal = ::mmap(0, len,
                    PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            if (MAP_FAILED == retVal) {
                int errsv = errno;
                ::close(fd);
                if (FileBacked == s_mmapworks) throw Exception("mmap", errsv);
            } else {
                assert(Unknown == s_mmapworks || FileBacked == s_mmapworks);
                s_mmapworks = FileBacked;
            }
            if (-1 == ::close(fd)) {
                int errsv = errno;
                ::munmap(retVal, len);
                throw Exception("close", errsv);
            }
            if (BidirMMapPipe::debugflag() && !msgprinted) {
                std::cerr << "   INFO: In " << __func__ << " (" << __FILE__ <<
                    ", line " << __LINE__ << "): mmapping temporary files "
                    "works, good!" << std::endl;
                msgprinted = true;
            }
            return retVal;
        }
        if (Copy == s_mmapworks || Unknown == s_mmapworks) {
            // fallback solution: mmap does not work on this OS (or does not
            // work for what we want to use it), so use a normal buffer of
            // memory instead, and collect data in that buffer - this needs an
            // additional write/read to/from the pipe(s), but there you go...
            if (BidirMMapPipe::debugflag() && !msgprinted) {
                std::cerr << "WARNING: In " << __func__ << " (" << __FILE__ <<
                    ", line " << __LINE__ << "): anonymous mmapping of "
                    "shared buffers failed, falling back to read/write on "
                    " pipes!" << std::endl;
                msgprinted = true;
            }
            s_mmapworks = Copy;
            void* retVal = std::malloc(len);
            if (!retVal) throw Exception("malloc", errno);
            return retVal;
        }
        // should never get here
        assert(false);
        return 0;
    }

    void PageChunk::domunmap(void* addr, unsigned len)
    {
        assert(len && 0 == (len % s_physpgsz));
        if (addr) {
            assert(Unknown != s_mmapworks);
            if (Copy != s_mmapworks) {
                if (-1 == ::munmap(addr, len))
                    throw Exception("munmap", errno);
            } else {
                std::free(addr);
            }
        }
    }

    void PageChunk::zap(Pages& p)
    {
        // try to mprotect the other bits of the pool with no access...
        // we'd really like a version of mremap here that can unmap all the
        // other pages in the chunk, but that does not exist, so we protect
        // the other pages in this chunk such that they may neither be read,
        // written nor executed, only the pages we're interested in for
        // communications stay readable and writable
        //
        // if an OS does not support changing the protection of a part of an
        // mmapped area, the mprotect calls below should just fail and not
        // change any protection, so we're a little less safe against
        // corruption, but everything should still work
        if (Copy != s_mmapworks) {
            unsigned char* p0 = reinterpret_cast<unsigned char*>(m_begin);
            unsigned char* p1 = reinterpret_cast<unsigned char*>(p[0u]);
            unsigned char* p2 = p1 + p.npages() * s_physpgsz;
            unsigned char* p3 = reinterpret_cast<unsigned char*>(m_end);
            if (p1 != p0) ::mprotect(p0, p1 - p0, PROT_NONE);
            if (p2 != p3) ::mprotect(p2, p3 - p2, PROT_NONE);
        }
        m_parent = 0;
        m_freelist.clear();
        m_nUsedGrp = 1;
        p.m_pimpl->m_parent = 0;
        m_begin = m_end = 0;
        // commit suicide
        delete this;
    }

    PagePool::PagePool(unsigned nPgPerGroup) :
        m_cursz(minsz), m_nPgPerGrp(nPgPerGroup)
    {
        // if logical and physical page size differ, we may have to adjust
        // m_nPgPerGrp to make things fit
        if (PageChunk::pagesize() != PageChunk::physPgSz()) {
            const unsigned mult =
                PageChunk::physPgSz() / PageChunk::pagesize();
            const unsigned desired = nPgPerGroup * PageChunk::pagesize();
            // round up to to next physical page boundary
            const unsigned actual = mult *
                (desired / mult + bool(desired % mult));
            const unsigned newPgPerGrp = actual / PageChunk::pagesize();
            if (BidirMMapPipe::debugflag()) {
                std::cerr << "   INFO: In " << __func__ << " (" <<
                    __FILE__ << ", line " << __LINE__ <<
                    "): physical page size " << PageChunk::physPgSz() <<
                    ", subdividing into logical pages of size " <<
                    PageChunk::pagesize() << ", adjusting nPgPerGroup " <<
                    m_nPgPerGrp << " -> " << newPgPerGrp <<
                    std::endl;
            }
            assert(newPgPerGrp >= m_nPgPerGrp);
            m_nPgPerGrp = newPgPerGrp;
        }
        std::fill(m_szmap, m_szmap + ((maxsz - minsz) / szincr), 0);
    }

    PagePool::~PagePool()
    {
        m_freelist.clear();
        for (ChunkList::iterator it = m_chunks.begin(); m_chunks.end() != it; ++it)
            delete *it;
        m_chunks.clear();
    }

    void PagePool::zap(Pages& p)
    {
        // unmap all pages but those pointed to by p
        m_freelist.clear();
        for (ChunkList::iterator it = m_chunks.begin(); m_chunks.end() != it; ++it) {
            if ((*it)->contains(p)) {
                (*it)->zap(p);
            } else {
                delete *it;
            }
        }
        m_chunks.clear();
        std::fill(m_szmap, m_szmap + ((maxsz - minsz) / szincr), 0);
        m_cursz = minsz;
    }

    void PagePool::clear_freelist()
    {
//        m_freelist.clear();
        for (auto chunk_it = m_freelist.begin(); chunk_it != m_freelist.end();) {
            // clear only the chunks which are fully unused
            if ( (*chunk_it)->empty()) {
                chunk_it = m_freelist.erase(chunk_it);
            } else {
                ++chunk_it;
            }
        }

    }

    Pages PagePool::pop()
    {
        if (m_freelist.empty()) {
            // allocate and register new chunk and put it on the freelist
            const int sz = nextChunkSz();
            Chunk *c = new Chunk(this,
                    sz * m_nPgPerGrp * pagesize(), m_nPgPerGrp);
            m_chunks.push_front(c);
            m_freelist.push_back(c);
            updateCurSz(sz, +1);
        }
        // get free element from first chunk on _freelist
        Chunk* c = m_freelist.front();
        Pages p(c->pop());
        // full chunks are removed from _freelist
        if (c->full()) m_freelist.pop_front();
        return p;
    }

    void PagePool::release(PageChunk* chunk)
    {
        assert(chunk->empty());
        // find chunk on freelist and remove
        ChunkList::iterator it = std::find(
                m_freelist.begin(), m_freelist.end(), chunk);
        std::stringstream exception_message;
      exception_message << "PagePool::release(PageChunk*) on PID " << getpid();
        if (m_freelist.end() == it)
            throw Exception(exception_message.str().c_str(), EINVAL);
        m_freelist.erase(it);
        // find chunk in m_chunks and remove
        it = std::find(m_chunks.begin(), m_chunks.end(), chunk);
        if (m_chunks.end() == it)
            throw Exception("PagePool::release(PageChunk*)", EINVAL);
        m_chunks.erase(it);
        const unsigned sz = chunk->len() / (pagesize() * m_nPgPerGrp);
        delete chunk;
        updateCurSz(sz, -1);
    }

    void PagePool::putOnFreeList(PageChunk* chunk)
    {
        assert(!chunk->full());
        m_freelist.push_back(chunk);
    }

    void PagePool::updateCurSz(int sz, int incr)
    {
        m_szmap[(sz - minsz) / szincr] += incr;
        m_cursz = minsz;
        for (int i = (maxsz - minsz) / szincr; i--; ) {
            if (m_szmap[i]) {
                m_cursz += i * szincr;
                break;
            }
        }
    }

    int PagePool::nextChunkSz() const
    {
        // no chunks with space available, figure out chunk size
        int sz = m_cursz;
        if (m_chunks.empty()) {
            // if we start allocating chunks, we start from minsz
            sz = minsz;
        } else {
            if (minsz >= sz) {
                // minimal sized chunks are always grown
                sz = minsz + szincr;
            } else {
                if (1 != m_chunks.size()) {
                    // if we have more than one completely filled chunk, grow
                    sz += szincr;
                } else {
                    // just one chunk left, try shrinking chunk size
                    sz -= szincr;
                }
            }
        }
        // clamp size to allowed range
        if (sz > maxsz) sz = maxsz;
        if (sz < minsz) sz = minsz;
        return sz;
    }
}

// static BidirMMapPipe members
pthread_mutex_t BidirMMapPipe::s_openpipesmutex = PTHREAD_MUTEX_INITIALIZER;
std::list<BidirMMapPipe*> BidirMMapPipe::s_openpipes;
BidirMMapPipe_impl::PagePool* BidirMMapPipe::s_pagepool = 0;
unsigned BidirMMapPipe::s_pagepoolrefcnt = 0;
int BidirMMapPipe::s_debugflag = 0;

BidirMMapPipe_impl::PagePool& BidirMMapPipe::pagepool()
{
    if (!s_pagepool) {
        s_pagepool = new BidirMMapPipe_impl::PagePool(TotPages);
    }
    return *s_pagepool;
}

void BidirMMapPipe::teardownall(void)
{
    pthread_mutex_lock(&s_openpipesmutex);
    while (!s_openpipes.empty()) {
        BidirMMapPipe *p = s_openpipes.front();
        pthread_mutex_unlock(&s_openpipesmutex);
        if (p->m_childPid) kill(p->m_childPid, SIGTERM);
        p->doClose(true, true);
        pthread_mutex_lock(&s_openpipesmutex);
    }
    pthread_mutex_unlock(&s_openpipesmutex);
}

BidirMMapPipe::BidirMMapPipe(const BidirMMapPipe&) :
    m_pages(pagepool().pop())
{
    // free pages again
    { BidirMMapPipe_impl::Pages p; p.swap(m_pages); }
    if (!s_pagepoolrefcnt) {
        delete s_pagepool;
        s_pagepool = 0;
    }
}


// When creating with keepLocal, the master (which loses its bipe pointers to
// the child) must manually wait for its children after they are closed!
// Otherwise zombies are created. The static function wait_for_child can be
// used for this by master.
BidirMMapPipe::BidirMMapPipe(bool useExceptions, bool useSocketpair, bool keepLocal) :
    m_pages(pagepool().pop()), m_busylist(0), m_freelist(0), m_dirtylist(0),
    m_inpipe(-1), m_outpipe(-1), m_flags(failbit), m_childPid(0),
    m_parentPid(::getpid()), kept_local(keepLocal)
{
    ++s_pagepoolrefcnt;
    assert(0 < TotPages && 0 == (TotPages & 1) && TotPages <= 256);
    int fds[4] = { -1, -1, -1, -1 };
    int myerrno;
    static bool firstcall = true;
    if (useExceptions) m_flags |= exceptionsbit;

    try {
        if (firstcall) {
            firstcall = false;
            // register a cleanup handler to make sure all BidirMMapPipes are torn
            // down, and child processes are sent a SIGTERM
            if (0 != atexit(BidirMMapPipe::teardownall))
                throw Exception("atexit", errno);
        }

        // build free lists
        for (unsigned i = 1; i < TotPages; ++i)
            m_pages[i - 1]->setNext(m_pages[i]);
        m_pages[PagesPerEnd - 1]->setNext(0);
        if (!useSocketpair) {
            // create pipes
            if (0 != ::pipe(&fds[0])) throw Exception("pipe", errno);
            if (0 != ::pipe(&fds[2])) throw Exception("pipe", errno);
        } else {
            if (0 != ::socketpair(AF_UNIX, SOCK_STREAM, 0, &fds[0]))
                throw Exception("socketpair", errno);
        }
        // fork the child
        pthread_mutex_lock(&s_openpipesmutex);
        char c;
        m_childPid = ::fork();
//#ifndef NDEBUG
//        if (m_childPid == 0) {
//          std::cerr << "ignoring SIGTRAP in child process PID " << getpid() << " to allow debugger breakpoints in parent; signal(SIGTRAP, SIG_IGN) return code: " << signal(SIGTRAP, SIG_IGN) << std::endl;
//          sleep(10);
//        }
//#endif  // NDEBUG
        switch (m_childPid) {
            case -1: // error in fork()
                myerrno = errno;
                pthread_mutex_unlock(&s_openpipesmutex);
                m_childPid = 0;
                throw Exception("fork", myerrno);
            case 0: // child
                // put the ends in the right place
                if (-1 != fds[2]) {
                    // pair of pipes
                    if (-1 == ::close(fds[0]) || (-1 == ::close(fds[3]))) {
                        myerrno = errno;
                        pthread_mutex_unlock(&s_openpipesmutex);
                        throw Exception("close", myerrno);
                    }
                    fds[0] = fds[3] = -1;
                    m_outpipe = fds[1];
                    m_inpipe = fds[2];
                } else {
                    // socket pair
                    if (-1 == ::close(fds[0])) {
                        myerrno = errno;
                        pthread_mutex_unlock(&s_openpipesmutex);
                        throw Exception("close", myerrno);
                    }
                    fds[0] = -1;
                    m_inpipe = m_outpipe = fds[1];
                }
                // choose whether to retain other pipes on parent or on children
                // (the latter can be useful in special cases)
                if(keepLocal) {
                    for (std::list<BidirMMapPipe*>::iterator it = s_openpipes.begin();
                            s_openpipes.end() != it; ) {
                        BidirMMapPipe* p = *it;
                        it = s_openpipes.erase(it);
                        p->doClose(true, true);
                    }
                    // if new pages are made independently this should allow 
                    // new BidirMMapPipes to be constructed in the child ??
                    //~ s_pagepoolrefcnt = 0;
                    //~ delete s_pagepool;
                    //~ s_pagepool = 0;
                    pagepool().zap(m_pages);
                } else {
                    pagepool().clear_freelist();
                }
                s_openpipes.push_front(this);
                pthread_mutex_unlock(&s_openpipesmutex);
                //~ // ok, put our pages on freelist
                m_freelist = m_pages[PagesPerEnd];
                // handshare with other end (to make sure it's alive)...
                c = 'C'; // ...hild
                if (1 != xferraw(m_outpipe, &c, 1, ::write))
                    throw Exception("handshake: xferraw write", EPIPE);
                if (1 != xferraw(m_inpipe, &c, 1, ::read))
                    throw Exception("handshake: xferraw read", EPIPE);
                if ('P' != c) throw Exception("handshake", EPIPE);
                break;
            default: // parent
                // put the ends in the right place
                if (-1 != fds[2]) {
                    // pair of pipes
                    if (-1 == ::close(fds[1]) || -1 == ::close(fds[2])) {
                        myerrno = errno;
                        pthread_mutex_unlock(&s_openpipesmutex);
                        throw Exception("close", myerrno);
                    }
                    fds[1] = fds[2] = -1;
                    m_outpipe = fds[3];
                    m_inpipe = fds[0];
                } else {
                    // socketpair
                    if (-1 == ::close(fds[1])) {
                        myerrno = errno;
                        pthread_mutex_unlock(&s_openpipesmutex);
                        throw Exception("close", myerrno);
                    }
                    fds[1] = -1;
                    m_inpipe = m_outpipe = fds[0];
                }
                // put on list of open pipes (so we can kill child processes
                // if things go wrong)
                if(keepLocal) {
                // 
                } else {
                    for (std::list<BidirMMapPipe*>::iterator it = s_openpipes.begin();
                            s_openpipes.end() != it; ) {
                        BidirMMapPipe* p = *it;
                        it = s_openpipes.erase(it);
                        p->doClose(true, true);
                    }
                    pagepool().zap(m_pages);
                }
                s_openpipes.push_front(this);
                pthread_mutex_unlock(&s_openpipesmutex);
                // ok, put our pages on freelist
                m_freelist = m_pages[0u];

                // handshare with other end (to make sure it's alive)...
                c = 'P'; // ...arent
                if (1 != xferraw(m_outpipe, &c, 1, ::write))
                    throw Exception("handshake: xferraw write", EPIPE);
                if (1 != xferraw(m_inpipe, &c, 1, ::read))
                    throw Exception("handshake: xferraw read", EPIPE);
                if ('C' != c) throw Exception("handshake", EPIPE);
                break;
        }
        // mark file descriptors for close on exec (we do not want to leak the
        // connection to anything we happen to exec)
        int fdflags = 0;
        if (-1 == ::fcntl(m_outpipe, F_GETFD, &fdflags))
            throw Exception("fcntl", errno);
        fdflags |= FD_CLOEXEC;
        if (-1 == ::fcntl(m_outpipe, F_SETFD, fdflags))
            throw Exception("fcntl", errno);
        if (m_inpipe != m_outpipe) {
            if (-1 == ::fcntl(m_inpipe, F_GETFD, &fdflags))
                throw Exception("fcntl", errno);
            fdflags |= FD_CLOEXEC;
            if (-1 == ::fcntl(m_inpipe, F_SETFD, fdflags))
                throw Exception("fcntl", errno);
        }
        // ok, finally, clear the failbit
        m_flags &= ~failbit;
        // all done
    } catch (BidirMMapPipe::Exception&) {
        if (0 != m_childPid) kill(m_childPid, SIGTERM);
        for (int i = 0; i < 4; ++i)
            if (-1 != fds[i] && 0 != fds[i]) ::close(fds[i]);
        {
            // free resources associated with mmapped pages
            BidirMMapPipe_impl::Pages p; p.swap(m_pages);
        }
        if (!--s_pagepoolrefcnt) {
            delete s_pagepool;
            s_pagepool = 0;
        }
        throw;
    }
}

int BidirMMapPipe::close()
{
    assert(!(m_flags & failbit));
    return doClose(false);
}


// static function
int BidirMMapPipe::wait_for_child(pid_t child_pid, bool may_throw) {
  int status = 0;
  int tmp;
  do {
    tmp = waitpid(child_pid, &status, WNOHANG);
  } while (-1 == tmp && EINTR == errno);
  if (-1 == tmp && may_throw) throw Exception("waitpid", errno);

  if (WIFEXITED(status)) {
    printf("exited, status=%d\n", WEXITSTATUS(status));
  } else if (WIFSIGNALED(status)) {
    printf("killed by signal %d\n", WTERMSIG(status));
  } else if (WIFSTOPPED(status)) {
    printf("stopped by signal %d\n", WSTOPSIG(status));
  } else if (WIFCONTINUED(status)) {
    printf("continued\n");
  }

  return status;
}



int BidirMMapPipe::doClose(bool force, bool holdlock)
{
    if (m_flags & failbit) return 0;
    // flush data to be written
    if (!force && -1 != m_outpipe && -1 != m_inpipe) flush();
    // shut down the write direction (no more writes from our side)
    if (m_inpipe == m_outpipe) {
        if (-1 != m_outpipe && !force && -1 == ::shutdown(m_outpipe, SHUT_WR))
            throw Exception("shutdown", errno);
        m_outpipe = -1;
    } else {
        if (-1 != m_outpipe && -1 == ::close(m_outpipe))
            if (!force) throw Exception("close", errno);
        m_outpipe = -1;
    }
    // shut down the write direction (no more writes from our side)
    // drain anything the other end might still want to send
    if (!force && -1 != m_inpipe) {
        // **************** THIS IS EXTREMELY UGLY: ****************
        // POLLHUP is not set reliably on pipe/socket shutdown on all
        // platforms, unfortunately, so we poll for readability here until
        // the other end closes, too
        //
        // the read loop below ensures that the other end sees the POLLIN that
        // is set on shutdown instead, and goes ahead to close its end
        //
        // if we don't do this, and close straight away, the other end
        // will catch a SIGPIPE or similar, and we don't want that
        int err;
        struct pollfd fds;
        fds.fd = m_inpipe;
        fds.events = POLLIN;
        fds.revents = 0;
        do {
            while ((err = ::poll(&fds, 1, 1 << 20)) >= 0) {
                if (fds.revents & (POLLERR | POLLHUP | POLLNVAL)) break;
                if (fds.revents & POLLIN) {
                    char c;
                    if (1 > ::read(m_inpipe, &c, 1)) break;
                }
            }
        } while (0 > err && EINTR == errno);
        // ignore all other poll errors
    }
    // close read end
    if (-1 != m_inpipe && -1 == ::close(m_inpipe))
        if (!force) throw Exception("close", errno);
    m_inpipe = -1;
    // unmap memory
    try {
        { BidirMMapPipe_impl::Pages p; p.swap(m_pages); }
        assert(s_pagepoolrefcnt!=0);
        if (!--s_pagepoolrefcnt) {
            delete s_pagepool;
            s_pagepool = 0;
        }
    } catch (std::exception&) {
        if (!force) throw;
    }
    m_busylist = m_freelist = m_dirtylist = 0;
    // wait for child process
    int retval = 0;
    if (!kept_local && isParent() && ::getpid() == m_parentPid) { // double check whether actually parent
      retval = BidirMMapPipe::wait_for_child(m_childPid, !force);
      m_childPid = 0; // feeling that m_childPid can become zero when it should not be
    }

    // When created with keepLocal (so kept_local is true), the master (which
    // lost its bipe pointers to the child) must manually wait for its children
    // after they are closed! Otherwise zombies are created.

    // remove from list of open pipes
    if (!holdlock) pthread_mutex_lock(&s_openpipesmutex);
    std::list<BidirMMapPipe*>::iterator it = std::find(
            s_openpipes.begin(), s_openpipes.end(), this);
    if (s_openpipes.end() != it) s_openpipes.erase(it);
    if (!holdlock) pthread_mutex_unlock(&s_openpipesmutex);
    m_flags |= failbit;

    return retval;
}

BidirMMapPipe::~BidirMMapPipe()
{ 
  doClose(false); }

BidirMMapPipe::size_type BidirMMapPipe::xferraw(
        int fd, void* addr, size_type len,
        ssize_t (*xferfn)(int, void*, std::size_t))
{
    size_type xferred = 0;
    unsigned char* buf = reinterpret_cast<unsigned char*>(addr);
    while (len) {
        ssize_t tmp = xferfn(fd, buf, len);
        if (tmp > 0) {
            xferred += tmp;
            len -= tmp;
            buf += tmp;
            continue;
        } else if (0 == tmp) {
            // check for end-of-file on pipe
            break;
        } else if (-1 == tmp) {
            // ok some error occurred, so figure out if we want to retry of throw
            switch (errno) {
                default:
                    // if anything was transferred, return number of bytes
                    // transferred so far, we can start throwing on the next
                    // transfer...
                    if (xferred) return xferred;
                    // else throw
                    throw Exception("xferraw", errno);
                case EAGAIN: // fallthrough intended
#if defined(EWOULDBLOCK) && EWOULDBLOCK != EAGAIN
                case EWOULDBLOCK: // fallthrough intended
#endif
                    std::cerr << "  ERROR: In " << __func__ << " (" <<
                        __FILE__ << ", line " << __LINE__ <<
                        "): expect transfer to block!" << std::endl;
                case EINTR:
                    break;
            }
            continue;
        } else {
            throw Exception("xferraw: unexpected return value from read/write",
                    errno);
        }
    }
    return xferred;
}

void BidirMMapPipe::sendpages(Page* plist)
{
    if (plist) {
        unsigned char pg = m_pages[plist];
        if (1 == xferraw(m_outpipe, &pg, 1, ::write)) {
            if (BidirMMapPipe_impl::PageChunk::Copy ==
                    BidirMMapPipe_impl::PageChunk::mmapVariety()) {
                // ok, have to copy pages through pipe
                for (Page* p = plist; p; p = p->next()) {
                    if (sizeof(Page) + p->size() !=
                            xferraw(m_outpipe, p, sizeof(Page) + p->size(),
                                ::write)) {
                        throw Exception("sendpages: short write", EPIPE);
                    }
                }
            }
        } else {
            throw Exception("sendpages: short write", EPIPE);
        }
    } else { assert(plist); }
}

unsigned BidirMMapPipe::recvpages()
{
    unsigned char pg;
    unsigned retVal = 0;
    Page *plisthead = 0, *plisttail = 0;
    if (1 == xferraw(m_inpipe, &pg, 1, ::read)) {
        plisthead = plisttail = m_pages[pg];
        // ok, have number of pages
        if (BidirMMapPipe_impl::PageChunk::Copy ==
                BidirMMapPipe_impl::PageChunk::mmapVariety()) {
            // ok, need to copy pages through pipe
            for (; plisttail; ++retVal) {
                Page* p = plisttail;
                if (sizeof(Page) == xferraw(m_inpipe, p, sizeof(Page),
                            ::read)) {
                    plisttail = p->next();
                    if (!p->size()) continue;
                    // break in case of read error
                    if (p->size() != xferraw(m_inpipe, p->begin(), p->size(),
                                ::read)) break;
                }
            }
        } else {
            retVal = lenPageList(plisthead);
        }
    }
    // put list of pages we just received into correct lists (busy/free)
    if (plisthead) feedPageLists(plisthead);
    // ok, retVal contains the number of pages read, so put them on the
    // correct lists
    return retVal;
}

unsigned BidirMMapPipe::recvpages_nonblock()
{
    struct pollfd fds;
    fds.fd = m_inpipe;
    fds.events = POLLIN;
    fds.revents = 0;
    unsigned retVal = 0;
    do {
        int rc = ::poll(&fds, 1, 0);
        if (0 > rc) {
            if (EINTR == errno) continue;
            break;
        }
        if (1 == retVal && fds.revents & POLLIN &&
                !(fds.revents & (POLLNVAL | POLLERR))) {
            // ok, we can read without blocking, so the other end has
            // something for us
            return recvpages();
        } else {
            break;
        }
    } while (true);
    return retVal;
}

unsigned BidirMMapPipe::lenPageList(const Page* p)
{
    unsigned n = 0;
    for ( ; p; p = p->next()) ++n;
    return n;
}

void BidirMMapPipe::feedPageLists(Page* plist)
{
    assert(plist);
    // get end of busy list
    Page *blend = m_busylist;
    while (blend && blend->next()) blend = blend->next();
    // ok, might have to send free pages to other end, and (if we do have to
    // send something to the other end) while we're at it, send any dirty
    // pages which are completely full, too
    Page *sendlisthead = 0, *sendlisttail = 0;
    // loop over plist
    while (plist) {
        Page* p = plist;
        plist = p->next();
        p->setNext(0);
        if (p->size()) {
            // busy page...
            p->pos() = 0;
            // put at end of busy list
            if (blend) blend->setNext(p);
            else m_busylist = p;
            blend = p;
        } else {
            // free page...
            // Very simple algorithm: once we're done with a page, we send it back
            // where it came from. If it's from our end, we put it on the free list, if
            // it's from the other end, we send it back.
            if ((isParent() && m_pages[p] >= PagesPerEnd) ||
                    (isChild() && m_pages[p] < PagesPerEnd)) {
                // page "belongs" to other end
                if (!sendlisthead) sendlisthead = p;
                if (sendlisttail) sendlisttail->setNext(p);
                sendlisttail = p;
            } else {
                // add page to freelist
                p->setNext(m_freelist);
                m_freelist = p;
            }
        }
    }
    // check if we have to send stuff to the other end
    if (sendlisthead) {
        // go through our list of dirty pages, and see what we can
        // send along
        Page* dp;
        while ((dp = m_dirtylist) && dp->full()) {
            Page* p = dp;
            // move head of dirty list
            m_dirtylist = p->next();
            // queue for sending
            p->setNext(0);
            sendlisttail->setNext(p);
            sendlisttail = p;
        }
        // poll if the other end is still alive - this needs that we first
        // close the write pipe of the other end when the remote end of the
        // connection is shutting down in doClose; we'll see that because we
        // get a POLLHUP on our inpipe
        const int nfds = (m_outpipe == m_inpipe) ? 1 : 2;
        struct pollfd fds[2];
        fds[0].fd = m_outpipe;
        fds[0].events = fds[0].revents = 0;
        if (m_outpipe != m_inpipe) {
            fds[1].fd = m_inpipe;
            fds[1].events = fds[1].revents = 0;
        } else {
            fds[0].events |= POLLIN;
        }
        int retVal = 0;
        do {
            retVal = ::poll(fds, nfds, 0);
            if (0 > retVal && EINTR == errno)
                continue;
            break;
        } while (true);
        if (0 <= retVal) {
            bool ok = !(fds[0].revents & (POLLERR | POLLNVAL | POLLHUP));
            if (m_outpipe != m_inpipe) {
                ok = ok && !(fds[1].revents & (POLLERR | POLLNVAL | POLLHUP));
            } else {
                if (ok && fds[0].revents & POLLIN) {
                    unsigned ret = recvpages();
                    if (!ret) ok = false;
                }
            }

            if (ok) sendpages(sendlisthead);
            // (if the pipe is dead already, we don't care that we leak the
            // contents of the pages on the send list here, so that is why
            // there's no else clause here)
        } else {
            throw Exception("feedPageLists: poll", errno);
        }
    }
}

void BidirMMapPipe::markPageDirty(Page* p)
{
    assert(p);
    assert(p == m_freelist);
    // remove from freelist
    m_freelist = p->next();
    p->setNext(0);
    // append to dirty list
    Page* dl = m_dirtylist;
    while (dl && dl->next()) dl = dl->next();
    if (dl) dl->setNext(p);
    else m_dirtylist = p;
}

BidirMMapPipe::Page* BidirMMapPipe::busypage()
{
    // queue any pages available for reading we can without blocking
    recvpages_nonblock();
    Page* p;
    // if there are no busy pages, try to get them from the other end,
    // block if we have to...
    while (!(p = m_busylist)) if (!recvpages()) return 0;
    return p;
}

BidirMMapPipe::Page* BidirMMapPipe::dirtypage()
{
    // queue any pages available for reading we can without blocking
    recvpages_nonblock();
    Page* p = m_dirtylist;
    // go to end of dirty list
    if (p) while (p->next()) p = p->next();
    if (!p || p->full()) {
        // need to append free page, so get one
        while (!(p = m_freelist)) if (!recvpages()) return 0;
        markPageDirty(p);
    }
    return p;
}

void BidirMMapPipe::flush()
{ return doFlush(true); }

void BidirMMapPipe::doFlush(bool forcePartialPages)
{
    assert(!(m_flags & failbit));
    // build a list of pages to flush
    Page *flushlisthead = 0, *flushlisttail = 0;
    while (m_dirtylist) {
        Page* p = m_dirtylist;
        if (!forcePartialPages && !p->full()) break;
        // remove dirty page from dirty list
        m_dirtylist = p->next();
        p->setNext(0);
        // and send it to other end
        if (!flushlisthead) flushlisthead = p;
        if (flushlisttail) flushlisttail->setNext(p);
        flushlisttail = p;
    }
    if (flushlisthead) sendpages(flushlisthead);
}

void BidirMMapPipe::purge()
{
    assert(!(m_flags & failbit));
    // join busy and dirty lists
    {
        Page *l = m_busylist;
        while (l && l->next()) l = l->next();
        if (l) l->setNext(m_dirtylist);
        else m_busylist = m_dirtylist;
    }
    // empty busy and dirty pages
    for (Page* p = m_busylist; p; p = p->next()) p->size() = 0;
    // put them on the free list
    if (m_busylist) feedPageLists(m_busylist);
    m_busylist = m_dirtylist = 0;
}

BidirMMapPipe::size_type BidirMMapPipe::bytesReadableNonBlocking()
{
    // queue all pages waiting for consumption in the pipe before we give an
    // answer
    recvpages_nonblock();
    size_type retVal = 0;
    for (Page* p = m_busylist; p; p = p->next())
        retVal += p->size() - p->pos();
    return retVal;
}

BidirMMapPipe::size_type BidirMMapPipe::bytesWritableNonBlocking()
{
    // queue all pages waiting for consumption in the pipe before we give an
    // answer
    recvpages_nonblock();
    // check if we could write to the pipe without blocking (we need to know
    // because we might need to check if flushing of dirty pages would block)
    bool couldwrite = false;
    {
        struct pollfd fds;
        fds.fd = m_outpipe;
        fds.events = POLLOUT;
        fds.revents = 0;
        int retVal = 0;
        do {
            retVal = ::poll(&fds, 1, 0);
            if (0 > retVal) {
                if (EINTR == errno) continue;
                throw Exception("bytesWritableNonBlocking: poll", errno);
            }
            if (1 == retVal && fds.revents & POLLOUT &&
                    !(fds.revents & (POLLNVAL | POLLERR | POLLHUP)))
                couldwrite = true;
            break;
        } while (true);
    }
    // ok, start counting bytes
    size_type retVal = 0;
    unsigned npages = 0;
    // go through the dirty list
    for (Page* p = m_dirtylist; p; p = p->next()) {
        ++npages;
        // if page only partially filled
        if (!p->full())
            retVal += p->free();
        if (npages >= FlushThresh && !couldwrite) break;
    }
    // go through the free list
    for (Page* p = m_freelist; p && (!m_dirtylist ||
                npages < FlushThresh || couldwrite); p = p->next()) {
        ++npages;
        retVal += Page::capacity();
    }
    return retVal;
}

BidirMMapPipe::size_type BidirMMapPipe::read(void* addr, size_type sz)
{
    assert(!(m_flags & failbit));
    size_type nread = 0;
    unsigned char *ap = reinterpret_cast<unsigned char*>(addr);
    try {
        while (sz) {
            // find next page to read from
            Page* p = busypage();
            if (!p) {
                m_flags |= eofbit;
                return nread;
            }
            unsigned char* pp = p->begin() + p->pos();
            size_type csz = std::min(size_type(p->remaining()), sz);
            std::copy(pp, pp + csz, ap);
            nread += csz;
            ap += csz;
            sz -= csz;
            p->pos() += csz;
            assert(p->size() >= p->pos());
            if (p->size() == p->pos()) {
                // if no unread data remains, page is free
                m_busylist = p->next();
                p->setNext(0);
                p->size() = 0;
                feedPageLists(p);
            }
        }
    } catch (Exception&) {
        m_flags |= rderrbit;
        if (m_flags & exceptionsbit) throw;
    }
    return nread;
}

BidirMMapPipe::size_type BidirMMapPipe::write(const void* addr, size_type sz)
{
    assert(!(m_flags & failbit));
    size_type written = 0;
    const unsigned char *ap = reinterpret_cast<const unsigned char*>(addr);
    try {
        while (sz) {
            // find next page to write to
            Page* p = dirtypage();
            if (!p) {
                m_flags |= eofbit;
                return written;
            }
            unsigned char* pp = p->begin() + p->size();
            size_type csz = std::min(size_type(p->free()), sz);
            std::copy(ap, ap + csz, pp);
            written += csz;
            ap += csz;
            p->size() += csz;
            sz -= csz;
            assert(p->capacity() >= p->size());
            if (p->full()) {
                // if page is full, see if we're above the flush threshold of
                // 3/4 of our pages
                if (lenPageList(m_dirtylist) >= FlushThresh)
                    doFlush(false);
            }
        }
    } catch (Exception&) {
        m_flags |= wrerrbit;
        if (m_flags & exceptionsbit) throw;
    }
    return written;
}

int BidirMMapPipe::poll(BidirMMapPipe::PollVector& pipes, int timeout)
{
    // go through pipes, and change flags where we already know without really
    // polling - stuff where we don't need poll to wait for its timeout in the
    // OS...
    bool canskiptimeout = false;
    std::vector<unsigned> masks(pipes.size(), ~(Readable | Writable));
    std::vector<unsigned>::iterator mit = masks.begin();
    for (PollVector::iterator it = pipes.begin(); pipes.end() != it;
            ++it, ++mit) {
        PollEntry& pe = *it;
        pe.revents = None;
        // null pipe is invalid
        if (!pe.pipe) {
           pe.revents |= Invalid;
           canskiptimeout = true;
           continue;
        }
        // closed pipe is invalid
        if (pe.pipe->closed()) pe.revents |= Invalid;
        // check for error
        if (pe.pipe->bad()) pe.revents |= Error;
        // check for end of file
        if (pe.pipe->eof()) pe.revents |= EndOfFile;
        // check if readable
        if (pe.events & Readable) {
            *mit |= Readable;
            if (pe.pipe->m_busylist) pe.revents |= Readable;
        }
        // check if writable
        if (pe.events & Writable) {
            *mit |= Writable;
            if (pe.pipe->m_freelist) {
                pe.revents |= Writable;
            } else {
                Page *dl = pe.pipe->m_dirtylist;
                while (dl && dl->next()) dl = dl->next();
                if (dl && dl->pos() < Page::capacity())
                    pe.revents |= Writable;
            }
        }
        if (pe.revents) canskiptimeout = true;
    }
    // set up the data structures required for the poll syscall
    std::vector<pollfd> fds;
    fds.reserve(2 * pipes.size());
    std::map<int, PollEntry*> fds2pipes;
    for (PollVector::const_iterator it = pipes.begin();
            pipes.end() != it; ++it) {
        const PollEntry& pe = *it;
        struct pollfd tmp;
        fds2pipes.insert(std::make_pair((tmp.fd = pe.pipe->m_inpipe),
                    const_cast<PollEntry*>(&pe)));
        tmp.events = tmp.revents = 0;
        // we always poll for readability; this allows us to queue pages
        // early
        tmp.events |= POLLIN;
        if (pe.pipe->m_outpipe != tmp.fd) {
            // ok, it's a pair of pipes
            fds.push_back(tmp);
            fds2pipes.insert(std::make_pair(
                        unsigned(tmp.fd = pe.pipe->m_outpipe),
                        const_cast<PollEntry*>(&pe)));
            tmp.events = 0;

        }
        if (pe.events & Writable) tmp.events |= POLLOUT;
        fds.push_back(tmp);
    }
    // poll
    int retVal = 0;
    do {
        retVal = ::poll(&fds[0], fds.size(), canskiptimeout ? 0 : timeout);
        if (0 > retVal) {
            if (EINTR == errno) continue;
            throw Exception("poll", errno);
        }
        break;
    } while (true);
    // fds may have changed state, so update...
    for (std::vector<pollfd>::iterator it = fds.begin();
            fds.end() != it; ++it) {
        pollfd& fe = *it;
        //if (!fe.revents) continue;
        --retVal;
        PollEntry& pe = *fds2pipes[fe.fd];
oncemore:
        if (fe.revents & POLLNVAL && fe.fd == pe.pipe->m_inpipe)
            pe.revents |= ReadInvalid;
        if (fe.revents & POLLNVAL && fe.fd == pe.pipe->m_outpipe)
            pe.revents |= WriteInvalid;
        if (fe.revents & POLLERR && fe.fd == pe.pipe->m_inpipe)
            pe.revents |= ReadError;
        if (fe.revents & POLLERR && fe.fd == pe.pipe->m_outpipe)
            pe.revents |= WriteError;
        if (fe.revents & POLLHUP && fe.fd == pe.pipe->m_inpipe)
            pe.revents |= ReadEndOfFile;
        if (fe.revents & POLLHUP && fe.fd == pe.pipe->m_outpipe)
            pe.revents |= WriteEndOfFile;
        if ((fe.revents & POLLIN) && fe.fd == pe.pipe->m_inpipe &&
                !(fe.revents & (POLLNVAL | POLLERR))) {
            // ok, there is at least one page for us to receive from the
            // other end
            if (0 == pe.pipe->recvpages()) continue;
            // more pages there?
            do {
                int tmp = ::poll(&fe, 1, 0);
                if (tmp > 0) goto oncemore; // yippie! I don't even feel bad!
                if (0 > tmp) {
                    if (EINTR == errno) continue;
                    throw Exception("poll", errno);
                }
                break;
            } while (true);
        }
        if (pe.pipe->m_busylist) pe.revents |= Readable;
        if (fe.revents & POLLOUT && fe.fd == pe.pipe->m_outpipe) {
            if (pe.pipe->m_freelist) {
                pe.revents |= Writable;
            } else {
                Page *dl = pe.pipe->m_dirtylist;
                while (dl && dl->next()) dl = dl->next();
                if (dl && dl->pos() < Page::capacity())
                    pe.revents |= Writable;
            }
        }
    }
    // apply correct masks, and count pipes with pending events
    int npipes = 0;
    mit = masks.begin();
    for (PollVector::iterator it = pipes.begin();
            pipes.end() != it; ++it, ++mit)
        if ((it->revents &= *mit)) ++npipes;
    return npipes;
}

BidirMMapPipe& BidirMMapPipe::operator<<(const char* str)
{
    size_t sz = std::strlen(str);
    *this << sz;
    if (sz) write(str, sz);
    return *this;
}

BidirMMapPipe& BidirMMapPipe::operator>>(char* (&str))
{
    size_t sz = 0;
    *this >> sz;
    if (good() && !eof()) {
        str = reinterpret_cast<char*>(std::realloc(str, sz + 1));
        if (!str) throw Exception("realloc", errno);
        if (sz) read(str, sz);
        str[sz] = 0;
    }
    return *this;
}

BidirMMapPipe& BidirMMapPipe::operator<<(const std::string& str)
{
    size_t sz = str.size();
    *this << sz;
    write(str.data(), sz);
    return *this;
}

BidirMMapPipe& BidirMMapPipe::operator>>(std::string& str)
{
    str.clear();
    size_t sz = 0;
    *this >> sz;
    if (good() && !eof()) {
        str.reserve(sz);
        for (unsigned char c; sz--; str.push_back(c)) *this >> c;
    }
    return *this;
}

END_NAMESPACE_ROOFIT

#ifdef TEST_BIDIRMMAPPIPE
using namespace RooFit;

int simplechild(BidirMMapPipe& pipe)
{
    // child does an echo loop
    while (pipe.good() && !pipe.eof()) {
        // read a string
        std::string str;
        pipe >> str;
        if (!pipe) return -1;
        if (pipe.eof()) break;
        if (!str.empty()) {
            std::cout << "[CHILD] :  read: " << str << std::endl;
            str = "... early in the morning?";
        }
        pipe << str << BidirMMapPipe::flush;
        // did our parent tell us to shut down?
        if (str.empty()) break;
        if (!pipe) return -1;
        if (pipe.eof()) break;
        std::cout << "[CHILD] : wrote: " << str << std::endl;
    }
    pipe.close();
    return 0;
}

#include <sstream>
int randomchild(BidirMMapPipe& pipe)
{
    // child sends out something at random intervals
    ::srand48(::getpid());
    {
        // wait for parent's go ahead signal
        std::string s;
        pipe >> s;
    }
    // no shutdown sequence needed on this side - we're producing the data,
    // and the parent can just read until we're done (when it'll get EOF)
    for (int i = 0; i < 5; ++i) {
        // sleep a random time between 0 and .9 seconds
        ::usleep(int(1e6 * ::drand48()));
        std::ostringstream buf;
        buf << "child pid " << ::getpid() << " sends message " << i;
        std::string str = buf.str();
        std::cout << "[CHILD] : " << str << std::endl;
        pipe << str << BidirMMapPipe::flush;
        if (!pipe) return -1;
        if (pipe.eof()) break;
    }
    // tell parent we're shutting down
    pipe << "" << BidirMMapPipe::flush;
    // wait for parent to acknowledge
    std::string s;
    pipe >> s;
    pipe.close();
    return 0;
}

int benchchildrtt(BidirMMapPipe& pipe)
{
    // child does the equivalent of listening for pings and sending the
    // packet back
    char* str = 0;
    while (pipe && !pipe.eof()) {
        pipe >> str;
        if (!pipe) {
            std::free(str);
            pipe.close();
            return -1;
        }
        if (pipe.eof()) break;
        pipe << str << BidirMMapPipe::flush;
        // if we have just completed the shutdown handshake, we break here
        if (!std::strlen(str)) break;
    }
    std::free(str);
    pipe.close();
    return 0;
}

int benchchildsink(BidirMMapPipe& pipe)
{
    // child behaves like a sink
    char* str = 0;
    while (pipe && !pipe.eof()) {
        pipe >> str;
        if (!std::strlen(str)) break;
    }
    pipe << "" << BidirMMapPipe::flush;
    std::free(str);
    pipe.close();
    return 0;
}

int benchchildsource(BidirMMapPipe& pipe)
{
    // child behaves like a source
    char* str = 0;
    for (unsigned i = 0; i <= 24; ++i) {
        str = reinterpret_cast<char*>(std::realloc(str, (1 << i) + 1));
        std::memset(str, '4', 1 << i);
        str[1 << i] = 0;
        for (unsigned j = 0; j < 1 << 7; ++j) {
            pipe << str;
            if (!pipe || pipe.eof()) {
                std::free(str);
                pipe.close();
                return -1;
            }
        }
        // tell parent we're done with this block size
        pipe << "" << BidirMMapPipe::flush;
    }
    // tell parent to shut down
    pipe << "" << BidirMMapPipe::flush;
    std::free(str);
    pipe.close();
    return 0;
}

BidirMMapPipe* spawnChild(int (*childexec)(BidirMMapPipe&))
{
    // create a pipe with the given child at the remote end
    BidirMMapPipe *p = new BidirMMapPipe();
    if (p->isChild()) {
        int retVal = childexec(*p);
        delete p;
        std::exit(retVal);
    }
    return p;
}

#include <sys/time.h>
#include <iomanip>
int main()
{
    // simple echo loop test
    {
        std::cout << "[PARENT]: simple challenge-response test, "
            "one child:" << std::endl;
        BidirMMapPipe* pipe = spawnChild(simplechild);
        for (int i = 0; i < 5; ++i) {
            std::string str("What shall we do with a drunken sailor...");
            *pipe << str << BidirMMapPipe::flush;
            if (!*pipe) return -1;
            std::cout << "[PARENT]: wrote: " << str << std::endl;
            *pipe >> str;
            if (!*pipe) return -1;
            std::cout << "[PARENT]:  read: " << str << std::endl;
        }
        // send shutdown string
        *pipe << "" << BidirMMapPipe::flush;
        // wait for shutdown handshake
        std::string s;
        *pipe >> s;
        int retVal = pipe->close();
        std::cout << "[PARENT]: exit status of child: " << retVal <<
            std::endl;
        if (retVal) return retVal;
        delete pipe;
    }

    // simple poll test - children send 5 results in random intervals
    {
        unsigned nch = 20;
        std::cout << std::endl << "[PARENT]: polling test, " << nch <<
            " children:" << std::endl;
        typedef BidirMMapPipe::PollEntry PollEntry;
        // poll data structure
        BidirMMapPipe::PollVector pipes;
        pipes.reserve(nch);
        // spawn children
        for (unsigned i = 0; i < nch; ++i) {
            std::cout << "[PARENT]: spawning child " << i << std::endl;
            pipes.push_back(PollEntry(spawnChild(randomchild),
                        BidirMMapPipe::Readable));
        }
        // wake children up
        std::cout << "[PARENT]: waking up children" << std::endl;
        for (unsigned i = 0; i < nch; ++i)
            *pipes[i].pipe << "" << BidirMMapPipe::flush;
        std::cout << "[PARENT]: waiting for events on children's pipes" << std::endl;
        // while at least some children alive
        while (!pipes.empty()) {
            // poll, wait until status change (infinite timeout)
            int npipes = BidirMMapPipe::poll(pipes, -1);
            // scan for pipes with changed status
            for (std::vector<PollEntry>::iterator it = pipes.begin();
                    npipes && pipes.end() != it; ) {
                if (!it->revents) {
                    // unchanged, next one
                    ++it;
                    continue;
                }
                --npipes; // maybe we can stop early...
                // read from pipes which are readable
                if (it->revents & BidirMMapPipe::Readable) {
                    std::string s;
                    *(it->pipe) >> s;
                    if (!s.empty()) {
                        std::cout << "[PARENT]: Read from pipe " << it->pipe <<
                            ": " << s << std::endl;
                        ++it;
                        continue;
                    } else {
                        // child is shutting down...
                        *(it->pipe) << "" << BidirMMapPipe::flush;
                        goto childcloses;
                    }
                }
                // retire pipes with error or end-of-file condition
                if (it->revents & (BidirMMapPipe::Error |
                            BidirMMapPipe::EndOfFile |
                            BidirMMapPipe::Invalid)) {
                    std::cerr << "[DEBUG]: Event on pipe " << it->pipe <<
                        " revents" <<
                        ((it->revents & BidirMMapPipe::Readable) ? " Readable" : "") <<
                        ((it->revents & BidirMMapPipe::Writable) ? " Writable" : "") <<
                        ((it->revents & BidirMMapPipe::ReadError) ? " ReadError" : "") <<
                        ((it->revents & BidirMMapPipe::WriteError) ? " WriteError" : "") <<
                        ((it->revents & BidirMMapPipe::ReadEndOfFile) ? " ReadEndOfFile" : "") <<
                        ((it->revents & BidirMMapPipe::WriteEndOfFile) ? " WriteEndOfFile" : "") <<
                        ((it->revents & BidirMMapPipe::ReadInvalid) ? " ReadInvalid" : "") <<
                        ((it->revents & BidirMMapPipe::WriteInvalid) ? " WriteInvalid" : "") <<
                        std::endl;
childcloses:
                    int retVal = it->pipe->close();
                    std::cout << "[PARENT]: child exit status: " <<
                        retVal << ", number of children still alive: " <<
                        (pipes.size() - 1) << std::endl;
                    if (retVal) return retVal;
                    delete it->pipe;
                    it = pipes.erase(it);
                    continue;
                }
            }
        }
    }
    // little benchmark - round trip time
    {
        std::cout << std::endl << "[PARENT]: benchmark: round-trip times vs block size" << std::endl;
        for (unsigned i = 0; i <= 24; ++i) {
            char *s = new char[1 + (1 << i)];
            std::memset(s, 'A', 1 << i);
            s[1 << i] = 0;
            const unsigned n = 1 << 7;
            double avg = 0., min = 1e42, max = -1e42;
            BidirMMapPipe *pipe = spawnChild(benchchildrtt);
            for (unsigned j = n; j--; ) {
                struct timeval t1;
                ::gettimeofday(&t1, 0);
                *pipe << s << BidirMMapPipe::flush;
                if (!*pipe || pipe->eof()) break;
                *pipe >> s;
                if (!*pipe || pipe->eof()) break;
                struct timeval t2;
                ::gettimeofday(&t2, 0);
                t2.tv_sec -= t1.tv_sec;
                t2.tv_usec -= t1.tv_usec;
                double dt = 1e-6 * double(t2.tv_usec) + double(t2.tv_sec);
                if (dt < min) min = dt;
                if (dt > max) max = dt;
                avg += dt;
            }
            // send a shutdown string
            *pipe << "" << BidirMMapPipe::flush;
            // get child's shutdown ok
            *pipe >> s;
            avg /= double(n);
            avg *= 1e6; min *= 1e6; max *= 1e6;
            int retVal = pipe->close();
            if (retVal) {
                std::cout << "[PARENT]: child exited with code " << retVal << std::endl;
                delete[] s;
                return retVal;
            }
            delete pipe;
            // there is a factor 2 in the formula for the transfer rate below,
            // because we transfer data of twice the size of the block - once
            // to the child, and once for the return trip
            std::cout << "block size " << std::setw(9) << (1 << i) <<
                " avg " << std::setw(7) << avg << " us min " <<
                std::setw(7) << min << " us max " << std::setw(7) << max <<
                "us speed " << std::setw(9) <<
                2. * (double(1 << i) / double(1 << 20) / (1e-6 * avg)) <<
                " MB/s" << std::endl;
            delete[] s;
        }
        std::cout << "[PARENT]: all children had exit code 0" << std::endl;
    }
    // little benchmark - child as sink
    {
        std::cout << std::endl << "[PARENT]: benchmark: raw transfer rate with child as sink" << std::endl;
        for (unsigned i = 0; i <= 24; ++i) {
            char *s = new char[1 + (1 << i)];
            std::memset(s, 'A', 1 << i);
            s[1 << i] = 0;
            const unsigned n = 1 << 7;
            double avg = 0., min = 1e42, max = -1e42;
            BidirMMapPipe *pipe = spawnChild(benchchildsink);
            for (unsigned j = n; j--; ) {
                struct timeval t1;
                ::gettimeofday(&t1, 0);
                // streaming mode - we do not flush here
                *pipe << s;
                if (!*pipe || pipe->eof()) break;
                struct timeval t2;
                ::gettimeofday(&t2, 0);
                t2.tv_sec -= t1.tv_sec;
                t2.tv_usec -= t1.tv_usec;
                double dt = 1e-6 * double(t2.tv_usec) + double(t2.tv_sec);
                if (dt < min) min = dt;
                if (dt > max) max = dt;
                avg += dt;
            }
            // send a shutdown string
            *pipe << "" << BidirMMapPipe::flush;
            // get child's shutdown ok
            *pipe >> s;
            avg /= double(n);
            avg *= 1e6; min *= 1e6; max *= 1e6;
            int retVal = pipe->close();
            if (retVal) {
                std::cout << "[PARENT]: child exited with code " << retVal << std::endl;
                return retVal;
            }
            delete pipe;
            std::cout << "block size " << std::setw(9) << (1 << i) <<
                " avg " << std::setw(7) << avg << " us min " <<
                std::setw(7) << min << " us max " << std::setw(7) << max <<
                "us speed " << std::setw(9) <<
                (double(1 << i) / double(1 << 20) / (1e-6 * avg)) <<
                " MB/s" << std::endl;
            delete[] s;
        }
        std::cout << "[PARENT]: all children had exit code 0" << std::endl;
    }
    // little benchmark - child as source
    {
        std::cout << std::endl << "[PARENT]: benchmark: raw transfer rate with child as source" << std::endl;
        char *s = 0;
        double avg = 0., min = 1e42, max = -1e42;
        unsigned n = 0, bsz = 0;
        BidirMMapPipe *pipe = spawnChild(benchchildsource);
        while (*pipe && !pipe->eof()) {
            struct timeval t1;
            ::gettimeofday(&t1, 0);
            // streaming mode - we do not flush here
            *pipe >> s;
            if (!*pipe || pipe->eof()) break;
            struct timeval t2;
            ::gettimeofday(&t2, 0);
            t2.tv_sec -= t1.tv_sec;
            t2.tv_usec -= t1.tv_usec;
            double dt = 1e-6 * double(t2.tv_usec) + double(t2.tv_sec);
            if (std::strlen(s)) {
                ++n;
                if (dt < min) min = dt;
                if (dt > max) max = dt;
                avg += dt;
                bsz = std::strlen(s);
            } else {
                if (!n) break;
                // next block size
                avg /= double(n);
                avg *= 1e6; min *= 1e6; max *= 1e6;

                std::cout << "block size " << std::setw(9) << bsz <<
                    " avg " << std::setw(7) << avg << " us min " <<
                    std::setw(7) << min << " us max " << std::setw(7) <<
                    max << "us speed " << std::setw(9) <<
                    (double(bsz) / double(1 << 20) / (1e-6 * avg)) <<
                    " MB/s" << std::endl;
                n = 0;
                avg = 0.;
                min = 1e42;
                max = -1e42;
            }
        }
        int retVal = pipe->close();
            std::cout << "[PARENT]: child exited with code " << retVal << std::endl;
        if (retVal) return retVal;
        delete pipe;
        std::free(s);
    }
    return 0;
}
#endif // TEST_BIDIRMMAPPIPE
#endif // _WIN32

// vim: ft=cpp:sw=4:tw=78:et
