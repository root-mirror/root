/** @file BidirMMapPipe.h
 *
 * header file for BidirMMapPipe, a class which forks off a child process and
 * serves as communications channel between parent and child
 *
 * @author Manuel Schiller <manuel.schiller@nikhef.nl>
 * @author Patrick Bos <p.bos@esciencecenter.nl>
 * @author Inti Pelupessy <i.pelupessy@esciencecenter.nl>
 * @date 2013-2018
 */

#ifndef BIDIRMMAPPIPE_H
#define BIDIRMMAPPIPE_H

#include <list>
#include <vector>
#include <cstring>
#include <unistd.h>
#include <chrono>
#include "RooTaskSpec.h"

// forward declarations
namespace RooFit {
  namespace MultiProcessV1 {
    enum class M2Q;
    enum class Q2M;
    enum class W2Q;
    enum class Q2W;
  }
}

#define BEGIN_NAMESPACE_ROOFIT namespace RooFit {
#define END_NAMESPACE_ROOFIT }

BEGIN_NAMESPACE_ROOFIT

/// namespace for implementation details of BidirMMapPipe
namespace BidirMMapPipe_impl {
    // forward declarations
    class Page;
    class PagePool;
    class Pages;

  /** @brief exception to throw if low-level OS calls go wrong
  *
  * @author Manuel Schiller <manuel.schiller@nikhef.nl>
  * @date 2013-07-07
  */
  class BidirMMapPipeException : public std::exception  {
  private:
    enum {
      s_sz = 256 ///< length of buffer
    };
    char m_buf[s_sz]; ///< buffer containing the error message

    /// for the POSIX version of strerror_r
    static int dostrerror_r(int err, char* buf, std::size_t sz,
                            int (*f)(int, char*, std::size_t));
    /// for the GNU version of strerror_r
    static int dostrerror_r(int, char*, std::size_t,
                            char* (*f)(int, char*, std::size_t));
  public:
    /// constructor taking error code, hint on operation (msg)
    BidirMMapPipeException(const char* msg, int err);
    /// return a destcription of what went wrong
    virtual const char* what() const noexcept;
  };


  /** @brief class representing a chunk of pages
     *
     * @author Manuel Schiller <manuel.schiller@nikhef.nl>
     * @date 2013-07-24
     *
     * allocating pages from the OS happens in chunks in order to not exhaust
     * the maximum allowed number of memory mappings per process; this class
     * takes care of such a chunk
     *
     * a page chunk allows callers to obtain or release pages in groups of
     * continuous pages of fixed size
     */
    class PageChunk {
        public:
            /// type of mmap support found
            typedef enum {
                Unknown,        ///< don't know yet what'll work
                Copy,           ///< mmap doesn't work, have to copy back and forth
                FileBacked,     ///< mmapping a temp file works
                DevZero,        ///< mmapping /dev/zero works
                Anonymous       ///< anonymous mmap works
            } MMapVariety;

        private:
            static unsigned s_physpgsz; ///< system physical page size
            static unsigned s_pagesize; ///< logical page size (run-time determined)
            /// mmap variety that works on this system
            static MMapVariety s_mmapworks;

            /// convenience typedef
            typedef BidirMMapPipeException Exception;

            void* m_begin;      ///< pointer to start of mmapped area
            void* m_end;        ///< pointer one behind end of mmapped area
            // FIXME: cannot keep freelist inline - other end may need that
            //        data, and we'd end up overwriting the page header
            std::list<void*> m_freelist; ///< free pages list
            PagePool* m_parent; ///< parent page pool
            unsigned m_nPgPerGrp; ///< number of pages per group
            unsigned m_nUsedGrp; ///< number of used page groups

            /// determine page size at run time
            static unsigned getPageSize();

            /// mmap pages, len is length of mmapped area in bytes
            static void* dommap(unsigned len);
            /// munmap pages p, len is length of mmapped area in bytes
            static void domunmap(void* p, unsigned len);
            /// forbid copying
            PageChunk(const PageChunk&) {}
            /// forbid assignment
            PageChunk& operator=(const PageChunk&) { return *this; }
        public:
            /// return the logical page size
            static unsigned pagesize() { return s_pagesize; }
            /// return the physical page size of the system
            static unsigned physPgSz() { return s_physpgsz; }
            /// return mmap variety support found
            static MMapVariety mmapVariety() { return s_mmapworks; }

            /// constructor
            PageChunk(PagePool* parent, unsigned length, unsigned nPgPerGroup);

            /// destructor
            ~PageChunk();

            /// return if p is contained in this PageChunk
            bool contains(const Pages& p) const;

            /// pop a group of pages off the free list
            Pages pop();

            /// push a group of pages onto the free list
            void push(const Pages& p);

            /// return length of chunk
            unsigned len() const
            {
                return reinterpret_cast<unsigned char*>(m_end) -
                    reinterpret_cast<unsigned char*>(m_begin);
            }
            /// return number of pages per page group
            unsigned nPagesPerGroup() const { return m_nPgPerGrp; }

            /// return true if no used page groups in this chunk
            bool empty() const { return !m_nUsedGrp; }

            /// return true if no free page groups in this chunk
            bool full() const { return m_freelist.empty(); }

            /// free all pages except for those pointed to by p
            void zap(Pages& p);
    };

    /** @brief handle class for a number of Pages
     *
     * @author Manuel Schiller <manuel.schiller@nikhef.nl>
     * @date 2013-07-24
     *
     * the associated pages are continuous in memory
     */
    class Pages {
        private:
            /// implementation
            typedef struct {
                PageChunk *m_parent; ///< pointer to parent pool
                Page* m_pages; ///< pointer to first page
                unsigned m_refcnt; ///< reference counter
                unsigned char m_npages; ///< length in pages
            } impl;
        public:
            /// default constructor
            Pages() : m_pimpl(0) { }

            /// destructor
            ~Pages();

            /** @brief copy constructor
             *
             * copy Pages handle to new object - old object loses ownership,
             * and becomes a dangling handle
             */
            Pages(const Pages& other);

            /** @brief assignment operator
             *
             * assign Pages handle to new object - old object loses ownership,
             * and becomes a dangling handle
             */
            Pages& operator=(const Pages& other);

            /// return page size
            static unsigned pagesize();

            /// return number of pages accessible
            unsigned npages() const { return m_pimpl->m_npages; }

            /// return page number pageno
            Page* page(unsigned pgno) const;

            /// return page number pageno
            Page* operator[](unsigned pgno) const { return page(pgno); }

            /// perform page to page number mapping
            unsigned pageno(Page* p) const;

            /// perform page to page number mapping
            unsigned operator[](Page* p) const { return pageno(p); }

            /// swap with other's contents
            void swap(Pages& other)
            {
                impl* tmp = other.m_pimpl;
                other.m_pimpl = m_pimpl;
                m_pimpl = tmp;
            }

        private:
            /// page pool is our friend - it's allowed to construct Pages
            friend class BidirMMapPipe_impl::PageChunk;

            /// pointer to implementation
            impl* m_pimpl;

            /// constructor
            Pages(PageChunk* parent, Page* pages, unsigned npg);
    };
}

/** @brief BidirMMapPipe creates a bidirectional channel between the current
 * process and a child it forks.
 *
 * @author Manuel Schiller <manuel.schiller@nikhef.nl>
 * @date 2013-07-07
 *
 * This class creates a bidirectional channel between this process and a child
 * it creates with fork().
 *
 * The channel is comrised of a small shared pool of buffer memory mmapped into
 * both process spaces, and two pipes to synchronise the exchange of data. The
 * idea behind using the pipes at all is to have some primitive which we can
 * block on without having to worry about atomic operations or polling, leaving
 * these tasks to the OS. In case the anonymous mmap cannot be performed on the
 * OS the code is running on (for whatever reason), the code falls back to
 * mmapping /dev/zero, mmapping a temporary file, or (if those all fail), a
 * dynamically allocated buffer which is then transmitted through the pipe(s),
 * a slightly slower alternative (because the data is copied more often).
 *
 * The channel supports five major operations: read(), write(), flush(),
 * purge() and close(). Reading and writing may block until the required buffer
 * space is available. Writes may queue up data to be sent to the other end
 * until either enough pages are full, or the user calls flush which forces
 * any unsent buffers to be sent to the other end. flush forces any data that
 * is to be sent to be sent. purge discards any buffered data waiting to be
 * read and/or sent. Closing the channel on the child returns zero, closing it
 * on the parent returns the child's exit status.
 *
 * The class also provides operator<< and operator>> for C++-style I/O for
 * basic data types (bool, char, short, int, long, long long, float, double
 * and their unsigned counterparts). Data is transmitted binary (i.e. no
 * formatting to strings like std::cout does). There are also overloads to
 * support C-style zero terminated strings and std::string. In terms of
 * performance, the former is to be preferred.
 *
 * If the caller needs to multiplex input and output to/from several pipes, the
 * class provides the poll() method which allows to block until an event occurs
 * on any of the polled pipes.
 *
 * After the BidirMMapPipe is closed, no further operations may be performed on
 * that object, save for the destructor which may still be called.
 *
 * If the BidirMMapPipe has not properly been closed, the destructor will call
 * close. However, the exit code of the child is lost in that case.
 *
 * Closing the object causes the mmapped memory to be unmapped and the two
 * pipes to be closed. We also install an atexit handler in the process of
 * creating BidirMMapPipes. This ensures that when the current process
 * terminates, a SIGTERM signal is sent to the child processes created for all
 * unclosed pipes to avoid leaving zombie processes in the OS's process table.
 *
 * BidirMMapPipe creation, closing and destruction are thread safe. If the
 * BidirMMapPipe is used in more than one thread, the other operations have to
 * be protected with a mutex (or something similar), though.
 *
 * End of file (other end closed its pipe, or died) is indicated with the eof()
 * method, serious I/O errors set a flags (bad(), fail()), and also throw
 * exceptions. For normal read/write operations, they can be suppressed (i.e.
 * error reporting only using flags) with a constructor argument.
 *
 * Technicalities:
 * - there is a pool of mmapped pages, half the pages are allocated to the
 *   parent process, half to the child
 * - when one side has accumulated enough data (or a flush forces dirty pages
 *   out to the other end), it sends these pages to the other end by writing a
 *   byte containing the page number into the pipe
 * - the other end (which has the pages mmapped, too) reads the page number(s)
 *   and puts the corresponding pages on its busy list
 * - as the other ends reads, it frees busy pages, and eventually tries to put
 *   them on the its list; if a page belongs to the other end of the
 *   connection, it is sent back
 * - lists of pages are sent across the pipe, not individual pages, in order
 *   to minimise the number of read/write operations needed
 * - when mmap works properly, only one bytes containing the page number of
 *   the page list head is sent back and forth; the contents of that page
 *   allow to access the rest of the page list sent, and page headers on the
 *   list tell the receiving end if the page is free or has to be added to the
 *   busy list
 * - when mmap does not work, we transfer one byte to indicate the head of the
 *   page list sent, and for each page on the list of sent pages, the page
 *   header and the page payload is sent (if the page is free, we only
 *   transmit the page header, and we never transmit more payload than
 *   the page actually contains)
 * - in the child, all open BidirMMapPipes but the current one are closed. this
 *   is done for two reasons: first, to conserve file descriptors and address
 *   space. second, if more than one process is meant to use such a
 *   BidirMMapPipe, synchronisation issues arise which can lead to bugs that
 *   are hard to find and understand. it's much better to come up with a design
 *   which does not need pipes to be shared among more than two processes.
 *
 * Here is a trivial example of a parent and a child talking to each other over
 * a BidirMMapPipe:
 * @code
 * #include <string>
 * #include <iostream>
 * #include <cstdlib>
 *
 * #include <BidirMMapPipe.h>
 *
 * int simplechild(BidirMMapPipe& pipe)
 * {
 *     // child does an echo loop
 *     while (pipe.good() && !pipe.eof()) {
 *         // read a string
 *         std::string str;
 *         pipe >> str;
 *         if (!pipe) return -1;
 *         if (pipe.eof()) break;
 *         // check if parent wants us to shut down
 *         if (!str.empty()) {
 *             std::cout << "[CHILD] :  read: " << str << std::endl;
 *             str = "... early in the morning?";
 *         }
 *         pipe << str << BidirMMapPipe::flush;
 *         if (str.empty()) break;
 *         if (!pipe) return -1;
 *         std::cout << "[CHILD] : wrote: " << str << std::endl;
 *     }
 *     // send shutdown request acknowledged
 *     pipe << "" << BidirMMapPipe::flush;
 *
 *     pipe.close();
 *     return 0;
 * }
 *
 * BidirMMapPipe* spawnChild(int (*childexec)(BidirMMapPipe&))
 * {
 *     BidirMMapPipe *p = new BidirMMapPipe();
 *     if (p->isChild()) {
 *         int retVal = childexec(*p);
 *         delete p;
 *         std::exit(retVal);
 *     }
 *     return p;
 * }
 *
 * int main()
 * {
 *     std::cout << "[PARENT]: simple challenge-response test, one child:" <<
 *             std::endl;
 *     BidirMMapPipe* pipe = spawnChild(simplechild);
 *     for (int i = 0; i < 5; ++i) {
 *         std::string str("What shall we do with a drunken sailor...");
 *         *pipe << str << BidirMMapPipe::flush;
 *         if (!*pipe) return -1;
 *         std::cout << "[PARENT]: wrote: " << str << std::endl;
 *         *pipe >> str;
 *         if (!*pipe) return -1;
 *         std::cout << "[PARENT]:  read: " << str << std::endl;
 *     }
 *     // ask child to shut down
 *     pipe << "" << BidirMMapPipe::flush;
 *     // wait for it to see the shutdown request
 *     std::string s;
 *     pipe >> s;
 *     std::cout << "[PARENT]: exit status of child: " << pipe->close() <<
 *             std::endl;
 *     delete pipe;
 *     return 0;
 * }
 * @endcode
 *
 * When designing your own protocols to use over the pipe, there are a few
 * things to bear in mind:
 * - Do as http does: When building a request, send all the options and
 *   properties of that request with the request itself in a single go (one
 *   flush). Then, the server has everything it needs, and hopefully, it'll
 *   shut up for a while and to let the client do something useful in the
 *   meantime... The same goes when the server replies to the request: include
 *   everything there is to know about the result of the request in the reply.
 * - The expensive operation should be the request that is made, all other
 *   operations should somehow be formulated as options or properties to that
 *   request.
 * - Include a shutdown handshake in whatever protocol you send over the
 *   pipe. That way, you can shut things down in a controlled way. Otherwise,
 *   and depending on your OS's scheduling quirks, you may catch a SIGPIPE if
 *   one end closes its pipe while the other is still trying to read.
 */
class BidirMMapPipe {
#ifndef _WIN32
    public:
        /// type used to represent sizes
        typedef std::size_t size_type;
        /// convenience typedef for BidirMMapPipeException
        typedef BidirMMapPipe_impl::BidirMMapPipeException Exception;
        /// flag bits for partial C++ iostream compatibility
        enum {
            eofbit = 1, ///< end of file reached
            failbit = 2, ///< logical failure (e.g. pipe closed)
            rderrbit = 4, ///< read error
            wrerrbit = 8, ///< write error
            badbit = rderrbit | wrerrbit, ///< general I/O error
            exceptionsbit = 16 ///< error reporting with exceptions
        };

        /** @brief constructor (forks!)
         *
         * Creates a bidirectional communications channel between this process
         * and a child the constructor forks. On return from the constructor,
         * isParent() and isChild() can be used to tell the parent end from the
         * child end of the pipe. In the child, all other open BidirMMapPipes
         * are closed.
         *
         * @param useExceptions read()/write() error reporting also done using
         *                      exceptions
         * @param useSocketpair use a socketpair instead of a pair or pipes
         * @param keepLocal     deetermines whether existing pipes are retained 
         *                      at the parent (default) or transferred to children
         *
         * Normally, exceptions are thrown for all serious I/O errors (apart
         * from end of file). Setting useExceptions to false will force the
         * read() and write() methods to only report serious I/O errors using
         * flags.
         *
         * When useSocketpair is true, use a pair of Unix domain sockets
         * created using socketpair instead a pair of pipes. The advantage is
         * that only one pair of file descriptors is needed instead of two
         * pairs which are needed for the pipe pair. Performance should very
         * similar on most platforms, especially if mmap works, since only
         * very little data is sent through the pipe(s)/socketpair.
         */
        BidirMMapPipe(bool useExceptions = true, bool useSocketpair = false, bool keepLocal = true);

        /** @brief destructor
         *
         * closes this end of pipe
         */
        ~BidirMMapPipe();

        /** @brief return the current setting of the debug flag
         *
         * @returns an integer with the debug Setting
         */
        static int debugflag() { return s_debugflag; }

        /** @brief set the debug flags
         *
         * @param flag  debug flags (if zero, no messages are printed)
         */
        static void setDebugflag(int flag) { s_debugflag = flag; }

        /** @brief read from pipe
         *
         * @param addr  where to put read data
         * @param sz    size of data to read (in bytes)
         * @returns     size of data read, or 0 in case of end-of-file
         *
         * read may block until data from other end is available. It will
         * return 0 if the other end closed the pipe.
         */
        size_type read(void* addr, size_type sz);

        /** @brief wirte to pipe
         *
         * @param addr  where to get data to write from
         * @param sz    size of data to write (in bytes)
         * @returns     size of data written, or 0 in case of end-of-file
         *
         * write may block until data can be written to other end (depends a
         * bit on available buffer space). It will return 0 if the other end
         * closed the pipe. The data is queued to be written on the next
         * convenient occasion, or it can be forced out with flush().
         */
        size_type write(const void* addr, size_type sz);

        /** @brief flush buffers with unwritten data
         *
         * This forces unwritten data to be written to the other end. The call
         * will block until this has been done (or the attempt failed with an
         * error).
         */
        void flush();

        /** @brief purge buffered data waiting to be read and/or written
         *
         * Discards all internal buffers.
         */
        void purge();

        /** @brief number of bytes that can be read without blocking
         *
         * @returns number of bytes that can be read without blocking
         */
        size_type bytesReadableNonBlocking();

        /** @brief number of bytes that can be written without blocking
         *
         * @returns number of bytes that can be written without blocking
         */
        size_type bytesWritableNonBlocking();

        /** @brief flush buffers, close pipe
         *
         * Flush buffers, discard unread data, closes the pipe. If the pipe is
         * in the parent process, it waits for the child.
         *
         * @returns exit code of child process in parent, zero in child
         */
        int close();

        /** @brief return PID of the process on the other end of the pipe
         *
         * @returns PID of the process running on the remote end
         */
        pid_t pidOtherEnd() const
        { return isChild() ? m_parentPid : m_childPid; }

        /// condition flags for poll
        enum PollFlags {
            None = 0,           ///< nothing special on this pipe
            Readable = 1,       ///< pipe has data for reading
            Writable = 2,       ///< pipe can be written to
            ReadError = 4,      ///< pipe error read end
            WriteError = 8,     ///< pipe error Write end
            Error = ReadError | WriteError, ///< pipe error
            ReadEndOfFile = 32, ///< read pipe in end-of-file state
            WriteEndOfFile = 64,///< write pipe in end-of-file state
            EndOfFile = ReadEndOfFile | WriteEndOfFile, ///< end of file
            ReadInvalid = 64,   ///< read end of pipe invalid
            WriteInvalid = 128, ///< write end of pipe invalid
            Invalid = ReadInvalid | WriteInvalid ///< invalid pipe
        };

        /// for poll() interface
        class PollEntry {
            public:
                BidirMMapPipe* pipe;    ///< pipe of interest
                unsigned events;        ///< events of interest (or'ed bitmask)
                unsigned revents;       ///< events that happened (or'ed bitmask)
                /// poll a pipe for all events
                PollEntry(BidirMMapPipe* _pipe) :
                    pipe(_pipe), events(None), revents(None) { }
                /// poll a pipe for specified events
                PollEntry(BidirMMapPipe* _pipe, int _events) :
                    pipe(_pipe), events(_events), revents(None) { }
        };
        /// convenience typedef for poll() interface
        typedef std::vector<PollEntry> PollVector;

        /** @brief poll a set of pipes for events (ready to read from, ready to
         * write to, error)
         *
         * @param pipes         set of pipes to check
         * @param timeout       timeout in milliseconds
         * @returns             positive number: number of pipes which have
         *                      status changes, 0: timeout, or no pipes with
         *                      status changed, -1 on error
         *
         * Timeout can be zero (check for specified events, and return), finite
         * (wait at most timeout milliseconds before returning), or -1
         * (infinite). The poll method returns when the timeout has elapsed,
         * or if an event occurs on one of the pipes being polled, whichever
         * happens earlier.
         *
         * Pipes is a vector of one or more PollEntries, which each list a pipe
         * and events to poll for. If events is left empty (zero), all
         * conditions are polled for, otherwise only the indicated ones. On
         * return, the revents fields contain the events that occurred for each
         * pipe; error Error, EndOfFile or Invalid events are always set,
         * regardless of wether they were in the set of requested events.
         *
         * poll may block slightly longer than specified by timeout due to OS
         * timer granularity and OS scheduling. Due to its implementation, the
         * poll call can also return early if the remote end of the page sends
         * a free page while polling (which is put on that pipe's freelist),
         * while that pipe is polled for e.g Reading. The status of the pipe is
         * indicated correctly in revents, and the caller can simply poll
         * again. (The reason this is done this way is because it helps to
         * replenish the pool of free pages and queue busy pages without
         * blocking.)
         *
         * Here's a piece of example code waiting on two pipes; if they become
         * readable they are read:
         * @code
         * #include <unistd.h>
         * #include <cstdlib>
         * #include <string>
         * #include <sstream>
         * #include <iostream>
         *
         * #include <BidirMMapPipe.h>
         *
         * // what to execute in the child
         * int randomchild(BidirMMapPipe& pipe)
         * {
         *     ::srand48(::getpid());
         *     for (int i = 0; i < 5; ++i) {
         *         // sleep a random time between 0 and .9 seconds
         *        ::usleep(int(1e6 * ::drand48()));
         *        std::ostringstream buf;
         *        buf << "child pid " << ::getpid() << " sends message " << i;
         *        std::cout << "[CHILD] : " << buf.str() << std::endl;
         *        pipe << buf.str() << BidirMMapPipe::flush;
         *        if (!pipe) return -1;
         *        if (pipe.eof()) break;
         *     }
         *     // tell parent we're done
         *     pipe << "" << BidirMMapPipe::flush;
         *     // wait for parent to acknowledge
         *     std::string s;
         *     pipe >> s;
         *     pipe.close();
         *     return 0;
         * }
         *
         * // function to spawn a child
         * BidirMMapPipe* spawnChild(int (*childexec)(BidirMMapPipe&))
         * {
         *     BidirMMapPipe *p = new BidirMMapPipe();
         *     if (p->isChild()) {
         *         int retVal = childexec(*p);
         *         delete p;
         *         std::exit(retVal);
         *     }
         *     return p;
         * }
         *
         * int main()
         * {
         *     typedef BidirMMapPipe::PollEntry PollEntry;
         *     // poll data structure
         *     BidirMMapPipe::PollVector pipes;
         *     pipes.reserve(3);
         *     // spawn children
         *     for (int i = 0; i < 3; ++i) {
         *         pipes.push_back(PollEntry(spawnChild(randomchild),
         *              BidirMMapPipe::Readable));
         *     }
         *     // while at least some children alive
         *     while (!pipes.empty()) {
         *         // poll, wait until status change (infinite timeout)
         *         int npipes = BidirMMapPipe::poll(pipes, -1);
         *         // scan for pipes with changed status
         *         for (std::vector<PollEntry>::iterator it = pipes.begin();
         *                 npipes && pipes.end() != it; ) {
         *             if (!it->revents) {
         *                 // unchanged, next one
         *                 ++it;
         *                 continue;
         *             }
         *             --npipes; // maybe we can stop early...
         *             // read from pipes which are readable
         *             if (it->revents & BidirMMapPipe::Readable) {
         *                 std::string s;
         *                 *(it->pipe) >> s;
         *                 if (!s.empty()) {
         *                 std::cout << "[PARENT]: Read from pipe " <<
         *                     it->pipe << ": " << s << std::endl;
         *                 ++it;
         *                 continue;
         *                 } else {
         *                     // child is shutting down...
         *                     *(it->pipe) << "" << BidirMMapPipe::flush;
         *                     goto childcloses;
         *                 }
         *             }
         *             // retire pipes with error or end-of-file condition
         *             if (it->revents & (BidirMMapPipe::Error |
         *                     BidirMMapPipe::EndOfFile |
         *                     BidirMMapPipe::Invalid)) {
         *                 std::cout << "[PARENT]: Error on pipe " <<
         *                     it->pipe << " revents " << it->revents <<
         *                     std::endl;
         * childcloses:
         *                 std::cout << "[PARENT]:\tchild exit status: " <<
         *                     it->pipe->close() << std::endl;
         *                 if (retVal) return retVal;
         *                 delete it->pipe;
         *                 it = pipes.erase(it);
         *                 continue;
         *             }
         *         }
         *     }
         *     return 0;
         * }
         * @endcode
         */
        static int poll(PollVector& pipes, int timeout);

        /** @brief return if this end of the pipe is the parent end
         *
         * @returns true if parent end of pipe
         */
        bool isParent() const { return m_childPid; }

        /** @brief return if this end of the pipe is the child end
         *
         * @returns true if child end of pipe
         */
        bool isChild() const { return !m_childPid; }

        /** @brief if BidirMMapPipe uses a socketpair for communications
         *
         * @returns true if BidirMMapPipe uses a socketpair for communications
         */
        bool usesSocketpair() const { return m_inpipe == m_outpipe; }

        /** @brief if BidirMMapPipe uses a pipe pair for communications
         *
         * @returns true if BidirMMapPipe uses a pipe pair for communications
         */
        bool usesPipepair() const { return m_inpipe != m_outpipe; }

        /** @brief return flags (end of file, BidirMMapPipe closed, ...)
         *
         * @returns flags (end of file, BidirMMapPipe closed, ...)
         */
        int rdstate() const { return m_flags; }

        /** @brief true if end-of-file
         *
         * @returns true if end-of-file
         */
        bool eof() const { return m_flags & eofbit; }

        /** @brief logical failure (e.g. I/O on closed BidirMMapPipe)
         *
         * @returns true in case of grave logical error (I/O on closed pipe,...)
         */
        bool fail() const { return m_flags & failbit; }

        /** @brief true on I/O error
         *
         * @returns true on I/O error
         */
        bool bad() const { return m_flags & badbit; }

        /** @brief status of stream is good
         *
         * @returns true if pipe is good (no errors, eof, ...)
         */
        bool good() const { return !(m_flags & (eofbit | failbit | badbit)); }

        /** @brief true if closed
         *
         * @returns true if stream is closed
         */
        bool closed() const { return m_flags & failbit; }

        /** @brief return true if not serious error (fail/bad)
         *
         * @returns true if stream is does not have serious error (fail/bad)
         *
         * (if EOF, this is still true)
         */
        operator bool() const { return !fail() && !bad(); }

        /** @brief return true if serious error (fail/bad)
         *
         * @returns true if stream has a serious error (fail/bad)
         */
        bool operator!() const { return fail() || bad(); }

#ifdef STREAMOP
#undef STREAMOP
#endif
#define STREAMOP(TYPE) \
        BidirMMapPipe& operator<<(const TYPE& val) \
        { write(&val, sizeof(TYPE)); return *this; } \
        BidirMMapPipe& operator>>(TYPE& val) \
        { read(&val, sizeof(TYPE)); return *this; }
        STREAMOP(bool); ///< C++ style stream operators for bool
        STREAMOP(char); ///< C++ style stream operators for char
        STREAMOP(short); ///< C++ style stream operators for short
        STREAMOP(int); ///< C++ style stream operators for int
        STREAMOP(long); ///< C++ style stream operators for long
        STREAMOP(long long); ///< C++ style stream operators for long long
        STREAMOP(unsigned char); ///< C++ style stream operators for unsigned char
        STREAMOP(unsigned short); ///< C++ style stream operators for unsigned short
        STREAMOP(unsigned int); ///< C++ style stream operators for unsigned int
        STREAMOP(unsigned long); ///< C++ style stream operators for unsigned long
        STREAMOP(unsigned long long); ///< C++ style stream operators for unsigned long long
        STREAMOP(float); ///< C++ style stream operators for float
        STREAMOP(double); ///< C++ style stream operators for double
#undef STREAMOP

        /** @brief write a C-style string
         *
         * @param str C-style string
         * @returns pipe written to
         */
        BidirMMapPipe& operator<<(const char* str);

        /** @brief read a C-style string
         *
         * @param str pointer to string (space allocated with malloc!)
         * @returns pipe read from
         *
         * since this is for C-style strings, we use malloc/realloc/free for
         * strings. passing in a NULL pointer is valid here, and the routine
         * will use realloc to allocate a chunk of memory of the right size.
         */
        BidirMMapPipe& operator>>(char* (&str));

        /** @brief write a std::string object
         *
         * @param str string to write
         * @returns pipe written to
         */
        BidirMMapPipe& operator<<(const std::string& str);

        /** @brief read a std::string object
         *
         * @param str string to be read
         * @returns pipe read from
         */
        BidirMMapPipe& operator>>(std::string& str);

        /** @brief write raw pointer to T to other side
         *
         * NOTE: This will not write the pointee! Only the value of the
         * pointer is transferred.
         *
         * @param tptr pointer to be written
         * @returns pipe written to
         */
        template<class T> BidirMMapPipe& operator<<(const T* tptr)
	  {   write(&tptr, sizeof(tptr)); return *this; }

        /** @brief read raw pointer to T from other side
         *
         * NOTE: This will not read the pointee! Only the value of the
         * pointer is transferred.
         *
         * @param tptr pointer to be read
         * @returns pipe read from
         */
        template<class T> BidirMMapPipe& operator>>(T* &tptr)
        { read(&tptr, sizeof(tptr)); return *this; }

        /** @brief I/O manipulator support
         *
         * @param manip manipulator
         * @returns pipe with manipulator applied
         *
         * example:
         * @code
         * pipe << BidirMMapPipe::flush;
         * @endcode
         */
        BidirMMapPipe& operator<<(BidirMMapPipe& (*manip)(BidirMMapPipe&))
        { return manip(*this); }

        /** @brief I/O manipulator support
         *
         * @param manip manipulator
         * @returns pipe with manipulator applied
         *
         * example:
         * @code
         * pipe >> BidirMMapPipe::purge;
         * @endcode
         */
        BidirMMapPipe& operator>>(BidirMMapPipe& (*manip)(BidirMMapPipe&))
        { return manip(*this); }

        /// for usage a la "pipe << flush;"
        static BidirMMapPipe& flush(BidirMMapPipe& pipe) { pipe.flush(); return pipe; }
        /// for usage a la "pipe << purge;"
        static BidirMMapPipe& purge(BidirMMapPipe& pipe) { pipe.purge(); return pipe; }

        static int wait_for_child(pid_t child_pid, bool may_throw);

    private:
        /// copy-construction forbidden
        BidirMMapPipe(const BidirMMapPipe&);
        /// assignment forbidden
        BidirMMapPipe& operator=(const BidirMMapPipe&) { return *this; }

        /// page is our friend
        friend class BidirMMapPipe_impl::Page;
        /// convenience typedef for Page
        typedef BidirMMapPipe_impl::Page Page;

        /// tuning constants
        enum {
            // TotPages = 16 will give 32k buffers at 4k page size for both
            // parent and child; if your average message to send is larger
            // than this, consider raising the value (max 256)
            TotPages = 16, ///< pages shared (child + parent)

            PagesPerEnd = TotPages / 2, ///< pages per pipe end

            // if FlushThresh pages are filled, the code forces a flush; 3/4
            // of the pages available seems to work quite well
            FlushThresh = (3 * PagesPerEnd) / 4 ///< flush threshold
        };

        // per-class members
        static pthread_mutex_t s_openpipesmutex; ///< protects s_openpipes
        /// list of open BidirMMapPipes
        static std::list<BidirMMapPipe*> s_openpipes;
        /// pool of mmapped pages
        static BidirMMapPipe_impl::PagePool* s_pagepool;
        /// page pool reference counter
        static unsigned s_pagepoolrefcnt;
        /// debug flag
        static int s_debugflag;

        /// return page pool
        static BidirMMapPipe_impl::PagePool& pagepool();

        // per-instance members
        BidirMMapPipe_impl::Pages m_pages; ///< mmapped pages
        Page* m_busylist; ///< linked list: busy pages (data to be read)
        Page* m_freelist; ///< linked list: free pages
        Page* m_dirtylist; ///< linked list: dirty pages (data to be sent)
        int m_inpipe; ///< pipe end from which data may be read
        int m_outpipe; ///< pipe end to which data may be written
        int m_flags; ///< flags (e.g. end of file)
        pid_t m_childPid; ///< pid of the child (zero if we're child)
        pid_t m_parentPid; ///< pid of the parent
        bool kept_local;


        /// cleanup routine - at exit, we want our children to get a SIGTERM...
        static void teardownall(void);

        /// return length of a page list
        static unsigned lenPageList(const Page* list);

        /** "feed" the busy and free lists with a list of pages
         *
         * @param plist linked list of pages
         *
         * goes through plist, puts free pages from plist onto the freelist
         * (or sends them to the remote end if they belong there), and puts
         * non-empty pages on plist onto the busy list
         */
        void feedPageLists(Page* plist);

        /// put on dirty pages list
        void markPageDirty(Page* p);

        /// transfer bytes through the pipe (reading, writing, may block)
        static size_type xferraw(int fd, void* addr, size_type len,
                ssize_t (*xferfn)(int, void*, std::size_t));
        /// transfer bytes through the pipe (reading, writing, may block)
        static size_type xferraw(int fd, void* addr, const size_type len,
                ssize_t (*xferfn)(int, const void*, std::size_t))
        {
            return xferraw(fd, addr, len,
                    reinterpret_cast<ssize_t (*)(
                        int, void*, std::size_t)>(xferfn));
        }

        /** @brief send page(s) to the other end (may block)
         *
         * @param plist linked list of pages to send
         *
         * the implementation gathers the different write(s) whereever
         * possible; if mmap works, this results in a single write to transfer
         * the list of pages sent, if we need to copy things through the pipe,
         * we have one write to transfer which pages are sent, and then one
         * write per page.
         */
        void sendpages(Page* plist);

        /** @brief receive a pages from the other end (may block), queue them
         *
         * @returns number of pages received
         *
         * this is an application-level scatter read, which gets the list of
         * pages to read from the pipe. if mmap works, it needs only one read
         * call (to get the head of the list of pages transferred). if we need
         * to copy pages through the pipe, we need to add one read for each
         * empty page, and two reads for each non-empty page.
         */
        unsigned recvpages();

        /** @brief receive pages from other end (non-blocking)
         *
         * @returns number of pages received
         *
         * like recvpages(), but does not block if nothing is available for
         * reading
         */
        unsigned recvpages_nonblock();

        /// get a busy page to read data from (may block)
        Page* busypage();
        /// get a dirty page to write data to (may block)
        Page* dirtypage();

        /// close the pipe (no flush if forced)
        int doClose(bool force, bool holdlock = false);
        /// perform the flush
        void doFlush(bool forcePartialPages = true);
#endif //_WIN32
};

END_NAMESPACE_ROOFIT

#undef BEGIN_NAMESPACE_ROOFIT
#undef END_NAMESPACE_ROOFIT

#endif // BIDIRMMAPPIPE_H

// vim: ft=cpp:sw=4:tw=78:et
