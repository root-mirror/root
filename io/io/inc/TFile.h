// @(#)root/io:$Id$
// Author: Rene Brun   28/11/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFile
#define ROOT_TFile


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFile                                                                //
//                                                                      //
// ROOT file.                                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <atomic>
#include <string>

#include "Compression.h"
#include "TDirectoryFile.h"
#include "TUrl.h"
#include "ROOT/RConcurrentHashColl.hxx"

// Not a part of TFile interface; provide a forward declaration instead of #include.
// #ifndef R__LESS_INCLUDES
// #include "TMap.h"
// #endif

#ifdef R__USE_IMT
#include "ROOT/TRWSpinLock.hxx"
#include <mutex>
#endif


class TMap;
class TFree;
class TArrayC;
class TArchiveFile;
class TFileOpenHandle;
class TFileCacheRead;
class TFileCacheWrite;
class TProcessID;
class TStopwatch;
class TFilePrefetch;

class TFile : public TDirectoryFile {
  friend class TDirectoryFile;
  friend class TFilePrefetch;
// TODO: We need to make sure only one TBasket is being written at a time
// if we are writing multiple baskets in parallel.
#ifdef R__USE_IMT
  friend class TBasket;
#endif

public:
   /// Asynchronous open request status
   enum EAsyncOpenStatus { kAOSNotAsync = -1,  kAOSFailure = 0,
                           kAOSInProgress = 1, kAOSSuccess = 2 };
   /// Open timeout constants
   enum EOpenTimeOut { kInstantTimeout = 0, kEternalTimeout = 999999999 };

   /// TTreeCache flushing semantics
   enum ECacheAction { kDisconnect = 0, kDoNotDisconnect = 1 };

protected:
   Double_t         fSumBuffer{0};            ///<Sum of buffer sizes of objects written so far
   Double_t         fSum2Buffer{0};           ///<Sum of squares of buffer sizes of objects written so far
   Long64_t         fBytesWrite{0};           ///<Number of bytes written to this file
   Long64_t         fBytesRead{0};            ///<Number of bytes read from this file
   Long64_t         fBytesReadExtra{0};       ///<Number of extra bytes (overhead) read by the readahead buffer
   Long64_t         fBEGIN{0};                ///<First used byte in file
   Long64_t         fEND{0};                  ///<Last used byte in file
   Long64_t         fSeekFree{0};             ///<Location on disk of free segments structure
   Long64_t         fSeekInfo{0};             ///<Location on disk of StreamerInfo record
   Int_t            fD{-1};                   ///<File descriptor
   Int_t            fVersion{0};              ///<File format version
   Int_t            fCompress{0};             ///<Compression level and algorithm
   Int_t            fNbytesFree{0};           ///<Number of bytes for free segments structure
   Int_t            fNbytesInfo{0};           ///<Number of bytes for StreamerInfo record
   Int_t            fWritten{0};              ///<Number of objects written so far
   Int_t            fNProcessIDs{0};          ///<Number of TProcessID written to this file
   Int_t            fReadCalls{0};            ///<Number of read calls ( not counting the cache calls )
   TString          fRealName;                ///<Effective real file name (not original url)
   TString          fOption;                  ///<File options
   Char_t           fUnits{0};                ///<Number of bytes for file pointers
   TList           *fFree{nullptr};           ///<Free segments linked list table
   TArrayC         *fClassIndex{nullptr};     ///<!Index of TStreamerInfo classes written to this file
   TObjArray       *fProcessIDs{nullptr};     ///<!Array of pointers to TProcessIDs
   Long64_t         fOffset{0};               ///<!Seek offset cache
   TArchiveFile    *fArchive{nullptr};        ///<!Archive file from which we read this file
   TFileCacheRead  *fCacheRead{nullptr};      ///<!Pointer to the read cache (if any)
   TMap            *fCacheReadMap{nullptr};   ///<!Pointer to the read cache (if any)
   TFileCacheWrite *fCacheWrite{nullptr};     ///<!Pointer to the write cache (if any)
   Long64_t         fArchiveOffset{0};        ///<!Offset at which file starts in archive
   Bool_t           fIsArchive{kFALSE};       ///<!True if this is a pure archive file
   Bool_t           fNoAnchorInName{kFALSE};  ///<!True if we don't want to force the anchor to be appended to the file name
   Bool_t           fIsRootFile{kTRUE};       ///<!True is this is a ROOT file, raw file otherwise
   Bool_t           fInitDone{kFALSE};        ///<!True if the file has been initialized
   Bool_t           fMustFlush{kTRUE};        ///<!True if the file buffers must be flushed
   Bool_t           fIsPcmFile{kFALSE};       ///<!True if the file is a ROOT pcm file.
   TFileOpenHandle *fAsyncHandle{nullptr};    ///<!For proper automatic cleanup
   EAsyncOpenStatus fAsyncOpenStatus{kAOSNotAsync}; ///<!Status of an asynchronous open request
   TUrl             fUrl;                     ///<!URL of file

   TList           *fInfoCache{nullptr};      ///<!Cached list of the streamer infos in this file
   TList           *fOpenPhases{nullptr};     ///<!Time info about open phases

#ifdef R__USE_IMT
   std::mutex                                 fWriteMutex;  ///<!Lock for writing baskets / keys into the file.
   static ROOT::Internal::RConcurrentHashColl fgTsSIHashes; ///<!TS Set of hashes built from read streamer infos
#endif

   static TList    *fgAsyncOpenRequests; //List of handles for pending open requests

   static TString   fgCacheFileDir;          ///<Directory where to locally stage files
   static Bool_t    fgCacheFileDisconnected; ///<Indicates, we trust in the files in the cache dir without stat on the cached file
   static Bool_t    fgCacheFileForce;        ///<Indicates, to force all READ to CACHEREAD
   static UInt_t    fgOpenTimeout;           ///<Timeout for open operations in ms  - 0 corresponds to blocking i/o
   static Bool_t    fgOnlyStaged ;           ///<Before the file is opened, it is checked, that the file is staged, if not, the open fails

   static std::atomic<Long64_t>  fgBytesWrite;            ///<Number of bytes written by all TFile objects
   static std::atomic<Long64_t>  fgBytesRead;             ///<Number of bytes read by all TFile objects
   static std::atomic<Long64_t>  fgFileCounter;           ///<Counter for all opened files
   static std::atomic<Int_t>     fgReadCalls;             ///<Number of bytes read from all TFile objects
   static Int_t     fgReadaheadSize;         ///<Readahead buffer size
   static Bool_t    fgReadInfo;              ///<if true (default) ReadStreamerInfo is called when opening a file

   virtual EAsyncOpenStatus GetAsyncOpenStatus() { return fAsyncOpenStatus; }
   virtual void        Init(Bool_t create);
           Bool_t      FlushWriteCache();
           Int_t       ReadBufferViaCache(char *buf, Int_t len);
           Int_t       WriteBufferViaCache(const char *buf, Int_t len);

   ////////////////////////////////////////////////////////////////////////////////
   /// \brief Simple struct of the return value of GetStreamerInfoListImpl
   struct InfoListRet {
      TList *fList;
      Int_t  fReturnCode;
      ROOT::Internal::RConcurrentHashColl::HashValue fHash;
   };

   virtual InfoListRet GetStreamerInfoListImpl(bool lookupSICache);

   // Creating projects
           Int_t       MakeProjectParMake(const char *packname, const char *filename);
           Int_t       MakeProjectParProofInf(const char *packname, const char *proofinfdir);

   // Interface to basic system I/O routines
   virtual Int_t       SysOpen(const char *pathname, Int_t flags, UInt_t mode);
   virtual Int_t       SysClose(Int_t fd);
   virtual Int_t       SysRead(Int_t fd, void *buf, Int_t len);
   virtual Int_t       SysWrite(Int_t fd, const void *buf, Int_t len);
   virtual Long64_t    SysSeek(Int_t fd, Long64_t offset, Int_t whence);
   virtual Int_t       SysStat(Int_t fd, Long_t *id, Long64_t *size, Long_t *flags, Long_t *modtime);
   virtual Int_t       SysSync(Int_t fd);

   // Interface for text-based TDirectory I/O
   virtual Long64_t    DirCreateEntry(TDirectory*) { return 0; }
   virtual Int_t       DirReadKeys(TDirectory*) { return 0; }
   virtual void        DirWriteKeys(TDirectory*) {}
   virtual void        DirWriteHeader(TDirectory*) {}

private:
   TFile(const TFile &) = delete;            //Files cannot be copied
   void operator=(const TFile &) = delete;

   static  void        CpProgress(Long64_t bytesread, Long64_t size, TStopwatch &watch);
   static  TFile      *OpenFromCache(const char *name, Option_t * = "",
                                     const char *ftitle = "", Int_t compress = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault,
                                     Int_t netopt = 0);

public:
   /// TFile status bits. BIT(13) is taken up by TObject
   enum EStatusBits {
      kRecovered     = BIT(10),
      kHasReferences = BIT(11),
      kDevNull       = BIT(12),
      kWriteError    = BIT(14),
      kBinaryFile    = BIT(15),
      kRedirected    = BIT(16),
      kReproducible  = BIT(17)
   };
   enum ERelativeTo { kBeg = 0, kCur = 1, kEnd = 2 };
   enum { kStartBigFile  = 2000000000 };
   /// File type
   enum EFileType { kDefault = 0, kLocal = 1, kNet = 2, kWeb = 3, kFile = 4, kMerge = 5 };

   TFile();
   TFile(const char *fname, Option_t *option="", const char *ftitle="", Int_t compress = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault);
   virtual ~TFile();

           void        Close(Option_t *option="") override; // *MENU*
           void        Copy(TObject &) const override { MayNotUse("Copy(TObject &)"); }
   virtual Bool_t      Cp(const char *dst, Bool_t progressbar = kTRUE,UInt_t buffersize = 1000000);
   virtual TKey*       CreateKey(TDirectory* mother, const TObject* obj, const char* name, Int_t bufsize);
   virtual TKey*       CreateKey(TDirectory* mother, const void* obj, const TClass* cl,
                                 const char* name, Int_t bufsize);
   static TFile      *&CurrentFile(); // Return the current file for this thread.
           void        Delete(const char *namecycle="") override;
           void        Draw(Option_t *option="") override;
   virtual void        DrawMap(const char *keys="*",Option_t *option=""); // *MENU*
           void        FillBuffer(char *&buffer) override;
   virtual void        Flush();
         TArchiveFile *GetArchive() const { return fArchive; }
           Long64_t    GetArchiveOffset() const { return fArchiveOffset; }
           Int_t       GetBestBuffer() const;
   virtual Int_t       GetBytesToPrefetch() const;
       TFileCacheRead *GetCacheRead(const TObject* tree = nullptr) const;
      TFileCacheWrite *GetCacheWrite() const;
           TArrayC    *GetClassIndex() const { return fClassIndex; }
           Int_t       GetCompressionAlgorithm() const;
           Int_t       GetCompressionLevel() const;
           Int_t       GetCompressionSettings() const;
           Float_t     GetCompressionFactor();
   virtual Long64_t    GetEND() const { return fEND; }
   virtual Int_t       GetErrno() const;
   virtual void        ResetErrno() const;
           Int_t       GetFd() const { return fD; }
   virtual const TUrl *GetEndpointUrl() const { return &fUrl; }
           TObjArray  *GetListOfProcessIDs() const {return fProcessIDs;}
           TList      *GetListOfFree() const { return fFree; }
   virtual Int_t       GetNfree() const { return fFree->GetSize(); }
   virtual Int_t       GetNProcessIDs() const { return fNProcessIDs; }
           Option_t   *GetOption() const override { return fOption.Data(); }
   virtual Long64_t    GetBytesRead() const { return fBytesRead; }
   virtual Long64_t    GetBytesReadExtra() const { return fBytesReadExtra; }
   virtual Long64_t    GetBytesWritten() const;
   virtual Int_t       GetReadCalls() const { return fReadCalls; }
           Int_t       GetVersion() const { return fVersion; }
           Int_t       GetRecordHeader(char *buf, Long64_t first, Int_t maxbytes,
                                       Int_t &nbytes, Int_t &objlen, Int_t &keylen);
   virtual Int_t       GetNbytesInfo() const {return fNbytesInfo;}
   virtual Int_t       GetNbytesFree() const {return fNbytesFree;}
   virtual TString     GetNewUrl() { return ""; }
           Long64_t    GetRelOffset() const { return fOffset - fArchiveOffset; }
   virtual Long64_t    GetSeekFree() const {return fSeekFree;}
   virtual Long64_t    GetSeekInfo() const {return fSeekInfo;}
   virtual Long64_t    GetSize() const;
   virtual TList      *GetStreamerInfoList() final; // Note: to override behavior, please override GetStreamerInfoListImpl
   const   TList      *GetStreamerInfoCache();
   virtual void        IncrementProcessIDs() { fNProcessIDs++; }
   virtual Bool_t      IsArchive() const { return fIsArchive; }
           Bool_t      IsBinary() const { return TestBit(kBinaryFile); }
           Bool_t      IsRaw() const { return !fIsRootFile; }
   virtual Bool_t      IsOpen() const;
           void        ls(Option_t *option="") const override;
   virtual void        MakeFree(Long64_t first, Long64_t last);
   virtual void        MakeProject(const char *dirname, const char *classes="*",
                                   Option_t *option="new"); // *MENU*
   virtual void        Map(Option_t *opt); // *MENU*
   virtual void        Map() { Map(""); }; // *MENU*
   virtual Bool_t      Matches(const char *name);
   virtual Bool_t      MustFlush() const {return fMustFlush;}
           void        Paint(Option_t *option="") override;
           void        Print(Option_t *option="") const override;
   virtual Bool_t      ReadBufferAsync(Long64_t offs, Int_t len);
   virtual Bool_t      ReadBuffer(char *buf, Int_t len);
   virtual Bool_t      ReadBuffer(char *buf, Long64_t pos, Int_t len);
   virtual Bool_t      ReadBuffers(char *buf, Long64_t *pos, Int_t *len, Int_t nbuf);
   virtual void        ReadFree();
   virtual TProcessID *ReadProcessID(UShort_t pidf);
   virtual void        ReadStreamerInfo();
   virtual Int_t       Recover();
   virtual Int_t       ReOpen(Option_t *mode);
   virtual void        Seek(Long64_t offset, ERelativeTo pos = kBeg);
   virtual void        SetCacheRead(TFileCacheRead *cache, TObject* tree = 0, ECacheAction action = kDisconnect);
   virtual void        SetCacheWrite(TFileCacheWrite *cache);
   virtual void        SetCompressionAlgorithm(Int_t algorithm = ROOT::RCompressionSetting::EAlgorithm::kUseGlobal);
   virtual void        SetCompressionLevel(Int_t level = ROOT::RCompressionSetting::ELevel::kUseMin);
   virtual void        SetCompressionSettings(Int_t settings = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault);
   virtual void        SetEND(Long64_t last) { fEND = last; }
   virtual void        SetOffset(Long64_t offset, ERelativeTo pos = kBeg);
   virtual void        SetOption(Option_t *option=">") { fOption = option; }
   virtual void        SetReadCalls(Int_t readcalls = 0) { fReadCalls = readcalls; }
   virtual void        ShowStreamerInfo();
           Int_t       Sizeof() const override;
           void        SumBuffer(Int_t bufsize);
   virtual Bool_t      WriteBuffer(const char *buf, Int_t len);
           Int_t       Write(const char *name=nullptr, Int_t opt=0, Int_t bufsiz=0) override;
           Int_t       Write(const char *name=nullptr, Int_t opt=0, Int_t bufsiz=0) const override;
   virtual void        WriteFree();
   virtual void        WriteHeader();
   virtual UShort_t    WriteProcessID(TProcessID *pid);
   virtual void        WriteStreamerInfo();

   static TFileOpenHandle
                      *AsyncOpen(const char *name, Option_t *option = "",
                                 const char *ftitle = "", Int_t compress = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault,
                                 Int_t netopt = 0);
   static TFile       *Open(const char *name, Option_t *option = "",
                            const char *ftitle = "", Int_t compress = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault,
                            Int_t netopt = 0);
   static TFile       *Open(TFileOpenHandle *handle);

   static EFileType    GetType(const char *name, Option_t *option = "", TString *prefix = nullptr);

   static EAsyncOpenStatus GetAsyncOpenStatus(const char *name);
   static EAsyncOpenStatus GetAsyncOpenStatus(TFileOpenHandle *handle);
   static const TUrl  *GetEndpointUrl(const char *name);

   static Long64_t     GetFileBytesRead();
   static Long64_t     GetFileBytesWritten();
   static Int_t        GetFileReadCalls();
   static Int_t        GetReadaheadSize();

   static void         SetFileBytesRead(Long64_t bytes = 0);
   static void         SetFileBytesWritten(Long64_t bytes = 0);
   static void         SetFileReadCalls(Int_t readcalls = 0);
   static void         SetReadaheadSize(Int_t bufsize = 256000);
   static void         SetReadStreamerInfo(Bool_t readinfo=kTRUE);
   static Bool_t       GetReadStreamerInfo();

   static Long64_t     GetFileCounter();
   static void         IncrementFileCounter();

   static Bool_t       SetCacheFileDir(ROOT::Internal::TStringView cacheDir, Bool_t operateDisconnected = kTRUE,
                                       Bool_t forceCacheread = kFALSE)
     { return SetCacheFileDir(std::string_view(cacheDir), operateDisconnected, forceCacheread); }
   static Bool_t       SetCacheFileDir(std::string_view cacheDir, Bool_t operateDisconnected = kTRUE,
                                       Bool_t forceCacheread = kFALSE);
   static const char  *GetCacheFileDir();
   static Bool_t       ShrinkCacheFileDir(Long64_t shrinkSize, Long_t cleanupInteval = 0);
   static Bool_t       Cp(const char *src, const char *dst, Bool_t progressbar = kTRUE,
                          UInt_t buffersize = 1000000);

   static UInt_t       SetOpenTimeout(UInt_t timeout);  // in ms
   static UInt_t       GetOpenTimeout(); // in ms
   static Bool_t       SetOnlyStaged(Bool_t onlystaged);
   static Bool_t       GetOnlyStaged();

   ClassDefOverride(TFile,8)  //ROOT file
};

#ifndef __CINT__
#define gFile (TFile::CurrentFile())

#elif defined(__MAKECINT__)
// To properly handle the use of gFile in header files (in static declarations)
R__EXTERN TFile   *gFile;
#endif

/**
\class TFileOpenHandle
\ingroup IO
Class holding info about the file being opened
*/
class TFileOpenHandle : public TNamed {

friend class TFile;
friend class TAlienFile;

private:
   TString  fOpt;            ///< Options
   Int_t    fCompress{0};    ///< Compression level and algorithm
   Int_t    fNetOpt{0};      ///< Network options
   TFile   *fFile{nullptr};  ///< TFile instance of the file being opened

   TFileOpenHandle(TFile *f) : TNamed("",""), fOpt(""), fCompress(ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault),
                               fNetOpt(0), fFile(f) { }
   TFileOpenHandle(const char *n, const char *o, const char *t, Int_t cmp,
                   Int_t no) : TNamed(n,t), fOpt(o), fCompress(cmp),
                               fNetOpt(no), fFile(nullptr) { }
   TFileOpenHandle(const TFileOpenHandle&) = delete;
   TFileOpenHandle& operator=(const TFileOpenHandle&) = delete;

   TFile      *GetFile() const { return fFile; }

public:
   ~TFileOpenHandle() { }

   Bool_t      Matches(const char *name);

   const char *GetOpt() const { return fOpt; }
   Int_t       GetCompress() const { return fCompress; }
   Int_t       GetNetOpt() const { return fNetOpt; }
};

//______________________________________________________________________________
inline Int_t TFile::GetCompressionAlgorithm() const
{
   return (fCompress < 0) ? -1 : fCompress / 100;
}

//______________________________________________________________________________
inline Int_t TFile::GetCompressionLevel() const
{
   return (fCompress < 0) ? -1 : fCompress % 100;
}

//______________________________________________________________________________
inline Int_t TFile::GetCompressionSettings() const
{
   return (fCompress < 0) ? -1 : fCompress;
}

#endif
