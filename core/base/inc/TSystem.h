// @(#)root/base:$Id$
// Author: Fons Rademakers   15/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TSystem
#define ROOT_TSystem


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TSystem                                                              //
//                                                                      //
// Abstract base class defining a generic interface to the underlying   //
// Operating System.                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include <cstdio>
#include <cctype>
#include <fcntl.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#endif

#include "TNamed.h"
#include "TInetAddress.h"
#include "TTimer.h"
#include <string>

class TSeqCollection;
class TFdSet;
class TVirtualMutex;

enum EAccessMode {
   kFileExists        = 0,
   kExecutePermission = 1,
   kWritePermission   = 2,
   kReadPermission    = 4
};

enum ELogOption {
   kLogPid            = 0x01,
   kLogCons           = 0x02
};

enum ELogLevel {
   kLogEmerg          = 0,
   kLogAlert          = 1,
   kLogCrit           = 2,
   kLogErr            = 3,
   kLogWarning        = 4,
   kLogNotice         = 5,
   kLogInfo           = 6,
   kLogDebug          = 7
};

enum ELogFacility {
   kLogLocal0,
   kLogLocal1,
   kLogLocal2,
   kLogLocal3,
   kLogLocal4,
   kLogLocal5,
   kLogLocal6,
   kLogLocal7
};

enum EFpeMask {
   kNoneMask         = 0x00,
   kInvalid          = 0x01,  // Invalid argument
   kDivByZero        = 0x02,  // Division by zero
   kOverflow         = 0x04,  // Overflow
   kUnderflow        = 0x08,  // Underflow
   kInexact          = 0x10,  // Inexact
   kDefaultMask      = 0x07,
   kAllMask          = 0x1F
};

enum EFileModeMask {
   kS_IFMT   = 0170000,   // bitmask for the file type bitfields
   kS_IFSOCK = 0140000,   // socket
   kS_IFLNK  = 0120000,   // symbolic link
   kS_IFOFF  = 0110000,   // offline file
   kS_IFREG  = 0100000,   // regular file
   kS_IFBLK  = 0060000,   // block device
   kS_IFDIR  = 0040000,   // directory
   kS_IFCHR  = 0020000,   // character device
   kS_IFIFO  = 0010000,   // fifo
   kS_ISUID  = 0004000,   // set UID bit
   kS_ISGID  = 0002000,   // set GID bit
   kS_ISVTX  = 0001000,   // sticky bit
   kS_IRWXU  = 00700,     // mask for file owner permissions
   kS_IRUSR  = 00400,     // owner has read permission
   kS_IWUSR  = 00200,     // owner has write permission
   kS_IXUSR  = 00100,     // owner has execute permission
   kS_IRWXG  = 00070,     // mask for group permissions
   kS_IRGRP  = 00040,     // group has read permission
   kS_IWGRP  = 00020,     // group has write permission
   kS_IXGRP  = 00010,     // group has execute permission
   kS_IRWXO  = 00007,     // mask for permissions for others (not in group)
   kS_IROTH  = 00004,     // others have read permission
   kS_IWOTH  = 00002,     // others have write permission
   kS_IXOTH  = 00001      // others have execute permission
};

inline Bool_t R_ISDIR(Int_t mode)  { return ((mode & kS_IFMT) == kS_IFDIR); }
inline Bool_t R_ISCHR(Int_t mode)  { return ((mode & kS_IFMT) == kS_IFCHR); }
inline Bool_t R_ISBLK(Int_t mode)  { return ((mode & kS_IFMT) == kS_IFBLK); }
inline Bool_t R_ISREG(Int_t mode)  { return ((mode & kS_IFMT) == kS_IFREG); }
inline Bool_t R_ISLNK(Int_t mode)  { return ((mode & kS_IFMT) == kS_IFLNK); }
inline Bool_t R_ISFIFO(Int_t mode) { return ((mode & kS_IFMT) == kS_IFIFO); }
inline Bool_t R_ISSOCK(Int_t mode) { return ((mode & kS_IFMT) == kS_IFSOCK); }
inline Bool_t R_ISOFF(Int_t mode)  { return ((mode & kS_IFMT) == kS_IFOFF); }

struct FileStat_t {
   Long_t   fDev;          // device id
   Long_t   fIno;          // inode
   Int_t    fMode;         // protection (combination of EFileModeMask bits)
   Int_t    fUid;          // user id of owner
   Int_t    fGid;          // group id of owner
   Long64_t fSize;         // total size in bytes
   Long_t   fMtime;        // modification date
   Bool_t   fIsLink;       // symbolic link
   TString  fUrl;          // end point url of file
   FileStat_t() : fDev(0), fIno(0), fMode(0), fUid(0), fGid(0), fSize(0),
                  fMtime(0), fIsLink(kFALSE), fUrl("") { }
};

struct UserGroup_t {
   Int_t    fUid;          // user id
   Int_t    fGid;          // group id
   TString  fUser;         // user name
   TString  fGroup;        // group name
   TString  fPasswd;       // password
   TString  fRealName;     // user full name
   TString  fShell;        // user preferred shell
   UserGroup_t() : fUid(0), fGid(0), fUser(), fGroup(), fPasswd(),
                   fRealName (), fShell() { }
};

struct SysInfo_t {
   TString   fOS;          // OS
   TString   fModel;       // computer model
   TString   fCpuType;     // type of cpu
   Int_t     fCpus;        // number of cpus
   Int_t     fCpuSpeed;    // cpu speed in MHz
   Int_t     fBusSpeed;    // bus speed in MHz
   Int_t     fL2Cache;     // level 2 cache size in KB
   Int_t     fPhysRam;     // physical RAM in MB
   SysInfo_t() : fOS(), fModel(), fCpuType(), fCpus(0), fCpuSpeed(0),
                 fBusSpeed(0), fL2Cache(0), fPhysRam(0) { }
   virtual ~SysInfo_t() { }
   ClassDef(SysInfo_t, 1); // System information - OS, CPU, RAM.
};

struct CpuInfo_t {
   Float_t   fLoad1m;      // cpu load average over 1 m
   Float_t   fLoad5m;      // cpu load average over 5 m
   Float_t   fLoad15m;     // cpu load average over 15 m
   Float_t   fUser;        // cpu user load in percentage
   Float_t   fSys;         // cpu sys load in percentage
   Float_t   fTotal;       // cpu user+sys load in percentage
   Float_t   fIdle;        // cpu idle percentage
   CpuInfo_t() : fLoad1m(0), fLoad5m(0), fLoad15m(0),
                 fUser(0), fSys(0), fTotal(0), fIdle(0) { }
   virtual ~CpuInfo_t() { }
   ClassDef(CpuInfo_t, 1); // CPU load information.
};

struct MemInfo_t {
   Int_t     fMemTotal;    // total RAM in MB
   Int_t     fMemUsed;     // used RAM in MB
   Int_t     fMemFree;     // free RAM in MB
   Int_t     fSwapTotal;   // total swap in MB
   Int_t     fSwapUsed;    // used swap in MB
   Int_t     fSwapFree;    // free swap in MB
   MemInfo_t() : fMemTotal(0), fMemUsed(0), fMemFree(0),
                 fSwapTotal(0), fSwapUsed(0), fSwapFree(0) { }
   virtual ~MemInfo_t() { }
   ClassDef(MemInfo_t, 1); // Memory utilization information.
};

struct ProcInfo_t {
   Float_t   fCpuUser;     // user time used by this process in seconds
   Float_t   fCpuSys;      // system time used by this process in seconds
   Long_t    fMemResident; // resident memory used by this process in KB
   Long_t    fMemVirtual;  // virtual memory used by this process in KB
   ProcInfo_t() : fCpuUser(0), fCpuSys(0), fMemResident(0),
                  fMemVirtual(0) { }
   virtual ~ProcInfo_t();
   ClassDef(ProcInfo_t, 1);// System resource usage of given process.
};

struct RedirectHandle_t {
   TString   fFile;            // File where the output was redirected
   TString   fStdOutTty;       // tty associated with stdout, if any (e.g. from ttyname(...))
   TString   fStdErrTty;       // tty associated with stderr, if any (e.g. from ttyname(...))
   Int_t     fStdOutDup{-1};   // Duplicated descriptor for stdout
   Int_t     fStdErrDup{-1};   // Duplicated descriptor for stderr
   Int_t     fReadOffSet{-1};  // Offset where to start reading the file (used by ShowOutput(...))
   RedirectHandle_t(const char *n = nullptr) : fFile(n) { }
   void Reset() { fFile = ""; fStdOutTty = ""; fStdErrTty = "";
                  fStdOutDup = -1; fStdErrDup = -1; fReadOffSet = -1; }
};

enum ESockOptions {
   kSendBuffer,        // size of send buffer
   kRecvBuffer,        // size of receive buffer
   kOobInline,         // OOB message inline
   kKeepAlive,         // keep socket alive
   kReuseAddr,         // allow reuse of local portion of address 5-tuple
   kNoDelay,           // send without delay
   kNoBlock,           // non-blocking I/O
   kProcessGroup,      // socket process group (used for SIGURG and SIGIO)
   kAtMark,            // are we at out-of-band mark (read only)
   kBytesToRead        // get number of bytes to read, FIONREAD (read only)
};

enum ESendRecvOptions {
   kDefault,           // default option (= 0)
   kOob,               // send or receive out-of-band data
   kPeek,              // peek at incoming message (receive only)
   kDontBlock          // send/recv as much data as possible without blocking
};

#ifdef __CINT__
typedef void *Func_t;
#else
typedef void ((*Func_t)());
#endif

R__EXTERN const char  *gRootDir;
R__EXTERN const char  *gProgName;
R__EXTERN const char  *gProgPath;
R__EXTERN TVirtualMutex *gSystemMutex;


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Asynchronous timer used for processing pending GUI and timer events  //
// every delay ms. Call in a tight computing loop                       //
// TProcessEventTimer::ProcessEvent(). If the timer did timeout this    //
// call will process the pending events and return kTRUE if the         //
// TROOT::IsInterrupted() flag is set (can be done by hitting key in    //
// canvas or selecting canvas menu item View/Interrupt.                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TProcessEventTimer : public TTimer {
public:
   TProcessEventTimer(Long_t delay);
   Bool_t Notify() { return kTRUE; }
   Bool_t ProcessEvents();
   ClassDef(TProcessEventTimer,0)  // Process pending events at fixed time intervals
};


class TSystem : public TNamed {

public:
   enum EAclicMode { kDefault, kDebug, kOpt };
   enum EAclicProperties {
      kFlatBuildDir = BIT(0)           // If set and a BuildDir is selected, then do not created sub-directories
   };

protected:
   TFdSet          *fReadmask{nullptr};         //!Files that should be checked for read events
   TFdSet          *fWritemask{nullptr};        //!Files that should be checked for write events
   TFdSet          *fReadready{nullptr};        //!Files with reads waiting
   TFdSet          *fWriteready{nullptr};       //!Files with writes waiting
   TFdSet          *fSignals{nullptr};          //!Signals that were trapped
   Int_t            fNfd{0};                    //Number of fd's in masks
   Int_t            fMaxrfd{-1};                //Largest fd in read mask
   Int_t            fMaxwfd{-1};                //Largest fd in write mask
   Int_t            fSigcnt{0};                 //Number of pending signals
   TString          fWdpath;                    //Working directory
   TString          fHostname;                  //Hostname
   Bool_t           fInsideNotify{kFALSE};      //Used by DispatchTimers()
   Int_t            fBeepFreq{0};               //Used by Beep()
   Int_t            fBeepDuration{0};           //Used by Beep()

   Bool_t           fInControl{kFALSE};         //True if in eventloop
   Bool_t           fDone{kFALSE};              //True if eventloop should be finished
   Int_t            fLevel{0};                  //Level of nested eventloops

   TSeqCollection  *fTimers{nullptr};           //List of timers
   TSeqCollection  *fSignalHandler{nullptr};    //List of signal handlers
   TSeqCollection  *fFileHandler{nullptr};      //List of file handlers
   TSeqCollection  *fStdExceptionHandler{nullptr}; //List of std::exception handlers
   TSeqCollection  *fOnExitList{nullptr};       //List of items to be cleaned-up on exit

   TString          fListLibs;                  //List shared libraries, cache used by GetLibraries

   TString          fBuildArch;                 //Architecture for which ROOT was built (passed to ./configure)
   TString          fBuildCompiler;             //Compiler used to build this ROOT
   TString          fBuildCompilerVersion;      //Compiler version used to build this ROOT
   TString          fBuildCompilerVersionStr;   //Compiler version identifier string used to build this ROOT
   TString          fBuildNode;                 //Detailed information where ROOT was built
   TString          fBuildDir;                  //Location where to build ACLiC shared library and use as scratch area.
   TString          fFlagsDebug;                //Flags for debug compilation
   TString          fFlagsOpt;                  //Flags for optimized compilation
   TString          fListPaths;                 //List of all include (fIncludePath + interpreter include path). Cache used by GetIncludePath
   TString          fIncludePath;               //Used to expand $IncludePath in the directives given to SetMakeSharedLib and SetMakeExe
   TString          fLinkedLibs;                //Used to expand $LinkedLibs in the directives given to SetMakeSharedLib and SetMakeExe
   TString          fSoExt;                     //Extension of shared library (.so, .sl, .a, .dll, etc.)
   TString          fObjExt;                    //Extension of object files (.o, .obj, etc.)
   EAclicMode       fAclicMode{kDefault};       //Whether the compilation should be done debug or opt
   TString          fMakeSharedLib;             //Directive used to build a shared library
   TString          fMakeExe;                   //Directive used to build an executable
   TString          fLinkdefSuffix;             //Default suffix for linkdef files to be used by ACLiC (see EACLiCProperties)
   Int_t            fAclicProperties{0};        //Various boolean flag for change ACLiC's behavior.
   TSeqCollection  *fCompiled{nullptr};         //List of shared libs from compiled macros to be deleted
   TSeqCollection  *fHelpers{nullptr};          //List of helper classes for alternative file/directory access

   TString &GetLastErrorString();             //Last system error message (thread local).
   const TString &GetLastErrorString() const; //Last system error message (thread local).

   TSystem               *FindHelper(const char *path, void *dirptr = nullptr);
   virtual Bool_t         ConsistentWith(const char *path, void *dirptr = nullptr);
   virtual const char    *ExpandFileName(const char *fname);
   virtual Bool_t         ExpandFileName(TString &fname);
   virtual void           SigAlarmInterruptsSyscalls(Bool_t) { }
   virtual const char    *GetLinkedLibraries();
   virtual void           DoBeep(Int_t /*freq*/=-1, Int_t /*duration*/=-1) const { printf("\a"); fflush(stdout); }

   static const char     *StripOffProto(const char *path, const char *proto);

private:
   TSystem(const TSystem&) = delete;
   TSystem& operator=(const TSystem&) = delete;
   Bool_t ExpandFileName(const char *fname, char *xname, const int kBufSize);

public:
   TSystem(const char *name = "Generic", const char *title = "Generic System");
   virtual ~TSystem();

   //---- Misc
   virtual Bool_t          Init();
   virtual void            SetProgname(const char *name);
   virtual void            SetDisplay();
   void                    SetErrorStr(const char *errstr);
   const char             *GetErrorStr() const { return GetLastErrorString(); }
   virtual const char     *GetError();
   void                    RemoveOnExit(TObject *obj);
   virtual const char     *HostName();
   virtual void            NotifyApplicationCreated();

   static Int_t            GetErrno();
   static void             ResetErrno();
   void                    Beep(Int_t freq=-1, Int_t duration=-1, Bool_t setDefault=kFALSE);
   void                    GetBeepDefaults(Int_t &freq, Int_t &duration) const { freq = fBeepFreq; duration = fBeepDuration; }

   //---- EventLoop
   virtual void            Run();
   virtual Bool_t          ProcessEvents();
   virtual void            DispatchOneEvent(Bool_t pendingOnly = kFALSE);
   virtual void            ExitLoop();
   Bool_t                  InControl() const { return fInControl; }
   virtual void            InnerLoop();
   virtual Int_t           Select(TList *active, Long_t timeout);
   virtual Int_t           Select(TFileHandler *fh, Long_t timeout);

   //---- Handling of system events
   virtual void            AddSignalHandler(TSignalHandler *sh);
   virtual TSignalHandler *RemoveSignalHandler(TSignalHandler *sh);
   virtual void            ResetSignal(ESignals sig, Bool_t reset = kTRUE);
   virtual void            ResetSignals();
   virtual void            IgnoreSignal(ESignals sig, Bool_t ignore = kTRUE);
   virtual void            IgnoreInterrupt(Bool_t ignore = kTRUE);
   virtual TSeqCollection *GetListOfSignalHandlers() const { return fSignalHandler; }
   virtual void            AddFileHandler(TFileHandler *fh);
   virtual TFileHandler   *RemoveFileHandler(TFileHandler *fh);
   virtual TSeqCollection *GetListOfFileHandlers() const { return fFileHandler; }
   virtual void            AddStdExceptionHandler(TStdExceptionHandler *eh);
   virtual TStdExceptionHandler *RemoveStdExceptionHandler(TStdExceptionHandler *eh);
   virtual TSeqCollection *GetListOfStdExceptionHandlers() const { return fStdExceptionHandler; }

   //---- Floating Point Exceptions Control
   virtual Int_t           GetFPEMask();
   virtual Int_t           SetFPEMask(Int_t mask = kDefaultMask);

   //---- Time & Date
   virtual TTime           Now();
   virtual TSeqCollection *GetListOfTimers() const { return fTimers; }
   virtual void            AddTimer(TTimer *t);
   virtual TTimer         *RemoveTimer(TTimer *t);
   virtual void            ResetTimer(TTimer *) { }
   virtual Long_t          NextTimeOut(Bool_t mode);
   virtual void            Sleep(UInt_t milliSec);

   //---- Processes
   virtual Int_t           Exec(const char *shellcmd);
   virtual FILE           *OpenPipe(const char *command, const char *mode);
   virtual int             ClosePipe(FILE *pipe);
   virtual TString         GetFromPipe(const char *command);
   virtual void            Exit(int code, Bool_t mode = kTRUE);
   virtual void            Abort(int code = 0);
   virtual int             GetPid();
   virtual void            StackTrace();

   //---- Directories
   virtual int             MakeDirectory(const char *name);
   virtual void           *OpenDirectory(const char *name);
   virtual void            FreeDirectory(void *dirp);
   virtual const char     *GetDirEntry(void *dirp);
   virtual void           *GetDirPtr() const { return nullptr; }
   virtual Bool_t          ChangeDirectory(const char *path);
   virtual const char     *WorkingDirectory();
   virtual std::string     GetWorkingDirectory() const;
   virtual const char     *HomeDirectory(const char *userName = nullptr);
   virtual std::string     GetHomeDirectory(const char *userName = nullptr) const;
   virtual int             mkdir(const char *name, Bool_t recursive = kFALSE);
   Bool_t                  cd(const char *path) { return ChangeDirectory(path); }
   const char             *pwd() { return WorkingDirectory(); }
   virtual const char     *TempDirectory() const;
   virtual FILE           *TempFileName(TString &base, const char *dir = nullptr);

   //---- Paths & Files
   virtual const char     *BaseName(const char *pathname);
   virtual const char     *DirName(const char *pathname);
   virtual TString         GetDirName(const char *pathname);
   virtual char           *ConcatFileName(const char *dir, const char *name);
   virtual Bool_t          IsAbsoluteFileName(const char *dir);
   virtual Bool_t          IsFileInIncludePath(const char *name, char **fullpath = nullptr);
   virtual const char     *PrependPathName(const char *dir, TString& name);
   virtual Bool_t          ExpandPathName(TString &path);
   virtual char           *ExpandPathName(const char *path);
   virtual Bool_t          AccessPathName(const char *path, EAccessMode mode = kFileExists);
   virtual Bool_t          IsPathLocal(const char *path);
   virtual int             CopyFile(const char *from, const char *to, Bool_t overwrite = kFALSE);
   virtual int             Rename(const char *from, const char *to);
   virtual int             Link(const char *from, const char *to);
   virtual int             Symlink(const char *from, const char *to);
   virtual int             Unlink(const char *name);
   int                     GetPathInfo(const char *path, Long_t *id, Long_t *size, Long_t *flags, Long_t *modtime);
   int                     GetPathInfo(const char *path, Long_t *id, Long64_t *size, Long_t *flags, Long_t *modtime);
   virtual int             GetPathInfo(const char *path, FileStat_t &buf);
   virtual int             GetFsInfo(const char *path, Long_t *id, Long_t *bsize, Long_t *blocks, Long_t *bfree);
   virtual int             Chmod(const char *file, UInt_t mode);
   virtual int             Umask(Int_t mask);
   virtual int             Utime(const char *file, Long_t modtime, Long_t actime);
   virtual const char     *UnixPathName(const char *unixpathname);
   virtual const char     *FindFile(const char *search, TString& file, EAccessMode mode = kFileExists);
   virtual char           *Which(const char *search, const char *file, EAccessMode mode = kFileExists);
   virtual TList          *GetVolumes(Option_t *) const { return 0; }

   //---- Users & Groups
   virtual Int_t           GetUid(const char *user = nullptr);
   virtual Int_t           GetGid(const char *group = nullptr);
   virtual Int_t           GetEffectiveUid();
   virtual Int_t           GetEffectiveGid();
   virtual UserGroup_t    *GetUserInfo(Int_t uid);
   virtual UserGroup_t    *GetUserInfo(const char *user = nullptr);
   virtual UserGroup_t    *GetGroupInfo(Int_t gid);
   virtual UserGroup_t    *GetGroupInfo(const char *group = nullptr);

   //---- Environment Manipulation
   virtual void            Setenv(const char *name, const char *value);
   virtual void            Unsetenv(const char *name);
   virtual const char     *Getenv(const char *env);

   //---- System Logging
   virtual void            Openlog(const char *name, Int_t options, ELogFacility facility);
   virtual void            Syslog(ELogLevel level, const char *mess);
   virtual void            Closelog();

   //---- Standard Output redirection
   virtual Int_t           RedirectOutput(const char *name, const char *mode = "a", RedirectHandle_t *h = nullptr);
   virtual void            ShowOutput(RedirectHandle_t *h);

   //---- Dynamic Loading
   virtual void            AddDynamicPath(const char *pathname);
   virtual const char     *GetDynamicPath();
   virtual void            SetDynamicPath(const char *pathname);
   char                   *DynamicPathName(const char *lib, Bool_t quiet = kFALSE);
   virtual const char     *FindDynamicLibrary(TString& lib, Bool_t quiet = kFALSE);
   virtual Func_t          DynFindSymbol(const char *module, const char *entry);
   virtual int             Load(const char *module, const char *entry = "", Bool_t system = kFALSE);
   virtual void            Unload(const char *module);
   virtual UInt_t          LoadAllLibraries();
   virtual void            ListSymbols(const char *module, const char *re = "");
   virtual void            ListLibraries(const char *regexp = "");
   virtual const char     *GetLibraries(const char *regexp = "",
                                        const char *option = "",
                                        Bool_t isRegexp = kTRUE);

   //---- RPC
   virtual TInetAddress    GetHostByName(const char *server);
   virtual TInetAddress    GetPeerName(int sock);
   virtual TInetAddress    GetSockName(int sock);
   virtual int             GetServiceByName(const char *service);
   virtual char           *GetServiceByPort(int port);
   virtual int             OpenConnection(const char *server, int port, int tcpwindowsize = -1, const char *protocol = "tcp");
   virtual int             AnnounceTcpService(int port, Bool_t reuse, int backlog, int tcpwindowsize = -1);
   virtual int             AnnounceUdpService(int port, int backlog);
   virtual int             AnnounceUnixService(int port, int backlog);
   virtual int             AnnounceUnixService(const char *sockpath, int backlog);
   virtual int             AcceptConnection(int sock);
   virtual void            CloseConnection(int sock, Bool_t force = kFALSE);
   virtual int             RecvRaw(int sock, void *buffer, int length, int flag);
   virtual int             SendRaw(int sock, const void *buffer, int length, int flag);
   virtual int             RecvBuf(int sock, void *buffer, int length);
   virtual int             SendBuf(int sock, const void *buffer, int length);
   virtual int             SetSockOpt(int sock, int kind, int val);
   virtual int             GetSockOpt(int sock, int kind, int *val);

   //---- System, CPU and Memory info
   virtual int             GetSysInfo(SysInfo_t *info) const;
   virtual int             GetCpuInfo(CpuInfo_t *info, Int_t sampleTime = 1000) const;
   virtual int             GetMemInfo(MemInfo_t *info) const;
   virtual int             GetProcInfo(ProcInfo_t *info) const;

   //---- ACLiC (Automatic Compiler of Shared Library for CINT)
   virtual void            AddIncludePath(const char *includePath);
   virtual void            AddLinkedLibs(const char *linkedLib);
   virtual int             CompileMacro(const char *filename, Option_t *opt = "", const char *library_name = "", const char *build_dir = "", UInt_t dirmode = 0);
   virtual Int_t           GetAclicProperties() const;
   virtual const char     *GetBuildArch() const;
   virtual const char     *GetBuildCompiler() const;
   virtual const char     *GetBuildCompilerVersion() const;
   virtual const char     *GetBuildCompilerVersionStr() const;
   virtual const char     *GetBuildNode() const;
   virtual const char     *GetBuildDir() const;
   virtual const char     *GetFlagsDebug() const;
   virtual const char     *GetFlagsOpt() const;
   virtual const char     *GetIncludePath();
   virtual const char     *GetLinkedLibs() const;
   virtual const char     *GetLinkdefSuffix() const;
   virtual EAclicMode      GetAclicMode() const;
   virtual const char     *GetMakeExe() const;
   virtual const char     *GetMakeSharedLib() const;
   virtual const char     *GetSoExt() const;
   virtual const char     *GetObjExt() const;
   virtual void            SetBuildDir(const char* build_dir, Bool_t isflat = kFALSE);
   virtual void            SetFlagsDebug(const char *);
   virtual void            SetFlagsOpt(const char *);
   virtual void            SetIncludePath(const char *includePath);
   virtual void            SetMakeExe(const char *directives);
   virtual void            SetAclicMode(EAclicMode mode);
   virtual void            SetMakeSharedLib(const char *directives);
   virtual void            SetLinkedLibs(const char *linkedLibs);
   virtual void            SetLinkdefSuffix(const char *suffix);
   virtual void            SetSoExt(const char *soExt);
   virtual void            SetObjExt(const char *objExt);
   virtual TString         SplitAclicMode(const char *filename, TString &mode, TString &args, TString &io) const;
   virtual void            CleanCompiledMacros();

   ClassDef(TSystem,0)  //ABC defining a generic interface to the OS
};

R__EXTERN TSystem *gSystem;
R__EXTERN TFileHandler *gXDisplay;  // Display server (X11) input event handler


#endif
