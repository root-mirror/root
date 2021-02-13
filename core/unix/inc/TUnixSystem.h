// @(#)root/unix:$Id$
// Author: Fons Rademakers   15/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TUnixSystem
#define ROOT_TUnixSystem


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TUnixSystem                                                          //
//                                                                      //
// Class providing an interface to the UNIX Operating System.           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TSystem.h"
#include "TSysEvtHandler.h"
#include "TTimer.h"
#include <string>

typedef void (*SigHandler_t)(ESignals);


class TUnixSystem : public TSystem {

private:
   void FillWithCwd(char *cwd) const;

protected:
   const char    *GetLinkedLibraries() override;

   // static functions providing semi-low level interface to raw Unix
   static int          UnixMakedir(const char *name);
   static void        *UnixOpendir(const char *name);
   static const char  *UnixGetdirentry(void *dir);
   static const char  *UnixHomedirectory(const char *user = nullptr);
   static const char  *UnixHomedirectory(const char *user, char *path, char *mydir);
   static Long64_t     UnixNow();
   static int          UnixWaitchild();
   static int          UnixSetitimer(Long_t ms);
   static int          UnixSelect(Int_t nfds, TFdSet *readready, TFdSet *writeready,
                                  Long_t timeout);
   static void         UnixSignal(ESignals sig, SigHandler_t h);
   static const char  *UnixSigname(ESignals sig);
   static void         UnixSigAlarmInterruptsSyscalls(Bool_t set);
   static void         UnixResetSignal(ESignals sig);
   static void         UnixResetSignals();
   static void         UnixIgnoreSignal(ESignals sig, Bool_t ignore);
   static int          UnixFilestat(const char *path, FileStat_t &buf);
   static int          UnixFSstat(const char *path, Long_t *id, Long_t *bsize,
                                  Long_t *blocks, Long_t *bfree);
   static int          UnixTcpConnect(const char *hostname, int port, int tcpwindowsize);
   static int          UnixUdpConnect(const char *hostname, int port);
   static int          UnixUnixConnect(int port);
   static int          UnixUnixConnect(const char *path);
   static int          UnixTcpService(int port, Bool_t reuse, int backlog,
                                      int tcpwindowsize);
   static int          UnixUdpService(int port, int backlog);
   static int          UnixUnixService(int port, int backlog);
   static int          UnixUnixService(const char *sockpath, int backlog);
   static int          UnixRecv(int sock, void *buf, int len, int flag);
   static int          UnixSend(int sock, const void *buf, int len, int flag);

public:
   TUnixSystem();
   virtual ~TUnixSystem();

   //---- Misc -------------------------------------------------
   Bool_t            Init() override;
   void              SetProgname(const char *name) override;
   void              SetDisplay() override;
   const char       *GetError() override;
   const char       *HostName() override;

   //---- EventLoop --------------------------------------------
   void              DispatchOneEvent(Bool_t pendingOnly = kFALSE) override;
   Int_t             Select(TList *active, Long_t timeout) override;
   Int_t             Select(TFileHandler *fh, Long_t timeout) override;

   //---- Handling of system events ----------------------------
   void              CheckChilds();
   Bool_t            CheckSignals(Bool_t sync);
   Bool_t            CheckDescriptors();
   void              DispatchSignals(ESignals sig);
   void              AddSignalHandler(TSignalHandler *sh) override;
   TSignalHandler   *RemoveSignalHandler(TSignalHandler *sh) override;
   void              ResetSignal(ESignals sig, Bool_t reset = kTRUE) override;
   void              ResetSignals() override;
   void              IgnoreSignal(ESignals sig, Bool_t ignore = kTRUE) override;
   void              SigAlarmInterruptsSyscalls(Bool_t set) override;
   void              AddFileHandler(TFileHandler *fh) override;
   TFileHandler     *RemoveFileHandler(TFileHandler *fh) override;

   //---- Floating Point Exceptions Control --------------------
   Int_t             GetFPEMask() override;
   Int_t             SetFPEMask(Int_t mask = kDefaultMask) override;

   //---- Time & Date ------------------------------------------
   TTime             Now() override;
   void              AddTimer(TTimer *ti) override;
   TTimer           *RemoveTimer(TTimer *ti) override;
   void              ResetTimer(TTimer *ti) override;
   Bool_t            DispatchTimers(Bool_t mode);
   void              Sleep(UInt_t milliSec) override;

   //---- Processes --------------------------------------------
   Int_t             Exec(const char *shellcmd) override;
   FILE             *OpenPipe(const char *shellcmd, const char *mode) override;
   int               ClosePipe(FILE *pipe) override;
   void              Exit(int code, Bool_t mode = kTRUE) override;
   void              Abort(int code = 0) override;
   int               GetPid() override;
   void              StackTrace() override;

   //---- Directories ------------------------------------------
   int               MakeDirectory(const char *name) override;
   void             *OpenDirectory(const char *name) override;
   void              FreeDirectory(void *dirp) override;
   const char       *GetDirEntry(void *dirp) override;
   Bool_t            ChangeDirectory(const char *path) override;
   const char       *WorkingDirectory() override;
   std::string       GetWorkingDirectory() const override;
   const char       *HomeDirectory(const char *userName = nullptr) override;
   std::string       GetHomeDirectory(const char *userName = nullptr) const override;
   const char       *TempDirectory() const override;
   FILE             *TempFileName(TString &base, const char *dir = nullptr) override;

   //---- Paths & Files ----------------------------------------
   const char       *PrependPathName(const char *dir, TString& name) override;
   Bool_t            ExpandPathName(TString &patbuf) override;
   char             *ExpandPathName(const char *path) override;
   Bool_t            AccessPathName(const char *path, EAccessMode mode = kFileExists) override;
   Bool_t            IsPathLocal(const char *path) override;
   int               CopyFile(const char *from, const char *to, Bool_t overwrite = kFALSE) override;
   int               Rename(const char *from, const char *to) override;
   int               Link(const char *from, const char *to) override;
   int               Symlink(const char *from, const char *to) override;
   int               Unlink(const char *name) override;
   int               GetPathInfo(const char *path, FileStat_t &buf) override;
   int               GetFsInfo(const char *path, Long_t *id, Long_t *bsize,
                               Long_t *blocks, Long_t *bfree) override;
   int               Chmod(const char *file, UInt_t mode) override;
   int               Umask(Int_t mask) override;
   int               Utime(const char *file, Long_t modtime, Long_t actime) override;
   const char       *FindFile(const char *search, TString& file, EAccessMode mode = kFileExists) override;

   //---- Users & Groups ---------------------------------------
   Int_t             GetUid(const char *user = nullptr) override;
   Int_t             GetGid(const char *group = nullptr) override;
   Int_t             GetEffectiveUid() override;
   Int_t             GetEffectiveGid() override;
   UserGroup_t      *GetUserInfo(Int_t uid) override;
   UserGroup_t      *GetUserInfo(const char *user = nullptr) override;
   UserGroup_t      *GetGroupInfo(Int_t gid) override;
   UserGroup_t      *GetGroupInfo(const char *group = nullptr) override;

   //---- Environment Manipulation -----------------------------
   const char       *Getenv(const char *name) override;
   void              Setenv(const char *name, const char *value) override;
   void              Unsetenv(const char *name) override;

   //---- System Logging ---------------------------------------
   void              Openlog(const char *name, Int_t options, ELogFacility facility) override;
   void              Syslog(ELogLevel level, const char *mess) override;
   void              Closelog() override;

   //---- Standard Output redirection --------------------------
   Int_t             RedirectOutput(const char *name, const char *mode = "a",
                                    RedirectHandle_t *h = nullptr) override;

   //---- Dynamic Loading --------------------------------------
   void              AddDynamicPath(const char *lib) override;
   const char       *GetDynamicPath() override;
   void              SetDynamicPath(const char *lib) override;
   const char       *FindDynamicLibrary(TString &lib, Bool_t quiet = kFALSE) override;
   Func_t            DynFindSymbol(const char *module, const char *entry) override;
   int               Load(const char *module, const char *entry = "", Bool_t system = kFALSE) override;
   void              Unload(const char *module) override;
   void              ListSymbols(const char *module, const char *re = "") override;
   void              ListLibraries(const char *regexp = "") override;

   //---- RPC --------------------------------------------------
   TInetAddress      GetHostByName(const char *server) override;
   TInetAddress      GetPeerName(int sock) override;
   TInetAddress      GetSockName(int sock) override;
   int               GetServiceByName(const char *service) override;
   char             *GetServiceByPort(int port) override;
   int               ConnectService(const char *server, int port, int tcpwindowsize, const char *protocol = "tcp");
   int               OpenConnection(const char *server, int port, int tcpwindowsize = -1, const char *protocol = "tcp") override;
   int               AnnounceTcpService(int port, Bool_t reuse, int backlog, int tcpwindowsize = -1) override;
   int               AnnounceUdpService(int port, int backlog) override;
   int               AnnounceUnixService(int port, int backlog) override;
   int               AnnounceUnixService(const char *sockpath, int backlog) override;
   int               AcceptConnection(int sock) override;
   void              CloseConnection(int sock, Bool_t force = kFALSE) override;
   int               RecvRaw(int sock, void *buffer, int length, int flag) override;
   int               SendRaw(int sock, const void *buffer, int length, int flag) override;
   int               RecvBuf(int sock, void *buffer, int length) override;
   int               SendBuf(int sock, const void *buffer, int length) override;
   int               SetSockOpt(int sock, int option, int val) override;
   int               GetSockOpt(int sock, int option, int *val) override;

   //---- System, CPU and Memory info
   int               GetSysInfo(SysInfo_t *info) const override;
   int               GetCpuInfo(CpuInfo_t *info, Int_t sampleTime = 1000) const override;
   int               GetMemInfo(MemInfo_t *info) const override;
   int               GetProcInfo(ProcInfo_t *info) const override;

   ClassDefOverride(TUnixSystem,0)  //Interface to Unix OS services
};

#endif
