// @(#)root/winnt:$Id$
// Author: Fons Rademakers   15/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef ROOT_TWinNTSystem
#define ROOT_TWinNTSystem

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TWinNTSystem                                                         //
//                                                                      //
// Class providing an interface to the Windows NT Operating System.     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TSystem.h"
#include <string>

#if !defined(__CINT__)
 #include "Windows4Root.h"
 #include <commctrl.h>
#else
 typedef void* HANDLE;
 struct WIN32_FIND_DATA;
 typedef void* HIMAGELIST;
 typedef void* HICON;
 typedef UChar_t BOOL;
 struct FILE;
#endif

#ifndef MAX_SID_SIZE
#define MAX_SID_SIZE   1024
#endif
#define MAX_NAME_STRING   1024

#define SID_GROUP         0
#define SID_MEMBER        1

struct passwd {
   char *pw_name;       // user name
   char *pw_passwd;     // user password
   int   pw_uid;        // user ID
   int   pw_gid;        // user's group ID
   int   pw_quota;      //
   char *pw_gecos;      // user's real (long) name
   char *pw_dir;        // user home directory
   char *pw_shell;      // shell command
   char *pw_group;      // user's group name
};

struct group {
   char   *gr_name;     // group name
   char   *gr_passwd;   // group password
   int    gr_gid;       // group id
   char   **gr_mem;     // group members
};


class TWinNTSystem : public TSystem {
public:
   // pointer to message handler func
   typedef Bool_t (*ThreadMsgFunc_t)(MSG*);

private:
   struct group     *fGroups{nullptr};           // Groups on local computer
   struct passwd    *fPasswords{nullptr};        // Users on local computer
   int               fNbUsers{0};                // Number of users on local computer
   int               fNbGroups{0};               // Number of groups on local computer
   int               fActUser{-1};               // Index of actual user in User list
   Bool_t            fGroupsInitDone{kFALSE};    // Flag used for Users and Groups initialization
   Bool_t            fFirstFile{kFALSE};         // Flag used by OpenDirectory/GetDirEntry

   HANDLE            fhProcess;                  // Handle of the current process
   void             *fGUIThreadHandle{nullptr};  // handle of GUI server (aka command) thread
   ULong_t           fGUIThreadId{0};            // id of GUI server (aka command) thread
   std::string       fDirNameBuffer;             // The string buffer to hold path name
   WIN32_FIND_DATA   fFindFileData;              // Structure to look for files (aka OpenDir under UNIX)

   Bool_t            DispatchTimers(Bool_t mode);
   Bool_t            CheckDescriptors();
   Bool_t            CheckSignals(Bool_t sync);
   Bool_t            CountMembers(const char *lpszGroupName);
   const char       *GetLinkedLibraries() override;
   Bool_t            GetNbGroups();
   Long_t            LookupSID (const char *lpszAccountName, int what, int &groupIdx, int &memberIdx);
   Bool_t            CollectMembers(const char *lpszGroupName, int &groupIdx, int &memberIdx);
   Bool_t            CollectGroups();
   Bool_t            InitUsersGroups();
   void              DoBeep(Int_t freq=-1, Int_t duration=-1) const override;

   static void       ThreadStub(void *Parameter) {((TWinNTSystem *)Parameter)->TimerThread();}
   void              TimerThread();
   void              FillWithHomeDirectory(const char *userName, char *mydir) const;
   char             *GetWorkingDirectory(char driveletter) const;


protected:
   static int        WinNTUnixConnect(int port);
   static int        WinNTUnixConnect(const char *path);
   static int        WinNTUdpConnect(const char *hostname, int port);

public:
   TWinNTSystem();
   virtual ~TWinNTSystem();

   //---- non-TSystem methods ----------------------------------
   HANDLE            GetProcess();
   Bool_t            HandleConsoleEvent();

   //---- Misc -------------------------------------------------
   Bool_t            Init() override;
   const char       *BaseName(const char *name) override;
   void              SetProgname(const char *name) override;
   const char       *GetError() override;
   const char       *HostName() override;
   void             *GetGUIThreadHandle() const {return fGUIThreadHandle;}
   ULong_t           GetGUIThreadId() const {return fGUIThreadId;}
   void              SetGUIThreadMsgHandler(ThreadMsgFunc_t func);
   void              NotifyApplicationCreated() override;


   //---- EventLoop --------------------------------------------
   Bool_t            ProcessEvents() override;
   void              DispatchOneEvent(Bool_t pendingOnly = kFALSE) override;
   void              ExitLoop() override;
   Int_t             Select(TList *active, Long_t timeout) override;
   Int_t             Select(TFileHandler *fh, Long_t timeout) override;

   //---- Handling of system events ----------------------------
   void              DispatchSignals(ESignals sig);
   void              AddSignalHandler(TSignalHandler *sh) override;
   TSignalHandler   *RemoveSignalHandler(TSignalHandler *sh) override;
   void              ResetSignal(ESignals sig, Bool_t reset = kTRUE) override;
   void              ResetSignals() override;
   void              IgnoreSignal(ESignals sig, Bool_t ignore = kTRUE) override;
   void              AddFileHandler(TFileHandler *fh) override;
   TFileHandler     *RemoveFileHandler(TFileHandler *fh) override;
   void              StackTrace() override;

   //---- Floating Point Exceptions Control --------------------
   Int_t             GetFPEMask() override;
   Int_t             SetFPEMask(Int_t mask = kDefaultMask) override;

   //---- Processes --------------------------------------------
   int               Exec(const char *shellcmd) override;
   FILE             *OpenPipe(const char *shellcmd, const char *mode) override;
   int               ClosePipe(FILE *pipe) override;
   void              Exit(int code, Bool_t mode = kTRUE) override;
   void              Abort(int code = 0) override;
   int               GetPid() override;

   //---- Environment manipulation -----------------------------
   const char       *Getenv(const char *name) override;
   void              Setenv(const char *name, const char *value) override;

   //---- Directories ------------------------------------------
   int               mkdir(const char *name, Bool_t recursive = kFALSE) override;
   int               MakeDirectory(const char *name) override;
   Bool_t            ChangeDirectory(const char *path) override;
   const char       *GetDirEntry(void *dirp) override;
   const char       *DirName(const char *pathname) override;
   TString           GetDirName(const char *pathname) override;
   void              FreeDirectory(void *dirp) override;
   void             *OpenDirectory(const char *name) override;
   const char       *WorkingDirectory(char driveletter);
   const char       *WorkingDirectory() override;
   std::string       GetWorkingDirectory() const override;
   const char       *HomeDirectory(const char *userName=0) override;
   std::string       GetHomeDirectory(const char *userName = nullptr) const override;
   const char       *TempDirectory() const override;
   FILE             *TempFileName(TString &base, const char *dir = nullptr) override;

   //---- Users & Groups ---------------------------------------
   Int_t             GetUid(const char *user = nullptr) override;
   Int_t             GetGid(const char *group = nullptr) override;
   Int_t             GetEffectiveUid() override;
   Int_t             GetEffectiveGid() override;
   UserGroup_t      *GetUserInfo(Int_t uid) override;
   UserGroup_t      *GetUserInfo(const char *user = nullptr) override;
   UserGroup_t      *GetGroupInfo(Int_t gid) override;
   UserGroup_t      *GetGroupInfo(const char *group = nullptr) override;

   //---- Paths & Files ----------------------------------------
   const char        DriveName(const char *pathname="/");
   const char       *PrependPathName(const char *dir, TString& name) override;
   Bool_t            ExpandPathName(TString &patbuf) override;
   char             *ExpandPathName(const char *path) override;
   Bool_t            AccessPathName(const char *path, EAccessMode mode = kFileExists) override;
   Bool_t            IsPathLocal(const char *path) override;
   Bool_t            IsAbsoluteFileName(const char *dir) override;
   int               CopyFile(const char *from, const char *to, Bool_t overwrite = kFALSE) override;
   int               Rename(const char *from, const char *to) override;
   int               Link(const char *from, const char *to) override;
   int               Symlink(const char *from, const char *to) override;
   int               Unlink(const char *name) override;
   int               SetNonBlock(int fd);
   int               GetPathInfo(const char *path, FileStat_t &buf) override;
   int               GetFsInfo(const char *path, Long_t *id, Long_t *bsize,
                                 Long_t *blocks, Long_t *bfree) override;
   int               Chmod(const char *file, UInt_t mode) override;
   int               Umask(Int_t mask) override;
   int               Utime(const char *file, Long_t modtime, Long_t actime) override;
   const char       *UnixPathName(const char *unixpathname) override;
   const char       *FindFile(const char *search, TString& file, EAccessMode mode = kFileExists) override;
   TList            *GetVolumes(Option_t *opt = "") const override;

   //---- Standard Output redirection --------------------------
   Int_t             RedirectOutput(const char *name, const char *mode = "a", RedirectHandle_t *h = nullptr) override;

   //---- Dynamic Loading --------------------------------------
   void              AddDynamicPath(const char *dir) override;
   const char       *GetDynamicPath() override;
   void              SetDynamicPath(const char *path) override;
   const char       *FindDynamicLibrary(TString &lib, Bool_t quiet = kFALSE) override;
   int               Load(const char *module, const char *entry = "", Bool_t system = kFALSE) override;
   const char       *GetLibraries(const char *regexp = "",
                                  const char *option = "",
                                  Bool_t isRegexp = kTRUE) override;

   //---- Time & Date -------------------------------------------
   TTime             Now() override;
   void              AddTimer(TTimer *ti) override;
   TTimer           *RemoveTimer(TTimer *ti) override;
   void              Sleep(UInt_t milliSec) override;
   Double_t          GetRealTime();
   Double_t          GetCPUTime();

   //---- RPC --------------------------------------------------
   int               ConnectService(const char *servername, int port, int tcpwindowsize, const char *protocol = "tcp");
   TInetAddress      GetHostByName(const char *server) override;
   TInetAddress      GetPeerName(int sock) override;
   TInetAddress      GetSockName(int sock) override;
   int               GetServiceByName(const char *service) override;
   char              *GetServiceByPort(int port) override;
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
   int               SetSockOpt(int sock, int opt, int val) override;
   int               GetSockOpt(int sock, int opt, int *val) override;

   //---- System, CPU and Memory info
   Int_t             GetSysInfo(SysInfo_t *info) const override;
   Int_t             GetCpuInfo(CpuInfo_t *info, Int_t sampleTime = 1000) const override;
   Int_t             GetMemInfo(MemInfo_t *info) const override;
   Int_t             GetProcInfo(ProcInfo_t *info) const override;

   ClassDefOverride(TWinNTSystem, 0)
};

R__EXTERN ULong_t gConsoleWindow;   // console window handle

#endif
