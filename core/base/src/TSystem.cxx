// @(#)root/base:$Id: 8944840ba34631ec28efc779647618db43c0eee5 $
// Author: Fons Rademakers   15/09/95

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TSystem
\ingroup Base

Abstract base class defining a generic interface to the underlying
Operating System.
This is not an ABC in the strict sense of the (C++) word. For
every member function there is an implementation (often not more
than a call to AbstractMethod() which prints a warning saying
that the method should be overridden in a derived class), which
allows a simple partial implementation for new OS'es.
*/

#include <ROOT/FoundationUtils.hxx>
#include "strlcpy.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TException.h"
#include "TROOT.h"
#include "TClass.h"
#include "TClassTable.h"
#include "TEnv.h"
#include "TOrdCollection.h"
#include "TObject.h"
#include "TInterpreter.h"
#include "TRegexp.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TError.h"
#include "TPluginManager.h"
#include "TUrl.h"
#include "TVirtualMutex.h"
#include "TVersionCheck.h"
#include "compiledata.h"
#include "RConfigure.h"
#include "THashList.h"
#include "ThreadLocalStorage.h"

#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <set>

#ifdef WIN32
#include <io.h>
#endif

const char *gRootDir = nullptr;
const char *gProgName = nullptr;
const char *gProgPath = nullptr;

TSystem      *gSystem   = nullptr;
TFileHandler *gXDisplay = nullptr;  // Display server event handler, set in TGClient

static Int_t *gLibraryVersion    = nullptr;   // Set in TVersionCheck, used in Load()
static Int_t  gLibraryVersionIdx = 0;   // Set in TVersionCheck, used in Load()
static Int_t  gLibraryVersionMax = 256;

// Pin vtable
ProcInfo_t::~ProcInfo_t() {}
ClassImp(TProcessEventTimer);

////////////////////////////////////////////////////////////////////////////////
/// Create async event processor timer. Delay is in milliseconds.

TProcessEventTimer::TProcessEventTimer(Long_t delay) : TTimer(delay, kFALSE)
{
   gROOT->SetInterrupt(kFALSE);
   TurnOn();
}

////////////////////////////////////////////////////////////////////////////////
/// Process events if timer did time out. Returns kTRUE if interrupt
/// flag is set (by hitting a key in the canvas or selecting the
/// Interrupt menu item in canvas or some other action).

Bool_t TProcessEventTimer::ProcessEvents()
{
   if (fTimeout) {
      if (gSystem->ProcessEvents()) {
         Remove();
         return kTRUE;
      } else {
         Reset();
         return kFALSE;
      }
   }
   return kFALSE;
}



ClassImp(TSystem);

TVirtualMutex* gSystemMutex = nullptr;



////////////////////////////////////////////////////////////////////////////////
/// Strip off protocol string from specified path

const char *TSystem::StripOffProto(const char *path, const char *proto)
{
   return !strncmp(path, proto, strlen(proto)) ? path + strlen(proto) : path;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new OS interface.

TSystem::TSystem(const char *name, const char *title) : TNamed(name, title)
{
   if (gSystem && name[0] != '-' && strcmp(name, "Generic"))
      Error("TSystem", "only one instance of TSystem allowed");

   if (!gLibraryVersion) {
      gLibraryVersion = new Int_t [gLibraryVersionMax];
      memset(gLibraryVersion, 0, gLibraryVersionMax*sizeof(Int_t));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete the OS interface.

TSystem::~TSystem()
{
   if (fOnExitList) {
      fOnExitList->Delete();
      SafeDelete(fOnExitList);
   }

   if (fSignalHandler) {
      fSignalHandler->Delete();
      SafeDelete(fSignalHandler);
   }

   if (fFileHandler) {
      fFileHandler->Delete();
      SafeDelete(fFileHandler);
   }

   if (fStdExceptionHandler) {
      fStdExceptionHandler->Delete();
      SafeDelete(fStdExceptionHandler);
   }

   if (fTimers) {
      fTimers->Delete();
      SafeDelete(fTimers);
   }

   if (fCompiled) {
      fCompiled->Delete();
      SafeDelete(fCompiled);
   }

   if (fHelpers) {
      fHelpers->Delete();
      SafeDelete(fHelpers);
   }

   if (gSystem == this)
      gSystem = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the OS interface.

Bool_t TSystem::Init()
{
   fNfd    = 0;
   fMaxrfd = -1;
   fMaxwfd = -1;

   fSigcnt = 0;
   fLevel  = 0;

   fSignalHandler       = new TOrdCollection;
   fFileHandler         = new TOrdCollection;
   fStdExceptionHandler = new TOrdCollection;
   fTimers              = new TOrdCollection;

   fBuildArch     = BUILD_ARCH;
   fBuildCompiler = COMPILER;
   fBuildCompilerVersion = COMPILERVERS;
   fBuildCompilerVersionStr = COMPILERVERSSTR;
   fBuildNode     = BUILD_NODE;
   fFlagsDebug    = CXXDEBUG;
   fFlagsOpt      = CXXOPT;
   fIncludePath   = INCLUDEPATH;
   fLinkedLibs    = LINKEDLIBS;
   fSoExt         = SOEXT;
   fObjExt        = OBJEXT;
   fAclicMode     = kDefault;
   fMakeSharedLib = MAKESHAREDLIB;
   fMakeExe       = MAKEEXE;
   fCompiled      = new TOrdCollection;

   if (gEnv && fBeepDuration == 0 && fBeepFreq == 0) {
      fBeepDuration = gEnv->GetValue("Root.System.BeepDuration", 100);
      fBeepFreq     = gEnv->GetValue("Root.System.BeepFreq", 440);
   }
   if (!fName.CompareTo("Generic")) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the application name (from command line, argv[0]) and copy it in
/// gProgName.

void TSystem::SetProgname(const char *name)
{
   delete [] gProgName;
   gProgName = StrDup(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Set DISPLAY environment variable based on utmp entry. Only for UNIX.

void TSystem::SetDisplay()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Set the system error string. This string will be used by GetError().
/// To be used in case one does not want or can use the system error
/// string (e.g. because error is generated by a third party POSIX like
/// library that does not use standard errno).

void TSystem::SetErrorStr(const char *errstr)
{
   ResetErrno();   // so GetError() uses the fLastErrorString
   GetLastErrorString() = errstr;
}

////////////////////////////////////////////////////////////////////////////////
/// Return system error string.

const char *TSystem::GetError()
{
   if (GetErrno() == 0 && !GetLastErrorString().IsNull())
      return GetLastErrorString().Data();
   return Form("errno: %d", GetErrno());
}

////////////////////////////////////////////////////////////////////////////////
/// Static function returning system error number.

Int_t TSystem::GetErrno()
{
#ifdef _REENTRANT
   return errno; // errno can be a macro if _REENTRANT is set
#else
#ifdef R__SOLARIS_CC50
   return ::errno;
#else
   return errno;
#endif
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Static function resetting system error number.

void TSystem::ResetErrno()
{
#ifdef _REENTRANT
   errno = 0; // errno can be a macro if _REENTRANT is set
#else
#ifdef R__SOLARIS_CC50
   ::errno = 0;
#else
   errno = 0;
#endif
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Objects that should be deleted on exit of the OS interface.

void TSystem::RemoveOnExit(TObject *obj)
{
   if (!fOnExitList)
      fOnExitList = new TOrdCollection;
   if (!fOnExitList->FindObject(obj))
      fOnExitList->Add(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the system's host name.

const char *TSystem::HostName()
{
   return "Local host";
}

////////////////////////////////////////////////////////////////////////////////
/// Hook to tell TSystem that the TApplication object has been created.

void TSystem::NotifyApplicationCreated()
{
   // Currently needed only for WinNT interface.
}

////////////////////////////////////////////////////////////////////////////////
/// Beep for duration milliseconds with a tone of frequency freq.
/// Defaults to printing the `\a` character to stdout.
/// If freq or duration is <0 respectively, use default value.
/// If setDefault is set, only set the frequency and duration as
/// new defaults, but don't beep.
/// If default freq or duration is <0, never beep (silence)

void TSystem::Beep(Int_t freq /*=-1*/, Int_t duration /*=-1*/,
                   Bool_t setDefault /*=kFALSE*/)
{
   if (setDefault) {
      fBeepFreq     = freq;
      fBeepDuration = duration;
      return;
   }
   if (fBeepDuration < 0 || fBeepFreq < 0) return; // silence
   if (freq < 0) freq = fBeepFreq;
   if (duration < 0) duration = fBeepDuration;
   DoBeep(freq, duration);
}

//---- EventLoop ---------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// System event loop.

void TSystem::Run()
{
   fInControl = kTRUE;
   fDone      = kFALSE;

loop_entry:
   try {
      RETRY {
         while (!fDone) {
            gApplication->StartIdleing();
            InnerLoop();
            gApplication->StopIdleing();
         }
      } ENDTRY;
   }
   catch (std::exception& exc) {
      TIter next(fStdExceptionHandler);
      TStdExceptionHandler* eh = 0;
      while ((eh = (TStdExceptionHandler*) next())) {
         switch (eh->Handle(exc))
         {
            case TStdExceptionHandler::kSEProceed:
               break;
            case TStdExceptionHandler::kSEHandled:
               goto loop_entry;
               break;
            case TStdExceptionHandler::kSEAbort:
               Warning("Run", "instructed to abort");
               goto loop_end;
               break;
         }
      }
      throw;
   }
   catch (const char *str) {
      printf("%s\n", str);
   }
   // handle every exception
   catch (...) {
      Warning("Run", "handle uncaugth exception, terminating");
   }

loop_end:
   fInControl = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Exit from event loop.

void TSystem::ExitLoop()
{
   fDone = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Inner event loop.

void TSystem::InnerLoop()
{
   fLevel++;
   DispatchOneEvent();
   fLevel--;
}

////////////////////////////////////////////////////////////////////////////////
/// Process pending events (GUI, timers, sockets). Returns the result of
/// TROOT::IsInterrupted(). The interrupt flag (TROOT::SetInterrupt())
/// can be set during the handling of the events. This mechanism allows
/// macros running in tight calculating loops to be interrupted by some
/// GUI event (depending on the interval with which this method is
/// called). For example hitting ctrl-c in a canvas will set the
/// interrupt flag.

Bool_t TSystem::ProcessEvents()
{
   gROOT->SetInterrupt(kFALSE);

   if (!gROOT->TestBit(TObject::kInvalidObject))
      DispatchOneEvent(kTRUE);

   return gROOT->IsInterrupted();
}

////////////////////////////////////////////////////////////////////////////////
/// Dispatch a single event.

void TSystem::DispatchOneEvent(Bool_t)
{
   AbstractMethod("DispatchOneEvent");
}

////////////////////////////////////////////////////////////////////////////////
/// Sleep milliSec milli seconds.

void TSystem::Sleep(UInt_t)
{
   AbstractMethod("Sleep");
}

////////////////////////////////////////////////////////////////////////////////
/// Select on active file descriptors (called by TMonitor).

Int_t TSystem::Select(TList *, Long_t)
{
   AbstractMethod("Select");
   return -1;
}
////////////////////////////////////////////////////////////////////////////////
/// Select on active file descriptors (called by TMonitor).

Int_t TSystem::Select(TFileHandler *, Long_t)
{
   AbstractMethod("Select");
   return -1;
}

//---- handling of system events -----------------------------------------------
////////////////////////////////////////////////////////////////////////////////
/// Get current time in milliseconds since 0:00 Jan 1 1995.

TTime TSystem::Now()
{
   return TTime(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Add timer to list of system timers.

void TSystem::AddTimer(TTimer *ti)
{
   if (ti && fTimers && (fTimers->FindObject(ti) == nullptr))
      fTimers->Add(ti);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove timer from list of system timers. Returns removed timer or 0
/// if timer was not active.

TTimer *TSystem::RemoveTimer(TTimer *ti)
{
   if (fTimers) {
      TTimer *tr = (TTimer*) fTimers->Remove(ti);
      return tr;
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Time when next timer of mode (synchronous=kTRUE or
/// asynchronous=kFALSE) will time-out (in ms).

Long_t TSystem::NextTimeOut(Bool_t mode)
{
   if (!fTimers) return -1;

   TOrdCollectionIter it((TOrdCollection*)fTimers);
   TTimer *t, *to = nullptr;
   Long64_t tt, tnow = Now();
   Long_t   timeout = -1;

   while ((t = (TTimer *) it.Next())) {
      if (t->IsSync() == mode) {
         tt = (Long64_t)t->GetAbsTime() - tnow;
         if (tt < 0) tt = 0;
         if (timeout == -1) {
            timeout = (Long_t)tt;
            to = t;
         }
         if (tt < timeout) {
            timeout = (Long_t)tt;
            to = t;
         }
      }
   }

   if (to && to->IsAsync() && timeout > 0) {
      if (to->IsInterruptingSyscalls())
         SigAlarmInterruptsSyscalls(kTRUE);
      else
         SigAlarmInterruptsSyscalls(kFALSE);
   }

   return timeout;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a signal handler to list of system signal handlers. Only adds
/// the handler if it is not already in the list of signal handlers.

void TSystem::AddSignalHandler(TSignalHandler *h)
{
   if (h && fSignalHandler && (fSignalHandler->FindObject(h) == nullptr))
      fSignalHandler->Add(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a signal handler from list of signal handlers. Returns
/// the handler or 0 if the handler was not in the list of signal handlers.

TSignalHandler *TSystem::RemoveSignalHandler(TSignalHandler *h)
{
   if (fSignalHandler)
      return (TSignalHandler *)fSignalHandler->Remove(h);

   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a file handler to the list of system file handlers. Only adds
/// the handler if it is not already in the list of file handlers.

void TSystem::AddFileHandler(TFileHandler *h)
{
   if (h && fFileHandler && !fFileHandler->FindObject(h))
      fFileHandler->Add(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a file handler from the list of file handlers. Returns
/// the handler or 0 if the handler was not in the list of file handlers.

TFileHandler *TSystem::RemoveFileHandler(TFileHandler *h)
{
   if (fFileHandler)
      return (TFileHandler *)fFileHandler->Remove(h);

   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// If reset is true reset the signal handler for the specified signal
/// to the default handler, else restore previous behaviour.

void TSystem::ResetSignal(ESignals /*sig*/, Bool_t /*reset*/)
{
   AbstractMethod("ResetSignal");
}

////////////////////////////////////////////////////////////////////////////////
/// Reset signals handlers to previous behaviour.

void TSystem::ResetSignals()
{
   AbstractMethod("ResetSignals");
}

////////////////////////////////////////////////////////////////////////////////
/// If ignore is true ignore the specified signal, else restore previous
/// behaviour.

void TSystem::IgnoreSignal(ESignals /*sig*/, Bool_t /*ignore*/)
{
   AbstractMethod("IgnoreSignal");
}

////////////////////////////////////////////////////////////////////////////////
/// If ignore is true ignore the interrupt signal, else restore previous
/// behaviour. Typically call ignore interrupt before writing to disk.

void TSystem::IgnoreInterrupt(Bool_t ignore)
{
   IgnoreSignal(kSigInterrupt, ignore);
}

////////////////////////////////////////////////////////////////////////////////
/// Add an exception handler to list of system exception handlers. Only adds
/// the handler if it is not already in the list of exception handlers.

void TSystem::AddStdExceptionHandler(TStdExceptionHandler *eh)
{
   if (eh && fStdExceptionHandler && !fStdExceptionHandler->FindObject(eh))
      fStdExceptionHandler->Add(eh);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove an exception handler from list of exception handlers. Returns
/// the handler or 0 if the handler was not in the list of exception handlers.

TStdExceptionHandler *TSystem::RemoveStdExceptionHandler(TStdExceptionHandler *eh)
{
   if (fStdExceptionHandler)
      return (TStdExceptionHandler *)fStdExceptionHandler->Remove(eh);

   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the bitmap of conditions that trigger a floating point exception.

Int_t TSystem::GetFPEMask()
{
   AbstractMethod("GetFPEMask");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set which conditions trigger a floating point exception.
/// Return the previous set of conditions.

Int_t TSystem::SetFPEMask(Int_t)
{
   AbstractMethod("SetFPEMask");
   return 0;
}

//---- Processes ---------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Execute a command.

int TSystem::Exec(const char *)
{
   AbstractMethod("Exec");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a pipe.

FILE *TSystem::OpenPipe(const char *, const char *)
{
   AbstractMethod("OpenPipe");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Close the pipe.

int TSystem::ClosePipe(FILE*)
{
   AbstractMethod("ClosePipe");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute command and return output in TString.

TString TSystem::GetFromPipe(const char *command)
{
   TString out;

   FILE *pipe = OpenPipe(command, "r");
   if (!pipe) {
      SysError("GetFromPipe", "cannot run command \"%s\"", command);
      return out;
   }

   TString line;
   while (line.Gets(pipe)) {
      if (out != "")
         out += "\n";
      out += line;
   }

   Int_t r = ClosePipe(pipe);
   if (r) {
      Error("GetFromPipe", "command \"%s\" returned %d", command, r);
   }
   return out;
}

////////////////////////////////////////////////////////////////////////////////
/// Get process id.

int TSystem::GetPid()
{
   AbstractMethod("GetPid");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Exit the application.

void TSystem::Exit(int, Bool_t)
{
   AbstractMethod("Exit");
}

////////////////////////////////////////////////////////////////////////////////
/// Abort the application.

void TSystem::Abort(int)
{
   AbstractMethod("Abort");
}

////////////////////////////////////////////////////////////////////////////////
/// Print a stack trace.

void TSystem::StackTrace()
{
   AbstractMethod("StackTrace");
}


//---- Directories -------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Create helper TSystem to handle file and directory operations that
/// might be special for remote file access.

TSystem *TSystem::FindHelper(const char *path, void *dirptr)
{
   TSystem *helper = nullptr;
   {
      R__READ_LOCKGUARD(ROOT::gCoreMutex);

      if (!fHelpers) {
         R__WRITE_LOCKGUARD(ROOT::gCoreMutex);
         fHelpers = new TOrdCollection;
      }

      if (path) {
         if (!GetDirPtr()) {
            TUrl url(path, kTRUE);
            if (!strcmp(url.GetProtocol(), "file"))
               return nullptr;
         }
      }

      // look for existing helpers
      TIter next(fHelpers);
      while ((helper = (TSystem*) next()))
         if (helper->ConsistentWith(path, dirptr))
            return helper;

      if (!path)
         return nullptr;
   }

   // create new helper
   TRegexp re("^root.*:");  // also roots, rootk, etc
   TString pname = path;
   TPluginHandler *h;
   if (pname.BeginsWith("xroot:") || pname.Index(re) != kNPOS) {
      // (x)rootd daemon ...
      if ((h = gROOT->GetPluginManager()->FindHandler("TSystem", path))) {
         if (h->LoadPlugin() == -1)
            return nullptr;
         helper = (TSystem*) h->ExecPlugin(2, path, kFALSE);
      }
   } else if ((h = gROOT->GetPluginManager()->FindHandler("TSystem", path))) {
      if (h->LoadPlugin() == -1)
         return nullptr;
      helper = (TSystem*) h->ExecPlugin(0);
   }

   if (helper) {
      R__WRITE_LOCKGUARD(ROOT::gCoreMutex);
      fHelpers->Add(helper);
   }

   return helper;
}

////////////////////////////////////////////////////////////////////////////////
/// Check consistency of this helper with the one required
/// by 'path' or 'dirptr'

Bool_t TSystem::ConsistentWith(const char *path, void *dirptr)
{
   Bool_t checkproto = kFALSE;
   if (path) {
      if (!GetDirPtr()) {
         TUrl url(path, kTRUE);
         if (!strncmp(url.GetProtocol(), GetName(), strlen(GetName())))
            checkproto = kTRUE;
      }
   }

   Bool_t checkdir = kFALSE;
   if (GetDirPtr() && GetDirPtr() == dirptr)
      checkdir = kTRUE;

   return (checkproto || checkdir);
}

////////////////////////////////////////////////////////////////////////////////
/// Make a directory. Returns 0 in case of success and
/// -1 if the directory could not be created (either already exists or
/// illegal path name).

int TSystem::MakeDirectory(const char *)
{
   AbstractMethod("MakeDirectory");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a directory. Returns 0 if directory does not exist.

void *TSystem::OpenDirectory(const char *)
{
   AbstractMethod("OpenDirectory");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Free a directory.

void TSystem::FreeDirectory(void *)
{
   AbstractMethod("FreeDirectory");
}

////////////////////////////////////////////////////////////////////////////////
/// Get a directory entry. Returns 0 if no more entries.

const char *TSystem::GetDirEntry(void *)
{
   AbstractMethod("GetDirEntry");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Change directory.

Bool_t TSystem::ChangeDirectory(const char *)
{
   AbstractMethod("ChangeDirectory");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return working directory.

const char *TSystem::WorkingDirectory()
{
   return  nullptr;
}

//////////////////////////////////////////////////////////////////////////////
/// Return working directory.

std::string TSystem::GetWorkingDirectory() const
{
   return std::string();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the user's home directory.

const char *TSystem::HomeDirectory(const char *)
{
   return nullptr;
}

//////////////////////////////////////////////////////////////////////////////
/// Return the user's home directory.

std::string TSystem::GetHomeDirectory(const char *) const
{
   return std::string();
}

////////////////////////////////////////////////////////////////////////////////
/// Make a file system directory. Returns 0 in case of success and
/// -1 if the directory could not be created (either already exists or
/// illegal path name).
/// If 'recursive' is true, makes parent directories as needed.

int TSystem::mkdir(const char *name, Bool_t recursive)
{
   if (recursive) {
      TString safeName = name; // local copy in case 'name' is output from
                               // TSystem::DirName as it uses static buffers
      TString dirname = GetDirName(safeName.Data());
      if (dirname.IsNull()) {
         // well we should not have to make the root of the file system!
         // (and this avoid infinite recursions!)
         return -1;
      }
      if (AccessPathName(dirname.Data(), kFileExists)) {
         int res = mkdir(dirname.Data(), kTRUE);
         if (res) return res;
      }
      if (!AccessPathName(safeName.Data(), kFileExists)) {
         return -1;
      }
   }

   return MakeDirectory(name);
}

//---- Paths & Files -----------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Base name of a file name. Base name of /user/root is root.

const char *TSystem::BaseName(const char *name)
{
   if (name) {
      if (name[0] == '/' && name[1] == '\0')
         return name;
      char *cp;
      if ((cp = (char *)strrchr(name, '/')))
         return ++cp;
      return name;
   }
   Error("BaseName", "name = 0");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if dir is an absolute pathname.

Bool_t TSystem::IsAbsoluteFileName(const char *dir)
{
   if (dir)
      return dir[0] == '/';
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if 'name' is a file that can be found in the ROOT include
/// path or the current directory.
/// If 'name' contains any ACLiC style information (e.g. trailing +[+][g|O]),
/// it will be striped off 'name'.
/// If fullpath is != 0, the full path to the file is returned in *fullpath,
/// which must be deleted by the caller.

Bool_t TSystem::IsFileInIncludePath(const char *name, char **fullpath)
{
   if (!name || !name[0]) return kFALSE;

   TString aclicMode;
   TString arguments;
   TString io;
   TString realname = SplitAclicMode(name, aclicMode, arguments, io);

   TString fileLocation = GetDirName(realname);

   TString incPath = gSystem->GetIncludePath(); // of the form -Idir1  -Idir2 -Idir3
   incPath.Append(":").Prepend(" ");
   incPath.ReplaceAll(" -I",":");       // of form :dir1 :dir2:dir3
   while ( incPath.Index(" :") != -1 ) {
      incPath.ReplaceAll(" :",":");
   }
   // Remove double quotes around path expressions.
   incPath.ReplaceAll("\":", ":");
   incPath.ReplaceAll(":\"", ":");

   incPath.Prepend(fileLocation+":.:");

   char *actual = Which(incPath,realname);

   if (!actual) {
      return kFALSE;
   } else {
      if (fullpath)
         *fullpath = actual;
      else
         delete [] actual;
      return kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the directory name in pathname. DirName of /user/root is /user.
/// In case no dirname is specified "." is returned.

const char *TSystem::DirName(const char *pathname)
{
   auto res = GetDirName(pathname);
   if (res.IsNull() || (res == "."))
      return ".";

   R__LOCKGUARD2(gSystemMutex);

   static Ssiz_t len = 0;
   static char *buf = nullptr;
   if (res.Length() >= len) {
      if (buf) delete [] buf;
      len = res.Length() + 50;
      buf = new char [len];
   }
   strncpy(buf, res.Data(), len);
   return buf;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the directory name in pathname.
/// DirName of /user/root is /user.
/// DirName of /user/root/ is also /user.
/// In case no dirname is specified "." is returned.

TString TSystem::GetDirName(const char *pathname)
{
   if (!pathname || !strchr(pathname, '/'))
      return ".";

   auto pathlen = strlen(pathname);

   const char *r = pathname + pathlen - 1;
   // First skip the trailing '/'
   while ((r > pathname) && (*r == '/'))
      --r;
   // Then find the next non slash
   while ((r > pathname) && (*r != '/'))
      --r;

   // Then skip duplicate slashes
   // Note the 'r>buf' is a strict comparison to allows '/topdir' to return '/'
   while ((r > pathname) && (*r == '/'))
      --r;
   // If all was cut away, we encountered a rel. path like 'subdir/'
   // and ended up at '.'.
   if ((r == pathname) && (*r != '/'))
      return ".";

   return TString(pathname, r + 1 - pathname);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert from a local pathname to a Unix pathname. E.g. from `\user\root` to
/// `/user/root`.

const char *TSystem::UnixPathName(const char *name)
{
   return name;
}

////////////////////////////////////////////////////////////////////////////////
/// Concatenate a directory and a file name. User must delete returned string.

char *TSystem::ConcatFileName(const char *dir, const char *name)
{
   TString nameString(name);
   PrependPathName(dir, nameString);
   return StrDup(nameString.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Concatenate a directory and a file name.

const char *TSystem::PrependPathName(const char *, TString&)
{
   AbstractMethod("PrependPathName");
   return nullptr;
}


//---- Paths & Files -----------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Expand a pathname getting rid of special shell characters like ~.$, etc.
/// For Unix/Win32 compatibility use $(XXX) instead of $XXX when using
/// environment variables in a pathname. If compatibility is not an issue
/// you can use on Unix directly $XXX. This is a protected function called
/// from the OS specific system classes, like TUnixSystem and TWinNTSystem.
/// Returns the expanded filename or 0 in case of error.

const char *TSystem::ExpandFileName(const char *fname)
{
   const int kBufSize = kMAXPATHLEN;
   TTHREAD_TLS_ARRAY(char, kBufSize, xname);

   Bool_t res = ExpandFileName(fname, xname, kBufSize);
   if (res)
      return nullptr;
   else
      return xname;
}

//////////////////////////////////////////////////////////////////////////////
/// Expand a pathname getting rid of special shell characters like ~.$, etc.
/// This function is analogous to ExpandFileName(const char *), except that
/// it receives a TString reference of the pathname to be expanded.
/// Returns kTRUE in case of error and kFALSE otherwise.

Bool_t TSystem::ExpandFileName(TString &fname)
{
   const int kBufSize = kMAXPATHLEN;
   char xname[kBufSize];

   Bool_t res = ExpandFileName(fname.Data(), xname, kBufSize);
   if (!res)
      fname = xname;

   return res;
}

////////////////////////////////////////////////////////////////////////////
/// Private method for pathname expansion.
/// Returns kTRUE in case of error and kFALSE otherwise.

Bool_t TSystem::ExpandFileName(const char *fname, char *xname, const int kBufSize)
{
   int n, ier, iter, lx, ncopy;
   char *inp, *out, *x, *t, *buff;
   const char *b, *c, *e;
   const char *p;
   buff = new char[kBufSize * 4];

   iter = 0; xname[0] = 0; inp = buff + kBufSize; out = inp + kBufSize;
   inp[-1] = ' '; inp[0] = 0; out[-1] = ' ';
   c = fname + strspn(fname, " \t\f\r");
   //VP  if (isalnum(c[0])) { strcpy(inp, WorkingDirectory()); strcat(inp, "/"); } // add $cwd

   strlcat(inp, c, kBufSize);

again:
   iter++; c = inp; ier = 0;
   x = out; x[0] = 0;

   p = 0; e = 0;
   if (c[0] == '~' && c[1] == '/') { // ~/ case
      std::string hd = GetHomeDirectory();
      p = hd.c_str();
      e = c + 1;
      if (p) {                         // we have smth to copy
         strlcpy(x, p, kBufSize);
         x += strlen(p);
         c = e;
      } else {
         ++ier;
         ++c;
      }
   } else if (c[0] == '~' && c[1] != '/') { // ~user case
      n = strcspn(c+1, "/ ");
      assert((n+1) < kBufSize && "This should have been prevented by the truncation 'strlcat(inp, c, kBufSize)'");
      // There is no overlap here as the buffer is segment in 4 strings of at most kBufSize
      (void)strlcpy(buff, c+1, n+1); // strlcpy copy 'size-1' characters.
      std::string hd = GetHomeDirectory(buff);
      e = c+1+n;
      if (!hd.empty()) {                   // we have smth to copy
         p = hd.c_str();
         strlcpy(x, p, kBufSize);
         x += strlen(p);
         c = e;
      } else {
         x++[0] = c[0];
         //++ier;
         ++c;
      }
   }

   for ( ; c[0]; c++) {

      p = 0; e = 0;

      if (c[0] == '.' && c[1] == '/' && c[-1] == ' ') { // $cwd
         std::string wd = GetWorkingDirectory();
         strlcpy(buff, wd.c_str(), kBufSize);
         p = buff;
         e = c + 1;
      }
      if (p) {                          // we have smth to copy */
         strlcpy(x, p, kBufSize); x += strlen(p); c = e-1; continue;
      }

      if (c[0] != '$') {                // not $, simple copy
         x++[0] = c[0];
      } else {                          // we have a $
         b = c+1;
         if (c[1] == '(') b++;
         if (c[1] == '{') b++;
         if (b[0] == '$')
            e = b+1;
         else
            for (e = b; isalnum(e[0]) || e[0] == '_'; e++) ;
         buff[0] = 0; strncat(buff, b, e-b);
         p = Getenv(buff);
         if (!p) {                      // too bad, try UPPER case
            for (t = buff; (t[0] = toupper(t[0])); t++) ;
            p = Getenv(buff);
         }
         if (!p) {                      // too bad, try Lower case
            for (t = buff; (t[0] = tolower(t[0])); t++) ;
            p = Getenv(buff);
         }
         if (!p && !strcmp(buff, "cwd")) { // it is $cwd
            std::string wd = GetWorkingDirectory();
            strlcpy(buff, wd.c_str(), kBufSize);
            p = buff;
         }
         if (!p && !strcmp(buff, "$")) { // it is $$ (replace by GetPid())
            snprintf(buff,kBufSize*4, "%d", GetPid());
            p = buff;
         }
         if (!p) {                      // too bad, nothing can help
#ifdef WIN32
            // if we're on windows, we can have \\SomeMachine\C$ - don't
            // complain about that, if '$' is followed by nothing or a
            // path delimiter.
            if (c[1] && c[1]!='\\' && c[1]!=';' && c[1]!='/')
               ier++;
#else
            ier++;
#endif
            x++[0] = c[0];
         } else {                       // It is OK, copy result
            int lp = strlen(p);
            if (lp >= kBufSize) {
               // make sure lx will be >= kBufSize (see below)
               strlcpy(x, p, kBufSize);
               x += kBufSize;
               break;
            }
            strcpy(x,p);
            x += lp;
            c = (b==c+1) ? e-1 : e;
         }
      }
   }

   x[0] = 0; lx = x - out;
   if (ier && iter < 3) { strlcpy(inp, out, kBufSize); goto again; }
   ncopy = (lx >= kBufSize) ? kBufSize-1 : lx;
   xname[0] = 0; strncat(xname, out, ncopy);

   delete[] buff;

   if (ier || ncopy != lx) {
      ::Error("TSystem::ExpandFileName", "input: %s, output: %s", fname, xname);
      return kTRUE;
   }

   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Expand a pathname getting rid of special shell characters like ~.$, etc.
/// For Unix/Win32 compatibility use $(XXX) instead of $XXX when using
/// environment variables in a pathname. If compatibility is not an issue
/// you can use on Unix directly $XXX.

Bool_t TSystem::ExpandPathName(TString&)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Expand a pathname getting rid of special shell characters like ~.$, etc.
/// For Unix/Win32 compatibility use $(XXX) instead of $XXX when using
/// environment variables in a pathname. If compatibility is not an issue
/// you can use on Unix directly $XXX. The user must delete returned string.

char *TSystem::ExpandPathName(const char *)
{
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns FALSE if one can access a file using the specified access mode.
/// The file name must not contain any special shell characters line ~ or $,
/// in those cases first call ExpandPathName().
/// Attention, bizarre convention of return value!!

Bool_t TSystem::AccessPathName(const char *, EAccessMode)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns TRUE if the url in 'path' points to the local file system.
/// This is used to avoid going through the NIC card for local operations.

Bool_t TSystem::IsPathLocal(const char *path)
{
   Bool_t localPath = kTRUE;

   TUrl url(path);
   if (strlen(url.GetHost()) > 0) {
      // Check locality
      localPath = kFALSE;
      TInetAddress a(gSystem->GetHostByName(url.GetHost()));
      TInetAddress b(gSystem->GetHostByName(gSystem->HostName()));
      if (!strcmp(a.GetHostName(), b.GetHostName()) ||
          !strcmp(a.GetHostAddress(), b.GetHostAddress())) {
         // Host OK
         localPath = kTRUE;
         // Check the user if specified
         if (strlen(url.GetUser()) > 0) {
            UserGroup_t *u = gSystem->GetUserInfo();
            if (u) {
               if (strcmp(u->fUser, url.GetUser()))
                  // Requested a different user
                  localPath = kFALSE;
               delete u;
            }
         }
      }
   }
   // Done
   return localPath;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a file. If overwrite is true and file already exists the
/// file will be overwritten. Returns 0 when successful, -1 in case
/// of file open failure, -2 in case the file already exists and overwrite
/// was false and -3 in case of error during copy.

int TSystem::CopyFile(const char *, const char *, Bool_t)
{
   AbstractMethod("CopyFile");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Rename a file.

int TSystem::Rename(const char *, const char *)
{
   AbstractMethod("Rename");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a link from file1 to file2.

int TSystem::Link(const char *, const char *)
{
   AbstractMethod("Link");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a symbolic link from file1 to file2.

int TSystem::Symlink(const char *, const char *)
{
   AbstractMethod("Symlink");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Unlink, i.e. remove, a file.
///
/// If the file is currently open by the current or another process, the behavior of this function is
/// implementation-defined (in particular, POSIX systems unlink the file name, while Windows does not allow the
/// file to be deleted and the operation is a no-op).

int TSystem::Unlink(const char *)
{
   AbstractMethod("Unlink");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get info about a file: id, size, flags, modification time.
///  - Id      is (statbuf.st_dev << 24) + statbuf.st_ino
///  - Size    is the file size
///  - Flags   is file type: 0 is regular file, bit 0 set executable,
///                          bit 1 set directory, bit 2 set special file
///                          (socket, fifo, pipe, etc.)
/// Modtime is modification time.
/// The function returns 0 in case of success and 1 if the file could
/// not be stat'ed.

int TSystem::GetPathInfo(const char *path, Long_t *id, Long_t *size,
                         Long_t *flags, Long_t *modtime)
{
   Long64_t lsize;

   int res = GetPathInfo(path, id, &lsize, flags, modtime);

   if (res == 0 && size) {
      if (sizeof(Long_t) == 4 && lsize > kMaxInt) {
         Error("GetPathInfo", "file %s > 2 GB, use GetPathInfo() with Long64_t size", path);
         *size = kMaxInt;
      } else {
         *size = (Long_t)lsize;
      }
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Get info about a file: id, size, flags, modification time.
///  - Id      is (statbuf.st_dev << 24) + statbuf.st_ino
///  - Size    is the file size
///  - Flags   is file type: 0 is regular file, bit 0 set executable,
///                          bit 1 set directory, bit 2 set special file
///                          (socket, fifo, pipe, etc.)
/// Modtime is modification time.
/// The function returns 0 in case of success and 1 if the file could
/// not be stat'ed.

int TSystem::GetPathInfo(const char *path, Long_t *id, Long64_t *size,
                         Long_t *flags, Long_t *modtime)
{
   FileStat_t buf;

   int res = GetPathInfo(path, buf);

   if (res == 0) {
      if (id)
         *id = (buf.fDev << 24) + buf.fIno;
      if (size)
         *size = buf.fSize;
      if (modtime)
         *modtime = buf.fMtime;
      if (flags) {
         *flags = 0;
         if (buf.fMode & (kS_IXUSR|kS_IXGRP|kS_IXOTH))
            *flags |= 1;
         if (R_ISDIR(buf.fMode))
            *flags |= 2;
         if (!R_ISREG(buf.fMode) && !R_ISDIR(buf.fMode))
            *flags |= 4;
      }
   }

   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Get info about a file. Info is returned in the form of a FileStat_t
/// structure (see TSystem.h).
/// The function returns 0 in case of success and 1 if the file could
/// not be stat'ed.

int TSystem::GetPathInfo(const char *, FileStat_t &)
{
   AbstractMethod("GetPathInfo(const char *, FileStat_t&)");
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get info about a file system: fs type, block size, number of blocks,
/// number of free blocks.

int TSystem::GetFsInfo(const char *, Long_t *, Long_t *, Long_t *, Long_t *)
{
   AbstractMethod("GetFsInfo");
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a user configured or systemwide directory to create
/// temporary files in.

const char *TSystem::TempDirectory() const
{
   AbstractMethod("TempDirectory");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a secure temporary file by appending a unique
/// 6 letter string to base. The file will be created in
/// a standard (system) directory or in the directory
/// provided in dir. The full filename is returned in base
/// and a filepointer is returned for safely writing to the file
/// (this avoids certain security problems). Returns 0 in case
/// of error.

FILE *TSystem::TempFileName(TString &, const char *)
{
   AbstractMethod("TempFileName");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the file permission bits. Returns -1 in case or error, 0 otherwise.

int TSystem::Chmod(const char *, UInt_t)
{
   AbstractMethod("Chmod");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the process file creation mode mask.

int TSystem::Umask(Int_t)
{
   AbstractMethod("Umask");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the a files modification and access times. If actime = 0 it will be
/// set to the modtime. Returns 0 on success and -1 in case of error.

int TSystem::Utime(const char *, Long_t, Long_t)
{
   AbstractMethod("Utime");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Find location of file in a search path. Return value points to TString for
/// compatibility with Which(const char *, const char *, EAccessMode).
/// Returns 0 in case file is not found.

const char *TSystem::FindFile(const char *, TString&, EAccessMode)
{
   AbstractMethod("FindFile");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Find location of file in a search path. User must delete returned string.
/// Returns 0 in case file is not found.

char *TSystem::Which(const char *search, const char *wfil, EAccessMode mode)
{
   TString wfilString(wfil);
   FindFile(search, wfilString, mode);
   if (wfilString.IsNull())
      return nullptr;
   return StrDup(wfilString.Data());
}

//---- Users & Groups ----------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Returns the user's id. If user = 0, returns current user's id.

Int_t TSystem::GetUid(const char * /*user*/)
{
   AbstractMethod("GetUid");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the effective user id. The effective id corresponds to the
/// set id bit on the file being executed.

Int_t TSystem::GetEffectiveUid()
{
   AbstractMethod("GetEffectiveUid");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the group's id. If group = 0, returns current user's group.

Int_t TSystem::GetGid(const char * /*group*/)
{
   AbstractMethod("GetGid");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the effective group id. The effective group id corresponds
/// to the set id bit on the file being executed.

Int_t TSystem::GetEffectiveGid()
{
   AbstractMethod("GetEffectiveGid");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns all user info in the UserGroup_t structure. The returned
/// structure must be deleted by the user. In case of error 0 is returned.

UserGroup_t *TSystem::GetUserInfo(Int_t /*uid*/)
{
   AbstractMethod("GetUserInfo");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns all user info in the UserGroup_t structure. If user = 0, returns
/// current user's id info. The returned structure must be deleted by the
/// user. In case of error 0 is returned.

UserGroup_t *TSystem::GetUserInfo(const char * /*user*/)
{
   AbstractMethod("GetUserInfo");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns all group info in the UserGroup_t structure. The only active
/// fields in the UserGroup_t structure for this call are:
///  - fGid and fGroup
/// The returned structure must be deleted by the user. In case of
/// error 0 is returned.

UserGroup_t *TSystem::GetGroupInfo(Int_t /*gid*/)
{
   AbstractMethod("GetGroupInfo");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns all group info in the UserGroup_t structure. The only active
/// fields in the UserGroup_t structure for this call are:
///  - fGid and fGroup
/// If group = 0, returns current user's group. The returned structure
/// must be deleted by the user. In case of error 0 is returned.

UserGroup_t *TSystem::GetGroupInfo(const char * /*group*/)
{
   AbstractMethod("GetGroupInfo");
   return nullptr;
}

//---- environment manipulation ------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Set environment variable.

void TSystem::Setenv(const char *, const char *)
{
   AbstractMethod("Setenv");
}

////////////////////////////////////////////////////////////////////////////////
/// Unset environment variable.

void TSystem::Unsetenv(const char *name)
{
   Setenv(name, "");
}

////////////////////////////////////////////////////////////////////////////////
/// Get environment variable.

const char *TSystem::Getenv(const char *)
{
   AbstractMethod("Getenv");
   return nullptr;
}

//---- System Logging ----------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Open connection to system log daemon. For the use of the options and
/// facility see the Unix openlog man page.

void TSystem::Openlog(const char *, Int_t, ELogFacility)
{
   AbstractMethod("Openlog");
}

////////////////////////////////////////////////////////////////////////////////
/// Send mess to syslog daemon. Level is the logging level and mess the
/// message that will be written on the log.

void TSystem::Syslog(ELogLevel, const char *)
{
   AbstractMethod("Syslog");
}

////////////////////////////////////////////////////////////////////////////////
/// Close connection to system log daemon.

void TSystem::Closelog()
{
   AbstractMethod("Closelog");
}

//---- Standard output redirection ---------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Redirect standard output (stdout, stderr) to the specified file.
/// If the file argument is 0 the output is set again to stderr, stdout.
/// The second argument specifies whether the output should be added to the
/// file ("a", default) or the file be truncated before ("w").
/// The implementations of this function save internally the current state into
/// a static structure.
///
/// The call can be made reentrant by specifying the opaque structure pointed
/// by 'h', which is filled with the relevant information. The handle 'h'
/// obtained on the first call must then be used in any subsequent call,
/// included ShowOutput, to display the redirected output.
/// Returns 0 on success, -1 in case of error.

Int_t TSystem::RedirectOutput(const char *, const char *, RedirectHandle_t *)
{
   AbstractMethod("RedirectOutput");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Display the content associated with the redirection described by the
/// opaque handle 'h'.

void TSystem::ShowOutput(RedirectHandle_t *h)
{
   // Check input ...
   if (!h) {
      Error("ShowOutput", "handle not specified");
      return;
   }

   // ... and file access
   if (gSystem->AccessPathName(h->fFile, kReadPermission)) {
      Error("ShowOutput", "file '%s' cannot be read", h->fFile.Data());
      return;
   }

   // Open the file
   FILE *f = nullptr;
   if (!(f = fopen(h->fFile.Data(), "r"))) {
      Error("ShowOutput", "file '%s' cannot be open", h->fFile.Data());
      return;
   }

   // Determine the number of bytes to be read from the file.
   off_t ltot = lseek(fileno(f), (off_t) 0, SEEK_END);
   Int_t begin = (h->fReadOffSet > 0 && h->fReadOffSet < ltot) ? h->fReadOffSet : 0;
   lseek(fileno(f), (off_t) begin, SEEK_SET);
   Int_t left = ltot - begin;

   // Now readout from file
   const Int_t kMAXBUF = 16384;
   char buf[kMAXBUF];
   Int_t wanted = (left > kMAXBUF-1) ? kMAXBUF-1 : left;
   Int_t len;
   do {
      while ((len = read(fileno(f), buf, wanted)) < 0 &&
               TSystem::GetErrno() == EINTR)
         TSystem::ResetErrno();

      if (len < 0) {
         SysError("ShowOutput", "error reading log file");
         break;
      }

      // Null-terminate
      buf[len] = 0;
      fprintf(stderr,"%s", buf);

      // Update counters
      left -= len;
      wanted = (left > kMAXBUF) ? kMAXBUF : left;

   } while (len > 0 && left > 0);

   // Do not display twice the same thing
   h->fReadOffSet = ltot;
   fclose(f);
}

//---- Dynamic Loading ---------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Add a new directory to the dynamic path.

void TSystem::AddDynamicPath(const char *)
{
   AbstractMethod("AddDynamicPath");
}

////////////////////////////////////////////////////////////////////////////////
/// Return the dynamic path (used to find shared libraries).

const char *TSystem::GetDynamicPath()
{
   AbstractMethod("GetDynamicPath");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the dynamic path to a new value.
/// If the value of 'path' is zero, the dynamic path is reset to its
/// default value.

void TSystem::SetDynamicPath(const char *)
{
   AbstractMethod("SetDynamicPath");
}


////////////////////////////////////////////////////////////////////////////////
/// Figure out if left and right points to the same
/// object in the file system.

static bool R__MatchFilename(const char *left, const char *right)
{
   if (left == right) return kTRUE;

   if (left==0 || right==0) return kFALSE;

   if ( (strcmp(right,left)==0) ) {
      return kTRUE;
   }

#ifdef G__WIN32

   char leftname[_MAX_PATH];
   char rightname[_MAX_PATH];
   _fullpath( leftname, left, _MAX_PATH );
   _fullpath( rightname, right, _MAX_PATH );
   return ((stricmp(leftname, rightname)==0));
#else
   struct stat rightBuf;
   struct stat leftBuf;
   return (   ( 0 == stat( left, & leftBuf ) )
       && ( 0 == stat( right, & rightBuf ) )
       && ( leftBuf.st_dev == rightBuf.st_dev )     // Files on same device
       && ( leftBuf.st_ino == rightBuf.st_ino )     // Files on same inode (but this is not unique on AFS so we need the next 2 test
       && ( leftBuf.st_size == rightBuf.st_size )   // Files of same size
       && ( leftBuf.st_mtime == rightBuf.st_mtime ) // Files modified at the same time
           );
#endif
}


////////////////////////////////////////////////////////////////////////////////
/// Load a shared library. Returns 0 on successful loading, 1 in
/// case lib was already loaded, -1 in case lib does not exist
/// or in case of error and -2 in case of version mismatch.
/// When entry is specified the loaded lib is
/// searched for this entry point (return -1 when entry does not exist,
/// 0 otherwise). When the system flag is kTRUE, the library is considered
/// a permanent system library that should not be unloaded during the
/// course of the session.

int TSystem::Load(const char *module, const char *entry, Bool_t system)
{
   // don't load libraries that have already been loaded
   TString libs( GetLibraries() );
   TString l(BaseName(module));

   Ssiz_t idx = l.Last('.');
   if (idx != kNPOS) {
      l.Remove(idx+1);
   }
   for (idx = libs.Index(l); idx != kNPOS; idx = libs.Index(l,idx+1)) {
      // The libs contains the sub-string 'l', let's make sure it is
      // not just part of a larger name.
      if (idx == 0 || libs[idx-1] == '/' || libs[idx-1] == '\\') {
         Ssiz_t len = libs.Length();
         idx += l.Length();
         if (!l.EndsWith(".") && libs[idx]=='.')
            idx++;
         // Skip the soversion.
         while (idx < len && isdigit(libs[idx])) {
            ++idx;
            // No need to test for len here, at worse idx==len and lib[idx]=='\0'
            if (libs[idx] == '.') {
               ++idx;
            }
         }
         while (idx < len && libs[idx] != '.') {
            if (libs[idx] == ' ' || idx+1 == len) {
               return 1;
            }
            ++idx;
         }
      }
   }
   if (l[l.Length()-1] == '.') {
      l.Remove(l.Length()-1);
   }
   if (l.BeginsWith("lib")) {
      l.Replace(0, 3, "-l");
      for(idx = libs.Index(l); idx != kNPOS; idx = libs.Index(l,idx+1)) {
         if ((idx == 0 || libs[idx-1] == ' ') &&
             (libs[idx+l.Length()] == ' ' || libs[idx+l.Length()] == 0)) {
            return 1;
         }
      }
   }

   char *path = DynamicPathName(module);

   int ret = -1;
   if (path) {
      // load any dependent libraries
      TString deplibs = gInterpreter->GetSharedLibDeps(path);
      if (!deplibs.IsNull()) {
         TString delim(" ");
         TObjArray *tokens = deplibs.Tokenize(delim);
         for (Int_t i = tokens->GetEntriesFast()-1; i > 0; i--) {
            const char *deplib = ((TObjString*)tokens->At(i))->GetName();
            if (strcmp(module,deplib)==0) {
               continue;
            }
            if (gDebug > 0)
               Info("Load", "loading dependent library %s for library %s",
                    deplib, ((TObjString*)tokens->At(0))->GetName());
            if ((ret = Load(deplib, "", system)) < 0) {
               delete tokens;
               delete [] path;
               return ret;
            }
         }
         delete tokens;
      }
      if (!system) {
         // Mark the library in $ROOTSYS/lib as system.
         TString dirname = GetDirName(path);
         system = R__MatchFilename(TROOT::GetLibDir(), dirname.Data());

         if (!system) {
            system = R__MatchFilename(TROOT::GetBinDir(), dirname.Data());
         }
      }

      gLibraryVersionIdx++;
      if (gLibraryVersionIdx == gLibraryVersionMax) {
         gLibraryVersionMax *= 2;
         gLibraryVersion = TStorage::ReAllocInt(gLibraryVersion, gLibraryVersionMax, gLibraryVersionIdx);
      }
      ret = gInterpreter->Load(path, system);
      if (ret < 0) ret = -1;
      if (gDebug > 0)
         Info("Load", "loaded library %s, status %d", path, ret);
      if (ret == 0 && gLibraryVersion[gLibraryVersionIdx]) {
         int v = TROOT::ConvertVersionCode2Int(gLibraryVersion[gLibraryVersionIdx]);
         Error("Load", "version mismatch, %s = %d, ROOT = %d",
               path, v, gROOT->GetVersionInt());
         ret = -2;
         gLibraryVersion[gLibraryVersionIdx] = 0;
      }
      gLibraryVersionIdx--;
      delete [] path;
   }

   if (!entry || !entry[0] || ret < 0) return ret;

   Func_t f = DynFindSymbol(module, entry);
   if (f) return 0;
   return -1;
}

///////////////////////////////////////////////////////////////////////////////
/// Load all libraries known to ROOT via the rootmap system.
/// Returns the number of top level libraries successfully loaded.

UInt_t TSystem::LoadAllLibraries()
{
   UInt_t nlibs = 0;

   TEnv* mapfile = gInterpreter->GetMapfile();
   if (!mapfile || !mapfile->GetTable()) return 0;

   std::set<std::string> loadedlibs;
   std::set<std::string> failedlibs;

   TEnvRec* rec = 0;
   TIter iEnvRec(mapfile->GetTable());
   while ((rec = (TEnvRec*) iEnvRec())) {
      TString libs = rec->GetValue();
      TString lib;
      Ssiz_t pos = 0;
      while (libs.Tokenize(lib, pos)) {
         // check that none of the libs failed to load
         if (failedlibs.find(lib.Data()) != failedlibs.end()) {
            // don't load it or any of its dependencies
            libs = "";
            break;
         }
      }
      pos = 0;
      while (libs.Tokenize(lib, pos)) {
         // ignore libCore - it's already loaded
         if (lib.BeginsWith("libCore"))
            continue;

         if (loadedlibs.find(lib.Data()) == loadedlibs.end()) {
            // just load the first library - TSystem will do the rest.
            auto res = gSystem->Load(lib);
            if (res >=0) {
               if (res == 0) ++nlibs;
               loadedlibs.insert(lib.Data());
            } else {
               failedlibs.insert(lib.Data());
            }
         }
      }
   }
   return nlibs;
}

////////////////////////////////////////////////////////////////////////////////
/// Find a dynamic library called lib using the system search paths.
/// Appends known extensions if needed. Returned string must be deleted
/// by the user!

char *TSystem::DynamicPathName(const char *lib, Bool_t quiet /*=kFALSE*/)
{
   TString sLib(lib);
   if (FindDynamicLibrary(sLib, quiet))
      return StrDup(sLib);
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Find a dynamic library using the system search paths. lib will be updated
/// to contain the absolute filename if found. Returns lib if found, or NULL
/// if a library called lib was not found.
/// This function does not open the library.

const char *TSystem::FindDynamicLibrary(TString&, Bool_t)
{
   AbstractMethod("FindDynamicLibrary");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Find specific entry point in specified library. Specify "*" for lib
/// to search in all libraries.

Func_t TSystem::DynFindSymbol(const char * /*lib*/, const char *entry)
{
   return (Func_t) gInterpreter->FindSym(entry);
}

////////////////////////////////////////////////////////////////////////////////
/// Unload a shared library.

void TSystem::Unload(const char *module)
{
   char *path;
   if ((path = DynamicPathName(module))) {
      gInterpreter->UnloadFile(path);
      delete [] path;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// List symbols in a shared library.

void TSystem::ListSymbols(const char *, const char *)
{
   AbstractMethod("ListSymbols");
}

////////////////////////////////////////////////////////////////////////////////
/// List all loaded shared libraries. Regexp is a wildcard expression,
/// see TRegexp::MakeWildcard.

void TSystem::ListLibraries(const char *regexp)
{
   TString libs = GetLibraries(regexp);
   TRegexp separator("[^ \\t\\s]+");
   TString s;
   Ssiz_t start = 0, index = 0, end = 0;
   int i = 0;

   Printf(" ");
   Printf("Loaded shared libraries");
   Printf("=======================");

   while ((start < libs.Length()) && (index != kNPOS)) {
      index = libs.Index(separator, &end, start);
      if (index >= 0) {
         s = libs(index, end);
         if (s.BeginsWith("-")) {
            if (s.BeginsWith("-l")) {
               Printf("%s", s.Data());
               i++;
            }
         } else {
            Printf("%s", s.Data());
            i++;
         }
      }
      start += end+1;
   }

   Printf("-----------------------");
   Printf("%d libraries loaded", i);
   Printf("=======================");
}

////////////////////////////////////////////////////////////////////////////////
/// Return the thread local storage for the custom last error message

TString &TSystem::GetLastErrorString()
{
   TTHREAD_TLS_DECL( TString, gLastErrorString);
   return gLastErrorString;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the thread local storage for the custom last error message

const TString &TSystem::GetLastErrorString() const
{
   return const_cast<TSystem*>(this)->GetLastErrorString();
}

////////////////////////////////////////////////////////////////////////////////
/// Get list of shared libraries loaded at the start of the executable.
/// Returns 0 in case list cannot be obtained or in case of error.

const char *TSystem::GetLinkedLibraries()
{
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a space separated list of loaded shared libraries.
/// Regexp is a wildcard expression, see TRegexp::MakeWildcard.
/// This list is of a format suitable for a linker, i.e it may contain
/// -Lpathname and/or -lNameOfLib.
/// Option can be any of:
///  - S: shared libraries loaded at the start of the executable, because
///       they were specified on the link line.
///  - D: shared libraries dynamically loaded after the start of the program.
///  - L: this option is ignored, and available for backward compatibility.

const char *TSystem::GetLibraries(const char *regexp, const char *options,
                                  Bool_t isRegexp)
{
   fListLibs.Clear();

   TString libs;
   TString opt(options);
   Bool_t so2dylib = (opt.First('L') != kNPOS);
   if (so2dylib)
      opt.ReplaceAll("L", "");

   if (opt.IsNull() || opt.First('D') != kNPOS)
      libs += gInterpreter->GetSharedLibs();

   // Cint currently register all libraries that
   // are loaded and have a dictionary in them, this
   // includes all the libraries that are included
   // in the list of (hard) linked libraries.

   TString slinked;
   const char *linked;
   if ((linked = GetLinkedLibraries())) {
      if (fLinkedLibs != LINKEDLIBS) {
         // This is not the default value, we need to keep the custom part.
         TString custom = fLinkedLibs;
         custom.ReplaceAll(LINKEDLIBS,linked);
         if (custom == fLinkedLibs) {
            // no replacement done, let's append linked
            slinked.Append(linked);
            slinked.Append(" ");
         }
         slinked.Append(custom);
      } else {
         slinked.Append(linked);
      }
   } else {
      slinked.Append(fLinkedLibs);
   }

   if (opt.IsNull() || opt.First('S') != kNPOS) {
      // We are done, the statically linked libraries are already included.
      if (libs.Length() == 0) {
         libs = slinked;
      } else {
         // We need to add the missing linked library

         static TString lastLinked;
         static TString lastAddMissing;
         if ( lastLinked != slinked ) {
            // Recalculate only if there was a change.
            static TRegexp separator("[^ \\t\\s]+");
            lastLinked = slinked;
            lastAddMissing.Clear();

            Ssiz_t start, index, end;
            start = index = end = 0;

            while ((start < slinked.Length()) && (index != kNPOS)) {
               index = slinked.Index(separator,&end,start);
               if (index >= 0) {
                  TString sub = slinked(index,end);
                  if (sub[0]=='-' && sub[1]=='L') {
                     lastAddMissing.Prepend(" ");
                     lastAddMissing.Prepend(sub);
                  } else {
                     if (libs.Index(sub) == kNPOS) {
                        lastAddMissing.Prepend(" ");
                        lastAddMissing.Prepend(sub);
                     }
                  }
               }
               start += end+1;
            }
         }
         libs.Prepend(lastAddMissing);
      }
   } else if (libs.Length() != 0) {
      // Let remove the statically linked library
      // from the list.
      static TRegexp separator("[^ \\t\\s]+");
      Ssiz_t start, index, end;
      start = index = end = 0;

      while ((start < slinked.Length()) && (index != kNPOS)) {
         index = slinked.Index(separator,&end,start);
         if (index >= 0) {
            TString sub = slinked(index,end);
            if (sub[0]!='-' && sub[1]!='L') {
               libs.ReplaceAll(sub,"");
            }
         }
         start += end+1;
      }
      libs = libs.Strip(TString::kBoth);
   }

   // Select according to regexp
   if (regexp && *regexp) {
      static TRegexp separator("[^ \\t\\s]+");
      TRegexp user_re(regexp, kTRUE);
      TString s;
      Ssiz_t start, index, end;
      start = index = end = 0;

      while ((start < libs.Length()) && (index != kNPOS)) {
         index = libs.Index(separator,&end,start);
         if (index >= 0) {
            s = libs(index,end);
            if ((isRegexp && s.Index(user_re) != kNPOS) ||
                (!isRegexp && s.Index(regexp) != kNPOS)) {
               if (!fListLibs.IsNull())
                  fListLibs.Append(" ");
               fListLibs.Append(s);
            }
         }
         start += end+1;
      }
   } else
      fListLibs = libs;

#if defined(R__MACOSX)
// We need to remove the libraries that are dynamically loaded and not linked
{
   TString libs2 = fListLibs;
   TString maclibs;

   static TRegexp separator("[^ \\t\\s]+");
   static TRegexp dynload("/lib-dynload/");

   Ssiz_t start, index, end;
   start = index = end = 0;

   while ((start < libs2.Length()) && (index != kNPOS)) {
      index = libs2.Index(separator, &end, start);
      if (index >= 0) {
         TString s = libs2(index, end);
         if (s.Index(dynload) == kNPOS) {
            if (!maclibs.IsNull()) maclibs.Append(" ");
            maclibs.Append(s);
         }
      }
      start += end+1;
   }
   fListLibs = maclibs;
}
#endif

   return fListLibs.Data();
}

//---- RPC ---------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Get Internet Protocol (IP) address of host.

TInetAddress TSystem::GetHostByName(const char *)
{
   AbstractMethod("GetHostByName");
   return TInetAddress();
}

////////////////////////////////////////////////////////////////////////////////
/// Get Internet Protocol (IP) address of remote host and port #.

TInetAddress TSystem::GetPeerName(int)
{
   AbstractMethod("GetPeerName");
   return TInetAddress();
}

////////////////////////////////////////////////////////////////////////////////
/// Get Internet Protocol (IP) address of host and port #.

TInetAddress TSystem::GetSockName(int)
{
   AbstractMethod("GetSockName");
   return TInetAddress();
}

////////////////////////////////////////////////////////////////////////////////
/// Get port # of internet service.

int TSystem::GetServiceByName(const char *)
{
   AbstractMethod("GetServiceByName");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get name of internet service.

char *TSystem::GetServiceByPort(int)
{
   AbstractMethod("GetServiceByPort");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a connection to another host.

int TSystem::OpenConnection(const char *, int, int, const char *)
{
   AbstractMethod("OpenConnection");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Announce TCP/IP service.

int TSystem::AnnounceTcpService(int, Bool_t, int, int)
{
   AbstractMethod("AnnounceTcpService");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Announce UDP service.

int TSystem::AnnounceUdpService(int, int)
{
   AbstractMethod("AnnounceUdpService");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Announce unix domain service.

int TSystem::AnnounceUnixService(int, int)
{
   AbstractMethod("AnnounceUnixService");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Announce unix domain service.

int TSystem::AnnounceUnixService(const char *, int)
{
   AbstractMethod("AnnounceUnixService");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Accept a connection.

int TSystem::AcceptConnection(int)
{
   AbstractMethod("AcceptConnection");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Close socket connection.

void TSystem::CloseConnection(int, Bool_t)
{
   AbstractMethod("CloseConnection");
}

////////////////////////////////////////////////////////////////////////////////
/// Receive exactly length bytes into buffer. Use opt to receive out-of-band
/// data or to have a peek at what is in the buffer (see TSocket).

int TSystem::RecvRaw(int, void *, int, int)
{
   AbstractMethod("RecvRaw");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Send exactly length bytes from buffer. Use opt to send out-of-band
/// data (see TSocket).

int TSystem::SendRaw(int, const void *, int, int)
{
   AbstractMethod("SendRaw");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Receive a buffer headed by a length indicator.

int TSystem::RecvBuf(int, void *, int)
{
   AbstractMethod("RecvBuf");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Send a buffer headed by a length indicator.

int TSystem::SendBuf(int, const void *, int)
{
   AbstractMethod("SendBuf");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set socket option.

int TSystem::SetSockOpt(int, int, int)
{
   AbstractMethod("SetSockOpt");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get socket option.

int TSystem::GetSockOpt(int, int, int*)
{
   AbstractMethod("GetSockOpt");
   return -1;
}

//---- System, CPU and Memory info ---------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Returns static system info, like OS type, CPU type, number of CPUs
/// RAM size, etc into the SysInfo_t structure. Returns -1 in case of error,
/// 0 otherwise.

int TSystem::GetSysInfo(SysInfo_t *) const
{
   AbstractMethod("GetSysInfo");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns cpu load average and load info into the CpuInfo_t structure.
/// Returns -1 in case of error, 0 otherwise. Use sampleTime to set the
/// interval over which the CPU load will be measured, in ms (default 1000).

int TSystem::GetCpuInfo(CpuInfo_t *, Int_t) const
{
   AbstractMethod("GetCpuInfo");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns ram and swap memory usage info into the MemInfo_t structure.
/// Returns -1 in case of error, 0 otherwise.

int TSystem::GetMemInfo(MemInfo_t *) const
{
   AbstractMethod("GetMemInfo");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns cpu and memory used by this process into the ProcInfo_t structure.
/// Returns -1 in case of error, 0 otherwise.

int TSystem::GetProcInfo(ProcInfo_t *) const
{
   AbstractMethod("GetProcInfo");
   return -1;
}

//---- Script Compiler ---------------------------------------------------------

void AssignAndDelete(TString& target, char *tobedeleted)
{
   // Assign the char* value to the TString and then delete it.

   target = tobedeleted;
   delete [] tobedeleted;
}

#ifdef WIN32

static TString R__Exec(const char *cmd)
{
   // Execute a command and return the stdout in a string.

   FILE * f = gSystem->OpenPipe(cmd,"r");
   if (!f) {
      return "";
   }
   TString result;

   char x;
   while ((x = fgetc(f))!=EOF ) {
      if (x=='\n' || x=='\r') break;
      result += x;
   }

   fclose(f);
   return result;
}

static void R__FixLink(TString &cmd)
{
   // Replace the call to 'link' by a full path name call based on where cl.exe is.
   // This prevents us from using inadvertently the link.exe provided by cygwin.

   // check if link is the microsoft one...
   TString res = R__Exec("link 2>&1");
   if (res.Length()) {
      if (res.Contains("Microsoft (R) Incremental Linker"))
         return;
   }
   // else check availability of cygpath...
   res = R__Exec("cygpath . 2>&1");
   if (res.Length()) {
      if (res != ".")
         return;
   }

   res = R__Exec("which cl.exe 2>&1|grep cl|sed 's,cl\\.exe$,link\\.exe,' 2>&1");
   if (res.Length()) {
      res = R__Exec(Form("cygpath -w '%s' 2>&1",res.Data()));
      if (res.Length()) {
         cmd.ReplaceAll(" link ",Form(" \"%s\" ",res.Data()));
      }
   }
}
#endif

#if defined(__CYGWIN__)
static void R__AddPath(TString &target, const TString &path) {
   if (path.Length() > 2 && path[1]==':') {
      target += TString::Format("/cygdrive/%c",path[0]) + path(2,path.Length()-2);
   } else {
      target += path;
   }
}
#else
static void R__AddPath(TString &target, const TString &path) {
   target += path;
}
#endif

static void R__WriteDependencyFile(const TString & build_loc, const TString &depfilename, const TString &filename, const TString &library, const TString &libname,
                                   const TString &extension, const char *version_var_prefix, const TString &includes, const TString &defines, const TString &incPath)
{
   // Generate the dependency via standard output, not searching the
   // standard include directories,

#ifndef WIN32
   const char * stderrfile = "/dev/null";
#else
   TString stderrfile;
   AssignAndDelete( stderrfile, gSystem->ConcatFileName(build_loc,"stderr.tmp") );
#endif
   TString bakdepfilename = depfilename + ".bak";

#ifdef WIN32
   TString touch = "echo # > "; touch += "\"" + depfilename + "\"";
#else
   TString touch = "echo > "; touch += "\"" + depfilename + "\"";
#endif
   TString builddep = "rmkdepend";
   gSystem->PrependPathName(TROOT::GetBinDir(), builddep);
   builddep += " \"-f";
   builddep += depfilename;
   builddep += "\" -o_" + extension + "." + gSystem->GetSoExt() + " ";
   if (build_loc.BeginsWith(gSystem->WorkingDirectory())) {
      Int_t len = strlen(gSystem->WorkingDirectory());
      if ( build_loc.Length() > (len+1) ) {
         builddep += " \"-p";
         if (build_loc[len] == '/' || build_loc[len+1] != '\\' ) {
            // Since the path is now ran through TSystem::ExpandPathName the single \ is also possible.
            R__AddPath(builddep, build_loc.Data() + len + 1 );
         } else {
            // Case of dir\\name
            R__AddPath(builddep, build_loc.Data() + len + 2 );
         }
         builddep += "/\" ";
      }
   } else {
      builddep += " \"-p";
      R__AddPath(builddep, build_loc);
      builddep += "/\" ";
   }
   builddep += " -Y -- ";
   TString rootsysInclude = TROOT::GetIncludeDir();
   builddep += " \"-I"+rootsysInclude+"\" "; // cflags
   builddep += includes;
   builddep += defines;
   builddep += " -- \"";
   builddep += filename;
   builddep += "\" ";
   TString targetname;
   if (library.BeginsWith(gSystem->WorkingDirectory())) {
      Int_t len = strlen(gSystem->WorkingDirectory());
      if ( library.Length() > (len+1) ) {
         if (library[len] == '/' || library[len+1] != '\\' ) {
            targetname = library.Data() + len + 1;
         } else {
            targetname = library.Data() + len + 2;
         }
      } else {
         targetname = library;
      }
   } else {
      targetname = library;
   }
   builddep += " \"";
   builddep += "-t";
   R__AddPath(builddep, targetname);
   builddep += "\" > ";
   builddep += stderrfile;
   builddep += " 2>&1 ";

   TString adddictdep = "echo ";
   R__AddPath(adddictdep,targetname);
   adddictdep += ": ";
#if defined(R__HAS_CLING_DICTVERSION)
   {
      char *clingdictversion = gSystem->Which(incPath,"clingdictversion.h");
      if (clingdictversion) {
         R__AddPath(adddictdep,clingdictversion);
         adddictdep += " ";
         delete [] clingdictversion;
      } else {
         R__AddPath(adddictdep,rootsysInclude+"/clingdictversion.h ");
      }
   }
#endif
   {
     const char *dictHeaders[] = { "RVersion.h", "RConfig.h", "TClass.h",
       "TDictAttributeMap.h","TInterpreter.h","TROOT.h","TBuffer.h",
       "TMemberInspector.h","TError.h","RtypesImp.h","TIsAProxy.h",
       "TFileMergeInfo.h","TCollectionProxyInfo.h"};

      for (unsigned int h=0; h < sizeof(dictHeaders)/sizeof(dictHeaders[0]); ++h)
      {
         char *rootVersion = gSystem->Which(incPath,dictHeaders[h]);
         if (rootVersion) {
            R__AddPath(adddictdep,rootVersion);
            delete [] rootVersion;
         } else {
            R__AddPath(adddictdep,rootsysInclude + "/" + dictHeaders[h]);
         }
         adddictdep += " ";
      }
   }
   {
      // Add dependency on rootcling.
      char *rootCling = gSystem->Which(gSystem->Getenv("PATH"),"rootcling");
      if (rootCling) {
         R__AddPath(adddictdep,rootCling);
         adddictdep += " ";
         delete [] rootCling;
      }
   }
   adddictdep += " >> \""+depfilename+"\"";

   TString addversiondep( "echo ");
   addversiondep += libname + version_var_prefix + " \"" + ROOT_RELEASE + "\" >> \""+depfilename+"\"";

   if (gDebug > 4)  {
      ::Info("ACLiC", "%s", touch.Data());
      ::Info("ACLiC", "%s", builddep.Data());
      ::Info("ACLiC", "%s", adddictdep.Data());
   }

   Int_t depbuilt = !gSystem->Exec(touch);
   if (depbuilt) depbuilt = !gSystem->Exec(builddep);
   if (depbuilt) depbuilt = !gSystem->Exec(adddictdep);
   if (depbuilt) depbuilt = !gSystem->Exec(addversiondep);

   if (!depbuilt) {
      ::Warning("ACLiC","Failed to generate the dependency file for %s",
                library.Data());
   } else {
#ifdef WIN32
      gSystem->Unlink(stderrfile);
#endif
      gSystem->Unlink(bakdepfilename);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// This method compiles and loads a shared library containing
/// the code from the file "filename".
///
/// The return value is true (1) in case of success and false (0)
/// in case of error.
///
/// The possible options are:
///    - k : keep the shared library after the session end.
///    - f : force recompilation.
///    - g : compile with debug symbol
///    - O : optimized the code
///    - c : compile only, do not attempt to load the library.
///    - s : silence all informational output
///    - v : output all information output
///    - d : debug ACLiC, keep all the output files.
///    - - : if buildir is set, use a flat structure (see buildir below)
///
/// If library_specified is specified, CompileMacro generates the file
/// "library_specified".soext where soext is the shared library extension for
/// the current platform.
///
/// If build_dir is specified, it is used as an alternative 'root' for the
/// generation of the shared library.  The library is stored in a sub-directories
/// of 'build_dir' including the full pathname of the script unless a flat
/// directory structure is requested ('-' option).  With the '-' option the libraries
/// are created directly in the directory 'build_dir'; in particular this means that
/// 2 scripts with the same name in different source directory will over-write each
/// other's library.
/// See also TSystem::SetBuildDir.
///
/// If dirmode is not zero and we need to create the target directory, the
/// file mode bit will be change to 'dirmode' using chmod.
///
/// If library_specified is not specified, CompileMacro generate a default name
/// for library by taking the name of the file "filename" but replacing the
/// dot before the extension by an underscore and by adding the shared
/// library extension for the current platform.
/// For example on most platform, hsimple.cxx will generate hsimple_cxx.so
///
/// It uses the directive fMakeSharedLibs to create a shared library.
/// If loading the shared library fails, it tries to output a list of missing
/// symbols by creating an executable (on some platforms like OSF, this does
/// not HAVE to be an executable) containing the script. It uses the
/// directive fMakeExe to do so.
/// For both directives, before passing them to TSystem::Exec, it expands the
/// variables $SourceFiles, $SharedLib, $LibName, $IncludePath, $LinkedLibs,
/// $DepLibs, $ExeName and $ObjectFiles. See SetMakeSharedLib() for more
/// information on those variables.
///
/// This method is used to implement the following feature:
///
/// Synopsis:
///
/// The purpose of this addition is to allow the user to use an external
/// compiler to create a shared library from its C++ macro (scripts).
/// Currently in order to execute a script, a user has to type at the root
/// prompt
/// ~~~ {.cpp}
///  .X myfunc.C(arg1,arg2)
/// ~~~
/// We allow them to type:
/// ~~~ {.cpp}
///  .X myfunc.C++(arg1,arg2)
/// ~~~
/// or
/// ~~~ {.cpp}
///  .X myfunc.C+(arg1,arg2)
/// ~~~
/// In which case an external compiler will be called to create a shared
/// library.  This shared library will then be loaded and the function
/// myfunc will be called with the two arguments.  With '++' the shared library
/// is always recompiled.  With '+' the shared library is recompiled only
/// if it does not exist yet or the macro file is newer than the shared
/// library.
///
/// Of course the + and ++ notation is supported in similar way for .x and .L.
///
/// Through the function TSystem::SetMakeSharedLib(), the user will be able to
/// indicate, with shell commands, how to build a shared library (a good
/// default will be provided). The most common change, namely where to find
/// header files, will be available through the function
/// TSystem::SetIncludePath().
/// A good default will be provided so that a typical user session should be at
/// most:
/// ~~~ {.cpp}
/// root[1] gSystem->SetIncludePath("-I$ROOTSYS/include
/// -I$HOME/mypackage/include");
/// root[2] .x myfunc.C++(10,20);
/// ~~~
/// The user may sometimes try to compile a script before it has loaded all the
/// needed shared libraries.  In this case we want to be helpful and output a
/// list of the unresolved symbols. So if the loading of the created shared
/// library fails, we will try to build a executable that contains the
/// script. The linker should then output a list of missing symbols.
///
/// To support this we provide a TSystem::SetMakeExe() function, that sets the
/// directive telling how to create an executable. The loader will need
/// to be informed of all the libraries available. The information about
/// the libraries that has been loaded by .L and TSystem::Load() is accessible
/// to the script compiler. However, the information about
/// the libraries that have been selected at link time by the application
/// builder (like the root libraries for root.exe) are not available and need
/// to be explicitly listed in fLinkedLibs (either by default or by a call to
/// TSystem::SetLinkedLibs()).
///
/// To simplify customization we could also add to the .rootrc support for the
/// variables
/// ~~~ {.cpp}
/// Unix.*.Root.IncludePath:     -I$ROOTSYS/include
/// WinNT.*.Root.IncludePath:    -I%ROOTSYS%/include
///
/// Unix.*.Root.LinkedLibs:      -L$ROOTSYS/lib -lBase ....
/// WinNT.*.Root.LinkedLibs:     %ROOTSYS%/lib/*.lib msvcrt.lib ....
/// ~~~
/// And also support for MakeSharedLibs() and MakeExe().
///
/// (the ... have to be replaced by the actual values and are here only to
/// shorten this comment).

int TSystem::CompileMacro(const char *filename, Option_t *opt,
                          const char *library_specified,
                          const char *build_dir,
                          UInt_t dirmode)
{
   static const char *version_var_prefix = "__ROOTBUILDVERSION=";

   // ======= Analyze the options
   Bool_t keep = kFALSE;
   Bool_t recompile = kFALSE;
   EAclicMode mode = fAclicMode;
   Bool_t loadLib = kTRUE;
   Bool_t withInfo = kTRUE;
   Bool_t verbose = kFALSE;
   Bool_t internalDebug = kFALSE;
   if (opt) {
      keep = (strchr(opt,'k')!=0);
      recompile = (strchr(opt,'f')!=0);
      if (strchr(opt,'O')!=0) {
         mode = kOpt;
      }
      if (strchr(opt,'g')!=0) {
         mode = kDebug;
      }
      if (strchr(opt,'c')!=0) {
         loadLib = kFALSE;
      }
      withInfo = strchr(opt, 's') == nullptr;
      verbose = strchr(opt, 'v') != nullptr;
      internalDebug = strchr(opt, 'd') != nullptr;
   }
   if (mode==kDefault) {
      TString rootbuild = ROOTBUILD;
      if (rootbuild.Index("debug",0,TString::kIgnoreCase)==kNPOS) {
         mode=kOpt;
      } else {
         mode=kDebug;
      }
   }
   UInt_t verboseLevel = verbose ? 7 : gDebug;
   Bool_t flatBuildDir = (fAclicProperties & kFlatBuildDir) || (strchr(opt,'-')!=0);

   // if non-zero, build_loc indicates where to build the shared library.
   TString build_loc = ExpandFileName(GetBuildDir());
   if (build_dir && strlen(build_dir)) build_loc = build_dir;
   if (build_loc == ".") {
      build_loc = WorkingDirectory();
   } else if (build_loc.Length() && (!IsAbsoluteFileName(build_loc)) ) {
      AssignAndDelete( build_loc , ConcatFileName( WorkingDirectory(), build_loc ) );
   }

   // Get the include directory list in the dir1:dir2:dir3 format
   // [Used for generating the .d file and to look for header files for
   // the linkdef file]
   TString incPath = GetIncludePath(); // of the form -Idir1  -Idir2 -Idir3
   incPath.Append(":").Prepend(" ");
   if (gEnv) {
      TString fromConfig = gEnv->GetValue("ACLiC.IncludePaths","");
      incPath.Append(fromConfig);
   }
   incPath.ReplaceAll(" -I",":");       // of form :dir1 :dir2:dir3
   auto posISysRoot = incPath.Index(" -isysroot \"");
   if (posISysRoot != kNPOS) {
      auto posISysRootEnd = incPath.Index('"', posISysRoot + 12);
      if (posISysRootEnd != kNPOS) {
         // NOTE: should probably just skip isysroot for dependency analysis.
         // (And will, in the future - once we rely on compiler-generated .d files.)
         incPath.Insert(posISysRootEnd - 1, "/usr/include/");
         incPath.Replace(posISysRoot, 12, ":\"");
      }
   }
   while ( incPath.Index(" :") != -1 ) {
      incPath.ReplaceAll(" :",":");
   }
   incPath.Prepend(":.:");
   incPath.Prepend(WorkingDirectory());

   // ======= Get the right file names for the dictionary and the shared library
   TString expFileName(filename);
   ExpandPathName( expFileName );
   expFileName = gSystem->UnixPathName(expFileName);
   TString library = expFileName;
   if (! IsAbsoluteFileName(library) )
   {
      const char *whichlibrary = Which(incPath,library);
      if (whichlibrary) {
         library = whichlibrary;
         delete [] whichlibrary;
      } else {
         ::Error("ACLiC","The file %s can not be found in the include path: %s",filename,incPath.Data());
         return kFALSE;
      }
   } else {
      if (gSystem->AccessPathName(library)) {
         ::Error("ACLiC","The file %s can not be found.",filename);
         return kFALSE;
      }
   }
   { // Remove multiple '/' characters, rootcling treats them as comments.
      Ssiz_t pos = 0;
      while ((pos = library.Index("//", 2, pos, TString::kExact)) != kNPOS) {
         library.Remove(pos, 1);
      }
   }
   library = gSystem->UnixPathName(library);
   TString filename_fullpath = library;

   TString file_dirname = GetDirName( filename_fullpath );
   // For some probably good reason, DirName on Windows returns the 'name' of
   // the directory, omitting the drive letter (even if there was one). In
   // consequence the result is not usable as a 'root directory', we need to
   // add the drive letter if there was one..
   if (library.Length()>1 && isalpha(library[0]) && library[1]==':') {
      file_dirname.Prepend(library(0,2));
   }
   TString file_location( file_dirname  ); // Location of the script.
   incPath.Prepend( file_location + ":" );

   Ssiz_t dot_pos = library.Last('.');
   TString extension, libname_noext = library;
   if (dot_pos >= 0) {
      libname_noext.Remove(dot_pos);
      extension = library(dot_pos+1, library.Length()-dot_pos-1);
   }

   // Extension of shared library is platform dependent!!
   TString suffix = TString("_") + extension + "." + fSoExt;
   if (dot_pos >= 0)
      library.Replace( dot_pos, library.Length()-dot_pos, suffix);
   else
      library.Append(suffix);

   TString libname ( BaseName( libname_noext ) );
   libname.Append("_").Append(extension);

   if (library_specified && strlen(library_specified) ) {
      // Use the specified name instead of the default
      libname = BaseName( library_specified );
      library = library_specified;
      ExpandPathName( library );
      if (! IsAbsoluteFileName(library) ) {
         AssignAndDelete( library , ConcatFileName( WorkingDirectory(), library ) );
      }
      library = TString(library) + "." + fSoExt;
   }
   library = gSystem->UnixPathName(library);

   TString libname_ext ( libname );
   libname_ext +=  "." + fSoExt;

   TString lib_dirname = GetDirName( library );
   // For some probably good reason, DirName on Windows returns the 'name' of
   // the directory, omitting the drive letter (even if there was one). In
   // consequence the result is not useable as a 'root directory', we need to
   // add the drive letter if there was one..
   if (library.Length()>1 && isalpha(library[0]) && library[1]==':') {
      lib_dirname.Prepend(library(0,2));
   }
   // Strip potential, somewhat redundant '/.' from the pathname ...
   if ( strncmp( &(lib_dirname[lib_dirname.Length()-2]), "/.", 2) == 0 ) {
      lib_dirname.Remove(lib_dirname.Length()-2);
   }
   if ( strncmp( &(lib_dirname[lib_dirname.Length()-2]), "\\.", 2) == 0 ) {
      lib_dirname.Remove(lib_dirname.Length()-2);
   }
   TString lib_location( lib_dirname );
   Bool_t mkdirFailed = kFALSE;

   if (build_loc.Length()==0) {
      build_loc = lib_location;
   } else {
      // Removes an existing disk specification from the names
      TRegexp disk_finder ("[A-z]:");
      Int_t pos = library.Index( disk_finder );
      if (pos==0) library.Remove(pos,3);
      pos = lib_location.Index( disk_finder );
      if (pos==0) lib_location.Remove(pos,3);

      if (flatBuildDir) {
         AssignAndDelete( library, ConcatFileName( build_loc, libname_ext) );
      } else {
         AssignAndDelete( library, ConcatFileName( build_loc, library) );
      }

      Bool_t canWriteBuild_loc = !gSystem->AccessPathName(build_loc,kWritePermission);
      TString build_loc_store( build_loc );
      if (!flatBuildDir) {
         AssignAndDelete( build_loc, ConcatFileName( build_loc, lib_location) );
      }

      if (gSystem->AccessPathName(build_loc,kFileExists)) {
         mkdirFailed = (0 != mkdir(build_loc, true));
         if (mkdirFailed && !canWriteBuild_loc) {
            // The mkdir failed __and__ we can not write to the target directory,
            // let make sure the error message will be about the target directory
            build_loc = build_loc_store;
            mkdirFailed = kFALSE;
         } else if (!mkdirFailed && dirmode!=0) {
            Chmod(build_loc,dirmode);
         }
      }
   }
   library = gSystem->UnixPathName(library);

   // ======= Check if the library need to loaded or compiled
   if (!gInterpreter->IsLibraryLoaded(library) && gInterpreter->IsLoaded(expFileName)) {
      // the script has already been loaded in interpreted mode
      // Let's warn the user and unload it.

      if (withInfo) {
         ::Info("ACLiC","script has already been loaded in interpreted mode");
         ::Info("ACLiC","unloading %s and compiling it", filename);
      }

      if ( gInterpreter->UnloadFile( expFileName ) != 0 ) {
         // We can not unload it.
         return kFALSE;
      }
   }

   // Calculate the -I lines
   TString includes = GetIncludePath();
   includes.Prepend(' ');

   {
      // I need to replace the -Isomerelativepath by -I../ (or -I..\ on NT)
      TRegexp rel_inc(" -I[^\"/\\$%-][^:-]+");
      Int_t len,pos;
      pos = rel_inc.Index(includes,&len);
      while( len != 0 ) {
         TString sub = includes(pos,len);
         sub.Remove(0,3); // Remove ' -I'
         AssignAndDelete( sub, ConcatFileName( WorkingDirectory(), sub ) );
         sub.Prepend(" -I\"");
         sub.Chop(); // Remove trailing space (i.e between the -Is ...
         sub.Append("\" ");
         includes.Replace(pos,len,sub);
         pos = rel_inc.Index(includes,&len);
      }
   }
   {
       // I need to replace the -I"somerelativepath" by -I"$cwd/ (or -I"$cwd\ on NT)
      TRegexp rel_inc(" -I\"[^/\\$%-][^:-]+");
      Int_t len,pos;
      pos = rel_inc.Index(includes,&len);
      while( len != 0 ) {
         TString sub = includes(pos,len);
         sub.Remove(0,4); // Remove ' -I"'
         AssignAndDelete( sub, ConcatFileName( WorkingDirectory(), sub ) );
         sub.Prepend(" -I\"");
         includes.Replace(pos,len,sub);
         pos = rel_inc.Index(includes,&len);
      }
   }
   //includes += " -I\"" + build_loc;
   //includes += "\" -I\"";
   //includes += WorkingDirectory();
//   if (includes[includes.Length()-1] == '\\') {
//      // The current directory is (most likely) the root of a windows drive and
//      // has a trailing \ which would espace the quote if left by itself.
//      includes += '\\';
//   }
//   includes += "\"";
   if (gEnv) {
      TString fromConfig = gEnv->GetValue("ACLiC.IncludePaths","");
      includes.Append(" ").Append(fromConfig).Append(" ");
   }

   // Extract the -D for the dependency generation.
   TString defines = " ";
   {
      TString cmd = GetMakeSharedLib();
      TRegexp rel_def("-D[^\\s\\t\\n\\r]*");
      Int_t len,pos;
      pos = rel_def.Index(cmd,&len);
      while( len != 0 ) {
         defines += cmd(pos,len);
         defines += " ";
         pos = rel_def.Index(cmd,&len,pos+1);
      }

   }

   TString emergency_loc;
   {
      UserGroup_t *ug = gSystem->GetUserInfo(gSystem->GetUid());
      if (ug) {
         AssignAndDelete( emergency_loc, ConcatFileName( TempDirectory(), ug->fUser ) );
         delete ug;
      } else {
         emergency_loc = TempDirectory();
      }
   }

   Bool_t canWrite = !gSystem->AccessPathName(build_loc,kWritePermission);

   Bool_t modified = kFALSE;

   // Generate the dependency filename
   TString depdir = build_loc;
   TString depfilename;
   AssignAndDelete( depfilename, ConcatFileName(depdir, BaseName(libname_noext)) );
   depfilename += "_" + extension + ".d";

   if ( !recompile ) {

      Long_t lib_time, file_time;

      if ((gSystem->GetPathInfo( library, 0, (Long_t*)0, 0, &lib_time ) != 0) ||
          (gSystem->GetPathInfo( expFileName, 0, (Long_t*)0, 0, &file_time ) == 0 &&
          (lib_time < file_time))) {

         // the library does not exist or is older than the script.
         recompile = kTRUE;
         modified  = kTRUE;

      } else {

         if ( gSystem->GetPathInfo( depfilename, 0,(Long_t*) 0, 0, &file_time ) != 0 ) {
            if (!canWrite) {
               depdir = emergency_loc;
               AssignAndDelete( depfilename, ConcatFileName(depdir, BaseName(libname_noext)) );
               depfilename += "_" + extension + ".d";
            }
            R__WriteDependencyFile(build_loc, depfilename, filename_fullpath, library, libname, extension, version_var_prefix, includes, defines, incPath);
         }
      }

      if (!modified) {

         // We need to check the dependencies
         FILE * depfile = fopen(depfilename.Data(),"r");
         if (depfile==0) {
            // there is no accessible dependency file, let's assume the library has been
            // modified
            modified = kTRUE;
            recompile = kTRUE;

         } else {

            TString version_var = libname + version_var_prefix;

            Int_t sz = 256;
            char *line = new char[sz];
            line[0] = 0;

            int c;
            Int_t current = 0;
            Int_t nested = 0;
            Bool_t hasversion = false;

            while ((c = fgetc(depfile)) != EOF) {
               if (c=='#') {
                  // skip comment
                  while ((c = fgetc(depfile)) != EOF) {
                     if (c=='\n') {
                        break;
                     }
                  }
                  continue;
               }
               if (current && line[current-1]=='=' && strncmp(version_var.Data(),line,current)==0) {

                  // The next word will be the version number.
                  hasversion = kTRUE;
                  line[0] = 0;
                  current = 0;
               } else if (isspace(c) && !nested) {
                  if (current) {
                     if (line[current-1]!=':') {
                        // ignore target
                        line[current] = 0;

                        Long_t filetime;
                        if (hasversion) {
                           modified |= strcmp(ROOT_RELEASE,line)!=0;
                           hasversion = kFALSE;
                        } else if ( gSystem->GetPathInfo( line, 0, (Long_t*)0, 0, &filetime ) == 0 ) {
                           modified |= ( lib_time <= filetime );
                        }
                     }
                  }
                  current = 0;
                  line[0] = 0;
               } else {
                  if (current==sz-1) {
                     sz = 2*sz;
                     char *newline = new char[sz];
                     memcpy(newline,line, current);
                     delete [] line;
                     line = newline;
                  }
                  if (c=='"') nested = !nested;
                  else {
                     line[current] = c;
                     current++;
                  }
               }
            }
            delete [] line;
            fclose(depfile);
            recompile = modified;

         }

      }
   }

   if ( gInterpreter->IsLibraryLoaded(library)
        || strlen(GetLibraries(library,"D",kFALSE)) != 0 ) {
      // The library has already been built and loaded.

      Bool_t reload = kFALSE;
      TNamed *libinfo = (TNamed*)fCompiled->FindObject(library);
      if (libinfo) {
         Long_t load_time = libinfo->GetUniqueID();
         Long_t lib_time;
         if ( gSystem->GetPathInfo( library, 0, (Long_t*)0, 0, &lib_time ) == 0
              && (lib_time>load_time)) {
            reload = kTRUE;
         }
      }

      if ( !recompile && reload ) {

         if (withInfo) {
            ::Info("ACLiC","%s has been modified and will be reloaded",
                   libname.Data());
         }
         if ( gInterpreter->UnloadFile( library.Data() ) != 0 ) {
            // The library is being used. We can not unload it.
            return kFALSE;
         }
         if (libinfo) {
            fCompiled->Remove(libinfo);
            delete libinfo;
            libinfo = 0;
         }
         TNamed *k = new TNamed(library,library);
         Long_t lib_time;
         gSystem->GetPathInfo( library, 0, (Long_t*)0, 0, &lib_time );
         k->SetUniqueID(lib_time);
         if (!keep) k->SetBit(kMustCleanup);
         fCompiled->Add(k);

         return !gSystem->Load(library);
      }

      if (withInfo) {
         ::Info("ACLiC","%s script has already been compiled and loaded",
                modified ? "modified" : "unmodified");
      }

      if ( !recompile ) {
         return kTRUE;
      } else {
         if (withInfo) {
            ::Info("ACLiC","it will be regenerated and reloaded!");
         }
         if ( gInterpreter->UnloadFile( library.Data() ) != 0 ) {
            // The library is being used. We can not unload it.
            return kFALSE;
         }
         if (libinfo) {
            fCompiled->Remove(libinfo);
            delete libinfo;
            libinfo = 0;
         }
         Unlink(library);
      }

   }

   TString libmapfilename;
   AssignAndDelete( libmapfilename, ConcatFileName( build_loc, libname ) );
   libmapfilename += ".rootmap";
#if (defined(R__MACOSX) && !defined(MAC_OS_X_VERSION_10_5)) || defined(R__WIN32)
   Bool_t produceRootmap = kTRUE;
#else
   Bool_t produceRootmap = kFALSE;
#endif
   Bool_t linkDepLibraries = !produceRootmap;
   if (gEnv) {
#if (defined(R__MACOSX) && !defined(MAC_OS_X_VERSION_10_5))
      Int_t linkLibs = gEnv->GetValue("ACLiC.LinkLibs",2);
#elif defined(R__WIN32)
      Int_t linkLibs = gEnv->GetValue("ACLiC.LinkLibs",3);
#else
      Int_t linkLibs = gEnv->GetValue("ACLiC.LinkLibs",1);
#endif
      produceRootmap = linkLibs & 0x2;
      linkDepLibraries = linkLibs & 0x1;
   }

   // FIXME: Triggers clang false positive warning -Wunused-lambda-capture.
   /*constexpr const*/ bool useCxxModules =
#ifdef R__USE_CXXMODULES
    true;
#else
    false;
#endif

   // FIXME: Switch to generic polymorphic when we make c++14 default.
   auto ForeachSharedLibDep = [](const char *lib, std::function<bool(const char *)> f) {
      using namespace std;
      string deps = gInterpreter->GetSharedLibDeps(lib, /*tryDyld*/ true);
      istringstream iss(deps);
      vector<string> libs{istream_iterator<std::string>{iss}, istream_iterator<string>{}};
      // Skip the first element: it is a relative path to `lib`.
      for (auto I = libs.begin() + 1, E = libs.end(); I != E; ++I)
         if (!f(I->c_str()))
            break;
   };
   auto LoadLibrary = [useCxxModules, produceRootmap, ForeachSharedLibDep](const TString &lib) {
      // We have no rootmap files or modules to construct `-l` flags enabling
      // explicit linking. We have to resolve the dependencies by ourselves
      // taking the job of the dyld.
      // FIXME: This is a rare case where we have rootcling running with
      // modules disabled. Remove this code once we fully switch to modules,
      // or implement a special flag in rootcling which selective enables
      // modules for dependent libraries and does not produce a module for
      // the ACLiC library.
      if (useCxxModules && !produceRootmap) {
         std::function<bool(const char *)> LoadLibF = [](const char *dep) {
            return gInterpreter->Load(dep, /*skipReload*/ true) >= 0;
         };
         ForeachSharedLibDep(lib, LoadLibF);
      }
      return !gSystem->Load(lib);
   };

   if (!recompile) {
      // The library already exist, let's just load it.
      if (loadLib) {
         TNamed *k = new TNamed(library,library);
         Long_t lib_time;
         gSystem->GetPathInfo( library, 0, (Long_t*)0, 0, &lib_time );
         k->SetUniqueID(lib_time);
         if (!keep) k->SetBit(kMustCleanup);
         fCompiled->Add(k);

         gInterpreter->GetSharedLibDeps(library);

         return LoadLibrary(library);
      }
      else return kTRUE;
   }

   if (!canWrite && recompile) {

      if (mkdirFailed) {
         ::Warning("ACLiC","Could not create the directory: %s",
                build_loc.Data());
      } else {
         ::Warning("ACLiC","%s is not writable!",
                   build_loc.Data());
      }
      if (emergency_loc == build_dir ) {
         ::Error("ACLiC","%s is the last resort location (i.e. temp location)",build_loc.Data());
         return kFALSE;
      }
      ::Warning("ACLiC","Output will be written to %s",
                emergency_loc.Data());
      return CompileMacro(expFileName, opt, library_specified, emergency_loc, dirmode);
   }

   if (withInfo) {
      Info("ACLiC","creating shared library %s",library.Data());
   }

   R__WriteDependencyFile(build_loc, depfilename, filename_fullpath, library, libname, extension, version_var_prefix, includes, defines, incPath);

   // ======= Select the dictionary name
   TString dict = libname + "_ACLiC_dict";

   // the file name end up in the file produced
   // by rootcling as a variable name so all character need to be valid!
   static const int maxforbidden = 27;
   static const char *forbidden_chars[maxforbidden] =
      { "+","-","*","/","&","%","|","^",">","<",
        "=","~",".","(",")","[","]","!",",","$",
        " ",":","'","#","@","\\","\"" };
   for( int ic = 0; ic < maxforbidden; ic++ ) {
      dict.ReplaceAll( forbidden_chars[ic],"_" );
   }
   if ( dict.Last('.')!=dict.Length()-1 ) dict.Append(".");
   AssignAndDelete( dict, ConcatFileName( build_loc, dict ) );
   TString dicth = dict;
   TString dictObj = dict;
   dict += "cxx"; //no need to keep the extension of the original file, any extension will do
   dicth += "h";
   dictObj += fObjExt;

   // ======= Generate a linkdef file

   TString linkdef;
   AssignAndDelete( linkdef, ConcatFileName( build_loc, libname ) );
   linkdef += "_ACLiC_linkdef.h";
   std::ofstream linkdefFile( linkdef, std::ios::out );
   linkdefFile << "// File Automatically generated by the ROOT Script Compiler "
               << std::endl;
   linkdefFile << std::endl;
   linkdefFile << "#ifdef __CINT__" << std::endl;
   linkdefFile << std::endl;
   linkdefFile << "#pragma link C++ nestedclasses;" << std::endl;
   linkdefFile << "#pragma link C++ nestedtypedefs;" << std::endl;
   linkdefFile << std::endl;

   // We want to look for a header file that has the same name as the macro

   const char * extensions[] = { ".h", ".hh", ".hpp", ".hxx",  ".hPP", ".hXX" };

   int i;
   for (i = 0; i < 6; i++ ) {
      char * name;
      TString extra_linkdef = BaseName( libname_noext );
      extra_linkdef.Append(GetLinkdefSuffix());
      extra_linkdef.Append(extensions[i]);
      name = Which(incPath,extra_linkdef);
      if (name) {
         if (verboseLevel>4 && withInfo) {
            Info("ACLiC","including extra linkdef file: %s",name);
         }
         linkdefFile << "#include \"" << name << "\"" << std::endl;
         delete [] name;
      }
   }

   if (verboseLevel>5 && withInfo) {
      Info("ACLiC","looking for header in: %s",incPath.Data());
   }
   for (i = 0; i < 6; i++ ) {
      char * name;
      TString lookup = BaseName( libname_noext );
      lookup.Append(extensions[i]);
      name = Which(incPath,lookup);
      if (name) {
         linkdefFile << "#pragma link C++ defined_in "<<gSystem->UnixPathName(name)<<";"<< std::endl;
         delete [] name;
      }
   }
   linkdefFile << "#pragma link C++ defined_in \""<<filename_fullpath << "\";" << std::endl;
   linkdefFile << std::endl;
   linkdefFile << "#endif" << std::endl;
   linkdefFile.close();
   // ======= Generate the list of rootmap files to be looked at

   TString mapfile;
   AssignAndDelete( mapfile, ConcatFileName( build_loc, libname ) );
   mapfile += "_ACLiC_map";
   TString mapfilein = mapfile + ".in";
   TString mapfileout = mapfile + ".out";

   Bool_t needLoadMap = kFALSE;
   if (!useCxxModules) {
      if (gInterpreter->GetSharedLibDeps(library) != nullptr) {
         gInterpreter->UnloadLibraryMap(libname);
         needLoadMap = kTRUE;
      }
   }

   std::ofstream mapfileStream( mapfilein, std::ios::out );
   {
      TString name = ".rootmap";
      TString sname = "system.rootmap";
      TString file;
      AssignAndDelete(file, ConcatFileName(TROOT::GetEtcDir(), sname) );
      if (gSystem->AccessPathName(file)) {
         // for backward compatibility check also $ROOTSYS/system<name> if
         // $ROOTSYS/etc/system<name> does not exist
         AssignAndDelete(file, ConcatFileName(TROOT::GetRootSys(), sname));
         if (gSystem->AccessPathName(file)) {
            // for backward compatibility check also $ROOTSYS/<name> if
            // $ROOTSYS/system<name> does not exist
            AssignAndDelete(file, ConcatFileName(TROOT::GetRootSys(), name));
         }
      }
      mapfileStream << file << std::endl;
      AssignAndDelete(file, ConcatFileName(gSystem->HomeDirectory(), name) );
      mapfileStream << file << std::endl;
      mapfileStream << name << std::endl;
      if (gInterpreter->GetRootMapFiles()) {
         for (i = 0; i < gInterpreter->GetRootMapFiles()->GetEntriesFast(); i++) {
            mapfileStream << ((TNamed*)gInterpreter->GetRootMapFiles()->At(i))->GetTitle() << std::endl;
         }
      }
   }
   mapfileStream.close();

   // ======= Generate the rootcling command line
   TString rcling = "rootcling";
   PrependPathName(TROOT::GetBinDir(), rcling);
   rcling += " -v0 \"--lib-list-prefix=";
   rcling += mapfile;
   rcling += "\" -f \"";
   rcling.Append(dict).Append("\" ");

   if (produceRootmap && !useCxxModules) {
      rcling += " -rml " + libname + " -rmf \"" + libmapfilename + "\" ";
      rcling.Append("-DR__ACLIC_ROOTMAP ");
   }
   rcling.Append(GetIncludePath()).Append(" -D__ACLIC__ ");
   if (gEnv) {
      TString fromConfig = gEnv->GetValue("ACLiC.IncludePaths","");
      rcling.Append(fromConfig);
      TString extraFlags = gEnv->GetValue("ACLiC.ExtraRootclingFlags","");
      if (!extraFlags.IsNull()) {
        extraFlags.Prepend(" ");
        extraFlags.Append(" ");
        rcling.Append(extraFlags);
      }
   }

   // Create a modulemap
   // FIXME: Merge the modulemap generation from cmake and here in rootcling.
   if (useCxxModules && produceRootmap) {
      rcling += " -cxxmodule ";
      // TString moduleMapFileName = file_dirname + "/" + libname + ".modulemap";
      TString moduleName = libname + "_ACLiC_dict";
      if (moduleName.BeginsWith("lib"))
          moduleName = moduleName.Remove(0, 3);
      TString moduleMapName = moduleName + ".modulemap";
      TString moduleMapFullPath = build_loc + "/" + moduleMapName;
      // A modulemap may exist from previous runs, overwrite it.
      if (verboseLevel > 3 && !AccessPathName(moduleMapFullPath))
         ::Info("ACLiC", "File %s already exists!", moduleMapFullPath.Data());

      std::string curDir = ROOT::FoundationUtils::GetCurrentDir();
      std::string relative_path = ROOT::FoundationUtils::MakePathRelative(filename_fullpath.Data(), curDir);
      std::ofstream moduleMapFile(moduleMapFullPath, std::ios::out);
      moduleMapFile << "module \"" << moduleName << "\" {" << std::endl;
      moduleMapFile << "  header \"" << relative_path << "\"" << std::endl;
      moduleMapFile << "  export *" << std::endl;
      moduleMapFile << "  link \"" << libname_ext << "\"" << std::endl;
      moduleMapFile << "}" << std::endl;
      moduleMapFile.close();
      gInterpreter->RegisterPrebuiltModulePath(build_loc.Data(), moduleMapName.Data());
      rcling.Append(" \"-moduleMapFile=" + moduleMapFullPath + "\" ");
   }

   rcling.Append(" \"").Append(filename_fullpath).Append("\" ");
   rcling.Append("\"").Append(linkdef).Append("\"");

   // ======= Run rootcling
   if (withInfo) {
      if (verboseLevel>3) {
         ::Info("ACLiC","creating the dictionary files");
         if (verboseLevel>4)  ::Info("ACLiC", "%s", rcling.Data());
      }
   }

   ///\returns true on success.
   auto ExecAndReport = [](TString cmd) -> bool {
      Int_t result = gSystem->Exec(cmd);
      if (result) {
         if (result == 139)
            ::Error("ACLiC", "Executing '%s' failed with a core dump!", cmd.Data());
         else
            ::Error("ACLiC", "Executing '%s' failed!", cmd.Data());
      }
      return !result;
   };

   Bool_t result = ExecAndReport(rcling);
   TString depLibraries;

   // ======= Load the library the script might depend on
   if (result) {
      TString linkedlibs = GetLibraries("", "S");
      TString libtoload;
      TString all_libtoload;
      std::ifstream liblist(mapfileout);

      while ( liblist >> libtoload ) {
         // Load the needed library except for the library we are currently building!
         if (libtoload == "#") {
            // The comment terminates the list of libraries.
            std::string toskipcomment;
            std::getline(liblist,toskipcomment);
            break;
         }
         if (libtoload != library && libtoload != libname && libtoload != libname_ext) {
            if (produceRootmap) {
               if (loadLib || linkDepLibraries /* For GetLibraries to Work */) {
                  result = gROOT->LoadClass("", libtoload) >= 0;
                  if (!result) {
                     // We failed to load one of the dependency.
                     break;
                  }
               }
               if (!linkedlibs.Contains(libtoload)) {
                  all_libtoload.Append(" ").Append(libtoload);
                  depLibraries.Append(" ");
                  depLibraries.Append(GetLibraries(libtoload,"DSL",kFALSE));
                  depLibraries = depLibraries.Strip(); // Remove any trailing spaces.
               }
            } else {
               gROOT->LoadClass("", libtoload);
            }
         }
         unsigned char c = liblist.peek();
         if (c=='\n' || c=='\r') {
            // Consume the character
            liblist.get();
            break;
         }
      }

//      depLibraries = all_libtoload;
//      depLibraries.ReplaceAll(" lib"," -l");
//      depLibraries.ReplaceAll(TString::Format(".%s",fSoExt.Data()),"");
   }

   // ======= Calculate the libraries for linking:
   TString linkLibraries;
   /*
     this is intentionally disabled until it can become useful
     if (gEnv) {
        linkLibraries =  gEnv->GetValue("ACLiC.Libraries","");
        linkLibraries.Prepend(" ");
     }
   */
   TString linkLibrariesNoQuotes(GetLibraries("","SDL"));
   // We need to enclose the single paths in quotes to account for paths with spaces
   TString librariesWithQuotes;
   TString singleLibrary;
   Bool_t collectingSingleLibraryNameTokens = kFALSE;
   for (auto tokenObj : *linkLibrariesNoQuotes.Tokenize(" ")) {
      singleLibrary = ((TObjString*)tokenObj)->GetString();
      if (singleLibrary[0]=='-' || !AccessPathName(singleLibrary)) {
         if (collectingSingleLibraryNameTokens) {
            librariesWithQuotes.Chop();
            librariesWithQuotes += "\" \"" + singleLibrary + "\"";
            collectingSingleLibraryNameTokens = kFALSE;
         } else {
            librariesWithQuotes += " \"" + singleLibrary + "\"";
         }
      } else {
         if (collectingSingleLibraryNameTokens) {
            librariesWithQuotes += singleLibrary + " ";
         } else {
            collectingSingleLibraryNameTokens = kTRUE;
            librariesWithQuotes += " \"" + singleLibrary + " ";
         }
      }
   }

#ifdef _MSC_VER
   linkLibraries.Prepend(linkLibrariesNoQuotes);
#else
   linkLibraries.Prepend(librariesWithQuotes);
#endif

   // ======= Generate the build command lines
   TString cmd = fMakeSharedLib;
   // we do not add filename because it is already included via the dictionary(in dicth) !
   // dict.Append(" ").Append(filename);
   cmd.ReplaceAll("$SourceFiles","-D__ACLIC__ \"$SourceFiles\"");
   cmd.ReplaceAll("$SourceFiles",dict);
   cmd.ReplaceAll("$ObjectFiles","\"$ObjectFiles\"");
   cmd.ReplaceAll("$ObjectFiles",dictObj);
   cmd.ReplaceAll("$IncludePath",includes);
   cmd.ReplaceAll("$SharedLib","\"$SharedLib\"");
   cmd.ReplaceAll("$SharedLib",library);
   if (linkDepLibraries) {
      if (produceRootmap) {
         cmd.ReplaceAll("$DepLibs",depLibraries);
      } else {
         cmd.ReplaceAll("$DepLibs",linkLibraries);
      }
   }
   cmd.ReplaceAll("$LinkedLibs",linkLibraries);
   cmd.ReplaceAll("$LibName",libname);
   cmd.ReplaceAll("\"$BuildDir","$BuildDir");
   cmd.ReplaceAll("$BuildDir","\"$BuildDir\"");
   cmd.ReplaceAll("$BuildDir",build_loc);
   if (mode==kDebug) {
      cmd.ReplaceAll("$Opt",fFlagsDebug);
   } else {
      cmd.ReplaceAll("$Opt",fFlagsOpt);
   }
#ifdef WIN32
   R__FixLink(cmd);
   cmd.ReplaceAll("-std=", "-std:");
#endif

   TString testcmd = fMakeExe;
   TString fakeMain;
   AssignAndDelete( fakeMain, ConcatFileName( build_loc, libname ) );
   fakeMain += "_ACLiC_main";
   fakeMain += extension;
   std::ofstream fakeMainFile( fakeMain, std::ios::out );
   fakeMainFile << "// File Automatically generated by the ROOT Script Compiler "
                << std::endl;
   fakeMainFile << "int main(char*argc,char**argvv) {};" << std::endl;
   fakeMainFile.close();
   // We could append this fake main routine to the compilation line.
   // But in this case compiler may output the name of the dictionary file
   // and of the fakeMain file while it compiles it. (this would be useless
   // confusing output).
   // We could also the fake main routine to the end of the dictionary file
   // however compilation would fail if a main is already there
   // (like stress.cxx)
   // dict.Append(" ").Append(fakeMain);
   TString exec;
   AssignAndDelete( exec, ConcatFileName( build_loc, libname ) );
   exec += "_ACLiC_exec";
   testcmd.ReplaceAll("$SourceFiles","-D__ACLIC__ \"$SourceFiles\"");
   testcmd.ReplaceAll("$SourceFiles",dict);
   testcmd.ReplaceAll("$ObjectFiles","\"$ObjectFiles\"");
   testcmd.ReplaceAll("$ObjectFiles",dictObj);
   testcmd.ReplaceAll("$IncludePath",includes);
   testcmd.ReplaceAll("$ExeName",exec);
   testcmd.ReplaceAll("$LinkedLibs",linkLibraries);
   testcmd.ReplaceAll("$BuildDir",build_loc);
   if (mode==kDebug)
      testcmd.ReplaceAll("$Opt",fFlagsDebug);
   else
      testcmd.ReplaceAll("$Opt",fFlagsOpt);

#ifdef WIN32
   R__FixLink(testcmd);
   testcmd.ReplaceAll("-std=", "-std:");
#endif

   // ======= Build the library
   if (result) {
      if (verboseLevel>3 && withInfo) {
         ::Info("ACLiC","compiling the dictionary and script files");
         if (verboseLevel>4)  ::Info("ACLiC", "%s", cmd.Data());
      }
      Int_t success = ExecAndReport(cmd);
      if (!success) {
         if (produceRootmap) {
            gSystem->Unlink(libmapfilename);
         }
      }
      result = success;
   }

   if ( result ) {
      if (linkDepLibraries) {
         // We may have unresolved symbols. Use dyld to resolve the dependent
         // libraries and relink.
         // FIXME: We will likely have duplicated libraries as we are appending
         // FIXME: This likely makes rootcling --lib-list-prefix redundant.
         TString depLibsFullPaths;
         std::function<bool(const char *)> CollectF = [&depLibsFullPaths](const char *dep) {
            TString LibFullPath(dep);
            if (!gSystem->FindDynamicLibrary(LibFullPath, /*quiet=*/true)) {
               ::Error("TSystem::CompileMacro", "Cannot find library '%s'", dep);
               return false; // abort
            }
            depLibsFullPaths += " " + LibFullPath;
            return true;
         };
         ForeachSharedLibDep(library, CollectF);

         TString relink_cmd = cmd.Strip(TString::kTrailing, ';');
         relink_cmd += depLibsFullPaths;
         result = ExecAndReport(relink_cmd);
      }

      TNamed *k = new TNamed(library,library);
      Long_t lib_time;
      gSystem->GetPathInfo( library, 0, (Long_t*)0, 0, &lib_time );
      k->SetUniqueID(lib_time);
      if (!keep) k->SetBit(kMustCleanup);
      fCompiled->Add(k);

      if (needLoadMap) {
          gInterpreter->LoadLibraryMap(libmapfilename);
      }
      if (verboseLevel>3 && withInfo)  ::Info("ACLiC","loading the shared library");
      if (loadLib)
         result = LoadLibrary(library);
      else
         result = kTRUE;

      if ( !result ) {
         if (verboseLevel>3 && withInfo) {
            ::Info("ACLiC","testing for missing symbols:");
            if (verboseLevel>4)  ::Info("ACLiC", "%s", testcmd.Data());
         }
         gSystem->Exec(testcmd);
         gSystem->Unlink( exec );
      }

   };

   if (verboseLevel<=5 && !internalDebug) {
      gSystem->Unlink( dict );
      gSystem->Unlink( dicth );
      gSystem->Unlink( dictObj );
      gSystem->Unlink( linkdef );
      gSystem->Unlink( mapfilein );
      gSystem->Unlink( mapfileout );
      gSystem->Unlink( fakeMain );
      gSystem->Unlink( exec );
   }
   if (verboseLevel>6) {
      rcling.Prepend("echo ");
      cmd.Prepend("echo \" ").Append(" \" ");
      testcmd.Prepend("echo \" ").Append(" \" ");
      gSystem->Exec(rcling);
      gSystem->Exec( cmd );
      gSystem->Exec(testcmd);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the ACLiC properties field.   See EAclicProperties for details
/// on the semantic of each bit.

Int_t TSystem::GetAclicProperties() const
{
   return fAclicProperties;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the build architecture.

const char *TSystem::GetBuildArch() const
{
   return fBuildArch;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the build compiler

const char *TSystem::GetBuildCompiler() const
{
   return fBuildCompiler;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the build compiler version

const char *TSystem::GetBuildCompilerVersion() const
{
   return fBuildCompilerVersion;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the build compiler version identifier string

const char *TSystem::GetBuildCompilerVersionStr() const
{
   return fBuildCompilerVersionStr;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the build node name.

const char *TSystem::GetBuildNode() const
{
   return fBuildNode;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the path of the build directory.

const char *TSystem::GetBuildDir() const
{
   if (fBuildDir.Length()==0) {
      if (!gEnv) return "";
      const_cast<TSystem*>(this)->fBuildDir = gEnv->GetValue("ACLiC.BuildDir","");
   }
   return fBuildDir;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the debug flags.

const char *TSystem::GetFlagsDebug() const
{
   return fFlagsDebug;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the optimization flags.

const char *TSystem::GetFlagsOpt() const
{
   return fFlagsOpt;
}

////////////////////////////////////////////////////////////////////////////////
/// AclicMode indicates whether the library should be built in
/// debug mode or optimized.  The values are:
///  - TSystem::kDefault : compile the same as the current ROOT
///  - TSystem::kDebug : compiled in debug mode
///  - TSystem::kOpt : optimized the library

TSystem::EAclicMode TSystem::GetAclicMode() const
{
   return fAclicMode;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the command line use to make a shared library.
/// See TSystem::CompileMacro for more details.

const char *TSystem::GetMakeSharedLib() const
{
   return fMakeSharedLib;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the command line use to make an executable.
/// See TSystem::CompileMacro for more details.

const char *TSystem::GetMakeExe() const
{
   return fMakeExe;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the list of include path.

const char *TSystem::GetIncludePath()
{
   fListPaths = fIncludePath;
#ifndef _MSC_VER
   // FIXME: This is a temporary fix for the following error with ACLiC
   // (and this is apparently not needed anyway):
   // 48: input_line_12:8:38: error: use of undeclared identifier 'IC'
   // 48: "C:/Users/bellenot/build/debug/etc" -IC:/Users/bellenot/build/debug/etc//cling -IC:/Users/bellenot/build/debug/include"",
   // 48:                                      ^
   // 48: Error in <ACLiC>: Dictionary generation failed!
   fListPaths.Append(" ").Append(gInterpreter->GetIncludePath());
#endif
   return fListPaths;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the list of library linked to this executable.
/// See TSystem::CompileMacro for more details.

const char *TSystem::GetLinkedLibs() const
{
   return fLinkedLibs;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the linkdef suffix chosen by the user for ACLiC.
/// See TSystem::CompileMacro for more details.

const char *TSystem::GetLinkdefSuffix() const
{
   if (fLinkdefSuffix.Length()==0) {
      if (!gEnv) return "_linkdef";
      const_cast<TSystem*>(this)->fLinkdefSuffix = gEnv->GetValue("ACLiC.Linkdef","_linkdef");
   }
   return fLinkdefSuffix;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the shared library extension.

const char *TSystem::GetSoExt() const
{
   return fSoExt;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the object file extension.

const char *TSystem::GetObjExt() const
{
   return fObjExt;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the location where ACLiC will create libraries and use as
/// a scratch area.
///
/// If 'isflat' is false, then the libraries are actually stored in
/// sub-directories of 'build_dir' including the full pathname of the
/// script.  If the script is location at /full/path/name/macro.C
/// the library will be located at 'build_dir+/full/path/name/macro_C.so'
/// If 'isflat' is true, then no subdirectory is created and the library
/// is created directly in the directory 'build_dir'.  Note that in this
/// mode there is a risk than 2 script of the same in different source
/// directory will over-write each other.

void TSystem::SetBuildDir(const char *build_dir, Bool_t isflat)
{
   fBuildDir = build_dir;
   if (isflat)
      fAclicProperties |= (kFlatBuildDir & kBitMask);
   else
      fAclicProperties &= ~(kFlatBuildDir & kBitMask);
}

////////////////////////////////////////////////////////////////////////////////
/// FlagsDebug should contain the options to pass to the C++ compiler
/// in order to compile the library in debug mode.

void TSystem::SetFlagsDebug(const char *flags)
{
   fFlagsDebug = flags;
}

////////////////////////////////////////////////////////////////////////////////
/// FlagsOpt should contain the options to pass to the C++ compiler
/// in order to compile the library in optimized mode.

void TSystem::SetFlagsOpt(const char *flags)
{
   fFlagsOpt = flags;
}

////////////////////////////////////////////////////////////////////////////////
/// AclicMode indicates whether the library should be built in
/// debug mode or optimized.  The values are:
///  - TSystem::kDefault : compile the same as the current ROOT
///  - TSystem::kDebug : compiled in debug mode
///  - TSystem::kOpt : optimized the library

void TSystem::SetAclicMode(EAclicMode mode)
{
   fAclicMode = mode;
}

////////////////////////////////////////////////////////////////////////////////
/// Directives has the same syntax as the argument of SetMakeSharedLib but is
/// used to create an executable. This creation is used as a means to output
/// a list of unresolved symbols, when loading a shared library has failed.
/// The required variable is $ExeName rather than $SharedLib, e.g.:
/// ~~~ {.cpp}
/// gSystem->SetMakeExe(
/// "g++ -Wall -fPIC $IncludePath $SourceFiles
///  -o $ExeName $LinkedLibs -L/usr/X11R6/lib -lX11 -lm -ldl -rdynamic");
/// ~~~

void TSystem::SetMakeExe(const char *directives)
{
   fMakeExe = directives;
   // NOTE: add verification that the directives has the required variables
}

////////////////////////////////////////////////////////////////////////////////
/// Directives should contain the description on how to compile and link a
/// shared lib. This description can be any valid shell command, including
/// the use of ';' to separate several instructions. However, shell specific
/// construct should be avoided. In particular this description can contain
/// environment variables, like $ROOTSYS (or %ROOTSYS% on windows).
/// ~~~ {.cpp}
/// Five special variables will be expanded before execution:
///   Variable name       Expands to
///   -------------       ----------
///   $SourceFiles        Name of source files to be compiled
///   $SharedLib          Name of the shared library being created
///   $LibName            Name of shared library without extension
///   $BuildDir           Directory where the files will be created
///   $IncludePath        value of fIncludePath
///   $LinkedLibs         value of fLinkedLibs
///   $DepLibs            libraries on which this library depends on
///   $ObjectFiles        Name of source files to be compiler with
///                       their extension changed to .o or .obj
///   $Opt                location of the optimization/debug options
///                       set fFlagsDebug and fFlagsOpt
/// ~~~
/// e.g.:
/// ~~~ {.cpp}
/// gSystem->SetMakeSharedLib(
/// "KCC -n32 --strict $IncludePath -K0 \$Opt $SourceFile
///  --no_exceptions --signed_chars --display_error_number
///  --diag_suppress 68 -o $SharedLib");
///
/// gSystem->setMakeSharedLib(
/// "Cxx $IncludePath -c $SourceFile;
///  ld  -L/usr/lib/cmplrs/cxx -rpath /usr/lib/cmplrs/cxx -expect_unresolved
///  \$Opt -shared /usr/lib/cmplrs/cc/crt0.o /usr/lib/cmplrs/cxx/_main.o
///  -o $SharedLib $ObjectFile -lcxxstd -lcxx -lexc -lots -lc"
///
/// gSystem->SetMakeSharedLib(
/// "$HOME/mygcc/bin/g++ \$Opt -Wall -fPIC $IncludePath $SourceFile
///  -shared -o $SharedLib");
///
/// gSystem->SetMakeSharedLib(
/// "cl -DWIN32  -D_WIN32 -D_MT -D_DLL -MD /O2 /G5 /MD -DWIN32
///  -D_WINDOWS $IncludePath $SourceFile
///  /link -PDB:NONE /NODEFAULTLIB /INCREMENTAL:NO /RELEASE /NOLOGO
///  $LinkedLibs -entry:_DllMainCRTStartup@12 -dll /out:$SharedLib")
/// ~~~

void TSystem::SetMakeSharedLib(const char *directives)
{
   fMakeSharedLib = directives;
   // NOTE: add verification that the directives has the required variables
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Add a directory to the already set include path.
/// \param[in] includePath The path to the directory.
/// \note This interface is mostly relevant for ACLiC and it does *not* inform
///       gInterpreter for this include path. If the TInterpreter needs to know
///       about the include path please use \c gInterpreter->AddIncludePath .
/// \warning The path should start with the \c -I prefix, i.e.
///          <tt>gSystem->AddIncludePath("-I /path/to/my/includes")</tt>.
void TSystem::AddIncludePath(const char *includePath)
{
   if (includePath) {
      fIncludePath += " ";
      fIncludePath += includePath;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add linkedLib to already set linked libs.

void TSystem::AddLinkedLibs(const char *linkedLib)
{
   if (linkedLib) {
      fLinkedLibs += " ";
      fLinkedLibs += linkedLib;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// IncludePath should contain the list of compiler flags to indicate where
/// to find user defined header files. It is used to expand $IncludePath in
/// the directives given to SetMakeSharedLib() and SetMakeExe(), e.g.:
/// ~~~ {.cpp}
///    gSystem->SetInclude("-I$ROOTSYS/include -Imydirectory/include");
/// ~~~
/// the default value of IncludePath on Unix is:
/// ~~~ {.cpp}
///    "-I$ROOTSYS/include "
/// ~~~
/// and on Windows:
/// ~~~ {.cpp}
///    "/I%ROOTSYS%/include "
/// ~~~

void TSystem::SetIncludePath(const char *includePath)
{
   fIncludePath = includePath;
}

////////////////////////////////////////////////////////////////////////////////
/// LinkedLibs should contain the library directory and list of libraries
/// needed to recreate the current executable. It is used to expand $LinkedLibs
/// in the directives given to SetMakeSharedLib() and SetMakeExe()
/// The default value on Unix is: `root-config --glibs`

void  TSystem::SetLinkedLibs(const char *linkedLibs)
{
   fLinkedLibs = linkedLibs;
}

////////////////////////////////////////////////////////////////////////////////
/// The 'suffix' will be appended to the name of a script loaded by ACLiC
/// and used to locate any eventual additional linkdef information that
/// ACLiC should used to produce the dictionary.
///
/// So by default, when doing .L MyScript.cxx, ACLiC will look
/// for a file name MyScript_linkdef and having one of the .h (.hpp,
/// etc.) extensions.  If such a file exist, it will be added to
/// the end of the linkdef file used to created the ACLiC dictionary.
/// This effectively enable the full customization of the creation
/// of the dictionary.  It should be noted that the file is intended
/// as a linkdef `fragment`, so usually you would not list the
/// typical:
/// ~~~ {.cpp}
///   #pragma link off ....
/// ~~~

void  TSystem::SetLinkdefSuffix(const char *suffix)
{
   fLinkdefSuffix = suffix;
}


////////////////////////////////////////////////////////////////////////////////
/// Set shared library extension, should be either .so, .sl, .a, .dll, etc.

void TSystem::SetSoExt(const char *SoExt)
{
   fSoExt = SoExt;
}

////////////////////////////////////////////////////////////////////////////////
/// Set object files extension, should be either .o, .obj, etc.

void TSystem::SetObjExt(const char *ObjExt)
{
   fObjExt = ObjExt;
}

////////////////////////////////////////////////////////////////////////////////
/// This method split a filename of the form:
/// ~~~ {.cpp}
///   [path/]macro.C[+|++[k|f|g|O|c|s|d|v|-]][(args)].
/// ~~~
/// It stores the ACliC mode [+|++[options]] in 'mode',
/// the arguments (including parenthesis) in arg
/// and the I/O indirection in io

TString TSystem::SplitAclicMode(const char *filename, TString &aclicMode,
                                TString &arguments, TString &io) const
{
   char *fname = Strip(filename);
   TString filenameCopy = fname;
   filenameCopy = filenameCopy.Strip();

   if (filenameCopy.EndsWith(";")) {
     filenameCopy.Remove(filenameCopy.Length() - 1);
     filenameCopy = filenameCopy.Strip();
   }
   if (filenameCopy.EndsWith(")")) {
      Ssiz_t posArgEnd = filenameCopy.Length() - 1;
      // There is an argument; find its start!
      int parenNestCount = 1;
      bool inString = false;
      Ssiz_t posArgBegin = posArgEnd - 1;
      for (; parenNestCount && posArgBegin >= 0; --posArgBegin) {
         // Escaped if the previous character is a `\` - but not if it
         // itself is preceded by a `\`!
         if (posArgBegin > 0 && filenameCopy[posArgBegin] == '\\' &&
             (posArgBegin == 1 || filenameCopy[posArgBegin - 1] != '\\')) {
            // skip escape.
            --posArgBegin;
            continue;
         }
         switch (filenameCopy[posArgBegin]) {
         case ')':
            if (!inString)
               ++parenNestCount;
            break;
         case '(':
            if (!inString)
               --parenNestCount;
            break;
         case '"': inString = !inString; break;
         }
      }
      if (parenNestCount || inString) {
         Error("SplitAclicMode", "Cannot parse argument in %s", filename);
      } else {
         arguments = filenameCopy(posArgBegin + 1, posArgEnd - 1);
         fname[posArgBegin + 1] = 0;
      }
   }

   // strip off I/O redirect tokens from filename
   {
      char *s2   = nullptr;
      char *s3;
      s2 = strstr(fname, ">>");
      if (!s2) s2 = strstr(fname, "2>");
      if (!s2) s2 = strchr(fname, '>');
      s3 = strchr(fname, '<');
      if (s2 && s3) s2 = s2<s3 ? s2 : s3;
      if (s3 && !s2) s2 = s3;
      if (s2==fname) {
         io = fname;
         aclicMode = "";
         arguments = "";
         delete []fname;
         return "";
      } else if (s2) {
         s2--;
         while (s2 && *s2 == ' ') s2--;
         s2++;
         io = s2; // ssave = *s2;
         *s2 = 0;
      } else
         io = "";
   }

   // remove the possible ACLiC + or ++ and g or O etc
   aclicMode.Clear();
   int len = strlen(fname);
   TString mode;
   while (len > 1) {
      if (strchr("kfgOcsdv-", fname[len - 1])) {
         mode += fname[len - 1];
         --len;
      } else {
         break;
      }
   }
   Bool_t compile = len && fname[len - 1] == '+';
   Bool_t remove  = compile && len > 1 && fname[len - 2] == '+';
   if (compile) {
      if (mode.Length()) {
         fname[len] = 0;
      }
      if (remove) {
         fname[strlen(fname)-2] = 0;
         aclicMode = "++";
      } else {
         fname[strlen(fname)-1] = 0;
         aclicMode = "+";
      }
      if (mode.Length())
         aclicMode += mode;
   }

   TString resFilename = fname;

   delete []fname;
   return resFilename;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the shared libs produced by the CompileMacro() function.

void TSystem::CleanCompiledMacros()
{
   TIter next(fCompiled);
   TNamed *lib;
   while ((lib = (TNamed*)next())) {
      if (lib->TestBit(kMustCleanup)) Unlink(lib->GetTitle());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Register version of plugin library.

TVersionCheck::TVersionCheck(int versionCode)
{
   if (versionCode != TROOT::RootVersionCode() && gLibraryVersion)
      gLibraryVersion[gLibraryVersionIdx] = versionCode;
}
