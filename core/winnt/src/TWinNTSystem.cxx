// @(#)root/winnt:$Id: db9b3139b1551a1b4e31a17f57866a276d5cd419 $
// Author: Fons Rademakers   15/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// TWinNTSystem                                                                 //
//                                                                              //
// Class providing an interface to the Windows NT/Windows 95 Operating Systems. //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////


#ifdef HAVE_CONFIG
#include "config.h"
#endif

#include "Windows4Root.h"
#include "ROOT/FoundationUtils.hxx"
#include "TWinNTSystem.h"
#include "TROOT.h"
#include "TError.h"
#include "TOrdCollection.h"
#include "TRegexp.h"
#include "TException.h"
#include "TEnv.h"
#include "TApplication.h"
#include "TWin32SplashThread.h"
#include "Win32Constants.h"
#include "TInterpreter.h"
#include "TVirtualX.h"
#include "TUrl.h"
#include "ThreadLocalStorage.h"
#include "snprintf.h"
#include "strlcpy.h"

#include <sys/utime.h>
#include <sys/timeb.h>
#include <process.h>
#include <io.h>
#include <direct.h>
#include <ctype.h>
#include <float.h>
#include <sys/stat.h>
#include <signal.h>
#include <stdio.h>
#include <errno.h>
#include <lm.h>
#include <dbghelp.h>
#include <Tlhelp32.h>
#include <sstream>
#include <iostream>
#include <list>
#include <shlobj.h>
#include <conio.h>

#if defined (_MSC_VER) && (_MSC_VER >= 1400)
   #include <intrin.h>
#elif defined (_M_IX86)
   static void __cpuid(int* cpuid_data, int info_type)
   {
      __asm {
         push ebx
         push edi
         mov edi, cpuid_data
         mov eax, info_type
         cpuid
         mov [edi], eax
         mov [edi + 4], ebx
         mov [edi + 8], ecx
         mov [edi + 12], edx
         pop edi
         pop ebx
      }
   }
   __int64 __rdtsc()
   {
      LARGE_INTEGER li;
      __asm {
         rdtsc
         mov li.LowPart, eax
         mov li.HighPart, edx
      }
      return li.QuadPart;
   }
#else
   static void __cpuid(int* cpuid_data, int) {
      cpuid_data[0] = 0x00000000;
      cpuid_data[1] = 0x00000000;
      cpuid_data[2] = 0x00000000;
      cpuid_data[3] = 0x00000000;
   }
   __int64 __rdtsc() { return (__int64)0; }
#endif

extern "C" {
   extern void Gl_setwidth(int width);
   void *_ReturnAddress(void);
}

//////////////////// Windows TFdSet ////////////////////////////////////////////////
class TFdSet {
private:
   fd_set *fds_bits; // file descriptors (according MSDN maximum is 64)
public:
   TFdSet() { fds_bits = new fd_set; fds_bits->fd_count = 0; }
   virtual ~TFdSet() { delete fds_bits; }
   void  Copy(TFdSet &fd) const { memcpy((void*)fd.fds_bits, fds_bits, sizeof(fd_set)); }
   TFdSet(const TFdSet& fd) { fd.Copy(*this); }
   TFdSet& operator=(const TFdSet& fd)  { fd.Copy(*this); return *this; }
   void  Zero() { fds_bits->fd_count = 0; }
   void  Set(Int_t fd)
   {
      if (fds_bits->fd_count < FD_SETSIZE-1) // protect out of bound access (64)
         fds_bits->fd_array[fds_bits->fd_count++] = (SOCKET)fd;
      else
         ::SysError("TFdSet::Set", "fd_count will exeed FD_SETSIZE");
   }
   void  Clr(Int_t fd)
   {
      int i;
      for (i=0; i<fds_bits->fd_count; i++) {
         if (fds_bits->fd_array[i]==(SOCKET)fd) {
            while (i<fds_bits->fd_count-1) {
               fds_bits->fd_array[i] = fds_bits->fd_array[i+1];
               i++;
            }
            fds_bits->fd_count--;
            break;
         }
      }
   }
   Int_t IsSet(Int_t fd) { return __WSAFDIsSet((SOCKET)fd, fds_bits); }
   Int_t *GetBits() { return fds_bits && fds_bits->fd_count ? (Int_t*)fds_bits : 0; }
   UInt_t GetCount() { return (UInt_t)fds_bits->fd_count; }
   Int_t GetFd(Int_t i) { return i<fds_bits->fd_count ? fds_bits->fd_array[i] : 0; }
};

namespace {
   const char *kProtocolName   = "tcp";
   typedef void (*SigHandler_t)(ESignals);
   static TWinNTSystem::ThreadMsgFunc_t gGUIThreadMsgFunc = 0;      // GUI thread message handler func

   static HANDLE gGlobalEvent;
   static HANDLE gTimerThreadHandle;
   typedef NET_API_STATUS (WINAPI *pfn1)(LPVOID);
   typedef NET_API_STATUS (WINAPI *pfn2)(LPCWSTR, LPCWSTR, DWORD, LPBYTE*);
   typedef NET_API_STATUS (WINAPI *pfn3)(LPCWSTR, LPCWSTR, DWORD, LPBYTE*,
                                       DWORD, LPDWORD, LPDWORD, PDWORD);
   typedef NET_API_STATUS (WINAPI *pfn4)(LPCWSTR, DWORD, LPBYTE*, DWORD, LPDWORD,
                                       LPDWORD, PDWORD);
   static pfn1 p2NetApiBufferFree;
   static pfn2 p2NetUserGetInfo;
   static pfn3 p2NetLocalGroupGetMembers;
   static pfn4 p2NetLocalGroupEnum;

   static struct signal_map {
      int code;
      SigHandler_t handler;
      const char *signame;
   } signal_map[kMAXSIGNALS] = {   // the order of the signals should be identical
      -1 /*SIGBUS*/,   0, "bus error",    // to the one in SysEvtHandler.h
      SIGSEGV,  0, "segmentation violation",
      -1 /*SIGSYS*/,   0, "bad argument to system call",
      -1 /*SIGPIPE*/,  0, "write on a pipe with no one to read it",
      SIGILL,   0, "illegal instruction",
      SIGABRT,  0, "abort",
      -1 /*SIGQUIT*/,  0, "quit",
      SIGINT,   0, "interrupt",
      -1 /*SIGWINCH*/, 0, "window size change",
      -1 /*SIGALRM*/,  0, "alarm clock",
      -1 /*SIGCHLD*/,  0, "death of a child",
      -1 /*SIGURG*/,   0, "urgent data arrived on an I/O channel",
      SIGFPE,   0, "floating point exception",
      SIGTERM,  0, "termination signal",
      -1 /*SIGUSR1*/,  0, "user-defined signal 1",
      -1 /*SIGUSR2*/,  0, "user-defined signal 2"
   };

   ////// static functions providing interface to raw WinNT ////////////////////

   //---- RPC -------------------------------------------------------------------
   //*-* Error codes set by the Windows Sockets implementation are not made available
   //*-* via the errno variable. Additionally, for the getXbyY class of functions,
   //*-* error codes are NOT made available via the h_errno variable. Instead, error
   //*-* codes are accessed by using the WSAGetLastError . This function is provided
   //*-* in Windows Sockets as a precursor (and eventually an alias) for the Win32
   //*-* function GetLastError. This is intended to provide a reliable way for a thread
   //*-* in a multithreaded process to obtain per-thread error information.

   /////////////////////////////////////////////////////////////////////////////
   /// Receive exactly length bytes into buffer. Returns number of bytes
   /// received. Returns -1 in case of error, -2 in case of MSG_OOB
   /// and errno == EWOULDBLOCK, -3 in case of MSG_OOB and errno == EINVAL
   /// and -4 in case of kNonBlock and errno == EWOULDBLOCK.
   /// Returns -5 if pipe broken or reset by peer (EPIPE || ECONNRESET).

   static int WinNTRecv(int socket, void *buffer, int length, int flag)
   {
      if (socket == -1) return -1;
      SOCKET sock = socket;

      int once = 0;
      if (flag == -1) {
         flag = 0;
         once = 1;
      }
      if (flag == MSG_PEEK) {
         once = 1;
      }

      int nrecv, n;
      char *buf = (char *)buffer;

      for (n = 0; n < length; n += nrecv) {
         if ((nrecv = ::recv(sock, buf+n, length-n, flag)) <= 0) {
            if (nrecv == 0) {
               break;        // EOF
            }
            if (flag == MSG_OOB) {
               if (::WSAGetLastError() == WSAEWOULDBLOCK) {
                  return -2;
               } else if (::WSAGetLastError() == WSAEINVAL) {
                  return -3;
               }
            }
            if (::WSAGetLastError() == WSAEWOULDBLOCK) {
               return -4;
            } else {
               if (::WSAGetLastError() != WSAEINTR)
                  ::SysError("TWinNTSystem::WinNTRecv", "recv");
               if (::WSAGetLastError() == EPIPE ||
                  ::WSAGetLastError() == WSAECONNRESET)
                  return -5;
               else
                  return -1;
            }
         }
         if (once) {
            return nrecv;
         }
      }
      return n;
   }

   /////////////////////////////////////////////////////////////////////////////
   /// Send exactly length bytes from buffer. Returns -1 in case of error,
   /// otherwise number of sent bytes. Returns -4 in case of kNoBlock and
   /// errno == EWOULDBLOCK. Returns -5 if pipe broken or reset by peer
   /// (EPIPE || ECONNRESET).

   static int WinNTSend(int socket, const void *buffer, int length, int flag)
   {
      if (socket < 0) return -1;
      SOCKET sock = socket;

      int once = 0;
      if (flag == -1) {
         flag = 0;
         once = 1;
      }

      int nsent, n;
      const char *buf = (const char *)buffer;

      for (n = 0; n < length; n += nsent) {
         if ((nsent = ::send(sock, buf+n, length-n, flag)) <= 0) {
            if (nsent == 0) {
               break;
            }
            if (::WSAGetLastError() == WSAEWOULDBLOCK) {
               return -4;
            } else {
               if (::WSAGetLastError() != WSAEINTR)
                  ::SysError("TWinNTSystem::WinNTSend", "send");
               if (::WSAGetLastError() == EPIPE ||
                  ::WSAGetLastError() == WSAECONNRESET)
                  return -5;
               else
                  return -1;
            }
         }
         if (once) {
            return nsent;
         }
      }
      return n;
   }

   /////////////////////////////////////////////////////////////////////////////
   /// Wait for events on the file descriptors specified in the readready and
   /// writeready masks or for timeout (in milliseconds) to occur.

   static int WinNTSelect(TFdSet *readready, TFdSet *writeready, Long_t timeout)
   {
      int retcode;
      fd_set* rbits = readready ? (fd_set*)readready->GetBits() : 0;
      fd_set* wbits = writeready ? (fd_set*)writeready->GetBits() : 0;

      if (timeout >= 0) {
         timeval tv;
         tv.tv_sec  = timeout / 1000;
         tv.tv_usec = (timeout % 1000) * 1000;

         retcode = ::select(0, rbits, wbits, 0, &tv);
      } else {
         retcode = ::select(0, rbits, wbits, 0, 0);
      }

      if (retcode == SOCKET_ERROR) {
         int errcode = ::WSAGetLastError();

         // if file descriptor is not a socket, assume it is the pipe used
         // by TXSocket
         if (errcode == WSAENOTSOCK) {
            struct __stat64 buf;
            int result = _fstat64( readready->GetFd(0), &buf );
            if ( result == 0 ) {
               if (buf.st_size > 0)
                  return 1;
            }
            // yield execution to another thread that is ready to run
            // if no other thread is ready, sleep 1 ms before to return
            if (gGlobalEvent) {
               ::WaitForSingleObject(gGlobalEvent, 1);
               ::ResetEvent(gGlobalEvent);
            }
            return 0;
         }

         if ( errcode == WSAEINTR) {
            TSystem::ResetErrno();  // errno is not self reseting
            return -2;
         }
         if (errcode == EBADF) {
            return -3;
         }
         return -1;
      }
      return retcode;
   }

   /////////////////////////////////////////////////////////////////////////////
   /// Get shared library search path.

   static const char *DynamicPath(const char *newpath = 0, Bool_t reset = kFALSE)
   {
      static TString dynpath;

      if (reset || newpath) {
         dynpath = "";
      }
      if (newpath) {

         dynpath = newpath;

      } else if (dynpath == "") {
         TString rdynpath = gEnv ? gEnv->GetValue("Root.DynamicPath", (char*)0) : "";
         rdynpath.ReplaceAll("; ", ";");  // in case DynamicPath was extended
         if (rdynpath == "") {
            rdynpath = ".;"; rdynpath += TROOT::GetBinDir();
         }
         TString path = gSystem->Getenv("PATH");
         if (path == "")
            dynpath = rdynpath;
         else {
            dynpath = path; dynpath += ";"; dynpath += rdynpath;
         }

      }

      if (!dynpath.Contains(TROOT::GetLibDir())) {
         dynpath += ";"; dynpath += TROOT::GetLibDir();
      }

      return dynpath;
   }

   /////////////////////////////////////////////////////////////////////////////
   /// Call the signal handler associated with the signal.

   static void sighandler(int sig)
   {
      for (int i = 0; i < kMAXSIGNALS; i++) {
         if (signal_map[i].code == sig) {
            (*signal_map[i].handler)((ESignals)i);
            return;
         }
      }
   }

   /////////////////////////////////////////////////////////////////////////////
   /// Set a signal handler for a signal.

   static void WinNTSignal(ESignals sig, SigHandler_t handler)
   {
      signal_map[sig].handler = handler;
      if (signal_map[sig].code != -1)
         (SigHandler_t)signal(signal_map[sig].code, sighandler);
   }

   /////////////////////////////////////////////////////////////////////////////
   /// Return the signal name associated with a signal.

   static const char *WinNTSigname(ESignals sig)
   {
      return signal_map[sig].signame;
   }

   /////////////////////////////////////////////////////////////////////////////
   /// WinNT signal handler.

   static BOOL ConsoleSigHandler(DWORD sig)
   {
      switch (sig) {
         case CTRL_C_EVENT:
            if (gSystem) {
               ((TWinNTSystem*)gSystem)->DispatchSignals(kSigInterrupt);
            }
            else {
               Break("TInterruptHandler::Notify", "keyboard interrupt");
               if (TROOT::Initialized()) {
                  gInterpreter->RewindDictionary();
               }
            }
            return kTRUE;
         case CTRL_BREAK_EVENT:
         case CTRL_LOGOFF_EVENT:
         case CTRL_SHUTDOWN_EVENT:
         case CTRL_CLOSE_EVENT:
         default:
            printf("\n *** Break *** keyboard interrupt - ROOT is terminated\n");
            gSystem->Exit(-1);
            return kTRUE;
      }
   }

   static CONTEXT *fgXcptContext = 0;
   /////////////////////////////////////////////////////////////////////////////

   static void SigHandler(ESignals sig)
   {
      if (gSystem)
         ((TWinNTSystem*)gSystem)->DispatchSignals(sig);
   }

   /////////////////////////////////////////////////////////////////////////////
   /// Function that's called when an unhandled exception occurs.
   /// Produces a stack trace, and lets the system deal with it
   /// as if it was an unhandled excecption (usually ::abort)

   LONG WINAPI ExceptionFilter(LPEXCEPTION_POINTERS pXcp)
   {
      fgXcptContext = pXcp->ContextRecord;
      gSystem->StackTrace();
      return EXCEPTION_CONTINUE_SEARCH;
   }


#pragma intrinsic(_ReturnAddress)
#pragma auto_inline(off)
   DWORD_PTR GetProgramCounter()
   {
      // Returns the current program counter.
      return (DWORD_PTR)_ReturnAddress();
   }
#pragma auto_inline(on)

   /////////////////////////////////////////////////////////////////////////////
   /// Message processing loop for the TGWin32 related GUI
   /// thread for processing windows messages (aka Main/Server thread).
   /// We need to start the thread outside the TGWin32 / GUI related
   /// dll, because starting threads at DLL init time does not work.
   /// Instead, we start an ideling thread at binary startup, and only
   /// call the "real" message processing function
   /// TGWin32::GUIThreadMessageFunc() once gVirtualX comes up.

   static DWORD WINAPI GUIThreadMessageProcessingLoop(void *p)
   {
      MSG msg;

      // force to create message queue
      ::PeekMessage(&msg, NULL, WM_USER, WM_USER, PM_NOREMOVE);

      Int_t erret = 0;
      Bool_t endLoop = kFALSE;
      while (!endLoop) {
         if (gGlobalEvent) ::SetEvent(gGlobalEvent);
         erret = ::GetMessage(&msg, NULL, NULL, NULL);
         if (erret <= 0) endLoop = kTRUE;
         if (gGUIThreadMsgFunc)
            endLoop = (*gGUIThreadMsgFunc)(&msg);
      }

      gVirtualX->CloseDisplay();

      // exit thread
      if (erret == -1) {
         erret = ::GetLastError();
         Error("MsgLoop", "Error in GetMessage");
         ::ExitThread(-1);
      } else {
         ::ExitThread(0);
      }
      return 0;
   }

   //=========================================================================
   // Load IMAGEHLP.DLL and get the address of functions in it that we'll use
   // by Microsoft, from http://www.microsoft.com/msj/0597/hoodtextfigs.htm#fig1
   //=========================================================================
   // Make typedefs for some IMAGEHLP.DLL functions so that we can use them
   // with GetProcAddress
   typedef BOOL (__stdcall *SYMINITIALIZEPROC)( HANDLE, LPSTR, BOOL );
   typedef BOOL (__stdcall *SYMCLEANUPPROC)( HANDLE );
   typedef BOOL (__stdcall *STACKWALK64PROC)
               ( DWORD, HANDLE, HANDLE, LPSTACKFRAME64, LPVOID,
               PREAD_PROCESS_MEMORY_ROUTINE,PFUNCTION_TABLE_ACCESS_ROUTINE,
               PGET_MODULE_BASE_ROUTINE, PTRANSLATE_ADDRESS_ROUTINE );
   typedef LPVOID (__stdcall *SYMFUNCTIONTABLEACCESS64PROC)( HANDLE, DWORD64 );
   typedef DWORD (__stdcall *SYMGETMODULEBASE64PROC)( HANDLE, DWORD64 );
   typedef BOOL (__stdcall *SYMGETMODULEINFO64PROC)(HANDLE, DWORD64, PIMAGEHLP_MODULE64);
   typedef BOOL (__stdcall *SYMGETSYMFROMADDR64PROC)( HANDLE, DWORD64, PDWORD64, PIMAGEHLP_SYMBOL64);
   typedef BOOL (__stdcall *SYMGETLINEFROMADDR64PROC)(HANDLE, DWORD64, PDWORD, PIMAGEHLP_LINE64);
   typedef DWORD (__stdcall *UNDECORATESYMBOLNAMEPROC)(PCSTR, PSTR, DWORD, DWORD);


   static SYMINITIALIZEPROC _SymInitialize = 0;
   static SYMCLEANUPPROC _SymCleanup = 0;
   static STACKWALK64PROC _StackWalk64 = 0;
   static SYMFUNCTIONTABLEACCESS64PROC _SymFunctionTableAccess64 = 0;
   static SYMGETMODULEBASE64PROC _SymGetModuleBase64 = 0;
   static SYMGETMODULEINFO64PROC _SymGetModuleInfo64 = 0;
   static SYMGETSYMFROMADDR64PROC _SymGetSymFromAddr64 = 0;
   static SYMGETLINEFROMADDR64PROC _SymGetLineFromAddr64 = 0;
   static UNDECORATESYMBOLNAMEPROC _UnDecorateSymbolName = 0;

   BOOL InitImagehlpFunctions()
   {
      // Fetches function addresses from IMAGEHLP.DLL at run-time, so we
      // don't need to link against its import library. These functions
      // are used in StackTrace; if they cannot be found (e.g. because
      // IMAGEHLP.DLL doesn't exist or has the wrong version) we cannot
      // produce a stack trace.

      HMODULE hModImagehlp = LoadLibrary( "IMAGEHLP.DLL" );
      if (!hModImagehlp)
         return FALSE;

      _SymInitialize = (SYMINITIALIZEPROC) GetProcAddress( hModImagehlp, "SymInitialize" );
      if (!_SymInitialize)
         return FALSE;

      _SymCleanup = (SYMCLEANUPPROC) GetProcAddress( hModImagehlp, "SymCleanup" );
      if (!_SymCleanup)
         return FALSE;

      _StackWalk64 = (STACKWALK64PROC) GetProcAddress( hModImagehlp, "StackWalk64" );
      if (!_StackWalk64)
         return FALSE;

      _SymFunctionTableAccess64 = (SYMFUNCTIONTABLEACCESS64PROC) GetProcAddress(hModImagehlp, "SymFunctionTableAccess64" );
      if (!_SymFunctionTableAccess64)
         return FALSE;

      _SymGetModuleBase64=(SYMGETMODULEBASE64PROC)GetProcAddress(hModImagehlp, "SymGetModuleBase64");
      if (!_SymGetModuleBase64)
         return FALSE;

      _SymGetModuleInfo64=(SYMGETMODULEINFO64PROC)GetProcAddress(hModImagehlp, "SymGetModuleInfo64");
      if (!_SymGetModuleInfo64)
         return FALSE;

      _SymGetSymFromAddr64=(SYMGETSYMFROMADDR64PROC)GetProcAddress(hModImagehlp, "SymGetSymFromAddr64");
      if (!_SymGetSymFromAddr64)
         return FALSE;

      _SymGetLineFromAddr64=(SYMGETLINEFROMADDR64PROC)GetProcAddress(hModImagehlp, "SymGetLineFromAddr64");
      if (!_SymGetLineFromAddr64)
         return FALSE;

      _UnDecorateSymbolName=(UNDECORATESYMBOLNAMEPROC)GetProcAddress(hModImagehlp, "UnDecorateSymbolName");
      if (!_UnDecorateSymbolName)
         return FALSE;

      if (!_SymInitialize(GetCurrentProcess(), 0, TRUE ))
         return FALSE;

      return TRUE;
   }

   // stack trace helpers getModuleName, getFunctionName by
   /**************************************************************************
   * VRS - The Virtual Rendering System
   * Copyright (C) 2000-2004 Computer Graphics Systems Group at the
   * Hasso-Plattner-Institute (HPI), Potsdam, Germany.
   * This library is free software; you can redistribute it and/or modify it
   * under the terms of the GNU Lesser General Public License as published by
   * the Free Software Foundation; either version 2.1 of the License, or
   * (at your option) any later version.
   ***************************************************************************/
   std::string GetModuleName(DWORD64 address)
   {
      // Return the name of the module that contains the function at address.
      // Used by StackTrace.
      std::ostringstream out;
      HANDLE process = ::GetCurrentProcess();

      DWORD lineDisplacement = 0;
      IMAGEHLP_LINE64 line;
      ::ZeroMemory(&line, sizeof(line));
      line.SizeOfStruct = sizeof(line);
      if(_SymGetLineFromAddr64(process, address, &lineDisplacement, &line)) {
            out << line.FileName << "(" << line.LineNumber << "): ";
      } else {
            IMAGEHLP_MODULE64 module;
            ::ZeroMemory(&module, sizeof(module));
            module.SizeOfStruct = sizeof(module);
            if(_SymGetModuleInfo64(process, address, &module)) {
               out << module.ModuleName << "!";
            } else {
               out << "0x" << std::hex << address << std::dec << " ";
            }
      }

      return out.str();
   }

   std::string GetFunctionName(DWORD64 address)
   {
      // Return the name of the function at address.
      // Used by StackTrace.
      DWORD64 symbolDisplacement = 0;
      HANDLE process = ::GetCurrentProcess();

      const unsigned int SYMBOL_BUFFER_SIZE = 8192;
      char symbolBuffer[SYMBOL_BUFFER_SIZE];
      PIMAGEHLP_SYMBOL64 symbol = reinterpret_cast<PIMAGEHLP_SYMBOL64>(symbolBuffer);
      ::ZeroMemory(symbol, SYMBOL_BUFFER_SIZE);
      symbol->SizeOfStruct = SYMBOL_BUFFER_SIZE;
      symbol->MaxNameLength = SYMBOL_BUFFER_SIZE - sizeof(IMAGEHLP_SYMBOL64);

      if(_SymGetSymFromAddr64(process, address, &symbolDisplacement, symbol)) {
            // Make the symbol readable for humans
            const unsigned int NAME_SIZE = 8192;
            char name[NAME_SIZE];
            _UnDecorateSymbolName(
               symbol->Name,
               name,
               NAME_SIZE,
               UNDNAME_COMPLETE             |
               UNDNAME_NO_THISTYPE          |
               UNDNAME_NO_SPECIAL_SYMS      |
               UNDNAME_NO_MEMBER_TYPE       |
               UNDNAME_NO_MS_KEYWORDS       |
               UNDNAME_NO_ACCESS_SPECIFIERS
            );

            std::string result;
            result += name;
            result += "()";
            return result;
      } else {
            return "??";
      }
   }

   ////// Shortcuts helper functions IsShortcut and ResolveShortCut ///////////

   /////////////////////////////////////////////////////////////////////////////
   /// Validates if a file name has extension '.lnk'. Returns true if file
   /// name have extension same as Window's shortcut file (.lnk).

   static BOOL IsShortcut(const char *filename)
   {
      //File extension for the Window's shortcuts (.lnk)
      const char *extLnk = ".lnk";
      if (filename != NULL) {
         //Validate extension
         TString strfilename(filename);
         if (strfilename.EndsWith(extLnk))
            return TRUE;
      }
      return FALSE;
   }

   /////////////////////////////////////////////////////////////////////////////
   /// Resolve a ShellLink (i.e. c:\path\shortcut.lnk) to a real path.

   static BOOL ResolveShortCut(LPCSTR pszShortcutFile, char *pszPath, int maxbuf)
   {
      HRESULT hres;
      IShellLink* psl;
      char szGotPath[MAX_PATH];
      WIN32_FIND_DATA wfd;

      *pszPath = 0;   // assume failure

      // Make typedefs for some ole32.dll functions so that we can use them
      // with GetProcAddress
      typedef HRESULT (__stdcall *COINITIALIZEPROC)( LPVOID );
      static COINITIALIZEPROC _CoInitialize = 0;
      typedef void (__stdcall *COUNINITIALIZEPROC)( void );
      static COUNINITIALIZEPROC _CoUninitialize = 0;
      typedef HRESULT (__stdcall *COCREATEINSTANCEPROC)( REFCLSID, LPUNKNOWN,
                       DWORD, REFIID, LPVOID );
      static COCREATEINSTANCEPROC _CoCreateInstance = 0;

      HMODULE hModImagehlp = LoadLibrary( "ole32.dll" );
      if (!hModImagehlp)
         return FALSE;

      _CoInitialize = (COINITIALIZEPROC) GetProcAddress( hModImagehlp, "CoInitialize" );
      if (!_CoInitialize)
         return FALSE;
      _CoUninitialize = (COUNINITIALIZEPROC) GetProcAddress( hModImagehlp, "CoUninitialize");
      if (!_CoUninitialize)
         return FALSE;
      _CoCreateInstance = (COCREATEINSTANCEPROC) GetProcAddress( hModImagehlp, "CoCreateInstance" );
      if (!_CoCreateInstance)
         return FALSE;

      _CoInitialize(NULL);

      hres = _CoCreateInstance(CLSID_ShellLink, NULL, CLSCTX_INPROC_SERVER,
                               IID_IShellLink, (void **) &psl);
      if (SUCCEEDED(hres)) {
         IPersistFile* ppf;

         hres = psl->QueryInterface(IID_IPersistFile, (void **) &ppf);
         if (SUCCEEDED(hres)) {
            WCHAR wsz[MAX_PATH];
            MultiByteToWideChar(CP_ACP, 0, pszShortcutFile, -1, wsz, MAX_PATH);

            hres = ppf->Load(wsz, STGM_READ);
            if (SUCCEEDED(hres)) {
               hres = psl->Resolve(HWND_DESKTOP, SLR_ANY_MATCH | SLR_NO_UI | SLR_UPDATE);
               if (SUCCEEDED(hres)) {
                  strlcpy(szGotPath, pszShortcutFile,MAX_PATH);
                  hres = psl->GetPath(szGotPath, MAX_PATH, (WIN32_FIND_DATA *)&wfd,
                                      SLGP_UNCPRIORITY | SLGP_RAWPATH);
                  strlcpy(pszPath,szGotPath, maxbuf);
                  if (maxbuf) pszPath[maxbuf-1] = 0;
               }
            }
            ppf->Release();
         }
         psl->Release();
      }
      _CoUninitialize();

      return SUCCEEDED(hres);
   }

   void UpdateRegistry(TWinNTSystem* sys, char* buf /* size of buffer: MAX_MODULE_NAME32 + 1 */) {
      // register ROOT as the .root file handler:
      GetModuleFileName(0, buf, MAX_MODULE_NAME32 + 1);
      if (strcmp(sys->TWinNTSystem::BaseName(buf), "root.exe"))
         return;
      HKEY regCUS;
      if (!::RegOpenKeyEx(HKEY_CURRENT_USER, "Software", 0, KEY_READ, &regCUS) == ERROR_SUCCESS)
         return;
      HKEY regCUSC;
      if (!::RegOpenKeyEx(regCUS, "Classes", 0, KEY_READ, &regCUSC) == ERROR_SUCCESS) {
         ::RegCloseKey(regCUS);
         return;
      }

      HKEY regROOT;
      bool regROOTwrite = false;
      TString iconloc(buf);
      iconloc += ",-101";

      if (::RegOpenKeyEx(regCUSC, "ROOTDEV.ROOT", 0, KEY_READ, &regROOT) != ERROR_SUCCESS) {
         ::RegCloseKey(regCUSC);
         if (::RegOpenKeyEx(regCUS, "Classes", 0, KEY_READ | KEY_WRITE, &regCUSC) == ERROR_SUCCESS &&
            ::RegCreateKeyEx(regCUSC, "ROOTDEV.ROOT", 0, NULL, 0, KEY_READ | KEY_WRITE,
            NULL, &regROOT, NULL) == ERROR_SUCCESS) {
            regROOTwrite = true;
         }
      } else {
         HKEY regROOTIcon;
         if (::RegOpenKeyEx(regROOT, "DefaultIcon", 0, KEY_READ, &regROOTIcon) == ERROR_SUCCESS) {
            char bufIconLoc[1024];
            DWORD dwType;
            DWORD dwSize = sizeof(bufIconLoc);

            if (::RegQueryValueEx(regROOTIcon, NULL, NULL, &dwType, (BYTE*)bufIconLoc, &dwSize) == ERROR_SUCCESS)
               regROOTwrite = (iconloc != bufIconLoc);
            else
               regROOTwrite = true;
            ::RegCloseKey(regROOTIcon);
         } else
            regROOTwrite = true;
         if (regROOTwrite) {
            // re-open for writing
            ::RegCloseKey(regCUSC);
            ::RegCloseKey(regROOT);
            if (::RegOpenKeyEx(regCUS, "Classes", 0, KEY_READ | KEY_WRITE, &regCUSC) != ERROR_SUCCESS) {
               // error opening key for writing:
               regROOTwrite = false;
            } else {
               if (::RegOpenKeyEx(regCUSC, "ROOTDEV.ROOT", 0, KEY_WRITE, &regROOT) != ERROR_SUCCESS) {
                  // error opening key for writing:
                  regROOTwrite = false;
                  ::RegCloseKey(regCUSC);
               }
            }
         }
      }

      // determine the fileopen.C file path:
      TString fileopen = "fileopen.C";
      TString rootmacrodir = "macros";
      sys->PrependPathName(getenv("ROOTSYS"), rootmacrodir);
      sys->PrependPathName(rootmacrodir.Data(), fileopen);

      if (regROOTwrite) {
         // only write to registry if fileopen.C is readable
         regROOTwrite = (::_access(fileopen, kReadPermission) == 0);
      }

      if (!regROOTwrite) {
         ::RegCloseKey(regROOT);
         ::RegCloseKey(regCUSC);
         ::RegCloseKey(regCUS);
         return;
      }

      static const char apptitle[] = "ROOT data file";
      ::RegSetValueEx(regROOT, NULL, 0, REG_SZ, (BYTE*)apptitle, sizeof(apptitle));
      DWORD editflags = /*FTA_OpenIsSafe*/ 0x00010000; // trust downloaded files
      ::RegSetValueEx(regROOT, "EditFlags", 0, REG_DWORD, (BYTE*)&editflags, sizeof(editflags));

      HKEY regROOTIcon;
      if (::RegCreateKeyEx(regROOT, "DefaultIcon", 0, NULL, 0, KEY_READ | KEY_WRITE,
                           NULL, &regROOTIcon, NULL) == ERROR_SUCCESS) {
         TString iconloc(buf);
         iconloc += ",-101";
         ::RegSetValueEx(regROOTIcon, NULL, 0, REG_SZ, (BYTE*)iconloc.Data(), iconloc.Length() + 1);
         ::RegCloseKey(regROOTIcon);
      }

      // "open" verb
      HKEY regROOTshell;
      if (::RegCreateKeyEx(regROOT, "shell", 0, NULL, 0, KEY_READ | KEY_WRITE,
                           NULL, &regROOTshell, NULL) == ERROR_SUCCESS) {
         HKEY regShellOpen;
         if (::RegCreateKeyEx(regROOTshell, "open", 0, NULL, 0, KEY_READ | KEY_WRITE,
                              NULL, &regShellOpen, NULL) == ERROR_SUCCESS) {
            HKEY regShellOpenCmd;
            if (::RegCreateKeyEx(regShellOpen, "command", 0, NULL, 0, KEY_READ | KEY_WRITE,
                                 NULL, &regShellOpenCmd, NULL) == ERROR_SUCCESS) {
               TString cmd(buf);
               cmd += " -l \"%1\" \"";
               cmd += fileopen;
               cmd += "\"";
               ::RegSetValueEx(regShellOpenCmd, NULL, 0, REG_SZ, (BYTE*)cmd.Data(), cmd.Length() + 1);
               ::RegCloseKey(regShellOpenCmd);
            }
            ::RegCloseKey(regShellOpen);
         }
         ::RegCloseKey(regROOTshell);
      }
      ::RegCloseKey(regROOT);

      if (::RegCreateKeyEx(regCUSC, ".root", 0, NULL, 0, KEY_READ | KEY_WRITE,
                           NULL, &regROOT, NULL) == ERROR_SUCCESS) {
         static const char appname[] = "ROOTDEV.ROOT";
         ::RegSetValueEx(regROOT, NULL, 0, REG_SZ, (BYTE*)appname, sizeof(appname));
      }
      ::RegCloseKey(regCUSC);
      ::RegCloseKey(regCUS);

      // tell Windows that the association was changed
      ::SHChangeNotify(SHCNE_ASSOCCHANGED, SHCNF_IDLIST, NULL, NULL);
   } // UpdateRegistry()

   /////////////////////////////////////////////////////////////////////////////
   /// return kFALSE if option "-l" was specified as main programm command arg

   bool NeedSplash()
   {
      static bool once = true;
      TString arg;

      if (!once || gROOT->IsBatch()) return false;
      TString cmdline(::GetCommandLine());
      Int_t i = 0, from = 0;
      while (cmdline.Tokenize(arg, from, " ")) {
         arg.Strip(TString::kBoth);
         if (i == 0 && ((arg != "root") && (arg != "rootn") &&
             (arg != "root.exe") && (arg != "rootn.exe"))) return false;
         else if ((arg == "-l") || (arg == "-b")) return false;
         ++i;
      }
      if (once) {
         once = false;
         return true;
      }
      return false;
   }

   /////////////////////////////////////////////////////////////////////////////

   static void SetConsoleWindowName()
   {
      char pszNewWindowTitle[1024]; // contains fabricated WindowTitle
      char pszOldWindowTitle[1024]; // contains original WindowTitle
      HANDLE hStdout;
      CONSOLE_SCREEN_BUFFER_INFO csbiInfo;

      if (!::GetConsoleTitle(pszOldWindowTitle, 1024))
         return;
      // format a "unique" NewWindowTitle
      wsprintf(pszNewWindowTitle,"%d/%d", ::GetTickCount(), ::GetCurrentProcessId());
      // change current window title
      if (!::SetConsoleTitle(pszNewWindowTitle))
         return;
      // ensure window title has been updated
      ::Sleep(40);
      // look for NewWindowTitle
      gConsoleWindow = (ULong_t)::FindWindow(0, pszNewWindowTitle);
      if (gConsoleWindow) {
         // restore original window title
         ::ShowWindow((HWND)gConsoleWindow, SW_RESTORE);
         //::SetForegroundWindow((HWND)gConsoleWindow);
         ::SetConsoleTitle("ROOT session");
      }
      hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
      // adding the ENABLE_VIRTUAL_TERMINAL_PROCESSING flag would enable the ANSI control
      // character sequences (e.g. `\033[39m`), but then it breaks the WRAP_AT_EOL_OUTPUT
      ::SetConsoleMode(hStdout, ENABLE_PROCESSED_OUTPUT |
                       ENABLE_WRAP_AT_EOL_OUTPUT);
      if (!::GetConsoleScreenBufferInfo(hStdout, &csbiInfo))
         return;
      Gl_setwidth(csbiInfo.dwMaximumWindowSize.X);
   }

} // end unnamed namespace


///////////////////////////////////////////////////////////////////////////////
ClassImp(TWinNTSystem);

ULong_t gConsoleWindow = 0;

////////////////////////////////////////////////////////////////////////////////
///

Bool_t TWinNTSystem::HandleConsoleEvent()
{
   TSignalHandler *sh;
   TIter next(fSignalHandler);
   ESignals s;

   while (sh = (TSignalHandler*)next()) {
      s = sh->GetSignal();
      if (s == kSigInterrupt) {
         sh->Notify();
         Throw(SIGINT);
         return kTRUE;
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// ctor

TWinNTSystem::TWinNTSystem() : TSystem("WinNT", "WinNT System")
{
   fhProcess = ::GetCurrentProcess();

   WSADATA WSAData;
   int initwinsock = 0;

   if (initwinsock = ::WSAStartup(MAKEWORD(2, 0), &WSAData)) {
      Error("TWinNTSystem()","Starting sockets failed");
   }

   // use ::MessageBeep by default for TWinNTSystem
   fBeepDuration = 1;
   fBeepFreq     = 0;
   if (gEnv) {
      fBeepDuration = gEnv->GetValue("Root.System.BeepDuration", 1);
      fBeepFreq     = gEnv->GetValue("Root.System.BeepFreq", 0);
   }

   char *buf = new char[MAX_MODULE_NAME32 + 1];

#ifdef ROOTPREFIX
   if (gSystem->Getenv("ROOTIGNOREPREFIX")) {
#endif
   // set ROOTSYS
   HMODULE hModCore = ::GetModuleHandle("libCore.dll");
   if (hModCore) {
      ::GetModuleFileName(hModCore, buf, MAX_MODULE_NAME32 + 1);
      char *pLibName = strstr(buf, "libCore.dll");
      if (pLibName) {
         --pLibName; // skip trailing \\ or /
         while (--pLibName >= buf && *pLibName != '\\' && *pLibName != '/');
         *pLibName = 0; // replace trailing \\ or / with 0
         TString check_path = buf;
         check_path += "\\etc";
         // look for $ROOTSYS (it should contain the "etc" subdirectory)
         while (buf[0] && GetFileAttributes(check_path.Data()) == INVALID_FILE_ATTRIBUTES) {
            while (--pLibName >= buf && *pLibName != '\\' && *pLibName != '/');
            *pLibName = 0;
            check_path = buf;
            check_path += "\\etc";
         }
         if (buf[0]) {
            Setenv("ROOTSYS", buf);
            TString path = buf;
            path += "\\bin;";
            path += Getenv("PATH");
            Setenv("PATH", path.Data());
         }
      }
   }
#ifdef ROOTPREFIX
   }
#endif

   UpdateRegistry(this, buf);

   delete [] buf;
}

////////////////////////////////////////////////////////////////////////////////
/// dtor

TWinNTSystem::~TWinNTSystem()
{
   // Revert back the accuracy of Sleep() without needing to link to winmm.lib
   typedef UINT (WINAPI* LPTIMEENDPERIOD)( UINT uPeriod );
   HINSTANCE hInstWinMM = LoadLibrary( "winmm.dll" );
   if( hInstWinMM ) {
      LPTIMEENDPERIOD pTimeEndPeriod = (LPTIMEENDPERIOD)GetProcAddress( hInstWinMM, "timeEndPeriod" );
      if( NULL != pTimeEndPeriod )
         pTimeEndPeriod(1);
      FreeLibrary(hInstWinMM);
   }
   // Clean up the WinSocket connectios
   ::WSACleanup();

   if (gGlobalEvent) {
      ::ResetEvent(gGlobalEvent);
      ::CloseHandle(gGlobalEvent);
      gGlobalEvent = 0;
   }
   if (gTimerThreadHandle) {
      ::TerminateThread(gTimerThreadHandle, 0);
      ::CloseHandle(gTimerThreadHandle);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize WinNT system interface.

Bool_t TWinNTSystem::Init()
{
   if (TSystem::Init())
      return kTRUE;

   fReadmask = new TFdSet;
   fWritemask = new TFdSet;
   fReadready = new TFdSet;
   fWriteready = new TFdSet;
   fSignals = new TFdSet;
   fNfd    = 0;

   //--- install default handlers
   // Actually: don't. If we want a stack trace we need a context for the
   // signal. Signals don't have one. If we don't handle them, Windows will
   // raise an exception, which has a context, and which is handled by
   // ExceptionFilter.
   //WinNTSignal(kSigChild,                 SigHandler);
   //WinNTSignal(kSigBus,                   SigHandler);
   WinNTSignal(kSigSegmentationViolation, SigHandler);
   WinNTSignal(kSigIllegalInstruction,    SigHandler);
   WinNTSignal(kSigAbort,                 SigHandler);
   //WinNTSignal(kSigSystem,                SigHandler);
   //WinNTSignal(kSigPipe,                  SigHandler);
   //WinNTSignal(kSigAlarm,                 SigHandler);
   WinNTSignal(kSigFloatingException,     SigHandler);
   ::SetUnhandledExceptionFilter(ExceptionFilter);

   fSigcnt = 0;

   // This is a fallback in case TROOT::GetRootSys() can't determine ROOTSYS
   gRootDir = ROOT::FoundationUtils::GetFallbackRootSys().c_str();

   // Increase the accuracy of Sleep() without needing to link to winmm.lib
   typedef UINT (WINAPI* LPTIMEBEGINPERIOD)( UINT uPeriod );
   HINSTANCE hInstWinMM = LoadLibrary( "winmm.dll" );
   if( hInstWinMM ) {
      LPTIMEBEGINPERIOD pTimeBeginPeriod = (LPTIMEBEGINPERIOD)GetProcAddress( hInstWinMM, "timeBeginPeriod" );
      if( NULL != pTimeBeginPeriod )
         pTimeBeginPeriod(1);
      FreeLibrary(hInstWinMM);
   }
   gTimerThreadHandle = ::CreateThread(NULL, NULL, (LPTHREAD_START_ROUTINE)ThreadStub,
                        this, NULL, NULL);

   gGlobalEvent = ::CreateEvent(NULL, TRUE, FALSE, NULL);
   fGUIThreadHandle = ::CreateThread( NULL, 0, &GUIThreadMessageProcessingLoop, 0, 0, &fGUIThreadId );

   char *buf = new char[MAX_MODULE_NAME32 + 1];
   HMODULE hModCore = ::GetModuleHandle("libCore.dll");
   if (hModCore) {
      ::GetModuleFileName(hModCore, buf, MAX_MODULE_NAME32 + 1);
      char *pLibName = strstr(buf, "libCore.dll");
      --pLibName; // remove trailing \\ or /
      *pLibName = 0;
      // add the directory containing libCore.dll in the dynamic search path
      if (buf[0]) AddDynamicPath(buf);
   }
   delete [] buf;
   SetConsoleWindowName();
   fGroupsInitDone = kFALSE;
   fFirstFile = kTRUE;

   return kFALSE;
}

//---- Misc --------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Base name of a file name. Base name of /user/root is root.
/// But the base name of '/' is '/'
///                      'c:\' is 'c:\'

const char *TWinNTSystem::BaseName(const char *name)
{
   // BB 28/10/05 : Removed (commented out) StrDup() :
   // - To get same behaviour on Windows and on Linux
   // - To avoid the need to use #ifdefs
   // - Solve memory leaks (mainly in TTF::SetTextFont())
   // No need for the calling routine to use free() anymore.

   if (name) {
      int idx = 0;
      const char *symbol=name;

      // Skip leading blanks
      while ( (*symbol == ' ' || *symbol == '\t') && *symbol) symbol++;

      if (*symbol) {
         if (isalpha(symbol[idx]) && symbol[idx+1] == ':') idx = 2;
         if ( (symbol[idx] == '/'  ||  symbol[idx] == '\\')  &&  symbol[idx+1] == '\0') {
            //return StrDup(symbol);
            return symbol;
         }
      } else {
         Error("BaseName", "name = 0");
         return nullptr;
      }
      char *cp;
      char *bslash = (char *)strrchr(&symbol[idx],'\\');
      char *rslash = (char *)strrchr(&symbol[idx],'/');
      if (cp = (std::max)(rslash, bslash)) {
         //return StrDup(++cp);
         return ++cp;
      }
      //return StrDup(&symbol[idx]);
      return &symbol[idx];
   }
   Error("BaseName", "name = 0");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the application name (from command line, argv[0]) and copy it in
/// gProgName. Copy the application pathname in gProgPath.

void TWinNTSystem::SetProgname(const char *name)
{
   ULong_t  idot = 0;
   char *dot = nullptr;
   char *progname;
   char *fullname = nullptr; // the program name with extension

  // On command prompt the progname can be supplied with no extension (under Windows)
   ULong_t namelen=name ? strlen(name) : 0;
   if (name && namelen > 0) {
      // Check whether the name contains "extention"
      fullname = new char[namelen+5];
      strlcpy(fullname, name,namelen+5);
      if ( !strrchr(fullname, '.') )
         strlcat(fullname, ".exe",namelen+5);

      progname = StrDup(BaseName(fullname));
      dot = strrchr(progname, '.');
      idot = dot ? (ULong_t)(dot - progname) : strlen(progname);

      char *which = nullptr;

      if (IsAbsoluteFileName(fullname) && !AccessPathName(fullname)) {
         which = StrDup(fullname);
      } else {
         which = Which(Form("%s;%s", WorkingDirectory(), Getenv("PATH")), progname);
      }

      if (which) {
         TString dirname;
         char driveletter = DriveName(which);
         TString d = GetDirName(which);

         if (driveletter) {
            dirname.Form("%c:%s", driveletter, d.Data());
         } else {
            dirname = d;
         }

         gProgPath = StrDup(dirname);
      } else {
         // Do not issue a warning - ROOT is not using gProgPath anyway.
         // Warning("SetProgname",
         //   "Cannot find this program named \"%s\" (Did you create a TApplication? Is this program in your %%PATH%%?)",
         //   fullname);
         gProgPath = StrDup(WorkingDirectory());
      }

      // Cut the extension for progname off
      progname[idot] = '\0';
      gProgName = StrDup(progname);
      if (which) delete [] which;
      delete[] fullname;
      delete[] progname;
   }
   if (::NeedSplash()) {
      new TWin32SplashThread(FALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return system error string.

const char *TWinNTSystem::GetError()
{
   Int_t err = GetErrno();
   if (err == 0 && GetLastErrorString() != "")
      return GetLastErrorString();
   if (err < 0 || err >= sys_nerr) {
      static TString error_msg;
      error_msg.Form("errno out of range %d", err);
      return error_msg;
   }
   return sys_errlist[err];
}

////////////////////////////////////////////////////////////////////////////////
/// Return the system's host name.

const char *TWinNTSystem::HostName()
{
   if (fHostname == "")
      fHostname = ::getenv("COMPUTERNAME");
   if (fHostname == "") {
      // This requires a DNS query - but we need it for fallback
      char hn[64];
      DWORD il = sizeof(hn);
      ::GetComputerName(hn, &il);
      fHostname = hn;
   }
   return fHostname;
}

////////////////////////////////////////////////////////////////////////////////
/// Beep. If freq==0 (the default for TWinNTSystem), use ::MessageBeep.
/// Otherwise ::Beep with freq and duration.

void TWinNTSystem::DoBeep(Int_t freq /*=-1*/, Int_t duration /*=-1*/) const
{
   if (freq == 0) {
      ::MessageBeep(-1);
      return;
   }
   if (freq < 37) freq = 440;
   if (duration < 0) duration = 100;
   ::Beep(freq, duration);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the (static part of) the event handler func for GUI messages.

void TWinNTSystem::SetGUIThreadMsgHandler(ThreadMsgFunc_t func)
{
   gGUIThreadMsgFunc = func;
}

////////////////////////////////////////////////////////////////////////////////
/// Hook to tell TSystem that the TApplication object has been created.

void TWinNTSystem::NotifyApplicationCreated()
{
   // send a dummy message to the GUI thread to kick it into life
   ::PostThreadMessage(fGUIThreadId, 0, NULL, 0L);
}


//---- EventLoop ---------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Add a file handler to the list of system file handlers. Only adds
/// the handler if it is not already in the list of file handlers.

void TWinNTSystem::AddFileHandler(TFileHandler *h)
{
   TSystem::AddFileHandler(h);
   if (h) {
      int fd = h->GetFd();
      if (!fd) return;

      if (h->HasReadInterest()) {
         fReadmask->Set(fd);
      }
      if (h->HasWriteInterest()) {
         fWritemask->Set(fd);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a file handler from the list of file handlers. Returns
/// the handler or 0 if the handler was not in the list of file handlers.

TFileHandler *TWinNTSystem::RemoveFileHandler(TFileHandler *h)
{
   if (!h) return nullptr;

   TFileHandler *oh = TSystem::RemoveFileHandler(h);
   if (oh) {       // found
      fReadmask->Clr(h->GetFd());
      fWritemask->Clr(h->GetFd());
   }
   return oh;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a signal handler to list of system signal handlers. Only adds
/// the handler if it is not already in the list of signal handlers.

void TWinNTSystem::AddSignalHandler(TSignalHandler *h)
{
   Bool_t set_console = kFALSE;
   ESignals  sig = h->GetSignal();

   if (sig == kSigInterrupt) {
      set_console = kTRUE;
      TSignalHandler *hs;
      TIter next(fSignalHandler);

      while ((hs = (TSignalHandler*) next())) {
         if (hs->GetSignal() == kSigInterrupt)
            set_console = kFALSE;
      }
   }
   TSystem::AddSignalHandler(h);

   // Add our handler to the list of the console handlers
   if (set_console)
      ::SetConsoleCtrlHandler((PHANDLER_ROUTINE)ConsoleSigHandler, TRUE);
   else
      WinNTSignal(h->GetSignal(), SigHandler);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a signal handler from list of signal handlers. Returns
/// the handler or 0 if the handler was not in the list of signal handlers.

TSignalHandler *TWinNTSystem::RemoveSignalHandler(TSignalHandler *h)
{
   if (!h) return nullptr;

   int sig = h->GetSignal();

   if (sig = kSigInterrupt) {
      Bool_t last = kTRUE;
      TSignalHandler *hs;
      TIter next(fSignalHandler);

      while ((hs = (TSignalHandler*) next())) {
         if (hs->GetSignal() == kSigInterrupt)
            last = kFALSE;
      }
      // Remove our handler from the list of the console handlers
      if (last)
         ::SetConsoleCtrlHandler((PHANDLER_ROUTINE)ConsoleSigHandler, FALSE);
   }
   return TSystem::RemoveSignalHandler(h);
}

////////////////////////////////////////////////////////////////////////////////
/// If reset is true reset the signal handler for the specified signal
/// to the default handler, else restore previous behaviour.

void TWinNTSystem::ResetSignal(ESignals sig, Bool_t reset)
{
   //FIXME!
}

////////////////////////////////////////////////////////////////////////////////
/// Reset signals handlers to previous behaviour.

void TWinNTSystem::ResetSignals()
{
   //FIXME!
}

////////////////////////////////////////////////////////////////////////////////
/// If ignore is true ignore the specified signal, else restore previous
/// behaviour.

void TWinNTSystem::IgnoreSignal(ESignals sig, Bool_t ignore)
{
   // FIXME!
}

////////////////////////////////////////////////////////////////////////////////
/// Print a stack trace, if gEnv entry "Root.Stacktrace" is unset or 1,
/// and if the image helper functions can be found (see InitImagehlpFunctions()).
/// The stack trace is printed for each thread; if fgXcptContext is set (e.g.
/// because there was an exception) use it to define the current thread's context.
/// For each frame in the stack, the frame's module name, the frame's function
/// name, and the frame's line number are printed.

void TWinNTSystem::StackTrace()
{
   if (!gEnv->GetValue("Root.Stacktrace", 1))
      return;

   HANDLE snapshot = ::CreateToolhelp32Snapshot(TH32CS_SNAPTHREAD,::GetCurrentProcessId());

   std::cerr.flush();
   fflush (stderr);

   if (!InitImagehlpFunctions()) {
      std::cerr << "No stack trace: cannot find (functions in) dbghelp.dll!" << std::endl;
      return;
   }

   // what system are we on?
   SYSTEM_INFO sysInfo;
   ::GetSystemInfo(&sysInfo);
   DWORD machineType = IMAGE_FILE_MACHINE_I386;
   switch (sysInfo.wProcessorArchitecture) {
      case PROCESSOR_ARCHITECTURE_AMD64:
         machineType = IMAGE_FILE_MACHINE_AMD64;
         break;
      case PROCESSOR_ARCHITECTURE_IA64:
         machineType = IMAGE_FILE_MACHINE_IA64;
         break;
   }

   DWORD currentThreadID = ::GetCurrentThreadId();
   DWORD currentProcessID = ::GetCurrentProcessId();

   if (snapshot == INVALID_HANDLE_VALUE) return;

   THREADENTRY32 threadentry;
   threadentry.dwSize = sizeof(THREADENTRY32);
   if (!::Thread32First(snapshot, &threadentry)) return;

   std::cerr << std::endl << "==========================================" << std::endl;
   std::cerr << "=============== STACKTRACE ===============" << std::endl;
   std::cerr << "==========================================" << std::endl << std::endl;
   UInt_t iThread = 0;
   do {
      if (threadentry.th32OwnerProcessID != currentProcessID)
         continue;
      HANDLE thread = ::OpenThread(THREAD_GET_CONTEXT|THREAD_SUSPEND_RESUME|THREAD_QUERY_INFORMATION,
         FALSE, threadentry.th32ThreadID);
      CONTEXT context;
      memset(&context, 0, sizeof(CONTEXT));

      if (threadentry.th32ThreadID != currentThreadID) {
         ::SuspendThread(thread);
         context.ContextFlags = CONTEXT_ALL;
         ::GetThreadContext(thread, &context);
         ::ResumeThread(thread);
      } else {
         if (fgXcptContext) {
            context = *fgXcptContext;
         } else {
            typedef void (WINAPI *RTLCCTXT)(PCONTEXT);
            RTLCCTXT p2RtlCCtxt = (RTLCCTXT) ::GetProcAddress(
               GetModuleHandle("kernel32.dll"), "RtlCaptureContext");
            if (p2RtlCCtxt) {
               context.ContextFlags = CONTEXT_ALL;
               p2RtlCCtxt(&context);
            }
         }
      }

      STACKFRAME64 frame;
      ::ZeroMemory(&frame, sizeof(frame));

      frame.AddrPC.Mode      = AddrModeFlat;
      frame.AddrFrame.Mode   = AddrModeFlat;
      frame.AddrStack.Mode   = AddrModeFlat;
#if defined(_M_IX86)
      frame.AddrPC.Offset    = context.Eip;
      frame.AddrFrame.Offset = context.Ebp;
      frame.AddrStack.Offset = context.Esp;
#elif defined(_M_X64)
      frame.AddrPC.Offset    = context.Rip;
      frame.AddrFrame.Offset = context.Rsp;
      frame.AddrStack.Offset = context.Rsp;
#elif defined(_M_IA64)
      frame.AddrPC.Offset    = context.StIIP;
      frame.AddrFrame.Offset = context.IntSp;
      frame.AddrStack.Offset = context.IntSp;
      frame.AddrBStore.Offset= context.RsBSP;
#else
      std::cerr << "Stack traces not supported on your architecture yet." << std::endl;
      return;
#endif

      Bool_t bFirst = kTRUE;
      while (_StackWalk64(machineType, (HANDLE)::GetCurrentProcess(), thread, (LPSTACKFRAME64)&frame,
         (LPVOID)&context, (PREAD_PROCESS_MEMORY_ROUTINE)NULL, (PFUNCTION_TABLE_ACCESS_ROUTINE)_SymFunctionTableAccess64,
         (PGET_MODULE_BASE_ROUTINE)_SymGetModuleBase64, NULL)) {
         if (bFirst)
            std::cerr << std::endl << "================ Thread " << iThread++ << " ================" << std::endl;
         if (!bFirst || threadentry.th32ThreadID != currentThreadID) {
            const std::string moduleName   = GetModuleName(frame.AddrPC.Offset);
            const std::string functionName = GetFunctionName(frame.AddrPC.Offset);
            std::cerr << "  " << moduleName << functionName << std::endl;
         }
         bFirst = kFALSE;
      }
      ::CloseHandle(thread);
   } while (::Thread32Next(snapshot, &threadentry));

   std::cerr << std::endl << "==========================================" << std::endl;
   std::cerr << "============= END STACKTRACE =============" << std::endl;
   std::cerr << "==========================================" << std::endl << std::endl;
   ::CloseHandle(snapshot);
   _SymCleanup(GetCurrentProcess());
}

////////////////////////////////////////////////////////////////////////////////
/// Return the bitmap of conditions that trigger a floating point exception.

Int_t TWinNTSystem::GetFPEMask()
{
   Int_t mask = 0;
   UInt_t oldmask = _statusfp( );

   if (oldmask & _EM_INVALID  )   mask |= kInvalid;
   if (oldmask & _EM_ZERODIVIDE)  mask |= kDivByZero;
   if (oldmask & _EM_OVERFLOW )   mask |= kOverflow;
   if (oldmask & _EM_UNDERFLOW)   mask |= kUnderflow;
   if (oldmask & _EM_INEXACT  )   mask |= kInexact;

   return mask;
}

////////////////////////////////////////////////////////////////////////////////
/// Set which conditions trigger a floating point exception.
/// Return the previous set of conditions.

Int_t TWinNTSystem::SetFPEMask(Int_t mask)
{
   Int_t old = GetFPEMask();

   UInt_t newm = 0;
   if (mask & kInvalid  )   newm |= _EM_INVALID;
   if (mask & kDivByZero)   newm |= _EM_ZERODIVIDE;
   if (mask & kOverflow )   newm |= _EM_OVERFLOW;
   if (mask & kUnderflow)   newm |= _EM_UNDERFLOW;
   if (mask & kInexact  )   newm |= _EM_INEXACT;

   UInt_t cm = ::_statusfp();
   cm &= ~newm;
   ::_controlfp(cm , _MCW_EM);

   return old;
}

////////////////////////////////////////////////////////////////////////////////
/// process pending events, i.e. DispatchOneEvent(kTRUE)

Bool_t TWinNTSystem::ProcessEvents()
{
   return TSystem::ProcessEvents();
}

////////////////////////////////////////////////////////////////////////////////
/// Dispatch a single event in TApplication::Run() loop

void TWinNTSystem::DispatchOneEvent(Bool_t pendingOnly)
{
   // check for keyboard events
   if (pendingOnly && gGlobalEvent) ::SetEvent(gGlobalEvent);

   Bool_t pollOnce = pendingOnly;

   while (1) {
      if (_kbhit()) {
         if (gROOT->GetApplication()) {
            gApplication->HandleTermInput();
            if (gSplash) {    // terminate splash window after first key press
               delete gSplash;
               gSplash = 0;
            }
            if (!pendingOnly) {
               return;
            }
         }
      }
      if (gROOT->IsLineProcessing() && (!gVirtualX || !gVirtualX->IsCmdThread())) {
         if (!pendingOnly) {
            // yield execution to another thread that is ready to run
            // if no other thread is ready, sleep 1 ms before to return
            if (gGlobalEvent) {
               ::WaitForSingleObject(gGlobalEvent, 1);
               ::ResetEvent(gGlobalEvent);
            }
            return;
         }
      }
      // first handle any GUI events
      if (gXDisplay && !gROOT->IsBatch()) {
         if (gXDisplay->Notify()) {
            if (!pendingOnly) {
               return;
            }
         }
      }

      // check for file descriptors ready for reading/writing
      if ((fNfd > 0) && fFileHandler && (fFileHandler->GetSize() > 0)) {
         if (CheckDescriptors()) {
            if (!pendingOnly) {
               return;
            }
         }
      }
      fNfd = 0;
      fReadready->Zero();
      fWriteready->Zero();

      if (pendingOnly && !pollOnce)
         return;

      // check synchronous signals
      if (fSigcnt > 0 && fSignalHandler->GetSize() > 0) {
         if (CheckSignals(kTRUE)) {
            if (!pendingOnly) {
               return;
            }
         }
      }
      fSigcnt = 0;
      fSignals->Zero();

      // handle past due timers
      Long_t nextto;
      if (fTimers && fTimers->GetSize() > 0) {
         if (DispatchTimers(kTRUE)) {
            // prevent timers from blocking the rest types of events
            nextto = NextTimeOut(kTRUE);
            if (nextto > (kItimerResolution>>1) || nextto == -1) {
               return;
            }
         }
      }

      // if in pendingOnly mode poll once file descriptor activity
      nextto = NextTimeOut(kTRUE);
      if (pendingOnly) {
         if (fFileHandler && fFileHandler->GetSize() == 0)
            return;
         nextto = 0;
         pollOnce = kFALSE;
      }

      if (fReadmask && !fReadmask->GetBits() &&
          fWritemask && !fWritemask->GetBits()) {
         // yield execution to another thread that is ready to run
         // if no other thread is ready, sleep 1 ms before to return
         if (!pendingOnly && gGlobalEvent) {
            ::WaitForSingleObject(gGlobalEvent, 1);
            ::ResetEvent(gGlobalEvent);
         }
         return;
      }

      *fReadready  = *fReadmask;
      *fWriteready = *fWritemask;

      fNfd = WinNTSelect(fReadready, fWriteready, nextto);

      // serious error has happened -> reset all file descrptors
      if ((fNfd < 0) && (fNfd != -2)) {
         int rc, i;

         for (i = 0; i < fReadmask->GetCount(); i++) {
            TFdSet t;
            Int_t fd = fReadmask->GetFd(i);
            t.Set(fd);
            if (fReadmask->IsSet(fd)) {
               rc = WinNTSelect(&t, 0, 0);
               if (rc < 0 && rc != -2) {
                  ::SysError("DispatchOneEvent", "select: read error on %d\n", fd);
                  fReadmask->Clr(fd);
               }
            }
         }

         for (i = 0; i < fWritemask->GetCount(); i++) {
            TFdSet t;
            Int_t fd = fWritemask->GetFd(i);
            t.Set(fd);

            if (fWritemask->IsSet(fd)) {
               rc = WinNTSelect(0, &t, 0);
               if (rc < 0 && rc != -2) {
                  ::SysError("DispatchOneEvent", "select: write error on %d\n", fd);
                  fWritemask->Clr(fd);
               }
            }
            t.Clr(fd);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Exit from event loop.

void TWinNTSystem::ExitLoop()
{
   TSystem::ExitLoop();
}

//---- handling of system events -----------------------------------------------
////////////////////////////////////////////////////////////////////////////////
/// Handle and dispatch signals.

void TWinNTSystem::DispatchSignals(ESignals sig)
{
   if (sig == kSigInterrupt) {
      fSignals->Set(sig);
      fSigcnt++;
   }
   else {
      if (gExceptionHandler) {
         //sig is ESignal, should it be mapped to the correct signal number?
         if (sig == kSigFloatingException) _fpreset();
         gExceptionHandler->HandleException(sig);
      } else {
         if (sig == kSigAbort)
            return;
         //map to the real signal code + set the
         //high order bit to indicate a signal (?)
         StackTrace();
         if (TROOT::Initialized()) {
             ::Throw(sig);
         }
      }
      Abort(-1);
   }

   // check a-synchronous signals
   if (fSigcnt > 0 && fSignalHandler->GetSize() > 0)
      CheckSignals(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if some signals were raised and call their Notify() member.

Bool_t TWinNTSystem::CheckSignals(Bool_t sync)
{
   TSignalHandler *sh;
   Int_t sigdone = -1;
   {
      TIter next(fSignalHandler);

      while (sh = (TSignalHandler*)next()) {
         if (sync == sh->IsSync()) {
            ESignals sig = sh->GetSignal();
            if ((fSignals->IsSet(sig) && sigdone == -1) || sigdone == sig) {
               if (sigdone == -1) {
                  fSignals->Clr(sig);
                  sigdone = sig;
                  fSigcnt--;
               }
               sh->Notify();
            }
         }
      }
   }
   if (sigdone != -1) return kTRUE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if there is activity on some file descriptors and call their
/// Notify() member.

Bool_t TWinNTSystem::CheckDescriptors()
{
   TFileHandler *fh;
   Int_t  fddone = -1;
   Bool_t read   = kFALSE;

   TOrdCollectionIter it((TOrdCollection*)fFileHandler);

   while ((fh = (TFileHandler*) it.Next())) {
      Int_t fd = fh->GetFd();
      if (!fd) continue; // ignore TTermInputHandler

      if ((fReadready->IsSet(fd) && fddone == -1) ||
          (fddone == fd && read)) {
         if (fddone == -1) {
            fReadready->Clr(fd);
            fddone = fd;
            read = kTRUE;
            fNfd--;
         }
         fh->ReadNotify();
      }
      if ((fWriteready->IsSet(fd) && fddone == -1) ||
          (fddone == fd && !read)) {
         if (fddone == -1) {
            fWriteready->Clr(fd);
            fddone = fd;
            read = kFALSE;
            fNfd--;
         }
         fh->WriteNotify();
      }
   }
   if (fddone != -1) return kTRUE;

   return kFALSE;
}

//---- Directories -------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Make a file system directory. Returns 0 in case of success and
/// -1 if the directory could not be created (either already exists or
/// illegal path name).
/// If 'recursive' is true, makes parent directories as needed.

int TWinNTSystem::mkdir(const char *name, Bool_t recursive)
{
   if (recursive) {
      TString dirname = GetDirName(name);
      if (dirname.Length() == 0) {
         // well we should not have to make the root of the file system!
         // (and this avoid infinite recursions!)
         return 0;
      }
      if (IsAbsoluteFileName(name)) {
         // For some good reason DirName strips off the drive letter
         // (if present), we need it to make the directory on the
         // right disk, so let's put it back!
         const char driveletter = DriveName(name);
         if (driveletter) {
            dirname.Prepend(":");
            dirname.Prepend(driveletter);
         }
      }
      if (AccessPathName(dirname, kFileExists)) {
         int res = this->mkdir(dirname, kTRUE);
         if (res) return res;
      }
      if (!AccessPathName(name, kFileExists)) {
         return -1;
      }
   }
   return MakeDirectory(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Make a WinNT file system directory. Returns 0 in case of success and
/// -1 if the directory could not be created (either already exists or
/// illegal path name).

int  TWinNTSystem::MakeDirectory(const char *name)
{
   TSystem *helper = FindHelper(name);
   if (helper) {
      return helper->MakeDirectory(name);
   }
   const char *proto = (strstr(name, "file:///")) ? "file://" : "file:";
#ifdef WATCOM
   // It must be as follows
   if (!name) return 0;
   return ::mkdir(StripOffProto(name, proto));
#else
   // but to be in line with TUnixSystem I did like this
   if (!name) return 0;
   return ::_mkdir(StripOffProto(name, proto));
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Close a WinNT file system directory.

void TWinNTSystem::FreeDirectory(void *dirp)
{
   TSystem *helper = FindHelper(0, dirp);
   if (helper) {
      helper->FreeDirectory(dirp);
      return;
   }

   if (dirp) {
      ::FindClose(dirp);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the next directory entry.

const char *TWinNTSystem::GetDirEntry(void *dirp)
{
   TSystem *helper = FindHelper(0, dirp);
   if (helper) {
      return helper->GetDirEntry(dirp);
   }

   if (dirp) {
      HANDLE searchFile = (HANDLE)dirp;
      if (fFirstFile) {
         // when calling TWinNTSystem::OpenDirectory(), the fFindFileData
         // structure is filled by a call to FindFirstFile().
         // So first returns this one, before calling FindNextFile()
         fFirstFile = kFALSE;
         return (const char *)fFindFileData.cFileName;
      }
      if (::FindNextFile(searchFile, &fFindFileData)) {
         return (const char *)fFindFileData.cFileName;
      }
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Change directory.

Bool_t TWinNTSystem::ChangeDirectory(const char *path)
{
   Bool_t ret = (Bool_t) (::chdir(path) == 0);
   if (fWdpath != "")
      fWdpath = "";   // invalidate path cache
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Inline function to check for a double-backslash at the
/// beginning of a string
///

__inline BOOL DBL_BSLASH(LPCTSTR psz)
{
   return (psz[0] == TEXT('\\') && psz[1] == TEXT('\\'));
}

////////////////////////////////////////////////////////////////////////////////
/// Returns TRUE if the given string is a UNC path.
///
/// TRUE
///      "\\foo\bar"
///      "\\foo"         <- careful
///      "\\"
/// FALSE
///      "\foo"
///      "foo"
///      "c:\foo"

BOOL PathIsUNC(LPCTSTR pszPath)
{
   return DBL_BSLASH(pszPath);
}

#pragma data_seg(".text", "CODE")
const TCHAR c_szColonSlash[] = TEXT(":\\");
#pragma data_seg()

////////////////////////////////////////////////////////////////////////////////
///
/// check if a path is a root
///
/// returns:
///  TRUE for "\" "X:\" "\\foo\asdf" "\\foo\"
///  FALSE for others
///

BOOL PathIsRoot(LPCTSTR pPath)
{
   if (!IsDBCSLeadByte(*pPath)) {
      if (!lstrcmpi(pPath + 1, c_szColonSlash))
         // "X:\" case
         return TRUE;
   }
   if ((*pPath == TEXT('\\')) && (*(pPath + 1) == 0))
      // "\" case
      return TRUE;
   if (DBL_BSLASH(pPath)) {
      // smells like UNC name
      LPCTSTR p;
      int cBackslashes = 0;
      for (p = pPath + 2; *p; p = CharNext(p)) {
         if (*p == TEXT('\\') && (++cBackslashes > 1))
            return FALSE;   // not a bare UNC name, therefore not a root dir
      }
      // end of string with only 1 more backslash
      // must be a bare UNC, which looks like a root dir
      return TRUE;
   }
   return FALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a directory. Returns 0 if directory does not exist.

void *TWinNTSystem::OpenDirectory(const char *fdir)
{
   TSystem *helper = FindHelper(fdir);
   if (helper) {
      return helper->OpenDirectory(fdir);
   }

   const char *proto = (strstr(fdir, "file:///")) ? "file://" : "file:";
   const char *sdir = StripOffProto(fdir, proto);

   char *dir = new char[MAX_PATH];
   if (IsShortcut(sdir)) {
      if (!ResolveShortCut(sdir, dir, MAX_PATH))
         strlcpy(dir, sdir,MAX_PATH);
   }
   else
      strlcpy(dir, sdir,MAX_PATH);

   int nche = strlen(dir)+3;
   char *entry = new char[nche];
   struct _stati64 finfo;

   if(PathIsUNC(dir)) {
      strlcpy(entry, dir,nche);
      if ((entry[strlen(dir)-1] == '/') || (entry[strlen(dir)-1] == '\\' )) {
         entry[strlen(dir)-1] = '\0';
      }
      if(PathIsRoot(entry)) {
         strlcat(entry,"\\",nche);
      }
      if (_stati64(entry, &finfo) < 0) {
         delete [] entry;
         delete [] dir;
         return nullptr;
      }
   } else {
      strlcpy(entry, dir,nche);
      if ((entry[strlen(dir)-1] == '/') || (entry[strlen(dir)-1] == '\\' )) {
         if(!PathIsRoot(entry))
            entry[strlen(dir)-1] = '\0';
      }
      if (_stati64(entry, &finfo) < 0) {
         delete [] entry;
         delete [] dir;
         return nullptr;
      }
   }

   if (finfo.st_mode & S_IFDIR) {
      strlcpy(entry, dir,nche);
      if (!(entry[strlen(dir)-1] == '/' || entry[strlen(dir)-1] == '\\' )) {
         strlcat(entry,"\\",nche);
      }
      if (entry[strlen(dir)-1] == ' ')
         entry[strlen(dir)-1] = '\0';
      strlcat(entry,"*",nche);

      HANDLE searchFile;
      searchFile = ::FindFirstFile(entry, &fFindFileData);
      if (searchFile == INVALID_HANDLE_VALUE) {
         ((TWinNTSystem *)gSystem)->Error( "Unable to find' for reading:", entry);
         delete [] entry;
         delete [] dir;
         return nullptr;
      }
      delete [] entry;
      delete [] dir;
      fFirstFile = kTRUE;
      return searchFile;
   }

   delete [] entry;
   delete [] dir;
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the working directory for the default drive

const char *TWinNTSystem::WorkingDirectory()
{
   return WorkingDirectory('\0');
}

//////////////////////////////////////////////////////////////////////////////
/// Return the working directory for the default drive

std::string TWinNTSystem::GetWorkingDirectory() const
{
   char *wdpath = GetWorkingDirectory('\0');
   std::string cwd;
   if (wdpath) {
      cwd = wdpath;
      free(wdpath);
   }
   return cwd;
}

////////////////////////////////////////////////////////////////////////////////
///  Return working directory for the selected drive
///  driveletter == 0 means return the working durectory for the default drive

const char *TWinNTSystem::WorkingDirectory(char driveletter)
{
   char *wdpath = GetWorkingDirectory(driveletter);
   if (wdpath) {
      fWdpath = wdpath;

      // Make sure the drive letter is upper case
      if (fWdpath[1] == ':')
         fWdpath[0] = toupper(fWdpath[0]);

      free(wdpath);
   }
   return fWdpath;
}

//////////////////////////////////////////////////////////////////////////////
///  Return working directory for the selected drive (helper function).
///  The caller must free the return value.

char *TWinNTSystem::GetWorkingDirectory(char driveletter) const
{
   char *wdpath = nullptr;
   char drive = driveletter ? toupper( driveletter ) - 'A' + 1 : 0;

   // don't use cache as user can call chdir() directly somewhere else
   //if (fWdpath != "" )
   //   return fWdpath;

   if (!(wdpath = ::_getdcwd( (int)drive, wdpath, kMAXPATHLEN))) {
      free(wdpath);
      Warning("WorkingDirectory", "getcwd() failed");
      return nullptr;
   }

   return wdpath;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the user's home directory.

const char *TWinNTSystem::HomeDirectory(const char *userName)
{
   static char mydir[kMAXPATHLEN] = "./";
   FillWithHomeDirectory(userName, mydir);
   return mydir;
}

//////////////////////////////////////////////////////////////////////////////
/// Return the user's home directory.

std::string TWinNTSystem::GetHomeDirectory(const char *userName) const
{
   char mydir[kMAXPATHLEN] = "./";
   FillWithHomeDirectory(userName, mydir);
   return std::string(mydir);
}

//////////////////////////////////////////////////////////////////////////////
/// Fill buffer with user's home directory.

void TWinNTSystem::FillWithHomeDirectory(const char *userName, char *mydir) const
{
   const char *h = nullptr;
   if (!(h = ::getenv("home"))) h = ::getenv("HOME");

   if (h) {
      strlcpy(mydir, h,kMAXPATHLEN);
   } else {
      // for Windows NT HOME might be defined as either $(HOMESHARE)/$(HOMEPATH)
      //                                         or     $(HOMEDRIVE)/$(HOMEPATH)
      h = ::getenv("HOMESHARE");
      if (!h)  h = ::getenv("HOMEDRIVE");
      if (h) {
         strlcpy(mydir, h,kMAXPATHLEN);
         h = ::getenv("HOMEPATH");
         if(h) strlcat(mydir, h,kMAXPATHLEN);
      }
      // on Windows Vista HOME is usually defined as $(USERPROFILE)
      if (!h) {
         h = ::getenv("USERPROFILE");
         if (h) strlcpy(mydir, h,kMAXPATHLEN);
      }
   }
   // Make sure the drive letter is upper case
   if (mydir[1] == ':')
      mydir[0] = toupper(mydir[0]);
}


////////////////////////////////////////////////////////////////////////////////
/// Return a user configured or systemwide directory to create
/// temporary files in.

const char *TWinNTSystem::TempDirectory() const
{
   const char *dir =  gSystem->Getenv("TEMP");
   if (!dir)   dir =  gSystem->Getenv("TEMPDIR");
   if (!dir)   dir =  gSystem->Getenv("TEMP_DIR");
   if (!dir)   dir =  gSystem->Getenv("TMP");
   if (!dir)   dir =  gSystem->Getenv("TMPDIR");
   if (!dir)   dir =  gSystem->Getenv("TMP_DIR");
   if (!dir) dir = "c:\\";

   return dir;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a secure temporary file by appending a unique
/// 6 letter string to base. The file will be created in
/// a standard (system) directory or in the directory
/// provided in dir. The full filename is returned in base
/// and a filepointer is returned for safely writing to the file
/// (this avoids certain security problems). Returns 0 in case
/// of error.

FILE *TWinNTSystem::TempFileName(TString &base, const char *dir)
{
   char tmpName[MAX_PATH];

   ::GetTempFileName(dir ? dir : TempDirectory(), base.Data(), 0, tmpName);
   base = tmpName;
   FILE *fp = fopen(tmpName, "w+");

   if (!fp) ::SysError("TempFileName", "error opening %s", tmpName);

   return fp;
}

//---- Paths & Files -----------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Get list of volumes (drives) mounted on the system.
/// The returned TList must be deleted by the user using "delete".

TList *TWinNTSystem::GetVolumes(Option_t *opt) const
{
   Int_t   curdrive;
   UInt_t  type;
   TString sDrive, sType;
   char    szFs[32];

   if (!opt || !opt[0]) {
      return 0;
   }

   // prevent the system dialog box to pop-up if a drive is empty
   UINT nOldErrorMode = ::SetErrorMode(SEM_FAILCRITICALERRORS);
   TList *drives = new TList();
   drives->SetOwner();
   // Save current drive
   curdrive = _getdrive();
   if (strstr(opt, "cur")) {
      *szFs='\0';
      sDrive.Form("%c:", (curdrive + 'A' - 1));
      sType.Form("Unknown Drive (%s)", sDrive.Data());
      ::GetVolumeInformation(Form("%s\\", sDrive.Data()), NULL, 0, NULL, NULL,
                             NULL, (LPSTR)szFs, 32);
      type = ::GetDriveType(sDrive.Data());
      switch (type) {
         case DRIVE_UNKNOWN:
         case DRIVE_NO_ROOT_DIR:
            break;
         case DRIVE_REMOVABLE:
            sType.Form("Removable Disk (%s)", sDrive.Data());
            break;
         case DRIVE_FIXED:
            sType.Form("Local Disk (%s)", sDrive.Data());
            break;
         case DRIVE_REMOTE:
            sType.Form("Network Drive (%s) (%s)", szFs, sDrive.Data());
            break;
         case DRIVE_CDROM:
            sType.Form("CD/DVD Drive (%s)", sDrive.Data());
            break;
         case DRIVE_RAMDISK:
            sType.Form("RAM Disk (%s)", sDrive.Data());
            break;
      }
      drives->Add(new TNamed(sDrive.Data(), sType.Data()));
   }
   else if (strstr(opt, "all")) {
      TCHAR szTemp[512];
      szTemp[0] = '\0';
      if (::GetLogicalDriveStrings(511, szTemp)) {
         TCHAR szDrive[3] = TEXT(" :");
         TCHAR* p = szTemp;
         do {
            // Copy the drive letter to the template string
            *szDrive = *p;
            *szFs='\0';
            sDrive.Form("%s", szDrive);
            // skip floppy drives, to avoid accessing them each time...
            if ((sDrive == "A:") || (sDrive == "B:")) {
               while (*p++);
               continue;
            }
            sType.Form("Unknown Drive (%s)", sDrive.Data());
            ::GetVolumeInformation(Form("%s\\", sDrive.Data()), NULL, 0, NULL,
                                   NULL, NULL, (LPSTR)szFs, 32);
            type = ::GetDriveType(sDrive.Data());
            switch (type) {
               case DRIVE_UNKNOWN:
               case DRIVE_NO_ROOT_DIR:
                  break;
               case DRIVE_REMOVABLE:
                  sType.Form("Removable Disk (%s)", sDrive.Data());
                  break;
               case DRIVE_FIXED:
                  sType.Form("Local Disk (%s)", sDrive.Data());
                  break;
               case DRIVE_REMOTE:
                  sType.Form("Network Drive (%s) (%s)", szFs, sDrive.Data());
                  break;
               case DRIVE_CDROM:
                  sType.Form("CD/DVD Drive (%s)", sDrive.Data());
                  break;
               case DRIVE_RAMDISK:
                  sType.Form("RAM Disk (%s)", sDrive.Data());
                  break;
            }
            drives->Add(new TNamed(sDrive.Data(), sType.Data()));
            // Go to the next NULL character.
            while (*p++);
         } while (*p); // end of string
      }
   }
   // restore previous error mode
   ::SetErrorMode(nOldErrorMode);
   return drives;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the directory name in pathname. DirName of c:/user/root is /user.
/// It creates output with 'new char []' operator. Returned string has to
/// be deleted.

const char *TWinNTSystem::DirName(const char *pathname)
{
   fDirNameBuffer = GetDirName(pathname);
   return fDirNameBuffer.c_str();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the directory name in pathname. DirName of c:/user/root is /user.
/// DirName of c:/user/root/ is /user/root.

TString TWinNTSystem::GetDirName(const char *pathname)
{
   // Create a buffer to keep the path name
   if (pathname) {
      if (strchr(pathname, '/') || strchr(pathname, '\\')) {
         const char *rslash = strrchr(pathname, '/');
         const char *bslash = strrchr(pathname, '\\');
         const char *r = std::max(rslash, bslash);
         const char *ptr = pathname;
         while (ptr <= r) {
            if (*ptr == ':') {
               // Windows path may contain a drive letter
               // For NTFS ":" may be a "stream" delimiter as well
               pathname =  ptr + 1;
               break;
            }
            ptr++;
         }
         int len =  r - pathname;
         if (len > 0)
            return TString(pathname, len);
      }
   }
   return "";
}

////////////////////////////////////////////////////////////////////////////////
/// Return the drive letter in pathname. DriveName of 'c:/user/root' is 'c'
///
///   Input:
///     - pathname - the string containing file name
///
///   Return:
///     - Letter representing the drive letter in the file name
///     - The current drive if the pathname has no drive assigment
///     - 0 if pathname is an empty string  or uses UNC syntax
///
///   Note:
///      It doesn't check whether pathname represents a 'real' filename.
///      This subroutine looks for 'single letter' followed by a ':'.

const char TWinNTSystem::DriveName(const char *pathname)
{
   if (!pathname)    return 0;
   if (!pathname[0]) return 0;

   const char *lpchar;
   lpchar = pathname;

   // Skip blanks
   while(*lpchar == ' ') lpchar++;

   if (isalpha((int)*lpchar) && *(lpchar+1) == ':') {
      return *lpchar;
   }
   // Test UNC syntax
   if ( (*lpchar == '\\' || *lpchar == '/' ) &&
        (*(lpchar+1) == '\\' || *(lpchar+1) == '/') ) return 0;

   // return the current drive
   return DriveName(WorkingDirectory());
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if dir is an absolute pathname.

Bool_t TWinNTSystem::IsAbsoluteFileName(const char *dir)
{
   if (dir) {
      int idx = 0;
      if (strchr(dir,':')) idx = 2;
      return  (dir[idx] == '/' || dir[idx] == '\\');
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert a pathname to a unix pathname. E.g. form \user\root to /user/root.
/// General rules for applications creating names for directories and files or
/// processing names supplied by the user include the following:
///
///  *  Use any character in the current code page for a name, but do not use
///     a path separator, a character in the range 0 through 31, or any character
///     explicitly disallowed by the file system. A name can contain characters
///     in the extended character set (128-255).
///  *  Use the backslash (\), the forward slash (/), or both to separate
///     components in a path. No other character is acceptable as a path separator.
///  *  Use a period (.) as a directory component in a path to represent the
///     current directory.
///  *  Use two consecutive periods (..) as a directory component in a path to
///     represent the parent of the current directory.
///  *  Use a period (.) to separate components in a directory name or filename.
///  *  Do not use the following characters in directory names or filenames, because
///     they are reserved for Windows:
///                      < > : " / \ |
///  *  Do not use reserved words, such as aux, con, and prn, as filenames or
///     directory names.
///  *  Process a path as a null-terminated string. The maximum length for a path
///     is given by MAX_PATH.
///  *  Do not assume case sensitivity. Consider names such as OSCAR, Oscar, and
///     oscar to be the same.

const char *TWinNTSystem::UnixPathName(const char *name)
{
   const int kBufSize = 1024;
   TTHREAD_TLS_ARRAY(char, kBufSize, temp);

   strlcpy(temp, name, kBufSize);
   char *currentChar = temp;

   // This can not change the size of the string.
   while (*currentChar != '\0') {
      if (*currentChar == '\\') *currentChar = '/';
      currentChar++;
   }
   return temp;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns FALSE if one can access a file using the specified access mode.
/// Mode is the same as for the WinNT access(2) function.
/// Attention, bizarre convention of return value!!

Bool_t TWinNTSystem::AccessPathName(const char *path, EAccessMode mode)
{
   TSystem *helper = FindHelper(path);
   if (helper)
      return helper->AccessPathName(path, mode);

   // prevent the system dialog box to pop-up if a drive is empty
   UINT nOldErrorMode = ::SetErrorMode(SEM_FAILCRITICALERRORS);
   if (mode==kExecutePermission)
      // cannot test on exe - use read instead
      mode=kReadPermission;
   const char *proto = (strstr(path, "file:///")) ? "file://" : "file:";
   if (::_access(StripOffProto(path, proto), mode) == 0) {
      // restore previous error mode
      ::SetErrorMode(nOldErrorMode);
      return kFALSE;
   }
   GetLastErrorString() = GetError();
   // restore previous error mode
   ::SetErrorMode(nOldErrorMode);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns TRUE if the url in 'path' points to the local file system.
/// This is used to avoid going through the NIC card for local operations.

Bool_t TWinNTSystem::IsPathLocal(const char *path)
{
   TSystem *helper = FindHelper(path);
   if (helper)
      return helper->IsPathLocal(path);

   return TSystem::IsPathLocal(path);
}

////////////////////////////////////////////////////////////////////////////////
/// Concatenate a directory and a file name.

const char *TWinNTSystem::PrependPathName(const char *dir, TString& name)
{
   if (name == ".") name = "";
   if (dir && dir[0]) {
      // Test whether the last symbol of the directory is a separator
      char last = dir[strlen(dir) - 1];
      if (last != '/' && last != '\\') {
         name.Prepend('\\');
      }
      name.Prepend(dir);
      name.ReplaceAll("/", "\\");
   }
   return name.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a file. If overwrite is true and file already exists the
/// file will be overwritten. Returns 0 when successful, -1 in case
/// of failure, -2 in case the file already exists and overwrite was false.

int TWinNTSystem::CopyFile(const char *f, const char *t, Bool_t overwrite)
{
   if (AccessPathName(f, kReadPermission)) return -1;
   if (!AccessPathName(t) && !overwrite) return -2;

   Bool_t ret = ::CopyFileA(f, t, kFALSE);

   if (!ret) return -1;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Rename a file. Returns 0 when successful, -1 in case of failure.

int TWinNTSystem::Rename(const char *f, const char *t)
{
   int ret = ::rename(f, t);
   GetLastErrorString() = GetError();
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Get info about a file. Info is returned in the form of a FileStat_t
/// structure (see TSystem.h).
/// The function returns 0 in case of success and 1 if the file could
/// not be stat'ed.

int TWinNTSystem::GetPathInfo(const char *path, FileStat_t &buf)
{
   TSystem *helper = FindHelper(path);
   if (helper)
      return helper->GetPathInfo(path, buf);

   struct _stati64 sbuf;

   // Remove trailing backslashes
   const char *proto = (strstr(path, "file:///")) ? "file://" : "file:";
   char *newpath = StrDup(StripOffProto(path, proto));
   int l = strlen(newpath);
   while (l > 1) {
      if (newpath[--l] != '\\' || newpath[--l] != '/') {
         break;
      }
      newpath[l] = '\0';
   }

   if (newpath && ::_stati64(newpath, &sbuf) >= 0) {

      buf.fDev    = sbuf.st_dev;
      buf.fIno    = sbuf.st_ino;
      buf.fMode   = sbuf.st_mode;
      buf.fUid    = sbuf.st_uid;
      buf.fGid    = sbuf.st_gid;
      buf.fSize   = sbuf.st_size;
      buf.fMtime  = sbuf.st_mtime;
      buf.fIsLink = IsShortcut(newpath); // kFALSE;

      char *lpath = new char[MAX_PATH];
      if (IsShortcut(newpath)) {
         struct _stati64 sbuf2;
         if (ResolveShortCut(newpath, lpath, MAX_PATH)) {
            if (::_stati64(lpath, &sbuf2) >= 0) {
               buf.fMode   = sbuf2.st_mode;
            }
         }
      }
      delete [] lpath;

      delete [] newpath;
      return 0;
   }
   delete [] newpath;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get info about a file system: id, bsize, bfree, blocks.
/// Id      is file system type (machine dependend, see statfs())
/// Bsize   is block size of file system
/// Blocks  is total number of blocks in file system
/// Bfree   is number of free blocks in file system
/// The function returns 0 in case of success and 1 if the file system could
/// not be stat'ed.

int TWinNTSystem::GetFsInfo(const char *path, Long_t *id, Long_t *bsize,
                            Long_t *blocks, Long_t *bfree)
{
   // address of root directory of the file system
   LPCTSTR lpRootPathName = path;

   // address of name of the volume
   LPTSTR  lpVolumeNameBuffer = 0;
   DWORD   nVolumeNameSize = 0;

   DWORD   volumeSerialNumber;     // volume serial number
   DWORD   maximumComponentLength; // system's maximum filename length

   // file system flags
   DWORD fileSystemFlags;

   // address of name of file system
   char  fileSystemNameBuffer[512];
   DWORD nFileSystemNameSize = sizeof(fileSystemNameBuffer);

   // prevent the system dialog box to pop-up if the drive is empty
   UINT nOldErrorMode = ::SetErrorMode(SEM_FAILCRITICALERRORS);
   if (!::GetVolumeInformation(lpRootPathName,
                               lpVolumeNameBuffer, nVolumeNameSize,
                               &volumeSerialNumber,
                               &maximumComponentLength,
                               &fileSystemFlags,
                               fileSystemNameBuffer, nFileSystemNameSize)) {
      // restore previous error mode
      ::SetErrorMode(nOldErrorMode);
      return 1;
   }

   const char *fsNames[] = { "FAT", "NTFS" };
   int i;
   for (i = 0; i < 2; i++) {
      if (!strncmp(fileSystemNameBuffer, fsNames[i], nFileSystemNameSize))
         break;
   }
   *id = i;

   DWORD sectorsPerCluster;      // # sectors per cluster
   DWORD bytesPerSector;         // # bytes per sector
   DWORD numberOfFreeClusters;   // # free clusters
   DWORD totalNumberOfClusters;  // # total of clusters

   if (!::GetDiskFreeSpace(lpRootPathName,
                           &sectorsPerCluster,
                           &bytesPerSector,
                           &numberOfFreeClusters,
                           &totalNumberOfClusters)) {
      // restore previous error mode
      ::SetErrorMode(nOldErrorMode);
      return 1;
   }
   // restore previous error mode
   ::SetErrorMode(nOldErrorMode);

   *bsize  = sectorsPerCluster * bytesPerSector;
   *blocks = totalNumberOfClusters;
   *bfree  = numberOfFreeClusters;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a link from file1 to file2.

int TWinNTSystem::Link(const char *from, const char *to)
{
   struct   _stati64 finfo;
   char     winDrive[256];
   char     winDir[256];
   char     winName[256];
   char     winExt[256];
   char     linkname[1024];
   LPTSTR   lpszFilePart;
   TCHAR    szPath[MAX_PATH];
   DWORD    dwRet = 0;

   typedef BOOL (__stdcall *CREATEHARDLINKPROC)( LPCTSTR, LPCTSTR, LPSECURITY_ATTRIBUTES );
   static CREATEHARDLINKPROC _CreateHardLink = 0;

   HMODULE hModImagehlp = LoadLibrary( "Kernel32.dll" );
   if (!hModImagehlp)
      return -1;

#ifdef _UNICODE
   _CreateHardLink = (CREATEHARDLINKPROC) GetProcAddress( hModImagehlp, "CreateHardLinkW" );
#else
   _CreateHardLink = (CREATEHARDLINKPROC) GetProcAddress( hModImagehlp, "CreateHardLinkA" );
#endif
   if (!_CreateHardLink)
      return -1;

   dwRet = GetFullPathName(from, sizeof(szPath) / sizeof(TCHAR),
                           szPath, &lpszFilePart);

   if (_stati64(szPath, &finfo) < 0)
      return -1;

   if (finfo.st_mode & S_IFDIR)
      return -1;

   snprintf(linkname,1024,"%s",to);
   _splitpath(linkname,winDrive,winDir,winName,winExt);
   if ((!winDrive[0] ) &&
       (!winDir[0] ))  {
      _splitpath(szPath,winDrive,winDir,winName,winExt);
      snprintf(linkname,1024,"%s\\%s\\%s", winDrive, winDir, to);
   }
   else if (!winDrive[0])  {
      _splitpath(szPath,winDrive,winDir,winName,winExt);
      snprintf(linkname,1024,"%s\\%s", winDrive, to);
   }

   if (!_CreateHardLink(linkname, szPath, NULL))
      return -1;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a symlink from file1 to file2. Returns 0 when successful,
/// -1 in case of failure.

int TWinNTSystem::Symlink(const char *from, const char *to)
{
   HRESULT        hRes;                  /* Returned COM result code */
   IShellLink*    pShellLink;            /* IShellLink object pointer */
   IPersistFile*  pPersistFile;          /* IPersistFile object pointer */
   WCHAR          wszLinkfile[MAX_PATH]; /* pszLinkfile as Unicode string */
   int            iWideCharsWritten;     /* Number of wide characters written */
   DWORD          dwRet = 0;
   LPTSTR         lpszFilePart;
   TCHAR          szPath[MAX_PATH];

   hRes = E_INVALIDARG;
   if ((from == NULL) || (!from[0]) || (to == NULL) ||
       (!to[0]))
      return -1;

   // Make typedefs for some ole32.dll functions so that we can use them
   // with GetProcAddress
   typedef HRESULT (__stdcall *COINITIALIZEPROC)( LPVOID );
   static COINITIALIZEPROC _CoInitialize = 0;
   typedef void (__stdcall *COUNINITIALIZEPROC)( void );
   static COUNINITIALIZEPROC _CoUninitialize = 0;
   typedef HRESULT (__stdcall *COCREATEINSTANCEPROC)( REFCLSID, LPUNKNOWN, DWORD, REFIID, LPVOID );
   static COCREATEINSTANCEPROC _CoCreateInstance = 0;

   HMODULE hModImagehlp = LoadLibrary( "ole32.dll" );
   if (!hModImagehlp)
      return -1;

   _CoInitialize = (COINITIALIZEPROC) GetProcAddress( hModImagehlp, "CoInitialize" );
   if (!_CoInitialize)
      return -1;
   _CoUninitialize = (COUNINITIALIZEPROC) GetProcAddress( hModImagehlp, "CoUninitialize" );
   if (!_CoUninitialize)
      return -1;
   _CoCreateInstance = (COCREATEINSTANCEPROC) GetProcAddress( hModImagehlp, "CoCreateInstance" );
   if (!_CoCreateInstance)
      return -1;

   TString linkname(to);
   if (!linkname.EndsWith(".lnk"))
      linkname.Append(".lnk");

   _CoInitialize(NULL);

   // Retrieve the full path and file name of a specified file
   dwRet = GetFullPathName(from, sizeof(szPath) / sizeof(TCHAR),
                           szPath, &lpszFilePart);
   hRes = _CoCreateInstance(CLSID_ShellLink, NULL, CLSCTX_INPROC_SERVER,
                           IID_IShellLink, (LPVOID *)&pShellLink);
   if (SUCCEEDED(hRes)) {
      // Set the fields in the IShellLink object
      hRes = pShellLink->SetPath(szPath);
      // Use the IPersistFile object to save the shell link
      hRes = pShellLink->QueryInterface(IID_IPersistFile, (void **)&pPersistFile);
      if (SUCCEEDED(hRes)){
         iWideCharsWritten = MultiByteToWideChar(CP_ACP, 0, linkname.Data(), -1,
                                                 wszLinkfile, MAX_PATH);
         hRes = pPersistFile->Save(wszLinkfile, TRUE);
         pPersistFile->Release();
      }
      pShellLink->Release();
   }
   _CoUninitialize();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Unlink, i.e. remove, a file or directory.
///
/// If the file is currently open by the current or another process Windows does not allow the file to be deleted and
/// the operation is a no-op.

int TWinNTSystem::Unlink(const char *name)
{
   TSystem *helper = FindHelper(name);
   if (helper)
      return helper->Unlink(name);

   struct _stati64 finfo;

   if (_stati64(name, &finfo) < 0) {
      return -1;
   }

   if (finfo.st_mode & S_IFDIR) {
      return ::_rmdir(name);
   } else {
      return ::_unlink(name);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Make descriptor fd non-blocking.

int TWinNTSystem::SetNonBlock(int fd)
{
   if (::ioctlsocket(fd, FIONBIO, (u_long *)1) == SOCKET_ERROR) {
      ::SysError("SetNonBlock", "ioctlsocket");
      return -1;
   }
   return 0;
}

// expand the metacharacters as in the shell

static const char
   *shellMeta  = "~*[]{}?$%",
   *shellStuff = "(){}<>\"'",
   shellEscape = '\\';

////////////////////////////////////////////////////////////////////////////////
/// Expand a pathname getting rid of special shell characaters like ~.$, etc.

Bool_t TWinNTSystem::ExpandPathName(TString &patbuf0)
{
   const char *patbuf = (const char *)patbuf0;
   const char *p;
   char   *cmd = nullptr;
   char  *q;

   Int_t old_level = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kFatal; // Explicitly remove all messages
   if (patbuf0.BeginsWith("\\")) {
      const char driveletter = DriveName(patbuf);
      if (driveletter) {
         patbuf0.Prepend(":");
         patbuf0.Prepend(driveletter);
      }
   }
   TUrl urlpath(patbuf0, kTRUE);
   TString proto = urlpath.GetProtocol();
   gErrorIgnoreLevel = old_level;
   if (!proto.EqualTo("file")) // don't expand urls!!!
      return kFALSE;

   // skip the "file:" protocol, if any
   if (patbuf0.BeginsWith("file:"))
      patbuf += 5;

   // skip leading blanks
   while (*patbuf == ' ') {
      patbuf++;
   }

   // skip leading ':'
   while (*patbuf == ':') {
      patbuf++;
   }

   // skip leading ';'
   while (*patbuf == ';') {
      patbuf++;
   }

   // Transform a Unix list of directories into a Windows list
   // by changing the separator from ':' into ';'
   for (q = (char*)patbuf; *q; q++) {
      if ( *q == ':' ) {
         // We are avoiding substitution in the case of
         // ....;c:.... and of ...;root:/... where root can be any url protocol
         if ( (((q-2)>patbuf) && ( (*(q-2)!=';') || !isalpha(*(q-1)) )) &&
              *(q+1)!='/' ) {
            *q=';';
         }
      }
   }
   // any shell meta characters ?
   for (p = patbuf; *p; p++) {
      if (strchr(shellMeta, *p)) {
         goto needshell;
      }
   }
   return kFALSE;

needshell:

   // Because (problably) we built with cygwin, the path name like:
   //     LOCALS~1\\Temp
   // gets extended to
   //     LOCALSc:\\Devel
   // The most likely cause is that '~' is used with Unix semantic of the
   // home directory (and it also cuts the path short after ... who knows why!)
   // So we need to detect this case and prevents its expansion :(.

   char replacement[4];

   // intentionally a non visible, unlikely character
   for (int k = 0; k<3; k++) replacement[k] = 0x1;

   replacement[3] = 0x0;
   Ssiz_t pos = 0;
   TRegexp TildaNum = "~[0-9]";

   while ( (pos = patbuf0.Index(TildaNum,pos)) != kNPOS ) {
      patbuf0.Replace(pos, 1, replacement);
   }

   // escape shell quote characters
   // EscChar(patbuf, stuffedPat, sizeof(stuffedPat), shellStuff, shellEscape);
   ExpandFileName(patbuf0);
   Int_t lbuf = ::ExpandEnvironmentStrings(
                                 patbuf0.Data(), // pointer to string with environment variables
                                 cmd,            // pointer to string with expanded environment variables
                                 0               // maximum characters in expanded string
                              );
   if (lbuf > 0) {
      cmd = new char[lbuf+1];
      ::ExpandEnvironmentStrings(
                               patbuf0.Data(), // pointer to string with environment variables
                               cmd,            // pointer to string with expanded environment variables
                               lbuf            // maximum characters in expanded string
                               );
      patbuf0 = cmd;
      patbuf0.ReplaceAll(replacement, "~");
      delete [] cmd;
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Expand a pathname getting rid of special shell characaters like ~.$, etc.
/// User must delete returned string.

char *TWinNTSystem::ExpandPathName(const char *path)
{
   char newpath[MAX_PATH];
   if (IsShortcut(path)) {
      if (!ResolveShortCut(path, newpath, MAX_PATH))
         strlcpy(newpath, path, MAX_PATH);
   }
   else
      strlcpy(newpath, path, MAX_PATH);
   TString patbuf = newpath;
   if (ExpandPathName(patbuf))
      return nullptr;

   return StrDup(patbuf.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Set the file permission bits. Returns -1 in case or error, 0 otherwise.
/// On windows mode can only be a combination of "user read" (0400),
/// "user write" (0200) or "user read | user write" (0600). Any other value
/// for mode are ignored.

int TWinNTSystem::Chmod(const char *file, UInt_t mode)
{
   return ::_chmod(file, mode);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the process file creation mode mask.

int TWinNTSystem::Umask(Int_t mask)
{
   return ::umask(mask);
}

////////////////////////////////////////////////////////////////////////////////
/// Set a files modification and access times. If actime = 0 it will be
/// set to the modtime. Returns 0 on success and -1 in case of error.

int TWinNTSystem::Utime(const char *file, Long_t modtime, Long_t actime)
{
   if (AccessPathName(file, kWritePermission)) {
      Error("Utime", "need write permission for %s to change utime", file);
      return -1;
   }
   if (!actime) actime = modtime;

   struct utimbuf t;
   t.actime  = (time_t)actime;
   t.modtime = (time_t)modtime;
   return ::utime(file, &t);
}

////////////////////////////////////////////////////////////////////////////////
/// Find location of file in a search path.
/// User must delete returned string. Returns 0 in case file is not found.

const char *TWinNTSystem::FindFile(const char *search, TString& infile, EAccessMode mode)
{
   // Windows cannot check on execution mode - all we can do is kReadPermission
   if (mode==kExecutePermission)
      mode=kReadPermission;

   // Expand parameters

   gSystem->ExpandPathName(infile);
   // Check whether this infile has the absolute path first
   if (IsAbsoluteFileName(infile.Data()) ) {
      if (!AccessPathName(infile.Data(), mode))
         return infile.Data();
      infile = "";
      return nullptr;
   }
   TString exsearch(search);
   gSystem->ExpandPathName(exsearch);

   // Need to use Windows delimiters
   Int_t lastDelim = -1;
   for(int i=0; i < exsearch.Length(); ++i) {
      switch( exsearch[i] ) {
         case ':':
            // Replace the ':' unless there are after a disk suffix (aka ;c:\mydirec...)
            if (i-lastDelim!=2) exsearch[i] = ';';
            lastDelim = i;
            break;
         case ';': lastDelim = i; break;
      }
   }

   // Check access
   struct stat finfo;
   char name[kMAXPATHLEN];
   char *lpFilePart = nullptr;
   if (::SearchPath(exsearch.Data(), infile.Data(), NULL, kMAXPATHLEN, name, &lpFilePart) &&
       ::access(name, mode) == 0 && stat(name, &finfo) == 0 &&
       finfo.st_mode & S_IFREG) {
      if (gEnv->GetValue("Root.ShowPath", 0)) {
         Printf("Which: %s = %s", infile, name);
      }
      infile = name;
      return infile.Data();
   }
   infile = "";
   return nullptr;
}

//---- Users & Groups ----------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Collect local users and groups accounts information

Bool_t TWinNTSystem::InitUsersGroups()
{
   // Net* API functions allowed and OS is Windows NT/2000/XP
   if ((gEnv->GetValue("WinNT.UseNetAPI", 0)) && (::GetVersion() < 0x80000000)) {
      fActUser = -1;
      fNbGroups = fNbUsers = 0;
      HINSTANCE netapi = ::LoadLibrary("netapi32.DLL");
      if (!netapi) return kFALSE;

      p2NetApiBufferFree  = (pfn1)::GetProcAddress(netapi, "NetApiBufferFree");
      p2NetUserGetInfo  = (pfn2)::GetProcAddress(netapi, "NetUserGetInfo");
      p2NetLocalGroupGetMembers  = (pfn3)::GetProcAddress(netapi, "NetLocalGroupGetMembers");
      p2NetLocalGroupEnum = (pfn4)::GetProcAddress(netapi, "NetLocalGroupEnum");

      if (!p2NetApiBufferFree || !p2NetUserGetInfo ||
          !p2NetLocalGroupGetMembers || !p2NetLocalGroupEnum) return kFALSE;

      GetNbGroups();

      fGroups = (struct group *)calloc(fNbGroups, sizeof(struct group));
      for(int i=0;i<fNbGroups;i++) {
        fGroups[i].gr_mem = (char **)calloc(fNbUsers, sizeof (char*));
      }
      fPasswords = (struct passwd *)calloc(fNbUsers, sizeof(struct passwd));

      CollectGroups();
      ::FreeLibrary(netapi);
   }
   fGroupsInitDone = kTRUE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TWinNTSystem::CountMembers(const char *lpszGroupName)
{
   NET_API_STATUS NetStatus = NERR_Success;
   LPBYTE Data = NULL;
   DWORD Index = 0, ResumeHandle = 0, Total = 0;
   LOCALGROUP_MEMBERS_INFO_1 *MemberInfo;
   WCHAR wszGroupName[256];
   int iRetOp = 0;
   DWORD dwLastError = 0;

   iRetOp = MultiByteToWideChar (
            (UINT)CP_ACP,                // code page
            (DWORD)MB_PRECOMPOSED,       // character-type options
            (LPCSTR)lpszGroupName,       // address of string to map
            (int)-1,                     // number of bytes in string
            (LPWSTR)wszGroupName,        // address of wide-character buffer
            (int)sizeof(wszGroupName) ); // size of buffer

   if (iRetOp == 0) {
      dwLastError = GetLastError();
      if (Data)
         p2NetApiBufferFree(Data);
      return FALSE;
   }

   // The NetLocalGroupGetMembers() API retrieves a list of the members
   // of a particular local group.
   NetStatus = p2NetLocalGroupGetMembers (NULL, wszGroupName, 1,
                            &Data, 8192, &Index, &Total, &ResumeHandle );

   if (NetStatus != NERR_Success || Data == NULL) {
      dwLastError = GetLastError();

      if (dwLastError == ERROR_ENVVAR_NOT_FOUND) {
         // This usually means that the current Group has no members.
         // We call NetLocalGroupGetMembers() again.
         // This time, we set the level to 0.
         // We do this just to confirm that the number of members in
         // this group is zero.
         NetStatus = p2NetLocalGroupGetMembers ( NULL, wszGroupName, 0,
                                  &Data, 8192, &Index, &Total, &ResumeHandle );
      }

      if (Data)
         p2NetApiBufferFree(Data);
      return FALSE;
   }

   fNbUsers += Total;
   MemberInfo = (LOCALGROUP_MEMBERS_INFO_1 *)Data;

   if (Data)
      p2NetApiBufferFree(Data);

   return TRUE;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TWinNTSystem::GetNbGroups()
{
   NET_API_STATUS NetStatus = NERR_Success;
   LPBYTE Data = NULL;
   DWORD Index = 0, ResumeHandle = 0, Total = 0, i;
   LOCALGROUP_INFO_0 *GroupInfo;
   char szAnsiName[256];
   DWORD dwLastError = 0;
   int  iRetOp = 0;

   NetStatus = p2NetLocalGroupEnum(NULL, 0, &Data, 8192, &Index,
                                    &Total, &ResumeHandle );

   if (NetStatus != NERR_Success || Data == NULL) {
      dwLastError = GetLastError();
      if (Data)
         p2NetApiBufferFree(Data);
      return FALSE;
   }

   fNbGroups = Total;
   GroupInfo = (LOCALGROUP_INFO_0 *)Data;
   for (i=0; i < Total; i++) {
      // Convert group name from UNICODE to ansi.
      iRetOp = WideCharToMultiByte (
               (UINT)CP_ACP,                    // code page
               (DWORD)0,                        // performance and mapping flags
               (LPCWSTR)(GroupInfo->lgrpi0_name), // address of wide-char string
               (int)-1,                        // number of characters in string
               (LPSTR)szAnsiName,            // address of buffer for new string
               (int)(sizeof(szAnsiName)),    // size of buffer
               (LPCSTR)NULL,     // address of default for unmappable characters
               (LPBOOL)NULL );     // address of flag set when default char used.

      // Now lookup all members of this group and record down their names and
      // SIDs into the output file.
      CountMembers((LPCTSTR)szAnsiName);

      GroupInfo++;
   }

   if (Data)
      p2NetApiBufferFree(Data);

   return TRUE;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Take the name and look up a SID so that we can get full
/// domain/user information
///

Long_t TWinNTSystem::LookupSID (const char *lpszAccountName, int what,
                                int &groupIdx, int &memberIdx)
{
   BOOL bRetOp = FALSE;
   PSID pSid = NULL;
   DWORD dwSidSize, dwDomainNameSize;
   BYTE bySidBuffer[MAX_SID_SIZE];
   char szDomainName[MAX_NAME_STRING];
   SID_NAME_USE sidType;
   PUCHAR puchar_SubAuthCount = NULL;
   SID_IDENTIFIER_AUTHORITY sid_identifier_authority;
   PSID_IDENTIFIER_AUTHORITY psid_identifier_authority = NULL;
   unsigned char j = 0;
   DWORD dwLastError = 0;

   pSid = (PSID)bySidBuffer;
   dwSidSize = sizeof(bySidBuffer);
   dwDomainNameSize = sizeof(szDomainName);

   bRetOp = LookupAccountName (
            (LPCTSTR)NULL,             // address of string for system name
            (LPCTSTR)lpszAccountName,  // address of string for account name
            (PSID)pSid,                // address of security identifier
            (LPDWORD)&dwSidSize,       // address of size of security identifier
            (LPTSTR)szDomainName,      // address of string for referenced domain
            (LPDWORD)&dwDomainNameSize,// address of size of domain string
            (PSID_NAME_USE)&sidType ); // address of SID-type indicator

   if (bRetOp == FALSE) {
      dwLastError = GetLastError();
      return -1;  // Unable to obtain Account SID.
   }

   bRetOp = IsValidSid((PSID)pSid);

   if (bRetOp == FALSE) {
      dwLastError = GetLastError();
      return -2;  // SID returned is invalid.
   }

   // Obtain via APIs the identifier authority value.
   psid_identifier_authority = GetSidIdentifierAuthority ((PSID)pSid);

   // Make a copy of it.
   memcpy (&sid_identifier_authority, psid_identifier_authority,
       sizeof(SID_IDENTIFIER_AUTHORITY));

   // Determine how many sub-authority values there are in the current SID.
   puchar_SubAuthCount = (PUCHAR)GetSidSubAuthorityCount((PSID)pSid);
   // Assign it to a more convenient variable.
   j = (unsigned char)(*puchar_SubAuthCount);
   // Now obtain all the sub-authority values from the current SID.
   DWORD dwSubAuth = 0;
   PDWORD pdwSubAuth = NULL;
   // Obtain the current sub-authority DWORD (referenced by a pointer)
   pdwSubAuth = (PDWORD)GetSidSubAuthority (
                (PSID)pSid,  // address of security identifier to query
                (DWORD)j-1); // index of subauthority to retrieve
   dwSubAuth = *pdwSubAuth;
   if(what == SID_MEMBER) {
       fPasswords[memberIdx].pw_uid = dwSubAuth;
       fPasswords[memberIdx].pw_gid = fGroups[groupIdx].gr_gid;
       fPasswords[memberIdx].pw_group = strdup(fGroups[groupIdx].gr_name);
   }
   else if(what == SID_GROUP) {
       fGroups[groupIdx].gr_gid = dwSubAuth;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///

Bool_t TWinNTSystem::CollectMembers(const char *lpszGroupName, int &groupIdx,
                                    int &memberIdx)
{

   NET_API_STATUS NetStatus = NERR_Success;
   LPBYTE Data = NULL;
   DWORD Index = 0, ResumeHandle = 0, Total = 0, i;
   LOCALGROUP_MEMBERS_INFO_1 *MemberInfo;
   char szAnsiMemberName[256];
   char szFullMemberName[256];
   char szMemberHomeDir[256];
   WCHAR wszGroupName[256];
   int iRetOp = 0;
   char  act_name[256];
   DWORD length = sizeof (act_name);
   DWORD dwLastError = 0;
   LPUSER_INFO_11  pUI11Buf = NULL;
   NET_API_STATUS  nStatus;

   iRetOp = MultiByteToWideChar (
            (UINT)CP_ACP,                // code page
            (DWORD)MB_PRECOMPOSED,       // character-type options
            (LPCSTR)lpszGroupName,       // address of string to map
            (int)-1,                     // number of bytes in string
            (LPWSTR)wszGroupName,        // address of wide-character buffer
            (int)sizeof(wszGroupName) ); // size of buffer

   if (iRetOp == 0) {
      dwLastError = GetLastError();
      if (Data)
         p2NetApiBufferFree(Data);
      return FALSE;
   }

   GetUserName (act_name, &length);

   // The NetLocalGroupGetMembers() API retrieves a list of the members
   // of a particular local group.
   NetStatus = p2NetLocalGroupGetMembers (NULL, wszGroupName, 1,
                            &Data, 8192, &Index, &Total, &ResumeHandle );

   if (NetStatus != NERR_Success || Data == NULL) {
      dwLastError = GetLastError();

      if (dwLastError == ERROR_ENVVAR_NOT_FOUND) {
         // This usually means that the current Group has no members.
         // We call NetLocalGroupGetMembers() again.
         // This time, we set the level to 0.
         // We do this just to confirm that the number of members in
         // this group is zero.
         NetStatus = p2NetLocalGroupGetMembers ( NULL, wszGroupName, 0,
                                  &Data, 8192, &Index, &Total, &ResumeHandle );
      }

      if (Data)
         p2NetApiBufferFree(Data);
      return FALSE;
   }

   MemberInfo = (LOCALGROUP_MEMBERS_INFO_1 *)Data;
   for (i=0; i < Total; i++) {
      iRetOp = WideCharToMultiByte (
               (UINT)CP_ACP,                     // code page
               (DWORD)0,                         // performance and mapping flags
               (LPCWSTR)(MemberInfo->lgrmi1_name), // address of wide-char string
               (int)-1,                         // number of characters in string
               (LPSTR)szAnsiMemberName,       // address of buffer for new string
               (int)(sizeof(szAnsiMemberName)), // size of buffer
               (LPCSTR)NULL,      // address of default for unmappable characters
               (LPBOOL)NULL );      // address of flag set when default char used.

      if (iRetOp == 0) {
         dwLastError = GetLastError();
      }

      fPasswords[memberIdx].pw_name = strdup(szAnsiMemberName);
      fPasswords[memberIdx].pw_passwd = strdup("");
      fGroups[groupIdx].gr_mem[i] = strdup(szAnsiMemberName);

      if(fActUser == -1 && !stricmp(fPasswords[memberIdx].pw_name,act_name))
                      fActUser = memberIdx;


      TCHAR szUserName[255]=TEXT("");
      MultiByteToWideChar(CP_ACP, 0, szAnsiMemberName, -1, (LPWSTR)szUserName, 255);
      //
      // Call the NetUserGetInfo function; specify level 10.
      //
      nStatus = p2NetUserGetInfo(NULL, (LPCWSTR)szUserName, 11, (LPBYTE *)&pUI11Buf);
      //
      // If the call succeeds, print the user information.
      //
      if (nStatus == NERR_Success) {
         if (pUI11Buf != NULL) {
            wsprintf(szFullMemberName,"%S",pUI11Buf->usri11_full_name);
            fPasswords[memberIdx].pw_gecos = strdup(szFullMemberName);
            wsprintf(szMemberHomeDir,"%S",pUI11Buf->usri11_home_dir);
            fPasswords[memberIdx].pw_dir = strdup(szMemberHomeDir);
         }
      }
      if((fPasswords[memberIdx].pw_gecos == NULL) || (strlen(fPasswords[memberIdx].pw_gecos) == 0))
         fPasswords[memberIdx].pw_gecos = strdup(fPasswords[memberIdx].pw_name);
      if((fPasswords[memberIdx].pw_dir == NULL) || (strlen(fPasswords[memberIdx].pw_dir) == 0))
         fPasswords[memberIdx].pw_dir = strdup("c:\\");
      //
      // Free the allocated memory.
      //
      if (pUI11Buf != NULL) {
         p2NetApiBufferFree(pUI11Buf);
         pUI11Buf = NULL;
      }

      /* Ensure SHELL is defined. */
      if (getenv("SHELL") == NULL)
         putenv ((GetVersion () & 0x80000000) ? "SHELL=command" : "SHELL=cmd");

      /* Set dir and shell from environment variables. */
      fPasswords[memberIdx].pw_shell = getenv("SHELL");

      // Find out the SID of the Member.
      LookupSID ((LPCTSTR)szAnsiMemberName, SID_MEMBER, groupIdx, memberIdx);
      memberIdx++;
      MemberInfo++;
   }
   if(fActUser == -1)  fActUser = 0;

   if (Data)
      p2NetApiBufferFree(Data);

   return TRUE;
}

////////////////////////////////////////////////////////////////////////////////
///

Bool_t TWinNTSystem::CollectGroups()
{
   NET_API_STATUS NetStatus = NERR_Success;
   LPBYTE Data = NULL;
   DWORD Index = 0, ResumeHandle = 0, Total = 0, i;
   LOCALGROUP_INFO_0 *GroupInfo;
   char szAnsiName[256];
   DWORD dwLastError = 0;
   int  iRetOp = 0, iGroupIdx = 0, iMemberIdx = 0;

   NetStatus = p2NetLocalGroupEnum(NULL, 0, &Data, 8192, &Index,
                                    &Total, &ResumeHandle );

   if (NetStatus != NERR_Success || Data == NULL) {
      dwLastError = GetLastError();
      if (Data)
         p2NetApiBufferFree(Data);
      return FALSE;
   }

   GroupInfo = (LOCALGROUP_INFO_0 *)Data;
   for (i=0; i < Total; i++) {
      // Convert group name from UNICODE to ansi.
      iRetOp = WideCharToMultiByte (
               (UINT)CP_ACP,                    // code page
               (DWORD)0,                        // performance and mapping flags
               (LPCWSTR)(GroupInfo->lgrpi0_name), // address of wide-char string
               (int)-1,                        // number of characters in string
               (LPSTR)szAnsiName,            // address of buffer for new string
               (int)(sizeof(szAnsiName)),    // size of buffer
               (LPCSTR)NULL,     // address of default for unmappable characters
               (LPBOOL)NULL );     // address of flag set when default char used.

      fGroups[iGroupIdx].gr_name = strdup(szAnsiName);
      fGroups[iGroupIdx].gr_passwd = strdup("");

      // Find out the SID of the Group.
      LookupSID ((LPCTSTR)szAnsiName, SID_GROUP, iGroupIdx, iMemberIdx);
      // Now lookup all members of this group and record down their names and
      // SIDs into the output file.
      CollectMembers((LPCTSTR)szAnsiName, iGroupIdx, iMemberIdx);

      iGroupIdx++;
      GroupInfo++;
   }

   if (Data)
      p2NetApiBufferFree(Data);

   return TRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the user's id. If user = 0, returns current user's id.

Int_t TWinNTSystem::GetUid(const char *user)
{
   if(!fGroupsInitDone)
      InitUsersGroups();

   // Net* API functions not allowed or OS not Windows NT/2000/XP
   if ((!gEnv->GetValue("WinNT.UseNetAPI", 0)) || (::GetVersion() >= 0x80000000)) {
      int   uid;
      char  name[256];
      DWORD length = sizeof (name);
      if (::GetUserName (name, &length)) {
         if (stricmp ("administrator", name) == 0)
            uid = 0;
         else
            uid = 123;
      }
      else {
         uid = 123;
      }
      return uid;
   }
   if (!user || !user[0])
      return fPasswords[fActUser].pw_uid;
   else {
      struct passwd *pwd = nullptr;
      for(int i=0;i<fNbUsers;i++) {
         if (!stricmp (user, fPasswords[i].pw_name)) {
            pwd = &fPasswords[i];
            break;
         }
      }
      if (pwd)
         return pwd->pw_uid;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the effective user id. The effective id corresponds to the
/// set id bit on the file being executed.

Int_t TWinNTSystem::GetEffectiveUid()
{
   if(!fGroupsInitDone)
      InitUsersGroups();

   // Net* API functions not allowed or OS not Windows NT/2000/XP
   if ((!gEnv->GetValue("WinNT.UseNetAPI", 0)) || (::GetVersion() >= 0x80000000)) {
      int   uid;
      char  name[256];
      DWORD length = sizeof (name);
      if (::GetUserName (name, &length)) {
         if (stricmp ("administrator", name) == 0)
            uid = 0;
         else
            uid = 123;
      }
      else {
         uid = 123;
      }
      return uid;
   }
   return fPasswords[fActUser].pw_uid;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the group's id. If group = 0, returns current user's group.

Int_t TWinNTSystem::GetGid(const char *group)
{
   if(!fGroupsInitDone)
      InitUsersGroups();

   // Net* API functions not allowed or OS not Windows NT/2000/XP
   if ((!gEnv->GetValue("WinNT.UseNetAPI", 0)) || (::GetVersion() >= 0x80000000)) {
      int   gid;
      char  name[256];
      DWORD length = sizeof (name);
      if (::GetUserName (name, &length)) {
         if (stricmp ("administrator", name) == 0)
            gid = 0;
         else
            gid = 123;
      }
      else {
         gid = 123;
      }
      return gid;
   }
   if (!group || !group[0])
      return fPasswords[fActUser].pw_gid;
   else {
      struct group *grp = nullptr;
      for(int i=0;i<fNbGroups;i++) {
         if (!stricmp (group, fGroups[i].gr_name)) {
            grp = &fGroups[i];
            break;
         }
      }
      if (grp)
         return grp->gr_gid;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the effective group id. The effective group id corresponds
/// to the set id bit on the file being executed.

Int_t TWinNTSystem::GetEffectiveGid()
{
   if(!fGroupsInitDone)
      InitUsersGroups();

   // Net* API functions not allowed or OS not Windows NT/2000/XP
   if ((!gEnv->GetValue("WinNT.UseNetAPI", 0)) || (::GetVersion() >= 0x80000000)) {
      int   gid;
      char  name[256];
      DWORD length = sizeof (name);
      if (::GetUserName (name, &length)) {
         if (stricmp ("administrator", name) == 0)
            gid = 0;
         else
            gid = 123;
      }
      else {
         gid = 123;
      }
      return gid;
   }
   return fPasswords[fActUser].pw_gid;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns all user info in the UserGroup_t structure. The returned
/// structure must be deleted by the user. In case of error 0 is returned.

UserGroup_t *TWinNTSystem::GetUserInfo(Int_t uid)
{
   if(!fGroupsInitDone)
      InitUsersGroups();

   // Net* API functions not allowed or OS not Windows NT/2000/XP
   if ((!gEnv->GetValue("WinNT.UseNetAPI", 0)) || (::GetVersion() >= 0x80000000)) {
      char  name[256];
      DWORD length = sizeof (name);
      UserGroup_t *ug = new UserGroup_t;
      if (::GetUserName (name, &length)) {
         ug->fUser = name;
         if (stricmp ("administrator", name) == 0) {
            ug->fUid = 0;
            ug->fGroup = "administrators";
         }
         else {
            ug->fUid = 123;
            ug->fGroup = "users";
         }
         ug->fGid = ug->fUid;
      }
      else {
         ug->fUser = "unknown";
         ug->fGroup = "unknown";
         ug->fUid = ug->fGid = 123;
      }
      ug->fPasswd = "";
      ug->fRealName = ug->fUser;
      ug->fShell = "command";
      return ug;
   }
   struct passwd *pwd = 0;
   if (uid == 0)
      pwd = &fPasswords[fActUser];
   else {
      for (int i = 0; i < fNbUsers; i++) {
         if (uid == fPasswords[i].pw_uid) {
            pwd = &fPasswords[i];
            break;
         }
      }
   }
   if (pwd) {
      UserGroup_t *ug = new UserGroup_t;
      ug->fUid      = pwd->pw_uid;
      ug->fGid      = pwd->pw_gid;
      ug->fUser     = pwd->pw_name;
      ug->fPasswd   = pwd->pw_passwd;
      ug->fRealName = pwd->pw_gecos;
      ug->fShell    = pwd->pw_shell;
      ug->fGroup    = pwd->pw_group;
      return ug;
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns all user info in the UserGroup_t structure. If user = 0, returns
/// current user's id info. The returned structure must be deleted by the
/// user. In case of error 0 is returned.

UserGroup_t *TWinNTSystem::GetUserInfo(const char *user)
{
   return GetUserInfo(GetUid(user));
}

////////////////////////////////////////////////////////////////////////////////
/// Returns all group info in the UserGroup_t structure. The only active
/// fields in the UserGroup_t structure for this call are:
///    fGid and fGroup
/// The returned structure must be deleted by the user. In case of
/// error 0 is returned.

UserGroup_t *TWinNTSystem::GetGroupInfo(Int_t gid)
{
   if(!fGroupsInitDone)
      InitUsersGroups();

   // Net* API functions not allowed or OS not Windows NT/2000/XP
   if ((!gEnv->GetValue("WinNT.UseNetAPI", 0)) || (::GetVersion() >= 0x80000000)) {
      char  name[256];
      DWORD length = sizeof (name);
      UserGroup_t *gr = new UserGroup_t;
      if (::GetUserName (name, &length)) {
         if (stricmp ("administrator", name) == 0) {
            gr->fGroup = "administrators";
            gr->fGid = 0;
         }
         else {
            gr->fGroup = "users";
            gr->fGid = 123;
         }
      }
      else {
         gr->fGroup = "unknown";
         gr->fGid = 123;
      }
      gr->fUid = 0;
      return gr;
   }
   struct group *grp = nullptr;
   for(int i=0;i<fNbGroups;i++) {
      if (gid == fGroups[i].gr_gid) {
         grp = &fGroups[i];
         break;
      }
   }
   if (grp) {
      UserGroup_t *gr = new UserGroup_t;
      gr->fUid   = 0;
      gr->fGid   = grp->gr_gid;
      gr->fGroup = grp->gr_name;
      return gr;
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns all group info in the UserGroup_t structure. The only active
/// fields in the UserGroup_t structure for this call are:
///    fGid and fGroup
/// If group = 0, returns current user's group. The returned structure
/// must be deleted by the user. In case of error 0 is returned.

UserGroup_t *TWinNTSystem::GetGroupInfo(const char *group)
{
   return GetGroupInfo(GetGid(group));
}

//---- environment manipulation ------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Set environment variable.

void TWinNTSystem::Setenv(const char *name, const char *value)
{
   ::_putenv(TString::Format("%s=%s", name, value));
}

////////////////////////////////////////////////////////////////////////////////
/// Get environment variable.

const char *TWinNTSystem::Getenv(const char *name)
{
   const char *env = ::getenv(name);
   if (!env) {
      if (::_stricmp(name,"home") == 0 ) {
        env = HomeDirectory();
      } else if (::_stricmp(name, "rootsys") == 0 ) {
        env = gRootDir;
      }
   }
   return env;
}

//---- Processes ---------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Execute a command.

int TWinNTSystem::Exec(const char *shellcmd)
{
   return ::system(shellcmd);
}

////////////////////////////////////////////////////////////////////////////////
/// Open a pipe.

FILE *TWinNTSystem::OpenPipe(const char *command, const char *mode)
{
  return ::_popen(command, mode);
}

////////////////////////////////////////////////////////////////////////////////
/// Close the pipe.

int TWinNTSystem::ClosePipe(FILE *pipe)
{
  return ::_pclose(pipe);
}

////////////////////////////////////////////////////////////////////////////////
/// Get process id.

int TWinNTSystem::GetPid()
{
   return ::getpid();
}

////////////////////////////////////////////////////////////////////////////////
/// Get current process handle

HANDLE TWinNTSystem::GetProcess()
{
  return fhProcess;
}

////////////////////////////////////////////////////////////////////////////////
/// Exit the application.

void TWinNTSystem::Exit(int code, Bool_t mode)
{
   // Insures that the files and sockets are closed before any library is unloaded
   // and before emptying CINT.
   // FIXME: Unify with TROOT::ShutDown.
   if (gROOT) {
      gROOT->CloseFiles();
      if (gROOT->GetListOfBrowsers()) {
         // GetListOfBrowsers()->Delete() creates problems when a browser is
         // created on the stack, calling CloseWindow() solves the problem
         if (gROOT->IsBatch())
            gROOT->GetListOfBrowsers()->Delete();
         else {
            TBrowser *b;
            TIter next(gROOT->GetListOfBrowsers());
            while ((b = (TBrowser*) next()))
               gROOT->ProcessLine(TString::Format("\
                  if (((TBrowser*)0x%lx)->GetBrowserImp() &&\
                      ((TBrowser*)0x%lx)->GetBrowserImp()->GetMainFrame()) \
                     ((TBrowser*)0x%lx)->GetBrowserImp()->GetMainFrame()->CloseWindow();\
                  else delete (TBrowser*)0x%lx", (ULong_t)b, (ULong_t)b, (ULong_t)b, (ULong_t)b));
         }
      }
   }
   TROOT::ShutDown();
   gVirtualX->CloseDisplay();

   if (mode) {
      ::exit(code);
   } else {
      ::_exit(code);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Abort the application.

void TWinNTSystem::Abort(int)
{
   IgnoreSignal(kSigAbort);
   ::abort();
}

//---- Standard output redirection ---------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Redirect standard output (stdout, stderr) to the specified file.
/// If the file argument is 0 the output is set again to stderr, stdout.
/// The second argument specifies whether the output should be added to the
/// file ("a", default) or the file be truncated before ("w").
/// This function saves internally the current state into a static structure.
/// The call can be made reentrant by specifying the opaque structure pointed
/// by 'h', which is filled with the relevant information. The handle 'h'
/// obtained on the first call must then be used in any subsequent call,
/// included ShowOutput, to display the redirected output.
/// Returns 0 on success, -1 in case of error.

Int_t TWinNTSystem::RedirectOutput(const char *file, const char *mode,
                                   RedirectHandle_t *h)
{
   FILE *fout, *ferr;
   static int fd1=0, fd2=0;
   static fpos_t pos1=0, pos2=0;
   // Instance to be used if the caller does not passes 'h'
   static RedirectHandle_t loch;
   Int_t rc = 0;

   // Which handle to use ?
   RedirectHandle_t *xh = (h) ? h : &loch;

   if (file) {
      // Make sure mode makes sense; default "a"
      const char *m = (mode[0] == 'a' || mode[0] == 'w') ? mode : "a";

      // Current file size
      xh->fReadOffSet = 0;
      if (m[0] == 'a') {
         // If the file exists, save the current size
         FileStat_t st;
         if (!gSystem->GetPathInfo(file, st))
            xh->fReadOffSet = (st.fSize > 0) ? st.fSize : xh->fReadOffSet;
      }
      xh->fFile = file;

      fflush(stdout);
      fgetpos(stdout, &pos1);
      fd1 = _dup(fileno(stdout));
      // redirect stdout & stderr
      if ((fout = freopen(file, m, stdout)) == 0) {
         SysError("RedirectOutput", "could not freopen stdout");
         if (fd1 > 0) {
            _dup2(fd1, fileno(stdout));
            close(fd1);
         }
         clearerr(stdout);
         fsetpos(stdout, &pos1);
         fd1 = fd2 = 0;
         return -1;
      }
      fflush(stderr);
      fgetpos(stderr, &pos2);
      fd2 = _dup(fileno(stderr));
      if ((ferr = freopen(file, m, stderr)) == 0) {
         SysError("RedirectOutput", "could not freopen stderr");
         if (fd1 > 0) {
            _dup2(fd1, fileno(stdout));
            close(fd1);
         }
         clearerr(stdout);
         fsetpos(stdout, &pos1);
         if (fd2 > 0) {
            _dup2(fd2, fileno(stderr));
            close(fd2);
         }
         clearerr(stderr);
         fsetpos(stderr, &pos2);
         fd1 = fd2 = 0;
         return -1;
      }
      if (m[0] == 'a') {
         fseek(fout, 0, SEEK_END);
         fseek(ferr, 0, SEEK_END);
      }
   } else {
      // Restore stdout & stderr
      fflush(stdout);
      if (fd1) {
         if (fd1 > 0) {
            if (_dup2(fd1, fileno(stdout))) {
               SysError("RedirectOutput", "could not restore stdout");
               rc = -1;
            }
            close(fd1);
         }
         clearerr(stdout);
         fsetpos(stdout, &pos1);
         fd1 = 0;
      }

      fflush(stderr);
      if (fd2) {
         if (fd2 > 0) {
            if (_dup2(fd2, fileno(stderr))) {
               SysError("RedirectOutput", "could not restore stderr");
               rc = -1;
            }
            close(fd2);
         }
         clearerr(stderr);
         fsetpos(stderr, &pos2);
         fd2 = 0;
      }

      // Reset the static instance, if using that
      if (xh == &loch)
         xh->Reset();
   }
   return rc;
}

//---- dynamic loading and linking ---------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Add a new directory to the dynamic path.

void TWinNTSystem::AddDynamicPath(const char *dir)
{
   if (dir) {
      TString oldpath = DynamicPath(0, kFALSE);
      oldpath.Append(";");
      oldpath.Append(dir);
      DynamicPath(oldpath);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the dynamic path (used to find shared libraries).

const char* TWinNTSystem::GetDynamicPath()
{
   return DynamicPath(0, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the dynamic path to a new value.
/// If the value of 'path' is zero, the dynamic path is reset to its
/// default value.

void TWinNTSystem::SetDynamicPath(const char *path)
{
   if (!path)
      DynamicPath(0, kTRUE);
   else
      DynamicPath(path);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns and updates sLib to the path of a dynamic library
///  (searches for library in the dynamic library search path).
/// If no file name extension is provided it tries .DLL.

const char *TWinNTSystem::FindDynamicLibrary(TString &sLib, Bool_t quiet)
{
   int len = sLib.Length();
   if (len > 4 && (!stricmp(sLib.Data()+len-4, ".dll"))) {
      if (gSystem->FindFile(GetDynamicPath(), sLib, kReadPermission))
         return sLib.Data();
   } else {
      TString sLibDll(sLib);
      sLibDll += ".dll";
      if (gSystem->FindFile(GetDynamicPath(), sLibDll, kReadPermission)) {
         sLibDll.Swap(sLib);
         return sLib.Data();
      }
   }

   if (!quiet) {
      Error("DynamicPathName",
            "%s does not exist in %s,\nor has wrong file extension (.dll)",
             sLib.Data(), GetDynamicPath());
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Load a shared library. Returns 0 on successful loading, 1 in
/// case lib was already loaded and -1 in case lib does not exist
/// or in case of error.

int TWinNTSystem::Load(const char *module, const char *entry, Bool_t system)
{
   return TSystem::Load(module, entry, system);
}

/* nonstandard extension used : zero-sized array in struct/union */
#pragma warning(push)
#pragma warning(disable:4200)
////////////////////////////////////////////////////////////////////////////////
/// Get list of shared libraries loaded at the start of the executable.
/// Returns 0 in case list cannot be obtained or in case of error.

const char *TWinNTSystem::GetLinkedLibraries()
{
   char winDrive[256];
   char winDir[256];
   char winName[256];
   char winExt[256];

   if (!gApplication) return nullptr;

   static Bool_t once = kFALSE;
   static TString linkedLibs;

   if (!linkedLibs.IsNull())
      return linkedLibs;

   if (once)
      return nullptr;

   char *exe = gSystem->Which(Getenv("PATH"), gApplication->Argv(0),
                              kExecutePermission);
   if (!exe) {
      once = kTRUE;
      return nullptr;
   }

   HANDLE hFile, hMapping;
   void *basepointer;

   if((hFile = CreateFile(exe,GENERIC_READ,FILE_SHARE_READ,0,OPEN_EXISTING,FILE_FLAG_SEQUENTIAL_SCAN,0))==INVALID_HANDLE_VALUE) {
      delete [] exe;
      return nullptr;
   }
   if(!(hMapping = CreateFileMapping(hFile,0,PAGE_READONLY|SEC_COMMIT,0,0,0))) {
      CloseHandle(hFile);
      delete [] exe;
      return nullptr;
   }
   if(!(basepointer = MapViewOfFile(hMapping,FILE_MAP_READ,0,0,0))) {
      CloseHandle(hMapping);
      CloseHandle(hFile);
      delete [] exe;
      return nullptr;
   }

   int sect;
   IMAGE_DOS_HEADER *dos_head = (IMAGE_DOS_HEADER *)basepointer;
   struct header {
      DWORD signature;
      IMAGE_FILE_HEADER _head;
      IMAGE_OPTIONAL_HEADER opt_head;
      IMAGE_SECTION_HEADER section_header[];  // actual number in NumberOfSections
   };
   struct header *pheader;
   const IMAGE_SECTION_HEADER * section_header;

   if(dos_head->e_magic!='ZM') {
      delete [] exe;
      return nullptr;
   }  // verify DOS-EXE-Header
   // after end of DOS-EXE-Header: offset to PE-Header
   pheader = (struct header *)((char*)dos_head + dos_head->e_lfanew);

   if(IsBadReadPtr(pheader,sizeof(struct header))) { // start of PE-Header
      delete [] exe;
      return nullptr;
   }
   if(pheader->signature!=IMAGE_NT_SIGNATURE) {      // verify PE format
      switch((unsigned short)pheader->signature) {
         case IMAGE_DOS_SIGNATURE:
            delete [] exe;
            return nullptr;
         case IMAGE_OS2_SIGNATURE:
            delete [] exe;
            return nullptr;
         case IMAGE_OS2_SIGNATURE_LE:
            delete [] exe;
            return nullptr;
         default: // unknown signature
            delete [] exe;
            return nullptr;
      }
   }
#define isin(address,start,length) ((address)>=(start) && (address)<(start)+(length))
   TString odump;
   // walk through sections
   for(sect=0,section_header=pheader->section_header;
       sect<pheader->_head.NumberOfSections;sect++,section_header++) {
      int directory;
      const void * const section_data =
            (char*)basepointer + section_header->PointerToRawData;
      for(directory=0;directory<IMAGE_NUMBEROF_DIRECTORY_ENTRIES;directory++) {
         if(isin(pheader->opt_head.DataDirectory[directory].VirtualAddress,
                 section_header->VirtualAddress,
                 section_header->SizeOfRawData)) {
            const IMAGE_IMPORT_DESCRIPTOR *stuff_start =
                 (IMAGE_IMPORT_DESCRIPTOR *)((char*)section_data +
                 (pheader->opt_head.DataDirectory[directory].VirtualAddress -
                  section_header->VirtualAddress));
            // (virtual address of stuff - virtual address of section) =
            // offset of stuff in section
            const unsigned stuff_length =
                  pheader->opt_head.DataDirectory[directory].Size;
            if(directory == IMAGE_DIRECTORY_ENTRY_IMPORT) {
               while(!IsBadReadPtr(stuff_start,sizeof(*stuff_start)) &&
                      stuff_start->Name) {
                  TString dll = (char*)section_data +
                               ((DWORD)(stuff_start->Name)) -
                                section_header->VirtualAddress;
                  if (dll.EndsWith(".dll")) {
                     char *dllPath = DynamicPathName(dll, kTRUE);
                     if (dllPath) {
                        char *winPath = getenv("windir");
                        _splitpath(winPath,winDrive,winDir,winName,winExt);
                        if(!strstr(dllPath, winDir)) {
                           if (!linkedLibs.IsNull())
                              linkedLibs += " ";
                           linkedLibs += dllPath;
                        }
                     }
                     delete [] dllPath;
                  }
                  stuff_start++;
               }
            }
         }
      }
   }

   UnmapViewOfFile(basepointer);
   CloseHandle(hMapping);
   CloseHandle(hFile);

   delete [] exe;

   once = kTRUE;

   if (linkedLibs.IsNull())
      return nullptr;

   return linkedLibs;
}
#pragma warning(pop)

////////////////////////////////////////////////////////////////////////////////
/// Return a space separated list of loaded shared libraries.
/// This list is of a format suitable for a linker, i.e it may contain
/// -Lpathname and/or -lNameOfLib.
/// Option can be any of:
///   S: shared libraries loaded at the start of the executable, because
///      they were specified on the link line.
///   D: shared libraries dynamically loaded after the start of the program.
///   L: list the .LIB rather than the .DLL (this is intended for linking)
///      [This options is not the default]

const char *TWinNTSystem::GetLibraries(const char *regexp, const char *options,
                                       Bool_t isRegexp)
{
   TString ntlibs;
   struct _stat buf;
   std::string str;
   char drive[_MAX_DRIVE], dir[_MAX_DIR], fname[_MAX_FNAME], ext[_MAX_EXT];
   TString libs(TSystem::GetLibraries(regexp, options, isRegexp));
   TString opt = options;
   std::vector<std::string> all_libs, libpaths;

   if ( (opt.First('L')!=kNPOS) ) {
      libs.ReplaceAll("/","\\");
      // get the %LIB% environment path list
      std::stringstream libenv(gSystem->Getenv("LIB"));
      while (getline(libenv, str, ';')) {
         libpaths.push_back(str);
      }
      // now get the list of libraries
      std::stringstream libraries(libs.Data());
      while (getline(libraries, str, ' ')) {
         std::string::size_type first, last;
         // if the line begins with "-L", it's a linker option
         // (e.g. -LIBPATH:%ROOTSYS%\\lib), so add it to the path list
         if (str.rfind("-L", 0) == 0) {
            first = str.find_first_of('%');
            last = str.find_last_of('%');
            if ((first != std::string::npos) && (last != std::string::npos) &&
                (first != last)) {
               // if there is a string between %%, this is an environment
               // variable (e.g. %ROOTSYS%), so let's try to resolve it
               // and replace it with the real path
               std::string var = str.substr(first+1, last-first-1);
               std::string env(gSystem->Getenv(var.c_str()));
               if (!env.empty()) {
                  // the environment variable exist and properly resolved
                  // so add the last part of the path and add it to the list
                  env += str.substr(last+1);
                  libpaths.push_back(env);
               }
            }
            // keep the linker instuction in the final list
            ntlibs.Append(str.c_str());
            ntlibs += " ";
            continue;
         }
         // replace the '.dll' or '.DLL' extension by '.lib'
         last = str.rfind(".dll");
         if (last != std::string::npos)
            str.replace(last, 4, ".lib");
         last = str.rfind(".DLL");
         if (last != std::string::npos)
            str.replace(last, 4, ".lib");
         if (str.rfind(".lib") != std::string::npos ||
             str.rfind(".LIB") != std::string::npos) {
            // check if the .lib with its full path exists
            if (_stat( str.c_str(), &buf ) == 0) {
               // file exists, so keep it with full path in our final list
               ntlibs.Append(str.c_str());
               ntlibs += " ";
               continue;
            }
         }
         // full path not found, so split it to extract the library name
         // only, set its extension to '.lib' and add it to the list of
         // libraries to search
         _splitpath(str.c_str(), drive, dir, fname, ext);
         std::string libname(fname);
         libname += ".lib";
         all_libs.push_back(libname);
      }
      for (auto lib : all_libs) {
         // loop over all libraries to check which one exists
         for (auto libpath : libpaths) {
            // check in each path of the %LIB% environment
            std::string path_lib(libpath);
            path_lib += "\\";
            path_lib += lib;
            if (_stat( path_lib.c_str(), &buf ) == 0) {
                // file exists, add it to the final list of libraries
               ntlibs.Append(lib.c_str());
               ntlibs += " ";
            }
         }
      }
   } else {
      ntlibs = libs;
   }

   fListLibs = ntlibs;
   fListLibs.ReplaceAll("/","\\");
   return fListLibs;
}

//---- Time & Date -------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Add timer to list of system timers.

void TWinNTSystem::AddTimer(TTimer *ti)
{
   TSystem::AddTimer(ti);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove timer from list of system timers.

TTimer *TWinNTSystem::RemoveTimer(TTimer *ti)
{
   if (!ti) return nullptr;

   TTimer *t = TSystem::RemoveTimer(ti);
   return t;
}

////////////////////////////////////////////////////////////////////////////////
/// Special Thread to check asynchronous timers.

void TWinNTSystem::TimerThread()
{
   while (1) {
      if (!fInsideNotify)
         DispatchTimers(kFALSE);
      ::Sleep(kItimerResolution/2);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle and dispatch timers. If mode = kTRUE dispatch synchronous
/// timers else a-synchronous timers.

Bool_t TWinNTSystem::DispatchTimers(Bool_t mode)
{
   if (!fTimers) return kFALSE;

   fInsideNotify = kTRUE;

   TOrdCollectionIter it((TOrdCollection*)fTimers);
   TTimer *t;
   Bool_t  timedout = kFALSE;

   while ((t = (TTimer *) it.Next())) {
      // NB: the timer resolution is added in TTimer::CheckTimer()
      TTime now = Now();
      if (mode && t->IsSync()) {
         if (t->CheckTimer(now)) {
            timedout = kTRUE;
         }
      } else if (!mode && t->IsAsync()) {
         if (t->CheckTimer(now)) {
            timedout = kTRUE;
         }
      }
   }
   fInsideNotify = kFALSE;

   return timedout;
}

const Double_t gTicks = 1.0e-7;
////////////////////////////////////////////////////////////////////////////////
///

Double_t TWinNTSystem::GetRealTime()
{
   union {
      FILETIME ftFileTime;
      __int64  ftInt64;
   } ftRealTime; // time the process has spent in kernel mode

   ::GetSystemTimeAsFileTime(&ftRealTime.ftFileTime);
   return (Double_t)ftRealTime.ftInt64 * gTicks;
}

////////////////////////////////////////////////////////////////////////////////
///

Double_t TWinNTSystem::GetCPUTime()
{
   OSVERSIONINFO OsVersionInfo;

//*-*         Value                      Platform
//*-*  ----------------------------------------------------
//*-*  VER_PLATFORM_WIN32s              Win32s on Windows 3.1
//*-*  VER_PLATFORM_WIN32_WINDOWS       Win32 on Windows 95
//*-*  VER_PLATFORM_WIN32_NT            Windows NT
//*-*

   OsVersionInfo.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
   GetVersionEx(&OsVersionInfo);
   if (OsVersionInfo.dwPlatformId == VER_PLATFORM_WIN32_NT) {
      DWORD       ret;
      FILETIME    ftCreate,       // when the process was created
                  ftExit;         // when the process exited

      union {
         FILETIME ftFileTime;
         __int64  ftInt64;
      } ftKernel; // time the process has spent in kernel mode

      union {
         FILETIME ftFileTime;
         __int64  ftInt64;
      } ftUser;   // time the process has spent in user mode

      HANDLE hThread = GetCurrentThread();
      ret = GetThreadTimes (hThread, &ftCreate, &ftExit,
                                     &ftKernel.ftFileTime,
                                     &ftUser.ftFileTime);
      if (ret != TRUE){
         ret = ::GetLastError();
         ::Error("GetCPUTime", " Error on GetProcessTimes 0x%lx", (int)ret);
      }

      // Process times are returned in a 64-bit structure, as the number of
      // 100 nanosecond ticks since 1 January 1601.  User mode and kernel mode
      // times for this process are in separate 64-bit structures.
      // To convert to floating point seconds, we will:
      //          Convert sum of high 32-bit quantities to 64-bit int

       return (Double_t) (ftKernel.ftInt64 + ftUser.ftInt64) * gTicks;
   } else {
      return GetRealTime();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get current time in milliseconds since 0:00 Jan 1 1995.

TTime TWinNTSystem::Now()
{
   static time_t jan95 = 0;
   if (!jan95) {
      struct tm tp;
      tp.tm_year  = 95;
      tp.tm_mon   = 0;
      tp.tm_mday  = 1;
      tp.tm_hour  = 0;
      tp.tm_min   = 0;
      tp.tm_sec   = 0;
      tp.tm_isdst = -1;

      jan95 = mktime(&tp);
      if ((int)jan95 == -1) {
         ::SysError("TWinNTSystem::Now", "error converting 950001 0:00 to time_t");
         return 0;
      }
   }

   _timeb now;
   _ftime(&now);
   return TTime((now.time-(Long_t)jan95)*1000 + now.millitm);
}

////////////////////////////////////////////////////////////////////////////////
/// Sleep milliSec milli seconds.
/// The Sleep function suspends the execution of the CURRENT THREAD for
/// a specified interval.

void TWinNTSystem::Sleep(UInt_t milliSec)
{
   ::Sleep(milliSec);
}

////////////////////////////////////////////////////////////////////////////////
/// Select on file descriptors. The timeout to is in millisec.

Int_t TWinNTSystem::Select(TList *act, Long_t to)
{
   Int_t rc = -4;

   TFdSet rd, wr;
   Int_t mxfd = -1;
   TIter next(act);
   TFileHandler *h = nullptr;
   while ((h = (TFileHandler *) next())) {
      Int_t fd = h->GetFd();
      if (h->HasReadInterest())
         rd.Set(fd);
      if (h->HasWriteInterest())
         wr.Set(fd);
      h->ResetReadyMask();
   }
   rc = WinNTSelect(&rd, &wr, to);

   // Set readiness bits
   if (rc > 0) {
      next.Reset();
      while ((h = (TFileHandler *) next())) {
         Int_t fd = h->GetFd();
         if (rd.IsSet(fd))
            h->SetReadReady();
         if (wr.IsSet(fd))
            h->SetWriteReady();
      }
   }

   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Select on the file descriptor related to file handler h.
/// The timeout to is in millisec.

Int_t TWinNTSystem::Select(TFileHandler *h, Long_t to)
{
   Int_t rc = -4;

   TFdSet rd, wr;
   Int_t fd = -1;
   if (h) {
      fd = h->GetFd();
      if (h->HasReadInterest())
         rd.Set(fd);
      if (h->HasWriteInterest())
         wr.Set(fd);
      h->ResetReadyMask();
      rc = WinNTSelect(&rd, &wr, to);
   }

   // Fill output lists, if required
   if (rc > 0) {
      if (rd.IsSet(fd))
         h->SetReadReady();
      if (wr.IsSet(fd))
         h->SetWriteReady();
   }

   return rc;
}

//---- RPC ---------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
/// Get port # of internet service.

int TWinNTSystem::GetServiceByName(const char *servicename)
{
   struct servent *sp;

   if ((sp = ::getservbyname(servicename, kProtocolName)) == 0) {
      Error("GetServiceByName", "no service \"%s\" with protocol \"%s\"\n",
             servicename, kProtocolName);
      return -1;
   }
   return ::ntohs(sp->s_port);
}

////////////////////////////////////////////////////////////////////////////////

char *TWinNTSystem::GetServiceByPort(int port)
{
   // Get name of internet service.

   struct servent *sp;

   if ((sp = ::getservbyport(::htons(port), kProtocolName)) == 0) {
      return Form("%d", port);
   }
   return sp->s_name;
}

////////////////////////////////////////////////////////////////////////////////
/// Get Internet Protocol (IP) address of host.

TInetAddress TWinNTSystem::GetHostByName(const char *hostname)
{
   struct hostent *host_ptr;
   const char     *host;
   int             type;
   UInt_t          addr;    // good for 4 byte addresses

   if ((addr = ::inet_addr(hostname)) != INADDR_NONE) {
      type = AF_INET;
      if ((host_ptr = ::gethostbyaddr((const char *)&addr,
                                      sizeof(addr), AF_INET))) {
         host = host_ptr->h_name;
         TInetAddress a(host, ntohl(addr), type);
         UInt_t addr2;
         Int_t  i;
         for (i = 1; host_ptr->h_addr_list[i]; i++) {
            memcpy(&addr2, host_ptr->h_addr_list[i], host_ptr->h_length);
            a.AddAddress(ntohl(addr2));
         }
         for (i = 0; host_ptr->h_aliases[i]; i++)
            a.AddAlias(host_ptr->h_aliases[i]);
         return a;
      } else {
         host = "UnNamedHost";
      }
   } else if ((host_ptr = ::gethostbyname(hostname))) {
      // Check the address type for an internet host
      if (host_ptr->h_addrtype != AF_INET) {
         Error("GetHostByName", "%s is not an internet host\n", hostname);
         return TInetAddress();
      }
      memcpy(&addr, host_ptr->h_addr, host_ptr->h_length);
      host = host_ptr->h_name;
      type = host_ptr->h_addrtype;
      TInetAddress a(host, ntohl(addr), type);
      UInt_t addr2;
      Int_t  i;
      for (i = 1; host_ptr->h_addr_list[i]; i++) {
         memcpy(&addr2, host_ptr->h_addr_list[i], host_ptr->h_length);
         a.AddAddress(ntohl(addr2));
      }
      for (i = 0; host_ptr->h_aliases[i]; i++)
         a.AddAlias(host_ptr->h_aliases[i]);
      return a;
   } else {
      if (gDebug > 0) Error("GetHostByName", "unknown host %s", hostname);
      return TInetAddress(hostname, 0, -1);
   }

   return TInetAddress(host, ::ntohl(addr), type);
}

////////////////////////////////////////////////////////////////////////////////
/// Get Internet Protocol (IP) address of remote host and port #.

TInetAddress TWinNTSystem::GetPeerName(int socket)
{
   SOCKET sock = socket;
   struct sockaddr_in addr;
   int len = sizeof(addr);

   if (::getpeername(sock, (struct sockaddr *)&addr, &len) == SOCKET_ERROR) {
      ::SysError("GetPeerName", "getpeername");
      return TInetAddress();
   }

   struct hostent *host_ptr;
   const char *hostname;
   int         family;
   UInt_t      iaddr;

   if ((host_ptr = ::gethostbyaddr((const char *)&addr.sin_addr,
                                   sizeof(addr.sin_addr), AF_INET))) {
      memcpy(&iaddr, host_ptr->h_addr, host_ptr->h_length);
      hostname = host_ptr->h_name;
      family   = host_ptr->h_addrtype;
   } else {
      memcpy(&iaddr, &addr.sin_addr, sizeof(addr.sin_addr));
      hostname = "????";
      family   = AF_INET;
   }

   return TInetAddress(hostname, ::ntohl(iaddr), family, ::ntohs(addr.sin_port));
}

////////////////////////////////////////////////////////////////////////////////
/// Get Internet Protocol (IP) address of host and port #.

TInetAddress TWinNTSystem::GetSockName(int socket)
{
   SOCKET sock = socket;
   struct sockaddr_in addr;
   int len = sizeof(addr);

   if (::getsockname(sock, (struct sockaddr *)&addr, &len) == SOCKET_ERROR) {
      ::SysError("GetSockName", "getsockname");
      return TInetAddress();
   }

   struct hostent *host_ptr;
   const char *hostname;
   int         family;
   UInt_t      iaddr;

   if ((host_ptr = ::gethostbyaddr((const char *)&addr.sin_addr,
                                   sizeof(addr.sin_addr), AF_INET))) {
      memcpy(&iaddr, host_ptr->h_addr, host_ptr->h_length);
      hostname = host_ptr->h_name;
      family   = host_ptr->h_addrtype;
   } else {
      memcpy(&iaddr, &addr.sin_addr, sizeof(addr.sin_addr));
      hostname = "????";
      family   = AF_INET;
   }

   return TInetAddress(hostname, ::ntohl(iaddr), family, ::ntohs(addr.sin_port));
}

////////////////////////////////////////////////////////////////////////////////
/// Announce unix domain service.

int TWinNTSystem::AnnounceUnixService(int port, int backlog)
{
   SOCKET sock;

   // Create socket
   if ((sock = ::socket(AF_INET, SOCK_STREAM, 0)) == INVALID_SOCKET) {
      ::SysError("TWinNTSystem::AnnounceUnixService", "socket");
      return -1;
   }

   struct sockaddr_in inserver;
   memset(&inserver, 0, sizeof(inserver));
   inserver.sin_family = AF_INET;
   inserver.sin_addr.s_addr = ::htonl(INADDR_LOOPBACK);
   inserver.sin_port = port;

   // Bind socket
   if (port > 0) {
      if (::bind(sock, (struct sockaddr*) &inserver, sizeof(inserver)) == SOCKET_ERROR) {
         ::SysError("TWinNTSystem::AnnounceUnixService", "bind");
         return -2;
      }
   }
   // Start accepting connections
   if (::listen(sock, backlog)) {
      ::SysError("TWinNTSystem::AnnounceUnixService", "listen");
      return -1;
   }
   return (int)sock;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a socket on path 'sockpath', bind to it and start listening for Unix
/// domain connections to it. Returns socket fd or -1.

int TWinNTSystem::AnnounceUnixService(const char *sockpath, int backlog)
{
   if (!sockpath || strlen(sockpath) <= 0) {
      ::SysError("TWinNTSystem::AnnounceUnixService", "socket path undefined");
      return -1;
   }

   struct sockaddr_in myaddr;
   FILE * fp;
   int len = sizeof myaddr;
   int rc;
   int sock;

   // Create socket
   if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
      ::SysError("TWinNTSystem::AnnounceUnixService", "socket");
      return -1;
   }

   memset(&myaddr, 0, sizeof(myaddr));
   myaddr.sin_port = 0;
   myaddr.sin_family = AF_INET;
   myaddr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);

   rc = bind(sock, (struct sockaddr *)&myaddr, len);
   if (rc) {
      ::SysError("TWinNTSystem::AnnounceUnixService", "bind");
      return rc;
   }
   rc = getsockname(sock, (struct sockaddr *)&myaddr, &len);
   if (rc) {
      ::SysError("TWinNTSystem::AnnounceUnixService", "getsockname");
      return rc;
   }
   TString socketpath = sockpath;
   socketpath.ReplaceAll("/", "\\");
   fp = fopen(socketpath, "wb");
   if (!fp) {
      ::SysError("TWinNTSystem::AnnounceUnixService", "fopen");
      return -1;
   }
   fprintf(fp, "%d", myaddr.sin_port);
   fclose(fp);

   // Start accepting connections
   if (listen(sock, backlog)) {
      ::SysError("TWinNTSystem::AnnounceUnixService", "listen");
      return -1;
   }

   return sock;
}

////////////////////////////////////////////////////////////////////////////////
/// Close socket.

void TWinNTSystem::CloseConnection(int socket, Bool_t force)
{
   if (socket == -1) return;
   SOCKET sock = socket;

   if (force) {
      ::shutdown(sock, 2);
   }
   struct linger linger = {0, 0};
   ::setsockopt(sock, SOL_SOCKET, SO_LINGER, (char *) &linger, sizeof(linger));
   while (::closesocket(sock) == SOCKET_ERROR && WSAGetLastError() == WSAEINTR) {
      TSystem::ResetErrno();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Receive a buffer headed by a length indicator. Length is the size of
/// the buffer. Returns the number of bytes received in buf or -1 in
/// case of error.

int TWinNTSystem::RecvBuf(int sock, void *buf, int length)
{
   Int_t header;

   if (WinNTRecv(sock, &header, sizeof(header), 0) > 0) {
      int count = ::ntohl(header);

      if (count > length) {
         Error("RecvBuf", "record header exceeds buffer size");
         return -1;
      } else if (count > 0) {
         if (WinNTRecv(sock, buf, count, 0) < 0) {
            Error("RecvBuf", "cannot receive buffer");
            return -1;
         }
      }
      return count;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Send a buffer headed by a length indicator. Returns length of sent buffer
/// or -1 in case of error.

int TWinNTSystem::SendBuf(int sock, const void *buf, int length)
{
   Int_t header = ::htonl(length);

   if (WinNTSend(sock, &header, sizeof(header), 0) < 0) {
      Error("SendBuf", "cannot send header");
      return -1;
   }
   if (length > 0) {
      if (WinNTSend(sock, buf, length, 0) < 0) {
         Error("SendBuf", "cannot send buffer");
         return -1;
      }
   }
   return length;
}

////////////////////////////////////////////////////////////////////////////////
/// Receive exactly length bytes into buffer. Use opt to receive out-of-band
/// data or to have a peek at what is in the buffer (see TSocket). Buffer
/// must be able to store at least length bytes. Returns the number of
/// bytes received (can be 0 if other side of connection was closed) or -1
/// in case of error, -2 in case of MSG_OOB and errno == EWOULDBLOCK, -3
/// in case of MSG_OOB and errno == EINVAL and -4 in case of kNoBlock and
/// errno == EWOULDBLOCK. Returns -5 if pipe broken or reset by peer
/// (EPIPE || ECONNRESET).

int TWinNTSystem::RecvRaw(int sock, void *buf, int length, int opt)
{
   int flag;

   switch (opt) {
   case kDefault:
      flag = 0;
      break;
   case kOob:
      flag = MSG_OOB;
      break;
   case kPeek:
      flag = MSG_PEEK;
      break;
   case kDontBlock:
      flag = -1;
      break;
   default:
      flag = 0;
      break;
   }

   int n;
   if ((n = WinNTRecv(sock, buf, length, flag)) <= 0) {
      if (n == -1) {
         Error("RecvRaw", "cannot receive buffer");
      }
      return n;
   }
   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Send exactly length bytes from buffer. Use opt to send out-of-band
/// data (see TSocket). Returns the number of bytes sent or -1 in case of
/// error. Returns -4 in case of kNoBlock and errno == EWOULDBLOCK.
/// Returns -5 if pipe broken or reset by peer (EPIPE || ECONNRESET).

int TWinNTSystem::SendRaw(int sock, const void *buf, int length, int opt)
{
   int flag;

   switch (opt) {
   case kDefault:
      flag = 0;
      break;
   case kOob:
      flag = MSG_OOB;
      break;
   case kDontBlock:
      flag = -1;
      break;
   case kPeek:            // receive only option (see RecvRaw)
   default:
      flag = 0;
      break;
   }

   int n;
   if ((n = WinNTSend(sock, buf, length, flag)) <= 0) {
      if (n == -1 && GetErrno() != EINTR) {
         Error("SendRaw", "cannot send buffer");
      }
      return n;
   }
   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Set socket option.

int  TWinNTSystem::SetSockOpt(int socket, int opt, int value)
{
   u_long val = value;
   if (socket == -1) return -1;
   SOCKET sock = socket;

   switch (opt) {
   case kSendBuffer:
      if (::setsockopt(sock, SOL_SOCKET, SO_SNDBUF, (char*)&val, sizeof(val)) == SOCKET_ERROR) {
         ::SysError("SetSockOpt", "setsockopt(SO_SNDBUF)");
         return -1;
      }
      break;
   case kRecvBuffer:
      if (::setsockopt(sock, SOL_SOCKET, SO_RCVBUF, (char*)&val, sizeof(val)) == SOCKET_ERROR) {
         ::SysError("SetSockOpt", "setsockopt(SO_RCVBUF)");
         return -1;
      }
      break;
   case kOobInline:
      if (::setsockopt(sock, SOL_SOCKET, SO_OOBINLINE, (char*)&val, sizeof(val)) == SOCKET_ERROR) {
         SysError("SetSockOpt", "setsockopt(SO_OOBINLINE)");
         return -1;
      }
      break;
   case kKeepAlive:
      if (::setsockopt(sock, SOL_SOCKET, SO_KEEPALIVE, (char*)&val, sizeof(val)) == SOCKET_ERROR) {
         ::SysError("SetSockOpt", "setsockopt(SO_KEEPALIVE)");
         return -1;
      }
      break;
   case kReuseAddr:
      if (::setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, (char*)&val, sizeof(val)) == SOCKET_ERROR) {
         ::SysError("SetSockOpt", "setsockopt(SO_REUSEADDR)");
         return -1;
      }
      break;
   case kNoDelay:
      if (::setsockopt(sock, IPPROTO_TCP, TCP_NODELAY, (char*)&val, sizeof(val)) == SOCKET_ERROR) {
         ::SysError("SetSockOpt", "setsockopt(TCP_NODELAY)");
         return -1;
      }
      break;
   case kNoBlock:
      if (::ioctlsocket(sock, FIONBIO, &val) == SOCKET_ERROR) {
         ::SysError("SetSockOpt", "ioctl(FIONBIO)");
         return -1;
      }
      break;
#if 0
   case kProcessGroup:
      if (::ioctl(sock, SIOCSPGRP, &val) == -1) {
         ::SysError("SetSockOpt", "ioctl(SIOCSPGRP)");
         return -1;
      }
      break;
#endif
   case kAtMark:       // read-only option (see GetSockOpt)
   case kBytesToRead:  // read-only option
   default:
      Error("SetSockOpt", "illegal option (%d)", opt);
      return -1;
      break;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get socket option.

int TWinNTSystem::GetSockOpt(int socket, int opt, int *val)
{
   if (socket == -1) return -1;
   SOCKET sock = socket;

   int optlen = sizeof(*val);

   switch (opt) {
   case kSendBuffer:
      if (::getsockopt(sock, SOL_SOCKET, SO_SNDBUF, (char*)val, &optlen) == SOCKET_ERROR) {
         ::SysError("GetSockOpt", "getsockopt(SO_SNDBUF)");
         return -1;
      }
      break;
   case kRecvBuffer:
      if (::getsockopt(sock, SOL_SOCKET, SO_RCVBUF, (char*)val, &optlen) == SOCKET_ERROR) {
         ::SysError("GetSockOpt", "getsockopt(SO_RCVBUF)");
         return -1;
      }
      break;
   case kOobInline:
      if (::getsockopt(sock, SOL_SOCKET, SO_OOBINLINE, (char*)val, &optlen) == SOCKET_ERROR) {
         ::SysError("GetSockOpt", "getsockopt(SO_OOBINLINE)");
         return -1;
      }
      break;
   case kKeepAlive:
      if (::getsockopt(sock, SOL_SOCKET, SO_KEEPALIVE, (char*)val, &optlen) == SOCKET_ERROR) {
         ::SysError("GetSockOpt", "getsockopt(SO_KEEPALIVE)");
         return -1;
      }
      break;
   case kReuseAddr:
      if (::getsockopt(sock, SOL_SOCKET, SO_REUSEADDR, (char*)val, &optlen) == SOCKET_ERROR) {
         ::SysError("GetSockOpt", "getsockopt(SO_REUSEADDR)");
         return -1;
      }
      break;
   case kNoDelay:
      if (::getsockopt(sock, IPPROTO_TCP, TCP_NODELAY, (char*)val, &optlen) == SOCKET_ERROR) {
         ::SysError("GetSockOpt", "getsockopt(TCP_NODELAY)");
         return -1;
      }
      break;
   case kNoBlock:
      {
         int flg = 0;
         if (sock == INVALID_SOCKET) {
            ::SysError("GetSockOpt", "INVALID_SOCKET");
         }
         *val = flg; //  & O_NDELAY;  It is not been defined for WIN32
         return -1;
      }
      break;
#if 0
   case kProcessGroup:
      if (::ioctlsocket(sock, SIOCGPGRP, (u_long*)val) == SOCKET_ERROR) {
         ::SysError("GetSockOpt", "ioctl(SIOCGPGRP)");
         return -1;
      }
      break;
#endif
   case kAtMark:
      if (::ioctlsocket(sock, SIOCATMARK, (u_long*)val) == SOCKET_ERROR) {
         ::SysError("GetSockOpt", "ioctl(SIOCATMARK)");
         return -1;
      }
      break;
   case kBytesToRead:
      if (::ioctlsocket(sock, FIONREAD, (u_long*)val) == SOCKET_ERROR) {
         ::SysError("GetSockOpt", "ioctl(FIONREAD)");
         return -1;
      }
      break;
   default:
      Error("GetSockOpt", "illegal option (%d)", opt);
      *val = 0;
      return -1;
      break;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Connect to service servicename on server servername.

int TWinNTSystem::ConnectService(const char *servername, int port,
                                 int tcpwindowsize, const char *protocol)
{
   short  sport;
   struct servent *sp;

   if (!strcmp(servername, "unix")) {
      return WinNTUnixConnect(port);
   }
   else if (!gSystem->AccessPathName(servername) || servername[0] == '/' ||
            (servername[1] == ':' && servername[2] == '/')) {
      return WinNTUnixConnect(servername);
   }

   if (!strcmp(protocol, "udp")){
      return WinNTUdpConnect(servername, port);
   }

   if ((sp = ::getservbyport(::htons(port), kProtocolName))) {
      sport = sp->s_port;
   } else {
      sport = ::htons(port);
   }

   TInetAddress addr = gSystem->GetHostByName(servername);
   if (!addr.IsValid()) return -1;
   UInt_t adr = ::htonl(addr.GetAddress());

   struct sockaddr_in server;
   memset(&server, 0, sizeof(server));
   memcpy(&server.sin_addr, &adr, sizeof(adr));
   server.sin_family = addr.GetFamily();
   server.sin_port   = sport;

   // Create socket
   SOCKET sock;
   if ((sock = ::socket(AF_INET, SOCK_STREAM, 0)) == INVALID_SOCKET) {
      ::SysError("TWinNTSystem::WinNTConnectTcp", "socket");
      return -1;
   }

   if (tcpwindowsize > 0) {
      gSystem->SetSockOpt((int)sock, kRecvBuffer, tcpwindowsize);
      gSystem->SetSockOpt((int)sock, kSendBuffer, tcpwindowsize);
   }

   if (::connect(sock, (struct sockaddr*) &server, sizeof(server)) == INVALID_SOCKET) {
      //::SysError("TWinNTSystem::UnixConnectTcp", "connect");
      ::closesocket(sock);
      return -1;
   }
   return (int) sock;
}

////////////////////////////////////////////////////////////////////////////////
/// Connect to a Unix domain socket.

int TWinNTSystem::WinNTUnixConnect(int port)
{
   struct sockaddr_in myaddr;
   int sock;

   memset(&myaddr, 0, sizeof(myaddr));
   myaddr.sin_family = AF_INET;
   myaddr.sin_port = port;
   myaddr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);

   // Open socket
   if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
      ::SysError("TWinNTSystem::WinNTUnixConnect", "socket");
      return -1;
   }

   while ((connect(sock, (struct sockaddr *)&myaddr, sizeof myaddr)) == -1) {
      if (GetErrno() == EINTR)
         ResetErrno();
      else {
         ::SysError("TWinNTSystem::WinNTUnixConnect", "connect");
         close(sock);
         return -1;
      }
   }
   return sock;
}

////////////////////////////////////////////////////////////////////////////////
/// Connect to a Unix domain socket. Returns -1 in case of error.

int TWinNTSystem::WinNTUnixConnect(const char *sockpath)
{
   FILE *fp;
   int port = 0;

   if (!sockpath || strlen(sockpath) <= 0) {
      ::SysError("TWinNTSystem::WinNTUnixConnect", "socket path undefined");
      return -1;
   }
   TString socketpath = sockpath;
   socketpath.ReplaceAll("/", "\\");
   fp = fopen(socketpath.Data(), "rb");
   if (!fp) {
      ::SysError("TWinNTSystem::WinNTUnixConnect", "fopen");
      return -1;
   }
   fscanf(fp, "%d", &port);
   fclose(fp);
   /* XXX: set errno in this case */
   if (port < 0 || port > 65535) {
      ::SysError("TWinNTSystem::WinNTUnixConnect", "invalid port");
      return -1;
   }
   return WinNTUnixConnect(port);
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a UDP socket connection
/// Is called via the TSocket constructor. Returns -1 in case of error.

int TWinNTSystem::WinNTUdpConnect(const char *hostname, int port)
{
   short  sport;
   struct servent *sp;

   if ((sp = getservbyport(htons(port), kProtocolName)))
      sport = sp->s_port;
   else
      sport = htons(port);

   TInetAddress addr = gSystem->GetHostByName(hostname);
   if (!addr.IsValid()) return -1;
   UInt_t adr = htonl(addr.GetAddress());

   struct sockaddr_in server;
   memset(&server, 0, sizeof(server));
   memcpy(&server.sin_addr, &adr, sizeof(adr));
   server.sin_family = addr.GetFamily();
   server.sin_port   = sport;

   // Create socket
   int sock;
   if ((sock = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP)) < 0) {
      ::SysError("TWinNTSystem::WinNTUdpConnect", "socket (%s:%d)",
                 hostname, port);
      return -1;
   }

   while (connect(sock, (struct sockaddr*) &server, sizeof(server)) == -1) {
      if (GetErrno() == EINTR)
         ResetErrno();
      else {
         ::SysError("TWinNTSystem::WinNTUdpConnect", "connect (%s:%d)",
                    hostname, port);
         close(sock);
         return -1;
      }
   }
   return sock;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a connection to a service on a server. Returns -1 in case
/// connection cannot be opened.
/// Use tcpwindowsize to specify the size of the receive buffer, it has
/// to be specified here to make sure the window scale option is set (for
/// tcpwindowsize > 65KB and for platforms supporting window scaling).
/// Is called via the TSocket constructor.

int TWinNTSystem::OpenConnection(const char *server, int port, int tcpwindowsize,
                                 const char *protocol)
{
   return ConnectService(server, port, tcpwindowsize, protocol);
}

////////////////////////////////////////////////////////////////////////////////
/// Announce TCP/IP service.
/// Open a socket, bind to it and start listening for TCP/IP connections
/// on the port. If reuse is true reuse the address, backlog specifies
/// how many sockets can be waiting to be accepted.
/// Use tcpwindowsize to specify the size of the receive buffer, it has
/// to be specified here to make sure the window scale option is set (for
/// tcpwindowsize > 65KB and for platforms supporting window scaling).
/// Returns socket fd or -1 if socket() failed, -2 if bind() failed
/// or -3 if listen() failed.

int TWinNTSystem::AnnounceTcpService(int port, Bool_t reuse, int backlog,
                                     int tcpwindowsize)
{
   short  sport;
   struct servent *sp;
   const short kSOCKET_MINPORT = 5000, kSOCKET_MAXPORT = 15000;
   short tryport = kSOCKET_MINPORT;

   if ((sp = ::getservbyport(::htons(port), kProtocolName))) {
      sport = sp->s_port;
   } else {
      sport = ::htons(port);
   }

   if (port == 0 && reuse) {
      ::Error("TWinNTSystem::WinNTTcpService", "cannot do a port scan while reuse is true");
      return -1;
   }

   if ((sp = ::getservbyport(::htons(port), kProtocolName))) {
      sport = sp->s_port;
   } else {
      sport = ::htons(port);
   }

   // Create tcp socket
   SOCKET sock;
   if ((sock = ::socket(AF_INET, SOCK_STREAM, 0)) < 0) {
      ::SysError("TWinNTSystem::WinNTTcpService", "socket");
      return -1;
   }

   if (reuse) {
      gSystem->SetSockOpt((int)sock, kReuseAddr, 1);
   }

   if (tcpwindowsize > 0) {
      gSystem->SetSockOpt((int)sock, kRecvBuffer, tcpwindowsize);
      gSystem->SetSockOpt((int)sock, kSendBuffer, tcpwindowsize);
   }

   struct sockaddr_in inserver;
   memset(&inserver, 0, sizeof(inserver));
   inserver.sin_family = AF_INET;
   inserver.sin_addr.s_addr = ::htonl(INADDR_ANY);
   inserver.sin_port = sport;

   // Bind socket
   if (port > 0) {
      if (::bind(sock, (struct sockaddr*) &inserver, sizeof(inserver)) == SOCKET_ERROR) {
         ::SysError("TWinNTSystem::WinNTTcpService", "bind");
         return -2;
      }
   } else {
      int bret;
      do {
         inserver.sin_port = ::htons(tryport);
         bret = ::bind(sock, (struct sockaddr*) &inserver, sizeof(inserver));
         tryport++;
      } while (bret == SOCKET_ERROR && WSAGetLastError() == WSAEADDRINUSE &&
               tryport < kSOCKET_MAXPORT);
      if (bret == SOCKET_ERROR) {
         ::SysError("TWinNTSystem::WinNTTcpService", "bind (port scan)");
         return -2;
      }
   }

   // Start accepting connections
   if (::listen(sock, backlog) == SOCKET_ERROR) {
      ::SysError("TWinNTSystem::WinNTTcpService", "listen");
      return -3;
   }
   return (int)sock;
}

////////////////////////////////////////////////////////////////////////////////
/// Announce UDP service.

int TWinNTSystem::AnnounceUdpService(int port, int backlog)
{
   // Open a socket, bind to it and start listening for UDP connections
   // on the port. If reuse is true reuse the address, backlog specifies
   // how many sockets can be waiting to be accepted. If port is 0 a port
   // scan will be done to find a free port. This option is mutual exlusive
   // with the reuse option.

   const short kSOCKET_MINPORT = 5000, kSOCKET_MAXPORT = 15000;
   short  sport, tryport = kSOCKET_MINPORT;
   struct servent *sp;

   if ((sp = getservbyport(htons(port), kProtocolName)))
      sport = sp->s_port;
   else
      sport = htons(port);

   // Create udp socket
   int sock;
   if ((sock = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP)) < 0) {
      ::SysError("TUnixSystem::UnixUdpService", "socket");
      return -1;
   }

   struct sockaddr_in inserver;
   memset(&inserver, 0, sizeof(inserver));
   inserver.sin_family = AF_INET;
   inserver.sin_addr.s_addr = htonl(INADDR_ANY);
   inserver.sin_port = sport;

   // Bind socket
   if (port > 0) {
      if (bind(sock, (struct sockaddr*) &inserver, sizeof(inserver))) {
         ::SysError("TWinNTSystem::AnnounceUdpService", "bind");
         return -2;
      }
   } else {
      int bret;
      do {
         inserver.sin_port = htons(tryport);
         bret = bind(sock, (struct sockaddr*) &inserver, sizeof(inserver));
         tryport++;
      } while (bret == SOCKET_ERROR && WSAGetLastError() == WSAEADDRINUSE &&
               tryport < kSOCKET_MAXPORT);
      if (bret < 0) {
         ::SysError("TWinNTSystem::AnnounceUdpService", "bind (port scan)");
         return -2;
      }
   }

   // Start accepting connections
   if (listen(sock, backlog)) {
      ::SysError("TWinNTSystem::AnnounceUdpService", "listen");
      return -3;
   }

   return sock;
}

////////////////////////////////////////////////////////////////////////////////
/// Accept a connection. In case of an error return -1. In case
/// non-blocking I/O is enabled and no connections are available
/// return -2.

int TWinNTSystem::AcceptConnection(int socket)
{
   int soc = -1;
   SOCKET sock = socket;

   while ((soc = ::accept(sock, 0, 0)) == INVALID_SOCKET &&
          (::WSAGetLastError() == WSAEINTR)) {
      TSystem::ResetErrno();
   }

   if (soc == -1) {
      if (::WSAGetLastError() == WSAEWOULDBLOCK) {
         return -2;
      } else {
         ::SysError("AcceptConnection", "accept");
         return -1;
      }
   }
   return soc;
}

//---- System, CPU and Memory info ---------------------------------------------

// !!! using undocumented functions and structures !!!

#define SystemBasicInformation         0
#define SystemPerformanceInformation   2

typedef struct
{
   DWORD dwUnknown1;
   ULONG uKeMaximumIncrement;
   ULONG uPageSize;
   ULONG uMmNumberOfPhysicalPages;
   ULONG uMmLowestPhysicalPage;
   ULONG UMmHighestPhysicalPage;
   ULONG uAllocationGranularity;
   PVOID pLowestUserAddress;
   PVOID pMmHighestUserAddress;
   ULONG uKeActiveProcessors;
   BYTE  bKeNumberProcessors;
   BYTE  bUnknown2;
   WORD  bUnknown3;
} SYSTEM_BASIC_INFORMATION;

typedef struct
{
   LARGE_INTEGER  liIdleTime;
   DWORD    dwSpare[76];
} SYSTEM_PERFORMANCE_INFORMATION;

typedef struct _PROCESS_MEMORY_COUNTERS {
   DWORD cb;
   DWORD PageFaultCount;
   SIZE_T PeakWorkingSetSize;
   SIZE_T WorkingSetSize;
   SIZE_T QuotaPeakPagedPoolUsage;
   SIZE_T QuotaPagedPoolUsage;
   SIZE_T QuotaPeakNonPagedPoolUsage;
   SIZE_T QuotaNonPagedPoolUsage;
   SIZE_T PagefileUsage;
   SIZE_T PeakPagefileUsage;
} PROCESS_MEMORY_COUNTERS, *PPROCESS_MEMORY_COUNTERS;

typedef LONG (WINAPI *PROCNTQSI) (UINT, PVOID, ULONG, PULONG);

#define Li2Double(x) ((double)((x).HighPart) * 4.294967296E9 + (double)((x).LowPart))

////////////////////////////////////////////////////////////////////////////////
/// Calculate the CPU clock speed using the 'rdtsc' instruction.
/// RDTSC: Read Time Stamp Counter.

static DWORD GetCPUSpeed()
{
   LARGE_INTEGER ulFreq, ulTicks, ulValue, ulStartCounter;

   // Query for high-resolution counter frequency
   // (this is not the CPU frequency):
   if (QueryPerformanceFrequency(&ulFreq)) {
      // Query current value:
      QueryPerformanceCounter(&ulTicks);
      // Calculate end value (one second interval);
      // this is (current + frequency)
      ulValue.QuadPart = ulTicks.QuadPart + ulFreq.QuadPart/10;
      ulStartCounter.QuadPart = __rdtsc();

      // Loop for one second (measured with the high-resolution counter):
      do {
         QueryPerformanceCounter(&ulTicks);
      } while (ulTicks.QuadPart <= ulValue.QuadPart);
      // Now again read CPU time-stamp counter:
      return (DWORD)((__rdtsc() - ulStartCounter.QuadPart)/100000);
   } else {
      // No high-resolution counter present:
      return 0;
   }
}

#define BUFSIZE 80
#define SM_SERVERR2 89
typedef void (WINAPI *PGNSI)(LPSYSTEM_INFO);

////////////////////////////////////////////////////////////////////////////////

static const char *GetWindowsVersion()
{
   OSVERSIONINFOEX osvi;
   SYSTEM_INFO si;
   PGNSI pGNSI;
   BOOL bOsVersionInfoEx;
   static char *strReturn = nullptr;
   char temp[512];

   if (!strReturn)
      strReturn = new char[2048];
   else
      return strReturn;

   ZeroMemory(&si, sizeof(SYSTEM_INFO));
   ZeroMemory(&osvi, sizeof(OSVERSIONINFOEX));

   // Try calling GetVersionEx using the OSVERSIONINFOEX structure.
   // If that fails, try using the OSVERSIONINFO structure.

   osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEX);

   if( !(bOsVersionInfoEx = GetVersionEx ((OSVERSIONINFO *) &osvi)) )
   {
      osvi.dwOSVersionInfoSize = sizeof (OSVERSIONINFO);
      if (! GetVersionEx ( (OSVERSIONINFO *) &osvi) )
         return "";
   }

   // Call GetNativeSystemInfo if supported or GetSystemInfo otherwise.
   pGNSI = (PGNSI) GetProcAddress( GetModuleHandle("kernel32.dll"),
                                   "GetNativeSystemInfo");
   if(NULL != pGNSI)
      pGNSI(&si);
   else GetSystemInfo(&si);

   switch (osvi.dwPlatformId)
   {
      // Test for the Windows NT product family.
      case VER_PLATFORM_WIN32_NT:

         // Test for the specific product.
         if ( osvi.dwMajorVersion == 6 && osvi.dwMinorVersion == 0 )
         {
            if( osvi.wProductType == VER_NT_WORKSTATION )
                strlcpy(strReturn, "Microsoft Windows Vista ",2048);
            else strlcpy(strReturn, "Windows Server \"Longhorn\" " ,2048);
         }
         if ( osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 2 )
         {
            if( GetSystemMetrics(SM_SERVERR2) )
               strlcpy(strReturn, "Microsoft Windows Server 2003 \"R2\" ",2048);
            else if( osvi.wProductType == VER_NT_WORKSTATION &&
                      si.wProcessorArchitecture==PROCESSOR_ARCHITECTURE_AMD64)
            {
               strlcpy(strReturn, "Microsoft Windows XP Professional x64 Edition ",2048);
            }
            else strlcpy(strReturn, "Microsoft Windows Server 2003, ",2048);
         }
         if ( osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 1 )
            strlcpy(strReturn, "Microsoft Windows XP ",2048);

         if ( osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 0 )
            strlcpy(strReturn, "Microsoft Windows 2000 ",2048);

         if ( osvi.dwMajorVersion <= 4 )
            strlcpy(strReturn, "Microsoft Windows NT ",2048);

         // Test for specific product on Windows NT 4.0 SP6 and later.
         if( bOsVersionInfoEx )
         {
            // Test for the workstation type.
            if ( osvi.wProductType == VER_NT_WORKSTATION &&
                 si.wProcessorArchitecture!=PROCESSOR_ARCHITECTURE_AMD64)
            {
               if( osvi.dwMajorVersion == 4 )
                  strlcat(strReturn, "Workstation 4.0 ",2048 );
               else if( osvi.wSuiteMask & VER_SUITE_PERSONAL )
                  strlcat(strReturn, "Home Edition " ,2048);
               else strlcat(strReturn, "Professional " ,2048);
            }
            // Test for the server type.
            else if ( osvi.wProductType == VER_NT_SERVER ||
                      osvi.wProductType == VER_NT_DOMAIN_CONTROLLER )
            {
               if(osvi.dwMajorVersion==5 && osvi.dwMinorVersion==2)
               {
                  if ( si.wProcessorArchitecture==PROCESSOR_ARCHITECTURE_IA64 )
                  {
                      if( osvi.wSuiteMask & VER_SUITE_DATACENTER )
                         strlcat(strReturn, "Datacenter Edition for Itanium-based Systems",2048 );
                      else if( osvi.wSuiteMask & VER_SUITE_ENTERPRISE )
                         strlcat(strReturn, "Enterprise Edition for Itanium-based Systems" ,2048);
                  }
                  else if ( si.wProcessorArchitecture==PROCESSOR_ARCHITECTURE_AMD64 )
                  {
                      if( osvi.wSuiteMask & VER_SUITE_DATACENTER )
                         strlcat(strReturn, "Datacenter x64 Edition ",2048 );
                      else if( osvi.wSuiteMask & VER_SUITE_ENTERPRISE )
                         strlcat(strReturn, "Enterprise x64 Edition ",2048 );
                      else strlcat(strReturn, "Standard x64 Edition ",2048 );
                  }
                  else
                  {
                      if( osvi.wSuiteMask & VER_SUITE_DATACENTER )
                         strlcat(strReturn, "Datacenter Edition ",2048 );
                      else if( osvi.wSuiteMask & VER_SUITE_ENTERPRISE )
                         strlcat(strReturn, "Enterprise Edition ",2048 );
                      else if ( osvi.wSuiteMask == VER_SUITE_BLADE )
                         strlcat(strReturn, "Web Edition " ,2048);
                      else strlcat(strReturn, "Standard Edition ",2048 );
                  }
               }
               else if(osvi.dwMajorVersion==5 && osvi.dwMinorVersion==0)
               {
                  if( osvi.wSuiteMask & VER_SUITE_DATACENTER )
                     strlcat(strReturn, "Datacenter Server ",2048 );
                  else if( osvi.wSuiteMask & VER_SUITE_ENTERPRISE )
                     strlcat(strReturn, "Advanced Server ",2048 );
                  else strlcat(strReturn, "Server ",2048 );
               }
               else  // Windows NT 4.0
               {
                  if( osvi.wSuiteMask & VER_SUITE_ENTERPRISE )
                     strlcat(strReturn, "Server 4.0, Enterprise Edition " ,2048);
                  else strlcat(strReturn, "Server 4.0 ",2048 );
               }
            }
         }
         // Test for specific product on Windows NT 4.0 SP5 and earlier
         else
         {
            HKEY hKey;
            TCHAR szProductType[BUFSIZE];
            DWORD dwBufLen=BUFSIZE*sizeof(TCHAR);
            LONG lRet;

            lRet = RegOpenKeyEx( HKEY_LOCAL_MACHINE,
                                 "SYSTEM\\CurrentControlSet\\Control\\ProductOptions",
                                 0, KEY_QUERY_VALUE, &hKey );
            if( lRet != ERROR_SUCCESS )
               return "";

            lRet = RegQueryValueEx( hKey, "ProductType", NULL, NULL,
                                   (LPBYTE) szProductType, &dwBufLen);
            RegCloseKey( hKey );

            if( (lRet != ERROR_SUCCESS) || (dwBufLen > BUFSIZE*sizeof(TCHAR)) )
               return "";

            if ( lstrcmpi( "WINNT", szProductType) == 0 )
               strlcat(strReturn, "Workstation " ,2048);
            if ( lstrcmpi( "LANMANNT", szProductType) == 0 )
               strlcat(strReturn, "Server " ,2048);
            if ( lstrcmpi( "SERVERNT", szProductType) == 0 )
               strlcat(strReturn, "Advanced Server " ,2048);
            snprintf(temp,512, "%d.%d ", osvi.dwMajorVersion, osvi.dwMinorVersion);
            strlcat(strReturn, temp,2048);
         }

         // Display service pack (if any) and build number.

         if( osvi.dwMajorVersion == 4 &&
             lstrcmpi( osvi.szCSDVersion, "Service Pack 6" ) == 0 )
         {
            HKEY hKey;
            LONG lRet;

            // Test for SP6 versus SP6a.
            lRet = RegOpenKeyEx( HKEY_LOCAL_MACHINE,
                                 "SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion\\Hotfix\\Q246009",
                                 0, KEY_QUERY_VALUE, &hKey );
            if( lRet == ERROR_SUCCESS ) {
               snprintf(temp, 512, "Service Pack 6a (Build %d)", osvi.dwBuildNumber & 0xFFFF );
               strlcat(strReturn, temp,2048 );
            }
            else // Windows NT 4.0 prior to SP6a
            {
               snprintf(temp,512, "%s (Build %d)", osvi.szCSDVersion, osvi.dwBuildNumber & 0xFFFF);
               strlcat(strReturn, temp,2048 );
            }

            RegCloseKey( hKey );
         }
         else // not Windows NT 4.0
         {
            snprintf(temp, 512,"%s (Build %d)", osvi.szCSDVersion, osvi.dwBuildNumber & 0xFFFF);
            strlcat(strReturn, temp,2048 );
         }

         break;

      // Test for the Windows Me/98/95.
      case VER_PLATFORM_WIN32_WINDOWS:

         if (osvi.dwMajorVersion == 4 && osvi.dwMinorVersion == 0)
         {
             strlcpy(strReturn, "Microsoft Windows 95 ",2048);
             if (osvi.szCSDVersion[1]=='C' || osvi.szCSDVersion[1]=='B')
                strlcat(strReturn, "OSR2 " ,2048);
         }

         if (osvi.dwMajorVersion == 4 && osvi.dwMinorVersion == 10)
         {
             strlcpy(strReturn, "Microsoft Windows 98 ",2048);
             if ( osvi.szCSDVersion[1]=='A' || osvi.szCSDVersion[1]=='B')
                strlcat(strReturn, "SE ",2048 );
         }

         if (osvi.dwMajorVersion == 4 && osvi.dwMinorVersion == 90)
         {
             strlcpy(strReturn, "Microsoft Windows Millennium Edition",2048);
         }
         break;

      case VER_PLATFORM_WIN32s:
         strlcpy(strReturn, "Microsoft Win32s",2048);
         break;
   }
   return strReturn;
}

////////////////////////////////////////////////////////////////////////////////
/// Use assembly to retrieve the L2 cache information ...

static int GetL2CacheSize()
{
   unsigned nHighestFeatureEx;
   int nBuff[4];

   __cpuid(nBuff, 0x80000000);
   nHighestFeatureEx = (unsigned)nBuff[0];
   // Get cache size
   if (nHighestFeatureEx >= 0x80000006) {
      __cpuid(nBuff, 0x80000006);
      return (((unsigned)nBuff[2])>>16);
   }
   else return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get system info for Windows NT.

static void GetWinNTSysInfo(SysInfo_t *sysinfo)
{
   SYSTEM_PERFORMANCE_INFORMATION   SysPerfInfo;
   SYSTEM_INFO sysInfo;
   MEMORYSTATUSEX statex;
   OSVERSIONINFO OsVersionInfo;
   HKEY hKey;
   char  szKeyValueString[80];
   DWORD szKeyValueDword;
   DWORD dwBufLen;
   LONG  status;
   PROCNTQSI  NtQuerySystemInformation;

   NtQuerySystemInformation = (PROCNTQSI)GetProcAddress(
         GetModuleHandle("ntdll"), "NtQuerySystemInformation");

   if (!NtQuerySystemInformation) {
      ::Error("GetWinNTSysInfo",
              "Error on GetProcAddress(NtQuerySystemInformation)");
      return;
   }

   status = NtQuerySystemInformation(SystemPerformanceInformation,
                                     &SysPerfInfo, sizeof(SysPerfInfo),
                                     NULL);
   OsVersionInfo.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
   GetVersionEx(&OsVersionInfo);
   GetSystemInfo(&sysInfo);
   statex.dwLength = sizeof(statex);
   if (!GlobalMemoryStatusEx(&statex)) {
      ::Error("GetWinNTSysInfo", "Error on GlobalMemoryStatusEx()");
      return;
   }
   sysinfo->fCpus     = sysInfo.dwNumberOfProcessors;
   sysinfo->fPhysRam  = (Int_t)(statex.ullTotalPhys  >> 20);
   sysinfo->fOS       = GetWindowsVersion();
   sysinfo->fModel    = "";
   sysinfo->fCpuType  = "";
   sysinfo->fCpuSpeed = GetCPUSpeed();
   sysinfo->fBusSpeed = 0;  // bus speed in MHz
   sysinfo->fL2Cache  = GetL2CacheSize();

   status = RegOpenKeyEx(HKEY_LOCAL_MACHINE, "HARDWARE\\DESCRIPTION\\System",
                         0, KEY_QUERY_VALUE, &hKey);
   if (status == ERROR_SUCCESS) {
      dwBufLen = sizeof(szKeyValueString);
      RegQueryValueEx(hKey, "Identifier", NULL, NULL,(LPBYTE)szKeyValueString,
                      &dwBufLen);
      sysinfo->fModel = szKeyValueString;
      RegCloseKey (hKey);
   }
   status = RegOpenKeyEx(HKEY_LOCAL_MACHINE,
                         "Hardware\\Description\\System\\CentralProcessor\\0",
                         0, KEY_QUERY_VALUE, &hKey);
   if (status == ERROR_SUCCESS) {
      dwBufLen = sizeof(szKeyValueString);
      status = RegQueryValueEx(hKey, "ProcessorNameString", NULL, NULL,
                               (LPBYTE)szKeyValueString, &dwBufLen);
      if (status == ERROR_SUCCESS)
         sysinfo->fCpuType = szKeyValueString;
      dwBufLen = sizeof(DWORD);
      status = RegQueryValueEx(hKey,"~MHz",NULL,NULL,(LPBYTE)&szKeyValueDword,
                               &dwBufLen);
      if ((status == ERROR_SUCCESS) && ((sysinfo->fCpuSpeed <= 0) ||
         (sysinfo->fCpuSpeed < (szKeyValueDword >> 1))))
         sysinfo->fCpuSpeed = (Int_t)szKeyValueDword;
      RegCloseKey (hKey);
   }
   sysinfo->fCpuType.Remove(TString::kBoth, ' ');
   sysinfo->fModel.Remove(TString::kBoth, ' ');
}

////////////////////////////////////////////////////////////////////////////////
/// Get CPU stat for Window. Use sampleTime to set the interval over which
/// the CPU load will be measured, in ms (default 1000).

static void GetWinNTCpuInfo(CpuInfo_t *cpuinfo, Int_t sampleTime)
{
   SYSTEM_INFO sysInfo;
   Float_t  idle_ratio, kernel_ratio, user_ratio, total_ratio;
   FILETIME ft_sys_idle, ft_sys_kernel, ft_sys_user, ft_fun_time;
   SYSTEMTIME st_fun_time;

   ULARGE_INTEGER ul_sys_idle, ul_sys_kernel, ul_sys_user;
   static ULARGE_INTEGER ul_sys_idleold = {0, 0};
   static ULARGE_INTEGER ul_sys_kernelold = {0, 0};
   static ULARGE_INTEGER ul_sys_userold = {0, 0};
   ULARGE_INTEGER ul_sys_idle_diff, ul_sys_kernel_diff, ul_sys_user_diff;

   ULARGE_INTEGER ul_fun_time;
   ULARGE_INTEGER ul_fun_timeold = {0, 0};
   ULARGE_INTEGER ul_fun_time_diff;

   typedef BOOL (__stdcall *GetSystemTimesProc)( LPFILETIME lpIdleTime,
                 LPFILETIME lpKernelTime, LPFILETIME lpUserTime );
   static GetSystemTimesProc pGetSystemTimes = 0;

   HMODULE hModImagehlp = LoadLibrary( "Kernel32.dll" );
   if (!hModImagehlp) {
      ::Error("GetWinNTCpuInfo", "Error on LoadLibrary(Kernel32.dll)");
      return;
   }

   pGetSystemTimes = (GetSystemTimesProc) GetProcAddress( hModImagehlp,
                      "GetSystemTimes" );
   if (!pGetSystemTimes) {
      ::Error("GetWinNTCpuInfo", "Error on GetProcAddress(GetSystemTimes)");
      return;
   }
   GetSystemInfo(&sysInfo);

again:
   pGetSystemTimes(&ft_sys_idle,&ft_sys_kernel,&ft_sys_user);
   GetSystemTime(&st_fun_time);
   SystemTimeToFileTime(&st_fun_time,&ft_fun_time);

   memcpy(&ul_sys_idle, &ft_sys_idle, sizeof(FILETIME));
   memcpy(&ul_sys_kernel, &ft_sys_kernel, sizeof(FILETIME));
   memcpy(&ul_sys_user, &ft_sys_user, sizeof(FILETIME));
   memcpy(&ul_fun_time, &ft_fun_time, sizeof(FILETIME));

   ul_sys_idle_diff.QuadPart   = ul_sys_idle.QuadPart -
                                 ul_sys_idleold.QuadPart;
   ul_sys_kernel_diff.QuadPart = ul_sys_kernel.QuadPart -
                                 ul_sys_kernelold.QuadPart;
   ul_sys_user_diff.QuadPart   = ul_sys_user.QuadPart -
                                 ul_sys_userold.QuadPart;

   ul_fun_time_diff.QuadPart = ul_fun_time.QuadPart -
                               ul_fun_timeold.QuadPart;

   ul_sys_idleold.QuadPart   = ul_sys_idle.QuadPart;
   ul_sys_kernelold.QuadPart = ul_sys_kernel.QuadPart;
   ul_sys_userold.QuadPart   = ul_sys_user.QuadPart;

   if (ul_fun_timeold.QuadPart == 0) {
      Sleep(sampleTime);
      ul_fun_timeold.QuadPart = ul_fun_time.QuadPart;
      goto again;
   }
   ul_fun_timeold.QuadPart = ul_fun_time.QuadPart;

   idle_ratio = (Float_t)(Li2Double(ul_sys_idle_diff)/
                          Li2Double(ul_fun_time_diff))*100.0;
   user_ratio = (Float_t)(Li2Double(ul_sys_user_diff)/
                          Li2Double(ul_fun_time_diff))*100.0;
   kernel_ratio = (Float_t)(Li2Double(ul_sys_kernel_diff)/
                            Li2Double(ul_fun_time_diff))*100.0;
   idle_ratio /= (Float_t)sysInfo.dwNumberOfProcessors;
   user_ratio /= (Float_t)sysInfo.dwNumberOfProcessors;
   kernel_ratio /= (Float_t)sysInfo.dwNumberOfProcessors;
   total_ratio = 100.0 - idle_ratio;

   cpuinfo->fLoad1m  = 0; // cpu load average over 1 m
   cpuinfo->fLoad5m  = 0; // cpu load average over 5 m
   cpuinfo->fLoad15m = 0; // cpu load average over 15 m
   cpuinfo->fUser    = user_ratio; // cpu user load in percentage
   cpuinfo->fSys     = kernel_ratio; // cpu sys load in percentage
   cpuinfo->fTotal   = total_ratio; // cpu user+sys load in percentage
   cpuinfo->fIdle    = idle_ratio; // cpu idle percentage
}

////////////////////////////////////////////////////////////////////////////////
/// Get VM stat for Windows NT.

static void GetWinNTMemInfo(MemInfo_t *meminfo)
{
   Long64_t total, used, free, swap_total, swap_used, swap_avail;
   MEMORYSTATUSEX statex;
   statex.dwLength = sizeof(statex);
   if (!GlobalMemoryStatusEx(&statex)) {
      ::Error("GetWinNTMemInfo", "Error on GlobalMemoryStatusEx()");
      return;
   }
   used  = (Long64_t)(statex.ullTotalPhys - statex.ullAvailPhys);
   free  = (Long64_t) statex.ullAvailPhys;
   total = (Long64_t) statex.ullTotalPhys;

   meminfo->fMemTotal  = (Int_t) (total >> 20); // divide by 1024 * 1024
   meminfo->fMemUsed   = (Int_t) (used >> 20);
   meminfo->fMemFree   = (Int_t) (free >> 20);

   swap_total = (Long64_t)(statex.ullTotalPageFile - statex.ullTotalPhys);
   swap_avail = (Long64_t)(statex.ullAvailPageFile - statex.ullAvailPhys);
   swap_used  = swap_total - swap_avail;

   meminfo->fSwapTotal = (Int_t) (swap_total >> 20);
   meminfo->fSwapUsed  = (Int_t) (swap_used >> 20);
   meminfo->fSwapFree  = (Int_t) (swap_avail >> 20);
}

////////////////////////////////////////////////////////////////////////////////
/// Get process info for this process on Windows NT.

static void GetWinNTProcInfo(ProcInfo_t *procinfo)
{
   PROCESS_MEMORY_COUNTERS pmc;
   FILETIME    starttime, exittime, kerneltime, usertime;
   timeval     ru_stime, ru_utime;
   ULARGE_INTEGER li;

   typedef BOOL (__stdcall *GetProcessMemoryInfoProc)( HANDLE Process,
                 PPROCESS_MEMORY_COUNTERS ppsmemCounters, DWORD cb );
   static GetProcessMemoryInfoProc pGetProcessMemoryInfo = 0;

   HMODULE hModImagehlp = LoadLibrary( "Psapi.dll" );
   if (!hModImagehlp) {
      ::Error("GetWinNTProcInfo", "Error on LoadLibrary(Psapi.dll)");
      return;
   }

   pGetProcessMemoryInfo = (GetProcessMemoryInfoProc) GetProcAddress(
                            hModImagehlp, "GetProcessMemoryInfo" );
   if (!pGetProcessMemoryInfo) {
      ::Error("GetWinNTProcInfo",
              "Error on GetProcAddress(GetProcessMemoryInfo)");
      return;
   }

   if ( pGetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc)) ) {
      procinfo->fMemResident = pmc.WorkingSetSize / 1024;
      procinfo->fMemVirtual  = pmc.PagefileUsage / 1024;
   }
   if ( GetProcessTimes(GetCurrentProcess(), &starttime, &exittime,
      &kerneltime, &usertime)) {

      /* Convert FILETIMEs (0.1 us) to struct timeval */
      memcpy(&li, &kerneltime, sizeof(FILETIME));
      li.QuadPart /= 10L;         /* Convert to microseconds */
      ru_stime.tv_sec = li.QuadPart / 1000000L;
      ru_stime.tv_usec = li.QuadPart % 1000000L;

      memcpy(&li, &usertime, sizeof(FILETIME));
      li.QuadPart /= 10L;         /* Convert to microseconds */
      ru_utime.tv_sec = li.QuadPart / 1000000L;
      ru_utime.tv_usec = li.QuadPart % 1000000L;

      procinfo->fCpuUser = (Float_t)(ru_utime.tv_sec) +
                           ((Float_t)(ru_utime.tv_usec) / 1000000.);
      procinfo->fCpuSys  = (Float_t)(ru_stime.tv_sec) +
                           ((Float_t)(ru_stime.tv_usec) / 1000000.);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns static system info, like OS type, CPU type, number of CPUs
/// RAM size, etc into the SysInfo_t structure. Returns -1 in case of error,
/// 0 otherwise.

Int_t TWinNTSystem::GetSysInfo(SysInfo_t *info) const
{
   if (!info) return -1;
   GetWinNTSysInfo(info);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns cpu load average and load info into the CpuInfo_t structure.
/// Returns -1 in case of error, 0 otherwise. Use sampleTime to set the
/// interval over which the CPU load will be measured, in ms (default 1000).

Int_t TWinNTSystem::GetCpuInfo(CpuInfo_t *info, Int_t sampleTime) const
{
   if (!info) return -1;
   GetWinNTCpuInfo(info, sampleTime);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns ram and swap memory usage info into the MemInfo_t structure.
/// Returns -1 in case of error, 0 otherwise.

Int_t TWinNTSystem::GetMemInfo(MemInfo_t *info) const
{
   if (!info) return -1;
   GetWinNTMemInfo(info);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns cpu and memory used by this process into the ProcInfo_t structure.
/// Returns -1 in case of error, 0 otherwise.

Int_t TWinNTSystem::GetProcInfo(ProcInfo_t *info) const
{
   if (!info) return -1;
   GetWinNTProcInfo(info);
   return 0;
}
