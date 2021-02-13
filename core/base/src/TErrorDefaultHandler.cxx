/// \file TErrorDefaultHandler.cxx
/// \date 2020-06-14

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef WIN32
#include <windows.h>
#endif

#include <TEnv.h>
#include <TError.h>
#include <ThreadLocalStorage.h>
#include <TSystem.h>
#include <Varargs.h>

#include <cstdio>
#include <cstdlib>
#include <cctype> // for tolower
#include <cstring> // for strdup
#include <mutex>

// Integrate with macOS crash reporter.
#ifdef __APPLE__
extern "C" {
static const char *__crashreporter_info__ = 0;
asm(".desc ___crashreporter_info__, 0x10");
}
#endif


/// Serializes error output, destructed by the gROOT destructor via ReleaseDefaultErrorHandler()
static std::mutex *GetErrorMutex() {
   static std::mutex *m = new std::mutex();
   return m;
}


namespace ROOT {
namespace Internal {

void ReleaseDefaultErrorHandler()
{
   delete GetErrorMutex();
}

} // Internal namespace
} // ROOT namespace


/// Print debugging message to stderr and, on Windows, to the system debugger.
static void DebugPrint(const char *fmt, ...)
{
   TTHREAD_TLS(Int_t) buf_size = 2048;
   TTHREAD_TLS(char*) buf = 0;

   va_list ap;
   va_start(ap, fmt);

again:
   if (!buf)
      buf = new char[buf_size];

   Int_t n = vsnprintf(buf, buf_size, fmt, ap);
   // old vsnprintf's return -1 if string is truncated new ones return
   // total number of characters that would have been written
   if (n == -1 || n >= buf_size) {
      if (n == -1)
         buf_size *= 2;
      else
         buf_size = n+1;
      delete [] buf;
      buf = 0;
      va_end(ap);
      va_start(ap, fmt);
      goto again;
   }
   va_end(ap);

   // Serialize the actual printing.
   std::lock_guard<std::mutex> guard(*GetErrorMutex());

   const char *toprint = buf; // Work around for older platform where we use TThreadTLSWrapper
   fprintf(stderr, "%s", toprint);

#ifdef WIN32
   ::OutputDebugString(buf);
#endif
}


/// The default error handler function. It prints the message on stderr and
/// if abort is set it aborts the application.  Replaces the minimal error handler
/// of TError.h as part of the gROOT construction.  TError's minimal handler is put
/// back in place during the gROOT destruction.
void DefaultErrorHandler(Int_t level, Bool_t abort_bool, const char *location, const char *msg)
{
   if (gErrorIgnoreLevel == kUnset) {
      std::lock_guard<std::mutex> guard(*GetErrorMutex());

      gErrorIgnoreLevel = 0;
      if (gEnv) {
         std::string slevel;
         auto cstrlevel = gEnv->GetValue("Root.ErrorIgnoreLevel", "Print");
         while (cstrlevel && *cstrlevel) {
            slevel.push_back(tolower(*cstrlevel));
            cstrlevel++;
         }

         if (slevel == "print")
            gErrorIgnoreLevel = kPrint;
         else if (slevel == "info")
            gErrorIgnoreLevel = kInfo;
         else if (slevel == "warning")
            gErrorIgnoreLevel = kWarning;
         else if (slevel == "error")
            gErrorIgnoreLevel = kError;
         else if (slevel == "break")
            gErrorIgnoreLevel = kBreak;
         else if (slevel == "syserror")
            gErrorIgnoreLevel = kSysError;
         else if (slevel == "fatal")
            gErrorIgnoreLevel = kFatal;
      }
   }

   if (level < gErrorIgnoreLevel)
      return;

   const char *type = 0;

   if (level >= kInfo)
      type = "Info";
   if (level >= kWarning)
      type = "Warning";
   if (level >= kError)
      type = "Error";
   if (level >= kBreak)
      type = "\n *** Break ***";
   if (level >= kSysError)
      type = "SysError";
   if (level >= kFatal)
      type = "Fatal";

   std::string smsg;
   if (level >= kPrint && level < kInfo)
      smsg = msg;
   else if (level >= kBreak && level < kSysError)
      smsg = std::string(type) + " " + msg;
   else if (!location || !location[0])
      smsg = std::string(type) + ": " + msg;
   else
      smsg = std::string(type) + " in <" + location + ">: " + msg;

   DebugPrint("%s\n", smsg.c_str());

   fflush(stderr);
   if (abort_bool) {

#ifdef __APPLE__
      if (__crashreporter_info__)
         delete [] __crashreporter_info__;
      __crashreporter_info__ = strdup(smsg.c_str());
#endif

      DebugPrint("aborting\n");
      fflush(stderr);
      if (gSystem) {
         gSystem->StackTrace();
         gSystem->Abort();
      } else {
         abort();
      }
   }
}
