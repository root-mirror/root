// @(#)root/main:$Id$
// Author: Fons Rademakers   15/02/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// PMain                                                                //
//                                                                      //
// Main program used to create PROOF server application.                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <errno.h>

#ifdef WIN32
#include <io.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif
#include <stdlib.h>
#include <sys/types.h>

#include <ROOT/RConfig.hxx>
#include "RConfigure.h"
#include "TApplication.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TSystem.h"


static Int_t gLogLevel = 0;

// Special type for the hook to the TXProofServ constructor, needed to avoid
// using the plugin manager
typedef TApplication *(*TProofServ_t)(Int_t *argc, char **argv, FILE *flog);

////////////////////////////////////////////////////////////////////////////////
/// Read envs from file 'envfile' and add them to the env space

static void ReadPutEnvs(const char *envfile)
{
   // Check inputs
   if (!envfile || strlen(envfile) <= 0) return;

   // Open the file
   FILE *fenv = fopen(envfile, "r");
   if (!fenv) return;

   // Read lines
   char ln[4096];
   while (fgets(ln, sizeof(ln), fenv)) {
      int l = strlen(ln);
      // Strip '\n'
      if (l > 0 && ln[l-1] == '\n') { ln[l-1] = '\0'; l--; }
      // Skip comments or empty line
      if (l <= 0 || ln[0] == '#') continue;
      // Skip lines not in the form '<name>=<value>'
      if (!strchr(ln, '=')) continue;
      // Good line
      char *ev = new char[l+1];
      strcpy(ev, ln);
      putenv(ev);
   }

   // Close the file
   fclose(fenv);
}

////////////////////////////////////////////////////////////////////////////////
/// Redirect stdout to 'logfile'. This log file will be flushed to the
/// client or master after each command.
/// If donotredir != 0 just reopen the file for usage in TProofServ (already redirected).
/// On success return a pointer to the open log file. Return 0 on failure.

static FILE *RedirectOutput(const char *logfile, const char *loc, Int_t donotredir)
{
   if (loc)
      fprintf(stderr, "%s: RedirectOutput: enter: %s (do-not-redir: %d)\n",
                      loc, logfile, donotredir);

   if (donotredir == 0) {
      if (!logfile || strlen(logfile) <= 0) {
         fprintf(stderr,"%s: RedirectOutput: logfile path undefined\n", loc);
         return 0;
      }

      if (loc)
         fprintf(stderr,"%s: RedirectOutput: reopen %s\n", loc, logfile);
      FILE *flog = freopen(logfile, "a", stdout);
      if (!flog) {
         fprintf(stderr,"%s: RedirectOutput: could not freopen stdout\n", loc);
         return 0;
      }

      if (loc)
         fprintf(stderr,"%s: RedirectOutput: dup2 ...\n", loc);
      if ((dup2(fileno(stdout), fileno(stderr))) < 0) {
         fprintf(stderr,"%s: RedirectOutput: could not redirect stderr\n", loc);
         return 0;
      }
   }

   if (loc)
      fprintf(stderr,"%s: RedirectOutput: read open ...\n", loc);
   FILE *fLog = fopen(logfile, "r");
   if (!fLog) {
      fprintf(stderr,"%s: RedirectOutput: could not open logfile %s\n", loc, logfile);
      return 0;
   }

   if (loc)
      fprintf(stderr,"%s: RedirectOutput: done!\n", loc);
   // We are done
   return fLog;
}

////////////////////////////////////////////////////////////////////////////////
/// Set limits on the address space (virtual memory) if required.

static void SetMaxMemLimits(const char *prog)
{
#ifndef WIN32
   const char *assoft = gSystem->Getenv("ROOTPROOFASSOFT");
   const char *ashard = gSystem->Getenv("ROOTPROOFASHARD");

   if (assoft || ashard) {
      struct rlimit aslim, aslimref;
      if (getrlimit(RLIMIT_AS, &aslimref) != 0) {
         fprintf(stderr,"%s: problems getting RLIMIT_AS values (errno: %d)\n", prog, errno);
         exit(1);
      }
      if (gLogLevel > 0)
         fprintf(stderr, "%s: memory limits currently set to %lld (soft) and %lld (hard) bytes\n",
                         prog, (Long64_t)aslimref.rlim_cur, (Long64_t)aslimref.rlim_max);
      aslim.rlim_cur = aslimref.rlim_cur;
      aslim.rlim_max = aslimref.rlim_max;
      if (assoft) {
         Long_t rlim_cur = strtol(assoft, 0, 10);
         if (rlim_cur < kMaxLong && rlim_cur > 0)
            aslim.rlim_cur = (rlim_t) rlim_cur * (1024 * 1024);
      }
      if (ashard) {
         Long_t rlim_max = strtol(ashard, 0, 10);
         if (rlim_max < kMaxLong && rlim_max > 0)
            aslim.rlim_max = (rlim_t) rlim_max * (1024 * 1024);
      }
      // Change the limits, if required
      if ((aslim.rlim_cur != aslimref.rlim_cur) || (aslim.rlim_max != aslimref.rlim_max)) {
         fprintf(stderr, "%s: setting memory limits to %lld (soft) and %lld (hard) bytes\n",
                         prog, (Long64_t)aslim.rlim_cur, (Long64_t)aslim.rlim_max);
         if (setrlimit(RLIMIT_AS, &aslim) != 0) {
            fprintf(stderr,"%s: problems setting RLIMIT_AS values (errno: %d)\n", prog, errno);
            exit(1);
         }
      }
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// PROOF server main program.

int main(int argc, char **argv)
{
#ifdef R__DEBUG
   int debug = 1;
   while (debug)
      ;
#endif
   if (argc >= 6) {
      // Read and put system envs
      ReadPutEnvs(argv[5]);
   }

   gLogLevel = (argc >= 5) ? strtol(argv[4], 0, 10) : -1;
   if (gLogLevel < 0 && gSystem->Getenv("ROOTPROOFLOGLEVEL"))
      gLogLevel = atoi(gSystem->Getenv("ROOTPROOFLOGLEVEL"));
   if (gLogLevel > 0)
      fprintf(stderr,"%s: starting %s\n", argv[1], argv[0]);

   // Redirect the output
   FILE *fLog = 0;
   const char *loc = 0;
   const char *logfile = gSystem->Getenv("ROOTPROOFLOGFILE");
   Int_t donotredir = 0;
   if (gSystem->Getenv("ROOTPROOFDONOTREDIR")) {
      donotredir++;
      TString anr(gSystem->Getenv("ROOTPROOFDONOTREDIR"));
      if (anr.IsDigit()) donotredir = anr.Atoi();
   }
   if (logfile && donotredir != 1) {
      loc = (gLogLevel > 0) ? argv[1] : 0;
      if (gLogLevel > 0)
         fprintf(stderr,"%s: redirecting output to %s\n", argv[1], logfile);
      if (!(fLog = RedirectOutput(logfile, loc, donotredir))) {
         fprintf(stderr,"%s: problems redirecting output to file %s\n", argv[1], logfile);
         exit(1);
      }
   }
   if (gLogLevel > 0)
      fprintf(stderr,"%s: output redirected to: %s\n",
             argv[1], (logfile ? logfile : "+++not redirected+++"));

   SetMaxMemLimits(argv[1]);

   gROOT->SetBatch();
   TApplication *theApp = 0;

   TString getter("GetTXProofServ");
   TString prooflib = "libProofx";
   if (argc > 2) {
      if (!strcmp(argv[2], "lite")) {
         // Lite version for local processing: this is in libProof and has its own getter
         prooflib = "libProof";
         getter = "GetTProofServLite";
      } else if (strcmp(argv[2], "xpd")) {
         // Not XPD: set the appropriate names for library and getter
         prooflib = "libProof";
         getter = "GetTProofServ";
      }
   }
   char *p = 0;
   if ((p = gSystem->DynamicPathName(prooflib, kTRUE))) {
      delete[] p;
      if (gSystem->Load(prooflib) == -1) {
         fprintf(stderr,"%s: can't load %s\n", argv[1], prooflib.Data());
         exit(1);
      }
   } else {
      fprintf(stderr,"%s: can't locate %s\n", argv[1], prooflib.Data());
      exit(1);
   }

   // Locate constructor
   Func_t f = gSystem->DynFindSymbol(prooflib, getter);
   if (f) {
      theApp = (TApplication *) (*((TProofServ_t)f))(&argc, argv, fLog);
   } else {
      fprintf(stderr,"%s: can't find %s\n", argv[1], getter.Data());
      exit(1);
   }

   // Ready to run
   if (gLogLevel > 0)
      fprintf(stderr,"%s: running the TProofServ application\n", argv[1]);

   theApp->Run();

   // We can exit now
   gSystem->Exit(0);
}
