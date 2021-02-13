// @(#)root/base:$Id$
// Author: Fons Rademakers   22/12/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TApplication
\ingroup Base

This class creates the ROOT Application Environment that interfaces
to the windowing system eventloop and eventhandlers.
This class must be instantiated exactly once in any given
application. Normally the specific application class inherits from
TApplication (see TRint).
*/

#include "RConfigure.h"
#include "TApplication.h"
#include "TException.h"
#include "TGuiFactory.h"
#include "TVirtualX.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TError.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TTimer.h"
#include "TInterpreter.h"
#include "TStyle.h"
#include "TVirtualPad.h"
#include "TEnv.h"
#include "TColor.h"
#include "TPluginManager.h"
#include "TClassTable.h"
#include "TBrowser.h"
#include "TUrl.h"
#include "TVirtualMutex.h"
#include "TClassEdit.h"
#include "TMethod.h"
#include "TDataMember.h"
#include "TApplicationCommandLineOptionsHelp.h"
#include "TPRegexp.h"
#include <cstdlib>
#include <iostream>
#include <fstream>

TApplication *gApplication = nullptr;
Bool_t TApplication::fgGraphNeeded = kFALSE;
Bool_t TApplication::fgGraphInit = kFALSE;
TList *TApplication::fgApplications = nullptr;  // List of available applications

////////////////////////////////////////////////////////////////////////////////

class TIdleTimer : public TTimer {
public:
   TIdleTimer(Long_t ms) : TTimer(ms, kTRUE) { }
   Bool_t Notify();
};

////////////////////////////////////////////////////////////////////////////////
/// Notify handler.

Bool_t TIdleTimer::Notify()
{
   gApplication->HandleIdleTimer();
   Reset();
   return kFALSE;
}


ClassImp(TApplication);

static void CallEndOfProcessCleanups()
{
   // Insure that the files, canvases and sockets are closed.

   // If we get here, the tear down has started.  We have no way to know what
   // has or has not yet been done.  In particular on Ubuntu, this was called
   // after the function static in TSystem.cxx has been destructed.  So we
   // set gROOT in its end-of-life mode which prevents executing code, like
   // autoloading libraries (!) that is pointless ...
   if (gROOT) {
      gROOT->SetBit(kInvalidObject);
      gROOT->EndOfProcessCleanups();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Default ctor. Can be used by classes deriving from TApplication.

TApplication::TApplication() :
   fArgc(0), fArgv(0), fAppImp(0), fIsRunning(kFALSE), fReturnFromRun(kFALSE),
   fNoLog(kFALSE), fNoLogo(kFALSE), fQuit(kFALSE), fUseMemstat(kFALSE),
   fFiles(0), fIdleTimer(0), fSigHandler(0), fExitOnException(kDontExit),
   fAppRemote(0)
{
   ResetBit(kProcessRemotely);
}

////////////////////////////////////////////////////////////////////////////////
/// Create an application environment. The application environment
/// provides an interface to the graphics system and eventloop
/// (be it X, Windows, MacOS or BeOS). After creating the application
/// object start the eventloop by calling its Run() method. The command
/// line options recognized by TApplication are described in the GetOptions()
/// method. The recognized options are removed from the argument array.
/// The original list of argument options can be retrieved via the Argc()
/// and Argv() methods. The appClassName "proofserv" is reserved for the
/// PROOF system. The "options" and "numOptions" arguments are not used,
/// except if you want to by-pass the argv processing by GetOptions()
/// in which case you should specify numOptions<0. All options will
/// still be available via the Argv() method for later use.

TApplication::TApplication(const char *appClassName, Int_t *argc, char **argv,
                           void * /*options*/, Int_t numOptions) :
   fArgc(0), fArgv(0), fAppImp(0), fIsRunning(kFALSE), fReturnFromRun(kFALSE),
   fNoLog(kFALSE), fNoLogo(kFALSE), fQuit(kFALSE), fUseMemstat(kFALSE),
   fFiles(0), fIdleTimer(0), fSigHandler(0), fExitOnException(kDontExit),
   fAppRemote(0)
{
   R__LOCKGUARD(gInterpreterMutex);

   // Create the list of applications the first time
   if (!fgApplications)
      fgApplications = new TList;

   // Add the new TApplication early, so that the destructor of the
   // default TApplication (if it is called in the block of code below)
   // will not destroy the files, socket or TColor that have already been
   // created.
   fgApplications->Add(this);

   if (gApplication && gApplication->TestBit(kDefaultApplication)) {
      // allow default TApplication to be replaced by a "real" TApplication
      delete gApplication;
      gApplication = 0;
      gROOT->SetBatch(kFALSE);
      fgGraphInit = kFALSE;
   }

   if (gApplication) {
      Error("TApplication", "only one instance of TApplication allowed");
      fgApplications->Remove(this);
      return;
   }

   if (!gROOT)
      ::Fatal("TApplication::TApplication", "ROOT system not initialized");

   if (!gSystem)
      ::Fatal("TApplication::TApplication", "gSystem not initialized");

   static Bool_t hasRegisterAtExit(kFALSE);
   if (!hasRegisterAtExit) {
      // If we are the first TApplication register the atexit)
      atexit(CallEndOfProcessCleanups);
      hasRegisterAtExit = kTRUE;
   }
   gROOT->SetName(appClassName);

   // copy command line arguments, can be later accessed via Argc() and Argv()
   if (argc && *argc > 0) {
      fArgc = *argc;
      fArgv = (char **)new char*[fArgc];
   }

   for (int i = 0; i < fArgc; i++)
      fArgv[i] = StrDup(argv[i]);

   if (numOptions >= 0)
      GetOptions(argc, argv);

   if (fArgv)
      gSystem->SetProgname(fArgv[0]);

   // Tell TSystem the TApplication has been created
   gSystem->NotifyApplicationCreated();

   fAppImp = gGuiFactory->CreateApplicationImp(appClassName, argc, argv);
   ResetBit(kProcessRemotely);

   // Initialize the graphics environment
   if (gClassTable->GetDict("TPad")) {
      fgGraphNeeded = kTRUE;
      InitializeGraphics();
   }

   // Save current interpreter context
   gInterpreter->SaveContext();
   gInterpreter->SaveGlobalsContext();

   // to allow user to interact with TCanvas's under WIN32
   gROOT->SetLineHasBeenProcessed();

   // activate TMemStat
   if (fUseMemstat || gEnv->GetValue("Root.TMemStat", 0)) {
      fUseMemstat = kTRUE;
      Int_t buffersize = gEnv->GetValue("Root.TMemStat.buffersize", 100000);
      Int_t maxcalls   = gEnv->GetValue("Root.TMemStat.maxcalls", 5000000);
      const char *ssystem = gEnv->GetValue("Root.TMemStat.system","gnubuiltin");
      if (maxcalls > 0) {
         gROOT->ProcessLine(Form("new TMemStat(\"%s\",%d,%d);",ssystem,buffersize,maxcalls));
      }
   }

   //Needs to be done last
   gApplication = this;
   gROOT->SetApplication(this);

}

////////////////////////////////////////////////////////////////////////////////
/// TApplication dtor.

TApplication::~TApplication()
{
   for (int i = 0; i < fArgc; i++)
      if (fArgv[i]) delete [] fArgv[i];
   delete [] fArgv;

   if (fgApplications)
      fgApplications->Remove(this);

   //close TMemStat
   if (fUseMemstat) {
      ProcessLine("TMemStat::Close()");
      fUseMemstat = kFALSE;
   }

   // Reduce the risk of the files or sockets being closed after the
   // end of 'main' (or more exactly before the library start being
   // unloaded).
   if (fgApplications == 0 || fgApplications->FirstLink() == 0 ) {
      TROOT::ShutDown();
   }

   // Now that all the canvases and files have been closed we can
   // delete the implementation.
   SafeDelete(fAppImp);
}

////////////////////////////////////////////////////////////////////////////////
/// Static method. This method should be called from static library
/// initializers if the library needs the low level graphics system.

void TApplication::NeedGraphicsLibs()
{
   fgGraphNeeded = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the graphics environment.

void TApplication::InitializeGraphics()
{
   if (fgGraphInit || !fgGraphNeeded) return;

   // Load the graphics related libraries
   LoadGraphicsLibs();

   // Try to load TrueType font renderer. Only try to load if not in batch
   // mode and Root.UseTTFonts is true and Root.TTFontPath exists. Abort silently
   // if libttf or libGX11TTF are not found in $ROOTSYS/lib or $ROOTSYS/ttf/lib.
   const char *ttpath = gEnv->GetValue("Root.TTFontPath",
                                       TROOT::GetTTFFontDir());
   char *ttfont = gSystem->Which(ttpath, "arialbd.ttf", kReadPermission);
   // Check for use of DFSG - fonts
   if (!ttfont)
      ttfont = gSystem->Which(ttpath, "FreeSansBold.ttf", kReadPermission);

#if !defined(R__WIN32)
   if (!gROOT->IsBatch() && !strcmp(gVirtualX->GetName(), "X11") &&
       ttfont && gEnv->GetValue("Root.UseTTFonts", 1)) {
      if (gClassTable->GetDict("TGX11TTF")) {
         // in principle we should not have linked anything against libGX11TTF
         // but with ACLiC this can happen, initialize TGX11TTF by hand
         // (normally this is done by the static library initializer)
         ProcessLine("TGX11TTF::Activate();");
      } else {
         TPluginHandler *h;
         if ((h = gROOT->GetPluginManager()->FindHandler("TVirtualX", "x11ttf")))
            if (h->LoadPlugin() == -1)
               Info("InitializeGraphics", "no TTF support");
      }
   }
#endif
   delete [] ttfont;

   // Create WM dependent application environment
   if (fAppImp)
      delete fAppImp;
   fAppImp = gGuiFactory->CreateApplicationImp(gROOT->GetName(), &fArgc, fArgv);
   if (!fAppImp) {
      MakeBatch();
      fAppImp = gGuiFactory->CreateApplicationImp(gROOT->GetName(), &fArgc, fArgv);
   }

   // Create the canvas colors early so they are allocated before
   // any color table expensive bitmaps get allocated in GUI routines (like
   // creation of XPM bitmaps).
   TColor::InitializeColors();

   // Hook for further initializing the WM dependent application environment
   Init();

   // Set default screen factor (if not disabled in rc file)
   if (gEnv->GetValue("Canvas.UseScreenFactor", 1)) {
      Int_t  x, y;
      UInt_t w, h;
      if (gVirtualX) {
         gVirtualX->GetGeometry(-1, x, y, w, h);
         if (h > 0 && h < 1000) gStyle->SetScreenFactor(0.0011*h);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clear list containing macro files passed as program arguments.
/// This method is called from TRint::Run() to ensure that the macro
/// files are only executed the first time Run() is called.

void TApplication::ClearInputFiles()
{
   if (fFiles) {
      fFiles->Delete();
      SafeDelete(fFiles);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return specified argument.

char *TApplication::Argv(Int_t index) const
{
   if (fArgv) {
      if (index >= fArgc) {
         Error("Argv", "index (%d) >= number of arguments (%d)", index, fArgc);
         return 0;
      }
      return fArgv[index];
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get and handle command line options. Arguments handled are removed
/// from the argument array. See CommandLineOptionsHelp.h for options.

void TApplication::GetOptions(Int_t *argc, char **argv)
{
   static char null[1] = { "" };

   fNoLog = kFALSE;
   fQuit  = kFALSE;
   fFiles = 0;

   if (!argc)
      return;

   int i, j;
   TString pwd;

   for (i = 1; i < *argc; i++) {
      if (!strcmp(argv[i], "-?") || !strncmp(argv[i], "-h", 2) ||
          !strncmp(argv[i], "--help", 6)) {
         fprintf(stderr, kCommandLineOptionsHelp);
         Terminate(0);
      } else if (!strncmp(argv[i], "--version", 9)) {
         fprintf(stderr, "ROOT Version: %s\n", gROOT->GetVersion());
         fprintf(stderr, "Built for %s on %s\n",
                 gSystem->GetBuildArch(),
                 gROOT->GetGitDate());

         fprintf(stderr, "From %s@%s\n",
                gROOT->GetGitBranch(),
                gROOT->GetGitCommit());

         Terminate(0);
      } else if (!strcmp(argv[i], "-config")) {
         fprintf(stderr, "ROOT ./configure options:\n%s\n", gROOT->GetConfigOptions());
         Terminate(0);
      } else if (!strcmp(argv[i], "-memstat")) {
         fUseMemstat = kTRUE;
         argv[i] = null;
      } else if (!strcmp(argv[i], "-b")) {
         MakeBatch();
         argv[i] = null;
      } else if (!strcmp(argv[i], "-n")) {
         fNoLog = kTRUE;
         argv[i] = null;
      } else if (!strcmp(argv[i], "-t")) {
         ROOT::EnableImplicitMT();
         // EnableImplicitMT() only enables thread safety if IMT was configured;
         // enable thread safety even with IMT off:
         ROOT::EnableThreadSafety();
         argv[i] = null;
      } else if (!strcmp(argv[i], "-q")) {
         fQuit = kTRUE;
         argv[i] = null;
      } else if (!strcmp(argv[i], "-l")) {
         // used by front-end program to not display splash screen
         fNoLogo = kTRUE;
         argv[i] = null;
      } else if (!strcmp(argv[i], "-x")) {
         fExitOnException = kExit;
         argv[i] = null;
      } else if (!strcmp(argv[i], "-splash")) {
         // used when started by front-end program to signal that
         // splash screen can be popped down (TRint::PrintLogo())
         argv[i] = null;
      } else if (strncmp(argv[i], "--web", 5) == 0) {
         // the web mode is requested
         const char *opt = argv[i] + 5;
         argv[i] = null;
         TString argw;
         if (gROOT->IsBatch()) argw = "batch";
         if (*opt == '=') argw.Append(opt+1);
         if (gSystem->Load("libROOTWebDisplay") >= 0) {
            gROOT->SetWebDisplay(argw.Data());
            gEnv->SetValue("Gui.Factory", "web");
         } else {
            Error("GetOptions", "--web option not supported, ROOT should be built with at least c++14 enabled");
         }
      } else if (!strcmp(argv[i], "-e")) {
         argv[i] = null;
         ++i;

         if ( i < *argc ) {
            if (!fFiles) fFiles = new TObjArray;
            TObjString *expr = new TObjString(argv[i]);
            expr->SetBit(kExpression);
            fFiles->Add(expr);
            argv[i] = null;
         } else {
            Warning("GetOptions", "-e must be followed by an expression.");
         }
      } else if (!strcmp(argv[i], "--")) {
         TObjString* macro = nullptr;
         bool warnShown = false;

         if (fFiles) {
            for (auto f: *fFiles) {
               TObjString* file = dynamic_cast<TObjString*>(f);
               if (!file) {
                  if (!dynamic_cast<TNamed*>(f)) {
                     Error("GetOptions()", "Inconsistent file entry (not a TObjString)!");
                     f->Dump();
                  } // else we did not find the file.
                  continue;
               }

               if (file->TestBit(kExpression))
                  continue;
               if (file->String().EndsWith(".root"))
                  continue;
               if (file->String().Contains('('))
                  continue;

               if (macro && !warnShown && (warnShown = true))
                  Warning("GetOptions", "-- is used with several macros. "
                                        "The arguments will be passed to the last one.");

               macro = file;
            }
         }

         if (macro) {
            argv[i] = null;
            ++i;
            TString& str = macro->String();

            str += '(';
            for (; i < *argc; i++) {
               str += argv[i];
               str += ',';
               argv[i] = null;
            }
            str.EndsWith(",") ? str[str.Length() - 1] = ')' : str += ')';
         } else {
            Warning("GetOptions", "no macro to pass arguments to was provided. "
                                  "Everything after the -- will be ignored.");
            for (; i < *argc; i++)
               argv[i] = null;
         }
      } else if (argv[i][0] != '-' && argv[i][0] != '+') {
         Long64_t size;
         Long_t id, flags, modtime;
         char *arg = strchr(argv[i], '(');
         if (arg) *arg = '\0';
         char *dir = gSystem->ExpandPathName(argv[i]);
         // ROOT-9959: we do not continue if we could not expand the path
         if (!dir) continue;
         TUrl udir(dir, kTRUE);
         // remove options and anchor to check the path
         TString sfx = udir.GetFileAndOptions();
         TString fln = udir.GetFile();
         sfx.Replace(sfx.Index(fln), fln.Length(), "");
         TString path = udir.GetFile();
         if (strcmp(udir.GetProtocol(), "file")) {
            path = udir.GetUrl();
            path.Replace(path.Index(sfx), sfx.Length(), "");
         }
         // 'path' is the full URL without suffices (options and/or anchor)
         if (arg) *arg = '(';
         if (!arg && !gSystem->GetPathInfo(path.Data(), &id, &size, &flags, &modtime)) {
            if ((flags & 2)) {
               // if directory set it in fWorkDir
               if (pwd == "") {
                  pwd = gSystem->WorkingDirectory();
                  fWorkDir = dir;
                  gSystem->ChangeDirectory(dir);
                  argv[i] = null;
               } else if (!strcmp(gROOT->GetName(), "Rint")) {
                  Warning("GetOptions", "only one directory argument can be specified (%s)", dir);
               }
            } else if (size > 0) {
               // if file add to list of files to be processed
               if (!fFiles) fFiles = new TObjArray;
               fFiles->Add(new TObjString(path.Data()));
               argv[i] = null;
            } else {
               Warning("GetOptions", "file %s has size 0, skipping", dir);
            }
         } else {
            if (TString(udir.GetFile()).EndsWith(".root")) {
               if (!strcmp(udir.GetProtocol(), "file")) {
                  // file ending on .root but does not exist, likely a typo
                  // warn user if plain root...
                  if (!strcmp(gROOT->GetName(), "Rint"))
                     Warning("GetOptions", "file %s not found", dir);
               } else {
                  // remote file, give it the benefit of the doubt and add it to list of files
                  if (!fFiles) fFiles = new TObjArray;
                  fFiles->Add(new TObjString(argv[i]));
                  argv[i] = null;
               }
            } else {
               TString mode,fargs,io;
               TString fname = gSystem->SplitAclicMode(dir,mode,fargs,io);
               char *mac;
               if (!fFiles) fFiles = new TObjArray;
               if ((mac = gSystem->Which(TROOT::GetMacroPath(), fname,
                                         kReadPermission))) {
                  // if file add to list of files to be processed
                  fFiles->Add(new TObjString(argv[i]));
                  argv[i] = null;
                  delete [] mac;
               } else {
                  // if file add an invalid entry to list of files to be processed
                  fFiles->Add(new TNamed("NOT FOUND!", argv[i]));
                  // only warn if we're plain root,
                  // other progs might have their own params
                  if (!strcmp(gROOT->GetName(), "Rint"))
                     Warning("GetOptions", "macro %s not found", fname.Data());
               }
            }
         }
         delete [] dir;
      }
      // ignore unknown options
   }

   // go back to startup directory
   if (pwd != "")
      gSystem->ChangeDirectory(pwd);

   // remove handled arguments from argument array
   j = 0;
   for (i = 0; i < *argc; i++) {
      if (strcmp(argv[i], "")) {
         argv[j] = argv[i];
         j++;
      }
   }

   *argc = j;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle idle timeout. When this timer expires the registered idle command
/// will be executed by this routine and a signal will be emitted.

void TApplication::HandleIdleTimer()
{
   if (!fIdleCommand.IsNull())
      ProcessLine(GetIdleCommand());

   Emit("HandleIdleTimer()");
}

////////////////////////////////////////////////////////////////////////////////
/// Handle exceptions (kSigBus, kSigSegmentationViolation,
/// kSigIllegalInstruction and kSigFloatingException) trapped in TSystem.
/// Specific TApplication implementations may want something different here.

void TApplication::HandleException(Int_t sig)
{
   if (TROOT::Initialized()) {
      if (gException) {
         gInterpreter->RewindDictionary();
         gInterpreter->ClearFileBusy();
      }
      if (fExitOnException == kExit)
         gSystem->Exit(128 + sig);
      else if (fExitOnException == kAbort)
         gSystem->Abort();
      else
         Throw(sig);
   }
   gSystem->Exit(128 + sig);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the exit on exception option. Setting this option determines what
/// happens in HandleException() in case an exception (kSigBus,
/// kSigSegmentationViolation, kSigIllegalInstruction or kSigFloatingException)
/// is trapped. Choices are: kDontExit (default), kExit or kAbort.
/// Returns the previous value.

TApplication::EExitOnException TApplication::ExitOnException(TApplication::EExitOnException opt)
{
   EExitOnException old = fExitOnException;
   fExitOnException = opt;
   return old;
}

/////////////////////////////////////////////////////////////////////////////////
/// The function generates and executes a command that loads the Doxygen URL in
/// a browser. It works for Mac, Windows and Linux. In the case of Linux, the
/// function also checks if the DISPLAY is set. If it isn't, a warning message
/// and the URL will be displayed on the terminal.
///
/// \param[in] url web page to be displayed in a browser

void TApplication::OpenInBrowser(const TString &url)
{
   // We check what operating system the user has.
#ifdef R__MACOSX
   // Command for opening a browser on Mac.
   TString cMac("open ");
   // We generate the full command and execute it.
   cMac.Append(url);
   gSystem->Exec(cMac);
#elif defined(R__WIN32)
   // Command for opening a browser on Windows.
   TString cWindows("start ");
   cWindows.Append(url);
   gSystem->Exec(cWindows);
#else
   // Command for opening a browser in Linux.
   TString cLinux("xdg-open ");
   // For Linux we check if the DISPLAY is set.
   if (gSystem->Getenv("DISPLAY")) {
      // If the DISPLAY is set it will open the browser.
      cLinux.Append(url);
      gSystem->Exec(cLinux);
   } else {
      // Else the user will have a warning and the URL in the terminal.
      Warning("OpenInBrowser", "The $DISPLAY is not set! Please open (e.g. Ctrl-click) %s\n", url.Data());
   }
#endif
}

namespace {
enum EUrl { kURLforClass, kURLforNameSpace, kURLforStruct };
////////////////////////////////////////////////////////////////////////////////
/// The function generates a URL address for class or namespace (scopeName).
/// This is the URL to the online reference guide, generated by Doxygen.
/// With the enumeration "EUrl" we pick which case we need - the one for
/// class (kURLforClass) or the one for namespace (kURLforNameSpace).
///
/// \param[in] scopeName the name of the class or the namespace
/// \param[in] scopeType the enumerator for class or namespace

static TString UrlGenerator(TString scopeName, EUrl scopeType)
{
   // We start the URL with a static part, the same for all scopes and members.
   TString url = "https://root.cern/doc/";
   // Then we check the ROOT version used.
   TPRegexp re4(R"(.*/v(\d)-(\d\d)-00-patches)");
   const char *branchName = gROOT->GetGitBranch();
   TObjArray *objarr = re4.MatchS(branchName);
   TString version;
   // We extract the correct version name for the URL.
   if (objarr && objarr->GetEntries() == 3) {
      // We have a valid version of ROOT and we will extract the correct name for the URL.
      version = ((TObjString *)objarr->At(1))->GetString() + ((TObjString *)objarr->At(2))->GetString();
   } else {
      // If it's not a supported version, we will go to "master" branch.
      version = "master";
   }
   delete objarr;
   url.Append(version);
   url.Append("/");
   // We will replace all "::" with "_1_1" and all "_" with "__" in the
   // classes definitions, due to Doxygen syntax requirements.
   scopeName.ReplaceAll("_", "__");
   scopeName.ReplaceAll("::", "_1_1");
   // We build the URL for the correct scope type and name.
   if (scopeType == kURLforClass) {
      url.Append("class");
   } else if (scopeType == kURLforStruct) {
      url.Append("struct");
   } else {
      url.Append("namespace");
   }
   url.Append(scopeName);
   url.Append(".html");
   return url;
}
} // namespace

namespace {
////////////////////////////////////////////////////////////////////////////////
/// The function returns a TString with the arguments of a method from the
/// scope (scopeName), but modified with respect to Doxygen syntax - spacing
/// around special symbols and adding the missing scopes ("std::").
/// "FormatMethodArgsForDoxygen" works for functions defined inside namespaces
/// as well. We avoid looking up twice for the TFunction by passing "func".
///
/// \param[in] scopeName the name of the class/namespace/struct
/// \param[in] func pointer to the method

static TString FormatMethodArgsForDoxygen(const TString &scopeName, TFunction *func)
{
   // With "GetSignature" we get the arguments of the method and put them in a TString.
   TString methodArguments = func->GetSignature();
   // "methodArguments" is modified with respect of Doxygen requirements.
   methodArguments.ReplaceAll(" = ", "=");
   methodArguments.ReplaceAll("* ", " *");
   methodArguments.ReplaceAll("*=", " *=");
   methodArguments.ReplaceAll("*)", " *)");
   methodArguments.ReplaceAll("*,", " *,");
   methodArguments.ReplaceAll("*& ", " *&");
   methodArguments.ReplaceAll("& ", " &");
   // TODO: prepend "std::" to all stdlib classes!
   methodArguments.ReplaceAll("ostream", "std::ostream");
   methodArguments.ReplaceAll("istream", "std::istream");
   methodArguments.ReplaceAll("map", "std::map");
   methodArguments.ReplaceAll("vector", "std::vector");
   // We need to replace the "currentClass::foo" with "foo" in the arguments.
   // TODO: protect the global functions.
   TString scopeNameRE("\\b");
   scopeNameRE.Append(scopeName);
   scopeNameRE.Append("::\\b");
   TPRegexp argFix(scopeNameRE);
   argFix.Substitute(methodArguments, "");
   return methodArguments;
}
} // namespace

namespace {
////////////////////////////////////////////////////////////////////////////////
/// The function checks if a member function of a scope is defined as inline.
/// If so, it also checks if it is virtual. Then the return type of "func" is
/// modified for the need of Doxygen and with respect to the function
/// definition. We pass pointer to the method (func) to not re-do the
/// TFunction lookup.
///
/// \param[in] scopeName the name of the class/namespace/struct
/// \param[in] func pointer to the method

static TString FormatReturnTypeForDoxygen(const TString &scopeName, TFunction *func)
{
   // We put the return type of "func" in a TString "returnType".
   TString returnType = func->GetReturnTypeName();
   // If the return type is a type nested in the current class, it will appear scoped (Class::Enumeration).
   // Below we make sure to remove the current class, because the syntax of Doxygen requires it.
   TString scopeNameRE("\\b");
   scopeNameRE.Append(scopeName);
   scopeNameRE.Append("::\\b");
   TPRegexp returnFix(scopeNameRE);
   returnFix.Substitute(returnType, "");
   // We check is if the method is defined as inline.
   if (func->ExtraProperty() & kIsInlined) {
      // We check if the function is defined as virtual.
      if (func->Property() & kIsVirtual) {
         // If the function is virtual, we append "virtual" before the return type.
         returnType.Prepend("virtual ");
      }
      returnType.ReplaceAll(" *", "*");
   } else {
      // If the function is not inline we only change the spacing in "returnType"
      returnType.ReplaceAll("*", " *");
   }
   // In any case (with no respect to virtual/inline check) we need to change
   // the return type as following.
   // TODO: prepend "std::" to all stdlib classes!
   returnType.ReplaceAll("istream", "std::istream");
   returnType.ReplaceAll("ostream", "std::ostream");
   returnType.ReplaceAll("map", "std::map");
   returnType.ReplaceAll("vector", "std::vector");
   returnType.ReplaceAll("&", " &");
   return returnType;
}
} // namespace

namespace {
////////////////////////////////////////////////////////////////////////////////
/// The function generates a URL for "dataMemberName" defined in "scopeName".
/// It returns a TString with the URL used in the online reference guide,
/// generated with Doxygen. For data members the URL consist of 2 parts -
/// URL for "scopeName" and a part for "dataMemberName".
/// For enumerator, the URL could be separated into 3 parts - URL for
/// "scopeName", part for the enumeration and a part for the enumerator.
///
/// \param[in] scopeName the name of the class/namespace/struct
/// \param[in] dataMemberName the name of the data member/enumerator
/// \param[in] dataMember pointer to the data member/enumerator
/// \param[in] scopeType enumerator to the scope type

static TString
GetUrlForDataMember(const TString &scopeName, const TString &dataMemberName, TDataMember *dataMember, EUrl scopeType)
{
   // We first check if the data member is not enumerator.
   if (!dataMember->IsEnum()) {
      // If we work with data members, we have to append a hashed with MD5 text, consisting of:
      // "Type ClassName::DataMemberNameDataMemberName(arguments)".
      // We first get the type of the data member.
      TString md5DataMember(dataMember->GetFullTypeName());
      md5DataMember.Append(" ");
      // We append the scopeName and "::".
      md5DataMember.Append(scopeName);
      md5DataMember.Append("::");
      // We append the dataMemberName twice.
      md5DataMember.Append(dataMemberName);
      md5DataMember.Append(dataMemberName);
      // We call UrlGenerator for the scopeName.
      TString urlForDataMember = UrlGenerator(scopeName, scopeType);
      // Then we append "#a" and the hashed text.
      urlForDataMember.Append("#a");
      urlForDataMember.Append(md5DataMember.MD5());
      return urlForDataMember;
   }
   // If the data member is enumerator, then we first have to check if the enumeration is anonymous.
   // Doxygen requires different syntax for anonymous enumeration ("scopeName::@1@1").
   // We create a TString with the name of the scope and the enumeration from which the enumerator is.
   TString scopeEnumeration = dataMember->GetTrueTypeName();
   TString md5EnumClass;
   if (scopeEnumeration.Contains("(anonymous)")) {
      // FIXME: need to investigate the numbering scheme.
      md5EnumClass.Append(scopeName);
      md5EnumClass.Append("::@1@1");
   } else {
      // If the enumeration is not anonymous we put "scopeName::Enumeration" in a TString,
      // which will be hashed with MD5 later.
      md5EnumClass.Append(scopeEnumeration);
      // We extract the part after "::" (this is the enumerator name).
      TString enumOnlyName = TClassEdit::GetUnqualifiedName(scopeEnumeration);
      // The syntax is "Class::EnumeratorEnumerator
      md5EnumClass.Append(enumOnlyName);
   }
   // The next part of the URL is hashed "@ scopeName::EnumeratorEnumerator".
   TString md5Enumerator("@ ");
   md5Enumerator.Append(scopeName);
   md5Enumerator.Append("::");
   md5Enumerator.Append(dataMemberName);
   md5Enumerator.Append(dataMemberName);
   // We make the URL for the "scopeName".
   TString url = UrlGenerator(scopeName, scopeType);
   // Then we have to append the hashed text for the enumerator.
   url.Append("#a");
   url.Append(md5EnumClass.MD5());
   // We append "a" and then the next hashed text.
   url.Append("a");
   url.Append(md5Enumerator.MD5());
   return url;
}
} // namespace

namespace {
////////////////////////////////////////////////////////////////////////////////
/// The function generates URL for enumeration. The hashed text consist of:
/// "Class::EnumerationEnumeration".
///
/// \param[in] scopeName the name of the class/namespace/struct
/// \param[in] enumeration the name of the enumeration
/// \param[in] scopeType enumerator for class/namespace/struct

static TString GetUrlForEnumeration(TString scopeName, const TString &enumeration, EUrl scopeType)
{
   // The URL consists of URL for the "scopeName", "#a" and hashed as MD5 text.
   // The text is "Class::EnumerationEnumeration.
   TString md5Enumeration(scopeName);
   md5Enumeration.Append("::");
   md5Enumeration.Append(enumeration);
   md5Enumeration.Append(enumeration);
   // We make the URL for the scope "scopeName".
   TString url(UrlGenerator(scopeName, scopeType));
   // Then we have to append "#a" and the hashed text.
   url.Append("#a");
   url.Append(md5Enumeration.MD5());
   return url;
}
} // namespace

namespace {
enum EMethodKind { kURLforMethod, kURLforStructor };
////////////////////////////////////////////////////////////////////////////////
/// The function generates URL for any member function (including Constructor/
/// Destructor) of "scopeName". Doxygen first generates the URL for the scope.
/// We do that with the help of "UrlGenerator". Then we append "#a" and a
/// hashed with MD5 text. It consists of:
/// "ReturnType ScopeName::MethodNameMethodName(Method arguments)".
/// For constructor/destructor of a class, the return type is not appended.
///
/// \param[in] scopeName the name of the class/namespace/struct
/// \param[in] methodName the name of the method from the scope
/// \param[in] func pointer to the method
/// \param[in] methodType enumerator for method or constructor
/// \param[in] scopeType enumerator for class/namespace/struct

static TString GetUrlForMethod(const TString &scopeName, const TString &methodName, TFunction *func,
                               EMethodKind methodType, EUrl scopeType)
{
   TString md5Text;
   if (methodType == kURLforMethod) {
      // In the case of method, we append the return type too.
      // "FormatReturnTypeForDoxygen" modifies the return type with respect to Doxygen's requirement.
      md5Text.Append((FormatReturnTypeForDoxygen(scopeName, func)));
      if (scopeType == kURLforNameSpace) {
         // We need to append "constexpr" if we work with constexpr functions in namespaces.
         if (func->Property() & kIsConstexpr) {
            md5Text.Prepend("constexpr ");
         }
      }
      md5Text.Append(" ");
   }
   // We append ScopeName::MethodNameMethodName.
   md5Text.Append(scopeName);
   md5Text.Append("::");
   md5Text.Append(methodName);
   md5Text.Append(methodName);
   // We use "FormatMethodArgsForDoxygen" to modify the arguments of Method with respect of Doxygen.
   md5Text.Append(FormatMethodArgsForDoxygen(scopeName, func));
   // We generate the URL for the class/namespace/struct.
   TString url = UrlGenerator(scopeName, scopeType);
   url.Append("#a");
   // We append the hashed text.
   url.Append(md5Text.MD5());
   return url;
}
} // namespace


////////////////////////////////////////////////////////////////////////////////
/// It opens the online reference guide, generated with Doxygen, for the
/// chosen scope (class/namespace/struct) or member (method/function/
/// data member/enumeration/enumerator. If the user types incorrect value,
/// it will return an error or warning.
///
/// \param[in] strippedClass the scope or scope::member

void TApplication::OpenReferenceGuideFor(const TString &strippedClass)
{
   // We check if the user is searching for a scope and if the scope exists.
   if (TClass *clas = TClass::GetClass(strippedClass)) {
      // We check what scope he is searching for (class/namespace/struct).
      // Enumerators will switch between the possible cases.
      EUrl scopeType;
      if (clas->Property() & kIsNamespace) {
         scopeType = kURLforNameSpace;
      } else if (clas->Property() & kIsStruct) {
         scopeType = kURLforStruct;
      } else {
         scopeType = kURLforClass;
      }
      // If the user search directly for a scope we open the URL for him with OpenInBrowser.
      OpenInBrowser(UrlGenerator(strippedClass, scopeType));
      return;
   }
   // Else we subtract the name of the method and remove it from the command.
   TString memberName = TClassEdit::GetUnqualifiedName(strippedClass);
   // Error out if "strippedClass" is un-scoped (and it's not a class, see `TClass::GetClass(strippedClass)` above).
   // TODO: Global functions.
   if (strippedClass == memberName) {
      Error("OpenReferenceGuideFor", "Unknown entity \"%s\" - global variables / functions not supported yet!",
            strippedClass.Data());
      return;
   }
   // Else we remove the member name to be left with the scope.
   TString scopeName = strippedClass(0, strippedClass.Length() - memberName.Length() - 2);
   // We check if the scope exists in ROOT.
   TClass *cl = TClass::GetClass(scopeName);
   if (!cl) {
      // That's a member of something ROOT doesn't know.
      Warning("OpenReferenceGuideFor", "\"%s\" does not exist in ROOT!", scopeName.Data());
      return;
   }
   // We have enumerators for the three available cases - class, namespace and struct.
   EUrl scopeType;
   if (cl->Property() & kIsNamespace) {
      scopeType = kURLforNameSpace;
   } else if (cl->Property() & kIsStruct) {
      scopeType = kURLforStruct;
   } else {
      scopeType = kURLforClass;
   }
   // If the user wants to search for a method, we take its name (memberName) and
   // modify it - we delete everything starting at the first "(" so the user won't have to
   // do it by hand when they use Tab.
   int bracket = memberName.First("(");
   if (bracket > 0) {
      memberName.Remove(bracket);
   }
   // We check if "memberName" is a member function of "cl" or any of its base classes.
   if (TFunction *func = cl->GetMethodAllAny(memberName)) {
      // If so we find the name of the class that it belongs to.
      TString baseClName = ((TMethod *)func)->GetClass()->GetName();
      // We define an enumerator to distinguish between structor and method.
      EMethodKind methodType;
      // We check if "memberName" is a constructor.
      if (baseClName == memberName) {
         methodType = kURLforStructor;
         // We check if "memberName" is a destructor.
      } else if (memberName[0] == '~') {
         methodType = kURLforStructor;
         // We check if "memberName" is a method.
      } else {
         methodType = kURLforMethod;
      }
      // We call "GetUrlForMethod" for the correct class and scope.
      OpenInBrowser(GetUrlForMethod(baseClName, memberName, func, methodType, scopeType));
      return;
   }
   // We check if "memberName" is an enumeration.
   if (cl->GetListOfEnums()->FindObject(memberName)) {
      // If so with OpenInBrowser we open the URL generated with GetUrlForEnumeration
      // with respect to the "scopeType".
      OpenInBrowser(GetUrlForEnumeration(scopeName, memberName, scopeType));
      return;
   }

   // We check if "memberName" is enumerator defined in one the base classes of "scopeName".
   if (auto enumerator = (TDataMember *)cl->GetListOfAllPublicDataMembers()->FindObject(memberName)) {
      // We find the actual scope (might be in a base) and open the URL in a browser.
      TString baseClName = ((TMethod *)enumerator->GetClass())->GetName();
      OpenInBrowser(GetUrlForDataMember(baseClName, memberName, enumerator, scopeType));
      return;
   }

   // Warning message will appear if the user types the function name incorrectly
   // or the function is not a member function of "cl" or any of its base classes.
   Warning("Help", "cannot find \"%s\" as member of %s or its base classes! Check %s\n", memberName.Data(),
           scopeName.Data(), UrlGenerator(scopeName, scopeType).Data());
}

////////////////////////////////////////////////////////////////////////////////
/// The function lists useful commands (".help") or opens the online reference
/// guide, generated with Doxygen (".help scope" or ".help scope::member").
///
/// \param[in] line command from the command line

void TApplication::Help(const char *line)
{
   // We first check if the user wants to print the help on the interpreter.
   TString strippedCommand = TString(line).Strip(TString::kBoth);
   // If the user chooses ".help" or ".?".
   if ((strippedCommand == ".help") || (strippedCommand == ".?")) {
      gInterpreter->ProcessLine(line);
      Printf("\nROOT special commands.");
      Printf("==========================================================================");
      Printf("   .pwd                : show current directory, pad and style");
      Printf("   .ls                 : list contents of current directory");
      Printf("   .which [file]       : shows path of macro file");
      Printf("   .help Class         : opens the reference guide for that class");
      Printf("   .help Class::Member : opens the reference guide for function/member");
      return;
   } else {
      // If the user wants to use the extended ".help scopeName" command to access
      // the online reference guide, we first check if the command starts correctly.
      if ((!strippedCommand.BeginsWith(".help ")) && (!strippedCommand.BeginsWith(".? "))) {
         Error("Help", "Unknown command!");
         return;
      }
      // We remove the command ".help" or ".?" from the TString.
      if (strippedCommand.BeginsWith(".? ")) {
         strippedCommand.Remove(0, 3);
      } else {
         strippedCommand.Remove(0, 5);
      }
      // We strip the command line after removing ".help" or ".?".
      strippedCommand = strippedCommand.Strip(TString::kBoth);
      // We call the function what handles the extended ".help scopeName" command.
      OpenReferenceGuideFor(strippedCommand);
   }
}

/// Load shared libs necessary for graphics. These libraries are only
/// loaded when gROOT->IsBatch() is kFALSE.

void TApplication::LoadGraphicsLibs()
{
   if (gROOT->IsBatch()) return;

   TPluginHandler *h;
   if ((h = gROOT->GetPluginManager()->FindHandler("TVirtualPad")))
      if (h->LoadPlugin() == -1)
         return;

   TString name;
   TString title1 = "ROOT interface to ";
   TString nativex, title;
   TString nativeg = "root";

#ifdef R__WIN32
   nativex = "win32gdk";
   name    = "Win32gdk";
   title   = title1 + "Win32gdk";
#elif defined(R__HAS_COCOA)
   nativex = "quartz";
   name    = "quartz";
   title   = title1 + "Quartz";
#else
   nativex = "x11";
   name    = "X11";
   title   = title1 + "X11";
#endif

   TString guiBackend(gEnv->GetValue("Gui.Backend", "native"));
   guiBackend.ToLower();
   if (guiBackend == "native") {
      guiBackend = nativex;
   } else {
      name  = guiBackend;
      title = title1 + guiBackend;
   }
   TString guiFactory(gEnv->GetValue("Gui.Factory", "native"));
   guiFactory.ToLower();
   if (guiFactory == "native")
      guiFactory = nativeg;

   if ((h = gROOT->GetPluginManager()->FindHandler("TVirtualX", guiBackend))) {
      if (h->LoadPlugin() == -1) {
         gROOT->SetBatch(kTRUE);
         return;
      }
      gVirtualX = (TVirtualX *) h->ExecPlugin(2, name.Data(), title.Data());
      fgGraphInit = kTRUE;
   }
   if ((h = gROOT->GetPluginManager()->FindHandler("TGuiFactory", guiFactory))) {
      if (h->LoadPlugin() == -1) {
         gROOT->SetBatch(kTRUE);
         return;
      }
      gGuiFactory = (TGuiFactory *) h->ExecPlugin(0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Switch to batch mode.

void TApplication::MakeBatch()
{
   gROOT->SetBatch();
   if (gGuiFactory != gBatchGuiFactory) delete gGuiFactory;
   gGuiFactory = gBatchGuiFactory;
#ifndef R__WIN32
   if (gVirtualX != gGXBatch) delete gVirtualX;
#endif
   gVirtualX = gGXBatch;
}

////////////////////////////////////////////////////////////////////////////////
/// Parse the content of a line starting with ".R" (already stripped-off)
/// The format is
/// ~~~ {.cpp}
///      [user@]host[:dir] [-l user] [-d dbg] [script]
/// ~~~
/// The variable 'dir' is the remote directory to be used as working dir.
/// The username can be specified in two ways, "-l" having the priority
/// (as in ssh).
/// A 'dbg' value > 0 gives increasing verbosity.
/// The last argument 'script' allows to specify an alternative script to
/// be executed remotely to startup the session.

Int_t TApplication::ParseRemoteLine(const char *ln,
                                   TString &hostdir, TString &user,
                                   Int_t &dbg, TString &script)
{
   if (!ln || strlen(ln) <= 0)
      return 0;

   Int_t rc = 0;
   Bool_t isHostDir = kTRUE;
   Bool_t isScript = kFALSE;
   Bool_t isUser = kFALSE;
   Bool_t isDbg = kFALSE;

   TString line(ln);
   TString tkn;
   Int_t from = 0;
   while (line.Tokenize(tkn, from, " ")) {
      if (tkn == "-l") {
         // Next is a user name
         isUser = kTRUE;
      } else if (tkn == "-d") {
         isDbg = kTRUE;
      } else if (tkn == "-close") {
         rc = 1;
      } else if (tkn.BeginsWith("-")) {
         ::Warning("TApplication::ParseRemoteLine","unknown option: %s", tkn.Data());
      } else {
         if (isUser) {
            user = tkn;
            isUser = kFALSE;
         } else if (isDbg) {
            dbg = tkn.Atoi();
            isDbg = kFALSE;
         } else if (isHostDir) {
            hostdir = tkn;
            hostdir.ReplaceAll(":","/");
            isHostDir = kFALSE;
            isScript = kTRUE;
         } else if (isScript) {
            // Add everything left
            script = tkn;
            script.Insert(0, "\"");
            script += "\"";
            isScript = kFALSE;
            break;
         }
      }
   }

   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Process the content of a line starting with ".R" (already stripped-off)
/// The format is
/// ~~~ {.cpp}
///      [user@]host[:dir] [-l user] [-d dbg] [script] | [host] -close
/// ~~~
/// The variable 'dir' is the remote directory to be used as working dir.
/// The username can be specified in two ways, "-l" having the priority
/// (as in ssh).
/// A 'dbg' value > 0 gives increasing verbosity.
/// The last argument 'script' allows to specify an alternative script to
/// be executed remotely to startup the session.

Long_t TApplication::ProcessRemote(const char *line, Int_t *)
{
   if (!line) return 0;

   if (!strncmp(line, "-?", 2) || !strncmp(line, "-h", 2) ||
       !strncmp(line, "--help", 6)) {
      Info("ProcessRemote", "remote session help:");
      Printf(".R [user@]host[:dir] [-l user] [-d dbg] [[<]script] | [host] -close");
      Printf("Create a ROOT session on the specified remote host.");
      Printf("The variable \"dir\" is the remote directory to be used as working dir.");
      Printf("The username can be specified in two ways, \"-l\" having the priority");
      Printf("(as in ssh). A \"dbg\" value > 0 gives increasing verbosity.");
      Printf("The last argument \"script\" allows to specify an alternative script to");
      Printf("be executed remotely to startup the session, \"roots\" being");
      Printf("the default. If the script is preceded by a \"<\" the script will be");
      Printf("sourced, after which \"roots\" is executed. The sourced script can be ");
      Printf("used to change the PATH and other variables, allowing an alternative");
      Printf("\"roots\" script to be found.");
      Printf("To close down a session do \".R host -close\".");
      Printf("To switch between sessions do \".R host\", to switch to the local");
      Printf("session do \".R\".");
      Printf("To list all open sessions do \"gApplication->GetApplications()->Print()\".");
      return 0;
   }

   TString hostdir, user, script;
   Int_t dbg = 0;
   Int_t rc = ParseRemoteLine(line, hostdir, user, dbg, script);
   if (hostdir.Length() <= 0) {
      // Close the remote application if required
      if (rc == 1) {
         TApplication::Close(fAppRemote);
         delete fAppRemote;
      }
      // Return to local run
      fAppRemote = 0;
      // Done
      return 1;
   } else if (rc == 1) {
      // close an existing remote application
      TApplication *ap = TApplication::Open(hostdir, 0, 0);
      if (ap) {
         TApplication::Close(ap);
         delete ap;
      }
   }
   // Attach or start a remote application
   if (user.Length() > 0)
      hostdir.Insert(0,Form("%s@", user.Data()));
   const char *sc = (script.Length() > 0) ? script.Data() : 0;
   TApplication *ap = TApplication::Open(hostdir, dbg, sc);
   if (ap) {
      fAppRemote = ap;
   }

   // Done
   return 1;
}

namespace {
   static int PrintFile(const char* filename) {
      TString sFileName(filename);
      gSystem->ExpandPathName(sFileName);
      if (gSystem->AccessPathName(sFileName)) {
         Error("ProcessLine()", "Cannot find file %s", filename);
         return 1;
      }
      std::ifstream instr(sFileName);
      TString content;
      content.ReadFile(instr);
      Printf("%s", content.Data());
      return 0;
   }
   } // namespace

////////////////////////////////////////////////////////////////////////////////
/// Process a single command line, either a C++ statement or an interpreter
/// command starting with a ".".
/// Return the return value of the command cast to a long.

Long_t TApplication::ProcessLine(const char *line, Bool_t sync, Int_t *err)
{
   if (!line || !*line) return 0;

   // If we are asked to go remote do it
   if (!strncmp(line, ".R", 2)) {
      Int_t n = 2;
      while (*(line+n) == ' ')
         n++;
      return ProcessRemote(line+n, err);
   }

   // Redirect, if requested
   if (fAppRemote && TestBit(kProcessRemotely)) {
      ResetBit(kProcessRemotely);
      return fAppRemote->ProcessLine(line, err);
   }

   if (!strncasecmp(line, ".qqqqqqq", 7)) {
      gSystem->Abort();
   } else if (!strncasecmp(line, ".qqqqq", 5)) {
      Info("ProcessLine", "Bye... (try '.qqqqqqq' if still running)");
      gSystem->Exit(1);
   } else if (!strncasecmp(line, ".exit", 4) || !strncasecmp(line, ".quit", 2)) {
      Terminate(0);
      return 0;
   }

   if (!strncmp(line, ".?", 2) || !strncmp(line, ".help", 5)) {
      Help(line);
      return 1;
   }

   if (!strncmp(line, ".demo", 5)) {
      if (gROOT->IsBatch()) {
         Error("ProcessLine", "Cannot show demos in batch mode!");
         return 1;
      }
      ProcessLine(".x " + TROOT::GetTutorialDir() + "/demos.C");
      return 0;
   }

   if (!strncmp(line, ".license", 8)) {
      return PrintFile(TROOT::GetDocDir() + "/LICENSE");
   }

   if (!strncmp(line, ".credits", 8)) {
      TString credits = TROOT::GetDocDir() + "/CREDITS";
      if (gSystem->AccessPathName(credits, kReadPermission))
         credits = TROOT::GetDocDir() + "/README/CREDITS";
      return PrintFile(credits);
   }

   if (!strncmp(line, ".pwd", 4)) {
      if (gDirectory)
         Printf("Current directory: %s", gDirectory->GetPath());
      if (gPad)
         Printf("Current pad:       %s", gPad->GetName());
      if (gStyle)
         Printf("Current style:     %s", gStyle->GetName());
      return 1;
   }

   if (!strncmp(line, ".ls", 3)) {
      const char *opt = 0;
      if (line[3]) opt = &line[3];
      if (gDirectory) gDirectory->ls(opt);
      return 1;
   }

   if (!strncmp(line, ".which", 6)) {
      char *fn  = Strip(line+7);
      char *s = strtok(fn, "+("); // this method does not need to be reentrant
      char *mac = gSystem->Which(TROOT::GetMacroPath(), s, kReadPermission);
      if (!mac)
         Printf("No macro %s in path %s", s, TROOT::GetMacroPath());
      else
         Printf("%s", mac);
      delete [] fn;
      delete [] mac;
      return mac ? 1 : 0;
   }

   if (!strncmp(line, ".L", 2) || !strncmp(line, ".U", 2)) {
      TString aclicMode;
      TString arguments;
      TString io;
      TString fname = gSystem->SplitAclicMode(line+3, aclicMode, arguments, io);

      char *mac = gSystem->Which(TROOT::GetMacroPath(), fname, kReadPermission);
      if (arguments.Length()) {
         Warning("ProcessLine", "argument(s) \"%s\" ignored with .%c", arguments.Data(),
                 line[1]);
      }
      Long_t retval = 0;
      if (!mac)
         Error("ProcessLine", "macro %s not found in path %s", fname.Data(),
               TROOT::GetMacroPath());
      else {
         TString cmd(line+1);
         Ssiz_t posSpace = cmd.Index(' ');
         if (posSpace == -1) cmd.Remove(1);
         else cmd.Remove(posSpace);
         TString tempbuf;
         if (sync) {
            tempbuf.Form(".%s %s%s%s", cmd.Data(), mac, aclicMode.Data(),io.Data());
            retval = gInterpreter->ProcessLineSynch(tempbuf,
                                                   (TInterpreter::EErrorCode*)err);
         } else {
            tempbuf.Form(".%s %s%s%s", cmd.Data(), mac, aclicMode.Data(),io.Data());
            retval = gInterpreter->ProcessLine(tempbuf,
                                              (TInterpreter::EErrorCode*)err);
         }
      }

      delete [] mac;

      InitializeGraphics();

      return retval;
   }

   if (!strncmp(line, ".X", 2) || !strncmp(line, ".x", 2)) {
      return ProcessFile(line+3, err, line[2] == 'k');
   }

   if (!strcmp(line, ".reset")) {
      // Do nothing, .reset disabled in CINT because too many side effects
      Printf("*** .reset not allowed, please use gROOT->Reset() ***");
      return 0;

#if 0
      // delete the ROOT dictionary since CINT will destroy all objects
      // referenced by the dictionary classes (TClass et. al.)
      gROOT->GetListOfClasses()->Delete();
      // fall through
#endif
   }

   if (sync)
      return gInterpreter->ProcessLineSynch(line, (TInterpreter::EErrorCode*)err);
   else
      return gInterpreter->ProcessLine(line, (TInterpreter::EErrorCode*)err);
}

////////////////////////////////////////////////////////////////////////////////
/// Process a file containing a C++ macro.

Long_t TApplication::ProcessFile(const char *file, Int_t *error, Bool_t keep)
{
   return ExecuteFile(file, error, keep);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute a file containing a C++ macro (static method). Can be used
/// while TApplication is not yet created.

Long_t TApplication::ExecuteFile(const char *file, Int_t *error, Bool_t keep)
{
   static const Int_t kBufSize = 1024;

   if (!file || !*file) return 0;

   TString aclicMode;
   TString arguments;
   TString io;
   TString fname = gSystem->SplitAclicMode(file, aclicMode, arguments, io);

   char *exnam = gSystem->Which(TROOT::GetMacroPath(), fname, kReadPermission);
   if (!exnam) {
      ::Error("TApplication::ExecuteFile", "macro %s not found in path %s", fname.Data(),
              TROOT::GetMacroPath());
      delete [] exnam;
      if (error)
         *error = (Int_t)TInterpreter::kRecoverable;
      return 0;
   }

   ::std::ifstream macro(exnam, std::ios::in);
   if (!macro.good()) {
      ::Error("TApplication::ExecuteFile", "%s no such file", exnam);
      if (error)
         *error = (Int_t)TInterpreter::kRecoverable;
      delete [] exnam;
      return 0;
   }

   char currentline[kBufSize];
   char dummyline[kBufSize];
   int tempfile = 0;
   int comment  = 0;
   int ifndefc  = 0;
   int ifdef    = 0;
   char *s      = 0;
   Bool_t execute = kFALSE;
   Long_t retval = 0;

   while (1) {
      bool res = (bool)macro.getline(currentline, kBufSize);
      if (macro.eof()) break;
      if (!res) {
         // Probably only read kBufSize, let's ignore the remainder of
         // the line.
         macro.clear();
         while (!macro.getline(dummyline, kBufSize) && !macro.eof()) {
            macro.clear();
         }
      }
      s = currentline;
      while (s && (*s == ' ' || *s == '\t')) s++;   // strip-off leading blanks

      // very simple minded pre-processor parsing, only works in case macro file
      // starts with "#ifndef __CINT__". In that case everything till next
      // "#else" or "#endif" will be skipped.
      if (*s == '#') {
         char *cs = Compress(currentline);
         if (strstr(cs, "#ifndef__CINT__") ||
             strstr(cs, "#if!defined(__CINT__)"))
            ifndefc = 1;
         else if (ifndefc && (strstr(cs, "#ifdef") || strstr(cs, "#ifndef") ||
                  strstr(cs, "#ifdefined") || strstr(cs, "#if!defined")))
            ifdef++;
         else if (ifndefc && strstr(cs, "#endif")) {
            if (ifdef)
               ifdef--;
            else
               ifndefc = 0;
         } else if (ifndefc && !ifdef && strstr(cs, "#else"))
            ifndefc = 0;
         delete [] cs;
      }
      if (!*s || *s == '#' || ifndefc || !strncmp(s, "//", 2)) continue;

      if (!comment && (!strncmp(s, ".X", 2) || !strncmp(s, ".x", 2))) {
         retval = ExecuteFile(s+3);
         execute = kTRUE;
         continue;
      }

      if (!strncmp(s, "/*", 2)) comment = 1;
      if (comment) {
         // handle slightly more complex cases like: /*  */  /*
again:
         s = strstr(s, "*/");
         if (s) {
            comment = 0;
            s += 2;

            while (s && (*s == ' ' || *s == '\t')) s++; // strip-off leading blanks
            if (!*s) continue;
            if (!strncmp(s, "//", 2)) continue;
            if (!strncmp(s, "/*", 2)) {
               comment = 1;
               goto again;
            }
         }
      }
      if (!comment && *s == '{') tempfile = 1;
      if (!comment) break;
   }
   macro.close();

   if (!execute) {
      TString exname = exnam;
      if (!tempfile) {
         // We have a script that does NOT contain an unnamed macro,
         // so we can call the script compiler on it.
         exname += aclicMode;
      }
      exname += arguments;
      exname += io;

      TString tempbuf;
      if (tempfile) {
         tempbuf.Form(".x %s", exname.Data());
      } else {
         tempbuf.Form(".X%s %s", keep ? "k" : " ", exname.Data());
      }
      retval = gInterpreter->ProcessLineSynch(tempbuf,(TInterpreter::EErrorCode*)error);
   }

   delete [] exnam;
   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Main application eventloop. Calls system dependent eventloop via gSystem.

void TApplication::Run(Bool_t retrn)
{
   SetReturnFromRun(retrn);

   fIsRunning = kTRUE;

   gSystem->Run();
   fIsRunning = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the command to be executed after the system has been idle for
/// idleTimeInSec seconds. Normally called via TROOT::Idle(...).

void TApplication::SetIdleTimer(UInt_t idleTimeInSec, const char *command)
{
   if (fIdleTimer) RemoveIdleTimer();
   fIdleCommand = command;
   fIdleTimer = new TIdleTimer(idleTimeInSec*1000);
   gSystem->AddTimer(fIdleTimer);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove idle timer. Normally called via TROOT::Idle(0).

void TApplication::RemoveIdleTimer()
{
   if (fIdleTimer) {
      // timers are removed from the gSystem timer list by their dtor
      SafeDelete(fIdleTimer);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called when system starts idleing.

void TApplication::StartIdleing()
{
   if (fIdleTimer) {
      fIdleTimer->Reset();
      gSystem->AddTimer(fIdleTimer);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called when system stops idleing.

void TApplication::StopIdleing()
{
   if (fIdleTimer)
      gSystem->RemoveTimer(fIdleTimer);
}

////////////////////////////////////////////////////////////////////////////////
/// What to do when tab is pressed. Re-implemented by TRint.
/// See TTabCom::Hook() for meaning of return values.

Int_t TApplication::TabCompletionHook(char* /*buf*/, int* /*pLoc*/, std::ostream& /*out*/)
{
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Terminate the application by call TSystem::Exit() unless application has
/// been told to return from Run(), by a call to SetReturnFromRun().

void TApplication::Terminate(Int_t status)
{
   Emit("Terminate(Int_t)", status);

   if (fReturnFromRun)
      gSystem->ExitLoop();
   else {
      //close TMemStat
      if (fUseMemstat) {
         ProcessLine("TMemStat::Close()");
         fUseMemstat = kFALSE;
      }

      gSystem->Exit(status);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Emit signal when a line has been processed.

void TApplication::LineProcessed(const char *line)
{
   Emit("LineProcessed(const char*)", line);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit signal when console keyboard key was pressed.

void TApplication::KeyPressed(Int_t key)
{
   Emit("KeyPressed(Int_t)", key);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit signal when return key was pressed.

void TApplication::ReturnPressed(char *text )
{
   Emit("ReturnPressed(char*)", text);
}

////////////////////////////////////////////////////////////////////////////////
/// Set console echo mode:
///
///  - mode = kTRUE  - echo input symbols
///  - mode = kFALSE - noecho input symbols

void TApplication::SetEchoMode(Bool_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Static function used to create a default application environment.

void TApplication::CreateApplication()
{
   R__LOCKGUARD(gROOTMutex);
   // gApplication is set at the end of 'new TApplication.
   if (!gApplication) {
      char *a = StrDup("RootApp");
      char *b = StrDup("-b");
      char *argv[2];
      Int_t argc = 2;
      argv[0] = a;
      argv[1] = b;
      new TApplication("RootApp", &argc, argv, 0, 0);
      if (gDebug > 0)
         Printf("<TApplication::CreateApplication>: "
                "created default TApplication");
      delete [] a; delete [] b;
      gApplication->SetBit(kDefaultApplication);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Static function used to attach to an existing remote application
/// or to start one.

TApplication *TApplication::Open(const char *url,
                                  Int_t debug, const char *script)
{
   TApplication *ap = 0;
   TUrl nu(url);
   Int_t nnew = 0;

   // Look among the existing ones
   if (fgApplications) {
      TIter nxa(fgApplications);
      while ((ap = (TApplication *) nxa())) {
         TString apn(ap->ApplicationName());
         if (apn == url) {
            // Found matching application
            return ap;
         } else {
            // Check if same machine and user
            TUrl au(apn);
            if (strlen(au.GetUser()) > 0 && strlen(nu.GetUser()) > 0 &&
                !strcmp(au.GetUser(), nu.GetUser())) {
               if (!strncmp(au.GetHost(), nu.GetHost(), strlen(nu.GetHost())))
                  // New session on a known machine
                  nnew++;
            }
         }
      }
   } else {
      ::Error("TApplication::Open", "list of applications undefined - protocol error");
      return ap;
   }

   // If new session on a known machine pass the number as option
   if (nnew > 0) {
      nnew++;
      nu.SetOptions(Form("%d", nnew));
   }

   // Instantiate the TApplication object to be run
   TPluginHandler *h = 0;
   if ((h = gROOT->GetPluginManager()->FindHandler("TApplication","remote"))) {
      if (h->LoadPlugin() == 0) {
         ap = (TApplication *) h->ExecPlugin(3, nu.GetUrl(), debug, script);
      } else {
         ::Error("TApplication::Open", "failed to load plugin for TApplicationRemote");
      }
   } else {
      ::Error("TApplication::Open", "failed to find plugin for TApplicationRemote");
   }

   // Add to the list
   if (ap && !(ap->TestBit(kInvalidObject))) {
      fgApplications->Add(ap);
      gROOT->GetListOfBrowsables()->Add(ap, ap->ApplicationName());
      TIter next(gROOT->GetListOfBrowsers());
      TBrowser *b;
      while ((b = (TBrowser*) next()))
         b->Add(ap, ap->ApplicationName());
      gROOT->RefreshBrowsers();
   } else {
      SafeDelete(ap);
      ::Error("TApplication::Open",
              "TApplicationRemote for %s could not be instantiated", url);
   }

   // Done
   return ap;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function used to close a remote application

void TApplication::Close(TApplication *app)
{
   if (app) {
      app->Terminate(0);
      fgApplications->Remove(app);
      gROOT->GetListOfBrowsables()->RecursiveRemove(app);
      TIter next(gROOT->GetListOfBrowsers());
      TBrowser *b;
      while ((b = (TBrowser*) next()))
         b->RecursiveRemove(app);
      gROOT->RefreshBrowsers();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Show available sessions

void TApplication::ls(Option_t *opt) const
{
   if (fgApplications) {
      TIter nxa(fgApplications);
      TApplication *a = 0;
      while ((a = (TApplication *) nxa())) {
         a->Print(opt);
      }
   } else {
      Print(opt);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the list of available applications

TList *TApplication::GetApplications()
{
   return fgApplications;
}
