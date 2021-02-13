// @(#)root/rint:$Id$
// Author: Rene Brun   17/02/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Rint                                                                 //
//                                                                      //
// Rint is the ROOT Interactive Interface. It allows interactive access //
// to the ROOT system via the Cling C/C++ interpreter.                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TClass.h"
#include "TClassEdit.h"
#include "TVirtualX.h"
#include "TObjectTable.h"
#include "TClassTable.h"
#include "TStopwatch.h"
#include "TBenchmark.h"
#include "TRint.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TSysEvtHandler.h"
#include "TSystemDirectory.h"
#include "TError.h"
#include "TException.h"
#include "TInterpreter.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TStorage.h" // ROOT::Internal::gMmallocDesc
#include "ThreadLocalStorage.h"
#include "TTabCom.h"
#include <cstdlib>
#include <algorithm>

#include "Getline.h"
#include "strlcpy.h"
#include "snprintf.h"

#ifdef R__UNIX
#include <signal.h>
#endif

////////////////////////////////////////////////////////////////////////////////

static Int_t Key_Pressed(Int_t key)
{
   gApplication->KeyPressed(key);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

static Int_t BeepHook()
{
   if (!gSystem) return 0;
   gSystem->Beep();
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Restore terminal to non-raw mode.

static void ResetTermAtExit()
{
   Getlinem(kCleanUp, 0);
}


//----- Interrupt signal handler -----------------------------------------------
////////////////////////////////////////////////////////////////////////////////

class TInterruptHandler : public TSignalHandler {
public:
   TInterruptHandler() : TSignalHandler(kSigInterrupt, kFALSE) { }
   Bool_t  Notify();
};

////////////////////////////////////////////////////////////////////////////////
/// TRint interrupt handler.

Bool_t TInterruptHandler::Notify()
{
   if (fDelay) {
      fDelay++;
      return kTRUE;
   }

   // make sure we use the sbrk heap (in case of mapped files)
   ROOT::Internal::gMmallocDesc = 0;

   if (TROOT::Initialized() && gROOT->IsLineProcessing()) {
      Break("TInterruptHandler::Notify", "keyboard interrupt");
      Getlinem(kInit, "Root > ");
      gCling->Reset();
#ifndef WIN32
   if (gException)
      Throw(GetSignal());
#endif
   } else {
      // Reset input.
      Getlinem(kClear, ((TRint*)gApplication)->GetPrompt());
   }

   return kTRUE;
}

//----- Terminal Input file handler --------------------------------------------
////////////////////////////////////////////////////////////////////////////////

class TTermInputHandler : public TFileHandler {
public:
   TTermInputHandler(Int_t fd) : TFileHandler(fd, 1) { }
   Bool_t Notify();
   Bool_t ReadNotify() { return Notify(); }
};

////////////////////////////////////////////////////////////////////////////////
/// Notify implementation.  Call the application interupt handler.

Bool_t TTermInputHandler::Notify()
{
   return gApplication->HandleTermInput();
}


ClassImp(TRint);

////////////////////////////////////////////////////////////////////////////////
/// Create an application environment. The TRint environment provides an
/// interface to the WM manager functionality and eventloop via inheritance
/// of TApplication and in addition provides interactive access to
/// the Cling C++ interpreter via the command line.

TRint::TRint(const char *appClassName, Int_t *argc, char **argv, void *options,
             Int_t numOptions, Bool_t noLogo):
   TApplication(appClassName, argc, argv, options, numOptions),
   fCaughtSignal(-1)
{
   fNcmd          = 0;
   fDefaultPrompt = "root [%d] ";
   fInterrupt     = kFALSE;

   gBenchmark = new TBenchmark();

   if (!noLogo && !NoLogoOpt()) {
      Bool_t lite = (Bool_t) gEnv->GetValue("Rint.WelcomeLite", 0);
      PrintLogo(lite);
   }

   // Explicitly load libMathCore it cannot be auto-loaded it when using one
   // of its freestanding functions. Once functions can trigger autoloading we
   // can get rid of this.
   if (!gClassTable->GetDict("TRandom"))
      gSystem->Load("libMathCore");

   if (!gInterpreter->HasPCMForLibrary("std")) {
      // Load some frequently used includes
      Int_t includes = gEnv->GetValue("Rint.Includes", 1);
      // When the interactive ROOT starts, it can automatically load some frequently
      // used includes. However, this introduces several overheads
      //   -The initialisation takes more time
      //   -Memory overhead when including <vector>
      // In $ROOTSYS/etc/system.rootrc, you can set the variable Rint.Includes to 0
      // to disable the loading of these includes at startup.
      // You can set the variable to 1 (default) to load only <iostream>, <string> and <DllImport.h>
      // You can set it to 2 to load in addition <vector> and <utility>
      // We strongly recommend setting the variable to 2 if your scripts include <vector>
      // and you execute your scripts multiple times.
      if (includes > 0) {
         TString code;
         code = "#include <iostream>\n"
            "#include <string>\n" // for std::string std::iostream.
            "#include <DllImport.h>\n";// Defined R__EXTERN
         if (includes > 1) {
            code += "#include <vector>\n"
               "#include <utility>";
         }
         ProcessLine(code, kTRUE);
      }
   }

   // Load user functions
   const char *logon;
   logon = gEnv->GetValue("Rint.Load", (char*)0);
   if (logon) {
      char *mac = gSystem->Which(TROOT::GetMacroPath(), logon, kReadPermission);
      if (mac)
         ProcessLine(Form(".L %s",logon), kTRUE);
      delete [] mac;
   }

   // Execute logon macro
   ExecLogon();

   // Save current interpreter context
   gCling->SaveContext();
   gCling->SaveGlobalsContext();

   // Install interrupt and terminal input handlers
   TInterruptHandler *ih = new TInterruptHandler();
   ih->Add();
   SetSignalHandler(ih);

   // Handle stdin events
   fInputHandler = new TTermInputHandler(0);
   fInputHandler->Add();

   // Goto into raw terminal input mode
   char defhist[kMAXPATHLEN];
   snprintf(defhist, sizeof(defhist), "%s/.root_hist", gSystem->HomeDirectory());
   logon = gEnv->GetValue("Rint.History", defhist);
   // In the code we had HistorySize and HistorySave, in the rootrc and doc
   // we have HistSize and HistSave. Keep the doc as it is and check
   // now also for HistSize and HistSave in case the user did not use
   // the History versions
   int hist_size = gEnv->GetValue("Rint.HistorySize", 500);
   if (hist_size == 500)
      hist_size = gEnv->GetValue("Rint.HistSize", 500);
   int hist_save = gEnv->GetValue("Rint.HistorySave", 400);
   if (hist_save == 400)
      hist_save = gEnv->GetValue("Rint.HistSave", 400);
   const char *envHist = gSystem->Getenv("ROOT_HIST");
   if (envHist) {
      hist_size = atoi(envHist);
      envHist = strchr(envHist, ':');
      if (envHist)
         hist_save = atoi(envHist+1);
   }
   Gl_histsize(hist_size, hist_save);
   Gl_histinit((char *)logon);

   // black on white or white on black?
   static const char* defaultColorsBW[] = {
      "bold blue", "magenta", "bold green", "bold red underlined", "default"
   };
   static const char* defaultColorsWB[] = {
      "yellow", "magenta", "bold green", "bold red underlined", "default"
   };

   const char** defaultColors = defaultColorsBW;
   TString revColor = gEnv->GetValue("Rint.ReverseColor", "no");
   if (revColor.Contains("yes", TString::kIgnoreCase)) {
      defaultColors = defaultColorsWB;
   }
   TString colorType = gEnv->GetValue("Rint.TypeColor", defaultColors[0]);
   TString colorTabCom = gEnv->GetValue("Rint.TabComColor", defaultColors[1]);
   TString colorBracket = gEnv->GetValue("Rint.BracketColor", defaultColors[2]);
   TString colorBadBracket = gEnv->GetValue("Rint.BadBracketColor", defaultColors[3]);
   TString colorPrompt = gEnv->GetValue("Rint.PromptColor", defaultColors[4]);
   Gl_setColors(colorType, colorTabCom, colorBracket, colorBadBracket, colorPrompt);

   Gl_windowchanged();

   atexit(ResetTermAtExit);

   // Setup for tab completion
   gTabCom      = new TTabCom;
   Gl_in_key    = &Key_Pressed;
   Gl_beep_hook = &BeepHook;

   // tell Cling to use our getline
   gCling->SetGetline(Getline, Gl_histadd);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TRint::~TRint()
{
   delete gTabCom;
   gTabCom = 0;
   Gl_in_key = 0;
   Gl_beep_hook = 0;
   fInputHandler->Remove();
   delete fInputHandler;
   // We can't know where the signal handler was changed since we started ...
   // so for now let's not delete it.
//   TSignalHandler *ih  = GetSignalHandler();
//   ih->Remove();
//   SetSignalHandler(0);
//   delete ih;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute logon macro's. There are three levels of logon macros that
/// will be executed: the system logon etc/system.rootlogon.C, the global
/// user logon ~/.rootlogon.C and the local ./.rootlogon.C. For backward
/// compatibility also the logon macro as specified by the Rint.Logon
/// environment setting, by default ./rootlogon.C, will be executed.
/// No logon macros will be executed when the system is started with
/// the -n option.

void TRint::ExecLogon()
{
   if (NoLogOpt()) return;

   TString name = ".rootlogon.C";
   TString sname = "system";
   sname += name;
   char *s = gSystem->ConcatFileName(TROOT::GetEtcDir(), sname);
   if (!gSystem->AccessPathName(s, kReadPermission)) {
      ProcessFile(s);
   }
   delete [] s;
   s = gSystem->ConcatFileName(gSystem->HomeDirectory(), name);
   if (!gSystem->AccessPathName(s, kReadPermission)) {
      ProcessFile(s);
   }
   delete [] s;
   // avoid executing ~/.rootlogon.C twice
   if (strcmp(gSystem->HomeDirectory(), gSystem->WorkingDirectory())) {
      if (!gSystem->AccessPathName(name, kReadPermission))
         ProcessFile(name);
   }

   // execute also the logon macro specified by "Rint.Logon"
   const char *logon = gEnv->GetValue("Rint.Logon", (char*)0);
   if (logon) {
      char *mac = gSystem->Which(TROOT::GetMacroPath(), logon, kReadPermission);
      if (mac)
         ProcessFile(logon);
      delete [] mac;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Main application eventloop. First process files given on the command
/// line and then go into the main application event loop, unless the -q
/// command line option was specified in which case the program terminates.
/// When return is true this method returns even when -q was specified.
///
/// When QuitOpt is true and return is false, terminate the application with
/// an error code equal to either the ProcessLine error (if any) or the
/// return value of the command casted to a long.

void TRint::Run(Bool_t retrn)
{
   if (!QuitOpt()) {
      // Promt prompt only if we are expecting / allowing input.
      Getlinem(kInit, GetPrompt());
   }

   Long_t retval = 0;
   Int_t  error = 0;
   volatile Bool_t needGetlinemInit = kFALSE;

   if (strlen(WorkingDirectory())) {
      // if directory specified as argument make it the working directory
      gSystem->ChangeDirectory(WorkingDirectory());
      TSystemDirectory *workdir = new TSystemDirectory("workdir", gSystem->WorkingDirectory());
      TObject *w = gROOT->GetListOfBrowsables()->FindObject("workdir");
      TObjLink *lnk = gROOT->GetListOfBrowsables()->FirstLink();
      while (lnk) {
         if (lnk->GetObject() == w) {
            lnk->SetObject(workdir);
            lnk->SetOption(gSystem->WorkingDirectory());
            break;
         }
         lnk = lnk->Next();
      }
      delete w;
   }

   // Process shell command line input files
   if (InputFiles()) {
      // Make sure that calls into the event loop
      // ignore end-of-file on the terminal.
      fInputHandler->DeActivate();
      TIter next(InputFiles());
      RETRY {
         retval = 0; error = 0;
         Int_t nfile = 0;
         while (TObject *fileObj = next()) {
            if (dynamic_cast<TNamed*>(fileObj)) {
               // A file that TApplication did not find. Note the error.
               retval = 1;
               continue;
            }
            TObjString *file = (TObjString *)fileObj;
            char cmd[kMAXPATHLEN+50];
            if (!fNcmd)
               printf("\n");
            Bool_t rootfile = kFALSE;

            if (file->TestBit(kExpression)) {
               snprintf(cmd, kMAXPATHLEN+50, "%s", (const char*)file->String());
            } else {
               if (file->String().EndsWith(".root") || file->String().BeginsWith("file:")) {
                  rootfile = kTRUE;
               } else {
                  rootfile = gROOT->IsRootFile(file->String());
               }
               if (rootfile) {
                  // special trick to be able to open files using UNC path names
                  if (file->String().BeginsWith("\\\\"))
                     file->String().Prepend("\\\\");
                  file->String().ReplaceAll("\\","/");
                  const char *rfile = (const char*)file->String();
                  Printf("Attaching file %s as _file%d...", rfile, nfile);
                  snprintf(cmd, kMAXPATHLEN+50, "TFile *_file%d = TFile::Open(\"%s\")", nfile++, rfile);
               } else {
                  Printf("Processing %s...", (const char*)file->String());
                  snprintf(cmd, kMAXPATHLEN+50, ".x %s", (const char*)file->String());
               }
            }
            Getlinem(kCleanUp, 0);
            Gl_histadd(cmd);
            fNcmd++;

            // The ProcessLine might throw an 'exception'.  In this case,
            // GetLinem(kInit,"Root >") is called and we are jump back
            // to RETRY ... and we have to avoid the Getlinem(kInit, GetPrompt());
            needGetlinemInit = kFALSE;
            retval = ProcessLineNr("ROOT_cli_", cmd, &error);
            gCling->EndOfLineAction();

            // The ProcessLine has successfully completed and we need
            // to call Getlinem(kInit, GetPrompt());
            needGetlinemInit = kTRUE;

            if (error != 0 || fCaughtSignal != -1) break;
         }
      } ENDTRY;

      if (QuitOpt()) {
         if (retrn) return;
         if (error) {
            retval = error;
         } else if (fCaughtSignal != -1) {
            retval = fCaughtSignal + 128;
         }
         // Bring retval into sensible range, 0..255.
         if (retval < 0 || retval > 255)
            retval = 255;
         Terminate(retval);
      }

      // Allow end-of-file on the terminal to be noticed
      // after we finish processing the command line input files.
      fInputHandler->Activate();

      ClearInputFiles();

      if (needGetlinemInit) Getlinem(kInit, GetPrompt());
   }

   if (QuitOpt()) {
      printf("\n");
      if (retrn) return;
      Terminate(fCaughtSignal != -1 ? fCaughtSignal + 128 : 0);
   }

   TApplication::Run(retrn);

   // Reset to happiness.
   fCaughtSignal = -1;

   Getlinem(kCleanUp, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Print the ROOT logo on standard output.

void TRint::PrintLogo(Bool_t lite)
{
   if (!lite) {
      // Fancy formatting: the content of lines are format strings; their %s is
      // replaced by spaces needed to make all lines as long as the longest line.
      std::vector<TString> lines;
      // Here, %%s results in %s after TString::Format():
      lines.emplace_back(TString::Format("Welcome to ROOT %s%%shttps://root.cern",
                                         gROOT->GetVersion()));
      lines.emplace_back(TString::Format("(c) 1995-2021, The ROOT Team; conception: R. Brun, F. Rademakers%%s"));
      lines.emplace_back(TString::Format("Built for %s on %s%%s", gSystem->GetBuildArch(), gROOT->GetGitDate()));
      if (!strcmp(gROOT->GetGitBranch(), gROOT->GetGitCommit())) {
         static const char *months[] = {"January","February","March","April","May",
                                        "June","July","August","September","October",
                                        "November","December"};
         Int_t idatqq = gROOT->GetVersionDate();
         Int_t iday   = idatqq%100;
         Int_t imonth = (idatqq/100)%100;
         Int_t iyear  = (idatqq/10000);

         lines.emplace_back(TString::Format("From tag %s, %d %s %4d%%s",
                                            gROOT->GetGitBranch(),
                                            iday,months[imonth-1],iyear));
      } else {
         // If branch and commit are identical - e.g. "v5-34-18" - then we have
         // a release build. Else specify the git hash this build was made from.
         lines.emplace_back(TString::Format("From %s@%s %%s",
                                            gROOT->GetGitBranch(),
                                            gROOT->GetGitCommit()));
      }
      lines.emplace_back(TString::Format("With %s %%s",
                                         gSystem->GetBuildCompilerVersionStr()));
      lines.emplace_back(TString("Try '.help', '.demo', '.license', '.credits', '.quit'/'.q'%s"));

      // Find the longest line and its length:
      auto itLongest = std::max_element(lines.begin(), lines.end(),
                                        [](const TString& left, const TString& right) {
                                           return left.Length() < right.Length(); });
      Ssiz_t lenLongest = itLongest->Length();


      Printf("   %s", TString('-', lenLongest).Data());
      for (const auto& line: lines) {
         // Print the line, expanded with the necessary spaces at %s, and
         // surrounded by some ASCII art.
         Printf("  | %s |",
                TString::Format(line.Data(),
                                TString(' ', lenLongest - line.Length()).Data()).Data());
      }
      Printf("   %s\n", TString('-', lenLongest).Data());
   }

#ifdef R__UNIX
   // Popdown X logo, only if started with -splash option
   for (int i = 0; i < Argc(); i++)
      if (!strcmp(Argv(i), "-splash"))
         kill(getppid(), SIGUSR1);
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Get prompt from interpreter. Either "root [n]" or "end with '}'".

char *TRint::GetPrompt()
{
   char *s = gCling->GetPrompt();
   if (s[0])
      strlcpy(fPrompt, s, sizeof(fPrompt));
   else
      snprintf(fPrompt, sizeof(fPrompt), fDefaultPrompt.Data(), fNcmd);

   return fPrompt;
}

////////////////////////////////////////////////////////////////////////////////
/// Set a new default prompt. It returns the previous prompt.
/// The prompt may contain a %d which will be replaced by the commend
/// number. The default prompt is "root [%d] ". The maximum length of
/// the prompt is 55 characters. To set the prompt in an interactive
/// session do:
/// root [0] ((TRint*)gROOT->GetApplication())->SetPrompt("aap> ")
/// aap>

const char *TRint::SetPrompt(const char *newPrompt)
{
   static TString op;
   op = fDefaultPrompt;

   if (newPrompt && strlen(newPrompt) <= 55)
      fDefaultPrompt = newPrompt;
   else
      Error("SetPrompt", "newPrompt too long (> 55 characters)");

   return op.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle input coming from terminal.

Bool_t TRint::HandleTermInput()
{
   static TStopwatch timer;
   const char *line;

   if ((line = Getlinem(kOneChar, 0))) {
      if (line[0] == 0 && Gl_eof())
         Terminate(0);

      gVirtualX->SetKeyAutoRepeat(kTRUE);

      Gl_histadd(line);

      TString sline = line;

      // strip off '\n' and leading and trailing blanks
      sline = sline.Chop();
      sline = sline.Strip(TString::kBoth);
      ReturnPressed((char*)sline.Data());

      fInterrupt = kFALSE;

      if (!gCling->GetMore() && !sline.IsNull()) fNcmd++;

      // prevent recursive calling of this input handler
      fInputHandler->DeActivate();

      if (gROOT->Timer()) timer.Start();

      TTHREAD_TLS(Bool_t) added;
      added = kFALSE; // reset on each call.

      // This is needed when working with remote sessions
      SetBit(kProcessRemotely);

      try {
         TRY {
            if (!sline.IsNull())
               LineProcessed(sline);
            ProcessLineNr("ROOT_prompt_", sline);
         } CATCH(excode) {
            // enable again input handler
            fInputHandler->Activate();
            added = kTRUE;
            Throw(excode);
         } ENDTRY;
      }
      // handle every exception
      catch (std::exception& e) {
         // enable again intput handler
         if (!added) fInputHandler->Activate();

         int err;
         char *demangledType_c = TClassEdit::DemangleTypeIdName(typeid(e), err);
         const char* demangledType = demangledType_c;
         if (err) {
            demangledType_c = nullptr;
            demangledType = "<UNKNOWN>";
         }
         Error("HandleTermInput()", "%s caught: %s", demangledType, e.what());
         free(demangledType_c);
      }
      catch (...) {
         // enable again intput handler
         if (!added) fInputHandler->Activate();
         Error("HandleTermInput()", "Exception caught!");
      }

      if (gROOT->Timer()) timer.Print("u");

      // enable again intput handler
      fInputHandler->Activate();

      if (!sline.BeginsWith(".reset"))
         gCling->EndOfLineAction();

      gTabCom->ClearAll();
      Getlinem(kInit, GetPrompt());
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle signals (kSigBus, kSigSegmentationViolation,
/// kSigIllegalInstruction and kSigFloatingException) trapped in TSystem.
/// Specific TApplication implementations may want something different here.

void TRint::HandleException(Int_t sig)
{
   fCaughtSignal = sig;
   if (TROOT::Initialized()) {
      if (gException) {
         Getlinem(kCleanUp, 0);
         Getlinem(kInit, "Root > ");
      }
   }
   TApplication::HandleException(sig);
}

////////////////////////////////////////////////////////////////////////////////
/// Terminate the application. Reset the terminal to sane mode and call
/// the logoff macro defined via Rint.Logoff environment variable.

void TRint::Terminate(Int_t status)
{
   Getlinem(kCleanUp, 0);

   if (ReturnFromRun()) {
      gSystem->ExitLoop();
   } else {
      delete gTabCom;
      gTabCom = 0;

      //Execute logoff macro
      const char *logoff;
      logoff = gEnv->GetValue("Rint.Logoff", (char*)0);
      if (logoff && !NoLogOpt()) {
         char *mac = gSystem->Which(TROOT::GetMacroPath(), logoff, kReadPermission);
         if (mac)
            ProcessFile(logoff);
         delete [] mac;
      }

      TApplication::Terminate(status);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set console mode:
///
///  mode = kTRUE  - echo input symbols
///  mode = kFALSE - noecho input symbols

void TRint::SetEchoMode(Bool_t mode)
{
   Gl_config("noecho", mode ? 0 : 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Process the content of a line starting with ".R" (already stripped-off)
/// The format is
///      [user@]host[:dir] [-l user] [-d dbg] [script]
/// The variable 'dir' is the remote directory to be used as working dir.
/// The username can be specified in two ways, "-l" having the priority
/// (as in ssh).
/// A 'dbg' value > 0 gives increasing verbosity.
/// The last argument 'script' allows to specify an alternative script to
/// be executed remotely to startup the session.

Long_t TRint::ProcessRemote(const char *line, Int_t *)
{
   Long_t ret = TApplication::ProcessRemote(line);

   if (ret == 1) {
      if (fAppRemote) {
         TString prompt; prompt.Form("%s:root [%%d] ", fAppRemote->ApplicationName());
         SetPrompt(prompt);
      } else {
         SetPrompt("root [%d] ");
      }
   }

   return ret;
}


////////////////////////////////////////////////////////////////////////////////
/// Calls ProcessLine() possibly prepending a #line directive for
/// better diagnostics. Must be called after fNcmd has been increased for
/// the next line.

Long_t  TRint::ProcessLineNr(const char* filestem, const char *line, Int_t *error /*= 0*/)
{
   Int_t err;
   if (!error)
      error = &err;
   if (line && line[0] != '.') {
      TString lineWithNr = TString::Format("#line 1 \"%s%d\"\n", filestem, fNcmd - 1);
      int res = ProcessLine(lineWithNr + line, kFALSE, error);
      if (*error == TInterpreter::kProcessing) {
         if (!fNonContinuePrompt.Length())
            fNonContinuePrompt = fDefaultPrompt;
         SetPrompt("root (cont'ed, cancel with .@) [%d]");
      } else if (fNonContinuePrompt.Length()) {
         SetPrompt(fNonContinuePrompt);
         fNonContinuePrompt.Clear();
      }
      return res;
   }
   if (line && line[0] == '.' && line[1] == '@') {
      ProcessLine(line, kFALSE, error);
      SetPrompt("root [%d] ");
   }
   return ProcessLine(line, kFALSE, error);
}


////////////////////////////////////////////////////////////////////////////////
/// Forward tab completion request to our TTabCom::Hook().

Int_t TRint::TabCompletionHook(char *buf, int *pLoc, std::ostream& out)
{
   if (gTabCom)
      return gTabCom->Hook(buf, pLoc, out);

   return -1;
}
