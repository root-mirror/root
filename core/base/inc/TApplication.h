// @(#)root/base:$Id$
// Author: Fons Rademakers   22/12/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TApplication
#define ROOT_TApplication


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TApplication                                                         //
//                                                                      //
// This class creates the ROOT Application Environment that interfaces  //
// to the windowing system eventloop and eventhandlers.                 //
// This class must be instantiated exactly once in any given            //
// application. Normally the specific application class inherits from   //
// TApplication (see TRint).                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"

#include "TQObject.h"

#include "TApplicationImp.h"

class TObjArray;
class TTimer;
class TSignalHandler;


class TApplication : public TObject, public TQObject {

public:
   // TApplication specific bits
   enum EStatusBits {
      kProcessRemotely    = BIT(15), // TRUE if this line has to be processed remotely
      kDefaultApplication = BIT(16)  // TRUE if created via CreateApplication()
   };
   // TApplication specific bits for fFiles
   enum EFileBits {
      kExpression = BIT(14)  // If the arguments is an expression rather than a file.
   };
   enum EExitOnException {
      kDontExit,
      kExit,
      kAbort
   };

private:
   Int_t              fArgc;            //Number of com   mand line arguments
   char             **fArgv;            //Command line arguments
   TApplicationImp   *fAppImp;          //!Window system specific application implementation
   Bool_t             fIsRunning;       //True when in event loop (Run() has been called)
   Bool_t             fReturnFromRun;   //When true return from Run()
   Bool_t             fNoLog;           //Do not process logon and logoff macros
   Bool_t             fNoLogo;          //Do not show splash screen and welcome message
   Bool_t             fQuit;            //Exit after having processed input files
   Bool_t             fUseMemstat;      //Run with TMemStat enabled
   TObjArray         *fFiles;           //Array of input files or C++ expression (TObjString's) specified via argv
   TString            fWorkDir;         //Working directory specified via argv
   TString            fIdleCommand;     //Command to execute while application is idle
   TTimer            *fIdleTimer;       //Idle timer
   TSignalHandler    *fSigHandler;      //Interrupt handler
   EExitOnException   fExitOnException; //Exit on exception option

   static Bool_t      fgGraphNeeded;    // True if graphics libs need to be initialized
   static Bool_t      fgGraphInit;      // True if graphics libs initialized

   TApplication(const TApplication&) = delete;
   TApplication& operator=(const TApplication&) = delete;

protected:
   TApplication      *fAppRemote;      //Current remote application, if defined

   static TList      *fgApplications;  //List of available applications

   TApplication();

   virtual Long_t     ProcessRemote(const char *line, Int_t *error = 0);
   virtual void       Help(const char *line);
   virtual void       LoadGraphicsLibs();
   virtual void       MakeBatch();
   void               SetSignalHandler(TSignalHandler *sh) { fSigHandler = sh; }

   static Int_t       ParseRemoteLine(const char *ln,
                                      TString &hostdir, TString &user,
                                      Int_t &dbg, TString &script);
   static TApplication *Open(const char *url, Int_t debug, const char *script);
   static void          Close(TApplication *app);

public:
   TApplication(const char *appClassName, Int_t *argc, char **argv,
                void *options = 0, Int_t numOptions = 0);
   virtual ~TApplication();

   void            InitializeGraphics();
   virtual void    GetOptions(Int_t *argc, char **argv);
   TSignalHandler *GetSignalHandler() const { return fSigHandler; }
   virtual void    SetEchoMode(Bool_t mode);
   void OpenInBrowser(const TString & url);
   void OpenReferenceGuideFor(const TString & strippedClass);
   virtual void    HandleException(Int_t sig);
   virtual void    HandleIdleTimer();   //*SIGNAL*
   virtual Bool_t  HandleTermInput() { return kFALSE; }
   virtual void    Init() { fAppImp->Init(); }
   virtual Long_t  ProcessLine(const char *line, Bool_t sync = kFALSE, Int_t *error = 0);
   virtual Long_t  ProcessFile(const char *file, Int_t *error = 0, Bool_t keep = kFALSE);
   virtual void    Run(Bool_t retrn = kFALSE);
   virtual void    SetIdleTimer(UInt_t idleTimeInSec, const char *command);
   virtual void    RemoveIdleTimer();
   const char     *GetIdleCommand() const { return fIdleCommand; }
   virtual void    StartIdleing();
   virtual void    StopIdleing();
   EExitOnException ExitOnException(EExitOnException opt = kExit);

   virtual const char *ApplicationName() const { return fAppImp->ApplicationName(); }
   virtual void    Show()    { fAppImp->Show(); }
   virtual void    Hide()    { fAppImp->Hide(); }
   virtual void    Iconify() { fAppImp->Iconify(); }
   virtual void    Open()    { fAppImp->Open(); }
   virtual void    Raise()   { fAppImp->Raise(); }
   virtual void    Lower()   { fAppImp->Lower(); }
   virtual Bool_t  IsCmdThread() { return fAppImp ? fAppImp->IsCmdThread() : kTRUE; }
   virtual TApplicationImp *GetApplicationImp() { return fAppImp;}

   virtual void    ls(Option_t *option="") const;

   Int_t           Argc() const  { return fArgc; }
   char          **Argv() const  { return fArgv; }
   char           *Argv(Int_t index) const;
   Bool_t          NoLogOpt() const { return fNoLog; }
   Bool_t          NoLogoOpt() const { return fNoLogo; }
   Bool_t          QuitOpt() const { return fQuit; }
   TObjArray      *InputFiles() const { return fFiles; }
   const char     *WorkingDirectory() const { return fWorkDir; }
   void            ClearInputFiles();

   TApplication   *GetAppRemote() const { return fAppRemote; }

   Bool_t          IsRunning() const { return fIsRunning; }
   Bool_t          ReturnFromRun() const { return fReturnFromRun; }
   void            SetReturnFromRun(Bool_t ret) { fReturnFromRun = ret; }

   virtual void    LineProcessed(const char *line);   //*SIGNAL*
   virtual void    Terminate(Int_t status = 0);       //*SIGNAL*
   virtual void    KeyPressed(Int_t key);             //*SIGNAL*
   virtual void    ReturnPressed(char *text );        //*SIGNAL*
   virtual Int_t   TabCompletionHook(char *buf, int *pLoc, std::ostream& out);

   static Long_t   ExecuteFile(const char *file, Int_t *error = 0, Bool_t keep = kFALSE);
   static TList   *GetApplications();
   static void     CreateApplication();
   static void     NeedGraphicsLibs();

   ClassDef(TApplication,0)  //GUI application singleton
};

R__EXTERN TApplication *gApplication;

#endif
