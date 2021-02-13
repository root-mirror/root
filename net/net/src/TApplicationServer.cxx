// @(#)root/net:$Id$
// Author: G. Ganis  10/5/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TApplicationServer                                                   //
//                                                                      //
// TApplicationServer is the remote application run by the roots main   //
// program. The input is taken from the socket connection to the client.//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RConfigure.h"
#include <ROOT/RConfig.hxx>
#include "snprintf.h"
#include <iostream>

#ifdef WIN32
   #include <io.h>
   typedef long off_t;
#endif
#include <cstdlib>
#include <cerrno>
#include <ctime>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#if (defined(__FreeBSD__) && (__FreeBSD__ < 4)) || \
    (defined(__APPLE__) && (!defined(MAC_OS_X_VERSION_10_3) || \
     (MAC_OS_X_VERSION_MAX_ALLOWED < MAC_OS_X_VERSION_10_3)))
#include <sys/file.h>
#define lockf(fd, op, sz)   flock((fd), (op))
#ifndef F_LOCK
#define F_LOCK             (LOCK_EX | LOCK_NB)
#endif
#ifndef F_ULOCK
#define F_ULOCK             LOCK_UN
#endif
#endif

#include "RRemoteProtocol.h"

#include "TApplicationServer.h"
#include "TBenchmark.h"
#include "TEnv.h"
#include "TError.h"
#include "TException.h"
#include "TInterpreter.h"
#include "TMD5.h"
#include "TMessage.h"
#include "TROOT.h"
#include "TSocket.h"
#include "TSystem.h"
#include "TRemoteObject.h"
#include "TUrl.h"
#include "compiledata.h"
#include "TClass.h"


//----- Interrupt signal handler -----------------------------------------------
////////////////////////////////////////////////////////////////////////////////

class TASInterruptHandler : public TSignalHandler {
   TApplicationServer  *fServ;
public:
   TASInterruptHandler(TApplicationServer *s)
      : TSignalHandler(kSigUrgent, kFALSE) { fServ = s; }
   Bool_t  Notify();
};

////////////////////////////////////////////////////////////////////////////////
/// Handle this interrupt

Bool_t TASInterruptHandler::Notify()
{
   fServ->HandleUrgentData();
   if (TROOT::Initialized()) {
      Throw(GetSignal());
   }
   return kTRUE;
}

//----- SigPipe signal handler -------------------------------------------------
////////////////////////////////////////////////////////////////////////////////

class TASSigPipeHandler : public TSignalHandler {
   TApplicationServer  *fServ;
public:
   TASSigPipeHandler(TApplicationServer *s) : TSignalHandler(kSigPipe, kFALSE)
      { fServ = s; }
   Bool_t  Notify();
};

////////////////////////////////////////////////////////////////////////////////
/// Handle this signal

Bool_t TASSigPipeHandler::Notify()
{
   fServ->HandleSigPipe();
   return kTRUE;
}

//----- Input handler for messages from client -----------------------
////////////////////////////////////////////////////////////////////////////////

class TASInputHandler : public TFileHandler {
   TApplicationServer  *fServ;
public:
   TASInputHandler(TApplicationServer *s, Int_t fd) : TFileHandler(fd, 1)
      { fServ = s; }
   Bool_t Notify();
   Bool_t ReadNotify() { return Notify(); }
};

////////////////////////////////////////////////////////////////////////////////
/// Handle this input

Bool_t TASInputHandler::Notify()
{
   fServ->HandleSocketInput();
   return kTRUE;
}

TString TASLogHandler::fgPfx = ""; // Default prefix to be prepended to messages
////////////////////////////////////////////////////////////////////////////////
/// Execute 'cmd' in a pipe and handle output messages from the related file

TASLogHandler::TASLogHandler(const char *cmd, TSocket *s, const char *pfx)
              : TFileHandler(-1, 1), fSocket(s), fPfx(pfx)
{
   ResetBit(kFileIsPipe);
   fFile = 0;
   if (s && cmd) {
      fFile = gSystem->OpenPipe(cmd, "r");
      if (fFile) {
         SetFd(fileno(fFile));
         // Notify what already in the file
         Notify();
         // Used in the destructor
         SetBit(kFileIsPipe);
      } else {
         fSocket = 0;
         Error("TASLogHandler", "executing command in pipe");
      }
   } else {
      Error("TASLogHandler",
            "undefined command (%p) or socket (%p)", (int *)cmd, s);
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Handle available message from the open file 'f'

TASLogHandler::TASLogHandler(FILE *f, TSocket *s, const char *pfx)
              : TFileHandler(-1, 1), fSocket(s), fPfx(pfx)
{
   ResetBit(kFileIsPipe);
   fFile = 0;
   if (s && f) {
      fFile = f;
      SetFd(fileno(fFile));
      // Notify what already in the file
      Notify();
   } else {
      Error("TASLogHandler", "undefined file (%p) or socket (%p)", f, s);
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Handle available message in the open file

TASLogHandler::~TASLogHandler()
{
   if (TestBit(kFileIsPipe) && fFile)
      gSystem->ClosePipe(fFile);
   fFile = 0;
   fSocket = 0;
   ResetBit(kFileIsPipe);
}
////////////////////////////////////////////////////////////////////////////////
/// Handle available message in the open file

Bool_t TASLogHandler::Notify()
{
   if (IsValid()) {
      TMessage m(kMESS_ANY);
      // Read buffer
      char line[4096];
      char *plf = 0;
      while (fgets(line, sizeof(line), fFile)) {
         if ((plf = strchr(line, '\n')))
            *plf = 0;
         // Send the message one level up
         m.Reset(kMESS_ANY);
         m << (Int_t)kRRT_Message;
         if (fPfx.Length() > 0) {
            // Prepend prefix specific to this instance
            m << TString(Form("%s: %s", fPfx.Data(), line));
         } else if (fgPfx.Length() > 0) {
            // Prepend default prefix
            m << TString(Form("%s: %s", fgPfx.Data(), line));
         } else {
            // Nothing to prepend
            m << TString(line);
         }
         fSocket->Send(m);
      }
   }
   return kTRUE;
}
////////////////////////////////////////////////////////////////////////////////
/// Static method to set the default prefix

void TASLogHandler::SetDefaultPrefix(const char *pfx)
{
   fgPfx = pfx;
}

////////////////////////////////////////////////////////////////////////////////
/// Init a guard for executing a command in a pipe

TASLogHandlerGuard::TASLogHandlerGuard(const char *cmd, TSocket *s,
                                       const char *pfx, Bool_t on)
{
   fExecHandler = 0;
   if (cmd && on) {
      fExecHandler = new TASLogHandler(cmd, s, pfx);
      if (fExecHandler->IsValid()) {
         gSystem->AddFileHandler(fExecHandler);
      } else {
         Error("TASLogHandlerGuard","invalid handler");
      }
   } else {
      if (on)
         Error("TASLogHandlerGuard","undefined command");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Init a guard for executing a command in a pipe

TASLogHandlerGuard::TASLogHandlerGuard(FILE *f, TSocket *s,
                                       const char *pfx, Bool_t on)
{
   fExecHandler = 0;
   if (f && on) {
      fExecHandler = new TASLogHandler(f, s, pfx);
      if (fExecHandler->IsValid()) {
         gSystem->AddFileHandler(fExecHandler);
      } else {
         Error("TASLogHandlerGuard","invalid handler");
      }
   } else {
      if (on)
         Error("TASLogHandlerGuard","undefined file");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Close a guard for executing a command in a pipe

TASLogHandlerGuard::~TASLogHandlerGuard()
{
   if (fExecHandler && fExecHandler->IsValid()) {
      gSystem->RemoveFileHandler(fExecHandler);
      SafeDelete(fExecHandler);
   }
}

ClassImp(TApplicationServer);

////////////////////////////////////////////////////////////////////////////////
/// Main constructor. Create an application environment. The TApplicationServer
/// environment provides an eventloop via inheritance of TApplication.

TApplicationServer::TApplicationServer(Int_t *argc, char **argv,
                                       FILE *flog, const char *logfile)
       : TApplication("server", argc, argv, 0, -1)
{
   // Parse options
   GetOptions(argc, argv);

   // Abort on higher than kSysError's and set error handler
   gErrorAbortLevel = kSysError + 1;
   SetErrorHandler(ErrorHandler);

   fInterrupt       = kFALSE;
   fSocket          = 0;
   fWorkingDir      = 0;

   fLogFilePath     = logfile;
   fLogFile         = flog;
   fLogFileDes      = -1;
   if (!fLogFile || (fLogFileDes = fileno(fLogFile)) < 0)
      // For some reason we failed setting a redirection; we cannot continue
      Terminate(0);
   fRealTimeLog     = kFALSE;
   fSentCanvases    = 0;

   // Default prefix for notifications
   TASLogHandler::SetDefaultPrefix(Form("roots:%s", gSystem->HostName()));

   // Now we contact back the client: if we fail we set ourselves
   // as invalid
   fIsValid = kFALSE;

   if (!(fSocket = new TSocket(GetHost(), GetPort()))) {
      Terminate(0);
      return;
   }
   Int_t sock = fSocket->GetDescriptor();

   if (Setup() != 0) {
      Error("TApplicationServer", "failed to setup - quitting");
      SendLogFile(-98);
      Terminate(0);
   }

   // Everybody expects std::iostream to be available, so load it...
   ProcessLine("#include <iostream>", kTRUE);
   ProcessLine("#include <string>",kTRUE); // for std::string std::iostream.

   // Load user functions
   const char *logon;
   logon = gEnv->GetValue("Rint.Load", (char *)0);
   if (logon) {
      char *mac = gSystem->Which(TROOT::GetMacroPath(), logon, kReadPermission);
      if (mac)
         ProcessLine(Form(".L %s", logon), kTRUE);
      delete [] mac;
   }

   // Execute logon macro
   ExecLogon();

   // Init benchmarking
   gBenchmark = new TBenchmark();

   // Save current interpreter context
   gInterpreter->SaveContext();
   gInterpreter->SaveGlobalsContext();

   // Install interrupt and message input handlers
   gSystem->AddSignalHandler(new TASInterruptHandler(this));
   gSystem->AddFileHandler(new TASInputHandler(this, sock));

   // We are done
   fIsValid = kTRUE;

   // Startup notification
   BrowseDirectory(0);
   SendLogFile();
}

////////////////////////////////////////////////////////////////////////////////
/// Print the Remote Server logo on standard output.
/// Return 0 on success, -1 on failure

Int_t TApplicationServer::Setup()
{
   char str[512];
   snprintf(str, 512, "**** Remote session @ %s started ****", gSystem->HostName());
   if (fSocket->Send(str) != 1+static_cast<Int_t>(strlen(str))) {
      Error("Setup", "failed to send startup message");
      return -1;
   }

   // Send our protocol level to the client
   if (fSocket->Send(kRRemote_Protocol, kROOTD_PROTOCOL) != 2*sizeof(Int_t)) {
      Error("Setup", "failed to send local protocol");
      return -1;
   }

   // Send the host name and full path to log file
   TMessage msg(kMESS_ANY);
   msg << TString(gSystem->HostName()) << fLogFilePath;
   fSocket->Send(msg);

   // Set working directory
   fWorkDir = gSystem->WorkingDirectory();
   if (strlen(fUrl.GetFile()) > 0) {
      fWorkDir = fUrl.GetFile();
      gSystem->ExpandPathName(fWorkDir);
   }

   // Go to working dir
   if (gSystem->AccessPathName(fWorkDir)) {
      gSystem->mkdir(fWorkDir, kTRUE);
      if (!gSystem->ChangeDirectory(fWorkDir)) {
         SysError("Setup", "can not change to directory %s",
                  fWorkDir.Data());
      }
   } else {
      if (!gSystem->ChangeDirectory(fWorkDir)) {
         gSystem->Unlink(fWorkDir);
         gSystem->mkdir(fWorkDir, kTRUE);
         if (!gSystem->ChangeDirectory(fWorkDir)) {
            SysError("Setup", "can not change to directory %s",
                     fWorkDir.Data());
         }
      }
   }

#if 0 // G.Ganis May 11, 2007
   // This needs to be fixed: we disable for the time being
   // Socket options: incoming OOB should generate a SIGURG
   if (fSocket->SetOption(kProcessGroup, (-1)*gSystem->GetPid()) != 0)
      SysWarning("Setup", "failed to enable SIGURG generation on incoming OOB");
#endif

   // Send messages off immediately to reduce latency
   if (fSocket->SetOption(kNoDelay, 1) != 0) {}
      //SysWarning("Setup", "failed to set no-delay option on input socket");

   // Check every two hours if client is still alive
   if (fSocket->SetOption(kKeepAlive, 1) != 0) {}
      //SysWarning("Setup", "failed to set keepalive option on input socket");

   // Install SigPipe handler to handle kKeepAlive failure
   gSystem->AddSignalHandler(new TASSigPipeHandler(this));

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup. Not really necessary since after this dtor there is no
/// live anyway.

TApplicationServer::~TApplicationServer()
{
   fSentCanvases->SetOwner(kFALSE);
   SafeDelete(fSentCanvases);
   SafeDelete(fSocket);
   close(fLogFileDes);
}

////////////////////////////////////////////////////////////////////////////////
/// Get and handle command line options. Fixed format:
/// "protocol  url"

void TApplicationServer::GetOptions(Int_t *argc, char **argv)
{
   if (*argc < 4) {
      Fatal("GetOptions", "must be started with 4 arguments");
      gSystem->Exit(1);
   }

   // Protocol run by the client
   fProtocol = TString(argv[1]).Atoi();

   // Client URL
   fUrl.SetUrl(argv[2]);

   // Debug level
   gDebug = 0;
   TString argdbg(argv[3]);
   if (argdbg.BeginsWith("-d=")) {
      argdbg.ReplaceAll("-d=","");
      gDebug = argdbg.Atoi();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Main server eventloop.

void TApplicationServer::Run(Bool_t retrn)
{
   // Setup the server
   if (fIsValid) {
      // Run the main event loop
      TApplication::Run(retrn);
   } else {
      Error("Run", "invalid instance: cannot Run()");
      gSystem->Exit(1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle input coming from the client or from the master server.

void TApplicationServer::HandleSocketInput()
{
   TMessage *mess;
   char      str[2048];
   Int_t     what;

   if (fSocket->Recv(mess) <= 0) {
      // Pending: do something more intelligent here
      // but at least get a message in the log file
      Error("HandleSocketInput", "retrieving message from input socket");
      Terminate(0);
      return;
   }

   what = mess->What();
   if (gDebug > 0)
      Info("HandleSocketInput", "got message of type %d", what);

   switch (what) {

      case kMESS_CINT:
         {  TASLogHandlerGuard hg(fLogFile, fSocket, "", fRealTimeLog);
            mess->ReadString(str, sizeof(str));
            if (gDebug > 1)
               Info("HandleSocketInput:kMESS_CINT", "processing: %s...", str);
            ProcessLine(str);
         }
         SendCanvases();
         SendLogFile();
         break;

      case kMESS_STRING:
         mess->ReadString(str, sizeof(str));
         break;

      case kMESS_OBJECT:
         mess->ReadObject(mess->GetClass());
         break;

      case kMESS_ANY:
         {
            Int_t type;
            (*mess) >> type;
            switch (type) {
               case kRRT_Reset:
                  mess->ReadString(str, sizeof(str));
                  Reset(str);
                  break;

               case kRRT_CheckFile:
                  // Handle file checking request
                  HandleCheckFile(mess);
                  break;

               case kRRT_File:
                  // A file follows
                  mess->ReadString(str, sizeof(str));
                  {  char   name[2048], i1[20], i2[40];
                     sscanf(str, "%2047s %19s %39s", name, i1, i2);
                     Int_t  bin = atoi(i1);
                     Long_t size = atol(i2);
                     ReceiveFile(name, bin ? kTRUE : kFALSE, size);
                  }
                  break;

               case kRRT_Terminate:
                  // Terminate the session (will not return from here)
                  Int_t status;
                  (*mess) >> status;
                  Terminate(status);
                  break;

               default:
                  break;
            }
         }
         SendLogFile();
         break;
      default:
         Warning("HandleSocketInput","message type unknown (%d)", what);
         SendLogFile();
         break;
   }

   delete mess;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle Out-Of-Band data sent by the master or client.

void TApplicationServer::HandleUrgentData()
{
   char  oob_byte;
   Int_t n, nch, wasted = 0;

   const Int_t kBufSize = 1024;
   char waste[kBufSize];

   // Real-time notification of messages
   TASLogHandlerGuard hg(fLogFile, fSocket, "", fRealTimeLog);

   Info("HandleUrgentData", "handling oob...");

   // Receive the OOB byte
   while ((n = fSocket->RecvRaw(&oob_byte, 1, kOob)) < 0) {
      if (n == -2) {   // EWOULDBLOCK
         //
         // The OOB data has not yet arrived: flush the input stream
         //
         // In some systems (Solaris) regular recv() does not return upon
         // receipt of the oob byte, which makes the below call to recv()
         // block indefinitely if there are no other data in the queue.
         // FIONREAD ioctl can be used to check if there are actually any
         // data to be flushed.  If not, wait for a while for the oob byte
         // to arrive and try to read it again.
         //
         fSocket->GetOption(kBytesToRead, nch);
         if (nch == 0) {
            gSystem->Sleep(1000);
            continue;
         }

         if (nch > kBufSize) nch = kBufSize;
         n = fSocket->RecvRaw(waste, nch);
         if (n <= 0) {
            Error("HandleUrgentData", "error receiving waste");
            break;
         }
         wasted = 1;
      } else {
         Error("HandleUrgentData", "error receiving OOB (n = %d)",n);
         return;
      }
   }

   Info("HandleUrgentData", "got OOB byte: %d\n", oob_byte);

   switch (oob_byte) {

      case kRRI_Hard:
         Info("HandleUrgentData", "*** Hard Interrupt");

         // Flush input socket
         while (1) {
            Int_t atmark;

            fSocket->GetOption(kAtMark, atmark);

            if (atmark) {
               // Send the OOB byte back so that the client knows where
               // to stop flushing its input stream of obsolete messages
               n = fSocket->SendRaw(&oob_byte, 1, kOob);
               if (n <= 0)
                  Error("HandleUrgentData", "error sending OOB");
               break;
            }

            // find out number of bytes to read before atmark
            fSocket->GetOption(kBytesToRead, nch);
            if (nch == 0) {
               gSystem->Sleep(1000);
               continue;
            }

            if (nch > kBufSize) nch = kBufSize;
            n = fSocket->RecvRaw(waste, nch);
            if (n <= 0) {
               Error("HandleUrgentData", "error receiving waste (2)");
               break;
            }
         }

         SendLogFile();

         break;

      case kRRI_Soft:
         Info("HandleUrgentData", "Soft Interrupt");

         if (wasted) {
            Error("HandleUrgentData", "soft interrupt flushed stream");
            break;
         }

         Interrupt();

         SendLogFile();

         break;

      case kRRI_Shutdown:
         Info("HandleUrgentData", "Shutdown Interrupt");

         Terminate(0);

         break;

      default:
         Error("HandleUrgentData", "unexpected OOB byte");
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called when the client is not alive anymore (i.e. when kKeepAlive
/// has failed).

void TApplicationServer::HandleSigPipe()
{
   // Real-time notification of messages
   TASLogHandlerGuard hg(fLogFile, fSocket, "", fRealTimeLog);

   Info("HandleSigPipe", "client died");
   Terminate(0);  // will not return from here....
}

////////////////////////////////////////////////////////////////////////////////
/// Reset environment to be ready for execution of next command.

void TApplicationServer::Reset(const char *dir)
{
   // First go to new directory.
   gDirectory->cd(dir);

   // Clear interpreter environment.
   gROOT->Reset();

   // Make sure current directory is empty (don't delete anything when
   // we happen to be in the ROOT memory only directory!?)
   if (gDirectory != gROOT) {
      gDirectory->Delete();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Receive a file, either sent by a client or a master server.
/// If bin is true it is a binary file, other wise it is an ASCII
/// file and we need to check for Windows \r tokens. Returns -1 in
/// case of error, 0 otherwise.

Int_t TApplicationServer::ReceiveFile(const char *file, Bool_t bin, Long64_t size)
{
   if (size <= 0) return 0;

   // open file, overwrite already existing file
   Int_t fd = open(file, O_CREAT | O_TRUNC | O_WRONLY, 0600);
   if (fd < 0) {
      SysError("ReceiveFile", "error opening file %s", file);
      return -1;
   }

   const Int_t kMAXBUF = 16384;  //32768  //16384  //65536;
   char buf[kMAXBUF], cpy[kMAXBUF];

   Int_t    left, r;
   Long64_t filesize = 0;

   while (filesize < size) {
      left = Int_t(size - filesize);
      if (left > kMAXBUF)
         left = kMAXBUF;
      r = fSocket->RecvRaw(&buf, left);
      if (r > 0) {
         char *p = buf;

         filesize += r;
         while (r) {
            Int_t w;

            if (!bin) {
               Int_t k = 0, i = 0, j = 0;
               char *q;
               while (i < r) {
                  if (p[i] == '\r') {
                     i++;
                     k++;
                  }
                  cpy[j++] = buf[i++];
               }
               q = cpy;
               r -= k;
               w = write(fd, q, r);
            } else {
               w = write(fd, p, r);
            }

            if (w < 0) {
               SysError("ReceiveFile", "error writing to file %s", file);
               close(fd);
               return -1;
            }
            r -= w;
            p += w;
         }
      } else if (r < 0) {
         Error("ReceiveFile", "error during receiving file %s", file);
         close(fd);
         return -1;
      }
   }

   close(fd);

   chmod(file, 0644);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Send log file to master.
/// If start > -1 send only bytes in the range from start to end,
/// if end <= start send everything from start.

void TApplicationServer::SendLogFile(Int_t status, Int_t start, Int_t end)
{
   // Determine the number of bytes left to be read from the log file.
   fflush(stdout);

   off_t ltot=0, lnow=0;
   Int_t left = -1;
   Bool_t adhoc = kFALSE;

   if (fLogFileDes > -1) {
      ltot = lseek(fileno(stdout), (off_t) 0, SEEK_END);
      lnow = lseek(fLogFileDes, (off_t) 0, SEEK_CUR);
      if (lnow == -1) {
         SysError("SendLogFile", "lseek failed");
         lnow = 0;
      }

      if (start > -1) {
         lseek(fLogFileDes, (off_t) start, SEEK_SET);
         if (end <= start || end > ltot)
            end = ltot;
         left = (Int_t)(end - start);
         if (end < ltot)
            left++;
         adhoc = kTRUE;
      } else {
         left = (Int_t)(ltot - lnow);
      }
   }

   TMessage m(kMESS_ANY);

   if (left > 0) {

      m << (Int_t)kRRT_LogFile << left;
      fSocket->Send(m);

      const Int_t kMAXBUF = 32768;  //16384  //65536;
      char buf[kMAXBUF];
      Int_t wanted = (left > kMAXBUF) ? kMAXBUF : left;
      Int_t len;
      do {
         while ((len = read(fLogFileDes, buf, wanted)) < 0 &&
                TSystem::GetErrno() == EINTR)
            TSystem::ResetErrno();

         if (len < 0) {
            SysError("SendLogFile", "error reading log file");
            break;
         }

         if (end == ltot && len == wanted)
            buf[len-1] = '\n';

         if (fSocket->SendRaw(buf, len) < 0) {
            SysError("SendLogFile", "error sending log file");
            break;
         }

         // Update counters
         left -= len;
         wanted = (left > kMAXBUF) ? kMAXBUF : left;

      } while (len > 0 && left > 0);
   }

   // Restore initial position if partial send
   if (adhoc)
      lseek(fLogFileDes, lnow, SEEK_SET);

   m.Reset();
   m << (Int_t)kRRT_LogDone << status;

   fSocket->Send(m);
}

////////////////////////////////////////////////////////////////////////////////
/// Send any created canvas to client

Int_t TApplicationServer::SendCanvases()
{
   Int_t nc = 0;

   // Send back new canvases
   TMessage mess(kMESS_OBJECT);
   TIter next(gROOT->GetListOfCanvases());
   TObject *o = 0;
   while ((o = next())) {
      if (!fSentCanvases)
         fSentCanvases = new TList;
      Bool_t sentalready = kFALSE;
      // We cannot use FindObject here because there may be invalid
      // objects in the send list (i.e. deleted canvases)
      TObjLink *lnk = fSentCanvases->FirstLink();
      while (lnk) {
         TObject *sc = lnk->GetObject();
         lnk = lnk->Next();
         if ((sc->TestBit(kNotDeleted)) && sc == o)
            sentalready = kTRUE;
      }
      if (!sentalready) {
         if (gDebug > 0)
            Info("SendCanvases","new canvas found: %p", o);
         mess.Reset(kMESS_OBJECT);
         mess.WriteObject(o);
         fSocket->Send(mess);
         nc++;
         fSentCanvases->Add(o);
      }
   }
   return nc;
}

////////////////////////////////////////////////////////////////////////////////
/// Browse directory and send back its content to client.

Int_t TApplicationServer::BrowseDirectory(const char *dirname)
{
   Int_t nc = 0;

   TMessage mess(kMESS_OBJECT);
   if (!fWorkingDir || !dirname || !*dirname) {
      if (!fWorkingDir)
         fWorkingDir = new TRemoteObject(fWorkDir, fWorkDir, "TSystemDirectory");
      TList *list = fWorkingDir->Browse();
      mess.Reset(kMESS_OBJECT);
      mess.WriteObject(fWorkingDir);
      fSocket->Send(mess);
      nc++;
      list->Delete();
      delete list;
   }
   else if (fWorkingDir) {
      TRemoteObject dir(dirname, dirname, "TSystemDirectory");
      TList *list = dir.Browse();
      mess.Reset(kMESS_OBJECT);
      mess.WriteObject(list);
      fSocket->Send(mess);
      nc++;
      list->Delete();
      delete list;
   }
   return nc;
}

////////////////////////////////////////////////////////////////////////////////
/// Browse root file and send back its content;
/// if fname is null, send the full list of files.

Int_t TApplicationServer::BrowseFile(const char *fname)
{
   Int_t nc = 0;

   TList *list = new TList;
   TMessage mess(kMESS_OBJECT);
   if (!fname || !*fname) {
      // fname is null, so send the list of files.
      TIter next(gROOT->GetListOfFiles());
      TNamed *fh = 0;
      TRemoteObject *robj;
      while ((fh = (TNamed *)next())) {
         robj = new TRemoteObject(fh->GetName(), fh->GetTitle(), "TFile");
         list->Add(robj);
      }
      if (list->GetEntries() > 0) {
         mess.Reset(kMESS_OBJECT);
         mess.WriteObject(list);
         fSocket->Send(mess);
         nc++;
      }
   }
   else {
      // get Root file content and send the list of objects
      TDirectory *fh = (TDirectory *)gROOT->GetListOfFiles()->FindObject(fname);
      if (fh) {
         fh->cd();
         TRemoteObject dir(fh->GetName(), fh->GetTitle(), "TFile");
         TList *keylist = (TList *)gROOT->ProcessLine(Form("((TFile *)0x%lx)->GetListOfKeys();", (ULong_t)fh));
         TIter nextk(keylist);
         TNamed *key = 0;
         TRemoteObject *robj;
         while ((key = (TNamed *)nextk())) {
            robj = new TRemoteObject(key->GetName(), key->GetTitle(), "TKey");
            const char *classname = (const char *)gROOT->ProcessLine(Form("((TKey *)0x%lx)->GetClassName();", (ULong_t)key));
            robj->SetKeyClassName(classname);
            Bool_t isFolder = (Bool_t)gROOT->ProcessLine(Form("((TKey *)0x%lx)->IsFolder();", (ULong_t)key));
            robj->SetFolder(isFolder);
            robj->SetRemoteAddress((Long_t) key);
            list->Add(robj);
         }
         if (list->GetEntries() > 0) {
            mess.Reset(kMESS_OBJECT);
            mess.WriteObject(list);
            fSocket->Send(mess);
            nc++;
         }
      }
   }

   list->Delete();
   delete list;

   return nc;
}

////////////////////////////////////////////////////////////////////////////////
/// Read key object and send it back to client.

Int_t TApplicationServer::BrowseKey(const char *keyname)
{
   Int_t nc = 0;

   TMessage mess(kMESS_OBJECT);
   TNamed *obj = (TNamed *)gROOT->ProcessLine(Form("gFile->GetKey(\"%s\")->ReadObj();", keyname));
   if (obj) {
      mess.Reset(kMESS_OBJECT);
      mess.WriteObject(obj);
      fSocket->Send(mess);
      nc++;
   }
   return nc;
}

////////////////////////////////////////////////////////////////////////////////
/// Terminate the proof server.

void TApplicationServer::Terminate(Int_t status)
{
   // Close and remove the log file; remove the cleanup script
   if (fLogFile) {
      fclose(fLogFile);
      // Delete the log file unless we are in debug mode
      if (gDebug <= 0)
         gSystem->Unlink(fLogFilePath);
      TString cleanup = fLogFilePath;
      cleanup.ReplaceAll(".log", ".cleanup");
      gSystem->Unlink(cleanup);
   }

   // Remove input handler to avoid spurious signals in socket
   // selection for closing activities executed upon exit()
   TIter next(gSystem->GetListOfFileHandlers());
   TObject *fh = 0;
   while ((fh = next())) {
      TASInputHandler *ih = dynamic_cast<TASInputHandler *>(fh);
      if (ih)
         gSystem->RemoveFileHandler(ih);
   }

   // Stop processing events
   gSystem->Exit(status);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle file checking request.

void TApplicationServer::HandleCheckFile(TMessage *mess)
{
   TString filenam;
   TMD5    md5;
   TMessage m(kMESS_ANY);

   // Parse message
   (*mess) >> filenam >> md5;

   // check file in working directory
   TMD5 *md5local = TMD5::FileChecksum(filenam);
   if (md5local && md5 == (*md5local)) {
      // We have an updated copy of the file
      m << (Int_t) kRRT_CheckFile << (Bool_t) kTRUE;
      fSocket->Send(m);
      if (gDebug > 0)
         Info("HandleCheckFile", "up-to-date version of %s available", filenam.Data());
   } else {
      m << (Int_t) kRRT_CheckFile << (Bool_t) kFALSE;
      fSocket->Send(m);
      if (gDebug > 0)
         Info("HandleCheckFile", "file %s needs to be uploaded", filenam.Data());
   }
   delete md5local;
}

////////////////////////////////////////////////////////////////////////////////
/// The error handler function. It prints the message on stderr and
/// if abort is set it aborts the application.

void TApplicationServer::ErrorHandler(Int_t level, Bool_t abort, const char *location,
                              const char *msg)
{
   if (gErrorIgnoreLevel == kUnset) {
      gErrorIgnoreLevel = 0;
      if (gEnv) {
         TString slevel = gEnv->GetValue("Root.ErrorIgnoreLevel", "Print");
         if (!slevel.CompareTo("Print", TString::kIgnoreCase))
            gErrorIgnoreLevel = kPrint;
         else if (!slevel.CompareTo("Info", TString::kIgnoreCase))
            gErrorIgnoreLevel = kInfo;
         else if (!slevel.CompareTo("Warning", TString::kIgnoreCase))
            gErrorIgnoreLevel = kWarning;
         else if (!slevel.CompareTo("Error", TString::kIgnoreCase))
            gErrorIgnoreLevel = kError;
         else if (!slevel.CompareTo("Break", TString::kIgnoreCase))
            gErrorIgnoreLevel = kBreak;
         else if (!slevel.CompareTo("SysError", TString::kIgnoreCase))
            gErrorIgnoreLevel = kSysError;
         else if (!slevel.CompareTo("Fatal", TString::kIgnoreCase))
            gErrorIgnoreLevel = kFatal;
      }
   }

   if (level < gErrorIgnoreLevel)
      return;

   static TString syslogService;

   if (syslogService.IsNull()) {
      syslogService = "server";
      gSystem->Openlog(syslogService, kLogPid | kLogCons, kLogLocal5);
   }

   const char *type   = 0;
   ELogLevel loglevel = kLogInfo;

   if (level >= kPrint) {
      loglevel = kLogInfo;
      type = "Print";
   }
   if (level >= kInfo) {
      loglevel = kLogInfo;
      type = "Info";
   }
   if (level >= kWarning) {
      loglevel = kLogWarning;
      type = "Warning";
   }
   if (level >= kError) {
      loglevel = kLogErr;
      type = "Error";
   }
   if (level >= kBreak) {
      loglevel = kLogErr;
      type = "*** Break ***";
   }
   if (level >= kSysError) {
      loglevel = kLogErr;
      type = "SysError";
   }
   if (level >= kFatal) {
      loglevel = kLogErr;
      type = "Fatal";
   }

   TString node = "server";
   TString buf;

   if (!location || !location[0] ||
       (level >= kPrint && level < kInfo) ||
       (level >= kBreak && level < kSysError)) {
      fprintf(stderr, "%s on %s: %s\n", type, node.Data(), msg);
      buf.Form("%s:%s:%s", node.Data(), type, msg);
   } else {
      fprintf(stderr, "%s in <%s> on %s: %s\n", type, location, node.Data(), msg);
      buf.Form("%s:%s:<%s>:%s", node.Data(), type, location, msg);
   }
   fflush(stderr);

   gSystem->Syslog(loglevel, buf);

   if (abort) {
      fprintf(stderr, "aborting\n");
      fflush(stderr);
      gSystem->StackTrace();
      gSystem->Abort();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Parse a command line received from the client, making sure that the files
/// needed for the execution, if any, are available. The line is either a C++
/// statement or an interpreter command starting with a ".".
/// Return the return value of the command casted to a long.

Long_t TApplicationServer::ProcessLine(const char *line, Bool_t, Int_t *)
{
   if (!line || !*line) return 0;

   // If load or execute request we must make sure that we have the files.
   // If not we ask the client to send them, blocking until we have everything.
   if (!strncmp(line, ".L", 2) || !strncmp(line, ".U", 2) ||
       !strncmp(line, ".X", 2) || !strncmp(line, ".x", 2)) {
      TString aclicMode;
      TString arguments;
      TString io;
      TString fname = gSystem->SplitAclicMode(line+3, aclicMode, arguments, io);

      char *imp = gSystem->Which(TROOT::GetMacroPath(), fname, kReadPermission);
      if (!imp) {

         // Make sure that we can write in the directory where we are
         if (gSystem->AccessPathName(gSystem->WorkingDirectory(), kWritePermission)) {
            Error("ProcessLine","no write permission in %s", gSystem->WorkingDirectory());
            return 0;
         }

         if (gDebug > 0)
            Info("ProcessLine", "macro %s not found in path %s: asking the client",
                                fname.Data(), TROOT::GetMacroPath());
         TMessage m(kMESS_ANY);
         m << (Int_t) kRRT_SendFile << TString(gSystem->BaseName(fname));
         fSocket->Send(m);

         // Wait for the reply(ies)
         Int_t type;
         Bool_t filefollows = kTRUE;

         while (filefollows) {

            // Get a message
            TMessage *rm = 0;
            if (fSocket->Recv(rm) <= 0) {
               Error("ProcessLine","ask-file: received empty message from client");
               return 0;
            }
            if (rm->What() != kMESS_ANY) {
               Error("ProcessLine","ask-file: wrong message received (what: %d)", rm->What());
               return 0;
            }
            (*rm) >> type;
            if (type != kRRT_SendFile) {
               Error("ProcessLine","ask-file: wrong sub-type received (type: %d)", type);
               return 0;
            }
            (*rm) >> filefollows;
            if (filefollows) {
               // Read the file specifications
               if (fSocket->Recv(rm) <= 0) {
                  Error("ProcessLine","file: received empty message from client");
                  return 0;
               }
               if (rm->What() != kMESS_ANY) {
                  Error("ProcessLine","file: wrong message received (what: %d)", rm->What());
                  return 0;
               }
               (*rm) >> type;
               if (type != kRRT_File) {
                  Error("ProcessLine","file: wrong sub-type received (type: %d)", type);
                  return 0;
               }
               // A file follows
               char str[2048];
               rm->ReadString(str, sizeof(str));
               char   name[2048], i1[20], i2[40];
               sscanf(str, "%2047s %19s %39s", name, i1, i2);
               Int_t  bin = atoi(i1);
               Long_t size = atol(i2);
               ReceiveFile(name, bin ? kTRUE : kFALSE, size);
            }
         }
      }
      delete [] imp;
   }

   // Process the line now
   return TApplication::ProcessLine(line);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute logon macro's. There are three levels of logon macros that
/// will be executed: the system logon etc/system.rootlogon.C, the global
/// user logon ~/.rootlogon.C and the local ./.rootlogon.C. For backward
/// compatibility also the logon macro as specified by the Rint.Logon
/// environment setting, by default ./rootlogon.C, will be executed.
/// No logon macros will be executed when the system is started with
/// the -n option.

void TApplicationServer::ExecLogon()
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
