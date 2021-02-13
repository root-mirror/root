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
// TApplicationRemote                                                   //
//                                                                      //
// TApplicationRemote maps a remote session. It starts a remote session //
// and takes care of redirecting the commands to be processed to the    //
// remote session, to collect the graphic output objects and to display //
// them locally.                                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <errno.h>
#include <random>

#include "TApplicationRemote.h"

#include "TBrowser.h"
#include "TDirectory.h"
#include "TError.h"
#include "THashList.h"
#include "TMonitor.h"
#include "TROOT.h"
#include "TServerSocket.h"
#include "TSystem.h"
#include "TRemoteObject.h"
#include "snprintf.h"
#ifdef WIN32
#include <io.h>
#include <sys/types.h>
#endif

//
// TApplicationRemote Interrupt signal handler
////////////////////////////////////////////////////////////////////////////////
/// TApplicationRemote interrupt handler.

Bool_t TARInterruptHandler::Notify()
{
   Info("Notify","Processing interrupt signal ...");

   // Handle interrupt condition on socket(s)
   fApplicationRemote->Interrupt(kRRI_Hard);

   return kTRUE;
}


ClassImp(TApplicationRemote);

static const char *gScript = "roots";
static const char *gScriptCmd = "\\\"%s %d localhost:%d/%s -d=%d\\\"";
#ifndef WIN32
static const char *gSshCmd = "ssh %s -f4 %s -R %d:localhost:%d sh -c \
   \"'(sh=\\`basename \'\\\\\\$SHELL\'\\`; \
   if test xbash = x\'\\\\\\$sh\' -o xsh = x\'\\\\\\$sh\' -o xzsh = x\'\\\\\\$sh\' -o xdash = x\'\\\\\\$sh\'; then \
      \'\\\\\\$SHELL\' -l -c %s; \
   elif test xcsh = x\'\\\\\\$sh\' -o xtcsh = x\'\\\\\\$sh\' -o xksh = x\'\\\\\\$sh\'; then \
      \'\\\\\\$SHELL\' -c %s; \
   else \
      echo \\\"Unknown shell \'\\\\\\$SHELL\'\\\"; \
   fi)'\"";
#else
static const char *gSshCmd = "ssh %s -f4 %s -R %d:localhost:%d sh -c \
   \"'(sh=`basename $SHELL`; \
   if test xbash = x$sh -o xsh = x$sh -o xzsh = x$sh -o xdash = x$sh; then \
      $SHELL -l -c %s; \
   elif test xcsh = x$sh -o xtcsh = x$sh -o xksh = x$sh; then \
      $SHELL -c %s; \
   else \
      echo \"Unknown shell $SHELL\"; \
   fi)'\"";
#endif

Int_t TApplicationRemote::fgPortAttempts = 100; // number of attempts to find a port
Int_t TApplicationRemote::fgPortLower =  49152; // lower bound for ports
Int_t TApplicationRemote::fgPortUpper =  65535; // upper bound for ports

////////////////////////////////////////////////////////////////////////////////
/// Main constructor: start a remote session at 'url' accepting callbacks
/// on local port 'port'; if port is already in use scan up to 'scan - 1'
/// ports starting from port + 1, i.e. port + 1, ... , port + scan - 1

TApplicationRemote::TApplicationRemote(const char *url, Int_t debug,
                                       const char *script)
                   : TApplication(), fUrl(url)
{
   // Unique name (used also in the prompt)
   fName = fUrl.GetHost();
   if (strlen(fUrl.GetOptions()) > 0)
      fName += Form("-%s", fUrl.GetOptions());
   UserGroup_t *pw = gSystem->GetUserInfo(gSystem->GetEffectiveUid());
   TString user = (pw) ? (const char*) pw->fUser : "";
   SafeDelete(pw);
   if (strlen(fUrl.GetUser()) > 0 && user != fUrl.GetUser())
      fName.Insert(0,Form("%s@", fUrl.GetUser()));

   fIntHandler = 0;
   fSocket = 0;
   fMonitor = 0;
   fFileList = 0;
   fWorkingDir = 0;
   fRootFiles = 0;
   fReceivedObject = 0;
   ResetBit(kCollecting);

   // Create server socket; generate randomly a port to find a free one
   Int_t port = -1;
   Int_t na = fgPortAttempts;
   Long64_t now = gSystem->Now();
   std::default_random_engine randomEngine(now);
   std::uniform_int_distribution<Int_t> randomPort(fgPortLower, fgPortUpper);
   TServerSocket *ss = 0;
   while (na--) {
      port = randomPort(randomEngine);
      ss = new TServerSocket(port);
      if (ss->IsValid())
         break;
   }
   if (!ss || !ss->IsValid()) {
      Error("TApplicationRemote","unable to find a free port for connections");
      SetBit(kInvalidObject);
      return;
   }

   // Create a monitor and add the socket to it
   TMonitor *mon = new TMonitor;
   mon->Add(ss);

   // Start the remote server
   Int_t rport = (port < fgPortUpper) ? port + 1 : port - 1;
   TString sc = gScript;
   if (script && *script) {
      // script is enclosed by " ", so ignore first " char
      if (script[1] == '<') {
         if (script[2])
            sc.Form("source %s; %s", script+2, gScript);
         else
            Error("TApplicationRemote", "illegal script name <");
      } else
         sc = script;
   }
   sc.ReplaceAll("\"","");
   TString userhost = fUrl.GetHost();
   if (strlen(fUrl.GetUser()) > 0)
      userhost.Insert(0, Form("%s@", fUrl.GetUser()));
   const char *verb = "";
   if (debug > 0)
      verb = "-v";
   TString scriptCmd;
   scriptCmd.Form(gScriptCmd, sc.Data(), kRRemote_Protocol, rport, fUrl.GetFile(), debug);
   TString cmd;
   cmd.Form(gSshCmd, verb, userhost.Data(), rport, port, scriptCmd.Data(), scriptCmd.Data());
#ifdef WIN32
   // make sure that the Gpad and GUI libs are loaded
   TApplication::NeedGraphicsLibs();
   gApplication->InitializeGraphics();
#endif
   if (gDebug > 0)
      Info("TApplicationRemote", "executing: %s", cmd.Data());
   if (gSystem->Exec(cmd) != 0) {
      Info("TApplicationRemote", "an error occured during SSH connection");
      mon->DeActivateAll();
      delete mon;
      delete ss;
      SafeDelete(fSocket);
      SetBit(kInvalidObject);
      return;
   }

   // Wait for activity on the socket
   mon->Select();

   // Get the connection
   if (!(fSocket = ss->Accept())) {
      Error("TApplicationRemote", "failed to open connection");
      SetBit(kInvalidObject);
      return;
   }

   // Cleanup the monitor and the server socket
   mon->DeActivateAll();
   delete mon;
   delete ss;

   // Receive the startup message
   Int_t what;
   char buf[512];
   if (fSocket->Recv(buf, sizeof(buf), what) <= 0) {
      Error("TApplicationRemote", "failed to receive startup message");
      SafeDelete(fSocket);
      SetBit(kInvalidObject);
      return;
   }
   Printf("%s", buf);

   // Receive the protocol version run remotely
   if (fSocket->Recv(fProtocol, what) != 2*sizeof(Int_t)) {
      Error("TApplicationRemote", "failed to receive remote server protocol");
      SafeDelete(fSocket);
      SetBit(kInvalidObject);
      return;
   }
   if (fProtocol != kRRemote_Protocol)
      Info("TApplicationRemote","server runs a different protocol version: %d (vs %d)",
                     fProtocol, kRRemote_Protocol);

   TMessage *msg = 0;
   // Receive the protocol version run remotely
   if (fSocket->Recv(msg) < 0 || msg->What() != kMESS_ANY) {
      Error("TApplicationRemote", "failed to receive server info - protocol error");
      SafeDelete(fSocket);
      SetBit(kInvalidObject);
      return;
   }

   // Real host name and full path to remote log
   TString hostname;
   (*msg) >> hostname >> fLogFilePath;
   fUrl.SetHost(hostname);

   // Monitor the socket
   fMonitor = new TMonitor;
   fMonitor->Add(fSocket);

   // Set interrupt handler from now on
   fIntHandler = new TARInterruptHandler(this);

   // To get the right cleaning sequence
   gROOT->GetListOfSockets()->Remove(fSocket);
   gROOT->GetListOfSockets()->Add(this);

   fRootFiles = new TList;
   fRootFiles->SetName("Files");

   // Collect startup notifications
   Collect();

   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TApplicationRemote::~TApplicationRemote()
{
   gROOT->GetListOfSockets()->Remove(this);
   Terminate(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast a message to the remote session.
/// Returns 0 on success, -1 in case of error.

Int_t TApplicationRemote::Broadcast(const TMessage &mess)
{
   if (!IsValid()) return -1;

   if (fSocket->Send(mess) == -1) {
      Error("Broadcast", "could not send message");
      return -1;
   }
   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast a character string buffer to the remote session.
/// Use kind to set the TMessage what field.
/// Returns 0 on success, -1 in case of error.

Int_t TApplicationRemote::Broadcast(const char *str, Int_t kind, Int_t type)
{
   TMessage mess(kind);
   if (kind == kMESS_ANY)
      mess << type;
   if (str) mess.WriteString(str);
   return Broadcast(mess);
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast an object to the remote session.
/// Use kind to set the TMessage what field.
/// Returns 0 on success, -1 in case of error.

Int_t TApplicationRemote::BroadcastObject(const TObject *obj, Int_t kind)
{
   TMessage mess(kind);
   mess.WriteObject(obj);
   return Broadcast(mess);
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast a raw buffer of specified length to the remote session.
/// Returns 0 on success, -1 in case of error.

Int_t TApplicationRemote::BroadcastRaw(const void *buffer, Int_t length)
{
   if (!IsValid()) return -1;

   if (fSocket->SendRaw(buffer, length) == -1) {
      Error("Broadcast", "could not send raw buffer");
      return -1;
   }
   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Collect responses from the remote server.
/// Returns the number of messages received.
/// If timeout >= 0, wait at most timeout seconds (timeout = -1 by default,
/// which means wait forever).

Int_t TApplicationRemote::Collect(Long_t timeout)
{
   // Activate monitoring
   fMonitor->ActivateAll();
   if (!fMonitor->GetActive())
       return 0;

   // Timeout counter
   Long_t nto = timeout;
   if (gDebug > 2)
      Info("Collect","active: %d", fMonitor->GetActive());

   // On clients, handle Ctrl-C during collection
   if (fIntHandler)
      fIntHandler->Add();

   // We are now going to collect from the server
   SetBit(kCollecting);

   Int_t rc = 0, cnt = 0;
   while (fMonitor->GetActive() && (nto < 0 || nto > 0)) {

      // Wait for a ready socket
      TSocket *s = fMonitor->Select(1000);

      if (s && s != (TSocket *)(-1)) {
         // Get and analyse the info it did receive
         if ((rc = CollectInput()) != 0) {
            // Deactivate it if we are done with it
            fMonitor->DeActivate(s);
            if (gDebug > 2)
               Info("Collect","deactivating %p", s);
         }

         // Update counter (if no error occured)
         if (rc >= 0)
            cnt++;

      } else {
         // If not timed-out, exit if not stopped or not aborted
         // (player exits status is finished in such a case); otherwise,
         // we still need to collect the partial output info
         if (!s)
            fMonitor->DeActivateAll();
         // Decrease the timeout counter if requested
         if (s == (TSocket *)(-1) && nto > 0)
            nto--;
      }
   }

   // Collection is over
   ResetBit(kCollecting);

   // If timed-out, deactivate everything
   if (nto == 0)
      fMonitor->DeActivateAll();

   // Deactivate Ctrl-C special handler
   if (fIntHandler)
      fIntHandler->Remove();

   return cnt;
}

////////////////////////////////////////////////////////////////////////////////
/// Collect and analyze available input from the socket.
/// Returns 0 on success, -1 if any failure occurs.

Int_t TApplicationRemote::CollectInput()
{
   TMessage *mess;
   Int_t rc = 0;

   char      str[512];
   TObject  *obj;
   Int_t     what;
   Bool_t    delete_mess = kTRUE;

   if (fSocket->Recv(mess) < 0) {
      SetBit(kInvalidObject);
      SafeDelete(fSocket);
      return -1;
   }
   if (!mess) {
      // we get here in case the remote server died
      SetBit(kInvalidObject);
      SafeDelete(fSocket);
      return -1;
   }

   what = mess->What();

   if (gDebug > 2)
      Info("CollectInput","what %d", what);

   switch (what) {

      case kMESS_OBJECT:
         {  // The server sent over an object: read it in memory
            TObject *o = mess->ReadObject(mess->GetClass());
            // If a canvas, draw it
            if (TString(o->ClassName()) == "TCanvas")
               o->Draw();
            else if (TString(o->ClassName()) == "TRemoteObject") {
               TRemoteObject *robj = (TRemoteObject *)o;
               if (TString(robj->GetClassName()) == "TSystemDirectory") {
                  if (fWorkingDir == 0) {
                     fWorkingDir = (TRemoteObject *)o;
                  }
               }
            }
            else if (TString(o->ClassName()) == "TList") {
               TList *list = (TList *)o;
               TRemoteObject *robj = (TRemoteObject *)list->First();
               if (robj && (TString(robj->GetClassName()) == "TFile")) {
                  TIter next(list);
                  while ((robj = (TRemoteObject *)next())) {
                     if (!fRootFiles->FindObject(robj->GetName()))
                        fRootFiles->Add(robj);
                  }
                  gROOT->RefreshBrowsers();
               }
            }
            fReceivedObject = o;
         }
         break;

      case kMESS_ANY:
         // Generic message: read out the type
         {  Int_t type;
            (*mess) >> type;

            if (gDebug > 2)
               Info("CollectInput","type %d", type);

            switch (type) {

               case kRRT_GetObject:
                  // send server the object it asks for
                  mess->ReadString(str, sizeof(str));
                  obj = gDirectory->Get(str);
                  if (obj) {
                     fSocket->SendObject(obj);
                  } else {
                     Warning("CollectInput",
                             "server requested an object that we do not have");
                     fSocket->Send(kMESS_NOTOK);
                  }
                  break;

               case kRRT_Fatal:
                  // Fatal error
                  SafeDelete(fSocket);
                  rc = -1;
                  break;

               case kRRT_LogFile:
                  {  Int_t size;
                     (*mess) >> size;
                     RecvLogFile(size);
                  }
                  break;

               case kRRT_LogDone:
                  {  Int_t st;
                    (*mess) >> st;
                     if (st < 0) {
                        // Problem: object should not be used
                        SetBit(kInvalidObject);
                     }
                     if (gDebug > 1)
                        Info("CollectInput","kRTT_LogDone: status %d", st);
                     rc = 1;
                  }
                  break;

               case kRRT_Message:
                  {  TString msg;
                     Bool_t lfeed;
                     (*mess) >> msg >> lfeed;
                     if (lfeed)
                        fprintf(stderr,"%s\n", msg.Data());
                     else
                        fprintf(stderr,"%s\r", msg.Data());
                  }
                  break;

               case kRRT_SendFile:
                  {  TString fname;
                     (*mess) >> fname;
                     // Prepare the reply
                     TMessage m(kMESS_ANY);
                     m << (Int_t) kRRT_SendFile;
                     // The server needs a file: we send also the related header
                     // if we have it.
                     char *imp = gSystem->Which(TROOT::GetMacroPath(), fname, kReadPermission);
                     if (!imp) {
                        Error("CollectInput", "file %s not found in path(s) %s",
                                         fname.Data(), TROOT::GetMacroPath());
                        m << (Bool_t) kFALSE;
                        Broadcast(m);
                     } else {
                        TString impfile = imp;
                        delete [] imp;
                        Int_t dot = impfile.Last('.');

                        // Is there any associated header file
                        Bool_t hasHeader = kTRUE;
                        TString headfile = impfile;
                        if (dot != kNPOS)
                           headfile.Remove(dot);
                        headfile += ".h";
                        if (gSystem->AccessPathName(headfile, kReadPermission)) {
                           TString h = headfile;
                           headfile.Remove(dot);
                           headfile += ".hh";
                           if (gSystem->AccessPathName(headfile, kReadPermission)) {
                              hasHeader = kFALSE;
                              if (gDebug > 0)
                                 Info("CollectInput", "no associated header file"
                                                 " found: tried: %s %s",
                                                 h.Data(), headfile.Data());
                           }
                        }

                        // Send files now;
                        m << (Bool_t) kTRUE;
                        Broadcast(m);
                        if (SendFile(impfile, kForce) == -1) {
                           Info("CollectInput", "problems sending file %s", impfile.Data());
                           return 0;
                        }
                        if (hasHeader) {
                           Broadcast(m);
                           if (SendFile(headfile, kForce) == -1) {
                              Info("CollectInput", "problems sending file %s", headfile.Data());
                              return 0;
                           }
                        }
                     }
                     // End of transmission
                     m.Reset(kMESS_ANY);
                     m << (Int_t) kRRT_SendFile;
                     m << (Bool_t) kFALSE;
                     Broadcast(m);
                  }
                  break;

               default:
                  Warning("CollectInput","unknown type received from server: %d", type);
                  break;

            }
         }
         break;

      default:
         Error("CollectInput", "unknown command received from server: %d", what);
         SetBit(kInvalidObject);
         SafeDelete(fSocket);
         rc = -1;
         break;
   }

   // Cleanup
   if (delete_mess)
      delete mess;

   // We are done successfully
   return rc;
}


////////////////////////////////////////////////////////////////////////////////
/// Receive the log file from the server

void TApplicationRemote::RecvLogFile(Int_t size)
{
   const Int_t kMAXBUF = 16384;  //32768  //16384  //65536;
   char buf[kMAXBUF];

   // Append messages to active logging unit
   Int_t fdout = fileno(stdout);
   if (fdout < 0) {
      Warning("RecvLogFile", "file descriptor for outputs undefined (%d):"
                             " will not log msgs", fdout);
      return;
   }
   lseek(fdout, (off_t) 0, SEEK_END);

   Int_t  left, rec, r;
   Long_t filesize = 0;

   while (filesize < size) {
      left = Int_t(size - filesize);
      if (left > kMAXBUF)
         left = kMAXBUF;
      rec = fSocket->RecvRaw(&buf, left);
      filesize = (rec > 0) ? (filesize + rec) : filesize;
      if (rec > 0) {

         char *p = buf;
         r = rec;
         while (r) {
            Int_t w;

            w = write(fdout, p, r);

            if (w < 0) {
               SysError("RecvLogFile", "error writing to unit: %d", fdout);
               break;
            }
            r -= w;
            p += w;
         }
      } else if (rec < 0) {
         Error("RecvLogFile", "error during receiving log file");
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Send object to server.
/// Return 0 on success, -1 in case of error.

Int_t TApplicationRemote::SendObject(const TObject *obj)
{
   if (!IsValid() || !obj) return -1;

   TMessage mess(kMESS_OBJECT);
   mess.WriteObject(obj);
   return Broadcast(mess);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if a file needs to be send to the server. Use the following
/// algorithm:
///   - check if file appears in file map
///     - if yes, get file's modtime and check against time in map,
///       if modtime not same get md5 and compare against md5 in map,
///       if not same return kTRUE.
///     - if no, get file's md5 and modtime and store in file map, ask
///       slave if file exists with specific md5, if yes return kFALSE,
///       if no return kTRUE.
/// Returns kTRUE in case file needs to be send, returns kFALSE in case
/// file is already on remote node.

Bool_t TApplicationRemote::CheckFile(const char *file, Long_t modtime)
{
   Bool_t sendto = kFALSE;

   if (!IsValid()) return -1;

   // The filename for the cache
   TString fn = gSystem->BaseName(file);

   // Check if the file is already in the cache
   TARFileStat *fs = 0;
   if (fFileList && (fs = (TARFileStat *) fFileList->FindObject(fn))) {
      // File in cache
      if (fs->fModtime != modtime) {
         TMD5 *md5 = TMD5::FileChecksum(file);
         if (md5) {
            if ((*md5) != fs->fMD5) {
               sendto       = kTRUE;
               fs->fMD5      = *md5;
               fs->fModtime  = modtime;
            }
            delete md5;
         } else {
            Error("CheckFile", "could not calculate local MD5 check sum - dont send");
            return kFALSE;
         }
      }
   } else {
      // file not in the cache
      TMD5 *md5 = TMD5::FileChecksum(file);
      if (md5) {
         fs = new TARFileStat(fn, md5, modtime);
         if (!fFileList)
            fFileList = new THashList;
         fFileList->Add(fs);
         delete md5;
      } else {
         Error("CheckFile", "could not calculate local MD5 check sum - dont send");
         return kFALSE;
      }
      TMessage mess(kMESS_ANY);
      mess << Int_t(kRRT_CheckFile) << TString(gSystem->BaseName(file)) << fs->fMD5;
      fSocket->Send(mess);

      TMessage *reply;
      fSocket->Recv(reply);
      if (reply) {
         if (reply->What() == kMESS_ANY) {
            // Check the type
            Int_t type;
            Bool_t uptodate;
            (*reply) >> type >> uptodate;
            if (type != kRRT_CheckFile) {
               // Protocol error
               Warning("CheckFile", "received wrong type:"
                                    " %d (expected %d): protocol error?",
                                    type, (Int_t)kRRT_CheckFile);
            }
            sendto = uptodate ? kFALSE : kTRUE;
         } else {
            // Protocol error
            Error("CheckFile", "received wrong message: %d (expected %d)",
                               reply->What(), kMESS_ANY);
         }
      } else {
         Error("CheckFile", "received empty message");
      }
      // Collect logs
      Collect();
   }

   // Done
   return sendto;
}

////////////////////////////////////////////////////////////////////////////////
/// Send a file to the server. Return 0 on success, -1 in case of error.
/// If defined, the full path of the remote path will be rfile.
/// The mask 'opt' is an or of ESendFileOpt:
///
///       kAscii  (0x0)      if set true ascii file transfer is used
///       kBinary (0x1)      if set true binary file transfer is used
///       kForce  (0x2)      if not set an attempt is done to find out
///                          whether the file really needs to be downloaded
///                          (a valid copy may already exist in the cache
///                          from a previous run)
///

Int_t TApplicationRemote::SendFile(const char *file, Int_t opt, const char *rfile)
{
   if (!IsValid()) return -1;

#ifndef R__WIN32
   Int_t fd = open(file, O_RDONLY);
#else
   Int_t fd = open(file, O_RDONLY | O_BINARY);
#endif
   if (fd < 0) {
      SysError("SendFile", "cannot open file %s", file);
      return -1;
   }

   // Get info about the file
   Long64_t size;
   Long_t id, flags, modtime;
   if (gSystem->GetPathInfo(file, &id, &size, &flags, &modtime) == 1) {
      Error("SendFile", "cannot stat file %s", file);
      close(fd);
      return -1;
   }
   if (size == 0) {
      Error("SendFile", "empty file %s", file);
      close(fd);
      return -1;
   }

   // Decode options
   Bool_t bin   = (opt & kBinary)  ? kTRUE : kFALSE;
   Bool_t force = (opt & kForce)   ? kTRUE : kFALSE;

   const Int_t kMAXBUF = 32768;  //16384  //65536;
   char buf[kMAXBUF];

   const char *fnam = (rfile) ? rfile : gSystem->BaseName(file);

   Bool_t sendto = force ? kTRUE : CheckFile(file, modtime);

   // The value of 'size' is used as flag remotely, so we need to
   // reset it to 0 if we are not going to send the file
   size = sendto ? size : 0;

   if (gDebug > 1 && size > 0)
      Info("SendFile", "sending file %s", file);

   snprintf(buf, kMAXBUF, "%s %d %lld", fnam, bin, size);
   if (Broadcast(buf, kMESS_ANY, kRRT_File) == -1) {
      SafeDelete(fSocket);
      close(fd);
      return -1;
   }

   if (sendto) {

      lseek(fd, 0, SEEK_SET);

      Int_t len;
      do {
         while ((len = read(fd, buf, kMAXBUF)) < 0 && TSystem::GetErrno() == EINTR)
            TSystem::ResetErrno();

         if (len < 0) {
            SysError("SendFile", "error reading from file %s", file);
            Interrupt();
            close(fd);
            return -1;
         }

         if (len > 0 && fSocket->SendRaw(buf, len) == -1) {
            SysError("SendFile", "error writing to server @ %s:%d (now offline)",
                     fUrl.GetHost(), fUrl.GetPort());
            SafeDelete(fSocket);
            break;
         }

      } while (len > 0);
   }
   close(fd);

   // Get the log (during collection this will be done at the end
   if (!TestBit(kCollecting))
      Collect();

   // Done
   return IsValid() ? 0 : -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Terminate this session

void TApplicationRemote::Terminate(Int_t status)
{
   TMessage mess(kMESS_ANY);
   mess << (Int_t)kRRT_Terminate << status;
   Broadcast(mess);

   SafeDelete(fRootFiles);
   SafeDelete(fMonitor);
   SafeDelete(fSocket);
}

////////////////////////////////////////////////////////////////////////////////
/// Set port parameters for tunnelling. A value of -1 means unchanged

void TApplicationRemote::SetPortParam(Int_t lower, Int_t upper, Int_t attempts)
{
   if (lower > -1)
      fgPortLower = lower;
   if (upper > -1)
      fgPortUpper = upper;
   if (attempts > -1)
      fgPortAttempts = attempts;

   ::Info("TApplicationRemote::SetPortParam","port scan: %d attempts in [%d,%d]",
          fgPortAttempts, fgPortLower, fgPortUpper);
}

////////////////////////////////////////////////////////////////////////////////
/// Parse a single command line and forward the request to the remote server
/// where it will be processed. The line is either a C++ statement or an
/// interpreter command starting with a ".".
/// Return the return value of the command casted to a long.

Long_t TApplicationRemote::ProcessLine(const char *line, Bool_t, Int_t *)
{
   if (!line || !*line) return 0;

   if (!strncasecmp(line, ".q", 2)) {
      // terminate the session
      gApplication->ProcessLine(".R -close");
      return 0;
   }

   if (!strncmp(line, "?", 1)) {
      Help(line);
      return 1;
   }

   fReceivedObject = 0;

   // Init graphics
   InitializeGraphics();

   // Ok, now we pack the command and we send it over for processing
   Broadcast(line, kMESS_CINT);

   // And collect the results
   Collect();

   // Done
   return (Long_t)fReceivedObject;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Print some info about this instance

void TApplicationRemote::Print(Option_t *opt) const
{
   TString s(Form("OBJ: TApplicationRemote     %s", fName.Data()));
   Printf("%s", s.Data());
   if (opt && opt[0] == 'F') {
      s = "    url: ";
      if (strlen(fUrl.GetUser()) > 0)
         s += Form("%s@", fUrl.GetUser());
      s += fUrl.GetHostFQDN();
      s += Form("  logfile: %s", fLogFilePath.Data());
      Printf("%s", s.Data());
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Send interrupt OOB byte to server.
/// Returns 0 if ok, -1 in case of error

void TApplicationRemote::Interrupt(Int_t type)
{
   if (!IsValid()) return;

   fInterrupt = kTRUE;

#if 1
   Info("Interrupt", "*** Ctrl-C not yet enabled *** (type= %d)", type);
   return;
#else

   char oobc = (char) type;
   const int kBufSize = 1024;
   char waste[kBufSize];

   // Send one byte out-of-band message to server
   if (fSocket->SendRaw(&oobc, 1, kOob) <= 0) {
      Error("Interrupt", "error sending oobc to server");
      return;
   }

   if (type == kRRI_Hard) {
      char  oob_byte;
      int   n, nch, nbytes = 0, nloop = 0;

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
               Error("Interrupt", "error receiving waste from server");
               break;
            }
            nbytes += n;
         } else if (n == -3) {   // EINVAL
            //
            // The OOB data has not arrived yet
            //
            gSystem->Sleep(100);
            if (++nloop > 100) {  // 10 seconds time-out
               Error("Interrupt", "server does not respond");
               break;
            }
         } else {
            Error("Interrupt", "error receiving OOB from server");
            break;
         }
      }

      //
      // Continue flushing the input socket stream until the OOB
      // mark is reached
      //
      while (1) {
         int atmark;

         fSocket->GetOption(kAtMark, atmark);

         if (atmark)
            break;

         // find out number of bytes to read before atmark
         fSocket->GetOption(kBytesToRead, nch);
         if (nch == 0) {
            gSystem->Sleep(1000);
            continue;
         }

         if (nch > kBufSize) nch = kBufSize;
         n = fSocket->RecvRaw(waste, nch);
         if (n <= 0) {
            Error("Interrupt", "error receiving waste (2) from server");
            break;
         }
         nbytes += n;
      }
      if (nbytes > 0)
         Info("Interrupt", "server synchronized: %d bytes discarded", nbytes);

      // Get log file from server after a hard interrupt
      Collect();

   } else if (type == kRRI_Soft) {

      // Get log file from server after a soft interrupt
      Collect();

   } else if (type == kRRI_Shutdown) {

      ; // nothing expected to be returned

   } else {

      // Unexpected message, just receive log file
      Collect();
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Browse remote application (working directory and ROOT files).

void TApplicationRemote::Browse(TBrowser *b)
{
   b->Add(fRootFiles, "ROOT Files");
   b->Add(fWorkingDir, fWorkingDir->GetTitle());
   gROOT->RefreshBrowsers();
}
