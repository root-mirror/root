// @(#)root/proofx:$Id$
// Author: Gerardo Ganis  12/12/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
/**
  \defgroup proofx XProofD client Library
  \ingroup proof

  The XProofD client library, libProofx, contain the classes providing 
  the client to interact with the XRootD-based xproofd daemon.

*/

/** \class TXProofMgr
\ingroup proofx

Implementation of the functionality provided by TProofMgr in the case of a xproofd-based session.

*/

#include <errno.h>
#include <memory>
#ifdef WIN32
#include <io.h>
#endif

#include "Getline.h"
#include "snprintf.h"
#include "TList.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TProof.h"
#include "TProofLog.h"
#include "TXProofMgr.h"
#include "TXSocket.h"
#include "TXSocketHandler.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TSysEvtHandler.h"
#include "XProofProtocol.h"

#include "XrdProofConn.h"

ClassImp(TXProofMgr);

//
//----- ProofMgr Interrupt signal handler
//
class TProofMgrInterruptHandler : public TSignalHandler {
private:
   TProofMgr *fMgr;

   TProofMgrInterruptHandler(const TProofMgrInterruptHandler&); // Not implemented
   TProofMgrInterruptHandler& operator=(const TProofMgrInterruptHandler&); // Not implemented
public:
   TProofMgrInterruptHandler(TProofMgr *mgr)
      : TSignalHandler(kSigInterrupt, kFALSE), fMgr(mgr) { }
   Bool_t Notify();
};

////////////////////////////////////////////////////////////////////////////////
/// TProofMgr interrupt handler.

Bool_t TProofMgrInterruptHandler::Notify()
{
   // Only on clients
   if (isatty(0) != 0 && isatty(1) != 0) {
      TString u = fMgr->GetUrl();
      Printf("Opening new connection to %s", u.Data());
      TXSocket *s = new TXSocket(u, 'C', kPROOF_Protocol,
                                 kXPROOF_Protocol, 0, -1, (TXHandler *)fMgr);
      if (s && s->IsValid()) {
         // Set the interrupt flag on the server
         s->CtrlC();
      }
   }
   return kTRUE;
}

// AutoLoading hooks.
// These are needed to avoid using the plugin manager which may create
// problems in multi-threaded environments.
TProofMgr *GetTXProofMgr(const char *url, Int_t l, const char *al)
{ return ((TProofMgr *) new TXProofMgr(url, l, al)); }

class TXProofMgrInit {
public:
   TXProofMgrInit() {
      TProofMgr::SetTXProofMgrHook(&GetTXProofMgr);
}};
static TXProofMgrInit gxproofmgr_init;

////////////////////////////////////////////////////////////////////////////////
/// Create a PROOF manager for the standard (old) environment.

TXProofMgr::TXProofMgr(const char *url, Int_t dbg, const char *alias)
          : TProofMgr(url, dbg, alias)
{
   // Set the correct servert type
   fServType = kXProofd;

   // Initialize
   if (Init(dbg) != 0) {
      // Failure: make sure the socket is deleted so that its lack of
      // validity is correctly transmitted
      SafeDelete(fSocket);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Do real initialization: open the connection and set the relevant
/// variables.
/// Login and authentication are dealt with at this level, if required.
/// Return 0 in case of success, 1 if the remote server is a 'proofd',
/// -1 in case of error.

Int_t TXProofMgr::Init(Int_t)
{
   // Here we make sure that the port is explicitly specified in the URL,
   // even when it matches the default value
   TString u = fUrl.GetUrl(kTRUE);

   fSocket = 0;
   if (!(fSocket = new TXSocket(u, 'C', kPROOF_Protocol,
                                kXPROOF_Protocol, 0, -1, this)) ||
       !(fSocket->IsValid())) {
      if (!fSocket || !(fSocket->IsServProofd()))
         if (gDebug > 0)
            Error("Init", "while opening the connection to %s - exit (error: %d)",
                          u.Data(), (fSocket ? fSocket->GetOpenError() : -1));
      if (fSocket && fSocket->IsServProofd())
         fServType = TProofMgr::kProofd;
      return -1;
   }

   // Protocol run by remote PROOF server
   fRemoteProtocol = fSocket->GetRemoteProtocol();

   // We add the manager itself for correct destruction
   {  R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfSockets()->Remove(fSocket);
   }

   // Set interrupt PROOF handler from now on
   fIntHandler = new TProofMgrInterruptHandler(this);

   // We are done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor: close the connection

TXProofMgr::~TXProofMgr()
{
   SetInvalid();
}

////////////////////////////////////////////////////////////////////////////////
/// Invalidate this manager by closing the connection

void TXProofMgr::SetInvalid()
{
   if (fSocket)
      fSocket->Close("P");
   SafeDelete(fSocket);

   // Avoid destroying twice
   {  R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfSockets()->Remove(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Dummy version provided for completeness. Just returns a pointer to
/// existing session 'id' (as shown by TProof::QuerySessions) or 0 if 'id' is
/// not valid. The boolena 'gui' should be kTRUE when invoked from the GUI.

TProof *TXProofMgr::AttachSession(TProofDesc *d, Bool_t gui)
{
   if (!IsValid()) {
      Warning("AttachSession","invalid TXProofMgr - do nothing");
      return 0;
   }
   if (!d) {
      Warning("AttachSession","invalid description object - do nothing");
      return 0;
   }

   if (d->GetProof())
      // Nothing to do if already in contact with proofserv
      return d->GetProof();

   // Re-compose url
   TString u(Form("%s/?%d", fUrl.GetUrl(kTRUE), d->GetRemoteId()));

   // We need this to set correctly the kUsingSessionGui bit before the first
   // feedback messages arrive
   if (gui)
      u += "GUI";

   // Attach
   TProof *p = new TProof(u, 0, 0, gDebug, 0, this);
   if (p && p->IsValid()) {

      // Set reference manager
      p->SetManager(this);

      // Save record about this session
      Int_t st = (p->IsIdle()) ? TProofDesc::kIdle
                                 : TProofDesc::kRunning;
      d->SetStatus(st);
      d->SetProof(p);

      // Set session tag
      p->SetName(d->GetName());

   } else {
      // Session creation failed
      Error("AttachSession", "attaching to PROOF session");
   }
   return p;
}

////////////////////////////////////////////////////////////////////////////////
/// Detach session with 'id' from its proofserv. The 'id' is the number
/// shown by QuerySessions. The correspondent TProof object is deleted.
/// If id == 0 all the known sessions are detached.
/// Option opt="S" or "s" forces session shutdown.

void TXProofMgr::DetachSession(Int_t id, Option_t *opt)
{
   if (!IsValid()) {
      Warning("DetachSession","invalid TXProofMgr - do nothing");
      return;
   }

   if (id > 0) {
      // Single session request
      TProofDesc *d = GetProofDesc(id);
      if (d) {
         if (fSocket)
            fSocket->DisconnectSession(d->GetRemoteId(), opt);
         TProof *p = d->GetProof();
         fSessions->Remove(d);
         SafeDelete(p);
         delete d;
      }
   } else if (id == 0) {

      // Requesto to destroy all sessions
      if (fSocket) {
         TString o = Form("%sA",opt);
         fSocket->DisconnectSession(-1, o);
      }
      if (fSessions) {
         // Delete PROOF sessions
         TIter nxd(fSessions);
         TProofDesc *d = 0;
         while ((d = (TProofDesc *)nxd())) {
            TProof *p = d->GetProof();
            SafeDelete(p);
         }
         fSessions->Delete();
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Detach session 'p' from its proofserv. The instance 'p' is invalidated
/// and should be deleted by the caller

void TXProofMgr::DetachSession(TProof *p, Option_t *opt)
{
   if (!IsValid()) {
      Warning("DetachSession","invalid TXProofMgr - do nothing");
      return;
   }

   if (p) {
      // Single session request
      TProofDesc *d = GetProofDesc(p);
      if (d) {
         if (fSocket)
            fSocket->DisconnectSession(d->GetRemoteId(), opt);
         fSessions->Remove(d);
         p->Close(opt);
         delete d;
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if 'url' refers to the same 'user@host:port' entity as the URL
/// in memory. TProofMgr::MatchUrl cannot be used here because of the
/// 'double' default port, implying an additional check on the port effectively
/// open.

Bool_t TXProofMgr::MatchUrl(const char *url)
{
   if (!IsValid()) {
      Warning("MatchUrl","invalid TXProofMgr - do nothing");
      return 0;
   }

   TUrl u(url);

   // Correct URL protocol
   if (!strcmp(u.GetProtocol(), TUrl("a").GetProtocol()))
      u.SetProtocol("proof");

   if (u.GetPort() == TUrl("a").GetPort()) {
      // Set default port
      Int_t port = gSystem->GetServiceByName("proofd");
      if (port < 0)
         port = 1093;
      u.SetPort(port);
   }

   // Now we can check
   if (!strcmp(u.GetHostFQDN(), fUrl.GetHost()))
      if (u.GetPort() == fUrl.GetPort() ||
          u.GetPort() == fSocket->GetPort())
         if (strlen(u.GetUser()) <= 0 || !strcmp(u.GetUser(),fUrl.GetUser()))
            return kTRUE;

   // Match failed
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Show available workers

void TXProofMgr::ShowWorkers()
{
   if (!IsValid()) {
      Warning("ShowWorkers","invalid TXProofMgr - do nothing");
      return;
   }

   // Send the request
   TObjString *os = fSocket->SendCoordinator(kQueryWorkers);
   if (os) {
      TObjArray *oa = TString(os->GetName()).Tokenize(TString("&"));
      if (oa) {
         TIter nxos(oa);
         TObjString *to = 0;
         while ((to = (TObjString *) nxos()))
            // Now parse them ...
            Printf("+  %s", to->GetName());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Gets the URL to be prepended to paths when accessing the MSS associated
/// with the connected cluster, if any. The information is retrieved from
/// the cluster the first time or if retrieve is true.

const char *TXProofMgr::GetMssUrl(Bool_t retrieve)
{
   if (fMssUrl.IsNull() || retrieve) {
      // Nothing to do if not in contact with proofserv
      if (!IsValid()) {
         Error("GetMssUrl", "invalid TXProofMgr - do nothing");
         return 0;
      }
      // Server may not support it
      if (fSocket->GetXrdProofdVersion() < 1007) {
         Error("GetMssUrl", "functionality not supported by server");
         return 0;
      }
      TObjString *os = fSocket->SendCoordinator(kQueryMssUrl);
      if (os) {
         Printf("os: '%s'", os->GetName());
         fMssUrl = os->GetName();
         SafeDelete(os);
      } else {
         Error("GetMssUrl", "problems retrieving the required information");
         return 0;
      }
   } else if (!IsValid()) {
      Warning("GetMssUrl", "TXProofMgr is now invalid: information may not be valid");
      return 0;
   }

   // Done
   return fMssUrl.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Get list of sessions accessible to this manager

TList *TXProofMgr::QuerySessions(Option_t *opt)
{
   if (opt && !strncasecmp(opt,"L",1))
      // Just return the existing list
      return fSessions;

   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Warning("QuerySessions","invalid TXProofMgr - do nothing");
      return 0;
   }

   // Create list if not existing
   if (!fSessions) {
      fSessions = new TList();
      fSessions->SetOwner();
   }

   // Send the request
   TList *ocl = new TList;
   TObjString *os = fSocket->SendCoordinator(kQuerySessions);
   if (os) {
      TObjArray *oa = TString(os->GetName()).Tokenize(TString("|"));
      if (oa) {
         TProofDesc *d = 0;
         TIter nxos(oa);
         TObjString *to = (TObjString *) nxos();
         if (to && to->GetString().IsDigit() && !strncasecmp(opt,"S",1))
            Printf("// +++ %s session(s) currently active +++", to->GetName());
         while ((to = (TObjString *) nxos())) {
            // Now parse them ...
            Int_t id = -1, st = -1;
            TString al, tg, tk;
            Ssiz_t from = 0;
            while (to->GetString()[from] == ' ') { from++; }
            if (!to->GetString().Tokenize(tk, from, " ") || !tk.IsDigit()) continue;
            id = tk.Atoi();
            if (!to->GetString().Tokenize(tg, from, " ")) continue;
            if (!to->GetString().Tokenize(al, from, " ")) continue;
            if (!to->GetString().Tokenize(tk, from, " ") || !tk.IsDigit()) continue;
            st = tk.Atoi();
            // Add to the list, if not already there
            if (!(d = (TProofDesc *) fSessions->FindObject(tg))) {
               Int_t locid = fSessions->GetSize() + 1;
               d = new TProofDesc(tg, al, GetUrl(), locid, id, st, 0);
               fSessions->Add(d);
            } else {
               // Set missing / update info
               d->SetStatus(st);
               d->SetRemoteId(id);
               d->SetTitle(al);
            }
            // Add to the list for final garbage collection
            ocl->Add(new TObjString(tg));
         }
         SafeDelete(oa);
      }
      SafeDelete(os);
   }

   // Printout and Garbage collection
   if (fSessions->GetSize() > 0) {
      TIter nxd(fSessions);
      TProofDesc *d = 0;
      while ((d = (TProofDesc *)nxd())) {
         if (ocl->FindObject(d->GetName())) {
            if (opt && !strncasecmp(opt,"S",1))
               d->Print("");
         } else {
            fSessions->Remove(d);
            SafeDelete(d);
         }
      }
   }

   // We are done
   return fSessions;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle asynchronous input on the socket

Bool_t TXProofMgr::HandleInput(const void *)
{
   if (fSocket && fSocket->IsValid()) {
      TMessage *mess;
      if (fSocket->Recv(mess) >= 0) {
         Int_t what = mess->What();
         if (gDebug > 0)
            Info("HandleInput", "%p: got message type: %d", this, what);
         switch (what) {
            case kPROOF_TOUCH:
               fSocket->RemoteTouch();
               break;
            default:
               Warning("HandleInput", "%p: got unknown message type: %d", this, what);
               break;
         }
      }
   } else {
      Warning("HandleInput", "%p: got message but socket is invalid!", this);
   }

   // We are done
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle error on the input socket

Bool_t TXProofMgr::HandleError(const void *in)
{
   XHandleErr_t *herr = in ? (XHandleErr_t *)in : 0;

   // Try reconnection
   if (fSocket && herr && (herr->fOpt == 1)) {
      fSocket->Reconnect();
      if (fSocket && fSocket->IsValid()) {
         if (gDebug > 0)
            Printf("ProofMgr: connection to coordinator at %s re-established",
                   fUrl.GetUrl());
         return kFALSE;
      }
   }
   Printf("TXProofMgr::HandleError: %p: got called ...", this);

   // Interrupt any PROOF session in Collect
   if (fSessions && fSessions->GetSize() > 0) {
      TIter nxd(fSessions);
      TProofDesc *d = 0;
      while ((d = (TProofDesc *)nxd())) {
         TProof *p = (TProof *) d->GetProof();
         if (p)
            p->InterruptCurrentMonitor();
      }
   }
   if (gDebug > 0)
      Printf("TXProofMgr::HandleError: %p: DONE ... ", this);

   // We are done
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Send a cleanup request for the sessions associated with the current user.
/// If 'hard' is true sessions are signalled for termination and moved to
/// terminate at all stages (top master, sub-master, workers). Otherwise
/// (default) only top-master sessions are asked to terminate, triggering
/// a gentle session termination. In all cases all sessions should be gone
/// after a few (2 or 3) session checking cycles.
/// A user with superuser privileges can also asks cleaning for an different
/// user, specified by 'usr', or for all users (usr = *)
/// Return 0 on success, -1 in case of error.

Int_t TXProofMgr::Reset(Bool_t hard, const char *usr)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Warning("Reset","invalid TXProofMgr - do nothing");
      return -1;
   }

   Int_t h = (hard) ? 1 : 0;
   fSocket->SendCoordinator(kCleanupSessions, usr, h);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get logs or log tails from last session associated with this manager
/// instance.
/// The arguments allow to specify a session different from the last one:
///      isess   specifies a position relative to the last one, i.e. 1
///              for the next to last session; the absolute value is taken
///              so -1 and 1 are equivalent.
///      stag    specifies the unique tag of the wanted session
/// The special value stag = "NR" allows to just initialize the TProofLog
/// object w/o retrieving the files; this may be useful when the number
/// of workers is large and only a subset of logs is required.
/// If 'stag' is specified 'isess' is ignored (unless stag = "NR").
/// If 'pattern' is specified only the lines containing it are retrieved
/// (remote grep functionality); to filter out a pattern 'pat' use
/// pattern = "-v pat".
/// If 'rescan' is TRUE, masters will rescan the worker sandboxes for the exact
/// paths, instead of using the save information; may be useful when the
/// ssave information looks wrong or incomplete.
/// Returns a TProofLog object (to be deleted by the caller) on success,
/// 0 if something wrong happened.

TProofLog *TXProofMgr::GetSessionLogs(Int_t isess, const char *stag,
                                      const char *pattern, Bool_t rescan)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Warning("GetSessionLogs","invalid TXProofMgr - do nothing");
      return 0;
   }

   TProofLog *pl = 0;

   // The absolute value of isess counts
   isess = (isess > 0) ? -isess : isess;

   // Special option in stag
   bool retrieve = 1;
   TString sesstag(stag);
   if (sesstag == "NR") {
      retrieve = 0;
      sesstag = "";
   }

   // Get the list of paths
   Int_t xrs = (rescan) ? 1 : 0;
   TObjString *os = fSocket->SendCoordinator(kQueryLogPaths, sesstag.Data(), isess, -1, xrs);

   // Analyse it now
   Int_t ii = 0;
   if (os) {
      TString rs(os->GetName());
      Ssiz_t from = 0;
      // The session tag
      TString tag;
      if (!rs.Tokenize(tag, from, "|")) {
         Warning("GetSessionLogs", "Session tag undefined: corruption?\n"
                                   " (received string: %s)", os->GetName());
         return (TProofLog *)0;
      }
      // The pool url
      TString purl;
      if (!rs.Tokenize(purl, from, "|")) {
         Warning("GetSessionLogs", "Pool URL undefined: corruption?\n"
                                   " (received string: %s)", os->GetName());
         return (TProofLog *)0;
      }
      // Create the instance now
      if (!pl)
         pl = new TProofLog(tag, GetUrl(), this);

      // Per-node info
      TString to;
      while (rs.Tokenize(to, from, "|")) {
         if (!to.IsNull()) {
            TString ord(to);
            ord.Strip(TString::kLeading, ' ');
            TString url(ord);
            if ((ii = ord.Index(" ")) != kNPOS)
               ord.Remove(ii);
            if ((ii = url.Index(" ")) != kNPOS)
               url.Remove(0, ii + 1);
            // Add to the list (special tag for valgrind outputs)
            if (url.Contains(".valgrind")) ord += "-valgrind";
            pl->Add(ord, url);
            // Notify
            if (gDebug > 1)
               Info("GetSessionLogs", "ord: %s, url: %s", ord.Data(), url.Data());
         }
      }
      // Cleanup
      SafeDelete(os);
      // Retrieve the default part if required
      if (pl && retrieve) {
         const char *pat = pattern ? pattern : "-v \"| SvcMsg\"";
         if (pat && strlen(pat) > 0)
            pl->Retrieve("*", TProofLog::kGrep, 0, pat);
         else
            pl->Retrieve();
      }
   }

   // Done
   return pl;
}

////////////////////////////////////////////////////////////////////////////////
/// Read, via the coordinator, 'len' bytes from offset 'ofs' of 'file'.
/// Returns a TObjString with the content or 0, in case of failure

TObjString *TXProofMgr::ReadBuffer(const char *fin, Long64_t ofs, Int_t len)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Warning("ReadBuffer","invalid TXProofMgr - do nothing");
      return (TObjString *)0;
   }

   // Send the request
   return fSocket->SendCoordinator(kReadBuffer, fin, len, ofs, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Read, via the coordinator, 'fin' filtered. If 'pattern' starts with '|',
/// it represents a command filtering the output. Elsewhere, it is a grep
/// pattern. Returns a TObjString with the content or 0 in case of failure

TObjString *TXProofMgr::ReadBuffer(const char *fin, const char *pattern)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Warning("ReadBuffer", "invalid TXProofMgr - do nothing");
      return (TObjString *)0;
   }

   const char *ptr;
   Int_t type;  // 1 = grep, 2 = grep -v, 3 = pipe through cmd
   if (*pattern == '|') {
      ptr = &pattern[1];  // strip first char if it is a command
      type = 3;
   }
   else {
      ptr = pattern;
      type = 1;
   }

   // Prepare the buffer
   Int_t plen = strlen(ptr);
   Int_t lfi = strlen(fin);
   char *buf = new char[lfi + plen + 1];
   memcpy(buf, fin, lfi);
   memcpy(buf+lfi, ptr, plen);
   buf[lfi+plen] = 0;

   // Send the request
   return fSocket->SendCoordinator(kReadBuffer, buf, plen, 0, type);
}

////////////////////////////////////////////////////////////////////////////////
/// Display what ROOT versions are available on the cluster

void TXProofMgr::ShowROOTVersions()
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Warning("ShowROOTVersions","invalid TXProofMgr - do nothing");
      return;
   }

   // Send the request
   TObjString *os = fSocket->SendCoordinator(kQueryROOTVersions);
   if (os) {
      // Display it
      Printf("----------------------------------------------------------\n");
      Printf("Available versions (tag ROOT-vers remote-path PROOF-version):\n");
      Printf("%s", os->GetName());
      Printf("----------------------------------------------------------");
      SafeDelete(os);
   }

   // We are done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the default ROOT version to be used

Int_t TXProofMgr::SetROOTVersion(const char *tag)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Warning("SetROOTVersion","invalid TXProofMgr - do nothing");
      return -1;
   }

   // Send the request
   fSocket->SendCoordinator(kROOTVersion, tag);

   // We are done
   return (fSocket->GetOpenError() != kXR_noErrorYet) ? -1 : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Send a message to connected users. Only superusers can do this.
/// The first argument specifies the message or the file from where to take
/// the message.
/// The second argument specifies the user to which to send the message: if
/// empty or null the message is send to all the connected users.
/// return 0 in case of success, -1 in case of error

Int_t TXProofMgr::SendMsgToUsers(const char *msg, const char *usr)
{
   Int_t rc = 0;

   // Check input
   if (!msg || strlen(msg) <= 0) {
      Error("SendMsgToUsers","no message to send - do nothing");
      return -1;
   }

   // Buffer (max 32K)
   const Int_t kMAXBUF = 32768;
   char buf[kMAXBUF] = {0};
   char *p = &buf[0];
   size_t space = kMAXBUF - 1;
   Int_t lusr = 0;

   // A specific user?
   if (usr && strlen(usr) > 0 && (strlen(usr) != 1 || usr[0] != '*')) {
      lusr = (strlen(usr) + 3);
      snprintf(buf, kMAXBUF, "u:%s ", usr);
      p += lusr;
      space -= lusr;
   }

   ssize_t len = 0;
   // Is it from file ?
   if (!gSystem->AccessPathName(msg, kFileExists)) {
      // From file: can we read it ?
      if (gSystem->AccessPathName(msg, kReadPermission)) {
         Error("SendMsgToUsers","request to read message from unreadable file '%s'", msg);
         return -1;
      }
      // Open the file
      FILE *f = 0;
      if (!(f = fopen(msg, "r"))) {
         Error("SendMsgToUsers", "file '%s' cannot be open", msg);
         return -1;
      }
      // Determine the number of bytes to be read from the file.
      size_t left = 0;
      off_t rcsk = lseek(fileno(f), (off_t) 0, SEEK_END);
      if ((rcsk != (off_t)(-1))) {
         left = (size_t) rcsk;
         if ((lseek(fileno(f), (off_t) 0, SEEK_SET) == (off_t)(-1))) {
            Error("SendMsgToUsers", "cannot rewind open file (seek to 0)");
            fclose(f);
            return -1;
         }
      } else {
         Error("SendMsgToUsers", "cannot get size of open file (seek to END)");
         fclose(f);
         return -1;
      }
      // Now readout from file
      size_t wanted = left;
      if (wanted > space) {
         wanted = space;
         Warning("SendMsgToUsers",
                 "requested to send %lld bytes: max size is %lld bytes: truncating",
                 (Long64_t)left, (Long64_t)space);
      }
      do {
         while ((len = read(fileno(f), p, wanted)) < 0 &&
                  TSystem::GetErrno() == EINTR)
            TSystem::ResetErrno();
         if (len < 0) {
            SysError("SendMsgToUsers", "error reading file");
            break;
         }

         // Update counters
         left = (len >= (ssize_t)left) ? 0 : left - len;
         p += len;
         wanted = (left > kMAXBUF-1) ? kMAXBUF-1 : left;

      } while (len > 0 && left > 0);
      // Close file
      fclose(f);
   } else {
      // Add the message to the buffer
      len = strlen(msg);
      if (len > (ssize_t)space) {
         Warning("SendMsgToUsers",
                 "requested to send %lld bytes: max size is %lld bytes: truncating",
                 (Long64_t)len, (Long64_t)space);
         len = space;
      }
      memcpy(p, msg, len);
   }

   // Null-terminate
   buf[len + lusr] = 0;

   // Send the request
   fSocket->SendCoordinator(kSendMsgToUser, buf);

   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Run 'grep' on the nodes

void TXProofMgr::Grep(const char *what, const char *how, const char *where)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Error("Grep","invalid TXProofMgr - do nothing");
      return;
   }
   // Server may not support it
   if (fSocket->GetXrdProofdVersion() < 1006) {
      Error("Grep", "functionality not supported by server");
      return;
   }

   // Send the request
   TObjString *os = Exec(kGrep, what, how, where);

   // Show the result, if any
   if (os) Printf("%s", os->GetName());

   // Cleanup
   SafeDelete(os);
}

////////////////////////////////////////////////////////////////////////////////
/// Run 'find' on the nodes

void TXProofMgr::Find(const char *what, const char *how, const char *where)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Error("Find","invalid TXProofMgr - do nothing");
      return;
   }
   // Server may not support it
   if (fSocket->GetXrdProofdVersion() < 1006) {
      Error("Find", "functionality not supported by server (XrdProofd version: %d)",
                     fSocket->GetXrdProofdVersion());
      return;
   }

   // Send the request
   TObjString *os = Exec(kFind, what, how, where);

   // Show the result, if any
   if (os) Printf("%s", os->GetName());

   // Cleanup
   SafeDelete(os);
}

////////////////////////////////////////////////////////////////////////////////
/// Run 'ls' on the nodes

void TXProofMgr::Ls(const char *what, const char *how, const char *where)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Error("Ls","invalid TXProofMgr - do nothing");
      return;
   }
   // Server may not support it
   if (fSocket->GetXrdProofdVersion() < 1006) {
      Error("Ls", "functionality not supported by server");
      return;
   }

   // Send the request
   TObjString *os = Exec(kLs, what, how, where);

   // Show the result, if any
   if (os) Printf("%s", os->GetName());

   // Cleanup
   SafeDelete(os);
}

////////////////////////////////////////////////////////////////////////////////
/// Run 'more' on the nodes

void TXProofMgr::More(const char *what, const char *how, const char *where)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Error("More","invalid TXProofMgr - do nothing");
      return;
   }
   // Server may not support it
   if (fSocket->GetXrdProofdVersion() < 1006) {
      Error("More", "functionality not supported by server");
      return;
   }

   // Send the request
   TObjString *os = Exec(kMore, what, how, where);

   // Show the result, if any
   if (os) Printf("%s", os->GetName());

   // Cleanup
   SafeDelete(os);
}

////////////////////////////////////////////////////////////////////////////////
/// Run 'rm' on the nodes. The user is prompted before removal, unless 'how'
/// contains "--force" or a combination of single letter options including 'f',
/// e.g. "-fv".

Int_t TXProofMgr::Rm(const char *what, const char *how, const char *where)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Error("Rm","invalid TXProofMgr - do nothing");
      return -1;
   }
   // Server may not support it
   if (fSocket->GetXrdProofdVersion() < 1006) {
      Error("Rm", "functionality not supported by server");
      return -1;
   }

   TString prompt, ans("Y"), opt(how);
   Bool_t force = kFALSE;
   if (!opt.IsNull()) {
      TString t;
      Int_t from = 0;
      while (!force && opt.Tokenize(t, from, " ")) {
         if (t == "--force") {
            force = kTRUE;
         } else if (t.BeginsWith("-") && !t.BeginsWith("--") && t.Contains("f")) {
            force = kTRUE;
         }
      }
   }

   if (!force && isatty(0) != 0 && isatty(1) != 0) {
      // Really remove the file?
      prompt.Form("Do you really want to remove '%s'? [N/y]", what);
      ans = "";
      while (ans != "N" && ans != "Y") {
         ans = Getline(prompt.Data());
         ans.Remove(TString::kTrailing, '\n');
         if (ans == "") ans = "N";
         ans.ToUpper();
         if (ans != "N" && ans != "Y")
            Printf("Please answer y, Y, n or N");
      }
   }

   if (ans == "Y") {
      // Send the request
      TObjString *os = Exec(kRm, what, how, where);
      // Show the result, if any
      if (os) {
         if (gDebug > 1) Printf("%s", os->GetName());
         // Cleanup
         SafeDelete(os);
         // Success
         return 0;
      }
      // Failure
      return -1;
   }
   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Run 'tail' on the nodes

void TXProofMgr::Tail(const char *what, const char *how, const char *where)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Error("Tail","invalid TXProofMgr - do nothing");
      return;
   }
   // Server may not support it
   if (fSocket->GetXrdProofdVersion() < 1006) {
      Error("Tail", "functionality not supported by server");
      return;
   }

   // Send the request
   TObjString *os = Exec(kTail, what, how, where);

   // Show the result, if any
   if (os) Printf("%s", os->GetName());

   // Cleanup
   SafeDelete(os);
}

////////////////////////////////////////////////////////////////////////////////
/// Run 'md5sum' on one of the nodes

Int_t TXProofMgr::Md5sum(const char *what, TString &sum, const char *where)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Error("Md5sum","invalid TXProofMgr - do nothing");
      return -1;
   }
   // Server may not support it
   if (fSocket->GetXrdProofdVersion() < 1006) {
      Error("Md5sum", "functionality not supported by server");
      return -1;
   }

   if (where && !strcmp(where, "all")) {
      Error("Md5sum","cannot run on all nodes at once: please specify one");
      return -1;
   }

   // Send the request
   TObjString *os = Exec(kMd5sum, what, 0, where);

   // Show the result, if any
   if (os) {
      if (gDebug > 1) Printf("%s", os->GetName());
      sum = os->GetName();
      // Cleanup
      SafeDelete(os);
      // Success
      return 0;
   }
   // Failure
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Run 'stat' on one of the nodes

Int_t TXProofMgr::Stat(const char *what, FileStat_t &st, const char *where)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Error("Stat","invalid TXProofMgr - do nothing");
      return -1;
   }
   // Server may not support it
   if (fSocket->GetXrdProofdVersion() < 1006) {
      Error("Stat", "functionality not supported by server");
      return -1;
   }

   if (where && !strcmp(where, "all")) {
      Error("Stat","cannot run on all nodes at once: please specify one");
      return -1;
   }

   // Send the request
   TObjString *os = Exec(kStat, what, 0, where);

   // Show the result, if any
   if (os) {
      if (gDebug > 1) Printf("%s", os->GetName());
#if 0
      Int_t    mode, uid, gid, islink;
      Long_t   dev, ino, mtime;
      Long64_t size;
#ifdef R__WIN32
      sscanf(os->GetName(), "%ld %ld %d %d %d %I64d %ld %d", &dev, &ino, &mode,
                            &uid, &gid, &size, &mtime, &islink);
#else
      sscanf(os->GetName(), "%ld %ld %d %d %d %lld %ld %d", &dev, &ino, &mode,
                            &uid, &gid, &size, &mtime, &islink);
#endif
      if (dev == -1)
         return -1;
      st.fDev    = dev;
      st.fIno    = ino;
      st.fMode   = mode;
      st.fUid    = uid;
      st.fGid    = gid;
      st.fSize   = size;
      st.fMtime  = mtime;
      st.fIsLink = (islink == 1);
#else
      TString tkn;
      Ssiz_t from = 0;
      if (!os->GetString().Tokenize(tkn, from, "[ ]+") || !tkn.IsDigit()) return -1;
      st.fDev = tkn.Atoi();
      if (st.fDev == -1) return -1;
      if (!os->GetString().Tokenize(tkn, from, "[ ]+") || !tkn.IsDigit()) return -1;
      st.fIno = tkn.Atoi();
      if (!os->GetString().Tokenize(tkn, from, "[ ]+") || !tkn.IsDigit()) return -1;
      st.fMode = tkn.Atoi();
      if (!os->GetString().Tokenize(tkn, from, "[ ]+") || !tkn.IsDigit()) return -1;
      st.fUid = tkn.Atoi();
      if (!os->GetString().Tokenize(tkn, from, "[ ]+") || !tkn.IsDigit()) return -1;
      st.fGid = tkn.Atoi();
      if (!os->GetString().Tokenize(tkn, from, "[ ]+") || !tkn.IsDigit()) return -1;
      st.fSize = tkn.Atoll();
      if (!os->GetString().Tokenize(tkn, from, "[ ]+") || !tkn.IsDigit()) return -1;
      st.fMtime = tkn.Atoi();
      if (!os->GetString().Tokenize(tkn, from, "[ ]+") || !tkn.IsDigit()) return -1;
      st.fIsLink = (tkn.Atoi() == 1) ? kTRUE : kFALSE;
#endif

      // Cleanup
      SafeDelete(os);
      // Success
      return 0;
   }
   // Failure
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute 'action' (see EAdminExecType in 'XProofProtocol.h') at 'where'
/// (default master), with options 'how', on 'what'. The option specified by
/// 'how' are typically unix option for the relate commands. In addition to
/// the unix authorizations, the limitations are:
///
///      action = kRm        limited to the sandbox (but basic dirs cannot be
///                          removed) and on files owned by the user in the
///                          allowed directories
///      action = kTail      option '-f' is not supported and will be ignored
///

TObjString *TXProofMgr::Exec(Int_t action,
                             const char *what, const char *how, const char *where)
{
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Error("Exec","invalid TXProofMgr - do nothing");
      return (TObjString *)0;
   }
   // Server may not support it
   if (fSocket->GetXrdProofdVersion() < 1006) {
      Error("Exec", "functionality not supported by server");
      return (TObjString *)0;
   }
   // Check 'what'
   if (!what || strlen(what) <= 0) {
      Error("Exec","specifying a path is mandatory");
      return (TObjString *)0;
   }
   // Check the options
   TString opt(how);
   if (action == kTail && !opt.IsNull()) {
      // Keep only static options: -c, --bytes=N, -n , --lines=N, -N
      TString opts(how), o;
      Int_t from = 0;
      Bool_t isc = kFALSE, isn = kFALSE;
      while (opts.Tokenize(o, from, " ")) {
         // Skip values not starting with '-' is not argument to '-c' or '-n'
         if (!o.BeginsWith("-") && !isc && isn) continue;
         if (isc) {
            opt.Form("-c %s", o.Data());
            isc = kFALSE;
         }
         if (isn) {
            opt.Form("-n %s", o.Data());
            isn = kFALSE;
         }
         if (o == "-c") {
            isc = kTRUE;
         } else if (o == "-n") {
            isn = kTRUE;
         } else if (o == "--bytes=" || o == "--lines=") {
            opt = o;
         } else if (o.BeginsWith("-")) {
            o.Remove(TString::kLeading,'-');
            if (o.IsDigit()) opt.Form("-%s", o.Data());
         }
      }
   }

   // Build the command line
   TString cmd(where);
   if (cmd.IsNull()) cmd.Form("%s:%d", fUrl.GetHost(), fUrl.GetPort());
   cmd += "|";
   cmd += what;
   cmd += "|";
   cmd += opt;

   // On clients, handle Ctrl-C during collection
   if (fIntHandler) fIntHandler->Add();

   // Send the request
   TObjString *os = fSocket->SendCoordinator(kExec, cmd.Data(), action);

   // On clients, handle Ctrl-C during collection
   if (fIntHandler) fIntHandler->Remove();

   // Done
   return os;
}

////////////////////////////////////////////////////////////////////////////////
/// Get file 'remote' into 'local' from the master.
/// If opt contains "force", the file, if it exists remotely, is copied in all cases,
/// otherwise a check is done on the MD5sum.
/// If opt contains "silent" standard notificatons are not printed (errors and
/// warnings and prompts still are).
/// Return 0 on success, -1 on error.

Int_t TXProofMgr::GetFile(const char *remote, const char *local, const char *opt)
{
   Int_t rc = -1;
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Error("GetFile", "invalid TXProofMgr - do nothing");
      return rc;
   }
   // Server may not support it
   if (fSocket->GetXrdProofdVersion() < 1006) {
      Error("GetFile", "functionality not supported by server");
      return rc;
   }

   // Check remote path name
   TString filerem(remote);
   if (filerem.IsNull()) {
      Error("GetFile", "remote file path undefined");
      return rc;
   }

   // Parse option
   TString oo(opt);
   oo.ToUpper();
   Bool_t force = (oo.Contains("FORCE")) ? kTRUE : kFALSE;
   Bool_t silent = (oo.Contains("SILENT")) ? kTRUE : kFALSE;

   // Check local path name
   TString fileloc(local);
   if (fileloc.IsNull()) {
      // Set the same as the remote one, in the working dir
      fileloc = gSystem->BaseName(filerem);
   }
   gSystem->ExpandPathName(fileloc);

   // Default open and mode flags
#ifdef WIN32
   UInt_t openflags =  O_WRONLY | O_BINARY;
#else
   UInt_t openflags =  O_WRONLY;
#endif
   UInt_t openmode = 0600;

   // Get information about the local file
   UserGroup_t *ugloc = 0;
   Int_t rcloc = 0;
   FileStat_t stloc;
   if ((rcloc = gSystem->GetPathInfo(fileloc, stloc)) == 0) {
      if (R_ISDIR(stloc.fMode)) {
         // Add the filename of the remote file and re-check
         if (!fileloc.EndsWith("/")) fileloc += "/";
         fileloc += gSystem->BaseName(filerem);
         // Get again the status of the path
         rcloc = gSystem->GetPathInfo(fileloc, stloc);
      }
      if (rcloc == 0) {
         // It exists already. If it is not a regular file we cannot continue
         if (!R_ISREG(stloc.fMode)) {
            if (!silent)
               Printf("[GetFile] local file '%s' exists and is not regular: cannot continue",
                      fileloc.Data());
            return rc;
         }
         // Get our info
         if (!(ugloc = gSystem->GetUserInfo(gSystem->GetUid()))) {
            Error("GetFile", "cannot get user info for additional checks");
            return rc;
         }
         // Can we delete or overwrite it ?
         Bool_t owner = (ugloc->fUid == stloc.fUid && ugloc->fGid == stloc.fGid) ? kTRUE : kFALSE;
         Bool_t group = (!owner && ugloc->fGid == stloc.fGid) ? kTRUE : kFALSE;
         Bool_t other = (!owner && !group) ? kTRUE : kFALSE;
         delete ugloc;
         if ((owner && !(stloc.fMode & kS_IWUSR)) ||
             (group && !(stloc.fMode & kS_IWGRP)) || (other && !(stloc.fMode & kS_IWOTH))) {
            if (!silent) {
               Printf("[GetFile] file '%s' exists: no permission to delete or overwrite the file", fileloc.Data());
               Printf("[GetFile] ownership: owner: %d, group: %d, other: %d", owner, group, other);
               Printf("[GetFile] mode: %x", stloc.fMode);
            }
            return rc;
         }
         // In case we open the file, we need to truncate it
         openflags |=  O_CREAT | O_TRUNC;
      } else {
         // In case we open the file, we need to create it
         openflags |=  O_CREAT;
      }
   } else {
      // In case we open the file, we need to create it
      openflags |=  O_CREAT;
   }

   // Check the remote file exists and get it check sum
   TString remsum;
   if (Md5sum(filerem, remsum) != 0) {
      if (!silent)
         Printf("[GetFile] remote file '%s' does not exists or cannot be read", filerem.Data());
      return rc;
   }

   // If the file exists already locally, check if it is different
   bool same = 0;
   if (rcloc == 0 && !force) {
      TMD5 *md5loc = TMD5::FileChecksum(fileloc);
      if (md5loc) {
         if (remsum == md5loc->AsString()) {
            if (!silent) {
               Printf("[GetFile] local file '%s' and remote file '%s' have the same MD5 check sum",
                      fileloc.Data(), filerem.Data());
               Printf("[GetFile] use option 'force' to override");
            }
            same = 1;
         }
         delete md5loc;
      }

      // If a different file with the same name exists already, ask what to do
      if (!same) {
         const char *a = Getline("Local file exists already: would you like to overwrite it? [N/y]");
         if (a[0] == 'n' || a[0] == 'N' || a[0] == '\0') return 0;
      } else {
         return 0;
      }
   }

   // Open the local file for writing
   Int_t fdout = open(fileloc, openflags, openmode);
   if (fdout < 0) {
      Error("GetFile", "could not open local file '%s' for writing: errno: %d", local, errno);
      return rc;
   }

   // Build the command line
   TString cmd(filerem);

   // Disable TXSocket handling while receiving the file (CpProgress processes
   // pending events and this could screw-up synchronization in the TXSocket pipe)
   gSystem->RemoveFileHandler(TXSocketHandler::GetSocketHandler());

   // Send the request
   TStopwatch watch;
   watch.Start();
   TObjString *os = fSocket->SendCoordinator(kGetFile, cmd.Data());

   if (os) {
      // The message contains the size
      TString ssz(os->GetName());
      ssz.ReplaceAll(" ", "");
      if (!ssz.IsDigit()) {
         Error("GetFile", "received non-digit size string: '%s' ('%s')", os->GetName(), ssz.Data());
         close(fdout);
         return rc;
      }
      Long64_t size = ssz.Atoll();
      if (size <= 0) {
         Error("GetFile", "received null or negative size: %lld", size);
         close(fdout);
         return rc;
      }

      // Receive the file
      const Int_t kMAXBUF = 16384;  //32768  //16384  //65536;
      char buf[kMAXBUF];

      rc = 0;
      Int_t rec, r;
      Long64_t filesize = 0, left = 0;
      while (rc == 0 && filesize < size) {
         left = size - filesize;
         if (left > kMAXBUF) left = kMAXBUF;
         rec = fSocket->RecvRaw(&buf, left);
         filesize = (rec > 0) ? (filesize + rec) : filesize;
         if (rec > 0) {
            char *p = buf;
            r = rec;
            while (r) {
               Int_t w = 0;
               while ((w = write(fdout, p, r)) < 0 && TSystem::GetErrno() == EINTR)
                  TSystem::ResetErrno();
               if (w < 0) {
                  SysError("GetFile", "error writing to unit: %d", fdout);
                  rc = -1;
                  break;
               }
               r -= w;
               p += w;
            }
            // Basic progress bar
            CpProgress("GetFile", filesize, size, &watch);
         } else if (rec < 0) {
            rc = -1;
            Error("GetFile", "error during receiving file");
            break;
         }
      }
      // Finalize the progress bar
      CpProgress("GetFile", filesize, size, &watch, kTRUE);

   } else {
      Error("GetFile", "size not received");
      rc = -1;
   }

   // Restore socket handling while receiving the file
   gSystem->AddFileHandler(TXSocketHandler::GetSocketHandler());

   // Close local file
   close(fdout);
   watch.Stop();
   watch.Reset();

   if (rc == 0) {
      // Check if everything went fine
      std::unique_ptr<TMD5> md5loc(TMD5::FileChecksum(fileloc));
      if (!(md5loc.get())) {
         Error("GetFile", "cannot get MD5 checksum of the new local file '%s'", fileloc.Data());
         rc = -1;
      } else if (remsum != md5loc->AsString()) {
         Error("GetFile", "checksums for the local copy and the remote file differ: {rem:%s,loc:%s}",
                           remsum.Data(), md5loc->AsString());
         rc = -1;
      }
   }
   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Put file 'local'to 'remote' to the master
/// If opt is "force", the file, if it exists remotely, is copied in all cases,
/// otherwise a check is done on the MD5sum.
/// Return 0 on success, -1 on error

Int_t TXProofMgr::PutFile(const char *local, const char *remote, const char *opt)
{
   Int_t rc = -1;
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Error("PutFile", "invalid TXProofMgr - do nothing");
      return rc;
   }
   // Server may not support it
   if (fSocket->GetXrdProofdVersion() < 1006) {
      Error("PutFile", "functionality not supported by server");
      return rc;
   }

   // Check local path name
   TString fileloc(local);
   if (fileloc.IsNull()) {
      Error("PutFile", "local file path undefined");
      return rc;
   }
   gSystem->ExpandPathName(fileloc);

   // Parse option
   TString oo(opt);
   oo.ToUpper();
   Bool_t force = (oo == "FORCE") ? kTRUE : kFALSE;

   // Check remote path name
   TString filerem(remote);
   if (filerem.IsNull()) {
      // Set the same as the local one, in the working dir
      filerem.Form("~/%s", gSystem->BaseName(fileloc));
   } else if (filerem.EndsWith("/")) {
      // Remote path is a directory: add the file name as in the local one
      filerem += gSystem->BaseName(fileloc);
   }

   // Default open flags
#ifdef WIN32
   UInt_t openflags =  O_RDONLY | O_BINARY;
#else
   UInt_t openflags =  O_RDONLY;
#endif

   // Get information about the local file
   Int_t rcloc = 0;
   FileStat_t stloc;
   if ((rcloc = gSystem->GetPathInfo(fileloc, stloc)) != 0 || !R_ISREG(stloc.fMode)) {
      // It dies not exists or it is not a regular file: we cannot continue
      const char *why = (rcloc == 0) ? "is not regular" : "does not exists";
      Printf("[PutFile] local file '%s' %s: cannot continue", fileloc.Data(), why);
      return rc;
   }
   // Get our info
   UserGroup_t *ugloc = 0;
   if (!(ugloc = gSystem->GetUserInfo(gSystem->GetUid()))) {
      Error("PutFile", "cannot get user info for additional checks");
      return rc;
   }
   // Can we read it ?
   Bool_t owner = (ugloc->fUid == stloc.fUid && ugloc->fGid == stloc.fGid) ? kTRUE : kFALSE;
   Bool_t group = (!owner && ugloc->fGid == stloc.fGid) ? kTRUE : kFALSE;
   Bool_t other = (!owner && !group) ? kTRUE : kFALSE;
   delete ugloc;
   if ((owner && !(stloc.fMode & kS_IRUSR)) ||
       (group && !(stloc.fMode & kS_IRGRP)) || (other && !(stloc.fMode & kS_IROTH))) {
      Printf("[PutFile] file '%s': no permission to read the file", fileloc.Data());
      Printf("[PutFile] ownership: owner: %d, group: %d, other: %d", owner, group, other);
      Printf("[PutFile] mode: %x", stloc.fMode);
      return rc;
   }

   // Local MD5 sum
   TString locsum;
   TMD5 *md5loc = TMD5::FileChecksum(fileloc);
   if (!md5loc) {
      Error("PutFile", "cannot calculate the check sum for '%s'", fileloc.Data());
      return rc;
   } else {
      locsum = md5loc->AsString();
      delete md5loc;
   }

   // Check the remote file exists and get it check sum
   Bool_t same = kFALSE;
   FileStat_t strem;
   TString remsum;
   if (Stat(filerem, strem) == 0) {
      if (Md5sum(filerem, remsum) != 0) {
         Printf("[PutFile] remote file exists but the check sum calculation failed");
         return rc;
      }
      // Check sums
      if (remsum == locsum) {
         if (!force) {
            Printf("[PutFile] local file '%s' and remote file '%s' have the same MD5 check sum",
                              fileloc.Data(), filerem.Data());
            Printf("[PutFile] use option 'force' to override");
         }
         same = kTRUE;
      }
      if (!force) {
         // If a different file with the same name exists already, ask what to do
         if (!same) {
            const char *a = Getline("Remote file exists already: would you like to overwrite it? [N/y]");
            if (a[0] == 'n' || a[0] == 'N' || a[0] == '\0') return 0;
            force = kTRUE;
         } else {
            return 0;
         }
      }
   }

   // Open the local file
   int fd = open(fileloc.Data(), openflags);
   if (fd < 0) {
      Error("PutFile", "cannot open file '%s': %d", fileloc.Data(), errno);
      return -1;
   }

   // Build the command line: 'path size [opt]'
   TString cmd;
   cmd.Form("%s %lld", filerem.Data(), stloc.fSize);
   if (force) cmd += " force";

   // Disable TXSocket handling while sending the file (CpProgress processes
   // pending events and this could screw-up synchronization in the TXSocket pipe)
   gSystem->RemoveFileHandler(TXSocketHandler::GetSocketHandler());

   // Send the request
   TStopwatch watch;
   watch.Start();
   TObjString *os = fSocket->SendCoordinator(kPutFile, cmd.Data());

   if (os) {

      // Send over the file
      const Int_t kMAXBUF = 16384;  //32768  //16384  //65536;
      char buf[kMAXBUF];

      Long64_t pos = 0;
      lseek(fd, pos, SEEK_SET);

      rc = 0;
      while (rc == 0 && pos < stloc.fSize) {
         Long64_t left = stloc.fSize - pos;
         if (left > kMAXBUF) left = kMAXBUF;
         Int_t siz;
         while ((siz = read(fd, &buf[0], left)) < 0 && TSystem::GetErrno() == EINTR)
            TSystem::ResetErrno();
         if (siz < 0 || siz != left) {
            Error("PutFile", "error reading from file: errno: %d", errno);
            rc = -1;
            break;
         }
         Int_t src = 0;
         if ((src = fSocket->fConn->WriteRaw((void *)&buf[0], left)) != left) {
            Error("PutFile", "error sending over: errno: %d (rc: %d)", TSystem::GetErrno(), src);
            rc = -1;
            break;
         }
         // Basic progress bar
         CpProgress("PutFile", pos, stloc.fSize, &watch);
         // Re-position
         pos += left;
      }
      // Finalize the progress bar
      CpProgress("PutFile", pos, stloc.fSize, &watch, kTRUE);

   } else {
      Error("PutFile", "command could not be executed");
      rc = -1;
   }

   // Restore TXSocket handling
   gSystem->AddFileHandler(TXSocketHandler::GetSocketHandler());

   // Close local file
   close(fd);
   watch.Stop();
   watch.Reset();

   if (rc == 0) {
      // Check if everything went fine
      if (Md5sum(filerem, remsum) != 0) {
         Printf("[PutFile] cannot get MD5 checksum of the new remote file '%s'", filerem.Data());
         rc = -1;
      } else if (remsum != locsum) {
         Printf("[PutFile] checksums for the local copy and the remote file differ: {rem:%s, loc:%s}",
                           remsum.Data(), locsum.Data());
         rc = -1;
      }
   }

   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Print file copy progress.

void TXProofMgr::CpProgress(const char *pfx, Long64_t bytes,
                            Long64_t size, TStopwatch *watch, Bool_t cr)
{
   // Protection
   if (!pfx || size == 0 || !watch) return;

   fprintf(stderr, "[%s] Total %.02f MB\t|", pfx, (Double_t)size/1048576);

   for (int l = 0; l < 20; l++) {
      if (size > 0) {
         if (l < 20*bytes/size)
            fprintf(stderr, "=");
         else if (l == 20*bytes/size)
            fprintf(stderr, ">");
         else if (l > 20*bytes/size)
            fprintf(stderr, ".");
      } else
         fprintf(stderr, "=");
   }
   // Allow to update the GUI while uploading files
   gSystem->ProcessEvents();
   watch->Stop();
   Double_t copytime = watch->RealTime();
   fprintf(stderr, "| %.02f %% [%.01f MB/s]\r",
                   100.0*bytes/size, bytes/copytime/1048576.);
   if (cr) fprintf(stderr, "\n");
   watch->Continue();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy files in/out of the sandbox. Either 'src' or 'dst' must be in the
/// sandbox.
/// Return 0 on success, -1 on error

Int_t TXProofMgr::Cp(const char *src, const char *dst, const char *fmt)
{
   Int_t rc = -1;
   // Nothing to do if not in contact with proofserv
   if (!IsValid()) {
      Error("Cp", "invalid TXProofMgr - do nothing");
      return rc;
   }
   // Server may not support it
   if (fSocket->GetXrdProofdVersion() < 1006) {
      Error("Cp", "functionality not supported by server");
      return rc;
   }

   // Check source path name
   TString filesrc(src);
   if (filesrc.IsNull()) {
      Error("Cp", "source file path undefined");
      return rc;
   }
   // Check destination path name
   TString filedst(dst);
   if (filedst.IsNull()) {
      filedst = gSystem->BaseName(TUrl(filesrc.Data()).GetFile());
   } else if (filedst.EndsWith("/")) {
      // Remote path is a directory: add the file name as in the local one
      filedst += gSystem->BaseName(filesrc);
   }

   // Make sure that local files are in the format file://host/<file> otherwise
   // the URL class in the server will not parse them correctly
   TUrl usrc = TUrl(filesrc.Data(), kTRUE).GetUrl();
   filesrc = usrc.GetUrl();
   if (!strcmp(usrc.GetProtocol(), "file"))
      filesrc.Form("file://host/%s", usrc.GetFileAndOptions());
   TUrl udst = TUrl(filedst.Data(), kTRUE).GetUrl();
   filedst = udst.GetUrl();
   if (!strcmp(udst.GetProtocol(), "file"))
      filedst.Form("file://host/%s", udst.GetFileAndOptions());

   // Prepare the command
   TString cmd;
   cmd.Form("%s %s %s", filesrc.Data(), filedst.Data(), (fmt ? fmt : ""));

   // On clients, handle Ctrl-C during collection
   if (fIntHandler) fIntHandler->Add();

   // Send the request
   TObjString *os = fSocket->SendCoordinator(kCpFile, cmd.Data());

   // On clients, handle Ctrl-C during collection
   if (fIntHandler) fIntHandler->Remove();

   // Show the result, if any
   if (os) {
      if (gDebug > 0) Printf("%s", os->GetName());
      rc = 0;
   }

   // Done
   return rc;
}
