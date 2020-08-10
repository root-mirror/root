// @(#)root/net:$Id$
// Author: Fons Rademakers   14/08/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TNetFile                                                             //
//                                                                      //
// A TNetFile is like a normal TFile except that it reads and writes    //
// its data via a rootd server (for more on the rootd daemon see the    //
// source files root/rootd/src/*.cxx). TNetFile file names are in       //
// standard URL format with protocol "root" or "roots". The following   //
// are valid TNetFile URL's:                                            //
//                                                                      //
//    root://hpbrun.cern.ch/root/hsimple.root                           //
//    root://pcna49a:5151/~na49/data/run821.root                        //
//    root://pcna49d.cern.ch:5050//v1/data/run810.root                  //
//                                                                      //
// The only difference with the well known httpd URL's is that the root //
// of the remote file tree is the user's home directory. Therefore an   //
// absolute pathname requires a // after the host or port specifier     //
// (see last example). Further the expansion of the standard shell      //
// characters, like ~, $, .., are handled as expected.                  //
// TNetFile (actually TUrl) uses 1094 as default port for rootd.        //
//                                                                      //
// Connecting to a rootd requires the remote user id and password.      //
// TNetFile allows three ways for you to provide your login:            //
//   1) Setting it globally via the static functions:                   //
//         TAuthenticate::SetGlobalUser() and                           //
//         TAuthenticate::SetGlobalPasswd()                             //
//   2) Getting it from the ~/.netrc file (same file as used by ftp)    //
//   3) Command line prompt                                             //
// The different methods will be tried in the order given above.        //
// On machines with AFS rootd will authenticate using AFS (if it was    //
// compiled with AFS support).                                          //
//                                                                      //
// If the protocol is specified as "rootk" kerberos5 will be used for   //
// authentication.                                                      //
//                                                                      //
// The rootd daemon lives in the directory $ROOTSYS/bin. It can be      //
// started either via inetd or by hand from the command line (no need   //
// to be super user). For more info about rootd see the web page:       //
// Begin_Html <a href=http://root.cern/root/NetFile.html>NetFile</a> //
// End_Html                                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <errno.h>

#include "Bytes.h"
#include "NetErrors.h"
#include "TApplication.h"
#include "TEnv.h"
#include "TNetFile.h"
#include "TPSocket.h"
#include "TROOT.h"
#include "TSysEvtHandler.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TVirtualPerfStats.h"

// fgClientProtocol is now in TAuthenticate

ClassImp(TNetFile);
ClassImp(TNetSystem);

////////////////////////////////////////////////////////////////////////////////
/// Create a TNetFile object. This is actually done inside Create(), so
/// for a description of the options and other arguments see Create().
/// Normally a TNetFile is created via TFile::Open().

TNetFile::TNetFile(const char *url, Option_t *option, const char *ftitle,
                   Int_t compress, Int_t netopt)
   : TFile(url, "NET", ftitle, compress), fEndpointUrl(url)
{
   fSocket = 0;
   Create(url, option, netopt);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TNetFile object. To be used by derived classes, that need
/// to initialize the TFile base class but not open a connection at this
/// moment.

TNetFile::TNetFile(const char *url, const char *ftitle, Int_t compress, Bool_t)
   : TFile(url, "NET", ftitle, compress), fEndpointUrl(url)
{
   fSocket    = 0;
   fProtocol  = 0;
   fErrorCode = 0;
   fNetopt    = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// TNetFile dtor. Send close message and close socket.

TNetFile::~TNetFile()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Open a remote file. Requires fOption to be set correctly.

Int_t TNetFile::SysOpen(const char * /*file*/, Int_t /*flags*/, UInt_t /*mode*/)
{
   if (!fSocket) {

      Create(fUrl.GetUrl(), fOption, fNetopt);
      if (!fSocket) return -1;

   } else {

      if (fProtocol > 15) {
         fSocket->Send(Form("%s %s", fUrl.GetFile(), ToLower(fOption).Data()),
                       kROOTD_OPEN);
      } else {
         // Old daemon versions expect an additional slash at beginning
         fSocket->Send(Form("/%s %s", fUrl.GetFile(), ToLower(fOption).Data()),
                       kROOTD_OPEN);
      }

      EMessageTypes kind;
      int stat;
      Recv(stat, kind);

      if (kind == kROOTD_ERR) {
         PrintError("SysOpen", stat);
         return -1;
      }
   }

   // This means ok for net files
   return -2;  // set as fD in ReOpen
}

////////////////////////////////////////////////////////////////////////////////
/// Close currently open file.

Int_t TNetFile::SysClose(Int_t /*fd*/)
{
   if (fSocket)
      fSocket->Send(kROOTD_CLOSE);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return file stat information. The interface and return value is
/// identical to TSystem::GetPathInfo().

Int_t TNetFile::SysStat(Int_t, Long_t *id, Long64_t *size, Long_t *flags, Long_t *modtime)
{
   if (fProtocol < 3) return 1;

   if (!fSocket) return 1;

   fSocket->Send(kROOTD_FSTAT);

   char  msg[1024];
   Int_t kind;
   fSocket->Recv(msg, sizeof(msg), kind);

   Int_t  mode, uid, gid, islink;
   Long_t dev, ino;

   if (fProtocol > 12) {
#ifdef R__WIN32
      sscanf(msg, "%ld %ld %d %d %d %I64d %ld %d", &dev, &ino, &mode,
            &uid, &gid, size, modtime, &islink);
#else
      sscanf(msg, "%ld %ld %d %d %d %lld %ld %d", &dev, &ino, &mode,
            &uid, &gid, size, modtime, &islink);
#endif
      if (dev == -1)
         return 1;
      if (id)
         *id = (dev << 24) + ino;
      if (flags) {
         *flags = 0;
         if (mode & (kS_IXUSR|kS_IXGRP|kS_IXOTH))
            *flags |= 1;
         if (R_ISDIR(mode))
            *flags |= 2;
         if (!R_ISREG(mode) && !R_ISDIR(mode))
            *flags |= 4;
      }
   } else {
#ifdef R__WIN32
      sscanf(msg, "%ld %I64d %ld %ld", id, size, flags, modtime);
#else
      sscanf(msg, "%ld %lld %ld %ld", id, size, flags, modtime);
#endif
      if (*id == -1)
         return 1;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Close remote file.

void TNetFile::Close(Option_t *opt)
{
   if (!fSocket) return;

   TFile::Close(opt);

   if (fProtocol > 6)
      fSocket->Send(kROOTD_BYE);

   SafeDelete(fSocket);

   fD = -1;  // so TFile::IsOpen() returns false when in TFile::~TFile
}

////////////////////////////////////////////////////////////////////////////////
/// Flush file to disk.

void TNetFile::Flush()
{
   FlushWriteCache();

   if (fSocket && fWritable)
      fSocket->Send(kROOTD_FLUSH);
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize a TNetFile object.

void TNetFile::Init(Bool_t create)
{
   Seek(0);

   TFile::Init(create);
   fD = -2;   // so TFile::IsOpen() returns true when in TFile::~TFile
}

////////////////////////////////////////////////////////////////////////////////
/// Retruns kTRUE if file is open, kFALSE otherwise.

Bool_t TNetFile::IsOpen() const
{
   return fSocket == 0 ? kFALSE : kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Print some info about the net file.

void TNetFile::Print(Option_t *) const
{
   const char *fname = fUrl.GetFile();
   Printf("URL:           %s",   ((TUrl*)&fUrl)->GetUrl());
   Printf("Remote file:   %s",   &fname[1]);
   Printf("Remote user:   %s",   fUser.Data());
   Printf("Title:         %s",   fTitle.Data());
   Printf("Option:        %s",   fOption.Data());
   Printf("Bytes written: %lld", fBytesWrite);
   Printf("Bytes read:    %lld", fBytesRead);
}

////////////////////////////////////////////////////////////////////////////////
/// Print error string depending on error code.

void TNetFile::PrintError(const char *where, Int_t err)
{
   fErrorCode = err;
   Error(where, "%s", gRootdErrStr[err]);
}

////////////////////////////////////////////////////////////////////////////////
/// Reopen a file with a different access mode, like from READ to
/// UPDATE or from NEW, CREATE, RECREATE, UPDATE to READ. Thus the
/// mode argument can be either "READ" or "UPDATE". The method returns
/// 0 in case the mode was successfully modified, 1 in case the mode
/// did not change (was already as requested or wrong input arguments)
/// and -1 in case of failure, in which case the file cannot be used
/// anymore.

Int_t TNetFile::ReOpen(Option_t *mode)
{
   if (fProtocol < 7) {
      Error("ReOpen", "operation not supported by remote rootd (protocol = %d)",
            fProtocol);
      return 1;
   }

   return TFile::ReOpen(mode);
}

////////////////////////////////////////////////////////////////////////////////
/// Read specified byte range from remote file via rootd daemon.
/// Returns kTRUE in case of error.

Bool_t TNetFile::ReadBuffer(char *buf, Int_t len)
{
   if (!fSocket) return kTRUE;
   if (len == 0)
      return kFALSE;

   Bool_t result = kFALSE;

   Int_t st;
   if ((st = ReadBufferViaCache(buf, len))) {
      if (st == 2)
         return kTRUE;
      return kFALSE;
   }

   if (gApplication && gApplication->GetSignalHandler())
      gApplication->GetSignalHandler()->Delay();

   Double_t start = 0;
   if (gPerfStats) start = TTimeStamp();

   if (fSocket->Send(Form("%lld %d", fOffset, len), kROOTD_GET) < 0) {
      Error("ReadBuffer", "error sending kROOTD_GET command");
      result = kTRUE;
      goto end;
   }

   Int_t         stat, n;
   EMessageTypes kind;

   fErrorCode = -1;
   if (Recv(stat, kind) < 0 || kind == kROOTD_ERR) {
      PrintError("ReadBuffer", stat);
      result = kTRUE;
      goto end;
   }

   while ((n = fSocket->RecvRaw(buf, len)) < 0 && TSystem::GetErrno() == EINTR)
      TSystem::ResetErrno();

   if (n != len) {
      Error("ReadBuffer", "error receiving buffer of length %d, got %d", len, n);
      result = kTRUE;
      goto end;
   }

   fOffset += len;

   fBytesRead  += len;
   fReadCalls++;
#ifdef R__WIN32
   SetFileBytesRead(GetFileBytesRead() + len);
   SetFileReadCalls(GetFileReadCalls() + 1);
#else
   fgBytesRead += len;
   fgReadCalls++;
#endif

end:

   if (gPerfStats)
      gPerfStats->FileReadEvent(this, len, start);

   if (gApplication && gApplication->GetSignalHandler())
      gApplication->GetSignalHandler()->HandleDelayedSignal();

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Read specified byte range from remote file via rootd daemon.
/// Returns kTRUE in case of error.

Bool_t TNetFile::ReadBuffer(char *buf, Long64_t pos, Int_t len)
{
   SetOffset(pos);
   return ReadBuffer(buf, len);
}

////////////////////////////////////////////////////////////////////////////////
/// Read a list of buffers given in pos[] and len[] and return it in a single
/// buffer.
/// Returns kTRUE in case of error.

Bool_t TNetFile::ReadBuffers(char *buf,  Long64_t *pos, Int_t *len, Int_t nbuf)
{
   if (!fSocket) return kTRUE;

   // If it's an old version of the protocol try the default TFile::ReadBuffers
   if (fProtocol < 17)
      return TFile::ReadBuffers(buf, pos, len, nbuf);

   Int_t   stat;
   Int_t   blockSize = 262144;  //Let's say we transfer 256KB at the time
   Bool_t  result = kFALSE;
   EMessageTypes kind;
   TString data_buf;      // buf to put the info

   if (gApplication && gApplication->GetSignalHandler())
      gApplication->GetSignalHandler()->Delay();

   Double_t start = 0;
   if (gPerfStats) start = TTimeStamp();

   // Make the string with a list of offsets and lengths
   Long64_t total_len = 0;
   Long64_t actual_pos;
   for(Int_t i = 0; i < nbuf; i++) {
      data_buf += pos[i] + fArchiveOffset;
      data_buf += "-";
      data_buf += len[i];
      data_buf += "/";
      total_len += len[i];
   }

   // Send the command with the length of the info and number of buffers
   if (fSocket->Send(Form("%d %d %d", nbuf, data_buf.Length(), blockSize),
                          kROOTD_GETS) < 0) {
      Error("ReadBuffers", "error sending kROOTD_GETS command");
      result = kTRUE;
      goto end;
   }
   // Send buffer with the list of offsets and lengths
   if (fSocket->SendRaw(data_buf, data_buf.Length()) < 0) {
      Error("ReadBuffers", "error sending buffer");
      result = kTRUE;
      goto end;
   }

   fErrorCode = -1;
   if (Recv(stat, kind) < 0 || kind == kROOTD_ERR) {
      PrintError("ReadBuffers", stat);
      result = kTRUE;
      goto end;
   }

   actual_pos = 0;
   while (actual_pos < total_len) {
      Long64_t left = total_len - actual_pos;
      if (left > blockSize)
         left = blockSize;

      Int_t n;
      while ((n = fSocket->RecvRaw(buf + actual_pos, Int_t(left))) < 0 &&
             TSystem::GetErrno() == EINTR)
         TSystem::ResetErrno();

      if (n != Int_t(left)) {
         Error("GetBuffers", "error receiving buffer of length %d, got %d",
               Int_t(left), n);
         result = kTRUE ;
         goto end;
      }
      actual_pos += left;
   }

   fBytesRead  += total_len;
   fReadCalls++;
#ifdef R__WIN32
   SetFileBytesRead(GetFileBytesRead() + total_len);
   SetFileReadCalls(GetFileReadCalls() + 1);
#else
   fgBytesRead += total_len;
   fgReadCalls++;
#endif

end:

   if (gPerfStats)
      gPerfStats->FileReadEvent(this, total_len, start);

   if (gApplication && gApplication->GetSignalHandler())
      gApplication->GetSignalHandler()->HandleDelayedSignal();

   // If found problems try the generic implementation
   if (result) {
      if (gDebug > 0)
         Info("ReadBuffers", "Couldnt use the specific implementation, calling TFile::ReadBuffers");
      return TFile::ReadBuffers(buf, pos, len, nbuf);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Write specified byte range to remote file via rootd daemon.
/// Returns kTRUE in case of error.

Bool_t TNetFile::WriteBuffer(const char *buf, Int_t len)
{
   if (!fSocket || !fWritable) return kTRUE;

   Bool_t result = kFALSE;

   Int_t st;
   if ((st = WriteBufferViaCache(buf, len))) {
      if (st == 2)
         return kTRUE;
      return kFALSE;
   }

   gSystem->IgnoreInterrupt();

   if (fSocket->Send(Form("%lld %d", fOffset, len), kROOTD_PUT) < 0) {
      SetBit(kWriteError);
      Error("WriteBuffer", "error sending kROOTD_PUT command");
      result = kTRUE;
      goto end;
   }
   if (fSocket->SendRaw(buf, len) < 0) {
      SetBit(kWriteError);
      Error("WriteBuffer", "error sending buffer");
      result = kTRUE;
      goto end;
   }

   Int_t         stat;
   EMessageTypes kind;

   fErrorCode = -1;
   if (Recv(stat, kind) < 0 || kind == kROOTD_ERR) {
      SetBit(kWriteError);
      PrintError("WriteBuffer", stat);
      result = kTRUE;
      goto end;
   }

   fOffset += len;

   fBytesWrite  += len;
#ifdef R__WIN32
   SetFileBytesWritten(GetFileBytesWritten() + len);
#else
   fgBytesWrite += len;
#endif

end:
   gSystem->IgnoreInterrupt(kFALSE);

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return status from rootd server and message kind. Returns -1 in
/// case of error otherwise 8 (sizeof 2 words, status and kind).

Int_t TNetFile::Recv(Int_t &status, EMessageTypes &kind)
{
   kind   = kROOTD_ERR;
   status = 0;

   if (!fSocket) return -1;

   Int_t what;
   Int_t n = fSocket->Recv(status, what);
   kind = (EMessageTypes) what;
   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Set position from where to start reading.

void TNetFile::Seek(Long64_t offset, ERelativeTo pos)
{
   SetOffset(offset, pos);
}

////////////////////////////////////////////////////////////////////////////////
/// Connect to remote rootd server.

void TNetFile::ConnectServer(Int_t *stat, EMessageTypes *kind, Int_t netopt,
                             Int_t tcpwindowsize, Bool_t forceOpen,
                             Bool_t forceRead)
{
   TString fn = fUrl.GetFile();

   // Create Authenticated socket
   Int_t sSize = netopt < -1 ? -netopt : 1;
   TString url(fUrl.GetProtocol());
   if (url.Contains("root")) {
      url.Insert(4,"dp");
   } else {
      url = "rootdp";
   }
   url += TString(Form("://%s@%s:%d",
                       fUrl.GetUser(), fUrl.GetHost(), fUrl.GetPort()));
   fSocket = TSocket::CreateAuthSocket(url, sSize, tcpwindowsize, fSocket, stat);
   if (!fSocket || (fSocket && !fSocket->IsAuthenticated())) {
      if (sSize > 1)
         Error("TNetFile", "can't open %d-stream connection to rootd on "
               "host %s at port %d", sSize, fUrl.GetHost(), fUrl.GetPort());
      else
         Error("TNetFile", "can't open connection to rootd on "
               "host %s at port %d", fUrl.GetHost(), fUrl.GetPort());
      *kind = kROOTD_ERR;
      goto zombie;
   }

   // Check if rootd supports new options
   fProtocol = fSocket->GetRemoteProtocol();
   if (forceRead && fProtocol < 5) {
      Warning("ConnectServer", "rootd does not support \"+read\" option");
      forceRead = kFALSE;
   }

   // Open the file
   if (fProtocol < 16)
      // For backward compatibility we need to add a '/' at the beginning
      fn.Insert(0,"/");
   if (forceOpen)
      fSocket->Send(Form("%s %s", fn.Data(),
                              ToLower("f"+fOption).Data()), kROOTD_OPEN);
   else if (forceRead)
      fSocket->Send(Form("%s %s", fn.Data(), "+read"), kROOTD_OPEN);
   else
      fSocket->Send(Form("%s %s", fn.Data(),
                              ToLower(fOption).Data()), kROOTD_OPEN);

   EMessageTypes tmpkind;
   int  tmpstat;
   Recv(tmpstat, tmpkind);
   *stat = tmpstat;
   *kind = tmpkind;

   return;

zombie:
   // error in file opening occured, make this object a zombie
   MakeZombie();
   SafeDelete(fSocket);
   gDirectory = gROOT;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a NetFile object. A net file is the same as a TFile
/// except that it is being accessed via a rootd server. The url
/// argument must be of the form: root[k]://host.dom.ain/file.root.
/// When protocol is "rootk" try using kerberos5 authentication.
/// If the file specified in the URL does not exist, is not accessable
/// or can not be created the kZombie bit will be set in the TNetFile
/// object. Use IsZombie() to see if the file is accessable.
/// If the remote daemon thinks the file is still connected, while you are
/// sure this is not the case you can force open the file by preceding the
/// option argument with an "-", e.g.: "-recreate". Do this only
/// in cases when you are very sure nobody else is using the file.
/// To bypass the writelock on a file, to allow the reading of a file
/// that is being written by another process, explicitly specify the
/// "+read" option ("read" being the default option).
/// The netopt argument can be used to specify the size of the tcp window in
/// bytes (for more info see: http://www.psc.edu/networking/perf_tune.html).
/// The default and minimum tcp window size is 65535 bytes.
/// If netopt < -1 then |netopt| is the number of parallel sockets that will
/// be used to connect to rootd. This option should be used on fat pipes
/// (i.e. high bandwidth, high latency links). The ideal number of parallel
/// sockets depends on the bandwidth*delay product. Generally 5-7 is a good
/// number.
/// For a description of the option and other arguments see the TFile ctor.
/// The preferred interface to this constructor is via TFile::Open().

void TNetFile::Create(const char * /*url*/, Option_t *option, Int_t netopt)
{
   Int_t tcpwindowsize = 65535;

   fErrorCode = -1;
   fNetopt    = netopt;
   fOption    = option;

   Bool_t forceOpen = kFALSE;
   if (option[0] == '-') {
      fOption   = &option[1];
      forceOpen = kTRUE;
   }
   // accept 'f', like 'frecreate' still for backward compatibility
   if (option[0] == 'F' || option[0] == 'f') {
      fOption   = &option[1];
      forceOpen = kTRUE;
   }

   Bool_t forceRead = kFALSE;
   if (!strcasecmp(option, "+read")) {
      fOption   = &option[1];
      forceRead = kTRUE;
   }

   fOption.ToUpper();

   if (fOption == "NEW")
      fOption = "CREATE";

   Bool_t create   = (fOption == "CREATE") ? kTRUE : kFALSE;
   Bool_t recreate = (fOption == "RECREATE") ? kTRUE : kFALSE;
   Bool_t update   = (fOption == "UPDATE") ? kTRUE : kFALSE;
   Bool_t read     = (fOption == "READ") ? kTRUE : kFALSE;
   if (!create && !recreate && !update && !read) {
      read    = kTRUE;
      fOption = "READ";
   }

   if (!fUrl.IsValid()) {
      Error("Create", "invalid URL specified: %s", fUrl.GetUrl());
      goto zombie;
   }

   if (netopt > tcpwindowsize)
      tcpwindowsize = netopt;

   // Open connection to remote rootd server
   EMessageTypes kind;
   Int_t stat;
   ConnectServer(&stat, &kind, netopt, tcpwindowsize, forceOpen, forceRead);
   if (gDebug > 2) Info("Create", "got from host %d %d", stat, kind);

   if (kind == kROOTD_ERR) {
      PrintError("Create", stat);
      Error("Create", "failing on file %s", fUrl.GetUrl());
      goto zombie;
   }

   if (recreate) {
      recreate = kFALSE;
      create   = kTRUE;
      fOption  = "CREATE";
   }

   if (update && stat > 1) {
      update = kFALSE;
      create = kTRUE;
      stat   = 1;
   }

   if (stat == 1)
      fWritable = kTRUE;
   else
      fWritable = kFALSE;

   Init(create);
   return;

zombie:
   // error in file opening occured, make this object a zombie
   MakeZombie();
   SafeDelete(fSocket);
   gDirectory = gROOT;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a NetFile object using an existing connection (socket s).
/// Provided for use in TXNetFile.
/// See:
///    TNetFile::Create(const char *url, Option_t *option, Int_t netopt)
/// for details about the arguments.

void TNetFile::Create(TSocket *s, Option_t *option, Int_t netopt)
{
   // Import socket
   fSocket = s;

   // Create the connection
   Create(s->GetUrl(), option, netopt);
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if 'url' matches the coordinates of this file.
/// Check the full URL, including port and FQDN.

Bool_t TNetFile::Matches(const char *url)
{
   // Run standard check on fUrl, first
   Bool_t rc = TFile::Matches(url);
   if (rc)
      // Ok, standard check enough
      return kTRUE;

   // Check also endpoint URL
   TUrl u(url);
   if (!strcmp(u.GetFile(),fEndpointUrl.GetFile())) {
      // Candidate info
      TString fqdn = u.GetHostFQDN();

      // Check ports
      if (u.GetPort() == fEndpointUrl.GetPort()) {
         TString fqdnref = fEndpointUrl.GetHostFQDN();
         if (fqdn == fqdnref)
            // Ok, coordinates match
            return kTRUE;
      }
   }

   // Default is not matching
   return kFALSE;
}

//
// TNetSystem: the directory handler for net files
//

////////////////////////////////////////////////////////////////////////////////
/// Create helper class that allows directory access via rootd.
/// Use ftpowner = TRUE (default) if this instance is responsible
/// for cleaning of the underlying TFTP connection; this allows
/// to have control on the order of the final cleaning.

TNetSystem::TNetSystem(Bool_t ftpowner)
           : TSystem("-root", "Net file Helper System")
{
   // name must start with '-' to bypass the TSystem singleton check
   SetName("root");

   fDir = kFALSE;
   fDirp = 0;
   fFTP  = 0;
   fFTPOwner = ftpowner;
   fUser = "";
   fHost = "";
   fPort = -1;
   fIsLocal = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create helper class that allows directory access via rootd.
/// Use ftpowner = TRUE (default) if this instance is responsible
/// for cleaning of the underlying TFTP connection; this allows
/// to have control on the order of the final cleaning.

TNetSystem::TNetSystem(const char *url, Bool_t ftpowner)
           : TSystem("-root", "Net file Helper System")
{
   // name must start with '-' to bypass the TSystem singleton check
   SetName("root");

   fFTPOwner = ftpowner;
   fIsLocal = kFALSE;
   Create(url);
}

////////////////////////////////////////////////////////////////////////////////
/// Parse and save coordinates of the remote entity (user, host, port, ...)

void TNetSystem::InitRemoteEntity(const char *url)
{
   TUrl turl(url);

   // Remote username: local as default
   fUser = turl.GetUser();
   if (!fUser.Length()) {
      UserGroup_t *u = gSystem->GetUserInfo();
      if (u)
         fUser = u->fUser;
      delete u;
   }

   // Check and save the host FQDN ...
   fHost = turl.GetHostFQDN();

   // Remote port: the default should be 1094 because we are here
   // only if the protocol is "root://"
   fPort = turl.GetPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TNetSystem object.

void TNetSystem::Create(const char *url, TSocket *sock)
{
   // If we got here protocol must be at least its short form "^root.*:" :
   // make sure that it is in the full form to avoid problems in TFTP
   TString surl(url);
   if (!surl.Contains("://")) {
      surl.Insert(surl.Index(":")+1,"//");
   }
   TUrl turl(surl);

   fDir  = kFALSE;
   fDirp = 0;
   fFTP  = 0;

   // Check locality, taking into account possible prefixes
   fLocalPrefix = "";
   fIsLocal = kFALSE;
   // We may have been asked explicitly to go through the daemon
   Bool_t forceRemote = gEnv->GetValue("Path.ForceRemote", 0);
   TString opts = TUrl(url).GetOptions();
   if (opts.Contains("remote=1"))
      forceRemote = kTRUE;
   else if (opts.Contains("remote=0"))
      forceRemote = kFALSE;
   if (!forceRemote) {
      if ((fIsLocal = TSystem::IsPathLocal(url))) {
         fLocalPrefix = gEnv->GetValue("Path.Localroot","");
         // Nothing more to do
         return;
      }
   }

   // Fill in user, host, port
   InitRemoteEntity(surl);

   // Build a TFTP url
   if (fHost.Length()) {
      TString eurl = "";
      // First the protocol
      if (strlen(turl.GetProtocol())) {
         eurl = turl.GetProtocol();
         eurl += "://";
      } else
         eurl = "root://";
      // Add user, if specified
      if (strlen(turl.GetUser())) {
         eurl += turl.GetUser();
         eurl += "@";
      }
      // Add host
      eurl += fHost;
      // Add port
      eurl += ":";
      eurl += turl.GetPort();

      fFTP  = new TFTP(eurl, 1, TFTP::kDfltWindowSize, sock);
      if (fFTP && fFTP->IsOpen()) {
         if (fFTP->GetSocket()->GetRemoteProtocol() < 12) {
            Error("Create",
                  "remote daemon does not support 'system' functionality");
            fFTP->Close();
            delete fFTP;
         } else {
            fUser = fFTP->GetSocket()->GetSecContext()->GetUser();
            fHost = fFTP->GetSocket()->GetSecContext()->GetHost();
            // If responsible for the TFTP connection, remove it from the
            // socket global list to avoid problems with double deletion
            // at final cleanup
            if (fFTPOwner)
               gROOT->GetListOfSockets()->Remove(fFTP);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TNetSystem::~TNetSystem()
{
   // Close FTP connection
   if (fFTPOwner) {
      if (fFTP) {
         if (fFTP->IsOpen()) {

            // Close remote directory if still open
            if (fDir) {
               fFTP->FreeDirectory(kFALSE);
               fDir = kFALSE;
            }
            fFTP->Close();
         }
         delete fFTP;
      }
   }
   fDirp = 0;
   fFTP  = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a directory via rootd.

Int_t TNetSystem::MakeDirectory(const char *dir)
{
   // If local, use the local TSystem
   if (fIsLocal) {
      TString edir = TUrl(dir).GetFile();
      if (!fLocalPrefix.IsNull())
         edir.Insert(0, fLocalPrefix);
      return gSystem->MakeDirectory(edir);
   }

   if (fFTP && fFTP->IsOpen()) {
      // Extract the directory name
      TString edir = TUrl(dir).GetFile();
      return fFTP->MakeDirectory(edir,kFALSE);
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a directory and return an opaque pointer to a dir structure.
/// Returns nullptr in case of error.

void *TNetSystem::OpenDirectory(const char *dir)
{
   // If local, use the local TSystem
   if (fIsLocal) {
      TString edir = TUrl(dir).GetFile();
      if (!fLocalPrefix.IsNull())
         edir.Insert(0, fLocalPrefix);
      return gSystem->OpenDirectory(edir);
   }

   if (!fFTP || !fFTP->IsOpen())
      return nullptr;

   if (fDir) {
      if (gDebug > 0)
         Info("OpenDirectory", "a directory is already open: close it first");
      fFTP->FreeDirectory(kFALSE);
      fDir = kFALSE;
   }

   // Extract the directory name
   TString edir = TUrl(dir).GetFile();

   if (fFTP->OpenDirectory(edir,kFALSE)) {
      fDir = kTRUE;
      fDirp = (void *)&fDir;
      return fDirp;
   } else
      return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Free directory via rootd.

void TNetSystem::FreeDirectory(void *dirp)
{
   // If local, use the local TSystem
   if (fIsLocal) {
      gSystem->FreeDirectory(dirp);
      return;
   }

   if (dirp != fDirp) {
      Error("FreeDirectory", "invalid directory pointer (should never happen)");
      return;
   }

   if (fFTP && fFTP->IsOpen()) {
      if (fDir) {
         fFTP->FreeDirectory(kFALSE);
         fDir = kFALSE;
         fDirp = 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get directory entry via rootd. Returns 0 in case no more entries.

const char *TNetSystem::GetDirEntry(void *dirp)
{
   // If local, use the local TSystem
   if (fIsLocal) {
      return gSystem->GetDirEntry(dirp);
   }

   if (dirp != fDirp) {
      Error("GetDirEntry", "invalid directory pointer (should never happen)");
      return 0;
   }

   if (fFTP && fFTP->IsOpen() && fDir) {
      return fFTP->GetDirEntry(kFALSE);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get info about a file. Info is returned in the form of a FileStat_t
/// structure (see TSystem.h).
/// The function returns 0 in case of success and 1 if the file could
/// not be stat'ed.

Int_t TNetSystem::GetPathInfo(const char *path, FileStat_t &buf)
{
   // If local, use the local TSystem
   if (fIsLocal) {
      TString epath = TUrl(path).GetFile();
      if (!fLocalPrefix.IsNull())
         epath.Insert(0, fLocalPrefix);
      return gSystem->GetPathInfo(epath, buf);
   }

   if (fFTP && fFTP->IsOpen()) {
      // Extract the directory name
      TString epath = TUrl(path).GetFile();
      fFTP->GetPathInfo(epath, buf, kFALSE);
      return 0;
   }
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns FALSE if one can access a file using the specified access mode.
/// Mode is the same as for the Unix access(2) function.
/// Attention, bizarre convention of return value!!

Bool_t TNetSystem::AccessPathName(const char *path, EAccessMode mode)
{
   // If local, use the local TSystem
   if (fIsLocal) {
      TString epath = TUrl(path).GetFile();
      if (!fLocalPrefix.IsNull())
         epath.Insert(0, fLocalPrefix);
      return gSystem->AccessPathName(epath, mode);
   }

   if (fFTP && fFTP->IsOpen()) {
      // Extract the directory name
      TString epath = TUrl(path).GetFile();
      return fFTP->AccessPathName(epath, mode, kFALSE);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check consistency of this helper with the one required
/// by 'path' or 'dirptr'.

Bool_t TNetSystem::ConsistentWith(const char *path, void *dirptr)
{
   // Standard check: only the protocol part of 'path' is required to match
   Bool_t checkstd = TSystem::ConsistentWith(path, dirptr);
   if (!checkstd) return kFALSE;

   // Require match of 'user' and 'host'
   Bool_t checknet = path ? kFALSE : kTRUE;
   if (path && strlen(path)) {

      // Get user name
      TUrl url(path);
      TString user = url.GetUser();
      if (user.IsNull() && !fUser.IsNull()) {
         UserGroup_t *u = gSystem->GetUserInfo();
         if (u)
            user = u->fUser;
         delete u;
      }

      // Get host name
      TString host = url.GetHostFQDN();

      // Get port
      Int_t port = url.GetPort();
      if (gDebug > 1)
         Info("ConsistentWith", "fUser:'%s' (%s), fHost:'%s' (%s), fPort:%d (%d)",
                                fUser.Data(), user.Data(), fHost.Data(), host.Data(),
                                fPort, port);

      if (user == fUser && host == fHost && port == fPort)
         checknet = kTRUE;
   }

   return (checkstd && checknet);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a path

Int_t TNetSystem::Unlink(const char *path)
{
   // If local, use the local TSystem
   if (fIsLocal) {
      TString epath = TUrl(path).GetFile();
      if (!fLocalPrefix.IsNull())
         epath.Insert(0, fLocalPrefix);
      return gSystem->Unlink(epath);
   }

   // Not implemented for rootd
   Warning("Unlink", "functionality not implemented - ignored (path: %s)", path);
   return -1;
}
