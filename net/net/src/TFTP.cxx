// @(#)root/net:$Id$
// Author: Fons Rademakers   13/02/2001

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFTP                                                                 //
//                                                                      //
// This class provides all infrastructure for a performant file         //
// transfer protocol. It works in conjuction with the rootd daemon      //
// and can use parallel sockets to improve performance over fat pipes.  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <ROOT/RConfig.hxx>

#include <fcntl.h>
#include <errno.h>
#include <sys/stat.h>
#ifndef R__WIN32
#   include <unistd.h>
#else
#   define ssize_t int
#   include <io.h>
#   include <sys/types.h>
#endif

#include "strlcpy.h"
#include "TFTP.h"
#include "TPSocket.h"
#include "TUrl.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TError.h"
#include "NetErrors.h"
#include "TRegexp.h"
#include "TVirtualMutex.h"

#if defined(R__UNIX) || defined(R__MACOSX)
#define HAVE_MMAP
#endif

#ifdef HAVE_MMAP
#   include <sys/mman.h>
#ifndef MAP_FILE
#define MAP_FILE 0           /* compatability flag */
#endif
#endif


Long64_t TFTP::fgBytesWrite = 0;
Long64_t TFTP::fgBytesRead  = 0;


ClassImp(TFTP);

////////////////////////////////////////////////////////////////////////////////
/// Open connection to host specified by the url using par parallel sockets.
/// The url has the form: [root[s,k]://]host[:port].
/// If port is not specified the default rootd port (1094) will be used.
/// Using wsize one can specify the tcp window size. Normally this is not
/// needed when using parallel sockets.
/// An existing connection (TSocket *sock) can also be used to establish
/// the FTP session.

TFTP::TFTP(const char *url, Int_t par, Int_t wsize, TSocket *sock)
{
   fSocket = sock;

   TString s = url;
   if (s.Contains("://")) {
      if (!s.BeginsWith("root")) {
         Error("TFTP",
               "url must be of the form \"[root[up,s,k,g,h,ug]://]host[:port]\"");
         MakeZombie();
         return;
      }
   } else
      s = "root://" + s;

   Init(s, par, wsize);
}

////////////////////////////////////////////////////////////////////////////////
/// Set up the actual connection.

void TFTP::Init(const char *surl, Int_t par, Int_t wsize)
{
   TUrl url(surl);
   TString hurl(url.GetProtocol());
   if (hurl.Contains("root")) {
      hurl.Insert(4,"dp");
   } else {
      hurl = "rootdp";
   }
   hurl += TString(Form("://%s@%s:%d",
                        url.GetUser(), url.GetHost(), url.GetPort()));
   fSocket = TSocket::CreateAuthSocket(hurl, par, wsize, fSocket);
   if (!fSocket || !fSocket->IsAuthenticated()) {
      if (par > 1)
         Error("TFTP", "can't open %d-stream connection to rootd on "
               "host %s at port %d", par, url.GetHost(), url.GetPort());
      else
         Error("TFTP", "can't open connection to rootd on "
               "host %s at port %d", url.GetHost(), url.GetPort());
      goto zombie;
   }

   fProtocol = fSocket->GetRemoteProtocol();
   fUser = fSocket->GetSecContext()->GetUser();

   fHost       = url.GetHost();
   fPort       = url.GetPort();
   fParallel   = par;
   fWindowSize = wsize;
   fLastBlock  = 0;
   fRestartAt  = 0;
   fBlockSize  = kDfltBlockSize;
   fMode       = kBinary;
   fBytesWrite = 0;
   fBytesRead  = 0;

   // Replace our socket in the list with this
   // for consistency during the final cleanup
   // (The socket will be delete by us when everything is ok remotely)
   {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfSockets()->Remove(fSocket);
      gROOT->GetListOfSockets()->Add(this);
   }
   return;

zombie:
   MakeZombie();
   SafeDelete(fSocket);
}

////////////////////////////////////////////////////////////////////////////////
/// TFTP dtor. Send close message and close socket.

TFTP::~TFTP()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Print some info about the FTP connection.

void TFTP::Print(Option_t *) const
{
   TString secCont;

   Printf("Local host:           %s", gSystem->HostName());
   Printf("Remote host:          %s [%d]", fHost.Data(), fPort);
   Printf("Remote user:          %s", fUser.Data());
   if (fSocket->IsAuthenticated())
      Printf("Security context:     %s",
                                      fSocket->GetSecContext()->AsString(secCont));
   Printf("Rootd protocol vers.: %d", fSocket->GetRemoteProtocol());
   if (fParallel > 1) {
      Printf("Parallel sockets:     %d", fParallel);
   }
   Printf("TCP window size:      %d",   fWindowSize);
   Printf("Rootd protocol:       %d",   fProtocol);
   Printf("Transfer block size:  %d",   fBlockSize);
   Printf("Transfer mode:        %s",   fMode ? "ascii" : "binary");
   Printf("Bytes sent:           %lld", fBytesWrite);
   Printf("Bytes received:       %lld", fBytesRead);
}

////////////////////////////////////////////////////////////////////////////////
/// Print error string depending on error code.

void TFTP::PrintError(const char *where, Int_t err) const
{
   Error(where, "%s", gRootdErrStr[err]);
}

////////////////////////////////////////////////////////////////////////////////
/// Return status from rootd server and message kind. Returns -1 in
/// case of error otherwise 8 (sizeof 2 words, status and kind).

Int_t TFTP::Recv(Int_t &status, EMessageTypes &kind) const
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
/// Make sure the block size is a power of two, with a minimum of 32768.

void TFTP::SetBlockSize(Int_t blockSize)
{
   if (blockSize < 32768) {
      fBlockSize = 32768;
      return;
   }

   int i;
   for (i = 0; i < int(sizeof(blockSize)*8); i++)
      if ((blockSize >> i) == 1)
         break;

   fBlockSize = 1 << i;
}

////////////////////////////////////////////////////////////////////////////////
/// Transfer file to remote host. Returns number of bytes
/// sent or < 0 in case of error. Error -1 connection is still
/// open, error -2 connection has been closed. In case of failure
/// fRestartAt is set to the number of bytes correclty transfered.
/// Calling PutFile() immediately afterwards will restart at fRestartAt.
/// If this is not desired call SetRestartAt(0) before calling PutFile().
/// If rootd reports that the file is locked, and you are sure this is not
/// the case (e.g. due to a crash), you can force unlock it by prepending
/// the remoteName with a '-'.

Long64_t TFTP::PutFile(const char *file, const char *remoteName)
{
   if (!IsOpen() || !file || !*file) return -1;

#if defined(R__WIN32) || defined(R__WINGCC)
   Int_t fd = open(file, O_RDONLY | O_BINARY);
#elif defined(R__SEEK64)
   Int_t fd = open64(file, O_RDONLY);
#else
   Int_t fd = open(file, O_RDONLY);
#endif
   if (fd < 0) {
      Error("PutFile", "cannot open %s in read mode", file);
      return -1;
   }

   Long64_t size;
   Long_t id, flags, modtime;
   if (gSystem->GetPathInfo(file, &id, &size, &flags, &modtime) == 0) {
      if (flags > 1) {
         Error("PutFile", "%s not a regular file (%ld)", file, flags);
         close(fd);
         return -1;
      }
   } else {
      Warning("PutFile", "could not stat %s", file);
      close(fd);
      return -1;
   }

   if (!remoteName)
      remoteName = file;

   Long64_t restartat = fRestartAt;

   // check if restartat value makes sense
   if (restartat && (restartat >= size))
      restartat = 0;

   if (fSocket->Send(Form("%s %d %d %lld %lld", remoteName, fBlockSize, fMode,
                     size, restartat), kROOTD_PUTFILE) < 0) {
      Error("PutFile", "error sending kROOTD_PUTFILE command");
      close(fd);
      return -2;
   }

   Int_t         stat;
   EMessageTypes kind;

   if (Recv(stat, kind) < 0 || kind == kROOTD_ERR) {
      PrintError("PutFile", stat);
      close(fd);
      return -1;
   }

   Info("PutFile", "sending file %s (%lld bytes, starting at %lld)",
        file, size, restartat);

   TStopwatch timer;
   timer.Start();

   Long64_t pos = restartat & ~(fBlockSize-1);
   Int_t skip = restartat - pos;

#ifndef HAVE_MMAP
   char *buf = new char[fBlockSize];
#if defined(R__SEEK64)
   lseek64(fd, pos, SEEK_SET);
#elif defined(R__WIN32)
   _lseeki64(fd, pos, SEEK_SET);
#else
   lseek(fd, pos, SEEK_SET);
#endif
#endif

   while (pos < size) {
      Long64_t left = Long64_t(size - pos);
      if (left > fBlockSize)
         left = fBlockSize;
#ifdef HAVE_MMAP
#if defined(R__SEEK64)
      char *buf = (char*) mmap64(0, left, PROT_READ, MAP_FILE | MAP_SHARED, fd, pos);
#else
      char *buf = (char*) mmap(0, left, PROT_READ, MAP_FILE | MAP_SHARED, fd, pos);
#endif
      if (buf == (char *) -1) {
         Error("PutFile", "mmap of file %s failed", file);
         close(fd);
         return -1;
      }
#else
      Int_t siz;
      while ((siz = read(fd, buf, left)) < 0 && TSystem::GetErrno() == EINTR)
         TSystem::ResetErrno();
      if (siz < 0 || siz != left) {
         Error("PutFile", "error reading from file %s", file);
         // Send urgent message to rootd to stop tranfer
         delete [] buf;
         close(fd);
         return -1;
      }
#endif

      if (fSocket->SendRaw(buf+skip, left-skip) < 0) {
         Error("PutFile", "error sending buffer");
         // Send urgent message to rootd to stop transfer
#ifdef HAVE_MMAP
         munmap(buf, left);
#else
         delete [] buf;
#endif
         close(fd);
         return -2;
      }

      fBytesWrite  += left-skip;
      fgBytesWrite += left-skip;

      fRestartAt = pos;   // bytes correctly sent up till now

      pos += left;
      skip = 0;

#ifdef HAVE_MMAP
      munmap(buf, left);
#endif
   }

#ifndef HAVE_MMAP
   delete [] buf;
#endif

   close(fd);

   fRestartAt = 0;

   // get acknowlegdement from server that file was stored correctly
   if (Recv(stat, kind) < 0 || kind == kROOTD_ERR) {
      PrintError("PutFile", stat);
      close(fd);
      return -1;
   }

   // provide timing numbers
   Double_t speed, t = timer.RealTime();
   if (t > 0)
      speed = Double_t(size - restartat) / t;
   else
      speed = 0.0;
   if (speed > 524288)
      Info("PutFile", "%.3f seconds, %.2f Mbytes per second",
           t, speed / 1048576);
   else if (speed > 512)
      Info("PutFile", "%.3f seconds, %.2f Kbytes per second",
           t, speed / 1024);
   else
      Info("PutFile", "%.3f seconds, %.2f bytes per second",
           t, speed);

   return Long64_t(size - restartat);
}

////////////////////////////////////////////////////////////////////////////////
/// Transfer file from remote host. Returns number of bytes
/// received or < 0 in case of error. Error -1 connection is still
/// open, error -2 connection has been closed. In case of failure
/// fRestartAt is set to the number of bytes correclty transfered.
/// Calling GetFile() immediately afterwards will restart at fRestartAt.
/// If this is not desired call SetRestartAt(0) before calling GetFile().
/// If rootd reports that the file is locked, and you are sure this is not
/// the case (e.g. due to a crash), you can force unlock it by prepending
/// the file name with a '-'.

Long64_t TFTP::GetFile(const char *file, const char *localName)
{
   if (!IsOpen() || !file || !*file) return -1;

   if (!localName) {
      if (file[0] == '-')
         localName = file+1;
      else
         localName = file;
   }

   Long64_t restartat = fRestartAt;

   if (fSocket->Send(Form("%s %d %d %lld", file, fBlockSize, fMode,
                     restartat), kROOTD_GETFILE) < 0) {
      Error("GetFile", "error sending kROOTD_GETFILE command");
      return -2;
   }

   Int_t         stat;
   EMessageTypes kind;

   if (Recv(stat, kind) < 0 || kind == kROOTD_ERR) {
      PrintError("GetFile", stat);
      return -1;
   }

   // get size of remote file
   Long64_t size;
   Int_t    what;
   char     mess[128];

   if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
      Error("GetFile", "error receiving remote file size");
      return -2;
   }
#ifdef R__WIN32
   sscanf(mess, "%I64d", &size);
#else
   sscanf(mess, "%lld", &size);
#endif

   // check if restartat value makes sense
   if (restartat && (restartat >= size))
      restartat = 0;

   // open local file
   Int_t fd;
   if (!restartat) {
#if defined(R__WIN32) || defined(R__WINGCC)
      if (fMode == kBinary)
         fd = open(localName, O_CREAT | O_TRUNC | O_WRONLY | O_BINARY,
                   S_IREAD | S_IWRITE);
      else
         fd = open(localName, O_CREAT | O_TRUNC | O_WRONLY,
                   S_IREAD | S_IWRITE);
#elif defined(R__SEEK64)
      fd = open64(localName, O_CREAT | O_TRUNC | O_WRONLY, 0600);
#else
      fd = open(localName, O_CREAT | O_TRUNC | O_WRONLY, 0600);
#endif
   } else {
#if defined(R__WIN32) || defined(R__WINGCC)
      if (fMode == kBinary)
         fd = open(localName, O_WRONLY | O_BINARY, S_IREAD | S_IWRITE);
      else
         fd = open(localName, O_WRONLY, S_IREAD | S_IWRITE);
#elif defined(R__SEEK64)
      fd = open64(localName, O_WRONLY, 0600);
#else
      fd = open(localName, O_WRONLY, 0600);
#endif
   }

   if (fd < 0) {
      Error("GetFile", "cannot open %s", localName);
      // send urgent message to rootd to stop tranfer
      return -1;
   }

   // check file system space
   if (strcmp(localName, "/dev/null")) {
      Long_t id, bsize, blocks, bfree;
      if (gSystem->GetFsInfo(localName, &id, &bsize, &blocks, &bfree) == 0) {
         Long64_t space = (Long64_t)bsize * (Long64_t)bfree;
         if (space < size - restartat) {
            Error("GetFile", "not enough space to store file %s", localName);
            // send urgent message to rootd to stop tranfer
            close(fd);
            return -1;
         }
      } else
         Warning("GetFile", "could not determine if there is enough free space to store file");
   }

   // seek to restartat position
   if (restartat) {
#if defined(R__SEEK64)
      if (lseek64(fd, restartat, SEEK_SET) < 0) {
#elif defined(R__WIN32)
      if (_lseeki64(fd, restartat, SEEK_SET) < 0) {
#else
      if (lseek(fd, restartat, SEEK_SET) < 0) {
#endif
         Error("GetFile", "cannot seek to position %lld in file %s",
               restartat, localName);
         // if cannot seek send urgent message to rootd to stop tranfer
         close(fd);
         return -1;
      }
   }

   Info("GetFile", "getting file %s (%lld bytes, starting at %lld)",
        localName, size, restartat);

   TStopwatch timer;
   timer.Start();

   char *buf = new char[fBlockSize];
   char *buf2 = 0;
   if (fMode == kAscii)
      buf2 = new char[fBlockSize];

   Long64_t pos = restartat & ~(fBlockSize-1);
   Int_t skip = restartat - pos;

   while (pos < size) {
      Long64_t left = size - pos;
      if (left > fBlockSize)
         left = fBlockSize;

      Int_t n;
      while ((n = fSocket->RecvRaw(buf, Int_t(left-skip))) < 0 &&
             TSystem::GetErrno() == EINTR)
         TSystem::ResetErrno();

      if (n != Int_t(left-skip)) {
         Error("GetFile", "error receiving buffer of length %d, got %d",
               Int_t(left-skip), n);
         close(fd);
         delete [] buf; delete [] buf2;
         return -2;
      }

      // in case of ascii file, loop over buffer and remove \r's
      ssize_t siz;
      if (fMode == kAscii) {
         Int_t i = 0, j = 0;
         while (i < n) {
            if (buf[i] == '\r')
               i++;
            else
               buf2[j++] = buf[i++];
         }
         n = j;
         while ((siz = write(fd, buf2, n)) < 0 && TSystem::GetErrno() == EINTR)
            TSystem::ResetErrno();
      } else {
         while ((siz = write(fd, buf, n)) < 0 && TSystem::GetErrno() == EINTR)
            TSystem::ResetErrno();
      }

      if (siz < 0) {
         SysError("GetFile", "error writing file %s", localName);
         // send urgent message to rootd to stop tranfer
         close(fd);
         delete [] buf; delete [] buf2;
         return -1;
      }

      if (siz != n) {
         Error("GetFile", "error writing all requested bytes to file %s, wrote %ld of %d",
               localName, (Long_t)siz, n);
         // send urgent message to rootd to stop tranfer
         close(fd);
         delete [] buf; delete [] buf2;
         return -1;
      }

      fBytesRead  += left-skip;
      fgBytesRead += left-skip;

      fRestartAt = pos;   // bytes correctly received up till now

      pos += left;
      skip = 0;
   }

   delete [] buf; delete [] buf2;

#ifndef R__WIN32
   fchmod(fd, 0644);
#endif

   close(fd);

   fRestartAt = 0;

   // provide timing numbers
   Double_t speed, t = timer.RealTime();
   if (t > 0)
      speed = Double_t(size - restartat) / t;
   else
      speed = 0.0;
   if (speed > 524288)
      Info("GetFile", "%.3f seconds, %.2f Mbytes per second",
           t, speed / 1048576);
   else if (speed > 512)
      Info("GetFile", "%.3f seconds, %.2f Kbytes per second",
           t, speed / 1024);
   else
      Info("GetFile", "%.3f seconds, %.2f bytes per second",
           t, speed);

   return Long64_t(size - restartat);
}

////////////////////////////////////////////////////////////////////////////////
/// Change the remote directory. If the remote directory contains a .message
/// file and it is < 1024 characters then the contents is echoed back.
/// Returns 0 in case of success and -1 in case of failure.

Int_t TFTP::ChangeDirectory(const char *dir) const
{
   if (!IsOpen()) return -1;

   if (!dir || !*dir) {
      Error("ChangeDirectory", "illegal directory name specified");
      return -1;
   }

   if (fSocket->Send(Form("%s", dir), kROOTD_CHDIR) < 0) {
      Error("ChangeDirectory", "error sending kROOTD_CHDIR command");
      return -1;
   }

   Int_t what;
   char  mess[1024];

   if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
      Error("ChangeDirectory", "error receiving chdir confirmation");
      return -1;
   }
   if (what == kMESS_STRING) {
      Printf("%s\n", mess);

      if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
         Error("ChangeDirectory", "error receiving chdir confirmation");
         return -1;
      }
   }

   Info("ChangeDirectory", "%s", mess);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a remote directory. Anonymous users may not create directories.
/// Returns 0 in case of success and -1 in case of failure.

Int_t TFTP::MakeDirectory(const char *dir, Bool_t print) const
{
   if (!IsOpen()) return -1;

   if (!dir || !*dir) {
      Error("MakeDirectory", "illegal directory name specified");
      return -1;
   }

   if (fSocket->Send(Form("%s", dir), kROOTD_MKDIR) < 0) {
      Error("MakeDirectory", "error sending kROOTD_MKDIR command");
      return -1;
   }

   Int_t what;
   char  mess[1024];

   if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
      Error("MakeDirectory", "error receiving mkdir confirmation");
      return -1;
   }

   if (print)
      Info("MakeDirectory", "%s", mess);

   if (!strncmp(mess,"OK:",3))
      return 1;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a remote directory. Anonymous users may not delete directories.
/// Returns 0 in case of success and -1 in case of failure.

Int_t TFTP::DeleteDirectory(const char *dir) const
{
   if (!IsOpen()) return -1;

   if (!dir || !*dir) {
      Error("DeleteDirectory", "illegal directory name specified");
      return -1;
   }

   if (fSocket->Send(Form("%s", dir), kROOTD_RMDIR) < 0) {
      Error("DeleteDirectory", "error sending kROOTD_RMDIR command");
      return -1;
   }

   Int_t what;
   char  mess[1024];

   if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
      Error("DeleteDirectory", "error receiving rmdir confirmation");
      return -1;
   }

   Info("DeleteDirectory", "%s", mess);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// List remote directory. With cmd you specify the options and directory
/// to be listed to ls. Returns 0 in case of success and -1 in case of
/// failure.

Int_t TFTP::ListDirectory(Option_t *cmd) const
{
   if (!IsOpen()) return -1;

   if (!cmd || !*cmd)
      cmd = "ls .";

   if (fSocket->Send(Form("%s", cmd), kROOTD_LSDIR) < 0) {
      Error("ListDirectory", "error sending kROOTD_LSDIR command");
      return -1;
   }

   Int_t what;
   char  mess[1024];

   do {
      if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
         Error("ListDirectory", "error receiving lsdir confirmation");
         return -1;
      }
      printf("%s", mess);
   } while (what == kMESS_STRING);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Print path of remote working directory. Returns 0 in case of succes and
/// -1 in cse of failure.

Int_t TFTP::PrintDirectory() const
{
   if (!IsOpen()) return -1;

   if (fSocket->Send("", kROOTD_PWD) < 0) {
      Error("DeleteDirectory", "error sending kROOTD_PWD command");
      return -1;
   }

   Int_t what;
   char  mess[1024];

   if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
      Error("PrintDirectory", "error receiving pwd confirmation");
      return -1;
   }

   Info("PrintDirectory", "%s", mess);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Rename a remote file. Anonymous users may not rename files.
/// Returns 0 in case of success and -1 in case of failure.

Int_t TFTP::RenameFile(const char *file1, const char *file2) const
{
   if (!IsOpen()) return -1;

   if (!file1 || !file2 || !*file1 || !*file2) {
      Error("RenameFile", "illegal file names specified");
      return -1;
   }

   if (fSocket->Send(Form("%s %s", file1, file2), kROOTD_MV) < 0) {
      Error("RenameFile", "error sending kROOTD_MV command");
      return -1;
   }

   Int_t what;
   char  mess[1024];

   if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
      Error("RenameFile", "error receiving mv confirmation");
      return -1;
   }

   Info("RenameFile", "%s", mess);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a remote file. Anonymous users may not delete files.
/// Returns 0 in case of success and -1 in case of failure.

Int_t TFTP::DeleteFile(const char *file) const
{
   if (!IsOpen()) return -1;

   if (!file || !*file) {
      Error("DeleteFile", "illegal file name specified");
      return -1;
   }

   if (fSocket->Send(Form("%s", file), kROOTD_RM) < 0) {
      Error("DeleteFile", "error sending kROOTD_RM command");
      return -1;
   }

   Int_t what;
   char  mess[1024];

   if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
      Error("DeleteFile", "error receiving rm confirmation");
      return -1;
   }

   Info("DeleteFile", "%s", mess);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Change permissions of a remote file. Anonymous users may not
/// chnage permissions. Returns 0 in case of success and -1 in case
/// of failure.

Int_t TFTP::ChangePermission(const char *file, Int_t mode) const
{
   if (!IsOpen()) return -1;

   if (!file || !*file) {
      Error("ChangePermission", "illegal file name specified");
      return -1;
   }

   if (fSocket->Send(Form("%s %d", file, mode), kROOTD_CHMOD) < 0) {
      Error("ChangePermission", "error sending kROOTD_CHMOD command");
      return -1;
   }

   Int_t what;
   char  mess[1024];

   if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
      Error("ChangePermission", "error receiving chmod confirmation");
      return -1;
   }

   Info("ChangePermission", "%s", mess);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Close ftp connection. Returns 0 in case of success and -1 in case of
/// failure.

Int_t TFTP::Close()
{
   if (!IsOpen()) return -1;

   if (fSocket->Send(kROOTD_CLOSE) < 0) {
      Error("Close", "error sending kROOTD_CLOSE command");
      return -1;
   }

   // Ask for remote shutdown
   if (fProtocol > 6)
      fSocket->Send(kROOTD_BYE);

   // Remove from the list of Sockets
   {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfSockets()->Remove(this);
   }

   // Delete socket here
   SafeDelete(fSocket);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a directory via rootd.
/// Returns kTRUE in case of success.
/// Returns kFALSE in case of error.

Bool_t TFTP::OpenDirectory(const char *dir, Bool_t print)
{
   fDir = kFALSE;

   if (!IsOpen()) return fDir;

   if (fProtocol < 12) {
      Error("OpenDirectory", "call not supported by remote rootd");
      return fDir;
   }

   if (!dir || !*dir) {
      Error("OpenDirectory", "illegal directory name specified");
      return fDir;
   }

   if (fSocket->Send(Form("%s", dir), kROOTD_OPENDIR) < 0) {
      Error("OpenDirectory", "error sending kROOTD_OPENDIR command");
      return fDir;
   }

   Int_t what;
   char  mess[1024];;

   if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
      Error("OpenDirectory", "error receiving opendir confirmation");
      return fDir;
   }

   if (print)
      Info("OpenDirectory", "%s", mess);

   if (!strncmp(mess,"OK:",3)) {
      fDir = kTRUE;
      return fDir;
   }
   return fDir;
}

////////////////////////////////////////////////////////////////////////////////
/// Free a remotely open directory via rootd.

void TFTP::FreeDirectory(Bool_t print)
{
   if (!IsOpen() || !fDir) return;

   if (fProtocol < 12) {
      Error("FreeDirectory", "call not supported by remote rootd");
      return;
   }

   if (fSocket->Send(kROOTD_FREEDIR) < 0) {
      Error("FreeDirectory", "error sending kROOTD_FREEDIR command");
      return;
   }

   Int_t what;
   char  mess[1024];;

   if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
      Error("FreeDirectory", "error receiving freedir confirmation");
      return;
   }

   if (print)
      Info("FreeDirectory", "%s", mess);

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Get directory entry via rootd.
/// Returns 0 in case no more entries or in case of error.

const char *TFTP::GetDirEntry(Bool_t print)
{
   static char dirent[1024] = {0};

   if (!IsOpen() || !fDir) return 0;

   if (fProtocol < 12) {
      Error("GetDirEntry", "call not supported by remote rootd");
      return 0;
   }

   if (fSocket->Send(kROOTD_DIRENTRY) < 0) {
      Error("GetDirEntry", "error sending kROOTD_DIRENTRY command");
      return 0;
   }

   Int_t what;
   char  mess[1024];;

   if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
      Error("GetDirEntry", "error receiving dir entry confirmation");
      return 0;
   }

   if (print)
      Info("GetDirEntry", "%s", mess);

   if (!strncmp(mess,"OK:",3)) {
      strlcpy(dirent,mess+3, sizeof(dirent));
      return (const char *)dirent;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get info about a file. Info is returned in the form of a FileStat_t
/// structure (see TSystem.h).
/// The function returns 0 in case of success and 1 if the file could
/// not be stat'ed.

Int_t TFTP::GetPathInfo(const char *path, FileStat_t &buf, Bool_t print)
{
   TUrl url(path);

   if (!IsOpen()) return 1;

   if (fProtocol < 12) {
      Error("GetPathInfo", "call not supported by remote rootd");
      return 1;
   }

   if (!path || !*path) {
      Error("GetPathInfo", "illegal path name specified");
      return 1;
   }

   if (fSocket->Send(Form("%s", path), kROOTD_FSTAT) < 0) {
      Error("GetPathInfo", "error sending kROOTD_FSTAT command");
      return 1;
   }

   Int_t what;
   char  mess[1024];;

   if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
      Error("GetPathInfo", "error receiving fstat confirmation");
      return 1;
   }
   if (print)
      Info("GetPathInfo", "%s", mess);

   Int_t    mode, uid, gid, islink;
   Long_t   id, flags, dev, ino, mtime;
   Long64_t size;
   if (fProtocol > 12) {
#ifdef R__WIN32
      sscanf(mess, "%ld %ld %d %d %d %I64d %ld %d", &dev, &ino, &mode,
             &uid, &gid, &size, &mtime, &islink);
#else
      sscanf(mess, "%ld %ld %d %d %d %lld %ld %d", &dev, &ino, &mode,
             &uid, &gid, &size, &mtime, &islink);
#endif
      if (dev == -1)
         return 1;
      buf.fDev    = dev;
      buf.fIno    = ino;
      buf.fMode   = mode;
      buf.fUid    = uid;
      buf.fGid    = gid;
      buf.fSize   = size;
      buf.fMtime  = mtime;
      buf.fIsLink = (islink == 1);
   } else {
#ifdef R__WIN32
      sscanf(mess, "%ld %I64d %ld %ld", &id, &size, &flags, &mtime);
#else
      sscanf(mess, "%ld %lld %ld %ld", &id, &size, &flags, &mtime);
#endif
      if (id == -1)
         return 1;
      buf.fDev    = (id >> 24);
      buf.fIno    = (id & 0x00FFFFFF);
      if (flags == 0)
         buf.fMode = kS_IFREG;
      if (flags & 1)
         buf.fMode = (kS_IFREG|kS_IXUSR|kS_IXGRP|kS_IXOTH);
      if (flags & 2)
         buf.fMode = kS_IFDIR;
      if (flags & 4)
         buf.fMode = kS_IFSOCK;
      buf.fSize   = size;
      buf.fMtime  = mtime;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kFALSE if one can access a file using the specified access mode.
/// Mode is the same as for the Unix access(2) function.
/// Attention, bizarre convention of return value!!

Bool_t TFTP::AccessPathName(const char *path, EAccessMode mode, Bool_t print)
{
   if (!IsOpen()) return kTRUE;

   if (fProtocol < 12) {
      Error("AccessPathName", "call not supported by remote rootd");
      return kTRUE;
   }

   if (!path || !*path) {
      Error("AccessPathName", "illegal path name specified");
      return kTRUE;
   }

   if (fSocket->Send(Form("%s %d", path, mode), kROOTD_ACCESS) < 0) {
      Error("AccessPathName", "error sending kROOTD_ACCESS command");
      return kTRUE;
   }

   Int_t what;
   char  mess[1024];;

   if (fSocket->Recv(mess, sizeof(mess), what) < 0) {
      Error("AccessPathName", "error receiving access confirmation");
      return kTRUE;
   }
   if (print)
      Info("AccessPathName", "%s", mess);

   if (!strncmp(mess,"OK",2))
      return kFALSE;
   else
      return kTRUE;
}
