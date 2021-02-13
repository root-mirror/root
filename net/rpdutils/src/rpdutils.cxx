// @(#)root/rpdutils:$Id$
// Author: Gerardo Ganis    7/4/2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// rpdutils                                                             //
//                                                                      //
// Set of utilities for rootd/proofd daemon authentication.             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RConfigure.h"
#include "TError.h"
#include <ROOT/RConfig.hxx>
#include "strlcpy.h"

#include <ctype.h>
#include <fcntl.h>
#include <pwd.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <errno.h>
#include <netdb.h>
#include <math.h>
#include "snprintf.h"

#if defined(__CYGWIN__) && defined(__GNUC__)
#   define cygwingcc
#endif

#if defined(__linux__) && !defined(linux)
# define linux
#endif

#if defined(linux) || defined(__sun) || defined(__sgi) || \
    defined(_AIX) || defined(__FreeBSD__) || defined(__OpenBSD__) || \
    defined(__APPLE__) || defined(__MACH__) || defined(cygwingcc)
#include <grp.h>
#include <sys/types.h>
#include <signal.h>
#endif

#ifdef _AIX
extern "C" int ruserok(char *, int, char *, char *);
#endif

#if defined(__APPLE__)
#include <sys/mount.h>
extern "C" int fstatfs(int file_descriptor, struct statfs *buffer);
#elif defined(linux) || defined(__hpux)
#include <sys/vfs.h>
#elif defined(__FreeBSD__) || defined(__OpenBSD__)
#include <sys/param.h>
#include <sys/mount.h>
#else
#include <sys/statfs.h>
#endif

#if defined(linux)
#   include <features.h>
#   if __GNU_LIBRARY__ == 6
#      ifndef R__GLIBC
#         define R__GLIBC
#      endif
#   endif
#endif
#if defined(cygwingcc) || (defined(__MACH__) && !defined(__APPLE__))
#   define R__GLIBC
#endif

#if (defined(__FreeBSD__) && (__FreeBSD__ < 4)) || defined(__OpenBSD__) || \
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

#if defined(cygwingcc) && !defined(F_LOCK) && !defined(F_ULOCK)
#define F_LOCK F_WRLCK
#define F_ULOCK F_UNLCK
int ruserok(const char *, int, const char *, const char *) {
   return 0;
}
static int fcntl_lockf(int fd, int op, off_t off)
{
   flock fl;
   fl.l_whence = SEEK_SET;
   fl.l_start  = off;
   fl.l_len    = 0;       // whole file
   fl.l_pid    = getpid();
   fl.l_type   = op;
   return fcntl(fd, F_SETLK, &fl);
}
#define lockf fcntl_lockf
#endif

#if defined(__sun) || defined(R__GLIBC)
#include <crypt.h>
#endif

#if defined(__osf__) || defined(__sgi) || defined(R__MACOSX)
extern "C" char *crypt(const char *, const char *);
#endif

#ifdef R__WIN32
#define R__NOCRYPT
#endif

#ifdef R__NOCRYPT
static std::string gRndmSalt = std::string("ABCDEFGH");
#endif

#if defined(__sun)
#ifndef R__SHADOWPW
#define R__SHADOWPW
#endif
#endif

#ifdef R__SHADOWPW
#include <shadow.h>
#endif

#ifdef R__SSL
// SSL specific headers for RSA keys
#include <openssl/bio.h>
#include <openssl/blowfish.h>
#include <openssl/err.h>
#include <openssl/pem.h>
#include <openssl/rand.h>
#include <openssl/rsa.h>
#include <openssl/ssl.h>
#endif

#include "rpdp.h"
#include "rsadef.h"
#include "rsalib.h"
//
// To improve error logging for UsrPwd on the client side
static ERootdErrors gUsrPwdErr[4][4] = {
   {kErrNoPasswd, kErrNoPassHEquNoFiles, kErrNoPassHEquBadFiles, kErrNoPassHEquFailed},
   {kErrBadPasswd, kErrBadPassHEquNoFiles, kErrBadPassHEquBadFiles, kErrBadPassHEquFailed},
   {kErrBadRtag, kErrBadRtagHEquNoFiles, kErrBadRtagHEquBadFiles, kErrBadRtagHEquFailed},
   {kErrBadPwdFile, kErrBadPwdFileHEquNoFiles, kErrBadPwdFileHEquBadFiles,
    kErrBadPwdFileHEquFailed}};

//--- Machine specific routines ------------------------------------------------

#if defined(__sgi) && !defined(__GNUG__) && (SGI_REL<62)
extern "C" {
   int seteuid(int euid);
   int setegid(int egid);
}
#endif

#if defined(_AIX)
extern "C" {
   int seteuid(uid_t euid);
   int setegid(gid_t egid);
}
#endif

#if !defined(__hpux) && !defined(linux) && !defined(__FreeBSD__) && \
    !defined(__OpenBSD__) || defined(cygwingcc)
static int setresgid(gid_t r, gid_t e, gid_t)
{
   if (setgid(r) == -1)
      return -1;
   return setegid(e);
}

static int setresuid(uid_t r, uid_t e, uid_t)
{
   if (setuid(r) == -1)
      return -1;
   return seteuid(e);
}
#else
#if defined(linux) && !defined(R__HAS_SETRESUID)
extern "C" {
   int setresgid(gid_t r, gid_t e, gid_t s);
   int setresuid(uid_t r, uid_t e, uid_t s);
}
#endif
#endif

#if defined(__sun)
#if defined(R__SUNGCC3)
extern "C" int gethostname(char *, unsigned int);
#endif
#endif

extern int gDebug;

namespace ROOT {

//
// rpdutils module globals
ErrorHandler_t gErrSys   = 0;
ErrorHandler_t gErrFatal = 0;
ErrorHandler_t gErr      = 0;
bool gSysLog = 0;
std::string gServName[3] = { "sockd", "rootd", "proofd" };

//
// Local global consts
static const int gAUTH_CLR_MSK = 0x1;     // Masks for authentication methods
static const int gMAXTABSIZE = 50000000;

static const std::string gAuthMeth[kMAXSEC] = { "UsrPwd", "Unsupported", "Unsupported",
                                                "Unsupported", "Unsupported", "Unsupported" };
static const std::string gAuthTab    = "/rpdauthtab";   // auth table
static const std::string gDaemonRc   = ".rootdaemonrc"; // daemon access rules
static const std::string gRootdPass  = ".rootdpass";    // special rootd passwd
static const std::string gKeyRoot    = "/rpk.";         // Root for key files

//
// RW dir for temporary files (needed by gRpdAuthTab: do not move)
static std::string gTmpDir = "/tmp";

//
// Local global vars
static int gAuthProtocol = -1;  // Protocol used fro a successful authentication
static char gBufOld[kMAXRECVBUF] = {0}; // msg sync for old client (<=3.05/07)
static bool gCheckHostsEquiv = 1;
static int gClientOld = 0;              // msg sync for old client (<=3.05/07)
static int gClientProtocol = -1;
static int gCryptRequired = -1;
static std::string gCryptToken;
static int gAllowMeth[kMAXSEC];
static int gAnon = 0;
static int gExistingAuth = 0;
static int gAuthListSent = 0;
static int gHaveMeth[kMAXSEC];
static EMessageTypes gKindOld;          // msg sync for old client (<=3.05/07)
static int gMethInit = 0;
static int gNumAllow = -1;
static int gNumLeft = -1;
static int gOffSet = -1;
static std::string gOpenHost = "????";
static int gParentId = -1;
static char gPasswd[kMAXUSERLEN] = { 0 };
static char gPubKey[kMAXPATHLEN] = { 0 };
static int gPubKeyLen = 0;
static int gRandInit = 0;
static int gRemPid = -1;
static bool gRequireAuth = 1;
static int gReUseAllow = 0x1F;  // define methods for which tokens can be asked
static int gReUseRequired = -1;
static int gDoLogin = 0;  // perform login
static std::string gRpdAuthTab = std::string(gTmpDir).append(gAuthTab);
static std::string gRpdKeyRoot = std::string(gTmpDir).append(gKeyRoot);
static rsa_NUMBER gRSA_d;
static rsa_NUMBER gRSA_n;
static int gRSAInit = 0;
static int gRSAKey = 0;
static rsa_KEY gRSAPriKey;
static rsa_KEY_export gRSAPubExport[2] = {{0,0},{0,0}};
static rsa_KEY gRSAPubKey;
#ifdef R__SSL
static BF_KEY gBFKey;            // Session symmetric key
static RSA *gRSASSLKey = 0;      // Local RSA SSL key
#endif
static int gSaltRequired = -1;
static int gSec = -1;
static int gServerProtocol = -1;
static EService gService = kROOTD;
static int gTriedMeth[kMAXSEC];
static char gUser[64] = { 0 };
static char *gUserAllow[kMAXSEC] = { 0 };          // User access control
static unsigned int gUserAlwLen[kMAXSEC] = { 0 };
static unsigned int gUserIgnLen[kMAXSEC] = { 0 };
static char *gUserIgnore[kMAXSEC] = { 0 };

////////////////////////////////////////////////////////////////////////////////
/// rand() implementation using /udev/random or /dev/random, if available

static int rpd_rand()
{
#ifndef WIN32
   int frnd = open("/dev/urandom", O_RDONLY);
   if (frnd < 0) frnd = open("/dev/random", O_RDONLY);
   int r;
   if (frnd >= 0) {
      ssize_t rs = read(frnd, (void *) &r, sizeof(int));
      close(frnd);
      if (r < 0) r = -r;
      if (rs == sizeof(int)) return r;
   }
   ErrorInfo("+++ERROR+++ : rpd_rand: neither /dev/urandom nor /dev/random are available or readable!");
   struct timeval tv;
   if (gettimeofday(&tv,0) == 0) {
      int t1, t2;
      memcpy((void *)&t1, (void *)&tv.tv_sec, sizeof(int));
      memcpy((void *)&t2, (void *)&tv.tv_usec, sizeof(int));
      r = t1 + t2;
      if (r < 0) r = -r;
      return r;
   }
   return -1;
#else
   // No special random device available: use rand()
   return rand();
#endif
}

////////////////////////////////////////////////////////////////////////////////
///  reads in at most one less than len characters from open
///  descriptor fd and stores them into the buffer pointed to by buf.
///  Reading stops after an EOF or a newline. If a newline is
///  read, it  is stored into the buffer.
///  A '\0' is stored after the last character in the buffer.
///  The number of characters read is returned (newline included).
///  Returns < 0 in case of error.

static int reads(int fd, char *buf, int len)
{
   int k = 0;
   int nread = -1;
   int nr = read(fd,buf,1);
   while (nr > 0 && buf[k] != '\n' && k < (len-1)) {
      k++;
      nr = read(fd,buf+k,1);
   }
   if (k >= len-1) {
      buf[k] = 0;
      nread = len-1;
   } else if (buf[k] == '\n'){
      if (k <= len-2) {
         buf[k+1] = 0;
         nread = k+1;
      } else {
         buf[k] = 0;
         nread = k;
      }
   } else if (nr == 0) {
      if (k > 0) {
         buf[k-1] = 0;
         nread = k-1;
      } else {
         buf[0] = 0;
         nread = 0;
      }
   } else if (nr < 0) {
      if (k > 0) {
         buf[k] = 0;
         nread = -(k-1);
      } else {
         buf[0] = 0;
         nread = -1;
      }
   }
   // Fix the lengths
   if (nread >= 0) buf[nread] = 0;

   return nread;
}

////////////////////////////////////////////////////////////////////////////////
/// Case insensitive string compare of n characters.

static int rpdstrncasecmp(const char *str1, const char *str2, int n)
{
   while (n > 0) {
      int c1 = *str1;
      int c2 = *str2;

      if (isupper(c1))
         c1 = tolower(c1);

      if (isupper(c2))
         c2 = tolower(c2);

      if (c1 != c2)
         return c1 - c2;

      str1++;
      str2++;
      n--;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Case insensitive string compare.

static int rpdstrcasecmp(const char *str1, const char *str2)
{
   return rpdstrncasecmp(str1, str2, strlen(str2) + 1);
}

////////////////////////////////////////////////////////////////////////////////
/// To avoid problems due to compiler optmization
/// Taken from Viega&Messier, "Secure Programming Cookbook", O'Really, #13.2
/// (see discussion there)

static volatile void *rpdmemset(volatile void *dst, int c, int len)
{
   volatile char *buf;

   for (buf = (volatile char *)dst; len; (buf[--len] = c)) { }
   return dst;
}

#ifdef R__NOCRYPT
////////////////////////////////////////////////////////////////////////////////
/// This applies simple nor encryption with sa to the first 64 bytes
/// pw. Returns the hex of the result (max length 128).
/// This is foreseen for systms where crypt is not available
/// (on windows ...), to provide some protection of tokens.

char *rpdcrypt(const char *pw, const char *sa)
{
   static char buf[129];
   char tbuf[64];
   int np = (strlen(pw) < 64) ? strlen(pw) : 64;
   int ns = strlen(sa);
   char c;

   int i = 0;
   for (i=0; i<np; i++) {
      // We use in turn all the salt bits; and re-use
      // if they are not enough
      int ks = i%ns;
      tbuf[i] = pw[i]^sa[ks];
      // Convert the result in two hexs: the first ...
      int j = 0xF & tbuf[i];
      if (j < 10)
         c = 48 + j;
      else
         c = 55 + j;
      int k = 2*i;
      buf[k] = c;
      // .. the second
      j = (0xF0 & tbuf[i]) >> 4;
      if (j < 10)
         c = 48 + j;
      else
         c = 55 + j;
      k = 2*i + 1;
      buf[k] = c;
   }
   // Null termination
   buf[np*2] = 0;

   return buf;
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// Change the value of the static gSysLog to syslog.
/// Recognized values:
///                       0      log to syslog (for root started daemons)
///                       1      log to stderr (for user started daemons)

void RpdSetSysLogFlag(int syslog)
{
   gSysLog = syslog;
   if (gDebug > 2)
      ErrorInfo("RpdSetSysLogFlag: gSysLog set to %d", gSysLog);
}

////////////////////////////////////////////////////////////////////////////////
/// Change the value of the static gMethInit to methinit.
/// Recognized values:
///                       0      reset
///                       1      initialized already

void RpdSetMethInitFlag(int methinit)
{
   gMethInit = methinit;
   if (gDebug > 2)
      ErrorInfo("RpdSetMethInitFlag: gMethInit set to %d", gMethInit);
}
////////////////////////////////////////////////////////////////////////////////
/// Return pointer to the root string for key files
/// Used by proofd.

const char *RpdGetKeyRoot()
{
   return (const char *)gRpdKeyRoot.c_str();
}

////////////////////////////////////////////////////////////////////////////////
/// Return protocol version run by the client.
/// Used by proofd.

int RpdGetClientProtocol()
{
   return gClientProtocol;
}

////////////////////////////////////////////////////////////////////////////////
/// Return authentication protocol used for the handshake.
/// Used by proofd.

int RpdGetAuthProtocol()
{
   return gAuthProtocol;
}

////////////////////////////////////////////////////////////////////////////////
/// Return offset in the authtab file.
/// Used by proofd.

int RpdGetOffSet()
{
   return gOffSet;
}

////////////////////////////////////////////////////////////////////////////////

int RpdGetAuthMethod(int kind)
{
   int method = -1;

   if (kind == kROOTD_USER)
      method = 0;

   return method;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete Public Key file
/// Returns: 0 if ok
///          1 if error unlinking (check errno);

int RpdDeleteKeyFile(int ofs)
{
   int retval = 0;

   std::string pukfile = gRpdKeyRoot;
   pukfile.append(ItoA(ofs));

   // Some debug info
   if (gDebug > 2) {
      ErrorInfo("RpdDeleteKeyFile: proc uid:%d gid:%d",
      getuid(),getgid());
   }

   // Unlink
   if (unlink(pukfile.c_str()) == -1) {
      if (gDebug > 0 && GetErrno() != ENOENT) {
         ErrorInfo("RpdDeleteKeyFile: problems unlinking pub"
                   " key file '%s' (errno: %d)",
                   pukfile.c_str(),GetErrno());
      }
      retval = 1;
   }
   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Update tab file.
/// If ilck <= 0 open and lock the file; if ilck > 0, use file
/// descriptor ilck, which should correspond to an open and locked file.
/// If opt = -1 : delete file (backup saved in <file>.bak);
/// If opt =  0 : eliminate all inactive entries
///               (if line="size" act only if size > gMAXTABSIZE)
/// if opt =  1 : append 'line'.
/// Returns -1 in case of error.
/// Returns offset for 'line' and token for opt = 1.
/// Returns new file size for opt = 0.

int RpdUpdateAuthTab(int opt, const char *line, char **token, int ilck)
{
   int retval = -1;
   int itab = 0;
   char fbuf[kMAXPATHLEN];

   if (gDebug > 2)
      ErrorInfo("RpdUpdateAuthTab: analyzing: opt: %d, line: %s, ilck: %d",
                opt, line, ilck);

   if (ilck <= 0) {

      // Open file for reading/ writing
      itab = open(gRpdAuthTab.c_str(), O_RDWR);
      if (itab == -1) {
         if (opt == 1 && GetErrno() == ENOENT) {
            // Try creating the file
            itab = open(gRpdAuthTab.c_str(), O_RDWR | O_CREAT, 0600);
            if (itab == -1) {
               ErrorInfo("RpdUpdateAuthTab: opt=%d: error opening %s"
                         "(errno: %d)",
                         opt, gRpdAuthTab.c_str(), GetErrno());
               return retval;
            }
         } else {
            ErrorInfo("RpdUpdateAuthTab: opt=%d: error opening %s"
                      " (errno: %d)",
                       opt, gRpdAuthTab.c_str(), GetErrno());
            return retval;
         }
      }

      // lock tab file
      if (lockf(itab, F_LOCK, (off_t) 1) == -1) {
         ErrorInfo("RpdUpdateAuthTab: opt=%d: error locking %s"
                   " (errno: %d)", opt, gRpdAuthTab.c_str(), GetErrno());
         close(itab);
         return retval;
      }
      if (gDebug > 0)
         ErrorInfo("RpdUpdateAuthTab: opt= %d - file LOCKED", opt);
   } else {
      itab = ilck;
   }

   // File size
   int fsize = 0;
   if ((fsize = lseek(itab, 0, SEEK_END)) == -1) {
      ErrorInfo("RpdUpdateAuthTab: opt=%d: lseek error (errno: %d)",
           opt, GetErrno());
      goto goingout;
   }

   // Set indicator to beginning
   if (lseek(itab, 0, SEEK_SET) == -1) {
      ErrorInfo("RpdUpdateAuthTab: opt=%d: lseek error (errno: %d)",
           opt, GetErrno());
      goto goingout;
   }

   if (opt == -1) {

      //
      // Save file in .bak and delete its content

      // Open backup file
      std::string bak = std::string(gRpdAuthTab).append(".bak");
      int ibak = open(bak.c_str(), O_RDWR | O_CREAT, 0600);
      if (ibak == -1) {
         ErrorInfo("RpdUpdateAuthTab: opt=%d: error opening/creating %s"
                   " (errno: %d)", opt, bak.c_str(), GetErrno());
         goto goingout;
      }

      // Truncate file to new length
      if (ftruncate(ibak, 0) == -1)
         ErrorInfo("RpdUpdateAuthTab: opt=%d: ftruncate error (%s)"
                   " (errno: %d)", opt, bak.c_str(), GetErrno());

      // Copy the content
      char buf[kMAXPATHLEN];
      int ofs = 0, nr = 0;
      while ((nr = reads(itab, buf, sizeof(buf)))) {
         int slen = strlen(buf);

         // Make sure there is a '\n' before writing
         if (buf[slen-1] != '\n') {
            if (slen >= kMAXPATHLEN -1)
               buf[slen-1] = '\n';
            else {
               buf[slen] = '\n';
               buf[slen+1] = '\0';
            }
         }
         if (slen) {
            while (write(ibak, buf, slen) < 0 && GetErrno() == EINTR)
               ResetErrno();
         }

         // Delete Public Key file
         RpdDeleteKeyFile(ofs);
         // Next OffSet
         ofs += slen;
      }
      close(ibak);

      // Truncate file to new length
      if (ftruncate(itab, 0) == -1)
         ErrorInfo("RpdUpdateAuthTab: opt=%d: ftruncate error (%s)"
                   " (errno: %d)", opt, gRpdAuthTab.c_str(), GetErrno());
      retval = 0;

   } else if (opt == 0) {

      //
      // Cleanup the file (remove inactive entries)

      // Now scan over entries
      int pr = 0, pw = 0;
      int lsec, act = 0, oldofs = 0, bytesread = 0;
      char ln[kMAXPATHLEN], dumm[kMAXPATHLEN];
      bool fwr = 0;

      while ((bytesread = reads(itab, ln, sizeof(ln)))) {

         bool ok = 1;
         // Current position
         if ((pr = lseek(itab,0,SEEK_CUR)) < 0) {
            ErrorInfo("RpdUpdateAuthTab: opt=%d: problems lseeking file %s"
                      " (errno: %d)", opt, gRpdAuthTab.c_str(), errno);
            fwr = 1;
            ok = 0;
         }

         // Check file corruption: length and number of items
         int slen = bytesread;
         if (ok && slen < 1) {
            ErrorInfo("RpdUpdateAuthTab: opt=%d: file %s seems corrupted"
                      " (slen: %d)", opt, gRpdAuthTab.c_str(), slen);
            fwr = 1;
            ok = 0;
         }
         if (ok) {
            // Check file corruption: number of items
            int ns = sscanf(ln, "%d %d %4095s", &lsec, &act, dumm);
            if (ns < 3 ) {
               ErrorInfo("RpdUpdateAuthTab: opt=%d: file %s seems corrupted"
                         " (ns: %d)", opt, gRpdAuthTab.c_str(), ns);
               fwr = 1;
               ok = 0;
            }
         }

         if (ok && act > 0) {
            if (fwr) {
               // We have to update the key file name
               int nr = 0;
               if ((nr = RpdRenameKeyFile(oldofs,pw)) == 0) {
                  // Write the entry at new position
                  lseek(itab, pw, SEEK_SET);

                  if (ln[slen-1] != '\n') {
                     if (slen >= kMAXPATHLEN -1)
                        ln[slen-1] = '\n';
                     else {
                        ln[slen] = '\n';
                        ln[slen+1] = '\0';
                     }
                  }
                  while (write(itab, ln, strlen(ln)) < 0
                         && GetErrno() == EINTR)
                     ResetErrno();
                  pw += strlen(ln);
               } else
                  RpdDeleteKeyFile(oldofs);
               lseek(itab, pr, SEEK_SET);
            } else
               pw += strlen(ln);
         } else {
            fwr = 1;
         }
         // Set old offset
         oldofs = pr;
      }

      // Truncate file to new length
      if (ftruncate(itab, pw) == -1)
         ErrorInfo("RpdUpdateAuthTab: opt=%d: ftruncate error (errno: %d)",
              opt, GetErrno());

      // Return new file size
      retval = pw;

   } else if (opt == 1) {

      //
      // Add 'line' at the end
      // (check size and cleanup/truncate if needed)

      // Check size ...
      if ((int)(fsize+strlen(line)) > gMAXTABSIZE) {

         // If it is going to be too big, cleanup or truncate first
         fsize = RpdUpdateAuthTab(0,(const char *)0,0,itab);

         // If still too big: delete everything
         if ((int)(fsize+strlen(line)) > gMAXTABSIZE)
            fsize = RpdUpdateAuthTab(-1,(const char *)0,0,itab);
      }
      // We are going to write at the end
      retval = lseek(itab, 0, SEEK_END);

      // Save first RSA public key into file for later use by the
      // same or other rootd/proofd; we will update the tab file
      // only if this operation is successful
      int ntry = 10;
      int rs = 0;
      while ((rs = RpdSavePubKey(gPubKey, retval, gUser)) == 2 && ntry--) {
         // We are here if a file with the same name exists already
         // and can not be deleted: we shift the offset with a
         // dummy entry
         char ltmp[256];
         SPrintf(ltmp, 256,
                 "0 0 %d %d %s error: pubkey file in use: shift offset\n",
                 gRSAKey, gRemPid, gOpenHost.c_str());

         // adds line
         while (write(itab, ltmp, strlen(ltmp)) < 0 && GetErrno() == EINTR)
            ResetErrno();

         // Set to the new end
         retval = lseek(itab, 0, SEEK_END);
      }

      if (rs > 0) {
         // Something wrong
         retval = -1;
         if (gDebug > 0)
            ErrorInfo("RpdUpdateAuthTab: pub key could not be saved (%d)",rs);
      } else {
         // Generate token
         *token = RpdGetRandString(3, 8);   // 8 crypt-like chars
#ifndef R__NOCRYPT
         char *cryptToken = crypt(*token, *token);
#else
         char *cryptToken = rpdcrypt(*token,gRndmSalt.c_str());
#endif
         SPrintf(fbuf, kMAXPATHLEN, "%s %s\n", line, cryptToken);
         if (gDebug > 2)
            ErrorInfo("RpdUpdateAuthTab: token: '%s'", cryptToken);
         // Save it for later use in kSOCKD servers
         gCryptToken = std::string(cryptToken);

         // adds line
         while (write(itab, fbuf, strlen(fbuf)) < 0 && GetErrno() == EINTR)
            ResetErrno();
      }

   } else {

      //
      // Unknown option
      ErrorInfo("RpdUpdateAuthTab: unrecognized option (opt= %d)", opt);
   }

   goingout:
   if (itab != ilck) {
      // unlock the file
      lseek(itab, 0, SEEK_SET);
      if (lockf(itab, F_ULOCK, (off_t) 1) == -1) {
         ErrorInfo("RpdUpdateAuthTab: error unlocking %s",
                   gRpdAuthTab.c_str());
      }

      // closing file ...
      close(itab);
   }

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// De-activates entry related to token with crypt crypttoken.
/// Returns: 0 if successful
///         -4 if entry not found or inactive
///         -1 problems opening auth tab file
///         -2 problems locking auth tab file
///         -3 auth tab file does not exists

int RpdCleanupAuthTab(const char *crypttoken)
{
   int retval = -4;

   if (gDebug > 2)
      ErrorInfo("RpdCleanupAuthTab: Crypt-token: '%s'",crypttoken);

   // Open file for update
   int itab = -1;
   if ((itab = open(gRpdAuthTab.c_str(), O_RDWR)) == -1) {
      if (GetErrno() == ENOENT) {
         if (gDebug > 0)
            ErrorInfo("RpdCleanupAuthTab: file %s does not exist",
                       gRpdAuthTab.c_str());
         return -3;
      } else {
         ErrorInfo("RpdCleanupAuthTab: error opening %s (errno: %d)",
                  gRpdAuthTab.c_str(), GetErrno());
         return -1;
      }
   }

   // lock tab file
   if (lockf(itab, F_LOCK, (off_t) 1) == -1) {
      ErrorInfo("RpdCleanupAuthTab: error locking %s (errno: %d)",
                gRpdAuthTab.c_str(), GetErrno());
      close(itab);
      return -2;
   }
   if (gDebug > 0)
      ErrorInfo("RpdCleanupAuthTab: file LOCKED (ctkn: '%s')",crypttoken);


   // Now access entry or scan over entries
   int pr = 0, pw = 0;
   int nw, lsec, act, remid, pkey;
   char line[kMAXPATHLEN];

   // Set indicators
   if ((pr = lseek(itab, 0, SEEK_SET)) < 0) {
      ErrorInfo("RpdCleanupAuthTab: error lseeking %s (errno: %d)",
                gRpdAuthTab.c_str(), GetErrno());
      close(itab);
      return -2;
   }
   pw = pr;
   while (reads(itab,line, sizeof(line))) {

      pr += strlen(line);
      if (gDebug > 2)
         ErrorInfo("RpdCleanupAuthTab: pr:%d pw:%d (line:%s) (pId:%d)",
                    pr, pw, line, gParentId);

      char dum1[kMAXPATHLEN] = {0}, host[kMAXUSERLEN] = {0}, user[kMAXUSERLEN] = {0},
           ctkn[30] = {0}, dum2[30] = {0};
      nw = sscanf(line, "%d %d %d %d %127s %127s %29s %4095s %29s",
                        &lsec, &act, &pkey, &remid, host, user, ctkn, dum1, dum2);

      int deactivate = 0;

      if (act > 0) {

         if (lsec == 3 && nw == 9) {
            if (!strncmp(dum2,crypttoken,strlen(crypttoken)))
               deactivate = 1;
         } else if (nw == 7) {
            if (!strncmp(ctkn,crypttoken,strlen(crypttoken)))
               deactivate = 1;
         }

         // Deactivate active entries: remote client has gone ...
         if (deactivate) {

            retval = 0;

            // Delete Public Key file
            RpdDeleteKeyFile(pw);

            // Locate 'act' ... skeep initial spaces, if any
            int slen = (int)strlen(line);
            int ka = 0;
            while (ka < slen && line[ka] == 32)
               ka++;
            // skeep method
            while (ka < slen && line[ka] != 32)
               ka++;
            // skeep spaces before 'act'
            while (ka < slen && line[ka] == 32)
               ka++;
            // This is 'act'
            line[ka] = '0';
            // Make sure there is a '\n' before writing
            int sl = strlen(line);
            if (line[sl-1] != '\n') {
               if (sl >= kMAXPATHLEN -1)
                  line[sl-1] = '\n';
               else {
                  line[sl] = '\n';
                  line[sl+1] = '\0';
               }
            }
            // Write it now
            lseek(itab, pw, SEEK_SET);
            while (write(itab, line, strlen(line)) < 0
                   && GetErrno() == EINTR)
               ResetErrno();
            // We are done
            lseek(itab,  0, SEEK_END);
         }
      }
      pw = pr;
   }

   // unlock the file
   lseek(itab, 0, SEEK_SET);
   if (lockf(itab, F_ULOCK, (off_t) 1) == -1) {
      ErrorInfo("RpdCleanupAuthTab: error unlocking %s", gRpdAuthTab.c_str());
   }
   // closing file ...
   close(itab);

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// In tab file, cleanup (set inactive) entry at offset
/// 'OffSet' from remote PiD 'RemId' at 'Host'.
/// If Host="all" or RemId=0 discard all entries.
/// Return number of entries not cleaned properly ...

int RpdCleanupAuthTab(const char *Host, int RemId, int OffSet)
{
   int retval = 0;

   if (gDebug > 2)
      ErrorInfo("RpdCleanupAuthTab: Host: '%s', RemId:%d, OffSet: %d",
                Host, RemId, OffSet);

   // Open file for update
   int itab = -1;
   if ((itab = open(gRpdAuthTab.c_str(), O_RDWR)) == -1) {
      if (GetErrno() == ENOENT) {
         if (gDebug > 0)
            ErrorInfo("RpdCleanupAuthTab: file %s does not exist",
                       gRpdAuthTab.c_str());
         return -3;
      } else {
         ErrorInfo("RpdCleanupAuthTab: error opening %s (errno: %d)",
                  gRpdAuthTab.c_str(), GetErrno());
         return -1;
      }
   }

   // lock tab file
   if (lockf(itab, F_LOCK, (off_t) 1) == -1) {
      ErrorInfo("RpdCleanupAuthTab: error locking %s (errno: %d)",
                gRpdAuthTab.c_str(), GetErrno());
      close(itab);
      //     return retval;
      return -2;
   }
   if (gDebug > 0)
      ErrorInfo("RpdCleanupAuthTab: file LOCKED"
                " (Host: '%s', RemId:%d, OffSet: %d)",
                  Host, RemId, OffSet);

   // Now access entry or scan over entries
   int pr = 0, pw = 0;
   int nw, lsec, act, remid, pkey;
   char line[kMAXPATHLEN];

   // Clean all flag
   int all = (!strcmp(Host, "all") || RemId == 0);

   // Set indicator
   if (all || OffSet < 0)
      pr = lseek(itab, 0, SEEK_SET);
   else
      pr = lseek(itab, OffSet, SEEK_SET);
   if (pr < 0) {
      ErrorInfo("RpdCleanupAuthTab: error lseeking %s (errno: %d)",
                gRpdAuthTab.c_str(), GetErrno());
      close(itab);
      //     return retval;
      return -2;
   }
   pw = pr;
   while (reads(itab,line, sizeof(line))) {

      pr += strlen(line);
      if (gDebug > 2)
         ErrorInfo("RpdCleanupAuthTab: pr:%d pw:%d (line:%s) (pId:%d)",
                    pr, pw, line, gParentId);

      char dumm[kMAXPATHLEN], host[kMAXUSERLEN], user[kMAXUSERLEN], shmbuf[30];
      nw = sscanf(line, "%d %d %d %d %127s %127s %29s %4095s",
                        &lsec, &act, &pkey, &remid, host, user, shmbuf, dumm);

      if (nw > 5) {
         if (all || OffSet > -1 ||
            (strstr(line,Host) && (RemId == remid))) {

            // Delete Public Key file
            RpdDeleteKeyFile(pw);

            // Deactivate active entries: remote client has gone ...
            if (act > 0) {

               // Locate 'act' ... skeep initial spaces, if any
               int slen = (int)strlen(line);
               int ka = 0;
               while (ka < slen && line[ka] == 32)
                  ka++;
               // skeep method
               while (ka < slen && line[ka] != 32)
                  ka++;
               // skeep spaces before 'act'
               while (ka < slen && line[ka] == 32)
                  ka++;
               // This is 'act'
               line[ka] = '0';
               // Make sure there is a '\n' before writing
               int sl = strlen(line);
               if (line[sl-1] != '\n') {
                  if (sl >= kMAXPATHLEN -1)
                     line[sl-1] = '\n';
                  else {
                     line[sl] = '\n';
                     line[sl+1] = '\0';
                  }
               }
               // Write it now
               lseek(itab, pw, SEEK_SET);
               while (write(itab, line, strlen(line)) < 0
                      && GetErrno() == EINTR)
                  ResetErrno();
               if (all || OffSet < 0)
                  lseek(itab, pr, SEEK_SET);
               else
                  lseek(itab,  0, SEEK_END);
            }
         }
      }
      pw = pr;
   }

   // unlock the file
   lseek(itab, 0, SEEK_SET);
   if (lockf(itab, F_ULOCK, (off_t) 1) == -1) {
      ErrorInfo("RpdCleanupAuthTab: error unlocking %s", gRpdAuthTab.c_str());
   }
   // closing file ...
   close(itab);

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Check authentication entry in tab file.

int RpdCheckAuthTab(int Sec, const char *User, const char *Host, int RemId,
                    int *OffSet)
{
   int retval = 0;
   if (gDebug > 2)
      ErrorInfo("RpdCheckAuthTab: analyzing: %d %s %s %d %d", Sec, User,
                Host, RemId, *OffSet);

   // Check OffSet first
   char *tkn = 0, *user =0;
   int shmid;
   bool goodOfs = RpdCheckOffSet(Sec,User,Host,RemId,
                                 OffSet,&tkn,&shmid,&user);
   if (gDebug > 2)
      ErrorInfo("RpdCheckAuthTab: goodOfs: %d", goodOfs);

   // Notify the result of the check
   int tag = 0;
   if (gClientProtocol >= 10) {
      if (goodOfs) {
         if (gClientProtocol > 11) {
            // Generate tag
            RpdInitRand();
            while ((tag = rpd_rand()) == 1) ; // .ne.1 for backward comptibility

            // We will receive the user token next
            NetSend(tag, kROOTD_AUTH);
         } else
            // We will receive the user token next
            NetSend(1, kROOTD_AUTH);
      } else {
         // No authentication available for re-use
         NetSend(0, kROOTD_AUTH);
         // Cleanup and return: we need a new one ...
         if (tkn) delete[] tkn;
         if (user) delete[] user;
         // ... no need to continue receiving the old token
         return retval;
      }
   }

   // Now Receive Token
   int ofs = *OffSet;
   char *token = 0;
   if (gRSAKey > 0) {
      if (RpdSecureRecv(&token) == -1) {
         ErrorInfo("RpdCheckAuthTab: problems secure-"
                   "receiving token %s",
                   "- may result in authentication failure ");
      }

   } else {
      EMessageTypes kind;
      int lenToken = 9;
      token = new char[lenToken];
      NetRecv(token, lenToken, kind);
      if (kind != kMESS_STRING)
         ErrorInfo
             ("RpdCheckAuthTab: got msg kind: %d instead of %d (kMESS_STRING)",
              kind, kMESS_STRING);
      // Invert Token
      for (int i = 0; i < (int) strlen(token); i++) {
         token[i] = ~token[i];
      }
   }
   if (gDebug > 2)
      ErrorInfo
          ("RpdCheckAuthTab: received from client: token: '%s' ",
           token);

   // Check tag, if there
   if (token && strlen(token) > 8) {
      // Create hex from tag
      char tagref[9] = {0};
      SPrintf(tagref,9,"%08x",tag);
      if (strncmp(token+8,tagref,8)) {
         ErrorInfo("RpdCheckAuthTab: token tag does not match - failure");
         goodOfs = 0;
      } else
         // Drop tag
         token[8] = 0;
   }

   // Now check Token validity
   if (goodOfs && token && RpdCheckToken(token, tkn)) {
      retval = 1;
      // Comunicate new offset to remote client
      *OffSet = ofs;
   }

   if (tkn) delete[] tkn;
   if (token) delete[] token;
   if (user) delete[] user;

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Check offset received from client entry in tab file.

int RpdCheckOffSet(int Sec, const char *User, const char *Host, int RemId,
                   int *OffSet, char **Token, int *ShmId, char **GlbsUser)
{
   int retval = 0;
   bool goodOfs = 0;
   int ofs = *OffSet >= 0 ? *OffSet : 0;

   if (gDebug > 2)
      ErrorInfo("RpdCheckOffSet: analyzing: %d %s %s %d %d", Sec, User,
                Host, RemId, *OffSet);

   // Open file
   int itab = open(gRpdAuthTab.c_str(), O_RDWR);
   if (itab == -1) {
      if (GetErrno() == ENOENT)
         ErrorInfo("RpcCheckOffSet: file %s does not exist",
                   gRpdAuthTab.c_str());
      else
         ErrorInfo("RpcCheckOffSet: error opening %s (errno: %d)",
                   gRpdAuthTab.c_str(), GetErrno());
      return retval;
   }
   // lock tab file
   if (lockf(itab, F_LOCK, (off_t) 1) == -1) {
      ErrorInfo("RpcCheckOffSet: error locking %s (errno: %d)",
                gRpdAuthTab.c_str(), GetErrno());
      close(itab);
      return retval;
   }
   if (gDebug > 0)
      ErrorInfo("RpdCheckOffSet: file LOCKED");

   // File is open: set position at wanted location
   if (lseek(itab, ofs, SEEK_SET) < 0) {
      ErrorInfo("RpcCheckOffSet: error lseeking %s (errno: %d)",
                gRpdAuthTab.c_str(), GetErrno());
      close(itab);
      return retval;
   }

   // Now read the entry
   char line[kMAXPATHLEN];
   if (reads(itab,line, sizeof(line)) < 0) {
      ErrorInfo("RpcCheckOffSet: error reading %d bytes from %s (errno: %d)",
                sizeof(line), gRpdAuthTab.c_str(), GetErrno());
      close(itab);
      return retval;
   }

   // and parse its content according to auth method
   int lsec, act, remid, shmid = -1;
   char host[kMAXPATHLEN], usr[kMAXPATHLEN], subj[kMAXPATHLEN],
       dumm[kMAXPATHLEN], tkn[20];
   int nw =
       sscanf(line, "%d %d %d %d %4095s %4095s %19s %4095s",
                    &lsec, &act, &gRSAKey, &remid, host, usr, tkn, dumm);
   if (gDebug > 2)
      ErrorInfo("RpdCheckOffSet: found line: %s", line);

   if (nw > 5 && act > 0) {
      if (lsec == Sec) {
         if (lsec == 3) {
            sscanf(line, "%d %d %d %d %4095s %4095s %d %4095s %19s %4095s",
                         &lsec, &act, &gRSAKey, &remid, host, usr, &shmid, subj, tkn, dumm);
            if ((remid == RemId)
                && !strcmp(host, Host) && !strcmp(subj, User))
               goodOfs = 1;
         } else {
            if ((remid == RemId) &&
                !strcmp(host, Host) && !strcmp(usr, User))
               goodOfs = 1;
         }
      }
   }
   if (!goodOfs) {
      // Tab may have been cleaned in the meantime ... try a scan
      lseek(itab, 0, SEEK_SET);
      ofs = 0;
      while (reads(itab, line, sizeof(line))) {

         nw = sscanf(line, "%d %d %d %d %4095s %4095s %19s %4095s",
                     &lsec, &act, &gRSAKey, &remid, host, usr, tkn, dumm);
         if (gDebug > 2)
            ErrorInfo("RpdCheckOffSet: found line: %s", line);

         if (nw > 5 && act > 0) {
            if (lsec == Sec) {
               if (lsec == 3) {
                  sscanf(line, "%d %d %d %d %4095s %4095s %d %4095s %19s %4095s",
                         &lsec, &act, &gRSAKey, &remid, host, usr, &shmid, subj, tkn, dumm);
                  if ((remid == RemId)
                      && !strcmp(host, Host) && !strcmp(subj, User)) {
                     goodOfs = 1;
                     goto found;
                  }
               } else {
                  if ((remid == RemId) &&
                      !strcmp(host, Host) && !strcmp(usr, User)) {
                     goodOfs = 1;
                     goto found;
                  }
               }
            }
         }
      }
   }

   found:
   // unlock the file
   lseek(itab, 0, SEEK_SET);
   if (lockf(itab, F_ULOCK, (off_t) 1) == -1) {
      ErrorInfo("RpcCheckOffSet: error unlocking %s",
                gRpdAuthTab.c_str());
   }
   // closing file ...
   close(itab);

   // Read public key
   std::string pukfile = gRpdKeyRoot;
   pukfile.append(ItoA(*OffSet));
   if (gDebug > 2)
      ErrorInfo("RpdCheckOffSet: RSAKey ofs file: %d %d '%s' ",
                gRSAKey, ofs, pukfile.c_str());

   struct passwd *pw = getpwnam(usr);
   if (pw) {
      uid_t fromUid = getuid();
      uid_t fromEUid = geteuid();
      // The check must be done with 'usr' UIDs to prevent
      // unauthorized users from forcing the server to read
      // manipulated key files.
      if (fromUid == 0)
         if (setresuid(pw->pw_uid, pw->pw_uid, fromEUid) == -1)
            // Since we could not set the user IDs, we will
            // not trust the client
            goodOfs = 0;

      // Get the key now
      if (goodOfs)
         if (RpdGetRSAKeys(pukfile.c_str(), 1) < 1)
            goodOfs = 0;

      // Reset original IDs
      if (getuid() != fromUid)
         if (setresuid(fromUid, fromEUid, pw->pw_uid) == -1)
            goodOfs = 0;

   } else {
      // Since we could not set the user IDs, we will
      // not trust the client
      goodOfs = 0;
      if (gDebug > 0)
         ErrorInfo("RpdCheckOffSet: error in getpwname(%s) (errno: %d)",
                   usr,GetErrno());
   }

   if (gDebug > 2)
      ErrorInfo("RpdCheckOffSet: goodOfs: %d (active: %d)",
                goodOfs, act);

   // Comunicate new offset to remote client
   if (goodOfs) {

      // Rename the key file, if needed
      if (*OffSet > 0 && *OffSet != ofs) {
         if (RpdRenameKeyFile(*OffSet,ofs) > 0) {
            goodOfs = 0;
            // Error: set entry inactive
            RpdCleanupAuthTab(Host,RemId,ofs);
         }
      }

      *OffSet = ofs;
      // return token if requested
      if (Token) {
         const size_t tokenSize = strlen(tkn)+1;
         *Token = new char[tokenSize];
         strlcpy(*Token,tkn,tokenSize);
      }
      if (Sec == 3) {
         if (GlbsUser) {
            const size_t glbsUserSize = strlen(usr)+1;
            *GlbsUser = new char[glbsUserSize];
            strlcpy(*GlbsUser,usr,glbsUserSize);
         }
         if (ShmId)
            *ShmId = shmid;
      }
   }

   return goodOfs;
}

////////////////////////////////////////////////////////////////////////////////
/// Rename public file with new offset
/// Returns: 0 if OK
///          1 if problems renaming

int RpdRenameKeyFile(int oldofs, int newofs)
{
   int retval = 0;

   // Old name
   std::string oldname = gRpdKeyRoot;
   oldname.append(ItoA(oldofs));
   // New name
   std::string newname = gRpdKeyRoot;
   newname.append(ItoA(newofs));

   if (rename(oldname.c_str(), newname.c_str()) == -1) {
      if (gDebug > 0)
         ErrorInfo("RpdRenameKeyFile: error renaming key file"
                   " %s to %s (errno: %d)",
                   oldname.c_str(),newname.c_str(),GetErrno());
      retval = 2;
   }

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Check token validity.

bool RpdCheckToken(char *token, char *tknref)
{
   // Get rid of '\n'
   char *s = strchr(token, '\n');
   if (s)
      *s = 0;
   s = strchr(tknref, '\n');
   if (s)
      *s = 0;

#ifndef R__NOCRYPT
   char *tkn_crypt = crypt(token, tknref);
   int tlen = 13;
#else
   char *tkn_crypt = rpdcrypt(token,gRndmSalt.c_str());
   int tlen = 16;
#endif

   if (gDebug > 2)
      ErrorInfo("RpdCheckToken: ref:'%s' crypt:'%s'", tknref, tkn_crypt);

   if (!strncmp(tkn_crypt, tknref, tlen))
      return 1;
   else
      return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the requiring subject has already authenticated during this session
/// and its 'ticket' is still valid.

int RpdReUseAuth(const char *sstr, int kind)
{
   int lenU, offset, opt;
   gOffSet = -1;
   gExistingAuth = 0;
   int auth= 0;

   if (gDebug > 2)
      ErrorInfo("RpdReUseAuth: analyzing: %s, %d", sstr, kind);

   char user[64];

   // kClear
   if (kind == kROOTD_USER) {
      if (!(gReUseAllow & gAUTH_CLR_MSK)) {
         return 0;              // re-authentication required by administrator
      }
      gSec = 0;
      // Decode subject string
      sscanf(sstr, "%d %d %d %d %63s", &gRemPid, &offset, &opt, &lenU, user);
      user[lenU] = '\0';
      if ((gReUseRequired = (opt & kAUTH_REUSE_MSK))) {
         gOffSet = offset;
         if (gRemPid > 0 && gOffSet > -1) {
            auth =
                RpdCheckAuthTab(gSec, user, gOpenHost.c_str(), gRemPid, &gOffSet);
         }
         if ((auth == 1) && (offset != gOffSet))
            auth = 2;
         // Fill gUser and free allocated memory
         strlcpy(gUser, user, sizeof(gUser));
      }
   }

   // Flag if existing token has been re-used
   if (auth > 0)
      gExistingAuth = 1;

   // Return value
   return auth;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if required auth method is allowed for 'Host'.
/// If 'yes', returns 0, if 'no', returns 1, the number of allowed
/// methods in NumAllow, and the codes of the allowed methods (in order
/// of preference) in AllowMeth. Memory for AllowMeth must be allocated
/// outside. Directives read from (in decreasing order of priority):
/// $ROOTDAEMONRC, $HOME/.rootdaemonrc (privately startd daemons only)
/// or $ROOTETCDIR/system.rootdaemonrc.

int RpdCheckAuthAllow(int Sec, const char *Host)
{
   int retval = 1, found = 0;

   std::string theDaemonRc;

   // Check if a non-standard file has been requested
   if (getenv("ROOTDAEMONRC"))
      theDaemonRc = getenv("ROOTDAEMONRC");

   if (theDaemonRc.length() <= 0) {
      if (getuid()) {
         // Check if user has a private daemon access file ...
         struct passwd *pw = getpwuid(getuid());
         if (pw != 0) {
            theDaemonRc = std::string(pw->pw_dir).append("/");
            theDaemonRc.append(gDaemonRc);
         } else {
            if (getenv("ROOTETCDIR")) {
               theDaemonRc = std::string(getenv("ROOTETCDIR")).append("/system");
               theDaemonRc.append(gDaemonRc);
            } else
               theDaemonRc = std::string("/etc/root/system").append(gDaemonRc);
         }
      } else {
         // If running as super-user, check system file only
         if (getenv("ROOTETCDIR")) {
            theDaemonRc = std::string(getenv("ROOTETCDIR")).append("/system");
            theDaemonRc.append(gDaemonRc);
         } else
            theDaemonRc = std::string("/etc/root/system").append(gDaemonRc);
      }
   }
   if (gDebug > 2)
      ErrorInfo("RpdCheckAuthAllow: Checking file: %s for meth:%d"
                " host:%s (gNumAllow: %d)",
                theDaemonRc.c_str(), Sec, Host, gNumAllow);

   // Check if info already loaded (not first call ...)
   if (gMethInit == 1) {

      // Look for the method in the allowed list and flag this method
      // as tried, if found ...
      int newtry = 0, i;
      for (i = 0; i < gNumAllow; i++) {
         if (gTriedMeth[i] == 0 && gAllowMeth[i] == Sec) {
            newtry = 1;
            gTriedMeth[i] = 1;
            gNumLeft--;
         }
      }
      if (newtry == 0) {
         ErrorInfo
             ("RpdCheckAuthAllow: new auth method proposed by %s",
              " client not in the list or already attempted");
         return retval;
      }
      retval = 0;

   } else {

      // Open file
      FILE *ftab = fopen(theDaemonRc.c_str(), "r");
      if (ftab == 0) {
         if (GetErrno() == ENOENT)
            ErrorInfo("RpdCheckAuthAllow: file %s does not exist",
                      theDaemonRc.c_str());
         else
            ErrorInfo("RpdCheckAuthAllow: error opening %s (errno: %d)",
                      theDaemonRc.c_str(), GetErrno());
      }
      // Now read the entry
      char line[kMAXPATHLEN], host[kMAXPATHLEN], rest[kMAXPATHLEN],
          cmth[kMAXPATHLEN];
      int nmet = 0, mth[6] = { 0 };

      int cont = 0, jm = -1;
      while (ftab && fgets(line, sizeof(line), ftab)) {
         int i;
         if (line[0] == '#')
            continue;           // skip comment lines
         if (line[strlen(line) - 1] == '\n')
            line[strlen(line) - 1] = '\0';   // get rid of '\n', if any ...
         // Analyze the line now ...
         int nw = 0;
         char *pstr = line;
         // Check if a continuation line
         if (cont == 1) {
            cont = 0;
            strlcpy(rest, pstr, kMAXPATHLEN);
         } else {
            jm = -1;
            // Get 'host' first ...
            nw = sscanf(pstr, "%4095s %4095s", host, rest);
            if (nw < 2)
               continue;        // no method defined for this host
            pstr = line + strlen(host) + 1;

            // Check if a service is specified
            char *pcol = strstr(host, ":");
            if (pcol) {
               if (!strstr(pcol+1, gServName[gService].c_str()))
                  continue;
               else
                  host[(int)(pcol-host)] = '\0';
            }
            if (!host[0])
               strlcpy(host, "default", kMAXPATHLEN);

            if (gDebug > 2)
               ErrorInfo("RpdCheckAuthAllow: found host: %s ", host);

            if (strcmp(host, "default")) {
               // now check validity of 'host' format
               if (!RpdCheckHost(Host,host)) {
                  goto next;
               }
            } else {
               // This is a default entry: ignore it if a host-specific entry was already
               // found, analyse it otherwise ...
               if (found == 1)
                  goto next;
            }

            // Reset mth[kMAXSEC]
            nmet = 0;
            for (i = 0; i < kMAXSEC; i++) {
               mth[i] = -1;
            }

         }

         // We are at the end and there will be a continuation line ...
         if (rest[0] == '\\') {
            cont = 1;
            continue;
         }

         while (pstr != 0) {
            int tmet = -1;
            char *pd = 0, *pd2 = 0;
            cmth[0] = '\0';
            rest[0] = '\0';
            nw = sscanf(pstr, "%4095s %4095s", cmth, rest);
            if (!strcmp(cmth, "none")) {
               nmet = 0;
               goto nexti;
            }
            pd = strchr(cmth, ':');
            // Parse the method
            char tmp[20];
            if (pd != 0) {
               int mlen = pd - cmth;
               strncpy(tmp, cmth, mlen);
               tmp[mlen] = '\0';
            } else {
               strlcpy(tmp, cmth, sizeof(tmp));
            }

            if (strlen(tmp) > 1) {

               for (tmet = 0; tmet < kMAXSEC; tmet++) {
                  if (!rpdstrcasecmp(gAuthMeth[tmet].c_str(), tmp))
                     break;
               }
               if (tmet < kMAXSEC) {
                  if (gDebug > 2)
                     ErrorInfo("RpdCheckAuthAllow: tmet %d", tmet);
               } else {
                  if (gDebug > 1)
                     ErrorInfo("RpdCheckAuthAllow: unknown methods"
                               " %s - ignore", tmp);
                  goto nexti;
               }

            } else {
               tmet = atoi(tmp);
            }
            jm = -1;
            if (gDebug > 2)
               ErrorInfo("RpdCheckAuthAllow: found method %d (have?:%d)",
                         tmet, (tmet >= 0 && tmet < kMAXSEC) ? gHaveMeth[tmet] : 0);
            if (tmet >= 0 && tmet < kMAXSEC) {
               if (gHaveMeth[tmet] == 1) {
                  int ii;
                  for (ii = 0; ii < nmet; ii++) {
                     if (mth[ii] == tmet) {
                        jm = ii;
                     }
                  }
               } else
                  goto nexti;
            } else
               goto nexti;
            if (jm == -1) {
               // New method ...
               mth[nmet] = tmet;
               jm = nmet;
               nmet++;
            }
            // Now parse users list, if any ...
            while (pd != 0 && (int) (pd[1]) != 32) {
               pd2 = strchr(pd + 1, ':');
               if (pd[1] == '-') {
                  pd += 2;
                  // Ignore
                  if (gUserIgnore[mth[jm]] == 0) {
                     gUserIgnLen[mth[jm]] = kMAXPATHLEN;
                     gUserIgnore[mth[jm]] = new char[gUserIgnLen[mth[jm]]];
                     gUserIgnore[mth[jm]][0] = '\0';
                  }
                  if (strlen(gUserIgnore[mth[jm]]) >
                      (gUserIgnLen[mth[jm]] - 10)) {
                     char *tmpUI = strdup(gUserIgnore[mth[jm]]);
                     free(gUserIgnore[mth[jm]]);
                     gUserIgnLen[mth[jm]] += kMAXPATHLEN;
                     gUserIgnore[mth[jm]] = new char[gUserIgnLen[mth[jm]]];
                     strlcpy(gUserIgnore[mth[jm]], tmpUI, sizeof(gUserIgnLen[mth[jm]]));
                     free(tmpUI);
                  }
                  char usr[256];
                  if (pd2 != 0) {
                     int ulen = pd2 - pd;
                     strncpy(usr, pd, ulen);
                     usr[ulen] = '\0';
                  } else {
                     strlcpy(usr, pd, sizeof(usr));
                  }
                  struct passwd *pw = getpwnam(usr);
                  if (pw != 0)
                     SPrintf(gUserIgnore[mth[jm]], gUserIgnLen[mth[jm]], "%s %d",
                             gUserIgnore[mth[jm]], (int)pw->pw_uid);
               } else {
                  pd += 1;
                  if (pd[1] == '+')
                     pd += 1;
                  // Keep
                  if (gUserAllow[mth[jm]] == 0) {
                     gUserAlwLen[mth[jm]] = kMAXPATHLEN;
                     gUserAllow[mth[jm]] = new char[gUserAlwLen[mth[jm]]];
                     gUserAllow[mth[jm]][0] = '\0';
                  }
                  if (strlen(gUserAllow[mth[jm]]) >
                      (gUserAlwLen[mth[jm]] - 10)) {
                     char *tmpUI = strdup(gUserAllow[mth[jm]]);
                     free(gUserAllow[mth[jm]]);
                     gUserAlwLen[mth[jm]] += kMAXPATHLEN;
                     gUserAllow[mth[jm]] = new char[gUserAlwLen[mth[jm]]];
                     strlcpy(gUserAllow[mth[jm]], tmpUI, sizeof(gUserAlwLen[mth[jm]]));
                     free(tmpUI);
                  }
                  char usr[256];
                  if (pd2 != 0) {
                     int ulen = pd2 - pd;
                     strncpy(usr, pd, ulen);
                     usr[ulen] = '\0';
                  } else {
                     strlcpy(usr, pd, sizeof(usr));
                  }
                  struct passwd *pw = getpwnam(usr);
                  if (pw != 0)
                     SPrintf(gUserAllow[mth[jm]], gUserAlwLen[mth[jm]], "%s %d",
                             gUserAllow[mth[jm]], (int)pw->pw_uid);
               }
               pd = pd2;
            }
            // Get next item
            nexti:
            if (nw > 1 && (int) rest[0] != 92) {
               pstr = strstr(pstr, rest);
            } else {
               if ((int) rest[0] == 92)
                  cont = 1;
               pstr = 0;
            }
         }
         if (gDebug > 2) {
            ErrorInfo("RpdCheckAuthAllow: for host %s found %d methods",
                      host, nmet);
            ErrorInfo("RpdCheckAuthAllow: %d %d %d %d %d %d", mth[0],
                      mth[1], mth[2], mth[3], mth[4], mth[5]);
         }
         // Found new entry matching: superseed previous result
         found = 1;
         retval = 1;
         gNumAllow = gNumLeft = nmet;
         for (i = 0; i < kMAXSEC; i++) {
            gAllowMeth[i] = -1;
            gTriedMeth[i] = 0;
            if (i < gNumAllow) {
               gAllowMeth[i] = mth[i];
               if (Sec == mth[i]) {
                  retval = 0;
                  gNumLeft--;
                  gTriedMeth[i] = 1;
               }
            }
         }
         next:
         continue;
      }

      // closing file ...
      if (ftab)
         fclose(ftab);

      // Host specific directives have been checked for ...
      gMethInit = 1;

      // Use defaults if nothing found
      if (!found) {
         if (gDebug > 2)
            ErrorInfo
            ("RpdCheckAuthAllow: no specific or 'default' entry found: %s",
             "using system defaults");
         int i;
         for (i = 0; i < gNumAllow; i++) {
            if (Sec == gAllowMeth[i]) {
               retval = 0;
               gNumLeft--;
               gTriedMeth[i] = 1;
            }
         }

      }

   }
   if (gDebug > 2) {
      ErrorInfo
          ("RpdCheckAuthAllow: returning: %d (gNumAllow: %d, gNumLeft:%d)",
           retval, gNumAllow, gNumLeft);
      int i, jm;
      for (i = 0; i < kMAXSEC; i++) {
         jm = gAllowMeth[i];
         if (gUserAlwLen[jm] > 0)
            ErrorInfo("RpdCheckAuthAllow: users allowed for method %d: %s",
                      jm, gUserAllow[jm]);
      }
      for (i = 0; i < kMAXSEC; i++) {
         jm = gAllowMeth[i];
         if (gUserIgnLen[jm] > 0)
            ErrorInfo("RpdCheckAuthAllow: users ignored for method %d: %s",
                      jm, gUserIgnore[jm]);
      }
   }

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if 'host' is compatible with 'Host' taking into account
/// wild cards in the host name
/// Returns 1 if successful, 0 otherwise ...

int RpdCheckHost(const char *Host, const char *host)
{
   int rc = 1;

   // Strings must be both defined
   if (!Host || !host)
      return 0;

   // If host is a just wild card accept it
   if (!strcmp(host,"*"))
      return 1;

   // Try now to understand whether it is an address or a name ...
   int name = 0, i = 0;
   for (i = 0; i < (int) strlen(host); i++) {
      if ((host[i] < 48 || host[i] > 57) &&
           host[i] != '*' && host[i] != '.') {
         name = 1;
         break;
      }
   }

   // If ref host is an IP, get IP of Host
   char *hh;
   if (!name) {
      hh = RpdGetIP(Host);
      if (gDebug > 2)
         ErrorInfo("RpdCheckHost: Checking Host IP: %s", hh);
   } else {
      const size_t hhSize = strlen(Host)+1;
      hh = new char[hhSize];
      strlcpy(hh,Host,hhSize);
      if (gDebug > 2)
         ErrorInfo("RpdCheckHost: Checking Host name: %s", hh);
   }

   // Check if starts with wild
   // Starting with '.' defines a full (sub)domain
   int sos = 0;
   if (host[0] == '*' || host[0] == '.')
      sos = 1;

   // Check if ends with wild
   // Ending with '.' defines a name
   int eos = 0, le = strlen(host);
   if (host[le-1] == '*' || host[le-1] == '.')
      eos = 1;

   int first= 1;
   int ends= 0;
   int starts= 0;
   const size_t hSize = strlen(host)+1;
   char *h = new char[hSize];
   strlcpy(h,host,hSize);
   char *tk = strtok(h,"*");
   while (tk) {

      char *ps = strstr(hh,tk);
      if (!ps) {
         rc = 0;
         break;
      }
      if (!sos && first && ps == hh)
         starts = 1;
      first = 0;

      if (ps == hh + strlen(hh) - strlen(tk))
         ends = 1;

      tk = strtok(0,"*");

   }
   delete[] h;
   delete[] hh;

   if ((!sos || !eos) && !starts && !ends)
      rc = 0;

   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Get IP address of 'host' as a string. String must be deleted by
/// the user.

char *RpdGetIP(const char *host)
{
   struct hostent *h;
   unsigned long ip;
   unsigned char ip_fld[4];

   // Check server name
   if ((h = gethostbyname(host)) == 0) {
      ErrorInfo("RpdGetIP: unknown host %s", host);
      return 0;
   }
   // Decode ...
   ip = ntohl(*(unsigned long *) h->h_addr_list[0]);
   ip_fld[0] = (unsigned char) ((0xFF000000 & ip) >> 24);
   ip_fld[1] = (unsigned char) ((0x00FF0000 & ip) >> 16);
   ip_fld[2] = (unsigned char) ((0x0000FF00 & ip) >> 8);
   ip_fld[3] = (unsigned char) ((0x000000FF & ip));

   // Prepare output
   char *output = new char[20];
   SPrintf(output, 20, "%d.%d.%d.%d",
           ip_fld[0], ip_fld[1], ip_fld[2], ip_fld[3]);

   // return
   return output;
}

////////////////////////////////////////////////////////////////////////////////
/// Send list of authentication methods not yet tried.

void RpdSendAuthList()
{
   if (gDebug > 2)
      ErrorInfo("RpdSendAuthList: analyzing (gNumLeft: %d)", gNumLeft);

   // Send Number of methods left
   NetSend(gNumLeft, kROOTD_NEGOTIA);

   if (gNumLeft > 0) {
      int i = 0;
      std::string alist;
      char cm[5];
      for (i = 0; i < gNumAllow; i++) {
         if (gDebug > 2)
            ErrorInfo("RpdSendAuthList: gTriedMeth[%d]: %d", i,
                      gTriedMeth[i]);
         if (gTriedMeth[i] == 0) {
            SPrintf(cm, 5, " %d",gAllowMeth[i]);
            alist.append(cm);
         }
      }
      NetSend(alist.c_str(), alist.length() + 1, kMESS_STRING);
      if (gDebug > 2)
         ErrorInfo("RpdSendAuthList: sent list: %s", alist.c_str());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Authenticate via Kerberos.

int RpdKrb5Auth(const char *)
{
   ::Error("RpdKrb5Auth", "Kerberos5 no longer supported by ROOT");
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Secure Remote Password protocol (no longer supported)

int RpdSRPUser(const char *)
{
   ::Error("RpdSRPUser", "SRP no longer supported by ROOT");
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the requesting {host,user} can be granted immediate
/// login on the base of the information found in /etc/hosts.equiv
/// and/or $HOME/.rhosts. The two files must be trustable, i.e. owned
/// and modifiable only by 'root' and by 'user', respectively (0600).
/// Returns 1 in case access can be granted, 0 in any other case
/// (errout contains a code for error logging on the client side)
///
/// NB: entries granting access in one of the two files cannot be
///     overriden in the other file; so, system admins cannot close
///     access from a host and user cannot stop access to their
///     account if the administrator has decided so; as an example,
///     if this entry is found in /etc/hosts.equiv
///
///     remote.host.dom auser
///
///     (allowing user named 'auser' from host 'remote.host.dom' to
///     login to any non-root local account without specifying a
///     password) the following entries in $home/.rhosts are ignored
///
///     remote.host.dom -auser
///     -remote.host.dom
///
///     and access to 'auser' is always granted. This is a "feature"
///     of ruserok.
///

int RpdCheckHostsEquiv(const char *host, const char *ruser,
                       const char *user, int &errout)
{
   int rc = 0;

   // Effective uid
   int rootuser = 0;
   if (!geteuid() && !getegid())
      rootuser = 1;

   // Check the files only if i) at least one exists; ii) those existing
   // have the right permission settings
   bool badfiles = 0;
   int  nfiles = 0;

   // Check system file /etc/hosts.equiv if non-root
   char hostsequiv[20] = { "/etc/hosts.equiv" };
   if (!rootuser) {

      // Get info about the file ...
      struct stat st;
      if (stat(hostsequiv,&st) == -1) {
         if (GetErrno() != ENOENT) {
            ErrorInfo("RpdCheckHostsEquiv: cannot stat /etc/hosts.equiv"
                      " (errno: %d)",GetErrno());
            badfiles = 1;
         } else
            if (gDebug > 1)
               ErrorInfo("RpdCheckHostsEquiv: %s does not exist",
                         hostsequiv);
      } else {

         // Require 'root' ownership
         if (st.st_uid || st.st_gid) {
            if (gDebug > 0)
               ErrorInfo("RpdCheckHostsEquiv: /etc/hosts.equiv not owned by"
                         " system (uid: %d, gid: %d)",st.st_uid,st.st_gid);
            badfiles = 1;
         } else {

            // Require WRITE permission only for owner
            if ((st.st_mode & S_IWGRP) || (st.st_mode & S_IWOTH)) {
               if (gDebug > 0)
                  ErrorInfo("RpdCheckHostsEquiv: group or others have write"
                            " permission on /etc/hosts.equiv: do not trust"
                            " it (g: %d, o: %d)",
                            (st.st_mode & S_IWGRP),(st.st_mode & S_IWOTH));
               badfiles = 1;
            } else
               // Good file
               nfiles++;
         }
      }
   }

   // Check local file
   char rhosts[kMAXPATHLEN] = {0};
   if (!badfiles) {

      struct passwd *pw = getpwnam(user);
      if (pw) {
         int ldir = strlen(pw->pw_dir);
         ldir = (ldir > kMAXPATHLEN - 9) ? (kMAXPATHLEN - 9) : ldir;
         memcpy(rhosts,pw->pw_dir,ldir);
         memcpy(rhosts+ldir,"/.rhosts",8);
         rhosts[ldir+8] = 0;
         if (gDebug > 2)
            ErrorInfo("RpdCheckHostsEquiv: checking for user file %s ...",rhosts);
      } else {
         if (gDebug > 0)
            ErrorInfo("RpdCheckHostsEquiv: cannot get user info with getpwnam"
                   " (errno: %d)",GetErrno());
         badfiles = 1;
      }

      if (!badfiles) {
         // Check the $HOME/.rhosts file ... ownership and protections
         struct stat st;
         if (stat(rhosts,&st) == -1) {
            if (GetErrno() != ENOENT) {
               ErrorInfo("RpdCheckHostsEquiv: cannot stat $HOME/.rhosts"
                      " (errno: %d)",GetErrno());
               badfiles = 1;
            } else
               ErrorInfo("RpdCheckHostsEquiv: %s/.rhosts does not exist",
                         pw->pw_dir);
         } else {

            // Only use file when its access rights are 0600
            if (!S_ISREG(st.st_mode) || S_ISDIR(st.st_mode) ||
                (st.st_mode & 0777) != (S_IRUSR | S_IWUSR)) {
               if (gDebug > 0)
                  ErrorInfo("RpdCheckHostsEquiv: unsecure permission setting"
                            " found for $HOME/.rhosts: 0%o (must be 0600)",
                            (st.st_mode & 0777));
               badfiles = 1;
            } else
               // Good file
               nfiles++;
         }
      }
   }

   // if files are not available or have wrong permissions or are
   // not accessible, give up
   if (!nfiles) {
      if (gDebug > 0)
         ErrorInfo("RpdCheckHostsEquiv: no files to check");
      errout = 1;
      if (badfiles) {
         if (gDebug > 0)
            ErrorInfo("RpdCheckHostsEquiv: config files cannot be used"
                      " (check permissions)");
         errout = 2;
      }
      return rc;
   }

   // Ok, now use ruserok to find out if {host,ruser,user}
   // is trusted
#if defined(__sgi) || defined(_AIX)
   if (ruserok((char*)host,rootuser,(char*)ruser,(char*)user) == 0) {
#else
   if (ruserok(host,rootuser,ruser,user) == 0) {
#endif
      if (gDebug > 0)
         ErrorInfo("RpdCheckHostsEquiv: remote user %s authorized to"
                   " access %s's area",ruser,user);
      rc = 1;
   } else {
      if (gDebug > 0)
         ErrorInfo("RpdCheckHostsEquiv: no special permission from"
                   " %s or %s",hostsequiv,rhosts);
      errout = 3;
   }

   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Check received user's password against password in $HOME/.rootdpass.
/// The password is retrieved in RpdUser and temporarly saved in gPasswd.
/// Returns 1 in case of success authentication, 0 otherwise.

int RpdCheckSpecialPass(const char *passwd)
{
   // Check inputs
   if (!passwd)
      return 0;

   // and the saved the password
   if (strlen(gPasswd) <= 0)
      return 0;

   // Ok, point to the saved passwd (retrieved in RpdUser)
   char *rootdpass = gPasswd;
   int n = 0;

   if (gClientProtocol > 8 && gSaltRequired > 0) {
      n = strlen(rootdpass);
      if (strncmp(passwd, rootdpass, n + 1) != 0) {
         if (gDebug > 0)
            ErrorInfo("RpdCheckSpecialPass: wrong password");
         rpdmemset((volatile void *)rootdpass,0,n);
         return 0;
      }
   } else {
#ifndef R__NOCRYPT
      char *pass_crypt = crypt(passwd, rootdpass);
#else
      char *pass_crypt = (char *)passwd;
#endif
      n = strlen(rootdpass);
      if (strncmp(pass_crypt, rootdpass, n+1) != 0) {
         if (gDebug > 0)
            ErrorInfo("RpdCheckSpecialPass: wrong password");
         rpdmemset((volatile void *)rootdpass,0,n);
         return 0;
      }
   }

   if (gDebug > 0)
      ErrorInfo
          ("RpdCheckSpecialPass: user %s authenticated via ~/.rootdpass",
           gUser);

   rpdmemset((volatile void *)rootdpass,0,n);
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Check user's password.

int RpdPass(const char *pass, int errheq)
{
   char passwd[128];
   char *passw;
   char *pass_crypt;
   struct passwd *pw;
#ifdef R__SHADOWPW
   struct spwd *spw;
#endif
   int afs_auth = 0;

   if (gDebug > 2)
      ErrorInfo("RpdPass: Enter (pass length: %d)", (int)strlen(pass));

   int auth = 0;
   errheq = (errheq > -1 && errheq < 4) ? errheq : 0;
   if (!*gUser) {
      if (gClientProtocol > 11)
         NetSend(gUsrPwdErr[0][errheq], kROOTD_ERR);
      else
         NetSend(kErrFatal, kROOTD_ERR);
      if (gDebug > 0)
         ErrorInfo("RpdPass: user needs to be specified first");
      return auth;
   }

   if (!pass) {
      if (gClientProtocol > 11)
         NetSend(gUsrPwdErr[1][errheq], kROOTD_ERR);
      else
         NetSend(kErrNoPasswd, kROOTD_ERR);
      if (gDebug > 0)
         ErrorInfo("RpdPass: no password specified");
      return auth;
   }
   int n = strlen(pass);
   // Passwd length should be in the correct range ...
   if (!n) {
      if (gClientProtocol > 11)
         NetSend(gUsrPwdErr[1][errheq], kROOTD_ERR);
      else
         NetSend(kErrBadPasswd, kROOTD_ERR);
      if (gDebug > 0)
         ErrorInfo("RpdPass: null passwd not allowed");
      return auth;
   }
   if (n > (int) sizeof(passwd)) {
      if (gClientProtocol > 11)
         NetSend(gUsrPwdErr[1][errheq], kROOTD_ERR);
      else
         NetSend(kErrBadPasswd, kROOTD_ERR);
      if (gDebug > 0)
         ErrorInfo("RpdPass: passwd too long");
      return auth;
   }
   // Inversion is done in RpdUser, if needed
   strlcpy(passwd, pass, sizeof(passwd));

   // Special treatment for anonimous ...
   if (gAnon) {
      strlcpy(gPasswd, passwd, sizeof(gPasswd));
      goto authok;
   }
   // ... and SpecialPass ...
   if (RpdCheckSpecialPass(passwd)) {
      goto authok;
   }
   // Get local passwd info for gUser
   if (!(pw = getpwnam(gUser))) {
      ErrorInfo("RpdPass: getpwnam failed!");
      return auth;
   }

#ifdef R__SHADOWPW
      // System V Rel 4 style shadow passwords
      if ((spw = getspnam(gUser)) == 0) {
         if (gDebug > 0)
            ErrorInfo("RpdPass: Shadow passwd not available for user %s",
                   gUser);
         passw = pw->pw_passwd;
      } else
         passw = spw->sp_pwdp;
#else
      passw = pw->pw_passwd;
#endif
#ifndef R__NOCRYPT
      if (gClientProtocol <= 8 || !gSaltRequired) {
         char salt[20] = {0};
         int lenS = 2;
         if (!strncmp(passw, "$1$", 3)) {
            // Shadow passwd
            char *pd = strstr(passw + 4, "$");
            lenS = (int) (pd - passw);
            strncpy(salt, passw, lenS);
         } else
            strncpy(salt, passw, lenS);
         salt[lenS] = 0;
         pass_crypt = crypt(passwd, salt);   // Comment this
      } else {
         pass_crypt = passwd;
      }
#else
      pass_crypt = passwd;
#endif
      n = strlen(passw);
      if (strncmp(pass_crypt, passw, n + 1) != 0) {
         if (gClientProtocol > 11)
            NetSend(gUsrPwdErr[1][errheq], kROOTD_ERR);
         else
            NetSend(kErrBadPasswd, kROOTD_ERR);
         if (gDebug > 0)
            ErrorInfo("RpdPass: invalid password for user %s", gUser);
         return auth;
      }
      if (gDebug > 2)
         ErrorInfo("RpdPass: valid password for user %s", gUser);

   authok:
   auth = afs_auth ? 5 : 1;
   gSec = 0;

   if (gClientProtocol > 8) {
      // Set an entry in the auth tab file for later (re)use, if required ...
      int offset = -1;
      char *token = 0;
      char line[kMAXPATHLEN];
      if ((gReUseAllow & gAUTH_CLR_MSK) && gReUseRequired) {

         SPrintf(line, kMAXPATHLEN, "0 1 %d %d %s %s",
                 gRSAKey, gRemPid, gOpenHost.c_str(), gUser);
         if (!afs_auth || gService == kPROOFD)
            offset = RpdUpdateAuthTab(1, line, &token);
         if (gDebug > 2)
            ErrorInfo("RpdPass: got offset %d", offset);

         // Comunicate login user name to client
         SPrintf(line, kMAXPATHLEN, "%s %d", gUser, offset);
         if (gDebug > 2)
            ErrorInfo("RpdPass: sending back line %s", line);
         NetSend(strlen(line), kROOTD_PASS);   // Send message length first
         NetSend(line, kMESS_STRING);

         if (offset > -1) {
            if (gDebug > 2)
               ErrorInfo("RpdPass: sending token %s (Crypt: %d)", token,
                         gCryptRequired);
            if (gCryptRequired) {
               // Send over the token
               if (RpdSecureSend(token) == -1) {
                  if (gDebug > 0)
                     ErrorInfo("RpdPass: problems secure-sending token"
                            " - may result in corrupted token");
               }
            } else {
               // Send token inverted
               for (int i = 0; i < (int) strlen(token); i++) {
                  token[i] = ~token[i];
               }
               NetSend(token, kMESS_STRING);
            }
            delete[] token;
         }
         gOffSet = offset;

      } else {
         // Comunicate login user name to client
         SPrintf(line, kMAXPATHLEN, "%s -1", gUser);
         if (gDebug > 2)
            ErrorInfo("RpdPass: sending back line %s", line);
         NetSend(strlen(line), kROOTD_PASS);   // Send message length first
         NetSend(line, kMESS_STRING);
      }
   }

   return auth;
}

////////////////////////////////////////////////////////////////////////////////
/// Prepare for globus authentication: check hostcer.conf and get
/// the credential handle. This is run once at daemon start-up

int RpdGlobusInit()
{
   ::Error("RpdGlobusInit", "Globus is no longer supported by ROOT");
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Authenticate via Globus.

int RpdGlobusAuth(const char *)
{
   ::Error("RpdGlobusInit", "Globus is no longer supported by ROOT");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// RFIO protocol (no longer supported by ROOT)

int RpdRfioAuth(const char *)
{
   ::Error("RpdRfioAuth", "RfioAuth no longer supported by ROOT");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Terminate correctly by cleaning up the auth table and closing the file.
/// Called upon receipt of a kROOTD_CLEANUP and on SIGPIPE.

void RpdAuthCleanup(const char *sstr, int opt)
{
   int rpid = 0, sec = -1, offs = -1, nw = 0;
   char usr[64] = {0};
   if (sstr)
      nw = sscanf(sstr, "%d %d %d %63s", &rpid, &sec, &offs, usr);

   // Turn back to superuser for cleaning, if the case
   if (getuid() == 0) {
      if (setresgid(0, 0, 0) == -1)
         if (gDebug > 0)
            ErrorInfo("RpdAuthCleanup: can't setgid to superuser");
      if (setresuid(0, 0, 0) == -1)
         if (gDebug > 0)
            ErrorInfo("RpdAuthCleanup: can't setuid to superuser");
   }
   if (opt == 0) {
      RpdCleanupAuthTab("all", 0, -1);            // Cleanup everything (SIGPIPE)
      ErrorInfo("RpdAuthCleanup: cleanup ('all',0) done");
   } else if (opt == 1) {
      if (nw == 1) {
         // host specific cleanup
         RpdCleanupAuthTab(gOpenHost.c_str(), rpid, -1);
         ErrorInfo("RpdAuthCleanup: cleanup ('%s',%d) done",
                   gOpenHost.c_str(), rpid);
      } else if (nw == 4) {
         // (host,usr,method) specific cleanup
         if (RpdCheckOffSet(sec,usr,gOpenHost.c_str(),rpid,&offs,0,0,0)) {
            RpdCleanupAuthTab(gOpenHost.c_str(), rpid, offs);
            ErrorInfo("RpdAuthCleanup: cleanup (%s,%d,%d,%d,%s) done",
                      gOpenHost.c_str(), rpid, sec, offs, usr);
         } else {
            ErrorInfo("RpdAuthCleanup: cleanup not done: %s",
                      "wrong offset or already cleaned up");
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

void RpdInitAuth()
{
   // Size check done in RpdUpdateAuthTab(1,...)

   // Reset
   int i;
   gNumAllow = gNumLeft = 0;
   for (i = 0; i < kMAXSEC; i++) {
      gAllowMeth[i] = -1;
      gHaveMeth[i] = 1;
   }

   // List of default authentication methods
   RpdDefaultAuthAllow();
}

////////////////////////////////////////////////////////////////////////////////
/// Check configuration options and running daemons to build a default list
/// of secure methods.

void RpdDefaultAuthAllow()
{
   if (gDebug > 2)
      ErrorInfo("RpdDefaultAuthAllow: Enter");

   // UsrPwdClear
   gAllowMeth[gNumAllow] = 0;
   gNumAllow++;
   gNumLeft++;

   // No SRP method
   gHaveMeth[1] = 0;

   // No Kerberos method
   gHaveMeth[2] = 0;

   // No Globus method
   gHaveMeth[3] = 0;

   if (gDebug > 2) {
      int i;
      std::string temp;
      char cm[5];
      if (gNumAllow == 0)
         temp.append("none");
      for (i = 0; i < gNumAllow; i++) {
         SPrintf(cm, 5, " %3d",gAllowMeth[i]);
         temp.append(cm);
      }
      ErrorInfo
          ("RpdDefaultAuthAllow: default list of secure methods available: %s",
           temp.c_str());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check the running of process 'daemon'.
/// Info got from 'ps ax'.

int RpdCheckDaemon(const char *daemon)
{
   char cmd[kMAXPATHLEN] = { 0 };
   int ch, i = 0, cnt = 0;

   if (gDebug > 2)
      ErrorInfo("RpdCheckDaemon: Enter ... %s", daemon);

   // Return if empty
   if (daemon == 0 || !daemon[0])
      return cnt;

   // Build command
   SPrintf(cmd, kMAXPATHLEN, "ps ax | grep %s 2>/dev/null", daemon);

   // Run it ...
   FILE *fp = popen(cmd, "r");
   if (fp != 0) {
      for (ch = fgetc(fp); ch != EOF; ch = fgetc(fp)) {
         if (ch != 10) {
            cmd[i++] = ch;
         } else {
            cmd[i] = '\0';
            if (strstr(cmd, "grep") == 0 && strstr(cmd, "rootd") == 0
                && strstr(cmd, "proofd") == 0) {
               cnt++;
            }
            i = 0;
         }
      }
      if (i > 0) {
         cmd[i] = '\0';
         cnt++;
      }
      pclose(fp);
      if (gDebug > 2)
         ErrorInfo("RpdCheckDaemon: found %d instances of daemon %s",
                   cnt, daemon);

   } else {
      ErrorInfo("RpdCheckDaemon: problems executing cmd ...");
   }
   return cnt;
}

////////////////////////////////////////////////////////////////////////////////
/// Check user id. If user id is not equal to rootd's effective uid, user
/// will not be allowed access, unless effective uid = 0 (i.e. root).

int RpdUser(const char *sstr)
{
   const int kMaxBuf = 256;
   char recvbuf[kMaxBuf];
   EMessageTypes kind;
   struct passwd *pw;
   if (gDebug > 2)
      ErrorInfo("RpdUser: Enter ... %s", sstr);

   int auth = 0;

   // Nothing can be done if empty message
   if (!*sstr) {
      NetSend(kErrBadUser, kROOTD_ERR);
      ErrorInfo("RpdUser: received empty string");
      return auth;
   }
   // Parse input message
   char user[kMAXUSERLEN] = {0};
   char ruser[kMAXUSERLEN] = {0};
   if (gClientProtocol > 8) {
      int ulen, ofs, opt, rulen;
      // Decode subject string
      int nw = sscanf(sstr, "%d %d %d %d %63s %d %63s",
                            &gRemPid, &ofs, &opt, &ulen, user, &rulen, ruser);
      ulen = (ulen >= kMAXUSERLEN) ? kMAXUSERLEN-1 : ulen;
      rulen = (rulen >= kMAXUSERLEN) ? kMAXUSERLEN-1 : rulen;
      user[ulen] = '\0';
      if (nw > 5)
         ruser[rulen] = '\0';
      gReUseRequired = (opt & kAUTH_REUSE_MSK);
      gCryptRequired = (opt & kAUTH_CRYPT_MSK);
      gSaltRequired  = (opt & kAUTH_SSALT_MSK);
      gOffSet = ofs;
#ifdef R__SSL
      if (gRSASSLKey) {
         // Determine type of RSA key required
         gRSAKey = (opt & kAUTH_RSATY_MSK) ? 2 : 1;
      } else
         gRSAKey = 1;
#else
      gRSAKey = 1;
#endif
   } else {
      SPrintf(user,kMAXUSERLEN,"%s",sstr);
   }
   if (gDebug > 2)
      ErrorInfo("RpdUser: gReUseRequired: %d gCryptRequired: %d gRSAKey: %d",
                gReUseRequired, gCryptRequired, gRSAKey);

   ERootdErrors err = kErrNoUser;
   if (gService == kROOTD) {
      // Default anonymous account ...
      if (!strcmp(user, "anonymous")) {
         user[0] = '\0';
         strlcpy(user, "rootd", sizeof(user));
      }
   }

   if ((pw = getpwnam(user)) == 0) {
      NetSend(err, kROOTD_ERR);
      ErrorInfo("RpdUser: user %s unknown", user);
      return auth;
   }

   // If server is not started as root and user is not same as the
   // one who started rootd then authetication is not ok.
   uid_t uid = getuid();
   if (uid && uid != pw->pw_uid) {
      NetSend(kErrBadUser, kROOTD_ERR);
      ErrorInfo("RpdUser: user not same as effective user of rootd");
      return auth;
   }

   // Check if the administrator allows authentication
   char cuid[20];
   SPrintf(cuid, 20, "%d", (int)pw->pw_uid);
   if (gUserIgnLen[0] > 0 && strstr(gUserIgnore[0], cuid) != 0) {
      NetSend(kErrNotAllowed, kROOTD_ERR);
      ErrorInfo
          ("RpdUser: user (%d,%s) not authorized to use UsrPwd method",
           uid, pw->pw_name);
      return auth;
   }
   if (gUserAlwLen[0] > 0 && strstr(gUserAllow[0], cuid) == 0) {
      NetSend(kErrNotAllowed, kROOTD_ERR);
      ErrorInfo
          ("RpdUser: user (%d,%s) not authorized to use UsrPwd method",
           uid, pw->pw_name);
      return auth;
   }

   // Check /etc/hosts.equiv and/or $HOME/.rhosts
   int errheq = 0;
   if (gCheckHostsEquiv && strlen(ruser)) {
      if (RpdCheckHostsEquiv(gOpenHost.c_str(),ruser,user,errheq)) {
         auth = 3;
         strlcpy(gUser, user, sizeof(gUser));
         return auth;
      }
   }

   // Check if of type anonymous ...
   if (!strcmp(pw->pw_shell, "/bin/false")) {
      err = kErrNoAnon;
      gAnon = 1;
      gReUseRequired = 0;
   }

   // Check if authorized
   // If not anonymous, try to get passwd
   // (if our system uses shadow passwds and we are not superuser
   // we cannot authenticate users ...)
   //   char *passw = 0;
   gPasswd[0] = 0;
   char *passw = gPasswd;
   int errrdp = 0;
   if (gAnon == 0) {

      // Check ROOT specific passwd first
      int rcsp = RpdRetrieveSpecialPass(user,gRootdPass.c_str(),
                                        gPasswd,sizeof(gPasswd));
      if (rcsp < 0)
         errrdp = (rcsp == -2) ? 3 : 0;

      if (!passw[0] || !strcmp(passw, "x")) {

#ifdef R__SHADOWPW
         struct spwd *spw = 0;
         // System V Rel 4 style shadow passwords
         if ((spw = getspnam(user)) == 0) {
            if (gDebug > 0) {
               ErrorInfo("RpdUser: Shadow passwd not accessible for user %s",user);
               ErrorInfo("RpdUser: trying normal system passwd");
            }
         } else
            passw = spw->sp_pwdp;
#else
         passw = pw->pw_passwd;
#endif
         // Check if successful
         if (!passw[0] || !strcmp(passw, "x")) {
            if (gClientProtocol > 11)
               NetSend(gUsrPwdErr[errrdp][errheq], kROOTD_ERR);
            else
               NetSend(kErrNotAllowed, kROOTD_ERR);
            ErrorInfo("RpdUser: passwd hash not available for user %s", user);
            ErrorInfo
                ("RpdUser: user %s cannot be authenticated with this method",
                 user);
            return auth;
         }
      }
   }
   // Ok: Save username and go to next steps
   strlcpy(gUser, user, sizeof(gUser));

   // Salt vars
   char salt[30] = { 0 };
   char ctag[11] = { 0 };
   int  rtag = 0;
   int lenS = 0;

   if (gClientProtocol > 8) {

      // Prepare status flag to send back
      if (gAnon == 1) {
         // Anonymous user: we will receive a text pass in the form
         // user@remote.host.dom
         NetSend(-1, kROOTD_AUTH);

      } else {

         if (gCryptRequired) {
            // Named user: first we receive a session public key
            // Ask for the RSA key
            NetSend(gRSAKey, kROOTD_RSAKEY);

            // Receive the key securely
            if (RpdRecvClientRSAKey()) {
               ErrorInfo("RpdUser: could not import a valid key -"
                         " switch off reuse for this session");
               gReUseRequired = 0;
            }

            // We get a random tag
            if (gClientProtocol > 11) {
               RpdInitRand();
               rtag = rpd_rand();
               SPrintf(ctag, 11, "#%08x#",rtag);
            }

            if (gSaltRequired) {
               // The crypt man page says that alternative salts can be in the form '$1$...$';
               // but on Ubuntu 10.04 the salt are in the form '$j$...$' where j is 6 or other.
               if (passw[0] == '$' && passw[2] == '$') {
                  // Shadow passwd
                  char *pd = strstr(passw + 4, "$");
                  lenS = (int) (pd - passw);
                  strncpy(salt, passw, lenS);
                  salt[lenS] = 0;
               } else {
                  lenS = 2;
                  strncpy(salt, passw, lenS);
                  salt[lenS] = 0;
               }

               if (gDebug > 2)
                  ErrorInfo("RpdUser: salt: '%s' ",salt);

               // We add the random tag here
               if (gClientProtocol > 11) {
                  strncat(salt, ctag, sizeof(salt) - sizeof(ctag) - 1);
                  salt[sizeof(salt) - 1] = '\0';
               }

               // Send it over encrypted
               if (RpdSecureSend(salt) == -1) {
                  ErrorInfo("RpdUser: problems secure-sending salt -"
                            " may result in corrupted salt");
               }
            } else {
               if (gClientProtocol > 11) {
                  // We send the random tag here
                  if (RpdSecureSend(ctag) == -1) {
                     ErrorInfo("RpdUser: problems secure-sending rndmtag -"
                               " may result in corrupted rndmtag");
                  }
               } else
                  NetSend(0, kMESS_ANY);
            }
         } else {
            // We continue the authentication process in clear
            NetSend(0, kROOTD_AUTH);
         }
      }

   } else {
      // If we are talking to a old client protocol
      NetSend(0, kROOTD_AUTH);
   }

   // Get the password hash or anonymous string
   if (NetRecv(recvbuf, kMaxBuf, kind) < 0) {
      NetSend(kErrFatal, kROOTD_ERR);
      ErrorInfo("RpdUser: error receiving message");
      return auth;
   }
   if (kind != kROOTD_PASS) {
      NetSend(kErrFatal, kROOTD_ERR);
      ErrorInfo("RpdUser: received wrong message type: %d (expecting: %d)",
                kind, (int) kROOTD_PASS);
      return auth;
   }
   if (!strncmp(recvbuf,"-1",2)) {
      if (gDebug > 0)
         ErrorInfo("RpdUser: client did not send a password - return");
      return auth;
   }
   // Get passwd
   char *passwd = 0;
   int lpwd = 0;
   if (gAnon == 0 && gClientProtocol > 8 && gCryptRequired) {

      // Receive encrypted pass or its hash
      if (RpdSecureRecv(&passwd) == -1) {
         ErrorInfo
             ("RpdUser: problems secure-receiving pass hash - %s",
              "may result in authentication failure");
      }
      // Length of the password buffer
      lpwd = strlen(passwd);

      // Check the random tag, if any
      if (strlen(ctag)) {

         // Check first that there is enough space for the tag
         int plen = lpwd;
         if (plen > 9 &&
             passwd[plen-1] == '#' && passwd[plen-10] == '#') {
            if (strncmp(ctag,&passwd[plen-10],10)) {
               // The tag does not match; failure
               if (gClientProtocol > 11)
                  NetSend(gUsrPwdErr[2][errheq], kROOTD_ERR);
               else
                  NetSend(kErrBadPasswd, kROOTD_ERR);
               ErrorInfo("RpdUser: rndm tag mis-match"
                         " (%s vs %s) - Failure",&passwd[plen-10],ctag);
               delete[] passwd;
               return auth;
            }

            // Tag ok: drop it
            plen -= 10;
            passwd[plen] = 0;

         } else {
            // The tag is not there or incomplete; failure
            if (gClientProtocol > 11)
               NetSend(gUsrPwdErr[2][errheq], kROOTD_ERR);
            else
               NetSend(kErrBadPasswd, kROOTD_ERR);
            ErrorInfo("RpdUser: rndm tag missing or incomplete"
                      " (pw length: %d) - Failure", plen);
            delete[] passwd;
            return auth;
         }
      }

      // If we required an hash check that we got it
      // (the client sends the passwd if the crypt version is different)
      if (gSaltRequired && lenS) {
         if (strncmp(passwd,salt,lenS))
            gSaltRequired = 0;
      }

   } else {

      // Receive clear or anonymous pass
      passwd = new char[strlen(recvbuf) + 1];

      // Re-invert pass
      int i, n = strlen(recvbuf);
      for (i = 0; i < n; i++)
         passwd[i] = ~recvbuf[i];
      passwd[i] = '\0';

      if (gDebug > 2 && gAnon)
         ErrorInfo("RpdUser: received anonymous pass: '%s'", passwd);
   }

   // Check the passwd and login if ok ...
   auth = RpdPass(passwd,errheq);

   // Erase memory used for password
   passwd = (char *)rpdmemset((volatile void *)passwd,0,lpwd);
   delete[] passwd;

   return auth;
}

////////////////////////////////////////////////////////////////////////////////
/// Try a guess of the client protocol from what they sent over
/// the net ...

int RpdGuessClientProt(const char *buf, EMessageTypes kind)
{
   if (gDebug > 2)
      ErrorInfo("RpdGuessClientProt: Enter: buf: '%s', kind: %d", buf,
                (int) kind);

   // Assume same version as us.
   int proto = 9;

   // Clear authentication
   if (kind == kROOTD_USER) {
      char usr[64], rest[256];
      int ns = sscanf(buf, "%63s %255s", usr, rest);
      if (ns == 1)
         proto = 8;
   }

   if (gDebug > 2)
      ErrorInfo("RpdGuessClientProt: guess for gClientProtocol is %d",
                proto);

   // Return the guess
   return proto;
}

////////////////////////////////////////////////////////////////////////////////
/// Allocates and Fills a NULL terminated buffer of length Len+1 with
/// Len random characters.
/// Return pointer to the buffer (to be deleted by the caller)
/// Opt = 0      any non dangerous char
///       1      letters and numbers  (upper and lower case)
///       2      hex characters       (upper and lower case)
///       3      crypt like           [a-zA-Z0-9./]

char *RpdGetRandString(int Opt, int Len)
{
   unsigned int iimx[4][4] = {
      { 0x0, 0xffffff08, 0xafffffff, 0x2ffffffe }, // Opt = 0
      { 0x0, 0x3ff0000,  0x7fffffe,  0x7fffffe },  // Opt = 1
      { 0x0, 0x3ff0000,  0x7e,       0x7e },       // Opt = 2
      { 0x0, 0x3ffc000,  0x7fffffe,  0x7fffffe }   // Opt = 3
   };

   const char *cOpt[4] = { "Any", "LetNum", "Hex", "Crypt" };

   //  Default option 0
   if (Opt < 0 || Opt > 3) {
      Opt = 0;
      if (gDebug > 2)
         ErrorInfo("RpdGetRandString: Unknown option: %d : assume 0", Opt);
   }
   if (gDebug > 2)
      ErrorInfo("RpdGetRandString: Enter ... Len: %d %s", Len, cOpt[Opt]);

   // Allocate buffer
   char *buf = new char[Len + 1];

   // Init Random machinery ...
   if (!gRandInit)
      RpdInitRand();

   // randomize
   int k = 0;
   int i, j, l, m, frnd;
   while (k < Len) {
      frnd = rpd_rand();
      for (m = 7; m < 32; m += 7) {
         i = 0x7F & (frnd >> m);
         j = i / 32;
         l = i - j * 32;
         if ((iimx[Opt][j] & (1 << l))) {
            buf[k] = i;
            k++;
         }
         if (k == Len)
            break;
      }
   }

   // NULL terminated
   buf[Len] = 0;
   if (gDebug > 2)
      ErrorInfo("RpdGetRandString: got '%s' ", buf);

   return buf;
}

////////////////////////////////////////////////////////////////////////////////
/// Get public key from file pubkey (Opt == 1) or string pubkey (Opt == 0).

int RpdGetRSAKeys(const char *pubkey, int Opt)
{
   char str[kMAXPATHLEN] = { 0 };
   int keytype = 0;

   if (gDebug > 2)
      ErrorInfo("RpdGetRSAKeys: enter: string len: %d, opt %d ",
                 gPubKeyLen, Opt);

   if (!pubkey)
      return keytype;

   char *theKey = 0;
   FILE *fKey = 0;
   // Parse input type
   if (Opt == 1) {

      // Ok, now open it
      fKey = fopen(pubkey, "r");
      if (!fKey) {
         if (GetErrno() == EACCES) {
            struct passwd *pw = getpwuid(getuid());
            char *usr = 0;
            if (pw)
               usr = pw->pw_name;
            ErrorInfo("RpdGetRSAKeys: access to key file %s denied"
                      " to user: %s", pubkey, (usr ? usr : (char *)"????"));
         } else
            ErrorInfo("RpdGetRSAKeys: cannot open key file"
                      " %s (errno: %d)", pubkey, GetErrno());
         return 0;
      }
      // Check first the permissions: should be 0600
      struct stat st;
      if (fstat(fileno(fKey), &st) == -1) {
         ErrorInfo("RpdGetRSAKeys: cannot stat descriptor %d"
                   " %s (errno: %d)", fileno(fKey), GetErrno());
         fclose(fKey);
         return 0;
      }
      if (!S_ISREG(st.st_mode) || S_ISDIR(st.st_mode) ||
          (st.st_mode & 0777) != (S_IRUSR | S_IWUSR)) {
         ErrorInfo("RpdGetRSAKeys: key file %s: wrong permissions"
                   " 0%o (should be 0600)", pubkey, (st.st_mode & 0777));
         fclose(fKey);
         return 0;
      }
      gPubKeyLen = fread((void *)str,1,sizeof(str),fKey);
      if (gDebug > 2)
         ErrorInfo("RpdGetRSAKeys: length of the read key: %d",gPubKeyLen);

      // This the key
      theKey = str;
   } else {
      // the key is the argument
      theKey = (char *)pubkey;
   }

   if (gPubKeyLen > 0) {

      // Skip spaces at beginning, if any
      int k = 0;
      while (theKey[k] == 32) k++;

      keytype = gRSAKey;

      // The format of keytype 1 is #<hex_n>#<hex_d>#
      char *pd1 = 0, *pd2 = 0, *pd3 = 0;
      pd1 = strstr(theKey, "#");
      if (pd1) pd2 = strstr(pd1 + 1, "#");
      if (pd2) pd3 = strstr(pd2 + 1, "#");
      if (keytype == 1) {
         if (!pd1 || !pd2 || !pd3) {
            if (gDebug > 0)
               ErrorInfo("RpdGetRSAKeys: bad format for keytype %d"
                         " - exit", keytype);
            keytype = 0;
         }
      }
      if (keytype == 1) {

         if (gDebug > 2)
            ErrorInfo("RpdGetRSAKeys: keytype %d ", keytype);

         // Get <hex_n> ...
         int l1 = (int) (pd2 - pd1 - 1);
         char *n_exp_RSA = new char[l1 + 1];
         strncpy(n_exp_RSA, pd1 + 1, l1);
         n_exp_RSA[l1] = 0;
         if (gDebug > 2)
            ErrorInfo("RpdGetRSAKeys: got %d bytes for n_exp_RSA",
                      strlen(n_exp_RSA));
         // Now <hex_d>
         int l2 = (int) (pd3 - pd2 - 1);
         char *d_exp_RSA = new char[l2 + 1];
         strncpy(d_exp_RSA, pd2 + 1, l2);
         d_exp_RSA[l2] = 0;
         if (gDebug > 2)
            ErrorInfo("RpdGetRSAKeys: got %d bytes for d_exp_RSA",
                      strlen(d_exp_RSA));

         rsa_num_sget(&gRSA_n, n_exp_RSA);
         rsa_num_sget(&gRSA_d, d_exp_RSA);

         delete[] n_exp_RSA;
         delete[] d_exp_RSA;

      } else if (keytype == 2){

#ifdef R__SSL
         // try SSL
         if (gDebug > 2)
            ErrorInfo("RpdGetRSAKeys: keytype %d ", keytype);

         // Now set the key locally in BF form
         BF_set_key(&gBFKey, gPubKeyLen, (const unsigned char *)theKey);
#else
         if (gDebug > 0) {
            ErrorInfo("RpdGetRSAKeys: not compiled with SSL support:"
                   " you should not have got here!");
         }
#endif
      }
   }

   if (fKey)
      fclose(fKey);

   return keytype;
}

////////////////////////////////////////////////////////////////////////////////
/// Save RSA public key into file for later use by other rootd/proofd.
/// Return: 0 if ok
///         1 if not ok
///         2 if not ok because file already exists and cannot be
///           overwritten

int RpdSavePubKey(const char *PubKey, int OffSet, char *user)
{
   int retval = 0;

   if (gRSAKey == 0 || OffSet < 0)
      return 1;

   std::string pukfile = gRpdKeyRoot;
   pukfile.append(ItoA(OffSet));

   // Unlink the file first
   if (unlink(pukfile.c_str()) == -1) {
      if (GetErrno() != ENOENT)
         // File exists and cannot overwritten by this process
         return 2;
   }

   // Create file
   int ipuk = -1;
   ipuk = open(pukfile.c_str(), O_WRONLY | O_CREAT, 0600);
   if (ipuk == -1) {
      ErrorInfo("RpdSavePubKey: cannot open file %s (errno: %d)",
                pukfile.c_str(),GetErrno());
      if (GetErrno() == ENOENT)
         return 2;
      else
         return 1;
   }

   // If root process set ownership of the pub key to the user
   if (getuid() == 0) {
      struct passwd *pw = getpwnam(user);
      if (pw) {
         if (fchown(ipuk,pw->pw_uid,pw->pw_gid) == -1) {
            ErrorInfo("RpdSavePubKey: cannot change ownership"
                      " of %s (errno: %d)",pukfile.c_str(),GetErrno());
            retval = 1;
         }
      } else {
         ErrorInfo("RpdSavePubKey: getpwnam failure (errno: %d)",GetErrno());
         retval = 1;
      }
   }

   // Write the key if no error occured
   if (retval == 0) {
      while (write(ipuk, PubKey, gPubKeyLen) < 0 && GetErrno() == EINTR)
         ResetErrno();
   }

   // close the file
   close(ipuk);

   // Over
   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Encode null terminated str using the session private key indcated by Key
/// and sends it over the network.
/// Returns number of bytes sent.or -1 in case of error.

int RpdSecureSend(char *str)
{
   char buftmp[kMAXSECBUF];
   char buflen[20];

   int ttmp = 0;
   int nsen = -1;

   if (gRSAKey == 1) {
      strncpy(buftmp, str, sizeof(buftmp) - 1);
      buftmp[kMAXSECBUF - 1] = '\0';
      ttmp = rsa_encode(buftmp, strlen(buftmp) + 1, gRSA_n, gRSA_d);
   } else if (gRSAKey == 2) {
#ifdef R__SSL
      ttmp = strlen(str);
      if ((ttmp % 8) > 0)            // It should be a multiple of 8!
         ttmp = ((ttmp + 8)/8) * 8;
      unsigned char iv[8];
      memset((void *)&iv[0],0,8);
      BF_cbc_encrypt((const unsigned char *)str, (unsigned char *)buftmp,
                     strlen(str), &gBFKey, iv, BF_ENCRYPT);
#else
      ErrorInfo("RpdSecureSend: Not compiled with SSL support:"
                " you should not have got here! - return");
#endif
   } else {
      ErrorInfo("RpdSecureSend: Unknown key option (%d) - return",
                gRSAKey);
   }

   // Send the buffer now
   SPrintf(buflen, 20, "%d", ttmp);
   NetSend(buflen, kROOTD_ENCRYPT);
   nsen = NetSendRaw(buftmp, ttmp);
   if (gDebug > 4)
      ErrorInfo("RpdSecureSend: sent %d bytes (expected: %d) - keytype: %d",
                 nsen, ttmp, gRSAKey);

   return nsen;
}

////////////////////////////////////////////////////////////////////////////////
/// Receive buffer and decode it in str using key indicated by Key type.
/// Return number of received bytes or -1 in case of error.

int RpdSecureRecv(char **str)
{
   char buftmp[kMAXSECBUF];
   char buflen[20];

   int nrec = -1;
   // We must get a pointer ...
   if (!str)
      return nrec;

   if (gDebug > 2)
      ErrorInfo("RpdSecureRecv: enter ... (key is %d)", gRSAKey);

   EMessageTypes kind;
   NetRecv(buflen, 20, kind);
   int len = atoi(buflen);
   if (gDebug > 4)
      ErrorInfo("RpdSecureRecv: got len '%s' %d ", buflen, len);
   if (!strncmp(buflen, "-1", 2))
      return nrec;

   // receive the buffer
   nrec = NetRecvRaw(buftmp,len);

   // decode it
   if (gRSAKey == 1) {
      rsa_decode(buftmp, len, gRSA_n, gRSA_d);
      if (gDebug > 2)
         ErrorInfo("RpdSecureRecv: Local: decoded string is %d bytes long",
                   strlen(buftmp));

      // Prepare output
      const size_t strSize = strlen(buftmp) + 1;
      *str = new char[strSize];
      strlcpy(*str, buftmp, strSize);
   } else if (gRSAKey == 2) {
#ifdef R__SSL
      unsigned char iv[8];
      memset((void *)&iv[0],0,8);
      *str = new char[nrec + 1];
      BF_cbc_encrypt((const unsigned char *)buftmp, (unsigned char *)(*str),
                      nrec, &gBFKey, iv, BF_DECRYPT);
      (*str)[nrec] = '\0';
#else
      ErrorInfo("RpdSecureRecv: Not compiled with SSL support:"
                " you should not have got here! - return");
#endif
   } else {
      ErrorInfo("RpdSecureRecv: Unknown key option (%d) - return",
                gRSAKey);
   }

   return nrec;
}

////////////////////////////////////////////////////////////////////////////////
/// Generate a valid pair of private/public RSA keys to protect for
/// authentication password and token exchange
/// Returns 1 if a good key pair is not found after kMAXRSATRIES attempts
/// Returns 0 if a good key pair is found
/// If setrndinit = 1, no futher init of the random engine

int RpdGenRSAKeys(int setrndinit)
{
   if (gDebug > 2)
      ErrorInfo("RpdGenRSAKeys: enter");

   // Init Random machinery ...
   if (!gRandInit)
      RpdInitRand();
   gRandInit = setrndinit;

#ifdef R__NOCRYPT
   // Generate a random salt
   char *rsalt = RpdGetRandString(3,8);
   if (rsalt) {
      gRndmSalt = std::string(rsalt);
      delete[] rsalt;
   } else {
      if (gDebug > 0)
         ErrorInfo("RpdGenRSAKeys: could not generate random salt");
   }
#endif

#ifdef R__SSL
   // Generate also the SSL key
   if (gDebug > 2)
      ErrorInfo("RpdGenRSAKeys: Generate RSA SSL keys");

   // Init SSL ...
   SSL_library_init();

   //  ... and its error strings
   SSL_load_error_strings();

   // Load Ciphers
   OpenSSL_add_all_ciphers();

   // Number of bits for key
   Int_t nbits = 1024;

   // Public exponent
   Int_t pubex = 17;

   // Init random engine
   char *rbuf = RpdGetRandString(0,40);
   RAND_seed(rbuf,strlen(rbuf));

   // Generate Key
   gRSASSLKey = RSA_generate_key(nbits,pubex,0,0);

   // Bio for exporting the pub key
   BIO *bkey = BIO_new(BIO_s_mem());

   // Write public key to BIO
   PEM_write_bio_RSAPublicKey(bkey,gRSASSLKey);

   // Read key from BIO to buf
   Int_t sbuf = 2*RSA_size(gRSASSLKey);
   char *kbuf = new char[sbuf];
   BIO_read(bkey,(void *)kbuf,sbuf);
   BIO_free(bkey);

   // Prepare export
   gRSAPubExport[1].len = sbuf;
   gRSAPubExport[1].keys = new char[gRSAPubExport[1].len + 2];
   strncpy(gRSAPubExport[1].keys,kbuf,gRSAPubExport[1].len);
   gRSAPubExport[1].keys[gRSAPubExport[1].len-1] = '\0';
   delete[] kbuf;
   if (gDebug > 2)
      ErrorInfo("RpdGenRSAKeys: SSL: export pub:\n%.*s",
           gRSAPubExport[1].len,gRSAPubExport[1].keys);

   // We have at least one key
   gRSAInit = 1;

#endif

   // Sometimes some bunch is not decrypted correctly
   // That's why we make retries to make sure that encryption/decryption
   // works as expected
   bool notOK = 1;
   rsa_NUMBER p1, p2, rsa_n, rsa_e, rsa_d;
   int l_n = 0, l_d = 0;
#if R__RSADEB
   Int_t l_e = 0;
   char buf[rsa_STRLEN];
#endif
   char buf_n[rsa_STRLEN], buf_e[rsa_STRLEN], buf_d[rsa_STRLEN];

   int nAttempts = 0;
   int thePrimeLen = kPRIMELENGTH;
   int thePrimeExp = kPRIMEEXP + 5;   // Prime probability = 1-0.5^thePrimeExp
   while (notOK && nAttempts < kMAXRSATRIES) {

      nAttempts++;
      if (gDebug > 2 && nAttempts > 1) {
            ErrorInfo("RpdGenRSAKeys: retry no. %d",nAttempts);
         srand(rpd_rand());
      }

      // Valid pair of primes
      p1 = rsa_genprim(thePrimeLen, thePrimeExp);
      p2 = rsa_genprim(thePrimeLen+1, thePrimeExp);

      // Retry if equal
      int nPrimes = 0;
      while (rsa_cmp(&p1, &p2) == 0 && nPrimes < kMAXRSATRIES) {
         nPrimes++;
         if (gDebug > 2)
            ErrorInfo("RpdGenRSAKeys: equal primes: regenerate (%d times)",nPrimes);
         srand(rpd_rand());
         p1 = rsa_genprim(thePrimeLen, thePrimeExp);
         p2 = rsa_genprim(thePrimeLen+1, thePrimeExp);
      }

#if R__RSADEB
      if (gDebug > 2) {
         rsa_num_sput(&p1, buf, rsa_STRLEN);
         ErrorInfo("RpdGenRSAKeys: local: p1: '%s' ", buf);
         rsa_num_sput(&p2, buf, rsa_STRLEN);
         ErrorInfo("RpdGenRSAKeys: local: p2: '%s' ", buf);
      }
#endif

      // Generate keys
      if (rsa_genrsa(p1, p2, &rsa_n, &rsa_e, &rsa_d)) {
         if (gDebug > 0)
            ErrorInfo("RpdGenRSAKeys: genrsa: attempt %d to generate"
                      " keys failed",nAttempts);
         continue;
      }

      // Determine their lengths
      rsa_num_sput(&rsa_n, buf_n, rsa_STRLEN);
      l_n = strlen(buf_n);
      rsa_num_sput(&rsa_e, buf_e, rsa_STRLEN);
#if R__RSADEB
      l_e = strlen(buf_e);
#endif
      rsa_num_sput(&rsa_d, buf_d, rsa_STRLEN);
      l_d = strlen(buf_d);

#if R__RSADEB
      if (gDebug > 2) {
         ErrorInfo("RpdGenRSAKeys: local: n: '%s' length: %d", buf_n, l_n);
         ErrorInfo("RpdGenRSAKeys: local: e: '%s' length: %d", buf_e, l_e);
         ErrorInfo("RpdGenRSAKeys: local: d: '%s' length: %d", buf_d, l_d);
      }
#endif
      if (rsa_cmp(&rsa_n, &rsa_e) <= 0)
         continue;
      if (rsa_cmp(&rsa_n, &rsa_d) <= 0)
         continue;

      // Now we try the keys
      char test[2 * rsa_STRLEN] = "ThisIsTheStringTest01203456-+/";
      Int_t lTes = 31;
      char *dumT = RpdGetRandString(0, lTes - 1);
      strncpy(test, dumT, lTes);
      delete[]dumT;
      char buf[2 * rsa_STRLEN];
      if (gDebug > 3)
         ErrorInfo("RpdGenRSAKeys: local: test string: '%s' ", test);

      // Private/Public
      strncpy(buf, test, lTes);
      buf[lTes] = 0;

      // Try encryption with private key
      int lout = rsa_encode(buf, lTes, rsa_n, rsa_e);
      if (gDebug > 3)
         ErrorInfo("GenRSAKeys: local: length of crypted string: %d bytes", lout);

      // Try decryption with public key
      rsa_decode(buf, lout, rsa_n, rsa_d);
      buf[lTes] = 0;
      if (gDebug > 3)
         ErrorInfo("RpdGenRSAKeys: local: after private/public : '%s' ", buf);

      if (strncmp(test, buf, lTes))
         continue;

      // Public/Private
      strncpy(buf, test, lTes);
      buf[lTes] = 0;

      // Try encryption with public key
      lout = rsa_encode(buf, lTes, rsa_n, rsa_d);
      if (gDebug > 3)
         ErrorInfo("RpdGenRSAKeys: local: length of crypted string: %d bytes ",
              lout);

      // Try decryption with private key
      rsa_decode(buf, lout, rsa_n, rsa_e);
      buf[lTes] = 0;
      if (gDebug > 3)
         ErrorInfo("RpdGenRSAKeys: local: after public/private : '%s' ", buf);

      if (strncmp(test, buf, lTes))
         continue;

      notOK = 0;
   }

   if (notOK) {
      ErrorInfo("RpdGenRSAKeys: unable to generate good RSA key pair"
                " (%d attempts)- return",kMAXRSATRIES);
      return 1;
   }

   // Save Private key
   rsa_assign(&gRSAPriKey.n, &rsa_n);
   rsa_assign(&gRSAPriKey.e, &rsa_e);

   // Save Public key
   rsa_assign(&gRSAPubKey.n, &rsa_n);
   rsa_assign(&gRSAPubKey.e, &rsa_d);

#if R__RSADEB
   if (gDebug > 2) {
      // Determine their lengths
      ErrorInfo("RpdGenRSAKeys: local: generated keys are:");
      ErrorInfo("RpdGenRSAKeys: local: n: '%s' length: %d", buf_n, l_n);
      ErrorInfo("RpdGenRSAKeys: local: e: '%s' length: %d", buf_e, l_e);
      ErrorInfo("RpdGenRSAKeys: local: d: '%s' length: %d", buf_d, l_d);
   }
#endif
   // Export form
   gRSAPubExport[0].len = l_n + l_d + 4;
   if (gRSAPubExport[0].keys)
      delete[] gRSAPubExport[0].keys;
   gRSAPubExport[0].keys = new char[gRSAPubExport[0].len];

   gRSAPubExport[0].keys[0] = '#';
   memcpy(gRSAPubExport[0].keys + 1, buf_n, l_n);
   gRSAPubExport[0].keys[l_n + 1] = '#';
   memcpy(gRSAPubExport[0].keys + l_n + 2, buf_d, l_d);
   gRSAPubExport[0].keys[l_n + l_d + 2] = '#';
   gRSAPubExport[0].keys[l_n + l_d + 3] = 0;
#if R__RSADEB
   if (gDebug > 2)
      ErrorInfo("RpdGenRSAKeys: local: export pub: '%s'",
                gRSAPubExport[0].keys);
#else
   if (gDebug > 2)
      ErrorInfo("RpdGenRSAKeys: local: export pub length: %d bytes",
                gRSAPubExport[0].len);
#endif

   gRSAInit = 1;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Generates local public/private RSA key pair
/// Send request for Client Public Key and Local public key
/// Receive encoded Client Key
/// Decode Client public key
/// NB: key is not saved to file here

int RpdRecvClientRSAKey()
{
   if (gRSAInit == 0) {
      // Generate Local RSA keys for the session
      if (RpdGenRSAKeys(1)) {
         ErrorInfo("RpdRecvClientRSAKey: unable to generate local keys");
         return 1;
      }
   }

   // Send server public key
   int key = gRSAKey - 1;
   NetSend(gRSAPubExport[key].keys, gRSAPubExport[key].len, kROOTD_RSAKEY);

   // Receive length of message with encode client public key
   EMessageTypes kind;
   char buflen[40];
   NetRecv(buflen, 20, kind);
   gPubKeyLen = atoi(buflen);
   if (gDebug > 3)
      ErrorInfo("RpdRecvClientRSAKey: got len '%s' %d ", buflen, gPubKeyLen);

   int nrec = 0;

   if (gRSAKey == 1) {

      // Receive and decode encoded public key
      nrec = NetRecvRaw(gPubKey, gPubKeyLen);

      rsa_decode(gPubKey, gPubKeyLen, gRSAPriKey.n, gRSAPriKey.e);
      if (gDebug > 2)
         ErrorInfo("RpdRecvClientRSAKey: Local: decoded string is %d bytes long ",
            strlen(gPubKey));
      gPubKeyLen = strlen(gPubKey);

   } else if (gRSAKey == 2) {
#ifdef R__SSL
      int ndec = 0;
      int lcmax = RSA_size(gRSASSLKey);
      char btmp[kMAXSECBUF];
      int nr = gPubKeyLen;
      int kd = 0;
      while (nr > 0) {
         // Receive and decode encoded public key
         nrec += NetRecvRaw(btmp, lcmax);
         if ((ndec = RSA_private_decrypt(lcmax,(unsigned char *)btmp,
                                    (unsigned char *)&gPubKey[kd],
                                    gRSASSLKey,
                                    RSA_PKCS1_PADDING)) < 0) {
            char errstr[120];
            ERR_error_string(ERR_get_error(), errstr);
            ErrorInfo("RpdRecvClientRSAKey: SSL: error: '%s' ",errstr);
         }
         nr -= lcmax;
         kd += ndec;
      }
      gPubKeyLen = kd;
#else
      if (gDebug > 0)
         ErrorInfo("RpdRecvClientRSAKey: not compiled with SSL support"
                   ": you should not have got here!");
      return 1;
#endif
   } else {
      if (gDebug > 0)
         ErrorInfo("RpdRecvClientRSAKey: unknown key type (%d)", gRSAKey);
   }


   // Import Key and Determine key type
   if (RpdGetRSAKeys(gPubKey, 0) != gRSAKey) {
      ErrorInfo("RpdRecvClientRSAKey:"
                " could not import a valid key (type %d)",gRSAKey);
      char *elogfile = new char[gRpdKeyRoot.length() + 11];
      SPrintf(elogfile, gRpdKeyRoot.length() + 11, "%.*serr.XXXXXX", (int)gRpdKeyRoot.length(), gRpdKeyRoot.c_str());
      mode_t oldumask = umask(0700);
      int ielog = mkstemp(elogfile);
      umask(oldumask);
      if (ielog != -1) {
         char line[kMAXPATHLEN] = {0};
         //
         SPrintf(line,kMAXPATHLEN,
                 " + RpdRecvClientRSAKey: error importing key\n + type: %d\n"
                 " + length: %d\n + key: %s\n + (%d bytes were received)",
                 gRSAKey, gPubKeyLen, gPubKey, nrec);
         while (write(ielog, line, strlen(line)) < 0 && GetErrno() == EINTR)
            ResetErrno();
         close (ielog);
      }
      delete [] elogfile;
      return 2;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Init random machine.

void RpdInitRand()
{
   const char *randdev = "/dev/urandom";

   int fd;
   unsigned int seed;
   if ((fd = open(randdev, O_RDONLY)) != -1) {
      if (gDebug > 2)
         ErrorInfo("RpdInitRand: taking seed from %s", randdev);
      if (read(fd, &seed, sizeof(seed))) {;}
      close(fd);
   } else {
      if (gDebug > 2)
         ErrorInfo("RpdInitRand: %s not available: using time()", randdev);
      seed = time(0);   //better use times() + win32 equivalent
   }
   srand(seed);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle user authentication.

int RpdAuthenticate()
{
   char buf[kMAXRECVBUF];
   EMessageTypes kind;

//#define R__DEBUG
#ifdef R__DEBUG
   int debug = 1;
   while (debug)
      ;
#endif

   // Reset gAuth (if we have been called this means that we need
   // to check at least that a valid authentication exists ...)
   int auth = 0;

   while (!auth) {

      // Receive next
      if (!gClientOld) {
         if (NetRecv(buf, kMAXRECVBUF, kind) < 0) {
            Error(gErr, -1, "RpdAuthenticate: error receiving message");
            return auth;
         }
      } else {
         strlcpy(buf,gBufOld, sizeof(buf));
         kind = gKindOld;
         gBufOld[0] = '\0';
         gClientOld = 0;
      }

      // If this is a rootd contacted via a TXNetFile we need to
      // receive again the buffer
      if (gService == kROOTD && kind == kROOTD_PROTOCOL) {
         if (NetRecv(buf, kMAXRECVBUF, kind) < 0) {
            Error(gErr, -1, "RpdAuthenticate: error receiving message");
            return auth;
         }
      }

      // Decode the method ...
      gAuthProtocol = RpdGetAuthMethod(kind);

      if (gDebug > 2) {
         if (kind != kROOTD_PASS) {
            ErrorInfo("RpdAuthenticate got: %d -- %s", kind, buf);
         } else {
            ErrorInfo("RpdAuthenticate got: %d ", kind);
         }
      }

      // Guess the client procotol if not received via Rootd/ProofdProtocol
      if (gClientProtocol == 0)
         gClientProtocol = RpdGuessClientProt(buf, kind);

      // If the client supports it check if we accept the method proposed;
      // if not send back the list of accepted methods, if any ...
      if (gAuthProtocol != -1 && gClientProtocol > 8) {

         // Check if accepted ...
         if (RpdCheckAuthAllow(gAuthProtocol, gOpenHost.c_str())) {
            if (gNumAllow>0) {
               if (gAuthListSent == 0) {
                  if (gDebug > 0)
                     ErrorInfo("Authenticate: %s method not"
                               " accepted from host: %s",
                                gAuthMeth[gAuthProtocol].c_str(),
                                gOpenHost.c_str());
                  NetSend(kErrNotAllowed, kROOTD_ERR);
                  RpdSendAuthList();
                  gAuthListSent = 1;
                  goto next;
               } else {
                  Error(gErr,kErrNotAllowed,"Authenticate: method not"
                       " in the list sent to the client");
                  return auth;
               }
            } else {
               Error(gErr,kErrConnectionRefused,"Authenticate:"
                       " connection refused from host %s", gOpenHost.c_str());
               return auth;
            }
         }

         // Then check if a previous authentication exists and is valid
         if ((auth = RpdReUseAuth(buf, kind)))
            goto next;
      }

      // Reset global variable
      auth = 0;

      switch (kind) {
         case kROOTD_USER:
            auth = RpdUser(buf);
            break;
         case kROOTD_PASS:
            auth = RpdPass(buf);
            break;
         case kROOTD_CLEANUP:
            RpdAuthCleanup(buf,1);
            ErrorInfo("RpdAuthenticate: authentication stuff cleaned - exit");
            // Fallthrough next case now to free the keys
         case kROOTD_BYE:
            RpdFreeKeys();
            return auth;
         default:
            Error(gErr,-1,"RpdAuthenticate: received bad opcode %d", kind);
            return auth;
      }

      if (gClientProtocol > 8) {

         // If failure prepare or continue negotiation
         int doneg = (gAuthProtocol != -1 || kind == kROOTD_PASS);
         if (gDebug > 2 && doneg)
            ErrorInfo("RpdAuthenticate: kind:%d meth:%d auth:%d gNumLeft:%d",
                      kind, gAuthProtocol, auth, gNumLeft);

         // If authentication failure, check if other methods could be tried ...
         if (auth == 0 && doneg) {
            if (gNumLeft > 0) {
               if (gAuthListSent == 0) {
                  RpdSendAuthList();
                  gAuthListSent = 1;
               } else
                  NetSend(-1, kROOTD_NEGOTIA);
            } else {
               NetSend(0, kROOTD_NEGOTIA);
               Error(gErr, -1, "RpdAuthenticate: authentication failed");
               return auth;
            }
         }
      }
next:
      continue;
   }

   return auth;
}
////////////////////////////////////////////////////////////////////////////////
/// Free space allocated for encryption keys

void RpdFreeKeys()
{
   if (gRSAPubExport[0].keys) delete[] gRSAPubExport[0].keys;
   if (gRSAPubExport[1].keys) delete[] gRSAPubExport[1].keys;
#ifdef R__SSL
   RSA_free(gRSASSLKey);
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Receives client protocol and returns daemon protocol.
/// Returns:  0 if ok
///          -1 if any error occured
///          -2 if special action (e.g. cleanup): no need to continue

int RpdProtocol(int ServType)
{
   int rc = 0;

//#define R__DEBUG
#ifdef R__DEBUG
   int debug = 1;
   while (debug)
      ;
#endif

   if (gDebug > 2)
      ErrorInfo("RpdProtocol: Enter: server type = %d", ServType);

   int readbuf = 1;
   EMessageTypes kind;
   char proto[kMAXRECVBUF];

   // For backward compatibility, for rootd we need to understand
   // whether we are talking to a OLD client: protocol information is
   // available only later on ...
   int lbuf[2];
   if (NetRecvRaw(lbuf, sizeof(lbuf)) < 0) {
      NetSend(kErrFatal, kROOTD_ERR);
      ErrorInfo("RpdProtocol: error receiving message");
      return -1;
   }

   // if kind is {kROOTD_PROTOCOL, kROOTD_CLEANUP}
   // receive the rest
   kind = (EMessageTypes) ntohl(lbuf[1]);
   int len = ntohl(lbuf[0]);
   if (gDebug > 1)
      ErrorInfo("RpdProtocol: kind: %d %d",kind,len);
   if (kind == kROOTD_PROTOCOL || kind == kROOTD_CLEANUP) {
      // Receive the rest
      char *buf = 0;
      len -= sizeof(int);
      if (gDebug > 1)
         ErrorInfo("RpdProtocol: len: %d",len);
      if (len) {
         buf = new char[len];
         if (NetRecvRaw(buf, len) < 0) {
            NetSend(kErrFatal, kROOTD_ERR);
            ErrorInfo("RpdProtocol: error receiving message");
            delete[] buf;
            return -1;
         }
         strlcpy(proto, buf, sizeof(proto));
      } else {
         // Empty buffer
         proto[0] = '\0';
      }
      if (gDebug > 1)
         ErrorInfo("RpdProtocol: proto buff: %s", buf ? buf : "---");
      // Copy buffer for later use
      readbuf = 0;
      if (buf) delete[] buf;
   } else if (ServType == kROOTD && kind == 0 && len == 0) {
      // TNetFile via TXNetFile: receive client protocol
      // read first next 12 bytes and discard them
      int llen = 12;
      char *buf = new char[llen];
      if (NetRecvRaw(buf, llen) < 0) {
         NetSend(kErrFatal, kROOTD_ERR);
         ErrorInfo("RpdProtocol: error receiving message");
         if (buf) delete[] buf;
         return -1;
      }
      if (buf) delete[] buf;
      // Send back the 'type'
      int type = htonl(8);
      if (NetSendRaw(&type,sizeof(type)) < 0) {
         NetSend(kErrFatal, kROOTD_ERR);
         ErrorInfo("RpdProtocol: error sending type to TXNetFile");
         return -1;
      }
      // Now read the client protocol
      llen = 4;
      buf = new char[llen];
      if (NetRecvRaw(buf,llen) < 0) {
         NetSend(kErrFatal, kROOTD_ERR);
         ErrorInfo("RpdProtocol: error receiving message");
         delete[] buf;
         return -1;
      }
      strlcpy(proto,buf, sizeof(proto));
      kind = kROOTD_PROTOCOL;
      readbuf = 0;
      delete[] buf;
   } else {
      // Need to open parallel sockets first
      int size = ntohl(lbuf[1]);
      // Read port
      int port;
      if (NetRecvRaw(&port, sizeof(int)) < 0) {
         NetSend(kErrFatal, kROOTD_ERR);
         ErrorInfo("RpdProtocol: error receiving message");
         return -1;
      }
      port = ntohl(port);
      if (gDebug > 0)
         ErrorInfo("RpdProtocol: port = %d, size = %d", port, size);
      if (size > 1)
         NetParOpen(port, size);
   }

   int done = 0;
   gClientOld = 0;
   while (!done) {

      // Receive next
      if (readbuf) {
         if (NetRecv(proto, kMAXRECVBUF, kind) < 0) {
            ErrorInfo("RpdProtocol: error receiving message");
            return -1;
         }
      }
      readbuf = 1;

      switch(kind) {

         case kROOTD_CLEANUP:
            RpdAuthCleanup(proto,1);
            ErrorInfo("RpdProtocol: authentication stuff cleaned");
            done = 1;
            rc = -2;
            break;
         case kROOTD_BYE:
            RpdFreeKeys();
            NetClose();
            done = 1;
            rc = -2;
            break;
         case kROOTD_PROTOCOL:

            if (strlen(proto) > 0) {
               gClientProtocol = atoi(proto);
            } else {
               if (ServType == kROOTD) {
                  // This is an old (TNetFile,TFTP) client:
                  // send our protocol first ...
                  if (NetSend(gServerProtocol, kROOTD_PROTOCOL) < 0) {
                     ErrorInfo("RpdProtocol: error sending kROOTD_PROTOCOL");
                     rc = -1;
                  }
                  // ... and receive protocol via kROOTD_PROTOCOL2
                  if (NetRecv(proto, kMAXRECVBUF, kind) < 0) {
                     ErrorInfo("RpdProtocol: error receiving message");
                     rc = -1;
                  }
                  if (kind != kROOTD_PROTOCOL2) {
                     strlcpy(gBufOld, proto, sizeof(gBufOld));
                     gKindOld = kind;
                     gClientOld = 1;
                     gClientProtocol = 0;
                  } else
                     gClientProtocol = atoi(proto);
               } else
                  gClientProtocol = 0;
            }
            if (!gClientOld) {
               // send our protocol
               // if we do not require authentication say it here
               Int_t protoanswer = gServerProtocol;
               if (!gRequireAuth && gClientProtocol > 10)
                  protoanswer += 1000;
               // Notify
               if (gDebug > 0) {
                  ErrorInfo("RpdProtocol: gClientProtocol = %d",
                            gClientProtocol);
                  ErrorInfo("RpdProtocol: Sending gServerProtocol = %d",
                            protoanswer);
               }
               if (NetSend(protoanswer, kROOTD_PROTOCOL) < 0) {
                  ErrorInfo("RpdProtocol: error sending kROOTD_PROTOCOL");
                  rc = -1;
               }
            }
            done = 1;
            break;
         default:
            ErrorInfo("RpdProtocol: received bad option (%d)",kind);
            rc = -1;
            done = 1;
            break;
      } // Switch

   } // done

   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Authentication was successful, set user environment.

int RpdLogin(int ServType, int auth)
{
//   if (gDebug > 2)
      ErrorInfo("RpdLogin: enter: Server: %d, gUser: %s, auth: %d",
                ServType, gUser, auth);

   // Login only if requested
   if (gDoLogin == 0)
      return -2;

   struct passwd *pw = getpwnam(gUser);

   if (!pw) {
      ErrorInfo("RpdLogin: user %s does not exist locally\n", gUser);
      return -1;
   }

   if (getuid() == 0) {
      // Anonymous users are confined to their corner
      if (gAnon) {
         // We need to do it before chroot, otherwise it does not work
         if (chdir(pw->pw_dir) == -1) {
            ErrorInfo("RpdLogin: can't change directory to %s (errno: %d)",
                      pw->pw_dir, errno);
            return -1;
         }
         if (chroot(pw->pw_dir) == -1) {
            ErrorInfo("RpdLogin: can't chroot to %s", pw->pw_dir);
            return -1;
         }
      }

      // set access control list from /etc/initgroup
      initgroups(gUser, pw->pw_gid);

      // set uid and gid
      if (setresgid(pw->pw_gid, pw->pw_gid, 0) == -1) {
         ErrorInfo("RpdLogin: can't setgid for user %s", gUser);
         return -1;
      }
      if (setresuid(pw->pw_uid, pw->pw_uid, 0) == -1) {
         ErrorInfo("RpdLogin: can't setuid for user %s", gUser);
         return -1;
      }
   }

   if (ServType == 2) {
      // set HOME env
      char *home = new char[8+strlen(pw->pw_dir)];
      SPrintf(home, 8+strlen(pw->pw_dir), "HOME=%s", pw->pw_dir);
      putenv(home);
   }

   // Change user's HOME, if required (for anon it is already done)
   if (gDoLogin == 2 && !gAnon) {
      if (chdir(pw->pw_dir) == -1) {
         ErrorInfo("RpdLogin: can't change directory to %s (errno: %d)",
                   pw->pw_dir, errno);
         return -1;
      }
   }

   umask(022);

   // Notify authentication to client ...
   NetSend(auth, kROOTD_AUTH);
   // Send also new offset if it changed ...
   if (auth == 2) NetSend(gOffSet, kROOTD_AUTH);

   if (gDebug > 0)
      ErrorInfo("RpdLogin: user %s logged in", gUser);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform the action needed to commence the new session:
/// Version called by TServerSocket.
///   - set debug flag
///   - check authentication table
///   - Inquire protocol
///   - authenticate the client
/// Returns logged-in user, the remote client procotol cproto,
/// the authentication protocol (ROOT internal) number is returned
/// in meth, type indicates the kind of authentication:
///       0 = new authentication
///       1 = existing authentication
///       2 = existing authentication with updated offset
/// and the crypted-token in ctoken (used later for cleaning).
/// Called just after opening the connection

int RpdInitSession(int servtype, std::string &user,
                   int &cproto, int &meth, int &type, std::string &ctoken)
{
   std::string pwd;
   int auth = RpdInitSession(servtype,user,cproto,meth,pwd);
   if (auth == 1)
      if (gExistingAuth)
         type = 1;
      else
         type = 0;
   else if (auth == 2)
      type = 2;
   ctoken = gCryptToken;

   return auth;
}
////////////////////////////////////////////////////////////////////////////////
/// Perform the action needed to commence the new session:
///   - set debug flag
///   - check authentication table
///   - Inquire protocol
///   - authenticate the client
///   - login the client
/// Returns 1 for a PROOF master server, 0 otherwise
/// Returns logged-in user, the remote client procotol cproto, the
/// client kind of user anon and, if anonymous user, the client passwd.
/// If TServerSocket (servtype==kSOCKD), the protocol number is returned
/// in anon.
/// Called just after opening the connection

int RpdInitSession(int servtype, std::string &user,
                   int &cproto, int &anon, std::string &passwd)
{
   if (gDebug > 2)
      ErrorInfo("RpdInitSession: %s", gServName[servtype].c_str());

   int retval = 0;

   // CleanUp authentication table, if needed or required ...
   RpdInitAuth();

   // Get Host name
   NetGetRemoteHost(gOpenHost);

   if (servtype == kPROOFD) {

      // find out if we are supposed to be a master or a slave server
      char  msg[80];
      if (NetRecv(msg, sizeof(msg)) < 0) {
         ErrorInfo("RpdInitSession: Cannot receive master/slave status");
         return -1;
      }

      retval = !strcmp(msg, "master") ? 1 : 0;

      if (gDebug > 0)
         ErrorInfo("RpdInitSession: PROOF master/slave = %s", msg);
   }

   // Get protocol first
   // Failure typically indicate special actions like cleanup
   // which do not need additional work
   // The calling program will then decide what to do
   int rcp = RpdProtocol(servtype);
   if (rcp != 0) {
      if (rcp == -1)
         ErrorInfo("RpdInitSession: error getting remote protocol");
      else if (rcp != -2)
         ErrorInfo("RpdInitSession: unknown error from RpdProtocol");
      return rcp;
   }

   // Check if authentication is required
   // Old clients do not support no authentication mode
   bool runAuth = (gClientProtocol < 11 || gRequireAuth) ? 1 : 0;

   // user authentication (does not return in case of failure)
   int auth = 0;
   if (runAuth) {
      auth = RpdAuthenticate();
      if (auth == 0) {
         ErrorInfo("RpdInitSession: unsuccessful authentication attempt");
         return -1;
      }
   } else {
      auth = RpdNoAuth(servtype);
   }

   // Login the user (if in rootd/proofd environment)
   if (gDoLogin > 0) {
      if (RpdLogin(servtype,auth) != 0) {
         ErrorInfo("RpdInitSession: unsuccessful login attempt");
         // Notify failure to client ...
         NetSend(0, kROOTD_AUTH);
         return -1;
      }
   } else {
      // Notify authentication to client ...
      NetSend(auth, kROOTD_AUTH);
      // Send also new offset if it changed ...
      if (auth == 2)
         NetSend(gOffSet, kROOTD_AUTH);
      if (gDebug > 0)
         ErrorInfo("RpdInitSession: User '%s' authenticated", gUser);
      retval = auth;
   }

   // Output vars
   user = std::string(gUser);
   cproto = gClientProtocol;
   if (servtype == kSOCKD)
      anon = gSec;
   else
      anon = gAnon;
   if (gAnon)
      passwd = std::string(gPasswd);

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform the action needed to commence the new session:
///   - set debug flag
///   - check authentication table
///   - Inquire protocol
///   - authenticate the client
///   - login the client
/// Returns 1 for a PROOF master server, 0 otherwise
/// Returns logged-in user and remote process id in rid
/// Called just after opening the connection

int RpdInitSession(int servtype, std::string &user, int &rid)
{
   int dum1 = 0, dum2 = 0;
   std::string dum3;
   rid = gRemPid;
   return RpdInitSession(servtype,user,dum1,dum2,dum3);

}


////////////////////////////////////////////////////////////////////////////////
/// Perform entrance formalities in case of no authentication
/// mode, i.e. get target user and check if authorized
/// Don't return if something goes wrong

int RpdNoAuth(int servtype)
{
   if (gDebug > 1)
      ErrorInfo("RpdNoAuth: no authentication required");

   // Special value for this case
   int auth = 0;

   // Receive target username
   if (servtype == kROOTD || servtype == kPROOFD) {

      char buf[kMAXPATHLEN];
      EMessageTypes kind;
      if (NetRecv(buf, kMAXPATHLEN, kind) < 0) {
         NetSend(kErrBadMess, kROOTD_ERR);
         ErrorInfo("RpdNoAuth: error receiving target user");
         goto quit;
      }

      if (kind == kROOTD_BYE)
         goto quit;

      if (kind != kROOTD_USER) {
         NetSend(kErrBadOp, kROOTD_ERR);
         ErrorInfo("RpdNoAuth: protocol error:"
          " received msg type: %d, expecting: %d", kind, kROOTD_USER);
         goto quit;
      }

      // Decode buffer
      char ruser[kMAXUSERLEN], user[kMAXUSERLEN];
      int nw = sscanf(buf,"%64s %64s",ruser,user);
      if (nw <= 0 || !strcmp(ruser,"-1")) {
         NetSend(kErrBadMess, kROOTD_ERR);
         ErrorInfo("RpdNoAuth: received uncorrect information: %s", buf);
         goto quit;
      }
      // If target user not send, assume user == ruser
      if (nw == 1)
         snprintf(user,kMAXUSERLEN,"%s",ruser);

      struct passwd *pw = 0;
      if ((pw = getpwnam(user)) == 0) {
         NetSend(kErrNoUser, kROOTD_ERR);
         ErrorInfo("RpdNoAuth: user %s unknown", user);
         goto quit;
      }

      // If server is not started as root, require user to be the same
      // as the one who started rootd
      uid_t uid = getuid();
      if (uid && uid != pw->pw_uid) {
         NetSend(kErrBadUser, kROOTD_ERR);
         ErrorInfo("RpdNoAuth: user not same as effective user of rootd");
         goto quit;
      }

      if (gDebug > 2)
         ErrorInfo("RpdNoAuth: remote user: %s, target user: %s",ruser,user);

      SPrintf(gUser, 63, "%s", user);
   }

   auth = 4;

   quit:
   return auth;
}

////////////////////////////////////////////////////////////////////////////////
/// Change current user id to uid (and gid).

int RpdSetUid(int uid)
{
   if (gDebug > 2)
      ErrorInfo("RpdSetUid: enter ...uid: %d", uid);

   struct passwd *pw = getpwuid(uid);

   if (!pw) {
      ErrorInfo("RpdSetUid: uid %d does not exist locally", uid);
      return -1;
   } else if (chdir(pw->pw_dir) == -1) {
      ErrorInfo("RpdSetUid: can't change directory to %s", pw->pw_dir);
      return -1;
   }

   if (getuid() == 0) {

      // set access control list from /etc/initgroup
      initgroups(pw->pw_name, pw->pw_gid);

      // set uid and gid
      if (setresgid(pw->pw_gid, pw->pw_gid, 0) == -1) {
         ErrorInfo("RpdSetUid: can't setgid for uid %d", uid);
         return -1;
      }
      if (setresuid(pw->pw_uid, pw->pw_uid, 0) == -1) {
         ErrorInfo("RpdSetUid: can't setuid for uid %d", uid);
         return -1;
      }
   }

   if (gDebug > 0)
      ErrorInfo("RpdSetUid: uid set (%d,%s)", uid, pw->pw_name);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Change defaults job control options.

void RpdInit(EService serv, int pid, int sproto, unsigned int options,
             int rumsk, int, const char *tmpd, const char * /* asrpp */, int login)
{
   gService        = serv;
   gParentId       = pid;
   gServerProtocol = sproto;
   gReUseAllow     = rumsk;
   gDoLogin        = login;

   // Parse options
   gCheckHostsEquiv= (bool)((options & kDMN_HOSTEQ) != 0);
   gRequireAuth    = (bool)((options & kDMN_RQAUTH) != 0);
   gSysLog         = (bool)((options & kDMN_SYSLOG) != 0);

   if (tmpd && strlen(tmpd)) {
      gTmpDir      = tmpd;
      gRpdAuthTab  = gTmpDir + gAuthTab;
      gRpdKeyRoot  = gTmpDir + gKeyRoot;
   }
   // Auth Tab and public key files are exclusive to this family
   gRpdAuthTab.append(".");
   gRpdAuthTab.append(ItoA(getuid()));
   gRpdKeyRoot.append(ItoA(getuid()));
   gRpdKeyRoot.append("_");

   if (gDebug > 0) {
      ErrorInfo("RpdInit: gService= %s, gSysLog= %d",
                 gServName[gService].c_str(), gSysLog);
      ErrorInfo("RpdInit: gParentId= %d", gParentId);
      ErrorInfo("RpdInit: gRequireAuth= %d, gCheckHostEquiv= %d",
                 gRequireAuth, gCheckHostsEquiv);
      ErrorInfo("RpdInit: gReUseAllow= 0x%x", gReUseAllow);
      ErrorInfo("RpdInit: gServerProtocol= %d", gServerProtocol);
      ErrorInfo("RpdInit: gDoLogin= %d", gDoLogin);
      if (tmpd)
         ErrorInfo("RpdInit: gTmpDir= %s", gTmpDir.c_str());
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Acts like snprintf with some printout in case of error if required
/// Returns number of  characters printed (excluding the trailing `\0').
/// Returns 0 is buf or size are not defined or inconsistent.
/// Returns -1 if the buffer is truncated.

int SPrintf(char *buf, size_t size, const char *va_(fmt), ...)
{
   // Check buf
   if (!buf) {
      if (gDebug > 0)
         ErrorInfo("SPrintf: buffer not allocated: do nothing");
      return 0;
   }

   // Check size
   if (size < 1) {
      if (gDebug > 0)
         ErrorInfo("SPrintf: cannot determine buffer size (%d): do nothing",size);
      return 0;
   }

   // Now fill buf
   va_list ap;
   va_start(ap,va_(fmt));
   int np = vsnprintf(buf, size, fmt, ap);
   va_end(ap);

   if (np == -1 && gDebug > 0)
      ErrorInfo("SPrintf: buffer truncated (%s)",buf);

   return np;
}


////////////////////////////////////////////////////////////////////////////////
/// Return pointer to a static string containing the string
/// version of integer 'i', up to a max of kMAXCHR (=30)
/// characters; returns "-1" if more chars are needed.

char *ItoA(int i)
{
   const int kMAXCHR = 30;
   static char str[kMAXCHR];

   // This is the number of characters we need
   int nchr = (int)log10(double(i)) + 1;
   if (nchr > kMAXCHR)
      strlcpy(str,"-1", sizeof(str));
   else
      snprintf(str,30,"%d",i);

   return str;
}

////////////////////////////////////////////////////////////////////////////////
/// Set global pointers to error handler functions

void RpdSetErrorHandler(ErrorHandler_t err, ErrorHandler_t sys, ErrorHandler_t fatal)
{
   gErr      = err;
   gErrSys   = sys;
   gErrFatal = fatal;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve specific ROOT password from $HOME/fpw, if any.
/// To avoid problems with NFS-root-squashing, if 'root' changes temporarly the
/// uid/gid to those of the target user (usr).
/// If OK, returns pass length and fill 'pass' with the password, null-terminated.
/// ('pass' is allocated externally to contain max lpwmax bytes).
/// If the file does not exists, return 0 and an empty pass.
/// If any problems with the file occurs, return a negative
/// code, -2 indicating wrong file permissions.
/// If any problem with changing ugid's occurs, prints a warning trying anyhow
/// to read the password hash.

int RpdRetrieveSpecialPass(const char *usr, const char *fpw, char *pass, int lpwmax)
{
   int rc = -1;
   int len = 0, n = 0, fid = -1;

   // Check inputs
   if (!usr || !pass) {
      if (gDebug > 0)
         ErrorInfo("RpdRetrieveSpecialPass: invalid arguments:"
                   " us:%p, sp:%p", usr, pass);
      return rc;
   }

   struct passwd *pw = getpwnam(usr);
   if (!pw) {
      if (gDebug > 0)
         ErrorInfo("RpdRetrieveSpecialPass: user '%s' does not exist", usr);
      return rc;
   }

   // target and actual uid
   int uid = pw->pw_uid;
   int ouid = getuid();

   // Temporary change to target user ID to avoid NFS squashing problems
   if (ouid == 0) {

      // set access control list from /etc/initgroup
      if (initgroups(pw->pw_name, pw->pw_gid) == -1)
         ErrorInfo("RpdRetrieveSpecialPass: can't initgroups for uid %d"
                   " (errno: %d)", uid, GetErrno());
      // set uid and gid
      if (setresgid(pw->pw_gid, pw->pw_gid, 0) == -1)
         ErrorInfo("RpdRetrieveSpecialPass: can't setgid for gid %d"
                   " (errno: %d)", pw->pw_gid, GetErrno());
      if (setresuid(pw->pw_uid, pw->pw_uid, 0) == -1)
         ErrorInfo("RpdRetrieveSpecialPass: can't setuid for uid %d"
                   " (errno: %d)", uid, GetErrno());
   }

   // The file now
   char rootdpass[kMAXPATHLEN];
   SPrintf(rootdpass, kMAXPATHLEN, "%s/%s", pw->pw_dir, fpw);

   if (gDebug > 0)
      ErrorInfo
         ("RpdRetrieveSpecialPass: checking file %s for user %s",rootdpass,
           pw->pw_name);


   if ((fid = open(rootdpass, O_RDONLY)) == -1) {
      ErrorInfo("RpdRetrieveSpecialPass: cannot open password file"
                " %s (errno: %d)", rootdpass, GetErrno());
      rc = -1;
      goto back;
   }
   // Check first the permissions: should be 0600
   struct stat st;
   if (fstat(fid, &st) == -1) {
      ErrorInfo("RpdRetrieveSpecialPass: cannot stat descriptor %d"
                  " %s (errno: %d)", fid, GetErrno());
      close(fid);
      rc = -1;
      goto back;
   }
   if (!S_ISREG(st.st_mode) || S_ISDIR(st.st_mode) ||
       (st.st_mode & (S_IWGRP | S_IWOTH | S_IRGRP | S_IROTH)) != 0) {
      ErrorInfo("RpdRetrieveSpecialPass: pass file %s: wrong permissions"
                " 0%o (should be 0600)", rootdpass, (st.st_mode & 0777));
      ErrorInfo("RpdRetrieveSpecialPass: %d %d",
                S_ISREG(st.st_mode),S_ISDIR(st.st_mode));
      close(fid);
      rc = -2;
      goto back;
   }

   if ((n = read(fid, pass, lpwmax - 1)) <= 0) {
      close(fid);
      ErrorInfo("RpdRetrieveSpecialPass: cannot read password file"
                " %s (errno: %d)", rootdpass, GetErrno());
      rc = -1;
      goto back;
   }
   close(fid);

   // Get rid of special trailing chars
   len = n;
   while (len-- && (pass[len] == '\n' || pass[len] == 32))
      pass[len] = 0;

   // Null-terminate
   pass[++len] = 0;
   rc = len;

   back:
   // Change back uid's
   if (ouid == 0) {
      // set uid and gid
      if (setresgid(0, 0, 0) == -1)
         ErrorInfo("RpdRetrieveSpecialPass: can't re-setgid for gid 0"
                   " (errno: %d)", GetErrno());
      if (setresuid(0, 0, 0) == -1)
         ErrorInfo("RpdRetrieveSpecialPass: can't re-setuid for uid 0"
                   " (errno: %d)", GetErrno());
   }

   // We are done
   return rc;
}

} // namespace ROOT
