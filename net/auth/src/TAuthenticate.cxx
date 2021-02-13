// @(#)root/auth:$Id: f2cfa663e232707e1201467b5805ff1d13575326 $
// Author: Fons Rademakers   26/11/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TAuthenticate                                                        //
//                                                                      //
// An authentication module for ROOT based network services, like rootd //
// and proofd.                                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RConfigure.h"

#include "TAuthenticate.h"
#include "TApplication.h"
#include "THostAuth.h"
#include "TRootSecContext.h"
#include "TPluginManager.h"
#include "TNetFile.h"
#include "TPSocket.h"
#include "TMessage.h"
#include "TSystem.h"
#include "TError.h"
#include "Getline.h"
#include "TROOT.h"
#include "TEnv.h"
#include "TList.h"
#include "NetErrors.h"
#include "TRegexp.h"
#include "TVirtualMutex.h"
#include "TTimer.h"
#include "TBase64.h"
#include "strlcpy.h"
#include "snprintf.h"

#include "rsafun.h"

#ifndef R__LYNXOS
#include <sys/stat.h>
#endif
#include <errno.h>
#include <sys/types.h>
#include <time.h>
#if !defined(R__WIN32) && !defined(R__MACOSX) && !defined(R__FBSD) && \
    !defined(R__OBSD)
#include <crypt.h>
#endif
#ifdef WIN32
#  include <io.h>
#endif /* WIN32 */
#if defined(R__LINUX) || defined(R__FBSD) || defined(R__OBSD)
#  include <unistd.h>
#endif
#include <stdlib.h>
#ifndef WIN32
#  include <sys/time.h>
#endif /* WIN32 */

#if defined(R__MACOSX)
extern "C" char *crypt(const char *, const char *);
#endif

#ifdef R__SSL
// SSL specific headers
#   include <openssl/bio.h>
#   include <openssl/err.h>
#   include <openssl/pem.h>
#   include <openssl/rand.h>
#   include <openssl/rsa.h>
#   include <openssl/ssl.h>
#   include <openssl/blowfish.h>
#endif

struct R__rsa_KEY: rsa_KEY { R__rsa_KEY(): rsa_KEY() {} };
struct R__rsa_KEY_export: rsa_KEY_export {};
struct R__rsa_NUMBER: rsa_NUMBER {};

#ifdef R__SSL
   static BF_KEY fgBFKey; // Blowfish symmetric key
#endif

// Statics initialization
TList          *TAuthenticate::fgAuthInfo = 0;
TString         TAuthenticate::fgAuthMeth[] = { "UsrPwd", "Unsupported", "Unsupported",
                                                "Unsupported", "Unsupported", "Unsupported" };
Bool_t          TAuthenticate::fgAuthReUse;
TString         TAuthenticate::fgDefaultUser;
TDatime         TAuthenticate::fgExpDate;
TDatime         TAuthenticate::fgLastAuthrc;    // Time of last reading of fgRootAuthrc
TString         TAuthenticate::fgPasswd;
TPluginHandler *TAuthenticate::fgPasswdDialog = (TPluginHandler *)(-1);
Bool_t          TAuthenticate::fgPromptUser;
TList          *TAuthenticate::fgProofAuthInfo = 0;
Bool_t          TAuthenticate::fgPwHash;
Bool_t          TAuthenticate::fgReadHomeAuthrc = kTRUE; // on/off search for $HOME/.rootauthrc
TString         TAuthenticate::fgRootAuthrc;    // Path to last rootauthrc-like file read
Int_t           TAuthenticate::fgRSAKey  = -1;  // Default RSA key type to be used
Int_t           TAuthenticate::fgRSAInit = 0;
R__rsa_KEY         TAuthenticate::fgRSAPriKey;
R__rsa_KEY_export R__fgRSAPubExport[2] = {{}, {}};
R__rsa_KEY_export* TAuthenticate::fgRSAPubExport = R__fgRSAPubExport;
R__rsa_KEY         TAuthenticate::fgRSAPubKey;
SecureAuth_t    TAuthenticate::fgSecAuthHook;
TString         TAuthenticate::fgUser;
Bool_t          TAuthenticate::fgUsrPwdCrypt;
Int_t           TAuthenticate::fgLastError = -1;
Int_t           TAuthenticate::fgAuthTO = -2;       // Timeout value

// ID of the main thread as unique identifier
Int_t           TAuthenticate::fgProcessID = -1;

TVirtualMutex *gAuthenticateMutex = 0;

// Standard version of Sec Context match checking
Int_t StdCheckSecCtx(const char *, TRootSecContext *);


ClassImp(TAuthenticate);

////////////////////////////////////////////////////////////////////////////////
/// rand() implementation using /udev/random or /dev/random, if available

static int auth_rand()
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
   Printf("+++ERROR+++ : auth_rand: neither /dev/urandom nor /dev/random are available or readable!");
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
/// Create authentication object.

TAuthenticate::TAuthenticate(TSocket *sock, const char *remote,
                             const char *proto, const char *user)
{
   if (gDebug > 2 && gAuthenticateMutex)
      Info("Authenticate", "locking mutex (pid:  %d)",gSystem->GetPid());
   R__LOCKGUARD2(gAuthenticateMutex);

   // In PROOF decode the buffer sent by the client, if any
   if (gROOT->IsProofServ())
      ProofAuthSetup();

   // Use the ID of the starting thread as unique identifier
   if (fgProcessID < 0)
      fgProcessID = gSystem->GetPid();

   if (fgAuthTO == -2)
      fgAuthTO = gEnv->GetValue("Auth.Timeout",-1);

   fSocket   = sock;
   fRemote   = remote;
   fHostAuth = 0;
   fVersion  = 5;                // The latest, by default
   fSecContext = 0;

   if (gDebug > 2)
      Info("TAuthenticate", "Enter: local host: %s, user is: %s (proto: %s)",
           gSystem->HostName(), user, proto);

   // Set protocol string.
   // Check if version should be different ...
   char *pdd;
   Int_t servtype = TSocket::kSOCKD;
   if (proto && strlen(proto) > 0) {
      char *sproto = StrDup(proto);
      if ((pdd = strstr(sproto, ":")) != 0) {
         int rproto = atoi(pdd + 1);
         *pdd = '\0';
         if (strstr(sproto, "root") != 0) {
            if (rproto < 12 ) {
               fVersion = 4;
               if (rproto < 11 ) {
                  fVersion = 3;
                  if (rproto < 9 ) {
                     fVersion = 2;
                     if (rproto < 8) {
                        fVersion = 1;
                        if (rproto < 6)
                           fVersion = 0;
                     }
                  }
               }
            }
            servtype = TSocket::kROOTD;
         }
         if (strstr(sproto, "proof") != 0) {
            if (rproto < 11) {
               fVersion = 4;
               if (rproto < 10) {
                  fVersion = 3;
                  if (rproto < 8) {
                     fVersion = 2;
                     if (rproto < 7)
                        fVersion = 1;
                  }
               }
            }
            servtype = TSocket::kPROOFD;
         }
         if (gDebug > 3)
            Info("TAuthenticate",
                 "service: %s (remote protocol: %d): fVersion: %d", sproto,
                 rproto, fVersion);
      }
      fProtocol = sproto;
      delete [] sproto;
   }

   // Check or get user name
   fUser = "";
   TString checkUser;
   if (user && strlen(user) > 0) {
      fUser = user;
      checkUser = user;
   } else {
      UserGroup_t *u = gSystem->GetUserInfo();
      if (u)
         checkUser = u->fUser;
      delete u;
   }
   fPasswd = "";
   fPwHash = kFALSE;

   // Type of RSA key
   if (fgRSAKey < 0) {
      fgRSAKey  = 0;                // Default key
#ifdef R__SSL
      // Another choice possible: check user preferences
      if (gEnv->GetValue("RSA.KeyType",0) == 1)
         fgRSAKey = 1;
#endif
   }
   // This is the key actually used: we propose the default
   // to the server, and behave according to its reply
   fRSAKey = fgRSAKey;
   if (gDebug > 3)
      Info("TAuthenticate","RSA key: default type %d", fgRSAKey);

   // RSA key generation (one per session)
   if (!fgRSAInit) {
      GenRSAKeys();
      fgRSAInit = 1;
   }

   // Check and save the host FQDN ...
   TString fqdn;
   TInetAddress addr = gSystem->GetHostByName(fRemote);
   if (addr.IsValid())
      fqdn = addr.GetHostName();
   TString fqdnsrv;
   fqdnsrv.Form("%s:%d",fqdn.Data(),servtype);

   // Read directives from files; re-read if files have changed
   TAuthenticate::ReadRootAuthrc();

   if (gDebug > 3) {
      Info("TAuthenticate",
           "number of HostAuth Instantiations in memory: %d",
           GetAuthInfo()->GetSize());
      TAuthenticate::Show("H");
      TAuthenticate::Show("P");
   }

   // Check the list of auth info for already loaded info about this host
   fHostAuth = GetHostAuth(fqdnsrv, checkUser);

   //
   // If generic THostAuth (i.e. with wild card or user == any)
   // make a personalized memory copy of this THostAuth
   if (strchr(fHostAuth->GetHost(),'*') || strchr(fHostAuth->GetHost(),'*') ||
       fHostAuth->GetServer() == -1 ) {
      fHostAuth = new THostAuth(*fHostAuth);
      fHostAuth->SetHost(fqdn);
      fHostAuth->SetUser(checkUser);
      fHostAuth->SetServer(servtype);
   }

   // If a specific method has been requested via the protocol
   // set it as first
   Int_t sec = -1;
   TString tmp = fProtocol;
   tmp.ReplaceAll("root",4,"",0);
   tmp.ReplaceAll("proof",5,"",0);
   tmp.ReplaceAll("sock",4,"",0);
   if (!strncmp(tmp.Data(),"up",2))
      sec = 0;
   else if (!strncmp(tmp.Data(),"s",1))
      sec = 1;
   else if (!strncmp(tmp.Data(),"k",1))
      sec = 2;
   else if (!strncmp(tmp.Data(),"g",1))
      sec = 3;
   else if (!strncmp(tmp.Data(),"h",1))
      sec = 4;
   else if (!strncmp(tmp.Data(),"ug",2))
      sec = 5;
   if (sec > -1 && sec < kMAXSEC) {
      if (fHostAuth->HasMethod(sec)) {
         fHostAuth->SetFirst(sec);
      } else {
         char *dtmp = GetDefaultDetails(sec, 1, checkUser);
         TString det(dtmp);
         fHostAuth->AddFirst(sec, det);
         if (dtmp)
            delete [] dtmp;
      }
   }

   // This is what we have in memory
   if (gDebug > 3) {
      TIter next(fHostAuth->Established());
      TRootSecContext *ctx;
      while ((ctx = (TRootSecContext *) next()))
         ctx->Print("0");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called in connection with a timer timeout

void TAuthenticate::CatchTimeOut()
{
   Info("CatchTimeOut", "%d sec timeout expired (protocol: %s)",
        fgAuthTO, fgAuthMeth[fSecurity].Data());

   fTimeOut = 1;
   if (fSocket)
      fSocket->Close("force");

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Authenticate to remote rootd or proofd server. Return kTRUE if
/// authentication succeeded.

Bool_t TAuthenticate::Authenticate()
{
   if (gDebug > 2 && gAuthenticateMutex)
      Info("Authenticate", "locking mutex (pid:  %d)",gSystem->GetPid());
   R__LOCKGUARD2(gAuthenticateMutex);

   Bool_t rc = kFALSE;
   Int_t st = -1;
   Int_t remMeth = 0, rMth[kMAXSEC], tMth[kMAXSEC] = {0};
   Int_t meth = 0;
   char noSupport[80] = { 0 };
   char triedMeth[80] = { 0 };
   Int_t ntry = 0;

   TString user, passwd;
   Bool_t pwhash;

   if (gDebug > 2)
      Info("Authenticate", "enter: fUser: %s", fUser.Data());

   //
   // Setup timeout timer, if required
   TTimer *alarm = 0;
   if (fgAuthTO > 0) {
      alarm = new TTimer(0, kFALSE);
      alarm->SetInterruptSyscalls();
      // The method CatchTimeOut will be called at timeout
      alarm->Connect("Timeout()", "TAuthenticate", this, "CatchTimeOut()");
   }

negotia:
   st = -1;
   tMth[meth] = 1;
   ntry++;
   if (gDebug > 2)
      Info("Authenticate", "try #: %d", ntry);

   user = "";
   passwd = "";
   pwhash = kFALSE;

   // Security level from the list (if not in cleanup mode ...)
   fSecurity = (ESecurity) fHostAuth->GetMethod(meth);
   fDetails = fHostAuth->GetDetails((Int_t) fSecurity);
   if (gDebug > 2)
      Info("Authenticate",
           "trying authentication: method:%d, default details:%s",
           fSecurity, fDetails.Data());

   // Keep track of tried methods in a list
   if (triedMeth[0] != '\0')
      (void) strlcat(triedMeth, " ", sizeof(triedMeth) - 1);

   (void) strlcat(triedMeth, fgAuthMeth[fSecurity].Data(), sizeof(triedMeth) - 1);

   // Set environments
   SetEnvironment();

   st = -1;

   //
   // Reset timeout variables and start timer
   fTimeOut = 0;
   if (fgAuthTO > 0 && alarm) {
      alarm->Start(fgAuthTO*1000, kTRUE);
   }

   // Auth calls depend of fSec
   if (fSecurity == kClear) {

      rc = kFALSE;

      // UsrPwd Authentication
      user = fgDefaultUser;
      if (user != "")
         CheckNetrc(user, passwd, pwhash, kFALSE);
      if (passwd == "") {
         if (fgPromptUser) {
            char *u = PromptUser(fRemote);
            user = u;
            delete[] u;
         }
         rc = GetUserPasswd(user, passwd, pwhash, kFALSE);
      }
      fUser = user;
      fPasswd = passwd;

      if (!rc) {

         if (fUser != "root")
            st = ClearAuth(user, passwd, pwhash);
      } else {
         Error("Authenticate",
               "unable to get user name for UsrPwd authentication");
      }

   }

   // Stop timer
   if (alarm) alarm->Stop();

   // Flag timeout condition
   st = (fTimeOut > 0) ? -3 : st;

   //
   // Analyse the result now ...
   // Type of action after the analysis:
   // 0 = return, 1 = negotiation, 2 = send kROOTD_BYE + 3,
   // 3 = print failure and return
   Int_t action = 0;
   Int_t nmet = fHostAuth->NumMethods();
   Int_t remloc = nmet - ntry;
   if (gDebug > 0)
      Info("Authenticate","remloc: %d, ntry: %d, meth: %d, fSecurity: %d",
                           remloc, ntry, meth, fSecurity);
   Int_t kind, stat;
   switch (st) {

   case 1:
      //
      // Success
      fHostAuth->CountSuccess((Int_t)fSecurity);
      if (gDebug > 2)
         fSecContext->Print();
      if (fSecContext->IsActive())
         fSecContext->AddForCleanup(fSocket->GetPort(),
                                    fSocket->GetRemoteProtocol(),fSocket->GetServType());
      rc = kTRUE;
      break;

   case 0:
      //
      // Failure
      fHostAuth->CountFailure((Int_t)fSecurity);
      if (fVersion < 2) {
         //
         // Negotiation not supported by old daemons ...
         if (gDebug > 2)
            Info("Authenticate",
                 "negotiation not supported remotely: try next method, if any");
         if (meth < nmet - 1) {
            meth++;
            action = 1;
         } else {
            action = 2;
         }
         rc = kFALSE;
         break;
      }
      //
      // Attempt negotiation ...
      if (fSocket->Recv(stat, kind) < 0) {
         action = 0;
         rc = kFALSE;
      }
      if (gDebug > 2)
         Info("Authenticate",
              "after failed attempt: kind= %d, stat= %d", kind, stat);
      if (kind == kROOTD_ERR) {
         action = 2;
         rc = kFALSE;
      } else if (kind == kROOTD_NEGOTIA) {
         if (stat > 0) {
            int len = 3 * stat;
            char *answer = new char[len];
            int nrec = fSocket->Recv(answer, len, kind);  // returns user
            if (nrec < 0) {
               delete[] answer; // delete buffer while it exit switch() scope
               action = 0;
               rc = kFALSE;
               break;
            }
            if (kind != kMESS_STRING)
               Warning("Authenticate",
                       "strings with accepted methods not received (%d:%d)",
                       kind, nrec);
            remMeth =
               sscanf(answer, "%d %d %d %d %d %d", &rMth[0], &rMth[1],
                      &rMth[2], &rMth[3], &rMth[4], &rMth[5]);
            if (gDebug > 0 && remloc > 0)
               Info("Authenticate",
                    "remotely allowed methods not yet tried: %s",
                    answer);
            delete[] answer;
         } else if (stat == 0) {
            Info("Authenticate",
                 "no more methods accepted remotely to be tried");
            action = 3;
            rc = kFALSE;
            break;
         }
         // If no more local methods, return
         if (remloc < 1) {
            action = 2;
            rc = kFALSE;
            break;
         }
         // Look if a non-tried method matches
         int i, j;
         std::string available{};
         Bool_t methfound = kFALSE;
         for (i = 0; i < remMeth; i++) {
            for (j = 0; j < nmet; j++) {
               if (fHostAuth->GetMethod(j) == rMth[i] && tMth[j] == 0) {
                  meth = j;
                  action = 1;
                  methfound = kTRUE;
                  break;
               }
               if (i == 0)
                  available += " " + std::to_string(fHostAuth->GetMethod(j));
            }
            if (methfound) break;
         }
         if (methfound) break;
         //
         // No method left to be tried: notify and exit
         if (gDebug > 0)
            Warning("Authenticate", "no match with those locally available: %s", available.c_str());
         action = 2;
         rc = kFALSE;
         break;
      } else {        // unknown message code at this stage
         action = 3;
         rc = kFALSE;
         break;
      }
      break;

   case -1:
      //
      // Method not supported
      fHostAuth->CountFailure((Int_t)fSecurity);
      if (gDebug > 2)
         Info("Authenticate",
              "method not even started: insufficient or wrong info: %s",
              "try with next method, if any");
      fHostAuth->RemoveMethod(fSecurity);
      nmet--;
      if (nmet > 0) {
         action = 1;
      } else
         action = 2;

      break;

   case -2:
      //
      // Remote host does not accepts connections from local host
      fHostAuth->CountFailure((Int_t)fSecurity);
      if (fVersion <= 2)
         if (gDebug > 2)
            Warning("Authenticate",
                    "status code -2 not expected from old daemons");
      rc = kFALSE;
      break;

   case -3:
      //
      // Timeout: we set the method as last one, should the caller
      // decide to retry, if it will attempt first something else.
      // (We can not retry directly, because the server will not be
      //  synchronized ...)
      fHostAuth->CountFailure((Int_t)fSecurity);
      if (gDebug > 2)
         Info("Authenticate", "got a timeout");
      fHostAuth->SetLast(fSecurity);
      if (meth < nmet - 1) {
         fTimeOut = 2;
      } else
         fTimeOut = 1;
      rc = kFALSE;
      break;

   default:
      fHostAuth->CountFailure((Int_t)fSecurity);
      if (gDebug > 2)
         Info("Authenticate", "unknown status code: %d - assume failure",st);
      rc = kFALSE;
      action = 0;
      break;
   }

   switch (action) {
   case 1:
      goto negotia;
      // No break but we go away anyhow
   case 2:
      fSocket->Send("0", kROOTD_BYE);
      // fallthrough
   case 3:
      if (strlen(noSupport) > 0)
         Info("Authenticate", "attempted methods %s are not supported"
              " by remote server version", noSupport);
      Info("Authenticate",
           "failure: list of attempted methods: %s", triedMeth);
      AuthError("Authenticate",-1);
      rc = kFALSE;
      break;
   default:
      break;
   }

   // Cleanup timer
   SafeDelete(alarm);

   return rc;

}

////////////////////////////////////////////////////////////////////////////////
/// Set default authentication environment. The values are inferred
/// from fSecurity and fDetails.

void TAuthenticate::SetEnvironment()
{
   R__LOCKGUARD2(gAuthenticateMutex);

   if (gDebug > 2)
      Info("SetEnvironment",
           "setting environment: fSecurity:%d, fDetails:%s", fSecurity,
           fDetails.Data());

   // Defaults
   fgDefaultUser = fgUser;
   fgAuthReUse = kTRUE;
   fgPromptUser = kFALSE;

   // Decode fDetails, is non empty ...
   if (fDetails != "") {
      char usdef[kMAXPATHLEN] = { 0 };
      char pt[5] = { 0 }, ru[5] = { 0 };
      Int_t hh = 0, mm = 0;
      char us[kMAXPATHLEN] = {0}, cp[kMAXPATHLEN] = {0};
      const char *ptr;

      TString usrPromptDef = TString(GetAuthMethod(fSecurity)) + ".LoginPrompt";
      if ((ptr = strstr(fDetails, "pt:")) != 0) {
         sscanf(ptr + 3, "%4s %8191s", pt, usdef);
      } else {
         if (!strncasecmp(gEnv->GetValue(usrPromptDef,""),"no",2) ||
             !strncmp(gEnv->GetValue(usrPromptDef,""),"0",1))
            strncpy(pt,"0",2);
         else
            strncpy(pt,"1",2);
      }
      TString usrReUseDef = TString(GetAuthMethod(fSecurity)) + ".ReUse";
      if ((ptr = strstr(fDetails, "ru:")) != 0) {
         sscanf(ptr + 3, "%4s %8191s", ru, usdef);
      } else {
         if (!strncasecmp(gEnv->GetValue(usrReUseDef,""),"no",2) ||
             !strncmp(gEnv->GetValue(usrReUseDef,""),"0",1))
            strncpy(ru,"0",2);
         else
            strncpy(ru,"1",2);
      }
      TString usrValidDef = TString(GetAuthMethod(fSecurity)) + ".Valid";
      TString hours(gEnv->GetValue(usrValidDef,"24:00"));
      Int_t pd = 0;
      if ((pd = hours.Index(":")) > -1) {
         TString minutes = hours;
         hours.Resize(pd);
         minutes.Replace(0,pd+1,"");
         hh = atoi(hours.Data());
         mm = atoi(minutes.Data());
      } else {
         hh = atoi(hours.Data());
         mm = 0;
      }

      // Now action depends on method ...
      if (fSecurity == kClear) {
         if ((ptr = strstr(fDetails, "us:")) != 0)
            sscanf(ptr + 3, "%8191s %8191s", us, usdef);
         if ((ptr = strstr(fDetails, "cp:")) != 0)
            sscanf(ptr + 3, "%8191s %8191s", cp, usdef);
         if (gDebug > 2)
            Info("SetEnvironment", "details:%s, pt:%s, ru:%s, us:%s cp:%s",
                 fDetails.Data(), pt, ru, us, cp);
      } else {
         if ((ptr = strstr(fDetails, "us:")) != 0)
            sscanf(ptr + 3, "%8191s %8191s", us, usdef);
         if (gDebug > 2)
            Info("SetEnvironment", "details:%s, pt:%s, ru:%s, us:%s",
                 fDetails.Data(), pt, ru, us);
      }

      // Set Prompt flag
      if (!strncasecmp(pt, "yes",3) || !strncmp(pt, "1", 1))
         fgPromptUser = kTRUE;

      // Set ReUse flag
      if (!gROOT->IsProofServ()) {
         fgAuthReUse = kTRUE;
         if (!strncasecmp(ru, "no",2) || !strncmp(ru, "0",1))
            fgAuthReUse = kFALSE;
      }

      // Set Expiring date
      fgExpDate = TDatime();
      fgExpDate.Set(fgExpDate.Convert() + hh*3600 + mm*60);

      // UnSet Crypt flag for UsrPwd, if requested
      if (fSecurity == kClear) {
         fgUsrPwdCrypt = kTRUE;
         if (!strncmp(cp, "no", 2) || !strncmp(cp, "0", 1))
            fgUsrPwdCrypt = kFALSE;
      }
      // Build UserDefaults
      usdef[0] = '\0';
      // give highest priority to command-line specification
      if (fUser == "") {
         if (strlen(us) > 0) snprintf(usdef, kMAXPATHLEN, "%s", us);
      } else {
         snprintf(usdef, kMAXPATHLEN, "%s", fUser.Data());
      }

      if (strlen(usdef) > 0) {
         fgDefaultUser = usdef;
      } else {
         if (fgUser != "") {
            fgDefaultUser = fgUser;
         } else {
            UserGroup_t *u = gSystem->GetUserInfo();
            if (u)
               fgDefaultUser = u->fUser;
            delete u;
         }
      }
      if (fgDefaultUser == "anonymous" || fgDefaultUser == "rootd" ||
          fgUser != "" || fUser != "") {
         // when set by user don't prompt for it anymore
         fgPromptUser = kFALSE;
      }

      if (gDebug > 2)
         Info("SetEnvironment", "usdef:%s", fgDefaultUser.Data());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Try to get user name and passwd from several sources.

Bool_t TAuthenticate::GetUserPasswd(TString &user, TString &passwd,
                                    Bool_t &pwhash, Bool_t srppwd)
{
   if (srppwd) {
      Error("GetUserPasswd", "SRP no longer supported by ROOT");
      return 1;
   }

   if (gDebug > 3)
      Info("GetUserPasswd", "Enter: User: '%s' Hash:%d SRP:%d",
           user.Data(),(Int_t)pwhash,(Int_t)false);

   // Get user and passwd set via static functions SetUser and SetPasswd.
   if (user == "" && fgUser != "")
      user = fgUser;

   if (fgUser != "" && user == fgUser) {
      if (passwd == "" && fgPasswd != "") {
         passwd = fgPasswd;
         pwhash = fgPwHash;
      }
   }

   if (gDebug > 3)
      Info("GetUserPasswd", "In memory: User: '%s' Hash:%d",
           user.Data(),(Int_t)pwhash);

   // Check system info for user if still not defined
   if (user == "") {
      UserGroup_t *u = gSystem->GetUserInfo();
      if (u)
         user = u->fUser;
      delete u;
      if (gDebug > 3)
         Info("GetUserPasswd", "In memory: User: '%s' Hash:%d",
              user.Data(),(Int_t)pwhash);
   }

   // Check ~/.rootnetrc and ~/.netrc files if user was not set via
   // the static SetUser() method.
   if (user == "" || passwd == "") {
      if (gDebug > 3)
         Info("GetUserPasswd", "Checking .netrc family ...");
      CheckNetrc(user, passwd, pwhash, /* srppwd */ false);
   }
   if (gDebug > 3)
      Info("GetUserPasswd", "From .netrc family: User: '%s' Hash:%d",
           user.Data(),(Int_t)pwhash);

   // If user also not set via  ~/.rootnetrc or ~/.netrc ask user.
   if (user == "") {
      char *p = PromptUser(fRemote);
      user = p;
      delete [] p;
      if (user == "") {
         Error("GetUserPasswd", "user name not set");
         return 1;
      }
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Try to get user name and passwd from the ~/.rootnetrc or
/// ~/.netrc files. For more info see the version with 4 arguments.
/// This version is maintained for backward compatability reasons.

Bool_t TAuthenticate::CheckNetrc(TString &user, TString &passwd)
{
   Bool_t hash = false;
   return CheckNetrc(user, passwd, hash, /* srppwd */ false);
}

////////////////////////////////////////////////////////////////////////////////
/// Try to get user name and passwd from the ~/.rootnetrc or
/// ~/.netrc files. First ~/.rootnetrc is tried, after that ~/.netrc.
/// These files will only be used when their access masks are 0600.
/// Returns kTRUE if user and passwd were found for the machine
/// specified in the URL. If kFALSE, user and passwd are "".
/// The boolean pwhash is set to kTRUE if the returned passwd is to
/// be understood as password hash, i.e. if the 'password-hash' keyword
/// is found in the 'machine' lines; not implemented for 'secure'
/// and the .netrc file.
/// The format of these files are:
///
/// # this is a comment line
/// machine <machine fqdn> login <user> password <passwd>
/// machine <machine fqdn> login <user> password-hash <passwd>
///
/// and in addition ~/.rootnetrc also supports:
///
/// secure <machine fqdn> login <user> password <passwd>
///
/// <machine fqdn> may be a domain name or contain the wild card '*'.
///
/// for the secure protocols. All lines must start in the first column.

Bool_t TAuthenticate::CheckNetrc(TString &user, TString &passwd,
                                 Bool_t &pwhash, Bool_t srppwd)
{
   if (srppwd) {
      Error("CheckNetrc", "SRP no longer supported by ROOT");
      return 1;
   }

   Bool_t result = kFALSE;
   Bool_t first = kTRUE;
   TString remote = fRemote;

   passwd = "";
   pwhash = kFALSE;

   char *net =
      gSystem->ConcatFileName(gSystem->HomeDirectory(), ".rootnetrc");

   // Determine FQDN of the host ...
   TInetAddress addr = gSystem->GetHostByName(fRemote);
   if (addr.IsValid())
      remote = addr.GetHostName();

again:
   // Only use file when its access rights are 0600
   FileStat_t buf;
   if (gSystem->GetPathInfo(net, buf) == 0) {
#ifdef WIN32
      // Since Win32 does not have proper protections use file always
      bool mode0600 = true;
#else
      bool mode0600 = (buf.fMode & 0777) == (kS_IRUSR | kS_IWUSR);
#endif
      if (R_ISREG(buf.fMode) && !R_ISDIR(buf.fMode) && mode0600) {
            FILE *fd = fopen(net, "r");
            char line[256];
            while (fgets(line, sizeof(line), fd) != 0) {
               if (line[0] == '#')
                  continue;
               char word[6][64];
               int nword = sscanf(line, "%63s %63s %63s %63s %63s %63s",
                                        word[0], word[1], word[2], word[3], word[4], word[5]);
               if (nword != 6)
                  continue;
               if (strcmp(word[0], "machine"))
                  continue;
               if (strcmp(word[2], "login"))
                  continue;
               if (strcmp(word[4], "password") && strcmp(word[4], "password-hash"))
                  continue;

               // Treat the host name found in file as a regular expression
               // with '*' as a wild card
               TString href(word[1]);
               href.ReplaceAll("*",".*");
               TRegexp rg(href);
               if (remote.Index(rg) != kNPOS) {
                  if (user == "") {
                     user = word[3];
                     passwd = word[5];
                     if (!strcmp(word[4], "password-hash"))
                        pwhash = kTRUE;
                     result = kTRUE;
                     break;
                  } else {
                     if (!strcmp(word[3], user.Data())) {
                        passwd = word[5];
                        if (!strcmp(word[4], "password-hash"))
                           pwhash = kTRUE;
                        result = kTRUE;
                        break;
                     }
                  }
               }
            }
            fclose(fd);
         } else
            Warning("CheckNetrc",
                    "file %s exists but has not 0600 permission", net);
      }
      delete [] net;

      if (first && !result) {
         net = gSystem->ConcatFileName(gSystem->HomeDirectory(), ".netrc");
         first = kFALSE;
         goto again;
      }

      return result;
   }

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the global user.

const char *TAuthenticate::GetGlobalUser()
{
   return fgUser;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the global password hash flag.

Bool_t TAuthenticate::GetGlobalPwHash()
{
   return fgPwHash;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the global SRP password flag.

Bool_t TAuthenticate::GetGlobalSRPPwd()
{
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning default expiring date for new validity contexts

TDatime TAuthenticate::GetGlobalExpDate()
{
   return fgExpDate;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the default user information.

const char *TAuthenticate::GetDefaultUser()
{
   return fgDefaultUser;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the principal to be used to init Krb5 tickets.

const char *TAuthenticate::GetKrb5Principal()
{
   ::Error("Krb5Auth", "Kerberos5 is no longer supported by ROOT");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the authentication reuse settings.

Bool_t TAuthenticate::GetAuthReUse()
{
   return fgAuthReUse;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the prompt user settings.

Bool_t TAuthenticate::GetPromptUser()
{
   return fgPromptUser;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the method corresponding to idx.

const char *TAuthenticate::GetAuthMethod(Int_t idx)
{
   R__LOCKGUARD2(gAuthenticateMutex);

   if (idx < 0 || idx > kMAXSEC-1) {
      ::Error("Authenticate::GetAuthMethod", "idx out of bounds (%d)", idx);
      idx = 0;
   }
   return fgAuthMeth[idx];
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the method index (which can be used to find
/// the method in GetAuthMethod()). Returns -1 in case meth is not found.

Int_t TAuthenticate::GetAuthMethodIdx(const char *meth)
{
   R__LOCKGUARD2(gAuthenticateMutex);

   if (meth && meth[0]) {
      for (Int_t i = 0; i < kMAXSEC; i++) {
         if (!fgAuthMeth[i].CompareTo(meth, TString::kIgnoreCase))
            return i;
      }
   }

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to prompt for the user name to be used for authentication
/// to rootd or proofd. User is asked to type user name.
/// Returns user name (which must be deleted by caller) or 0.
/// If non-interactive run (eg ProofServ) returns default user.

char *TAuthenticate::PromptUser(const char *remote)
{
   R__LOCKGUARD2(gAuthenticateMutex);

   const char *user;
   if (fgDefaultUser != "")
      user = fgDefaultUser;
   else
      user = gSystem->Getenv("USER");
#ifdef R__WIN32
   if (!user)
      user = gSystem->Getenv("USERNAME");
#endif
   if (isatty(0) == 0 || isatty(1) == 0) {
      ::Warning("TAuthenticate::PromptUser",
                "not tty: cannot prompt for user, returning default");
      if (strlen(user))
         return StrDup(user);
      else
         return StrDup("None");
   }

   const char *usrIn = Getline(Form("Name (%s:%s): ", remote, user));
   if (usrIn[0]) {
      TString usr(usrIn);
      usr.Remove(usr.Length() - 1); // get rid of \n
      if (!usr.IsNull())
         return StrDup(usr);
      else
         return StrDup(user);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to prompt for the user's passwd to be used for
/// authentication to rootd or proofd. Uses non-echoing command line
/// to get passwd. Returns passwd (which must de deleted by caller) or 0.
/// If non-interactive run (eg ProofServ) returns -1

char *TAuthenticate::PromptPasswd(const char *prompt)
{
   if (isatty(0) == 0 || isatty(1) == 0) {
      ::Warning("TAuthenticate::PromptPasswd",
                "not tty: cannot prompt for passwd, returning -1");
      static char noint[4] = {"-1"};
      return StrDup(noint);
   }

   char buf[128];
   const char *pw = buf;
   // Get the plugin for the passwd dialog box, if needed
   if (!gROOT->IsBatch() && (fgPasswdDialog == (TPluginHandler *)(-1)) &&
       gEnv->GetValue("Auth.UsePasswdDialogBox", 1) == 1) {
      if ((fgPasswdDialog =
           gROOT->GetPluginManager()->FindHandler("TGPasswdDialog"))) {
         if (fgPasswdDialog->LoadPlugin() == -1) {
            fgPasswdDialog = 0;
            ::Warning("TAuthenticate",
                      "could not load plugin for the password dialog box");
         }
      }
   }
   if (fgPasswdDialog && (fgPasswdDialog != (TPluginHandler *)(-1))) {

      // Use graphic dialog
      fgPasswdDialog->ExecPlugin(3, prompt, buf, 128);

      // Wait until the user is done
      while (gROOT->IsInterrupted())
         gSystem->DispatchOneEvent(kFALSE);

   } else {
      Gl_config("noecho", 1);
      pw = Getline(prompt);
      Gl_config("noecho", 0);
   }

   // Final checks
   if (pw[0]) {
      TString spw(pw);
      if (spw.EndsWith("\n"))
         spw.Remove(spw.Length() - 1);   // get rid of \n
      char *rpw = StrDup(spw);
      return rpw;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the globus authorization hook (no longer supported)

GlobusAuth_t TAuthenticate::GetGlobusAuthHook()
{
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the RSA public keys.

const char *TAuthenticate::GetRSAPubExport(Int_t key)
{
   key = (key >= 0 && key <= 1) ? key : 0;
   return fgRSAPubExport[key].keys;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the RSA initialization flag.

Int_t TAuthenticate::GetRSAInit()
{
   return fgRSAInit;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method setting the default type of RSA key.

void TAuthenticate::SetDefaultRSAKeyType(Int_t key)
{
   if (key >= 0 && key <= 1)
      fgRSAKey = key;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method setting RSA initialization flag.

void TAuthenticate::SetRSAInit(Int_t init)
{
   fgRSAInit = init;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the list with authentication details.

TList *TAuthenticate::GetAuthInfo()
{
   R__LOCKGUARD2(gAuthenticateMutex);

   if (!fgAuthInfo)
      fgAuthInfo = new TList;
   return fgAuthInfo;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the list with authentication directives
/// to be sent to proof.

TList *TAuthenticate::GetProofAuthInfo()
{
   R__LOCKGUARD2(gAuthenticateMutex);

   if (!fgProofAuthInfo)
      fgProofAuthInfo = new TList;
   return fgProofAuthInfo;
}

////////////////////////////////////////////////////////////////////////////////
/// Print error string depending on error code.

void TAuthenticate::AuthError(const char *where, Int_t err)
{
   R__LOCKGUARD2(gAuthenticateMutex);

   // Make sure it is in range
   err = (err < kErrError) ? ((err > -1) ? err : -1) : kErrError;

   Int_t erc = err;
   Bool_t forceprint = kFALSE;
   TString lasterr = "";
   if (err == -1) {
      forceprint = kTRUE;
      erc = fgLastError;
      lasterr = "(last error only; re-run with gDebug > 0 for more details)";
   }

   if (erc > -1)
      if (gDebug > 0 || forceprint) {
         if (gRootdErrStr[erc])
            ::Error(Form("TAuthenticate::%s", where), "%s %s",
                    gRootdErrStr[erc], lasterr.Data());
         else
            ::Error(Form("TAuthenticate::%s", where),
                    "unknown error code: server must be running a newer ROOT version %s",
                    lasterr.Data());
      }

   // Update last error code
   fgLastError = err;
}

////////////////////////////////////////////////////////////////////////////////
/// Set global user name to be used for authentication to rootd or proofd.

void TAuthenticate::SetGlobalUser(const char *user)
{
   R__LOCKGUARD2(gAuthenticateMutex);

   if (fgUser != "")
      fgUser = "";

   if (user && user[0])
      fgUser = user;
}

////////////////////////////////////////////////////////////////////////////////
/// Set global passwd to be used for authentication to rootd or proofd.

void TAuthenticate::SetGlobalPasswd(const char *passwd)
{
   R__LOCKGUARD2(gAuthenticateMutex);

   if (fgPasswd != "")
      fgPasswd = "";

   if (passwd && passwd[0])
      fgPasswd = passwd;
}

////////////////////////////////////////////////////////////////////////////////
/// Set global passwd hash flag to be used for authentication to rootd or proofd.

void TAuthenticate::SetGlobalPwHash(Bool_t pwhash)
{
   fgPwHash = pwhash;
}

////////////////////////////////////////////////////////////////////////////////
/// Set global SRP passwd flag to be used for authentication to rootd or proofd.

void TAuthenticate::SetGlobalSRPPwd(Bool_t)
{
   ::Error("SetGlobalSRPPwd", "SRP no longer supported by ROOT");
}

////////////////////////////////////////////////////////////////////////////////
/// Set flag controlling the reading of $HOME/.rootauthrc.
/// In PROOF the administrator may want to switch off private settings.
/// Always true, may only be set false via option to proofd.

void TAuthenticate::SetReadHomeAuthrc(Bool_t readhomeauthrc)
{
   fgReadHomeAuthrc = readhomeauthrc;
}

////////////////////////////////////////////////////////////////////////////////
/// Set default expiring date for new validity contexts

void TAuthenticate::SetGlobalExpDate(TDatime expdate)
{
   fgExpDate = expdate;
}

////////////////////////////////////////////////////////////////////////////////
/// Set default user name.

void TAuthenticate::SetDefaultUser(const char *defaultuser)
{
   if (fgDefaultUser != "")
      fgDefaultUser = "";

   if (defaultuser && defaultuser[0])
      fgDefaultUser = defaultuser;
}

////////////////////////////////////////////////////////////////////////////////
/// Set timeout (active if > 0)

void TAuthenticate::SetTimeOut(Int_t to)
{
   fgAuthTO = (to <= 0) ? -1 : to;
}

////////////////////////////////////////////////////////////////////////////////
/// Set global AuthReUse flag

void TAuthenticate::SetAuthReUse(Bool_t authreuse)
{
   fgAuthReUse = authreuse;
}

////////////////////////////////////////////////////////////////////////////////
/// Set global PromptUser flag

void TAuthenticate::SetPromptUser(Bool_t promptuser)
{
   fgPromptUser = promptuser;
}

////////////////////////////////////////////////////////////////////////////////
/// Set secure authorization function.

void TAuthenticate::SetSecureAuthHook(SecureAuth_t func)
{
   fgSecAuthHook = func;
}

////////////////////////////////////////////////////////////////////////////////
/// Set kerberos5 authorization function. Automatically called when
/// libKrb5Auth is loaded.

void TAuthenticate::SetKrb5AuthHook(Krb5Auth_t)
{
   ::Error("Krb5Auth", "Kerberos5 is no longer supported by ROOT");
}

////////////////////////////////////////////////////////////////////////////////
/// Set Globus authorization function. Automatically called when
/// libGlobusAuth is loaded.

void TAuthenticate::SetGlobusAuthHook(GlobusAuth_t)
{
   ::Error("GlobusAuth", "Globus is no longer supported by ROOT");
}

////////////////////////////////////////////////////////////////////////////////
/// SSH client authentication code (no longer supported)

Int_t TAuthenticate::SshAuth(TString & /* user */)
{
   ::Error("SshAuth", "SSH is no longer supported by ROOT");
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Method returning the user to be used for the ssh login (no longer supported)

const char *TAuthenticate::GetSshUser(TString /* user */) const
{
   ::Error("GetSshUser", "SSH is no longer supported by ROOT");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if 'host' matches 'href':
/// this means either equal or "containing" it, even with wild cards *
/// in the first field (in the case 'href' is a name, ie not IP address)
/// Returns kTRUE if the two matches.

Bool_t TAuthenticate::CheckHost(const char *host, const char *href)
{
   R__LOCKGUARD2(gAuthenticateMutex);

   Bool_t retval = kTRUE;

   // Both strings should have been defined
   if (!host || !href)
      return kFALSE;

   // 'href' == '*' indicates any 'host' ...
   if (!strcmp(href,"*"))
      return kTRUE;

   // If 'href' contains at a letter or an hyphen it is assumed to be
   // a host name. Otherwise a name.
   // Check also for wild cards
   Bool_t name = kFALSE;
   TRegexp rename("[+a-zA-Z]");
   Int_t len;
   if (rename.Index(href,&len) != -1 || strstr(href,"-"))
      name = kTRUE;

   // Check also for wild cards
   Bool_t wild = kFALSE;
   if (strstr(href,"*"))
      wild = kTRUE;

   // Now build the regular expression for final checking
   TRegexp rehost(href,wild);

   // host to check
   TString theHost(host);
   if (!name) {
      TInetAddress addr = gSystem->GetHostByName(host);
      theHost = addr.GetHostAddress();
      if (gDebug > 2)
         ::Info("TAuthenticate::CheckHost", "checking host IP: %s", theHost.Data());
   }

   // Check 'host' against 'rehost'
   Ssiz_t pos = rehost.Index(theHost,&len);
   if (pos == -1)
      retval = kFALSE;

   // If IP and no wilds, it should match either
   // the beginning or the end of the string
   if (!wild) {
      if (pos > 0 && pos != (Ssiz_t)(theHost.Length()-strlen(href)))
         retval = kFALSE;
   }

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// RFIO authentication (no longer supported)

Int_t TAuthenticate::RfioAuth(TString &)
{
   ::Error("RfioAuth", "RfioAuth is no longer supported by ROOT");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// UsrPwd client authentication code.
/// Returns 0 in case authentication failed
///         1 in case of success

Int_t TAuthenticate::ClearAuth(TString &user, TString &passwd, Bool_t &pwdhash)
{
   R__LOCKGUARD2(gAuthenticateMutex);

   if (gDebug > 2)
      Info("ClearAuth", "enter: user: %s (passwd hashed?: %d)",
           user.Data(),(Int_t)pwdhash);

   Int_t reuse    = fgAuthReUse;
   Int_t prompt   = fgPromptUser;
   Int_t cryptopt = fgUsrPwdCrypt;
   Int_t needsalt = 1;
   if (pwdhash)
      needsalt = 0;
   fDetails = TString::Format("pt:%d ru:%d cp:%d us:",
                              fgPromptUser, fgAuthReUse, fgUsrPwdCrypt) + user;
   if (gDebug > 2)
      Info("ClearAuth", "ru:%d pt:%d cp:%d ns:%d rk:%d",
           fgAuthReUse,fgPromptUser,fgUsrPwdCrypt,needsalt,fgRSAKey);
#ifdef R__WIN32
   needsalt = 0;
#endif
   Int_t stat, kind;

   if (fVersion > 1) {

      //
      // New protocol
      //
      Int_t anon = 0;
      TString salt = "";
      TString pashash = "";

      // Get effective user (fro remote checks in $HOME/.rhosts)
      UserGroup_t *pw = gSystem->GetUserInfo(gSystem->GetEffectiveUid());
      TString effUser;
      if (pw) {
         effUser = TString(pw->fUser);
         delete pw;
      } else
         effUser = user;

      // Create options string
      int opt = (reuse * kAUTH_REUSE_MSK) + (cryptopt * kAUTH_CRYPT_MSK) +
         (needsalt * kAUTH_SSALT_MSK) + (fRSAKey * kAUTH_RSATY_MSK);
      TString options;
      options.Form("%d %ld %s %ld %s", opt,
                   (Long_t)user.Length(), user.Data(),
                   (Long_t)effUser.Length(), effUser.Data());

      // Check established authentications
      kind = kROOTD_USER;
      stat = reuse;
      Int_t rc = 0;
      if ((rc = AuthExists(user, (Int_t) TAuthenticate::kClear, options,
                           &kind, &stat, &StdCheckSecCtx)) == 1) {
         // A valid authentication exists: we are done ...
         return 1;
      }
      if (rc == -2) {
         return rc;
      }
      if (stat == kErrNotAllowed && kind == kROOTD_ERR) {
         return 0;
      }

      if (kind == kROOTD_AUTH && stat == -1) {
         if (gDebug > 3)
            Info("ClearAuth", "anonymous user");
         anon  = 1;
         cryptopt = 0;
         reuse = 0;
         needsalt = 0;
      }

      // The random tag in hex representation
      // Protection against reply attacks
      char ctag[11] = {0};
      if (anon == 0 && cryptopt == 1) {

         // Check that we got the right thing ..
         if (kind != kROOTD_RSAKEY || stat < 1 || stat > 2 ) {
            // Check for errors
            if (kind != kROOTD_ERR) {
               Warning("ClearAuth",
                       "problems recvn RSA key flag: got message %d, flag: %d",
                       kind, stat);
            }
            return 0;
         }
         if (gDebug > 3)
            Info("ClearAuth", "get key request ...");

         // Save type of key
         fRSAKey = stat - 1;

         // Send the key securely
         if (SendRSAPublicKey(fSocket,fRSAKey) < 0)
            return 0;

         int slen = 0;
         if (needsalt) {
            // Receive password salt
            char *tmpsalt = 0;
            if ((slen = SecureRecv(fSocket, 1, fRSAKey, &tmpsalt)) == -1) {
               Warning("ClearAuth", "problems secure-receiving salt -"
                       " may result in corrupted salt");
               Warning("ClearAuth", "switch off reuse for this session");
               needsalt = 0;
               return 0;
            }
            if (slen) {
               // Extract random tag, if there
               if (slen > 9) {
                  int ltmp = slen;
                  while (ltmp && tmpsalt[ltmp-1] != '#') ltmp--;
                  if (ltmp) {
                     if (tmpsalt[ltmp-1] == '#' &&
                         tmpsalt[ltmp-10] == '#') {
                        strlcpy(ctag,&tmpsalt[ltmp-10],11);
                        // We drop the random tag
                        ltmp -= 10;
                        tmpsalt[ltmp] = 0;
                        // Update salt length
                        slen -= 10;
                     }
                  }
                  if (!tmpsalt[0]) {
                     // No salt left
                     needsalt = 0;
                     slen = 0;
                  }
               }
               if (slen)
                  salt = TString(tmpsalt);
               delete [] tmpsalt;
            }
            if (gDebug > 2)
               Info("ClearAuth", "got salt: '%s' (len: %d)", salt.Data(), slen);
         } else {
            if (gDebug > 2)
               Info("ClearAuth", "Salt not required");
            char *tmptag = 0;
            if (SecureRecv(fSocket, 1, fRSAKey, &tmptag) == -1) {
               Warning("ClearAuth", "problems secure-receiving rndmtag -"
                       " may result in corrupted rndmtag");
            }
            if (tmptag) {
               strlcpy(ctag, tmptag, 11);
               delete [] tmptag;
            }
         }
         // We may not have got a salt (if the server may not access it
         // or if it needs the full password, like for AFS ...)
         if (!slen)
            needsalt = 0;
      }
      // Now get the password either from prompt or from memory, if saved already
      if (anon == 1) {

         if (fgPasswd.Contains("@")) {
            // Anonymous like login with user chosen passwd ...
            passwd = fgPasswd;
         } else {
            // Anonymous like login with automatic passwd generation ...
            TString localuser;
            pw = gSystem->GetUserInfo();
            if (pw) {
               char *u = StrDup(pw->fUser);
               localuser = u;
               delete[] u;
            }
            delete pw;
            static TString localFQDN;
            if (localFQDN == "") {
               TInetAddress addr = gSystem->GetHostByName(gSystem->HostName());
               if (addr.IsValid())
                  localFQDN = addr.GetHostName();
            }
            passwd.Form("%s@%s", localuser.Data(), localFQDN.Data());
            if (gDebug > 2)
               Info("ClearAuth",
                    "automatically generated anonymous passwd: %s",
                    passwd.Data());
         }

      } else {

         if (prompt == 1 || pashash.Length() == 0) {

            if (passwd == "") {
               TString xp;
               xp.Form("%s@%s password: ", user.Data(),fRemote.Data());
               char *pwd = PromptPasswd(xp);
               passwd = TString(pwd);
               delete [] pwd;
               if (passwd == "") {
                  Error("ClearAuth", "password not set");
                  fSocket->Send("-1", kROOTD_PASS);  // Needs this for consistency
                  return 0;
               }
            }
            if (needsalt && !pwdhash) {
#ifndef R__WIN32
               pashash = TString(crypt(passwd, salt));
               if (!pashash.BeginsWith(salt)) {
                  // not the right version of the crypt function:
                  // do not send hash
                  pashash = passwd;
               }
#else
               pashash = passwd;
#endif
            } else {
               pashash = passwd;
            }
         }

      }

      // Store password for later use
      fgUser = fUser;
      fgPwHash = kFALSE;
      fPwHash = kFALSE;
      fgPasswd = passwd;
      fPasswd = passwd;

      // Send it to server
      if (anon == 0 && cryptopt == 1) {

         // Needs to send this for consistency
         if (fSocket->Send("\0", kROOTD_PASS) < 0)
            return 0;

         // Add the random tag received from the server
         // (if any); makes packets non re-usable
         if (strlen(ctag))
            pashash += ctag;

         if (SecureSend(fSocket, 1, fRSAKey, pashash.Data()) == -1) {
            Warning("ClearAuth", "problems secure-sending pass hash"
                    " - may result in authentication failure");
            return 0;
         }
      } else {

         // Standard technique: invert passwd
         if (passwd != "") {
            for (int i = 0; i < passwd.Length(); i++) {
               char inv = ~passwd(i);
               passwd.Replace(i, 1, inv);
            }
         }
         if (fSocket->Send(passwd.Data(), kROOTD_PASS) < 0)
            return 0;
      }

      Int_t nrec = 0;
      // Receive username used for login
      if ((nrec = fSocket->Recv(stat, kind)) < 0 )  // returns user
         return 0;
      if (gDebug > 3)
         Info("ClearAuth", "after kROOTD_PASS: kind= %d, stat= %d", kind,
              stat);

      // Check for errors
      if (kind == kROOTD_ERR) {
         AuthError("ClearAuth", stat);
         fgPasswd = "";
         return 0;
      }

      if (kind != kROOTD_PASS || stat < 1)
         Warning("ClearAuth",
                 "problems recvn (user,offset) length (%d:%d bytes:%d)",
                 kind, stat, nrec);

      // Get user and offset
      char answer[256];
      int reclen = (stat+1 > 256) ? 256 : stat+1;
      if ((nrec = fSocket->Recv(answer, reclen, kind)) < 0)
         return 0;
      if (kind != kMESS_STRING)
         Warning("ClearAuth",
                 "username and offset not received (%d:%d)", kind,
                 nrec);

      // Parse answer
      char lUser[128];
      Int_t offset = -1;
      sscanf(answer, "%127s %d", lUser, &offset);
      if (gDebug > 3)
         Info("ClearAuth",
              "received from server: user: %s, offset: %d (%s)", lUser,
              offset, answer);

      // Return username
      user = lUser;

      char *token = 0;
      if (reuse == 1 && offset > -1) {
         // Receive token
         if (cryptopt == 1) {
            if (SecureRecv(fSocket, 1, fRSAKey, &token) == -1) {
               Warning("ClearAuth",
                       "problems secure-receiving token -"
                       " may result in corrupted token");
               return 0;
            }
         } else {
            Int_t tlen = 9;
            token = new char[tlen];
            if (fSocket->Recv(token, tlen, kind) < 0) {
               delete [] token;
               return 0;
            }
            if (kind != kMESS_STRING)
               Warning("ClearAuth", "token not received (%d:%d)", kind,
                       nrec);
            // Invert token
            for (int i = 0; i < (int) strlen(token); i++) {
               token[i] = ~token[i];
            }

         }
         if (gDebug > 3)
            Info("ClearAuth", "received from server: token: '%s' ",
                 token);
      }
      TPwdCtx *pwdctx = new TPwdCtx(fPasswd,fPwHash);
      // Create SecContext object
      fSecContext = fHostAuth->CreateSecContext((const char *)lUser, fRemote,
                                                kClear, offset, fDetails, (const char *)token,
                                                fgExpDate, (void *)pwdctx, fRSAKey);

      // Release allocated memory ...
      if (token)
         delete [] token;

      // This from remote login
      if (fSocket->Recv(stat, kind) < 0)
         return 0;


      if (kind == kROOTD_AUTH && stat >= 1) {
         if (stat == 5 && fSocket->GetServType() == TSocket::kPROOFD)
            // AFS: we cannot reuse the token because remotely the
            // daemon token must be re-initialized; for PROOF, we
            // just flag the entry as AFS; this allows to skip reusing
            // but to keep the session key for password forwarding
            fSecContext->SetID("AFS authentication");
         return 1;
      } else {
         fgPasswd = "";
         if (kind == kROOTD_ERR)
            AuthError("ClearAuth", stat);
         return 0;
      }

   } else {

      // Old Protocol

      // Send username
      if (fSocket->Send(user.Data(), kROOTD_USER) < 0)
         return 0;

      // Get replay from server
      if (fSocket->Recv(stat, kind) < 0)
         return 0;

      // This check should guarantee backward compatibility with a private
      // version of rootd used by CDF
      if (kind == kROOTD_AUTH && stat == 1) {
         fSecContext =
            fHostAuth->CreateSecContext(user,fRemote,kClear,-1,fDetails,0);
         return 1;
      }

      if (kind == kROOTD_ERR) {
         TString server = "sockd";
         if (fProtocol.Contains("root"))
            server = "rootd";
         if (fProtocol.Contains("proof"))
            server = "proofd";
         if (stat == kErrConnectionRefused) {
            if (gDebug > 0)
               Error("ClearAuth",
                     "%s@%s does not accept connections from %s@%s",
                     server.Data(),fRemote.Data(),
                     fUser.Data(),gSystem->HostName());
            return -2;
         } else if (stat == kErrNotAllowed) {
            if (gDebug > 0)
               Error("ClearAuth",
                     "%s@%s does not accept %s authentication from %s@%s",
                     server.Data(),fRemote.Data(),
                     TAuthenticate::fgAuthMeth[0].Data(),
                     fUser.Data(),gSystem->HostName());
         } else
            AuthError("ClearAuth", stat);
         return 0;
      }
      // Prepare passwd to send
   badpass1:
      if (passwd == "") {
         TString xp;
         xp.Form("%s@%s password: ", user.Data(),fRemote.Data());
         char *p = PromptPasswd(xp);
         passwd = p;
         delete [] p;
         if (passwd == "")
            Error("ClearAuth", "password not set");
      }
      if (fUser == "anonymous" || fUser == "rootd") {
         if (!passwd.Contains("@")) {
            Warning("ClearAuth",
                    "please use passwd of form: user@host.do.main");
            passwd = "";
            goto badpass1;
         }
      }

      fgPasswd = passwd;
      fPasswd = passwd;

      // Invert passwd
      if (passwd != "") {
         for (int i = 0; i < passwd.Length(); i++) {
            char inv = ~passwd(i);
            passwd.Replace(i, 1, inv);
         }
      }
      // Send it over the net
      if (fSocket->Send(passwd, kROOTD_PASS) < 0)
         return 0;

      // Get result of attempt
      if (fSocket->Recv(stat, kind) < 0)  // returns user
         return 0;
      if (gDebug > 3)
         Info("ClearAuth", "after kROOTD_PASS: kind= %d, stat= %d", kind,
              stat);

      if (kind == kROOTD_AUTH && stat == 1) {
         fSecContext =
            fHostAuth->CreateSecContext(user,fRemote,kClear,-1,fDetails,0);
         return 1;
      } else {
         if (kind == kROOTD_ERR)
            AuthError("ClearAuth", stat);
         return 0;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets fUser=user and search fgAuthInfo for the entry pertaining to
/// (host,user), setting fHostAuth accordingly.
/// If opt = "P" use fgProofAuthInfo list instead
/// If no entry is found fHostAuth is not changed

THostAuth *TAuthenticate::GetHostAuth(const char *host, const char *user,
                                      Option_t *opt, Int_t *exact)
{
   if (exact)
      *exact = 0;

   if (gDebug > 2)
      ::Info("TAuthenticate::GetHostAuth", "enter ... %s ... %s", host, user);

   // Strip off the servertype, if any
   Int_t srvtyp = -1;
   TString hostname = host;
   if (hostname.Contains(":")) {
      char *ps = (char *)strstr(host,":");
      if (ps)
         srvtyp = atoi(ps+1);
      hostname.Remove(hostname.Index(":"));
   }
   TString hostFQDN = hostname;
   if (strncmp(host,"default",7) && !hostFQDN.Contains("*")) {
      TInetAddress addr = gSystem->GetHostByName(hostFQDN);
      if (addr.IsValid())
         hostFQDN = addr.GetHostName();
   }
   TString usr = user;
   if (!usr.Length())
      usr = "*";
   THostAuth *rHA = 0;

   // Check list of auth info for already loaded info about this host
   TIter *next = new TIter(GetAuthInfo());
   if (!strncasecmp(opt,"P",1)) {
      SafeDelete(next);
      next = new TIter(GetProofAuthInfo());
   }

   THostAuth *ai;
   Bool_t notFound = kTRUE;
   Bool_t serverOK = kTRUE;
   while ((ai = (THostAuth *) (*next)())) {
      if (gDebug > 3)
         ai->Print("Authenticate::GetHostAuth");

      // server
      if (!(serverOK = (ai->GetServer() == -1) ||
            (ai->GetServer() == srvtyp)))
         continue;

      // Use default entry if existing and nothing more specific is found
      if (!strcmp(ai->GetHost(),"default") && serverOK && notFound)
         rHA = ai;

      // Check
      if (CheckHost(hostFQDN,ai->GetHost()) &&
          CheckHost(usr,ai->GetUser())     && serverOK) {
         rHA = ai;
         notFound = kFALSE;
      }

      if (hostFQDN == ai->GetHost() &&
          usr == ai->GetUser()     && srvtyp == ai->GetServer() ) {
         rHA = ai;
         if (exact)
            *exact = 1;
         break;
      }
   }
   SafeDelete(next);
   return rHA;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if a THostAuth with exact match for {host,user} exists
/// in the fgAuthInfo list
/// If opt = "P" use ProofAuthInfo list instead
/// Returns pointer to it or 0

THostAuth *TAuthenticate::HasHostAuth(const char *host, const char *user,
                                      Option_t *opt)
{
   if (gDebug > 2)
      ::Info("TAuthenticate::HasHostAuth", "enter ... %s ... %s", host, user);

   // Strip off the servertype, if any
   Int_t srvtyp = -1;
   TString hostFQDN = host;
   if (hostFQDN.Contains(":")) {
      char *ps = (char *)strstr(host,":");
      if (ps)
         srvtyp = atoi(ps+1);
      hostFQDN.Remove(hostFQDN.Index(":"));
   }
   if (strncmp(host,"default",7) && !hostFQDN.Contains("*")) {
      TInetAddress addr = gSystem->GetHostByName(hostFQDN);
      if (addr.IsValid())
         hostFQDN = addr.GetHostName();
   }

   TIter *next = new TIter(GetAuthInfo());
   if (!strncasecmp(opt,"P",1)) {
      SafeDelete(next);
      next = new TIter(GetProofAuthInfo());
   }
   THostAuth *ai;
   while ((ai = (THostAuth *) (*next)())) {

      if (hostFQDN == ai->GetHost() &&
          !strcmp(user, ai->GetUser()) && srvtyp == ai->GetServer()) {
         SafeDelete(next);
         return ai;
      }
   }
   SafeDelete(next);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Expands include directives found in fexp files
/// The expanded, temporary file, is pointed to by 'ftmp'
/// and should be already open. To be called recursively.

void TAuthenticate::FileExpand(const char *fexp, FILE *ftmp)
{
   FILE *fin;
   char line[kMAXPATHLEN];
   char cinc[20], fileinc[kMAXPATHLEN];

   if (gDebug > 2)
      ::Info("TAuthenticate::FileExpand", "enter ... '%s' ... 0x%lx", fexp, (Long_t)ftmp);

   fin = fopen(fexp, "r");
   if (fin == 0)
      return;

   while (fgets(line, sizeof(line), fin) != 0) {
      // Skip comment lines
      if (line[0] == '#')
         continue;
      if (line[strlen(line) - 1] == '\n')
         line[strlen(line) - 1] = '\0';
      if (gDebug > 2)
         ::Info("TAuthenticate::FileExpand", "read line ... '%s'", line);
      int nw = sscanf(line, "%19s %8191s", cinc, fileinc);
      if (nw < 1)
         continue;              // Not enough info in this line
      if (strcmp(cinc, "include") != 0) {
         // copy line in temporary file
         fprintf(ftmp, "%s\n", line);
      } else {

         // Drop quotes or double quotes, if any
         TString ln(line);
         ln.ReplaceAll("\"",1,"",0);
         ln.ReplaceAll("'",1,"",0);
         sscanf(ln.Data(), "%19s %8191s", cinc, fileinc);

         // support environment directories ...
         if (fileinc[0] == '$') {
            TString finc(fileinc);
            TString edir(fileinc);
            if (edir.Contains("/")) {
               edir.Remove(edir.Index("/"));
               edir.Remove(0,1);
               if (gSystem->Getenv(edir.Data())) {
                  finc.Remove(0,1);
                  finc.ReplaceAll(edir.Data(),gSystem->Getenv(edir.Data()));
                  fileinc[0] = '\0';
                  strncpy(fileinc,finc.Data(),kMAXPATHLEN);
                  fileinc[kMAXPATHLEN-1] = '\0';
               }
            }
         }

         // open (expand) file in temporary file ...
         if (fileinc[0] == '~') {
            // needs to expand
            int flen =
               strlen(fileinc) + strlen(gSystem->HomeDirectory()) + 10;
            char *ffull = new char[flen];
            snprintf(ffull, flen, "%s/%s", gSystem->HomeDirectory(), fileinc + 1);
            if (strlen(ffull) < kMAXPATHLEN - 1) strlcpy(fileinc, ffull,kMAXPATHLEN);
            delete [] ffull;
         }
         // Check if file exist and can be read ... ignore if not ...
         if (!gSystem->AccessPathName(fileinc, kReadPermission)) {
            FileExpand(fileinc, ftmp);
         } else {
            ::Warning("TAuthenticate::FileExpand",
                      "file specified by 'include' cannot be open or read (%s)",
                      fileinc);
         }
      }
   }
   fclose(fin);
}

////////////////////////////////////////////////////////////////////////////////
/// Determine default authentication details for method 'sec' and user 'usr'.
/// Checks .rootrc family files. Returned string must be deleted by the user.

char *TAuthenticate::GetDefaultDetails(int sec, int opt, const char *usr)
{
   char temp[kMAXPATHLEN] = { 0 };
   const char copt[2][5] = { "no", "yes" };

   if (gDebug > 2)
      ::Info("TAuthenticate::GetDefaultDetails",
             "enter ... %d ...pt:%d ... '%s'", sec, opt, usr);

   if (opt < 0 || opt > 1)
      opt = 1;

   // UsrPwd
   if (sec == TAuthenticate::kClear) {
      if (!usr[0] || !strncmp(usr,"*",1))
         usr = gEnv->GetValue("UsrPwd.Login", "");
      snprintf(temp, kMAXPATHLEN, "pt:%s ru:%s cp:%s us:%s",
               gEnv->GetValue("UsrPwd.LoginPrompt", copt[opt]),
               gEnv->GetValue("UsrPwd.ReUse", "1"),
               gEnv->GetValue("UsrPwd.Crypt", "1"), usr);
   }

   if (gDebug > 2)
      ::Info("TAuthenticate::GetDefaultDetails", "returning ... %s", temp);

   return StrDup(temp);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove THostAuth instance from the list

void TAuthenticate::RemoveHostAuth(THostAuth * ha, Option_t *opt)
{
   if (!strncasecmp(opt,"P",1))
      GetProofAuthInfo()->Remove(ha);
   else
      GetAuthInfo()->Remove(ha);
   // ... destroy it
   delete ha;
}

////////////////////////////////////////////////////////////////////////////////
/// Print info about the authentication sector.
/// If 'opt' contains 's' or 'S' prints information about established TSecContext,
/// else prints information about THostAuth (if 'opt' is 'p' or 'P', prints
/// Proof related information)

void TAuthenticate::Show(Option_t *opt)
{
   TString sopt(opt);

   if (sopt.Contains("s",TString::kIgnoreCase)) {

      // Print established security contexts
      TIter next(gROOT->GetListOfSecContexts());
      TSecContext *sc = 0;
      while ((sc = (TSecContext *)next()))
         sc->Print();

   } else {

      ::Info("::Print",
             " +--------------------------- BEGIN --------------------------------+");
      ::Info("::Print",
             " +                                                                  +");
      if (sopt.Contains("p",TString::kIgnoreCase)) {
         ::Info("::Print",
                " + List fgProofAuthInfo has %4d members                            +",
                GetProofAuthInfo()->GetSize());
         ::Info("::Print",
                " +                                                                  +");
         ::Info("::Print",
                " +------------------------------------------------------------------+");
         TIter next(GetProofAuthInfo());
         THostAuth *ai;
         while ((ai = (THostAuth *) next())) {
            ai->Print();
         }
      } else {
         ::Info("::Print",
                " + List fgAuthInfo has %4d members                                 +",
                GetAuthInfo()->GetSize());
         ::Info("::Print",
                " +                                                                  +");
         ::Info("::Print",
                " +------------------------------------------------------------------+");
         TIter next(GetAuthInfo());
         THostAuth *ai;
         while ((ai = (THostAuth *) next())) {
            ai->Print();
            ai->PrintEstablished();
         }
      }
      ::Info("::Print",
             " +---------------------------- END ---------------------------------+");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check if we have a valid established sec context in memory
/// Retrieves relevant info and negotiates with server.
/// options = "Opt,strlen(username),username.Data()"
/// message = kROOTD_USER, ...

Int_t TAuthenticate::AuthExists(TString username, Int_t method, const char *options,
                                Int_t *message, Int_t *rflag,
                                CheckSecCtx_t checksecctx)
{
   // Welcome message, if requested ...
   if (gDebug > 2)
      Info("AuthExists","%d: enter: msg: %d options: '%s'",
           method,*message, options);

   // Look for an existing security context matching this request
   Bool_t notHA = kFALSE;

   // First in the local list
   TIter next(fHostAuth->Established());
   TRootSecContext *secctx;
   while ((secctx = (TRootSecContext *)next())) {
      if (secctx->GetMethod() == method) {
         if (fRemote == secctx->GetHost()) {
            if (checksecctx &&
                (*checksecctx)(username,secctx) == 1)
               break;
         }
      }
   }

   // If nothing found, try the all list
   if (!secctx) {
      next = TIter(gROOT->GetListOfSecContexts());
      while ((secctx = (TRootSecContext *)next())) {
         if (secctx->GetMethod() == method) {
            if (fRemote == secctx->GetHost()) {
               if (checksecctx &&
                   (*checksecctx)(username,secctx) == 1) {
                  notHA = kTRUE;
                  break;
               }
            }
         }
      }
   }

   // If we have been given a valid sec context retrieve some info
   Int_t offset = -1;
   TString token;
   if (secctx) {
      offset = secctx->GetOffSet();
      token = secctx->GetToken();
      if (gDebug > 2)
         Info("AuthExists",
              "found valid TSecContext: offset: %d token: '%s'",
              offset, token.Data());
   }

   // Prepare string to be sent to the server
   TString sstr;
   sstr.Form("%d %d %s", fgProcessID, offset, options);

   // Send message
   if (fSocket->Send(sstr, *message) < 0)
      return -2;

   Int_t reuse = *rflag;
   if (reuse == 1 && offset > -1) {

      // Receive result of checking offset
      // But only for recent servers
      // NB: not backward compatible with dev version 4.00.02: switch
      // off 'reuse' for such servers to avoid hanging at this point.
      Int_t rproto = fSocket->GetRemoteProtocol();
      Bool_t oldsrv = ((fProtocol.BeginsWith("root") && rproto == 9) ||
                       (fProtocol.BeginsWith("proof") && rproto == 8));
      Int_t stat = 1, kind;
      if (!oldsrv) {
         if (fSocket->Recv(stat, kind) < 0)
            return -2;
         if (kind != kROOTD_AUTH)
            Warning("AuthExists","protocol error: expecting %d got %d"
                    " (value: %d)",kROOTD_AUTH,kind,stat);
      }

      if (stat > 0) {
         if (gDebug > 2)
            Info("AuthExists","offset OK");

         Int_t rsaKey = secctx->GetRSAKey();
         if (gDebug > 2)
            Info("AuthExists", "key type: %d", rsaKey);

         if (rsaKey > -1) {

            // Recent servers send a random tag in stat
            // It has to be signed too
            if (stat > 1) {
               // Create hex from tag
               char tag[9] = {0};
               snprintf(tag, 9, "%08x",stat);
               // Add to token
               token += tag;
            }

            // Send token encrypted
            if (SecureSend(fSocket, 1, rsaKey, token) == -1) {
               Warning("AuthExists", "problems secure-sending token %s",
                       "- may trigger problems in proofing Id ");
               return -2;
            }
         } else {
            // Send inverted
            for (int i = 0; i < token.Length(); i++) {
               char inv = ~token(i);
               token.Replace(i, 1, inv);
            }
            if (fSocket->Send(token, kMESS_STRING) < 0)
               return -2;
         }
      } else {
         if (gDebug > 0)
            Info("AuthExists","offset not OK - rerun authentication");
         // If the sec context was not valid, deactivate it ...
         if (secctx)
            secctx->DeActivate("");
      }
   }

   Int_t stat, kind;
   if (fSocket->Recv(stat, kind) < 0)
      return -2;
   if (gDebug > 3)
      Info("AuthExists","%d: after msg %d: kind= %d, stat= %d",
           method,*message, kind, stat);

   // Return flags
   *message = kind;
   *rflag = stat;

   if (kind == kROOTD_ERR) {
      TString server = "sockd";
      if (fSocket->GetServType() == TSocket::kROOTD)
         server = "rootd";
      if (fSocket->GetServType() == TSocket::kPROOFD)
         server = "proofd";
      if (stat == kErrConnectionRefused) {
         Error("AuthExists","%s@%s does not accept connections from %s@%s",
               server.Data(),fRemote.Data(),fUser.Data(),gSystem->HostName());
         return -2;
      } else if (stat == kErrNotAllowed) {
         if (gDebug > 0)
            Info("AuthExists",
                 "%s@%s does not accept %s authentication from %s@%s",
                 server.Data(),fRemote.Data(), fgAuthMeth[method].Data(),
                 fUser.Data(),gSystem->HostName());
      } else
         AuthError("AuthExists", stat);

      // If the sec context was not valid, deactivate it ...
      if (secctx)
         secctx->DeActivate("");
      return 0;
   }

   if (kind == kROOTD_AUTH && stat >= 1) {
      if (!secctx)
         secctx =
            fHostAuth->CreateSecContext(fUser,fRemote,method,-stat,fDetails,0);
      if (gDebug > 3) {
         if (stat == 1)
            Info("AuthExists", "valid authentication exists");
         if (stat == 2)
            Info("AuthExists", "valid authentication exists: offset changed");
         if (stat == 3)
            Info("AuthExists", "remote access authorized by /etc/hosts.equiv");
         if (stat == 4)
            Info("AuthExists", "no authentication required remotely");
      }

      if (stat == 2) {
         int newOffSet;
         // Receive new offset ...
         if (fSocket->Recv(newOffSet, kind) < 0)
            return -2;
         // ... and save it
         secctx->SetOffSet(newOffSet);
      }

      fSecContext = secctx;
      // Add it to local list for later use (if not already there)
      if (notHA)
         fHostAuth->Established()->Add(secctx);
      return 1;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize random machine using seed from /dev/urandom
/// (or current time if /dev/urandom not available).

void TAuthenticate::InitRandom()
{
   static Bool_t notinit = kTRUE;

   if (notinit) {
      const char *randdev = "/dev/urandom";
      Int_t fd;
      UInt_t seed;
      if ((fd = open(randdev, O_RDONLY)) != -1) {
         if (gDebug > 2)
            ::Info("InitRandom", "taking seed from %s", randdev);
         if (read(fd, &seed, sizeof(seed)) != sizeof(seed))
            ::Warning("InitRandom", "could not read seed from %s", randdev);
         close(fd);
      } else {
         if (gDebug > 2)
            ::Info("InitRandom", "%s not available: using time()", randdev);
         seed = time(0);   //better use times() + win32 equivalent
      }
      srand(seed);
      notinit = kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Generate a valid pair of private/public RSA keys to protect for
/// authentication token exchange

Int_t TAuthenticate::GenRSAKeys()
{
   if (gDebug > 2)
      Info("GenRSAKeys", "enter");

   if (fgRSAInit == 1) {
      if (gDebug > 2)
         Info("GenRSAKeys", "Keys prviously generated - return");
   }

   // This is for dynamic loads ...
   TString lib = "libRsa";

   // This is the local RSA implementation
   if (!TRSA_fun::RSA_genprim()) {
      char *p;
      if ((p = gSystem->DynamicPathName(lib, kTRUE))) {
         delete [] p;
         gSystem->Load(lib);
      }
   }

   // Init random machine
   TAuthenticate::InitRandom();

#ifdef R__SSL
   if (fgRSAKey == 1) {
      // Generate also the SSL key
      if (gDebug > 2)
         Info("GenRSAKeys","SSL: Generate Blowfish key");

      // Init SSL ...
      SSL_library_init();

      //  ... and its error strings
      SSL_load_error_strings();

      // Load Ciphers
      OpenSSL_add_all_ciphers();

      // Number of bits for key
      Int_t nbits = gEnv->GetValue("SSL.BFBits",256);

      // Minimum is 128
      nbits = (nbits >= 128) ? nbits : 128;

      // Max to limit size of buffers to 15912 (internal limitation)
      nbits = (nbits <= 15912) ? nbits : 15912;

      // Closer Number of chars
      Int_t klen = nbits / 8 ;

      // Init random engine
      char *rbuf = GetRandString(0,klen);
      RAND_seed(rbuf,strlen(rbuf));

      // This is what we export
      fgRSAPubExport[1].len = klen;
      fgRSAPubExport[1].keys = rbuf;
      if (gDebug > 2)
         Info("GenRSAKeys","SSL: BF key length: %d", fgRSAPubExport[1].len);

      // Now set the key locally in BF form
      BF_set_key(&fgBFKey, klen, (const unsigned char *)rbuf);
   }
#endif

   // Sometimes some bunch is not decrypted correctly
   // That's why we make retries to make sure that encryption/decryption
   // works as expected
   Bool_t notOk = 1;
   rsa_NUMBER p1, p2, rsa_n, rsa_e, rsa_d;
   Int_t l_n = 0, l_d = 0;
   char buf_n[rsa_STRLEN], buf_e[rsa_STRLEN], buf_d[rsa_STRLEN];
#if R__RSADE
   Int_t l_e;
   char buf[rsa_STRLEN];
#endif

   Int_t nAttempts = 0;
   Int_t thePrimeLen = kPRIMELENGTH;
   Int_t thePrimeExp = kPRIMEEXP;   // Prime probability = 1-0.5^thePrimeExp
   while (notOk && nAttempts < kMAXRSATRIES) {

      nAttempts++;
      if (gDebug > 2 && nAttempts > 1) {
         Info("GenRSAKeys", "retry no. %d",nAttempts);
         srand(auth_rand());
      }

      // Valid pair of primes
      p1 = TRSA_fun::RSA_genprim()(thePrimeLen, thePrimeExp);
      p2 = TRSA_fun::RSA_genprim()(thePrimeLen+1, thePrimeExp);

      // Retry if equal
      Int_t nPrimes = 0;
      while (TRSA_fun::RSA_cmp()(&p1, &p2) == 0 && nPrimes < kMAXRSATRIES) {
         nPrimes++;
         if (gDebug > 2)
            Info("GenRSAKeys", "equal primes: regenerate (%d times)",nPrimes);
         srand(auth_rand());
         p1 = TRSA_fun::RSA_genprim()(thePrimeLen, thePrimeExp);
         p2 = TRSA_fun::RSA_genprim()(thePrimeLen+1, thePrimeExp);
      }
#if R__RSADEB
      if (gDebug > 3) {
         TRSA_fun::RSA_num_sput()(&p1, buf, rsa_STRLEN);
         Info("GenRSAKeys", "local: p1: '%s' ", buf);
         TRSA_fun::RSA_num_sput()(&p2, buf, rsa_STRLEN);
         Info("GenRSAKeys", "local: p2: '%s' ", buf);
      }
#endif
      // Generate keys
      if (TRSA_fun::RSA_genrsa()(p1, p2, &rsa_n, &rsa_e, &rsa_d)) {
         if (gDebug > 2 && nAttempts > 1)
            Info("GenRSAKeys"," genrsa: unable to generate keys (%d)",
                 nAttempts);
         continue;
      }

      // Get equivalent strings and determine their lengths
      TRSA_fun::RSA_num_sput()(&rsa_n, buf_n, rsa_STRLEN);
      l_n = strlen(buf_n);
      TRSA_fun::RSA_num_sput()(&rsa_e, buf_e, rsa_STRLEN);
#if R__RSADEB
      l_e = strlen(buf_e);
#endif
      TRSA_fun::RSA_num_sput()(&rsa_d, buf_d, rsa_STRLEN);
      l_d = strlen(buf_d);

#if R__RSADEB
      if (gDebug > 3) {
         Info("GenRSAKeys", "local: n: '%s' length: %d", buf_n, l_n);
         Info("GenRSAKeys", "local: e: '%s' length: %d", buf_e, l_e);
         Info("GenRSAKeys", "local: d: '%s' length: %d", buf_d, l_d);
      }
#endif
      if (TRSA_fun::RSA_cmp()(&rsa_n, &rsa_e) <= 0)
         continue;
      if (TRSA_fun::RSA_cmp()(&rsa_n, &rsa_d) <= 0)
         continue;

      // Now we try the keys
      char test[2 * rsa_STRLEN] = "ThisIsTheStringTest01203456-+/";
      Int_t lTes = 31;
      char *tdum = GetRandString(0, lTes - 1);
      strlcpy(test, tdum, lTes+1);
      delete [] tdum;
      char buf[2 * rsa_STRLEN];
      if (gDebug > 3)
         Info("GenRSAKeys", "local: test string: '%s' ", test);

      // Private/Public
      strlcpy(buf, test, lTes+1);

      // Try encryption with private key
      int lout = TRSA_fun::RSA_encode()(buf, lTes, rsa_n, rsa_e);
      if (gDebug > 3)
         Info("GenRSAKeys",
              "local: length of crypted string: %d bytes", lout);

      // Try decryption with public key
      TRSA_fun::RSA_decode()(buf, lout, rsa_n, rsa_d);
      buf[lTes] = 0;
      if (gDebug > 3)
         Info("GenRSAKeys", "local: after private/public : '%s' ", buf);

      if (strncmp(test, buf, lTes))
         continue;

      // Public/Private
      strlcpy(buf, test, lTes+1);

      // Try encryption with public key
      lout = TRSA_fun::RSA_encode()(buf, lTes, rsa_n, rsa_d);
      if (gDebug > 3)
         Info("GenRSAKeys", "local: length of crypted string: %d bytes ",
              lout);

      // Try decryption with private key
      TRSA_fun::RSA_decode()(buf, lout, rsa_n, rsa_e);
      buf[lTes] = 0;
      if (gDebug > 3)
         Info("GenRSAKeys", "local: after public/private : '%s' ", buf);

      if (strncmp(test, buf, lTes))
         continue;

      notOk = 0;
   }

   // Save Private key
   TRSA_fun::RSA_assign()(&fgRSAPriKey.n, &rsa_n);
   TRSA_fun::RSA_assign()(&fgRSAPriKey.e, &rsa_e);

   // Save Public key
   TRSA_fun::RSA_assign()(&fgRSAPubKey.n, &rsa_n);
   TRSA_fun::RSA_assign()(&fgRSAPubKey.e, &rsa_d);

#if R__RSADEB
   if (gDebug > 2) {
      // Determine their lengths
      Info("GenRSAKeys", "local: generated keys are:");
      Info("GenRSAKeys", "local: n: '%s' length: %d", buf_n, l_n);
      Info("GenRSAKeys", "local: e: '%s' length: %d", buf_e, l_e);
      Info("GenRSAKeys", "local: d: '%s' length: %d", buf_d, l_d);
   }
#endif
   // Export form
   if (fgRSAPubExport[0].keys) {
      delete [] fgRSAPubExport[0].keys;
      fgRSAPubExport[0].len = 0;
   }
   fgRSAPubExport[0].len = l_n + l_d + 4;
   fgRSAPubExport[0].keys = new char[fgRSAPubExport[0].len];

   fgRSAPubExport[0].keys[0] = '#';
   memcpy(fgRSAPubExport[0].keys + 1, buf_n, l_n);
   fgRSAPubExport[0].keys[l_n + 1] = '#';
   memcpy(fgRSAPubExport[0].keys + l_n + 2, buf_d, l_d);
   fgRSAPubExport[0].keys[l_n + l_d + 2] = '#';
   fgRSAPubExport[0].keys[l_n + l_d + 3] = 0;
#if R__RSADEB
   if (gDebug > 2)
      Info("GenRSAKeys", "local: export pub: '%s'", fgRSAPubExport[0].keys);
#else
   if (gDebug > 2)
      Info("GenRSAKeys", "local: export pub length: %d bytes", fgRSAPubExport[0].len);
#endif

   // Set availability flag
   fgRSAInit = 1;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Allocates and fills a 0 terminated buffer of length len+1 with
/// len random characters.
/// Returns pointer to the buffer (to be deleted by the caller)
/// opt = 0      any non dangerous char
///       1      letters and numbers  (upper and lower case)
///       2      hex characters       (upper and lower case)

char *TAuthenticate::GetRandString(Int_t opt, Int_t len)
{
   unsigned int iimx[4][4] = {
      {0x0, 0xffffff08, 0xafffffff, 0x2ffffffe}, // opt = 0
      {0x0, 0x3ff0000, 0x7fffffe, 0x7fffffe},    // opt = 1
      {0x0, 0x3ff0000, 0x7e, 0x7e},              // opt = 2
      {0x0, 0x3ffc000, 0x7fffffe, 0x7fffffe}     // opt = 3
   };

   const char *cOpt[4] = { "Any", "LetNum", "Hex", "Crypt" };

   //  Default option 0
   if (opt < 0 || opt > 2) {
      opt = 0;
      if (gDebug > 2)
         Info("GetRandString", "unknown option: %d : assume 0", opt);
   }
   if (gDebug > 2)
      Info("GetRandString", "enter ... len: %d %s", len, cOpt[opt]);

   // Allocate buffer
   char *buf = new char[len + 1];

   // Init random machine (if needed)
   TAuthenticate::InitRandom();

   // randomize
   Int_t k = 0;
   Int_t i, j, l, m, frnd;
   while (k < len) {
      frnd = auth_rand();
      for (m = 7; m < 32; m += 7) {
         i = 0x7F & (frnd >> m);
         j = i / 32;
         l = i - j * 32;
         if ((iimx[opt][j] & (1 << l))) {
            buf[k] = i;
            k++;
         }
         if (k == len)
            break;
      }
   }

   // null terminated
   buf[len] = 0;
   if (gDebug > 3)
      Info("GetRandString", "got '%s' ", buf);

   return buf;
}

////////////////////////////////////////////////////////////////////////////////
/// Encode null terminated str using the session private key indicated by enc
/// and sends it over the network
/// Returns number of bytes sent, or -1 in case of error.
/// enc = 1 for private encoding, enc = 2 for public encoding

Int_t TAuthenticate::SecureSend(TSocket *sock, Int_t enc,
                                Int_t key, const char *str)
{
   char buftmp[kMAXSECBUF];
   char buflen[20];

   if (gDebug > 2)
      ::Info("TAuthenticate::SecureSend", "local: enter ... (enc: %d)", enc);

   Int_t slen = strlen(str) + 1;
   Int_t ttmp = 0;
   Int_t nsen = -1;

   if (key == 0) {
      strlcpy(buftmp, str, slen+1);

      if (enc == 1)
         ttmp = TRSA_fun::RSA_encode()(buftmp, slen, fgRSAPriKey.n,
                                       fgRSAPriKey.e);
      else if (enc == 2)
         ttmp = TRSA_fun::RSA_encode()(buftmp, slen, fgRSAPubKey.n,
                                       fgRSAPubKey.e);
      else
         return nsen;
   } else if (key == 1) {

#ifdef R__SSL
      ttmp = strlen(str);
      if ((ttmp % 8) > 0)           // It should be a multiple of 8!
         ttmp = ((ttmp + 8)/8) * 8;
      unsigned char iv[8];
      memset((void *)&iv[0],0,8);
      BF_cbc_encrypt((const unsigned char *)str, (unsigned char *)buftmp,
                     strlen(str), &fgBFKey, iv, BF_ENCRYPT);
#else
      if (gDebug > 0)
         ::Info("TAuthenticate::SecureSend","not compiled with SSL support:"
                " you should not have got here!");
#endif
   } else {
      if (gDebug > 0)
         ::Info("TAuthenticate::SecureSend","unknown key type (%d)",key);
      return nsen;
   }

   snprintf(buflen,20,"%d",ttmp);
   if (sock->Send(buflen, kROOTD_ENCRYPT) < 0)
      return -1;
   nsen = sock->SendRaw(buftmp, ttmp);
   if (gDebug > 3)
      ::Info("TAuthenticate::SecureSend",
             "local: sent %d bytes (expected: %d)", nsen,ttmp);

   return nsen;
}

////////////////////////////////////////////////////////////////////////////////
/// Receive str from sock and decode it using key indicated by key type
/// Return number of received bytes or -1 in case of error.
/// dec = 1 for private decoding, dec = 2 for public decoding

Int_t TAuthenticate::SecureRecv(TSocket *sock, Int_t dec, Int_t key, char **str)
{

   char buftmp[kMAXSECBUF];
   char buflen[20];

   Int_t nrec = -1;
   // We must get a pointer ...
   if (!str)
      return nrec;

   Int_t kind;
   if (sock->Recv(buflen, 20, kind) < 0)
      return -1;
   Int_t len = atoi(buflen);
   if (gDebug > 3)
      ::Info("TAuthenticate::SecureRecv", "got len '%s' %d (msg kind: %d)",
             buflen, len, kind);
   if (len == 0) {
      return len;
   }
   if (!strncmp(buflen, "-1", 2))
      return nrec;

   // Receive buffer
   if ((nrec = sock->RecvRaw(buftmp, len)) < 0)
      return nrec;
   if (key == 0) {
      if (dec == 1)
         TRSA_fun::RSA_decode()(buftmp, len, fgRSAPriKey.n, fgRSAPriKey.e);
      else if (dec == 2)
         TRSA_fun::RSA_decode()(buftmp, len, fgRSAPubKey.n, fgRSAPubKey.e);
      else
         return -1;

      // Prepare output
      const size_t strSize = strlen(buftmp) + 1;
      *str = new char[strSize];
      strlcpy(*str, buftmp, strSize);

   } else if (key == 1) {
#ifdef R__SSL
      unsigned char iv[8];
      memset((void *)&iv[0],0,8);
      *str = new char[nrec + 1];
      BF_cbc_encrypt((const unsigned char *)buftmp, (unsigned char *)(*str),
                     nrec, &fgBFKey, iv, BF_DECRYPT);
      (*str)[nrec] = '\0';
#else
      if (gDebug > 0)
         ::Info("TAuthenticate::SecureRecv","not compiled with SSL support:"
                " you should not have got here!");
#endif
   } else {
      if (gDebug > 0)
         ::Info("TAuthenticate::SecureRecv","unknown key type (%d)",key);
      return -1;
   }

   nrec= strlen(*str);

   return nrec;
}

////////////////////////////////////////////////////////////////////////////////
/// Store RSA public keys from export string rsaPubExport.

Int_t TAuthenticate::DecodeRSAPublic(const char *rsaPubExport, R__rsa_NUMBER &rsa_n,
                                     R__rsa_NUMBER &rsa_d, char **rsassl)
{
   if (!rsaPubExport)
      return -1;

   if (gDebug > 2)
      ::Info("TAuthenticate::DecodeRSAPublic",
             "enter: string length: %ld bytes", (Long_t)strlen(rsaPubExport));

   char str[kMAXPATHLEN] = { 0 };
   Int_t klen = strlen(rsaPubExport);
   if (klen > kMAXPATHLEN - 1) {
      ::Info("TAuthenticate::DecodeRSAPublic",
             "key too long (%d): truncate to %d",klen,kMAXPATHLEN);
      klen = kMAXPATHLEN - 1;
   }
   memcpy(str, rsaPubExport, klen);
   str[klen] ='\0';

   Int_t keytype = -1;

   if (klen > 0) {

      // Skip spaces at beginning, if any
      int k = 0;
      while (str[k] == 32) k++;

      if (str[k] == '#') {

         keytype = 0;

         // The format is #<hex_n>#<hex_d>#
         char *pd1 = strstr(str, "#");
         char *pd2 = pd1 ? strstr(pd1 + 1, "#") : (char *)0;
         char *pd3 = pd2 ? strstr(pd2 + 1, "#") : (char *)0;
         if (pd1 && pd2 && pd3) {
            // Get <hex_n> ...
            int l1 = (int) (pd2 - pd1 - 1);
            char *rsa_n_exp = new char[l1 + 1];
            strlcpy(rsa_n_exp, pd1 + 1, l1+1);
            if (gDebug > 2)
               ::Info("TAuthenticate::DecodeRSAPublic",
                      "got %ld bytes for rsa_n_exp", (Long_t)strlen(rsa_n_exp));
            // Now <hex_d>
            int l2 = (int) (pd3 - pd2 - 1);
            char *rsa_d_exp = new char[l2 + 1];
            strlcpy(rsa_d_exp, pd2 + 1, 13);
            if (gDebug > 2)
               ::Info("TAuthenticate::DecodeRSAPublic",
                      "got %ld bytes for rsa_d_exp", (Long_t)strlen(rsa_d_exp));

            TRSA_fun::RSA_num_sget()(&rsa_n, rsa_n_exp);
            TRSA_fun::RSA_num_sget()(&rsa_d, rsa_d_exp);

            delete[] rsa_n_exp;
            delete[] rsa_d_exp;

         } else
            ::Info("TAuthenticate::DecodeRSAPublic","bad format for input string");
#ifdef R__SSL
      } else {
         // try SSL
         keytype = 1;

         RSA *rsatmp;

         // Bio for exporting the pub key
         BIO *bpub = BIO_new(BIO_s_mem());

         // Write key from kbuf to BIO
         BIO_write(bpub,(void *)str,strlen(str));

         // Read pub key from BIO
         if (!(rsatmp = PEM_read_bio_RSAPublicKey(bpub, 0, 0, 0))) {
            if (gDebug > 0)
               ::Info("TAuthenticate::DecodeRSAPublic",
                        "unable to read pub key from bio");
         } else
            if (rsassl)
               *rsassl = (char *)rsatmp;
            else
               ::Info("TAuthenticate::DecodeRSAPublic",
                        "no space allocated for output variable");
         BIO_free(bpub);
      }
#else
      } else {
         if (rsassl) { }   // To avoid compiler complains
         if (gDebug > 0)
            ::Info("TAuthenticate::DecodeRSAPublic","not compiled with SSL support:"
                   " you should not have got here!");
      }
#endif
   }

   return keytype;
}

////////////////////////////////////////////////////////////////////////////////
/// Store RSA public keys from export string rsaPubExport.
/// Returns type of stored key, or -1 is not recognized

Int_t TAuthenticate::SetRSAPublic(const char *rsaPubExport, Int_t klen)
{
   if (gDebug > 2)
      ::Info("TAuthenticate::SetRSAPublic",
             "enter: string length %ld bytes", (Long_t)strlen(rsaPubExport));

   Int_t rsakey = -1;
   if (!rsaPubExport)
      return rsakey;

   if (klen > 0) {

      // Skip spaces at beginning, if any
      int k0 = 0;
      while (rsaPubExport[k0] == 32) k0++;
      int k2 = klen - 1;

      // Parse rsaPubExport
      // Type 0 is in the form
      //
      //   #< gt 10 exa chars >#< gt 10 exa chars >#
      //
      rsakey = 1;
      if (rsaPubExport[k0] == '#' && rsaPubExport[k2] == '#') {
         char *p0 = (char *)&rsaPubExport[k0];
         char *p2 = (char *)&rsaPubExport[k2];
         char *p1 = strchr(p0+1,'#');
         if (p1 > p0 && p1 < p2) {
            Int_t l01 = (Int_t)(p1-p0)-1;
            Int_t l12 = (Int_t)(p2-p1)-1;
            if (l01 >= kPRIMELENGTH*2 && l12 >= kPRIMELENGTH*2) {
               // Require exadecimal chars in between
               char *c = p0+1;
               while (c < p1 && ((*c < 58 && *c > 47) || (*c < 91 && *c > 64)))
                  c++;
               if (c == p1) {
                  c++;
                  while (c < p2 && ((*c < 58 && *c > 47) || (*c < 91 && *c > 64)))
                     c++;
                  if (c == p2)
                     rsakey = 0;
               }
            }
         }
      }
      if (gDebug > 3)
         ::Info("TAuthenticate::SetRSAPublic"," Key type: %d",rsakey);
      if (rsakey == 0) {

         // Decode input string
         R__rsa_NUMBER rsa_n, rsa_d;
         rsakey = TAuthenticate::DecodeRSAPublic(rsaPubExport,rsa_n,rsa_d);

         // Save Public key
         TRSA_fun::RSA_assign()(&fgRSAPubKey.n, &rsa_n);
         TRSA_fun::RSA_assign()(&fgRSAPubKey.e, &rsa_d);

      } else {
         rsakey = 1;
#ifdef R__SSL
         // Now set the key locally in BF form
         BF_set_key(&fgBFKey, klen, (const unsigned char *)rsaPubExport);
#else
         if (gDebug > 0)
            ::Info("TAuthenticate::SetRSAPublic",
                   "not compiled with SSL support:"
                   " you should not have got here!");
#endif
      }
   }

   return rsakey;
}

////////////////////////////////////////////////////////////////////////////////
/// Receives server RSA Public key
/// Sends local RSA public key encoded

Int_t TAuthenticate::SendRSAPublicKey(TSocket *socket, Int_t key)
{
   // Receive server public key
   char serverPubKey[kMAXSECBUF];
   int kind, nr = 0;
   if ((nr = socket->Recv(serverPubKey, kMAXSECBUF, kind)) < 0)
      return nr;
   if (gDebug > 3)
      ::Info("TAuthenticate::SendRSAPublicKey",
             "received key from server %ld bytes", (Long_t)strlen(serverPubKey));

   // Decode it
   R__rsa_NUMBER rsa_n, rsa_d;
#ifdef R__SSL
   char *tmprsa = 0;
   if (TAuthenticate::DecodeRSAPublic(serverPubKey,rsa_n,rsa_d,
                                      &tmprsa) != key) {
      if (tmprsa)
         RSA_free((RSA *)tmprsa);
      return -1;
   }
   RSA *RSASSLServer = (RSA *)tmprsa;
#else
   if (TAuthenticate::DecodeRSAPublic(serverPubKey,rsa_n,rsa_d) != key)
      return -1;
#endif

   // Send local public key, encodes
   char buftmp[kMAXSECBUF] = {0};
   char buflen[20] = {0};
   Int_t slen = fgRSAPubExport[key].len;
   Int_t ttmp = 0;
   if (key == 0) {
      strlcpy(buftmp,fgRSAPubExport[key].keys,slen+1);
      ttmp = TRSA_fun::RSA_encode()(buftmp, slen, rsa_n, rsa_d);
      snprintf(buflen, 20, "%d", ttmp);
   } else if (key == 1) {
#ifdef R__SSL
      Int_t lcmax = RSA_size(RSASSLServer) - 11;
      Int_t kk = 0;
      Int_t ke = 0;
      Int_t ns = slen;
      while (ns > 0) {
         Int_t lc = (ns > lcmax) ? lcmax : ns ;
         if ((ttmp = RSA_public_encrypt(lc,
                                        (unsigned char *)&fgRSAPubExport[key].keys[kk],
                                        (unsigned char *)&buftmp[ke],
                                        RSASSLServer,RSA_PKCS1_PADDING)) < 0) {
            char errstr[120];
            ERR_error_string(ERR_get_error(), errstr);
            ::Info("TAuthenticate::SendRSAPublicKey","SSL: error: '%s' ",errstr);
         }
         kk += lc;
         ke += ttmp;
         ns -= lc;
      }
      ttmp = ke;
      snprintf(buflen, 20, "%d", ttmp);
#else
      if (gDebug > 0)
         ::Info("TAuthenticate::SendRSAPublicKey","not compiled with SSL support:"
                " you should not have got here!");
      return -1;
#endif
   } else {
      if (gDebug > 0)
         ::Info("TAuthenticate::SendRSAPublicKey","unknown key type (%d)",key);
#ifdef R__SSL
      if (RSASSLServer)
         RSA_free(RSASSLServer);
#endif
      return -1;
   }

   // Send length first
   if ((nr = socket->Send(buflen, kROOTD_ENCRYPT)) < 0)
      return nr;
   // Send Key. second ...
   Int_t nsen = socket->SendRaw(buftmp, ttmp);
   if (gDebug > 3)
      ::Info("TAuthenticate::SendRSAPublicKey",
             "local: sent %d bytes (expected: %d)", nsen,ttmp);
#ifdef R__SSL
   if (RSASSLServer)
      RSA_free(RSASSLServer);
#endif
   return nsen;
}

////////////////////////////////////////////////////////////////////////////////
/// Read authentication directives from $ROOTAUTHRC, $HOME/.rootauthrc or
/// <Root_etc_dir>/system.rootauthrc and create related THostAuth objects.
/// Files are read only if they changed since last reading
/// If 'proofconf' is defined, check also file proofconf for directives

Int_t TAuthenticate::ReadRootAuthrc()
{
   // rootauthrc family
   char *authrc = 0;
   if (gSystem->Getenv("ROOTAUTHRC") != 0) {
      authrc = StrDup(gSystem->Getenv("ROOTAUTHRC"));
   } else {
      if (fgReadHomeAuthrc)
         authrc = gSystem->ConcatFileName(gSystem->HomeDirectory(), ".rootauthrc");
   }
   if (authrc && gDebug > 2)
      ::Info("TAuthenticate::ReadRootAuthrc", "Checking file: %s", authrc);
   if (!authrc || gSystem->AccessPathName(authrc, kReadPermission)) {
      if (authrc && gDebug > 1)
         ::Info("TAuthenticate::ReadRootAuthrc",
                "file %s cannot be read (errno: %d)", authrc, errno);
      delete [] authrc;
      authrc = gSystem->ConcatFileName(TROOT::GetEtcDir(), "system.rootauthrc");
      if (gDebug > 2)
         ::Info("TAuthenticate::ReadRootAuthrc", "Checking system file: %s", authrc);
      if (gSystem->AccessPathName(authrc, kReadPermission)) {
         if (gDebug > 1)
            ::Info("TAuthenticate::ReadRootAuthrc",
                   "file %s cannot be read (errno: %d)", authrc, errno);
         delete [] authrc;
         return 0;
      }
   }

   // Check if file has changed since last read
   TString tRootAuthrc = authrc;
   if (tRootAuthrc == fgRootAuthrc) {
      struct stat si;
      stat(tRootAuthrc, &si);
      if ((UInt_t)si.st_mtime < fgLastAuthrc.Convert()) {
         if (gDebug > 1)
            ::Info("TAuthenticate::ReadRootAuthrc",
                   "file %s already read", authrc);
         delete [] authrc;
         return 0;
      }
   }

   // Save filename in static variable
   fgRootAuthrc = tRootAuthrc;
   fgLastAuthrc = TDatime();

   // THostAuth lists
   TList *authinfo = TAuthenticate::GetAuthInfo();
   TList *proofauthinfo = TAuthenticate::GetProofAuthInfo();

   // Expand File into temporary file name and open it
   int expand = 1;
   TString filetmp = "rootauthrc";
   FILE *ftmp = gSystem->TempFileName(filetmp);
   if (gDebug > 2)
      ::Info("TAuthenticate::ReadRootAuthrc", "got tmp file: %s open at 0x%lx",
             filetmp.Data(), (Long_t)ftmp);
   if (ftmp == 0)
      expand = 0;  // Problems opening temporary file: ignore 'include's ...

   FILE *fd = 0;
   // If the temporary file is open, copy everything to the new file ...
   if (expand == 1) {
      TAuthenticate::FileExpand(authrc, ftmp);
      fd = ftmp;
      rewind(fd);
   } else {
      // Open file
      fd = fopen(authrc, "r");
      if (fd == 0) {
         if (gDebug > 2)
            ::Info("TAuthenticate::ReadRootAuthrc",
                   "file %s cannot be open (errno: %d)", authrc, errno);
         delete [] authrc;
         return 0;
      }
   }

   // Now scan file for meaningful directives
   TList tmpAuthInfo;
   char line[kMAXPATHLEN];
   Bool_t cont = kFALSE;
   TString proofserv;
   while (fgets(line, sizeof(line), fd) != 0) {

      // Skip comment lines
      if (line[0] == '#')
         continue;

      // Get rid of end of line '\n', if there ...
      if (line[strlen(line) - 1] == '\n')
         line[strlen(line) - 1] = '\0';

      // Skip empty lines
      if (!line[0])
         continue;

      // Now scan
      const size_t tmpSize = strlen(line) + 1;
      char *tmp = new char[tmpSize];
      if (!tmp) {
         ::Error("TAuthenticate::ReadRootAuthrc",
                 "could not allocate temporary buffer");
         fclose(fd);
         return 0;
      }
      strlcpy(tmp, line, tmpSize);
      char *nxt = strtok(tmp," ");

      if (!strcmp(nxt, "proofserv") || cont) {

         // Building the list of data servers for proof (analyzed at the end)
         char *ph = 0;
         if (cont)
            ph = nxt;
         else
            ph = strtok(0," ");
         while (ph) {
            if (*ph != 92) {
               proofserv += TString((const char *)ph);
               proofserv += TString(" ");
               cont = kFALSE;
            } else {
               cont = kTRUE;
            }
            ph = strtok(0," ");
         }

      } else {

         TString hostsrv = nxt;
         TString host   = hostsrv;
         TString server = "";
         if (hostsrv.Contains(":")) {
            server = hostsrv;
            host.Remove(host.Index(":"));
            server.Remove(0,server.Index(":")+1);
         }
         Int_t srvtyp = -1;
         if (server.Length()) {
            if (server == "0" || server.BeginsWith("sock"))
               srvtyp = TSocket::kSOCKD;
            else if (server == "1" || server.BeginsWith("root"))
               srvtyp = TSocket::kROOTD;
            else if (server == "2" || server.BeginsWith("proof"))
               srvtyp = TSocket::kPROOFD;
         }

         // Line with host info directives
         TString user = "*";

         nxt = strtok(0," ");
         if (!strncmp(nxt,"user",4)) {
            nxt = strtok(0," ");
            if (strncmp(nxt,"list",4) && strncmp(nxt,"method",6)) {
               user = TString(nxt);
               nxt = strtok(0," ");
            }
         }

         // Get related THostAuth, if exists in the tmp list,
         TIter next(&tmpAuthInfo);
         THostAuth *ha;
         while ((ha = (THostAuth *)next())) {
            if (host == ha->GetHost() && user == ha->GetUser() &&
                srvtyp == ha->GetServer())
               break;
         }
         if (!ha) {
            // Create a new one
            ha = new THostAuth(host,srvtyp,user);
            tmpAuthInfo.Add(ha);
         }

         if (!strncmp(nxt,"list",4)) {
            // list of methods for {host,usr}
            Int_t nm = 0, me[kMAXSEC] = {0};
            char *mth = strtok(0," ");
            while (mth) {
               Int_t met = -1;
               if (strlen(mth) > 1) {
                  // Method passed as string: translate it to number
                  met = GetAuthMethodIdx(mth);
                  if (met == -1 && gDebug > 2)
                     ::Info("TAuthenticate::ReadRootAuthrc",
                            "unrecognized method (%s): ", mth);
               } else {
                  met = atoi(mth);
               }
               if (met > -1 && met < kMAXSEC)
                  me[nm++] = met;
               mth = strtok(0," ");
            }
            if (nm)
               ha->ReOrder(nm,me);

         } else if (!strncmp(nxt,"method",6)) {

            // details for {host,usr,method}
            char *mth = strtok(0," ");
            Int_t met = -1;
            if (strlen(mth) > 1) {
               // Method passed as string: translate it to number
               met = GetAuthMethodIdx(mth);
               if (met == -1 && gDebug > 2)
                  ::Info("TAuthenticate::ReadRootAuthrc",
                         "unrecognized method (%s): ", mth);
            } else {
               met = atoi(mth);
            }
            if (met > -1 && met < kMAXSEC) {
               const char *det = 0;
               nxt = strtok(0," ");
               if (nxt) {
                  det = (const char *)strstr(line,nxt);
               }
               if (ha->HasMethod(met))
                  ha->SetDetails(met,det);
               else
                  ha->AddMethod(met,det);
            }
         }
      }
      if (tmp) delete [] tmp;
   }
   // Close file and remove it if temporary
   fclose(fd);
   if (expand == 1)
      gSystem->Unlink(filetmp);
   // Cleanup allocated memory
   delete [] authrc;

   // Update authinfo with new info found
   TAuthenticate::MergeHostAuthList(authinfo,&tmpAuthInfo);

   // Print those left, if requested ...
   if (gDebug > 2)
      TAuthenticate::Show();

   // Now create the list of THostAuth to be sent over to
   // the Master/Slaves, if requested ...
   TList tmpproofauthinfo;
   if (proofserv.Length() > 0) {
      char *tmps = new char[proofserv.Length()+1];
      strlcpy(tmps,proofserv.Data(),proofserv.Length()+1);
      char *nxt = strtok(tmps," ");
      while (nxt) {
         TString tmp((const char *)nxt);
         Int_t pdd = -1;
         // host
         TString host;
         if ((pdd = tmp.Index(":")) == -1) {
            host = tmp;
         } else {
            host = tmp;
            host.Resize(pdd);
            if (!host.Length())
               host = "*";
            tmp.Remove(0,pdd+1);
         }
         // user
         TString user;
         if ((pdd = tmp.Index(":")) == -1) {
            user = tmp;
         } else {
            user = tmp;
            user.Resize(pdd);
            if (!user.Length())
               user = "*";
            tmp.Remove(0,pdd+1);
         }
         // method(s)
         TString meth;
         Int_t nm = 0, me[kMAXSEC] = {0}, met = -1;
         while (tmp.Length() > 0) {
            meth = tmp;
            if ((pdd = tmp.Index(":")) > -1)
               meth.Resize(pdd);
            if (meth.Length() > 1) {
               // Method passed as string: translate it to number
               met = GetAuthMethodIdx(meth.Data());
               if (met == -1 && gDebug > 2)
                  ::Info("TAuthenticate::ReadRootAuthrc",
                         "unrecognized method (%s): ",meth.Data());
            } else if (meth.Length() == 1) {
               met = atoi(meth.Data());
               if (met > -1 && met < kMAXSEC)
                  me[nm++] = met;
            }
            if (pdd > -1)
               tmp.Remove(0,pdd+1);
            else
               tmp.Resize(0);
         }

         // Get related THostAuth, if exists, or create a new one
         THostAuth *ha = 0;
         THostAuth *hatmp = TAuthenticate::GetHostAuth(host,user);
         if (!hatmp) {
            ha = new THostAuth(host,user,nm,me,0);
         } else {
            // Create an empty THostAuth
            ha = new THostAuth(host,user);
            // Update with hatmp info
            ha->Update(hatmp);
            // ReOrder following new directives
            ha->ReOrder(nm,me);
         }
         // Add to the tmp list
         tmpproofauthinfo.Add(ha);
         // Go to next
         nxt = strtok(0," ");
      }
      delete [] tmps;
   }

   // Update proofauthinfo with new info found
   TAuthenticate::MergeHostAuthList(proofauthinfo,&tmpproofauthinfo,"P");
   // Print those, if requested ...
   if (gDebug > 2)
      TAuthenticate::Show("P");

   return authinfo->GetSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the authentication method can be attempted for the client.

Bool_t TAuthenticate::CheckProofAuth(Int_t cSec, TString &out)
{
   Bool_t rc = kFALSE;
   const char netrc[2][20] = { "/.netrc", "/.rootnetrc" };
   TString user;

   // Get user logon name
   UserGroup_t *pw = gSystem->GetUserInfo();
   if (pw) {
      user = TString(pw->fUser);
      delete pw;
   } else {
      ::Info("CheckProofAuth",
             "not properly logged on (getpwuid unable to find relevant info)!");
      out = "";
      return rc;
   }

   // UsrPwd
   if (cSec == (Int_t) TAuthenticate::kClear) {
      Int_t i = 0;
      for (; i < 2; i++) {
         TString infofile = TString(gSystem->HomeDirectory())+TString(netrc[i]);
         if (!gSystem->AccessPathName(infofile, kReadPermission))
            rc = kTRUE;
      }
      if (rc)
         out.Form("pt:0 ru:1 us:%s",user.Data());
   }

   if (gDebug > 3) {
      if (strlen(out) > 0)
         ::Info("CheckProofAuth",
                "meth: %d ... is available: details: %s", cSec, out.Data());
      else
         ::Info("CheckProofAuth",
                "meth: %d ... is NOT available", cSec);
   }

   // return
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Standard version of CheckSecCtx to be passed to TAuthenticate::AuthExists
/// Check if User is matches the one in Ctx
/// Returns: 1 if ok, 0 if not
/// Deactivates Ctx is not valid

Int_t StdCheckSecCtx(const char *user, TRootSecContext *ctx)
{
   Int_t rc = 0;

   if (ctx->IsActive()) {
      if (!strcmp(user,ctx->GetUser()) &&
          strncmp("AFS", ctx->GetID(), 3))
         rc = 1;
   }
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Tool for updating fgAuthInfo or fgProofAuthInfo
/// 'nin' contains list of last input information through (re)reading
/// of a rootauthrc-alike file. 'nin' info has priority.
/// 'std' is cleaned from inactive members.
/// 'nin' members used to update existing members in 'std' are
/// removed from 'nin', do that they do not leak
/// opt = "P" for proofauthinfo.

void TAuthenticate::MergeHostAuthList(TList *std, TList *nin, Option_t *opt)
{
   // Remove inactive from the 'std'
   TIter nxstd(std);
   THostAuth *ha;
   while ((ha = (THostAuth *) nxstd())) {
      if (!ha->IsActive()) {
         std->Remove(ha);
         SafeDelete(ha);
      }
   }

   // Merge 'nin' info in 'std'
   TIter nxnew(nin);
   THostAuth *hanew;
   while ((hanew = (THostAuth *)nxnew())) {
      if (hanew->NumMethods()) {
         TString hostsrv;
         hostsrv.Form("%s:%d",hanew->GetHost(),hanew->GetServer());
         THostAuth *hastd =
            TAuthenticate::HasHostAuth(hostsrv,hanew->GetUser(),opt);
         if (hastd) {
            // Update with new info
            hastd->Update(hanew);
            // Flag for removal
            hanew->DeActivate();
         } else {
            // Add new ThostAuth to std
            std->Add(hanew);
         }
      } else
         // Flag for removal empty objects
         hanew->DeActivate();
   }

   // Cleanup memory before quitting
   nxnew.Reset();
   while ((hanew = (THostAuth *)nxnew())) {
      if (!hanew->IsActive()) {
         nin->Remove(hanew);
         SafeDelete(hanew);
      }
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Tool for removing SecContext ctx from THostAuth listed in
/// fgAuthInfo or fgProofAuthInfo

void TAuthenticate::RemoveSecContext(TRootSecContext *ctx)
{
   THostAuth *ha = 0;

   // authinfo first
   TIter nxai(GetAuthInfo());
   while ((ha = (THostAuth *)nxai())) {
      TIter next(ha->Established());
      TRootSecContext *lctx = 0;
      while ((lctx = (TRootSecContext *) next())) {
         if (lctx == ctx) {
            ha->Established()->Remove(ctx);
            break;
         }
      }
   }

   // proofauthinfo second
   TIter nxpa(GetProofAuthInfo());
   while ((ha = (THostAuth *)nxpa())) {
      TIter next(ha->Established());
      TRootSecContext *lctx = 0;
      while ((lctx = (TRootSecContext *) next())) {
         if (lctx == ctx) {
            ha->Established()->Remove(ctx);
            break;
         }
      }
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Authentication related stuff setup in TProofServ.
/// This is the place where the buffer send by the client / master is
/// decoded. It contains also password information, if the case requires.
/// Return 0 on success, -1 on failure.

Int_t TAuthenticate::ProofAuthSetup()
{
   static Bool_t done = kFALSE;

   // Only once
   if (done)
      return 0;
   done = kTRUE;

   // Localise the buffer and decode it
   const char *p = gSystem->Getenv("ROOTPROOFAUTHSETUP");
   if (!p) {
      if (gDebug > 2)
         Info("ProofAuthSetup","Buffer not found: nothing to do");
      return 0;
   }
   TString mbuf = TBase64::Decode(p);

   // Create the message
   TMessage *mess = new TMessage((void*)mbuf.Data(), mbuf.Length()+sizeof(UInt_t));

   // Extract the information
   TString user = "";
   TString passwd = "";
   Bool_t  pwhash = kFALSE;
   Bool_t  srppwd = kFALSE;
   Int_t  rsakey = -1;
   *mess >> user >> passwd >> pwhash >> srppwd >> rsakey;

   // Set Globals for later use
   TAuthenticate::SetGlobalUser(user);
   TAuthenticate::SetGlobalPasswd(passwd);
   TAuthenticate::SetGlobalPwHash(pwhash);
   TAuthenticate::SetDefaultRSAKeyType(rsakey);
   const char *h = gSystem->Getenv("ROOTHOMEAUTHRC");
   if (h) {
      Bool_t rha = (Bool_t)(strtol(h, (char **)0, 10));
      TAuthenticate::SetReadHomeAuthrc(rha);
   }

   // Extract the list of THostAuth
   TList *pha = (TList *)mess->ReadObject(TList::Class());
   if (!pha) {
      if (gDebug > 0)
         Info("ProofAuthSetup","List of THostAuth not found");
      return 0;
   }

   Bool_t master = gROOT->IsProofServ();
   TIter next(pha);
   THostAuth *ha = 0;
   while ((ha = (THostAuth *)next())) {

      // Check if there is already one compatible
      Int_t kExact = 0;
      THostAuth *haex = 0;
      Bool_t fromProofAI = kFALSE;
      if (master) {
         // Look first in the proof list
         haex = TAuthenticate::GetHostAuth(ha->GetHost(),ha->GetUser(),"P",&kExact);
         // If nothing found, look also in the standard list
         if (!haex) {
            haex =
               TAuthenticate::GetHostAuth(ha->GetHost(),ha->GetUser(),"R",&kExact);
         } else
            fromProofAI = kTRUE;
      } else {
         // For slaves look first in the standard list only
         haex = TAuthenticate::GetHostAuth(ha->GetHost(),ha->GetUser(),"R",&kExact);
      }

      if (haex) {
         // If yes, action depends on whether it matches exactly or not
         if (kExact == 1) {
            // Update info in authinfo if Slave or in proofauthinfo
            // if Master and the entry was already in proofauthinfo
            if (!master || fromProofAI) {
               // update this existing one with the information found in
               // in the new one, if needed
               haex->Update(ha);
               // Delete temporary THostAuth
               SafeDelete(ha);
            } else
               // Master, entry not already in proofauthinfo,
               // Add it to the list
               TAuthenticate::GetProofAuthInfo()->Add(ha);
         } else {
            // update this new one with the information found in
            // in the existing one (if needed) and ...
            Int_t i = 0;
            for (; i < haex->NumMethods(); i++) {
               Int_t met = haex->GetMethod(i);
               if (!ha->HasMethod(met))
                  ha->AddMethod(met,haex->GetDetails(met));
            }
            if (master)
               // ... add the new one to the list
               TAuthenticate::GetProofAuthInfo()->Add(ha);
            else
               // We add this one to the standard list
               TAuthenticate::GetAuthInfo()->Add(ha);
         }
      } else {
         if (master)
            // We add this one to the list for forwarding
            TAuthenticate::GetProofAuthInfo()->Add(ha);
         else
            // We add this one to the standard list
            TAuthenticate::GetAuthInfo()->Add(ha);
      }
   }

   // We are done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Setup of authetication related stuff in PROOF run after a
/// successful authentication.
/// Return 0 on success, -1 on failure.

Int_t TAuthenticate::ProofAuthSetup(TSocket *sock, Bool_t /* client */)
{
   // Fill some useful info
   TSecContext *sc    = sock->GetSecContext();
   TString user       = sc->GetUser();
   Int_t remoteOffSet = sc->GetOffSet();

   // send user name to remote host
   // for UsrPwd method send also passwd, rsa encoded
   TMessage pubkey;
   TString passwd = "";
   Bool_t  pwhash = kFALSE;
   Bool_t  srppwd = kFALSE;

   Bool_t upwd = sc->IsA("UsrPwd");

   TPwdCtx *pwdctx = 0;
   if (remoteOffSet > -1 && upwd)
      pwdctx = (TPwdCtx *)(sc->GetContext());

   if (upwd && pwdctx) {
      passwd = pwdctx->GetPasswd();
      pwhash = pwdctx->IsPwHash();
   }

   Int_t keytyp = ((TRootSecContext *)sc)->GetRSAKey();

   // Prepare buffer
   TMessage mess;
   mess << user << passwd << pwhash << srppwd << keytyp;

   // Add THostAuth info
   mess.WriteObject(TAuthenticate::GetProofAuthInfo());

   // Get buffer as a base 64 string
   char *mbuf = mess.Buffer();
   Int_t mlen = mess.Length();
   TString messb64 = TBase64::Encode(mbuf, mlen);

   if (gDebug > 2)
      ::Info("ProofAuthSetup","sending %d bytes", messb64.Length());

   // Send it over
   if (remoteOffSet > -1) {
      if (TAuthenticate::SecureSend(sock, 1, keytyp, messb64.Data()) == -1) {
         ::Error("ProofAuthSetup","problems secure-sending message buffer");
         return -1;
      }
   } else {
      // There is no encryption key: send it plain
      char buflen[20];
      snprintf(buflen,20, "%d", messb64.Length());
      if (sock->Send(buflen, kMESS_ANY) < 0) {
         ::Error("ProofAuthSetup","plain: problems sending message length");
         return -1;
      }
      if (sock->SendRaw(messb64.Data(), messb64.Length()) < 0) {
         ::Error("ProofAuthSetup","problems sending message buffer");
         return -1;
      }
   }

   // We are done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning supported client protocol.

Int_t TAuthenticate::GetClientProtocol()
{
   return TSocket::GetClientProtocol();
}

//
// The code below is needed by TSlave and TProofServ for backward
// compatibility.
//

////////////////////////////////////////////////////////////////////////////////
/// Sends the list of the relevant THostAuth objects to the master or
/// to the active slaves, typically data servers external to the proof
/// cluster. The list is of THostAuth to be sent is specified by
/// TAuthenticate::fgProofAuthInfo after directives found in the
/// .rootauthrc family files ('proofserv' key)
/// Returns -1 if a problem sending THostAuth has occured, -2 in case
/// of problems closing the transmission.

static Int_t SendHostAuth(TSocket *s)
{
   Int_t retval = 0, ns = 0;

   if (!s) {
      Error("SendHostAuth","invalid input: socket undefined");
      return -1;
   }


   TIter next(TAuthenticate::GetProofAuthInfo());
   THostAuth *ha;
   while ((ha = (THostAuth *)next())) {
      TString buf;
      ha->AsString(buf);
      if((ns = s->Send(buf, kPROOF_HOSTAUTH)) < 1) {
         retval = -1;
         break;
      }
      if (gDebug > 2)
         Info("SendHostAuth","sent %d bytes (%s)",ns,buf.Data());
   }

   // End of transmission ...
   if ((ns = s->Send("END", kPROOF_HOSTAUTH)) < 1)
      retval = -2;
   if (gDebug > 2)
      Info("SendHostAuth","sent %d bytes for closing",ns);

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// Receive from client/master directives for authentications, create
/// related THostAuth and add them to the TAuthenticate::ProofAuthInfo
/// list. Opt = "M" or "m" if Master, "S" or "s" if Proof slave.
/// The 'proofconf' file is read only if Master

static Int_t RecvHostAuth(TSocket *s, Option_t *opt)
{
   if (!s) {
      Error("RecvHostAuth","invalid input: socket undefined");
      return -1;
   }

   // Check if Master
   Bool_t master = !strncasecmp(opt,"M",1) ? kTRUE : kFALSE;

   // First read directives from <rootauthrc>, <proofconf> and alike files
   TAuthenticate::ReadRootAuthrc();

   // Receive buffer
   Int_t kind;
   char buf[kMAXSECBUF];
   Int_t nr = s->Recv(buf, kMAXSECBUF, kind);
   if (nr < 0 || kind != kPROOF_HOSTAUTH) {
      Error("RecvHostAuth", "received: kind: %d (%d bytes)", kind, nr);
      return -1;
   }
   if (gDebug > 2)
      Info("RecvHostAuth","received %d bytes (%s)",nr,buf);

   while (strcmp(buf, "END")) {
      // Clean buffer
      Int_t nc = (nr >= kMAXSECBUF) ? kMAXSECBUF - 1 : nr ;
      buf[nc] = '\0';

      // Create THostAuth
      THostAuth *ha = new THostAuth((const char *)&buf);

      // Check if there is already one compatible
      Int_t kExact = 0;
      THostAuth *haex = 0;
      Bool_t fromProofAI = kFALSE;
      if (master) {
         // Look first in the proof list
         haex = TAuthenticate::GetHostAuth(ha->GetHost(),ha->GetUser(),"P",&kExact);
         // If nothing found, look also in the standard list
         if (!haex) {
            haex =
               TAuthenticate::GetHostAuth(ha->GetHost(),ha->GetUser(),"R",&kExact);
         } else
            fromProofAI = kTRUE;
      } else {
         // For slaves look first in the standard list only
         haex = TAuthenticate::GetHostAuth(ha->GetHost(),ha->GetUser(),"R",&kExact);
      }

      if (haex) {
         // If yes, action depends on whether it matches exactly or not
         if (kExact == 1) {
            // Update info in authinfo if Slave or in proofauthinfo
            // if master and the entry was already in proofauthinfo
            if (!master || fromProofAI) {
               // update this existing one with the information found in
               // in the new one, if needed
               haex->Update(ha);
               // Delete temporary THostAuth
               SafeDelete(ha);
            } else
               // master, entry not already in proofauthinfo,
               // Add it to the list
               TAuthenticate::GetProofAuthInfo()->Add(ha);
         } else {
            // update this new one with the information found in
            // in the existing one (if needed) and ...
            Int_t i = 0;
            for (; i < haex->NumMethods(); i++) {
               Int_t met = haex->GetMethod(i);
               if (!ha->HasMethod(met))
                  ha->AddMethod(met,haex->GetDetails(met));
            }
            if (master)
               // ... add the new one to the list
               TAuthenticate::GetProofAuthInfo()->Add(ha);
            else
               // We add this one to the standard list
               TAuthenticate::GetAuthInfo()->Add(ha);
         }
      } else {
         if (master)
            // We add this one to the list for forwarding
            TAuthenticate::GetProofAuthInfo()->Add(ha);
         else
            // We add this one to the standard list
            TAuthenticate::GetAuthInfo()->Add(ha);
      }


      // Get the next one
      nr = s->Recv(buf, kMAXSECBUF, kind);
      if (nr < 0 || kind != kPROOF_HOSTAUTH) {
         Info("RecvHostAuth","Error: received: kind: %d (%d bytes)", kind, nr);
         return -1;
      }
      if (gDebug > 2)
         Info("RecvHostAuth","received %d bytes (%s)",nr,buf);
   }

   return 0;
}

extern "C" {

////////////////////////////////////////////////////////////////////////////////
/// Setup of authetication in PROOF run after successful opening
/// of the socket. Provided for backward compatibility.
/// Return 0 on success, -1 on failure.

Int_t OldSlaveAuthSetup(TSocket *sock,
                        Bool_t /* master */, TString ord, TString conf)
{

   // Fill some useful info
   TSecContext *sc    = sock->GetSecContext();
   TString user       = sc->GetUser();
   Int_t proofdProto  = sock->GetRemoteProtocol();
   Int_t remoteOffSet = sc->GetOffSet();

   // send user name to remote host
   // for UsrPwd method send also passwd, rsa encoded
   TMessage pubkey;
   TString passwd = "";
   Bool_t  pwhash = kFALSE;
   Bool_t  srppwd = kFALSE;

   Bool_t upwd = sc->IsA("UsrPwd");

   TPwdCtx *pwdctx = 0;
   if (remoteOffSet > -1 && upwd)
      pwdctx = (TPwdCtx *)(sc->GetContext());

   if (upwd && pwdctx) {

      // Send offset to identify remotely the public part of RSA key
      if (sock->Send(remoteOffSet, kROOTD_RSAKEY) != 2*sizeof(Int_t)) {
         Error("OldAuthSetup", "failed to send offset in RSA key");
         return -1;
      }

      if (pwdctx) {
         passwd = pwdctx->GetPasswd();
         pwhash = pwdctx->IsPwHash();
      }

      Int_t keytyp = ((TRootSecContext *)sc)->GetRSAKey();
      if (TAuthenticate::SecureSend(sock, 1, keytyp, passwd.Data()) == -1) {
         if (remoteOffSet > -1)
            Warning("OldAuthSetup","problems secure-sending pass hash %s",
                    "- may result in failures");
         // If non RSA encoding available try passwd inversion
         if (upwd) {
            for (int i = 0; i < passwd.Length(); i++) {
               char inv = ~passwd(i);
               passwd.Replace(i, 1, inv);
            }
            TMessage mess;
            mess << passwd;
            if (sock->Send(mess) < 0) {
               Error("OldAuthSetup", "failed to send inverted password");
               return -1;
            }
         }
      }

   } else {

      // Send notification of no offset to be sent ...
      if (sock->Send(-2, kROOTD_RSAKEY) != 2*sizeof(Int_t)) {
         Error("OldAuthSetup", "failed to send no offset notification in RSA key");
         return -1;
      }
   }

   // Send ordinal (and config) info to slave (or master)
   TMessage mess;
   mess << user << pwhash << srppwd << ord << conf;

   if (sock->Send(mess) < 0) {
      Error("OldAuthSetup", "failed to send ordinal and config info");
      return -1;
   }

   if (proofdProto > 6) {
      // Now we send authentication details to access, e.g., data servers
      // not in the proof cluster and to be propagated to slaves.
      // This is triggered by the 'proofserv <dserv1> <dserv2> ...'
      // line in .rootauthrc
      if (SendHostAuth(sock) < 0) {
         Error("OldAuthSetup", "failed to send HostAuth info");
         return -1;
      }
   }

   // We are done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Authentication related setup in TProofServ run after successful
/// startup. Provided for backward compatibility.
/// Return 0 on success, -1 on failure.

Int_t OldProofServAuthSetup(TSocket *sock, Bool_t master, Int_t protocol,
                            TString &user, TString &ord, TString &conf)
{
   // First receive, decode and store the public part of RSA key
   Int_t retval, kind;
   if (sock->Recv(retval, kind) != 2*sizeof(Int_t)) {
      //other side has closed connection
      Info("OldProofServAuthSetup",
           "socket has been closed due to protocol mismatch - Exiting");
      return -1;
   }

   Int_t rsakey = 0;
   TString passwd;
   if (kind == kROOTD_RSAKEY) {

      if (retval > -1) {
         if (gSystem->Getenv("ROOTKEYFILE")) {

            TString keyfile = gSystem->Getenv("ROOTKEYFILE");
            keyfile += retval;

            FILE *fKey = 0;
            char pubkey[kMAXPATHLEN] = { 0 };
            if (!gSystem->AccessPathName(keyfile.Data(), kReadPermission)) {
               if ((fKey = fopen(keyfile.Data(), "r"))) {
                  Int_t klen = fread((void *)pubkey,1,sizeof(pubkey),fKey);
                  if (klen <= 0) {
                     Error("OldProofServAuthSetup",
                           "failed to read public key from '%s'", keyfile.Data());
                     fclose(fKey);
                     return -1;
                  }
                  pubkey[klen] = 0;
                  // Set RSA key
                  rsakey = TAuthenticate::SetRSAPublic(pubkey,klen);
                  fclose(fKey);
               } else {
                  Error("OldProofServAuthSetup", "failed to open '%s'", keyfile.Data());
                  return -1;
               }
            }
         }

         // Receive passwd
         char *pwd = 0;
         if (TAuthenticate::SecureRecv(sock, 2, rsakey, &pwd) < 0) {
            Error("OldProofServAuthSetup", "failed to receive password");
            return -1;
         }
         passwd = pwd;
         delete[] pwd;

      } else if (retval == -1) {

         // Receive inverted passwd
         TMessage *mess;
         if ((sock->Recv(mess) <= 0) || !mess) {
            Error("OldProofServAuthSetup", "failed to receive inverted password");
            return -1;
         }
         (*mess) >> passwd;
         delete mess;

         for (Int_t i = 0; i < passwd.Length(); i++) {
            char inv = ~passwd(i);
            passwd.Replace(i, 1, inv);
         }

      }
   }

   // Receive final information
   TMessage *mess;
   if ((sock->Recv(mess) <= 0) || !mess) {
      Error("OldProofServAuthSetup", "failed to receive ordinal and config info");
      return -1;
   }

   // Decode it
   Bool_t pwhash, srppwd;
   if (master) {
      if (protocol < 4) {
         (*mess) >> user >> pwhash >> srppwd >> conf;
         ord = "0";
      } else {
         (*mess) >> user >> pwhash >> srppwd >> ord >> conf;
      }
   } else {
      if (protocol < 4) {
         Int_t iord;
         (*mess) >> user >> pwhash >> srppwd >> iord;
         ord = "0.";
         ord += iord;
      } else {
         (*mess) >> user >> pwhash >> srppwd >> ord >> conf;
      }
   }
   delete mess;

   // Set Globals for later use
   TAuthenticate::SetGlobalUser(user);
   TAuthenticate::SetGlobalPasswd(passwd);
   TAuthenticate::SetGlobalPwHash(pwhash);
   TAuthenticate::SetDefaultRSAKeyType(rsakey);
   const char *h = gSystem->Getenv("ROOTHOMEAUTHRC");
   if (h) {
      Bool_t rha = (Bool_t)(strtol(h, (char **)0, 10));
      TAuthenticate::SetReadHomeAuthrc(rha);
   }

   // Read user or system authentication directives and
   // receive auth info transmitted from the client
   Int_t harc = master ? RecvHostAuth(sock, "M") : RecvHostAuth(sock, "S");

   if (harc < 0) {
      Error("OldProofServAuthSetup", "failed to receive HostAuth info");
      return -1;
   }

   // We are done
   return 0;
}

}  // extern "C"
