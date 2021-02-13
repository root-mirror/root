// @(#)root/base:$Id$
// Author: Fons Rademakers   17/01/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TUrl
\ingroup Base

This class represents a WWW compatible URL.
It provides member functions to return the different parts of
an URL. The supported url format is:
~~~ {.cpp}
  [proto://][user[:passwd]@]host[:port]/file.ext[#anchor][?options]
~~~
*/

#include <stdlib.h>
#include "TUrl.h"
#include "THashList.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TMap.h"
#include "TROOT.h"

#include <atomic>

TObjArray *TUrl::fgSpecialProtocols = nullptr;
THashList *TUrl::fgHostFQDNs = nullptr;

#ifdef R__COMPLETE_MEM_TERMINATION
namespace {
   class TUrlCleanup {
      TObjArray **fSpecialProtocols;
      THashList **fHostFQDNs;
   public:
      TUrlCleanup(TObjArray **protocols, THashList **hosts) : fSpecialProtocols(protocols),fHostFQDNs(hosts) {}
      ~TUrlCleanup() {
         if (*fSpecialProtocols) (*fSpecialProtocols)->Delete();
         delete *fSpecialProtocols;
         *fSpecialProtocols = 0;
         if (*fHostFQDNs) (*fHostFQDNs)->Delete();
         delete *fHostFQDNs;
         *fHostFQDNs = 0;
      }
   };
}
#endif

ClassImp(TUrl);

////////////////////////////////////////////////////////////////////////////////
/// Parse url character string and split in its different subcomponents.
/// Use IsValid() to check if URL is legal.
/// ~~~ {.cpp}
/// url: [proto://][user[:passwd]@]host[:port]/file.ext[?options][#anchor]
/// ~~~
/// Known protocols: http, root, proof, ftp, news and any special protocols
/// defined in the rootrc Url.Special key.
/// The default protocol is "http", unless defaultIsFile is true in which
/// case the url is assumed to be of type "file".
/// If a passwd contains a @ it must be escaped by a \\, e.g.
/// "pip@" becomes "pip\\@".
///
/// Default ports: http=80, root=1094, proof=1093, ftp=20, news=119.
/// Port #1093 has been assigned by IANA (www.iana.org) to proofd.
/// Port #1094 has been assigned by IANA (www.iana.org) to rootd.

TUrl::TUrl(const char *url, Bool_t defaultIsFile)
{
   SetUrl(url, defaultIsFile);

#ifdef R__COMPLETE_MEM_TERMINATION
   static TUrlCleanup cleanup(&fgSpecialProtocols,&fgHostFQDNs);
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup.

TUrl::~TUrl()
{
   delete fOptionsMap;
}

////////////////////////////////////////////////////////////////////////////////
/// Parse url character string and split in its different subcomponents.
/// Use IsValid() to check if URL is legal.
///~~~ {.cpp}
/// url: [proto://][user[:passwd]@]host[:port]/file.ext[?options][#anchor]
///~~~
/// Known protocols: http, root, proof, ftp, news and any special protocols
/// defined in the rootrc Url.Special key.
/// The default protocol is "http", unless defaultIsFile is true in which
/// case the url is assumed to be of type "file".
/// If a passwd contains a @ it must be escaped by a \\, e.g.
/// "pip@" becomes "pip\\@".
///
/// Default ports: http=80, root=1094, proof=1093, ftp=20, news=119.
/// Port #1093 has been assigned by IANA (www.iana.org) to proofd.
/// Port #1094 has been assigned by IANA (www.iana.org) to rootd.

void TUrl::SetUrl(const char *url, Bool_t defaultIsFile)
{
   delete fOptionsMap;
   fOptionsMap = nullptr;

   if (!url || !url[0]) {
      fPort = -1;
      return;
   }

   // Set defaults
   fUrl        = "";
   fProtocol   = "http";
   fUser       = "";
   fPasswd     = "";
   fHost       = "";
   fPort       = 80;
   fFile       = "";
   fAnchor     = "";
   fOptions    = "";
   fFileOA     = "";
   fHostFQ     = "";

   // if url starts with a / consider it as a file url
   if (url[0] == '/')
      defaultIsFile = kTRUE;

   // Find protocol
   char *s, sav;

   char *u, *u0 = Strip(url);
tryfile:
   u = u0;

   // Handle special protocol cases: "file:", etc.
   for (int i = 0; i < GetSpecialProtocols()->GetEntriesFast(); i++) {
      TObjString *os = (TObjString*) GetSpecialProtocols()->UncheckedAt(i);
      TString s1 = os->GetString();
      int l = s1.Length();
      Bool_t stripoff = kFALSE;
      if (s1.EndsWith("/-")) {
         stripoff = kTRUE;
         s1 = s1.Strip(TString::kTrailing, '-');
         l--;
      }
      if (!strncmp(u, s1, l)) {
         if (s1(0) == '/' && s1(l-1) == '/') {
            // case with file namespace like: /alien/user/file.root
            fProtocol = s1(1, l-2);
            if (stripoff)
               l--;    // strip off namespace prefix from file name
            else
               l = 0;  // leave namespace prefix as part of file name
         } else {
            // case with protocol, like: file:/data/file.root
            fProtocol = s1(0, l-1);
         }
         if (!strncmp(u+l, "//", 2))
            u += l+2;
         else
            u += l;
         fPort = 0;

         FindFile(u, kFALSE);

         delete [] u0;
         return;
      }
   }

   u = u0;

   char *x, *t, *s2;
   // allow x:/path as Windows filename
   if ((s = strstr(u, ":/")) && u+1 != s) {
      if (*(s+2) != '/') {
         Error("TUrl", "%s malformed, URL must contain \"://\"", u0);
         fPort = -1;
         goto cleanup;
      }
      sav = *s;
      *s = 0;
      SetProtocol(u, kTRUE);
      *s = sav;
      s += 3;
      // allow url of form: "proto://"
   } else {
      if (defaultIsFile) {
         char *newu = new char [strlen("file:") + strlen(u0) + 1];
         sprintf(newu, "file:%s", u0);
         delete [] u0;
         u0 = newu;
         goto tryfile;
      }
      s = u;
   }

   // Find user and passwd
   u = s;
   t = s;
again:
   if ((s = strchr(t, '@')) && (
       ((x = strchr(t, '/')) && s < x) ||
       ((x = strchr(t, '?')) && s < x) ||
       ((x = strchr(t, '#')) && s < x) ||
       (!strchr(t, '/'))
      )) {
      if (*(s-1) == '\\') {
         t = s+1;
         goto again;
      }
      sav = *s;
      *s = 0;
      if ((s2 = strchr(u, ':'))) {
         *s2 = 0;
         fUser = u;
         *s2 = ':';
         s2++;
         if (*s2) {
            fPasswd = s2;
            fPasswd.ReplaceAll("\\@", "@");
         }
      } else
         fUser = u;
      *s = sav;
      s++;
   } else
      s = u;

   // Find host
   u = s;
   if ((s = strchr(u, ':')) || (s = strchr(u, '/')) || (s = strchr(u, '?')) || (s = strchr(u, '#'))) {
      if ((strchr (u, ':') > strchr(u, '/')) && (strchr (u, '/')))
         s = strchr(u, '/');
      sav = *s;
      *s = 0;
      fHost = u;
      *s = sav;
      if (sav == ':') {
         s++;
         // Get port #
         if (!*s) {
            fPort = -1;
            goto cleanup;
         }
         u = s;
         if ((s = strchr(u, '/')) || (s = strchr(u, '?')) || (s = strchr(u, '#'))) {
            sav = *s;
            *s = 0;
            fPort = atoi(u);
            *s = sav;
         } else {
            fPort = atoi(u);
            goto cleanup;
         }
      }
   } else {
      fHost = u;
      goto cleanup;
   }

   if (!*s) goto cleanup;

   // Find file
   u = s;
   if (*u == '/' && fHost.Length())
      u++;

   FindFile(u);

cleanup:
   delete [] u0;
}

////////////////////////////////////////////////////////////////////////////////
/// Find file and optionally anchor and options.

void TUrl::FindFile(char *u, Bool_t stripDoubleSlash)
{
   char *s, sav;

   // Locate anchor and options, if any
   char *opt = strchr(u, '?');
   char *anc = strchr(u, '#');

   // URL invalid if anchor is coming before the options
   if (opt && anc && opt > anc) {
      fPort = -1;
      return;
   }

   if ((s = opt) || (s = anc)) {
      sav = *s;
      *s = 0;
      fFile = u;
      if (stripDoubleSlash)
         fFile.ReplaceAll("//", "/");
      *s = sav;
      s++;
      if (sav == '?') {
         // Get options
         if (!*s) {
            // options string is empty
            return;
         }
         u = s;
         if ((s = strchr(u, '#'))) {
            sav = *s;
            *s = 0;
            fOptions = u;
            *s = sav;
            s++;
         } else {
            fOptions = u;
            return;
         }
      }
      if (!*s) {
         // anchor string is empty
         return;
      }
   } else {
      fFile = u;
      if (stripDoubleSlash)
         fFile.ReplaceAll("//", "/");
      return;
   }

   // Set anchor
   fAnchor = s;
}

////////////////////////////////////////////////////////////////////////////////
/// TUrl copy ctor.

TUrl::TUrl(const TUrl &url) : TObject(url)
{
   fUrl        = url.fUrl;
   fProtocol   = url.fProtocol;
   fUser       = url.fUser;
   fPasswd     = url.fPasswd;
   fHost       = url.fHost;
   fFile       = url.fFile;
   fAnchor     = url.fAnchor;
   fOptions    = url.fOptions;
   fPort       = url.fPort;
   fFileOA     = url.fFileOA;
   fHostFQ     = url.fHostFQ;
   fOptionsMap = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// TUrl assignment operator.

TUrl &TUrl::operator=(const TUrl &rhs)
{
   if (this != &rhs) {
      TObject::operator=(rhs);
      fUrl        = rhs.fUrl;
      fProtocol   = rhs.fProtocol;
      fUser       = rhs.fUser;
      fPasswd     = rhs.fPasswd;
      fHost       = rhs.fHost;
      fFile       = rhs.fFile;
      fAnchor     = rhs.fAnchor;
      fOptions    = rhs.fOptions;
      fPort       = rhs.fPort;
      fFileOA     = rhs.fFileOA;
      fHostFQ     = rhs.fHostFQ;
      delete fOptionsMap;
      fOptionsMap = nullptr;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Return full URL. If withDflt is kTRUE, explicitly add the port even
/// if it matches the default value for the URL protocol.

const char *TUrl::GetUrl(Bool_t withDeflt) const
{
   if (((TestBit(kUrlWithDefaultPort) && !withDeflt) ||
       (!TestBit(kUrlWithDefaultPort) && withDeflt)) &&
       TestBit(kUrlHasDefaultPort))
      fUrl = "";

   if (IsValid() && fUrl == "") {
      // Handle special protocol cases: file:, etc.
      for (int i = 0; i < GetSpecialProtocols()->GetEntriesFast(); i++) {
         TObjString *os = (TObjString*) GetSpecialProtocols()->UncheckedAt(i);
         TString &s = os->String();
         int l = s.Length();
         if (fProtocol == s(0, l-1)) {
            if (fFile[0] == '/')
               fUrl = fProtocol + "://" + fFile;
            else
               fUrl = fProtocol + ":" + fFile;
            if (fOptions != "") {
               fUrl += "?";
               fUrl += fOptions;
            }
            if (fAnchor != "") {
               fUrl += "#";
               fUrl += fAnchor;
            }
            return fUrl;
         }
      }

      Bool_t deflt = kFALSE;
      if ((!fProtocol.CompareTo("http")  && fPort == 80)   ||
          (fProtocol.BeginsWith("proof") && fPort == 1093) ||
          (fProtocol.BeginsWith("root")  && fPort == 1094) ||
          (!fProtocol.CompareTo("ftp")   && fPort == 20)   ||
          (!fProtocol.CompareTo("news")  && fPort == 119)  ||
          (!fProtocol.CompareTo("https") && fPort == 443)  ||
          fPort == 0) {
         deflt = kTRUE;
         ((TUrl *)this)->SetBit(kUrlHasDefaultPort);
      }

      fUrl = fProtocol + "://";
      if (fUser != "") {
         fUrl += fUser;
         if (fPasswd != "") {
            fUrl += ":";
            TString passwd = fPasswd;
            passwd.ReplaceAll("@", "\\@");
            fUrl += passwd;
         }
         fUrl += "@";
      }
      if (withDeflt)
         ((TUrl*)this)->SetBit(kUrlWithDefaultPort);
      else
         ((TUrl*)this)->ResetBit(kUrlWithDefaultPort);

      if (!deflt || withDeflt) {
         char p[10];
         sprintf(p, "%d", fPort);
         fUrl = fUrl + fHost + ":" + p + "/" + fFile;
      } else
         fUrl = fUrl + fHost + "/" + fFile;
      if (fOptions != "") {
         fUrl += "?";
         fUrl += fOptions;
      }
      if (fAnchor != "") {
         fUrl += "#";
         fUrl += fAnchor;
      }
   }

   fUrl.ReplaceAll("////", "///");
   return fUrl;
}

////////////////////////////////////////////////////////////////////////////////
/// Return fully qualified domain name of url host. If host cannot be
/// resolved or not valid return the host name as originally specified.

const char *TUrl::GetHostFQDN() const
{
   if (fHostFQ == "") {
      // Check if we already resolved it
      TNamed *fqdn = fgHostFQDNs ? (TNamed *) fgHostFQDNs->FindObject(fHost) : 0;
      if (!fqdn) {
         TInetAddress adr(gSystem->GetHostByName(fHost));
         if (adr.IsValid()) {
            fHostFQ = adr.GetHostName();
         } else
            fHostFQ = "-";
         R__LOCKGUARD(gROOTMutex);
         if (!fgHostFQDNs) {
            fgHostFQDNs = new THashList;
            fgHostFQDNs->SetOwner();
         }
         if (fgHostFQDNs && !fgHostFQDNs->FindObject(fHost))
            fgHostFQDNs->Add(new TNamed(fHost,fHostFQ));
      } else {
         fHostFQ = fqdn->GetTitle();
      }
   }
   if (fHostFQ == "-")
      return fHost;
   return fHostFQ;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the file and its options (the string specified behind the ?).
/// Convenience function useful when the option is used to pass
/// authentication/access information for the specified file.

const char *TUrl::GetFileAndOptions() const
{
   if (fFileOA == "") {
      fFileOA = fFile;
      if (fOptions != "") {
         fFileOA += "?";
         fFileOA += fOptions;
      }
      if (fAnchor != "") {
         fFileOA += "#";
         fFileOA += fAnchor;
      }
   }
   return fFileOA;
}

////////////////////////////////////////////////////////////////////////////////
/// Set protocol and, optionally, change the port accordingly.

void TUrl::SetProtocol(const char *proto, Bool_t setDefaultPort)
{
   fProtocol = proto;
   if (setDefaultPort) {
      if (!fProtocol.CompareTo("http"))
         fPort = 80;
      else if (!fProtocol.CompareTo("https"))
         fPort = 443;
      else if (fProtocol.BeginsWith("proof"))  // can also be proofs or proofk
         fPort = 1093;
      else if (fProtocol.BeginsWith("root"))   // can also be roots or rootk
         fPort = 1094;
      else if (!fProtocol.CompareTo("ftp"))
         fPort = 20;
      else if (!fProtocol.CompareTo("news"))
         fPort = 119;
      else {
         // generic protocol (no default port)
         fPort = 0;
      }
   }
   fUrl = "";
}

////////////////////////////////////////////////////////////////////////////////
/// Compare two urls as strings.

Int_t TUrl::Compare(const TObject *obj) const
{
   if (this == obj) return 0;
   if (TUrl::Class() != obj->IsA()) return -1;
   return TString(GetUrl()).CompareTo(((TUrl*)obj)->GetUrl(), TString::kExact);
}

////////////////////////////////////////////////////////////////////////////////
/// Print URL on stdout.

void TUrl::Print(Option_t *) const
{
   if (fPort == -1)
      Printf("Illegal URL");

   Printf("%s", GetUrl());
}

////////////////////////////////////////////////////////////////////////////////
/// Read the list of special protocols from the rootrc files.
/// These protocols will be parsed in a protocol and a file part,
/// no host or other info will be determined. This is typically
/// used for legacy file descriptions like: file:/path/file.root.

TObjArray *TUrl::GetSpecialProtocols()
{
   static std::atomic_bool usedEnv = ATOMIC_VAR_INIT(false);

   if (!gEnv) {
      R__LOCKGUARD(gROOTMutex);
      if (!fgSpecialProtocols)
         fgSpecialProtocols = new TObjArray;
      if (fgSpecialProtocols->GetEntriesFast() == 0)
         fgSpecialProtocols->Add(new TObjString("file:"));
      return fgSpecialProtocols;
   }

   if (usedEnv)
      return fgSpecialProtocols;

   R__LOCKGUARD(gROOTMutex);

   // Some other thread might have set it up in the meantime.
   if (usedEnv)
      return fgSpecialProtocols;

   if (fgSpecialProtocols)
      fgSpecialProtocols->Delete();

   if (!fgSpecialProtocols)
      fgSpecialProtocols = new TObjArray;

   const char *protos = gEnv->GetValue("Url.Special", "file: hpss: dcache: dcap:");

   if (protos) {
      Int_t cnt = 0;
      char *p = StrDup(protos);
      while (1) {
         TObjString *proto = new TObjString(strtok(!cnt ? p : 0, " "));
         if (proto->String().IsNull()) {
            delete proto;
            break;
         }
         fgSpecialProtocols->Add(proto);
         cnt++;
      }
      delete [] p;
   }
   usedEnv = true;
   return fgSpecialProtocols;
}


////////////////////////////////////////////////////////////////////////////////
/// Parse URL options into a key/value map.

void TUrl::ParseOptions() const
{
   if (fOptionsMap) return;

   TString urloptions = GetOptions();
   if (urloptions.IsNull())
      return;

   TObjArray *objOptions = urloptions.Tokenize("&");
   for (Int_t n = 0; n < objOptions->GetEntries(); n++) {
      TString loption = ((TObjString *) objOptions->At(n))->GetName();
      TObjArray *objTags = loption.Tokenize("=");
      if (!fOptionsMap) {
         fOptionsMap = new TMap;
         fOptionsMap->SetOwnerKeyValue();
      }
      if (objTags->GetEntries() == 2) {
         TString key = ((TObjString *) objTags->At(0))->GetName();
         TString value = ((TObjString *) objTags->At(1))->GetName();
         fOptionsMap->Add(new TObjString(key), new TObjString(value));
      } else {
         TString key = ((TObjString *) objTags->At(0))->GetName();
         fOptionsMap->Add(new TObjString(key), nullptr);
      }
      delete objTags;
   }
   delete objOptions;
}


////////////////////////////////////////////////////////////////////////////////
/// Return a value for a given key from the URL options.
/// Returns 0 in case key is not found.

const char *TUrl::GetValueFromOptions(const char *key) const
{
   if (!key) return nullptr;
   ParseOptions();
   TObject *option = fOptionsMap ? fOptionsMap->GetValue(key) : nullptr;
   return option ? option->GetName() : nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a value for a given key from the URL options as an Int_t,
/// a missing key returns -1.

Int_t TUrl::GetIntValueFromOptions(const char *key) const
{
   if (!key) return -1;
   ParseOptions();
   TObject *option = fOptionsMap ? fOptionsMap->GetValue(key) : nullptr;
   return option ? atoi(option->GetName()) : -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if the given key appears in the URL options list.

Bool_t TUrl::HasOption(const char *key) const
{
   if (!key) return kFALSE;
   ParseOptions();

   return fOptionsMap && fOptionsMap->FindObject(key);
}

////////////////////////////////////////////////////////////////////////////////
/// Recompute the path removing all relative directory jumps via '..'.

void TUrl::CleanRelativePath()
{
   Ssiz_t slash = 0;
   while ( (slash = fFile.Index("/..") ) != kNPOS) {
      // find backwards the next '/'
      Bool_t found = kFALSE;
      for (int l = slash-1; l >=0; l--) {
         if (fFile[l] == '/') {
            // found previous '/'
            fFile.Remove(l, slash+3-l);
            found = kTRUE;
            break;
         }
      }
      if (!found)
        break;
   }
}
