// @(#)root/net:$Id$
// Author: Adrien Devresse and Tigran Mkrtchyan

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TDavixFile                                                           //
//                                                                      //
// A TDavixFile is like a normal TFile except that it uses              //
// libdavix to read/write remote files.                                 //
// It supports HTTP and HTTPS in a number of dialects and options       //
//  e.g. S3 is one of them                                              //
// Other caracteristics come from the full support of Davix,            //
//  e.g. full redirection support in any circumstance                   //
//                                                                      //
// Authors:     Adrien Devresse (CERN IT/SDC)                           //
//              Tigran Mkrtchyan (DESY)                                 //
//                                                                      //
// Checks and ROOT5 porting:                                            //
//              Fabrizio Furano (CERN IT/SDC)                           //
//                                                                      //
// September 2013                                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "ROOT/RLogger.hxx"
#include "TDavixFile.h"
#include "TROOT.h"
#include "TSocket.h"
#include "Bytes.h"
#include "TError.h"
#include "TEnv.h"
#include "TBase64.h"
#include "TVirtualPerfStats.h"
#include "TDavixFileInternal.h"
#include "snprintf.h"

#include <cerrno>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <davix.hpp>
#include <sstream>
#include <string>
#include <cstring>


static const std::string VERSION = "0.2.0";

static const std::string gUserAgent = "ROOT/" + std::string(gROOT->GetVersion()) +
" TDavixFile/" + VERSION + " davix/" + Davix::version();

// The prefix that is used to find the variables in the gEnv
#define ENVPFX "Davix."

ClassImp(TDavixFile);

using namespace Davix;

const char* grid_mode_opt = "grid_mode=yes";
const char* ca_check_opt = "ca_check=no";
const char* s3_seckey_opt = "s3seckey=";
const char* s3_acckey_opt = "s3acckey=";
const char* s3_region_opt = "s3region=";
const char* s3_token_opt = "s3token=";
const char* s3_alternate_opt = "s3alternate=";
const char* open_mode_read = "READ";
const char* open_mode_create = "CREATE";
const char* open_mode_new = "NEW";
const char* open_mode_update = "UPDATE";

static TMutex createLock;
static Context* davix_context_s = NULL;

////////////////////////////////////////////////////////////////////////////////

ROOT::Experimental::RLogChannel &TDavixLogChannel()
{
   static ROOT::Experimental::RLogChannel sLog("ROOT.TDavix");
   return sLog;
}

////////////////////////////////////////////////////////////////////////////////

bool isno(const char *str)
{
   if (!str) return false;

   if (!strcmp(str, "n") || !strcmp(str, "no") || !strcmp(str, "0") || !strcmp(str, "false")) return true;

   return false;

}

bool strToBool(const char *str, bool defvalue) {
    if(!str) return defvalue;

    if(strcmp(str, "n") == 0 || strcmp(str, "no") == 0  || strcmp(str, "0") == 0 || strcmp(str, "false") == 0) return false;
    if(strcmp(str, "y") == 0 || strcmp(str, "yes") == 0 || strcmp(str, "1") == 0 || strcmp(str, "true") == 0)  return true;

    return defvalue;
}

////////////////////////////////////////////////////////////////////////////////

int configure_open_flag(const std::string &str, int old_flag)
{
   if (strcasecmp(str.c_str(), open_mode_read) == 0)
      old_flag |= O_RDONLY;
   if ((strcasecmp(str.c_str(), open_mode_create) == 0)
         || (strcasecmp(str.c_str(), open_mode_new) == 0)) {
      old_flag |= (O_CREAT | O_WRONLY | O_TRUNC);
   }
   if ((strcasecmp(str.c_str(), open_mode_update) == 0)) {
      old_flag |= (O_RDWR);
   }
   return old_flag;
}

////////////////////////////////////////////////////////////////////////////////

static void ConfigureDavixLogLevel()
{
   Int_t log_level = (gEnv) ? gEnv->GetValue("Davix.Debug", 0) : 0;

   switch (log_level) {
      case 0:
         davix_set_log_level(0);
         break;
      case 1:
         davix_set_log_level(DAVIX_LOG_WARNING);
         break;
      case 2:
         davix_set_log_level(DAVIX_LOG_VERBOSE);
         break;
      case 3:
         davix_set_log_level(DAVIX_LOG_DEBUG);
         break;
      default:
         davix_set_log_level(DAVIX_LOG_ALL);
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////

bool normalizeToken(const std::string &input_token, std::string &output_token)
{
   static const std::string whitespace = " \t\f\n\v\r";
   // Tokens are meant for use in an HTTP header; these are forbidden
   // characters.
   static const std::string nonheader_whitespace = "\r\n";
   auto begin = input_token.find_first_not_of(whitespace);
   // No token here, but not an error.
   if (begin == std::string::npos) {
       output_token = "";
       return true;
   }

   std::string token = input_token.substr(begin);
   auto end = token.find_last_not_of(whitespace);
   token = token.substr(0, end + 1);

   // If non-permitted header characters are present ("\r\n"),
   // then this token is not permitted.
   if (token.find(nonheader_whitespace) != std::string::npos) {
       output_token = "";
       R__LOG_ERROR(TDavixLogChannel()) << "Token discovery failure: token contains non-permitted character sequence (\\r\\n)";
       return false;
   }
   output_token = token;
   return true;
}

bool findTokenInFile(const std::string &token_file, std::string &output_token)
{
   R__LOG_INFO(TDavixLogChannel()) << "Looking for token in file " << token_file.c_str();
   int fd = open(token_file.c_str(), O_RDONLY);
   if (fd == -1) {
      output_token = "";
      if (errno == ENOENT) {
         return true;
      }
      R__LOG_ERROR(TDavixLogChannel()) << "Cannot open '" << token_file << "', error: " << strerror(errno);
      return false;
   }

   // As a pragmatic matter, we limit tokens to 16KB; this is often the limit
   // in size of HTTP headers for many web servers.
   // This also avoids unbounded memory use in case if the user points us to
   // a large file.
   static const size_t max_size = 16384;
   std::vector<char> input_buffer;
   input_buffer.resize(max_size);
   ssize_t retval = read(fd, &input_buffer[0], max_size);
   close(fd);

   if (retval == -1) {
      output_token = "";
      R__LOG_ERROR(TDavixLogChannel()) << "Token discovery failure: failed to read file " << token_file.c_str() << "', error: " << strerror(errno);
      return false;
   }
   if (retval == 16384) {
      // Token may have been truncated!  Erring on the side of caution.
      R__LOG_ERROR(TDavixLogChannel()) << "Token discovery failure: token was larger than 16KB limit.";
      return false;
   }

   std::string token(&input_buffer[0], retval);

   return normalizeToken(token, output_token);
}

// Token discovery process is based on WLCG specification:
// https://github.com/WLCG-AuthZ-WG/bearer-token-discovery/blob/master/specification.md
std::string DiscoverToken()
{
   const char *bearer_token = getenv("BEARER_TOKEN");
   std::string token;
   if (bearer_token && *bearer_token){
      if (!normalizeToken(bearer_token, token)) {return "";}
      if (!token.empty()) {return token;}
   }

   const char *bearer_token_file = getenv("BEARER_TOKEN_FILE");
   if (bearer_token_file) {
      if (!findTokenInFile(bearer_token_file, token)) {return "";}
      if (!token.empty()) {return token;}
   }

#ifndef WIN32
   uid_t euid = geteuid();
   std::string fname = "/bt_u";
   fname += std::to_string(euid);

   const char *xdg_runtime_dir = getenv("XDG_RUNTIME_DIR");
   if (xdg_runtime_dir) {
      std::string xdg_token_file = std::string(xdg_runtime_dir) + fname;
      if (!findTokenInFile(xdg_token_file, token)) {return "";}
      if (!token.empty()) {return token;}
   }

   if (!findTokenInFile("/tmp" + fname, token)) {return "";}
   return token;
#else
   // WLCG profile doesn't define search paths on Windows;
   // skip this for now.
   return "";
#endif
}

///////////////////////////////////////////////////////////////////
// Authn implementation, Locate and get VOMS cred if exist

////////////////////////////////////////////////////////////////////////////////

static void TDavixFile_http_get_ucert(std::string &ucert, std::string &ukey)
{
   char default_proxy[64];
   const char *genvvar = 0, *genvvar1 = 0;
   // The gEnv has higher priority, let's look for a proxy cert
   genvvar = gEnv->GetValue("Davix.GSI.UserProxy", (const char *) NULL);
   if (genvvar) {
      ucert = ukey = genvvar;
      if (gDebug > 0)
         Info("TDavixFile_http_get_ucert", "Found proxy in gEnv");
      return;
   }

   // Try explicit environment for proxy
   if (getenv("X509_USER_PROXY")) {
      if (gDebug > 0)
         Info("TDavixFile_http_get_ucert", "Found proxy in X509_USER_PROXY");
      ucert = ukey = getenv("X509_USER_PROXY");
      return;
   }

   // Try with default location
   snprintf(default_proxy, sizeof(default_proxy), "/tmp/x509up_u%d",
            geteuid());

   if (access(default_proxy, R_OK) == 0) {
      if (gDebug > 0)
         Info("TDavixFile_http_get_ucert", "Found proxy in /tmp");
      ucert = ukey = default_proxy;
      return;
   }

   // It seems we got no proxy, let's try to gather the keys
   genvvar = gEnv->GetValue("Davix.GSI.UserCert", (const char *) NULL);
   genvvar1 = gEnv->GetValue("Davix.GSI.UserKey", (const char *) NULL);
   if (genvvar || genvvar1) {
      if (gDebug > 0)
         Info("TDavixFile_http_get_ucert", "Found cert and key in gEnv");

      ucert = genvvar;
      ukey = genvvar1;
      return;
   }

   // try with X509_* environment
   if (getenv("X509_USER_CERT"))
      ucert = getenv("X509_USER_CERT");
   if (getenv("X509_USER_KEY"))
      ukey = getenv("X509_USER_KEY");

   if ((ucert.size() > 0) || (ukey.size() > 0)) {
      if (gDebug > 0)
         Info("TDavixFile_http_get_ucert", "Found cert and key in gEnv");
   }
   return;

}

////////////////////////////////////////////////////////////////////////////////

static int TDavixFile_http_authn_cert_X509(void *userdata, const Davix::SessionInfo &info,
      Davix::X509Credential *cert, Davix::DavixError **err)
{
   (void) userdata; // keep quiete compilation warnings
   (void) info;
   std::string ucert, ukey;
   TDavixFile_http_get_ucert(ucert, ukey);

   if (ucert.empty() || ukey.empty()) {
      Davix::DavixError::setupError(err, "TDavixFile",
                                    Davix::StatusCode::AuthentificationError,
                                    "Could not set the user's proxy or certificate");
      return -1;
   }
   return cert->loadFromFilePEM(ukey, ucert, "", err);
}
/////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TDavixFileInternal::~TDavixFileInternal()
{
   delete davixPosix;
   delete davixParam;
}

////////////////////////////////////////////////////////////////////////////////

Context *TDavixFileInternal::getDavixInstance()
{
   if (davix_context_s == NULL) {
      TLockGuard guard(&createLock);
      if (davix_context_s == NULL) {
         davix_context_s = new Context();
      }
   }
   return davix_context_s;
}

////////////////////////////////////////////////////////////////////////////////

Davix_fd *TDavixFileInternal::Open()
{
   DavixError *davixErr = NULL;
   Davix_fd *fd = davixPosix->open(davixParam, fUrl.GetUrl(), oflags, &davixErr);
   if (fd == NULL) {
       // An error has occurred.. We might be able to recover with metalinks.
       // Try to populate the replicas vector. If successful, TFile will try
       // the replicas one by one

       replicas.clear();
       DavixError *davixErr2 = NULL;
       try {
           DavFile file(*davixContext, Davix::Uri(fUrl.GetUrl()));
           std::vector<DavFile> replicasLocal = file.getReplicas(NULL, &davixErr2);
           for(size_t i = 0; i < replicasLocal.size(); i++) {
             replicas.push_back(replicasLocal[i].getUri().getString());
           }
       }
       catch(...) {}
       DavixError::clearError(&davixErr2);

       if(replicas.empty()) {
           // I was unable to retrieve a list of replicas: propagate the original
           // error.
           Error("DavixOpen", "can not open file \"%s\" with davix: %s (%d)",
                 fUrl.GetUrl(),
                 davixErr->getErrMsg().c_str(), davixErr->getStatus());
        }
        DavixError::clearError(&davixErr);
   } else {
      // setup ROOT style read
      davixPosix->fadvise(fd, 0, 300, Davix::AdviseRandom);
   }

   return fd;
}

////////////////////////////////////////////////////////////////////////////////

void TDavixFileInternal::Close()
{
   DavixError *davixErr = NULL;
   if (davixFd != NULL && davixPosix->close(davixFd, &davixErr)) {
      Error("DavixClose", "can not to close file with davix: %s (%d)",
            davixErr->getErrMsg().c_str(), davixErr->getStatus());
      DavixError::clearError(&davixErr);
   }
}

////////////////////////////////////////////////////////////////////////////////

void TDavixFileInternal::enableGridMode()
{
   const char *env_var = NULL;

   if (gDebug > 1)
      Info("enableGridMode", " grid mode enabled !");

   if( ( env_var = getenv("X509_CERT_DIR")) == NULL){
      env_var= "/etc/grid-security/certificates/";
   }
   davixParam->addCertificateAuthorityPath(env_var);
   if (gDebug > 0)
      Info("enableGridMode", "Adding CAdir %s", env_var);
}

////////////////////////////////////////////////////////////////////////////////

// Only newer versions of davix support setting the S3 region and STS tokens.
// But it's only possible to check the davix version through a #define starting from
// 0.6.4.
// I have no way to check if setAwsRegion is available, so let's use SFINAE. :-)
// The first overload will always take priority - if "substitution" fails, meaning
// setAwsRegion is not there, the compiler will pick the second overload with
// the ellipses. (...)

template<typename TRequestParams = Davix::RequestParams>
static auto awsRegion(TRequestParams *parameters, const char *region)
  -> decltype(parameters->setAwsRegion(region), void())
{
   if (gDebug > 1) Info("awsRegion", "Setting S3 Region to '%s' - v4 signature will be used", region);
   parameters->setAwsRegion(region);
}

template<typename TRequestParams = Davix::RequestParams>
static void awsRegion(...) {
   Warning("setAwsRegion", "Unable to set AWS region, not supported by this version of davix");
}

// Identical SFINAE trick as above for setAwsToken
template<typename TRequestParams = Davix::RequestParams>
static auto awsToken(TRequestParams *parameters, const char *token)
  -> decltype(parameters->setAwsToken(token), void())
{
   if (gDebug > 1) Info("awsToken", "Setting S3 STS temporary credentials");
   parameters->setAwsToken(token);
}

template<typename TRequestParams = Davix::RequestParams>
static void awsToken(...) {
   Warning("awsToken", "Unable to set AWS token, not supported by this version of davix");
}

// Identical SFINAE trick as above for setAwsAlternate
template<typename TRequestParams = Davix::RequestParams>
static auto awsAlternate(TRequestParams *parameters, bool option)
  -> decltype(parameters->setAwsAlternate(option), void())
{
   if (gDebug > 1) Info("awsAlternate", "Setting S3 path-based bucket option (s3alternate)");
   parameters->setAwsAlternate(option);
}

template<typename TRequestParams = Davix::RequestParams>
static void awsAlternate(...) {
   Warning("awsAlternate", "Unable to set AWS path-based bucket option (s3alternate), not supported by this version of davix");
}

void TDavixFileInternal::setAwsRegion(const std::string & region) {
   if(!region.empty()) {
      awsRegion(davixParam, region.c_str());
   }
}

void TDavixFileInternal::setAwsToken(const std::string & token) {
   if(!token.empty()) {
      awsToken(davixParam, token.c_str());
   }
}

void TDavixFileInternal::setAwsAlternate(const bool & option) {
   awsAlternate(davixParam, option);
}


void TDavixFileInternal::setS3Auth(const std::string &secret, const std::string &access,
                                   const std::string &region, const std::string &token)
{
   if (gDebug > 1) {
      Info("setS3Auth", " Aws S3 tokens configured");
   }
   davixParam->setAwsAuthorizationKeys(secret, access);
   davixParam->setProtocol(RequestProtocol::AwsS3);

   setAwsRegion(region);
   setAwsToken(token);
}

////////////////////////////////////////////////////////////////////////////////

void TDavixFileInternal::parseConfig()
{
   const char *env_var = NULL, *env_var2 = NULL;
   // default opts
   davixParam->setTransparentRedirectionSupport(true);
   davixParam->setClientCertCallbackX509(&TDavixFile_http_authn_cert_X509, NULL);

   // setup CADIR
   env_var = gEnv->GetValue("Davix.GSI.CAdir", (const char *) NULL);
   if (env_var) {
      davixParam->addCertificateAuthorityPath(env_var);
      if (gDebug > 0)
         Info("parseConfig", "Add CAdir: %s", env_var);
   }

   // CA Check
   bool ca_check_local = !isno(gEnv->GetValue("Davix.GSI.CACheck", (const char *)"y"));
   davixParam->setSSLCAcheck(ca_check_local);
   if (gDebug > 0)
      Info("parseConfig", "Setting CAcheck to %s", ((ca_check_local) ? ("true") : ("false")));

   // WLCG Bearer tokens check
   std::string prefix = "Bearer ";
   auto token = DiscoverToken();
   if (token.c_str()) {
      // header: "Authorization: Bearer mytoken"
      davixParam->addHeader("Authorization", prefix + token);
   }

   // S3 Auth
   if (((env_var = gEnv->GetValue("Davix.S3.SecretKey", getenv("S3_SECRET_KEY"))) != NULL)
         && ((env_var2 = gEnv->GetValue("Davix.S3.AccessKey", getenv("S3_ACCESS_KEY"))) != NULL)) {
      Info("parseConfig", "Setting S3 SecretKey and AccessKey. Access Key : %s ", env_var2);
      davixParam->setAwsAuthorizationKeys(env_var, env_var2);

      // need to set region?
      if ( (env_var = gEnv->GetValue("Davix.S3.Region", getenv("S3_REGION"))) != NULL) {
         setAwsRegion(env_var);
      }
      // need to set STS token?
      if( (env_var = gEnv->GetValue("Davix.S3.Token", getenv("S3_TOKEN"))) != NULL) {
         setAwsToken(env_var);
      }
      // need to set aws alternate?
      if( (env_var = gEnv->GetValue("Davix.S3.Alternate", getenv("S3_ALTERNATE"))) != NULL) {
         setAwsAlternate(strToBool(env_var, false));
      }
   }

   env_var = gEnv->GetValue("Davix.GSI.GridMode", (const char *)"y");
   if (!isno(env_var))
      enableGridMode();
}

////////////////////////////////////////////////////////////////////////////////
/// intput params

void TDavixFileInternal::parseParams(Option_t *option)
{
   std::stringstream ss(option);
   std::string item;
   std::vector<std::string> parsed_options;
   // parameters
   std::string s3seckey, s3acckey, s3region, s3token;

   while (std::getline(ss, item, ' ')) {
      parsed_options.push_back(item);
   }

   for (std::vector<std::string>::iterator it = parsed_options.begin(); it < parsed_options.end(); ++it) {
      // grid mode option
      if ((strcasecmp(it->c_str(), grid_mode_opt)) == 0) {
         enableGridMode();
      }
      // ca check option
      if ((strcasecmp(it->c_str(), ca_check_opt)) == 0) {
         davixParam->setSSLCAcheck(false);
      }
      // s3 sec key
      if (strncasecmp(it->c_str(), s3_seckey_opt, strlen(s3_seckey_opt)) == 0) {
         s3seckey = std::string(it->c_str() + strlen(s3_seckey_opt));
      }
      // s3 access key
      if (strncasecmp(it->c_str(), s3_acckey_opt, strlen(s3_acckey_opt)) == 0) {
         s3acckey = std::string(it->c_str() + strlen(s3_acckey_opt));
      }
      // s3 region
      if (strncasecmp(it->c_str(), s3_region_opt, strlen(s3_region_opt)) == 0) {
         s3region = std::string(it->c_str() + strlen(s3_region_opt));
      }
      // s3 sts token
      if (strncasecmp(it->c_str(), s3_token_opt, strlen(s3_token_opt)) == 0) {
         s3token = std::string(it->c_str() + strlen(s3_token_opt));
      }
      // s3 alternate option
      if (strncasecmp(it->c_str(), s3_alternate_opt, strlen(s3_alternate_opt)) == 0) {
         setAwsAlternate(strToBool(it->c_str() + strlen(s3_alternate_opt), false));
      }
      // open mods
      oflags = configure_open_flag(*it, oflags);
   }

   if (s3seckey.size() > 0) {
      setS3Auth(s3seckey, s3acckey, s3region, s3token);
   }

   if (oflags == 0) // default open mode
      oflags = O_RDONLY;
}

////////////////////////////////////////////////////////////////////////////////

void TDavixFileInternal::init()
{
   davixPosix = new DavPosix(davixContext);
   davixParam = new RequestParams();
   davixParam->setUserAgent(gUserAgent);
   davixParam->setMetalinkMode(Davix::MetalinkMode::Disable);
   ConfigureDavixLogLevel();
   parseConfig();
   parseParams(opt);
}

////////////////////////////////////////////////////////////////////////////////

Int_t TDavixFileInternal::DavixStat(const char *url, struct stat *st)
{
   DavixError *davixErr = NULL;

   if (davixPosix->stat(davixParam, url, st, &davixErr) < 0) {

      Error("DavixStat", "can not stat the file with davix: %s (%d)",
            davixErr->getErrMsg().c_str(), davixErr->getStatus());
      DavixError::clearError(&davixErr);
      return 0;
   }
   return 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TDavixFile::TDavixFile(const char *url, Option_t *opt, const char *ftitle, Int_t compress) : TFile(url, "WEB"),
   d_ptr(new TDavixFileInternal(fUrl, opt))
{
   (void) ftitle;
   (void) compress;
   Init(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////

TDavixFile::~TDavixFile()
{
   d_ptr->Close();
   delete d_ptr;
}

////////////////////////////////////////////////////////////////////////////////

void TDavixFile::Init(Bool_t init)
{
   (void) init;
   //initialize davix
   d_ptr->init();
   // pre-open file
   if ((d_ptr->getDavixFileInstance()) == NULL){
         MakeZombie();
         gDirectory = gROOT;
         return;
   }
   TFile::Init(kFALSE);
   fOffset = 0;
   fD = -2; // so TFile::IsOpen() will return true when in TFile::~TFi */
}

TString TDavixFile::GetNewUrl() {
   std::vector<std::string> replicas = d_ptr->getReplicas();
   TString newUrl;
   if(!replicas.empty()) {
      std::stringstream ss;
      for(size_t i = 0; i < replicas.size(); i++) {
         ss << replicas[i];
         if(i != replicas.size()-1) ss << "|";
      }
      newUrl = ss.str();
   }
   return newUrl;
}

////////////////////////////////////////////////////////////////////////////////
/// Set position from where to start reading.

void TDavixFile::Seek(Long64_t offset, ERelativeTo pos)
{
   TLockGuard guard(&(d_ptr->positionLock));
   switch (pos) {
      case kBeg:
         fOffset = offset + fArchiveOffset;
         break;
      case kCur:
         fOffset += offset;
         break;
      case kEnd:
         // this option is not used currently in the ROOT code
         if (fArchiveOffset)
            Error("Seek", "seeking from end in archive is not (yet) supported");
         fOffset = fEND - offset; // is fEND really EOF or logical EOF?
         break;
   }

   if (gDebug > 1)
      Info("Seek", " move cursor to %lld"
           , fOffset);
}

////////////////////////////////////////////////////////////////////////////////
/// Read specified byte range from remote file via HTTP.
/// Returns kTRUE in case of error.

Bool_t TDavixFile::ReadBuffer(char *buf, Int_t len)
{
   TLockGuard guard(&(d_ptr->positionLock));
   Davix_fd *fd;
   if ((fd = d_ptr->getDavixFileInstance()) == NULL)
      return kTRUE;
   Long64_t ret = DavixReadBuffer(fd, buf, len);
   if (ret < 0)
      return kTRUE;

   if (gDebug > 1)
      Info("ReadBuffer", "%lld bytes of data read sequentially"
           " (%d requested)", ret, len);

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TDavixFile::ReadBuffer(char *buf, Long64_t pos, Int_t len)
{
   Davix_fd *fd;
   if ((fd = d_ptr->getDavixFileInstance()) == NULL)
      return kTRUE;

   Long64_t ret = DavixPReadBuffer(fd, buf, pos, len);
   if (ret < 0)
      return kTRUE;

   if (gDebug > 1)
      Info("ReadBuffer", "%lld bytes of data read from offset"
           " %lld (%d requested)", ret, pos, len);
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TDavixFile::ReadBufferAsync(Long64_t offs, Int_t len)
{
   Davix_fd *fd;
   if ((fd = d_ptr->getDavixFileInstance()) == NULL)
      return kFALSE;

   d_ptr->davixPosix->fadvise(fd, static_cast<dav_off_t>(offs), static_cast<dav_size_t>(len), Davix::AdviseRandom);

   if (gDebug > 1)
      Info("ReadBufferAsync", "%d bytes of data prefected from offset"
           " %lld ",  len, offs);
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TDavixFile::ReadBuffers(char *buf, Long64_t *pos, Int_t *len, Int_t nbuf)
{
   Davix_fd *fd;
   if ((fd = d_ptr->getDavixFileInstance()) == NULL)
      return kTRUE;

   Long64_t ret = DavixReadBuffers(fd, buf, pos, len, nbuf);
   if (ret < 0)
      return kTRUE;

   if (gDebug > 1)
      Info("ReadBuffers", "%lld bytes of data read from a list of %d buffers",
           ret, nbuf);

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TDavixFile::WriteBuffer(const char *buf, Int_t len)
{
   Davix_fd *fd;
   if ((fd = d_ptr->getDavixFileInstance()) == NULL)
      return kTRUE;

   Long64_t ret = DavixWriteBuffer(fd, buf, len);
   if (ret < 0)
      return kTRUE;

   if (gDebug > 1)
      Info("WriteBuffer", "%lld bytes of data write"
           " %d requested", ret, len);
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

void TDavixFile::setCACheck(Bool_t check)
{
   d_ptr->davixParam->setSSLCAcheck((bool)check);
}

////////////////////////////////////////////////////////////////////////////////

void TDavixFile::enableGridMode()
{
   d_ptr->enableGridMode();
}

////////////////////////////////////////////////////////////////////////////////

bool TDavixFileInternal::isMyDird(void *fd)
{
   TLockGuard l(&(openLock));
   std::vector<void *>::iterator f = std::find(dirdVec.begin(), dirdVec.end(), fd);
   return (f != dirdVec.end());
}

////////////////////////////////////////////////////////////////////////////////

void TDavixFileInternal::addDird(void *fd)
{
   TLockGuard l(&(openLock));
   dirdVec.push_back(fd);
}

////////////////////////////////////////////////////////////////////////////////

void TDavixFileInternal::removeDird(void *fd)
{
   TLockGuard l(&(openLock));
   std::vector<void *>::iterator f = std::find(dirdVec.begin(), dirdVec.end(), fd);
   if (f != dirdVec.end())
      dirdVec.erase(f);
}

////////////////////////////////////////////////////////////////////////////////

Long64_t TDavixFile::GetSize() const
{
   struct stat st;
   Int_t ret = d_ptr->DavixStat(fUrl.GetUrl(), &st);
   if (ret) {
      if (gDebug > 1)
         Info("GetSize", "file size requested:  %lld", (Long64_t)st.st_size);
      return st.st_size;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////

Double_t TDavixFile::eventStart()
{
   if (gPerfStats)
      return TTimeStamp();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// set TFile state info

void TDavixFile::eventStop(Double_t t_start, Long64_t len, bool read)
{
  if(read) {
   fBytesRead += len;
   fReadCalls += 1;

   SetFileBytesRead(GetFileBytesRead() + len);
   SetFileReadCalls(GetFileReadCalls() + 1);

   if (gPerfStats)
      gPerfStats->FileReadEvent(this, (Int_t) len, t_start);
  } else {
    fBytesWrite += len;
    SetFileBytesWritten(GetFileBytesWritten() + len);
  }
}

////////////////////////////////////////////////////////////////////////////////

Long64_t TDavixFile::DavixReadBuffer(Davix_fd *fd, char *buf, Int_t len)
{
   DavixError *davixErr = NULL;
   Double_t start_time = eventStart();

   Long64_t ret = d_ptr->davixPosix->pread(fd, buf, len, fOffset, &davixErr);
   if (ret < 0) {
      Error("DavixReadBuffer", "can not read data with davix: %s (%d)",
            davixErr->getErrMsg().c_str(), davixErr->getStatus());
      DavixError::clearError(&davixErr);
   } else {
      fOffset += ret;
      eventStop(start_time, ret);
   }

   return ret;
}

////////////////////////////////////////////////////////////////////////////////

Long64_t TDavixFile::DavixWriteBuffer(Davix_fd *fd, const char *buf, Int_t len)
{
   DavixError *davixErr = NULL;
   Double_t start_time = eventStart();

   Long64_t ret = d_ptr->davixPosix->pwrite(fd, buf, len, fOffset, &davixErr);
   if (ret < 0) {
      Error("DavixWriteBuffer", "can not write data with davix: %s (%d)",
            davixErr->getErrMsg().c_str(), davixErr->getStatus());
      DavixError::clearError(&davixErr);
   } else {
      fOffset += ret;
      eventStop(start_time, ret, false);
   }

   return ret;
}

////////////////////////////////////////////////////////////////////////////////

Long64_t TDavixFile::DavixPReadBuffer(Davix_fd *fd, char *buf, Long64_t pos, Int_t len)
{
   DavixError *davixErr = NULL;
   Double_t start_time = eventStart();

   Long64_t ret = d_ptr->davixPosix->pread(fd, buf, len, pos, &davixErr);
   if (ret < 0) {
      Error("DavixPReadBuffer", "can not read data with davix: %s (%d)",
            davixErr->getErrMsg().c_str(), davixErr->getStatus());
      DavixError::clearError(&davixErr);
   } else {
      eventStop(start_time, ret);
   }


   return ret;
}

////////////////////////////////////////////////////////////////////////////////

Long64_t TDavixFile::DavixReadBuffers(Davix_fd *fd, char *buf, Long64_t *pos, Int_t *len, Int_t nbuf)
{
   DavixError *davixErr = NULL;
   Double_t start_time = eventStart();
   DavIOVecInput in[nbuf];
   DavIOVecOuput out[nbuf];

   int lastPos = 0;
   for (Int_t i = 0; i < nbuf; ++i) {
      in[i].diov_buffer = &buf[lastPos];
      in[i].diov_offset = pos[i];
      in[i].diov_size = len[i];
      lastPos += len[i];
   }

   Long64_t ret = d_ptr->davixPosix->preadVec(fd, in, out, nbuf, &davixErr);
   if (ret < 0) {
      Error("DavixReadBuffers", "can not read data with davix: %s (%d)",
            davixErr->getErrMsg().c_str(), davixErr->getStatus());
      DavixError::clearError(&davixErr);
   } else {
      eventStop(start_time, ret);
   }

   return ret;
}
