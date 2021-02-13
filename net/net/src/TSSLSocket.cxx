// @(#)root/net:$Id: TSSLSocket.cxx
// Author: Alejandro Alvarez 16/09/2011

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TSSLSocket                                                           //
//                                                                      //
// A TSocket wrapped in by SSL.                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <openssl/ssl.h>
#include <cstdio>
#include "strlcpy.h"
#include "TSSLSocket.h"
#include "TSystem.h"

// Static properties
char TSSLSocket::fgSSLCAFile[FILENAME_MAX] = "";
char TSSLSocket::fgSSLCAPath[FILENAME_MAX] = "";
char TSSLSocket::fgSSLUCert[FILENAME_MAX]  = "";
char TSSLSocket::fgSSLUKey[FILENAME_MAX]   = "";

////////////////////////////////////////////////////////////////////////////////
/// Wraps the socket with OpenSSL.

void TSSLSocket::WrapWithSSL(void)
{
   SSL_library_init();

   // New context
   if (!(fSSLCtx = SSL_CTX_new(SSLv23_method()))) {
      Error("WrapWithSSL", "the context could not be created");
      goto wrapFailed;
   }

   if ((fgSSLCAFile[0] || fgSSLCAPath[0]) && SSL_CTX_load_verify_locations(fSSLCtx, fgSSLCAFile, fgSSLCAPath) == 0) {
      Error("WrapWithSSL", "could not set the CA file and/or the CA path");
      goto wrapFailed;
   }

   if (fgSSLUCert[0] && SSL_CTX_use_certificate_chain_file(fSSLCtx, fgSSLUCert) == 0) {
      Error("WrapWithSSL", "could not set the client certificate");
      goto wrapFailed;
   }

   if (fgSSLUKey[0] && SSL_CTX_use_PrivateKey_file(fSSLCtx, fgSSLUKey, SSL_FILETYPE_PEM) == 0) {
      Error("WrapWithSSL", "could not set the client private key");
      goto wrapFailed;
   }

   // New SSL structure
   if (!(fSSL = SSL_new(fSSLCtx))) {
      Error("WrapWithSSL", "cannot create the ssl struct");
      goto wrapFailed;
   }

   // Bind to the socket
   if (SSL_set_fd(fSSL, fSocket) != 1) {
      Error("WrapWithSSL", "cannot bind to the socket %d", fSocket);
      goto wrapFailed;
   }

   // Open connection
   if (SSL_connect(fSSL) != 1) {
      Error("WrapWithSSL", "cannot connect");
      goto wrapFailed;
   }

   return;

wrapFailed:
   Close();
   return;
}

////////////////////////////////////////////////////////////////////////////////

ClassImp(TSSLSocket);

////////////////////////////////////////////////////////////////////////////////

TSSLSocket::TSSLSocket(TInetAddress addr, const char *service, Int_t tcpwindowsize)
   : TSocket(addr, service, tcpwindowsize)
{
   WrapWithSSL();
}

////////////////////////////////////////////////////////////////////////////////

TSSLSocket::TSSLSocket(TInetAddress addr, Int_t port, Int_t tcpwindowsize)
   : TSocket(addr, port, tcpwindowsize)
{
   WrapWithSSL();
}

////////////////////////////////////////////////////////////////////////////////

TSSLSocket::TSSLSocket(const char *host, const char *service, Int_t tcpwindowsize)
   : TSocket(host, service, tcpwindowsize)
{
   WrapWithSSL();
}

////////////////////////////////////////////////////////////////////////////////

TSSLSocket::TSSLSocket(const char *url, Int_t port, Int_t tcpwindowsize)
   : TSocket(url, port, tcpwindowsize)
{
   WrapWithSSL();
}

////////////////////////////////////////////////////////////////////////////////

TSSLSocket::TSSLSocket(const char *sockpath) : TSocket(sockpath)
{
   WrapWithSSL();
}

////////////////////////////////////////////////////////////////////////////////

TSSLSocket::TSSLSocket(Int_t desc) : TSocket(desc)
{
   WrapWithSSL();
}

////////////////////////////////////////////////////////////////////////////////

TSSLSocket::TSSLSocket(Int_t desc, const char *sockpath) : TSocket(desc, sockpath)
{
   WrapWithSSL();
}

////////////////////////////////////////////////////////////////////////////////

TSSLSocket::TSSLSocket(const TSSLSocket &s) : TSocket(s)
{
   WrapWithSSL();
}

////////////////////////////////////////////////////////////////////////////////
/// Close gracefully the connection, and free SSL structures.

TSSLSocket::~TSSLSocket()
{
   Close();
   if (fSSL)
     SSL_free(fSSL);
   if (fSSLCtx)
     SSL_CTX_free(fSSLCtx);
}

////////////////////////////////////////////////////////////////////////////////
/// Close the SSL connection.

void TSSLSocket::Close(Option_t *option)
{
   if (fSSL)
      SSL_shutdown(fSSL);
   TSocket::Close(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Set up the static configuration variables.

void TSSLSocket::SetUpSSL(const char *cafile, const char *capath,
                          const char *ucert,  const char *ukey)
{
   if (cafile)
      strlcpy(fgSSLCAFile, cafile, FILENAME_MAX);
   if (capath)
      strlcpy(fgSSLCAPath, capath, FILENAME_MAX);
   if (ucert)
      strlcpy(fgSSLUCert,  ucert,  FILENAME_MAX);
   if (ukey)
      strlcpy(fgSSLUKey,   ukey,   FILENAME_MAX);
}

////////////////////////////////////////////////////////////////////////////////

Int_t TSSLSocket::Recv(TMessage *& /*mess */)
{
   Error("Recv", "not implemented");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Receive a raw buffer of specified length bytes.

Int_t TSSLSocket::RecvRaw(void *buffer, Int_t length, ESendRecvOptions opt)
{
   TSystem::ResetErrno();

   if (fSocket == -1) return -1;
   if (length == 0)   return 0;

   ResetBit(TSocket::kBrokenConn);

   Int_t n;
   Int_t offset = 0;
   Int_t remain = length;

   // SSL_read/SSL_peek may not return the total length at once
   while (remain > 0) {
     if (opt == kPeek)
        n = SSL_peek(fSSL, (char*)buffer + offset, (int)remain);
     else
        n = SSL_read(fSSL, (char*)buffer + offset, (int)remain);

     if (n <= 0) {
        if (gDebug > 0)
           Error("RecvRaw", "failed to read from the socket");

        if (SSL_get_error(fSSL, n) == SSL_ERROR_ZERO_RETURN || SSL_get_error(fSSL, n) == SSL_ERROR_SYSCALL) {
           // Connection closed, reset or broken
           SetBit(TSocket::kBrokenConn);
           SSL_set_quiet_shutdown(fSSL, 1); // Socket is gone, sending "close notify" will fail
           Close();
        }
        return n;
     }

     // When peeking, just return the available data, don't loop. Otherwise,
     // we may copy the same chunk of data multiple times into the
     // output buffer, for instance when there is no more recent data
     // in the socket's internal reception buffers.
     // Note that in this case we don't update the counters of data received
     // through this socket. They will be updated when the data is actually
     // read. This avoids double counting.
     if (opt == kPeek) return n;

     offset += n;
     remain -= n;
   }

   fBytesRecv  += length;
   fgBytesRecv += length;

   Touch();  // update usage timestamp

   return offset;
}

////////////////////////////////////////////////////////////////////////////////

Int_t TSSLSocket::Send(const TMessage & /* mess */)
{
   Error("Send", "not implemented");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Send a raw buffer of specified length.

Int_t TSSLSocket::SendRaw(const void *buffer, Int_t length, ESendRecvOptions /* opt */)
{
   TSystem::ResetErrno();

   if (fSocket == -1) return -1;

   ResetBit(TSocket::kBrokenConn);

   Int_t nsent;
   if ((nsent = SSL_write(fSSL, buffer, (int)length)) <= 0) {
      if (SSL_get_error(fSSL, nsent) == SSL_ERROR_ZERO_RETURN) {
         // Connection reset or broken: close
         SetBit(TSocket::kBrokenConn);
         Close();
      }
      return nsent;
   }

   fBytesSent  += nsent;
   fgBytesSent += nsent;

   Touch(); // update usage timestamp

   return nsent;
}
