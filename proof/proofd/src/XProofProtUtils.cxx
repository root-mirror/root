// @(#)root/proofd:$Id$
// Author: Gerardo Ganis  12/12/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProofProtUtils.cxx                                                  //
//                                                                      //
// Authors: G. Ganis, CERN 2005                                         //
//                                                                      //
// Utility functions prototypes for client-to-server                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifdef __APPLE__
#   ifndef __macos__
#      define __macos__
#   endif
#endif
#ifdef __sun
#   ifndef __solaris__
#      define __solaris__
#   endif
#endif
#ifndef WIN32
#  include <sys/types.h>
#  ifndef ROOT_XrdFour
#     include <netinet/in.h>
#  endif
#endif
#include "XrdSys/XrdSysPlatform.hh"
#include "XProofProtocol.h"
#include "XProofProtUtils.h"
#include "Bytes.h"

#include <cstdio>

namespace XPD {

////////////////////////////////////////////////////////////////////////////////
/// This function applies the network byte order on those
/// parts of the 16-bytes buffer, only if it is composed
/// by some binary part
/// Return 0 if OK, -1 in case the ID is unknown

int clientMarshall(XPClientRequest* str)
{
   switch(str->header.requestid) {

   case kXP_login:
      str->login.pid = htonl(str->login.pid);
      break;
   case kXP_auth:
      // no swap on ASCII fields
      break;
   case kXP_create:
      // no swap on ASCII fields
      str->proof.int1 = htonl(str->proof.int1);
      break;
   case kXP_destroy:
      str->proof.sid = htonl(str->proof.sid);
      break;
   case kXP_attach:
      str->proof.sid = htonl(str->proof.sid);
      break;
   case kXP_detach:
      str->proof.sid = htonl(str->proof.sid);
      break;
   case kXP_cleanup:
      str->proof.sid = htonl(str->proof.sid);
      str->proof.int1 = htonl(str->proof.int1);
      str->proof.int2 = htonl(str->proof.int2);
      break;
   case kXP_sendmsg:
      str->sendrcv.sid = htonl(str->sendrcv.sid);
      str->sendrcv.opt = htonl(str->sendrcv.opt);
      str->sendrcv.cid = htonl(str->sendrcv.cid);
      break;
   case kXP_admin:
      str->proof.sid = htonl(str->proof.sid);
      str->proof.int1 = htonl(str->proof.int1);
      str->proof.int2 = htonl(str->proof.int2);
      str->proof.int3 = htonl(str->proof.int3);
      break;
   case kXP_readbuf:
      str->readbuf.ofs = htonll(str->readbuf.ofs);
      str->readbuf.len = htonl(str->readbuf.len);
      str->readbuf.int1 = htonl(str->readbuf.int1);
      break;
   case kXP_interrupt:
      str->interrupt.sid = htonl(str->interrupt.sid);
      str->interrupt.type = htonl(str->interrupt.type);
      break;
   case kXP_ping:
      str->sendrcv.sid = htonl(str->sendrcv.sid);
      str->sendrcv.opt = htonl(str->sendrcv.opt);
      break;
   case kXP_urgent:
      str->proof.sid = htonl(str->proof.sid);
      str->proof.int1 = htonl(str->proof.int1);
      str->proof.int2 = htonl(str->proof.int2);
      str->proof.int3 = htonl(str->proof.int3);
      break;
   case kXP_touch:
      str->sendrcv.sid = htonl(str->sendrcv.sid);
      break;
   case kXP_ctrlc:
      str->proof.sid = htonl(str->sendrcv.sid);
      break;
   default:
      fprintf(stderr,"clientMarshall: unknown req ID: %d (0x%x)\n",
                      str->header.requestid, str->header.requestid);
      return -1;
      break;
   }

   str->header.requestid = htons(str->header.requestid);
   str->header.dlen      = htonl(str->header.dlen);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////

void clientUnmarshall(struct ServerResponseHeader* str)
{
   str->status = ntohs(str->status);
   str->dlen   = ntohl(str->dlen);
}

////////////////////////////////////////////////////////////////////////////////

void ServerResponseHeader2NetFmt(struct ServerResponseHeader *srh)
{
   srh->status = htons(srh->status);
   srh->dlen   = htonl(srh->dlen);
}

////////////////////////////////////////////////////////////////////////////////

void ServerInitHandShake2HostFmt(struct ServerInitHandShake *srh)
{
   srh->msglen  = ntohl(srh->msglen);
   srh->protover = ntohl(srh->protover);
   srh->msgval  = ntohl(srh->msgval);
}

////////////////////////////////////////////////////////////////////////////////
/// This procedure convert the request code id (an integer defined in
/// XProtocol.hh) in the ascii label (human readable)

char *convertRequestIdToChar(kXR_int16 requestid)
{
   switch(requestid) {

   case kXP_login:
      return (char *)"kXP_login";
   case kXP_auth:
      return (char *)"kXP_auth";
   case kXP_create:
      return (char *)"kXP_create";
   case kXP_destroy:
      return (char *)"kXP_destroy";
   case kXP_attach:
      return (char *)"kXP_attach";
   case kXP_detach:
      return (char *)"kXP_detach";
   case kXP_sendmsg:
      return (char *)"kXP_sendmsg";
   case kXP_admin:
      return (char *)"kXP_admin";
   case kXP_readbuf:
      return (char *)"kXP_readbuf";
   case kXP_interrupt:
      return (char *)"kXP_interrupt";
   case kXP_ping:
      return (char *)"kXP_ping";
   case kXP_cleanup:
      return (char *)"kXP_cleanup";
   case kXP_urgent:
      return (char *)"kXP_urgent";
   case kXP_touch:
      return (char *)"kXP_touch";
   case kXP_ctrlc:
      return (char *)"kXP_ctrlc";
   default:
      return (char *)"kXP_UNKNOWN";
   }
}

////////////////////////////////////////////////////////////////////////////////

char *convertRespStatusToChar(kXR_int16 status)
{
   switch( status) {
   case kXP_ok:
      return (char *)"kXP_ok";
      break;
   case kXP_oksofar:
      return (char *)"kXP_oksofar";
      break;
   case kXP_attn:
      return (char *)"kXP_attn";
      break;
   case kXP_authmore:
      return (char *)"kXP_authmore";
      break;
   case kXP_error:
      return (char *)"kXP_error";
      break;
   case kXP_wait:
      return (char *)"kXP_wait";
      break;
   default:
      return (char *)"kXP_UNKNOWN";
      break;
   }
}

////////////////////////////////////////////////////////////////////////////////

void smartPrintClientHeader(XPClientRequest* hdr)
{
   printf("\n\n================= DUMPING CLIENT REQUEST HEADER =================\n");

   printf("%40s0x%.2x 0x%.2x\n", "ClientHeader.streamid = ",
          hdr->header.streamid[0], hdr->header.streamid[1]);

   printf("%40s%s (%d)\n", "ClientHeader.requestid = ",
          convertRequestIdToChar(hdr->header.requestid), hdr->header.requestid);

   void *tmp;
   switch(hdr->header.requestid) {

   case kXP_login:
      printf("%40s%d \n", "ClientHeader.login.pid = ", hdr->login.pid);
      printf("%40s%s\n", "ClientHeader.login_body.username = ", hdr->login.username);
      tmp = &hdr->login.reserved[0];
      printf("%40s0 repeated %d times\n", "ClientHeader.login.reserved = ",
             *((kXR_int16 *)tmp)); // use tmp to avoid type punned warning
      printf("%40s%d\n", "ClientHeader.login.role = ", (kXR_int32)hdr->login.role[0]);
      break;
   case kXP_auth:
      printf("%40s0 repeated %d times\n", "ClientHeader.auth.reserved = ",
             (kXR_int32)sizeof(hdr->auth.reserved));
      printf("  ClientHeader.auth.credtype= 0x%.2x 0x%.2x 0x%.2x 0x%.2x \n",
             hdr->auth.credtype[0], hdr->auth.credtype[1],
             hdr->auth.credtype[2], hdr->auth.credtype[3]);
      break;
   case kXP_create:
      break;
   case kXP_destroy:
      printf("%40s%d \n", "ClientHeader.proof.sid = ", hdr->proof.sid);
      break;
   case kXP_attach:
      printf("%40s%d \n", "ClientHeader.proof.sid = ", hdr->proof.sid);
      break;
   case kXP_detach:
      printf("%40s%d \n", "ClientHeader.proof.sid = ", hdr->proof.sid);
      break;
   case kXP_ctrlc:
      printf("%40s%d \n", "ClientHeader.proof.sid = ", hdr->proof.sid);
      break;
   case kXP_cleanup:
      printf("%40s%d \n", "ClientHeader.proof.sid = ", hdr->proof.sid);
      printf("%40s%d \n", "ClientHeader.proof.int1 = ", hdr->proof.int1);
      printf("%40s%d \n", "ClientHeader.proof.int2 = ", hdr->proof.int2);
      break;
   case kXP_sendmsg:
      printf("%40s%d \n", "ClientHeader.sendrcv.sid = ", hdr->sendrcv.sid);
      printf("%40s%d \n", "ClientHeader.sendrcv.opt = ", hdr->sendrcv.opt);
      printf("%40s%d \n", "ClientHeader.sendrcv.cid = ", hdr->sendrcv.cid);
      break;
   case kXP_interrupt:
      printf("%40s%d \n", "ClientHeader.interrupt.sid = ", hdr->interrupt.sid);
      printf("%40s%d \n", "ClientHeader.interrupt.type = ", hdr->interrupt.type);
      break;
   case kXP_ping:
      printf("%40s%d \n", "ClientHeader.sendrcv.sid = ", hdr->sendrcv.sid);
      printf("%40s%d \n", "ClientHeader.sendrcv.opt = ", hdr->sendrcv.opt);
      break;
   case kXP_touch:
      printf("%40s%d \n", "ClientHeader.sendrcv.sid = ", hdr->sendrcv.sid);
      break;
   case kXP_admin:
   case kXP_urgent:
      printf("%40s%d \n", "ClientHeader.proof.sid = ", hdr->proof.sid);
      printf("%40s%d \n", "ClientHeader.proof.int1 = ", hdr->proof.int1);
      printf("%40s%d \n", "ClientHeader.proof.int2 = ", hdr->proof.int2);
      printf("%40s%d \n", "ClientHeader.proof.int3 = ", hdr->proof.int3);
      break;
    case kXP_readbuf:
      printf("%40s%lld \n", "ClientHeader.readbuf.ofs = ", hdr->readbuf.ofs);
      printf("%40s%d \n", "ClientHeader.readbuf.len = ", hdr->readbuf.len);
      break;
    default:
      printf("Unknown request ID: %d ! \n", hdr->header.requestid);
  }

   printf("%40s%d", "ClientHeader.header.dlen = ", hdr->header.dlen);
   printf("\n=================== END CLIENT HEADER DUMPING ===================\n\n");
}

////////////////////////////////////////////////////////////////////////////////

void smartPrintServerHeader(struct ServerResponseHeader* hdr)
{
   printf("\n\n======== DUMPING SERVER RESPONSE HEADER ========\n");
   printf("%30s0x%.2x 0x%.2x\n", "ServerHeader.streamid = ",
          hdr->streamid[0], hdr->streamid[1]);
   switch(hdr->status) {
   case kXP_ok:
      printf("%30skXP_ok", "ServerHeader.status = ");
      break;
   case kXP_attn:
      printf("%30skXP_attn", "ServerHeader.status = ");
      break;
   case kXP_authmore:
      printf("%30skXP_authmore", "ServerHeader.status = ");
      break;
   case kXP_error:
      printf("%30skXP_error", "ServerHeader.status = ");
      break;
   case kXP_oksofar:
      printf("%30skXP_oksofar", "ServerHeader.status = ");
      break;
   case kXP_wait:
      printf("%30skXP_wait", "ServerHeader.status = ");
      break;
   }
   printf(" (%d)\n", hdr->status);
   printf("%30s%d", "ServerHeader.dlen = ", hdr->dlen);
   printf("\n========== END DUMPING SERVER HEADER ===========\n\n");
}

} // namespace ROOT
