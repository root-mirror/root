/* @(#)root/multiproc:$Id$ */
// Author: Enrico Guiraud July 2015
// Modified: G Ganis Jan 2017

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "MPCode.h"
#include "MPSendRecv.h"
#include "TError.h"
#include "TMPWorker.h"
#include "TSystem.h"
#include <memory> //unique_ptr
#include <string>

//////////////////////////////////////////////////////////////////////////
///
/// \class TMPWorker
///
/// This class works in conjuction with TMPClient, reacting to messages
/// received from it as specified by the Notify and HandleInput methods.
/// When TMPClient::Fork is called, a TMPWorker instance is passed to it
/// which will take control of the ROOT session in the children processes.
///
/// After forking, every time a message is sent or broadcast to the workers,
/// TMPWorker::Notify is called and the message is retrieved.
/// Messages exchanged between TMPClient and TMPWorker should be sent with
/// the MPSend() standalone function.\n
/// If the code of the message received is above 1000 (i.e. it is an MPCode)
/// the qualified TMPWorker::HandleInput method is called, that takes care
/// of handling the most generic type of messages. Otherwise the unqualified
/// (possibly overridden) version of HandleInput is called, allowing classes
/// that inherit from TMPWorker to manage their own protocol.\n
/// An application's worker class should inherit from TMPWorker and implement
/// a HandleInput method that overrides TMPWorker's.\n
///
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
/// This method is called by children processes right after forking.
/// Initialization of worker properties that must be delayed until after
/// forking must be done here.\n
/// For example, Init saves the pid into fPid, and adds the TMPWorker to
/// the main eventloop (as a TFileHandler).\n
/// Make sure this operations are performed also by overriding implementations,
/// e.g. by calling TMPWorker::Init explicitly.
void TMPWorker::Init(int fd, unsigned workerN)
{
   fS.reset(new TSocket(fd, "MPsock")); //TSocket's constructor with this signature seems much faster than TSocket(int fd)
   fPid = getpid();
   fNWorker = workerN;
   fId = "W" + std::to_string(GetNWorker()) + "|P" + std::to_string(GetPid());
}


void TMPWorker::Run()
{
   while(true) {
      MPCodeBufPair msg = MPRecv(fS.get());
      if (msg.first == MPCode::kRecvError) {
         Error("TMPWorker::Run", "Lost connection to client\n");
         gSystem->Exit(0);
      }

      if (msg.first < 1000)
         HandleInput(msg); //call overridden method
      else
         TMPWorker::HandleInput(msg); //call this class' method
  }
}


//////////////////////////////////////////////////////////////////////////
/// Handle a message with an EMPCode.
/// This method is called upon receiving a message with a code >= 1000 (i.e.
/// EMPCode). It handles the most generic types of messages.\n
/// Classes inheriting from TMPWorker should implement their own HandleInput
/// function, that should be able to handle codes specific to that application.\n
/// The appropriate version of the HandleInput method (TMPWorker's or the
/// overriding version) is automatically called depending on the message code.
void TMPWorker::HandleInput(MPCodeBufPair &msg)
{
   unsigned code = msg.first;

   std::string reply = fId;
   if (code == MPCode::kMessage) {
      //general message, ignore it
      reply += ": ok";
      MPSend(fS.get(), MPCode::kMessage, reply.c_str());
   } else if (code == MPCode::kError) {
      //general error, ignore it
      reply += ": ko";
      MPSend(fS.get(), MPCode::kMessage, reply.c_str());
   } else if (code == MPCode::kShutdownOrder || code == MPCode::kFatalError) {
      //client is asking the server to shutdown or client is dying
      MPSend(fS.get(), MPCode::kShutdownNotice, reply.c_str());
      gSystem->Exit(0);
   } else {
      reply += ": unknown code received. code=" + std::to_string(code);
      MPSend(fS.get(), MPCode::kError, reply.c_str());
   }
}

//////////////////////////////////////////////////////////////////////////
/// Error sender

void TMPWorker::SendError(const std::string& errmsg, unsigned int errcode)
{
   std::string reply = fId + ": " + errmsg;
   MPSend(GetSocket(), errcode, reply.c_str());
}
