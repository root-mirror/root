/* @(#)root/multiproc:$Id$ */
// Author: Enrico Guiraud July 2015

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "MPCode.h"
#include "TGuiFactory.h" //gGuiFactory
#include "TError.h" //gErrorIgnoreLevel
#include "TMPClient.h"
#include "TMPWorker.h"
#include "TROOT.h" //gROOT
#include "TSocket.h"
#include "TSystem.h" //gSystem
#include "TVirtualX.h" //gVirtualX
#include <errno.h> //errno, used by socketpair
#include <memory> //unique_ptr
#include <sys/socket.h> //socketpair
#include <sys/wait.h> // waitpid
#include <unistd.h> // close, fork
#include <dlfcn.h>

//////////////////////////////////////////////////////////////////////////
///
/// \class TMPClient
///
/// Base class for multiprocess applications' clients. It provides a
/// simple interface to fork a ROOT session into server/worker sessions
/// and exchange messages with them. Multiprocessing applications can build
/// on TMPClient and TMPWorker: the class providing multiprocess
/// functionalities to users should inherit (possibly privately) from
/// TMPClient, and the workers executing tasks should inherit from TMPWorker.
///
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
/// Class constructor.
/// \param nWorkers
/// \parblock
/// the number of children processes that will be created by
/// Fork, i.e. the number of workers that will be available after this call.
/// The default value (0) means that a number of workers equal to the number
/// of cores of the machine is going to be spawned. If that information is
/// not available, 2 workers are created instead.
/// \endparblock
TMPClient::TMPClient(unsigned nWorkers) : fIsParent(true), fWorkerPids(), fMon(), fNWorkers(0)
{
   // decide on number of workers
   if (nWorkers) {
      fNWorkers = nWorkers;
   } else {
      SysInfo_t si;
      if (gSystem->GetSysInfo(&si) == 0)
         fNWorkers = si.fCpus;
      else
         fNWorkers = 2;
   }
}


//////////////////////////////////////////////////////////////////////////
/// Class destructor.
/// This method is in charge of shutting down any remaining worker,
/// closing off connections and reap the terminated children processes.
TMPClient::~TMPClient()
{
   Broadcast(MPCode::kShutdownOrder);
   TList *l = fMon.GetListOfActives();
   l->Delete();
   delete l;
   l = fMon.GetListOfDeActives();
   l->Delete();
   delete l;
   fMon.RemoveAll();
   ReapWorkers();
}

namespace ROOT {
   namespace Internal {
      /// Class to acquire and release the Python GIL where it applies, i.e.
      /// if libPython is loaded and the interpreter is initialized.
      class TGILRAII {
         using Py_IsInitialized_type = int (*)(void);
         using PyGILState_Ensure_type = void* (*)(void);
         using PyGILState_Release_type = void (*)(void*);
         void* fPyGILState_STATE = nullptr;
         template<class FPTYPE>
         FPTYPE GetSymT(const char* name) {return (FPTYPE) dlsym(nullptr,name);}
      public:
         TGILRAII()
         {
            auto Py_IsInitialized = GetSymT<Py_IsInitialized_type>("Py_IsInitialized");
            if (!Py_IsInitialized || !Py_IsInitialized()) return;
            auto PyGILState_Ensure = GetSymT<PyGILState_Ensure_type>("PyGILState_Ensure");
            if (PyGILState_Ensure) fPyGILState_STATE = PyGILState_Ensure();
         }

         ~TGILRAII()
         {
            auto PyGILState_Release = GetSymT<PyGILState_Release_type>("PyGILState_Release");
            if (fPyGILState_STATE && PyGILState_Release) PyGILState_Release(fPyGILState_STATE);
         }
      };
   }
}

//////////////////////////////////////////////////////////////////////////
/// This method forks the ROOT session into fNWorkers children processes.
/// The ROOT sessions spawned in this way will not have graphical
/// capabilities and will not read from standard input, but will be
/// connected to the original (interactive) session through TSockets.
/// The children processes' PIDs are added to the fWorkerPids vector.
/// The parent session can then communicate with the children using the
/// Broadcast and MPSend methods, and receive messages through MPRecv.\n
/// \param server
/// \parblock
/// A pointer to an instance of the class that will take control
/// of the subprocesses after forking. Applications should implement their
/// own class inheriting from TMPWorker. Behaviour can be customized
/// overriding TMPWorker::HandleInput.
/// \endparblock
/// \return true if Fork succeeded, false otherwise
bool TMPClient::Fork(TMPWorker &server)
{
   std::string basePath = "/tmp/ROOTMP-";

   //fork as many times as needed and save pids
   pid_t pid = 1; //must be positive to handle the case in which fNWorkers is 0
   int sockets[2]; //sockets file descriptors
   unsigned nWorker = 0;
   for (; nWorker < fNWorkers; ++nWorker) {
      //create socket pair
      int ret = socketpair(AF_UNIX, SOCK_STREAM, 0, sockets);
      if (ret != 0) {
         Error("TMPClient::Fork", "[E][C] Could not create socketpair. Error n. . Now retrying.\n%d", errno);
         --nWorker;
         continue;
      }

      //fork
      {
         ROOT::Internal::TGILRAII tgilraai;
         pid = fork();
      }

      if (!pid) {
         //child process, exit loop. sockets[1] is the fd that should be used
         break;
      } else {
         //parent process, create TSocket with current value of sockets[0]
         close(sockets[1]); //we don't need this
         TSocket *s = new TSocket(sockets[0], (std::to_string(pid)).c_str()); //TSocket's constructor with this signature seems much faster than TSocket(int fd)
         if (s && s->IsValid()) {
            fMon.Add(s);
            fWorkerPids.push_back(pid);
         } else {
            Error("TMPClient::Fork","[E][C] Could not connect to worker with pid %d. Giving up.\n", pid);
            delete s;
         }
      }
   }

   if (pid) {
      //parent returns here
      return true;
   } else {
      //CHILD/WORKER
      fIsParent = false;
      close(sockets[0]); //we don't need this

      //override signal handler (make the servers exit on SIGINT)
      TSeqCollection *signalHandlers = gSystem->GetListOfSignalHandlers();
      TSignalHandler *sh = nullptr;
      if (signalHandlers && signalHandlers->GetSize() > 0)
         sh = (TSignalHandler *)signalHandlers->First();
      if (sh)
         gSystem->RemoveSignalHandler(sh);

      //remove stdin from eventloop and close it
      TSeqCollection *fileHandlers = gSystem->GetListOfFileHandlers();
      if (fileHandlers) {
         for (auto h : *fileHandlers) {
            if (h && ((TFileHandler *)h)->GetFd() == 0) {
               gSystem->RemoveFileHandler((TFileHandler *)h);
               break;
            }
         }
      }
      close(0);

      std::unique_ptr<TList> lofact(fMon.GetListOfActives());
      while (lofact && (lofact->GetSize() > 0)) {
         TSocket *s = (TSocket *) lofact->First();
         lofact.reset(nullptr); // delete list before object is destroyed
         fMon.Remove(s);
         delete s;
         lofact.reset(fMon.GetListOfActives());
      }

      std::unique_ptr<TList> lofdeact(fMon.GetListOfDeActives());
      while (lofdeact && (lofdeact->GetSize() > 0)) {
         TSocket *s = (TSocket *) lofdeact->First();
         lofdeact.reset(nullptr); // delete list before object is destroyed
         fMon.Remove(s);
         delete s;
         lofdeact.reset(fMon.GetListOfDeActives());
      }

      //disable graphics
      //these instructions were copied from TApplication::MakeBatch
      gROOT->SetBatch();
      if (gGuiFactory != gBatchGuiFactory)
         delete gGuiFactory;
      gGuiFactory = gBatchGuiFactory;
#ifndef R__WIN32
      if (gVirtualX != gGXBatch)
         delete gVirtualX;
#endif
      gVirtualX = gGXBatch;

      //prepare server and add it to eventloop
      server.Init(sockets[1], nWorker);

      //enter worker loop
      server.Run();
   }

   //control should never reach here
   return true;
}


//////////////////////////////////////////////////////////////////////////
/// Send a message with the specified code to at most nMessages workers.
/// Sockets can either be in an "active" or "non-active" state. This method
/// activates all the sockets through which the client is connected to the
/// workers, and deactivates them when a message is sent to the corresponding
/// worker. This way the sockets pertaining to workers who have been left
/// idle will be the only ones in the active list
/// (TSocket::GetMonitor()->GetListOfActives()) after execution.
/// \param code the code to send (e.g. EMPCode)
/// \param nMessages
/// \parblock
/// the maximum number of messages to send.
/// If `nMessages == 0 || nMessage > fNWorkers`, send a message to every worker.
/// \endparblock
/// \return the number of messages successfully sent
unsigned TMPClient::Broadcast(unsigned code, unsigned nMessages)
{
   if (nMessages == 0)
      nMessages = fNWorkers;
   unsigned count = 0;
   fMon.ActivateAll();

   //send message to all sockets
   std::unique_ptr<TList> lp(fMon.GetListOfActives());
   for (auto s : *lp) {
      if (count == nMessages)
         break;
      if (MPSend((TSocket *)s, code)) {
         fMon.DeActivate((TSocket *)s);
         ++count;
      } else {
         Error("TMPClient:Broadcast", "[E] Could not send message to server\n");
      }
   }

   return count;
}


//////////////////////////////////////////////////////////////////////////
/// DeActivate a certain socket.
/// This does not remove it from the monitor: it will be reactivated by
/// the next call to Broadcast() (or possibly other methods that are
/// specified to do so).\n
/// A socket should be DeActivated when the corresponding
/// worker is done *for now* and we want to stop listening to this worker's
/// socket. If the worker is done *forever*, Remove() should be used instead.
/// \param s the socket to be deactivated
void TMPClient::DeActivate(TSocket *s)
{
   fMon.DeActivate(s);
}


//////////////////////////////////////////////////////////////////////////
/// Remove a certain socket from the monitor.
/// A socket should be Removed from the monitor when the
/// corresponding worker is done *forever*. For example HandleMPCode()
/// calls this method on sockets pertaining to workers which sent an
/// MPCode::kShutdownNotice.\n
/// If the worker is done *for now*, DeActivate should be used instead.
/// \param s the socket to be removed from the monitor fMon
void TMPClient::Remove(TSocket *s)
{
   fMon.Remove(s);
   delete s;
}


//////////////////////////////////////////////////////////////////////////
/// Wait on worker processes and remove their pids from fWorkerPids.
/// A blocking waitpid is called, but this should actually not block
/// execution since ReapWorkers should only be called when all workers
/// have already quit. ReapWorkers is then called not to leave zombie
/// processes hanging around, and to clean-up fWorkerPids.
void TMPClient::ReapWorkers()
{
   for (auto &pid : fWorkerPids) {
      waitpid(pid, nullptr, 0);
   }
   fWorkerPids.clear();
}


//////////////////////////////////////////////////////////////////////////
/// Handle messages containing an EMPCode.
/// This method should be called upon receiving a message with a code >= 1000
/// (i.e. EMPCode). It handles the most generic types of messages.\n
/// Classes inheriting from TMPClient should implement a similar method
/// to handle message codes specific to the application they're part of.\n
/// \param msg the MPCodeBufPair returned by a MPRecv call
/// \param s
/// \parblock
/// a pointer to the socket from which the message has been received is passed.
/// This way HandleMPCode knows which socket to reply on.
/// \endparblock
void TMPClient::HandleMPCode(MPCodeBufPair &msg, TSocket *s)
{
   unsigned code = msg.first;
   //message contains server's pid. retrieve it
   const char *str = ReadBuffer<const char*>(msg.second.get());

   if (code == MPCode::kMessage) {
      Error("TMPClient::HandleMPCode", "[I][C] message received: %s\n", str);
   } else if (code == MPCode::kError) {
      Error("TMPClient::HandleMPCode", "[E][C] error message received: %s\n", str);
   } else if (code == MPCode::kShutdownNotice || code == MPCode::kFatalError) {
      if (gDebug > 0) //generally users don't want to know this
         Error("TMPClient::HandleMPCode", "[I][C] shutdown notice received from %s\n", str);
      Remove(s);
   } else
       Error("TMPClient::HandleMPCode", "[W][C] unknown code received. code=%d\n", code);
   delete [] str;
}
