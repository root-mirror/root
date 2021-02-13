// Author: Sergey Linev <s.linev@gsi.de>
// Date: 2017-10-16
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RWebWindowsManager
#define ROOT7_RWebWindowsManager

#include <memory>
#include <string>
#include <thread>
#include <mutex>

#include "THttpEngine.h"

#include <ROOT/RWebWindow.hxx>

class THttpServer;
class THttpWSHandler;

namespace ROOT {
namespace Experimental {

class RWebWindowsManager {

   friend class RWebWindow;

private:
   std::unique_ptr<THttpServer> fServer; ///<! central communication with the all used displays
   std::string fAddr;                    ///<! HTTP address of the server
   std::recursive_mutex fMutex;          ///<! main mutex, used for window creations
   unsigned fIdCnt{0};                   ///<! counter for identifiers
   bool fUseHttpThrd{false};             ///<! use special thread for THttpServer
   bool fUseSenderThreads{false};        ///<! use extra threads for sending data from RWebWindow to clients
   float fLaunchTmout{30.};              ///<! timeout in seconds to start browser process, default 30s

   /// Returns true if http server use special thread for requests processing (default off)
   bool IsUseHttpThread() const { return fUseHttpThrd; }

   /// Returns true if extra threads to send data via websockets will be used (default off)
   bool IsUseSenderThreads() const { return fUseSenderThreads; }

   /// Returns timeout for launching new browser process
   float GetLaunchTmout() const { return fLaunchTmout; }

   void Unregister(RWebWindow &win);

   /// Show window in specified location, see Show() method for more details
   unsigned ShowWindow(RWebWindow &win, bool batch_mode, const RWebDisplayArgs &args);

   int WaitFor(RWebWindow &win, WebWindowWaitFunc_t check, bool timed = false, double tm = -1);

   static bool IsMainThrd();

   std::string GetUrl(const RWebWindow &win, bool remote = false);

   bool CreateServer(bool with_http = false);

public:
   RWebWindowsManager();

   ~RWebWindowsManager();

   /// Returns THttpServer instance
   THttpServer *GetServer() const { return fServer.get(); }

   /// Returns http address of the server, empty string when not available
   std::string GetServerAddr() const { return fAddr; }

   static std::shared_ptr<RWebWindowsManager> &Instance();

   std::shared_ptr<RWebWindow> CreateWindow();

   void Terminate();
};

} // namespace Experimental
} // namespace ROOT

#endif
