// Author: Sergey Linev <s.linev@gsi.de>
// Date: 2018-10-24
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RWebDisplayArgs
#define ROOT7_RWebDisplayArgs

#include <string>
#include <memory>

class THttpServer;

namespace ROOT {
namespace Experimental {

class RLogChannel;
/// Log channel for WebGUI diagnostics.
RLogChannel &WebGUILog();

class RWebWindow;

class RWebDisplayArgs {

friend class RWebWindow;

public:
   enum EBrowserKind {
      kChrome,   ///< Google Chrome browser
      kFirefox,  ///< Mozilla Firefox browser
      kNative,   ///< either Chrome or Firefox - both support major functionality
      kCEF,      ///< Chromium Embedded Framework - local display with CEF libs
      kQt5,      ///< QWebEngine libraries - Chrome code packed in qt5
      kLocal,    ///< either CEF or Qt5 - both runs on local display without real http server
      kStandard, ///< standard system web browser, not recognized by ROOT, without batch mode
      kEmbedded, ///< window will be embedded into other, no extra browser need to be started
      kCustom    ///< custom web browser, execution string should be provided
   };

protected:
   EBrowserKind fKind{kNative};   ///<! id of web browser used for display
   std::string fUrl;              ///<! URL to display
   std::string fExtraArgs;        ///<! extra arguments which will be append to exec string
   std::string fPageContent;      ///<! HTML page content
   std::string fRedirectOutput;   ///<! filename where browser output should be redirected
   bool fHeadless{false};         ///<! is browser runs in headless mode
   bool fStandalone{true};        ///<! indicates if browser should run isolated from other browser instances
   THttpServer *fServer{nullptr}; ///<! http server which handle all requests
   int fWidth{0};                 ///<! custom window width, when not specified - used RWebWindow geometry
   int fHeight{0};                ///<! custom window height, when not specified - used RWebWindow geometry
   int fX{-1};                    ///<! custom window x position, negative is default
   int fY{-1};                    ///<! custom window y position, negative is default
   std::string fUrlOpt;           ///<! extra URL options, which are append to window URL
   std::string fExec;             ///<! string to run browser, used with kCustom type
   void *fDriverData{nullptr};    ///<! special data delivered to driver, can be used for QWebEngine

   std::shared_ptr<RWebWindow> fMaster; ///<!  master window
   int fMasterChannel{-1};              ///<!  used master channel

   bool SetSizeAsStr(const std::string &str);
   bool SetPosAsStr(const std::string &str);

public:
   RWebDisplayArgs();

   RWebDisplayArgs(const std::string &browser);

   RWebDisplayArgs(const char *browser);

   RWebDisplayArgs(int width, int height, int x = -1, int y = -1, const std::string &browser = "");

   RWebDisplayArgs(std::shared_ptr<RWebWindow> master, int channel = -1);

   virtual ~RWebDisplayArgs();

   RWebDisplayArgs &SetBrowserKind(const std::string &kind);
   /// set browser kind, see EBrowserKind for allowed values
   RWebDisplayArgs &SetBrowserKind(EBrowserKind kind) { fKind = kind; return *this; }
   /// returns configured browser kind, see EBrowserKind for supported values
   EBrowserKind GetBrowserKind() const { return fKind; }
   std::string GetBrowserName() const;

   void SetMasterWindow(std::shared_ptr<RWebWindow> master, int channel = -1);

   /// returns true if local display like CEF or Qt5 QWebEngine should be used
   bool IsLocalDisplay() const
   {
      return (GetBrowserKind() == kLocal) || (GetBrowserKind() == kCEF) || (GetBrowserKind() == kQt5);
   }

   /// returns true if browser supports headless (batch) mode, used for image production
   bool IsSupportHeadless() const
   {
      return (GetBrowserKind() == kNative) || (GetBrowserKind() == kChrome) || (GetBrowserKind() == kCEF) || (GetBrowserKind() == kQt5);
   }

   /// set window url
   RWebDisplayArgs &SetUrl(const std::string &url) { fUrl = url; return *this; }
   /// returns window url
   const std::string &GetUrl() const { return fUrl; }

   /// set window url
   RWebDisplayArgs &SetPageContent(const std::string &cont) { fPageContent = cont; return *this; }
   /// returns window url
   const std::string &GetPageContent() const { return fPageContent; }

   /// Set standalone mode for running browser, default on
   /// When disabled, normal browser window (or just tab) will be started
   void SetStandalone(bool on = true) { fStandalone = on; }
   /// Return true if browser should runs in standalone mode
   bool IsStandalone() const { return fStandalone; }

   /// set window url options
   RWebDisplayArgs &SetUrlOpt(const std::string &opt) { fUrlOpt = opt; return *this; }
   /// returns window url options
   const std::string &GetUrlOpt() const { return fUrlOpt; }

   /// append extra url options, add "&" as separator if required
   void AppendUrlOpt(const std::string &opt);

   /// returns window url with append options
   std::string GetFullUrl() const;

   /// set headless mode
   void SetHeadless(bool on = true) { fHeadless = on; }
   /// returns headless mode
   bool IsHeadless() const { return fHeadless; }

   /// set preferable web window width
   RWebDisplayArgs &SetWidth(int w = 0) { fWidth = w; return *this; }
   /// set preferable web window height
   RWebDisplayArgs &SetHeight(int h = 0) { fHeight = h; return *this; }
   RWebDisplayArgs &SetSize(int w, int h) { fWidth = w; fHeight = h; return *this; }

   /// set preferable web window x position, negative is default
   RWebDisplayArgs &SetX(int x = -1) { fX = x; return *this; }
   /// set preferable web window y position, negative is default
   RWebDisplayArgs &SetY(int y = -1) { fY = y; return *this; }
   RWebDisplayArgs &SetPos(int x = -1, int y = -1) { fX = x; fY = y; return *this; }

   /// returns preferable web window width
   int GetWidth() const { return fWidth; }
   /// returns preferable web window height
   int GetHeight() const { return fHeight; }
   /// set preferable web window x position
   int GetX() const { return fX; }
   /// set preferable web window y position
   int GetY() const { return fY; }

   void SetExtraArgs(const std::string &args) { fExtraArgs = args; }
   const std::string &GetExtraArgs() const { return fExtraArgs; }

   void SetRedirectOutput(const std::string &fname = "") { fRedirectOutput = fname; }
   const std::string &GetRedirectOutput() const { return fRedirectOutput; }

   /// set custom executable to start web browser
   void SetCustomExec(const std::string &exec);
   /// returns custom executable to start web browser
   std::string GetCustomExec() const;

   /// set http server instance, used for window display
   void SetHttpServer(THttpServer *serv) { fServer = serv; }
   /// returns http server instance, used for window display
   THttpServer *GetHttpServer() const { return fServer; }

   /// [internal] set web-driver data, used to start window
   void SetDriverData(void *data) { fDriverData = data; }
   /// [internal] returns web-driver data, used to start window
   void *GetDriverData() const { return fDriverData; }

   static std::string GetQt5EmbedQualifier(const void *qparent, const std::string &urlopt = "");
};

}
}

#endif
