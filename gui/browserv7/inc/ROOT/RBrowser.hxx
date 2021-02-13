// Authors: Bertrand Bellenot <bertrand.bellenot@cern.ch> Sergey Linev <S.Linev@gsi.de>
// Date: 2019-02-28
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RBrowser
#define ROOT7_RBrowser

#include <ROOT/RWebWindow.hxx>
#include <ROOT/RBrowserData.hxx>

#include <vector>
#include <memory>

namespace ROOT {
namespace Experimental {

class RBrowserWidget;

class RBrowser {

protected:

   std::string fTitle;  ///<! title
   unsigned fConnId{0}; ///<! default connection id

   bool fUseRCanvas{false};             ///<!  which canvas should be used
   std::string fActiveWidgetName;        ///<! name of active widget
   std::vector<std::shared_ptr<RBrowserWidget>> fWidgets; ///<!  all browser widgets
   int fWidgetCnt{0};                                     ///<! counter for created widgets

   std::shared_ptr<RWebWindow> fWebWindow;   ///<! web window to browser

   RBrowserData  fBrowsable;                   ///<! central browsing element

   std::shared_ptr<RBrowserWidget> AddWidget(const std::string &kind);
   std::shared_ptr<RBrowserWidget> FindWidget(const std::string &name) const;
   std::shared_ptr<RBrowserWidget> GetActiveWidget() const { return FindWidget(fActiveWidgetName); }

   void CloseTab(const std::string &name);

   std::string ProcessBrowserRequest(const std::string &msg);
   std::string ProcessDblClick(std::vector<std::string> &args);
   std::string NewWidgetMsg(std::shared_ptr<RBrowserWidget> &widget);
   long ProcessRunMacro(const std::string &file_path);
   void ProcessSaveFile(const std::string &fname, const std::string &content);
   std::string GetCurrentWorkingDirectory();

   std::vector<std::string> GetRootHistory();
   std::vector<std::string> GetRootLogs();

   void SendInitMsg(unsigned connid);
   void ProcessMsg(unsigned connid, const std::string &arg);

public:
   RBrowser(bool use_rcanvas = true);
   virtual ~RBrowser();

   bool GetUseRCanvas() const { return fUseRCanvas; }
   void SetUseRCanvas(bool on = true) { fUseRCanvas = on; }

   /// show Browser in specified place
   void Show(const RWebDisplayArgs &args = "", bool always_start_new_browser = false);

   /// hide Browser
   void Hide();

};

} // namespace Experimental
} // namespace ROOT

#endif
