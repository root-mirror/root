// Author: Sergey Linev <S.Linev@gsi.de>
// Date: 2021-01-22
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RBrowserWidget
#define ROOT7_RBrowserWidget

#include <memory>
#include <map>
#include <string>

#include <ROOT/Browsable/RElement.hxx>

namespace ROOT {
namespace Experimental {

/** \class ROOT::Experimental::RBrowserWidget
\ingroup rbrowser
Abstract Web-based widget, which can be used in the RBrowser
Used to embed canvas, geometry viewer and potentially any other widgets
*/

class RBrowserWidget {

   std::string fName;   ///<!  widget name

   Browsable::RElementPath_t  fPath;  ///<! path of drawn element

public:

   explicit RBrowserWidget(const std::string &name) : fName(name) {};
   virtual ~RBrowserWidget() = default;

   virtual void Show(const std::string &) = 0;

   virtual void ResetConn() {}

   void SetPath(const Browsable::RElementPath_t &path) { fPath = path; }
   const Browsable::RElementPath_t &GetPath() const { return fPath; }

   const std::string &GetName() const { return fName; }
   virtual std::string GetKind() const = 0;
   virtual std::string GetUrl() = 0;
   virtual std::string GetTitle() { return ""; }

   virtual bool DrawElement(std::shared_ptr<Browsable::RElement> &, const std::string &) { return false; }
   virtual std::string SendWidgetContent() { return ""; }
};

class RBrowserWidgetProvider {
protected:
   using ProvidersMap_t = std::map<std::string, RBrowserWidgetProvider*>;

   virtual std::shared_ptr<RBrowserWidget> Create(const std::string &) = 0;

   static ProvidersMap_t& GetMap();

public:

   explicit RBrowserWidgetProvider(const std::string &kind);
   virtual ~RBrowserWidgetProvider();

   static std::shared_ptr<RBrowserWidget> CreateWidget(const std::string &kind, const std::string &name);
};


} // namespace Experimental
} // namespace ROOT

#endif
