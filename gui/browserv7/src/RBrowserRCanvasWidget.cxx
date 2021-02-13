// Author: Sergey Linev <S.Linev@gsi.de>
// Date: 2021-01-25
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "RBrowserWidget.hxx"

#include <ROOT/Browsable/RProvider.hxx>

#include <ROOT/RCanvas.hxx>

using namespace ROOT::Experimental;

using namespace std::string_literals;


class RBrowserRCanvasWidget : public RBrowserWidget {

   std::shared_ptr<RCanvas> fCanvas; ///<! drawn canvas

public:

   RBrowserRCanvasWidget(const std::string &name) : RBrowserWidget(name)
   {
      fCanvas = RCanvas::Create(name);
   }

   virtual ~RBrowserRCanvasWidget() = default;

   std::string GetKind() const override { return "rcanvas"s; }

   void Show(const std::string &arg) override
   {
      fCanvas->Show(arg);
   }

   std::string GetUrl() override
   {
      return "../"s + fCanvas->GetWindowAddr() + "/"s;
   }

   bool DrawElement(std::shared_ptr<Browsable::RElement> &elem, const std::string &opt) override
   {
      if (!elem->IsCapable(Browsable::RElement::kActDraw7))
         return false;

      auto obj = elem->GetObject();
      if (!obj)
         return false;

      std::shared_ptr<RPadBase> subpad = fCanvas;

      if (obj && Browsable::RProvider::Draw7(subpad, obj, opt)) {
         fCanvas->Modified();
         fCanvas->Update(true);
         return true;
      }

      return false;
   }

};

// ======================================================================

class RBrowserRCanvasProvider : public RBrowserWidgetProvider {
protected:
   std::shared_ptr<RBrowserWidget> Create(const std::string &name) final
   {
      return std::make_shared<RBrowserRCanvasWidget>(name);
   }
public:
   RBrowserRCanvasProvider() : RBrowserWidgetProvider("rcanvas") {}
   ~RBrowserRCanvasProvider() = default;
} sRBrowserRCanvasProvider;
