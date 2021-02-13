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

#include "TCanvas.h"
#include "TWebCanvas.h"

using namespace ROOT::Experimental;

using namespace std::string_literals;


class RBrowserTCanvasWidget : public RBrowserWidget {

   std::unique_ptr<TCanvas> fCanvas; ///<! drawn canvas
   TWebCanvas *fWebCanvas{nullptr};  ///<! web implementation, owned by TCanvas

   std::unique_ptr<Browsable::RHolder> fObject; // TObject drawing in the TCanvas

public:

   RBrowserTCanvasWidget(const std::string &name) : RBrowserWidget(name)
   {
      fCanvas = std::make_unique<TCanvas>(kFALSE);
      fCanvas->SetName(name.c_str());
      fCanvas->SetTitle(name.c_str());
      fCanvas->ResetBit(TCanvas::kShowEditor);
      fCanvas->ResetBit(TCanvas::kShowToolBar);
      fCanvas->SetCanvas(fCanvas.get());
      fCanvas->SetBatch(kTRUE); // mark canvas as batch
      fCanvas->SetEditable(kTRUE); // ensure fPrimitives are created

      // create implementation
      fWebCanvas = new TWebCanvas(fCanvas.get(), "title", 0, 0, 800, 600);

      // assign implementation
      fCanvas->SetCanvasImp(fWebCanvas);
   }

   virtual ~RBrowserTCanvasWidget() = default;

   std::string GetKind() const override { return "tcanvas"s; }

   void Show(const std::string &arg) override
   {
      fWebCanvas->ShowWebWindow(arg);
   }

   std::string GetUrl() override
   {
      return "../"s + fWebCanvas->GetWebWindow()->GetAddr() + "/"s;
   }

   bool DrawElement(std::shared_ptr<Browsable::RElement> &elem, const std::string &opt) override
   {
      if (!elem->IsCapable(Browsable::RElement::kActDraw6))
         return false;

      fObject = elem->GetObject();

      if (fObject && Browsable::RProvider::Draw6(fCanvas.get(), fObject, opt)) {
         fCanvas->ForceUpdate();
         return true;
      }

      return false;
   }

};

// ======================================================================

class RBrowserTCanvasProvider : public RBrowserWidgetProvider {
protected:
   std::shared_ptr<RBrowserWidget> Create(const std::string &name) final
   {
      return std::make_shared<RBrowserTCanvasWidget>(name);
   }
public:
   RBrowserTCanvasProvider() : RBrowserWidgetProvider("tcanvas") {}
   ~RBrowserTCanvasProvider() = default;
} sRBrowserTCanvasProvider;

