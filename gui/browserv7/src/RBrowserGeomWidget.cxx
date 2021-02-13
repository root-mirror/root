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

#include "RBrowserWidget.hxx"

#include <ROOT/REveGeomViewer.hxx>

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"

using namespace ROOT::Experimental;

using namespace std::string_literals;


class RBrowserGeomWidget : public RBrowserWidget {
   REveGeomViewer fViewer;

   std::unique_ptr<Browsable::RHolder> fObject; // geometry object

   /** Create dummy geometry - when nothing else is there */
   TGeoManager *MakeDummy()
   {
      auto prev = gGeoManager;
      gGeoManager = nullptr;

      auto mgr = new TGeoManager("box", "poza1");
      TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
      TGeoMedium *med = new TGeoMedium("MED",1,mat);
      TGeoVolume *top = mgr->MakeBox("TOP",med,100,100,100);
      mgr->SetTopVolume(top);
      TGeoVolume *vol = mgr->MakeBox("BOX",med, 20,30,40);
      vol->SetLineColor(kRed);
      vol->SetLineWidth(2);
      top->AddNode(vol,1);
      mgr->CloseGeometry();

      gGeoManager = prev;
      return mgr;
   }

public:

   RBrowserGeomWidget(const std::string &name) : RBrowserWidget(name), fViewer()
   {
      fViewer.SetTitle(name);
      fViewer.SetShowHierarchy(false);

      // fViewer.SetGeometry(MakeDummy());
   }

   virtual ~RBrowserGeomWidget() = default;

   std::string GetKind() const override { return "geom"s; }

   void Show(const std::string &arg) override
   {
      fViewer.Show(arg);
   }

   std::string GetUrl() override
   {
      return "../"s + fViewer.GetWindowAddr() + "/"s;
   }

   bool DrawElement(std::shared_ptr<Browsable::RElement> &elem, const std::string &) override
   {
      if (!elem->IsCapable(Browsable::RElement::kActGeom))
         return false;

      fObject = elem->GetObject();
      if (!fObject)
         return false;

      auto vol = fObject->Get<TGeoVolume>();
      if (vol) {
         fViewer.SetGeometry(vol->GetGeoManager(), vol->GetName());
         return true;
      }

      // only handle TGeoManager now
      auto mgr = fObject->Get<TGeoManager>();
      if (!mgr) {
         fObject.release();
         return false;
      }

      fViewer.SetGeometry(const_cast<TGeoManager *>(mgr));

      return true;
   }

};

// ======================================================================

class RBrowserGeomProvider : public RBrowserWidgetProvider {
protected:
   std::shared_ptr<RBrowserWidget> Create(const std::string &name) final
   {
      return std::make_shared<RBrowserGeomWidget>(name);
   }
public:
   RBrowserGeomProvider() : RBrowserWidgetProvider("geom") {}
   ~RBrowserGeomProvider() = default;
} sRBrowserGeomProvider;
