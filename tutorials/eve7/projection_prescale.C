/// \file
/// \ingroup tutorial_eve7
///  This example display projection prescale
///
/// \macro_code
///


#include <sstream>
#include <iostream>

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TSystem.h"
#include "TFile.h"
#include "TRandom.h"

#include <ROOT/REveGeoShape.hxx>
#include <ROOT/REveTrans.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveViewer.hxx>
#include <ROOT/REveElement.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REvePointSet.hxx>
#include <ROOT/REveLine.hxx>
#include <ROOT/REveProjectionManager.hxx>
#include <ROOT/REveProjectionBases.hxx>

namespace REX = ROOT::Experimental;

REX::REveManager *eveMng = nullptr;

REX::REvePointSet* getPointSet(int npoints = 2, float s=2, int color=28)
{
   TRandom &r = *gRandom;

   auto ps = new REX::REvePointSet("testPnts", "title", npoints);

   for (Int_t i=0; i<npoints; ++i)
       ps->SetNextPoint(r.Uniform(-s,s), r.Uniform(-s,s), r.Uniform(-s,s));

   ps->SetMarkerColor(color);
   ps->SetMarkerSize(3+r.Uniform(1, 2));
   // ps->SetMarkerStyle(4);
   return ps;
}

void makeProjectedViewsAndScene(REX::REveProjection::EPType_e type, bool scale)
{
   auto rPhiGeomScene  = eveMng->SpawnNewScene(Form("Project%s Geo", scale ? "PreScaled" : ""));
   auto rPhiEventScene = eveMng->SpawnNewScene(Form("Project%s Event", scale ? "PreScaled" : ""));

   auto mngRhoPhi = new REX::REveProjectionManager(type);
   if (scale) {
      REX::REveProjection* p = mngRhoPhi->GetProjection();
      p->AddPreScaleEntry(0, 0,   4);    // r scale 4 from 0
      p->AddPreScaleEntry(0, 45,  1);    // r scale 1 from 45
      p->AddPreScaleEntry(0, 310, 0.5);
      p->SetUsePreScale(kTRUE);
   }
   auto rphiView = eveMng->SpawnNewViewer("Projected View", "");
   rphiView->AddScene(rPhiGeomScene);
   rphiView->AddScene(rPhiEventScene);

   for (auto &ie : eveMng->GetGlobalScene()->RefChildren())
      mngRhoPhi->ImportElements(ie, rPhiGeomScene);

   for (auto &ie : eveMng->GetEventScene()->RefChildren())
      mngRhoPhi->ImportElements(ie, rPhiEventScene);
}

TGeoNode* getNodeFromPath( TGeoNode* top, std::string path)
{
   TGeoNode* node = top;
   istringstream f(path);
   string s;
   while (getline(f, s, '/'))
      node = node->GetVolume()->FindNode(s.c_str());

   return node;
}


void projection_prescale(std::string type = "RhPhi")
{
   eveMng = REX::REveManager::Create();

   // static scene
   TFile::SetCacheFileDir(".");
   auto geoManager = eveMng->GetGeometry("http://root.cern/files/cms.root");
   TGeoNode* top = geoManager->GetTopVolume()->FindNode("CMSE_1");
   auto holder = new REX::REveElement("MUON");
   eveMng->GetGlobalScene()->AddElement(holder);
   auto n = getNodeFromPath(top, "MUON_1/MB_1");
   auto m = new REX::REveGeoShape("MB_1");
   m->SetShape(n->GetVolume()->GetShape());
   m->SetMainColor(kOrange);
   holder->AddElement(m);

   auto bv =  n->GetVolume();
   for (int i = 1; i < 5; ++i ) {

      auto n = bv->FindNode(Form("MBXC_%d",i));
      auto gss = n->GetVolume()->GetShape();
      auto b1s = new REX::REveGeoShape(Form("Arc %d", i));
      b1s->InitMainTrans();
      const double* move = n->GetMatrix()->GetTranslation();
      b1s->RefMainTrans().SetFrom( *(n->GetMatrix()));
      b1s->SetShape(gss);
      b1s->SetMainColor(kBlue);
      holder->AddElement(b1s);
   }

   // event scene
   auto line = new REX::REveLine();
   line->SetNextPoint(0, 0, 0);
   float a = 300;
   line->SetNextPoint(a, a, a);
   eveMng->GetEventScene()->AddElement(line);

   auto line2 = new REX::REveLine();
   line2->SetNextPoint(0, 0, 0);
   float b = 30;
   line2->SetNextPoint(b, b+5, b);
   line2->SetMainColor(kRed);
   eveMng->GetEventScene()->AddElement(line2);

   auto points = getPointSet(10, 30);
   eveMng->GetEventScene()->AddElement(points);

   // make scaled and plain projected views
   if (type == "RPhi") {
      makeProjectedViewsAndScene(REX::REveProjection::kPT_RPhi, true);
      makeProjectedViewsAndScene(REX::REveProjection::kPT_RPhi, false);
   }
   else {
      makeProjectedViewsAndScene(REX::REveProjection::kPT_RhoZ, true);
      makeProjectedViewsAndScene(REX::REveProjection::kPT_RhoZ, false);
   }

   eveMng->Show();
}

