#include <ROOT/REveElement.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveViewer.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveCaloData.hxx>
#include <ROOT/REveCalo.hxx>
#include <ROOT/REveJetCone.hxx>
#include <ROOT/REveGeoShape.hxx>

#include "TFile.h"
#include "TGeoTube.h"

const char* histFile =
   "http://amraktad.web.cern.ch/amraktad/cms_calo_hist.root";


const Double_t kR_min = 128;
const Double_t kR_max = 129;
const Double_t kZ_d   = 268.39;

using namespace ROOT::Experimental;


void makeCalo2D(REveCalo3D* calo3d, const char* pname, REveProjection::EPType_e t)
{
   auto eveMng = ROOT::Experimental::gEve;

   auto geomScene  = eveMng->SpawnNewScene(Form("%s Geometry", pname), pname);
   auto eventScene = eveMng->SpawnNewScene(Form("%s Event Data", pname), pname);

   auto mng = new REveProjectionManager();
   mng->SetProjection(t);

   // project Calo3D
   REveCalo2D* calo2d = (REveCalo2D*) mng->ImportElements(calo3d, eventScene);

   // project reference geometry
   REveElement* geo2d =  mng->ImportElements(eveMng->GetGlobalScene()->FirstChild(), geomScene);

   auto view = eveMng->SpawnNewViewer("RPhi View", "");
   view->AddScene(geomScene);
   view->AddScene(eventScene);
}

void add_jet(REveElement* parent, const char* name,
             Float_t eta, Float_t phi,
             Float_t deta, Float_t dphi)
{
   auto jet = new REveJetCone(name, name);
   jet->SetMainTransparency(60);
   jet->SetLineColor(kRed);
   jet->SetCylinder(129 - 10, 268.36 - 10);
   jet->AddEllipticCone(eta, phi, deta, dphi);
   jet->SetPickable(kTRUE);
   jet->SetHighlightFrame(kFALSE);
   parent->AddElement(jet);
}

void calorimeters()
{
   auto eveMng = REveManager::Create();

   TFile::SetCacheFileDir(".");
   auto hf = TFile::Open(histFile, "CACHEREAD");
   auto ecalHist = (TH2F*)hf->Get("ecalLego");
   auto hcalHist = (TH2F*)hf->Get("hcalLego");
   auto data = new REveCaloDataHist();
   data->AddHistogram(ecalHist);
   data->RefSliceInfo(0).Setup("ECAL", 0.f, kBlue);
   data->AddHistogram(hcalHist);
   data->RefSliceInfo(1).Setup("HCAL", 0.1, kRed);
   eveMng->GetEventScene()->AddElement(data);

   auto b1 = new REveGeoShape("Barrel 1");
   b1->SetShape(new TGeoTube(kR_min, kR_max, kZ_d));
   b1->SetMainColor(kCyan);
   eveMng->GetGlobalScene()->AddElement(b1);

   auto calo3d = new REveCalo3D(data);
   calo3d->SetBarrelRadius(kR_max);
   calo3d->SetEndCapPos(kZ_d);
   calo3d->SetMaxTowerH(300);
   eveMng->GetEventScene()->AddElement(calo3d);

   add_jet(calo3d, "JetCone Lojz",  1.4,  1.0, 0.4, 0.2);
   add_jet(calo3d, "JetCone Mici", -2.0, -2.1, 0.2, 0.4);

   makeCalo2D(calo3d, "RPhi", REveProjection::kPT_RPhi);
   makeCalo2D(calo3d, "RhoZ", REveProjection::kPT_RhoZ);

   eveMng->Show();
}


