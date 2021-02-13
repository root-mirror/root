/// \file
/// \ingroup tutorial_eve7
///  This example display geometry, tracks and hits in web browser
///
/// \macro_code
///

#include <vector>
#include <string>
#include <iostream>

#include "TClass.h"
#include "TRandom.h"
#include "TGeoTube.h"
#include "TGeoSphere.h"
#include "TParticle.h"
#include "TApplication.h"
#include "TMatrixDSym.h"
#include "TVector.h"
#include "TMatrixDEigen.h"

#include <ROOT/REveGeoShape.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveViewer.hxx>
#include <ROOT/REveElement.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveUtil.hxx>
#include <ROOT/REveGeoShape.hxx>
#include <ROOT/REveProjectionManager.hxx>
#include <ROOT/REveProjectionBases.hxx>
#include <ROOT/REvePointSet.hxx>
#include <ROOT/REveJetCone.hxx>
#include <ROOT/REveTrans.hxx>

#include <ROOT/REveTrack.hxx>
#include <ROOT/REveTrackPropagator.hxx>

namespace REX = ROOT::Experimental;

// globals
REX::REveManager *eveMng = nullptr;
REX::REveProjectionManager *mngRhoPhi = nullptr;
REX::REveProjectionManager *mngRhoZ   = nullptr;
REX::REveScene  *rPhiGeomScene = nullptr, *rPhiEventScene = nullptr;
REX::REveScene  *rhoZGeomScene = nullptr, *rhoZEventScene = nullptr;
REX::REveViewer *rphiView = nullptr;
REX::REveViewer *rhoZView = nullptr;

const Double_t kR_min = 240;
const Double_t kR_max = 250;
const Double_t kZ_d   = 300;


REX::REvePointSet *getPointSet(int npoints = 2, float s=2, int color=28)
{
   TRandom &r = *gRandom;

   auto ps = new REX::REvePointSet("fu", "", npoints);

   for (Int_t i=0; i<npoints; ++i)
       ps->SetNextPoint(r.Uniform(-s,s), r.Uniform(-s,s), r.Uniform(-s,s));

   ps->SetMarkerColor(color);
   ps->SetMarkerSize(3+r.Uniform(1, 7));
   ps->SetMarkerStyle(4);
   return ps;
}

void addPoints()
{
   REX::REveElement* event = eveMng->GetEventScene();

   auto pntHolder = new REX::REveElement("Hits");

   auto ps1 = getPointSet(20, 100);
   ps1->SetName("Points_1");
   ps1->SetTitle("Points_1 title"); // used as tooltip

   pntHolder->AddElement(ps1);

   auto ps2 = getPointSet(10, 200, 4);
   ps2->SetName("Points_2");
   ps2->SetTitle("Points_2 title"); // used as tooltip
   pntHolder->AddElement(ps2);

   event->AddElement(pntHolder);
}

void addTracks()
{
   TRandom &r = *gRandom;

   REX::REveElement* event = eveMng->GetEventScene();
   auto prop = new REX::REveTrackPropagator();
   prop->SetMagFieldObj(new REX::REveMagFieldDuo(350, 3.5, -2.0));
   prop->SetMaxR(300);
   prop->SetMaxZ(600);
   prop->SetMaxOrbs(6);

   auto trackHolder = new REX::REveElement("Tracks");

   double v = 0.2;
   double m = 5;

   int N_Tracks = 10 + r.Integer(20);
   for (int i = 0; i < N_Tracks; i++)
   {
      TParticle* p = new TParticle();

      int pdg = 11 * (r.Integer(2) > 0 ? 1 : -1);
      p->SetPdgCode(pdg);

      p->SetProductionVertex(r.Uniform(-v,v), r.Uniform(-v,v), r.Uniform(-v,v), 1);
      p->SetMomentum(r.Uniform(-m,m), r.Uniform(-m,m), r.Uniform(-m,m)*r.Uniform(1, 3), 1);
      auto track = new REX::REveTrack(p, 1, prop);
      track->MakeTrack();
      if (i % 4 == 3) track->SetLineStyle(2); // enabled dashed style for some tracks
      track->SetMainColor(kBlue);
      track->SetName(Form("RandomTrack_%d", i));
      track->SetTitle(Form("RandomTrack_%d title", i)); // used as tooltip
      trackHolder->AddElement(track);
   }

   event->AddElement(trackHolder);
}

void addJets()
{
   TRandom &r = *gRandom;

   REX::REveElement *event = eveMng->GetEventScene();
   auto jetHolder = new REX::REveElement("Jets");

   int N_Jets = 5 + r.Integer(5);
   for (int i = 0; i < N_Jets; i++)
   {
      auto jet = new REX::REveJetCone(Form("Jet_%d", i));
      jet->SetTitle(Form("Jet_%d title", i)); // used as tooltip
      jet->SetCylinder(2*kR_max, 2*kZ_d);
      jet->AddEllipticCone(r.Uniform(-3.5, 3.5), r.Uniform(0, TMath::TwoPi()),
                           r.Uniform(0.02, 0.2), r.Uniform(0.02, 0.3));
      jet->SetFillColor(kPink - 8);
      jet->SetLineColor(kViolet - 7);

      jetHolder->AddElement(jet);
   }
   event->AddElement(jetHolder);
}

void makeEventScene()
{
   addPoints();
   addTracks();
   addJets();
}

void makeGeometryScene()
{
   auto b1 = new REX::REveGeoShape("Barrel 1");
   b1->SetShape(new TGeoTube(kR_min, kR_max, kZ_d));
   b1->SetMainColor(kCyan);
   eveMng->GetGlobalScene()->AddElement(b1);

   // Debug of surface fill in RPhi (index buffer screwed).
   // b1->SetNSegments(3);
   b1->SetNSegments(40);
}


void createProjectionStuff()
{
   // project RhoPhi
   rPhiGeomScene  = eveMng->SpawnNewScene("RPhi Geometry","RPhi");
   rPhiEventScene = eveMng->SpawnNewScene("RPhi Event Data","RPhi");

   mngRhoPhi = new REX::REveProjectionManager(REX::REveProjection::kPT_RPhi);

   rphiView = eveMng->SpawnNewViewer("RPhi View", "");
   rphiView->AddScene(rPhiGeomScene);
   rphiView->AddScene(rPhiEventScene);

   // ----------------------------------------------------------------

   rhoZGeomScene  = eveMng->SpawnNewScene("RhoZ Geometry", "RhoZ");
   rhoZEventScene = eveMng->SpawnNewScene("RhoZ Event Data","RhoZ");

   mngRhoZ = new REX::REveProjectionManager(REX::REveProjection::kPT_RhoZ);

   rhoZView = eveMng->SpawnNewViewer("RhoZ View", "");
   rhoZView->AddScene(rhoZGeomScene);
   rhoZView->AddScene(rhoZEventScene);
}

void projectScenes(bool geomp, bool eventp)
{
   if (geomp)
   {
      for (auto &ie : eveMng->GetGlobalScene()->RefChildren())
      {
         mngRhoPhi->ImportElements(ie, rPhiGeomScene);
         mngRhoZ  ->ImportElements(ie, rhoZGeomScene);
      }
   }
   if (eventp)
   {
      for (auto &ie : eveMng->GetEventScene()->RefChildren())
      {
         mngRhoPhi->ImportElements(ie, rPhiEventScene);
         mngRhoZ  ->ImportElements(ie, rhoZEventScene);
      }
   }

   // auto t0 = eveMng->GetEventScene()->FindChild("Tracks")->FirstChild();
   // printf("t0=%p, %s %s\n", t0, t0->GetElementName(), t0->IsA()->GetName());
   // dynamic_cast<REX::REveTrack*>(t0)->Print("all");

   // auto t1 = rPhiEventScene->FindChild("Tracks [P]")->FirstChild();
   // printf("t1=%p, %s %s\n", t1, t1->GetElementName(), t1->IsA()->GetName());
   // dynamic_cast<REX::REveTrack*>(t1)->Print("all");
}

//==============================================================================

class EventManager : public REX::REveElement
{
public:
   EventManager() = default;

   virtual ~EventManager() {}

   virtual void NextEvent()
   {
      eveMng->DisableRedraw();
      auto scene =  eveMng->GetEventScene();
      scene->DestroyElements();
      makeEventScene();
      for (auto &ie : scene->RefChildren())
      {
         if (mngRhoPhi)
         mngRhoPhi->ImportElements(ie, rPhiEventScene);
         if (mngRhoZ)
         mngRhoZ  ->ImportElements(ie, rhoZEventScene);
      }
      eveMng->EnableRedraw();
      eveMng->DoRedraw3D();
   }

   virtual void QuitRoot()
   {
      printf("Quit ROOT\n");
      if (gApplication) gApplication->Terminate();
   }

};

void event_demo()
{
   // disable browser cache - all scripts and html files will be loaded every time, useful for development
   // gEnv->SetValue("WebGui.HttpMaxAge", 0);

   gRandom->SetSeed(0); // make random seed

   eveMng = REX::REveManager::Create();

   auto eventMng = new EventManager();
   eventMng->SetName("EventManager");
   eveMng->GetWorld()->AddElement(eventMng);

   eveMng->GetWorld()->AddCommand("QuitRoot", "sap-icon://log", eventMng, "QuitRoot()");

   eveMng->GetWorld()->AddCommand("NextEvent", "sap-icon://step", eventMng, "NextEvent()");

   makeGeometryScene();
   makeEventScene();

   if (1) {
      createProjectionStuff();
      projectScenes(true, true);
   }

   eveMng->Show();
}
