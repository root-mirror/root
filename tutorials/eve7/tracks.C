/// \file
/// \ingroup tutorial_eve7
///  This example display only points in web browser
///
/// \macro_code
///

#include "TRandom.h"
#include "TParticle.h"
#include <ROOT/REveElement.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveTrack.hxx>
#include <ROOT/REveTrackPropagator.hxx>

namespace REX = ROOT::Experimental;

void makeTracks(int N_Tracks, REX::REveElement* trackHolder)
{
   TRandom &r = *gRandom;
   auto prop = new REX::REveTrackPropagator();
   prop->SetMagFieldObj(new REX::REveMagFieldDuo(350, 3.5, -2.0));
   prop->SetMaxR(300);
   prop->SetMaxZ(600);
   prop->SetMaxOrbs(6);
   // Default is kHelix propagator.
   // prop->SetStepper(REX::REveTrackPropagator::kRungeKutta);

   double v = 0.5;
   double m = 5;

   for (int i = 0; i < N_Tracks; i++)
   {
      auto p = new TParticle();

      int pdg = 11 * (r.Integer(2) > 0 ? 1 : -1);
      p->SetPdgCode(pdg);

      p->SetProductionVertex(r.Uniform(-v,v), r.Uniform(-v,v), r.Uniform(-v,v), 1);
      p->SetMomentum(r.Uniform(-m,m), r.Uniform(-m,m), r.Uniform(-m,m)*r.Uniform(1, 3), 1);
      auto track = new REX::REveTrack(p, 1, prop);
      track->MakeTrack();
      track->SetMainColor(kBlue);
      track->SetName(Form("RandomTrack_%d", i));
      trackHolder->AddElement(track);
   }
}

void tracks()
{
   auto eveMng = REX::REveManager::Create();

   auto trackHolder = new REX::REveElement("Tracks");
   eveMng->GetEventScene()->AddElement(trackHolder);
   makeTracks(100, trackHolder);

   eveMng->Show();
}
