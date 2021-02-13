#include <ROOT/REveElement.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveViewer.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REvePointSet.hxx>


using namespace ROOT::Experimental;


void makeTracks(int N_Tracks, REveElement* trackHolder, float* pos)
{
   TRandom &r = *gRandom;
   auto prop = new REveTrackPropagator();
   prop->SetMagFieldObj(new REveMagFieldDuo(5, 3.5, -2.0));
   prop->SetMaxR(9);
   prop->SetMaxZ(15);

   double v = 0.1;
   double m = 5;
   for (int i = 0; i < N_Tracks; i++)
   {
      auto p = new TParticle();
      int pdg = 11 * (r.Integer(2) > 0 ? 1 : -1);
      p->SetPdgCode(pdg);
      p->SetProductionVertex(r.Uniform(-v,v) + pos[0], r.Uniform(-v,v) + pos[1], r.Uniform(-v,v) + pos[2], 1);
      p->SetMomentum(r.Uniform(-m,m), r.Uniform(-m,m), r.Uniform(-m,m)*r.Uniform(1, 3), 1);
      auto track = new REveTrack(p, 1, prop);
      track->MakeTrack();
      track->SetMainColor(kBlue);
      track->SetName(Form("Track_%d", i));
      trackHolder->AddElement(track);
   }
}

void makeProjected(REveElement* el, const char* pname, REveProjection::EPType_e t)
{
   auto eveMng = ROOT::Experimental::gEve;
   auto eventScene = eveMng->SpawnNewScene(Form("%s Event Data", pname), pname);
   auto mng = new REveProjectionManager();
   mng->SetProjection(t);
   mng->ImportElements(el, eventScene);
   auto view = eveMng->SpawnNewViewer(pname);
   view->AddScene(eventScene);
}

void error_ellipse()
{
   auto  eveMng = REveManager::Create();
   float pos[3] = {1.46589e-06,-1.30522e-05,-1.98267e-05};

   // symnetric matrix
   double a[16] = {1.46589e-01,-1.30522e-02,-1.98267e-02, 0,
      -1.30522e-02, 4.22955e-02,-5.86628e-03, 0,
      -1.98267e-02,-5.86628e-03, 2.12836e-01, 0,
      0, 0, 0, 1};

   REveTrans t;
   t.SetFrom(a);
   TMatrixDSym xxx(3);
   for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
      {
         xxx(i,j) = t(i+1,j+1);
      }

   TMatrixDEigen eig(xxx);
   TVectorD xxxEig ( eig.GetEigenValues() );
   xxxEig = xxxEig.Sqrt();

   TMatrixD vecEig = eig.GetEigenVectors();
   REveVector v[3]; int ei = 0;
   for (int i = 0; i < 3; ++i)
   {
      v[i].Set(vecEig(0,i), vecEig(1,i), vecEig(2,i));
      v[i] *=  xxxEig(i);
   }
   REveElement *event = eveMng->GetEventScene();

   REveEllipsoid* ellipse = new REveEllipsoid("VertexError");
   ellipse->InitMainTrans();
   ellipse->SetMainColor(kGreen + 10);
   ellipse->SetLineWidth(2);
   ellipse->SetBaseVectors(v[0], v[1], v[2]);
   ellipse->Outline();
   event->AddElement(ellipse);
   
   auto ps = new REvePointSet("Vertices");
   ps->SetMainColor(kYellow);
   ps->SetNextPoint(pos[0], pos[1], pos[2]);
   ps->SetMarkerStyle(4);
   ps->SetMarkerSize(4);
   float rng = 1;
   for(int i=0; i < 6; ++i)
      ps->SetNextPoint(i*(rng/3) - rng, 0, 0);
   event->AddElement(ps);
   
   auto trackHolder = new REveElement("Tracks");
   eveMng->GetEventScene()->AddElement(trackHolder);
   makeTracks(10, trackHolder, pos);
   
   makeProjected(eveMng->GetEventScene(), "RPhi", REveProjection::kPT_RPhi);
   makeProjected(eveMng->GetEventScene(), "RhoZ", REveProjection::kPT_RhoZ);
   eveMng->Show();
}
