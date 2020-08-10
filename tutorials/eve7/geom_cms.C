#include <sstream>
#include <iostream>

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TSystem.h"
#include "TFile.h"
#include <ROOT/REveGeoShape.hxx>
#include <ROOT/REveTrans.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveViewer.hxx>
#include <ROOT/REveElement.hxx>
#include <ROOT/REveManager.hxx>

R__LOAD_LIBRARY(libGeom);
R__LOAD_LIBRARY(libROOTEve);

namespace REX = ROOT::Experimental;

void makeEveGeoShape(TGeoNode* n, REX::REveTrans& trans, REX::REveElement* holder)
{
   auto gss = n->GetVolume()->GetShape();
   auto b1s = new REX::REveGeoShape(n->GetName());
   b1s->InitMainTrans();
   b1s->RefMainTrans().SetFrom(trans.Array());
   b1s->SetShape(gss);
   b1s->SetMainColor(kCyan);
   holder->AddElement(b1s);
}


void filterChildNodes(TGeoNode* pn, REX::REveTrans& trans,  REX::REveElement* holder, std::string mat, int maxlevel, int level)
{
   ++level;
   if (level > maxlevel)
      return;

   for (int i = 0; i < pn->GetNdaughters(); ++i)
   {
      TGeoNode* n = pn->GetDaughter(i);
      TGeoMaterial* material = n->GetVolume()->GetMaterial();
      REX::REveTrans ctrans;
      ctrans.SetFrom(trans.Array());

      {
         TGeoMatrix     *gm = n->GetMatrix();
         const Double_t *rm = gm->GetRotationMatrix();
         const Double_t *tv = gm->GetTranslation();
         REX::REveTrans t;
         t(1,1) = rm[0]; t(1,2) = rm[1]; t(1,3) = rm[2];
         t(2,1) = rm[3]; t(2,2) = rm[4]; t(2,3) = rm[5];
         t(3,1) = rm[6]; t(3,2) = rm[7]; t(3,3) = rm[8];
         t(1,4) = tv[0]; t(2,4) = tv[1]; t(3,4) = tv[2];
         ctrans *= t;
      }

      std::string mn = material->GetName();
      if (mn == mat) {
         n->ls();
         makeEveGeoShape(n, ctrans, holder);
      }
      filterChildNodes(n, ctrans,holder,  mat, maxlevel, level);
   }
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

void geom_cms()
{

   auto eveMng = REX::REveManager::Create();

   TFile::SetCacheFileDir(".");

   auto geoManager = eveMng->GetGeometry("http://root.cern/files/cms.root");
   TGeoNode* top = geoManager->GetTopVolume()->FindNode("CMSE_1");

   // tracker
   {
      auto holder = new REX::REveElement("Tracker");
      eveMng->GetGlobalScene()->AddElement(holder);
      TGeoNode* n = getNodeFromPath(top, "TRAK_1/SVTX_1/TGBX_1/GAW1_1");
      REX::REveTrans trans;
      std::string material = "TOB_Silicon";
      filterChildNodes(n, trans, holder, material,  6, 0);
   }

   // muon
   {
      auto holder = new REX::REveElement("MUON");
      eveMng->GetGlobalScene()->AddElement(holder);

      auto n = getNodeFromPath(top, "MUON_1/MB_1");

      std::string material = "M_B_Air";
      REX::REveTrans trans;
      filterChildNodes(n, trans, holder, material, 1, 0);

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
   }

   eveMng->Show();
}
