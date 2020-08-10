/// \file
/// \ingroup tutorial_eve
/// Shows CMS geometry.
///
/// \image html eve_geom_cms.png
/// \macro_code
///
/// \author Matevz Tadel

void geom_cms()
{
   TEveManager::Create();

   TFile::SetCacheFileDir(".");
   gGeoManager = gEve->GetGeometry("http://root.cern/files/cms.root");
   gGeoManager->DefaultColors();

   auto top = gGeoManager->GetTopVolume()->FindNode("CMSE_1")->GetVolume();

   TEveGeoTopNode* trk = new TEveGeoTopNode(gGeoManager, top->FindNode("TRAK_1"));
   trk->SetVisLevel(6);
   gEve->AddGlobalElement(trk);

   auto calo = new TEveGeoTopNode(gGeoManager, top->FindNode("CALO_1"));
   calo->SetVisLevel(3);
   gEve->AddGlobalElement(calo);

   auto muon = new TEveGeoTopNode(gGeoManager, top->FindNode("MUON_1"));
   muon->SetVisLevel(4);
   gEve->AddGlobalElement(muon);

   gEve->FullRedraw3D(kTRUE);

   // EClipType not exported to CINT (see TGLUtil.h):
   // 0 - no clip, 1 - clip plane, 2 - clip box
   auto v = gEve->GetDefaultGLViewer();
   v->GetClipSet()->SetClipType(TGLClip::EType(1));
   v->ColorSet().Background().SetColor(kMagenta+4);
   v->SetGuideState(TGLUtil::kAxesEdge, kTRUE, kFALSE, 0);
   v->RefreshPadEditor(v);

   v->CurrentCamera().RotateRad(-1.2, 0.5);
   v->DoDraw();
}
