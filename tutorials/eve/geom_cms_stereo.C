/// \file
/// \ingroup tutorial_eve
/// Shows CMS geometry in stereo mode.
/// This requires quad-buffer support in the OpenGL hardware / driver,
/// otherwise a fatal error occurs.
///
/// \image html eve_geom_cms_stereo.png
/// \macro_code
///
/// \author Matevz Tadel

void geom_cms_stereo(Bool_t quad_buf=kTRUE)
{
   TEveManager::Create();

   TFile::SetCacheFileDir(".");
   gGeoManager = gEve->GetGeometry("http://root.cern/files/cms.root");
   gGeoManager->DefaultColors();

   auto top = gGeoManager->GetTopVolume()->FindNode("CMSE_1")->GetVolume();

   auto trk = new TEveGeoTopNode(gGeoManager, top->FindNode("TRAK_1"));
   trk->SetVisLevel(6);
   gEve->AddGlobalElement(trk);

   auto calo = new TEveGeoTopNode(gGeoManager, top->FindNode("CALO_1"));
   calo->SetVisLevel(3);
   gEve->AddGlobalElement(calo);

   auto muon = new TEveGeoTopNode(gGeoManager, top->FindNode("MUON_1"));
   muon->SetVisLevel(4);
   gEve->AddGlobalElement(muon);

   // --- Stereo ---

   TEveWindowSlot* slot = 0;
   slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());

   auto  sv = new TEveViewer("Stereo GL", "Stereoscopic view");
   sv->SpawnGLViewer(gEve->GetEditor(), kTRUE, quad_buf);
   sv->AddScene(gEve->GetGlobalScene());

   slot->ReplaceWindow(sv);

   gEve->GetViewers()->AddElement(sv);

   gEve->GetBrowser()->GetTabRight()->SetTab(1);

   // --- Redraw ---

   gEve->FullRedraw3D(kTRUE);
   gEve->EditElement(sv);

   // --- Fix editor ---

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
