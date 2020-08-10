/// \file
/// \ingroup tutorial_eve
/// Shows LHCB geometry.
///
/// \image html eve_geom_lhcb.png
/// \macro_code
///
/// \author Matevz Tadel

void geom_lhcb()
{
   TEveManager::Create();

   TFile::SetCacheFileDir(".");
   gGeoManager = gEve->GetGeometry("http://root.cern/files/lhcbfull.root");
   gGeoManager->DefaultColors();

   auto tn = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
   tn->SetVisLevel(4);
   gEve->AddGlobalElement(tn);

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
