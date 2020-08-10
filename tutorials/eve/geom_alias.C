/// \file
/// \ingroup tutorial_eve
/// Demonstrates usage of geometry aliases - merge ALICE ITS with ATLAS MUON.
///
/// \image html eve_geom_alias.png
/// \macro_code
///
/// \author Matevz Tadel

void geom_alias()
{
   TEveManager::Create();

   gEve->RegisterGeometryAlias("ALICE", "http://root.cern/files/alice.root");
   gEve->RegisterGeometryAlias("ATLAS", "http://root.cern/files/atlas.root");


   gGeoManager = gEve->GetGeometryByAlias("ALICE");

   auto node1 = gGeoManager->GetTopVolume()->FindNode("ITSV_1");
   auto its   = new TEveGeoTopNode(gGeoManager, node1);
   gEve->AddGlobalElement(its);

   gGeoManager = gEve->GetGeometryByAlias("ATLAS");

   auto node2 = gGeoManager->GetTopVolume()->FindNode("OUTE_1");
   auto atlas = new TEveGeoTopNode(gGeoManager, node2);
   gEve->AddGlobalElement(atlas);

   gEve->FullRedraw3D(kTRUE);

   // EClipType not exported to CINT (see TGLUtil.h):
   // 0 - no clip, 1 - clip plane, 2 - clip box
   auto v = gEve->GetDefaultGLViewer();
   v->GetClipSet()->SetClipType(TGLClip::EType(2));
   v->RefreshPadEditor(v);

   v->CurrentCamera().RotateRad(-0.5, -2.4);
   v->DoDraw();
}
