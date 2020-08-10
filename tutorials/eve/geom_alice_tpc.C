/// \file
/// \ingroup tutorial_eve
/// Shows geometry of ALICE TPC.
///
/// \image html eve_geom_alice_tpc.png
/// \macro_code
///
/// \author Matevz Tadel

void geom_alice_tpc()
{
   TEveManager::Create();

   gGeoManager = gEve->GetGeometry("http://root.cern/files/alice.root");

   auto node = gGeoManager->GetTopVolume()->FindNode("TPC_M_1");
   auto tpc  = new TEveGeoTopNode(gGeoManager, node);
   gEve->AddGlobalElement(tpc);

   gEve->Redraw3D(kTRUE);
}
