/// \file
/// \ingroup tutorial_geom
/// Script drawing a detector geometry (here D0).
///
/// by default the geometry is drawn using the GL viewer
/// Using the TBrowser, you can select other components
/// if the file containing the geometry is not found in the local
/// directory, it is automatically read from the ROOT web site.
///   - run with `.x geomD0.C`    top level detectors are transparent
///   - or       `.x geomD0.C(1)` top level detectors are visible
///
/// \image html geom_geomD0.png width=800px
/// \macro_code
///
/// \authors Bertrand Bellenot, Rene Brun

void RecursiveInvisible(TGeoVolume *vol);
void RecursiveTransparency(TGeoVolume *vol, Int_t transp);

void geomD0(Int_t allVisible=0) {
   TGeoManager::Import("http://root.cern.ch/files/d0.root");
   gGeoManager->DefaultColors();
   gGeoManager->SetMaxVisNodes(40000);
   //gGeoManager->SetVisLevel(4);
   if (!allVisible) {
      RecursiveInvisible(gGeoManager->GetVolume("D0-"));
      RecursiveInvisible(gGeoManager->GetVolume("D0+"));
      RecursiveInvisible(gGeoManager->GetVolume("D0WZ"));
      RecursiveInvisible(gGeoManager->GetVolume("D0WL"));
      RecursiveTransparency(gGeoManager->GetVolume("MUON"), 90);
   }

   gGeoManager->GetVolume("D0")->Draw("ogl");
}

void RecursiveInvisible(TGeoVolume *vol)
{
   vol->InvisibleAll();
   Int_t nd = vol->GetNdaughters();
   for (Int_t i=0; i<nd; i++) {
      RecursiveInvisible(vol->GetNode(i)->GetVolume());
   }
}

void RecursiveTransparency(TGeoVolume *vol, Int_t transp)
{
   vol->SetTransparency(transp);
   Int_t nd = vol->GetNdaughters();
   for (Int_t i=0; i<nd; i++) {
      RecursiveTransparency(vol->GetNode(i)->GetVolume(), transp);
   }
}
