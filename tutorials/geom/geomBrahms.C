/// \file
/// \ingroup tutorial_geom
/// Script drawing a detector geometry (here BRAHMS).
///
/// by default the geometry is drawn using the GL viewer
/// Using the TBrowser, you can select other components
/// if the file containing the geometry is not found in the local
/// directory, it is automatically read from the ROOT web site.
///
/// \image html geom_geomBrahms.png width=800px
/// \macro_code
///
/// \author Rene Brun

void geomBrahms() {
   TGeoManager::Import("http://root.cern.ch/files/brahms.root");
   gGeoManager->GetVolume("CAVE")->Draw("ogl");
   new TBrowser;
}
