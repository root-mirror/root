/// \file
/// \ingroup tutorial_geom
/// Macro allowing to vizualize tessellations from Wavefront's .obj format.
///
/// \image html geom_visualizeWavefrontObj.png width=500px
/// \macro_code
///
/// \author Andrei Gheata

#include <TROOT.h>
#include <TColor.h>
#include <TDatime.h>
#include <TRandom3.h>
#include <TGeoManager.h>
#include <TGeoTessellated.h>

//______________________________________________________________________________
int randomColor()
{
   gRandom = new TRandom3();
   TDatime dt;
   gRandom->SetSeed(dt.GetTime());
   int ci = TColor::GetFreeColorIndex();
   TColor *color = new TColor(ci, gRandom->Rndm(), gRandom->Rndm(), gRandom->Rndm());
   return ci;
}

//______________________________________________________________________________
void visualizeWavefrontObj(const char *dot_obj_file="", bool check = false)
{
   // Input a file in .obj format (https://en.wikipedia.org/wiki/Wavefront_.obj_file)
   // The file should have a single object inside, only vertex and faces information is used

   TString name = dot_obj_file;
   TString sfile = dot_obj_file;
   if (sfile.IsNull()) {
      sfile = gROOT->GetTutorialsDir();
      sfile += "/geom/teddy.obj";
   }
   name.ReplaceAll(".obj", "");
   gROOT->GetListOfCanvases()->Delete();
   if (gGeoManager)
      delete gGeoManager;
   new TGeoManager(name, "Imported from .obj file");
   TGeoMaterial *mat = new TGeoMaterial("Al", 26.98, 13, 2.7);
   TGeoMedium *med = new TGeoMedium("MED", 1, mat);
   TGeoVolume *top = gGeoManager->MakeBox("TOP", med, 10, 10, 10);
   gGeoManager->SetTopVolume(top);

   auto tsl = TGeoTessellated::ImportFromObjFormat(sfile.Data(), check);
   if (!tsl) return;
   tsl->ResizeCenter(5.);

   TGeoVolume *vol = new TGeoVolume(name, tsl, med);
   vol->SetLineColor(randomColor());
   vol->SetLineWidth(2);
   top->AddNode(vol, 1);
   gGeoManager->CloseGeometry();
   if (!gROOT->IsBatch()) top->Draw("ogl");
}
