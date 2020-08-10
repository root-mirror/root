/// \file
/// \ingroup tutorial_eve
/// Demonstrates usage of automatic 2D projections - class TEveProjectionManager.
///
/// \image html eve_projection.png
/// \macro_code
///
/// \author Matevz Tadel

const char* esd_geom_file_name =
   "http://root.cern/files/alice_ESDgeometry.root";

void projection()
{
   TFile::SetCacheFileDir(".");
   TEveManager::Create();

   // camera
   auto s = gEve->SpawnNewScene("Projected Event");
   gEve->GetDefaultViewer()->AddScene(s);
   auto v = gEve->GetDefaultGLViewer();
   v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   TGLOrthoCamera& cam = (TGLOrthoCamera&) v->CurrentCamera();
   cam.SetZoomMinMax(0.2, 20);

   // projections
   auto mng = new TEveProjectionManager(TEveProjection::kPT_RPhi);
   s->AddElement(mng);
   auto axes = new TEveProjectionAxes(mng);
   axes->SetTitle("TEveProjections demo");
   s->AddElement(axes);
   gEve->AddToListTree(axes, kTRUE);
   gEve->AddToListTree(mng, kTRUE);

   // Simple geometry
   auto geom = TFile::Open(esd_geom_file_name, "CACHEREAD");
   if (!geom)
      return;

   auto gse = (TEveGeoShapeExtract*) geom->Get("Gentle");
   auto gsre = TEveGeoShape::ImportShapeExtract(gse, 0);
   geom->Close();
   delete geom;
   gsre->SetPickableRecursively(kTRUE);
   gEve->AddGlobalElement(gsre);
   gEve->GetGlobalScene()->SetRnrState(kFALSE);
   mng->ImportElements(gsre);

   auto line = new TEveLine;
   line->SetMainColor(kGreen);
   for (Int_t i=0; i<160; ++i)
      line->SetNextPoint(120*sin(0.2*i), 120*cos(0.2*i), 80-i);
   gEve->AddElement(line);
   mng->ImportElements(line);
   line->SetRnrSelf(kFALSE);

   gEve->Redraw3D(kTRUE);
}
