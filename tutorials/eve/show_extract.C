/// \file
/// \ingroup tutorial_eve
/// Helper script for showing of extracted / simplified geometries.
/// By default shows a simplified ALICE geometry.
///
/// \image html eve_show_extract.png
/// \macro_code
///
/// \author Matevz Tadel

void show_extract(const char* file="http://root.cern/files/alice_ESDgeometry.root")
{
  TEveManager::Create();

  TFile::Open(file);

  TIter next(gDirectory->GetListOfKeys());
  TKey* key;
  TString xxx("TEveGeoShapeExtract");

  while ((key = (TKey*) next()))
  {
    if (xxx == key->GetClassName())
    {
      auto gse = (TEveGeoShapeExtract*) key->ReadObj();
      auto gs  = TEveGeoShape::ImportShapeExtract(gse, 0);
      gEve->AddGlobalElement(gs);
    }
  }

  gEve->Redraw3D(kTRUE);
}
