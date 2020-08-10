/// \file
/// \ingroup tutorial_eve
/// Shows geometry of ALICE ITS.
///
/// \image html eve_geom_alice_its.png
/// \macro_code
///
/// \author Matevz Tadel

#include "TEveManager.h"
#include "TEveGeoNode.h"

#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"

#include "TString.h"

void geom_alice_its()
{
   TEveManager::Create();

   gGeoManager = gEve->GetGeometry("http://root.cern/files/alice.root");

   TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("ITSV_1");
   TEveGeoTopNode* its = new TEveGeoTopNode(gGeoManager, node);
   gEve->AddGlobalElement(its);

   gEve->Redraw3D(kTRUE);
}


//==============================================================================
// Demonstrate extraction of volumes matching certain criteria.
//==============================================================================

// Should be run in compiled mode -- CINT has issues with recursion.
//
// 1. Creation:
//    root
//      .L geom_alice_its.C+
//      extract_ssd_modules()
//      .q
//    This creates file "test-extract.root" in current dir.
//
// 2. Viewing:
//    root
//      .x show_extract.C("test-extract.root")

TEveGeoNode* descend_extract(TGeoNode* node)
{
   // We only return something if:
   // - this is a node of interest;
   // - one of the daughters returns something of interest.

   const TString material("ITS_SI$");

   TEveGeoNode *res = 0;

   auto medium = node->GetVolume()->GetMedium();
   if (medium && material == medium->GetName()) {
      // Node of interest - instantiate eve representation and return.
      res = new TEveGeoNode(node);
      return res;
   }

   Int_t nd = node->GetNdaughters();
   for (Int_t i = 0; i < nd; ++i) {
      auto ed = descend_extract(node->GetDaughter(i));

      if (ed) {
         if (res == 0) res = new TEveGeoNode(node);
         res->AddElement(ed);
      }
   }

   return res;
}

void extract_ssd_modules()
{
   const TString kEH("extract_ssd_modules");

   TEveManager::Create();

   gGeoManager = gEve->GetGeometry("http://root.cern/files/alice.root");

   Bool_t s = gGeoManager->cd("/ITSV_1/ITSD_1/IT56_1");
   if (!s) {
      Error(kEH, "Start node not found.");
      return;
   }

   auto node = gGeoManager->GetCurrentNode();

   TEveGeoNode *egn = descend_extract(node);

   if (egn == 0) {
      Warning(kEH, "No matching nodes found.");
      return;
   }

   egn->SaveExtract("test-extract.root", "AliSDD", kTRUE);
}
