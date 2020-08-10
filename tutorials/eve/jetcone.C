/// \file
/// \ingroup tutorial_eve
/// Demonstrates usage of TEveJetCone class.
///
/// \image html eve_jetcone.png
/// \macro_code
///
/// \author Jochen Thaeder

const char* esd_geom_file_name =
   "http://root.cern/files/alice_ESDgeometry.root";

TEveVector GetTEveVector(Float_t eta, Float_t phi);
void       geomGentleTPC();

void jetcone()
{
   TEveManager::Create();

   using namespace TMath;

   TRandom r(0);

   // -- Set Constants
   Int_t nCones  = 10;
   Int_t nTracks = 200;

   Float_t coneRadius = 0.4;
   Float_t length = 300.;

   // -- Define palette
   auto pal = new TEveRGBAPalette(0, 500);

   // -----------------------------------------------------------------------
   // -- Line sets
   // -----------------------------------------------------------------------

   // -- Define cone center
   auto axis = new TEveStraightLineSet("Cone Axis");
   axis->SetLineColor(kGreen);
   axis->SetLineWidth(2);

   auto tracksXYZ = new TEveStraightLineSet("StraightLinesXYZ");
   tracksXYZ->SetLineColor(kRed);
   tracksXYZ->SetLineWidth(2);

   auto tracksEtaPhi = new TEveStraightLineSet("StraightLinesEtaPhi");
   tracksEtaPhi->SetLineColor(kYellow);
   tracksEtaPhi->SetLineWidth(2);

   auto tracksSeedEtaPhi = new TEveStraightLineSet("StraightLinesEtaPhiSeed");
   tracksSeedEtaPhi->SetLineColor(kBlue);
   tracksSeedEtaPhi->SetLineWidth(2);

   // -----------------------------------------------------------------------
   // -- Draw track distribution in XYZ in TPC Volume +/-250, +/-250, +/-250
   // -----------------------------------------------------------------------

   for ( Int_t track=0; track < nTracks ; track++ ) {

      Float_t trackX = r.Uniform(-250.0, 250.0);
      Float_t trackY = r.Uniform(-250.0, 250.0);
      Float_t trackZ = r.Uniform(-250.0, 250.0);
      Float_t trackR = Sqrt(trackX*trackX + trackY*trackY + trackZ*trackZ);

      TEveVector trackDir(trackX/trackR, trackY/trackR ,trackZ/trackR);
      TEveVector trackEnd = trackDir * length;
      tracksXYZ->AddLine(0., 0., 0., trackEnd.fX, trackEnd.fY, trackEnd.fZ );
   }

   // -----------------------------------------------------------------------
   // -- Draw track distribution in eta phi in TPC Volume +/-0.9, {0, 2Pi}
   // -----------------------------------------------------------------------

   for ( Int_t track=0; track < nTracks ; track++ ) {

      Float_t trackEta = r.Uniform(-0.9, 0.9);
      Float_t trackPhi = r.Uniform(0.0, TwoPi());

      TEveVector trackDir( GetTEveVector(trackEta, trackPhi) );
      TEveVector trackEnd = trackDir * length;

      if ( trackEta > coneRadius || trackEta < -coneRadius )
         tracksEtaPhi->AddLine(0., 0., 0.,
                               trackEnd.fX, trackEnd.fY, trackEnd.fZ);
      else
         tracksSeedEtaPhi->AddLine(0., 0., 0.,
                                   trackEnd.fX, trackEnd.fY, trackEnd.fZ);
   }

   // -----------------------------------------------------------------------
   // -- Draw cones
   // -----------------------------------------------------------------------

   for ( Int_t iter = 0; iter < nCones; ++iter ) {

      // -- Get Random ( eta ,phi )
      Float_t coneEta = r.Uniform(-0.9, 0.9);
      Float_t conePhi = r.Uniform(0.0, TwoPi() );

      // -- Primary vertex as origin
      TEveVector coneOrigin(0.0,0.0,0.0);

      // -- Get Cone Axis - axis line 10% longer than cone height
      TEveVector coneAxis ( GetTEveVector( coneEta, conePhi) );
      coneAxis *= length * 1.1;

      axis->AddLine( 0., 0., 0., coneAxis.fX, coneAxis.fY, coneAxis.fZ );

      // -- Draw jet cone
      auto jetCone = new TEveJetCone("JetCone");
      jetCone->SetPickable(kTRUE);
      jetCone->SetCylinder( 250., 250. );
      if ( (jetCone->AddCone( coneEta, conePhi, coneRadius   ) ) != -1)
         gEve->AddElement( jetCone );
   }

   // -----------------------------------------------------------------------

   // -- Add cone axis
   gEve->AddElement(axis);

   // -- Add lines
   //  gEve->AddElement(tracksXYZ);
   gEve->AddElement(tracksSeedEtaPhi);
   gEve->AddElement(tracksEtaPhi);

   // -- Load TPC geometry
   geomGentleTPC();

   gEve->Redraw3D(kTRUE);

   return;
}

//___________________________________________________________________________
TEveVector GetTEveVector(Float_t eta, Float_t phi)
{
  using namespace TMath;

  TEveVector vec( (Float_t) Cos ( (Double_t) phi)/ CosH( (Double_t) eta ),
                 (Float_t) Sin ( (Double_t) phi)/ CosH( (Double_t) eta ),
                 (Float_t) TanH( (Double_t) eta ) );
  return vec;
}

//__________________________________________________________________________
void geomGentleTPC()
{
   // Simple geometry
   TFile::SetCacheFileDir(".");
   TFile* geom = TFile::Open(esd_geom_file_name, "CACHEREAD");
   if (!geom)
      return;

   TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) geom->Get("Gentle");
   TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse, 0);
   geom->Close();
   delete geom;

   gEve->AddGlobalElement(gsre);

   TEveElement* elTRD = gsre->FindChild("TRD+TOF");
   elTRD->SetRnrState(kFALSE);

   TEveElement* elPHOS = gsre->FindChild("PHOS");
   elPHOS->SetRnrState(kFALSE);

   TEveElement* elHMPID = gsre->FindChild("HMPID");
   elHMPID->SetRnrState(kFALSE);
}
