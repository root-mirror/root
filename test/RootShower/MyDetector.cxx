// Author: Bertrand Bellenot   22/08/02

/*************************************************************************
 * Copyright (C) 1995-2002, Bertrand Bellenot.                           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see the LICENSE file.                         *
 *************************************************************************/

#include <cstdlib>
#include "TMath.h"

#include "MyDetector.h"

//______________________________________________________________________________
//
// MyDetector class implementation
//______________________________________________________________________________


ClassImp(MyDetector);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

MyDetector::MyDetector()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

MyDetector::~MyDetector()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize detector with material and dimensions

void MyDetector::Init()
{
   Int_t           i;
   Double_t        x;
   TGeoMaterial    *material;
   TGeoMaterial    *lead;
   TGeoMaterial    *polystyrene;
   TGeoMaterial    *nai;
   TGeoMedium      *scintillator[4];
   TGeoMedium      *discriminator;
   TGeoMedium      *calorimeter;
   TGeoVolume      *volume[7];
   TGeoTranslation *trans[6];

   fTotalELoss = 0.0;
   fDimX = 40.0;
   fDimY = 62.0;
   fDimZ = 40.0;
   fMaxX = fDimX / 2.;
   fMinX = -fMaxX;

   fMinY = 0.0;
   fMaxY = fDimY+fMinY;

   fMaxZ = fDimZ / 2.;
   fMinZ = -fMaxZ;

   new TGeoManager("MyDetector", "MyDetector");

   trans[0] = new TGeoTranslation(0., 2., 0.);
   trans[1] = new TGeoTranslation(0., 5.0, 0.);
   trans[2] = new TGeoTranslation(0., 8.0, 0.);
   trans[3] = new TGeoTranslation(0., 12.0, 0.);
   trans[4] = new TGeoTranslation(0., 36.0, 0.);
   trans[5] = new TGeoTranslation(0., 60.0, 0.);

   new TGeoMaterial("Iron",55.85f,26,7.87f,1.76f,16.7598f);
   lead = new TGeoMaterial("Lead",207.2f,82,11.35f,0.56f,17.0925f);
   polystyrene = new TGeoMaterial("Polystyrene",13.01f,7,1.032f,42.4f,79.36f);
   new TGeoMaterial("BGO",175.92f,74,7.1f,1.12f,22.11f);
   new TGeoMaterial("CsI",129.90f,54,4.53f,1.85f,36.8653f);
   nai = new TGeoMaterial("NaI",117.10f,50,3.67f,2.59f,41.1444f);
   new TGeoMaterial("Al",26.981539f,13,2.7f,8.9f,39.407407f);

   scintillator[0] = new TGeoMedium("SCINT0",1, polystyrene);
   discriminator   = new TGeoMedium("DISCR", 2, lead);
   scintillator[1] = new TGeoMedium("SCINT1",3, polystyrene);
   scintillator[2] = new TGeoMedium("SCINT2",4, polystyrene);
   calorimeter     = new TGeoMedium("CALOR", 5, nai);
   scintillator[3] = new TGeoMedium("SCINT3",6, polystyrene);

   TGeoMaterial *mat = new TGeoMaterial("VOID");
   TGeoMedium *med = new TGeoMedium("MED", 0, mat);
   volume[0] = gGeoManager->MakeBox("TOP", med, 40, 62, 40);
   volume[0]->SetVisibility(kFALSE);
   gGeoManager->SetTopVolume(volume[0]);

   volume[1] = gGeoManager->MakeBox("BOX0",scintillator[0], 20.0,2.00,20.0);
   volume[1]->SetLineColor(7);
   volume[1]->SetTransparency(50);
   volume[1]->SetLineWidth(1);
   volume[0]->AddNode(volume[1],0,trans[0]);

   volume[2] = gGeoManager->MakeBox("BOX1",discriminator, 20.0,1.00,20.0);
   volume[2]->SetLineColor(15);
   volume[2]->SetTransparency(50);
   volume[2]->SetLineWidth(1);
   volume[0]->AddNode(volume[2],1,trans[1]);

   volume[3] = gGeoManager->MakeBox("BOX2",scintillator[1], 20.0,2.00,20.0);
   volume[3]->SetLineColor(7);
   volume[3]->SetTransparency(50);
   volume[3]->SetLineWidth(1);
   volume[0]->AddNode(volume[3],2,trans[2]);

   volume[4] = gGeoManager->MakeBox("BOX3",scintillator[2], 20.0,2.00,20.0);
   volume[4]->SetLineColor(7);
   volume[4]->SetTransparency(50);
   volume[4]->SetLineWidth(1);
   volume[0]->AddNode(volume[4],3,trans[3]);

   volume[5] = gGeoManager->MakeBox("BOX4",calorimeter, 20.0,22.00,20.0);
   volume[5]->SetLineColor(38);
   volume[5]->SetTransparency(50);
   volume[5]->SetLineWidth(1);
   volume[0]->AddNode(volume[5],4,trans[4]);

   volume[6] = gGeoManager->MakeBox("BOX5",scintillator[3], 20.0,2.00,20.0);
   volume[6]->SetLineColor(7);
   volume[6]->SetTransparency(50);
   volume[6]->SetLineWidth(1);
   volume[0]->AddNode(volume[6],5,trans[5]);

   gGeoManager->CloseGeometry();

   for (i=0;i<6;i++) {
      material = volume[i+1]->GetMaterial();
      // Ionisation constant in MeV/cm
      fI[i] = 16.0e-06 * (TMath::Power(material->GetZ(),0.9));
      fI[i] *= 1.0e-03; // in GeV/cm...

      // first factor in the Bethe-Bloch equation in MeV/cm
      fPreconst[i] = 0.3071 * material->GetDensity() *
                     (material->GetZ() / material->GetA());
      fPreconst[i] *= 1.0e-03; // in GeV/cm...

      // Critical Energy in MeV
      fEc[i] = 800.0 / (material->GetZ() + 1.2);
      fEc[i] *= 1.0e-03; // in GeV...

      // Radiation Length in cm
      fX0[i] = material->GetRadLen();

      // Time step dT in ms about 0.015 times X_0/c
      fdT[i] = 0.005 * (fX0[i] / CSpeed);

      // Scatter angle
      x = fdT[i] * CSpeed;
      fTheta0[i] = TMath::Sqrt(2.0) * 13.6 * TMath::Sqrt(x / fX0[i]) *
                   (1.0 + 0.038 * TMath::Log(x / fX0[i]));
      fTheta0[i] *= 1.0e-03;
   }
}

