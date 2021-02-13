// Author: Sergey Linev, GSI  13/01/2021

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ExampleWidget.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"

#include <ROOT/RCanvas.hxx>
#include <ROOT/TObjectDrawable.hxx>
#include <ROOT/REveGeomViewer.hxx>

#include <QMessageBox>

ExampleWidget::ExampleWidget(QWidget *parent, const char* name) :
   QWidget(parent)
{
   setupUi(this);

   setAttribute(Qt::WA_DeleteOnClose);

   setObjectName(name);

   // create sample histogram

   fHisto = new TH1F("gaus1","Example of TH1 drawing in TCanvas", 100, -5, 5);
   fHisto->FillRandom("gaus", 10000);
   fHisto->SetDirectory(nullptr);

   gPad = fxTCanvasWidget->getCanvas();
   fHisto->Draw();

   static constexpr int nth2points = 40;
   fHisto2 = std::make_shared<TH2I>("gaus2", "Example of TH2 drawing in RCanvas", nth2points, -5, 5, nth2points, -5, 5);
   fHisto2->SetDirectory(nullptr);
   for (int n=0;n<nth2points;++n) {
      for (int k=0;k<nth2points;++k) {
         double x = 10.*n/nth2points-5.;
         double y = 10.*k/nth2points-5.;
         fHisto2->SetBinContent(fHisto2->GetBin(n+1, k+1), (int) (1000*TMath::Gaus(x)*TMath::Gaus(y)));
      }
   }

   fxRCanvasWidget->getCanvas()->Draw<ROOT::Experimental::TObjectDrawable>(fHisto2, "col");

   CreateDummyGeometry();
}

ExampleWidget::~ExampleWidget()
{
}


void ExampleWidget::CreateDummyGeometry()
{
   auto viewer = fxGeomViewerWidget->getGeomViewer();

   new TGeoManager("tubeseg", "poza3");
   TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
   TGeoMedium *med = new TGeoMedium("MED",1,mat);
   TGeoVolume *top = gGeoManager->MakeBox("TOP",med,100,100,100);
   gGeoManager->SetTopVolume(top);
   TGeoVolume *vol = gGeoManager->MakeTubs("TUBESEG",med, 20,30,40,-30,270);
   vol->SetLineColor(kRed);
   vol->SetLineWidth(2);
   top->AddNode(vol,1);
   gGeoManager->CloseGeometry();
//   gGeoManager->SetNsegments(40);
   gGeoManager->SetNsegments(80);

   viewer->SetGeometry(gGeoManager);
}

void ExampleWidget::ImportCmsGeometry()
{
   TFile::SetCacheFileDir(".");

   if (!TGeoManager::Import("https://root.cern/files/cms.root")) {
      CreateDummyGeometry();
      return;
   }

   auto viewer = fxGeomViewerWidget->getGeomViewer();

   gGeoManager->DefaultColors();
   gGeoManager->SetVisLevel(4);
   gGeoManager->GetVolume("TRAK")->InvisibleAll();
   gGeoManager->GetVolume("HVP2")->SetTransparency(20);
   gGeoManager->GetVolume("HVEQ")->SetTransparency(20);
   gGeoManager->GetVolume("YE4")->SetTransparency(10);
   gGeoManager->GetVolume("YE3")->SetTransparency(20);
   gGeoManager->GetVolume("RB2")->SetTransparency(99);
   gGeoManager->GetVolume("RB3")->SetTransparency(99);
   gGeoManager->GetVolume("COCF")->SetTransparency(99);
   gGeoManager->GetVolume("HEC1")->SetLineColor(7);
   gGeoManager->GetVolume("EAP1")->SetLineColor(7);
   gGeoManager->GetVolume("EAP2")->SetLineColor(7);
   gGeoManager->GetVolume("EAP3")->SetLineColor(7);
   gGeoManager->GetVolume("EAP4")->SetLineColor(7);
   gGeoManager->GetVolume("HTC1")->SetLineColor(2);

   viewer->SetGeometry(gGeoManager);

   // select volume to draw
   viewer->SelectVolume("CMSE");

   // specify JSROOT draw options - here clipping on X,Y,Z axes
   viewer->SetDrawOptions("clipxyz");

   // set default limits for number of visible nodes and faces
   // when viewer created, initial values exported from TGeoManager
   viewer->SetLimits();
}

void ExampleWidget::InfoButton_clicked()
{
   QMessageBox::information(this,"QtRoot standalone example","Demo how QRootCanvas can be inserted into the QWidget");
}


void ExampleWidget::CmsButton_clicked()
{
   ImportCmsGeometry();
   fxTabWidget->setCurrentIndex(3);
}

void ExampleWidget::ExitButton_clicked()
{
   // when widget closed, application automatically exit
   close();
}
