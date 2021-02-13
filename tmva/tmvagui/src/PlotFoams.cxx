#include "TMVA/PlotFoams.h"

#include "TControlBar.h"
#include "TVectorT.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TMVA/PDEFoamKernelBase.h"
#include "TMVA/PDEFoamKernelTrivial.h"

#include <sstream>
#include <string>
#include <cfloat>

#include "TMVA/PDEFoam.h"

void TMVA::PlotFoams( TString fileName,
                      bool useTMVAStyle )
{
   cout << "read file: " << fileName << endl;
   cout << "kValue = " << kValue << endl;
   TFile *file = TFile::Open(fileName);

   // set style and remove existing canvas'
   TMVAGlob::Initialize( useTMVAStyle );

   // create control bar
   TControlBar* cbar = new TControlBar( "vertical", "Choose cell value for plot:", 50, 50 );
   if ((gDirectory->Get("SignalFoam") && gDirectory->Get("BgFoam")) ||
       gDirectory->Get("MultiTargetRegressionFoam")) {
      TString macro = Form( "TMVA::Plot(\"%s\",%s, \"Event density\", %s)",
                            fileName.Data(), "TMVA::kValueDensity", (useTMVAStyle ? "kTRUE" : "kFALSE") );
      cbar->AddButton( "Event density", macro, "Plot event density", "button" );
   } else if (gDirectory->Get("DiscrFoam") || gDirectory->Get("MultiClassFoam0")){
      TString macro = Form( "TMVA::Plot(\"%s\", %s, \"Discriminator\", %s)",
                            fileName.Data(), "TMVA::kValue", (useTMVAStyle ? "kTRUE" : "kFALSE") );
      cbar->AddButton( "Discriminator", macro, "Plot discriminator", "button" );
   } else if (gDirectory->Get("MonoTargetRegressionFoam")){
      TString macro = Form( "TMVA::Plot(\"%s\", %s, \"Target\", %s)",
                            fileName.Data(), "TMVA::kValue",  (useTMVAStyle ? "kTRUE" : "kFALSE") );
      cbar->AddButton( "Target", macro, "Plot target", "button" );
   } else {
      cout << "Error: no foams found in file: " << fileName << endl;
      return;
   }

   TString macro_rms = Form( "TMVA::Plot(\"%s\", %s, \"Variance\", %s)",
                             fileName.Data(), "TMVA::kRms", (useTMVAStyle ? "kTRUE" : "kFALSE") );
   cbar->AddButton( "Variance", macro_rms, "Plot variance", "button" );
   TString macro_rms_ov_mean = Form( "TMVA::Plot(\"%s\", %s, \"Variance/Mean\", %s)",
                                     fileName.Data(), "TMVA::kRmsOvMean", (useTMVAStyle ? "kTRUE" : "kFALSE") );
   cbar->AddButton( "Variance/Mean", macro_rms_ov_mean, "Plot variance over mean", "button" );
   TString macro_cell_tree = Form( "TMVA::PlotCellTree(\"%s\", \"Cell tree\", %s)",
                                   fileName.Data(), (useTMVAStyle ? "kTRUE" : "kFALSE") );
   cbar->AddButton( "Cell tree", macro_cell_tree, "Plot cell tree", "button" );

   cbar->Show();
   file->Close();
}

// foam plotting macro
void TMVA::Plot(TString fileName, TMVA::ECellValue cv, TString cv_long, bool useTMVAStyle )
{
   cout << "read file: " << fileName << endl;
   TFile *file = TFile::Open(fileName);

   gStyle->SetNumberContours(999);
   if (useTMVAStyle)  TMVAGlob::SetTMVAStyle();

   // fileNamed foams and foam type
   TMVA::PDEFoam* SignalFoam = (TMVA::PDEFoam*) gDirectory->Get("SignalFoam");
   TMVA::PDEFoam* BgFoam     = (TMVA::PDEFoam*) gDirectory->Get("BgFoam");
   TMVA::PDEFoam* DiscrFoam  = (TMVA::PDEFoam*) gDirectory->Get("DiscrFoam");
   TMVA::PDEFoam* MultiClassFoam0 = (TMVA::PDEFoam*) gDirectory->Get("MultiClassFoam0");
   TMVA::PDEFoam* MonoTargetRegressionFoam = (TMVA::PDEFoam*) gDirectory->Get("MonoTargetRegressionFoam");
   TMVA::PDEFoam* MultiTargetRegressionFoam = (TMVA::PDEFoam*) gDirectory->Get("MultiTargetRegressionFoam");
   TList foam_list; // the foams and their captions
   if (SignalFoam && BgFoam) {
      foam_list.Add(new TPair(SignalFoam, new TObjString("Signal Foam")));
      foam_list.Add(new TPair(BgFoam, new TObjString("Background Foam")));
   } else if (DiscrFoam) {
      foam_list.Add(new TPair(DiscrFoam, new TObjString("Discriminator Foam")));
   } else if (MultiClassFoam0) {
      UInt_t cls = 0;
      TMVA::PDEFoam *fm = NULL;
      while ((fm = (TMVA::PDEFoam*) gDirectory->Get(Form("MultiClassFoam%u", cls)))) {
         foam_list.Add(new TPair(fm, new TObjString(Form("Discriminator Foam %u",cls))));
         cls++;
      }
   } else if (MonoTargetRegressionFoam) {
      foam_list.Add(new TPair(MonoTargetRegressionFoam,
                              new TObjString("MonoTargetRegression Foam")));
   } else if (MultiTargetRegressionFoam) {
      foam_list.Add(new TPair(MultiTargetRegressionFoam,
                              new TObjString("MultiTargetRegression Foam")));
   } else {
      cout << "ERROR: no Foams found in file: " << fileName << endl;
      return;
   }

   // loop over all foams and print out a debug message
   TListIter foamIter(&foam_list);
   TPair *fm_pair = NULL;
   Int_t kDim = 0; // foam dimensions
   while ((fm_pair = (TPair*) foamIter())) {
      kDim = ((TMVA::PDEFoam*) fm_pair->Key())->GetTotDim();
      cout << "Foam loaded: " << ((TObjString*) fm_pair->Value())->String()
           << " (dimension = " << kDim << ")" << endl;
   }

   // kernel to use for the projection
   TMVA::PDEFoamKernelBase *kernel = new TMVA::PDEFoamKernelTrivial();

   // plot foams
   if (kDim == 1) {
      Plot1DimFoams(foam_list, cv, cv_long, kernel);
   } else {
      PlotNDimFoams(foam_list, cv, cv_long, kernel);
   }

   file->Close();
}


void TMVA::Plot1DimFoams(TList& foam_list, TMVA::ECellValue cell_value,
                         const TString& cell_value_description,
                         TMVA::PDEFoamKernelBase* kernel)
{
   // visualize a 1 dimensional PDEFoam via a histogram
   TCanvas* canvas = NULL;
   TH1D* projection = NULL;

   // loop over all foams and draw the histogram
   TListIter it(&foam_list);
   TPair* fm_pair = NULL;    // the (foam, caption) pair
   while ((fm_pair = (TPair*) it())) {
      TMVA::PDEFoam* foam = (TMVA::PDEFoam*) fm_pair->Key();
      if (!foam) continue;
      TString foam_caption(((TObjString*) fm_pair->Value())->String());
      TString variable_name(foam->GetVariableName(0)->String());

      canvas = new TCanvas(Form("canvas_%p",foam),
                           "1-dimensional PDEFoam", 400, 400);

      projection = foam->Draw1Dim(cell_value, 100, kernel);
      projection->SetTitle(cell_value_description + " of " + foam_caption
                           + ";" + variable_name);
      projection->Draw();
      projection->SetDirectory(0);

      canvas->Update();
   }
}


void TMVA::PlotNDimFoams(TList& foam_list, TMVA::ECellValue cell_value,
                         const TString& cell_value_description,
                         TMVA::PDEFoamKernelBase* kernel)
{
   // draw 2 dimensional PDEFoam projections
   TCanvas* canvas = NULL;
   TH2D* projection = NULL;

   // loop over all foams and draw the projection
   TListIter it(&foam_list);
   TPair* fm_pair = NULL;    // the (foam, caption) pair
   while ((fm_pair = (TPair*) it())) {
      TMVA::PDEFoam* foam = (TMVA::PDEFoam*) fm_pair->Key();
      if (!foam) continue;
      TString foam_caption(((TObjString*) fm_pair->Value())->String());
      const Int_t kDim = ((TMVA::PDEFoam*) fm_pair->Key())->GetTotDim();

      // draw all possible projections (kDim*(kDim-1)/2)
      for (Int_t i = 0; i < kDim; ++i) {
         for (Int_t k = i + 1; k < kDim; ++k) {

            canvas = new TCanvas(Form("canvas_%p_%i:%i", foam, i, k),
                                 Form("Foam projections %i:%i", i, k),
                                 (Int_t)(400/(1.-0.2)), 400);
            canvas->SetRightMargin(0.2);

            TString title = Form("%s of %s: Projection %s:%s;%s;%s",
                                 cell_value_description.Data(),
                                 foam_caption.Data(),
                                 foam->GetVariableName(i)->String().Data(),
                                 foam->GetVariableName(k)->String().Data(),
                                 foam->GetVariableName(i)->String().Data(),
                                 foam->GetVariableName(k)->String().Data());

            projection = foam->Project2(i, k, cell_value, kernel);
            projection->SetTitle(title);
            projection->Draw("COLZ");
            projection->SetDirectory(0);

            canvas->Update();
         }
      }
   } // loop over foams
}


void TMVA::PlotCellTree(TString fileName, TString cv_long, bool useTMVAStyle )
{
   // Draw the PDEFoam cell tree

   cout << "read file: " << fileName << endl;
   TFile *file = TFile::Open(fileName);

   if (useTMVAStyle)  TMVAGlob::SetTMVAStyle();

   // find foams
   TListIter foamIter(gDirectory->GetListOfKeys());
   TKey *foam_key = NULL; // the foam key
   TCanvas *canv = NULL;  // the canvas
   while ((foam_key = (TKey*) foamIter())) {
      TString name(foam_key->GetName());
      TString class_name(foam_key->GetClassName());
      if (!class_name.Contains("PDEFoam"))
         continue;
      cout << "PDEFoam found: " << class_name << " " << name << endl;

      // read the foam
      TMVA::PDEFoam *foam = (TMVA::PDEFoam*) foam_key->ReadObj();
      canv = new TCanvas(Form("canvas_%s",name.Data()),
                         Form("%s of %s",cv_long.Data(),name.Data()), 640, 480);
      canv->cd();
      // get cell tree depth
      const UInt_t   depth = foam->GetRootCell()->GetTreeDepth();
      const Double_t ystep = 1.0 / depth;
      DrawCell(foam->GetRootCell(), foam, 0.5, 1.-0.5*ystep, 0.25, ystep);
   }

   file->Close();
}

void TMVA::DrawCell( TMVA::PDEFoamCell *cell, TMVA::PDEFoam *foam,
                     Double_t x, Double_t y,
                     Double_t xscale,  Double_t yscale )
{
   // recursively draw cell and it's daughters

   Float_t xsize = xscale*1.5;
   Float_t ysize = yscale/3;
   if (xsize > 0.15) xsize=0.1; //xscale/2;
   if (cell->GetDau0() != NULL) {
      TLine *a1 = new TLine(x-xscale/4, y-ysize, x-xscale, y-ysize*2);
      a1->SetLineWidth(2);
      a1->Draw();
      DrawCell(cell->GetDau0(), foam, x-xscale, y-yscale, xscale/2, yscale);
   }
   if (cell->GetDau1() != NULL){
      TLine *a1 = new TLine(x+xscale/4, y-ysize, x+xscale, y-ysize*2);
      a1->SetLineWidth(2);
      a1->Draw();
      DrawCell(cell->GetDau1(), foam, x+xscale, y-yscale, xscale/2, yscale);
   }

   TPaveText *t = new TPaveText(x-xsize, y-ysize, x+xsize, y+ysize, "NDC");

   t->SetBorderSize(1);
   t->SetFillStyle(1);

   // draw all cell elements
   t->AddText( Form("Intg=%.5f", cell->GetIntg()) );
   t->AddText( Form("Var=%.5f", cell->GetDriv()) );
   TVectorD *vec = (TVectorD*) cell->GetElement();
   if (vec) {
      for (Int_t i = 0; i < vec->GetNrows(); ++i)
         t->AddText( Form("E[%i]=%.5f", i, (*vec)[i]) );
   }

   if (cell->GetStat() != 1) {
      // cell is inactive --> draw split point
      t->SetFillColor( TColor::GetColor("#BBBBBB") );
      t->SetTextColor( TColor::GetColor("#000000") );

      // cell position and size
      TMVA::PDEFoamVect cellPosi(foam->GetTotDim()), cellSize(foam->GetTotDim());
      cell->GetHcub(cellPosi, cellSize);
      Int_t    kBest = cell->GetBest(); // best division variable
      Double_t xBest = cell->GetXdiv(); // best division point
      t->AddText( Form("dim=%i", kBest) );
      t->AddText( Form("cut=%.5g", foam->VarTransformInvers(kBest,cellPosi[kBest] + xBest*cellSize[kBest])) );
   } else {
      t->SetFillColor( TColor::GetColor("#DD0033") );
      t->SetTextColor( TColor::GetColor("#FFFFFF") );
   }

   t->Draw();
}
