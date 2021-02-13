#include "TMVA/correlations.h"
#include "TMVA/Config.h"


#include "TH2.h"
#include "TPaletteAxis.h"

// this macro plots the correlation matrix of the various input
// variables used in TMVA (e.g. running TMVAnalysis.C).  Signal and
// Background are plotted separately

// input: - Input file (result from TMVA),
//        - use of colors or grey scale
//        - use of TMVA plotting TStyle
void TMVA::correlations(TString dataset, TString fin , Bool_t isRegression , 
                        Bool_t /* greyScale */ , Bool_t useTMVAStyle  )
{

   // set style and remove existing canvas'
   TMVAGlob::Initialize( useTMVAStyle );

   // checks if file with name "fin" is already open, and if not opens one
   TFile* file = TMVAGlob::OpenFile( fin );  

   // signal and background or regression problem
   Int_t ncls = (isRegression ? 1 : 2 );
   TString hName[2] = { "CorrelationMatrixS", "CorrelationMatrixB" };
   if (isRegression) hName[0]= "CorrelationMatrix";
   const Int_t width = 600;
   for (Int_t ic=0; ic<ncls; ic++) {

      TH2* h2 = dynamic_cast<TH2*> (file->GetDirectory(dataset.Data())->Get( hName[ic] ));
      if(!h2) {
         cout << "Did not find histogram " << hName[ic] << " in " << fin << endl;
         continue;
      }

      TCanvas* c = new TCanvas( hName[ic], 
                                Form("Correlations between MVA input variables (%s)", 
                                     (isRegression ? "" : (ic==0 ? "signal" : "background"))), 
                                ic*(width+5)+200, 0, width, width ); 
      Float_t newMargin1 = 0.13;
      Float_t newMargin2 = 0.15;
      if (gConfig().fVariablePlotting.fUsePaperStyle) newMargin2 = 0.13;

      c->SetGrid();
      c->SetTicks();
      c->SetLeftMargin  ( newMargin2 );
      c->SetBottomMargin( newMargin2 );
      c->SetRightMargin ( newMargin1 );
      c->SetTopMargin   ( newMargin1 );
      gStyle->SetPalette( 1, 0 );


      gStyle->SetPaintTextFormat( "3g" );

      h2->SetMarkerSize( 1.5 );
      h2->SetMarkerColor( 0 );
      Float_t labelSize = 0.040;
      h2->GetXaxis()->SetLabelSize( labelSize );
      h2->GetYaxis()->SetLabelSize( labelSize );
      h2->LabelsOption( "d" );
      h2->SetLabelOffset( 0.011 );// label offset on x axis    

      h2->Draw("colz"); // color pads   
      c->Update();

      // modify properties of paletteAxis
      TPaletteAxis* paletteAxis = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject( "palette" );
      paletteAxis->SetLabelSize( 0.03 );
      paletteAxis->SetX1NDC( paletteAxis->GetX1NDC() + 0.02 );

      h2->Draw("textsame");  // add text

      // add comment    
      TText* t = new TText( 0.53, 0.88, "Linear correlation coefficients in %" );
      t->SetNDC();
      t->SetTextSize( 0.026 );
      t->AppendPad();    

      // TMVAGlob::plot_logo( );
      c->Update();

      TString fname = dataset+"/plots/";
      fname += hName[ic];
      TMVAGlob::imgconv( c, fname );
   }
}
