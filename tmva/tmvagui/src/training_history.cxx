#include "TMVA/training_history.h"

#include "TH2F.h"
#include "TFile.h"
#include "TIterator.h"
#include "TKey.h"

void TMVA::plot_training_history(TString dataset, TFile* /*file*/, TDirectory* BinDir)
{

   Bool_t __PLOT_LOGO__  = kTRUE;
   Bool_t __SAVE_IMAGE__ = kTRUE;

   // the coordinates
   Float_t x1 = 999;
   Float_t x2 = -999;
   Float_t y1 = 999.;
   Float_t y2 = -999;

   // create canvas
   TCanvas* c = new TCanvas( "c", "the canvas", 200, 0, 650, 500 );

   // global style settings
   c->SetGrid();
   c->SetTicks();

   // legend
   Float_t x0L = 0.107,     y0H = 0.899;
   Float_t dxL = 0.557-x0L, dyH = 0.22;

   TLegend *legend = new TLegend( x0L, y0H-dyH, x0L+dxL, y0H );
   //legend->SetTextSize( 0.05 );
   legend->SetHeader( "MVA Method:" );
   legend->SetMargin( 0.4 );

   TString xtit = "Training Step";
   TString ytit = "";

   TString ftit = "Training History";

   TString hNameRef = "TrainingHistory";

   TList xhists;
   TList xmethods;
   UInt_t xnm = TMVAGlob::GetListOfMethods( xmethods ,BinDir);
   if (xnm==0){
      cout << "ups .. no methods found in to plot training history for ... give up"  << endl;
      return;
   }
   TIter xnext(&xmethods);
   // loop over all methods
   TKey *xkey;
   while ((xkey = (TKey*)xnext())) {
      TDirectory * mDir = (TDirectory*)xkey->ReadObj();
      TList titles;
      UInt_t ninst = TMVAGlob::GetListOfTitles(mDir,titles);
      if (ninst==0) cout << "hmm... sorry, but this printout was supposed to be only to keep the compiler quite.. never supposed to happen :(" << endl;
      TIter nextTitle(&titles);
      TKey *titkey;
      TDirectory *titDir;
      while ((titkey = TMVAGlob::NextKey(nextTitle,"TDirectory"))) {
         titDir = (TDirectory *)titkey->ReadObj();
         TString methodTitle;
         TMVAGlob::GetMethodTitle(methodTitle,titDir);
         TIter nextKey( titDir->GetListOfKeys() );
         TKey *hkey2;
         while ((hkey2 = TMVAGlob::NextKey(nextKey,"TH1"))) {
            TH1 *h = (TH1*)hkey2->ReadObj();
            TString hname = h->GetName();
            if (hname.Contains( hNameRef ) && hname.BeginsWith( "TrainingHistory_" )) {

               if (h->GetMaximum() > y2) y2 = h->GetMaximum()*1.1;
               if (h->GetMinimum() < y1) y1 = h->GetMinimum();
               if (h->GetBinLowEdge(0) < x1 ) x1 = h->GetBinLowEdge(0);
               if (h->GetBinLowEdge(h->GetNbinsX()+1 ) > x2 ) x2 = h->GetBinLowEdge(h->GetNbinsX()+1 );

            }
         }
      }
   }


   // draw empty frame
   if(gROOT->FindObject("frame")!=0) gROOT->FindObject("frame")->Delete();
   TH2F* frame = new TH2F( "frame", ftit, 500, x1, x2, 500, y1, y2 );
   frame->GetXaxis()->SetTitle( xtit );
   frame->GetYaxis()->SetTitle( ytit );
   TMVAGlob::SetFrameStyle( frame, 1.0 );

   frame->Draw();

   Int_t color = 1;
   Int_t nmva  = 0;
   TKey *key;

   TList hists;
   TList methods;
   UInt_t nm = TMVAGlob::GetListOfMethods( methods,BinDir );
   if (nm==0){
      cout << "ups .. no methods found in to plot ROC curve for ... give up"  << endl;
      return;
   }
   //   TIter next(file->GetListOfKeys());
   TIter next(&methods);

   // loop over all methods
   while ((key = (TKey*)next())) {
      TDirectory * mDir = (TDirectory*)key->ReadObj();
      TList titles;
      UInt_t ninst = TMVAGlob::GetListOfTitles(mDir,titles);
      if (ninst==0) cout << "hmm...  sorry, but this printout was supposed to be only to keep the compiler quite.. never supposed to happen :(" << endl;
      TIter nextTitle(&titles);
      TKey *titkey;
      TDirectory *titDir;
      while ((titkey = TMVAGlob::NextKey(nextTitle,"TDirectory"))) {
         titDir = (TDirectory *)titkey->ReadObj();
         TString methodTitle;
         TMVAGlob::GetMethodTitle(methodTitle,titDir);
         TIter nextKey( titDir->GetListOfKeys() );
         TKey *hkey2;
         while ((hkey2 = TMVAGlob::NextKey(nextKey,"TH1"))) {
            TH1 *h = (TH1*)hkey2->ReadObj();
            TString hname = h->GetName();
            if (hname.Contains( hNameRef ) && hname.BeginsWith( "TrainingHistory_" )) {
               h->SetLineWidth(3);
               h->SetLineColor(color);
               color++; if (color == 5 || color == 10 || color == 11) color++;
               h->Draw("csame");
               hists.Add(h);
               nmva++;
            }
         }
      }
   }

   while (hists.GetSize()) {
      TListIter hIt(&hists);
      TH1* hist(0);
      Double_t largestInt=-1;
      TH1* histWithLargestInt(0);
      while ((hist = (TH1*)hIt())!=0) {
         Double_t integral = hist->Integral(1,hist->FindBin(0.9999));
         if (integral>largestInt) {
            largestInt = integral;
            histWithLargestInt = hist;
         }
      }
      if (histWithLargestInt == 0) {
         cout << "ERROR - unknown hist \"histWithLargestInt\" --> serious problem in ROOT file" << endl;
         break;
      }
      TString histLabel = TString(histWithLargestInt->GetTitle()).ReplaceAll("MVA_", "");
      histLabel.ReplaceAll("TrainingHistory_", "");
      legend->AddEntry(histWithLargestInt, histLabel, "l");
      hists.Remove(histWithLargestInt);
   }

   // rescale legend box size
   // current box size has been tuned for 3 MVAs + 1 title
   dyH *= (1. + (Float_t(TMath::Min(10,nmva) - 3.0)/4.0) );
   legend->SetY1( y0H - dyH);

   // redraw axes
   frame->Draw("sameaxis");
   legend->Draw("same");

   // ============================================================

   if (__PLOT_LOGO__) TMVAGlob::plot_logo();

   // ============================================================

   c->Update();

   TString fname = dataset+"/plots/" + hNameRef;
   if (TString(BinDir->GetName()).Contains("multicut")){
      TString fprepend(BinDir->GetName());
      fprepend.ReplaceAll("multicutMVA_","");
      fname = dataset+"plots/" + fprepend + "_" + hNameRef;
   }
   if (__SAVE_IMAGE__) TMVAGlob::imgconv( c, fname );

   return;
}

void TMVA::training_history(TString dataset, TString fin , Bool_t useTMVAStyle )
{

   // set style and remove existing canvas'
   TMVAGlob::Initialize( useTMVAStyle );

   // checks if file with name "fin" is already open, and if not opens one
   TFile* file = TMVAGlob::OpenFile( fin );

   plot_training_history(dataset, file, file->GetDirectory(dataset.Data()));

   return;
}
