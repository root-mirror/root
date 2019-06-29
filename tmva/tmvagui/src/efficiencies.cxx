#include "TMVA/efficiencies.h"

#include "TH2F.h"
#include "TFile.h"
#include "TIterator.h"
#include "TKey.h"

void TMVA::plot_efficiencies(TString dataset, TFile* /*file*/, Int_t type , TDirectory* BinDir)
{
   // input:   - Input file (result from TMVA),
   //          - type = 1 --> plot efficiency(B) versus eff(S)
   //                 = 2 --> plot rejection (B) versus efficiency (S)
   //                 = 3 --> plot 1/eff(B) versus efficiency (S)

   Bool_t __PLOT_LOGO__  = kTRUE;
   Bool_t __SAVE_IMAGE__ = kTRUE;

   // the coordinates
   Float_t x1 = 0;
   Float_t x2 = 1;
   Float_t y1 = 0;
   Float_t y2 = 1;

   // reverse order if "rejection"
   if (type == 2) {
      Float_t z = y1;
      y1 = 1 - y2;
      y2 = 1 - z;
      //      cout << "--- type==2: plot background rejection versus signal efficiency" << endl;
   } else if (type == 3) {
      y1 = 0;
      y2 = -1; // will be set to the max found in the histograms

   } else {
      //  cout << "--- type==1: plot background efficiency versus signal efficiency" << endl;
   }
   // create canvas
   TCanvas* c = new TCanvas( "c", "the canvas", 200, 0, 650, 500 );

   // global style settings
   c->SetGrid();
   c->SetTicks();

   // legend
   Float_t x0L = 0.107,     y0H = 0.899;
   Float_t dxL = 0.457-x0L, dyH = 0.22;
   if (type == 2) {
      x0L = 0.15;
      y0H = 1 - y0H + dyH + 0.07;
   }
   TLegend *legend = new TLegend( x0L, y0H-dyH, x0L+dxL, y0H );
   //legend->SetTextSize( 0.05 );
   legend->SetHeader( "MVA Method:" );
   legend->SetMargin( 0.4 );

   TString xtit = "Signal efficiency";
   TString ytit = "Background efficiency";
   if (type == 2) ytit = "Background rejection";
   if (type == 3) ytit = "1/(Background eff.)";
   TString ftit = ytit + " versus " + xtit;

   TString hNameRef = "effBvsS";
   if (type == 2) hNameRef = "rejBvsS";
   if (type == 3) hNameRef = "invBeffvsSeff";


   if (TString(BinDir->GetName()).Contains("multicut")){
      ftit += "  Bin: ";
      ftit += (BinDir->GetTitle());
   }

   TList xhists;
   TList xmethods;
   UInt_t xnm = TMVAGlob::GetListOfMethods( xmethods ,BinDir);
   if (xnm==0){
      cout << "ups .. no methods found in to plot ROC curve for ... give up"  << endl;
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
            if (hname.Contains( hNameRef ) && hname.BeginsWith( "MVA_" )) {
               if (type==3 && h->GetMaximum() > y2) y2 = h->GetMaximum()*1.1;
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
            if (hname.Contains( hNameRef ) && hname.BeginsWith( "MVA_" )) {
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
      legend->AddEntry(histWithLargestInt,TString(histWithLargestInt->GetTitle()).ReplaceAll("MVA_",""),"l");
      hists.Remove(histWithLargestInt);
   }

   // rescale legend box size
   // current box size has been tuned for 3 MVAs + 1 title
   if (type == 1) {
      dyH *= (1.0 + Float_t(nmva - 3.0)/4.0);
      legend->SetY1( y0H - dyH );
   }
   else {
      dyH *= (Float_t(TMath::Min(10,nmva) - 3.0)/4.0);
      legend->SetY2( y0H + dyH);
   }

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

void TMVA::efficiencies(TString dataset, TString fin , Int_t type , Bool_t useTMVAStyle )
{
   // argument: type = 1 --> plot efficiency(B) versus eff(S)
   //           type = 2 --> plot rejection (B) versus efficiency (S)

   // set style and remove existing canvas'
   TMVAGlob::Initialize( useTMVAStyle );

   // checks if file with name "fin" is already open, and if not opens one
   TFile* file = TMVAGlob::OpenFile( fin );

   plot_efficiencies(dataset, file, type, file->GetDirectory(dataset.Data()));

   return;
}

