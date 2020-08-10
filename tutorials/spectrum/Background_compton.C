/// \file
/// \ingroup tutorial_spectrum
/// \notebook
/// Example to illustrate the background estimator (class TSpectrum) including
/// Compton edges.
///
/// \macro_image
/// \macro_code
///
/// \authors Miroslav Morhac, Olivier Couet

void Background_compton() {
   Int_t i;
   const Int_t nbins = 256;
   Double_t xmin     = 0;
   Double_t xmax     = nbins;
   Double_t source[nbins];
   gROOT->ForceStyle();

   TH1F *d1 = new TH1F("d1","",nbins,xmin,xmax);

   TString dir  = gROOT->GetTutorialDir();
   TString file = dir+"/spectrum/TSpectrum.root";
   TFile *f     = new TFile(file.Data());
   TH1F *back = (TH1F*) f->Get("back3");
   back->SetTitle("Estimation of background with Compton edges under peaks");
   back->GetXaxis()->SetRange(1,nbins);
   back->Draw("L");

   TSpectrum *s = new TSpectrum();

   for (i = 0; i < nbins; i++) source[i]=back->GetBinContent(i + 1);
   s->Background(source,nbins,10,TSpectrum::kBackDecreasingWindow,
                 TSpectrum::kBackOrder8,kTRUE,
                 TSpectrum::kBackSmoothing5,kTRUE);
   for (i = 0; i < nbins; i++) d1->SetBinContent(i + 1,source[i]);
   d1->SetLineColor(kRed);
   d1->Draw("SAME L");
}
