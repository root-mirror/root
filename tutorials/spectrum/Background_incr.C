/// \file
/// \ingroup tutorial_spectrum
/// \notebook
/// Example to illustrate the background estimator (class TSpectrum).
///
/// \macro_image
/// \macro_code
///
/// \authors Miroslav Morhac, Olivier Couet

void Background_incr() {
   Int_t i;
   const Int_t nbins = 1024;
   Double_t xmin     = 0;
   Double_t xmax     = nbins;
   Double_t source[nbins];
   gROOT->ForceStyle();

   TH1F *d    = new TH1F("d","",nbins,xmin,xmax);

   TString dir  = gROOT->GetTutorialDir();
   TString file = dir+"/spectrum/TSpectrum.root";
   TFile *f     = new TFile(file.Data());
   TH1F *back = (TH1F*) f->Get("back1");
   back->GetXaxis()->SetRange(1,nbins);
   back->SetTitle("Estimation of background with increasing window");
   back->Draw("L");

   TSpectrum *s = new TSpectrum();

   for (i = 0; i < nbins; i++) source[i] = back->GetBinContent(i + 1);

   // Estimate the background
   s->Background(source,nbins,6,TSpectrum::kBackIncreasingWindow,
                 TSpectrum::kBackOrder2,kFALSE,
                 TSpectrum::kBackSmoothing3,kFALSE);

   // Draw the estimated background
   for (i = 0; i < nbins; i++) d->SetBinContent(i + 1,source[i]);
   d->SetLineColor(kRed);
   d->Draw("SAME L");
}
