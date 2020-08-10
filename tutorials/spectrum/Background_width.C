/// \file
/// \ingroup tutorial_spectrum
/// \notebook
/// Example to illustrate the influence of the clipping window width on the
/// estimated background.
///
/// \macro_image
/// \macro_code
///
/// \authors Miroslav Morhac, Olivier Couet

void Background_width() {
   Int_t i;
   const Int_t nbins = 1024;
   Double_t xmin     = 0;
   Double_t xmax     = nbins;
   Double_t source[nbins];
   gROOT->ForceStyle();

   TString dir  = gROOT->GetTutorialDir();
   TString file = dir+"/spectrum/TSpectrum.root";
   TFile *f     = new TFile(file.Data());
   TH1F *back = (TH1F*) f->Get("back1");
   TH1F *d1 = new TH1F("d1","",nbins,xmin,xmax);
   TH1F *d2 = new TH1F("d2","",nbins,xmin,xmax);
   TH1F *d3 = new TH1F("d3","",nbins,xmin,xmax);

   back->GetXaxis()->SetRange(1,nbins);
   back->SetTitle("Influence of clipping window width on the estimated background");
   back->Draw("L");

   TSpectrum *s = new TSpectrum();

   for (i = 0; i < nbins; i++) source[i]=back->GetBinContent(i + 1);
   s->Background(source,nbins,4,TSpectrum::kBackDecreasingWindow,
                 TSpectrum::kBackOrder2,kFALSE,
                 TSpectrum::kBackSmoothing3,kFALSE);
   for (i = 0; i < nbins; i++) d1->SetBinContent(i + 1,source[i]);
   d1->SetLineColor(kRed);
   d1->Draw("SAME L");

   for (i = 0; i < nbins; i++) source[i]=back->GetBinContent(i + 1);
   s->Background(source,nbins,6,TSpectrum::kBackDecreasingWindow,
                 TSpectrum::kBackOrder2,kFALSE,
                 TSpectrum::kBackSmoothing3,kFALSE);
   for (i = 0; i < nbins; i++) d2->SetBinContent(i + 1,source[i]);
   d2->SetLineColor(kOrange);
   d2->Draw("SAME L");

   for (i = 0; i < nbins; i++) source[i]=back->GetBinContent(i + 1);
   s->Background(source,nbins,8,TSpectrum::kBackDecreasingWindow,
                 TSpectrum::kBackOrder2,kFALSE,
                 TSpectrum::kBackSmoothing3,kFALSE);
   for (i = 0; i < nbins; i++) d3->SetBinContent(i + 1,source[i]);
   d3->SetLineColor(kGreen);
   d3->Draw("SAME L");
}
