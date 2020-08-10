/// \file
/// \ingroup tutorial_spectrum
/// \notebook
/// Example to illustrate deconvolution function (class TSpectrum).
///
/// \macro_image
/// \macro_code
///
/// \authors Miroslav Morhac, Olivier Couet

void Deconvolution_wide() {
   Int_t i;
   const Int_t nbins = 256;
   Double_t xmin     = 0;
   Double_t xmax     = nbins;
   Double_t source[nbins];
   Double_t response[nbins];
   gROOT->ForceStyle();

   TString dir  = gROOT->GetTutorialDir();
   TString file = dir+"/spectrum/TSpectrum.root";
   TFile *f     = new TFile(file.Data());
   TH1F *h = (TH1F*) f->Get("decon3");
   h->SetTitle("Deconvolution of closely positioned overlapping peaks using Gold deconvolution method");
   TH1F *d = (TH1F*) f->Get("decon_response_wide");

   for (i = 0; i < nbins; i++) source[i]   = h->GetBinContent(i + 1);
   for (i = 0; i < nbins; i++) response[i] = d->GetBinContent(i + 1);

   h->SetMaximum(50000);
   h->Draw("L");
   TSpectrum *s = new TSpectrum();
   s->Deconvolution(source,response,256,10000,1,1);

   for (i = 0; i < nbins; i++) d->SetBinContent(i + 1,source[i]);
   d->SetLineColor(kRed);
   d->Draw("SAME L");
}
