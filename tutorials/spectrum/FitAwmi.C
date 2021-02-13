/// \file
/// \ingroup tutorial_spectrum
/// \notebook
/// This macro fits the source spectrum using the AWMI algorithm
/// from the "TSpectrumFit" class ("TSpectrum" class is used to find peaks).
///
/// To try this macro, in a ROOT (5 or 6) prompt, do:
///
/// ~~~{.cpp}
///  root > .x FitAwmi.C
/// ~~~
///
/// or:
///
/// ~~~{.cpp}
///  root > .x FitAwmi.C++
///  root > FitAwmi(); // re-run with another random set of peaks
/// ~~~
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author

#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TSpectrumFit.h"
#include "TPolyMarker.h"
#include "TList.h"

#include <iostream>

TH1F *FitAwmi_Create_Spectrum(void) {
   Int_t nbins = 1000;
   Double_t xmin = -10., xmax = 10.;
   delete gROOT->FindObject("h"); // prevent "memory leak"
   TH1F *h = new TH1F("h", "simulated spectrum", nbins, xmin, xmax);
   h->SetStats(kFALSE);
   TF1 f("f", "TMath::Gaus(x, [0], [1], 1)", xmin, xmax);
   // f.SetParNames("mean", "sigma");
   gRandom->SetSeed(0); // make it really random
   // create well separated peaks with exactly known means and areas
   // note: TSpectrumFit assumes that all peaks have the same sigma
   Double_t sigma = (xmax - xmin) / Double_t(nbins) * Int_t(gRandom->Uniform(2., 6.));
   Int_t npeaks = 0;
   while (xmax > (xmin + 6. * sigma)) {
      npeaks++;
      xmin += 3. * sigma; // "mean"
      f.SetParameters(xmin, sigma);
      Double_t area = 1. * Int_t(gRandom->Uniform(1., 11.));
      h->Add(&f, area, ""); // "" ... or ... "I"
      std::cout << "created "
                << xmin << " "
                << (area / sigma / TMath::Sqrt(TMath::TwoPi())) << " "
                << area << std::endl;
      xmin += 3. * sigma;
   }
   std::cout << "the total number of created peaks = " << npeaks
             << " with sigma = " << sigma << std::endl;
   return h;
}

void FitAwmi(void) {

   TH1F *h = FitAwmi_Create_Spectrum();

   TCanvas *cFit = ((TCanvas *)(gROOT->GetListOfCanvases()->FindObject("cFit")));
   if (!cFit) cFit = new TCanvas("cFit", "cFit", 10, 10, 1000, 700);
   else cFit->Clear();
   h->Draw("L");
   Int_t i, nfound, bin;
   Int_t nbins = h->GetNbinsX();

   Double_t *source = new Double_t[nbins];
   Double_t *dest = new Double_t[nbins];

   for (i = 0; i < nbins; i++) source[i] = h->GetBinContent(i + 1);
   TSpectrum *s = new TSpectrum(); // note: default maxpositions = 100
   // searching for candidate peaks positions
   nfound = s->SearchHighRes(source, dest, nbins, 2., 2., kFALSE, 10000, kFALSE, 0);
   // filling in the initial estimates of the input parameters
   Bool_t *FixPos = new Bool_t[nfound];
   Bool_t *FixAmp = new Bool_t[nfound];
   for(i = 0; i < nfound; i++) FixAmp[i] = FixPos[i] = kFALSE;

   Double_t *Pos, *Amp = new Double_t[nfound]; // ROOT 6

   Pos = s->GetPositionX(); // 0 ... (nbins - 1)
   for (i = 0; i < nfound; i++) {
      bin = 1 + Int_t(Pos[i] + 0.5); // the "nearest" bin
      Amp[i] = h->GetBinContent(bin);
   }
   TSpectrumFit *pfit = new TSpectrumFit(nfound);
   pfit->SetFitParameters(0, (nbins - 1), 1000, 0.1, pfit->kFitOptimChiCounts,
                          pfit->kFitAlphaHalving, pfit->kFitPower2,
                          pfit->kFitTaylorOrderFirst);
   pfit->SetPeakParameters(2., kFALSE, Pos, FixPos, Amp, FixAmp);
   // pfit->SetBackgroundParameters(source[0], kFALSE, 0., kFALSE, 0., kFALSE);
   pfit->FitAwmi(source);
   Double_t *Positions = pfit->GetPositions();
   Double_t *PositionsErrors = pfit->GetPositionsErrors();
   Double_t *Amplitudes = pfit->GetAmplitudes();
   Double_t *AmplitudesErrors = pfit->GetAmplitudesErrors();
   Double_t *Areas = pfit->GetAreas();
   Double_t *AreasErrors = pfit->GetAreasErrors();
   delete gROOT->FindObject("d"); // prevent "memory leak"
   TH1F *d = new TH1F(*h); d->SetNameTitle("d", ""); d->Reset("M");
   for (i = 0; i < nbins; i++) d->SetBinContent(i + 1, source[i]);
   Double_t x1 = d->GetBinCenter(1), dx = d->GetBinWidth(1);
   Double_t sigma, sigmaErr;
   pfit->GetSigma(sigma, sigmaErr);

   // current TSpectrumFit needs a sqrt(2) correction factor for sigma
   sigma /= TMath::Sqrt2(); sigmaErr /= TMath::Sqrt2();
   // convert "bin numbers" into "x-axis values"
   sigma *= dx; sigmaErr *= dx;

   std::cout << "the total number of found peaks = " << nfound
             << " with sigma = " << sigma << " (+-" << sigmaErr << ")"
             << std::endl;
   std::cout << "fit chi^2 = " << pfit->GetChi() << std::endl;
   for (i = 0; i < nfound; i++) {
      bin = 1 + Int_t(Positions[i] + 0.5); // the "nearest" bin
      Pos[i] = d->GetBinCenter(bin);
      Amp[i] = d->GetBinContent(bin);

      // convert "bin numbers" into "x-axis values"
      Positions[i] = x1 + Positions[i] * dx;
      PositionsErrors[i] *= dx;
      Areas[i] *= dx;
      AreasErrors[i] *= dx;

      std::cout << "found "
                << Positions[i] << " (+-" << PositionsErrors[i] << ") "
                << Amplitudes[i] << " (+-" << AmplitudesErrors[i] << ") "
                << Areas[i] << " (+-" << AreasErrors[i] << ")"
                << std::endl;
   }
   d->SetLineColor(kRed); d->SetLineWidth(1);
   d->Draw("SAME L");
   TPolyMarker *pm = ((TPolyMarker*)(h->GetListOfFunctions()->FindObject("TPolyMarker")));
   if (pm) {
      h->GetListOfFunctions()->Remove(pm);
      delete pm;
   }
   pm = new TPolyMarker(nfound, Pos, Amp);
   h->GetListOfFunctions()->Add(pm);
   pm->SetMarkerStyle(23);
   pm->SetMarkerColor(kRed);
   pm->SetMarkerSize(1);
   // cleanup
   delete pfit;
   delete [] Amp;
   delete [] FixAmp;
   delete [] FixPos;
   delete s;
   delete [] dest;
   delete [] source;
   return;
}
