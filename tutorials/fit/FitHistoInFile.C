/// \file
/// \notebook
/// \brief Example for fitting a signal + background model to a histogram found in a file.
///
/// This example can be executed as:
///     root > .x FitHistoInFile.C
///
/// \macro_image
/// \macro_code
/// \macro_output
/// \author Author E. von Toerne
/// Based on FittingDemo.C by Rene Brun

#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"

// Function parameters are passed as an array to TF1. Here, we
// define the position of each parameter in this array.
// Note: N_PAR will give us the total number of parameters. Make
// sure it is always the last entry!
enum ParIndex_t {
   Bkg0=0, Bkg1=1, Bkg2,
   SigScale, SigSigma, SigMean,
   N_PAR};
// Use this map to (re-)name parameters for the plot
const std::map<ParIndex_t,std::string> parNames{
   {Bkg0, "Bkg0"}, {Bkg1, "Bkg1"}, {Bkg2, "Bkg2"},
   {SigScale, "Gauss scale"}, {SigSigma, "Gauss #sigma"}, {SigMean, "Gauss #mu"}
};


// Quadratic background function
Double_t background(Double_t *x, Double_t *par) {
   return par[Bkg0] + par[Bkg1]*x[0] + par[Bkg2]*x[0]*x[0];
}

// Gauss Peak function
Double_t signal(Double_t *x, Double_t *par) {
   return par[SigScale]*TMath::Gaus(x[0], par[SigMean], par[SigSigma], true);
}

// Sum of background and peak function. We pass x and the fit parameters
// down to the signal and background functions.
Double_t fitFunction(Double_t *x, Double_t *par) {
   return background(x, par) + signal(x, par);
}

// Fit "fitFunction" to the histogram, and draw results on the canvas `c1`.
void FitRoutine(TCanvas* c1, TH1* histo, float fitxmin, float fitxmax, TString filename){
   c1->cd();
   // create a TF1 with the range from 0 to 3 and N_PAR parameters (six by default)
   TF1 fitFcn("fitFcn",fitFunction,fitxmin,fitxmax,N_PAR);
   fitFcn.SetNpx(500);
   fitFcn.SetLineWidth(2);
   fitFcn.SetLineColor(kBlue);

   // Assign the names from the map "parNames". Optional, but makes a nicer plot.
   for (auto& idx_name : parNames) {
     fitFcn.SetParName(idx_name.first, idx_name.second.c_str());
   }

   // Fit. First set ok-ish starting values for the parameters
   fitFcn.SetParameters(30,0,0,50.,0.1,1.);
   histo->GetXaxis()->SetRange(2,40);
   histo->Fit("fitFcn","VR+","ep");

   // improve the picture:
   // Draw signal and background functions separately
   TF1 backFcn("backFcn",background,fitxmin,fitxmax,N_PAR);
   backFcn.SetLineColor(kRed);
   TF1 signalFcn("signalFcn",signal,fitxmin,fitxmax,N_PAR);
   signalFcn.SetLineColor(kBlue);
   signalFcn.SetNpx(500);

   // Retrieve fit parameters, and copy them to the signal and background functions
   Double_t par[N_PAR];
   fitFcn.GetParameters(par);

   backFcn.SetParameters(par);
   backFcn.DrawCopy("same");

   signalFcn.SetParameters(par);
   signalFcn.SetLineColor(kGreen);
   signalFcn.DrawCopy("same");

   const double binwidth = histo->GetXaxis()->GetBinWidth(1);
   const double integral = signalFcn.Integral(0.,3.);
   cout << "number of signal events = " << integral/binwidth << " " << binwidth<< endl;

   // draw the legend
   TLegend legend(0.15,0.7,0.28,0.85);
   legend.SetTextFont(72);
   legend.SetTextSize(0.03);
   legend.AddEntry(histo,"Data","lpe");
   legend.AddEntry(&backFcn,"Bgd","l");
   legend.AddEntry(&signalFcn,"Sig","l");
   legend.AddEntry(&fitFcn,"Sig+Bgd","l");
   legend.DrawClone();
   histo->Draw("esame");
   c1->SaveAs(filename);
}

// Create a file with example data
void CreateRootFile(){
   // The data in array form
   const int nBins = 60;
   Double_t data[nBins] = { 6, 1,10,12, 6,13,23,22,15,21,
                           23,26,36,25,27,35,40,44,66,81,
                           75,57,43,37,36,31,35,36,43,32,
                           40,37,38,33,36,44,42,37,32,32,
                           43,44,35,33,33,39,29,41,32,44,
                           26,39,29,35,32,21,21,15,25,15};

   TFile* f = new TFile("exampleRootFile.root","RECREATE");
   TH1D *histo = new TH1D("histo", "Gauss Peak on Quadratic Background;x;Events/0.05",60,0,3);
   for(int i=0; i < nBins;  i++) histo->SetBinContent(i+1,data[i]);
   f->Write();
   f->Close();
}

// Show how to fit a histogram that's in a file.
void FitHistoInFile() {
   gStyle->SetOptFit(1111);
   gStyle->SetOptStat(0);

   // fit range from fitxmin to fitxmax
   float fitxmin=0.2;
   float fitxmax=2.7;

   TCanvas *c1 = new TCanvas("c1","Fitting Demo of Histogram in File",10,10,700,500);
   CreateRootFile();
   TFile* f = new TFile("exampleRootFile.root");
   TH1D* histo= nullptr;
   f->GetObject("histo",histo);
   if (!histo){
      cout << "histo not found"<<endl;
      return;
   }
   histo->SetMarkerStyle(21);
   histo->SetMarkerSize(0.8);
   // now call the fit routine
   FitRoutine(c1,histo, fitxmin, fitxmax,"FitHistoInFile.pdf");
}
