/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Data and categories: advanced options for importing data from ROOT TTree and THx histograms
///
/// Basic import options are demonstrated in rf102_dataimport.C
///
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TH1.h"
#include "TTree.h"
#include "TRandom.h"
#include <map>

using namespace RooFit;

TH1 *makeTH1(const char *name, Double_t mean, Double_t sigma);
TTree *makeTTree();

void rf401_importttreethx()
{
   // I m p o r t  m u l t i p l e   T H 1   i n t o   a   R o o D a t a H i s t
   // --------------------------------------------------------------------------

   // Create thee ROOT TH1 histograms
   TH1 *hh_1 = makeTH1("hh1", 0, 3);
   TH1 *hh_2 = makeTH1("hh2", -3, 1);
   TH1 *hh_3 = makeTH1("hh3", +3, 4);

   // Declare observable x
   RooRealVar x("x", "x", -10, 10);

   // Create category observable c that serves as index for the ROOT histograms
   RooCategory c("c", "c", {{"SampleA",0}, {"SampleB",1}, {"SampleC",2}});

   // Create a binned dataset that imports contents of all TH1 mapped by index category c
   RooDataHist *dh = new RooDataHist("dh", "dh", x, Index(c), Import("SampleA", *hh_1), Import("SampleB", *hh_2),
                                     Import("SampleC", *hh_3));
   dh->Print();

   // Alternative constructor form for importing multiple histograms
   std::map<std::string, TH1 *> hmap;
   hmap["SampleA"] = hh_1;
   hmap["SampleB"] = hh_2;
   hmap["SampleC"] = hh_3;
   RooDataHist *dh2 = new RooDataHist("dh", "dh", x, c, hmap);
   dh2->Print();

   // I m p o r t i n g   a   T T r e e   i n t o   a   R o o D a t a S e t   w i t h   c u t s
   // -----------------------------------------------------------------------------------------

   TTree *tree = makeTTree();

   // Define observables y,z
   RooRealVar y("y", "y", -10, 10);
   RooRealVar z("z", "z", -10, 10);

   // Import only observables (y,z)
   RooDataSet ds("ds", "ds", RooArgSet(x, y), Import(*tree));
   ds.Print();

   // Import observables (x,y,z) but only event for which (y+z<0) is true
   RooDataSet ds2("ds2", "ds2", RooArgSet(x, y, z), Import(*tree), Cut("y+z<0"));
   ds2.Print();

   // I m p o r t i n g   i n t e g e r   T T r e e   b r a n c h e s
   // ---------------------------------------------------------------

   // Import integer tree branch as RooRealVar
   RooRealVar i("i", "i", 0, 5);
   RooDataSet ds3("ds3", "ds3", RooArgSet(i, x), Import(*tree));
   ds3.Print();

   // Define category i
   RooCategory icat("i", "i");
   icat.defineType("State0", 0);
   icat.defineType("State1", 1);

   // Import integer tree branch as RooCategory (only events with i==0 and i==1
   // will be imported as those are the only defined states)
   RooDataSet ds4("ds4", "ds4", RooArgSet(icat, x), Import(*tree));
   ds4.Print();

   // I m p o r t  m u l t i p l e   R o o D a t a S e t s   i n t o   a   R o o D a t a S e t
   // ----------------------------------------------------------------------------------------

   // Create three RooDataSets in (y,z)
   RooDataSet *dsA = (RooDataSet *)ds2.reduce(RooArgSet(x, y), "z<-5");
   RooDataSet *dsB = (RooDataSet *)ds2.reduce(RooArgSet(x, y), "abs(z)<5");
   RooDataSet *dsC = (RooDataSet *)ds2.reduce(RooArgSet(x, y), "z>5");

   // Create a dataset that imports contents of all the above datasets mapped by index category c
   RooDataSet *dsABC = new RooDataSet("dsABC", "dsABC", RooArgSet(x, y), Index(c), Import("SampleA", *dsA),
                                      Import("SampleB", *dsB), Import("SampleC", *dsC));

   dsABC->Print();
}

TH1 *makeTH1(const char *name, Double_t mean, Double_t sigma)
{
   // Create ROOT TH1 filled with a Gaussian distribution

   TH1D *hh = new TH1D(name, name, 100, -10, 10);
   for (int i = 0; i < 1000; i++) {
      hh->Fill(gRandom->Gaus(mean, sigma));
   }
   return hh;
}

TTree *makeTTree()
{
   // Create ROOT TTree filled with a Gaussian distribution in x and a uniform distribution in y

   TTree *tree = new TTree("tree", "tree");
   Double_t *px = new Double_t;
   Double_t *py = new Double_t;
   Double_t *pz = new Double_t;
   Int_t *pi = new Int_t;
   tree->Branch("x", px, "x/D");
   tree->Branch("y", py, "y/D");
   tree->Branch("z", pz, "z/D");
   tree->Branch("i", pi, "i/I");
   for (int i = 0; i < 100; i++) {
      *px = gRandom->Gaus(0, 3);
      *py = gRandom->Uniform() * 30 - 15;
      *pz = gRandom->Gaus(0, 5);
      *pi = i % 3;
      tree->Fill();
   }
   return tree;
}
