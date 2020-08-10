/// \file
/// \ingroup tutorial_roostats
/// \notebook -js
/// \brief High Level Factory: creation of a simple model
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Danilo Piparo

#include <fstream>
#include "TString.h"
#include "TROOT.h"
#include "RooGlobalFunc.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooStats/HLFactory.h"

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using namespace std;

void rs601_HLFactoryexample()
{

   // --- Build the datacard and dump to file---

   TString card_name("HLFavtoryexample.rs");
   ofstream ofile(card_name);
   ofile << "// The simplest card\n\n"
         << "gauss = Gaussian(mes[5.20,5.30],mean[5.28,5.2,5.3],width[0.0027,0.001,1]);\n"
         << "argus = ArgusBG(mes,5.291,argpar[-20,-100,-1]);\n"
         << "sum = SUM(nsig[200,0,10000]*gauss,nbkg[800,0,10000]*argus);\n\n";

   ofile.close();

   HLFactory hlf("HLFavtoryexample", card_name, false);

   // --- Take elements out of the internal workspace ---

   auto w = hlf.GetWs();

   auto mes = dynamic_cast<RooRealVar *>(w->arg("mes"));
   auto sum = dynamic_cast<RooAbsPdf *>(w->pdf("sum"));
   auto argus = dynamic_cast<RooAbsPdf *>(w->pdf("argus"));

   // --- Generate a toyMC sample from composite PDF ---
   auto data = sum->generate(*mes, 2000);

   // --- Perform extended ML fit of composite PDF to toy data ---
   sum->fitTo(*data);

   // --- Plot toy data and composite PDF overlaid ---
   auto mesframe = mes->frame();
   data->plotOn(mesframe);
   sum->plotOn(mesframe);
   sum->plotOn(mesframe, Components(*argus), LineStyle(kDashed));

   gROOT->SetStyle("Plain");
   mesframe->Draw();
}
