/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Data and categories: tools for manipulation of (un)binned datasets
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
using namespace RooFit;

// WVE Add reduction by range

void rf402_datahandling()
{

   // Binned (RooDataHist) and unbinned datasets (RooDataSet) share
   // many properties and inherit from a common abstract base class
   // (RooAbsData), that provides an interface for all operations
   // that can be performed regardless of the data format

   RooRealVar x("x", "x", -10, 10);
   RooRealVar y("y", "y", 0, 40);
   RooCategory c("c", "c");
   c.defineType("Plus", +1);
   c.defineType("Minus", -1);

   // B a s i c   O p e r a t i o n s   o n   u n b i n n e d   d a t a s e t s
   // --------------------------------------------------------------

   // RooDataSet is an unbinned dataset (a collection of points in N-dimensional space)
   RooDataSet d("d", "d", RooArgSet(x, y, c));

   // Unlike RooAbsArgs (RooAbsPdf,RooFormulaVar,....) datasets are not attached to
   // the variables they are constructed from. Instead they are attached to an internal
   // clone of the supplied set of arguments

   // Fill d with dummy values
   Int_t i;
   for (i = 0; i < 1000; i++) {
      x = i / 50 - 10;
      y = sqrt(1.0 * i);
      c.setLabel((i % 2) ? "Plus" : "Minus");

      // We must explicitly refer to x,y,c here to pass the values because
      // d is not linked to them (as explained above)
      d.add(RooArgSet(x, y, c));
   }
   d.Print("v");
   cout << endl;

   // The get() function returns a pointer to the internal copy of the RooArgSet(x,y,c)
   // supplied in the constructor
   const RooArgSet *row = d.get();
   row->Print("v");
   cout << endl;

   // Get with an argument loads a specific data point in row and returns
   // a pointer to row argset. get() always returns the same pointer, unless
   // an invalid row number is specified. In that case a null ptr is returned
   d.get(900)->Print("v");
   cout << endl;

   // R e d u c i n g ,   A p p e n d i n g   a n d   M e r g i n g
   // -------------------------------------------------------------

   // The reduce() function returns a new dataset which is a subset of the original
   cout << endl << ">> d1 has only columns x,c" << endl;
   RooDataSet *d1 = (RooDataSet *)d.reduce(RooArgSet(x, c));
   d1->Print("v");

   cout << endl << ">> d2 has only column y" << endl;
   RooDataSet *d2 = (RooDataSet *)d.reduce(RooArgSet(y));
   d2->Print("v");

   cout << endl << ">> d3 has only the points with y>5.17" << endl;
   RooDataSet *d3 = (RooDataSet *)d.reduce("y>5.17");
   d3->Print("v");

   cout << endl << ">> d4 has only columns x,c for data points with y>5.17" << endl;
   RooDataSet *d4 = (RooDataSet *)d.reduce(RooArgSet(x, c), "y>5.17");
   d4->Print("v");

   // The merge() function adds two data set column-wise
   cout << endl << ">> merge d2(y) with d1(x,c) to form d1(x,c,y)" << endl;
   d1->merge(d2);
   d1->Print("v");

   // The append() function addes two datasets row-wise
   cout << endl << ">> append data points of d3 to d1" << endl;
   d1->append(*d3);
   d1->Print("v");

   // O p e r a t i o n s   o n   b i n n e d   d a t a s e t s
   // ---------------------------------------------------------

   // A binned dataset can be constructed empty, from an unbinned dataset, or
   // from a ROOT native histogram (TH1,2,3)

   cout << ">> construct dh (binned) from d(unbinned) but only take the x and y dimensions," << endl
        << ">> the category 'c' will be projected in the filling process" << endl;

   // The binning of real variables (like x,y) is done using their fit range
   // 'get/setRange()' and number of specified fit bins 'get/setBins()'.
   // Category dimensions of binned datasets get one bin per defined category state
   x.setBins(10);
   y.setBins(10);
   RooDataHist dh("dh", "binned version of d", RooArgSet(x, y), d);
   dh.Print("v");

   RooPlot *yframe = y.frame(Bins(10), Title("Operations on binned datasets"));
   dh.plotOn(yframe); // plot projection of 2D binned data on y

   // Examine the statistics of a binned dataset
   cout << ">> number of bins in dh   : " << dh.numEntries() << endl;
   cout << ">> sum of weights in dh   : " << dh.sum(kFALSE) << endl;
   cout << ">> integral over histogram: " << dh.sum(kTRUE) << endl; // accounts for bin volume

   // Locate a bin from a set of coordinates and retrieve its properties
   x = 0.3;
   y = 20.5;
   cout << ">> retrieving the properties of the bin enclosing coordinate (x,y) = (0.3,20.5) " << endl;
   cout << " bin center:" << endl;
   dh.get(RooArgSet(x, y))->Print("v");         // load bin center coordinates in internal buffer
   cout << " weight = " << dh.weight() << endl; // return weight of last loaded coordinates

   // Reduce the 2-dimensional binned dataset to a 1-dimensional binned dataset
   //
   // All reduce() methods are interfaced in RooAbsData. All reduction techniques
   // demonstrated on unbinned datasets can be applied to binned datasets as well.
   cout << ">> Creating 1-dimensional projection on y of dh for bins with x>0" << endl;
   RooDataHist *dh2 = (RooDataHist *)dh.reduce(y, "x>0");
   dh2->Print("v");

   // Add dh2 to yframe and redraw
   dh2->plotOn(yframe, LineColor(kRed), MarkerColor(kRed));

   // S a v i n g   a n d   l o a d i n g   f r o m   f i l e
   // -------------------------------------------------------

   // Datasets can be persisted with ROOT I/O
   cout << endl << ">> Persisting d via ROOT I/O" << endl;
   TFile f("rf402_datahandling.root", "RECREATE");
   d.Write();
   f.ls();

   // To read back in future session:
   // > TFile f("rf402_datahandling.root") ;
   // > RooDataSet* d = (RooDataSet*) f.FindObject("d") ;

   new TCanvas("rf402_datahandling", "rf402_datahandling", 600, 600);
   gPad->SetLeftMargin(0.15);
   yframe->GetYaxis()->SetTitleOffset(1.4);
   yframe->Draw();
}
