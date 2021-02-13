// @(#)root/hist:$Id$
// Author: Rene Brun   12/10/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TEnv.h"
#include "TBrowser.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TVirtualPad.h"
#include "TVirtualFitter.h"
#include "TPluginManager.h"
#include "TMath.h"
#include "TF1.h"
#include "strlcpy.h"

#include "HFitInterface.h"
#include "Fit/DataRange.h"
#include "Math/MinimizerOptions.h"

#include <iostream>
#include <cstdlib>
#include <cctype>

extern void H1LeastSquareSeqnd(Int_t n, Double_t *a, Int_t idim, Int_t &ifail, Int_t k, Double_t *b);

ClassImp(TMultiGraph);


////////////////////////////////////////////////////////////////////////////////

/** \class TMultiGraph
    \ingroup Hist
     \brief A TMultiGraph is a collection of TGraph (or derived) objects.

- [Introduction](#MG00)
- [MultiGraphs' drawing](#MG01)
    - [Setting drawing options](#MG01a)
    - [Titles setting](#MG01b)
    - [The option "3D"](#MG01c)
    - [Legend drawing](#MG01d)
    - [Automatic coloring](#MG01e)
    - [Reverse axis](#MG01f)
- [MultiGraphs' fitting](#MG02)
    - [Fit box position](#MG02a)
- [Axis' limits setting](#MG03)


### <a name="MG00"></a> Introduction

A TMultiGraph allows to manipulate a set of graphs as a single entity. In particular,
when drawn, the X and Y axis ranges are automatically computed such as all the graphs
will be visible.

`TMultiGraph::Add` should be used to add a new graph to the list.

The TMultiGraph owns the objects in the list.

The number of graphs in a multigraph can be retrieve with:
~~~ {.cpp}
mg->GetListOfGraphs()->GetEntries();
~~~

### <a name="MG00"></a> MultiGraphs' Drawing

The drawing options are the same as for TGraph.
Like for TGraph, the painting is performed thanks to the TGraphPainter
class. All details about the various painting options are given in this class.

Example:
~~~ {.cpp}
     TGraph *gr1 = new TGraph(...
     TGraphErrors *gr2 = new TGraphErrors(...
     TMultiGraph *mg = new TMultiGraph();
     mg->Add(gr1,"lp");
     mg->Add(gr2,"cp");
     mg->Draw("a");
~~~

#### <a name="MG01a"></a> Setting drawing options

The drawing option for each TGraph may be specified as an optional
second argument of the `Add` function.

If a draw option is specified, it will be used to draw the graph,
otherwise the graph will be drawn with the option specified in
`TMultiGraph::Draw`

#### <a name="MG01b"></a> Titles setting

The global title and the axis titles can be modified the following way:

~~~ {.cpp}
   [...]
   auto mg = new TMultiGraph;
   mg->SetTitle("title;xaxis title; yaxis title");
   mg->Add(g1);
   mg->Add(g2);
   mg->Draw("apl");
~~~

#### <a name="MG01c"></a> The option "3D"

A special option `3D` allows to draw the graphs in a 3D space. See the
following example:

Begin_Macro(source)
{
   auto c0 = new TCanvas("c1","multigraph L3",200,10,700,500);

   auto mg = new TMultiGraph();

   auto gr1 = new TGraph(); gr1->SetLineColor(kBlue);
   auto gr2 = new TGraph(); gr2->SetLineColor(kRed);
   auto gr3 = new TGraph(); gr3->SetLineColor(kGreen);
   auto gr4 = new TGraph(); gr4->SetLineColor(kOrange);

   Double_t dx = 6.28/1000;
   Double_t x  = -3.14;

   for (int i=0; i<=1000; i++) {
      x = x+dx;
      gr1->SetPoint(i,x,2.*TMath::Sin(x));
      gr2->SetPoint(i,x,TMath::Cos(x));
      gr3->SetPoint(i,x,TMath::Cos(x*x));
      gr4->SetPoint(i,x,TMath::Cos(x*x*x));
   }

   mg->Add(gr4); gr4->SetTitle("Cos(x*x*x)"); gr4->SetLineWidth(3);
   mg->Add(gr3); gr3->SetTitle("Cos(x*x)")  ; gr3->SetLineWidth(3);
   mg->Add(gr2); gr2->SetTitle("Cos(x)")    ; gr2->SetLineWidth(3);
   mg->Add(gr1); gr1->SetTitle("2*Sin(x)")  ; gr1->SetLineWidth(3);

   mg->SetTitle("Multi-graph Title; X-axis Title; Y-axis Title");

   mg->Draw("a fb l3d");

   mg->GetHistogram()->GetXaxis()->SetRangeUser(0.,2.5);
   gPad->Modified();
   gPad->Update();
}
End_Macro

#### <a name="MG01d"></a> Legend drawing

The method TPad::BuildLegend is able to extract the graphs inside a
multigraph. The following example demonstrate this.

Begin_Macro(source)
{
   auto c3 = new TCanvas("c3","c3",600, 400);

   auto mg = new TMultiGraph("mg","mg");

   const Int_t size = 10;

   double px[size];
   double py1[size];
   double py2[size];
   double py3[size];

   for ( int i = 0; i <  size ; ++i ) {
      px[i] = i;
      py1[i] = size - i;
      py2[i] = size - 0.5 * i;
      py3[i] = size - 0.6 * i;
   }

   auto gr1 = new TGraph( size, px, py1 );
   gr1->SetName("gr1");
   gr1->SetTitle("graph 1");
   gr1->SetMarkerStyle(21);
   gr1->SetDrawOption("AP");
   gr1->SetLineColor(2);
   gr1->SetLineWidth(4);
   gr1->SetFillStyle(0);

   auto gr2 = new TGraph( size, px, py2 );
   gr2->SetName("gr2");
   gr2->SetTitle("graph 2");
   gr2->SetMarkerStyle(22);
   gr2->SetMarkerColor(2);
   gr2->SetDrawOption("P");
   gr2->SetLineColor(3);
   gr2->SetLineWidth(4);
   gr2->SetFillStyle(0);

   auto gr3 = new TGraph( size, px, py3 );
   gr3->SetName("gr3");
   gr3->SetTitle("graph 3");
   gr3->SetMarkerStyle(23);
   gr3->SetLineColor(4);
   gr3->SetLineWidth(4);
   gr3->SetFillStyle(0);

   mg->Add( gr1 );
   mg->Add( gr2 );

   gr3->Draw("ALP");
   mg->Draw("LP");
   c3->BuildLegend();
}
End_Macro

#### <a name="MG01e"></a> Automatic coloring

Automatic coloring according to the current palette is available as shown in the
following example:

Begin_Macro(source)
../../../tutorials/graphs/multigraphpalettecolor.C
End_Macro

#### <a name="MG01f"></a> Reverse axis

\since **ROOT version 6.19/02**

When a TMultiGraph is drawn, the X-axis is drawn with increasing values from left to
right and the Y-axis from bottom to top. The two options RX and RY allow to change
this order. The option RX allows to draw the X-axis with increasing values from
right to left and the RY option allows to draw the Y-axis with increasing values
from top to bottom. The following example illustrate how to use these options.

Begin_Macro(source)
{
   auto *c = new TCanvas();
   c->Divide(2,1);

   auto *g1 = new TGraphErrors();
   g1->SetPoint(0,-4,-3);
   g1->SetPoint(1,1,1);
   g1->SetPoint(2,2,1);
   g1->SetPoint(3,3,4);
   g1->SetPoint(4,5,5);
   g1->SetPointError(0,1.,2.);
   g1->SetPointError(1,2,1);
   g1->SetPointError(2,2,3);
   g1->SetPointError(3,3,2);
   g1->SetPointError(4,4,5);
   g1->SetMarkerStyle(21);

   auto *g2 = new TGraph();
   g2->SetPoint(0,4,8);
   g2->SetPoint(1,5,9);
   g2->SetPoint(2,6,10);
   g2->SetPoint(3,10,11);
   g2->SetPoint(4,15,12);
   g2->SetLineColor(kRed);
   g2->SetLineWidth(5);

   auto mg = new TMultiGraph();
   mg->Add(g1,"P");
   mg->Add(g2,"L");

   c->cd(1); gPad->SetGrid(1,1);
   mg->Draw("A");

   c->cd(2); gPad->SetGrid(1,1);
   mg->Draw("A RX RY");
}
End_Macro

### <a name="MG02"></a> MultiGraphs' fitting

The following example shows how to fit a TMultiGraph.

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","c1",600,400);

   Double_t px1[2] = {2.,4.};
   Double_t dx1[2] = {0.1,0.1};
   Double_t py1[2] = {2.1,4.0};
   Double_t dy1[2] = {0.3,0.2};

   Double_t px2[2] = {3.,5.};
   Double_t dx2[2] = {0.1,0.1};
   Double_t py2[2] = {3.2,4.8};
   Double_t dy2[2] = {0.3,0.2};

   gStyle->SetOptFit(0001);

   auto g1 = new TGraphErrors(2,px1,py1,dx1,dy1);
   g1->SetMarkerStyle(21);
   g1->SetMarkerColor(2);

   auto g2 = new TGraphErrors(2,px2,py2,dx2,dy2);
   g2->SetMarkerStyle(22);
   g2->SetMarkerColor(3);

   auto g = new TMultiGraph();
   g->Add(g1);
   g->Add(g2);

   g->Draw("AP");

   g->Fit("pol1","FQ");
}
End_Macro

#### <a name="MG02a"></a> Fit box position

When the graphs in a TMultiGraph are fitted, the fit parameters boxes
overlap. The following example shows how to make them all visible.


Begin_Macro(source)
../../../tutorials/graphs/multigraph.C
End_Macro

### <a name="MG03"></a> Axis' limits setting

The axis limits can be changed the like for TGraph. The same methods apply on
the multigraph.
Note the two differents ways to change limits on X and Y axis.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",600,400);

   TGraph *g[3];
   Double_t x[10] = {0,1,2,3,4,5,6,7,8,9};
   Double_t y[10] = {1,2,3,4,5,5,4,3,2,1};
   auto mg = new TMultiGraph();
   for (int i=0; i<3; i++) {
      g[i] = new TGraph(10, x, y);
      g[i]->SetMarkerStyle(20);
      g[i]->SetMarkerColor(i+2);
      for (int j=0; j<10; j++) y[j] = y[j]-1;
      mg->Add(g[i]);
   }
   mg->Draw("APL");
   mg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
   mg->GetYaxis()->SetTitle("Coefficients");

   // Change the axis limits
   gPad->Modified();
   mg->GetXaxis()->SetLimits(1.5,7.5);
   mg->SetMinimum(0.);
   mg->SetMaximum(10.);
}
End_Macro
*/


////////////////////////////////////////////////////////////////////////////////
/// TMultiGraph default constructor.

TMultiGraph::TMultiGraph(): TNamed()
{
   fGraphs    = 0;
   fFunctions = 0;
   fHistogram = 0;
   fMaximum   = -1111;
   fMinimum   = -1111;
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor with name and title.

TMultiGraph::TMultiGraph(const char *name, const char *title)
       : TNamed(name,title)
{
   fGraphs    = 0;
   fFunctions = 0;
   fHistogram = 0;
   fMaximum   = -1111;
   fMinimum   = -1111;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TMultiGraph::TMultiGraph(const TMultiGraph& mg) :
  TNamed (mg),
  fGraphs(mg.fGraphs),
  fFunctions(mg.fFunctions),
  fHistogram(mg.fHistogram),
  fMaximum(mg.fMaximum),
  fMinimum(mg.fMinimum)
{
}


////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TMultiGraph& TMultiGraph::operator=(const TMultiGraph& mg)
{
   if (this!=&mg) {
      TNamed::operator=(mg);
      fGraphs=mg.fGraphs;
      fFunctions=mg.fFunctions;
      fHistogram=mg.fHistogram;
      fMaximum=mg.fMaximum;
      fMinimum=mg.fMinimum;
   }
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// TMultiGraph destructor.

TMultiGraph::~TMultiGraph()
{
   if (!fGraphs) return;
   TGraph *g;
   TIter   next(fGraphs);
   while ((g = (TGraph*) next())) {
      g->ResetBit(kMustCleanup);
   }
   fGraphs->Delete();
   delete fGraphs;
   fGraphs = 0;
   delete fHistogram;
   fHistogram = 0;
   if (fFunctions) {
      fFunctions->SetBit(kInvalidObject);
      //special logic to support the case where the same object is
      //added multiple times in fFunctions.
      //This case happens when the same object is added with different
      //drawing modes
      TObject *obj;
      while ((obj  = fFunctions->First())) {
         while (fFunctions->Remove(obj)) { }
         delete obj;
      }
      delete fFunctions;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Add a new graph to the list of graphs.
/// Note that the graph is now owned by the TMultigraph.
/// Deleting the TMultiGraph object will automatically delete the graphs.
/// You should not delete the graphs when the TMultigraph is still active.

void TMultiGraph::Add(TGraph *graph, Option_t *chopt)
{
   if (!fGraphs) fGraphs = new TList();
   graph->SetBit(kMustCleanup);
   fGraphs->Add(graph,chopt);
}


////////////////////////////////////////////////////////////////////////////////
/// Add all the graphs in "multigraph" to the list of graphs.
///
///   - If "chopt" is defined all the graphs in "multigraph" will be added with
///     the "chopt" option.
///   - If "chopt" is undefined each graph will be added with the option it had
///     in "multigraph".

void TMultiGraph::Add(TMultiGraph *multigraph, Option_t *chopt)
{
   TList *graphlist = multigraph->GetListOfGraphs();
   if (!graphlist) return;

   if (!fGraphs) fGraphs = new TList();

   TObjOptLink *lnk = (TObjOptLink*)graphlist->FirstLink();
   TObject *obj = 0;

   while (lnk) {
      obj = lnk->GetObject();
      if (!strlen(chopt)) fGraphs->Add(obj,lnk->GetOption());
      else                fGraphs->Add(obj,chopt);
      lnk = (TObjOptLink*)lnk->Next();
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Browse multigraph.

void TMultiGraph::Browse(TBrowser *b)
{
   TString opt = gEnv->GetValue("TGraph.BrowseOption", "");
   if (opt.IsNull()) {
      opt = b ? b->GetDrawOption() : "alp";
      opt = (opt == "") ? "alp" : opt.Data();
   }
   Draw(opt.Data());
   gPad->Update();
}


////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to each graph.

Int_t TMultiGraph::DistancetoPrimitive(Int_t px, Int_t py)
{
   // Are we on the axis?
   const Int_t kMaxDiff = 10;
   Int_t distance = 9999;
   if (fHistogram) {
      distance = fHistogram->DistancetoPrimitive(px,py);
      if (distance <= 0) return distance;
   }

   // Loop on the list of graphs
   if (!fGraphs) return distance;
   TGraph *g;
   TIter   next(fGraphs);
   while ((g = (TGraph*) next())) {
      Int_t dist = g->DistancetoPrimitive(px,py);
      if (dist <= 0) return 0;
      if (dist < kMaxDiff) {gPad->SetSelected(g); return dist;}
   }
   return distance;
}


////////////////////////////////////////////////////////////////////////////////
/// Draw this multigraph with its current attributes.
///
///  Options to draw a graph are described in TGraphPainter.
///
///  The drawing option for each TGraph may be specified as an optional
///  second argument of the Add function. You can use GetGraphDrawOption
///  to return this option.
///
///  If a draw option is specified, it will be used to draw the graph,
///  otherwise the graph will be drawn with the option specified in
///  TMultiGraph::Draw. Use GetDrawOption to return the option specified
///  when drawing the TMultiGraph.

void TMultiGraph::Draw(Option_t *option)
{
   TString opt = option;
   opt.ToLower();

   if (gPad) {
      if (!gPad->IsEditable()) gROOT->MakeDefCanvas();
      if (opt.Contains("a")) gPad->Clear();
   }
   AppendPad(option);
}


////////////////////////////////////////////////////////////////////////////////
/// Fit this graph with function with name fname.
///
///  interface to TF1::Fit(TF1 *f1...

TFitResultPtr TMultiGraph::Fit(const char *fname, Option_t *option, Option_t *, Axis_t xmin, Axis_t xmax)
{
   char *linear;
   linear= (char*)strstr(fname, "++");
   if (linear) {
      TF1 f1(fname, fname, xmin, xmax);
      return Fit(&f1,option,"",xmin,xmax);
   }
   TF1 * f1 = (TF1*)gROOT->GetFunction(fname);
   if (!f1) { Printf("Unknown function: %s",fname); return -1; }

   return Fit(f1,option,"",xmin,xmax);
}


////////////////////////////////////////////////////////////////////////////////
/// Fit this multigraph with function f1.
///
///  In this function all graphs of the multigraph are fitted simultaneously
///
///  f1 is an already predefined function created by TF1.
///  Predefined functions such as gaus, expo and poln are automatically
///  created by ROOT.
///
///  The list of fit options is given in parameter `option`which may takes the
///  following values:
///
///   - "W" Ignore all the point errors
///   - "U" Use a User specified fitting algorithm (via SetFCN)
///   - "Q" Quiet mode (minimum printing)
///   - "V" Verbose mode (default is between Q and V)
///   - "B" Use this option when you want to fix one or more parameters
///                   and the fitting function is like "gaus","expo","poln","landau".
///   - "R" Use the Range specified in the function range
///   - "N" Do not store the graphics function, do not draw
///   - "0" Do not plot the result of the fit. By default the fitted function
///     is drawn unless the option"N" above is specified.
///   - "+" Add this new fitted function to the list of fitted functions
///     (by default, any previous function is deleted)
///   - "C" In case of linear fitting, not calculate the chisquare (saves time)
///   - "F" If fitting a polN, switch to minuit fitter
///   - "ROB" In case of linear fitting, compute the LTS regression
///      coefficients (robust(resistant) regression), using
///      the default fraction of good points
///   - "ROB=0.x" - compute the LTS regression coefficients, using
///     0.x as a fraction of good points
///
///  When the fit is drawn (by default), the parameter goption may be used
///  to specify a list of graphics options. See TGraph::Paint for a complete
///  list of these options.
///
///  In order to use the Range option, one must first create a function
///  with the expression to be fitted. For example, if your graph
///  has a defined range between -4 and 4 and you want to fit a gaussian
///  only in the interval 1 to 3, you can do:
/// ~~~ {.cpp}
///        TF1 *f1 = new TF1("f1","gaus",1,3);
///        graph->Fit("f1","R");
/// ~~~
///
///  ### Who is calling this function ?
///
///  Note that this function is called when calling TGraphErrors::Fit
///  or TGraphAsymmErrors::Fit ot TGraphBentErrors::Fit
///  see the discussion below on the errors calculation.
///
///  ### Setting initial conditions
///
///  Parameters must be initialized before invoking the Fit function.
///  The setting of the parameter initial values is automatic for the
///  predefined functions : poln, expo, gaus, landau. One can however disable
///  this automatic computation by specifying the option "B".
///  You can specify boundary limits for some or all parameters via
/// ~~~ {.cpp}
///        f1->SetParLimits(p_number, parmin, parmax);
/// ~~~
///  if `parmin>=parmax`, the parameter is fixed
///  Note that you are not forced to fix the limits for all parameters.
///  For example, if you fit a function with 6 parameters, you can do:
/// ~~~ {.cpp}
///     func->SetParameters(0,3.1,1.e-6,0.1,-8,100);
///     func->SetParLimits(4,-10,-4);
///     func->SetParLimits(5, 1,1);
/// ~~~
///  With this setup, parameters 0->3 can vary freely
///  Parameter 4 has boundaries [-10,-4] with initial value -8
///  Parameter 5 is fixed to 100.
///
///  ### Fit range
///
///  The fit range can be specified in two ways:
///
///   - specify rxmax > rxmin (default is rxmin=rxmax=0)
///   - specify the option "R". In this case, the function will be taken
///     instead of the full graph range.
///
///  ### Changing the fitting function
///
///   By default a chi2 fitting function is used for fitting the TGraphs's.
///   The function is implemented in `FitUtil::EvaluateChi2`.
///   In case of TGraphErrors an effective chi2 is used
///   (see TGraphErrors fit in TGraph::Fit) and is implemented in
///   `FitUtil::EvaluateChi2Effective`
///   To specify a User defined fitting function, specify option "U" and
///   call the following function:
/// ~~~ {.cpp}
///   TVirtualFitter::Fitter(mygraph)->SetFCN(MyFittingFunction)
/// ~~~
///   where MyFittingFunction is of type:
/// ~~~ {.cpp}
///   extern void MyFittingFunction(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag);
/// ~~~
///
///  ### Access to the fit result
///
///  The function returns a TFitResultPtr which can hold a  pointer to a TFitResult object.
///  By default the TFitResultPtr contains only the status of the fit and it converts
///  automatically to an integer. If the option "S" is instead used, TFitResultPtr contains
///  the TFitResult and behaves as a smart pointer to it. For example one can do:
/// ~~~ {.cpp}
///     TFitResultPtr r = graph->Fit("myFunc","S");
///     TMatrixDSym cov = r->GetCovarianceMatrix();  //  to access the covariance matrix
///     Double_t par0   = r->Parameter(0); // retrieve the value for the parameter 0
///     Double_t err0   = r->ParError(0); // retrieve the error for the parameter 0
///     r->Print("V");     // print full information of fit including covariance matrix
///     r->Write();        // store the result in a file
/// ~~~
///
///   The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also
///   from the fitted function.
///
///  ### Associated functions
///
///  One or more object (typically a TF1*) can be added to the list
///  of functions (fFunctions) associated to each graph.
///  When TGraph::Fit is invoked, the fitted function is added to this list.
///  Given a graph gr, one can retrieve an associated function
///  with:
/// ~~~ {.cpp}
///   TF1 *myfunc = gr->GetFunction("myfunc");
/// ~~~
///
///  If the graph is made persistent, the list of
///  associated functions is also persistent. Given a pointer (see above)
///  to an associated function myfunc, one can retrieve the function/fit
///  parameters with calls such as:
/// ~~~ {.cpp}
///    Double_t chi2 = myfunc->GetChisquare();
///    Double_t par0 = myfunc->GetParameter(0); //value of 1st parameter
///    Double_t err0 = myfunc->GetParError(0);  //error on first parameter
/// ~~~
///
///  ### Fit Statistics
///
///  You can change the statistics box to display the fit parameters with
///  the TStyle::SetOptFit(mode) method. This mode has four digits.
///  mode = pcev  (default = 0111)
///
///   - v = 1;  print name/values of parameters
///   - e = 1;  print errors (if e=1, v must be 1)
///   - c = 1;  print Chisquare/Number of degrees of freedom
///   - p = 1;  print Probability
///
///  For example: `gStyle->SetOptFit(1011);`
///  prints the fit probability, parameter names/values, and errors.
///  You can change the position of the statistics box with these lines
///  (where g is a pointer to the TGraph):
///
/// ~~~ {.cpp}
///  Root > TPaveStats *st = (TPaveStats*)g->GetListOfFunctions()->FindObject("stats")
///  Root > st->SetX1NDC(newx1); //new x start position
///  Root > st->SetX2NDC(newx2); //new x end position
/// ~~~

TFitResultPtr TMultiGraph::Fit(TF1 *f1, Option_t *option, Option_t *goption, Axis_t rxmin, Axis_t rxmax)
{
   // internal multigraph fitting methods
   Foption_t fitOption;
   ROOT::Fit::FitOptionsMake(ROOT::Fit::kGraph,option,fitOption);

   // create range and minimizer options with default values
   ROOT::Fit::DataRange range(rxmin,rxmax);
   ROOT::Math::MinimizerOptions minOption;
   return ROOT::Fit::FitObject(this, f1 , fitOption , minOption, goption, range);

}

////////////////////////////////////////////////////////////////////////////////
/// Display a panel with all histogram fit options.
/// See class TFitPanel for example

void TMultiGraph::FitPanel()
{
   if (!gPad)
      gROOT->MakeDefCanvas();

   if (!gPad) {
      Error("FitPanel", "Unable to create a default canvas");
      return;
   }

   // use plugin manager to create instance of TFitEditor
   TPluginHandler *handler = gROOT->GetPluginManager()->FindHandler("TFitEditor");
   if (handler && handler->LoadPlugin() != -1) {
      if (handler->ExecPlugin(2, gPad, this) == 0)
         Error("FitPanel", "Unable to crate the FitPanel");
   }
   else
         Error("FitPanel", "Unable to find the FitPanel plug-in");
}

////////////////////////////////////////////////////////////////////////////////
/// Return the draw option for the TGraph `gr` in this TMultiGraph.
/// The return option is the one specified when calling TMultiGraph::Add(gr,option).

Option_t *TMultiGraph::GetGraphDrawOption(const TGraph *gr) const
{
   if (!fGraphs || !gr) return "";
   TListIter next(fGraphs);
   TObject *obj;
   while ((obj = next())) {
      if (obj == (TObject*)gr) return next.GetOption();
   }
   return "";
}


////////////////////////////////////////////////////////////////////////////////
/// Compute Initial values of parameters for a gaussian.

void TMultiGraph::InitGaus(Double_t xmin, Double_t xmax)
{
   Double_t allcha, sumx, sumx2, x, val, rms, mean;
   Int_t bin;
   const Double_t sqrtpi = 2.506628;

   // Compute mean value and RMS of the graph in the given range
   Int_t np = 0;
   allcha = sumx = sumx2 = 0;
   TGraph *g;
   TIter next(fGraphs);
   Double_t *px, *py;
   Int_t npp; //number of points in each graph
   while ((g = (TGraph*) next())) {
      px=g->GetX();
      py=g->GetY();
      npp=g->GetN();
      for (bin=0; bin<npp; bin++) {
         x=px[bin];
         if (x<xmin || x>xmax) continue;
         np++;
         val=py[bin];
         sumx+=val*x;
         sumx2+=val*x*x;
         allcha+=val;
      }
   }
   if (np == 0 || allcha == 0) return;
   mean = sumx/allcha;
   rms  = TMath::Sqrt(sumx2/allcha - mean*mean);

   Double_t binwidx = TMath::Abs((xmax-xmin)/np);
   if (rms == 0) rms = 1;
   TVirtualFitter *grFitter = TVirtualFitter::GetFitter();
   TF1 *f1 = (TF1*)grFitter->GetUserFunc();
   f1->SetParameter(0,binwidx*allcha/(sqrtpi*rms));
   f1->SetParameter(1,mean);
   f1->SetParameter(2,rms);
   f1->SetParLimits(2,0,10*rms);
}


////////////////////////////////////////////////////////////////////////////////
/// Compute Initial values of parameters for an exponential.

void TMultiGraph::InitExpo(Double_t xmin, Double_t xmax)
{
   Double_t constant, slope;
   Int_t ifail;

   LeastSquareLinearFit(-1, constant, slope, ifail, xmin, xmax);

   TVirtualFitter *grFitter = TVirtualFitter::GetFitter();
   TF1 *f1 = (TF1*)grFitter->GetUserFunc();
   f1->SetParameter(0,constant);
   f1->SetParameter(1,slope);
}


////////////////////////////////////////////////////////////////////////////////
/// Compute Initial values of parameters for a polynom.

void TMultiGraph::InitPolynom(Double_t xmin, Double_t xmax)
{
   Double_t fitpar[25];

   TVirtualFitter *grFitter = TVirtualFitter::GetFitter();
   TF1 *f1 = (TF1*)grFitter->GetUserFunc();
   Int_t npar   = f1->GetNpar();

   LeastSquareFit(npar, fitpar, xmin, xmax);

   for (Int_t i=0;i<npar;i++) f1->SetParameter(i, fitpar[i]);
}


////////////////////////////////////////////////////////////////////////////////
/// Least squares lpolynomial fitting without weights.
///
///   - m     number of parameters
///   - a     array of parameters
///   - first 1st point number to fit (default =0)
///   - last  last point number to fit (default=fNpoints-1)
///
///   based on CERNLIB routine LSQ: Translated to C++ by Rene Brun

void TMultiGraph::LeastSquareFit(Int_t m, Double_t *a, Double_t xmin, Double_t xmax)
{
   const Double_t zero = 0.;
   const Double_t one = 1.;
   const Int_t idim = 20;

   Double_t  b[400]        /* was [20][20] */;
   Int_t i, k, l, ifail, bin;
   Double_t power;
   Double_t da[20], xk, yk;


   //count the total number of points to fit
   TGraph *g;
   TIter next(fGraphs);
   Double_t *px, *py;
   Int_t n=0;
   Int_t npp;
   while ((g = (TGraph*) next())) {
      px=g->GetX();
      npp=g->GetN();
      for (bin=0; bin<npp; bin++) {
         xk=px[bin];
         if (xk < xmin || xk > xmax) continue;
         n++;
      }
   }
   if (m <= 2) {
      LeastSquareLinearFit(n, a[0], a[1], ifail, xmin, xmax);
      return;
   }
   if (m > idim || m > n) return;
   da[0] = zero;
   for (l = 2; l <= m; ++l) {
      b[l-1]           = zero;
      b[m + l*20 - 21] = zero;
      da[l-1]          = zero;
   }
   Int_t np = 0;

   next.Reset();
   while ((g = (TGraph*) next())) {
      px=g->GetX();
      py=g->GetY();
      npp=g->GetN();

      for (k = 0; k <= npp; ++k) {
         xk     = px[k];
         if (xk < xmin || xk > xmax) continue;
         np++;
         yk     = py[k];
         power  = one;
         da[0] += yk;
         for (l = 2; l <= m; ++l) {
            power   *= xk;
            b[l-1]  += power;
            da[l-1] += power*yk;
         }
         for (l = 2; l <= m; ++l) {
            power            *= xk;
            b[m + l*20 - 21] += power;
         }
      }
   }
   b[0]  = Double_t(np);
   for (i = 3; i <= m; ++i) {
      for (k = i; k <= m; ++k) {
         b[k - 1 + (i-1)*20 - 21] = b[k + (i-2)*20 - 21];
      }
   }
   H1LeastSquareSeqnd(m, b, idim, ifail, 1, da);

   if (ifail < 0) {
      //a[0] = fY[0];
      py=((TGraph *)fGraphs->First())->GetY();
      a[0]=py[0];
      for (i=1; i<m; ++i) a[i] = 0;
      return;
   }
   for (i=0; i<m; ++i) a[i] = da[i];
}


////////////////////////////////////////////////////////////////////////////////
/// Least square linear fit without weights.
///
///  Fit a straight line (a0 + a1*x) to the data in this graph.
///
///   - ndata:  number of points to fit
///   - first:  first point number to fit
///   - last:   last point to fit O(ndata should be last-first
///   - ifail:  return parameter indicating the status of the fit (ifail=0, fit is OK)
///
///   extracted from CERNLIB LLSQ: Translated to C++ by Rene Brun

void TMultiGraph::LeastSquareLinearFit(Int_t ndata, Double_t &a0, Double_t &a1,
                                       Int_t &ifail, Double_t xmin, Double_t xmax)
{
   Double_t xbar, ybar, x2bar;
   Int_t i;
   Double_t xybar;
   Double_t fn, xk, yk;
   Double_t det;

   ifail = -2;
   xbar  = ybar = x2bar = xybar = 0;
   Int_t np = 0;
   TGraph *g;
   TIter next(fGraphs);
   Double_t *px, *py;
   Int_t npp;
   while ((g = (TGraph*) next())) {
      px=g->GetX();
      py=g->GetY();
      npp=g->GetN();
      for (i = 0; i < npp; ++i) {
         xk = px[i];
         if (xk < xmin || xk > xmax) continue;
         np++;
         yk = py[i];
         if (ndata < 0) {
            if (yk <= 0) yk = 1e-9;
            yk = TMath::Log(yk);
         }
         xbar  += xk;
         ybar  += yk;
         x2bar += xk*xk;
         xybar += xk*yk;
      }
   }
   fn    = Double_t(np);
   det   = fn*x2bar - xbar*xbar;
   ifail = -1;
   if (det <= 0) {
      if (fn > 0) a0 = ybar/fn;
      else        a0 = 0;
      a1 = 0;
      return;
   }
   ifail = 0;
   a0 = (x2bar*ybar - xbar*xybar) / det;
   a1 = (fn*xybar - xbar*ybar) / det;
}


////////////////////////////////////////////////////////////////////////////////
/// Return 1 if the point (x,y) is inside one of the graphs 0 otherwise.

Int_t TMultiGraph::IsInside(Double_t x, Double_t y) const
{
   Int_t in = 0;
   if (!fGraphs) return in;
   TGraph *g;
   TIter next(fGraphs);
   while ((g = (TGraph*) next())) {
      in = g->IsInside(x, y);
      if (in) return in;
   }
   return in;
}


////////////////////////////////////////////////////////////////////////////////
/// Returns a pointer to the histogram used to draw the axis.
/// Takes into account following cases.
///
///    1. if `fHistogram` exists it is returned
///    2. if `fHistogram` doesn't exists and `gPad` exists `gPad` is updated. That
///       may trigger the creation of `fHistogram`. If `fHistogram` still does not
///       exit but `hframe` does (if user called `TPad::DrawFrame`) the pointer to
///       `hframe` histogram is returned
///    3. after the two previous steps, if `fHistogram` still doesn't exist, then
///       it is created.

TH1F *TMultiGraph::GetHistogram()
{
   if (fHistogram) return fHistogram;

   if (gPad) {
      gPad->Modified();
      gPad->Update();
      if (fHistogram) return fHistogram;
      TH1F *h1 = (TH1F*)gPad->FindObject("hframe");
      if (h1) return h1;
   }

   Bool_t initialrangeset = kFALSE;
   Double_t rwxmin = 0.,rwxmax = 0.,rwymin = 0.,rwymax = 0.;
   TGraph *g;
   Int_t npt = 100 ;
   TIter   next(fGraphs);
   while ((g = (TGraph*) next())) {
      if (g->GetN() <= 0) continue;
      if (initialrangeset) {
         Double_t rx1,ry1,rx2,ry2;
         g->ComputeRange(rx1, ry1, rx2, ry2);
         if (rx1 < rwxmin) rwxmin = rx1;
         if (ry1 < rwymin) rwymin = ry1;
         if (rx2 > rwxmax) rwxmax = rx2;
         if (ry2 > rwymax) rwymax = ry2;
      } else {
         g->ComputeRange(rwxmin, rwymin, rwxmax, rwymax);
         initialrangeset = kTRUE;
      }
      if (g->GetN() > npt) npt = g->GetN();
   }
   if (rwxmin == rwxmax) rwxmax += 1.;
   if (rwymin == rwymax) rwymax += 1.;
   double dx = 0.05*(rwxmax-rwxmin);
   double dy = 0.05*(rwymax-rwymin);
   rwxmin = rwxmin - dx;
   rwxmax = rwxmax + dx;
   rwymin = rwymin - dy;
   rwymax = rwymax + dy;
   fHistogram = new TH1F(GetName(),GetTitle(),npt,rwxmin,rwxmax);
   if (!fHistogram) return 0;
   fHistogram->SetMinimum(rwymin);
   fHistogram->SetBit(TH1::kNoStats);
   fHistogram->SetMaximum(rwymax);
   fHistogram->GetYaxis()->SetLimits(rwymin,rwymax);
   fHistogram->SetDirectory(0);
   return fHistogram;
}


////////////////////////////////////////////////////////////////////////////////
/// Return pointer to function with name.
///
/// Functions such as TGraph::Fit store the fitted function in the list of
/// functions of this graph.

TF1 *TMultiGraph::GetFunction(const char *name) const
{
   if (!fFunctions) return 0;
   return (TF1*)fFunctions->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to list of functions.
/// If pointer is null create the list

TList *TMultiGraph::GetListOfFunctions()
{
   if (!fFunctions) fFunctions = new TList();
   return fFunctions;
}


////////////////////////////////////////////////////////////////////////////////
/// Get x axis of the graph.
/// This method returns a valid axis only after the TMultigraph has been drawn.

TAxis *TMultiGraph::GetXaxis()
{
   TH1 *h = GetHistogram();
   if (!h) return 0;
   return h->GetXaxis();
}


////////////////////////////////////////////////////////////////////////////////
/// Get y axis of the graph.
/// This method returns a valid axis only after the TMultigraph has been drawn.

TAxis *TMultiGraph::GetYaxis()
{
   TH1 *h = GetHistogram();
   if (!h) return 0;
   return h->GetYaxis();
}


////////////////////////////////////////////////////////////////////////////////
/// Paint all the graphs of this multigraph.

void TMultiGraph::Paint(Option_t *choptin)
{
   const TPickerStackGuard pushGuard(this);

   if (!fGraphs) return;
   if (fGraphs->GetSize() == 0) return;

   char option[128];
   strlcpy(option,choptin,128);
   Int_t nch = strlen(choptin);
   for (Int_t i=0;i<nch;i++) option[i] = toupper(option[i]);

   // Automatic color
   char *l1 = strstr(option,"PFC"); // Automatic Fill Color
   char *l2 = strstr(option,"PLC"); // Automatic Line Color
   char *l3 = strstr(option,"PMC"); // Automatic Marker Color
   if (l1 || l2 || l3) {
      TString opt1 = option; opt1.ToLower();
      if (l1) memcpy(l1,"   ",3);
      if (l2) memcpy(l2,"   ",3);
      if (l3) memcpy(l3,"   ",3);
      TObjOptLink *lnk = (TObjOptLink*)fGraphs->FirstLink();
      TGraph* gAti;
      Int_t ngraphs = fGraphs->GetSize();
      Int_t ic;
      gPad->IncrementPaletteColor(ngraphs, opt1);
      for (Int_t i=0;i<ngraphs;i++) {
         ic = gPad->NextPaletteColor();
         gAti = (TGraph*)(fGraphs->At(i));
         if (l1) gAti->SetFillColor(ic);
         if (l2) gAti->SetLineColor(ic);
         if (l3) gAti->SetMarkerColor(ic);
         lnk = (TObjOptLink*)lnk->Next();
      }
   }

   char *l;

   TString chopt = option;

   l = (char*)strstr(chopt.Data(),"3D");
   if (l) {
      l = (char*)strstr(chopt.Data(),"L");
      if (l) PaintPolyLine3D(chopt.Data());
      return;
   }

   l = (char*)strstr(chopt.Data(),"PADS");
   if (l) {
      chopt.ReplaceAll("PADS","");
      PaintPads(chopt.Data());
      return;
   }

   char *lrx = (char *)strstr(chopt.Data(), "RX"); // Reverse graphs along X axis
   char *lry = (char *)strstr(chopt.Data(), "RY"); // Reverse graphs along Y axis
   if (lrx || lry) {
      PaintReverse(chopt.Data());
      return;
   }

   TGraph *g;

   l = (char*)strstr(chopt.Data(),"A");
   if (l) {
      *l = ' ';
      TIter   next(fGraphs);
      Int_t npt = 100;
      Double_t maximum, minimum, rwxmin, rwxmax, rwymin, rwymax, uxmin, uxmax, dx, dy;
      rwxmin    = gPad->GetUxmin();
      rwxmax    = gPad->GetUxmax();
      rwymin    = gPad->GetUymin();
      rwymax    = gPad->GetUymax();
      char *xtitle = 0;
      char *ytitle = 0;
      Int_t firstx = 0;
      Int_t lastx  = 0;
      Bool_t timedisplay = kFALSE;
      char *timeformat = 0;

      if (fHistogram) {
         //cleanup in case of a previous unzoom and in case one of the TGraph has changed
         TObjOptLink *lnk = (TObjOptLink*)fGraphs->FirstLink();
         TGraph* gAti;
         Int_t ngraphs = fGraphs->GetSize();
         Bool_t reset_hist = kFALSE;
         for (Int_t i=0;i<ngraphs;i++) {
            gAti = (TGraph*)(fGraphs->At(i));
            if(gAti->TestBit(TGraph::kResetHisto)) {reset_hist = kTRUE; break;}
            lnk = (TObjOptLink*)lnk->Next();
         }
         if (fHistogram->GetMinimum() >= fHistogram->GetMaximum() || reset_hist) {
            nch = strlen(fHistogram->GetXaxis()->GetTitle());
            firstx = fHistogram->GetXaxis()->GetFirst();
            lastx  = fHistogram->GetXaxis()->GetLast();
            timedisplay = fHistogram->GetXaxis()->GetTimeDisplay();
            if (nch) {
               xtitle = new char[nch+1];
               strlcpy(xtitle,fHistogram->GetXaxis()->GetTitle(),nch+1);
            }
            nch = strlen(fHistogram->GetYaxis()->GetTitle());
            if (nch) {
               ytitle = new char[nch+1];
               strlcpy(ytitle,fHistogram->GetYaxis()->GetTitle(),nch+1);
            }
            nch = strlen(fHistogram->GetXaxis()->GetTimeFormat());
            if (nch) {
              timeformat = new char[nch+1];
              strlcpy(timeformat,fHistogram->GetXaxis()->GetTimeFormat(),nch+1);
            }
            delete fHistogram;
            fHistogram = 0;
         }
      }
      if (fHistogram) {
         minimum = fHistogram->GetYaxis()->GetXmin();
         maximum = fHistogram->GetYaxis()->GetXmax();
         uxmin   = gPad->PadtoX(rwxmin);
         uxmax   = gPad->PadtoX(rwxmax);
      } else {
         Bool_t initialrangeset = kFALSE;
         while ((g = (TGraph*) next())) {
            if (g->GetN() <= 0) continue;
            if (initialrangeset) {
               Double_t rx1,ry1,rx2,ry2;
               g->ComputeRange(rx1, ry1, rx2, ry2);
               if (rx1 < rwxmin) rwxmin = rx1;
               if (ry1 < rwymin) rwymin = ry1;
               if (rx2 > rwxmax) rwxmax = rx2;
               if (ry2 > rwymax) rwymax = ry2;
            } else {
               g->ComputeRange(rwxmin, rwymin, rwxmax, rwymax);
               initialrangeset = kTRUE;
            }
            if (g->GetN() > npt) npt = g->GetN();
         }
         if (rwxmin == rwxmax) rwxmax += 1.;
         if (rwymin == rwymax) rwymax += 1.;
         dx = 0.05*(rwxmax-rwxmin);
         dy = 0.05*(rwymax-rwymin);
         uxmin    = rwxmin - dx;
         uxmax    = rwxmax + dx;
         if (gPad->GetLogy()) {
            if (rwymin <= 0) rwymin = 0.001*rwymax;
            minimum = rwymin/(1+0.5*TMath::Log10(rwymax/rwymin));
            maximum = rwymax*(1+0.2*TMath::Log10(rwymax/rwymin));
         } else {
            minimum  = rwymin - dy;
            maximum  = rwymax + dy;
         }
         if (minimum < 0 && rwymin >= 0) minimum = 0;
         if (maximum > 0 && rwymax <= 0) maximum = 0;
      }

      if (fMinimum != -1111) rwymin = minimum = fMinimum;
      if (fMaximum != -1111) rwymax = maximum = fMaximum;
      if (uxmin < 0 && rwxmin >= 0) {
         if (gPad->GetLogx()) uxmin = 0.9*rwxmin;
         //else                 uxmin = 0;
      }
      if (uxmax > 0 && rwxmax <= 0) {
         if (gPad->GetLogx()) uxmax = 1.1*rwxmax;
         //else                 uxmax = 0;
      }
      if (minimum < 0 && rwymin >= 0) {
         if (gPad->GetLogy()) minimum = 0.9*rwymin;
         //else                minimum = 0;
      }
      if (maximum > 0 && rwymax <= 0) {
         if (gPad->GetLogy()) maximum = 1.1*rwymax;
         //else                maximum = 0;
      }
      if (minimum <= 0 && gPad->GetLogy()) minimum = 0.001*maximum;
      if (uxmin <= 0 && gPad->GetLogx()) {
         if (uxmax > 1000) uxmin = 1;
         else              uxmin = 0.001*uxmax;
      }
      rwymin = minimum;
      rwymax = maximum;
      if (fHistogram) {
         fHistogram->GetYaxis()->SetLimits(rwymin,rwymax);
      }

      // Create a temporary histogram to draw the axis
      if (!fHistogram) {
         // the graph is created with at least as many channels as there are points
         // to permit zooming on the full range
         rwxmin = uxmin;
         rwxmax = uxmax;
         fHistogram = new TH1F(GetName(),GetTitle(),npt,rwxmin,rwxmax);
         if (!fHistogram) return;
         fHistogram->SetMinimum(rwymin);
         fHistogram->SetBit(TH1::kNoStats);
         fHistogram->SetMaximum(rwymax);
         fHistogram->GetYaxis()->SetLimits(rwymin,rwymax);
         fHistogram->SetDirectory(0);
         if (xtitle) {fHistogram->GetXaxis()->SetTitle(xtitle); delete [] xtitle;}
         if (ytitle) {fHistogram->GetYaxis()->SetTitle(ytitle); delete [] ytitle;}
         if (firstx != lastx) fHistogram->GetXaxis()->SetRange(firstx,lastx);
         if (timedisplay) {fHistogram->GetXaxis()->SetTimeDisplay(timedisplay);}
         if (timeformat) {fHistogram->GetXaxis()->SetTimeFormat(timeformat); delete [] timeformat;}
      }
      fHistogram->Paint("0");
   }

   TGraph *gfit = nullptr;
   if (fGraphs) {
      TObjOptLink *lnk = (TObjOptLink*)fGraphs->FirstLink();
      TObject *obj = 0;

      chopt.ReplaceAll("A","");

      while (lnk) {

         obj = lnk->GetObject();

         gPad->PushSelectableObject(obj);

         if (!gPad->PadInHighlightMode() || (gPad->PadInHighlightMode() && obj == gPad->GetSelected())) {
            TString opt = lnk->GetOption();
            if (!opt.IsWhitespace())
               obj->Paint(opt.ReplaceAll("A","").Data());
            else {
               if (!chopt.IsWhitespace()) obj->Paint(chopt.Data());
               else                       obj->Paint("L");
            }
         }

         lnk = (TObjOptLink*)lnk->Next();
      }

      gfit = (TGraph*)obj; // pick one TGraph in the list to paint the fit parameters.
   }

   TObject *f;
   TF1 *fit = nullptr;
   if (fFunctions) {
      TIter   next(fFunctions);
      while ((f = (TObject*) next())) {
         if (f->InheritsFrom(TF1::Class())) {
            if (f->TestBit(TF1::kNotDraw) == 0) f->Paint("lsame");
            fit = (TF1*)f;
         } else  {
            f->Paint();
         }
      }
   }

   if (gfit && fit) gfit->PaintStats(fit);
}


////////////////////////////////////////////////////////////////////////////////
/// Divides the active pad and draws all Graphs in the Multigraph separately.

void TMultiGraph::PaintPads(Option_t *option)
{
   TIter next(fGraphs);
   Int_t neededPads = fGraphs->GetSize();
   Int_t existingPads = 0;
   TString opt = (TString)option;

   TVirtualPad *curPad = gPad;
   TObject *obj;
   TIter nextPad(curPad->GetListOfPrimitives());

   while ((obj = nextPad())) {
      if (obj->InheritsFrom(TVirtualPad::Class())) existingPads++;
   }
   if (existingPads < neededPads) {
      curPad->Clear();
      Int_t nx = (Int_t)TMath::Sqrt((Double_t)neededPads);
      if (nx*nx < neededPads) nx++;
      Int_t ny = nx;
      if (((nx*ny)-nx) >= neededPads) ny--;
      curPad->Divide(nx,ny);
   }
   Int_t i = 0;
   TGraph *g;

   TObjOptLink *lnk = (TObjOptLink*)fGraphs->FirstLink();
   obj = 0;

   while (lnk) {
      g = (TGraph*)lnk->GetObject();
      i++;
      curPad->cd(i);
      TString apopt = lnk->GetOption();
      if (strlen(apopt)) {
         g->Draw((apopt.Append("A")).Data());
      } else {
         if (strlen(opt)) g->Draw(opt.Append("A"));
         else             g->Draw("LA");
      }
      lnk = (TObjOptLink*)lnk->Next();
   }

   curPad->cd();
}


////////////////////////////////////////////////////////////////////////////////
/// Paint all the graphs of this multigraph as 3D lines.

void TMultiGraph::PaintPolyLine3D(Option_t *option)
{
   Int_t i, npt=0;
   char *l;
   Double_t rwxmin=0., rwxmax=0., rwymin=0., rwymax=0.;
   TIter next(fGraphs);
   TGraph *g;

   g = (TGraph*) next();
   if (g) {
      g->ComputeRange(rwxmin, rwymin, rwxmax, rwymax);
      npt = g->GetN();
   }

   if (!fHistogram) {
      fHistogram = new TH1F(GetName(),GetTitle(),npt,rwxmin,rwxmax);
   }

   while ((g = (TGraph*) next())) {
      Double_t rx1,ry1,rx2,ry2;
      g->ComputeRange(rx1, ry1, rx2, ry2);
      if (rx1 < rwxmin) rwxmin = rx1;
      if (ry1 < rwymin) rwymin = ry1;
      if (rx2 > rwxmax) rwxmax = rx2;
      if (ry2 > rwymax) rwymax = ry2;
      if (g->GetN() > npt) npt = g->GetN();
   }

   Int_t ndiv = fGraphs->GetSize();

   TH2F* frame = new TH2F("frame","", ndiv, 0., (Double_t)(ndiv), npt, rwxmin, rwxmax);
   if (fHistogram) {
      frame->SetTitle(fHistogram->GetTitle());
      frame->GetYaxis()->SetTitle(fHistogram->GetXaxis()->GetTitle());
      frame->GetYaxis()->SetRange(fHistogram->GetXaxis()->GetFirst(), fHistogram->GetXaxis()->GetLast());
      frame->GetZaxis()->SetTitle(fHistogram->GetYaxis()->GetTitle());
   }

   TAxis *Xaxis = frame->GetXaxis();
   Xaxis->SetNdivisions(-ndiv);
   next.Reset();
   for (i=ndiv; i>=1; i--) {
      g = (TGraph*) next();
      Xaxis->SetBinLabel(i, g->GetTitle());
   }

   frame->SetStats(kFALSE);
   if (fMinimum != -1111) frame->SetMinimum(fMinimum);
   else                   frame->SetMinimum(rwymin);
   if (fMaximum != -1111) frame->SetMaximum(fMaximum);
   else                   frame->SetMaximum(rwymax);

   l = (char*)strstr(option,"A");
   if (l) frame->Paint("lego9,fb,bb");
   l = (char*)strstr(option,"BB");
   if (!l) frame->Paint("lego9,fb,a,same");

   Double_t *x, *y;
   Double_t xyz1[3], xyz2[3];

   Double_t xl = frame->GetYaxis()->GetBinLowEdge(frame->GetYaxis()->GetFirst());
   Double_t xu = frame->GetYaxis()->GetBinUpEdge(frame->GetYaxis()->GetLast());
   Double_t yl = frame->GetMinimum();
   Double_t yu = frame->GetMaximum();
   Double_t xc[2],yc[2];
   next.Reset();
   Int_t j = ndiv;

   while ((g = (TGraph*) next())) {
      npt = g->GetN();
      x   = g->GetX();
      y   = g->GetY();
      gPad->SetLineColor(g->GetLineColor());
      gPad->SetLineWidth(g->GetLineWidth());
      gPad->SetLineStyle(g->GetLineStyle());
      gPad->TAttLine::Modify();
      for (i=0; i<npt-1; i++) {
         xc[0] = x[i];
         xc[1] = x[i+1];
         yc[0] = y[i];
         yc[1] = y[i+1];
         if (gPad->Clip(&xc[0], &yc[0], xl, yl, xu, yu)<2) {
            xyz1[0] = j-0.5;
            xyz1[1] = xc[0];
            xyz1[2] = yc[0];
            xyz2[0] = j-0.5;
            xyz2[1] = xc[1];
            xyz2[2] = yc[1];
            gPad->PaintLine3D(xyz1, xyz2);
         }
      }
      j--;
   }

   l = (char*)strstr(option,"FB");
   if (!l) frame->Paint("lego9,bb,a,same");
   delete frame;
}


////////////////////////////////////////////////////////////////////////////////
/// Paint all the graphs of this multigraph reverting values along X and/or Y axis.
/// New graphs are created.

void TMultiGraph::PaintReverse(Option_t *option)
{
   auto *h = GetHistogram();
   TH1F *hg = nullptr;
   TGraph *fg = nullptr;
   if (!h)
      return;
   TString mgopt = option;
   mgopt.ToLower();

   TIter next(fGraphs);
   TGraph *g;
   Bool_t first = kTRUE;
   TString gopt;
   while ((g = (TGraph *)next())) {
      gopt = GetGraphDrawOption(g);
      gopt.Append(mgopt);
      if (first) {
         fg = g;
         hg = fg->GetHistogram();
         fg->SetHistogram(h);
         fg->Paint(gopt.Data());
         first = kFALSE;
      } else {
         g->Paint(gopt.ReplaceAll("a", "").Data());
      }
   }
   if (fg)
      fg->SetHistogram(hg);
}


////////////////////////////////////////////////////////////////////////////////
/// Print the list of graphs.

void TMultiGraph::Print(Option_t *option) const
{
   TGraph *g;
   if (fGraphs) {
      TIter   next(fGraphs);
      while ((g = (TGraph*) next())) {
         g->Print(option);
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Recursively remove this object from a list. Typically implemented
/// by classes that can contain multiple references to a same object.

void TMultiGraph::RecursiveRemove(TObject *obj)
{
   if (!fGraphs) return;
   TObject *objr = fGraphs->Remove(obj);
   if (!objr) return;
   delete fHistogram; fHistogram = 0;
   if (gPad) gPad->Modified();
}


////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out.

void TMultiGraph::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';
   out<<"   "<<std::endl;
   if (gROOT->ClassSaved(TMultiGraph::Class())) {
      out<<"   ";
   } else {
      out<<"   TMultiGraph *";
   }
   out<<"multigraph = new TMultiGraph();"<<std::endl;
   out<<"   multigraph->SetName("<<quote<<GetName()<<quote<<");"<<std::endl;
   out<<"   multigraph->SetTitle("<<quote<<GetTitle()<<quote<<");"<<std::endl;

   if (fGraphs) {
      TObjOptLink *lnk = (TObjOptLink*)fGraphs->FirstLink();
      TObject *g;

      while (lnk) {
         g = lnk->GetObject();
         g->SavePrimitive(out, Form("multigraph%s",lnk->GetOption()));
         lnk = (TObjOptLink*)lnk->Next();
      }
   }
   const char *l = strstr(option,"th2poly");
   if (l) {
      out<<"   "<<l+7<<"->AddBin(multigraph);"<<std::endl;
   } else {
      out<<"   multigraph->Draw(" <<quote<<option<<quote<<");"<<std::endl;
   }
   TAxis *xaxis = GetXaxis();
   TAxis *yaxis = GetYaxis();

   if (xaxis) {
     out<<"   multigraph->GetXaxis()->SetLimits("<<xaxis->GetXmin()<<", "<<xaxis->GetXmax()<<");"<<std::endl;
     xaxis->SaveAttributes(out, "multigraph","->GetXaxis()");
   }
   if (yaxis) yaxis->SaveAttributes(out, "multigraph","->GetYaxis()");
   if (fMinimum != -1111) out<<"   multigraph->SetMinimum("<<fMinimum<<");"<<std::endl;
   if (fMaximum != -1111) out<<"   multigraph->SetMaximum("<<fMaximum<<");"<<std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Set multigraph maximum.

void TMultiGraph::SetMaximum(Double_t maximum)
{
   fMaximum = maximum;
   if (fHistogram)  fHistogram->SetMaximum(maximum);
}


////////////////////////////////////////////////////////////////////////////////
/// Set multigraph minimum.

void TMultiGraph::SetMinimum(Double_t minimum)
{
   fMinimum = minimum;
   if (fHistogram) fHistogram->SetMinimum(minimum);
}


////////////////////////////////////////////////////////////////////////////////
/// Get iterator over internal graphs list.

TIter TMultiGraph::begin() const
{
  return TIter(fGraphs);
}
