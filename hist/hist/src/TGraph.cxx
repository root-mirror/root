// @(#)root/hist:$Id$
// Author: Rene Brun, Olivier Couet   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#include "TROOT.h"
#include "TBuffer.h"
#include "TEnv.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVectorD.h"
#include "Foption.h"
#include "TRandom.h"
#include "TSpline.h"
#include "TVirtualFitter.h"
#include "TVirtualPad.h"
#include "TVirtualGraphPainter.h"
#include "TBrowser.h"
#include "TSystem.h"
#include "TPluginManager.h"
#include "strtok.h"

#include <cstdlib>
#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include <cstring>

#include "HFitInterface.h"
#include "Fit/DataRange.h"
#include "Math/MinimizerOptions.h"

extern void H1LeastSquareSeqnd(Int_t n, Double_t *a, Int_t idim, Int_t &ifail, Int_t k, Double_t *b);

ClassImp(TGraph);

////////////////////////////////////////////////////////////////////////////////

/** \class TGraph
    \ingroup Hist
A TGraph is an object made of two arrays X and Y with npoints each.
The TGraph painting is performed thanks to the TGraphPainter
class. All details about the various painting options are given in this class.

#### Notes

  - Unlike histogram or tree (or even TGraph2D), TGraph objects
    are not automatically attached to the current TFile, in order to keep the
    management and size of the TGraph as small as possible.
  - The TGraph constructors do not have the TGraph title and name as parameters.
    A TGraph has the default title and name "Graph". To change the default title
    and name `SetTitle` and `SetName` should be called on the TGraph after its creation.
    TGraph was a light weight object to start with, like TPolyline or TPolyMarker.
    That’s why it did not have any title and name parameters in the constructors.

The picture below gives an example:

Begin_Macro(source)
{
   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,500,300);
   Double_t x[100], y[100];
   Int_t n = 20;
   for (Int_t i=0;i<n;i++) {
     x[i] = i*0.1;
     y[i] = 10*sin(x[i]+0.2);
   }
   TGraph* gr = new TGraph(n,x,y);
   gr->Draw("AC*");
}
End_Macro
*/

////////////////////////////////////////////////////////////////////////////////
/// Graph default constructor.

TGraph::TGraph(): TNamed(), TAttLine(), TAttFill(0, 1000), TAttMarker()
{
   fNpoints = -1;  //will be reset to 0 in CtorAllocate
   if (!CtorAllocate()) return;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with only the number of points set
/// the arrays x and y will be set later

TGraph::TGraph(Int_t n)
   : TNamed("Graph", "Graph"), TAttLine(), TAttFill(0, 1000), TAttMarker()
{
   fNpoints = n;
   if (!CtorAllocate()) return;
   FillZero(0, fNpoints);
}

////////////////////////////////////////////////////////////////////////////////
/// Graph normal constructor with ints.

TGraph::TGraph(Int_t n, const Int_t *x, const Int_t *y)
   : TNamed("Graph", "Graph"), TAttLine(), TAttFill(0, 1000), TAttMarker()
{
   if (!x || !y) {
      fNpoints = 0;
   } else {
      fNpoints = n;
   }
   if (!CtorAllocate()) return;
   for (Int_t i = 0; i < n; i++) {
      fX[i] = (Double_t)x[i];
      fY[i] = (Double_t)y[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Graph normal constructor with floats.

TGraph::TGraph(Int_t n, const Float_t *x, const Float_t *y)
   : TNamed("Graph", "Graph"), TAttLine(), TAttFill(0, 1000), TAttMarker()
{
   if (!x || !y) {
      fNpoints = 0;
   } else {
      fNpoints = n;
   }
   if (!CtorAllocate()) return;
   for (Int_t i = 0; i < n; i++) {
      fX[i] = x[i];
      fY[i] = y[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Graph normal constructor with doubles.

TGraph::TGraph(Int_t n, const Double_t *x, const Double_t *y)
   : TNamed("Graph", "Graph"), TAttLine(), TAttFill(0, 1000), TAttMarker()
{
   if (!x || !y) {
      fNpoints = 0;
   } else {
      fNpoints = n;
   }
   if (!CtorAllocate()) return;
   n = fNpoints * sizeof(Double_t);
   memcpy(fX, x, n);
   memcpy(fY, y, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor for this graph

TGraph::TGraph(const TGraph &gr)
   : TNamed(gr), TAttLine(gr), TAttFill(gr), TAttMarker(gr)
{
   fNpoints = gr.fNpoints;
   fMaxSize = gr.fMaxSize;
   if (gr.fFunctions) fFunctions = (TList*)gr.fFunctions->Clone();
   else fFunctions = new TList;
   if (gr.fHistogram) fHistogram = (TH1F*)gr.fHistogram->Clone();
   else fHistogram = 0;
   fMinimum = gr.fMinimum;
   fMaximum = gr.fMaximum;
   if (!fMaxSize) {
      fX = fY = 0;
      return;
   } else {
      fX = new Double_t[fMaxSize];
      fY = new Double_t[fMaxSize];
   }

   Int_t n = gr.GetN() * sizeof(Double_t);
   memcpy(fX, gr.fX, n);
   memcpy(fY, gr.fY, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Equal operator for this graph

TGraph& TGraph::operator=(const TGraph &gr)
{
   if (this != &gr) {
      TNamed::operator=(gr);
      TAttLine::operator=(gr);
      TAttFill::operator=(gr);
      TAttMarker::operator=(gr);

      fNpoints = gr.fNpoints;
      fMaxSize = gr.fMaxSize;

      // delete list of functions and their contents before copying it
      if (fFunctions) {
         // delete previous lists of functions
         if (!fFunctions->IsEmpty()) {
            fFunctions->SetBit(kInvalidObject);
            // use TList::Remove to take into account the case the same object is
            // added multiple times in the list
            TObject *obj;
            while ((obj  = fFunctions->First())) {
               while (fFunctions->Remove(obj)) { }
               delete obj;
            }
         }
         delete fFunctions;
      }

      if (gr.fFunctions) fFunctions = (TList*)gr.fFunctions->Clone();
      else fFunctions = new TList;

      if (fHistogram) delete fHistogram;
      if (gr.fHistogram) fHistogram = new TH1F(*(gr.fHistogram));
      else fHistogram = 0;

      fMinimum = gr.fMinimum;
      fMaximum = gr.fMaximum;
      if (fX) delete [] fX;
      if (fY) delete [] fY;
      if (!fMaxSize) {
         fX = fY = 0;
         return *this;
      } else {
         fX = new Double_t[fMaxSize];
         fY = new Double_t[fMaxSize];
      }

      Int_t n = gr.GetN() * sizeof(Double_t);
      if (n > 0) {
         memcpy(fX, gr.fX, n);
         memcpy(fY, gr.fY, n);
      }
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Graph constructor with two vectors of floats in input
/// A graph is build with the X coordinates taken from vx and Y coord from vy
/// The number of points in the graph is the minimum of number of points
/// in vx and vy.

TGraph::TGraph(const TVectorF &vx, const TVectorF &vy)
   : TNamed("Graph", "Graph"), TAttLine(), TAttFill(0, 1000), TAttMarker()
{
   fNpoints = TMath::Min(vx.GetNrows(), vy.GetNrows());
   if (!CtorAllocate()) return;
   Int_t ivxlow  = vx.GetLwb();
   Int_t ivylow  = vy.GetLwb();
   for (Int_t i = 0; i < fNpoints; i++) {
      fX[i]  = vx(i + ivxlow);
      fY[i]  = vy(i + ivylow);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Graph constructor with two vectors of doubles in input
/// A graph is build with the X coordinates taken from vx and Y coord from vy
/// The number of points in the graph is the minimum of number of points
/// in vx and vy.

TGraph::TGraph(const TVectorD &vx, const TVectorD &vy)
   : TNamed("Graph", "Graph"), TAttLine(), TAttFill(0, 1000), TAttMarker()
{
   fNpoints = TMath::Min(vx.GetNrows(), vy.GetNrows());
   if (!CtorAllocate()) return;
   Int_t ivxlow  = vx.GetLwb();
   Int_t ivylow  = vy.GetLwb();
   for (Int_t i = 0; i < fNpoints; i++) {
      fX[i]  = vx(i + ivxlow);
      fY[i]  = vy(i + ivylow);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Graph constructor importing its parameters from the TH1 object passed as argument

TGraph::TGraph(const TH1 *h)
   : TNamed("Graph", "Graph"), TAttLine(), TAttFill(0, 1000), TAttMarker()
{
   if (!h) {
      Error("TGraph", "Pointer to histogram is null");
      fNpoints = 0;
      return;
   }
   if (h->GetDimension() != 1) {
      Error("TGraph", "Histogram must be 1-D; h %s is %d-D", h->GetName(), h->GetDimension());
      fNpoints = 0;
   } else {
      fNpoints = ((TH1*)h)->GetXaxis()->GetNbins();
   }

   if (!CtorAllocate()) return;

   TAxis *xaxis = ((TH1*)h)->GetXaxis();
   for (Int_t i = 0; i < fNpoints; i++) {
      fX[i] = xaxis->GetBinCenter(i + 1);
      fY[i] = h->GetBinContent(i + 1);
   }
   h->TAttLine::Copy(*this);
   h->TAttFill::Copy(*this);
   h->TAttMarker::Copy(*this);

   std::string gname = "Graph_from_" + std::string(h->GetName());
   SetName(gname.c_str());
   SetTitle(h->GetTitle());
}

////////////////////////////////////////////////////////////////////////////////
/// Graph constructor importing its parameters from the TF1 object passed as argument
/// - if option =="" (default), a TGraph is created with points computed
///                at the fNpx points of f.
/// - if option =="d", a TGraph is created with points computed with the derivatives
///                at the fNpx points of f.
/// - if option =="i", a TGraph is created with points computed with the integral
///                at the fNpx points of f.
/// - if option =="I", a TGraph is created with points computed with the integral
///                at the fNpx+1 points of f and the integral is normalized to 1.

TGraph::TGraph(const TF1 *f, Option_t *option)
   : TNamed("Graph", "Graph"), TAttLine(), TAttFill(0, 1000), TAttMarker()
{
   char coption = ' ';
   if (!f) {
      Error("TGraph", "Pointer to function is null");
      fNpoints = 0;
   } else {
      fNpoints   = f->GetNpx();
      if (option) coption = *option;
      if (coption == 'i' || coption == 'I') fNpoints++;
   }
   if (!CtorAllocate()) return;

   Double_t xmin = f->GetXmin();
   Double_t xmax = f->GetXmax();
   Double_t dx   = (xmax - xmin) / fNpoints;
   Double_t integ = 0;
   Int_t i;
   for (i = 0; i < fNpoints; i++) {
      if (coption == 'i' || coption == 'I') {
         fX[i] = xmin + i * dx;
         if (i == 0) fY[i] = 0;
         else        fY[i] = integ + ((TF1*)f)->Integral(fX[i] - dx, fX[i]);
         integ = fY[i];
      } else if (coption == 'd' || coption == 'D') {
         fX[i] = xmin + (i + 0.5) * dx;
         fY[i] = ((TF1*)f)->Derivative(fX[i]);
      } else {
         fX[i] = xmin + (i + 0.5) * dx;
         fY[i] = ((TF1*)f)->Eval(fX[i]);
      }
   }
   if (integ != 0 && coption == 'I') {
      for (i = 1; i < fNpoints; i++) fY[i] /= integ;
   }

   f->TAttLine::Copy(*this);
   f->TAttFill::Copy(*this);
   f->TAttMarker::Copy(*this);

   SetName(f->GetName());
   SetTitle(f->GetTitle());
}

////////////////////////////////////////////////////////////////////////////////
/// Graph constructor reading input from filename.
///
/// `filename` is assumed to contain at least two columns of numbers.
/// the string format is by default `"%lg %lg"`.
/// This is a standard c formatting for `scanf()`.
///
/// If columns of numbers should be skipped, a `"%*lg"` or `"%*s"` for each column
/// can be added,  e.g. `"%lg %*lg %lg"` would read x-values from the first and
/// y-values from the third column.
///
/// For files separated by a specific delimiter different from ' ' and '\t' (e.g.
/// ';' in csv files) you can avoid using `%*s` to bypass this delimiter by explicitly
/// specify the `option` argument,
/// e.g. option=`" \t,;"` for columns of figures separated by any of these characters
/// (' ', '\t', ',', ';')
/// used once (e.g. `"1;1"`) or in a combined way (`" 1;,;;  1"`).
/// Note in that case, the instantiation is about two times slower.

TGraph::TGraph(const char *filename, const char *format, Option_t *option)
   : TNamed("Graph", filename), TAttLine(), TAttFill(0, 1000), TAttMarker()
{
   Double_t x, y;
   TString fname = filename;
   gSystem->ExpandPathName(fname);

   std::ifstream infile(fname.Data());
   if (!infile.good()) {
      MakeZombie();
      Error("TGraph", "Cannot open file: %s, TGraph is Zombie", filename);
      fNpoints = 0;
      return;
   } else {
      fNpoints = 100;  //initial number of points
   }
   if (!CtorAllocate()) return;
   std::string line;
   Int_t np = 0;

   // No delimiters specified (standard constructor).
   if (strcmp(option, "") == 0) {

      while (std::getline(infile, line, '\n')) {
         if (2 != sscanf(line.c_str(), format, &x, &y)) {
            continue; //skip empty and ill-formed lines
         }
         SetPoint(np, x, y);
         np++;
      }
      Set(np);

      // A delimiter has been specified in "option"
   } else {

      // Checking format and creating its boolean counterpart
      TString format_ = TString(format) ;
      format_.ReplaceAll(" ", "") ;
      format_.ReplaceAll("\t", "") ;
      format_.ReplaceAll("lg", "") ;
      format_.ReplaceAll("s", "") ;
      format_.ReplaceAll("%*", "0") ;
      format_.ReplaceAll("%", "1") ;
      if (!format_.IsDigit()) {
         Error("TGraph", "Incorrect input format! Allowed formats are {\"%%lg\",\"%%*lg\" or \"%%*s\"}");
         return;
      }
      Int_t ntokens = format_.Length() ;
      if (ntokens < 2) {
         Error("TGraph", "Incorrect input format! Only %d tag(s) in format whereas 2 \"%%lg\" tags are expected!", ntokens);
         return;
      }
      Int_t ntokensToBeSaved = 0 ;
      Bool_t * isTokenToBeSaved = new Bool_t [ntokens] ;
      for (Int_t idx = 0; idx < ntokens; idx++) {
         isTokenToBeSaved[idx] = TString::Format("%c", format_[idx]).Atoi() ; //atoi(&format_[idx]) does not work for some reason...
         if (isTokenToBeSaved[idx] == 1) {
            ntokensToBeSaved++ ;
         }
      }
      if (ntokens >= 2 && ntokensToBeSaved != 2) { //first condition not to repeat the previous error message
         Error("TGraph", "Incorrect input format! There are %d \"%%lg\" tag(s) in format whereas 2 and only 2 are expected!", ntokensToBeSaved);
         delete [] isTokenToBeSaved ;
         return;
      }

      // Initializing loop variables
      Bool_t isLineToBeSkipped = kFALSE ; //empty and ill-formed lines
      char * token = NULL ;
      TString token_str = "" ;
      Int_t token_idx = 0 ;
      Double_t * value = new Double_t [2] ; //x,y buffers
      Int_t value_idx = 0 ;

      // Looping
      char *rest;
      while (std::getline(infile, line, '\n')) {
         if (line != "") {
            if (line[line.size() - 1] == char(13)) {  // removing DOS CR character
               line.erase(line.end() - 1, line.end()) ;
            }
            //token = R__STRTOK_R(const_cast<char *>(line.c_str()), option, rest);
            token = R__STRTOK_R(const_cast<char *>(line.c_str()), option, &rest);
            while (token != NULL && value_idx < 2) {
               if (isTokenToBeSaved[token_idx]) {
                  token_str = TString(token) ;
                  token_str.ReplaceAll("\t", "") ;
                  if (!token_str.IsFloat()) {
                     isLineToBeSkipped = kTRUE ;
                     break ;
                  } else {
                     value[value_idx] = token_str.Atof() ;
                     value_idx++ ;
                  }
               }
               token = R__STRTOK_R(NULL, option, &rest); // next token
               token_idx++ ;
            }
            if (!isLineToBeSkipped && value_idx == 2) {
               x = value[0] ;
               y = value[1] ;
               SetPoint(np, x, y) ;
               np++ ;
            }
         }
         isLineToBeSkipped = kFALSE ;
         token = NULL ;
         token_idx = 0 ;
         value_idx = 0 ;
      }
      Set(np) ;

      // Cleaning
      delete [] isTokenToBeSaved ;
      delete [] value ;
      delete token ;
   }
   infile.close();
}

////////////////////////////////////////////////////////////////////////////////
/// Graph default destructor.

TGraph::~TGraph()
{
   delete [] fX;
   delete [] fY;
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
      fFunctions = 0; //to avoid accessing a deleted object in RecursiveRemove
   }
   delete fHistogram;
}

////////////////////////////////////////////////////////////////////////////////
/// Allocate internal data structures for `newsize` points.

Double_t **TGraph::Allocate(Int_t newsize) {
   return AllocateArrays(2, newsize);
}

////////////////////////////////////////////////////////////////////////////////
/// Allocate arrays.

Double_t** TGraph::AllocateArrays(Int_t Narrays, Int_t arraySize)
{
   if (arraySize < 0) {
      arraySize = 0;
   }
   Double_t **newarrays = new Double_t*[Narrays];
   if (!arraySize) {
      for (Int_t i = 0; i < Narrays; ++i)
         newarrays[i] = 0;
   } else {
      for (Int_t i = 0; i < Narrays; ++i)
         newarrays[i] = new Double_t[arraySize];
   }
   fMaxSize = arraySize;
   return newarrays;
}

////////////////////////////////////////////////////////////////////////////////
/// Apply function f to all the data points
/// f may be a 1-D function TF1 or 2-d function TF2
/// The Y values of the graph are replaced by the new values computed
/// using the function

void TGraph::Apply(TF1 *f)
{
   if (fHistogram) SetBit(kResetHisto);

   for (Int_t i = 0; i < fNpoints; i++) {
      fY[i] = f->Eval(fX[i], fY[i]);
   }
   if (gPad) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Browse

void TGraph::Browse(TBrowser *b)
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
/// Return the chisquare of this graph with respect to f1.
/// The chisquare is computed as the sum of the quantity below at each point:
/// \f[
///   \frac{(y-f1(x))^{2}}{ey^{2}+(\frac{1}{2}(exl+exh)f1'(x))^{2}}
/// \f]
/// where x and y are the graph point coordinates and f1'(x) is the derivative of function f1(x).
/// This method to approximate the uncertainty in y because of the errors in x, is called
/// "effective variance" method.
/// In case of a pure TGraph, the denominator is 1.
/// In case of a TGraphErrors or TGraphAsymmErrors the errors are taken
/// into account.
/// By default the range of the graph is used whatever function range.
/// Use option "R" to use the function range

Double_t TGraph::Chisquare(TF1 *func, Option_t * option) const
{
   if (!func) {
      Error("Chisquare","Function pointer is Null - return -1");
      return -1;
   }

   TString opt(option); opt.ToUpper();
   bool useRange = opt.Contains("R");

   return ROOT::Fit::Chisquare(*this, *func,useRange);
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if point number "left"'s argument (angle with respect to positive
/// x-axis) is bigger than that of point number "right". Can be used by Sort.

Bool_t TGraph::CompareArg(const TGraph* gr, Int_t left, Int_t right)
{
   Double_t xl = 0, yl = 0, xr = 0, yr = 0;
   gr->GetPoint(left, xl, yl);
   gr->GetPoint(right, xr, yr);
   return (TMath::ATan2(yl, xl) > TMath::ATan2(yr, xr));
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if fX[left] > fX[right]. Can be used by Sort.

Bool_t TGraph::CompareX(const TGraph* gr, Int_t left, Int_t right)
{
   return gr->fX[left] > gr->fX[right];
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if fY[left] > fY[right]. Can be used by Sort.

Bool_t TGraph::CompareY(const TGraph* gr, Int_t left, Int_t right)
{
   return gr->fY[left] > gr->fY[right];
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if point number "left"'s distance to origin is bigger than
/// that of point number "right". Can be used by Sort.

Bool_t TGraph::CompareRadius(const TGraph* gr, Int_t left, Int_t right)
{
   return gr->fX[left] * gr->fX[left] + gr->fY[left] * gr->fY[left]
          > gr->fX[right] * gr->fX[right] + gr->fY[right] * gr->fY[right];
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the x/y range of the points in this graph

void TGraph::ComputeRange(Double_t &xmin, Double_t &ymin, Double_t &xmax, Double_t &ymax) const
{
   if (fNpoints <= 0) {
      xmin = xmax = ymin = ymax = 0;
      return;
   }
   xmin = xmax = fX[0];
   ymin = ymax = fY[0];

   Double_t xminl = 0; // Positive minimum. Used in case of log scale along X axis.
   Double_t yminl = 0; // Positive minimum. Used in case of log scale along Y axis.

   for (Int_t i = 1; i < fNpoints; i++) {
      if (fX[i] < xmin) xmin = fX[i];
      if (fX[i] > xmax) xmax = fX[i];
      if (fY[i] < ymin) ymin = fY[i];
      if (fY[i] > ymax) ymax = fY[i];
      if (ymin>0 && (yminl==0 || ymin<yminl)) yminl = ymin;
      if (xmin>0 && (xminl==0 || xmin<xminl)) xminl = xmin;
   }

   if (gPad && gPad->GetLogy() && yminl>0) ymin = yminl;
   if (gPad && gPad->GetLogx() && xminl>0) xmin = xminl;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy points from fX and fY to arrays[0] and arrays[1]
/// or to fX and fY if arrays == 0 and ibegin != iend.
/// If newarrays is non null, replace fX, fY with pointers from newarrays[0,1].
/// Delete newarrays, old fX and fY

void TGraph::CopyAndRelease(Double_t **newarrays, Int_t ibegin, Int_t iend,
                            Int_t obegin)
{
   CopyPoints(newarrays, ibegin, iend, obegin);
   if (newarrays) {
      delete[] fX;
      fX = newarrays[0];
      delete[] fY;
      fY = newarrays[1];
      delete[] newarrays;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy points from fX and fY to arrays[0] and arrays[1]
/// or to fX and fY if arrays == 0 and ibegin != iend.

Bool_t TGraph::CopyPoints(Double_t **arrays, Int_t ibegin, Int_t iend,
                          Int_t obegin)
{
   if (ibegin < 0 || iend <= ibegin || obegin < 0) { // Error;
      return kFALSE;
   }
   if (!arrays && ibegin == obegin) { // No copying is needed
      return kFALSE;
   }
   Int_t n = (iend - ibegin) * sizeof(Double_t);
   if (arrays) {
      memmove(&arrays[0][obegin], &fX[ibegin], n);
      memmove(&arrays[1][obegin], &fY[ibegin], n);
   } else {
      memmove(&fX[obegin], &fX[ibegin], n);
      memmove(&fY[obegin], &fY[ibegin], n);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// In constructors set fNpoints than call this method.
/// Return kFALSE if the graph will contain no points.
///Note: This function should be called only from the constructor
/// since it does not delete previously existing arrays

Bool_t TGraph::CtorAllocate()
{
   fHistogram = 0;
   fMaximum = -1111;
   fMinimum = -1111;
   SetBit(kClipFrame);
   fFunctions = new TList;
   if (fNpoints <= 0) {
      fNpoints = 0;
      fMaxSize   = 0;
      fX         = 0;
      fY         = 0;
      return kFALSE;
   } else {
      fMaxSize   = fNpoints;
      fX = new Double_t[fMaxSize];
      fY = new Double_t[fMaxSize];
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this graph with its current attributes.
///
/// The options to draw a graph are described in TGraphPainter class.

void TGraph::Draw(Option_t *option)
{
   TString opt = option;
   opt.ToLower();

   if (opt.Contains("same")) {
      opt.ReplaceAll("same", "");
   }

   // in case of option *, set marker style to 3 (star) and replace
   // * option by option P.
   Ssiz_t pos;
   if ((pos = opt.Index("*")) != kNPOS) {
      SetMarkerStyle(3);
      opt.Replace(pos, 1, "p");
   }

   // If no option is specified, it is defined as "alp" in case there
   // no current pad or if the current pad as no axis defined.
   if (!strlen(option)) {
      if (gPad) {
         if (!gPad->GetListOfPrimitives()->FindObject("TFrame")) opt = "alp";
      } else {
         opt = "alp";
      }
   }

   if (gPad) {
      if (!gPad->IsEditable()) gROOT->MakeDefCanvas();
      if (opt.Contains("a")) gPad->Clear();
   }

   AppendPad(opt);

   gPad->IncrementPaletteColor(1, opt);

}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a graph.
///
///  Compute the closest distance of approach from point px,py to this line.
///  The distance is computed in pixels units.

Int_t TGraph::DistancetoPrimitive(Int_t px, Int_t py)
{
   TVirtualGraphPainter *painter = TVirtualGraphPainter::GetPainter();
   if (painter) return painter->DistancetoPrimitiveHelper(this, px, py);
   else return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this graph with new attributes.

void TGraph::DrawGraph(Int_t n, const Int_t *x, const Int_t *y, Option_t *option)
{
   TGraph *newgraph = new TGraph(n, x, y);
   TAttLine::Copy(*newgraph);
   TAttFill::Copy(*newgraph);
   TAttMarker::Copy(*newgraph);
   newgraph->SetBit(kCanDelete);
   newgraph->AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this graph with new attributes.

void TGraph::DrawGraph(Int_t n, const Float_t *x, const Float_t *y, Option_t *option)
{
   TGraph *newgraph = new TGraph(n, x, y);
   TAttLine::Copy(*newgraph);
   TAttFill::Copy(*newgraph);
   TAttMarker::Copy(*newgraph);
   newgraph->SetBit(kCanDelete);
   newgraph->AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this graph with new attributes.

void TGraph::DrawGraph(Int_t n, const Double_t *x, const Double_t *y, Option_t *option)
{
   const Double_t *xx = x;
   const Double_t *yy = y;
   if (!xx) xx = fX;
   if (!yy) yy = fY;
   TGraph *newgraph = new TGraph(n, xx, yy);
   TAttLine::Copy(*newgraph);
   TAttFill::Copy(*newgraph);
   TAttMarker::Copy(*newgraph);
   newgraph->SetBit(kCanDelete);
   newgraph->AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Display a panel with all graph drawing options.

void TGraph::DrawPanel()
{
   TVirtualGraphPainter *painter = TVirtualGraphPainter::GetPainter();
   if (painter) painter->DrawPanelHelper(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Interpolate points in this graph at x using a TSpline.
///
///  - if spline==0 and option="" a linear interpolation between the two points
///    close to x is computed. If x is outside the graph range, a linear
///    extrapolation is computed.
///  - if spline==0 and option="S" a TSpline3 object is created using this graph
///    and the interpolated value from the spline is returned.
///    the internally created spline is deleted on return.
///  - if spline is specified, it is used to return the interpolated value.
///
///   If the points are sorted in X a binary search is used (significantly faster)
///   One needs to set the bit  TGraph::SetBit(TGraph::kIsSortedX) before calling
///   TGraph::Eval to indicate that the graph is sorted in X.

Double_t TGraph::Eval(Double_t x, TSpline *spline, Option_t *option) const
{

   if (spline) {
      //spline interpolation using the input spline
      return spline->Eval(x);
   }

   if (fNpoints == 0) return 0;
   if (fNpoints == 1) return fY[0];

   if (option && *option) {
      TString opt = option;
      opt.ToLower();
      // create a TSpline every time when using option "s" and no spline pointer is given
      if (opt.Contains("s")) {

         // points must be sorted before using a TSpline
         std::vector<Double_t> xsort(fNpoints);
         std::vector<Double_t> ysort(fNpoints);
         std::vector<Int_t> indxsort(fNpoints);
         TMath::Sort(fNpoints, fX, &indxsort[0], false);
         for (Int_t i = 0; i < fNpoints; ++i) {
            xsort[i] = fX[ indxsort[i] ];
            ysort[i] = fY[ indxsort[i] ];
         }

         // spline interpolation creating a new spline
         TSpline3 s("", &xsort[0], &ysort[0], fNpoints);
         Double_t result = s.Eval(x);
         return result;
      }
   }
   //linear interpolation
   //In case x is < fX[0] or > fX[fNpoints-1] return the extrapolated point

   //find points in graph around x assuming points are not sorted
   // (if point are sorted use a binary search)
   Int_t low  = -1;
   Int_t up  = -1;
   if (TestBit(TGraph::kIsSortedX) ) {
      low = TMath::BinarySearch(fNpoints, fX, x);
      if (low == -1)  {
         // use first two points for doing an extrapolation
         low = 0;
      }
      if (fX[low] == x) return fY[low];
      if (low == fNpoints-1) low--; // for extrapolating
      up = low+1;
   }
   else {
      // case TGraph is not sorted

   // find neighbours simply looping  all points
   // and find also the 2 adjacent points: (low2 < low < x < up < up2 )
   // needed in case x is outside the graph ascissa interval
      Int_t low2 = -1;
      Int_t up2 = -1;

      for (Int_t i = 0; i < fNpoints; ++i) {
         if (fX[i] < x) {
            if (low == -1 || fX[i] > fX[low])  {
               low2 = low;
               low = i;
            } else if (low2 == -1) low2 = i;
         } else if (fX[i] > x) {
            if (up  == -1 || fX[i] < fX[up])  {
               up2 = up;
               up = i;
            } else if (up2 == -1) up2 = i;
         } else // case x == fX[i]
            return fY[i]; // no interpolation needed
      }

      // treat cases when x is outside graph min max abscissa
      if (up == -1)  {
         up  = low;
         low = low2;
      }
      if (low == -1) {
         low = up;
         up  = up2;
      }
   }
   // do now the linear interpolation
   assert(low != -1 && up != -1);

   if (fX[low] == fX[up]) return fY[low];
   Double_t yn = fY[up] + (x - fX[up]) * (fY[low] - fY[up]) / (fX[low] - fX[up]);
   return yn;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.
///
///  This member function is called when a graph is clicked with the locator
///
///  If Left button clicked on one of the line end points, this point
///     follows the cursor until button is released.
///
///  if Middle button clicked, the line is moved parallel to itself
///     until the button is released.

void TGraph::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   TVirtualGraphPainter *painter = TVirtualGraphPainter::GetPainter();
   if (painter) painter->ExecuteEventHelper(this, event, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// If array sizes <= newsize, expand storage to 2*newsize.

void TGraph::Expand(Int_t newsize)
{
   Double_t **ps = ExpandAndCopy(newsize, fNpoints);
   CopyAndRelease(ps, 0, 0, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// If graph capacity is less than newsize points then make array sizes
/// equal to least multiple of step to contain newsize points.

void TGraph::Expand(Int_t newsize, Int_t step)
{
   if (newsize <= fMaxSize) {
      return;
   }
   Double_t **ps = Allocate(step * (newsize / step + (newsize % step ? 1 : 0)));
   CopyAndRelease(ps, 0, fNpoints, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// if size > fMaxSize allocate new arrays of 2*size points and copy iend first
/// points.
/// Return pointer to new arrays.

Double_t **TGraph::ExpandAndCopy(Int_t size, Int_t iend)
{
   if (size <= fMaxSize) {
      return 0;
   }
   Double_t **newarrays = Allocate(2 * size);
   CopyPoints(newarrays, 0, iend, 0);
   return newarrays;
}

////////////////////////////////////////////////////////////////////////////////
/// Set zero values for point arrays in the range [begin, end)
/// Should be redefined in descendant classes

void TGraph::FillZero(Int_t begin, Int_t end, Bool_t)
{
   memset(fX + begin, 0, (end - begin)*sizeof(Double_t));
   memset(fY + begin, 0, (end - begin)*sizeof(Double_t));
}

////////////////////////////////////////////////////////////////////////////////
/// Search object named name in the list of functions

TObject *TGraph::FindObject(const char *name) const
{
   if (fFunctions) return fFunctions->FindObject(name);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Search object obj in the list of functions

TObject *TGraph::FindObject(const TObject *obj) const
{
   if (fFunctions) return fFunctions->FindObject(obj);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Fit this graph with function with name fname.
///
/// interface to TGraph::Fit(TF1 *f1...
///
/// fname is the name of an already predefined function created by TF1 or TF2
/// Predefined functions such as gaus, expo and poln are automatically
/// created by ROOT.
///
/// fname can also be a formula, accepted by the linear fitter (linear parts divided
/// by "++" sign), for example "x++sin(x)" for fitting "[0]*x+[1]*sin(x)"

TFitResultPtr TGraph::Fit(const char *fname, Option_t *option, Option_t *, Axis_t xmin, Axis_t xmax)
{
   char *linear;
   linear = (char*) strstr(fname, "++");
   if (linear) {
      TF1 f1(fname, fname, xmin, xmax);
      return Fit(&f1, option, "", xmin, xmax);
   }
   TF1 * f1 = (TF1*)gROOT->GetFunction(fname);
   if (!f1) {
      Printf("Unknown function: %s", fname);
      return -1;
   }
   return Fit(f1, option, "", xmin, xmax);
}

////////////////////////////////////////////////////////////////////////////////
/// Fit this graph with function f1.
///
/// f1 is an already predefined function created by TF1.
/// Predefined functions such as gaus, expo and poln are automatically
/// created by ROOT.
///
/// The list of fit options is given in parameter option.
///
/// option | description
/// -------|------------
/// "W" | Ignore all point errors when fitting a TGraphErrors or TGraphAsymmErrors
/// "U" | Use a User specified fitting algorithm (via SetFCN)
/// "Q" | Quiet mode (minimum printing)
/// "V" | Verbose mode (default is between Q and V)
/// "E" | Perform better Errors estimation using Minos technique
/// "B" | User defined parameter settings are used for predefined functions like "gaus", "expo", "poln", "landau". Use this option when you want to fix one or more parameters for these functions.
/// "M" | More. Improve fit results. It uses the IMPROVE command of TMinuit (see TMinuit::mnimpr). This algorithm attempts to improve the found local minimum by searching for a better one.
/// "R" | Use the Range specified in the function range
/// "N" | Do not store the graphics function, do not draw
/// "0" | Do not plot the result of the fit. By default the fitted function is drawn unless the option "N" above is specified.
/// "+" | Add this new fitted function to the list of fitted functions (by default, any previous function is deleted)
/// "C" | In case of linear fitting, do not calculate the chisquare (saves time)
/// "F" | If fitting a polN, use the minuit fitter
/// "EX0" | When fitting a TGraphErrors or TGraphAsymErrors do not consider errors in the X coordinates
/// "ROB" | In case of linear fitting, compute the LTS regression coefficients (robust (resistant) regression), using the default fraction of good points "ROB=0.x" - compute the LTS regression coefficients, using 0.x as a fraction of good points
/// "S" |  The result of the fit is returned in the TFitResultPtr (see below Access to the Fit Result)
///
/// When the fit is drawn (by default), the parameter goption may be used
/// to specify a list of graphics options. See TGraphPainter for a complete
/// list of these options.
///
/// In order to use the Range option, one must first create a function
/// with the expression to be fitted. For example, if your graph
/// has a defined range between -4 and 4 and you want to fit a gaussian
/// only in the interval 1 to 3, you can do:
///
///        TF1 *f1 = new TF1("f1","gaus",1,3);
///        graph->Fit("f1","R");
///
/// Who is calling this function:
///
/// Note that this function is called when calling TGraphErrors::Fit
/// or TGraphAsymmErrors::Fit ot TGraphBentErrors::Fit
/// See the discussion below on error calculation.
///
/// ### Linear fitting:
///   When the fitting function is linear (contains the "++" sign) or the fitting
///   function is a polynomial, a linear fitter is initialised.
///   To create a linear function, use the following syntax: linear parts
///   separated by "++" sign.
///   Example: to fit the parameters of "[0]*x + [1]*sin(x)", create a
///    TF1 *f1=new TF1("f1", "x++sin(x)", xmin, xmax);
///   For such a TF1 you don't have to set the initial conditions.
///   Going via the linear fitter for functions, linear in parameters, gives a
///   considerable advantage in speed.
///
/// ### Setting initial conditions:
///
///   Parameters must be initialized before invoking the Fit function.
///   The setting of the parameter initial values is automatic for the
///   predefined functions : poln, expo, gaus, landau. One can however disable
///   this automatic computation by specifying the option "B".
///   You can specify boundary limits for some or all parameters via
///
///        f1->SetParLimits(p_number, parmin, parmax);
///   If parmin>=parmax, the parameter is fixed
///   Note that you are not forced to fix the limits for all parameters.
///   For example, if you fit a function with 6 parameters, you can do:
///
///     func->SetParameters(0,3.1,1.e-6,0.1,-8,100);
///     func->SetParLimits(4,-10,-4);
///     func->SetParLimits(5, 1,1);
///   With this setup, parameters 0->3 can vary freely.
///   Parameter 4 has boundaries [-10,-4] with initial value -8.
///   Parameter 5 is fixed to 100.
///
/// ### Fit range:
///
///   The fit range can be specified in two ways:
///     - specify rxmax > rxmin (default is rxmin=rxmax=0)
///     - specify the option "R". In this case, the function will be taken
///       instead of the full graph range.
///
/// ### Changing the fitting function:
///
///   By default a chi2 fitting function is used for fitting a TGraph.
///   The function is implemented in FitUtil::EvaluateChi2.
///   In case of TGraphErrors an effective chi2 is used (see below TGraphErrors fit)
///   To specify a User defined fitting function, specify option "U" and
///   call the following functions:
///
///      TVirtualFitter::Fitter(mygraph)->SetFCN(MyFittingFunction)
///   where MyFittingFunction is of type:
///   extern void MyFittingFunction(Int_t &npar, Double_t *gin, Double_t &f,
///                                 Double_t *u, Int_t flag);
///
///
/// ### TGraphErrors fit:
///
///   In case of a TGraphErrors object, when x errors are present, the error along x,
///   is projected along the y-direction by calculating the function at the points x-exlow and
///   x+exhigh. The chisquare is then computed as the sum of the quantity below at each point:
///
/// \f[
///   \frac{(y-f(x))^{2}}{ey^{2}+(\frac{1}{2}(exl+exh)f'(x))^{2}}
/// \f]
///
///   where x and y are the point coordinates, and f'(x) is the derivative of the
///   function f(x).
///
///   In case the function lies below (above) the data point, ey is ey_low (ey_high).
///
///   thanks to Andy Haas (haas@yahoo.com) for adding the case with TGraphAsymmErrors
///             University of Washington
///
///   The approach used to approximate the uncertainty in y because of the
///   errors in x is to make it equal the error in x times the slope of the line.
///   The improvement, compared to the first method (f(x+ exhigh) - f(x-exlow))/2
///   is of (error of x)**2 order. This approach is called "effective variance method".
///   This improvement has been made in version 4.00/08 by Anna Kreshuk.
///   The implementation is provided in the function FitUtil::EvaluateChi2Effective
///
/// NOTE:
/// 1. By using the "effective variance" method a simple linear regression
///    becomes a non-linear case, which takes several iterations
///    instead of 0 as in the linear case.
/// 2. The effective variance technique assumes that there is no correlation
///    between the x and y coordinate.
/// 3. The standard chi2 (least square) method without error in the coordinates (x) can
///    be forced by using option "EX0"
/// 4. The linear fitter doesn't take into account the errors in x. When fitting a
///    TGraphErrors with a linear functions the errors in x will not be considered.
///    If errors in x are important, go through minuit (use option "F" for polynomial fitting).
/// 5. When fitting a TGraph (i.e. no errors associated with each point),
///    a correction is applied to the errors on the parameters with the following
///    formula: errorp *= sqrt(chisquare/(ndf-1))
///
///  ## Access to the fit result
///  The function returns a TFitResultPtr which can hold a  pointer to a TFitResult object.
///  By default the TFitResultPtr contains only the status of the fit which is return by an
///  automatic conversion of the TFitResultPtr to an integer. One can write in this case
///  directly:
///
///      Int_t fitStatus =  h->Fit(myFunc)
///
///  If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves
///  as a smart pointer to it. For example one can do:
///
///      TFitResultPtr r = h->Fit(myFunc,"S");
///      TMatrixDSym cov = r->GetCovarianceMatrix();  //  to access the covariance matrix
///      Double_t chi2   = r->Chi2(); // to retrieve the fit chi2
///      Double_t par0   = r->Value(0); // retrieve the value for the parameter 0
///      Double_t err0   = r->ParError(0); // retrieve the error for the parameter 0
///      r->Print("V");     // print full information of fit including covariance matrix
///      r->Write();        // store the result in a file
///
///  The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also
///  from the fitted function.
///  If the histogram is made persistent, the list of
///  associated functions is also persistent. Given a pointer (see above)
///  to an associated function myfunc, one can retrieve the function/fit
///  parameters with calls such as:
///
///      Double_t chi2 = myfunc->GetChisquare();
///      Double_t par0 = myfunc->GetParameter(0); //value of 1st parameter
///      Double_t err0 = myfunc->GetParError(0);  //error on first parameter
///
///
/// ### Access to the fit status
///  The status of the fit can be obtained converting the TFitResultPtr to an integer
///  independently if the fit option "S" is used or not:
///
///       TFitResultPtr r = h->Fit(myFunc,opt);
///       Int_t fitStatus = r;
///
///  The fitStatus is 0 if the fit is OK (i.e. no error occurred).
///  The value of the fit status code is negative in case of an error not connected with the
///  minimization procedure, for example when a wrong function is used.
///  Otherwise the return value is the one returned from the minimization procedure.
///  When TMinuit (default case) or Minuit2 are used as minimizer the status returned is :
///  fitStatus =  migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult.
///  TMinuit will return 0 (for migrad, minos, hesse or improve) in case of success and 4 in
///  case of error (see the documentation of TMinuit::mnexcm). So for example, for an error
///  only in Minos but not in Migrad a fitStatus of 40 will be returned.
///  Minuit2 will return also 0 in case of success and different values in migrad, minos or
///  hesse depending on the error.   See in this case the documentation of
///  Minuit2Minimizer::Minimize for the migradResult, Minuit2Minimizer::GetMinosError for the
///  minosResult and Minuit2Minimizer::Hesse for the hesseResult.
///  If other minimizers are used see their specific documentation for the status code
///  returned. For example in the case of Fumili, for the status returned see TFumili::Minimize.
///
/// ### Associated functions:
///   One or more object (typically a TF1*) can be added to the list
///   of functions (fFunctions) associated with each graph.
///   When TGraph::Fit is invoked, the fitted function is added to this list.
///   Given a graph gr, one can retrieve an associated function
///   with:  TF1 *myfunc = gr->GetFunction("myfunc");
///
///   If the graph is made persistent, the list of associated functions is also
///   persistent. Given a pointer (see above) to an associated function myfunc,
///   one can retrieve the function/fit parameters with calls such as:
///
///       Double_t chi2 = myfunc->GetChisquare();
///       Double_t par0 = myfunc->GetParameter(0); //value of 1st parameter
///       Double_t err0 = myfunc->GetParError(0);  //error on first parameter
///
/// ### Fit Statistics
///   You can change the statistics box to display the fit parameters with
///   the TStyle::SetOptFit(mode) method. This mode has four digits.
///   mode = pcev  (default = 0111)
///
///       v = 1;  print name/values of parameters
///       e = 1;  print errors (if e=1, v must be 1)
///       c = 1;  print Chisquare/Number of degrees of freedom
///       p = 1;  print Probability
///
///   For example: gStyle->SetOptFit(1011);
///   prints the fit probability, parameter names/values, and errors.
///   You can change the position of the statistics box with these lines
///   (where g is a pointer to the TGraph):
///
///       Root > TPaveStats *st = (TPaveStats*)g->GetListOfFunctions()->FindObject("stats")
///       Root > st->SetX1NDC(newx1); //new x start position
///       Root > st->SetX2NDC(newx2); //new x end position
///

TFitResultPtr TGraph::Fit(TF1 *f1, Option_t *option, Option_t *goption, Axis_t rxmin, Axis_t rxmax)
{
   Foption_t fitOption;
   ROOT::Fit::FitOptionsMake(ROOT::Fit::kGraph, option, fitOption);
   // create range and minimizer options with default values
   ROOT::Fit::DataRange range(rxmin, rxmax);
   ROOT::Math::MinimizerOptions minOption;
   return ROOT::Fit::FitObject(this, f1 , fitOption , minOption, goption, range);
}

////////////////////////////////////////////////////////////////////////////////
/// Display a GUI panel with all graph fit options.
///
/// See class TFitEditor for example

void TGraph::FitPanel()
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
   } else
      Error("FitPanel", "Unable to find the FitPanel plug-in");
}

////////////////////////////////////////////////////////////////////////////////
/// Return graph correlation factor

Double_t TGraph::GetCorrelationFactor() const
{
   Double_t rms1 = GetRMS(1);
   if (rms1 == 0) return 0;
   Double_t rms2 = GetRMS(2);
   if (rms2 == 0) return 0;
   return GetCovariance() / rms1 / rms2;
}

////////////////////////////////////////////////////////////////////////////////
/// Return covariance of vectors x,y

Double_t TGraph::GetCovariance() const
{
   if (fNpoints <= 0) return 0;
   Double_t sum = fNpoints, sumx = 0, sumy = 0, sumxy = 0;

   for (Int_t i = 0; i < fNpoints; i++) {
      sumx  += fX[i];
      sumy  += fY[i];
      sumxy += fX[i] * fY[i];
   }
   return sumxy / sum - sumx / sum * sumy / sum;
}

////////////////////////////////////////////////////////////////////////////////
/// Return mean value of X (axis=1)  or Y (axis=2)

Double_t TGraph::GetMean(Int_t axis) const
{
   if (axis < 1 || axis > 2) return 0;
   if (fNpoints <= 0) return 0;
   Double_t sumx = 0;
   for (Int_t i = 0; i < fNpoints; i++) {
      if (axis == 1) sumx += fX[i];
      else           sumx += fY[i];
   }
   return sumx / fNpoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Return RMS of X (axis=1)  or Y (axis=2)

Double_t TGraph::GetRMS(Int_t axis) const
{
   if (axis < 1 || axis > 2) return 0;
   if (fNpoints <= 0) return 0;
   Double_t sumx = 0, sumx2 = 0;
   for (Int_t i = 0; i < fNpoints; i++) {
      if (axis == 1) {
         sumx += fX[i];
         sumx2 += fX[i] * fX[i];
      } else           {
         sumx += fY[i];
         sumx2 += fY[i] * fY[i];
      }
   }
   Double_t x = sumx / fNpoints;
   Double_t rms2 = TMath::Abs(sumx2 / fNpoints - x * x);
   return TMath::Sqrt(rms2);
}

////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It always returns a negative value. Real implementation in TGraphErrors

Double_t TGraph::GetErrorX(Int_t) const
{
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It always returns a negative value. Real implementation in TGraphErrors

Double_t TGraph::GetErrorY(Int_t) const
{
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It always returns a negative value. Real implementation in TGraphErrors
/// and TGraphAsymmErrors

Double_t TGraph::GetErrorXhigh(Int_t) const
{
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It always returns a negative value. Real implementation in TGraphErrors
/// and TGraphAsymmErrors

Double_t TGraph::GetErrorXlow(Int_t) const
{
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It always returns a negative value. Real implementation in TGraphErrors
/// and TGraphAsymmErrors

Double_t TGraph::GetErrorYhigh(Int_t) const
{
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It always returns a negative value. Real implementation in TGraphErrors
/// and TGraphAsymmErrors

Double_t TGraph::GetErrorYlow(Int_t) const
{
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to function with name.
///
/// Functions such as TGraph::Fit store the fitted function in the list of
/// functions of this graph.

TF1 *TGraph::GetFunction(const char *name) const
{
   if (!fFunctions) return 0;
   return (TF1*)fFunctions->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a pointer to the histogram used to draw the axis
/// Takes into account the two following cases.
///  1. option 'A' was specified in TGraph::Draw. Return fHistogram
///  2. user had called TPad::DrawFrame. return pointer to hframe histogram

TH1F *TGraph::GetHistogram() const
{
   Double_t rwxmin, rwxmax, rwymin, rwymax, maximum, minimum, dx, dy;
   Double_t uxmin, uxmax;

   ComputeRange(rwxmin, rwymin, rwxmax, rwymax);  //this is redefined in TGraphErrors

   // (if fHistogram exist) && (if the log scale is on) &&
   // (if the computed range minimum is > 0) && (if the fHistogram minimum is zero)
   // then it means fHistogram limits have been computed in linear scale
   // therefore they might be too strict and cut some points. In that case the
   // fHistogram limits should be recomputed ie: the existing fHistogram
   // should not be returned.
   TH1F *historg = 0;
   if (fHistogram) {
      if (!TestBit(kResetHisto)) {
         if (gPad && gPad->GetLogx()) {
            if (rwxmin <= 0 || fHistogram->GetXaxis()->GetXmin() != 0) return fHistogram;
         } else if (gPad && gPad->GetLogy()) {
            if (rwymin <= 0 || fHistogram->GetMinimum() != 0) return fHistogram;
         } else {
            return fHistogram;
         }
      } else {
         const_cast <TGraph*>(this)->ResetBit(kResetHisto);
      }
      historg = fHistogram;
   }

   if (rwxmin == rwxmax) rwxmax += 1.;
   if (rwymin == rwymax) rwymax += 1.;
   dx = 0.1 * (rwxmax - rwxmin);
   dy = 0.1 * (rwymax - rwymin);
   uxmin    = rwxmin - dx;
   uxmax    = rwxmax + dx;
   minimum  = rwymin - dy;
   maximum  = rwymax + dy;

   if (fMinimum != -1111) minimum = fMinimum;
   if (fMaximum != -1111) maximum = fMaximum;

   // the graph is created with at least as many channels as there are points
   // to permit zooming on the full range
   if (uxmin < 0 && rwxmin >= 0) {
      if (gPad && gPad->GetLogx()) uxmin = 0.9 * rwxmin;
      else                         uxmin = 0;
   }
   if (uxmax > 0 && rwxmax <= 0) {
      if (gPad && gPad->GetLogx()) uxmax = 1.1 * rwxmax;
      else                         uxmax = 0;
   }

   if (minimum < 0 && rwymin >= 0) minimum = 0.9 * rwymin;

   if (minimum <= 0 && gPad && gPad->GetLogy()) minimum = 0.001 * maximum;
   if (uxmin <= 0 && gPad && gPad->GetLogx()) {
      if (uxmax > 1000) uxmin = 1;
      else              uxmin = 0.001 * uxmax;
   }

   rwxmin = uxmin;
   rwxmax = uxmax;
   Int_t npt = 100;
   if (fNpoints > npt) npt = fNpoints;
   const char *gname = GetName();
   if (!gname[0]) gname = "Graph";
   // do not add the histogram to gDirectory
   // use local TDirectory::TContect that will set temporarly gDirectory to a nullptr and
   // will avoid that histogram is added in the global directory
   {
      TDirectory::TContext ctx(nullptr);
      ((TGraph*)this)->fHistogram = new TH1F(gname, GetTitle(), npt, rwxmin, rwxmax);
   }
   if (!fHistogram) return 0;
   fHistogram->SetMinimum(minimum);
   fHistogram->SetBit(TH1::kNoStats);
   fHistogram->SetMaximum(maximum);
   fHistogram->GetYaxis()->SetLimits(minimum, maximum);
   // Restore the axis attributes if needed
   if (historg) {
      fHistogram->GetXaxis()->SetTitle(historg->GetXaxis()->GetTitle());
      fHistogram->GetXaxis()->CenterTitle(historg->GetXaxis()->GetCenterTitle());
      fHistogram->GetXaxis()->RotateTitle(historg->GetXaxis()->GetRotateTitle());
      fHistogram->GetXaxis()->SetNoExponent(historg->GetXaxis()->GetNoExponent());
      fHistogram->GetXaxis()->SetNdivisions(historg->GetXaxis()->GetNdivisions());
      fHistogram->GetXaxis()->SetLabelFont(historg->GetXaxis()->GetLabelFont());
      fHistogram->GetXaxis()->SetLabelOffset(historg->GetXaxis()->GetLabelOffset());
      fHistogram->GetXaxis()->SetLabelSize(historg->GetXaxis()->GetLabelSize());
      fHistogram->GetXaxis()->SetTitleSize(historg->GetXaxis()->GetTitleSize());
      fHistogram->GetXaxis()->SetTitleOffset(historg->GetXaxis()->GetTitleOffset());
      fHistogram->GetXaxis()->SetTitleFont(historg->GetXaxis()->GetTitleFont());
      fHistogram->GetXaxis()->SetTimeDisplay(historg->GetXaxis()->GetTimeDisplay());
      fHistogram->GetXaxis()->SetTimeFormat(historg->GetXaxis()->GetTimeFormat());

      fHistogram->GetYaxis()->SetTitle(historg->GetYaxis()->GetTitle());
      fHistogram->GetYaxis()->CenterTitle(historg->GetYaxis()->GetCenterTitle());
      fHistogram->GetYaxis()->RotateTitle(historg->GetYaxis()->GetRotateTitle());
      fHistogram->GetYaxis()->SetNoExponent(historg->GetYaxis()->GetNoExponent());
      fHistogram->GetYaxis()->SetNdivisions(historg->GetYaxis()->GetNdivisions());
      fHistogram->GetYaxis()->SetLabelFont(historg->GetYaxis()->GetLabelFont());
      fHistogram->GetYaxis()->SetLabelOffset(historg->GetYaxis()->GetLabelOffset());
      fHistogram->GetYaxis()->SetLabelSize(historg->GetYaxis()->GetLabelSize());
      fHistogram->GetYaxis()->SetTitleSize(historg->GetYaxis()->GetTitleSize());
      fHistogram->GetYaxis()->SetTitleOffset(historg->GetYaxis()->GetTitleOffset());
      fHistogram->GetYaxis()->SetTitleFont(historg->GetYaxis()->GetTitleFont());
      fHistogram->GetYaxis()->SetTimeDisplay(historg->GetYaxis()->GetTimeDisplay());
      fHistogram->GetYaxis()->SetTimeFormat(historg->GetYaxis()->GetTimeFormat());
      delete historg;
   }
   return fHistogram;
}

////////////////////////////////////////////////////////////////////////////////
/// Get x and y values for point number i.
/// The function returns -1 in case of an invalid request or the point number otherwise

Int_t TGraph::GetPoint(Int_t i, Double_t &x, Double_t &y) const
{
   if (i < 0 || i >= fNpoints || !fX || !fY) return -1;
   x = fX[i];
   y = fY[i];
   return i;
}

////////////////////////////////////////////////////////////////////////////////
/// Get x value for point i.

Double_t TGraph::GetPointX(Int_t i) const
{
   if (i < 0 || i >= fNpoints || !fX)
      return -1.;

   return fX[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Get y value for point i.

Double_t TGraph::GetPointY(Int_t i) const
{
   if (i < 0 || i >= fNpoints || !fY)
      return -1.;

   return fY[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Get x axis of the graph.

TAxis *TGraph::GetXaxis() const
{
   TH1 *h = GetHistogram();
   if (!h) return 0;
   return h->GetXaxis();
}

////////////////////////////////////////////////////////////////////////////////
/// Get y axis of the graph.

TAxis *TGraph::GetYaxis() const
{
   TH1 *h = GetHistogram();
   if (!h) return 0;
   return h->GetYaxis();
}

////////////////////////////////////////////////////////////////////////////////
/// Implementation to get information on point of graph at cursor position
/// Adapted from class TH1

char *TGraph::GetObjectInfo(Int_t px, Int_t py) const
{
   // localize point
   Int_t ipoint = -2;
   Int_t i;
   // start with a small window (in case the mouse is very close to one point)
   for (i = 0; i < fNpoints; i++) {
      Int_t dpx = px - gPad->XtoAbsPixel(gPad->XtoPad(fX[i]));
      Int_t dpy = py - gPad->YtoAbsPixel(gPad->YtoPad(fY[i]));

      if (dpx * dpx + dpy * dpy < 25) {
         ipoint = i;
         break;
      }
   }

   Double_t x = gPad->PadtoX(gPad->AbsPixeltoX(px));
   Double_t y = gPad->PadtoY(gPad->AbsPixeltoY(py));

   if (ipoint == -2)
      return Form("x=%g, y=%g", x, y);

   Double_t xval = fX[ipoint];
   Double_t yval = fY[ipoint];

   return Form("x=%g, y=%g, point=%d, xval=%g, yval=%g", x, y, ipoint, xval, yval);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute Initial values of parameters for a gaussian.

void TGraph::InitGaus(Double_t xmin, Double_t xmax)
{
   Double_t allcha, sumx, sumx2, x, val, rms, mean;
   Int_t bin;
   const Double_t sqrtpi = 2.506628;

   // Compute mean value and RMS of the graph in the given range
   if (xmax <= xmin) {
      xmin = fX[0];
      xmax = fX[fNpoints-1];
   }
   Int_t np = 0;
   allcha = sumx = sumx2 = 0;
   for (bin = 0; bin < fNpoints; bin++) {
      x       = fX[bin];
      if (x < xmin || x > xmax) continue;
      np++;
      val     = fY[bin];
      sumx   += val * x;
      sumx2  += val * x * x;
      allcha += val;
   }
   if (np == 0 || allcha == 0) return;
   mean = sumx / allcha;
   rms  = TMath::Sqrt(sumx2 / allcha - mean * mean);
   Double_t binwidx = TMath::Abs((xmax - xmin) / np);
   if (rms == 0) rms = 1;
   TVirtualFitter *grFitter = TVirtualFitter::GetFitter();
   TF1 *f1 = (TF1*)grFitter->GetUserFunc();
   f1->SetParameter(0, binwidx * allcha / (sqrtpi * rms));
   f1->SetParameter(1, mean);
   f1->SetParameter(2, rms);
   f1->SetParLimits(2, 0, 10 * rms);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute Initial values of parameters for an exponential.

void TGraph::InitExpo(Double_t xmin, Double_t xmax)
{
   Double_t constant, slope;
   Int_t ifail;
   if (xmax <= xmin) {
      xmin = fX[0];
      xmax = fX[fNpoints-1];
   }
   Int_t nchanx = fNpoints;

   LeastSquareLinearFit(-nchanx, constant, slope, ifail, xmin, xmax);

   TVirtualFitter *grFitter = TVirtualFitter::GetFitter();
   TF1 *f1 = (TF1*)grFitter->GetUserFunc();
   f1->SetParameter(0, constant);
   f1->SetParameter(1, slope);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute Initial values of parameters for a polynom.

void TGraph::InitPolynom(Double_t xmin, Double_t xmax)
{
   Double_t fitpar[25];

   TVirtualFitter *grFitter = TVirtualFitter::GetFitter();
   TF1 *f1 = (TF1*)grFitter->GetUserFunc();
   Int_t npar   = f1->GetNpar();
   if (xmax <= xmin) {
      xmin = fX[0];
      xmax = fX[fNpoints-1];
   }

   LeastSquareFit(npar, fitpar, xmin, xmax);

   for (Int_t i = 0; i < npar; i++) f1->SetParameter(i, fitpar[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Insert a new point at the mouse position

Int_t TGraph::InsertPoint()
{
   Int_t px = gPad->GetEventX();
   Int_t py = gPad->GetEventY();

   //localize point where to insert
   Int_t ipoint = -2;
   Int_t i, d = 0;
   // start with a small window (in case the mouse is very close to one point)
   for (i = 0; i < fNpoints - 1; i++) {
      d = DistancetoLine(px, py, gPad->XtoPad(fX[i]), gPad->YtoPad(fY[i]), gPad->XtoPad(fX[i+1]), gPad->YtoPad(fY[i+1]));
      if (d < 5) {
         ipoint = i + 1;
         break;
      }
   }
   if (ipoint == -2) {
      //may be we are far from one point, try again with a larger window
      for (i = 0; i < fNpoints - 1; i++) {
         d = DistancetoLine(px, py, gPad->XtoPad(fX[i]), gPad->YtoPad(fY[i]), gPad->XtoPad(fX[i+1]), gPad->YtoPad(fY[i+1]));
         if (d < 10) {
            ipoint = i + 1;
            break;
         }
      }
   }
   if (ipoint == -2) {
      //distinguish between first and last point
      Int_t dpx = px - gPad->XtoAbsPixel(gPad->XtoPad(fX[0]));
      Int_t dpy = py - gPad->YtoAbsPixel(gPad->XtoPad(fY[0]));
      if (dpx * dpx + dpy * dpy < 25) ipoint = 0;
      else                      ipoint = fNpoints;
   }


   InsertPointBefore(ipoint, gPad->AbsPixeltoX(px), gPad->AbsPixeltoY(py));

   gPad->Modified();
   return ipoint;
}


////////////////////////////////////////////////////////////////////////////////
/// Insert a new point with coordinates (x,y) before the point number `ipoint`.

void TGraph::InsertPointBefore(Int_t ipoint, Double_t x, Double_t y)
{
   if (ipoint < 0) {
      Error("TGraph", "Inserted point index should be >= 0");
      return;
   }

   if (ipoint > fNpoints) {
      Error("TGraph", "Inserted point index should be <= %d", fNpoints);
      return;
   }

   if (ipoint == fNpoints) {
       SetPoint(ipoint, x, y);
       return;
   }

   Double_t **ps = ExpandAndCopy(fNpoints + 1, ipoint);
   CopyAndRelease(ps, ipoint, fNpoints++, ipoint + 1);

   // To avoid redefinitions in descendant classes
   FillZero(ipoint, ipoint + 1);

   fX[ipoint] = x;
   fY[ipoint] = y;
}


////////////////////////////////////////////////////////////////////////////////
/// Integrate the TGraph data within a given (index) range.
/// Note that this function computes the area of the polygon enclosed by the points of the TGraph.
/// The polygon segments, which are defined by the points of the TGraph, do not need to form a closed polygon,
/// since the last polygon segment, which closes the polygon, is taken as the line connecting the last TGraph point
/// with the first one. It is clear that the order of the point is essential in defining the polygon.
/// Also note that the segments should not intersect.
///
/// NB:
///  - if last=-1 (default) last is set to the last point.
///  - if (first <0) the first point (0) is taken.
///
/// ### Method:
///
/// There are many ways to calculate the surface of a polygon. It all depends on what kind of data
/// you have to deal with. The most evident solution would be to divide the polygon in triangles and
/// calculate the surface of them. But this can quickly become complicated as you will have to test
/// every segments of every triangles and check if they are intersecting with a current polygon's
/// segment or if it goes outside the polygon. Many calculations that would lead to many problems...
///
/// ### The solution (implemented by R.Brun)
/// Fortunately for us, there is a simple way to solve this problem, as long as the polygon's
/// segments don't intersect.
/// It takes the x coordinate of the current vertex and multiply it by the y coordinate of the next
/// vertex. Then it subtracts from it the result of the y coordinate of the current vertex multiplied
/// by the x coordinate of the next vertex. Then divide the result by 2 to get the surface/area.
///
/// ### Sources
///  - http://forums.wolfram.com/mathgroup/archive/1998/Mar/msg00462.html
///  - http://stackoverflow.com/questions/451426/how-do-i-calculate-the-surface-area-of-a-2d-polygon

Double_t TGraph::Integral(Int_t first, Int_t last) const
{
   if (first < 0) first = 0;
   if (last < 0) last = fNpoints - 1;
   if (last >= fNpoints) last = fNpoints - 1;
   if (first >= last) return 0;
   Int_t np = last - first + 1;
   Double_t sum = 0.0;
   //for(Int_t i=first;i<=last;i++) {
   //   Int_t j = first + (i-first+1)%np;
   //   sum += TMath::Abs(fX[i]*fY[j]);
   //   sum -= TMath::Abs(fY[i]*fX[j]);
   //}
   for (Int_t i = first; i <= last; i++) {
      Int_t j = first + (i - first + 1) % np;
      sum += (fY[i] + fY[j]) * (fX[j] - fX[i]);
   }
   return 0.5 * TMath::Abs(sum);
}

////////////////////////////////////////////////////////////////////////////////
/// Return 1 if the point (x,y) is inside the polygon defined by
/// the graph vertices 0 otherwise.
///
/// Algorithm:
///
/// The loop is executed with the end-point coordinates of a line segment
/// (X1,Y1)-(X2,Y2) and the Y-coordinate of a horizontal line.
/// The counter inter is incremented if the line (X1,Y1)-(X2,Y2) intersects
/// the horizontal line. In this case XINT is set to the X-coordinate of the
/// intersection point. If inter is an odd number, then the point x,y is within
/// the polygon.

Int_t TGraph::IsInside(Double_t x, Double_t y) const
{
   return (Int_t)TMath::IsInside(x, y, fNpoints, fX, fY);
}

////////////////////////////////////////////////////////////////////////////////
/// Least squares polynomial fitting without weights.
///
/// \param [in] m     number of parameters
/// \param [in] a     array of parameters
/// \param [in] xmin  1st point number to fit (default =0)
/// \param [in] xmax  last point number to fit (default=fNpoints-1)
///
/// based on CERNLIB routine LSQ: Translated to C++ by Rene Brun

void TGraph::LeastSquareFit(Int_t m, Double_t *a, Double_t xmin, Double_t xmax)
{
   const Double_t zero = 0.;
   const Double_t one = 1.;
   const Int_t idim = 20;

   Double_t  b[400]        /* was [20][20] */;
   Int_t i, k, l, ifail;
   Double_t power;
   Double_t da[20], xk, yk;
   Int_t n = fNpoints;
   if (xmax <= xmin) {
      xmin = fX[0];
      xmax = fX[fNpoints-1];
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
   for (k = 0; k < fNpoints; ++k) {
      xk     = fX[k];
      if (xk < xmin || xk > xmax) continue;
      np++;
      yk     = fY[k];
      power  = one;
      da[0] += yk;
      for (l = 2; l <= m; ++l) {
         power   *= xk;
         b[l-1]  += power;
         da[l-1] += power * yk;
      }
      for (l = 2; l <= m; ++l) {
         power            *= xk;
         b[m + l*20 - 21] += power;
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
      a[0] = fY[0];
      for (i = 1; i < m; ++i) a[i] = 0;
      return;
   }
   for (i = 0; i < m; ++i) a[i] = da[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Least square linear fit without weights.
///
///  Fit a straight line (a0 + a1*x) to the data in this graph.
///
/// \param [in] ndata        if ndata<0, fits the logarithm of the graph (used in InitExpo() to set
///                          the initial parameter values for a fit with exponential function.
/// \param [in] a0           constant
/// \param [in] a1           slope
/// \param [in] ifail        return parameter indicating the status of the fit (ifail=0, fit is OK)
/// \param [in] xmin, xmax   fitting range
///
///  extracted from CERNLIB LLSQ: Translated to C++ by Rene Brun

void TGraph::LeastSquareLinearFit(Int_t ndata, Double_t &a0, Double_t &a1, Int_t &ifail, Double_t xmin, Double_t xmax)
{
   Double_t xbar, ybar, x2bar;
   Int_t i;
   Double_t xybar;
   Double_t fn, xk, yk;
   Double_t det;
   if (xmax <= xmin) {
      xmin = fX[0];
      xmax = fX[fNpoints-1];
   }

   ifail = -2;
   xbar  = ybar = x2bar = xybar = 0;
   Int_t np = 0;
   for (i = 0; i < fNpoints; ++i) {
      xk = fX[i];
      if (xk < xmin || xk > xmax) continue;
      np++;
      yk = fY[i];
      if (ndata < 0) {
         if (yk <= 0) yk = 1e-9;
         yk = TMath::Log(yk);
      }
      xbar  += xk;
      ybar  += yk;
      x2bar += xk * xk;
      xybar += xk * yk;
   }
   fn    = Double_t(np);
   det   = fn * x2bar - xbar * xbar;
   ifail = -1;
   if (det <= 0) {
      if (fn > 0) a0 = ybar / fn;
      else        a0 = 0;
      a1 = 0;
      return;
   }
   ifail = 0;
   a0 = (x2bar * ybar - xbar * xybar) / det;
   a1 = (fn * xybar - xbar * ybar) / det;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this graph with its current attributes.

void TGraph::Paint(Option_t *option)
{
   TVirtualGraphPainter *painter = TVirtualGraphPainter::GetPainter();
   if (painter) painter->PaintHelper(this, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the (x,y) as a graph.

void TGraph::PaintGraph(Int_t npoints, const Double_t *x, const Double_t *y, Option_t *chopt)
{
   TVirtualGraphPainter *painter = TVirtualGraphPainter::GetPainter();
   if (painter) painter->PaintGraph(this, npoints, x, y, chopt);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the (x,y) as a histogram.

void TGraph::PaintGrapHist(Int_t npoints, const Double_t *x, const Double_t *y, Option_t *chopt)
{
   TVirtualGraphPainter *painter = TVirtualGraphPainter::GetPainter();
   if (painter) painter->PaintGrapHist(this, npoints, x, y, chopt);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the stats

void TGraph::PaintStats(TF1 *fit)
{
   TVirtualGraphPainter *painter = TVirtualGraphPainter::GetPainter();
   if (painter) painter->PaintStats(this, fit);
}

////////////////////////////////////////////////////////////////////////////////
/// Print graph values.

void TGraph::Print(Option_t *) const
{
   for (Int_t i = 0; i < fNpoints; i++) {
      printf("x[%d]=%g, y[%d]=%g\n", i, fX[i], i, fY[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively remove object from the list of functions

void TGraph::RecursiveRemove(TObject *obj)
{
   if (fFunctions) {
      if (!fFunctions->TestBit(kInvalidObject)) fFunctions->RecursiveRemove(obj);
   }
   if (fHistogram == obj) fHistogram = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete point close to the mouse position

Int_t TGraph::RemovePoint()
{
   Int_t px = gPad->GetEventX();
   Int_t py = gPad->GetEventY();

   //localize point to be deleted
   Int_t ipoint = -2;
   Int_t i;
   // start with a small window (in case the mouse is very close to one point)
   for (i = 0; i < fNpoints; i++) {
      Int_t dpx = px - gPad->XtoAbsPixel(gPad->XtoPad(fX[i]));
      Int_t dpy = py - gPad->YtoAbsPixel(gPad->YtoPad(fY[i]));
      if (dpx * dpx + dpy * dpy < 100) {
         ipoint = i;
         break;
      }
   }
   return RemovePoint(ipoint);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete point number ipoint

Int_t TGraph::RemovePoint(Int_t ipoint)
{
   if (ipoint < 0) return -1;
   if (ipoint >= fNpoints) return -1;

   Double_t **ps = ShrinkAndCopy(fNpoints - 1, ipoint);
   CopyAndRelease(ps, ipoint + 1, fNpoints--, ipoint);
   if (gPad) gPad->Modified();
   return ipoint;
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TGraph::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';
   out << "   " << std::endl;
   static Int_t frameNumber = 0;
   frameNumber++;

   if (fNpoints >= 1) {
      Int_t i;
      TString fXName = TString(GetName()) + Form("_fx%d",frameNumber);
      TString fYName = TString(GetName()) + Form("_fy%d",frameNumber);
      out << "   Double_t " << fXName << "[" << fNpoints << "] = {" << std::endl;
      for (i = 0; i < fNpoints-1; i++) out << "   " << fX[i] << "," << std::endl;
      out << "   " << fX[fNpoints-1] << "};" << std::endl;
      out << "   Double_t " << fYName << "[" << fNpoints << "] = {" << std::endl;
      for (i = 0; i < fNpoints-1; i++) out << "   " << fY[i] << "," << std::endl;
      out << "   " << fY[fNpoints-1] << "};" << std::endl;
      if (gROOT->ClassSaved(TGraph::Class())) out << "   ";
      else out << "   TGraph *";
      out << "graph = new TGraph(" << fNpoints << "," << fXName << "," << fYName << ");" << std::endl;
   } else {
      if (gROOT->ClassSaved(TGraph::Class())) out << "   ";
      else out << "   TGraph *";
      out << "graph = new TGraph();" << std::endl;
   }

   out << "   graph->SetName(" << quote << GetName() << quote << ");" << std::endl;
   out << "   graph->SetTitle(" << quote << GetTitle() << quote << ");" << std::endl;

   SaveFillAttributes(out, "graph", 0, 1001);
   SaveLineAttributes(out, "graph", 1, 1, 1);
   SaveMarkerAttributes(out, "graph", 1, 1, 1);

   if (fHistogram) {
      TString hname = fHistogram->GetName();
      hname += frameNumber;
      fHistogram->SetName(Form("Graph_%s", hname.Data()));
      fHistogram->SavePrimitive(out, "nodraw");
      out << "   graph->SetHistogram(" << fHistogram->GetName() << ");" << std::endl;
      out << "   " << std::endl;
   }

   // save list of functions
   TIter next(fFunctions);
   TObject *obj;
   while ((obj = next())) {
      obj->SavePrimitive(out, Form("nodraw #%d\n",++frameNumber));
      if (obj->InheritsFrom("TPaveStats")) {
         out << "   graph->GetListOfFunctions()->Add(ptstats);" << std::endl;
         out << "   ptstats->SetParent(graph->GetListOfFunctions());" << std::endl;
      } else {
         TString objname;
         objname.Form("%s%d",obj->GetName(),frameNumber);
         if (obj->InheritsFrom("TF1")) {
            out << "   " << objname << "->SetParent(graph);\n";
         }
         out << "   graph->GetListOfFunctions()->Add("
             << objname << ");" << std::endl;
      }
   }

   const char *l;
   l = strstr(option, "multigraph");
   if (l) {
      out << "   multigraph->Add(graph," << quote << l + 10 << quote << ");" << std::endl;
      return;
   }
   l = strstr(option, "th2poly");
   if (l) {
      out << "   " << l + 7 << "->AddBin(graph);" << std::endl;
      return;
   }
   out << "   graph->Draw(" << quote << option << quote << ");" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Set number of points in the graph
/// Existing coordinates are preserved
/// New coordinates above fNpoints are preset to 0.

void TGraph::Set(Int_t n)
{
   if (n < 0) n = 0;
   if (n == fNpoints) return;
   Double_t **ps = Allocate(n);
   CopyAndRelease(ps, 0, TMath::Min(fNpoints, n), 0);
   if (n > fNpoints) {
      FillZero(fNpoints, n, kFALSE);
   }
   fNpoints = n;
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if kNotEditable bit is not set, kFALSE otherwise.

Bool_t TGraph::GetEditable() const
{
   return TestBit(kNotEditable) ? kFALSE : kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// if editable=kFALSE, the graph cannot be modified with the mouse
///  by default a TGraph is editable

void TGraph::SetEditable(Bool_t editable)
{
   if (editable) ResetBit(kNotEditable);
   else          SetBit(kNotEditable);
}

////////////////////////////////////////////////////////////////////////////////
/// Set highlight (enable/disble) mode for the graph
/// by default highlight mode is disable

void TGraph::SetHighlight(Bool_t set)
{
   if (IsHighlight() == set) return;

   TVirtualGraphPainter *painter = TVirtualGraphPainter::GetPainter();
   if (!painter) return;
   SetBit(kIsHighlight, set);
   painter->SetHighlight(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the maximum of the graph.

void TGraph::SetMaximum(Double_t maximum)
{
   fMaximum = maximum;
   GetHistogram()->SetMaximum(maximum);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the minimum of the graph.

void TGraph::SetMinimum(Double_t minimum)
{
   fMinimum = minimum;
   GetHistogram()->SetMinimum(minimum);
}

////////////////////////////////////////////////////////////////////////////////
/// Set x and y values for point number i.

void TGraph::SetPoint(Int_t i, Double_t x, Double_t y)
{
   if (i < 0) return;
   if (fHistogram) SetBit(kResetHisto);

   if (i >= fMaxSize) {
      Double_t **ps = ExpandAndCopy(i + 1, fNpoints);
      CopyAndRelease(ps, 0, 0, 0);
   }
   if (i >= fNpoints) {
      // points above i can be not initialized
      // set zero up to i-th point to avoid redefinition
      // of this method in descendant classes
      FillZero(fNpoints, i + 1);
      fNpoints = i + 1;
   }
   fX[i] = x;
   fY[i] = y;
   if (gPad) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Set x value for point i.

void TGraph::SetPointX(Int_t i, Double_t x)
{
    SetPoint(i, x, GetPointY(i));
}

////////////////////////////////////////////////////////////////////////////////
/// Set y value for point i.

void TGraph::SetPointY(Int_t i, Double_t y)
{
    SetPoint(i, GetPointX(i), y);
}

////////////////////////////////////////////////////////////////////////////////
/// Set graph name.
void TGraph::SetName(const char *name)
{
   fName = name;
   if (fHistogram) fHistogram->SetName(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Change (i.e. set) the title
///
/// if title is in the form `stringt;stringx;stringy;stringz`
/// the graph title is set to `stringt`, the x axis title to `stringx`,
/// the y axis title to `stringy`, and the z axis title to `stringz`.
///
/// To insert the character `;` in one of the titles, one should use `#;`
/// or `#semicolon`.

void TGraph::SetTitle(const char* title)
{
   fTitle = title;
   fTitle.ReplaceAll("#;",2,"#semicolon",10);
   Int_t p = fTitle.Index(";");

   if (p>0) {
      if (!fHistogram) GetHistogram();
      fHistogram->SetTitle(title);
      Int_t n = fTitle.Length()-p;
      if (p>0) fTitle.Remove(p,n);
      fTitle.ReplaceAll("#semicolon",10,"#;",2);
   } else {
      if (fHistogram) fHistogram->SetTitle(title);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set graph name and title

void TGraph::SetNameTitle(const char *name, const char *title)
{
   SetName(name);
   SetTitle(title);
}

////////////////////////////////////////////////////////////////////////////////
/// if size*2 <= fMaxSize allocate new arrays of size points,
/// copy points [0,oend).
/// Return newarray (passed or new instance if it was zero
/// and allocations are needed)

Double_t **TGraph::ShrinkAndCopy(Int_t size, Int_t oend)
{
   if (size * 2 > fMaxSize || !fMaxSize) {
      return 0;
   }
   Double_t **newarrays = Allocate(size);
   CopyPoints(newarrays, 0, oend, 0);
   return newarrays;
}

////////////////////////////////////////////////////////////////////////////////
/// Sorts the points of this TGraph using in-place quicksort (see e.g. older glibc).
/// To compare two points the function parameter greaterfunc is used (see TGraph::CompareX for an
/// example of such a method, which is also the default comparison function for Sort). After
/// the sort, greaterfunc(this, i, j) will return kTRUE for all i>j if ascending == kTRUE, and
/// kFALSE otherwise.
///
/// The last two parameters are used for the recursive quick sort, stating the range to be sorted
///
/// Examples:
/// ~~~ {.cpp}
///   // sort points along x axis
///   graph->Sort();
///   // sort points along their distance to origin
///   graph->Sort(&TGraph::CompareRadius);
///
///   Bool_t CompareErrors(const TGraph* gr, Int_t i, Int_t j) {
///     const TGraphErrors* ge=(const TGraphErrors*)gr;
///     return (ge->GetEY()[i]>ge->GetEY()[j]); }
///   // sort using the above comparison function, largest errors first
///   graph->Sort(&CompareErrors, kFALSE);
/// ~~~

void TGraph::Sort(Bool_t (*greaterfunc)(const TGraph*, Int_t, Int_t) /*=TGraph::CompareX()*/,
                  Bool_t ascending /*=kTRUE*/, Int_t low /* =0 */, Int_t high /* =-1111 */)
{

   // set the bit in case of an ascending =sort in X
   if (greaterfunc == TGraph::CompareX && ascending  && low == 0 && high == -1111)
      SetBit(TGraph::kIsSortedX);

   if (high == -1111) high = GetN() - 1;
   //  Termination condition
   if (high <= low) return;

   int left, right;
   left = low; // low is the pivot element
   right = high;
   while (left < right) {
      // move left while item < pivot
      while (left <= high && greaterfunc(this, left, low) != ascending)
         left++;
      // move right while item > pivot
      while (right > low && greaterfunc(this, right, low) == ascending)
         right--;
      if (left < right && left < high && right > low)
         SwapPoints(left, right);
   }
   // right is final position for the pivot
   if (right > low)
      SwapPoints(low, right);
   Sort(greaterfunc, ascending, low, right - 1);
   Sort(greaterfunc, ascending, right + 1, high);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TGraph.

void TGraph::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         b.ReadClassBuffer(TGraph::Class(), this, R__v, R__s, R__c);
         if (fHistogram) fHistogram->SetDirectory(0);
         TIter next(fFunctions);
         TObject *obj;
         while ((obj = next())) {
            if (obj->InheritsFrom(TF1::Class())) {
               TF1 *f1 = (TF1*)obj;
               f1->SetParent(this);
            }
         }
         fMaxSize = fNpoints;
         return;
      }
      //====process old versions before automatic schema evolution
      TNamed::Streamer(b);
      TAttLine::Streamer(b);
      TAttFill::Streamer(b);
      TAttMarker::Streamer(b);
      b >> fNpoints;
      fMaxSize = fNpoints;
      fX = new Double_t[fNpoints];
      fY = new Double_t[fNpoints];
      if (R__v < 2) {
         Float_t *x = new Float_t[fNpoints];
         Float_t *y = new Float_t[fNpoints];
         b.ReadFastArray(x, fNpoints);
         b.ReadFastArray(y, fNpoints);
         for (Int_t i = 0; i < fNpoints; i++) {
            fX[i] = x[i];
            fY[i] = y[i];
         }
         delete [] y;
         delete [] x;
      } else {
         b.ReadFastArray(fX, fNpoints);
         b.ReadFastArray(fY, fNpoints);
      }
      b >> fFunctions;
      b >> fHistogram;
      if (fHistogram) fHistogram->SetDirectory(0);
      if (R__v < 2) {
         Float_t mi, ma;
         b >> mi;
         b >> ma;
         fMinimum = mi;
         fMaximum = ma;
      } else {
         b >> fMinimum;
         b >> fMaximum;
      }
      b.CheckByteCount(R__s, R__c, TGraph::IsA());
      //====end of old versions

   } else {
      b.WriteClassBuffer(TGraph::Class(), this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Swap points.

void TGraph::SwapPoints(Int_t pos1, Int_t pos2)
{
   SwapValues(fX, pos1, pos2);
   SwapValues(fY, pos1, pos2);
}

////////////////////////////////////////////////////////////////////////////////
/// Swap values.

void TGraph::SwapValues(Double_t* arr, Int_t pos1, Int_t pos2)
{
   Double_t tmp = arr[pos1];
   arr[pos1] = arr[pos2];
   arr[pos2] = tmp;
}

////////////////////////////////////////////////////////////////////////////////
/// Set current style settings in this graph
/// This function is called when either TCanvas::UseCurrentStyle
/// or TROOT::ForceStyle have been invoked.

void TGraph::UseCurrentStyle()
{
   if (gStyle->IsReading()) {
      SetFillColor(gStyle->GetHistFillColor());
      SetFillStyle(gStyle->GetHistFillStyle());
      SetLineColor(gStyle->GetHistLineColor());
      SetLineStyle(gStyle->GetHistLineStyle());
      SetLineWidth(gStyle->GetHistLineWidth());
      SetMarkerColor(gStyle->GetMarkerColor());
      SetMarkerStyle(gStyle->GetMarkerStyle());
      SetMarkerSize(gStyle->GetMarkerSize());
   } else {
      gStyle->SetHistFillColor(GetFillColor());
      gStyle->SetHistFillStyle(GetFillStyle());
      gStyle->SetHistLineColor(GetLineColor());
      gStyle->SetHistLineStyle(GetLineStyle());
      gStyle->SetHistLineWidth(GetLineWidth());
      gStyle->SetMarkerColor(GetMarkerColor());
      gStyle->SetMarkerStyle(GetMarkerStyle());
      gStyle->SetMarkerSize(GetMarkerSize());
   }
   if (fHistogram) fHistogram->UseCurrentStyle();

   TIter next(GetListOfFunctions());
   TObject *obj;

   while ((obj = next())) {
      obj->UseCurrentStyle();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Adds all graphs from the collection to this graph.
/// Returns the total number of poins in the result or -1 in case of an error.

Int_t TGraph::Merge(TCollection* li)
{
   TIter next(li);
   while (TObject* o = next()) {
      TGraph *g = dynamic_cast<TGraph*>(o);
      if (!g) {
         Error("Merge",
               "Cannot merge - an object which doesn't inherit from TGraph found in the list");
         return -1;
      }
      DoMerge(g);
   }
   return GetN();
}

////////////////////////////////////////////////////////////////////////////////
///  protected function to perform the merge operation of a graph

Bool_t TGraph::DoMerge(const TGraph* g)
{
   Double_t x = 0, y = 0;
   for (Int_t i = 0 ; i < g->GetN(); i++) {
      g->GetPoint(i, x, y);
      SetPoint(GetN(), x, y);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Move all graph points on specified values dx,dy
/// If log argument specified, calculation done in logarithmic scale like:
///  new_value = exp( log(old_value) + delta );

void TGraph::MovePoints(Double_t dx, Double_t dy, Bool_t logx, Bool_t logy)
{
   Double_t x = 0, y = 0;
   for (Int_t i = 0 ; i < GetN(); i++) {
      GetPoint(i, x, y);
      if (!logx) {
         x += dx;
      } else if (x > 0) {
         x = TMath::Exp(TMath::Log(x) + dx);
      }
      if (!logy) {
         y += dy;
      } else if (y > 0) {
         y = TMath::Exp(TMath::Log(y) + dy);
      }
      SetPoint(i, x, y);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Find zero of a continuous function.
/// This function finds a real zero of the continuous real
/// function Y(X) in a given interval (A,B). See accompanying
/// notes for details of the argument list and calling sequence

void TGraph::Zero(Int_t &k, Double_t AZ, Double_t BZ, Double_t E2, Double_t &X, Double_t &Y
                  , Int_t maxiterations)
{
   static Double_t a, b, ya, ytest, y1, x1, h;
   static Int_t j1, it, j3, j2;
   Double_t yb, x2;
   yb = 0;

   //       Calculate Y(X) at X=AZ.
   if (k <= 0) {
      a  = AZ;
      b  = BZ;
      X  = a;
      j1 = 1;
      it = 1;
      k  = j1;
      return;
   }

   //       Test whether Y(X) is sufficiently small.

   if (TMath::Abs(Y) <= E2) {
      k = 2;
      return;
   }

   //       Calculate Y(X) at X=BZ.

   if (j1 == 1) {
      ya = Y;
      X  = b;
      j1 = 2;
      return;
   }
   //       Test whether the signs of Y(AZ) and Y(BZ) are different.
   //       if not, begin the binary subdivision.

   if (j1 != 2) goto L100;
   if (ya * Y < 0) goto L120;
   x1 = a;
   y1 = ya;
   j1 = 3;
   h  = b - a;
   j2 = 1;
   x2 = a + 0.5 * h;
   j3 = 1;
   it++;      //*-*-   Check whether (maxiterations) function values have been calculated.
   if (it >= maxiterations) k = j1;
   else                     X = x2;
   return;

   //      Test whether a bracket has been found .
   //      If not,continue the search

L100:
   if (j1 > 3) goto L170;
   if (ya*Y >= 0) {
      if (j3 >= j2) {
         h  = 0.5 * h;
         j2 = 2 * j2;
         a  = x1;
         ya = y1;
         x2 = a + 0.5 * h;
         j3 = 1;
      } else {
         a  = X;
         ya = Y;
         x2 = X + h;
         j3++;
      }
      it++;
      if (it >= maxiterations) k = j1;
      else                     X = x2;
      return;
   }

   //       The first bracket has been found.calculate the next X by the
   //       secant method based on the bracket.

L120:
   b  = X;
   yb = Y;
   j1 = 4;
L130:
   if (TMath::Abs(ya) > TMath::Abs(yb)) {
      x1 = a;
      y1 = ya;
      X  = b;
      Y  = yb;
   } else                                 {
      x1 = b;
      y1 = yb;
      X  = a;
      Y  = ya;
   }

   //       Use the secant method based on the function values y1 and Y.
   //       check that x2 is inside the interval (a,b).

L150:
   x2    = X - Y * (X - x1) / (Y - y1);
   x1    = X;
   y1    = Y;
   ytest = 0.5 * TMath::Min(TMath::Abs(ya), TMath::Abs(yb));
   if ((x2 - a)*(x2 - b) < 0) {
      it++;
      if (it >= maxiterations) k = j1;
      else                     X = x2;
      return;
   }

   //       Calculate the next value of X by bisection . Check whether
   //       the maximum accuracy has been achieved.

L160:
   x2    = 0.5 * (a + b);
   ytest = 0;
   if ((x2 - a)*(x2 - b) >= 0) {
      k = 2;
      return;
   }
   it++;
   if (it >= maxiterations) k = j1;
   else                     X = x2;
   return;


   //       Revise the bracket (a,b).

L170:
   if (j1 != 4) return;
   if (ya * Y < 0) {
      b  = X;
      yb = Y;
   } else          {
      a  = X;
      ya = Y;
   }

   //       Use ytest to decide the method for the next value of X.

   if (ytest <= 0) goto L130;
   if (TMath::Abs(Y) - ytest <= 0) goto L150;
   goto L160;
}
