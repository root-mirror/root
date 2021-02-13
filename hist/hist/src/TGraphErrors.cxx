// @(#)root/hist:$Id$
// Author: Rene Brun   15/09/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#include "TROOT.h"
#include "TBuffer.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVirtualPad.h"
#include "TH1.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TSystem.h"
#include "strtok.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>

ClassImp(TGraphErrors);


////////////////////////////////////////////////////////////////////////////////

/** \class TGraphErrors
    \ingroup Hist
A TGraphErrors is a TGraph with error bars.

The TGraphErrors painting is performed thanks to the TGraphPainter
class. All details about the various painting options are given in this class.

The picture below gives an example:

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
   c1->SetFillColor(42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   const Int_t n = 10;
   Double_t x[n]  = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
   Double_t y[n]  = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
   Double_t ex[n] = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};
   Double_t ey[n] = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
   auto gr = new TGraphErrors(n,x,y,ex,ey);
   gr->SetTitle("TGraphErrors Example");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("ALP");
}
End_Macro
*/


////////////////////////////////////////////////////////////////////////////////
/// TGraphErrors default constructor.

TGraphErrors::TGraphErrors(): TGraph()
{
   if (!CtorAllocate()) return;
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphErrors normal constructor.
///
///  the arrays are preset to zero

TGraphErrors::TGraphErrors(Int_t n)
   : TGraph(n)
{
   if (!CtorAllocate()) return;
   FillZero(0, fNpoints);
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphErrors normal constructor.
///
///  if ex or ey are null, the corresponding arrays are preset to zero

TGraphErrors::TGraphErrors(Int_t n, const Float_t *x, const Float_t *y, const Float_t *ex, const Float_t *ey)
   : TGraph(n, x, y)
{
   if (!CtorAllocate()) return;

   for (Int_t i = 0; i < n; i++) {
      if (ex) fEX[i] = ex[i];
      else    fEX[i] = 0;
      if (ey) fEY[i] = ey[i];
      else    fEY[i] = 0;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphErrors normal constructor.
///
///  if ex or ey are null, the corresponding arrays are preset to zero

TGraphErrors::TGraphErrors(Int_t n, const Double_t *x, const Double_t *y, const Double_t *ex, const Double_t *ey)
   : TGraph(n, x, y)
{
   if (!CtorAllocate()) return;

   n = sizeof(Double_t) * fNpoints;
   if (ex) memcpy(fEX, ex, n);
   else    memset(fEX, 0, n);
   if (ey) memcpy(fEY, ey, n);
   else    memset(fEY, 0, n);
}


////////////////////////////////////////////////////////////////////////////////
/// constructor with four vectors of floats in input
/// A grapherrors is built with the X coordinates taken from vx and Y coord from vy
/// and the errors from vectors vex and vey.
/// The number of points in the graph is the minimum of number of points
/// in vx and vy.

TGraphErrors::TGraphErrors(const TVectorF &vx, const TVectorF &vy, const TVectorF &vex, const TVectorF &vey)
   : TGraph(TMath::Min(vx.GetNrows(), vy.GetNrows()), vx.GetMatrixArray(), vy.GetMatrixArray() )
{
   if (!CtorAllocate()) return;
   Int_t ivexlow = vex.GetLwb();
   Int_t iveylow = vey.GetLwb();
   for (Int_t i = 0; i < fNpoints; i++) {
      fEX[i]  = vex(i + ivexlow);
      fEY[i]  = vey(i + iveylow);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// constructor with four vectors of doubles in input
/// A grapherrors is built with the X coordinates taken from vx and Y coord from vy
/// and the errors from vectors vex and vey.
/// The number of points in the graph is the minimum of number of points
/// in vx and vy.

TGraphErrors::TGraphErrors(const TVectorD  &vx, const TVectorD  &vy, const TVectorD  &vex, const TVectorD  &vey)
   : TGraph(TMath::Min(vx.GetNrows(), vy.GetNrows()), vx.GetMatrixArray(), vy.GetMatrixArray() )
{
   if (!CtorAllocate()) return;
   Int_t ivexlow = vex.GetLwb();
   Int_t iveylow = vey.GetLwb();
   for (Int_t i = 0; i < fNpoints; i++) {
      fEX[i]  = vex(i + ivexlow);
      fEY[i]  = vey(i + iveylow);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphErrors copy constructor

TGraphErrors::TGraphErrors(const TGraphErrors &gr)
   : TGraph(gr)
{
   if (!CtorAllocate()) return;

   Int_t n = sizeof(Double_t) * fNpoints;
   memcpy(fEX, gr.fEX, n);
   memcpy(fEY, gr.fEY, n);
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphErrors assignment operator

TGraphErrors& TGraphErrors::operator=(const TGraphErrors &gr)
{
   if (this != &gr) {
      TGraph::operator=(gr);
      // N.B CtorAllocate does not delete arrays
      if (fEX) delete [] fEX;
      if (fEY) delete [] fEY;
      if (!CtorAllocate()) return *this;

      Int_t n = sizeof(Double_t) * fNpoints;
      memcpy(fEX, gr.fEX, n);
      memcpy(fEY, gr.fEY, n);
   }
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphErrors constructor importing its parameters from the TH1 object passed as argument

TGraphErrors::TGraphErrors(const TH1 *h)
   : TGraph(h)
{
   if (!CtorAllocate()) return;

   for (Int_t i = 0; i < fNpoints; i++) {
      fEX[i] = h->GetBinWidth(i + 1) * gStyle->GetErrorX();
      fEY[i] = h->GetBinError(i + 1);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// GraphErrors constructor reading input from filename
/// filename is assumed to contain at least 2 columns of numbers
///
/// Convention for format (default="%lg %lg %lg %lg)
///   - format = "%lg %lg"         read only 2 first columns into X,Y
///   - format = "%lg %lg %lg"     read only 3 first columns into X,Y and EY
///   - format = "%lg %lg %lg %lg" read only 4 first columns into X,Y,EX,EY.
///
/// For files separated by a specific delimiter different from ' ' and '\t' (e.g. ';' in csv files)
/// you can avoid using %*s to bypass this delimiter by explicitly specify the "option" argument,
/// e.g. option=" \t,;" for columns of figures separated by any of these characters (' ', '\t', ',', ';')
/// used once (e.g. "1;1") or in a combined way (" 1;,;;  1").
/// Note in that case, the instantiation is about 2 times slower.
/// In case a delimiter is specified, the format "%lg %lg %lg" will read X,Y,EX.

TGraphErrors::TGraphErrors(const char *filename, const char *format, Option_t *option)
   : TGraph(100)
{
   if (!CtorAllocate()) return;
   Double_t x, y, ex, ey;
   TString fname = filename;
   gSystem->ExpandPathName(fname);
   std::ifstream infile(fname.Data());
   if (!infile.good()) {
      MakeZombie();
      Error("TGraphErrors", "Cannot open file: %s, TGraphErrors is Zombie", filename);
      fNpoints = 0;
      return;
   }
   std::string line;
   Int_t np = 0;

   if (strcmp(option, "") == 0) { // No delimiters specified (standard constructor).

      Int_t ncol = CalculateScanfFields(format);  //count number of columns in format
      Int_t res;
      while (std::getline(infile, line, '\n')) {
         ex = ey = 0;
         if (ncol < 3) {
            res = sscanf(line.c_str(), format, &x, &y);
         } else if (ncol < 4) {
            res = sscanf(line.c_str(), format, &x, &y, &ey);
         } else {
            res = sscanf(line.c_str(), format, &x, &y, &ex, &ey);
         }
         if (res < 2) {
            continue; //skip empty and ill-formed lines
         }
         SetPoint(np, x, y);
         SetPointError(np, ex, ey);
         np++;
      }
      Set(np);

   } else { // A delimiter has been specified in "option"

      // Checking format and creating its boolean equivalent
      TString format_ = TString(format) ;
      format_.ReplaceAll(" ", "") ;
      format_.ReplaceAll("\t", "") ;
      format_.ReplaceAll("lg", "") ;
      format_.ReplaceAll("s", "") ;
      format_.ReplaceAll("%*", "0") ;
      format_.ReplaceAll("%", "1") ;
      if (!format_.IsDigit()) {
         Error("TGraphErrors", "Incorrect input format! Allowed format tags are {\"%%lg\",\"%%*lg\" or \"%%*s\"}");
         return ;
      }
      Int_t ntokens = format_.Length() ;
      if (ntokens < 2) {
         Error("TGraphErrors", "Incorrect input format! Only %d tag(s) in format whereas at least 2 \"%%lg\" tags are expected!", ntokens);
         return ;
      }
      Int_t ntokensToBeSaved = 0 ;
      Bool_t * isTokenToBeSaved = new Bool_t [ntokens] ;
      for (Int_t idx = 0; idx < ntokens; idx++) {
         isTokenToBeSaved[idx] = TString::Format("%c", format_[idx]).Atoi() ; //atoi(&format_[idx]) does not work for some reason...
         if (isTokenToBeSaved[idx] == 1) {
            ntokensToBeSaved++ ;
         }
      }
      if (ntokens >= 2 && (ntokensToBeSaved < 2 || ntokensToBeSaved > 4)) { //first condition not to repeat the previous error message
         Error("TGraphErrors", "Incorrect input format! There are %d \"%%lg\" tag(s) in format whereas 2,3 or 4 are expected!", ntokensToBeSaved);
         delete [] isTokenToBeSaved ;
         return ;
      }

      // Initializing loop variables
      Bool_t isLineToBeSkipped = kFALSE ; //empty and ill-formed lines
      char * token = NULL ;
      TString token_str = "" ;
      Int_t token_idx = 0 ;
      Double_t * value = new Double_t [4] ; //x,y,ex,ey buffers
      for (Int_t k = 0; k < 4; k++) {
         value[k] = 0. ;
      }
      Int_t value_idx = 0 ;

      // Looping
      char *rest;
      while (std::getline(infile, line, '\n')) {
         if (line != "") {
            if (line[line.size() - 1] == char(13)) {  // removing DOS CR character
               line.erase(line.end() - 1, line.end()) ;
            }
            token = R__STRTOK_R(const_cast<char *>(line.c_str()), option, &rest);
            while (token != NULL && value_idx < ntokensToBeSaved) {
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
            if (!isLineToBeSkipped && value_idx > 1) { //i.e. 2,3 or 4
               x = value[0] ;
               y = value[1] ;
               ex = value[2] ;
               ey = value[3] ;
               SetPoint(np, x, y) ;
               SetPointError(np, ex, ey);
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
/// TGraphErrors default destructor.

TGraphErrors::~TGraphErrors()
{
   delete [] fEX;
   delete [] fEY;
}


////////////////////////////////////////////////////////////////////////////////
/// apply function to all the data points
/// y = f(x,y)
///
/// The error is calculated as ey=(f(x,y+ey)-f(x,y-ey))/2
/// This is the same as error(fy) = df/dy * ey for small errors
///
/// For generic functions the symmetric errors might become non-symmetric
/// and are averaged here. Use TGraphAsymmErrors if desired.
///
/// error on x doesn't change
/// function suggested/implemented by Miroslav Helbich <helbich@mail.desy.de>

void TGraphErrors::Apply(TF1 *f)
{
   Double_t x, y, ex, ey;

   if (fHistogram) {
      delete fHistogram;
      fHistogram = 0;
   }
   for (Int_t i = 0; i < GetN(); i++) {
      GetPoint(i, x, y);
      ex = GetErrorX(i);
      ey = GetErrorY(i);

      SetPoint(i, x, f->Eval(x, y));
      SetPointError(i, ex, TMath::Abs(f->Eval(x, y + ey) - f->Eval(x, y - ey)) / 2.);
   }
   if (gPad) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// apply function to all the data points
/// x = f(x,y)
///
/// The error is calculated as ex=(f(x+ex,y)-f(x-ex,y))/2
/// This is the same as error(fx) = df/dx * ex for small errors
///
/// For generic functions the symmetric errors might become non-symmetric
/// and are averaged here. Use TGraphAsymmErrors if desired.
///
/// error on y doesn't change

void TGraphErrors::ApplyX(TF1 *f)
{
   Double_t x, y, ex, ey;

   if (fHistogram) {
      delete fHistogram;
      fHistogram = 0;
   }
   for (Int_t i = 0; i < GetN(); i++) {
      GetPoint(i, x, y);
      ex = GetErrorX(i);
      ey = GetErrorY(i);

      SetPoint(i, f->Eval(x,y), y);
      SetPointError(i, TMath::Abs(f->Eval(x + ex, y) - f->Eval(x - ex, y)) / 2. , ey);
   }
   if (gPad) gPad->Modified();
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate scan fields.

Int_t TGraphErrors::CalculateScanfFields(const char *fmt)
{
   Int_t fields = 0;
   while ((fmt = strchr(fmt, '%'))) {
      Bool_t skip = kFALSE;
      while (*(++fmt)) {
         if ('[' == *fmt) {
            if (*++fmt && '^' == *fmt) ++fmt; // "%[^]a]"
            if (*++fmt && ']' == *fmt) ++fmt; // "%[]a]" or "%[^]a]"
            while (*fmt && *fmt != ']')
               ++fmt;
            if (!skip) ++fields;
            break;
         }
         if ('%' == *fmt) break; // %% literal %
         if ('*' == *fmt) {
            skip = kTRUE; // %*d -- skip a number
         } else if (strchr("dDiouxXxfegEscpn", *fmt)) {
            if (!skip) ++fields;
            break;
         }
         // skip modifiers & field width
      }
   }
   return fields;
}


////////////////////////////////////////////////////////////////////////////////
/// Compute range.

void TGraphErrors::ComputeRange(Double_t &xmin, Double_t &ymin, Double_t &xmax, Double_t &ymax) const
{
   TGraph::ComputeRange(xmin, ymin, xmax, ymax);

   for (Int_t i = 0; i < fNpoints; i++) {
      if (fX[i] - fEX[i] < xmin) {
         if (gPad && gPad->GetLogx()) {
            if (fEX[i] < fX[i]) xmin = fX[i] - fEX[i];
            else                xmin = TMath::Min(xmin, fX[i] / 3);
         } else {
            xmin = fX[i] - fEX[i];
         }
      }
      if (fX[i] + fEX[i] > xmax) xmax = fX[i] + fEX[i];
      if (fY[i] - fEY[i] < ymin) {
         if (gPad && gPad->GetLogy()) {
            if (fEY[i] < fY[i]) ymin = fY[i] - fEY[i];
            else                ymin = TMath::Min(ymin, fY[i] / 3);
         } else {
            ymin = fY[i] - fEY[i];
         }
      }
      if (fY[i] + fEY[i] > ymax) ymax = fY[i] + fEY[i];
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Copy and release.

void TGraphErrors::CopyAndRelease(Double_t **newarrays,
                                  Int_t ibegin, Int_t iend, Int_t obegin)
{
   CopyPoints(newarrays, ibegin, iend, obegin);
   if (newarrays) {
      delete[] fX;
      fX = newarrays[2];
      delete[] fY;
      fY = newarrays[3];
      delete[] fEX;
      fEX = newarrays[0];
      delete[] fEY;
      fEY = newarrays[1];
      delete[] newarrays;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Copy errors from fEX and fEY to arrays[0] and arrays[1]
/// or to fX and fY. Copy points.

Bool_t TGraphErrors::CopyPoints(Double_t **arrays, Int_t ibegin, Int_t iend,
                                Int_t obegin)
{
   if (TGraph::CopyPoints(arrays ? arrays + 2 : 0, ibegin, iend, obegin)) {
      Int_t n = (iend - ibegin) * sizeof(Double_t);
      if (arrays) {
         memmove(&arrays[0][obegin], &fEX[ibegin], n);
         memmove(&arrays[1][obegin], &fEY[ibegin], n);
      } else {
         memmove(&fEX[obegin], &fEX[ibegin], n);
         memmove(&fEY[obegin], &fEY[ibegin], n);
      }
      return kTRUE;
   } else {
      return kFALSE;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor allocate.
///Note: This function should be called only from the constructor
/// since it does not delete previously existing arrays

Bool_t TGraphErrors::CtorAllocate()
{

   if (!fNpoints) {
      fEX = fEY = 0;
      return kFALSE;
   } else {
      fEX = new Double_t[fMaxSize];
      fEY = new Double_t[fMaxSize];
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///  protected function to perform the merge operation of a graph with errors

Bool_t TGraphErrors::DoMerge(const TGraph *g)
{
   if (g->GetN() == 0) return kFALSE;

   Double_t * ex = g->GetEX();
   Double_t * ey = g->GetEY();
   if (ex == 0 || ey == 0 ) {
      if (g->IsA() != TGraph::Class() )
         Warning("DoMerge","Merging a %s is not compatible with a TGraphErrors - errors will be ignored",g->IsA()->GetName());
      return TGraph::DoMerge(g);
   }
   for (Int_t i = 0 ; i < g->GetN(); i++) {
      Int_t ipoint = GetN();
      Double_t x = g->GetX()[i];
      Double_t y = g->GetY()[i];
      SetPoint(ipoint, x, y);
      SetPointError( ipoint, ex[i], ey[i] );
   }
   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Set zero values for point arrays in the range [begin, end]

void TGraphErrors::FillZero(Int_t begin, Int_t end, Bool_t from_ctor)
{
   if (!from_ctor) {
      TGraph::FillZero(begin, end, from_ctor);
   }
   Int_t n = (end - begin) * sizeof(Double_t);
   memset(fEX + begin, 0, n);
   memset(fEY + begin, 0, n);
}


////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It returns the error along X at point i.

Double_t TGraphErrors::GetErrorX(Int_t i) const
{
   if (i < 0 || i >= fNpoints) return -1;
   if (fEX) return fEX[i];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It returns the error along Y at point i.

Double_t TGraphErrors::GetErrorY(Int_t i) const
{
   if (i < 0 || i >= fNpoints) return -1;
   if (fEY) return fEY[i];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It returns the error along X at point i.

Double_t TGraphErrors::GetErrorXhigh(Int_t i) const
{
   if (i < 0 || i >= fNpoints) return -1;
   if (fEX) return fEX[i];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It returns the error along X at point i.

Double_t TGraphErrors::GetErrorXlow(Int_t i) const
{
   if (i < 0 || i >= fNpoints) return -1;
   if (fEX) return fEX[i];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It returns the error along X at point i.

Double_t TGraphErrors::GetErrorYhigh(Int_t i) const
{
   if (i < 0 || i >= fNpoints) return -1;
   if (fEY) return fEY[i];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It returns the error along X at point i.

Double_t TGraphErrors::GetErrorYlow(Int_t i) const
{
   if (i < 0 || i >= fNpoints) return -1;
   if (fEY) return fEY[i];
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Adds all graphs with errors from the collection to this graph.
/// Returns the total number of poins in the result or -1 in case of an error.

Int_t TGraphErrors::Merge(TCollection* li)
{
   TIter next(li);
   while (TObject* o = next()) {
      TGraph *g = dynamic_cast<TGraph*>(o);
      if (!g) {
         Error("Merge",
               "Cannot merge - an object which doesn't inherit from TGraph found in the list");
         return -1;
      }
      int n0 = GetN();
      int n1 = n0+g->GetN();
      Set(n1);
      Double_t * x = g->GetX();
      Double_t * y = g->GetY();
      Double_t * ex = g->GetEX();
      Double_t * ey = g->GetEY();
      for (Int_t i = 0 ; i < g->GetN(); i++) {
         SetPoint(n0+i, x[i], y[i]);
         if (ex) fEX[n0+i] = ex[i];
         if (ey) fEY[n0+i] = ey[i];
      }
   }
   return GetN();
}

////////////////////////////////////////////////////////////////////////////////
/// Print graph and errors values.

void TGraphErrors::Print(Option_t *) const
{
   for (Int_t i = 0; i < fNpoints; i++) {
      printf("x[%d]=%g, y[%d]=%g, ex[%d]=%g, ey[%d]=%g\n", i, fX[i], i, fY[i], i, fEX[i], i, fEY[i]);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TGraphErrors::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';
   out << "   " << std::endl;
   static Int_t frameNumber = 1000;
   frameNumber++;

   Int_t i;
   TString fXName  = TString(GetName()) + Form("_fx%d",frameNumber);
   TString fYName  = TString(GetName()) + Form("_fy%d",frameNumber);
   TString fEXName = TString(GetName()) + Form("_fex%d",frameNumber);
   TString fEYName = TString(GetName()) + Form("_fey%d",frameNumber);
   out << "   Double_t " << fXName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fX[i] << "," << std::endl;
   out << "   " << fX[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fYName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fY[i] << "," << std::endl;
   out << "   " << fY[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fEXName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEX[i] << "," << std::endl;
   out << "   " << fEX[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fEYName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEY[i] << "," << std::endl;
   out << "   " << fEY[fNpoints-1] << "};" << std::endl;

   if (gROOT->ClassSaved(TGraphErrors::Class())) out << "   ";
   else out << "   TGraphErrors *";
   out << "gre = new TGraphErrors(" << fNpoints << ","
                                    << fXName   << ","  << fYName  << ","
                                    << fEXName  << ","  << fEYName << ");"
                                    << std::endl;

   out << "   gre->SetName(" << quote << GetName() << quote << ");" << std::endl;
   out << "   gre->SetTitle(" << quote << GetTitle() << quote << ");" << std::endl;

   SaveFillAttributes(out, "gre", 0, 1001);
   SaveLineAttributes(out, "gre", 1, 1, 1);
   SaveMarkerAttributes(out, "gre", 1, 1, 1);

   if (fHistogram) {
      TString hname = fHistogram->GetName();
      hname += frameNumber;
      fHistogram->SetName(Form("Graph_%s", hname.Data()));
      fHistogram->SavePrimitive(out, "nodraw");
      out << "   gre->SetHistogram(" << fHistogram->GetName() << ");" << std::endl;
      out << "   " << std::endl;
   }

   // save list of functions
   TIter next(fFunctions);
   TObject *obj;
   while ((obj = next())) {
      obj->SavePrimitive(out, Form("nodraw #%d\n",++frameNumber));
      if (obj->InheritsFrom("TPaveStats")) {
         out << "   gre->GetListOfFunctions()->Add(ptstats);" << std::endl;
         out << "   ptstats->SetParent(gre->GetListOfFunctions());" << std::endl;
      } else {
         TString objname;
         objname.Form("%s%d",obj->GetName(),frameNumber);
         if (obj->InheritsFrom("TF1")) {
            out << "   " << objname << "->SetParent(gre);\n";
         }
         out << "   gre->GetListOfFunctions()->Add("
             << objname << ");" << std::endl;
      }
   }

   const char *l = strstr(option, "multigraph");
   if (l) {
      out << "   multigraph->Add(gre," << quote << l + 10 << quote << ");" << std::endl;
   } else {
      out << "   gre->Draw(" << quote << option << quote << ");" << std::endl;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Set ex and ey values for point pointed by the mouse.

void TGraphErrors::SetPointError(Double_t ex, Double_t ey)
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
      if (dpx * dpx + dpy * dpy < 25) {
         ipoint = i;
         break;
      }
   }
   if (ipoint == -2) return;

   fEX[ipoint] = ex;
   fEY[ipoint] = ey;
   gPad->Modified();
}


////////////////////////////////////////////////////////////////////////////////
/// Set ex and ey values for point number i.

void TGraphErrors::SetPointError(Int_t i, Double_t ex, Double_t ey)
{
   if (i < 0) return;
   if (i >= fNpoints) {
      // re-allocate the object
      TGraphErrors::SetPoint(i, 0, 0);
   }
   fEX[i] = ex;
   fEY[i] = ey;
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TGraphErrors.

void TGraphErrors::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         b.ReadClassBuffer(TGraphErrors::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TGraph::Streamer(b);
      fEX = new Double_t[fNpoints];
      fEY = new Double_t[fNpoints];
      if (R__v < 2) {
         Float_t *ex = new Float_t[fNpoints];
         Float_t *ey = new Float_t[fNpoints];
         b.ReadFastArray(ex, fNpoints);
         b.ReadFastArray(ey, fNpoints);
         for (Int_t i = 0; i < fNpoints; i++) {
            fEX[i] = ex[i];
            fEY[i] = ey[i];
         }
         delete [] ey;
         delete [] ex;
      } else {
         b.ReadFastArray(fEX, fNpoints);
         b.ReadFastArray(fEY, fNpoints);
      }
      b.CheckByteCount(R__s, R__c, TGraphErrors::IsA());
      //====end of old versions

   } else {
      b.WriteClassBuffer(TGraphErrors::Class(), this);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Swap points

void TGraphErrors::SwapPoints(Int_t pos1, Int_t pos2)
{
   SwapValues(fEX, pos1, pos2);
   SwapValues(fEY, pos1, pos2);
   TGraph::SwapPoints(pos1, pos2);
}
