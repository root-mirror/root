// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss, Jan Therhaag

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Tools                                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch> - CERN, Switzerland           *
 *      Jan Therhaag       <Jan.Therhaag@cern.ch>     - U of Bonn, Germany        *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::Tools
\ingroup TMVA
Global auxiliary applications and data treatment routines.
*/

#include "TMVA/Tools.h"

#include "TMVA/Config.h"
#include "TMVA/Event.h"
#include "TMVA/Version.h"
#include "TMVA/PDF.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Types.h"

#include "TObjString.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TSpline.h"
#include "TVector.h"
#include "TMatrixD.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "TTreeFormula.h"
#include "TXMLEngine.h"
#include "TROOT.h"

#include <algorithm>
#include <cstdlib>
#include <iomanip>

using namespace std;

#if __cplusplus > 199711L
std::atomic<TMVA::Tools*> TMVA::Tools::fgTools{0};
#else
TMVA::Tools* TMVA::Tools::fgTools = 0;
#endif

TMVA::Tools& TMVA::gTools()                 { return TMVA::Tools::Instance(); }
TMVA::Tools& TMVA::Tools::Instance()        {
#if __cplusplus > 199711L
   if(!fgTools) {
      Tools* tmp = new Tools();
      Tools* expected = 0;
      if(! fgTools.compare_exchange_strong(expected,tmp)) {
         //another thread beat us
         delete tmp;
      }
   }
   return *fgTools;
#else
   return fgTools?*(fgTools): *(fgTools = new Tools());
#endif
}
void         TMVA::Tools::DestroyInstance() {
   //NOTE: there is no thread safe way to do this so
   // one must only call this method ones in an executable
#if __cplusplus > 199711L
   if (fgTools != 0) { delete fgTools.load(); fgTools=0; }
#else
   if (fgTools != 0) { delete fgTools; fgTools=0; }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::Tools::Tools() :
   fRegexp("$&|!%^&()'<>?= "),
   fLogger(new MsgLogger("Tools")),
   fXMLEngine(new TXMLEngine())
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::Tools::~Tools()
{
   delete fLogger;
   delete fXMLEngine;
}

////////////////////////////////////////////////////////////////////////////////
/// normalise to output range: [-1, 1]

Double_t TMVA::Tools::NormVariable( Double_t x, Double_t xmin, Double_t xmax )
{
   return 2*(x - xmin)/(xmax - xmin) - 1.0;
}

////////////////////////////////////////////////////////////////////////////////
/// compute "separation" defined as
/// \f[
/// <s2> = \frac{1}{2} \int_{-\infty}^{+\infty} \frac{(S(x) - B(x))^2}{(S(x) + B(x))} dx
/// \f]

Double_t TMVA::Tools::GetSeparation( TH1* S, TH1* B ) const
{
   Double_t separation = 0;

   // sanity checks
   // signal and background histograms must have same number of bins and
   // same limits
   if ((S->GetNbinsX() != B->GetNbinsX()) || (S->GetNbinsX() <= 0)) {
      Log() << kFATAL << "<GetSeparation> signal and background"
            << " histograms have different number of bins: "
            << S->GetNbinsX() << " : " << B->GetNbinsX() << Endl;
   }

   if (S->GetXaxis()->GetXmin() != B->GetXaxis()->GetXmin() ||
       S->GetXaxis()->GetXmax() != B->GetXaxis()->GetXmax() ||
       S->GetXaxis()->GetXmax() <= S->GetXaxis()->GetXmin()) {
      Log() << kINFO << S->GetXaxis()->GetXmin() << " " << B->GetXaxis()->GetXmin()
            << " " << S->GetXaxis()->GetXmax() << " " << B->GetXaxis()->GetXmax()
            << " " << S->GetXaxis()->GetXmax() << " " << S->GetXaxis()->GetXmin() << Endl;
      Log() << kFATAL << "<GetSeparation> signal and background"
            << " histograms have different or invalid dimensions:" << Endl;
   }

   Int_t    nstep  = S->GetNbinsX();
   Double_t intBin = (S->GetXaxis()->GetXmax() - S->GetXaxis()->GetXmin())/nstep;
   Double_t nS     = S->GetSumOfWeights()*intBin;
   Double_t nB     = B->GetSumOfWeights()*intBin;

   if (nS > 0 && nB > 0) {
      for (Int_t bin=0; bin<nstep; bin++) {
         Double_t s = S->GetBinContent( bin+1 )/Double_t(nS);
         Double_t b = B->GetBinContent( bin+1 )/Double_t(nB);
         // separation
         if (s + b > 0) separation += (s - b)*(s - b)/(s + b);
      }
      separation *= (0.5*intBin);
   }
   else {
      Log() << kWARNING << "<GetSeparation> histograms with zero entries: "
            << nS << " : " << nB << " cannot compute separation"
            << Endl;
      separation = 0;
   }

   return separation;
}

////////////////////////////////////////////////////////////////////////////////
/// compute "separation" defined as
/// \f[
/// <s2> = \frac{1}{2} \int_{-\infty}^{+\infty} \frac{(S(x) - B(x))^2}{(S(x) + B(x))} dx
/// \f]

Double_t TMVA::Tools::GetSeparation( const PDF& pdfS, const PDF& pdfB ) const
{
   Double_t xmin = pdfS.GetXmin();
   Double_t xmax = pdfS.GetXmax();
   // sanity check
   if (xmin != pdfB.GetXmin() || xmax != pdfB.GetXmax()) {
      Log() << kFATAL << "<GetSeparation> Mismatch in PDF limits: "
            << xmin << " " << pdfB.GetXmin() << xmax << " " << pdfB.GetXmax()  << Endl;
   }

   Double_t separation = 0;
   Int_t    nstep      = 100;
   Double_t intBin     = (xmax - xmin)/Double_t(nstep);
   for (Int_t bin=0; bin<nstep; bin++) {
      Double_t x = (bin + 0.5)*intBin + xmin;
      Double_t s = pdfS.GetVal( x );
      Double_t b = pdfB.GetVal( x );
      // separation
      if (s + b > 0) separation += (s - b)*(s - b)/(s + b);
   }
   separation *= (0.5*intBin);

   return separation;
}

////////////////////////////////////////////////////////////////////////////////
/// sanity check

void TMVA::Tools::ComputeStat( const std::vector<TMVA::Event*>& events, std::vector<Float_t>* valVec,
                               Double_t& meanS, Double_t& meanB,
                               Double_t& rmsS,  Double_t& rmsB,
                               Double_t& xmin,  Double_t& xmax,
                               Int_t signalClass, Bool_t  norm )
{
   if (0 == valVec)
      Log() << kFATAL << "<Tools::ComputeStat> value vector is zero pointer" << Endl;

   if ( events.size() != valVec->size() )
      Log() << kWARNING << "<Tools::ComputeStat> event and value vector have different lengths "
            << events.size() << "!=" << valVec->size() << Endl;

   Long64_t entries = valVec->size();

   // first fill signal and background in arrays before analysis
   Double_t* varVecS  = new Double_t[entries];
   Double_t* varVecB  = new Double_t[entries];
   Double_t* wgtVecS  = new Double_t[entries];
   Double_t* wgtVecB  = new Double_t[entries];
   xmin               = +DBL_MAX;
   xmax               = -DBL_MAX;
   Long64_t nEventsS  = 0;
   Long64_t nEventsB  = 0;
   Double_t xmin_ = 0, xmax_ = 0;

   if (norm) {
      xmin_ = *std::min( valVec->begin(), valVec->end() );
      xmax_ = *std::max( valVec->begin(), valVec->end() );
   }

   for (Int_t ievt=0; ievt<entries; ievt++) {
      Double_t theVar = (*valVec)[ievt];
      if (norm) theVar = Tools::NormVariable( theVar, xmin_, xmax_ );

      if (Int_t(events[ievt]->GetClass()) == signalClass ){
         wgtVecS[nEventsS]   = events[ievt]->GetWeight(); // this is signal
         varVecS[nEventsS++] = theVar; // this is signal
      }
      else {
         wgtVecB[nEventsB]   = events[ievt]->GetWeight(); // this is signal
         varVecB[nEventsB++] = theVar; // this is background
      }

      if (theVar > xmax) xmax = theVar;
      if (theVar < xmin) xmin = theVar;
   }
   // ++nEventsS;
   // ++nEventsB;

   // basic statistics
   // !!! TMath::Mean allows for weights, but NOT for negative weights
   //     and TMath::RMS doesn't allow for weights all together...
   meanS = TMVA::Tools::Mean( nEventsS, varVecS, wgtVecS );
   meanB = TMVA::Tools::Mean( nEventsB, varVecB, wgtVecB );
   rmsS  = TMVA::Tools::RMS ( nEventsS, varVecS, wgtVecS );
   rmsB  = TMVA::Tools::RMS ( nEventsB, varVecB, wgtVecB );

   delete [] varVecS;
   delete [] varVecB;
   delete [] wgtVecS;
   delete [] wgtVecB;
}

////////////////////////////////////////////////////////////////////////////////
/// square-root of symmetric matrix
/// of course the resulting sqrtMat is also symmetric, but it's easier to
/// treat it as a general matrix

TMatrixD* TMVA::Tools::GetSQRootMatrix( TMatrixDSym* symMat )
{
   Int_t n = symMat->GetNrows();

   // compute eigenvectors
   TMatrixDSymEigen* eigen = new TMatrixDSymEigen( *symMat );

   // D = ST C S
   TMatrixD* si = new TMatrixD( eigen->GetEigenVectors() );
   TMatrixD* s  = new TMatrixD( *si ); // copy
   si->Transpose( *si ); // invert (= transpose)

   // diagonal matrices
   TMatrixD* d = new TMatrixD( n, n);
   d->Mult( (*si), (*symMat) ); (*d) *= (*s);

   // sanity check: matrix must be diagonal and positive definit
   Int_t i, j;
   Double_t epsilon = 1.0e-8;
   for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
         if ((i != j && TMath::Abs((*d)(i,j))/TMath::Sqrt((*d)(i,i)*(*d)(j,j)) > epsilon) ||
             (i == j && (*d)(i,i) < 0)) {
            //d->Print();
            Log() << kWARNING << "<GetSQRootMatrix> error in matrix diagonalization; printed S and B" << Endl;
         }
      }
   }

   // make exactly diagonal
   for (i=0; i<n; i++) for (j=0; j<n; j++) if (j != i) (*d)(i,j) = 0;

   // compute the square-root C' of covariance matrix: C = C'*C'
   for (i=0; i<n; i++) (*d)(i,i) = TMath::Sqrt((*d)(i,i));

   TMatrixD* sqrtMat = new TMatrixD( n, n );
   sqrtMat->Mult( (*s), (*d) );
   (*sqrtMat) *= (*si);

   // invert square-root matrices
   sqrtMat->Invert();

   delete eigen;
   delete s;
   delete si;
   delete d;

   return sqrtMat;
}

////////////////////////////////////////////////////////////////////////////////
/// turns covariance into correlation matrix

const TMatrixD* TMVA::Tools::GetCorrelationMatrix( const TMatrixD* covMat )
{

   if (covMat == 0) return 0;
   // sanity check
   Int_t nvar = covMat->GetNrows();
   if (nvar != covMat->GetNcols())
      Log() << kFATAL << "<GetCorrelationMatrix> input matrix not quadratic" << Endl;

   Log() << kWARNING;
   TMatrixD* corrMat = new TMatrixD( nvar, nvar );
   for (Int_t ivar=0; ivar<nvar; ivar++) {
      for (Int_t jvar=0; jvar<nvar; jvar++) {
         if (ivar != jvar) {
            Double_t d = (*covMat)(ivar, ivar)*(*covMat)(jvar, jvar);
            if (d > 1E-20) {
               (*corrMat)(ivar, jvar) = (*covMat)(ivar, jvar)/TMath::Sqrt(d);
            } else {
               Log() <<  "<GetCorrelationMatrix> zero variances for variables "
                     << "(" << ivar << ", " << jvar << ")" << Endl;
               (*corrMat)(ivar, jvar) = 0;
            }
            if (TMath::Abs( (*corrMat)(ivar,jvar))  > 1){
               Log() << kWARNING
                     <<  " Element  corr("<<ivar<<","<<ivar<<")=" << (*corrMat)(ivar,jvar)
                     << " sigma2="<<d
                     << " cov("<<ivar<<","<<ivar<<")=" <<(*covMat)(ivar, ivar)
                     << " cov("<<jvar<<","<<jvar<<")=" <<(*covMat)(jvar, jvar)
                     << Endl;

            }
         }
         else (*corrMat)(ivar, ivar) = 1.0;
      }
   }
   Log() << Endl;
   return corrMat;
}

////////////////////////////////////////////////////////////////////////////////
/// projects variable from tree into normalised histogram

TH1* TMVA::Tools::projNormTH1F( TTree* theTree, const TString& theVarName,
                                const TString& name, Int_t nbins,
                                Double_t xmin, Double_t xmax, const TString& cut )
{
   // needed because of ROOT bug (feature) that excludes events that have value == xmax
   xmax += 0.00001;

   TH1* hist = new TH1F( name, name, nbins, xmin, xmax );
   hist->Sumw2(); // enable quadratic errors
   theTree->Project( name, theVarName, cut );
   NormHist( hist );
   return hist;
}

////////////////////////////////////////////////////////////////////////////////
/// normalises histogram

Double_t TMVA::Tools::NormHist( TH1* theHist, Double_t norm )
{
   if (!theHist) return 0;

   if (theHist->GetSumw2N() == 0) theHist->Sumw2();
   if (theHist->GetSumOfWeights() != 0) {
      Double_t   w  = ( theHist->GetSumOfWeights()
                        *(theHist->GetXaxis()->GetXmax() - theHist->GetXaxis()->GetXmin())/theHist->GetNbinsX() );
      if (w > 0) theHist->Scale( norm/w );
      return w;
   }

   return 1.0;
}

////////////////////////////////////////////////////////////////////////////////
/// Parse the string and cut into labels separated by ":"

TList* TMVA::Tools::ParseFormatLine( TString formatString, const char* sep )
{
   TList*   labelList = new TList();
   labelList->SetOwner();
   while (formatString.First(sep)==0) formatString.Remove(0,1); // remove initial separators

   while (formatString.Length()>0) {
      if (formatString.First(sep) == -1) { // no more separator
         labelList->Add(new TObjString(formatString.Data()));
         formatString="";
         break;
      }

      Ssiz_t posSep = formatString.First(sep);
      labelList->Add(new TObjString(TString(formatString(0,posSep)).Data()));
      formatString.Remove(0,posSep+1);

      while (formatString.First(sep)==0) formatString.Remove(0,1); // remove additional separators

   }
   return labelList;
}

////////////////////////////////////////////////////////////////////////////////
/// parse option string for ANN methods
/// default settings (should be defined in theOption string)
///
/// format and syntax of option string: "3000:N:N+2:N-3:6"
///
/// where:
///      -  3000 - number of training cycles (epochs)
///      -  N    - number of nodes in first hidden layer, where N is the number
///                of discriminating variables used (note that the first ANN
///                layer necessarily has N nodes, and hence is not given).
///      -  N+2  - number of nodes in 2nd hidden layer (2 nodes more than
///                number of variables)
///      -  N-3  - number of nodes in 3rd hidden layer (3 nodes less than
///                number of variables)
///      -  6    - 6 nodes in last (4th) hidden layer (note that the last ANN
///                layer in MVA has 2 nodes, each one for signal and background
///                classes)

vector<Int_t>* TMVA::Tools::ParseANNOptionString( TString theOptions, Int_t nvar,
                                                  vector<Int_t>* nodes )
{
   TList* list  = TMVA::Tools::ParseFormatLine( theOptions, ":" );


   // sanity check
   if (list->GetSize() < 1) {
      Log() << kFATAL << "<ParseANNOptionString> unrecognized option string: " << theOptions << Endl;
   }

   // add number of cycles
   nodes->push_back( atoi( ((TObjString*)list->At(0))->GetString() ) );

   Int_t a;
   if (list->GetSize() > 1) {
      for (Int_t i=1; i<list->GetSize(); i++) {
         TString s = ((TObjString*)list->At(i))->GetString();
         s.ToUpper();
         if (s(0) == 'N')  {
            if (s.Length() > 1) nodes->push_back( nvar + atoi(&s[1]) );
            else                nodes->push_back( nvar );
         }
         else if ((a = atoi( s )) > 0) nodes->push_back( atoi(s ) );
         else {
            Log() << kFATAL << "<ParseANNOptionString> unrecognized option string: " << theOptions << Endl;
         }
      }
   }

   return nodes;
}

////////////////////////////////////////////////////////////////////////////////
/// check quality of splining by comparing splines and histograms in each bin

Bool_t TMVA::Tools::CheckSplines( const TH1* theHist, const TSpline* theSpline )
{
   const Double_t sanityCrit = 0.01; // relative deviation

   Bool_t retval = kTRUE;
   for (Int_t ibin=1; ibin<=theHist->GetNbinsX(); ibin++) {
      Double_t x  = theHist->GetBinCenter( ibin );
      Double_t yh = theHist->GetBinContent( ibin ); // the histogram output
      Double_t ys = theSpline->Eval( x );           // the spline output

      if (ys + yh > 0) {
         Double_t dev = 0.5*(ys - yh)/(ys + yh);
         if (TMath::Abs(dev) > sanityCrit) {
            Log() << kFATAL << "<CheckSplines> Spline failed sanity criterion; "
                  << " relative deviation from histogram: " << dev
                  << " in (bin, value): (" << ibin << ", " << x << ")" << Endl;
            retval = kFALSE;
         }
      }
   }

   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// computes difference between two vectors

std::vector<Double_t> TMVA::Tools::MVADiff( std::vector<Double_t>& a, std::vector<Double_t>& b )
{
   if (a.size() != b.size()) {
      throw;
   }
   vector<Double_t> result(a.size());
   for (UInt_t i=0; i<a.size();i++) result[i]=a[i]-b[i];
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// scales double vector

void TMVA::Tools::Scale( std::vector<Double_t>& v, Double_t f )
{
   for (UInt_t i=0; i<v.size();i++) v[i]*=f;
}

////////////////////////////////////////////////////////////////////////////////
/// scales float vector

void TMVA::Tools::Scale( std::vector<Float_t>& v, Float_t f )
{
   for (UInt_t i=0; i<v.size();i++) v[i]*=f;
}

////////////////////////////////////////////////////////////////////////////////
/// sort 2D vector (AND in parallel a TString vector) in such a way
/// that the "first vector is sorted" and the other vectors are reshuffled
/// in the same way as necessary to have the first vector sorted.
/// I.e. the correlation between the elements is kept.

void TMVA::Tools::UsefulSortAscending( std::vector<vector<Double_t> >& v, std::vector<TString>* vs ){
   UInt_t nArrays=v.size();
   Double_t temp;
   if (nArrays > 0) {
      UInt_t sizeofarray=v[0].size();
      for (UInt_t i=0; i<sizeofarray; i++) {
         for (UInt_t j=sizeofarray-1; j>i; j--) {
            if (v[0][j-1] > v[0][j]) {
               for (UInt_t k=0; k< nArrays; k++) {
                  temp = v[k][j-1]; v[k][j-1] = v[k][j]; v[k][j] = temp;
               }
               if (nullptr != vs) {
                  TString temps = (*vs)[j-1]; (*vs)[j-1] = (*vs)[j]; (*vs)[j] = temps;
               }
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// sort 2D vector (AND in parallel a TString vector) in such a way
/// that the "first vector is sorted" and the other vectors are reshuffled
/// in the same way as necessary to have the first vector sorted.
/// I.e. the correlation between the elements is kept.

void TMVA::Tools::UsefulSortDescending( std::vector<std::vector<Double_t> >& v, std::vector<TString>* vs )
{
   UInt_t nArrays=v.size();
   Double_t temp;
   if (nArrays > 0) {
      UInt_t sizeofarray=v[0].size();
      for (UInt_t i=0; i<sizeofarray; i++) {
         for (UInt_t j=sizeofarray-1; j>i; j--) {
            if (v[0][j-1] < v[0][j]) {
               for (UInt_t k=0; k< nArrays; k++) {
                  temp = v[k][j-1]; v[k][j-1] = v[k][j]; v[k][j] = temp;
               }
               if (nullptr != vs) {
                  TString temps = (*vs)[j-1]; (*vs)[j-1] = (*vs)[j]; (*vs)[j] = temps;
               }
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Mutual Information method for non-linear correlations estimates in 2D histogram
/// Author: Moritz Backes, Geneva (2009)

Double_t TMVA::Tools::GetMutualInformation( const TH2F& h_ )
{
   Double_t hi = h_.Integral();
   if (hi == 0) return -1;

   // copy histogram and rebin to speed up procedure
   TH2F h( h_ );
   h.RebinX(2);
   h.RebinY(2);

   Double_t mutualInfo = 0.;
   Int_t maxBinX = h.GetNbinsX();
   Int_t maxBinY = h.GetNbinsY();
   for (Int_t x = 1; x <= maxBinX; x++) {
      for (Int_t y = 1; y <= maxBinY; y++) {
         Double_t p_xy = h.GetBinContent(x,y)/hi;
         Double_t p_x  = h.Integral(x,x,1,maxBinY)/hi;
         Double_t p_y  = h.Integral(1,maxBinX,y,y)/hi;
         if (p_x > 0. && p_y > 0. && p_xy > 0.){
            mutualInfo += p_xy*TMath::Log(p_xy / (p_x * p_y));
         }
      }
   }

   return mutualInfo;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute Correlation Ratio of 2D histogram to estimate functional dependency between two variables
/// Author: Moritz Backes, Geneva (2009)

Double_t TMVA::Tools::GetCorrelationRatio( const TH2F& h_ )
{
   Double_t hi = h_.Integral();
   if (hi == 0.) return -1;

   // copy histogram and rebin to speed up procedure
   TH2F h( h_ );
   h.RebinX(2);
   h.RebinY(2);

   Double_t corrRatio = 0.;
   Double_t y_mean = h.ProjectionY()->GetMean();
   for (Int_t ix=1; ix<=h.GetNbinsX(); ix++) {
      corrRatio += (h.Integral(ix,ix,1,h.GetNbinsY())/hi)*pow((GetYMean_binX(h,ix)-y_mean),2);
   }
   corrRatio /= pow(h.ProjectionY()->GetRMS(),2);
   return corrRatio;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the mean in Y for a given bin X of a 2D histogram

Double_t TMVA::Tools::GetYMean_binX( const TH2& h, Int_t bin_x )
{
   if (h.Integral(bin_x,bin_x,1,h.GetNbinsY()) == 0.) {return 0;}
   Double_t y_bin_mean = 0.;
   TH1* py = h.ProjectionY();
   for (Int_t y = 1; y <= h.GetNbinsY(); y++){
      y_bin_mean += h.GetBinContent(bin_x,y)*py->GetBinCenter(y);
   }
   y_bin_mean /= h.Integral(bin_x,bin_x,1,h.GetNbinsY());
   return y_bin_mean;
}

////////////////////////////////////////////////////////////////////////////////
/// Transpose quadratic histogram

TH2F* TMVA::Tools::TransposeHist( const TH2F& h )
{
   // sanity check
   if (h.GetNbinsX() != h.GetNbinsY()) {
      Log() << kFATAL << "<TransposeHist> cannot transpose non-quadratic histogram" << Endl;
   }

   TH2F *transposedHisto = new TH2F( h );
   for (Int_t ix=1; ix <= h.GetNbinsX(); ix++){
      for (Int_t iy=1; iy <= h.GetNbinsY(); iy++){
         transposedHisto->SetBinContent(iy,ix,h.GetBinContent(ix,iy));
      }
   }

   // copy stats (thanks to Swagato Banerjee for pointing out the missing stats information)
   Double_t stats_old[7];
   Double_t stats_new[7];

   h.GetStats(stats_old);
   stats_new[0] = stats_old[0];
   stats_new[1] = stats_old[1];
   stats_new[2] = stats_old[4];
   stats_new[3] = stats_old[5];
   stats_new[4] = stats_old[2];
   stats_new[5] = stats_old[3];
   stats_new[6] = stats_old[6];
   transposedHisto->PutStats(stats_new);

   return transposedHisto; // ownership returned
}

////////////////////////////////////////////////////////////////////////////////
/// check for "silence" option in configuration option string

Bool_t TMVA::Tools::CheckForSilentOption( const TString& cs ) const
{
   Bool_t isSilent = kFALSE;

   TString s( cs );
   s.ToLower();
   s.ReplaceAll(" ","");
   if (s.Contains("silent") && !s.Contains("silent=f")) {
      if (!s.Contains("!silent") || s.Contains("silent=t")) isSilent = kTRUE;
   }

   return isSilent;
}

////////////////////////////////////////////////////////////////////////////////
/// check if verbosity "V" set in option

Bool_t TMVA::Tools::CheckForVerboseOption( const TString& cs ) const
{
   Bool_t isVerbose = kFALSE;

   TString s( cs );
   s.ToLower();
   s.ReplaceAll(" ","");
   std::vector<TString> v = SplitString( s, ':' );
   for (std::vector<TString>::iterator it = v.begin(); it != v.end(); ++it) {
      if ((*it == "v" || *it == "verbose") && !it->Contains("!")) isVerbose = kTRUE;
   }

   return isVerbose;
}

////////////////////////////////////////////////////////////////////////////////
/// sort vector

void TMVA::Tools::UsefulSortDescending( std::vector<Double_t>& v )
{
   vector< vector<Double_t> > vtemp;
   vtemp.push_back(v);
   UsefulSortDescending(vtemp);
   v = vtemp[0];
}

////////////////////////////////////////////////////////////////////////////////
/// sort vector

void TMVA::Tools::UsefulSortAscending( std::vector<Double_t>& v )
{
   vector<vector<Double_t> > vtemp;
   vtemp.push_back(v);
   UsefulSortAscending(vtemp);
   v = vtemp[0];
}

////////////////////////////////////////////////////////////////////////////////
/// find index of maximum entry in vector

Int_t TMVA::Tools::GetIndexMaxElement( std::vector<Double_t>& v )
{
   if (v.empty()) return -1;

   Int_t pos=0; Double_t mx=v[0];
   for (UInt_t i=0; i<v.size(); i++){
      if (v[i] > mx){
         mx=v[i];
         pos=i;
      }
   }
   return pos;
}

////////////////////////////////////////////////////////////////////////////////
/// find index of minimum entry in vector

Int_t TMVA::Tools::GetIndexMinElement( std::vector<Double_t>& v )
{
   if (v.empty()) return -1;

   Int_t pos=0; Double_t mn=v[0];
   for (UInt_t i=0; i<v.size(); i++){
      if (v[i] < mn){
         mn=v[i];
         pos=i;
      }
   }
   return pos;
}


////////////////////////////////////////////////////////////////////////////////
/// check if regular expression
/// helper function to search for "$!%^&()'<>?= " in a string

Bool_t TMVA::Tools::ContainsRegularExpression( const TString& s )
{
   Bool_t  regular = kFALSE;
   for (Int_t i = 0; i < Tools::fRegexp.Length(); i++)
      if (s.Contains( Tools::fRegexp[i] )) { regular = kTRUE; break; }

   return regular;
}

////////////////////////////////////////////////////////////////////////////////
/// replace regular expressions
/// helper function to remove all occurrences "$!%^&()'<>?= " from a string
/// and replace all ::,$,*,/,+,- with _M_,_S_,_T_,_D_,_P_,_M_ respectively

TString TMVA::Tools::ReplaceRegularExpressions( const TString& s, const TString& r )
{
   TString snew = s;
   for (Int_t i = 0; i < Tools::fRegexp.Length(); i++)
      snew.ReplaceAll( Tools::fRegexp[i], r );

   snew.ReplaceAll( "::", r );
   snew.ReplaceAll( "$", "_S_" );
   snew.ReplaceAll( "&", "_A_" );
   snew.ReplaceAll( "%", "_MOD_" );
   snew.ReplaceAll( "|", "_O_" );
   snew.ReplaceAll( "*", "_T_" );
   snew.ReplaceAll( "/", "_D_" );
   snew.ReplaceAll( "+", "_P_" );
   snew.ReplaceAll( "-", "_M_" );
   snew.ReplaceAll( " ", "_" );
   snew.ReplaceAll( "[", "_" );
   snew.ReplaceAll( "]", "_" );
   snew.ReplaceAll( "=", "_E_" );
   snew.ReplaceAll( ">", "_GT_" );
   snew.ReplaceAll( "<", "_LT_" );
   snew.ReplaceAll( "(", "_" );
   snew.ReplaceAll( ")", "_" );

   return snew;
}

////////////////////////////////////////////////////////////////////////////////
/// human readable color strings

const TString& TMVA::Tools::Color( const TString& c )
{
   static const TString gClr_none         = "" ;
   static const TString gClr_white        = "\033[1;37m";  // white
   static const TString gClr_black        = "\033[30m";    // black
   static const TString gClr_blue         = "\033[34m";    // blue
   static const TString gClr_red          = "\033[1;31m" ; // red
   static const TString gClr_yellow       = "\033[1;33m";  // yellow
   static const TString gClr_darkred      = "\033[31m";    // dark red
   static const TString gClr_darkgreen    = "\033[32m";    // dark green
   static const TString gClr_darkyellow   = "\033[33m";    // dark yellow

   static const TString gClr_bold         = "\033[1m"    ; // bold
   static const TString gClr_black_b      = "\033[30m"   ; // bold black
   static const TString gClr_lblue_b      = "\033[1;34m" ; // bold light blue
   static const TString gClr_cyan_b       = "\033[0;36m" ; // bold cyan
   static const TString gClr_lgreen_b     = "\033[1;32m";  // bold light green

   static const TString gClr_blue_bg      = "\033[44m";    // blue background
   static const TString gClr_red_bg       = "\033[1;41m";  // white on red background
   static const TString gClr_whiteonblue  = "\033[1;44m";  // white on blue background
   static const TString gClr_whiteongreen = "\033[1;42m";  // white on green background
   static const TString gClr_grey_bg      = "\033[47m";    // grey background

   static const TString gClr_reset  = "\033[0m";     // reset

   if (!gConfig().UseColor()) return gClr_none;

   if (c == "white" )         return gClr_white;
   if (c == "blue"  )         return gClr_blue;
   if (c == "black"  )        return gClr_black;
   if (c == "lightblue")      return gClr_cyan_b;
   if (c == "yellow")         return gClr_yellow;
   if (c == "red"   )         return gClr_red;
   if (c == "dred"  )         return gClr_darkred;
   if (c == "dgreen")         return gClr_darkgreen;
   if (c == "lgreenb")        return gClr_lgreen_b;
   if (c == "dyellow")        return gClr_darkyellow;

   if (c == "bold")           return gClr_bold;
   if (c == "bblack")         return gClr_black_b;

   if (c == "blue_bgd")       return gClr_blue_bg;
   if (c == "red_bgd" )       return gClr_red_bg;

   if (c == "white_on_blue" ) return gClr_whiteonblue;
   if (c == "white_on_green") return gClr_whiteongreen;

   if (c == "reset") return gClr_reset;

   std::cout << "Unknown color " << c << std::endl;
   exit(1);

   return gClr_none;
}

////////////////////////////////////////////////////////////////////////////////
/// formatted output of simple table

void TMVA::Tools::FormattedOutput( const std::vector<Double_t>& values, const std::vector<TString>& V,
                                   const TString titleVars, const TString titleValues, MsgLogger& logger,
                                   TString format )
{
   // sanity check
   UInt_t nvar = V.size();
   if ((UInt_t)values.size() != nvar) {
      logger << kFATAL << "<FormattedOutput> fatal error with dimensions: "
             << values.size() << " OR " << " != " << nvar << Endl;
   }

   // find maximum length in V (and column title)
   UInt_t maxL = 7;
   std::vector<UInt_t> vLengths;
   for (UInt_t ivar=0; ivar<nvar; ivar++) maxL = TMath::Max( (UInt_t)V[ivar].Length(), maxL );
   maxL = TMath::Max( (UInt_t)titleVars.Length(), maxL );

   // column length
   UInt_t maxV = 7;
   maxV = TMath::Max( (UInt_t)titleValues.Length() + 1, maxL );

   // full column length
   UInt_t clen = maxL + maxV + 3;

   // bar line
   for (UInt_t i=0; i<clen; i++) logger << "-";
   logger << Endl;

   // title bar
   logger << setw(maxL) << titleVars << ":";
   logger << setw(maxV+1) << titleValues << ":";
   logger << Endl;
   for (UInt_t i=0; i<clen; i++) logger << "-";
   logger << Endl;

   // the numbers
   for (UInt_t irow=0; irow<nvar; irow++) {
      logger << setw(maxL) << V[irow] << ":";
      logger << setw(maxV+1) << Form( format.Data(), values[irow] );
      logger << Endl;
   }

   // bar line
   for (UInt_t i=0; i<clen; i++) logger << "-";
   logger << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// formatted output of matrix (with labels)

void TMVA::Tools::FormattedOutput( const TMatrixD& M, const std::vector<TString>& V, MsgLogger& logger )
{
   // sanity check: matrix must be quadratic
   UInt_t nvar = V.size();
   if ((UInt_t)M.GetNcols() != nvar || (UInt_t)M.GetNrows() != nvar) {
      logger << kFATAL << "<FormattedOutput> fatal error with dimensions: "
             << M.GetNcols() << " OR " << M.GetNrows() << " != " << nvar << " ==> abort" << Endl;
   }

   // get length of each variable, and maximum length
   UInt_t minL = 7;
   UInt_t maxL = minL;
   std::vector<UInt_t> vLengths;
   for (UInt_t ivar=0; ivar<nvar; ivar++) {
      vLengths.push_back(TMath::Max( (UInt_t)V[ivar].Length(), minL ));
      maxL = TMath::Max( vLengths.back(), maxL );
   }

   // count column length
   UInt_t clen = maxL+1;
   for (UInt_t icol=0; icol<nvar; icol++) clen += vLengths[icol]+1;

   // bar line
   for (UInt_t i=0; i<clen; i++) logger << "-";
   logger << Endl;

   // title bar
   logger << setw(maxL+1) << " ";
   for (UInt_t icol=0; icol<nvar; icol++) logger << setw(vLengths[icol]+1) << V[icol];
   logger << Endl;

   // the numbers
   for (UInt_t irow=0; irow<nvar; irow++) {
      logger << setw(maxL) << V[irow] << ":";
      for (UInt_t icol=0; icol<nvar; icol++) {
         logger << setw(vLengths[icol]+1) << Form( "%+1.3f", M(irow,icol) );
      }
      logger << Endl;
   }

   // bar line
   for (UInt_t i=0; i<clen; i++) logger << "-";
   logger << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// formatted output of matrix (with labels)

void TMVA::Tools::FormattedOutput( const TMatrixD& M,
                                   const std::vector<TString>& vert, const std::vector<TString>& horiz,
                                   MsgLogger& logger )
{
   // sanity check: matrix must be quadratic
   UInt_t nvvar = vert.size();
   UInt_t nhvar = horiz.size();

   // get length of each variable, and maximum length
   UInt_t minL = 7;
   UInt_t maxL = minL;
   std::vector<UInt_t> vLengths;
   for (UInt_t ivar=0; ivar<nvvar; ivar++) {
      vLengths.push_back(TMath::Max( (UInt_t)vert[ivar].Length(), minL ));
      maxL = TMath::Max( vLengths.back(), maxL );
   }

   // count column length
   UInt_t minLh = 7;
   UInt_t maxLh = minLh;
   std::vector<UInt_t> hLengths;
   for (UInt_t ivar=0; ivar<nhvar; ivar++) {
      hLengths.push_back(TMath::Max( (UInt_t)horiz[ivar].Length(), minL ));
      maxLh = TMath::Max( hLengths.back(), maxLh );
   }

   UInt_t clen = maxLh+1;
   for (UInt_t icol=0; icol<nhvar; icol++) clen += hLengths[icol]+1;

   // bar line
   for (UInt_t i=0; i<clen; i++) logger << "-";
   logger << Endl;

   // title bar
   logger << setw(maxL+1) << " ";
   for (UInt_t icol=0; icol<nhvar; icol++) logger << setw(hLengths[icol]+1) << horiz[icol];
   logger << Endl;

   // the numbers
   for (UInt_t irow=0; irow<nvvar; irow++) {
      logger << setw(maxL) << vert[irow] << ":";
      for (UInt_t icol=0; icol<nhvar; icol++) {
         logger << setw(hLengths[icol]+1) << Form( "%+1.3f", M(irow,icol) );
      }
      logger << Endl;
   }

   // bar line
   for (UInt_t i=0; i<clen; i++) logger << "-";
   logger << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// histogramming utility

TString TMVA::Tools::GetXTitleWithUnit( const TString& title, const TString& unit )
{
   return ( unit == "" ? title : ( title + "  [" + unit + "]" ) );
}

////////////////////////////////////////////////////////////////////////////////
/// histogramming utility

TString TMVA::Tools::GetYTitleWithUnit( const TH1& h, const TString& unit, Bool_t normalised )
{
   TString retval = ( normalised ? "(1/N) " : "" );
   retval += Form( "dN_{ }/^{ }%.3g %s", h.GetXaxis()->GetBinWidth(1), unit.Data() );
   return retval;
}

////////////////////////////////////////////////////////////////////////////////
/// writes a float value with the available precision to a stream

void TMVA::Tools::WriteFloatArbitraryPrecision( Float_t val, ostream& os )
{
   os << val << " :: ";
   void * c = &val;
   for (int i=0; i<4; i++) {
      Int_t ic = *((char*)c+i)-'\0';
      if (ic<0) ic+=256;
      os << ic << " ";
   }
   os << ":: ";
}

////////////////////////////////////////////////////////////////////////////////
/// reads a float value with the available precision from a stream

void TMVA::Tools::ReadFloatArbitraryPrecision( Float_t& val, istream& is )
{
   Float_t a = 0;
   is >> a;
   TString dn;
   is >> dn;
   Int_t c[4];
   void * ap = &a;
   for (int i=0; i<4; i++) {
      is >> c[i];
      *((char*)ap+i) = '\0'+c[i];
   }
   is >> dn;
   val = a;
}

// XML file reading/writing helper functions

////////////////////////////////////////////////////////////////////////////////
/// add attribute from xml

Bool_t TMVA::Tools::HasAttr( void* node, const char* attrname )
{
   return xmlengine().HasAttr(node, attrname);
}

////////////////////////////////////////////////////////////////////////////////
/// add attribute from xml

void TMVA::Tools::ReadAttr( void* node, const char* attrname, TString& value )
{
   if (!HasAttr(node, attrname)) {
      const char * nodename = xmlengine().GetNodeName(node);
      Log() << kFATAL << "Trying to read non-existing attribute '" << attrname << "' from xml node '" << nodename << "'" << Endl;
   }
   const char* val = xmlengine().GetAttr(node, attrname);
   value = TString(val);
}

////////////////////////////////////////////////////////////////////////////////
/// add attribute to node

void TMVA::Tools::AddAttr( void* node, const char* attrname, const char* value )
{
   if( node == 0 ) return;
   gTools().xmlengine().NewAttr(node, 0, attrname, value );
}

////////////////////////////////////////////////////////////////////////////////
/// add child node

void* TMVA::Tools::AddChild( void* parent, const char* childname, const char* content, bool isRootNode )
{
   if( !isRootNode && parent == 0 ) return 0;
   return gTools().xmlengine().NewChild(parent, 0, childname, content);
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TMVA::Tools::AddComment( void* node, const char* comment ) {
   if( node == 0 ) return kFALSE;
   return gTools().xmlengine().AddComment(node, comment);
}

////////////////////////////////////////////////////////////////////////////////
/// get parent node

void* TMVA::Tools::GetParent( void* child)
{
   void* par = xmlengine().GetParent(child);

   return par;
}

////////////////////////////////////////////////////////////////////////////////
/// get child node

void* TMVA::Tools::GetChild( void* parent, const char* childname )
{
   void* ch = xmlengine().GetChild(parent);
   if (childname != 0) {
      while (ch!=0 && strcmp(xmlengine().GetNodeName(ch),childname) != 0) ch = xmlengine().GetNext(ch);
   }
   return ch;
}

////////////////////////////////////////////////////////////////////////////////
/// XML helpers

void* TMVA::Tools::GetNextChild( void* prevchild, const char* childname )
{
   void* ch = xmlengine().GetNext(prevchild);
   if (childname != 0) {
      while (ch!=0 && strcmp(xmlengine().GetNodeName(ch),childname)!=0) ch = xmlengine().GetNext(ch);
   }
   return ch;
}

////////////////////////////////////////////////////////////////////////////////
/// XML helpers

const char* TMVA::Tools::GetContent( void* node )
{
   return xmlengine().GetNodeContent(node);
}

////////////////////////////////////////////////////////////////////////////////
/// XML helpers

const char* TMVA::Tools::GetName( void* node )
{
   return xmlengine().GetNodeName(node);
}

////////////////////////////////////////////////////////////////////////////////
/// XML helpers

Bool_t TMVA::Tools::AddRawLine( void* node, const char * raw )
{
   return xmlengine().AddRawLine( node, raw );
}

////////////////////////////////////////////////////////////////////////////////
/// splits the option string at 'separator' and fills the list
/// 'splitV' with the primitive strings

std::vector<TString> TMVA::Tools::SplitString(const TString& theOpt, const char separator ) const
{
   std::vector<TString> splitV;
   TString splitOpt(theOpt);
   splitOpt.ReplaceAll("\n"," ");
   splitOpt = splitOpt.Strip(TString::kBoth,separator);
   while (splitOpt.Length()>0) {
      if ( !splitOpt.Contains(separator) ) {
         splitV.push_back(splitOpt);
         break;
      }
      else {
         TString toSave = splitOpt(0,splitOpt.First(separator));
         splitV.push_back(toSave);
         splitOpt = splitOpt(splitOpt.First(separator),splitOpt.Length());
      }
      splitOpt = splitOpt.Strip(TString::kLeading,separator);
   }
   return splitV;
}

////////////////////////////////////////////////////////////////////////////////
/// string tools

TString TMVA::Tools::StringFromInt( Long_t i )
{
   std::stringstream s;
   s << i;
   return TString(s.str().c_str());
}

////////////////////////////////////////////////////////////////////////////////
/// string tools

TString TMVA::Tools::StringFromDouble( Double_t d )
{
   std::stringstream s;
   s << Form( "%5.8e", d );
   return TString(s.str().c_str());
}

////////////////////////////////////////////////////////////////////////////////
/// XML helpers

void TMVA::Tools::WriteTMatrixDToXML( void* node, const char* name, TMatrixD* mat )
{
   void* matnode = xmlengine().NewChild(node, 0, name);
   xmlengine().NewAttr(matnode,0,"Rows", StringFromInt(mat->GetNrows()) );
   xmlengine().NewAttr(matnode,0,"Columns", StringFromInt(mat->GetNcols()) );
   std::stringstream s;
   for (Int_t row = 0; row<mat->GetNrows(); row++) {
      for (Int_t col = 0; col<mat->GetNcols(); col++) {
         s << Form( "%5.15e ", (*mat)[row][col] );
      }
   }
   xmlengine().AddRawLine( matnode, s.str().c_str() );
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::Tools::WriteTVectorDToXML( void* node, const char* name, TVectorD* vec )
{
   TMatrixD mat(1,vec->GetNoElements(),&((*vec)[0]));
   WriteTMatrixDToXML( node, name, &mat );
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::Tools::ReadTVectorDFromXML( void* node, const char* name, TVectorD* vec )
{
   TMatrixD mat(1,vec->GetNoElements(),&((*vec)[0]));
   ReadTMatrixDFromXML( node, name, &mat );
   for (int i=0;i<vec->GetNoElements();++i) (*vec)[i] = mat[0][i];
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::Tools::ReadTMatrixDFromXML( void* node, const char* name, TMatrixD* mat )
{
   if (strcmp(xmlengine().GetNodeName(node),name)!=0){
      Log() << kWARNING << "Possible Error: Name of matrix in weight file"
            << " does not match name of matrix passed as argument!" << Endl;
   }
   Int_t nrows, ncols;
   ReadAttr( node, "Rows",    nrows );
   ReadAttr( node, "Columns", ncols );
   if (mat->GetNrows() != nrows || mat->GetNcols() != ncols){
      Log() << kWARNING << "Possible Error: Dimension of matrix in weight file"
            << " does not match dimension of matrix passed as argument!" << Endl;
      mat->ResizeTo(nrows,ncols);
   }
   const char* content = xmlengine().GetNodeContent(node);
   std::stringstream s(content);
   for (Int_t row = 0; row<nrows; row++) {
      for (Int_t col = 0; col<ncols; col++) {
         s >> (*mat)[row][col];
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// direct output, eg, when starting ROOT session -> no use of Logger here

void TMVA::Tools::TMVAWelcomeMessage()
{
   std::cout << std::endl;
   std::cout << Color("bold") << "TMVA -- Toolkit for Multivariate Data Analysis" << Color("reset") << std::endl;
   std::cout << "        " << "Version " << TMVA_RELEASE << ", " << TMVA_RELEASE_DATE << std::endl;
   std::cout << "        " << "Copyright (C) 2005-2010 CERN, MPI-K Heidelberg, Us of Bonn and Victoria" << std::endl;
   std::cout << "        " << "Home page:     http://tmva.sf.net" << std::endl;
   std::cout << "        " << "Citation info: http://tmva.sf.net/citeTMVA.html" << std::endl;
   std::cout << "        " << "License:       http://tmva.sf.net/LICENSE" << std::endl << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// prints the TMVA release number and date

void TMVA::Tools::TMVAVersionMessage( MsgLogger& logger )
{
   logger << "___________TMVA Version " << TMVA_RELEASE << ", " << TMVA_RELEASE_DATE
  << "" << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// prints the ROOT release number and date

void TMVA::Tools::ROOTVersionMessage( MsgLogger& logger )
{
   static const char * const months[] = { "Jan","Feb","Mar","Apr","May",
                                          "Jun","Jul","Aug","Sep","Oct",
                                          "Nov","Dec" };
   Int_t   idatqq = gROOT->GetVersionDate();
   Int_t   iday   = idatqq%100;
   Int_t   imonth = (idatqq/100)%100;
   Int_t   iyear  = (idatqq/10000);
   TString versionDate = Form("%s %d, %4d",months[imonth-1],iday,iyear);

   logger << kHEADER ;
   logger << "You are running ROOT Version: " << gROOT->GetVersion() << ", " << versionDate << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// various kinds of welcome messages
/// ASCII text generated by this site: http://www.network-science.de/ascii/

void TMVA::Tools::TMVAWelcomeMessage( MsgLogger& logger, EWelcomeMessage msgType )
{
   switch (msgType) {

   case kStandardWelcomeMsg:
      logger << Color("white") << "TMVA -- Toolkit for Multivariate Analysis" << Color("reset") << Endl;
      logger << "Copyright (C) 2005-2006 CERN, LAPP & MPI-K Heidelberg and Victoria U." << Endl;
      logger << "Home page http://tmva.sourceforge.net" << Endl;
      logger << "All rights reserved, please read http://tmva.sf.net/license.txt" << Endl << Endl;
      break;

   case kIsometricWelcomeMsg:
      logger << "   ___           ___           ___           ___      " << Endl;
      logger << "  /\\  \\         /\\__\\         /\\__\\         /\\  \\     " << Endl;
      logger << "  \\:\\  \\       /::|  |       /:/  /        /::\\  \\    " << Endl;
      logger << "   \\:\\  \\     /:|:|  |      /:/  /        /:/\\:\\  \\   " << Endl;
      logger << "   /::\\  \\   /:/|:|__|__   /:/__/  ___   /::\\~\\:\\  \\  " << Endl;
      logger << "  /:/\\:\\__\\ /:/ |::::\\__\\  |:|  | /\\__\\ /:/\\:\\ \\:\\__\\ " << Endl;
      logger << " /:/  \\/__/ \\/__/~~/:/  /  |:|  |/:/  / \\/__\\:\\/:/  / " << Endl;
      logger << "/:/  /            /:/  /   |:|__/:/  /       \\::/  /  " << Endl;
      logger << "\\/__/            /:/  /     \\::::/__/        /:/  /   " << Endl;
      logger << "                /:/  /       ~~~~           /:/  /    " << Endl;
      logger << "                \\/__/                       \\/__/     " << Endl << Endl;
      break;

   case kBlockWelcomeMsg:
      logger << Endl;
      logger << "_|_|_|_|_|  _|      _|  _|      _|    _|_|    " << Endl;
      logger << "    _|      _|_|  _|_|  _|      _|  _|    _|  " << Endl;
      logger << "    _|      _|  _|  _|  _|      _|  _|_|_|_|  " << Endl;
      logger << "    _|      _|      _|    _|  _|    _|    _|  " << Endl;
      logger << "    _|      _|      _|      _|      _|    _|  " << Endl << Endl;
      break;

   case kLeanWelcomeMsg:
      logger << Endl;
      logger << "_/_/_/_/_/  _/      _/  _/      _/    _/_/   " << Endl;
      logger << "   _/      _/_/  _/_/  _/      _/  _/    _/  " << Endl;
      logger << "  _/      _/  _/  _/  _/      _/  _/_/_/_/   " << Endl;
      logger << " _/      _/      _/    _/  _/    _/    _/    " << Endl;
      logger << "_/      _/      _/      _/      _/    _/     " << Endl << Endl;
      break;

   case kLogoWelcomeMsg:
      logger << Endl;
      logger << "_/_/_/_/_/ _|      _|  _|      _|    _|_|   " << Endl;
      logger << "   _/      _|_|  _|_|  _|      _|  _|    _| " << Endl;
      logger << "  _/       _|  _|  _|  _|      _|  _|_|_|_| " << Endl;
      logger << " _/        _|      _|    _|  _|    _|    _| " << Endl;
      logger << "_/         _|      _|      _|      _|    _| " << Endl << Endl;
      break;

   case kSmall1WelcomeMsg:
      logger << " _____ __  ____   ___   " << Endl;
      logger << "|_   _|  \\/  \\ \\ / /_\\  " << Endl;
      logger << "  | | | |\\/| |\\ V / _ \\ " << Endl;
      logger << "  |_| |_|  |_| \\_/_/ \\_\\" << Endl << Endl;
      break;

   case kSmall2WelcomeMsg:
      logger << " _____ __  ____     ___     " << Endl;
      logger << "|_   _|  \\/  \\ \\   / / \\    " << Endl;
      logger << "  | | | |\\/| |\\ \\ / / _ \\   " << Endl;
      logger << "  | | | |  | | \\ V / ___ \\  " << Endl;
      logger << "  |_| |_|  |_|  \\_/_/   \\_\\ " << Endl << Endl;
      break;

   case kOriginalWelcomeMsgColor:
      logger << kINFO << "" << Color("red")
             << "_______________________________________" << Color("reset") << Endl;
      logger << kINFO << "" << Color("blue")
             << Color("red_bgd") << Color("bwhite") << " // " << Color("reset")
             << Color("white") << Color("blue_bgd")
             << "|\\  /|| \\  //  /\\\\\\\\\\\\\\\\\\\\\\\\ \\ \\ \\ " << Color("reset") << Endl;
      logger << kINFO << ""<< Color("blue")
             << Color("red_bgd") << Color("white") << "//  " << Color("reset")
             << Color("white") << Color("blue_bgd")
             << "| \\/ ||  \\//  /--\\\\\\\\\\\\\\\\\\\\\\\\ \\ \\ \\" << Color("reset") << Endl;
      break;

   case kOriginalWelcomeMsgBW:
      logger << kINFO << ""
             << "_______________________________________" << Endl;
      logger << kINFO << " // "
             << "|\\  /|| \\  //  /\\\\\\\\\\\\\\\\\\\\\\\\ \\ \\ \\ " << Endl;
      logger << kINFO << "//  "
             << "| \\/ ||  \\//  /--\\\\\\\\\\\\\\\\\\\\\\\\ \\ \\ \\" << Endl;
      break;

   default:
      logger << kFATAL << "unknown message type: " << msgType << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// kinds of TMVA citation

void TMVA::Tools::TMVACitation( MsgLogger& logger, ECitation citType )
{
   switch (citType) {

   case kPlainText:
      logger << "A. Hoecker, P. Speckmayer, J. Stelzer, J. Therhaag, E. von Toerne, H. Voss" << Endl;
      logger << "\"TMVA - Toolkit for Multivariate Data Analysis\" PoS ACAT:040,2007. e-Print: physics/0703039" << Endl;
      break;

   case kBibTeX:
      logger << "@Article{TMVA2007," << Endl;
      logger << "     author    = \"Hoecker, Andreas and Speckmayer, Peter and Stelzer, Joerg " << Endl;
      logger << "                   and Therhaag, Jan and von Toerne, Eckhard and Voss, Helge\"," << Endl;
      logger << "     title     = \"{TMVA: Toolkit for multivariate data analysis}\"," << Endl;
      logger << "     journal   = \"PoS\"," << Endl;
      logger << "     volume    = \"ACAT\"," << Endl;
      logger << "     year      = \"2007\"," << Endl;
      logger << "     pages     = \"040\"," << Endl;
      logger << "     eprint    = \"physics/0703039\"," << Endl;
      logger << "     archivePrefix = \"arXiv\"," << Endl;
      logger << "     SLACcitation  = \"%%CITATION = PHYSICS/0703039;%%\"" << Endl;
      logger << "}" << Endl;
      break;

   case kLaTeX:
      logger << "%\\cite{TMVA2007}" << Endl;
      logger << "\\bibitem{TMVA2007}" << Endl;
      logger << "  A.~Hoecker, P.~Speckmayer, J.~Stelzer, J.~Therhaag, E.~von Toerne, H.~Voss" << Endl;
      logger << "  %``TMVA: Toolkit for multivariate data analysis,''" << Endl;
      logger << "  PoS A {\\bf CAT} (2007) 040" << Endl;
      logger << "  [arXiv:physics/0703039]." << Endl;
      logger << "  %%CITATION = POSCI,ACAT,040;%%" << Endl;
      break;

   case kHtmlLink:
      //  logger << kINFO << "  " << Endl;
      logger << kHEADER << gTools().Color("bold")
             << "Thank you for using TMVA!" << gTools().Color("reset") << Endl;
      logger << kINFO << gTools().Color("bold")
             << "For citation information, please visit: http://tmva.sf.net/citeTMVA.html"
             << gTools().Color("reset") << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TMVA::Tools::HistoHasEquidistantBins(const TH1& h)
{
   return !(h.GetXaxis()->GetXbins()->fN);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<TMatrixDSym*>*
TMVA::Tools::CalcCovarianceMatrices( const std::vector<const Event*>& events, Int_t maxCls, VariableTransformBase* transformBase )
{
   std::vector<Event*> eventVector;
   for (std::vector<const Event*>::const_iterator it = events.begin(), itEnd = events.end(); it != itEnd; ++it)
      {
         eventVector.push_back (new Event(*(*it)));
      }
   std::vector<TMatrixDSym*>* returnValue = CalcCovarianceMatrices (eventVector, maxCls, transformBase);
   for (std::vector<Event*>::const_iterator it = eventVector.begin(), itEnd = eventVector.end(); it != itEnd; ++it)
      {
         delete (*it);
      }
   return returnValue;
}

////////////////////////////////////////////////////////////////////////////////
/// compute covariance matrices

std::vector<TMatrixDSym*>*
TMVA::Tools::CalcCovarianceMatrices( const std::vector<Event*>& events, Int_t maxCls, VariableTransformBase* transformBase )
{
   if (events.empty()) {
      Log() << kWARNING << " Asked to calculate a covariance matrix for an empty event vectors.. sorry cannot do that -> return NULL"<<Endl;
      return 0;
   }

   UInt_t nvars=0, ntgts=0, nspcts=0;
   if (transformBase)
      transformBase->CountVariableTypes( nvars, ntgts, nspcts );
   else {
      nvars =events.at(0)->GetNVariables ();
      ntgts =events.at(0)->GetNTargets   ();
      nspcts=events.at(0)->GetNSpectators();
   }


   // init matrices
   Int_t matNum = maxCls;
   if (maxCls > 1 ) matNum++; // if more than one classes, then produce one matrix for all events as well (beside the matrices for each class)

   std::vector<TVectorD*>* vec = new std::vector<TVectorD*>(matNum);
   std::vector<TMatrixD*>* mat2 = new std::vector<TMatrixD*>(matNum);
   std::vector<Double_t> count(matNum);
   count.assign(matNum,0);

   Int_t cls = 0;
   TVectorD* v;
   TMatrixD* m;
   UInt_t ivar=0, jvar=0;
   for (cls = 0; cls < matNum ; cls++) {
      vec->at(cls) = new TVectorD(nvars);
      mat2->at(cls) = new TMatrixD(nvars,nvars);
      v = vec->at(cls);
      m = mat2->at(cls);

      for (ivar=0; ivar<nvars; ivar++) {
         (*v)(ivar) = 0;
         for (jvar=0; jvar<nvars; jvar++) {
            (*m)(ivar, jvar) = 0;
         }
      }
   }

   // perform event loop
   for (UInt_t i=0; i<events.size(); i++) {

      // fill the event
      const Event * ev = events[i];
      cls = ev->GetClass();
      Double_t weight = ev->GetWeight();

      std::vector<Float_t> input;
      std::vector<Char_t> mask; // entries with kTRUE must not be transformed
      // Bool_t hasMaskedEntries = kFALSE;
      if (transformBase) {
         /* hasMaskedEntries = */ transformBase->GetInput (ev, input, mask);
      } else {
         for (ivar=0; ivar<nvars; ++ivar) {
            input.push_back (ev->GetValue(ivar));
         }
      }

      if (maxCls > 1) {
         v = vec->at(matNum-1);
         m = mat2->at(matNum-1);

         count.at(matNum-1)+=weight; // count used events
         for (ivar=0; ivar<nvars; ivar++) {

            Double_t xi = input.at (ivar);
            (*v)(ivar) += xi*weight;
            (*m)(ivar, ivar) += (xi*xi*weight);

            for (jvar=ivar+1; jvar<nvars; jvar++) {
               Double_t xj = input.at (jvar);
               (*m)(ivar, jvar) += (xi*xj*weight);
               (*m)(jvar, ivar) = (*m)(ivar, jvar); // symmetric matrix
            }
         }
      }

      count.at(cls)+=weight; // count used events
      v = vec->at(cls);
      m = mat2->at(cls);
      for (ivar=0; ivar<nvars; ivar++) {
         Double_t xi = input.at (ivar);
         (*v)(ivar) += xi*weight;
         (*m)(ivar, ivar) += (xi*xi*weight);

         for (jvar=ivar+1; jvar<nvars; jvar++) {
            Double_t xj = input.at (jvar);
            (*m)(ivar, jvar) += (xi*xj*weight);
            (*m)(jvar, ivar) = (*m)(ivar, jvar); // symmetric matrix
         }
      }
   }

   // variance-covariance
   std::vector<TMatrixDSym*>* mat = new std::vector<TMatrixDSym*>(matNum);
   for (cls = 0; cls < matNum; cls++) {
      v = vec->at(cls);
      m = mat2->at(cls);

      mat->at(cls) = new TMatrixDSym(nvars);

      Double_t n = count.at(cls);
      for (ivar=0; ivar<nvars; ivar++) {
         for (jvar=0; jvar<nvars; jvar++) {
            (*(mat->at(cls)))(ivar, jvar) = (*m)(ivar, jvar)/n - (*v)(ivar)*(*v)(jvar)/(n*n);
         }
      }
      delete v;
      delete m;
   }

   delete mat2;
   delete vec;

   return mat;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the weighted mean of an array defined by the first and
/// last iterators. The w iterator should point to the first element
/// of a vector of weights of the same size as the main array.

template <typename Iterator, typename WeightIterator>
Double_t TMVA::Tools::Mean ( Iterator first,  Iterator last,  WeightIterator w)
{
   Double_t sum = 0;
   Double_t sumw = 0;
   int i = 0;
   if (w==NULL)
      {
         while ( first != last )
            {
               // if ( *w < 0) {
               //    ::Error("TMVA::Tools::Mean","w[%d] = %.4e < 0 ?!",i,*w);
               //    return 0;
               // } // SURE, why wouldn't you allow for negative event weights here ?? :)
               sum  += (*first);
               sumw += 1.0 ;
               ++first;
               ++i;
            }
         if (sumw <= 0) {
            ::Error("TMVA::Tools::Mean","sum of weights <= 0 ?! that's a bit too much of negative event weights :) ");
            return 0;
         }
      }
   else
      {
         while ( first != last )
            {
               // if ( *w < 0) {
               //    ::Error("TMVA::Tools::Mean","w[%d] = %.4e < 0 ?!",i,*w);
               //    return 0;
               // } // SURE, why wouldn't you allow for negative event weights here ?? :)
               sum  += (*w) * (*first);
               sumw += (*w) ;
               ++w;
               ++first;
               ++i;
            }
         if (sumw <= 0) {
            ::Error("TMVA::Tools::Mean","sum of weights <= 0 ?! that's a bit too much of negative event weights :) ");
            return 0;
         }
      }
   return sum/sumw;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the weighted mean of an array a with length n.

template <typename T>
Double_t TMVA::Tools::Mean(Long64_t n, const T *a, const Double_t *w)
{
   if (w) {
      return TMVA::Tools::Mean(a, a+n, w);
   } else {
      return TMath::Mean(a, a+n);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the Standard Deviation of an array defined by the iterators.
/// Note that this function returns the sigma(standard deviation) and
/// not the root mean square of the array.

template <typename Iterator, typename WeightIterator>
Double_t TMVA::Tools::RMS(Iterator first, Iterator last, WeightIterator w)
{

   Double_t sum = 0;
   Double_t sum2 = 0;
   Double_t sumw = 0;

   Double_t adouble;
   if (w==NULL)
      {
         while ( first != last ) {
            adouble=Double_t(*first);
            sum  += adouble;
            sum2 += adouble*adouble;
            sumw += 1.0;
            ++first;
         }
      }
   else
      {
         while ( first != last ) {
            adouble=Double_t(*first);
            sum  += adouble * (*w);
            sum2 += adouble*adouble * (*w);
            sumw += (*w);
            ++first;
            ++w;
         }
      }
   Double_t norm = 1./sumw;
   Double_t mean = sum*norm;
   Double_t rms = TMath::Sqrt(TMath::Abs(sum2*norm -mean*mean));
   return rms;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the Standard Deviation of an array a with length n.
/// Note that this function returns the sigma(standard deviation) and
/// not the root mean square of the array.

template <typename T>
Double_t TMVA::Tools::RMS(Long64_t n, const T *a, const Double_t *w)
{

   if (w) {
      return TMVA::Tools::RMS(a, a+n, w);
   } else {
      return TMath::RMS(a, a+n);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// get the cumulative distribution of a histogram

TH1* TMVA::Tools::GetCumulativeDist( TH1* h)
{
   TH1* cumulativeDist= (TH1*) h->Clone(Form("%sCumul",h->GetTitle()));
   //cumulativeDist->Smooth(5); // with this, I get less beautiful ROC curves, hence out!

   Float_t partialSum = 0;
   Float_t inverseSum = 0.;

   Float_t val;
   for (Int_t ibinEnd=1, ibin=cumulativeDist->GetNbinsX(); ibin >=ibinEnd ; ibin--){
      val = cumulativeDist->GetBinContent(ibin);
      if (val>0) inverseSum += val;
   }
   inverseSum = 1/inverseSum; // as I learned multiplications are much faster than division, and later I need one per bin. Well, not that it would really matter here I guess :)

   for (Int_t ibinEnd=1, ibin=cumulativeDist->GetNbinsX(); ibin >=ibinEnd ; ibin--){
      val = cumulativeDist->GetBinContent(ibin);
      if (val>0) partialSum += val;
      cumulativeDist->SetBinContent(ibin,partialSum*inverseSum);
   }
   return cumulativeDist;
}

void TMVA::Tools::ReadAttr(void *node, const char *attrname, float &value)
{
   // read attribute from xml
   const char *val = xmlengine().GetAttr(node, attrname);
   if (val == nullptr) {
      const char *nodename = xmlengine().GetNodeName(node);
      Log() << kFATAL << "Trying to read non-existing attribute '" << attrname << "' from xml node '" << nodename << "'"
            << Endl;
   } else
      value = atof(val);
}

void TMVA::Tools::ReadAttr(void *node, const char *attrname, int &value)
{
   // read attribute from xml
   const char *val = xmlengine().GetAttr(node, attrname);
   if (val == nullptr) {
      const char *nodename = xmlengine().GetNodeName(node);
      Log() << kFATAL << "Trying to read non-existing attribute '" << attrname << "' from xml node '" << nodename << "'"
            << Endl;
   } else
      value = atoi(val);
}

void TMVA::Tools::ReadAttr(void *node, const char *attrname, short &value)
{
   // read attribute from xml
   const char *val = xmlengine().GetAttr(node, attrname);
   if (val == nullptr) {
      const char *nodename = xmlengine().GetNodeName(node);
      Log() << kFATAL << "Trying to read non-existing attribute '" << attrname << "' from xml node '" << nodename << "'"
            << Endl;
   } else
      value = atoi(val);
}
