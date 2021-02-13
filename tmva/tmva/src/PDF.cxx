// @(#)root/tmva $Id$
// Author: Asen Christov, Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss, Jan Therhaag, Eckhard von Toerne

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : PDF                                                                   *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Asen Christov   <christov@physik.uni-freiburg.de> - Freiburg U., Germany  *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Jan Therhaag       <Jan.Therhaag@cern.ch>     - U of Bonn, Germany        *
 *      Eckhard von Toerne <evt@physik.uni-bonn.de>   - U of Bonn, Germany        *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      Freiburg U., Germany                                                      *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::PDF
\ingroup TMVA
PDF wrapper for histograms; uses user-defined spline interpolation.
*/

#include "TMVA/PDF.h"

#include "TMVA/Configurable.h"
#include "TMVA/KDEKernel.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Types.h"
#include "TMVA/Tools.h"
#include "TMVA/TSpline1.h"
#include "TMVA/TSpline2.h"
#include "TMVA/Version.h"

#include "TF1.h"
#include "TH1F.h"
#include "TMath.h"
#include "TVectorD.h"
#include "TGraph.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

// static configuration settings
const Int_t    TMVA::PDF::fgNbin_PdfHist      = 10000;
const Bool_t   TMVA::PDF::fgManualIntegration = kTRUE;
const Double_t TMVA::PDF::fgEpsilon           = 1.0e-12;

ClassImp(TMVA::PDF);

////////////////////////////////////////////////////////////////////////////////
/// default constructor needed for ROOT I/O

TMVA::PDF::PDF( const TString& name, Bool_t norm )
: Configurable   (""),
   fUseHistogram  ( kFALSE ),
   fPDFName       ( name ),
   fNsmooth       ( 0 ),
   fMinNsmooth    (-1 ),
   fMaxNsmooth    (-1 ),
   fNSmoothHist   ( 0 ),
   fInterpolMethod( PDF::kSpline2 ),
   fSpline        ( 0 ),
   fPDFHist       ( 0 ),
   fHist          ( 0 ),
   fHistOriginal  ( 0 ),
   fGraph         ( 0 ),
   fIGetVal       ( 0 ),
   fHistAvgEvtPerBin  ( 0 ),
   fHistDefinedNBins  ( 0 ),
   fKDEtypeString     ( 0 ),
   fKDEiterString     ( 0 ),
   fBorderMethodString( 0 ),
   fInterpolateString ( 0 ),
   fKDEtype       ( KDEKernel::kNone ),
   fKDEiter       ( KDEKernel::kNonadaptiveKDE ),
   fKDEborder     ( KDEKernel::kNoTreatment ),
   fFineFactor    ( 0. ),
   fReadingVersion( 0 ),
   fCheckHist     ( kFALSE ),
   fNormalize     ( norm ),
   fSuffix        ( "" ),
   fLogger        ( 0 )
{
   fLogger   = new MsgLogger(this);
   GetThisPdfThreadLocal() = this;
}

////////////////////////////////////////////////////////////////////////////////
/// constructor of spline based PDF:

TMVA::PDF::PDF( const TString& name,
                const TH1 *hist,
                PDF::EInterpolateMethod method,
                Int_t minnsmooth,
                Int_t maxnsmooth,
                Bool_t checkHist,
                Bool_t norm) :
   Configurable   (""),
   fUseHistogram  ( kFALSE ),
   fPDFName       ( name ),
   fMinNsmooth    ( minnsmooth ),
   fMaxNsmooth    ( maxnsmooth ),
   fNSmoothHist   ( 0 ),
   fInterpolMethod( method ),
   fSpline        ( 0 ),
   fPDFHist       ( 0 ),
   fHist          ( 0 ),
   fHistOriginal  ( 0 ),
   fGraph         ( 0 ),
   fIGetVal       ( 0 ),
   fHistAvgEvtPerBin  ( 0 ),
   fHistDefinedNBins  ( 0 ),
   fKDEtypeString     ( 0 ),
   fKDEiterString     ( 0 ),
   fBorderMethodString( 0 ),
   fInterpolateString ( 0 ),
   fKDEtype       ( KDEKernel::kNone ),
   fKDEiter       ( KDEKernel::kNonadaptiveKDE ),
   fKDEborder     ( KDEKernel::kNoTreatment ),
   fFineFactor    ( 0. ),
   fReadingVersion( 0 ),
   fCheckHist     ( checkHist ),
   fNormalize     ( norm ),
   fSuffix        ( "" ),
   fLogger        ( 0 )
{
   fLogger   = new MsgLogger(this);
   BuildPDF( hist );
}

////////////////////////////////////////////////////////////////////////////////
/// constructor of kernel based PDF:

TMVA::PDF::PDF( const TString& name,
                const TH1* hist,
                KDEKernel::EKernelType ktype,
                KDEKernel::EKernelIter kiter,
                KDEKernel::EKernelBorder kborder,
                Float_t FineFactor,
                Bool_t norm) :
   Configurable   (""),
   fUseHistogram  ( kFALSE ),
   fPDFName       ( name ),
   fNsmooth       ( 0 ),
   fMinNsmooth    (-1 ),
   fMaxNsmooth    (-1 ),
   fNSmoothHist   ( 0 ),
   fInterpolMethod( PDF::kKDE ),
   fSpline        ( 0 ),
   fPDFHist       ( 0 ),
   fHist          ( 0 ),
   fHistOriginal  ( 0 ),
   fGraph         ( 0 ),
   fIGetVal       ( 0 ),
   fHistAvgEvtPerBin  ( 0 ),
   fHistDefinedNBins  ( 0 ),
   fKDEtypeString     ( 0 ),
   fKDEiterString     ( 0 ),
   fBorderMethodString( 0 ),
   fInterpolateString ( 0 ),
   fKDEtype       ( ktype ),
   fKDEiter       ( kiter ),
   fKDEborder     ( kborder ),
   fFineFactor    ( FineFactor),
   fReadingVersion( 0 ),
   fCheckHist     ( kFALSE ),
   fNormalize     ( norm ),
   fSuffix        ( "" ),
   fLogger        ( 0 )
{
   fLogger   = new MsgLogger(this);
   BuildPDF( hist );
}

////////////////////////////////////////////////////////////////////////////////

TMVA::PDF::PDF( const TString& name,
                const TString& options,
                const TString& suffix,
                PDF* defaultPDF,
                Bool_t norm) :
   Configurable   (options),
   fUseHistogram  ( kFALSE ),
   fPDFName       ( name ),
   fNsmooth       ( 0 ),
   fMinNsmooth    ( -1 ),
   fMaxNsmooth    ( -1 ),
   fNSmoothHist   ( 0 ),
   fInterpolMethod( PDF::kSpline0 ),
   fSpline        ( 0 ),
   fPDFHist       ( 0 ),
   fHist          ( 0 ),
   fHistOriginal  ( 0 ),
   fGraph         ( 0 ),
   fIGetVal       ( 0 ),
   fHistAvgEvtPerBin  ( 50 ),
   fHistDefinedNBins  ( 0 ),
   fKDEtypeString     ( "Gauss" ),
   fKDEiterString     ( "Nonadaptive" ),
   fBorderMethodString( "None" ),
   fInterpolateString ( "Spline2" ),
   fKDEtype       ( KDEKernel::kNone ),
   fKDEiter       ( KDEKernel::kNonadaptiveKDE ),
   fKDEborder     ( KDEKernel::kNoTreatment ),
   fFineFactor    ( 1. ),
   fReadingVersion( 0 ),
   fCheckHist     ( kFALSE ),
   fNormalize     ( norm ),
   fSuffix        ( suffix ),
   fLogger        ( 0 )
{
   fLogger   = new MsgLogger(this);
   if (defaultPDF != 0) {
      fNsmooth            = defaultPDF->fNsmooth;
      fMinNsmooth         = defaultPDF->fMinNsmooth;
      fMaxNsmooth         = defaultPDF->fMaxNsmooth;
      fHistAvgEvtPerBin   = defaultPDF->fHistAvgEvtPerBin;
      fHistAvgEvtPerBin   = defaultPDF->fHistAvgEvtPerBin;
      fInterpolateString  = defaultPDF->fInterpolateString;
      fKDEtypeString      = defaultPDF->fKDEtypeString;
      fKDEiterString      = defaultPDF->fKDEiterString;
      fFineFactor         = defaultPDF->fFineFactor;
      fBorderMethodString = defaultPDF->fBorderMethodString;
      fCheckHist          = defaultPDF->fCheckHist;
      fHistDefinedNBins   = defaultPDF->fHistDefinedNBins;
   }
}

//___________________fNSmoothHist____________________________________________________
TMVA::PDF::~PDF()
{
   // destructor
   if (fSpline       != NULL) delete fSpline;
   if (fHist         != NULL) delete fHist;
   if (fPDFHist      != NULL) delete fPDFHist;
   if (fHistOriginal != NULL) delete fHistOriginal;
   if (fIGetVal      != NULL) delete fIGetVal;
   if (fGraph        != NULL) delete fGraph;
   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::PDF::BuildPDF( const TH1* hist )
{
   GetThisPdfThreadLocal() = this;
   // sanity check
   if (hist == NULL) Log() << kFATAL << "Called without valid histogram pointer!" << Endl;

   // histogram should be non empty
   if (hist->GetEntries() <= 0)
      Log() << kFATAL << "Number of entries <= 0 (" << hist->GetEntries() << " in histogram: " << hist->GetTitle() << ")" << Endl;

   if (fInterpolMethod == PDF::kKDE) {
     Log()<< kDEBUG << "Create "
            << ((fKDEiter == KDEKernel::kNonadaptiveKDE) ? "nonadaptive " :
                (fKDEiter == KDEKernel::kAdaptiveKDE)    ? "adaptive "    : "??? ")
            << ((fKDEtype == KDEKernel::kGauss)          ? "Gauss "       : "??? ")
            << "type KDE kernel for histogram: \"" << hist->GetName() << "\""
            << Endl;
   }
   else {
      // another sanity check (nsmooth<0 indicated build with KDE)
      if (fMinNsmooth<0)
         Log() << kFATAL << "PDF construction called with minnsmooth<0" << Endl;
      else if (fMaxNsmooth<=0)
         fMaxNsmooth = fMinNsmooth;
      else if (fMaxNsmooth<fMinNsmooth)
         Log() << kFATAL << "PDF construction called with maxnsmooth<minnsmooth" << Endl;
   }

   fHistOriginal = (TH1F*)hist->Clone( TString(hist->GetName()) + "_original" );
   fHist         = (TH1F*)hist->Clone( TString(hist->GetName()) + "_smoothed" );
   fHistOriginal->SetTitle( fHistOriginal->GetName() ); // reset to new title as well
   fHist        ->SetTitle( fHist->GetName() );

   // do not store in current target file
   fHistOriginal->SetDirectory(0);
   fHist        ->SetDirectory(0);
   fUseHistogram = kFALSE;

   if (fInterpolMethod == PDF::kKDE) BuildKDEPDF();
   else                              BuildSplinePDF();
}

////////////////////////////////////////////////////////////////////////////////

Int_t TMVA::PDF::GetHistNBins ( Int_t evtNum )
{
   Int_t ResolutionFactor = (fInterpolMethod == PDF::kKDE) ? 5 : 1;
   if (evtNum == 0 && fHistDefinedNBins == 0)
      Log() << kFATAL << "No number of bins set for PDF" << Endl;
   else if (fHistDefinedNBins > 0)
      return fHistDefinedNBins * ResolutionFactor;
   else if ( evtNum > 0 && fHistAvgEvtPerBin > 0 )
      return evtNum / fHistAvgEvtPerBin * ResolutionFactor;
   else
      Log() << kFATAL << "No number of bins or average event per bin set for PDF" << fHistAvgEvtPerBin << Endl;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// build the PDF from the original histograms

void TMVA::PDF::BuildSplinePDF()
{
   // (not useful for discrete distributions, or if no splines are requested)
   if (fInterpolMethod != PDF::kSpline0 && fCheckHist) CheckHist();
   // use ROOT TH1 smooth method
   fNSmoothHist = 0;
   if (fMaxNsmooth > 0 && fMinNsmooth >=0 ) SmoothHistogram();

   // fill histogramm to graph
   FillHistToGraph();

   //   fGraph->Print();
   switch (fInterpolMethod) {

   case kSpline0:
      // use original histogram as reference
      // this is useful, eg, for discrete variables
      fUseHistogram = kTRUE;
      break;

   case kSpline1:
      fSpline = new TMVA::TSpline1( "spline1", new TGraph(*fGraph) );
      break;

   case kSpline2:
      fSpline = new TMVA::TSpline2( "spline2", new TGraph(*fGraph) );
      break;

   case kSpline3:
      fSpline = new TSpline3( "spline3", new TGraph(*fGraph) );
      break;

   case kSpline5:
      fSpline = new TSpline5( "spline5", new TGraph(*fGraph) );
      break;

   default:
      Log() << kWARNING << "No valid interpolation method given! Use Spline2" << Endl;
      fSpline = new TMVA::TSpline2( "spline2", new TGraph(*fGraph) );
      Log() << kFATAL << " Well.. .thinking about it, I better quit so you notice you are forced to fix the mistake " << Endl;
      std::exit(1);
   }
   // fill into histogram
   FillSplineToHist();

   if (!UseHistogram()) {
      fSpline->SetTitle( (TString)fHist->GetTitle() + fSpline->GetTitle() );
      fSpline->SetName ( (TString)fHist->GetName()  + fSpline->GetName()  );
   }


   // sanity check
   Double_t integral = GetIntegral();
   if (integral < 0) Log() << kFATAL << "Integral: " << integral << " <= 0" << Endl;

   // normalize
   if (fNormalize)
      if (integral>0) fPDFHist->Scale( 1.0/integral );

   fPDFHist->SetDirectory(0);
}

////////////////////////////////////////////////////////////////////////////////
/// creates high-binned reference histogram to be used instead of the
/// PDF for speed reasons

void TMVA::PDF::BuildKDEPDF()
{
   fPDFHist = new TH1F( "", "", fgNbin_PdfHist, GetXmin(), GetXmax() );
   fPDFHist->SetTitle( (TString)fHist->GetTitle() + "_hist from_KDE" );
   fPDFHist->SetName ( (TString)fHist->GetName()  + "_hist_from_KDE" );

   // create the kernel object
   Float_t histoLowEdge   = fHist->GetBinLowEdge(1);
   Float_t histoUpperEdge = fPDFHist->GetBinLowEdge(fPDFHist->GetNbinsX()) + fPDFHist->GetBinWidth(fPDFHist->GetNbinsX());
   KDEKernel *kern = new TMVA::KDEKernel( fKDEiter,
                                          fHist, histoLowEdge, histoUpperEdge,
                                          fKDEborder, fFineFactor );
   kern->SetKernelType(fKDEtype);

   for (Int_t i=1;i<fHist->GetNbinsX();i++) {
      // loop over the bins of the original histo
      for (Int_t j=1;j<fPDFHist->GetNbinsX();j++) {
         // loop over the bins of the new histo and fill it
         fPDFHist->AddBinContent(j,fHist->GetBinContent(i)*
                                 kern->GetBinKernelIntegral(fPDFHist->GetBinLowEdge(j),
                                                            fPDFHist->GetBinLowEdge(j+1),
                                                            fHist->GetBinCenter(i),
                                                            i)
                                 );
      }
      if (fKDEborder == 3) { // mirror the samples and fill them again
         // in order to save time do the mirroring only for the first (the lower) 1/5 of the histo to the left;
         // and the last (the higher) 1/5 of the histo to the right.
         // the middle part of the histo, which is not mirrored, has no influence on the border effects anyway ...
         if (i < fHist->GetNbinsX()/5  ) {  // the first (the lower) 1/5 of the histo
            for (Int_t j=1;j<fPDFHist->GetNbinsX();j++) {
               // loop over the bins of the PDF histo and fill it
               fPDFHist->AddBinContent(j,fHist->GetBinContent(i)*
                                       kern->GetBinKernelIntegral(fPDFHist->GetBinLowEdge(j),
                                                                  fPDFHist->GetBinLowEdge(j+1),
                                                                  2*histoLowEdge-fHist->GetBinCenter(i), //  mirroring to the left
                                                                  i)
                                       );
            }
         }
         if (i > 4*fHist->GetNbinsX()/5) { // the last (the higher) 1/5 of the histo
            for (Int_t j=1;j<fPDFHist->GetNbinsX();j++) {
               // loop over the bins of the PDF histo and fill it
               fPDFHist->AddBinContent( j,fHist->GetBinContent(i)*
                                        kern->GetBinKernelIntegral(fPDFHist->GetBinLowEdge(j),
                                                                   fPDFHist->GetBinLowEdge(j+1),
                                                                   2*histoUpperEdge-fHist->GetBinCenter(i), i) );
            }
         }
      }
   }

   fPDFHist->SetEntries(fHist->GetEntries());

   delete kern;

   // sanity check
   Double_t integral = GetIntegral();
   if (integral < 0) Log() << kFATAL << "Integral: " << integral << " <= 0" << Endl;

   // normalize
   if (fNormalize)
      if (integral>0) fPDFHist->Scale( 1.0/integral );
   fPDFHist->SetDirectory(0);
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::PDF::SmoothHistogram()
{
   if(fHist->GetNbinsX()==1) return;
   if (fMaxNsmooth == fMinNsmooth) {
      fHist->Smooth( fMinNsmooth );
      return;
   }

   //calculating Mean, RMS of the relative errors and using them to set
   //the boundaries of the linear transformation
   Float_t Err=0, ErrAvg=0, ErrRMS=0 ; Int_t num=0, smooth;
   for (Int_t bin=0; bin<fHist->GetNbinsX(); bin++) {
      if (fHist->GetBinContent(bin+1) <= fHist->GetBinError(bin+1)) continue;
      Err = fHist->GetBinError(bin+1) / fHist->GetBinContent(bin+1);
      ErrAvg += Err; ErrRMS += Err*Err; num++;
   }
   ErrAvg /= num;
   ErrRMS = TMath::Sqrt(ErrRMS/num-ErrAvg*ErrAvg) ;

   //linearly convent the histogram to a vector of smoothnum
   Float_t MaxErr=ErrAvg+ErrRMS, MinErr=ErrAvg-ErrRMS;
   fNSmoothHist = new TH1I("","",fHist->GetNbinsX(),0,fHist->GetNbinsX());
   fNSmoothHist->SetTitle( (TString)fHist->GetTitle() + "_Nsmooth" );
   fNSmoothHist->SetName ( (TString)fHist->GetName()  + "_Nsmooth" );
   for (Int_t bin=0; bin<fHist->GetNbinsX(); bin++) {
      if (fHist->GetBinContent(bin+1) <= fHist->GetBinError(bin+1))
         smooth=fMaxNsmooth;
      else {
         Err = fHist->GetBinError(bin+1) / fHist->GetBinContent(bin+1);
         smooth=(Int_t)((Err-MinErr) /(MaxErr-MinErr) * (fMaxNsmooth-fMinNsmooth)) + fMinNsmooth;
      }
      smooth = TMath::Max(fMinNsmooth,TMath::Min(fMaxNsmooth,smooth));
      fNSmoothHist->SetBinContent(bin+1,smooth);
   }

   //find regions of constant smoothnum, starting from the highest amount of
   //smoothing. So the last iteration all the histogram will be smoothed as a whole
   for (Int_t n=fMaxNsmooth; n>=0; n--) {
      //all the histogram has to be smoothed at least fMinNsmooth
      if (n <= fMinNsmooth) { fHist->Smooth(); continue; }
      Int_t MinBin=-1,MaxBin =-1;
      for (Int_t bin=0; bin < fHist->GetNbinsX(); bin++) {
         if (fNSmoothHist->GetBinContent(bin+1) >= n) {
            if (MinBin==-1) MinBin = bin;
            else MaxBin=bin;
         }
         else if (MaxBin >= 0) {
            fHist->Smooth(1,"R");
            MinBin=MaxBin=-1;
         }
         else     // can't smooth a single bin
            MinBin = -1;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Simple conversion

void TMVA::PDF::FillHistToGraph()
{
   fGraph=new TGraph(fHist);
}

////////////////////////////////////////////////////////////////////////////////
/// creates high-binned reference histogram to be used instead of the
/// PDF for speed reasons

void TMVA::PDF::FillSplineToHist()
{
   if (UseHistogram()) {
      // no spline given, use the original histogram
      fPDFHist = (TH1*)fHist->Clone();
      fPDFHist->SetTitle( (TString)fHist->GetTitle() + "_hist from_spline0" );
      fPDFHist->SetName ( (TString)fHist->GetName()  + "_hist_from_spline0" );
   }
   else {
      // create new reference histogram
      fPDFHist = new TH1F( "", "", fgNbin_PdfHist, GetXmin(), GetXmax() );
      fPDFHist->SetTitle( (TString)fHist->GetTitle() + "_hist from_" + fSpline->GetTitle() );
      fPDFHist->SetName ( (TString)fHist->GetName()  + "_hist_from_" + fSpline->GetTitle() );

      for (Int_t bin=1; bin <= fgNbin_PdfHist; bin++) {
         Double_t x = fPDFHist->GetBinCenter( bin );
         Double_t y = fSpline->Eval( x );
         // sanity correction: in cases where strong slopes exist, accidentally, the
         // splines can go to zero; in this case we set the corresponding bin content
         // equal to the bin content of the original histogram
         if (y <= fgEpsilon) y = fHist->GetBinContent( fHist->FindBin( x ) );
         fPDFHist->SetBinContent( bin, TMath::Max(y, fgEpsilon) );
      }
   }
   fPDFHist->SetDirectory(0);
}

////////////////////////////////////////////////////////////////////////////////
/// sanity check: compare PDF with original histogram

void TMVA::PDF::CheckHist() const
{
   if (fHist == NULL) {
      Log() << kFATAL << "<CheckHist> Called without valid histogram pointer!" << Endl;
   }

   Int_t nbins = fHist->GetNbinsX();

   Int_t emptyBins=0;
   // count number of empty bins
   for (Int_t bin=1; bin<=nbins; bin++)
      if (fHist->GetBinContent(bin) == 0) emptyBins++;

   if (((Float_t)emptyBins/(Float_t)nbins) > 0.5) {
      Log() << kWARNING << "More than 50% (" << (((Float_t)emptyBins/(Float_t)nbins)*100)
            <<"%) of the bins in hist '"
            << fHist->GetName() << "' are empty!" << Endl;
      Log() << kWARNING << "X_min=" << GetXmin()
            << " mean=" << fHist->GetMean() << " X_max= " << GetXmax() << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// comparison of original histogram with reference PDF

void TMVA::PDF::ValidatePDF( TH1* originalHist ) const
{
   // if no histogram is given, use the original one (the one the PDF was made from)
   if (!originalHist) originalHist = fHistOriginal;

   Int_t    nbins = originalHist->GetNbinsX();

   // treat errors properly
   if (originalHist->GetSumw2()->GetSize() == 0) originalHist->Sumw2();

   // ---- first validation: simple(st) possible chi2 test
   // count number of empty bins
   Double_t chi2 = 0;
   Int_t    ndof = 0;
   Int_t    nc1  = 0; // deviation counters
   Int_t    nc2  = 0; // deviation counters
   Int_t    nc3  = 0; // deviation counters
   Int_t    nc6  = 0; // deviation counters
   for (Int_t bin=1; bin<=nbins; bin++) {
      Double_t x  = originalHist->GetBinCenter( bin );
      Double_t y  = originalHist->GetBinContent( bin );
      Double_t ey = originalHist->GetBinError( bin );

      Int_t binPdfHist = fPDFHist->FindBin( x );
      if (binPdfHist<0) continue; // happens only if hist-dim>3

      Double_t yref = GetVal( x );
      Double_t rref = ( originalHist->GetSumOfWeights()/fPDFHist->GetSumOfWeights() *
                        originalHist->GetBinWidth( bin )/fPDFHist->GetBinWidth( binPdfHist ) );

      if (y > 0) {
         ndof++;
         Double_t d = TMath::Abs( (y - yref*rref)/ey );
         //         std::cout << "bin: " << bin << "  val: " << x << "  data(err): " << y << "(" << ey << ")   pdf: "
         //              << yref << "  dev(chi2): " << d << "(" << chi2 << ")  rref: " << rref << std::endl;
         chi2 += d*d;
         if (d > 1) { nc1++; if (d > 2) { nc2++; if (d > 3) { nc3++; if (d > 6) nc6++; } } }
      }
   }

   Log() << kDEBUG << "Validation result for PDF \"" << originalHist->GetTitle() << "\"" << ": " << Endl;
   Log() << kDEBUG << Form( "    chi2/ndof(!=0) = %.1f/%i = %.2f (Prob = %.2f)",
                  chi2, ndof, chi2/ndof, TMath::Prob( chi2, ndof ) ) << Endl;
   if ((1.0 - TMath::Prob( chi2, ndof )) > 0.9999994) {
      Log() << kDEBUG << "Comparison of the original histogram \"" << originalHist->GetTitle() << "\"" << Endl;
      Log() << kDEBUG << "with the corresponding PDF gave a chi2/ndof of " << chi2/ndof << "," << Endl;
      Log() << kDEBUG << "which corresponds to a deviation of more than 5 sigma! Please check!" << Endl;
   }
   Log() << kDEBUG << Form( "    #bins-found(#expected-bins) deviating > [1,2,3,6] sigmas: " \
                  "[%i(%i),%i(%i),%i(%i),%i(%i)]",
                  nc1, Int_t(TMath::Prob(1.0,1)*ndof), nc2, Int_t(TMath::Prob(4.0,1)*ndof),
                  nc3, Int_t(TMath::Prob(9.0,1)*ndof), nc6, Int_t(TMath::Prob(36.0,1)*ndof) ) << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// computes normalisation

Double_t TMVA::PDF::GetIntegral() const
{
   Double_t integral = fPDFHist->GetSumOfWeights();
   integral *= GetPdfHistBinWidth();

   return integral;
}

////////////////////////////////////////////////////////////////////////////////
/// static external auxiliary function (integrand)

Double_t TMVA::PDF::IGetVal( Double_t* x, Double_t* )
{
   return ThisPDF()->GetVal( x[0] );
}

////////////////////////////////////////////////////////////////////////////////
/// computes PDF integral within given ranges

Double_t TMVA::PDF::GetIntegral( Double_t xmin, Double_t xmax )
{
   Double_t  integral = 0;

   if (fgManualIntegration) {

      // compute integral by summing over bins
      Int_t imin = fPDFHist->FindBin(xmin);
      Int_t imax = fPDFHist->FindBin(xmax);
      if (imin < 1)                     imin = 1;
      if (imax > fPDFHist->GetNbinsX()) imax = fPDFHist->GetNbinsX();

      for (Int_t bini = imin; bini <= imax; bini++) {
         Float_t dx = fPDFHist->GetBinWidth(bini);
         // correct for bin fractions in first and last bin
         if      (bini == imin) dx = fPDFHist->GetBinLowEdge(bini+1) - xmin;
         else if (bini == imax) dx = xmax - fPDFHist->GetBinLowEdge(bini);
         if (dx < 0 && dx > -1.0e-8) dx = 0; // protect against float->double numerical effects
         if (dx<0) {
            Log() << kWARNING
                  << "dx   = " << dx   << std::endl
                  << "bini = " << bini << std::endl
                  << "xmin = " << xmin << std::endl
                  << "xmax = " << xmax << std::endl
                  << "imin = " << imin << std::endl
                  << "imax = " << imax << std::endl
                  << "low edge of imin" << fPDFHist->GetBinLowEdge(imin) << std::endl
                  << "low edge of imin+1" << fPDFHist->GetBinLowEdge(imin+1) << Endl;
            Log() << kFATAL << "<GetIntegral> dx = " << dx << " < 0" << Endl;
         }
         integral += fPDFHist->GetBinContent(bini)*dx;
      }

   }
   else {

      // compute via Gaussian quadrature (C++ version of CERNLIB function DGAUSS)
      if (fIGetVal == 0) fIGetVal = new TF1( "IGetVal", PDF::IGetVal, GetXmin(), GetXmax(), 0 );
      integral = fIGetVal->Integral( xmin, xmax );
   }

   return integral;
}

////////////////////////////////////////////////////////////////////////////////
/// returns value PDF(x)

Double_t TMVA::PDF::GetVal( Double_t x ) const
{
   // check which is filled
   Int_t bin = fPDFHist->FindBin(x);
   bin = TMath::Max(bin,1);
   bin = TMath::Min(bin,fPDFHist->GetNbinsX());

   Double_t retval = 0;

   if (UseHistogram()) {
      // use directly histogram bins (this is for discrete PDFs)
      retval = fPDFHist->GetBinContent( bin );
   }
   else {
      // linear interpolation
      Int_t nextbin = bin;
      if ((x > fPDFHist->GetBinCenter(bin) && bin != fPDFHist->GetNbinsX()) || bin == 1)
         nextbin++;
      else
         nextbin--;

      // linear interpolation between adjacent bins
      Double_t dx = fPDFHist->GetBinCenter( bin )  - fPDFHist->GetBinCenter( nextbin );
      Double_t dy = fPDFHist->GetBinContent( bin ) - fPDFHist->GetBinContent( nextbin );
      retval = fPDFHist->GetBinContent( bin ) + (x - fPDFHist->GetBinCenter( bin ))*dy/dx;
   }

   return TMath::Max( retval, fgEpsilon );
}

////////////////////////////////////////////////////////////////////////////////
/// returns value \f$ PDF^{-1}(y) \f$

Double_t TMVA::PDF::GetValInverse( Double_t y, Bool_t isMonotonouslyIncreasingFunction ) const
{
   Int_t    lowerBin=0,      higherBin=0;
   Double_t lowerBinValue=0, higherBinValue=0;
   FindBinInverse(fPDFHist,lowerBin,higherBin,lowerBinValue,higherBinValue,y,isMonotonouslyIncreasingFunction);

   Double_t xValueLowerBin =fPDFHist->GetBinCenter (lowerBin);
   Double_t xValueHigherBin=fPDFHist->GetBinCenter (higherBin);

   Double_t length  =(higherBinValue-lowerBinValue);
   Double_t fraction=lowerBinValue;
   if (length>1.0e-10)
      fraction=(y-lowerBinValue)/length;

   Double_t lengthX =xValueHigherBin-xValueLowerBin;
   Double_t x       =xValueLowerBin+lengthX*fraction;

   // comparison for test purposes
   //   std::cout << "lb " << lowerBin << "  hb " << higherBin << "  lbv " << lowerBinValue << "  hbv " << higherBinValue << "  frac " << fraction << std::endl;
   //   std::cout << "y " << y << "  inv x " << x << "  straight y " << GetVal(x) << std::endl;

   return x;
}

////////////////////////////////////////////////////////////////////////////////
/// find bin from value on ordinate

void TMVA::PDF::FindBinInverse( const TH1* histogram, Int_t& lowerBin, Int_t& higherBin, Double_t& lowerBinValue, Double_t& higherBinValue,
                                Double_t y, Bool_t isMonotonouslyIncreasingFunction ) const
{
   if (isMonotonouslyIncreasingFunction) {
      higherBin=histogram->GetNbinsX();
      lowerBin =0;

      Int_t bin=higherBin/2;

      while (bin>lowerBin && bin<higherBin) {
         Double_t binContent=histogram->GetBinContent(bin);

         if (y<binContent) {
            higherBin     =bin;
            higherBinValue=binContent;
         }
         else if (y>=binContent){
            lowerBin      =bin;
            lowerBinValue =binContent;
         }
         bin=lowerBin+(higherBin-lowerBin)/2;
      }
      return;
   }
   // search sequentially
   for (Int_t bin=0, binEnd=histogram->GetNbinsX(); bin<binEnd; ++bin) {
      Double_t binContent=histogram->GetBinContent(bin);
      if (binContent<y) {
         lowerBin =bin;
         higherBin=bin;
         lowerBinValue =binContent;
         higherBinValue=binContent;
      }
      else {
         higherBin=bin;
         higherBinValue=binContent;
         break;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// define the options (their key words) that can be set in the option string
///
/// know options:
///
/// PDFInterpol[ivar] <string>   Spline0, Spline1, Spline2 <default>, Spline3, Spline5, KDE  used to interpolate reference histograms
///             if no variable index is given, it is valid for ALL the variables
///
///  - NSmooth           <int>    how often the input histos are smoothed
///  - MinNSmooth        <int>    min number of smoothing iterations, for bins with most data
///  - MaxNSmooth        <int>    max number of smoothing iterations, for bins with least data
///  - NAvEvtPerBin      <int>    minimum average number of events per PDF bin
///  - TransformOutput   <bool>   transform (often strongly peaked) likelihood output through sigmoid inversion
///  - fKDEtype          <KernelType>   type of the Kernel to use (1 is Gaussian)
///  - fKDEiter          <KerneIter>    number of iterations (1 --> "static KDE", 2 --> "adaptive KDE")
///  - fBorderMethod     <KernelBorder> the method to take care about "border" effects (1=no treatment , 2=kernel renormalization, 3=sample mirroring)

void TMVA::PDF::DeclareOptions()
{
   DeclareOptionRef( fNsmooth, Form("NSmooth%s",fSuffix.Data()),
                     "Number of smoothing iterations for the input histograms" );
   DeclareOptionRef( fMinNsmooth, Form("MinNSmooth%s",fSuffix.Data()),
                     "Min number of smoothing iterations, for bins with most data" );

   DeclareOptionRef( fMaxNsmooth, Form("MaxNSmooth%s",fSuffix.Data()),
                     "Max number of smoothing iterations, for bins with least data" );

   DeclareOptionRef( fHistAvgEvtPerBin, Form("NAvEvtPerBin%s",fSuffix.Data()),
                     "Average number of events per PDF bin" );

   DeclareOptionRef( fHistDefinedNBins, Form("Nbins%s",fSuffix.Data()),
                     "Defined number of bins for the histogram from which the PDF is created" );

   DeclareOptionRef( fCheckHist, Form("CheckHist%s",fSuffix.Data()),
                     "Whether or not to check the source histogram of the PDF" );

   DeclareOptionRef( fInterpolateString, Form("PDFInterpol%s",fSuffix.Data()),
                     "Interpolation method for reference histograms (e.g. Spline2 or KDE)" );
   AddPreDefVal(TString("Spline0")); // take histogram
   AddPreDefVal(TString("Spline1")); // linear interpolation between bins
   AddPreDefVal(TString("Spline2")); // quadratic interpolation
   AddPreDefVal(TString("Spline3")); // cubic interpolation
   AddPreDefVal(TString("Spline5")); // fifth order polynom interpolation
   AddPreDefVal(TString("KDE"));     // use kernel density estimator

   DeclareOptionRef( fKDEtypeString, Form("KDEtype%s",fSuffix.Data()), "KDE kernel type (1=Gauss)" );
   AddPreDefVal(TString("Gauss"));

   DeclareOptionRef( fKDEiterString, Form("KDEiter%s",fSuffix.Data()), "Number of iterations (1=non-adaptive, 2=adaptive)" );
   AddPreDefVal(TString("Nonadaptive"));
   AddPreDefVal(TString("Adaptive"));

   DeclareOptionRef( fFineFactor , Form("KDEFineFactor%s",fSuffix.Data()),
                     "Fine tuning factor for Adaptive KDE: Factor to multiply the width of the kernel");

   DeclareOptionRef( fBorderMethodString, Form("KDEborder%s",fSuffix.Data()),
                     "Border effects treatment (1=no treatment , 2=kernel renormalization, 3=sample mirroring)" );
   AddPreDefVal(TString("None"));
   AddPreDefVal(TString("Renorm"));
   AddPreDefVal(TString("Mirror"));

   SetConfigName( GetName() );
   SetConfigDescription( "Configuration options for the PDF class" );
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::PDF::ProcessOptions()
{
   if (fNsmooth < 0) fNsmooth = 0; // set back to default

   if (fMaxNsmooth < 0 || fMinNsmooth < 0) { // use "Nsmooth" variable
      fMinNsmooth = fMaxNsmooth = fNsmooth;
   }

   if (fMaxNsmooth < fMinNsmooth && fMinNsmooth >= 0) { // sanity check
      Log() << kFATAL << "ERROR: MaxNsmooth = "
            << fMaxNsmooth << " < MinNsmooth = " << fMinNsmooth << Endl;
   }

   if (fMaxNsmooth < 0 || fMinNsmooth < 0) {
      Log() << kFATAL << "ERROR: MaxNsmooth = "
            << fMaxNsmooth << " or MinNsmooth = " << fMinNsmooth << " smaller than zero" << Endl;
   }

   //processing the options
   if      (fInterpolateString == "Spline0") fInterpolMethod = TMVA::PDF::kSpline0;
   else if (fInterpolateString == "Spline1") fInterpolMethod = TMVA::PDF::kSpline1;
   else if (fInterpolateString == "Spline2") fInterpolMethod = TMVA::PDF::kSpline2;
   else if (fInterpolateString == "Spline3") fInterpolMethod = TMVA::PDF::kSpline3;
   else if (fInterpolateString == "Spline5") fInterpolMethod = TMVA::PDF::kSpline5;
   else if (fInterpolateString == "KDE"    ) fInterpolMethod = TMVA::PDF::kKDE;
   else if (fInterpolateString != ""       ) {
      Log() << kFATAL << "unknown setting for option 'InterpolateMethod': " << fKDEtypeString << ((fSuffix=="")?"":Form(" for pdf with suffix %s",fSuffix.Data())) << Endl;
   }

   // init KDE options
   if      (fKDEtypeString == "Gauss"      ) fKDEtype = KDEKernel::kGauss;
   else if (fKDEtypeString != ""           )
      Log() << kFATAL << "unknown setting for option 'KDEtype': " << fKDEtypeString << ((fSuffix=="")?"":Form(" for pdf with suffix %s",fSuffix.Data())) << Endl;
   if      (fKDEiterString == "Nonadaptive") fKDEiter = KDEKernel::kNonadaptiveKDE;
   else if (fKDEiterString == "Adaptive"   ) fKDEiter = KDEKernel::kAdaptiveKDE;
   else if (fKDEiterString != ""           )// nothing more known
      Log() << kFATAL << "unknown setting for option 'KDEiter': " << fKDEtypeString << ((fSuffix=="")?"":Form(" for pdf with suffix %s",fSuffix.Data())) << Endl;

   if       ( fBorderMethodString == "None"   ) fKDEborder= KDEKernel::kNoTreatment;
   else if  ( fBorderMethodString == "Renorm" ) fKDEborder= KDEKernel::kKernelRenorm;
   else if  ( fBorderMethodString == "Mirror" ) fKDEborder= KDEKernel::kSampleMirror;
   else if  ( fKDEiterString != ""            ) { // nothing more known
      Log() << kFATAL << "unknown setting for option 'KDEBorder': " << fKDEtypeString << ((fSuffix=="")?"":Form(" for pdf with suffix %s",fSuffix.Data())) << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// XML file writing

void TMVA::PDF::AddXMLTo( void* parent )
{
   void* pdfxml = gTools().AddChild(parent, "PDF");
   gTools().AddAttr(pdfxml, "Name",           fPDFName );
   gTools().AddAttr(pdfxml, "MinNSmooth",     fMinNsmooth );
   gTools().AddAttr(pdfxml, "MaxNSmooth",     fMaxNsmooth );
   gTools().AddAttr(pdfxml, "InterpolMethod", fInterpolMethod );
   gTools().AddAttr(pdfxml, "KDE_type",       fKDEtype );
   gTools().AddAttr(pdfxml, "KDE_iter",       fKDEiter );
   gTools().AddAttr(pdfxml, "KDE_border",     fKDEborder );
   gTools().AddAttr(pdfxml, "KDE_finefactor", fFineFactor );
   void* pdfhist = gTools().AddChild(pdfxml,"Histogram" );
   TH1*  histToWrite = GetOriginalHist();
   Bool_t hasEquidistantBinning = gTools().HistoHasEquidistantBins(*histToWrite);
   gTools().AddAttr(pdfhist, "Name",  histToWrite->GetName() );
   gTools().AddAttr(pdfhist, "NBins", histToWrite->GetNbinsX() );
   gTools().AddAttr(pdfhist, "XMin",  histToWrite->GetXaxis()->GetXmin() );
   gTools().AddAttr(pdfhist, "XMax",  histToWrite->GetXaxis()->GetXmax() );
   gTools().AddAttr(pdfhist, "HasEquidistantBins", hasEquidistantBinning );

   TString bincontent("");
   for (Int_t i=0; i<histToWrite->GetNbinsX(); i++) {
      bincontent += gTools().StringFromDouble(histToWrite->GetBinContent(i+1));
      bincontent += " ";
   }
   gTools().AddRawLine(pdfhist, bincontent );

   if (!hasEquidistantBinning) {
      void* pdfhistbins = gTools().AddChild(pdfxml,"HistogramBinning" );
      gTools().AddAttr(pdfhistbins, "NBins", histToWrite->GetNbinsX() );
      TString binns("");
      for (Int_t i=1; i<=histToWrite->GetNbinsX()+1; i++) {
         binns += gTools().StringFromDouble(histToWrite->GetXaxis()->GetBinLowEdge(i));
         binns += " ";
      }
      gTools().AddRawLine(pdfhistbins, binns );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// XML file reading

void TMVA::PDF::ReadXML( void* pdfnode )
{
   UInt_t enumint;

   gTools().ReadAttr(pdfnode, "MinNSmooth",     fMinNsmooth );
   gTools().ReadAttr(pdfnode, "MaxNSmooth",     fMaxNsmooth );
   gTools().ReadAttr(pdfnode, "InterpolMethod", enumint ); fInterpolMethod = EInterpolateMethod(enumint);
   gTools().ReadAttr(pdfnode, "KDE_type",       enumint ); fKDEtype        = KDEKernel::EKernelType(enumint);
   gTools().ReadAttr(pdfnode, "KDE_iter",       enumint ); fKDEiter        = KDEKernel::EKernelIter(enumint);
   gTools().ReadAttr(pdfnode, "KDE_border",     enumint ); fKDEborder      = KDEKernel::EKernelBorder(enumint);
   gTools().ReadAttr(pdfnode, "KDE_finefactor", fFineFactor );

   TString hname;
   UInt_t nbins;
   Double_t xmin, xmax;
   Bool_t hasEquidistantBinning;

   void* histch = gTools().GetChild(pdfnode);
   gTools().ReadAttr( histch, "Name",  hname );
   gTools().ReadAttr( histch, "NBins", nbins );
   gTools().ReadAttr( histch, "XMin",  xmin );
   gTools().ReadAttr( histch, "XMax",  xmax );
   gTools().ReadAttr( histch, "HasEquidistantBins", hasEquidistantBinning );

   // recreate the original hist
   TH1* newhist = 0;
   if (hasEquidistantBinning) {
      newhist = new TH1F( hname, hname, nbins, xmin, xmax );
      newhist->SetDirectory(0);
      const char* content = gTools().GetContent(histch);
      std::stringstream s(content);
      Double_t val;
      for (UInt_t i=0; i<nbins; i++) {
         s >> val;
         newhist->SetBinContent(i+1,val);
      }
   }
   else {
      const char* content = gTools().GetContent(histch);
      std::stringstream s(content);
      Double_t val;
      void* binch = gTools().GetNextChild(histch);
      UInt_t nbinning;
      gTools().ReadAttr( binch, "NBins", nbinning );
      TVectorD binns(nbinning+1);
      if (nbinning != nbins) {
         Log() << kFATAL << "Number of bins in content and binning array differs"<<Endl;
      }
      const char* binString = gTools().GetContent(binch);
      std::stringstream sb(binString);
      for (UInt_t i=0; i<=nbins; i++) sb >> binns[i];
      newhist =  new TH1F( hname, hname, nbins, binns.GetMatrixArray() );
      newhist->SetDirectory(0);
      for (UInt_t i=0; i<nbins; i++) {
         s >> val;
         newhist->SetBinContent(i+1,val);
      }
   }

   TString hnameSmooth = hname;
   hnameSmooth.ReplaceAll( "_original", "_smoothed" );

   if (fHistOriginal != 0) delete fHistOriginal;
   fHistOriginal = newhist;
   fHist = (TH1F*)fHistOriginal->Clone( hnameSmooth );
   fHist->SetTitle( hnameSmooth );
   fHist->SetDirectory(0);

   if (fInterpolMethod == PDF::kKDE) BuildKDEPDF();
   else                              BuildSplinePDF();
}

////////////////////////////////////////////////////////////////////////////////
/// write the pdf

std::ostream& TMVA::operator<< ( std::ostream& os, const PDF& pdf )
{
   Int_t dp = os.precision();
   os << "MinNSmooth      " << pdf.fMinNsmooth << std::endl;
   os << "MaxNSmooth      " << pdf.fMaxNsmooth << std::endl;
   os << "InterpolMethod  " << pdf.fInterpolMethod << std::endl;
   os << "KDE_type        " << pdf.fKDEtype << std::endl;
   os << "KDE_iter        " << pdf.fKDEiter << std::endl;
   os << "KDE_border      " << pdf.fKDEborder << std::endl;
   os << "KDE_finefactor  " << pdf.fFineFactor << std::endl;

   TH1* histToWrite = pdf.GetOriginalHist();

   const Int_t nBins = histToWrite->GetNbinsX();

   // note: this is a schema change introduced for v3.7.3
   os << "Histogram       "
      << histToWrite->GetName()
      << "   " << nBins                                 // nbins
      << "   " << std::setprecision(12) << histToWrite->GetXaxis()->GetXmin()    // x_min
      << "   " << std::setprecision(12) << histToWrite->GetXaxis()->GetXmax()    // x_max
      << std::endl;

   // write the smoothed hist
   os << "Weights " << std::endl;
   os << std::setprecision(8);
   for (Int_t i=0; i<nBins; i++) {
      os << std::setw(15) << std::left << histToWrite->GetBinContent(i+1) << std::right << " ";
      if ((i+1)%5==0) os << std::endl;
   }

   os << std::setprecision(dp);
   return os; // Return the output stream.
}

////////////////////////////////////////////////////////////////////////////////
/// read the tree from an std::istream

std::istream& TMVA::operator>> ( std::istream& istr, PDF& pdf )
{
   TString devnullS;
   Int_t   valI;
   Int_t   nbins=-1; // default binning will cause an exit
   Float_t xmin=-1., xmax=-1.;
   TString hname="_original";
   Bool_t doneReading = kFALSE;
   while (!doneReading) {
      istr >> devnullS;
      if (devnullS=="NSmooth")
         {istr >> pdf.fMinNsmooth; pdf.fMaxNsmooth=pdf.fMinNsmooth;}
      else if (devnullS=="MinNSmooth") istr >> pdf.fMinNsmooth;
      else if (devnullS=="MaxNSmooth") istr >> pdf.fMaxNsmooth;
      // have to do this with strings to be more stable with developing code
      else if (devnullS == "InterpolMethod") { istr >> valI; pdf.fInterpolMethod = PDF::EInterpolateMethod(valI);}
      else if (devnullS == "KDE_type")       { istr >> valI; pdf.fKDEtype        = KDEKernel::EKernelType(valI); }
      else if (devnullS == "KDE_iter")       { istr >> valI; pdf.fKDEiter        = KDEKernel::EKernelIter(valI);}
      else if (devnullS == "KDE_border")     { istr >> valI; pdf.fKDEborder      = KDEKernel::EKernelBorder(valI);}
      else if (devnullS == "KDE_finefactor") {
         istr  >> pdf.fFineFactor;
         if (pdf.GetReadingVersion() != 0 && pdf.GetReadingVersion() < TMVA_VERSION(3,7,3)) { // here we expect the histogram limits if the version is below 3.7.3. When version == 0, the newest TMVA version is assumed.
            // coverity[tainted_data_argument]
            istr  >> nbins >> xmin >> xmax;
            doneReading = kTRUE;
         }
      }
      else if (devnullS == "Histogram")     { istr  >> hname >> nbins >> xmin >> xmax; }
      else if (devnullS == "Weights")       { doneReading = kTRUE; }
   }

   TString hnameSmooth = hname;
   hnameSmooth.ReplaceAll( "_original", "_smoothed" );

   // recreate the original hist
   if (nbins==-1) {
      std::cout << "PDF, trying to create a histogram without defined binning"<< std::endl;
      std::exit(1);
   }
   TH1* newhist = new TH1F( hname,hname, nbins, xmin, xmax );
   newhist->SetDirectory(0);
   Float_t val;
   for (Int_t i=0; i<nbins; i++) {
      istr >> val;
      newhist->SetBinContent(i+1,val);
   }

   if (pdf.fHistOriginal != 0) delete pdf.fHistOriginal;
   pdf.fHistOriginal = newhist;
   pdf.fHist = (TH1F*)pdf.fHistOriginal->Clone( hnameSmooth );
   pdf.fHist->SetTitle( hnameSmooth );
   pdf.fHist->SetDirectory(0);

   if (pdf.fMinNsmooth>=0) pdf.BuildSplinePDF();
   else {
      pdf.fInterpolMethod = PDF::kKDE;
      pdf.BuildKDEPDF();
   }

   return istr;
}

TMVA::PDF*  TMVA::PDF::ThisPDF( void )
{
   // return global "this" pointer of PDF
   return GetThisPdfThreadLocal();
}
