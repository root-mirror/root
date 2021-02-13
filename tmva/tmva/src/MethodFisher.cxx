// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Xavier Prudent, Joerg Stelzer, Helge Voss, Kai Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate Data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodFisher                                                          *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Original author of this Fisher-Discriminant implementation:                    *
 *      Andre Gaidot, CEA-France;                                                 *
 *      (Translation from FORTRAN)                                                *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Xavier Prudent  <prudent@lapp.in2p3.fr>  - LAPP, France                   *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      LAPP, Annecy, France                                                      *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::MethodFisher
\ingroup TMVA

Fisher and Mahalanobis Discriminants (Linear Discriminant Analysis)

In the method of Fisher discriminants event selection is performed
in a transformed variable space with zero linear correlations, by
distinguishing the mean values of the signal and background
distributions.

The linear discriminant analysis determines an axis in the (correlated)
hyperspace of the input variables
such that, when projecting the output classes (signal and background)
upon this axis, they are pushed as far as possible away from each other,
while events of a same class are confined in a close vicinity.
The linearity property of this method is reflected in the metric with
which "far apart" and "close vicinity" are determined: the covariance
matrix of the discriminant variable space.

The classification of the events in signal and background classes
relies on the following characteristics (only): overall sample means, \f$ x_i \f$,
for each input variable, \f$ i \f$,
class-specific sample means, \f$ x_{S(B),i}\f$,
and total covariance matrix \f$ T_{ij} \f$. The covariance matrix
can be decomposed into the sum of a _within_ (\f$ W_{ij} \f$)
and a _between-class_ (\f$ B_{ij} \f$) class matrix. They describe
the dispersion of events relative to the means of their own class (within-class
matrix), and relative to the overall sample means (between-class matrix).
The Fisher coefficients, \f$ F_i \f$, are then given by

\f[
F_i = \frac{\sqrt{N_s N_b}}{N_s + N_b} \sum_{j=1}^{N_{SB}} W_{ij}^{-1} (\bar{X}_{Sj} - \bar{X}_{Bj})
\f]

where in TMVA is set \f$ N_S = N_B \f$, so that the factor
in front of the sum simplifies to \f$ \frac{1}{2}\f$.
The Fisher discriminant then reads

\f[
X_{Fi} = F_0 + \sum_{i=1}^{N_{SB}} F_i X_i
\f]

The offset \f$ F_0 \f$ centers the sample mean of \f$ x_{Fi} \f$
at zero. Instead of using the within-class matrix, the Mahalanobis variant
determines the Fisher coefficients as follows:

\f[
F_i = \frac{\sqrt{N_s N_b}}{N_s + N_b} \sum_{j=1}^{N_{SB}} (W + B)_{ij}^{-1} (\bar{X}_{Sj} - \bar{X}_{Bj})
\f]

with resulting \f$ x_{Ma} \f$ that are very similar to the \f$ x_{Fi} \f$.

TMVA provides two outputs for the ranking of the input variables:

  -  __Fisher test:__ the Fisher analysis aims at simultaneously maximising
the between-class separation, while minimising the within-class dispersion.
A useful measure of the discrimination power of a variable is hence given
by the diagonal quantity: \f$ \frac{B_{ii}}{W_{ii}} \f$ .

  -  __Discrimination power:__ the value of the Fisher coefficient is a
measure of the discriminating power of a variable. The discrimination power
of set of input variables can therefore be measured by the scalar

\f[
\lambda = \frac{\sqrt{N_s N_b}}{N_s + N_b} \sum_{j=1}^{N_{SB}} F_i (\bar{X}_{Sj} - \bar{X}_{Bj})
\f]

The corresponding numbers are printed on standard output.
*/

#include "TMVA/MethodFisher.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/Configurable.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/Event.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Ranking.h"
#include "TMVA/Tools.h"
#include "TMVA/TransformationHandler.h"
#include "TMVA/Types.h"
#include "TMVA/VariableTransformBase.h"

#include "TMath.h"
#include "TMatrix.h"
#include "TList.h"

#include <iostream>
#include <iomanip>
#include <cassert>

REGISTER_METHOD(Fisher)

ClassImp(TMVA::MethodFisher);

////////////////////////////////////////////////////////////////////////////////
/// standard constructor for the "Fisher"

TMVA::MethodFisher::MethodFisher( const TString& jobName,
                                  const TString& methodTitle,
                                  DataSetInfo& dsi,
                                  const TString& theOption ) :
   MethodBase( jobName, Types::kFisher, methodTitle, dsi, theOption),
   fMeanMatx     ( 0 ),
   fTheMethod    ( "Fisher" ),
   fFisherMethod ( kFisher ),
   fBetw         ( 0 ),
   fWith         ( 0 ),
   fCov          ( 0 ),
   fSumOfWeightsS( 0 ),
   fSumOfWeightsB( 0 ),
   fDiscrimPow   ( 0 ),
   fFisherCoeff  ( 0 ),
   fF0           ( 0 )
{
}

////////////////////////////////////////////////////////////////////////////////
/// constructor from weight file

TMVA::MethodFisher::MethodFisher( DataSetInfo& dsi,
                                  const TString& theWeightFile) :
   MethodBase( Types::kFisher, dsi, theWeightFile),
   fMeanMatx     ( 0 ),
   fTheMethod    ( "Fisher" ),
   fFisherMethod ( kFisher ),
   fBetw         ( 0 ),
   fWith         ( 0 ),
   fCov          ( 0 ),
   fSumOfWeightsS( 0 ),
   fSumOfWeightsB( 0 ),
   fDiscrimPow   ( 0 ),
   fFisherCoeff  ( 0 ),
   fF0           ( 0 )
{
}

////////////////////////////////////////////////////////////////////////////////
/// default initialization called by all constructors

void TMVA::MethodFisher::Init( void )
{
   // allocate Fisher coefficients
   fFisherCoeff = new std::vector<Double_t>( GetNvar() );

   // the minimum requirement to declare an event signal-like
   SetSignalReferenceCut( 0.0 );

   // this is the preparation for training
   InitMatrices();
}

////////////////////////////////////////////////////////////////////////////////
/// MethodFisher options:
/// format and syntax of option string: "type"
/// where type is "Fisher" or "Mahalanobis"

void TMVA::MethodFisher::DeclareOptions()
{
   DeclareOptionRef( fTheMethod = "Fisher", "Method", "Discrimination method" );
   AddPreDefVal(TString("Fisher"));
   AddPreDefVal(TString("Mahalanobis"));
}

////////////////////////////////////////////////////////////////////////////////
/// process user options

void TMVA::MethodFisher::ProcessOptions()
{
   if (fTheMethod ==  "Fisher" ) fFisherMethod = kFisher;
   else                          fFisherMethod = kMahalanobis;

   // this is the preparation for training
   InitMatrices();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodFisher::~MethodFisher( void )
{
   if (fBetw       ) { delete fBetw; fBetw = 0; }
   if (fWith       ) { delete fWith; fWith = 0; }
   if (fCov        ) { delete fCov;  fCov = 0; }
   if (fDiscrimPow ) { delete fDiscrimPow; fDiscrimPow = 0; }
   if (fFisherCoeff) { delete fFisherCoeff; fFisherCoeff = 0; }
}

////////////////////////////////////////////////////////////////////////////////
/// Fisher can only handle classification with 2 classes

Bool_t TMVA::MethodFisher::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t /*numberTargets*/ )
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// computation of Fisher coefficients by series of matrix operations

void TMVA::MethodFisher::Train( void )
{
   // get mean value of each variables for signal, backgd and signal+backgd
   GetMean();

   // get the matrix of covariance 'within class'
   GetCov_WithinClass();

   // get the matrix of covariance 'between class'
   GetCov_BetweenClass();

   // get the matrix of covariance 'between class'
   GetCov_Full();

   //--------------------------------------------------------------

   // get the Fisher coefficients
   GetFisherCoeff();

   // get the discriminating power of each variables
   GetDiscrimPower();

   // nice output
   PrintCoefficients();

   ExitFromTraining();
}

////////////////////////////////////////////////////////////////////////////////
/// returns the Fisher value (no fixed range)

Double_t TMVA::MethodFisher::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   const Event * ev = GetEvent();
   Double_t result = fF0;
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++)
      result += (*fFisherCoeff)[ivar]*ev->GetValue(ivar);

   // cannot determine error
   NoErrorCalc(err, errUpper);

   return result;

}

////////////////////////////////////////////////////////////////////////////////
/// initialization method; creates global matrices and vectors

void TMVA::MethodFisher::InitMatrices( void )
{
   // average value of each variables for S, B, S+B
   fMeanMatx = new TMatrixD( GetNvar(), 3 );

   // the covariance 'within class' and 'between class' matrices
   fBetw = new TMatrixD( GetNvar(), GetNvar() );
   fWith = new TMatrixD( GetNvar(), GetNvar() );
   fCov  = new TMatrixD( GetNvar(), GetNvar() );

   // discriminating power
   fDiscrimPow = new std::vector<Double_t>( GetNvar() );
}

////////////////////////////////////////////////////////////////////////////////
/// compute mean values of variables in each sample, and the overall means

void TMVA::MethodFisher::GetMean( void )
{
   // initialize internal sum-of-weights variables
   fSumOfWeightsS = 0;
   fSumOfWeightsB = 0;

   const UInt_t nvar = DataInfo().GetNVariables();

   // init vectors
   Double_t* sumS = new Double_t[nvar];
   Double_t* sumB = new Double_t[nvar];
   for (UInt_t ivar=0; ivar<nvar; ivar++) { sumS[ivar] = sumB[ivar] = 0; }

   // compute sample means
   for (Int_t ievt=0; ievt<Data()->GetNEvents(); ievt++) {

      // read the Training Event into "event"
      const Event * ev = GetEvent(ievt);

      // sum of weights
      Double_t weight = ev->GetWeight();
      if (DataInfo().IsSignal(ev)) fSumOfWeightsS += weight;
      else                         fSumOfWeightsB += weight;

      Double_t* sum = DataInfo().IsSignal(ev) ? sumS : sumB;

      for (UInt_t ivar=0; ivar<nvar; ivar++) sum[ivar] += ev->GetValue( ivar )*weight;
   }

   for (UInt_t ivar=0; ivar<nvar; ivar++) {
      (*fMeanMatx)( ivar, 2 ) = sumS[ivar];
      (*fMeanMatx)( ivar, 0 ) = sumS[ivar]/fSumOfWeightsS;

      (*fMeanMatx)( ivar, 2 ) += sumB[ivar];
      (*fMeanMatx)( ivar, 1 ) = sumB[ivar]/fSumOfWeightsB;

      // signal + background
      (*fMeanMatx)( ivar, 2 ) /= (fSumOfWeightsS + fSumOfWeightsB);
   }

   //   fMeanMatx->Print();
   delete [] sumS;
   delete [] sumB;
}

////////////////////////////////////////////////////////////////////////////////
/// the matrix of covariance 'within class' reflects the dispersion of the
/// events relative to the center of gravity of their own class

void TMVA::MethodFisher::GetCov_WithinClass( void )
{
   // assert required
   assert( fSumOfWeightsS > 0 && fSumOfWeightsB > 0 );

   // product matrices (x-<x>)(y-<y>) where x;y are variables

   // init
   const Int_t nvar  = GetNvar();
   const Int_t nvar2 = nvar*nvar;
   Double_t *sumSig  = new Double_t[nvar2];
   Double_t *sumBgd  = new Double_t[nvar2];
   Double_t *xval    = new Double_t[nvar];
   memset(sumSig,0,nvar2*sizeof(Double_t));
   memset(sumBgd,0,nvar2*sizeof(Double_t));

   // 'within class' covariance
   for (Int_t ievt=0; ievt<Data()->GetNEvents(); ievt++) {

      // read the Training Event into "event"
      const Event* ev = GetEvent(ievt);

      Double_t weight = ev->GetWeight(); // may ignore events with negative weights

      for (Int_t x=0; x<nvar; x++) xval[x] = ev->GetValue( x );
      Int_t k=0;
      for (Int_t x=0; x<nvar; x++) {
         for (Int_t y=0; y<nvar; y++) {
            if (DataInfo().IsSignal(ev)) {
               Double_t v = ( (xval[x] - (*fMeanMatx)(x, 0))*(xval[y] - (*fMeanMatx)(y, 0)) )*weight;
               sumSig[k] += v;
            }else{
               Double_t v = ( (xval[x] - (*fMeanMatx)(x, 1))*(xval[y] - (*fMeanMatx)(y, 1)) )*weight;
               sumBgd[k] += v;
            }
            k++;
         }
      }
   }
   Int_t k=0;
   for (Int_t x=0; x<nvar; x++) {
      for (Int_t y=0; y<nvar; y++) {
         //(*fWith)(x, y) = (sumSig[k] + sumBgd[k])/(fSumOfWeightsS + fSumOfWeightsB);
         // HHV: I am still convinced that THIS is how it should be (below) However, while
         // the old version corresponded so nicely with LD, the FIXED version does not, unless
         // we agree to change LD. For LD, it is not "defined" to my knowledge how the weights
         // are weighted, while it is clear how the "Within" matrix for Fisher should be calculated
         // (i.e. as seen below). In order to agree with the Fisher classifier, one would have to
         // weigh signal and background such that they correspond to the same number of effective
         // (weighted) events.
         // THAT is NOT done currently, but just "event weights" are used.
         (*fWith)(x, y) = sumSig[k]/fSumOfWeightsS + sumBgd[k]/fSumOfWeightsB;
         k++;
      }
   }

   delete [] sumSig;
   delete [] sumBgd;
   delete [] xval;
}

////////////////////////////////////////////////////////////////////////////////
/// the matrix of covariance 'between class' reflects the dispersion of the
/// events of a class relative to the global center of gravity of all the class
/// hence the separation between classes

void TMVA::MethodFisher::GetCov_BetweenClass( void )
{
   // assert required
   assert( fSumOfWeightsS > 0 && fSumOfWeightsB > 0);

   Double_t prodSig, prodBgd;

   for (UInt_t x=0; x<GetNvar(); x++) {
      for (UInt_t y=0; y<GetNvar(); y++) {

         prodSig = ( ((*fMeanMatx)(x, 0) - (*fMeanMatx)(x, 2))*
                     ((*fMeanMatx)(y, 0) - (*fMeanMatx)(y, 2)) );
         prodBgd = ( ((*fMeanMatx)(x, 1) - (*fMeanMatx)(x, 2))*
                     ((*fMeanMatx)(y, 1) - (*fMeanMatx)(y, 2)) );

         (*fBetw)(x, y) = (fSumOfWeightsS*prodSig + fSumOfWeightsB*prodBgd) / (fSumOfWeightsS + fSumOfWeightsB);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// compute full covariance matrix from sum of within and between matrices

void TMVA::MethodFisher::GetCov_Full( void )
{
   for (UInt_t x=0; x<GetNvar(); x++)
      for (UInt_t y=0; y<GetNvar(); y++)
         (*fCov)(x, y) = (*fWith)(x, y) + (*fBetw)(x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Fisher = Sum { [coeff]*[variables] }
///
/// let Xs be the array of the mean values of variables for signal evts
/// let Xb be the array of the mean values of variables for backgd evts
/// let InvWith be the inverse matrix of the 'within class' correlation matrix
///
/// then the array of Fisher coefficients is
/// [coeff] =sqrt(fNsig*fNbgd)/fNevt*transpose{Xs-Xb}*InvWith

void TMVA::MethodFisher::GetFisherCoeff( void )
{
   // assert required
   assert( fSumOfWeightsS > 0 && fSumOfWeightsB > 0);

   // invert covariance matrix
   TMatrixD* theMat = 0;
   switch (GetFisherMethod()) {
   case kFisher:
      theMat = fWith;
      break;
   case kMahalanobis:
      theMat = fCov;
      break;
   default:
      Log() << kFATAL << "<GetFisherCoeff> undefined method" << GetFisherMethod() << Endl;
   }

   TMatrixD invCov( *theMat );

   if ( TMath::Abs(invCov.Determinant()) < 10E-24 ) {
      Log() << kWARNING << "<GetFisherCoeff> matrix is almost singular with determinant="
            << TMath::Abs(invCov.Determinant())
            << " did you use the variables that are linear combinations or highly correlated?"
            << Endl;
   }
   if ( TMath::Abs(invCov.Determinant()) < 10E-120 ) {
      theMat->Print();
      Log() << kFATAL << "<GetFisherCoeff> matrix is singular with determinant="
            << TMath::Abs(invCov.Determinant())
            << " did you use the variables that are linear combinations? \n"
            << " do you any clue as to what went wrong in above printout of the covariance matrix? "
            << Endl;
   }

   invCov.Invert();

   // apply rescaling factor
   Double_t xfact = TMath::Sqrt( fSumOfWeightsS*fSumOfWeightsB ) / (fSumOfWeightsS + fSumOfWeightsB);

   // compute difference of mean values
   std::vector<Double_t> diffMeans( GetNvar() );
   UInt_t ivar, jvar;
   for (ivar=0; ivar<GetNvar(); ivar++) {
      (*fFisherCoeff)[ivar] = 0;

      for (jvar=0; jvar<GetNvar(); jvar++) {
         Double_t d = (*fMeanMatx)(jvar, 0) - (*fMeanMatx)(jvar, 1);
         (*fFisherCoeff)[ivar] += invCov(ivar, jvar)*d;
      }
      // rescale
      (*fFisherCoeff)[ivar] *= xfact;
   }


   // offset correction
   fF0 = 0.0;
   for (ivar=0; ivar<GetNvar(); ivar++){
      fF0 += (*fFisherCoeff)[ivar]*((*fMeanMatx)(ivar, 0) + (*fMeanMatx)(ivar, 1));
   }
   fF0 /= -2.0;
}

////////////////////////////////////////////////////////////////////////////////
/// computation of discrimination power indicator for each variable
/// small values of "fWith" indicates little compactness of sig & of backgd
/// big values of "fBetw" indicates large separation between sig & backgd
///
/// we want signal & backgd classes as compact and separated as possible
/// the discriminating power is then defined as the ration "fBetw/fWith"

void TMVA::MethodFisher::GetDiscrimPower( void )
{
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      if ((*fCov)(ivar, ivar) != 0)
         (*fDiscrimPow)[ivar] = (*fBetw)(ivar, ivar)/(*fCov)(ivar, ivar);
      else
         (*fDiscrimPow)[ivar] = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// computes ranking of input variables

const TMVA::Ranking* TMVA::MethodFisher::CreateRanking()
{
   // create the ranking object
   fRanking = new Ranking( GetName(), "Discr. power" );

   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      fRanking->AddRank( Rank( GetInputLabel(ivar), (*fDiscrimPow)[ivar] ) );
   }

   return fRanking;
}

////////////////////////////////////////////////////////////////////////////////
/// display Fisher coefficients and discriminating power for each variable
/// check maximum length of variable name

void TMVA::MethodFisher::PrintCoefficients( void )
{
   Log() << kHEADER << "Results for Fisher coefficients:" << Endl;

   if (GetTransformationHandler().GetTransformationList().GetSize() != 0) {
      Log() << kINFO << "NOTE: The coefficients must be applied to TRANFORMED variables" << Endl;
      Log() << kINFO << "  List of the transformation: " << Endl;
      TListIter trIt(&GetTransformationHandler().GetTransformationList());
      while (VariableTransformBase *trf = (VariableTransformBase*) trIt()) {
         Log() << kINFO << "  -- " << trf->GetName() << Endl;
      }
   }
   std::vector<TString>  vars;
   std::vector<Double_t> coeffs;
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      vars  .push_back( GetInputLabel(ivar) );
      coeffs.push_back(  (*fFisherCoeff)[ivar] );
   }
   vars  .push_back( "(offset)" );
   coeffs.push_back( fF0 );
   TMVA::gTools().FormattedOutput( coeffs, vars, "Variable" , "Coefficient", Log() );

   // for (int i=0; i<coeffs.size(); i++)
   //    std::cout << "fisher coeff["<<i<<"]="<<coeffs[i]<<std::endl;

   if (IsNormalised()) {
      Log() << kINFO << "NOTE: You have chosen to use the \"Normalise\" booking option. Hence, the" << Endl;
      Log() << kINFO << "      coefficients must be applied to NORMALISED (') variables as follows:" << Endl;
      Int_t maxL = 0;
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++) if (GetInputLabel(ivar).Length() > maxL) maxL = GetInputLabel(ivar).Length();

      // Print normalisation expression (see Tools.cxx): "2*(x - xmin)/(xmax - xmin) - 1.0"
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
         Log() << kINFO
               << std::setw(maxL+9) << TString("[") + GetInputLabel(ivar) + "]' = 2*("
               << std::setw(maxL+2) << TString("[") + GetInputLabel(ivar) + "]"
               << std::setw(3) << (GetXmin(ivar) > 0 ? " - " : " + ")
               << std::setw(6) << TMath::Abs(GetXmin(ivar)) << std::setw(3) << ")/"
               << std::setw(6) << (GetXmax(ivar) -  GetXmin(ivar) )
               << std::setw(3) << " - 1"
               << Endl;
      }
      Log() << kINFO << "The TMVA Reader will properly account for this normalisation, but if the" << Endl;
      Log() << kINFO << "Fisher classifier is applied outside the Reader, the transformation must be" << Endl;
      Log() << kINFO << "implemented -- or the \"Normalise\" option is removed and Fisher retrained." << Endl;
      Log() << kINFO << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// read Fisher coefficients from weight file

void TMVA::MethodFisher::ReadWeightsFromStream( std::istream& istr )
{
   istr >> fF0;
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) istr >> (*fFisherCoeff)[ivar];
}

////////////////////////////////////////////////////////////////////////////////
/// create XML description of Fisher classifier

void TMVA::MethodFisher::AddWeightsXMLTo( void* parent ) const
{
   void* wght = gTools().AddChild(parent, "Weights");
   gTools().AddAttr( wght, "NCoeff", GetNvar()+1 );
   void* coeffxml = gTools().AddChild(wght, "Coefficient");
   gTools().AddAttr( coeffxml, "Index", 0   );
   gTools().AddAttr( coeffxml, "Value", fF0 );
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      coeffxml = gTools().AddChild( wght, "Coefficient" );
      gTools().AddAttr( coeffxml, "Index", ivar+1 );
      gTools().AddAttr( coeffxml, "Value", (*fFisherCoeff)[ivar] );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// read Fisher coefficients from xml weight file

void TMVA::MethodFisher::ReadWeightsFromXML( void* wghtnode )
{
   UInt_t ncoeff, coeffidx;
   gTools().ReadAttr( wghtnode, "NCoeff", ncoeff );
   fFisherCoeff->resize(ncoeff-1);

   void* ch = gTools().GetChild(wghtnode);
   Double_t coeff;
   while (ch) {
      gTools().ReadAttr( ch, "Index", coeffidx );
      gTools().ReadAttr( ch, "Value", coeff    );
      if (coeffidx==0) fF0 = coeff;
      else             (*fFisherCoeff)[coeffidx-1] = coeff;
      ch = gTools().GetNextChild(ch);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// write Fisher-specific classifier response

void TMVA::MethodFisher::MakeClassSpecific( std::ostream& fout, const TString& className ) const
{
   Int_t dp = fout.precision();
   fout << "   double              fFisher0;" << std::endl;
   fout << "   std::vector<double> fFisherCoefficients;" << std::endl;
   fout << "};" << std::endl;
   fout << "" << std::endl;
   fout << "inline void " << className << "::Initialize() " << std::endl;
   fout << "{" << std::endl;
   fout << "   fFisher0 = " << std::setprecision(12) << fF0 << ";" << std::endl;
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      fout << "   fFisherCoefficients.push_back( " << std::setprecision(12) << (*fFisherCoeff)[ivar] << " );" << std::endl;
   }
   fout << std::endl;
   fout << "   // sanity check" << std::endl;
   fout << "   if (fFisherCoefficients.size() != fNvars) {" << std::endl;
   fout << "      std::cout << \"Problem in class \\\"\" << fClassName << \"\\\"::Initialize: mismatch in number of input values\"" << std::endl;
   fout << "                << fFisherCoefficients.size() << \" != \" << fNvars << std::endl;" << std::endl;
   fout << "      fStatusIsClean = false;" << std::endl;
   fout << "   }         " << std::endl;
   fout << "}" << std::endl;
   fout << std::endl;
   fout << "inline double " << className << "::GetMvaValue__( const std::vector<double>& inputValues ) const" << std::endl;
   fout << "{" << std::endl;
   fout << "   double retval = fFisher0;" << std::endl;
   fout << "   for (size_t ivar = 0; ivar < fNvars; ivar++) {" << std::endl;
   fout << "      retval += fFisherCoefficients[ivar]*inputValues[ivar];" << std::endl;
   fout << "   }" << std::endl;
   fout << std::endl;
   fout << "   return retval;" << std::endl;
   fout << "}" << std::endl;
   fout << std::endl;
   fout << "// Clean up" << std::endl;
   fout << "inline void " << className << "::Clear() " << std::endl;
   fout << "{" << std::endl;
   fout << "   // clear coefficients" << std::endl;
   fout << "   fFisherCoefficients.clear(); " << std::endl;
   fout << "}" << std::endl;
   fout << std::setprecision(dp);
}

////////////////////////////////////////////////////////////////////////////////
/// get help message text
///
/// typical length of text line:
///         "|--------------------------------------------------------------|"

void TMVA::MethodFisher::GetHelpMessage() const
{
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "Fisher discriminants select events by distinguishing the mean " << Endl;
   Log() << "values of the signal and background distributions in a trans- " << Endl;
   Log() << "formed variable space where linear correlations are removed." << Endl;
   Log() << Endl;
   Log() << "   (More precisely: the \"linear discriminator\" determines" << Endl;
   Log() << "    an axis in the (correlated) hyperspace of the input " << Endl;
   Log() << "    variables such that, when projecting the output classes " << Endl;
   Log() << "    (signal and background) upon this axis, they are pushed " << Endl;
   Log() << "    as far as possible away from each other, while events" << Endl;
   Log() << "    of a same class are confined in a close vicinity. The  " << Endl;
   Log() << "    linearity property of this classifier is reflected in the " << Endl;
   Log() << "    metric with which \"far apart\" and \"close vicinity\" are " << Endl;
   Log() << "    determined: the covariance matrix of the discriminating" << Endl;
   Log() << "    variable space.)" << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance optimisation:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "Optimal performance for Fisher discriminants is obtained for " << Endl;
   Log() << "linearly correlated Gaussian-distributed variables. Any deviation" << Endl;
   Log() << "from this ideal reduces the achievable separation power. In " << Endl;
   Log() << "particular, no discrimination at all is achieved for a variable" << Endl;
   Log() << "that has the same sample mean for signal and background, even if " << Endl;
   Log() << "the shapes of the distributions are very different. Thus, Fisher " << Endl;
   Log() << "discriminants often benefit from suitable transformations of the " << Endl;
   Log() << "input variables. For example, if a variable x in [-1,1] has a " << Endl;
   Log() << "a parabolic signal distributions, and a uniform background" << Endl;
   Log() << "distributions, their mean value is zero in both cases, leading " << Endl;
   Log() << "to no separation. The simple transformation x -> |x| renders this " << Endl;
   Log() << "variable powerful for the use in a Fisher discriminant." << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance tuning via configuration options:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "<None>" << Endl;
}
