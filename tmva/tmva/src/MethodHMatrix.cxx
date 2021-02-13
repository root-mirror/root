// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate Data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMVA::MethodHMatrix                                                   *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header file for description)                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch> - CERN, Switzerland           *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#include "TMVA/MethodHMatrix.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include "TMatrix.h"
#include "TVectorT.h"
#include "TList.h"

#include <iostream>
#include <algorithm>

REGISTER_METHOD(HMatrix)

ClassImp(TMVA::MethodHMatrix);

/*! \class TMVA::MethodHMatrix
\ingroup TMVA

  H-Matrix method, which is implemented as a simple comparison of
  chi-squared estimators for signal and background, taking into
  account the linear correlations between the input variables

  This MVA approach is used by the \f$D\emptyset \f$ collaboration (FNAL) for the
  purpose of electron identification (see, eg.,
  [hep-ex/9507007](http://arxiv.org/abs/hep-ex/9507007)).
  As it is implemented in TMVA, it is usually equivalent or worse than
  the Fisher-Mahalanobis discriminant, and it has only been added for
  the purpose of completeness.
  Two chi^2 estimators are computed for an event, each one
  for signal and background, using the estimates for the means and
  covariance matrices obtained from the training sample:<br>

\f[
\chi^2_\eta = (x_\eta(i) - \bar{x}_\eta)^T C_\eta^{-1} (x_\eta(i) - \bar{x}_\eta), \eta = S,B
\f]

  TMVA then uses as normalised analyser for event \f$ (i) \f$ the ratio:
\f[
\frac{(chi_S(i)^2 - chi_B^2(i))}{(chi_S^2(i) + chi_B^2(i))}
\f]
*/



////////////////////////////////////////////////////////////////////////////////
/// standard constructor for the H-Matrix method

   TMVA::MethodHMatrix::MethodHMatrix( const TString& jobName,
                                       const TString& methodTitle,
                                       DataSetInfo& theData,
                                       const TString& theOption )
   : TMVA::MethodBase( jobName, Types::kHMatrix, methodTitle, theData, theOption)
   ,fInvHMatrixS(0)
   ,fInvHMatrixB(0)
   ,fVecMeanS(0)
   ,fVecMeanB(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// constructor from weight file

TMVA::MethodHMatrix::MethodHMatrix( DataSetInfo& theData,
                                    const TString& theWeightFile)
   : TMVA::MethodBase( Types::kHMatrix, theData, theWeightFile)
   ,fInvHMatrixS(0)
   ,fInvHMatrixB(0)
   ,fVecMeanS(0)
   ,fVecMeanB(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// default initialization called by all constructors

void TMVA::MethodHMatrix::Init( void )
{
   //SetNormalised( kFALSE ); obsolete!

   fInvHMatrixS = new TMatrixD( GetNvar(), GetNvar() );
   fInvHMatrixB = new TMatrixD( GetNvar(), GetNvar() );
   fVecMeanS    = new TVectorD( GetNvar() );
   fVecMeanB    = new TVectorD( GetNvar() );

   // the minimum requirement to declare an event signal-like
   SetSignalReferenceCut( 0.0 );
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodHMatrix::~MethodHMatrix( void )
{
   if (NULL != fInvHMatrixS) delete fInvHMatrixS;
   if (NULL != fInvHMatrixB) delete fInvHMatrixB;
   if (NULL != fVecMeanS   ) delete fVecMeanS;
   if (NULL != fVecMeanB   ) delete fVecMeanB;
}

////////////////////////////////////////////////////////////////////////////////
/// FDA can handle classification with 2 classes and regression with one regression-target

Bool_t TMVA::MethodHMatrix::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t /*numberTargets*/ )
{
   if( type == Types::kClassification && numberClasses == 2 ) return kTRUE;
   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// MethodHMatrix options: none (apart from those implemented in MethodBase)

void TMVA::MethodHMatrix::DeclareOptions()
{
}

////////////////////////////////////////////////////////////////////////////////
/// process user options

void TMVA::MethodHMatrix::ProcessOptions()
{
}

////////////////////////////////////////////////////////////////////////////////
/// computes H-matrices for signal and background samples

void TMVA::MethodHMatrix::Train( void )
{
   // covariance matrices for signal and background
   ComputeCovariance( kTRUE,  fInvHMatrixS );
   ComputeCovariance( kFALSE, fInvHMatrixB );

   // sanity checks
   if (TMath::Abs(fInvHMatrixS->Determinant()) < 10E-24) {
      Log() << kWARNING << "<Train> H-matrix  S is almost singular with determinant= "
            << TMath::Abs(fInvHMatrixS->Determinant())
            << " did you use the variables that are linear combinations or highly correlated ???"
            << Endl;
   }
   if (TMath::Abs(fInvHMatrixB->Determinant()) < 10E-24) {
      Log() << kWARNING << "<Train> H-matrix  B is almost singular with determinant= "
            << TMath::Abs(fInvHMatrixB->Determinant())
            << " did you use the variables that are linear combinations or highly correlated ???"
            << Endl;
   }

   if (TMath::Abs(fInvHMatrixS->Determinant()) < 10E-120) {
      Log() << kFATAL << "<Train> H-matrix  S is singular with determinant= "
            << TMath::Abs(fInvHMatrixS->Determinant())
            << " did you use the variables that are linear combinations ???"
            << Endl;
   }
   if (TMath::Abs(fInvHMatrixB->Determinant()) < 10E-120) {
      Log() << kFATAL << "<Train> H-matrix  B is singular with determinant= "
            << TMath::Abs(fInvHMatrixB->Determinant())
            << " did you use the variables that are linear combinations ???"
            << Endl;
   }

   // invert matrix
   fInvHMatrixS->Invert();
   fInvHMatrixB->Invert();
   ExitFromTraining();
}

////////////////////////////////////////////////////////////////////////////////
/// compute covariance matrix

void TMVA::MethodHMatrix::ComputeCovariance( Bool_t isSignal, TMatrixD* mat )
{
   Data()->SetCurrentType(Types::kTraining);

   const UInt_t nvar = DataInfo().GetNVariables();
   UInt_t ivar, jvar;

   // init matrices
   TVectorD vec(nvar);        vec  *= 0;
   TMatrixD mat2(nvar, nvar); mat2 *= 0;

   // initialize internal sum-of-weights variables
   Double_t sumOfWeights = 0;
   Double_t *xval = new Double_t[nvar];

   // perform event loop
   for (Int_t i=0, iEnd=Data()->GetNEvents(); i<iEnd; ++i) {

      // retrieve the original (not transformed) event
      const Event* origEvt = Data()->GetEvent(i);
      Double_t weight = origEvt->GetWeight();

      // in case event with neg weights are to be ignored
      if (IgnoreEventsWithNegWeightsInTraining() && weight <= 0) continue;

      if (DataInfo().IsSignal(origEvt) != isSignal) continue;

      // transform the event
      GetTransformationHandler().SetTransformationReferenceClass( origEvt->GetClass() );
      const Event* ev = GetTransformationHandler().Transform( origEvt );

      // event is of good type
      sumOfWeights += weight;

      // mean values
      for (ivar=0; ivar<nvar; ivar++) xval[ivar] = ev->GetValue(ivar);

      // covariance matrix
      for (ivar=0; ivar<nvar; ivar++) {

         vec(ivar)        += xval[ivar]*weight;
         mat2(ivar, ivar) += (xval[ivar]*xval[ivar])*weight;

         for (jvar=ivar+1; jvar<nvar; jvar++) {
            mat2(ivar, jvar) += (xval[ivar]*xval[jvar])*weight;
            mat2(jvar, ivar) = mat2(ivar, jvar); // symmetric matrix
         }
      }
   }

   // variance-covariance
   for (ivar=0; ivar<nvar; ivar++) {

      if (isSignal) (*fVecMeanS)(ivar) = vec(ivar)/sumOfWeights;
      else          (*fVecMeanB)(ivar) = vec(ivar)/sumOfWeights;

      for (jvar=0; jvar<nvar; jvar++) {
         (*mat)(ivar, jvar) = mat2(ivar, jvar)/sumOfWeights - vec(ivar)*vec(jvar)/(sumOfWeights*sumOfWeights);
      }
   }

   delete [] xval;
}

////////////////////////////////////////////////////////////////////////////////
/// returns the H-matrix signal estimator

Double_t TMVA::MethodHMatrix::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   Double_t s = GetChi2( Types::kSignal     );
   Double_t b = GetChi2( Types::kBackground );

   if (s+b < 0) Log() << kFATAL << "big trouble: s+b: " << s+b << Endl;

   // cannot determine error
   NoErrorCalc(err, errUpper);

   return (b - s)/(s + b);
}

////////////////////////////////////////////////////////////////////////////////
/// compute chi2-estimator for event according to type (signal/background)

Double_t TMVA::MethodHMatrix::GetChi2( Types::ESBType type )
{
   // get original (not transformed) event

   const Event* origEvt = fTmpEvent ? fTmpEvent:Data()->GetEvent();

   // loop over variables
   UInt_t ivar(0), jvar(0), nvar(GetNvar());
   std::vector<Double_t> val( nvar );

   // transform the event according to the given type (signal/background)
   if (type==Types::kSignal)
      GetTransformationHandler().SetTransformationReferenceClass( fSignalClass     );
   else
      GetTransformationHandler().SetTransformationReferenceClass( fBackgroundClass );

   const Event* ev = GetTransformationHandler().Transform( origEvt );

   for (ivar=0; ivar<nvar; ivar++) val[ivar] = ev->GetValue( ivar );

   Double_t chi2 = 0;
   for (ivar=0; ivar<nvar; ivar++) {
      for (jvar=0; jvar<nvar; jvar++) {
         if (type == Types::kSignal)
            chi2 += ( (val[ivar] - (*fVecMeanS)(ivar))*(val[jvar] - (*fVecMeanS)(jvar))
                      * (*fInvHMatrixS)(ivar,jvar) );
         else
            chi2 += ( (val[ivar] - (*fVecMeanB)(ivar))*(val[jvar] - (*fVecMeanB)(jvar))
                      * (*fInvHMatrixB)(ivar,jvar) );
      }
   }

   // sanity check
   if (chi2 < 0) Log() << kFATAL << "<GetChi2> negative chi2: " << chi2 << Endl;

   return chi2;
}

////////////////////////////////////////////////////////////////////////////////
/// create XML description for HMatrix classification

void TMVA::MethodHMatrix::AddWeightsXMLTo( void* parent ) const
{
   void* wght = gTools().AddChild(parent, "Weights");
   gTools().WriteTVectorDToXML( wght, "VecMeanS", fVecMeanS    );
   gTools().WriteTVectorDToXML( wght, "VecMeanB", fVecMeanB    );
   gTools().WriteTMatrixDToXML( wght, "InvHMatS", fInvHMatrixS );
   gTools().WriteTMatrixDToXML( wght, "InvHMatB", fInvHMatrixB );
}

////////////////////////////////////////////////////////////////////////////////
/// read weights from XML file

void TMVA::MethodHMatrix::ReadWeightsFromXML( void* wghtnode )
{
   void* descnode = gTools().GetChild(wghtnode);
   gTools().ReadTVectorDFromXML( descnode, "VecMeanS", fVecMeanS    );
   descnode = gTools().GetNextChild(descnode);
   gTools().ReadTVectorDFromXML( descnode, "VecMeanB", fVecMeanB    );
   descnode = gTools().GetNextChild(descnode);
   gTools().ReadTMatrixDFromXML( descnode, "InvHMatS", fInvHMatrixS );
   descnode = gTools().GetNextChild(descnode);
   gTools().ReadTMatrixDFromXML( descnode, "InvHMatB", fInvHMatrixB );
}

////////////////////////////////////////////////////////////////////////////////
/// read variable names and min/max
/// NOTE: the latter values are mandatory for the normalisation
/// in the reader application !!!

void  TMVA::MethodHMatrix::ReadWeightsFromStream( std::istream& istr )
{
   UInt_t ivar,jvar;
   TString var, dummy;
   istr >> dummy;
   //this->SetMethodName(dummy);

   // mean vectors
   for (ivar=0; ivar<GetNvar(); ivar++)
      istr >> (*fVecMeanS)(ivar) >> (*fVecMeanB)(ivar);

   // inverse covariance matrices (signal)
   for (ivar=0; ivar<GetNvar(); ivar++)
      for (jvar=0; jvar<GetNvar(); jvar++)
         istr >> (*fInvHMatrixS)(ivar,jvar);

   // inverse covariance matrices (background)
   for (ivar=0; ivar<GetNvar(); ivar++)
      for (jvar=0; jvar<GetNvar(); jvar++)
         istr >> (*fInvHMatrixB)(ivar,jvar);
}

////////////////////////////////////////////////////////////////////////////////
/// write Fisher-specific classifier response

void TMVA::MethodHMatrix::MakeClassSpecific( std::ostream& fout, const TString& className ) const
{
   fout << "   // arrays of input evt vs. variable " << std::endl;
   fout << "   double fInvHMatrixS[" << GetNvar() << "][" << GetNvar() << "]; // inverse H-matrix (signal)" << std::endl;
   fout << "   double fInvHMatrixB[" << GetNvar() << "][" << GetNvar() << "]; // inverse H-matrix (background)" << std::endl;
   fout << "   double fVecMeanS[" << GetNvar() << "];    // vector of mean values (signal)" << std::endl;
   fout << "   double fVecMeanB[" << GetNvar() << "];    // vector of mean values (background)" << std::endl;
   fout << "   " << std::endl;
   fout << "   double GetChi2( const std::vector<double>& inputValues, int type ) const;" << std::endl;
   fout << "};" << std::endl;
   fout << "   " << std::endl;
   fout << "void " << className << "::Initialize() " << std::endl;
   fout << "{" << std::endl;
   fout << "   // init vectors with mean values" << std::endl;
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      fout << "   fVecMeanS[" << ivar << "] = " << (*fVecMeanS)(ivar) << ";" << std::endl;
      fout << "   fVecMeanB[" << ivar << "] = " << (*fVecMeanB)(ivar) << ";" << std::endl;
   }
   fout << "   " << std::endl;
   fout << "   // init H-matrices" << std::endl;
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      for (UInt_t jvar=0; jvar<GetNvar(); jvar++) {
         fout << "   fInvHMatrixS[" << ivar << "][" << jvar << "] = "
              << (*fInvHMatrixS)(ivar,jvar) << ";" << std::endl;
         fout << "   fInvHMatrixB[" << ivar << "][" << jvar << "] = "
              << (*fInvHMatrixB)(ivar,jvar) << ";" << std::endl;
      }
   }
   fout << "}" << std::endl;
   fout << "   " << std::endl;
   fout << "inline double " << className << "::GetMvaValue__( const std::vector<double>& inputValues ) const" << std::endl;
   fout << "{" << std::endl;
   fout << "   // returns the H-matrix signal estimator" << std::endl;
   fout << "   std::vector<double> inputValuesSig = inputValues;" << std::endl;
   fout << "   std::vector<double> inputValuesBgd = inputValues;" << std::endl;
   if (GetTransformationHandler().GetTransformationList().GetSize() != 0) {

      UInt_t signalClass    =DataInfo().GetClassInfo("Signal")->GetNumber();
      UInt_t backgroundClass=DataInfo().GetClassInfo("Background")->GetNumber();

      fout << "   Transform(inputValuesSig," << signalClass << ");" << std::endl;
      fout << "   Transform(inputValuesBgd," << backgroundClass << ");" << std::endl;
   }

   //   fout << "   for(uint i=0; i<GetNvar(); ++i) std::cout << inputValuesSig.at(i) << \"  \" << inputValuesBgd.at(i) << std::endl; " << std::endl;

   fout << "   double s = GetChi2( inputValuesSig, " << Types::kSignal << " );" << std::endl;
   fout << "   double b = GetChi2( inputValuesBgd, " << Types::kBackground << " );" << std::endl;

   //   fout << "   std::cout << s << \"  \" << b << std::endl; " << std::endl;

   fout << "   " << std::endl;
   fout << "   if (s+b <= 0) std::cout << \"Problem in class " << className << "::GetMvaValue__: s+b = \"" << std::endl;
   fout << "                           << s+b << \" <= 0 \"  << std::endl;" << std::endl;
   fout << "   " << std::endl;
   fout << "   return (b - s)/(s + b);" << std::endl;
   fout << "}" << std::endl;
   fout << "   " << std::endl;
   fout << "inline double " << className << "::GetChi2( const std::vector<double>& inputValues, int type ) const" << std::endl;
   fout << "{" << std::endl;
   fout << "   // compute chi2-estimator for event according to type (signal/background)" << std::endl;
   fout << "   " << std::endl;
   fout << "   size_t ivar,jvar;" << std::endl;
   fout << "   double chi2 = 0;" << std::endl;
   fout << "   for (ivar=0; ivar<GetNvar(); ivar++) {" << std::endl;
   fout << "      for (jvar=0; jvar<GetNvar(); jvar++) {" << std::endl;
   fout << "         if (type == " << Types::kSignal << ") " << std::endl;
   fout << "            chi2 += ( (inputValues[ivar] - fVecMeanS[ivar])*(inputValues[jvar] - fVecMeanS[jvar])" << std::endl;
   fout << "                      * fInvHMatrixS[ivar][jvar] );" << std::endl;
   fout << "         else" << std::endl;
   fout << "            chi2 += ( (inputValues[ivar] - fVecMeanB[ivar])*(inputValues[jvar] - fVecMeanB[jvar])" << std::endl;
   fout << "                      * fInvHMatrixB[ivar][jvar] );" << std::endl;
   fout << "      }" << std::endl;
   fout << "   }   // loop over variables   " << std::endl;
   fout << "   " << std::endl;
   fout << "   // sanity check" << std::endl;
   fout << "   if (chi2 < 0) std::cout << \"Problem in class " << className << "::GetChi2: chi2 = \"" << std::endl;
   fout << "                           << chi2 << \" < 0 \"  << std::endl;" << std::endl;
   fout << "   " << std::endl;
   fout << "   return chi2;" << std::endl;
   fout << "}" << std::endl;
   fout << "   " << std::endl;
   fout << "// Clean up" << std::endl;
   fout << "inline void " << className << "::Clear() " << std::endl;
   fout << "{" << std::endl;
   fout << "   // nothing to clear" << std::endl;
   fout << "}" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// get help message text
///
/// typical length of text line:
///         "|--------------------------------------------------------------|"

void TMVA::MethodHMatrix::GetHelpMessage() const
{
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "The H-Matrix classifier discriminates one class (signal) of a feature" << Endl;
   Log() << "vector from another (background). The correlated elements of the" << Endl;
   Log() << "vector are assumed to be Gaussian distributed, and the inverse of" << Endl;
   Log() << "the covariance matrix is the H-Matrix. A multivariate chi-squared" << Endl;
   Log() << "estimator is built that exploits differences in the mean values of" << Endl;
   Log() << "the vector elements between the two classes for the purpose of" << Endl;
   Log() << "discrimination." << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance optimisation:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "The TMVA implementation of the H-Matrix classifier has been shown" << Endl;
   Log() << "to underperform in comparison with the corresponding Fisher discriminant," << Endl;
   Log() << "when using similar assumptions and complexity. Its use is therefore" << Endl;
   Log() << "depreciated. Only in cases where the background model is strongly" << Endl;
   Log() << "non-Gaussian, H-Matrix may perform better than Fisher. In such" << Endl;
   Log() << "occurrences the user is advised to employ non-linear classifiers. " << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance tuning via configuration options:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "None" << Endl;
}
