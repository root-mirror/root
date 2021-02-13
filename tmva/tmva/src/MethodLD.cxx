// @(#)root/tmva $Id$
// Author: Krzysztof Danielowski, Kamil Kraszewski, Maciej Kruk, Jan Therhaag

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodLD                                                              *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Linear Discriminant - Simple Linear Regression and Classification         *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Krzysztof Danielowski   <danielow@cern.ch>      - IFJ PAN & AGH, Poland   *
 *      Kamil Kraszewski        <kalq@cern.ch>          - IFJ PAN & UJ, Poland    *
 *      Maciej Kruk             <mkruk@cern.ch>         - IFJ PAN & AGH, Poland   *
 *      Jan Therhaag            <therhaag@physik.uni-bonn.de> - Uni Bonn, Germany *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      PAN, Poland                                                               *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 *                                                                                *
 **********************************************************************************/

/*! \class TMVA::MethodLD
\ingroup TMVA
Linear Discriminant.

Can compute multidimensional output for regression
(although it computes every dimension separately)
*/

#include "TMVA/MethodLD.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/Configurable.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/PDF.h"
#include "TMVA/Ranking.h"
#include "TMVA/Tools.h"
#include "TMVA/TransformationHandler.h"
#include "TMVA/Types.h"
#include "TMVA/VariableTransformBase.h"

#include "TMath.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "TList.h"

#include <iostream>
#include <iomanip>

using std::vector;

REGISTER_METHOD(LD)

ClassImp(TMVA::MethodLD);

////////////////////////////////////////////////////////////////////////////////
/// standard constructor for the LD

   TMVA::MethodLD::MethodLD( const TString& jobName,
                             const TString& methodTitle,
                             DataSetInfo& dsi,
                             const TString& theOption ) :
   MethodBase( jobName, Types::kLD, methodTitle, dsi, theOption),
   fNRegOut   ( 0 ),
   fSumMatx   ( 0 ),
   fSumValMatx( 0 ),
   fCoeffMatx ( 0 ),
   fLDCoeff   ( 0 )
{
}

////////////////////////////////////////////////////////////////////////////////
/// constructor from weight file

TMVA::MethodLD::MethodLD( DataSetInfo& theData, const TString& theWeightFile)
   : MethodBase( Types::kLD, theData, theWeightFile),
     fNRegOut   ( 0 ),
     fSumMatx   ( 0 ),
     fSumValMatx( 0 ),
     fCoeffMatx ( 0 ),
     fLDCoeff   ( 0 )
{
}

////////////////////////////////////////////////////////////////////////////////
/// default initialization called by all constructors

void TMVA::MethodLD::Init( void )
{
   if(DataInfo().GetNTargets()!=0) fNRegOut = DataInfo().GetNTargets();
   else                fNRegOut = 1;

   fLDCoeff = new vector< vector< Double_t >* >(fNRegOut);
   for (Int_t iout = 0; iout<fNRegOut; iout++){
      (*fLDCoeff)[iout] = new std::vector<Double_t>( GetNvar()+1 );
   }

   // the minimum requirement to declare an event signal-like
   SetSignalReferenceCut( 0.0 );
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodLD::~MethodLD( void )
{
   if (fSumMatx)    { delete fSumMatx;    fSumMatx    = 0; }
   if (fSumValMatx) { delete fSumValMatx; fSumValMatx = 0; }
   if (fCoeffMatx)  { delete fCoeffMatx;  fCoeffMatx  = 0; }
   if (fLDCoeff) {
      for (vector< vector< Double_t >* >::iterator vi=fLDCoeff->begin(); vi!=fLDCoeff->end(); ++vi){
         if (*vi) { delete *vi; *vi = 0; }
      }
      delete fLDCoeff; fLDCoeff = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// LD can handle classification with 2 classes and regression with one regression-target

Bool_t TMVA::MethodLD::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets )
{
   if      (type == Types::kClassification && numberClasses == 2) return kTRUE;
   else if (type == Types::kRegression     && numberTargets == 1) {
      Log() << "regression with " << numberTargets << " targets.";
      return kTRUE;
   }
   else return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// compute fSumMatx

void TMVA::MethodLD::Train( void )
{
   GetSum();

   // compute fSumValMatx
   GetSumVal();

   // compute fCoeffMatx and fLDCoeff
   GetLDCoeff();

   // nice output
   PrintCoefficients();

   ExitFromTraining();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the MVA classification output

Double_t TMVA::MethodLD::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   const Event* ev = GetEvent();

   if (fRegressionReturnVal == NULL) fRegressionReturnVal = new vector< Float_t >();
   fRegressionReturnVal->resize( fNRegOut );

   for (Int_t iout = 0; iout<fNRegOut; iout++) {
      (*fRegressionReturnVal)[iout] = (*(*fLDCoeff)[iout])[0] ;

      int icoeff=0;
      for (std::vector<Float_t>::const_iterator it = ev->GetValues().begin();it!=ev->GetValues().end();++it){
         (*fRegressionReturnVal)[iout] += (*(*fLDCoeff)[iout])[++icoeff] * (*it);
      }
   }

   // cannot determine error
   NoErrorCalc(err, errUpper);

   return (*fRegressionReturnVal)[0];
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates the regression output

const std::vector< Float_t >& TMVA::MethodLD::GetRegressionValues()
{
   const Event* ev = GetEvent();

   if (fRegressionReturnVal == NULL) fRegressionReturnVal = new vector< Float_t >();
   fRegressionReturnVal->resize( fNRegOut );

   for (Int_t iout = 0; iout<fNRegOut; iout++) {
      (*fRegressionReturnVal)[iout] = (*(*fLDCoeff)[iout])[0] ;

      int icoeff = 0;
      for (std::vector<Float_t>::const_iterator it = ev->GetValues().begin();it!=ev->GetValues().end();++it){
         (*fRegressionReturnVal)[iout] += (*(*fLDCoeff)[iout])[++icoeff] * (*it);
      }
   }

   // perform inverse transformation
   Event* evT = new Event(*ev);
   for (Int_t iout = 0; iout<fNRegOut; iout++) evT->SetTarget(iout,(*fRegressionReturnVal)[iout]);

   const Event* evT2 = GetTransformationHandler().InverseTransform( evT );
   fRegressionReturnVal->clear();
   for (Int_t iout = 0; iout<fNRegOut; iout++) fRegressionReturnVal->push_back(evT2->GetTarget(iout));

   delete evT;
   return (*fRegressionReturnVal);
}

////////////////////////////////////////////////////////////////////////////////
/// Initialization method; creates global matrices and vectors

void TMVA::MethodLD::InitMatrices( void )
{
   fSumMatx    = new TMatrixD( GetNvar()+1, GetNvar()+1 );
   fSumValMatx = new TMatrixD( GetNvar()+1, fNRegOut );
   fCoeffMatx  = new TMatrixD( GetNvar()+1, fNRegOut );

}

////////////////////////////////////////////////////////////////////////////////
/// Calculates the matrix transposed(X)*W*X with W being the diagonal weight matrix
/// and X the coordinates values

void TMVA::MethodLD::GetSum( void )
{
   const UInt_t nvar = DataInfo().GetNVariables();

   for (UInt_t ivar = 0; ivar<=nvar; ivar++){
      for (UInt_t jvar = 0; jvar<=nvar; jvar++) (*fSumMatx)( ivar, jvar ) = 0;
   }

   // compute sample means
   Long64_t nevts = Data()->GetNEvents();
   for (Int_t ievt=0; ievt<nevts; ievt++) {
      const Event * ev = GetEvent(ievt);
      Double_t weight = ev->GetWeight();

      if (IgnoreEventsWithNegWeightsInTraining() && weight <= 0) continue;

      // Sum of weights
      (*fSumMatx)( 0, 0 ) += weight;

      // Sum of coordinates
      for (UInt_t ivar=0; ivar<nvar; ivar++) {
         (*fSumMatx)( ivar+1, 0 ) += ev->GetValue( ivar ) * weight;
         (*fSumMatx)( 0, ivar+1 ) += ev->GetValue( ivar ) * weight;
      }

      // Sum of products of coordinates
      for (UInt_t ivar=0; ivar<nvar; ivar++){
         for (UInt_t jvar=0; jvar<nvar; jvar++){
            (*fSumMatx)( ivar+1, jvar+1 ) += ev->GetValue( ivar ) * ev->GetValue( jvar ) * weight;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates the vector transposed(X)*W*Y with Y being the target vector

void TMVA::MethodLD::GetSumVal( void )
{
   const UInt_t nvar = DataInfo().GetNVariables();

   for (Int_t ivar = 0; ivar<fNRegOut; ivar++){
      for (UInt_t jvar = 0; jvar<=nvar; jvar++){
         (*fSumValMatx)(jvar,ivar) = 0;
      }
   }

   // Sum of coordinates multiplied by values
   for (Int_t ievt=0; ievt<Data()->GetNEvents(); ievt++) {

      // retrieve the event
      const Event* ev = GetEvent(ievt);
      Double_t weight = ev->GetWeight();

      // in case event with neg weights are to be ignored
      if (IgnoreEventsWithNegWeightsInTraining() && weight <= 0) continue;

      for (Int_t ivar=0; ivar<fNRegOut; ivar++) {

         Double_t val = weight;

         if (!DoRegression()){
            val *= DataInfo().IsSignal(ev); // yes it works.. but I'm still surprised (Helge).. would have not set y_B to zero though..
         }else {//for regression
            val *= ev->GetTarget( ivar );
         }
         (*fSumValMatx)( 0,ivar ) += val;
         for (UInt_t jvar=0; jvar<nvar; jvar++) {
            (*fSumValMatx)(jvar+1,ivar ) += ev->GetValue(jvar) * val;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates the coefficients used for classification/regression

void TMVA::MethodLD::GetLDCoeff( void )
{
   const UInt_t nvar = DataInfo().GetNVariables();

   for (Int_t ivar = 0; ivar<fNRegOut; ivar++){
      TMatrixD invSum( *fSumMatx );
      if ( TMath::Abs(invSum.Determinant()) < 10E-24 ) {
         Log() << kWARNING << "<GetCoeff> matrix is almost singular with determinant="
               << TMath::Abs(invSum.Determinant())
               << " did you use the variables that are linear combinations or highly correlated?"
               << Endl;
      }
      if ( TMath::Abs(invSum.Determinant()) < 10E-120 ) {
         Log() << kFATAL << "<GetCoeff> matrix is singular with determinant="
               << TMath::Abs(invSum.Determinant())
               << " did you use the variables that are linear combinations?"
               << Endl;
      }
      invSum.Invert();

      fCoeffMatx = new TMatrixD( invSum * (*fSumValMatx));
      for (UInt_t jvar = 0; jvar<nvar+1; jvar++) {
         (*(*fLDCoeff)[ivar])[jvar] = (*fCoeffMatx)(jvar, ivar );
      }
      if (!DoRegression()) {
         (*(*fLDCoeff)[ivar])[0]=0.0;
         for (UInt_t jvar = 1; jvar<nvar+1; jvar++){
            (*(*fLDCoeff)[ivar])[0]+=(*fCoeffMatx)(jvar,ivar)*(*fSumMatx)(0,jvar)/(*fSumMatx)( 0, 0 );
         }
         (*(*fLDCoeff)[ivar])[0]/=-2.0;
      }

   }
}

////////////////////////////////////////////////////////////////////////////////
/// read LD coefficients from weight file

void  TMVA::MethodLD::ReadWeightsFromStream( std::istream& istr )
{
   for (Int_t iout=0; iout<fNRegOut; iout++){
      for (UInt_t icoeff=0; icoeff<GetNvar()+1; icoeff++){
         istr >> (*(*fLDCoeff)[iout])[icoeff];
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// create XML description for LD classification and regression
/// (for arbitrary number of output classes/targets)

void TMVA::MethodLD::AddWeightsXMLTo( void* parent ) const
{
   void* wght = gTools().AddChild(parent, "Weights");
   gTools().AddAttr( wght, "NOut",   fNRegOut    );
   gTools().AddAttr( wght, "NCoeff", GetNvar()+1 );
   for (Int_t iout=0; iout<fNRegOut; iout++) {
      for (UInt_t icoeff=0; icoeff<GetNvar()+1; icoeff++) {
         void* coeffxml = gTools().AddChild( wght, "Coefficient" );
         gTools().AddAttr( coeffxml, "IndexOut",   iout   );
         gTools().AddAttr( coeffxml, "IndexCoeff", icoeff );
         gTools().AddAttr( coeffxml, "Value",      (*(*fLDCoeff)[iout])[icoeff] );
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// read coefficients from xml weight file

void TMVA::MethodLD::ReadWeightsFromXML( void* wghtnode )
{
   UInt_t ncoeff;
   gTools().ReadAttr( wghtnode, "NOut",   fNRegOut );
   gTools().ReadAttr( wghtnode, "NCoeff", ncoeff   );

   // sanity checks
   if (ncoeff != GetNvar()+1) Log() << kFATAL << "Mismatch in number of output variables/coefficients: "
                                    << ncoeff << " != " << GetNvar()+1 << Endl;

   // create vector with coefficients (double vector due to arbitrary output dimension)
   if (fLDCoeff) {
      for (vector< vector< Double_t >* >::iterator vi=fLDCoeff->begin(); vi!=fLDCoeff->end(); ++vi){
         if (*vi) { delete *vi; *vi = 0; }
      }
      delete fLDCoeff; fLDCoeff = 0;
   }
   fLDCoeff = new vector< vector< Double_t >* >(fNRegOut);
   for (Int_t ivar = 0; ivar<fNRegOut; ivar++) (*fLDCoeff)[ivar] = new std::vector<Double_t>( ncoeff );

   void* ch = gTools().GetChild(wghtnode);
   Double_t coeff;
   Int_t iout, icoeff;
   while (ch) {
      gTools().ReadAttr( ch, "IndexOut",   iout   );
      gTools().ReadAttr( ch, "IndexCoeff", icoeff );
      gTools().ReadAttr( ch, "Value",      coeff  );

      (*(*fLDCoeff)[iout])[icoeff] = coeff;

      ch = gTools().GetNextChild(ch);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// write LD-specific classifier response

void TMVA::MethodLD::MakeClassSpecific( std::ostream& fout, const TString& className ) const
{
   fout << "   std::vector<double> fLDCoefficients;" << std::endl;
   fout << "};" << std::endl;
   fout << "" << std::endl;
   fout << "inline void " << className << "::Initialize() " << std::endl;
   fout << "{" << std::endl;
   for (UInt_t ivar=0; ivar<GetNvar()+1; ivar++) {
      Int_t dp = fout.precision();
      fout << "   fLDCoefficients.push_back( "
           << std::setprecision(12) << (*(*fLDCoeff)[0])[ivar]
           << std::setprecision(dp) << " );" << std::endl;
   }
   fout << std::endl;
   fout << "   // sanity check" << std::endl;
   fout << "   if (fLDCoefficients.size() != fNvars+1) {" << std::endl;
   fout << "      std::cout << \"Problem in class \\\"\" << fClassName << \"\\\"::Initialize: mismatch in number of input values\"" << std::endl;
   fout << "                << fLDCoefficients.size() << \" != \" << fNvars+1 << std::endl;" << std::endl;
   fout << "      fStatusIsClean = false;" << std::endl;
   fout << "   }         " << std::endl;
   fout << "}" << std::endl;
   fout << std::endl;
   fout << "inline double " << className << "::GetMvaValue__( const std::vector<double>& inputValues ) const" << std::endl;
   fout << "{" << std::endl;
   fout << "   double retval = fLDCoefficients[0];" << std::endl;
   fout << "   for (size_t ivar = 1; ivar < fNvars+1; ivar++) {" << std::endl;
   fout << "      retval += fLDCoefficients[ivar]*inputValues[ivar-1];" << std::endl;
   fout << "   }" << std::endl;
   fout << std::endl;
   fout << "   return retval;" << std::endl;
   fout << "}" << std::endl;
   fout << std::endl;
   fout << "// Clean up" << std::endl;
   fout << "inline void " << className << "::Clear() " << std::endl;
   fout << "{" << std::endl;
   fout << "   // clear coefficients" << std::endl;
   fout << "   fLDCoefficients.clear(); " << std::endl;
   fout << "}" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////
/// computes ranking of input variables

const TMVA::Ranking* TMVA::MethodLD::CreateRanking()
{
   // create the ranking object
   fRanking = new Ranking( GetName(), "Discr. power" );

   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      fRanking->AddRank( Rank( GetInputLabel(ivar), TMath::Abs((* (*fLDCoeff)[0])[ivar+1] )) );
   }

   return fRanking;
}

////////////////////////////////////////////////////////////////////////////////
/// MethodLD options

void TMVA::MethodLD::DeclareOptions()
{
   AddPreDefVal(TString("LD"));
}

////////////////////////////////////////////////////////////////////////////////
/// this is the preparation for training

void TMVA::MethodLD::ProcessOptions()
{
   if (HasTrainingTree()) InitMatrices();
}

////////////////////////////////////////////////////////////////////////////////
/// Display the classification/regression coefficients for each variable

void TMVA::MethodLD::PrintCoefficients( void )
{
   Log() << kHEADER << "Results for LD coefficients:" << Endl;

   if (GetTransformationHandler().GetTransformationList().GetSize() != 0) {
      Log() << kINFO << "NOTE: The coefficients must be applied to TRANFORMED variables" << Endl;
      Log() << kINFO << "      List of the transformation: " << Endl;
      TListIter trIt(&GetTransformationHandler().GetTransformationList());
      while (VariableTransformBase *trf = (VariableTransformBase*) trIt() ) {
         Log() << kINFO << "  -- " << trf->GetName() << Endl;
      }
   }
   std::vector<TString>  vars;
   std::vector<Double_t> coeffs;
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      vars  .push_back( GetInputLabel(ivar) );
      coeffs.push_back( (* (*fLDCoeff)[0])[ivar+1] );
   }
   vars  .push_back( "(offset)" );
   coeffs.push_back((* (*fLDCoeff)[0])[0] );
   TMVA::gTools().FormattedOutput( coeffs, vars, "Variable" , "Coefficient", Log() );
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
      Log() << kINFO << "LD classifier is applied outside the Reader, the transformation must be" << Endl;
      Log() << kINFO << "implemented -- or the \"Normalise\" option is removed and LD retrained." << Endl;
      Log() << kINFO << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// get help message text
///
/// typical length of text line:
///         "|--------------------------------------------------------------|"

void TMVA::MethodLD::GetHelpMessage() const
{
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "Linear discriminants select events by distinguishing the mean " << Endl;
   Log() << "values of the signal and background distributions in a trans- " << Endl;
   Log() << "formed variable space where linear correlations are removed." << Endl;
   Log() << "The LD implementation here is equivalent to the \"Fisher\" discriminant" << Endl;
   Log() << "for classification, but also provides linear regression." << Endl;
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
   Log() << "Optimal performance for the linear discriminant is obtained for " << Endl;
   Log() << "linearly correlated Gaussian-distributed variables. Any deviation" << Endl;
   Log() << "from this ideal reduces the achievable separation power. In " << Endl;
   Log() << "particular, no discrimination at all is achieved for a variable" << Endl;
   Log() << "that has the same sample mean for signal and background, even if " << Endl;
   Log() << "the shapes of the distributions are very different. Thus, the linear " << Endl;
   Log() << "discriminant often benefits from a suitable transformation of the " << Endl;
   Log() << "input variables. For example, if a variable x in [-1,1] has a " << Endl;
   Log() << "a parabolic signal distributions, and a uniform background" << Endl;
   Log() << "distributions, their mean value is zero in both cases, leading " << Endl;
   Log() << "to no separation. The simple transformation x -> |x| renders this " << Endl;
   Log() << "variable powerful for the use in a linear discriminant." << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance tuning via configuration options:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "<None>" << Endl;
}
