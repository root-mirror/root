// @(#)root/mathcore:$Id$
// Author: L. Moneta Wed Aug 30 11:05:34 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Implementation file for class FitResult

#include "Fit/FitResult.h"

#include "Fit/FitConfig.h"

#include "Fit/BinData.h"

//#include "Fit/Chi2FCN.h"

#include "Math/Minimizer.h"

#include "Math/IParamFunction.h"
#include "Math/OneDimFunctionAdapter.h"

#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"

#include "TMath.h"
#include "Math/RichardsonDerivator.h"
#include "Math/Error.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace ROOT {

   namespace Fit {


const int gInitialResultStatus = -99; // use this special convention to flag it when printing result

FitResult::FitResult() :
   fValid(false), fNormalized(false), fNFree(0), fNdf(0), fNCalls(0),
   fStatus(-1), fCovStatus(0), fVal(0), fEdm(-1), fChi2(-1)
{
   // Default constructor implementation.
}

FitResult::FitResult(const FitConfig & fconfig) :
   fValid(false),
   fNormalized(false),
   fNFree(0),
   fNdf(0),
   fNCalls(0),
   fStatus(gInitialResultStatus),
   fCovStatus(0),
   fVal(0),
   fEdm(-1),
   fChi2(-1),
   fFitFunc(0),
   fParams(std::vector<double>( fconfig.NPar() ) ),
   fErrors(std::vector<double>( fconfig.NPar() ) ),
   fParNames(std::vector<std::string> ( fconfig.NPar() ) )
{
   // create a Fit result from a fit config (i.e. with initial parameter values
   // and errors equal to step values
   // The model function is NULL in this case

   // set minimizer type and algorithm
   fMinimType = fconfig.MinimizerType();
   // append algorithm name for minimizer that support it
   if ( (fMinimType.find("Fumili") == std::string::npos) &&
        (fMinimType.find("GSLMultiFit") == std::string::npos)
      ) {
      if (fconfig.MinimizerAlgoType() != "") fMinimType += " / " + fconfig.MinimizerAlgoType();
   }

   // get parameter values and errors (step sizes)
   unsigned int npar = fconfig.NPar();
   for (unsigned int i = 0; i < npar; ++i ) {
      const ParameterSettings & par = fconfig.ParSettings(i);
      fParams[i]   =  par.Value();
      fErrors[i]   =  par.StepSize();
      fParNames[i] =  par.Name();
      if (par.IsFixed() ) fFixedParams[i] = true;
      else fNFree++;
      if (par.IsBound() ) {
         double lower = (par.HasLowerLimit()) ? par.LowerLimit() : - std::numeric_limits<double>::infinity() ;
         double upper = (par.HasUpperLimit()) ? par.UpperLimit() :   std::numeric_limits<double>::infinity() ;
         fBoundParams[i] = fParamBounds.size();
         fParamBounds.push_back(std::make_pair(lower,upper));
      }
   }
   std::cout << "create fit result from config - nfree " << fNFree << std::endl;
}

void FitResult::FillResult(const std::shared_ptr<ROOT::Math::Minimizer> & min, const FitConfig & fconfig, const std::shared_ptr<IModelFunction> & func,
                     bool isValid,  unsigned int sizeOfData, bool binnedFit, const  ROOT::Math::IMultiGenFunction * chi2func, unsigned int ncalls )
{
   // Fill the FitResult after minimization using result from Minimizers

   // minimizer must exist
   assert(min);

   fValid = isValid;
   fNFree= min->NFree();
   fNCalls = min->NCalls();
   fStatus = min->Status();
   fCovStatus= min->CovMatrixStatus();
   fVal  =  min->MinValue();
   fEdm = min->Edm();

   fMinimizer= min;
   fFitFunc = func;

   fMinimType = fconfig.MinimizerName();

   // replace ncalls if minimizer does not support it (they are taken then from the FitMethodFunction)
   if (fNCalls == 0) fNCalls = ncalls;

   const unsigned int npar = min->NDim();
   if (npar == 0) return;

   if (min->X() )
      fParams = std::vector<double>(min->X(), min->X() + npar);
   else {
      // case minimizer does not provide minimum values (it failed) take from configuration
      fParams.resize(npar);
      for (unsigned int i = 0; i < npar; ++i ) {
         fParams[i] = ( fconfig.ParSettings(i).Value() );
      }
   }

   if (sizeOfData >  min->NFree() ) fNdf = sizeOfData - min->NFree();


   // set right parameters in function (in case minimizer did not do before)
   // do also when fit is not valid
   if (func ) {
      // I think we can avoid cloning the model function
      //fFitFunc = dynamic_cast<IModelFunction *>( func->Clone() );
      //assert(fFitFunc);
      fFitFunc->SetParameters(&fParams.front());
   }
   else {
      // when no fFitFunc is present take parameters from FitConfig
      fParNames.reserve( npar );
      for (unsigned int i = 0; i < npar; ++i ) {
         fParNames.push_back( fconfig.ParSettings(i).Name() );
      }
   }


   // check for fixed or limited parameters
   unsigned int nfree = 0;
   for (unsigned int ipar = 0; ipar < npar; ++ipar) {
      const ParameterSettings & par = fconfig.ParSettings(ipar);
      if (par.IsFixed() ) fFixedParams[ipar] = true;
      else nfree++;
      if (par.IsBound() ) {
         double lower = (par.HasLowerLimit()) ? par.LowerLimit() : - std::numeric_limits<double>::infinity() ;
         double upper = (par.HasUpperLimit()) ? par.UpperLimit() :   std::numeric_limits<double>::infinity() ;
         fBoundParams[ipar] = fParamBounds.size();
         fParamBounds.push_back(std::make_pair(lower,upper));
      }
   }
   // check if nfree (from FitConfig) and fNFree (from minimizer) are consistent
   if (nfree != fNFree ) {
      MATH_ERROR_MSG("FitResult","FitConfiguration and Minimizer result are not consistent");
      std::cout << "Number of free parameters from FitConfig = " << nfree << std::endl;
      std::cout << "Number of free parameters from Minimizer = " << fNFree << std::endl;
   }

   // if flag is binned compute a chi2 when a chi2 function is given
   if (binnedFit) {
      if (chi2func == 0)
         fChi2 = fVal;
      else {
         // compute chi2 equivalent for likelihood fits
         // NB: empty bins are considered
         fChi2 = (*chi2func)(&fParams[0]);
      }
   }

   // fill error matrix
   // if minimizer provides error provides also error matrix
   if (min->Errors() != 0) {

      fErrors = std::vector<double>(min->Errors(), min->Errors() + npar ) ;

      if (fCovStatus != 0) {
         unsigned int r = npar * (  npar + 1 )/2;
         fCovMatrix.reserve(r);
         for (unsigned int i = 0; i < npar; ++i)
            for (unsigned int j = 0; j <= i; ++j)
               fCovMatrix.push_back(min->CovMatrix(i,j) );
      }
      // minos errors are set separetly when calling Fitter::CalculateMinosErrors()

      // globalCC
      fGlobalCC.reserve(npar);
      for (unsigned int i = 0; i < npar; ++i) {
         double globcc = min->GlobalCC(i);
         if (globcc < 0) break; // it is not supported by that minimizer
         fGlobalCC.push_back(globcc);
      }

   }

}

FitResult::~FitResult() {
   // destructor. FitResult manages the fit Function pointer
   //if (fFitFunc) delete fFitFunc;
}

FitResult::FitResult(const FitResult &rhs) :
   fFitFunc(0)
{
   // Implementation of copy constructor
   (*this) = rhs;
}

FitResult & FitResult::operator = (const FitResult &rhs) {
   // Implementation of assignment operator.
   if (this == &rhs) return *this;  // time saving self-test

   // Manages the fitted function
   // if (fFitFunc) delete fFitFunc;
   // fFitFunc = 0;
   // if (rhs.fFitFunc != 0 ) {
   //    fFitFunc = dynamic_cast<IModelFunction *>( (rhs.fFitFunc)->Clone() );
   //    assert(fFitFunc != 0);
   // }

   // copy all other data members
   fValid = rhs.fValid;
   fNormalized = rhs.fNormalized;
   fNFree = rhs.fNFree;
   fNdf = rhs.fNdf;
   fNCalls = rhs.fNCalls;
   fCovStatus = rhs.fCovStatus;
   fStatus = rhs.fStatus;
   fVal = rhs.fVal;
   fEdm = rhs.fEdm;
   fChi2 = rhs.fChi2;
   fMinimizer = rhs.fMinimizer;
   fObjFunc = rhs.fObjFunc;
   fFitFunc = rhs.fFitFunc;

   fFixedParams = rhs.fFixedParams;
   fBoundParams = rhs.fBoundParams;
   fParamBounds = rhs.fParamBounds;
   fParams = rhs.fParams;
   fErrors = rhs.fErrors;
   fCovMatrix = rhs.fCovMatrix;
   fGlobalCC = rhs.fGlobalCC;
   fMinosErrors = rhs.fMinosErrors;

   fMinimType = rhs.fMinimType;
   fParNames = rhs.fParNames;

   return *this;

}

bool FitResult::Update(const std::shared_ptr<ROOT::Math::Minimizer> & min, const ROOT::Fit::FitConfig & fconfig, bool isValid, unsigned int ncalls) {
   // update fit result with new status from minimizer
   // ncalls if it is not zero is used instead of value from minimizer

   fMinimizer = min;

   // in case minimizer changes
   fMinimType = fconfig.MinimizerName();

   const unsigned int npar = fParams.size();
   if (min->NDim() != npar ) {
      MATH_ERROR_MSG("FitResult::Update","Wrong minimizer status ");
      return false;
   }
   if (min->X() == 0 ) {
      MATH_ERROR_MSG("FitResult::Update","Invalid minimizer status ");
      return false;
   }
   //fNFree = min->NFree();
   if (fNFree != min->NFree() ) {
      MATH_ERROR_MSG("FitResult::Update","Configuration has changed ");
      return false;
   }

   fValid = isValid;
   // update minimum value
   fVal = min->MinValue();
   fEdm = min->Edm();
   fStatus = min->Status();
   fCovStatus = min->CovMatrixStatus();

   // update number of function calls
   if ( min->NCalls() > 0)   fNCalls = min->NCalls();
   else fNCalls = ncalls;

   // copy parameter value and errors
   std::copy(min->X(), min->X() + npar, fParams.begin());


   // set parameters  in fit model function
   if (fFitFunc) fFitFunc->SetParameters(&fParams.front());

   if (min->Errors() != 0)  {

      if (fErrors.size() != npar) fErrors.resize(npar);

      std::copy(min->Errors(), min->Errors() + npar, fErrors.begin() ) ;

      if (fCovStatus != 0) {

         // update error matrix
         unsigned int r = npar * (  npar + 1 )/2;
         if (fCovMatrix.size() != r) fCovMatrix.resize(r);
         unsigned int l = 0;
         for (unsigned int i = 0; i < npar; ++i) {
            for (unsigned int j = 0; j <= i; ++j)
               fCovMatrix[l++] = min->CovMatrix(i,j);
         }
      }

      // update global CC
      if (fGlobalCC.size() != npar) fGlobalCC.resize(npar);
      for (unsigned int i = 0; i < npar; ++i) {
         double globcc = min->GlobalCC(i);
         if (globcc < 0) {
            fGlobalCC.clear();
            break; // it is not supported by that minimizer
         }
         fGlobalCC[i] = globcc;
      }

   }
   return true;
}

void FitResult::NormalizeErrors() {
   // normalize errors and covariance matrix according to chi2 value
   if (fNdf == 0 || fChi2 <= 0) return;
   double s2 = fChi2/fNdf;
   double s = std::sqrt(fChi2/fNdf);
   for (unsigned int i = 0; i < fErrors.size() ; ++i)
      fErrors[i] *= s;
   for (unsigned int i = 0; i < fCovMatrix.size() ; ++i)
      fCovMatrix[i] *= s2;

   fNormalized = true;
}


double FitResult::Prob() const {
   // fit probability
   return ROOT::Math::chisquared_cdf_c(fChi2, static_cast<double>(fNdf) );
}

bool FitResult::HasMinosError(unsigned int i) const {
   // query if the parameter i has the Minos error
   std::map<unsigned int, std::pair<double,double> >::const_iterator itr = fMinosErrors.find(i);
   return (itr !=  fMinosErrors.end() );
}


double FitResult::LowerError(unsigned int i) const {
   // return lower Minos error for parameter i
   //  return the parabolic error if Minos error has not been calculated for the parameter i
   std::map<unsigned int, std::pair<double,double> >::const_iterator itr = fMinosErrors.find(i);
   return ( itr != fMinosErrors.end() ) ? itr->second.first : Error(i) ;
}

double FitResult::UpperError(unsigned int i) const {
   // return upper Minos error for parameter i
   //  return the parabolic error if Minos error has not been calculated for the parameter i
   std::map<unsigned int, std::pair<double,double> >::const_iterator itr = fMinosErrors.find(i);
   return ( itr != fMinosErrors.end() ) ? itr->second.second : Error(i) ;
}

void FitResult::SetMinosError(unsigned int i, double elow, double eup) {
   // set the Minos error for parameter i
   fMinosErrors[i] = std::make_pair(elow,eup);
}

int FitResult::Index(const std::string & name) const {
   // find index for given parameter name
   if (! fFitFunc) return -1;
   unsigned int npar = fParams.size();
   for (unsigned int i = 0; i < npar; ++i)
      if ( fFitFunc->ParameterName(i) == name) return i;

   return -1; // case name is not found
}

bool FitResult::IsParameterBound(unsigned int ipar) const {
   return fBoundParams.find(ipar) != fBoundParams.end();
}

bool FitResult::IsParameterFixed(unsigned int ipar) const {
   return fFixedParams.find(ipar) != fFixedParams.end();
}

bool FitResult::ParameterBounds(unsigned int ipar, double & lower, double & upper) const {
   std::map<unsigned int, unsigned int>::const_iterator itr =  fBoundParams.find(ipar);
   if (itr ==  fBoundParams.end() ) {
      lower =  -std::numeric_limits<Double_t>::infinity();
      upper =  std::numeric_limits<Double_t>::infinity();
      return false;
   }
   assert(itr->second < fParamBounds.size() );
   lower = fParamBounds[itr->second].first;
   upper = fParamBounds[itr->second].second;
   return false;
}

std::string FitResult::ParName(unsigned int ipar) const {
   // return parameter name
   if (fFitFunc) return fFitFunc->ParameterName(ipar);
   else if (ipar < fParNames.size() ) return fParNames[ipar];
   return "param_" + ROOT::Math::Util::ToString(ipar);
}

void FitResult::Print(std::ostream & os, bool doCovMatrix) const {
   // print the result in the given stream
   // need to add also minos errors , globalCC, etc..
   unsigned int npar = fParams.size();
   if (npar == 0) {
      os << "<Empty FitResult>\n";
      return;
   }
   os << "\n****************************************\n";
   if (!fValid) {
      if (fStatus != gInitialResultStatus) {
         os << "         Invalid FitResult";
         os << "  (status = " << fStatus << " )";
      }
      else {
         os << "      FitResult before fitting";
      }
      os << "\n****************************************\n";
   }

   //os << "            FitResult                   \n\n";
   os << "Minimizer is " << fMinimType << std::endl;
   const unsigned int nw = 25; // spacing for text
   const unsigned int nn = 12; // spacing for numbers
   const std::ios_base::fmtflags prFmt = os.setf(std::ios::left,std::ios::adjustfield); // set left alignment

   if (fVal != fChi2 || fChi2 < 0)
      os << std::left << std::setw(nw) << "MinFCN" << " = " << std::right << std::setw(nn) << fVal << std::endl;
   if (fChi2 >= 0)
      os << std::left << std::setw(nw) <<  "Chi2"         << " = " << std::right << std::setw(nn) << fChi2 << std::endl;
   os << std::left << std::setw(nw) << "NDf"              << " = " << std::right << std::setw(nn) << fNdf << std::endl;
   if (fMinimType.find("Linear") == std::string::npos) {  // no need to print this for linear fits
      if (fEdm >=0) os << std::left << std::setw(nw) << "Edm"    << " = " << std::right << std::setw(nn) << fEdm << std::endl;
      os << std::left << std::setw(nw) << "NCalls" << " = " << std::right << std::setw(nn) << fNCalls << std::endl;
   }
   for (unsigned int i = 0; i < npar; ++i) {
      os << std::left << std::setw(nw) << GetParameterName(i);
      os << " = " << std::right << std::setw(nn) << fParams[i];
      if (IsParameterFixed(i) )
         os << std::setw(9) << " "  << std::setw(nn) << " " << " \t (fixed)";
      else {
         if (fErrors.size() != 0)
            os << "   +/-   " << std::left << std::setw(nn) << fErrors[i] << std::right;
         if (HasMinosError(i))
            os << "  " << std::left  << std::setw(nn) << LowerError(i) << " +" << std::setw(nn) << UpperError(i)
               << " (Minos) ";
         if (IsParameterBound(i))
            os << " \t (limited)";
      }
      os << std::endl;
   }

   // restore stremam adjustfield
   if (prFmt != os.flags() ) os.setf(prFmt, std::ios::adjustfield);

   if (doCovMatrix) PrintCovMatrix(os);
}

void FitResult::PrintCovMatrix(std::ostream &os) const {
   // print the covariance and correlation matrix
   if (!fValid) return;
   if (fCovMatrix.size() == 0) return;
//   os << "****************************************\n";
   os << "\nCovariance Matrix:\n\n";
   unsigned int npar = fParams.size();
   const int kPrec = 5;
   const int kWidth = 8;
   const int parw = 12;
   const int matw = kWidth+4;

   // query previous precision and format flags
   int prevPrec = os.precision(kPrec);
   const std::ios_base::fmtflags prevFmt = os.flags();

   os << std::setw(parw) << " " << "\t";
   for (unsigned int i = 0; i < npar; ++i) {
      if (!IsParameterFixed(i) ) {
         os << std::right  << std::setw(matw)  << GetParameterName(i) ;
      }
   }
   os << std::endl;
   for (unsigned int i = 0; i < npar; ++i) {
      if (!IsParameterFixed(i) ) {
         os << std::left << std::setw(parw) << GetParameterName(i) << "\t";
         for (unsigned int j = 0; j < npar; ++j) {
            if (!IsParameterFixed(j) ) {
               os.precision(kPrec); os.width(kWidth);  os << std::right << std::setw(matw) << CovMatrix(i,j);
            }
         }
         os << std::endl;
      }
   }
//   os << "****************************************\n";
   os << "\nCorrelation Matrix:\n\n";
   os << std::setw(parw) << " " << "\t";
   for (unsigned int i = 0; i < npar; ++i) {
      if (!IsParameterFixed(i) ) {
         os << std::right << std::setw(matw)  << GetParameterName(i) ;
      }
   }
   os << std::endl;
   for (unsigned int i = 0; i < npar; ++i) {
      if (!IsParameterFixed(i) ) {
         os << std::left << std::setw(parw) << std::left << GetParameterName(i) << "\t";
         for (unsigned int j = 0; j < npar; ++j) {
            if (!IsParameterFixed(j) ) {
               os.precision(kPrec); os.width(kWidth);  os << std::right << std::setw(matw) << Correlation(i,j);
            }
         }
         os << std::endl;
      }
   }
   // restore alignment and precision
   os.setf(prevFmt, std::ios::adjustfield);
   os.precision(prevPrec);
}

void FitResult::GetConfidenceIntervals(unsigned int n, unsigned int stride1, unsigned int stride2, const double * x, double * ci, double cl, bool norm ) const {
   // stride1 stride in coordinate  stride2 stride in dimension space
   // i.e. i-th point in k-dimension is x[ stride1 * i + stride2 * k]
   // compute the confidence interval of the fit on the given data points
   // the dimension of the data points must match the dimension of the fit function
   // confidence intervals are returned in array ci

   if (!fFitFunc) {
      // check if model function exists
      MATH_ERROR_MSG("FitResult::GetConfidenceIntervals","Cannot compute Confidence Intervals without fit model function");
      return;
   }
   assert(fFitFunc);

   // use student quantile in case of normalized errors
   double corrFactor = 1;
   if (fChi2 <= 0 || fNdf == 0) norm = false;
   if (norm)
      corrFactor = TMath::StudentQuantile(0.5 + cl/2, fNdf) * std::sqrt( fChi2/fNdf );
   else
      // correction to apply to the errors given a CL different than 1 sigma (cl=0.683) 
      corrFactor = ROOT::Math::normal_quantile(0.5 + cl/2, 1);



   unsigned int ndim = fFitFunc->NDim();
   unsigned int npar = fFitFunc->NPar();

   std::vector<double> xpoint(ndim);
   std::vector<double> grad(npar);
   std::vector<double> vsum(npar);

   // loop on the points
   for (unsigned int ipoint = 0; ipoint < n; ++ipoint) {

      for (unsigned int kdim = 0; kdim < ndim; ++kdim) {
         unsigned int i = ipoint * stride1 + kdim * stride2;
         assert(i < ndim*n);
         xpoint[kdim] = x[i];
      }

      // calculate gradient of fitted function w.r.t the parameters
      ROOT::Math::RichardsonDerivator d;
      for (unsigned int ipar = 0; ipar < npar; ++ipar) {
         if (!IsParameterFixed(ipar)) {
            ROOT::Math::OneDimParamFunctionAdapter<const ROOT::Math::IParamMultiFunction &> fadapter(*fFitFunc,&xpoint.front(),&fParams.front(),ipar);
            d.SetFunction(fadapter);
            // compute step size as a small fraction of the error
            // (see numerical recipes in C 5.7.8)   1.E-5 is ~ (eps)^1/3
            if ( fErrors[ipar] > 0 )
               d.SetStepSize( std::max( fErrors[ipar]*1.E-5, 1.E-15) );
            else
               d.SetStepSize( std::min(std::max(fParams[ipar]*1.E-5, 1.E-15), 0.0001 ) );

            grad[ipar] = d(fParams[ipar] ); // evaluate df/dp
         }
         else
            grad[ipar] = 0.;  // for fixed parameters
      }

      // multiply covariance matrix with gradient
      vsum.assign(npar,0.0);
      for (unsigned int ipar = 0; ipar < npar; ++ipar) {
         for (unsigned int jpar = 0; jpar < npar; ++jpar) {
             vsum[ipar] += CovMatrix(ipar,jpar) * grad[jpar];
         }
      }
      // multiply gradient by vsum
      double r2 = 0;
      for (unsigned int ipar = 0; ipar < npar; ++ipar) {
         r2 += grad[ipar] * vsum[ipar];
      }
      double r = std::sqrt(r2);
      ci[ipoint] = r * corrFactor;
   }
}

void FitResult::GetConfidenceIntervals(const BinData & data, double * ci, double cl, bool norm ) const {
   // implement confidence intervals from a given bin data sets
   // currently copy the data from Bindata.
   // could implement otherwise directly
   unsigned int ndim = data.NDim();
   unsigned int np = data.NPoints();
   std::vector<double> xdata( ndim * np );
   for (unsigned int i = 0; i < np ; ++i) {
      const double * x = data.Coords(i);
      std::vector<double>::iterator itr = xdata.begin()+ ndim * i;
      std::copy(x,x+ndim,itr);
   }
   // points are arraned as x0,y0,z0, ....xN,yN,zN  (stride1=ndim, stride2=1)
   GetConfidenceIntervals(np,ndim,1,&xdata.front(),ci,cl,norm);
}

std::vector<double> FitResult::GetConfidenceIntervals(double cl, bool norm ) const {
   // implement confidence intervals using stored data sets (if can be retrieved from objective function)
   // it works only in case of chi2 or binned likelihood fits
    const BinData * data = FittedBinData();
    std::vector<double> result;
    if (data) {
       result.resize(data->NPoints() );
       GetConfidenceIntervals(*data, result.data(), cl, norm);
    }
    else {
      MATH_ERROR_MSG("FitResult::GetConfidenceIntervals","Cannot compute Confidence Intervals without the fit bin data");
    }
    return result;
}

// const BinData * GetFitBinData() const {
//    // return a pointer to the binned data used in the fit
//    // works only for chi2 or binned likelihood fits
//    // thus when the objective function stored is a Chi2Func or a PoissonLikelihood
//    ROOT::Math::IMultiGenFunction * f = fObjFunc->get();
//    Chi2Function * chi2func = dynamic_cast<Chi2Function*>(f);
//    if (chi2func) return &(chi2func->Data());
//    PoissonLLFunction * pllfunc = dynamic_cast<PoissonLLFunction*>(f);
//    if (pllfunc) return &(pllfunc->Data());
//    Chi2GradFunction * chi2gradfunc = dynamic_cast<Chi2GradFunction*>(f);
//    if (chi2gradfunc) return &(chi2gradfunc->Data());
//    PoissonLLGradFunction * pllgradfunc = dynamic_cast<PoissonLLFunction*>(f);
//    if (pllgradfunc) return &(pllgradfunc->Data());
//    MATH_WARN_MSG("FitResult::GetFitBinData","Cannot retrun fit bin data set if objective function is not of a known type");
//    return nullptr;
// }

const BinData * FitResult::FittedBinData() const {
   return dynamic_cast<const BinData*> ( fFitData.get() );
}

////////////////////////////////////////////////////////////////////////////////
///  Scan parameter ipar between value of xmin and xmax
///  A array for x and y points should be provided

bool FitResult::Scan(unsigned int ipar, unsigned int &npoints, double *pntsx, double *pntsy, double xmin, double xmax)
{
   if (!pntsx || !pntsy || !npoints)
      return false;

   if (!fMinimizer) {
      MATH_ERROR_MSG("FitResult::Scan", "Minimizer is not available - cannot Scan");
      return false;
   }

   return fMinimizer->Scan(ipar, npoints, pntsx, pntsy, xmin, xmax);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 2D contour around the minimum for the parameter ipar and jpar
/// if a minimum does not exist or is invalid it will return false
/// A array for x and y points should be provided
/// Pass optionally the confidence level, default is 0.683
/// it is assumed that ErrorDef() defines the right error definition
/// (i.e 1 sigma error for one parameter). If not the confidence level are scaled to new level

bool FitResult::Contour(unsigned int ipar, unsigned int jpar, unsigned int &npoints, double *pntsx, double *pntsy, double confLevel)
{
   if (!pntsx || !pntsy || !npoints)
      return false;

   if (!fMinimizer) {
      MATH_ERROR_MSG("FitResult::Contour", "Minimizer is not available - cannot produce Contour");
      return false;
   }

   // get error level used for fitting
   double upScale = fMinimizer->ErrorDef();

   double upVal = TMath::ChisquareQuantile(confLevel, 2); // 2 is number of parameter we do the contour

   // set required error definition in minimizer
   fMinimizer->SetErrorDef(upScale * upVal);

   bool ret = fMinimizer->Contour(ipar, jpar, npoints, pntsx, pntsy);

   // restore the error level used for fitting
   fMinimizer->SetErrorDef(upScale);

   return ret;
}

   } // end namespace Fit

} // end namespace ROOT
