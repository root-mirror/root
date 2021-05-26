// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


/** \class RooStats::BayesianCalculator
    \ingroup Roostats

BayesianCalculator is a concrete implementation of IntervalCalculator, providing the computation
of a credible interval using a Bayesian method.
The class works only for one single parameter of interest and it integrates the likelihood function with the given prior
probability density function to compute the posterior probability. The result of the class is a one dimensional interval
(class SimpleInterval ), which is obtained from inverting the cumulative posterior distribution.
This calculator works then only for model with a single parameter of interest.
The model can instead have several nuisance parameters which are integrated (marginalized) in the computation of the posterior function.
The integration and normalization of the posterior is computed using numerical integration methods provided by ROOT.
See the MCMCCalculator for model with multiple parameters of interest.

The interface allows one to construct the class by passing the data set, probability density function for the model, the prior
functions and then the parameter of interest to scan. The nuisance parameters can also be passed to be marginalized when
computing the posterior. Alternatively, the class can be constructed by passing the data and the ModelConfig containing
all the needed information (model pdf, prior pdf, parameter of interest, nuisance parameters, etc..)

After configuring the calculator, one only needs to ask GetInterval(), which
will return an SimpleInterval object. By default the extrem of the integral are obtained by inverting directly the
cumulative posterior distribution. By using the method SetScanOfPosterior(nbins) the interval is then obtained by
scanning  the posterior function in the given number of points. The first method is in general faster but it requires an
integration one extra dimension  ( in the poi in addition to the nuisance parameters), therefore in some case it can be
less robust.

The class can also return the posterior function (method GetPosteriorFunction) or if needed the normalized
posterior function (the posterior pdf) (method GetPosteriorPdf). A posterior plot is also obtained using
the GetPosteriorPlot method.

The class allows to use different integration methods for integrating in (marginalizing) the nuisances and in the poi. All the numerical
integration methods of ROOT can be used via the method SetIntegrationType (see more in the documentation of
this method).

Calculator estimating a credible interval using the Bayesian procedure.
The calculator computes given the model the posterior distribution and estimates the
credible interval from the given function.
*/


// include other header files

#include "RooAbsFunc.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooBrentRootFinder.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"

// include header file of this class
#include "RooStats/BayesianCalculator.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/RooStatsUtils.h"

#include "Math/IFunction.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/Integrator.h"
#include "Math/RootFinder.h"
#include "Math/BrentMinimizer1D.h"
#include "RooFunctor.h"
#include "RooFunctor1DBinding.h"
#include "RooTFnBinding.h"
#include "RooMsgService.h"

#include "TAxis.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"

#include <map>
#include <cmath>

#include "RConfigure.h"

ClassImp(RooStats::BayesianCalculator);

using namespace std;

namespace RooStats {


// first some utility classes and functions

#ifdef R__HAS_MATHMORE
   const ROOT::Math::RootFinder::EType kRootFinderType = ROOT::Math::RootFinder::kGSL_BRENT;
#else
   const ROOT::Math::RootFinder::EType kRootFinderType = ROOT::Math::RootFinder::kBRENT;
#endif




struct  LikelihoodFunction {
   LikelihoodFunction(RooFunctor & f, RooFunctor * prior = 0, double offset = 0) :
      fFunc(f), fPrior(prior),
      fOffset(offset), fMaxL(0) {
      fFunc.binding().resetNumCall();
   }

   void SetPrior(RooFunctor * prior) { fPrior = prior; }

   double operator() (const double *x ) const {
      double nll = fFunc(x) - fOffset;
      double likelihood =  std::exp(-nll);

      if (fPrior) likelihood *= (*fPrior)(x);

      int nCalls = fFunc.binding().numCall();
      if (nCalls > 0 && nCalls % 1000 == 0) {
         ooccoutD((TObject*)0,Eval) << "Likelihood evaluation ncalls = " << nCalls
                                    << " x0 " << x[0] << "  nll = " << nll+fOffset;
         if (fPrior) ooccoutD((TObject*)0,Eval) << " prior(x) = " << (*fPrior)(x);
         ooccoutD((TObject*)0,Eval) << " likelihood " << likelihood
                                    << " max Likelihood " << fMaxL << std::endl;
      }

      if  (likelihood > fMaxL ) {
         fMaxL = likelihood;
         if ( likelihood > 1.E10) {
            ooccoutW((TObject*)0,Eval) << "LikelihoodFunction::()  WARNING - Huge likelihood value found for  parameters ";
            for (int i = 0; i < fFunc.nObs(); ++i)
               ooccoutW((TObject*)0,Eval) << " x[" << i << " ] = " << x[i];
            ooccoutW((TObject*)0,Eval) << "  nll = " << nll << " L = " << likelihood << std::endl;
         }
      }

      return likelihood;
   }

   // for the 1D case
   double operator() (double x) const {
      // just call the previous method
      assert(fFunc.nObs() == 1); // check nobs = 1
      double tmp = x;
      return (*this)(&tmp);
   }

   RooFunctor & fFunc;     // functor representing the nll function
   RooFunctor * fPrior;     // functor representing the prior function
   double fOffset;         //  offset used to bring the nll in a reasonable range for computing the exponent
   mutable double fMaxL;
};


// Posterior CDF Function class
// for integral of posterior function in nuisance and POI
// 1-Dim function as function of the poi

class PosteriorCdfFunction : public ROOT::Math::IGenFunction {

public:

   PosteriorCdfFunction(RooAbsReal & nll,  RooArgList & bindParams, RooAbsReal * prior = 0, const char * integType = 0, double nllMinimum = 0) :
      fFunctor(nll, bindParams, RooArgList() ),              // functor
      fPriorFunc(nullptr),
      fLikelihood(fFunctor, 0, nllMinimum),         // integral of exp(-nll) function
      fIntegrator(ROOT::Math::IntegratorMultiDim::GetType(integType) ),  // integrator
      fXmin(bindParams.getSize() ),               // vector of parameters (min values)
      fXmax(bindParams.getSize() ),               // vector of parameter (max values)
      fNorm(1.0), fNormErr(0.0), fOffset(0), fMaxPOI(0),
      fHasNorm(false),  fUseOldValues(true), fError(false)
   {

      if (prior) {
         fPriorFunc = std::make_shared<RooFunctor>(*prior, bindParams, RooArgList());
         fLikelihood.SetPrior(fPriorFunc.get() );
      }

      fIntegrator.SetFunction(fLikelihood, bindParams.getSize() );

      ooccoutD((TObject*)0,NumIntegration) << "PosteriorCdfFunction::Compute integral of posterior in nuisance and poi. "
                                           << " nllMinimum is " << nllMinimum << std::endl;

      std::vector<double> par(bindParams.getSize());
      for (unsigned int i = 0; i < fXmin.size(); ++i) {
         RooRealVar & var = (RooRealVar &) bindParams[i];
         fXmin[i] = var.getMin();
         fXmax[i] = var.getMax();
         par[i] = var.getVal();
         ooccoutD((TObject*)0,NumIntegration) << "PosteriorFunction::Integrate" << var.GetName()
                                              << " in interval [ " <<  fXmin[i] << " , " << fXmax[i] << " ] " << std::endl;
      }

      fIntegrator.Options().Print(ooccoutD((TObject*)0,NumIntegration));

      // store max POI value because it will be changed when evaluating the function
      fMaxPOI = fXmax[0];

      // compute first the normalization with  the poi
      fNorm = (*this)( fMaxPOI );
      if (fError) ooccoutE((TObject*)0,NumIntegration) << "PosteriorFunction::Error computing normalization - norm = " << fNorm << std::endl;
      fHasNorm = true;
      fNormCdfValues.insert(std::make_pair(fXmin[0], 0) );
      fNormCdfValues.insert(std::make_pair(fXmax[0], 1.0) );

   }

   // copy constructor (needed for Cloning the object)
   // need special treatment because integrator
   // has no copy constructor
   PosteriorCdfFunction(const PosteriorCdfFunction & rhs) :
      ROOT::Math::IGenFunction(),
      fFunctor(rhs.fFunctor),
      //fPriorFunc(std::shared_ptr<RooFunctor>((RooFunctor*)0)),
      fPriorFunc(rhs.fPriorFunc),
      fLikelihood(fFunctor, fPriorFunc.get(), rhs.fLikelihood.fOffset),
      fIntegrator(ROOT::Math::IntegratorMultiDim::GetType( rhs.fIntegrator.Name().c_str() ) ),  // integrator
      fXmin( rhs.fXmin),
      fXmax( rhs.fXmax),
      fNorm( rhs.fNorm),
      fNormErr( rhs.fNormErr),
      fOffset(rhs.fOffset),
      fMaxPOI(rhs.fMaxPOI),
      fHasNorm(rhs.fHasNorm),
      fUseOldValues(rhs.fUseOldValues),
      fError(rhs.fError),
      fNormCdfValues(rhs.fNormCdfValues)
   {
      fIntegrator.SetFunction(fLikelihood, fXmin.size() );
      // need special treatment for the smart pointer
      // if (rhs.fPriorFunc.get() ) {
      //    fPriorFunc = std::shared_ptr<RooFunctor>(new RooFunctor(*(rhs.fPriorFunc) ) );
      //    fLikelihood.SetPrior( fPriorFunc.get() );
      // }
   }


   bool HasError() const { return fError; }


   ROOT::Math::IGenFunction * Clone() const {
      ooccoutD((TObject*)0,NumIntegration) << " cloning function .........." << std::endl;
      return new PosteriorCdfFunction(*this);
   }

   // offset value for computing the root
   void SetOffset(double offset) { fOffset = offset; }

private:

   // make assignment operator private
   PosteriorCdfFunction& operator=(const PosteriorCdfFunction &) {
      return *this;
   }

   double DoEval (double x) const {

      // evaluate cdf at poi value x by integrating poi from [xmin,x] and all the nuisances
      fXmax[0] = x;
      if (x <= fXmin[0] ) return -fOffset;
      // could also avoid a function evaluation at maximum
      if (x >= fMaxPOI && fHasNorm) return 1. - fOffset;  // cdf is bound to these values

      // computes the integral using a previous cdf estimate
      double  normcdf0 = 0;
      if (fHasNorm && fUseOldValues) {
         // look in the map of the stored cdf values the closes one
         std::map<double,double>::iterator itr = fNormCdfValues.upper_bound(x);
         --itr;   // upper bound returns a position 1 up of the value we want
         if (itr != fNormCdfValues.end() ) {
            fXmin[0] = itr->first;
            normcdf0 = itr->second;
            // ooccoutD((TObject*)0,NumIntegration) << "PosteriorCdfFunction:   computing integral between in poi interval : "
            //                                      << fXmin[0] << " -  " << fXmax[0] << std::endl;
         }
      }

      fFunctor.binding().resetNumCall();  // reset number of calls for debug

      double cdf = fIntegrator.Integral(&fXmin[0],&fXmax[0]);
      double error = fIntegrator.Error();
      double normcdf =  cdf/fNorm;  // normalize the cdf

      ooccoutD((TObject*)0,NumIntegration) << "PosteriorCdfFunction: poi = [" << fXmin[0] << " , "
                                           << fXmax[0] << "] integral =  " << cdf << " +/- " << error
                                           << "  norm-integ = " << normcdf << " cdf(x) = " << normcdf+normcdf0
                                           << " ncalls = " << fFunctor.binding().numCall() << std::endl;

      if (TMath::IsNaN(cdf) || cdf > std::numeric_limits<double>::max()) {
         ooccoutE((TObject*)0,NumIntegration) << "PosteriorFunction::Error computing integral - cdf = "
                                              << cdf << std::endl;
         fError = true;
      }

      if (cdf != 0 && error/cdf > 0.2 )
         oocoutW((TObject*)0,NumIntegration) << "PosteriorCdfFunction: integration error  is larger than 20 %   x0 = " << fXmin[0]
                                              << " x = " << x << " cdf(x) = " << cdf << " +/- " << error << std::endl;

      if (!fHasNorm) {
         oocoutI((TObject*)0,NumIntegration) << "PosteriorCdfFunction - integral of posterior = "
                                             << cdf << " +/- " << error << std::endl;
         fNormErr = error;
         return cdf;
      }

      normcdf += normcdf0;

      // store values in the map
      if (fUseOldValues) {
         fNormCdfValues.insert(std::make_pair(x, normcdf) );
      }

      double errnorm = sqrt( error*error + normcdf*normcdf * fNormErr * fNormErr )/fNorm;
      if (normcdf > 1. + 3 * errnorm) {
         oocoutW((TObject*)0,NumIntegration) << "PosteriorCdfFunction: normalized cdf values is larger than 1"
                                              << " x = " << x << " normcdf(x) = " << normcdf << " +/- " << error/fNorm << std::endl;
      }

      return normcdf - fOffset;  // apply an offset (for finding the roots)
   }

   mutable RooFunctor fFunctor;                   // functor binding nll
   mutable std::shared_ptr<RooFunctor> fPriorFunc;  // functor binding the prior
   LikelihoodFunction fLikelihood;               // likelihood function
   mutable ROOT::Math::IntegratorMultiDim  fIntegrator; // integrator  (mutable because Integral() is not const
   mutable std::vector<double> fXmin;    // min value of parameters (poi+nuis) -
   mutable std::vector<double> fXmax;   // max value of parameters (poi+nuis) - max poi changes so it is mutable
   double fNorm;      // normalization value (computed in ctor)
   mutable double fNormErr;    // normalization error value (computed in ctor)
   double fOffset;   // offset for computing the root
   double fMaxPOI;  // maximum value of POI
   bool fHasNorm; // flag to control first call to the function
   bool fUseOldValues;  // use old cdf values
   mutable bool fError;     // flag to indicate if a numerical evaluation error occurred
   mutable std::map<double,double> fNormCdfValues;
};

//__________________________________________________________________
// Posterior Function class
// 1-Dim function as function of the poi
// and it integrated all the nuisance parameters

class PosteriorFunction : public ROOT::Math::IGenFunction {

public:


   PosteriorFunction(RooAbsReal & nll, RooRealVar & poi, RooArgList & nuisParams, RooAbsReal * prior = 0, const char * integType = 0, double
                     norm = 1.0,  double nllOffset = 0, int niter = 0) :
      fFunctor(nll, nuisParams, RooArgList() ),
      fPriorFunc(nullptr),
      fLikelihood(fFunctor, 0, nllOffset),
      fPoi(&poi),
      fXmin(nuisParams.getSize() ),
      fXmax(nuisParams.getSize() ),
      fNorm(norm),
      fError(0)
   {

      if (prior) {
         fPriorFunc = std::make_shared<RooFunctor>(*prior, nuisParams, RooArgList());
         fLikelihood.SetPrior(fPriorFunc.get() );
      }

      ooccoutD((TObject*)0,NumIntegration) << "PosteriorFunction::Evaluate the posterior function by integrating the nuisances: " << std::endl;
      for (unsigned int i = 0; i < fXmin.size(); ++i) {
         RooRealVar & var = (RooRealVar &) nuisParams[i];
         fXmin[i] = var.getMin();
         fXmax[i] = var.getMax();
         ooccoutD((TObject*)0,NumIntegration) << "PosteriorFunction::Integrate " << var.GetName()
                                              << " in interval [" <<  fXmin[i] << " , " << fXmax[i] << " ] " << std::endl;
      }
      if (fXmin.size() == 1) { // 1D case
         fIntegratorOneDim.reset( new ROOT::Math::Integrator(ROOT::Math::IntegratorOneDim::GetType(integType) ) );

         fIntegratorOneDim->SetFunction(fLikelihood);
         // interested only in relative tolerance
         //fIntegratorOneDim->SetAbsTolerance(1.E-300);
         fIntegratorOneDim->Options().Print(ooccoutD((TObject*)0,NumIntegration) );
      }
      else if (fXmin.size() > 1) { // multiDim case
         fIntegratorMultiDim.reset(new ROOT::Math::IntegratorMultiDim(ROOT::Math::IntegratorMultiDim::GetType(integType) ) );
         fIntegratorMultiDim->SetFunction(fLikelihood, fXmin.size());
         ROOT::Math::IntegratorMultiDimOptions opt = fIntegratorMultiDim->Options();
         if (niter > 0) {
            opt.SetNCalls(niter);
            fIntegratorMultiDim->SetOptions(opt);
         }
         //fIntegratorMultiDim->SetAbsTolerance(1.E-300);
         // print the options
         opt.Print(ooccoutD((TObject*)0,NumIntegration) );
      }
   }


   ROOT::Math::IGenFunction * Clone() const {
      assert(1);
      return 0; // cannot clone this function for integrator
   }

   double Error() const { return fError;}


private:
   double DoEval (double x) const {

      // evaluate posterior function at a poi value x by integrating all nuisance parameters

      fPoi->setVal(x);
      fFunctor.binding().resetNumCall();  // reset number of calls for debug

      double f = 0;
      double error = 0;
      if (fXmin.size() == 1) { // 1D case
         f = fIntegratorOneDim->Integral(fXmin[0],fXmax[0]);
         error = fIntegratorOneDim->Error();
      }
      else if (fXmin.size() > 1) { // multi-dim case
         f = fIntegratorMultiDim->Integral(&fXmin[0],&fXmax[0]);
         error = fIntegratorMultiDim->Error();
      } else {
         // no integration to be done
         f = fLikelihood(x);
      }

      // debug
      ooccoutD((TObject*)0,NumIntegration) << "PosteriorFunction:  POI value  =  "
                                           << x << "\tf(x) =  " << f << " +/- " << error
                                           << "  norm-f(x) = " << f/fNorm
                                           << " ncalls = " << fFunctor.binding().numCall() << std::endl;




      if (f != 0 && error/f > 0.2 )
         ooccoutW((TObject*)0,NumIntegration) << "PosteriorFunction::DoEval - Error from integration in "
                                              << fXmin.size() <<  " Dim is larger than 20 % "
                                              << "x = " << x << " p(x) = " << f << " +/- " << error << std::endl;

      fError = error / fNorm;
      return f / fNorm;
   }

   mutable RooFunctor fFunctor;
   mutable std::shared_ptr<RooFunctor> fPriorFunc;  // functor binding the prior
   LikelihoodFunction fLikelihood;
   RooRealVar * fPoi;
   std::unique_ptr<ROOT::Math::Integrator>  fIntegratorOneDim;
   std::unique_ptr<ROOT::Math::IntegratorMultiDim>  fIntegratorMultiDim;
   std::vector<double> fXmin;
   std::vector<double> fXmax;
   double fNorm;
   mutable double fError;
};

////////////////////////////////////////////////////////////////////////////////
/// Posterior function obtaining sampling  toy MC for the nuisance according to their pdf

class PosteriorFunctionFromToyMC : public ROOT::Math::IGenFunction {

public:


   PosteriorFunctionFromToyMC(RooAbsReal & nll, RooAbsPdf & pdf, RooRealVar & poi, RooArgList & nuisParams, RooAbsReal * prior = 0, double
                              nllOffset = 0, int niter = 0, bool redoToys = true ) :
      fFunctor(nll, nuisParams, RooArgList() ),
      fPriorFunc(nullptr),
      fLikelihood(fFunctor, 0, nllOffset),
      fPdf(&pdf),
      fPoi(&poi),
      fNuisParams(nuisParams),
      fGenParams(0),
      fNumIterations(niter),
      fError(-1),
      fRedoToys(redoToys)
   {
      if (niter == 0) fNumIterations = 100; // default value

      if (prior) {
         fPriorFunc = std::make_shared<RooFunctor>(*prior, nuisParams, RooArgList());
         fLikelihood.SetPrior(fPriorFunc.get() );
      }

      ooccoutI((TObject*)0,InputArguments) << "PosteriorFunctionFromToyMC::Evaluate the posterior function by randomizing the nuisances:  niter " << fNumIterations << std::endl;

      ooccoutI((TObject*)0,InputArguments) << "PosteriorFunctionFromToyMC::Pdf used for randomizing the nuisance is " << fPdf->GetName() << std::endl;
      // check that pdf contains  the nuisance
      RooArgSet * vars = fPdf->getVariables();
      for (int i = 0; i < fNuisParams.getSize(); ++i) {
         if (!vars->find( fNuisParams[i].GetName() ) ) {
            ooccoutW((TObject*)0,InputArguments) << "Nuisance parameter " << fNuisParams[i].GetName()
                                                 << " is not part of sampling pdf. "
                                                 << "they will be treated as constant " << std::endl;
         }
      }
      delete vars;

      if (!fRedoToys) {
         ooccoutI((TObject*)0,InputArguments) << "PosteriorFunctionFromToyMC::Generate nuisance toys only one time (for all POI points)" << std::endl;
         GenerateToys();
      }
   }

   virtual ~PosteriorFunctionFromToyMC() { if (fGenParams) delete fGenParams; }

   // generate first n-samples of the nuisance parameters
   void GenerateToys() const {
      if (fGenParams) delete fGenParams;
      fGenParams = fPdf->generate(fNuisParams, fNumIterations);
      if(fGenParams==0) {
         ooccoutE((TObject*)0,InputArguments) << "PosteriorFunctionFromToyMC - failed to generate nuisance parameters" << std::endl;
      }
   }

   double Error() const { return fError;}

   ROOT::Math::IGenFunction * Clone() const {
      // use default copy constructor
      //return new PosteriorFunctionFromToyMC(*this);
      //  clone not implemented
      assert(1);
      return 0;
   }

private:
   // evaluate the posterior at the poi value x
   double DoEval( double x) const {

      int npar = fNuisParams.getSize();
      assert (npar > 0);


      // generate the toys
      if (fRedoToys) GenerateToys();
      if (!fGenParams) return 0;

      // evaluate posterior function at a poi value x by integrating all nuisance parameters

      fPoi->setVal(x);

      // loop over all of the generate data
      double sum = 0;
      double sum2 = 0;

      for(int iter=0; iter<fNumIterations; ++iter) {

         // get the set of generated parameters and set the nuisance parameters to the generated values
         std::vector<double> p(npar);
         for (int i = 0; i < npar; ++i) {
            const RooArgSet* genset=fGenParams->get(iter);
            RooAbsArg * arg = genset->find( fNuisParams[i].GetName() );
            RooRealVar * var = dynamic_cast<RooRealVar*>(arg);
            assert(var != 0);
            p[i] = var->getVal();
            ((RooRealVar &) fNuisParams[i]).setVal(p[i]);
         }

         // evaluate now the likelihood function
         double fval =  fLikelihood( &p.front() );

         // likelihood already must contained the pdf we have sampled
         // so we must divided by it. The value must be normalized on all
         // other parameters
         RooArgSet arg(fNuisParams);
         double nuisPdfVal = fPdf->getVal(&arg);
         fval /= nuisPdfVal;


         if( fval > std::numeric_limits<double>::max()  ) {
            ooccoutE((TObject*)0,Eval) <<  "BayesianCalculator::EvalPosteriorFunctionFromToy : "
                        << "Likelihood evaluates to infinity " << std::endl;
            ooccoutE((TObject*)0,Eval) <<  "poi value =  " << x << std::endl;
            ooccoutE((TObject*)0,Eval) <<  "Nuisance  parameter values :  ";
            for (int i = 0; i < npar; ++i)
               ooccoutE((TObject*)0,Eval) << fNuisParams[i].GetName() << " = " << p[i] << " ";
            ooccoutE((TObject*)0,Eval) <<  " - return 0   " << std::endl;

            fError = 1.E30;
            return 0;
         }
         if(  TMath::IsNaN(fval) ) {
            ooccoutE((TObject*)0,Eval) <<  "BayesianCalculator::EvalPosteriorFunctionFromToy : "
                        << "Likelihood is a NaN " << std::endl;
            ooccoutE((TObject*)0,Eval) <<  "poi value =  " << x << std::endl;
            ooccoutE((TObject*)0,Eval) <<  "Nuisance  parameter values :  ";
            for (int i = 0; i < npar; ++i)
               ooccoutE((TObject*)0,Eval) << fNuisParams[i].GetName() << " = " << p[i] << " ";
            ooccoutE((TObject*)0,Eval) <<  " - return 0   " << std::endl;
            fError = 1.E30;
            return 0;
         }



         sum += fval;
         sum2 += fval*fval;
      }

      // compute the average and variance
      double val = sum/double(fNumIterations);
      double dval2 = std::max( sum2/double(fNumIterations) - val*val, 0.0);
      fError = std::sqrt( dval2 / fNumIterations);

      // debug
      ooccoutD((TObject*)0,NumIntegration) << "PosteriorFunctionFromToyMC:  POI value  =  "
                                           << x << "\tp(x) =  " << val << " +/- " << fError << std::endl;


      if (val != 0 && fError/val > 0.2 ) {
         ooccoutW((TObject*)0,NumIntegration) << "PosteriorFunctionFromToyMC::DoEval"
                                              << " - Error in estimating posterior is larger than 20% ! "
                                              << "x = " << x << " p(x) = " << val << " +/- " << fError << std::endl;
      }


      return val;
   }

   mutable RooFunctor fFunctor;
   mutable std::shared_ptr<RooFunctor> fPriorFunc;  // functor binding the prior
   LikelihoodFunction fLikelihood;
   mutable RooAbsPdf * fPdf;
   RooRealVar * fPoi;
   RooArgList fNuisParams;
   mutable RooDataSet * fGenParams;
   int fNumIterations;
   mutable double fError;
   bool fRedoToys;                    // do toys every iteration

};

////////////////////////////////////////////////////////////////////////////////
// Implementation of BayesianCalculator
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// default constructor

BayesianCalculator::BayesianCalculator() :
   fData(0),
   fPdf(0),
   fPriorPdf(0),
   fNuisancePdf(0),
   fProductPdf (0), fLogLike(0), fLikelihood (0), fIntegratedLikelihood (0), fPosteriorPdf(0),
   fPosteriorFunction(0), fApproxPosterior(0),
   fLower(0), fUpper(0),
   fNLLMin(0),
   fSize(0.05), fLeftSideFraction(0.5),
   fBrfPrecision(0.00005),
   fNScanBins(-1),
   fNumIterations(0),
   fValidInterval(false)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from data set, model pdf, parameter of interests and prior pdf
/// If nuisance parameters are given they will be integrated according either to the prior or
/// their constraint term included in the model

BayesianCalculator::BayesianCalculator( /* const char* name,  const char* title, */
                      RooAbsData& data,
                                                    RooAbsPdf& pdf,
                      const RooArgSet& POI,
                      RooAbsPdf& priorPdf,
                      const RooArgSet* nuisanceParameters ) :
   fData(&data),
   fPdf(&pdf),
   fPOI(POI),
   fPriorPdf(&priorPdf),
   fNuisancePdf(0),
   fProductPdf (0), fLogLike(0), fLikelihood (0), fIntegratedLikelihood (0), fPosteriorPdf(0),
   fPosteriorFunction(0), fApproxPosterior(0),
   fLower(0), fUpper(0),
   fNLLMin(0),
   fSize(0.05), fLeftSideFraction(0.5),
   fBrfPrecision(0.00005),
   fNScanBins(-1),
   fNumIterations(0),
   fValidInterval(false)
{

   if (nuisanceParameters) fNuisanceParameters.add(*nuisanceParameters);
   // remove constant nuisance parameters
   RemoveConstantParameters(&fNuisanceParameters);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a data set and a ModelConfig
/// model pdf, poi and nuisances will be taken from the ModelConfig

BayesianCalculator::BayesianCalculator( RooAbsData& data,
                       ModelConfig & model) :
   fData(&data),
   fPdf(model.GetPdf()),
   fPriorPdf( model.GetPriorPdf()),
   fNuisancePdf(0),
   fProductPdf (0), fLogLike(0), fLikelihood (0), fIntegratedLikelihood (0), fPosteriorPdf(0),
   fPosteriorFunction(0), fApproxPosterior(0),
   fLower(0), fUpper(0),
   fNLLMin(0),
   fSize(0.05), fLeftSideFraction(0.5),
   fBrfPrecision(0.00005),
   fNScanBins(-1),
   fNumIterations(0),
   fValidInterval(false)
{
   SetModel(model);
}


BayesianCalculator::~BayesianCalculator()
{
   // destructor
   ClearAll();
}

////////////////////////////////////////////////////////////////////////////////
/// clear all cached pdf objects

void BayesianCalculator::ClearAll() const {
   if (fProductPdf) delete fProductPdf;
   if (fLogLike) delete fLogLike;
   if (fLikelihood) delete fLikelihood;
   if (fIntegratedLikelihood) delete fIntegratedLikelihood;
   if (fPosteriorPdf) delete fPosteriorPdf;
   if (fPosteriorFunction) delete fPosteriorFunction;
   if (fApproxPosterior) delete fApproxPosterior;
   fPosteriorPdf = 0;
   fPosteriorFunction = 0;
   fProductPdf = 0;
   fLogLike = 0;
   fLikelihood = 0;
   fIntegratedLikelihood = 0;
   fLower = 0;
   fUpper = 0;
   fNLLMin = 0;
   fValidInterval = false;
}

////////////////////////////////////////////////////////////////////////////////
/// set the model to use
/// The model pdf, prior pdf, parameter of interest and nuisances
/// will be taken according to the model

void BayesianCalculator::SetModel(const ModelConfig & model) {

   fPdf = model.GetPdf();
   fPriorPdf =  model.GetPriorPdf();
   // assignment operator = does not do a real copy the sets (must use add method)
   fPOI.removeAll();
   fNuisanceParameters.removeAll();
   fConditionalObs.removeAll();
   fGlobalObs.removeAll();
   if (model.GetParametersOfInterest()) fPOI.add( *(model.GetParametersOfInterest()) );
   if (model.GetNuisanceParameters())  fNuisanceParameters.add( *(model.GetNuisanceParameters() ) );
   if (model.GetConditionalObservables())  fConditionalObs.add( *(model.GetConditionalObservables() ) );
   if (model.GetGlobalObservables())  fGlobalObs.add( *(model.GetGlobalObservables() ) );
   // remove constant nuisance parameters
   RemoveConstantParameters(&fNuisanceParameters);

   // invalidate the cached pointers
   ClearAll();
}

////////////////////////////////////////////////////////////////////////////////
/// Build and return the posterior function (not normalized) as a RooAbsReal
/// the posterior is obtained from the product of the likelihood function and the
/// prior pdf which is then integrated in the nuisance parameters (if existing).
/// A prior function for the nuisance can be specified either in the prior pdf object
/// or in the model itself. If no prior nuisance is specified, but prior parameters are then
/// the integration is performed assuming a flat prior for the nuisance parameters.
///
/// NOTE: the return object is managed by the BayesianCalculator class, users do not need to delete it,
///       but the object will be deleted when the BayesiabCalculator object is deleted

RooAbsReal* BayesianCalculator::GetPosteriorFunction() const
{

   if (fIntegratedLikelihood) return fIntegratedLikelihood;
   if (fLikelihood) return fLikelihood;

   // run some sanity checks
   if (!fPdf ) {
      coutE(InputArguments) << "BayesianCalculator::GetPosteriorPdf - missing pdf model" << std::endl;
      return 0;
   }
   if (fPOI.getSize() == 0) {
      coutE(InputArguments) << "BayesianCalculator::GetPosteriorPdf - missing parameter of interest" << std::endl;
      return 0;
   }
   if (fPOI.getSize() > 1) {
      coutE(InputArguments) << "BayesianCalculator::GetPosteriorPdf - current implementation works only on 1D intervals" << std::endl;
      return 0;
   }


   RooArgSet* constrainedParams = fPdf->getParameters(*fData);
   // remove the constant parameters
   RemoveConstantParameters(constrainedParams);

   //constrainedParams->Print("V");

   // use RooFit::Constrain() to be sure constraints terms are taken into account
   fLogLike = fPdf->createNLL(*fData, RooFit::Constrain(*constrainedParams), RooFit::ConditionalObservables(fConditionalObs), RooFit::GlobalObservables(fGlobalObs) );



   ccoutD(Eval) <<  "BayesianCalculator::GetPosteriorFunction : "
                << " pdf value " <<  fPdf->getVal()
                << " neglogLikelihood = " << fLogLike->getVal() << std::endl;

   if (fPriorPdf)
      ccoutD(Eval)  << "\t\t\t priorPOI value " << fPriorPdf->getVal() << std::endl;

   // check that likelihood evaluation is not infinity
   double nllVal =  fLogLike->getVal();
   if ( nllVal > std::numeric_limits<double>::max() ) {
      coutE(Eval) <<  "BayesianCalculator::GetPosteriorFunction : "
                  << " Negative log likelihood evaluates to infinity " << std::endl
                  << " Non-const Parameter values : ";
      RooArgList p(*constrainedParams);
      for (int i = 0; i < p.getSize(); ++i) {
         RooRealVar * v = dynamic_cast<RooRealVar *>(&p[i] );
         if (v!=0) ccoutE(Eval) << v->GetName() << " = " << v->getVal() << "   ";
      }
      ccoutE(Eval) << std::endl;
      ccoutE(Eval) << "--  Perform a full likelihood fit of the model before or set more reasonable parameter values"
                   << std::endl;
      coutE(Eval) << "BayesianCalculator::GetPosteriorFunction : " << " cannot compute posterior function "  << std::endl;
      return 0;
   }



   // need do find minimum of log-likelihood in the range to shift function
   // to avoid numerical errors when we compute the likelihood (overflows in the exponent)
   // N.B.: this works for only 1 parameter of interest otherwise Minuit should be used for finding the minimum
   RooFunctor * nllFunc = fLogLike->functor(fPOI);
   assert(nllFunc);
   ROOT::Math::Functor1D wnllFunc(*nllFunc);
   RooRealVar* poi = dynamic_cast<RooRealVar*>( fPOI.first() );
   assert(poi);

   // try to reduce some error messages
   //Bool_t silentMode = (RooMsgService::instance().globalKillBelow() >= RooFit::ERROR || RooMsgService::instance().silentMode()) ;
   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);


   coutI(Eval) <<  "BayesianCalculator::GetPosteriorFunction : "
               << " nll value " <<  nllVal << " poi value = " << poi->getVal() << std::endl;


   ROOT::Math::BrentMinimizer1D minim;
   minim.SetFunction(wnllFunc,poi->getMin(),poi->getMax() );
   bool ret  = minim.Minimize(100,1.E-3,1.E-3);
   fNLLMin = 0;
   if (ret) fNLLMin = minim.FValMinimum();

   coutI(Eval) << "BayesianCalculator::GetPosteriorFunction : minimum of NLL vs POI for POI =  "
          << poi->getVal() << " min NLL = " << fNLLMin << std::endl;

   delete nllFunc;

   delete constrainedParams;


   if ( fNuisanceParameters.getSize() == 0 ||  fIntegrationType.Contains("ROOFIT") ) {

      ccoutD(Eval) << "BayesianCalculator::GetPosteriorFunction : use ROOFIT integration  "
                   << std::endl;

#ifdef DOLATER // (not clear why this does not work)
      // need to make in this case a likelihood from the nll and make the product with the prior
      TString likeName = TString("likelihood_times_prior_") + TString(fPriorPdf->GetName());
      TString formula;
      formula.Form("exp(-@0+%f+log(@1))",fNLLMin);
      fLikelihood = new RooFormulaVar(likeName,formula,RooArgList(*fLogLike,*fPriorPdf));
#else
      // here use RooProdPdf (not very nice) but working

      if (fLogLike) delete fLogLike;
      if (fProductPdf) {
         delete fProductPdf;
         fProductPdf = 0;
      }

      // // create a unique name for the product pdf
      RooAbsPdf *  pdfAndPrior = fPdf;
      if (fPriorPdf) {
         TString prodName = TString("product_") + TString(fPdf->GetName()) + TString("_") + TString(fPriorPdf->GetName() );
         // save this as data member since it needs to be deleted afterwards
         fProductPdf = new RooProdPdf(prodName,"",RooArgList(*fPdf,*fPriorPdf));
         pdfAndPrior = fProductPdf;
      }

      RooArgSet* constrParams = fPdf->getParameters(*fData);
      // remove the constant parameters
      RemoveConstantParameters(constrParams);
      fLogLike = pdfAndPrior->createNLL(*fData, RooFit::Constrain(*constrParams),RooFit::ConditionalObservables(fConditionalObs),RooFit::GlobalObservables(fGlobalObs) );
      delete constrParams;

      TString likeName = TString("likelihood_times_prior_") + TString(pdfAndPrior->GetName());
      TString formula;
      formula.Form("exp(-@0+%f)",fNLLMin);
      fLikelihood = new RooFormulaVar(likeName,formula,RooArgList(*fLogLike));
#endif


      // if no nuisance parameter we can just return the likelihood function
      if (fNuisanceParameters.getSize() == 0) {
         fIntegratedLikelihood = fLikelihood;
         fLikelihood = 0;
      }
      else
         // case of using RooFit for the integration
         fIntegratedLikelihood = fLikelihood->createIntegral(fNuisanceParameters);


   }

   else if ( fIntegrationType.Contains("TOYMC") ) {
      // compute the posterior as expectation values of the likelihood function
      // sampling on the nuisance parameters

      RooArgList nuisParams(fNuisanceParameters);

      bool doToysEveryIteration = true;
      // if type is 1-TOYMC or TOYMC-1
      if ( fIntegrationType.Contains("1") || fIntegrationType.Contains("ONE")  ) doToysEveryIteration = false;

      RooAbsPdf * samplingPdf = (fNuisancePdf) ? fNuisancePdf : fPdf;
      if (!fNuisancePdf) {
         ccoutI(Eval) << "BayesianCalculator::GetPosteriorFunction : no nuisance pdf is provided, try using global pdf (this will be slower)"
                      << std::endl;
      }
      fPosteriorFunction = new PosteriorFunctionFromToyMC(*fLogLike, *samplingPdf, *poi, nuisParams, fPriorPdf, fNLLMin,
                                                          fNumIterations, doToysEveryIteration );

      TString name = "toyposteriorfunction_from_";
      name += fLogLike->GetName();
      fIntegratedLikelihood = new RooFunctor1DBinding(name,name,*fPosteriorFunction,*poi);

      // need to scan likelihood in this case
      if (fNScanBins <= 0) fNScanBins = 100;

   }

   else  {

      // use ROOT integration method if there are nuisance parameters

      RooArgList nuisParams(fNuisanceParameters);
      fPosteriorFunction = new PosteriorFunction(*fLogLike, *poi, nuisParams, fPriorPdf, fIntegrationType, 1., fNLLMin, fNumIterations );

      TString name = "posteriorfunction_from_";
      name += fLogLike->GetName();
      fIntegratedLikelihood = new RooFunctor1DBinding(name,name,*fPosteriorFunction,*poi);

   }


   if (RooAbsReal::numEvalErrors() > 0)
      coutW(Eval) << "BayesianCalculator::GetPosteriorFunction : " << RooAbsReal::numEvalErrors() << " errors reported in evaluating log-likelihood function "
                   << std::endl;
   RooAbsReal::clearEvalErrorLog();
   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);

   return fIntegratedLikelihood;

}

////////////////////////////////////////////////////////////////////////////////
/// Build and return the posterior pdf (i.e posterior function normalized to all range of poi)
/// Note that an extra integration in the POI is required for the normalization
/// NOTE: user must delete the returned object

RooAbsPdf* BayesianCalculator::GetPosteriorPdf() const
{

   RooAbsReal * plike = GetPosteriorFunction();
   if (!plike) return 0;


   // create a unique name on the posterior from the names of the components
   TString posteriorName = this->GetName() + TString("_posteriorPdf_") + plike->GetName();

   RooAbsPdf * posteriorPdf = new RooGenericPdf(posteriorName,"@0",*plike);

   return posteriorPdf;
}

////////////////////////////////////////////////////////////////////////////////
/// When an approximate posterior is computed binning the parameter of interest (poi) range
/// (see SetScanOfPosteriors) an histogram is created and can be returned to the user
///  A nullptr is instead returned when the posterior is computed without binning the poi.
///
/// NOTE: the returned object is managed by the BayesianCalculator class,
///  if the user wants to take ownership of the returned histogram, he needs to clone
///  or copy the return object.

TH1 *  BayesianCalculator::GetPosteriorHistogram() const
{
   return  (fApproxPosterior) ? fApproxPosterior->GetHistogram() : nullptr;
}


////////////////////////////////////////////////////////////////////////////////
/// return a RooPlot with the posterior  and the credibility region
/// NOTE: User takes ownership of the returned object

RooPlot* BayesianCalculator::GetPosteriorPlot(bool norm, double precision ) const
{

   GetPosteriorFunction();

   // if a scan is requested approximate the posterior
   if (fNScanBins > 0)
      ApproximatePosterior();

   RooAbsReal * posterior = fIntegratedLikelihood;
   if (norm) {
      // delete and re-do always posterior pdf (could be invalid after approximating it)
      if (fPosteriorPdf) delete fPosteriorPdf;
      fPosteriorPdf = GetPosteriorPdf();
      posterior = fPosteriorPdf;
   }
   if (!posterior) return 0;

   if (!fValidInterval) GetInterval();

   RooAbsRealLValue* poi = dynamic_cast<RooAbsRealLValue*>( fPOI.first() );
   assert(poi);


   RooPlot* plot = poi->frame();
   if (!plot) return 0;

   // try to reduce some error messages
   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);

   plot->SetTitle(TString("Posterior probability of parameter \"")+TString(poi->GetName())+TString("\""));
   posterior->plotOn(plot,RooFit::Range(fLower,fUpper,kFALSE),RooFit::VLines(),RooFit::DrawOption("F"),RooFit::MoveToBack(),RooFit::FillColor(kGray),RooFit::Precision(precision));
   posterior->plotOn(plot);
   plot->GetYaxis()->SetTitle("posterior function");

   // reset the counts and default mode
   RooAbsReal::clearEvalErrorLog();
   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);

   return plot;
}

////////////////////////////////////////////////////////////////////////////////
/// set the integration type (possible type are) :
///
///  - 1D integration ( used when only one nuisance and when the posterior is scanned):
///    adaptive , gauss, nonadaptive
///  -  multidim:
///      - ADAPTIVE,   adaptive numerical integration
///                    The parameter numIters (settable with SetNumIters) is  the max number of function calls.
///                    It can be reduced to make the integration faster but it will be difficult to reach the required tolerance
///      - VEGAS   MC integration method based on importance sampling - numIters is number of function calls
///                Extra Vegas parameter can be set using  IntegratorMultiDimOptions class
///      - MISER    MC integration method based on stratified sampling
///                See also http://en.wikipedia.org/wiki/Monte_Carlo_integration for VEGAS and MISER description
///      - PLAIN    simple MC integration method, where the max  number of calls can be specified using SetNumIters(numIters)
///
/// Extra integration types are:
///
///   - TOYMC:
///       evaluate posterior by generating toy MC for the nuisance parameters. It is a MC
///       integration, where the function is sampled according to the nuisance. It is convenient to use when all
///       the nuisance are uncorrelated and it is efficient to generate them
///       The toy are generated by default for each  poi values
///       (this method has been proposed and provided by J.P Chou)
///   - 1-TOYMC  : same method as before but in this case the toys are generated only one time and then used for
///       each poi value. It can be convenient when the generation time is much larger than the evaluation time,
///       otherwise it is recommended to re-generate the toy for each poi scanned point of the posterior function
///   - ROOFIT:
///       use roofit default integration methods which will produce a nested integral (not recommended for more
///       than 1 nuisance parameters)

void BayesianCalculator::SetIntegrationType(const char * type) {
   // if type = 0 use default specified via class IntegratorMultiDimOptions::SetDefaultIntegrator
   fIntegrationType = TString(type);
   fIntegrationType.ToUpper();
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the interval. By Default a central interval is computed
/// and the result is a SimpleInterval object.
///
/// Using the method (to be called before SetInterval) SetLeftSideTailFraction the user can choose the type of interval.
/// By default the returned interval is a central interval with the confidence level specified
/// previously in the constructor ( LeftSideTailFraction = 0.5).
///  - For lower limit use SetLeftSideTailFraction = 1
///  - For upper limit use SetLeftSideTailFraction = 0
///  - for shortest intervals use SetLeftSideTailFraction = -1 or call the method SetShortestInterval()
///
/// NOTE: The BayesianCalculator covers only the case with one
/// single parameter of interest
///
/// NOTE: User takes ownership of the returned object

SimpleInterval* BayesianCalculator::GetInterval() const
{

   if (fValidInterval)
      coutW(Eval) << "BayesianCalculator::GetInterval - recomputing interval for the same CL and same model" << std::endl;

   RooRealVar* poi = dynamic_cast<RooRealVar*>( fPOI.first() );
   if (!poi) {
      coutE(Eval) << "BayesianCalculator::GetInterval - no parameter of interest is set " << std::endl;
      return 0;
   }



   // get integrated likelihood (posterior function)
   GetPosteriorFunction();

   //Bool_t silentMode = (RooMsgService::instance().globalKillBelow() >= RooFit::ERROR || RooMsgService::instance().silentMode()) ;
   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);

   if (fLeftSideFraction < 0 ) {
      // compute short intervals
      ComputeShortestInterval();
   }
   else {
      // compute the other intervals

      double lowerCutOff = fLeftSideFraction * fSize;
      double upperCutOff = 1. - (1.- fLeftSideFraction) * fSize;


      if (fNScanBins > 0) {
         ComputeIntervalFromApproxPosterior(lowerCutOff, upperCutOff);
      }

      else {
         // use integration method if there are nuisance parameters
         if (fNuisanceParameters.getSize() > 0) {
            ComputeIntervalFromCdf(lowerCutOff, upperCutOff);
         }
         else {
            // case of no nuisance - just use createCdf from roofit
            ComputeIntervalUsingRooFit(lowerCutOff, upperCutOff);
         }
         // case cdf failed (scan then the posterior)
         if (!fValidInterval) {
            fNScanBins = 100;
            coutW(Eval) << "BayesianCalculator::GetInterval - computing integral from cdf failed - do a scan in "
                        << fNScanBins << " nbins " << std::endl;
            ComputeIntervalFromApproxPosterior(lowerCutOff, upperCutOff);
         }
      }
   }


   // reset the counts and default mode
   if (RooAbsReal::numEvalErrors() > 0)
      coutW(Eval) << "BayesianCalculator::GetInterval : " << RooAbsReal::numEvalErrors() << " errors reported in evaluating log-likelihood function "
                   << std::endl;

   RooAbsReal::clearEvalErrorLog();
   RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);

   if (!fValidInterval) {
      fLower = 1; fUpper = 0;
      coutE(Eval) << "BayesianCalculator::GetInterval - cannot compute a valid interval - return a dummy [1,0] interval"
      <<  std::endl;
   }
   else {
      coutI(Eval) << "BayesianCalculator::GetInterval - found a valid interval : [" << fLower << " , "
                << fUpper << " ]" << std::endl;
   }

   TString interval_name = TString("BayesianInterval_a") + TString(this->GetName());
   SimpleInterval * interval = new SimpleInterval(interval_name,*poi,fLower,fUpper,ConfidenceLevel());
   interval->SetTitle("SimpleInterval from BayesianCalculator");

   return interval;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the value of the parameter for the point in
/// parameter-space that is the most likely.
///  How do we do if there are points that are equi-probable?
/// use approximate posterior
/// t.b.d use real function to find the mode

double BayesianCalculator::GetMode() const {

   ApproximatePosterior();
   TH1 * h = fApproxPosterior->GetHistogram();
   return h->GetBinCenter(h->GetMaximumBin() );
   //return  fApproxPosterior->GetMaximumX();
}

////////////////////////////////////////////////////////////////////////////////
/// internal function compute the interval using RooFit

void BayesianCalculator::ComputeIntervalUsingRooFit(double lowerCutOff, double upperCutOff ) const {

   coutI(Eval) <<  "BayesianCalculator: Compute interval using RooFit:  posteriorPdf + createCdf + RooBrentRootFinder " << std::endl;

   RooRealVar* poi = dynamic_cast<RooRealVar*>( fPOI.first() );
   assert(poi);

   fValidInterval = false;
   if (!fPosteriorPdf) fPosteriorPdf = (RooAbsPdf*) GetPosteriorPdf();
   if (!fPosteriorPdf) return;

   RooAbsReal* cdf = fPosteriorPdf->createCdf(fPOI,RooFit::ScanNoCdf());
   if (!cdf) return;

   RooAbsFunc* cdf_bind = cdf->bindVars(fPOI,&fPOI);
   if (!cdf_bind) return;

   RooBrentRootFinder brf(*cdf_bind);
   brf.setTol(fBrfPrecision); // set the brf precision

   double tmpVal = poi->getVal();  // patch used because findRoot changes the value of poi

   bool ret = true;
   if (lowerCutOff > 0) {
      double y = lowerCutOff;
      ret &= brf.findRoot(fLower,poi->getMin(),poi->getMax(),y);
   }
   else
      fLower = poi->getMin();

   if (upperCutOff < 1.0) {
      double y=upperCutOff;
      ret &= brf.findRoot(fUpper,poi->getMin(),poi->getMax(),y);
   }
   else
      fUpper = poi->getMax();
   if (!ret)  coutE(Eval) << "BayesianCalculator::GetInterval "
                           << "Error returned from Root finder, estimated interval is not fully correct"
                           << std::endl;
   else
      fValidInterval = true;


   poi->setVal(tmpVal); // patch: restore the original value of poi

   delete cdf_bind;
   delete cdf;
}

////////////////////////////////////////////////////////////////////////////////
/// internal function compute the interval using Cdf integration

void BayesianCalculator::ComputeIntervalFromCdf(double lowerCutOff, double upperCutOff ) const {

   fValidInterval = false;

   coutI(InputArguments) <<  "BayesianCalculator:GetInterval Compute the interval from the posterior cdf " << std::endl;

   RooRealVar* poi = dynamic_cast<RooRealVar*>( fPOI.first() );
   assert(poi);
   if (GetPosteriorFunction() == 0) {
      coutE(InputArguments) <<  "BayesianCalculator::GetInterval() cannot make posterior Function " << std::endl;
      return;
   }

   // need to remove the constant parameters
   RooArgList bindParams;
   bindParams.add(fPOI);
   bindParams.add(fNuisanceParameters);

   // this code could be put inside the PosteriorCdfFunction

   //bindParams.Print("V");

   PosteriorCdfFunction cdf(*fLogLike, bindParams, fPriorPdf, fIntegrationType, fNLLMin );
   if( cdf.HasError() ) {
      coutE(Eval) <<  "BayesianCalculator: Numerical error computing CDF integral - try a different method " << std::endl;
      return;
   }

   //find the roots

   ROOT::Math::RootFinder rf(kRootFinderType);

   ccoutD(Eval) << "BayesianCalculator::GetInterval - finding roots of posterior using RF " << rf.Name()
                << " with precision = " << fBrfPrecision;

   if (lowerCutOff > 0) {
      cdf.SetOffset(lowerCutOff);
      ccoutD(NumIntegration) << "Integrating posterior to get cdf and search lower limit at p =" << lowerCutOff << std::endl;
      bool ok = rf.Solve(cdf, poi->getMin(),poi->getMax() , 200,fBrfPrecision, fBrfPrecision);
      if( cdf.HasError() )
         coutW(Eval) <<  "BayesianCalculator: Numerical error integrating the  CDF   " << std::endl;
      if (!ok) {
         coutE(NumIntegration) << "BayesianCalculator::GetInterval - Error from root finder when searching lower limit !" << std::endl;
         return;
      }
      fLower = rf.Root();
   }
   else {
      fLower = poi->getMin();
   }
   if (upperCutOff < 1.0) {
      cdf.SetOffset(upperCutOff);
      ccoutD(NumIntegration) << "Integrating posterior to get cdf and search upper interval limit at p =" << upperCutOff << std::endl;
      bool ok = rf.Solve(cdf, fLower,poi->getMax() , 200, fBrfPrecision, fBrfPrecision);
      if( cdf.HasError() )
         coutW(Eval) <<  "BayesianCalculator: Numerical error integrating the  CDF   " << std::endl;
      if (!ok)  {
         coutE(NumIntegration) << "BayesianCalculator::GetInterval - Error from root finder when searching upper limit !" << std::endl;
         return;
      }
      fUpper = rf.Root();
   }
   else {
      fUpper = poi->getMax();
   }

   fValidInterval = true;
}

////////////////////////////////////////////////////////////////////////////////
/// approximate posterior in nbins using a TF1
/// scan the poi values and evaluate the posterior at each point
/// and save the result in a cloned TF1
/// For each point the posterior is evaluated by integrating the nuisance
/// parameters

void BayesianCalculator::ApproximatePosterior() const {

   if (fApproxPosterior) {
      // if number of bins of existing function is >= requested one - no need to redo the scan
      if (fApproxPosterior->GetNpx() >= fNScanBins) return;
      // otherwise redo the scan
      delete fApproxPosterior;
      fApproxPosterior = 0;
   }


   RooAbsReal * posterior = GetPosteriorFunction();
   if (!posterior) return;


   TF1 * tmp = posterior->asTF(fPOI);
   assert(tmp != 0);
   // binned the function in nbins and evaluate at those points
   if (fNScanBins > 0)  tmp->SetNpx(fNScanBins);  // if not use default of TF1 (which is 100)

   coutI(Eval) << "BayesianCalculator - scan posterior function in nbins = " << tmp->GetNpx() << std::endl;

   fApproxPosterior = (TF1*) tmp->Clone();
   // save this function for future reuse
   // I can delete now original posterior and use this approximated copy
   delete tmp;
   TString name = posterior->GetName() + TString("_approx");
   TString title = posterior->GetTitle() + TString("_approx");
   RooAbsReal * posterior2 = new RooTFnBinding(name,title,fApproxPosterior,fPOI);
   if (posterior == fIntegratedLikelihood) {
      delete fIntegratedLikelihood;
      fIntegratedLikelihood = posterior2;
   }
   else if (posterior == fLikelihood) {
      delete fLikelihood;
      fLikelihood = posterior2;
   }
   else {
      assert(1); // should never happen this case
   }
}

////////////////////////////////////////////////////////////////////////////////
/// compute the interval using the approximate posterior function

void BayesianCalculator::ComputeIntervalFromApproxPosterior(double lowerCutOff, double upperCutOff ) const {

   ccoutD(Eval) <<  "BayesianCalculator: Compute interval from the approximate posterior " << std::endl;

   ApproximatePosterior();
   if (!fApproxPosterior) return;

   double prob[2];
   double limits[2] = {0,0};
   prob[0] = lowerCutOff;
   prob[1] = upperCutOff;
   fApproxPosterior->GetQuantiles(2,limits,prob);
   fLower = limits[0];
   fUpper = limits[1];
   fValidInterval = true;
}

////////////////////////////////////////////////////////////////////////////////
/// compute the shortest interval from the histogram representing the posterior


void BayesianCalculator::ComputeShortestInterval( ) const {
   coutI(Eval) << "BayesianCalculator - computing shortest interval with CL = " << 1.-fSize << std::endl;

   // compute via the approx posterior function
   ApproximatePosterior();
   if (!fApproxPosterior) return;

   TH1D * h1 = dynamic_cast<TH1D*>(fApproxPosterior->GetHistogram() );
   assert(h1 != 0);
   h1->SetName(fApproxPosterior->GetName());
   // get bins and sort them
   double * bins = h1->GetArray();
   // exclude under/overflow bins
   int n = h1->GetSize()-2;
   std::vector<int> index(n);
   //  exclude bins[0] (the underflow bin content)
   TMath::Sort(n, bins+1, &index[0]);
   // find cut off as test size
   double sum = 0;
   double actualCL = 0;
   double upper =  h1->GetXaxis()->GetXmin();
   double lower =  h1->GetXaxis()->GetXmax();
   double norm = h1->GetSumOfWeights();

   for (int i = 0; i < n; ++i)  {
      int idx = index[i];
      double p = bins[ idx] / norm;
      sum += p;
      if (sum > 1.-fSize ) {
         actualCL = sum - p;
         break;
      }

      // histogram bin content starts from 1
      if ( h1->GetBinLowEdge(idx+1) < lower)
         lower = h1->GetBinLowEdge(idx+1);
      if ( h1->GetXaxis()->GetBinUpEdge(idx+1) > upper)
         upper = h1->GetXaxis()->GetBinUpEdge(idx+1);
   }

   ccoutD(Eval) << "BayesianCalculator::ComputeShortestInterval - actual interval CL = "
                << actualCL << " difference from requested is " << (actualCL-(1.-fSize))/fSize*100. << "%  "
                << " limits are [ " << lower << " , " << " upper ] " << std::endl;


   if (lower < upper) {
      fLower = lower;
      fUpper = upper;



      if ( std::abs(actualCL-(1.-fSize)) > 0.1*(1.-fSize) )
         coutW(Eval) << "BayesianCalculator::ComputeShortestInterval - actual interval CL = "
                     << actualCL << " differs more than 10% from desired CL value - must increase nbins "
                     << n << " to an higher value " << std::endl;
   }
   else
      coutE(Eval) << "BayesianCalculator::ComputeShortestInterval " << n << " bins are not sufficient " << std::endl;

   fValidInterval = true;

}




} // end namespace RooStats
