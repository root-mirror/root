// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef RooStats_LikelihoodInterval
#define RooStats_LikelihoodInterval

#include "RooStats/ConfInterval.h"

#include "RooArgSet.h"

#include "RooAbsReal.h"

#include "Math/IFunctionfwd.h"

#include <map>
#include <memory>
#include <string>

namespace ROOT {
   namespace Math {
      class Minimizer;
   }
}

namespace RooStats {

   class LikelihoodInterval : public ConfInterval {

   public:

      /// default constructor
      explicit LikelihoodInterval(const char* name = 0);

      //// construct the interval from a Profile Likelihood object, parameter of interest and optionally a snapshot of
      //// POI with their best fit values
      LikelihoodInterval(const char* name, RooAbsReal*, const RooArgSet*,  RooArgSet * = 0);

      /// destructor
      virtual ~LikelihoodInterval();

      /// check if given point is in the interval
      virtual Bool_t IsInInterval(const RooArgSet&) const;

      /// set the confidence level for the interval (e.g 0.682 for a 1-sigma interval)
      virtual void SetConfidenceLevel(Double_t cl) {fConfidenceLevel = cl; ResetLimits(); }

      /// return confidence level
      virtual Double_t ConfidenceLevel() const {return fConfidenceLevel;}

      /// return a cloned list of parameters of interest.  User manages the return object
      virtual  RooArgSet* GetParameters() const;

      /// check if parameters are correct (i.e. they are the POI of this interval)
      Bool_t CheckParameters(const RooArgSet&) const ;


      /// return the lower bound of the interval on a given parameter
      Double_t LowerLimit(const RooRealVar& param) { bool ok; return LowerLimit(param,ok); }
      Double_t LowerLimit(const RooRealVar& param, bool & status) ;

      /// return the upper bound of the interval on a given parameter
      Double_t UpperLimit(const RooRealVar& param) { bool ok; return UpperLimit(param,ok); }
      Double_t UpperLimit(const RooRealVar& param, bool & status) ;

      /// find both lower and upper interval boundaries for a given parameter
      /// return false if the bounds have not been found
      Bool_t FindLimits(const RooRealVar & param, double & lower, double &upper);

      /// return the 2D-contour points for the given subset of parameters
      /// by default make the contour using 30 points. The User has to preallocate the x and y array which will return
      /// the set of x and y points defining the contour.
      /// The return value of the function specify the number of contour point found.
      /// In case of error a zero is returned
      Int_t GetContourPoints(const RooRealVar & paramX, const RooRealVar & paramY, Double_t * x, Double_t *y, Int_t npoints = 30);

      /// return the profile log-likelihood ratio function
      RooAbsReal* GetLikelihoodRatio() {return fLikelihoodRatio;}

      /// return a pointer to a snapshot with best fit parameter of interest
      const RooArgSet * GetBestFitParameters() const { return fBestFitParams; }

   protected:

      /// reset the cached limit values
      void ResetLimits();

      /// internal function to create the minimizer for finding the contours
      bool CreateMinimizer();

   private:

      RooArgSet   fParameters; /// parameters of interest for this interval
      RooArgSet * fBestFitParams; /// snapshot of the model parameters with best fit value (managed internally)
      RooAbsReal* fLikelihoodRatio; /// likelihood ratio function used to make contours (managed internally)
      Double_t fConfidenceLevel; /// Requested confidence level (eg. 0.95 for 95% CL)
      std::map<std::string, double> fLowerLimits; /// map with cached lower bound values
      std::map<std::string, double> fUpperLimits; /// map with cached upper bound values
      std::shared_ptr<ROOT::Math::Minimizer > fMinimizer; //! transient pointer to minimizer class used to find limits and contour
      std::shared_ptr<RooFunctor>           fFunctor;   //! transient pointer to functor class used by the minimizer
      std::shared_ptr<ROOT::Math::IMultiGenFunction> fMinFunc; //! transient pointer to the minimization function

      ClassDef(LikelihoodInterval,1)  /// Concrete implementation of a ConfInterval based on a likelihood ratio

   };
}

#endif
