// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_IntervalCalculator
#define ROOSTATS_IntervalCalculator

#include "Rtypes.h"

class RooAbsData;
class RooWorkspace;

namespace RooStats {

   class ConfInterval;

   class ModelConfig;

/** \class IntervalCalculator
    \ingroup Roostats

IntervalCalculator is an interface class for a tools which produce RooStats
ConfIntervals. The interface currently assumes that any interval calculator can
be configured by specifying:

  - a model,
  - a data set,
  - a set of parameters of interest,
  - a set of nuisance parameters (eg. parameters on which the model depends, but are not of interest), and
  - a confidence level or size of the test (eg. rate of Type I error).

The interface allows one to pass the model, data, and parameters via a workspace
and then specify them with names. The interface will be extended so that one does
not need to use a workspace.

After configuring the calculator, one only needs to ask GetInterval, which will
return a ConfInterval pointer.

The concrete implementations of this interface should deal with the details of
how the nuisance parameters are dealt with (eg. integration vs. profiling) and
which test-statistic is used (perhaps this should be added to the interface).

The motivation for this interface is that we hope to be able to specify the
problem in a common way for several concrete calculators.

*/


   class IntervalCalculator {

   public:

      virtual ~IntervalCalculator() {}

      /// Main interface to get a ConfInterval, pure virtual
      virtual ConfInterval* GetInterval() const = 0;

      /// Get the size of the test (eg. rate of Type I error)
      virtual Double_t Size() const = 0;

      /// Get the Confidence level for the test
      virtual Double_t ConfidenceLevel()  const = 0;

      /// Set the DataSet ( add to the the workspace if not already there ?)
      virtual void SetData(RooAbsData&) = 0;

      /// Set the Model
      virtual void SetModel(const ModelConfig & /* model */) = 0;

      /// set the size of the test (rate of Type I error) ( e.g. 0.05 for a 95% Confidence Interval)
      virtual void SetTestSize(Double_t size) = 0;

      /// set the confidence level for the interval (e.g. 0.95 for a 95% Confidence Interval)
      virtual void SetConfidenceLevel(Double_t cl) = 0;

   protected:
      ClassDef(IntervalCalculator,1)   // Interface for tools setting limits (producing confidence intervals)
   };
}


#endif
