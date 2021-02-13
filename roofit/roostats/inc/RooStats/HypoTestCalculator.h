// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_HypoTestCalculator
#define ROOSTATS_HypoTestCalculator

#include "Rtypes.h"

// class RooAbsPdf;
class RooAbsData;
// class RooArgSet;
class RooWorkspace;


namespace RooStats {

   class HypoTestResult;

   class ModelConfig;


/** \class HypoTestCalculator
    \ingroup Roostats

HypoTestCalculator is an interface class for a tools which produce RooStats
HypoTestResults. The interface currently assumes that any hypothesis test
calculator can be configured by specifying:

  - a model for the null,
  - a model for the alternate,
  - a data set,
  - a set of parameters of which specify the null (including values and const/non-const status), and
  - a set of parameters of which specify the alternate (including values and const/non-const status).

The interface allows one to pass the model, data, and parameters via a workspace
and then specify them with names. The interface will be extended so that one does
not need to use a workspace.

After configuring the calculator, one only needs to ask GetHypoTest, which will
return a HypoTestResult pointer.

The concrete implementations of this interface should deal with the details of
how the nuisance parameters are dealt with (eg. integration vs. profiling) and
which test-statistic is used (perhaps this should be added to the interface).

The motivation for this interface is that we hope to be able to specify the
problem in a common way for several concrete calculators.

*/


   class HypoTestCalculator {

   public:


      virtual ~HypoTestCalculator() {}

      // main interface to get a HypoTestResult, pure virtual
      virtual HypoTestResult* GetHypoTest() const = 0;

      // Set a common model for both the null and alternate, add to the the workspace if not already there
      virtual void SetCommonModel(const ModelConfig& model) {
         SetNullModel(model);
         SetAlternateModel(model);
      }

      // Set the model for the null hypothesis
      virtual void SetNullModel(const ModelConfig& model) = 0;
      // Set the model for the alternate hypothesis
      virtual void SetAlternateModel(const ModelConfig& model) = 0;
      // Set the DataSet
      virtual void SetData(RooAbsData& data) = 0;


   protected:
      ClassDef(HypoTestCalculator,1)  // Interface for tools doing hypothesis tests
   };
}


#endif
