// @(#)root/roostats:$Id$

/*************************************************************************
 * Project: RooStats                                                     *
 * Package: RooFit/RooStats                                              *
 * Authors:                                                              *
 *   Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke       *
 *************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_HybridResult
#define ROOSTATS_HybridResult

#include "RooStats/HypoTestResult.h"

#include <vector>

namespace RooStats {

   class HybridPlot;

   class HybridResult : public HypoTestResult {

   public:

      /// Default constructor
      explicit HybridResult(const char *name = 0);

      /// Constructor for HybridResult
      HybridResult(const char *name,
         const std::vector<double>& testStat_sb_vals,
                   const std::vector<double>& testStat_b_vals,
         bool sumLargerValues=true);


      /// Destructor of HybridResult
      virtual ~HybridResult();

      void SetDataTestStatistics(double testStat_data_val);

      void Add(HybridResult* other);

      HybridPlot* GetPlot(const char* name,const char* title, int n_bins);

      void PrintMore(const char* options);

      /// Get test statistics values for the sb model
      std::vector<double> GetTestStat_sb(){return fTestStat_sb;}

      /// Get test statistics values for the b model
      std::vector<double> GetTestStat_b(){return fTestStat_b;}

      /// Get test statistics value for data
      double GetTestStat_data(){ return fTestStat_data;}

      // Return p-value for null hypothesis
      Double_t NullPValue() const;

      // Return p-value for alternate hypothesis
      Double_t AlternatePValue() const;

      /// The error on the "confidence level" of the null hypothesis
      Double_t CLbError() const;

      /// The error on the "confidence level" of the alternative hypothesis
      Double_t CLsplusbError() const;

      /// The error on the ratio \f$CL_{s+b}/CL_{b}\f$
      Double_t CLsError() const;

   private:

      std::vector<double> fTestStat_b; // vector of results for B-only toy-MC
      std::vector<double> fTestStat_sb; // vector of results for S+B toy-MC
      double fTestStat_data; // results (test statistics) evaluated for data

      mutable bool fComputationsNulDoneFlag; // flag if the fNullPValue computation have been already done or not (ie need to be refreshed)
      mutable bool fComputationsAltDoneFlag; // flag if the fAlternatePValue computation have been already done or not (ie need to be refreshed)
      bool fSumLargerValues; // p-value for velues of testStat >= testStat_data (or testStat <= testStat_data)

   protected:

      ClassDef(HybridResult,1)  // Class containing the results of the HybridCalculator
   };
}

#endif
