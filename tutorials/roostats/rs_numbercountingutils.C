/// \file
/// \ingroup tutorial_roostats
/// \notebook -nodraw
/// 'Number Counting Utils' RooStats tutorial
///
/// This tutorial shows an example of the RooStats standalone
/// utilities that calculate the p-value or Z value (eg. significance in
/// 1-sided Gaussian standard deviations) for a number counting experiment.
/// This is a hypothesis test between background only and signal-plus-background.
/// The background estimate has uncertainty derived from an auxiliary or sideband
/// measurement.
///
/// Documentation for these utilities can be found here:
/// http://root.cern/root/html/RooStats__NumberCountingUtils.html
///
///
/// This problem is often called a proto-type problem for high energy physics.
/// In some references it is referred to as the on/off problem.
///
/// The problem is treated in a fully frequentist fashion by
/// interpreting the relative background uncertainty as
/// being due to an auxiliary or sideband observation
/// that is also Poisson distributed with only background.
/// Finally, one considers the test as a ratio of Poisson means
/// where an interval is well known based on the conditioning on the total
/// number of events and the binomial distribution.
/// For more on this, see
///  - http://arxiv.org/abs/0905.3831
///  - http://arxiv.org/abs/physics/physics/0702156
///  - http://arxiv.org/abs/physics/0511028
///
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Kyle Cranmer

#include "RooStats/RooStatsUtils.h"
#include <iostream>

using namespace RooFit;
using namespace RooStats; // the utilities are in the RooStats namespace
using namespace std;

void rs_numbercountingutils()
{

   // From the root prompt, you can see the full list of functions by using tab-completion
   // ~~~{.bash}
   // root [0] RooStats::NumberCountingUtils::  <tab>
   // BinomialExpZ
   // BinomialWithTauExpZ
   // BinomialObsZ
   // BinomialWithTauObsZ
   // BinomialExpP
   // BinomialWithTauExpP
   // BinomialObsP
   // BinomialWithTauObsP
   // ~~~

   // For each of the utilities you can inspect the arguments by tab completion
   // ~~~{.bash}
   // root [1] NumberCountingUtils::BinomialExpZ( <tab>
   // Double_t BinomialExpZ(Double_t sExp, Double_t bExp, Double_t fractionalBUncertainty)
   // ~~~

   // -------------------------------------------------
   // Here we see common usages where the experimenter
   // has a relative background uncertainty, without
   // explicit reference to the auxiliary or sideband
   // measurement

   // -------------------------------------------------------------
   // Expected p-values and significance with background uncertainty
   double sExpected = 50;
   double bExpected = 100;
   double relativeBkgUncert = 0.1;

   double pExp = NumberCountingUtils::BinomialExpP(sExpected, bExpected, relativeBkgUncert);
   double zExp = NumberCountingUtils::BinomialExpZ(sExpected, bExpected, relativeBkgUncert);
   cout << "expected p-value =" << pExp << "  Z value (Gaussian sigma) = " << zExp << endl;

   // -------------------------------------------------
   // Expected p-values and significance with background uncertainty
   double observed = 150;
   double pObs = NumberCountingUtils::BinomialObsP(observed, bExpected, relativeBkgUncert);
   double zObs = NumberCountingUtils::BinomialObsZ(observed, bExpected, relativeBkgUncert);
   cout << "observed p-value =" << pObs << "  Z value (Gaussian sigma) = " << zObs << endl;

   // ---------------------------------------------------------
   // Here we see usages where the experimenter has knowledge
   // about the properties of the auxiliary or sideband
   // measurement.  In particular, the ratio tau of background
   // in the auxiliary measurement to the main measurement.
   // Large values of tau mean small background uncertainty
   // because the sideband is very constraining.

   // Usage:
   // ~~~{.bash}
   // root [0] RooStats::NumberCountingUtils::BinomialWithTauExpP(
   // Double_t BinomialWithTauExpP(Double_t sExp, Double_t bExp, Double_t tau)
   // ~~~

   // --------------------------------------------------------------
   // Expected p-values and significance with background uncertainty
   double tau = 1;

   double pExpWithTau = NumberCountingUtils::BinomialWithTauExpP(sExpected, bExpected, tau);
   double zExpWithTau = NumberCountingUtils::BinomialWithTauExpZ(sExpected, bExpected, tau);
   cout << "expected p-value =" << pExpWithTau << "  Z value (Gaussian sigma) = " << zExpWithTau << endl;

   // ---------------------------------------------------------------
   // Expected p-values and significance with background uncertainty
   double pObsWithTau = NumberCountingUtils::BinomialWithTauObsP(observed, bExpected, tau);
   double zObsWithTau = NumberCountingUtils::BinomialWithTauObsZ(observed, bExpected, tau);
   cout << "observed p-value =" << pObsWithTau << "  Z value (Gaussian sigma) = " << zObsWithTau << endl;
}
