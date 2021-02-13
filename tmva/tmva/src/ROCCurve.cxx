// @(#)root/tmva $Id$
// Author: Omar Zapata, Lorenzo Moneta, Sergei Gleyzer, Simon Pfreundschuh and Kim Albertsson

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : ROCCurve                                                              *
 *                                                                                *
 * Description:                                                                   *
 *      This is class to compute ROC Integral (AUC)                               *
 *                                                                                *
 * Authors :                                                                      *
 *      Omar Zapata     <Omar.Zapata@cern.ch>    - UdeA/ITM Colombia              *
 *      Lorenzo Moneta  <Lorenzo.Moneta@cern.ch> - CERN, Switzerland              *
 *      Sergei Gleyzer  <Sergei.Gleyzer@cern.ch> - U of Florida & CERN            *
 *      Kim Albertsson  <kim.albertsson@cern.ch> - LTU & CERN                     *
 *                                                                                *
 * Copyright (c) 2015:                                                            *
 *      CERN, Switzerland                                                         *
 *      UdeA/ITM, Colombia                                                        *
 *      U. of Florida, USA                                                        *
 **********************************************************************************/

/*! \class TMVA::ROCCurve
\ingroup TMVA

*/
#include "TMVA/Tools.h"
#include "TMVA/TSpline1.h"
#include "TMVA/ROCCurve.h"
#include "TMVA/Config.h"
#include "TMVA/Version.h"
#include "TMVA/MsgLogger.h"
#include "TGraph.h"

#include <algorithm>
#include <vector>
#include <cassert>

using namespace std;

auto tupleSort = [](std::tuple<Float_t, Float_t, Bool_t> _a, std::tuple<Float_t, Float_t, Bool_t> _b) {
   return std::get<0>(_a) < std::get<0>(_b);
};

//_______________________________________________________________________
TMVA::ROCCurve::ROCCurve(const std::vector<std::tuple<Float_t, Float_t, Bool_t>> &mvas)
   : fLogger(new TMVA::MsgLogger("ROCCurve")), fGraph(NULL), fMva(mvas)
{
}

////////////////////////////////////////////////////////////////////////////////
///

TMVA::ROCCurve::ROCCurve(const std::vector<Float_t> &mvaValues, const std::vector<Bool_t> &mvaTargets,
                         const std::vector<Float_t> &mvaWeights)
   : fLogger(new TMVA::MsgLogger("ROCCurve")), fGraph(NULL)
{
   assert(mvaValues.size() == mvaTargets.size());
   assert(mvaValues.size() == mvaWeights.size());

   for (UInt_t i = 0; i < mvaValues.size(); i++) {
      fMva.emplace_back(mvaValues[i], mvaWeights[i], mvaTargets[i]);
   }

   std::sort(fMva.begin(), fMva.end(), tupleSort);
}

////////////////////////////////////////////////////////////////////////////////
///

TMVA::ROCCurve::ROCCurve(const std::vector<Float_t> &mvaValues, const std::vector<Bool_t> &mvaTargets)
   : fLogger(new TMVA::MsgLogger("ROCCurve")), fGraph(NULL)
{
   assert(mvaValues.size() == mvaTargets.size());

   for (UInt_t i = 0; i < mvaValues.size(); i++) {
      fMva.emplace_back(mvaValues[i], 1, mvaTargets[i]);
   }

   std::sort(fMva.begin(), fMva.end(), tupleSort);
}

////////////////////////////////////////////////////////////////////////////////
///

TMVA::ROCCurve::ROCCurve(const std::vector<Float_t> &mvaSignal, const std::vector<Float_t> &mvaBackground)
   : fLogger(new TMVA::MsgLogger("ROCCurve")), fGraph(NULL)
{
   for (UInt_t i = 0; i < mvaSignal.size(); i++) {
      fMva.emplace_back(mvaSignal[i], 1, kTRUE);
   }

   for (UInt_t i = 0; i < mvaBackground.size(); i++) {
      fMva.emplace_back(mvaBackground[i], 1, kFALSE);
   }

   std::sort(fMva.begin(), fMva.end(), tupleSort);
}

////////////////////////////////////////////////////////////////////////////////
///

TMVA::ROCCurve::ROCCurve(const std::vector<Float_t> &mvaSignal, const std::vector<Float_t> &mvaBackground,
                         const std::vector<Float_t> &mvaSignalWeights, const std::vector<Float_t> &mvaBackgroundWeights)
   : fLogger(new TMVA::MsgLogger("ROCCurve")), fGraph(NULL)
{
   assert(mvaSignal.size() == mvaSignalWeights.size());
   assert(mvaBackground.size() == mvaBackgroundWeights.size());

   for (UInt_t i = 0; i < mvaSignal.size(); i++) {
      fMva.emplace_back(mvaSignal[i], mvaSignalWeights[i], kTRUE);
   }

   for (UInt_t i = 0; i < mvaBackground.size(); i++) {
      fMva.emplace_back(mvaBackground[i], mvaBackgroundWeights[i], kFALSE);
   }

   std::sort(fMva.begin(), fMva.end(), tupleSort);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::ROCCurve::~ROCCurve() {
   delete fLogger;
   if(fGraph) delete fGraph;
}

TMVA::MsgLogger &TMVA::ROCCurve::Log() const
{
   if (!fLogger)
      fLogger = new TMVA::MsgLogger("ROCCurve");
   return *fLogger;
}

////////////////////////////////////////////////////////////////////////////////
///

std::vector<Double_t> TMVA::ROCCurve::ComputeSpecificity(const UInt_t num_points)
{
   if (num_points <= 2) {
      return {0.0, 1.0};
   }

   std::vector<Double_t> specificity_vector;
   std::vector<Double_t> true_negatives;
   specificity_vector.reserve(fMva.size());
   true_negatives.reserve(fMva.size());

   Double_t true_negatives_sum = 0.0;
   for (auto &ev : fMva) {
      // auto value = std::get<0>(ev);
      auto weight = std::get<1>(ev);
      auto isSignal = std::get<2>(ev);

      true_negatives_sum += weight * (!isSignal);
      true_negatives.push_back(true_negatives_sum);
   }

   specificity_vector.push_back(0.0);
   Double_t total_background = true_negatives_sum;
   for (auto &tn : true_negatives) {
      Double_t specificity =
         (total_background <= std::numeric_limits<Double_t>::min()) ? (0.0) : (tn / total_background);
      specificity_vector.push_back(specificity);
   }
   specificity_vector.push_back(1.0);

   return specificity_vector;
}

////////////////////////////////////////////////////////////////////////////////
///

std::vector<Double_t> TMVA::ROCCurve::ComputeSensitivity(const UInt_t num_points)
{
   if (num_points <= 2) {
      return {1.0, 0.0};
   }

   std::vector<Double_t> sensitivity_vector;
   std::vector<Double_t> true_positives;
   sensitivity_vector.reserve(fMva.size());
   true_positives.reserve(fMva.size());

   Double_t true_positives_sum = 0.0;
   for (auto it = fMva.rbegin(); it != fMva.rend(); ++it) {
      // auto value = std::get<0>(*it);
      auto weight = std::get<1>(*it);
      auto isSignal = std::get<2>(*it);

      true_positives_sum += weight * (isSignal);
      true_positives.push_back(true_positives_sum);
   }
   std::reverse(true_positives.begin(), true_positives.end());

   sensitivity_vector.push_back(1.0);
   Double_t total_signal = true_positives_sum;
   for (auto &tp : true_positives) {
      Double_t sensitivity = (total_signal <= std::numeric_limits<Double_t>::min()) ? (0.0) : (tp / total_signal);
      sensitivity_vector.push_back(sensitivity);
   }
   sensitivity_vector.push_back(0.0);

   return sensitivity_vector;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the signal efficiency (sensitivity) for a given background
/// efficiency (sensitivity).
///
/// @param effB         Background efficiency for which to calculate signal
///                     efficiency.
/// @param num_points   Number of points used for the underlying histogram.
///                     The number of bins will be num_points - 1.
///

Double_t TMVA::ROCCurve::GetEffSForEffB(Double_t effB, const UInt_t num_points)
{
   assert(0.0 <= effB && effB <= 1.0);

   auto effS_vec = ComputeSensitivity(num_points);
   auto effB_vec = ComputeSpecificity(num_points);

   // Specificity is actually rejB, so we need to transform it.
   auto complement = [](Double_t x) { return 1 - x; };
   std::transform(effB_vec.begin(), effB_vec.end(), effB_vec.begin(), complement);

   // Since TSpline1 uses binary search (and assumes ascending sorting) we must ensure this.
   std::reverse(effS_vec.begin(), effS_vec.end());
   std::reverse(effB_vec.begin(), effB_vec.end());

   TGraph *graph = new TGraph(effS_vec.size(), &effB_vec[0], &effS_vec[0]);

   // TSpline1 does linear interpolation of ROC curve
   TSpline1 rocSpline = TSpline1("", graph);
   return rocSpline.Eval(effB);
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates the ROC integral (AUC)
///
/// @param num_points Granularity of the resulting curve used for integration.
///                     The curve will be subdivided into num_points - 1 regions
///                     where the performance of the classifier is sampled.
///                     Larger number means more accurate, but more costly,
///                     evaluation.

Double_t TMVA::ROCCurve::GetROCIntegral(const UInt_t num_points)
{
   auto sensitivity = ComputeSensitivity(num_points);
   auto specificity = ComputeSpecificity(num_points);

   Double_t integral = 0.0;
   for (UInt_t i = 0; i < sensitivity.size() - 1; i++) {
      // FNR, false negatigve rate = 1 - Sensitivity
      Double_t currFnr = 1 - sensitivity[i];
      Double_t nextFnr = 1 - sensitivity[i + 1];
      // Trapezodial integration
      integral += 0.5 * (nextFnr - currFnr) * (specificity[i] + specificity[i + 1]);
   }

   return integral;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a new TGraph containing the ROC curve. Sensitivity is on the x-axis,
/// specificity on the y-axis.
///
/// @param num_points Granularity of the resulting curve. The curve will be subdivided
///                     into num_points - 1 regions where the performance of the
///                     classifier is sampled. Larger number means more accurate,
///                     but more costly, evaluation.

TGraph *TMVA::ROCCurve::GetROCCurve(const UInt_t num_points)
{
   if (fGraph != nullptr) {
      delete fGraph;
   }

   auto sensitivity = ComputeSensitivity(num_points);
   auto specificity = ComputeSpecificity(num_points);

   fGraph = new TGraph(sensitivity.size(), &sensitivity[0], &specificity[0]);

   return fGraph;
}
