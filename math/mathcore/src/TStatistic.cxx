// @(#)root/base:$Id$
// Author: G. Ganis 2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TStatistic.h"

#include "TROOT.h"
#include "TList.h"

// clang-format off
/**
* \class TStatistic
* \ingroup MathCore
* \brief Statistical variable, defined by its mean and variance (RMS). Named, streamable, storable and mergeable.
*/
// clang-format on

templateClassImp(TStatistic);

////////////////////////////////////////////////////////////////////////////
/// \brief Constructor from a vector of values
/// \param[in] name The name given to the object
/// \param[in] n The total number of entries
/// \param[in] val The vector of values
/// \param[in] w The vector of weights for the values
///
/// Recursively calls the TStatistic::Fill() function to fill the object.
TStatistic::TStatistic(const char *name, Int_t n, const Double_t *val, const Double_t *w)
         : fName(name), fN(0), fW(0.), fW2(0.), fM(0.), fM2(0.), fMin(TMath::Limits<Double_t>::Max()), fMax(TMath::Limits<Double_t>::Min())
{
   if (n > 0) {
      for (Int_t i = 0; i < n; i++) {
         if (w) {
            Fill(val[i], w[i]);
         } else {
            Fill(val[i]);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TStatistic destructor.
TStatistic::~TStatistic()
{
   // Required since we overload TObject::Hash.
   ROOT::CallRecursiveRemoveIfNeeded(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Increment the entries in the object by one value-weight pair.
/// \param[in] val Value to fill the Tstatistic with
/// \param[in] w The weight of the value
///
/// Also updates statistics in the object. For number of entries, sum of weights,
/// sum of squared weights and sum of (value*weight), one extra value is added to the
/// statistic. For the sum of squared (value*weight) pairs, the function uses formula 1.4
/// in Chan-Golub, LeVeque : Algorithms for computing the Sample Variance (1983),
/// genralized by LM for the case of weights:
/// \f[
///   \frac{w_j}{\sum_{i=0}^{j} w_i \cdot \sum_{i=0}^{j-1} w_i} \cdot
///   \left(
///         \sum_{i=0}^{j} w_i \cdot val_i -
///         \sum_{i=0}^{j} \left(w_i \cdot val_i\right)
///   \right)
/// \f]
///
/// The minimum(maximum) is computed by checking that the fill value
/// is either less(greater) than the current minimum(maximum)
void TStatistic::Fill(Double_t val, Double_t w) {


   if (w == 0) return;
   // increase data count
   fN++;

   // update sum of weights
   Double_t tW = fW + w;

   // update sum of (value * weight) pairs
   fM += w * val;

   // update minimum and maximum values
   fMin = (val < fMin) ? val : fMin;
   fMax = (val > fMax) ? val : fMax;

//      Double_t dt = val - fM ;
   if (tW == 0) {
      Warning("Fill","Sum of weights is zero - ignore current data point");
      fN--;
      return;
   }

   if (fW != 0) {  // from the second time
      Double_t rr = ( tW * val - fM);
      fM2 += w * rr * rr / (tW * fW);
   }
   fW = tW;
   fW2 += w*w;
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Print the content of the object
///
/// Prints the statistics held by the object in one line. These include the mean,
/// mean error, RMS, the total number of values, the minimum and the maximum.
void TStatistic::Print(Option_t *) const {
   TROOT::IndentLevel();
   Printf(" OBJ: TStatistic\t %s \t Mean = %.5g +- %.4g \t RMS = %.5g \t Count = %lld \t Min = %.5g \t Max = %.5g",
          fName.Data(), GetMean(), GetMeanErr(), GetRMS(), GetN(), GetMin(), GetMax());
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Merge implementation of TStatistic
/// \param[in] in Other TStatistic objects to be added to the current one
///
/// The function merges the statistics of all objects together to form a new one.
/// Merging quantities is done via simple addition for the following class data members:
/// - number of entries fN
/// - the sum of weights fW
/// - the sum of squared weights fW2
/// - the sum of (value*weight) fM
///
/// The sum of squared (value*weight) pairs fM2 is updated using the same formula as in
/// TStatistic::Fill() function.
///
/// The minimum(maximum) is updated by checking that the minimum(maximum) of
/// the next TStatistic object in the queue is either less(greater) than the current minimum(maximum).
Int_t TStatistic::Merge(TCollection *in) {

   // Let's organise the list of objects to merge excluding the empty ones
   std::vector<TStatistic*> statPtrs;
   if (this->fN != 0LL) statPtrs.push_back(this);
   TStatistic *statPtr;
   for (auto o : *in) {
      if ((statPtr = dynamic_cast<TStatistic *>(o)) && statPtr->fN != 0LL) {
         statPtrs.push_back(statPtr);
      }
   }

   // No object included this has entries
   const auto nStatsPtrs = statPtrs.size();

   // Early return possible in case nothing has been filled
   if (nStatsPtrs == 0) return 0;

   // Merge the statistic quantities into local variables to then
   // update the data members of this object
   auto firstStatPtr = statPtrs[0];
   auto N = firstStatPtr->fN;
   auto M = firstStatPtr->fM;
   auto M2 = firstStatPtr->fM2;
   auto W = firstStatPtr->fW;
   auto W2 = firstStatPtr->fW2;
   auto Min = firstStatPtr->fMin;
   auto Max = firstStatPtr->fMax;
   for (auto i = 1U; i < nStatsPtrs; ++i) {
      auto c = statPtrs[i];
      double temp = (c->fW) / (W)*M - c->fM;
      M2 += c->fM2 + W / (c->fW * (c->fW + W)) * temp * temp;
      M += c->fM;
      W += c->fW;
      W2 += c->fW2;
      N += c->fN;
      Min = (c->fMin < Min) ? c->fMin : Min;
      Max = (c->fMax > Max) ? c->fMax : Max;
   }

   // Now update the data members of this object
   fN = N;
   fW = W;
   fW2 = W2;
   fM = M;
   fM2 = M2;
   fMin = Min;
   fMax = Max;

   return nStatsPtrs;

}
