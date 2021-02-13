// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FunctionMinimum
#define ROOT_Minuit2_FunctionMinimum

#include "Minuit2/BasicFunctionMinimum.h"

#include <vector>
#include <memory>

#ifdef G__DICTIONARY
typedef ROOT::Minuit2::MinimumState MinimumState;
#endif

namespace ROOT {

namespace Minuit2 {

//______________________________________________________________________________________________
/**
    class holding the full result of the minimization;
    both internal and external (MnUserParameterState) representation available
    for the parameters at the Minimum
 */

class FunctionMinimum {

public:
   using MnAboveMaxEdm = BasicFunctionMinimum::MnAboveMaxEdm;
   using MnReachedCallLimit = BasicFunctionMinimum::MnReachedCallLimit;

public:
   /// constructor from only MinimumSeed. Minimum is only from seed result not full minimization
   FunctionMinimum(const MinimumSeed &seed, double up) : fData(std::make_shared<BasicFunctionMinimum>(seed, up)) {}

   /// constructor at the end of a successfull minimization from seed and vector of states
   FunctionMinimum(const MinimumSeed &seed, const std::vector<MinimumState> &states, double up)
      : fData(std::make_shared<BasicFunctionMinimum>(seed, states, up))
   {
   }

   /// constructor at the end of a failed minimization due to exceeding function call limit
   FunctionMinimum(const MinimumSeed &seed, const std::vector<MinimumState> &states, double up, MnReachedCallLimit)
      : fData(std::make_shared<BasicFunctionMinimum>(seed, states, up, MnReachedCallLimit{}))
   {
   }

   /// constructor at the end of a failed minimization due to edm above maximum value
   FunctionMinimum(const MinimumSeed &seed, const std::vector<MinimumState> &states, double up, MnAboveMaxEdm)
      : fData(std::make_shared<BasicFunctionMinimum>(seed, states, up, MnAboveMaxEdm{}))
   {
   }

   // add new state
   void Add(const MinimumState &state) { fData->Add(state); }

   // add new state
   void Add(const MinimumState &state, MnAboveMaxEdm) { fData->Add(state, MnAboveMaxEdm{}); }

   const MinimumSeed &Seed() const { return fData->Seed(); }
   const std::vector<ROOT::Minuit2::MinimumState> &States() const { return fData->States(); }

   // user representation of state at Minimum
   const MnUserParameterState &UserState() const { return fData->UserState(); }
   const MnUserParameters &UserParameters() const { return fData->UserParameters(); }
   const MnUserCovariance &UserCovariance() const { return fData->UserCovariance(); }

   // forward interface of last state
   const MinimumState &State() const { return fData->State(); }
   const MinimumParameters &Parameters() const { return fData->Parameters(); }
   const MinimumError &Error() const { return fData->Error(); }
   const FunctionGradient &Grad() const { return fData->Grad(); }
   double Fval() const { return fData->Fval(); }
   double Edm() const { return fData->Edm(); }
   int NFcn() const { return fData->NFcn(); }

   double Up() const { return fData->Up(); }
   bool IsValid() const { return fData->IsValid(); }
   bool HasValidParameters() const { return fData->HasValidParameters(); }
   bool HasValidCovariance() const { return fData->HasValidCovariance(); }
   bool HasAccurateCovar() const { return fData->HasAccurateCovar(); }
   bool HasPosDefCovar() const { return fData->HasPosDefCovar(); }
   bool HasMadePosDefCovar() const { return fData->HasMadePosDefCovar(); }
   bool HesseFailed() const { return fData->HesseFailed(); }
   bool HasCovariance() const { return fData->HasCovariance(); }
   bool IsAboveMaxEdm() const { return fData->IsAboveMaxEdm(); }
   bool HasReachedCallLimit() const { return fData->HasReachedCallLimit(); }

   void SetErrorDef(double up) { return fData->SetErrorDef(up); }

private:
   std::shared_ptr<BasicFunctionMinimum> fData;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_FunctionMinimum
