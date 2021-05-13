// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FunctionMinimum
#define ROOT_Minuit2_FunctionMinimum

#include "Minuit2/MinimumSeed.h"
#include "Minuit2/MinimumState.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnUserTransformation.h"

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
   enum Status {
      MnValid,
      MnReachedCallLimit,
      MnAboveMaxEdm,
   };

public:
   /// Constructor from only MinimumSeed. Minimum is only from seed result not the full minimization
   FunctionMinimum(const MinimumSeed &seed, double up)
      : FunctionMinimum(seed,
                        std::vector<MinimumState>(1, MinimumState(seed.Parameters(), seed.Error(), seed.Gradient(),
                                                                  seed.Parameters().Fval(), seed.NFcn())),
                        up)
   {
   }

   /// Constructor at the end of a minimization from seed and vector of states
   FunctionMinimum(const MinimumSeed &seed, const std::vector<MinimumState> &states, double up, Status status = MnValid)
      : fPtr{new Data{seed, states, up, status == MnAboveMaxEdm, status == MnReachedCallLimit, {}}}
   {
   }

   /// add latest minimization state (for example add Hesse result after Migrad)
   void Add(const MinimumState &state, Status status = MnValid)
   {
      fPtr->fStates.push_back(state);
      // LM : update also the user state
      fPtr->fUserState = MnUserParameterState(State(), Up(), Seed().Trafo());
      fPtr->fAboveMaxEdm = status == MnAboveMaxEdm;
      fPtr->fReachedCallLimit = status == MnReachedCallLimit;
   }

   const MinimumSeed &Seed() const { return fPtr->fSeed; }
   const std::vector<MinimumState> &States() const { return fPtr->fStates; }

   // user representation of state at Minimum
   const MnUserParameterState &UserState() const
   {
      if (!fPtr->fUserState.IsValid())
         fPtr->fUserState = MnUserParameterState(State(), Up(), Seed().Trafo());
      return fPtr->fUserState;
   }
   const MnUserParameters &UserParameters() const { return UserState().Parameters(); }
   const MnUserCovariance &UserCovariance() const { return UserState().Covariance(); }

   // forward interface of last state
   const MinimumState &State() const { return States().back(); }
   const MinimumParameters &Parameters() const { return States().back().Parameters(); }
   const MinimumError &Error() const { return States().back().Error(); }
   const FunctionGradient &Grad() const { return States().back().Gradient(); }
   double Fval() const { return States().back().Fval(); }
   double Edm() const { return States().back().Edm(); }
   int NFcn() const { return States().back().NFcn(); }

   double Up() const { return fPtr->fErrorDef; }
   bool IsValid() const { return State().IsValid() && !IsAboveMaxEdm() && !HasReachedCallLimit(); }
   bool HasValidParameters() const { return State().Parameters().IsValid(); }
   bool HasValidCovariance() const { return State().Error().IsValid(); }
   bool HasAccurateCovar() const { return State().Error().IsAccurate(); }
   bool HasPosDefCovar() const { return State().Error().IsPosDef(); }
   bool HasMadePosDefCovar() const { return State().Error().IsMadePosDef(); }
   bool HesseFailed() const { return State().Error().HesseFailed(); }
   bool HasCovariance() const { return State().Error().IsAvailable(); }
   bool IsAboveMaxEdm() const { return fPtr->fAboveMaxEdm; }
   bool HasReachedCallLimit() const { return fPtr->fReachedCallLimit; }

   void SetErrorDef(double up)
   {
      fPtr->fErrorDef = up;
      // update user state for new valeu of up (scaling of errors)
      fPtr->fUserState = MnUserParameterState(State(), up, Seed().Trafo());
   }

private:
   struct Data {
      MinimumSeed fSeed;
      std::vector<MinimumState> fStates;
      double fErrorDef;
      bool fAboveMaxEdm;
      bool fReachedCallLimit;
      mutable MnUserParameterState fUserState;
   };

   std::shared_ptr<Data> fPtr;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_FunctionMinimum
