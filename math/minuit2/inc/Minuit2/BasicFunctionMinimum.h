// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_BasicFunctionMinimum
#define ROOT_Minuit2_BasicFunctionMinimum

#include "Minuit2/MinimumSeed.h"
#include "Minuit2/MinimumState.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnUserTransformation.h"

#include "Minuit2/StackAllocator.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

// extern StackAllocator gStackAllocator;

/** result of the minimization;
    both internal and external (MnUserParameterState) representation available
    For the parameters at the Minimum
 */

class BasicFunctionMinimum {

public:
   class MnReachedCallLimit {
   };
   class MnAboveMaxEdm {
   };

public:
   /// Constructor from only MinimumSeed. Minimum is only from seed result not the full minimization
   BasicFunctionMinimum(const MinimumSeed &seed, double up)
      : fSeed(seed), fStates(std::vector<MinimumState>(1, MinimumState(seed.Parameters(), seed.Error(), seed.Gradient(),
                                                                       seed.Parameters().Fval(), seed.NFcn()))),
        fErrorDef(up), fAboveMaxEdm(false), fReachedCallLimit(false), fUserState(MnUserParameterState())
   {
   }

   /// Constructor at the end of a successfull minimization from seed and vector of states
   BasicFunctionMinimum(const MinimumSeed &seed, const std::vector<MinimumState> &states, double up)
      : fSeed(seed), fStates(states), fErrorDef(up), fAboveMaxEdm(false), fReachedCallLimit(false),
        fUserState(MnUserParameterState())
   {
   }

   /// Constructor at the end of a failed minimization due to exceeding function call limit
   BasicFunctionMinimum(const MinimumSeed &seed, const std::vector<MinimumState> &states, double up, MnReachedCallLimit)
      : fSeed(seed), fStates(states), fErrorDef(up), fAboveMaxEdm(false), fReachedCallLimit(true),
        fUserState(MnUserParameterState())
   {
   }

   /// Constructor at the end of a failed minimization due to edm above maximum value
   BasicFunctionMinimum(const MinimumSeed &seed, const std::vector<MinimumState> &states, double up, MnAboveMaxEdm)
      : fSeed(seed), fStates(states), fErrorDef(up), fAboveMaxEdm(true), fReachedCallLimit(false),
        fUserState(MnUserParameterState())
   {
   }

   /// Copy constructor
   BasicFunctionMinimum(const BasicFunctionMinimum &min)
      : fSeed(min.fSeed), fStates(min.fStates), fErrorDef(min.fErrorDef), fAboveMaxEdm(min.fAboveMaxEdm),
        fReachedCallLimit(min.fReachedCallLimit), fUserState(min.fUserState)
   {
   }

   BasicFunctionMinimum &operator=(const BasicFunctionMinimum &min)
   {
      fSeed = min.fSeed;
      fStates = min.fStates;
      fErrorDef = min.fErrorDef;
      fAboveMaxEdm = min.fAboveMaxEdm;
      fReachedCallLimit = min.fReachedCallLimit;
      fUserState = min.fUserState;
      return *this;
   }

   ~BasicFunctionMinimum() {}

   /// add latest minimization state (for example add Hesse result after Migrad)
   void Add(const MinimumState &state)
   {
      fStates.push_back(state);
      // LM : update also the user state
      fUserState = MnUserParameterState(State(), Up(), Seed().Trafo());
      // reset maxedm flag. If new state has edm over max other method must be used
      fAboveMaxEdm = false;
   }

   /// Add a new state and flag that edm is above maximum
   void Add(const MinimumState &state, MnAboveMaxEdm)
   {
      Add(state);
      fAboveMaxEdm = true;
   }

   const MinimumSeed &Seed() const { return fSeed; }
   const std::vector<MinimumState> &States() const { return fStates; }

   // user representation of state at Minimum
   const MnUserParameterState &UserState() const
   {
      if (!fUserState.IsValid())
         fUserState = MnUserParameterState(State(), Up(), Seed().Trafo());
      return fUserState;
   }
   const MnUserParameters &UserParameters() const
   {
      if (!fUserState.IsValid())
         fUserState = MnUserParameterState(State(), Up(), Seed().Trafo());
      return fUserState.Parameters();
   }
   const MnUserCovariance &UserCovariance() const
   {
      if (!fUserState.IsValid())
         fUserState = MnUserParameterState(State(), Up(), Seed().Trafo());
      return fUserState.Covariance();
   }

   void *operator new(size_t nbytes) { return StackAllocatorHolder::Get().Allocate(nbytes); }

   void operator delete(void *p, size_t /*nbytes */) { StackAllocatorHolder::Get().Deallocate(p); }

   // forward interface of last state
   const MinimumState &State() const { return fStates.back(); }
   const MinimumParameters &Parameters() const { return fStates.back().Parameters(); }
   const MinimumError &Error() const { return fStates.back().Error(); }
   const FunctionGradient &Grad() const { return fStates.back().Gradient(); }
   double Fval() const { return fStates.back().Fval(); }
   double Edm() const { return fStates.back().Edm(); }
   int NFcn() const { return fStates.back().NFcn(); }

   double Up() const { return fErrorDef; }
   bool IsValid() const { return State().IsValid() && !IsAboveMaxEdm() && !HasReachedCallLimit(); }
   bool HasValidParameters() const { return State().Parameters().IsValid(); }
   bool HasValidCovariance() const { return State().Error().IsValid(); }
   bool HasAccurateCovar() const { return State().Error().IsAccurate(); }
   bool HasPosDefCovar() const { return State().Error().IsPosDef(); }
   bool HasMadePosDefCovar() const { return State().Error().IsMadePosDef(); }
   bool HesseFailed() const { return State().Error().HesseFailed(); }
   bool HasCovariance() const { return State().Error().IsAvailable(); }
   bool IsAboveMaxEdm() const { return fAboveMaxEdm; }
   bool HasReachedCallLimit() const { return fReachedCallLimit; }

   void SetErrorDef(double up)
   {
      fErrorDef = up;
      // update user state for new valeu of up (scaling of errors)
      fUserState = MnUserParameterState(State(), up, Seed().Trafo());
   }

private:
   MinimumSeed fSeed;
   std::vector<MinimumState> fStates;
   double fErrorDef;
   bool fAboveMaxEdm;
   bool fReachedCallLimit;
   mutable MnUserParameterState fUserState;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_BasicFunctionMinimum
