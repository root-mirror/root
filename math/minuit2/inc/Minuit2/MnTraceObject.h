// @(#)root/minuit2:$Id$
// Author:  L. Moneta 2012

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2012 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnTraceObject
#define ROOT_Minuit2_MnTraceObject

namespace ROOT {

namespace Minuit2 {

class MinimumState;
class MnUserParameterState;

class MnTraceObject {

public:
   MnTraceObject(int parNumber = -1) : fUserState(0), fParNumber(parNumber) {}

   virtual ~MnTraceObject() {}

   virtual void Init(const MnUserParameterState &state) { fUserState = &state; }

   virtual void operator()(int i, const MinimumState &state);

   const MnUserParameterState &UserState() const { return *fUserState; }

   void SetParNumber(int number) { fParNumber = number; }

   int ParNumber() const { return fParNumber; }

private:
   const MnUserParameterState *fUserState;
   int fParNumber;
};

} // namespace Minuit2
} // namespace ROOT

#endif // ROOT_Minuit2_MnTraceIter
