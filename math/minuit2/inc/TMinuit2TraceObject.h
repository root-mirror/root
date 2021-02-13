// @(#)root/minuit2:$Id$
// Author:  L. Moneta 2012

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2012 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_TMinuit2TraceObject
#define ROOT_TMinuit2TraceObject

#include "TNamed.h"
#include "Minuit2/MnTraceObject.h"

class TH1;
class TVirtualPad;
class TList;

namespace ROOT {

namespace Minuit2 {

class MinimumState;
class MnUserParameterState;

} // namespace Minuit2
} // namespace ROOT

class TMinuit2TraceObject : public ROOT::Minuit2::MnTraceObject, public TNamed {

public:
   TMinuit2TraceObject(int parNumber = -1);

   virtual ~TMinuit2TraceObject();

   virtual void Init(const ROOT::Minuit2::MnUserParameterState &state);

   virtual void operator()(int i, const ROOT::Minuit2::MinimumState &state);

   ClassDef(TMinuit2TraceObject, 0) // Example Trace Object for Minuit2

      private :

      int fIterOffset;      // offset in iteration in case of combined minimizers
   TH1 *fHistoFval;         // Function value histogram
   TH1 *fHistoEdm;          // Edm histogram
   TList *fHistoParList;    // list of parameter values histograms
   TVirtualPad *fOldPad;    // old existing current pad
   TVirtualPad *fMinuitPad; // new pad with trace histograms
};

#endif // ROOT_TMinuit2TraceObject
