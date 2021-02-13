// @(#)root/hist:$Id$
// Author: Christophe.Delaere@cern.ch   21/08/2002

#ifndef ROOT_TLimit
#define ROOT_TLimit

#include "Rtypes.h"

#include "TVectorDfwd.h"

class TConfidenceLevel;
class TRandom;
class TLimitDataSource;
class TArrayD;
class TOrdCollection;
class TH1;
class TObjArray;

class TLimit {
 protected:
   static bool Fluctuate(TLimitDataSource * input, TLimitDataSource * output, bool init,TRandom *, bool stat=false);
   static Double_t LogLikelihood(Double_t s, Double_t b, Double_t b2, Double_t d);

public:
   TLimit() {}
   virtual ~TLimit() {}
   static TConfidenceLevel *ComputeLimit(TLimitDataSource * data,
                                         Int_t nmc =50000,
                                         bool stat = false,
                                         TRandom * generator = 0);
   static TConfidenceLevel *ComputeLimit(Double_t s, Double_t b, Int_t d,
                                         Int_t nmc =50000,
                                         bool stat = false,
                                         TRandom * generator = 0);
   static TConfidenceLevel *ComputeLimit(Double_t s, Double_t b, Int_t d,
                                         TVectorD* se, TVectorD* be, TObjArray*,
                                         Int_t nmc =50000,
                                         bool stat = false,
                                         TRandom * generator = 0);
   static TConfidenceLevel *ComputeLimit(TH1* s, TH1* b, TH1* d,
                                         Int_t nmc =50000,
                                         bool stat = false,
                                         TRandom * generator = 0);
   static TConfidenceLevel *ComputeLimit(TH1* s, TH1* b, TH1* d,
                                         TVectorD* se, TVectorD* be, TObjArray*,
                                         Int_t nmc =50000,
                                         bool stat = false,
                                         TRandom * generator = 0);
 private:
   static TArrayD *fgTable;              // a log table... just to speed up calculation
   static TOrdCollection *fgSystNames;   // Collection of systematics names
   ClassDef(TLimit, 2)          // Class to compute 95% CL limits
};

#endif

