/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2016-2020, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
#ifndef ROOT_ROOFIT_TESTSTATISTICS_RooRealL
#define ROOT_ROOFIT_TESTSTATISTICS_RooRealL

#include "Rtypes.h"  // ClassDef, ClassImp
#include <memory>  // shared_ptr
#include <RooAbsReal.h>
#include "RooSetProxy.h"

namespace RooFit {
namespace TestStatistics {

// forward declaration
class RooAbsL;

class RooRealL : public RooAbsReal {
public:
   RooRealL(const char *name, const char *title, std::shared_ptr<RooAbsL> likelihood);
   RooRealL(const RooRealL& other, const char* name=0);

   // pure virtual overrides:
   Double_t evaluate() const override;
   TObject* clone(const char* newname) const override;
   // virtual overrides:
   double globalNormalization() const;

   double get_carry() const;

   Double_t defaultErrorLevel() const override;
private:
   std::shared_ptr<RooAbsL> likelihood_;
   mutable double eval_carry = 0;
   RooSetProxy vars_proxy_;    // sets up client-server connections

   ClassDefOverride(RooRealL, 0);
};

}
}

#endif // ROOT_ROOFIT_TESTSTATISTICS_RooRealL
