/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2016-2021, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
#ifndef ROOT_ROOFIT_TESTSTATISTICS_optional_parameter_types
#define ROOT_ROOFIT_TESTSTATISTICS_optional_parameter_types

#include <RooArgSet.h>

namespace RooFit {
namespace TestStatistics {

// strongly named container types for use as optional parameters in the test statistics constructors

struct ConstrainedParameters {
   ConstrainedParameters() = default;
   explicit ConstrainedParameters(const RooArgSet &parameters);
   RooArgSet set;
};

struct ExternalConstraints {
   ExternalConstraints() = default;
   explicit ExternalConstraints(const RooArgSet &constraints);
   RooArgSet set;
};

struct GlobalObservables {
   GlobalObservables() = default;
   explicit GlobalObservables(const RooArgSet &observables);
   RooArgSet set;
};

}
}

#endif // ROOT_ROOFIT_TESTSTATISTICS_optional_parameter_types
