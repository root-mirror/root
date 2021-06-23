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

#include <TestStatistics/optional_parameter_types.h>

namespace RooFit {
namespace TestStatistics {

ConstrainedParameters::ConstrainedParameters(const RooArgSet &parameters) : set(parameters) {}

ExternalConstraints::ExternalConstraints(const RooArgSet &constraints) : set(constraints) {}

GlobalObservables::GlobalObservables(const RooArgSet &observables) : set(observables) {}

}
}