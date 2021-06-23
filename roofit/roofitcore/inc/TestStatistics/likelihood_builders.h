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
#ifndef ROOT_ROOFIT_TESTSTATISTICS_likelihood_builders
#define ROOT_ROOFIT_TESTSTATISTICS_likelihood_builders

#include <TestStatistics/RooAbsL.h>
#include <TestStatistics/optional_parameter_types.h>

// forward declarations
class RooAbsPdf;
class RooAbsData;

namespace RooFit {
namespace TestStatistics {

std::shared_ptr<RooAbsL> build_simultaneous_likelihood(RooAbsPdf *pdf, RooAbsData *data, RooAbsL::Extended extended = RooAbsL::Extended::Auto,
                                                       ConstrainedParameters constrained_parameters = {}, ExternalConstraints external_constraints = {},
                                                       GlobalObservables global_observables = {}, std::string global_observables_tag = {});

// delegating builder calls, for more convenient "optional" parameter passing
std::shared_ptr<RooAbsL> build_simultaneous_likelihood(RooAbsPdf* pdf, RooAbsData* data, ConstrainedParameters constrained_parameters);
std::shared_ptr<RooAbsL> build_simultaneous_likelihood(RooAbsPdf* pdf, RooAbsData* data, ExternalConstraints external_constraints);
std::shared_ptr<RooAbsL> build_simultaneous_likelihood(RooAbsPdf* pdf, RooAbsData* data, GlobalObservables global_observables);
std::shared_ptr<RooAbsL> build_simultaneous_likelihood(RooAbsPdf* pdf, RooAbsData* data, std::string global_observables_tag);
std::shared_ptr<RooAbsL> build_simultaneous_likelihood(RooAbsPdf* pdf, RooAbsData* data, ConstrainedParameters constrained_parameters, GlobalObservables global_observables);

std::shared_ptr<RooAbsL> build_unbinned_constrained_likelihood(RooAbsPdf *pdf, RooAbsData *data, RooAbsL::Extended extended = RooAbsL::Extended::Auto,
                                                               ConstrainedParameters constrained_parameters = {}, ExternalConstraints external_constraints = {},
                                                               GlobalObservables global_observables = {}, std::string global_observables_tag = {});

// delegating builder calls, for more convenient "optional" parameter passing
std::shared_ptr<RooAbsL> build_unbinned_constrained_likelihood(RooAbsPdf* pdf, RooAbsData* data, ConstrainedParameters constrained_parameters);
std::shared_ptr<RooAbsL> build_unbinned_constrained_likelihood(RooAbsPdf* pdf, RooAbsData* data, ExternalConstraints external_constraints);
std::shared_ptr<RooAbsL> build_unbinned_constrained_likelihood(RooAbsPdf* pdf, RooAbsData* data, GlobalObservables global_observables);
std::shared_ptr<RooAbsL> build_unbinned_constrained_likelihood(RooAbsPdf* pdf, RooAbsData* data, std::string global_observables_tag);

}
}

#endif // ROOT_ROOFIT_TESTSTATISTICS_likelihood_builders
