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
#ifndef ROOT_ROOFIT_TESTSTATISTICS_RooSubsidiaryL
#define ROOT_ROOFIT_TESTSTATISTICS_RooSubsidiaryL

#include <TestStatistics/RooAbsL.h>
#include <RooArgList.h>
#include <RooArgSet.h>

namespace RooFit {
namespace TestStatistics {

/// Gathers all subsidiary PDF terms from the component PDFs of RooSumL likelihoods.
/// These are summed separately for increased numerical stability, since these terms are often
/// small and cause numerical variances in their original PDFs, whereas by summing as one
/// separate subsidiary collective term, it is numerically very stable.
/// Note that when a subsidiary PDF is part of multiple component PDFs, it will only be summed
/// once in this class! This doesn't change the derivative of the log likelihood (which is what
/// matters in fitting the likelihood), but does change the value of the (log-)likelihood itself.
class RooSubsidiaryL : public RooAbsL {
public:
   RooSubsidiaryL(const std::string & parent_pdf_name, const RooArgSet & pdfs, const RooArgSet & parameter_set);

   double evaluate_partition(Section events, std::size_t components_begin, std::size_t components_end) override;
   RooArgSet * getParameters() override;
   std::string GetName() const override;
   std::string GetTitle() const override;
   std::size_t numDataEntries() const override;

   void constOptimizeTestStatistic(RooAbsArg::ConstOpCode opcode, bool doAlsoTrackingOpt) override;

private:
   std::string parent_pdf_name_;
   RooArgList subsidiary_pdfs_{"subsidiary_pdfs"};    // Set of subsidiary PDF or "constraint" terms
   RooArgSet parameter_set_{"parameter_set"};       // Set of parameters to which constraints apply
};

}
}

#endif // ROOT_ROOFIT_TESTSTATISTICS_RooSubsidiaryL
