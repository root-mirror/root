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
#include <TestStatistics/RooSubsidiaryL.h>
#include <TestStatistics/kahan_sum.h>
#include <RooAbsPdf.h> // for dynamic cast
#include <RooErrorHandler.h>

namespace RooFit {
namespace TestStatistics {

RooSubsidiaryL::RooSubsidiaryL(const std::string& parent_pdf_name, const RooArgSet &pdfs, const RooArgSet &parameter_set)
   : RooAbsL(nullptr, nullptr, 0, 0, RooAbsL::Extended::No), parent_pdf_name_(parent_pdf_name)
{
   std::unique_ptr<TIterator> inputIter{pdfs.createIterator()};
   RooAbsArg *comp;
   while ((comp = (RooAbsArg *)inputIter->Next())) {
      if (!dynamic_cast<RooAbsPdf *>(comp)) {
         oocoutE((TObject *)0, InputArguments) << "RooSubsidiaryL::ctor(" << GetName() << ") ERROR: component "
                                               << comp->GetName() << " is not of type RooAbsPdf" << std::endl;
         RooErrorHandler::softAbort();
      }
      subsidiary_pdfs_.add(*comp);
   }
   parameter_set_.add(parameter_set);
}

double
RooSubsidiaryL::evaluate_partition(RooAbsL::Section events, std::size_t /*components_begin*/, std::size_t /*components_end*/)
{
   if (events.begin_fraction != 0 || events.end_fraction != 1) {
      oocoutW((TObject *)0, InputArguments) << "RooSubsidiaryL::evaluate_partition can only calculate everything, so "
                                               "section should be {0,1}, but it's not!"
                                            << std::endl;
   }

   double sum = 0, carry = 0;
   RooAbsReal *comp;
   RooFIter setIter1 = subsidiary_pdfs_.fwdIterator();

   while ((comp = (RooAbsReal *)setIter1.next())) {
      double term = -((RooAbsPdf *)comp)->getLogVal(&parameter_set_);
      std::tie(sum, carry) = kahan_add(sum, term, carry);
   }
   eval_carry_ = carry;

   return sum;
}

RooArgSet *RooSubsidiaryL::getParameters()
{
   return &parameter_set_;
}
std::string RooSubsidiaryL::GetName() const
{
   return std::string("subsidiary_pdf_of_") + parent_pdf_name_;
}
std::string RooSubsidiaryL::GetTitle() const
{
   return std::string("Subsidiary PDF set of simultaneous PDF ") + parent_pdf_name_;
}
std::size_t RooSubsidiaryL::numDataEntries() const
{
   // function only used in LikelihoodJob::evaluate, but this class must always be evaluated over Section(0,1), so not useful
   return 0;
}

void RooSubsidiaryL::constOptimizeTestStatistic(RooAbsArg::ConstOpCode /*opcode*/, bool /*doAlsoTrackingOpt*/) {
   // do nothing, there's no dataset here to cache
}

} // namespace TestStatistics
} // namespace RooFit
