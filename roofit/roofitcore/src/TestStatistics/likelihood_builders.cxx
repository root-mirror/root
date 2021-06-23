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

#include <TestStatistics/likelihood_builders.h>

#include <RooSimultaneous.h>
#include <TestStatistics/RooBinnedL.h>
#include <TestStatistics/RooUnbinnedL.h>
#include <TestStatistics/RooSubsidiaryL.h>
#include <TestStatistics/RooSumL.h>
#include <RooAbsPdf.h>
#include <RooAbsData.h>
#include <RooRealSumPdf.h>
#include <RooProdPdf.h>
#include <ROOT/RMakeUnique.hxx>
#include <memory>


namespace RooFit {
namespace TestStatistics {

std::unique_ptr<RooSubsidiaryL> build_constraints(RooAbsPdf *pdf, RooAbsData *data,
                                                  ConstrainedParameters constrained_parameters, ExternalConstraints external_constraints,
                                                  GlobalObservables global_observables, std::string global_observables_tag)
{
   // BEGIN CONSTRAINT COLLECTION; copied from RooAbsPdf::createNLL

   Bool_t doStripDisconnected = kFALSE;
   // If no explicit list of parameters to be constrained is specified apply default algorithm
   // All terms of RooProdPdfs that do not contain observables and share a parameters with one or more
   // terms that do contain observables are added as constraints.
#ifndef NDEBUG
   bool did_default_constraint_algo = false;
   std::size_t N_default_constraints = 0;
#endif
   if (constrained_parameters.set.getSize() == 0) {
      std::unique_ptr<RooArgSet> default_constraints{pdf->getParameters(*data, kFALSE)};
      constrained_parameters.set.add(*default_constraints);
      doStripDisconnected = kTRUE;
#ifndef NDEBUG
      did_default_constraint_algo = true;
      N_default_constraints = default_constraints->getSize();
#endif
   }
#ifndef NDEBUG
   if (did_default_constraint_algo) {
      assert(N_default_constraints == constrained_parameters.set.getSize());
   }
#endif

   // Collect internal and external constraint specifications
   RooArgSet allConstraints;

   if (!global_observables_tag.empty()) {
      std::cout << "DEBUG: global_observables_tag > 0" << std::endl;
      if (global_observables.set.getSize() > 0) {
         global_observables.set.removeAll();
      }
      std::unique_ptr<RooArgSet> allVars {pdf->getVariables()};
      global_observables.set.add(*dynamic_cast<RooArgSet *>(allVars->selectByAttrib(global_observables_tag.c_str(), kTRUE)));
      oocoutI((TObject*)nullptr, Minimization) << "User-defined specification of global observables definition with tag named '" <<  global_observables_tag << "'" << std::endl;
   } else if (global_observables.set.getSize() == 0) {
      // neither global_observables nor global_observables_tag was given - try if a default tag is defined in the head node
      const char* defGlobObsTag = pdf->getStringAttribute("DefaultGlobalObservablesTag");
      if (defGlobObsTag) {
         oocoutI((TObject*)nullptr, Minimization) << "p.d.f. provides built-in specification of global observables definition with tag named '" <<  defGlobObsTag << "'" << std::endl;
         std::unique_ptr<RooArgSet> allVars {pdf->getVariables()};
         global_observables.set.add(*dynamic_cast<RooArgSet *>(allVars->selectByAttrib(defGlobObsTag, kTRUE)));
      }
   }

   // EGP: removed workspace (RooAbsPdf::_myws) based stuff for now; TODO: reconnect this class to workspaces

   if (constrained_parameters.set.getSize() > 0) {
      std::unique_ptr<RooArgSet> constraints{pdf->getAllConstraints(*data->get(), constrained_parameters.set, doStripDisconnected)};
      allConstraints.add(*constraints);
   }
   if (external_constraints.set.getSize() > 0) {
      allConstraints.add(external_constraints.set);
   }

   std::unique_ptr<RooSubsidiaryL> subsidiary_likelihood;
   // Include constraints, if any, in likelihood
   if (allConstraints.getSize() > 0) {

      oocoutI((TObject*) nullptr, Minimization) << " Including the following contraint terms in minimization: " << allConstraints << std::endl;
      if (global_observables.set.getSize() > 0) {
         oocoutI((TObject*) nullptr, Minimization) << "The following global observables have been defined: " << global_observables.set << std::endl;
      }
      std::string name("likelihood for pdf ");
      name += pdf->GetName();
      subsidiary_likelihood = std::make_unique<RooSubsidiaryL>(name, allConstraints,
                                                               (global_observables.set.getSize() > 0) ? global_observables.set : constrained_parameters.set);
   }

   // END CONSTRAINT COLLECTION; copied from RooAbsPdf::createNLL

   return subsidiary_likelihood;
}


std::shared_ptr<RooAbsL> build_simultaneous_likelihood(RooAbsPdf *pdf, RooAbsData *data, RooAbsL::Extended extended,
                                                       ConstrainedParameters constrained_parameters, ExternalConstraints external_constraints,
                                                       GlobalObservables global_observables, std::string global_observables_tag) {
   auto sim_pdf = dynamic_cast<RooSimultaneous *>(pdf);
   if (sim_pdf == nullptr) {
      throw std::logic_error("Can only build RooSumL from RooSimultaneous pdf!");
   }

   // the rest of this function is an adaptation of RooAbsTestStatistic::initSimMode:

   RooAbsCategoryLValue &simCat = (RooAbsCategoryLValue &)sim_pdf->indexCat();

   // FIXME: process_empty_data_sets was previously member function processEmptyDataSets(), which returned extended_ for RooSumL (previously named RooSimultaneousL). See RooAbsL and RooSumL/RooUnbinnedL implementations. I think this should work here, but be careful in other subclasses!
   bool process_empty_data_sets = RooAbsL::is_extended(pdf, extended);

   TString simCatName(simCat.GetName());
   // Note: important not to use cloned dataset here (possible when this code is run in Roo[...]L ctor), use the
   // original one (which is data_ in Roo[...]L ctors, but data here)
   std::unique_ptr<TList> dsetList{data->split(simCat, process_empty_data_sets)};
   if (!dsetList) {
      throw std::logic_error("build_simultaneous_likelihood ERROR, index category of simultaneous pdf is missing in dataset, aborting");
   }

   // Count number of used states
   std::size_t N_components = 0;

   RooCatType *type;
   std::unique_ptr<TIterator> catIter{simCat.typeIterator()};
   while ((type = (RooCatType *)catIter->Next())) {
      // Retrieve the PDF for this simCat state
      RooAbsPdf *component_pdf = sim_pdf->getPdf(type->GetName());
      auto dset = (RooAbsData *)dsetList->FindObject(type->GetName());

      if (component_pdf && dset && (0. != dset->sumEntries() || process_empty_data_sets)) {
         ++N_components;
      }
   }

   // Allocate arrays
   std::vector<std::unique_ptr<RooAbsL>> components;
   components.reserve(N_components);
   //   _gofSplitMode.resize(N_components);  // not used, Hybrid mode only, see below

   // Create array of regular fit contexts, containing subset of data and single fitCat PDF
   catIter->Reset();
   std::size_t n = 0;
   while ((type = (RooCatType *)catIter->Next())) {
      // Retrieve the PDF for this simCat state
      RooAbsPdf *component_pdf = sim_pdf->getPdf(type->GetName());
      auto dset = (RooAbsData *)dsetList->FindObject(type->GetName());

      if (component_pdf && dset && (0. != dset->sumEntries() || process_empty_data_sets)) {
         ooccoutI((TObject *)nullptr, Fitting)
         << "RooSumL: creating slave calculator #" << n << " for state " << type->GetName() << " ("
                                                   << dset->numEntries() << " dataset entries)" << std::endl;

         // *** START HERE
         // WVE HACK determine if we have a RooRealSumPdf and then treat it like a binned likelihood
         RooAbsPdf *binnedPdf = 0;
         Bool_t binnedL = kFALSE;
         if (component_pdf->getAttribute("BinnedLikelihood") &&
             component_pdf->IsA()->InheritsFrom(RooRealSumPdf::Class())) {
            // Simplest case: top-level of component is a RRSP
            binnedPdf = component_pdf;
            binnedL = kTRUE;
         } else if (component_pdf->IsA()->InheritsFrom(RooProdPdf::Class())) {
            // Default case: top-level pdf is a product of RRSP and other pdfs
            RooFIter iter = ((RooProdPdf *)component_pdf)->pdfList().fwdIterator();
            RooAbsArg *component;
            while ((component = iter.next())) {
               if (component->getAttribute("BinnedLikelihood") &&
                   component->IsA()->InheritsFrom(RooRealSumPdf::Class())) {
                  binnedPdf = (RooAbsPdf *)component;
                  binnedL = kTRUE;
               }
               if (component->getAttribute("MAIN_MEASUREMENT")) {
                  // not really a binned pdf, but this prevents a (potentially) long list of subsidiary measurements to
                  // be passed to the slave calculator
                  binnedPdf = (RooAbsPdf *)component;
               }
            }
         }
         // WVE END HACK
         // Below here directly pass binnedPdf instead of PROD(binnedPdf,constraints) as constraints are evaluated
         // elsewhere anyway and omitting them reduces model complexity and associated handling/cloning times
         //         if (_splitRange && rangeName) {
         //            _gofArray[n] =
         //               create(type->GetName(), type->GetName(), (binnedPdf ? *binnedPdf : *component_pdf), *dset,
         //               *projDeps,
         //                      Form("%s_%s", rangeName, type->GetName()), addCoefRangeName, _nCPU * (_mpinterl ? -1 :
         //                      1), _mpinterl, _CPUAffinity, _verbose, _splitRange, binnedL);
         //         } else {
         //            _gofArray[n] =
         //               create(type->GetName(), type->GetName(), (binnedPdf ? *binnedPdf : *component_pdf), *dset,
         //               *projDeps, rangeName,
         //                      addCoefRangeName, _nCPU, _mpinterl, _CPUAffinity, _verbose, _splitRange, binnedL);
         if (binnedL) {
            components.push_back(std::make_unique<RooBinnedL>((binnedPdf ? binnedPdf : component_pdf), dset));
         } else {
            components.push_back(std::make_unique<RooUnbinnedL>((binnedPdf ? binnedPdf : component_pdf), dset));
         }
         //         }
         components.back()->set_sim_count(N_components);
         // *** END HERE

         // TODO: left out Hybrid mode for now, evaluate later whether to reinclude (also then change
         // evaluate_partition)
         //         // Fill per-component split mode with Bulk Partition for now so that Auto will map to bulk-splitting
         //         of all components if (_mpinterl==RooFit::Hybrid) {
         //            if (dset->numEntries()<10) {
         //               //cout << "RAT::initSim("<< GetName() << ") MP mode is auto, setting split mode for component
         //               "<< n << " to SimComponents"<< endl ; _gofSplitMode[n] = RooFit::SimComponents;
         //               _gofArray[n]->_mpinterl = RooFit::SimComponents;
         //            } else {
         //               //cout << "RAT::initSim("<< GetName() << ") MP mode is auto, setting split mode for component
         //               "<< n << " to BulkPartition"<< endl ; _gofSplitMode[n] = RooFit::BulkPartition;
         //               _gofArray[n]->_mpinterl = RooFit::BulkPartition;
         //            }
         //         }
         //
         // Servers may have been redirected between instantiation and (deferred) initialization

         std::unique_ptr<RooArgSet> actualParams{binnedPdf ? binnedPdf->getParameters(dset)
                                                           : component_pdf->getParameters(dset)};
         std::unique_ptr<RooArgSet> selTargetParams{(RooArgSet *)pdf->getParameters(*data)->selectCommon(*actualParams)};

         // TODO: I don't think we have to redirect servers, because our classes make no use of those, but we should
         // make sure. Do we need to reset the parameter set instead?
         //         components_.back()->recursiveRedirectServers(*selTargetParams);
         assert(selTargetParams->equals(*components.back()->getParameters()));

         ++n;
      } else {
         if ((!dset || (0. != dset->sumEntries() && !process_empty_data_sets)) && component_pdf) {
            ooccoutD((TObject *)nullptr, Fitting) << "RooSumL: state " << type->GetName()
                                                                       << " has no data entries, no slave calculator created" << std::endl;
         }
      }
   }
   oocoutI((TObject *)nullptr, Fitting) << "RooSumL: created " << n << " slave calculators." << std::endl;

   std::unique_ptr<RooAbsL> subsidiary = build_constraints(pdf, data, constrained_parameters, external_constraints, global_observables, global_observables_tag);
   if (subsidiary) {
      components.push_back(std::move(subsidiary));
   }

   return std::make_shared<RooSumL>(pdf, data, std::move(components), extended);
}

// delegating convenience overloads
std::shared_ptr<RooAbsL> build_simultaneous_likelihood(RooAbsPdf* pdf, RooAbsData* data, ConstrainedParameters constrained_parameters)
{
   return build_simultaneous_likelihood(pdf, data, RooAbsL::Extended::Auto, constrained_parameters);
}
std::shared_ptr<RooAbsL> build_simultaneous_likelihood(RooAbsPdf* pdf, RooAbsData* data, ExternalConstraints external_constraints)
{
   return build_simultaneous_likelihood(pdf, data, RooAbsL::Extended::Auto, {}, external_constraints);
}
std::shared_ptr<RooAbsL> build_simultaneous_likelihood(RooAbsPdf* pdf, RooAbsData* data, GlobalObservables global_observables)
{
   return build_simultaneous_likelihood(pdf, data, RooAbsL::Extended::Auto, {}, {}, global_observables);
}
std::shared_ptr<RooAbsL> build_simultaneous_likelihood(RooAbsPdf* pdf, RooAbsData* data, std::string global_observables_tag)
{
   return build_simultaneous_likelihood(pdf, data, RooAbsL::Extended::Auto, {}, {}, {}, global_observables_tag);
}
std::shared_ptr<RooAbsL> build_simultaneous_likelihood(RooAbsPdf* pdf, RooAbsData* data, ConstrainedParameters constrained_parameters, GlobalObservables global_observables)
{
   return build_simultaneous_likelihood(pdf, data, RooAbsL::Extended::Auto, constrained_parameters, {}, global_observables);
}


std::shared_ptr<RooAbsL> build_unbinned_constrained_likelihood(RooAbsPdf *pdf, RooAbsData *data, RooAbsL::Extended extended,
                                                               ConstrainedParameters constrained_parameters, ExternalConstraints external_constraints,
                                                               GlobalObservables global_observables, std::string global_observables_tag) {
   std::vector<std::unique_ptr<RooAbsL>> components;
   components.reserve(2);
   components.push_back(std::make_unique<RooUnbinnedL>(pdf, data, extended));
   components.push_back(build_constraints(pdf, data, constrained_parameters, external_constraints, global_observables, global_observables_tag));
   return std::make_shared<RooSumL>(pdf, data, std::move(components), extended);
}

// delegating convenience overloads
std::shared_ptr<RooAbsL> build_unbinned_constrained_likelihood(RooAbsPdf* pdf, RooAbsData* data, ConstrainedParameters constrained_parameters)
{
   return build_unbinned_constrained_likelihood(pdf, data, RooAbsL::Extended::Auto, constrained_parameters);
}
std::shared_ptr<RooAbsL> build_unbinned_constrained_likelihood(RooAbsPdf* pdf, RooAbsData* data, ExternalConstraints external_constraints)
{
   return build_unbinned_constrained_likelihood(pdf, data, RooAbsL::Extended::Auto, {}, external_constraints);
}
std::shared_ptr<RooAbsL> build_unbinned_constrained_likelihood(RooAbsPdf* pdf, RooAbsData* data, GlobalObservables global_observables)
{
   return build_unbinned_constrained_likelihood(pdf, data, RooAbsL::Extended::Auto, {}, {}, global_observables);
}
std::shared_ptr<RooAbsL> build_unbinned_constrained_likelihood(RooAbsPdf* pdf, RooAbsData* data, std::string global_observables_tag)
{
   return build_unbinned_constrained_likelihood(pdf, data, RooAbsL::Extended::Auto, {}, {}, {}, global_observables_tag);
}

}
}

