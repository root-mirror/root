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

#include <TestStatistics/optimization.h>

#include <RooMsgService.h>
#include <RooVectorDataStore.h> // complete type for dynamic cast
#include <RooAbsReal.h>
#include <RooArgSet.h>
#include <RooAbsData.h>

namespace RooFit {
namespace TestStatistics {

RooArgSet ConstantTermsOptimizer::requiredExtraObservables()
{
   // TODO: the RooAbsOptTestStatistics::requiredExtraObservables() call this code was copied
   //       from was overloaded for RooXYChi2Var only; implement different options when necessary
   return RooArgSet();
}

void ConstantTermsOptimizer::enable_constant_terms_optimization(RooAbsReal *function, RooArgSet *norm_set, RooAbsData *dataset,
                                        bool applyTrackingOpt)
{
   // Trigger create of all object caches now in nodes that have deferred object creation
   // so that cache contents can be processed immediately
   function->getVal(norm_set);

   // Apply tracking optimization here. Default strategy is to track components
   // of RooAddPdfs and RooRealSumPdfs. If these components are a RooProdPdf
   // or a RooProduct respectively, track the components of these products instead
   // of the product term
   RooArgSet trackNodes;

   // Add safety check here - applyTrackingOpt will only be applied if present
   // dataset is constructed in terms of a RooVectorDataStore
   if (applyTrackingOpt) {
      if (!dynamic_cast<RooVectorDataStore *>(dataset->store())) {
         oocoutW((TObject *)nullptr, Optimization)
            << "enable_constant_terms_optimization(function: " << function->GetName()
            << ", dataset: " << dataset->GetName()
            << ") WARNING Cache-and-track optimization (Optimize level 2) is only available for datasets"
            << " implemented in terms of RooVectorDataStore - ignoring this option for current dataset" << std::endl;
         applyTrackingOpt = kFALSE;
      }
   }

   if (applyTrackingOpt) {
      RooArgSet branches;
      function->branchNodeServerList(&branches);
      RooFIter iter = branches.fwdIterator();
      RooAbsArg *arg;
      while ((arg = iter.next())) {
         arg->setCacheAndTrackHints(trackNodes);
      }
      // Do not set CacheAndTrack on constant expressions
      auto constNodes = (RooArgSet *)trackNodes.selectByAttrib("Constant", kTRUE);
      trackNodes.remove(*constNodes);
      delete constNodes;

      // Set CacheAndTrack flag on all remaining nodes
      trackNodes.setAttribAll("CacheAndTrack", kTRUE);
   }

   // Find all nodes that depend exclusively on constant parameters
   RooArgSet cached_nodes;

   function->findConstantNodes(*dataset->get(), cached_nodes);

   // Cache constant nodes with dataset - also cache entries corresponding to zero-weights in data when using
   // BinnedLikelihood
   // NOTE: we pass nullptr as cache-owner here, because we don't use the cacheOwner() anywhere in TestStatistics
   // TODO: make sure this (nullptr) is always correct
   dataset->cacheArgs(nullptr, cached_nodes, norm_set, !function->getAttribute("BinnedLikelihood"));

   // Put all cached nodes in AClean value caching mode so that their evaluate() is never called
   TIterator *cIter = cached_nodes.createIterator();
   RooAbsArg *cacheArg;
   while ((cacheArg = (RooAbsArg *)cIter->Next())) {
      cacheArg->setOperMode(RooAbsArg::AClean);
   }
   delete cIter;

   RooArgSet *constNodes = (RooArgSet *)cached_nodes.selectByAttrib("ConstantExpressionCached", kTRUE);
   RooArgSet actualTrackNodes(cached_nodes);
   actualTrackNodes.remove(*constNodes);
   if (constNodes->getSize() > 0) {
      if (constNodes->getSize() < 20) {
         oocoutI((TObject*)nullptr, Minimization)
            << " The following expressions have been identified as constant and will be precalculated and cached: "
            << *constNodes << std::endl;
      } else {
         oocoutI((TObject*)nullptr, Minimization) << " A total of " << constNodes->getSize()
                             << " expressions have been identified as constant and will be precalculated and cached."
                             << std::endl;
      }
   }
   if (actualTrackNodes.getSize() > 0) {
      if (actualTrackNodes.getSize() < 20) {
         oocoutI((TObject*)nullptr, Minimization) << " The following expressions will be evaluated in cache-and-track mode: "
                             << actualTrackNodes << std::endl;
      } else {
         oocoutI((TObject*)nullptr, Minimization) << " A total of " << constNodes->getSize()
                             << " expressions will be evaluated in cache-and-track-mode." << std::endl;
      }
   }
   delete constNodes;

   // Disable reading of observables that are no longer used
   dataset->optimizeReadingWithCaching(*function, cached_nodes, requiredExtraObservables());

   //   _optimized = kTRUE;
}

void ConstantTermsOptimizer::disable_constant_terms_optimization(RooAbsReal *function, RooArgSet *norm_set, RooArgSet *observables,
                                                                 RooAbsData *dataset)
{
   // Delete the cache
   dataset->resetCache();

   // Reactivate all tree branches
   dataset->setArgStatus(*dataset->get(), kTRUE);

   // Reset all nodes to ADirty
   optimize_caching(function, norm_set, observables, dataset);

   // Disable propagation of dirty state flags for observables
   dataset->setDirtyProp(kFALSE);

   //   _cachedNodes.removeAll();

   //   _optimized = kFALSE;
}

void ConstantTermsOptimizer::optimize_caching(RooAbsReal *function, RooArgSet *norm_set, RooArgSet *observables, RooAbsData *dataset)
{
   // Trigger create of all object caches now in nodes that have deferred object creation
   // so that cache contents can be processed immediately
   function->getVal(norm_set);

   // Set value caching mode for all nodes that depend on any of the observables to ADirty
   function->optimizeCacheMode(*observables);

   // Disable propagation of dirty state flags for observables
   dataset->setDirtyProp(kFALSE);

   // Disable reading of observables that are not used
   dataset->optimizeReadingWithCaching(*function, RooArgSet(), requiredExtraObservables()) ;
}

} // namespace TestStatistics
} // namespace RooFit