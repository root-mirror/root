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
#include <TestStatistics/RooAbsL.h>
#include <TestStatistics/optimization.h>
#include "RooAbsPdf.h"
#include "RooAbsData.h"

#include <ROOT/RMakeUnique.hxx>

// for dynamic casts in init_clones:
#include "RooAbsRealLValue.h"
#include "RooRealVar.h"
#include "RooDataHist.h"

// other stuff in init_clones:
#include "RooErrorHandler.h"
#include "RooMsgService.h"

// concrete classes in getParameters (testing, remove later)
#include "RooRealSumPdf.h"
#include "RooProdPdf.h"

namespace RooFit {
namespace TestStatistics {

// static function
bool RooAbsL::is_extended(RooAbsPdf* pdf, Extended extended)
{
   bool return_value;
   switch (extended) {
   case RooAbsL::Extended::No: {
      return_value = false;
      break;
   }
   case RooAbsL::Extended::Yes: {
      return_value = true;
      break;
   }
   case RooAbsL::Extended::Auto: {
      return_value = ((pdf->extendMode() == RooAbsPdf::CanBeExtended || pdf->extendMode() == RooAbsPdf::MustBeExtended))
                        ? true
                        : false;
      break;
   }
   }
   return return_value;
}

// private ctor
RooAbsL::RooAbsL(std::shared_ptr<RooAbsPdf> pdf, std::shared_ptr<RooAbsData> data,
                 std::size_t N_events, std::size_t N_components, Extended extended)
   : pdf_(std::move(pdf)), data_(std::move(data)), N_events(N_events), N_components(N_components)
{
   //   std::unique_ptr<RooArgSet> obs {pdf->getObservables(*data)};
   //   data->attachBuffers(*obs);
   extended_ = is_extended(pdf_.get(), extended);
   if (extended == Extended::Auto) {
      if (extended_) {
         oocoutI((TObject *)nullptr, Minimization)
            << "in RooAbsL ctor: p.d.f. provides expected number of events, including extended term in likelihood."
            << std::endl;
      }
   }
}

// this constructor clones the pdf/data and owns those cloned copies
RooAbsL::RooAbsL(RooAbsL::ClonePdfData in, std::size_t N_events, std::size_t N_components, Extended extended)
  : RooAbsL(std::shared_ptr<RooAbsPdf>(static_cast<RooAbsPdf *>(in.pdf->cloneTree())),
     std::shared_ptr<RooAbsData>(static_cast<RooAbsData *>(in.data->Clone())), N_events, N_components, extended)
{
   init_clones(*in.pdf, *in.data);
}

// this constructor does not clone pdf/data and uses the shared_ptr aliasing constructor to make it non-owning
RooAbsL::RooAbsL(RooAbsPdf *inpdf, RooAbsData *indata, std::size_t N_events, std::size_t N_components,
                 Extended extended)
   : RooAbsL({std::shared_ptr<RooAbsPdf>(nullptr), inpdf}, {std::shared_ptr<RooAbsData>(nullptr), indata}, N_events, N_components, extended)
{}


RooAbsL::RooAbsL(const RooAbsL &other)
   : pdf_(other.pdf_), data_(other.data_), N_events(other.N_events),
     N_components(other.N_components), extended_(other.extended_), sim_count_(other.sim_count_), eval_carry_(other.eval_carry_)
{
   // it can never be one, since we just copied the shared_ptr; if it is, something really weird is going on; also they must be equal (usually either zero or two)
   assert((pdf_.use_count() != 1) && (data_.use_count() != 1) && (pdf_.use_count() == data_.use_count()));
   // TODO: use aliasing ctor in initialization list, and then check in body here whether pdf and data were clones; if so, they need to be cloned again (and init_clones called on them)
   if ((pdf_.use_count() > 1) && (data_.use_count() > 1)) {
      pdf_.reset(static_cast<RooAbsPdf *>(other.pdf_->cloneTree()));
      data_.reset(static_cast<RooAbsData *>(other.data_->Clone()));
      init_clones(*other.pdf_, *other.data_);
   }
}

void RooAbsL::init_clones(RooAbsPdf &inpdf, RooAbsData &indata)
{
   //   RooArgSet obs(*indata.get()) ;
   //   obs.remove(projDeps,kTRUE,kTRUE) ;

   // ******************************************************************
   // *** PART 1 *** Clone incoming pdf, attach to each other *
   // ******************************************************************

   // moved to ctor
   //   pdf = static_cast<RooAbsPdf *>(inpdf->cloneTree());

   // Attach FUNC to data set
   auto _funcObsSet = pdf_->getObservables(indata);

   if (pdf_->getAttribute("BinnedLikelihood")) {
      pdf_->setAttribute("BinnedLikelihoodActive");
   }

   // Reattach FUNC to original parameters
   std::unique_ptr<RooArgSet> origParams{inpdf.getParameters(indata)};
   pdf_->recursiveRedirectServers(*origParams);

   // Mark all projected dependents as such
   //   if (projDeps.getSize()>0) {
   //      RooArgSet *projDataDeps = (RooArgSet*) _funcObsSet->selectCommon(projDeps) ;
   //      projDataDeps->setAttribAll("projectedDependent") ;
   //      delete projDataDeps ;
   //   }

   // TODO: do we need this here? Or in RooSumL?
   //   // If PDF is a RooProdPdf (with possible constraint terms)
   //   // analyze pdf for actual parameters (i.e those in unconnected constraint terms should be
   //   // ignored as here so that the test statistic will not be recalculated if those
   //   // are changed
   //   RooProdPdf* pdfWithCons = dynamic_cast<RooProdPdf*>(pdf) ;
   //   if (pdfWithCons) {
   //
   //      RooArgSet* connPars = pdfWithCons->getConnectedParameters(*indata.get()) ;
   //      // Add connected parameters as servers
   //      _paramSet.removeAll() ;
   //      _paramSet.add(*connPars) ;
   //      delete connPars ;
   //
   //   } else {
   //      // Add parameters as servers
   //      _paramSet.add(*origParams) ;
   //   }

   // Store normalization set
   _normSet.reset((RooArgSet *)indata.get()->snapshot(kFALSE));

   // Expand list of observables with any observables used in parameterized ranges
   RooAbsArg *realDep;
   RooFIter iter = _funcObsSet->fwdIterator();
   while ((realDep = iter.next())) {
      RooAbsRealLValue *realDepRLV = dynamic_cast<RooAbsRealLValue *>(realDep);
      if (realDepRLV && realDepRLV->isDerived()) {
         RooArgSet tmp2;
         realDepRLV->leafNodeServerList(&tmp2, 0, kTRUE);
         _funcObsSet->add(tmp2, kTRUE);
      }
   }

   // ******************************************************************
   // *** PART 2 *** Clone and adjust incoming data, attach to PDF     *
   // ******************************************************************

   // Check if the fit ranges of the dependents in the data and in the FUNC are consistent
   const RooArgSet *dataDepSet = indata.get();
   iter = _funcObsSet->fwdIterator();
   RooAbsArg *arg;
   while ((arg = iter.next())) {

      // Check that both dataset and function argument are of type RooRealVar
      RooRealVar *realReal = dynamic_cast<RooRealVar *>(arg);
      if (!realReal) {
         continue;
      }
      RooRealVar *datReal = dynamic_cast<RooRealVar *>(dataDepSet->find(realReal->GetName()));
      if (!datReal) {
         continue;
      }

      // Check that range of observables in pdf is equal or contained in range of observables in data

      if (!realReal->getBinning().lowBoundFunc() && realReal->getMin() < (datReal->getMin() - 1e-6)) {
         oocoutE((TObject *)0, InputArguments) << "RooAbsL: ERROR minimum of FUNC observable " << arg->GetName() << "("
                                               << realReal->getMin() << ") is smaller than that of " << arg->GetName()
                                               << " in the dataset (" << datReal->getMin() << ")" << std::endl;
         RooErrorHandler::softAbort();
         return;
      }

      if (!realReal->getBinning().highBoundFunc() && realReal->getMax() > (datReal->getMax() + 1e-6)) {
         oocoutE((TObject *)0, InputArguments)
            << "RooAbsL: ERROR maximum of FUNC observable " << arg->GetName() << " is larger than that of "
            << arg->GetName() << " in the dataset" << std::endl;
         RooErrorHandler::softAbort();
         return;
      }
   }

   //   // Copy data and strip entries lost by adjusted fit range, data ranges will be copied from realDepSet ranges
   //   if (rangeName && strlen(rangeName)) {
   //      data = ((RooAbsData &)indata).reduce(RooFit::SelectVars(*_funcObsSet), RooFit::CutRange(rangeName));
   //      //     cout << "RooAbsOptTestStatistic: reducing dataset to fit in range named " << rangeName << " resulting
   //      //     dataset has " << data->sumEntries() << " events" << endl ;
   //   } else {

   // moved to ctor
   //      data = static_cast<RooAbsData *>(indata.Clone());

   //   }
   //   _ownData = kTRUE;

   // ******************************************************************
   // *** PART 3 *** Make adjustments for fit ranges, if specified     *
   // ******************************************************************

   //   RooArgSet *origObsSet = inpdf.getObservables(indata);
   //   RooArgSet *dataObsSet = (RooArgSet *)data->get();
   //   if (rangeName && strlen(rangeName)) {
   //      cxcoutI(Fitting) << "RooAbsOptTestStatistic::ctor(" << GetName()
   //                       << ") constructing test statistic for sub-range named " << rangeName << endl;
   //      // cout << "now adjusting observable ranges to requested fit range" << endl ;
   //
   //      // Adjust FUNC normalization ranges to requested fitRange, store original ranges for RooAddPdf coefficient
   //      // interpretation
   //      iter = _funcObsSet->fwdIterator();
   //      while ((arg = iter.next())) {
   //
   //         RooRealVar *realObs = dynamic_cast<RooRealVar *>(arg);
   //         if (realObs) {
   //
   //            // If no explicit range is given for RooAddPdf coefficients, create explicit named range equivalent to
   //            // original observables range
   //            if (!(addCoefRangeName && strlen(addCoefRangeName))) {
   //               realObs->setRange(Form("NormalizationRangeFor%s", rangeName), realObs->getMin(), realObs->getMax());
   //               // 	  cout << "RAOTS::ctor() setting range " << Form("NormalizationRangeFor%s",rangeName) << " on
   //               // observable "
   //               // 	       << realObs->GetName() << " to [" << realObs->getMin() << "," << realObs->getMax() << "]"
   //               <<
   //               // endl ;
   //            }
   //
   //            // Adjust range of function observable to those of given named range
   //            realObs->setRange(realObs->getMin(rangeName), realObs->getMax(rangeName));
   //            //  	cout << "RAOTS::ctor() setting normalization range on observable "
   //            //  	     << realObs->GetName() << " to [" << realObs->getMin() << "," << realObs->getMax() << "]" <<
   //            endl
   //            //  ;
   //
   //            // Adjust range of data observable to those of given named range
   //            RooRealVar *dataObs = (RooRealVar *)dataObsSet->find(realObs->GetName());
   //            dataObs->setRange(realObs->getMin(rangeName), realObs->getMax(rangeName));
   //
   //            // Keep track of list of fit ranges in string attribute fit range of original p.d.f.
   //            if (!_splitRange) {
   //               const char *origAttrib = inpdf.getStringAttribute("fitrange");
   //               if (origAttrib) {
   //                  inpdf.setStringAttribute("fitrange", Form("%s,fit_%s", origAttrib, GetName()));
   //               } else {
   //                  inpdf.setStringAttribute("fitrange", Form("fit_%s", GetName()));
   //               }
   //               RooRealVar *origObs = (RooRealVar *)origObsSet->find(arg->GetName());
   //               if (origObs) {
   //                  origObs->setRange(Form("fit_%s", GetName()), realObs->getMin(rangeName),
   //                  realObs->getMax(rangeName));
   //               }
   //            }
   //         }
   //      }
   //   }
   //   delete origObsSet;

   // If dataset is binned, activate caching of bins that are invalid because they're outside the
   // updated range definition (WVE need to add virtual interface here)
   RooDataHist *tmph = dynamic_cast<RooDataHist *>(data_.get());
   if (tmph) {
      tmph->cacheValidEntries();
   }

   //   // Fix RooAddPdf coefficients to original normalization range
   //   if (rangeName && strlen(rangeName)) {
   //
   //      // WVE Remove projected dependents from normalization
   //      pdf->fixAddCoefNormalization(*data->get(), kFALSE);
   //
   //      if (addCoefRangeName && strlen(addCoefRangeName)) {
   //         cxcoutI(Fitting) << "RooAbsOptTestStatistic::ctor(" << GetName()
   //                          << ") fixing interpretation of coefficients of any RooAddPdf component to range "
   //                          << addCoefRangeName << endl;
   //         pdf->fixAddCoefRange(addCoefRangeName, kFALSE);
   //      } else {
   //         cxcoutI(Fitting) << "RooAbsOptTestStatistic::ctor(" << GetName()
   //                          << ") fixing interpretation of coefficients of any RooAddPdf to full domain of
   //                          observables "
   //                          << endl;
   //         pdf->fixAddCoefRange(Form("NormalizationRangeFor%s", rangeName), kFALSE);
   //      }
   //   }

   // This is deferred from part 2 - but must happen after part 3 - otherwise invalid bins cannot be properly marked in
   // cacheValidEntries
   data_->attachBuffers(*_funcObsSet);
   // TODO: we pass event count to the ctor in the subclasses currently, because it's split into components and events
   // now
   //   setEventCount(data->numEntries());

   // *********************************************************************
   // *** PART 4 *** Adjust normalization range for projected observables *
   // *********************************************************************

   //   // Remove projected dependents from normalization set
   //   if (projDeps.getSize() > 0) {
   //
   //      _projDeps = (RooArgSet *)projDeps.snapshot(kFALSE);
   //
   //      // RooArgSet* tobedel = (RooArgSet*) _normSet->selectCommon(*_projDeps) ;
   //      _normSet->remove(*_projDeps, kTRUE, kTRUE);
   //
   //      //     // Delete owned projected dependent copy in _normSet
   //      //     TIterator* ii = tobedel->createIterator() ;
   //      //     RooAbsArg* aa ;
   //      //     while((aa=(RooAbsArg*)ii->Next())) {
   //      //       delete aa ;
   //      //     }
   //      //     delete ii ;
   //      //     delete tobedel ;
   //
   //      // Mark all projected dependents as such
   //      RooArgSet *projDataDeps = (RooArgSet *)_funcObsSet->selectCommon(*_projDeps);
   //      projDataDeps->setAttribAll("projectedDependent");
   //      delete projDataDeps;
   //   }

   //   coutI(Optimization)
   //      << "RooAbsOptTestStatistic::ctor(" << GetName()
   //      << ") optimizing internal clone of p.d.f for likelihood evaluation."
   //      << "Lazy evaluation and associated change tracking will disabled for all nodes that depend on observables"
   //      << endl;

   // *********************************************************************
   // *** PART 4 *** Finalization and activation of optimization          *
   // *********************************************************************

   //_origFunc = _func ;
   //_origData = _data ;

   //   // Redirect pointers of base class to clone
   //   _func = pdf ;
   //   _data = data ;

   // TODO: why this call?
   //   pdf->getVal(_normSet);

   //   cout << "ROATS::ctor(" << GetName() << ") funcClone structure dump BEFORE opt" << endl ;
   //   pdf->Print("t") ;

   //   optimizeCaching() ;

   //   cout << "ROATS::ctor(" << GetName() << ") funcClone structure dump AFTER opt" << endl ;
   //   pdf->Print("t") ;

   // optimization steps (copied from ROATS::optimizeCaching)

   pdf_->getVal(_normSet.get());
   // Set value caching mode for all nodes that depend on any of the observables to ADirty
   pdf_->optimizeCacheMode(*_funcObsSet);
   // Disable propagation of dirty state flags for observables
   data_->setDirtyProp(kFALSE);

   // Disable reading of observables that are not used
   data_->optimizeReadingWithCaching(*pdf_, RooArgSet(), RooArgSet()) ;

}

RooArgSet *RooAbsL::getParameters()
{
   auto ding = pdf_->getParameters(*data_);
   return ding;

//   // *** START HERE
//   // WVE HACK determine if we have a RooRealSumPdf and then treat it like a binned likelihood
//   RooAbsPdf *binnedPdf = 0;
//   if (pdf_->getAttribute("BinnedLikelihood") && pdf_->IsA()->InheritsFrom(RooRealSumPdf::Class())) {
//      // Simplest case: top-level of component is a RRSP
//      binnedPdf = pdf_.get();
//   } else if (pdf_->IsA()->InheritsFrom(RooProdPdf::Class())) {
//      // Default case: top-level pdf is a product of RRSP and other pdfs
//      RooFIter iter = ((RooProdPdf *)pdf_.get())->pdfList().fwdIterator();
//      RooAbsArg *component;
//      while ((component = iter.next())) {
//         if (component->getAttribute("BinnedLikelihood") &&
//             component->IsA()->InheritsFrom(RooRealSumPdf::Class())) {
//            binnedPdf = (RooAbsPdf *)component;
//         }
//         if (component->getAttribute("MAIN_MEASUREMENT")) {
//            // not really a binned pdf, but this prevents a (potentially) long list of subsidiary measurements to
//            // be passed to the slave calculator
//            binnedPdf = (RooAbsPdf *)component;
//         }
//      }
//   }
//   // WVE END HACK
//
//   std::unique_ptr<RooArgSet> actualParams {binnedPdf ? binnedPdf->getParameters(data_.get()) : pdf_->getParameters(data_.get())};
//   RooArgSet* selTargetParams = (RooArgSet *)ding->selectCommon(*actualParams);
//
//   std::cout << "RooAbsL::getParameters:" << std::endl;
//   selTargetParams->Print("v");
//
//   return selTargetParams;
}

void RooAbsL::constOptimizeTestStatistic(RooAbsArg::ConstOpCode opcode, bool doAlsoTrackingOpt)
{
   // to be further implemented, this is just a first test implementation
   if (opcode == RooAbsArg::Activate) {
      ConstantTermsOptimizer::enable_constant_terms_optimization(pdf_.get(), _normSet.get(), data_.get(), doAlsoTrackingOpt);
   }
}

std::string RooAbsL::GetName() const
{
   std::string output("likelihood of pdf ");
   output.append(pdf_->GetName());
   return output;
}

std::string RooAbsL::GetTitle() const
{
   std::string output("likelihood of pdf ");
   output.append(pdf_->GetTitle());
   return output;
}

double RooAbsL::defaultErrorLevel() const
{
   return 0.5;
}
std::size_t RooAbsL::numDataEntries() const
{
   return static_cast<std::size_t>(data_->numEntries());
}

void RooAbsL::optimize_pdf()
{
   // TODO: implement
}

std::size_t RooAbsL::get_N_events() const
{
   return N_events;
}

std::size_t RooAbsL::get_N_components() const
{
   return N_components;
}

double RooAbsL::get_carry() const
{
   return eval_carry_;
}

bool RooAbsL::is_extended() const
{
   return extended_;
}

void RooAbsL::set_sim_count(std::size_t value)
{
   sim_count_ = value;
}


} // namespace TestStatistics
} // namespace RooFit
