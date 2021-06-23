/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooNLLVar.cxx
\class RooNLLVar
\ingroup Roofitcore

Class RooNLLVar implements a -log(likelihood) calculation from a dataset
and a PDF. The NLL is calculated as
\f[
 \sum_\mathrm{data} -\log( \mathrm{pdf}(x_\mathrm{data}))
\f]
In extended mode, a
\f$ N_\mathrm{expect} - N_\mathrm{observed}*log(N_\mathrm{expect}) \f$ term is added.
**/

#include "RooNLLVar.h"

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooCmdConfig.h"
#include "RooMsgService.h"
#include "RooAbsDataStore.h"
#include "RooRealMPFE.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooProdPdf.h"
#include "RooNaNPacker.h"
#include "RunContext.h"

#ifdef ROOFIT_CHECK_CACHED_VALUES
#include <iomanip>
#endif

#include "TMath.h"
#include "Math/Util.h"

#include <algorithm>

namespace {
  template<class ...Args>
  RooAbsTestStatistic::Configuration makeRooAbsTestStatisticCfg(Args const& ... args) {
    RooAbsTestStatistic::Configuration cfg;
    cfg.rangeName = RooCmdConfig::decodeStringOnTheFly("RooNLLVar::RooNLLVar","RangeWithName",0,"",args...);
    cfg.addCoefRangeName = RooCmdConfig::decodeStringOnTheFly("RooNLLVar::RooNLLVar","AddCoefRange",0,"",args...);
    cfg.nCPU = RooCmdConfig::decodeIntOnTheFly("RooNLLVar::RooNLLVar","NumCPU",0,1,args...);
    cfg.interleave = RooFit::BulkPartition;
    cfg.CPUAffinity = static_cast<bool>(RooCmdConfig::decodeIntOnTheFly("RooNLLVar::RooNLLVar","CPUAffinity",0,1,args...));
    cfg.verbose = static_cast<bool>(RooCmdConfig::decodeIntOnTheFly("RooNLLVar::RooNLLVar","Verbose",0,1,args...));
    cfg.splitCutRange = static_cast<bool>(RooCmdConfig::decodeIntOnTheFly("RooNLLVar::RooNLLVar","SplitRange",0,0,args...));
    cfg.cloneInputData = static_cast<bool>(RooCmdConfig::decodeIntOnTheFly("RooNLLVar::RooNLLVar","CloneData",0,1,args...));
    cfg.integrateOverBinsPrecision = RooCmdConfig::decodeDoubleOnTheFly("RooNLLVar::RooNLLVar", "IntegrateBins", 0, -1., {args...});
    return cfg;
  }
}

ClassImp(RooNLLVar)

RooArgSet RooNLLVar::_emptySet ;

RooNLLVar::RooNLLVar()
{ }

////////////////////////////////////////////////////////////////////////////////
/// Construct likelihood from given p.d.f and (binned or unbinned dataset)
///
///  Argument                 | Description
///  -------------------------|------------
///  Extended()               | Include extended term in calculation
///  NumCPU()                 | Activate parallel processing feature
///  Range()                  | Fit only selected region
///  SumCoefRange()           | Set the range in which to interpret the coefficients of RooAddPdf components
///  SplitRange()             | Fit range is split by index category of simultaneous PDF
///  ConditionalObservables() | Define conditional observables
///  Verbose()                | Verbose output of GOF framework classes
///  CloneData()              | Clone input dataset for internal use (default is kTRUE)
///  BatchMode()              | Evaluate batches of data events (faster if PDFs support it)
///  IntegrateBins() | Integrate PDF within each bin. This sets the desired precision. Only useful for binned fits.
RooNLLVar::RooNLLVar(const char *name, const char* title, RooAbsPdf& pdf, RooAbsData& indata,
		     const RooCmdArg& arg1, const RooCmdArg& arg2,const RooCmdArg& arg3,
		     const RooCmdArg& arg4, const RooCmdArg& arg5,const RooCmdArg& arg6,
		     const RooCmdArg& arg7, const RooCmdArg& arg8,const RooCmdArg& arg9) :
  RooAbsOptTestStatistic(name,title,pdf,indata,
                         *static_cast<const RooArgSet*>(RooCmdConfig::decodeObjOnTheFly(
                             "RooNLLVar::RooNLLVar","ProjectedObservables",0,&_emptySet,
                             arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)),
                         makeRooAbsTestStatisticCfg(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9))
{
  RooCmdConfig pc("RooNLLVar::RooNLLVar") ;
  pc.allowUndefined() ;
  pc.defineInt("extended","Extended",0,kFALSE) ;
  pc.defineInt("BatchMode", "BatchMode", 0, false);

  pc.process(arg1) ;  pc.process(arg2) ;  pc.process(arg3) ;
  pc.process(arg4) ;  pc.process(arg5) ;  pc.process(arg6) ;
  pc.process(arg7) ;  pc.process(arg8) ;  pc.process(arg9) ;

  _extended = pc.getInt("extended") ;
  _batchEvaluations = pc.getInt("BatchMode");
  _weightSq = kFALSE ;
  _first = kTRUE ;
  _offsetSaveW2 = 0.;

  _binnedPdf = 0 ;
}



////////////////////////////////////////////////////////////////////////////////
/// Construct likelihood from given p.d.f and (binned or unbinned dataset)
/// For internal use.

RooNLLVar::RooNLLVar(const char *name, const char *title, RooAbsPdf& pdf, RooAbsData& indata,
                     RooAbsTestStatistic::Configuration const& cfg, bool extended) :
  RooAbsOptTestStatistic(name,title,pdf,indata,RooArgSet(),cfg),
  _extended(extended),
  _weightSq(kFALSE),
  _first(kTRUE)
{
  // If binned likelihood flag is set, pdf is a RooRealSumPdf representing a yield vector
  // for a binned likelihood calculation
  _binnedPdf = cfg.binnedL ? (RooRealSumPdf*)_funcClone : 0 ;

  // Retrieve and cache bin widths needed to convert un-normalized binnedPdf values back to yields
  if (_binnedPdf) {

    // The Active label will disable pdf integral calculations
    _binnedPdf->setAttribute("BinnedLikelihoodActive") ;

    RooArgSet* obs = _funcClone->getObservables(_dataClone) ;
    if (obs->getSize()!=1) {
      _binnedPdf = 0 ;
    } else {
      RooRealVar* var = (RooRealVar*) obs->first() ;
      std::list<Double_t>* boundaries = _binnedPdf->binBoundaries(*var,var->getMin(),var->getMax()) ;
      std::list<Double_t>::iterator biter = boundaries->begin() ;
      _binw.resize(boundaries->size()-1) ;
      Double_t lastBound = (*biter) ;
      ++biter ;
      int ibin=0 ;
      while (biter!=boundaries->end()) {
	_binw[ibin] = (*biter) - lastBound ;
	lastBound = (*biter) ;
	ibin++ ;
	++biter ;
      }
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Construct likelihood from given p.d.f and (binned or unbinned dataset)
/// For internal use.

RooNLLVar::RooNLLVar(const char *name, const char *title, RooAbsPdf& pdf, RooAbsData& indata,
                     const RooArgSet& projDeps,
                     RooAbsTestStatistic::Configuration const& cfg, bool extended) :
  RooAbsOptTestStatistic(name,title,pdf,indata,projDeps, cfg),
  _extended(extended),
  _weightSq(kFALSE),
  _first(kTRUE)
{
  // If binned likelihood flag is set, pdf is a RooRealSumPdf representing a yield vector
  // for a binned likelihood calculation
  _binnedPdf = cfg.binnedL ? (RooRealSumPdf*)_funcClone : 0 ;

  // Retrieve and cache bin widths needed to convert un-normalized binnedPdf values back to yields
  if (_binnedPdf) {

    RooArgSet* obs = _funcClone->getObservables(_dataClone) ;
    if (obs->getSize()!=1) {
      _binnedPdf = 0 ;
    } else {
      RooRealVar* var = (RooRealVar*) obs->first() ;
      std::list<Double_t>* boundaries = _binnedPdf->binBoundaries(*var,var->getMin(),var->getMax()) ;
      std::list<Double_t>::iterator biter = boundaries->begin() ;
      _binw.resize(boundaries->size()-1) ;
      Double_t lastBound = (*biter) ;
      ++biter ;
      int ibin=0 ;
      while (biter!=boundaries->end()) {
	_binw[ibin] = (*biter) - lastBound ;
	lastBound = (*biter) ;
	ibin++ ;
	++biter ;
      }
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooNLLVar::RooNLLVar(const RooNLLVar& other, const char* name) :
  RooAbsOptTestStatistic(other,name),
  _extended(other._extended),
  _batchEvaluations(other._batchEvaluations),
  _weightSq(other._weightSq),
  _first(kTRUE),
  _offsetSaveW2(other._offsetSaveW2),
  _binw(other._binw) {
  _binnedPdf = other._binnedPdf ? (RooRealSumPdf*)_funcClone : 0 ;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a test statistic using several properties of the current instance. This is used to duplicate
/// the test statistic in multi-processing scenarios.
RooAbsTestStatistic* RooNLLVar::create(const char *name, const char *title, RooAbsReal& pdf, RooAbsData& adata,
            const RooArgSet& projDeps, RooAbsTestStatistic::Configuration const& cfg) {
  RooAbsPdf & thePdf = dynamic_cast<RooAbsPdf&>(pdf);
  // check if pdf can be extended
  bool extendedPdf = _extended && thePdf.canBeExtended();

  auto testStat = new RooNLLVar(name, title, thePdf, adata, projDeps, cfg, extendedPdf);
  testStat->batchMode(_batchEvaluations);
  return testStat;
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooNLLVar::~RooNLLVar()
{
}




////////////////////////////////////////////////////////////////////////////////

void RooNLLVar::applyWeightSquared(Bool_t flag)
{
  if (_gofOpMode==Slave) {
    if (flag != _weightSq) {
      _weightSq = flag;
      std::swap(_offset, _offsetSaveW2);
    }
    setValueDirty();
  } else if ( _gofOpMode==MPMaster) {
    for (Int_t i=0 ; i<_nCPU ; i++)
      _mpfeArray[i]->applyNLLWeightSquared(flag);
  } else if ( _gofOpMode==SimMaster) {
    for (Int_t i=0 ; i<_nGof ; i++)
      ((RooNLLVar*)_gofArray[i])->applyWeightSquared(flag);
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate and return likelihood on subset of data.
/// \param[in] firstEvent First event to be processed.
/// \param[in] lastEvent  First event not to be processed, any more.
/// \param[in] stepSize   Steps between events.
/// \note For batch computations, the step size **must** be one.
///
/// If this an extended likelihood, the extended term is added to the return likelihood
/// in the batch that encounters the event with index 0.

Double_t RooNLLVar::evaluatePartition(std::size_t firstEvent, std::size_t lastEvent, std::size_t stepSize) const
{
  // Throughout the calculation, we use Kahan's algorithm for summing to
  // prevent loss of precision - this is a factor four more expensive than
  // straight addition, but since evaluating the PDF is usually much more
  // expensive than that, we tolerate the additional cost...
  ROOT::Math::KahanSum<double> result{0.0};
  double sumWeight{0.0};

  RooAbsPdf* pdfClone = (RooAbsPdf*) _funcClone ;

  // cout << "RooNLLVar::evaluatePartition(" << GetName() << ") projDeps = " << (_projDeps?*_projDeps:RooArgSet()) << endl ;

  _dataClone->store()->recalculateCache( _projDeps, firstEvent, lastEvent, stepSize, (_binnedPdf?kFALSE:kTRUE) ) ;



  // If pdf is marked as binned - do a binned likelihood calculation here (sum of log-Poisson for each bin)
  if (_binnedPdf) {
    ROOT::Math::KahanSum<double> sumWeightKahanSum{0.0};
    for (auto i=firstEvent ; i<lastEvent ; i+=stepSize) {

      _dataClone->get(i) ;

      if (!_dataClone->valid()) continue;

      Double_t eventWeight = _dataClone->weight();


      // Calculate log(Poisson(N|mu) for this bin
      Double_t N = eventWeight ;
      Double_t mu = _binnedPdf->getVal()*_binw[i] ;
      //cout << "RooNLLVar::binnedL(" << GetName() << ") N=" << N << " mu = " << mu << endl ;

      if (mu<=0 && N>0) {

        // Catch error condition: data present where zero events are predicted
        logEvalError(Form("Observed %f events in bin %lu with zero event yield",N,(unsigned long)i)) ;

      } else if (fabs(mu)<1e-10 && fabs(N)<1e-10) {

        // Special handling of this case since log(Poisson(0,0)=0 but can't be calculated with usual log-formula
        // since log(mu)=0. No update of result is required since term=0.

      } else {

        result += -1*(-mu + N*log(mu) - TMath::LnGamma(N+1));
        sumWeightKahanSum += eventWeight;

      }
    }

    sumWeight = sumWeightKahanSum.Sum();

  } else { //unbinned PDF

    if (_batchEvaluations) {
      std::tie(result, sumWeight) = computeBatched(stepSize, firstEvent, lastEvent);
#ifdef ROOFIT_CHECK_CACHED_VALUES

      ROOT::Math::KahanSum<double> resultScalar, sumWeightScalar;
      std::tie(resultScalar, sumWeightScalar) = computeScalar(stepSize, firstEvent, lastEvent);
      double carryScalar = resultScalar.Carry();

      constexpr bool alwaysPrint = false;

      if (alwaysPrint || fabs(result - resultScalar)/resultScalar > 5.E-15) {
        std::cerr << "RooNLLVar: result is off\n\t" << std::setprecision(15) << result
            << "\n\t" << resultScalar << std::endl;
      }

      if (alwaysPrint || fabs(carry - carryScalar)/carryScalar > 500.) {
        std::cerr << "RooNLLVar: carry is far off\n\t" << std::setprecision(15) << carry
            << "\n\t" << carryScalar << std::endl;
      }

      if (alwaysPrint || fabs(sumWeight - sumWeightScalar)/sumWeightScalar > 1.E-15) {
        std::cerr << "RooNLLVar: sumWeight is off\n\t" << std::setprecision(15) << sumWeight
            << "\n\t" << sumWeightScalar << std::endl;
      }

#endif
    } else { //scalar mode
      std::tie(result, sumWeight) = computeScalar(stepSize, firstEvent, lastEvent);
    }

    // include the extended maximum likelihood term, if requested
    if(_extended && _setNum==_extSet) {
      if (_weightSq) {


        // Calculate sum of weights-squared here for extended term
        Double_t sumW2;
        if (_batchEvaluations) {
          const RooSpan<const double> eventWeights = _dataClone->getWeightBatch(0, _nEvents);
          if (eventWeights.empty()) {
            sumW2 = (lastEvent - firstEvent) * _dataClone->weightSquared();
          } else {
            ROOT::Math::KahanSum<double, 4u> kahanWeight;
            for (std::size_t i = 0; i < eventWeights.size(); ++i) {
              kahanWeight.AddIndexed(eventWeights[i] * eventWeights[i], i);
            }
            sumW2 = kahanWeight.Sum();
          }
        } else { // scalar mode
          ROOT::Math::KahanSum<double> sumW2KahanSum;
          for (decltype(_dataClone->numEntries()) i = 0; i < _dataClone->numEntries() ; i++) {
            _dataClone->get(i);
            sumW2KahanSum += _dataClone->weightSquared();
          }
          sumW2 = sumW2KahanSum.Sum();
        }

        Double_t expected= pdfClone->expectedEvents(_dataClone->get());

        // Adjust calculation of extended term with W^2 weighting: adjust poisson such that
        // estimate of Nexpected stays at the same value, but has a different variance, rescale
        // both the observed and expected count of the Poisson with a factor sum[w] / sum[w^2] which is
        // the effective weight of the Poisson term.
        // i.e. change Poisson(Nobs = sum[w]| Nexp ) --> Poisson( sum[w] * sum[w] / sum[w^2] | Nexp * sum[w] / sum[w^2] )
        // weighted by the effective weight  sum[w^2]/ sum[w] in the likelihood.
        // Since here we compute the likelihood with the weight square we need to multiply by the
        // square of the effective weight
        // expectedW = expected * sum[w] / sum[w^2]   : effective expected entries
        // observedW =  sum[w]  * sum[w] / sum[w^2]   : effective observed entries
        // The extended term for the likelihood weighted by the square of the weight will be then:
        //  (sum[w^2]/ sum[w] )^2 * expectedW -  (sum[w^2]/ sum[w] )^2 * observedW * log (expectedW)  and this is
        //  using the previous expressions for expectedW and observedW
        //  sum[w^2] / sum[w] * expected - sum[w^2] * log (expectedW)
        //  and since the weights are constants in the likelihood we can use log(expected) instead of log(expectedW)

        Double_t expectedW2 = expected * sumW2 / _dataClone->sumEntries() ;
        Double_t extra= expectedW2 - sumW2*log(expected );

        // Double_t extra = pdfClone->extendedTerm(sumW2, _dataClone->get());

        result += extra;
      } else {
        result += pdfClone->extendedTerm(_dataClone->sumEntries(), _dataClone->get());
      }
    }
  } //unbinned PDF


  // If part of simultaneous PDF normalize probability over
  // number of simultaneous PDFs: -sum(log(p/n)) = -sum(log(p)) + N*log(n)
  if (_simCount>1) {
    result += sumWeight * log(1.0*_simCount);
  }


  // At the end of the first full calculation, wire the caches
  if (_first) {
    _first = kFALSE ;
    _funcClone->wireAllCaches() ;
  }


  // Check if value offset flag is set.
  if (_doOffset) {

    // If no offset is stored enable this feature now
    if (_offset==0 && result !=0 ) {
      coutI(Minimization) << "RooNLLVar::evaluatePartition(" << GetName() << ") first = "<< firstEvent << " last = " << lastEvent << " Likelihood offset now set to " << result << std::endl ;
      _offset = result ;
    }

    // Subtract offset
    result -= _offset;
  }

  _evalCarry = result.Carry();
  return result.Sum() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Compute probabilites of all data events. Use faster batch interface.
/// \param[in] stepSize Stride when moving through the dataset.
///   \note For batch computations, the step size **must** be one.
/// \param[in] firstEvent  First event to be processed.
/// \param[in] lastEvent   First event not to be processed.
/// \return Tuple with (Kahan sum of probabilities, carry of kahan sum, sum of weights)
RooNLLVar::ComputeResult RooNLLVar::computeBatched(std::size_t stepSize, std::size_t firstEvent, std::size_t lastEvent) const
{
  const auto nEvents = lastEvent - firstEvent;

  if (stepSize != 1) {
    throw std::invalid_argument(std::string("Error in ") + __FILE__ + ": Step size for batch computations can only be 1.");
  }

  auto pdfClone = static_cast<const RooAbsPdf*>(_funcClone);

  // Create a RunContext that will own the memory where computation results are stored.
  // Holding on to this struct in between function calls will make sure that the memory
  // is only allocated once.
  if (!_evalData) {
    _evalData.reset(new RooBatchCompute::RunContext);
  }
  _evalData->clear();
  _dataClone->getBatches(*_evalData, firstEvent, nEvents);

  auto results = pdfClone->getLogProbabilities(*_evalData, _normSet);

#ifdef ROOFIT_CHECK_CACHED_VALUES

  for (std::size_t evtNo = firstEvent; evtNo < std::min(lastEvent, firstEvent + 10); ++evtNo) {
    _dataClone->get(evtNo);
    if (_dataClone->weight() == 0.) // 0-weight events are not cached, so cannot compare against them.
      continue;

    assert(_dataClone->valid());
    try {
      // Cross check results with strict tolerance and complain
      BatchInterfaceAccessor::checkBatchComputation(*pdfClone, *_evalData, evtNo-firstEvent, _normSet, 1.E-13);
    } catch (std::exception& e) {
      std::cerr << __FILE__ << ":" << __LINE__ << " ERROR when checking batch computation for event " << evtNo << ":\n"
          << e.what() << std::endl;

      // It becomes a real problem if it's very wrong. We fail in this case:
      try {
        BatchInterfaceAccessor::checkBatchComputation(*pdfClone, *_evalData, evtNo-firstEvent, _normSet, 1.E-9);
      } catch (std::exception& e2) {
        assert(false);
      }
    }
  }

#endif


  // Compute sum of event weights. First check if we need squared weights
  const RooSpan<const double> eventWeights = _dataClone->getWeightBatch(firstEvent, nEvents);
  //Capture member for lambda:
  const bool retrieveSquaredWeights = _weightSq;
  auto retrieveWeight = [&eventWeights, retrieveSquaredWeights](std::size_t i) {
    return retrieveSquaredWeights ? eventWeights[i] * eventWeights[i] : eventWeights[i];
  };

  //Sum the event weights and probabilities
  ROOT::Math::KahanSum<double, 4u> kahanProb;
  double uniformSingleEventWeight{0.0};
  double sumOfWeights;
  if (eventWeights.empty()) {
    uniformSingleEventWeight = retrieveSquaredWeights ? _dataClone->weightSquared() : _dataClone->weight();
    sumOfWeights = nEvents * uniformSingleEventWeight;
    for (std::size_t i = 0; i < results.size(); ++i) { //CHECK_VECTORISE
      kahanProb.AddIndexed(-uniformSingleEventWeight * results[i], i);
    }
  } else {
    assert(results.size() == eventWeights.size());
    ROOT::Math::KahanSum<double, 4u> kahanWeight;
    for (std::size_t i = 0; i < results.size(); ++i) { //CHECK_VECTORISE
      const double weight = retrieveWeight(i);
      kahanProb.AddIndexed(-weight * results[i], i);
      kahanWeight.AddIndexed(weight, i);
    }
    sumOfWeights = kahanWeight.Sum();
  }

  if (std::isnan(kahanProb.Sum())) {
    // Special handling of evaluation errors.
    // We can recover if the bin/event that results in NaN has a weight of zero:
    ROOT::Math::KahanSum<double, 4u> kahanSanitised;
    RooNaNPacker nanPacker;
    for (std::size_t i = 0; i < results.size(); ++i) {
      double weight = eventWeights.empty() ? uniformSingleEventWeight : retrieveWeight(i);

      if (weight == 0.)
        continue;

      if (std::isnan(results[i])) {
        nanPacker.accumulate(results[i]);
      } else {
        kahanSanitised += -weight * results[i];
      }
    }

    // Some events with evaluation errors. Return "badness" of errors.
    if (nanPacker.getPayload() > 0.) {
      return {{nanPacker.getNaNWithPayload()}, sumOfWeights};
    } else {
      return {kahanSanitised, sumOfWeights};
    }
  }

  return {kahanProb, sumOfWeights};
}


RooNLLVar::ComputeResult RooNLLVar::computeScalar(std::size_t stepSize, std::size_t firstEvent, std::size_t lastEvent) const {
  auto pdfClone = static_cast<const RooAbsPdf*>(_funcClone);

  ROOT::Math::KahanSum<double> kahanWeight;
  ROOT::Math::KahanSum<double> kahanProb;
  RooNaNPacker packedNaN(0.f);

  for (auto i=firstEvent; i<lastEvent; i+=stepSize) {
    _dataClone->get(i) ;

    if (!_dataClone->valid()) continue;

    Double_t eventWeight = _dataClone->weight(); //FIXME
    if (0. == eventWeight * eventWeight) continue ;
    if (_weightSq) eventWeight = _dataClone->weightSquared() ;

    const double term = -eventWeight * pdfClone->getLogVal(_normSet);

    kahanWeight.Add(eventWeight);
    kahanProb.Add(term);
    packedNaN.accumulate(term);
  }

  if (packedNaN.getPayload() != 0.) {
    // Some events with evaluation errors. Return "badness" of errors.
    return {{packedNaN.getNaNWithPayload()}, kahanWeight.Sum()};
  }

  return {kahanProb, kahanWeight.Sum()};
}
