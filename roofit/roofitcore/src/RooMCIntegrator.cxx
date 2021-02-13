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
\file RooMCIntegrator.cxx
\class RooMCIntegrator
\ingroup Roofitcore

RooMCIntegrator implements an adaptive multi-dimensional Monte Carlo
numerical integration, following the VEGAS algorithm originally described
in G. P. Lepage, J. Comp. Phys. 27, 192(1978). This implementation is
based on a C version from the 0.9 beta release of the GNU scientific library.
**/

#include "RooFit.h"
#include "Riostream.h"

#include "TMath.h"
#include "TClass.h"
#include "RooMCIntegrator.h"
#include "RooArgSet.h"
#include "RooNumber.h"
#include "RooAbsArg.h"
#include "RooNumIntFactory.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooMsgService.h"

#include <math.h>



using namespace std;

ClassImp(RooMCIntegrator);
;

// Register this class with RooNumIntFactory


////////////////////////////////////////////////////////////////////////////////
/// This function registers class RooMCIntegrator, its configuration options
/// and its capabilities with RooNumIntFactory

void RooMCIntegrator::registerIntegrator(RooNumIntFactory& fact)
{
  RooCategory samplingMode("samplingMode","Sampling Mode") ;
  samplingMode.defineType("Importance",RooMCIntegrator::Importance) ;
  samplingMode.defineType("ImportanceOnly",RooMCIntegrator::ImportanceOnly) ;
  samplingMode.defineType("Stratified",RooMCIntegrator::Stratified) ;
  samplingMode.setIndex(RooMCIntegrator::Importance) ;

  RooCategory genType("genType","Generator Type") ;
  genType.defineType("QuasiRandom",RooMCIntegrator::QuasiRandom) ;
  genType.defineType("PseudoRandom",RooMCIntegrator::PseudoRandom) ;
  genType.setIndex(RooMCIntegrator::QuasiRandom) ;

  RooCategory verbose("verbose","Verbose flag") ;
  verbose.defineType("true",1) ;
  verbose.defineType("false",0) ;
  verbose.setIndex(0) ;

  RooRealVar alpha("alpha","Grid structure constant",1.5) ;
  RooRealVar nRefineIter("nRefineIter","Number of refining iterations",5) ;
  RooRealVar nRefinePerDim("nRefinePerDim","Number of refining samples (per dimension)",1000) ;
  RooRealVar nIntPerDim("nIntPerDim","Number of integration samples (per dimension)",5000) ;
  
  // Create prototype integrator
  RooMCIntegrator* proto = new RooMCIntegrator() ;

  // Register prototype and default config with factory
  fact.storeProtoIntegrator(proto,RooArgSet(samplingMode,genType,verbose,alpha,nRefineIter,nRefinePerDim,nIntPerDim)) ;

  // Make this method the default for all N>2-dim integrals
  RooNumIntConfig::defaultConfig().methodND().setLabel(proto->IsA()->GetName()) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Default constructor 
/// 
/// coverity[UNINIT_CTOR] 

 RooMCIntegrator::RooMCIntegrator()
{
}



////////////////////////////////////////////////////////////////////////////////
/// Construct an integrator over 'function' with given sampling mode
/// and generator type.  The sampling mode can be 'Importance'
/// (default), 'ImportanceOnly' and 'Stratified'. The generator type
/// can be 'QuasiRandom' (default) and 'PseudoRandom'. Consult the original
/// VEGAS documentation on details of the mode and type parameters.

RooMCIntegrator::RooMCIntegrator(const RooAbsFunc& function, SamplingMode mode,
				 GeneratorType genType, Bool_t verbose) :
  RooAbsIntegrator(function), _grid(function), _verbose(verbose),
  _alpha(1.5),  _mode(mode), _genType(genType),
  _nRefineIter(5),_nRefinePerDim(1000),_nIntegratePerDim(5000)
{
  // coverity[UNINIT_CTOR]
  if(!(_valid= _grid.isValid())) return;
  if(_verbose) _grid.Print();
} 



////////////////////////////////////////////////////////////////////////////////
/// Construct an integrator over 'function' where the configuration details
/// are taken from 'config'

RooMCIntegrator::RooMCIntegrator(const RooAbsFunc& function, const RooNumIntConfig& config) :
  RooAbsIntegrator(function), _grid(function)
{ 
  const RooArgSet& configSet = config.getConfigSection(IsA()->GetName()) ;
  _verbose = (Bool_t) configSet.getCatIndex("verbose",0) ;
  _alpha = configSet.getRealValue("alpha",1.5) ;
  _mode = (SamplingMode) configSet.getCatIndex("samplingMode",Importance) ;
  _genType = (GeneratorType) configSet.getCatIndex("genType",QuasiRandom) ;
  _nRefineIter = (Int_t) configSet.getRealValue("nRefineIter",5) ;
  _nRefinePerDim = (Int_t) configSet.getRealValue("nRefinePerDim",1000) ;
  _nIntegratePerDim = (Int_t) configSet.getRealValue("nIntPerDim",5000) ;

  // check that our grid initialized without errors
  if(!(_valid= _grid.isValid())) return;
  if(_verbose) _grid.Print();
} 



////////////////////////////////////////////////////////////////////////////////
/// Return clone of this generator operating on given function with given configuration
/// Needed to support RooNumIntFactory

RooAbsIntegrator* RooMCIntegrator::clone(const RooAbsFunc& function, const RooNumIntConfig& config) const
{
  return new RooMCIntegrator(function,config) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooMCIntegrator::~RooMCIntegrator() 
{
}



////////////////////////////////////////////////////////////////////////////////
/// Check if we can integrate over the current domain. If return value
/// is kTRUE we cannot handle the current limits (e.g. where the domain
/// of one or more observables is open ended.

Bool_t RooMCIntegrator::checkLimits() const 
{
  return _grid.initialize(*integrand());
}



////////////////////////////////////////////////////////////////////////////////
/// Evaluate the integral using a fixed number of calls to evaluate the integrand
/// equal to about 10k per dimension. Use the first 5k calls to refine the grid
/// over 5 iterations of 1k calls each, and the remaining 5k calls for a single
/// high statistics integration.

Double_t RooMCIntegrator::integral(const Double_t* /*yvec*/) 
{
  _timer.Start(kTRUE);
  vegas(AllStages,_nRefinePerDim*_grid.getDimension(),_nRefineIter);
  Double_t ret = vegas(ReuseGrid,_nIntegratePerDim*_grid.getDimension(),1);
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Perform one step of Monte Carlo integration using the specified number of iterations
/// with (approximately) the specified number of integrand evaluation calls per iteration.
/// Use the VEGAS algorithm, starting from the specified stage. Returns the best estimate
/// of the integral. Also sets *absError to the estimated absolute error of the integral
/// estimate if absError is non-zero.

Double_t RooMCIntegrator::vegas(Stage stage, UInt_t calls, UInt_t iterations, Double_t *absError) 
{
  //cout << "VEGAS stage = " << stage << " calls = " << calls << " iterations = " << iterations << endl ;

  // reset the grid to its initial state if we are starting from scratch
  if(stage == AllStages) _grid.initialize(*_function);

  // reset the results of previous calculations on this grid, but reuse the grid itself.
  if(stage <= ReuseGrid) {
    _wtd_int_sum = 0;
    _sum_wgts = 0;
    _chi_sum = 0;
    _it_num = 1;
    _samples = 0;
  }

  // refine the results of previous calculations on the current grid.
  if(stage <= RefineGrid) {
    UInt_t bins = RooGrid::maxBins;
    UInt_t boxes = 1;
    UInt_t dim(_grid.getDimension());

    // select the sampling mode for the next step
    if(_mode != ImportanceOnly) {
      // calculate the largest number of equal subdivisions ("boxes") along each
      // axis that results in an average of no more than 2 integrand calls per cell
      boxes = (UInt_t)floor(TMath::Power(calls/2.0,1.0/dim));
      // use stratified sampling if we are allowed enough calls (or equivalently,
      // if the dimension is low enough)
      _mode = Importance;
      if (2*boxes >= RooGrid::maxBins) {
	_mode = Stratified;
	// adjust the number of bins and boxes to give an integral number >= 1 of boxes per bin
	Int_t box_per_bin= (boxes > RooGrid::maxBins) ? boxes/RooGrid::maxBins : 1;
	bins= boxes/box_per_bin;
	if(bins > RooGrid::maxBins) bins= RooGrid::maxBins;
	boxes = box_per_bin * bins;	
	oocxcoutD((TObject*)0,Integration) << "RooMCIntegrator: using stratified sampling with " << bins << " bins and "
					 << box_per_bin << " boxes/bin" << endl;
      }
      else {
	oocxcoutD((TObject*)0,Integration) << "RooMCIntegrator: using importance sampling with " << bins << " bins and "
					 << boxes << " boxes" << endl;
      }
    }

    // calculate the total number of n-dim boxes for this step
    Double_t tot_boxes = TMath::Power((Double_t)boxes,(Double_t)dim);

    // increase the total number of calls to get at least 2 calls per box, if necessary
    _calls_per_box = (UInt_t)(calls/tot_boxes);
    if(_calls_per_box < 2) _calls_per_box= 2;
    calls= (UInt_t)(_calls_per_box*tot_boxes);

    // calculate the Jacobean factor: volume/(avg # of calls/bin)
    _jac = _grid.getVolume()*TMath::Power((Double_t)bins,(Double_t)dim)/calls;

    // setup our grid to use the calculated number of boxes and bins
    _grid.setNBoxes(boxes);
    if(bins != _grid.getNBins()) _grid.resize(bins);
  }

  // allocate memory for some book-keeping arrays
  UInt_t *box= _grid.createIndexVector();
  UInt_t *bin= _grid.createIndexVector();
  Double_t *x= _grid.createPoint();

  // loop over iterations for this step
  Double_t cum_int(0),cum_sig(0);
  _it_start = _it_num;
  _chisq = 0.0;
  for (UInt_t it = 0; it < iterations; it++) {
    Double_t intgrl(0),intgrl_sq(0),sig(0),jacbin(_jac);
    
    _it_num = _it_start + it;
    
    // reset the values associated with each grid cell
    _grid.resetValues();

    // loop over grid boxes
    _grid.firstBox(box);
    do {
      Double_t m(0),q(0);
      // loop over integrand evaluations within this grid box
      for(UInt_t k = 0; k < _calls_per_box; k++) {
        // generate a random point in this box
        Double_t bin_vol(0);
        _grid.generatePoint(box, x, bin, bin_vol, _genType == QuasiRandom ? kTRUE : kFALSE);
        // evaluate the integrand at the generated point
        Double_t fval= jacbin*bin_vol*integrand(x);
        // update mean and variance calculations
        Double_t d = fval - m;
        m+= d / (k + 1.0);
        q+= d * d * (k / (k + 1.0));
        // accumulate the results of this evaluation (importance sampling only)
        if (_mode != Stratified) _grid.accumulate(bin, fval*fval);
      }
      intgrl += m * _calls_per_box;
      Double_t f_sq_sum = q * _calls_per_box ;
      sig += f_sq_sum ;

      // accumulate the results for this grid box (stratified sampling only)      
      if (_mode == Stratified) _grid.accumulate(bin, f_sq_sum);

      // print occasional progress messages
      if(_timer.RealTime() > 30) {
        std::size_t index = 0;
        std::size_t sizeOfDim = 1;

        for (unsigned int i=0; i < _grid.getDimension(); ++i) {
          index += box[i] * sizeOfDim;
          sizeOfDim *= _grid.getNBoxes();
        }
        oocoutP(this,Integration) << "RooMCIntegrator: still working ... iteration "
            << it << '/' << iterations << "  box " << index << "/"<< std::pow(_grid.getNBoxes(), _grid.getDimension()) << endl;
        _timer.Start(kTRUE);
      }
      else {
        _timer.Start(kFALSE);
      }

    } while(_grid.nextBox(box));

    // compute final results for this iteration
    Double_t wgt;
    sig = sig / (_calls_per_box - 1.0)  ;    
    if (sig > 0) {
      wgt = 1.0 / sig;
    }
    else if (_sum_wgts > 0) {
      wgt = _sum_wgts / _samples;
    }
    else {
      wgt = 0.0;
    }
    intgrl_sq = intgrl * intgrl;
    _result = intgrl;
    _sigma  = sqrt(sig);
    
    if (wgt > 0.0) {
      _samples++ ;
      _sum_wgts += wgt;
      _wtd_int_sum += intgrl * wgt;
      _chi_sum += intgrl_sq * wgt;
      
      cum_int = _wtd_int_sum / _sum_wgts;
      cum_sig = sqrt (1 / _sum_wgts);
      
      if (_samples > 1) {
	_chisq = (_chi_sum - _wtd_int_sum * cum_int)/(_samples - 1.0);
      }
    }
    else {
      cum_int += (intgrl - cum_int) / (it + 1.0);
      cum_sig = 0.0;
    }         
    oocxcoutD((TObject*)0,Integration) << "=== Iteration " << _it_num << " : I = " << intgrl << " +/- " << sqrt(sig) << endl
				     << "    Cumulative : I = " << cum_int << " +/- " << cum_sig << "( chi2 = " << _chisq
				     << ")" << endl;
    // print the grid after the final iteration
    if (oodologD((TObject*)0,Integration)) {
      if(it + 1 == iterations) _grid.Print("V");
    }
    _grid.refine(_alpha);
  }

  // cleanup
  delete[] bin;
  delete[] box;
  delete[] x;

  if(absError) *absError = cum_sig;
  return cum_int;
}
