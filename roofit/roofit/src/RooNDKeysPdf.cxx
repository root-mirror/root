/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooNDKeysPdf.cxx 31258 2009-11-17 22:41:06Z wouter $
 * Authors:                                                                  *
 *   Max Baak, CERN, mbaak@cern.ch *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/** \class RooNDKeysPdf
    \ingroup Roofit

Generic N-dimensional implementation of a kernel estimation p.d.f. This p.d.f. models the distribution
of an arbitrary input dataset as a superposition of Gaussian kernels, one for each data point,
each contributing 1/N to the total integral of the p.d.f.
If the 'adaptive mode' is enabled, the width of the Gaussian is adaptively calculated from the
local density of events, i.e. narrow for regions with high event density to preserve details and
wide for regions with log event density to promote smoothness. The details of the general algorithm
are described in the following paper:
Cranmer KS, Kernel Estimation in High-Energy Physics.
            Computer Physics Communications 136:198-207,2001 - e-Print Archive: hep ex/0011057
For multi-dimensional datasets, the kernels are modeled by multidimensional Gaussians. The kernels are
constructed such that they reflect the correlation coefficients between the observables
in the input dataset.
**/

#include <iostream>
#include <algorithm>
#include <string>

#include "TMath.h"
#include "TMatrixDSymEigen.h"
#include "RooNDKeysPdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooHist.h"
#include "RooMsgService.h"
#include "RooChangeTracker.h"

#include "TError.h"

using namespace std;

ClassImp(RooNDKeysPdf)

////////////////////////////////////////////////////////////////////////////////
/// Construct N-dimensional kernel estimation p.d.f. in observables 'varList'
/// from dataset 'data'. Options can be
///
///  - 'a' = Use adaptive kernels (width varies with local event density)
///  - 'm' = Mirror data points over observable boundaries. Improves modeling
///         behavior at edges for distributions that are not close to zero
///         at edge
///  - 'd' = Debug flag
///  - 'v' = Verbose flag
///
/// The parameter rho (default = 1) provides an overall scale factor that can
/// be applied to the bandwith calculated for each kernel. The nSigma parameter
/// determines the size of the box that is used to search for contributing kernels
/// around a given point in observable space. The nSigma parameters is used
/// in case of non-adaptive bandwidths and for the 1st non-adaptive pass for
/// the calculation of adaptive keys p.d.f.s.
///
/// The optional weight arguments allows to specify an observable or function
/// expression in observables that specifies the weight of each event.

RooNDKeysPdf::RooNDKeysPdf(const char *name, const char *title,
			   const RooArgList& varList, const RooAbsData& data,
			   TString options, Double_t rho, Double_t nSigma, Bool_t rotate, Bool_t sortInput) : 
  RooAbsPdf(name,title),
  _varList("varList","List of variables",this),
  _rhoList("rhoList","List of rho parameters",this),
  _dataP(0),
  _data(data),
  _options(options),
  _widthFactor(rho),
  _nSigma(nSigma),
  _weights(&_weights0),
  _rotate(rotate),
  _sortInput(sortInput),
  _nAdpt(1),
  _tracker(0)
{
  // Constructor
  _varItr    = _varList.createIterator() ;
  _rhoItr    = _rhoList.createIterator() ;

  TIterator* varItr = varList.createIterator() ;
  RooAbsArg* var ;
  for (Int_t i=0; (var = (RooAbsArg*)varItr->Next()); ++i) {
    if (!dynamic_cast<RooAbsReal*>(var)) {
      coutE(InputArguments) << "RooNDKeysPdf::ctor(" << GetName() << ") ERROR: variable " << var->GetName()
             << " is not of type RooAbsReal" << endl ;
      R__ASSERT(0) ;
    }
    _varList.add(*var) ;
    _varName.push_back(var->GetName());
  }
  delete varItr ;

  createPdf();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooNDKeysPdf::RooNDKeysPdf(const char *name, const char *title,
			   const RooArgList& varList, const TH1& hist,
			   TString options, Double_t rho, Double_t nSigma, Bool_t rotate, Bool_t sortInput) : 
  RooAbsPdf(name,title),
  _varList("varList","List of variables",this),
  _rhoList("rhoList","List of rho parameters",this),
  _dataP( createDatasetFromHist(varList,hist) ), 
  _data( *_dataP ),
  _options(options),
  _widthFactor(rho),
  _nSigma(nSigma),
  _weights(&_weights0),
  _rotate(rotate),
  _sortInput(sortInput),
  _nAdpt(1),
  _tracker(0)
{
  _varItr    = _varList.createIterator() ;
  _rhoItr    = _rhoList.createIterator() ;

  TIterator* varItr = varList.createIterator() ;
  RooAbsArg* var ;
  for (Int_t i=0; (var = (RooAbsArg*)varItr->Next()); ++i) {
    if (!dynamic_cast<RooAbsReal*>(var)) {
      coutE(InputArguments) << "RooNDKeysPdf::ctor(" << GetName() << ") ERROR: variable " << var->GetName() 
			    << " is not of type RooAbsReal" << endl ;
      assert(0) ;
    }
    _varList.add(*var) ;
    _varName.push_back(var->GetName());
  }
  delete varItr ;

  createPdf();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooNDKeysPdf::RooNDKeysPdf(const char *name, const char *title,
			   const RooArgList& varList, const RooAbsData& data, const TVectorD& rho,
			   TString options, Double_t nSigma, Bool_t rotate, Bool_t sortInput) : 
  RooAbsPdf(name,title),
  _varList("varList","List of variables",this),
  _rhoList("rhoList","List of rho parameters",this),
  _dataP(0),
  _data(data),
  _options(options),
  _widthFactor(-1.0),
  _nSigma(nSigma),
  _weights(&_weights0),
  _rotate(rotate),
  _sortInput(sortInput),
  _nAdpt(1),
  _tracker(0)
{
  _varItr    = _varList.createIterator() ;
  _rhoItr    = _rhoList.createIterator() ;

  TIterator* varItr = varList.createIterator() ;
  RooAbsArg* var ;
  for (Int_t i=0; (var = (RooAbsArg*)varItr->Next()); ++i) {
    if (!dynamic_cast<RooAbsReal*>(var)) {
      coutE(InputArguments) << "RooNDKeysPdf::ctor(" << GetName() << ") ERROR: variable " << var->GetName() 
			    << " is not of type RooAbsReal" << endl ;
      R__ASSERT(0) ;
    }
    _varList.add(*var) ;
    _varName.push_back(var->GetName());
  }
  delete varItr ;

  // copy rho widths
  if( _varList.getSize() != rho.GetNrows() ) {
    coutE(InputArguments) << "ERROR:  RooNDKeysPdf::RooNDKeysPdf() : The vector-size of rho is different from that of varList."
			  << "Unable to create the PDF." << endl;
    R__ASSERT ( _varList.getSize()==rho.GetNrows() );
  }

  // negative width factor will serve as a switch in initialize()
  // negative value means that a vector has been provided as input,
  // and that _rho has already been set ...
  _rho.resize( rho.GetNrows() );
  for  (Int_t j=0; j<rho.GetNrows(); j++) { _rho[j]=rho[j]; /*cout<<"RooNDKeysPdf ctor, _rho["<<j<<"]="<<_rho[j]<<endl;*/ }

  createPdf(); // calls initialize ...
}

////////////////////////////////////////////////////////////////////////////////
/// Backward compatibility constructor for (1-dim) RooKeysPdf. If you are a new user,
/// please use the first constructor form.

RooNDKeysPdf::RooNDKeysPdf(const char *name, const char *title,
			   const RooArgList& varList, const RooAbsData& data, const RooArgList& rhoList,
			   TString options, Double_t nSigma, Bool_t rotate, Bool_t sortInput) : 
  RooAbsPdf(name,title),
  _varList("varList","List of variables",this),
  _rhoList("rhoList","List of rho parameters",this),
  _dataP(0),
  _data(data),
  _options(options),
  _widthFactor(-1.0),
  _nSigma(nSigma),
  _weights(&_weights0),
  _rotate(rotate),
  _sortInput(sortInput),
  _nAdpt(1)
{
  _varItr    = _varList.createIterator() ;
  _rhoItr = _rhoList.createIterator() ;

  TIterator* varItr = varList.createIterator() ;
  RooAbsArg* var ;
  for (Int_t i=0; (var = (RooAbsArg*)varItr->Next()); ++i) {
    if (!dynamic_cast<RooAbsReal*>(var)) {
      coutE(InputArguments) << "RooNDKeysPdf::ctor(" << GetName() << ") ERROR: variable " << var->GetName() 
			    << " is not of type RooAbsReal" << endl ;
      assert(0) ;
    }
    _varList.add(*var) ;
    _varName.push_back(var->GetName());
  }
  delete varItr ;

  TIterator* rhoItr = rhoList.createIterator() ;
  RooAbsArg* rho ;
  _rho.resize( rhoList.getSize(), 1.0 );

  for (Int_t i=0; (rho = (RooAbsArg*)rhoItr->Next()); ++i) {
    if (!dynamic_cast<RooAbsReal*>(rho)) {
      coutE(InputArguments) << "RooNDKeysPdf::ctor(" << GetName() << ") ERROR: parameter " << rho->GetName() 
			    << " is not of type RooRealVar" << endl ;
      assert(0) ;
    }
    _rhoList.add(*rho) ;
    _rho[i] = (dynamic_cast<RooAbsReal*>(rho))->getVal();
  }
  delete rhoItr ;

  // copy rho widths
  if( (_varList.getSize() != _rhoList.getSize()) ) {
    coutE(InputArguments) << "ERROR:  RooNDKeysPdf::RooNDKeysPdf() : The size of rhoList is different from varList."
			  << "Unable to create the PDF." << endl;
    assert ( _varList.getSize()==_rhoList.getSize() );
  }

  // keep track of changes in rho parameters
  _tracker = new RooChangeTracker("tracker","track rho parameters", _rhoList, true); // check for value updates
  (void) _tracker->hasChanged(true); // first evaluation always true for new parameters (?)

  createPdf();
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooNDKeysPdf::RooNDKeysPdf(const char *name, const char *title,
			   const RooArgList& varList, const TH1& hist, const RooArgList& rhoList,
			   TString options, Double_t nSigma, Bool_t rotate, Bool_t sortInput) : 
  RooAbsPdf(name,title),
  _varList("varList","List of variables",this),
  _rhoList("rhoList","List of rho parameters",this),
  _dataP( createDatasetFromHist(varList,hist) ), 
  _data( *_dataP ),
  _options(options),
  _widthFactor(-1),
  _nSigma(nSigma),
  _weights(&_weights0),
  _rotate(rotate),
  _sortInput(sortInput),
  _nAdpt(1)
{
  _varItr    = _varList.createIterator() ;
  _rhoItr    = _rhoList.createIterator() ;

  TIterator* varItr = varList.createIterator() ;
  RooAbsArg* var ;
  for (Int_t i=0; (var = (RooAbsArg*)varItr->Next()); ++i) {
    if (!dynamic_cast<RooAbsReal*>(var)) {
      coutE(InputArguments) << "RooNDKeysPdf::ctor(" << GetName() << ") ERROR: variable " << var->GetName() 
			    << " is not of type RooAbsReal" << endl ;
      assert(0) ;
    }
    _varList.add(*var) ;
    _varName.push_back(var->GetName());
  }
  delete varItr ;

  // copy rho widths
  TIterator* rhoItr = rhoList.createIterator() ;
  RooAbsArg* rho ;
  _rho.resize( rhoList.getSize(), 1.0 );

  for (Int_t i=0; (rho = (RooAbsArg*)rhoItr->Next()); ++i) {
    if (!dynamic_cast<RooAbsReal*>(rho)) {
      coutE(InputArguments) << "RooNDKeysPdf::ctor(" << GetName() << ") ERROR: parameter " << rho->GetName() 
			    << " is not of type RooRealVar" << endl ;
      assert(0) ;
    }
    _rhoList.add(*rho) ;
    _rho[i] = (dynamic_cast<RooAbsReal*>(rho))->getVal();
  }
  delete rhoItr ;

  if( (_varList.getSize() != _rhoList.getSize()) ) {
    coutE(InputArguments) << "ERROR:  RooNDKeysPdf::RooNDKeysPdf() : The size of rhoList is different from varList."
			  << "Unable to create the PDF." << endl;
    assert ( _varList.getSize()==_rhoList.getSize() );
  }

  // keep track of changes in rho parameters
  _tracker = new RooChangeTracker("tracker","track rho parameters", _rhoList, true); // check for value updates
  (void) _tracker->hasChanged(true); // first evaluation always true for new parameters (?)

  createPdf();
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooNDKeysPdf::RooNDKeysPdf(const char *name, const char *title,
                           RooAbsReal& x, const RooAbsData& data,
                           Mirror mirror, Double_t rho, Double_t nSigma, Bool_t rotate, Bool_t sortInput) : 
  RooAbsPdf(name,title),
  _varList("varList","List of variables",this),
  _rhoList("rhoList","List of rho parameters",this),
  _dataP(0),
  _data(data),
  _options("a"),
  _widthFactor(rho),
  _nSigma(nSigma),
  _weights(&_weights0),
  _rotate(rotate),
  _sortInput(sortInput),
  _nAdpt(1),
  _tracker(0)
{ 
  _varItr = _varList.createIterator() ;
  _rhoItr = _rhoList.createIterator() ;
  
  _varList.add(x) ;
  _varName.push_back(x.GetName());

  if (mirror!=NoMirror) {
    if (mirror!=MirrorBoth)
      coutW(InputArguments) << "RooNDKeysPdf::RooNDKeysPdf() : Warning : asymmetric mirror(s) no longer supported." << endl;
    _options="m";
  }

  createPdf();
}

////////////////////////////////////////////////////////////////////////////////
/// Backward compatibility constructor for Roo2DKeysPdf. If you are a new user,
/// please use the first constructor form.

RooNDKeysPdf::RooNDKeysPdf(const char *name, const char *title, RooAbsReal& x, RooAbsReal & y,
                           const RooAbsData& data, TString options, Double_t rho, Double_t nSigma, Bool_t rotate, Bool_t sortInput) : 
  RooAbsPdf(name,title),
  _varList("varList","List of variables",this),
  _rhoList("rhoList","List of rho parameters",this),
  _dataP(0),
  _data(data),
  _options(options),
  _widthFactor(rho),
  _nSigma(nSigma),
  _weights(&_weights0),
  _rotate(rotate),
  _sortInput(sortInput),
  _nAdpt(1),
  _tracker(0)
{ 
  _varItr = _varList.createIterator() ;
  _rhoItr = _rhoList.createIterator() ;

  _varList.add(RooArgSet(x,y)) ;
  _varName.push_back(x.GetName());
  _varName.push_back(y.GetName());

  createPdf();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

RooNDKeysPdf::RooNDKeysPdf(const RooNDKeysPdf& other, const char* name) :
  RooAbsPdf(other,name),
  _varList("varList",this,other._varList),
  _rhoList("rhoList",this,other._rhoList),
  _dataP( other._dataP!=NULL ? new RooDataSet(*other._dataP) : NULL ),
  _data(  other._dataP!=NULL ? *_dataP : other._data ), 
  _options(other._options),
  _widthFactor(other._widthFactor),
  _nSigma(other._nSigma),
  _weights(&_weights0),
  _rotate(other._rotate),
  _sortInput(other._sortInput),
  _nAdpt(other._nAdpt)
{
  _tracker = ( other._tracker!=NULL ? new RooChangeTracker(*other._tracker) : NULL );
  //if (_tracker!=NULL) { _tracker->hasChanged(true); }

  _varItr      = _varList.createIterator() ;
  _rhoItr      = _rhoList.createIterator() ;

  _fixedShape  = other._fixedShape;
  _mirror      = other._mirror;
  _debug       = other._debug;
  _verbose     = other._verbose;
  _sqrt2pi     = other._sqrt2pi;
  _nDim        = other._nDim;
  _nEvents     = other._nEvents;
  _nEventsM    = other._nEventsM;
  _nEventsW    = other._nEventsW;
  _d           = other._d;
  _n           = other._n;
  _dataPts     = other._dataPts;
  _dataPtsR    = other._dataPtsR;
  _weights0    = other._weights0;
  _weights1    = other._weights1;
  if (_options.Contains("a")) { _weights = &_weights1; }
  //_sortIdcs    = other._sortIdcs;
  _sortTVIdcs  = other._sortTVIdcs;
  _varName     = other._varName;
  _rho         = other._rho;
  _x           = other._x;
  _x0          = other._x0 ;
  _x1          = other._x1 ;
  _x2          = other._x2 ;
  _xDatLo      = other._xDatLo;
  _xDatHi      = other._xDatHi;
  _xDatLo3s    = other._xDatLo3s;
  _xDatHi3s    = other._xDatHi3s;
  _mean        = other._mean;
  _sigma       = other._sigma;

  // BoxInfo
  _netFluxZ    = other._netFluxZ;
  _nEventsBW   = other._nEventsBW;
  _nEventsBMSW = other._nEventsBMSW;
  _xVarLo      = other._xVarLo;
  _xVarHi      = other._xVarHi;
  _xVarLoM3s   = other._xVarLoM3s;
  _xVarLoP3s   = other._xVarLoP3s;
  _xVarHiM3s   = other._xVarHiM3s;
  _xVarHiP3s   = other._xVarHiP3s;
  _bpsIdcs     = other._bpsIdcs;
  _ibNoSort    = other._ibNoSort;
  _sIdcs       = other._sIdcs;
  _bIdcs       = other._bIdcs;
  _bmsIdcs     = other._bmsIdcs;

  _rangeBoxInfo= other._rangeBoxInfo ;
  _fullBoxInfo = other._fullBoxInfo ;

  _idx         = other._idx;
  _minWeight   = other._minWeight;
  _maxWeight   = other._maxWeight;
  _wMap        = other._wMap;

  _covMat      = new TMatrixDSym(*other._covMat);
  _corrMat     = new TMatrixDSym(*other._corrMat);
  _rotMat      = new TMatrixD(*other._rotMat);
  _sigmaR      = new TVectorD(*other._sigmaR);
  _dx          = new TVectorD(*other._dx);
  _sigmaAvgR   = other._sigmaAvgR;
}

////////////////////////////////////////////////////////////////////////////////

RooNDKeysPdf::~RooNDKeysPdf()
{
  if (_varItr)    delete _varItr;
  if (_rhoItr)    delete _rhoItr;
  if (_covMat)    delete _covMat;
  if (_corrMat)   delete _corrMat;
  if (_rotMat)    delete _rotMat;
  if (_sigmaR)    delete _sigmaR;
  if (_dx)        delete _dx;
  if (_dataP)     delete _dataP;
  if (_tracker)   delete _tracker;

  // delete all the boxinfos map
  while ( !_rangeBoxInfo.empty() ) {
    map<pair<string,int>,BoxInfo*>::iterator iter = _rangeBoxInfo.begin();
    BoxInfo* box= (*iter).second;
    _rangeBoxInfo.erase(iter);
    delete box;
  }

  _dataPts.clear();
  _dataPtsR.clear();
  _weights0.clear();
  _weights1.clear();
  //_sortIdcs.clear();
  _sortTVIdcs.clear();
}

////////////////////////////////////////////////////////////////////////////////
/// evaluation order of constructor.

void RooNDKeysPdf::createPdf(Bool_t firstCall) const
{
  if (firstCall) {
    // set options
    setOptions();
    // initialization
    initialize();
  }


  // copy dataset, calculate sigma_i's, determine min and max event weight
  loadDataSet(firstCall);

  // mirror dataset around dataset boundaries -- does not depend on event weights
  if (_mirror) mirrorDataSet();

  // store indices and weights of events with high enough weights
  loadWeightSet();

  // store indices of events in variable boundaries and box shell.
//calculateShell(&_fullBoxInfo);
  // calculate normalization needed in analyticalIntegral()
//calculatePreNorm(&_fullBoxInfo);

  // lookup table for determining which events contribute to a certain coordinate
  sortDataIndices();

  // determine static and/or adaptive bandwidth
  calculateBandWidth();
}

////////////////////////////////////////////////////////////////////////////////
/// set the configuration

void RooNDKeysPdf::setOptions() const
{
  _options.ToLower();

  if( _options.Contains("a") ) { _weights = &_weights1; }
  else                         { _weights = &_weights0; }
  if( _options.Contains("m") )   _mirror = true;
  else                           _mirror = false;
  if( _options.Contains("d") )   _debug  = true;
  else                           _debug  = false;
  if( _options.Contains("v") ) { _debug = true;  _verbose = true; }
  else                         { _debug = false; _verbose = false; }

  cxcoutD(InputArguments) << "RooNDKeysPdf::setOptions()    options = " << _options
         << "\n\tbandWidthType    = " << _options.Contains("a")
         << "\n\tmirror           = " << _mirror
         << "\n\tdebug            = " << _debug
         << "\n\tverbose          = " << _verbose
         << endl;

  if (_nSigma<2.0) {
    coutW(InputArguments) << "RooNDKeysPdf::setOptions() : Warning : nSigma = " << _nSigma << " < 2.0. "
           << "Calculated normalization could be too large."
           << endl;
  }

  // number of adaptive width iterations. Default is 1.
  if( _options.Contains("a") ) {
    if ( !sscanf( _options.Data(),"%d%*s", &_nAdpt ) ) { _nAdpt=1; }
  }
}

////////////////////////////////////////////////////////////////////////////////
/// initialization

void RooNDKeysPdf::initialize() const
{
  _sqrt2pi   = sqrt(2.0*TMath::Pi()) ;
  _nDim      = _varList.getSize();
  _nEvents   = (Int_t)_data.numEntries();
  _nEventsM  = _nEvents;
  _fixedShape= kFALSE;

  _netFluxZ = kFALSE;
  _nEventsBW = 0;
  _nEventsBMSW = 0;

  if(_nDim==0) {
    coutE(InputArguments) << "ERROR:  RooNDKeysPdf::initialize() : The observable list is empty. "
           << "Unable to begin generating the PDF." << endl;
    R__ASSERT (_nDim!=0);
  }

  if(_nEvents==0) {
    coutE(InputArguments) << "ERROR:  RooNDKeysPdf::initialize() : The input data set is empty. "
           << "Unable to begin generating the PDF." << endl;
    R__ASSERT (_nEvents!=0);
  }

  _d         = static_cast<Double_t>(_nDim);

  vector<Double_t> dummy(_nDim,0.);
  _dataPts.resize(_nEvents,dummy);
  _weights0.resize(_nEvents,dummy);
  //_sortIdcs.resize(_nDim);
  _sortTVIdcs.resize(_nDim);

  //rdh _rho.resize(_nDim,_widthFactor);

  if (_widthFactor>0) { _rho.resize(_nDim,_widthFactor); }
  // else: _rho has been provided as external input

  _x.resize(_nDim,0.);
  _x0.resize(_nDim,0.);
  _x1.resize(_nDim,0.);
  _x2.resize(_nDim,0.);

  _mean.resize(_nDim,0.);
  _sigma.resize(_nDim,0.);

  _xDatLo.resize(_nDim,0.);
  _xDatHi.resize(_nDim,0.);
  _xDatLo3s.resize(_nDim,0.);
  _xDatHi3s.resize(_nDim,0.);

  boxInfoInit(&_fullBoxInfo,0,0xFFFF);

  _minWeight=0;
  _maxWeight=0;
  _wMap.clear();

  _covMat = 0;
  _corrMat= 0;
  _rotMat = 0;
  _sigmaR = 0;
  _dx = new TVectorD(_nDim); _dx->Zero();
  _dataPtsR.resize(_nEvents,*_dx);

  _varItr->Reset() ;
  RooRealVar* var ;
  for(Int_t j=0; (var=(RooRealVar*)_varItr->Next()); ++j) {
    _xDatLo[j] = var->getMin();
    _xDatHi[j] = var->getMax();
  }
}

////////////////////////////////////////////////////////////////////////////////
/// copy the dataset and calculate some useful variables

void RooNDKeysPdf::loadDataSet(Bool_t firstCall) const
{
  // first some initialization
  _nEventsW=0.;

  TMatrixD mat(_nDim,_nDim);
  if (!_covMat)  _covMat = new TMatrixDSym(_nDim);
  if (!_corrMat) _corrMat= new TMatrixDSym(_nDim);
  if (!_rotMat)  _rotMat = new TMatrixD(_nDim,_nDim);
  if (!_sigmaR)  _sigmaR = new TVectorD(_nDim);

  mat.Zero();
  _covMat->Zero();
  _corrMat->Zero();
  _rotMat->Zero();
  _sigmaR->Zero();

  const RooArgSet* values= _data.get();
  vector<RooRealVar*> dVars(_nDim);
  for  (Int_t j=0; j<_nDim; j++) {
    dVars[j] = (RooRealVar*)values->find(_varName[j].c_str());
    _x0[j]=_x1[j]=_x2[j]=0.;
  }

  _idx.clear();
  for (Int_t i=0; i<_nEvents; i++) {

    _data.get(i); // fills dVars
    _idx.push_back(i);
    vector<Double_t>& point  = _dataPts[i];
    TVectorD& pointV = _dataPtsR[i];

    Double_t myweight = _data.weight(); // default is one?
    if ( TMath::Abs(myweight)>_maxWeight ) { _maxWeight = TMath::Abs(myweight); }
    _nEventsW += myweight;

    for (Int_t j=0; j<_nDim; j++) {
      for (Int_t k=0; k<_nDim; k++) {
   mat(j,k) += dVars[j]->getVal() * dVars[k]->getVal() * myweight;
      }
      // only need to do once
      if (firstCall)
   point[j] = pointV[j] = dVars[j]->getVal();

      _x0[j] += 1. * myweight;
      _x1[j] += point[j] * myweight ;
      _x2[j] += point[j] * point[j] * myweight ;
      if (_x2[j]!=_x2[j]) exit(3);

      // only need to do once
      if (firstCall) {
   if (point[j]<_xDatLo[j]) { _xDatLo[j]=point[j]; }
   if (point[j]>_xDatHi[j]) { _xDatHi[j]=point[j]; }
      }
    }
  }

  _n = TMath::Power(4./(_nEventsW*(_d+2.)), 1./(_d+4.)) ;
  // = (4/[n(dim(R) + 2)])^1/(dim(R)+4); dim(R) = 2
  _minWeight = (0.5 - TMath::Erf(_nSigma/sqrt(2.))/2.) * _maxWeight;

  for (Int_t j=0; j<_nDim; j++) {
    _mean[j]  = _x1[j]/_x0[j];
    _sigma[j] = sqrt(_x2[j]/_x0[j]-_mean[j]*_mean[j]);
  }

  for (Int_t j=0; j<_nDim; j++) {
    for (Int_t k=0; k<_nDim; k++) {
      (*_covMat)(j,k) = mat(j,k)/_x0[j] - _mean[j]*_mean[k];
    }
  }

  for (Int_t j=0; j<_nDim; j++) {
    for (Int_t k=0; k<_nDim; k++)
      (*_corrMat)(j,k) = (*_covMat)(j,k)/(_sigma[j]*_sigma[k]) ;
  }

  // use raw sigmas (without rho) for sigmaAvgR
  TMatrixDSymEigen evCalculator(*_covMat);
  TVectorD sigmaRraw = evCalculator.GetEigenValues();
  for (Int_t j=0; j<_nDim; j++) { sigmaRraw[j] = sqrt(sigmaRraw[j]); }

  _sigmaAvgR=1.;
  for (Int_t j=0; j<_nDim; j++) { _sigmaAvgR *= sigmaRraw[j]; }
  _sigmaAvgR = TMath::Power(_sigmaAvgR, 1./_d) ;

  // find decorrelation matrix and eigenvalues (R)
  if (_nDim>1 && _rotate) {
    // new: rotation matrix now independent of rho evaluation
    *_rotMat = evCalculator.GetEigenVectors();
    *_rotMat = _rotMat->T(); // transpose
  } else {
    TMatrixD haar(_nDim,_nDim);
    TMatrixD unit(TMatrixD::kUnit,haar);
    *_rotMat = unit;
  }

  // update sigmas (rho dependent)
  updateRho();

  //// rho no longer used after this.
  //// Now set rho = 1 because sigmaR now contains rho
  //for (Int_t j=0; j<_nDim; j++) { _rho[j] = 1.; }  // reset: important!

  if (_verbose) {
    //_covMat->Print();
    _rotMat->Print();
    _corrMat->Print();
    _sigmaR->Print();
  }

  if (_nDim>1 && _rotate) {
    // apply rotation
    for (Int_t i=0; i<_nEvents; i++) {
      TVectorD& pointR = _dataPtsR[i];
      pointR *= *_rotMat;
    }
  }

  coutI(Contents) << "RooNDKeysPdf::loadDataSet(" << this << ")" 
		  << "\n Number of events in dataset: " << _nEvents 
		  << "\n Weighted number of events in dataset: " << _nEventsW 
		  << endl; 
}

////////////////////////////////////////////////////////////////////////////////
/// determine mirror dataset.
/// mirror points are added around the physical boundaries of the dataset
/// Two steps:
/// 1. For each entry, determine if it should be mirrored (the mirror configuration).
/// 2. For each mirror configuration, make the mirror points.

void RooNDKeysPdf::mirrorDataSet() const
{
  for (Int_t j=0; j<_nDim; j++) {
    _xDatLo3s[j] = _xDatLo[j] + _nSigma * (_n * _sigma[j]);
    _xDatHi3s[j] = _xDatHi[j] - _nSigma * (_n * _sigma[j]);

    //cout<<"xDatLo3s["<<j<<"]="<<_xDatLo3s[j]<<endl;
    //cout<<"xDatHi3s["<<j<<"]="<<_xDatHi3s[j]<<endl;
  }

  vector<Double_t> dummy(_nDim,0.);

  // 1.
  for (Int_t i=0; i<_nEvents; i++) {
    vector<Double_t>& x = _dataPts[i];

    Int_t size = 1;
    vector<vector<Double_t> > mpoints(size,dummy);
    vector<vector<Int_t> > mjdcs(size);

    // create all mirror configurations for event i
    for (Int_t j=0; j<_nDim; j++) {

      vector<Int_t>& mjdxK = mjdcs[0];
      vector<Double_t>& mpointK = mpoints[0];

      // single mirror *at physical boundaries*
      if ((x[j]>_xDatLo[j] && x[j]<_xDatLo3s[j]) && x[j]<(_xDatLo[j]+_xDatHi[j])/2.) {
   mpointK[j] = 2.*_xDatLo[j]-x[j];
   mjdxK.push_back(j);
      } else if ((x[j]<_xDatHi[j] && x[j]>_xDatHi3s[j]) && x[j]>(_xDatLo[j]+_xDatHi[j])/2.) {
   mpointK[j] = 2.*_xDatHi[j]-x[j];
   mjdxK.push_back(j);
      }
    }

    vector<Int_t>& mjdx0 = mjdcs[0];
    // no mirror point(s) for this event
    if (size==1 && mjdx0.size()==0) continue;

    // 2.
    // generate all mirror points for event i
    vector<Int_t>& mjdx = mjdcs[0];
    vector<Double_t>& mpoint = mpoints[0];

    // number of mirror points for this mirror configuration
    Int_t eMir = 1 << mjdx.size();
    vector<vector<Double_t> > epoints(eMir,x);

    for (Int_t m=0; m<Int_t(mjdx.size()); m++) {
      Int_t size1 = 1 << m;
      Int_t size2 = 1 << (m+1);
      // copy all previous mirror points
      for (Int_t l=size1; l<size2; ++l) {
   epoints[l] = epoints[l-size1];
   // fill high mirror points
   vector<Double_t>& epoint = epoints[l];
   epoint[mjdx[Int_t(mjdx.size()-1)-m]] = mpoint[mjdx[Int_t(mjdx.size()-1)-m]];
      }
    }

    // remove duplicate mirror points
    // note that: first epoint == x
    epoints.erase(epoints.begin());

    // add mirror points of event i to total dataset
    TVectorD pointR(_nDim);

    for (Int_t m=0; m<Int_t(epoints.size()); m++) {
      _idx.push_back(i);
      _dataPts.push_back(epoints[m]);
      //_weights0.push_back(_weights0[i]);
      for (Int_t j=0; j<_nDim; j++) { pointR[j] = (epoints[m])[j]; }
      if (_nDim>1 && _rotate) { pointR *= *_rotMat; }
      _dataPtsR.push_back(pointR);
    }

    epoints.clear();
    mpoints.clear();
    mjdcs.clear();
  } // end of event loop

  _nEventsM = Int_t(_dataPts.size());
}

////////////////////////////////////////////////////////////////////////////////

void RooNDKeysPdf::loadWeightSet() const
{
  _wMap.clear();

  for (Int_t i=0; i<_nEventsM; i++) {
    _data.get(_idx[i]);
    Double_t myweight = _data.weight();
    //if ( TMath::Abs(myweight)>_minWeight ) {
      _wMap[i] = myweight;
    //}
  }

  coutI(Contents) << "RooNDKeysPdf::loadWeightSet(" << this << ") : Number of weighted events : " << _wMap.size() << endl;
}

////////////////////////////////////////////////////////////////////////////////
/// determine points in +/- nSigma shell around the box determined by the variable
/// ranges. These points are needed in the normalization, to determine probability
/// leakage in and out of the box.

void RooNDKeysPdf::calculateShell(BoxInfo* bi) const
{
  for (Int_t j=0; j<_nDim; j++) {
    if (bi->xVarLo[j]==_xDatLo[j] && bi->xVarHi[j]==_xDatHi[j]) {
      bi->netFluxZ = bi->netFluxZ && kTRUE;
    } else { bi->netFluxZ = kFALSE; }

    bi->xVarLoM3s[j] = bi->xVarLo[j] - _nSigma * (_n * _sigma[j]);
    bi->xVarLoP3s[j] = bi->xVarLo[j] + _nSigma * (_n * _sigma[j]);
    bi->xVarHiM3s[j] = bi->xVarHi[j] - _nSigma * (_n * _sigma[j]);
    bi->xVarHiP3s[j] = bi->xVarHi[j] + _nSigma * (_n * _sigma[j]);

    //cout<<"bi->xVarLoM3s["<<j<<"]="<<bi->xVarLoM3s[j]<<endl;
    //cout<<"bi->xVarLoP3s["<<j<<"]="<<bi->xVarLoP3s[j]<<endl;
    //cout<<"bi->xVarHiM3s["<<j<<"]="<<bi->xVarHiM3s[j]<<endl;
    //cout<<"bi->xVarHiM3s["<<j<<"]="<<bi->xVarHiM3s[j]<<endl;
  }

  map<Int_t,Double_t>::iterator wMapItr = _wMap.begin();

  //for (Int_t i=0; i<_nEventsM; i++) {
  for (; wMapItr!=_wMap.end(); ++wMapItr) {
    Int_t i = (*wMapItr).first;

    const vector<Double_t>& x = _dataPts[i];
    Bool_t inVarRange(kTRUE);
    Bool_t inVarRangePlusShell(kTRUE);

    for (Int_t j=0; j<_nDim; j++) {

      if (x[j]>bi->xVarLo[j] && x[j]<bi->xVarHi[j]) {
   inVarRange = inVarRange && kTRUE;
      } else { inVarRange = inVarRange && kFALSE; }

      if (x[j]>bi->xVarLoM3s[j] && x[j]<bi->xVarHiP3s[j]) {
   inVarRangePlusShell = inVarRangePlusShell && kTRUE;
      } else { inVarRangePlusShell = inVarRangePlusShell && kFALSE; }
    }

    // event in range?
    if (inVarRange) {
      bi->bIdcs.push_back(i);
    }

    // event in shell?
    if (inVarRangePlusShell) {
      bi->bpsIdcs[i] = kTRUE;
      Bool_t inShell(kFALSE);
      for (Int_t j=0; j<_nDim; j++) {
   if ((x[j]>bi->xVarLoM3s[j] && x[j]<bi->xVarLoP3s[j]) && x[j]<(bi->xVarLo[j]+bi->xVarHi[j])/2.) {
     inShell = kTRUE;
   } else if ((x[j]>bi->xVarHiM3s[j] && x[j]<bi->xVarHiP3s[j]) && x[j]>(bi->xVarLo[j]+bi->xVarHi[j])/2.) {
     inShell = kTRUE;
   }
      }
      if (inShell) bi->sIdcs.push_back(i); // needed for normalization
      else {
   bi->bmsIdcs.push_back(i);          // idem
      }
    }
  }

  coutI(Contents) << "RooNDKeysPdf::calculateShell() : "
        << "\n Events in shell " << bi->sIdcs.size()
        << "\n Events in box " << bi->bIdcs.size()
        << "\n Events in box and shell " << bi->bpsIdcs.size()
        << endl;
}

////////////////////////////////////////////////////////////////////////////////
///bi->nEventsBMSW=0.;
///bi->nEventsBW=0.;

void RooNDKeysPdf::calculatePreNorm(BoxInfo* bi) const
{
  // box minus shell
  for (Int_t i=0; i<Int_t(bi->bmsIdcs.size()); i++)
    bi->nEventsBMSW += _wMap[bi->bmsIdcs[i]];

  // box
  for (Int_t i=0; i<Int_t(bi->bIdcs.size()); i++)
    bi->nEventsBW += _wMap[bi->bIdcs[i]];

  cxcoutD(Eval) << "RooNDKeysPdf::calculatePreNorm() : "
         << "\n nEventsBMSW " << bi->nEventsBMSW
         << "\n nEventsBW " << bi->nEventsBW
         << endl;
}

////////////////////////////////////////////////////////////////////////////////
/// sort entries, as needed for loopRange()

void RooNDKeysPdf::sortDataIndices(BoxInfo* bi) const
{
  // will loop over all events by default
  if (!_sortInput) { 
    _ibNoSort.clear();
    for (unsigned int i=0; i<_dataPtsR.size(); ++i) {
      _ibNoSort[i] = kTRUE;
    }
    return;
  }
  
  itVec itrVecR;
  vector<TVectorD>::iterator dpRItr = _dataPtsR.begin();
  for (Int_t i=0; dpRItr!=_dataPtsR.end(); ++dpRItr, ++i) {
    if (bi) {
      if (bi->bpsIdcs.find(i)!=bi->bpsIdcs.end())
      //if (_wMap.find(i)!=_wMap.end())
   itrVecR.push_back(itPair(i,dpRItr));
    } else itrVecR.push_back(itPair(i,dpRItr));
  }

  for (Int_t j=0; j<_nDim; j++) {
    _sortTVIdcs[j].clear();
    sort(itrVecR.begin(),itrVecR.end(),SorterTV_L2H(j));
    _sortTVIdcs[j] = itrVecR;
  }

  for (Int_t j=0; j<_nDim; j++) {
    cxcoutD(Eval) << "RooNDKeysPdf::sortDataIndices() : Number of sorted events : " << _sortTVIdcs[j].size() << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

void RooNDKeysPdf::calculateBandWidth() const
{
  cxcoutD(Eval) << "RooNDKeysPdf::calculateBandWidth()" << endl;

  // non-adaptive bandwidth
  // (default, and needed to calculate adaptive bandwidth)

  if(!_options.Contains("a")) {
      cxcoutD(Eval) << "RooNDKeysPdf::calculateBandWidth() Using static bandwidth." << endl;
  }    
  
  // fixed width approximation
  for (Int_t i=0; i<_nEvents; i++) {
    vector<Double_t>& weight = _weights0[i];
    for (Int_t j=0; j<_nDim; j++) { 
      weight[j] = _n * (*_sigmaR)[j] ; 
      //cout<<"j: "<<j<<", _n: "<<_n<<", sigmaR="<<(*_sigmaR)[j]<<", weight="<<weight[j]<<endl;
    }
  }

  // adaptive width
  if( _options.Contains("a") ) {
    cxcoutD(Eval) << "RooNDKeysPdf::calculateBandWidth() Using adaptive bandwidth." << endl;

    double sqrt12=sqrt(12.);
    double sqrtSigmaAvgR=sqrt(_sigmaAvgR);

    vector<Double_t> dummy(_nDim,0.);
    _weights1.resize(_nEvents,dummy);

    std::vector<std::vector<Double_t> >* weights_prev(0);
    std::vector<std::vector<Double_t> >* weights_new(0);

    //cout << "Number of adaptive iterations: " << _nAdpt << endl;

    for (Int_t k=1; k<=_nAdpt; ++k) {  

      //cout << "  Cycle: " << k << endl;

      // if multiple adaptive iterations, need to swap weight sets 
      if (k%2) {
	weights_prev = &_weights0;
	weights_new  = &_weights1;
      } else {
	weights_prev = &_weights1;
	weights_new  = &_weights0;
      }

      for(Int_t i=0; i<_nEvents; ++i) {	
	vector<Double_t>& x = _dataPts[i];
	Double_t f =  TMath::Power( gauss(x,*weights_prev)/_nEventsW , -1./(2.*_d) ) ;
	
	vector<Double_t>& weight = (*weights_new)[i];
	for (Int_t j=0; j<_nDim; j++) {
	  Double_t norm = (_n*(*_sigmaR)[j]) / sqrtSigmaAvgR ; 
	  weight[j] = norm * f / sqrt12 ;  //  note additional factor of sqrt(12) compared with HEP-EX/0011057
	}
      }
    }
    // this is the latest updated weights set
    _weights = weights_new;
  } 
}

////////////////////////////////////////////////////////////////////////////////
/// loop over all closest point to x, as determined by loopRange()

Double_t RooNDKeysPdf::gauss(vector<Double_t>& x, vector<vector<Double_t> >& weights) const
{
  if(_nEvents==0) return 0.;

  Double_t z=0.;
  map<Int_t,Bool_t> ibMap;

  // determine input loop range for event x
  if (_sortInput) { loopRange(x,ibMap); }

  map<Int_t,Bool_t>::iterator ibMapItr, ibMapEnd;
  ibMapItr = ( _sortInput ? ibMap.begin() : _ibNoSort.begin() ) ;
  ibMapEnd = ( _sortInput ? ibMap.end()   : _ibNoSort.end() ) ;

  for (; ibMapItr!=ibMapEnd; ++ibMapItr) {
    Int_t i = (*ibMapItr).first;

    Double_t g(1.);

    if(i>=(Int_t)_idx.size()) {continue;} //---> 1.myline

    const vector<Double_t>& point  = _dataPts[i];
    const vector<Double_t>& weight = weights[_idx[i]];

    for (Int_t j=0; j<_nDim; j++) {
      (*_dx)[j] = x[j]-point[j];
    }

    if (_nDim>1 && _rotate) {
      *_dx *= *_rotMat; // rotate to decorrelated frame!
    }

    for (Int_t j=0; j<_nDim; j++) {
      Double_t r = (*_dx)[j];  // x[j] - point[j];
      Double_t c = 1./(2.*weight[j]*weight[j]);

      //cout << "j = " << j << " x[j] = " << point[j] << " w = " << weight[j] << endl;

      g *= exp( -c*r*r );
      g *= 1./(_sqrt2pi*weight[j]);
    }
    z += (g*_wMap[_idx[i]]);
  }
  return z;
}

////////////////////////////////////////////////////////////////////////////////
/// determine closest points to x, to loop over in evaluate()

void RooNDKeysPdf::loopRange(vector<Double_t>& x, map<Int_t,Bool_t>& ibMap) const
{
  ibMap.clear();
  TVectorD xRm(_nDim);
  TVectorD xRp(_nDim);

  for (Int_t j=0; j<_nDim; j++) { xRm[j] = xRp[j] = x[j]; }

  if (_nDim>1 && _rotate) {
    xRm *= *_rotMat;
    xRp *= *_rotMat;
  }
  for (Int_t j=0; j<_nDim; j++) {
    xRm[j] -= _nSigma * (_n * (*_sigmaR)[j]);
    xRp[j] += _nSigma * (_n * (*_sigmaR)[j]);
    //cout<<"xRm["<<j<<"]="<<xRm[j]<<endl;
    //cout<<"xRp["<<j<<"]="<<xRp[j]<<endl;
  }

  vector<TVectorD> xvecRm(1,xRm);
  vector<TVectorD> xvecRp(1,xRp);

  map<Int_t,Bool_t> ibMapRT;

  for (Int_t j=0; j<_nDim; j++) {
    ibMap.clear();
    itVec::iterator lo = lower_bound(_sortTVIdcs[j].begin(), _sortTVIdcs[j].end(),
                 itPair(0,xvecRm.begin()), SorterTV_L2H(j));
    itVec::iterator hi = upper_bound(_sortTVIdcs[j].begin(), _sortTVIdcs[j].end(),
				     itPair(0,xvecRp.begin()), SorterTV_L2H(j));
    itVec::iterator it = lo;

    if (j==0) {
      if (_nDim==1) { for (it=lo; it!=hi; ++it) ibMap[(*it).first] = kTRUE; }
      else { for (it=lo; it!=hi; ++it) ibMapRT[(*it).first] = kTRUE; }
      continue;
    }

    for (it=lo; it!=hi; ++it)
      if (ibMapRT.find((*it).first)!=ibMapRT.end()) { ibMap[(*it).first] = kTRUE; }

    ibMapRT.clear();
    if (j!=_nDim-1) { ibMapRT = ibMap; }
  }
}

////////////////////////////////////////////////////////////////////////////////

void RooNDKeysPdf::boxInfoInit(BoxInfo* bi, const char* rangeName, Int_t /*code*/) const
{
  vector<Bool_t> doInt(_nDim,kTRUE);

  bi->filled = kFALSE;

  bi->xVarLo.resize(_nDim,0.);
  bi->xVarHi.resize(_nDim,0.);
  bi->xVarLoM3s.resize(_nDim,0.);
  bi->xVarLoP3s.resize(_nDim,0.);
  bi->xVarHiM3s.resize(_nDim,0.);
  bi->xVarHiP3s.resize(_nDim,0.);

  bi->netFluxZ = kTRUE;
  bi->bpsIdcs.clear();
  bi->bIdcs.clear();
  bi->sIdcs.clear();
  bi->bmsIdcs.clear();

  bi->nEventsBMSW=0.;
  bi->nEventsBW=0.;

  _varItr->Reset() ;
  RooRealVar* var ;
  for(Int_t j=0; (var=(RooRealVar*)_varItr->Next()); ++j) {
    if (doInt[j]) {
      bi->xVarLo[j] = var->getMin(rangeName);
      bi->xVarHi[j] = var->getMax(rangeName);
    } else {
      bi->xVarLo[j] = var->getVal() ;
      bi->xVarHi[j] = var->getVal() ;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooNDKeysPdf::evaluate() const
{
  if ( _tracker!=NULL && _tracker->hasChanged(kTRUE) ) {
    updateRho(); // update internal rho parameters
    // redetermine static and/or adaptive bandwidth
    calculateBandWidth();
  } 

  _varItr->Reset() ;
  RooAbsReal* var ;
  const RooArgSet* nset = _varList.nset() ;
  for(Int_t j=0; (var=(RooAbsReal*)_varItr->Next()); ++j) {
    _x[j] = var->getVal(nset);
  }

  Double_t val = gauss(_x,*_weights);
  //cout<<"returning "<<val<<endl;

  if (val>=1E-20)
    return val ;
  else
    return (1E-20) ;
}

////////////////////////////////////////////////////////////////////////////////

Int_t RooNDKeysPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  if (rangeName) return 0 ;

  Int_t code=0;
  if (matchArgs(allVars,analVars,RooArgSet(_varList))) { code=1; }

  return code;

}

////////////////////////////////////////////////////////////////////////////////

Double_t RooNDKeysPdf::analyticalIntegral(Int_t code, const char* rangeName) const
{
  cxcoutD(Eval) << "Calling RooNDKeysPdf::analyticalIntegral(" << GetName() << ") with code " << code
         << " and rangeName " << (rangeName?rangeName:"<none>") << endl;

  // determine which observables need to be integrated over ...
  Int_t nComb = 1 << (_nDim);
  R__ASSERT(code>=1 && code<nComb) ;

  vector<Bool_t> doInt(_nDim,kTRUE);

  // get BoxInfo
  BoxInfo* bi(0);

  if (rangeName) {
    string rangeNameStr(rangeName) ;
    bi = _rangeBoxInfo[make_pair(rangeNameStr,code)] ;
    if (!bi) {
      bi = new BoxInfo ;
      _rangeBoxInfo[make_pair(rangeNameStr,code)] = bi ;
      boxInfoInit(bi,rangeName,code);
    }
  } else bi= &_fullBoxInfo ;

  // have boundaries changed?
  Bool_t newBounds(kFALSE);
  _varItr->Reset() ;
  RooRealVar* var ;
  for(Int_t j=0; (var=(RooRealVar*)_varItr->Next()); ++j) {
    if ((var->getMin(rangeName)-bi->xVarLo[j]!=0) ||
   (var->getMax(rangeName)-bi->xVarHi[j]!=0)) {
      newBounds = kTRUE;
    }
  }

  // reset
  if (newBounds) {
    cxcoutD(Eval) << "RooNDKeysPdf::analyticalIntegral() : Found new boundaries ... " << (rangeName?rangeName:"<none>") << endl;
    boxInfoInit(bi,rangeName,code);
  }

  // recalculates netFluxZero and nEventsIR
  if (!bi->filled || newBounds) {
    // Fill box info with contents
    calculateShell(bi);
    calculatePreNorm(bi);
    bi->filled = kTRUE;
    sortDataIndices(bi);
  }

  // first guess
  Double_t norm=bi->nEventsBW;

  if (_mirror && bi->netFluxZ) {
    // KEYS expression is self-normalized
    cxcoutD(Eval) << "RooNDKeysPdf::analyticalIntegral() : Using mirrored normalization : " << bi->nEventsBW << endl;
    return bi->nEventsBW;
  }
  // calculate leakage in and out of variable range box
  else
  {
    norm = bi->nEventsBMSW;
    if (norm<0.) norm=0.;

    for (Int_t i=0; i<Int_t(bi->sIdcs.size()); ++i) {
      Double_t prob=1.;
      const vector<Double_t>& x = _dataPts[bi->sIdcs[i]];
      const vector<Double_t>& weight = (*_weights)[_idx[bi->sIdcs[i]]];

      vector<Double_t> chi(_nDim,100.);

      for (Int_t j=0; j<_nDim; j++) {
   if(!doInt[j]) continue;

   if ((x[j]>bi->xVarLoM3s[j] && x[j]<bi->xVarLoP3s[j]) && x[j]<(bi->xVarLo[j]+bi->xVarHi[j])/2.)
     chi[j] = (x[j]-bi->xVarLo[j])/weight[j];
   else if ((x[j]>bi->xVarHiM3s[j] && x[j]<bi->xVarHiP3s[j]) && x[j]>(bi->xVarLo[j]+bi->xVarHi[j])/2.)
     chi[j] = (bi->xVarHi[j]-x[j])/weight[j];

   if (chi[j]>0) // inVarRange
     prob *= (0.5 + TMath::Erf(fabs(chi[j])/sqrt(2.))/2.);
   else // outside Var range
     prob *= (0.5 - TMath::Erf(fabs(chi[j])/sqrt(2.))/2.);
      }

      norm += prob * _wMap[_idx[bi->sIdcs[i]]];
    }

    cxcoutD(Eval) << "RooNDKeysPdf::analyticalIntegral() : Final normalization : " << norm << " " << bi->nEventsBW << endl;
    return norm;
  }
}



RooDataSet* 
////////////////////////////////////////////////////////////////////////////////

RooNDKeysPdf::createDatasetFromHist( const RooArgList& varList, const TH1& hist ) const
{
  std::vector<RooRealVar*> varVec;
  RooArgSet varsAndWeightSet;

  TIterator* varItr = varList.createIterator() ;
  RooAbsArg* var ;
  for (Int_t i=0; (var = (RooAbsArg*)varItr->Next()); ++i) {
    if (!dynamic_cast<RooRealVar*>(var)) {
      coutE(InputArguments) << "RooNDKeysPdf::createDatasetFromHist(" << GetName() << ") WARNING: variable " << var->GetName() 
			    << " is not of type RooRealVar. Skip." << endl ;
      continue;
    }
    varsAndWeightSet.add(*var) ; // used for dataset creation
    varVec.push_back( static_cast<RooRealVar*>(var) ); // used for setting the variables.
  }
  delete varItr ;

  /// Add weight 
  RooRealVar weight("weight", "event weight", 0);
  varsAndWeightSet.add(weight);

  /// determine histogram dimensionality 
  unsigned int histndim(0);
  std::string classname = hist.ClassName();
  if      (classname.find("TH1")==0) { histndim=1; }
  else if (classname.find("TH2")==0) { histndim=2; }
  else if (classname.find("TH3")==0) { histndim=3; }
  assert ( histndim==varVec.size() );

  if ( histndim>3 || histndim<=0 ) {
    coutE(InputArguments) << "RooNDKeysPdf::createDatasetFromHist(" << GetName() 
			  << ") ERROR: input histogram dimension not between [1-3]: " << histndim
			  << endl ;
    assert(0);
  }

  /// dataset creation
  RooDataSet* dataFromHist = new RooDataSet( "datasetFromHist", "datasetFromHist", varsAndWeightSet, weight.GetName() );

  /// dataset filling
  for (int i=1; i<=hist.GetXaxis()->GetNbins(); ++i) { 
    // 1 or more dimension
    
    Double_t xval = hist.GetXaxis()->GetBinCenter(i);
    varVec[0]->setVal( xval );
    
    if(varVec.size()==1) {
      Double_t fval = hist.GetBinContent(i);
      weight.setVal(fval);
      dataFromHist->add( varsAndWeightSet, fval );
    } else { // 2 or more dimensions
      
      for(int j=1; j<=hist.GetYaxis()->GetNbins(); ++j) {
	Double_t yval = hist.GetYaxis()->GetBinCenter(j);
	varVec[1]->setVal( yval );
	
	if(varVec.size()==2) { 
	  Double_t fval = hist.GetBinContent(i,j);
	  weight.setVal(fval);
	  dataFromHist->add( varsAndWeightSet, fval );
	} else { // 3 dimensions 
	  
	  for(int k=1; k<=hist.GetZaxis()->GetNbins(); ++k) {
	    Double_t zval = hist.GetZaxis()->GetBinCenter(k);
	    varVec[2]->setVal( zval );

	    Double_t fval = hist.GetBinContent(i,j,k);
	    weight.setVal(fval);
	    dataFromHist->add( varsAndWeightSet, fval );
	  }
	}
      }
    }
  }

  return dataFromHist;
}



TMatrixD 
////////////////////////////////////////////////////////////////////////////////
/// Return evaluated weights

RooNDKeysPdf::getWeights( const int& k ) const
{
  TMatrixD mref( _nEvents, _nDim+1 );

  cxcoutD(Eval) << "RooNDKeysPdf::getWeights() Return evaluated weights." << endl;
  
  for(Int_t i=0; i<_nEvents; ++i) {    
    vector<Double_t>& x = _dataPts[i];
    for (Int_t j=0; j<_nDim; j++) {
      mref(i,j) = x[j];
    }
    
    vector<Double_t>& weight = (*_weights)[i];
    mref(i,_nDim) = weight[k];
  }
  
  return mref;
}



void 
////////////////////////////////////////////////////////////////////////////////

RooNDKeysPdf::updateRho() const
{
  _rhoItr->Reset();
  RooAbsReal* rho(0) ;
  for(Int_t j=0; (rho=(RooAbsReal*)_rhoItr->Next()); ++j) {
    _rho[j] = rho->getVal();
  }

  if (_nDim>1 && _rotate) {
    TMatrixDSym covMatRho(_nDim); // covariance matrix times rho parameters
    for (Int_t j=0; j<_nDim; j++) {
      for (Int_t k=0; k<_nDim; k++) { 
	covMatRho(j,k) = (*_covMat)(j,k) * _rho[j] * _rho[k];
      }
    }
    // find decorrelation matrix and eigenvalues (R)
    TMatrixDSymEigen evCalculatorRho(covMatRho);
    *_sigmaR = evCalculatorRho.GetEigenValues();
    for (Int_t j=0; j<_nDim; j++) { (*_sigmaR)[j] = sqrt((*_sigmaR)[j]); }
  } else {
    for (Int_t j=0; j<_nDim; j++) { (*_sigmaR)[j] = (_sigma[j] * _rho[j]); } // * rho
  }
}


