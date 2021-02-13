// @(#)root/roostats:$Id:  cranmer $
// Author: Kyle Cranmer, George Lewis 
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////

/** \class ParamHistFunc
 * \ingroup HistFactory
 *   A class which maps the current values of a RooRealVar
 *  (or a set of RooRealVars) to one of a number of RooRealVars:
 *
 *  `ParamHistFunc: {val1, val2, ...} -> {gamma (RooRealVar)}`
 *
 *  The intended interpretation is that each parameter in the
 *  range represent the height of a bin over the domain
 *  space.
 *
 *  The 'createParamSet' is an easy way to create these
 *  parameters from a set of observables. They are
 *  stored using the "TH1" ordering convention (as compared
 *  to the RooDataHist convention, which is used internally
 *  and one must map between the two).
 *
 *  All indices include '0':<br>
 *  \f$ \gamma_{i,j} \f$ = `paramSet[ size(i)*j + i ]`
 *
 *  ie assuming the dimensions are 5*5:<br>
 *  \f$ \gamma_{2,1} \f$ = `paramSet[ 5*1 + 2 ] = paramSet[7]`
 */


#include <sstream>
#include <math.h>
#include <stdexcept>
#include <iostream>

#include "TH1.h"

#include "RooFit.h"
#include "RooStats/HistFactory/ParamHistFunc.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"

#include "RooConstVar.h"
#include "RooBinning.h"
#include "RooErrorHandler.h"

#include "RooGaussian.h"
#include "RooHistFunc.h"
#include "RooArgSet.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMsgService.h"

// Forward declared:
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooWorkspace.h"

//using namespace std;

ClassImp(ParamHistFunc);


////////////////////////////////////////////////////////////////////////////////

ParamHistFunc::ParamHistFunc() : _numBins(0)
{
  _dataSet.removeSelfFromDir(); // files must not delete _dataSet.
}


////////////////////////////////////////////////////////////////////////////////
/// Create a function which returns binewise-values
/// This class contains N RooRealVar's, one for each
/// bin from the given RooRealVar.
///
/// The value of the function in the ith bin is 
/// given by:
///
/// F(i) = gamma_i * nominal(i)
///
/// Where the nominal values are simply fixed
/// numbers (default = 1.0 for all i)
ParamHistFunc::ParamHistFunc(const char* name, const char* title, 
			     const RooArgList& vars, const RooArgList& paramSet) :
  RooAbsReal(name, title),
  _dataVars("!dataVars","data Vars",       this),
  _paramSet("!paramSet","bin parameters",  this),
  _numBins(0),
  _dataSet( (std::string(name)+"_dataSet").c_str(), "", vars)
{

  // Create the dataset that stores the binning info:
  
  //  _dataSet = RooDataSet("

  _dataSet.removeSelfFromDir(); // files must not delete _dataSet.

  // Set the binning
  // //_binning = var.getBinning().clone() ;
  
  // Create the set of parameters
  // controlling the height of each bin

  // Get the number of bins
  _numBins = GetNumBins( vars );

  // Add the parameters (with checking)
  addVarSet( vars );
  addParamSet( paramSet );
}


////////////////////////////////////////////////////////////////////////////////
/// Create a function which returns bin-wise values.
/// This class allows to multiply bin contents of histograms
/// with the values of a set of RooRealVars.
///
/// The value of the function in the ith bin is 
/// given by:
/// \f[
///   F(i) = \gamma_{i} * \mathrm{nominal}(i)
/// \f]
///
/// Where the nominal values are taken from the histogram,
/// and the \f$ \gamma_{i} \f$ can be set from the outside.
ParamHistFunc::ParamHistFunc(const char* name, const char* title, 
			     const RooArgList& vars, const RooArgList& paramSet,
			     const TH1* Hist ) :
  RooAbsReal(name, title),
  //  _dataVar("!dataVar","data Var", this, (RooRealVar&) var),
  _dataVars("!dataVars","data Vars",       this),
  _paramSet("!paramSet","bin parameters",  this),
  _numBins(0),
  _dataSet( (std::string(name)+"_dataSet").c_str(), "", vars, Hist)
{

  _dataSet.removeSelfFromDir(); // files must not delete _dataSet.

  // Get the number of bins
  _numBins = GetNumBins( vars );

  // Add the parameters (with checking)
  addVarSet( vars );
  addParamSet( paramSet );
 
}


Int_t ParamHistFunc::GetNumBins( const RooArgSet& vars ) {
  
  // A helper method to get the number of bins
  
  if( vars.getSize() == 0 ) return 0;
    
  Int_t numBins = 1;

  RooFIter varIter = vars.fwdIterator() ;
  RooAbsArg* comp ;
  while((comp = (RooAbsArg*) varIter.next())) {
    if (!dynamic_cast<RooRealVar*>(comp)) {
      std::cout << "ParamHistFunc::GetNumBins" << vars.GetName() << ") ERROR: component " 
	   << comp->GetName() 
	   << " in vars list is not of type RooRealVar" << std::endl ;
      RooErrorHandler::softAbort() ;
      return -1;
    }
    RooRealVar* var = (RooRealVar*) comp;

    Int_t varNumBins = var->numBins();
    numBins *= varNumBins;
  }
    
  return numBins;

}


////////////////////////////////////////////////////////////////////////////////

ParamHistFunc::ParamHistFunc(const ParamHistFunc& other, const char* name) :
  RooAbsReal(other, name), 
  _dataVars("!dataVars", this, other._dataVars ),
  _paramSet("!paramSet", this, other._paramSet),
  _numBins( other._numBins ),
  _binMap( other._binMap ),
  _dataSet( other._dataSet )
{
  _dataSet.removeSelfFromDir(); // files must not delete _dataSet.

  // Copy constructor
  // Member _ownedList is intentionally not copy-constructed -- ownership is not transferred
}


////////////////////////////////////////////////////////////////////////////////

ParamHistFunc::~ParamHistFunc() 
{
  ;
}


////////////////////////////////////////////////////////////////////////////////
/// Get the index of the gamma parameter associated
/// with the current bin.
/// This number is the "RooDataSet" style index
/// and it must be because it uses the RooDataSet method directly
/// This is intended to be fed into the getParameter(Int_t) method:
///
/// RooRealVar currentParam = getParameter( getCurrentBin() );
Int_t ParamHistFunc::getCurrentBin() const {
  Int_t dataSetIndex = _dataSet.getIndex( _dataVars ); // calcTreeIndex();
  return dataSetIndex;

}


////////////////////////////////////////////////////////////////////////////////
/// Get the parameter associate with the the
/// input RooDataHist style index
/// It uses the binMap to convert the RooDataSet style index
/// into the TH1 style index (which is how they are stored
/// internally in the '_paramSet' vector
RooRealVar& ParamHistFunc::getParameter( Int_t index ) const {
  Int_t gammaIndex = -1;
  if( _binMap.find( index ) != _binMap.end() ) {
    gammaIndex = _binMap[ index ];
  }
  else {
    std::cout << "Error: ParamHistFunc internal bin index map "
	      << "not properly configured" << std::endl;
    throw -1;
  }

  return (RooRealVar&) _paramSet[gammaIndex];
}


////////////////////////////////////////////////////////////////////////////////

RooRealVar& ParamHistFunc::getParameter() const {
  Int_t index = getCurrentBin();
  return getParameter( index );
}


void ParamHistFunc::setParamConst( Int_t index, Bool_t varConst ) {
  RooRealVar& var = getParameter( index );
  var.setConstant( varConst );
}


void ParamHistFunc::setConstant( bool constant ) {
  for( int i=0; i < numBins(); ++i) {
    setParamConst(i, constant);
  }
}


////////////////////////////////////////////////////////////////////////////////

void ParamHistFunc::setShape( TH1* shape ) {
  int num_hist_bins = shape->GetNbinsX()*shape->GetNbinsY()*shape->GetNbinsZ();

  if( num_hist_bins != numBins() ) {
    std::cout << "Error - ParamHistFunc: cannot set Shape of ParamHistFunc: " << GetName()
	      << " using histogram: " << shape->GetName()
	      << ". Bins don't match" << std::endl;
    throw std::runtime_error("setShape");
  }


  Int_t TH1BinNumber = 0;
  for( Int_t i = 0; i < numBins(); ++i) {
    
    TH1BinNumber++;
    
    while( shape->IsBinUnderflow(TH1BinNumber) || shape->IsBinOverflow(TH1BinNumber) ){
      TH1BinNumber++;
    }

    //RooRealVar& var = dynamic_cast<RooRealVar&>(getParameter(i));
    RooRealVar& var = dynamic_cast<RooRealVar&>(_paramSet[i]);
    var.setVal( shape->GetBinContent(TH1BinNumber) );
  }

}


////////////////////////////////////////////////////////////////////////////////
/// Create the list of RooRealVar
/// parameters which represent the
/// height of the histogram bins.
/// The list 'vars' represents the 
/// observables (corresponding to histogram bins)
/// that these newly created parameters will 
/// be mapped to. (ie, we create one parameter
/// per observable in vars and per bin in each observable)

/// Store them in a list using:
/// _paramSet.add( createParamSet() );
/// This list is stored in the "TH1" index order
RooArgList ParamHistFunc::createParamSet(RooWorkspace& w, const std::string& Prefix, 
					 const RooArgList& vars) {


  // Get the number of bins
  // in the nominal histogram

  RooArgList paramSet;

  Int_t numVars = vars.getSize();
  Int_t numBins = GetNumBins( vars );

  if( numVars == 0 ) {
    std::cout << "Warning - ParamHistFunc::createParamSet() :"
	 << " No Variables provided.  Not making constraint terms." 
	 << std::endl;
    return paramSet;
  }

  else if( numVars == 1 ) {
 
    // For each bin, create a RooRealVar
    for( Int_t i = 0; i < numBins; ++i) {

      std::stringstream VarNameStream;
      VarNameStream << Prefix << "_bin_" << i;
      std::string VarName = VarNameStream.str();

      RooRealVar gamma( VarName.c_str(), VarName.c_str(), 1.0 ); 
      // "Hard-Code" a minimum of 0.0
      gamma.setMin( 0.0 );
      gamma.setConstant( false );

      w.import( gamma, RooFit::RecycleConflictNodes() );
      RooRealVar* gamma_wspace = (RooRealVar*) w.var( VarName.c_str() );

      paramSet.add( *gamma_wspace );

    }
  }

  else if( numVars == 2 ) {
 
    // Create a vector of indices
    // all starting at 0
    std::vector< Int_t > Indices(numVars, 0);

    RooRealVar* varx = (RooRealVar*) vars.at(0);
    RooRealVar* vary = (RooRealVar*) vars.at(1);
    
    // For each bin, create a RooRealVar
    for( Int_t j = 0; j < vary->numBins(); ++j) {
      for( Int_t i = 0; i < varx->numBins(); ++i) {

	// Ordering is important:
	// To match TH1, list goes over x bins
	// first, then y

	std::stringstream VarNameStream;
	VarNameStream << Prefix << "_bin_" << i << "_" << j;
	std::string VarName = VarNameStream.str();

	RooRealVar gamma( VarName.c_str(), VarName.c_str(), 1.0 ); 
	// "Hard-Code" a minimum of 0.0
	gamma.setMin( 0.0 );
	gamma.setConstant( false );
	  
	w.import( gamma, RooFit::RecycleConflictNodes() );
	RooRealVar* gamma_wspace = (RooRealVar*) w.var( VarName.c_str() );
	  
	paramSet.add( *gamma_wspace );
	  
      }
    }
  }

  else if( numVars == 3 ) {
 
    // Create a vector of indices
    // all starting at 0
    std::vector< Int_t > Indices(numVars, 0);

    RooRealVar* varx = (RooRealVar*) vars.at(0);
    RooRealVar* vary = (RooRealVar*) vars.at(1);
    RooRealVar* varz = (RooRealVar*) vars.at(2);
    
    // For each bin, create a RooRealVar
    for( Int_t k = 0; k < varz->numBins(); ++k) {
      for( Int_t j = 0; j < vary->numBins(); ++j) {
	for( Int_t i = 0; i < varx->numBins(); ++i) {
	  
	  // Ordering is important:
	  // To match TH1, list goes over x bins
	  // first, then y, then z

	  std::stringstream VarNameStream;
	  VarNameStream << Prefix << "_bin_" << i << "_" << j << "_" << k;
	  std::string VarName = VarNameStream.str();
	
	  RooRealVar gamma( VarName.c_str(), VarName.c_str(), 1.0 ); 
	  // "Hard-Code" a minimum of 0.0
	  gamma.setMin( 0.0 );
	  gamma.setConstant( false );
	  
	  w.import( gamma, RooFit::RecycleConflictNodes() );
	  RooRealVar* gamma_wspace = (RooRealVar*) w.var( VarName.c_str() );
	  
	  paramSet.add( *gamma_wspace );
	
	}
      }
    }
  }

  else {
    std::cout << " Error: ParamHistFunc doesn't support dimensions > 3D " <<  std::endl;
  }

  return paramSet;  

}


////////////////////////////////////////////////////////////////////////////////
/// Create the list of RooRealVar parameters which scale the
/// height of histogram bins.
/// The list `vars` represents the observables (corresponding to histogram bins)
/// that these newly created parameters will
/// be mapped to. *I.e.*, we create one parameter
/// per observable in `vars` and per bin in each observable.
///
/// The new parameters are initialised to 1 with an uncertainty of +/- 1.,
/// their range is set to the function arguments.
///
/// Store the parameters in a list using:
/// ```
/// _paramSet.add( createParamSet() );
/// ```
/// This list is stored in the "TH1" index order.
RooArgList ParamHistFunc::createParamSet(RooWorkspace& w, const std::string& Prefix, 
					 const RooArgList& vars, 
					 Double_t gamma_min, Double_t gamma_max) {



  RooArgList params = ParamHistFunc::createParamSet( w, Prefix, vars );

  for (auto comp : params) {
    auto var = static_cast<RooRealVar*>(comp);

    var->setMin( gamma_min );
    var->setMax( gamma_max );
  }

  return params;

}


////////////////////////////////////////////////////////////////////////////////
/// Create the list of RooRealVar
/// parameters which represent the
/// height of the histogram bins.
/// Store them in a list
RooArgList ParamHistFunc::createParamSet(const std::string& Prefix, Int_t numBins, 
					 Double_t gamma_min, Double_t gamma_max) {


  // _paramSet.add( createParamSet() );

  // Get the number of bins
  // in the nominal histogram

  RooArgList paramSet;

  if( gamma_max <= gamma_min ) {

    std::cout << "Warning: gamma_min <= gamma_max: Using default values (0, 10)" << std::endl;

    gamma_min = 0.0;
    gamma_max = 10.0;

  }

  Double_t gamma_nominal = 1.0;

  if( gamma_nominal < gamma_min ) {
    gamma_nominal = gamma_min;
  }

  if( gamma_nominal > gamma_max ) {
    gamma_nominal = gamma_max;
  }

  // For each bin, create a RooRealVar
  for( Int_t i = 0; i < numBins; ++i) {

    std::stringstream VarNameStream;
    VarNameStream << Prefix << "_bin_" << i;
    std::string VarName = VarNameStream.str();

    RooRealVar* gamma = new RooRealVar( VarName.c_str(), VarName.c_str(), 
					gamma_nominal, gamma_min, gamma_max );
    gamma->setConstant( false );
    paramSet.add( *gamma );

  }

  return paramSet;

}


////////////////////////////////////////////////////////////////////////////////
/// return 0 for success
/// return 1 for failure
/// Check that the elements 
/// are actually RooRealVar's
/// If so, add them to the 
/// list of vars
Int_t ParamHistFunc::addVarSet( const RooArgList& vars ) {


  int numVars = 0;

  RooFIter varIter = vars.fwdIterator() ;
  RooAbsArg* comp ;
  while((comp = (RooAbsArg*) varIter.next())) {
    if (!dynamic_cast<RooRealVar*>(comp)) {
      coutE(InputArguments) << "ParamHistFunc::(" << GetName() << ") ERROR: component " 
			    << comp->GetName() << " in variables list is not of type RooRealVar" 
			    << std::endl;
      RooErrorHandler::softAbort() ;
      return 1;
    }

    _dataVars.add( *comp );
    numVars++;

  }

  Int_t numBinsX = 1;
  Int_t numBinsY = 1;
  Int_t numBinsZ = 1;

  if( numVars == 1 ) {
    RooRealVar* varX = (RooRealVar*) _dataVars.at(0);
    numBinsX = varX->numBins();
    numBinsY = 1;
    numBinsZ = 1;
  } else  if( numVars == 2 ) {
    RooRealVar* varX = (RooRealVar*) _dataVars.at(0);
    RooRealVar* varY = (RooRealVar*) _dataVars.at(1);
    numBinsX = varX->numBins();
    numBinsY = varY->numBins();
    numBinsZ = 1;
  } else  if( numVars == 3 ) {
    RooRealVar* varX = (RooRealVar*) _dataVars.at(0);
    RooRealVar* varY = (RooRealVar*) _dataVars.at(1);
    RooRealVar* varZ = (RooRealVar*) _dataVars.at(2);
    numBinsX = varX->numBins();
    numBinsY = varY->numBins();
    numBinsZ = varZ->numBins();
  } else {
    std::cout << "ParamHistFunc() - Only works for 1-3 variables (1d-3d)" << std::endl;
    throw -1;  
  }

  // Fill the mapping between
  // RooDataHist bins and TH1 Bins:

  // Clear the map
  _binMap.clear();

  // Fill the map
  for( Int_t i = 0; i < numBinsX; ++i ) {
    for( Int_t j = 0; j < numBinsY; ++j ) {
      for( Int_t k = 0; k < numBinsZ; ++k ) {
	
	Int_t RooDataSetBin = k + j*numBinsZ + i*numBinsY*numBinsZ; 
	Int_t TH1HistBin    = i + j*numBinsX + k*numBinsX*numBinsY; 
	  
	_binMap[RooDataSetBin] = TH1HistBin;
	
      }
    }
  }
  
  return 0;

}


////////////////////////////////////////////////////////////////////////////////

Int_t ParamHistFunc::addParamSet( const RooArgList& params ) {
  // return 0 for success
  // return 1 for failure

  // Check that the supplied list has
  // the right number of arguments:

  Int_t numVarBins  = _numBins;
  Int_t numElements = params.getSize();

  if( numVarBins != numElements ) {
    std::cout << "ParamHistFunc::addParamSet - ERROR - " 
	      << "Supplied list of parameters " << params.GetName()
	      << " has " << numElements << " elements but the ParamHistFunc"
	      << GetName() << " has " << numVarBins << " bins."
	      << std::endl;
    return 1;

  }
  
  // Check that the elements 
  // are actually RooRealVar's
  // If so, add them to the 
  // list of params

  RooFIter paramIter = params.fwdIterator() ;
  RooAbsArg* comp ;
  while((comp = (RooAbsArg*) paramIter.next())) {
    if (!dynamic_cast<RooRealVar*>(comp)) {
      coutE(InputArguments) << "ParamHistFunc::(" << GetName() << ") ERROR: component " 
			    << comp->GetName() << " in paramater list is not of type RooRealVar" 
			    << std::endl;
      RooErrorHandler::softAbort() ;
      return 1;
    }

    _paramSet.add( *comp );

  }
  
  return 0;

}


////////////////////////////////////////////////////////////////////////////////

Double_t ParamHistFunc::evaluate() const 
{
  // Find the bin cooresponding to the current
  // value of the RooRealVar:

  RooRealVar* param = (RooRealVar*) &(getParameter());
  Double_t value = param->getVal();
  return value;
  
}


////////////////////////////////////////////////////////////////////////////////
/// Advertise that all integrals can be handled internally.

Int_t ParamHistFunc::getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, 
					     const RooArgSet* normSet, 
					     const char* /*rangeName*/) const 
{
  // Handle trivial no-integration scenario
  if (allVars.getSize()==0) return 0 ;
  if (_forceNumInt) return 0 ;


  // Select subset of allVars that are actual dependents
  analVars.add(allVars) ;

  // Check if this configuration was created before
  Int_t sterileIdx(-1) ;
  CacheElem* cache = (CacheElem*) _normIntMgr.getObj(normSet,&analVars,&sterileIdx,(const char*)0) ;
  if (cache) {
    return _normIntMgr.lastIndex()+1 ;
  }
  
  // Create new cache element
  cache = new CacheElem ;

  // Store cache element
  Int_t code = _normIntMgr.setObj(normSet,&analVars,(RooAbsCacheElement*)cache,0) ;

  return code+1 ; 

}


////////////////////////////////////////////////////////////////////////////////
/// Implement analytical integrations by doing appropriate weighting from  component integrals
/// functions to integrators of components

Double_t ParamHistFunc::analyticalIntegralWN(Int_t /*code*/, const RooArgSet* /*normSet2*/,
					     const char* /*rangeName*/) const 
{
  Double_t value(0) ;

  // Simply loop over bins, 
  // get the height, and
  // multiply by the bind width
  
  RooFIter paramIter = _paramSet.fwdIterator();
  RooRealVar* param = NULL;
  Int_t nominalItr = 0;
  while((param = (RooRealVar*) paramIter.next())) {

    // Get the gamma's value
    Double_t paramVal  = (*param).getVal();
    
    // Get the bin volume
    _dataSet.get( nominalItr );
    Double_t binVolumeDS  = _dataSet.binVolume(); //_binning->binWidth( nominalItr );
    
    // Finally, get the subtotal
    value += paramVal*binVolumeDS;

    ++nominalItr;

    /*
    std::cout << "Integrating : "
	      << " bin: "  << nomValue
	      << " binVolume:  "  << binVolumeDS
	      << " paramValue:  "  << paramVal
	      << " nomValue:  "  << nomValue
	      << " subTotal:  "  << value
	      << std::endl;
    */

  }

  return value;

}



////////////////////////////////////////////////////////////////////////////////
/// Return sampling hint for making curves of (projections) of this function
/// as the recursive division strategy of RooCurve cannot deal efficiently
/// with the vertical lines that occur in a non-interpolated histogram

std::list<Double_t>* ParamHistFunc::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, 
						Double_t xhi) const
{
  // copied and edited from RooHistFunc
  RooAbsLValue* lvarg = &obs;

  // Retrieve position of all bin boundaries
  const RooAbsBinning* binning = lvarg->getBinningPtr(0) ;
  Double_t* boundaries = binning->array() ;

  std::list<Double_t>* hint = new std::list<Double_t> ;

  // Widen range slighty
  xlo = xlo - 0.01*(xhi-xlo) ;
  xhi = xhi + 0.01*(xhi-xlo) ;

  Double_t delta = (xhi-xlo)*1e-8 ;
 
  // Construct array with pairs of points positioned epsilon to the left and
  // right of the bin boundaries
  for (Int_t i=0 ; i<binning->numBoundaries() ; i++) {
    if (boundaries[i]>=xlo && boundaries[i]<=xhi) {
      hint->push_back(boundaries[i]-delta) ;
      hint->push_back(boundaries[i]+delta) ;
    }
  }
  return hint ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return sampling hint for making curves of (projections) of this function
/// as the recursive division strategy of RooCurve cannot deal efficiently
/// with the vertical lines that occur in a non-interpolated histogram

std::list<Double_t>* ParamHistFunc::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, 
						  Double_t xhi) const 
{
  // copied and edited from RooHistFunc
  RooAbsLValue* lvarg = &obs;

  // Retrieve position of all bin boundaries
  const RooAbsBinning* binning = lvarg->getBinningPtr(0) ;
  Double_t* boundaries = binning->array() ;

  std::list<Double_t>* hint = new std::list<Double_t> ;

  // Construct array with pairs of points positioned epsilon to the left and
  // right of the bin boundaries
  for (Int_t i=0 ; i<binning->numBoundaries() ; i++) {
    if (boundaries[i]>=xlo && boundaries[i]<=xhi) {
      hint->push_back(boundaries[i]) ;
    }
  }

  return hint ;
}
