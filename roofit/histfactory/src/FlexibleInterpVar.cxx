// @(#)root/roostats:$Id:  cranmer $
// Author: Kyle Cranmer, Akira Shibata
// Author: Giovanni Petrucciani (UCSD) (log-interpolation)
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////

/** \class RooStats::HistFactory::FlexibleInterpVar
 *  \ingroup HistFactory 
 */

#include "RooFit.h"

#include "Riostream.h"
#include <math.h>
#include "TMath.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooMsgService.h"
#include "RooTrace.h"

#include "RooStats/HistFactory/FlexibleInterpVar.h"

using namespace std;

ClassImp(RooStats::HistFactory::FlexibleInterpVar);

using namespace RooStats;
using namespace HistFactory;

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

FlexibleInterpVar::FlexibleInterpVar()
{
  _paramIter = _paramList.createIterator() ;
  _nominal = 0;
  _interpBoundary=1.;
  _logInit = kFALSE ;
  TRACE_CREATE
}


////////////////////////////////////////////////////////////////////////////////

FlexibleInterpVar::FlexibleInterpVar(const char* name, const char* title, 
		       const RooArgList& paramList, 
		       Double_t argNominal, vector<double> lowVec, vector<double> highVec) :
  RooAbsReal(name, title),
  _paramList("paramList","List of paramficients",this),
  _nominal(argNominal), _low(lowVec), _high(highVec), _interpBoundary(1.)
{
  _logInit = kFALSE ;
  _paramIter = _paramList.createIterator() ;


  TIterator* paramIter = paramList.createIterator() ;
  RooAbsArg* param ;
  while((param = (RooAbsArg*)paramIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(param)) {
      coutE(InputArguments) << "FlexibleInterpVar::ctor(" << GetName() << ") ERROR: paramficient " << param->GetName() 
			    << " is not of type RooAbsReal" << endl ;
      R__ASSERT(0) ;
    }
    _paramList.add(*param) ;
    _interpCode.push_back(0); // default code
  }
  if (int(_low.size() ) != _paramList.getSize() || _low.size() != _high.size()) {
     coutE(InputArguments) << "FlexibleInterpVar::ctor(" << GetName() << ") invalid input low/high vectors " << endl;
     R__ASSERT(int(_low.size() ) == _paramList.getSize());
     R__ASSERT(_low.size() == _high.size());
  }

  delete paramIter ;
  TRACE_CREATE

}

////////////////////////////////////////////////////////////////////////////////

FlexibleInterpVar::FlexibleInterpVar(const char* name, const char* title, 
		       const RooArgList& paramList, 
		       double argNominal, const RooArgList& lowList, const RooArgList& highList) :
  RooAbsReal(name, title),
  _paramList("paramList","List of paramficients",this),
  _nominal(argNominal), _interpBoundary(1.)
{
  RooFIter lowIter = lowList.fwdIterator() ;
  RooAbsReal* val ; 
  while ((val = (RooAbsReal*) lowIter.next())) {
    _low.push_back(val->getVal()) ;
  }

  RooFIter highIter = highList.fwdIterator() ;
  while ((val = (RooAbsReal*) highIter.next())) {
    _high.push_back(val->getVal()) ;
  }
  
  
  _logInit = kFALSE ;
  _paramIter = _paramList.createIterator() ;


  TIterator* paramIter = paramList.createIterator() ;
  RooAbsArg* param ;
  while((param = (RooAbsArg*)paramIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(param)) {
      coutE(InputArguments) << "FlexibleInterpVar::ctor(" << GetName() << ") ERROR: paramficient " << param->GetName() 
			    << " is not of type RooAbsReal" << endl ;
      R__ASSERT(0) ;
    }
    _paramList.add(*param) ;
    _interpCode.push_back(0); // default code
  }
  if (int(_low.size() ) != _paramList.getSize() || _low.size() != _high.size()) {
     coutE(InputArguments) << "FlexibleInterpVar::ctor(" << GetName() << ") invalid input low/high lists " << endl;
     R__ASSERT(int(_low.size() ) == _paramList.getSize());
     R__ASSERT(_low.size() == _high.size());
  }

  delete paramIter ;
  TRACE_CREATE

}




////////////////////////////////////////////////////////////////////////////////

FlexibleInterpVar::FlexibleInterpVar(const char* name, const char* title, 
				     const RooArgList& paramList, 
				     double argNominal, vector<double> lowVec, vector<double> highVec,
				     vector<int> code) :
  RooAbsReal(name, title),
  _paramList("paramList","List of paramficients",this),
  _nominal(argNominal), _low(lowVec), _high(highVec), _interpCode(code), _interpBoundary(1.)
{
  _logInit = kFALSE ;
  _paramIter = _paramList.createIterator() ;


  TIterator* paramIter = paramList.createIterator() ;
  RooAbsArg* param ;
  while((param = (RooAbsArg*)paramIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(param)) {
      coutE(InputArguments) << "FlexibleInterpVar::ctor(" << GetName() << ") ERROR: paramficient " << param->GetName()
			    << " is not of type RooAbsReal" << endl ;
      // use R__ASSERT which remains also in release mode
      R__ASSERT(0) ;
    }
    _paramList.add(*param) ;
  } 
  if (int(_low.size() ) != _paramList.getSize() || _low.size() != _high.size() || _low.size() != _interpCode.size()) {
     coutE(InputArguments) << "FlexibleInterpVar::ctor(" << GetName() << ") invalid input vectors " << endl;
     R__ASSERT(int(_low.size() ) == _paramList.getSize());
     R__ASSERT(_low.size() == _high.size());
     R__ASSERT(_low.size() == _interpCode.size());
  }
  delete paramIter ;
  TRACE_CREATE

}

////////////////////////////////////////////////////////////////////////////////
/// Constructor of flat polynomial function

FlexibleInterpVar::FlexibleInterpVar(const char* name, const char* title) :
  RooAbsReal(name, title),
  _paramList("paramList","List of coefficients",this), 
  _nominal(0), _interpBoundary(1.)
{
  _logInit = kFALSE ;
  _paramIter = _paramList.createIterator() ;
  TRACE_CREATE
}

////////////////////////////////////////////////////////////////////////////////

FlexibleInterpVar::FlexibleInterpVar(const FlexibleInterpVar& other, const char* name) :
  RooAbsReal(other, name), 
  _paramList("paramList",this,other._paramList),
  _nominal(other._nominal), _low(other._low), _high(other._high), _interpCode(other._interpCode), _interpBoundary(other._interpBoundary)
  
{
  // Copy constructor
  _logInit = kFALSE ;
  _paramIter = _paramList.createIterator() ;
  TRACE_CREATE
  
}


////////////////////////////////////////////////////////////////////////////////
/// Destructor

FlexibleInterpVar::~FlexibleInterpVar() 
{
  delete _paramIter ;
  TRACE_DESTROY
}


////////////////////////////////////////////////////////////////////////////////

void FlexibleInterpVar::setInterpCode(RooAbsReal& param, int code){
  int index = _paramList.index(&param);
  if(index<0){
      coutE(InputArguments) << "FlexibleInterpVar::setInterpCode ERROR:  " << param.GetName() 
			    << " is not in list" << endl ;
  } else {
      coutW(InputArguments) << "FlexibleInterpVar::setInterpCode :  " << param.GetName() 
			    << " is now " << code << endl ;
    _interpCode.at(index) = code;
  }
  // GHL: Adding suggestion by Swagato:
  _logInit = kFALSE ;
  setValueDirty();
}

////////////////////////////////////////////////////////////////////////////////

void FlexibleInterpVar::setAllInterpCodes(int code){
  for(unsigned int i=0; i<_interpCode.size(); ++i){
    _interpCode.at(i) = code;
  }
  // GHL: Adding suggestion by Swagato:
  _logInit = kFALSE ;
  setValueDirty();

}

////////////////////////////////////////////////////////////////////////////////

void FlexibleInterpVar::setNominal(Double_t newNominal){
  coutW(InputArguments) << "FlexibleInterpVar::setNominal : nominal is now " << newNominal << endl ;
  _nominal = newNominal;

  _logInit = kFALSE ;

  setValueDirty();
}

////////////////////////////////////////////////////////////////////////////////

void FlexibleInterpVar::setLow(RooAbsReal& param, Double_t newLow){
  int index = _paramList.index(&param);
  if(index<0){
      coutE(InputArguments) << "FlexibleInterpVar::setLow ERROR:  " << param.GetName() 
			    << " is not in list" << endl ;
  } else {
      coutW(InputArguments) << "FlexibleInterpVar::setLow :  " << param.GetName() 
			    << " is now " << newLow << endl ;
    _low.at(index) = newLow;
  }

  _logInit = kFALSE ;

  setValueDirty();
}

////////////////////////////////////////////////////////////////////////////////

void FlexibleInterpVar::setHigh(RooAbsReal& param, Double_t newHigh){
  int index = _paramList.index(&param);
  if(index<0){
      coutE(InputArguments) << "FlexibleInterpVar::setHigh ERROR:  " << param.GetName() 
			    << " is not in list" << endl ;
  } else {
      coutW(InputArguments) << "FlexibleInterpVar::setHigh :  " << param.GetName() 
			    << " is now " << newHigh << endl ;
    _high.at(index) = newHigh;
  }

  _logInit = kFALSE ;
  setValueDirty();
}

////////////////////////////////////////////////////////////////////////////////

void FlexibleInterpVar::printAllInterpCodes(){
  for(unsigned int i=0; i<_interpCode.size(); ++i){
    coutI(InputArguments) <<"interp code for " << _paramList.at(i)->GetName() << " = " << _interpCode.at(i) <<endl;
    // GHL: Adding suggestion by Swagato:
    if( _low.at(i) <= 0.001 ) coutE(InputArguments) << GetName() << ", " << _paramList.at(i)->GetName() << ": low value = " << _low.at(i) << endl;
    if( _high.at(i) <= 0.001 ) coutE(InputArguments) << GetName() << ", " << _paramList.at(i)->GetName() << ": high value = " << _high.at(i) << endl;
  }

}

////////////////////////////////////////////////////////////////////////////////

double FlexibleInterpVar::PolyInterpValue(int i, double x) const { 
   // code for polynomial interpolation used when interpCode=4

   double boundary = _interpBoundary;

   double x0 = boundary;


   // cache the polynomial coefficient values
   // which do not depend on x but on the boundaries values
   if (!_logInit) {
      
      _logInit=kTRUE ;

      unsigned int n = _low.size(); 
      assert(n == _high.size() );
      
      _polCoeff.resize(n*6) ;
      
      for (unsigned int j = 0; j < n ; j++) {
         
         // location of the 6 coefficient for the j-th variable
         double * coeff = &_polCoeff[j * 6]; 
         
         // GHL: Swagato's suggestions
         double pow_up       =  std::pow(_high[j]/_nominal, x0);
         double pow_down     =  std::pow(_low[j]/_nominal,  x0);
         double logHi        =  std::log(_high[j]) ; 
         double logLo        =  std::log(_low[j] );
         double pow_up_log   = _high[j] <= 0.0 ? 0.0 : pow_up      * logHi;
         double pow_down_log = _low[j] <= 0.0 ? 0.0 : -pow_down    * logLo;
         double pow_up_log2  = _high[j] <= 0.0 ? 0.0 : pow_up_log  * logHi;
         double pow_down_log2= _low[j] <= 0.0 ? 0.0 : -pow_down_log* logLo;

         double S0 = (pow_up+pow_down)/2;
         double A0 = (pow_up-pow_down)/2;
         double S1 = (pow_up_log+pow_down_log)/2;
         double A1 = (pow_up_log-pow_down_log)/2;
         double S2 = (pow_up_log2+pow_down_log2)/2;
         double A2 = (pow_up_log2-pow_down_log2)/2;
         
         //fcns+der+2nd_der are eq at bd
         
         // cache  coefficient of the polynomial 
         coeff[0] = 1./(8*x0)        *(      15*A0 -  7*x0*S1 + x0*x0*A2);
         coeff[1] = 1./(8*x0*x0)     *(-24 + 24*S0 -  9*x0*A1 + x0*x0*S2);
         coeff[2] = 1./(4*pow(x0, 3))*(    -  5*A0 +  5*x0*S1 - x0*x0*A2);
         coeff[3] = 1./(4*pow(x0, 4))*( 12 - 12*S0 +  7*x0*A1 - x0*x0*S2);
         coeff[4] = 1./(8*pow(x0, 5))*(    +  3*A0 -  3*x0*S1 + x0*x0*A2);
         coeff[5] = 1./(8*pow(x0, 6))*( -8 +  8*S0 -  5*x0*A1 + x0*x0*S2);
         
      }
      
   }
   
   // GHL: Swagato's suggestions
   // if( _low[i] == 0 ) _low[i] = 0.0001;
   // if( _high[i] == 0 ) _high[i] = 0.0001;
   
   // get pointer to location of coefficients in the vector

   assert(int(_polCoeff.size()) > i ); 
   const double * coefficients = &_polCoeff.front() + 6*i;  
   
   double a = coefficients[0];
   double b = coefficients[1];
   double c = coefficients[2];
   double d = coefficients[3];
   double e = coefficients[4];
   double f = coefficients[5];
   

   // evaluate the 6-th degree polynomial using Horner's method
   double value = 1. + x * (a + x * ( b + x * ( c + x * ( d + x * ( e + x * f ) ) ) ) );
   return value; 
}

////////////////////////////////////////////////////////////////////////////////
/// Const getters

const RooListProxy& FlexibleInterpVar::variables() const { return _paramList; }
double FlexibleInterpVar::nominal() const { return _nominal; }        
const std::vector<double>& FlexibleInterpVar::low() const { return _low; }
const std::vector<double>& FlexibleInterpVar::high() const { return _high; }

////////////////////////////////////////////////////////////////////////////////
/// Calculate and return value of polynomial

Double_t FlexibleInterpVar::evaluate() const 
{
  Double_t total(_nominal) ;
  _paramIter->Reset() ;

  RooAbsReal* param ;
  //const RooArgSet* nset = _paramList.nset() ;
  int i=0;

  // TString name = GetName();
  // if (name == TString("ZHWW_ll12_vzll_epsilon") )
  //    //    std::cout << "evaluate flexible interp var - init flag is " << _logInit << std::endl;

  while((param=(RooAbsReal*)_paramIter->Next())) {
    //    param->Print("v");


    Int_t icode = _interpCode[i] ;

    switch(icode) {

    case 0: {
      // piece-wise linear
      if(param->getVal()>0)
	total +=  param->getVal()*(_high[i] - _nominal );
      else
	total += param->getVal()*(_nominal - _low[i]);
      break ;
    }
    case 1: {
      // pice-wise log
      if(param->getVal()>=0)
	total *= pow(_high[i]/_nominal, +param->getVal());
      else
	total *= pow(_low[i]/_nominal,  -param->getVal());
      break ;
    }
    case 2: {
      // parabolic with linear
      double a = 0.5*(_high[i]+_low[i])-_nominal;
      double b = 0.5*(_high[i]-_low[i]);
      double c = 0;
      if(param->getVal()>1 ){
	total += (2*a+b)*(param->getVal()-1)+_high[i]-_nominal;
      } else if(param->getVal()<-1 ) {
	total += -1*(2*a-b)*(param->getVal()+1)+_low[i]-_nominal;
      } else {
	total +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
      }
      break ;
    }
    case 3: {
      //parabolic version of log-normal
      double a = 0.5*(_high[i]+_low[i])-_nominal;
      double b = 0.5*(_high[i]-_low[i]);
      double c = 0;
      if(param->getVal()>1 ){
	total += (2*a+b)*(param->getVal()-1)+_high[i]-_nominal;
      } else if(param->getVal()<-1 ) {
	total += -1*(2*a-b)*(param->getVal()+1)+_low[i]-_nominal;
      } else {
	total +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
      }
      break ;
    }

    case 4: {
      double boundary = _interpBoundary;
      double x = param->getVal(); 
      //std::cout << icode << " param " << param->GetName() << "  " << param->getVal() << " boundary " << boundary << std::endl;

      if(x >= boundary)
      {
         total *= std::pow(_high[i]/_nominal, +param->getVal());
      }
      else if (x <= -boundary)
      {
         total *= std::pow(_low[i]/_nominal, -param->getVal());
      }
      else if (x != 0)
      {
         total *= PolyInterpValue(i, x);
      }
      break ;
    }
    default: {
      coutE(InputArguments) << "FlexibleInterpVar::evaluate ERROR:  " << param->GetName() 
			    << " with unknown interpolation code" << endl ;
    }
    }
    ++i;
  }

  if(total<=0) {
     total= TMath::Limits<double>::Min();
  }    

  return total;
}

void FlexibleInterpVar::printMultiline(ostream& os, Int_t contents, 
				       Bool_t verbose, TString indent) const
{
  RooAbsReal::printMultiline(os,contents,verbose,indent);
  os << indent << "--- FlexibleInterpVar ---" << endl;
  printFlexibleInterpVars(os);
}

void FlexibleInterpVar::printFlexibleInterpVars(ostream& os) const
{
  _paramIter->Reset();
  for (int i=0;i<(int)_low.size();i++) {
    RooAbsReal* param=(RooAbsReal*)_paramIter->Next();
    os << setw(36) << param->GetName()<<": "<<setw(7) << _low[i]<<"  "<<setw(7) << _high[i]
       <<endl;
  }
}


