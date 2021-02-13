
/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Author:                                                                   *
 * Tristan du Pree, Nikhef, Amsterdam, tdupree@nikhef.nl                     *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/** \class RooMultiBinomial
    \ingroup Roofit

RooMultiBinomial is an efficiency function which makes all combinations of
efficiencies given as input different efficiency functions for different categories.

Given a dataset with a category C that determines if a given
event is accepted (1) or rejected (0) for the efficiency to be measured,
this class evaluates as F if C is 'accept' and as (1-F) if
C is 'reject'. Values of F below 0 and above 1 are clipped.
F may have an arbitrary number of dependents and parameters

The combination only 'reject' can be chosen to be visible or not visible
(and hence this efficiency is then equal to zero).
**/

#include "RooFit.h"

#include "RooMultiBinomial.h"
#include "RooStreamParser.h"
#include "RooArgList.h"
#include "RooAbsCategory.h"
#include "RooMsgService.h"
#include <string>
#include <vector>

using namespace std ;

ClassImp(RooMultiBinomial);

////////////////////////////////////////////////////////////////////////////////
/// Construct the efficiency functions from a list of efficiency functions
/// and a list of categories cat with two states (0,1) that indicate if a given
/// event should be counted as rejected or accepted respectively

RooMultiBinomial::RooMultiBinomial(const char *name, const char *title,
               const RooArgList& effFuncList,
               const RooArgList& catList,
               Bool_t ignoreNonVisible) :
  RooAbsReal(name,title),
  _catList("catList","list of cats", this),
  _effFuncList("effFuncList","list of eff funcs",this),
  _ignoreNonVisible(ignoreNonVisible)
{
  _catList.add(catList);
  _effFuncList.add(effFuncList);

  if (_catList.getSize() != effFuncList.getSize()) {
    coutE(InputArguments) << "RooMultiBinomial::ctor(" << GetName() << ") ERROR: Wrong input, should have equal number of categories and efficiencies." << endl;
    throw string("RooMultiBinomial::ctor() ERROR: Wrong input, should have equal number of categories and efficiencies") ;
  }

}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooMultiBinomial::RooMultiBinomial(const RooMultiBinomial& other, const char* name) :
  RooAbsReal(other, name),
  _catList("catList",this,other._catList),
  _effFuncList("effFuncList",this,other._effFuncList),
  _ignoreNonVisible(other._ignoreNonVisible)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooMultiBinomial::~RooMultiBinomial()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the raw value of the function which is the effFunc
/// value if cat==1 and it is (1-effFunc) if cat==0

Double_t RooMultiBinomial::evaluate() const
{
  Int_t effFuncListSize = _effFuncList.getSize();

  // Get efficiency function for category i

  vector<Double_t> effFuncVal(effFuncListSize);
  for (int i=0; i<effFuncListSize; ++i) {
    effFuncVal[i] = ((RooAbsReal&)_effFuncList[i]).getVal() ;
  }

  // Truncate efficiency functions in range 0.0-1.0

  for (int i=0; i<effFuncListSize; ++i) {
    if (effFuncVal[i]>1) {
      coutW(Eval) << "WARNING: Efficiency >1 (equal to " << effFuncVal[i]
        << " ), for i = " << i << "...TRUNCATED" << endl;
      effFuncVal[i] = 1.0 ;
    } else if (effFuncVal[i]<0) {
      effFuncVal[i] = 0.0 ;
      coutW(Eval) << "WARNING: Efficiency <0 (equal to " << effFuncVal[i]
        << " ), for i = " << i << "...TRUNCATED" << endl;
    }
  }

  vector<Double_t> effValue(effFuncListSize);
  Bool_t notVisible = true;

  // Calculate efficiency per accept/reject decision

  for (int i=0; i<effFuncListSize; ++i) {
    if ( ((RooAbsCategory&)_catList[i]).getCurrentIndex() == 1) {
      // Accept case
      effValue[i] = effFuncVal[i] ;
      notVisible = false;
    } else if ( ((RooAbsCategory&)_catList[i]).getCurrentIndex() == 0){
      // Reject case
      effValue[i] = 1 - effFuncVal[i] ;
    } else {
      coutW(Eval) << "WARNING: WRONG CATEGORY NAMES GIVEN!, label = " << ((RooAbsCategory&)_catList[i]).getCurrentIndex() << endl;
      effValue[i] = 0;
    }
  }

  Double_t _effVal = 1.;

  // Calculate efficiency for combination of accept/reject categories
  // put equal to zero if combination of only zeros AND chosen to be invisible

  for (int i=0; i<effFuncListSize; ++i) {
    _effVal=_effVal*effValue[i];
    if (notVisible && _ignoreNonVisible){
      _effVal=0;
    }
  }

  return _effVal;

}
