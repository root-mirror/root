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
\file RooGenProdProj.cxx
\class RooGenProdProj
\ingroup Roofitcore


RooGenProdProj is an auxiliary class for RooProdPdf that calculates
a general normalised projection of a product of non-factorising PDFs, e.g.
\f[
 P_{x,xy} = \frac{\int ( P1 * P2 * \ldots) \mathrm{d}x}{\int ( P1 * P2 * \ldots ) \mathrm{d}x \mathrm{d}y}
\f]

Partial integrals, which factorise and can be calculated, are calculated
analytically. Remaining non-factorising observables are integrated numerically.
**/


#include "RooFit.h"

#include "Riostream.h"
#include <math.h>

#include "RooGenProdProj.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooErrorHandler.h"
#include "RooProduct.h"

using namespace std;

ClassImp(RooGenProdProj);



////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooGenProdProj::RooGenProdProj() :
  _compSetOwnedN(0), 
  _compSetOwnedD(0),
  _haveD(kFALSE)
{
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor for a normalization projection of the product of p.d.f.s _prodSet
/// integrated over _intSet in range isetRangeName while normalized over _normSet

RooGenProdProj::RooGenProdProj(const char *name, const char *title, const RooArgSet& _prodSet, const RooArgSet& _intSet, 
			       const RooArgSet& _normSet, const char* isetRangeName, const char* normRangeName, Bool_t doFactorize) :
  RooAbsReal(name, title),
  _compSetOwnedN(0), 
  _compSetOwnedD(0),
  _compSetN("compSetN","Set of integral components owned by numerator",this,kFALSE),
  _compSetD("compSetD","Set of integral components owned by denominator",this,kFALSE),
  _intList("intList","List of integrals",this,kTRUE),
  _haveD(kFALSE)
{
  // Set expensive object cache to that of first item in prodSet
  setExpensiveObjectCache(_prodSet.first()->expensiveObjectCache()) ;

  // Create owners of components created in ctor
  _compSetOwnedN = new RooArgSet ;
  _compSetOwnedD = new RooArgSet ;

  RooAbsReal* numerator = makeIntegral("numerator",_prodSet,_intSet,*_compSetOwnedN,isetRangeName,doFactorize) ;
  RooAbsReal* denominator = makeIntegral("denominator",_prodSet,_normSet,*_compSetOwnedD,normRangeName,doFactorize) ;

//   cout << "RooGenProdPdf::ctor(" << GetName() << ") numerator = " << numerator->GetName() << endl ;
//   numerator->printComponentTree() ;
//   cout << "RooGenProdPdf::ctor(" << GetName() << ") denominator = " << denominator->GetName() << endl ;
//   denominator->printComponentTree() ;

  // Copy all components in (non-owning) set proxy
  _compSetN.add(*_compSetOwnedN) ;
  _compSetD.add(*_compSetOwnedD) ;
  
  _intList.add(*numerator) ;
  if (denominator) {
    _intList.add(*denominator) ;
    _haveD = kTRUE ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooGenProdProj::RooGenProdProj(const RooGenProdProj& other, const char* name) :
  RooAbsReal(other, name), 
  _compSetOwnedN(0), 
  _compSetOwnedD(0),
  _compSetN("compSetN","Set of integral components owned by numerator",this),
  _compSetD("compSetD","Set of integral components owned by denominator",this),
  _intList("intList","List of integrals",this)
{
  // Explicitly remove all server links at this point
  TIterator* iter = serverIterator() ;
  RooAbsArg* server ;
  while((server=(RooAbsArg*)iter->Next())) {
    removeServer(*server,kTRUE) ;
  }
  delete iter ;

  // Copy constructor
  _compSetOwnedN = (RooArgSet*) other._compSetN.snapshot() ;
  _compSetN.add(*_compSetOwnedN) ;

  _compSetOwnedD = (RooArgSet*) other._compSetD.snapshot() ;
  _compSetD.add(*_compSetOwnedD) ;

  RooAbsArg* arg ;
  TIterator* nIter = _compSetOwnedN->createIterator() ;  
  while((arg=(RooAbsArg*)nIter->Next())) {
//     cout << "ownedN elem " << arg->GetName() << "(" << arg << ")" << endl ;
    arg->setOperMode(_operMode) ;
  }
  delete nIter ;
  TIterator* dIter = _compSetOwnedD->createIterator() ;
  while((arg=(RooAbsArg*)dIter->Next())) {
//     cout << "ownedD elem " << arg->GetName() << "(" << arg << ")" << endl ;
    arg->setOperMode(_operMode) ;
  }
  delete dIter ;

  // Fill _intList
  _haveD = other._haveD ;
  _intList.add(*_compSetN.find(other._intList.at(0)->GetName())) ;
  if (other._haveD) {
    _intList.add(*_compSetD.find(other._intList.at(1)->GetName())) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooGenProdProj::~RooGenProdProj()
{
  if (_compSetOwnedN) delete _compSetOwnedN ;
  if (_compSetOwnedD) delete _compSetOwnedD ;
}



////////////////////////////////////////////////////////////////////////////////
/// Utility function to create integral for product over certain observables.
/// \param[in] name Name of integral to be created.
/// \param[in] compSet All components of the product.
/// \param[in] intSet Observables to be integrated.
/// \param[out] saveSet All component objects needed to represent the product integral are added as owned members to saveSet.
/// \note The set owns new components that are created for the integral.
/// \param[in] isetRangeName Integral range.
/// \param[in] doFactorize
///
/// \return A RooAbsReal object representing the requested integral. The object is owned by `saveSet`.
///
/// The integration is factorized into components as much as possible and done analytically as far as possible.
RooAbsReal* RooGenProdProj::makeIntegral(const char* name, const RooArgSet& compSet, const RooArgSet& intSet, 
					 RooArgSet& saveSet, const char* isetRangeName, Bool_t doFactorize) 
{
  RooArgSet anaIntSet, numIntSet ;

  // First determine subset of observables in intSet that are factorizable
  for (const auto arg : intSet) {
    auto count = std::count_if(compSet.begin(), compSet.end(), [arg](const RooAbsArg* pdfAsArg){
      auto pdf = static_cast<const RooAbsPdf*>(pdfAsArg);
      return (pdf->dependsOn(*arg));
    });

    if (count==1) {
      anaIntSet.add(*arg) ;
    }    
  }

  // Determine which of the factorizable integrals can be done analytically
  RooArgSet prodSet ;
  numIntSet.add(intSet) ;
  
  for (const auto pdfAsArg : compSet) {
    auto pdf = static_cast<const RooAbsPdf*>(pdfAsArg);

    if (doFactorize && pdf->dependsOn(anaIntSet)) {
      RooArgSet anaSet ;
      Int_t code = pdf->getAnalyticalIntegralWN(anaIntSet,anaSet,0,isetRangeName) ;
      if (code!=0) {
        // Analytical integral, create integral object
        RooAbsReal* pai = pdf->createIntegral(anaSet,isetRangeName) ;
        pai->setOperMode(_operMode) ;

        // Add to integral to product
        prodSet.add(*pai) ;

        // Remove analytically integratable observables from numeric integration list
        numIntSet.remove(anaSet) ;

        // Declare ownership of integral
        saveSet.addOwned(*pai) ;
      } else {
        // Analytic integration of factorizable observable not possible, add straight pdf to product
        prodSet.add(*pdf) ;
      }      
    } else {
      // Non-factorizable observables, add straight pdf to product
      prodSet.add(*pdf) ;
    }
  }

  // Create product of (partial) analytical integrals
  TString prodName ;
  if (isetRangeName) {
    prodName = Form("%s_%s_Range[%s]",GetName(),name,isetRangeName) ;
  } else {
    prodName = Form("%s_%s",GetName(),name) ;
  }

  // Create clones of the elements in prodSet. These need to be cloned
  // because when caching optimisation lvl 2 is activated, pre-computed
  // values are side-loaded into the elements.
  // Those pre-cached values already contain normalisation constants, so
  // the integral comes out wrongly. Therefore, we create here nodes that
  // don't participate in any caching, which are used to compute integrals.
  RooArgSet prodSetClone;
  prodSet.snapshot(prodSetClone, false);
  saveSet.addOwned(prodSetClone);
  prodSetClone.releaseOwnership();

  RooProduct* prod = new RooProduct(prodName, "product", prodSetClone);
  prod->setExpensiveObjectCache(expensiveObjectCache()) ;
  prod->setOperMode(_operMode) ;

  // Declare ownership of product
  saveSet.addOwned(*prod);

  // Create integral performing remaining numeric integration over (partial) analytic product
  RooAbsReal* ret = prod->createIntegral(numIntSet,isetRangeName) ;
  ret->setOperMode(_operMode) ;
  saveSet.addOwned(*ret) ;

  // Caller owners returned master integral object
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate and return value of normalization projection

Double_t RooGenProdProj::evaluate() const 
{  
  Double_t nom = ((RooAbsReal*)_intList.at(0))->getVal() ;

  if (!_haveD) return nom ;

  Double_t den = ((RooAbsReal*)_intList.at(1))->getVal() ;

  //cout << "RooGenProdProj::eval(" << GetName() << ") nom = " << nom << " den = " << den << endl ;

  return nom / den ;
}



////////////////////////////////////////////////////////////////////////////////
/// Intercept cache mode operation changes and propagate them to the components

void RooGenProdProj::operModeHook() 
{
  // WVE use cache manager here!

  RooAbsArg* arg ;
  TIterator* nIter = _compSetOwnedN->createIterator() ;  
  while((arg=(RooAbsArg*)nIter->Next())) {
    arg->setOperMode(_operMode) ;
  }
  delete nIter ;

  TIterator* dIter = _compSetOwnedD->createIterator() ;
  while((arg=(RooAbsArg*)dIter->Next())) {
    arg->setOperMode(_operMode) ;
  }
  delete dIter ;

  _intList.at(0)->setOperMode(_operMode) ;
  if (_haveD) _intList.at(1)->setOperMode(Auto) ; // Denominator always stays in Auto mode (normalization integral)
}







