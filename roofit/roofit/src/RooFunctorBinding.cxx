/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, NIKHEF, verkerke@nikhef.nl                         *
 *                                                                           *
 * Copyright (c) 2000-2008, NIKHEF, Regents of the University of California  *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 *****************************************************************************/

/** \class RooFunctorBinding
    \ingroup Roofit

RooFunctorBinding makes math functions from ROOT usable in RooFit. It takes
a ROOT::Math::IBaseFunctionMultiDim, and binds the variables of this function to
the RooFit variables passed in the constructor.

Instances of function binding
classes are fully functional RooFit function objects with one exception:
if the bound function is *not* a standard TMath or MathMore function the
class cannot be persisted in a RooWorkspace without registering the function
pointer first using RooCFunction1Binding<T1,T2>::register().
**/

/** \class RooFunctorPdfBinding
    \ingroup Roofit
RooFunctorPdfBinding makes math functions from ROOT usable as PDFs in RooFit. It takes
a ROOT::Math::IBaseFunctionMultiDim, and binds the variables of this function to
the RooFit variables passed in the constructor.
When the PDF is evaluated, the bound function is evaluated, and also integrated numerically
to normalise it to unity over the range of its observables.

Instances of function binding
classes are fully functional RooFit function objects with one exception:
if the bound function is *not* a standard TMath or MathMore function the
class cannot be persisted in a RooWorkspace without registering the function
pointer first using RooCFunction1Binding<T1,T2>::register().
**/

#include "Riostream.h"
#include "RooFunctorBinding.h"

using namespace std ;

ClassImp(RooFunctorBinding);
ClassImp(RooFunctorPdfBinding);

////////////////////////////////////////////////////////////////////////////////
/// Create a RooFit function that makes `ftor` usable in RooFit.
/// \param name Name of the object.
/// \param title Title (e.g. for plotting)
/// \param ftor Functor instance to be evaluated.
/// \param v RooFit variables to be passed to the function.
RooFunctorBinding::RooFunctorBinding(const char *name, const char *title, const ROOT::Math::IBaseFunctionMultiDim& ftor, const RooArgList& v) :
  RooAbsReal(name,title),
  func(&ftor),
  vars("vars","vars",this)
{
  // Check that function dimension and number of variables match
  if (ftor.NDim()!=UInt_t(v.getSize())) {
    coutE(InputArguments) << "RooFunctorBinding::ctor(" << GetName() << ") ERROR number of provided variables (" << v.getSize()
           << ") does not match dimensionality of function (" << ftor.NDim() << ")" << endl ;
    throw string("RooFunctor::ctor ERROR") ;
  }
  x = new Double_t[func->NDim()] ;
  vars.add(v) ;
}

////////////////////////////////////////////////////////////////////////////////
RooFunctorBinding::RooFunctorBinding(const RooFunctorBinding& other, const char* name) :
  RooAbsReal(other,name),
  func(other.func),
  vars("vars",this,other.vars)
{
  // Copy constructor
  x = new Double_t[func->NDim()] ;
}

////////////////////////////////////////////////////////////////////////////////
void RooFunctorBinding::printArgs(ostream& os) const {
  // Print object arguments and name/address of function pointer
  os << "[ function=" << func << " " ;
  for (Int_t i=0 ; i<numProxies() ; i++) {
    RooAbsProxy* p = getProxy(i) ;
    if (!TString(p->name()).BeginsWith("!")) {
      p->print(os) ;
      os << " " ;
    }
  }
  os << "]" ;
}

////////////////////////////////////////////////////////////////////////////////
Double_t RooFunctorBinding::evaluate() const {
    // Return value of embedded function using value of referenced variable x
    for (int i=0 ; i<vars.getSize() ; i++) {
      x[i] = ((RooAbsReal*)vars.at(i))->getVal() ;
    }
    return (*func)(x) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a RooFit PDF that makes `ftor` usable as a PDF in RooFit.
/// \param name Name of the object.
/// \param title Title (e.g. for plotting)
/// \param ftor Functor instance to be evaluated and normalised.
/// \param v RooFit variables to be passed to the function.
RooFunctorPdfBinding::RooFunctorPdfBinding(const char *name, const char *title, const ROOT::Math::IBaseFunctionMultiDim& ftor, const RooArgList& v) :
  RooAbsPdf(name,title),
  func(&ftor),
  vars("vars","vars",this)
{
  // Check that function dimension and number of variables match
  if (ftor.NDim()!=UInt_t(v.getSize())) {
    coutE(InputArguments) << "RooFunctorPdfBinding::ctor(" << GetName() << ") ERROR number of provided variables (" << v.getSize()
           << ") does not match dimensionality of function (" << ftor.NDim() << ")" << endl ;
    throw string("RooFunctor::ctor ERROR") ;
  }
  x = new Double_t[func->NDim()] ;
  vars.add(v) ;
}

////////////////////////////////////////////////////////////////////////////////
RooFunctorPdfBinding::RooFunctorPdfBinding(const RooFunctorPdfBinding& other, const char* name) :
  RooAbsPdf(other,name),
  func(other.func),
  vars("vars",this,other.vars)
{
  // Copy constructor
  x = new Double_t[func->NDim()] ;
}

////////////////////////////////////////////////////////////////////////////////
void RooFunctorPdfBinding::printArgs(ostream& os) const {
  // Print object arguments and name/address of function pointer
  os << "[ function=" << func << " " ;
  for (Int_t i=0 ; i<numProxies() ; i++) {
    RooAbsProxy* p = getProxy(i) ;
    if (!TString(p->name()).BeginsWith("!")) {
      p->print(os) ;
      os << " " ;
    }
  }
  os << "]" ;
}

////////////////////////////////////////////////////////////////////////////////
Double_t RooFunctorPdfBinding::evaluate() const {
    // Return value of embedded function using value of referenced variable x
    for (int i=0 ; i<vars.getSize() ; i++) {
      x[i] = ((RooAbsReal*)vars.at(i))->getVal() ;
    }
    return (*func)(x) ;
  }

namespace RooFit {

  RooAbsReal* bindFunction(const char* name, const ROOT::Math::IBaseFunctionMultiDim& ftor,const RooArgList& vars) {
    return new RooFunctorBinding(name,name,ftor,vars) ;
  }

  RooAbsPdf*  bindPdf(const char* name, const ROOT::Math::IBaseFunctionMultiDim& ftor, const RooArgList& vars) {
    return new RooFunctorPdfBinding(name,name,ftor,vars) ;
  }

}
