/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooIntegratorBinding.h,v 1.4 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_INTEGRATOR_BINDING
#define ROO_INTEGRATOR_BINDING

#include "RooAbsFunc.h"
#include "RooAbsIntegrator.h"

class RooIntegratorBinding : public RooAbsFunc {
public:
  RooIntegratorBinding(RooAbsIntegrator& integrator) :
    RooAbsFunc(integrator.integrand()->getDimension()-1), _integrator(&integrator) {} ;
  virtual ~RooIntegratorBinding() {} ;

  inline virtual Double_t operator()(const Double_t xvector[]) const { _ncall++ ; return _integrator->integral(xvector) ; }
  inline virtual Double_t getMinLimit(UInt_t index) const { return _integrator->integrand()->getMinLimit(index+1); }
  inline virtual Double_t getMaxLimit(UInt_t index) const { return _integrator->integrand()->getMaxLimit(index+1); }

protected:
  RooAbsIntegrator* _integrator ;  // Numeric integrator


  ClassDef(RooIntegratorBinding,0) // Function binding representing output of numeric integrator
};

#endif

