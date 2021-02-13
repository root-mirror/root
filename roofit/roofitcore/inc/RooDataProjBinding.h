/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooDataProjBinding.h,v 1.6 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_DATA_PROJ_BINDING
#define ROO_DATA_PROJ_BINDING

#include "RooRealBinding.h"
class RooAbsReal ;
class RooAbsData ;
class RooSuperCategory ;
class Roo1DTable ;

class RooDataProjBinding : public RooRealBinding {
public:
  RooDataProjBinding(const RooAbsReal &real, const RooAbsData& data, const RooArgSet &vars, const RooArgSet* normSet=0) ;
  virtual ~RooDataProjBinding() ;

  virtual Double_t operator()(const Double_t xvector[]) const;

  RooSpan<const double> getValues(std::vector<RooSpan<const double>> coordinates) const;

protected:

  mutable Bool_t _first   ;  // Bit indicating if operator() has been called yet
  const RooAbsReal* _real ;  // Real function to be projected
  const RooAbsData* _data ;  // Dataset used for projection
  const RooArgSet*  _nset ;  // Normalization set for real function

  RooSuperCategory* _superCat ;  // Supercategory constructed from _data's category variables
  Roo1DTable* _catTable ;        // Supercategory table generated from _data
  mutable std::unique_ptr<std::vector<double>> _batchBuffer; //! Storage for handing out spans.

  ClassDef(RooDataProjBinding,0) // RealFunc/Dataset binding for data projection of a real function
};

#endif

