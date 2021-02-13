/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooAICRegistry.h,v 1.11 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_AIC_REGISTRY
#define ROO_AIC_REGISTRY

#include <vector>
#include "Rtypes.h"

class RooArgSet;

typedef RooArgSet* pRooArgSet ;

class RooAICRegistry {

public:
  RooAICRegistry(UInt_t size = 10) ;
  RooAICRegistry(const RooAICRegistry& other) ;
  virtual ~RooAICRegistry() ;

  Int_t store(const std::vector<Int_t>& codeList, RooArgSet* set1 = 0, RooArgSet* set2 = 0,
              RooArgSet* set3 = 0, RooArgSet* set4 = 0);
  const std::vector<Int_t>& retrieve(Int_t masterCode) const ;
  const std::vector<Int_t>&  retrieve(Int_t masterCode, pRooArgSet& set1) const ;
  const std::vector<Int_t>&  retrieve(Int_t masterCode, pRooArgSet& set1, pRooArgSet& set2) const ;
  const std::vector<Int_t>&  retrieve(Int_t masterCode, pRooArgSet& set1,
                                      pRooArgSet& set2, pRooArgSet& set3, pRooArgSet& set4) const ;

protected:

  std::vector<std::vector<Int_t> > _clArr;       //! Array of array of code lists
  std::vector<pRooArgSet> _asArr1;  //! Array of 1st RooArgSet pointers
  std::vector<pRooArgSet> _asArr2;  //! Array of 2nd RooArgSet pointers
  std::vector<pRooArgSet> _asArr3;  //! Array of 3rd RooArgSet pointers
  std::vector<pRooArgSet> _asArr4;  //! Array of 4th RooArgSet pointers

  ClassDef(RooAICRegistry,2) // Registry for analytical integration codes
} ;

#endif 
