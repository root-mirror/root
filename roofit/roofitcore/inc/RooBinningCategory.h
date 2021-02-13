/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
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
#ifndef ROO_BINNING_CATEGORY
#define ROO_BINNING_CATEGORY

#include "RooAbsCategory.h"
#include "RooTemplateProxy.h"
#include "TString.h"

class RooBinningCategory : public RooAbsCategory {

public:
  // Constructors etc.
  inline RooBinningCategory() { }
  RooBinningCategory(const char *name, const char *title, RooAbsRealLValue& inputVar, const char* binningName=0, const char* catTypeName=0);
  RooBinningCategory(const RooBinningCategory& other, const char *name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooBinningCategory(*this, newname); }
  virtual ~RooBinningCategory();

  // Printing interface (human readable)
  virtual void printMultiline(std::ostream& os, Int_t content, Bool_t verbose=kFALSE, TString indent="") const ;

protected:
  
  void initialize(const char* catTypeName=0) ;

  RooTemplateProxy<RooAbsRealLValue> _inputVar; // Input variable that is mapped
  TString _bname ;         // Name of the binning specification to be used to perform the mapping

  virtual value_type evaluate() const;
  /// The shape of this category does not need to be recomputed, as it creates states on the fly.
  void recomputeShape() { }

  ClassDef(RooBinningCategory,1) // RealVar-to-Category function defined by bin boundaries on input var
};

#endif
