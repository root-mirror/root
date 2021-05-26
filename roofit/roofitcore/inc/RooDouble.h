/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooDouble.h,v 1.8 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_DOUBLE
#define ROO_DOUBLE

#include "Rtypes.h"
#include "TNamed.h"

class RooDouble : public TNamed {
public:

  RooDouble() : _value(0) {
    // Default constructor
  } ;
  RooDouble(Double_t value) ;
  RooDouble(const RooDouble& other) : TNamed(other), _value(other._value) {}
  virtual ~RooDouble() {
    // Destructor
  } ;

  // Double_t cast operator
  inline operator Double_t() const {
    // Return value of contained double
    return _value ;
  }
  RooDouble& operator=(Double_t value) {
    // Return true if contained double equals value
    _value = value ; return *this ;
  }

  // Sorting interface ;
  Int_t Compare(const TObject* other) const ;
  virtual Bool_t IsSortable() const {
    // We are a sortable object
    return kTRUE ;
  }

protected:

  Double_t _value ; // Value payload
  ClassDef(RooDouble,1) // Container class for Double_t
};

#endif
