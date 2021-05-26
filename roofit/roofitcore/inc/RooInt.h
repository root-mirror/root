/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooInt.h,v 1.6 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_INT
#define ROO_INT

#include "Rtypes.h"
#include "TNamed.h"

class RooInt : public TNamed {
public:

  RooInt() : _value(0) {} ;
  RooInt(Int_t value) : TNamed(), _value(value) {} ;
  RooInt(const RooInt& other) : TNamed(other), _value(other._value) {}
  virtual ~RooInt() {} ;

  // Double_t cast operator
  inline operator Int_t() const { return _value ; }
  RooInt& operator=(Int_t value) { _value = value ; return *this ; }

  // Sorting interface ;
  Int_t Compare(const TObject* other) const ;
  virtual Bool_t IsSortable() const { return kTRUE ; }

protected:

  Int_t _value ; // Payload
  ClassDef(RooInt,1) // Container class for Int_t
};

#endif
