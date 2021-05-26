/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooTObjWrap.h,v 1.7 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_TOBJ_WRAP
#define ROO_TOBJ_WRAP

#include "Rtypes.h"
#include "TNamed.h"
#include "RooLinkedList.h"

class RooTObjWrap : public TNamed {
public:

  RooTObjWrap(Bool_t isArray=kFALSE) : _isArray(isArray), _owning(kFALSE) {} ;
  RooTObjWrap(TObject* inObj, Bool_t isArray=kFALSE) : TNamed(), _isArray(isArray), _owning(kFALSE) { if (inObj) _list.Add(inObj) ; }
  RooTObjWrap(const RooTObjWrap& other) : TNamed(other),  _isArray(other._isArray), _owning(kFALSE), _list(other._list) {}
  virtual ~RooTObjWrap() { if (_owning) _list.Delete() ; } ;

  void setOwning(Bool_t flag) { _owning = flag ; }
  TObject* obj() const { return _list.At(0) ; }
  const RooLinkedList& objList() const { return _list ; }

  void setObj(TObject* inObj) {
     if (!_isArray) {
         _list.Clear() ;
     }
    if (inObj) _list.Add(inObj) ;
   }

protected:

  Bool_t _isArray ;
  Bool_t _owning ;
  RooLinkedList _list ;
  ClassDef(RooTObjWrap,2) // Container class for Int_t
};

#endif
