/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooCategorySharedProperties.h,v 1.2 2007/05/11 09:11:30 verkerke Exp $
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
#ifndef ROO_CATEGORY_SHARED_PROPERTY
#define ROO_CATEGORY_SHARED_PROPERTY

#include "RooSharedProperties.h"
#include "RooLinkedList.h"

class RooCategorySharedProperties : public RooSharedProperties {
public:

  RooCategorySharedProperties() ;
  RooCategorySharedProperties(const char* uuidstr) ;
  RooCategorySharedProperties(const RooCategorySharedProperties& other) ;
  virtual ~RooCategorySharedProperties() ;

  RooSharedProperties* clone() { return new RooCategorySharedProperties(*this)  ; }

protected:

  friend class RooCategory ;

  RooLinkedList _altRanges ;  // Optional alternative ranges

  ClassDef(RooCategorySharedProperties,1) // Shared properties of a RooCategory clone set
};


#endif
