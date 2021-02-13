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
#ifndef HISTFACTORY_SIMULTANEOUS
#define HISTFACTORY_SIMULTANEOUS

#include "RooSimultaneous.h"

#include <string>
#include <map>

namespace RooStats{
namespace HistFactory{


class HistFactorySimultaneous : public RooSimultaneous {
public:

  // Constructors, assignment etc
  inline HistFactorySimultaneous() : RooSimultaneous() {} //_plotCoefNormRange(0) { }
  HistFactorySimultaneous(const char *name, const char *title, RooAbsCategoryLValue& indexCat) ;
  HistFactorySimultaneous(const char *name, const char *title, std::map<std::string,RooAbsPdf*> pdfMap, RooAbsCategoryLValue& inIndexCat) ;
  HistFactorySimultaneous(const char *name, const char *title, const RooArgList& pdfList, RooAbsCategoryLValue& indexCat) ;
  HistFactorySimultaneous(const HistFactorySimultaneous& other, const char* name=0);
  HistFactorySimultaneous(const RooSimultaneous& other, const char* name=0);
  ~HistFactorySimultaneous();

  virtual TObject* clone(const char* newname) const { return new HistFactorySimultaneous(*this,newname) ; }

  virtual RooAbsReal* createNLL(RooAbsData& data, const RooLinkedList& cmdList);

  virtual RooAbsReal* createNLL(RooAbsData& data, 
				const RooCmdArg& arg1 = RooCmdArg::none(), const RooCmdArg& arg2 = RooCmdArg::none(), 
				const RooCmdArg& arg3 = RooCmdArg::none(), const RooCmdArg& arg4 = RooCmdArg::none(), 
				const RooCmdArg& arg5 = RooCmdArg::none(), const RooCmdArg& arg6 = RooCmdArg::none(), 
				const RooCmdArg& arg7 = RooCmdArg::none(), const RooCmdArg& arg8 = RooCmdArg::none());
    
  
protected:

  ClassDef(RooStats::HistFactory::HistFactorySimultaneous,2)  // Simultaneous operator p.d.f, functions like C++  'switch()' on input p.d.fs operating on index category5A
};

}
}

#endif
