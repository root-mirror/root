/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooTreeData.h,v 1.41 2007/07/16 21:04:28 wouter Exp $
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
#ifndef ROO_TREE_DATA
#define ROO_TREE_DATA

#include "RooAbsData.h"
#include "TString.h"
#include "RooArgSet.h"
class TTree ;


class RooTreeData : public RooAbsData {
public:

  // Constructors, factory methods etc.
  RooTreeData() : _tree(0) {} ;
  const TTree *GetTree() const { return _tree; }

private:

  TTree *_tree ;           // TTree holding the data points
  RooArgSet _truth;        // Truth variables
  TString _blindString ;   // Blinding string (optionally read from ASCII files)

  ClassDef(RooTreeData,1) // Dummy class for legacy RooDataSet support
};


#endif
