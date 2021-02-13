/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
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

/**
\file RooRangeBinning.cxx
\class RooRangeBinning
\ingroup Roofitcore

RooRangeBinning is binning/range definition that only defines a range
but no binning. It it used to store named ranges created by
the RooRealVar::setRange() method
**/

#include "RooFit.h"

#include "RooNumber.h"
#include "RooMsgService.h"
#include "Riostream.h"

#include "RooRangeBinning.h"

using namespace std;

ClassImp(RooRangeBinning); 
;



////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooRangeBinning::RooRangeBinning(const char* name) :
  RooAbsBinning(name)
{
  _range[0] = -RooNumber::infinity() ;
  _range[1] = +RooNumber::infinity() ;

}


////////////////////////////////////////////////////////////////////////////////
/// Construct binning with range [xmin,xmax] with no binning substructure

RooRangeBinning::RooRangeBinning(Double_t xmin, Double_t xmax, const char* name) :
  RooAbsBinning(name)
{
  _range[0] = xmin ;
  _range[1] = xmax ;
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooRangeBinning::RooRangeBinning(const RooRangeBinning& other, const char* name) :
  RooAbsBinning(name)
{
  _range[0] = other._range[0] ;
  _range[1] = other._range[1] ;
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor 

RooRangeBinning::~RooRangeBinning() 
{
}



////////////////////////////////////////////////////////////////////////////////
/// Change limits of the binning to [xlo,xhi]

void RooRangeBinning::setRange(Double_t xlo, Double_t xhi) 
{
  if (xlo>xhi) {
    oocoutE((TObject*)0,InputArguments) << "RooRangeBinning::setRange: ERROR low bound > high bound" << endl ;
    return ;
  }

  _range[0] = xlo ;
  _range[1] = xhi ;
}
