/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
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

/** \class RooUnblindCPAsymVar
    \ingroup Roofit

Implementation of BlindTools' CP asymmetry blinding method
A RooUnblindCPAsymVar object is a real valued function
object, constructed from a blind value holder and a
set of unblinding parameters. When supplied to a PDF
in lieu of a regular parameter, the blind value holder
supplied to the unblinded objects will in a fit be minimized
to blind value corresponding to the actual minimum of the
parameter. The transformation is chosen such that the
the error on the blind parameters is identical to that
of the unblind parameter
**/

#include "RooFit.h"

#include "RooArgSet.h"
#include "RooUnblindCPAsymVar.h"


using namespace std;

ClassImp(RooUnblindCPAsymVar);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooUnblindCPAsymVar::RooUnblindCPAsymVar()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a given RooAbsReal (to hold the blind value) and a set of blinding parameters

RooUnblindCPAsymVar::RooUnblindCPAsymVar(const char *name, const char *title,
                const char *blindString, RooAbsReal& cpasym)
  : RooAbsHiddenReal(name,title),
  _asym("asym","CP Asymmetry",this,cpasym),
  _blindEngine(blindString)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a given RooAbsReal (to hold the blind value) and a set of blinding parameters

RooUnblindCPAsymVar::RooUnblindCPAsymVar(const char *name, const char *title,
                const char *blindString, RooAbsReal& cpasym, RooAbsCategory& blindState)
  : RooAbsHiddenReal(name,title,blindState),
  _asym("asym","CP Asymmetry",this,cpasym),
  _blindEngine(blindString)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooUnblindCPAsymVar::RooUnblindCPAsymVar(const RooUnblindCPAsymVar& other, const char* name) :
  RooAbsHiddenReal(other, name),
  _asym("asym",this,other._asym),
  _blindEngine(other._blindEngine)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooUnblindCPAsymVar::~RooUnblindCPAsymVar()
{
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooUnblindCPAsymVar::evaluate() const
{
  if (isHidden()) {
    // Blinding active for this event
    return _blindEngine.UnHideAsym(_asym);
  } else {
    // Blinding not active for this event
    return _asym ;
  }
}
