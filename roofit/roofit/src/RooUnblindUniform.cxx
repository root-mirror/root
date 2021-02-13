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

/** \class RooUnblindUniform
    \ingroup Roofit

Implementation of BlindTools' offset blinding method.
A RooUnblindUniform object is a real-valued function
object, constructed from a parameter to be blinded and a
set of config parameters to change the blinding method.
When supplied to a PDF in lieu of the regular parameter,
a transformation will be applied such that the likelihood is computed with the actual
value of the parameter, but RooFit (, the user, MINUIT) see only
the transformed (blinded) value. The transformation is chosen such that
the error of the blind parameter is identical to that
of the original parameter.
**/

#include "RooFit.h"

#include "RooArgSet.h"
#include "RooUnblindUniform.h"


using namespace std;

ClassImp(RooUnblindUniform);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooUnblindUniform::RooUnblindUniform()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a given RooAbsReal (to hold the blinded value) and a set
/// of blinding parameters.
/// \param name Name of this transformation
/// \param title Title (for plotting)
/// \param blindString String to initialise the random generator
/// \param scale Scale the offset. High values make the blinding more violent.
/// \param blindValue The parameter to be blinded. After the fit, this parameter will
/// only hold the blinded values.

RooUnblindUniform::RooUnblindUniform(const char *name, const char *title,
                const char *blindString, Double_t scale, RooAbsReal& blindValue)
  : RooAbsHiddenReal(name,title),
  _value("value","Uniform blinded value",this,blindValue),
  _blindEngine(blindString,RooBlindTools::full,0.,scale)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooUnblindUniform::RooUnblindUniform(const RooUnblindUniform& other, const char* name) :
  RooAbsHiddenReal(other, name),
  _value("asym",this,other._value),
  _blindEngine(other._blindEngine)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooUnblindUniform::~RooUnblindUniform()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate RooBlindTools unhide-offset method on blind value

Double_t RooUnblindUniform::evaluate() const
{
  return _blindEngine.UnHideUniform(_value);
}
