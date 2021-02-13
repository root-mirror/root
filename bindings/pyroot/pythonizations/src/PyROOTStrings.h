// Author: Enric Tejedor CERN  06/2018
// Original PyROOT code by Wim Lavrijsen, LBL

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef PYROOT_PYSTRINGS_H
#define PYROOT_PYSTRINGS_H

#include "Python.h"
#include "DllImport.h"

namespace PyROOT {

// python strings kept for performance reasons

namespace PyStrings {

R__EXTERN PyObject *gBranch;
R__EXTERN PyObject *gFitFCN;
R__EXTERN PyObject *gROOTns;
R__EXTERN PyObject *gSetBranchAddress;
R__EXTERN PyObject *gSetFCN;
R__EXTERN PyObject *gTClassDynCast;

R__EXTERN PyObject* gClass;

} // namespace PyStrings

bool CreatePyStrings();
PyObject *DestroyPyStrings();

} // namespace PyROOT

#endif // !PYROOT_PYSTRINGS_H
