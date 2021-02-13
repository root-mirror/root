// Author: Enric Tejedor CERN  06/2018
// Original PyROOT code by Wim Lavrijsen, LBL

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef PYROOT_ROOTWRAPPER_H
#define PYROOT_ROOTWRAPPER_H

#include "Python.h"

namespace PyROOT {

// initialize ROOT
void Init();

// clean up all objects controlled by TMemoryRegulator
PyObject *ClearProxiedObjects(PyObject *self, PyObject *args);

} // namespace PyROOT

#endif // !PYROOT_ROOTWRAPPER_H
