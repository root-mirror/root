// Author: Enric Tejedor CERN  04/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef PYROOT_FACADEHELPERS_H
#define PYROOT_FACADEHELPERS_H

#include "Python.h"

namespace PyROOT {

// Create an indexable buffer starting at the received address
PyObject *CreateBufferFromAddress(PyObject *self, PyObject *addr);

} // namespace PyROOT

#endif // !PYROOT_FACADEHELPERS_H
