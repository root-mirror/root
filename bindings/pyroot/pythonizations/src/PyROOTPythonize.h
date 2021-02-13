// Author: Enric Tejedor CERN  06/2018
// Original PyROOT code by Wim Lavrijsen, LBL

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef PYROOT_PYTHONIZE_H
#define PYROOT_PYTHONIZE_H

#include "Python.h"

namespace PyROOT {

PyObject *AddCPPInstancePickling(PyObject *self, PyObject *args);

PyObject *AddPrettyPrintingPyz(PyObject *self, PyObject *args);

PyObject *AddDirectoryGetAttrPyz(PyObject *self, PyObject *args);
PyObject *AddDirectoryWritePyz(PyObject *self, PyObject *args);
PyObject *AddTDirectoryFileGetPyz(PyObject *self, PyObject *args);
PyObject *AddFileOpenPyz(PyObject *self, PyObject *args);

PyObject *AddBranchAttrSyntax(PyObject *self, PyObject *args);
PyObject *BranchPyz(PyObject *self, PyObject *args);
PyObject *SetBranchAddressPyz(PyObject *self, PyObject *args);

PyObject *AddTClassDynamicCastPyz(PyObject *self, PyObject *args);

PyObject *AddTObjectEqNePyz(PyObject *self, PyObject *args);

PyObject *AddSetItemTCAPyz(PyObject *self, PyObject *args);

PyObject *AddUsingToClass(PyObject *self, PyObject *args);

PyObject *AsRVec(PyObject *self, PyObject *obj);
PyObject *AsRTensor(PyObject *self, PyObject *obj);

PyObject *CPPInstanceExpand(PyObject *self, PyObject *args);

PyObject *GetEndianess(PyObject *self, PyObject *args);
PyObject *GetDataPointer(PyObject *self, PyObject *args);
PyObject *GetSizeOfType(PyObject *self, PyObject *args);

PyObject *MakeNumpyDataFrame(PyObject *self, PyObject *obj);

} // namespace PyROOT

#endif // !PYROOT_PYTHONIZE_H
