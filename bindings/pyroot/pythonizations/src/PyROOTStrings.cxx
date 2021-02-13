// Author: Enric Tejedor CERN  06/2018
// Original PyROOT code by Wim Lavrijsen, LBL

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// Bindings
#include "CPyCppyy.h"
#include "PyROOTStrings.h"

// Define cached python strings
PyObject *PyROOT::PyStrings::gBranch = nullptr;
PyObject *PyROOT::PyStrings::gFitFCN = nullptr;
PyObject *PyROOT::PyStrings::gROOTns = nullptr;
PyObject *PyROOT::PyStrings::gSetBranchAddress = nullptr;
PyObject *PyROOT::PyStrings::gSetFCN = nullptr;
PyObject *PyROOT::PyStrings::gTClassDynCast = nullptr;
PyObject *PyROOT::PyStrings::gClass = nullptr;

#define PYROOT_INITIALIZE_STRING(var, str)                                    \
   if (!(PyStrings::var = CPyCppyy_PyText_InternFromString((char *)#str))) \
   return false

bool PyROOT::CreatePyStrings()
{
   // Build cache of commonly used python strings (the cache is python intern, so
   // all strings are shared python-wide, not just in PyROOT).
   PYROOT_INITIALIZE_STRING(gBranch, Branch);
   PYROOT_INITIALIZE_STRING(gFitFCN, FitFCN);
   PYROOT_INITIALIZE_STRING(gROOTns, ROOT);
   PYROOT_INITIALIZE_STRING(gSetBranchAddress, SetBranchAddress);
   PYROOT_INITIALIZE_STRING(gSetFCN, SetFCN);
   PYROOT_INITIALIZE_STRING(gTClassDynCast, _TClass__DynamicCast);
   PYROOT_INITIALIZE_STRING(gClass, __class__);

   return true;
}

/// Remove all cached python strings.

PyObject *PyROOT::DestroyPyStrings()
{
   Py_DECREF(PyStrings::gBranch);
   PyStrings::gBranch = nullptr;
   Py_DECREF(PyStrings::gFitFCN);
   PyStrings::gFitFCN = nullptr;
   Py_DECREF(PyStrings::gROOTns);
   PyStrings::gROOTns = nullptr;
   Py_DECREF(PyStrings::gSetBranchAddress);
   PyStrings::gSetBranchAddress = nullptr;
   Py_DECREF(PyStrings::gSetFCN);
   PyStrings::gSetFCN = nullptr;
   Py_DECREF(PyStrings::gTClassDynCast);
   PyStrings::gTClassDynCast = nullptr;
   Py_DECREF(PyStrings::gClass);
   PyStrings::gClass = nullptr;

   Py_INCREF(Py_None);
   return Py_None;
}
