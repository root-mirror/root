// Author: Danilo Piparo CERN  08/2018
// Original PyROOT code by Wim Lavrijsen, LBL

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "CPyCppyy.h"
#include "CPPOverload.h"
#include "PyROOTPythonize.h"
#include "Python.h"

using namespace CPyCppyy;

////////////////////////////////////////////////////////////////////////////
/// \brief Make TFile::Open equivalent to a constructor
/// \param[in] self Always null, since this is a module function.
/// \param[in] args Pointer to a Python tuple object containing the arguments
PyObject *PyROOT::AddFileOpenPyz(PyObject * /* self */, PyObject *args)
{
   PyObject *pyclass = PyTuple_GetItem(args, 0);
   // TFile::Open really is a constructor
   PyObject *attr = PyObject_GetAttrString(pyclass, (char *)"Open");
   if (CPPOverload_Check(attr)) {
      ((CPPOverload *)attr)->fMethodInfo->fFlags |= CallContext::kIsCreator;
   }
   Py_XDECREF(attr);
   Py_RETURN_NONE;
}
