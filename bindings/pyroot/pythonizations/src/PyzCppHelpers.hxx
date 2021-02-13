// Author: Danilo Piparo CERN  08/2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef PYROOT_PYZCPPHELPERS
#define PYROOT_PYZCPPHELPERS

#include "CPyCppyy.h"
#include "CPPInstance.h"
#include "TClass.h"
#include "RConfig.h"

#include <string>

PyObject *CallPyObjMethod(PyObject *obj, const char *meth);
PyObject *CallPyObjMethod(PyObject *obj, const char *meth, PyObject *arg1);
PyObject *BoolNot(PyObject *value);
TClass *GetTClass(const CPyCppyy::CPPInstance *pyobj);
std::string GetCppTypeFromNumpyType(const std::string& dtype);
PyObject *GetArrayInterface(PyObject *obj);
unsigned long long GetDataPointerFromArrayInterface(PyObject *obj);
std::string GetTypestrFromArrayInterface(PyObject *obj);
unsigned int GetDatatypeSizeFromTypestr(const std::string& typestr);
bool CheckEndianessFromTypestr(const std::string& typestr);

#endif // PYROOT_PYZCPPHELPERS
