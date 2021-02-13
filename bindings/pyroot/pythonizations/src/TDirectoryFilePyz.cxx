// Author: Massimiliano Galli CERN  08/2019
// Original PyROOT code by Wim Lavrijsen, LBL

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "CPyCppyy.h"
#include "CPPInstance.h"
#include "CPPOverload.h"
#include "PyROOTPythonize.h"
#include "ProxyWrappers.h"
#include "Python.h"
#include "Utility.h"
#include "PyzCppHelpers.hxx"
#include "TClass.h"
#include "TKey.h"
#include "TDirectoryFile.h"

using namespace CPyCppyy;

////////////////////////////////////////////////////////////////////////////
/// \brief Allow access to objects through the method Get()
/// This concerns both TDirectoryFile and TFile, since the latter
/// inherits the Get method from the former.
/// We decided not to inject this behavior directly in TDirectory
/// because this one already has a templated method Get which, when
/// invoked from Python, returns an object of the derived class (e.g. TH1F)
/// and not a generic TObject.
/// In case the object is not found, a null pointer is returned.
PyObject *TDirectoryFileGetPyz(CPPInstance *self, PyObject *pynamecycle)
{
   // Pythonization of TDirectoryFile::Get that handles non-TObject deriveds
   if (!CPPInstance_Check(self)) {
      PyErr_SetString(PyExc_TypeError,
                      "T(Directory)File::Get must be called with a T(Directory)File instance as first argument");
      return nullptr;
   }
   auto dirf = (TDirectoryFile *)GetTClass(self)->DynamicCast(TDirectoryFile::Class(), self->GetObject());
   if (!dirf) {
      PyErr_SetString(PyExc_ReferenceError, "attempt to access a null-pointer");
      return nullptr;
   }
   const char *namecycle = CPyCppyy_PyText_AsString(pynamecycle);
   if (!namecycle)
      return nullptr; // TypeError already set
   auto key = dirf->GetKey(namecycle);
   if (key) {
      void *addr = dirf->GetObjectChecked(namecycle, key->GetClassName());
      return BindCppObjectNoCast(addr, (Cppyy::TCppType_t)Cppyy::GetScope(key->GetClassName()), kFALSE);
   }
   // no key? for better or worse, call normal Get()
   void *addr = dirf->Get(namecycle);
   return BindCppObject(addr, (Cppyy::TCppType_t)Cppyy::GetScope("TObject"), kFALSE);
}

////////////////////////////////////////////////////////////////////////////
/// \brief Add pythonisation of TDirectoryFile::Get()
/// \param[in] self Always null, since this is a module function.
/// \param[in] args Pointer to a Python tuple object containing the arguments
/// received from Python.
///
/// Allow access to objects through the Get() method.
/// (e.g. dirfile.Get(key))
PyObject *PyROOT::AddTDirectoryFileGetPyz(PyObject * /* self */, PyObject *args)
{
   PyObject *pyclass = PyTuple_GetItem(args, 0);
   Utility::AddToClass(pyclass, "Get", (PyCFunction)TDirectoryFileGetPyz, METH_O);
   Py_RETURN_NONE;
}
