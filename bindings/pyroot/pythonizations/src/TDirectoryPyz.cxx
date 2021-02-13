// Author: Danilo Piparo, Stefan Wunsch, Massimiliano Galli CERN  08/2018
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
#include "PyROOTPythonize.h"
#include "ProxyWrappers.h"
#include "Utility.h"
#include "PyzCppHelpers.hxx"
#include "TClass.h"
#include "TDirectory.h"
#include "TKey.h"
#include "Python.h"

using namespace CPyCppyy;

////////////////////////////////////////////////////////////////////////////////
/// \brief Implements the WriteObject method of TDirectory
/// This method allows to write objects into TDirectory instances with this
/// syntax:
/// ~~~{.python}
/// myDir.WriteObject(myObj, "myKeyName")
/// ~~~
PyObject *TDirectoryWriteObject(CPPInstance *self, PyObject *args)
{
   CPPInstance *wrt = nullptr;
   PyObject *name = nullptr;
   PyObject *option = nullptr;
   Int_t bufsize = 0;
   if (!PyArg_ParseTuple(args, const_cast<char *>("O!O!|O!i:TDirectory::WriteObject"),
                         &CPPInstance_Type, &wrt,
                         &CPyCppyy_PyText_Type, &name,
                         &CPyCppyy_PyText_Type, &option,
                         &bufsize))
      return nullptr;
   auto dir = (TDirectory *)GetTClass(self)->DynamicCast(TDirectory::Class(), self->GetObject());
   if (!dir) {
      PyErr_SetString(PyExc_TypeError,
                      "TDirectory::WriteObject must be called with a TDirectory instance as first argument");
      return nullptr;
   }
   Int_t result = 0;
   if (option != nullptr) {
      result = dir->WriteObjectAny(wrt->GetObject(), GetTClass(wrt), CPyCppyy_PyText_AsString(name),
                                   CPyCppyy_PyText_AsString(option), bufsize);
   } else {
      result = dir->WriteObjectAny(wrt->GetObject(), GetTClass(wrt), CPyCppyy_PyText_AsString(name));
   }
   return PyInt_FromLong((Long_t)result);
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Implements a getter to assign to TDirectory.__getattr__
/// Method that is assigned to TDirectory.__getattr__. It relies on Get to
/// obtain the object from the TDirectory and adds on top:
/// - Raising an AttributeError if the object does not exist
/// - Caching the result of a successful get for future re-attempts.
/// Once cached, the same object is retrieved every time.
/// This pythonisation is inherited by TDirectoryFile and TFile.
PyObject *TDirectoryGetAttr(PyObject *self, PyObject *attr)
{
   // Injection of TDirectory.__getattr__ that raises AttributeError on failure.
   PyObject *result = CallPyObjMethod(self, "Get", attr);
   if (!result)
      return result;

   if (!PyObject_IsTrue(result)) {
      PyObject *astr = PyObject_Str(attr);
      PyObject *stypestr = PyObject_Str(PyObject_Type(self));
      PyErr_Format(PyExc_AttributeError, "%s object has no attribute \'%s\'", CPyCppyy_PyText_AsString(stypestr),
                   CPyCppyy_PyText_AsString(astr));
      Py_DECREF(astr);
      Py_DECREF(result);
      return nullptr;
   }

   // Caching behavior seems to be more clear to the user; can always override said
   // behavior (i.e. re-read from file) with an explicit Get() call
   PyObject_SetAttr(self, attr, result);
   return result;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Add attr syntax to TDirectory
/// \param[in] self Always null, since this is a module function.
/// \param[in] args Pointer to a Python tuple object containing the arguments
/// This allows to use TDirectory and daughters (such as TDirectoryFile and TFile)
/// as follows
/// ~~~{.python}
/// myfile.mydir.mysubdir.myHist.Draw()
/// ~~~
PyObject *PyROOT::AddDirectoryGetAttrPyz(PyObject * /* self */, PyObject *args)
{
   PyObject *pyclass = PyTuple_GetItem(args, 0);
   Utility::AddToClass(pyclass, "__getattr__", (PyCFunction)TDirectoryGetAttr, METH_O);
   Py_RETURN_NONE;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Add pythonisation of TDirectory::WriteObject
/// \param[in] self Always null, since this is a module function.
/// \param[in] args Pointer to a Python tuple object containing the arguments
PyObject *PyROOT::AddDirectoryWritePyz(PyObject * /* self */, PyObject *args)
{
   PyObject *pyclass = PyTuple_GetItem(args, 0);
   Utility::AddToClass(pyclass, "WriteObject", (PyCFunction)TDirectoryWriteObject);
   Py_RETURN_NONE;
}
