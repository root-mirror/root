// Author: Stefan Wunsch CERN  02/2019
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
#include "ProxyWrappers.h"
#include "PyROOTPythonize.h"
#include "TInterpreter.h"
#include "PyzCppHelpers.hxx"

#include <sstream>

////////////////////////////////////////////////////////////////////////////
/// \brief Adopt memory of a Python object with array interface using an RVec
/// \param[in] obj PyObject with array interface
///
/// This function returns an RVec which adopts the memory of the given
/// PyObject. The RVec takes the data pointer and the size from the array
/// interface dictionary.
PyObject *PyROOT::AsRVec(PyObject * /*self*/, PyObject * obj)
{
   if (!obj) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: Invalid Python object.");
      return NULL;
   }

   // Get array interface of object
   auto pyinterface = GetArrayInterface(obj);
   if (pyinterface == NULL)
      return NULL;

   // Get the data-pointer
   const auto data = GetDataPointerFromArrayInterface(pyinterface);
   if (data == 0)
      return NULL;

   // Get the size of the contiguous memory
   auto pyshape = PyDict_GetItemString(pyinterface, "shape");
   if (!pyshape) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: __array_interface__['shape'] does not exist.");
      return NULL;
   }
   long size = 0;
   for (int i = 0; i < PyTuple_Size(pyshape); i++) {
      if (size == 0) size = 1;
      size *= PyLong_AsLong(PyTuple_GetItem(pyshape, i));
   }

   // Get the typestring and properties thereof
   const auto typestr = GetTypestrFromArrayInterface(pyinterface);
   if (typestr.compare("") == 0)
      return NULL;
   if (!CheckEndianessFromTypestr(typestr))
      return NULL;

   const auto dtype = typestr.substr(1, typestr.size());
   std::string cppdtype = GetCppTypeFromNumpyType(dtype);
   if (cppdtype.compare("") == 0)
      return NULL;

   // Construct an RVec of the correct data-type
   const std::string klassname = "ROOT::VecOps::RVec<" + cppdtype + ">";
   std::stringstream prefix;
#ifdef _MSC_VER
   prefix << "0x";
#endif
   auto address = (void*) gInterpreter->Calc("new " + klassname + "(reinterpret_cast<" + cppdtype + "*>(" + prefix.str() + data + ")," + size + ")");

   // Bind the object to a Python-side proxy
   auto klass = (Cppyy::TCppType_t)Cppyy::GetScope(klassname);
   auto pyobj = CPyCppyy::BindCppObject(address, klass);

   // Give Python the ownership of the underlying C++ object
   ((CPyCppyy::CPPInstance*)pyobj)->PythonOwns();

   // Bind pyobject holding adopted memory to the RVec
   if (PyObject_SetAttrString(pyobj, "__adopted__", obj)) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: Failed to set Python object as attribute __adopted__.");
      return NULL;
   }

   // Clean-up and return
   Py_DECREF(pyinterface);
   return pyobj;
}
