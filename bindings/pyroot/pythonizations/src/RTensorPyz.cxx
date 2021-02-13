// Author: Stefan Wunsch CERN  07/2019
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
/// \brief Adopt memory of a Python object with array interface using an RTensor
/// \param[in] obj PyObject with array interface
///
/// This function returns an RTensor which adopts the memory of the given
/// PyObject. The RTensor takes the data pointer and the shape from the array
/// interface dictionary.
PyObject *PyROOT::AsRTensor(PyObject * /*self*/, PyObject * obj)
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
   std::vector<std::size_t> shape;
   for (Py_ssize_t i = 0; i < PyTuple_Size(pyshape); i++) {
      const auto s = PyLong_AsLong(PyTuple_GetItem(pyshape, i));
      shape.push_back(s);
   }

   // Get the typestring and properties thereof
   const auto typestr = GetTypestrFromArrayInterface(pyinterface);
   if (typestr.compare("") == 0)
      return NULL;
   const auto dtypesize = GetDatatypeSizeFromTypestr(typestr);
   if (!CheckEndianessFromTypestr(typestr))
      return NULL;

   const auto dtype = typestr.substr(1, typestr.size());
   std::string cppdtype = GetCppTypeFromNumpyType(dtype);
   if (cppdtype.compare("") == 0)
      return NULL;

   // Get strides
   if (!PyObject_HasAttrString(obj, "strides")) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: Object does not have method 'strides'.");
      return NULL;
   }
   auto pystrides = PyObject_GetAttrString(obj, "strides");
   std::vector<std::size_t> strides;
   for (Py_ssize_t i = 0; i < PyTuple_Size(pystrides); i++) {
      strides.push_back(PyInt_AsLong(PyTuple_GetItem(pystrides, i)) / dtypesize);
   }
   Py_DECREF(pystrides);

   // Infer memory layout from strides
   bool rowMajor = true;
   if (strides.size() > 1) {
      if (strides.front() < strides.back()) rowMajor = false;
   }

   // Construct an RTensor of the correct data-type
   const std::string klassname = "TMVA::Experimental::RTensor<" + cppdtype + ",std::vector<" + cppdtype + ">>";
   std::stringstream code;
   code << "new " << klassname << "(reinterpret_cast<" << cppdtype << "*>(" << data << "),{";
   for (auto s: shape) code << s << ",";
   code << "},{";
   for (auto s: strides) code << s << ",";
   code << "},";
   if (rowMajor) {
      code << "TMVA::Experimental::MemoryLayout::RowMajor";
   }
   else {
      code << "TMVA::Experimental::MemoryLayout::ColumnMajor";
   }
   code << ")";
   const auto codestr = code.str();
   auto address = (void*) gInterpreter->Calc(codestr.c_str());

   // Bind the object to a Python-side proxy
   auto klass = (Cppyy::TCppType_t)Cppyy::GetScope(klassname);
   auto pyobj = CPyCppyy::BindCppObject(address, klass);

   // Give Python the ownership of the underlying C++ object
   ((CPyCppyy::CPPInstance*)pyobj)->PythonOwns();

   // Bind pyobject holding adopted memory to the RTensor
   if (PyObject_SetAttrString(pyobj, "__adopted__", obj)) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: Failed to set Python object as attribute __adopted__.");
      return NULL;
   }

   // Clean-up and return
   Py_DECREF(pyinterface);
   return pyobj;
}
