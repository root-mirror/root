// Author: Danilo Piparo CERN  08/2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**

Set of helper functions that are invoked from the C++ implementation of
pythonizations.

*/
#include "PyzCppHelpers.hxx"

// Call method with signature: obj->meth()
PyObject *CallPyObjMethod(PyObject *obj, const char *meth)
{
   return PyObject_CallMethod(obj, const_cast<char *>(meth), const_cast<char *>(""));
}

// Call method with signature: obj->meth(arg1)
PyObject *CallPyObjMethod(PyObject *obj, const char *meth, PyObject *arg1)
{
   return PyObject_CallMethod(obj, const_cast<char *>(meth), const_cast<char *>("O"), arg1);
}

// Convert generic python object into a boolean value
PyObject *BoolNot(PyObject *value)
{
   if (PyObject_IsTrue(value) == 1) {
      Py_DECREF(value);
      Py_RETURN_FALSE;
   } else {
      Py_XDECREF(value);
      Py_RETURN_TRUE;
   }
}

// Get the TClass of the C++ object proxied by pyobj
TClass *GetTClass(const CPyCppyy::CPPInstance *pyobj)
{
   return TClass::GetClass(Cppyy::GetScopedFinalName(pyobj->ObjectIsA()).c_str());
}

////////////////////////////////////////////////////////////////////////////
/// \brief Convert Numpy data-type string to the according C++ data-type string
/// \param[in] dtype Numpy data-type string
/// \return C++ data-type string
///
/// If the input data-tyep is not known, the function returns an empty string.
std::string GetCppTypeFromNumpyType(const std::string& dtype) {
   if (dtype == "i4") {
      return "int";
   } else if (dtype == "u4") {
      return "unsigned int";
   } else if (dtype == "i8") {
      return "Long64_t";
   } else if (dtype == "u8") {
      return "ULong64_t";
   } else if (dtype == "f4") {
      return "float";
   } else if (dtype == "f8") {
      return "double";
   } else {
      PyErr_SetString(PyExc_RuntimeError, ("Object not convertible: Python object has unknown data-type '" + dtype + "'.").c_str());
      return "";
   }
}

////////////////////////////////////////////////////////////////////////////
/// \brief Get Numpy array interface and perform error handling
/// \param[in] obj PyObject with array interface dictionary
/// \return Array interface dictionary
PyObject *GetArrayInterface(PyObject *obj)
{
   auto pyinterface = PyObject_GetAttrString(obj, "__array_interface__");
   if (!pyinterface) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: __array_interface__ does not exist.");
      return NULL;
   }
   if (!PyDict_Check(pyinterface)) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: __array_interface__ is not a dictionary.");
      return NULL;
   }
   return pyinterface;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Get data pointer from Numpy array interface and perform error handling
/// \param[in] obj Array interface dictionary
/// \return Data pointer
unsigned long long GetDataPointerFromArrayInterface(PyObject *obj)
{
   auto pydata = PyDict_GetItemString(obj, "data");
   if (!pydata) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: __array_interface__['data'] does not exist.");
      return 0;
   }
   return PyLong_AsLong(PyTuple_GetItem(pydata, 0));
}

////////////////////////////////////////////////////////////////////////////
/// \brief Get type string from Numpy array interface and perform error handling
/// \param[in] obj Array interface dictionary
/// \return Type string
std::string GetTypestrFromArrayInterface(PyObject *obj)
{
   auto pytypestr = PyDict_GetItemString(obj, "typestr");
   if (!pytypestr) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: __array_interface__['typestr'] does not exist.");
      return "";
   }
   std::string typestr = CPyCppyy_PyText_AsString(pytypestr);
   const auto length = typestr.length();
   if(length != 3) {
      PyErr_SetString(PyExc_RuntimeError,
              ("Object not convertible: __array_interface__['typestr'] returned '" + typestr + "' with invalid length unequal 3.").c_str());
      return "";
   }
   return typestr;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Get size of data type in bytes from Numpy type string
/// \param[in] typestr Numpy type string
/// \return Size in bytes
unsigned int GetDatatypeSizeFromTypestr(const std::string& typestr)
{
   const auto length = typestr.size();
   const auto dtypesizestr = typestr.substr(length - 1, length);
   return std::stoi(dtypesizestr);
}

////////////////////////////////////////////////////////////////////////////
/// \brief Check whether endianess in type string matches the endianess of ROOT
/// \param[in] typestr Numpy type string
/// \return Boolean indicating whether they match
bool CheckEndianessFromTypestr(const std::string& typestr)
{
   const auto endianess = typestr.substr(1, 2);
#ifdef R__BYTESWAP
   const auto byteswap = "<";
#else
   const auto byteswap = ">";
#endif
   if (!endianess.compare(byteswap)) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: Endianess of __array_interface__['typestr'] does not match endianess of ROOT.");
      return false;
   }
   return true;
}
