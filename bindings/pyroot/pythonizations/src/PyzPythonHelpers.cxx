// Author: Stefan Wunsch CERN  08/2018
// Original PyROOT code by Wim Lavrijsen, LBL

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**

Set of helper functions that are invoked from the pythonizors, on the
Python side. For that purpose, they are included in the interface of the
PyROOT extension module.

*/

#include "CPyCppyy.h"
#include "CPPInstance.h"
#include "CPPOverload.h"

#include "PyROOTPythonize.h"

#include "RConfig.h"
#include "TInterpreter.h"

#include <sstream>

// needed to properly resolve (dllimport) symbols on Windows
namespace CPyCppyy {
   namespace PyStrings {
      R__EXTERN PyObject *gMRO;
   }
}

////////////////////////////////////////////////////////////////////////////
/// \brief Get size of C++ data-type
/// \param[in] self Always null, since this is a module function.
/// \param[in] args C++ data-type as Python string
///
/// This function returns the length of a C++ data-type in bytes
/// as a Python integer.
PyObject *PyROOT::GetSizeOfType(PyObject * /*self*/, PyObject *args)
{
   // Get name of data-type
   PyObject *pydtype = PyTuple_GetItem(args, 0);
   std::string dtype = CPyCppyy_PyText_AsString(pydtype);

   // Call interpreter to get size of data-type using `sizeof`
   size_t size = 0;
   std::stringstream code;
   code << "*((size_t*)" << std::showbase << (uintptr_t)&size << ") = (size_t)sizeof(" << dtype << ")";
   gInterpreter->Calc(code.str().c_str());

   // Return size of data-type as integer
   PyObject *pysize = PyInt_FromLong(size);
   return pysize;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Get pointer to the data of an object
/// \param[in] self Always null, since this is a module function.
/// \param[in] args[0] Python representation of the C++ object.
/// \param[in] args[1] Data-type of the C++ object as Python string
/// \param[in] args[2] Method to be called on the C++ object to get the data pointer as Python string
///
/// This function returns the pointer to the data of an object as an Python
/// integer retrieved by the given method.
PyObject *PyROOT::GetDataPointer(PyObject * /*self*/, PyObject *args)
{
   // Get pointer of C++ object
   PyObject *pyobj = PyTuple_GetItem(args, 0);
   auto instance = (CPyCppyy::CPPInstance *)(pyobj);
   auto cppobj = instance->GetObject();

   // Get name of C++ object as string
   PyObject *pycppname = PyTuple_GetItem(args, 1);
   std::string cppname = CPyCppyy_PyText_AsString(pycppname);

   // Get name of method to be called to get the data pointer
   PyObject *pymethodname = PyTuple_GetItem(args, 2);
   std::string methodname = CPyCppyy_PyText_AsString(pymethodname);

   // Call interpreter to get pointer to data
   uintptr_t pointer = 0;
   std::stringstream code;
   code << "*((intptr_t*)" << std::showbase << (uintptr_t)&pointer << ") = reinterpret_cast<uintptr_t>(reinterpret_cast<"
        << cppname << "*>(" << std::showbase << (uintptr_t)cppobj << ")->" << methodname << "())";
   gInterpreter->Calc(code.str().c_str());

   // Return pointer as integer
   PyObject *pypointer = PyLong_FromUnsignedLongLong(pointer);
   return pypointer;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Get endianess of the system
/// \param[in] self Always null, since this is a module function.
/// \param[in] args Pointer to an empty Python tuple.
/// \return Endianess as Python string
///
/// This function returns endianess of the system as a Python integer. The
/// return value is either '<' or '>' for little or big endian, respectively.
PyObject *PyROOT::GetEndianess(PyObject * /* self */, PyObject * /* args */)
{
#ifdef R__BYTESWAP
   return CPyCppyy_PyText_FromString("<");
#else
   return CPyCppyy_PyText_FromString(">");
#endif
}

using namespace CPyCppyy;

// Helper to add base class methods to the derived class one
static bool AddUsingToClass(PyObject *pyclass, const char *method)
{
   CPPOverload *derivedMethod = (CPPOverload *)PyObject_GetAttrString(pyclass, const_cast<char *>(method));
   if (!CPPOverload_Check(derivedMethod)) {
      Py_XDECREF(derivedMethod);
      return false;
   }

   PyObject *mro = PyObject_GetAttr(pyclass, PyStrings::gMRO);
   if (!mro || !PyTuple_Check(mro)) {
      Py_XDECREF(mro);
      Py_DECREF(derivedMethod);
      return false;
   }

   CPPOverload *baseMethod = nullptr;
   for (int i = 1; i < PyTuple_GET_SIZE(mro); ++i) {
      baseMethod = (CPPOverload *)PyObject_GetAttrString(PyTuple_GET_ITEM(mro, i), const_cast<char *>(method));

      if (!baseMethod) {
         PyErr_Clear();
         continue;
      }

      if (CPPOverload_Check(baseMethod))
         break;

      Py_DECREF(baseMethod);
      baseMethod = nullptr;
   }

   Py_DECREF(mro);

   if (!CPPOverload_Check(baseMethod)) {
      Py_XDECREF(baseMethod);
      Py_DECREF(derivedMethod);
      return false;
   }

   for (PyCallable *pc : baseMethod->fMethodInfo->fMethods) {
      derivedMethod->AdoptMethod(pc->Clone());
   }

   Py_DECREF(baseMethod);
   Py_DECREF(derivedMethod);

   return true;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Add base class overloads of a given method to a derived class
/// \param[in] self Always null, since this is a module function.
/// \param[in] args[0] Derived class.
/// \param[in] args[1] Name of the method whose base class overloads to
///                    inject in the derived class.
///
/// This function adds base class overloads to a derived class for a given
/// method. This covers the 'using' case, which is not supported by default
/// by the bindings.
PyObject *PyROOT::AddUsingToClass(PyObject * /* self */, PyObject *args)
{
   // Get derived class to pythonize
   PyObject *pyclass = PyTuple_GetItem(args, 0);

   // Get method name where to add overloads
   PyObject *pyname = PyTuple_GetItem(args, 1);
   auto cppname = CPyCppyy_PyText_AsString(pyname);

   AddUsingToClass(pyclass, cppname);

   Py_RETURN_NONE;
}
