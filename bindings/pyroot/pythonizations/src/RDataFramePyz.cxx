// Author: Stefan Wunsch CERN  04/2019
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
#include "RConfig.h"
#include "TInterpreter.h"
#include "CPyCppyy/API.h"

#include <utility> // std::pair
#include <sstream> // std::stringstream

////////////////////////////////////////////////////////////////////////////
/// \brief Make an RDataFrame from a dictionary of numpy arrays
/// \param[in] pydata Dictionary with numpy arrays
///
/// This function takes a dictionary of numpy arrays and creates an RDataFrame
/// using the keys as column names and the numpy arrays as data.
PyObject *PyROOT::MakeNumpyDataFrame(PyObject * /*self*/, PyObject * pydata)
{
   if (!pydata) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: Invalid Python object.");
      return NULL;
   }

   if (!PyDict_Check(pydata)) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: Python object is not a dictionary.");
      return NULL;
   }

   if (PyDict_Size(pydata) == 0) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: Dictionary is empty.");
      return NULL;
   }


   // Add PyObject (dictionary) holding RVecs to data source
   std::stringstream code;
   code << "ROOT::Internal::RDF::MakeNumpyDataFrame(";
   std::stringstream pyaddress;
   auto pyvecs = PyDict_New();
#ifdef _MSC_VER
   pyaddress << "0x";
#endif
   pyaddress << pyvecs;
   code << "reinterpret_cast<PyObject*>(" << pyaddress.str() << "), ";

   // Iterate over dictionary, convert numpy arrays to RVecs and put together interpreter code
   PyObject *key, *value;
   Py_ssize_t pos = 0;
   const auto size = PyObject_Size(pydata);
   auto counter = 0u;
   while (PyDict_Next(pydata, &pos, &key, &value)) {
      // Get name of key
      if (!CPyCppyy_PyText_Check(key)) {
         PyErr_SetString(PyExc_RuntimeError, "Object not convertible: Dictionary key is not convertible to a string.");
         return NULL;
      }
      std::string keystr = CPyCppyy_PyText_AsString(key);

      // Convert value to RVec and attach to dictionary
      auto pyvec = PyROOT::AsRVec(NULL, value);
      if (pyvec == NULL) {
         PyErr_SetString(PyExc_RuntimeError,
                         ("Object not convertible: Dictionary entry " + keystr + " is not convertible with AsRVec.").c_str());
         return NULL;
      }
      PyDict_SetItem(pyvecs, key, pyvec);
      Py_DECREF(pyvec);

      // Add pairs of column name and associated RVec to signature
      std::string vectype = Cppyy::GetScopedFinalName(((CPyCppyy::CPPInstance*)pyvec)->ObjectIsA());
      std::stringstream vecaddress;
#ifdef _MSC_VER
      vecaddress << "0x";
#endif
      vecaddress << ((CPyCppyy::CPPInstance*)pyvec)->GetObject();
      code << "std::pair<std::string, " << vectype <<  "*>(\"" + keystr
           << "\", reinterpret_cast<" << vectype+ "*>(" << vecaddress.str() << "))";
      if (counter != size - 1) {
         code << ",";
      } else {
         code << ");";
      }
      counter++;
   }

   // Create RDataFrame and build Python proxy
   const auto err = gInterpreter->Declare("#include \"ROOT/RNumpyDS.hxx\"");
   if (!err) {
      PyErr_SetString(PyExc_RuntimeError, "Failed to find \"ROOT/RNumpyDS.hxx\".");
      return NULL;
   }
   const auto codeStr = code.str();
   auto address = (void*) gInterpreter->Calc(codeStr.c_str());
   const auto pythonOwns = true;
   auto pyobj = CPyCppyy::Instance_FromVoidPtr(address, "ROOT::RDataFrame", pythonOwns);

   // Bind pyobject holding adopted memory to the RVec
   if (PyObject_SetAttrString(pyobj, "__data__", pyvecs)) {
      PyErr_SetString(PyExc_RuntimeError, "Object not convertible: Failed to set dictionary as attribute __data__.");
      return NULL;
   }
   Py_DECREF(pyvecs);

   return pyobj;
}
