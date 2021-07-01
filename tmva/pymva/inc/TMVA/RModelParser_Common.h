// @(#)root/tmva/pymva $Id$

#ifndef TMVA_PYMVA_RMODELPARSER_COMMON
#define TMVA_PYMVA_RMODELPARSER_COMMON

#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "TMVA/RTensor.hxx"
#include "TString.h"


#include "TMVA/RModel.hxx"
#include "TMVA/SOFIE_common.hxx"
#include<algorithm>

namespace TMVA{

static PyObject *fPyReturn = NULL;

///////////////////////////////////////////////////////////////////////////////
/// Execute Python code from string
///
/// \param[in] code Python code as string
/// \param[in] fGlobalNS Global Namespace dictionary
/// \param[in] fLocalNS  Local  Namespace dictionary
/// \param[in] errorMessage Error message which shall be shown if the execution fails
/// \param[in] start Start symbol
///
/// Helper function to run python code from string in local namespace with
/// error handling
/// `start` defines the start symbol defined in PyRun_String (Py_eval_input,
/// Py_single_input, Py_file_input)

void PyRunString(TString code, PyObject *fGlobalNS, PyObject *fLocalNS, TString errorMessage="Failed to run python code", int start=Py_single_input) {
   fPyReturn = PyRun_String(code, start, fGlobalNS, fLocalNS);
   if (!fPyReturn) {
      std::cout<<"Failed to run python code: " << code <<"\n";
      std::cout<< "Python error message:\n";
      PyErr_Print();
      std::cout<<errorMessage;
   }
 }

namespace Experimental{
   namespace SOFIE{
      RTensor<float> getArray(PyObject* value){
         //Check and modify the function signature
         PyArrayObject* weightArray = (PyArrayObject*)value;
         std::vector<std::size_t>shape;
         std::vector<std::size_t>strides;

         //Preparing the shape vector
         for(npy_intp* j=PyArray_SHAPE(weightArray); j<PyArray_SHAPE(weightArray)+PyArray_NDIM(weightArray); ++j){
            shape.push_back((std::size_t)(*j));
            }

         //Preparing the strides vector
         for(npy_intp* k=PyArray_STRIDES(weightArray); k<PyArray_STRIDES(weightArray)+PyArray_NDIM(weightArray); ++k){
            strides.push_back((std::size_t)(*k));
            }

         //Declaring the RTensor object for storing weights values.
         RTensor<float>x((float*)PyArray_DATA(weightArray),shape,strides);
         return x;
         }
         }
         }

const char* PyStringAsString(PyObject* str){
   #if PY_MAJOR_VERSION < 3   // for Python2
      const char *stra_name = PyBytes_AsString(str);
      // need to add string delimiter for Python2
      TString sname = TString::Format("'%s'",stra_name);
      const char * name = sname.Data();
#else   // for Python3
      PyObject* repr = PyObject_Repr(str);
      PyObject* stra = PyUnicode_AsEncodedString(repr, "utf-8", "~E~");
      const char *name = PyBytes_AsString(stra);
#endif
return name;
}

std::string toLower(std::string name){
          std::transform(name.begin(), name.end(), name.begin(), ::tolower);
          return name;
   }

}

#endif //TMVA_PYMVA_RMODELPARSER_COMMON
