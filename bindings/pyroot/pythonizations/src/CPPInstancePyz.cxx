// Author: Massimiliano Galli CERN  07/2019
// Original PyROOT code by Wim Lavrijsen, LBL

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// Bindings
#include "CPyCppyy.h"
#include "PyROOTPythonize.h"
#include "CPPInstance.h"
#include "ProxyWrappers.h"
#include "Converters.h"
#include "Utility.h"
#include "PyzCppHelpers.hxx"
#include "TBufferFile.h"
#include "CustomPyTypes.h"

using namespace CPyCppyy;

namespace PyROOT {
extern PyObject *gRootModule;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Deserialize pickled objects
/// \param[in] self Always null, since this is a module function.
/// \param[in] args Pointer to a Python tuple object containing the arguments
/// received from Python.
///
/// Helper function that deserializes pickled objects. It needs to be
/// included in the extension module API because otherwise it is not
/// callable from Python. This is important because it will be Python
/// itself calling it when trying to expand a serialized object.
PyObject *PyROOT::CPPInstanceExpand(PyObject * /*self*/, PyObject *args)
{
   PyObject *pybuf = 0, *pyname = 0;
   if (!PyArg_ParseTuple(args, const_cast<char *>("O!O!:__expand__"), &PyBytes_Type, &pybuf, &PyBytes_Type, &pyname))
      return 0;
   const char *clname = PyBytes_AS_STRING(pyname);
   // TBuffer and its derived classes can't write themselves, but can be created
   // directly from the buffer, so handle them in a special case
   void *newObj = 0;
   if (strcmp(clname, "TBufferFile") == 0) {
      TBufferFile *buf = new TBufferFile(TBuffer::kWrite);
      buf->WriteFastArray(PyBytes_AS_STRING(pybuf), PyBytes_GET_SIZE(pybuf));
      newObj = buf;
   } else {
      // use the PyString macro's to by-pass error checking; do not adopt the buffer,
      // as the local TBufferFile can go out of scope (there is no copying)
      TBufferFile buf(TBuffer::kRead, PyBytes_GET_SIZE(pybuf), PyBytes_AS_STRING(pybuf), kFALSE);
      newObj = buf.ReadObjectAny(0);
   }
   PyObject *result = BindCppObject(newObj, Cppyy::GetScope(clname));
   if (result) {
      // this object is to be owned by the Python interpreter, assuming that the call
      // originated from there
      ((CPPInstance *)result)->PythonOwns();
   }
   return result;
}

/// PyROOT object proxy pickle support
/// Turn the object proxy instance into a character stream and return for
/// pickle, together with the callable object that can restore the stream
/// into the object proxy instance.
PyObject *op_reduce(CPPInstance *self, PyObject * /*args*/)
{
   // keep a borrowed reference around to the callable function for expanding;
   // because it is borrowed, it means that there can be no pickling during the
   // shutdown of the libPyROOT module
   static PyObject *s_expand =
      PyDict_GetItemString(PyModule_GetDict(PyROOT::gRootModule), const_cast<char *>("_CPPInstance__expand__"));

   // TBuffer and its derived classes can't write themselves, but can be created
   // directly from the buffer, so handle them in a special case
   static Cppyy::TCppType_t s_bfClass = Cppyy::GetScope("TBufferFile");
   TBufferFile *buff = 0;
   if (s_bfClass == self->ObjectIsA()) {
      buff = (TBufferFile *)self->GetObject();
   } else {
      // no cast is needed, but WriteObject taking a TClass argument is protected,
      // so use WriteObjectAny()
      static TBufferFile s_buff(TBuffer::kWrite);
      s_buff.Reset();
      // to delete
      if (s_buff.WriteObjectAny(self->GetObject(),
                                TClass::GetClass(Cppyy::GetScopedFinalName(self->ObjectIsA()).c_str())) != 1) {
         PyErr_Format(PyExc_IOError, "could not stream object of type %s",
                      Cppyy::GetScopedFinalName(self->ObjectIsA()).c_str());
         return 0;
      }
      buff = &s_buff;
   }
   // use a string for the serialized result, as a python buffer will not copy
   // the buffer contents; use a string for the class name, used when casting
   // on reading back in (see CPPInstanceExpand defined above)
   PyObject *res2 = PyTuple_New(2);
   PyTuple_SET_ITEM(res2, 0, PyBytes_FromStringAndSize(buff->Buffer(), buff->Length()));
   PyTuple_SET_ITEM(res2, 1, PyBytes_FromString(Cppyy::GetScopedFinalName(self->ObjectIsA()).c_str()));

   PyObject *result = PyTuple_New(2);
   Py_INCREF(s_expand);
   PyTuple_SET_ITEM(result, 0, s_expand);
   PyTuple_SET_ITEM(result, 1, res2);

   return result;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Set __reduce__ attribute for CPPInstance objects
/// \param[in] self Always null, since this is a module function.
/// \param[in] args Pointer to a Python tuple object containing the arguments
/// received from Python.
///
/// The C++ function op_reduce defined above is wrapped in a Python method
/// so that it can be injected in CPPInstance
PyObject *PyROOT::AddCPPInstancePickling(PyObject * /*self*/, PyObject *args)
{
   PyObject *pyclass = PyTuple_GetItem(args, 0);

   const char *attr = "__reduce__";

   PyMethodDef *pdef = new PyMethodDef();
   pdef->ml_name = attr;
   pdef->ml_meth = (PyCFunction)op_reduce;
   pdef->ml_flags = METH_NOARGS;
   pdef->ml_doc = nullptr;

   PyObject *func = PyCFunction_New(pdef, nullptr);
   PyObject *method = CustomInstanceMethod_New(func, nullptr, pyclass);

   // here PyObject_GenericSetAttr is used because CPPInstance does not allow
   // attribute assignment using PyObject_SetAttr
   // for more info refer to:
   // https://bitbucket.org/wlav/cppyy/issues/110/user-defined-classes-in-c-dont-seem-to-be
   PyObject_GenericSetAttr(pyclass, CPyCppyy_PyText_FromString(attr), method);
   Py_DECREF(method);
   Py_DECREF(func);

   Py_RETURN_NONE;
}
