// Author: Enric Tejedor CERN  04/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// Bindings
#include "FacadeHelpers.hxx"

// Cppyy
#include "LowLevelViews.h"

// ROOT
#include "RtypesCore.h"
#include "ROOT/RConfig.hxx"

//////////////////////////////////////////////////////////////////////////
/// \brief Get a buffer starting at a given address
/// \param[in] self Always null, since this is a module function.
/// \param[in] addr Address to create buffer from
///
/// \return A cppyy LowLevelView object on the received address, i.e. an
///         indexable buffer starting at that address.
PyObject *PyROOT::CreateBufferFromAddress(PyObject * /* self */, PyObject *addr)
{
   if (!addr) {
      PyErr_SetString(PyExc_RuntimeError, "Unable to create buffer from invalid address");
      return NULL;
   }

   Long64_t cAddr = PyLong_AsLongLong(addr);
   if (cAddr == -1 && PyErr_Occurred()) {
      PyErr_SetString(PyExc_RuntimeError, "Unable to create buffer: address is not a valid integer");
      return NULL;
   }

#ifdef R__B64
   return CPyCppyy::CreateLowLevelView((Long64_t*)cAddr);
#else
   return CPyCppyy::CreateLowLevelView((Int_t*)cAddr);
#endif
}
