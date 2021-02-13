// Author: Enric Tejedor CERN  06/2018
// Original PyROOT code by Wim Lavrijsen, LBL

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// Bindings
#include "PyROOTWrapper.h"
#include "TMemoryRegulator.h"

// Cppyy
#include "CPyCppyy.h"
#include "ProxyWrappers.h"

// ROOT
#include "TROOT.h"
#include "TSystem.h"
#include "TClass.h"
#include "TInterpreter.h"
#include "DllImport.h"

namespace PyROOT {
R__EXTERN PyObject *gRootModule;
}

using namespace PyROOT;

namespace {

static void AddToGlobalScope(const char *label, const char * /* hdr */, TObject *obj, Cppyy::TCppType_t klass)
{
   // Bind the given object with the given class in the global scope with the
   // given label for its reference.
   PyModule_AddObject(gRootModule, const_cast<char *>(label), CPyCppyy::BindCppObjectNoCast(obj, klass));
}

} // unnamed namespace

static TMemoryRegulator &GetMemoryRegulator()
{
   static TMemoryRegulator m;
   return m;
}

void PyROOT::Init()
{
   // Initialize and acquire the GIL to allow for threading in ROOT
   PyEval_InitThreads();

   // Memory management
   gROOT->GetListOfCleanups()->Add(&GetMemoryRegulator());

   // Bind ROOT globals that will be needed in ROOT.py
   AddToGlobalScope("gROOT", "TROOT.h", gROOT, Cppyy::GetScope(gROOT->IsA()->GetName()));
   AddToGlobalScope("gSystem", "TSystem.h", gSystem, Cppyy::GetScope(gSystem->IsA()->GetName()));
   AddToGlobalScope("gInterpreter", "TInterpreter.h", gInterpreter, Cppyy::GetScope(gInterpreter->IsA()->GetName()));
}

PyObject *PyROOT::ClearProxiedObjects(PyObject * /* self */, PyObject * /* args */)
{
   // Delete all memory-regulated objects
   GetMemoryRegulator().ClearProxiedObjects();
   Py_RETURN_NONE;
}
