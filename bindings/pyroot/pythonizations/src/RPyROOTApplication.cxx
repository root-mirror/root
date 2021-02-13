// Author: Enric Tejedor CERN  04/2019
// Original PyROOT code by Wim Lavrijsen, LBL

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// Bindings
#include "Python.h"
#include "CPyCppyy.h"
#include "RPyROOTApplication.h"

// ROOT
#include "TInterpreter.h"
#include "TSystem.h"
#include "TBenchmark.h"
#include "TStyle.h"
#include "TError.h"
#include "Getline.h"
#include "TVirtualMutex.h"

////////////////////////////////////////////////////////////////////////////
/// \brief Create an RPyROOTApplication.
/// \param[in] ignoreCmdLineOpts True if Python command line options should
///            be ignored.
/// \return false if gApplication is not null, true otherwise.
///
/// If ignoreCmdLineOpts is false, this method processes the command line
/// arguments from sys.argv. A distinction between arguments for
/// TApplication and user arguments can be made by using "-" or "--" as a
/// separator on the command line.
///
/// For example, to enable batch mode from the command line:
/// > python script_name.py -b -- user_arg1 ... user_argn
/// or, if the user script receives no arguments:
/// > python script_name.py -b
bool PyROOT::RPyROOTApplication::CreateApplication(int ignoreCmdLineOpts)
{
   if (!gApplication) {
      int argc = 1;
      char **argv = nullptr;

      if (ignoreCmdLineOpts) {
         argv = new char *[argc];
      } else {
         // Retrieve sys.argv list from Python
         PyObject *argl = PySys_GetObject(const_cast<char *>("argv"));

         if (argl && 0 < PyList_Size(argl))
            argc = (int)PyList_GET_SIZE(argl);

         argv = new char *[argc];
         for (int i = 1; i < argc; ++i) {
            char *argi = const_cast<char *>(CPyCppyy_PyText_AsString(PyList_GET_ITEM(argl, i)));
            if (strcmp(argi, "-") == 0 || strcmp(argi, "--") == 0) {
               // Stop collecting options, the remaining are for the Python script
               argc = i; // includes program name
               break;
            }
            argv[i] = argi;
         }
      }

#if PY_VERSION_HEX < 0x03000000
      if (Py_GetProgramName() && strlen(Py_GetProgramName()) != 0)
         argv[0] = Py_GetProgramName();
      else
         argv[0] = (char *)"python";
#else
      argv[0] = (char *)"python";
#endif

      gApplication = new RPyROOTApplication("PyROOT", &argc, argv);
      delete[] argv; // TApplication ctor has copied argv, so done with it

      return true;
   }

   return false;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Setup the basic ROOT globals gBenchmark, gStyle and gProgname,
/// if not already set.
void PyROOT::RPyROOTApplication::InitROOTGlobals()
{
   if (!gBenchmark)
      gBenchmark = new TBenchmark();
   if (!gStyle)
      gStyle = new TStyle();

   if (!gProgName) // should have been set by TApplication
#if PY_VERSION_HEX < 0x03000000
      gSystem->SetProgname(Py_GetProgramName());
#else
      gSystem->SetProgname("python");
#endif
}

////////////////////////////////////////////////////////////////////////////
/// \brief Translate ROOT error/warning to Python.
static void ErrMsgHandler(int level, Bool_t abort, const char *location, const char *msg)
{
   // Initialization from gEnv (the default handler will return w/o msg b/c level too low)
   if (gErrorIgnoreLevel == kUnset)
      ::DefaultErrorHandler(kUnset - 1, kFALSE, "", "");

   if (level < gErrorIgnoreLevel)
      return;

   // Turn warnings into Python warnings
   if (level >= kError) {
      ::DefaultErrorHandler(level, abort, location, msg);
   } else if (level >= kWarning) {
      static const char *emptyString = "";
      if (!location)
         location = emptyString;
      // This warning might be triggered while holding the ROOT lock, while
      // some other thread is holding the GIL and waiting for the ROOT lock.
      // That will trigger a deadlock.
      // So if ROOT is in MT mode, use ROOT's error handler that doesn't take
      // the GIL.
      if (!gGlobalMutex) {
         // Either printout or raise exception, depending on user settings
         PyErr_WarnExplicit(NULL, (char *)msg, (char *)location, 0, (char *)"ROOT", NULL);
      } else {
         ::DefaultErrorHandler(level, abort, location, msg);
      }
   } else {
      ::DefaultErrorHandler(level, abort, location, msg);
   }
}

////////////////////////////////////////////////////////////////////////////
/// \brief Install the ROOT message handler which will turn ROOT error
/// messages into Python exceptions.
void PyROOT::RPyROOTApplication::InitROOTMessageCallback()
{
   SetErrorHandler((ErrorHandlerFunc_t)&ErrMsgHandler);
}

////////////////////////////////////////////////////////////////////////////
/// \brief Initialize an RPyROOTApplication.
/// \param[in] self Always null, since this is a module function.
/// \param[in] ignoreCmdLineOpts True if Python command line options should
///            be ignored.
PyObject *PyROOT::RPyROOTApplication::InitApplication(PyObject * /*self*/, PyObject *args)
{
   int argc = PyTuple_GET_SIZE(args);
   if (argc == 1) { 
      PyObject *ignoreCmdLineOpts = PyTuple_GetItem(args, 0); 
      
      if (!PyBool_Check(ignoreCmdLineOpts)) {
         PyErr_SetString(PyExc_TypeError, "Expected boolean type as argument.");
         return nullptr;
      }

      if (CreateApplication(PyObject_IsTrue(ignoreCmdLineOpts))) {
         InitROOTGlobals();
         InitROOTMessageCallback();
      }
   } else {
      PyErr_Format(PyExc_TypeError, "Expected 1 argument, %d passed.", argc);
      return nullptr;
   }

   Py_RETURN_NONE;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Construct a TApplication for PyROOT.
/// \param[in] acn Application class name.
/// \param[in] argc Number of arguments.
/// \param[in] argv Arguments.
PyROOT::RPyROOTApplication::RPyROOTApplication(const char *acn, int *argc, char **argv) : TApplication(acn, argc, argv)
{
   // Save current interpreter context
   gInterpreter->SaveContext();
   gInterpreter->SaveGlobalsContext();

   // Prevent crashes on accessing history
   Gl_histinit((char *)"-");

   // Prevent ROOT from exiting python
   SetReturnFromRun(true);
}

namespace {
static int (*sOldInputHook)() = nullptr;
static PyThreadState *sInputHookEventThreadState = nullptr;

static int EventInputHook()
{
   // This method is supposed to be called from CPython's command line and
   // drives the GUI
   PyEval_RestoreThread(sInputHookEventThreadState);
   gSystem->ProcessEvents();
   PyEval_SaveThread();

   if (sOldInputHook)
      return sOldInputHook();

   return 0;
}

} // unnamed namespace

////////////////////////////////////////////////////////////////////////////
/// \brief Install a method hook for sending events to the GUI.
/// \param[in] self Always null, since this is a module function.
/// \param[in] args Pointer to an empty Python tuple.
PyObject *PyROOT::RPyROOTApplication::InstallGUIEventInputHook(PyObject * /* self */, PyObject * /* args */)
{
   if (PyOS_InputHook && PyOS_InputHook != &EventInputHook)
      sOldInputHook = PyOS_InputHook;

   sInputHookEventThreadState = PyThreadState_Get();

   PyOS_InputHook = (int (*)()) & EventInputHook;

   Py_RETURN_NONE;
}
