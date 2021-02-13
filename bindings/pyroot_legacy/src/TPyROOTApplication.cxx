// Author: Wim Lavrijsen, February 2006

// Bindings
#include "PyROOT.h"
#include "TPyROOTApplication.h"
#include "Utility.h"

// ROOT
#include "TROOT.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TBenchmark.h"
#include "TStyle.h"
#include "TError.h"
#include "Getline.h"
#ifdef R__WIN32
#include "TVirtualX.h"
#endif

// Standard
#include <string.h>


//______________________________________________________________________________
//                   Setup interactive application for python
//                   ========================================
//
// The TPyROOTApplication sets up the nuts and bolts for interactive ROOT use
// from python, closely following TRint. Note that not everything is done here,
// some bits (such as e.g. the use of exception hook for shell escapes) are more
// easily done in python and you'll thus find them ROOT.py
//
// The intended use of this class is from python only. It is used by default in
// ROOT.py, so if you do not want to have a TApplication derived object created
// for you, you'll need to load libPyROOT.so instead.
//
// The static InitXYZ functions are used in conjunction with TPyROOTApplication
// in ROOT.py, but they can be used independently.
//
// NOTE: This class will receive the command line arguments from sys.argv. A
// distinction between arguments for TApplication and user arguments can be
// made by using "-" or "--" as a separator on the command line.


//- data ---------------------------------------------------------------------
ClassImp(PyROOT::TPyROOTApplication);


//- constructors/destructor --------------------------------------------------
PyROOT::TPyROOTApplication::TPyROOTApplication(
   const char* acn, int* argc, char** argv, Bool_t /*bLoadLibs*/ ) :
      TApplication( acn, argc, argv )
{
// The following code is redundant with ROOT6 and the PCH: the headers are
// available to the interpreter.
// // Create a TApplication derived for use with interactive ROOT from python. A
// // set of standard, often used libs is loaded if bLoadLibs is true (default).
//
//    if ( bLoadLibs )   // note that this section could be programmed in python
//    {
//    // follow TRint to minimize differences with root.exe (note: changed <pair>
//    // to <utility> for Cling, which is correct)
//       ProcessLine( "#include <iostream>", kTRUE );
//       ProcessLine( "#include <string>",   kTRUE ); // for std::string iostream.
//       ProcessLine( "#include <vector>",   kTRUE ); // needed because they're used within the
//       ProcessLine( "#include <utility>",  kTRUE ); //  core ROOT dicts and CINT won't be able
//                                                    //  to properly unload these files
//    }

#ifdef WIN32
   // switch win32 proxy main thread id
   if (gVirtualX)
      ProcessLine("((TGWin32 *)gVirtualX)->SetUserThreadId(0);", kTRUE);
#endif

// save current interpreter context
   gInterpreter->SaveContext();
   gInterpreter->SaveGlobalsContext();

// prevent crashes on accessing history
   Gl_histinit( (char*)"-" );

// prevent ROOT from exiting python
   SetReturnFromRun( kTRUE );
}


//- static public members ----------------------------------------------------
Bool_t PyROOT::TPyROOTApplication::CreatePyROOTApplication( Bool_t bLoadLibs )
{
// Create a TPyROOTApplication. Returns false if gApplication is not null.

   if ( ! gApplication ) {
   // retrieve arg list from python, translate to raw C, pass on
      PyObject* argl = PySys_GetObject( const_cast< char* >( "argv" ) );

      int argc = 1;
      if ( argl && 0 < PyList_Size( argl ) ) argc = (int)PyList_GET_SIZE( argl );
      char** argv = new char*[ argc ];
      for ( int i = 1; i < argc; ++i ) {
         char* argi = const_cast< char* >( PyROOT_PyUnicode_AsString( PyList_GET_ITEM( argl, i ) ) );
         if ( strcmp( argi, "-" ) == 0 || strcmp( argi, "--" ) == 0 ) {
         // stop collecting options, the remaining are for the python script
            argc = i;    // includes program name
            break;
         }
         argv[ i ] = argi;
      }
#if PY_VERSION_HEX < 0x03000000
      if ( Py_GetProgramName() && strlen( Py_GetProgramName() ) != 0 )
         argv[ 0 ] = Py_GetProgramName();
      else
         argv[ 0 ] = (char*)"python";
#else
// TODO: convert the wchar_t*
      argv[ 0 ] = (char*)"python";
#endif

      gApplication = new TPyROOTApplication( "PyROOT", &argc, argv, bLoadLibs );
      delete[] argv;     // TApplication ctor has copied argv, so done with it

      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Setup the basic ROOT globals gBenchmark, gStyle, gProgname, if not already
/// set. Always returns true.

Bool_t PyROOT::TPyROOTApplication::InitROOTGlobals()
{
   if ( ! gBenchmark ) gBenchmark = new TBenchmark();
   if ( ! gStyle ) gStyle = new TStyle();

   if ( ! gProgName )              // should have been set by TApplication
#if PY_VERSION_HEX < 0x03000000
      gSystem->SetProgname( Py_GetProgramName() );
#else
// TODO: convert the wchar_t*
      gSystem->SetProgname( "python" );
#endif

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Install ROOT message handler which will turn ROOT error message into
/// python exceptions. Always returns true.

Bool_t PyROOT::TPyROOTApplication::InitROOTMessageCallback()
{
   SetErrorHandler( (ErrorHandlerFunc_t)&Utility::ErrMsgHandler );
   return kTRUE;
}
