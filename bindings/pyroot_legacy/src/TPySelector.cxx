// Author: Wim Lavrijsen   March 2008

// Bindings
#include "PyROOT.h"
#include "TPySelector.h"
#include "ObjectProxy.h"
#include "MethodProxy.h"
#include "RootWrapper.h"
#include "TMemoryRegulator.h"

//- ROOT
#include "TPython.h"
#include "TString.h"

//______________________________________________________________________________
//                      Python equivalent PROOF base class
//                      ==================================
//
// The problem with deriving a python class from a PyROOT bound class and then
// handing it back to a C++ framework, is that the virtual function dispatching
// of C++ is completely oblivious to the methods overridden in python. To work
// within the PROOF (C++) framework, a python class should derive from the class
// TPySelector. This class provides the proper overrides on the C++ side, and
// then forwards them, as apropriate, to python.
//
// This is an example set of scripts:
//
// ### PROOF running script, very close to equivalent .C (prooftest.py)
// import time
// from ROOT import *
//
// dataset = TDSet( 'TTree', 'h42' )
// dataset.Add( 'root:// .... myfile.root' )
//
// proof = TProof.Open('')
// time.sleep(1)                     # needed for GUI to settle
// print dataset.Process( 'TPySelector', 'aapje' )
// ### EOF
//
// ### selector module (aapje.py, name has to match as per above)
// from ROOT import TPySelector
//
// class MyPySelector( TPySelector ):
//    def Begin( self ):
//       print 'py: beginning'
//
//    def SlaveBegin( self, tree ):
//       print 'py: slave beginning'
//
//    def Process( self, entry ):
//       if self.fChain.GetEntry( entry ) <= 0:
//          return 0
//       print 'py: processing', self.fChain.MyVar
//       return 1
//
//   def SlaveTerminate( self ):
//       print 'py: slave terminating'
//
//   def Terminate( self ):
//       print 'py: terminating'
// ### EOF


//- data ---------------------------------------------------------------------
ClassImp(TPySelector);


//- private helpers ----------------------------------------------------------
void TPySelector::SetupPySelf()
{
// Install the python side identity of the TPySelector.
   if ( fPySelf && fPySelf != Py_None )
      return;                      // already created ...

// split option as needed for the module part and the (optional) user part
   std::string opt = GetOption();
   std::string::size_type pos = opt.find( '#' );
   std::string module = opt.substr( 0, pos );
   std::string user = (pos == std::string::npos) ? "" : opt.substr( pos+1, std::string::npos );

   TString impst = TString::Format( "import %s", module.c_str() );

// reset user option
   SetOption( user.c_str() );

// use TPython to ensure that the interpreter is initialized
   if ( ! TPython::Exec( (const char*)impst ) ) {
      Abort( "failed to load provided python module" );  // Exec already printed error trace
      return;
   }

// get the TPySelector python class
   PyObject* tpysel = PyObject_GetAttrString(
      PyImport_AddModule( const_cast< char* >( "libPyROOT" ) ),
      const_cast< char* >( "TPySelector" ) );

// get handle to the module
   PyObject* pymod = PyImport_AddModule( const_cast< char* >( module.c_str() ) );

// get the module dictionary to loop over
   PyObject* dict = PyModule_GetDict( pymod );
   Py_INCREF( dict );

// locate the TSelector derived class
   PyObject* allvalues = PyDict_Values( dict );

   PyObject* pyclass = 0;
   for ( int i = 0; i < PyList_GET_SIZE( allvalues ); ++i ) {
      PyObject* value = PyList_GET_ITEM( allvalues, i );
      Py_INCREF( value );

      if ( PyType_Check( value ) && PyObject_IsSubclass( value, tpysel ) ) {
         if ( PyObject_RichCompareBool( value, tpysel, Py_NE ) ) {   // i.e., if not equal
            pyclass = value;
            break;
         }
      }

      Py_DECREF( value );
   }

   Py_DECREF( allvalues );
   Py_DECREF( dict );
   Py_DECREF( tpysel );

   if ( ! pyclass ) {
      Abort( "no TSelector derived class available in provided module" );
      return;
   }

   PyObject* args = PyTuple_New( 0 );
   PyObject* self = PyObject_Call( pyclass, args, 0 );
   Py_DECREF( args );
   Py_DECREF( pyclass );

// final check before declaring success ...
   if ( ! self || ! PyROOT::ObjectProxy_Check( self ) ) {
      if ( ! PyErr_Occurred() )
         PyErr_SetString( PyExc_RuntimeError, "could not create python selector" );
      Py_XDECREF( self );
      Abort( 0 );
      return;
   }

// steal reference to new self, since the deletion will come from the C++ side
   Py_XDECREF( fPySelf );
   fPySelf = self;

// inject ourselves into the base of self; destroy old identity if need be (which happens
// if the user calls the default ctor unnecessarily)
   TPySelector* oldselector = (TPySelector*)((PyROOT::ObjectProxy*)fPySelf)->fObject;
   ((PyROOT::ObjectProxy*)fPySelf)->fObject = this;
   if ( oldselector ) {
      PyROOT::TMemoryRegulator::UnregisterObject( oldselector );
      delete oldselector;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Forward <method> to python.

PyObject* TPySelector::CallSelf( const char* method, PyObject* pyobject )
{
   if ( ! fPySelf || fPySelf == Py_None ) {
      Py_INCREF( Py_None );
      return Py_None;
   }

   PyObject* result = 0;

// get the named method and check for python side overload by not accepting the
// binding's methodproxy
   PyObject* pymethod = PyObject_GetAttrString( fPySelf, const_cast< char* >( method ) );
   if ( ! PyROOT::MethodProxy_CheckExact( pymethod ) ) {
      if ( pyobject )
         result = PyObject_CallFunction( pymethod, const_cast< char* >( "O" ), pyobject );
      else
         result = PyObject_CallFunction( pymethod, const_cast< char* >( "" ) );
   } else {
   // silently ignore if method not overridden (note that the above can't lead
   // to a python exception, since this (TPySelector) class contains the method
   // so it is always to be found)
      Py_INCREF( Py_None );
      result = Py_None;
   }

   Py_XDECREF( pymethod );

   if ( ! result )
      Abort( 0 );

   return result;
}


//- constructors/destructor --------------------------------------------------
TPySelector::TPySelector( TTree*, PyObject* self ) : fChain( 0 ), fPySelf( 0 )
{
// Construct a TSelector derived with <self> as the underlying, which is
// generally 0 to start out with in the current PROOF framework.
   if ( self ) {
   // steal reference as this is us, as seen from python
      fPySelf = self;
   } else {
      Py_INCREF( Py_None );        // using None allows clearer diagnostics
      fPySelf = Py_None;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. Only deref if still holding on to Py_None (circular otherwise).

TPySelector::~TPySelector()
{
   if ( fPySelf == Py_None ) {
      Py_DECREF( fPySelf );
   }
}


//- public functions ---------------------------------------------------------
Int_t TPySelector::Version() const {
// Return version number of this selector. First forward; if not overridden, then
// yield an obvious "undefined" number,
   PyObject* result = const_cast< TPySelector* >( this )->CallSelf( "Version" );
   if ( result && result != Py_None ) {
      Int_t ires = (Int_t)PyLong_AsLong( result );
      Py_DECREF( result );
      return ires;
   } else if ( result == Py_None ) {
      Py_DECREF( result );
   }
   return -99;
}

////////////////////////////////////////////////////////////////////////////////
/// Boilerplate get entry; same as for generated code; not forwarded.

Int_t TPySelector::GetEntry( Long64_t entry, Int_t getall )
{
   return fChain ? fChain->GetTree()->GetEntry( entry, getall ) : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize with the current tree to be used; not forwarded (may be called
/// multiple times, and is called from Begin() and SlaveBegin() ).

void TPySelector::Init( TTree* tree )
{
   if ( ! tree )
      return;

// set the fChain beforehand so that the python side may correct if needed
   fChain = tree;

// forward call
   PyObject* pytree = PyROOT::BindCppObject( (void*)tree, tree->IsA()->GetName() );
   PyObject* result = CallSelf( "Init", pytree );
   Py_DECREF( pytree );

   if ( ! result )
      Abort( 0 );

   Py_XDECREF( result );
}

////////////////////////////////////////////////////////////////////////////////
/// Forward call to derived Notify() if available.

Bool_t TPySelector::Notify()
{
   PyObject* result = CallSelf( "Notify" );

   if ( ! result )
      Abort( 0 );

   Py_XDECREF( result );

// by default, return kTRUE, b/c the Abort will stop the processing anyway on
// a real error, so if we get here it usually means that there is no Notify()
// override on the python side of things
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// First function called, and used to setup the python self; forward call.

void TPySelector::Begin( TTree* )
{
   SetupPySelf();

// As per the generated code: the tree argument is deprecated (on PROOF 0 is
// passed), and hence not forwarded.
   PyObject* result = CallSelf( "Begin" );

   if ( ! result )
      Abort( 0 );

   Py_XDECREF( result );
}

////////////////////////////////////////////////////////////////////////////////
/// First function called on worker node, needs to make sure python self is setup,
/// then store the tree to be used, initialize client, and forward call.

void TPySelector::SlaveBegin( TTree* tree )
{
   SetupPySelf();
   Init( tree );

   PyObject* result = 0;
   if ( tree ) {
      PyObject* pytree = PyROOT::BindCppObject( (void*)tree, tree->IsA()->GetName() );
      result = CallSelf( "SlaveBegin", pytree );
      Py_DECREF( pytree );
   } else {
      result = CallSelf( "SlaveBegin", Py_None );
   }

   if ( ! result )
      Abort( 0 );

   Py_XDECREF( result );
}

////////////////////////////////////////////////////////////////////////////////
/// Actual processing; call is forwarded to python self.

Bool_t TPySelector::Process( Long64_t entry )
{
   if ( ! fPySelf || fPySelf == Py_None ) {
   // would like to set a python error, but can't risk that in case of a
   // configuration problem, as it would be absorbed ...

   // simply returning kFALSE will not stop processing; need to set abort
      Abort( "no python selector instance available" );
      return kFALSE;
   }

   PyObject* result = PyObject_CallMethod( fPySelf,
      const_cast< char* >( "Process" ), const_cast< char* >( "L" ), entry );
   if ( ! result ) {
      Abort( 0 );
      return kFALSE;
   }

   Bool_t bresult = (Bool_t)PyLong_AsLong( result );
   Py_DECREF( result );
   return bresult;
}

////////////////////////////////////////////////////////////////////////////////
/// End of client; call is forwarded to python self.

void TPySelector::SlaveTerminate()
{
   PyObject* result = CallSelf( "SlaveTerminate" );

   if ( ! result )
      Abort( 0 );

   Py_XDECREF( result );
}

////////////////////////////////////////////////////////////////////////////////
/// End of job; call is forwarded to python self.

void TPySelector::Terminate()
{
   PyObject* result = CallSelf( "Terminate" );

   if ( ! result )
      Abort( 0 );

   Py_XDECREF( result );
}

////////////////////////////////////////////////////////////////////////////////
/// If no 'why' given, read from python error

void TPySelector::Abort( const char* why, EAbort what )
{
   if ( ! why && PyErr_Occurred() ) {
      PyObject *pytype = 0, *pyvalue = 0, *pytrace = 0;
      PyErr_Fetch( &pytype, &pyvalue, &pytrace );

   // abort is delayed (done at end of loop, message is current)
      PyObject* pystr = PyObject_Str( pyvalue );
      Abort( PyROOT_PyUnicode_AsString( pystr ), what );
      Py_DECREF( pystr );

      PyErr_Restore( pytype, pyvalue, pytrace );
   } else
      TSelector::Abort( why ? why : "", what );
}
