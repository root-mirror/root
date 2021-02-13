// Author: Wim Lavrijsen, Nov 2008

// Bindings
#include "PyROOT.h"
#include "PyStrings.h"


//- data _____________________________________________________________________
PyObject* PyROOT::PyStrings::gBases = 0;
PyObject* PyROOT::PyStrings::gBase = 0;
PyObject* PyROOT::PyStrings::gClass = 0;
PyObject* PyROOT::PyStrings::gCppEq = 0;
PyObject* PyROOT::PyStrings::gCppNe = 0;
PyObject* PyROOT::PyStrings::gDeref = 0;
PyObject* PyROOT::PyStrings::gDict = 0;
PyObject* PyROOT::PyStrings::gEmptyString = 0;
PyObject* PyROOT::PyStrings::gEq = 0;
PyObject* PyROOT::PyStrings::gFollow = 0;
PyObject* PyROOT::PyStrings::gGetItem = 0;
PyObject* PyROOT::PyStrings::gInit = 0;
PyObject* PyROOT::PyStrings::gIter = 0;
PyObject* PyROOT::PyStrings::gLen = 0;
PyObject* PyROOT::PyStrings::gLifeLine = 0;
PyObject* PyROOT::PyStrings::gModule = 0;
PyObject* PyROOT::PyStrings::gMRO = 0;
PyObject* PyROOT::PyStrings::gName = 0;
PyObject* PyROOT::PyStrings::gCppName = 0;
PyObject* PyROOT::PyStrings::gCppNameNew = 0;
PyObject* PyROOT::PyStrings::gNe = 0;
PyObject* PyROOT::PyStrings::gTypeCode = 0;

PyObject* PyROOT::PyStrings::gAdd = 0;
PyObject* PyROOT::PyStrings::gSub = 0;
PyObject* PyROOT::PyStrings::gMul = 0;
PyObject* PyROOT::PyStrings::gDiv = 0;

PyObject* PyROOT::PyStrings::gAt = 0;
PyObject* PyROOT::PyStrings::gBegin = 0;
PyObject* PyROOT::PyStrings::gEnd = 0;
PyObject* PyROOT::PyStrings::gFirst = 0;
PyObject* PyROOT::PyStrings::gSecond = 0;
PyObject* PyROOT::PyStrings::gSize = 0;
PyObject* PyROOT::PyStrings::gGetSize = 0;
PyObject* PyROOT::PyStrings::ggetSize = 0;
PyObject* PyROOT::PyStrings::gTemplate = 0;
PyObject* PyROOT::PyStrings::gVectorAt = 0;

PyObject* PyROOT::PyStrings::gBranch = 0;
PyObject* PyROOT::PyStrings::gFitFCN = 0;
PyObject* PyROOT::PyStrings::gROOTns = 0;
PyObject* PyROOT::PyStrings::gSetBranchAddress = 0;
PyObject* PyROOT::PyStrings::gSetFCN = 0;
PyObject* PyROOT::PyStrings::gTClassDynCast = 0;


////////////////////////////////////////////////////////////////////////////////

#define PYROOT_INITIALIZE_STRING( var, str )                                 \
   if ( ! ( PyStrings::var = PyROOT_PyUnicode_InternFromString( (char*)#str ) ) )    \
      return kFALSE

Bool_t PyROOT::CreatePyStrings() {
// Build cache of commonly used python strings (the cache is python intern, so
// all strings are shared python-wide, not just in PyROOT).
   PYROOT_INITIALIZE_STRING( gBases, __bases__ );
   PYROOT_INITIALIZE_STRING( gBase, __base__ );
   PYROOT_INITIALIZE_STRING( gClass, __class__ );
   PYROOT_INITIALIZE_STRING( gCppEq, __cpp_eq__ );
   PYROOT_INITIALIZE_STRING( gCppNe, __cpp_ne__ );
   PYROOT_INITIALIZE_STRING( gDeref, __deref__ );
   PYROOT_INITIALIZE_STRING( gDict, __dict__ );
   if ( ! ( PyStrings::gEmptyString = PyROOT_PyUnicode_FromString( (char*)"" ) ) )
      return kFALSE;
   PYROOT_INITIALIZE_STRING( gEq, __eq__ );
   PYROOT_INITIALIZE_STRING( gFollow, __follow__ );
   PYROOT_INITIALIZE_STRING( gGetItem, __getitem__ );
   PYROOT_INITIALIZE_STRING( gInit, __init__ );
   PYROOT_INITIALIZE_STRING( gIter, __iter__ );
   PYROOT_INITIALIZE_STRING( gLen, __len__ );
   PYROOT_INITIALIZE_STRING( gLifeLine, __lifeline );
   PYROOT_INITIALIZE_STRING( gModule, __module__ );
   PYROOT_INITIALIZE_STRING( gMRO, __mro__ );
   PYROOT_INITIALIZE_STRING( gName, __name__ );
   PYROOT_INITIALIZE_STRING( gCppName, __cppname__ );
   PYROOT_INITIALIZE_STRING( gCppNameNew, __cpp_name__ );
   PYROOT_INITIALIZE_STRING( gNe, __ne__ );
   PYROOT_INITIALIZE_STRING( gTypeCode, typecode );

   PYROOT_INITIALIZE_STRING( gAdd, __add__ );
   PYROOT_INITIALIZE_STRING( gSub, __sub__ );
   PYROOT_INITIALIZE_STRING( gMul, __mul__ );
   PYROOT_INITIALIZE_STRING( gDiv, PYROOT__div__ );

   PYROOT_INITIALIZE_STRING( gAt, at );
   PYROOT_INITIALIZE_STRING( gBegin, begin );
   PYROOT_INITIALIZE_STRING( gEnd, end );
   PYROOT_INITIALIZE_STRING( gFirst, first );
   PYROOT_INITIALIZE_STRING( gSecond, second );
   PYROOT_INITIALIZE_STRING( gSize, size );
   PYROOT_INITIALIZE_STRING( gGetSize, GetSize );
   PYROOT_INITIALIZE_STRING( ggetSize, getSize );
   PYROOT_INITIALIZE_STRING( gTemplate, Template );
   PYROOT_INITIALIZE_STRING( gVectorAt, _vector__at );

   PYROOT_INITIALIZE_STRING( gBranch, Branch );
   PYROOT_INITIALIZE_STRING( gFitFCN, FitFCN );
   PYROOT_INITIALIZE_STRING( gROOTns, ROOT );
   PYROOT_INITIALIZE_STRING( gSetBranchAddress, SetBranchAddress );
   PYROOT_INITIALIZE_STRING( gSetFCN, SetFCN );
   PYROOT_INITIALIZE_STRING( gTClassDynCast, _TClass__DynamicCast );

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all cached python strings.

PyObject* PyROOT::DestroyPyStrings() {
   Py_DECREF( PyStrings::gBases ); PyStrings::gBases = 0;
   Py_DECREF( PyStrings::gBase ); PyStrings::gBase = 0;
   Py_DECREF( PyStrings::gClass ); PyStrings::gClass = 0;
   Py_DECREF( PyStrings::gCppEq ); PyStrings::gCppEq = 0;
   Py_DECREF( PyStrings::gCppNe ); PyStrings::gCppNe = 0;
   Py_DECREF( PyStrings::gDeref ); PyStrings::gDeref = 0;
   Py_DECREF( PyStrings::gDict ); PyStrings::gDict = 0;
   Py_DECREF( PyStrings::gEmptyString ); PyStrings::gEmptyString = 0;
   Py_DECREF( PyStrings::gEq ); PyStrings::gEq = 0;
   Py_DECREF( PyStrings::gFollow ); PyStrings::gFollow = 0;
   Py_DECREF( PyStrings::gGetItem ); PyStrings::gGetItem = 0;
   Py_DECREF( PyStrings::gInit ); PyStrings::gInit = 0;
   Py_DECREF( PyStrings::gIter ); PyStrings::gIter = 0;
   Py_DECREF( PyStrings::gLen ); PyStrings::gLen = 0;
   Py_DECREF( PyStrings::gLifeLine ); PyStrings::gLifeLine = 0;
   Py_DECREF( PyStrings::gModule ); PyStrings::gModule = 0;
   Py_DECREF( PyStrings::gMRO ); PyStrings::gMRO = 0;
   Py_DECREF( PyStrings::gName ); PyStrings::gName = 0;
   Py_DECREF( PyStrings::gCppName ); PyStrings::gCppName = 0;
   Py_DECREF( PyStrings::gCppNameNew ); PyStrings::gCppNameNew = 0;
   Py_DECREF( PyStrings::gNe ); PyStrings::gNe = 0;
   Py_DECREF( PyStrings::gTypeCode ); PyStrings::gTypeCode = 0;

   Py_DECREF( PyStrings::gAdd ); PyStrings::gAdd = 0;
   Py_DECREF( PyStrings::gSub ); PyStrings::gSub = 0;
   Py_DECREF( PyStrings::gMul ); PyStrings::gMul = 0;
   Py_DECREF( PyStrings::gDiv ); PyStrings::gDiv = 0;

   Py_DECREF( PyStrings::gAt ); PyStrings::gAt = 0;
   Py_DECREF( PyStrings::gBegin ); PyStrings::gBegin = 0;
   Py_DECREF( PyStrings::gEnd ); PyStrings::gEnd = 0;
   Py_DECREF( PyStrings::gFirst ); PyStrings::gFirst = 0;
   Py_DECREF( PyStrings::gSecond ); PyStrings::gSecond = 0;
   Py_DECREF( PyStrings::gSize ); PyStrings::gSize = 0;
   Py_DECREF( PyStrings::gGetSize ); PyStrings::gGetSize = 0;
   Py_DECREF( PyStrings::ggetSize ); PyStrings::ggetSize = 0;
   Py_DECREF( PyStrings::gTemplate ); PyStrings::gTemplate = 0;
   Py_DECREF( PyStrings::gVectorAt ); PyStrings::gVectorAt = 0;

   Py_DECREF( PyStrings::gBranch ); PyStrings::gBranch = 0;
   Py_DECREF( PyStrings::gFitFCN ); PyStrings::gFitFCN = 0;
   Py_DECREF( PyStrings::gROOTns ); PyStrings::gROOTns = 0;
   Py_DECREF( PyStrings::gSetBranchAddress ); PyStrings::gSetBranchAddress = 0;
   Py_DECREF( PyStrings::gSetFCN ); PyStrings::gSetFCN = 0;
   Py_DECREF( PyStrings::gTClassDynCast ); PyStrings::gTClassDynCast = 0;

   Py_INCREF( Py_None );
   return Py_None;
}
