// Bindings
#include "CPyCppyy.h"
#define CPYCPPYY_INTERNAL 1
#include "CPyCppyy/PyException.h"
#undef CPYCPPYY_INTERNAL


//______________________________________________________________________________
//                 C++ exception for throwing python exceptions
//                 ============================================
// Purpose: A C++ exception class for throwing python exceptions
//          through C++ code.
// Created: Apr, 2004, Scott Snyder, from the version in D0's python_util.
//
// Note: Don't be tempted to declare the virtual functions defined here
//       as inline.
//       If you do, you may not be able to properly throw these
//       exceptions across shared libraries.


//- constructors/destructor --------------------------------------------------
CPyCppyy::PyException::PyException()
{
// default constructor
}

CPyCppyy::PyException::~PyException() noexcept
{
// destructor
}


//- public members -----------------------------------------------------------
const char* CPyCppyy::PyException::what() const noexcept
{
// Return reason for throwing this exception: a python exception was raised.
    return "python exception";
}
