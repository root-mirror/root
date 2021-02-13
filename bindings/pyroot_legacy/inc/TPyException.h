// @(#)root/pyroot:$Name:  $:$Id$
// Author: Scott Snyder, Apr 2004

#ifndef ROOT_TPyException
#define ROOT_TPyException

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TPyException                                                             //
//                                                                          //
// Purpose: A C++ exception class for throwing python exceptions            //
//          through C++ code.                                               //
// Created: Apr, 2004, sss, from the version in D0's python_util.           //
//                                                                          //
// The situation is:                                                        //
//   - We're calling ROOT C++ code from python.                             //
//   - The C++ code can call back to python.                                //
//   - What to do then if the python callback throws an exception?          //
//                                                                          //
// We need to get the control flow back to where PyROOT makes the ROOT call.//
// To do that we throw a TPyException.                                      //
// We can then catch this exception when we do the ROOT call.               //
//                                                                          //
// Note that we don't need to save any state in the exception -- it's       //
// already in the python error info variables.                              //
// (??? Actually, if the program is multithreaded, this is dangerous        //
// if the code has released and reacquired the lock along the call chain.   //
// Punt on this for now, though.)                                           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

// ROOT
#include "DllImport.h"
#include "Rtypes.h"

// Standard
#include <exception>


namespace PyROOT {

class TPyException : public std::exception {
public:
// default constructor
   TPyException();

// destructor
   virtual ~TPyException() noexcept;

// give reason for raised exception
   virtual const char* what() const noexcept;

   ClassDef(TPyException,0)   //C++ exception for throwing python exceptions
};

} // namespace PyROOT

#if defined(G__DICTIONARY) && defined(R__SOLARIS)
// Force the inclusion of rw/math.h
#include <limits>
// Work around interaction between a struct named exception in math.h,
// std::exception and the use of using namespace std;
#if (__SUNPRO_CC < 0x5050)
#define exception std::exception
#endif
#endif
#endif
