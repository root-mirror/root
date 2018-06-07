#ifndef CPYCPPYY_TPyException
#define CPYCPPYY_TPyException

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TPyException                                                             //
//                                                                          //
// Purpose: A C++ exception class for throwing python exceptions            //
//          through C++ code.                                               //
// Created: Apr, 2004, Scott Snyder, from the version in D0's python_util.  //
//                                                                          //
// The situation is:                                                        //
//   - We're calling C++ code from python.                                  //
//   - The C++ code can call back to python.                                //
//   - What to do then if the python callback throws an exception?          //
//                                                                          //
// We need to get the control flow back to where CPyCppyy calls C++.        //
// To do that we throw a TPyException.                                      //
// We can then catch this exception when we do the C++ call.                //
//                                                                          //
// Note that we don't need to save any state in the exception -- it's       //
// already in the python error info variables.                              //
// (??? Actually, if the program is multithreaded, this is dangerous        //
// if the code has released and reacquired the lock along the call chain.   //
// Punt on this for now, though.)                                           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

// Standard
#include <exception>


namespace CPyCppyy {

class TPyException : public std::exception {
public:
// default constructor
    TPyException();

// destructor
    virtual ~TPyException() noexcept;

// give reason for raised exception
    virtual const char* what() const noexcept;
};

} // namespace CPyCppyy

#endif // !CPYCPPYY_TPyException
