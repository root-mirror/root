// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen   Aug 2013

#ifndef ROOT_TPyArg
#define ROOT_TPyArg

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TPyArg                                                                   //
//                                                                          //
// Morphing argument type from evaluating python expressions.               //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


// ROOT
#include "Rtypes.h"

// Python
struct _object;
typedef _object PyObject;

// Standard
#include <vector>


class TPyArg {
public:
// converting constructors
   TPyArg( PyObject* );
   TPyArg( Int_t );
   TPyArg( Long_t );
   TPyArg( Double_t );
   TPyArg( const char* );

   TPyArg( const TPyArg& );
   TPyArg& operator=( const TPyArg& );
   virtual ~TPyArg();

// "extractor"
   operator PyObject*() const;

// constructor and generic dispatch
   static void CallConstructor( PyObject*& pyself, PyObject* pyclass, const std::vector<TPyArg>& args );
   static void CallConstructor( PyObject*& pyself, PyObject* pyclass );   // default ctor
   static PyObject* CallMethod( PyObject* pymeth, const std::vector<TPyArg>& args );
   static void CallDestructor( PyObject*& pyself, PyObject* pymeth, const std::vector<TPyArg>& args );
   static void CallDestructor( PyObject*& pyself );

   ClassDef(TPyArg,1)   //Python morphing argument type

private:
   mutable PyObject* fPyObject;        //! converted C++ value as python object
};

#endif
