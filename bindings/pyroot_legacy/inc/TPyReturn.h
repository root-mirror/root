// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen   May 2004

#ifndef ROOT_TPyReturn
#define ROOT_TPyReturn

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TPyReturn                                                                //
//                                                                          //
// Morphing return type from evaluating python expressions.                 //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


// ROOT
#include "Rtypes.h"

// Python
struct _object;
typedef _object PyObject;


class TPyReturn {
public:
   TPyReturn();
   TPyReturn( PyObject* pyobject );
   TPyReturn( const TPyReturn& );
   TPyReturn& operator=( const TPyReturn& );
   virtual ~TPyReturn();

// conversions to standard types, may fail if unconvertible
   operator char*() const;
   operator const char*() const;
   operator Char_t() const;

   operator Long_t() const;
   operator Int_t() const { return (Int_t)operator Long_t(); }
   operator Short_t() const { return (Short_t)operator Long_t(); }

   operator ULong_t() const;
   operator UInt_t() const { return (UInt_t)operator ULong_t(); }
   operator UShort_t() const { return (UShort_t)operator ULong_t(); }

   operator Double_t() const;
   operator Float_t() const { return (Float_t)operator Double_t(); }

// used for both TObject and PyObject conversions
   operator void*() const;

   template<class T>
   operator T*() const { return (T*)(void*)*this; }

// used strictly for PyObject conversions
   operator PyObject*() const;

   ClassDef(TPyReturn,1)   //Python morphing return object

private:
   PyObject* fPyObject;            //! actual python object
};

#endif
