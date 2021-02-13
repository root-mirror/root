// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Apr 2004

#ifndef PYROOT_TPYBUFFERFACTORY_H
#define PYROOT_TPYBUFFERFACTORY_H


namespace PyROOT {

/** Factory for python buffers of non-string type
      @author  WLAV
      @date    10/28/2004
      @version 1.5
*/

class TPyBufferFactory {
public:
   static TPyBufferFactory* Instance();

   PyObject* PyBuffer_FromMemory( Bool_t* buf, Py_ssize_t size = -1 );
   PyObject* PyBuffer_FromMemory( Bool_t* buf, PyObject* sizeCallback );
   PyObject* PyBuffer_FromMemory( Char_t* buf, Py_ssize_t size = -1 );
   PyObject* PyBuffer_FromMemory( Char_t* buf, PyObject* sizeCallback );
   PyObject* PyBuffer_FromMemory( UChar_t* buf, Py_ssize_t size = -1 );
   PyObject* PyBuffer_FromMemory( UChar_t* buf, PyObject* sizeCallback );
   PyObject* PyBuffer_FromMemory( Short_t* buf, Py_ssize_t size = -1 );
   PyObject* PyBuffer_FromMemory( Short_t* buf, PyObject* sizeCallback );
   PyObject* PyBuffer_FromMemory( UShort_t* buf, Py_ssize_t size = -1 );
   PyObject* PyBuffer_FromMemory( UShort_t* buf, PyObject* sizeCallback );
   PyObject* PyBuffer_FromMemory( Int_t* buf, Py_ssize_t size = -1 );
   PyObject* PyBuffer_FromMemory( Int_t* buf, PyObject* sizeCallback );
   PyObject* PyBuffer_FromMemory( UInt_t* buf, Py_ssize_t size = -1 );
   PyObject* PyBuffer_FromMemory( UInt_t* buf, PyObject* sizeCallback );
   PyObject* PyBuffer_FromMemory( Long_t* buf, Py_ssize_t size = -1 );
   PyObject* PyBuffer_FromMemory( Long_t* buf, PyObject* sizeCallback );
   PyObject* PyBuffer_FromMemory( ULong_t* buf, Py_ssize_t size = -1 );
   PyObject* PyBuffer_FromMemory( ULong_t* buf, PyObject* sizeCallback );
   PyObject* PyBuffer_FromMemory( Float_t* buf, Py_ssize_t size = -1 );
   PyObject* PyBuffer_FromMemory( Float_t* buf, PyObject* sizeCallback );
   PyObject* PyBuffer_FromMemory( Double_t* buf, Py_ssize_t size = -1 );
   PyObject* PyBuffer_FromMemory( Double_t* buf, PyObject* sizeCallback );

protected:
   TPyBufferFactory();
   ~TPyBufferFactory();
};

typedef TPyBufferFactory BufFac_t;

} // namespace PyROOT


#endif // !PYROOT_TPYBUFFERFACTORY_H
