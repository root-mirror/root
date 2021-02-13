// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Jan 2005
#ifndef PYROOT_EXECUTORS_H
#define PYROOT_EXECUTORS_H

// Bindings
#include "TCallContext.h"

// Standard
#include <string>


namespace PyROOT {

   class TExecutor {
   public:
      virtual ~TExecutor() {}
      virtual PyObject* Execute(
         Cppyy::TCppMethod_t, Cppyy::TCppObject_t, TCallContext* ) = 0;
   };

#define PYROOT_DECLARE_BASIC_EXECUTOR( name )                                 \
   class T##name##Executor : public TExecutor {                               \
   public:                                                                    \
      virtual PyObject* Execute(                                              \
         Cppyy::TCppMethod_t, Cppyy::TCppObject_t, TCallContext* );           \
   }

// executors for built-ins
   PYROOT_DECLARE_BASIC_EXECUTOR( Bool );
   PYROOT_DECLARE_BASIC_EXECUTOR( BoolConstRef );
   PYROOT_DECLARE_BASIC_EXECUTOR( Char );
   PYROOT_DECLARE_BASIC_EXECUTOR( CharConstRef );
   PYROOT_DECLARE_BASIC_EXECUTOR( UChar );
   PYROOT_DECLARE_BASIC_EXECUTOR( UCharConstRef );
   PYROOT_DECLARE_BASIC_EXECUTOR( Short );
   PYROOT_DECLARE_BASIC_EXECUTOR( Int );
   PYROOT_DECLARE_BASIC_EXECUTOR( Long );
   PYROOT_DECLARE_BASIC_EXECUTOR( ULong );
   PYROOT_DECLARE_BASIC_EXECUTOR( LongLong );
   PYROOT_DECLARE_BASIC_EXECUTOR( ULongLong );
   PYROOT_DECLARE_BASIC_EXECUTOR( Float );
   PYROOT_DECLARE_BASIC_EXECUTOR( Double );
   PYROOT_DECLARE_BASIC_EXECUTOR( LongDouble );
   PYROOT_DECLARE_BASIC_EXECUTOR( Void );
   PYROOT_DECLARE_BASIC_EXECUTOR( CString );

// pointer/array executors
   PYROOT_DECLARE_BASIC_EXECUTOR( VoidArray );
   PYROOT_DECLARE_BASIC_EXECUTOR( BoolArray );
   PYROOT_DECLARE_BASIC_EXECUTOR( ShortArray );
   PYROOT_DECLARE_BASIC_EXECUTOR( UShortArray );
   PYROOT_DECLARE_BASIC_EXECUTOR( CharArray );
   PYROOT_DECLARE_BASIC_EXECUTOR( UCharArray );
   PYROOT_DECLARE_BASIC_EXECUTOR( IntArray );
   PYROOT_DECLARE_BASIC_EXECUTOR( UIntArray );
   PYROOT_DECLARE_BASIC_EXECUTOR( LongArray );
   PYROOT_DECLARE_BASIC_EXECUTOR( ULongArray );
   PYROOT_DECLARE_BASIC_EXECUTOR( FloatArray );
   PYROOT_DECLARE_BASIC_EXECUTOR( DoubleArray );

// special cases
   PYROOT_DECLARE_BASIC_EXECUTOR( STLString );
   PYROOT_DECLARE_BASIC_EXECUTOR( TGlobal );

   class TCppObjectExecutor : public TExecutor {
   public:
      TCppObjectExecutor( Cppyy::TCppType_t klass ) : fClass( klass ) {}
      virtual PyObject* Execute(
         Cppyy::TCppMethod_t, Cppyy::TCppObject_t,TCallContext* );

   protected:
      Cppyy::TCppType_t fClass;
   };

   class TCppObjectByValueExecutor : public TCppObjectExecutor {
   public:
      using TCppObjectExecutor::TCppObjectExecutor;
      virtual PyObject* Execute(
         Cppyy::TCppMethod_t, Cppyy::TCppObject_t,TCallContext* );
   };

   class TRefExecutor : public TExecutor {
   public:
      TRefExecutor() : fAssignable( 0 ) {}

   public:
      virtual Bool_t SetAssignable( PyObject* );

   protected:
      PyObject* fAssignable;
   };

   PYROOT_DECLARE_BASIC_EXECUTOR( Constructor );
   PYROOT_DECLARE_BASIC_EXECUTOR( PyObject );

#define PYROOT_DECLARE_BASIC_REFEXECUTOR( name )                              \
   class T##name##RefExecutor : public TRefExecutor {                         \
   public:                                                                    \
      virtual PyObject* Execute(                                              \
         Cppyy::TCppMethod_t, Cppyy::TCppObject_t, TCallContext* );           \
   }

   PYROOT_DECLARE_BASIC_REFEXECUTOR( Bool );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( Char );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( UChar );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( Short );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( UShort );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( Int );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( UInt );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( Long );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( ULong );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( LongLong );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( ULongLong );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( Float );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( Double );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( LongDouble );
   PYROOT_DECLARE_BASIC_REFEXECUTOR( STLString );

// special cases
   class TCppObjectRefExecutor : public TRefExecutor {
   public:
      TCppObjectRefExecutor( Cppyy::TCppType_t klass ) : fClass( klass ) {}
      virtual PyObject* Execute(
         Cppyy::TCppMethod_t, Cppyy::TCppObject_t, TCallContext* );

   protected:
      Cppyy::TCppType_t fClass;
   };

   class TCppObjectPtrPtrExecutor : public TCppObjectExecutor {
   public:
      using TCppObjectExecutor::TCppObjectExecutor;
      virtual PyObject* Execute(
         Cppyy::TCppMethod_t, Cppyy::TCppObject_t, TCallContext* );
   };

   class TCppObjectPtrRefExecutor : public TCppObjectExecutor {
   public:
      using TCppObjectExecutor::TCppObjectExecutor;
      virtual PyObject* Execute(
         Cppyy::TCppMethod_t, Cppyy::TCppObject_t, TCallContext* );
   };

   class TCppObjectArrayExecutor : public TCppObjectExecutor {
   public:
      TCppObjectArrayExecutor( Cppyy::TCppType_t klass, Py_ssize_t array_size )
         : TCppObjectExecutor ( klass ), fArraySize( array_size ) {}
      virtual PyObject* Execute(
         Cppyy::TCppMethod_t, Cppyy::TCppObject_t, TCallContext* );

   protected:
      Py_ssize_t fArraySize;
   };

// smart pointer executors
   class TCppObjectBySmartPtrExecutor : public TExecutor {
   public:
      TCppObjectBySmartPtrExecutor( Cppyy::TCppType_t klass, Cppyy::TCppType_t rawPtrType,
         Cppyy::TCppMethod_t deref ) : fClass( klass ), fRawPtrType( rawPtrType ), fDereferencer( deref ) {}

      virtual PyObject* Execute(
         Cppyy::TCppMethod_t, Cppyy::TCppObject_t, TCallContext* );

   protected:
      Cppyy::TCppType_t   fClass;
      Cppyy::TCppType_t   fRawPtrType;
      Cppyy::TCppMethod_t fDereferencer;
   };

   class TCppObjectBySmartPtrPtrExecutor : public TCppObjectBySmartPtrExecutor {
   public:
      using TCppObjectBySmartPtrExecutor::TCppObjectBySmartPtrExecutor;

      virtual PyObject* Execute(
         Cppyy::TCppMethod_t, Cppyy::TCppObject_t, TCallContext* );
   };

   class TCppObjectBySmartPtrRefExecutor : public TRefExecutor {
   public:
     TCppObjectBySmartPtrRefExecutor( Cppyy::TCppType_t klass, Cppyy::TCppType_t rawPtrType,
        Cppyy::TCppMethod_t deref ) : fClass( klass ), fRawPtrType( rawPtrType ), fDereferencer( deref ) {}

      virtual PyObject* Execute(
         Cppyy::TCppMethod_t, Cppyy::TCppObject_t,TCallContext* );

   protected:
      Cppyy::TCppType_t fClass;
      Cppyy::TCppType_t fRawPtrType;
      Cppyy::TCppMethod_t fDereferencer;
   };

// create executor from fully qualified type
   TExecutor* CreateExecutor( const std::string& fullType,
                              Bool_t manage_smart_ptr = kTRUE );

} // namespace PyROOT

#endif // !PYROOT_EXECUTORS_H
