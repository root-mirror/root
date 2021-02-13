// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Jan 2005
#ifndef PYROOT_CONVERTERS_H
#define PYROOT_CONVERTERS_H

// ROOT
#include "TString.h"

// Standard
#include <limits.h>
#include <string>


namespace PyROOT {

   class ObjectProxy;
   struct TParameter;
   struct TCallContext;

   class TConverter {
   public:
      virtual ~TConverter() {}

   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 ) = 0;
      virtual PyObject* FromMemory( void* address );
      virtual Bool_t ToMemory( PyObject* value, void* address );
   };

#define PYROOT_DECLARE_BASIC_CONVERTER( name )                                \
   class T##name##Converter : public TConverter {                             \
   public:                                                                    \
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );\
      virtual PyObject* FromMemory( void* );                                  \
      virtual Bool_t ToMemory( PyObject*, void* );                            \
   };                                                                         \
                                                                              \
   class TConst##name##RefConverter : public TConverter {                     \
   public:                                                                    \
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );\
   }


#define PYROOT_DECLARE_BASIC_CONVERTER2( name, base )                         \
   class T##name##Converter : public T##base##Converter {                     \
   public:                                                                    \
      virtual PyObject* FromMemory( void* );                                  \
      virtual Bool_t ToMemory( PyObject*, void* );                            \
   };                                                                         \
                                                                              \
   class TConst##name##RefConverter : public TConverter {                     \
   public:                                                                    \
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );\
   }

#define PYROOT_DECLARE_REF_CONVERTER( name )                                  \
   class T##name##RefConverter : public TConverter {                          \
   public:                                                                    \
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );\
   };

#define PYROOT_DECLARE_ARRAY_CONVERTER( name )                                \
   class T##name##Converter : public TConverter {                             \
   public:                                                                    \
      T##name##Converter( Py_ssize_t size = -1 ) { fSize = size; }            \
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );\
      virtual PyObject* FromMemory( void* );                                  \
      virtual Bool_t ToMemory( PyObject*, void* );                            \
   private:                                                                   \
      Py_ssize_t fSize;                                                       \
   };                                                                         \
                                                                              \
   class T##name##RefConverter : public T##name##Converter {                  \
   public:                                                                    \
      using T##name##Converter::T##name##Converter;                           \
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );\
   }

// converters for built-ins
   PYROOT_DECLARE_BASIC_CONVERTER( Long );
   PYROOT_DECLARE_BASIC_CONVERTER( Bool );
   PYROOT_DECLARE_BASIC_CONVERTER( Char );
   PYROOT_DECLARE_BASIC_CONVERTER( UChar );
   PYROOT_DECLARE_BASIC_CONVERTER( Short );
   PYROOT_DECLARE_BASIC_CONVERTER( UShort );
   PYROOT_DECLARE_BASIC_CONVERTER( Int );
   PYROOT_DECLARE_BASIC_CONVERTER( ULong );
   PYROOT_DECLARE_BASIC_CONVERTER2( UInt, ULong );
   PYROOT_DECLARE_BASIC_CONVERTER( LongLong );
   PYROOT_DECLARE_BASIC_CONVERTER( ULongLong );
   PYROOT_DECLARE_BASIC_CONVERTER( Double );
   PYROOT_DECLARE_BASIC_CONVERTER( Float );
   PYROOT_DECLARE_BASIC_CONVERTER( LongDouble );

   PYROOT_DECLARE_REF_CONVERTER( Int );
   PYROOT_DECLARE_REF_CONVERTER( Long );
   PYROOT_DECLARE_REF_CONVERTER( Double );

   class TVoidConverter : public TConverter {
   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
   };

   class TCStringConverter : public TConverter {
   public:
      TCStringConverter( UInt_t maxSize = UINT_MAX ) : fMaxSize( maxSize ) {}

   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
      virtual PyObject* FromMemory( void* address );
      virtual Bool_t ToMemory( PyObject* value, void* address );

   protected:
      std::string fBuffer;
      UInt_t fMaxSize;
   };

   class TNonConstCStringConverter : public TCStringConverter {
   public:
      TNonConstCStringConverter( UInt_t maxSize = UINT_MAX ) : TCStringConverter( maxSize ) {}

   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
      virtual PyObject* FromMemory( void* address );
   };

   class TNonConstUCStringConverter : public TNonConstCStringConverter {
   public:
      TNonConstUCStringConverter( UInt_t maxSize = UINT_MAX ) : TNonConstCStringConverter( maxSize ) {}

   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
   };

// pointer/array conversions
   class TVoidArrayConverter : public TConverter {
   public:
      TVoidArrayConverter( Bool_t keepControl = kTRUE ) { fKeepControl = keepControl; }
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
      virtual PyObject* FromMemory( void* address );
      virtual Bool_t ToMemory( PyObject* value, void* address );

   protected:
      virtual Bool_t GetAddressSpecialCase( PyObject* pyobject, void*& address );

   protected:
      Bool_t KeepControl() { return fKeepControl; }

   private:
      Bool_t fKeepControl;
   };

   PYROOT_DECLARE_ARRAY_CONVERTER( BoolArray );
   PYROOT_DECLARE_ARRAY_CONVERTER( ShortArray );
   PYROOT_DECLARE_ARRAY_CONVERTER( UShortArray );
   PYROOT_DECLARE_ARRAY_CONVERTER( IntArray );
   PYROOT_DECLARE_ARRAY_CONVERTER( UIntArray );
   PYROOT_DECLARE_ARRAY_CONVERTER( LongArray );
   PYROOT_DECLARE_ARRAY_CONVERTER( ULongArray );
   PYROOT_DECLARE_ARRAY_CONVERTER( FloatArray );
   PYROOT_DECLARE_ARRAY_CONVERTER( DoubleArray );

   class TLongLongArrayConverter : public TVoidArrayConverter {
   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
   };

// converters for special cases
   class TCppObjectConverter : public TVoidArrayConverter {
   public:
      TCppObjectConverter( Cppyy::TCppType_t klass, Bool_t keepControl = kFALSE ) :
         TVoidArrayConverter( keepControl ), fClass( klass ), fObjProxy(nullptr) {}

   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
      virtual PyObject* FromMemory( void* address );
      virtual Bool_t ToMemory( PyObject* value, void* address );

   protected:
      Cppyy::TCppType_t fClass;
      ObjectProxy* fObjProxy;
   };

   class TStrictCppObjectConverter : public TCppObjectConverter {
   public:
      using TCppObjectConverter::TCppObjectConverter;

   protected:
      virtual Bool_t GetAddressSpecialCase( PyObject*, void*& ) { return kFALSE; }
   };

   class TValueCppObjectConverter : public TStrictCppObjectConverter {
   public:
      using TStrictCppObjectConverter::TStrictCppObjectConverter;

   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
   };

   class TRefCppObjectConverter : public TConverter  {
   public:
      TRefCppObjectConverter( Cppyy::TCppType_t klass ) : fClass( klass ), fObjProxy(nullptr) {}

   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );

   protected:
      Cppyy::TCppType_t fClass;
      ObjectProxy* fObjProxy;
   };

   template <bool ISREFERENCE>
   class TCppObjectPtrConverter : public TCppObjectConverter {
   public:
      using TCppObjectConverter::TCppObjectConverter;

   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
      virtual PyObject* FromMemory( void* address );
      virtual Bool_t ToMemory( PyObject* value, void* address );
   };

   extern template class TCppObjectPtrConverter<true>;
   extern template class TCppObjectPtrConverter<false>;

   class TCppObjectArrayConverter : public TCppObjectConverter {
   public:
      TCppObjectArrayConverter( Cppyy::TCppType_t klass, size_t size, Bool_t keepControl = kFALSE ) :
         TCppObjectConverter( klass, keepControl ), m_size( size ) {}

   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
      virtual PyObject* FromMemory( void* address );
      virtual Bool_t ToMemory( PyObject* value, void* address );

   protected:
      size_t m_size;
   };

// CLING WORKAROUND -- classes for STL iterators are completely undefined in that
// they come in a bazillion different guises, so just do whatever
   class TSTLIteratorConverter : public TConverter {
   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
   };
// -- END CLING WORKAROUND

   class TVoidPtrRefConverter : public TConverter {
   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
   };

   class TVoidPtrPtrConverter : public TConverter {
   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
      virtual PyObject* FromMemory( void* address );
   };

   PYROOT_DECLARE_BASIC_CONVERTER( PyObject );

#define PYROOT_DECLARE_STRING_CONVERTER( name, strtype )                      \
   class T##name##Converter : public TCppObjectConverter {                    \
   public:                                                                    \
      T##name##Converter( Bool_t keepControl = kTRUE );                       \
   public:                                                                    \
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );\
      virtual PyObject* FromMemory( void* address );                          \
      virtual Bool_t ToMemory( PyObject* value, void* address );              \
   private:                                                                   \
      strtype fBuffer;                                                        \
   }

   PYROOT_DECLARE_STRING_CONVERTER( TString,   TString );
   PYROOT_DECLARE_STRING_CONVERTER( STLString, std::string );
   PYROOT_DECLARE_STRING_CONVERTER( STLStringView, std::string_view );

   class TNotImplementedConverter : public TConverter {
   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* = 0 );
   };

// smart pointer converter
   class TSmartPtrCppObjectConverter : public TConverter  {
   public:
      TSmartPtrCppObjectConverter( Cppyy::TCppType_t klass,
                                   Cppyy::TCppType_t rawPtrType,
                                   Cppyy::TCppMethod_t deref,
                                   Bool_t keepControl = kFALSE,
                                   Bool_t handlePtr = kFALSE )
         : fClass( klass ), fRawPtrType( rawPtrType ), fDereferencer( deref ),
           fKeepControl( keepControl ), fHandlePtr( handlePtr ) {}

   public:
      virtual Bool_t SetArg( PyObject*, TParameter&, TCallContext* ctxt = 0 );
      virtual PyObject* FromMemory( void* address );
      //virtual Bool_t ToMemory( PyObject* value, void* address );

   protected:
      virtual Bool_t GetAddressSpecialCase( PyObject*, void*& ) { return kFALSE; }

      Cppyy::TCppType_t   fClass;
      Cppyy::TCppType_t   fRawPtrType;
      Cppyy::TCppMethod_t fDereferencer;
      Bool_t              fKeepControl;
      Bool_t              fHandlePtr;
   };

// create converter from fully qualified type
   TConverter* CreateConverter( const std::string& fullType, Long_t size = -1 );

} // namespace PyROOT

#endif // !PYROOT_CONVERTERS_H
