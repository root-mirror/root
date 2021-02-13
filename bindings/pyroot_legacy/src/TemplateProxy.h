// Author: Wim Lavrijsen, Jan 2008

#ifndef PYROOT_TEMPLATEPROXY_H
#define PYROOT_TEMPLATEPROXY_H

// ROOT
#include "DllImport.h"

// Standard
#include <string>


namespace PyROOT {

   class PyCallable;
   class MethodProxy;

/** Template proxy object to return functions and methods
      @author  WLAV
      @date    01/15/2008
      @version 1.0
 */

   class TemplateProxy {
   public:
      void Set( const std::string& name, PyObject* pyclass );

   public:               // public, as the python C-API works with C structs
      PyObject_HEAD
      PyObject* fSelf;             // must be first (same layout as MethodProxy)
      PyObject* fPyClass;
      PyObject* fPyName;
      MethodProxy* fNonTemplated;  // holder for non-template overloads
      MethodProxy* fTemplated;     // holder for templated overloads

   public:
      void AddOverload( MethodProxy* mp );
      void AddOverload( PyCallable* pc );
      void AddTemplate( PyCallable* pc );

   private:              // private, as the python C-API will handle creation
      TemplateProxy() {}
   };


//- template proxy type and type verification --------------------------------
   R__EXTERN PyTypeObject TemplateProxy_Type;

   template< typename T >
   inline Bool_t TemplateProxy_Check( T* object )
   {
      return object && PyObject_TypeCheck( object, &TemplateProxy_Type );
   }

   template< typename T >
   inline Bool_t TemplateProxy_CheckExact( T* object )
   {
      return object && Py_TYPE(object) == &TemplateProxy_Type;
   }

//- creation -----------------------------------------------------------------
   inline TemplateProxy* TemplateProxy_New( const std::string& name, PyObject* pyclass )
   {
   // Create and initialize a new template method proxy for the class.
      TemplateProxy* pytmpl = (TemplateProxy*)TemplateProxy_Type.tp_new( &TemplateProxy_Type, 0, 0 );
      pytmpl->Set( name, pyclass );
      return pytmpl;
   }

} // namespace PyROOT

#endif // !PYROOT_TEMPLATEPROXY_H
