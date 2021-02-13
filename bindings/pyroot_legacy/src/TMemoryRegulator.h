// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Apr 2004

#ifndef PYROOT_TMEMORYREGULATOR_H
#define PYROOT_TMEMORYREGULATOR_H

// ROOT
#include "TObject.h"

// Standard
#include <unordered_map>


namespace PyROOT {

/** Communicate object destruction across ROOT/CINT/PyROOT/
      @author  WLAV
      @date    11/23/2004
      @version 2.2
 */

   class ObjectProxy;

   class TMemoryRegulator : public TObject {
   public:
      TMemoryRegulator();
      ~TMemoryRegulator();

   // callback for ROOT/CINT
      virtual void RecursiveRemove( TObject* object );

   // cleanup of all tracked objects
      void ClearProxiedObjects();

   // add a python object to the table of managed objects
      static Bool_t RegisterObject( ObjectProxy* pyobj, TObject* object );

   // remove a python object from the table of managed objects, w/o notification
      static Bool_t UnregisterObject( TObject* object );

   // new reference to python object corresponding to object, or 0 on failure
      static PyObject* RetrieveObject( TObject* object, Cppyy::TCppType_t klass );

   // callback when weak refs to managed objects are destroyed
      static PyObject* ObjectEraseCallback( PyObject*, PyObject* pyref );

   private:
      typedef std::unordered_map< TObject*, PyObject* > ObjectMap_t;
      typedef std::unordered_map< PyObject*, ObjectMap_t::iterator > WeakRefMap_t;

      static ObjectMap_t*  fgObjectTable;
      static WeakRefMap_t* fgWeakRefTable;
   };

} // namespace PyROOT

#endif // !PYROOT_TMEMORYREGULATOR_H
