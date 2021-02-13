
// Author: Enric Tejedor CERN  08/2019

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef PYROOT_TMEMORYREGULATOR_H
#define PYROOT_TMEMORYREGULATOR_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMemoryRegulator                                                     //
//                                                                      //
// Sets hooks in Cppyy's MemoryRegulator to keep track of the TObjects  //
// that are constructed and destructed. For those objects, a map is     //
// filled, where the key is the address of the object and the value is  //
// the class to which the object belongs.                               //
//                                                                      //
// The TMemoryRegulator object, created in PyROOTWrapper.cxx, is added  //
// to the list of cleanups and its RecursiveRemove method is called by  //
// ROOT to manage the memory of TObjects being deleted.                 //
// In RecursiveRemove, the object being deleted is already a TNamed, so //
// the information about its actual class is not available anymore.     //
// To solve the problem, the map above is used to know the class of the //
// object, so that Cppyy's RecursiveRemove can be called passing the    //
// class as argument.                                                   //
//////////////////////////////////////////////////////////////////////////

// Bindings
// CPyCppyy.h must be go first, since it includes Python.h, which must be
// included before any standard header
#include "CPyCppyy.h"
#include "MemoryRegulator.h"

// ROOT
#include "TObject.h"
#include "TClass.h"

// Stl
#include <unordered_map>

namespace PyROOT {

typedef std::unordered_map<Cppyy::TCppObject_t, Cppyy::TCppType_t> ObjectMap_t;

class TMemoryRegulator : public TObject {
private:
   static ObjectMap_t fObjectMap; // key: object address; value: object class id

   static std::pair<bool, bool> RegisterHook(Cppyy::TCppObject_t, Cppyy::TCppType_t);

   static std::pair<bool, bool> UnregisterHook(Cppyy::TCppObject_t, Cppyy::TCppType_t);

public:
   TMemoryRegulator();

   virtual void RecursiveRemove(TObject *);

   void ClearProxiedObjects();
};

} // namespace PyROOT

#endif // !PYROOT_TMEMORYREGULATOR_H
