
// Author: Enric Tejedor CERN  08/2019

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TMemoryRegulator.h"

#include "ProxyWrappers.h"
#include "CPPInstance.h"
#include "CPPInstance.h"

using namespace CPyCppyy;

PyROOT::ObjectMap_t PyROOT::TMemoryRegulator::fObjectMap = PyROOT::ObjectMap_t();

////////////////////////////////////////////////////////////////////////////
/// \brief Constructor. Registers the hooks to run on Cppyy's object
///        construction and destruction
PyROOT::TMemoryRegulator::TMemoryRegulator()
{
   MemoryRegulator::SetRegisterHook(PyROOT::TMemoryRegulator::RegisterHook);
   MemoryRegulator::SetUnregisterHook(PyROOT::TMemoryRegulator::UnregisterHook);
}

////////////////////////////////////////////////////////////////////////////
/// \brief Register a hook that Cppyy runs when constructing an object.
/// \param[in] cppobj Address of the object.
/// \param[in] klass Class id of the object.
/// \return Pair of two booleans. First indicates success, second tells
///         Cppyy if we want to continue running RegisterPyObject
std::pair<bool, bool> PyROOT::TMemoryRegulator::RegisterHook(Cppyy::TCppObject_t cppobj, Cppyy::TCppType_t klass)
{
   static Cppyy::TCppType_t tobjectTypeID = (Cppyy::TCppType_t)Cppyy::GetScope("TObject");

   if (Cppyy::IsSubtype(klass, tobjectTypeID)) {
      ObjectMap_t::iterator ppo = fObjectMap.find(cppobj);
      if (ppo == fObjectMap.end()) {
         fObjectMap.insert({cppobj, klass});
      }
   }

   return {true, true};
}

////////////////////////////////////////////////////////////////////////////
/// \brief Register a hook that Cppyy runs when deleting an object.
/// \param[in] cppobj Address of the object.
/// \param[in] klass Class id of the object.
/// \return Pair of two booleans. First indicates success, second tells
///         Cppyy if we want to continue running UnRegisterPyObject
std::pair<bool, bool> PyROOT::TMemoryRegulator::UnregisterHook(Cppyy::TCppObject_t cppobj, Cppyy::TCppType_t klass)
{
   static Cppyy::TCppType_t tobjectTypeID = (Cppyy::TCppType_t)Cppyy::GetScope("TObject");

   if (Cppyy::IsSubtype(klass, tobjectTypeID)) {
      ObjectMap_t::iterator ppo = fObjectMap.find(cppobj);
      if (ppo != fObjectMap.end()) {
         fObjectMap.erase(ppo);
      }
   }

   return {true, true};
}

////////////////////////////////////////////////////////////////////////////
/// \brief Get the class id of the TObject being deleted and run Cppyy's
///        RecursiveRemove.
/// \param[in] object Object being destructed.
void PyROOT::TMemoryRegulator::RecursiveRemove(TObject *object)
{
   auto cppobj = (Cppyy::TCppObject_t)object;
   Cppyy::TCppType_t klass = 0;

   ObjectMap_t::iterator ppo = fObjectMap.find(cppobj);
   if (ppo != fObjectMap.end()) {
      klass = ppo->second;
      MemoryRegulator::RecursiveRemove(cppobj, klass);
      fObjectMap.erase(ppo);
   }
}

////////////////////////////////////////////////////////////////////////////
/// \brief Clean up all tracked objects.
void PyROOT::TMemoryRegulator::ClearProxiedObjects()
{
   while (!fObjectMap.empty()) {
      auto elem = fObjectMap.begin();
      auto cppobj = elem->first;
      auto klassid = elem->second;
      auto pyclass = CreateScopeProxy(klassid);
      auto pyobj = (CPPInstance *)MemoryRegulator::RetrievePyObject(cppobj, pyclass);

      if (pyobj && (pyobj->fFlags & CPPInstance::kIsOwner)) {
         // Only delete the C++ object if the Python proxy owns it.
         // If it is a value, cppyy deletes it in RecursiveRemove as part of
         // the proxy cleanup.
         auto o = static_cast<TObject *>(cppobj);
         bool isValue = pyobj->fFlags & CPPInstance::kIsValue;
         RecursiveRemove(o);
         if (!isValue)
            delete o;
      }
      else {
         // Non-owning proxy, just unregister to clean tables.
         // The proxy deletion by Python will have no effect on C++, so all good
         MemoryRegulator::UnregisterPyObject(pyobj, pyclass);
      }
   }
}
