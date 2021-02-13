// Author: Stefan Wunsch, Enrico Guiraud CERN  09/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RDF_RRESULTHANDLE
#define ROOT_RDF_RRESULTHANDLE

#include "ROOT/RResultPtr.hxx"
#include "ROOT/RDF/RLoopManager.hxx"
#include "ROOT/RDF/RActionBase.hxx"
#include "ROOT/RDF/Utils.hxx" // TypeID2TypeName

#include <memory>
#include <sstream>
#include <typeinfo>
#include <stdexcept> // std::runtime_error

namespace ROOT {
namespace RDF {

class RResultHandle {
   ROOT::Detail::RDF::RLoopManager* fLoopManager; //< Pointer to the loop manager
   /// Owning pointer to the action that will produce this result.
   /// Ownership is shared with RResultPtrs and RResultHandles that refer to the same result.
   std::shared_ptr<ROOT::Internal::RDF::RActionBase> fActionPtr;
   std::shared_ptr<void> fObjPtr; ///< Type erased shared pointer encapsulating the wrapped result
   const std::type_info& fType; ///< Type of the wrapped result

   // The ROOT::RDF::RunGraphs helper has to access the loop manager to check whether two RResultHandles belong to the same computation graph
   friend void RunGraphs(std::vector<RResultHandle>);

   /// Get the pointer to the encapsulated result.
   /// Ownership is not transferred to the caller.
   /// Triggers event loop and execution of all actions booked in the associated RLoopManager.
   void* Get()
   {
      if (!fActionPtr->HasRun())
         fLoopManager->Run();
      return fObjPtr.get();
   }

   /// Compare given type to the type of the wrapped result and throw if the types don't match.
   void CheckType(const std::type_info& type)
   {
      if (fType != type) {
         std::stringstream ss;
         ss << "Got the type " << ROOT::Internal::RDF::TypeID2TypeName(type)
            << " but the RResultHandle refers to a result of type " << ROOT::Internal::RDF::TypeID2TypeName(fType)
            << ".";
         throw std::runtime_error(ss.str());
      }
   }

public:
   template <class T>
   RResultHandle(const RResultPtr<T> &resultPtr) : fLoopManager(resultPtr.fLoopManager), fActionPtr(resultPtr.fActionPtr), fObjPtr(resultPtr.fObjPtr), fType(typeid(T)) {}

   RResultHandle(const RResultHandle&) = default;
   RResultHandle(RResultHandle&&) = default;

   /// Get the pointer to the encapsulated object.
   /// Triggers event loop and execution of all actions booked in the associated RLoopManager.
   /// \tparam T Type of the action result
   template <class T>
   T* GetPtr()
   {
      CheckType(typeid(T));
      return static_cast<T*>(Get());
   }

   /// Get a const reference to the encapsulated object.
   /// Triggers event loop and execution of all actions booked in the associated RLoopManager.
   /// \tparam T Type of the action result
   template <class T>
   const T& GetValue()
   {
      CheckType(typeid(T));
      return *static_cast<T*>(Get());
   }

   /// Check whether the result has already been computed
   ///
   /// ~~~{.cpp}
   /// std::vector<RResultHandle> results;
   /// results.emplace_back(df.Mean<double>("var"));
   /// res.IsReady(); // false, access will trigger event loop
   /// std::cout << res.GetValue<double>() << std::endl; // triggers event loop
   /// res.IsReady(); // true
   /// ~~~
   bool IsReady() const {
      return fActionPtr->HasRun();
   }

   bool operator==(const RResultHandle &rhs) const
   {
      return fObjPtr == rhs.fObjPtr;
   }

   bool operator!=(const RResultHandle &rhs) const
   {
      return !(fObjPtr == rhs.fObjPtr);
   }
};

} // namespace RDF
} // namespace ROOT

#endif // ROOT_RDF_RRESULTHANDLE
