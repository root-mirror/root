// Author: Enrico Guiraud, CERN  08/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RDF/RDefineReader.hxx>
#include <ROOT/RDF/RDefineBase.hxx>
#include <ROOT/RDF/Utils.hxx> // TypeID2TypeName
#include <TClass.h>

#include <stdexcept> // std::runtime_error
#include <string>
#include <typeinfo>

void ROOT::Internal::RDF::CheckDefineType(RDefineBase &define, const std::type_info &tid)
{
   const auto &colTId = define.GetTypeId();

   // Here we compare names and not typeinfos since they may come from two different contexts: a compiled
   // and a jitted one.
   const auto diffTypes = (0 != std::strcmp(colTId.name(), tid.name()));
   auto inheritedType = [&]() {
      auto colTClass = TClass::GetClass(colTId);
      return colTClass && colTClass->InheritsFrom(TClass::GetClass(tid));
   };

   if (diffTypes && !inheritedType()) {
      const auto tName = TypeID2TypeName(tid);
      const auto colTypeName = TypeID2TypeName(colTId);
      std::string errMsg = "RDefineReader: column \"" + define.GetName() + "\" is being used as ";
      if (tName.empty()) {
         errMsg += tid.name();
         errMsg += " (extracted from type info)";
      } else {
         errMsg += tName;
      }
      errMsg += " but defined column has type ";
      if (colTypeName.empty()) {
         auto &id = colTId;
         errMsg += id.name();
         errMsg += " (extracted from type info)";
      } else {
         errMsg += colTypeName;
      }
      throw std::runtime_error(errMsg);
   }
}
