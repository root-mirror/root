// Author: Enrico Guiraud, Danilo Piparo, Massimo Tumolo CERN  06/2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RDFBOOKEDCUSTOMCOLUMNS
#define ROOT_RDFBOOKEDCUSTOMCOLUMNS

#include <memory>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include "TString.h"

namespace ROOT {
namespace Detail {
namespace RDF {
class RDefineBase;
}
}

namespace Internal {
namespace RDF {

namespace RDFDetail = ROOT::Detail::RDF;

/**
 * \class ROOT::Internal::RDF::RBookedDefines
 * \ingroup dataframe
 * \brief Encapsulates the columns defined by the user
 */

class RBookedDefines {
   using RDefineBasePtrMap_t = std::map<std::string, std::shared_ptr<RDFDetail::RDefineBase>>;
   using ColumnNames_t = std::vector<std::string>;

   // Since RBookedDefines is meant to be an immutable, copy-on-write object, the actual values are set as const
   using RDefineBasePtrMapPtr_t = std::shared_ptr<const RDefineBasePtrMap_t>;
   using ColumnNamesPtr_t = std::shared_ptr<const ColumnNames_t>;

private:
   RDefineBasePtrMapPtr_t fDefines;
   ColumnNamesPtr_t fDefinesNames;  // also abused to keep track of aliases for each branch of the computation graph

public:
   ////////////////////////////////////////////////////////////////////////////
   /// \brief Copy-ctor for RBookedDefines.
   RBookedDefines(const RBookedDefines &) = default;

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Move-ctor for RBookedDefines.
   RBookedDefines(RBookedDefines &&) = default;

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Copy-assignment operator for RBookedDefines.
   RBookedDefines &operator=(const RBookedDefines &) = default;

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Creates the object starting from the provided maps
   RBookedDefines(RDefineBasePtrMapPtr_t defines, ColumnNamesPtr_t defineNames)
      : fDefines(defines), fDefinesNames(defineNames)
   {
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Creates a new wrapper with empty maps
   RBookedDefines()
      : fDefines(std::make_shared<RDefineBasePtrMap_t>()),
        fDefinesNames(std::make_shared<ColumnNames_t>())
   {
   }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Returns the list of the names of the defined columns
   ColumnNames_t GetNames() const { return *fDefinesNames; }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Returns the list of the pointers to the defined columns
   const RDefineBasePtrMap_t &GetColumns() const { return *fDefines; }

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Check if the provided name is tracked in the names list
   bool HasName(std::string_view name) const;

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Add a new booked column.
   /// Internally it recreates the map with the new column, and swaps it with the old one.
   void AddColumn(const std::shared_ptr<RDFDetail::RDefineBase> &column, std::string_view name);

   ////////////////////////////////////////////////////////////////////////////
   /// \brief Add a new name to the list returned by `GetNames` without booking a new column.
   ///
   /// This is needed because we abuse fDefinesNames to also keep track of the aliases defined
   /// in each branch of the computation graph.
   /// Internally it recreates the vector with the new name, and swaps it with the old one.
   void AddName(std::string_view name);
};

} // Namespace RDF
} // Namespace Internal
} // Namespace ROOT

#endif
