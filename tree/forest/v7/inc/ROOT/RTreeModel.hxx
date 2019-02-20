/// \file ROOT/RTreeModel.hxx
/// \ingroup Forest ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-10-04
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RTreeModel
#define ROOT7_RTreeModel

#include <ROOT/RStringView.hxx>
#include <ROOT/RTreeEntry.hxx>
#include <ROOT/RTreeField.hxx>
#include <ROOT/RTreeValue.hxx>

#include <TError.h>

#include <memory>
#include <utility>

namespace ROOT {
namespace Experimental {

// clang-format off
/**
\class ROOT::Experimental::RTreeModel
\ingroup Forest
\brief The RTreeModel encapulates the schema of a tree.

The tree model comprises a collection of hierarchically organized fields. From a frozen model, "entries"
can be extracted. For convenience, the model provides a default entry. Models have a unique model identifier
that faciliates checking whether entries are compatible with it (i.e.: have been extracted from that model).
A model needs to be frozen before it can be used to create an RTree.
*/
// clang-format on
class RTreeModel {
   /// Hierarchy of fields consiting of simple types and collections (sub trees)
   RTreeFieldRoot fRootField;
   /// Contains tree values corresponding to the created fields
   RTreeEntry fDefaultEntry;

public:
   /// Adds a field whose type is not known at compile time.  Thus there is no shared pointer returned.
   void AddField(std::unique_ptr<Detail::RTreeFieldBase> field);

   /// Creates a new field and a corresponding tree value that is managed by a shared pointer.
   template <typename T, typename... ArgsT>
   std::shared_ptr<T> AddField(std::string_view fieldName, ArgsT&&... args) {
      auto field = std::make_unique<RTreeField<T>>(fieldName);
      auto ptr = fDefaultEntry.AddValue<T>(field.get(), std::forward<ArgsT>(args)...);
      fRootField.Attach(std::move(field));
      return ptr;
   }

   RTreeFieldRoot* GetRootField() { return &fRootField; }
   RTreeEntry* GetDefaultEntry() { return &fDefaultEntry; }
};

} // namespace Exerimental
} // namespace ROOT

#endif
