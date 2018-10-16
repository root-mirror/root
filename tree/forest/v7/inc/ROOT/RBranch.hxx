/// \file ROOT/RBranch.hxx
/// \ingroup Forest ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-10-09
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RBranch
#define ROOT7_RBranch

#include <ROOT/RCargo.hxx>
#include <ROOT/RColumn.hxx>
#include <ROOT/RStringView.hxx>
#include <ROOT/RTreeUtil.hxx>

#include <memory>
#include <vector>

namespace ROOT {
namespace Experimental {

namespace Detail {

class RCargoBase;
class RPageStorage;

// clang-format off
/**
\class ROOT::Experimental::RBranchBase
\ingroup Forest
\brief An RBranchBase translates read and write calls from/to underlying columns to/from cargo objects

The RBranchBase and its type-safe descendants provide the object to column mapper. They map C++ objects
to primitive columns, where the mapping is trivial for simple types such as 'double'. The branch knows
based on its type and the branch name the type(s) and name(s) of the columns.
*/
// clang-format on
class RBranchBase {
private:
   /// A branch on a trivial type that maps as-is to a single column
   bool fIsSimple;
   /// All branches have a main column. For nested branches, the main column is the index branch. Points into fColumns.
   RColumn *fPrincipalColumn;
   /// The columns are connected either to a sink or to a source (not to both); they are owned by the branch.
   std::vector<RColumn> fColumns;

protected:
   /// Operations on values of complex types, e.g. ones that involve multiple columns or for which no direct
   /// column type exists.
   virtual void DoAppend(const RCargoBase &cargo) = 0;
   virtual void DoRead(TreeIndex_t index, const RCargoBase &cargo) = 0;
   virtual void DoReadV(TreeIndex_t index, TreeIndex_t count, void *dst) = 0;

public:
   /// The constructor creates the underlying column objects and connects them to either a sink or a source.
   RBranchBase(std::string_view name);
   virtual ~RBranchBase();

   /// Registeres the backing columns with the physical storage
   virtual void GenerateColumns(Detail::RPageStorage &storage) = 0;

   /// Generates a cargo object of the branch type.
   virtual std::unique_ptr<RCargoBase> GenerateCargo() = 0;

   /// Write the value stored in cargo to a tree. The cargo object has to be of the same type as the branch.
   void Append(const RCargoBase &cargo) {
     if (!fIsSimple) {
        DoAppend(cargo);
        return;
     }
     fPrincipalColumn->Append(*(cargo.fPrincipalElement));
   }

   /// Populate a single cargo object with data from the tree, which needs to be of the fitting type.
   /// Reading copies data into the memory given by cargo.
   void Read(TreeIndex_t index, RCargoBase &cargo) {
      if (!fIsSimple) {
         DoRead(index, cargo);
         return;
      }
      fPrincipalColumn->Read(index, cargo.fPrincipalElement);
   }

   /// Type unsafe bulk read interface; dst must point to a vector of objects of the branch type.
   /// TODO(jblomer): can this be type safe?
   void ReadV(TreeIndex_t index, TreeIndex_t count, void *dst)
   {
      if (!fIsSimple) {
         DoReadV(index, count, dst);
         return;
      }
      fPrincipalColumn->ReadV(index, count, dst);
   }

   /// Only for simple types, let the memory enclosed in cargo simply point into the page buffer.
   /// The resulting cargo object may only be used for as long as no request to another item of this branch is made
   /// because another index might trigger a swap of the page buffer.
   /// The dst location must be an object of the branch type.
   void Map(TreeIndex_t index, void **dst) {
      if (!fIsSimple) {
         // TODO(jblomer)
      }
      fPrincipalColumn->Map(index, dst);
   }

   /// The number of elements in the principal column. For top level branches, the number of entries.
   TreeIndex_t GetNItems();

   /// Ensure that all received items are written from page buffers to the storage.
   void Flush();

   RPageSource* GetSource();
};

} // namespace Detail

/// A Branch representing a collection
class RBranchSubtree : public Detail::RBranchBase {
protected:
   void DoAppend(const Detail::RCargoBase &cargo) final;
   void DoRead(TreeIndex_t index, const Detail::RCargoBase &cargo) final;
   void DoReadV(TreeIndex_t index, TreeIndex_t count, void *dst) final;

public:
   RBranchSubtree(std::string_view name);
   ~RBranchSubtree();

   void GenerateColumns(Detail::RPageStorage &storage) final;
   std::unique_ptr<Detail::RCargoBase> GenerateCargo() final;
   void Attach(RBranchBase* child);
};


/// Supported types are implemented as template specializations
template <typename T>
class RBranch : public Detail::RBranchBase {
};

/// TODO(jblomer): template specializations for simple types and TClass

} // namespace Experimental
} // namespace ROOT

#endif
