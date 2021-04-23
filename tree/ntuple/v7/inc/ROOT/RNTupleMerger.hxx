/// \file ROOT/RNTupleMerger.hxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch> & Max Orok <maxwellorok@gmail.com>
/// \date 2020-07-08
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RNTupleMerger
#define ROOT7_RNTupleMerger

#include <ROOT/RError.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleUtil.hxx>

#include <map>
#include <memory>
#include <vector>

namespace ROOT {
namespace Experimental {

// clang-format off
/**
\class ROOT::Experimental::RFieldMerger
\ingroup NTuple
\brief Two-way merge between NTuple fields

Provides a mapping of field IDs and column IDs between structurally equivalent field trees.
The merger is initialized with a reference field and its sub fields. It can then merge (map) multiple
input fields. Successful merges return an input ID that can be subsequently used to map IDs of the
input field to the reference field and vice versa. The reference field has input ID 0.

In order to merge successfully, the input field must have at least the same type and subfields than
the reference field. The input field can have additional sub fields, which are ignored.

Normally, the constuctor and the Merge() method are called with zero fields, which would merge
complete NTuple schemas. But it can be called with fields deeper in a hierarchy, too.
*/
// clang-format on
class RFieldMerger {
public:
   /// In order to inpect a field, its subfields and the corresponding columns, both a field id and an RNTuple
   /// descriptor are required
   struct RFieldHandle {
      const RNTupleDescriptor &desc;
      DescriptorId_t fid = kInvalidDescriptorId;
   };

private:
   /// Helper class to map descriptor IDs in both directions, from the reference to an input and vice versa.
   struct RBiMap {
      std::unordered_map<int, int> fInput2Reference;
      std::unordered_map<int, int> fReference2Input;

      void Insert(int input, int reference)
      {
         fInput2Reference[input] = reference;
         fReference2Input[reference] = input;
      }
   };

   /// The field ID of the top-most reference field
   DescriptorId_t fReferenceId;

   /// A clone of the initial field passed in the constructor and all its sub fields.
   std::unordered_map<DescriptorId_t, RFieldDescriptor> fReferenceFields;
   /// Maps from parent field ID to sub field IDs in the reference field tree
   std::unordered_map<DescriptorId_t, std::vector<DescriptorId_t>> fSubFieldLinks;
   /// Maps field ID to its column IDs in the reference field tree. The order of the columns is preserved,
   /// so that the vector index is identical to the column index (RColumnDescriptor::fIndex).
   std::unordered_map<DescriptorId_t, std::vector<DescriptorId_t>> fColumnLinks;

   /// Bi-directional mapping between input field IDs and reference field IDs.
   /// The input ID is used as an index.
   std::vector<RBiMap> fFieldIdMaps;
   /// Bi-directional mapping between input column IDs and reference column IDs.
   /// The input ID is used as an index.
   std::vector<RBiMap> fColumnIdMaps;

   void MakeReference(const RFieldHandle &referenceHandle);

   RResult<void> MergeImpl(DescriptorId_t referenceId, const RFieldHandle &inputHandle);

public:
   explicit RFieldMerger(const RFieldHandle &referenceHandle);

   /// Returns a unique input ID on success, i.e. if the input field in the input handle is structurally
   /// equivalent to the reference field.
   RResult<int> Merge(const RFieldHandle &inputHandle);

   DescriptorId_t GetReferenceFieldId(DescriptorId_t fieldId, int inputId) const {
      return fFieldIdMaps[inputId].fInput2Reference.at(fieldId);
   }
   DescriptorId_t GetReferenceColumnId(DescriptorId_t columnId, int inputId) const {
      return fColumnIdMaps[inputId].fInput2Reference.at(columnId);
   }

   DescriptorId_t GetInputFieldId(DescriptorId_t fieldId, int inputId) const {
      return fFieldIdMaps[inputId].fReference2Input.at(fieldId);
   }
   DescriptorId_t GetInputColumnId(DescriptorId_t columnId, int inputId) const {
      return fColumnIdMaps[inputId].fReference2Input.at(columnId);
   }
};

} // namespace Experimental
} // namespace ROOT

#endif
