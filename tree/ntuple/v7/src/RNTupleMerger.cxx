/// \file RNTupleMerger.cxx
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

#include <ROOT/RError.hxx>
#include <ROOT/RMiniFile.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleMerger.hxx>
#include <ROOT/RNTupleUtil.hxx>

#include <TError.h>

Long64_t ROOT::Experimental::RNTuple::Merge(TCollection* inputs, TFileMergeInfo* mergeInfo) {
   if (inputs == nullptr || mergeInfo == nullptr) {
      return -1;
   }
   return -1;
}


////////////////////////////////////////////////////////////////////////////////


void ROOT::Experimental::RFieldMerger::MakeReference(const RFieldHandle &referenceHandle)
{
   fReferenceFields[referenceHandle.fid] = referenceHandle.desc.GetFieldDescriptor(referenceHandle.fid).Clone();
   fFieldIdMaps[0].Insert(referenceHandle.fid, referenceHandle.fid);
   fSubFieldLinks[referenceHandle.fid] = std::vector<DescriptorId_t>();
   fColumnLinks[referenceHandle.fid] = std::vector<DescriptorId_t>();

   // Register columns connected to the field
   for (const auto &c : referenceHandle.desc.GetColumnRange(referenceHandle.fid)) {
      fColumnIdMaps[0].Insert(c.GetId(), c.GetId());
      fColumnLinks[referenceHandle.fid].emplace_back(c.GetId());
   }

   // Recurse into the sub fields
   for (const auto &f : referenceHandle.desc.GetFieldRange(referenceHandle.fid)) {
      MakeReference({referenceHandle.desc, f.GetId()});
      fSubFieldLinks[referenceHandle.fid].emplace_back(f.GetId());
   }
}


ROOT::Experimental::RFieldMerger::RFieldMerger(const RFieldHandle &referenceHandle)
   : fReferenceId(referenceHandle.fid)
{
   fFieldIdMaps.emplace_back(RBiMap());
   fColumnIdMaps.emplace_back(RBiMap());
   MakeReference(referenceHandle);
}


ROOT::Experimental::RResult<void>
ROOT::Experimental::RFieldMerger::MergeImpl(DescriptorId_t referenceId, const RFieldHandle &inputHandle)
{
   const RFieldDescriptor &lhs = fReferenceFields.at(referenceId);
   const RFieldDescriptor &rhs = inputHandle.desc.GetFieldDescriptor(inputHandle.fid);
   auto qualifiedName = inputHandle.desc.GetQualifiedFieldName(inputHandle.fid);

   // TODO(jblomer): compare _normalized_ type name
   if (lhs.GetTypeName() != rhs.GetTypeName())
      return R__FAIL("field merge error: type mismatch for " + qualifiedName);
   if (lhs.GetFieldVersion() != rhs.GetFieldVersion())
      return R__FAIL("field merge error: field version mismatch for " + qualifiedName);
   if (lhs.GetTypeVersion() != rhs.GetTypeVersion())
      return R__FAIL("field merge error: type version mismatch for " + qualifiedName);
   if (lhs.GetNRepetitions() != rhs.GetNRepetitions())
      return R__FAIL("field merge error: fixed-size array mismatch for " + qualifiedName);
   R__ASSERT(lhs.GetStructure() == rhs.GetStructure());

   // Map the field and column IDs between input and reference
   fFieldIdMaps.rbegin()->Insert(inputHandle.fid, referenceId);
   const auto &referenceColumnIds = fColumnLinks.at(referenceId);
   for (size_t i = 0; i < referenceColumnIds.size(); ++i) {
      auto inputColumnId = inputHandle.desc.FindColumnId(inputHandle.fid, i);
      R__ASSERT(inputColumnId != kInvalidDescriptorId);
      fColumnIdMaps.rbegin()->Insert(inputColumnId, referenceColumnIds[i]);
   }

   // Recurse into the sub fields in the reference field tree
   for (auto subId : fSubFieldLinks.at(referenceId)) {
      const auto &subField = fReferenceFields.at(subId);
      auto matchingInput = inputHandle.desc.FindFieldId(subField.GetFieldName(), inputHandle.fid);
      if (matchingInput == kInvalidDescriptorId)
         return R__FAIL("field merge error: missing field " + qualifiedName + "." + subField.GetFieldName());

      auto success = MergeImpl(subId, {inputHandle.desc, matchingInput});
      if (!success)
         return R__FORWARD_ERROR(success);
   }

   return RResult<void>::Success();
}


ROOT::Experimental::RResult<int> ROOT::Experimental::RFieldMerger::Merge(const RFieldHandle &inputHandle)
{
   int inputId = fFieldIdMaps.size();
   fFieldIdMaps.emplace_back(RBiMap());
   fColumnIdMaps.emplace_back(RBiMap());

   auto result = MergeImpl(fReferenceId, inputHandle);
   if (!result) {
      fFieldIdMaps.pop_back();
      fColumnIdMaps.pop_back();
      return R__FORWARD_ERROR(result);
   }

   return inputId;
}
