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

#include <iostream>

namespace {

// legacy bridge interface type
using RNTupleList = std::vector<std::pair<
   std::string,  // path
   std::string   // ntuple name
>>;

} // anonymous namespace

// legacy bridge interface to hadd via TFileMerger
Long64_t ROOT::Experimental::RNTuple::Merge(TCollection* inputs, TFileMergeInfo* mergeInfo) {
   if (inputs == nullptr || mergeInfo == nullptr) {
      return -1;
   }
   std::string& output_file = *(static_cast<std::string*>(static_cast<void*>(mergeInfo)));
   if (output_file.empty()) {
      return -1;
   }
   std::cout << "RNTuple merger output file is " << output_file << "\n";

   RNTupleList& ntuples = *(static_cast<RNTupleList*>(static_cast<void*>(inputs)));
   if (ntuples.size() == 0) {
      std::cout << "got empty list of RNTuples to merge\n";
      return -1;
   }
   for (const auto& ntpl: ntuples) {
      std::cout << "merging ntuple from file '" << ntpl.first << " named '" << ntpl.second << "'\n";
   }
   // todo(max)
   // 1. open the RNTuple page sources with the ntuple names and paths
   // std::vector<std::unique_ptr<RPageSource>> merge_inputs;
   // for (const auto& ntpl: ntuples) {
   //    merge_inputs.emplace_back(RPageSource::Create(ntpl.second, ntpl.first));
   // }
   //
   // 2. open the output file as a page sink
   // - use the first ntuple name as the output ntuple name (?)
   // - can we assume they're all the same?
   // std::unique_ptr<RPageSink> merge_output = RPageSink::Create(ntuples.front().second, output_file);
   //
   // 3. call the merge function on the page sources and write it to the page sink
   // auto res = MergeRNTuples(merge_inputs, merge_output);
   // if (!res) {
   //    return -1;
   // }
   std::cout << "RNTuple merging is unimplemented\n";
   return -1;
}


////////////////////////////////////////////////////////////////////////////////


ROOT::Experimental::RResult<ROOT::Experimental::RFieldMerger>
ROOT::Experimental::RFieldMerger::Merge(const ROOT::Experimental::RFieldDescriptor &lhs,
   const ROOT::Experimental::RFieldDescriptor &rhs)
{
   return R__FAIL("couldn't merge field " + lhs.GetFieldName() + " with field "
      + rhs.GetFieldName() + " (unimplemented!)");
}
