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

Long64_t ROOT::Experimental::RNTuple::Merge(TCollection* inputs, TFileMergeInfo* mergeInfo) {
   if (inputs == nullptr || mergeInfo == nullptr) {
      return -1;
   }
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
