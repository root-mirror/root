/// \file RPageStorageFriend.hxx
/// \ingroup NTuple ROOT7
/// \author Simon Leisibach <simon.leisibach@gmail.com>
/// \date 2019-10-04
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RPageStorageFriend
#define ROOT7_RPageStorageFriend

#include <ROOT/RNTupleOptions.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleMetrics.hxx>
#include <ROOT/RPageStorage.hxx>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace ROOT {
namespace Experimental {
namespace Detail {

// clang-format off
/**
\class ROOT::Experimental::Detail::RPageSourceFriend
\ingroup NTuple
\brief A pagesource generated from multiple files with the same clusters and number of entries. It acts like a PageSource for a file where multiple files were merged into one.
 
 An instance of RPageSourceFriend is created in RNTupleReader::Open() when a std::vector of filenames is passed as an argument instead of a filename-string and in RNTupleReader::ChainReader(). It first creates a RPageSource (including its descriptor) for each filename passed and initalizes its members. Later it merges the information from all the descriptors to create a single descriptor containing all information. After that its main job is to assign PopulatePage() and ReleasePage() to the correct RPageSources in fSources.
 
 The total number of columns of a reader with this RPageSource is the sum of columns for all the RPageSources in fSources.
*/
// clang-format on

// Note(lesimon): Currently is will fail if the size of the i-th cluster in a file is not the same for all files. Should
// it be changed to deal with clusters of varying sizes across files?
class RPageSourceFriend : public RPageSource {
private:
   /// Holds a RPageSource pointer for each file.
   std::vector<std::unique_ptr<RPageSource>> fSources;
   /// Holds the cumulative number of columns per file. It's size is number of files + 1.
   /// fNColumnPerSource[i] holds the number of columns of the i-th file, fNColumnPerSource[0] is always = 0
   std::vector<std::size_t> fNColumnPerSource;
   /// Keeps track to which RPageSource a populated page belongs to and how often the same page was populated but not
   /// released yet.
   struct PageInfoFriend {
      /// Tells that the RPage belongs to the RPageSource fSources.at(fSourceId).
      std::size_t fSourceId;
      /// Tells how often many instances of the same page were populated but not released yet.
      std::size_t fNSamePagePopulated;
   };
   /// Maps the buffer of a RPage (void*) to its RPageSource.
   std::unordered_map<void *, PageInfoFriend> fPageMapper;
   /// Is set to true when the meta-data of the clusters don't match.
   /// Getting pages from a unsafe RPageStorageFrame can lead to undefined behaviour.
   bool fUnsafe = false;
   // TODO(jblomer): implement fMetrics
   RNTupleMetrics fMetrics /*= RNTupleMetrics(std::string{"RPageStorageFriend"})*/;

protected:
   RNTupleDescriptor DoAttach();
   /// Checks if the field and column meta-data of all files are the same. If the meta-data doesn't match
   /// RNTupleReader::Open and RNTupleReader::ChainReader will return a nullptr.
   bool CompareFileMetaData();
   /// initializes fNColumnPerSource
   void InitializeMemberVariables();

public:
   RPageSourceFriend(std::string_view ntupleName, std::vector<std::string> locationVec,
                     const RNTupleReadOptions &options = RNTupleReadOptions());
   RPageSourceFriend(std::string_view ntupleName, std::vector<RPageSource *> sources,
                     const RNTupleReadOptions &options = RNTupleReadOptions());
   RPageSourceFriend(std::string_view ntupleName, std::vector<std::unique_ptr<RPageSource>> &&sources,
                     const RNTupleReadOptions &options = RNTupleReadOptions());
   ~RPageSourceFriend() = default;
   std::unique_ptr<RPageSource> Clone() const final;

   RPage PopulatePage(ColumnHandle_t columnHandle, NTupleSize_t globalIndex) final;
   RPage PopulatePage(ColumnHandle_t columnHandle, const RClusterIndex &clusterIndex) final;
   void ReleasePage(RPage &page) final;

   bool IsUnsafe() const { return fUnsafe; }
   RNTupleMetrics &GetMetrics() final { return fMetrics; }
};

} // namespace Detail
} // namespace Experimental
} // namespace ROOT

#endif
