/// \file ROOT/RNTupleOptions.hxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2019-08-25
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RNTupleOptions
#define ROOT7_RNTupleOptions

#include <Compression.h>
#include <ROOT/RNTupleUtil.hxx>

namespace ROOT {
namespace Experimental {

// clang-format off
/**
\class ROOT::Experimental::ENTupleContainerFormat
\ingroup NTuple
\brief Describes the options for wrapping RNTuple data in files
*/
// clang-format on
enum class ENTupleContainerFormat {
   kTFile, // ROOT TFile
   kBare, // A thin envelope supporting a single RNTuple only
};


// clang-format off
/**
\class ROOT::Experimental::RNTupleWriteOptions
\ingroup NTuple
\brief Common user-tunable settings for storing ntuples

All page sink classes need to support the common options.
*/
// clang-format on
class RNTupleWriteOptions {
   int fCompression{RCompressionSetting::EDefaults::kUseAnalysis};
   ENTupleContainerFormat fContainerFormat{ENTupleContainerFormat::kTFile};
   /// Whether the user has explicitly set these options, e.g. to override a
   /// TFile's compression setting.
   bool fIsCompressionOverride = false;
   NTupleSize_t fNEntriesPerCluster = 64000;
   NTupleSize_t fNElementsPerPage = 10000;
   bool fUseBufferedWrite = true;

public:
   virtual ~RNTupleWriteOptions() = default;
   virtual std::unique_ptr<RNTupleWriteOptions> Clone() const
   { return std::make_unique<RNTupleWriteOptions>(*this); }

   int GetCompression() const { return fCompression; }
   void SetCompression(int val) {
      fIsCompressionOverride = true;
      fCompression = val;
   }
   void SetCompression(RCompressionSetting::EAlgorithm algorithm, int compressionLevel) {
      SetCompression(CompressionSettings(algorithm, compressionLevel));
   }
   /// Returns true if the compression settings were explicitly set by the user.
   bool IsCompressionOverride() const { return fIsCompressionOverride; }

   ENTupleContainerFormat GetContainerFormat() const { return fContainerFormat; }
   void SetContainerFormat(ENTupleContainerFormat val) { fContainerFormat = val; }

   NTupleSize_t GetNElementsPerPage() const { return fNElementsPerPage; }
   void SetNElementsPerPage(NTupleSize_t val) { fNElementsPerPage = val; }

   NTupleSize_t GetNEntriesPerCluster() const { return fNEntriesPerCluster; }
   void SetNEntriesPerCluster(NTupleSize_t val) { fNEntriesPerCluster = val; }

   bool GetUseBufferedWrite() const { return fUseBufferedWrite; }
   void SetUseBufferedWrite(bool val) { fUseBufferedWrite = val; }
};

// clang-format off
/**
\class ROOT::Experimental::RNTupleWriteOptionsDaos
\ingroup NTuple
\brief DAOS-specific user-tunable settings for storing ntuples
*/
// clang-format on
class RNTupleWriteOptionsDaos : public RNTupleWriteOptions {
  std::string fObjectClass{"SX"};

public:
   virtual ~RNTupleWriteOptionsDaos() = default;
   std::unique_ptr<RNTupleWriteOptions> Clone() const override
   { return std::make_unique<RNTupleWriteOptionsDaos>(*this); }

   const std::string &GetObjectClass() const { return fObjectClass; }
   /// Set the object class used to generate OIDs that relate to user data. Any
   /// `OC_xxx` constant defined in `daos_obj_class.h` may be used here without
   /// the OC_ prefix.
   void SetObjectClass(const std::string &val) { fObjectClass = val; }
};

// clang-format off
/**
\class ROOT::Experimental::RNTupleReadOptions
\ingroup NTuple
\brief Common user-tunable settings for reading ntuples

All page source classes need to support the common options.
*/
// clang-format on
class RNTupleReadOptions {
public:
   enum EClusterCache {
      kOff,
      kOn,
      kDefault = kOn,
   };

private:
   EClusterCache fClusterCache = EClusterCache::kDefault;

public:
   EClusterCache GetClusterCache() const { return fClusterCache; }
   void SetClusterCache(EClusterCache val) { fClusterCache = val; }
};

} // namespace Experimental
} // namespace ROOT

#endif
