// @(#)root/zip:$Id$
// Author: David Dagenhart   May 2011

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_Compression
#define ROOT_Compression

#include "RtypesCore.h"

namespace ROOT {

/// The global settings depend on a global variable named R__ZipMode which can be
/// modified by a global function named R__SetZipMode. Both are defined in Bits.h.
///
///  - The default is to use the global setting and the default of the global
///    setting is to use the ZLIB compression algorithm.
///  - The LZMA algorithm (from the XZ package) is also available. The LZMA
///    compression usually results in greater compression factors, but takes
///    more CPU time and memory when compressing. LZMA memory usage is particularly
///    high for compression levels 8 and 9.
///  - Finally, the LZ4 package results in worse compression ratios
///    than ZLIB but achieves much faster decompression rates.
///
/// The current algorithms support level 1 to 9. The higher the level the greater
/// the compression and more CPU time and memory resources used during compression.
/// Level 0 means no compression.
///
/// Recommendation for the compression algorithm's levels:
///  - ZLIB is recommended to be used with compression level 1 [101]
///  - LZMA is recommended to be used with compression level 7-8 (higher is better,
///   since in the case of LZMA we don't care about compression/decompression speed)
///   [207 - 208]
///  - LZ4 is recommended to be used with compression level 4 [404]
///  - ZSTD is recommended to be used with compression level 5 [505]

struct RCompressionSetting {
   struct EDefaults { /// Note: this is only temporarily a struct and will become a enum class hence the name convention
                      /// used.
      enum EValues {
         /// Use the global compression setting for this process; may be affected by rootrc.
         kUseGlobal = 0,
         /// Use the compile-time default setting
         kUseCompiledDefault = 101,
         /// Use the default analysis setting; fast reading but poor compression ratio
         kUseAnalysis = 404,
         /// Use the new recommended general-purpose setting; it is a best trade-off between compression ratio/decompression speed
         kUseGeneralPurpose = 505,
         /// Use the setting that results in the smallest files; very slow read and write
         kUseSmallest = 207,
      };
   };
   struct ELevel { /// Note: this is only temporarily a struct and will become a enum class hence the name convention
                   /// used.
      enum EValues {
         /// Some objects use this value to denote that the compression algorithm
         /// should be inherited from the parent object
         kInherit = -1,
         /// Compression level reserved for "uncompressed state"
         kUncompressed = 0,
         /// Compression level reserved when we are not sure what to use (1 is for the fastest compression)
         kUseMin = 1,
         /// Compression level reserved for ZLIB compression algorithm (fastest compression)
         kDefaultZLIB = 1,
         /// Compression level reserved for LZ4 compression algorithm (trade-off between file ratio/decompression speed)
         kDefaultLZ4 = 4,
         /// Compression level reserved for ZSYD compression algorithm (trade-off between file ratio/decompression speed)
         kDefaultZSTD = 5,
         /// Compression level reserved for old ROOT compression algorithm
         kDefaultOld = 6,
         /// Compression level reserved for LZMA compression algorithm (slowest compression with smallest files)
         kDefaultLZMA = 7
      };
   };
   struct EAlgorithm { /// Note: this is only temporarily a struct and will become a enum class hence the name
                        /// convention used.
      enum EValues {
         /// Some objects use this value to denote that the compression algorithm
         /// should be inherited from the parent object (e.g., TBranch should get the algorithm from the TTree)
         kInherit = -1,
         /// Use the global compression algorithm
         kUseGlobal = 0,
         /// Use ZLIB compression
         kZLIB,
         /// Use LZMA compression
         kLZMA,
         /// Use the old compression algorithm
         kOldCompressionAlgo,
         /// Use LZ4 compression
         kLZ4,
         /// Use ZSTD compression
         kZSTD,
         /// Undefined compression algorithm (must be kept the last of the list in case a new algorithm is added).
         kUndefined
      };
   };
};

enum ECompressionAlgorithm {
   /// Deprecated name, do *not* use:
   kUseGlobalCompressionSetting = RCompressionSetting::EAlgorithm::kUseGlobal,
   /// Deprecated name, do *not* use:
   kUseGlobalSetting = RCompressionSetting::EAlgorithm::kUseGlobal,
   /// Deprecated name, do *not* use:
   kZLIB = RCompressionSetting::EAlgorithm::kZLIB,
   /// Deprecated name, do *not* use:
   kLZMA = RCompressionSetting::EAlgorithm::kLZMA,
   /// Deprecated name, do *not* use:
   kOldCompressionAlgo = RCompressionSetting::EAlgorithm::kOldCompressionAlgo,
   /// Deprecated name, do *not* use:
   kLZ4 = RCompressionSetting::EAlgorithm::kLZ4,
   /// Deprecated name, do *not* use:
   kZSTD = RCompressionSetting::EAlgorithm::kZSTD,
   /// Deprecated name, do *not* use:
   kUndefinedCompressionAlgorithm = RCompressionSetting::EAlgorithm::kUndefined
};

int CompressionSettings(RCompressionSetting::EAlgorithm algorithm, int compressionLevel);
/// Deprecated name, do *not* use:
int CompressionSettings(ROOT::ECompressionAlgorithm algorithm, int compressionLevel);
} // namespace ROOT

#endif
