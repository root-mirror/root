// Author: Brian Bockelman UNL 09/2017

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TIO_FEATURES
#define ROOT_TIO_FEATURES

#include "Rtypes.h"

class TBasket;
class TBranch;
class TTree;

// keep it here to have a note that was removed
// #ifndef R__LESS_INCLUDES
// #include "TBasket.h"
// #include <memory>
// #include <sstream>
// #include <initializer_list>
// #include <typeinfo>
// #include <type_traits> // is_same, enable_if
// #endif

namespace ROOT {


// These are the known, supported, and enabled-by-default features for ROOT IO.
//
// Note that the `kSupported` members for EIOFeatures, Experimental::EIOFeatures, and
// Experiment::EIOUnsupportedFeatures should have no intersection and a union of equal
// to BITS(kIOFeatureCount).
//
enum class EIOFeatures {
   kSupported = 0  // Union of all known, supported, and enabled-by-default features (currently none).
};


namespace Experimental {

// These are the known and supported "experimental" features, not enabled by default.
// When these are enabled by default, they will move to `ROOT::EIOFeatures`.
//
// Note that these all show up in TBasket::EIOBits, but it is desired to have the enum be at
// the "ROOT-IO-wide" level and not restricted to TBasket -- even if all the currently-foreseen
// usage of this mechanism somehow involves baskets currently.
enum class EIOFeatures {
   kGenerateOffsetMap = BIT(0),
   kSupported = kGenerateOffsetMap  // Union of all features in this enum.
};


// These are previous experimental features that are not supported in this series.
// NOTE: the intent is that there is never an IO feature that goes into the ROOT:: namespace
// but is unsupported.
enum class EIOUnsupportedFeatures {
   kUnsupported = 0  // Union of all features in this enum.
};


}  // namespace Experimental


class TIOFeatures {
friend class ::TTree;
friend class ::TBranch;
friend class ::TBasket;

public:
   TIOFeatures() {}

   void Clear(EIOFeatures bits);
   void Clear(Experimental::EIOUnsupportedFeatures bits);
   void Clear(Experimental::EIOFeatures bits);
   bool Set(EIOFeatures bits);
   bool Set(Experimental::EIOFeatures bits);
   bool Set(const std::string &);
   bool Test(EIOFeatures bits) const;
   bool Test(Experimental::EIOFeatures bits) const;
   bool Test(Experimental::EIOUnsupportedFeatures bits) const;
   void Print() const;

   // The number of known, defined IO features (supported / unsupported / experimental).
   static constexpr int kIOFeatureCount = 1;

private:
   // These methods allow access to the raw bitset underlying
   // this object, breaking type safety.  They are necessary for
   // efficient interaction with TTree / TBranch / TBasket, but left
   // private to prevent users from interacting with the raw bits.
   TIOFeatures(UChar_t IOBits) : fIOBits(IOBits) {}
   UChar_t GetFeatures() const;
   void Set(UChar_t newBits) {fIOBits = newBits;}

   UChar_t fIOBits{0};
};

}  // namespace ROOT

#endif  // ROOT_TIO_FEATURES
