/// \file ROOT/RColumnModel.hxx
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

#ifndef ROOT7_RColumnModel
#define ROOT7_RColumnModel

#include <ROOT/RStringView.hxx>

namespace ROOT {
namespace Experimental {

// clang-format off
/**
\class ROOT::Experimental::EColumnType
\ingroup Forest
\brief The possible trivial, native content types of a column
*/
// clang-format on
enum class EColumnType {
   kUnknown,
   kIndex, // type for root columns of nested collections
   kByte,
   kReal64,
   kReal32,
   kReal16,
   kReal8,
   kInt64,
   kInt32,
   kInt16,
   //...
};

// clang-format off
/**
\class ROOT::Experimental::RColumnModel
\ingroup Forest
\brief Holds the static meta-data of a column in a tree
*/
// clang-format on
class RColumnModel {
public:
   RColumnModel(std::string_view name, EColumnType type, bool is_sorted);
};

} // namespace Experimental
} // namespace ROOT

#endif
