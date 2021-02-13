/// \file ROOT/RAxisConfig.hxx
/// \ingroup Hist ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2020-02-05
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAxisConfig
#define ROOT7_RAxisConfig

#include <string>
#include <utility>
#include <vector>

#include "ROOT/RStringView.hxx"

namespace ROOT {
namespace Experimental {

/**
\class RAxisConfig
Objects used to configure the different axis types. It can store the
properties of all ROOT-provided axis types, together with the type of the axis.

TODO: that's what a variant will be invented for!
*/
class RAxisConfig {
public:
   enum EKind {
      kEquidistant, ///< represents a RAxisEquidistant
      kGrow,        ///< represents a RAxisGrow
      kIrregular,   ///< represents a RAxisIrregular
      kLabels,      ///< represents a RAxisLabels
      kNumKinds
   };

private:
   std::string fTitle;
   int fNBinsNoOver;
   EKind fKind;                      ///< The kind of axis represented by this configuration
   std::vector<double> fBinBorders;  ///< Bin borders of the RAxisIrregular
   std::vector<std::string> fLabels; ///< Bin labels for a RAxisLabels

   /// Represents a `RAxisEquidistant` or `RAxisGrow` with `nbins` (excluding over- and
   /// underflow bins) from `from` to `to`, with an axis title.
   explicit RAxisConfig(std::string_view title, int nbins, double from, double to, EKind kind)
      : fTitle(title), fNBinsNoOver(nbins), fKind(kind), fBinBorders(2)
   {
      if (from > to)
         std::swap(to, from);

      fBinBorders[0] = from;
      fBinBorders[1] = to;
   }

public:
   /// Tag type signalling that an axis should be able to grow; used for calling
   /// the appropriate constructor.
   struct Grow_t {
   };
   /// Tag signalling that an axis should be able to grow; used for calling the
   /// appropriate constructor like so:
   ///     RAxisConfig ac(RAxisConfig::Grow, 10, 0., 1.);
   constexpr static const Grow_t Grow{};

   /// Represents a `RAxisEquidistant` with `nbins` from `from` to `to`, and
   /// axis title.
   RAxisConfig(std::string_view title, int nbins, double from, double to)
      : RAxisConfig(title, nbins, from, to, kEquidistant)
   {}

   /// Represents a `RAxisEquidistant` with `nbins` from `from` to `to`.
   RAxisConfig(int nbins, double from, double to): RAxisConfig("", nbins, from, to, kEquidistant) {}

   /// Represents a `RAxisGrow` with `nbins` from `from` to `to`, and axis title.
   RAxisConfig(std::string_view title, Grow_t, int nbins, double from, double to)
      : RAxisConfig(title, nbins, from, to, kGrow)
   {}

   /// Represents a `RAxisGrow` with `nbins` from `from` to `to`.
   RAxisConfig(Grow_t, int nbins, double from, double to): RAxisConfig("", nbins, from, to, kGrow) {}

   /// Represents a `RAxisIrregular` with `binborders` and title.
   RAxisConfig(std::string_view title, const std::vector<double> &binborders)
      : fTitle(title), fNBinsNoOver(binborders.size() - 1), fKind(kIrregular), fBinBorders(binborders)
   {}

   /// Represents a `RAxisIrregular` with `binborders`.
   RAxisConfig(const std::vector<double> &binborders): RAxisConfig("", binborders) {}

   /// Represents a `RAxisIrregular` with `binborders` and title.
   RAxisConfig(std::string_view title, std::vector<double> &&binborders) noexcept
      : fTitle(title), fNBinsNoOver(binborders.size() - 1), fKind(kIrregular),
        fBinBorders(std::move(binborders))
   {}

   /// Represents a `RAxisIrregular` with `binborders`.
   RAxisConfig(std::vector<double> &&binborders) noexcept: RAxisConfig("", std::move(binborders)) {}

   /// Represents a `RAxisLabels` with `labels` and title.
   RAxisConfig(std::string_view title, const std::vector<std::string_view> &labels)
      : fTitle(title), fNBinsNoOver(labels.size()), fKind(kLabels), fLabels(labels.begin(), labels.end())
   {}

   /// Represents a `RAxisLabels` with `labels`.
   RAxisConfig(const std::vector<std::string_view> &labels): RAxisConfig("", labels) {}

   /// Represents a `RAxisLabels` with `labels` and title.
   RAxisConfig(std::string_view title, std::vector<std::string> &&labels)
      : fTitle(title), fNBinsNoOver(labels.size()), fKind(kLabels), fLabels(std::move(labels))
   {}

   /// Represents a `RAxisLabels` with `labels`.
   RAxisConfig(std::vector<std::string> &&labels): RAxisConfig("", std::move(labels)) {}

   /// Get the axis's title
   const std::string &GetTitle() const { return fTitle; }

   /// Get the axis kind represented by this `RAxisConfig`.
   EKind GetKind() const noexcept { return fKind; }

   /// Get the number of bins, excluding under- and overflow.
   int GetNBinsNoOver() const noexcept { return fNBinsNoOver; }

   /// Get the bin borders; non-empty if the GetKind() == kIrregular.
   const std::vector<double> &GetBinBorders() const noexcept { return fBinBorders; }

   /// Get the bin labels; non-empty if the GetKind() == kLabels.
   const std::vector<std::string> &GetBinLabels() const noexcept { return fLabels; }
};

namespace Internal {

/// Converts a RAxisConfig of whatever kind to the corresponding RAxisBase-derived
/// object.
template <RAxisConfig::EKind>
struct AxisConfigToType; // Only specializations are defined.

} // namespace Internal
} // namespace Experimental
} // namespace ROOT

#endif // ROOT7_RAxisConfig header guard