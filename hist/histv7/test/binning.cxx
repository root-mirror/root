#include "gtest/gtest.h"
#include "ROOT/RAxis.hxx"
#include "ROOT/RHistData.hxx"
#include "ROOT/RHistImpl.hxx"
#include <array>
#include <sstream>
#include <vector>

using namespace ROOT::Experimental;

// Convenience tool for enumerating bins in row-major order along an axis,
// collecting some useful bin properties along the way.
struct BinProperties {
   int index;
   double from, center, to;
};

std::vector<BinProperties> GetAllBinProperties(const RAxisBase& axis) {
   std::vector<BinProperties> result;
   result.reserve(axis.GetNBins());

   auto push_bin = [&axis, &result](int bin) {
      result.push_back({ bin,
                         axis.GetBinFrom(bin),
                         axis.GetBinCenter(bin),
                         axis.GetBinTo(bin) });
   };

   if (!axis.CanGrow()) push_bin(-1);
   for (int bin = 1; bin <= axis.GetNBinsNoOver(); ++bin) push_bin(bin);
   if (!axis.CanGrow()) push_bin(-2);

   return result;
}

// Generic test that an N-dimensional histogram with a certain axis
// configuration is binned correctly.
template <typename... Axes>
void TestHistogramBinning(Axes&&... axes) {
   // Query everything there is to know about the axis bins
   static constexpr std::size_t NDIMS = sizeof...(axes);
   std::array<std::vector<BinProperties>, NDIMS>
      axis_bin_properties { GetAllBinProperties(axes)... };

   // Build an RHistImpl with the provided axis configuration
   Detail::RHistImpl<Detail::RHistData<NDIMS,
                                       double,
                                       std::vector<double>,
                                       RHistStatContent>,
                     Axes...>
      hist(std::forward<Axes>(axes)...);

   // Track the index of the last regular and overflow bin that was seen
   int last_regular_bin = 0;
   int last_overflow_bin = 0;

   // Prepare to iterate over bins
   using BinPropertiesIter = std::vector<BinProperties>::const_iterator;
   std::array<BinPropertiesIter, NDIMS> local_bins_iters, local_bins_ends;
   for (std::size_t axis = 0; axis < NDIMS; ++axis) {
      local_bins_iters[axis] = axis_bin_properties[axis].cbegin();
      local_bins_ends[axis] = axis_bin_properties[axis].cend();
   }

   // Do the bin iteration in row-major order
   while (local_bins_iters.back() != local_bins_ends.back()) {
      // Check various properties of the current bin
      std::array<double, NDIMS> bin_center, bin_from, bin_to;
      std::array<int, NDIMS> local_bin_indices;
      for (std::size_t axis = 0; axis < NDIMS; ++axis) {
         local_bin_indices[axis] = local_bins_iters[axis]->index;
         bin_center[axis] = local_bins_iters[axis]->center;
         bin_from[axis] = local_bins_iters[axis]->from;
         bin_to[axis] = local_bins_iters[axis]->to;
      }

      // Make googletest report the current bin index on failure
      std::stringstream sstream("On histogram bin [ ",
                                std::ios_base::out | std::ios_base::ate);
      for (std::size_t axis = 0; axis < NDIMS; ++axis) {
         sstream << local_bin_indices[axis] << ' ';
      }
      sstream << ']';
      SCOPED_TRACE(sstream.str());

      // Check the global bin index
      const int global_bin = hist.GetBinIndexFromLocalBins(local_bin_indices);
      if (std::all_of(local_bin_indices.cbegin(), local_bin_indices.cend(),
                      [](int bin) { return bin >= 1; })) {
         // This is a regular bin, it should be attributed to the next positive
         // index after the previously observed regular bin index since we're
         // iterating in row-major order.
         EXPECT_EQ(global_bin, ++last_regular_bin);
      } else {
         // This is an overflow bin, it should be attributed to the next
         // negative index after the previously observed overflow bin since
         // we're iterating in row-major order.
         EXPECT_EQ(global_bin, --last_overflow_bin);
      }

      // Check the local bin coordinates
      EXPECT_EQ(hist.GetLocalBins(global_bin), local_bin_indices);

      EXPECT_EQ(hist.GetBinCenter(global_bin), bin_center);
      EXPECT_EQ(hist.GetBinFrom(global_bin), bin_from);
      EXPECT_EQ(hist.GetBinTo(global_bin), bin_to);

      // Go to the next bin in row-major order
      for (std::size_t axis = 0; axis < NDIMS; ++axis) {
         ++local_bins_iters[axis];
         if (local_bins_iters[axis] != local_bins_ends[axis]) {
            break;
         } else if (axis != NDIMS-1) {
            local_bins_iters[axis] = axis_bin_properties[axis].cbegin();
            continue;
         }
      }
   }
}

// ===

// 1D histograms should be binned just like their underlying axis

//   UF  Reg1  Reg2  ...  RegN  OF
//  -------------------------------
//   -1    1     2          N   -2
TEST(HistImplBinning, Eq) {
   SCOPED_TRACE("1D histogram with equidistant axis");
   TestHistogramBinning(RAxisEquidistant(6, -7.5, 5.8));
}

//   Reg1  Reg2  ...  RegN
//  -----------------------
//     1     2          N 
TEST(HistImplBinning, Grow) {
   SCOPED_TRACE("1D histogram with growable axis");
   TestHistogramBinning(RAxisGrow(3, 3.0, 5.3));
}

// ===

// 2D+ histograms should be binned in row-major order

//                     Axis 0
//               UF  Reg1  Reg2  OF
//        ---------------------------
//     A   UF  | -1   -2    -3   -4
//     x  Reg1 | -5    1     2   -6
//     .  Reg2 | -7    3     4   -8
//        Reg3 | -9    5     6   -10
//     1   OF  | -11 -12   -13   -14
TEST(HistImplBinning, EqEq) {
   SCOPED_TRACE("2D histogram with equidistant-equidistant axes");
   TestHistogramBinning(RAxisEquidistant(8, -9.5, 4.7),
                        RAxisEquidistant(5, -3.2, -2.5));
}

//                     Axis 0
//               UF  Reg1  Reg2  OF
//        --------------------------
//     A  Reg1 | -1    1     2   -2
//     x  Reg2 | -3    3     4   -4
//     1  Reg3 | -5    5     6   -6
TEST(HistImplBinning, EqGrow) {
   SCOPED_TRACE("2D histogram with equidistant-growable axes");
   TestHistogramBinning(RAxisEquidistant(3, -5.2, 3.1),
                        RAxisGrow(3, 3.9, 4.4));
}

//                Axis 0
//              Reg1  Reg2
//        ----------------
//     A   UF  | -1    -2
//     x  Reg1 |  1     2
//     .  Reg2 |  3     4
//        Reg3 |  5     6
//     1   OF  | -3    -4 
TEST(HistImplBinning, GrowEq) {
   SCOPED_TRACE("2D histogram with growable-equidistant axes");
   TestHistogramBinning(RAxisGrow(2, -7.9, 7.5),
                        RAxisEquidistant(8, 3.1, 3.3));
}

//                 Axis 0
//               Reg1  Reg2
//        ------------------
//     A  Reg1 |   1     2
//     x  Reg2 |   3     4
//     1  Reg3 |   5     6
TEST(HistImplBinning, GrowGrow) {
   SCOPED_TRACE("2D histogram with growable-growable axes");
   TestHistogramBinning(RAxisGrow(5, -7.2, -2.1),
                        RAxisGrow(5, 2.9, 9.6));
}

// ===

// For 3D histograms, comment ASCII art illustrates the global binning of
// successive planes as one iterates over axis 2. Each plane is represented
// using the same convention that was used for 2D histograms.

//   -1   -2   -3   -4    |   -21  -22  -23  -24   |   -35  -36  -37  -38
//   -5   -6   -7   -8    |   -25   1    2   -26   |   -39  -40  -41  -42
//   -9   -10  -11  -12   |   -27   3    4   -28   |   -43  -44  -45  -46
//   -13  -14  -15  -16   |   -29   5    6   -30   |   -47  -48  -49  -50
//   -17  -18  -19  -20   |   -31  -32  -33  -34   |   -51  -52  -53  -54
TEST(HistImplBinning, EqEqEq) {
   SCOPED_TRACE("3D histogram with equidistant-equidistant-equidistant axes");
   TestHistogramBinning(RAxisEquidistant(6, -2.2, 9.3),
                        RAxisEquidistant(4, -9.3, -7.4),
                        RAxisEquidistant(5, -5.7, 3.2));
}

//   -1   -2   -3   -4
//   -5    1    2   -6
//   -7    3    4   -8
//   -9    5    6   -10
//   -11  -12  -13  -14
TEST(HistImplBinning, EqEqGrow) {
   SCOPED_TRACE("3D histogram with equidistant-equidistant-growable axes");
   TestHistogramBinning(RAxisEquidistant(7, -7.8, -2.4),
                        RAxisEquidistant(6, 2.0, 2.5),
                        RAxisGrow(7, -4.5, -3.1));
}

//   -1   -2   -3   -4    |   -13   1    2   -14   |   -19  -20  -21  -22
//   -5   -6   -7   -8    |   -15   3    4   -16   |   -23  -24  -25  -26
//   -9   -10  -11  -12   |   -17   5    6   -18   |   -27  -28  -29  -30
TEST(HistImplBinning, EqGrowEq) {
   SCOPED_TRACE("3D histogram with equidistant-growable-equidistant axes");
   TestHistogramBinning(RAxisEquidistant(9, -4.5, 2.1),
                        RAxisGrow(5, -7.3, -5.5),
                        RAxisEquidistant(3, -8.8, 3.6));
}

//   -1    1    2   -2
//   -3    3    4   -4
//   -5    5    6   -6
TEST(HistImplBinning, EqGrowGrow) {
   SCOPED_TRACE("3D histogram with equidistant-growable-growable axes");
   TestHistogramBinning(RAxisEquidistant(7, 4.8, 7.8),
                        RAxisGrow(2, -3.7, 4.8),
                        RAxisGrow(9, 4.0, 6.7));
}

//   -1   -2    |   -11  -12   |  -15  -16
//   -3   -4    |    1    2    |  -17  -18
//   -5   -6    |    3    4    |  -19  -20
//   -7   -8    |    5    6    |  -21  -22
//   -9   -10   |   -13  -14   |  -23  -24
TEST(HistImplBinning, GrowEqEq) {
   SCOPED_TRACE("3D histogram with growable-equidistant-equidistant axes");
   TestHistogramBinning(RAxisGrow(2, -7.8, 8.5),
                        RAxisEquidistant(3, -8.7, -3.4),
                        RAxisEquidistant(9, 1.7, 3.3));
}

//   -1   -2
//    1    2
//    3    4
//    5    6
//   -3   -4
TEST(HistImplBinning, GrowEqGrow) {
   SCOPED_TRACE("3D histogram with growable-equidistant-growable axes");
   TestHistogramBinning(RAxisGrow(8, 0.6, 1.0),
                        RAxisEquidistant(2, -1.8, 2.5),
                        RAxisGrow(4, -1.9, 4.0));
}

//   -1   -2    |    1    2    |  -7   -8
//   -3   -4    |    3    4    |  -9   -10
//   -5   -6    |    5    6    |  -11  -12
TEST(HistImplBinning, GrowGrowEq) {
   SCOPED_TRACE("3D histogram with growable-growable-equidistant axes");
   TestHistogramBinning(RAxisGrow(3, -8.2, 0.0),
                        RAxisGrow(6, -4.8, 2.5),
                        RAxisEquidistant(6, -3.9, -2.6));
}

//   1    2
//   3    4
//   5    6
TEST(HistImplBinning, GrowGrowGrow) {
   SCOPED_TRACE("3D histogram with growable-growable-growable axes");
   TestHistogramBinning(RAxisGrow(5, -1.7, 9.6),
                        RAxisGrow(9, -6.1, 8.7),
                        RAxisGrow(9, -4.9, 7.6));
}
