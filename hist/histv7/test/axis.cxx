#include "gtest/gtest.h"
#include "ROOT/RAxis.hxx"

using namespace ROOT::Experimental;

// FIXME: Passing this directly as an initializer list is ambiguous because
//        the compiler doesn't know if it should convert the inner const char*
//        literals to std::string_view or std::string.
std::vector<std::string_view> labels{"abc", "de", "fghi", "j", "klmno"};

// Test RAxisConfig and conversion to concrete axis types
TEST(AxisTest, Config) {
  // Equidistant
  {
    auto test = [](const RAxisConfig& cfg, std::string_view title) {
      EXPECT_EQ(cfg.GetTitle(), title);
      EXPECT_EQ(cfg.GetNBinsNoOver(), 10);
      EXPECT_EQ(cfg.GetKind(), RAxisConfig::kEquidistant);
      EXPECT_EQ(cfg.GetBinBorders().size(), 2u);
      EXPECT_EQ(cfg.GetBinBorders()[0], 1.2);
      EXPECT_EQ(cfg.GetBinBorders()[1], 3.4);
      EXPECT_EQ(cfg.GetBinLabels().size(), 0u);

      RAxisEquidistant axis = Internal::AxisConfigToType<RAxisConfig::kEquidistant>()(cfg);
      EXPECT_EQ(axis.GetTitle(), title);
      EXPECT_EQ(axis.GetNBinsNoOver(), 10);
      EXPECT_EQ(axis.GetMinimum(), 1.2);
      EXPECT_DOUBLE_EQ(axis.GetMaximum(), 3.4);
    };

    {
      SCOPED_TRACE("Equidistant axis config w/o title");
      test({10, 1.2, 3.4}, "");
    }

    {
      SCOPED_TRACE("Equidistant axis config with title");
      test({"RITLE_E", 10, 1.2, 3.4}, "RITLE_E");
    }
  }

  // Growable
  {
    auto test = [](const RAxisConfig& cfg, std::string_view title) {
      EXPECT_EQ(cfg.GetTitle(), title);
      EXPECT_EQ(cfg.GetNBinsNoOver(), 10);
      EXPECT_EQ(cfg.GetKind(), RAxisConfig::kGrow);
      EXPECT_EQ(cfg.GetBinBorders().size(), 2u);
      EXPECT_EQ(cfg.GetBinBorders()[0], 1.2);
      EXPECT_EQ(cfg.GetBinBorders()[1], 3.4);
      EXPECT_EQ(cfg.GetBinLabels().size(), 0u);

      RAxisGrow axis = Internal::AxisConfigToType<RAxisConfig::kGrow>()(cfg);
      EXPECT_EQ(axis.GetTitle(), title);
      EXPECT_EQ(axis.GetNBinsNoOver(), 10);
      EXPECT_EQ(axis.GetMinimum(), 1.2);
      EXPECT_DOUBLE_EQ(axis.GetMaximum(), 3.4);
    };

    {
      SCOPED_TRACE("Growable axis config w/o title");
      test({RAxisConfig::Grow, 10, 1.2, 3.4}, "");
    }

    {
      SCOPED_TRACE("Growable axis config with title");
      test({"RITLE_G", RAxisConfig::Grow, 10, 1.2, 3.4}, "RITLE_G");
    }
  }

  // Irregular
  {
    auto test = [](const RAxisConfig& cfg, std::string_view title) {
      EXPECT_EQ(cfg.GetTitle(), title);
      EXPECT_EQ(cfg.GetNBinsNoOver(), 3);
      EXPECT_EQ(cfg.GetKind(), RAxisConfig::kIrregular);
      EXPECT_EQ(cfg.GetBinBorders().size(), 4u);
      EXPECT_EQ(cfg.GetBinBorders()[0], 2.3);
      EXPECT_EQ(cfg.GetBinBorders()[1], 5.7);
      EXPECT_EQ(cfg.GetBinBorders()[2], 11.13);
      EXPECT_EQ(cfg.GetBinBorders()[3], 17.19);
      EXPECT_EQ(cfg.GetBinLabels().size(), 0u);

      RAxisIrregular axis = Internal::AxisConfigToType<RAxisConfig::kIrregular>()(cfg);
      EXPECT_EQ(axis.GetTitle(), title);
      EXPECT_EQ(axis.GetNBinsNoOver(), 3);
      EXPECT_EQ(axis.GetBinBorders().size(), 4u);
      EXPECT_EQ(axis.GetBinBorders()[0], 2.3);
      EXPECT_EQ(axis.GetBinBorders()[1], 5.7);
      EXPECT_EQ(axis.GetBinBorders()[2], 11.13);
      EXPECT_EQ(axis.GetBinBorders()[3], 17.19);
    };

    {
      SCOPED_TRACE("Irregular axis config w/o title");
      test({{2.3, 5.7, 11.13, 17.19}}, "");
    }

    {
      SCOPED_TRACE("Irregular axis config with title");
      test({"RITLE_I", {2.3, 5.7, 11.13, 17.19}}, "RITLE_I");
    }
  }

  // Labels
  {
    auto test = [](const RAxisConfig& cfg, std::string_view title) {
      EXPECT_EQ(cfg.GetTitle(), title);
      EXPECT_EQ(cfg.GetNBinsNoOver(), 5);
      EXPECT_EQ(cfg.GetKind(), RAxisConfig::kLabels);
      EXPECT_EQ(cfg.GetBinBorders().size(), 0u);
      EXPECT_EQ(cfg.GetBinLabels().size(), 5u);
      EXPECT_EQ(cfg.GetBinLabels()[0], "abc");
      EXPECT_EQ(cfg.GetBinLabels()[1], "de");
      EXPECT_EQ(cfg.GetBinLabels()[2], "fghi");
      EXPECT_EQ(cfg.GetBinLabels()[3], "j");
      EXPECT_EQ(cfg.GetBinLabels()[4], "klmno");

      RAxisLabels axis = Internal::AxisConfigToType<RAxisConfig::kLabels>()(cfg);
      EXPECT_EQ(axis.GetTitle(), title);
      EXPECT_EQ(axis.GetNBinsNoOver(), 5);
      EXPECT_EQ(axis.GetBinLabels().size(), 5u);
      EXPECT_EQ(axis.GetBinLabels()[0], "abc");
      EXPECT_EQ(axis.GetBinLabels()[1], "de");
      EXPECT_EQ(axis.GetBinLabels()[2], "fghi");
      EXPECT_EQ(axis.GetBinLabels()[3], "j");
      EXPECT_EQ(axis.GetBinLabels()[4], "klmno");
    };

    {
      SCOPED_TRACE("Labeled axis config w/o title");
      test({labels}, "");
    }

    {
      SCOPED_TRACE("Labeled axis config with title");
      test({"RITLE_L", labels}, "RITLE_L");
    }
  }
}

TEST(AxisTest, Iterator) {
  auto it = RAxisBase::const_iterator(42);
  EXPECT_EQ(*it, 42);

  {
    auto it2 = ++it;
    EXPECT_EQ(*it, 43);
    EXPECT_EQ(*it2, 43);
    auto it3 = --it;
    EXPECT_EQ(*it, 42);
    EXPECT_EQ(*it3, 42);
  }

  {
    auto it2 = it++;
    EXPECT_EQ(*it, 43);
    EXPECT_EQ(*it2, 42);
    auto it3 = it--;
    EXPECT_EQ(*it, 42);
    EXPECT_EQ(*it3, 43);
  }

  {
    auto it2 = (it += 7);
    EXPECT_EQ(*it, 49);
    EXPECT_EQ(*it2, 49);
    auto it3 = (it -= 7);
    EXPECT_EQ(*it, 42);
    EXPECT_EQ(*it3, 42);
  }

  {
    auto it2 = it + 7;
    EXPECT_EQ(*it, 42);
    EXPECT_EQ(*it2, 49);
    auto it3 = 7 + it;
    EXPECT_EQ(*it, 42);
    EXPECT_EQ(*it3, 49);
    auto it4 = it - 7;
    EXPECT_EQ(*it, 42);
    EXPECT_EQ(*it4, 35);
  }

  {
    auto it2 = RAxisBase::const_iterator(54);
    EXPECT_EQ(it2 - it, 12);
    EXPECT_EQ(it[12], 54);
  }

  {
    auto it_m1 = RAxisBase::const_iterator(41);
    auto it_p1 = RAxisBase::const_iterator(43);

    EXPECT_EQ(it < it_m1, false);
    EXPECT_EQ(it < it, false);
    EXPECT_EQ(it < it_p1, true);

    EXPECT_EQ(it > it_m1, true);
    EXPECT_EQ(it > it, false);
    EXPECT_EQ(it > it_p1, false);

    EXPECT_EQ(it <= it_m1, false);
    EXPECT_EQ(it <= it, true);
    EXPECT_EQ(it <= it_p1, true);

    EXPECT_EQ(it >= it_m1, true);
    EXPECT_EQ(it >= it, true);
    EXPECT_EQ(it >= it_p1, false);

    EXPECT_EQ(it == it_m1, false);
    EXPECT_EQ(it == it, true);
    EXPECT_EQ(it == it_p1, false);

    EXPECT_EQ(it != it_m1, true);
    EXPECT_EQ(it != it, false);
    EXPECT_EQ(it != it_p1, true);
  }
}

// Common test items for RAxisBase child classes
void test_axis_base(const RAxisBase& axis,
                    std::string_view title,
                    bool can_grow,
                    int n_bins_no_over,
                    double minimum,
                    double maximum) {
  EXPECT_EQ(axis.GetTitle(), title);
  EXPECT_EQ(axis.CanGrow(), can_grow);
  EXPECT_EQ(axis.GetNBinsNoOver(), n_bins_no_over);

  const int n_overflow_bins = can_grow ? 0 : 2;
  EXPECT_EQ(axis.GetNOverflowBins(), n_overflow_bins);
  EXPECT_EQ(axis.GetNBins(), n_bins_no_over + n_overflow_bins);

  const int kInvalidBin = RAxisBase::kInvalidBin;
  const int underflow_bin = can_grow ? kInvalidBin : -1;
  EXPECT_EQ(axis.GetUnderflowBin(), underflow_bin);

  const int overflow_bin = can_grow ? kInvalidBin : -2;
  EXPECT_EQ(axis.GetOverflowBin(), overflow_bin);

  EXPECT_EQ(axis.GetFirstBin(), 1);
  EXPECT_EQ(axis.GetLastBin(), n_bins_no_over);

  EXPECT_EQ(*axis.begin(), can_grow ? underflow_bin+1 : underflow_bin+2);
  EXPECT_EQ(*axis.end(), n_bins_no_over + 1);

  int nBins = 0;
  for (auto iter = axis.begin(); iter != axis.end(); ++iter) {
    ++nBins;
  }
  EXPECT_EQ(nBins, n_bins_no_over);

  EXPECT_DOUBLE_EQ(axis.GetMinimum(), minimum);
  EXPECT_DOUBLE_EQ(axis.GetMaximum(), maximum);
}

// Common test items for RAxisEquidistant child classes
void test_axis_equidistant(const RAxisEquidistant& axis,
                           std::string_view title,
                           bool can_grow,
                           int n_bins_no_over,
                           double minimum,
                           double maximum) {
  test_axis_base(axis, title, can_grow, n_bins_no_over, minimum, maximum);

  const double bin_width = (maximum - minimum) / n_bins_no_over;
  EXPECT_DOUBLE_EQ(axis.GetBinWidth(), bin_width);
  EXPECT_DOUBLE_EQ(axis.GetInverseBinWidth(), 1.0/bin_width);

  const int kInvalidBin = RAxisBase::kInvalidBin;
  const int underflow_findbin_res = can_grow ? kInvalidBin : -1;
  EXPECT_EQ(axis.FindBin(std::numeric_limits<double>::lowest()),
            underflow_findbin_res);
  EXPECT_EQ(axis.FindBin(minimum-0.01*bin_width), underflow_findbin_res);
  const int first_bin = 1;
  EXPECT_EQ(axis.FindBin(minimum+0.01*bin_width), first_bin);
  EXPECT_EQ(axis.FindBin(minimum+0.99*bin_width), first_bin);
  EXPECT_EQ(axis.FindBin(minimum+1.01*bin_width), first_bin+1);
  const int last_bin = first_bin + n_bins_no_over - 1;
  EXPECT_EQ(axis.FindBin(maximum-0.01*bin_width), last_bin);
  const int overflow_findbin_res = can_grow ? kInvalidBin : -2;
  EXPECT_EQ(axis.FindBin(maximum+0.01*bin_width), overflow_findbin_res);
  EXPECT_EQ(axis.FindBin(std::numeric_limits<double>::max()),
            overflow_findbin_res);

  // NOTE: Result of GetBinFrom on underflow bins, GetBinTo on overflow bins and
  //       GetBinCenter on either is considered unspecified for now. If we do
  //       ultimately decide to specify this behavior, please add a test here.
  if (!can_grow) {
    EXPECT_DOUBLE_EQ(axis.GetBinTo(0), minimum);
  }
  EXPECT_DOUBLE_EQ(axis.GetBinFrom(first_bin), minimum);
  EXPECT_DOUBLE_EQ(axis.GetBinCenter(first_bin), minimum+0.5*bin_width);
  EXPECT_DOUBLE_EQ(axis.GetBinTo(first_bin), minimum+bin_width);
  EXPECT_DOUBLE_EQ(axis.GetBinFrom(first_bin+1), minimum+bin_width);
  EXPECT_DOUBLE_EQ(axis.GetBinCenter(first_bin+1), minimum+1.5*bin_width);
  EXPECT_DOUBLE_EQ(axis.GetBinTo(first_bin+1), minimum+2*bin_width);
  EXPECT_DOUBLE_EQ(axis.GetBinFrom(last_bin), maximum-bin_width);
  EXPECT_DOUBLE_EQ(axis.GetBinCenter(last_bin), maximum-0.5*bin_width);
  EXPECT_DOUBLE_EQ(axis.GetBinTo(last_bin), maximum);
  if (!can_grow) {
    EXPECT_DOUBLE_EQ(axis.GetBinFrom(n_bins_no_over+1), maximum);
  }

  EXPECT_EQ(axis.GetBinIndexForLowEdge(std::numeric_limits<double>::lowest()),
            kInvalidBin);
  EXPECT_EQ(axis.GetBinIndexForLowEdge(minimum-bin_width), kInvalidBin);
  EXPECT_EQ(axis.GetBinIndexForLowEdge(minimum-0.5*bin_width), kInvalidBin);
  EXPECT_EQ(axis.GetBinIndexForLowEdge(minimum), first_bin);
  EXPECT_EQ(axis.GetBinIndexForLowEdge(minimum+0.5*bin_width), kInvalidBin);
  EXPECT_EQ(axis.GetBinIndexForLowEdge(minimum+bin_width), first_bin+1);
  EXPECT_EQ(axis.GetBinIndexForLowEdge(maximum-1.5*bin_width), kInvalidBin);
  EXPECT_EQ(axis.GetBinIndexForLowEdge(maximum-bin_width), last_bin);
  EXPECT_EQ(axis.GetBinIndexForLowEdge(maximum-0.5*bin_width), kInvalidBin);
  EXPECT_EQ(axis.GetBinIndexForLowEdge(maximum), last_bin+1);
  EXPECT_EQ(axis.GetBinIndexForLowEdge(maximum+0.5*bin_width), kInvalidBin);
  EXPECT_EQ(axis.GetBinIndexForLowEdge(maximum+bin_width), kInvalidBin);
  EXPECT_EQ(axis.GetBinIndexForLowEdge(std::numeric_limits<double>::max()),
            kInvalidBin);
}

TEST(AxisTest, Equidistant) {
  auto test = [](const RAxisEquidistant& axis,
                 std::string_view title,
                 int nbins,
                 double min,
                 double max) {
    test_axis_equidistant(axis, title, false, nbins, min, max);

    RAxisConfig cfg(axis);
    EXPECT_EQ(cfg.GetTitle(), title);
    EXPECT_EQ(cfg.GetNBinsNoOver(), nbins);
    EXPECT_EQ(cfg.GetKind(), RAxisConfig::kEquidistant);
    EXPECT_EQ(cfg.GetBinBorders().size(), 2u);
    EXPECT_EQ(cfg.GetBinBorders()[0], min);
    EXPECT_DOUBLE_EQ(cfg.GetBinBorders()[1], max);
    EXPECT_EQ(cfg.GetBinLabels().size(), 0u);
  };

  {
    SCOPED_TRACE("Equidistant axis w/o title, normal binning");
    test(RAxisEquidistant(10, 1.2, 3.4), "", 10, 1.2, 3.4);
  }

  {
    SCOPED_TRACE("Equidistant axis with title, normal binning");
    test(RAxisEquidistant("RITLE_E2", 5, -8.9, -6.7),
         "RITLE_E2",
         5,
         -8.9,
         -6.7);
  }
}

TEST(AxisTest, Growable) {
  auto test = [](RAxisGrow& axis, std::string_view title) {
    const RAxisGrow& caxis = axis;

    test_axis_equidistant(caxis, title, true, 10, 1.2, 3.4);

    RAxisConfig cfg(caxis);
    EXPECT_EQ(cfg.GetTitle(), title);
    EXPECT_EQ(cfg.GetNBinsNoOver(), 10);
    EXPECT_EQ(cfg.GetKind(), RAxisConfig::kGrow);
    EXPECT_EQ(cfg.GetBinBorders().size(), 2u);
    EXPECT_EQ(cfg.GetBinBorders()[0], 1.2);
    EXPECT_DOUBLE_EQ(cfg.GetBinBorders()[1], 3.4);
    EXPECT_EQ(cfg.GetBinLabels().size(), 0u);

    // FIXME: Can't test RAxisGrow::Grow() as this method is not implemented.
    //        Once it's implemented, please factor out commonalities with
    //        the RAxisLabels test.
  };

  {
    SCOPED_TRACE("Growable axis w/o title");
    RAxisGrow grow1(10, 1.2, 3.4);
    test(grow1, "");
  }

  {
    SCOPED_TRACE("Growable axis with title");
    RAxisGrow grow2("RITLE_G2", 10, 1.2, 3.4);
    test(grow2, "RITLE_G2");
  }
}

TEST(AxisTest, Irregular) {
  auto test = [](const RAxisIrregular& axis,
                 std::string_view title,
                 const std::vector<double>& bin_borders) {
    const int n_bins_no_over = bin_borders.size() - 1;
    test_axis_base(axis,
                   title,
                   false,
                   n_bins_no_over,
                   bin_borders.front(),
                   bin_borders.back());

    const int underflow_findbin_res = -1;
    const int overflow_findbin_res = -2;
    EXPECT_EQ(axis.FindBin(std::numeric_limits<double>::lowest()), underflow_findbin_res);
    EXPECT_EQ(axis.FindBin(bin_borders.front() - 0.01), underflow_findbin_res);
    for (int bin = 1; bin <= n_bins_no_over; ++bin) {
      const double bin_width = bin_borders[bin] - bin_borders[bin-1];
      EXPECT_EQ(axis.FindBin(bin_borders[bin-1] + 0.01 * bin_width), bin);
      EXPECT_EQ(axis.FindBin(bin_borders[bin] - 0.01 * bin_width), bin);
    }
    EXPECT_EQ(axis.FindBin(bin_borders.back() + 0.01), overflow_findbin_res);
    EXPECT_EQ(axis.FindBin(std::numeric_limits<double>::max()), overflow_findbin_res);

    EXPECT_DOUBLE_EQ(axis.GetBinTo(underflow_findbin_res), bin_borders[0]);
    for (int bin = 1; bin <= n_bins_no_over; ++bin) {
      const double left_border = bin_borders[bin-1];
      const double right_border = bin_borders[bin];
      EXPECT_DOUBLE_EQ(axis.GetBinCenter(bin), (left_border + right_border) / 2.0);
      EXPECT_DOUBLE_EQ(axis.GetBinFrom(bin), left_border);
      EXPECT_DOUBLE_EQ(axis.GetBinTo(bin), right_border);
    }
    EXPECT_DOUBLE_EQ(axis.GetBinFrom(overflow_findbin_res), bin_borders.back());

    EXPECT_EQ(axis.GetBinBorders(), bin_borders);

    const int kInvalidBin = RAxisBase::kInvalidBin;
    EXPECT_EQ(axis.GetBinIndexForLowEdge(std::numeric_limits<double>::lowest()),
              kInvalidBin);
    for (int iborder = 0; iborder < (int) bin_borders.size(); ++iborder) {
      const double border = bin_borders[iborder];
      EXPECT_EQ(axis.GetBinIndexForLowEdge(border - 0.01), kInvalidBin);
      EXPECT_EQ(axis.GetBinIndexForLowEdge(border), iborder + 1);
      EXPECT_EQ(axis.GetBinIndexForLowEdge(border + 0.01), kInvalidBin);
    }
    EXPECT_EQ(axis.GetBinIndexForLowEdge(std::numeric_limits<double>::max()),
              kInvalidBin);

    RAxisConfig cfg(axis);
    EXPECT_EQ(cfg.GetTitle(), title);
    EXPECT_EQ(cfg.GetNBinsNoOver(), n_bins_no_over);
    EXPECT_EQ(cfg.GetKind(), RAxisConfig::kIrregular);
    EXPECT_EQ(cfg.GetBinBorders(), bin_borders);
    EXPECT_EQ(cfg.GetBinLabels().size(), 0u);
  };

  {
    SCOPED_TRACE("Irregular axis w/o title");
    test(RAxisIrregular({-1.2, 3.4, 5.6}), "", {-1.2, 3.4, 5.6});
  }

  {
    SCOPED_TRACE("Irregular axis with title");
    test(RAxisIrregular("RITLE_I2", {2.3, 5.7, 11.13, 17.19}),
         "RITLE_I2",
         {2.3, 5.7, 11.13, 17.19});
  }
}

TEST(AxisTest, Labels) {
  // Test that an RAxisLabels has the expected properties
  //
  // This function was extracted from the test lambda to reduce nesting. It
  // assumes that the input axis was constructed as RAxisLabels(expected_title,
  // labels), and received some number of operations that may have inserted new
  // labels from that point, ultimately leading to the expected_labels set.
  //
  const auto check_labeled_axis =
    [](const RAxisLabels& axis,
       std::string_view expected_title,
       const std::vector<std::string_view>& expected_labels) {
    // Notice that the RAxisBase configuration is _not_ updated when new
    // labels are added. This is by design, according to the RAxisLabels docs.
    // The configuration would be updated on Grow(), but we can't test Grow()
    // right now since it isn't implemented yet...
    test_axis_equidistant(axis,
                          expected_title,
                          true,
                          labels.size(),
                          0.0,
                          static_cast<double>(labels.size()));

    EXPECT_EQ(axis.GetBinLabels().size(), expected_labels.size());
    for (size_t i = 0; i < expected_labels.size(); ++i) {
      EXPECT_EQ(axis.GetBinLabels()[i], expected_labels[i]);
    }

    // Compare the RAxisLabels with various variations of itself
    using BinningCompat = RAxisBase::BinningCompatibility;
    using LabeledCompat = RAxisBase::LabeledBinningCompatibility;
    using CompatFlags = LabeledCompat::Flags;
    auto checkLabeledCompat =
      [&axis](const auto& sourceLabels,
              int expectedCompatFlags,
              auto&& expectedExtraSourceLabels) {
        const RAxisLabels source(sourceLabels);
        const LabeledCompat expected(CompatFlags(expectedCompatFlags),
                                     std::move(expectedExtraSourceLabels));
        EXPECT_EQ(axis.CheckLabeledBinningCompat(source), expected);
        EXPECT_EQ(axis.CheckBinningCompat(source), BinningCompat(expected));
      };
    const int uncommittedTargetLabels =
      (expected_labels.size() - axis.GetNBinsNoOver());
    {
      SCOPED_TRACE("Compare with newly created axis, same labels");
      checkLabeledCompat(
        expected_labels,
        (uncommittedTargetLabels > 0) * CompatFlags::kTargetMustGrow,
        std::vector<std::string_view>{}
      );
    }
    {
      const std::vector<std::string_view> missing_last_label(
        expected_labels.cbegin(), expected_labels.cend() - 1);
      SCOPED_TRACE("Compare with newly created axis, missing last label");
      checkLabeledCompat(
        missing_last_label,
        (uncommittedTargetLabels > 1) * CompatFlags::kTargetMustGrow
        + (uncommittedTargetLabels == 0) * CompatFlags::kExtraTargetBins,
        std::vector<std::string_view>{}
      );
    }
    {
      auto one_extra_label = expected_labels;
      one_extra_label.push_back("I AM ROOT");
      SCOPED_TRACE("Compare with newly created axis, one extra label");
      checkLabeledCompat(
        one_extra_label,
        CompatFlags::kTargetMustGrow,
        std::vector<std::string_view>{"I AM ROOT"}
      );
    }
    auto swapped_labels = expected_labels;
    std::swap(swapped_labels[0], swapped_labels[expected_labels.size()-1]);
    {
      SCOPED_TRACE("Compare with newly created axis, swapped label");
      checkLabeledCompat(
        swapped_labels,
        (uncommittedTargetLabels > 0) * CompatFlags::kTargetMustGrow
        + CompatFlags::kLabelOrderDiffers,
        std::vector<std::string_view>{}
      );
    }
    {
      auto changed_one_label = expected_labels;
      changed_one_label[0] = "I AM ROOT";
      SCOPED_TRACE("Compare with newly created axis, changed one label");
      checkLabeledCompat(
        changed_one_label,
        CompatFlags::kTargetMustGrow
        + CompatFlags::kLabelOrderDiffers
        + CompatFlags::kExtraTargetBins,
        std::vector<std::string_view>{"I AM ROOT"}
      );
    }
    {
      auto missing_first_label = expected_labels;
      missing_first_label.erase(missing_first_label.cbegin());
      SCOPED_TRACE("Compare with newly created axis, missing first label");
      checkLabeledCompat(
        missing_first_label,
        (uncommittedTargetLabels > 0) * CompatFlags::kTargetMustGrow
        + CompatFlags::kLabelOrderDiffers
        + CompatFlags::kExtraTargetBins,
        std::vector<std::string_view>{}
      );
    }
    swapped_labels.push_back("I AM ROOT");
    {
      SCOPED_TRACE("Compare with newly created axis, swapped label + one extra");
      checkLabeledCompat(
        swapped_labels,
        CompatFlags::kTargetMustGrow
        + CompatFlags::kLabelOrderDiffers,
        std::vector<std::string_view>{"I AM ROOT"}
      );
    }

    // RAxisLabels's labeled binning scheme is not considered compatible with
    // any sort of numerical binning
    EXPECT_EQ(axis.CheckBinningCompat(RAxisEquidistant(4, 0.1, 2.3)),
              BinningCompat());
    EXPECT_EQ(RAxisEquidistant(4, 0.1, 2.3).CheckBinningCompat(axis),
              BinningCompat());
    EXPECT_EQ(axis.CheckBinningCompat(RAxisIrregular({1.2, 3.5, 7.9})),
              BinningCompat());
    EXPECT_EQ(RAxisIrregular({1.2, 3.5, 7.9}).CheckBinningCompat(axis),
              BinningCompat());
    EXPECT_EQ(axis.CheckBinningCompat(RAxisGrow(4, 0.1, 2.3)),
              BinningCompat());
    EXPECT_EQ(RAxisGrow(4, 0.1, 2.3).CheckBinningCompat(axis),
              BinningCompat());

    RAxisConfig cfg(axis);
    EXPECT_EQ(cfg.GetTitle(), expected_title);
    EXPECT_EQ(cfg.GetNBinsNoOver(), static_cast<int>(expected_labels.size()));
    EXPECT_EQ(cfg.GetKind(), RAxisConfig::kLabels);
    EXPECT_EQ(cfg.GetBinBorders().size(), 0u);
    EXPECT_EQ(cfg.GetBinLabels().size(), expected_labels.size());
    for (size_t i = 0; i < expected_labels.size(); ++i) {
      EXPECT_EQ(cfg.GetBinLabels()[i], expected_labels[i]);
    }
  };

  auto test = [&](RAxisLabels& axis, std::string_view title) {
    {
      SCOPED_TRACE("Original labels configuration");
      check_labeled_axis(axis, title, labels);
    }

    // Bin queries aren't const in general, but should effectively be when
    // querying bins which already exist.
    EXPECT_EQ(axis.FindBinByName("abc"), 0);
    EXPECT_EQ(axis.FindBinByName("de"), 1);
    EXPECT_EQ(axis.FindBinByName("fghi"), 2);
    EXPECT_EQ(axis.FindBinByName("j"), 3);
    EXPECT_EQ(axis.FindBinByName("klmno"), 4);
    EXPECT_EQ(axis.GetBinCenterByName("abc"), 0.5);
    EXPECT_EQ(axis.GetBinCenterByName("de"), 1.5);
    EXPECT_EQ(axis.GetBinCenterByName("fghi"), 2.5);
    EXPECT_EQ(axis.GetBinCenterByName("j"), 3.5);
    EXPECT_EQ(axis.GetBinCenterByName("klmno"), 4.5);
    {
      SCOPED_TRACE("After querying existing labels");
      check_labeled_axis(axis, title, labels);
    }

    // FIXME: Can't test RAxisGrow::Grow() as this method is not implemented.
    //        Once it's implemented, please factor out commonalities with
    //        the RAxisGrow test.

    // Now let's add some new bins
    auto new_labels = labels;
    EXPECT_EQ(axis.FindBinByName("pq"), 5);
    new_labels.push_back("pq");
    {
      SCOPED_TRACE("After querying a first new label");
      check_labeled_axis(axis, title, new_labels);
    }
    EXPECT_EQ(axis.GetBinCenterByName("pq"), 5.5);
    {
      SCOPED_TRACE("After querying the first new label's center");
      check_labeled_axis(axis, title, new_labels);
    }
    EXPECT_EQ(axis.GetBinCenterByName("rst"), 6.5);
    new_labels.push_back("rst");
    {
      SCOPED_TRACE("After querying a second new label's center");
      check_labeled_axis(axis, title, new_labels);
    }
    EXPECT_EQ(axis.FindBinByName("rst"), 6);
    {
      SCOPED_TRACE("After querying the second new label again");
      check_labeled_axis(axis, title, new_labels);
    }
  };

  {
    SCOPED_TRACE("Labeled axis w/o title");
    RAxisLabels axis(labels);
    test(axis, "");
  }

  {
    SCOPED_TRACE("Labeled axis with title");
    RAxisLabels axis("RITLE_L2", labels);
    test(axis, "RITLE_L2");
  }
}

// NOTE: Labeled binning comparisons are tested in RAxisLabels tests
TEST(AxisTest, NumericBinningCompatibility) {
  using BinningCompat = RAxisBase::BinningCompatibility;
  using NumericCompat = RAxisBase::NumericBinningCompatibility;
  using CompatFlags = NumericCompat::Flags;

  // Check numerix axis binning compatibility with minimal boilerplate
  const auto checkNumericCompat = [](const RAxisBase& target,
                                     const RAxisBase& source,
                                     int expectedCompatFlags) {
      const NumericCompat expected{CompatFlags(expectedCompatFlags)};
      EXPECT_EQ(target.CheckBinningCompat(source), BinningCompat(expected));
    };

  // Make an equidistantly binned axis using various RAxis types
  const auto makeEquidistant = [](int numBins, float min, float max) {
    return RAxisEquidistant(numBins, min, max);
  };
  const auto makeGrowable = [](int numBins, float min, float max) {
    return RAxisGrow(numBins, min, max);
  };
  const auto makeEqBinnedIrregular = [](int numBins, float min, float max) {
    std::vector<double> binBorders(numBins+1);
    for (int i = 0; i <= numBins; ++i) {
      binBorders[i] = min + static_cast<double>(i) / numBins * (max - min);
    }
    return RAxisIrregular(std::move(binBorders));
  };

  // Test scenarios where an equidistantly binned axis is merged into another
  // equidistantly binned axis, which may be growable or not.
  const auto testEqBinnedToEqBinned = [&](const auto& makeTarget,
                                          const auto& makeSource) {
    const auto target = makeTarget(6, 1.2, 4.2);
    const bool fixedSource = !(makeSource(6, 1.2, 4.2).CanGrow());
    const bool needEmptyUnderOver = fixedSource && target.CanGrow();
    {
      SCOPED_TRACE("Source axis has the same binning");
      const auto source = makeSource(6, 1.2, 4.2);
      const bool sameGrowability = source.CanGrow() == target.CanGrow();
      checkNumericCompat(target,
                         source,
                         CompatFlags::kTrivialRegularBinMapping
                         + CompatFlags::kRegularBinBijection
                         + sameGrowability * CompatFlags::kFullBinBijection
                         + needEmptyUnderOver * CompatFlags::kNeedEmptyUnderflow
                         + needEmptyUnderOver * CompatFlags::kNeedEmptyOverflow);
    }
    {
      SCOPED_TRACE("Source axis has one less bin on the left");
      checkNumericCompat(target,
                         makeSource(5, 1.7, 4.2),
                         fixedSource * CompatFlags::kNeedEmptyUnderflow
                         + needEmptyUnderOver * CompatFlags::kNeedEmptyOverflow);
    }
    {
      SCOPED_TRACE("Source axis has one less bin on the right");
      checkNumericCompat(target,
                         makeSource(5, 1.2, 3.7),
                         CompatFlags::kTrivialRegularBinMapping
                         + needEmptyUnderOver * CompatFlags::kNeedEmptyUnderflow
                         + fixedSource * CompatFlags::kNeedEmptyOverflow);
    }
    {
      SCOPED_TRACE("Source axis bins are larger by an integer factor");
      checkNumericCompat(target,
                         makeSource(3, 1.2, 4.2),
                         CompatFlags::kRegularBinAliasing
                         + needEmptyUnderOver * CompatFlags::kNeedEmptyUnderflow
                         + needEmptyUnderOver * CompatFlags::kNeedEmptyOverflow);
    }
    {
      SCOPED_TRACE("Source axis bins are larger by a non-integer factor");
      checkNumericCompat(target,
                         makeSource(4, 1.2, 4.2),
                         CompatFlags::kMergingIsLossy
                         + CompatFlags::kRegularBinAliasing
                         + needEmptyUnderOver * CompatFlags::kNeedEmptyUnderflow
                         + needEmptyUnderOver * CompatFlags::kNeedEmptyOverflow);
    }
    {
      SCOPED_TRACE("Source axis bins are smaller by an integer factor");
      checkNumericCompat(target,
                         makeSource(12, 1.2, 4.2),
                         CompatFlags::kMergingIsLossy
                         + needEmptyUnderOver * CompatFlags::kNeedEmptyUnderflow
                         + needEmptyUnderOver * CompatFlags::kNeedEmptyOverflow);
    }
    {
      SCOPED_TRACE("Source axis bins are smaller by a non-integer factor");
      checkNumericCompat(target,
                         makeSource(11, 1.2, 4.2),
                         CompatFlags::kMergingIsLossy
                         + CompatFlags::kRegularBinAliasing
                         + needEmptyUnderOver * CompatFlags::kNeedEmptyUnderflow
                         + needEmptyUnderOver * CompatFlags::kNeedEmptyOverflow);
    }
  };

  // Test scenarios where an equidistantly binned axis is merged into another
  // non-growable, equidistantly binned axis.
  const auto testEqBinnedToNonGrowable = [&](const auto& makeTarget,
                                             const auto& makeSource) {
    testEqBinnedToEqBinned(makeTarget, makeSource);
    const auto target = makeTarget(6, 1.2, 4.2);
    const bool fixedSource = !(makeSource(6, 1.2, 4.2).CanGrow());
    {
      SCOPED_TRACE("Source axis has one more bin on the left");
      checkNumericCompat(target,
                         makeSource(7, 0.7, 4.2),
                         CompatFlags::kMergingIsLossy);
    }
    {
      SCOPED_TRACE("Source axis has one more bin on the right");
      checkNumericCompat(target,
                         makeSource(7, 1.2, 4.7),
                         CompatFlags::kMergingIsLossy);
    }
    {
      SCOPED_TRACE("Source axis is shifted forward by 0.2 bins");
      checkNumericCompat(target,
                         makeSource(6, 1.3, 4.3),
                         CompatFlags::kTrivialRegularBinMapping
                         + CompatFlags::kRegularBinBijection
                         + fixedSource * CompatFlags::kFullBinBijection
                         + CompatFlags::kMergingIsLossy
                         + CompatFlags::kRegularBinAliasing
                         + fixedSource * CompatFlags::kNeedEmptyUnderflow);
    }
    {
      SCOPED_TRACE("Source axis is shifted forward by 1 bin");
      checkNumericCompat(target,
                         makeSource(6, 1.7, 4.7),
                         CompatFlags::kMergingIsLossy
                         + fixedSource * CompatFlags::kNeedEmptyUnderflow);
    }
    {
      SCOPED_TRACE("Source axis is shifted forward by 1.2 bins");
      checkNumericCompat(target,
                         makeSource(6, 1.8, 4.8),
                         CompatFlags::kMergingIsLossy
                         + CompatFlags::kRegularBinAliasing
                         + fixedSource * CompatFlags::kNeedEmptyUnderflow);
    }
    {
      SCOPED_TRACE("Source axis is shifted backward by 0.2 bins");
      checkNumericCompat(target,
                         makeSource(6, 1.1, 4.1),
                         CompatFlags::kMergingIsLossy
                         + CompatFlags::kRegularBinAliasing
                         + fixedSource * CompatFlags::kNeedEmptyOverflow);
    }
    {
      SCOPED_TRACE("Source axis is shifted backward by 1 bin");
      checkNumericCompat(target,
                         makeSource(6, 0.7, 3.7),
                         CompatFlags::kMergingIsLossy
                         + fixedSource * CompatFlags::kNeedEmptyOverflow);
    }
    {
      SCOPED_TRACE("Source axis is shifted backward by 1.2 bins");
      checkNumericCompat(target,
                         makeSource(6, 0.6, 3.6),
                         CompatFlags::kMergingIsLossy
                         + CompatFlags::kRegularBinAliasing
                         + fixedSource * CompatFlags::kNeedEmptyOverflow);
    }
  };

  // Test scenarios where an irregularly binned axis is merged into an
  // equidistantly binned axis, which may be growable or not.
  const auto testIrregularToEqBinned = [&](const auto& makeTarget) {
    const auto target = makeTarget(6, 1.2, 4.2);
    const bool growableTarget = target.CanGrow();
    const bool fixedTarget = !growableTarget;
    // NOTE: No need to test removing the first or last bin border, that's
    //       morally equivalent to removing a bin which we test separately.
    {
      SCOPED_TRACE("Source axis has an extra inner bin border");
      checkNumericCompat(target,
                         RAxisIrregular({1.2, 1.4, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2}),
                         CompatFlags::kMergingIsLossy
                         + growableTarget * CompatFlags::kNeedEmptyUnderflow
                         + growableTarget * CompatFlags::kNeedEmptyOverflow);
    }
    {
      SCOPED_TRACE("Source axis has one less inner bin border");
      checkNumericCompat(target,
                         RAxisIrregular({1.2, 2.2, 2.7, 3.2, 3.7, 4.2}),
                         CompatFlags::kRegularBinAliasing
                         + growableTarget * CompatFlags::kNeedEmptyUnderflow
                         + growableTarget * CompatFlags::kNeedEmptyOverflow);
    }
    {
      SCOPED_TRACE("First source border is shifted forward");
      checkNumericCompat(target,
                         RAxisIrregular({1.3, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2}),
                         CompatFlags::kTrivialRegularBinMapping
                         + CompatFlags::kRegularBinBijection
                         + fixedTarget * CompatFlags::kFullBinBijection
                         + CompatFlags::kMergingIsLossy
                         + CompatFlags::kNeedEmptyUnderflow
                         + growableTarget * CompatFlags::kNeedEmptyOverflow);
    }
    {
      SCOPED_TRACE("Second source border is shifted forward");
      checkNumericCompat(target,
                         RAxisIrregular({1.2, 1.8, 2.2, 2.7, 3.2, 3.7, 4.2}),
                         CompatFlags::kTrivialRegularBinMapping
                         + CompatFlags::kRegularBinBijection
                         + fixedTarget * CompatFlags::kFullBinBijection
                         + CompatFlags::kMergingIsLossy
                         + CompatFlags::kRegularBinAliasing
                         + growableTarget * CompatFlags::kNeedEmptyUnderflow
                         + growableTarget * CompatFlags::kNeedEmptyOverflow);
    }
    {
      SCOPED_TRACE("Second source border is shifted backward");
      checkNumericCompat(target,
                         RAxisIrregular({1.2, 1.6, 2.2, 2.7, 3.2, 3.7, 4.2}),
                         CompatFlags::kMergingIsLossy
                         + CompatFlags::kRegularBinAliasing
                         + growableTarget * CompatFlags::kNeedEmptyUnderflow
                         + growableTarget * CompatFlags::kNeedEmptyOverflow);
    }
    {
      SCOPED_TRACE("Last source border is shifted backward");
      checkNumericCompat(target,
                         RAxisIrregular({1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.1}),
                         CompatFlags::kTrivialRegularBinMapping
                         + CompatFlags::kRegularBinBijection
                         + fixedTarget * CompatFlags::kFullBinBijection
                         + CompatFlags::kMergingIsLossy
                         + growableTarget * CompatFlags::kNeedEmptyUnderflow
                         + CompatFlags::kNeedEmptyOverflow);
    }
  };

  // Test binning compatibility when merging into an equidistant axis
  {
    SCOPED_TRACE("Target axis is equidistant");
    {
      SCOPED_TRACE("Source axis is equidistant");
      testEqBinnedToNonGrowable(makeEquidistant, makeEquidistant);
    }
    {
      SCOPED_TRACE("Source axis is growable");
      testEqBinnedToNonGrowable(makeEquidistant, makeGrowable);
    }
    {
      SCOPED_TRACE("Source axis is irregular");
      testEqBinnedToNonGrowable(makeEquidistant, makeEqBinnedIrregular);
      testIrregularToEqBinned(makeEquidistant);

      // Outcomes which are specific to the Eq<-Irr scenario
      //
      // If the target axis is not growable, the distance by which we shift
      // extremal bin borders, or the size of the bins that we add on the sides,
      // do not matter. It's all underflow and overflow range anyway. So we do
      // not need to test creating bins on non-standard width on the sides of
      // the source axis, or to test shifting the first/last borders by
      // different amounts, as we'll need to do for growable axes.
      //
      const RAxisEquidistant target(6, 1.2, 4.2);
      {
        SCOPED_TRACE("First source border is shifted backward");
        checkNumericCompat(target,
                           RAxisIrregular({1.1, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2}),
                           CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing);
      }
      {
        SCOPED_TRACE("Last source border is shifted forward");
        checkNumericCompat(target,
                           RAxisIrregular({1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.3}),
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kRegularBinBijection
                           + CompatFlags::kFullBinBijection
                           + CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing);
      }
    }
  }

  // Test binning compatibility when merging into an irregular axis
  {
    SCOPED_TRACE("Target axis is irregular");

    // Deduplicated tests for merging an equidistantly binned source into an
    // irregularly binned target.
    //
    // This is highly symmetrical to the case of merging an irregularly binned
    // source into an equidistantly binned target, but I didn't find a good way
    // to deduplicate it without sacrificing code readability yet...
    //
    const auto testEqBinnedToIrregular = [&](const auto& makeSource) {
      testEqBinnedToEqBinned(makeEqBinnedIrregular, makeSource);
      const auto source = makeSource(6, 1.2, 4.2);
      const bool fixedSource = !source.CanGrow();
      {
        SCOPED_TRACE("Target axis has an extra inner bin border");
        checkNumericCompat(RAxisIrregular({1.2, 1.4, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2}),
                           source,
                           CompatFlags::kRegularBinAliasing);
      }
      {
        SCOPED_TRACE("Target axis has one less inner bin border");
        checkNumericCompat(RAxisIrregular({1.2, 2.2, 2.7, 3.2, 3.7, 4.2}),
                           source,
                           CompatFlags::kMergingIsLossy);
      }
      {
        SCOPED_TRACE("First target border is shifted forward");
        checkNumericCompat(RAxisIrregular({1.3, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2}),
                           source,
                           fixedSource * CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing);
      }
      {
        SCOPED_TRACE("First target border is shifted backward");
        checkNumericCompat(RAxisIrregular({1.1, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2}),
                           source,
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kRegularBinBijection
                           + fixedSource * CompatFlags::kFullBinBijection
                           + CompatFlags::kMergingIsLossy
                           + fixedSource * CompatFlags::kNeedEmptyUnderflow);
      }
      {
        SCOPED_TRACE("Second target border is shifted forward");
        checkNumericCompat(RAxisIrregular({1.2, 1.8, 2.2, 2.7, 3.2, 3.7, 4.2}),
                           source,
                           CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing);
      }
      {
        SCOPED_TRACE("Second target border is shifted backward");
        checkNumericCompat(RAxisIrregular({1.2, 1.6, 2.2, 2.7, 3.2, 3.7, 4.2}),
                           source,
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kRegularBinBijection
                           + fixedSource * CompatFlags::kFullBinBijection
                           + CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing);
      }
      {
        SCOPED_TRACE("Last target border is shifted forward");
        checkNumericCompat(RAxisIrregular({1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.3}),
                           source,
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kRegularBinBijection
                           + fixedSource * CompatFlags::kFullBinBijection
                           + CompatFlags::kMergingIsLossy
                           + fixedSource * CompatFlags::kNeedEmptyOverflow);
      }
      {
        SCOPED_TRACE("Last target border is shifted backward");
        checkNumericCompat(RAxisIrregular({1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.1}),
                           source,
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kRegularBinBijection
                           + fixedSource * CompatFlags::kFullBinBijection
                           + CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing);
      }
    };

    {
      SCOPED_TRACE("Source axis is equidistant");
      testEqBinnedToIrregular(makeEquidistant);
    }
    {
      SCOPED_TRACE("Source axis is growable");
      testEqBinnedToIrregular(makeGrowable);
    }
    {
      SCOPED_TRACE("Source axis is irregular");
      testEqBinnedToIrregular(makeEqBinnedIrregular);
      // NOTE: There are Irr<->Irr specific scenarios, but I did not find one
      //       which is _qualitatively_ different from the Irr<->EqBinned ones.
      //       Please add some here as needed.
    }
  }

  // Test binning compatibility when merging into a growable axis
  {
    SCOPED_TRACE("Target axis is growable");

    // Test scenarios where an equidistantly binned axis is merged into another
    // non-growable, equidistantly binned axis.
    const auto testEqBinnedToGrowable = [&](const auto& makeSource) {
      testEqBinnedToEqBinned(makeGrowable, makeSource);
      const auto target = RAxisGrow(6, 1.2, 4.2);
      const bool growableSource = makeSource(6, 1.2, 4.2).CanGrow();
      const bool fixedSource = !growableSource;
      {
        SCOPED_TRACE("Source axis has one more bin on the left");
        checkNumericCompat(target,
                           makeSource(7, 0.7, 4.2),
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kRegularBinBijection
                           + growableSource * CompatFlags::kFullBinBijection
                           + fixedSource * CompatFlags::kNeedEmptyUnderflow
                           + fixedSource * CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
      {
        SCOPED_TRACE("Source axis has one more bin on the right");
        checkNumericCompat(target,
                           makeSource(7, 1.2, 4.7),
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kRegularBinBijection
                           + growableSource * CompatFlags::kFullBinBijection
                           + fixedSource * CompatFlags::kNeedEmptyUnderflow
                           + fixedSource * CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
      {
        SCOPED_TRACE("Source axis is shifted forward by 0.2 bins");
        checkNumericCompat(target,
                           makeSource(6, 1.3, 4.3),
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing
                           + fixedSource * CompatFlags::kNeedEmptyUnderflow
                           + fixedSource * CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
      {
        SCOPED_TRACE("Source axis is shifted forward by 1 bin");
        checkNumericCompat(target,
                           makeSource(6, 1.7, 4.7),
                           fixedSource * CompatFlags::kNeedEmptyUnderflow
                           + fixedSource * CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
      {
        SCOPED_TRACE("Source axis is shifted forward by 1.2 bins");
        checkNumericCompat(target,
                           makeSource(6, 1.8, 4.8),
                           CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing
                           + fixedSource * CompatFlags::kNeedEmptyUnderflow
                           + fixedSource * CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
      {
        SCOPED_TRACE("Source axis is shifted backward by 0.2 bins");
        checkNumericCompat(target,
                           makeSource(6, 1.1, 4.1),
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing
                           + fixedSource * CompatFlags::kNeedEmptyUnderflow
                           + fixedSource * CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
      {
        SCOPED_TRACE("Source axis is shifted backward by 1 bin");
        checkNumericCompat(target,
                           makeSource(6, 0.7, 3.7),
                           CompatFlags::kTrivialRegularBinMapping
                           + fixedSource * CompatFlags::kNeedEmptyUnderflow
                           + fixedSource * CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
      {
        SCOPED_TRACE("Source axis is shifted backward by 1.2 bins");
        checkNumericCompat(target,
                           makeSource(6, 0.6, 3.6),
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing
                           + fixedSource * CompatFlags::kNeedEmptyUnderflow
                           + fixedSource * CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
    };

    {
      SCOPED_TRACE("Source axis is equidistant");
      testEqBinnedToGrowable(makeEquidistant);
    }
    {
      SCOPED_TRACE("Source axis is growable");
      testEqBinnedToGrowable(makeGrowable);
    }
    {
      SCOPED_TRACE("Source axis is irregular");
      testEqBinnedToGrowable(makeEqBinnedIrregular);
      testIrregularToEqBinned(makeGrowable);

      // Outcomes which are specific to the Grow<-Irr scenario
      const RAxisGrow target(6, 1.2, 4.2);
      {
        SCOPED_TRACE("Creating a left source border at -0.2 bins");
        checkNumericCompat(target,
                           RAxisIrregular({1.1, 1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2}),
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kRegularBinBijection
                           + CompatFlags::kMergingIsLossy
                           + CompatFlags::kNeedEmptyUnderflow
                           + CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
      {
        SCOPED_TRACE("Creating a left source border at -1.2 bin");
        checkNumericCompat(target,
                           RAxisIrregular({0.6, 1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2}),
                           CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing
                           + CompatFlags::kNeedEmptyUnderflow
                           + CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
      {
        SCOPED_TRACE("Creating a right source border at +0.2 bins");
        checkNumericCompat(target,
                           RAxisIrregular({1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2, 4.3}),
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kRegularBinBijection
                           + CompatFlags::kMergingIsLossy
                           + CompatFlags::kNeedEmptyUnderflow
                           + CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
      {
        SCOPED_TRACE("Creating a right source border at +1.2 bin");
        checkNumericCompat(target,
                           RAxisIrregular({1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2, 4.8}),
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing
                           + CompatFlags::kNeedEmptyUnderflow
                           + CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
      {
        SCOPED_TRACE("Shifting the left source border by -0.2 bins");
        checkNumericCompat(target,
                           RAxisIrregular({1.1, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2}),
                           CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing
                           + CompatFlags::kNeedEmptyUnderflow
                           + CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
      {
        SCOPED_TRACE("Shifting the right source border by +0.2 bins");
        checkNumericCompat(target,
                           RAxisIrregular({1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.3}),
                           CompatFlags::kTrivialRegularBinMapping
                           + CompatFlags::kMergingIsLossy
                           + CompatFlags::kRegularBinAliasing
                           + CompatFlags::kNeedEmptyUnderflow
                           + CompatFlags::kNeedEmptyOverflow
                           + CompatFlags::kTargetMustGrow);
      }
    }
  }
}

TEST(AxisTest, ReverseBinLimits) {
  {
    RAxisConfig cfg(10, 3.4, 1.2);
    EXPECT_EQ(cfg.GetBinBorders().size(), 2u);
    EXPECT_DOUBLE_EQ(cfg.GetBinBorders()[0], 1.2);
    EXPECT_DOUBLE_EQ(cfg.GetBinBorders()[1], 3.4);
    EXPECT_EQ(cfg.GetNBinsNoOver(), 10);

    // NOTE: This auto-reversal does _not_ happen when using the explicit
    //       RAxisEquidistant constructor, at the time of writing.
    //
    // RAxisEquidistant axis(10, 3.4, 1.2);
    // EXPECT_DOUBLE_EQ(axis.GetMinimum(), 1.2);
    // EXPECT_DOUBLE_EQ(axis.GetMaximum(), 3.4);
    // EXPECT_DOUBLE_EQ(axis.GetBinFrom(1), 1.2);
    // EXPECT_DOUBLE_EQ(axis.GetBinTo(10), 3.4);
    // EXPECT_EQ(axis.GetNBinsNoOver(), 10);
  }

  {
    RAxisConfig cfg(RAxisConfig::Grow, 10, 3.4, 1.2);
    EXPECT_EQ(cfg.GetBinBorders().size(), 2u);
    EXPECT_DOUBLE_EQ(cfg.GetBinBorders()[0], 1.2);
    EXPECT_DOUBLE_EQ(cfg.GetBinBorders()[1], 3.4);
    EXPECT_EQ(cfg.GetNBinsNoOver(), 10);

    // NOTE: This auto-reversal does _not_ happen when using the explicit
    //       RAxisGrow constructor, at the time of writing.
    //
    // RAxisGrow axis(10, 3.4, 1.2);
    // EXPECT_DOUBLE_EQ(axis.GetMinimum(), 1.2);
    // EXPECT_DOUBLE_EQ(axis.GetMaximum(), 3.4);
    // EXPECT_DOUBLE_EQ(axis.GetBinFrom(0), 1.2);
    // EXPECT_DOUBLE_EQ(axis.GetBinTo(9), 3.4);
    // EXPECT_EQ(axis.GetNBinsNoOver(), 10);
  }
}
