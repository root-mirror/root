#include "ROOT/RDataFrame.hxx"
#include "ROOT/TSeq.hxx"

#include "gtest/gtest.h"

using namespace ROOT::RDF;

template <typename COLL>
void CheckBins(const TAxis *axis, const COLL &v)
{
   auto nBins = axis->GetNbins();
   auto nBinsp1 = nBins + 1;
   auto nBinsm1 = nBins - 1;
   EXPECT_EQ(nBinsp1, (int)v.size());
   for (auto i : ROOT::TSeqI(1, nBinsp1)) {
      EXPECT_DOUBLE_EQ(axis->GetBinLowEdge(i), (double)v[i - 1]);
   }
   EXPECT_DOUBLE_EQ(axis->GetBinUpEdge(nBinsm1), (double)v[nBinsm1]);
}

/*
These tests help us checking that the binning of the final histograms/profiles is correct.
We test construction of models with full blown histos/profiles and also using constructors
that mimick the ones of those very classes.
*/

TEST(RDataFrameHistoModels, Histo1D)
{
   ROOT::RDataFrame tdf(10);
   auto x = 0.;
   auto d = tdf.Define("x", [&x]() { return x++; }).Define("w", [&x]() { return x + 1.; });
   auto h1 = d.Histo1D(::TH1D("h1", "h1", 10, 0, 10), "x");
   auto h2 = d.Histo1D({"h2", "h2", 10, 0, 10}, "x");
   auto h1w = d.Histo1D(::TH1D("h0w", "h0w", 10, 0, 10), "x", "w");
   auto h2w = d.Histo1D({"h2w", "h2w", 10, 0, 10}, "x", "w");
   std::vector<float> edgesf{1, 2, 3, 4, 5, 6, 10};
   auto h1edgesf = d.Histo1D(::TH1D("h1edgesf", "h1edgesf", (int)edgesf.size() - 1, edgesf.data()), "x");
   auto h2edgesf = d.Histo1D({"h2edgesf", "h2edgesf", (int)edgesf.size() - 1, edgesf.data()}, "x");
   std::vector<double> edgesd{1, 2, 3, 4, 5, 6, 10};
   auto h1edgesd = d.Histo1D(::TH1D("h1edgesd", "h1edgesd", (int)edgesd.size() - 1, edgesd.data()), "x");
   auto h2edgesd = d.Histo1D({"h2edgesd", "h2edgesd", (int)edgesd.size() - 1, edgesd.data()}, "x");

   TH1DModel m0("m0", "m0", 10, 0, 10);
   TH1DModel m1(::TH1D("m1", "m1", 10, 0, 10));

   auto hm0 = d.Histo1D(m0, "x");
   auto hm1 = d.Histo1D(m1, "x");
   auto hm0w = d.Histo1D(m0, "x", "w");
   auto hm1w = d.Histo1D(m1, "x", "w");

   std::vector<double> ref({0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.});

   CheckBins(h1->GetXaxis(), ref);
   CheckBins(h2->GetXaxis(), ref);
   CheckBins(hm0->GetXaxis(), ref);
   CheckBins(hm1->GetXaxis(), ref);

   CheckBins(h1edgesf->GetXaxis(), edgesf);
   CheckBins(h2edgesf->GetXaxis(), edgesf);
   CheckBins(h1edgesd->GetXaxis(), edgesd);
   CheckBins(h2edgesd->GetXaxis(), edgesd);
}

TEST(RDataFrameHistoModels, Prof1D)
{
   ROOT::RDataFrame tdf(10);
   auto x = 0.;
   auto d = tdf.Define("x", [&x]() { return x++; }).Define("y", [&x]() { return x - .3; }).Define("w", [&x]() {
      return x + 1.;
   });
   auto h1 = d.Profile1D(::TProfile("p1", "p1", 10, 0, 10, 0, 12), "x", "y");
   auto h2 = d.Profile1D({"p1", "p1", 10, 0, 10, 0, 12}, "x", "y");
   auto h1w = d.Profile1D(::TProfile("p0w", "p0w", 10, 0, 10, 0, 12), "x", "y", "w");
   auto h2w = d.Profile1D({"p2w", "p2w", 10, 0, 10, 0, 12}, "x", "y", "w");
   std::vector<float> edgesf{1.1f, 2.f, 3.f, 4.f, 5.f, 6.f, 10.f};
   auto h1edgesf = d.Profile1D(::TProfile("h1edgesf", "h1edgesf", (int)edgesf.size() - 1, edgesf.data()), "x", "y");
   auto h2edgesf = d.Profile1D({"h2edgesf", "h2edgesf", (int)edgesf.size() - 1, edgesf.data()}, "x", "y");
   std::vector<double> edgesd{1.1, 2, 3, 4, 5, 6, 10};
   auto h1edgesd = d.Profile1D(::TProfile("h1edgesd", "h1edgesd", (int)edgesd.size() - 1, edgesd.data()), "x", "y");
   auto h2edgesd = d.Profile1D({"h2edgesd", "h2edgesd", (int)edgesd.size() - 1, edgesd.data()}, "x", "y");

   TProfile1DModel m0("m0", "m0", 10, 0, 10, 0, 12);
   TProfile1DModel m1(::TProfile("m1", "m1", 10, 0, 10, 0, 12));

   auto hm0 = d.Profile1D(m0, "x", "y");
   auto hm1 = d.Profile1D(m1, "x", "y");
   auto hm0w = d.Profile1D(m0, "x", "y", "w");
   auto hm1w = d.Profile1D(m1, "x", "y", "w");

   std::vector<double> ref({0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.});

   CheckBins(h1->GetXaxis(), ref);
   CheckBins(h2->GetXaxis(), ref);
   CheckBins(hm0->GetXaxis(), ref);
   CheckBins(hm1->GetXaxis(), ref);

   CheckBins(h1edgesf->GetXaxis(), edgesf);
   CheckBins(h2edgesf->GetXaxis(), edgesf);
   CheckBins(h1edgesd->GetXaxis(), edgesd);
   CheckBins(h2edgesd->GetXaxis(), edgesd);
}

TEST(RDataFrameHistoModels, Histo2D)
{
   ROOT::RDataFrame tdf(10);
   auto x = 0.;
   auto d = tdf.Define("x", [&x]() { return x++; }).Define("y", [&x]() { return x + .1; });
   auto h1 = d.Histo2D(::TH2D("h1", "h1", 10, 0, 10, 5, 0, 10), "x", "y");
   auto h2 = d.Histo2D({"h2", "h2", 10, 0, 10, 5, 0, 10}, "x", "y");
   std::vector<double> edgesX{1, 2, 3, 4, 5, 6, 10};
   std::vector<double> edgesY{1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 10.1};
   auto h1eX = d.Histo2D(::TH2D("h1eX", "h1eX", (int)edgesX.size() - 1, edgesX.data(), 5, 0, 10), "x", "y");
   auto h2eX = d.Histo2D({"h2eX", "h2eX", (int)edgesX.size() - 1, edgesX.data(), 5, 0, 10}, "x", "y");
   auto h1eY = d.Histo2D(::TH2D("h1eY", "h1eY", 5, 0, 10, (int)edgesY.size() - 1, edgesY.data()), "x", "y");
   auto h2eY = d.Histo2D({"h2eY", "h2eY", 5, 0, 10, (int)edgesY.size() - 1, edgesY.data()}, "x", "y");
   auto h1eXeY = d.Histo2D(
      ::TH2D("h1eXeY", "h1eXeY", (int)edgesX.size() - 1, edgesX.data(), (int)edgesY.size() - 1, edgesY.data()), "x",
      "y");
   auto h2eXeY = d.Histo2D(
      {"h2eXeY", "h2eXeY", (int)edgesX.size() - 1, edgesX.data(), (int)edgesY.size() - 1, edgesY.data()}, "x", "y");
   std::vector<float> edgesXf{1, 2, 3, 4, 5, 6, 10};
   std::vector<float> edgesYf{1.1f, 2.1f, 3.1f, 4.1f, 5.1f, 6.1f, 10.1f};
   auto h1eXeYf = d.Histo2D(
      ::TH2D("h1eXeYf", "h1eXeYf", (int)edgesXf.size() - 1, edgesXf.data(), (int)edgesYf.size() - 1, edgesYf.data()),
      "x", "y");
   auto h2eXeYf = d.Histo2D(
      {"h2eXeYf", "h2eXeYf", (int)edgesXf.size() - 1, edgesXf.data(), (int)edgesYf.size() - 1, edgesYf.data()}, "x",
      "y");

   TH2DModel m0("m0", "m0", 10, 0, 10, 5, 0, 10);
   TH2DModel m1(::TH2D("m1", "m1", 10, 0, 10, 5, 0, 10));

   auto hm0 = d.Histo2D(m0, "x", "y");
   auto hm1 = d.Histo2D(m1, "x", "y");
   auto hm0w = d.Histo2D(m0, "x", "y", "x");
   auto hm1w = d.Histo2D(m1, "x", "y", "x");

   std::vector<double> ref0({0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.});
   std::vector<double> ref1({0., 2., 4., 6., 8., 10.});

   CheckBins(h1->GetXaxis(), ref0);
   CheckBins(h1->GetYaxis(), ref1);
   CheckBins(h2->GetXaxis(), ref0);
   CheckBins(h2->GetYaxis(), ref1);
   CheckBins(h1eX->GetXaxis(), edgesX);
   CheckBins(h1eX->GetYaxis(), ref1);
   CheckBins(h2eX->GetXaxis(), edgesX);
   CheckBins(h2eX->GetYaxis(), ref1);
   CheckBins(h1eY->GetXaxis(), ref1);
   CheckBins(h1eY->GetYaxis(), edgesY);
   CheckBins(h2eY->GetXaxis(), ref1);
   CheckBins(h2eY->GetYaxis(), edgesY);
   CheckBins(hm0->GetXaxis(), ref0);
   CheckBins(hm0->GetYaxis(), ref1);
   CheckBins(hm1->GetXaxis(), ref0);
   CheckBins(hm1->GetYaxis(), ref1);

   CheckBins(h1eXeY->GetXaxis(), edgesX);
   CheckBins(h1eXeY->GetYaxis(), edgesY);
   CheckBins(h2eXeY->GetXaxis(), edgesX);
   CheckBins(h2eXeY->GetYaxis(), edgesY);
   CheckBins(h1eXeYf->GetXaxis(), edgesXf);
   CheckBins(h1eXeYf->GetYaxis(), edgesYf);
   CheckBins(h2eXeYf->GetXaxis(), edgesXf);
   CheckBins(h2eXeYf->GetYaxis(), edgesYf);
}

TEST(RDataFrameHistoModels, Profile2D)
{
   ROOT::RDataFrame tdf(10);
   auto x = 0.;
   auto d = tdf.Define("x", [&x]() { return x++; }).Define("y", [&x]() { return x + .1; }).Define("z", [&x]() {
      return x + 1.2;
   });
   auto h1 = d.Profile2D(::TProfile2D("h1", "h1", 10, 0, 10, 5, 0, 10), "x", "y", "z");
   auto h2 = d.Profile2D({"h2", "h2", 10, 0, 10, 5, 0, 10}, "x", "y", "z");
   std::vector<double> edgesX{1, 2, 3, 4, 5, 6, 10};
   std::vector<double> edgesY{1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 10.1};
   auto h1eX =
      d.Profile2D(::TProfile2D("h1eX", "h1eX", (int)edgesX.size() - 1, edgesX.data(), 5, 0, 10), "x", "y", "z");
   auto h2eX = d.Profile2D({"h2eX", "h2eX", (int)edgesX.size() - 1, edgesX.data(), 5, 0, 10}, "x", "y", "z");
   auto h1eY =
      d.Profile2D(::TProfile2D("h1eY", "h1eY", 5, 0, 10, (int)edgesY.size() - 1, edgesY.data()), "x", "y", "z");
   auto h2eY = d.Profile2D({"h2eY", "h2eY", 5, 0, 10, (int)edgesY.size() - 1, edgesY.data()}, "x", "y", "z");
   auto h1eXeY = d.Profile2D(
      ::TProfile2D("h1eXeY", "h1eXeY", (int)edgesX.size() - 1, edgesX.data(), (int)edgesY.size() - 1, edgesY.data()),
      "x", "y", "z");
   auto h2eXeY =
      d.Profile2D({"h2eXeY", "h2eXeY", (int)edgesX.size() - 1, edgesX.data(), (int)edgesY.size() - 1, edgesY.data()},
                  "x", "y", "z");

   TProfile2DModel m0("m0", "m0", 10, 0, 10, 5, 0, 10);
   TProfile2DModel m1(::TProfile2D("m1", "m1", 10, 0, 10, 5, 0, 10));

   auto hm0 = d.Profile2D(m0, "x", "y", "z");
   auto hm1 = d.Profile2D(m1, "x", "y", "z");
   auto hm0w = d.Profile2D(m0, "x", "y", "z", "x");
   auto hm1w = d.Profile2D(m1, "x", "y", "z", "x");

   std::vector<double> ref0({0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.});
   std::vector<double> ref1({0., 2., 4., 6., 8., 10.});

   CheckBins(h1->GetXaxis(), ref0);
   CheckBins(h1->GetYaxis(), ref1);
   CheckBins(h2->GetXaxis(), ref0);
   CheckBins(h2->GetYaxis(), ref1);
   CheckBins(h1eX->GetXaxis(), edgesX);
   CheckBins(h1eX->GetYaxis(), ref1);
   CheckBins(h2eX->GetXaxis(), edgesX);
   CheckBins(h2eX->GetYaxis(), ref1);
   CheckBins(h1eY->GetXaxis(), ref1);
   CheckBins(h1eY->GetYaxis(), edgesY);
   CheckBins(h2eY->GetXaxis(), ref1);
   CheckBins(h2eY->GetYaxis(), edgesY);
   CheckBins(hm0->GetXaxis(), ref0);
   CheckBins(hm0->GetYaxis(), ref1);
   CheckBins(hm1->GetXaxis(), ref0);
   CheckBins(hm1->GetYaxis(), ref1);

   CheckBins(h1eXeY->GetXaxis(), edgesX);
   CheckBins(h1eXeY->GetYaxis(), edgesY);
   CheckBins(h2eXeY->GetXaxis(), edgesX);
   CheckBins(h2eXeY->GetYaxis(), edgesY);
}

TEST(RDataFrameHistoModels, Histo3D)
{
   ROOT::RDataFrame tdf(10);
   auto x = 0.;
   auto d = tdf.Define("x", [&x]() { return x++; }).Define("y", [&x]() { return x + .1; }).Define("z", [&x]() {
      return x + .1;
   });
   auto h1 = d.Histo3D(::TH3D("h1", "h1", 10, 0, 10, 5, 0, 10, 2, 0, 10), "x", "y", "z");
   auto h2 = d.Histo3D({"h2", "h2", 10, 0, 10, 5, 0, 10, 2, 0, 10}, "x", "y", "z");

   std::vector<double> edgesXd{1, 2, 3, 4, 5, 6, 10};
   std::vector<double> edgesYd{1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 10.1};
   std::vector<double> edgesZd{1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 10.2};
   auto h1e = d.Histo3D(::TH3D("h1e", "h1e", (int)edgesXd.size() - 1, edgesXd.data(), (int)edgesYd.size() - 1,
                               edgesYd.data(), (int)edgesZd.size() - 1, edgesZd.data()),
                        "x", "y", "z");
   auto h2e = d.Histo3D({"h2e", "h2e", (int)edgesXd.size() - 1, edgesXd.data(), (int)edgesYd.size() - 1, edgesYd.data(),
                         (int)edgesZd.size() - 1, edgesZd.data()},
                        "x", "y", "z");

   TH3DModel m0("m0", "m0", 2, 0, 10, 5, 0, 10, 10, 0, 10);
   TH3DModel m1(::TH3D("m1", "m1", 2, 0, 10, 5, 0, 10, 10, 0, 10));

   auto hm0 = d.Histo3D(m0, "x", "y", "z");
   auto hm1 = d.Histo3D(m1, "x", "y", "z");
   auto hm0w = d.Histo3D(m0, "x", "y", "z", "z");
   auto hm1w = d.Histo3D(m1, "x", "y", "z", "z");

   std::vector<double> ref0({0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.});
   std::vector<double> ref1({0., 2., 4., 6., 8., 10.});
   std::vector<double> ref2({0., 5., 10.});

   CheckBins(h1e->GetXaxis(), edgesXd);
   CheckBins(h1e->GetYaxis(), edgesYd);
   CheckBins(h1e->GetZaxis(), edgesZd);
   CheckBins(h2e->GetXaxis(), edgesXd);
   CheckBins(h2e->GetYaxis(), edgesYd);
   CheckBins(h2e->GetZaxis(), edgesZd);
   CheckBins(h1->GetXaxis(), ref0);
   CheckBins(h1->GetYaxis(), ref1);
   CheckBins(h1->GetZaxis(), ref2);
   CheckBins(h2->GetXaxis(), ref0);
   CheckBins(h2->GetYaxis(), ref1);
   CheckBins(h2->GetZaxis(), ref2);
   CheckBins(hm0->GetXaxis(), ref2);
   CheckBins(hm0->GetYaxis(), ref1);
   CheckBins(hm0->GetZaxis(), ref0);
   CheckBins(hm0w->GetXaxis(), ref2);
   CheckBins(hm0w->GetYaxis(), ref1);
   CheckBins(hm0w->GetZaxis(), ref0);
}

TEST(RDataFrameHisto, FillVecBool)
{
    const auto n = 10u;
    ROOT::RDataFrame df(n);
    auto h = df.Define("x", [] { return ROOT::RVec<bool>{true, false}; }).Histo1D<ROOT::RVec<bool>>({"h", "h", 2, -0.5, 1.5}, "x");
    EXPECT_EQ(h->GetBinContent(0), 0u);
    EXPECT_EQ(h->GetBinContent(1), n);
    EXPECT_EQ(h->GetBinContent(2), n);
    EXPECT_EQ(h->GetBinContent(3), 0u);
}
