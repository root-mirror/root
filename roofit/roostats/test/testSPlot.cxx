

#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooLinearVar.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"

#include "gtest/gtest.h"

TEST(SPlot, UseWithRooLinearVar) {
  RooRealVar x("x", "observable", 0, 0, 20);
  RooRealVar m("m", "mean", 5., -10, 10);
  RooRealVar s("s", "sigma", 2., 0, 10);
  RooGaussian gaus("gaus", "gaus", x, m, s);

  RooRealVar a("a", "exp", -0.2, -10., 0.);
  RooExponential ex("ex", "ex", x, a);

  RooRealVar common("common", "common scale", 3., 0, 10);
  RooRealVar r1("r1", "ratio 1", 0.3, 0, 10);
  RooRealVar r2("r2", "ratio 2", 0.5, 0, 10);
  RooLinearVar c1("c1", "c1", r1, common, RooFit::RooConst(0.));
  RooLinearVar c2("c2", "c2", r2, common, RooFit::RooConst(0.));
  RooAddPdf sum("sum", "sum", RooArgSet(gaus, ex), RooArgSet(c1, c2));

  auto data = sum.generate(x, 1000);

  RooStats::SPlot splot("splot", "splot", *data, &sum, RooArgSet(c1, c2));
  EXPECT_EQ(splot.GetNumSWeightVars(), 2);
  EXPECT_NE(data->get(0)->find("c1_sw"), nullptr);
}
