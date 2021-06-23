/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *   IP, Inti Pelupessy,  NL eScience Center, i.pelupessy@esciencecenter.nl  *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#include "gtest/gtest.h"
#include "../test_lib.h"

TEST(MPFEnll, getVal) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  // check whether MPFE produces the same results when using different NumCPU or mode.
  // this defines the baseline against which we compare our MP NLL
  RooRandom::randomGenerator()->SetSeed(3);
  // N.B.: it passes on seeds 1 and 2

  RooWorkspace w;
  w.factory("Gaussian::g(x[-5,5],mu[0,-3,3],sigma[1])");
  auto x = w.var("x");
  RooAbsPdf *pdf = w.pdf("g");
  RooRealVar *mu = w.var("mu");
  RooDataSet *data = pdf->generate(RooArgSet(*x), 10000);
  double results[4];

  RooArgSet values = RooArgSet(*mu, *pdf);

  auto nll1 = pdf->createNLL(*data, RooFit::NumCPU(1));
  results[0] = nll1->getVal();
  delete nll1;
  auto nll2 = pdf->createNLL(*data, RooFit::NumCPU(2));
  results[1] = nll2->getVal();
  delete nll2;
  auto nll3 = pdf->createNLL(*data, RooFit::NumCPU(3));
  results[2] = nll3->getVal();
  delete nll3;
  auto nll4 = pdf->createNLL(*data, RooFit::NumCPU(4));
  results[3] = nll4->getVal();
  delete nll4;
  auto nll1b = pdf->createNLL(*data, RooFit::NumCPU(1));
  auto result1b = nll1b->getVal();
  delete nll1b;
  auto nll2b = pdf->createNLL(*data, RooFit::NumCPU(2));
  auto result2b = nll2b->getVal();
  delete nll2b;

  auto nll1_mpfe = pdf->createNLL(*data, RooFit::NumCPU(-1));
  auto result1_mpfe = nll1_mpfe->getVal();
  delete nll1_mpfe;

  auto nll1_interleave = pdf->createNLL(*data, RooFit::NumCPU(1, 1));
  auto result_interleave1 = nll1_interleave->getVal();
  delete nll1_interleave;
  auto nll2_interleave = pdf->createNLL(*data, RooFit::NumCPU(2, 1));
  auto result_interleave2 = nll2_interleave->getVal();
  delete nll2_interleave;
  auto nll3_interleave = pdf->createNLL(*data, RooFit::NumCPU(3, 1));
  auto result_interleave3 = nll3_interleave->getVal();
  delete nll3_interleave;

  EXPECT_DOUBLE_EQ(Hex(results[0]), Hex(results[1]));
  EXPECT_DOUBLE_EQ(Hex(results[0]), Hex(results[2]));
  EXPECT_DOUBLE_EQ(Hex(results[0]), Hex(results[3]));

  EXPECT_DOUBLE_EQ(Hex(results[0]), Hex(result1b));
  EXPECT_DOUBLE_EQ(Hex(results[1]), Hex(result2b));
  EXPECT_DOUBLE_EQ(Hex(results[0]), Hex(result1_mpfe));

  EXPECT_DOUBLE_EQ(Hex(results[0]), Hex(result_interleave1));
  EXPECT_DOUBLE_EQ(Hex(results[0]), Hex(result_interleave2));
  EXPECT_DOUBLE_EQ(Hex(results[0]), Hex(result_interleave3));
}
