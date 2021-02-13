// Tests for the RooAbsCollection and derived classes
// Authors: Stephan Hageboeck, CERN  05/2020
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooHelpers.h"

#include "gtest/gtest.h"

/// ROOT-10845 IsOnHeap() always returned false.
TEST(RooArgSet, IsOnHeap) {
  auto setp = new RooArgSet();
  EXPECT_TRUE(setp->IsOnHeap());

  RooArgSet setStack;
  EXPECT_FALSE(setStack.IsOnHeap());
}

TEST(RooArgSet_List, VariadicTemplateConstructor) {
  RooRealVar x("x", "x", 0.);
  RooRealVar y("y", "y", 0.);
  RooRealVar z("z", "z", 0.);

  RooArgSet theSet(x, "Hallo");
  EXPECT_EQ(theSet.size(), 1);
  EXPECT_STREQ(theSet.GetName(), "Hallo");

  RooArgSet theSet2(x, y, "Hallo2");
  EXPECT_EQ(theSet2.size(), 2);
  EXPECT_STREQ(theSet2.GetName(), "Hallo2");

  RooArgSet theSet3(x, y, z, "Hallo3");
  EXPECT_EQ(theSet3.size(), 3);
  EXPECT_STREQ(theSet3.GetName(), "Hallo3");


  RooArgList theList(x, "Hallo");
  EXPECT_EQ(theList.size(), 1);
  EXPECT_STREQ(theList.GetName(), "Hallo");

  RooArgList theList2(x, y, "Hallo2");
  EXPECT_EQ(theList2.size(), 2);
  EXPECT_STREQ(theList2.GetName(), "Hallo2");

  RooArgList theList3(x, y, z, "Hallo3");
  EXPECT_EQ(theList3.size(), 3);
  EXPECT_STREQ(theList3.GetName(), "Hallo3");

}

TEST(RooArgSet, SubscriptOperator) {
  RooRealVar x("x", "x", 0.);
  RooRealVar y("y", "y", 0.);
  RooArgSet theSet(x, y);

  EXPECT_EQ(theSet[0], &x);
  EXPECT_EQ(theSet[1], &y);

  const RooAbsArg& xRef = theSet["x"];
  const RooAbsArg& yRef = theSet["y"];
  EXPECT_EQ(&xRef, static_cast<RooAbsArg*>(&x));
  EXPECT_EQ(&yRef, static_cast<RooAbsArg*>(&y));

  // Silence error to come next:
  RooHelpers::HijackMessageStream hijack(RooFit::ERROR, RooFit::InputArguments);
  EXPECT_THROW(theSet[nullptr], std::invalid_argument);
}

TEST(RooArgSet, FromVector) {
  std::vector<RooRealVar> vars;
  vars.emplace_back("x", "x", 0.);
  vars.emplace_back("y", "y", 0.);

  RooArgSet theSet(vars.begin(), vars.end(), "Hallo");
  EXPECT_EQ(theSet.size(), 2);
  EXPECT_STREQ(theSet.GetName(), "Hallo");
}
