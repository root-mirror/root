// Tests for the RooProxy
// Author: Stephan Hageboeck, CERN  01/2019

#include "RooProxy.h"
#include "RooCategoryProxy.h"

#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooGenericPdf.h"

#include "gtest/gtest.h"

struct DummyClass : public RooAbsPdf {
    DummyClass(RooAbsCategory& theCat, RooRealVar& theVar, RooAbsPdf* thePdf = nullptr) :
      cat("catProxy", "Stores categories", this, theCat),
      var("varProxy", "Stores variables", this, theVar),
      pdf("pdfProxy", "Stores pdfs", this) {
      if (thePdf) {
        pdf.setArg(*thePdf);
      }
    }

    double evaluate() const {
      return 1.;
    }

    void clear() {
      clearValueAndShapeDirty();
    }

    TObject* clone(const char*) const {
      return new TObject();
    }

    RooCategoryProxy cat;
    RooProxy<RooRealVar> var;
    RooProxy<RooAbsPdf>  pdf;
};


TEST(RooProxy, CategoryProxy) {
  RooCategory myCat("myCat", "A category");
  myCat.defineType("A", 1);
  myCat.defineType("B", 2);
  myCat.defineType("NotA", -1);
  std::string longStr(500, '*');
  myCat.defineType(longStr.c_str(), 500);

  RooRealVar x("x", "x", -10, 10);

  DummyClass dummy(myCat, x);
  dummy.clear();

  dummy.cat = 2;
  EXPECT_TRUE(dummy.isValueDirty());
  dummy.clear();
  EXPECT_TRUE(dummy.cat == 2);
  EXPECT_STREQ(dummy.cat.label(), "B");

  dummy.cat = longStr.c_str();
  EXPECT_TRUE(dummy.isValueDirty());
  dummy.clear();
  EXPECT_TRUE(dummy.cat == 500);
  EXPECT_STREQ(dummy.cat.label(), longStr.c_str());

  dummy.cat = std::string("NotA");
  EXPECT_TRUE(dummy.isValueDirty());
  dummy.clear();
  EXPECT_TRUE(dummy.cat == -1);
  EXPECT_STREQ(dummy.cat.label(), "NotA");

  dummy.var = 2.;
  EXPECT_TRUE(dummy.isValueDirty());
  dummy.clear();
  EXPECT_TRUE(dummy.var == 2.);
}


TEST(RooProxy, CategoryProxyBatchAccess) {
  RooCategory myCat("myCat", "A category");
  myCat.defineType("A", 1);
  myCat.defineType("B", 2);
  myCat.defineType("NotA", -1);
  std::string longStr(500, '*');
  myCat.defineType(longStr.c_str(), 500);

  RooRealVar x("x", "x", -10, 10);

  DummyClass dummy(myCat, x);

  RooDataSet data("data", "data", RooArgSet(x, myCat));
  auto& xData = dynamic_cast<RooRealVar&>((*data.get())["x"]);
  auto& catData = dynamic_cast<RooAbsCategoryLValue&>((*data.get())["myCat"]);
  for (unsigned int i=0; i < 9; ++i) {
    xData = i;
    catData = i%2 + 1;
    data.fill();
  }
  xData = 9;
  catData = -1;
  data.fill();

  data.attachBuffers(*dummy.getVariables());

  auto theBatch = dummy.cat->getValBatch(0, 10);
  ASSERT_FALSE(theBatch.empty());
  EXPECT_EQ(theBatch.size(), 10ul);
  EXPECT_EQ(theBatch[0], 1);
  EXPECT_EQ(theBatch[1], 2);
  EXPECT_EQ(theBatch[8], 1);
  EXPECT_EQ(theBatch[9], -1);
}


TEST(RooProxy, RealProxy) {
  RooCategory myCat("myCat", "A category");
  RooRealVar x("x", "x", -10, 10);
  DummyClass dummy(myCat, x);

  RooDataSet data("data", "data", x);
  RooRealVar& xData = dynamic_cast<RooRealVar&>((*data.get())["x"]);
  for (unsigned int i=0; i < 10; ++i) {
    xData = i;
    data.fill();
  }
  data.attachBuffers(dummy);

  auto theBatch = dummy.var->getValBatch(0, 100);
  ASSERT_FALSE(theBatch.empty());
  EXPECT_EQ(theBatch.size(), 10ul);
  EXPECT_EQ(theBatch[2], 2.);
}


TEST(RooProxy, PdfProxy) {
  RooCategory myCat("myCat", "A category");
  RooRealVar x("x", "x", -10, 10);
  RooGenericPdf generic("generic", "generic", "1.+x", x);
  DummyClass dummy(myCat, x, &generic);
  dummy.Print("T");

  RooDataSet data("data", "data", x);
  for (unsigned int i=0; i < 10; ++i) {
    x = i;
    data.fill();
  }
  data.attachBuffers(dummy);

  auto theBatch = dummy.pdf->getValBatch(0, 100);
  ASSERT_FALSE(theBatch.empty());
  EXPECT_EQ(theBatch.size(), 10ul);
  EXPECT_EQ(theBatch[2], 3.);
}
