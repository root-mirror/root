// Tests for the RooTemplateProxy and RooCategory, and if they work together
// Author: Stephan Hageboeck, CERN  01/2019

#include "RooTemplateProxy.h"
#include "RooCategoryProxy.h"

#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooGenericPdf.h"
#include "RooWorkspace.h"

#include "TFile.h"
#include "TMemFile.h"

#include "gtest/gtest.h"


TEST(RooCategory, CategoryDefineMultiState) {
  RooCategory myCat("myCat", "A category", { {"0Lep", 0}, {"1Lep", 1}, {"2Lep", 2}, {"3Lep", 3} });

  std::vector<std::string> targets{"0Lep", "1Lep", "2Lep", "3Lep"};
  int i = 0;
  for (const auto& nameAndIdx : myCat) {
    EXPECT_EQ(nameAndIdx.second, i);
    EXPECT_EQ(nameAndIdx.first, targets[i]);
    i++;
  }
  EXPECT_EQ(myCat.lookupName(1), "1Lep");
  EXPECT_EQ(myCat.lookupIndex("2Lep"), 2);
}


TEST(RooCategory, WriteAndReadBack) {
  RooCategory myCat("myCat", "A category", {
      {"0Lep", 0}, {"1Lep", 1}, {"2Lep", 2}, {"3Lep", 3} });
  myCat.setIndex(2);

  TMemFile file("memfile", "RECREATE");
  file.WriteObject(&myCat, "myCat");

  RooCategory* readCat = nullptr;
  file.GetObject("myCat", readCat);
  ASSERT_NE(readCat, nullptr);
  EXPECT_EQ(readCat->getCurrentIndex(), myCat.getCurrentIndex());
  EXPECT_STREQ(readCat->getLabel(), myCat.getLabel());
  EXPECT_EQ(readCat->size(), 4u);

  auto readIt = std::next(readCat->begin());
  auto origIt = std::next(myCat.begin());
  EXPECT_EQ(readIt->first, origIt->first);
  EXPECT_EQ(readIt->second, origIt->second);
}


class RooCategoryIO : public testing::TestWithParam<const char*> {
public:
  /* Test the reading of a simple mock category that has the states
   * one = 0
   * two = 1
   * three = 2
   * four = 3
   * The ranges "evens" and "odds" for even and odd state names are defined.
   * Now, we check that set ranges are read and written properly, and that
   * sharing of those ranges works even after reading back.
   * A mock file can be created as follows:
      RooCategory cat("cat", "a category")
      cat.defineType("one")
      cat.defineType("two")
      cat.defineType("three")
      cat.defineType("four")
      cat.addToRange("evens", "two,four")
      cat.addToRange("odds", "one,three")
      RooDataSet data("data", "a dataset with a category", RooArgSet(cat))
      data.fill()
      TFile outfile("/tmp/testCategories.root", "RECREATE")
      outfile.WriteObject(&cat, "catOrig")
      outfile.WriteObject(&data, "data")
   */
  void SetUp() override {
    TFile file(GetParam(), "READ");
    ASSERT_TRUE(file.IsOpen());

    file.GetObject("catOrig", cat);
    ASSERT_NE(cat, nullptr);

    file.GetObject("data", data);
    ASSERT_NE(data, nullptr);

    catFromDataset = dynamic_cast<RooCategory*>(data->get(0)->find("cat"));
    ASSERT_NE(catFromDataset, nullptr);
  }

  void TearDown() override {
    delete cat;
    delete data;
  }

protected:
  enum State_t {one = 0, two = 1, three = 2, four = 3};
  RooCategory* cat{nullptr};
  RooDataSet* data{nullptr};
  RooCategory* catFromDataset{nullptr};
};

TEST_P(RooCategoryIO, ReadWithRanges) {
  for (RooCategory* theCat : {cat, catFromDataset}) {
    ASSERT_TRUE(theCat->hasRange("odds"));
    ASSERT_TRUE(theCat->hasRange("evens"));

    EXPECT_TRUE(theCat->isStateInRange("odds",  one));
    EXPECT_TRUE(theCat->isStateInRange("odds",  "three"));
    EXPECT_TRUE(theCat->isStateInRange("evens", two));
    EXPECT_TRUE(theCat->isStateInRange("evens", "four"));

    EXPECT_FALSE(theCat->isStateInRange("odds",  "two"));
    EXPECT_FALSE(theCat->isStateInRange("odds",  four));
    EXPECT_FALSE(theCat->isStateInRange("evens", "one"));
    EXPECT_FALSE(theCat->isStateInRange("evens", three));
  }
}

TEST_P(RooCategoryIO, TestThatRangesAreShared) {
  cat->addToRange("evens", "three");
  EXPECT_TRUE(catFromDataset->isStateInRange("evens", "three"));
}

INSTANTIATE_TEST_SUITE_P(IO_SchemaEvol, RooCategoryIO,
    testing::Values("categories_v620.root", "categories_v621.root", "categories_v622.root"));



TEST(RooCategory, BracketOperator) {
  RooCategory myCat;
  myCat["0Lep"] = 0;
  myCat["1Lep"];
  myCat["Negative"] = -1;
  myCat["2Lep"];

  std::map<int, std::string> targets{{-1,"Negative"}, {0,"0Lep"}, {1,"1Lep"}, {2,"2Lep"}};
  for (const auto& nameAndIndex : myCat) {
    ASSERT_NE(targets.find(nameAndIndex.second), targets.end());
    EXPECT_EQ(nameAndIndex.first, targets[nameAndIndex.second]);
  }
}


TEST(RooCategory, OverwriteActiveState) {
  RooCategory myCat;
  myCat["0Lep"] = 1;
  myCat["1Lep"] = 2;

  EXPECT_EQ(myCat.getCurrentIndex(), 1);

  RooCategory otherCat;
  otherCat["test1"] = 1;
  otherCat["test2"] = 2;

  EXPECT_STREQ(otherCat.getCurrentLabel(), "test1");
}


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
    RooTemplateProxy<RooRealVar> var;
    RooTemplateProxy<RooAbsPdf>  pdf;
};


TEST(RooTemplateProxy, CategoryProxy) {
  RooCategory myCat("myCat", "A category");
  myCat.defineType("A", 1);
  myCat.defineType("B", 2);
  myCat.defineType("NotA", -1);
  std::string longStr(500, '*');
  myCat.defineType(longStr, 500);

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

// Read a simple v6.20 workspace to test proxy schema evolution
TEST(RooProxy, Read6_20) {
  TFile file("testProxiesAndCategories_1.root", "READ");
  ASSERT_TRUE(file.IsOpen());

  RooWorkspace* ws = nullptr;
  file.GetObject("ws", ws);
  ASSERT_NE(ws, nullptr);

  auto pdf = ws->pdf("gaus");
  EXPECT_NE(pdf, nullptr);
  const char* names[3] = {"x", "m", "s"};
  for (int i=0; i<3; ++i) {
    ASSERT_NE(pdf->findServer(names[i]), nullptr);
  }
}
