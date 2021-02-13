#include "gtest/gtest.h"

// test iterating histogram bins and using the new THistRange and
// THBinIterator classes

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom.h"
#include "THistRange.h"
#include "TStopwatch.h"
#include "TKDTreeBinning.h"
#include "TH2Poly.h"

// global variables
int printLevel = 0;
bool fastMode = true;

double testAxisIter(TH1 *h)
{

   TStopwatch w;
   w.Start();


   int xmin = h->GetXaxis()->GetFirst();
   int xmax = h->GetXaxis()->GetLast();
   int ymin = h->GetYaxis()->GetFirst();
   int ymax = h->GetYaxis()->GetLast();
   int zmin = h->GetZaxis()->GetFirst();
   int zmax = h->GetZaxis()->GetLast();

   double sum = 0;
   for (int i = xmin; i <= xmax; ++i) {
      for (int j = ymin; j <= ymax; ++j) {
         for (int k = zmin; k <= zmax; ++k) {
            sum += h->GetBinContent(i, j, k); // this is faster than GetBin
            // int iglob = h->GetBin(i,j,k);
            // sum += h->GetArray()[iglob];
         }
      }
   }
   w.Stop();
   if (printLevel) {
      w.Print();
      std::cout << "Axis iteration: integral is : " << sum << std::endl;
   }
   return sum;
}

double testBinGlobalIter(TH1 *h)
{

   TStopwatch w;
   w.Start();

   //std::cout << "Use GetBinXYZ and global bin iteration " << std::endl;


   int xmin = h->GetXaxis()->GetFirst();
   int xmax = h->GetXaxis()->GetLast();
   int ymin = h->GetYaxis()->GetFirst();
   int ymax = h->GetYaxis()->GetLast();
   int zmin = h->GetZaxis()->GetFirst();
   int zmax = h->GetZaxis()->GetLast();

   int n = h->GetNcells();

   // look at case when one uses all in range bins
   bool useAllInBins = false;
   bool useTHPoly = (h->IsA() == TH2Poly::Class());// || h->IsA() == TProfile2Poly::Class());
   if (!h->GetXaxis()->TestBit(TAxis::kAxisRange) && !h->GetYaxis()->TestBit(TAxis::kAxisRange) &&
       !h->GetZaxis()->TestBit(TAxis::kAxisRange))
      useAllInBins = true;

   if (h->GetDimension() < 3) {
      zmin = 0;
      zmax = 0;
      if (h->GetDimension() < 2) {
         ymin = 0;
         ymax = 0;
      }
   }

   double sum = 0;
   if (useTHPoly) {
      for (int i = 1; i <= n; ++i) {
         sum += h->GetBinContent(i);
      }
   } else if (!useAllInBins) {
      int ix, iy, iz;
      for (int i = 0; i < n; ++i) {
         // use GetBinXYZ to get axis bin numbers
         h->GetBinXYZ(i, ix, iy, iz);
         if (ix >= xmin && ix <= xmax && iy >= ymin && iy <= ymax && iz >= zmin && iz <= zmax) {
            sum += h->GetBinContent(i);
         }
      }
   } else {
      // exclude only underflow/overflow
      for (int i = 0; i < n; ++i) {
         if (h->IsBinOverflow(i) || h->IsBinUnderflow(i))
            continue;
         sum += h->GetBinContent(i);
      }
   }
   if (printLevel) {
      w.Print();
      std::cout << "GlobalBinIteration: integral is : " << sum << std::endl;
   }
   return sum;
}

double testTHBinIterator(TH1 * h) {

   //std::cout << "Use new TBinIterator class " << std::endl;

   TStopwatch w;
   w.Start();

   THistRange r(h);

   double sum = 0;

   if (printLevel > 1 && h->IsA() == TH2Poly::Class()) {
      for (auto &bin : r) {
         sum += h->GetBinContent(bin);
         std::cout << "iterating on bin " << bin << " content " << h->GetBinContent(bin) << std::endl;
      }
   } else {
      // standard iteration
      for (auto &bin : r) {
         sum += h->GetBinContent(bin);
      }
   }

   if (printLevel) {
      w.Print();
      std::cout << "THBinIterator: integral is : " << sum << std::endl;
   }
   return sum;
}

// creation of histograms
int nEvents = 1000000;
TH1 *createTH1()
{
   int nbins = (fastMode) ? 1000000 : 100000000;
   if (printLevel)
      std::cout << "TH1 tests: create and fill a 1d histogram with " << std::pow(nbins,1) << " nbins" << std::endl;

   auto h1 = new TH1D("h1", "h1", nbins, 0, 10);
   const int n = nEvents;
   for (int i = 0; i < n;  ++i) {
      double x = gRandom->Gaus(5, 2);
      h1->Fill(x);
   }
   return h1;
}
TH1 *createTH2()
{
   int nbins = (fastMode) ? 1000 : 10000;
   if (printLevel)
      std::cout << "TH2 tests: create and fill a 2d histogram with " << std::pow(nbins, 2) << " nbins" << std::endl;
   auto h1 = new TH2D("h2", "h2", nbins, 0, 10, nbins, 0, 10);
   const int n = nEvents;
   for (int i = 0; i < n; ++i) {
      double x = gRandom->Gaus(5, 2);
      double y = gRandom->Gaus(1, 3);
      h1->Fill(x, y);
   }
   return h1;
}
TH1 *createTH3()
{
   int nbins = (fastMode) ? 100 : 500;
   if (printLevel)
      std::cout << "TH3 tests : create and fill 3d histogram with " << std::pow(nbins, 3) << " nbins" << std::endl;
   auto h1 = new TH3D("h3", "h3", nbins, 0, 10, nbins, 0, 10, nbins, 0, 10);
   const int n = nEvents;
   for (int i = 0; i < n; ++i) {
      double x = gRandom->Gaus(5, 2);
      double y = gRandom->Gaus(1, 3);
      double z = gRandom->Gaus(0, 5);
      h1->Fill(x, y, z);
   }
   return h1;
}

// create a TH2Poly using teh TKDTree binning class
TH1 * createTH2Poly() {
   // generate multidim data
   const int n = nEvents;
   std::vector<double> data(2 * n);
   unsigned int nbins = (fastMode) ? 1000 : 10000;

   if (printLevel)
      std::cout << "TH2Poly tests : create and fill TH2Poly histogram with " << nbins << " nbins" << std::endl;

   for (int i = 0; i < n; ++i) {
      double x = gRandom->Gaus(5, 2);
      double y = gRandom->Gaus(1, 3);
      data[i] = x;
      data[n + i] = y;
   }

   TKDTreeBinning* bins = new TKDTreeBinning(n, 2, data, nbins);
   R__ASSERT(bins->GetNBins() == nbins);


   TH2Poly *h2pol = new TH2Poly("h2Poly", "KDTree binning", bins->GetDataMin(0), bins->GetDataMax(0),
                                bins->GetDataMin(1), bins->GetDataMax(1));

   const Double_t *binsMinEdges = bins->GetBinsMinEdges();
   const Double_t *binsMaxEdges = bins->GetBinsMaxEdges();
   for (UInt_t i = 0; i < nbins; ++i) {
      UInt_t edgeDim = i * bins->GetDim();
      h2pol->AddBin(binsMinEdges[edgeDim], binsMinEdges[edgeDim + 1], binsMaxEdges[edgeDim], binsMaxEdges[edgeDim + 1]);
   }

   for (UInt_t i = 1; i <= nbins; ++i)
      h2pol->SetBinContent(i, bins->GetBinContent(i - 1));

   return h2pol;
}
// test classes
class THIterationTestBase : public ::testing::Test {
protected:


   static TH1 *h1;
   static double refValue;
   static void TearDownTestSuite() {
      //GOOD:  this is called after running same test suite and not at
      // end of main program.
      if (printLevel) std::cout << "Delete histogram " << h1->GetName() << std::endl;
      delete h1;
   }
};
TH1 * THIterationTestBase::h1 = nullptr;
double THIterationTestBase::refValue = 0;

// TH1 test classes -
class TestTH1DFull : public THIterationTestBase {
protected:
   static void SetUpTestSuite() {
      h1 = createTH1();
      refValue = h1->Integral();
   }
};
class TestTH1DRange : public THIterationTestBase {
protected:
   static void SetUpTestSuite()
   {
      h1 = createTH1();
      // test 60% of histogram
      int xmin = 0.2 * h1->GetNbinsX();
      int xmax = 0.8 * h1->GetNbinsX();
      h1->GetXaxis()->SetRange(xmin, xmax);
      refValue = h1->Integral(xmin,xmax);
   }
};
// TH2 classes

class TestTH2DFull : public THIterationTestBase {
protected:
   static void SetUpTestSuite()
   {
      h1 = createTH2();
      refValue = h1->Integral();
   }
};

class TestTH2DRange : public THIterationTestBase {
protected:
   static void SetUpTestSuite()
   {
      h1 = createTH2();
      // test 60% of histogram
      int xmin = 0.2 * h1->GetNbinsX();
      int xmax = 0.8 * h1->GetNbinsX();
      h1->GetXaxis()->SetRange(xmin, xmax);
      int ymin = 0.2 * h1->GetNbinsY();
      int ymax = 0.8 * h1->GetNbinsY();
      h1->GetYaxis()->SetRange(ymin, ymax);
      refValue = ((TH2* )h1)->Integral(xmin, xmax, ymin, ymax);
   }
};
// TH3 classes
class TestTH3DFull : public THIterationTestBase {
protected:
   static void SetUpTestSuite()
   {
      h1 = createTH3();
      refValue = h1->Integral();
   }
};

class TestTH3DRange : public THIterationTestBase {
protected:
   static void SetUpTestSuite()
   {
      h1 = createTH3();
      // createTH3();
      // test 60% of histogram
      // test 60% of histogram
      int xmin = 0.2 * h1->GetNbinsX();
      int xmax = 0.8 * h1->GetNbinsX();
      h1->GetXaxis()->SetRange(xmin, xmax);
      int ymin = 0.2 * h1->GetNbinsY();
      int ymax = 0.8 * h1->GetNbinsY();
      h1->GetYaxis()->SetRange(ymin, ymax);
      int zmin = 0.2 * h1->GetNbinsZ();
      int zmax = 0.8 * h1->GetNbinsZ();
      h1->GetZaxis()->SetRange(zmin, zmax);
      refValue = ((TH3*) h1)->Integral(xmin, xmax, ymin, ymax,zmin,zmax);
   }
};
// TH2Pol
class TestTH2PolyFull : public THIterationTestBase {
protected:
   static void SetUpTestSuite()
   {
      h1 = createTH2Poly();
      refValue = h1->Integral();
      if (printLevel)
         std::cout << "Ref TH2Poly integral is " << refValue << std::endl;
   }
};

// test with full range

TEST_F(TestTH1DFull, AxisIteration)
{
   EXPECT_EQ(testAxisIter(h1), refValue);
}
TEST_F(TestTH1DFull, GlobalIteration)
{
   EXPECT_EQ(testBinGlobalIter(h1), refValue);
}
TEST_F(TestTH1DFull, BinIterator)
{
   EXPECT_EQ(testTHBinIterator(h1), refValue);
}

// test with restricted range
TEST_F(TestTH1DRange, AxisIteration)
{
   EXPECT_EQ(testAxisIter(h1), refValue);
}
TEST_F(TestTH1DRange, GlobalIteration)
{
   EXPECT_EQ(testBinGlobalIter(h1), refValue);
}
TEST_F(TestTH1DRange, BinIterator)
{
   EXPECT_EQ(testTHBinIterator(h1), refValue);
}

// TH2D tests
TEST_F(TestTH2DFull, AxisIteration)
{
   EXPECT_EQ(testAxisIter(h1), refValue);
}
TEST_F(TestTH2DFull, GlobalIteration)
{
   EXPECT_EQ(testBinGlobalIter(h1), refValue);
}
TEST_F(TestTH2DFull, BinIterator)
{
   EXPECT_EQ(testTHBinIterator(h1), refValue);
}
// 2D tests with restricted range
TEST_F(TestTH2DRange, AxisIteration)
{
   EXPECT_EQ(testAxisIter(h1), refValue);
}
TEST_F(TestTH2DRange, GlobalIteration)
{
   EXPECT_EQ(testBinGlobalIter(h1), refValue);
}
TEST_F(TestTH2DRange, BinIterator)
{
   EXPECT_EQ(testTHBinIterator(h1), refValue);
}

// TH3D tests
TEST_F(TestTH3DFull, AxisIteration)
{
   EXPECT_EQ(testAxisIter(h1), refValue);
}
TEST_F(TestTH3DFull, GlobalIteration)
{
   EXPECT_EQ(testBinGlobalIter(h1), refValue);
}
TEST_F(TestTH3DFull, BinIterator)
{
   EXPECT_EQ(testTHBinIterator(h1), refValue);
}
// 3D tests with restricted range
TEST_F(TestTH3DRange, AxisIteration)
{
   EXPECT_EQ(testAxisIter(h1), refValue);
}
TEST_F(TestTH3DRange, GlobalIteration)
{
   EXPECT_EQ(testBinGlobalIter(h1), refValue);
}
TEST_F(TestTH3DRange, BinIterator)
{
   EXPECT_EQ(testTHBinIterator(h1), refValue);
}
// th2poly
TEST_F(TestTH2PolyFull, GlobalIteration)
{
   EXPECT_EQ(testBinGlobalIter(h1), refValue);
}
TEST_F(TestTH2PolyFull, BinIterator)
{
   EXPECT_EQ(testTHBinIterator(h1), refValue);
}

int main(int argc, char **argv)
{

   // Disables elapsed time by default.
   //::testing::GTEST_FLAG(print_time) = false;


   // Parse command line arguments
   for (Int_t i = 1; i < argc; i++) {
      std::string arg = argv[i];
      if (arg == "-v") {
         std::cout << "---running in verbose mode" << std::endl;
         printLevel = 1;
      } else if (arg == "-b") {
         std::cout << "---running in benchmark mode to measure time" << std::endl;
         fastMode = false;
      }
   }


   // This allows the user to override the flag on the command line.
   ::testing::InitGoogleTest(&argc, argv);


   return RUN_ALL_TESTS();
}