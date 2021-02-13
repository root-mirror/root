#include <ROOT/TSeq.hxx>
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TSystem.h"

#include "gtest/gtest.h"

TEST(TTreeReaderArray, Vector)
{
   TTree *tree = new TTree("TTreeReaderArrayTree", "In-memory test tree");
   std::vector<float> vecf{17.f, 18.f, 19.f, 20.f, 21.f};
   tree->Branch("vec", &vecf);

   tree->Fill();
   tree->Fill();
   tree->Fill();

   tree->ResetBranchAddresses();

   TTreeReader tr(tree);
   TTreeReaderArray<float> vec(tr, "vec");

   tr.SetEntry(1);

   EXPECT_EQ(5u, vec.GetSize());
   EXPECT_FLOAT_EQ(19.f, vec[2]);
   EXPECT_FLOAT_EQ(17.f, vec[0]);
}

TEST(TTreeReaderArray, MultiReaders)
{
   // See https://root.cern.ch/phpBB3/viewtopic.php?f=3&t=22790
   TTree *tree = new TTree("TTreeReaderArrayTree", "In-memory test tree");
   double Double[6] = {42.f, 43.f, 44.f, 45.f, 46.f, 47.f};
   tree->Branch("D", &Double, "D[4]/D");

   tree->Fill();
   tree->Fill();
   tree->Fill();

   tree->ResetBranchAddresses();

   TTreeReader TR(tree);
   TTreeReaderArray<double> trDouble0(TR, "D");
   TTreeReaderArray<double> trDouble1(TR, "D");
   TTreeReaderArray<double> trDouble2(TR, "D");
   TTreeReaderArray<double> trDouble3(TR, "D");
   TTreeReaderArray<double> trDouble4(TR, "D");
   TTreeReaderArray<double> trDouble5(TR, "D");

   TR.SetEntry(1);

   EXPECT_EQ(4u, trDouble0.GetSize());
   EXPECT_EQ(4u, trDouble1.GetSize());
   EXPECT_EQ(4u, trDouble2.GetSize());
   EXPECT_EQ(4u, trDouble3.GetSize());
   EXPECT_EQ(4u, trDouble4.GetSize());
   EXPECT_EQ(4u, trDouble5.GetSize());
   for (int i = 0; i < 4; ++i) {
      EXPECT_DOUBLE_EQ(Double[i], trDouble0[i]);
      EXPECT_DOUBLE_EQ(Double[i], trDouble1[i]);
      EXPECT_DOUBLE_EQ(Double[i], trDouble2[i]);
      EXPECT_DOUBLE_EQ(Double[i], trDouble3[i]);
      EXPECT_DOUBLE_EQ(Double[i], trDouble4[i]);
      EXPECT_DOUBLE_EQ(Double[i], trDouble5[i]);
   }
}

void checkRV(TTreeReaderArray<bool> &rval, const std::vector<bool> &arr, std::string_view compType)
{
   EXPECT_EQ(rval.GetSize(), arr.size());
   for (auto i : ROOT::TSeqI(rval.GetSize())) {
      const auto bRead = rval[i];
      const auto bExpected = arr[i];
      EXPECT_EQ(bRead, bExpected) << "In the case of " << compType << " element " << i << " of the read collection is "
                                  << std::boolalpha << bRead << " which differs from the expected one, which is "
                                  << bExpected;
   }
}

TEST(TTreeReaderArray, BoolCollections)
{

   // References
   std::vector<bool> va_ref0{false, true, true};
   std::vector<bool> a_ref0{true, false, false, false, true, true, false};
   std::vector<bool> v_ref0{true, false, false, false, true, true};

   std::vector<bool> va_ref1{false, true, false, false, true};
   std::vector<bool> a_ref1{true, true, false, false, false, true, true};
   std::vector<bool> v_ref1{true, false, false, true, true, true, true, false, false, true, true};

   // Tree Setup
   auto fileName = "TTreeReaderArray_BoolCollections.root";
   TFile f(fileName, "recreate");
   TTree t("t", "t");

   int n;
   bool va[5];
   bool a[7];
   std::vector<bool> v;

   t.Branch("v", &v);
   t.Branch("a", a, "a[7]/O");
   t.Branch("n", &n);
   t.Branch("va", va, "va[n]/O");

   // Filling the tree
   n = va_ref0.size();
   std::copy(va_ref0.begin(), va_ref0.end(), va);
   std::copy(a_ref0.begin(), a_ref0.end(), a);
   v = v_ref0;
   t.Fill();

   n = va_ref1.size();
   std::copy(va_ref1.begin(), va_ref1.end(), va);
   std::copy(a_ref1.begin(), a_ref1.end(), a);
   v = v_ref1;
   t.Fill();

   t.Write();
   f.Close();

   // Set up the ttree reader
   TFile f2(fileName);
   TTreeReader r("t", &f2);

   TTreeReaderArray<bool> rv(r, "v");
   TTreeReaderArray<bool> ra(r, "a");
   TTreeReaderArray<bool> rva(r, "va");

   // Loop by hand
   r.Next();
   checkRV(rv, v_ref0, "vector<bool> ev 0");
   checkRV(ra, a_ref0, "fixed size array of bools ev 0");
   checkRV(rva, va_ref0, "variable size array of bools ev 0");

   r.Next();
   checkRV(rv, v_ref1, "vector<bool> ev 1");
   checkRV(ra, a_ref1, "fixed size array of bools ev 1");
   checkRV(rva, va_ref1, "variable size array of bools ev 1");

   gSystem->Unlink(fileName);
}

TEST(TTreeReaderArray, Double32_t)
{
   TTree t("t", "t");

   int n;
   Double32_t arr[64];
   t.Branch("n", &n);
   t.Branch("arr", arr, "arr[n]/d[0,0,10]");
   t.Branch("arr2", arr, "arr2[n]/D");
   std::vector<int> sizes{20, 30};
   float globalIndex = 1.f;
   for (auto ievt : {0, 1}) {
      n = sizes[ievt];
      for (auto inumb = 0; inumb < n; ++inumb) {
         arr[inumb] = 1024 * globalIndex++;
      }
      t.Fill();
   }
   TTreeReader r(&t);
   TTreeReaderArray<double> arrra(r, "arr");
   TTreeReaderArray<double> arr2ra(r, "arr2");
   while (r.Next()) {
      const auto arr_size = arrra.GetSize();
      EXPECT_EQ(arr_size, arr2ra.GetSize()) << "The size of the collections differ!";
      for (auto i = 0U; i < arr_size; ++i) {
         EXPECT_DOUBLE_EQ(arrra[i], arr2ra[i]) << "The content of the element at index " << i
                                               << " in the collections differs!";
      }
   }
}

TEST(TTreeReaderArray, Float16_t)
{
   TTree t("t", "t");

   int n;
   Float16_t arr[64];
   t.Branch("n", &n);
   t.Branch("arr", arr, "arr[n]/f[0,0,10]");
   t.Branch("arr2", arr, "arr2[n]/F");
   std::vector<int> sizes{20, 30};
   float globalIndex = 1.f;
   for (auto ievt : {0, 1}) {
      n = sizes[ievt];
      for (auto inumb = 0; inumb < n; ++inumb) {
         arr[inumb] = 1024 * globalIndex++;
      }
      t.Fill();
   }
   TTreeReader r(&t);
   TTreeReaderArray<float> arrra(r, "arr");
   TTreeReaderArray<float> arr2ra(r, "arr2");
   while (r.Next()) {
      const auto arr_size = arrra.GetSize();
      EXPECT_EQ(arr_size, arr2ra.GetSize()) << "The size of the collections differ!";
      for (auto i = 0U; i < arr_size; ++i) {
         EXPECT_FLOAT_EQ(arrra[i], arr2ra[i]) << "The content of the element at index " << i
                                              << " in the collections differs!";
      }
   }
}

TEST(TTreeReaderArray, ROOT10397)
{
   TTree t("t", "t");
   float x[10];
   int n;
   struct {
      int n = 10;
      float z[10];
   } z;
   t.Branch("n", &n, "n/I");
   t.Branch("x", &x, "y[n]/F");
   t.Branch("z", &z, "n/I:z[n]/F");
   for (int i = 7; i < 10; i++) {
      n = i;
      for (int j = 0; j < n; j++) {
         x[j] = j;
      }
      z.n = 13 - i;
      for (int j = 0; j < 10; ++j)
         z.z[j] = z.n;
      t.Fill();
   };

   TTreeReader r(&t);
   TTreeReaderArray<float> xr(r, "x.y");
   TTreeReaderArray<float> zr(r, "z.z");
   r.Next();
   EXPECT_EQ(xr.GetSize(), 7);
   EXPECT_EQ(zr.GetSize(), 13 - 7);
}

TEST(TTreeReaderArray, LongIntArray)
{
   const auto fname = "TTreeReaderArrayLongIntArray.root";
   {
      TFile f(fname, "recreate");
      long int G[3] = {std::numeric_limits<long int>::min(), 42, std::numeric_limits<long int>::max()};
      int size = 2;
      unsigned long int *g = new unsigned long int[size];
      g[0] = 42;
      g[1] = std::numeric_limits<unsigned long int>::max();
      TTree t("t", "t");
      t.Branch("G", G, "G[3]/G");
      t.Branch("n", &size);
      t.Branch("g", g, "g[n]/g");
      t.Fill();
      t.Write();
   }

   TFile f(fname);
   TTreeReader r("t", &f);
   TTreeReaderArray<long int> rG(r, "G");
   TTreeReaderArray<unsigned long int> rg(r, "g");
   EXPECT_TRUE(r.Next());
   ASSERT_EQ(rG.GetSize(), 3);
   EXPECT_EQ(rG[0], std::numeric_limits<long int>::min());
   EXPECT_EQ(rG[1], 42);
   ASSERT_EQ(rg.GetSize(), 2);
   EXPECT_EQ(rg[0], 42);
   EXPECT_EQ(rg[1], std::numeric_limits<unsigned long int>::max());
   EXPECT_FALSE(r.Next());
}
