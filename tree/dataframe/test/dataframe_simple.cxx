/****** Run RDataFrame tests both with and without IMT enabled *******/
#include <gtest/gtest.h>
#include <ROOTUnitTestSupport.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/TSeq.hxx>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TInterpreter.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm> // std::sort
#include <array>
#include <chrono>
#include <thread>
#include <set>
#include <random>

using namespace ROOT;
using namespace ROOT::RDF;
using namespace ROOT::VecOps;

// Fixture for all tests in this file. If parameter is true, run with implicit MT, else run sequentially
class RDFSimpleTests : public ::testing::TestWithParam<bool> {
protected:
   RDFSimpleTests() : NSLOTS(GetParam() ? 4u : 1u)
   {
      if (GetParam())
         ROOT::EnableImplicitMT(NSLOTS);
   }
   ~RDFSimpleTests()
   {
      if (GetParam())
         ROOT::DisableImplicitMT();
   }
   const unsigned int NSLOTS;
};

// Create file `filename` containing a test tree `treeName` with `nevents` events
void FillTree(const char *filename, const char *treeName, int nevents = 0)
{
   TFile f(filename, "RECREATE");
   TTree t(treeName, treeName);
   t.SetAutoFlush(1); // yes, one event per cluster: to make MT more meaningful
   double b1;
   int b2;
   double b3[2];
   unsigned int n;
   int b4[2] = {21, 42};
   t.Branch("b1", &b1);
   t.Branch("b2", &b2);
   t.Branch("b3", b3, "b3[2]/D");
   t.Branch("n", &n);
   t.Branch("b4", b4, "b4[n]/I");
   for (int i = 0; i < nevents; ++i) {
      b1 = i;
      b2 = i * i;
      b3[0] = b1;
      b3[1] = -b1;
      n = i % 2 + 1;
      t.Fill();
   }
   t.Write();
   f.Close();
}

TEST_P(RDFSimpleTests, CreateEmpty)
{
   RDataFrame tdf(10);
   auto c = tdf.Count();
   EXPECT_EQ(10U, *c);
}

TEST_P(RDFSimpleTests, CreateZeroEntries)
{
   RDataFrame tdf(0);
   auto c = tdf.Count();
   EXPECT_EQ(0U, *c);
}

TEST_P(RDFSimpleTests, CreateZeroEntriesWithBranches)
{
   auto filename = "dataframe_simple_0.root";
   auto treename = "t";
   // create input file (at most once per execution of the parametrized gtest)
   static bool hasFile = false;
   if (!hasFile) {
      FillTree(filename, treename);
      hasFile = true;
   }
   RDataFrame tdf(treename, filename);
   auto c = tdf.Count();
   auto m = tdf.Mean("b1");
   EXPECT_EQ(0U, *c);
   EXPECT_EQ(0., *m);
}

TEST_P(RDFSimpleTests, BuildWithTDirectory)
{
   auto filename = "dataframe_simple_1.root";
   auto treename = "t";
   // create input file (at most once per execution of the parametrized gtest)
   static bool hasFile = false;
   if (!hasFile) {
      FillTree(filename, treename, 50);
      hasFile = true;
   }
   TFile f(filename);
   RDataFrame tdf(treename, &f);
   auto c = tdf.Count();
   EXPECT_EQ(50U, *c);
}

// Jitting of column types
TEST_P(RDFSimpleTests, TypeGuessing)
{
   auto filename = "dataframe_simple_2.root";
   auto treename = "t";
   // create input file (at most once per execution of the parametrized gtest)
   static bool hasFile = false;
   if (!hasFile) {
      FillTree(filename, treename, 50);
      hasFile = true;
   }
   RDataFrame tdf(treename, filename, {"b1"});
   auto hcompiled = tdf.Histo1D<double>();
   auto hjitted = tdf.Histo1D();
   EXPECT_EQ(50, hcompiled->GetEntries());
   EXPECT_EQ(50, hjitted->GetEntries());
   EXPECT_DOUBLE_EQ(hcompiled->GetMean(), hjitted->GetMean());
}

// Define

TEST_P(RDFSimpleTests, Define_lambda)
{
   RDataFrame tdf(10);
   auto d = tdf.Define("i", []() { return 1; });
   auto m = d.Mean("i");
   EXPECT_DOUBLE_EQ(1., *m);
}

int DefineFunction()
{
   return 1;
}

TEST_P(RDFSimpleTests, Define_function)
{
   RDataFrame tdf(10);
   auto d = tdf.Define("i", DefineFunction);
   auto m = d.Mean("i");
   EXPECT_DOUBLE_EQ(1., *m);
}

struct DefineStruct {
   int operator()() { return 1; }
};

TEST_P(RDFSimpleTests, Define_functor)
{
   RDataFrame tdf(10);
   DefineStruct def;
   auto d = tdf.Define("i", def);
   auto m = d.Mean("i");
   EXPECT_DOUBLE_EQ(1., *m);
}

TEST_P(RDFSimpleTests, Define_jitted)
{
   RDataFrame tdf(10);
   auto d = tdf.Define("i", "1");
   auto m = d.Mean("i");
   EXPECT_DOUBLE_EQ(1., *m);
}

struct RFoo {};

TEST_P(RDFSimpleTests, Define_jitted_type_unknown_to_interpreter)
{
   RDataFrame tdf(10);
   auto d = tdf.Define("foo", [](){return RFoo();});
   auto d2 = tdf.Define("foo2", [](){return std::array<RFoo, 2>();});

   // We check that if nothing is done with RFoo in jitted strings everything works fine
   EXPECT_EQ(10U, *d.Count());

   EXPECT_ANY_THROW(d.Define("foo3", "foo*2"));
}

TEST_P(RDFSimpleTests, Define_jitted_complex)
{
   // this test case (as all others) is usually run twice, in IMT and non-IMT mode,
   // but we only want to create the TRandom object once.
   static bool hasJittedTRandom = false;
   if (!hasJittedTRandom) {
      gInterpreter->ProcessLine("TRandom r;");
      hasJittedTRandom = true;
   }
   gInterpreter->ProcessLine("r.SetSeed(1);");
   RDataFrame tdf(50);
   auto d = tdf.Define("i", "r.Uniform(0.,8.)");
   auto m = d.Max("i");
   EXPECT_EQ(7.867497533559811628, *m);
}

TEST_P(RDFSimpleTests, Define_jitted_complex_array_sum)
{
   RDataFrame tdf(10);
   auto d = tdf.Define("x", "3.0")
               .Define("y", "4.0")
               .Define("z", "12.0")
               .Define("v", "std::array<double, 3> v{x, y, z}; return v;")
               .Define("r", "double r2 = 0.0; for (auto&& w : v) r2 += w*w; return sqrt(r2);");
   auto m = d.Max("r");
   EXPECT_DOUBLE_EQ(13.0, *m);
}

TEST_P(RDFSimpleTests, Define_jitted_defines_with_return)
{
   RDataFrame tdf(10);
   auto d = tdf.Define("my_return_x", "return 3.0")
               .Define("return_y", "4.0 // with a comment")
               .Define("v", "std::array<double, 2> v{my_return_x, return_y}; return v; // also with comment")
               .Define("r", "double r2 = 0.0; for (auto&& w : v) r2 += w*w; return sqrt(r2);");
   auto m = d.Max("r");
   EXPECT_DOUBLE_EQ(5.0, *m);
}

// Define + Filters
TEST_P(RDFSimpleTests, Define_Filter)
{
   TRandom r(1);
   RDataFrame tdf(50);
   auto d = tdf.Define("r", [&r]() { return r.Uniform(0., 8.); });
   auto df = d.Filter([](double x) { return x > 5; }, {"r"});
   auto m = df.Max("r");
   EXPECT_DOUBLE_EQ(7.867497533559811628, *m);
}

TEST_P(RDFSimpleTests, Define_Filter_jitted)
{
   TRandom r(1);
   RDataFrame tdf(50);
   auto d = tdf.Define("r", [&r]() { return r.Uniform(0., 8.); });
   auto df = d.Filter("r>5");
   auto m = df.Max("r");
   EXPECT_EQ(7.867497533559811628, *m);
}

TEST_P(RDFSimpleTests, Define_Filter_named)
{
   TRandom r(1);
   RDataFrame tdf(50);
   auto d = tdf.Define("r", [&r]() { return r.Uniform(0., 8.); });
   auto df = d.Filter([](double x) { return x > 5; }, {"r"}, "myFilter");
   auto m = df.Max("r");
   EXPECT_EQ(7.867497533559811628, *m);
}

TEST_P(RDFSimpleTests, Define_Filter_named_jitted)
{
   TRandom r(1);
   RDataFrame tdf(50);
   auto d = tdf.Define("r", [&r]() { return r.Uniform(0., 8.); });
   auto df = d.Filter("r>5", "myFilter");
   auto m = df.Max("r");
   EXPECT_EQ(7.867497533559811628, *m);
}

// jitted Define + Filters
TEST_P(RDFSimpleTests, Define_jitted_Filter)
{
   gInterpreter->ProcessLine("r.SetSeed(1);");
   RDataFrame tdf(50);
   auto d = tdf.Define("r", "r.Uniform(0.,8.)");
   auto df = d.Filter([](double x) { return x > 5; }, {"r"});
   auto m = df.Max("r");
   EXPECT_EQ(7.867497533559811628, *m);
}

TEST_P(RDFSimpleTests, Define_jitted_Filter_jitted)
{
   gInterpreter->ProcessLine("r.SetSeed(1);");
   RDataFrame tdf(50);
   auto d = tdf.Define("r", "r.Uniform(0.,8.)");
   auto df = d.Filter("r>5");
   auto m = df.Max("r");
   EXPECT_EQ(7.867497533559811628, *m);
}

TEST_P(RDFSimpleTests, Define_jitted_Filter_named)
{
   gInterpreter->ProcessLine("r.SetSeed(1);");
   RDataFrame tdf(50);
   auto d = tdf.Define("r", "r.Uniform(0.,8.)");
   auto df = d.Filter([](double x) { return x > 5; }, {"r"}, "myFilter");
   auto m = df.Max("r");
   EXPECT_EQ(7.867497533559811628, *m);
}

TEST_P(RDFSimpleTests, Define_jitted_Filter_named_jitted)
{
   gInterpreter->ProcessLine("r.SetSeed(1);");
   RDataFrame tdf(50);
   auto d = tdf.Define("r", "r.Uniform(0.,8.)");
   auto df = d.Filter("r>5", "myFilter");
   auto m = df.Max("r");
   EXPECT_EQ(7.867497533559811628, *m);
}

TEST_P(RDFSimpleTests, Define_jitted_Filter_complex_array)
{
   gInterpreter->ProcessLine("r.SetSeed(1);");
   RDataFrame tdf(50);
   auto d = tdf.Define("x", "r.Uniform(0.0, 1.0)")
               .Define("y", "r.Uniform(0.0, 1.0)")
               .Define("z", "r.Uniform(0.0, 1.0)")
               .Define("v", "std::array<double, 3> v{x, y, z}; return v;")
               .Define("r", "double r2 = 0.0; for (auto&& w : v) r2 += w*w; return sqrt(r2);");
   auto dfin = d.Filter("r <= 1.0", "inside");
   auto dfout = d.Filter("bool out = r > 1.0; return out;", "outside");
   auto in = dfin.Count();
   auto out = dfout.Count();

   EXPECT_TRUE(*in < 50U);
   EXPECT_TRUE(*out < 50U);
   EXPECT_EQ(50U, *in + *out);
}

TEST_P(RDFSimpleTests, DefineSlotConsistency)
{
   RDataFrame df(8);
   auto m = df.DefineSlot("x", [](unsigned int) { return 1.; }).Max("x");
   EXPECT_EQ(1., *m);
}

TEST_P(RDFSimpleTests, DefineSlot)
{
   std::vector<int> values(NSLOTS);
   for (auto i = 0u; i < NSLOTS; ++i)
      values[i] = i;
   RDataFrame df(NSLOTS);
   auto ddf = df.DefineSlot("s", [values](unsigned int slot) { return values[slot]; });
   auto m = ddf.Max("s");
   EXPECT_EQ(*m, NSLOTS - 1); // no matter the order of processing, the higher slot number is always taken at least once
}

TEST_P(RDFSimpleTests, DefineSlotCheckMT)
{
   std::vector<unsigned int> ids(NSLOTS, 0u);
   RDataFrame d(NSLOTS);
   auto m = d.DefineSlot("x", [&](unsigned int slot) {
                ids[slot] = 1u;
                return 1;
             }).Max("x");
   EXPECT_EQ(1, *m); // just in case

   const auto nUsedSlots = std::accumulate(ids.begin(), ids.end(), 0u);
   EXPECT_GT(nUsedSlots, 0u);
   EXPECT_LE(nUsedSlots, NSLOTS);
}

TEST_P(RDFSimpleTests, DefineSlotEntry)
{
   const auto nEntries = 8u;
   RDataFrame df(nEntries);
   auto es = df.DefineSlotEntry("e", [](unsigned int, ULong64_t e) { return e; }).Take<ULong64_t>("e");
   auto entries = *es;
   std::sort(entries.begin(), entries.end());
   for (auto i = 0u; i < nEntries; ++i) {
      EXPECT_EQ(i, entries[i]);
   }
}

TEST_P(RDFSimpleTests, Define_Multiple)
{
   RDataFrame tdf(3);
   auto root = tdf.Define("root", "0");
   auto branch1 = root.Define("b", []() { return 1; });
   auto branch2 = root.Define("b", "2");
   auto branch3 = root.Define("b", "4.2"); // check a second jitted branch (checks for ROOT-9754)

   auto rootMax1 = branch1.Max("root");
   auto rootMax2 = branch2.Max<int>("root");
   auto branch1Max = branch1.Max("b");
   auto branch2Max = branch2.Max<int>("b");
   auto branch3Max = branch3.Max<double>("b");

   // Checking that both branches see the same root column
   EXPECT_EQ(*rootMax1, *rootMax2);
   EXPECT_EQ(*rootMax1, 0);

   // Name collision must not represent a problem
   EXPECT_EQ(*branch1Max, 1);
   EXPECT_EQ(*branch2Max, 2);
   EXPECT_EQ(*branch3Max, 4.2);
}

TEST_P(RDFSimpleTests, Define_Multiple_Filter)
{
   RDataFrame tdf(3);
   auto notFilteringFilter = [](int b1) { return b1; }; // also tests returning a type convertible to bool
   auto filteringFilter = [](int b2) { return b2 < 1; };

   auto root = tdf.Define("root", "0");
   auto branch1 = root.Define("b", []() { return 1; }).Filter(notFilteringFilter, {"b"});
   auto branch2 = root.Define("b", "2").Filter(filteringFilter, {"b"});

   auto branch1Count = branch1.Count();
   auto branch2Count = branch2.Count();

   EXPECT_EQ(*branch1Count, 3U);
   EXPECT_EQ(*branch2Count, 0U);
}

TEST_P(RDFSimpleTests, GetNSlots)
{
   EXPECT_EQ(NSLOTS, ROOT::Internal::RDF::GetNSlots());
}

TEST_P(RDFSimpleTests, GetNRuns)
{
   RDataFrame df(3);
   EXPECT_EQ(df.GetNRuns(), 0u);

   auto sum1 = df.Sum("rdfentry_");
   sum1.GetValue();
   EXPECT_EQ(df.GetNRuns(), 1u);

   auto sum2 = df.Sum("rdfentry_");
   sum2.GetValue();
   EXPECT_EQ(df.GetNRuns(), 2u);
}

TEST_P(RDFSimpleTests, CArraysFromTree)
{
   auto filename = "dataframe_simple_3.root";
   auto treename = "t";
   // create input file (at most once per execution of the parametrized gtest)
   static bool hasFile = false;
   if (!hasFile) {
      FillTree(filename, treename, 10);
      hasFile = true;
   }
   RDataFrame df(treename, filename);

   // no jitting
   auto h = df.Filter([](double b1, unsigned int n, RVec<double> &b3,
                         RVec<int> &b4) { return b3[0] == b1 && b4[0] == 21 && b4.size() == n; },
                      {"b1", "n", "b3", "b4"})
               .Histo1D<RVec<double>>("b3");
   EXPECT_EQ(20, h->GetEntries());

   // jitting
   auto h_jit = df.Filter(/*"b3[0] == b1"*/ "b4[0] == 21" /*"b4.size() == n"*/).Histo1D("b3");
   EXPECT_EQ(20, h_jit->GetEntries());
}

TEST_P(RDFSimpleTests, TakeCarrays)
{
   auto treeName = "t";
   auto fileName = "TakeCarrays.root";

   {
      TFile f(fileName, "RECREATE");
      TTree t(treeName, treeName);
      float arr[4];
      t.Branch("arr", arr, "arr[4]/F");
      for (auto i : ROOT::TSeqU(4)) {
         for (auto j : ROOT::TSeqU(4)) {
            arr[j] = i + j;
         }
         t.Fill();
      }
      t.Write();
   }

   RDataFrame tdf(treeName, fileName);
   // no auto here: we check that the type is a COLL<vector<float>>!
   using ColType_t = ROOT::VecOps::RVec<float>;
   std::vector<ColType_t> v = *tdf.Take<ColType_t>("arr");
   std::deque<ColType_t> d = *tdf.Take<ColType_t, std::deque<ColType_t>>("arr");
   std::list<ColType_t> l = *tdf.Take<ColType_t, std::list<ColType_t>>("arr");

   auto lit = l.begin();
   auto ifloat = 0.f;
   for (auto i : ROOT::TSeqU(4)) {
      const auto &vv = v[i];
      const auto &dv = d[i];
      const auto &lv = *lit;
      for (auto j : ROOT::TSeqU(4)) {
         const auto ref = ifloat + j;
         EXPECT_EQ(ref, vv[j]);
         EXPECT_EQ(ref, dv[j]);
         EXPECT_EQ(ref, lv[j]);
      }
      ifloat++;
      lit++;
   }

   // Now we check that the tvecs are not adopting
   EXPECT_TRUE(v[0].data() != v[1].data());
   EXPECT_TRUE(v[1].data() != v[2].data());
   EXPECT_TRUE(v[2].data() != v[3].data());

   gSystem->Unlink(fileName);
}

TEST_P(RDFSimpleTests, Reduce)
{
   auto d = RDataFrame(5).DefineSlotEntry("x", [](unsigned int, ULong64_t e) { return static_cast<int>(e) + 1; });
   auto r1 = d.Reduce([](int x, int y) { return x + y; }, "x");
   auto r2 = d.Reduce([](int x, int y) { return x * y; }, "x", 1);
   auto r3 = d.Define("b", [] { return true; }).Reduce([](bool a, bool b) { return a == b; }, "b", true);
   EXPECT_EQ(*r1, 15);
   EXPECT_EQ(*r2, 120);
   EXPECT_EQ(*r3, true);
}

TEST_P(RDFSimpleTests, Aggregate)
{
   auto d = RDataFrame(5).DefineSlotEntry("x", [](unsigned int, ULong64_t e) { return static_cast<int>(e) + 1; });
   // acc U(U,T), merge U(U,U), default initValue
   auto r1 = d.Aggregate([](int x, int y) { return x + y; }, [](int x, int y) { return x + y; }, "x");
   // acc U(U,T), merge U(U,U), initValue
   auto r2 = d.Aggregate([](int x, int y) { return x * y; }, [](int x, int y) { return x * y; }, "x", 1);
   EXPECT_EQ(*r1, 15);
   EXPECT_EQ(*r2, 120);
}



TEST_P(RDFSimpleTests, AggregateGraph)
{
   auto d = RDataFrame(20).DefineSlotEntry("x", [](unsigned int, ULong64_t e) { return static_cast<double>(e); });
   auto graph = d.Aggregate([](TGraph &g, double x) { g.SetPoint(g.GetN(), x, x * x); },
                            [](std::vector<TGraph> &graphs) {
                               TList l;
                               for (auto g = graphs.begin() + 1; g != graphs.end(); ++g)
                                  l.Add(&(*g));
                               graphs[0].Merge(&l);
                            },
                            "x");
   EXPECT_EQ(graph->GetN(), 20);
   // collect data-points, sorted by x values (they can be scrambled in MT executions)
   using Point_t = std::pair<double, double>;
   std::vector<Point_t> points;
   Point_t p;
   for (int i = 0; i < 20; ++i) {
      graph->GetPoint(i, p.first, p.second);
      points.emplace_back(p);
   }
   // check that all data-points are retrieved correctly
   std::sort(points.begin(), points.end(), [](const Point_t &p1, const Point_t &p2) { return p1.first < p2.first; });
   for (int i = 0; i < 20; ++i) {
      EXPECT_DOUBLE_EQ(points[i].first, i);
      EXPECT_DOUBLE_EQ(points[i].second, i * i);
   }
}

TEST_P(RDFSimpleTests, Graph)
{
   static const int NR_ELEMENTS = 20;

   // Define the source for the graph
   std::vector<int> source(NR_ELEMENTS);
   for (int i = 0; i < NR_ELEMENTS; ++i)
      source[i] = i;

   // Create the graph from the Dataframe
   ROOT::RDataFrame d(NR_ELEMENTS);
   auto dd = d.DefineSlotEntry("x1",
                               [&source](unsigned int slot, ULong64_t entry) {
                                  (void)slot;
                                  return source[entry];
                               })
                .DefineSlotEntry("x2", [&source](unsigned int slot, ULong64_t entry) {
                   (void)slot;
                   return source[entry];
                });

   auto dfGraph = dd.Graph("x1", "x2");
   EXPECT_EQ(dfGraph->GetN(), NR_ELEMENTS);
   EXPECT_STREQ(dfGraph->GetXaxis()->GetTitle(), "x1");
   EXPECT_STREQ(dfGraph->GetYaxis()->GetTitle(), "x2");

   //To perform the test, it's easier to sort
   dfGraph->Sort();

   Double_t x, y;
   for (int i = 0; i < NR_ELEMENTS; ++i) {
      dfGraph->GetPoint(i, x, y);
      EXPECT_EQ(i, x);
      EXPECT_EQ(i, y);
   }
}

class MaxSlotHelper : public ROOT::Detail::RDF::RActionImpl<MaxSlotHelper> {
   const std::shared_ptr<unsigned int> fMaxSlot; // final result
   std::vector<unsigned int> fMaxSlots;          // per-thread partial results
public:
   MaxSlotHelper(unsigned int nSlots)
      : fMaxSlot(std::make_shared<unsigned int>(std::numeric_limits<unsigned int>::lowest())),
        fMaxSlots(nSlots, std::numeric_limits<unsigned int>::lowest())
   {
   }
   MaxSlotHelper(MaxSlotHelper &&) = default;
   MaxSlotHelper(const MaxSlotHelper &) = delete;
   using Result_t = unsigned int;
   std::shared_ptr<unsigned int> GetResultPtr() const { return fMaxSlot; }
   void Initialize() {}
   void InitTask(TTreeReader *, unsigned int) {}
   void Exec(unsigned int slot, unsigned int /*slot2*/) { fMaxSlots[slot] = std::max(fMaxSlots[slot], slot); }
   void Finalize() { *fMaxSlot = *std::max_element(fMaxSlots.begin(), fMaxSlots.end()); }

   std::string GetActionName() { return "MaxSlot"; }
};

TEST_P(RDFSimpleTests, BookCustomAction)
{
   RDataFrame d(1);
   const auto nWorkers = std::max(1u, ROOT::GetThreadPoolSize());
   const auto expected = nWorkers-1;

   auto maxSlot0 = d.Book<unsigned int>(MaxSlotHelper(nWorkers), {"tdfslot_"});
   auto maxSlot1 = d.Book<unsigned int>(MaxSlotHelper(nWorkers), {"rdfslot_"});
   EXPECT_EQ(*maxSlot0, expected);
   EXPECT_EQ(*maxSlot1, expected);
}

class StdDevTestHelper {
private:
   std::default_random_engine fGenerator;
   std::normal_distribution<double> fDistribution;
   std::vector<double> samples;

public:
   void GenerateNumbers(int n)
   {
      std::vector<double> numbers;
      for (int i = 0; i < n; ++i)
         numbers.push_back(fDistribution(fGenerator));
      samples = numbers;
   }

   double stdDevFromDefinition()
   {
      // Calculating the Variance using the definition
      int nSamples = samples.size();
      double mean = 0;
      for (int i = 0; i < nSamples; ++i) {
         mean += samples[i];
      }
      mean = mean / nSamples;

      double varianceRight = 0;

      for (int i = 0; i < nSamples; ++i) {
         varianceRight += std::pow((samples[i] - mean), 2);
      }
      varianceRight = varianceRight / (nSamples - 1);
      return std::sqrt(varianceRight);
   }

   double stdDevFromWelford()
   {
      ROOT::RDataFrame d(samples.size());
      return *d.DefineSlotEntry("x",
                                [this](unsigned int slot, ULong64_t entry) {
                                   (void)slot;
                                   return samples[entry];
                                })
                 .StdDev("x");
   }
};

TEST_P(RDFSimpleTests, StandardDeviation)
{
   const auto expected = 2.4494897427831779;
   RDataFrame rd(8);
   auto stdDev0 = rd.StdDev<ULong64_t>("tdfentry_");
   auto stdDev1 = rd.StdDev<ULong64_t>("rdfentry_");
   EXPECT_DOUBLE_EQ(*stdDev0, expected);
   EXPECT_DOUBLE_EQ(*stdDev1, expected);
}

TEST_P(RDFSimpleTests, StandardDeviationPrecision)
{
   const int maxNSamples = 100;
   const int step = 10;
   const int nTrials = 1;

   StdDevTestHelper helper;

   for (int j = 2; j < maxNSamples; j += step) {
      for (int i = 0; i < nTrials; ++i) {
         auto varianceFromDef = helper.stdDevFromDefinition();
         auto varianceFromWel = helper.stdDevFromWelford();
         EXPECT_DOUBLE_EQ(varianceFromDef, varianceFromWel);
         helper.GenerateNumbers(j);
      }
   }
}

TEST_P(RDFSimpleTests, StandardDeviationCollections)
{
   RDataFrame tdf(3);
   auto stdDev = tdf.Define("vector",
                            []() {
                               std::vector<int> v(3);
                               v[0] = 0;
                               v[1] = 1;
                               v[2] = 2;
                               return v;
                            })
                    .StdDev<std::vector<int>>("vector");
   EXPECT_DOUBLE_EQ(*stdDev, 0.86602540378443871);
}

TEST_P(RDFSimpleTests, StandardDeviationZero)
{
   RDataFrame rd1(8);
   auto stdDev = rd1.Define("b1", []() { return 0; }).StdDev("b1");
   EXPECT_DOUBLE_EQ(*stdDev, 0);
}

TEST_P(RDFSimpleTests, StandardDeviationOne)
{
   RDataFrame rd1(1);
   auto stdDev = rd1.Define("b1", []() { return 1; }).StdDev("b1");
   EXPECT_DOUBLE_EQ(*stdDev, 0);
}

TEST_P(RDFSimpleTests, StandardDeviationEmpty)
{
   RDataFrame rd1(0);
   auto stdDev = rd1.Define("b1", []() { return 0; }).StdDev("b1");
   EXPECT_DOUBLE_EQ(*stdDev, 0);
}

TEST(RDFSimpleTests, SumOfStrings)
{
   auto df = RDataFrame(2).Define("str", []() -> std::string { return "bla"; });
   EXPECT_EQ(*df.Sum<std::string>("str"), "blabla");
}


TEST(RDFSimpleTests, GenVector)
{
   // The leading underscore of "_hh" tests against ROOT-10305.
   ROOT::RDataFrame t(1);
   auto aa = t.Define("_hh", "ROOT::Math::PtEtaPhiMVector(1,1,1,1)").Define("h", "_hh.Rapidity()");
   auto m = aa.Mean("h");
   EXPECT_TRUE(0 != *m);
}

TEST(RDFSimpleTests, AutomaticNamesOfHisto1DAndGraph)
{
   auto df = RDataFrame(1).Define("x", [](){return 1;})
                          .Define("y", [](){return 1;});
   auto hx = df.Histo1D("x");
   auto hxy = df.Histo1D("x", "y");
   auto gxy = df.Graph("x", "y");

   EXPECT_STREQ(hx->GetName(), "x");
   EXPECT_STREQ(hx->GetTitle(), "x");
   EXPECT_STREQ(hx->GetXaxis()->GetTitle(), "x");
   EXPECT_STREQ(hx->GetYaxis()->GetTitle(), "count");
   EXPECT_STREQ(hxy->GetName(), "x_weighted_y");
   EXPECT_STREQ(hxy->GetTitle(), "x, weights: y");
   EXPECT_STREQ(hxy->GetXaxis()->GetTitle(), "x");
   EXPECT_STREQ(hxy->GetYaxis()->GetTitle(), "count * y");
   EXPECT_STREQ(gxy->GetName(), "x_vs_y");
   EXPECT_STREQ(gxy->GetTitle(), "x vs y");
   EXPECT_STREQ(gxy->GetXaxis()->GetTitle(), "x");
   EXPECT_STREQ(gxy->GetYaxis()->GetTitle(), "y");

}

TEST_P(RDFSimpleTests, DifferentTreesInDifferentThreads)
{
   const auto filename = "DifferentTreesInDifferentThreads.root";
   const auto treename = "mytree";
   {
      auto df = RDataFrame(64)
                   .Define("x", []() { return 1; })
                   .Define("y", []() { return 1; })
                   .Snapshot<int, int>(treename, filename, {"x", "y"}, {"RECREATE", ROOT::kZLIB, 4, 2, 99, false});
   }

   TFile f(filename);
   auto t = f.Get<TTree>(treename);
   RDataFrame df(*t);
   *df.Define("xy", [](int x, int y) { return x * y; }, {"x", "y"})
       .Filter([](int xy) { return xy > 0; }, {"xy"})
       .Count();

   gSystem->Unlink(filename);
}

TEST_P(RDFSimpleTests, HistosOneWeightPerEvent)
{
   using floats = std::vector<float>;
   auto df = RDataFrame(1);
   auto d = df.Define("v0", [](){floats v({1,2,3});return v;})
              .Define("v1", [](){floats v({4,5,6});return v;})
              .Define("v2", [](){floats v({7,8,9});return v;})
              .Define("w",[](){return 3;});

   auto h1 = d.Histo1D<floats, int>("v0","w");
   EXPECT_DOUBLE_EQ(h1->GetMean(), 2.);
   auto h2 = d.Histo2D<floats, floats, int>({"","",16,0,16,16,0,16}, "v0", "v1", "w");
   EXPECT_DOUBLE_EQ(h2->GetMean(), 2.);
   auto h3 = d.Histo3D<floats, floats, floats, int>({"","",16,0,16,16,0,16,16,0,16},"v0", "v1", "v2", "w");
   EXPECT_DOUBLE_EQ(h3->GetMean(), 2.);
}

TEST_P(RDFSimpleTests, ManyRangesPerWorker)
{
   auto filename = "ManyRangesPerWorker_file.root";
   {
      ROOT::RDataFrame(184).Define("i",[](){return 0;})
        .Snapshot<int>("t",filename,{"i"},{"RECREATE", ROOT::kZLIB, 1, 1, 99, false});
   }
   ROOT::RDataFrame("t",filename).Mean<int>("i");
   gSystem->Unlink(filename);
}

// ROOT-9736
TEST_P(RDFSimpleTests, NonExistingFile)
{
   ROOT::RDataFrame r("myTree", "nonexistingfile.root");

   TString expecteddiag;
   expecteddiag.Form("file %s/nonexistingfile.root does not exist", gSystem->pwd());

   // We try to use the tree for jitting: an exception is thrown
   ROOT_EXPECT_ERROR(EXPECT_ANY_THROW(r.Filter("inventedVar > 0")), "TFile::TFile", expecteddiag.Data());
}

// ROOT-10549: check we throw if a file is unreadable
TEST_P(RDFSimpleTests, NonExistingFileInChain)
{
   const auto filename = "rdf_nonexistingfileinchain.root";
   ROOT::RDataFrame(1).Define("x", [] { return 10; }).Snapshot<int>("t", filename, {"x"});

   ROOT::RDataFrame df("t", {filename, "doesnotexist.root"});

   const auto errmsg = "file %s/doesnotexist.root does not exist";
   TString expecteddiag;
   expecteddiag.Form(errmsg, gSystem->pwd());
   // in the single-thread case the error happens when TTreeReader is calling LoadTree the first time
   // otherwise we notice the file does not exist beforehand, e.g. in TTreeProcessorMT
   if (!ROOT::IsImplicitMTEnabled())
      expecteddiag += "\nWarning in <TTreeReader::SetEntryBase()>: There was an issue opening the last file associated "
                      "to the TChain being processed.";

   bool exceptionCaught = false;
   try {
      ROOT_EXPECT_ERROR(df.Count().GetValue(), "TFile::TFile", expecteddiag.Data());
   } catch (const std::runtime_error &e) {
      const std::string expected_msg =
         ROOT::IsImplicitMTEnabled()
            ? "TTreeProcessorMT::Process: an error occurred while opening file \"doesnotexist.root\""
            : "An error was encountered while processing the data. TTreeReader status code is: 5";
      EXPECT_EQ(e.what(), expected_msg);
      exceptionCaught = true;
   }
   EXPECT_TRUE(exceptionCaught);

   gSystem->Unlink(filename);
}

TEST_P(RDFSimpleTests, Stats)
{
   ROOT::RDataFrame r(256);
   auto rr = r.Define("v", [](ULong64_t e){return e;}, {"rdfentry_"})
              .Define("vec_v", [](ULong64_t e){return std::vector<ULong64_t>({e, e+1, e+2});}, {"v"})
              .Define("w", [](ULong64_t e){return 1./(e+1);}, {"v"})
              .Define("vec_w", [](double w){return std::vector<double>({w, w+1, w+2});}, {"w"})
              .Define("one", [](){return 1.;})
              .Define("ones", [](){return std::vector<double>({1.,1.,1.});});

   auto s0 = rr.Stats("v");
   auto s0c = rr.Stats<ULong64_t>("v");
   auto m0 = rr.Mean<ULong64_t>("v");
   auto v0 = rr.StdDev<ULong64_t>("v");
   auto s0prime = rr.Stats("v", "one");
   auto s0primec = rr.Stats<ULong64_t, double>("v", "one");
   auto s0w = rr.Stats("v", "w");
   auto s0wc = rr.Stats<ULong64_t, double>("v", "w");
   auto s1 = rr.Stats("vec_v");
   auto s1c = rr.Stats<std::vector<ULong64_t>>("vec_v");
   auto m1 = rr.Mean<std::vector<ULong64_t>>("vec_v");
   auto v1 = rr.StdDev<std::vector<ULong64_t>>("vec_v");
   auto s1w = rr.Stats("vec_v", "vec_w");
   auto s1wc = rr.Stats<std::vector<ULong64_t>, std::vector<double>>("vec_v", "vec_w");
   auto s1prime0 = rr.Stats("vec_v", "one");
   auto s1prime1 = rr.Stats("vec_v", "ones");
   auto s1prime0c = rr.Stats<std::vector<ULong64_t>, double>("vec_v", "one");
   auto s1prime1c = rr.Stats<std::vector<ULong64_t>, std::vector<double>>("vec_v", "ones");

   // Checks
   EXPECT_FLOAT_EQ(s0->GetMean(), 127.5);
   EXPECT_FLOAT_EQ(s0->GetMean(), s0c->GetMean());
   EXPECT_FLOAT_EQ(s0w->GetMean(), 40.800388);
   EXPECT_FLOAT_EQ(s0w->GetMean(), s0wc->GetMean());
   EXPECT_FLOAT_EQ(s0->GetMean(), s0prime->GetMean());
   EXPECT_FLOAT_EQ(s0->GetMean(), *m0);
   EXPECT_FLOAT_EQ(s0->GetRMS(), *v0);
   EXPECT_FLOAT_EQ(s1->GetMean(), 128.5);
   EXPECT_FLOAT_EQ(s1->GetMean(), s1c->GetMean());
   EXPECT_FLOAT_EQ(s1w->GetMean(), 127.12541);
   EXPECT_FLOAT_EQ(s1w->GetMean(), s1wc->GetMean());
   EXPECT_FLOAT_EQ(s1->GetMean(), *m1);
   EXPECT_FLOAT_EQ(s1->GetRMS(), *v1);
   EXPECT_FLOAT_EQ(s1->GetMean(), s1prime0->GetMean());
   EXPECT_FLOAT_EQ(s1->GetMean(), s1prime1->GetMean());
   EXPECT_FLOAT_EQ(s1->GetMean(), s1prime0c->GetMean());
   EXPECT_FLOAT_EQ(s1->GetMean(), s1prime1c->GetMean());

   // Check for the unsupported case
   EXPECT_ANY_THROW(rr.Stats<ULong64_t>("v", "one"));
}

// ROOT-10092
TEST(RDFSimpleTests, ScalarValuesCollectionWeights)
{
   ROOT::RDataFrame r(1);
   auto h = r.Define("x", [](){return 10;})
             .Define("y", [](){return ROOT::RVec<int>{1,2,3}; })
             .Histo1D<int, ROOT::RVec<int>>("x","y");

   // Check that the exception is thrown
   EXPECT_ANY_THROW(*h);
}

TEST_P(RDFSimpleTests, ChainWithDifferentTreeNames)
{
   const auto fname1 = "test_chainwithdifferenttreenames_1.root";
   const auto fname2 = "test_chainwithdifferenttreenames_2.root";
   {
      ROOT::RDataFrame(10).Define("x", [] { return 1; }).Snapshot<int>("t1", fname1, {"x"});
      ROOT::RDataFrame(10).Define("x", [] { return 3; }).Snapshot<int>("t2", fname2, {"x"});
   }

   // add trees to chain
   TChain c("t1");
   c.Add(fname1);
   c.Add((std::string(fname2) + "/t2").c_str());

   // pass chain to RDF and process trees with different names
   ROOT::RDataFrame df2(c);
   EXPECT_DOUBLE_EQ(*df2.Mean("x"), 2);

   gSystem->Unlink(fname1);
   gSystem->Unlink(fname2);
}

TEST_P(RDFSimpleTests, WritingToFundamentalType)
{
   EXPECT_THROW(ROOT::RDataFrame(1).Define("x", [] { return 1; }).Filter("x = 42"), std::runtime_error);
}

// run single-thread tests
INSTANTIATE_TEST_SUITE_P(Seq, RDFSimpleTests, ::testing::Values(false));

// run multi-thread tests
#ifdef R__USE_IMT
   INSTANTIATE_TEST_SUITE_P(MT, RDFSimpleTests, ::testing::Values(true));
#endif
