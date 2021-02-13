#include "ROOTUnitTestSupport.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/TSeq.hxx"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include <TInterpreter.h>
#include "TTree.h"
#include "TChain.h"
#include "gtest/gtest.h"
#include <limits>
#include <memory>
using namespace ROOT;         // RDataFrame
using namespace ROOT::RDF;    // RInterface
using namespace ROOT::VecOps; // RVec
using namespace ROOT::Detail::RDF;          // RLoopManager

/********* FIXTURES *********/
// fixture that provides a RDF with no data-source and a single integer column "ans" with value 42
class RDFSnapshot : public ::testing::Test {
protected:
   const ULong64_t nEvents = 100ull; // must be initialized before fLoopManager

private:
   RDataFrame fTdf;
   RInterface<RLoopManager> DefineAns()
   {
      return fTdf.Define("ans", []() { return 42; });
   }

protected:
   RDFSnapshot() : fTdf(nEvents), tdf(DefineAns()) {}
   RInterface<RLoopManager> tdf;
};

#ifdef R__USE_IMT
struct TIMTEnabler {
   TIMTEnabler(unsigned int nSlots) { ROOT::EnableImplicitMT(nSlots); }
   ~TIMTEnabler() { ROOT::DisableImplicitMT(); }
};

// fixture that enables implicit MT and provides a RDF with no data-source and a single column "x" containing
// normal-distributed doubles
class RDFSnapshotMT : public ::testing::Test {
protected:
   const ULong64_t kNEvents = 100ull; // must be initialized before fLoopManager
   const unsigned int kNSlots = 4u;

private:
   TIMTEnabler fIMTEnabler;
   RDataFrame fTdf;
   RInterface<RLoopManager> DefineAns()
   {
      return fTdf.Define("ans", []() { return 42; });
   }

protected:
   RDFSnapshotMT() : fIMTEnabler(kNSlots), fTdf(kNEvents), tdf(DefineAns()) {}
   RInterface<RLoopManager> tdf;
};
#endif // R__USE_IMT

// fixture that provides fixed and variable sized arrays as RDF columns
class RDFSnapshotArrays : public ::testing::Test {
protected:
   const static unsigned int kNEvents = 10u;
   static const std::vector<std::string> kFileNames;

   static void SetUpTestCase()
   {
      // write files containing c-arrays
      const auto eventsPerFile = kNEvents / kFileNames.size();
      auto curEvent = 0u;
      for (const auto &fname : kFileNames) {
         TFile f(fname.c_str(), "RECREATE");
         TTree t("arrayTree", "arrayTree");

         // doubles, floats
         const unsigned int fixedSize = 4u;
         float fixedSizeArr[fixedSize];
         t.Branch("fixedSizeArr", fixedSizeArr, ("fixedSizeArr[" + std::to_string(fixedSize) + "]/F").c_str());
         unsigned int size = 0u;
         t.Branch("size", &size);
         double *varSizeArr = new double[eventsPerFile * 100u];
         t.Branch("varSizeArr", varSizeArr, "varSizeArr[size]/D");

         // bools. std::vector<bool> makes bool treatment in RDF special
         bool fixedSizeBoolArr[fixedSize];
         t.Branch("fixedSizeBoolArr", fixedSizeBoolArr,
                  ("fixedSizeBoolArr[" + std::to_string(fixedSize) + "]/O").c_str());
         bool *varSizeBoolArr = new bool[eventsPerFile * 100u];
         t.Branch("varSizeBoolArr", varSizeBoolArr, "varSizeBoolArr[size]/O");

         // for each event, fill array elements
         for (auto i : ROOT::TSeqU(eventsPerFile)) {
            for (auto j : ROOT::TSeqU(4)) {
               fixedSizeArr[j] = curEvent * j;
               fixedSizeBoolArr[j] = j % 2 == 0;
            }
            size = (i + 1) * 100u;
            for (auto j : ROOT::TSeqU(size)) {
               varSizeArr[j] = curEvent * j;
               varSizeBoolArr[j] = j % 2 == 0;
            }
            t.Fill();
            ++curEvent;
         }
         t.Write();

         delete[] varSizeArr;
         delete[] varSizeBoolArr;
      }
   }

   static void TearDownTestCase()
   {
      for (const auto &fname : kFileNames)
         gSystem->Unlink(fname.c_str());
   }
};
const std::vector<std::string> RDFSnapshotArrays::kFileNames = {"test_snapshotarray1.root", "test_snapshotarray2.root"};

/********* SINGLE THREAD TESTS ***********/

TEST_F(RDFSnapshot, SnapshotCallAmbiguities)
{
   auto filename = "Snapshot_interface.root";

   tdf.Snapshot("t", filename, "an.*");
   tdf.Snapshot("t", filename, {"ans"});
   tdf.Snapshot("t", filename, {{"ans"}});

   gSystem->Unlink(filename);
}

// Test for ROOT-9210
TEST_F(RDFSnapshot, Snapshot_aliases)
{
   const auto alias0 = "myalias0";
   const auto alias1 = "myalias1";
   auto tdfa = tdf.Alias(alias0, "ans");
   auto tdfb = tdfa.Define("vec", [] { return RVec<int>{1,2,3}; }).Alias(alias1, "vec");
   testing::internal::CaptureStderr();
   auto snap = tdfb.Snapshot<int, RVec<int>>("mytree", "Snapshot_aliases.root", {alias0, alias1});
   std::string err = testing::internal::GetCapturedStderr();
   EXPECT_TRUE(err.empty()) << err;
   auto names = snap->GetColumnNames();
   EXPECT_EQ(2U, names.size());
   EXPECT_EQ(names, std::vector<std::string>({alias0, alias1}));

   auto takenCol = snap->Alias("a", alias0).Take<int>("a");
   for (auto i : takenCol) {
      EXPECT_EQ(42, i);
   }
}

// Test for ROOT-9122
TEST_F(RDFSnapshot, Snapshot_nocolumnmatch)
{
   const auto fname = "snapshotnocolumnmatch.root";
   RDataFrame d(1);
   auto op = [&](){
      d.Snapshot("t", fname, "x");
   };
   EXPECT_ANY_THROW(op());
   gSystem->Unlink(fname);
}

void TestSnapshotUpdate(RInterface<RLoopManager> &tdf, const std::string &outfile, const std::string &tree1,
                        const std::string &tree2, bool overwriteIfExists)
{
   // test snapshotting two trees to the same file opened in "UPDATE" mode
   auto df = tdf.Define("x", [] { return 10; });
   auto s1 = df.Snapshot<int>(tree1, outfile, {"x"});

   auto c1 = s1->Count();
   auto mean1 = s1->Mean<int>("x");
   EXPECT_EQ(100ull, *c1);
   EXPECT_DOUBLE_EQ(10., *mean1);

   RSnapshotOptions opts;
   opts.fMode = "UPDATE";
   opts.fOverwriteIfExists = overwriteIfExists;
   auto s2 = ROOT::RDataFrame(50ull).Define("x", [] { return 10; })
                                    .Snapshot<int>(tree2, outfile, {"x"}, opts);

   auto c2 = s2->Count();
   auto mean2 = s2->Mean<int>("x");
   EXPECT_EQ(50ull, *c2);
   EXPECT_DOUBLE_EQ(10., *mean2);

   // check that the output file contains both trees
   std::unique_ptr<TFile> f(TFile::Open(outfile.c_str()));
   EXPECT_NE(nullptr, f->Get<TTree>(tree1.c_str()));
   EXPECT_NE(nullptr, f->Get<TTree>(tree2.c_str()));

   // clean-up
   gSystem->Unlink(outfile.c_str());
}

TEST_F(RDFSnapshot, Snapshot_update_diff_treename)
{
   // test snapshotting two trees with different names
   TestSnapshotUpdate(tdf, "snap_update_difftreenames.root", "t1", "t2", false);
}

TEST_F(RDFSnapshot, Snapshot_update_same_treename)
{
   bool exceptionCaught = false;
   try {
      // test snapshotting two trees with same name
      TestSnapshotUpdate(tdf, "snap_update_sametreenames.root", "t", "t", false);
   } catch (const std::invalid_argument &e) {
      const std::string msg =
         "Snapshot: tree \"t\" already present in file \"snap_update_sametreenames.root\". If you want to delete the "
         "original tree and write another, please set RSnapshotOptions::fOverwriteIfExists to true.";
      EXPECT_EQ(e.what(), msg);
      exceptionCaught = true;
   }
   EXPECT_TRUE(exceptionCaught);
}

TEST_F(RDFSnapshot, Snapshot_update_overwrite)
{
   // test snapshotting two trees with different names
   TestSnapshotUpdate(tdf, "snap_update_overwrite.root", "t", "t", true);
}

void test_snapshot_options(RInterface<RLoopManager> &tdf)
{
   RSnapshotOptions opts;
   opts.fAutoFlush = 10;
   opts.fMode = "RECREATE";
   opts.fCompressionLevel = 6;

   const auto outfile = "snapshot_test_opts.root";
   for (auto algorithm : {ROOT::kZLIB, ROOT::kLZMA, ROOT::kLZ4, ROOT::kZSTD}) {
      opts.fCompressionAlgorithm = algorithm;

      auto s = tdf.Snapshot<int>("t", outfile, {"ans"}, opts);

      auto c = s->Count();
      auto min = s->Min<int>("ans");
      auto max = s->Max<int>("ans");
      auto mean = s->Mean<int>("ans");
      EXPECT_EQ(100ull, *c);
      EXPECT_EQ(42, *min);
      EXPECT_EQ(42, *max);
      EXPECT_EQ(42, *mean);

      std::unique_ptr<TFile> f(TFile::Open("snapshot_test_opts.root"));

      EXPECT_EQ(algorithm, f->GetCompressionAlgorithm());
      EXPECT_EQ(6, f->GetCompressionLevel());
   }

   // clean-up
   gSystem->Unlink(outfile);
}

TEST_F(RDFSnapshot, Snapshot_action_with_options)
{
   test_snapshot_options(tdf);
}

void checkSnapshotArrayFile(RResultPtr<RInterface<RLoopManager>> &df, unsigned int kNEvents)
{
   // fixedSizeArr and varSizeArr are RResultPtr<vector<vector<T>>>
   auto fixedSizeArr = df->Take<RVec<float>>("fixedSizeArr");
   auto varSizeArr = df->Take<RVec<double>>("varSizeArr");
   auto fixedSizeBoolArr = df->Take<RVec<bool>>("fixedSizeBoolArr");
   auto varSizeBoolArr = df->Take<RVec<bool>>("varSizeBoolArr");
   auto size = df->Take<unsigned int>("size");

   // check contents of fixed sized arrays
   const auto nEvents = fixedSizeArr->size();
   const auto fixedSizeSize = fixedSizeArr->front().size();
   EXPECT_EQ(nEvents, kNEvents);
   EXPECT_EQ(fixedSizeSize, 4u);
   for (auto i = 0u; i < nEvents; ++i) {
      for (auto j = 0u; j < fixedSizeSize; ++j) {
         EXPECT_DOUBLE_EQ(fixedSizeArr->at(i).at(j), i * j);
         EXPECT_EQ(fixedSizeBoolArr->at(i).at(j), j % 2 == 0);
      }
   }

   // check contents of variable sized arrays
   for (auto i = 0u; i < nEvents; ++i) {
      const auto thisSize = size->at(i);
      const auto &dv = varSizeArr->at(i);
      const auto &bv = varSizeBoolArr->at(i);
      EXPECT_EQ(thisSize, dv.size());
      EXPECT_EQ(thisSize, bv.size());
      for (auto j = 0u; j < thisSize; ++j) {
         EXPECT_DOUBLE_EQ(dv[j], i * j);
         EXPECT_EQ(bv[j], j % 2 == 0);
      }
   }
}

TEST_F(RDFSnapshotArrays, SingleThread)
{
   RDataFrame tdf("arrayTree", kFileNames);
   // template Snapshot
   // "size" _must_ be listed before "varSizeArr"!
   auto dt = tdf.Snapshot<RVec<float>, unsigned int, RVec<double>, RVec<bool>, RVec<bool>>(
      "outTree", "test_snapshotRVecoutST.root",
      {"fixedSizeArr", "size", "varSizeArr", "varSizeBoolArr", "fixedSizeBoolArr"});

   checkSnapshotArrayFile(dt, kNEvents);
}

TEST_F(RDFSnapshotArrays, SingleThreadJitted)
{
   RDataFrame tdf("arrayTree", kFileNames);
   // jitted Snapshot
   // "size" _must_ be listed before "varSizeArr"!
   auto dj = tdf.Snapshot("outTree", "test_snapshotRVecoutSTJitted.root",
                          {"fixedSizeArr", "size", "varSizeArr", "varSizeBoolArr", "fixedSizeBoolArr"});

   checkSnapshotArrayFile(dj, kNEvents);
}

void WriteColsWithCustomTitles(const std::string &tname, const std::string &fname)
{
   int i;
   float f;
   int a[2];
   TFile file(fname.c_str(), "RECREATE");
   TTree t(tname.c_str(), tname.c_str());
   auto b = t.Branch("float", &f);
   b->SetTitle("custom title");
   b = t.Branch("i", &i);
   b->SetTitle("custom title");
   b = t.Branch("arrint", &a, "arrint[2]/I");
   b->SetTitle("custom title");
   b = t.Branch("vararrint", &a, "vararrint[i]/I");
   b->SetTitle("custom title");

   i = 1;
   a[0] = 42;
   a[1] = 84;
   f = 4.2;
   t.Fill();

   i = 2;
   f = 8.4;
   t.Fill();

   t.Write();
}

void CheckColsWithCustomTitles(unsigned long long int entry, int i, const RVec<int> &arrint,
                               const RVec<int> &vararrint, float f)
{
   if (entry == 0) {
      EXPECT_EQ(i, 1);
      EXPECT_EQ(arrint.size(), 2u);
      EXPECT_EQ(arrint[0], 42);
      EXPECT_EQ(arrint[1], 84);
      EXPECT_EQ(vararrint.size(), 1u);
      EXPECT_EQ(vararrint[0], 42);
      EXPECT_FLOAT_EQ(f, 4.2f);
   } else if (entry == 1) {
      EXPECT_EQ(i, 2);
      EXPECT_EQ(arrint.size(), 2u);
      EXPECT_EQ(arrint[0], 42);
      EXPECT_EQ(arrint[1], 84);
      EXPECT_EQ(vararrint.size(), 2u);
      EXPECT_EQ(vararrint[0], 42);
      EXPECT_EQ(vararrint[1], 84);
      EXPECT_FLOAT_EQ(f, 8.4f);
   } else
      throw std::runtime_error("tree has more entries than expected");
}

TEST(RDFSnapshotMore, ColsWithCustomTitles)
{
   const auto fname = "colswithcustomtitles.root";
   const auto tname = "t";

   // write test tree
   WriteColsWithCustomTitles(tname, fname);

   // read and write test tree with RDF
   RDataFrame d(tname, fname);
   const std::string prefix = "snapshotted_";
   auto res_tdf =
      d.Snapshot<int, float, RVec<int>, RVec<int>>(tname, prefix + fname, {"i", "float", "arrint", "vararrint"});

   // check correct results have been written out
   res_tdf->Foreach(CheckColsWithCustomTitles, {"tdfentry_", "i", "arrint", "vararrint", "float"});

   // clean-up
   gSystem->Unlink(fname);
   gSystem->Unlink((prefix + fname).c_str());
}

TEST(RDFSnapshotMore, ReadWriteStdVec)
{
   // write a TFile containing a std::vector
   const auto fname = "readwritestdvec.root";
   const auto treename = "t";
   TFile f(fname, "RECREATE");
   TTree t(treename, treename);
   std::vector<int> v({42});
   std::vector<bool> vb({true, false, true}); // std::vector<bool> is special, not in a good way
   t.Branch("v", &v);
   t.Branch("vb", &vb);
   t.Fill();
   // as an extra test, make sure that the vector reallocates between first and second entry
   v = std::vector<int>(100000, 84);
   vb = std::vector<bool>(100000, true);
   t.Fill();
   t.Write();
   f.Close();

   auto outputChecker = [&treename](const char* filename){
      // check snapshot output
      TFile f2(filename);
      TTreeReader r(treename, &f2);
      TTreeReaderArray<int> rv(r, "v");
      TTreeReaderArray<bool> rvb(r, "vb");
      r.Next();
      EXPECT_EQ(rv.GetSize(), 1u);
      EXPECT_EQ(rv[0], 42);
      EXPECT_EQ(rvb.GetSize(), 3u);
      EXPECT_TRUE(rvb[0]);
      EXPECT_FALSE(rvb[1]);
      EXPECT_TRUE(rvb[2]);
      r.Next();
      EXPECT_EQ(rv.GetSize(), 100000u);
      EXPECT_EQ(rvb.GetSize(), 100000u);
      for (auto e : rv)
         EXPECT_EQ(e, 84);
      for (auto e : rvb)
         EXPECT_TRUE(e);
   };

   // read and write using RDataFrame

   const auto outfname1 = "out_readwritestdvec1.root";
   RDataFrame(treename, fname).Snapshot<std::vector<int>, std::vector<bool>>(treename, outfname1, {"v", "vb"});
   outputChecker(outfname1);

   const auto outfname2 = "out_readwritestdvec2.root";
   RDataFrame(treename, fname).Snapshot(treename, outfname2);
   outputChecker(outfname2);

   const auto outfname3 = "out_readwritestdvec3.root";
   RDataFrame(treename, fname).Snapshot<RVec<int>, RVec<bool>>(treename, outfname3, {"v", "vb"});
   outputChecker(outfname3);

   gSystem->Unlink(fname);
   gSystem->Unlink(outfname1);
   gSystem->Unlink(outfname2);
   gSystem->Unlink(outfname3);
}

void ReadWriteCarray(const char *outFileNameBase)
{
   // write a TFile containing a arrays
   std::string outFileNameBaseStr = outFileNameBase;
   const auto fname = outFileNameBaseStr + ".root";
   const auto treename = "t";
   TFile f(fname.c_str(), "RECREATE");
   TTree t(treename, treename);
   const auto maxArraySize = 100000U;
   auto size = 0;
   int v[maxArraySize];
   bool vb[maxArraySize];
   long int vl[maxArraySize];
   t.Branch("size", &size, "size/I");
   t.Branch("v", v, "v[size]/I");
   t.Branch("vb", vb, "vb[size]/O");
   t.Branch("vl", vl, "vl[size]/G");

   // Size 1
   size = 1;
   v[0] = 12;
   vb[0] = true;
   vl[0] = 8589934592; // 2**33
   t.Fill();

   // Size 0 (see ROOT-9860)
   size = 0;
   t.Fill();

   // Size 100k: this reallocates!
   size = maxArraySize;
   for (auto i : ROOT::TSeqU(size)) {
      v[i] = 84;
      vb[i] = true;
      vl[i] = 42;
   }
   t.Fill();

   // Size 3
   size = 3;
   v[0] = 42;
   v[1] = 43;
   v[2] = 44;
   vb[0] = true;
   vb[1] = false;
   vb[2] = true;
   vl[0] = -1;
   vl[1] = 0;
   vl[2] = 1;
   t.Fill();

   t.Write();
   f.Close();

   auto outputChecker = [&treename](const char *filename) {
      // check snapshot output
      TFile f2(filename);
      TTreeReader r(treename, &f2);
      TTreeReaderArray<int> rv(r, "v");
      TTreeReaderArray<bool> rvb(r, "vb");
      TTreeReaderArray<long int> rvl(r, "vl");

      // Size 1
      EXPECT_TRUE(r.Next());
      EXPECT_EQ(rv.GetSize(), 1u);
      EXPECT_EQ(rv[0], 12);
      EXPECT_EQ(rvb.GetSize(), 1u);
      EXPECT_TRUE(rvb[0]);
      EXPECT_EQ(rvl.GetSize(), 1u);
      EXPECT_EQ(rvl[0], 8589934592);

      // Size 0
      EXPECT_TRUE(r.Next());
      EXPECT_EQ(rv.GetSize(), 0u);
      EXPECT_EQ(rvb.GetSize(), 0u);
      EXPECT_EQ(rvl.GetSize(), 0u);

      // Size 100k
      EXPECT_TRUE(r.Next());
      EXPECT_EQ(rv.GetSize(), 100000u);
      EXPECT_EQ(rvb.GetSize(), 100000u);
      for (auto e : rv)
         EXPECT_EQ(e, 84);
      for (auto e : rvb)
         EXPECT_TRUE(e);
      for (auto e : rvl)
         EXPECT_EQ(e, 42);

      // Size 3
      EXPECT_TRUE(r.Next());
      EXPECT_EQ(rv.GetSize(), 3u);
      EXPECT_EQ(rv[0], 42);
      EXPECT_EQ(rv[1], 43);
      EXPECT_EQ(rv[2], 44);
      EXPECT_EQ(rvb.GetSize(), 3u);
      EXPECT_TRUE(rvb[0]);
      EXPECT_FALSE(rvb[1]);
      EXPECT_TRUE(rvb[2]);
      EXPECT_EQ(rvl.GetSize(), 3u);
      EXPECT_EQ(rvl[0], -1);
      EXPECT_EQ(rvl[1], 0);
      EXPECT_EQ(rvl[2], 1);

      EXPECT_FALSE(r.Next());
   };

   // read and write using RDataFrame
   const auto outfname1 = outFileNameBaseStr + "_out1.root";
   RDataFrame(treename, fname).Snapshot(treename, outfname1);
   outputChecker(outfname1.c_str());

   const auto outfname2 = outFileNameBaseStr + "_out2.root";
   RDataFrame(treename, fname)
      .Snapshot<int, RVec<int>, RVec<bool>, RVec<long int>>(treename, outfname2, {"size", "v", "vb", "vl"});
   outputChecker(outfname2.c_str());

   gSystem->Unlink(fname.c_str());
   gSystem->Unlink(outfname1.c_str());
   gSystem->Unlink(outfname2.c_str());
}

TEST(RDFSnapshotMore, ReadWriteCarray)
{
   ReadWriteCarray("ReadWriteCarray");
}

struct TwoInts {
   int a, b;
};

void WriteTreeWithLeaves(const std::string &treename, const std::string &fname)
{
   TFile f(fname.c_str(), "RECREATE");
   TTree t(treename.c_str(), treename.c_str());

   TwoInts ti{1, 2};
   t.Branch("v", &ti, "a/I:b/I");

   // TODO add checks for reading of multiple nested levels ("w.v.a")
   // when ROOT-9312 is solved and RDF supports "w.v.a" nested notation

   t.Fill();
   t.Write();
}

TEST(RDFSnapshotMore, ReadWriteNestedLeaves)
{
   const auto treename = "t";
   const auto fname = "readwritenestedleaves.root";
   WriteTreeWithLeaves(treename, fname);
   RDataFrame d(treename, fname);
   const auto outfname = "out_readwritenestedleaves.root";
   ROOT::RDF::RNode d2(d);
   ROOT_EXPECT_INFO((d2 = *d.Snapshot<int, int>(treename, outfname, {"v.a", "v.b"})), "Snapshot", "Column v.a will be saved as v_a\nInfo in <Snapshot>: Column v.b will be saved as v_b");
   EXPECT_EQ(d2.GetColumnNames(), std::vector<std::string>({"v_a", "v_b"}));
   auto check_a_b = [](int a, int b) {
      EXPECT_EQ(a, 1);
      EXPECT_EQ(b, 2);
   };
   d2.Foreach(check_a_b, {"v_a", "v_b"});
   gSystem->Unlink(fname);
   gSystem->Unlink(outfname);

   try {
      d.Define("v_a", [] { return 0; }).Snapshot<int, int>(treename, outfname, {"v.a", "v_a"});
   } catch (std::runtime_error &e) {
      const auto error_msg = "Column v.a would be written as v_a but this column already exists. Please use Alias to "
                             "select a new name for v.a";
      EXPECT_STREQ(e.what(), error_msg);
   }
}

TEST(RDFSnapshotMore, Lazy)
{
   const auto treename = "t";
   const auto fname0 = "lazy0.root";
   const auto fname1 = "lazy1.root";
   // make sure the file is not here beforehand
   gSystem->Unlink(fname0);
   RDataFrame d(1);
   auto v = 0U;
   auto genf = [&v](){++v;return 42;};
   RSnapshotOptions opts = {"RECREATE", ROOT::kZLIB, 0, 0, 99, true};
   auto ds = d.Define("c0", genf).Snapshot<int>(treename, fname0, {"c0"}, opts);
   EXPECT_EQ(v, 0U);
   EXPECT_TRUE(gSystem->AccessPathName(fname0)); // This returns FALSE if the file IS there
   auto ds2 = ds->Define("c1", genf).Snapshot<int>(treename, fname1, {"c1"}, opts);
   EXPECT_EQ(v, 1U);
   EXPECT_FALSE(gSystem->AccessPathName(fname0));
   EXPECT_TRUE(gSystem->AccessPathName(fname1));
   *ds2;
   EXPECT_EQ(v, 2U);
   EXPECT_FALSE(gSystem->AccessPathName(fname1));
   gSystem->Unlink(fname0);
   gSystem->Unlink(fname1);
}

TEST(RDFSnapshotMore, LazyNotTriggered)
{
   {
      auto d = ROOT::RDataFrame(1);
      ROOT::RDF::RSnapshotOptions opts;
      opts.fLazy = true;
      d.Snapshot<ULong64_t>("t", "foo.root", {"tdfentry_"}, opts);
   }
}

void CheckTClonesArrayOutput(const RVec<TH1D> &hvec)
{
   ASSERT_EQ(hvec.size(), 3);
   for (int i = 0; i < 3; ++i) {
      EXPECT_EQ(hvec[i].GetEntries(), 1);
      EXPECT_DOUBLE_EQ(hvec[i].GetMean(), i);
   }
}

void ReadWriteTClonesArray()
{
   {
      TClonesArray arr("TH1D", 3);
      for (int i = 0; i < 3; ++i) {
         auto *h = static_cast<TH1D *>(arr.ConstructedAt(i));
         h->SetBins(25, 0, 10);
         h->Fill(i);
      }
      TFile f("df_readwritetclonesarray.root", "recreate");
      TTree t("t", "t");
      t.Branch("arr", &arr);
      t.Fill();
      t.Write();
      f.Close();
   }

   {
      // write as TClonesArray
      auto out_df = ROOT::RDataFrame("t", "df_readwritetclonesarray.root")
                       .Snapshot<TClonesArray>("t", "df_readwriteclonesarray1.root", {"arr"});
      RVec<TH1D> hvec;

#ifndef NDEBUG
      ROOT_EXPECT_WARNING(
         hvec = out_df->Take<RVec<TH1D>>("arr")->at(0), "RTreeColumnReader::Get",
         "Branch arr hangs from a non-split branch. A copy is being performed in order to properly read the content.");
#else
      ROOT_EXPECT_NODIAG(hvec = out_df->Take<RVec<TH1D>>("arr")->at(0));
#endif
      CheckTClonesArrayOutput(hvec);
   }

   // FIXME uncomment when ROOT-10801 is solved
   //{
   //   gInterpreter->GenerateDictionary("vector<TH1D,ROOT::Detail::VecOps::RAdoptAllocator<TH1D>>",
   //                                    "vector;TH1D.h;ROOT/RVec.hxx");
   //   // write as RVecs
   //   auto out_df = ROOT::RDataFrame("t", "df_readwritetclonesarray.root")
   //                    .Snapshot<RVec<TH1D>>("t", "df_readwriteclonesarray2.root", {"arr"});
   //   const auto hvec = out_df->Take<RVec<TH1D>>("arr")->at(0);
   //   CheckTClonesArrayOutput(hvec);
   //}

   {
      // write as Snapshot wants
      auto out_df =
         ROOT::RDataFrame("t", "df_readwritetclonesarray.root").Snapshot("t", "df_readwriteclonesarray3.root", {"arr"});
      RVec<TH1D> hvec;
#ifndef NDEBUG
      ROOT_EXPECT_WARNING(
         hvec = out_df->Take<RVec<TH1D>>("arr")->at(0), "RTreeColumnReader::Get",
         "Branch arr hangs from a non-split branch. A copy is being performed in order to properly read the content.");
#else
      ROOT_EXPECT_NODIAG(hvec = out_df->Take<RVec<TH1D>>("arr")->at(0));
#endif
      CheckTClonesArrayOutput(hvec);
   }

   gSystem->Unlink("df_readwritetclonesarray.root");
   gSystem->Unlink("df_readwriteclonesarray1.root");
   gSystem->Unlink("df_readwriteclonesarray2.root");
   gSystem->Unlink("df_readwriteclonesarray3.root");
}

TEST(RDFSnapshotMore, TClonesArray)
{
   ReadWriteTClonesArray();
}

// ROOT-10702
TEST(RDFSnapshotMore, CompositeTypeWithNameClash)
{
   const auto fname = "snap_compositetypewithnameclash.root";
   gInterpreter->Declare("struct Int { int x; };");
   ROOT::RDataFrame df(3);
   auto snap_df = df.Define("i", "Int{-1};").Define("x", [] { return 1; }).Snapshot("t", fname);
   EXPECT_EQ(snap_df->Sum<int>("x").GetValue(), 3); // prints -3 if the wrong "x" is written out
   EXPECT_EQ(snap_df->Sum<int>("i.x").GetValue(), -3);

   gSystem->Unlink(fname);
}

// Test that we error out gracefully in case the output file specified for a Snapshot cannot be opened
TEST(RDFSnapshotMore, ForbiddenOutputFilename)
{
   ROOT::RDataFrame df(4);
   const auto out_fname = "/definitely/not/a/valid/path/f.root";

   // Compiled
   try {
      ROOT_EXPECT_SYSERROR(df.Snapshot<unsigned int>("t", out_fname, {"rdfslot_"}), "TFile::TFile",
                        "file /definitely/not/a/valid/path/f.root can not be opened No such file or directory")
   } catch (const std::runtime_error &e) {
      EXPECT_STREQ(e.what(), "Snapshot: could not create output file /definitely/not/a/valid/path/f.root");
   }

   // Jitted
   // If some other test case called EnableThreadSafety, the error printed here is of the form
   // "SysError in <TFile::TFile>: file /definitely/not/a/valid/path/f.root can not be opened No such file or directory\nError in <TReentrantRWLock::WriteUnLock>: Write lock already released for 0x55f179989378\n"
   // but the address printed changes every time
   EXPECT_THROW(df.Snapshot("t", out_fname, {"rdfslot_"}), std::runtime_error);
}

/********* MULTI THREAD TESTS ***********/
#ifdef R__USE_IMT
TEST_F(RDFSnapshotMT, Snapshot_update_diff_treename)
{
   // test snapshotting two trees with different names
   TestSnapshotUpdate(tdf, "snap_update_difftreenames.root", "t1", "t2", false);
}

TEST_F(RDFSnapshotMT, Snapshot_update_same_treename)
{
   bool exceptionCaught = false;
   try {
      // test snapshotting two trees with same name
      TestSnapshotUpdate(tdf, "snap_update_sametreenames.root", "t", "t", false);
   } catch (const std::invalid_argument &e) {
      const std::string msg =
         "Snapshot: tree \"t\" already present in file \"snap_update_sametreenames.root\". If you want to delete the "
         "original tree and write another, please set RSnapshotOptions::fOverwriteIfExists to true.";
      EXPECT_EQ(e.what(), msg);
      exceptionCaught = true;
   }
   EXPECT_TRUE(exceptionCaught);
}

TEST_F(RDFSnapshotMT, Snapshot_update_overwrite)
{
   // test snapshotting two trees with different names
   TestSnapshotUpdate(tdf, "snap_update_overwrite.root", "t", "t", true);
}

TEST_F(RDFSnapshotMT, Snapshot_action_with_options)
{
   test_snapshot_options(tdf);
}

TEST_F(RDFSnapshotMT, Reshuffled_friends)
{
   const auto fname = "snapshot_reshuffled_friends.root";
   tdf.Snapshot("t", fname);

   {
      // add reshuffled tree as friend
      TFile f(fname);
      TTree *t = f.Get<TTree>("t");
      TTree t2("t2", "t2");
      const auto expected = "Tree 't' has the kEntriesReshuffled bit set, and cannot be used as "
                            "friend nor can be added as a friend unless the main tree has a TTreeIndex on the friend "
                            "tree 't'. You can also unset the bit manually if you know what you are doing.";
      ROOT_EXPECT_ERROR(t2.AddFriend(t), "AddFriend", expected);
   }

   {
      // add friend to reshuffled tree
      TFile f(fname);
      TTree *t = f.Get<TTree>("t");
      TTree t2("t2", "t2");
      const auto expected = "Tree 't' has the kEntriesReshuffled bit set, and cannot be used as "
                            "friend nor can be added as a friend unless the main tree has a TTreeIndex on the friend "
                            "tree 't2'. You can also unset the bit manually if you know what you are doing.";
      ROOT_EXPECT_ERROR(t->AddFriend(&t2);, "AddFriend", expected);
   }
}

TEST(RDFSnapshotMore, ManyTasksPerThread)
{
   const auto nSlots = 4u;
   ROOT::EnableImplicitMT(nSlots);
   // make sure the file is not here beforehand
   gSystem->Unlink("snapshot_manytasks_out.root");

   // easiest way to be sure reading requires spawning of several tasks: create several input files
   const std::string inputFilePrefix = "snapshot_manytasks_";
   const auto tasksPerThread = 8u;
   const auto nInputFiles = nSlots * tasksPerThread;
   ROOT::RDataFrame d(1);
   auto dd = d.Define("x", []() { return 42; });
   for (auto i = 0u; i < nInputFiles; ++i)
      dd.Snapshot<int>("t", inputFilePrefix + std::to_string(i) + ".root", {"x"});

   // test multi-thread Snapshotting from many tasks per worker thread
   const auto outputFile = "snapshot_manytasks_out.root";
   ROOT::RDataFrame tdf("t", (inputFilePrefix + "*.root").c_str());
   tdf.Snapshot<int>("t", outputFile, {"x"});

   // check output contents
   ROOT::RDataFrame checkTdf("t", outputFile);
   auto c = checkTdf.Count();
   auto t = checkTdf.Take<int>("x");
   for (auto v : t)
      EXPECT_EQ(v, 42);
   EXPECT_EQ(*c, nInputFiles);

   // clean-up input files
   for (auto i = 0u; i < nInputFiles; ++i)
      gSystem->Unlink((inputFilePrefix + std::to_string(i) + ".root").c_str());
   gSystem->Unlink(outputFile);

   ROOT::DisableImplicitMT();
}

void checkSnapshotArrayFileMT(RResultPtr<RInterface<RLoopManager>> &df, unsigned int kNEvents)
{
   // fixedSizeArr and varSizeArr are RResultPtr<vector<vector<T>>>
   auto fixedSizeArr = df->Take<RVec<float>>("fixedSizeArr");
   auto varSizeArr = df->Take<RVec<double>>("varSizeArr");
   auto size = df->Take<unsigned int>("size");

   // multi-thread execution might have scrambled events w.r.t. the original file, so we just check overall properties
   const auto nEvents = fixedSizeArr->size();
   EXPECT_EQ(nEvents, kNEvents);
   // TODO check more!
}

TEST_F(RDFSnapshotArrays, MultiThread)
{
   ROOT::EnableImplicitMT(4);

   RDataFrame tdf("arrayTree", kFileNames);
   auto dt = tdf.Snapshot<RVec<float>, unsigned int, RVec<double>, RVec<bool>, RVec<bool>>(
      "outTree", "test_snapshotRVecoutMT.root",
      {"fixedSizeArr", "size", "varSizeArr", "varSizeBoolArr", "fixedSizeBoolArr"});

   checkSnapshotArrayFileMT(dt, kNEvents);

   ROOT::DisableImplicitMT();
}

TEST_F(RDFSnapshotArrays, MultiThreadJitted)
{
   ROOT::EnableImplicitMT(4);

   RDataFrame tdf("arrayTree", kFileNames);
   auto dj = tdf.Snapshot("outTree", "test_snapshotRVecoutMTJitted.root",
                          {"fixedSizeArr", "size", "varSizeArr", "varSizeBoolArr", "fixedSizeBoolArr"});

   checkSnapshotArrayFileMT(dj, kNEvents);

   ROOT::DisableImplicitMT();
}

TEST(RDFSnapshotMore, ColsWithCustomTitlesMT)
{
   const auto fname = "colswithcustomtitlesmt.root";
   const auto tname = "t";

   // write test tree
   WriteColsWithCustomTitles(tname, fname);

   // read and write test tree with RDF (in parallel)
   ROOT::EnableImplicitMT(4);
   RDataFrame d(tname, fname);
   const std::string prefix = "snapshotted_";
   auto res_tdf =
      d.Snapshot<int, float, RVec<int>, RVec<int>>(tname, prefix + fname, {"i", "float", "arrint", "vararrint"});

   // check correct results have been written out
   res_tdf->Foreach(CheckColsWithCustomTitles, {"tdfentry_", "i", "arrint", "vararrint", "float"});
   res_tdf->Foreach(CheckColsWithCustomTitles, {"rdfentry_", "i", "arrint", "vararrint", "float"});

   // clean-up
   gSystem->Unlink(fname);
   gSystem->Unlink((prefix + fname).c_str());
   ROOT::DisableImplicitMT();
}

TEST(RDFSnapshotMore, TreeWithFriendsMT)
{
   const auto fname = "treewithfriendsmt.root";
   RDataFrame(10).Define("x", []() { return 0; }).Snapshot<int>("t", fname, {"x"});

   ROOT::EnableImplicitMT();

   TFile file(fname);
   auto tree = file.Get<TTree>("t");
   TFile file2(fname);
   auto tree2 = file2.Get<TTree>("t");
   tree->AddFriend(tree2);

   const auto outfname = "out_treewithfriendsmt.root";
   RDataFrame df(*tree);
   df.Snapshot<int>("t", outfname, {"x"});
   ROOT::DisableImplicitMT();

   gSystem->Unlink(fname);
   gSystem->Unlink(outfname);
}

TEST(RDFSnapshotMore, JittedSnapshotAndAliasedColumns)
{
   ROOT::RDataFrame df(1);
   const auto fname = "out_aliasedcustomcolumn.root";
   // aliasing a custom column
   auto df2 = df.Define("x", [] { return 42; }).Alias("y", "x").Snapshot("t", fname, "y"); // must be jitted!
   EXPECT_EQ(df2->GetColumnNames(), std::vector<std::string>({"y"}));
   EXPECT_EQ(df2->Take<int>("y")->at(0), 42);

   // aliasing a column from a file
   const auto fname2 = "out_aliasedcustomcolumn2.root";
   df2->Alias("z", "y").Snapshot("t", fname2, "z");

   gSystem->Unlink(fname);
   gSystem->Unlink(fname2);
}


TEST(RDFSnapshotMore, LazyNotTriggeredMT)
{
   ROOT::EnableImplicitMT(4);
   const auto fname = "lazynottriggeredmt.root";
   {
      auto d = ROOT::RDataFrame(8);
      ROOT::RDF::RSnapshotOptions opts;
      opts.fLazy = true;
      d.Snapshot<ULong64_t, ULong64_t>("t", fname, {"tdfentry_", "rdfentry_"}, opts);
   }

   gSystem->Unlink(fname);
   ROOT::DisableImplicitMT();
}

TEST(RDFSnapshotMore, EmptyBuffersMT)
{
   const auto fname = "emptybuffersmt.root";
   const auto treename = "t";
   ROOT::EnableImplicitMT(4);
   ROOT::RDataFrame d(10);
   auto dd = d.DefineSlot("x", [](unsigned int s) { return s == 3 ? 0 : 1; })
               .Filter([](int x) { return x == 0; }, {"x"}, "f");
   auto r = dd.Report();
   dd.Snapshot<int>(treename, fname, {"x"});

   // check test sanity
   const auto passed = r->At("f").GetPass();
   EXPECT_GT(passed, 0u);

   // check result
   TFile f(fname);
   auto t = f.Get<TTree>(treename);
   EXPECT_EQ(t->GetListOfBranches()->GetEntries(), 1);
   EXPECT_EQ(t->GetEntries(), Long64_t(passed));

   ROOT::DisableImplicitMT();
   gSystem->Unlink(fname);
}

TEST(RDFSnapshotMore, ReadWriteCarrayMT)
{
   ROOT::EnableImplicitMT(4);
   ReadWriteCarray("ReadWriteCarrayMT");
   ROOT::DisableImplicitMT();
}

TEST(RDFSnapshotMore, TClonesArrayMT)
{
   TIMTEnabler _(4);
   ReadWriteTClonesArray();
}

// Test that we error out gracefully in case the output file specified for a Snapshot cannot be opened
TEST(RDFSnapshotMore, ForbiddenOutputFilenameMT)
{
   TIMTEnabler _(4);
   ROOT::RDataFrame df(4);
   const auto out_fname = "/definitely/not/a/valid/path/f.root";

   // Compiled
   try {
      const auto expected = "file /definitely/not/a/valid/path/f.root can not be opened No such file or directory";
      ROOT_EXPECT_SYSERROR(df.Snapshot<unsigned int>("t", out_fname, {"rdfslot_"}), "TFile::TFile", expected);
   } catch (const std::runtime_error &e) {
      EXPECT_STREQ(e.what(), "Snapshot: could not create output file /definitely/not/a/valid/path/f.root");
   }

   // Jitted
   // the error printed here is
   // "SysError in <TFile::TFile>: file /definitely/not/a/valid/path/f.root can not be opened No such file or directory\nError in <TReentrantRWLock::WriteUnLock>: Write lock already released for 0x55f179989378\n"
   // but the address printed changes every time
   EXPECT_THROW(df.Snapshot("t", out_fname, {"rdfslot_"}), std::runtime_error);
}

/**
 * Test against issue #6523 and #6640
 * Try to force `TTree::ChangeFile` behaviour. Within RDataFrame, this should
 * not happen and both sequential and multithreaded Snapshot should only create
 * one file.
 */
TEST(RDFSnapshotMore, SetMaxTreeSizeMT)
{
   // Set TTree max size to a low number. Normally this would trigger the
   // behaviour of TTree::ChangeFile, but not within RDataFrame.
   auto old_maxtreesize{TTree::GetMaxTreeSize()};
   TTree::SetMaxTreeSize(1000);

   // Create TTree, fill it and Snapshot (should create one single file).
   {
      TTree t{"T", "SetMaxTreeSize(1000)"};
      int x{};
      auto nentries{20000};

      t.Branch("x", &x, "x/I");

      for (auto i = 0; i < nentries; i++) {
         x = i;
         t.Fill();
      }

      ROOT::RDataFrame df{t};
      df.Snapshot<Int_t>("T", "rdfsnapshot_ttree_sequential_setmaxtreesize.root", {"x"});
   }

   // Create an RDF from the previously snapshotted file, then Snapshot again
   // with IMT enabled.
   {
      ROOT::EnableImplicitMT();

      ROOT::RDataFrame df{"T", "rdfsnapshot_ttree_sequential_setmaxtreesize.root"};
      df.Snapshot<Int_t>("T", "rdfsnapshot_imt_setmaxtreesize.root", {"x"});

      ROOT::DisableImplicitMT();
   }

   // Check the file for data integrity.
   {
      TFile f{"rdfsnapshot_imt_setmaxtreesize.root"};
      std::unique_ptr<TTree> t{f.Get<TTree>("T")};

      EXPECT_EQ(t->GetEntries(), 20000);

      int sum{0};
      int x{0};
      t->SetBranchAddress("x", &x);

      for (auto i = 0; i < t->GetEntries(); i++) {
         t->GetEntry(i);
         sum += x;
      }

      // sum(range(20000)) == 199990000
      EXPECT_EQ(sum, 199990000);
   }

   gSystem->Unlink("rdfsnapshot_ttree_sequential_setmaxtreesize.root");
   gSystem->Unlink("rdfsnapshot_imt_setmaxtreesize.root");

   // Reset TTree max size to its old value
   TTree::SetMaxTreeSize(old_maxtreesize);
}

#endif // R__USE_IMT

