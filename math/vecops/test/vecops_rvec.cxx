#include <gtest/gtest.h>
#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>
#include <Math/Vector4Dfwd.h>
#include <ROOT/RVec.hxx>
#include <ROOT/TSeq.hxx>
#include <TFile.h>
#include <TInterpreter.h>
#include <TTree.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <vector>
#include <sstream>
#include <cmath>

using namespace ROOT::VecOps;

void CheckEqual(const RVec<float> &a, const RVec<float> &b, std::string_view msg = "")
{
   const auto asize = a.size();
   const auto bsize = b.size();
   EXPECT_EQ(asize, bsize);
   for (unsigned int i = 0; i < asize; ++i) {
      EXPECT_FLOAT_EQ(a[i], b[i]) << msg;
   }
}

void CheckEqual(const RVec<double> &a, const RVec<double> &b, std::string_view msg = "")
{
   const auto asize = a.size();
   const auto bsize = b.size();
   EXPECT_EQ(asize, bsize);
   for (unsigned int i = 0; i < asize; ++i) {
      EXPECT_DOUBLE_EQ(a[i], b[i]) << msg;
   }
}

template <typename T, typename V>
void CheckEqual(const T &a, const V &b, std::string_view msg = "")
{
   const auto asize = a.size();
   const auto bsize = b.size();
   EXPECT_EQ(asize, bsize);
   for (unsigned int i = 0; i < asize; ++i) {
      EXPECT_EQ(a[i], b[i]) << msg;
   }
}

TEST(VecOps, DefaultCtor)
{
   ROOT::VecOps::RVec<int> v;
   EXPECT_EQ(v.size(), 0u);
}

TEST(VecOps, InitListCtor)
{
   ROOT::VecOps::RVec<int> v{1, 2, 3};
   EXPECT_EQ(v.size(), 3u);
   EXPECT_EQ(v[0], 1);
   EXPECT_EQ(v[1], 2);
   EXPECT_EQ(v[2], 3);
}

TEST(VecOps, CopyCtor)
{
   ROOT::VecOps::RVec<int> v1{1, 2, 3};
   ROOT::VecOps::RVec<int> v2(v1);
   EXPECT_EQ(v1.size(), 3u);
   EXPECT_EQ(v2.size(), 3u);
   EXPECT_EQ(v2[0], 1);
   EXPECT_EQ(v2[1], 2);
   EXPECT_EQ(v2[2], 3);
}

class TLeakChecker {
public:
   static bool fgDestroyed;
   ~TLeakChecker(){
      fgDestroyed = true;
   }
};
bool TLeakChecker::fgDestroyed = false;

// This test is here to check that emplace_back works
// also for T==bool. Indeed, on some platform, vector<bool>
// has no emplace_back. Notable examples are osx 10.14 and gcc 4.8
TEST(VecOps, EmplaceBack)
{
   ROOT::RVec<bool> vb; vb.emplace_back(true);
   ROOT::RVec<int> vi; vi.emplace_back(1);
}

TEST(VecOps, CopyCtorCheckNoLeak)
{
   ROOT::VecOps::RVec<TLeakChecker> ref;
   ref.emplace_back(TLeakChecker());
   ROOT::VecOps::RVec<TLeakChecker> proxy(ref.data(), ref.size());
   TLeakChecker::fgDestroyed = false;
   {
      auto v = proxy;
   }
   EXPECT_TRUE(TLeakChecker::fgDestroyed);

   TLeakChecker::fgDestroyed = false;
   ref.clear();
   EXPECT_TRUE(TLeakChecker::fgDestroyed);

}

TEST(VecOps, MoveCtor)
{
   ROOT::VecOps::RVec<int> v1{1, 2, 3};
   ROOT::VecOps::RVec<int> v2(std::move(v1));
   EXPECT_EQ(v2.size(), 3u);
}

TEST(VecOps, Conversion)
{
   ROOT::VecOps::RVec<float> fvec{1.0f, 2.0f, 3.0f};
   ROOT::VecOps::RVec<unsigned> uvec{1u, 2u, 3u};

   ROOT::VecOps::RVec<int>  ivec = uvec;
   ROOT::VecOps::RVec<long> lvec = ivec;

   EXPECT_EQ(1, ivec[0]);
   EXPECT_EQ(2, ivec[1]);
   EXPECT_EQ(3, ivec[2]);
   EXPECT_EQ(3u, ivec.size());
   EXPECT_EQ(1l, lvec[0]);
   EXPECT_EQ(2l, lvec[1]);
   EXPECT_EQ(3l, lvec[2]);
   EXPECT_EQ(3u, lvec.size());

   auto dvec1 = ROOT::VecOps::RVec<double>(fvec);
   auto dvec2 = ROOT::VecOps::RVec<double>(uvec);

   EXPECT_EQ(1.0, dvec1[0]);
   EXPECT_EQ(2.0, dvec1[1]);
   EXPECT_EQ(3.0, dvec1[2]);
   EXPECT_EQ(3u, dvec1.size());
   EXPECT_EQ(1.0, dvec2[0]);
   EXPECT_EQ(2.0, dvec2[1]);
   EXPECT_EQ(3.0, dvec2[2]);
   EXPECT_EQ(3u, dvec2.size());
}

TEST(VecOps, ArithmeticsUnary)
{
   ROOT::VecOps::RVec<int> ivec{1, 2, 3};
   ROOT::VecOps::RVec<int> pvec = +ivec;
   ROOT::VecOps::RVec<int> nvec = -ivec;
   ROOT::VecOps::RVec<int> tvec = ~ivec;

   EXPECT_EQ(1, pvec[0]);
   EXPECT_EQ(2, pvec[1]);
   EXPECT_EQ(3, pvec[2]);
   EXPECT_EQ(3u, pvec.size());

   EXPECT_EQ(-1, nvec[0]);
   EXPECT_EQ(-2, nvec[1]);
   EXPECT_EQ(-3, nvec[2]);
   EXPECT_EQ(3u, nvec.size());

   EXPECT_EQ(-2, tvec[0]);
   EXPECT_EQ(-3, tvec[1]);
   EXPECT_EQ(-4, tvec[2]);
   EXPECT_EQ(3u, tvec.size());
}

TEST(VecOps, MathScalar)
{
   ROOT::VecOps::RVec<double> ref{1, 2, 3};
   ROOT::VecOps::RVec<double> v(ref);
   int scalar = 3;
   auto plus = v + scalar;
   auto minus = v - scalar;
   auto mult = v * scalar;
   auto div = v / scalar;

   CheckEqual(plus, ref + scalar);
   CheckEqual(minus, ref - scalar);
   CheckEqual(mult, ref * scalar);
   CheckEqual(div, ref / scalar);

   // The same with views
   ROOT::VecOps::RVec<double> w(ref.data(), ref.size());
   plus = w + scalar;
   minus = w - scalar;
   mult = w * scalar;
   div = w / scalar;

   CheckEqual(plus, ref + scalar);
   CheckEqual(minus, ref - scalar);
   CheckEqual(mult, ref * scalar);
   CheckEqual(div, ref / scalar);
}

TEST(VecOps, MathScalarInPlace)
{
   ROOT::VecOps::RVec<double> ref{1, 2, 3};
   const ROOT::VecOps::RVec<double> v(ref);
   int scalar = 3;
   auto plus = v;
   plus += scalar;
   auto minus = v;
   minus -= scalar;
   auto mult = v;
   mult *= scalar;
   auto div = v;
   div /= scalar;

   CheckEqual(plus, ref + scalar);
   CheckEqual(minus, ref - scalar);
   CheckEqual(mult, ref * scalar);
   CheckEqual(div, ref / scalar);
}

TEST(VecOps, MathVector)
{
   ROOT::VecOps::RVec<double> ref{1, 2, 3};
   ROOT::VecOps::RVec<double> vec{3, 4, 5};
   ROOT::VecOps::RVec<double> v(ref);
   auto plus = v + vec;
   auto minus = v - vec;
   auto mult = v * vec;
   auto div = v / vec;

   CheckEqual(plus, ref + vec);
   CheckEqual(minus, ref - vec);
   CheckEqual(mult, ref * vec);
   CheckEqual(div, ref / vec);

   // The same with 1 view
   ROOT::VecOps::RVec<double> w(ref.data(), ref.size());
   plus = w + vec;
   minus = w - vec;
   mult = w * vec;
   div = w / vec;

   CheckEqual(plus, ref + vec);
   CheckEqual(minus, ref - vec);
   CheckEqual(mult, ref * vec);
   CheckEqual(div, ref / vec);

   // The same with 2 views
   ROOT::VecOps::RVec<double> w2(ref.data(), ref.size());
   plus = w + w2;
   minus = w - w2;
   mult = w * w2;
   div = w / w2;

   CheckEqual(plus, ref + w2);
   CheckEqual(minus, ref - w2);
   CheckEqual(mult, ref * w2);
   CheckEqual(div, ref / w2);
}

TEST(VecOps, MathVectorInPlace)
{
   ROOT::VecOps::RVec<double> ref{1, 2, 3};
   ROOT::VecOps::RVec<double> vec{3, 4, 5};
   ROOT::VecOps::RVec<double> v(ref);
   auto plus = v;
   plus += vec;
   auto minus = v;
   minus -= vec;
   auto mult = v;
   mult *= vec;
   auto div = v;
   div /= vec;

   CheckEqual(plus, ref + vec);
   CheckEqual(minus, ref - vec);
   CheckEqual(mult, ref * vec);
   CheckEqual(div, ref / vec);
}

TEST(VecOps, Filter)
{
   ROOT::VecOps::RVec<int> v{0, 1, 2, 3, 4, 5};
   const std::vector<int> vEvenRef{0, 2, 4};
   const std::vector<int> vOddRef{1, 3, 5};
   auto vEven = v[v % 2 == 0];
   auto vOdd = v[v % 2 == 1];
   CheckEqual(vEven, vEvenRef, "Even check");
   CheckEqual(vOdd, vOddRef, "Odd check");

   // now with the helper function
   vEven = Filter(v, [](int i) { return 0 == i % 2; });
   vOdd = Filter(v, [](int i) { return 1 == i % 2; });
   CheckEqual(vEven, vEvenRef, "Even check");
   CheckEqual(vOdd, vOddRef, "Odd check");
}

template <typename T, typename V>
std::string PrintRVec(ROOT::VecOps::RVec<T> v, V w)
{
   std::stringstream ss;
   ss << v << " " << w << std::endl;
   ss << v + w << std::endl;
   ss << v - w << std::endl;
   ss << v * w << std::endl;
   ss << v / w << std::endl;
   ss << (v > w) << std::endl;
   ss << (v >= w) << std::endl;
   ss << (v == w) << std::endl;
   ss << (v <= w) << std::endl;
   ss << (v < w) << std::endl;
   ss << w + v << std::endl;
   ss << w - v << std::endl;
   ss << w * v << std::endl;
   ss << w / v << std::endl;
   ss << (w > v) << std::endl;
   ss << (w >= v) << std::endl;
   ss << (w == v) << std::endl;
   ss << (w <= v) << std::endl;
   ss << (w < v) << std::endl;
   return ss.str();
}

TEST(VecOps, PrintOps)
{
   ROOT::VecOps::RVec<int> ref{1, 2, 3};
   ROOT::VecOps::RVec<int> v(ref);

   auto ref0 = R"ref0({ 1, 2, 3 } 2
{ 3, 4, 5 }
{ -1, 0, 1 }
{ 2, 4, 6 }
{ 0.5, 1, 1.5 }
{ 0, 0, 1 }
{ 0, 1, 1 }
{ 0, 1, 0 }
{ 1, 1, 0 }
{ 1, 0, 0 }
{ 3, 4, 5 }
{ 1, 0, -1 }
{ 2, 4, 6 }
{ 2, 1, 0.666667 }
{ 1, 0, 0 }
{ 1, 1, 0 }
{ 0, 1, 0 }
{ 0, 1, 1 }
{ 0, 0, 1 }
)ref0";
   auto t0 = PrintRVec(v, 2.);
   EXPECT_STREQ(t0.c_str(), ref0);
   auto ref1 = R"ref1({ 1, 2, 3 } { 3, 4, 5 }
{ 4, 6, 8 }
{ -2, -2, -2 }
{ 3, 8, 15 }
{ 0, 0, 0 }
{ 0, 0, 0 }
{ 0, 0, 0 }
{ 0, 0, 0 }
{ 1, 1, 1 }
{ 1, 1, 1 }
{ 4, 6, 8 }
{ 2, 2, 2 }
{ 3, 8, 15 }
{ 3, 2, 1 }
{ 1, 1, 1 }
{ 1, 1, 1 }
{ 0, 0, 0 }
{ 0, 0, 0 }
{ 0, 0, 0 }
)ref1";
   auto t1 = PrintRVec(v, ref + 2);
   EXPECT_STREQ(t1.c_str(), ref1);

   ROOT::VecOps::RVec<int> w(ref.data(), ref.size());

   auto ref2 = R"ref2({ 1, 2, 3 } 2
{ 3, 4, 5 }
{ -1, 0, 1 }
{ 2, 4, 6 }
{ 0.5, 1, 1.5 }
{ 0, 0, 1 }
{ 0, 1, 1 }
{ 0, 1, 0 }
{ 1, 1, 0 }
{ 1, 0, 0 }
{ 3, 4, 5 }
{ 1, 0, -1 }
{ 2, 4, 6 }
{ 2, 1, 0.666667 }
{ 1, 0, 0 }
{ 1, 1, 0 }
{ 0, 1, 0 }
{ 0, 1, 1 }
{ 0, 0, 1 }
)ref2";
   auto t2 = PrintRVec(v, 2.);
   EXPECT_STREQ(t2.c_str(), ref2);

   auto ref3 = R"ref3({ 1, 2, 3 } { 3, 4, 5 }
{ 4, 6, 8 }
{ -2, -2, -2 }
{ 3, 8, 15 }
{ 0, 0, 0 }
{ 0, 0, 0 }
{ 0, 0, 0 }
{ 0, 0, 0 }
{ 1, 1, 1 }
{ 1, 1, 1 }
{ 4, 6, 8 }
{ 2, 2, 2 }
{ 3, 8, 15 }
{ 3, 2, 1 }
{ 1, 1, 1 }
{ 1, 1, 1 }
{ 0, 0, 0 }
{ 0, 0, 0 }
{ 0, 0, 0 }
)ref3";
   auto t3 = PrintRVec(v, ref + 2);
   EXPECT_STREQ(t3.c_str(), ref3);
}

#ifdef R__HAS_VDT
#include <vdt/vdtMath.h>
#endif

TEST(VecOps, MathFuncs)
{
   ROOT::VecOps::RVec<double> u{1, 1, 1};
   ROOT::VecOps::RVec<double> v{1, 2, 3};
   ROOT::VecOps::RVec<double> w{1, 4, 27};

   CheckEqual(pow(1,v), u, " error checking math function pow");
   CheckEqual(pow(v,1), v, " error checking math function pow");
   CheckEqual(pow(v,v), w, " error checking math function pow");

   CheckEqual(sqrt(v), Map(v, [](double x) { return std::sqrt(x); }), " error checking math function sqrt");
   CheckEqual(log(v), Map(v, [](double x) { return std::log(x); }), " error checking math function log");
   CheckEqual(sin(v), Map(v, [](double x) { return std::sin(x); }), " error checking math function sin");
   CheckEqual(cos(v), Map(v, [](double x) { return std::cos(x); }), " error checking math function cos");
   CheckEqual(tan(v), Map(v, [](double x) { return std::tan(x); }), " error checking math function tan");
   CheckEqual(atan(v), Map(v, [](double x) { return std::atan(x); }), " error checking math function atan");
   CheckEqual(sinh(v), Map(v, [](double x) { return std::sinh(x); }), " error checking math function sinh");
   CheckEqual(cosh(v), Map(v, [](double x) { return std::cosh(x); }), " error checking math function cosh");
   CheckEqual(tanh(v), Map(v, [](double x) { return std::tanh(x); }), " error checking math function tanh");
   CheckEqual(asinh(v), Map(v, [](double x) { return std::asinh(x); }), " error checking math function asinh");
   CheckEqual(acosh(v), Map(v, [](double x) { return std::acosh(x); }), " error checking math function acosh");
   v /= 10.;
   CheckEqual(asin(v), Map(v, [](double x) { return std::asin(x); }), " error checking math function asin");
   CheckEqual(acos(v), Map(v, [](double x) { return std::acos(x); }), " error checking math function acos");
   CheckEqual(atanh(v), Map(v, [](double x) { return std::atanh(x); }), " error checking math function atanh");

#ifdef R__HAS_VDT
   #define CHECK_VDT_FUNC(F) \
   CheckEqual(fast_##F(v), Map(v, [](double x) { return vdt::fast_##F(x); }), "error checking vdt function " #F);

   CHECK_VDT_FUNC(exp)
   CHECK_VDT_FUNC(log)

   CHECK_VDT_FUNC(sin)
   CHECK_VDT_FUNC(sin)
   CHECK_VDT_FUNC(cos)
   CHECK_VDT_FUNC(atan)
   CHECK_VDT_FUNC(acos)
   CHECK_VDT_FUNC(atan)
#endif
}

TEST(VecOps, PhysicsSelections)
{
   // We emulate 8 muons
   ROOT::VecOps::RVec<short> mu_charge{1, 1, -1, -1, -1, 1, 1, -1};
   ROOT::VecOps::RVec<float> mu_pt{56.f, 45.f, 32.f, 24.f, 12.f, 8.f, 7.f, 6.2f};
   ROOT::VecOps::RVec<float> mu_eta{3.1f, -.2f, -1.1f, 1.f, 4.1f, 1.6f, 2.4f, -.5f};

   // Pick the pt of the muons with a pt greater than 10, an eta between -2 and 2 and a negative charge
   // or the ones with a pt > 20, outside the eta range -2:2 and with positive charge
   auto goodMuons_pt = mu_pt[(mu_pt > 10.f && abs(mu_eta) <= 2.f && mu_charge == -1) ||
                             (mu_pt > 15.f && abs(mu_eta) > 2.f && mu_charge == 1)];
   ROOT::VecOps::RVec<float> goodMuons_pt_ref = {56.f, 32.f, 24.f};
   CheckEqual(goodMuons_pt, goodMuons_pt_ref, "Muons quality cut");
}

template<typename T0>
void CheckEq(const T0 &v, const T0 &ref)
{
   auto vsize = v.size();
   auto refsize = ref.size();
   EXPECT_EQ(vsize, refsize) << "Sizes are: " << vsize << " " << refsize << std::endl;
   for (auto i : ROOT::TSeqI(vsize)) {
      EXPECT_EQ(v[i], ref[i]) << "RVecs differ" << std::endl;
   }
}

TEST(VecOps, InputOutput)
{
   auto filename = "vecops_inputoutput.root";
   auto treename = "t";

   const ROOT::VecOps::RVec<double> dref {1., 2., 3.};
   const ROOT::VecOps::RVec<float> fref {1.f, 2.f, 3.f};
   const ROOT::VecOps::RVec<UInt_t> uiref {1, 2, 3};
   const ROOT::VecOps::RVec<ULong_t> ulref {1UL, 2UL, 3UL};
   const ROOT::VecOps::RVec<ULong64_t> ullref {1ULL, 2ULL, 3ULL};
   const ROOT::VecOps::RVec<UShort_t> usref {1, 2, 3};
   const ROOT::VecOps::RVec<UChar_t> ucref {1, 2, 3};
   const ROOT::VecOps::RVec<Int_t> iref {1, 2, 3};;
   const ROOT::VecOps::RVec<Long_t> lref {1UL, 2UL, 3UL};;
   const ROOT::VecOps::RVec<Long64_t> llref {1ULL, 2ULL, 3ULL};
   const ROOT::VecOps::RVec<Short_t> sref {1, 2, 3};
   const ROOT::VecOps::RVec<Char_t> cref {1, 2, 3};
   const ROOT::VecOps::RVec<bool> bref {true, false, true};

   {
      auto d = dref;
      auto f = fref;
      auto ui = uiref;
      auto ul = ulref;
      auto ull = ullref;
      auto us = usref;
      auto uc = ucref;
      auto i = iref;
      auto l = lref;
      auto ll = llref;
      auto s = sref;
      auto c = cref;
      auto b = bref;
      TFile file(filename, "RECREATE");
      TTree t(treename, treename);
      t.Branch("d", &d);
      t.Branch("f", &f);
      t.Branch("ui", &ui);
      t.Branch("ul", &ul);
      t.Branch("ull", &ull);
      t.Branch("us", &us);
      t.Branch("uc", &uc);
      t.Branch("i", &i);
      t.Branch("l", &l);
      t.Branch("ll", &ll);
      t.Branch("s", &s);
      t.Branch("c", &c);
      t.Branch("b", &b);
      t.Fill();
      t.Write();
   }

   auto d = new ROOT::VecOps::RVec<double>();
   auto f = new ROOT::VecOps::RVec<float>;
   auto ui = new ROOT::VecOps::RVec<UInt_t>();
   auto ul = new ROOT::VecOps::RVec<ULong_t>();
   auto ull = new ROOT::VecOps::RVec<ULong64_t>();
   auto us = new ROOT::VecOps::RVec<UShort_t>();
   auto uc = new ROOT::VecOps::RVec<UChar_t>();
   auto i = new ROOT::VecOps::RVec<Int_t>();
   auto l = new ROOT::VecOps::RVec<Long_t>();
   auto ll = new ROOT::VecOps::RVec<Long64_t>();
   auto s = new ROOT::VecOps::RVec<Short_t>();
   auto c = new ROOT::VecOps::RVec<Char_t>();
   auto b = new ROOT::VecOps::RVec<bool>();

   TFile file(filename);
   TTree *tp;
   file.GetObject(treename, tp);
   auto &t = *tp;

   t.SetBranchAddress("d", &d);
   t.SetBranchAddress("f", &f);
   t.SetBranchAddress("ui", &ui);
   t.SetBranchAddress("ul", &ul);
   t.SetBranchAddress("ull", &ull);
   t.SetBranchAddress("us", &us);
   t.SetBranchAddress("uc", &uc);
   t.SetBranchAddress("i", &i);
   t.SetBranchAddress("l", &l);
   t.SetBranchAddress("ll", &ll);
   t.SetBranchAddress("s", &s);
   t.SetBranchAddress("c", &c);
   t.SetBranchAddress("b", &b);

   t.GetEntry(0);
   CheckEq(*d, dref);
   CheckEq(*f, fref);
   CheckEq(*d, dref);
   CheckEq(*f, fref);
   CheckEq(*d, dref);
   CheckEq(*f, fref);
   CheckEq(*d, dref);
   CheckEq(*f, fref);
   CheckEq(*d, dref);
   CheckEq(*f, fref);
   CheckEq(*b, bref);

   gSystem->Unlink(filename);

}

TEST(VecOps, SimpleStatOps)
{
   ROOT::VecOps::RVec<double> v0 {};
   ASSERT_DOUBLE_EQ(Sum(v0), 0.);
   ASSERT_DOUBLE_EQ(Mean(v0), 0.);
   ASSERT_DOUBLE_EQ(StdDev(v0), 0.);
   ASSERT_DOUBLE_EQ(Var(v0), 0.);

   ROOT::VecOps::RVec<double> v1 {42.};
   ASSERT_DOUBLE_EQ(Sum(v1), 42.);
   ASSERT_DOUBLE_EQ(Mean(v1), 42.);
   ASSERT_DOUBLE_EQ(Max(v1), 42.);
   ASSERT_DOUBLE_EQ(Min(v1), 42.);
   ASSERT_DOUBLE_EQ(ArgMax(v1), 0);
   ASSERT_DOUBLE_EQ(ArgMin(v1), 0);
   ASSERT_DOUBLE_EQ(StdDev(v1), 0.);
   ASSERT_DOUBLE_EQ(Var(v1), 0.);

   ROOT::VecOps::RVec<double> v2 {1., 2., 3.};
   ASSERT_DOUBLE_EQ(Sum(v2), 6.);
   ASSERT_DOUBLE_EQ(Mean(v2), 2.);
   ASSERT_DOUBLE_EQ(Max(v2), 3.);
   ASSERT_DOUBLE_EQ(Min(v2), 1.);
   ASSERT_DOUBLE_EQ(ArgMax(v2), 2);
   ASSERT_DOUBLE_EQ(ArgMin(v2), 0);
   ASSERT_DOUBLE_EQ(Var(v2), 1.);
   ASSERT_DOUBLE_EQ(StdDev(v2), 1.);

   ROOT::VecOps::RVec<double> v3 {10., 20., 32.};
   ASSERT_DOUBLE_EQ(Sum(v3), 62.);
   ASSERT_DOUBLE_EQ(Mean(v3), 20.666666666666668);
   ASSERT_DOUBLE_EQ(Max(v3), 32.);
   ASSERT_DOUBLE_EQ(Min(v3), 10.);
   ASSERT_DOUBLE_EQ(ArgMax(v3), 2);
   ASSERT_DOUBLE_EQ(ArgMin(v3), 0);
   ASSERT_DOUBLE_EQ(Var(v3), 121.33333333333337);
   ASSERT_DOUBLE_EQ(StdDev(v3), 11.015141094572206);

   ROOT::VecOps::RVec<int> v4 {2, 3, 1};
   ASSERT_DOUBLE_EQ(Sum(v4), 6.);
   ASSERT_DOUBLE_EQ(Mean(v4), 2.);
   ASSERT_DOUBLE_EQ(Max(v4), 3);
   ASSERT_DOUBLE_EQ(Min(v4), 1);
   ASSERT_DOUBLE_EQ(ArgMax(v4), 1);
   ASSERT_DOUBLE_EQ(ArgMin(v4), 2);
   ASSERT_DOUBLE_EQ(Var(v4), 1.);
   ASSERT_DOUBLE_EQ(StdDev(v4), 1.);
}

TEST(VecOps, Any)
{
   ROOT::VecOps::RVec<int> vi {0, 1, 2};
   EXPECT_TRUE(Any(vi));
   vi = {0, 0, 0};
   EXPECT_FALSE(Any(vi));
   vi = {1, 1};
   EXPECT_TRUE(Any(vi));
}

TEST(VecOps, All)
{
   ROOT::VecOps::RVec<int> vi {0, 1, 2};
   EXPECT_FALSE(All(vi));
   vi = {0, 0, 0};
   EXPECT_FALSE(All(vi));
   vi = {1, 1};
   EXPECT_TRUE(All(vi));
}

TEST(VecOps, Argsort)
{
   ROOT::VecOps::RVec<int> v{2, 0, 1};
   using size_type = typename ROOT::VecOps::RVec<int>::size_type;
   auto i = Argsort(v);
   ROOT::VecOps::RVec<size_type> ref{1, 2, 0};
   CheckEqual(i, ref);
}

TEST(VecOps, ArgsortWithComparisonOperator)
{
   ROOT::VecOps::RVec<int> v{2, 0, 1};
   using size_type = typename ROOT::VecOps::RVec<int>::size_type;

   auto i1 = Argsort(v, [](int x, int y){ return x < y; });
   ROOT::VecOps::RVec<size_type> ref1{1, 2, 0};
   CheckEqual(i1, ref1);

   auto i2 = Argsort(v, [](int x, int y){ return x > y; });
   ROOT::VecOps::RVec<size_type> ref2{0, 2, 1};
   CheckEqual(i2, ref2);
}

TEST(VecOps, TakeIndices)
{
   ROOT::VecOps::RVec<int> v0{2, 0, 1};
   ROOT::VecOps::RVec<typename ROOT::VecOps::RVec<int>::size_type> i{1, 2, 0, 0, 0};
   auto v1 = Take(v0, i);
   ROOT::VecOps::RVec<int> ref{0, 1, 2, 2, 2};
   CheckEqual(v1, ref);
}

TEST(VecOps, TakeFirst)
{
   ROOT::VecOps::RVec<int> v0{0, 1, 2};

   auto v1 = Take(v0, 2);
   ROOT::VecOps::RVec<int> ref{0, 1};
   CheckEqual(v1, ref);

   // Corner-case: Take zero entries
   auto v2 = Take(v0, 0);
   ROOT::VecOps::RVec<int> none{};
   CheckEqual(v2, none);
}

TEST(VecOps, TakeLast)
{
   ROOT::VecOps::RVec<int> v0{0, 1, 2};

   auto v1 = Take(v0, -2);
   ROOT::VecOps::RVec<int> ref{1, 2};
   CheckEqual(v1, ref);

   // Corner-case: Take zero entries
   auto v2 = Take(v0, 0);
   ROOT::VecOps::RVec<int> none{};
   CheckEqual(v2, none);
}

TEST(VecOps, Reverse)
{
   ROOT::VecOps::RVec<int> v0{0, 1, 2};

   auto v1 = Reverse(v0);
   ROOT::VecOps::RVec<int> ref{2, 1, 0};
   CheckEqual(v1, ref);

   // Corner-case: Empty vector
   ROOT::VecOps::RVec<int> none{};
   auto v2 = Reverse(none);
   CheckEqual(v2, none);
}

TEST(VecOps, Sort)
{
   ROOT::VecOps::RVec<int> v{2, 0, 1};

   // Sort in ascending order
   auto v1 = Sort(v);
   ROOT::VecOps::RVec<int> ref{0, 1, 2};
   CheckEqual(v1, ref);

   // Corner-case: Empty vector
   ROOT::VecOps::RVec<int> none{};
   auto v2 = Sort(none);
   CheckEqual(v2, none);
}

TEST(VecOps, SortWithComparisonOperator)
{
   ROOT::VecOps::RVec<int> v{2, 0, 1};

   // Sort with comparison operator
   auto v1 = Sort(v, std::greater<int>());
   ROOT::VecOps::RVec<int> ref{2, 1, 0};
   CheckEqual(v1, ref);

   // Corner-case: Empty vector
   ROOT::VecOps::RVec<int> none{};
   auto v2 = Sort(none, std::greater<int>());
   CheckEqual(v2, none);
}

TEST(VecOps, RVecBool)
{
   // std::vector<bool> is special, so RVec<bool> is special
   ROOT::VecOps::RVec<bool> v{true, false};
   auto v2 = v;
   EXPECT_EQ(v[0], true);
   EXPECT_EQ(v[1], false);
   EXPECT_EQ(v.size(), 2u);
   CheckEqual(v2, v);
}

TEST(VecOps, CombinationsTwoVectors)
{
   ROOT::VecOps::RVec<int> v1{1, 2, 3};
   ROOT::VecOps::RVec<int> v2{-4, -5};

   auto idx = Combinations(v1, v2);
   auto c1 = Take(v1, idx[0]);
   auto c2 = Take(v2, idx[1]);
   auto v3 = c1 * c2;

   RVec<int> ref{-4, -5, -8, -10, -12, -15};
   CheckEqual(v3, ref);

   // Test overload with size as input
   auto idx3 = Combinations(v1.size(), v2.size());
   auto c3 = Take(v1, idx3[0]);
   auto c4 = Take(v2, idx3[1]);
   auto v4 = c3 * c4;

   CheckEqual(v4, ref);

   // Corner-case: One collection is empty
   RVec<int> empty_int{};
   auto idx2 = Combinations(v1, empty_int);
   RVec<size_t> empty_size{};
   CheckEqual(idx2[0], empty_size);
   CheckEqual(idx2[1], empty_size);
}

TEST(VecOps, UniqueCombinationsSingleVector)
{
   // Doubles: x + y
   ROOT::VecOps::RVec<int> v1{1, 2, 3};
   auto idx1 = Combinations(v1, 2);
   auto c1 = Take(v1, idx1[0]);
   auto c2 = Take(v1, idx1[1]);
   auto v2 = c1 + c2;
   RVec<int> ref1{
      3,  // 1+2
      4,  // 1+3
      5}; // 2+3
   CheckEqual(v2, ref1);

   // Triples: x * y * z
   ROOT::VecOps::RVec<int> v3{1, 2, 3, 4};
   auto idx2 = Combinations(v3, 3);
   auto c3 = Take(v3, idx2[0]);
   auto c4 = Take(v3, idx2[1]);
   auto c5 = Take(v3, idx2[2]);
   auto v4 = c3 * c4 * c5;
   RVec<int> ref2{
      6,  // 1*2*3
      8,  // 1*2*3
      12, // 1*3*4
      24};// 2*3*4
   CheckEqual(v4, ref2);

   // Corner-case: Single combination
   ROOT::VecOps::RVec<int> v5{1};
   auto idx3 = Combinations(v5, 1);
   EXPECT_EQ(idx3.size(), 1u);
   EXPECT_EQ(idx3[0].size(), 1u);
   EXPECT_EQ(idx3[0][0], 0u);

   // Corner-case: Insert empty vector
   ROOT::VecOps::RVec<int> empty_int{};
   auto idx4 = Combinations(empty_int, 0);
   EXPECT_EQ(idx4.size(), 0u);

   // Corner-case: Request "zero-tuples"
   auto idx5 = Combinations(v1, 0);
   EXPECT_EQ(idx5.size(), 0u);
}

TEST(VecOps, PrintCollOfNonPrintable)
{
   auto code = "class A{};ROOT::VecOps::RVec<A> v(1);v";
   auto ret = gInterpreter->ProcessLine(code);
   EXPECT_TRUE(0 != ret) << "Error in printing an RVec collection of non printable objects.";
}

TEST(VecOps, Nonzero)
{
   ROOT::VecOps::RVec<int> v1{0, 1, 0, 3, 4, 0, 6};
   ROOT::VecOps::RVec<float> v2{0, 1, 0, 3, 4, 0, 6};
   auto v3 = Nonzero(v1);
   auto v4 = Nonzero(v2);
   ROOT::VecOps::RVec<size_t> ref1{1, 3, 4, 6};
   CheckEqual(v3, ref1);
   CheckEqual(v4, ref1);

   auto v5 = v1[v1<2];
   auto v6 = Nonzero(v5);
   ROOT::VecOps::RVec<size_t> ref2{1};
   CheckEqual(v6, ref2);
}

TEST(VecOps, Intersect)
{
   ROOT::VecOps::RVec<int> v1{0, 1, 2, 3};
   ROOT::VecOps::RVec<int> v2{2, 3, 4, 5};
   auto v3 = Intersect(v1, v2);
   ROOT::VecOps::RVec<int> ref1{2, 3};
   CheckEqual(v3, ref1);

   ROOT::VecOps::RVec<int> v4{4, 5, 3, 2};
   auto v5 = Intersect(v1, v4);
   CheckEqual(v5, ref1);

   // v2 already sorted
   auto v6 = Intersect(v1, v2, true);
   CheckEqual(v6, ref1);
}

TEST(VecOps, Where)
{
   // Use two vectors as arguments
   RVec<float> v0{1, 2, 3, 4};
   RVec<float> v1{-1, -2, -3, -4};
   auto v3 = Where(v0 > 1 && v0 < 4, v0, v1);
   RVec<float> ref1{-1, 2, 3, -4};
   CheckEqual(v3, ref1);

   // Broadcast false argument
   auto v4 = Where(v0 == 2 || v0 == 4, v0, -1.0f);
   RVec<float> ref2{-1, 2, -1, 4};
   CheckEqual(v4, ref2);

   // Broadcast true argument
   auto v5 = Where(v0 == 2 || v0 == 4, -1.0f, v1);
   RVec<float> ref3{-1, -1, -3, -1};
   CheckEqual(v5, ref3);

   // Broadcast both arguments
   auto v6 = Where(v0 == 2 || v0 == 4, -1.0f, 1.0f);
   RVec<float> ref4{1, -1, 1, -1};
   CheckEqual(v6, ref4);
}

TEST(VecOps, AtWithFallback)
{
   ROOT::VecOps::RVec<float> v({1.f, 2.f, 3.f});
   EXPECT_FLOAT_EQ(v.at(7, 99.f), 99.f);
}

TEST(VecOps, Concatenate)
{
   RVec<float> rvf {0.f, 1.f, 2.f};
   RVec<int> rvi {7, 8, 9};
   const auto res = Concatenate(rvf, rvi);
   const RVec<float> ref { 0.00000f, 1.00000f, 2.00000f, 7.00000f, 8.00000f, 9.00000f };
   CheckEqual(res, ref);
}

TEST(VecOps, DeltaPhi)
{
   // Two scalars (radians)
   // NOTE: These tests include the checks of the poundary effects
   const float c1 = M_PI;
   EXPECT_EQ(DeltaPhi(0.f, 2.f), 2.f);
   EXPECT_EQ(DeltaPhi(1.f, 0.f), -1.f);
   EXPECT_EQ(DeltaPhi(-0.5f, 0.5f), 1.f);
   EXPECT_EQ(DeltaPhi(0.f, 2.f * c1 - 1.f), -1.f);
   EXPECT_EQ(DeltaPhi(0.f, 4.f * c1 - 1.f), -1.f);
   EXPECT_EQ(DeltaPhi(0.f, -2.f * c1 + 1.f), 1.f);
   EXPECT_EQ(DeltaPhi(0.f, -4.f * c1 + 1.f), 1.f);

   // Two scalars (degrees)
   const float c2 = 180.f;
   EXPECT_EQ(DeltaPhi(0.f, 2.f, c2), 2.f);
   EXPECT_EQ(DeltaPhi(1.f, 0.f, c2), -1.f);
   EXPECT_EQ(DeltaPhi(-0.5f, 0.5f, c2), 1.f);
   EXPECT_EQ(DeltaPhi(0.f, 2.f * c2 - 1.f, c2), -1.f);
   EXPECT_EQ(DeltaPhi(0.f, 4.f * c2 - 1.f, c2), -1.f);
   EXPECT_EQ(DeltaPhi(0.f, -2.f * c2 + 1.f, c2), 1.f);
   EXPECT_EQ(DeltaPhi(0.f, -4.f * c2 + 1.f, c2), 1.f);

   // Two vectors
   RVec<float> v1 = {0.f, 1.f, -0.5f, 0.f, 0.f, 0.f, 0.f};
   RVec<float> v2 = {2.f, 0.f, 0.5f, 2.f * c1 - 1.f, 4.f * c1 - 1.f, -2.f * c1 + 1.f, -4.f * c1 + 1.f};
   auto dphi1 = DeltaPhi(v1, v2);
   RVec<float> r1 = {2.f, -1.f, 1.f, -1.f, -1.f, 1.f, 1.f};
   CheckEqual(dphi1, r1);

   auto dphi2 = DeltaPhi(v2, v1);
   auto r2 = r1 * -1.f;
   CheckEqual(dphi2, r2);

   // Check against TLorentzVector
   for (std::size_t i = 0; i < v1.size(); i++) {
      TLorentzVector p1, p2;
      p1.SetPtEtaPhiM(1.f, 1.f, v1[i], 1.f);
      p2.SetPtEtaPhiM(1.f, 1.f, v2[i], 1.f);
      EXPECT_NEAR(float(p2.DeltaPhi(p1)), dphi1[i], 1e-6);
   }

   // Vector and scalar
   RVec<float> v3 = {0.f, 1.f, c1, c1 + 1.f, 2.f * c1, 2.f * c1 + 1.f, -1.f, -c1, -c1 + 1.f, -2.f * c1, -2.f * c1 + 1.f};
   float v4 = 1.f;
   auto dphi3 = DeltaPhi(v3, v4);
   RVec<float> r3 = {1.f, 0.f, 1.f - c1, c1, 1.f, 0.f, 2.f, -c1 + 1.f, c1, 1.f, 0.f};
   CheckEqual(dphi3, r3);

   auto dphi4 = DeltaPhi(v4, v3);
   auto r4 = -1.f * r3;
   CheckEqual(dphi4, r4);
}

TEST(VecOps, InvariantMass)
{
   // Dummy particle collections
   RVec<double> mass1 = {50,  50,  50,   50,   100};
   RVec<double> pt1 =   {0,   5,   5,    10,   10};
   RVec<double> eta1 =  {0.0, 0.0, -1.0, 0.5,  2.5};
   RVec<double> phi1 =  {0.0, 0.0, 0.0,  -0.5, -2.4};

   RVec<double> mass2 = {40,  40,  40,  40,  30};
   RVec<double> pt2 =   {0,   5,   5,   10,  2};
   RVec<double> eta2 =  {0.0, 0.0, 0.5, 0.4, 1.2};
   RVec<double> phi2 =  {0.0, 0.0, 0.0, 0.5, 2.4};

   // Compute invariant mass of two particle system using both collections
   const auto invMass = InvariantMasses(pt1, eta1, phi1, mass1, pt2, eta2, phi2, mass2);

   for(size_t i=0; i<mass1.size(); i++) {
      TLorentzVector p1, p2;
      p1.SetPtEtaPhiM(pt1[i], eta1[i], phi1[i], mass1[i]);
      p2.SetPtEtaPhiM(pt2[i], eta2[i], phi2[i], mass2[i]);
      // NOTE: The accuracy of the optimized trigonometric functions is relatively
      // low and the test start to fail with an accuracy of 1e-5.
      EXPECT_NEAR((p1 + p2).M(), invMass[i], 1e-4);
   }

   // Compute invariant mass of multiple-particle system using a single collection
   const auto invMass2 = InvariantMass(pt1, eta1, phi1, mass1);

   TLorentzVector p3;
   p3.SetPtEtaPhiM(pt1[0], eta1[0], phi1[0], mass1[0]);
   for(size_t i=1; i<mass1.size(); i++) {
      TLorentzVector p4;
      p4.SetPtEtaPhiM(pt1[i], eta1[i], phi1[i], mass1[i]);
      p3 += p4;
   }

   EXPECT_NEAR(p3.M(), invMass2, 1e-4);

   const auto invMass3 = InvariantMass(pt2, eta2, phi2, mass2);

   TLorentzVector p5;
   p5.SetPtEtaPhiM(pt2[0], eta2[0], phi2[0], mass2[0]);
   for(size_t i=1; i<mass2.size(); i++) {
      TLorentzVector p6;
      p6.SetPtEtaPhiM(pt2[i], eta2[i], phi2[i], mass2[i]);
      p5 += p6;
   }

   EXPECT_NEAR(p5.M(), invMass3, 1e-4);
}

TEST(VecOps, DeltaR)
{
   RVec<double> eta1 =  {0.1, -1.0, -1.0, 0.5,  -2.5};
   RVec<double> eta2 =  {0.0, 0.0, 0.5, 2.4, 1.2};
   RVec<double> phi1 =  {1.0, 5.0, -1.0,  -0.5, -2.4};
   RVec<double> phi2 =  {0.0, 3.0, 6.0, 1.5, 1.4};

   auto dr = DeltaR(eta1, eta2, phi1, phi2);
   auto dr2 = DeltaR(eta2, eta1, phi2, phi1);

   for (std::size_t i = 0; i < eta1.size(); i++) {
      // Check against TLorentzVector
      TLorentzVector p1, p2;
      p1.SetPtEtaPhiM(1.f, eta1[i], phi1[i], 1.f);
      p2.SetPtEtaPhiM(1.f, eta2[i], phi2[i], 1.f);
      auto dr3 = p2.DeltaR(p1);
      EXPECT_NEAR(dr3, dr[i], 1e-6);
      EXPECT_NEAR(dr3, dr2[i], 1e-6);

      // Check scalar implementation
      auto dr4 = DeltaR(eta1[i], eta2[i], phi1[i], phi2[i]);
      EXPECT_NEAR(dr3, dr4, 1e-6);
   }
}

TEST(VecOps, Map)
{
   RVec<float> a({1.f, 2.f, 3.f});
   RVec<float> b({4.f, 5.f, 6.f});
   RVec<float> c({7.f, 8.f, 9.f});

   auto mod = [](float x, int y, double z) { return sqrt(x * x + y * y + z * z); };

   auto res = Map(a, c, c, mod);

   ROOT::VecOps::RVec<double> ref{9.9498743710661994, 11.489125293076057, 13.076696830622021};
   CheckEqual(res, ref);
}

TEST(VecOps, Construct)
{
   RVec<float> pts {15.5f, 34.32f, 12.95f};
   RVec<float> etas {.3f, 2.2f, 1.32f};
   RVec<float> phis {.1f, 3.02f, 2.2f};
   RVec<float> masses {105.65f, 105.65f, 105.65f};
   auto fourVects = Construct<ROOT::Math::PtEtaPhiMVector>(pts, etas, phis, masses);
   const ROOT::Math::PtEtaPhiMVector ref0 {15.5f, .3f, .1f, 105.65f};
   const ROOT::Math::PtEtaPhiMVector ref1 {34.32f, 2.2f, 3.02f, 105.65f};
   const ROOT::Math::PtEtaPhiMVector ref2 {12.95f, 1.32f, 2.2f, 105.65f};
   EXPECT_TRUE(fourVects[0] == ref0);
   EXPECT_TRUE(fourVects[1] == ref1);
   EXPECT_TRUE(fourVects[2] == ref2);
}

bool IsSmall(const RVec<int> &v)
{
   // the first array element is right after the 3 data members of SmallVectorBase
   return reinterpret_cast<std::uintptr_t>(v.begin()) - reinterpret_cast<std::uintptr_t>(&v) ==
          sizeof(void *) + 2 * sizeof(int);
}

// this is a regression test for https://github.com/root-project/root/issues/6796
TEST(VecOps, MemoryAdoptionAndClear)
{
   ROOT::RVec<int> v{1, 2, 3};
   EXPECT_TRUE(IsSmall(v));
   ROOT::RVec<int> v2(v.data(), v.size());
   v2[0] = 0;
   v2.clear();
   EXPECT_TRUE(All(v == RVec<int>{0, 2, 3}));
   v2.push_back(42);
   EXPECT_FALSE(IsSmall(v2)); // currently RVec does not go back to a small state after `clear()`
   EXPECT_TRUE(All(v2 == RVec<int>{42}));
}

// interaction between small buffer optimization and memory adoption
TEST(VecOps, MemoryAdoptionAndSBO)
{
   int *values = new int[3]{1, 2, 3};
   ROOT::RVec<int> v(values, 3);
   auto check = [](const RVec<int> &mv) {
      EXPECT_EQ(mv.size(), 3u);
      EXPECT_EQ(mv[0], 1);
      EXPECT_EQ(mv[1], 2);
      EXPECT_EQ(mv[2], 3);
   };
   check(v);
   EXPECT_FALSE(IsSmall(v));
   ROOT::RVec<int> v2 = std::move(v);
   EXPECT_TRUE(v.empty());
   check(v2);
   // this changes the RVec from memory adoption mode to "long" mode, even if the size is small
   // currently we don't allow going from memory adoption to small buffer mode directly, it could be future optimization
   v2.push_back(4);
   EXPECT_FALSE(IsSmall(v2));
   v2.clear();
   v2.push_back(1);
   v2.push_back(2);
   v2.push_back(3);
   check(v2);
   delete[] values;
   check(v2);
}

struct ThrowingCtor {
   ThrowingCtor() { throw std::runtime_error("This exception should have been caught."); }
};

struct ThrowingMove {
   ThrowingMove() {}
   ThrowingMove(const ThrowingMove &) {}
   ThrowingMove &operator=(const ThrowingMove &) { return *this; }
   ThrowingMove(ThrowingMove &&) { throw std::runtime_error("This exception should have been caught."); }
   ThrowingMove &operator=(const ThrowingMove &&) { return *this; }
};

// RVec does not guarantee exception safety, but we still want to test
// that we don't segfault or otherwise crash if element construction or move throws.
TEST(VecOps, NoExceptionSafety)
{
   EXPECT_NO_THROW(ROOT::RVec<ThrowingCtor>());

   EXPECT_THROW(ROOT::RVec<ThrowingCtor>(1), std::runtime_error);
   EXPECT_THROW(ROOT::RVec<ThrowingCtor>(42), std::runtime_error);
   ROOT::RVec<ThrowingCtor> v1;
   EXPECT_THROW(v1.push_back(ThrowingCtor{}), std::runtime_error);
   ROOT::RVec<ThrowingCtor> v2;
   EXPECT_THROW(v2.emplace_back(ThrowingCtor{}), std::runtime_error);

   ROOT::RVec<ThrowingMove> v3(2);
   ROOT::RVec<ThrowingMove> v4(42);
   EXPECT_THROW(std::swap(v3, v4), std::runtime_error);
   ThrowingMove tm;
   EXPECT_THROW(v3.emplace_back(std::move(tm)), std::runtime_error);

   // now with memory adoption
   ThrowingCtor *p1 = new ThrowingCtor[0];
   ROOT::RVec<ThrowingCtor> v5(p1, 0);
   EXPECT_THROW(v5.push_back(ThrowingCtor{}), std::runtime_error);
   delete[] p1;

   ThrowingMove *p2 = new ThrowingMove[2];
   ROOT::RVec<ThrowingMove> v6(p2, 2);
   EXPECT_THROW(std::swap(v6, v3), std::runtime_error);
   delete[] p2;
}
