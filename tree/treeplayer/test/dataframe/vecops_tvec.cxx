#include <gtest/gtest.h>
#include <ROOT/TVec.hxx>
#include <vector>
#include <sstream>

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
   ROOT::Experimental::VecOps::TVec<int> v;
   EXPECT_EQ(v.size(), 0u);
}

TEST(VecOps, InitListCtor)
{
   ROOT::Experimental::VecOps::TVec<int> v{1, 2, 3};
   EXPECT_EQ(v.size(), 3u);
   EXPECT_EQ(v[0], 1);
   EXPECT_EQ(v[1], 2);
   EXPECT_EQ(v[2], 3);
}

TEST(VecOps, CopyCtor)
{
   ROOT::Experimental::VecOps::TVec<int> v1{1, 2, 3};
   ROOT::Experimental::VecOps::TVec<int> v2(v1);
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

TEST(VecOps, CopyCtorCheckNoLeak)
{
   ROOT::Experimental::VecOps::TVec<TLeakChecker> ref;
   ref.emplace_back(TLeakChecker());
   ROOT::Experimental::VecOps::TVec<TLeakChecker> proxy(ref.data(), ref.size());
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
   ROOT::Experimental::VecOps::TVec<int> v1{1, 2, 3};
   ROOT::Experimental::VecOps::TVec<int> v2(std::move(v1));
   EXPECT_EQ(v1.size(), 0u);
   EXPECT_EQ(v2.size(), 3u);
}

TEST(VecOps, MathScalar)
{
   ROOT::Experimental::VecOps::TVec<double> ref{1, 2, 3};
   ROOT::Experimental::VecOps::TVec<double> v(ref);
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
   ROOT::Experimental::VecOps::TVec<double> w(ref.data(), ref.size());
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
   ROOT::Experimental::VecOps::TVec<double> ref{1, 2, 3};
   const ROOT::Experimental::VecOps::TVec<double> v(ref);
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
   ROOT::Experimental::VecOps::TVec<double> ref{1, 2, 3};
   ROOT::Experimental::VecOps::TVec<double> vec{3, 4, 5};
   ROOT::Experimental::VecOps::TVec<double> v(ref);
   auto plus = v + vec;
   auto minus = v - vec;
   auto mult = v * vec;
   auto div = v / vec;

   CheckEqual(plus, ref + vec);
   CheckEqual(minus, ref - vec);
   CheckEqual(mult, ref * vec);
   CheckEqual(div, ref / vec);

   // The same with 1 view
   ROOT::Experimental::VecOps::TVec<double> w(ref.data(), ref.size());
   plus = w + vec;
   minus = w - vec;
   mult = w * vec;
   div = w / vec;

   CheckEqual(plus, ref + vec);
   CheckEqual(minus, ref - vec);
   CheckEqual(mult, ref * vec);
   CheckEqual(div, ref / vec);

   // The same with 2 views
   ROOT::Experimental::VecOps::TVec<double> w2(ref.data(), ref.size());
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
   ROOT::Experimental::VecOps::TVec<double> ref{1, 2, 3};
   ROOT::Experimental::VecOps::TVec<double> vec{3, 4, 5};
   ROOT::Experimental::VecOps::TVec<double> v(ref);
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
   ROOT::Experimental::VecOps::TVec<int> v{0, 1, 2, 3, 4, 5};
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
std::string PrintTVec(ROOT::Experimental::VecOps::TVec<T> v, V w)
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
   ROOT::Experimental::VecOps::TVec<int> ref{1, 2, 3};
   ROOT::Experimental::VecOps::TVec<int> v(ref);

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
   auto t0 = PrintTVec(v, 2.);
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
   auto t1 = PrintTVec(v, ref + 2);
   EXPECT_STREQ(t1.c_str(), ref1);

   ROOT::Experimental::VecOps::TVec<int> w(ref.data(), ref.size());

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
   auto t2 = PrintTVec(v, 2.);
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
   auto t3 = PrintTVec(v, ref + 2);
   EXPECT_STREQ(t3.c_str(), ref3);
}

TEST(VecOps, MathFuncs)
{
   ROOT::Experimental::VecOps::TVec<double> v{1, 2, 3};
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
}

TEST(VecOps, PhysicsSelections)
{
   // We emulate 8 muons
   ROOT::Experimental::VecOps::TVec<short> mu_charge{1, 1, -1, -1, -1, 1, 1, -1};
   ROOT::Experimental::VecOps::TVec<float> mu_pt{56.f, 45.f, 32.f, 24.f, 12.f, 8.f, 7.f, 6.2f};
   ROOT::Experimental::VecOps::TVec<float> mu_eta{3.1f, -.2f, -1.1f, 1.f, 4.1f, 1.6f, 2.4f, -.5f};

   // Pick the pt of the muons with a pt greater than 10, an eta between -2 and 2 and a negative charge
   // or the ones with a pt > 20, outside the eta range -2:2 and with positive charge
   auto goodMuons_pt = mu_pt[(mu_pt > 10.f && abs(mu_eta) <= 2.f && mu_charge == -1) ||
                             (mu_pt > 15.f && abs(mu_eta) > 2.f && mu_charge == 1)];
   ROOT::Experimental::VecOps::TVec<float> goodMuons_pt_ref = {56.f, 32.f, 24.f};
   CheckEqual(goodMuons_pt, goodMuons_pt_ref, "Muons quality cut");
}
