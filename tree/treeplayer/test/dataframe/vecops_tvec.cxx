#include <gtest/gtest.h>
#include <ROOT/TVec.hxx>
#include <vector>

TEST(VecOps, DefaultCtor)
{
   ROOT::Experimental::VecOps::TVec<int> v;
   EXPECT_EQ(v.size(), 0u);
}

TEST(VecOps, InitListCtor)
{
   ROOT::Experimental::VecOps::TVec<int> v{1,2,3};
   EXPECT_EQ(v.size(), 3u);
   EXPECT_EQ(v[0], 1);
   EXPECT_EQ(v[1], 2);
   EXPECT_EQ(v[2], 3);
}

TEST(VecOps, CopyCtor)
{
   ROOT::Experimental::VecOps::TVec<int> v1{1,2,3};
   ROOT::Experimental::VecOps::TVec<int> v2(v1);
   EXPECT_EQ(v1.size(), 3u);
   EXPECT_EQ(v2.size(), 3u);
   EXPECT_EQ(v2[0], 1);
   EXPECT_EQ(v2[1], 2);
   EXPECT_EQ(v2[2], 3);
}

TEST(VecOps, MoveCtor)
{
   ROOT::Experimental::VecOps::TVec<int> v1{1,2,3};
   ROOT::Experimental::VecOps::TVec<int> v2(std::move(v1));
   EXPECT_EQ(v1.size(), 0u);
   EXPECT_EQ(v2.size(), 3u);
}

TEST(VecOps, MathScalar)
{
   using namespace ROOT::Experimental::VecOps;
   TVec<double> ref {1,2,3};
   TVec<double> v (ref);
   int scalar = 3;
   auto plus = v + scalar;
   auto minus = v - scalar;
   auto mult = v * scalar;
   auto div = v / scalar;

   for (unsigned int i = 0; i< v.size(); ++i) {
      EXPECT_EQ(plus[i], ref[i] + scalar);
      EXPECT_EQ(minus[i], ref[i] - scalar);
      EXPECT_EQ(mult[i], ref[i] * scalar);
      EXPECT_EQ(div[i], ref[i] / scalar);
   }

   // The same with views
   ROOT::Experimental::VecOps::TVec<double> w(ref.data(), ref.size());
   plus = w + scalar;
   minus = w - scalar;
   mult = w * scalar;
   div = w / scalar;

   for (unsigned int i = 0; i< v.size(); ++i) {
      EXPECT_EQ(plus[i], ref[i] + scalar);
      EXPECT_EQ(minus[i], ref[i] - scalar);
      EXPECT_EQ(mult[i], ref[i] * scalar);
      EXPECT_EQ(div[i], ref[i] / scalar);
   }
}

TEST(VecOps, MathVector)
{
   using namespace ROOT::Experimental::VecOps;
   TVec<double> ref {1,2,3};
   TVec<double> vec {3,4,5};
   TVec<double> v (ref);
   auto plus = v + vec;
   auto minus = v - vec;
   auto mult = v * vec;
   auto div = v / vec;

   for (unsigned int i = 0; i< v.size(); ++i) {
      EXPECT_EQ(plus[i], ref[i] + vec[i]);
      EXPECT_EQ(minus[i], ref[i] - vec[i]);
      EXPECT_EQ(mult[i], ref[i] * vec[i]);
      EXPECT_EQ(div[i], ref[i] / vec[i]);
   }

   // The same with 1 view
   ROOT::Experimental::VecOps::TVec<double> w(ref.data(), ref.size());
   plus = w + vec;
   minus = w - vec;
   mult = w * vec;
   div = w / vec;

   for (unsigned int i = 0; i< v.size(); ++i) {
      EXPECT_EQ(plus[i], ref[i] + vec[i]);
      EXPECT_EQ(minus[i], ref[i] - vec[i]);
      EXPECT_EQ(mult[i], ref[i] * vec[i]);
      EXPECT_EQ(div[i], ref[i] / vec[i]);
   }

   // The same with 2 views
   ROOT::Experimental::VecOps::TVec<double> w2(ref.data(), ref.size());
   plus = w + w2;
   minus = w - w2;
   mult = w * w2;
   div = w / w2;

   for (unsigned int i = 0; i< v.size(); ++i) {
      EXPECT_EQ(plus[i], ref[i] + w2[i]);
      EXPECT_EQ(minus[i], ref[i] - w2[i]);
      EXPECT_EQ(mult[i], ref[i] * w2[i]);
      EXPECT_EQ(div[i], ref[i] / w2[i]);
   }
}

TEST(VecOps, Filter)
{
   using namespace ROOT::Experimental::VecOps;
   TVec<int> v {0,1,2,3,4,5};
   std::vector<int> vEvenRef {0,2,4};
   std::vector<int> vOddRef {1,3,5};
   auto vEven = v[ v % 2 == 0];
   auto vOdd = v[ v % 2 == 1];
   for (int i = 0 ; i < 3; ++i) {
      EXPECT_EQ(vEven[i], vEvenRef[i]);
      EXPECT_EQ(vOdd[i], vOddRef[i]);
   }

}

template<typename T, typename V>
void CompAndPrintTVec(ROOT::Experimental::VecOps::TVec<T> v, V w)
{
   using namespace ROOT::Experimental::VecOps;
   std::cout << v << " " << w << std::endl;
   std::cout << v + w << std::endl;
   std::cout << v - w << std::endl;
   std::cout << v * w << std::endl;
   std::cout << v / w << std::endl;
   std::cout << (v > w) << std::endl;
   std::cout << (v >= w) << std::endl;
   std::cout << (v == w) << std::endl;
   std::cout << (v <= w) << std::endl;
   std::cout << (v < w) << std::endl;

}

TEST(VecOps, PrintOps)
{
   using namespace ROOT::Experimental::VecOps;
   TVec<int> ref {1,2,3};
   TVec<int> v (ref);
   CompAndPrintTVec(v, 2.);
   CompAndPrintTVec(v, ref+2);

   ROOT::Experimental::VecOps::TVec<int> w(ref.data(), ref.size());

   CompAndPrintTVec(v, 2.);
   CompAndPrintTVec(v, ref+2);

}

void CheckEq(std::string_view name, const ROOT::Experimental::VecOps::TVec<double> v, const ROOT::Experimental::VecOps::TVec<double> w)
{
   auto wsize = w.size();
   EXPECT_EQ(wsize, v.size());
   for (unsigned int i = 0 ; i < wsize; i++) {
      EXPECT_DOUBLE_EQ(v[i], w[i]) << " error while checking math function " << name;
   }
}

TEST(VecOps, MathFuncs)
{
   using namespace ROOT::Experimental::VecOps;
   TVec<double> v {1,2,3};
   CheckEq("sqrt", sqrt(v), Map(v, [](double x){ return std::sqrt(x);}));
   CheckEq("log",log(v), Map(v, [](double x){ return std::log(x);}));
   CheckEq("sin",sin(v), Map(v, [](double x){ return std::sin(x);}));
   CheckEq("cos",cos(v), Map(v, [](double x){ return std::cos(x);}));
   CheckEq("tan",tan(v), Map(v, [](double x){ return std::tan(x);}));
   CheckEq("atan",atan(v), Map(v, [](double x){ return std::atan(x);}));
   CheckEq("sinh",sinh(v), Map(v, [](double x){ return std::sinh(x);}));
   CheckEq("cosh",cosh(v), Map(v, [](double x){ return std::cosh(x);}));
   CheckEq("tanh",tanh(v), Map(v, [](double x){ return std::tanh(x);}));
   CheckEq("asinh",asinh(v), Map(v, [](double x){ return std::asinh(x);}));
   CheckEq("acosh",acosh(v), Map(v, [](double x){ return std::acosh(x);}));
   v = v / 10.;
   CheckEq("asin",asin(v), Map(v, [](double x){ return std::asin(x);}));
   CheckEq("acos",acos(v), Map(v, [](double x){ return std::acos(x);}));
   CheckEq("atanh",atanh(v), Map(v, [](double x){ return std::atanh(x);}));
}
