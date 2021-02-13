// @(#)root/test:$Id$
// Author: Lorenzo Moneta       06/2005
//         Martin Storø Nyfløtt 05/2017

#include "gtest/gtest.h"

#include "StatFunction.h"
#include "TestHelper.h"
#include "VectorTest.h"
#include "Rep.h"

using namespace ROOT::Math;

template <typename V1, typename V2, int Dim>
struct VectorTestTypeWrapper {
   typedef V1 V1_t;
   typedef V2 V2_t;
   static constexpr int GetDim() { return Dim; }
};

template <class T>
class GenVectorTest : public testing::Test {
protected:
   const int fNGen = 10000;
   VectorTest<T::GetDim()> fVectorTest;
   std::vector<typename T::V1_t> fV1;
   std::vector<typename T::V2_t> fV2;
   const int fDim = T::GetDim();

   virtual void SetUp()
   {
      fVectorTest.GenData();
      fV1.reserve(fNGen);
      fV2.reserve(fNGen);
   }

   std::string V1Name() { return VecType<typename T::V1_t>::name(); }

   std::string V2Name() { return VecType<typename T::V2_t>::name(); }

   std::string V1Name32() { return VecType<typename T::V1_t>::name32(); }

public:
   GenVectorTest() : fVectorTest(fNGen) {}
};

TYPED_TEST_SUITE_P(GenVectorTest);

// Test of Physics Vector (GenVector package)
TYPED_TEST_P(GenVectorTest, TestGenVectors)
{
   double s1, s2 = 0;
   double scale = 1;
   double sref1, sref2 = 0;

   this->fVectorTest.TestCreate(this->fV1);
   EXPECT_TRUE(IsNear(this->V1Name() + " creation", this->fV1.size(), this->fNGen, 1));
   s1 = this->fVectorTest.TestAddition(this->fV1);
   EXPECT_TRUE(IsNear(this->V1Name() + " addition", s1, this->fVectorTest.Sum(), this->fDim * 20));
   sref1 = s1;
   this->fV1.clear();
   EXPECT_TRUE(this->fV1.size() == 0);
   this->fVectorTest.TestCreateAndSet(this->fV1);
   EXPECT_TRUE(IsNear(this->V1Name() + " creation", this->fV1.size(), this->fNGen, 1));
   s2 = this->fVectorTest.TestAddition(this->fV1);
   EXPECT_TRUE(IsNear(this->V1Name() + " setting", s2, s1, 1));

   this->fVectorTest.TestConversion(this->fV1, this->fV2);
   EXPECT_TRUE(IsNear(this->V1Name() + " -> " + this->V2Name(), this->fV1.size(), this->fV2.size(), 1));

   scale = 1000;
   if (this->fDim == 2) scale = 1.E12; // to be understood
   if (this->fDim == 3) scale = 1.E4;  // problem with RhoEtaPhiVector
   s2 = this->fVectorTest.TestAddition(this->fV2);
   EXPECT_TRUE(IsNear("Vector conversion", s2, s1, scale));

   sref2 = s2;
   s1 = this->fVectorTest.TestOperations(this->fV1);
   scale = this->fDim * 20;
   if (this->fDim == 3 && this->V2Name() == "RhoEtaPhiVector") scale *= 12; // for problem with RhoEtaPhi
   if (this->fDim == 4 && ( this->V2Name() == "PtEtaPhiMVector"  || this->V2Name() == "PxPyPzMVector")) {
#if (defined(__arm__) || defined(__arm64__) || defined(__aarch64__))
      scale *= 1.E7;
#else
      scale *= 10;
#endif
#if defined(__FAST_MATH__) && defined(__clang__)
      scale *= 1.E6;
#endif  
   }

#if defined(R__LINUX) && !defined(R__B64)
   // problem of precision on linux 32
   if (this->fDim == 4) scale = 1000000000;
#endif
   // for problem with PtEtaPhiE
   if (this->fDim == 4 && this->V2Name() == "PtEtaPhiEVector") scale = 0.01 / (std::numeric_limits<double>::epsilon());
   s2 = this->fVectorTest.TestOperations(this->fV2);
   EXPECT_TRUE(IsNear(this->V2Name() + " operations", s2, s1, scale));

   s1 = this->fVectorTest.TestDelta(this->fV1);

   scale = this->fDim * 16;
   if (this->fDim == 4) scale *= 100; // for problem with PtEtaPhiE
   s2 = this->fVectorTest.TestDelta(this->fV2);
   EXPECT_TRUE(IsNear(this->V2Name() + " delta values", s2, s1, scale));

   double fsize = 0;
   int ir = 0;

   double estSize = this->fNGen * 8 * this->fDim + 10000; // add extra bytes
   scale = 0.1 / std::numeric_limits<double>::epsilon();
   fsize = this->fVectorTest.TestWrite(this->fV1);
   EXPECT_TRUE(IsNear(this->V1Name() + " write", fsize, estSize, scale));
   ir = this->fVectorTest.TestRead(this->fV1);
   EXPECT_TRUE(IsNear(this->V1Name() + " read", ir, 0));
   s1 = this->fVectorTest.TestAddition(this->fV1);
   EXPECT_TRUE(IsNear(this->V1Name() + " after read", s1, sref1));

   // test io vector 2
   fsize = this->fVectorTest.TestWrite(this->fV2);
   EXPECT_TRUE(IsNear(this->V2Name() + " write", fsize, estSize, scale));
   ir = this->fVectorTest.TestRead(this->fV2);
   EXPECT_TRUE(IsNear(this->V2Name() + " read", ir, 0));
   scale = 100; // gcc 4.3.2 gives and error for RhoEtaPhiVector for 32 bits
   s2 = this->fVectorTest.TestAddition(this->fV2);
   EXPECT_TRUE(IsNear(this->V2Name() + " after read", s2, sref2, scale));

   // test io of double 32 for vector 1
   if (this->fDim != 2) {
      estSize = this->fNGen * 4 * this->fDim + 10000; // add extra bytes
      scale = 0.1 / std::numeric_limits<double>::epsilon();

      fsize = this->fVectorTest.TestWrite(this->fV1, this->V1Name32());
      EXPECT_TRUE(IsNear(this->V1Name() + "_D32 write", fsize, estSize, scale));
      ir = this->fVectorTest.TestRead(this->fV1);
      EXPECT_TRUE(IsNear(this->V1Name() + "_D32 read", ir, 0));
      s1 = this->fVectorTest.TestAddition(this->fV1);
      EXPECT_TRUE(IsNear(this->V1Name() + "_D32 after read", s1, sref1, 1.E9));
   }
}

REGISTER_TYPED_TEST_SUITE_P(GenVectorTest, TestGenVectors);

typedef testing::Types<
   VectorTestTypeWrapper<XYVector, Polar2DVector, 2>, VectorTestTypeWrapper<XYZVector, Polar3DVector, 3>,
   VectorTestTypeWrapper<XYZVector, RhoEtaPhiVector, 3>, VectorTestTypeWrapper<XYZVector, RhoZPhiVector, 3>,
   VectorTestTypeWrapper<XYZTVector, PtEtaPhiEVector, 4>, VectorTestTypeWrapper<XYZTVector, PtEtaPhiMVector, 4>,
   VectorTestTypeWrapper<XYZTVector, PxPyPzMVector, 4>>
   VectorGenTypes_t;

INSTANTIATE_TYPED_TEST_SUITE_P(StressMathCore, GenVectorTest, VectorGenTypes_t);
