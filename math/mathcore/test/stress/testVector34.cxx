// @(#)root/test:$Id$
// Author: Lorenzo Moneta       06/2005
//         Martin Storø Nyfløtt 05/2017

#include "gtest/gtest.h"

#include "StatFunction.h"
#include "VectorTest.h"

using namespace ROOT::Math;

template <int Dim>
struct Vector32TestDimWrapper {
   static constexpr int GetDim() { return Dim; }
};

template <typename T>
class Vector34Test : public testing::Test {
protected:
   const int fNGen = 10000;
   VectorTest<T::GetDim()> fVectorTest;
   typedef SVector<double, T::GetDim()> SV_t;
   std::vector<SV_t> fV1;

   virtual void SetUp()
   {
      fVectorTest.GenData();
      fV1.reserve(fNGen);
   }

   int GetDim() { return T::GetDim(); }

   std::string VecTypeName() { return VecType<SV_t>::name(); }

   std::string VecTypeName32() { return VecType<SV_t>::name32(); }

public:
   Vector34Test() : fVectorTest(fNGen) {}
};

TYPED_TEST_SUITE_P(Vector34Test);

// test of Svector of dim 3 or 4
TYPED_TEST_P(Vector34Test, TestVector34)
{
   double s1 = 0;
   double sref1 = 0;

   std::string name = "SVector<double," + Util::ToString(this->GetDim()) + ">";
   this->fVectorTest.TestCreate(this->fV1);
   EXPECT_TRUE(IsNear(name + " creation", this->fV1.size(), this->fNGen, 1));
   s1 = this->fVectorTest.TestAddition(this->fV1);
   EXPECT_TRUE(IsNear(name + " addition", s1, this->fVectorTest.Sum(), this->GetDim() * 20));
   sref1 = s1;

   // test the io
   double fsize = 0;
   int ir = 0;

   double estSize = this->fNGen * 8 * this->GetDim() + 47000; // add extra bytes (why so much ? ) this is for ngen=10000
   double scale = 0.1 / std::numeric_limits<double>::epsilon();
   fsize = this->fVectorTest.TestWrite(this->fV1, "ROOT::Math::" + name);
   EXPECT_TRUE(IsNear(name + " write", fsize, estSize, scale));
   ir = this->fVectorTest.TestRead(this->fV1);
   EXPECT_TRUE(IsNear(name + " read", ir, 0, 1));
   s1 = this->fVectorTest.TestAddition(this->fV1);
   EXPECT_TRUE(IsNear(name + " after read", s1, sref1, 10));

   // test Double32
   estSize = this->fNGen * 4 * this->GetDim() + 47000; // add extra bytes
   scale = 0.1 / std::numeric_limits<double>::epsilon();

   fsize = this->fVectorTest.TestWrite(this->fV1, this->VecTypeName32());
   EXPECT_TRUE(IsNear(this->VecTypeName() + "_D32 write", fsize, estSize, scale));
   ir = this->fVectorTest.TestRead(this->fV1);
   EXPECT_TRUE(IsNear(this->VecTypeName() + "_D32 read", ir, 0, 1));
   s1 = this->fVectorTest.TestAddition(this->fV1);
   EXPECT_TRUE(IsNear(this->VecTypeName() + "_D32 after read", s1, sref1, 1.E9));
}

REGISTER_TYPED_TEST_SUITE_P(Vector34Test, TestVector34);

typedef testing::Types<Vector32TestDimWrapper<3>, Vector32TestDimWrapper<4>> SMatrixVectorTypes_t;

INSTANTIATE_TYPED_TEST_SUITE_P(SMatrix, Vector34Test, SMatrixVectorTypes_t);
