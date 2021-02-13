// @(#)root/test:$Id$
// Author: Lorenzo Moneta       06/2005
//         Martin Storø Nyfløtt 05/2017

#include "gtest/gtest.h"

#include "StatFunction.h"
#include "TestHelper.h"
#include "VectorTest.h"

template <class T>
class VectorTestFixture : public testing::Test {
protected:
   const int fNGen = 10000;
   const int fDim = T::kSize;
   std::string fName = T::Type();
   VectorTest<T::kSize> fVectorTest;
   std::vector<T> fV1;

   virtual void SetUp()
   {
      fVectorTest.GenDataN();
      fV1.reserve(fNGen);
   }

   bool IsD32() { return T::IsD32(); }

public:
   VectorTestFixture() : fVectorTest(fNGen) {}
};

TYPED_TEST_SUITE_P(VectorTestFixture);

// Test of a Composite Object (containing Vector's and Matrices)
TYPED_TEST_P(VectorTestFixture, TestCompositeObj)
{
   double s1 = 0;
   double sref1 = 0;

   this->fVectorTest.TestCreateSV(this->fV1);
   EXPECT_TRUE(IsNear(this->fName + " creation", this->fV1.size(), this->fNGen, 1));

   s1 = this->fVectorTest.TestAdditionTR(this->fV1);
   EXPECT_TRUE(IsNear(this->fName + " addition", s1, this->fVectorTest.Sum(), this->fDim * 4));
   sref1 = s1;

   double fsize = 0;
   int ir = 0;

   // the full name is needed for sym matrices
   auto typeName = this->fName;

   int wsize = 8;
   if (this->IsD32()) wsize = 4;

   double estSize = this->fNGen * wsize * this->fDim + 10000;
   double scale = 0.2 / std::numeric_limits<double>::epsilon();

   fsize = this->fVectorTest.TestWrite(this->fV1, typeName);
   EXPECT_TRUE(IsNear(this->fName + " write", fsize, estSize, scale));

   ir = this->fVectorTest.TestRead(this->fV1);
   EXPECT_TRUE(IsNear(this->fName + " read", ir, 0, 1));

   scale = 1;
   if (this->IsD32()) scale = 1.E9; // use float numbers

   s1 = this->fVectorTest.TestAdditionTR(this->fV1);
   EXPECT_TRUE(IsNear(this->fName + " after read", s1, sref1, scale));

   if (this->IsD32()) {
      // check at double precision type must fail otherwise Double's are stored
      ASSERT_FALSE(IsNear("Double32 test", s1, sref1));
   }
}

REGISTER_TYPED_TEST_SUITE_P(VectorTestFixture, TestCompositeObj);

typedef testing::Types<TrackD, TrackD32, TrackErrD, TrackErrD32, VecTrack<TrackD>, VecTrack<TrackErrD>> TrackTypes_t;

INSTANTIATE_TYPED_TEST_SUITE_P(StressMathCore, VectorTestFixture, TrackTypes_t);
