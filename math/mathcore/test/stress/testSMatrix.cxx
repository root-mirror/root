// @(#)root/test:$Id$
// Author: Lorenzo Moneta       06/2005
//         Martin Storø Nyfløtt 05/2017

#include "Math/Util.h"
#include "Math/SMatrix.h"
#include "TestHelper.h"
#include "gtest/gtest.h"

#include "VectorTest.h"
#include "TROOT.h"
#include "TSystem.h"
#include "Rep.h"

using namespace ROOT::Math;

template <int D1, int D2, typename Rep>
struct GenericSMatrixTypeWrapper {
   static constexpr int GetD1() { return D1; }
   static constexpr int GetD2() { return D2; }
   typedef Rep Rep_t;
};

template <typename T>
class SMatrixTest : public testing::Test {
protected:
   const int fNGen = 10000;

   typedef typename T::Rep_t::R_t R_t;
   typedef SMatrix<double, T::GetD1(), T::GetD2(), R_t> SM_t;
   std::vector<SM_t> v1;
   // in the case of sym matrices SM::kSize is different than R::kSize
   // need to use the R::kSize for Dim
   const int fDim = R_t::kSize;

   VectorTest<R_t::kSize> fVectorTest;

   virtual void SetUp()
   {
      fVectorTest.GenDataN();
      v1.reserve(fNGen);
   }

   std::string GetRepName() { return T::Rep_t::Name(); }

   std::string GetRepSName() { return T::Rep_t::SName(); }

   std::string GetRepName32() { return T::Rep_t::Name32(); }

   std::string GetD1AsString() { return Util::ToString(T::GetD1()); }

   std::string GetD2AsString() { return Util::ToString(T::GetD2()); }

   bool ShouldTestDouble32_t()
   {
      if (T::GetD1() != T::GetD2()) return false;
      if (T::GetD1() < 3 || T::GetD1() > 6) return false;
      return true;
   }

public:
   SMatrixTest() : fVectorTest(fNGen) {}
};

TYPED_TEST_SUITE_P(SMatrixTest);

// test of generic SMatrix
TYPED_TEST_P(SMatrixTest, TestSMatrix)
{
   double s1 = 0;
   double sref1 = 0;

   std::string name0 = "SMatrix<double," + this->GetD1AsString() + "," + this->GetD2AsString();
   std::string name = name0 + "," + this->GetRepSName() + ">";

   this->fVectorTest.TestCreateSV(this->v1);
   EXPECT_TRUE(IsNear(name + " creation", this->v1.size(), this->fNGen, 1));
   s1 = this->fVectorTest.TestAdditionSV(this->v1);
   EXPECT_TRUE(IsNear(name + " addition", s1, this->fVectorTest.Sum(), this->fDim * 4));
   sref1 = s1;

   // test the io
   double fsize = 0;
   int ir = 0;

   // the full name is needed for sym matrices
   std::string typeName = "ROOT::Math::" + name0 + "," + this->GetRepName() + ">";

   double estSize = this->fNGen * 8 * this->fDim + 10000;
   double scale = 0.1 / std::numeric_limits<double>::epsilon();
   fsize = this->fVectorTest.TestWrite(this->v1, typeName);
   EXPECT_TRUE(IsNear(name + " write", fsize, estSize, scale));
   ir = this->fVectorTest.TestRead(this->v1);
   EXPECT_TRUE(IsNear(name + " read", ir, 0, 1));
   s1 = this->fVectorTest.TestAdditionSV(this->v1);
   EXPECT_TRUE(IsNear(name + " after read", s1, sref1, 1));

   // test storing as Double32_t
   // dictionary exist only for square matrices between 3 and 6
   if (this->ShouldTestDouble32_t()) {
      double fsize32 = 0;

      name0 = "SMatrix<Double32_t," + this->GetD1AsString() + "," + this->GetD2AsString();
      name = name0 + "," + this->GetRepSName() + ">";
      typeName = "ROOT::Math::" + name0 + "," + this->GetRepName32() + ">";

      estSize = this->fNGen * 4 * this->fDim + 60158;
      scale = 0.1 / std::numeric_limits<double>::epsilon();
      fsize32 = this->fVectorTest.TestWrite(this->v1, typeName);
      EXPECT_TRUE(IsNear(name + " write", fsize32, estSize, scale));
      ir = this->fVectorTest.TestRead(this->v1);
      EXPECT_TRUE(IsNear(name + " read", ir, 0, 1));
      // we read back float (scale errors then)
      s1 = this->fVectorTest.TestAdditionSV(this->v1);
      EXPECT_TRUE(IsNear(name + " after read", s1, sref1, 1.E9));
   }
}

REGISTER_TYPED_TEST_SUITE_P(SMatrixTest, TestSMatrix);

typedef testing::Types<GenericSMatrixTypeWrapper<3, 4, RepStd<3, 4>>, GenericSMatrixTypeWrapper<4, 3, RepStd<4, 3>>,
                       GenericSMatrixTypeWrapper<3, 3, RepStd<3, 3>>,
                       GenericSMatrixTypeWrapper<5, 5, RepSym<5>>> // sym matrix
   SMatrixTestingTypes_t;

INSTANTIATE_TYPED_TEST_SUITE_P(SMatrix, SMatrixTest, SMatrixTestingTypes_t);
