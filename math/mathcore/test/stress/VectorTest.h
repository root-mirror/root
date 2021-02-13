#ifndef ROOT_VECTORTEST
#define ROOT_VECTORTEST

#include "TTree.h"
#include "TRandom3.h"
#include "gtest/gtest.h"
#include "TFile.h"

#include "TSystem.h"

#include "VecType.h"
#include "VecOp.h"

#include <vector>
#include <string>

//*******************************************************************************************************************
// GenVector tests
//*******************************************************************************************************************

// trait for getting vector name

template <int Dim>
class VectorTest {

private:
   // global data variables
   std::vector<double> fDataX;
   std::vector<double> fDataY;
   std::vector<double> fDataZ;
   std::vector<double> fDataE;

   int fNGen;
   int fN2Loop;

   double fSum; // total sum of x,y,z,t (for testing first addition)

public:
   VectorTest(int n1, int n2 = 0) : fNGen(n1), fN2Loop(n2) {}

   double Sum() const { return fSum; }

   void GenData()
   {
      // generate for all 4 d data
      TRandom3 r(111); // use a fixed seed to be able to reproduce tests
      fSum = 0;
      for (int i = 0; i < fNGen; ++i) {

         // generate a 4D vector and stores only the interested dimensions
         double phi = r.Rndm() * 3.1415926535897931;
         double eta = r.Uniform(-5., 5.);
         double pt = r.Exp(10.);
         double m = r.Uniform(0, 10.);
         if (i % 50 == 0) m = r.BreitWigner(1., 0.01);
         double E = sqrt(m * m + pt * pt * std::cosh(eta) * std::cosh(eta));

         // fill vectors
         ROOT::Math::PtEtaPhiEVector q(pt, eta, phi, E);
         fDataX.push_back(q.x());
         fDataY.push_back(q.y());
         fSum += q.x() + q.y();
         if (Dim >= 3) {
            fDataZ.push_back(q.z());
            fSum += q.z();
         }
         if (Dim >= 4) {
            fDataE.push_back(q.t());
            fSum += q.t();
         }
      }
      ASSERT_EQ(int(fDataX.size()), fNGen);
      ASSERT_EQ(int(fDataY.size()), fNGen);
      if (Dim >= 3) { ASSERT_EQ(int(fDataZ.size()), fNGen); }
      if (Dim >= 4) { ASSERT_EQ(int(fDataE.size()), fNGen); }
   }

   // gen data for a Ndim matrix or vector
   void GenDataN()
   {
      // generate for all 4 d data
      TRandom3 r(111); // use a fixed seed to be able to reproduce tests
      fSum = 0;
      fDataX.reserve(fNGen * Dim);
      for (int i = 0; i < fNGen * Dim; ++i) {

         // generate random data between [0,1]
         double x = r.Rndm();
         fSum += x;
         fDataX.push_back(x);
      }
   }

   typedef std::vector<double>::const_iterator DataIt_t;

   // test methods
   template <class V>
   void TestCreate(std::vector<V> &dataV)
   {
      DataIt_t x = fDataX.begin();
      DataIt_t y = fDataY.begin();
      DataIt_t z = fDataZ.begin();
      DataIt_t t = fDataE.begin();
      while (x != fDataX.end()) {
         dataV.push_back(VecOp<V, Dim>::Create(x, y, z, t));
         ASSERT_TRUE(int(dataV.size()) <= fNGen);
      }
   }

   template <class V>
   void TestCreateAndSet(std::vector<V> &dataV)
   {
      DataIt_t x = fDataX.begin();
      DataIt_t y = fDataY.begin();
      DataIt_t z = fDataZ.begin();
      DataIt_t t = fDataE.begin();
      while (x != fDataX.end()) {
         V v;
         VecOp<V, Dim>::Set(v, x, y, z, t);
         dataV.push_back(v);
         ASSERT_TRUE(int(dataV.size()) <= fNGen);
      }
   }

   template <class V>
   double TestAddition(const std::vector<V> &dataV)
   {
      V v0;
      for (int i = 0; i < fNGen; ++i) {
         v0 += dataV[i];
      }
      return VecOp<V, Dim>::Add(v0);
   }

   template <class V>
   double TestOperations(const std::vector<V> &dataV)
   {
      double tot = 0;
      for (int i = 0; i < fNGen - 1; ++i) {
         const V &v1 = dataV[i];
         const V &v2 = dataV[i + 1];
         double a = v1.R();
         double b = v2.mag2(); // mag2 is defined for all dimensions;
         double c = 1. / v1.Dot(v2);
         V v3 = c * (v1 / a + v2 / b);
         tot += VecOp<V, Dim>::Add(v3);
      }
      return tot;
   }

   // mantain loop in gen otherwise is proportional to N**@
   template <class V>
   double TestDelta(const std::vector<V> &dataV)
   {
      double tot = 0;
      for (int i = 0; i < fNGen - 1; ++i) {
         const V &v1 = dataV[i];
         const V &v2 = dataV[i + 1];
         tot += VecOp<V, Dim>::Delta(v1, v2);
      }
      return tot;
   }

   template <class V1, class V2>
   void TestConversion(std::vector<V1> &dataV1, std::vector<V2> &dataV2)
   {
      for (int i = 0; i < fNGen; ++i) {
         dataV2.push_back(V2(dataV1[i]));
      }
   }

   // rotation
   template <class V, class R>
   double TestRotation(std::vector<V> &dataV)
   {
      double sum = 0;
      double rotAngle = 1;
      for (unsigned int i = 0; i < fNGen; ++i) {
         V &v1 = dataV[i];
         V v2 = v1;
         v2.Rotate(rotAngle);
         sum += VecOp<V, Dim>::Add(v2);
      }
      return sum;
   }

   template <class V>
   double TestWrite(const std::vector<V> &dataV, std::string typeName = "", bool compress = false)
   {
      std::string fname = VecType<V>::name() + ".root";
      // replace < character with _
      TFile file(fname.c_str(), "RECREATE", "", compress);

      // create tree
      std::string tree_name = "Tree with" + VecType<V>::name();

      TTree tree("VectorTree", tree_name.c_str());

      V *v1 = new V();

      // need to add namespace to full type name
      if (typeName == "") {
         typeName = "ROOT::Math::" + VecType<V>::name();
      }

      TBranch *br = tree.Branch("Vector_branch", typeName.c_str(), &v1);
      if (br == 0) {
         std::cout << "Error creating branch for" << typeName << "\n\t typeid is " << typeid(*v1).name() << std::endl;
         return -1;
      }

      for (int i = 0; i < fNGen; ++i) {
         *v1 = dataV[i];
         tree.Fill();
      }

      tree.Print();

      file.Write();
      file.Close();
      return file.GetSize();
   }

   template <class V>
   int TestRead(std::vector<V> &dataV)
   {
      dataV.clear();
      dataV.reserve(fNGen);

      std::string fname = VecType<V>::name() + ".root";

      TFile f1(fname.c_str());
      if (f1.IsZombie()) {
         std::cout << " Error opening file " << fname << std::endl;
         return -1;
      }

      // create tree
      TTree *tree = dynamic_cast<TTree *>(f1.Get("VectorTree"));
      if (tree == 0) {
         std::cout << " Error reading file " << fname << std::endl;
         return -1;
      }

      V *v1 = 0;

      // cast to void * to avoid a warning
      tree->SetBranchAddress("Vector_branch", (void *)&v1);

      int n = (int)tree->GetEntries();
      if (n != fNGen) {
         std::cout << "wrong tree entries from file" << fname << std::endl;
         return -1;
      }

      for (int i = 0; i < n; ++i) {
         tree->GetEntry(i);
         dataV.push_back(*v1);
      }

      gSystem->Unlink(fname.c_str());

      return 0;
   }

   // test of SVEctor's or SMatrix
   template <class V>
   void TestCreateSV(std::vector<V> &dataV)
   {
      double *x = &fDataX.front();
      double *end = x + fDataX.size();
      // SVector cannot be created from a generic iterator (should be fixed)
      while (x != end) {
         dataV.push_back(V(x, x + Dim));
         assert(int(dataV.size()) <= fNGen);
         x += Dim;
      }
   }

   template <class V>
   double TestAdditionSV(const std::vector<V> &dataV)
   {
      V v0;
      for (int i = 0; i < fNGen; ++i) {
         v0 += dataV[i];
      }
      double tot = 0;
      typedef typename V::const_iterator It_t;
      for (It_t itr = v0.begin(); itr != v0.end(); ++itr) tot += *itr;

      return tot;
   }

   template <class V>
   double TestAdditionTR(const std::vector<V> &dataV)
   {
      V v0;
      for (int i = 0; i < fNGen; ++i) {
         v0 += dataV[i];
      }
      v0.Print();
      return v0.Sum();
   }
};

#endif // ROOT_VECTORTEST
