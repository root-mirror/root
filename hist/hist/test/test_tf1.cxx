#include "TF1.h"
#include "TF1NormSum.h"
#include "TObjString.h"
#include "TObjArray.h"

#include "gtest/gtest.h"

#include <iostream>

using namespace std;

class MyClass {
public:
   double operator()(double *x, double *p) { return *x + *p; };
};

class MyClassConst {
public:
   double operator()(const double *x, const double *p) { return *x + *p; };
};

double func(const double *x, const double *p)
{
   return *x + *p;
}
double functionConst(const double *x, const double *p)
{
   return *x + *p;
}

Float_t delta = 0.00000000001;

void coeffNamesGeneric(TString &formula, TObjArray *coeffNames)
{
   TF1 cn0("cn0", formula, 0, 1);
   ASSERT_EQ(cn0.GetNpar(), coeffNames->GetEntries());
   for (int i = 0; i < coeffNames->GetEntries(); i++) {
      TObjString *coeffObj = (TObjString *)coeffNames->At(i);
      TString coeffName = coeffObj->GetString();
      EXPECT_EQ(coeffName, TString(cn0.GetParName(i)));
   }
}

// Test that the NSUM names are copied correctly
void test_nsumCoeffNames()
{
   TObjArray *coeffNames = new TObjArray();
   coeffNames->SetOwner(kTRUE);
   coeffNames->Add(new TObjString("sg"));
   coeffNames->Add(new TObjString("bg"));
   coeffNames->Add(new TObjString("Mean"));
   coeffNames->Add(new TObjString("Sigma"));
   coeffNames->Add(new TObjString("Slope"));
   TString formula("NSUM([sg] * gaus, [bg] * expo)");

   coeffNamesGeneric(formula, coeffNames);

   delete coeffNames;
}

// Test that the NSUM is normalized as we'd expect
void test_normalization() {
   double xmin = -5;
   double xmax = 5;
   TF1 n0("n0", "NSUM(.5 * gaus, .5 * (x+[0])**2)", xmin, xmax);
   EXPECT_NEAR(n0.Integral(xmin, xmax), 1, delta);
   n0.SetParameter(4, 1); // should not affect integral
   EXPECT_NEAR(n0.Integral(xmin, xmax), 1, delta);
   n0.SetParameter(0, 0);
   EXPECT_NEAR(n0.Integral(xmin, xmax), .5, delta);

   TF1 n1("n1", "NSUM([sg] * gaus, [bg] * (x+[0])**2)", xmin, xmax);
   n1.SetParameter(0, .5);
   n1.SetParameter(1, .5);
   EXPECT_NEAR(n1.Integral(xmin, xmax), 1, delta);
   n0.SetParameter(0, 0);
   EXPECT_NEAR(n0.Integral(xmin, xmax), .5, delta);

   TF1 n2("n2", "NSUM([sg] * gaus, -0.5 * (x+[0])**2)", xmin, xmax);
   n2.SetParameter(0, .5);
   EXPECT_NEAR(n2.GetParameter(1), -.5, delta);
   EXPECT_NEAR(n2.Integral(xmin, xmax), 0, delta);
   n2.SetParameter(0, 0);
   EXPECT_NEAR(n2.Integral(xmin, xmax), -.5, delta);
}

void voigtHelper(double sigma, double lg)
{
   TF1 lor("lor", "breitwigner", -20, 20);
   lor.SetParameters(1, 0, lg);
   TF1 mygausn("mygausn", "gausn", -20, 20);
   mygausn.SetParameters(1, 0, sigma);

   TF1 conv("conv", "CONV(lor, mygausn)", -20, 20);

   // Voigt should just be the convolution of the gaussian and lorentzian
   TF1 myvoigt("myvoigt", "TMath::Voigt(x, [0], [1])", -20, 20);
   myvoigt.SetParameters(sigma, lg);

   for (double x = -19.5; x < 20; x += .5)
      EXPECT_NEAR(conv.Eval(x), myvoigt.Eval(x), .01 * conv.Eval(x));
}

// Test that we can change the range of TF1NormSum and TF1Convolution
void test_setRange()
{
   // Define TF1 using NSUM
   TF1 f1("f1", "NSUM([sg] * gaus, [bg] * expo)", 0, 1);
   f1.SetParameters(1, 1, 0, 1, 1);
   EXPECT_NEAR(f1.Integral(0, 1), 2, delta);

   // set range and check value of integral again
   f1.SetRange(-5, 5);
   EXPECT_NEAR(f1.Integral(-5, 5), 2, delta);

   // now same thing with CONV
   TF1 f2("f2", "CONV(gaus, breitwigner)", 0, 1);
   f2.SetParameters(1, 1, 1, 1, 1, 1);
   f2.SetRange(-100, 100); // making our convolution much more accurate
   EXPECT_NEAR(f2.Integral(-20, 20), 2.466, .005);
}

// Test that we can copy and clone TF1 objects based on NSUM and CONV
void test_copyClone()
{
   // Define original TF1 using NSUM
   TF1 *f1 = new TF1("f1", "NSUM(gaus, expo)", -5, 5);
   f1->SetParameters(1, 1, 0, 1, 1);

   // Make copy and test
   TF1 *f2 = new TF1();
   f1->Copy(*f2);
   for (double x = -5; x < 5; x += .5)
      EXPECT_NEAR(f1->Eval(x), f2->Eval(x), delta);

   // Make clone and test
   TF1 *f3 = (TF1 *)f1->Clone();
   for (double x = -5; x < 5; x += .5)
      EXPECT_NEAR(f1->Eval(x), f3->Eval(x), delta);

   delete f1;
   delete f2;
   delete f3;

   // Same thing for CONV

   // Define original TF1 using NSUM
   TF1 *f4 = new TF1("f4", "CONV(gaus, breitwigner)", -15, 15);
   f4->SetParameters(1, 1, 1, 1, 1, 1);

   // Make copy
   TF1 *f5 = new TF1();
   f4->Copy(*f5);
   for (double x = -5; x < 5; x += .5)
      EXPECT_NEAR(f4->Eval(x), f5->Eval(x), delta);

   // Make clone
   TF1 *f6 = (TF1 *)f4->Clone();
   for (double x = -5; x < 5; x += .5)
      EXPECT_NEAR(f4->Eval(x), f6->Eval(x), delta);

   delete f4;
   delete f5;
   delete f6;
}

// Test that the voigt can be expressed as a convolution of a gaussian and lorentzian
// Check that the values match to within 1%
void test_convVoigt()
{
   voigtHelper(.1, 1);
   voigtHelper(1, .1);
   voigtHelper(1, 1);
}

TEST(TF1, NsumCoeffNames)
{
   test_nsumCoeffNames();
}

TEST(TF1, Normalization)
{
   test_normalization();
}

TEST(TF1, ConvVoigt)
{
   test_convVoigt();
}

TEST(TF1, SetRange)
{
   test_setRange();
}

TEST(TF1, CopyClone)
{
   test_copyClone();
}

TEST(TF1, Constructors)
{

   MyClass testObj, testObjConst;
   double x = 1;
   double p = 1;

   std::function<double(const double *data, const double *param)> stdFunction = functionConst;
   std::vector<TF1> vtf1;
   vtf1.emplace_back(TF1("Functor", testObj, 0, 1, 0));
   vtf1.emplace_back(TF1("FunctorConst", testObjConst, 0, 1, 0));
   vtf1.emplace_back(TF1("FunctorPtr", &testObj, 0, 1, 0));
   vtf1.emplace_back(TF1("FunctorConstPtr", &testObjConst, 0, 1, 0));
   vtf1.emplace_back(TF1("function", func));
   vtf1.emplace_back(TF1("functionConst", functionConst));
   vtf1.emplace_back(TF1("lambda", [](double *x1, double *x2) { return *x1 + *x2; }));
   vtf1.emplace_back(TF1("lambdaConst", [](double *x1, double *x2) { return *x1 + *x2; }));
   vtf1.emplace_back(TF1("stdFunction", stdFunction));

   for (auto tf1 : vtf1)
      EXPECT_EQ(tf1(&x, &p), 2);
}
