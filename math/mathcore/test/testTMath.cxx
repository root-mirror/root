#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <typeinfo>

#include <TMath.h>

using namespace std;
using namespace TMath;

bool showVector = true;

template <typename T>
void testNormCross()
{
   //Float_t NormCross(const Float_t v1[3],const Float_t v2[3],Float_t out[3])

   T fv1[] = {1,2,3};
   T fv2[] = {4,5,6};
   T fout[3];

   std::string type;
   if ( strcmp( typeid(T).name(), "f" ) == 0 )
      type = "Float_t";
   else if ( strcmp( typeid(T).name(), "d" ) == 0 )
      type = "Double_t";
   else
      type = typeid(T).name();

   TMath::NormCross(fv1,fv2,fout);

   cout << "NormCross(const " << type << " v1[3],const "
   << type << " v2[3]," << type << " out[3]): out = ["
   << fout[0] << ", " << fout[1] << ", " << fout[2] << "]"
   << endl;
}


template <typename T, typename U>
void testArrayFunctions()
{
   const U n = 10;
   const U k = 3;
   U index[n];
   U is;

   T sa[n] = { 2, 55 ,23, 57, -9, 24, 6, 82, -4, 10};

   if ( showVector )
   {
      cout << "Vector a[] = {" << sa[0];
      for ( Int_t i = 1; i < n; ++i )
         cout << ", " << sa[i];
      cout << "}\n" << endl;
      showVector = false;
   }

   cout << "Min: a[" << LocMin(n, sa) << "] = " << MinElement(n, sa)
        << " Max: a[" << LocMax(n, sa) << "] = " << MaxElement(n, sa)
        << " Mean: " << Mean(n, sa)
        << " GeomMean: " << GeomMean(n, sa)
        << " RMS: " << RMS(n, sa)
        << " Median: " << Median(n, sa)
        << " KOrdStat(3): " << KOrdStat(n, sa, k)
        << endl;

   Sort(n, sa, index, kFALSE);
   cout << "Sorted a[] = {" << sa[index[0]];
   for ( Int_t i = 1; i < n; ++i )
      cout << ", " << sa[index[i]];
   cout << "}" << endl;

   sort(sa, sa+n);
   is = BinarySearch(n, sa, (T) 57);
   cout << "BinarySearch(n, a, 57) = " << is << "\n" << endl;
}

template <typename T>
void testIteratorFunctions()
{
   const Long64_t n = 10;
   vector<Int_t> index(n);
   Long64_t is;

   T tsa[n] = { 2, 55 ,23, 57, -9, 24, 6, 82, -4, 10};
   vector<T> sa(n);
   for ( int i = 0;  i < n; ++i ) sa[i] = tsa[i];

   if ( showVector )
   {
      cout << "\nVector a[] = {" << sa[0];
      for ( Int_t i = 1; i < n; ++i )
         cout << ", " << sa[i];
      cout << "}\n" << endl;
      showVector = false;
   }

   cout << "Min: " << *LocMin(sa.begin(), sa.end())
        << " Max: " << *LocMax(sa.begin(), sa.end())
        << " Mean: " << Mean(sa.begin(), sa.end())
        << " GeomMean: " << GeomMean(sa.begin(), sa.end())
        << " RMS: " << RMS(sa.begin(), sa.end())
        << endl;

   TMath::SortItr(sa.begin(), sa.end(), index.begin(), kFALSE);
   cout << "Sorted a[] = {" << sa[ index[0] ];
   for ( Int_t i = 1; i < n; ++i )
      cout << ", " << sa[ index[i] ];
   cout << "}" << endl;

   sort(&sa[0], &sa[0]+n);
   is = BinarySearch(n, &sa[0], (T) 57);
   cout << "BinarySearch(n, a, 57) = " << is << "\n" << endl;
}

template <typename T>
void testPoints(T x, T y)
{
   const Int_t n = 4;

   T dx[4] = {0, 0, 2, 2};
   T dy[4] = {0, 2, 2, 0};
   cout << "Point(" << x << "," << y << ") IsInside?: "
        << IsInside( x, y, n, dx, dy) << endl;

}

template <typename T>
void testPlane()
{
   T dp1[3] = {0,0,0};
   T dp2[3] = {1,0,0};
   T dp3[3] = {0,1,0};
   T dn[3];
   Normal2Plane(dp1, dp2, dp3, dn);
   cout << "Normal: ("
        << dn[0] << ", "
        << dn[1] << ", "
        << dn[2] << ")"
        << endl;
}

void testTMath()
{
   cout << "Starting tests on TMath..." << endl;

   cout << "\nNormCross tests: " << endl;

   testNormCross<Float_t>();
   testNormCross<Double_t>();

   cout << "\nArray functions tests: " << endl;

   testArrayFunctions<Short_t,Long64_t>();
   testArrayFunctions<Int_t,Long64_t>();
   testArrayFunctions<Float_t,Long64_t>();
   testArrayFunctions<Double_t,Long64_t>();
   testArrayFunctions<Double_t,Int_t>();
   testArrayFunctions<Long_t,Long64_t>();
   testArrayFunctions<Long64_t,Long64_t>();

   cout << "\nIterator functions tests: " << endl;

   testIteratorFunctions<Short_t>();
   testIteratorFunctions<Int_t>();
   testIteratorFunctions<Float_t>();
   testIteratorFunctions<Double_t>();
   testIteratorFunctions<Long_t>();
   testIteratorFunctions<Long64_t>();

   cout << "\nPoint functions tests: " << endl;

   testPoints<Double_t>(1.3, 0.5);
   testPoints<Float_t>(-0.2, 1.7);
   testPoints<Int_t>(1, 1);

   cout << "\nPLane functions tests: " << endl;

   testPlane<Double_t>();
   testPlane<Float_t>();
}

int main()
{
   testTMath();

   return 0;
}
