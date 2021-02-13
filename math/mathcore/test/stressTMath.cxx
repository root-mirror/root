#include <iostream>
#include <algorithm>
#include <typeinfo>

#include <TRandom2.h>
#include <TStopwatch.h>
#include <TMath.h>

using std::cout;
using std::endl;

const unsigned int NUMTEST = 500;

// #define DEBUG

template <typename T> T randD() {
   // use default seed to get same sequence
   static TRandom2 r;
   return (T) r.Uniform(-500,500);
}

double Time(TStopwatch & w) {
   //return w.CpuTime();
   return w.RealTime();
}


template <typename T> double stressVector(unsigned int size, const char* type)
{
   std::cout << "Generating random vector of '"
        << type << "' and size " << size << " ..." << std::endl << std::endl;

   double totalTime = 0;
   double totalUnitTime = 0;

   T *vector = new T[size];
   std::generate(vector, &vector[size], randD<T>);

#ifdef DEBUG
   for ( unsigned int i = 0; i < size; ++i )
      std::cout << vector[i] << " " << std::endl;
#endif

   TStopwatch w;
   std::cout.precision(6);

   unsigned int ntest = 3 * NUMTEST;
   w.Start( kTRUE );
   volatile int iresult = 0;
   for ( unsigned int i = 0; i < ntest; ++i )
      iresult = TMath::MinElement(size, vector);
   w.Stop();
   std::cout << "MinMaxElement() \tTotal Time: " << Time(w) << "  (s)\t\t"
        << " Time/call: " << Time(w)/(ntest)*1.E6 << "   (microsec)" << std::endl;
   totalUnitTime += Time(w)/ntest;
   totalTime += Time(w);

   ntest = 3 * NUMTEST;
   w.Start( kTRUE );
   for ( unsigned int i = 0; i < ntest; ++i )
      iresult  = TMath::LocMin(size, vector);
   w.Stop();
   std::cout << "LocMin/Max() \t\tTotal Time: " << Time(w) << "  (s)\t\t"
        << " Time/call: " << Time(w)/(ntest)*1.E6 << "   (microsec)" << std::endl;
   totalUnitTime += Time(w)/ntest;
   totalTime += Time(w);

   ntest = 10 * NUMTEST;
   w.Start( kTRUE );
   volatile double result;
   for ( unsigned int i = 0; i < ntest; ++i )
      result = TMath::Mean(size, vector);
   w.Stop();
   std::cout << "Mean() \t\t\tTotal Time: " << Time(w) << "  (s)\t\t"
        << " Time/call: " << Time(w)/(ntest)*1.E6 << "   (microsec)" << std::endl;
   totalUnitTime += Time(w)/ntest;
   totalTime += Time(w);

   ntest = (unsigned int) ( NUMTEST/2.5 );
   w.Start( kTRUE );
   for ( unsigned int i = 0; i < ntest; ++i )
      result = TMath::Median(size, vector);
   w.Stop();
   std::cout << "Median() \t\tTotal Time: " << Time(w) << "  (s)\t\t"
        << " Time/call: " << Time(w)/(ntest)*1.E6 << "   (microsec)" << std::endl;
   totalUnitTime += Time(w)/ntest;
   totalTime += Time(w);

   ntest = (unsigned int) ( 10 * NUMTEST );
   w.Start( kTRUE );
   for ( unsigned int i = 0; i < ntest; ++i )
      result = TMath::RMS(size, vector);
   w.Stop();
   std::cout << "RMS() \t\t\tTotal Time: " << Time(w) << "  (s)\t\t"
        << " Time/call: " << Time(w)/(ntest)*1.E6 << "   (microsec)" << std::endl;
   totalUnitTime += Time(w)/ntest;
   totalTime += Time(w);

   ntest = (unsigned int) ( NUMTEST/2.5 );
   w.Start( kTRUE );
   for ( unsigned int i = 0; i < ntest; ++i )
      result = TMath::GeomMean(size, vector);
   w.Stop();
   std::cout << "GeomMean() \t\tTotal Time: " << Time(w) << "  (s)\t\t"
        << " Time/call: " << Time(w)/(ntest)*1.E6 << "   (microsec)" << std::endl;
   totalUnitTime += Time(w)/ntest;
   totalTime += Time(w);

   UInt_t * index =new UInt_t[size];
   ntest = NUMTEST/10;
   w.Start( kTRUE );
   for ( unsigned int i = 0; i < ntest; ++i )
      TMath::Sort(size, vector, index, kFALSE);
   w.Stop();
   std::cout << "Sort() \t\t\tTotal Time: " << Time(w) << "  (s)\t\t"
        << " Time/call: " << Time(w)/(ntest)*1.E6 << "   (microsec)" << std::endl;
   totalUnitTime += Time(w)/ntest;
   totalTime += Time(w);

   std::sort(vector, vector + size);
#ifdef DEBUG
   for ( unsigned int i = 0; i < size; ++i )
      std::cout << vector[i] << " " << std::endl;
#endif
   ntest = 20000*NUMTEST;
   w.Start( kTRUE );
   for ( unsigned int i = 0; i < ntest; ++i )
      iresult = TMath::BinarySearch(size, vector, vector[ i % size ]);
   w.Stop();
   std::cout << "BinarySearch() \t\tTotal Time: " << Time(w) << "  (s)\t\t"
        << " Time/call: " << Time(w)/(ntest)*1.E6 << "   (microsec)" << std::endl;
   totalUnitTime += Time(w)/ntest;
   totalTime += Time(w);

   std::cout << "\nTotal Time :       "      << totalTime     << "  (s)\n"
        <<   "Total Time/call :  " << totalUnitTime*1.E3 << "  (ms)\n" << std::endl;

   // to avoid the warning for un-used variables
   (void)result;
   (void)iresult;

   delete [] vector;
   delete [] index;

   return totalUnitTime;
}

void stressTMath(unsigned int size, const char * type)
{
   double totalTime = 0;

   std::cout << "Stress Test Start..." << std::endl;

   if ( strcmp(type, "Short_t") == 0 )
      totalTime += stressVector<Short_t>(size, type);
   else if ( strcmp(type, "Int_t") == 0 )
      totalTime += stressVector<Int_t>(size, type);
   else if ( strcmp(type, "Float_t") == 0 )
      totalTime += stressVector<Float_t>(size, type);
   else if ( strcmp(type, "Long_t") == 0 )
      totalTime += stressVector<Long_t>(size, type);
   else if ( strcmp(type, "Long64_t") == 0 )
      totalTime += stressVector<Long64_t>(size, type);
   else
      totalTime += stressVector<Double_t>(size, "Double_t");

   //std::cout << "Total Test Time: " << totalTime << "\n" << std::endl;

   std::cout << "End of Stress Test..." << std::endl;

   return;
}


int main(int argc, char* argv[])
{
   // Default size and data type
   unsigned int size = 100000;
   const char *  type = "Double_t";

   if ( argc > 1 ) {
      if (strcmp(argv[1], "-h") == 0) {
         std::cout << "Usage: " << argv[0]
              << " [TYPE OF ARRAY] [SIZE OF ARRAY]\n\n"
              << "where [TYPE OF ARRAY] is one of the following:\n"
              << "\t\"Short_t\"\n"
              << "\t\"Int_t\"\n"
              << "\t\"Float_t\"\n"
              << "\t\"Long_t\"\n"
              << "\t\"Long64_t\"\n"
              << "\t \"Double_t\"\n"
              << std::endl;
         return 1;
      }
      type = argv[1];
   }


   if ( argc > 2 )
      size = (unsigned int) atoi(argv[2]);

   stressTMath(size, type);

   return 0;
}
