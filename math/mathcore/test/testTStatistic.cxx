
#include "TRandom3.h"
#include "TStatistic.h"
#include "TList.h"
#include "TStopwatch.h"
#include <iostream>

bool gVerbose = false;

void testTStatistic(Int_t n = 10000)
{

   // Array
   std::vector<Double_t> xx(n);
   std::vector<Double_t> ww(n);

   double true_mean = 10;
   double true_sigma = 1;
   double eps1 = 5*true_sigma/sqrt(double(n));  // 5 times error on the mean
   double eps2 = 5*true_sigma/sqrt(double(2.*n));  // 5 times error on the RMS

   // Gaussian first
   TRandom3 rand3;
   for (Int_t i = 0; i < n; i++) {
      xx[i] = rand3.Gaus(true_mean, true_sigma);
      ww[i] = rand3.Uniform(0.01, 3.00);
   }

   // Find minimum and maximum
   double min = xx[0]; // the minimum value in the array
   double max = xx[0]; // the maximum value in the array
   double eps_min_max = 0.000001; // epsilon to check for GetMin and GetMax functions
   for (Int_t i = 0; i < n; ++i) {
      if (xx[i] < min) {
         min = xx[i];
      } else if (xx[i] > max) {
         max = xx[i];
      }
   }

   TStopwatch stp;
   bool error = false;

   printf("\nTest without using weights :        ");

   TStatistic st0("st0", n, xx.data());
   stp.Stop();
   if (!TMath::AreEqualAbs(st0.GetMean(),true_mean, eps1) )   { Error("TestTStatistic-GetMean","Different value obtained for the unweighted data"); error = true; }
   if (!TMath::AreEqualAbs(st0.GetRMS(),true_sigma, eps2) )   { Error("TestTStatistic-GetRMS","Different value obtained for the unweighted data"); error = true; }
   if (!TMath::AreEqualAbs(st0.GetMin(), min, eps_min_max) )   { Error("TestTStatistic-GetMin","Different value obtained for the unweighted data"); error = true; }
   if (!TMath::AreEqualAbs(st0.GetMax(), max, eps_min_max) )   { Error("TestTStatistic-GetMax","Different value obtained for the unweighted data"); error = true; }

   if (error) printf("Failed\n");
   else printf("OK\n");
   if (error || gVerbose) {
      stp.Print();
      st0.Print();
   }

   // test with TMath
   printf("\nTest using TMath:                   ");
   error = false;
   stp.Start();
   double mean = TMath::Mean(xx.begin(), xx.end() );
   double rms  = TMath::RMS(xx.begin(), xx.end() );
   stp.Stop();
   if (!TMath::AreEqualAbs(mean,true_mean, eps1) )  {  Error("TestTStatistic::TMath::Mean","Different value obtained for the unweighted data"); error = true; }
   if (!TMath::AreEqualAbs(rms,true_sigma, eps2) )  {  Error("TestTStatistic::TMath::RMS","Different value obtained for the unweighted data"); error = true; }

   if (error) printf("Failed\n");
   else printf("OK\n");
   if (error || gVerbose) {
      stp.Print();
      printf("  TMATH         mu =  %.5g +- %.4g \t RMS = %.5g \n",mean, rms/sqrt(double(xx.size()) ), rms);
   }


   printf("\nTest using Weights :                ");
   error = false;
   stp.Start();
   TStatistic st1("st1", n, xx.data(), ww.data());
   stp.Stop();

   if (!TMath::AreEqualAbs(st1.GetMean(),true_mean, eps1) )  {  Error("TestTStatistic-GetMean","Different value obtained for the weighted data"); error = true; }
   if (!TMath::AreEqualAbs(st1.GetRMS(),true_sigma, eps2) )  {  Error("TestTStatistic-GetRMS","Different value obtained for the weighted data"); error = true; }

   if (error) printf("Failed\n");
   else printf("OK\n");
   if (error || gVerbose) {
      stp.Print();
      st1.Print();
   }


   // Incremental test
   printf("\nTest incremental filling :          ");
   error = false;
   TStatistic st2("st2");
   stp.Start();
   for (Int_t i = 0; i < n; i++) {
      st2.Fill(xx[i], ww[i]);
   }
   stp.Stop();

   if (!TMath::AreEqualRel(st1.GetMean(),st2.GetMean(), 1.E-15) )  {  Error("TestTStatistic-GetMean","2 Different values obtained for the weighted data"); error = true; }
   if (!TMath::AreEqualRel(st1.GetRMS(),st2.GetRMS(), 1.E-15) )    { Error("TestTStatistic-GetRMS","2 Different values obtained for the weighted data"); error = true; }

   if (error) printf("Failed\n");
   else printf("OK\n");
   if (error || gVerbose) {
      stp.Print();
      st2.Print();
   }


   // test merge
   int n1 = rand3.Uniform(10,n-10);

   // sort the data to have then two biased samples
   std::sort(xx.begin(), xx.end() );


   printf("\nTest merge :                        ");
   error = false;
   TStatistic sta("sta");
   TStatistic stb("stb");
   stp.Start();
   for (int i = 0; i < n ; ++i) {
      if (i < n1) sta.Fill(xx[i],ww[i] );
      else   stb.Fill(xx[i],ww[i] );
   }

   TList l; l.Add(&stb);
   sta.Merge(&l);
   stp.Stop();

   if (!TMath::AreEqualAbs(sta.GetMean(),true_mean, eps1) )   { Error("TestTStatistic-GetMean","Different value obtained for the merged data"); error = true; }
   if (!TMath::AreEqualAbs(sta.GetRMS(),true_sigma, eps2) )   { Error("TestTStatistic-GetRMS","Different value obtained for the merged data"); error = true; }

   if (error) printf("Failed\n");
   else printf("OK\n");
   if (error || gVerbose) {
      stp.Print();
      sta.Print();
   }

   printf("\nTest sorted data :                  ");
   error = false;
   stp.Start();
   TStatistic st3("st3", n, xx.data(), ww.data());
   stp.Stop();

   if (!TMath::AreEqualAbs(st3.GetMean(), sta.GetMean(), 1.E-10 ) ) {  Error("TestTStatistic-GetMean","Different value obtained for the sorted data");  error = true; }
   if (!TMath::AreEqualAbs(st3.GetRMS(), sta.GetRMS() , 1.E-10 ) )  {  Error("TestTStatistic-GetRMS","Different value obtained for the sorted data");  error = true; }

   if (error) printf("Failed\n");
   else printf("OK\n");
   if (error || gVerbose) {
      stp.Print();
      st3.Print();
   }


   // test with TMath
   printf("\nTest TMath with weights :           ");
   error = false;
   stp.Start();
   double meanw = TMath::Mean(xx.begin(), xx.end(), ww.begin() );
   double rmsw  = TMath::RMS(xx.begin(), xx.end(), ww.begin() );
   double neff = st2.GetW() * st2.GetW() / st2.GetW2();
   stp.Stop();

   if (!TMath::AreEqualAbs(meanw,true_mean, eps1) )  {  Error("TestTStatistic::TMath::Mean","Different value obtained for the weighted data"); error = true; }
   if (!TMath::AreEqualAbs(rmsw,true_sigma, eps2) )  {  Error("TestTStatistic::TMath::RMS","Different value obtained for the weighted data"); error = true; }

   if (error) printf("Failed\n");
   else printf("OK\n");
   if (error || gVerbose) {
      stp.Print();
      printf("  TMATH         mu =  %.5g +- %.4g \t RMS = %.5g \n",meanw, rmsw/sqrt(neff ), rmsw);
   }


}

int main(int argc, char **argv)
{
  // Parse command line arguments
  for (Int_t i=1 ;  i<argc ; i++) {
     std::string arg = argv[i] ;
     if (arg == "-v") {
      gVerbose = true;
     }
     if (arg == "-h") {
        std::cout << "Usage: " << argv[0] << " [-v]\n";
        std::cout << "  where:\n";
        std::cout << "     -v : verbose  mode";
        std::cout << std::endl;
        return -1;
     }
   }

   testTStatistic();

   return 0;

}
