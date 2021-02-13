#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include "Math/VavilovAccurate.h"
#include "Math/VavilovFast.h"
#include "VavilovTest.cxx"

int testVavilov() {

   ROOT::Math::Vavilov *v;
   //bool fast = false;
   int result = 0;
   int result_pdf[2], result_cdf[2], result_quant[2];

   for (int i = 0; i < 2; ++i) {

      if (i) {
        v = new ROOT::Math::VavilovFast (1, 1);
      }
      else {
        v = new ROOT::Math::VavilovAccurate (1, 1);
      }
      result += (result_pdf[i] = ROOT::Math::VavilovTest::PdfTest (*v, std::cout));
      result += (result_cdf[i] = ROOT::Math::VavilovTest::CdfTest (*v, std::cout));
      result += (result_quant[i] = ROOT::Math::VavilovTest::QuantileTest (*v, std::cout));
   }

   for (int i = 0; i < 2; ++i) {
      if (i) {
        std::cout << "\nResults for VavilovFast:\n";
      }
      else {
        std::cout << "\nResults for VavilovAccurate:\n";
      }
      std::cout << "PdfTest:      ";
      if (result_pdf[i] == 0)
        std::cout << "PASS\n";
      else
        std::cout << "FAIL: " << result_pdf[i] << " / 10\n";
      std::cout << "CdfTest:      ";
      if (result_cdf[i] == 0)
        std::cout << "PASS\n";
      else
        std::cout << "FAIL: " << result_cdf[i] << " / 10\n";
      std::cout << "QuantileTest: ";
      if (result_quant[i] == 0)
        std::cout << "PASS\n";
      else
        std::cout << "FAIL: " << result_quant[i] << " / 20\n";
   }
   std::cout << "\n\nOverall:      ";
   if (result == 0)
     std::cout << "PASS\n";
   else
     std::cout << "FAIL: " << result << " / 80\n";

   return result;
}

int main() {
   return testVavilov();
}
