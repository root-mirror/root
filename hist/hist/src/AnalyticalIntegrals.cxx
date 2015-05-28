//
//  AnalyticalIntegrals.cxx
//  
//
//  Created by Aurélie Flandi on 09.09.14.
//
//

#include <stdio.h>
#include "TROOT.h"
#include "TMath.h"
#include "AnalyticalIntegrals.h"
#include "Math/DistFuncMathCore.h" //for cdf


using namespace std;

Double_t AnalyticalIntegral(TF1 *f, Double_t a, Double_t b)
{

   Double_t xmin = a;
   Double_t xmax = b;
   Int_t    num  = f->GetNumber();
   Double_t *p   = f->GetParameters();
   Double_t result = 0.;

   TFormula * formula = f->GetFormula();
   if (!formula) {
      Error("TF1::AnalyticalIntegral","Invalid formula number - return a NaN"); 
      result = TMath::QuietNaN();
   }

   if   (num == 200)//expo: exp(p0+p1*x)
   {
      result = (exp(p[0])/p[1])*(exp(-p[1]*xmin)-exp(-p[1]*xmax));
   }
   else if (num == 100)//gaus: [0]*exp(-0.5*((x-[1])/[2])^2))
   {
      double amp   = p[0];
      double mean  = p[1];
      double sigma = p[2];
      if (formula->TestBit(TFormula::kNormalized) ) 
         result =  amp*( ROOT::Math::gaussian_cdf(xmax, sigma, mean)- ROOT::Math::gaussian_cdf(xmin, sigma, mean) );
      else
         result =  amp*sqrt(2*M_PI)*sigma*(ROOT::Math::gaussian_cdf(xmax, sigma, mean)- ROOT::Math::gaussian_cdf(xmin, sigma, mean));//
   }
   else if (num == 400)//landau: root::math::landau(x,mpv=0,sigma=1,bool norm=false)
   {

      double amp   = p[0];
      double mean  = p[1];
      double sigma = p[2];
      //printf("computing integral for landau in [%f,%f] for m=%f s = %f \n",xmin,xmax,mean,sigma);
      if (formula->TestBit(TFormula::kNormalized) ) 
         result = amp*(ROOT::Math::landau_cdf(xmax,sigma,mean) - ROOT::Math::landau_cdf(xmin,sigma,mean));
      else 
         result = amp*sigma*(ROOT::Math::landau_cdf(xmax,sigma,mean) - ROOT::Math::landau_cdf(xmin,sigma,mean));
   }
   else if (num >= 300 && num < 400)//polN
   {
      Int_t n = num - 300;
      for (int i=0;i<n+1;i++)
      {
         result += p[i]/(i+1)*(std::pow(xmax,i+1)-std::pow(xmin,i+1));
      }
   }
   else
      result = TMath::QuietNaN();
   
   return result;
}
