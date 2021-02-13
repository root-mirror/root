/// \file
/// \ingroup tutorial_fit
/// \notebook -nodraw
/// Multi-Dimensional Parametrisation and Fitting
///
/// \macro_output
/// \macro_code
///
/// \authors Rene Brun, Christian Holm Christensen

#include "Riostream.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TSystem.h"
#include "TBrowser.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMultiDimFit.h"
#include "TVectorD.h"
#include "TMath.h"


//____________________________________________________________________
void makeData(Double_t* x, Double_t& d, Double_t& e)
{
  // Make data points
  Double_t upp[5] = { 10, 10, 10, 10,  1 };
  Double_t low[5] = {  0,  0,  0,  0, .1 };
  for (int i = 0; i < 4; i++)
    x[i] = (upp[i] - low[i]) * gRandom->Rndm() + low[i];

  d = x[0] * TMath::Sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);

  e = gRandom->Gaus(upp[4],low[4]);
}

//____________________________________________________________________
int CompareResults(TMultiDimFit *fit, bool doFit)
{
      //Compare results with reference run


      // the right coefficients (before fit)
   double GoodCoeffsNoFit[] = {
   -4.37056,
   43.1468,
   13.432,
   13.4632,
   13.3964,
   13.328,
   13.3016,
   13.3519,
   4.49724,
   4.63876,
   4.89036,
   -3.69982,
   -3.98618,
   -3.86195,
   4.36054,
   -4.02597,
   4.57037,
   4.69845,
   2.83819,
   -3.48855,
   -3.97612
   };

      // the right coefficients (after fit)
   double GoodCoeffs[] = {
      -4.399,
      43.15,
      13.41,
      13.49,
      13.4,
      13.23,
      13.34,
      13.29,
      4.523,
      4.659,
      4.948,
      -4.026,
      -4.045,
      -3.939,
      4.421,
      -4.006,
      4.626,
      4.378,
      3.516,
      -4.111,
      -3.823,
   };

   // Good Powers
   int GoodPower[] = {
   1,  1,  1,  1,
   2,  1,  1,  1,
   1,  1,  1,  2,
   1,  1,  2,  1,
   1,  2,  1,  1,
   2,  2,  1,  1,
   2,  1,  1,  2,
   2,  1,  2,  1,
   1,  1,  1,  3,
   1,  3,  1,  1,
   1,  1,  5,  1,
   1,  1,  2,  2,
   1,  2,  1,  2,
   1,  2,  2,  1,
   2,  1,  1,  3,
   2,  2,  1,  2,
   2,  1,  3,  1,
   2,  3,  1,  1,
   1,  2,  2,  2,
   2,  1,  2,  2,
   2,  2,  2,  1
   };

   Int_t nc = fit->GetNCoefficients();
   Int_t nv = fit->GetNVariables();
   const Int_t *powers = fit->GetPowers();
   const Int_t *pindex = fit->GetPowerIndex();
   if (nc != 21) return 1;
   const TVectorD *coeffs = fit->GetCoefficients();
   int k = 0;
   for (Int_t i=0;i<nc;i++) {
      if (doFit) {
         if (!TMath::AreEqualRel((*coeffs)[i],GoodCoeffs[i],1e-3)) return 2;
      }
      else {
         if (TMath::Abs((*coeffs)[i] - GoodCoeffsNoFit[i]) > 5e-5) return 2;
      }
      for (Int_t j=0;j<nv;j++) {
         if (powers[pindex[i]*nv+j] != GoodPower[k]) return 3;
         k++;
      }
   }

   // now test the result of the generated function
   gROOT->ProcessLine(".L MDF.C");

   Double_t refMDF = (doFit) ? 43.95 : 43.98;
   // this does not work in CLing since the function is not defined
   //Double_t x[]    = {5,5,5,5};
   //Double_t rMDF   = MDF(x);
   //LM:  need to return the address of the result since it is casted to a long (this should not be in a tutorial !)
   Long_t iret = gROOT->ProcessLine(" Double_t xvalues[] = {5,5,5,5}; double result=MDF(xvalues); &result;");
   Double_t rMDF = * ( (Double_t*)iret);
   //printf("%f\n",rMDF);
   if (TMath::Abs(rMDF -refMDF) > 1e-2) return 4;
   return 0;
}

//____________________________________________________________________
Int_t multidimfit(bool doFit = true)
{

   cout << "*************************************************" << endl;
   cout << "*             Multidimensional Fit              *" << endl;
   cout << "*                                               *" << endl;
   cout << "* By Christian Holm <cholm@nbi.dk> 14/10/00     *" << endl;
   cout << "*************************************************" << endl;
   cout << endl;

   // Initialize global TRannom object.
   gRandom = new TRandom();

   // Open output file
   TFile* output = new TFile("mdf.root", "RECREATE");

   // Global data parameters
   Int_t nVars       = 4;
   Int_t nData       = 500;
   Double_t x[4];

   // make fit object and set parameters on it.
   TMultiDimFit* fit = new TMultiDimFit(nVars, TMultiDimFit::kMonomials,"v");

   Int_t mPowers[]   = { 6 , 6, 6, 6 };
   fit->SetMaxPowers(mPowers);
   fit->SetMaxFunctions(1000);
   fit->SetMaxStudy(1000);
   fit->SetMaxTerms(30);
   fit->SetPowerLimit(1);
   fit->SetMinAngle(10);
   fit->SetMaxAngle(10);
   fit->SetMinRelativeError(.01);

   // variables to hold the temporary input data
   Double_t d;
   Double_t e;

   // Print out the start parameters
   fit->Print("p");

   printf("======================================\n");

   // Create training sample
   Int_t i;
   for (i = 0; i < nData ; i++) {

      // Make some data
      makeData(x,d,e);

      // Add the row to the fit object
      fit->AddRow(x,d,e);
   }

   // Print out the statistics
   fit->Print("s");

   // Book histograms
   fit->MakeHistograms();

   // Find the parameterization
   fit->FindParameterization();

   // Print coefficents
   fit->Print("rc");

   // Get the min and max of variables from the training sample, used
   // for cuts in test sample.
   Double_t *xMax = new Double_t[nVars];
   Double_t *xMin = new Double_t[nVars];
   for (i = 0; i < nVars; i++) {
      xMax[i] = (*fit->GetMaxVariables())(i);
      xMin[i] = (*fit->GetMinVariables())(i);
   }

   nData = fit->GetNCoefficients() * 100;
   Int_t j;

   // Create test sample
   for (i = 0; i < nData ; i++) {
      // Make some data
      makeData(x,d,e);

      for (j = 0; j < nVars; j++)
         if (x[j] < xMin[j] || x[j] > xMax[j])
      break;

      // If we get through the loop above, all variables are in range
      if (j == nVars)
         // Add the row to the fit object
         fit->AddTestRow(x,d,e);
      else
         i--;
   }
   //delete gRandom;

   // Test the parameterizatio and coefficents using the test sample.
   if (doFit)
      fit->Fit("M");

   // Print result
   fit->Print("fc v");

   // Write code to file
   fit->MakeCode();

   // Write histograms to disk, and close file
   output->Write();
   output->Close();
   delete output;

   // Compare results with reference run
   Int_t compare = CompareResults(fit, doFit);
   if (!compare) {
      printf("\nmultidimfit ..............................................  OK\n");
   } else {
      printf("\nmultidimfit ..............................................  fails case %d\n",compare);
   }

   // We're done
   delete fit;
   delete [] xMin;
   delete [] xMax;
   return compare;
}
