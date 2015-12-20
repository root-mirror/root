// @(#)root/roostats:$Id$

/*************************************************************************
 * Project: RooStats                                                     *
 * Package: RooFit/RooStats                                              *
 * Authors:                                                              *
 *   Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke       *
 *************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////



#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h" // for RooFit::Extended()
#include "RooNLLVar.h"
#include "RooRealVar.h"
#include "RooAbsData.h"

#include "RooStats/HybridResult.h"
#include "RooStats/HybridPlot.h"


/// ClassImp for building the THtml documentation of the class 
using namespace std;

ClassImp(RooStats::HybridResult)

using namespace RooStats;

///////////////////////////////////////////////////////////////////////////

HybridResult::HybridResult( const char *name) :
   HypoTestResult(name),
   fTestStat_data(-999.),
   fComputationsNulDoneFlag(false),
   fComputationsAltDoneFlag(false),
   fSumLargerValues(false)
{
   // HybridResult default constructor (with name )
}

///////////////////////////////////////////////////////////////////////////

HybridResult::HybridResult( const char *name,
                            const std::vector<double>& testStat_sb_vals,
                            const std::vector<double>& testStat_b_vals,
			    bool sumLargerValues ) :
   HypoTestResult(name,0,0),
   fTestStat_data(-999.),
   fComputationsNulDoneFlag(false),
   fComputationsAltDoneFlag(false),
   fSumLargerValues(sumLargerValues)
{
   // HybridResult constructor (with name, title and vectors of S+B and B values)

   int vector_size_sb = testStat_sb_vals.size();
   assert(vector_size_sb>0);

   int vector_size_b = testStat_b_vals.size();
   assert(vector_size_b>0);

   fTestStat_sb.reserve(vector_size_sb);
   for (int i=0;i<vector_size_sb;++i)
      fTestStat_sb.push_back(testStat_sb_vals[i]);

   fTestStat_b.reserve(vector_size_b);
   for (int i=0;i<vector_size_b;++i)
      fTestStat_b.push_back(testStat_b_vals[i]);
}


///////////////////////////////////////////////////////////////////////////

HybridResult::~HybridResult()
{
   // HybridResult destructor

   fTestStat_sb.clear();
   fTestStat_b.clear();
}

///////////////////////////////////////////////////////////////////////////

void HybridResult::SetDataTestStatistics(double testStat_data_val)
{
   // set the value of the test statistics on data

   fComputationsAltDoneFlag = false;
   fComputationsNulDoneFlag = false;
   fTestStat_data = testStat_data_val;
   return;
}

///////////////////////////////////////////////////////////////////////////

double HybridResult::NullPValue() const
{
   // return 1-CL_b : the B p-value

   if (fComputationsNulDoneFlag==false) {
      int nToys = fTestStat_b.size();
      if (nToys==0) {
         std::cout << "Error: no toy data present. Returning -1.\n";
         return -1;
      }

      double larger_than_measured=0;
      if (fSumLargerValues) {
	for (int iToy=0;iToy<nToys;++iToy)
	  if ( fTestStat_b[iToy] >= fTestStat_data ) ++larger_than_measured;
      } else {
	for (int iToy=0;iToy<nToys;++iToy)
	  if ( fTestStat_b[iToy] <= fTestStat_data ) ++larger_than_measured;
      }

      if (larger_than_measured==0) std::cout << "Warning: CLb = 0 ... maybe more toys are needed!\n";

      fComputationsNulDoneFlag = true;
      fNullPValue = 1-larger_than_measured/nToys;
   }

   return fNullPValue;
}

///////////////////////////////////////////////////////////////////////////

double HybridResult::AlternatePValue() const
{
   // return CL_s+b : the S+B p-value

   if (fComputationsAltDoneFlag==false) {
      int nToys = fTestStat_b.size();
      if (nToys==0) {
         std::cout << "Error: no toy data present. Returning -1.\n";
         return -1;
      }

      double larger_than_measured=0;
      if (fSumLargerValues) {
	for (int iToy=0;iToy<nToys;++iToy)
	  if ( fTestStat_sb[iToy] >= fTestStat_data ) ++larger_than_measured;
      } else {
	for (int iToy=0;iToy<nToys;++iToy)
	  if ( fTestStat_sb[iToy] <= fTestStat_data ) ++larger_than_measured;
      }

      if (larger_than_measured==0) std::cout << "Warning: CLsb = 0 ... maybe more toys are needed!\n";

      fComputationsAltDoneFlag = true;
      fAlternatePValue = larger_than_measured/nToys;
   }

   return fAlternatePValue;
}

///////////////////////////////////////////////////////////////////////////

Double_t HybridResult::CLbError() const
{
  // Returns an estimate of the error on CLb assuming a binomial error on
  // CLb:
  // BEGIN_LATEX
  // #sigma_{CL_{b}} &=& #sqrt{CL_{b} #left( 1 - CL_{b} #right) / n_{toys}}
  // END_LATEX
  unsigned const int n = fTestStat_b.size();
  return TMath::Sqrt(CLb() * (1. - CLb()) / n);
}

///////////////////////////////////////////////////////////////////////////

Double_t HybridResult::CLsplusbError() const
{
  // Returns an estimate of the error on CLsplusb assuming a binomial
  // error on CLsplusb:
  // BEGIN_LATEX
  // #sigma_{CL_{s+b}} &=& #sqrt{CL_{s+b} #left( 1 - CL_{s+b} #right) / n_{toys}}
  // END_LATEX
  unsigned const int n = fTestStat_sb.size();
  return TMath::Sqrt(CLsplusb() * (1. - CLsplusb()) / n);
}

///////////////////////////////////////////////////////////////////////////

Double_t HybridResult::CLsError() const
{
  // Returns an estimate of the error on CLs through combination of the
  // errors on CLb and CLsplusb:
  // BEGIN_LATEX
  // #sigma_{CL_s} &=& CL_s #sqrt{#left( #frac{#sigma_{CL_{s+b}}}{CL_{s+b}} #right)^2 + #left( #frac{#sigma_{CL_{b}}}{CL_{b}} #right)^2}
  // END_LATEX
  unsigned const int n_b = fTestStat_b.size();
  unsigned const int n_sb = fTestStat_sb.size();
  
  if (CLb() == 0 || CLsplusb() == 0)
    return 0;
  
  double cl_b_err = (1. - CLb()) / (n_b * CLb());
  double cl_sb_err = (1. - CLsplusb()) / (n_sb * CLsplusb());
  
  return CLs() * TMath::Sqrt(cl_b_err + cl_sb_err);
}

///////////////////////////////////////////////////////////////////////////

void HybridResult::Add(HybridResult* other)
{
   // add additional toy-MC experiments to the current results
   // use the data test statistics of the added object if none is already present (otherwise, ignore the new one)

   int other_size_sb = other->GetTestStat_sb().size();
   for (int i=0;i<other_size_sb;++i)
      fTestStat_sb.push_back(other->GetTestStat_sb()[i]);

   int other_size_b = other->GetTestStat_b().size();
   for (int i=0;i<other_size_b;++i)
      fTestStat_b.push_back(other->GetTestStat_b()[i]);

   // if no data is present use the other's HybridResult's data
   if (fTestStat_data==-999.)
      fTestStat_data = other->GetTestStat_data();

   fComputationsAltDoneFlag = false;
   fComputationsNulDoneFlag = false;

   return;
}

///////////////////////////////////////////////////////////////////////////

HybridPlot* HybridResult::GetPlot(const char* name,const char* title, int n_bins)
{
   // prepare a plot showing a result and return a pointer to a HybridPlot object
   // the needed arguments are: an object name, a title and the number of bins in the plot

   // default plot name
   TString plot_name;
   if ( TString(name)=="" ) {
      plot_name += GetName();
      plot_name += "_plot";
   } else plot_name = name;

   // default plot title
   TString plot_title;
   if ( TString(title)=="" ) {
      plot_title += GetTitle();
      plot_title += "_plot (";
      plot_title += fTestStat_b.size();
      plot_title += " toys)";
   } else plot_title = title;

   HybridPlot* plot = new HybridPlot( plot_name.Data(),
                                      plot_title.Data(),
                                      fTestStat_sb,
                                      fTestStat_b,
                                      fTestStat_data,
                                      n_bins,
                                      true );
   return plot;
}

///////////////////////////////////////////////////////////////////////////

void HybridResult::PrintMore(const char* /*options */)
{
   // Print out some information about the results

   std::cout << "\nResults " << GetName() << ":\n"
             << " - Number of S+B toys: " << fTestStat_b.size() << std::endl
             << " - Number of B toys: " << fTestStat_sb.size() << std::endl
             << " - test statistics evaluated on data: " << fTestStat_data << std::endl
             << " - CL_b " << CLb() << std::endl
             << " - CL_s+b " << CLsplusb() << std::endl
             << " - CL_s " << CLs() << std::endl;

   return;
}










