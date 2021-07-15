// @(#)root/unuran:$Id$
// Authors: L. Moneta, J. Leydold Wed Feb 28 2007

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2010  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Implementation file for class TUnuranSampler
#include "TUnuranSampler.h"

#include "TUnuranContDist.h"
#include "TUnuranDiscrDist.h"
#include "TUnuranMultiContDist.h"
#include "TUnuran.h"
#include "Math/OneDimFunctionAdapter.h"
#include "Math/DistSamplerOptions.h"
#include "Math/GenAlgoOptions.h"
#include "Fit/DataRange.h"
//#include "Math/WrappedTF1.h"

#include "TRandom.h"
#include "TError.h"

#include "TF1.h"
#include <cassert>
#include <cmath>

ClassImp(TUnuranSampler);

TUnuranSampler::TUnuranSampler() : ROOT::Math::DistSampler(),
   fOneDim(false),
   fDiscrete(false),
   fHasMode(false), fHasArea(false),
   fMode(0), fArea(0),
   fUnuran(new TUnuran()  )
{
   fLevel = ROOT::Math::DistSamplerOptions::DefaultPrintLevel();
}

TUnuranSampler::~TUnuranSampler() {
   assert(fUnuran != 0);
   delete fUnuran;
}

bool TUnuranSampler::Init(const char * algo) {
   // initialize unuran classes using the given algorithm
   assert (fUnuran != 0 );
   bool ret = false;
   //case distribution has not been set
   // Maybe we are using the Unuran string API which contains also distribution string
   // try to initialize Unuran
   if (NDim() == 0)  {
      ret = fUnuran->Init(algo,"");
      if (!ret) { 
         Error("TUnuranSampler::Init",
         "Unuran initialization string is invalid or the Distribution function has not been set and one needs to call SetFunction first.");
         return false;
      }
      int ndim = fUnuran->GetDimension();
      assert(ndim > 0);
      fOneDim = (ndim == 1);
      fDiscrete = fUnuran->IsDistDiscrete(); 
      DoSetDimension(ndim);
      return true;
   }

   if (fLevel < 0) fLevel =  ROOT::Math::DistSamplerOptions::DefaultPrintLevel();

   TString method(algo);
   if (method.IsNull() ) {
      if (NDim() == 1) method = ROOT::Math::DistSamplerOptions::DefaultAlgorithm1D();
      else  method = ROOT::Math::DistSamplerOptions::DefaultAlgorithmND();
   }
   method.ToUpper();

   if (NDim() == 1) {
       // check if distribution is discrete by
      // using first string in the method name is "D"
      if (method.First("D") == 0) {
         if (fLevel>1) Info("TUnuranSampler::Init","Initialize one-dim discrete distribution with method %s",method.Data());
         ret =  DoInitDiscrete1D(method);
      }
      else {
         if (fLevel>1) Info("TUnuranSampler::Init","Initialize one-dim continuous distribution with method %s",method.Data());
         ret =  DoInit1D(method);
      }
   }
   else {
      if (fLevel>1) Info("TUnuranSampler::Init","Initialize multi-dim continuous distribution with method %s",method.Data());
      ret = DoInitND(method);
   }
   // set print level in UNURAN (must be done after having initialized) -
   if (fLevel>0) {
      //fUnuran->SetLogLevel(fLevel); ( seems not to work  disable for the time being)
      if (ret) Info("TUnuranSampler::Init","Successfully initailized Unuran with method %s",method.Data() );
      else Error("TUnuranSampler::Init","Failed to  initailize Unuran with method %s",method.Data() );
      // seems not to work in UNURAN (call only when level > 0 )
   }
   return ret;
}


bool TUnuranSampler::Init(const ROOT::Math::DistSamplerOptions & opt ) {
   // default initialization with algorithm name
   SetPrintLevel(opt.PrintLevel() );
   // check if there are extra options
   std::string optionStr = opt.Algorithm();
   auto extraOpts = opt.ExtraOptions();
   if (extraOpts) {
      ROOT::Math::GenAlgoOptions * opts = dynamic_cast<ROOT::Math::GenAlgoOptions*>(extraOpts);
      auto appendOption = [&](const std::string & key, const std::string & val) {
         optionStr += "; ";
         optionStr += key;
         if (!val.empty()) {
            optionStr += "=";
            optionStr += val;
         }
      };
      auto names = opts->GetAllNamedKeys();
      for ( auto & name : names) {
         std::string value = opts->NamedValue(name.c_str());
         appendOption(name,value);
      } 
      names = opts->GetAllIntKeys();
      for ( auto & name : names) {
         std::string value = ROOT::Math::Util::ToString(opts->IValue(name.c_str()));
         appendOption(name,value);
      } 
      names = opts->GetAllRealKeys();
      for ( auto & name : names) {
         std::string value = ROOT::Math::Util::ToString(opts->RValue(name.c_str()));
         appendOption(name,value);
      } 
   }
   Info("Init","Initialize UNU.RAN with Method option string: %s",optionStr.c_str());
   return Init(optionStr.c_str() );
}


bool TUnuranSampler::DoInit1D(const char * method) {
   // initialize for 1D sampling
   // need to create 1D interface from Multidim one
   // (to do: use directly 1D functions ??)
   // to do : add possibility for String API of UNURAN
   fOneDim = true;
   TUnuranContDist * dist = 0;
   if (fFunc1D == 0) {
      if (HasParentPdf()) {
         ROOT::Math::OneDimMultiFunctionAdapter<> function(ParentPdf() );
         dist = new TUnuranContDist(&function,fDPDF,fCDF,fUseLogPdf,true);
      }
      else {
         if (!fDPDF && !fCDF) {
            Error("DoInit1D", "No PDF, CDF or DPDF function has been set");
            return false;
         }
         dist = new TUnuranContDist(nullptr, fDPDF, fCDF, fUseLogPdf, true);
      }
   }
   else {
      dist = new TUnuranContDist(fFunc1D, fDPDF, fCDF, fUseLogPdf, true); // no need to copy the function
   }
   // set range in distribution (support only one range)
   const ROOT::Fit::DataRange & range = PdfRange();
   if (range.Size(0) > 0) {
      double xmin, xmax;
      range.GetRange(0,xmin,xmax);
      dist->SetDomain(xmin,xmax);
   }
   if (fHasMode) dist->SetMode(fMode);
   if (fHasArea) dist->SetPdfArea(fArea);

   bool ret = false;
   if (method) ret =  fUnuran->Init(*dist, method);
   else ret =  fUnuran->Init(*dist);
   delete dist;
   return ret;
}

bool TUnuranSampler::DoInitDiscrete1D(const char * method) {
   // initialize for 1D sampling of discrete distributions
   fOneDim = true;
   fDiscrete = true;
   TUnuranDiscrDist * dist = 0;
   if (fFunc1D == 0) {
      if (!HasParentPdf()) {
         Error("DoInitDiscrete1D", "No PMF has been defined");
         return false;
      }
      // need to copy the passed function pointer in this case
      ROOT::Math::OneDimMultiFunctionAdapter<> function(ParentPdf() );
      dist = new TUnuranDiscrDist(function,true);
   }
   else {
      // no need to copy the function since fFunc1D is managed outside
      dist = new TUnuranDiscrDist(*fFunc1D, false);
   }
   // set CDF if available
   if (fCDF) dist->SetCdf(*fCDF);
   // set range in distribution (support only one range)
   // otherwise 0, inf is assumed
   const ROOT::Fit::DataRange & range = PdfRange();
   if (range.Size(0) > 0) {
      double xmin, xmax;
      range.GetRange(0,xmin,xmax);
      if (xmin < 0) {
         Warning("DoInitDiscrete1D","range starts from negative values - set minimum to zero");
         xmin = 0;
      }
      dist->SetDomain(int(xmin+0.1),int(xmax+0.1));
   }
   if (fHasMode) dist->SetMode(int(fMode+0.1));
   if (fHasArea) dist->SetProbSum(fArea);

   bool ret =  fUnuran->Init(*dist, method);
   delete dist;
   return ret;
}


bool TUnuranSampler::DoInitND(const char * method) {
   // initialize for ND sampling
   if (!HasParentPdf()) {
      Error("DoInitND", "No PDF has been defined");
      return false;
   }
   TUnuranMultiContDist dist(ParentPdf());
   // set range in distribution (support only one range)
   const ROOT::Fit::DataRange & range = PdfRange();
   if (range.IsSet()) {
      std::vector<double> xmin(range.NDim() );
      std::vector<double> xmax(range.NDim() );
      range.GetRange(&xmin[0],&xmax[0]);
      dist.SetDomain(&xmin.front(),&xmax.front());
//       std::cout << " range is min = ";
//       for (int j = 0; j < NDim(); ++j) std::cout << xmin[j] << "   ";
//       std::cout << " max = ";
//       for (int j = 0; j < NDim(); ++j) std::cout << xmax[j] << "   ";
//       std::cout << std::endl;
   }
   fOneDim = false;
   if (fHasMode && fNDMode.size() == dist.NDim())
      dist.SetMode(fNDMode.data());

   if (method) return fUnuran->Init(dist, method);
   return fUnuran->Init(dist);
}

void TUnuranSampler::SetFunction(TF1 * pdf) {
   // set function from a TF1 pointer
   SetFunction<TF1>(*pdf, pdf->GetNdim());
}

void TUnuranSampler::SetRandom(TRandom * r) {
   // set random generator (must be called before Init to have effect)
   fUnuran->SetRandom(r);
}

void TUnuranSampler::SetSeed(unsigned int seed) {
   // set random generator seed (must be called before Init to have effect)
   fUnuran->SetSeed(seed);
}

TRandom * TUnuranSampler::GetRandom() {
   // get random generator used
   return  fUnuran->GetRandom();
}

double TUnuranSampler::Sample1D() {
   // sample 1D distributions
   return (fDiscrete) ? (double) fUnuran->SampleDiscr() : fUnuran->Sample();
}

bool TUnuranSampler::Sample(double * x) {
   // sample multi-dim distributions
   if (!fOneDim) return fUnuran->SampleMulti(x);
   x[0] = Sample1D();
   return true;
}


bool TUnuranSampler::SampleBin(double prob, double & value, double *error) {
   // sample a bin according to Poisson statistics
   TRandom * r = fUnuran->GetRandom();
   if (!r) return false;
   value = r->Poisson(prob);
   if (error) *error = std::sqrt(prob);
   return true;
}

void TUnuranSampler::SetMode(const std::vector<double> &mode)
{
   // set modes for multidim distribution
   if (mode.size() == ParentPdf().NDim()) {
      if (mode.size() == 1)
         fMode = mode[0];
      else 
         fNDMode = mode;

      fHasMode = true;
   }
   else {
      Error("SetMode", "modes vector is not compatible with function dimension of %d", (int)ParentPdf().NDim());
      fHasMode = false;
      fNDMode.clear();
   }
}

void TUnuranSampler::SetCdf(const ROOT::Math::IGenFunction &cdf) {
   fCDF = &cdf;
   // in case dimension has not been defined ( a pdf is not provided)
   if (NDim() == 0) DoSetDimension(1);
}

void TUnuranSampler::SetDPdf(const ROOT::Math::IGenFunction &dpdf) { 
   fDPDF = &dpdf;
   // in case dimension has not been defined ( a pdf is not provided)
   if (NDim() == 0) DoSetDimension(1);
}