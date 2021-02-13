// @(#)root/unuran:$Id$
// Authors: L. Moneta, J. Leydold Tue Sep 26 16:25:09 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Implementation file for class TUnuran

#include "TUnuran.h"

#include "TUnuranContDist.h"
#include "TUnuranMultiContDist.h"
#include "TUnuranDiscrDist.h"
#include "TUnuranEmpDist.h"

#include "UnuranRng.h"
#include "UnuranDistrAdapter.h"

#include "TRandom.h"

#include <cassert>

#include <unuran.h>

#include "TError.h"


TUnuran::TUnuran(TRandom * r, unsigned int debugLevel) :
   fGen(0),
   fUdistr(0),
   fUrng(0),
   fRng(r)
{
   // constructor implementation with a ROOT random generator
   // if no generator is given the ROOT default is used
   if (fRng == 0) fRng = gRandom;
   // set debug level at global level
   // (should be in a static  initialization function of the library ? )
   if ( debugLevel > 1)
      unur_set_default_debug(UNUR_DEBUG_ALL);
   else if (debugLevel == 1)
      unur_set_default_debug(UNUR_DEBUG_INIT);
   else
      unur_set_default_debug(UNUR_DEBUG_OFF);

}


TUnuran::~TUnuran()
{
   // Destructor implementation
   if (fGen != 0) unur_free(fGen);
   if (fUrng != 0) unur_urng_free(fUrng);
  // we can delete now the distribution object
   if (fUdistr != 0) unur_distr_free(fUdistr);
}

//private (no impl.)
TUnuran::TUnuran(const TUnuran &)
{
   // Implementation of copy constructor.
}

TUnuran & TUnuran::operator = (const TUnuran &rhs)
{
   // Implementation of assignment operator.
   if (this == &rhs) return *this;  // time saving self-test
   return *this;
}

bool  TUnuran::Init(const std::string & dist, const std::string & method)
{
   // initialize with a string
   std::string s = dist + " & " + method;
   fGen = unur_str2gen(s.c_str() );
   if (fGen == 0) {
      Error("Init","Cannot create generator object");
      return false;
   }
   if (! SetRandomGenerator() ) return false;

   return true;
}

bool TUnuran::Init(const TUnuranContDist & distr, const std::string  & method)
{
   // initialization with a distribution and and generator
   // the distribution object is copied in and managed by this class
   // use std::unique_ptr to manage previously existing distribution objects
   TUnuranContDist * distNew = distr.Clone();
   fDist.reset(distNew);

   fMethod = method;
   if (! SetContDistribution(*distNew) ) return false;
   if (! SetMethodAndInit() ) return false;
   if (! SetRandomGenerator() ) return false;
   return true;
}


bool TUnuran::Init(const TUnuranMultiContDist & distr, const std::string  & method)
{
   //  initialization with a distribution and method
   // the distribution object is copied in and managed by this class
   // use std::unique_ptr to manage previously existing distribution objects
   TUnuranMultiContDist * distNew = distr.Clone();
   fDist.reset(distNew);

   fMethod = method;
   if (! SetMultiDistribution(*distNew) ) return false;
   if (! SetMethodAndInit() ) return false;
   if (! SetRandomGenerator() ) return false;
   return true;
}


bool TUnuran::Init(const TUnuranDiscrDist & distr, const std::string & method ) {
   //   initialization with a distribution and and generator
   // the distribution object is copied in and managed by this class
   // use std::unique_ptr to manage previously existing distribution objects
   TUnuranDiscrDist * distNew = distr.Clone();
   fDist.reset(distNew);

   fMethod = method;
   if (! SetDiscreteDistribution(*distNew) ) return false;
   if (! SetMethodAndInit() ) return false;
   if (! SetRandomGenerator() ) return false;
   return true;
}

bool TUnuran::Init(const TUnuranEmpDist & distr, const std::string & method ) {
   //   initialization with a distribution and and generator
   // the distribution object is copied in and managed by this class
   // use std::unique_ptr to manage previously existing distribution objects
   TUnuranEmpDist * distNew = distr.Clone();
   fDist.reset(distNew);

   fMethod = method;
   if (distr.IsBinned()) fMethod = "hist";
   else if (distr.NDim() > 1) fMethod = "vempk";
   if (! SetEmpiricalDistribution(*distNew) ) return false;
   if (! SetMethodAndInit() ) return false;
   if (! SetRandomGenerator() ) return false;
   return true;
}


bool  TUnuran::SetRandomGenerator()
{
   // set an external random generator
   if (fRng == 0) return false;
   if (fGen == 0) return false;

   fUrng = unur_urng_new(&UnuranRng<TRandom>::Rndm, fRng );
   if (fUrng == 0) return false;
   unsigned int ret = 0;
   ret |= unur_urng_set_delete(fUrng, &UnuranRng<TRandom>::Delete);
   ret |= unur_urng_set_seed(fUrng, &UnuranRng<TRandom>::Seed);
   if (ret != 0) return false;

   unur_chg_urng( fGen, fUrng);
   return true;
}

bool  TUnuran::SetContDistribution(const TUnuranContDist & dist )
{
   // internal method to set in unuran the function pointer for a continuous univariate distribution
   if (fUdistr != 0)  unur_distr_free(fUdistr);
   fUdistr = unur_distr_cont_new();
   if (fUdistr == 0) return false;
   unsigned int ret = 0;
   ret = unur_distr_set_extobj(fUdistr, &dist);
   if ( ! dist.IsLogPdf() ) {
      ret |= unur_distr_cont_set_pdf(fUdistr, &ContDist::Pdf);
      ret |= unur_distr_cont_set_dpdf(fUdistr, &ContDist::Dpdf);
      if (dist.HasCdf() ) ret |= unur_distr_cont_set_cdf(fUdistr, &ContDist::Cdf);
   }
   else {
      // case user provides log of pdf
      ret |= unur_distr_cont_set_logpdf(fUdistr, &ContDist::Pdf);
      ret |= unur_distr_cont_set_dlogpdf(fUdistr, &ContDist::Dpdf);
   }

   double xmin, xmax = 0;
   if (dist.GetDomain(xmin,xmax) ) {
      ret = unur_distr_cont_set_domain(fUdistr,xmin,xmax);
      if (ret != 0)  {
         Error("SetContDistribution","invalid domain xmin = %g xmax = %g ",xmin,xmax);
         return false;
      }
   }
   if (dist.HasMode() ) {
      ret = unur_distr_cont_set_mode(fUdistr, dist.Mode());
      if (ret != 0)  {
         Error("SetContDistribution","invalid mode given,  mode = %g ",dist.Mode());
         return false;
      }
   }
   if (dist.HasPdfArea() ) {
      ret = unur_distr_cont_set_pdfarea(fUdistr, dist.PdfArea());
      if (ret != 0)  {
         Error("SetContDistribution","invalid area given,  area = %g ",dist.PdfArea());
         return false;
      }
   }

   return (ret ==0) ? true : false;
}


bool  TUnuran::SetMultiDistribution(const TUnuranMultiContDist & dist )
{
   // internal method to set in unuran the function pointer for a multivariate distribution
   if (fUdistr != 0)  unur_distr_free(fUdistr);
   fUdistr = unur_distr_cvec_new(dist.NDim() );
   if (fUdistr == 0) return false;
   unsigned int ret = 0;
   ret |= unur_distr_set_extobj(fUdistr, &dist );
   if ( ! dist.IsLogPdf() ) {
      ret |= unur_distr_cvec_set_pdf(fUdistr, &MultiDist::Pdf);
      ret |= unur_distr_cvec_set_dpdf(fUdistr, &MultiDist::Dpdf);
      ret |= unur_distr_cvec_set_pdpdf(fUdistr, &MultiDist::Pdpdf);
   }
   else {
      ret |= unur_distr_cvec_set_logpdf(fUdistr, &MultiDist::Pdf);
      ret |= unur_distr_cvec_set_dlogpdf(fUdistr, &MultiDist::Dpdf);
      ret |= unur_distr_cvec_set_pdlogpdf(fUdistr, &MultiDist::Pdpdf);
   }

   const double * xmin = dist.GetLowerDomain();
   const double * xmax = dist.GetUpperDomain();
   if ( xmin != 0 || xmax != 0 ) {
      ret = unur_distr_cvec_set_domain_rect(fUdistr,xmin,xmax);
      if (ret != 0)  {
         Error("SetMultiDistribution","invalid domain");
         return false;
      }
#ifdef OLDVERS
      Error("SetMultiDistribution","domain setting not available in UNURAN 0.8.1");
#endif

   }

   const double * xmode = dist.GetMode();
   if (xmode != 0) {
      ret = unur_distr_cvec_set_mode(fUdistr, xmode);
      if (ret != 0)  {
         Error("SetMultiDistribution","invalid mode");
         return false;
      }
   }
   return (ret ==0) ? true : false;
}

bool TUnuran::SetEmpiricalDistribution(const TUnuranEmpDist & dist) {

   // internal method to set in unuran the function pointer for am empiral distribution (from histogram)
   if (fUdistr != 0)  unur_distr_free(fUdistr);
   if (dist.NDim() == 1)
      fUdistr = unur_distr_cemp_new();
   else
      fUdistr = unur_distr_cvemp_new(dist.NDim() );

   if (fUdistr == 0) return false;
   unsigned int ret = 0;


   // get info from histogram
   if (dist.IsBinned() ) {
      int nbins = dist.Data().size();
      double min = dist.LowerBin();
      double max = dist.UpperBin();
      const double * pv = &(dist.Data().front());
      ret |= unur_distr_cemp_set_hist(fUdistr, pv, nbins, min, max);
#ifdef OLDVERS
      Error("SetEmpiricalDistribution","hist method not available in UNURAN 0.8.1");
#endif
   }
   else {
      const double * pv = &dist.Data().front();
      // n is number of points (size/ndim)
      int n = dist.Data().size()/dist.NDim();
      if (dist.NDim() == 1)
         ret |= unur_distr_cemp_set_data(fUdistr, pv, n);
      else
         ret |= unur_distr_cvemp_set_data(fUdistr, pv, n);
   }
   if (ret != 0) {
      Error("SetEmpiricalDistribution","invalid distribution object");
      return false;
   }
   return true;
}


bool  TUnuran::SetDiscreteDistribution(const TUnuranDiscrDist & dist)
{
   // internal method to set in unuran the function pointer for a discrete univariate distribution
   if (fUdistr != 0)  unur_distr_free(fUdistr);
   fUdistr = unur_distr_discr_new();
   if (fUdistr == 0) return false;
   unsigned int ret = 0;
   // if a probability mesh function is provided
   if (dist.ProbVec().size() == 0) {
      ret = unur_distr_set_extobj(fUdistr, &dist );
      ret |= unur_distr_discr_set_pmf(fUdistr, &DiscrDist::Pmf);
      if (dist.HasCdf() ) ret |= unur_distr_discr_set_cdf(fUdistr, &DiscrDist::Cdf);

   }
   else {
      // case user provides vector of probabilities
      ret |= unur_distr_discr_set_pv(fUdistr, &dist.ProbVec().front(), dist.ProbVec().size() );
   }

   int xmin, xmax = 0;
   if (dist.GetDomain(xmin,xmax) ) {
      ret = unur_distr_discr_set_domain(fUdistr,xmin,xmax);
      if (ret != 0)  {
         Error("SetDiscrDistribution","invalid domain xmin = %d xmax = %d ",xmin,xmax);
         return false;
      }
   }
   if (dist.HasMode() ) {
      ret = unur_distr_discr_set_mode(fUdistr, dist.Mode());
      if (ret != 0)  {
         Error("SetContDistribution","invalid mode given,  mode = %d ",dist.Mode());
         return false;
      }
   }
   if (dist.HasProbSum() ) {
      ret = unur_distr_discr_set_pmfsum(fUdistr, dist.ProbSum());
      if (ret != 0)  {
         Error("SetContDistribution","invalid sum given,  mode = %g ",dist.ProbSum());
         return false;
      }
   }

   return (ret ==0) ? true : false;
}


//bool TUnuran::SetMethodAndInit(const std::string & s) {
bool TUnuran::SetMethodAndInit() {

   // internal function to set a method from a distribution and
   // initialize unuran with the given distribution.
   if (fUdistr == 0) return false;

   struct unur_slist *mlist = NULL;

   UNUR_PAR * par = _unur_str2par(fUdistr, fMethod.c_str(), &mlist);
   if (par == 0) {
      Error("SetMethod","missing distribution information or syntax error");
      if (mlist != 0)  _unur_slist_free(mlist);
      return false;
   }


   // set unuran to not use a private copy of the distribution object
   unur_set_use_distr_privatecopy (par, false);

   // need to free fGen if already existing ?
   if (fGen != 0 )  unur_free(fGen);
   fGen = unur_init(par);
   _unur_slist_free(mlist);
   if (fGen == 0) {
      Error("SetMethod","initializing Unuran: condition for method violated");
      return false;
   }
   return true;
 }


int TUnuran::SampleDiscr()
{
   // sample one-dimensional distribution
   assert(fGen != 0);
   return unur_sample_discr(fGen);
}

double TUnuran::Sample()
{
   // sample one-dimensional distribution
   assert(fGen != 0);
   return unur_sample_cont(fGen);
}

bool TUnuran::SampleMulti(double * x)
{
   // sample multidimensional distribution
   if (fGen == 0) return false;
   unur_sample_vec(fGen,x);
   return true;
}

void TUnuran::SetSeed(unsigned int seed) {
   return fRng->SetSeed(seed);
}

bool  TUnuran::SetLogLevel(unsigned int debugLevel)
{
   if (fGen == 0) return false;
   int ret = 0;
   if ( debugLevel > 1)
      ret |= unur_chg_debug(fGen, UNUR_DEBUG_ALL);
   else if (debugLevel == 1)
      ret |= unur_chg_debug(fGen, UNUR_DEBUG_ALL);
   else
      ret |= unur_chg_debug(fGen, UNUR_DEBUG_OFF);

   return (ret ==0) ? true : false;

}

bool TUnuran::InitPoisson(double mu, const std::string & method) {
   // initializaton for a Poisson
   double p[1];
   p[0] = mu;

   fUdistr = unur_distr_poisson(p,1);

   fMethod = method;
   if (fUdistr == 0) return false;
   if (! SetMethodAndInit() ) return false;
   if (! SetRandomGenerator() ) return false;
   return true;
}


bool TUnuran::InitBinomial(unsigned int ntot, double prob, const std::string & method ) {
   // initializaton for a Binomial
   double par[2];
   par[0] = ntot;
   par[1] = prob;
   fUdistr = unur_distr_binomial(par,2);

   fMethod = method;
   if (fUdistr == 0) return false;
   if (! SetMethodAndInit() ) return false;
   if (! SetRandomGenerator() ) return false;
   return true;
}


bool TUnuran::ReInitDiscrDist(unsigned int npar, double * par) {
   // re-initialization of UNURAN without freeing and creating a new fGen object
   // works only for pre-defined distribution by changing their parameters
   if (!fGen ) return false;
   if (!fUdistr) return false;
   unur_distr_discr_set_pmfparams(fUdistr,par,npar);
   int iret = unur_reinit(fGen);
   if (iret) Warning("ReInitDiscrDist","re-init failed - a full initizialization must be performed");
   return (!iret);
}

