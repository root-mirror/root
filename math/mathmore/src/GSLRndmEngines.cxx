// @(#)root/mathmore:$Id$
// Authors: L. Moneta, A. Zsenei   08/2005

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2004 ROOT Foundation,  CERN/PH-SFT                   *
  *                                                                    *
  * This library is free software; you can redistribute it and/or      *
  * modify it under the terms of the GNU General Public License        *
  * as published by the Free Software Foundation; either version 2     *
  * of the License, or (at your option) any later version.             *
  *                                                                    *
  * This library is distributed in the hope that it will be useful,    *
  * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU   *
  * General Public License for more details.                           *
  *                                                                    *
  * You should have received a copy of the GNU General Public License  *
  * along with this library (see file COPYING); if not, write          *
  * to the Free Software Foundation, Inc., 59 Temple Place, Suite      *
  * 330, Boston, MA 02111-1307 USA, or contact the author.             *
  *                                                                    *
  **********************************************************************/

// Header file for class GSLRandom
//
// Created by: moneta  at Sun Nov 21 16:26:03 2004
//
// Last update: Sun Nov 21 16:26:03 2004
//



// need to be included later
#include <ctime>
#include <cassert>

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"


#include "Math/GSLRndmEngines.h"
#include "GSLRngWrapper.h"
// for wrapping in GSL ROOT engines
#include "GSLRngROOTWrapper.h"

extern double gsl_ran_gaussian_acr(  const gsl_rng * r, const double sigma);

namespace ROOT {
namespace Math {





  // default constructor (need to call set type later)
   GSLRandomEngine::GSLRandomEngine() :
      fRng(nullptr),
      fCurTime(0)
  { }

   // constructor from external rng
   // internal generator will be managed or not depending on
   // how the GSLRngWrapper is created
   GSLRandomEngine::GSLRandomEngine( GSLRngWrapper * rng) :
      fRng(new GSLRngWrapper(*rng) ),
      fCurTime(0)
   {}

   // copy constructor
   GSLRandomEngine::GSLRandomEngine(const GSLRandomEngine & eng) :
      fRng(new GSLRngWrapper(*eng.fRng) ),
      fCurTime(0)
   {}

   GSLRandomEngine::~GSLRandomEngine() {
      // destructor : call terminate if not yet called
      if (fRng) Terminate();
   }

   // assignment operator
   GSLRandomEngine & GSLRandomEngine::operator=(const GSLRandomEngine & eng) {
      if (this == &eng) return *this;
      if (fRng)
         *fRng = *eng.fRng;
      else
         fRng = new GSLRngWrapper(*eng.fRng);
      fCurTime = eng.fCurTime;
      return *this;
   }


   void GSLRandomEngine::Initialize() {
      // initialize the generator by allocating the GSL object
      // if type was not passed create with default generator
      if (!fRng) fRng = new GSLRngWrapper();
      fRng->Allocate();
   }

   void GSLRandomEngine::Terminate() {
      // terminate the generator by freeing the GSL object
      if (!fRng) return;
      fRng->Free();
      delete fRng;
      fRng = 0;
   }


   double GSLRandomEngine::operator() () const {
      // generate random between 0 and 1.
      // 0 is excluded
      return gsl_rng_uniform_pos(fRng->Rng() );
   }


   unsigned long GSLRandomEngine::RndmInt(unsigned long max) const {
      // generate a random integer number between 0  and MAX
      return gsl_rng_uniform_int( fRng->Rng(), max );
   }

   unsigned long GSLRandomEngine::MinInt() const {
      // return minimum integer value used in RndmInt
      return gsl_rng_min( fRng->Rng() );
   }

   unsigned long GSLRandomEngine::MaxInt() const {
      // return maximum integr value used in RndmInt
      return gsl_rng_max( fRng->Rng() );
   }

   void GSLRandomEngine::RandomArray(double * begin, double * end )  const {
      // generate array of randoms betweeen 0 and 1. 0 is excluded
      // specialization for double * (to be faster)
      for ( double * itr = begin; itr != end; ++itr ) {
         *itr = gsl_rng_uniform_pos(fRng->Rng() );
      }
   }

   void GSLRandomEngine::SetSeed(unsigned int seed) const  {
      // set the seed, if = 0then the seed is set randomly using an std::rand()
      // seeded with the current time. Be carefuk in case the current time is
      // the same in consecutive calls
      if (seed == 0) {
         // use like in root (use time)
         time_t curtime;
         time(&curtime);
         unsigned int ct = static_cast<unsigned int>(curtime);
         if (ct != fCurTime) {
            fCurTime = ct;
            // set the seed for rand
            srand(ct);
         }
         seed = rand();
      }

      assert(fRng);
      gsl_rng_set(fRng->Rng(), seed );
   }

   std::string GSLRandomEngine::Name() const {
      //////////////////////////////////////////////////////////////////////////

      assert ( fRng != 0);
      assert ( fRng->Rng() != 0 );
      return std::string( gsl_rng_name( fRng->Rng() ) );
   }

   unsigned int GSLRandomEngine::Size() const {
      //////////////////////////////////////////////////////////////////////////

      assert (fRng != 0);
      return gsl_rng_size( fRng->Rng() );
   }


   // Random distributions

   double GSLRandomEngine::GaussianZig(double sigma)  const
   {
      // Gaussian distribution. Use fast ziggurat algorithm implemented since GSL 1.8
      return gsl_ran_gaussian_ziggurat(  fRng->Rng(), sigma);
   }

   double GSLRandomEngine::Gaussian(double sigma)  const
   {
      // Gaussian distribution. Use default Box-Muller method
      return gsl_ran_gaussian(  fRng->Rng(), sigma);
   }

   double GSLRandomEngine::GaussianRatio(double sigma)  const
   {
      // Gaussian distribution. Use ratio method
      return gsl_ran_gaussian_ratio_method(  fRng->Rng(), sigma);
   }


   double GSLRandomEngine::GaussianTail(double a , double sigma) const
   {
      // Gaussian Tail distribution: eeturn values larger than a distributed
      // according to the gaussian
      return gsl_ran_gaussian_tail(  fRng->Rng(), a, sigma);
   }

   void GSLRandomEngine::Gaussian2D(double sigmaX, double sigmaY, double rho, double &x, double &y) const
   {
      // Gaussian Bivariate distribution, with correlation coefficient rho
      gsl_ran_bivariate_gaussian(  fRng->Rng(), sigmaX, sigmaY, rho, &x, &y);
   }

   double GSLRandomEngine::Exponential(double mu)  const
   {
      // Exponential distribution
      return gsl_ran_exponential(  fRng->Rng(), mu);
   }

   double GSLRandomEngine::Cauchy(double a) const
   {
      // Cauchy distribution
      return gsl_ran_cauchy(  fRng->Rng(), a);
   }

   double GSLRandomEngine::Landau() const
   {
      // Landau distribution
      return gsl_ran_landau(  fRng->Rng());
   }

   double GSLRandomEngine::Beta(double a, double b) const
   {
      // Beta distribution
      return gsl_ran_beta(  fRng->Rng(), a, b);
   }

   double GSLRandomEngine::Gamma(double a, double b) const
   {
      // Gamma distribution
      return gsl_ran_gamma(  fRng->Rng(), a, b);
   }

   double GSLRandomEngine::LogNormal(double zeta, double sigma) const
   {
      // Log normal distribution
      return gsl_ran_lognormal(  fRng->Rng(), zeta, sigma);
   }

   double GSLRandomEngine::ChiSquare(double nu) const
   {
      // Chi square distribution
      return gsl_ran_chisq(  fRng->Rng(), nu);
   }


   double GSLRandomEngine::FDist(double nu1, double nu2)  const
   {
      // F distribution
      return gsl_ran_fdist(  fRng->Rng(), nu1, nu2);
   }

   double GSLRandomEngine::tDist(double nu)  const
   {
      // t distribution
      return gsl_ran_tdist(  fRng->Rng(), nu);
   }

   double GSLRandomEngine::Rayleigh(double sigma)  const
   {
      // Rayleigh distribution
      return gsl_ran_rayleigh(  fRng->Rng(), sigma);
   }

   double GSLRandomEngine::Logistic(double a)  const
   {
      // Logistic distribution
      return gsl_ran_logistic(  fRng->Rng(), a);
   }

   double GSLRandomEngine::Pareto(double a, double b)  const
   {
      // Pareto distribution
      return gsl_ran_pareto(  fRng->Rng(), a, b);
   }

   void GSLRandomEngine::Dir2D(double &x, double &y) const
   {
      // generate random numbers in a 2D circle of radious 1
      gsl_ran_dir_2d(  fRng->Rng(), &x, &y);
   }

   void GSLRandomEngine::Dir3D(double &x, double &y, double &z) const
   {
      // generate random numbers in a 3D sphere of radious 1
      gsl_ran_dir_3d(  fRng->Rng(), &x, &y, &z);
   }

   unsigned int GSLRandomEngine::Poisson(double mu) const
   {
      // Poisson distribution
      return gsl_ran_poisson(  fRng->Rng(), mu);
   }

   unsigned int GSLRandomEngine::Binomial(double p, unsigned int n) const
   {
      // Binomial distribution
      return gsl_ran_binomial(  fRng->Rng(), p, n);
   }

   unsigned int GSLRandomEngine::NegativeBinomial(double p, double n) const
   {
      // Negative Binomial distribution
      return gsl_ran_negative_binomial(  fRng->Rng(), p, n);
   }


   std::vector<unsigned int>  GSLRandomEngine::Multinomial( unsigned int ntot, const std::vector<double> & p ) const
   {
      // Multinomial distribution  return vector of integers which sum is ntot
      std::vector<unsigned int> ival( p.size());
      gsl_ran_multinomial(  fRng->Rng(), p.size(), ntot, &p.front(), &ival[0]);
      return ival;
   }



   //----------------------------------------------------
   // generators
   //----------------------------------------------------

   /////////////////////////////////////////////////////////////////////////////

   GSLRngMT::GSLRngMT() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_mt19937));
      Initialize(); 
   }


   // old ranlux - equivalent to TRandom1
   GSLRngRanLux::GSLRngRanLux() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_ranlux) );
      Initialize(); 
   }

   // second generation of Ranlux (single precision version - luxury 1)
   GSLRngRanLuxS1::GSLRngRanLuxS1() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_ranlxs1) );
      Initialize(); 
   }

   // second generation of Ranlux (single precision version - luxury 2)
   GSLRngRanLuxS2::GSLRngRanLuxS2() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_ranlxs2) );
      Initialize(); 
   }

   // double precision  version - luxury 1
   GSLRngRanLuxD1::GSLRngRanLuxD1() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_ranlxd1) );
      Initialize(); 
   }

   // double precision  version - luxury 2
   GSLRngRanLuxD2::GSLRngRanLuxD2() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_ranlxd2) );
      Initialize(); 
   }

   /////////////////////////////////////////////////////////////////////////////

   GSLRngTaus::GSLRngTaus() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_taus2) );
      Initialize(); 
   }

   /////////////////////////////////////////////////////////////////////////////

   GSLRngGFSR4::GSLRngGFSR4() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_gfsr4) );
      Initialize(); 
   }

   /////////////////////////////////////////////////////////////////////////////

   GSLRngCMRG::GSLRngCMRG() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_cmrg) );
      Initialize(); 
   }

   /////////////////////////////////////////////////////////////////////////////

   GSLRngMRG::GSLRngMRG() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_mrg) );
      Initialize(); 
   }


   /////////////////////////////////////////////////////////////////////////////

   GSLRngRand::GSLRngRand() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_rand) );
      Initialize(); 
   }

   /////////////////////////////////////////////////////////////////////////////

   GSLRngRanMar::GSLRngRanMar() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_ranmar) );
      Initialize(); 
   }

   /////////////////////////////////////////////////////////////////////////////

   GSLRngMinStd::GSLRngMinStd() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_minstd) );
      Initialize(); 
   }


   // for extra engines based on ROOT
   GSLRngMixMax::GSLRngMixMax() : GSLRandomEngine()
   {
      SetType(new GSLRngWrapper(gsl_rng_mixmax) );
      Initialize(); // this creates the gsl_rng structure
      //  no real need to call CreateEngine since the underlined MIXMAX engine is created
      // by calling GSLMixMaxWrapper::Seed(gsl_default_seed) that is called 
      // when gsl_rng is allocated (in Initialize) 
      GSLMixMaxWrapper::CreateEngine(Engine()->Rng());
   }
   GSLRngMixMax::~GSLRngMixMax() {
      // we need to explicitly delete the ROOT wrapper class
      GSLMixMaxWrapper::FreeEngine(Engine()->Rng());
   }

} // namespace Math
} // namespace ROOT



