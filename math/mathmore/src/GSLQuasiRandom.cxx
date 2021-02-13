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
//
//
// Created by: moneta  at Sun Nov 21 :26:03 2004
//
// Last update: Sun Nov 21 16:26:03 2004
//

// need to be included later
#include <ctime>
#include <vector>
#include <cstdlib>
#include <cassert>

#include "Math/GSLQuasiRandom.h"
#include "GSLQRngWrapper.h"


namespace ROOT {
namespace Math {





  // default constructor (need to call set type later)
   GSLQuasiRandomEngine::GSLQuasiRandomEngine() :
      fQRng(0 )
  { }

   // constructor from external rng
   // internal generator will be managed or not depending on
   // how the GSLQRngWrapper is created
   GSLQuasiRandomEngine::GSLQuasiRandomEngine( GSLQRngWrapper * rng) :
      fQRng(new GSLQRngWrapper(*rng) )
   {}

   // copy constructor
   GSLQuasiRandomEngine::GSLQuasiRandomEngine(const GSLQuasiRandomEngine & eng) :
      fQRng(new GSLQRngWrapper(*eng.fQRng) )
   {}

   GSLQuasiRandomEngine::~GSLQuasiRandomEngine() {
      // destructor : call terminate if not yet called
      if (fQRng) Terminate();
   }

   // assignment operator
   GSLQuasiRandomEngine & GSLQuasiRandomEngine::operator=(const GSLQuasiRandomEngine & eng) {
      if (this == &eng) return *this;
      if (fQRng)
         *fQRng = *eng.fQRng;
      else
         fQRng = new GSLQRngWrapper(*eng.fQRng);
      return *this;
   }


   void GSLQuasiRandomEngine::Initialize(unsigned int dimension) {
      // initialize the generator by allocating the GSL object
      // if type was not passed create with default generator
      if (!fQRng) fQRng = new GSLQRngWrapper();
      fQRng->Allocate(dimension);
   }

   void GSLQuasiRandomEngine::Terminate() {
      // terminate the generator by freeing the GSL object
      if (!fQRng) return;
      fQRng->Free();
      delete fQRng;
      fQRng = 0;
   }


   double GSLQuasiRandomEngine::operator() () const {
      // generate next point in the quasi random sequence. The generate number x is 0 < x < 1
      // with 0 and 1 excluded
      // This method should be called only if dimension == 1
      assert(fQRng->Dimension() == 1);
      double x;
      gsl_qrng_get(fQRng->Rng(), &x );
      return x;
   }

   bool GSLQuasiRandomEngine::operator() (double * x) const {
      // generate next point in the quasi random sequence. The generate number x is 0 < x < 1
      // with 0 and 1 excluded
      int status = gsl_qrng_get(fQRng->Rng(), x );
      return (status == 0);
   }

   bool GSLQuasiRandomEngine::Skip(unsigned int n) const {
      // throw away the next n random numbers
      std::vector<double> xtmp(fQRng->Dimension() );
      int status = 0;
      for (unsigned int i = 0; i < n; ++i ) {
         status |= gsl_qrng_get(fQRng->Rng(), &xtmp[0] );
      }
      return status == 0;
   }

   bool GSLQuasiRandomEngine::GenerateArray(double * begin, double * end )  const {
      // generate array of randoms betweeen 0 and 1. 0 is excluded
      // specialization for double * (to be faster)
      int status = 0;
      for ( double * itr = begin; itr != end; itr+=fQRng->Dimension() ) {
         status |= gsl_qrng_get(fQRng->Rng(), itr );
      }
      return status == 0;
   }


   std::string GSLQuasiRandomEngine::Name() const {
      //////////////////////////////////////////////////////////////////////////

      assert (fQRng != 0);
      assert(fQRng->Rng() != 0);
      const char * name = gsl_qrng_name( fQRng->Rng() );
      if (!name)  return std::string();
      return std::string( name);
   }

   unsigned int GSLQuasiRandomEngine::Size() const {
      //////////////////////////////////////////////////////////////////////////

      assert (fQRng != 0);
      return gsl_qrng_size( fQRng->Rng() );
   }

   unsigned int GSLQuasiRandomEngine::NDim() const {
      //////////////////////////////////////////////////////////////////////////

      assert (fQRng != 0);
      return fQRng->Dimension();
   }




   //----------------------------------------------------
   // generators
   //----------------------------------------------------

   //----------------------------------------------------
   // generator based on Sobol sequence
   GSLQRngSobol::GSLQRngSobol() : GSLQuasiRandomEngine()
   {
      SetType(new GSLQRngWrapper(gsl_qrng_sobol));
   }


   // generator based on Bratley-Fox,Niederreiter sequence
   GSLQRngNiederreiter2::GSLQRngNiederreiter2() : GSLQuasiRandomEngine()
   {
      SetType(new GSLQRngWrapper(gsl_qrng_niederreiter_2) );
   }





} // namespace Math
} // namespace ROOT



