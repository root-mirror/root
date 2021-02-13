// wrapper class used to wrap ROOT random number engines with GSL interface

// @(#)root/mathmore:$Id$
// Author: L. Moneta Fri Aug 24 17:20:45 2007

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class GSLRngWrapper

#ifndef ROOT_Math_GSLRngROOTWrapper
#define ROOT_Math_GSLRngROOTWrapper

#include "gsl/gsl_rng.h"

#include <string>

namespace ROOT {

   namespace Math {

   /**
    * class for wrapping ROOT Engines in gsl_rng types which can be used as extra
    * GSL random number generators
    * For this we need to implment functions which will be called by gsl_rng.
    * The functions (Seed, Rndm, IntRndm) are passed in the gsl_rng_type and used to build a gsl_rng object.
    * When gsl_rng is alloacated, only the memory state is allocated using calloc(1,size), which gives a memory 
    * block of the given bytes and it initializes to zero. Therefore no constructor of GSLRngROOTWrapper can be called 
    * and also we cannot call non-static member funciton of the class. 
    * The underlined ROOT engine is then built and deleted using the functions CreateEngine() and FreeEngine(), 
    * called by the specific GSLRandomEngine class that instantiates for the  the generator (e.g. GSLRngMixMax)
    *
    **/
   template <class Engine>
   struct GSLRngROOTWrapper {

      Engine *fEngine = nullptr;

      // non need to have specific ctor/dtor since we cannot call them

      // function called by the specific GSLRndmEngine to create the ROOT engine
      static void CreateEngine(gsl_rng *r)
      {
         // the engine pointer is retrieved from r
         GSLRngROOTWrapper *wrng = ((GSLRngROOTWrapper *)r->state);
         //printf("calling create engine - engine pointer : %p\n", wrng->fEngine);
         if (!wrng->fEngine)
            wrng->fEngine = new Engine();
         // do nothing in case it is already created (e.g. when calling with default_seed)
      }

      static double Rndm(void *p) { return ((GSLRngROOTWrapper *)p)->fEngine->operator()(); }
      static unsigned long IntRndm(void *p) { return ((GSLRngROOTWrapper *)p)->fEngine->IntRndm(); }

      static void Seed(void *p, unsigned long seed)
      {
         auto wrng = ((GSLRngROOTWrapper *)p);
         // (GSL calls at the beginning with the defaul seed (typically zero))
         //printf("calling the seed function with %d on %p and engine %p\n", seed, p, wrng->fEngine);
         if (seed == gsl_rng_default_seed) {
            seed = 111; // avoid using 0 that for ROOT means a specific seed
            if (!wrng->fEngine) wrng->fEngine = new Engine();
         }
         assert(wrng->fEngine != nullptr);
         wrng->fEngine->SetSeed(seed);
      }
      static void FreeEngine(gsl_rng *r)
      {
         auto wrng = ((GSLRngROOTWrapper *)r->state);
         if (wrng->fEngine)
            delete wrng->fEngine;
         wrng->fEngine = nullptr;
      }

      static unsigned long Max() { return Engine::MaxInt(); }
      static unsigned long Min() { return Engine::MinInt(); }
      static size_t Size() { return sizeof(GSLRngROOTWrapper<Engine>); }
      static std::string Name() { return std::string("GSL_") + Engine::Name(); }
   };

   }  // end namespace Math
} // end namespace ROOT

#include "Math/MixMaxEngine.h"

// now define and implement the specific types    

typedef ROOT::Math::GSLRngROOTWrapper<ROOT::Math::MixMaxEngine<17,0>> GSLMixMaxWrapper;

static const std::string gsl_mixmax_name =  GSLMixMaxWrapper::Name();
static const gsl_rng_type mixmax_type =
{
   gsl_mixmax_name.c_str(), 
   GSLMixMaxWrapper::Max(), 
   GSLMixMaxWrapper::Min(),
   GSLMixMaxWrapper::Size(),
   &GSLMixMaxWrapper::Seed,
   &GSLMixMaxWrapper::IntRndm,
   &GSLMixMaxWrapper::Rndm
};

const gsl_rng_type *gsl_rng_mixmax = &mixmax_type;

#endif

