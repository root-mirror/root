// @(#)root/mathcore:$Id$
// Author: L. Moneta Tue Aug 4 2015

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2015  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// random engines based on ROOT 

#ifndef ROOT_Math_StdEngine
#define ROOT_Math_StdEngine

#include <random>

namespace ROOT {

   namespace Math {


      class StdRandomEngine {};

      template<class Generator>
      struct StdEngineType {
         static const char *  Name() { return "std_random_eng";}
      };
      template<>
         struct StdEngineType<std::minstd_rand> {
         static const char *  Name() { return "std_minstd_rand";}
      };
      template<>
      struct StdEngineType<std::mt19937> {
         static const char *  Name() { return "std_mt19937";}
      };
      template<>
      struct StdEngineType<std::mt19937_64> {
         static const char *  Name() { return "std_mt19937_64";}
      };
      template<>
      struct StdEngineType<std::ranlux24> {
         static const char *  Name() { return "std_ranlux24";}
      };
      template<>
      struct StdEngineType<std::ranlux48> {
         static const char *  Name() { return "std_ranlux48";}
      };
      template<>
      struct StdEngineType<std::knuth_b> {
         static const char *  Name() { return "std_knuth_b";}
      };
      template<>
      struct StdEngineType<std::random_device> {
         static const char *  Name() { return "std_random_device";}
      };
      
      
      /** 
          @ingroup Random
          Class to wrap engines fron the C++ standard random library in 
          the ROOT Random interface. 
          This casess is then by used by the generic TRandoGen class 
          to provide TRandom interrace generators for the C++ random generators.

          See for examples the TRandomMT64 and TRandomRanlux48 generators 
          which are typede's to TRandomGen instaniated with some 
          random engine from the C++ standard library. 

      */

      template <class Generator> 
      class StdEngine {


      public:

         typedef  StdRandomEngine BaseType; 
         typedef  typename Generator::result_type Result_t;

         StdEngine() : fGen() {
            fCONS = 1./fGen.max(); 
         }


         void SetSeed(Result_t seed) { fGen.seed(seed);}

         double Rndm() {
            Result_t rndm = fGen(); // generate integer number according to the type
            if (rndm != 0) return  fCONS*rndm;
            return Rndm();
         }

         Result_t IntRndm() {
            return fGen();
         }

         double operator() () {
            return Rndm(); 
         }

         static const char * Name()  {
            return StdEngineType<Generator>::Name(); 
         }

         static uint64_t MaxInt() { return Generator::max(); }


      private:
         Generator fGen;
         double fCONS;   //! cached value of maximum integer value generated
      };


      extern template class StdEngine<std::mt19937_64>;
      extern template class StdEngine<std::ranlux48>;

   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_StdEngine */
