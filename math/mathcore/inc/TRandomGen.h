// @(#)root/mathcore:$Id$
// Author: Rene Brun   04/03/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TRandomGen
#define ROOT_TRandomGen



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRandomGen
// @ingroup Random
//                                                                      //
// Generic random number generator class which is template on the type  //
// of engine. Using this class different random number generator all    //
// implementing the TRandom interface can be built.                     //
// The available random number engine that can be presently used are    //
//    * ROOT::Math::MixMaxEngine to create random number generators     //
//       based on the MIXMAX family of generators. Different generators //
//       can be created for different state N.                          //
//    * ROOT::MATH::StdEngine to create genersators based on engines    //
//      provided by the C++ standard libraries
//
//  Convenient typedef are defines to define the different types of
//  generators. These typedef are
//   * TRandomMixMax for the MixMaxEngine<240,0>  (MIXMAX with state N=240)
//   * TRandomMixMax17 for the MixMaxEngine<17,0>  (MIXMAX with state N=17)
//   * TRandomMixMax256 for the MixMaxEngine<256,2> (MIXMAX with state N=256 )
//   * TRandomMT64 for the  StdEngine<std::mt19937_64> ( MersenneTwister 64 bits)
//   * TRandomRanlux48 for the  StdEngine<std::ranlux48> (Ranlux 48 bits)
//
//                                                                     //
//////////////////////////////////////////////////////////////////////////

#include "TRandom.h"

#include <string>

template<class Engine>
class TRandomGen : public TRandom {

protected:

   Engine  fEngine;   // random number generator engine
public:
   
   TRandomGen(ULong_t seed=1) {
      fEngine.SetSeed(seed);
      SetName(TString::Format("Random_%s", std::string(fEngine.Name()).c_str()));
      SetTitle(TString::Format("Random number generator: %s", std::string(fEngine.Name()).c_str()));
   }
   virtual ~TRandomGen() {}
   using TRandom::Rndm; 
   virtual  Double_t Rndm( ) { return fEngine(); }
   virtual  void     RndmArray(Int_t n, Float_t *array) {
      for (int i = 0; i < n; ++i) array[i] = fEngine(); 
   }
   virtual  void     RndmArray(Int_t n, Double_t *array) {
            for (int i = 0; i < n; ++i) array[i] = fEngine(); 
   }
   virtual  void     SetSeed(ULong_t seed=0) {
      fEngine.SetSeed(seed);
   }

   ClassDef(TRandomGen,1)  //Generic Random number generator template on the Engine type
};

// some useful typedef
#include "Math/StdEngine.h"
#include "Math/MixMaxEngine.h"
#include "Math/RanluxppEngine.h"

// not working wight now for this classes
//#define  DEFINE_TEMPL_INSTANCE
#ifdef DEFINE_TEMPL_INSTANCE

extern template class  TRandomGen<ROOT::Math::MixMaxEngine<240,0>>;
extern template class TRandomGen<ROOT::Math::MixMaxEngine<256,2>>;
extern template class TRandomGen<ROOT::Math::MixMaxEngine<256,4>>;
extern template class TRandomGen<ROOT::Math::MixMaxEngine<17,0>>;
extern template class TRandomGen<ROOT::Math::MixMaxEngine<17,1>>;

extern template class TRandomGen<ROOT::Math::RanluxppEngine2048>;

extern template class  TRandomGen<ROOT::Math::StdEngine<std::mt19937_64> >;
extern template class  TRandomGen<ROOT::Math::StdEngine<std::ranlux48> >;

#endif
/**
  @ingroup Random
  MIXMAX generator based on a state of N=240.  
  This generator is described in this paper:

   K. Savvidy and G. Savvidy, *Spectrum and Entropy of C-systems. MIXMAX random number generator*,
  Chaos, Solitons & Fractals, Volume 91, (2016) pp. 33–38 http://dx.doi.org/10.1016/j.chaos.2016.05.003
 */
typedef TRandomGen<ROOT::Math::MixMaxEngine<240,0>> TRandomMixMax;

/**
  @ingroup Random
  MIXMAX generator based on a state of N=17. This generator has a fast seeding time
  compared to N=240. 
   This generator is described in this paper:

   K. Savvidy and G. Savvidy, *Spectrum and Entropy of C-systems. MIXMAX random number generator*,
  Chaos, Solitons & Fractals, Volume 91, (2016) pp. 33–38 http://dx.doi.org/10.1016/j.chaos.2016.05.003
 */
typedef TRandomGen<ROOT::Math::MixMaxEngine<17,0>> TRandomMixMax17;

/**
  @ingroup Random
  MIXMAX generator based on a state of N=256, based on the generator descrived in this 
  paper:

   K. Savvidy, *The MIXMAX random number generator*, Comp. Phys. Commun. 196 (2015), pp 161–165
  http://dx.doi.org/10.1016/j.cpc.2015.06.003

   This generator has been implemented with a skipping value of 2 iterations (so retaining one 
   matrix iteration every 3). 
  
 */
typedef TRandomGen<ROOT::Math::MixMaxEngine<256,2>> TRandomMixMax256;

typedef TRandomGen<ROOT::Math::RanluxppEngine2048> TRandomRanluxpp;

/**
  @ingroup Random
  Generator based on a the Mersenne-Twister generator with 64 bits, 
  using the implementation provided by the standard library, 
  std::mt19937_64 (see  http://www.cplusplus.com/reference/random/mt19937_64/ )

 */
typedef TRandomGen<ROOT::Math::StdEngine<std::mt19937_64> > TRandomMT64;
/**
  @ingroup Random
  Generator based on a the RanLux generator with 48 bits, 
  using the implementation provided by the standard library, 
  std::ranlux48 (see  http://www.cplusplus.com/reference/random/ranlux48/ )

 */
typedef TRandomGen<ROOT::Math::StdEngine<std::ranlux48> > TRandomRanlux48;


#endif
