// @(#)root/mathcore:$Id$
// Author: Peter Malzacher   31/08/99

/**

\class TRandom3

Random number generator class based on
  M. Matsumoto and T. Nishimura,
  Mersenne Twister: A 623-diminsionally equidistributed
  uniform pseudorandom number generator
  ACM Transactions on Modeling and Computer Simulation,
  Vol. 8, No. 1, January 1998, pp 3--30.

For more information see the Mersenne Twister homepage
  [http://www.math.keio.ac.jp/~matumoto/emt.html]

Advantage:

-  large period 2**19937 -1
-  relativly fast (slightly slower than TRandom2 but much faster than TRandom1)

Drawback:
-  a relative large internal state of 624 integers 
- generate only 32 random bits 
- not passing all the random generator tests. It fails some tests in TestU01
 (see [http://simul.iro.umontreal.ca/testu01/tu01.html])

An altenativly excellent generator passing all tests of TestU01, having 61 random bits and 
fast as Mersenne and Twister is MIXMAX (TRandomMixMax). 

@ingroup Random

*/

//////////////////////////////////////////////////////////////////////
// Aug.99 ROOT implementation based on CLHEP by P.Malzacher
//
// the original code contains the following copyright notice:
/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */
/* 02111-1307  USA                                                 */
/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */
/////////////////////////////////////////////////////////////////////

#include "TRandom3.h"
#include "TBuffer.h"
#include "TRandom2.h"
#include "TUUID.h"

TRandom *gRandom = new TRandom3();
#ifdef R__COMPLETE_MEM_TERMINATION
namespace {
   struct TRandomCleanup {
      ~TRandomCleanup() { delete gRandom; gRandom = nullptr; }
   };
   static TRandomCleanup gCleanupRandom;
}
#endif

ClassImp(TRandom3);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor
/// If seed is 0, the seed is automatically computed via a TUUID object.
/// In this case the seed is guaranteed to be unique in space and time.

TRandom3::TRandom3(UInt_t seed)
{
   SetName("Random3");
   SetTitle("Random number generator: Mersenne Twister");
   SetSeed(seed);
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor

TRandom3::~TRandom3()
{
}

////////////////////////////////////////////////////////////////////////////////
///  Machine independent random number generator.
///  Produces uniformly-distributed floating points in (0,1)
///  Method: Mersenne Twister

Double_t TRandom3::Rndm()
{
   UInt_t y;

   const Int_t  kM = 397;
   const Int_t  kN = 624;
   const UInt_t kTemperingMaskB =  0x9d2c5680;
   const UInt_t kTemperingMaskC =  0xefc60000;
   const UInt_t kUpperMask =       0x80000000;
   const UInt_t kLowerMask =       0x7fffffff;
   const UInt_t kMatrixA =         0x9908b0df;

   if (fCount624 >= kN) {
      Int_t i;

      for (i=0; i < kN-kM; i++) {
         y = (fMt[i] & kUpperMask) | (fMt[i+1] & kLowerMask);
         fMt[i] = fMt[i+kM] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
      }

      for (   ; i < kN-1    ; i++) {
         y = (fMt[i] & kUpperMask) | (fMt[i+1] & kLowerMask);
         fMt[i] = fMt[i+kM-kN] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
      }

      y = (fMt[kN-1] & kUpperMask) | (fMt[0] & kLowerMask);
      fMt[kN-1] = fMt[kM-1] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
      fCount624 = 0;
   }

   y = fMt[fCount624++];
   y ^=  (y >> 11);
   y ^= ((y << 7 ) & kTemperingMaskB );
   y ^= ((y << 15) & kTemperingMaskC );
   y ^=  (y >> 18);

   // 2.3283064365386963e-10 == 1./(max<UINt_t>+1)  -> then returned value cannot be = 1.0
   if (y) return ( (Double_t) y * 2.3283064365386963e-10); // * Power(2,-32)
   return Rndm();
}

////////////////////////////////////////////////////////////////////////////////
/// Return an array of n random numbers uniformly distributed in ]0,1]

void TRandom3::RndmArray(Int_t n, Float_t *array)
{
  for(Int_t i=0; i<n; i++) array[i]=(Float_t)Rndm();
}

////////////////////////////////////////////////////////////////////////////////
/// Return an array of n random numbers uniformly distributed in ]0,1]

void TRandom3::RndmArray(Int_t n, Double_t *array)
{
   Int_t k = 0;

   UInt_t y;

   const Int_t  kM = 397;
   const Int_t  kN = 624;
   const UInt_t kTemperingMaskB =  0x9d2c5680;
   const UInt_t kTemperingMaskC =  0xefc60000;
   const UInt_t kUpperMask =       0x80000000;
   const UInt_t kLowerMask =       0x7fffffff;
   const UInt_t kMatrixA =         0x9908b0df;

   while (k < n) {
      if (fCount624 >= kN) {
         Int_t i;

         for (i=0; i < kN-kM; i++) {
            y = (fMt[i] & kUpperMask) | (fMt[i+1] & kLowerMask);
            fMt[i] = fMt[i+kM] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
         }

         for (   ; i < kN-1    ; i++) {
            y = (fMt[i] & kUpperMask) | (fMt[i+1] & kLowerMask);
            fMt[i] = fMt[i+kM-kN] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
         }

         y = (fMt[kN-1] & kUpperMask) | (fMt[0] & kLowerMask);
         fMt[kN-1] = fMt[kM-1] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
         fCount624 = 0;
      }

      y = fMt[fCount624++];
      y ^=  (y >> 11);
      y ^= ((y << 7 ) & kTemperingMaskB );
      y ^= ((y << 15) & kTemperingMaskC );
      y ^=  (y >> 18);

      if (y) {
         array[k] = Double_t( y * 2.3283064365386963e-10); // * Power(2,-32)
         k++;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///  Set the random generator sequence
/// if seed is 0 (default value) a TUUID is generated and used to fill
/// the first 8 integers of the seed array.
/// In this case the seed is guaranteed to be unique in space and time.
/// Use upgraded seeding procedure to fix a known problem when seeding with values
/// with many zero in the bit pattern (like 2**28).
/// see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html

void TRandom3::SetSeed(ULong_t seed)
{
   TRandom::SetSeed(seed);
   fCount624 = 624;
   if (seed > 0) {
      fMt[0] = fSeed;

      // use multipliers from  Knuth's "Art of Computer Programming" Vol. 2, 3rd Ed. p.106
      for(Int_t i=1; i<624; i++) {
         fMt[i] = (1812433253 * ( fMt[i-1]  ^ ( fMt[i-1] >> 30)) + i );
      }

   } else {

      // use TRandom2 (which is based on TUUId to generate the seed
      // TRandom2 works fairly well  and has been tested against example
      // layout in https://savannah.cern.ch/bugs/?99516
      TRandom2 r(0);
      for (Int_t i = 0; i< 624; i++) {
         fMt[i]   = static_cast<UInt_t> (4294967296.*r.Rndm());
      }
      // warm up the generator calling it 10 times
      for (Int_t i = 0; i < 10; ++i) Rndm();
   }


}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TRandom3.

void TRandom3::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TRandom3::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TRandom::Streamer(R__b);
      R__b.ReadStaticArray(fMt);
      R__b >> fCount624;
      R__b.CheckByteCount(R__s, R__c, TRandom3::IsA());
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TRandom3::Class(),this);
   }
}
