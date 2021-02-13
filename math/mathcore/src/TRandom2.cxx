// @(#)root/mathcore:$Id$
// Author: Rene Brun, Lorenzo Moneta  17/05/2006

/**

\class TRandom2

Random number generator class based on the maximally quidistributed combined
Tausworthe generator by L'Ecuyer.

The period of the generator is 2**88 (about 10**26) and it uses only 3 words
for the state.

For more information see:
P. L'Ecuyer, Mathematics of Computation, 65, 213 (1996)
P. L'Ecuyer, Mathematics of Computation, 68, 225 (1999)

The publication are available online at
 [http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme.ps]
 [http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme2.ps]

@ingroup Random

*/

#include "TRandom2.h"
#include "TRandom3.h"
#include "TUUID.h"


ClassImp(TRandom2);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TRandom2::TRandom2(UInt_t seed)
{
   SetName("Random2");
   SetTitle("Random number generator with period of about  10**26");
   SetSeed(seed);
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor

TRandom2::~TRandom2()
{
}

////////////////////////////////////////////////////////////////////////////////
///  TausWorth generator from L'Ecuyer, uses as seed 3x32bits integers
///  Use a mask of 0xffffffffUL to make in work on 64 bit machines
///  Periodicity of about  10**26
///  Generate number in interval (0,1)  : 0 and 1 are not included in the interval

Double_t TRandom2::Rndm()
{
#define TAUSWORTHE(s,a,b,c,d) (((s &c) <<d) & 0xffffffffUL ) ^ ((((s <<a) & 0xffffffffUL )^s) >>b)

   // scale by 1./(Max<UINT> + 1) = 1./4294967296
   const double kScale = 2.3283064365386963e-10;    // range in 32 bit ( 1/(2**32)

   fSeed  = TAUSWORTHE (fSeed, 13, 19, 4294967294UL, 12);
   fSeed1 = TAUSWORTHE (fSeed1, 2, 25, 4294967288UL, 4);
   fSeed2 = TAUSWORTHE (fSeed2, 3, 11, 4294967280UL, 17);

   UInt_t iy = fSeed ^ fSeed1 ^ fSeed2;
   if (iy) return  kScale*static_cast<Double_t>(iy);
   return Rndm();
}

////////////////////////////////////////////////////////////////////////////////
/// Return an array of n random numbers uniformly distributed in ]0,1]

void TRandom2::RndmArray(Int_t n, Float_t *array)
{
   const double kScale = 2.3283064365386963e-10;    // range in 32 bit ( 1/(2**32)

   UInt_t iy;

   for(Int_t i=0; i<n; i++) {
      fSeed  = TAUSWORTHE (fSeed, 13, 19, 4294967294UL, 12);
      fSeed1 = TAUSWORTHE (fSeed1, 2, 25, 4294967288UL, 4);
      fSeed2 = TAUSWORTHE (fSeed2, 3, 11, 4294967280UL, 17);

      iy = fSeed ^ fSeed1 ^ fSeed2;
      if (iy) array[i] = (Float_t)(kScale*static_cast<Double_t>(iy));
      else    array[i] = Rndm();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return an array of n random numbers uniformly distributed in ]0,1]

void TRandom2::RndmArray(Int_t n, Double_t *array)
{
   const double kScale = 2.3283064365386963e-10;    // range in 32 bit ( 1/(2**32)

   UInt_t iy;
   for(Int_t i=0; i<n; i++) {
      fSeed  = TAUSWORTHE (fSeed, 13, 19, 4294967294UL, 12);
      fSeed1 = TAUSWORTHE (fSeed1, 2, 25, 4294967288UL, 4);
      fSeed2 = TAUSWORTHE (fSeed2, 3, 11, 4294967280UL, 17);

      iy = fSeed ^ fSeed1 ^ fSeed2;
      if (iy) array[i] = kScale*static_cast<Double_t>(iy);
      else    array[i] = Rndm();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the generator seed.
/// If the seed given is zero, generate automatically seed values which
/// are different every time by using TRandom3  and TUUID
/// If a seed is given generate the other two needed for the generator state using
/// a linear congruential generator
/// The only condition, stated at the end of the 1999 L'Ecuyer paper is that the seeds
/// must be greater than 1,7 and 15.
/// Note that after seeting the seed the generator is warmed up by calling it internally few
/// times therefore the returned seed in TRandom2::GetSeed will be a different value.
/// For this generator the user will have to store the provided seed externally
/// if he wants to reproduce the random sequence
void TRandom2::SetSeed(ULong_t seed)
{
#define LCG(n) ((69069 * n) & 0xffffffffUL)  // linear congurential generator

   if (seed > 0) {
      fSeed = LCG (seed);
      if (fSeed < 2) fSeed += 2UL;
      fSeed1 = LCG (fSeed);
      if (fSeed1 < 8) fSeed1 += 8UL;
      fSeed2 = LCG (fSeed1);
      if (fSeed2 < 16) fSeed2 += 16UL;
   } else {
      // initialize using a TUUID
      TUUID u;
      UChar_t uuid[16];
      u.GetUUID(uuid);
      fSeed  =  UInt_t(uuid[3])*16777216 + UInt_t(uuid[2])*65536 + UInt_t(uuid[1])*256 + UInt_t(uuid[0]);
      fSeed1  =  UInt_t(uuid[7])*16777216 + UInt_t(uuid[6])*65536 + UInt_t(uuid[5])*256 + UInt_t(uuid[4]);
      fSeed2  =  UInt_t(uuid[11])*16777216 + UInt_t(uuid[10])*65536 + UInt_t(uuid[9])*256 + UInt_t(uuid[8]);
      // use also the other bytes
      UInt_t seed3 = UInt_t(uuid[15])*16777216 + UInt_t(uuid[14])*65536 + UInt_t(uuid[13])*256 + UInt_t(uuid[12]);
      fSeed2 += seed3;


      //    TRandom r3(0);
      // fSeed   = static_cast<UInt_t> (4294967296.*r3.Rndm());
      // fSeed1  = static_cast<UInt_t> (4294967296.*r3.Rndm());
      // fSeed2  = static_cast<UInt_t> (4294967296.*r3.Rndm());

      if (fSeed < 2)   fSeed  += 2UL;
      if (fSeed1 < 8)  fSeed1 += 8UL;
      if (fSeed2 < 16) fSeed2 += 16UL;
   }

   // "warm it up" by calling it 6 times
   for (int i = 0; i < 6; ++i)
      Rndm();

   return;
}
