// @(#)root/mathcore:$Id$
// Author: Rene Brun, Lorenzo Moneta   15/12/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**

\class TRandom

@ingroup Random

This is the base class for the ROOT Random number generators.
This class defines the ROOT Random number interface and it should not be instantiated directly but used via its derived
classes. The generator provided in TRandom itself is a LCG (Linear Congruential Generator), the <a
href="https://www.gnu.org/software/gsl/manual/html_node/Unix-random-number-generators.html">BSD `rand` generator</a>,
that it should not be used because its period is only 2**31, i.e. approximatly 2 billion events, that can be generated
in just few seconds.

To generate random numbers, one should use the derived class, which  are :
- TRandom3: it is based on the "Mersenne Twister generator",
it is fast and a very long period of about \f$10^{6000}\f$. However it fails some of the most stringent tests of the
<a href="http://simul.iro.umontreal.ca/testu01/tu01.html">TestU01 suite</a>.
In addition this generator provide only numbers with 32 random bits, which might be not sufficient for some application
based on double or extended precision. This generator is however used in ROOT used to instantiate the global pointer to
the ROOT generator, *gRandom*.
- ::TRandomRanluxpp : New implementation of the Ranlux generator algorithm based on a fast modular multiplication of
576 bits. This new implementation is built on the idea and the original code of Alexei Sibidanov, described in his
<a href="https://arxiv.org/abs/1705.03123">paper </a>. It generates random numbers with 52 bit precision (double
precision) and it has an higher luxury level than the original Ranlux generator (`p = 2048` instead of `p=794`).
- ::TRandomMixMax: Generator based on the family of the MIXMAX matrix generators (see the
<a href="https://mixmax.hepforge.org">MIXMAX HEPFORGE Web page</a> and the the documentation of the class
ROOT::Math::MixMaxEngine for more information), that are base on the Asanov dynamical C systems. This generator has a
state of N=240 64 bit integers, proof random properties, it provides 61 random bits and it has a very large period
(\f$10^{4839}\f$). Furthermore, it provides the capability to be seeded with the guarantee that, for each given
different seed, a different sequence of random numbers will be generated. The only drawback is that the seeding time is
time consuming, of the order of 0.1 ms, while the time to generate a number is few ns (more than 10000 faster).
- ::TRandomMixMax17: Another MixMax generator, but with a smaller state, N=17, and this results in a smaller entropy
than the generator with N=240. However, it has the same seeding capabilities, with a much faster seeding time (about 200
times less than TRandomMixMax240 and comparable to TRandom3).
- ::TRandomMixMax256 : A variant of the MIXMAX generators, based on a state of N=256, and described in the
        <a  href="http://arxiv.org/abs/1403.5355">2015 paper</a>. This implementation has been modified with respect to
the paper, by skipping 2 internal interations, to provide improved random properties.
- ::TRandomMT64 :  Generator based on a the Mersenne-Twister generator with 64 bits,
  using the implementation provided by the standard library ( <a
href="http://www.cplusplus.com/reference/random/mt19937_64/">std::mt19937_64</a> )
- TRandom1  based on the RANLUX algorithm, has mathematically proven random proprieties
  and a period of about \f$10{171}\f$. It is however much slower than the others and it has only 24 random bits. It can
be constructed with different luxury levels.
- ::TRandomRanlux48 : Generator based on a the RanLux generator with 48 bits and highest luxury level
  using the implementation provided by the standard library (<a
href="http://www.cplusplus.com/reference/random/ranlux48/">std::ranlux48</a>). The drawback of this generator is its
slow generation time.
- TRandom2  is based on the Tausworthe generator of L'Ecuyer, and it has the advantage
of being fast and using only 3 words (of 32 bits) for the state. The period however is not impressively long, it is
10**26.

Using the template TRandomGen class (template on the contained Engine type), it is possible to add any generator based
on the standard C++ random library (see the C++ <a href="http://www.cplusplus.com/reference/random/">random</a>
documentation.) or different variants of the MIXMAX generator using the ROOT::Math::MixMaxEngine. Some of the listed
generator above (e.g. TRandomMixMax256 or TRandomMT64) are convenient typedef's of generator built using the template
TRandomGen class.

Please note also that this class (TRandom) implements also a very simple generator (linear congruential) with period =
\f$10^9\f$, known to have defects (the lower random bits are correlated) and it is failing the majority of the random
number generator tests. Therefore it should NOT be used in any statistical study.

The following table shows some timings (in nanoseconds/call)
for the random numbers obtained using a macbookpro 2.6 GHz Intel Core i7 CPU:


-   TRandom            3   ns/call     (but this is a very BAD Generator, not to be used)
-   TRandom2           5   ns/call
-   TRandom3           5   ns/call
-   ::TRandomMixMax      6   ns/call
-   ::TRandomMixMax17    6   ns/call
-   ::TRandomMT64        9   ns/call
-   ::TRandomMixMax256  10   ns/call
-   ::TRandomRanluxpp   14   ns/call
-   ::TRandom1          80   ns/call
-   ::TRandomRanlux48  250  ns/call

The following methods are provided to generate random numbers disctributed according to some basic distributions:

- Exp(Double_t tau)
- Integer(UInt_t imax)
- Gaus(Double_t mean, Double_t sigma)
- Rndm()
- Uniform(Double_t)
- Landau(Double_t mean, Double_t sigma)
- Poisson(Double_t mean)
- Binomial(Int_t ntot, Double_t prob)

Random numbers distributed according to 1-d, 2-d or 3-d distributions contained in TF1, TF2 or TF3 objects can also be
generated. For example, to get a random number distributed following abs(sin(x)/x)*sqrt(x) you can do : \code{.cpp} TF1
*f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,10); double r = f1->GetRandom(); \endcode or you can use the UNURAN
package. You need in this case to initialize UNURAN to the function you would like to generate. \code{.cpp} TUnuran u;
  u.Init(TUnuranDistrCont(f1));
  double r = u.Sample();
\endcode

The techniques of using directly a TF1,2 or 3 function is powerful and
can be used to generate numbers in the defined range of the function.
Getting a number from a TF1,2,3 function is also quite fast.
UNURAN is a  powerful and flexible tool which containes various methods for
generate random numbers for continuous distributions of one and multi-dimension.
It requires some set-up (initialization) phase and can be very fast when the distribution
parameters are not changed for every call.

The following table shows some timings (in nanosecond/call)
for basic functions,  TF1 functions and using UNURAN obtained running
the tutorial math/testrandom.C
Numbers have been obtained on an Intel Xeon Quad-core Harpertown (E5410) 2.33 GHz running
Linux SLC4 64 bit and compiled with gcc 3.4

~~~~
Distribution            nanoseconds/call
                    TRandom  TRandom1 TRandom2 TRandom3
Rndm..............    5.000  105.000    7.000   10.000
RndmArray.........    4.000  104.000    6.000    9.000
Gaus..............   36.000  180.000   40.000   48.000
Rannor............  118.000  220.000  120.000  124.000
Landau............   22.000  123.000   26.000   31.000
Exponential.......   93.000  198.000   98.000  104.000
Binomial(5,0.5)...   30.000  548.000   46.000   65.000
Binomial(15,0.5)..   75.000 1615.000  125.000  178.000
Poisson(3)........   96.000  494.000  109.000  125.000
Poisson(10).......  138.000 1236.000  165.000  203.000
Poisson(70).......  818.000 1195.000  835.000  844.000
Poisson(100)......  837.000 1218.000  849.000  864.000
GausTF1...........   83.000  180.000   87.000   88.000
LandauTF1.........   80.000  180.000   83.000   86.000
GausUNURAN........   40.000  139.000   41.000   44.000
PoissonUNURAN(10).   85.000  271.000   92.000  102.000
PoissonUNURAN(100)   62.000  256.000   69.000   78.000
~~~~

Note that the time to generate a number from an arbitrary TF1 function
using TF1::GetRandom or using TUnuran is  independent of the complexity of the function.

TH1::FillRandom(TH1 *) or TH1::FillRandom(const char *tf1name)
can be used to fill an histogram (1-d, 2-d, 3-d from an existing histogram
or from an existing function.

Note this interesting feature when working with objects.
 You can use several TRandom objects, each with their "independent"
 random sequence. For example, one can imagine
~~~~
    TRandom *eventGenerator = new TRandom();
    TRandom *tracking       = new TRandom();
~~~~
 `eventGenerator` can be used to generate the event kinematics.
 tracking can be used to track the generated particles with random numbers
 independent from eventGenerator.
 This very interesting feature gives the possibility to work with simple
 and very fast random number generators without worrying about
 random number periodicity as it was the case with Fortran.
 One can use TRandom::SetSeed to modify the seed of one generator.

A TRandom object may be written to a Root file

- as part of another object
- or with its own key (example: `gRandom->Write("Random")` )  ;

*/

#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "Math/QuantFuncMathCore.h"
#include "TUUID.h"

ClassImp(TRandom);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor. For seed see SetSeed().

TRandom::TRandom(UInt_t seed): TNamed("Random","Default Random number generator")
{
   SetSeed(seed);
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor. Can reset gRandom to 0 if gRandom points to this
/// generator.

TRandom::~TRandom()
{
   if (gRandom == this) gRandom = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Generates a random integer N according to the binomial law.
/// Coded from Los Alamos report LA-5061-MS.
///
/// N is binomially distributed between 0 and ntot inclusive
/// with mean prob*ntot and prob is between 0 and 1.
///
/// Note: This function should not be used when ntot is large (say >100).
/// The normal approximation is then recommended instead
/// (with mean =*ntot+0.5 and standard deviation sqrt(ntot*prob*(1-prob)).

Int_t TRandom::Binomial(Int_t ntot, Double_t prob)
{
   if (prob < 0 || prob > 1) return 0;
   Int_t n = 0;
   for (Int_t i=0;i<ntot;i++) {
      if (Rndm() > prob) continue;
      n++;
   }
   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a number distributed following a BreitWigner function with mean and gamma.

Double_t TRandom::BreitWigner(Double_t mean, Double_t gamma)
{
   Double_t rval, displ;
   rval = 2*Rndm() - 1;
   displ = 0.5*gamma*TMath::Tan(rval*TMath::PiOver2());

   return (mean+displ);
}

////////////////////////////////////////////////////////////////////////////////
/// Generates random vectors, uniformly distributed over a circle of given radius.
///   Input : r = circle radius
///   Output: x,y a random 2-d vector of length r

void TRandom::Circle(Double_t &x, Double_t &y, Double_t r)
{
   Double_t phi = Uniform(0,TMath::TwoPi());
   x = r*TMath::Cos(phi);
   y = r*TMath::Sin(phi);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns an exponential deviate.
///
///          exp( -t/tau )

Double_t TRandom::Exp(Double_t tau)
{
   Double_t x = Rndm();              // uniform on ] 0, 1 ]
   Double_t t = -tau * TMath::Log( x ); // convert to exponential distribution
   return t;
}

////////////////////////////////////////////////////////////////////////////////
/// Samples a random number from the standard Normal (Gaussian) Distribution
/// with the given mean and sigma.
/// Uses the Acceptance-complement ratio from W. Hoermann and G. Derflinger
/// This is one of the fastest existing method for generating normal random variables.
/// It is a factor 2/3 faster than the polar (Box-Muller) method used in the previous
/// version of TRandom::Gaus. The speed is comparable to the Ziggurat method (from Marsaglia)
/// implemented for example in GSL and available in the MathMore library.
///
/// REFERENCE:  - W. Hoermann and G. Derflinger (1990):
///              The ACR Method for generating normal random variables,
///              OR Spektrum 12 (1990), 181-185.
///
/// Implementation taken from
/// UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien

Double_t TRandom::Gaus(Double_t mean, Double_t sigma)
{
   const Double_t kC1 = 1.448242853;
   const Double_t kC2 = 3.307147487;
   const Double_t kC3 = 1.46754004;
   const Double_t kD1 = 1.036467755;
   const Double_t kD2 = 5.295844968;
   const Double_t kD3 = 3.631288474;
   const Double_t kHm = 0.483941449;
   const Double_t kZm = 0.107981933;
   const Double_t kHp = 4.132731354;
   const Double_t kZp = 18.52161694;
   const Double_t kPhln = 0.4515827053;
   const Double_t kHm1 = 0.516058551;
   const Double_t kHp1 = 3.132731354;
   const Double_t kHzm = 0.375959516;
   const Double_t kHzmp = 0.591923442;
   /*zhm 0.967882898*/

   const Double_t kAs = 0.8853395638;
   const Double_t kBs = 0.2452635696;
   const Double_t kCs = 0.2770276848;
   const Double_t kB  = 0.5029324303;
   const Double_t kX0 = 0.4571828819;
   const Double_t kYm = 0.187308492 ;
   const Double_t kS  = 0.7270572718 ;
   const Double_t kT  = 0.03895759111;

   Double_t result;
   Double_t rn,x,y,z;

   do {
      y = Rndm();

      if (y>kHm1) {
         result = kHp*y-kHp1; break; }

      else if (y<kZm) {
         rn = kZp*y-1;
         result = (rn>0) ? (1+rn) : (-1+rn);
         break;
      }

      else if (y<kHm) {
         rn = Rndm();
         rn = rn-1+rn;
         z = (rn>0) ? 2-rn : -2-rn;
         if ((kC1-y)*(kC3+TMath::Abs(z))<kC2) {
            result = z; break; }
         else {
            x = rn*rn;
            if ((y+kD1)*(kD3+x)<kD2) {
               result = rn; break; }
            else if (kHzmp-y<exp(-(z*z+kPhln)/2)) {
               result = z; break; }
            else if (y+kHzm<exp(-(x+kPhln)/2)) {
               result = rn; break; }
         }
      }

      while (1) {
         x = Rndm();
         y = kYm * Rndm();
         z = kX0 - kS*x - y;
         if (z>0)
            rn = 2+y/x;
         else {
            x = 1-x;
            y = kYm-y;
            rn = -(2+y/x);
         }
         if ((y-kAs+x)*(kCs+x)+kBs<0) {
            result = rn; break; }
         else if (y<x+kT)
            if (rn*rn<4*(kB-log(x))) {
               result = rn; break; }
      }
   } while(0);

   return mean + sigma * result;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a random integer uniformly distributed on the interval [ 0, imax-1 ].
/// Note that the interval contains the values of 0 and imax-1 but not imax.

UInt_t TRandom::Integer(UInt_t imax)
{
   UInt_t ui;
   ui = (UInt_t)(imax*Rndm());
   return ui;
}

////////////////////////////////////////////////////////////////////////////////
/// Generate a random number following a Landau distribution
/// with location parameter mu and scale parameter sigma:
///      Landau( (x-mu)/sigma )
/// Note that mu is not the mpv(most probable value) of the Landa distribution
/// and sigma is not the standard deviation of the distribution which is not defined.
/// For mu  =0 and sigma=1, the mpv = -0.22278
///
/// The Landau random number generation is implemented using the
/// function landau_quantile(x,sigma), which provides
/// the inverse of the landau cumulative distribution.
/// landau_quantile has been converted from CERNLIB ranlan(G110).

Double_t TRandom::Landau(Double_t mu, Double_t sigma)
{
   if (sigma <= 0) return 0;
   Double_t x = Rndm();
   Double_t res = mu + ROOT::Math::landau_quantile(x, sigma);
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Generates a random integer N according to a Poisson law.
/// Prob(N) = exp(-mean)*mean^N/Factorial(N)
///
/// Use a different procedure according to the mean value.
/// The algorithm is the same used by CLHEP.
/// For lower value (mean < 25) use the rejection method based on
/// the exponential.
/// For higher values use a rejection method comparing with a Lorentzian
/// distribution, as suggested by several authors.
/// This routine since is returning 32 bits integer will not work for values
/// larger than 2*10**9.
/// One should then use the Trandom::PoissonD for such large values.

Int_t TRandom::Poisson(Double_t mean)
{
   Int_t n;
   if (mean <= 0) return 0;
   if (mean < 25) {
      Double_t expmean = TMath::Exp(-mean);
      Double_t pir = 1;
      n = -1;
      while(1) {
         n++;
         pir *= Rndm();
         if (pir <= expmean) break;
      }
      return n;
   }
   // for large value we use inversion method
   else if (mean < 1E9) {
      Double_t em, t, y;
      Double_t sq, alxm, g;
      Double_t pi = TMath::Pi();

      sq = TMath::Sqrt(2.0*mean);
      alxm = TMath::Log(mean);
      g = mean*alxm - TMath::LnGamma(mean + 1.0);

      do {
         do {
            y = TMath::Tan(pi*Rndm());
            em = sq*y + mean;
         } while( em < 0.0 );

         em = TMath::Floor(em);
         t = 0.9*(1.0 + y*y)* TMath::Exp(em*alxm - TMath::LnGamma(em + 1.0) - g);
      } while( Rndm() > t );

      return static_cast<Int_t> (em);

   }
   else {
      // use Gaussian approximation vor very large values
      n = Int_t(Gaus(0,1)*TMath::Sqrt(mean) + mean +0.5);
      return n;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Generates a random number according to a Poisson law.
/// Prob(N) = exp(-mean)*mean^N/Factorial(N)
///
/// This function is a variant of TRandom::Poisson returning a double
/// instead of an integer.

Double_t TRandom::PoissonD(Double_t mean)
{
   Int_t n;
   if (mean <= 0) return 0;
   if (mean < 25) {
      Double_t expmean = TMath::Exp(-mean);
      Double_t pir = 1;
      n = -1;
      while(1) {
         n++;
         pir *= Rndm();
         if (pir <= expmean) break;
      }
      return static_cast<Double_t>(n);
   }
   // for large value we use inversion method
   else if (mean < 1E9) {
      Double_t em, t, y;
      Double_t sq, alxm, g;
      Double_t pi = TMath::Pi();

      sq = TMath::Sqrt(2.0*mean);
      alxm = TMath::Log(mean);
      g = mean*alxm - TMath::LnGamma(mean + 1.0);

      do {
         do {
            y = TMath::Tan(pi*Rndm());
            em = sq*y + mean;
         } while( em < 0.0 );

         em = TMath::Floor(em);
         t = 0.9*(1.0 + y*y)* TMath::Exp(em*alxm - TMath::LnGamma(em + 1.0) - g);
      } while( Rndm() > t );

      return em;

   } else {
      // use Gaussian approximation vor very large values
      return Gaus(0,1)*TMath::Sqrt(mean) + mean +0.5;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return 2 numbers distributed following a gaussian with mean=0 and sigma=1.

void TRandom::Rannor(Float_t &a, Float_t &b)
{
   Double_t r, x, y, z;

   y = Rndm();
   z = Rndm();
   x = z * 6.28318530717958623;
   r = TMath::Sqrt(-2*TMath::Log(y));
   a = (Float_t)(r * TMath::Sin(x));
   b = (Float_t)(r * TMath::Cos(x));
}

////////////////////////////////////////////////////////////////////////////////
/// Return 2 numbers distributed following a gaussian with mean=0 and sigma=1.

void TRandom::Rannor(Double_t &a, Double_t &b)
{
   Double_t r, x, y, z;

   y = Rndm();
   z = Rndm();
   x = z * 6.28318530717958623;
   r = TMath::Sqrt(-2*TMath::Log(y));
   a = r * TMath::Sin(x);
   b = r * TMath::Cos(x);
}

////////////////////////////////////////////////////////////////////////////////
/// Reads saved random generator status from filename.

void TRandom::ReadRandom(const char *filename)
{
   if (!gDirectory) return;
   char *fntmp = gSystem->ExpandPathName(filename);
   TDirectory *file = (TDirectory*)gROOT->ProcessLine(Form("TFile::Open(\"%s\");",fntmp));
   delete [] fntmp;
   if(file && file->GetFile()) {
      gDirectory->ReadTObject(this,GetName());
      delete file;
   }
}

////////////////////////////////////////////////////////////////////////////////
///  Machine independent random number generator.
///  Based on the BSD Unix (Rand) Linear congrential generator.
///  Produces uniformly-distributed floating points between 0 and 1.
///  Identical sequence on all machines of >= 32 bits.
///  Periodicity = 2**31, generates a number in (0,1).
///  Note that this is a generator which is known to have defects
///  (the lower random bits are correlated) and therefore should NOT be
///  used in any statistical study).

Double_t TRandom::Rndm( )
{
#ifdef OLD_TRANDOM_IMPL
   const Double_t kCONS = 4.6566128730774E-10;
   const Int_t kMASK24  = 2147483392;

   fSeed *= 69069;
   UInt_t jy = (fSeed&kMASK24); // Set lower 8 bits to zero to assure exact float
   if (jy) return kCONS*jy;
   return Rndm();
#endif

   // kCONS = 1./2147483648 = 1./(RAND_MAX+1) and RAND_MAX= 0x7fffffffUL
   const Double_t kCONS = 4.6566128730774E-10; // (1/pow(2,31)
   fSeed = (1103515245 * fSeed + 12345) & 0x7fffffffUL;

   if (fSeed) return  kCONS*fSeed;
   return Rndm();
}

////////////////////////////////////////////////////////////////////////////////
/// Return an array of n random numbers uniformly distributed in ]0,1].

void TRandom::RndmArray(Int_t n, Double_t *array)
{
   const Double_t kCONS = 4.6566128730774E-10; // (1/pow(2,31))
   Int_t i=0;
   while (i<n) {
      fSeed = (1103515245 * fSeed + 12345) & 0x7fffffffUL;
      if (fSeed) {array[i] = kCONS*fSeed; i++;}
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return an array of n random numbers uniformly distributed in ]0,1].

void TRandom::RndmArray(Int_t n, Float_t *array)
{
   const Double_t kCONS = 4.6566128730774E-10; // (1/pow(2,31))
   Int_t i=0;
   while (i<n) {
      fSeed = (1103515245 * fSeed + 12345) & 0x7fffffffUL;
      if (fSeed) {array[i] = Float_t(kCONS*fSeed); i++;}
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the random generator seed. Note that default value is zero, which is
/// different than the default value used when constructing the class.
/// If the seed is zero the seed is set to a random value
/// which in case of TRandom depends on the lowest 4 bytes of TUUID
/// The UUID will be identical if SetSeed(0) is called with time smaller than 100 ns
/// Instead if a different generator implementation is used (TRandom1, 2 or 3)
/// the seed is generated using a 128 bit UUID. This results in different seeds
/// and then random sequence for every SetSeed(0) call.

void TRandom::SetSeed(ULong_t seed)
{
   if( seed==0 ) {
      TUUID u;
      UChar_t uuid[16];
      u.GetUUID(uuid);
      fSeed  =  UInt_t(uuid[3])*16777216 + UInt_t(uuid[2])*65536 + UInt_t(uuid[1])*256 + UInt_t(uuid[0]);
   } else {
      fSeed = seed;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the random generator seed.
/// Note that this function returns the given seed only when using
/// as random generator engine TRandom itself, which is an LCG generator
/// and it has as seed (state) only one 32 bit word.
/// In case of the other generators GetSeed will return one of the state elements and not the
/// given seed. See the documentation of the corresponding generator used
/// (for example TRandom3::GetSeed() when using TRandom3 or gRandom.
/// If one needs to save the generator seed in order to be used later for obtaining reproducible
/// numbers, one should store the full generator, either in a file or in memory in a separate TRandom
/// object. Here is an example on how to store reproducible states:
/// ```
/// // set a unique seed
///  gRandom->SetSeed(0);
///  // save generator state in a different TRandom instance
///  TRandom* rngSaved = static_cast<TRandom*>(gRandom->Clone());
///  // now both rngSaved and gRandom will produce the same sequence of numbers
///  for (int i = 0; i < 10; ++i )
///     std::cout << "genrated number from gRandom : " << gRandom->Rndm() << "  from saved generator " <<
///     rngSaved->Rndm() << std::endl;
/// ```
UInt_t TRandom::GetSeed() const
{
   return fSeed;
}

////////////////////////////////////////////////////////////////////////////////
/// Generates random vectors, uniformly distributed over the surface
/// of a sphere of given radius.
///   Input : r = sphere radius
///   Output: x,y,z a random 3-d vector of length r
/// Method: (based on algorithm suggested by Knuth and attributed to Robert E Knop)
///         which uses less random numbers than the CERNLIB RN23DIM algorithm

void TRandom::Sphere(Double_t &x, Double_t &y, Double_t &z, Double_t r)
{
   Double_t a=0,b=0,r2=1;
   while (r2 > 0.25) {
      a  = Rndm() - 0.5;
      b  = Rndm() - 0.5;
      r2 =  a*a + b*b;
   }
   z = r* ( -1. + 8.0 * r2 );

   Double_t scale = 8.0 * r * TMath::Sqrt(0.25 - r2);
   x = a*scale;
   y = b*scale;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a uniform deviate on the interval  (0, x1).

Double_t TRandom::Uniform(Double_t x1)
{
   Double_t ans = Rndm();
   return x1*ans;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a uniform deviate on the interval (x1, x2).

Double_t TRandom::Uniform(Double_t x1, Double_t x2)
{
   Double_t ans= Rndm();
   return x1 + (x2-x1)*ans;
}

////////////////////////////////////////////////////////////////////////////////
/// Writes random generator status to filename.

void TRandom::WriteRandom(const char *filename) const
{
   if (!gDirectory) return;
   char *fntmp = gSystem->ExpandPathName(filename);
   TDirectory *file = (TDirectory*)gROOT->ProcessLine(Form("TFile::Open(\"%s\",\"recreate\");",fntmp));
   delete [] fntmp;
   if(file && file->GetFile()) {
      gDirectory->WriteTObject(this,GetName());
      delete file;
   }
}
