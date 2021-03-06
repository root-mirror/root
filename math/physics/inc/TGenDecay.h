/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/*!
  @class TGenDecay

@brief  Generates spherical decay and returns the weight of event equals to LIPS (Lorentz Invariant Phase Space). The
generator decays a central blob of a given 4-momentum into particles of specific masses.

### Authors

Radosław Kycia (kycia.radoslaw@gmail.com), Piotr Lebiedowicz, Antoni Szczurek

People who helped in the development of the project: Jacek Turnau, Janusz Chwastowski, Rafał Staszewski, Maciej
Trzebiński.


### Details

Generator implements TGenInterface interface and works with it is along with this interface.

The general way of working with the generator is as follows:
1. Prepare 4-momenta of decaying particle P and mass[nt] array of final products.
2. Initialize generator:  SetDecay();
3. Generate decay: Generate();
4. For each of particles enumerated by 0..nt-1 get their 4-momentum using GetDecay(): pfi = GetDecay(i);
5. Repeat 3 and 4 for another decay.

Generate() method returns the LIPS:

\f$ dLIPS = (2\pi)^4 \delta^{(4)}(P-\sum_{i=1}^{n}p_{i}) \prod_{i=1}^{n} \frac{d^{3}p_{i}}{(2\pi)^3 2E_{i}}\f$

Comparing to TDecay, it contains a random number generator.

The generator returns weighted events.

This class is not recommended for non-advanced users. If you want a robust and adaptive MC generator and integrator, use
TGenFoamDecay instead.


### Further details

- R.A. Kycia, P. Lebiedowicz, A. Szczurek, 'Decay: A Monte Carlo library for the decay of a particle with ROOT
compatibility', arXiv:2011.14750 [hep-ph]
- R. A. Kycia, J. Chwastowski, R. Staszewski, and J. Turnau, 'GenEx: A simple generator structure for exclusive
processes in high energy collisions', Commun. Comput. Phys.24no. 3, (2018) 860, arXiv:1411.6035 [hep-ph]
- R. A. Kycia, J. Turnau, J. J. Chwastowski, R. Staszewski, and M. Trzebinski, 'The adaptiveMonte Carlo toolbox for
phase space integration and generation', Commun. Comput. Phys.25no. 5, (2019) 1547, arXiv:1711.06087 [hep-ph]
- R. Hagedorn, 'Relativistic Kinematics: A guide to the kinematic problems of high-energy physics'. W.A. Benjamin, New
York, Amsterdam, 1964
- H.Pilkhun, The Interactions of Hadrons North-Holland 1967
People who helped in the development of the project:

*/

#ifndef ROOT_TGenDecay
#define ROOT_TGenDecay

#include <TGenInterface.h>
#include <TDecay.h>

#include <TLorentzVector.h>
#include <TRandom3.h>

#include <iostream> // std::cin, std::cout
#include <queue>    // std::queue
#include <assert.h>

using namespace std;

class TGenDecay : public TObject, public TGenInterface {
private:
   Int_t fNt;                  // number of decay particles
   Double_t fMass[18];         // masses of particles
   Double_t fBeta[3];          // betas of decaying particle
   Double_t fTeCmTm;           // total energy in the C.M. minus the total mass
   TLorentzVector fDecPro[18]; // kinematics of the generated particles
   Int_t fSeed;                 // seed for pseudorandom number generator

   Int_t kMAXP; // max number of particles (relict from TGenPhaseSpace)

   TDecay fDecay;    // decay engine
   TRandom3 fPseRan; // pseudorandom numbers

   queue<double> fRndQueue; // queue for random numbers

public:
   /// constructor
   TGenDecay() : fNt(0), fMass(), fBeta(), fTeCmTm(0.), fSeed(4357) {}
   /// copy constructor
   TGenDecay(const TGenDecay &gen);
   /// desctuctor
   virtual ~TGenDecay() {}
   /// assignment
   TGenDecay &operator=(const TGenDecay &gen);

   /// Sets up configuration of decay
   /// @param P 4-momentum of decaying particle (Momentum, Energy units are Gev/C, GeV)
   /// @param nt number of products of decay
   /// @param mass mass matrix of products of decay mass[nt]
   /// @returns kTRUE - the decay is permitted by kinematics; kFALSE - the decay is forbidden by kinematics
   /// @warning This should be first method to call since it sets up decay configuration.
   /// @warning The method also initialize FOAM.
   Bool_t SetDecay(TLorentzVector &P, Int_t nt, const Double_t *mass);

   /// Generate a single decay
   Double_t Generate(void);

   /// Collect 4-vector of products of decay
   ///  @param n  number of final particle to get from range 1...nt
   /// @warning You should call Generate() in first place.
   TLorentzVector *GetDecay(Int_t n);

   /// @returns 4-vector of n-th product of decay
   ///  @param n  number of final particle to get from range 1...nt
   Int_t GetNt() const { return fNt; }

   /// sets seed for pseudorandom number generator
   void setSeed(Int_t seed)
   {
      fSeed = seed;
      fPseRan.SetSeed(seed);
   };
};

#endif
