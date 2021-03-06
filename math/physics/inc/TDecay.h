/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/*!
  @class TDecay

  @brief  The generator decays a central blob of a given 4-momentum into particles of specific masses. Original GENBOD
algorithm with w few modifications.

### Authors

Radosław Kycia (kycia.radoslaw@gmail.com), Piotr Lebiedowicz, Antoni Szczurek

People who helped in the development of the project: Jacek Turnau, Janusz Chwastowski, Rafał Staszewski, Maciej
Trzebiński.

### Description

This is a helper class that can be base to construct MC generators that simulate decay. It returns weight equals LIPS
(Lorentz Invariant Phase Space):

\f$ dLIPS = (2\pi)^4 \delta^{(4)}(P-\sum_{i=1}^{n}p_{i}) \prod_{i=1}^{n} \frac{d^{3}p_{i}}{(2\pi)^3 2E_{i}}\f$

The generator can be used as a building block in more complicated decays.

The generator returns weighted events.

The class is adapted from the TGenPhaseSpace ROOT package.

This class is not recommended for non-advanced users. If you want a robust and adaptive MC generator and integrator, use
TGenFoamDecay instead.

#### The scheme of use

1. Prepare 4-momenta of decaying particle P and mass[nt] array of final products.
2. Initialize generator:  SetDecay();
3. Generate decay: Generate(); The class requires random numbers given to Generate(). They are used to make an event.
4. For each of particles enumerated by 0..nt-1 get their 4-momentum using GetDecay(): pfi = GetDecay(i);
5. Repeat 3 and 4 for another decay.


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

*/

#ifndef ROOT_TDecay
#define ROOT_TDecay

#include <TLorentzVector.h>

#include <iostream> // std::cin, std::cout
#include <queue>    // std::queue
#include <assert.h>

using namespace std;

class TDecay : public TObject {
private:
   Int_t fNt;                  // number of decay particles
   Double_t fMass[18];         // masses of particles
   Double_t fBeta[3];          // betas of decaying particle
   Double_t fTeCmTm;           // total energy in the C.M. minus the total mass
   TLorentzVector fDecPro[18]; // kinematics of the generated particles

   Int_t kMAXP; // max number of particles (relict from TGenPhaseSpace)

   Double_t PDK(Double_t a, Double_t b, Double_t c);

   /// factorial function
   int factorial(int n);

public:
   /// constructor
   TDecay() : fNt(0), fMass(), fBeta(), fTeCmTm(0.){};

   /// copy constructor
   TDecay(const TDecay &gen);

   /// desctructor
   virtual ~TDecay(){};

   /// assignment
   TDecay &operator=(const TDecay &gen);

   /// Sets up configuration of decay
   /// @param P 4-momentum of decaying particle
   /// @param nt number of products of decay
   /// @param mass mass matrix of products of decay mass[nt]
   /// @returns kTRUE - the decay is permitted by kinematics; kFALSE - the decay is forbidden by kinematics
   /// @warning This should be first method to call since it sets up decay configuration.
   Bool_t SetDecay(TLorentzVector &P, Int_t nt, const Double_t *mass);

   /// Generate a single decay
   /// @param rnd a queue of 3*nt-4 random numbers from external source.
   Double_t Generate(std::queue<double> &rnd);

   /// @returns 4-vector of n-th product of decay
   /// @param n  number of final particle to get from range 1...nt
   /// @warning You should call Generate() in first place.
   TLorentzVector *GetDecay(Int_t n);

   /// @returns number of final particles
   Int_t GetNt() const { return fNt; }
};

#endif
