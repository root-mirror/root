/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/*!
  @class TGenFoamDecay


  @brief  Adaptive MC generator of decay and integrator over LIPS. The generator decays a central blob of a given
4-momentum into particles of specific masses.

### Authors

Radosław Kycia (kycia.radoslaw@gmail.com), Piotr Lebiedowicz, Antoni Szczurek

People who helped in development of the project: Jacek Turnau, Janusz Chwastowski, Rafał Staszewski, Maciej Trzebiński.

### Details

This class is a multipurpose MC generator for decay, which also integrates functions over phase space. It uses FOAM
adaptive integrator.

To use this class, make your own generator class that inherits from TGenFoamDecay and redefine Integrand() method,
filling it with an integrand function that will be integrated over LIPS (Lorent Invariant Phase Space)  and possibly
kinematics cuts (return 0.0). If you do not have any specific integrand, then the best choice is to make Integrand()
return 1.0.

The volume element of LIPS is:

\f$ dLIPS = (2\pi)^4 \delta^{(4)}(P-\sum_{i=1}^{n}p_{i}) \prod_{i=1}^{n} \frac{d^{3}p_{i}}{(2\pi)^3 2E_{i}}\f$

Generator implements interface TGenInterface and work with it is along with this interface.

The general way of working with the generator is as follows:

1. Create your generator class derived from TGenFoamDecay and implement Integrand() method - return 1.0 if you do not
have any specific integrand.
2. Prepare 4-momenta of decaying particle P and mass[nt] array of final products.
3. Initialize generator:  SetDecay();
4. Generate decay: Generate();
5. For each of particles enumerated by 0..nt-1 get their 4-momentum using GetDecay(): pfi = GetDecay(i);
6. Repeat 3 and 4 for another decay.

If you prepared the integrand to be integrated over LIPS, you can print the integral final value by calling the
Finalize() method or getting the numerical value by GetIntegMC().

### Troubleshooting and suggestions

- By default, the weight of the event is 1.0; however, see OptRej parameter. Sometimes FOAM, even for OptRej=1 weight
are not precisely 1.0, so you can remove them or use weights, e.g., in histograms.
- If you do not have an integrand function (e.g., matrix element), then the optimal choice for Integrand() method is to
return 1.0 (uniform distribution in the phase space).
- If the integrand has 'small' support in phase space and you get some nonsense result/errors, you probably should
increase (say 10x or more) nSampl and nCells parameters of FOAM. It will increase the probability that adaptive MC will
have a chance to spot support of the function. However, the exploration phase will last longer. If FOAM still cannot
spot the support of the integral, then you can think of using a special Monte Carlo generator for non-spherical decay.
- FOAM has many various setting that can be tested when you deal with demanding integrand functions, see TFoam
documentation or the Foam paper:
 * S. Jadach, Foam: A General-Purpose Cellular Monte Carlo Event Generator, Comput.Phys.Commun. 152 (2003) 55-100
- Kinematic cuts can be effectively implemented in Integrand() method - just simply return 0.0 when kinematic cuts are
met.
- Foam produces a lot of diagnostic output. It is good to read it at the beginning of implementation. However, in
production code or when the generator is embedded in more extensive software, it is better to turn them off. To suppress
them set Chat to 0.


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
- S. Jadach, Foam: A General-Purpose Cellular Monte Carlo Event Generator, Comput.Phys.Commun. 152 (2003) 55-100

*/

#ifndef ROOT_TGenFoamDecay
#define ROOT_TGenFoamDecay

#include <TGenInterface.h>

#include <TLorentzVector.h>
#include <TFoam.h>
#include <TFoamIntegrand.h>
#include <TRandom3.h>

#include <stdlib.h>
#include <iostream> // std::cin, std::cout
#include <queue>    // std::queue
#include <assert.h>

#include <TDecay.h>

using namespace std;

class TGenFoamDecay : public TFoamIntegrand, public TGenInterface {
private:
   Int_t fNt;                  // number of decay particles
   Double_t fMass[18];         // masses of particles
   Double_t fBeta[3];          // betas of decaying particle
   Double_t fTeCmTm;           // total energy in the C.M. minus the total mass
   TLorentzVector fDecPro[18]; // kinematics of the generated particles
   Int_t fSeed;                // seed for pseudorandom generator

   Int_t kMAXP; // max number of particles (relict from TGenPhaseSpace)

   TDecay fDecay;    // decay engine
   TFoam *fFoam;     // adaptive integrator
   TRandom3 fPseRan; // pseudorandom number generator

   // FOAM parameters
   Int_t fNCells;   // Number of Cells
   Int_t fNSampl;   // Number of MC events per cell in build-up
   Int_t fNBin;     // Number of bins in build-up
   Int_t fOptRej;   // Wted events for OptRej=0; wt=1 for OptRej=1 (default)
   Int_t fOptDrive; // (D=2) Option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
   Int_t fEvPerBin; // Maximum events (equiv.) per bin in buid-up
   Int_t fChat;     // Chat level

   /// @returns weight of the process
   /// @param nDim  dimension of integration
   /// @param Xarg vector of probablilistic variables from [0;1] of dim nDim
   /// @warning it is required by Foam integrator
   Double_t Density(int nDim, Double_t *Xarg);

public:
   /// constructor
   TGenFoamDecay()
      : fNt(0), fMass(), fBeta(), fTeCmTm(0.), fSeed(4357), fNCells(1000), fNSampl(1000), fNBin(8), fOptRej(1), fOptDrive(2),
        fEvPerBin(25), fChat(1)
   {
      fFoam = new TFoam("FoamX");
   };
   /// copy constructor
   TGenFoamDecay(const TGenFoamDecay &gen);
   /// desctructor
   virtual ~TGenFoamDecay() { delete fFoam; }
   /// assignemt
   TGenFoamDecay &operator=(const TGenFoamDecay &gen);

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
   Int_t GetNt() const { return fNt; };

   /// @returns the function under integral over LIPS (Lorentz Invariant Phase Space)
   /// @param fNt  number of outgoing particles
   /// @param pf  array of TLorentzVectors of outgoing particles
   /// @warning It is set to 1.0. You should to redefine it (in your derived class, after inheritance) when you want to
   /// use full adaptation features.
   virtual Double_t Integrand(int fNt, TLorentzVector *pf);

   /// sets seed for pseudorandom number generator
   virtual void setSeed(Int_t seed)
   {
      fSeed = seed;
      fPseRan.SetSeed(seed);
   };

   /// finalize Foam printing out MC integral and statistics
   virtual void Finalize(void);

   ///@returns integral +- error of MC integration of function in Integrand() over LIPS
   virtual void GetIntegMC(Double_t &inetgral, Double_t &error);

   /// Sets FOAM number of Cells
   /// @warning It should be done before Generate()
   void SetnCells(Int_t nc) { fNCells = nc; };

   /// Sets FOAM number of MC events per cell in build-up
   /// @warning It should be done before Generate()
   void SetnSampl(Int_t ns) { fNSampl = ns; };

   /// Sets FOAM number of bins in build-up
   /// @warning It should be done before Generate()
   void SetnBin(Int_t nb) { fNBin = nb; };

   /// Sets FOAM Weigh events for OptRej=0; wt=1 for OptRej=1 (default)
   /// @warning It determines if the events will be weighted of not (weight=1).
   /// @warning It should be done before Generate()
   void SetOptRej(Int_t OptR) { fOptRej = OptR; };

   /// Sets FOAM (D=2) option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
   /// @warning It should be done before Generate()
   void SetOptDrive(Int_t OptD) { fOptDrive = OptD; };

   /// Sets FOAM maximum events (equiv.) per bin in buid-up
   /// @warning It should be done before Generate()
   void SetEvPerBin(Int_t Ev) { fEvPerBin = Ev; };

   /// Sets FOAM chat level
   /// @warning It should be done before Generate()
   void SetChat(Int_t Ch) { fChat = Ch; };
};

#endif
