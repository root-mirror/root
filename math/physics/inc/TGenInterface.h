/*************************************************************************
* Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/ 

/*!
  @class TGenInterface

  @brief  Minimal interface implemented by MC generators for the decay of a central blob of a given 4-momentum into particles of specific masses. 
  
### Authors
	
Radosław Kycia (kycia.radoslaw@gmail.com), Piotr Lebiedowicz, Antoni Szczurek

People who helped in development of the project: Jacek Turnau, Janusz Chwastowski, Rafał Staszewski, Maciej Trzebiński.

### Details
The general concept of working with MC generators is as follows:

1. Prepare 4-momenta of decaying particle P and mass[nt] array of final products.
2. Initialize generator:  SetDecay();
3. Generate decay: Generate();
4. For each of particles enumerated by 0..nt-1 get their 4-momentum using GetDecay(): pfi = GetDecay(i);
5. Repeat 3 and 4 for another decay.
 
@warning The class is an interface for integration of some distribution over phase space of decaying particle. In the end, you would probably like to normalize distribution properly.
 
The interface is implemented in:
- TGenFoamDecay - adaptive Monte Carlo generator and integrator; Probably contains all you need.
- TGenDecay - non-adaptive Monte Carlo generator; returns proper LIPS (Lorentz Invariant Phase Space) weight for an event.
- TDecay - the base class that returns LISP for an event; Can be used in more advanced MC generators as a part that makes decay.
- TGenPhaseSpace - legacy MC generator.

### Further details

- R.A. Kycia, P. Lebiedowicz, A. Szczurek, 'Decay: A Monte Carlo library for the decay of a particle with ROOT compatibility', arXiv:2011.14750 [hep-ph]
- R. A. Kycia, J. Chwastowski, R. Staszewski, and J. Turnau, 'GenEx: A simple generator structure for exclusive processes in high energy collisions', Commun. Comput. Phys.24no. 3, (2018) 860, arXiv:1411.6035 [hep-ph]
- R. A. Kycia, J. Turnau, J. J. Chwastowski, R. Staszewski, and M. Trzebinski, 'The adaptiveMonte Carlo toolbox for phase space integration and generation', Commun. Comput. Phys.25no. 5, (2019) 1547, arXiv:1711.06087 [hep-ph]
- R. Hagedorn, 'Relativistic Kinematics: A guide to the kinematic problems of high-energy physics'. W.A. Benjamin, New York, Amsterdam, 1964
- H.Pilkhun, The Interactions of Hadrons North-Holland 1967 
 
*/



#ifndef ROOT_TGenInterface
#define ROOT_TGenInterface

#include <TLorentzVector.h>


using namespace std;


class TGenInterface {
private:  

   Int_t  fNt;  // number of decay particles
   
public:

	///constructor
	TGenInterface(): fNt(0) {};
	
	///destructor
	virtual ~TGenInterface() {};

	/// Sets up configuration of decay
	/// @param P 4-momentum of decaying particle
	/// @param nt number of products of decay
	/// @param mass mass matrix of products of decay mass[nt]
	/// @returns kTRUE - the decay is permitted by kinematics; kFALSE - the decay is forbidden by kinematics
	/// @warning This should be first method to call since it sets up decay configuration.
	virtual Bool_t SetDecay(TLorentzVector &P, Int_t nt, const Double_t *mass)=0;
	
	/// Generate a single decay
	virtual Double_t  Generate(void)=0;
	
	/// @returns 4-vector of n-th product of decay
	/// @param n  number of final particle to get from range 1...nt
	/// @warning You should call Generate() in first place.
	virtual TLorentzVector *GetDecay(Int_t n)=0; 

	/// @returns number of final particles
	virtual Int_t GetNt()const { return fNt;};
   
	
};

#endif

