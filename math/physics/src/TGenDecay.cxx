/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGenDecay.h"

// const Int_t kMAXP = 18;

//__________________________________________________________________________________________________
TGenDecay::TGenDecay(const TGenDecay &gen) : TObject(gen)
{
   // copy constructor
   fNt = gen.fNt;
   fTeCmTm = gen.fTeCmTm;
   fBeta[0] = gen.fBeta[0];
   fBeta[1] = gen.fBeta[1];
   fBeta[2] = gen.fBeta[2];
   fDecay = gen.fDecay;
   fPseRan = gen.fPseRan;
   for (Int_t i = 0; i < fNt; i++) {
      fMass[i] = gen.fMass[i];
      fDecPro[i] = gen.fDecPro[i];
   }
}

//__________________________________________________________________________________________________
TGenDecay &TGenDecay::operator=(const TGenDecay &gen)
{
   // Assignment operator
   TObject::operator=(gen);
   fNt = gen.fNt;
   fTeCmTm = gen.fTeCmTm;
   fBeta[0] = gen.fBeta[0];
   fBeta[1] = gen.fBeta[1];
   fBeta[2] = gen.fBeta[2];
   fDecay = gen.fDecay;
   fPseRan = gen.fPseRan;
   for (Int_t i = 0; i < fNt; i++) {
      fMass[i] = gen.fMass[i];
      fDecPro[i] = gen.fDecPro[i];
   }
   return *this;
}

//__________________________________________________________________________________________________
Double_t TGenDecay::Generate(void)
{

   // clear queue - in case there was a previous run
   while (!fRndQueue.empty()) {
      fRndQueue.pop();
   }

   // number of degrees of freedom in the decay
   int nDim = 3 * fNt - 4;

   // put rnd numbers into queue
   for (int i = 0; i < nDim; i++) {
      fRndQueue.push(fPseRan.Rndm());
   }

   return fDecay.Generate(fRndQueue);
}

//__________________________________________________________________________________
TLorentzVector *TGenDecay::GetDecay(Int_t n)
{
   // return Lorentz vector corresponding to decay n
   if (n > fNt)
      return 0;

   return fDecay.GetDecay(n);
}

//_____________________________________________________________________________________
Bool_t TGenDecay::SetDecay(TLorentzVector &P, Int_t nt, const Double_t *mass)
{

   kMAXP = nt;

   Int_t n;
   fNt = nt;
   if (fNt < 2 || fNt > 18)
      return kFALSE; // no more then 18 particle

   //
   fTeCmTm = P.Mag(); // total energy in C.M. minus the sum of the masses
   for (n = 0; n < fNt; n++) {
      fMass[n] = mass[n];
      fTeCmTm -= mass[n];
   }

   if (fTeCmTm <= 0)
      return kFALSE; // not enough energy for this decay

   fPseRan.SetSeed(fSeed); // set seed

   fDecay.SetDecay(P, fNt, fMass); // set decay to TDecay

   return kTRUE;
}
