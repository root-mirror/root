/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGenFoamDecay.h"

// const Int_t kMAXP = 18;

Double_t TGenFoamDecay::Density(int /*nDim*/, Double_t *Xarg)
{

   // queue for random numbers
   queue<double> rndQueue;

   // put rnd numbers into queue
   for (int i = 0; i < 3 * fNt - 4; i++) {
      rndQueue.push(Xarg[i]);
   }

   // make decay and take d(LIPS)
   double wtdecay = fDecay.Generate(rndQueue);

   // get out particles
   TLorentzVector pf[fNt];
   for (int i = 0; i < fNt; i++) {
      pf[i] = *(fDecay.GetDecay(i));
   }

   // calculate integrand
   double integrand = Integrand(fNt, pf);

   return wtdecay * integrand;
}

//__________________________________________________________________________________________________
Double_t TGenFoamDecay::Integrand(int /*nt*/, TLorentzVector * /*pf*/)
{
   return 1.0; // default and probably overloaded for matrix element
}

//__________________________________________________________________________________________________
TGenFoamDecay::TGenFoamDecay(const TGenFoamDecay &gen)
{
   // copy constructor
   fNt = gen.fNt;
   fTeCmTm = gen.fTeCmTm;
   fBeta[0] = gen.fBeta[0];
   fBeta[1] = gen.fBeta[1];
   fBeta[2] = gen.fBeta[2];
   fDecay = gen.fDecay;
   fFoam = gen.fFoam;
   fPseRan = gen.fPseRan;
   for (Int_t i = 0; i < fNt; i++) {
      fMass[i] = gen.fMass[i];
      fDecPro[i] = gen.fDecPro[i];
   }
}

//__________________________________________________________________________________________________
TGenFoamDecay &TGenFoamDecay::operator=(const TGenFoamDecay &gen)
{
   // Assignment operator
   TObject::operator=(gen);
   fNt = gen.fNt;
   fTeCmTm = gen.fTeCmTm;
   fBeta[0] = gen.fBeta[0];
   fBeta[1] = gen.fBeta[1];
   fBeta[2] = gen.fBeta[2];
   fDecay = gen.fDecay;
   fFoam = gen.fFoam;
   fPseRan = gen.fPseRan;
   for (Int_t i = 0; i < fNt; i++) {
      fMass[i] = gen.fMass[i];
      fDecPro[i] = gen.fDecPro[i];
   }
   return *this;
}

//__________________________________________________________________________________________________
Double_t TGenFoamDecay::Generate(void)
{

   fFoam->MakeEvent();

   return fFoam->GetMCwt();
}

//__________________________________________________________________________________
TLorentzVector *TGenFoamDecay::GetDecay(Int_t n)
{

   if (n > fNt)
      return 0;

   // return Lorentz vector corresponding to decay of n-th particle
   return fDecay.GetDecay(n);
}

//_____________________________________________________________________________________
Bool_t TGenFoamDecay::SetDecay(TLorentzVector &P, Int_t nt, const Double_t *mass)
{

   kMAXP = nt;

   Int_t n;
   fNt = nt;
   if (fNt < 2 || fNt > 18)
      return kFALSE; // no more then 18 particle

   fTeCmTm = P.Mag(); // total energy in C.M. minus the sum of the masses
   for (n = 0; n < fNt; n++) {
      fMass[n] = mass[n];
      fTeCmTm -= mass[n];
   }

   if (fTeCmTm <= 0)
      return kFALSE; // not enough energy for this decay

   fDecay.SetDecay(P, fNt, fMass); // set decay to TDecay

   // initialize FOAM
   //=========================================================
   if (fChat > 0) {
      cout << "*****   Foam version " << fFoam->GetVersion() << "    *****" << endl;
   }
   fFoam->SetkDim(3 * fNt - 4);  // Mandatory!!!
   fFoam->SetnCells(fNCells);     // optional
   fFoam->SetnSampl(fNSampl);     // optional
   fFoam->SetnBin(fNBin);         // optional
   fFoam->SetOptRej(fOptRej);     // optional
   fFoam->SetOptDrive(fOptDrive); // optional
   fFoam->SetEvPerBin(fEvPerBin); // optional
   fFoam->SetChat(fChat);         // optional
   //===============================
   fFoam->SetRho(this);
   fFoam->SetPseRan(&fPseRan);

   // Initialize simulator
   fFoam->Initialize();

   return kTRUE;
}

//_____________________________________________________________________________________
void TGenFoamDecay::Finalize(void)
{
   Double_t MCresult, MCerror;
   Double_t eps = 0.0005;
   Double_t Effic, WtMax, AveWt, Sigma;
   Double_t IntNorm, Errel;
   fFoam->Finalize(IntNorm, Errel);              // final printout
   fFoam->GetIntegMC(MCresult, MCerror);         // get MC intnegral
   fFoam->GetWtParams(eps, AveWt, WtMax, Sigma); // get MC wt parameters
   long nCalls = fFoam->GetnCalls();
   Effic = 0;
   if (WtMax > 0)
      Effic = AveWt / WtMax;
   cout << "================================================================" << endl;
   cout << " MCresult= " << MCresult << " +- " << MCerror << " RelErr= " << MCerror / MCresult << endl;
   cout << " Dispersion/<wt>= " << Sigma / AveWt << endl;
   cout << "      <wt>/WtMax= " << Effic << ",    for epsilon = " << eps << endl;
   cout << " nCalls (initialization only) =   " << nCalls << endl;
   cout << "================================================================" << endl;
}

//_____________________________________________________________________________________
void TGenFoamDecay::GetIntegMC(Double_t &integral, Double_t &error)
{
   fFoam->GetIntegMC(integral, error);
}
