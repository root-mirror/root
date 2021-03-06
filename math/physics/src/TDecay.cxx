/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TDecay.h"

// const Int_t kMAXP = 18;

//_____________________________________________________________________________________
Double_t TDecay::PDK(Double_t a, Double_t b, Double_t c)
{
   // the PDK function
   Double_t x = (a - b - c) * (a + b + c) * (a - b + c) * (a + b - c);
   x = TMath::Sqrt(x) / (2 * a);
   return x;
}

namespace {
//_____________________________________________________________________________________
Int_t DoubleMax(const void *a, const void *b)
{
   // special max function
   Double_t aa = *((Double_t *)a);
   Double_t bb = *((Double_t *)b);
   if (aa > bb)
      return 1;
   if (aa < bb)
      return -1;
   return 0;
}
} // namespace

//__________________________________________________________________________________________________
TDecay::TDecay(const TDecay &gen) : TObject(gen)
{
   // copy constructor
   fNt = gen.fNt;
   fTeCmTm = gen.fTeCmTm;
   fBeta[0] = gen.fBeta[0];
   fBeta[1] = gen.fBeta[1];
   fBeta[2] = gen.fBeta[2];
   for (Int_t i = 0; i < fNt; i++) {
      fMass[i] = gen.fMass[i];
      fDecPro[i] = gen.fDecPro[i];
   }
}

//__________________________________________________________________________________________________
TDecay &TDecay::operator=(const TDecay &gen)
{
   // Assignment operator
   TObject::operator=(gen);
   fNt = gen.fNt;
   fTeCmTm = gen.fTeCmTm;
   fBeta[0] = gen.fBeta[0];
   fBeta[1] = gen.fBeta[1];
   fBeta[2] = gen.fBeta[2];
   for (Int_t i = 0; i < fNt; i++) {
      fMass[i] = gen.fMass[i];
      fDecPro[i] = gen.fDecPro[i];
   }

   return *this;
}

//__________________________________________________________________________________________________
int TDecay::factorial(int n)
{
   // the implementation is sufficient for small n
   return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

//__________________________________________________________________________________________________
Double_t TDecay::Generate(std::queue<double> &rnd)
{

   Double_t rno[kMAXP];
   rno[0] = 0;
   Int_t n;

   // check if phase space dim is ok
   assert((unsigned int)(3 * fNt - 4) == rnd.size());

   double weight = 1.0; // additional weight

   if (fNt > 2) {
      for (n = 1; n < fNt - 1; n++) {
         rno[n] = rnd.front(); // fNt-2 random numbers
         rnd.pop();
      }

      qsort(rno + 1, fNt - 2, sizeof(Double_t), DoubleMax); // sort them

      weight /= double(factorial(fNt - 2));
   }
   rno[fNt - 1] = 1;

   Double_t invMas[kMAXP], sum = 0;
   for (n = 0; n < fNt; n++) {
      sum += fMass[n];
      invMas[n] = rno[n] * fTeCmTm + sum;

      if ((n == 0) || (n == fNt - 1))
         continue;
      weight *= 2.0 * invMas[n] * fTeCmTm;
   }

   //
   //-----> compute the weight of the current event
   //
   Double_t wt = 1.0;
   Double_t pd[kMAXP];
   for (n = 0; n < fNt - 1; n++) {
      pd[n] = PDK(invMas[n + 1], invMas[n], fMass[n + 1]);
      wt *= pd[n];
      wt *= M_PI;
      wt /= invMas[n + 1] * pow(2.0 * M_PI, 3);
   }
   weight *= 2.0 * M_PI;
   wt *= weight;

   // New boost method - adjusted to spherical coordinates
   fDecPro[0].SetPxPyPzE(0, 0, pd[0], TMath::Sqrt(pd[0] * pd[0] + fMass[0] * fMass[0]));

   Int_t i = 1;
   Int_t j;
   while (1) {
      fDecPro[i].SetPxPyPzE(0, 0, -pd[i - 1], TMath::Sqrt(pd[i - 1] * pd[i - 1] + fMass[i] * fMass[i]));

      Double_t cY = 2 * rnd.front() - 1;
      rnd.pop();
      Double_t sY = TMath::Sqrt(1 - cY * cY);
      Double_t angZ = 2 * TMath::Pi() * rnd.front();
      rnd.pop();
      Double_t cZ = TMath::Cos(angZ);
      Double_t sZ = TMath::Sin(angZ);
      for (j = 0; j <= i; j++) {
         TLorentzVector *v = fDecPro + j;
         Double_t x = v->Px();
         Double_t z = v->Pz();
         v->SetPz(cY * z - sY * x);
         v->SetPx(sY * z + cY * x); // rotation around Y
         x = v->Px();
         Double_t y = v->Py();
         v->SetPx(cZ * x - sZ * y);
         v->SetPy(sZ * x + cZ * y); // rotation around Z
      }

      if (i == (fNt - 1))
         break;

      Double_t beta = pd[i] / sqrt(pd[i] * pd[i] + invMas[i] * invMas[i]);
      for (j = 0; j <= i; j++)
         fDecPro[j].Boost(0, 0, beta);
      i++;
   }

   //
   //---> final boost of all particles
   //
   for (n = 0; n < fNt; n++)
      fDecPro[n].Boost(fBeta[0], fBeta[1], fBeta[2]);

   //
   //---> return the weigth of event
   //
   return wt;
}

//__________________________________________________________________________________
TLorentzVector *TDecay::GetDecay(Int_t n)
{
   // return Lorentz vector corresponding to decay n
   if (n > fNt)
      return 0;
   return fDecPro + n;
}

//_____________________________________________________________________________________
Bool_t TDecay::SetDecay(TLorentzVector &P, Int_t nt, const Double_t *mass)
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

   //---->  save the betas of the decaying particle
   //
   if (P.Beta()) {
      Double_t w = P.Beta() / P.Rho();
      fBeta[0] = P(0) * w;
      fBeta[1] = P(1) * w;
      fBeta[2] = P(2) * w;
   } else
      fBeta[0] = fBeta[1] = fBeta[2] = 0;

   return kTRUE;
}
