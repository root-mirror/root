/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
/*
Unit tests for TGenFoamDecay
Author: Radosław Kycia (kycia.radoslaw@gmail.com)
*/

#include <iostream>
#include <fstream>

#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>

#include "TGenFoamDecay.h"

#include <gtest/gtest.h>

using namespace std;

// decay into 2 particles in final state
TEST(PhaseSpace2Dim, Integrand1)
{

   // Properites of particles:
   TDatabasePDG *PDGDatabese = TDatabasePDG::Instance();

   // Create your own generator class derived from TGenFoamDecay
   class Generator : public TGenFoamDecay {
      Double_t Integrand(int /*nt*/, TLorentzVector * /*pf*/)
      {

         // Here put your function to integrate over LIPS (Lorentz Invariant Phase Space)
         // You can use: pf[0]....pf[Nop-1] - 4 momenta of outgoing particles
         return 1.0;
      };
   };

   // make an object of your generator
   Generator generator;

   // silent mode
   generator.SetChat(0);

   // set up initial particle/blob that decays:

   // center of mass energy
   const double mmu = 0.1057; // muon rest mass GeV

   // CM 4-momentum - blob that initiates decay
   TLorentzVector pbCM;

   // set blob CM energy to mass of the muon with zero momentum - rest frame of particle
   pbCM.SetPxPyPzE(0.0, 0.0, 0.0, mmu);

   // set up outgoing particles configuration:

   // number of outgoing particles
   const int Nop = 2;

   // PDGID outgoing paricles (masses are taken from PDG table)
   // full list: https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
   int idOut[Nop] = {11, 11}; // e, e

   // out particles array for storing their 4-vectors
   TLorentzVector pf[Nop];

   // masses of products - from PDGDatabase
   double mass[Nop];

   // fill in the table with masses of particles using their codes
   for (int i = 0; i < Nop; i++) {
      mass[i] = PDGDatabese->GetParticle(idOut[i])->Mass();
   }

   // set Generator to decay configuration
   generator.SetDecay(pbCM, Nop, mass);

   long NevTot = 10e4; // Total number of events to generate

   // GENRATION LOOP
   for (long loop = 0; loop < NevTot; loop++) {

      // Generate event

      // make decay
      generator.Generate();

      // get out particles into pf[] array - order as in mass[] array
      for (int i = 0; i < Nop; i++) {
         pf[i] = *(generator.GetDecay(i));
      }
   }
   // END of GENRATION LOOP

   // generator.Finalize();

   // get integral and error from Generator
   Double_t integral, error;
   generator.GetIntegMC(integral, error);

   // Theoretical prediction //Theoretical prediction
   double th_integral = M_PI *
                        sqrt((pow(mmu, 2) - pow(mass[0] + mass[1], 2)) * (pow(mmu, 2) - pow(mass[0] - mass[1], 2))) /
                        (2.0 * pow(mmu, 2) * pow(2.0 * M_PI, 2));

   ASSERT_NEAR(integral, th_integral, 1e-15);
};

// decay into 3 particles in final state - 4-Fermi theory of muon decay
TEST(PhaseSpace3Dim, Fermi4ModelOfDecay)
{

   // Properites of particles:
   TDatabasePDG *PDGDatabese = TDatabasePDG::Instance();

   // Create your own generator class derived from TGenFoamDecay
   class Generator : public TGenFoamDecay {
      Double_t Integrand(int /*nt*/, TLorentzVector *pf)
      {

         // Here put your function to integrate over LIPS (Lorentz Invariant Phase Space)
         // You can use: pf[0]....pf[Nop-1] - 4 momenta of outgoing particles
         const double G = 1.166e-5;
         const double mmu = 0.1057;
         double M2 = 32.0 * G * G * (mmu * mmu - 2.0 * mmu * pf[2].E()) * mmu * pf[2].E();
         double coeff = 1.0 / (2.0 * mmu);

         return M2 * coeff;
      };
   };

   // make an object of your generator
   Generator generator;

   // silent mode
   generator.SetChat(0);

   // set up initial particle/blob that decays:

   // center of mass energy
   const double mmu = 0.1057; // muon rest mass GeV

   // CM 4-momentum - blob that initiates decay
   TLorentzVector pbCM;

   // set blob CM energy to mass of the muon with zero momentum - rest frame of particle
   pbCM.SetPxPyPzE(0.0, 0.0, 0.0, mmu);

   // set up outgoing particles configuration:

   // number of outgoing particles
   const int Nop = 3;

   // PDGID outgoing paricles (masses are taken from PDG table)
   // full list: https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
   int idOut[Nop] = {11, 14, 12}; // e, \nu_\mi, \nu_e

   // out particles array for storing their 4-vectors
   TLorentzVector pf[Nop];

   // masses of products - from PDGDatabase
   double mass[Nop];

   // fill in the table with masses of particles using their codes
   for (int i = 0; i < Nop; i++) {
      mass[i] = PDGDatabese->GetParticle(idOut[i])->Mass();
   }

   // set Generator to decay configuration
   generator.SetDecay(pbCM, Nop, mass);

   long NevTot = 10e4; // Total number of events to generate

   // GENRATION LOOP
   for (long loop = 0; loop < NevTot; loop++) {

      // Generate event

      // make decay
      generator.Generate();

      // get out particles into pf[] array - order as in mass[] array
      for (int i = 0; i < Nop; i++) {
         pf[i] = *(generator.GetDecay(i));
      }
   }
   // END of GENRATION LOOP

   // generator.Finalize();

   // get integral and error from Generator
   Double_t integral, error;
   generator.GetIntegMC(integral, error);

   // Theoretical prediction from (see Problem 5.3 in M.D. Schwatz, 'Quantum Field Theory and The Standard Model',
   // Cambridge  2014)
   const double G = 1.166e-5;
   double th_integral = (G * G * pow(mmu, 5)) / (192.0 * pow(M_PI, 3));

   ASSERT_TRUE(abs(integral - th_integral) < error);
};
