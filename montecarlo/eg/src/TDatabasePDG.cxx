// @(#)root/eg:$Id$
// Author: Pasha Murat   12/02/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "RConfigure.h"
#include "TROOT.h"
#include "TEnv.h"
#include "THashList.h"
#include "TExMap.h"
#include "TSystem.h"
#include "TDatabasePDG.h"
#include "TDecayChannel.h"
#include "TParticlePDG.h"
#include <stdlib.h>


/** \class TDatabasePDG
    \ingroup eg

Particle database manager class

This manager creates a list of particles which by default is
initialised from with the constants used by PYTHIA6 (plus some
other particles added). See definition and the format of the default
particle list in $ROOTSYS/etc/pdg_table.txt

There are 2 ways of redefining the name of the file containing the
particle properties

-# One can define the name in .rootrc file:

Root.DatabasePDG: $(HOME)/my_pdg_table.txt

-# One can use TDatabasePDG::ReadPDGTable method explicitly:

   - TDatabasePDG *pdg = new TDatabasePDG();
   - pdg->ReadPDGtable(filename)

See TParticlePDG for the description of a static particle properties.
See TParticle    for the description of a dynamic particle particle.

The current default pdg_table file displays lifetime 0 for some unstable particles.

*/

ClassImp(TDatabasePDG);

////////////////////////////////////////////////////////////////////////////////
/// Static function holding the instance.

TDatabasePDG** GetInstancePtr()
{
   static TDatabasePDG* fgInstance = nullptr;
   return &fgInstance;
}

////////////////////////////////////////////////////////////////////////////////
/// Create PDG database. Initialization of the DB has to be done via explicit
/// call to ReadDataBasePDG (also done by GetParticle methods)

TDatabasePDG::TDatabasePDG(): TNamed("PDGDB","The PDG particle data base")
{
   fParticleList  = 0;
   fPdgMap        = 0;
   fListOfClasses = 0;
   auto fgInstance = GetInstancePtr();
   if (*fgInstance != nullptr) {
      Warning("TDatabasePDG", "object already instantiated");
   } else {
      *fgInstance = this;
      gROOT->GetListOfSpecials()->Add(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup the PDG database.

TDatabasePDG::~TDatabasePDG()
{
   if (fParticleList) {
      fParticleList->Delete();
      delete fParticleList;    // this deletes all objects in the list
      if (fPdgMap) delete fPdgMap;
   }
                                // classes do not own particles...
   if (fListOfClasses) {
      fListOfClasses->Delete();
      delete fListOfClasses;
   }
   if (gROOT && !gROOT->TestBit(TObject::kInvalidObject))
      gROOT->GetListOfSpecials()->Remove(this);
   auto fgInstance = GetInstancePtr();
   *fgInstance = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
///static function

TDatabasePDG* TDatabasePDG::Instance()
{
   auto fgInstance = GetInstancePtr();
   if (*fgInstance == nullptr) {
      // Constructor creates a new instance, inits fgInstance.
      new TDatabasePDG();
   }
   return *fgInstance;
}

////////////////////////////////////////////////////////////////////////////////
/// Build fPdgMap mapping pdg-code to particle.
///
/// Initial size is set so as to be able to hold at least 600
/// particles: 521 in default table, ALICE adds 54 more.
/// To be revisited after LHC discovers SUSY.

void TDatabasePDG::BuildPdgMap() const
{
   fPdgMap = new TExMap(4*TMath::Max(600, fParticleList->GetEntries())/3 + 3);
   TIter next(fParticleList);
   TParticlePDG *p;
   while ((p = (TParticlePDG*)next())) {
      fPdgMap->Add((Long_t)p->PdgCode(), (Long_t)p);
   }
}

////////////////////////////////////////////////////////////////////////////////
///
///  Particle definition normal constructor. If the particle is set to be
///  stable, the decay width parameter does have no meaning and can be set to
///  any value. The parameters granularity, LowerCutOff and HighCutOff are
///  used for the construction of the mean free path look up tables. The
///  granularity will be the number of logwise energy points for which the
///  mean free path will be calculated.
///

TParticlePDG* TDatabasePDG::AddParticle(const char *name, const char *title,
                                        Double_t mass, Bool_t stable,
                                        Double_t width, Double_t charge,
                                        const char* ParticleClass,
                                        Int_t PDGcode,
                                        Int_t Anti,
                                        Int_t TrackingCode)
{
   TParticlePDG* old = GetParticle(PDGcode);

   if (old) {
      printf(" *** TDatabasePDG::AddParticle: particle with PDGcode=%d already defined\n",PDGcode);
      return 0;
   }

   TParticlePDG* p = new TParticlePDG(name, title, mass, stable, width,
                                     charge, ParticleClass, PDGcode, Anti,
                                     TrackingCode);
   fParticleList->Add(p);
   if (fPdgMap)
      fPdgMap->Add((Long_t)PDGcode, (Long_t)p);

   TParticleClassPDG* pclass = GetParticleClass(ParticleClass);

   if (!pclass) {
      pclass = new TParticleClassPDG(ParticleClass);
      fListOfClasses->Add(pclass);
   }

   pclass->AddParticle(p);

   return p;
}

////////////////////////////////////////////////////////////////////////////////
/// assuming particle has already been defined

TParticlePDG* TDatabasePDG::AddAntiParticle(const char* Name, Int_t PdgCode)
{
   TParticlePDG* old = GetParticle(PdgCode);

   if (old) {
      printf(" *** TDatabasePDG::AddAntiParticle: can't redefine parameters\n");
      return NULL;
   }

   Int_t pdg_code  = abs(PdgCode);
   TParticlePDG* p = GetParticle(pdg_code);

   if (!p) {
      printf(" *** TDatabasePDG::AddAntiParticle: particle with pdg code %d not known\n", pdg_code);
      return NULL;
   }

   TParticlePDG* ap = AddParticle(Name,
                                  Name,
                                  p->Mass(),
                                  1,
                                  p->Width(),
                                  -p->Charge(),
                                  p->ParticleClass(),
                                  PdgCode,
                                  1,
                                  p->TrackingCode());
   return ap;
}


////////////////////////////////////////////////////////////////////////////////
///
///  Get a pointer to the particle object according to the name given
///

TParticlePDG *TDatabasePDG::GetParticle(const char *name) const
{
   if (fParticleList == 0)  ((TDatabasePDG*)this)->ReadPDGTable();

   TParticlePDG *def = (TParticlePDG *)fParticleList->FindObject(name);
//     if (!def) {
//        Error("GetParticle","No match for %s exists!",name);
//     }
   return def;
}

////////////////////////////////////////////////////////////////////////////////
///
///  Get a pointer to the particle object according to the MC code number
///

TParticlePDG *TDatabasePDG::GetParticle(Int_t PDGcode) const
{
   if (fParticleList == 0)  ((TDatabasePDG*)this)->ReadPDGTable();
   if (fPdgMap       == 0)  BuildPdgMap();

   return (TParticlePDG*) (Long_t)fPdgMap->GetValue((Long_t)PDGcode);
}

////////////////////////////////////////////////////////////////////////////////
/// Print contents of PDG database.

void TDatabasePDG::Print(Option_t *option) const
{
   if (fParticleList == 0)  ((TDatabasePDG*)this)->ReadPDGTable();

   TIter next(fParticleList);
   TParticlePDG *p;
   while ((p = (TParticlePDG *)next())) {
      p->Print(option);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Converts Geant3 particle codes to PDG convention. (Geant4 uses
/// PDG convention already)
/// Source: BaBar User Guide, Neil I. Geddes,
///
/// see <A href="http://www.slac.stanford.edu/BFROOT/www/Computing/Environment/NewUser/htmlbug/node51.html"> Conversion table</A>
///
/// with some fixes by PB, marked with (PB) below. Checked against
/// PDG listings from 2000.
///
/// Paul Balm, Nov 19, 2001

Int_t TDatabasePDG::ConvertGeant3ToPdg(Int_t Geant3number) const
{


   switch(Geant3number) {

      case 1   : return 22;       // photon
      case 25  : return -2112;    // anti-neutron
      case 2   : return -11;      // e+
      case 26  : return -3122;    // anti-Lambda
      case 3   : return 11;       // e-
      case 27  : return -3222;    // Sigma-
      case 4   : return 12;       // e-neutrino (NB: flavour undefined by Geant)
      case 28  : return -3212;    // Sigma0
      case 5   : return -13;      // mu+
      case 29  : return -3112;    // Sigma+ (PB)*/
      case 6   : return 13;       // mu-
      case 30  : return -3322;    // Xi0
      case 7   : return 111;      // pi0
      case 31  : return -3312;    // Xi+
      case 8   : return 211;      // pi+
      case 32  : return -3334;    // Omega+ (PB)
      case 9   : return -211;     // pi-
      case 33  : return -15;      // tau+
      case 10  : return 130;      // K long
      case 34  : return 15;       // tau-
      case 11  : return 321;      // K+
      case 35  : return 411;      // D+
      case 12  : return -321;     // K-
      case 36  : return -411;     // D-
      case 13  : return 2112;     // n
      case 37  : return 421;      // D0
      case 14  : return 2212;     // p
      case 38  : return -421;     // D0
      case 15  : return -2212;    // anti-proton
      case 39  : return 431;      // Ds+
      case 16  : return 310;      // K short
      case 40  : return -431;     // anti Ds-
      case 17  : return 221;      // eta
      case 41  : return 4122;     // Lamba_c+
      case 18  : return 3122;     // Lambda
      case 42  : return 24;       // W+
      case 19  : return 3222;     // Sigma+
      case 43  : return -24;      // W-
      case 20  : return 3212;     // Sigma0
      case 44  : return 23;       // Z
      case 21  : return 3112;     // Sigma-
      case 45  : return 0;        // deuteron
      case 22  : return 3322;     // Xi0
      case 46  : return 0;        // triton
      case 23  : return 3312;     // Xi-
      case 47  : return 0;        // alpha
      case 24  : return 3334;     // Omega- (PB)
      case 48  : return 0;        // G nu ? PDG ID 0 is undefined

      default  : return 0;

   }
}

////////////////////////////////////////////////////////////////////////////////
/// Converts pdg code to geant3 id

Int_t TDatabasePDG::ConvertPdgToGeant3(Int_t pdgNumber) const
{
   switch(pdgNumber) {

      case   22     : return  1;    // photon
      case   -2112  : return  25;   // anti-neutron
      case   -11    : return  2;    // e+
      case   -3122  : return  26;   // anti-Lambda
      case   11     : return  3;    // e-
      case   -3222  : return  27;   // Sigma-
      case   12     : return  4;    // e-neutrino (NB: flavour undefined by Geant)
      case   -3212  : return  28;   // Sigma0
      case   -13    : return  5;    // mu+
      case   -3112  : return  29;   // Sigma+ (PB)*/
      case   13     : return  6;    // mu-
      case   -3322  : return  30;   // Xi0
      case   111    : return  7;    // pi0
      case   -3312  : return  31;   // Xi+
      case   211    : return  8;    // pi+
      case   -3334  : return  32;   // Omega+ (PB)
      case   -211   : return  9;    // pi-
      case   -15    : return  33;   // tau+
      case   130    : return  10;   // K long
      case   15     : return  34;   // tau-
      case   321    : return  11;   // K+
      case   411    : return  35;   // D+
      case   -321   : return  12;   // K-
      case   -411   : return  36;   // D-
      case   2112   : return  13;   // n
      case   421    : return  37;   // D0
      case   2212   : return  14;   // p
      case   -421   : return  38;   // D0
      case   -2212  : return  15;   // anti-proton
      case   431    : return  39;   // Ds+
      case   310    : return  16;   // K short
      case   -431   : return  40;   // anti Ds-
      case   221    : return  17;   // eta
      case   4122   : return  41;   // Lamba_c+
      case   3122   : return  18;   // Lambda
      case   24     : return  42;   // W+
      case   3222   : return  19;   // Sigma+
      case   -24    : return  43;   // W-
      case   3212   : return  20;   // Sigma0
      case   23     : return  44;   // Z
      case   3112   : return  21;   // Sigma-
      case   3322   : return  22;   // Xi0
      case   3312   : return  23;   // Xi-
      case   3334   : return  24;   // Omega- (PB)

      default  : return 0;

   }
}

////////////////////////////////////////////////////////////////////////////////
///
///  Converts the ISAJET Particle number into the PDG MC number
///

Int_t TDatabasePDG::ConvertIsajetToPdg(Int_t isaNumber) const
{
   switch (isaNumber) {
      case     1 : return     2; //     UP        .30000E+00       .67
      case    -1 : return    -2; //     UB        .30000E+00      -.67
      case     2 : return     1; //     DN        .30000E+00      -.33
      case    -2 : return    -1; //     DB        .30000E+00       .33
      case     3 : return     3; //     ST        .50000E+00      -.33
      case    -3 : return    -3; //     SB        .50000E+00       .33
      case     4 : return     4; //     CH        .16000E+01       .67
      case    -4 : return    -4; //     CB        .16000E+01      -.67
      case     5 : return     5; //     BT        .49000E+01      -.33
      case    -5 : return    -5; //     BB        .49000E+01       .33
      case     6 : return     6; //     TP        .17500E+03       .67
      case    -6 : return    -6; //     TB        .17500E+03      -.67
      case     9 : return    21; //     GL       0.               0.00
      case    80 : return    24; //     W+        SIN2W=.23       1.00
      case   -80 : return   -24; //     W-        SIN2W=.23      -1.00
      case    90 : return    23; //     Z0        SIN2W=.23       0.00
      case   230 : return   311; //     K0        .49767E+00      0.00
      case  -230 : return  -311; //     AK0       .49767E+00      0.00
      case   330 : return   331; //     ETAP      .95760E+00      0.00
      case   340 : return     0; //     F-        .20300E+01     -1.00
      case  -340 : return     0; //     F+        .20300E+01      1.00
      case   440 : return   441; //     ETAC      .29760E+01      0.00
      case   111 : return   113; //     RHO0      .77000E+00      0.00
      case   121 : return   213; //     RHO+      .77000E+00      1.00
      case  -121 : return  -213; //     RHO-      .77000E+00     -1.00
      case   221 : return   223; //     OMEG      .78260E+00      0.00
      case   131 : return   323; //     K*+       .88810E+00      1.00
      case  -131 : return  -323; //     K*-       .88810E+00     -1.00
      case   231 : return   313; //     K*0       .89220E+00      0.00
      case  -231 : return  -313; //     AK*0      .89220E+00      0.00
      case   331 : return   333; //     PHI       .10196E+01      0.00
      case  -140 : return   421; //     D0
      case   140 : return  -421; //     D0 bar
      case   141 : return  -423; //     AD*0      .20060E+01      0.00
      case  -141 : return   423; //     D*0       .20060E+01      0.00
      case  -240 : return  -411; //     D+
      case   240 : return   411; //     D-
      case   241 : return  -413; //     D*-       .20086E+01     -1.00
      case  -241 : return   413; //     D*+       .20086E+01      1.00
      case   341 : return     0; //     F*-       .21400E+01     -1.00
      case  -341 : return     0; //     F*+       .21400E+01      1.00
      case   441 : return   443; //     JPSI      .30970E+01      0.00

                                        // B-mesons, Bc still missing
      case   250 : return   511; // B0
      case  -250 : return  -511; // B0 bar
      case   150 : return   521; // B+
      case  -150 : return  -521; // B-
      case   350 : return   531; // Bs  0
      case  -350 : return  -531; // Bs  bar
      case   351 : return   533; // Bs* 0
      case  -351 : return  -533; // Bs* bar
      case   450 : return   541; // Bc  +
      case  -450 : return  -541; // Bc  bar

      case  1140 : return  4222; //     SC++      .24300E+01      2.00
      case -1140 : return -4222; //     ASC--     .24300E+01     -2.00
      case  1240 : return  4212; //     SC+       .24300E+01      1.00
      case -1240 : return -4212; //     ASC-      .24300E+01     -1.00
      case  2140 : return  4122; //     LC+       .22600E+01      1.00
      case -2140 : return -4122; //     ALC-      .22600E+01     -1.00
      case  2240 : return  4112; //     SC0       .24300E+01      0.00
      case -2240 : return -4112; //     ASC0      .24300E+01      0.00
      case  1340 : return     0; //     USC.      .25000E+01      1.00
      case -1340 : return     0; //     AUSC.     .25000E+01     -1.00
      case  3140 : return     0; //     SUC.      .24000E+01      1.00
      case -3140 : return     0; //     ASUC.     .24000E+01     -1.00
      case  2340 : return     0; //     DSC.      .25000E+01      0.00
      case -2340 : return     0; //     ADSC.     .25000E+01      0.00
      case  3240 : return     0; //     SDC.      .24000E+01      0.00
      case -3240 : return     0; //     ASDC.     .24000E+01      0.00
      case  3340 : return     0; //     SSC.      .26000E+01      0.00
      case -3340 : return     0; //     ASSC.     .26000E+01      0.00
      case  1440 : return     0; //     UCC.      .35500E+01      2.00
      case -1440 : return     0; //     AUCC.     .35500E+01     -2.00
      case  2440 : return     0; //     DCC.      .35500E+01      1.00
      case -2440 : return     0; //     ADCC.     .35500E+01     -1.00
      case  3440 : return     0; //     SCC.      .37000E+01      1.00
      case -3440 : return     0; //     ASCC.     .37000E+01     -1.00
      case  1111 : return  2224; //     DL++      .12320E+01      2.00
      case -1111 : return -2224; //     ADL--     .12320E+01     -2.00
      case  1121 : return  2214; //     DL+       .12320E+01      1.00
      case -1121 : return -2214; //     ADL-      .12320E+01     -1.00
      case  1221 : return  2114; //     DL0       .12320E+01      0.00
      case -1221 : return -2114; //     ADL0      .12320E+01      0.00
      case  2221 : return   1114; //     DL-       .12320E+01     -1.00
      case -2221 : return -1114; //     ADL+      .12320E+01      1.00
      case  1131 : return  3224; //     S*+       .13823E+01      1.00
      case -1131 : return -3224; //     AS*-      .13823E+01     -1.00
      case  1231 : return  3214; //     S*0       .13820E+01      0.00
      case -1231 : return -3214; //     AS*0      .13820E+01      0.00
      case  2231 : return  3114; //     S*-       .13875E+01     -1.00
      case -2231 : return -3114; //     AS*+      .13875E+01      1.00
      case  1331 : return  3324; //     XI*0      .15318E+01      0.00
      case -1331 : return -3324; //     AXI*0     .15318E+01      0.00
      case  2331 : return  3314; //     XI*-      .15350E+01     -1.00
      case -2331 : return -3314; //     AXI*+     .15350E+01      1.00
      case  3331 : return  3334; //     OM-       .16722E+01     -1.00
      case -3331 : return -3334; //     AOM+      .16722E+01      1.00
      case  1141 : return     0; //     UUC*      .26300E+01      2.00
      case -1141 : return     0; //     AUUC*     .26300E+01     -2.00
      case  1241 : return     0; //     UDC*      .26300E+01      1.00
      case -1241 : return     0; //     AUDC*     .26300E+01     -1.00
      case  2241 : return     0; //     DDC*      .26300E+01      0.00
      case -2241 : return     0; //     ADDC*     .26300E+01      0.00
      case  1341 : return     0; //     USC*      .27000E+01      1.00
      case -1341 : return     0; //     AUSC*     .27000E+01     -1.00
      case  2341 : return     0; //     DSC*      .27000E+01      0.00
      case -2341 : return     0; //     ADSC*     .27000E+01      0.00
      case  3341 : return     0; //     SSC*      .28000E+01      0.00
      case -3341 : return     0; //     ASSC*     .28000E+01      0.00
      case  1441 : return     0; //     UCC*      .37500E+01      2.00
      case -1441 : return     0; //     AUCC*     .37500E+01     -2.00
      case  2441 : return     0; //     DCC*      .37500E+01      1.00
      case -2441 : return     0; //     ADCC*     .37500E+01     -1.00
      case  3441 : return     0; //     SCC*      .39000E+01      1.00
      case -3441 : return     0; //     ASCC*     .39000E+01     -1.00
      case  4441 : return     0; //     CCC*      .48000E+01      2.00
      case -4441 : return     0; //     ACCC*     .48000E+01     -2.00
      case    10 : return    22; // Photon
      case    12 : return    11; // Electron
      case   -12 : return   -11; // Positron
      case    14 : return    13; // Muon-
      case   -14 : return   -13; // Muon+
      case    16 : return    15; // Tau-
      case   -16 : return   -15; // Tau+
      case    11 : return    12; // Neutrino e
      case   -11 : return   -12; // Anti Neutrino e
      case    13 : return    14; // Neutrino Muon
      case   -13 : return   -14; // Anti Neutrino Muon
      case    15 : return    16; // Neutrino Tau
      case   -15 : return   -16; // Anti Neutrino Tau
      case   110 : return   111; // Pion0
      case   120 : return   211; // Pion+
      case  -120 : return  -211; // Pion-
      case   220 : return   221; // Eta
      case   130 : return   321; // Kaon+
      case  -130 : return  -321; // Kaon-
      case   -20 : return   130; // Kaon Long
      case    20 : return   310; // Kaon Short

                                        // baryons
      case  1120 : return  2212; // Proton
      case -1120 : return -2212; // Anti Proton
      case  1220 : return  2112; // Neutron
      case -1220 : return -2112; // Anti Neutron
      case  2130 : return  3122; // Lambda
      case -2130 : return -3122; // Lambda bar
      case  1130 : return  3222; // Sigma+
      case -1130 : return -3222; // Sigma bar -
      case  1230 : return  3212; // Sigma0
      case -1230 : return -3212; // Sigma bar 0
      case  2230 : return  3112; // Sigma-
      case -2230 : return -3112; // Sigma bar +
      case  1330 : return  3322; // Xi0
      case -1330 : return -3322; // Xi bar 0
      case  2330 : return  3312; // Xi-
      case -2330 : return -3312; // Xi bar +
      default :    return 0;      // isajet or pdg number does not exist
   }
}

////////////////////////////////////////////////////////////////////////////////
/// read list of particles from a file
/// if the particle list does not exist, it is created, otherwise
/// particles are added to the existing list
/// See $ROOTSYS/etc/pdg_table.txt to see the file format

void TDatabasePDG::ReadPDGTable(const char *FileName)
{
   if (fParticleList == 0) {
      fParticleList  = new THashList;
      fListOfClasses = new TObjArray;
   }

   TString default_name;
   const char *fn;

   if (!FileName[0]) {
      default_name = "pdg_table.txt";
      gSystem->PrependPathName(TROOT::GetEtcDir(), default_name);
      fn = gEnv->GetValue("Root.DatabasePDG", default_name.Data());
   } else {
      fn = FileName;
   }

   FILE* file = fopen(fn,"r");
   if (file == 0) {
      Error("ReadPDGTable","Could not open PDG particle file %s",fn);
      return;
   }

   char      c[512];
   Int_t     class_number, anti, isospin, i3, spin, tracking_code;
   Int_t     ich, kf, nch, charge;
   char      name[30], class_name[30];
   Double_t  mass, width, branching_ratio;
   Int_t     dau[20];

   Int_t     idecay, decay_type, flavor, ndau, stable;

   Int_t input;
   while ( (input=getc(file)) != EOF) {
      c[0] = input;
      if (c[0] != '#') {
         ungetc(c[0],file);
         // read channel number
         // coverity [secure_coding : FALSE]
         if (fscanf(file,"%i",&ich)) {;}
         // coverity [secure_coding : FALSE]
         if (fscanf(file,"%s",name  )) {;}
         // coverity [secure_coding : FALSE]
         if (fscanf(file,"%i",&kf   )) {;}
         // coverity [secure_coding : FALSE]
         if (fscanf(file,"%i",&anti )) {;}

         if (kf < 0) {
            AddAntiParticle(name,kf);
            // nothing more on this line
            if (fgets(c,200,file)) {;}
         } else {
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&class_number)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%s",class_name)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&charge)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%le",&mass)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%le",&width)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&isospin)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&i3)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&spin)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&flavor)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&tracking_code)) {;}
            // coverity [secure_coding : FALSE]
            if (fscanf(file,"%i",&nch)) {;}
            // nothing more on this line
            if (fgets(c,200,file)) {;}
            if (width > 1e-10) stable = 0;
            else               stable = 1;

            // create particle

            TParticlePDG* part = AddParticle(name,
                                             name,
                                             mass,
                                             stable,
                                             width,
                                             charge,
                                             class_name,
                                             kf,
                                             anti,
                                             tracking_code);

            if (nch) {
               // read in decay channels
               ich = 0;
               Int_t c_input = 0;
               while ( ((c_input=getc(file)) != EOF) && (ich <nch)) {
                  c[0] = c_input;
                  if (c[0] != '#') {
                     ungetc(c[0],file);

                     // coverity [secure_coding : FALSE]
                     if (fscanf(file,"%i",&idecay)) {;}
                     // coverity [secure_coding : FALSE]
                     if (fscanf(file,"%i",&decay_type)) {;}
                     // coverity [secure_coding : FALSE]
                     if (fscanf(file,"%le",&branching_ratio)) {;}
                     // coverity [secure_coding : FALSE]
                     if (fscanf(file,"%i",&ndau)) {;}
                     for (int idau=0; idau<ndau; idau++) {
                        // coverity [secure_coding : FALSE]
                        if (fscanf(file,"%i",&dau[idau])) {;}
                     }
                     // add decay channel

                     if (part) part->AddDecayChannel(decay_type,branching_ratio,ndau,dau);
                     ich++;
                  }
                  // skip end of line
                  if (fgets(c,200,file)) {;}
               }
            }
         }
      } else {
         // skip end of line
         if (fgets(c,200,file)) {;}
      }
   }
   // in the end loop over the antiparticles and
   // define their decay lists
   TIter it(fParticleList);

   Int_t code[20];
   TParticlePDG  *ap, *p, *daughter;
   TDecayChannel *dc;

   while ((p = (TParticlePDG*) it.Next())) {

      // define decay channels for antiparticles
      if (p->PdgCode() < 0) {
         ap = GetParticle(-p->PdgCode());
         if (!ap) continue;
         nch = ap->NDecayChannels();
         for (ich=0; ich<nch; ich++) {
            dc = ap->DecayChannel(ich);
            if (!dc) continue;
            ndau = dc->NDaughters();
            for (int i=0; i<ndau; i++) {
               // conserve CPT

               code[i] = dc->DaughterPdgCode(i);
               daughter = GetParticle(code[i]);
               if (daughter && daughter->AntiParticle()) {
                  // this particle does have an
                  // antiparticle
                  code[i] = -code[i];
               }
            }
            p->AddDecayChannel(dc->MatrixElementCode(),
                               dc->BranchingRatio(),
                               dc->NDaughters(),
                               code);
         }
         p->SetAntiParticle(ap);
         ap->SetAntiParticle(p);
      }
   }

   fclose(file);
   return;
}


////////////////////////////////////////////////////////////////////////////////
///browse data base

void TDatabasePDG::Browse(TBrowser* b)
{
   if (fListOfClasses ) fListOfClasses->Browse(b);
}


////////////////////////////////////////////////////////////////////////////////
/// write contents of the particle DB into a file

Int_t TDatabasePDG::WritePDGTable(const char *filename)
{
   if (fParticleList == 0) {
      Error("WritePDGTable","Do not have a valid PDG particle list;"
                            " consider loading it with ReadPDGTable first.");
      return -1;
   }

   FILE *file = fopen(filename,"w");
   if (file == 0) {
      Error("WritePDGTable","Could not open PDG particle file %s",filename);
      return -1;
   }

   fprintf(file,"#--------------------------------------------------------------------\n");
   fprintf(file,"#    i   NAME.............  KF AP   CLASS      Q        MASS     WIDTH  2*I+1 I3 2*S+1 FLVR TrkCod N(dec)\n");
   fprintf(file,"#--------------------------------------------------------------------\n");

   Int_t nparts=fParticleList->GetEntries();
   for(Int_t i=0;i<nparts;++i) {
      TParticlePDG *p = dynamic_cast<TParticlePDG*>(fParticleList->At(i));
      if(!p) continue;

      Int_t ich=i+1;
      Int_t kf=p->PdgCode();
      fprintf(file,"%5i %-20s %- 6i ", ich, p->GetName(), kf);

      Int_t anti=p->AntiParticle() ? 1:0;
      if(kf<0) {
         for(Int_t j=0;j<nparts;++j) {
            TParticlePDG *dummy = dynamic_cast<TParticlePDG*>(fParticleList->At(j));
            if(dummy==p->AntiParticle()) {
               anti=j+1;
               break;
            }
         }
         fprintf(file,"%i 0\n",anti);
         continue;
      }

      fprintf(file,"%i ",anti);
      fprintf(file,"%i ",100);
      fprintf(file,"%s ",p->ParticleClass());
      fprintf(file,"% i ",(Int_t)p->Charge());
      fprintf(file,"%.5le ",p->Mass());
      fprintf(file,"%.5le ",p->Width());
      fprintf(file,"%i ",(Int_t)p->Isospin());
      fprintf(file,"%i ",(Int_t)p->I3());
      fprintf(file,"%i ",(Int_t)p->Spin());
      fprintf(file,"%i ",-1);
      fprintf(file,"%i ",p->TrackingCode());
      Int_t nch=p->NDecayChannels();
      fprintf(file,"%i\n",nch);
      if(nch==0) {
         continue;
      }
      fprintf(file,"#----------------------------------------------------------------------\n");
      fprintf(file,"#    decay  type(PY6)    BR     Nd         daughters(codes, then names)\n");
      fprintf(file,"#----------------------------------------------------------------------\n");
      for(Int_t j=0;j<nch;++j) {
         TDecayChannel *dc=p->DecayChannel(j);
         if (!dc) continue;
         fprintf(file,"%9i   ",dc->Number()+1);
         fprintf(file,"%3i   ",dc->MatrixElementCode());
         fprintf(file,"%.5le  ",dc->BranchingRatio());
         Int_t ndau=dc->NDaughters();
         fprintf(file,"%3i       ",ndau);
         for (int idau=0; idau<ndau; idau++) {
            fprintf(file,"%- 6i ",dc->DaughterPdgCode(idau));
         }
         for (int idau=0; idau<ndau; idau++) {
            TParticlePDG *dummy=GetParticle(dc->DaughterPdgCode(idau));
            if(dummy)
               fprintf(file,"%-10s ",dummy->GetName());
            else
               fprintf(file,"%-10s ","???");
         }
         fprintf(file,"\n");
      }
   }
   fclose(file);
   return nparts;
}
