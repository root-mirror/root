// @(#)root/geom:$Id$
// Author: Andrei Gheata   17/06/04

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoElement
\ingroup Geometry_classes
Base class for chemical elements
*/

/** \class TGeoElementRN
\ingroup Geometry_classes
Class representing a radionuclidevoid TGeoManager::SetDefaultRootUnits()
{
   if ( fgDefaultUnits == kRootUnits )   {
      return;
   }
   else if ( gGeometryLocked )    {
      TError::Fatal("TGeoManager","The system of units may only be changed once BEFORE any elements and materials are created!");
   }
   fgDefaultUnits = kRootUnits;
}

*/

/** \class TGeoElemIter
\ingroup Geometry_classes
Iterator for decay branches
*/

/** \class TGeoDecayChannel
\ingroup Geometry_classes
A decay channel for a radionuclide
*/

/** \class TGeoElementTable
\ingroup Geometry_classes
Table of elements
*/

#include "RConfigure.h"

#include <fstream>
#include <iomanip>

#include "TSystem.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TVirtualGeoPainter.h"
#include "TGeoManager.h"
#include "TGeoElement.h"
#include "TMath.h"
#include "TGeoPhysicalConstants.h"
#include "TGeant4PhysicalConstants.h"

// statics and globals
static const Int_t gMaxElem  = 110;
static const Int_t gMaxLevel = 8;
static const Int_t gMaxDecay = 15;

static const char gElName[gMaxElem][3] = {
          "H ","He","Li","Be","B ","C ","N ","O ","F ","Ne","Na","Mg",
          "Al","Si","P ","S ","Cl","Ar","K ","Ca","Sc","Ti","V ","Cr",
          "Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
          "Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
          "In","Sn","Sb","Te","I ","Xe","Cs","Ba","La","Ce","Pr","Nd",
          "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",
          "Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po",
          "At","Rn","Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm",
          "Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs",
          "Mt","Ds" };

static const char *gDecayName[gMaxDecay+1] = {
   "2BetaMinus", "BetaMinus", "NeutronEm", "ProtonEm", "Alpha", "ECF",
   "ElecCapt", "IsoTrans", "I", "SpontFiss", "2ProtonEm", "2NeutronEm",
   "2Alpha", "Carbon12", "Carbon14", "Stable" };

static const Int_t gDecayDeltaA[gMaxDecay] = {
              0,           0,          -1,         -1,       -4,
            -99,           0,           0,        -99,      -99,
             -2,          -2,          -8,        -12,      -14 };

static const Int_t gDecayDeltaZ[gMaxDecay] = {
              2,           1,           0,         -1,       -2,
            -99,          -1,           0,        -99,      -99,
             -2,           0,          -4,         -6,       -6 };
static const char gLevName[gMaxLevel]=" mnpqrs";

ClassImp(TGeoElement);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoElement::TGeoElement()
{
   TGeoManager::SetDefaultUnits(TGeoManager::GetDefaultUnits()); // Ensure nobody changes the units afterwards
   SetDefined(kFALSE);
   SetUsed(kFALSE);
   fZ = 0;
   fN = 0;
   fNisotopes = 0;
   fA = 0.0;
   fIsotopes = NULL;
   fAbundances = NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// Obsolete constructor

TGeoElement::TGeoElement(const char *name, const char *title, Int_t z, Double_t a)
            :TNamed(name, title)
{
   TGeoManager::SetDefaultUnits(TGeoManager::GetDefaultUnits()); // Ensure nobody changes the units afterwards
   SetDefined(kFALSE);
   SetUsed(kFALSE);
   fZ = z;
   fN = Int_t(a);
   fNisotopes = 0;
   fA = a;
   fIsotopes = NULL;
   fAbundances = NULL;
   ComputeDerivedQuantities();
}

////////////////////////////////////////////////////////////////////////////////
/// Element having isotopes.

TGeoElement::TGeoElement(const char *name, const char *title, Int_t nisotopes)
            :TNamed(name, title)
{
   TGeoManager::SetDefaultUnits(TGeoManager::GetDefaultUnits()); // Ensure nobody changes the units afterwards
   SetDefined(kFALSE);
   SetUsed(kFALSE);
   fZ = 0;
   fN = 0;
   fNisotopes = nisotopes;
   fA = 0.0;
   fIsotopes = new TObjArray(nisotopes);
   fAbundances = new Double_t[nisotopes];
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TGeoElement::TGeoElement(const char *name, const char *title, Int_t z, Int_t n, Double_t a)
            :TNamed(name, title)
{
   TGeoManager::SetDefaultUnits(TGeoManager::GetDefaultUnits()); // Ensure nobody changes the units afterwards
   SetDefined(kFALSE);
   SetUsed(kFALSE);
   fZ = z;
   fN = n;
   fNisotopes = 0;
   fA = a;
   fIsotopes = NULL;
   fAbundances = NULL;
   ComputeDerivedQuantities();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoElement::~TGeoElement()
{
   delete fIsotopes;
   delete [] fAbundances;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate properties for an atomic number

void TGeoElement::ComputeDerivedQuantities()
{
   // Radiation Length
   ComputeCoulombFactor();
   ComputeLradTsaiFactor();
}
////////////////////////////////////////////////////////////////////////////////
/// Compute Coulomb correction factor (Phys Rev. D50 3-1 (1994) page 1254)

void TGeoElement::ComputeCoulombFactor()
{
   static constexpr Double_t k1 = 0.0083 , k2 = 0.20206 ,k3 = 0.0020 , k4 = 0.0369;
   Double_t fsc = TGeoManager::kRootUnits == TGeoManager::GetDefaultUnits()
     ? TGeoUnit::fine_structure_const : TGeant4Unit::fine_structure_const;
   Double_t az2 = (fsc*fZ)*(fsc*fZ);
   Double_t az4 = az2 * az2;

   fCoulomb = (k1*az4 + k2 + 1./(1.+az2))*az2 - (k3*az4 + k4)*az4;
}
////////////////////////////////////////////////////////////////////////////////
/// Compute Tsai's Expression for the Radiation Length (Phys Rev. D50 3-1 (1994) page 1254)

void TGeoElement::ComputeLradTsaiFactor()
{
   static constexpr Double_t Lrad_light[]  = {5.31  , 4.79  , 4.74 ,  4.71} ;
   static constexpr Double_t Lprad_light[] = {6.144 , 5.621 , 5.805 , 5.924} ;

   fRadTsai = 0.0;
   if (fZ == 0) return;
   const Double_t logZ3 = TMath::Log(fZ)/3.;

   Double_t Lrad, Lprad;
   Double_t alpha_rcl2 = TGeoManager::kRootUnits == TGeoManager::GetDefaultUnits()
     ? TGeoUnit::alpha_rcl2 : TGeant4Unit::alpha_rcl2;
   Int_t iz = static_cast<Int_t>(fZ+0.5) - 1 ; // The static cast comes from G4lrint
   static const Double_t log184  = TMath::Log(184.15);
   static const Double_t log1194 = TMath::Log(1194.);
   if (iz <= 3) { Lrad = Lrad_light[iz] ;  Lprad = Lprad_light[iz] ; }
   else { Lrad = log184 - logZ3 ; Lprad = log1194 - 2*logZ3;}

   fRadTsai = 4*alpha_rcl2*fZ*(fZ*(Lrad-fCoulomb) + Lprad);
}
////////////////////////////////////////////////////////////////////////////////
/// Print this isotope

void TGeoElement::Print(Option_t *option) const
{
   printf("Element: %s      Z=%d   N=%f   A=%f [g/mole]\n", GetName(), fZ,Neff(),fA);
   if (HasIsotopes()) {
      for (Int_t i=0; i<fNisotopes; i++) {
         TGeoIsotope *iso = GetIsotope(i);
         printf("=>Isotope %s, abundance=%f :\n", iso->GetName(), fAbundances[i]);
         iso->Print(option);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns pointer to the table.

TGeoElementTable *TGeoElement::GetElementTable()
{
   if (!gGeoManager) {
      ::Error("TGeoElementTable::GetElementTable", "Create a geometry manager first");
      return NULL;
   }
   return gGeoManager->GetElementTable();
}

////////////////////////////////////////////////////////////////////////////////
/// Add an isotope for this element. All isotopes have to be isotopes of the same element.

void TGeoElement::AddIsotope(TGeoIsotope *isotope, Double_t relativeAbundance)
{
   if (!fIsotopes) {
      Fatal("AddIsotope", "Cannot add isotopes to normal elements. Use constructor with number of isotopes.");
      return;
   }
   Int_t ncurrent = 0;
   TGeoIsotope *isocrt;
   for (ncurrent=0; ncurrent<fNisotopes; ncurrent++)
      if (!fIsotopes->At(ncurrent)) break;
   if (ncurrent==fNisotopes) {
      Error("AddIsotope", "All %d isotopes of element %s already defined", fNisotopes, GetName());
      return;
   }
   // Check Z of the new isotope
   if ((fZ!=0) && (isotope->GetZ()!=fZ)) {
      Fatal("AddIsotope", "Trying to add isotope %s with different Z to the same element %s",
            isotope->GetName(), GetName());
      return;
   } else {
      fZ = isotope->GetZ();
   }
   fIsotopes->Add(isotope);
   fAbundances[ncurrent] = relativeAbundance;
   if (ncurrent==fNisotopes-1) {
      Double_t weight = 0.0;
      Double_t aeff = 0.0;
      Double_t neff = 0.0;
      for (Int_t i=0; i<fNisotopes; i++) {
         isocrt = (TGeoIsotope*)fIsotopes->At(i);
         aeff += fAbundances[i]*isocrt->GetA();
         neff += fAbundances[i]*isocrt->GetN();
         weight += fAbundances[i];
      }
      aeff /= weight;
      neff /= weight;
      fN = (Int_t)neff;
      fA = aeff;
   }
   ComputeDerivedQuantities();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns effective number of nucleons.

Double_t TGeoElement::Neff() const
{
   if (!fNisotopes) return fN;
   TGeoIsotope *isocrt;
   Double_t weight = 0.0;
   Double_t neff = 0.0;
   for (Int_t i=0; i<fNisotopes; i++) {
      isocrt = (TGeoIsotope*)fIsotopes->At(i);
      neff += fAbundances[i]*isocrt->GetN();
      weight += fAbundances[i];
   }
   neff /= weight;
   return neff;
}

////////////////////////////////////////////////////////////////////////////////
/// Return i-th isotope in the element.

TGeoIsotope *TGeoElement::GetIsotope(Int_t i) const
{
   if (i>=0 && i<fNisotopes) {
      return (TGeoIsotope*)fIsotopes->At(i);
   }
   return NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// Return relative abundance of i-th isotope in this element.

Double_t TGeoElement::GetRelativeAbundance(Int_t i) const
{
   if (i>=0 && i<fNisotopes) return fAbundances[i];
   return 0.0;
}

ClassImp(TGeoIsotope);

////////////////////////////////////////////////////////////////////////////////
/// Dummy I/O constructor

TGeoIsotope::TGeoIsotope()
            :TNamed(),
             fZ(0),
             fN(0),
             fA(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TGeoIsotope::TGeoIsotope(const char *name, Int_t z, Int_t n, Double_t a)
            :TNamed(name,""),
             fZ(z),
             fN(n),
             fA(a)
{
   if (z<1) Fatal("ctor", "Not allowed Z=%d (<1) for isotope: %s", z,name);
   if (n<z) Fatal("ctor", "Not allowed Z=%d < N=%d for isotope: %s", z,n,name);
   TGeoElement::GetElementTable()->AddIsotope(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Find existing isotope by name.

TGeoIsotope *TGeoIsotope::FindIsotope(const char *name)
{
   TGeoElementTable *elTable = TGeoElement::GetElementTable();
   if (!elTable) return 0;
   return elTable->FindIsotope(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Print this isotope

void TGeoIsotope::Print(Option_t *) const
{
   printf("Isotope: %s     Z=%d   N=%d   A=%f [g/mole]\n", GetName(), fZ,fN,fA);
}

ClassImp(TGeoElementRN);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoElementRN::TGeoElementRN()
{
   TObject::SetBit(kElementChecked,kFALSE);
   fENDFcode = 0;
   fIso      = 0;
   fLevel    = 0;
   fDeltaM   = 0;
   fHalfLife = 0;
   fNatAbun  = 0;
   fTH_F     = 0;
   fTG_F     = 0;
   fTH_S     = 0;
   fTG_S     = 0;
   fStatus   = 0;
   fRatio    = 0;
   fDecays   = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGeoElementRN::TGeoElementRN(Int_t aa, Int_t zz, Int_t iso, Double_t level,
               Double_t deltaM, Double_t halfLife, const char* JP,
               Double_t natAbun, Double_t th_f, Double_t tg_f, Double_t th_s,
               Double_t tg_s, Int_t status)
              :TGeoElement("", JP, zz, aa)
{
   TObject::SetBit(kElementChecked,kFALSE);
   fENDFcode = ENDF(aa,zz,iso);
   fIso      = iso;
   fLevel    = level;
   fDeltaM   = deltaM;
   fHalfLife = halfLife;
   fTitle    = JP;
   if (!fTitle.Length()) fTitle = "?";
   fNatAbun  = natAbun;
   fTH_F     = th_f;
   fTG_F     = tg_f;
   fTH_S     = th_s;
   fTG_S     = tg_s;
   fStatus   = status;
   fDecays   = 0;
   fRatio    = 0;
   MakeName(aa,zz,iso);
   if ((TMath::Abs(fHalfLife)<1.e-30) || fHalfLife<-1) Warning("ctor","Element %s has T1/2=%g [s]", fName.Data(), fHalfLife);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoElementRN::~TGeoElementRN()
{
   if (fDecays) {
      fDecays->Delete();
      delete fDecays;
   }
   if (fRatio) delete fRatio;
}

////////////////////////////////////////////////////////////////////////////////
/// Adds a decay mode for this element.

void TGeoElementRN::AddDecay(Int_t decay, Int_t diso, Double_t branchingRatio, Double_t qValue)
{
   if (branchingRatio<1E-20) {
      TString decayName;
      TGeoDecayChannel::DecayName(decay, decayName);
      Warning("AddDecay", "Decay %s of %s has BR=0. Not added.", decayName.Data(),fName.Data());
      return;
   }
   TGeoDecayChannel *dc = new TGeoDecayChannel(decay,diso,branchingRatio, qValue);
   dc->SetParent(this);
   if (!fDecays) fDecays = new TObjArray(5);
   fDecays->Add(dc);
}

////////////////////////////////////////////////////////////////////////////////
/// Adds a decay channel to the list of decays.

void TGeoElementRN::AddDecay(TGeoDecayChannel *dc)
{
   dc->SetParent(this);
   if (!fDecays) fDecays = new TObjArray(5);
   fDecays->Add(dc);
}

////////////////////////////////////////////////////////////////////////////////
/// Get number of decay channels of this element.

Int_t TGeoElementRN::GetNdecays() const
{
   if (!fDecays) return 0;
   return fDecays->GetEntriesFast();
}

////////////////////////////////////////////////////////////////////////////////
/// Get the activity in Bq of a gram of material made from this element.

Double_t TGeoElementRN::GetSpecificActivity() const
{
   static const Double_t ln2 = TMath::Log(2.);
   Double_t sa = (fHalfLife>0 && fA>0)?(ln2*TMath::Na()/fHalfLife/fA):0.;
   return sa;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if all decay chain of the element is OK.

Bool_t TGeoElementRN::CheckDecays() const
{
   if (TObject::TestBit(kElementChecked)) return kTRUE;
   TObject *oelem = (TObject*)this;
   TGeoDecayChannel *dc;
   TGeoElementRN *elem;
   TGeoElementTable *table = GetElementTable();
   TString decayName;
   if (!table) {
      Error("CheckDecays", "Element table not present");
      return kFALSE;
   }
   Bool_t resultOK = kTRUE;
   if (!fDecays) {
      oelem->SetBit(kElementChecked,kTRUE);
      return resultOK;
   }
   Double_t br = 0.;
   Int_t decayResult = 0;
   TIter next(fDecays);
   while ((dc=(TGeoDecayChannel*)next())) {
      br += dc->BranchingRatio();
      decayResult = DecayResult(dc);
      if (decayResult) {
         elem = table->GetElementRN(decayResult);
         if (!elem) {
            TGeoDecayChannel::DecayName(dc->Decay(),decayName);
            Error("CheckDecays", "Element after decay %s of %s not found in DB", decayName.Data(),fName.Data());
            return kFALSE;
         }
         dc->SetDaughter(elem);
         resultOK = elem->CheckDecays();
      }
   }
   if (TMath::Abs(br-100) > 1.E-3) {
      Warning("CheckDecays", "BR for decays of element %s sum-up = %f", fName.Data(), br);
      resultOK = kFALSE;
   }
   oelem->SetBit(kElementChecked,kTRUE);
   return resultOK;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns ENDF code of decay result.

Int_t TGeoElementRN::DecayResult(TGeoDecayChannel *dc) const
{
   Int_t da, dz, diso;
   dc->DecayShift(da, dz, diso);
   if (da == -99 || dz == -99) return 0;
   return ENDF(Int_t(fA)+da,fZ+dz,fIso+diso);
}

////////////////////////////////////////////////////////////////////////////////
/// Fills the input array with the set of RN elements resulting from the decay of
/// this one. All element in the list will contain the time evolution of their
/// proportion by number with respect to this element. The proportion can be
/// retrieved via the method TGeoElementRN::Ratio().
/// The precision represent the minimum cumulative branching ratio for
/// which decay products are still taken into account.

void TGeoElementRN::FillPopulation(TObjArray *population, Double_t precision, Double_t factor)
{
   TGeoElementRN *elem;
   TGeoElemIter next(this, precision);
   TGeoBatemanSol s(this);
   s.Normalize(factor);
   AddRatio(s);
   if (!population->FindObject(this)) population->Add(this);
   while ((elem=next())) {
      TGeoBatemanSol ratio(next.GetBranch());
      ratio.Normalize(factor);
      elem->AddRatio(ratio);
      if (!population->FindObject(elem)) population->Add(elem);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Generate a default name for the element.

void TGeoElementRN::MakeName(Int_t a, Int_t z, Int_t iso)
{
   fName = "";
   if (z==0 && a==1) {
      fName = "neutron";
      return;
   }
   if (z>=1 && z<= gMaxElem) fName += TString::Format("%3d-%s-",z,gElName[z-1]);
   else fName = "?? -?? -";
   if (a>=1 && a<=999) fName += TString::Format("%3.3d",a);
   else fName += "??";
   if (iso>0 && iso<gMaxLevel) fName += TString::Format("%c", gLevName[iso]);
   fName.ReplaceAll(" ","");
}

////////////////////////////////////////////////////////////////////////////////
/// Print info about the element;

void TGeoElementRN::Print(Option_t *option) const
{
   printf("\n%-12s ",fName.Data());
   printf("ENDF=%d; ",fENDFcode);
   printf("A=%d; ",(Int_t)fA);
   printf("Z=%d; ",fZ);
   printf("Iso=%d; ",fIso);
   printf("Level=%g[MeV]; ",fLevel);
   printf("Dmass=%g[MeV]; ",fDeltaM);
   if (fHalfLife>0) printf("Hlife=%g[s]\n",fHalfLife);
   else printf("Hlife=INF\n");
   printf("%13s"," ");
   printf("J/P=%s; ",fTitle.Data());
   printf("Abund=%g; ",fNatAbun);
   printf("Htox=%g; ",fTH_F);
   printf("Itox=%g; ",fTG_F);
   printf("Stat=%d\n",fStatus);
   if(!fDecays) return;
   printf("Decay modes:\n");
   TIter next(fDecays);
   TGeoDecayChannel *dc;
   while ((dc=(TGeoDecayChannel*)next())) dc->Print(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Create element from line record.

TGeoElementRN *TGeoElementRN::ReadElementRN(const char *line, Int_t &ndecays)
{
   Int_t a,z,iso,status;
   Double_t level, deltaM, halfLife, natAbun, th_f, tg_f, th_s, tg_s;
   char name[20],jp[20];
   sscanf(&line[0], "%s%d%d%d%lg%lg%lg%s%lg%lg%lg%lg%lg%d%d", name,&a,&z,&iso,&level,&deltaM,
          &halfLife,jp,&natAbun,&th_f,&tg_f,&th_s,&tg_s,&status,&ndecays);
   TGeoElementRN *elem = new TGeoElementRN(a,z,iso,level,deltaM,halfLife,
                                           jp,natAbun,th_f,tg_f,th_s,tg_s,status);
   return elem;
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive for RN elements.

void TGeoElementRN::SavePrimitive(std::ostream &out, Option_t *option)
{
   if (!strcmp(option,"h")) {
      // print a header if requested
      out << "#====================================================================================================================================" << std::endl;
      out << "#   Name      A    Z   ISO    LEV[MeV]  DM[MeV]   T1/2[s]        J/P     ABUND[%]    HTOX      ITOX      HTOX      ITOX    STAT NDCY" << std::endl;
      out << "#====================================================================================================================================" << std::endl;
   }
   out << std::setw(11) << fName.Data();
   out << std::setw(5) << (Int_t)fA;
   out << std::setw(5) << fZ;
   out << std::setw(5) << fIso;
   out << std::setw(10) << std::setiosflags(std::ios::fixed) << std::setprecision(5) << fLevel;
   out << std::setw(10) << std::setiosflags(std::ios::fixed) << std::setprecision(5) << fDeltaM;
   out << std::setw(10) << std::setiosflags(std::ios::scientific) << std::setprecision(3) << fHalfLife;
   out << std::setw(13) << fTitle.Data();
   out << std::setw(10) << std::setiosflags(std::ios::fixed) << std::setprecision(5) << fNatAbun;
   out << std::setw(10) << std::setiosflags(std::ios::fixed) << std::setprecision(5) << fTH_F;
   out << std::setw(10) << std::setiosflags(std::ios::fixed) << std::setprecision(5) << fTG_F;
   out << std::setw(10) << std::setiosflags(std::ios::fixed) << std::setprecision(5) << fTH_S;
   out << std::setw(10) << std::setiosflags(std::ios::fixed) << std::setprecision(5) << fTG_S;
   out << std::setw(5) << fStatus;
   Int_t ndecays = 0;
   if (fDecays) ndecays = fDecays->GetEntries();
   out << std::setw(5) << ndecays;
   out << std::endl;
   if (fDecays) {
      TIter next(fDecays);
      TGeoDecayChannel *dc;
      while ((dc=(TGeoDecayChannel*)next())) dc->SavePrimitive(out);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Adds a proportion ratio to the existing one.

void TGeoElementRN::AddRatio(TGeoBatemanSol &ratio)
{
   if (!fRatio) fRatio = new TGeoBatemanSol(ratio);
   else         *fRatio += ratio;
}

////////////////////////////////////////////////////////////////////////////////
/// Clears the existing ratio.

void TGeoElementRN::ResetRatio()
{
   if (fRatio) {
      delete fRatio;
      fRatio = 0;
   }
}

ClassImp(TGeoDecayChannel);

////////////////////////////////////////////////////////////////////////////////
/// Assignment.
///assignment operator

TGeoDecayChannel& TGeoDecayChannel::operator=(const TGeoDecayChannel& dc)
{
   if(this!=&dc) {
      TObject::operator=(dc);
      fDecay          = dc.fDecay;
      fDiso           = dc.fDiso;
      fBranchingRatio = dc.fBranchingRatio;
      fQvalue         = dc.fQvalue;
      fParent         = dc.fParent;
      fDaughter       = dc.fDaughter;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns name of decay.

const char *TGeoDecayChannel::GetName() const
{
   static TString name = "";
   name = "";
   if (!fDecay) return gDecayName[gMaxDecay];
   for (Int_t i=0; i<gMaxDecay; i++) {
      if (1<<i & fDecay) {
         if (name.Length()) name += "+";
         name += gDecayName[i];
      }
   }
   return name.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns decay name.

void TGeoDecayChannel::DecayName(UInt_t decay, TString &name)
{
   if (!decay) {
      name = gDecayName[gMaxDecay];
      return;
   }
   name = "";
   for (Int_t i=0; i<gMaxDecay; i++) {
      if (1<<i & decay) {
         if (name.Length()) name += "+";
         name += gDecayName[i];
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get index of this channel in the list of decays of the parent nuclide.

Int_t TGeoDecayChannel::GetIndex() const
{
   return fParent->Decays()->IndexOf(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Prints decay info.

void TGeoDecayChannel::Print(Option_t *) const
{
   TString name;
   DecayName(fDecay, name);
   printf("%-20s Diso: %3d BR: %9.3f%% Qval: %g\n", name.Data(),fDiso,fBranchingRatio,fQvalue);
}

////////////////////////////////////////////////////////////////////////////////
/// Create element from line record.

TGeoDecayChannel *TGeoDecayChannel::ReadDecay(const char *line)
{
   char name[80];
   Int_t decay,diso;
   Double_t branchingRatio, qValue;
   sscanf(&line[0], "%s%d%d%lg%lg", name,&decay,&diso,&branchingRatio,&qValue);
   TGeoDecayChannel *dc = new TGeoDecayChannel(decay,diso,branchingRatio,qValue);
   return dc;
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive for decays.

void TGeoDecayChannel::SavePrimitive(std::ostream &out, Option_t *)
{
   TString decayName;
   DecayName(fDecay, decayName);
   out << std::setw(50) << decayName.Data();
   out << std::setw(10) << fDecay;
   out << std::setw(10) << fDiso;
   out << std::setw(12) << std::setiosflags(std::ios::fixed) << std::setprecision(6) << fBranchingRatio;
   out << std::setw(12) << std::setiosflags(std::ios::fixed) << std::setprecision(6) << fQvalue;
   out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns variation in A, Z and Iso after decay.

void TGeoDecayChannel::DecayShift(Int_t &dA, Int_t &dZ, Int_t &dI) const
{
   dA=dZ=0;
   dI=fDiso;
   for(Int_t i=0; i<gMaxDecay; ++i) {
      if(1<<i & fDecay) {
         if(gDecayDeltaA[i] == -99 || gDecayDeltaZ[i] == -99 ) {
            dA=dZ=-99;
            return;
         }
         dA += gDecayDeltaA[i];
         dZ += gDecayDeltaZ[i];
      }
   }
}

ClassImp(TGeoElemIter);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TGeoElemIter::TGeoElemIter(TGeoElementRN *top, Double_t limit)
          : fTop(top), fElem(top), fBranch(0), fLevel(0), fLimitRatio(limit), fRatio(1.)
{
   fBranch = new TObjArray(10);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor.

TGeoElemIter::TGeoElemIter(const TGeoElemIter &iter)
             :fTop(iter.fTop),
              fElem(iter.fElem),
              fBranch(0),
              fLevel(iter.fLevel),
              fLimitRatio(iter.fLimitRatio),
              fRatio(iter.fRatio)
{
   if (iter.fBranch) {
      fBranch = new TObjArray(10);
      for (Int_t i=0; i<fLevel; i++) fBranch->Add(iter.fBranch->At(i));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoElemIter::~TGeoElemIter()
{
   if (fBranch) delete fBranch;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment.

TGeoElemIter &TGeoElemIter::operator=(const TGeoElemIter &iter)
{
   if (&iter == this) return *this;
   fTop   = iter.fTop;
   fElem  = iter.fElem;
   fLevel = iter.fLevel;
   if (iter.fBranch) {
      fBranch = new TObjArray(10);
      for (Int_t i=0; i<fLevel; i++) fBranch->Add(iter.fBranch->At(i));
   }
   fLimitRatio = iter.fLimitRatio;
   fRatio = iter.fRatio;
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// () operator.

TGeoElementRN *TGeoElemIter::operator()()
{
   return Next();
}

////////////////////////////////////////////////////////////////////////////////
/// Go upwards from the current location until the next branching, then down.

TGeoElementRN *TGeoElemIter::Up()
{
   TGeoDecayChannel *dc;
   Int_t ind, nd;
   while (fLevel) {
      // Current decay channel
      dc = (TGeoDecayChannel*)fBranch->At(fLevel-1);
      ind = dc->GetIndex();
      nd = dc->Parent()->GetNdecays();
      fRatio /= 0.01*dc->BranchingRatio();
      fElem = dc->Parent();
      fBranch->RemoveAt(--fLevel);
      ind++;
      while (ind<nd) {
         if (Down(ind++)) return (TGeoElementRN*)fElem;
      }
   }
   fElem = NULL;
   return NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// Go downwards from current level via ibranch as low in the tree as possible.
/// Return value flags if the operation was successful.

TGeoElementRN *TGeoElemIter::Down(Int_t ibranch)
{
   TGeoDecayChannel *dc = (TGeoDecayChannel*)fElem->Decays()->At(ibranch);
   if (!dc->Daughter()) return NULL;
   Double_t br = 0.01*fRatio*dc->BranchingRatio();
   if (br < fLimitRatio) return NULL;
   fLevel++;
   fRatio = br;
   fBranch->Add(dc);
   fElem = dc->Daughter();
   return (TGeoElementRN*)fElem;
}

////////////////////////////////////////////////////////////////////////////////
/// Return next element.

TGeoElementRN *TGeoElemIter::Next()
{
   if (!fElem) return NULL;
   // Check if this is the first iteration.
   Int_t nd = fElem->GetNdecays();
   for (Int_t i=0; i<nd; i++) if (Down(i)) return (TGeoElementRN*)fElem;
   return Up();
}

////////////////////////////////////////////////////////////////////////////////
/// Print info about the current decay branch.

void TGeoElemIter::Print(Option_t * /*option*/) const
{
   TGeoElementRN *elem;
   TGeoDecayChannel *dc;
   TString indent = "";
   printf("=== Chain with %g %%\n", 100*fRatio);
   for (Int_t i=0; i<fLevel; i++) {
      dc = (TGeoDecayChannel*)fBranch->At(i);
      elem = dc->Parent();
      printf("%s%s (%g%% %s) T1/2=%g\n", indent.Data(), elem->GetName(),dc->BranchingRatio(),dc->GetName(),elem->HalfLife());
      indent += " ";
      if (i==fLevel-1) {
         elem = dc->Daughter();
         printf("%s%s\n", indent.Data(), elem->GetName());
      }
   }
}

ClassImp(TGeoElementTable);

////////////////////////////////////////////////////////////////////////////////
/// default constructor

TGeoElementTable::TGeoElementTable()
{
   fNelements = 0;
   fNelementsRN = 0;
   fNisotopes = 0;
   fList      = 0;
   fListRN    = 0;
   fIsotopes = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoElementTable::TGeoElementTable(Int_t /*nelements*/)
{
   fNelements = 0;
   fNelementsRN = 0;
   fNisotopes = 0;
   fList = new TObjArray(128);
   fListRN    = 0;
   fIsotopes = 0;
   BuildDefaultElements();
//   BuildElementsRN();
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TGeoElementTable::TGeoElementTable(const TGeoElementTable& get) :
  TObject(get),
  fNelements(get.fNelements),
  fNelementsRN(get.fNelementsRN),
  fNisotopes(get.fNisotopes),
  fList(get.fList),
  fListRN(get.fListRN),
  fIsotopes(0)
{
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoElementTable& TGeoElementTable::operator=(const TGeoElementTable& get)
{
   if(this!=&get) {
      TObject::operator=(get);
      fNelements=get.fNelements;
      fNelementsRN=get.fNelementsRN;
      fNisotopes=get.fNisotopes;
      fList=get.fList;
      fListRN=get.fListRN;
      fIsotopes = 0;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoElementTable::~TGeoElementTable()
{
   if (fList) {
      fList->Delete();
      delete fList;
   }
   if (fListRN) {
      fListRN->Delete();
      delete fListRN;
   }
   if (fIsotopes) {
      fIsotopes->Delete();
      delete fIsotopes;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add an element to the table. Obsolete.

void TGeoElementTable::AddElement(const char *name, const char *title, Int_t z, Double_t a)
{
   if (!fList) fList = new TObjArray(128);
   fList->AddAtAndExpand(new TGeoElement(name,title,z,a), fNelements++);
}

////////////////////////////////////////////////////////////////////////////////
/// Add an element to the table.

void TGeoElementTable::AddElement(const char *name, const char *title, Int_t z, Int_t n, Double_t a)
{
   if (!fList) fList = new TObjArray(128);
   fList->AddAtAndExpand(new TGeoElement(name,title,z,n,a), fNelements++);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a custom element to the table.

void TGeoElementTable::AddElement(TGeoElement *elem)
{
   if (!fList) fList = new TObjArray(128);
   TGeoElement *orig = FindElement(elem->GetName());
   if (orig) {
      Error("AddElement", "Found element with same name: %s (%s). Cannot add to table.",
             orig->GetName(), orig->GetTitle());
      return;
   }
   fList->AddAtAndExpand(elem, fNelements++);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a radionuclide to the table and map it.

void TGeoElementTable::AddElementRN(TGeoElementRN *elem)
{
   if (!fListRN) fListRN = new TObjArray(3600);
   if (HasRNElements() && GetElementRN(elem->ENDFCode())) return;
//   elem->Print();
   fListRN->Add(elem);
   fNelementsRN++;
   fElementsRN.insert(ElementRNMap_t::value_type(elem->ENDFCode(), elem));
}

////////////////////////////////////////////////////////////////////////////////
/// Add isotope to the table.

void TGeoElementTable::AddIsotope(TGeoIsotope *isotope)
{
   if (FindIsotope(isotope->GetName())) {
      Error("AddIsotope", "Isotope with the same name: %s already in table. Not adding.",isotope->GetName());
      return;
   }
   if (!fIsotopes) fIsotopes = new TObjArray();
   fIsotopes->Add(isotope);
}

////////////////////////////////////////////////////////////////////////////////
/// Creates the default element table

void TGeoElementTable::BuildDefaultElements()
{
   if (HasDefaultElements()) return;
   AddElement("VACUUM","VACUUM"   ,0,   0, 0.0);
   AddElement("H"   ,"HYDROGEN"   ,1,   1, 1.00794);
   AddElement("HE"  ,"HELIUM"     ,2,   4, 4.002602);
   AddElement("LI"  ,"LITHIUM"    ,3,   7, 6.941);
   AddElement("BE"  ,"BERYLLIUM"  ,4,   9, 9.01218);
   AddElement("B"   ,"BORON"      ,5,  11, 10.811);
   AddElement("C"   ,"CARBON"     ,6,  12, 12.0107);
   AddElement("N"   ,"NITROGEN"   ,7,  14, 14.00674);
   AddElement("O"   ,"OXYGEN"     ,8,  16, 15.9994);
   AddElement("F"   ,"FLUORINE"   ,9,  19, 18.9984032);
   AddElement("NE"  ,"NEON"       ,10, 20, 20.1797);
   AddElement("NA"  ,"SODIUM"     ,11, 23, 22.989770);
   AddElement("MG"  ,"MAGNESIUM"  ,12, 24, 24.3050);
   AddElement("AL"  ,"ALUMINIUM"  ,13, 27, 26.981538);
   AddElement("SI"  ,"SILICON"    ,14, 28, 28.0855);
   AddElement("P"   ,"PHOSPHORUS" ,15, 31, 30.973761);
   AddElement("S"   ,"SULFUR"     ,16, 32, 32.066);
   AddElement("CL"  ,"CHLORINE"   ,17, 35, 35.4527);
   AddElement("AR"  ,"ARGON"      ,18, 40, 39.948);
   AddElement("K"   ,"POTASSIUM"  ,19, 39, 39.0983);
   AddElement("CA"  ,"CALCIUM"    ,20, 40, 40.078);
   AddElement("SC"  ,"SCANDIUM"   ,21, 45, 44.955910);
   AddElement("TI"  ,"TITANIUM"   ,22, 48, 47.867);
   AddElement("V"   ,"VANADIUM"   ,23, 51, 50.9415);
   AddElement("CR"  ,"CHROMIUM"   ,24, 52, 51.9961);
   AddElement("MN"  ,"MANGANESE"  ,25, 55, 54.938049);
   AddElement("FE"  ,"IRON"       ,26, 56, 55.845);
   AddElement("CO"  ,"COBALT"     ,27, 59, 58.933200);
   AddElement("NI"  ,"NICKEL"     ,28, 59, 58.6934);
   AddElement("CU"  ,"COPPER"     ,29, 64, 63.546);
   AddElement("ZN"  ,"ZINC"       ,30, 65, 65.39);
   AddElement("GA"  ,"GALLIUM"    ,31, 70, 69.723);
   AddElement("GE"  ,"GERMANIUM"  ,32, 73, 72.61);
   AddElement("AS"  ,"ARSENIC"    ,33, 75, 74.92160);
   AddElement("SE"  ,"SELENIUM"   ,34, 79, 78.96);
   AddElement("BR"  ,"BROMINE"    ,35, 80, 79.904);
   AddElement("KR"  ,"KRYPTON"    ,36, 84, 83.80);
   AddElement("RB"  ,"RUBIDIUM"   ,37, 85, 85.4678);
   AddElement("SR"  ,"STRONTIUM"  ,38, 88, 87.62);
   AddElement("Y"   ,"YTTRIUM"    ,39, 89, 88.90585);
   AddElement("ZR"  ,"ZIRCONIUM"  ,40, 91, 91.224);
   AddElement("NB"  ,"NIOBIUM"    ,41, 93, 92.90638);
   AddElement("MO"  ,"MOLYBDENUM" ,42, 96, 95.94);
   AddElement("TC"  ,"TECHNETIUM" ,43, 98, 98.0);
   AddElement("RU"  ,"RUTHENIUM"  ,44, 101, 101.07);
   AddElement("RH"  ,"RHODIUM"    ,45, 103, 102.90550);
   AddElement("PD"  ,"PALLADIUM"  ,46, 106, 106.42);
   AddElement("AG"  ,"SILVER"     ,47, 108, 107.8682);
   AddElement("CD"  ,"CADMIUM"    ,48, 112, 112.411);
   AddElement("IN"  ,"INDIUM"     ,49, 115, 114.818);
   AddElement("SN"  ,"TIN"        ,50, 119, 118.710);
   AddElement("SB"  ,"ANTIMONY"   ,51, 122, 121.760);
   AddElement("TE"  ,"TELLURIUM"  ,52, 128, 127.60);
   AddElement("I"   ,"IODINE"     ,53, 127, 126.90447);
   AddElement("XE"  ,"XENON"      ,54, 131, 131.29);
   AddElement("CS"  ,"CESIUM"     ,55, 133, 132.90545);
   AddElement("BA"  ,"BARIUM"     ,56, 137, 137.327);
   AddElement("LA"  ,"LANTHANUM"  ,57, 139, 138.9055);
   AddElement("CE"  ,"CERIUM"     ,58, 140, 140.116);
   AddElement("PR"  ,"PRASEODYMIUM" ,59, 141, 140.90765);
   AddElement("ND"  ,"NEODYMIUM"  ,60, 144, 144.24);
   AddElement("PM"  ,"PROMETHIUM" ,61, 145, 145.0);
   AddElement("SM"  ,"SAMARIUM"   ,62, 150, 150.36);
   AddElement("EU"  ,"EUROPIUM"   ,63, 152, 151.964);
   AddElement("GD"  ,"GADOLINIUM" ,64, 157, 157.25);
   AddElement("TB"  ,"TERBIUM"    ,65, 159, 158.92534);
   AddElement("DY"  ,"DYSPROSIUM" ,66, 162, 162.50);
   AddElement("HO"  ,"HOLMIUM"    ,67, 165, 164.93032);
   AddElement("ER"  ,"ERBIUM"     ,68, 167, 167.26);
   AddElement("TM"  ,"THULIUM"    ,69, 169, 168.93421);
   AddElement("YB"  ,"YTTERBIUM"  ,70, 173, 173.04);
   AddElement("LU"  ,"LUTETIUM"   ,71, 175, 174.967);
   AddElement("HF"  ,"HAFNIUM"    ,72, 178, 178.49);
   AddElement("TA"  ,"TANTALUM"   ,73, 181, 180.9479);
   AddElement("W"   ,"TUNGSTEN"   ,74, 184, 183.84);
   AddElement("RE"  ,"RHENIUM"    ,75, 186, 186.207);
   AddElement("OS"  ,"OSMIUM"     ,76, 190, 190.23);
   AddElement("IR"  ,"IRIDIUM"    ,77, 192, 192.217);
   AddElement("PT"  ,"PLATINUM"   ,78, 195, 195.078);
   AddElement("AU"  ,"GOLD"       ,79, 197, 196.96655);
   AddElement("HG"  ,"MERCURY"    ,80, 200, 200.59);
   AddElement("TL"  ,"THALLIUM"   ,81, 204, 204.3833);
   AddElement("PB"  ,"LEAD"       ,82, 207, 207.2);
   AddElement("BI"  ,"BISMUTH"    ,83, 209, 208.98038);
   AddElement("PO"  ,"POLONIUM"   ,84, 209, 209.0);
   AddElement("AT"  ,"ASTATINE"   ,85, 210, 210.0);
   AddElement("RN"  ,"RADON"      ,86, 222, 222.0);
   AddElement("FR"  ,"FRANCIUM"   ,87, 223, 223.0);
   AddElement("RA"  ,"RADIUM"     ,88, 226, 226.0);
   AddElement("AC"  ,"ACTINIUM"   ,89, 227, 227.0);
   AddElement("TH"  ,"THORIUM"    ,90, 232, 232.0381);
   AddElement("PA"  ,"PROTACTINIUM" ,91, 231, 231.03588);
   AddElement("U"   ,"URANIUM"    ,92, 238, 238.0289);
   AddElement("NP"  ,"NEPTUNIUM"  ,93, 237, 237.0);
   AddElement("PU"  ,"PLUTONIUM"  ,94, 244, 244.0);
   AddElement("AM"  ,"AMERICIUM"  ,95, 243, 243.0);
   AddElement("CM"  ,"CURIUM"     ,96, 247, 247.0);
   AddElement("BK"  ,"BERKELIUM"  ,97, 247, 247.0);
   AddElement("CF"  ,"CALIFORNIUM",98, 251, 251.0);
   AddElement("ES"  ,"EINSTEINIUM",99, 252, 252.0);
   AddElement("FM"  ,"FERMIUM"    ,100, 257, 257.0);
   AddElement("MD"  ,"MENDELEVIUM",101, 258, 258.0);
   AddElement("NO"  ,"NOBELIUM"   ,102, 259, 259.0);
   AddElement("LR"  ,"LAWRENCIUM" ,103, 262, 262.0);
   AddElement("RF"  ,"RUTHERFORDIUM",104, 261, 261.0);
   AddElement("DB"  ,"DUBNIUM" ,105, 262, 262.0);
   AddElement("SG"  ,"SEABORGIUM" ,106, 263, 263.0);
   AddElement("BH"  ,"BOHRIUM"    ,107, 262, 262.0);
   AddElement("HS"  ,"HASSIUM"    ,108, 265, 265.0);
   AddElement("MT"  ,"MEITNERIUM" ,109, 266, 266.0);
   AddElement("UUN" ,"UNUNNILIUM" ,110, 269, 269.0);
   AddElement("UUU" ,"UNUNUNIUM"  ,111, 272, 272.0);
   AddElement("UUB" ,"UNUNBIUM"   ,112, 277, 277.0);

   TObject::SetBit(kETDefaultElements,kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Creates the list of radionuclides.

void TGeoElementTable::ImportElementsRN()
{
   if (HasRNElements()) return;
   TGeoElementRN *elem;
   TString rnf = "RadioNuclides.txt";
   gSystem->PrependPathName(TROOT::GetEtcDir(), rnf);
   FILE *fp = fopen(rnf, "r");
   if (!fp) {
      Error("ImportElementsRN","File RadioNuclides.txt not found");
      return;
   }
   char line[150];
   Int_t ndecays = 0;
   Int_t i;
   while (fgets(&line[0],140,fp)) {
      if (line[0]=='#') continue;
      elem = TGeoElementRN::ReadElementRN(line, ndecays);
      for (i=0; i<ndecays; i++) {
         if (!fgets(&line[0],140,fp)) {
            Error("ImportElementsRN", "Error parsing RadioNuclides.txt file");
            fclose(fp);
            return;
         }
         TGeoDecayChannel *dc = TGeoDecayChannel::ReadDecay(line);
         elem->AddDecay(dc);
      }
      AddElementRN(elem);
//      elem->Print();
   }
   TObject::SetBit(kETRNElements,kTRUE);
   CheckTable();
   fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////
/// Checks status of element table.

Bool_t TGeoElementTable::CheckTable() const
{
   if (!HasRNElements()) return HasDefaultElements();
   TGeoElementRN *elem;
   Bool_t result = kTRUE;
   TIter next(fListRN);
   while ((elem=(TGeoElementRN*)next())) {
      if (!elem->CheckDecays()) result = kFALSE;
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Export radionuclides in a file.

void TGeoElementTable::ExportElementsRN(const char *filename)
{
   if (!HasRNElements()) return;
   TString sname = filename;
   if (!sname.Length()) sname = "RadioNuclides.txt";
   std::ofstream out;
   out.open(sname.Data(), std::ios::out);
   if (!out.good()) {
      Error("ExportElementsRN", "Cannot open file %s", sname.Data());
      return;
   }

   TGeoElementRN *elem;
   TIter next(fListRN);
   Int_t i=0;
   while ((elem=(TGeoElementRN*)next())) {
      if ((i%48)==0) elem->SavePrimitive(out,"h");
      else elem->SavePrimitive(out);
      i++;
   }
   out.close();
}

////////////////////////////////////////////////////////////////////////////////
/// Search an element by symbol or full name
/// Exact matching

TGeoElement *TGeoElementTable::FindElement(const char *name) const
{
   TGeoElement *elem;
   elem = (TGeoElement*)fList->FindObject(name);
   if (elem) return elem;
   // Search case insensitive by element name
   TString s(name);
   s.ToUpper();
   elem = (TGeoElement*)fList->FindObject(s.Data());
   if (elem) return elem;
   // Search by full name
   TIter next(fList);
   while ((elem=(TGeoElement*)next())) {
      if (s == elem->GetTitle()) return elem;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Find existing isotope by name. Not optimized for a big number of isotopes.

TGeoIsotope *TGeoElementTable::FindIsotope(const char *name) const
{
   if (!fIsotopes) return NULL;
   return (TGeoIsotope*)fIsotopes->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve a radionuclide by ENDF code.

TGeoElementRN *TGeoElementTable::GetElementRN(Int_t ENDFcode) const
{
   if (!HasRNElements()) {
      TGeoElementTable *table = (TGeoElementTable*)this;
      table->ImportElementsRN();
      if (!fListRN) return 0;
   }
   ElementRNMap_t::const_iterator it = fElementsRN.find(ENDFcode);
   if (it != fElementsRN.end()) return it->second;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve a radionuclide by a, z, and isomeric state.

TGeoElementRN *TGeoElementTable::GetElementRN(Int_t a, Int_t z, Int_t iso) const
{
   return GetElementRN(TGeoElementRN::ENDF(a,z,iso));
}

////////////////////////////////////////////////////////////////////////////////
/// Print table of elements. The accepted options are:
///  ""  - prints everything by default
///  "D" - prints default elements only
///  "I" - prints isotopes
///  "R" - prints radio-nuclides only if imported
///  "U" - prints user-defined elements only

void TGeoElementTable::Print(Option_t *option) const
{
   TString opt(option);
   opt.ToUpper();
   Int_t induser = HasDefaultElements() ? 113 : 0;
   // Default elements
   if (opt=="" || opt=="D") {
      if (induser) printf("================\nDefault elements\n================\n");
      for (Int_t iel=0; iel<induser; ++iel) fList->At(iel)->Print();
   }
   // Isotopes
   if (opt=="" || opt=="I") {
      if (fIsotopes) {
         printf("================\nIsotopes\n================\n");
         fIsotopes->Print();
      }
   }
   // Radio-nuclides
   if (opt=="" || opt=="R") {
      if (HasRNElements()) {
         printf("================\nRadio-nuclides\n================\n");
         fListRN->Print();
      }
   }
   // User-defined elements
   if (opt=="" || opt=="U") {
      if (fNelements>induser) printf("================\nUser elements\n================\n");
      for (Int_t iel=induser; iel<fNelements; ++iel) fList->At(iel)->Print();
   }
}

ClassImp(TGeoBatemanSol);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TGeoBatemanSol::TGeoBatemanSol(TGeoElementRN *elem)
               :TObject(), TAttLine(), TAttFill(), TAttMarker(),
                fElem(elem),
                fElemTop(elem),
                fCsize(10),
                fNcoeff(0),
                fFactor(1.),
                fTmin(0.),
                fTmax(0.),
                fCoeff(NULL)
{
   fCoeff = new BtCoef_t[fCsize];
   fNcoeff = 1;
   fCoeff[0].cn = 1.;
   Double_t t12 = elem->HalfLife();
   if (t12 == 0.) t12 = 1.e-30;
   if (elem->Stable()) fCoeff[0].lambda = 0.;
   else                fCoeff[0].lambda = TMath::Log(2.)/t12;
}

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TGeoBatemanSol::TGeoBatemanSol(const TObjArray *chain)
               :TObject(), TAttLine(), TAttFill(), TAttMarker(),
                fElem(NULL),
                fElemTop(NULL),
                fCsize(0),
                fNcoeff(0),
                fFactor(1.),
                fTmin(0.),
                fTmax(0.),
                fCoeff(NULL)
{
   TGeoDecayChannel *dc = (TGeoDecayChannel*)chain->At(0);
   if (dc) fElemTop = dc->Parent();
   dc = (TGeoDecayChannel*)chain->At(chain->GetEntriesFast()-1);
   if (dc) {
      fElem = dc->Daughter();
      fCsize = chain->GetEntriesFast()+1;
      fCoeff = new BtCoef_t[fCsize];
      FindSolution(chain);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TGeoBatemanSol::TGeoBatemanSol(const TGeoBatemanSol& other)
               :TObject(other), TAttLine(other), TAttFill(other), TAttMarker(other),
                fElem(other.fElem),
                fElemTop(other.fElemTop),
                fCsize(other.fCsize),
                fNcoeff(other.fNcoeff),
                fFactor(other.fFactor),
                fTmin(other.fTmin),
                fTmax(other.fTmax),
                fCoeff(NULL)
{
   if (fCsize) {
      fCoeff = new BtCoef_t[fCsize];
      for (Int_t i=0; i<fNcoeff; i++) {
         fCoeff[i].cn = other.fCoeff[i].cn;
         fCoeff[i].lambda = other.fCoeff[i].lambda;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoBatemanSol::~TGeoBatemanSol()
{
   if (fCoeff) delete [] fCoeff;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment.

TGeoBatemanSol& TGeoBatemanSol::operator=(const TGeoBatemanSol& other)
{
   if (this == &other) return *this;
   TObject::operator=(other);
   TAttLine::operator=(other);
   TAttFill::operator=(other);
   TAttMarker::operator=(other);
   fElem = other.fElem;
   fElemTop = other.fElemTop;
   if (fCoeff) delete [] fCoeff;
   fCoeff = 0;
   fCsize = other.fCsize;
   fNcoeff = other.fNcoeff;
   fFactor = other.fFactor;
   fTmin = other.fTmin;
   fTmax = other.fTmax;
   if (fCsize) {
      fCoeff = new BtCoef_t[fCsize];
      for (Int_t i=0; i<fNcoeff; i++) {
         fCoeff[i].cn = other.fCoeff[i].cn;
         fCoeff[i].lambda = other.fCoeff[i].lambda;
      }
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Addition of other solution.

TGeoBatemanSol& TGeoBatemanSol::operator+=(const TGeoBatemanSol& other)
{
   if (other.GetElement() != fElem) {
      Error("operator+=", "Cannot add 2 solutions for different elements");
      return *this;
   }
   Int_t i,j;
   BtCoef_t *coeff = fCoeff;
   Int_t ncoeff = fNcoeff + other.fNcoeff;
   if (ncoeff > fCsize) {
      fCsize = ncoeff;
      coeff = new BtCoef_t[ncoeff];
      for (i=0; i<fNcoeff; i++) {
         coeff[i].cn = fCoeff[i].cn;
         coeff[i].lambda = fCoeff[i].lambda;
      }
      delete [] fCoeff;
      fCoeff = coeff;
   }
   ncoeff = fNcoeff;
   for (j=0; j<other.fNcoeff; j++) {
      for (i=0; i<fNcoeff; i++) {
         if (coeff[i].lambda == other.fCoeff[j].lambda) {
            coeff[i].cn += fFactor * other.fCoeff[j].cn;
            break;
         }
      }
      if (i == fNcoeff) {
         coeff[ncoeff].cn = fFactor * other.fCoeff[j].cn;
         coeff[ncoeff].lambda = other.fCoeff[j].lambda;
         ncoeff++;
      }
   }
   fNcoeff = ncoeff;
   return *this;
}
////////////////////////////////////////////////////////////////////////////////
/// Find concentration of the element at a given time.

Double_t TGeoBatemanSol::Concentration(Double_t time) const
{
   Double_t conc = 0.;
   for (Int_t i=0; i<fNcoeff; i++)
      conc += fCoeff[i].cn * TMath::Exp(-fCoeff[i].lambda * time);
   return conc;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the solution of Bateman equation versus time.

void TGeoBatemanSol::Draw(Option_t *option)
{
   if (!gGeoManager) return;
   gGeoManager->GetGeomPainter()->DrawBatemanSol(this, option);
}

////////////////////////////////////////////////////////////////////////////////
/// Find the solution for the Bateman equations corresponding to the decay
/// chain described by an array ending with element X.
/// A->B->...->X
/// Cn = SUM [Ain * exp(-LMBDi*t)];
///      Cn    - concentration Nx/Na
///      n     - order of X in chain (A->B->X => n=3)
///      LMBDi - decay constant for element of order i in the chain
///      Ain = LMBD1*...*LMBD(n-1) * br1*...*br(n-1)/(LMBD1-LMBDi)...(LMBDn-LMBDi)
///      bri   - branching ratio for decay Ei->Ei+1

void TGeoBatemanSol::FindSolution(const TObjArray *array)
{
   fNcoeff = 0;
   if (!array || !array->GetEntriesFast()) return;
   Int_t n = array->GetEntriesFast();
   TGeoDecayChannel *dc = (TGeoDecayChannel*)array->At(n-1);
   TGeoElementRN *elem = dc->Daughter();
   if (elem != fElem) {
      Error("FindSolution", "Last element in the list must be %s\n", fElem->GetName());
      return;
   }
   Int_t i,j;
   Int_t order = n+1;
   if (!fCoeff) {
      fCsize = order;
      fCoeff = new BtCoef_t[fCsize];
   }
   if (fCsize < order) {
      delete [] fCoeff;
      fCsize = order;
      fCoeff = new BtCoef_t[fCsize];
   }

   Double_t *lambda = new Double_t[order];
   Double_t *br     = new Double_t[n];
   Double_t halflife;
   for (i=0; i<n; i++) {
      dc = (TGeoDecayChannel*)array->At(i);
      elem = dc->Parent();
      br[i] = 0.01 * dc->BranchingRatio();
      halflife = elem->HalfLife();
      if (halflife==0.) halflife = 1.e-30;
      if (elem->Stable()) lambda[i] = 0.;
      else                lambda[i] = TMath::Log(2.)/halflife;
      if (i==n-1) {
         elem = dc->Daughter();
         halflife = elem->HalfLife();
         if (halflife==0.) halflife = 1.e-30;
         if (elem->Stable()) lambda[n] = 0.;
         else                lambda[n] = TMath::Log(2.)/halflife;
      }
   }
   // Check if we have equal lambdas
   for (i=0; i<order-1; i++) {
      for (j=i+1; j<order; j++) {
         if (lambda[j] == lambda[i]) lambda[j] += 0.001*lambda[j];
      }
   }
   Double_t ain;
   Double_t pdlambda, plambdabr=1.;
   for (j=0; j<n; j++) plambdabr *= lambda[j]*br[j];
   for (i=0; i<order; i++) {
      pdlambda = 1.;
      for (j=0; j<n+1; j++) {
         if (j == i) continue;
         pdlambda *= lambda[j] - lambda[i];
      }
      if (pdlambda == 0.) {
         Error("FindSolution", "pdlambda=0 !!!");
         delete [] lambda;
         delete [] br;
         return;
      }
      ain = plambdabr/pdlambda;
      fCoeff[i].cn = ain;
      fCoeff[i].lambda = lambda[i];
   }
   fNcoeff = order;
   Normalize(fFactor);
   delete [] lambda;
   delete [] br;
}

////////////////////////////////////////////////////////////////////////////////
/// Normalize all coefficients with a given factor.

void TGeoBatemanSol::Normalize(Double_t factor)
{
   for (Int_t i=0; i<fNcoeff; i++) fCoeff[i].cn *= factor;
}

////////////////////////////////////////////////////////////////////////////////
/// Print concentration evolution.

void TGeoBatemanSol::Print(Option_t * /*option*/) const
{
   TString formula;
   formula.Form("N[%s]/N[%s] = ", fElem->GetName(), fElemTop->GetName());
   for (Int_t i=0; i<fNcoeff; i++) {
      if (i == fNcoeff-1) formula += TString::Format("%g*exp(-%g*t)", fCoeff[i].cn, fCoeff[i].lambda);
      else                formula += TString::Format("%g*exp(-%g*t) + ", fCoeff[i].cn, fCoeff[i].lambda);
   }
   printf("%s\n", formula.Data());
}

