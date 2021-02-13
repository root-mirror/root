// @(#)root/geom:$Id$
// Author: Andrei Gheata   25/10/01

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoMaterial
\ingroup Geometry_classes

Base class describing materials.

\image html geom_material.jpg
*/

#include <iostream>
#include "TMath.h"
#include "TObjArray.h"
#include "TGeoElement.h"
#include "TGeoManager.h"
#include "TGeoExtension.h"
#include "TGeoMaterial.h"
#include "TGeoPhysicalConstants.h"
#include "TGeant4PhysicalConstants.h"
#include "TGDMLMatrix.h"

// statics and globals

ClassImp(TGeoMaterial);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoMaterial::TGeoMaterial()
             :TNamed(), TAttFill(),
              fIndex(0),
              fA(0.),
              fZ(0.),
              fDensity(0.),
              fRadLen(0.),
              fIntLen(0.),
              fTemperature(0.),
              fPressure(0.),
              fState(kMatStateUndefined),
              fShader(NULL),
              fCerenkov(NULL),
              fElement(NULL),
              fUserExtension(0),
              fFWExtension(0)
{
   TGeoManager::SetDefaultUnits(TGeoManager::GetDefaultUnits()); // Ensure nobody changes the units afterwards
   SetUsed(kFALSE);
   fIndex    = -1;
   fTemperature = STP_temperature;
   fPressure = STP_pressure;
   fState = kMatStateUndefined;
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoMaterial::TGeoMaterial(const char *name)
             :TNamed(name, ""), TAttFill(),
              fIndex(0),
              fA(0.),
              fZ(0.),
              fDensity(0.),
              fRadLen(0.),
              fIntLen(0.),
              fTemperature(0.),
              fPressure(0.),
              fState(kMatStateUndefined),
              fShader(NULL),
              fCerenkov(NULL),
              fElement(NULL),
              fUserExtension(0),
              fFWExtension(0)
{
   TGeoManager::SetDefaultUnits(TGeoManager::GetDefaultUnits()); // Ensure nobody changes the units afterwards
   fName = fName.Strip();
   SetUsed(kFALSE);
   fIndex    = -1;
   fTemperature = STP_temperature;
   fPressure = STP_pressure;
   fState = kMatStateUndefined;

   if (!gGeoManager) {
      gGeoManager = new TGeoManager("Geometry", "default geometry");
   }
   gGeoManager->AddMaterial(this);
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoMaterial::TGeoMaterial(const char *name, Double_t a, Double_t z,
                Double_t rho, Double_t radlen, Double_t intlen)
             :TNamed(name, ""), TAttFill(),
              fIndex(0),
              fA(a),
              fZ(z),
              fDensity(rho),
              fRadLen(0.),
              fIntLen(0.),
              fTemperature(0.),
              fPressure(0.),
              fState(kMatStateUndefined),
              fShader(NULL),
              fCerenkov(NULL),
              fElement(NULL),
              fUserExtension(0),
              fFWExtension(0)
{
   TGeoManager::SetDefaultUnits(TGeoManager::GetDefaultUnits()); // Ensure nobody changes the units afterwards
   fName = fName.Strip();
   SetUsed(kFALSE);
   fIndex    = -1;
   fA        = a;
   fZ        = z;
   fDensity  = rho;
   fTemperature = STP_temperature;
   fPressure = STP_pressure;
   fState = kMatStateUndefined;
   SetRadLen(radlen, intlen);
   if (!gGeoManager) {
      gGeoManager = new TGeoManager("Geometry", "default geometry");
   }
   if (fZ - Int_t(fZ) > 1E-3)
      Warning("ctor", "Material %s defined with fractional Z=%f", GetName(), fZ);
   if (GetElement()) GetElement()->SetUsed();
   gGeoManager->AddMaterial(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with state, temperature and pressure.

TGeoMaterial::TGeoMaterial(const char *name, Double_t a, Double_t z, Double_t rho,
                EGeoMaterialState state, Double_t temperature, Double_t pressure)
             :TNamed(name, ""), TAttFill(),
              fIndex(0),
              fA(a),
              fZ(z),
              fDensity(rho),
              fRadLen(0.),
              fIntLen(0.),
              fTemperature(temperature),
              fPressure(pressure),
              fState(state),
              fShader(NULL),
              fCerenkov(NULL),
              fElement(NULL),
              fUserExtension(0),
              fFWExtension(0)
{
   TGeoManager::SetDefaultUnits(TGeoManager::GetDefaultUnits()); // Ensure nobody changes the units afterwards
   fName = fName.Strip();
   SetUsed(kFALSE);
   fIndex    = -1;
   SetRadLen(0,0);
   if (!gGeoManager) {
      gGeoManager = new TGeoManager("Geometry", "default geometry");
   }
   if (fZ - Int_t(fZ) > 1E-3)
      Warning("ctor", "Material %s defined with fractional Z=%f", GetName(), fZ);
   if (GetElement()) GetElement()->SetUsed();
   gGeoManager->AddMaterial(this);
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoMaterial::TGeoMaterial(const char *name, TGeoElement *elem, Double_t rho)
             :TNamed(name, ""), TAttFill(),
              fIndex(0),
              fA(0.),
              fZ(0.),
              fDensity(rho),
              fRadLen(0.),
              fIntLen(0.),
              fTemperature(0.),
              fPressure(0.),
              fState(kMatStateUndefined),
              fShader(NULL),
              fCerenkov(NULL),
              fElement(elem),
              fUserExtension(0),
              fFWExtension(0)
{
   TGeoManager::SetDefaultUnits(TGeoManager::GetDefaultUnits()); // Ensure nobody changes the units afterwards
   fName = fName.Strip();
   SetUsed(kFALSE);
   fIndex    = -1;
   fA        = elem->A();
   fZ        = elem->Z();
   SetRadLen(0,0);
   fTemperature = STP_temperature;
   fPressure = STP_pressure;
   fState = kMatStateUndefined;
   if (!gGeoManager) {
      gGeoManager = new TGeoManager("Geometry", "default geometry");
   }
   if (fZ - Int_t(fZ) > 1E-3)
      Warning("ctor", "Material %s defined with fractional Z=%f", GetName(), fZ);
   if (GetElement()) GetElement()->SetUsed();
   gGeoManager->AddMaterial(this);
}

////////////////////////////////////////////////////////////////////////////////

TGeoMaterial::TGeoMaterial(const TGeoMaterial& gm) :
              TNamed(gm),
              TAttFill(gm),
              fIndex(gm.fIndex),
              fA(gm.fA),
              fZ(gm.fZ),
              fDensity(gm.fDensity),
              fRadLen(gm.fRadLen),
              fIntLen(gm.fIntLen),
              fTemperature(gm.fTemperature),
              fPressure(gm.fPressure),
              fState(gm.fState),
              fShader(gm.fShader),
              fCerenkov(gm.fCerenkov),
              fElement(gm.fElement),
              fUserExtension(gm.fUserExtension->Grab()),
              fFWExtension(gm.fFWExtension->Grab())

{
   //copy constructor
   TGeoManager::SetDefaultUnits(TGeoManager::GetDefaultUnits()); // Ensure nobody changes the units afterwards
   fProperties.SetOwner();
   TIter next(&fProperties);
   TNamed *property;
   while ((property = (TNamed*)next())) fProperties.Add(new TNamed(*property));
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TGeoMaterial& TGeoMaterial::operator=(const TGeoMaterial& gm)
{
   if(this!=&gm) {
      TNamed::operator=(gm);
      TAttFill::operator=(gm);
      fIndex=gm.fIndex;
      fA=gm.fA;
      fZ=gm.fZ;
      fDensity=gm.fDensity;
      fRadLen=gm.fRadLen;
      fIntLen=gm.fIntLen;
      fTemperature=gm.fTemperature;
      fPressure=gm.fPressure;
      fState=gm.fState;
      fShader=gm.fShader;
      fCerenkov=gm.fCerenkov;
      fElement=gm.fElement;
      fUserExtension = gm.fUserExtension->Grab();
      fFWExtension = gm.fFWExtension->Grab();
      fProperties.SetOwner();
      TIter next(&fProperties);
      TNamed *property;
      while ((property = (TNamed*)next())) fProperties.Add(new TNamed(*property));
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoMaterial::~TGeoMaterial()
{
   if (fUserExtension) {fUserExtension->Release(); fUserExtension=0;}
   if (fFWExtension) {fFWExtension->Release(); fFWExtension=0;}
}

////////////////////////////////////////////////////////////////////////////////
/// Connect user-defined extension to the material. The material "grabs" a copy, so
/// the original object can be released by the producer. Release the previously
/// connected extension if any.
///
/// NOTE: This interface is intended for user extensions and is guaranteed not
/// to be used by TGeo

void TGeoMaterial::SetUserExtension(TGeoExtension *ext)
{
   if (fUserExtension) fUserExtension->Release();
   fUserExtension = 0;
   if (ext) fUserExtension = ext->Grab();
}

//_____________________________________________________________________________
const char *TGeoMaterial::GetPropertyRef(const char *property) const
{
   // Find reference for a given property
   TNamed *prop = (TNamed*)fProperties.FindObject(property);
   return (prop) ? prop->GetTitle() : nullptr;
}

//_____________________________________________________________________________
TGDMLMatrix *TGeoMaterial::GetProperty(const char *property) const
{
   // Find reference for a given property
   TNamed *prop = (TNamed*)fProperties.FindObject(property);
   if ( !prop ) return nullptr;
   return gGeoManager->GetGDMLMatrix(prop->GetTitle());
}

//_____________________________________________________________________________
TGDMLMatrix *TGeoMaterial::GetProperty(Int_t i) const
{
   // Find reference for a given property
   TNamed *prop = (TNamed*)fProperties.At(i);
   if ( !prop ) return nullptr;
   return gGeoManager->GetGDMLMatrix(prop->GetTitle());
}

//_____________________________________________________________________________
const char *TGeoMaterial::GetConstPropertyRef(const char *property) const
{
   // Find reference for a given constant property
   TNamed *prop = (TNamed*)fConstProperties.FindObject(property);
   return (prop) ? prop->GetTitle() : nullptr;
}

//_____________________________________________________________________________
Double_t TGeoMaterial::GetConstProperty(const char *property, Bool_t *err) const
{
   // Find reference for a given constant property
   TNamed *prop = (TNamed*)fConstProperties.FindObject(property);
   if (!prop) {
      if (err) *err = kTRUE;
      return 0.;
   }
   return gGeoManager->GetProperty(prop->GetTitle(), err);
}

//_____________________________________________________________________________
Double_t TGeoMaterial::GetConstProperty(Int_t i, Bool_t *err) const
{
   // Find reference for a given constant property
   TNamed *prop = (TNamed*)fConstProperties.At(i);
   if (!prop) {
      if (err) *err = kTRUE;
      return 0.;
   }
   return gGeoManager->GetProperty(prop->GetTitle(), err);
}

//_____________________________________________________________________________
bool TGeoMaterial::AddProperty(const char *property, const char *ref)
{
   fProperties.SetOwner();
   if (GetPropertyRef(property)) {
      Error("AddProperty", "Property %s already added to material %s",
         property, GetName());
      return false;
   }
   fProperties.Add(new TNamed(property, ref));
   return true;
}

//_____________________________________________________________________________
bool TGeoMaterial::AddConstProperty(const char *property, const char *ref)
{
   fConstProperties.SetOwner();
   if (GetConstPropertyRef(property)) {
      Error("AddConstProperty", "Constant property %s already added to material %s",
         property, GetName());
      return false;
   }
   fConstProperties.Add(new TNamed(property, ref));
   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Connect framework defined extension to the material. The material "grabs" a copy,
/// so the original object can be released by the producer. Release the previously
/// connected extension if any.
///
/// NOTE: This interface is intended for the use by TGeo and the users should
///       NOT connect extensions using this method

void TGeoMaterial::SetFWExtension(TGeoExtension *ext)
{
   if (fFWExtension) fFWExtension->Release();
   fFWExtension = 0;
   if (ext) fFWExtension = ext->Grab();
}

////////////////////////////////////////////////////////////////////////////////
/// Get a copy of the user extension pointer. The user must call Release() on
/// the copy pointer once this pointer is not needed anymore (equivalent to
/// delete() after calling new())

TGeoExtension *TGeoMaterial::GrabUserExtension() const
{
   if (fUserExtension) return fUserExtension->Grab();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get a copy of the framework extension pointer. The user must call Release() on
/// the copy pointer once this pointer is not needed anymore (equivalent to
/// delete() after calling new())

TGeoExtension *TGeoMaterial::GrabFWExtension() const
{
   if (fFWExtension) return fFWExtension->Grab();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Provide a pointer name containing uid.

char *TGeoMaterial::GetPointerName() const
{
   static TString name;
   name = TString::Format("pMat%d", GetUniqueID());
   return (char*)name.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Set radiation/absorption lengths. If the values are negative, their absolute value
/// is taken, otherwise radlen is recomputed using G3 formula.

void TGeoMaterial::SetRadLen(Double_t radlen, Double_t intlen)
{
   fRadLen = TMath::Abs(radlen);
   fIntLen = TMath::Abs(intlen);
   // Check for vacuum
   if (fA<0.9 || fZ<0.9) {
      if (radlen<-1e5 || intlen<-1e-5) {
         Error("SetRadLen","Material %s: user values taken for vacuum: radlen=%g or intlen=%g - too small", GetName(),fRadLen, fIntLen);
         return;
      }
      // Ignore positive values and take big numbers
      if (radlen>=0) fRadLen = 1.E30;
      if (intlen>=0) fIntLen = 1.E30;
      return;
   }
   TGeoManager::EDefaultUnits typ = TGeoManager::GetDefaultUnits();
   // compute radlen systematically with G3 formula for a valid material
   if ( typ == TGeoManager::kRootUnits && radlen>=0 ) {
      //taken grom Geant3 routine GSMATE
      constexpr Double_t alr2av = 1.39621E-03*TGeoUnit::cm2;
      constexpr Double_t al183  = 5.20948;
      fRadLen = fA/(alr2av*fDensity*fZ*(fZ +TGeoMaterial::ScreenFactor(fZ))*
                   (al183-TMath::Log(fZ)/3-TGeoMaterial::Coulomb(fZ)));
      fRadLen *= TGeoUnit::cm;
   }
   else if ( typ == TGeoManager::kG4Units && radlen>=0 ) {
      //taken grom Geant3 routine GSMATE
      constexpr Double_t alr2av = 1.39621E-03*TGeant4Unit::cm2;
      constexpr Double_t al183  = 5.20948;
      fRadLen = fA/(alr2av*fDensity*fZ*(fZ +TGeoMaterial::ScreenFactor(fZ))*
                   (al183-TMath::Log(fZ)/3-TGeoMaterial::Coulomb(fZ)));
      fRadLen *= TGeant4Unit::cm;
   }
   // Compute interaction length using the same formula as in GEANT4
   if ( typ == TGeoManager::kRootUnits && intlen>=0 ) {
      constexpr Double_t lambda0 = 35.*TGeoUnit::g/TGeoUnit::cm2;  // [g/cm^2]
      Double_t nilinv = 0.0;
      TGeoElement *elem = GetElement();
      if (!elem) {
         Fatal("SetRadLen", "Element not found for material %s", GetName());
         return;
      }
      Double_t nbAtomsPerVolume = TGeoUnit::Avogadro*fDensity/elem->A();
      nilinv += nbAtomsPerVolume*TMath::Power(elem->Neff(), 0.6666667);
      nilinv *= TGeoUnit::amu/lambda0;
      fIntLen = (nilinv<=0) ? TGeoShape::Big() : (TGeoUnit::cm/nilinv);
   }
   else if ( typ == TGeoManager::kG4Units && intlen>=0 ) {
      constexpr Double_t lambda0 = 35.*TGeant4Unit::g/TGeant4Unit::cm2;  // [g/cm^2]
      Double_t nilinv = 0.0;
      TGeoElement *elem = GetElement();
      if (!elem) {
         Fatal("SetRadLen", "Element not found for material %s", GetName());
         return;
      }
      Double_t nbAtomsPerVolume = TGeant4Unit::Avogadro*fDensity/elem->A();
      nilinv += nbAtomsPerVolume*TMath::Power(elem->Neff(), 0.6666667);
      nilinv *= TGeant4Unit::amu/lambda0;
      fIntLen = (nilinv<=0) ? TGeoShape::Big() : (TGeant4Unit::cm/nilinv);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// static function
///  Compute Coulomb correction for pair production and Brem
///  REFERENCE : EGS MANUAL SLAC 210 - UC32 - JUNE 78
///                        FORMULA 2.7.17

Double_t TGeoMaterial::Coulomb(Double_t z)
{
   Double_t az    = TGeoManager::kRootUnits == TGeoManager::GetDefaultUnits()
     ? TGeoUnit::fine_structure_const*z : TGeant4Unit::fine_structure_const*z;
   Double_t az2   = az*az;
   Double_t az4   = az2 * az2;
   Double_t fp    = ( 0.0083*az4 + 0.20206 + 1./(1.+az2) ) * az2;
   Double_t fm    = ( 0.0020*az4 + 0.0369  ) * az4;
   return fp - fm;
}

////////////////////////////////////////////////////////////////////////////////
/// return true if the other material has the same physical properties

Bool_t TGeoMaterial::IsEq(const TGeoMaterial *other) const
{
   if (other==this) return kTRUE;
   if (other->IsMixture()) return kFALSE;
   if (TMath::Abs(fA-other->GetA())>1E-3) return kFALSE;
   if (TMath::Abs(fZ-other->GetZ())>1E-3) return kFALSE;
   if (TMath::Abs(fDensity-other->GetDensity())>1E-6) return kFALSE;
   if (GetCerenkovProperties() != other->GetCerenkovProperties()) return kFALSE;
//   if (fRadLen != other->GetRadLen()) return kFALSE;
//   if (fIntLen != other->GetIntLen()) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// print characteristics of this material

void TGeoMaterial::Print(const Option_t * /*option*/) const
{
   printf("Material %s %s   A=%g Z=%g rho=%g radlen=%g intlen=%g index=%i\n", GetName(), GetTitle(),
          fA,fZ,fDensity, fRadLen, fIntLen, fIndex);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoMaterial::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TestBit(TGeoMaterial::kMatSavePrimitive)) return;
   char *name = GetPointerName();
   out << "// Material: " << GetName() << std::endl;
   out << "   a       = " << fA << ";" << std::endl;
   out << "   z       = " << fZ << ";" << std::endl;
   out << "   density = " << fDensity << ";" << std::endl;
   out << "   radl    = " << fRadLen << ";" << std::endl;
   out << "   absl    = " << fIntLen << ";" << std::endl;

   out << "   " << name << " = new TGeoMaterial(\"" << GetName() << "\", a,z,density,radl,absl);" << std::endl;
   out << "   " << name << "->SetIndex(" << GetIndex() << ");" << std::endl;
   SetBit(TGeoMaterial::kMatSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Get some default color related to this material.

Int_t TGeoMaterial::GetDefaultColor() const
{
   Int_t id = 1+ gGeoManager->GetListOfMaterials()->IndexOf(this);
   return (2+id%6);
}

////////////////////////////////////////////////////////////////////////////////
/// Get a pointer to the element this material is made of.
/// This second call is to avaoid warnings to not call a virtual
/// method from the constructor

TGeoElement *TGeoMaterial::GetElement() const
{
   if (fElement) return fElement;
   TGeoElementTable *table = gGeoManager->GetElementTable();
   return table->GetElement(Int_t(fZ));
}

////////////////////////////////////////////////////////////////////////////////
/// Get a pointer to the element this material is made of.

TGeoElement *TGeoMaterial::GetElement(Int_t) const
{
   if (fElement) return fElement;
   TGeoElementTable *table = gGeoManager->GetElementTable();
   return table->GetElement(Int_t(fZ));
}

////////////////////////////////////////////////////////////////////////////////
/// Single interface to get element properties.

void TGeoMaterial::GetElementProp(Double_t &a, Double_t &z, Double_t &w, Int_t)
{
   a = fA;
   z = fZ;
   w = 1.;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve material index in the list of materials

Int_t TGeoMaterial::GetIndex()
{
   if (fIndex>=0) return fIndex;
   TList *matlist = gGeoManager->GetListOfMaterials();
   fIndex = matlist->IndexOf(this);
   return fIndex;
}

////////////////////////////////////////////////////////////////////////////////
/// Create the material representing the decay product of this material at a
/// given time. The precision represent the minimum cumulative branching ratio for
/// which decay products are still taken into account.

TGeoMaterial *TGeoMaterial::DecayMaterial(Double_t time, Double_t precision)
{
   TObjArray *pop = new TObjArray();
   if (!fElement || !fElement->IsRadioNuclide()) return this;
   FillMaterialEvolution(pop, precision);
   Int_t ncomp = pop->GetEntriesFast();
   if (!ncomp) return this;
   TGeoElementRN *el;
   Double_t *weight = new Double_t[ncomp];
   Double_t amed = 0.;
   Int_t i;
   for (i=0; i<ncomp; i++) {
      el = (TGeoElementRN *)pop->At(i);
      weight[i] = el->Ratio()->Concentration(time) * el->A();
      amed += weight[i];
   }
   Double_t rho = fDensity*amed/fA;
   TGeoMixture *mix = 0;
   Int_t ncomp1 = ncomp;
   for (i=0; i<ncomp; i++) {
      if ((weight[i]/amed)<precision) {
         amed -= weight[i];
         ncomp1--;
      }
   }
   if (ncomp1<2) {
      el = (TGeoElementRN *)pop->At(0);
      delete [] weight;
      delete pop;
      if (ncomp1==1) return new TGeoMaterial(TString::Format("%s-evol",GetName()), el, rho);
      return NULL;
   }
   mix = new TGeoMixture(TString::Format("%s-evol",GetName()), ncomp, rho);
   for (i=0; i<ncomp; i++) {
      weight[i] /= amed;
      if (weight[i]<precision) continue;
      el = (TGeoElementRN *)pop->At(i);
      mix->AddElement(el, weight[i]);
   }
   delete [] weight;
   delete pop;
   return mix;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a user array with all the elements deriving from the possible
/// decay of the top element composing the mixture. Each element contained
/// by <population> may be a radionuclide having a Bateman solution attached.
/// The precision represent the minimum cumulative branching ratio for
/// which decay products are still taken into account.
/// To visualize the time evolution of each decay product one can use:
/// ~~~ {.cpp}
///    TGeoElement *elem = population->At(index);
///    TGeoElementRN *elemrn = 0;
///    if (elem->IsRadioNuclide()) elemrn = (TGeoElementRN*)elem;
/// ~~~
/// One can get Ni/N1(t=0) at any moment of time. Ni is the number of atoms
/// of one of the decay products, N1(0) is the number of atoms of the top
/// element at t=0.
/// ~~~ {.cpp}
///    Double_t fraction_weight = elemrn->Ratio()->Concentration(time);
/// ~~~
/// One can also display the time evolution of the fractional weight:
/// ~~~ {.cpp}
///    elemrn->Ratio()->Draw(option);
/// ~~~

void TGeoMaterial::FillMaterialEvolution(TObjArray *population, Double_t precision)
{
   if (population->GetEntriesFast()) {
      Error("FillMaterialEvolution", "Provide an empty array !");
      return;
   }
   TGeoElementTable *table = gGeoManager->GetElementTable();
   TGeoElement *elem;
   TGeoElementRN *elemrn;
   TIter next(table->GetElementsRN());
   while ((elemrn=(TGeoElementRN*)next())) elemrn->ResetRatio();
   elem = GetElement();
   if (!elem) {
      Fatal("FillMaterialEvolution", "Element not found for material %s", GetName());
      return;
   }
   if (!elem->IsRadioNuclide()) {
      population->Add(elem);
      return;
   }
   elemrn = (TGeoElementRN*)elem;
   elemrn->FillPopulation(population, precision);
}

/** \class TGeoMixture
\ingroup Geometry_classes

Mixtures of elements.

*/

ClassImp(TGeoMixture);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoMixture::TGeoMixture()
{
   fNelements = 0;
   fZmixture  = 0;
   fAmixture  = 0;
   fWeights   = 0;
   fNatoms    = 0;
   fVecNbOfAtomsPerVolume = 0;
   fElements  = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoMixture::TGeoMixture(const char *name, Int_t /*nel*/, Double_t rho)
            :TGeoMaterial(name)
{
   fZmixture   = 0;
   fAmixture   = 0;
   fWeights    = 0;
   fNelements  = 0;
   fNatoms     = 0;
   fVecNbOfAtomsPerVolume = 0;
   fDensity = rho;
   fElements   = 0;
   if (fDensity < 0) fDensity = 0.001;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoMixture::~TGeoMixture()
{
   if (fZmixture) delete[] fZmixture;
   if (fAmixture) delete[] fAmixture;
   if (fWeights)  delete[] fWeights;
   if (fNatoms)   delete[] fNatoms;
   if (fVecNbOfAtomsPerVolume) delete[] fVecNbOfAtomsPerVolume;
   if (fElements) delete fElements;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute effective A/Z and radiation length

void TGeoMixture::AverageProperties()
{
   TGeoManager::EDefaultUnits typ = TGeoManager::GetDefaultUnits();
   const Double_t cm   = (typ==TGeoManager::kRootUnits) ? TGeoUnit::cm   : TGeant4Unit::cm;
   const Double_t cm2  = (typ==TGeoManager::kRootUnits) ? TGeoUnit::cm2  : TGeant4Unit::cm2;
   const Double_t amu  = (typ==TGeoManager::kRootUnits) ? TGeoUnit::amu  : TGeant4Unit::amu; // [MeV/c^2]
   const Double_t gram = (typ==TGeoManager::kRootUnits) ? TGeoUnit::gram : TGeant4Unit::gram;
   const Double_t na   = (typ==TGeoManager::kRootUnits) ? TGeoUnit::Avogadro : TGeant4Unit::Avogadro;
   const Double_t alr2av  = 1.39621E-03 * cm2;
   const Double_t al183   = 5.20948;
   const Double_t lambda0 = 35.*gram/cm2;  // [g/cm^2]
   Double_t radinv = 0.0;
   Double_t nilinv = 0.0;
   Double_t nbAtomsPerVolume;
   fA = 0;
   fZ = 0;
   for (Int_t j=0;j<fNelements;j++) {
      if (fWeights[j] <= 0) continue;
      fA += fWeights[j]*fAmixture[j];
      fZ += fWeights[j]*fZmixture[j];
      nbAtomsPerVolume = na*fDensity*fWeights[j]/GetElement(j)->A();
      nilinv += nbAtomsPerVolume*TMath::Power(GetElement(j)->Neff(), 0.6666667);
      Double_t zc = fZmixture[j];
      Double_t alz = TMath::Log(zc)/3.;
      Double_t xinv = zc*(zc+TGeoMaterial::ScreenFactor(zc))*
         (al183-alz-TGeoMaterial::Coulomb(zc))/fAmixture[j];
      radinv += xinv*fWeights[j];
   }
   radinv *= alr2av*fDensity;
   if (radinv > 0) fRadLen = cm/radinv;
   // Compute interaction length
   nilinv *= amu/lambda0;
   fIntLen = (nilinv<=0) ? TGeoShape::Big() : (cm/nilinv);
}

////////////////////////////////////////////////////////////////////////////////
/// add an element to the mixture using fraction by weight
/// Check if the element is already defined

void TGeoMixture::AddElement(Double_t a, Double_t z, Double_t weight)
{
   TGeoElementTable *table = gGeoManager->GetElementTable();
   if (z<1 || z>table->GetNelements()-1)
      Fatal("AddElement", "Cannot add element having Z=%d to mixture %s", (Int_t)z, GetName());
   Int_t i;
   for (i=0; i<fNelements; i++) {
      if (TMath::Abs(z-fZmixture[i])<1.e-6  && TMath::Abs(a-fAmixture[i])<1.e-6) {
         fWeights[i] += weight;
         AverageProperties();
         return;
      }
   }
   if (!fNelements) {
      fZmixture = new Double_t[1];
      fAmixture = new Double_t[1];
      fWeights  = new Double_t[1];
   } else {
      Int_t nelements = fNelements+1;
      Double_t *zmixture = new Double_t[nelements];
      Double_t *amixture = new Double_t[nelements];
      Double_t *weights  = new Double_t[nelements];
      for (Int_t j=0; j<fNelements; j++) {
         zmixture[j] = fZmixture[j];
         amixture[j] = fAmixture[j];
         weights[j]  = fWeights[j];
      }
      delete [] fZmixture;
      delete [] fAmixture;
      delete [] fWeights;
      fZmixture = zmixture;
      fAmixture = amixture;
      fWeights  = weights;
   }

   fNelements++;
   i = fNelements - 1;
   fZmixture[i] = z;
   fAmixture[i] = a;
   fWeights[i]  = weight;
   if (z - Int_t(z) > 1E-3)
      Warning("DefineElement", "Mixture %s has element defined with fractional Z=%f", GetName(), z);
   GetElement(i)->SetDefined();
   table->GetElement((Int_t)z)->SetDefined();

   //compute equivalent radiation length (taken from Geant3/GSMIXT)
   AverageProperties();
}

////////////////////////////////////////////////////////////////////////////////
/// Define one component of the mixture as an existing material/mixture.

void TGeoMixture::AddElement(TGeoMaterial *mat, Double_t weight)
{
   TGeoElement *elnew, *elem;
   Double_t a,z;
   if (!mat->IsMixture()) {
      elem = mat->GetBaseElement();
      if (elem) {
         AddElement(elem, weight);
      } else {
         a = mat->GetA();
         z = mat->GetZ();
         AddElement(a, z, weight);
      }
      return;
   }
   // The material is a mixture.
   TGeoMixture *mix = (TGeoMixture*)mat;
   Double_t wnew;
   Int_t nelem = mix->GetNelements();
   Bool_t elfound;
   Int_t i,j;
   // loop the elements of the daughter mixture
   for (i=0; i<nelem; i++) {
      elfound = kFALSE;
      elnew = mix->GetElement(i);
      if (!elnew) continue;
      // check if we have the element already defined in the parent mixture
      for (j=0; j<fNelements; j++) {
         if (fWeights[j]<=0) continue;
         elem = GetElement(j);
         if (elem == elnew) {
            // element found, compute new weight
            fWeights[j] += weight * (mix->GetWmixt())[i];
            elfound = kTRUE;
            break;
         }
      }
      if (elfound) continue;
      // element not found, define it
      wnew = weight * (mix->GetWmixt())[i];
      AddElement(elnew, wnew);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// add an element to the mixture using fraction by weight

void TGeoMixture::AddElement(TGeoElement *elem, Double_t weight)
{
   TGeoElement *elemold;
   TGeoElementTable *table = gGeoManager->GetElementTable();
   if (!fElements) fElements = new TObjArray(128);
   Bool_t exist = kFALSE;
   // If previous elements were defined by A/Z, add corresponding TGeoElements
   for (Int_t i=0; i<fNelements; i++) {
      elemold = (TGeoElement*)fElements->At(i);
      if (!elemold) fElements->AddAt(elemold = table->GetElement((Int_t)fZmixture[i]), i);
      if (elemold == elem) exist = kTRUE;
   }
   if (!exist) fElements->AddAtAndExpand(elem, fNelements);
   AddElement(elem->A(), elem->Z(), weight);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a mixture element by number of atoms in the chemical formula.

void TGeoMixture::AddElement(TGeoElement *elem, Int_t natoms)
{
   Int_t i,j;
   Double_t amol;
   TGeoElement *elemold;
   TGeoElementTable *table = gGeoManager->GetElementTable();
   if (!fElements) fElements = new TObjArray(128);
   // Check if the element is already defined
   for (i=0; i<fNelements; i++) {
      elemold = (TGeoElement*)fElements->At(i);
      if (!elemold) fElements->AddAt(table->GetElement((Int_t)fZmixture[i]), i);
      else if (elemold != elem) continue;
      if ((elem==elemold) ||
          (TMath::Abs(elem->Z()-fZmixture[i])<1.e-6 && TMath::Abs(elem->A()-fAmixture[i])<1.e-6)) {
         fNatoms[i] += natoms;
         amol = 0.;
         for (j=0; j<fNelements; j++) amol += fAmixture[j]*fNatoms[j];
         for (j=0; j<fNelements; j++) fWeights[j] = fNatoms[j]*fAmixture[j]/amol;
         AverageProperties();
         return;
      }
   }
   // New element
   if (!fNelements) {
      fZmixture = new Double_t[1];
      fAmixture = new Double_t[1];
      fWeights  = new Double_t[1];
      fNatoms   = new Int_t[1];
   } else {
      if (!fNatoms) {
         Fatal("AddElement", "Cannot add element by natoms in mixture %s after defining elements by weight",
               GetName());
         return;
      }
      Int_t nelements = fNelements+1;
      Double_t *zmixture = new Double_t[nelements];
      Double_t *amixture = new Double_t[nelements];
      Double_t *weights  = new Double_t[nelements];
      Int_t *nnatoms  = new Int_t[nelements];
      for (j=0; j<fNelements; j++) {
         zmixture[j] = fZmixture[j];
         amixture[j] = fAmixture[j];
         weights[j]  = fWeights[j];
         nnatoms[j]  = fNatoms[j];
      }
      delete [] fZmixture;
      delete [] fAmixture;
      delete [] fWeights;
      delete [] fNatoms;
      fZmixture = zmixture;
      fAmixture = amixture;
      fWeights  = weights;
      fNatoms   = nnatoms;
   }
   fNelements++;
   Int_t iel = fNelements-1;
   fZmixture[iel] = elem->Z();
   fAmixture[iel] = elem->A();
   fNatoms[iel]  = natoms;
   fElements->AddAtAndExpand(elem, iel);
   amol = 0.;
   for (i=0; i<fNelements; i++) {
      if (fNatoms[i]<=0) return;
      amol += fAmixture[i]*fNatoms[i];
   }
   for (i=0; i<fNelements; i++) fWeights[i] = fNatoms[i]*fAmixture[i]/amol;
   table->GetElement(elem->Z())->SetDefined();
   AverageProperties();
}

////////////////////////////////////////////////////////////////////////////////
/// Define the mixture element at index iel by number of atoms in the chemical formula.

void TGeoMixture::DefineElement(Int_t /*iel*/, Int_t z, Int_t natoms)
{
   TGeoElementTable *table = gGeoManager->GetElementTable();
   TGeoElement *elem = table->GetElement(z);
   if (!elem) {
      Fatal("DefineElement", "In mixture %s, element with Z=%i not found",GetName(),z);
      return;
   }
   AddElement(elem, natoms);
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve the pointer to the element corresponding to component I.

TGeoElement *TGeoMixture::GetElement(Int_t i) const
{
   if (i<0 || i>=fNelements) {
      Error("GetElement", "Mixture %s has only %d elements", GetName(), fNelements);
      return 0;
   }
   TGeoElement *elem = 0;
   if (fElements) elem = (TGeoElement*)fElements->At(i);
   if (elem) return elem;
   TGeoElementTable *table = gGeoManager->GetElementTable();
   return table->GetElement(Int_t(fZmixture[i]));
}

////////////////////////////////////////////////////////////////////////////////
/// Get specific activity (in Bq/gram) for the whole mixture (no argument) or
/// for a given component.

Double_t TGeoMixture::GetSpecificActivity(Int_t i) const
{
   if (i>=0 && i<fNelements) return fWeights[i]*GetElement(i)->GetSpecificActivity();
   Double_t sa = 0;
   for (Int_t iel=0; iel<fNelements; iel++) {
      sa += fWeights[iel]*GetElement(iel)->GetSpecificActivity();
   }
   return sa;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the other material has the same physical properties

Bool_t TGeoMixture::IsEq(const TGeoMaterial *other) const
{
   if (other->IsEqual(this)) return kTRUE;
   if (!other->IsMixture()) return kFALSE;
   TGeoMixture *mix = (TGeoMixture*)other;
   if (!mix) return kFALSE;
   if (fNelements != mix->GetNelements()) return kFALSE;
   if (TMath::Abs(fA-other->GetA())>1E-3) return kFALSE;
   if (TMath::Abs(fZ-other->GetZ())>1E-3) return kFALSE;
   if (TMath::Abs(fDensity-other->GetDensity())>1E-6) return kFALSE;
   if (GetCerenkovProperties() != other->GetCerenkovProperties()) return kFALSE;
//   if (fRadLen != other->GetRadLen()) return kFALSE;
//   if (fIntLen != other->GetIntLen()) return kFALSE;
   for (Int_t i=0; i<fNelements; i++) {
      if (TMath::Abs(fZmixture[i]-(mix->GetZmixt())[i])>1E-3) return kFALSE;
      if (TMath::Abs(fAmixture[i]-(mix->GetAmixt())[i])>1E-3) return kFALSE;
      if (TMath::Abs(fWeights[i]-(mix->GetWmixt())[i])>1E-3) return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// print characteristics of this material

void TGeoMixture::Print(const Option_t * /*option*/) const
{
   printf("Mixture %s %s   Aeff=%g Zeff=%g rho=%g radlen=%g intlen=%g index=%i\n", GetName(), GetTitle(),
          fA,fZ,fDensity, fRadLen, fIntLen, fIndex);
   for (Int_t i=0; i<fNelements; i++) {
      if (fNatoms) printf("   Element #%i : %s  Z=%6.2f A=%6.2f w=%6.3f natoms=%d\n", i, GetElement(i)->GetName(),fZmixture[i],
             fAmixture[i], fWeights[i], fNatoms[i]);
      else printf("   Element #%i : %s  Z=%6.2f A=%6.2f w=%6.3f\n", i, GetElement(i)->GetName(),fZmixture[i],
             fAmixture[i], fWeights[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoMixture::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TestBit(TGeoMaterial::kMatSavePrimitive)) return;
   char *name = GetPointerName();
   out << "// Mixture: " << GetName() << std::endl;
   out << "   nel     = " << fNelements << ";" << std::endl;
   out << "   density = " << fDensity << ";" << std::endl;
   out << "   " << name << " = new TGeoMixture(\"" << GetName() << "\", nel,density);" << std::endl;
   for (Int_t i=0; i<fNelements; i++) {
      TGeoElement *el = GetElement(i);
      out << "      a = " << fAmixture[i] << ";   z = "<< fZmixture[i] << ";   w = " << fWeights[i] << ";  // " << el->GetName() << std::endl;
      out << "   " << name << "->DefineElement(" << i << ",a,z,w);" << std::endl;
   }
   out << "   " << name << "->SetIndex(" << GetIndex() << ");" << std::endl;
   SetBit(TGeoMaterial::kMatSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Create the mixture representing the decay product of this material at a
/// given time. The precision represent the minimum cumulative branching ratio for
/// which decay products are still taken into account.

TGeoMaterial *TGeoMixture::DecayMaterial(Double_t time, Double_t precision)
{
   TObjArray *pop = new TObjArray();
   FillMaterialEvolution(pop, precision);
   Int_t ncomp = pop->GetEntriesFast();
   if (!ncomp) return this;
   TGeoElement *elem;
   TGeoElementRN *el;
   Double_t *weight = new Double_t[ncomp];
   Double_t amed = 0.;
   Int_t i, j;
   for (i=0; i<ncomp; i++) {
      elem = (TGeoElement *)pop->At(i);
      if (!elem->IsRadioNuclide()) {
         j = fElements->IndexOf(elem);
         weight[i] = fWeights[j]*fAmixture[0]/fWeights[0];
      } else {
         el = (TGeoElementRN*)elem;
         weight[i] = el->Ratio()->Concentration(time) * el->A();
      }
      amed += weight[i];
   }
   Double_t rho = fDensity * fWeights[0] * amed/fAmixture[0];
   TGeoMixture *mix = 0;
   Int_t ncomp1 = ncomp;
   for (i=0; i<ncomp; i++) {
      if ((weight[i]/amed)<precision) {
         amed -= weight[i];
         ncomp1--;
      }
   }
   if (ncomp1<2) {
      el = (TGeoElementRN *)pop->At(0);
      delete [] weight;
      delete pop;
      if (ncomp1==1) return new TGeoMaterial(TString::Format("%s-evol",GetName()), el, rho);
      return NULL;
   }
   mix = new TGeoMixture(TString::Format("%s-evol",GetName()), ncomp, rho);
   for (i=0; i<ncomp; i++) {
      weight[i] /= amed;
      if (weight[i]<precision) continue;
      el = (TGeoElementRN *)pop->At(i);
      mix->AddElement(el, weight[i]);
   }
   delete [] weight;
   delete pop;
   return mix;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a user array with all the elements deriving from the possible
/// decay of the top elements composing the mixture. Each element contained
/// by <population> may be a radionuclide having a Bateman solution attached.
/// The precision represent the minimum cumulative branching ratio for
/// which decay products are still taken into account.
/// To visualize the time evolution of each decay product one can use:
/// ~~~ {.cpp}
///    TGeoElement *elem = population->At(index);
///    TGeoElementRN *elemrn = 0;
///    if (elem->IsRadioNuclide()) elemrn = (TGeoElementRN*)elem;
/// ~~~
/// One can get Ni/N1(t=0) at any moment of time. Ni is the number of atoms
/// of one of the decay products, N1(0) is the number of atoms of the first top
/// element at t=0.
/// ~~~ {.cpp}
///    Double_t fraction_weight = elemrn->Ratio()->Concentration(time);
/// ~~~
/// One can also display the time evolution of the fractional weight:
/// ~~~ {.cpp}
///    elemrn->Ratio()->Draw(option);
/// ~~~

void TGeoMixture::FillMaterialEvolution(TObjArray *population, Double_t precision)
{
   if (population->GetEntriesFast()) {
      Error("FillMaterialEvolution", "Provide an empty array !");
      return;
   }
   TGeoElementTable *table = gGeoManager->GetElementTable();
   TGeoElement *elem;
   TGeoElementRN *elemrn;
   TIter next(table->GetElementsRN());
   while ((elemrn=(TGeoElementRN*)next())) elemrn->ResetRatio();
   Double_t factor;
   for (Int_t i=0; i<fNelements; i++) {
      elem = GetElement(i);
      if (!elem->IsRadioNuclide()) {
         population->Add(elem);
         continue;
      }
      elemrn = (TGeoElementRN*)elem;
      factor = fWeights[i]*fAmixture[0]/(fWeights[0]*fAmixture[i]);
      elemrn->FillPopulation(population, precision, factor);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// static function
///  Compute screening factor for pair production and Bremsstrahlung
///  REFERENCE : EGS MANUAL SLAC 210 - UC32 - JUNE 78
///                        FORMULA 2.7.22

Double_t TGeoMaterial::ScreenFactor(Double_t z)
{
   const Double_t al183= 5.20948 , al1440 = 7.27239;
   Double_t alz  = TMath::Log(z)/3.;
   Double_t factor = (al1440 - 2*alz) / (al183 - alz - TGeoMaterial::Coulomb(z));
   return factor;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute Derived Quantities as in Geant4

void TGeoMixture::ComputeDerivedQuantities()
{
   const Double_t Na = (TGeoManager::GetDefaultUnits()==TGeoManager::kRootUnits)
     ? TGeoUnit::Avogadro : TGeant4Unit::Avogadro;

   if ( fVecNbOfAtomsPerVolume ) delete [] fVecNbOfAtomsPerVolume;

   fVecNbOfAtomsPerVolume = new Double_t[fNelements];

   // Formula taken from G4Material.cxx L312
   for (Int_t i=0; i<fNelements; ++i) {
      fVecNbOfAtomsPerVolume[i] = Na*fDensity*fWeights[i]/((TGeoElement*)fElements->At(i))->A();
   }
   ComputeRadiationLength();
   ComputeNuclearInterLength();
}


////////////////////////////////////////////////////////////////////////////////
/// Compute Radiation Length based on Geant4 formula

void TGeoMixture::ComputeRadiationLength()
{
   // Formula taken from G4Material.cxx L556
   const Double_t cm = (TGeoManager::GetDefaultUnits()==TGeoManager::kRootUnits) ? TGeoUnit::cm : TGeant4Unit::cm;
   Double_t radinv = 0.0 ;
   for (Int_t i=0;i<fNelements;++i) {
     radinv += fVecNbOfAtomsPerVolume[i]*((TGeoElement*)fElements->At(i))->GetfRadTsai();
   }
   fRadLen = (radinv <= 0.0 ? DBL_MAX : cm/radinv);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute Nuclear Interaction Length based on Geant4 formula
void TGeoMixture::ComputeNuclearInterLength()
{
   // Formula taken from G4Material.cxx L567
   TGeoManager::EDefaultUnits typ = TGeoManager::GetDefaultUnits();
   const Double_t g   = (typ==TGeoManager::kRootUnits) ? TGeoUnit::g   : TGeant4Unit::g;
   const Double_t cm  = (typ==TGeoManager::kRootUnits) ? TGeoUnit::cm  : TGeant4Unit::cm;
   const Double_t amu = (typ==TGeoManager::kRootUnits) ? TGeoUnit::amu : TGeant4Unit::amu;
   const Double_t lambda0  = 35*g/(cm*cm);
   const Double_t twothird = 2.0/3.0;
   Double_t NILinv = 0.0;
   for (Int_t i=0; i<fNelements; ++i) {
      Int_t Z = static_cast<Int_t>(((TGeoElement*)fElements->At(i))->Z()+0.5);
      Double_t A = ((TGeoElement*)fElements->At(i))->Neff();
      if(1 == Z) {
         NILinv += fVecNbOfAtomsPerVolume[i]*A;
      } else {
         NILinv += fVecNbOfAtomsPerVolume[i]*TMath::Exp(twothird*TMath::Log(A));
      }
   }
   NILinv *= amu/lambda0;
   fIntLen = (NILinv <= 0.0 ? DBL_MAX : cm/NILinv);
}
