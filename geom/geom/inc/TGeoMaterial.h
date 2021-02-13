// @(#)root/geom:$Id$
// Author: Andrei Gheata   25/10/01

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoMaterial
#define ROOT_TGeoMaterial

#include <TNamed.h>
#include <TAttFill.h>
#include <TList.h>

#ifdef R__LESS_INCLUDES
class TGeoElement;
class TGeoElementTable;
#else
#include "TGeoElement.h"
#endif

// forward declarations
class TGeoExtension;
class TGDMLMatrix;

// Some units used in G4
static const Double_t STP_temperature = 273.15;     // [K]
static const Double_t STP_pressure    = 6.32420e+8; // [MeV/mm3]

class TGeoMaterial : public TNamed,
                     public TAttFill
{
public:
   enum EGeoMaterial {
      kMatUsed   =   BIT(17),
      kMatSavePrimitive = BIT(18)
   };
   enum EGeoMaterialState {
      kMatStateUndefined,
      kMatStateSolid,
      kMatStateLiquid,
      kMatStateGas
   };

protected:
   Int_t                    fIndex;      // material index
   Double_t                 fA;          // A of material
   Double_t                 fZ;          // Z of material
   Double_t                 fDensity;    // density of material
   Double_t                 fRadLen;     // radiation length
   Double_t                 fIntLen;     // interaction length
   Double_t                 fTemperature; // temperature
   Double_t                 fPressure;   // pressure
   EGeoMaterialState        fState;      // material state
   TObject                 *fShader;     // shader with optical properties
   TObject                 *fCerenkov;   // pointer to class with Cerenkov properties
   TGeoElement             *fElement;    // pointer to element composing the material
   TList                    fProperties; // user-defined properties
   TList                    fConstProperties; // user-defined constant properties
   TGeoExtension           *fUserExtension;  //! Transient user-defined extension to materials
   TGeoExtension           *fFWExtension;    //! Transient framework-defined extension to materials

// methods
   TGeoMaterial(const TGeoMaterial&);
   TGeoMaterial& operator=(const TGeoMaterial&);


public:
   // constructors
   TGeoMaterial();
   TGeoMaterial(const char *name);
   TGeoMaterial(const char *name, Double_t a, Double_t z,
                Double_t rho, Double_t radlen=0, Double_t intlen=0);
   TGeoMaterial(const char *name, Double_t a, Double_t z, Double_t rho,
                EGeoMaterialState state, Double_t temperature=STP_temperature, Double_t pressure=STP_pressure);
   TGeoMaterial(const char *name, TGeoElement *elem, Double_t rho);

   // destructor
   virtual ~TGeoMaterial();
   // methods
   static  Double_t         Coulomb(Double_t z);
   // radioactive mixture evolution
   virtual TGeoMaterial    *DecayMaterial(Double_t time, Double_t precision=0.001);
   virtual void             FillMaterialEvolution(TObjArray *population, Double_t precision=0.001);
   // getters & setters
   bool                     AddProperty(const char *property, const char *ref);
   bool                     AddConstProperty(const char *property, const char *ref);
   Int_t                    GetNproperties() const { return fProperties.GetSize(); }
   Int_t                    GetNconstProperties() const { return fConstProperties.GetSize(); }
   const char              *GetPropertyRef(const char *property) const;
   const char              *GetPropertyRef(Int_t i) const { return (fProperties.At(i) ? fProperties.At(i)->GetTitle() : nullptr); }
   Double_t                 GetConstProperty(const char *property, Bool_t *error = nullptr) const;
   Double_t                 GetConstProperty(Int_t i, Bool_t *error = nullptr) const;
   const char              *GetConstPropertyRef(const char *property) const;
   const char              *GetConstPropertyRef(Int_t i) const { return (fConstProperties.At(i) ? fConstProperties.At(i)->GetTitle() : nullptr); }
   TList const             &GetProperties() const { return fProperties; }
   TList const             &GetConstProperties() const { return fConstProperties; }
   TGDMLMatrix*             GetProperty(const char* name)  const;
   TGDMLMatrix*             GetProperty(Int_t i)  const;
   virtual Int_t            GetByteCount() const {return sizeof(*this);}
   virtual Double_t         GetA() const       {return fA;}
   virtual Double_t         GetZ()  const      {return fZ;}
   virtual Int_t            GetDefaultColor() const;
   virtual Double_t         GetDensity() const {return fDensity;}
   virtual Int_t            GetNelements() const {return 1;}
   TGeoElement             *GetElement() const;
   virtual TGeoElement     *GetElement(Int_t i) const;
   virtual void             GetElementProp(Double_t &a, Double_t &z, Double_t &w, Int_t i=0);
   TGeoElement             *GetBaseElement() const {return fElement;}
   char                    *GetPointerName() const;
   virtual Double_t         GetRadLen() const  {return fRadLen;}
   virtual Double_t         GetIntLen() const  {return fIntLen;}
   Int_t                    GetIndex();
   virtual TObject         *GetCerenkovProperties() const {return fCerenkov;}
   Char_t                   GetTransparency() const {return (fFillStyle<3000 || fFillStyle>3100)?0:Char_t(fFillStyle-3000);}
   Double_t                 GetTemperature() const {return fTemperature;}
   Double_t                 GetPressure() const {return fPressure;}
   EGeoMaterialState        GetState() const {return fState;}
   virtual Double_t         GetSpecificActivity(Int_t) const {return 0.;}
   TGeoExtension           *GetUserExtension() const {return fUserExtension;}
   TGeoExtension           *GetFWExtension() const   {return fFWExtension;}
   TGeoExtension           *GrabUserExtension() const;
   TGeoExtension           *GrabFWExtension() const;
   virtual Bool_t           IsEq(const TGeoMaterial *other) const;
   Bool_t                   IsUsed() const {return TObject::TestBit(kMatUsed);}
   virtual Bool_t           IsMixture() const {return kFALSE;}
   virtual void             Print(const Option_t *option="") const;
   virtual void             SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void             SetA(Double_t a) {fA = a; SetRadLen(0);}
   virtual void             SetZ(Double_t z) {fZ = z; SetRadLen(0);}
   virtual void             SetDensity(Double_t density) {fDensity = density; SetRadLen(0);}
   void                     SetIndex(Int_t index) {fIndex=index;}
   virtual void             SetCerenkovProperties(TObject* cerenkov) {fCerenkov = cerenkov;}
   void                     SetRadLen(Double_t radlen, Double_t intlen=0.);
   void                     SetUsed(Bool_t flag=kTRUE) {TObject::SetBit(kMatUsed, flag);}
   void                     SetTransparency(Char_t transparency=0) {fFillStyle = 3000+transparency;}
   void                     SetTemperature(Double_t temperature) {fTemperature = temperature;}
   void                     SetPressure(Double_t pressure) {fPressure = pressure;}
   void                     SetState(EGeoMaterialState state) {fState = state;}
   void                     SetUserExtension(TGeoExtension *ext);
   void                     SetFWExtension(TGeoExtension *ext);
   static  Double_t         ScreenFactor(Double_t z);



   ClassDef(TGeoMaterial, 7)              // base material class

//***** Need to add classes and globals to LinkDef.h *****
};


class TGeoMixture : public TGeoMaterial
{
protected :
// data members
   Int_t                    fNelements;  // number of elements
   Double_t                *fZmixture;   // [fNelements] array of Z of the elements
   Double_t                *fAmixture;   // [fNelements] array of A of the elements
   Double_t                *fWeights;    // [fNelements] array of relative proportions by mass
   Int_t                   *fNatoms;     // [fNelements] array of numbers of atoms
   Double_t                *fVecNbOfAtomsPerVolume; //[fNelements] array of numbers of atoms per unit volume
   TObjArray               *fElements;   // array of elements composing the mixture
// methods
   TGeoMixture(const TGeoMixture&) = delete;
   TGeoMixture& operator=(const TGeoMixture&) = delete;
   void                     AverageProperties();

public:
   // constructors
   TGeoMixture();
   TGeoMixture(const char *name, Int_t nel, Double_t rho=-1);
   // destructor
   virtual ~TGeoMixture();
   // methods for adding elements
   void                     AddElement(Double_t a, Double_t z, Double_t weight);
   void                     AddElement(TGeoMaterial *mat, Double_t weight);
   void                     AddElement(TGeoElement *elem, Double_t weight);
   void                     AddElement(TGeoElement *elem, Int_t natoms);
   // backward compatibility for defining elements
   void                     DefineElement(Int_t iel, Double_t a, Double_t z, Double_t weight);
   void                     DefineElement(Int_t iel, TGeoElement *elem, Double_t weight);
   void                     DefineElement(Int_t iel, Int_t z, Int_t natoms);
   // radioactive mixture evolution
   virtual TGeoMaterial    *DecayMaterial(Double_t time, Double_t precision=0.001);
   virtual void             FillMaterialEvolution(TObjArray *population, Double_t precision=0.001);
   // getters
   virtual Int_t            GetByteCount() const {return 48+12*fNelements;}
   virtual TGeoElement     *GetElement(Int_t i=0) const;
   virtual void             GetElementProp(Double_t &a, Double_t &z, Double_t &w, Int_t i=0) {a=fAmixture[i]; z=fZmixture[i]; w=fWeights[i];}
   virtual Int_t            GetNelements() const {return fNelements;}
   Double_t                *GetZmixt() const     {return fZmixture;}
   Double_t                *GetAmixt() const     {return fAmixture;}
   Double_t                *GetWmixt() const     {return fWeights;}
   Int_t                   *GetNmixt() const     {return fNatoms;}
   virtual Double_t         GetSpecificActivity(Int_t i=-1) const;
   // utilities
   virtual Bool_t           IsEq(const TGeoMaterial *other) const;
   virtual Bool_t           IsMixture() const {return kTRUE;}
   virtual void             Print(const Option_t *option="") const;
   virtual void             SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void             SetA(Double_t a) {fA = a;}
   virtual void             SetZ(Double_t z) {fZ = z;}
   virtual void             SetDensity(Double_t density) {fDensity = density; AverageProperties();}
   void                     ComputeDerivedQuantities();
   void                     ComputeRadiationLength();
   void                     ComputeNuclearInterLength();

   ClassDef(TGeoMixture, 3)              // material mixtures
};

inline void TGeoMixture::DefineElement(Int_t, Double_t a, Double_t z, Double_t weight)
   {return AddElement(a,z,weight);}
inline void TGeoMixture::DefineElement(Int_t, TGeoElement *elem, Double_t weight)
   {return AddElement(elem,weight);}


#endif

