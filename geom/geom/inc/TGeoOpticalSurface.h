// @(#)root/geom:$Id$
// Author: Andrei Gheata   05/12/18

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoOpticalSurface
#define ROOT_TGeoOpticalSurface

#include <TNamed.h>
#include <TList.h>

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoOpticalSurface - class describing surface properties for           //
//                      compatibility with Geant4                         //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGDMLMatrix;

class TGeoOpticalSurface : public TNamed {
public:
   enum ESurfaceFinish {
      kFpolished,             // smooth perfectly polished surface
      kFpolishedfrontpainted, // smooth top-layer (front) paint
      kFpolishedbackpainted,  // same is 'polished' but with a back-paint

      kFground,             // rough surface
      kFgroundfrontpainted, // rough top-layer (front) paint
      kFgroundbackpainted,  // same as 'ground' but with a back-paint

      kFpolishedlumirrorair,  // mechanically polished surface, with lumirror
      kFpolishedlumirrorglue, // mechanically polished surface, with lumirror & meltmount
      kFpolishedair,          // mechanically polished surface
      kFpolishedteflonair,    // mechanically polished surface, with teflon
      kFpolishedtioair,       // mechanically polished surface, with tio paint
      kFpolishedtyvekair,     // mechanically polished surface, with tyvek
      kFpolishedvm2000air,    // mechanically polished surface, with esr film
      kFpolishedvm2000glue,   // mechanically polished surface, with esr film & meltmount

      kFetchedlumirrorair,  // chemically etched surface, with lumirror
      kFetchedlumirrorglue, // chemically etched surface, with lumirror & meltmount
      kFetchedair,          // chemically etched surface
      kFetchedteflonair,    // chemically etched surface, with teflon
      kFetchedtioair,       // chemically etched surface, with tio paint
      kFetchedtyvekair,     // chemically etched surface, with tyvek
      kFetchedvm2000air,    // chemically etched surface, with esr film
      kFetchedvm2000glue,   // chemically etched surface, with esr film & meltmount

      kFgroundlumirrorair,  // rough-cut surface, with lumirror
      kFgroundlumirrorglue, // rough-cut surface, with lumirror & meltmount
      kFgroundair,          // rough-cut surface
      kFgroundteflonair,    // rough-cut surface, with teflon
      kFgroundtioair,       // rough-cut surface, with tio paint
      kFgroundtyvekair,     // rough-cut surface, with tyvek
      kFgroundvm2000air,    // rough-cut surface, with esr film
      kFgroundvm2000glue,   // rough-cut surface, with esr film & meltmount

      // for DAVIS model
      kFRough_LUT,             // rough surface
      kFRoughTeflon_LUT,       // rough surface wrapped in Teflon tape
      kFRoughESR_LUT,          // rough surface wrapped with ESR
      kFRoughESRGrease_LUT,    // rough surface wrapped with ESR and coupled with opical grease
      kFPolished_LUT,          // polished surface
      kFPolishedTeflon_LUT,    // polished surface wrapped in Teflon tape
      kFPolishedESR_LUT,       // polished surface wrapped with ESR
      kFPolishedESRGrease_LUT, // polished surface wrapped with ESR and coupled with opical grease
      kFDetector_LUT           // polished surface with optical grease
   };

   enum ESurfaceModel {
      kMglisur,  // original GEANT3 model
      kMunified, // UNIFIED model
      kMLUT,     // Look-Up-Table model
      kMDAVIS,   // DAVIS model
      kMdichroic // dichroic filter
   };

   enum ESurfaceType {
      kTdielectric_metal,      // dielectric-metal interface
      kTdielectric_dielectric, // dielectric-dielectric interface
      kTdielectric_LUT,        // dielectric-Look-Up-Table interface
      kTdielectric_LUTDAVIS,   // dielectric-Look-Up-Table DAVIS interface
      kTdielectric_dichroic,   // dichroic filter interface
      kTfirsov,                // for Firsov Process
      kTx_ray                  // for x-ray mirror process
   };

private:
   std::string fName = "";                  // Surface name
   ESurfaceType fType = kTdielectric_metal; // Surface type
   ESurfaceModel fModel = kMglisur;         // Surface model
   ESurfaceFinish fFinish = kFpolished;     // Surface finish

   Double_t fValue = 0.0;      // The value used to determine sigmaalpha and polish
   Double_t fSigmaAlpha = 0.0; // The sigma of micro-facet polar angle
   Double_t fPolish = 0.0;     // Polish parameter in glisur model

   TList fProperties; // List of surface properties

   // No copy
   TGeoOpticalSurface(const TGeoOpticalSurface &) = delete;
   TGeoOpticalSurface &operator=(const TGeoOpticalSurface &) = delete;

public:
   // constructors
   TGeoOpticalSurface() {}

   TGeoOpticalSurface(const char *name, ESurfaceModel model = kMglisur, ESurfaceFinish finish = kFpolished,
                      ESurfaceType type = kTdielectric_dielectric, Double_t value = 1.0);

   virtual ~TGeoOpticalSurface() {}

   // Accessors
   bool AddProperty(const char *property, const char *ref);
   const char *GetPropertyRef(const char *property);
   TList const &GetProperties() const { return fProperties; }
   Int_t GetNproperties() const { return fProperties.GetSize(); }
   TGDMLMatrix* GetProperty(const char* name)  const;
   TGDMLMatrix* GetProperty(Int_t i)  const;
   ESurfaceType GetType() const { return fType; }
   ESurfaceModel GetModel() const { return fModel; }
   ESurfaceFinish GetFinish() const { return fFinish; }
   Double_t GetPolish() const { return fPolish; }
   Double_t GetValue() const { return fValue; }
   Double_t GetSigmaAlpha() const { return fSigmaAlpha; }

   void SetType(ESurfaceType type) { fType = type; }
   void SetModel(ESurfaceModel model) { fModel = model; }
   void SetFinish(ESurfaceFinish finish) { fFinish = finish; }
   void SetPolish(Double_t polish) { fPolish = polish; }
   void SetValue(Double_t value) { fValue = value; }
   void SetSigmaAlpha(Double_t sigmaalpha) { fSigmaAlpha = sigmaalpha; }

   void Print(Option_t *option = "") const;

   static ESurfaceType StringToType(const char *type);
   static const char *TypeToString(ESurfaceType type);
   static ESurfaceModel StringToModel(const char *model);
   static const char *ModelToString(ESurfaceModel model);
   static ESurfaceFinish StringToFinish(const char *finish);
   static const char *FinishToString(ESurfaceFinish finish);

   ClassDef(TGeoOpticalSurface, 1) // Class representing an optical surface
};

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoSkinSurface - class describing a surface having optical properties //
//                      surrounding a volume                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoVolume;

class TGeoSkinSurface : public TNamed {
private:
   TGeoOpticalSurface const *fSurface = nullptr; // Referenced optical surface
   TGeoVolume const *fVolume = nullptr;          // Referenced volume
public:
   TGeoSkinSurface() {}
   TGeoSkinSurface(const char *name, const char *ref, TGeoOpticalSurface const *surf, TGeoVolume const *vol)
      : TNamed(name, ref), fSurface(surf), fVolume(vol)
   {
   }
   virtual ~TGeoSkinSurface() {}

   TGeoOpticalSurface const *GetSurface() const { return fSurface; }
   TGeoVolume const *GetVolume() const { return fVolume; }

   void Print(Option_t *option = "") const;

   ClassDef(TGeoSkinSurface, 1) // A surface with optical properties surrounding a volume
};

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoBorderSurface - class describing a surface having optical          //
//                      properties between 2 touching volumes             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoNode;

class TGeoBorderSurface : public TNamed {
private:
   TGeoOpticalSurface const *fSurface = nullptr; // Referenced optical surface
   TGeoNode const *fNode1 = nullptr;             // Referenced node 1
   TGeoNode const *fNode2 = nullptr;             // Referenced node 2
public:
   TGeoBorderSurface() {}
   TGeoBorderSurface(const char *name, const char *ref, TGeoOpticalSurface const *surf, TGeoNode const *node1,
                     TGeoNode const *node2)
      : TNamed(name, ref), fSurface(surf), fNode1(node1), fNode2(node2)
   {
   }
   virtual ~TGeoBorderSurface() {}

   TGeoOpticalSurface const *GetSurface() const { return fSurface; }
   TGeoNode const *GetNode1() const { return fNode1; }
   TGeoNode const *GetNode2() const { return fNode2; }

   void Print(Option_t *option = "") const;

   ClassDef(TGeoBorderSurface, 1) // A surface with optical properties betwqeen 2 touching volumes
};

#endif // ROOT_TGeoOpticalSurface
