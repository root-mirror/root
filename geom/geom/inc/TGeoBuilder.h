// @(#)root/geom:$Id$
// Author: Mihaela Gheata   30/05/07

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoBuilder
#define ROOT_TGeoBuilder

#include "TObject.h"

class TGeoMaterial;
class TGeoMatrix;
class TGeoMedium;
class TGeoShape;
class TGeoVolume;
class TGeoVolumeAssembly;
class TGeoVolumeMulti;
class TGeoManager;

class TGeoBuilder : public TObject
{
protected:
   static TGeoBuilder    *fgInstance;            //! static pointer to singleton

   TGeoBuilder();
   TGeoBuilder(const TGeoBuilder&) = delete;
   TGeoBuilder& operator=(const TGeoBuilder&) = delete;

private :
   TGeoManager           *fGeometry;             //! current geometry

   void                   SetGeometry(TGeoManager *geom) {fGeometry = geom;}

public :
   virtual ~TGeoBuilder();

   static TGeoBuilder    *Instance(TGeoManager *geom);

   Int_t                  AddMaterial(TGeoMaterial *material);
   Int_t                  AddTransformation(TGeoMatrix *matrix);
   Int_t                  AddShape(TGeoShape *shape);
   void                   RegisterMatrix(TGeoMatrix *matrix);

   TGeoVolume            *MakeArb8(const char *name, TGeoMedium *medium,
                                     Double_t dz, Double_t *vertices=0);
   TGeoVolume            *MakeBox(const char *name, TGeoMedium *medium,
                                     Double_t dx, Double_t dy, Double_t dz);
   TGeoVolume            *MakeCone(const char *name, TGeoMedium *medium,
                                      Double_t dz, Double_t rmin1, Double_t rmax1,
                                      Double_t rmin2, Double_t rmax2);
   TGeoVolume            *MakeCons(const char *name, TGeoMedium *medium,
                                      Double_t dz, Double_t rmin1, Double_t rmax1,
                                      Double_t rmin2, Double_t rmax2,
                                      Double_t phi1, Double_t phi2);
   TGeoVolume            *MakeCtub(const char *name, TGeoMedium *medium,
                                      Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2,
                                      Double_t lx, Double_t ly, Double_t lz, Double_t tx, Double_t ty, Double_t tz);
   TGeoVolume            *MakeEltu(const char *name, TGeoMedium *medium,
                                      Double_t a, Double_t b, Double_t dz);
   TGeoVolume            *MakeGtra(const char *name, TGeoMedium *medium,
                                   Double_t dz, Double_t theta, Double_t phi, Double_t twist, Double_t h1,
                                   Double_t bl1, Double_t tl1, Double_t alpha1, Double_t h2, Double_t bl2,
                                   Double_t tl2, Double_t alpha2);
   TGeoVolume            *MakePara(const char *name, TGeoMedium *medium,
                                     Double_t dx, Double_t dy, Double_t dz,
                                     Double_t alpha, Double_t theta, Double_t phi);
   TGeoVolume            *MakePcon(const char *name, TGeoMedium *medium,
                                      Double_t phi, Double_t dphi, Int_t nz);
   TGeoVolume            *MakeParaboloid(const char *name, TGeoMedium *medium,
                                      Double_t rlo, Double_t rhi, Double_t dz);
   TGeoVolume            *MakeHype(const char *name, TGeoMedium *medium,
                                      Double_t rin, Double_t stin, Double_t rout, Double_t stout, Double_t dz);
   TGeoVolume            *MakePgon(const char *name, TGeoMedium *medium,
                                      Double_t phi, Double_t dphi, Int_t nedges, Int_t nz);
   TGeoVolume            *MakeSphere(const char *name, TGeoMedium *medium,
                                     Double_t rmin, Double_t rmax,
                                     Double_t themin=0, Double_t themax=180,
                                     Double_t phimin=0, Double_t phimax=360);
   TGeoVolume            *MakeTorus(const char *name, TGeoMedium *medium, Double_t r,
                                    Double_t rmin, Double_t rmax, Double_t phi1=0, Double_t dphi=360);
   TGeoVolume            *MakeTrap(const char *name, TGeoMedium *medium,
                                   Double_t dz, Double_t theta, Double_t phi, Double_t h1,
                                   Double_t bl1, Double_t tl1, Double_t alpha1, Double_t h2, Double_t bl2,
                                   Double_t tl2, Double_t alpha2);
   TGeoVolume            *MakeTrd1(const char *name, TGeoMedium *medium,
                                      Double_t dx1, Double_t dx2, Double_t dy, Double_t dz);
   TGeoVolume            *MakeTrd2(const char *name, TGeoMedium *medium,
                                      Double_t dx1, Double_t dx2, Double_t dy1, Double_t dy2,
                                      Double_t dz);
   TGeoVolume            *MakeTube(const char *name, TGeoMedium *medium,
                                      Double_t rmin, Double_t rmax, Double_t dz);
   TGeoVolume            *MakeTubs(const char *name, TGeoMedium *medium,
                                      Double_t rmin, Double_t rmax, Double_t dz,
                                      Double_t phi1, Double_t phi2);
   TGeoVolume            *MakeXtru(const char *name, TGeoMedium *medium,
                                   Int_t nz);
   TGeoVolumeAssembly    *MakeVolumeAssembly(const char *name);
   TGeoVolumeMulti       *MakeVolumeMulti(const char *name, TGeoMedium *medium);

   //--- GEANT3-like geometry creation
   TGeoVolume            *Division(const char *name, const char *mother, Int_t iaxis, Int_t ndiv,
                                         Double_t start, Double_t step, Int_t numed=0, Option_t *option="");
   void                   Matrix(Int_t index, Double_t theta1, Double_t phi1,
                                       Double_t theta2, Double_t phi2,
                                       Double_t theta3, Double_t phi3);
   TGeoMaterial          *Material(const char *name, Double_t a, Double_t z, Double_t dens, Int_t uid, Double_t radlen=0, Double_t intlen=0);
   TGeoMaterial          *Mixture(const char *name, Float_t *a, Float_t *z, Double_t dens,
                                        Int_t nelem, Float_t *wmat, Int_t uid);
   TGeoMaterial          *Mixture(const char *name, Double_t *a, Double_t *z, Double_t dens,
                                        Int_t nelem, Double_t *wmat, Int_t uid);
   TGeoMedium            *Medium(const char *name, Int_t numed, Int_t nmat, Int_t isvol,
                                       Int_t ifield, Double_t fieldm, Double_t tmaxfd,
                                       Double_t stemax, Double_t deemax, Double_t epsil,
                                       Double_t stmin);
   void                   Node(const char *name, Int_t nr, const char *mother,
                                     Double_t x, Double_t y, Double_t z, Int_t irot,
                                     Bool_t isOnly, Float_t *upar, Int_t npar=0);
   void                   Node(const char *name, Int_t nr, const char *mother,
                                     Double_t x, Double_t y, Double_t z, Int_t irot,
                                     Bool_t isOnly, Double_t *upar, Int_t npar=0);
   TGeoVolume            *Volume(const char *name, const char *shape, Int_t nmed,
                                       Float_t *upar, Int_t npar=0);
   TGeoVolume            *Volume(const char *name, const char *shape, Int_t nmed,
                                       Double_t *upar, Int_t npar=0);

   ClassDef(TGeoBuilder, 0)          // geometry builder singleton
};

#endif

