// @(#)root/geom:$Id$
// Author: Mihaela Gheata   30/05/07

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoBuilder
\ingroup Geometry_classes

Utility class for creating geometry objects.These will be associated
with the current selected geometry manager object:

`TGeoBuilder::Instance()->SetGeometry(gGeoManager);`

The geometry builder is a singleton that may be used to build one or more
geometries.
*/

#include "TList.h"
#include "TObjArray.h"

#include "TGeoManager.h"
#include "TGeoElement.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoMatrix.h"
#include "TGeoPara.h"
#include "TGeoParaboloid.h"
#include "TGeoTube.h"
#include "TGeoEltu.h"
#include "TGeoHype.h"
#include "TGeoCone.h"
#include "TGeoSphere.h"
#include "TGeoArb8.h"
#include "TGeoPgon.h"
#include "TGeoTrd1.h"
#include "TGeoTrd2.h"
#include "TGeoTorus.h"
#include "TGeoXtru.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"

#include "TGeoBuilder.h"

ClassImp(TGeoBuilder);

TGeoBuilder *TGeoBuilder::fgInstance = NULL;

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TGeoBuilder::TGeoBuilder()
            :fGeometry(nullptr)
{
   fgInstance = this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoBuilder::~TGeoBuilder()
{
   fgInstance = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to singleton.

TGeoBuilder *TGeoBuilder::Instance(TGeoManager *geom)
{
   if (!geom) {
      printf("ERROR: Cannot create geometry builder with NULL geometry\n");
      return NULL;
   }
   if (!fgInstance) fgInstance = new TGeoBuilder();
   fgInstance->SetGeometry(geom);
   return fgInstance;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a material to the list. Returns index of the material in list.

Int_t TGeoBuilder::AddMaterial(TGeoMaterial *material)
{
   if (!material) return -1;
   TList *materials = fGeometry->GetListOfMaterials();
   Int_t index = materials->GetSize();
   material->SetIndex(index);
   materials->Add(material);
   return index;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a matrix to the list. Returns index of the matrix in list.

Int_t TGeoBuilder::AddTransformation(TGeoMatrix *matrix)
{
   Int_t index = -1;
   if (!matrix) return -1;
   TObjArray *matrices = fGeometry->GetListOfMatrices();
   index = matrices->GetEntriesFast();
   matrices->AddAtAndExpand(matrix,index);
   return index;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a shape to the list. Returns index of the shape in list.

Int_t TGeoBuilder::AddShape(TGeoShape *shape)
{
   Int_t index = -1;
   if (!shape) return -1;
   TObjArray *shapes = fGeometry->GetListOfShapes();
   if (shape->IsRunTimeShape()) shapes =  fGeometry->GetListOfGShapes();
   index = shapes->GetEntriesFast();
   shapes->AddAtAndExpand(shape,index);
   return index;
}

////////////////////////////////////////////////////////////////////////////////
/// Register a matrix to the list of matrices. It will be cleaned-up at the
/// destruction TGeoManager.

void TGeoBuilder::RegisterMatrix(TGeoMatrix *matrix)
{
   if (matrix->IsRegistered()) return;
   TObjArray *matrices = fGeometry->GetListOfMatrices();
   Int_t nmat = matrices->GetEntriesFast();
   matrices->AddAtAndExpand(matrix, nmat);
}

////////////////////////////////////////////////////////////////////////////////
/// Make an TGeoArb8 volume.

TGeoVolume *TGeoBuilder::MakeArb8(const char *name, TGeoMedium *medium,
                                  Double_t dz, Double_t *vertices)
{
   TGeoArb8 *arb = new TGeoArb8(name, dz, vertices);
   TGeoVolume *vol = new TGeoVolume(name, arb, medium);
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a box shape with given medium.

TGeoVolume *TGeoBuilder::MakeBox(const char *name, TGeoMedium *medium,
                                    Double_t dx, Double_t dy, Double_t dz)
{
   TGeoBBox *box = new TGeoBBox(name, dx, dy, dz);
   TGeoVolume *vol = 0;
   if (box->IsRunTimeShape()) {
      vol = fGeometry->MakeVolumeMulti(name, medium);
      vol->SetShape(box);
   } else {
      vol = new TGeoVolume(name, box, medium);
   }
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a parallelepiped shape with given medium.

TGeoVolume *TGeoBuilder::MakePara(const char *name, TGeoMedium *medium,
                                    Double_t dx, Double_t dy, Double_t dz,
                                    Double_t alpha, Double_t theta, Double_t phi)
{
   if (TMath::Abs(alpha)<TGeoShape::Tolerance() && TMath::Abs(theta)<TGeoShape::Tolerance()) {
      Warning("MakePara","parallelepiped %s having alpha=0, theta=0 -> making box instead", name);
      return MakeBox(name, medium, dx, dy, dz);
   }
   TGeoPara *para=0;
   para = new TGeoPara(name, dx, dy, dz, alpha, theta, phi);
   TGeoVolume *vol = 0;
   if (para->IsRunTimeShape()) {
      vol = fGeometry->MakeVolumeMulti(name, medium);
      vol->SetShape(para);
   } else {
      vol = new TGeoVolume(name, para, medium);
   }
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a sphere shape with given medium

TGeoVolume *TGeoBuilder::MakeSphere(const char *name, TGeoMedium *medium,
                                    Double_t rmin, Double_t rmax, Double_t themin, Double_t themax,
                                    Double_t phimin, Double_t phimax)
{
   TGeoSphere *sph = new TGeoSphere(name, rmin, rmax, themin, themax, phimin, phimax);
   TGeoVolume *vol = new TGeoVolume(name, sph, medium);
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a torus shape with given medium.

TGeoVolume *TGeoBuilder::MakeTorus(const char *name, TGeoMedium *medium, Double_t r,
                                   Double_t rmin, Double_t rmax, Double_t phi1, Double_t dphi)
{
   TGeoTorus *tor = new TGeoTorus(name,r,rmin,rmax,phi1,dphi);
   TGeoVolume *vol = new TGeoVolume(name, tor, medium);
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a tube shape with given medium.

TGeoVolume *TGeoBuilder::MakeTube(const char *name, TGeoMedium *medium,
                                     Double_t rmin, Double_t rmax, Double_t dz)
{
   if (rmin>rmax) {
      Error("MakeTube", "tube %s, Rmin=%g greater than Rmax=%g", name,rmin,rmax);
   }
   TGeoTube *tube = new TGeoTube(name, rmin, rmax, dz);
   TGeoVolume *vol = 0;
   if (tube->IsRunTimeShape()) {
      vol = fGeometry->MakeVolumeMulti(name, medium);
      vol->SetShape(tube);
   } else {
      vol = new TGeoVolume(name, tube, medium);
   }
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a tube segment shape with given medium.

TGeoVolume *TGeoBuilder::MakeTubs(const char *name, TGeoMedium *medium,
                                     Double_t rmin, Double_t rmax, Double_t dz,
                                     Double_t phi1, Double_t phi2)
{
   TGeoTubeSeg *tubs = new TGeoTubeSeg(name, rmin, rmax, dz, phi1, phi2);
   TGeoVolume *vol = 0;
   if (tubs->IsRunTimeShape()) {
      vol = fGeometry->MakeVolumeMulti(name, medium);
      vol->SetShape(tubs);
   } else {
      vol = new TGeoVolume(name, tubs, medium);
   }
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a tube shape with given medium

TGeoVolume *TGeoBuilder::MakeEltu(const char *name, TGeoMedium *medium,
                                     Double_t a, Double_t b, Double_t dz)
{
   TGeoEltu *eltu = new TGeoEltu(name, a, b, dz);
   TGeoVolume *vol = 0;
   if (eltu->IsRunTimeShape()) {
      vol = fGeometry->MakeVolumeMulti(name, medium);
      vol->SetShape(eltu);
   } else {
      vol = new TGeoVolume(name, eltu, medium);
   }
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a tube shape with given medium

TGeoVolume *TGeoBuilder::MakeHype(const char *name, TGeoMedium *medium,
                                        Double_t rin, Double_t stin, Double_t rout, Double_t stout, Double_t dz)
{
   TGeoHype * hype =  new TGeoHype(name, rin,stin,rout,stout,dz);
   TGeoVolume *vol = 0;
   if (hype->IsRunTimeShape()) {
      vol = fGeometry->MakeVolumeMulti(name, medium);
      vol->SetShape(hype);
   } else {
      vol = new TGeoVolume(name, hype, medium);
   }
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a tube shape with given medium

TGeoVolume *TGeoBuilder::MakeParaboloid(const char *name, TGeoMedium *medium,
                                        Double_t rlo, Double_t rhi, Double_t dz)
{
   TGeoParaboloid *parab = new TGeoParaboloid(name, rlo, rhi, dz);
   TGeoVolume *vol = 0;
   if (parab->IsRunTimeShape()) {
      vol = fGeometry->MakeVolumeMulti(name, medium);
      vol->SetShape(parab);
   } else {
      vol = new TGeoVolume(name, parab, medium);
   }
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a tube segment shape with given medium

TGeoVolume *TGeoBuilder::MakeCtub(const char *name, TGeoMedium *medium,
                                     Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2,
                                     Double_t lx, Double_t ly, Double_t lz, Double_t tx, Double_t ty, Double_t tz)
{
   TGeoCtub *ctub = new TGeoCtub(name, rmin, rmax, dz, phi1, phi2, lx, ly, lz, tx, ty, tz);
   TGeoVolume *vol = new TGeoVolume(name, ctub, medium);
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a cone shape with given medium.

TGeoVolume *TGeoBuilder::MakeCone(const char *name, TGeoMedium *medium,
                                     Double_t dz, Double_t rmin1, Double_t rmax1,
                                     Double_t rmin2, Double_t rmax2)
{
   TGeoCone *cone = new TGeoCone(name, dz, rmin1, rmax1, rmin2, rmax2);
   TGeoVolume *vol = 0;
   if (cone->IsRunTimeShape()) {
      vol = fGeometry->MakeVolumeMulti(name, medium);
      vol->SetShape(cone);
   } else {
      vol = new TGeoVolume(name, cone, medium);
   }
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a cone segment shape with given medium

TGeoVolume *TGeoBuilder::MakeCons(const char *name, TGeoMedium *medium,
                                     Double_t dz, Double_t rmin1, Double_t rmax1,
                                     Double_t rmin2, Double_t rmax2,
                                     Double_t phi1, Double_t phi2)
{
   TGeoConeSeg *cons = new TGeoConeSeg(name, dz, rmin1, rmax1, rmin2, rmax2, phi1, phi2);
   TGeoVolume *vol = 0;
   if (cons->IsRunTimeShape()) {
      vol = fGeometry->MakeVolumeMulti(name, medium);
      vol->SetShape(cons);
   } else {
      vol = new TGeoVolume(name, cons, medium);
   }
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a polycone shape with given medium.

TGeoVolume *TGeoBuilder::MakePcon(const char *name, TGeoMedium *medium,
                                     Double_t phi, Double_t dphi, Int_t nz)
{
   TGeoPcon *pcon = new TGeoPcon(name, phi, dphi, nz);
   TGeoVolume *vol = new TGeoVolume(name, pcon, medium);
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a polygone shape with given medium.

TGeoVolume *TGeoBuilder::MakePgon(const char *name, TGeoMedium *medium,
                                     Double_t phi, Double_t dphi, Int_t nedges, Int_t nz)
{
   TGeoPgon *pgon = new TGeoPgon(name, phi, dphi, nedges, nz);
   TGeoVolume *vol = new TGeoVolume(name, pgon, medium);
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a TGeoTrd1 shape with given medium.

TGeoVolume *TGeoBuilder::MakeTrd1(const char *name, TGeoMedium *medium,
                                  Double_t dx1, Double_t dx2, Double_t dy, Double_t dz)
{
   TGeoTrd1 *trd1 = new TGeoTrd1(name, dx1, dx2, dy, dz);
   TGeoVolume *vol = 0;
   if (trd1->IsRunTimeShape()) {
      vol = fGeometry->MakeVolumeMulti(name, medium);
      vol->SetShape(trd1);
   } else {
      vol = new TGeoVolume(name, trd1, medium);
   }
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a TGeoTrd2 shape with given medium.

TGeoVolume *TGeoBuilder::MakeTrd2(const char *name, TGeoMedium *medium,
                                  Double_t dx1, Double_t dx2, Double_t dy1, Double_t dy2,
                                  Double_t dz)
{
   TGeoTrd2 *trd2 = new TGeoTrd2(name, dx1, dx2, dy1, dy2, dz);
   TGeoVolume *vol = 0;
   if (trd2->IsRunTimeShape()) {
      vol = fGeometry->MakeVolumeMulti(name, medium);
      vol->SetShape(trd2);
   } else {
      vol = new TGeoVolume(name, trd2, medium);
   }
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a trapezoid shape with given medium.

TGeoVolume *TGeoBuilder::MakeTrap(const char *name, TGeoMedium *medium,
                                  Double_t dz, Double_t theta, Double_t phi, Double_t h1,
                                  Double_t bl1, Double_t tl1, Double_t alpha1, Double_t h2, Double_t bl2,
                                  Double_t tl2, Double_t alpha2)
{
   TGeoTrap *trap = new TGeoTrap(name, dz, theta, phi, h1, bl1, tl1, alpha1, h2, bl2,
                                 tl2, alpha2);
   TGeoVolume *vol = new TGeoVolume(name, trap, medium);
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make in one step a volume pointing to a twisted trapezoid shape with given medium.

TGeoVolume *TGeoBuilder::MakeGtra(const char *name, TGeoMedium *medium,
                                  Double_t dz, Double_t theta, Double_t phi, Double_t twist, Double_t h1,
                                  Double_t bl1, Double_t tl1, Double_t alpha1, Double_t h2, Double_t bl2,
                                  Double_t tl2, Double_t alpha2)
{
   TGeoGtra *gtra = new TGeoGtra(name, dz, theta, phi, twist, h1, bl1, tl1, alpha1, h2, bl2,
                                 tl2, alpha2);
   TGeoVolume *vol = new TGeoVolume(name, gtra, medium);
   return vol;
}
////////////////////////////////////////////////////////////////////////////////
/// Make a TGeoXtru-shaped volume with nz planes

TGeoVolume *TGeoBuilder::MakeXtru(const char *name, TGeoMedium *medium, Int_t nz)
{
   TGeoXtru *xtru = new TGeoXtru(nz);
   xtru->SetName(name);
   TGeoVolume *vol = new TGeoVolume(name, xtru, medium);
   return vol;
}

////////////////////////////////////////////////////////////////////////////////
/// Make an assembly of volumes.

TGeoVolumeAssembly *TGeoBuilder::MakeVolumeAssembly(const char *name)
{
   return (new TGeoVolumeAssembly(name));
}

////////////////////////////////////////////////////////////////////////////////
/// Make a TGeoVolumeMulti handling a list of volumes.

TGeoVolumeMulti *TGeoBuilder::MakeVolumeMulti(const char *name, TGeoMedium *medium)
{
   return (new TGeoVolumeMulti(name, medium));
}


////////////////////////////////////////////////////////////////////////////////
/// Create a new volume by dividing an existing one (GEANT3 like)
///
/// Divides MOTHER into NDIV divisions called NAME
/// along axis IAXIS starting at coordinate value START
/// and having size STEP. The created volumes will have tracking
/// media ID=NUMED (if NUMED=0 -> same media as MOTHER)
///    The behavior of the division operation can be triggered using OPTION (case insensitive):
///
///  - N  - divide all range in NDIV cells (same effect as STEP<=0) (GSDVN in G3)
///  - NX - divide range starting with START in NDIV cells          (GSDVN2 in G3)
///  - S  - divide all range with given STEP. NDIV is computed and divisions will be centered
///        in full range (same effect as NDIV<=0)                (GSDVS, GSDVT in G3)
///  - SX - same as DVS, but from START position.                   (GSDVS2, GSDVT2 in G3)

TGeoVolume *TGeoBuilder::Division(const char *name, const char *mother, Int_t iaxis,
                                  Int_t ndiv, Double_t start, Double_t step, Int_t numed, Option_t *option)
{
   TGeoVolume *amother;
   TString sname = name;
   sname = sname.Strip();
   const char *vname = sname.Data();
   TString smname = mother;
   smname = smname.Strip();
   const char *mname = smname.Data();

   amother = (TGeoVolume*)fGeometry->GetListOfGVolumes()->FindObject(mname);
   if (!amother) amother = fGeometry->GetVolume(mname);
   if (amother) return amother->Divide(vname,iaxis,ndiv,start,step,numed, option);

   Error("Division","VOLUME: \"%s\" not defined",mname);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create rotation matrix named 'mat<index>'.
///
///  - index    rotation matrix number
///  - theta1   polar angle for axis X
///  - phi1     azimuthal angle for axis X
///  - theta2   polar angle for axis Y
///  - phi2     azimuthal angle for axis Y
///  - theta3   polar angle for axis Z
///  - phi3     azimuthal angle for axis Z
///

void TGeoBuilder::Matrix(Int_t index, Double_t theta1, Double_t phi1,
                         Double_t theta2, Double_t phi2,
                         Double_t theta3, Double_t phi3)
{
   TGeoRotation * rot = new TGeoRotation("",theta1,phi1,theta2,phi2,theta3,phi3);
   rot->SetUniqueID(index);
   rot->RegisterYourself();
}

////////////////////////////////////////////////////////////////////////////////
/// Create material with given A, Z and density, having an unique id.

TGeoMaterial *TGeoBuilder::Material(const char *name, Double_t a, Double_t z, Double_t dens, Int_t uid,Double_t radlen, Double_t intlen)
{
   TGeoMaterial *material = new TGeoMaterial(name,a,z,dens,radlen,intlen);
   material->SetUniqueID(uid);
   return material;
}

////////////////////////////////////////////////////////////////////////////////
/// Create mixture OR COMPOUND IMAT as composed by THE BASIC nelem
/// materials defined by arrays A,Z and WMAT, having an unique id.

TGeoMaterial *TGeoBuilder::Mixture(const char *name, Float_t *a, Float_t *z, Double_t dens,
                                   Int_t nelem, Float_t *wmat, Int_t uid)
{
   TGeoMixture *mix = new TGeoMixture(name,nelem,dens);
   mix->SetUniqueID(uid);
   Int_t i;
   for (i=0;i<nelem;i++) {
      mix->DefineElement(i,a[i],z[i],wmat[i]);
   }
   return (TGeoMaterial*)mix;
}

////////////////////////////////////////////////////////////////////////////////
/// Create mixture OR COMPOUND IMAT as composed by THE BASIC nelem
/// materials defined by arrays A,Z and WMAT, having an unique id.

TGeoMaterial *TGeoBuilder::Mixture(const char *name, Double_t *a, Double_t *z, Double_t dens,
                                   Int_t nelem, Double_t *wmat, Int_t uid)
{
   TGeoMixture *mix = new TGeoMixture(name,nelem,dens);
   mix->SetUniqueID(uid);
   Int_t i;
   for (i=0;i<nelem;i++) {
      mix->DefineElement(i,a[i],z[i],wmat[i]);
   }
   return (TGeoMaterial*)mix;
}

////////////////////////////////////////////////////////////////////////////////
/// Create tracking medium
///
///  - numed      tracking medium number assigned
///  - name      tracking medium name
///  - nmat      material number
///  - isvol     sensitive volume flag
///  - ifield    magnetic field
///  - fieldm    max. field value (kilogauss)
///  - tmaxfd    max. angle due to field (deg/step)
///  - stemax    max. step allowed
///  - deemax    max. fraction of energy lost in a step
///  - epsil     tracking precision (cm)
///  - stmin     min. step due to continuous processes (cm)
///
///  - ifield = 0 if no magnetic field; ifield = -1 if user decision in guswim;
///  - ifield = 1 if tracking performed with g3rkuta; ifield = 2 if tracking
///    performed with g3helix; ifield = 3 if tracking performed with g3helx3.
///

TGeoMedium *TGeoBuilder::Medium(const char *name, Int_t numed, Int_t nmat, Int_t isvol,
                                Int_t ifield, Double_t fieldm, Double_t tmaxfd,
                                Double_t stemax, Double_t deemax, Double_t epsil,
                                Double_t stmin)
{
   return new TGeoMedium(name,numed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a node called <name_nr> pointing to the volume called <name>
/// as daughter of the volume called <mother> (gspos). The relative matrix is
/// made of : a translation (x,y,z) and a rotation matrix named <matIROT>.
/// In case npar>0, create the volume to be positioned in mother, according
/// its actual parameters (gsposp).
///  - NAME   Volume name
///  - NUMBER Copy number of the volume
///  - MOTHER Mother volume name
///  - X      X coord. of the volume in mother ref. sys.
///  - Y      Y coord. of the volume in mother ref. sys.
///  - Z      Z coord. of the volume in mother ref. sys.
///  - IROT   Rotation matrix number w.r.t. mother ref. sys.
///  - ISONLY ONLY/MANY flag

void TGeoBuilder::Node(const char *name, Int_t nr, const char *mother,
                       Double_t x, Double_t y, Double_t z, Int_t irot,
                       Bool_t isOnly, Float_t *upar, Int_t npar)
{
   TGeoVolume *amother= 0;
   TGeoVolume *volume = 0;

   // look into special volume list first
   amother = fGeometry->FindVolumeFast(mother,kTRUE);
   if (!amother) amother = fGeometry->FindVolumeFast(mother);
   if (!amother) {
      TString mname = mother;
      mname = mname.Strip();
      Error("Node","Mother VOLUME \"%s\" not defined",mname.Data());
      return;
   }
   Int_t i;
   if (npar<=0) {
   //---> acting as G3 gspos
      if (gDebug > 0) Info("Node","Calling gspos, mother=%s, name=%s, nr=%d, x=%g, y=%g, z=%g, irot=%d, konly=%i",mother,name,nr,x,y,z,irot,(Int_t)isOnly);
      // look into special volume list first
      volume  = fGeometry->FindVolumeFast(name,kTRUE);
      if (!volume) volume = fGeometry->FindVolumeFast(name);
      if (!volume) {
         TString vname = name;
         vname = vname.Strip();
         Error("Node","VOLUME: \"%s\" not defined",vname.Data());
         return;
      }
      if (((TObject*)volume)->TestBit(TGeoVolume::kVolumeMulti) && !volume->GetShape()) {
         Error("Node", "cannot add multiple-volume object %s as node", volume->GetName());
         return;
      }
   } else {
   //---> acting as G3 gsposp
      TGeoVolumeMulti *vmulti  = (TGeoVolumeMulti*)fGeometry->FindVolumeFast(name, kTRUE);
      if (!vmulti) {
         volume = fGeometry->FindVolumeFast(name);
         if (volume) {
            Warning("Node", "volume: %s is defined as single -> ignoring shape parameters", volume->GetName());
            Node(name,nr,mother,x,y,z,irot,isOnly, upar);
            return;
         }
         TString vname = name;
         vname = vname.Strip();
         Error("Node","VOLUME: \"%s\" not defined ",vname.Data());
         return;
      }
      TGeoMedium *medium = vmulti->GetMedium();
      TString sh    = vmulti->GetTitle();
      sh.ToLower();
      if (sh.Contains("box")) {
         volume = MakeBox(name,medium,upar[0],upar[1],upar[2]);
      } else if (sh.Contains("trd1")) {
         volume = MakeTrd1(name,medium,upar[0],upar[1],upar[2],upar[3]);
      } else if (sh.Contains("trd2")) {
         volume = MakeTrd2(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
      } else if (sh.Contains("trap")) {
         volume = MakeTrap(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6],upar[7],upar[8],upar[9],upar[10]);
      } else if (sh.Contains("gtra")) {
         volume = MakeGtra(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6],upar[7],upar[8],upar[9],upar[10],upar[11]);
      } else if (sh.Contains("tube")) {
         volume = MakeTube(name,medium,upar[0],upar[1],upar[2]);
      } else if (sh.Contains("tubs")) {
         volume = MakeTubs(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
      } else if (sh.Contains("cone")) {
         volume = MakeCone(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
      } else if (sh.Contains("cons")) {
         volume = MakeCons(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6]);
      } else if (sh.Contains("pgon")) {
         volume = MakePgon(name,medium,upar[0],upar[1],(Int_t)upar[2],(Int_t)upar[3]);
         Int_t nz = (Int_t)upar[3];
         for (i=0;i<nz;i++) {
            ((TGeoPgon*)volume->GetShape())->DefineSection(i,upar[3*i+4],upar[3*i+5],upar[3*i+6]);
         }
      } else if (sh.Contains("pcon")) {
         volume = MakePcon(name,medium,upar[0],upar[1],(Int_t)upar[2]);
         Int_t nz = (Int_t)upar[2];
         for (i=0;i<nz;i++) {
            ((TGeoPcon*)volume->GetShape())->DefineSection(i,upar[3*i+3],upar[3*i+4],upar[3*i+5]);
         }
      } else if (sh.Contains("eltu")) {
         volume = MakeEltu(name,medium,upar[0],upar[1],upar[2]);
      } else if (sh.Contains("sphe")) {
         volume = MakeSphere(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5]);
      } else if (sh.Contains("ctub")) {
         volume = MakeCtub(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6],upar[7],upar[8],upar[9],upar[10]);
      } else if (sh.Contains("para")) {
         volume = MakePara(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5]);
      } else {
         Error("Node","cannot create shape %s",sh.Data());
      }

      if (!volume) return;
      vmulti->AddVolume(volume);
   }
   if (irot) {
      TGeoRotation *matrix = 0;
      TGeoMatrix *mat;
      TIter next(fGeometry->GetListOfMatrices());
      while ((mat=(TGeoMatrix*)next())) {
         if (mat->GetUniqueID()==UInt_t(irot)) {
            matrix = dynamic_cast<TGeoRotation*>(mat);
            break;
         }
      }
      if (!matrix) {
         Fatal("Node", "Node %s/%s_%d rotation %i not found",mother, name, nr ,irot);
         return;
      }
      if (isOnly) amother->AddNode(volume,nr,new TGeoCombiTrans(x,y,z,matrix));
      else        amother->AddNodeOverlap(volume,nr,new TGeoCombiTrans(x,y,z,matrix));
   } else {
      if (TMath::Abs(x)<TGeoShape::Tolerance() &&
          TMath::Abs(y)<TGeoShape::Tolerance() &&
          TMath::Abs(z)<TGeoShape::Tolerance()) {
         if (isOnly) amother->AddNode(volume,nr);
         else        amother->AddNodeOverlap(volume,nr);
      } else {
         if (isOnly) amother->AddNode(volume,nr,new TGeoTranslation(x,y,z));
         else        amother->AddNodeOverlap(volume,nr,new TGeoTranslation(x,y,z));
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a node called <name_nr> pointing to the volume called <name>
/// as daughter of the volume called <mother> (gspos). The relative matrix is
/// made of : a translation (x,y,z) and a rotation matrix named <matIROT>.
/// In case npar>0, create the volume to be positioned in mother, according
/// its actual parameters (gsposp).
///  - NAME   Volume name
///  - NUMBER Copy number of the volume
///  - MOTHER Mother volume name
///  - X      X coord. of the volume in mother ref. sys.
///  - Y      Y coord. of the volume in mother ref. sys.
///  - Z      Z coord. of the volume in mother ref. sys.
///  - IROT   Rotation matrix number w.r.t. mother ref. sys.
///  - ISONLY ONLY/MANY flag

void TGeoBuilder::Node(const char *name, Int_t nr, const char *mother,
                       Double_t x, Double_t y, Double_t z, Int_t irot,
                       Bool_t isOnly, Double_t *upar, Int_t npar)
{
   TGeoVolume *amother= 0;
   TGeoVolume *volume = 0;

   // look into special volume list first
   amother = fGeometry->FindVolumeFast(mother,kTRUE);
   if (!amother) amother = fGeometry->FindVolumeFast(mother);
   if (!amother) {
      TString mname = mother;
      mname = mname.Strip();
      Error("Node","Mother VOLUME \"%s\" not defined",mname.Data());
      return;
   }
   Int_t i;
   if (npar<=0) {
   //---> acting as G3 gspos
      if (gDebug > 0) Info("Node","Calling gspos, mother=%s, name=%s, nr=%d, x=%g, y=%g, z=%g, irot=%d, konly=%i",mother,name,nr,x,y,z,irot,(Int_t)isOnly);
      // look into special volume list first
      volume  = fGeometry->FindVolumeFast(name,kTRUE);
      if (!volume) volume = fGeometry->FindVolumeFast(name);
      if (!volume) {
         TString vname = name;
         vname = vname.Strip();
         Error("Node","VOLUME: \"%s\" not defined",vname.Data());
         return;
      }
      if (((TObject*)volume)->TestBit(TGeoVolume::kVolumeMulti) && !volume->GetShape()) {
         Error("Node", "cannot add multiple-volume object %s as node", volume->GetName());
         return;
      }
   } else {
   //---> acting as G3 gsposp
      TGeoVolumeMulti *vmulti  = (TGeoVolumeMulti*)fGeometry->FindVolumeFast(name, kTRUE);
      if (!vmulti) {
         volume = fGeometry->FindVolumeFast(name);
         if (volume) {
            Warning("Node", "volume: %s is defined as single -> ignoring shape parameters", volume->GetName());
            Node(name,nr,mother,x,y,z,irot,isOnly, upar);
            return;
         }
         TString vname = name;
         vname = vname.Strip();
         Error("Node","VOLUME: \"%s\" not defined ",vname.Data());
         return;
      }
      TGeoMedium *medium = vmulti->GetMedium();
      TString sh    = vmulti->GetTitle();
      sh.ToLower();
      if (sh.Contains("box")) {
         volume = MakeBox(name,medium,upar[0],upar[1],upar[2]);
      } else if (sh.Contains("trd1")) {
         volume = MakeTrd1(name,medium,upar[0],upar[1],upar[2],upar[3]);
      } else if (sh.Contains("trd2")) {
         volume = MakeTrd2(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
      } else if (sh.Contains("trap")) {
         volume = MakeTrap(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6],upar[7],upar[8],upar[9],upar[10]);
      } else if (sh.Contains("gtra")) {
         volume = MakeGtra(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6],upar[7],upar[8],upar[9],upar[10],upar[11]);
      } else if (sh.Contains("tube")) {
         volume = MakeTube(name,medium,upar[0],upar[1],upar[2]);
      } else if (sh.Contains("tubs")) {
         volume = MakeTubs(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
      } else if (sh.Contains("cone")) {
         volume = MakeCone(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
      } else if (sh.Contains("cons")) {
         volume = MakeCons(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6]);
      } else if (sh.Contains("pgon")) {
         volume = MakePgon(name,medium,upar[0],upar[1],(Int_t)upar[2],(Int_t)upar[3]);
         Int_t nz = (Int_t)upar[3];
         for (i=0;i<nz;i++) {
            ((TGeoPgon*)volume->GetShape())->DefineSection(i,upar[3*i+4],upar[3*i+5],upar[3*i+6]);
         }
      } else if (sh.Contains("pcon")) {
         volume = MakePcon(name,medium,upar[0],upar[1],(Int_t)upar[2]);
         Int_t nz = (Int_t)upar[2];
         for (i=0;i<nz;i++) {
            ((TGeoPcon*)volume->GetShape())->DefineSection(i,upar[3*i+3],upar[3*i+4],upar[3*i+5]);
         }
      } else if (sh.Contains("eltu")) {
         volume = MakeEltu(name,medium,upar[0],upar[1],upar[2]);
      } else if (sh.Contains("sphe")) {
         volume = MakeSphere(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5]);
      } else if (sh.Contains("ctub")) {
         volume = MakeCtub(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6],upar[7],upar[8],upar[9],upar[10]);
      } else if (sh.Contains("para")) {
         volume = MakePara(name,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5]);
      } else {
         Error("Node","cannot create shape %s",sh.Data());
      }

      if (!volume) return;
      vmulti->AddVolume(volume);
   }
   if (irot) {
      TGeoRotation *matrix = 0;
      TGeoMatrix *mat;
      TIter next(fGeometry->GetListOfMatrices());
      while ((mat=(TGeoMatrix*)next())) {
         if (mat->GetUniqueID()==UInt_t(irot)) {
            matrix = dynamic_cast<TGeoRotation*>(mat);
            break;
         }
      }
      if (!matrix) {
         Fatal("Node", "Node %s/%s_%d rotation %i not found",mother, name, nr ,irot);
         return;
      }
      if (isOnly) amother->AddNode(volume,nr,new TGeoCombiTrans(x,y,z,matrix));
      else        amother->AddNodeOverlap(volume,nr,new TGeoCombiTrans(x,y,z,matrix));
   } else {
      if (TMath::Abs(x)<TGeoShape::Tolerance() &&
          TMath::Abs(y)<TGeoShape::Tolerance() &&
          TMath::Abs(z)<TGeoShape::Tolerance()) {
         if (isOnly) amother->AddNode(volume,nr);
         else        amother->AddNodeOverlap(volume,nr);
      } else {
         if (isOnly) amother->AddNode(volume,nr,new TGeoTranslation(x,y,z));
         else        amother->AddNodeOverlap(volume,nr,new TGeoTranslation(x,y,z));
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a volume in GEANT3 style.
///  - NAME   Volume name
///  - SHAPE  Volume type
///  - NMED   Tracking medium number
///  - NPAR   Number of shape parameters
///  - UPAR   Vector containing shape parameters

TGeoVolume *TGeoBuilder::Volume(const char *name, const char *shape, Int_t nmed,
                                Float_t *upar, Int_t npar)
{
   Int_t i;
   TGeoVolume *volume = 0;
   TGeoMedium *medium = fGeometry->GetMedium(nmed);
   if (!medium) {
      Error("Volume","cannot create volume: %s, medium: %d is unknown",name,nmed);
      return 0;
   }
   TString sh = shape;
   TString sname = name;
   sname = sname.Strip();
   const char *vname = sname.Data();
   if (npar <= 0) {
      //--- create a TGeoVolumeMulti
      volume = MakeVolumeMulti(vname,medium);
      volume->SetTitle(shape);
      TGeoVolumeMulti *vmulti  = (TGeoVolumeMulti*)fGeometry->GetListOfGVolumes()->FindObject(vname);
      if (!vmulti) {
         Error("Volume","volume multi: %s not created",vname);
         return 0;
      }
      return vmulti;
   }
   //---> create a normal volume
   sh.ToLower();
   if (sh.Contains("box")) {
      volume = MakeBox(vname,medium,upar[0],upar[1],upar[2]);
   } else if (sh.Contains("trd1")) {
      volume = MakeTrd1(vname,medium,upar[0],upar[1],upar[2],upar[3]);
   } else if (sh.Contains("trd2")) {
      volume = MakeTrd2(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
   } else if (sh.Contains("trap")) {
      volume = MakeTrap(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6],upar[7],upar[8],upar[9],upar[10]);
   } else if (sh.Contains("gtra")) {
      volume = MakeGtra(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6],upar[7],upar[8],upar[9],upar[10],upar[11]);
   } else if (sh.Contains("tube")) {
      volume = MakeTube(vname,medium,upar[0],upar[1],upar[2]);
   } else if (sh.Contains("tubs")) {
      volume = MakeTubs(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
   } else if (sh.Contains("cone")) {
      volume = MakeCone(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
   } else if (sh.Contains("cons")) {
      volume = MakeCons(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6]);
   } else if (sh.Contains("pgon")) {
      volume = MakePgon(vname,medium,upar[0],upar[1],(Int_t)upar[2],(Int_t)upar[3]);
      Int_t nz = (Int_t)upar[3];
      for (i=0;i<nz;i++) {
         ((TGeoPgon*)volume->GetShape())->DefineSection(i,upar[3*i+4],upar[3*i+5],upar[3*i+6]);
      }
   } else if (sh.Contains("pcon")) {
      volume = MakePcon(vname,medium,upar[0],upar[1],(Int_t)upar[2]);
      Int_t nz = (Int_t)upar[2];
      for (i=0;i<nz;i++) {
         ((TGeoPcon*)volume->GetShape())->DefineSection(i,upar[3*i+3],upar[3*i+4],upar[3*i+5]);
      }
   } else if (sh.Contains("eltu")) {
      volume = MakeEltu(vname,medium,upar[0],upar[1],upar[2]);
   } else if (sh.Contains("sphe")) {
      volume = MakeSphere(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5]);
   } else if (sh.Contains("ctub")) {
      volume = MakeCtub(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6],upar[7],upar[8],upar[9],upar[10]);
   } else if (sh.Contains("para")) {
      volume = MakePara(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5]);
   } else if (sh.Contains("tor")) {
      volume = MakeTorus(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
   }

   if (!volume) {
      Error("Volume","volume: %s not created",vname);
      return 0;
   }
   return volume;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a volume in GEANT3 style.
///  - NAME   Volume name
///  - SHAPE  Volume type
///  - NMED   Tracking medium number
///  - NPAR   Number of shape parameters
///  - UPAR   Vector containing shape parameters

TGeoVolume *TGeoBuilder::Volume(const char *name, const char *shape, Int_t nmed,
                                Double_t *upar, Int_t npar)
{
   Int_t i;
   TGeoVolume *volume = 0;
   TGeoMedium *medium = fGeometry->GetMedium(nmed);
   if (!medium) {
      Error("Volume","cannot create volume: %s, medium: %d is unknown",name,nmed);
      return 0;
   }
   TString sh = shape;
   TString sname = name;
   sname = sname.Strip();
   const char *vname = sname.Data();
   if (npar <= 0) {
      //--- create a TGeoVolumeMulti
      volume = MakeVolumeMulti(vname,medium);
      volume->SetTitle(shape);
      TGeoVolumeMulti *vmulti  = (TGeoVolumeMulti*)fGeometry->GetListOfGVolumes()->FindObject(vname);
      if (!vmulti) {
         Error("Volume","volume multi: %s not created",vname);
         return 0;
      }
      return vmulti;
   }
   //---> create a normal volume
   sh.ToLower();
   if (sh.Contains("box")) {
      volume = MakeBox(vname,medium,upar[0],upar[1],upar[2]);
   } else if (sh.Contains("trd1")) {
      volume = MakeTrd1(vname,medium,upar[0],upar[1],upar[2],upar[3]);
   } else if (sh.Contains("trd2")) {
      volume = MakeTrd2(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
   } else if (sh.Contains("trap")) {
      volume = MakeTrap(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6],upar[7],upar[8],upar[9],upar[10]);
   } else if (sh.Contains("gtra")) {
      volume = MakeGtra(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6],upar[7],upar[8],upar[9],upar[10],upar[11]);
   } else if (sh.Contains("tube")) {
      volume = MakeTube(vname,medium,upar[0],upar[1],upar[2]);
   } else if (sh.Contains("tubs")) {
      volume = MakeTubs(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
   } else if (sh.Contains("cone")) {
      volume = MakeCone(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
   } else if (sh.Contains("cons")) {
      volume = MakeCons(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6]);
   } else if (sh.Contains("pgon")) {
      volume = MakePgon(vname,medium,upar[0],upar[1],(Int_t)upar[2],(Int_t)upar[3]);
      Int_t nz = (Int_t)upar[3];
      for (i=0;i<nz;i++) {
         ((TGeoPgon*)volume->GetShape())->DefineSection(i,upar[3*i+4],upar[3*i+5],upar[3*i+6]);
      }
   } else if (sh.Contains("pcon")) {
      volume = MakePcon(vname,medium,upar[0],upar[1],(Int_t)upar[2]);
      Int_t nz = (Int_t)upar[2];
      for (i=0;i<nz;i++) {
         ((TGeoPcon*)volume->GetShape())->DefineSection(i,upar[3*i+3],upar[3*i+4],upar[3*i+5]);
      }
   } else if (sh.Contains("eltu")) {
      volume = MakeEltu(vname,medium,upar[0],upar[1],upar[2]);
   } else if (sh.Contains("sphe")) {
      volume = MakeSphere(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5]);
   } else if (sh.Contains("ctub")) {
      volume = MakeCtub(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5],upar[6],upar[7],upar[8],upar[9],upar[10]);
   } else if (sh.Contains("para")) {
      volume = MakePara(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4],upar[5]);
   } else if (sh.Contains("tor")) {
      volume = MakeTorus(vname,medium,upar[0],upar[1],upar[2],upar[3],upar[4]);
   }

   if (!volume) {
      Error("Volume","volume: %s not created",vname);
      return 0;
   }
   return volume;
}

