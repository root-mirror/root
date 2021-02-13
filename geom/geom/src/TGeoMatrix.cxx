// @(#)root/geom:$Id$
// Author: Andrei Gheata   25/10/01

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGeoMatrix
\ingroup Geometry_classes

Geometrical transformation package.

  All geometrical transformations handled by the modeller are provided as a
built-in package. This was designed to minimize memory requirements and
optimize performance of point/vector master-to-local and local-to-master
computation. We need to have in mind that a transformation in TGeo has 2
major use-cases. The first one is for defining the placement of a volume
with respect to its container reference frame. This frame will be called
'master' and the frame of the positioned volume - 'local'. If T is a
transformation used for positioning volume daughters, then:

~~~ {.cpp}
         MASTER = T * LOCAL
~~~

  Therefore a local-to-master conversion will be performed by using T, while
a master-to-local by using its inverse. The second use case is the computation
of the global transformation of a given object in the geometry. Since the
geometry is built as 'volumes-inside-volumes', this global transformation
represent the pile-up of all local transformations in the corresponding
branch. The conversion from the global reference frame and the given object
is also called master-to-local, but it is handled by the manager class.
  A general homogenous transformation is defined as a 4x4 matrix embedding
a rotation, a translation and a scale. The advantage of this description
is that each basic transformation can be represented as a homogenous matrix,
composition being performed as simple matrix multiplication.

  Rotation:                      Inverse rotation:

~~~ {.cpp}
        r11  r12  r13   0              r11  r21  r31   0
        r21  r22  r23   0              r12  r22  r32   0
        r31  r32  r33   0              r13  r23  r33   0
         0    0    0    1               0    0    0    1
~~~

  Translation:                   Inverse translation:

~~~ {.cpp}
         1    0    0    tx               1    0    0   -tx
         0    1    0    ty               0    1    0   -ty
         0    0    1    tz               0    0    1   -tz
         0    0    0    1                0    0    0   1
~~~

  Scale:                         Inverse scale:

~~~ {.cpp}
         sx   0    0    0              1/sx  0    0    0
         0    sy   0    0               0   1/sy  0    0
         0    0    sz   0               0    0   1/sz  0
         0    0    0    1               0    0    0    1
~~~

 where:
       - `rij` are the 3x3 rotation matrix components,
       - `tx`, `ty`, `tz` are the translation components
       - `sx`, `sy`, `sz` are arbitrary scale constants on each axis,

  The disadvantage in using this approach is that computation for 4x4 matrices
is expensive. Even combining two translation would become a multiplication
of their corresponding matrices, which is quite an undesired effect. On the
other hand, it is not a good idea to store a translation as a block of 16
numbers. We have therefore chosen to implement each basic transformation type
as a class deriving from the same basic abstract class and handling its specific
data and point/vector transformation algorithms.

\image html geom_transf.jpg

### The base class TGeoMatrix defines abstract metods for:

#### translation, rotation and scale getters. Every derived class stores only
  its specific data, e.g. a translation stores an array of 3 doubles and a
  rotation an array of 9. However, asking which is the rotation array of a
  TGeoTranslation through the base TGeoMatrix interface is a legal operation.
  The answer in this case is a pointer to a global constant array representing
  an identity rotation.

~~~ {.cpp}
     Double_t *TGeoMatrix::GetTranslation()
     Double_t *TGeoMatrix::GetRotation()
     Double_t *TGeoMatrix::GetScale()
~~~

#### MasterToLocal() and LocalToMaster() point and vector transformations :

~~~ {.cpp}
     void      TGeoMatrix::MasterToLocal(const Double_t *master, Double_t *local)
     void      TGeoMatrix::LocalToMaster(const Double_t *local, Double_t *master)
     void      TGeoMatrix::MasterToLocalVect(const Double_t *master, Double_t *local)
     void      TGeoMatrix::LocalToMasterVect(const Double_t *local, Double_t *master)
~~~

  These allow correct conversion also for reflections.

#### Transformation type getters :

~~~ {.cpp}
     Bool_t    TGeoMatrix::IsIdentity()
     Bool_t    TGeoMatrix::IsTranslation()
     Bool_t    TGeoMatrix::IsRotation()
     Bool_t    TGeoMatrix::IsScale()
     Bool_t    TGeoMatrix::IsCombi() (translation + rotation)
     Bool_t    TGeoMatrix::IsGeneral() (translation + rotation + scale)
~~~

  Combinations of basic transformations are represented by specific classes
deriving from TGeoMatrix. In order to define a matrix as a combination of several
others, a special class TGeoHMatrix is provided. Here is an example of matrix
creation :

### Matrix creation example:

~~~ {.cpp}
  root[0] TGeoRotation r1,r2;
          r1.SetAngles(90,0,30);        // rotation defined by Euler angles
          r2.SetAngles(90,90,90,180,0,0); // rotation defined by GEANT3 angles
          TGeoTranslation t1(-10,10,0);
          TGeoTranslation t2(10,-10,5);
          TGeoCombiTrans c1(t1,r1);
          TGeoCombiTrans c2(t2,r2);
          TGeoHMatrix h = c1 * c2; // composition is done via TGeoHMatrix class
  root[7] TGeoHMatrix *ph = new TGeoHMatrix(hm); // this is the one we want to
                                               // use for positioning a volume
  root[8] ph->Print();
          ...
          pVolume->AddNode(pVolDaughter,id,ph) // now ph is owned by the manager
~~~

### Rule for matrix creation:
 Unless explicitly used for positioning nodes (TGeoVolume::AddNode()) all
matrices deletion have to be managed by users. Matrices passed to geometry
have to be created by using new() operator and their deletion is done by
TGeoManager class.

### Available geometrical transformations

#### TGeoTranslation
Represent a (dx,dy,dz) translation. Data members:
   Double_t fTranslation[3]. Translations can be added/subtracted.

~~~ {.cpp}
   TGeoTranslation t1;
   t1->SetTranslation(-5,10,4);
   TGeoTranslation *t2 = new TGeoTranslation(4,3,10);
   t2->Subtract(&t1);
~~~

#### Rotations
 Represent a pure rotation. Data members: Double_t fRotationMatrix[3*3].
   Rotations can be defined either by Euler angles, either, by GEANT3 angles :

~~~ {.cpp}
   TGeoRotation *r1 = new TGeoRotation();
   r1->SetAngles(phi, theta, psi); // all angles in degrees
~~~

   This represent the composition of : first a rotation about Z axis with
   angle phi, then a rotation with theta about the rotated X axis, and
   finally a rotation with psi about the new Z axis.

~~~ {.cpp}
   r1->SetAngles(th1,phi1, th2,phi2, th3,phi3)
~~~

   This is a rotation defined in GEANT3 style. Theta and phi are the spherical
   angles of each axis of the rotated coordinate system with respect to the
   initial one. This construction allows definition of malformed rotations,
   e.g. not orthogonal. A check is performed and an error message is issued
   in this case.

   Specific utilities : determinant, inverse.

#### Scale transformations
     Represent a scale shrinking/enlargement. Data
     members :Double_t fScale[3]. Not fully implemented yet.

#### Combined transformations
Represent a rotation followed by a translation.
Data members: Double_t fTranslation[3], TGeoRotation *fRotation.

~~~ {.cpp}
   TGeoRotation *rot = new TGeoRotation("rot",10,20,30);
   TGeoTranslation trans;
   ...
   TGeoCombiTrans *c1 = new TGeoCombiTrans(trans, rot);
   TGeoCombiTrans *c2 = new TGeoCombiTrans("somename",10,20,30,rot)
~~~


#### TGeoGenTrans
Combined transformations including a scale. Not implemented.

#### TGeoIdentity
A generic singleton matrix representing a identity transformation
   NOTE: identified by the global variable gGeoIdentity.
*/

#include <iostream>
#include "TObjArray.h"

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TMath.h"

TGeoIdentity *gGeoIdentity = nullptr;
const Int_t kN3 = 3*sizeof(Double_t);
const Int_t kN9 = 9*sizeof(Double_t);

// statics and globals

ClassImp(TGeoMatrix);

////////////////////////////////////////////////////////////////////////////////
/// dummy constructor

TGeoMatrix::TGeoMatrix()
{
  ResetBit(kGeoMatrixBits);
}

////////////////////////////////////////////////////////////////////////////////
/// copy constructor

TGeoMatrix::TGeoMatrix(const TGeoMatrix &other)
           :TNamed(other)
{
   ResetBit(kGeoRegistered);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TGeoMatrix::TGeoMatrix(const char *name)
           :TNamed(name, "")
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoMatrix::~TGeoMatrix()
{
   if (IsRegistered() && gGeoManager) {
      if (!gGeoManager->IsCleaning()) {
         gGeoManager->GetListOfMatrices()->Remove(this);
         Warning("dtor", "Registered matrix %s was removed", GetName());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if no rotation or the rotation is about Z axis

Bool_t TGeoMatrix::IsRotAboutZ() const
{
   if (IsIdentity()) return kTRUE;
   const Double_t *rot = GetRotationMatrix();
   if (TMath::Abs(rot[6])>1E-9) return kFALSE;
   if (TMath::Abs(rot[7])>1E-9) return kFALSE;
   if ((1.-TMath::Abs(rot[8]))>1E-9) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get total size in bytes of this

Int_t TGeoMatrix::GetByteCount() const
{
   Int_t count = 4+28+strlen(GetName())+strlen(GetTitle()); // fId + TNamed
   if (IsTranslation()) count += 12;
   if (IsScale()) count += 12;
   if (IsCombi() || IsGeneral()) count += 4 + 36;
   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// Provide a pointer name containing uid.

char *TGeoMatrix::GetPointerName() const
{
   static TString name;
   name = TString::Format("pMatrix%d", GetUniqueID());
   return (char*)name.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// The homogenous matrix associated with the transformation is used for
/// piling up's and visualization. A homogenous matrix is a 4*4 array
/// containing the translation, the rotation and the scale components
/// ~~~ {.cpp}
///          | R00*sx  R01    R02    dx |
///          | R10    R11*sy  R12    dy |
///          | R20     R21   R22*sz  dz |
///          |  0       0      0      1 |
/// ~~~
///   where Rij is the rotation matrix, (sx, sy, sz) is the scale
/// transformation and (dx, dy, dz) is the translation.

void TGeoMatrix::GetHomogenousMatrix(Double_t *hmat) const
{
   Double_t *hmatrix = hmat;
   const Double_t *mat = GetRotationMatrix();
   for (Int_t i=0; i<3; i++) {
      memcpy(hmatrix, mat, kN3);
      mat     += 3;
      hmatrix += 3;
      *hmatrix = 0.0;
      hmatrix++;
   }
   memcpy(hmatrix, GetTranslation(), kN3);
   hmatrix = hmat;
   if (IsScale()) {
      for (Int_t i=0; i<3; i++) {
         *hmatrix *= GetScale()[i];
         hmatrix  += 5;
      }
   }
   hmatrix[15] = 1.;
}

////////////////////////////////////////////////////////////////////////////////
/// convert a point by multiplying its column vector (x, y, z, 1) to matrix inverse

void TGeoMatrix::LocalToMaster(const Double_t *local, Double_t *master) const
{
   if (IsIdentity()) {
      memcpy(master, local, kN3);
      return;
   }
   Int_t i;
   const Double_t *tr = GetTranslation();
   if (!IsRotation()) {
      for (i=0; i<3; i++) master[i] = tr[i] + local[i];
      return;
   }
   const Double_t *rot = GetRotationMatrix();
   for (i=0; i<3; i++) {
      master[i] = tr[i]
                + local[0]*rot[3*i]
                + local[1]*rot[3*i+1]
                + local[2]*rot[3*i+2];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// convert a vector by multiplying its column vector (x, y, z, 1) to matrix inverse

void TGeoMatrix::LocalToMasterVect(const Double_t *local, Double_t *master) const
{
   if (!IsRotation()) {
      memcpy(master, local, kN3);
      return;
   }
   const Double_t *rot = GetRotationMatrix();
   for (Int_t i=0; i<3; i++) {
      master[i] = local[0]*rot[3*i]
                + local[1]*rot[3*i+1]
                + local[2]*rot[3*i+2];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// convert a point by multiplying its column vector (x, y, z, 1) to matrix inverse

void TGeoMatrix::LocalToMasterBomb(const Double_t *local, Double_t *master) const
{
   if (IsIdentity()) {
      memcpy(master, local, kN3);
      return;
   }
   Int_t i;
   const Double_t *tr = GetTranslation();
   Double_t bombtr[3] = {0.,0.,0.};
   gGeoManager->BombTranslation(tr, &bombtr[0]);
   if (!IsRotation()) {
      for (i=0; i<3; i++) master[i] = bombtr[i] + local[i];
      return;
   }
   const Double_t *rot = GetRotationMatrix();
   for (i=0; i<3; i++) {
      master[i] = bombtr[i]
                + local[0]*rot[3*i]
                + local[1]*rot[3*i+1]
                + local[2]*rot[3*i+2];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// convert a point by multiplying its column vector (x, y, z, 1) to matrix

void TGeoMatrix::MasterToLocal(const Double_t *master, Double_t *local) const
{
   if (IsIdentity()) {
      memcpy(local, master, kN3);
      return;
   }
   const Double_t *tr  = GetTranslation();
   Double_t mt0  = master[0]-tr[0];
   Double_t mt1  = master[1]-tr[1];
   Double_t mt2  = master[2]-tr[2];
   if (!IsRotation()) {
      local[0] = mt0;
      local[1] = mt1;
      local[2] = mt2;
      return;
   }
   const Double_t *rot = GetRotationMatrix();
   local[0] = mt0*rot[0] + mt1*rot[3] + mt2*rot[6];
   local[1] = mt0*rot[1] + mt1*rot[4] + mt2*rot[7];
   local[2] = mt0*rot[2] + mt1*rot[5] + mt2*rot[8];
}

////////////////////////////////////////////////////////////////////////////////
/// convert a point by multiplying its column vector (x, y, z, 1) to matrix

void TGeoMatrix::MasterToLocalVect(const Double_t *master, Double_t *local) const
{
   if (!IsRotation()) {
      memcpy(local, master, kN3);
      return;
   }
   const Double_t *rot = GetRotationMatrix();
   for (Int_t i=0; i<3; i++) {
      local[i] = master[0]*rot[i]
               + master[1]*rot[i+3]
               + master[2]*rot[i+6];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// convert a point by multiplying its column vector (x, y, z, 1) to matrix

void TGeoMatrix::MasterToLocalBomb(const Double_t *master, Double_t *local) const
{
   if (IsIdentity()) {
      memcpy(local, master, kN3);
      return;
   }
   const Double_t *tr = GetTranslation();
   Double_t bombtr[3] = {0.,0.,0.};
   Int_t i;
   gGeoManager->UnbombTranslation(tr, &bombtr[0]);
   if (!IsRotation()) {
      for (i=0; i<3; i++) local[i] = master[i]-bombtr[i];
      return;
   }
   const Double_t *rot = GetRotationMatrix();
   for (i=0; i<3; i++) {
      local[i] = (master[0]-bombtr[0])*rot[i]
               + (master[1]-bombtr[1])*rot[i+3]
               + (master[2]-bombtr[2])*rot[i+6];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Normalize a vector.

void TGeoMatrix::Normalize(Double_t *vect)
{
   Double_t normfactor = vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2];
   if (normfactor <= 1E-10) return;
   normfactor = 1./TMath::Sqrt(normfactor);
   vect[0] *= normfactor;
   vect[1] *= normfactor;
   vect[2] *= normfactor;
}

////////////////////////////////////////////////////////////////////////////////
/// print the matrix in 4x4 format

void TGeoMatrix::Print(Option_t *) const
{
   const Double_t *rot = GetRotationMatrix();
   const Double_t *tr  = GetTranslation();
   printf("matrix %s - tr=%d  rot=%d  refl=%d  scl=%d shr=%d reg=%d own=%d\n", GetName(),(Int_t)IsTranslation(),
          (Int_t)IsRotation(), (Int_t)IsReflection(), (Int_t)IsScale(), (Int_t)IsShared(), (Int_t)IsRegistered(),
          (Int_t)IsOwned());
   printf("%10.6f%12.6f%12.6f    Tx = %10.6f\n", rot[0], rot[1], rot[2], tr[0]);
   printf("%10.6f%12.6f%12.6f    Ty = %10.6f\n", rot[3], rot[4], rot[5], tr[1]);
   printf("%10.6f%12.6f%12.6f    Tz = %10.6f\n", rot[6], rot[7], rot[8], tr[2]);
   if (IsScale()) {
      const Double_t *scl  = GetScale();
      printf("Sx=%10.6fSy=%12.6fSz=%12.6f\n", scl[0], scl[1], scl[2]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply by a reflection respect to YZ.

void TGeoMatrix::ReflectX(Bool_t, Bool_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply by a reflection respect to ZX.

void TGeoMatrix::ReflectY(Bool_t, Bool_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply by a reflection respect to XY.

void TGeoMatrix::ReflectZ(Bool_t, Bool_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Register the matrix in the current manager, which will become the owner.

void TGeoMatrix::RegisterYourself()
{
   if (!gGeoManager) {
      Warning("RegisterYourself", "cannot register without geometry");
      return;
   }
   if (!IsRegistered()) {
      gGeoManager->RegisterMatrix(this);
      SetBit(kGeoRegistered);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If no name was supplied in the ctor, the type of transformation is checked.
/// A letter will be prepended to the name :
///  - t - translation
///  - r - rotation
///  - s - scale
///  - c - combi (translation + rotation)
///  - g - general (tr+rot+scale)
/// The index of the transformation in gGeoManager list of transformations will
/// be appended.

void TGeoMatrix::SetDefaultName()
{
   if (!gGeoManager) return;
   if (strlen(GetName())) return;
   char type = 'n';
   if (IsTranslation()) type = 't';
   if (IsRotation()) type = 'r';
   if (IsScale()) type = 's';
   if (IsCombi()) type = 'c';
   if (IsGeneral()) type = 'g';
   TObjArray *matrices = gGeoManager->GetListOfMatrices();
   Int_t index = 0;
   if (matrices) index =matrices->GetEntriesFast() - 1;
   TString name = TString::Format("%c%d", type, index);
   SetName(name);
}

/** \class TGeoTranslation
\ingroup Geometry_classes

Class describing translations. A translation is
basically an array of 3 doubles matching the positions 12, 13
and 14 in the homogenous matrix description.
*/

ClassImp(TGeoTranslation);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoTranslation::TGeoTranslation()
{
   for (Int_t i=0; i<3; i++) fTranslation[i] = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor.

TGeoTranslation::TGeoTranslation(const TGeoTranslation &other)
                :TGeoMatrix(other)
{
   SetTranslation(other);
}

////////////////////////////////////////////////////////////////////////////////
/// Ctor. based on a general matrix

TGeoTranslation::TGeoTranslation(const TGeoMatrix &other)
                :TGeoMatrix(other)
{
   ResetBit(kGeoRotation);
   ResetBit(kGeoScale);
   SetTranslation(other);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor defining the translation

TGeoTranslation::TGeoTranslation(Double_t dx, Double_t dy, Double_t dz)
                :TGeoMatrix("")
{
   if (dx || dy || dz) SetBit(kGeoTranslation);
   SetTranslation(dx, dy, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor defining the translation

TGeoTranslation::TGeoTranslation(const char *name, Double_t dx, Double_t dy, Double_t dz)
                :TGeoMatrix(name)
{
   if (dx || dy || dz) SetBit(kGeoTranslation);
   SetTranslation(dx, dy, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment from a general matrix

TGeoTranslation& TGeoTranslation::operator = (const TGeoMatrix &matrix)
{
   if (&matrix == this) return *this;
   Bool_t registered = TestBit(kGeoRegistered);
   TNamed::operator=(matrix);
   SetTranslation(matrix);
   SetBit(kGeoRegistered,registered);
   ResetBit(kGeoRotation);
   ResetBit(kGeoScale);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Translation composition

TGeoTranslation &TGeoTranslation::operator*=(const TGeoTranslation &right)
{
   const Double_t *tr = right.GetTranslation();
   fTranslation[0] += tr[0];
   fTranslation[1] += tr[1];
   fTranslation[2] += tr[2];
   if (!IsTranslation()) SetBit(kGeoTranslation, right.IsTranslation());
   return *this;
}

TGeoTranslation TGeoTranslation::operator*(const TGeoTranslation &right) const
{
   TGeoTranslation t = *this;
   t *= right;
   return t;
}

TGeoHMatrix TGeoTranslation::operator*(const TGeoMatrix &right) const
{
   TGeoHMatrix t = *this;
   t *= right;
   return t;
}

////////////////////////////////////////////////////////////////////////////////
/// Is-equal operator

Bool_t TGeoTranslation::operator ==(const TGeoTranslation &other) const
{
   if (&other == this) return kTRUE;
   const Double_t *tr = GetTranslation();
   const Double_t *otr = other.GetTranslation();
   for (auto i=0; i<3; i++)
      if (TMath::Abs(tr[i]-otr[i])>1.E-10) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a temporary inverse of this.

TGeoHMatrix TGeoTranslation::Inverse() const
{
   TGeoHMatrix h;
   h = *this;
   Double_t tr[3];
   tr[0] = -fTranslation[0];
   tr[1] = -fTranslation[1];
   tr[2] = -fTranslation[2];
   h.SetTranslation(tr);
   return h;
}

////////////////////////////////////////////////////////////////////////////////
/// Adding a translation to this one

void TGeoTranslation::Add(const TGeoTranslation *other)
{
   const Double_t *trans = other->GetTranslation();
   for (Int_t i=0; i<3; i++)
      fTranslation[i] += trans[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Make a clone of this matrix.

TGeoMatrix *TGeoTranslation::MakeClone() const
{
   TGeoMatrix *matrix = new TGeoTranslation(*this);
   return matrix;
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate about X axis of the master frame with angle expressed in degrees.

void TGeoTranslation::RotateX(Double_t /*angle*/)
{
   Warning("RotateX", "Not implemented. Use TGeoCombiTrans instead");
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate about Y axis of the master frame with angle expressed in degrees.

void TGeoTranslation::RotateY(Double_t /*angle*/)
{
   Warning("RotateY", "Not implemented. Use TGeoCombiTrans instead");
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate about Z axis of the master frame with angle expressed in degrees.

void TGeoTranslation::RotateZ(Double_t /*angle*/)
{
   Warning("RotateZ", "Not implemented. Use TGeoCombiTrans instead");
}

////////////////////////////////////////////////////////////////////////////////
/// Subtracting a translation from this one

void TGeoTranslation::Subtract(const TGeoTranslation *other)
{
   const Double_t *trans = other->GetTranslation();
   for (Int_t i=0; i<3; i++)
      fTranslation[i] -= trans[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Set translation components

void TGeoTranslation::SetTranslation(Double_t dx, Double_t dy, Double_t dz)
{
   fTranslation[0] = dx;
   fTranslation[1] = dy;
   fTranslation[2] = dz;
   if (dx || dy || dz) SetBit(kGeoTranslation);
   else                ResetBit(kGeoTranslation);
}

////////////////////////////////////////////////////////////////////////////////
/// Set translation components

void TGeoTranslation::SetTranslation(const TGeoMatrix &other)
{
   SetBit(kGeoTranslation, other.IsTranslation());
   const Double_t *transl = other.GetTranslation();
   memcpy(fTranslation, transl, kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// convert a point by multiplying its column vector (x, y, z, 1) to matrix inverse

void TGeoTranslation::LocalToMaster(const Double_t *local, Double_t *master) const
{
   const Double_t *tr = GetTranslation();
   for (Int_t i=0; i<3; i++)
      master[i] = tr[i] + local[i];
}

////////////////////////////////////////////////////////////////////////////////
/// convert a vector to MARS

void TGeoTranslation::LocalToMasterVect(const Double_t *local, Double_t *master) const
{
   memcpy(master, local, kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// convert a point by multiplying its column vector (x, y, z, 1) to matrix inverse

void TGeoTranslation::LocalToMasterBomb(const Double_t *local, Double_t *master) const
{
   const Double_t *tr = GetTranslation();
   Double_t bombtr[3] = {0.,0.,0.};
   gGeoManager->BombTranslation(tr, &bombtr[0]);
   for (Int_t i=0; i<3; i++)
      master[i] = bombtr[i] + local[i];
}

////////////////////////////////////////////////////////////////////////////////
/// convert a point by multiplying its column vector (x, y, z, 1) to matrix

void TGeoTranslation::MasterToLocal(const Double_t *master, Double_t *local) const
{
   const Double_t *tr = GetTranslation();
   for (Int_t i=0; i<3; i++)
      local[i] =  master[i]-tr[i];
}

////////////////////////////////////////////////////////////////////////////////
/// convert a vector from MARS to local

void TGeoTranslation::MasterToLocalVect(const Double_t *master, Double_t *local) const
{
   memcpy(local, master, kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// convert a point by multiplying its column vector (x, y, z, 1) to matrix

void TGeoTranslation::MasterToLocalBomb(const Double_t *master, Double_t *local) const
{
   const Double_t *tr = GetTranslation();
   Double_t bombtr[3] = {0.,0.,0.};
   gGeoManager->UnbombTranslation(tr, &bombtr[0]);
   for (Int_t i=0; i<3; i++)
      local[i] =  master[i]-bombtr[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoTranslation::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TestBit(kGeoSavePrimitive)) return;
   out << "   // Translation: " << GetName() << std::endl;
   out << "   dx = " << fTranslation[0] << ";" << std::endl;
   out << "   dy = " << fTranslation[1] << ";" << std::endl;
   out << "   dz = " << fTranslation[2] << ";" << std::endl;
   out << "   TGeoTranslation *" << GetPointerName() << " = new TGeoTranslation(\"" << GetName() << "\",dx,dy,dz);" << std::endl;
   TObject::SetBit(kGeoSavePrimitive);
}

/** \class TGeoRotation
\ingroup Geometry_classes
Class describing rotations. A rotation is a 3*3 array
Column vectors has to be orthogonal unit vectors.
*/

ClassImp(TGeoRotation);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TGeoRotation::TGeoRotation()
{
   for (Int_t i=0; i<9; i++) {
      if (i%4) fRotationMatrix[i] = 0;
      else fRotationMatrix[i] = 1.0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor.

TGeoRotation::TGeoRotation(const TGeoRotation &other)
             :TGeoMatrix(other)
{
   SetRotation(other);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor.

TGeoRotation::TGeoRotation(const TGeoMatrix &other)
             :TGeoMatrix(other)
{
   ResetBit(kGeoTranslation);
   ResetBit(kGeoScale);
   SetRotation(other);
}

////////////////////////////////////////////////////////////////////////////////
/// Named rotation constructor

TGeoRotation::TGeoRotation(const char *name)
             :TGeoMatrix(name)
{
   for (Int_t i=0; i<9; i++) {
      if (i%4) fRotationMatrix[i] = 0;
      else fRotationMatrix[i] = 1.0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Default rotation constructor with Euler angles. Phi is the rotation angle about
/// Z axis  and is done first, theta is the rotation about new X and is done
/// second, psi is the rotation angle about new Z and is done third. All angles are in
/// degrees.

TGeoRotation::TGeoRotation(const char *name, Double_t phi, Double_t theta, Double_t psi)
             :TGeoMatrix(name)
{
   SetAngles(phi, theta, psi);
}

////////////////////////////////////////////////////////////////////////////////
/// Rotation constructor a la GEANT3. Angles theta(i), phi(i) are the polar and azimuthal
/// angles of the (i) axis of the rotated system with respect to the initial non-rotated
/// system.
///   Example : the identity matrix (no rotation) is composed by
///      theta1=90, phi1=0, theta2=90, phi2=90, theta3=0, phi3=0
///   SetBit(kGeoRotation);

TGeoRotation::TGeoRotation(const char *name, Double_t theta1, Double_t phi1, Double_t theta2, Double_t phi2,
                           Double_t theta3, Double_t phi3)
             :TGeoMatrix(name)
{
   SetAngles(theta1, phi1, theta2, phi2, theta3, phi3);
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment from a general matrix

TGeoRotation& TGeoRotation::operator = (const TGeoMatrix &other)
{
   if (&other == this) return *this;
   Bool_t registered = TestBit(kGeoRegistered);
   TNamed::operator=(other);
   SetRotation(other);
   SetBit(kGeoRegistered,registered);
   ResetBit(kGeoTranslation);
   ResetBit(kGeoScale);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Composition

TGeoRotation &TGeoRotation::operator*=(const TGeoRotation &right)
{
   if (!right.IsRotation()) return *this;
   MultiplyBy(&right, true);
   return *this;
}

TGeoRotation TGeoRotation::operator*(const TGeoRotation &right) const
{
   TGeoRotation r = *this;
   r *= right;
   return r;
}

TGeoHMatrix TGeoRotation::operator*(const TGeoMatrix &right) const
{
   TGeoHMatrix t = *this;
   t *= right;
   return t;
}

////////////////////////////////////////////////////////////////////////////////
/// Is-equal operator

Bool_t TGeoRotation::operator ==(const TGeoRotation &other) const
{
   if (&other == this) return kTRUE;
   const Double_t *rot = GetRotationMatrix();
   const Double_t *orot = other.GetRotationMatrix();
   for (auto i=0; i<9; i++)
      if (TMath::Abs(rot[i]-orot[i])>1.E-10) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a temporary inverse of this.

TGeoHMatrix TGeoRotation::Inverse() const
{
   TGeoHMatrix h;
   h = *this;
   Double_t newrot[9];
   newrot[0] = fRotationMatrix[0];
   newrot[1] = fRotationMatrix[3];
   newrot[2] = fRotationMatrix[6];
   newrot[3] = fRotationMatrix[1];
   newrot[4] = fRotationMatrix[4];
   newrot[5] = fRotationMatrix[7];
   newrot[6] = fRotationMatrix[2];
   newrot[7] = fRotationMatrix[5];
   newrot[8] = fRotationMatrix[8];
   h.SetRotation(newrot);
   return h;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform orthogonality test for rotation.

Bool_t TGeoRotation::IsValid() const
{
   const Double_t *r = fRotationMatrix;
   Double_t cij;
   for (Int_t i=0; i<2; i++) {
      for (Int_t j=i+1; j<3; j++) {
         // check columns
         cij = TMath::Abs(r[i]*r[j]+r[i+3]*r[j+3]+r[i+6]*r[j+6]);
         if (cij>1E-4) return kFALSE;
         // check rows
         cij = TMath::Abs(r[3*i]*r[3*j]+r[3*i+1]*r[3*j+1]+r[3*i+2]*r[3*j+2]);
         if (cij>1E-4) return kFALSE;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// reset data members

void TGeoRotation::Clear(Option_t *)
{
   memcpy(fRotationMatrix,kIdentityMatrix,kN9);
   ResetBit(kGeoRotation);
}

////////////////////////////////////////////////////////////////////////////////
/// Perform a rotation about Z having the sine/cosine of the rotation angle.

void TGeoRotation::FastRotZ(const Double_t *sincos)
{
   fRotationMatrix[0] = sincos[1];
   fRotationMatrix[1] = -sincos[0];
   fRotationMatrix[3] = sincos[0];
   fRotationMatrix[4] = sincos[1];
   SetBit(kGeoRotation);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns rotation angle about Z axis in degrees. If the rotation is a pure
/// rotation about Z, fixX parameter does not matter, otherwise its meaning is:
///    - fixX = true  : result is the phi angle of the projection of the rotated X axis in the un-rotated XY
///    - fixX = false : result is the phi angle of the projection of the rotated Y axis - 90 degrees

Double_t TGeoRotation::GetPhiRotation(Bool_t fixX) const
{
   Double_t phi;
   if (fixX) phi = 180.*TMath::ATan2(-fRotationMatrix[1],fRotationMatrix[4])/TMath::Pi();
   else      phi = 180.*TMath::ATan2(fRotationMatrix[3], fRotationMatrix[0])/TMath::Pi();
   return phi;
}

////////////////////////////////////////////////////////////////////////////////
/// convert a point by multiplying its column vector (x, y, z, 1) to matrix inverse

void TGeoRotation::LocalToMaster(const Double_t *local, Double_t *master) const
{
   const Double_t *rot = GetRotationMatrix();
   for (Int_t i=0; i<3; i++) {
      master[i] = local[0]*rot[3*i]
                + local[1]*rot[3*i+1]
                + local[2]*rot[3*i+2];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// convert a point by multiplying its column vector (x, y, z, 1) to matrix

void TGeoRotation::MasterToLocal(const Double_t *master, Double_t *local) const
{
   const Double_t *rot = GetRotationMatrix();
   for (Int_t i=0; i<3; i++) {
      local[i] = master[0]*rot[i]
               + master[1]*rot[i+3]
               + master[2]*rot[i+6];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Make a clone of this matrix.

TGeoMatrix *TGeoRotation::MakeClone() const
{
   TGeoMatrix *matrix = new TGeoRotation(*this);
   return matrix;
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate about X axis of the master frame with angle expressed in degrees.

void TGeoRotation::RotateX(Double_t angle)
{
   SetBit(kGeoRotation);
   Double_t phi = angle*TMath::DegToRad();
   Double_t c = TMath::Cos(phi);
   Double_t s = TMath::Sin(phi);
   Double_t v[9];
   v[0] = fRotationMatrix[0];
   v[1] = fRotationMatrix[1];
   v[2] = fRotationMatrix[2];
   v[3] = c*fRotationMatrix[3]-s*fRotationMatrix[6];
   v[4] = c*fRotationMatrix[4]-s*fRotationMatrix[7];
   v[5] = c*fRotationMatrix[5]-s*fRotationMatrix[8];
   v[6] = s*fRotationMatrix[3]+c*fRotationMatrix[6];
   v[7] = s*fRotationMatrix[4]+c*fRotationMatrix[7];
   v[8] = s*fRotationMatrix[5]+c*fRotationMatrix[8];

   memcpy(fRotationMatrix, v, kN9);
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate about Y axis of the master frame with angle expressed in degrees.

void TGeoRotation::RotateY(Double_t angle)
{
   SetBit(kGeoRotation);
   Double_t phi = angle*TMath::DegToRad();
   Double_t c = TMath::Cos(phi);
   Double_t s = TMath::Sin(phi);
   Double_t v[9];
   v[0] = c*fRotationMatrix[0]+s*fRotationMatrix[6];
   v[1] = c*fRotationMatrix[1]+s*fRotationMatrix[7];
   v[2] = c*fRotationMatrix[2]+s*fRotationMatrix[8];
   v[3] = fRotationMatrix[3];
   v[4] = fRotationMatrix[4];
   v[5] = fRotationMatrix[5];
   v[6] = -s*fRotationMatrix[0]+c*fRotationMatrix[6];
   v[7] = -s*fRotationMatrix[1]+c*fRotationMatrix[7];
   v[8] = -s*fRotationMatrix[2]+c*fRotationMatrix[8];

   memcpy(fRotationMatrix, v, kN9);
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate about Z axis of the master frame with angle expressed in degrees.

void TGeoRotation::RotateZ(Double_t angle)
{
   SetBit(kGeoRotation);
   Double_t phi = angle*TMath::DegToRad();
   Double_t c = TMath::Cos(phi);
   Double_t s = TMath::Sin(phi);
   Double_t v[9];
   v[0] = c*fRotationMatrix[0]-s*fRotationMatrix[3];
   v[1] = c*fRotationMatrix[1]-s*fRotationMatrix[4];
   v[2] = c*fRotationMatrix[2]-s*fRotationMatrix[5];
   v[3] = s*fRotationMatrix[0]+c*fRotationMatrix[3];
   v[4] = s*fRotationMatrix[1]+c*fRotationMatrix[4];
   v[5] = s*fRotationMatrix[2]+c*fRotationMatrix[5];
   v[6] = fRotationMatrix[6];
   v[7] = fRotationMatrix[7];
   v[8] = fRotationMatrix[8];

   memcpy(&fRotationMatrix[0],v,kN9);
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply by a reflection respect to YZ.

void TGeoRotation::ReflectX(Bool_t leftside, Bool_t)
{
   if (leftside) {
      fRotationMatrix[0]=-fRotationMatrix[0];
      fRotationMatrix[1]=-fRotationMatrix[1];
      fRotationMatrix[2]=-fRotationMatrix[2];
   } else {
      fRotationMatrix[0]=-fRotationMatrix[0];
      fRotationMatrix[3]=-fRotationMatrix[3];
      fRotationMatrix[6]=-fRotationMatrix[6];
   }
   SetBit(kGeoRotation);
   SetBit(kGeoReflection, !IsReflection());
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply by a reflection respect to ZX.

void TGeoRotation::ReflectY(Bool_t leftside, Bool_t)
{
   if (leftside) {
      fRotationMatrix[3]=-fRotationMatrix[3];
      fRotationMatrix[4]=-fRotationMatrix[4];
      fRotationMatrix[5]=-fRotationMatrix[5];
   } else {
      fRotationMatrix[1]=-fRotationMatrix[1];
      fRotationMatrix[4]=-fRotationMatrix[4];
      fRotationMatrix[7]=-fRotationMatrix[7];
   }
   SetBit(kGeoRotation);
   SetBit(kGeoReflection, !IsReflection());
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply by a reflection respect to XY.

void TGeoRotation::ReflectZ(Bool_t leftside, Bool_t)
{
   if (leftside) {
      fRotationMatrix[6]=-fRotationMatrix[6];
      fRotationMatrix[7]=-fRotationMatrix[7];
      fRotationMatrix[8]=-fRotationMatrix[8];
   } else {
      fRotationMatrix[2]=-fRotationMatrix[2];
      fRotationMatrix[5]=-fRotationMatrix[5];
      fRotationMatrix[8]=-fRotationMatrix[8];
   }
   SetBit(kGeoRotation);
   SetBit(kGeoReflection, !IsReflection());
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoRotation::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TestBit(kGeoSavePrimitive)) return;
   out << "   // Rotation: " << GetName() << std::endl;
   Double_t th1,ph1,th2,ph2,th3,ph3;
   GetAngles(th1,ph1,th2,ph2,th3,ph3);
   out << "   thx = " << th1 << ";    phx = " << ph1 << ";" << std::endl;
   out << "   thy = " << th2 << ";    phy = " << ph2 << ";" << std::endl;
   out << "   thz = " << th3 << ";    phz = " << ph3 << ";" << std::endl;
   out << "   TGeoRotation *" << GetPointerName() << " = new TGeoRotation(\"" << GetName() << "\",thx,phx,thy,phy,thz,phz);" << std::endl;
   TObject::SetBit(kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy rotation elements from other rotation matrix.

void TGeoRotation::SetRotation(const TGeoMatrix &other)
{
   SetBit(kGeoRotation, other.IsRotation());
   SetMatrix(other.GetRotationMatrix());
}

////////////////////////////////////////////////////////////////////////////////
/// Set matrix elements according to Euler angles. Phi is the rotation angle about
/// Z axis  and is done first, theta is the rotation about new X and is done
/// second, psi is the rotation angle about new Z and is done third. All angles are in
/// degrees.

void TGeoRotation::SetAngles(Double_t phi, Double_t theta, Double_t psi)
{
   Double_t degrad = TMath::Pi()/180.;
   Double_t sinphi = TMath::Sin(degrad*phi);
   Double_t cosphi = TMath::Cos(degrad*phi);
   Double_t sinthe = TMath::Sin(degrad*theta);
   Double_t costhe = TMath::Cos(degrad*theta);
   Double_t sinpsi = TMath::Sin(degrad*psi);
   Double_t cospsi = TMath::Cos(degrad*psi);

   fRotationMatrix[0] =  cospsi*cosphi - costhe*sinphi*sinpsi;
   fRotationMatrix[1] = -sinpsi*cosphi - costhe*sinphi*cospsi;
   fRotationMatrix[2] =  sinthe*sinphi;
   fRotationMatrix[3] =  cospsi*sinphi + costhe*cosphi*sinpsi;
   fRotationMatrix[4] = -sinpsi*sinphi + costhe*cosphi*cospsi;
   fRotationMatrix[5] = -sinthe*cosphi;
   fRotationMatrix[6] =  sinpsi*sinthe;
   fRotationMatrix[7] =  cospsi*sinthe;
   fRotationMatrix[8] =  costhe;

   if (!IsValid()) Error("SetAngles", "invalid rotation (Euler angles : phi=%f theta=%f psi=%f)",phi,theta,psi);
   CheckMatrix();
}

////////////////////////////////////////////////////////////////////////////////
/// Set matrix elements in the GEANT3 way

void TGeoRotation::SetAngles(Double_t theta1, Double_t phi1, Double_t theta2, Double_t phi2,
                             Double_t theta3, Double_t phi3)
{
   Double_t degrad = TMath::Pi()/180.;
   fRotationMatrix[0] = TMath::Cos(degrad*phi1)*TMath::Sin(degrad*theta1);
   fRotationMatrix[3] = TMath::Sin(degrad*phi1)*TMath::Sin(degrad*theta1);
   fRotationMatrix[6] = TMath::Cos(degrad*theta1);
   fRotationMatrix[1] = TMath::Cos(degrad*phi2)*TMath::Sin(degrad*theta2);
   fRotationMatrix[4] = TMath::Sin(degrad*phi2)*TMath::Sin(degrad*theta2);
   fRotationMatrix[7] = TMath::Cos(degrad*theta2);
   fRotationMatrix[2] = TMath::Cos(degrad*phi3)*TMath::Sin(degrad*theta3);
   fRotationMatrix[5] = TMath::Sin(degrad*phi3)*TMath::Sin(degrad*theta3);
   fRotationMatrix[8] = TMath::Cos(degrad*theta3);
   // do the trick to eliminate most of the floating point errors
   for (Int_t i=0; i<9; i++) {
      if (TMath::Abs(fRotationMatrix[i])<1E-15) fRotationMatrix[i] = 0;
      if (TMath::Abs(fRotationMatrix[i]-1)<1E-15) fRotationMatrix[i] = 1;
      if (TMath::Abs(fRotationMatrix[i]+1)<1E-15) fRotationMatrix[i] = -1;
   }
   if (!IsValid()) Error("SetAngles", "invalid rotation (G3 angles, th1=%f phi1=%f, th2=%f ph2=%f, th3=%f phi3=%f)",
                          theta1,phi1,theta2,phi2,theta3,phi3);
   CheckMatrix();
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve rotation angles

void TGeoRotation::GetAngles(Double_t &theta1, Double_t &phi1, Double_t &theta2, Double_t &phi2,
                             Double_t &theta3, Double_t &phi3) const
{
   Double_t raddeg = 180./TMath::Pi();
   theta1 = raddeg*TMath::ACos(fRotationMatrix[6]);
   theta2 = raddeg*TMath::ACos(fRotationMatrix[7]);
   theta3 = raddeg*TMath::ACos(fRotationMatrix[8]);
   if (TMath::Abs(fRotationMatrix[0])<1E-6 && TMath::Abs(fRotationMatrix[3])<1E-6) phi1=0.;
   else phi1 = raddeg*TMath::ATan2(fRotationMatrix[3],fRotationMatrix[0]);
   if (phi1<0) phi1+=360.;
   if (TMath::Abs(fRotationMatrix[1])<1E-6 && TMath::Abs(fRotationMatrix[4])<1E-6) phi2=0.;
   else phi2 = raddeg*TMath::ATan2(fRotationMatrix[4],fRotationMatrix[1]);
   if (phi2<0) phi2+=360.;
   if (TMath::Abs(fRotationMatrix[2])<1E-6 && TMath::Abs(fRotationMatrix[5])<1E-6) phi3=0.;
   else phi3 = raddeg*TMath::ATan2(fRotationMatrix[5],fRotationMatrix[2]);
   if (phi3<0) phi3+=360.;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve Euler angles.

void TGeoRotation::GetAngles(Double_t &phi, Double_t &theta, Double_t &psi) const
{
   const Double_t *m = fRotationMatrix;
   // Check if theta is 0 or 180.
   if (TMath::Abs(1.-TMath::Abs(m[8]))<1.e-9) {
      theta = TMath::ACos(m[8])*TMath::RadToDeg();
      phi = TMath::ATan2(-m[8]*m[1],m[0])*TMath::RadToDeg();
      psi = 0.; // convention, phi+psi matters
      return;
   }
   // sin(theta) != 0
   phi = TMath::ATan2(m[2],-m[5]);
   Double_t sphi = TMath::Sin(phi);
   if (TMath::Abs(sphi)<1.e-9) theta = -TMath::ASin(m[5]/TMath::Cos(phi))*TMath::RadToDeg();
   else theta = TMath::ASin(m[2]/sphi)*TMath::RadToDeg();
   phi *= TMath::RadToDeg();
   psi = TMath::ATan2(m[6],m[7])*TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// computes determinant of the rotation matrix

Double_t TGeoRotation::Determinant() const
{
   Double_t
   det = fRotationMatrix[0]*fRotationMatrix[4]*fRotationMatrix[8] +
         fRotationMatrix[3]*fRotationMatrix[7]*fRotationMatrix[2] +
         fRotationMatrix[6]*fRotationMatrix[1]*fRotationMatrix[5] -
         fRotationMatrix[2]*fRotationMatrix[4]*fRotationMatrix[6] -
         fRotationMatrix[5]*fRotationMatrix[7]*fRotationMatrix[0] -
         fRotationMatrix[8]*fRotationMatrix[1]*fRotationMatrix[3];
   return det;
}

////////////////////////////////////////////////////////////////////////////////
/// performes an orthogonality check and finds if the matrix is a reflection
///   Warning("CheckMatrix", "orthogonality check not performed yet");

void TGeoRotation::CheckMatrix()
{
   if (Determinant() < 0) SetBit(kGeoReflection);
   Double_t dd = fRotationMatrix[0] + fRotationMatrix[4] + fRotationMatrix[8] - 3.;
   if (TMath::Abs(dd) < 1.E-12) ResetBit(kGeoRotation);
   else                         SetBit(kGeoRotation);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the inverse rotation matrix (which is simply the transpose)

void TGeoRotation::GetInverse(Double_t *invmat) const
{
   if (!invmat) {
      Error("GetInverse", "no place to store the inverse matrix");
      return;
   }
   for (Int_t i=0; i<3; i++) {
      for (Int_t j=0; j<3; j++) {
         invmat[3*i+j] = fRotationMatrix[3*j+i];
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply this rotation with the one specified by ROT.
/// -   after=TRUE (default): THIS*ROT
/// -   after=FALSE         : ROT*THIS

void TGeoRotation::MultiplyBy(const TGeoRotation *rot, Bool_t after)
{
   const Double_t *matleft, *matright;
   SetBit(kGeoRotation);
   Double_t  newmat[9] = {0};
   if (after) {
      matleft  = &fRotationMatrix[0];
      matright = rot->GetRotationMatrix();
   } else {
      matleft  = rot->GetRotationMatrix();
      matright = &fRotationMatrix[0];
   }
   for (Int_t i=0; i<3; i++) {
      for (Int_t j=0; j<3; j++) {
         for (Int_t k=0; k<3; k++) {
            newmat[3*i+j] += matleft[3*i+k] * matright[3*k+j];
         }
      }
   }
   memcpy(&fRotationMatrix[0], &newmat[0], kN9);
}

/** \class TGeoScale
\ingroup Geometry_classes
Class describing scale transformations. A scale is an
array of 3 doubles (sx, sy, sz) multiplying elements 0, 5 and 10
of the homogenous matrix. A scale is normalized : sx*sy*sz = 1
*/

ClassImp(TGeoScale);

////////////////////////////////////////////////////////////////////////////////
/// default constructor

TGeoScale::TGeoScale()
{
   SetBit(kGeoScale);
   for (Int_t i=0; i<3; i++) fScale[i] = 1.;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

TGeoScale::TGeoScale(const TGeoScale &other)
          :TGeoMatrix(other)
{
   SetScale(other);
}

////////////////////////////////////////////////////////////////////////////////
/// Ctor. based on a general matrix

TGeoScale::TGeoScale(const TGeoMatrix &other)
                :TGeoMatrix(other)
{
   ResetBit(kGeoTranslation);
   ResetBit(kGeoRotation);
   SetScale(other);
}

////////////////////////////////////////////////////////////////////////////////
/// default constructor

TGeoScale::TGeoScale(Double_t sx, Double_t sy, Double_t sz)
          :TGeoMatrix("")
{
   SetBit(kGeoScale);
   SetScale(sx, sy, sz);
}

////////////////////////////////////////////////////////////////////////////////
/// default constructor

TGeoScale::TGeoScale(const char *name, Double_t sx, Double_t sy, Double_t sz)
          :TGeoMatrix(name)
{
   SetBit(kGeoScale);
   SetScale(sx, sy, sz);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoScale::~TGeoScale()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment from a general matrix

TGeoScale& TGeoScale::operator = (const TGeoMatrix &matrix)
{
   if (&matrix == this) return *this;
   Bool_t registered = TestBit(kGeoRegistered);
   TNamed::operator=(matrix);
   SetScale(matrix);
   SetBit(kGeoRegistered,registered);
   ResetBit(kGeoTranslation);
   ResetBit(kGeoRotation);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Scale composition

TGeoScale &TGeoScale::operator*=(const TGeoScale &right)
{
   const Double_t *scl = right.GetScale();
   fScale[0] *= scl[0];
   fScale[1] *= scl[1];
   fScale[2] *= scl[2];
   SetBit(kGeoReflection, fScale[0] * fScale[1] * fScale[2] < 0);
   if (!IsScale()) SetBit(kGeoScale, right.IsScale());
   return *this;
}

TGeoScale TGeoScale::operator*(const TGeoScale &right) const
{
   TGeoScale s = *this;
   s *= right;
   return s;
}

TGeoHMatrix TGeoScale::operator*(const TGeoMatrix &right) const
{
   TGeoHMatrix t = *this;
   t *= right;
   return t;
}

////////////////////////////////////////////////////////////////////////////////
/// Is-equal operator

Bool_t TGeoScale::operator ==(const TGeoScale &other) const
{
   if (&other == this) return kTRUE;
   const Double_t *scl = GetScale();
   const Double_t *oscl = other.GetScale();
   for (auto i=0; i<3; i++)
      if (TMath::Abs(scl[i]-oscl[i])>1.E-10) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a temporary inverse of this.

TGeoHMatrix TGeoScale::Inverse() const
{
   TGeoHMatrix h;
   h = *this;
   Double_t scale[3];
   scale[0] = 1./fScale[0];
   scale[1] = 1./fScale[1];
   scale[2] = 1./fScale[2];
   h.SetScale(scale);
   return h;
}

////////////////////////////////////////////////////////////////////////////////
/// scale setter

void TGeoScale::SetScale(Double_t sx, Double_t sy, Double_t sz)
{
   if (TMath::Abs(sx*sy*sz) < 1.E-10) {
      Error("SetScale", "Invalid scale %f, %f, %f for transformation %s",sx,sy,sx,GetName());
      return;
   }
   fScale[0] = sx;
   fScale[1] = sy;
   fScale[2] = sz;
   if (sx*sy*sz<0) SetBit(kGeoReflection);
   else            SetBit(kGeoReflection, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set scale from other transformation

void TGeoScale::SetScale(const TGeoMatrix &other)
{
   SetBit(kGeoScale, other.IsScale());
   SetBit(kGeoReflection, other.IsReflection());
   memcpy(fScale, other.GetScale(), kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert a local point to the master frame.

void TGeoScale::LocalToMaster(const Double_t *local, Double_t *master) const
{
   master[0] = local[0]*fScale[0];
   master[1] = local[1]*fScale[1];
   master[2] = local[2]*fScale[2];
}

////////////////////////////////////////////////////////////////////////////////
/// Convert the local distance along unit vector DIR to master frame. If DIR
/// is not specified perform a conversion such as the returned distance is the
/// the minimum for all possible directions.

Double_t TGeoScale::LocalToMaster(Double_t dist, const Double_t *dir) const
{
   Double_t scale;
   if (!dir) {
      scale = TMath::Abs(fScale[0]);
      if (TMath::Abs(fScale[1])<scale) scale = TMath::Abs(fScale[1]);
      if (TMath::Abs(fScale[2])<scale) scale = TMath::Abs(fScale[2]);
   } else {
      scale = fScale[0]*fScale[0]*dir[0]*dir[0] +
              fScale[1]*fScale[1]*dir[1]*dir[1] +
              fScale[2]*fScale[2]*dir[2]*dir[2];
      scale = TMath::Sqrt(scale);
   }
   return scale*dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a clone of this matrix.

TGeoMatrix *TGeoScale::MakeClone() const
{
   TGeoMatrix *matrix = new TGeoScale(*this);
   return matrix;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert a global point to local frame.

void TGeoScale::MasterToLocal(const Double_t *master, Double_t *local) const
{
   local[0] = master[0]/fScale[0];
   local[1] = master[1]/fScale[1];
   local[2] = master[2]/fScale[2];
}

////////////////////////////////////////////////////////////////////////////////
/// Convert the distance along unit vector DIR to local frame. If DIR
/// is not specified perform a conversion such as the returned distance is the
/// the minimum for all possible directions.

Double_t TGeoScale::MasterToLocal(Double_t dist, const Double_t *dir) const
{
   Double_t scale;
   if (!dir) {
      scale = TMath::Abs(fScale[0]);
      if (TMath::Abs(fScale[1])>scale) scale = TMath::Abs(fScale[1]);
      if (TMath::Abs(fScale[2])>scale) scale = TMath::Abs(fScale[2]);
      scale = 1./scale;
   } else {
      scale = (dir[0]*dir[0])/(fScale[0]*fScale[0]) +
              (dir[1]*dir[1])/(fScale[1]*fScale[1]) +
              (dir[2]*dir[2])/(fScale[2]*fScale[2]);
      scale = TMath::Sqrt(scale);
   }
   return scale*dist;
}

/** \class TGeoCombiTrans
\ingroup Geometry_classes
Class describing rotation + translation. Most frequently used in the description
of TGeoNode 's
*/

ClassImp(TGeoCombiTrans);

////////////////////////////////////////////////////////////////////////////////
/// dummy ctor

TGeoCombiTrans::TGeoCombiTrans()
{
   for (Int_t i=0; i<3; i++) fTranslation[i] = 0.0;
   fRotation = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor from generic matrix.

TGeoCombiTrans::TGeoCombiTrans(const TGeoMatrix &other)
               :TGeoMatrix(other)
{
   ResetBit(kGeoScale);
   if (other.IsTranslation()) {
      SetBit(kGeoTranslation);
      memcpy(fTranslation,other.GetTranslation(),kN3);
   } else {
      for (Int_t i=0; i<3; i++) fTranslation[i] = 0.0;
   }
   if (other.IsRotation()) {
      SetBit(kGeoRotation);
      SetBit(kGeoMatrixOwned);
      fRotation = new TGeoRotation(other);
   }
   else fRotation = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a translation and a rotation.

TGeoCombiTrans::TGeoCombiTrans(const TGeoTranslation &tr, const TGeoRotation &rot)
{
   if (tr.IsTranslation()) {
      SetBit(kGeoTranslation);
      const Double_t *trans = tr.GetTranslation();
      memcpy(fTranslation, trans, kN3);
   } else {
      for (Int_t i=0; i<3; i++) fTranslation[i] = 0.0;
   }
   if (rot.IsRotation()) {
      SetBit(kGeoRotation);
      SetBit(kGeoMatrixOwned);
      fRotation = new TGeoRotation(rot);
      SetBit(kGeoReflection, rot.TestBit(kGeoReflection));
   }
   else fRotation = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Named ctor.

TGeoCombiTrans::TGeoCombiTrans(const char *name)
               :TGeoMatrix(name)
{
   for (Int_t i=0; i<3; i++) fTranslation[i] = 0.0;
   fRotation = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a translation specified by X,Y,Z and a pointer to a rotation. The rotation will not be owned by this.

TGeoCombiTrans::TGeoCombiTrans(Double_t dx, Double_t dy, Double_t dz, TGeoRotation *rot)
               :TGeoMatrix("")
{
   SetTranslation(dx, dy, dz);
   fRotation = 0;
   SetRotation(rot);
}

////////////////////////////////////////////////////////////////////////////////
/// Named ctor

TGeoCombiTrans::TGeoCombiTrans(const char *name, Double_t dx, Double_t dy, Double_t dz, TGeoRotation *rot)
               :TGeoMatrix(name)
{
   SetTranslation(dx, dy, dz);
   fRotation = 0;
   SetRotation(rot);
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator with generic matrix.

TGeoCombiTrans &TGeoCombiTrans::operator=(const TGeoMatrix &matrix)
{
   if (&matrix == this) return *this;
   Bool_t registered = TestBit(kGeoRegistered);
   Clear();
   TNamed::operator=(matrix);

   if (matrix.IsTranslation()) {
      memcpy(fTranslation,matrix.GetTranslation(),kN3);
   }
   if (matrix.IsRotation()) {
      if (!fRotation) {
         fRotation = new TGeoRotation();
         SetBit(kGeoMatrixOwned);
      } else {
         if (!TestBit(kGeoMatrixOwned)) {
            fRotation = new TGeoRotation();
            SetBit(kGeoMatrixOwned);
         }
      }
      fRotation->SetMatrix(matrix.GetRotationMatrix());
      fRotation->SetBit(kGeoReflection, matrix.TestBit(kGeoReflection));
      fRotation->SetBit(kGeoRotation);
   } else {
      if (fRotation && TestBit(kGeoMatrixOwned)) delete fRotation;
      ResetBit(kGeoMatrixOwned);
      fRotation = 0;
   }
   SetBit(kGeoRegistered,registered);
   ResetBit(kGeoScale);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Is-equal operator

Bool_t TGeoCombiTrans::operator==(const TGeoMatrix &other) const
{
   if (&other == this) return kTRUE;
   const Double_t *tr = GetTranslation();
   const Double_t *otr = other.GetTranslation();
   for (auto i=0; i<3; i++) if (TMath::Abs(tr[i]-otr[i])>1.E-10) return kFALSE;
   const Double_t *rot = GetRotationMatrix();
   const Double_t *orot = other.GetRotationMatrix();
   for (auto i=0; i<9; i++) if (TMath::Abs(rot[i]-orot[i])>1.E-10) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Composition

TGeoCombiTrans &TGeoCombiTrans::operator*=(const TGeoMatrix &right)
{
   Multiply(&right);
   return *this;
}

TGeoCombiTrans TGeoCombiTrans::operator*(const TGeoMatrix &right) const
{
   TGeoHMatrix h = *this;
   h *= right;
   return h;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoCombiTrans::~TGeoCombiTrans()
{
   if (fRotation) {
      if(TestBit(TGeoMatrix::kGeoMatrixOwned) && !fRotation->IsRegistered()) delete fRotation;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset translation/rotation to identity

void TGeoCombiTrans::Clear(Option_t *)
{
   if (IsTranslation()) {
      ResetBit(kGeoTranslation);
      memset(fTranslation, 0, kN3);
   }
   if (fRotation) {
      if (TestBit(kGeoMatrixOwned)) delete fRotation;
      fRotation = 0;
   }
   ResetBit(kGeoRotation);
   ResetBit(kGeoTranslation);
   ResetBit(kGeoMatrixOwned);
}

////////////////////////////////////////////////////////////////////////////////
/// Return a temporary inverse of this.

TGeoHMatrix TGeoCombiTrans::Inverse() const
{
   TGeoHMatrix h;
   h = *this;
   Bool_t is_tr = IsTranslation();
   Bool_t is_rot = IsRotation();
   Double_t tr[3];
   Double_t newrot[9];
   const Double_t *rot = GetRotationMatrix();
   tr[0] = -fTranslation[0]*rot[0] - fTranslation[1]*rot[3] - fTranslation[2]*rot[6];
   tr[1] = -fTranslation[0]*rot[1] - fTranslation[1]*rot[4] - fTranslation[2]*rot[7];
   tr[2] = -fTranslation[0]*rot[2] - fTranslation[1]*rot[5] - fTranslation[2]*rot[8];
   h.SetTranslation(tr);
   newrot[0] = rot[0];
   newrot[1] = rot[3];
   newrot[2] = rot[6];
   newrot[3] = rot[1];
   newrot[4] = rot[4];
   newrot[5] = rot[7];
   newrot[6] = rot[2];
   newrot[7] = rot[5];
   newrot[8] = rot[8];
   h.SetRotation(newrot);
   h.SetBit(kGeoTranslation,is_tr);
   h.SetBit(kGeoRotation,is_rot);
   return h;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a clone of this matrix.

TGeoMatrix *TGeoCombiTrans::MakeClone() const
{
   TGeoMatrix *matrix = new TGeoCombiTrans(*this);
   return matrix;
}

////////////////////////////////////////////////////////////////////////////////
/// multiply to the right with an other transformation
/// if right is identity matrix, just return

void TGeoCombiTrans::Multiply(const TGeoMatrix *right)
{
   if (right->IsIdentity()) return;
   TGeoHMatrix h = *this;
   h.Multiply(right);
   operator=(h);
}

////////////////////////////////////////////////////////////////////////////////
/// Register the matrix in the current manager, which will become the owner.

void TGeoCombiTrans::RegisterYourself()
{
   TGeoMatrix::RegisterYourself();
   if (fRotation && fRotation->IsRotation()) fRotation->RegisterYourself();
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate about X axis with angle expressed in degrees.

void TGeoCombiTrans::RotateX(Double_t angle)
{
   if (!fRotation || !TestBit(kGeoMatrixOwned)) {
      if (fRotation) fRotation = new TGeoRotation(*fRotation);
      else fRotation = new TGeoRotation();
      SetBit(kGeoMatrixOwned);
   }
   SetBit(kGeoRotation);
   const Double_t *rot = fRotation->GetRotationMatrix();
   Double_t phi = angle*TMath::DegToRad();
   Double_t c = TMath::Cos(phi);
   Double_t s = TMath::Sin(phi);
   Double_t v[9];
   v[0] = rot[0];
   v[1] = rot[1];
   v[2] = rot[2];
   v[3] = c*rot[3]-s*rot[6];
   v[4] = c*rot[4]-s*rot[7];
   v[5] = c*rot[5]-s*rot[8];
   v[6] = s*rot[3]+c*rot[6];
   v[7] = s*rot[4]+c*rot[7];
   v[8] = s*rot[5]+c*rot[8];
   fRotation->SetMatrix(v);
   fRotation->SetBit(kGeoRotation);
   if (!IsTranslation()) return;
   v[0] = fTranslation[0];
   v[1] = c*fTranslation[1]-s*fTranslation[2];
   v[2] = s*fTranslation[1]+c*fTranslation[2];
   memcpy(fTranslation,v,kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate about Y axis with angle expressed in degrees.

void TGeoCombiTrans::RotateY(Double_t angle)
{
   if (!fRotation || !TestBit(kGeoMatrixOwned)) {
      if (fRotation) fRotation = new TGeoRotation(*fRotation);
      else fRotation = new TGeoRotation();
      SetBit(kGeoMatrixOwned);
   }
   SetBit(kGeoRotation);
   const Double_t *rot = fRotation->GetRotationMatrix();
   Double_t phi = angle*TMath::DegToRad();
   Double_t c = TMath::Cos(phi);
   Double_t s = TMath::Sin(phi);
   Double_t v[9];
   v[0] = c*rot[0]+s*rot[6];
   v[1] = c*rot[1]+s*rot[7];
   v[2] = c*rot[2]+s*rot[8];
   v[3] = rot[3];
   v[4] = rot[4];
   v[5] = rot[5];
   v[6] = -s*rot[0]+c*rot[6];
   v[7] = -s*rot[1]+c*rot[7];
   v[8] = -s*rot[2]+c*rot[8];
   fRotation->SetMatrix(v);
   fRotation->SetBit(kGeoRotation);
   if (!IsTranslation()) return;
   v[0] = c*fTranslation[0]+s*fTranslation[2];
   v[1] = fTranslation[1];
   v[2] = -s*fTranslation[0]+c*fTranslation[2];
   memcpy(fTranslation,v,kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate about Z axis with angle expressed in degrees.

void TGeoCombiTrans::RotateZ(Double_t angle)
{
   if (!fRotation || !TestBit(kGeoMatrixOwned)) {
      if (fRotation) fRotation = new TGeoRotation(*fRotation);
      else fRotation = new TGeoRotation();
      SetBit(kGeoMatrixOwned);
   }
   SetBit(kGeoRotation);
   const Double_t *rot = fRotation->GetRotationMatrix();
   Double_t phi = angle*TMath::DegToRad();
   Double_t c = TMath::Cos(phi);
   Double_t s = TMath::Sin(phi);
   Double_t v[9];
   v[0] = c*rot[0]-s*rot[3];
   v[1] = c*rot[1]-s*rot[4];
   v[2] = c*rot[2]-s*rot[5];
   v[3] = s*rot[0]+c*rot[3];
   v[4] = s*rot[1]+c*rot[4];
   v[5] = s*rot[2]+c*rot[5];
   v[6] = rot[6];
   v[7] = rot[7];
   v[8] = rot[8];
   fRotation->SetMatrix(v);
   fRotation->SetBit(kGeoRotation);
   if (!IsTranslation()) return;
   v[0] = c*fTranslation[0]-s*fTranslation[1];
   v[1] = s*fTranslation[0]+c*fTranslation[1];
   v[2] = fTranslation[2];
   memcpy(fTranslation,v,kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply by a reflection respect to YZ.

void TGeoCombiTrans::ReflectX(Bool_t leftside, Bool_t rotonly)
{
   if (leftside && !rotonly) fTranslation[0] = -fTranslation[0];
   if (!fRotation || !TestBit(kGeoMatrixOwned)) {
      if (fRotation) fRotation = new TGeoRotation(*fRotation);
      else fRotation = new TGeoRotation();
      SetBit(kGeoMatrixOwned);
   }
   SetBit(kGeoRotation);
   fRotation->ReflectX(leftside);
   SetBit(kGeoReflection, !IsReflection());
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply by a reflection respect to ZX.

void TGeoCombiTrans::ReflectY(Bool_t leftside, Bool_t rotonly)
{
   if (leftside && !rotonly) fTranslation[1] = -fTranslation[1];
   if (!fRotation || !TestBit(kGeoMatrixOwned)) {
      if (fRotation) fRotation = new TGeoRotation(*fRotation);
      else fRotation = new TGeoRotation();
      SetBit(kGeoMatrixOwned);
   }
   SetBit(kGeoRotation);
   fRotation->ReflectY(leftside);
   SetBit(kGeoReflection, !IsReflection());
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply by a reflection respect to XY.

void TGeoCombiTrans::ReflectZ(Bool_t leftside, Bool_t rotonly)
{
   if (leftside && !rotonly) fTranslation[2] = -fTranslation[2];
   if (!fRotation || !TestBit(kGeoMatrixOwned)) {
      if (fRotation) fRotation = new TGeoRotation(*fRotation);
      else fRotation = new TGeoRotation();
      SetBit(kGeoMatrixOwned);
   }
   SetBit(kGeoRotation);
   fRotation->ReflectZ(leftside);
   SetBit(kGeoReflection, !IsReflection());
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoCombiTrans::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (TestBit(kGeoSavePrimitive)) return;
   out << "   // Combi transformation: " << GetName() << std::endl;
   out << "   dx = " << fTranslation[0] << ";" << std::endl;
   out << "   dy = " << fTranslation[1] << ";" << std::endl;
   out << "   dz = " << fTranslation[2] << ";" << std::endl;
   if (fRotation && fRotation->IsRotation()) {
      fRotation->SavePrimitive(out,option);
      out << "   " << GetPointerName() << " = new TGeoCombiTrans(\"" << GetName() << "\", dx,dy,dz,";
      out << fRotation->GetPointerName() << ");" << std::endl;
   } else {
      out << "   " << GetPointerName() << " = new TGeoCombiTrans(\"" << GetName() << "\");" << std::endl;
      out << "   " << GetPointerName() << "->SetTranslation(dx,dy,dz);" << std::endl;
   }
   TObject::SetBit(kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Assign a foreign rotation to the combi. The rotation is NOT owned by this.

void TGeoCombiTrans::SetRotation(const TGeoRotation *rot)
{
   if (fRotation && TestBit(kGeoMatrixOwned)) delete fRotation;
   fRotation = 0;
   ResetBit(TGeoMatrix::kGeoMatrixOwned);
   ResetBit(kGeoRotation);
   ResetBit(kGeoReflection);
   if (!rot) return;
   if (!rot->IsRotation()) return;

   SetBit(kGeoRotation);
   SetBit(kGeoReflection, rot->TestBit(kGeoReflection));
   TGeoRotation *rr = (TGeoRotation*)rot;
   fRotation = rr;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the rotation from another one.

void TGeoCombiTrans::SetRotation(const TGeoRotation &rot)
{
   if (fRotation && TestBit(kGeoMatrixOwned)) delete fRotation;
   fRotation = 0;
   if (!rot.IsRotation()) {
      ResetBit(kGeoRotation);
      ResetBit(kGeoReflection);
      ResetBit(TGeoMatrix::kGeoMatrixOwned);
      return;
   }

   SetBit(kGeoRotation);
   SetBit(kGeoReflection, rot.TestBit(kGeoReflection));
   fRotation = new TGeoRotation(rot);
   SetBit(kGeoMatrixOwned);
}

////////////////////////////////////////////////////////////////////////////////
/// copy the translation component

void TGeoCombiTrans::SetTranslation(const TGeoTranslation &tr)
{
   if (tr.IsTranslation()) {
      SetBit(kGeoTranslation);
      const Double_t *trans = tr.GetTranslation();
      memcpy(fTranslation, trans, kN3);
   } else {
      if (!IsTranslation()) return;
      memset(fTranslation, 0, kN3);
      ResetBit(kGeoTranslation);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// set the translation component

void TGeoCombiTrans::SetTranslation(Double_t dx, Double_t dy, Double_t dz)
{
   fTranslation[0] = dx;
   fTranslation[1] = dy;
   fTranslation[2] = dz;
   if (fTranslation[0] || fTranslation[1] || fTranslation[2]) SetBit(kGeoTranslation);
   else ResetBit(kGeoTranslation);
}

////////////////////////////////////////////////////////////////////////////////
/// set the translation component

void TGeoCombiTrans::SetTranslation(Double_t *vect)
{
   fTranslation[0] = vect[0];
   fTranslation[1] = vect[1];
   fTranslation[2] = vect[2];
   if (fTranslation[0] || fTranslation[1] || fTranslation[2]) SetBit(kGeoTranslation);
   else ResetBit(kGeoTranslation);
}

////////////////////////////////////////////////////////////////////////////////
/// get the rotation array

const Double_t *TGeoCombiTrans::GetRotationMatrix() const
{
   if (!fRotation) return kIdentityMatrix;
   return fRotation->GetRotationMatrix();
}

/** \class TGeoGenTrans
\ingroup Geometry_classes
Most general transformation, holding a translation, a rotation and a scale
*/

ClassImp(TGeoGenTrans);

////////////////////////////////////////////////////////////////////////////////
/// dummy ctor

TGeoGenTrans::TGeoGenTrans()
{
   SetBit(kGeoGenTrans);
   for (Int_t i=0; i<3; i++) fTranslation[i] = 0.0;
   for (Int_t j=0; j<3; j++) fScale[j] = 1.0;
   fRotation = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoGenTrans::TGeoGenTrans(const char *name)
             :TGeoCombiTrans(name)
{
   SetBit(kGeoGenTrans);
   for (Int_t i=0; i<3; i++) fTranslation[i] = 0.0;
   for (Int_t j=0; j<3; j++) fScale[j] = 1.0;
   fRotation = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoGenTrans::TGeoGenTrans(Double_t dx, Double_t dy, Double_t dz,
                           Double_t sx, Double_t sy, Double_t sz, TGeoRotation *rot)
             :TGeoCombiTrans("")
{
   SetBit(kGeoGenTrans);
   SetTranslation(dx, dy, dz);
   SetScale(sx, sy, sz);
   SetRotation(rot);
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoGenTrans::TGeoGenTrans(const char *name, Double_t dx, Double_t dy, Double_t dz,
                           Double_t sx, Double_t sy, Double_t sz, TGeoRotation *rot)
             :TGeoCombiTrans(name)
{
   SetBit(kGeoGenTrans);
   SetTranslation(dx, dy, dz);
   SetScale(sx, sy, sz);
   SetRotation(rot);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoGenTrans::~TGeoGenTrans()
{
}

////////////////////////////////////////////////////////////////////////////////
/// clear the fields of this transformation

void TGeoGenTrans::Clear(Option_t *)
{
   memset(&fTranslation[0], 0, kN3);
   memset(&fScale[0], 0, kN3);
   if (fRotation) fRotation->Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// set the scale

void TGeoGenTrans::SetScale(Double_t sx, Double_t sy, Double_t sz)
{
   if (sx<1.E-5 || sy<1.E-5 || sz<1.E-5) {
      Error("ctor", "Invalid scale");
      return;
   }
   fScale[0] = sx;
   fScale[1] = sy;
   fScale[2] = sz;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a temporary inverse of this.

TGeoHMatrix TGeoGenTrans::Inverse() const
{
   TGeoHMatrix h = *this;
   return h;
}

////////////////////////////////////////////////////////////////////////////////
/// A scale transformation should be normalized by sx*sy*sz factor

Bool_t TGeoGenTrans::Normalize()
{
   Double_t normfactor = fScale[0]*fScale[1]*fScale[2];
   if (normfactor <= 1E-5) return kFALSE;
   for (Int_t i=0; i<3; i++)
      fScale[i] /= normfactor;
   return kTRUE;
}

/** \class TGeoIdentity
\ingroup Geometry_classes
An identity transformation. It holds no data member
and returns pointers to static null translation and identity
transformations for rotation and scale
*/

ClassImp(TGeoIdentity);

////////////////////////////////////////////////////////////////////////////////
/// dummy ctor

TGeoIdentity::TGeoIdentity()
{
   if (!gGeoIdentity) gGeoIdentity = this;
   RegisterYourself();
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoIdentity::TGeoIdentity(const char *name)
             :TGeoMatrix(name)
{
   if (!gGeoIdentity) gGeoIdentity = this;
   RegisterYourself();
}

////////////////////////////////////////////////////////////////////////////////
/// Return a temporary inverse of this.

TGeoHMatrix TGeoIdentity::Inverse() const
{
   TGeoHMatrix h = *gGeoIdentity;
   return h;
}

/** \class TGeoHMatrix
\ingroup Geometry_classes

Matrix class used for computing global transformations
Should NOT be used for node definition. An instance of this class
is generally used to pile-up local transformations starting from
the top level physical node, down to the current node.
*/

ClassImp(TGeoHMatrix);

////////////////////////////////////////////////////////////////////////////////
/// dummy ctor

TGeoHMatrix::TGeoHMatrix()
{
   memset(&fTranslation[0], 0, kN3);
   memcpy(fRotationMatrix,kIdentityMatrix,kN9);
   memcpy(fScale,kUnitScale,kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoHMatrix::TGeoHMatrix(const char* name)
            :TGeoMatrix(name)
{
   memset(&fTranslation[0], 0, kN3);
   memcpy(fRotationMatrix,kIdentityMatrix,kN9);
   memcpy(fScale,kUnitScale,kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// assignment

TGeoHMatrix::TGeoHMatrix(const TGeoMatrix &matrix)
            :TGeoMatrix(matrix)
{
   memset(&fTranslation[0], 0, kN3);
   memcpy(fRotationMatrix,kIdentityMatrix,kN9);
   memcpy(fScale,kUnitScale,kN3);
   if (matrix.IsIdentity()) return;
   if (matrix.IsTranslation())
      SetTranslation(matrix.GetTranslation());
   if (matrix.IsRotation())
      memcpy(fRotationMatrix,matrix.GetRotationMatrix(),kN9);
   if (matrix.IsScale())
      memcpy(fScale,matrix.GetScale(),kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoHMatrix::~TGeoHMatrix()
{
}

////////////////////////////////////////////////////////////////////////////////
/// assignment

TGeoHMatrix &TGeoHMatrix::operator=(const TGeoMatrix *matrix)
{
   return TGeoHMatrix::operator=(*matrix);
}

////////////////////////////////////////////////////////////////////////////////
/// assignment

TGeoHMatrix &TGeoHMatrix::operator=(const TGeoMatrix &matrix)
{
   if (&matrix == this) return *this;
   Clear();
   Bool_t registered = TestBit(kGeoRegistered);
   TNamed::operator=(matrix);
   if (matrix.IsIdentity()) return *this;
   if (matrix.IsTranslation())
      memcpy(fTranslation,matrix.GetTranslation(),kN3);
   if (matrix.IsRotation())
      memcpy(fRotationMatrix,matrix.GetRotationMatrix(),kN9);
   if (matrix.IsScale())
      memcpy(fScale,matrix.GetScale(),kN3);
   SetBit(kGeoRegistered,registered);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Composition

TGeoHMatrix &TGeoHMatrix::operator*=(const TGeoMatrix &right)
{
   Multiply(&right);
   return *this;
}

TGeoHMatrix TGeoHMatrix::operator*(const TGeoMatrix &right) const
{
   TGeoHMatrix h = *this;
   h *= right;
   return h;
}

////////////////////////////////////////////////////////////////////////////////
/// Is-equal operator

Bool_t TGeoHMatrix::operator==(const TGeoMatrix &other) const
{
   if (&other == this) return kTRUE;
   const Double_t *tr = GetTranslation();
   const Double_t *otr = other.GetTranslation();
   for (auto i=0; i<3; i++) if (TMath::Abs(tr[i]-otr[i])>1.E-10) return kFALSE;
   const Double_t *rot = GetRotationMatrix();
   const Double_t *orot = other.GetRotationMatrix();
   for (auto i=0; i<9; i++) if (TMath::Abs(rot[i]-orot[i])>1.E-10) return kFALSE;
   const Double_t *scl = GetScale();
   const Double_t *oscl = other.GetScale();
   for (auto i=0; i<3; i++) if (TMath::Abs(scl[i]-oscl[i])>1.E-10) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Fast copy method.

void TGeoHMatrix::CopyFrom(const TGeoMatrix *other)
{
   SetBit(kGeoTranslation, other->IsTranslation());
   SetBit(kGeoRotation, other->IsRotation());
   SetBit(kGeoReflection, other->IsReflection());
   memcpy(fTranslation,other->GetTranslation(),kN3);
   memcpy(fRotationMatrix,other->GetRotationMatrix(),kN9);
}

////////////////////////////////////////////////////////////////////////////////
/// clear the data for this matrix

void TGeoHMatrix::Clear(Option_t *)
{
   SetBit(kGeoReflection, kFALSE);
   if (IsIdentity()) return;
   ResetBit(kGeoTranslation);
   ResetBit(kGeoRotation);
   ResetBit(kGeoScale);
   memcpy(fTranslation,kNullVector,kN3);
   memcpy(fRotationMatrix,kIdentityMatrix,kN9);
   memcpy(fScale,kUnitScale,kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// Make a clone of this matrix.

TGeoMatrix *TGeoHMatrix::MakeClone() const
{
   TGeoMatrix *matrix = new TGeoHMatrix(*this);
   return matrix;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform a rotation about Z having the sine/cosine of the rotation angle.

void TGeoHMatrix::FastRotZ(const Double_t *sincos)
{
   fRotationMatrix[0] = sincos[1];
   fRotationMatrix[1] = -sincos[0];
   fRotationMatrix[3] = sincos[0];
   fRotationMatrix[4] = sincos[1];
   SetBit(kGeoRotation);
}

////////////////////////////////////////////////////////////////////////////////
/// Return a temporary inverse of this.

TGeoHMatrix TGeoHMatrix::Inverse() const
{
   TGeoHMatrix h;
   h = *this;
   if (IsTranslation()) {
      Double_t tr[3];
      tr[0] = -fTranslation[0]*fRotationMatrix[0] - fTranslation[1]*fRotationMatrix[3] - fTranslation[2]*fRotationMatrix[6];
      tr[1] = -fTranslation[0]*fRotationMatrix[1] - fTranslation[1]*fRotationMatrix[4] - fTranslation[2]*fRotationMatrix[7];
      tr[2] = -fTranslation[0]*fRotationMatrix[2] - fTranslation[1]*fRotationMatrix[5] - fTranslation[2]*fRotationMatrix[8];
      h.SetTranslation(tr);
   }
   if (IsRotation()) {
      Double_t newrot[9];
      newrot[0] = fRotationMatrix[0];
      newrot[1] = fRotationMatrix[3];
      newrot[2] = fRotationMatrix[6];
      newrot[3] = fRotationMatrix[1];
      newrot[4] = fRotationMatrix[4];
      newrot[5] = fRotationMatrix[7];
      newrot[6] = fRotationMatrix[2];
      newrot[7] = fRotationMatrix[5];
      newrot[8] = fRotationMatrix[8];
      h.SetRotation(newrot);
   }
   if (IsScale()) {
      Double_t sc[3];
      sc[0] = 1./fScale[0];
      sc[1] = 1./fScale[1];
      sc[2] = 1./fScale[2];
      h.SetScale(sc);
   }
   return h;
}

////////////////////////////////////////////////////////////////////////////////
/// computes determinant of the rotation matrix

Double_t TGeoHMatrix::Determinant() const
{
   Double_t
   det = fRotationMatrix[0]*fRotationMatrix[4]*fRotationMatrix[8] +
         fRotationMatrix[3]*fRotationMatrix[7]*fRotationMatrix[2] +
         fRotationMatrix[6]*fRotationMatrix[1]*fRotationMatrix[5] -
         fRotationMatrix[2]*fRotationMatrix[4]*fRotationMatrix[6] -
         fRotationMatrix[5]*fRotationMatrix[7]*fRotationMatrix[0] -
         fRotationMatrix[8]*fRotationMatrix[1]*fRotationMatrix[3];
   return det;
}

////////////////////////////////////////////////////////////////////////////////
/// multiply to the right with an other transformation
/// if right is identity matrix, just return

void TGeoHMatrix::Multiply(const TGeoMatrix *right)
{
   if (right->IsIdentity()) return;
   const Double_t *r_tra = right->GetTranslation();
   const Double_t *r_rot = right->GetRotationMatrix();
   const Double_t *r_scl = right->GetScale();
   if (IsIdentity()) {
      if (right->IsRotation()) {
         SetBit(kGeoRotation);
         memcpy(fRotationMatrix,r_rot,kN9);
         if (right->IsReflection()) SetBit(kGeoReflection, !TestBit(kGeoReflection));
      }
      if (right->IsScale()) {
         SetBit(kGeoScale);
         memcpy(fScale,r_scl,kN3);
      }
      if (right->IsTranslation()) {
         SetBit(kGeoTranslation);
         memcpy(fTranslation,r_tra,kN3);
      }
      return;
   }
   Int_t i, j;
   Double_t new_rot[9];

   if (right->IsRotation())    {
      SetBit(kGeoRotation);
      if (right->IsReflection()) SetBit(kGeoReflection, !TestBit(kGeoReflection));
   }
   if (right->IsScale())       SetBit(kGeoScale);
   if (right->IsTranslation()) SetBit(kGeoTranslation);

   // new translation
   if (IsTranslation()) {
      for (i=0; i<3; i++) {
         fTranslation[i] += fRotationMatrix[3*i]*r_tra[0]
                         + fRotationMatrix[3*i+1]*r_tra[1]
                         + fRotationMatrix[3*i+2]*r_tra[2];
      }
   }
   if (IsRotation()) {
      // new rotation
      for (i=0; i<3; i++) {
         for (j=0; j<3; j++) {
               new_rot[3*i+j] = fRotationMatrix[3*i]*r_rot[j] +
                                fRotationMatrix[3*i+1]*r_rot[3+j] +
                                fRotationMatrix[3*i+2]*r_rot[6+j];
         }
      }
      memcpy(fRotationMatrix,new_rot,kN9);
   }
   // new scale
   if (IsScale()) {
      for (i=0; i<3; i++) fScale[i] *= r_scl[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// multiply to the left with an other transformation
/// if right is identity matrix, just return

void TGeoHMatrix::MultiplyLeft(const TGeoMatrix *left)
{
   if (left == gGeoIdentity) return;
   const Double_t *l_tra = left->GetTranslation();
   const Double_t *l_rot = left->GetRotationMatrix();
   const Double_t *l_scl = left->GetScale();
   if (IsIdentity()) {
      if (left->IsRotation()) {
         if (left->IsReflection()) SetBit(kGeoReflection, !TestBit(kGeoReflection));
         SetBit(kGeoRotation);
         memcpy(fRotationMatrix,l_rot,kN9);
      }
      if (left->IsScale()) {
         SetBit(kGeoScale);
         memcpy(fScale,l_scl,kN3);
      }
      if (left->IsTranslation()) {
         SetBit(kGeoTranslation);
         memcpy(fTranslation,l_tra,kN3);
      }
      return;
   }
   Int_t i, j;
   Double_t new_tra[3];
   Double_t new_rot[9];

   if (left->IsRotation()) {
      SetBit(kGeoRotation);
      if (left->IsReflection()) SetBit(kGeoReflection, !TestBit(kGeoReflection));
   }
   if (left->IsScale())       SetBit(kGeoScale);
   if (left->IsTranslation()) SetBit(kGeoTranslation);

   // new translation
   if (IsTranslation()) {
      for (i=0; i<3; i++) {
         new_tra[i] = l_tra[i]
                      + l_rot[3*i]*  fTranslation[0]
                      + l_rot[3*i+1]*fTranslation[1]
                      + l_rot[3*i+2]*fTranslation[2];
      }
      memcpy(fTranslation,new_tra,kN3);
   }
   if (IsRotation()) {
      // new rotation
      for (i=0; i<3; i++) {
         for (j=0; j<3; j++) {
               new_rot[3*i+j] = l_rot[3*i]*fRotationMatrix[j] +
                                l_rot[3*i+1]*fRotationMatrix[3+j] +
                                l_rot[3*i+2]*fRotationMatrix[6+j];
         }
      }
      memcpy(fRotationMatrix,new_rot,kN9);
   }
   // new scale
   if (IsScale()) {
      for (i=0; i<3; i++) fScale[i] *= l_scl[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate about X axis with angle expressed in degrees.

void TGeoHMatrix::RotateX(Double_t angle)
{
   SetBit(kGeoRotation);
   Double_t phi = angle*TMath::DegToRad();
   Double_t c = TMath::Cos(phi);
   Double_t s = TMath::Sin(phi);
   Double_t v[9];
   v[0] = fRotationMatrix[0];
   v[1] = fRotationMatrix[1];
   v[2] = fRotationMatrix[2];
   v[3] = c*fRotationMatrix[3]-s*fRotationMatrix[6];
   v[4] = c*fRotationMatrix[4]-s*fRotationMatrix[7];
   v[5] = c*fRotationMatrix[5]-s*fRotationMatrix[8];
   v[6] = s*fRotationMatrix[3]+c*fRotationMatrix[6];
   v[7] = s*fRotationMatrix[4]+c*fRotationMatrix[7];
   v[8] = s*fRotationMatrix[5]+c*fRotationMatrix[8];
   memcpy(fRotationMatrix, v, kN9);

   v[0] = fTranslation[0];
   v[1] = c*fTranslation[1]-s*fTranslation[2];
   v[2] = s*fTranslation[1]+c*fTranslation[2];
   memcpy(fTranslation,v,kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate about Y axis with angle expressed in degrees.

void TGeoHMatrix::RotateY(Double_t angle)
{
   SetBit(kGeoRotation);
   Double_t phi = angle*TMath::DegToRad();
   Double_t c = TMath::Cos(phi);
   Double_t s = TMath::Sin(phi);
   Double_t v[9];
   v[0] = c*fRotationMatrix[0]+s*fRotationMatrix[6];
   v[1] = c*fRotationMatrix[1]+s*fRotationMatrix[7];
   v[2] = c*fRotationMatrix[2]+s*fRotationMatrix[8];
   v[3] = fRotationMatrix[3];
   v[4] = fRotationMatrix[4];
   v[5] = fRotationMatrix[5];
   v[6] = -s*fRotationMatrix[0]+c*fRotationMatrix[6];
   v[7] = -s*fRotationMatrix[1]+c*fRotationMatrix[7];
   v[8] = -s*fRotationMatrix[2]+c*fRotationMatrix[8];
   memcpy(fRotationMatrix, v, kN9);

   v[0] = c*fTranslation[0]+s*fTranslation[2];
   v[1] = fTranslation[1];
   v[2] = -s*fTranslation[0]+c*fTranslation[2];
   memcpy(fTranslation,v,kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate about Z axis with angle expressed in degrees.

void TGeoHMatrix::RotateZ(Double_t angle)
{
   SetBit(kGeoRotation);
   Double_t phi = angle*TMath::DegToRad();
   Double_t c = TMath::Cos(phi);
   Double_t s = TMath::Sin(phi);
   Double_t v[9];
   v[0] = c*fRotationMatrix[0]-s*fRotationMatrix[3];
   v[1] = c*fRotationMatrix[1]-s*fRotationMatrix[4];
   v[2] = c*fRotationMatrix[2]-s*fRotationMatrix[5];
   v[3] = s*fRotationMatrix[0]+c*fRotationMatrix[3];
   v[4] = s*fRotationMatrix[1]+c*fRotationMatrix[4];
   v[5] = s*fRotationMatrix[2]+c*fRotationMatrix[5];
   v[6] = fRotationMatrix[6];
   v[7] = fRotationMatrix[7];
   v[8] = fRotationMatrix[8];
   memcpy(&fRotationMatrix[0],v,kN9);

   v[0] = c*fTranslation[0]-s*fTranslation[1];
   v[1] = s*fTranslation[0]+c*fTranslation[1];
   v[2] = fTranslation[2];
   memcpy(fTranslation,v,kN3);
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply by a reflection respect to YZ.

void TGeoHMatrix::ReflectX(Bool_t leftside, Bool_t rotonly)
{
   if (leftside && !rotonly) fTranslation[0] = -fTranslation[0];
   if (leftside) {
      fRotationMatrix[0]=-fRotationMatrix[0];
      fRotationMatrix[1]=-fRotationMatrix[1];
      fRotationMatrix[2]=-fRotationMatrix[2];
   } else {
      fRotationMatrix[0]=-fRotationMatrix[0];
      fRotationMatrix[3]=-fRotationMatrix[3];
      fRotationMatrix[6]=-fRotationMatrix[6];
   }
   SetBit(kGeoRotation);
   SetBit(kGeoReflection, !IsReflection());
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply by a reflection respect to ZX.

void TGeoHMatrix::ReflectY(Bool_t leftside, Bool_t rotonly)
{
   if (leftside && !rotonly) fTranslation[1] = -fTranslation[1];
   if (leftside) {
      fRotationMatrix[3]=-fRotationMatrix[3];
      fRotationMatrix[4]=-fRotationMatrix[4];
      fRotationMatrix[5]=-fRotationMatrix[5];
   } else {
      fRotationMatrix[1]=-fRotationMatrix[1];
      fRotationMatrix[4]=-fRotationMatrix[4];
      fRotationMatrix[7]=-fRotationMatrix[7];
   }
   SetBit(kGeoRotation);
   SetBit(kGeoReflection, !IsReflection());
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply by a reflection respect to XY.

void TGeoHMatrix::ReflectZ(Bool_t leftside, Bool_t rotonly)
{
   if (leftside && !rotonly) fTranslation[2] = -fTranslation[2];
   if (leftside) {
      fRotationMatrix[6]=-fRotationMatrix[6];
      fRotationMatrix[7]=-fRotationMatrix[7];
      fRotationMatrix[8]=-fRotationMatrix[8];
   } else {
      fRotationMatrix[2]=-fRotationMatrix[2];
      fRotationMatrix[5]=-fRotationMatrix[5];
      fRotationMatrix[8]=-fRotationMatrix[8];
   }
   SetBit(kGeoRotation);
   SetBit(kGeoReflection, !IsReflection());
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoHMatrix::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TestBit(kGeoSavePrimitive)) return;
   const Double_t *tr = fTranslation;
   const Double_t *rot = fRotationMatrix;
   out << "   // HMatrix: " << GetName() << std::endl;
   out << "   tr[0]  = " << tr[0] << ";    " << "tr[1] = " << tr[1] << ";    " << "tr[2] = " << tr[2] << ";" << std::endl;
   out << "   rot[0] =" << rot[0] << ";    " << "rot[1] = " << rot[1] << ";    " << "rot[2] = " << rot[2] << ";" << std::endl;
   out << "   rot[3] =" << rot[3] << ";    " << "rot[4] = " << rot[4] << ";    " << "rot[5] = " << rot[5] << ";" << std::endl;
   out << "   rot[6] =" << rot[6] << ";    " << "rot[7] = " << rot[7] << ";    " << "rot[8] = " << rot[8] << ";" << std::endl;
   char *name = GetPointerName();
   out << "   TGeoHMatrix *" << name << " = new TGeoHMatrix(\"" << GetName() << "\");" << std::endl;
   out << "   " << name << "->SetTranslation(tr);" << std::endl;
   out << "   " << name << "->SetRotation(rot);" << std::endl;
   if (IsTranslation()) out << "   " << name << "->SetBit(TGeoMatrix::kGeoTranslation);" << std::endl;
   if (IsRotation()) out << "   " << name << "->SetBit(TGeoMatrix::kGeoRotation);" << std::endl;
   if (IsReflection()) out << "   " << name << "->SetBit(TGeoMatrix::kGeoReflection);" << std::endl;
   TObject::SetBit(kGeoSavePrimitive);
}
