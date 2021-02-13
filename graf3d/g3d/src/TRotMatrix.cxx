// @(#)root/g3d:$Id$
// Author: Rene Brun   14/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TRotMatrix.h"
#include "TBuffer.h"
#include "TGeometry.h"
#include "TMath.h"

ClassImp(TRotMatrix);

/** \class TRotMatrix
\ingroup g3d
Manages a detector rotation matrix. See class TGeometry.
*/

////////////////////////////////////////////////////////////////////////////////
/// RotMatrix default constructor.

TRotMatrix::TRotMatrix()
{
   for (int i=0;i<9;i++) fMatrix[i] = 0;
   fNumber = 0;
   fPhi    = 0;
   fPsi    = 0;
   fTheta  = 0;
   fType   = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// RotMatrix normal constructor.

TRotMatrix::TRotMatrix(const char *name, const char *title, Double_t *matrix)
           :TNamed(name,title)
{
   fNumber = 0;
   fPhi    = 0;
   fPsi    = 0;
   fTheta  = 0;
   fType   = 0;

   if (!matrix) { Error("ctor","No rotation is supplied"); return; }

   SetMatrix(matrix);
   if (!gGeometry) gGeometry = new TGeometry();
   fNumber = gGeometry->GetListOfMatrices()->GetSize();
   gGeometry->GetListOfMatrices()->Add(this);
}

////////////////////////////////////////////////////////////////////////////////
/// RotMatrix normal constructor.

TRotMatrix::TRotMatrix(const char *name, const char *title, Double_t theta, Double_t phi, Double_t psi)
           :TNamed(name,title)
{
   printf("ERROR: This form of TRotMatrix constructor not implemented yet\n");

   Int_t i;
   fTheta  = theta;
   fPhi    = phi;
   fPsi    = psi;
   fType   = 2;
   for (i=0;i<9;i++) fMatrix[i] = 0;
   fMatrix[0] = 1;   fMatrix[4] = 1;   fMatrix[8] = 1;

   if (!gGeometry) gGeometry = new TGeometry();
   fNumber = gGeometry->GetListOfMatrices()->GetSize();
   gGeometry->GetListOfMatrices()->Add(this);
}

////////////////////////////////////////////////////////////////////////////////
/// RotMatrix normal constructor defined a la GEANT.
///
/// The TRotMatrix constructor with six angles uses the GEANT convention:
///
/// theta1 is the polar angle of the x-prim axis in the main reference system
/// (MRS), theta2 and theta3 have the same meaning for the y-prim and z-prim
/// axis.
///
/// Phi1 is the azimuthal angle of the x-prim in the MRS and phi2 and phi3
/// have the same meaning for y-prim and z-prim.
///
///
/// for example, the unit matrix is defined in the following way.
/// ~~~ {.cpp}
///     x-prim || x, y-prim || y, z-prim || z
///
///     means:  theta1=90, theta2=90, theta3=0, phi1=0, phi2=90, phi3=0
/// ~~~

TRotMatrix::TRotMatrix(const char *name, const char *title, Double_t theta1, Double_t phi1
                                                  , Double_t theta2, Double_t phi2
                                                  , Double_t theta3, Double_t phi3)
                :TNamed(name,title)
{
   SetAngles(theta1,phi1,theta2,phi2,theta3,phi3);

   if (!gGeometry) gGeometry = new TGeometry();
   fNumber = gGeometry->GetListOfMatrices()->GetSize();
   gGeometry->GetListOfMatrices()->Add(this);
}

////////////////////////////////////////////////////////////////////////////////
/// RotMatrix default destructor.

TRotMatrix::~TRotMatrix()
{
   if (gGeometry) gGeometry->GetListOfMatrices()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the value of the determinant of this matrix

Double_t  TRotMatrix::Determinant() const
{
   return
      fMatrix[0] * (fMatrix[4]*fMatrix[8] - fMatrix[7]*fMatrix[5])
    - fMatrix[3] * (fMatrix[1]*fMatrix[8] - fMatrix[7]*fMatrix[2])
    + fMatrix[6] * (fMatrix[1]*fMatrix[5] - fMatrix[4]*fMatrix[2]);
}

////////////////////////////////////////////////////////////////////////////////
///  Convert this matrix to the OpenGL [4x4]
///
/// ~~~ {.cpp}
///  [  fMatrix[0]   fMatrix[1]   fMatrix[2]    0  ]
///  [  fMatrix[3]   fMatrix[4]   fMatrix[5]    0  ]
///  [  fMatrix[6]   fMatrix[7]   fMatrix[8]    0  ]
///  [     0             0           0          1  ]
/// ~~~
///
///  Input:
///
///  Double_t *rGLMatrix: pointer to Double_t 4x4 buffer array
///
///  Return:
///
///  Double_t*: pointer to the input buffer

Double_t* TRotMatrix::GetGLMatrix(Double_t *rGLMatrix) const
{
   Double_t *glmatrix = rGLMatrix;
   const Double_t *matrix   = fMatrix;
   if (rGLMatrix)
   {
      for (Int_t i=0;i<3;i++) {
         for (Int_t j=0;j<3;j++) memcpy(glmatrix,matrix,3*sizeof(Double_t));
         matrix   += 3;
         glmatrix += 3;
         *glmatrix = 0.0;
         glmatrix++;
      }
      for (Int_t j=0;j<3;j++) {
         *glmatrix = 0.0;
         glmatrix++;
      }
      *glmatrix = 1.0;
   }
   return rGLMatrix;
}

////////////////////////////////////////////////////////////////////////////////
/// theta1 is the polar angle of the x-prim axis in the main reference system
/// (MRS), theta2 and theta3 have the same meaning for the y-prim and z-prim
/// axis.
///
/// Phi1 is the azimuthal angle of the x-prim in the MRS and phi2 and phi3
/// have the same meaning for y-prim and z-prim.
///
///
/// for example, the unit matrix is defined in the following way.
///
/// ~~~ {.cpp}
///     x-prim || x, y-prim || y, z-prim || z
///
///     means:  theta1=90, theta2=90, theta3=0, phi1=0, phi2=90, phi3=0
/// ~~~

const Double_t* TRotMatrix::SetAngles(Double_t theta1, Double_t phi1,
                Double_t theta2, Double_t phi2,Double_t theta3, Double_t phi3)
{
   const Double_t degrad = 0.0174532925199432958;

   fTheta  = theta1;
   fPhi    = phi1;
   fPsi    = theta2;

   fType   = 2;
   if (!strcmp(GetName(),"Identity")) fType = 0;

   fMatrix[0] = TMath::Sin(theta1*degrad)*TMath::Cos(phi1*degrad);
   fMatrix[1] = TMath::Sin(theta1*degrad)*TMath::Sin(phi1*degrad);
   fMatrix[2] = TMath::Cos(theta1*degrad);
   fMatrix[3] = TMath::Sin(theta2*degrad)*TMath::Cos(phi2*degrad);
   fMatrix[4] = TMath::Sin(theta2*degrad)*TMath::Sin(phi2*degrad);
   fMatrix[5] = TMath::Cos(theta2*degrad);
   fMatrix[6] = TMath::Sin(theta3*degrad)*TMath::Cos(phi3*degrad);
   fMatrix[7] = TMath::Sin(theta3*degrad)*TMath::Sin(phi3*degrad);
   fMatrix[8] = TMath::Cos(theta3*degrad);

   SetReflection();
   return fMatrix;
}

////////////////////////////////////////////////////////////////////////////////
/// copy predefined 3x3 matrix into TRotMatrix object

void TRotMatrix::SetMatrix(const Double_t *matrix)
{
   fTheta  = 0;
   fPhi    = 0;
   fPsi    = 0;
   fType   = 0;
   if (!matrix) return;
   fType   = 2;
   memcpy(fMatrix,matrix,9*sizeof(Double_t));
   SetReflection();
}

////////////////////////////////////////////////////////////////////////////////
/// Checks whether the determinant of this
/// matrix defines the reflection transformation
/// and set the "reflection" flag if any

void TRotMatrix::SetReflection()
{
   ResetBit(kReflection);
   if (Determinant() < 0) { fType=1; SetBit(kReflection);}
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TRotMatrix.

void TRotMatrix::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         R__b.ReadClassBuffer(TRotMatrix::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TNamed::Streamer(R__b);
      R__b >> fNumber;
      R__b >> fType;
      R__b >> fTheta;
      R__b >> fPhi;
      R__b >> fPsi;
      R__b.ReadStaticArray(fMatrix);
      R__b.CheckByteCount(R__s, R__c, TRotMatrix::IsA());
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TRotMatrix::Class(),this);
   }
}
