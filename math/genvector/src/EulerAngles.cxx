// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2005 , LCG ROOT FNAL MathLib Team                    *
  *                                                                    *
  *                                                                    *
  **********************************************************************/

// Implementation file for rotation in 3 dimensions, represented by EulerAngles
//
// Created by: Mark Fischler Thurs June 9  2005
//
// Last update: $Id$
//
#include "Math/GenVector/EulerAngles.h"

#include <cmath>

#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/Quaternion.h"
#include "Math/GenVector/RotationX.h"
#include "Math/GenVector/RotationY.h"
#include "Math/GenVector/RotationZ.h"

#include "Math/GenVector/AxisAnglefwd.h"

namespace ROOT {

namespace Math {

// ========== Constructors and Assignment =====================

void EulerAngles::Rectify()
{
   // rectify
   if ( fTheta < 0 || fTheta > Pi() ) {
      Scalar t = fTheta - std::floor( fTheta/(2*Pi()) ) * 2*Pi();
      if ( t <= Pi() ) {
         fTheta = t;
      } else {
         fTheta = 2*Pi() - t;
         fPhi =  fPhi + Pi();
         fPsi =  fPsi + Pi();
      }
   }

   if ( fPhi <= -Pi()|| fPhi > Pi() ) {
      fPhi = fPhi - std::floor( fPhi/(2*Pi()) +.5 ) * 2*Pi();
   }

   if ( fPsi <= -Pi()|| fPsi > Pi() ) {
      fPsi = fPsi - std::floor( fPsi/(2*Pi()) +.5 ) * 2*Pi();
   }

} // Rectify()


// ========== Operations =====================

// DisplacementVector3D< Cartesian3D<double> >
// EulerAngles::
// operator() (const DisplacementVector3D< Cartesian3D<double> > & v) const
// {
//   return Rotation3D(*this)(v);
// }


EulerAngles EulerAngles::operator * (const Rotation3D  & r) const {
   // combine with a Rotation3D
   return EulerAngles ( Rotation3D(*this) * r );
}

EulerAngles EulerAngles::operator * (const AxisAngle   & a) const {
   // combine with a AxisAngle
   return EulerAngles ( Quaternion(*this) * Quaternion(a) );
}

EulerAngles EulerAngles::operator * (const EulerAngles & e) const {
   // combine with a EulerAngles
   return EulerAngles ( Quaternion(*this) * Quaternion(e) );
}
EulerAngles EulerAngles::operator * (const Quaternion & q) const {
   // combination with a Quaternion
   return EulerAngles ( Quaternion(*this) * q );
}

EulerAngles EulerAngles::operator * (const RotationX  & r) const {
   // combine with a RotationX
   return EulerAngles ( Quaternion(*this) * r );
}

EulerAngles EulerAngles::operator * (const RotationY  & r) const {
   // combine with a RotationY
   return EulerAngles ( Quaternion(*this) * r );
}

EulerAngles EulerAngles::operator * (const RotationZ  & r) const {
   // combine with a RotationZ
   // TODO -- this can be made much faster because it merely adds
   //         the r.Angle() to phi.
   Scalar newPhi = fPhi + r.Angle();
   if ( newPhi <= -Pi()|| newPhi > Pi() ) {
      newPhi = newPhi - std::floor( newPhi/(2*Pi()) +.5 ) * 2*Pi();
   }
   return EulerAngles ( newPhi, fTheta, fPsi );
}

EulerAngles operator * ( RotationX const & r, EulerAngles const & e )  {
   return EulerAngles(r) * e;  // TODO: improve performance
}

EulerAngles operator * ( RotationY const & r, EulerAngles const & e )  {
   return EulerAngles(r) * e;  // TODO: improve performance
}

EulerAngles
operator * ( RotationZ const & r, EulerAngles const & e )  {
   return EulerAngles(r) * e;  // TODO: improve performance
}

// ========== I/O =====================

std::ostream & operator<< (std::ostream & os, const EulerAngles & e) {
   // TODO - this will need changing for machine-readable issues
   //        and even the human readable form may need formatiing improvements
   os << "\n{phi: " << e.Phi() << "   theta: " << e.Theta()
   << "   psi: " << e.Psi() << "}\n";
   return os;
}


} //namespace Math
} //namespace ROOT
