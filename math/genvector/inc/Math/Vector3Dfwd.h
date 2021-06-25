// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2005 , LCG ROOT MathLib Team                         *
  *                                                                    *
  *                                                                    *
  **********************************************************************/

// Header file Vector3Dfwd
//
// Created by: Lorenzo Moneta  at Mon May 30 18:08:35 2005
//
// Last update: Mon May 30 18:08:35 2005
//
#ifndef ROOT_Math_Vector3Dfwd
#define ROOT_Math_Vector3Dfwd  1

// forward declarations of displacement vectors (Vectors) and type defs definitions

namespace ROOT {

   namespace Math {

      template<class CoordSystem, class Tag> class DisplacementVector3D;

      template<typename T> class Cartesian3D;
      template<typename T> class CylindricalEta3D;
      template<typename T> class Polar3D;
      template<typename T> class Cylindrical3D;

      class DefaultCoordinateSystemTag;

      typedef DisplacementVector3D< Cartesian3D<double>, DefaultCoordinateSystemTag > XYZVector;
      typedef DisplacementVector3D< Cartesian3D<float>, DefaultCoordinateSystemTag > XYZVectorF;
      typedef XYZVector XYZVectorD;
      typedef DisplacementVector3D< CylindricalEta3D<double>, DefaultCoordinateSystemTag > RhoEtaPhiVector;
      typedef DisplacementVector3D< CylindricalEta3D<float>, DefaultCoordinateSystemTag > RhoEtaPhiVectorF;
      typedef RhoEtaPhiVector RhoEtaPhiVectorD;
      typedef DisplacementVector3D< Polar3D<double>, DefaultCoordinateSystemTag > Polar3DVector;
      typedef DisplacementVector3D< Polar3D<float>, DefaultCoordinateSystemTag > Polar3DVectorF;
      typedef Polar3DVector Polar3DVectorD;
      typedef DisplacementVector3D< Cylindrical3D<double>, DefaultCoordinateSystemTag > RhoZPhiVector;
      typedef DisplacementVector3D< Cylindrical3D<float>, DefaultCoordinateSystemTag > RhoZPhiVectorF;
      typedef RhoZPhiVector RhoZPhiVectorD;

   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_Vector3Dfwd  */
