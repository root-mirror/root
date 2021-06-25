// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

#ifndef ROOT_ROOT_Math_Vector3D
#define ROOT_ROOT_Math_Vector3D

// Defines typedefs to specific vectors and forward declarations.

namespace ROOT {

   namespace Math {

      template<class CoordSystem, class Tag> class DisplacementVector3D;

      template<typename T> class Cartesian3D;
      template<typename T> class CylindricalEta3D;
      template<typename T> class Polar3D;
      template<typename T> class Cylindrical3D;

      class DefaultCoordinateSystemTag;

      /// 3D Vector based on the cartesian coordinates x,y,z in double precision
      typedef DisplacementVector3D< Cartesian3D<double>, DefaultCoordinateSystemTag > XYZVector;

      /// 3D Vector based on the cartesian coordinates x,y,z in single precision
      typedef DisplacementVector3D< Cartesian3D<float>, DefaultCoordinateSystemTag > XYZVectorF;
      typedef XYZVector XYZVectorD;

      /// 3D Vector based on the eta based cylindrical coordinates rho, eta, phi in double precision.
      typedef DisplacementVector3D< CylindricalEta3D<double>, DefaultCoordinateSystemTag > RhoEtaPhiVector;

      /// 3D Vector based on the eta based cylindrical coordinates rho, eta, phi in single precision.
      typedef DisplacementVector3D< CylindricalEta3D<float>, DefaultCoordinateSystemTag > RhoEtaPhiVectorF;
      typedef RhoEtaPhiVector RhoEtaPhiVectorD;

      /// 3D Vector based on the polar coordinates rho, theta, phi in double precision.
      typedef DisplacementVector3D< Polar3D<double>, DefaultCoordinateSystemTag > Polar3DVector;

      /// 3D Vector based on the polar coordinates rho, theta, phi in single precision.
      typedef DisplacementVector3D< Polar3D<float>, DefaultCoordinateSystemTag > Polar3DVectorF;
      typedef Polar3DVector Polar3DVectorD;

      /// 3D Vector based on the cylindrical coordinates rho, z, phi in double precision.
      typedef DisplacementVector3D< Cylindrical3D<double>, DefaultCoordinateSystemTag > RhoZPhiVector;

      /// 3D Vector based on the cylindrical coordinates rho, z, phi in single precision.
      typedef DisplacementVector3D< Cylindrical3D<float>, DefaultCoordinateSystemTag > RhoZPhiVectorF;
      typedef RhoZPhiVector RhoZPhiVectorD;

   } // end namespace Math

} // end namespace ROOT

// Coordinate system types.
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/CylindricalEta3D.h"
#include "Math/GenVector/Polar3D.h"
#include "Math/GenVector/Cylindrical3D.h"

// Generic Vector3D class definition.
#include "Math/GenVector/DisplacementVector3D.h"

#endif
