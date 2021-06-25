// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

#ifndef ROOT_ROOT_Math_Vector2D
#define ROOT_ROOT_Math_Vector2D

// Defines typedefs to specific vectors and forward declarations.
//
namespace ROOT {

   namespace Math {

      template<class CoordSystem, class Tag> class DisplacementVector2D;

      template<typename T> class Cartesian2D;
      template<typename T> class Polar2D;

      class DefaultCoordinateSystemTag;

      /// 2D Vector based on the cartesian coordinates x,y in double precision
      typedef DisplacementVector2D< Cartesian2D<double>, DefaultCoordinateSystemTag > XYVector;
      typedef XYVector XYVectorD;

      /// 2D Vector based on the cartesian coordinates x,y,z in single precision
      typedef DisplacementVector2D< Cartesian2D<float>, DefaultCoordinateSystemTag > XYVectorF;

      /// 2D Vector based on the polar coordinates rho, phi in double precision.
      typedef DisplacementVector2D< Polar2D<double>, DefaultCoordinateSystemTag > Polar2DVector;
      typedef Polar2DVector Polar2DVectorD;

      /// 2D Vector based on the polar coordinates rho, phi in single precision.
      typedef DisplacementVector2D< Polar2D<float>, DefaultCoordinateSystemTag > Polar2DVectorF;

   } // end namespace Math

} // end namespace ROOT


// Coordinate system types.
#include "Math/GenVector/Cartesian2D.h"
#include "Math/GenVector/Polar2D.h"

// Generic Vector2D class definition.
#include "Math/GenVector/DisplacementVector2D.h"

#endif
