// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2005 , LCG ROOT MathLib Team                         *
  *                                                                    *
  *                                                                    *
  **********************************************************************/

// Header file Vector2Dfwd
//
// Created by: Lorenzo Moneta  at Mon Apr 16 2007
//
#ifndef ROOT_Math_Vector2Dfwd
#define ROOT_Math_Vector2Dfwd  1

// forward declarations of displacement vectors (Vectors) and type defs definitions

namespace ROOT {

   namespace Math {

      template<class CoordSystem, class Tag> class DisplacementVector2D;

      template<typename T> class Cartesian2D;
      template<typename T> class Polar2D;

      class DefaultCoordinateSystemTag;

      typedef DisplacementVector2D< Cartesian2D<double>, DefaultCoordinateSystemTag > XYVector;
      typedef XYVector XYVectorD;
      typedef DisplacementVector2D< Cartesian2D<float>, DefaultCoordinateSystemTag > XYVectorF;
      typedef DisplacementVector2D< Polar2D<double>, DefaultCoordinateSystemTag > Polar2DVector;
      typedef Polar2DVector Polar2DVectorD;
      typedef DisplacementVector2D< Polar2D<float>, DefaultCoordinateSystemTag > Polar2DVectorF;

   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_Vector2Dfwd  */
