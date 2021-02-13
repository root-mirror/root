\page Vector2DPage 2D Point and Vector Classes

Similar to the \ref Vector3DPage and \ref Point3DPage , typedefs are defined to avoid
exposing templated parameter to the users, for all 2D vectors based an double's and
float's. To use them, one must include the header file _Math/Vector2D.h_ or _Math/Point2D.h_.
The following typedef's, defined in the header file _Math/Vector2Dfwd.h_, are available for
the different instantiations of the template class ROOT::Math::DisplacementVector2D:

*   ROOT::Math::XYVector vector based on x,y coordinates (cartesian) in double precision
*   ROOT::Math::XYVectorF vector based on x,y coordinates (cartesian) in float precision
*   ROOT::Math::Polar2DVector vector based on r,phi coordinates (polar) in double precision
*   ROOT::Math::Polar2DVectorF vector based on r,phi coordinates (polar) in float precision

The typedef's, defined in the header file _Math/Point2Dfwd.h_, available for the different
instantiations of the template class ROOT::Math::PoistionVector2D are:

*   ROOT::Math::XYPoint vector based on x,y coordinates (cartesian) in double precision
*   ROOT::Math::XYPointF vector based on x,y coordinates (cartesian) in float precision
*   ROOT::Math::Polar2DPoint vector based on r,phi coordinates (polar) in double precision
*   ROOT::Math::Polar2DPointF vector based on r,phi coordinates (polar) in float precision

Similar constructs, functions and operations available for the 3D vectors and points
(see \ref Vector3DPage and \ref Point3DPage ) are available also for the 2D vector and
points. No transformations or rotation classes are available for the 2D vectors.
