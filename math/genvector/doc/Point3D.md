\page Point3DPage Point3D Classes

To avoid exposing templated parameter to the users, typedefs are defined for all types of
vectors based an double's and float's. To use them, one must include the header file _Math/Point3D.h_.
The following typedef's, defined in the header file _Math/Point3Dfwd.h_, are available for
the different instantiations of the template class ROOT::Math::PositionVector3D:

*   ROOT::Math::XYZPoint point based on x,y,z coordinates (cartesian) in double precision
*   ROOT::Math::XYZPointF point based on x,y,z coordinates (cartesian) in float precision
*   ROOT::Math::Polar3DPoint point based on r,theta,phi coordinates (polar) in double precision
*   ROOT::Math::Polar3DPointF point based on r,theta,phi coordinates (polar) in float precision
*   ROOT::Math::RhoZPhiPoint point based on rho,z,phi coordinates (cylindrical using z) in double precision
*   ROOT::Math::RhoZPhiPointF point based on rho,z,phi coordinates (cylindrical using z) in float precision
*   ROOT::Math::RhoEtaPhiPoint point based on rho,eta,phi coordinates (cylindrical using eta instead of z) in double precision
*   ROOT::Math::RhoEtaPhiPointF point based on rho,eta,phi coordinates (cylindrical using eta instead of z) in float precision

#### Constructors and Assignment

The following declarations are available:

```
 XYZPoint         p1;                     // create an empty vector (x = 0, y = 0, z = 0)
 XYZPoint         p2( 1,2,3);             // create a vector with x=1, y = 2, z = 3;
 Polar3DPoint     p3( 1, PI/2, PI);       // create a vector with r = 1, theta = PI/2 and phi=PI
 RhoEtaPHiPoint   p4( 1, 2, PI)           // create a vector with rho= 1, eta = 2, phi = PI
```

Note that each type of vector is constructed by passing its coordinates representations, so
a XYZPoint(1,2,3) is different from a Polar3DPoint(1,2,3).

In addition the Point classes can be constructed by any vector, which implements the accessors
x(), y() and z(). This con be another Point3D based on a different coordinate system types
or even any vector of a different package, like the CLHEP HepThreePoint that implements
the required signatures.

```
  XYZPoint             p1(1,2,3);
  RhoEtaPHiPoint       r2(v1);
  CLHEP::HepThreePoint q(1,2,3);
  XYZPoint             p3(q)
```

#### Coordinate Accessors and Setter Methods

For the Points classes we have the same getter and setter methods as for the Vector classes.
See the examples for the \ref Vector3DPage.

#### Point-Vector Operations

The following operations are possible between points and vector classes: ( p1 ,p2 and p3 are
instantiations of the ROOT::Math::PositionVector3D class, p1 and p3 of the same type v1 and
v2 are a ROOT::Math::DisplacementVector3D class )

```
p1 += v1;
p1 -= v1;
p3 = p1 + v1;      // p1 and p3 are the same type
p3 = v1 + p1;      // p3 is based on the same coordinate system as v1
p3 = p1 - v1;
p3 = v1 - p1;
v2 = p1 - p2;    // difference between points returns a vector v2 based on the same coordinate system as p1
```

Note that additions between two points is NOT possible and the difference between points returns a vector.

#### Other Operations

Exactly as for the 3D Vectors, the following operations are allowed:

*   comparison of points
*   scaling and division of points with a scalar
*   dot and cross product with any type of vector
