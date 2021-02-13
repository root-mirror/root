  \page Vector3DPage  Vector3D Classes


To avoid exposing templated parameter to the users, typedefs are defined for all types of
vectors based an double's and float's. To use them, one must include the header file _Math/Vector3D.h_.
The following typedef's, defined in the header file _Math/Vector3Dfwd.h_, are available
for the different instantiations of the template class ROOT::Math::DisplacementVector3D:

*   ROOT::Math::XYZVector vector based on x,y,z coordinates (cartesian) in double precision
*   ROOT::Math::XYZVectorF vector based on x,y,z coordinates (cartesian) in float precision
*   ROOT::Math::Polar3DVector vector based on r,theta,phi coordinates (polar) in double precision
*   ROOT::Math::Polar3DVectorF vector based on r,theta,phi coordinates (polar) in float precision
*   ROOT::Math::RhoZPhiVector vector based on rho, z,phi coordinates (cylindrical) in double precision
*   ROOT::Math::RhoZPhiVectorF vector based on rho, z,phi coordinates (cylindrical) in float precision
*   ROOT::Math::RhoEtaPhiVector vector based on rho,eta,phi coordinates (cylindrical using eta instead of z) in double precision
*   ROOT::Math::RhoEtaPhiVectorF vector based on rho,eta,phi coordinates (cylindrical using eta instead of z) in float precision

#### Constructors and Assignment

The following declarations are available:

```
 XYZVector         v1;                     // create an empty vector (x = 0, y = 0, z = 0)
 XYZVector         v2( 1,2,3);             // create a vector with x=1, y = 2, z = 3;
 Polar3DVector     v3( 1, PI/2, PI);       // create a vector with r = 1, theta = PI/2 and phi=PI
 RhoEtaPHiVector   v4( 1, 2, PI)           // create a vector with rho= 1, eta = 2, phi = PI
```

Note that each type of vector is constructed by passing its coordinates representations, so
a XYZVector(1,2,3) is different from a Polar3DVector(1,2,3).

In addition the Vector classes can be constructed by any vector, which implements the
accessors x(), y() and z(). This con be another Vector3D based on a different coordinate
system types or even any vector of a different package, like the CLHEP HepThreeVector that
implements the required signatures.

```
  XYZVector v1(1,2,3);
  RhoEtaPhiVector r2(v1);
  CLHEP::HepThreeVector q(1,2,3);
  XYZVector v3(q)
```

#### Coordinate Accessors

All the same coordinate accessors are available through the interface of the
class ROOT::Math::DisplacementVector3D. For example:

```
v1.X(); v1.X(); v1.Z()                     // returns cartesian components for the cartesian vector v1
v1.Rho(); v1.Eta(); v1.Phi()               // returns cylindrical components for the cartesian vector v1
r2.X(); r2.Y(); r2.Z()                     // returns cartesian components for the cylindrical vector r2
```

In addition, all the 3 coordinates of the vector can be retrieved with the GetCoordinates method:

```
double d[3];
v1.GetCoordinates(d);                     // fill d array with (x,y,z) components of v1
r2.GetCoordinates(d);                     // fill d array with (r,eta,phi) components of r2
std::vector <double>vc(3);
v1.GetCoordinates(vc.begin(),vc.end());   // fill std::vector with (x,y,z) components of v1</double>
```

To get more information on all the coordinate accessors see the reference documentation of
ROOT::Math::DisplacementVector3D.

#### Setter Methods

One can set only all the three coordinates via:

```
v1.SetCoordinates(c1,c2,c3);               // sets the (x,y,z) for a XYZVector
r2.SetCoordinates(c1,c2,c3);               // sets r,theta,phi for a Polar3DVector
r2.SetXYZ(x,y,z);                          // sets the three cartesian components for the Polar3DVector
```

Single coordinate setter methods are available for the basic vector coordinates, like SetX()
for a XYZVector or SetR() for a polar vector. Attempting to do a SetX() on a polar vector
will not compile.

```
XYZVector v1;      v1.SetX(1)             // OK setting x for a Cartesian vector
Polar3DVector v2;  v2.SetX(1)             // ERROR: cannot set  X for a Polar vector. Method will not compile
v2.SetR(1)                                // OK setting r for a Polar vector
```

In addition there are setter methods from C arrays or iterators.

```
double d[3] = {1.,2.,3.};
XYZVector v;
v.SetCoordinates(d);                      // set (x,y,z) components of v using values from d
```

or for example from an std::vector using the iterator

```
std::vector <double>w(3);
v.SetCoordinates(w.begin(),w.end());      // set (x,y,z) components of v using values from w</double>
```

#### Arithmetic Operations

The following operations are possible between Vector classes, even of different coordinate
system types: ( v1,v2 are any type of ROOT::Math::DisplacementVector3D classes, v3 is the
same type of v1; _a_ is a scalar value)

```
v1 += v2;
v1 -= v2;
v1 = - v2;
v1 *= a;
v1 /= a;
v2 = a * v1;
v2 = v1 / a;
v2 = v1 * a;
v3 = v1 + v2;
v3 = v1 - v2;
```

#### Comparison

For v1 and v2 of the same type (same coordinate system and same scalar type):

```
v1 == v2;
v1 != v2;
```

#### Dot and Cross Product

We support the dot and cross products, through the Dot() and Cross() method, with any
Vector (q) implementing x(), y() and z()

```
XYZVector v1(x,y,z);
double s = v1.Dot(q);
XYZVector v2 = v1.Cross(q);
```

Note that the multiplication between two vectors using the operator * is not supported
because is ambiguous.

#### Other Methods

```
XYZVector u = v1.Unit();               //  return unit vector parallel to v1
```
