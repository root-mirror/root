\page LorentzVectorPage LorentzVector Classes


To avoid exposing templated parameter to the users, typedefs are defined for all types of
vectors based an double's and float's. To use them, one must include the header
file _Math/Vector4D.h_. The following typedef's, defined in the header file
_Math/Vector4Dfwd.h_, are available for the different instantiations of the template class
ROOT::Math::LorentzVector:

*   ROOT::Math::XYZTVector vector based on x,y,z,t coordinates (cartesian) in double precision
*   ROOT::Math::XYZTVectorF vector based on x,y,z,t coordinates (cartesian) in float precision
*   ROOT::Math::PtEtaPhiEVector vector based on pt (rho),eta,phi and E (t) coordinates in double precision
*   ROOT::Math::PtEtaPhiMVector vector based on pt (rho),eta,phi and M (t) coordinates in double precision
*   ROOT::Math::PxPyPzMVector vector based on px,py,pz and M (mass) coordinates in double precision
*   ROOT::Math::PxPyPzEVector vector based on px,py,pz and E (energy) coordinates in double precision

The metric used for all the LorentzVector's is (-,-,-,+)

#### Constructors and Assignment

The following declarations are available:


```
XYZTVector               v1;               // create an empty vector (x = 0, y = 0, z = 0, t = 0)
XYZTVector               v2(1,2,3,4);      // create a vector with x=1, y = 2, z = 3, t = 4
PtEtaPhiEVector          v3(1,2,PI,5);     // create a vector with pt = 1, eta = 2, phi = PI, E = 5
```

Note that each type of vector is constructed by passing its coordinates representations,
so a XYZTVector(1,2,3,4) is different from a PtEtaPhiEVector(1,2,3,4).

In addition the Vector classes can be constructed by any vector, which implements the
accessors x(), y() and z() and t(). This cann be another ROOT::Math::LorentzVector based
on a different coordinate system or even any vector of a different package, like the CLHEP
HepLorentzVector that implements the required signature.

```
XYZTVector               v1(1,2,3,4);
PtEtaPhiEVector          v2(v1);

CLHEP::HepLorentzVector  q(1,2,3,4);
XYZTVector               v3(q)
```

#### Coordinate Accessors

All the same coordinate accessors are available through the interface of the class
ROOT::Math::LorentzVector. For example:

```
v1.X(); v1.X(); v1.Z(); v1.T()             // returns cartesian components for the cartesian vector v1
v2.Px(); v2.Py(); v2.Pz(); v2.E()          // returns cartesian components for the cylindrical vector v2
v1.Pt(); v1.Eta(); v1.Phi(); v1.M()        // returns other components for the cartesian vector v1
```

In addition, all the 4 coordinates of the vector can be retrieved with the GetCoordinates method:

```
double d[4];
v1.GetCoordinates(d);                      // fill d array with (x,y,z,t) components of v1
v2.GetCoordinates(d);                      // fill d array with (pt,eta,phi,e) components of v2
std::vector w(4);
v1.GetCoordinates(w.begin(),w.end());      // fill std::vector with (x,y,z,t) components of v1
```

To get information on all the coordinate accessors see the reference documentation of
ROOT::Math::LorentzVector

#### Setter Methods

One can set only all the three coordinates via:

```
v1.SetCoordinates(c1,c2,c3,c4);            // sets the (x,y,z,t) for a XYZTVector
v2.SetCoordinates(c1,c2,c3,c4);            // sets pt,eta,phi,e for a PtEtaPhiEVector
v2.SetXYZ(x,y,z,t);                        // sets the 4 cartesian components for the PtEtaPhiEVector
```

Single coordinate setter methods are available for the basic vector coordinates, like SetX()
for a XYZTVector or SetPt() for a PtEtaPhiEVector. Attempting to do a SetX() on a non
cartesian vector will not compile.

```
XYZTVector v1;      v1.SetX(1)             // OK setting x for a cartesian vector
PtEtaPhiEVector v2; v2.SetX(1)             // ERROR: cannot set X for a non-cartesian vector. Method will not compile
v2.SetR(1)                                 // OK setting Pt for a  PtEtaPhiEVector vector
```

In addition there are setter methods from C arrays or iterators.

```
double d[4] = {1.,2.,3.,4.};
XYZTVector v;
v.SetCoordinates(d);                      // set (x,y,z,t) components of v using values from d
```

or for example from an std::vector using the iterators

```
std::vector w(4);
v.SetCoordinates(w.begin(),w.end());      // set (x,y,z,t) components of v using values from w
```

#### Arithmetic Operations

The following operations are possible between LorentzVectors classes, even of different
coordinate system types: ( v and w are two LorentzVector of the same type, q is a generic
LorentzVector implementing x(), y(), z() and t() and a is a generic scalar
type: double, flot, int, etc.... )

```
v += q;
v -= q;
v  = -q;
v *= a;
v /= a;
w = v + q;
w = v - q;
w = v * a;
w = a * v;
w = v / a;
```

#### Comparison

```
v == w;
v != w;
```

#### Other Methods

```
a =  v.Dot(q);                           // dot product in metric (+,+,+,-) of two LorentzVector's
XYZVector s = v.Vect()                   // return the spatial components (x,y,z)
v.Beta();                                // return beta and gamma value
v.Gamma()                                // (vector must be time-like otherwise result is meaningless)
XYZVector b = v.BoostToCM()              // return boost vector which will bring the Vector in its mas frame (P=0)
```
