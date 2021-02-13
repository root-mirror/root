\page ExtUsagePage Examples with External Packages

### Connection to Linear Algebra classes

It is possible to use the Vector and Rotation classes together with the Linear Algebra
classes. It is possible to set and get the contents of any 3D or 4D Vector from a Linear
Algebra Vector class which implements an iterator or something which behaves like an iterator.
For example a pointer to a C array (double *) behaves like an iterator. It is then assumed
that the coordinates, like (x,y,z) will be stored contiguously.

```
TVectorD       r2(N)                        // %ROOT Linear Algebra Vector containing many vectors
XYZVector      v2;
v2.SetCoordinates(&r2[INDEX],&r2[index]+3); // construct vector from x=r[INDEX], y=r[INDEX+1], z=r[INDEX+2]
```

To fill a Linear Algebra Vector from a 3D or 4D Vector, with GetCoordinates() one can get the internal coordinate data.

```
HepVector      c(3);                        // CLHEP Linear algebra vector
v2.GetCoordinates(&c[0],&c[index]+3 )       // fill HepVector c with c[0] = x, c[1] = y , c[2] = z
```

Or using TVectorD

```
double * data[3];
v2.GetCoordinates(data,data+3);
TVectorD       r1(3,data);                  // create a new Linear Algebra vector copying the data
```

In the case of transformations, constructor and method to set/get components exist with Linear Algebra matrices. The requisite is that the matrix data are stored, for example in the cse of a LorentzRotation, from (0,0) thru (3,3)

```
TMatrixD(4,4) m;
LorentzRotation r(m)                        // create LorentzRotation from matrix m
r.GetComponents(m)                          // fill matrix m with LorentzRotation components
```

### Connection to Other Vector classes

The 3D and 4D vectors of the GenVector package can be constructed and assigned from any Vector, which satisfies the following requisites:

*   for 3D Vectors and Points must implement the x(), y() ans z() methods
*   for LorentzVectors must implement the x(), y(), z() and t() methods.

```
CLHEP::Hep3Vector hv;
XYZVector  v1(hv);                          //  create  3D Vector from  CLHEP 3D Vector

HepGeom::Point3D <double>hp;
XYZPoint p1(hp);                            // create a 3D Point from CLHEP geom Point

CLHEP::HepLorentzVector  hq;
XYZTVector    q(hq);                        // create a LorentzVector  from CLHEP L.V.</double>
```
