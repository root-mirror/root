\page TransformPage Vector Transformations

Transformations classes are grouped in Rotations (in 3 dimensions), Lorentz transformations
and Poincaré transformations, which are Translation/Rotation combinations. Each group
has several members which may model physically equivalent trasformations but with different
internal representations. All the classes are non-template and use double precision as the
scalar type The following types of transformation classes are defined:

*   3D Rotations:
    *   ROOT::Math::Rotation3D, rotation described by a 3x3 matrix of doubles
    *   ROOT::Math::EulerAngles rotation described by the three Euler angles (phi, theta and psi) following the GoldStein [definition](http://mathworld.wolfram.com/EulerAngles.html).
    *   ROOT::Math::RotationZYX rotation described by three angles defining a rotation first along the Z axis, then along the rotated Y' axis and then along the rotated X'' axis.
    *   ROOT::Math::AxisAngle, rotation described by a vector (axis) and an angle
    *   ROOT::Math::Quaternion, rotation described by a quaternion (4 numbers)
    *   ROOT::Math::RotationX, specialized rotation along the X axis
    *   ROOT::Math::RotationY, specialized rotation along the Y axis
    *   ROOT::Math::RotationZ, specialized rotation along the Z axis
*   3D Transformations (Rotations + Translations)
    *   ROOT::Math::Translation3D, (only translation) described by a 3D Vector
    *   ROOT::Math::Transform3D, (rotations and then translation) described by a 3x4 matrix (12 numbers)
*   Lorentz Rotations and Boost
    *   ROOT::Math::LorentzRotation , 4D rotation (3D rotation plus a boost) described by a 4x4 matrix
    *   ROOT::Math::Boost, a Lorentz boost in an arbitrary direction and described by a 4x4 symmetric matrix (10 numbers)
    *   ROOT::Math::BoostX, a boost in the X axis direction
    *   ROOT::Math::BoostY, a boost in the Y axis direction
    *   ROOT::Math::BoostZ, a boost in the Z axis direction

#### Constructors

All rotations and transformations are default constructible (giving the identity transformation).
All rotations are constructible taking a number of scalar arguments matching the number
(and order of components)

```
Rotation3D    rI;                                // create a summy rotation (Identity matrix)
RotationX     rX(M_PI);                          // create a rotationX with an angle PI
EulerAngles   rE(phi, theta, psi);               // create a Euler rotation with phi,theta,psi angles
XYZVector     u(ux,uy,uz);
AxisAngle     rA(u, delta);                      // create a rotation based on direction u with delta angle
```

In addition, all rotations and transformations (other than the axial rotations) and
transformations are constructible from (begin,end) iterators or from pointers which behave
like iterators.

```
double          data[9];
Rotation3D      r(data, data+9);                 // create a rotation from a rotation matrix
std::vector <double>w(12);
Transform3D     t(w.begin(),w.end());            // create a Transform3D from the content of a std::vector</double>
```

All rotations, except the axial rotations, are constructible and assigned from any other
type of rotation (including the axial):

```
Rotation3D    r(ROOT::Math::RotationX(PI));      // create a rotation 3D from a rotation along X axis of angle PI
EulerAngles   r2(r);                             // construct an Euler Rotation from A Rotation3D
AxisAngle     r3; r3 = r2;                       // assign an Axis Rotation from an Euler Rotation;
```

Transform3D (rotation + translation) can be constructed from a rotation and a translation vector

```
Rotation3D r; XYZVector v;
Transform3D   t1(r,v);                           // construct from rotation and then translation
Transform3D   t2(v,r);                           // construct inverse from first translation then rotation
Transform3D   t3(r);                             // construct from only a rotation (zero translation)
Transform3D   t4(v);                             // construct from only translation (identity rotation)
```

#### Operations

All transformations can be applied to vector and points using the _operator *_ or using the _operator()_

```
XYZVector  v1(...);
Rotation3D r(...);
XYZVector v2 = r*v1;                             // rotate vector v1 using r
v2 = r(v1)                                       // equivalent
```

Transformations can be combined using the operator * . Note that the rotations are not
commutative ans therefore the order is important

```
Rotation3D     r1(...);
Rotation3D     r2(...);
Rotation3D  r3 = r2*r1;                          // obtain a combine rotation r3 by applying first r1 then r2
```

We can combine rotations of different types, like Rotation3D with any other type of rotations.
The product of two different axial rotations return a Rotation3D:

```
RotationX        rx(1.);
RotationY        ry(2.);
Rotation3D  r = ry * rx;                         // rotation along X and then Y axis
```

It is also possible to invert all the transformation or return the inverse of a transformation

```
Rotation3D           r1(...);
                 r1.Invert();                    // invert the rotation modifying its content
Rotation3D  r2 =r1.Inverse();                    // return the inverse in a new rotation class
```

We have used rotation as examples, but all these operations can be applied to all the
transformation classes. Rotation3D, Transform3D and Translation3D classes can all be combined
via the _operator *_.

```
Rotation3D     r(AxisAngle(phi,ux,uy,uz));   // rotation of an angle phi around u.
Translation3D  d(dx,dy,dz);                  // translation of a vector d
Transform3D    t1 = d * r;                   // transformation obtained applying first the rotation
Transform3D    t2 = r * d;                   // transformation obtained applying first the translation
```

#### Set/GetComponents methods

Common methods to all the transformations are the Get and SetComponents. They can be used
to retrieve all the scalar values on which the trasformation is based. They can be used with
a signature based iterators or by using any foreign matrix which implements the _operator(i,j)_
or a different signatures depending on the transformation type.

```
RotationX  rx;  rx.SetComponents(1.)           // set agle of the X rotation
double d[9] = {........}
Rotation3D r;   r.SetComponents(d,d+9);        // set 9 components of 3D rotation
double d[16];
LorentzRotation lr;
lr.GetComponents( d, d+16);                    // get 16 components of a LorentzRotation
TMatrixD(3,4) m;
Transform3D t;  t.GetComponens(m);             // fill matrix of size 3x4 with components of the transform3D t
```

For more detailed documentation on all methods see the reference doc for the specific
transformation class.
