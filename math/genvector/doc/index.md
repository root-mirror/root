\defgroup GenVector GenVector
\ingroup Math
\brief Generic 2D, 3D and 4D vectors classes and their transformations (rotations).

**GenVector**, is a new package intended to represent vectors and their operations and
transformations, such as rotations and Lorentz transformations, in 2, 3 and 4 dimensions.
The 2D and 3D space are used to describe the geometry vectors and points, while the 4D
space-time is used for physics vectors representing relativistic particles.
These 2D,3D and 4D vectors are different from vectors of the Linear Algebra package which
describe generic N-dimensional vectors. Similar functionality is currently provided by the
CLHEP [Vector](http://proj-clhep.web.cern.ch/proj-clhep/manual/UserGuide/VectorDefs/index.html)
and [Geometry](http://www.hep.phy.cam.ac.uk/lhcb/doc/CLHEP/1.9.1.2/html/namespaceHepGeom.html)
packages and the %ROOT [Physics Vector](http://root.cern.ch/root/html/PHYSICS_Index.html)
classes (TVector2, TVector3 and TLorentzVector). It is also re-uses concepts and ideas from
the CMS [Common Vector package](http://lcgapp.cern.ch/doxygen/SEAL/snapshot/html/dir_000007.html).

In contrast to CLHEP or the %ROOT physics libraries, GenVector provides class templates for
modelling the vectors. There is a user-controlled freedom on how the vector is internally
represented. This is expressed by a choice of coordinate system which is supplied as a
template parameter when the vector is constructed. Furthermore each coordinate system is
itself a template, so that the user can specify the underlying scalar type.
In more detail, the main characteristics of GenVector are:

*   **Optimal runtime performances**

    We try to minimize any overhead in the run-time performances. We have deliberately
    avoided to have any virtual function and even virtual destructors in the classes and
    have inlined as much as possible all the functions. For this reason, we have chosen to
    use template classes to implement the GenVector concepts instead of abstract or base
    classes and virtual functions.

*   **Points and Vector concept**

    Mathematically vectors and points are two distinct concepts. They have different
    transformations, like vectors only rotate while points rotate and translate. You can
    add two vectors but not two points and the difference between two points is a vector.
    We then distinguish for the 2 and 3 dimensional case, between points and vectors,
    modeling them with different classes:
    *   ROOT::Math::DisplacementVector3D and ROOT::Math::DisplacementVector2D template
        classes describing 3 and 2 component direction and magnitude vectors, not rooted
        at any particular point;
    *   ROOT::Math::PositionVector3D template and ROOT::Math::PositionVector3D class
        modeling the points in 3 and 2 dimensions. For the 4D space-time vectors, we use
        the same class to model them, ROOT::Math::LorentzVector, since we have recognized
        a limited needs for modeling the functionality of a 4D point.

*   **Generic Coordinate System**

    The vector classes are based on a generic type of coordinate system, expressed as a
    template parameter of the class. Various classes exist to describe the various
    coordinates systems:
    *   **2D coordinate system** classes:
        *   ROOT::Math::Cartesian2D, based on <em>(x,y)</em> ;
        *   ROOT::Math::Polar2D, based on <em>(r, phi)</em> ;
    *   **3D coordinate system** classes:
        *   ROOT::Math::Cartesian3D, based on <em>(x,y,z)</em>;
        *   ROOT::Math::Polar3D, based on <em>(r, theta, phi)</em>;
        *   ROOT::Math::Cylindrical3D, based on <em>(rho, z, phi)</em>
        *   ROOT::Math::CylindricalEta3D, based on <em>(rho, eta, phi)</em>, where eta is the pseudo-rapidity;
    *   **4D coordinate system** classes:
        *   ROOT::Math::PxPyPzE4D, based on based on <em>(px,py,pz,E)</em>;
        *   ROOT::Math::PxPyPzM4D, based on based on <em>(px,py,pz,M)</em>;
        *   ROOT::Math::PtEtaPhiE4D, based on based on <em>(pt,eta,phi,E)</em>;
        *   ROOT::Math::PtEtaPhiM4D, based on based on <em>(pt,eta,phi,M)</em>;

    The angle _eta_ is defined between [0,\f$\pi\f$] and _phi_ between [-\f$\pi\f$,\f$\pi\f$].
    The angles are expressed in radians.

    Users can define the Vectors according to the coordinate type which
    is most efficient for their use. Transformations between the various coordinate
    systems are available through copy constructors or the assignment <em>(=)</em> operator.
    For maximum flexibility and minimize in some use case memory allocation, the
    coordinate system classes are templated on the scalar type. To avoid exposing
    templated parameter to the users, typedefs are defined for all types of vectors
    based an double's. See the \ref Vector3DPage "3D Vector", \ref Point3DPage
    "3D Point", \ref Vector2DPage "2D Vector and Point", and \ref LorentzVectorPage
    "LorentzVector" classes for all the possible types of vector classes which can
    be constructed by the user with the available coordinate system types.

*   **Coordinate System Tag**

    The 2D and 3D points and vector classes can be associated to a tag defining the
    coordinate system. This can be used to distinguish between vectors of different
    coordinate systems like global or local vectors. The coordinate system tag is a
    template parameter of the ROOT::Math::DisplacementVector3D
    (and ROOT::Math::DisplacementVector2D) and ROOT::Math::PositionVector3D
    (and ROOT::Math::PositionVector2D) classes. A default tag,
    ROOT::Math::DefaultCoordinateSystemTag, exists for users who don't need this
    functionality.

*   **Transformations**

    The transformations are modeled using simple (non-template) classes, using double as
    the scalar type to avoid too large numerical errors. The transformations are grouped
    in Rotations (in 3 dimensions), Lorentz transformations and Poincaré transformations,
    which are Translation/Rotation combinations. Each group has several members which may
    model physically equivalent trasformations but with different internal representations.
    Transformation classes can operate on all type of vectors using the  <em>operator()</em>
    or the _operator *_ and the transformations can also be combined via the _operator *_.
    In more detail the transformations available are:

    *   **3D Rotations classes**

        *   Rotation described by a 3x3 matrix (ROOT::Math::Rotation3D)
        *   Rotation described by Euler angles, Goldstein representation, (ROOT::Math::EulerAngles)
        *   Rotation described by 3-2-1 Euler angles (ROOT::Math::RotationZYX)
        *   Rotation described by a direction axis and an angle (ROOT::Math::AxisAngle)
        *   Rotation described by a quaternion (ROOT::Math::Quaternion)
        *   Optimized rotation around the x (ROOT::Math::RotationX), y (ROOT::Math::RotationY) and z
            (ROOT::Math::RotationZ) axis and described by just one angle.
    *   **3D Transformation**

        We describe the transformations defined as a composition between a rotation and a
        translation using the class ROOT::Math::Transform3D. It is important to note that
        transformations act differently on Vectors and Points. The Vectors only rotate,
        therefore when applying a transformation (rotation + translation) on a Vector,
        only the rotation operates while the translation has no effect. The interface for
        Transformations is similar to the one used in the CLHEP Geometry package
        (class [Transform3D](http://www.hep.phy.cam.ac.uk/lhcb/doc/CLHEP/1.9.1.2/html/classHepGeom_1_1Transform3D.html)).
        A class, ROOT::Math::Translation3D. describe transformations consisting of only a
        translation. Translation can be applied only on Points, applying them on Vector
        objects has no effect. The Translation3D class can be combined with both
        ROOT::Math::Rotation3D and ROOT::Math::Transform3D using the _operator *_ to
        obtain a new transformation as an instance of a Transform3D class.

    *   **Lorentz Rotation**

        *   Generic Lorentz Rotation described by a 4x4 matrix containing a 3D rotation part
            and a boost part (class ROOT::Math::LorentzRotation)
        *   A pure boost in an arbitrary direction and described by a 4x4 symmetric matrix or
            10 numbers (class ROOT::Math::Boost)
        *   Boost along the x (ROOT::Math::BoostX), y (ROOT::Math::BoostY) and z (ROOT::Math::BoostZ) axis.

Other main characteristics of the GenVector classes are:

*   **Minimal Vector classes interface**

    We have tried to keep the interface to a minimal level:
    *   We try to avoid methods providing the same functionality but with different names
        ( like getX() and x() ).
    *   we minimize the number of setter methods, avoiding methods which can be ambiguous
        and set the Vector classes in an inconsistent state. We provide only methods which
        set all the coordinates at the same time or set only the coordinates on which the
        vector is based, for example SetX() for a cartesian vector. We then enforce the use
        of transformations as rotations or translations (additions) for modifying the vector
        contents.
    *   The majority of the functionality, which is present in the CLHEP package, involving
        operations on two vectors, is moved in separated helper functions (see
        ROOT::Math::VectorUtil). This has the advantage that the basic interface will remain
        more stable with time while additional functions can be added easily.

*   **Naming Convention**

    As part of %ROOT, the GenVector package adheres to the prescribed ROOT naming convention,
    with some (approved) exceptions, as described here:
    *   Every class and function is in the _ROOT::Math_ namespace
    *   Member function names starts with upper-case letter, apart some exceptions (see later
        CLHEP compatibility)

*   **Compatibility with CLHEP Vector classes**

    *   For backward compatibility with CLHEP the Vector classes can be constructed easily
        from a CLHEP HepVector or HepLorentzVector, by using a template constructor, which
        requires only that the classes implement the accessors x(), y() and z() (and t()
        for the 4D).
    *   we have decided to provide Vector member function with the same naming convention
        as CLHEP for the most used functions like  <em>x()</em>,  <em>y()</em> and <em>z()</em>.

*   **Connection to Linear Algebra package**

    In some use cases, like in track reconstruction, it is needed to use the content of the
    vector and rotation classes in conjunction with linear algebra operations. We prefer to
    avoid any direct dependency to any Linear algebra package. However, we provide some hooks
    to convert to and from Linear Algebra classes.
    *   The vector and the transformation classes have methods which allow to get and set
        their data members (like SetCoordinates and GetCoordinates ) passing either a generic
        iterator or a pointer to a contiguous set of data, like a C array. This allows a
        easy connection with linear algebra package which allows creation of matrices using
        C arrays (like the %ROOT TMatrix classes) or iterators ( SMatrix classes )
    *   Multiplication between Linear Algebra matrix and GenVector Vectors is possible by
        using the template free functions ROOT::Math::VectorUtil::Mult. This works for any
        Linear Algebra matrix which implement the <em>operator(i,j)</em> and with first matrix element at _i=j=0_.

## Example of Usage

*   \ref Vector3DPage
*   \ref Point3DPage
*   \ref LorentzVectorPage
*   \ref TransformPage
*   ROOT::Math::VectorUtil (Helper functions)
*   \ref ExtUsagePage

## Packaging

This GenVector package is part of the \ref index and it can be built as an independent
package. A tar file can be downloaded from [here](../GenVector.tar.gz).

## Additional Documentation

A more detailed description of all the GenVector classes is available in this [document](http://seal.cern.ch/documents/mathlib/GenVector.pdf).

## References

1.  CLHEP Vector package ([User guide](http://proj-clhep.web.cern.ch/proj-clhep/manual/UserGuide/VectorDefs/index.html) and [reference doc](http://www.hep.phy.cam.ac.uk/lhcb/doc/CLHEP/1.9.1.2/html/dir_000027.html))
2.  [CLHEP Geometry package](http://www.hep.phy.cam.ac.uk/lhcb/doc/CLHEP/1.9.1.2/html/namespaceH)
3.  [%ROOT Physics Vector classes](http://root.cern.ch/root/html/PHYSICS_Index.html)
4.  [CMS Vector package](http://lcgapp.cern.ch/doxygen/SEAL/snapshot/html/dir_000007.html)
