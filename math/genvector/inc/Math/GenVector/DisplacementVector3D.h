// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2005 , LCG ROOT MathLib Team and                     *
  *                      FNAL LCG ROOT MathLib Team                    *
  *                                                                    *
  *                                                                    *
  **********************************************************************/

// Header source file for class DisplacementVector3D
//
// Created by: Lorenzo Moneta  at Mon May 30 12:21:43 2005
// Major rewrite: M. FIschler  at Wed Jun  8  2005
//
// Last update: $Id$
//

#ifndef ROOT_Math_GenVector_DisplacementVector3D
#define ROOT_Math_GenVector_DisplacementVector3D  1

#include "Math/GenVector/Cartesian3D.h"

#include "Math/GenVector/PositionVector3Dfwd.h"

#include "Math/GenVector/GenVectorIO.h"

#include "Math/GenVector/BitReproducible.h"

#include "Math/GenVector/CoordinateSystemTags.h"

#include <cassert>


namespace ROOT {

  namespace Math {


//__________________________________________________________________________________________
    /**
     Class describing a generic displacement vector in 3 dimensions.
     This class is templated on the type of Coordinate system.
     One example is the XYZVector which is a vector based on
     double precision x,y,z data members by using the
     ROOT::Math::Cartesian3D<double> Coordinate system.
     The class is having also an extra template parameter, the coordinate system tag,
     to be able to identify (tag) vector described in different reference coordinate system,
     like global or local coordinate systems.

     @ingroup GenVector
    */

    template <class CoordSystem, class Tag = DefaultCoordinateSystemTag >
    class DisplacementVector3D {

    public:

      typedef typename CoordSystem::Scalar Scalar;
      typedef CoordSystem CoordinateType;
      typedef Tag  CoordinateSystemTag;

      // ------ ctors ------

      /**
          Default constructor. Construct an empty object with zero values
      */
      DisplacementVector3D ( ) :   fCoordinates()  { }


      /**
         Construct from three values of type <em>Scalar</em>.
         In the case of a XYZVector the values are x,y,z
         In the case of  a polar vector they are r,theta, phi
      */
      DisplacementVector3D(Scalar a, Scalar b, Scalar c) :
        fCoordinates ( a , b,  c )  { }

     /**
          Construct from a displacement vector expressed in different
          coordinates, or using a different Scalar type, but with same coordinate system tag
      */
      template <class OtherCoords>
      explicit DisplacementVector3D( const DisplacementVector3D<OtherCoords, Tag> & v) :
        fCoordinates ( v.Coordinates() ) { }


      /**
         Construct from a position vector expressed in different coordinates
         but with the same coordinate system tag
      */
      template <class OtherCoords>
      explicit DisplacementVector3D( const PositionVector3D<OtherCoords,Tag> & p) :
        fCoordinates ( p.Coordinates() ) { }


      /**
          Construct from a foreign 3D vector type, for example, Hep3Vector
          Precondition: v must implement methods x(), y() and z()
      */
      template <class ForeignVector>
      explicit DisplacementVector3D( const ForeignVector & v) :
        fCoordinates ( Cartesian3D<Scalar>( v.x(), v.y(), v.z() ) ) { }


#ifdef LATER
      /**
         construct from a generic linear algebra  vector of at least size 3
         implementing operator [].
         \par v  LAVector
         \par index0   index where coordinates starts (typically zero)
         It works for all Coordinates types,
         ( x= v[index0] for Cartesian and r=v[index0] for Polar )
      */
      template <class LAVector>
      DisplacementVector3D(const LAVector & v, size_t index0 ) {
        fCoordinates = CoordSystem ( v[index0], v[index0+1], v[index0+2] );
      }
#endif

      // compiler-generated copy ctor and dtor are fine.

      // ------ assignment ------

      /**
          Assignment operator from a displacement vector of arbitrary type
      */
      template <class OtherCoords>
      DisplacementVector3D & operator=
                        ( const DisplacementVector3D<OtherCoords, Tag> & v) {
        fCoordinates = v.Coordinates();
        return *this;
      }

      /**
         Assignment operator from a position vector
         (not necessarily efficient unless one or the other is Cartesian)
      */
      template <class OtherCoords>
      DisplacementVector3D & operator=
                        ( const PositionVector3D<OtherCoords,Tag> & rhs) {
         SetXYZ(rhs.x(), rhs.y(), rhs.z());
         return *this;
      }


      /**
          Assignment from a foreign 3D vector type, for example, Hep3Vector
          Precondition: v must implement methods x(), y() and z()
      */
      template <class ForeignVector>
      DisplacementVector3D & operator= ( const ForeignVector & v) {
        SetXYZ( v.x(),  v.y(), v.z() );
        return *this;
      }


#ifdef LATER
      /**
         assign from a generic linear algebra  vector of at least size 3
         implementing operator []. This could be also a C array
         \par v  LAVector
         \par index0   index where coordinates starts (typically zero)
         It works for all Coordinates types,
         ( x= v[index0] for Cartesian and r=v[index0] for Polar )
      */
      template <class LAVector>
      DisplacementVector3D & assignFrom(const LAVector & v, size_t index0 = 0) {
        fCoordinates = CoordSystem  ( v[index0], v[index0+1], v[index0+2] );
        return *this;
      }
#endif

      // ------ Set, Get, and access coordinate data ------

      /**
          Retrieve a copy of the coordinates object
      */
      CoordSystem Coordinates() const {
        return fCoordinates;
      }

      /**
         Set internal data based on a C-style array of 3 Scalar numbers
       */
      DisplacementVector3D<CoordSystem, Tag>& SetCoordinates( const Scalar src[] )
       { fCoordinates.SetCoordinates(src); return *this; }

      /**
         Set internal data based on 3 Scalar numbers
       */
      DisplacementVector3D<CoordSystem, Tag>& SetCoordinates( Scalar a, Scalar b, Scalar c )
       { fCoordinates.SetCoordinates(a, b, c); return *this; }

      /**
         Set internal data based on 3 Scalars at *begin to *end
       */
      template <class IT>
      DisplacementVector3D<CoordSystem, Tag>& SetCoordinates( IT begin, IT end  )
      { IT a = begin; IT b = ++begin; IT c = ++begin;
        (void)end;
        assert (++begin==end);
        SetCoordinates (*a,*b,*c);
        return *this;
      }

      /**
         get internal data into 3 Scalar numbers
       */
      void GetCoordinates( Scalar& a, Scalar& b, Scalar& c ) const
                            { fCoordinates.GetCoordinates(a, b, c);  }

      /**
         get internal data into a C-style array of 3 Scalar numbers
       */
      void GetCoordinates( Scalar dest[] ) const
                            { fCoordinates.GetCoordinates(dest);  }

      /**
         get internal data into 3 Scalars at *begin to *end (3 past begin)
       */
      template <class IT>
      void GetCoordinates( IT begin, IT end ) const
      { IT a = begin; IT b = ++begin; IT c = ++begin;
        (void)end;
        assert (++begin==end);
        GetCoordinates (*a,*b,*c);
      }
      /**
         get internal data into 3 Scalars starting at *begin
       */
      template <class IT>
      void GetCoordinates( IT begin) const {
         Scalar a = Scalar(0);
         Scalar b = Scalar(0);
         Scalar c = Scalar(0);
         GetCoordinates(a, b, c);
         *begin++ = a;
         *begin++ = b;
         *begin   = c;
      }

      /**
         set the values of the vector from the cartesian components (x,y,z)
         (if the vector is held in polar or cylindrical eta coordinates,
         then (x, y, z) are converted to that form)
       */
      DisplacementVector3D<CoordSystem, Tag>& SetXYZ (Scalar a, Scalar b, Scalar c) {
         fCoordinates.SetXYZ(a, b, c);
         return *this;
      }

      // ------------------- Equality -----------------

      /**
        Exact equality
       */
      bool operator==(const DisplacementVector3D & rhs) const {
        return fCoordinates==rhs.fCoordinates;
      }
      bool operator!= (const DisplacementVector3D & rhs) const {
        return !(operator==(rhs));
      }

      // ------ Individual element access, in various coordinate systems ------

      /**
          Cartesian X, converting if necessary from internal coordinate system.
      */
      Scalar X() const { return fCoordinates.X(); }

      /**
          Cartesian Y, converting if necessary from internal coordinate system.
      */
      Scalar Y() const { return fCoordinates.Y(); }

      /**
          Cartesian Z, converting if necessary from internal coordinate system.
      */
      Scalar Z() const { return fCoordinates.Z(); }

      /**
          Polar R, converting if necessary from internal coordinate system.
      */
      Scalar R() const { return fCoordinates.R(); }

      /**
          Polar theta, converting if necessary from internal coordinate system.
      */
      Scalar Theta() const { return fCoordinates.Theta(); }

      /**
          Polar phi, converting if necessary from internal coordinate system.
      */
      Scalar Phi() const { return fCoordinates.Phi(); }

      /**
          Polar eta, converting if necessary from internal coordinate system.
      */
      Scalar Eta() const { return fCoordinates.Eta(); }

      /**
          Cylindrical transverse component rho
      */
      Scalar Rho() const { return fCoordinates.Rho(); }

      // ----- Other fundamental properties -----

      /**
          Magnitute squared ( r^2 in spherical coordinate)
      */
      Scalar Mag2() const { return fCoordinates.Mag2();}

      /**
         Transverse component squared (rho^2 in cylindrical coordinates.
      */
      Scalar Perp2() const { return fCoordinates.Perp2();}

      /**
         return unit vector parallel to this (scalar)
      */
      template <typename SCALAR = Scalar, typename std::enable_if<std::is_arithmetic<SCALAR>::value>::type * = nullptr>
      DisplacementVector3D Unit() const
      {
         const auto tot = R();
         return tot == 0 ? *this : DisplacementVector3D(*this) / tot;
      }

      /**
         return unit vector parallel to this (vector)
      */
      template <typename SCALAR = Scalar, typename std::enable_if<!std::is_arithmetic<SCALAR>::value>::type * = nullptr>
      DisplacementVector3D Unit() const
      {
         SCALAR tot            = R();
         tot(tot == SCALAR(0)) = SCALAR(1);
         return DisplacementVector3D(*this) / tot;
      }

      // ------ Setting of individual elements present in coordinate system ------

      /**
         Change X - Cartesian3D coordinates only
      */
       DisplacementVector3D<CoordSystem, Tag>& SetX (Scalar xx) { fCoordinates.SetX(xx); return *this;}

      /**
         Change Y - Cartesian3D coordinates only
      */
       DisplacementVector3D<CoordSystem, Tag>& SetY (Scalar yy) { fCoordinates.SetY(yy); return *this;}

      /**
         Change Z - Cartesian3D coordinates only
      */
       DisplacementVector3D<CoordSystem, Tag>& SetZ (Scalar zz) { fCoordinates.SetZ(zz); return *this;}

      /**
         Change R - Polar3D coordinates only
      */
       DisplacementVector3D<CoordSystem, Tag>& SetR (Scalar rr) { fCoordinates.SetR(rr); return *this;}

      /**
         Change Theta - Polar3D coordinates only
      */
       DisplacementVector3D<CoordSystem, Tag>& SetTheta (Scalar ang) { fCoordinates.SetTheta(ang); return *this;}

      /**
         Change Phi - Polar3D or CylindricalEta3D coordinates
      */
       DisplacementVector3D<CoordSystem, Tag>& SetPhi (Scalar ang) { fCoordinates.SetPhi(ang); return *this;}

      /**
         Change Rho - CylindricalEta3D coordinates only
      */
       DisplacementVector3D<CoordSystem, Tag>& SetRho (Scalar rr) { fCoordinates.SetRho(rr); return *this;}

      /**
         Change Eta - CylindricalEta3D coordinates only
      */
       DisplacementVector3D<CoordSystem, Tag>& SetEta (Scalar etaval) { fCoordinates.SetEta(etaval); return *this;}


      // ------ Operations combining two vectors ------
      // -- need to have the specialized version in order to avoid

      /**
       Return the scalar (dot) product of two displacement vectors.
       It is possible to perform the product for any type of vector coordinates,
       but they must have the same coordinate system tag
      */
      template< class OtherCoords >
      Scalar Dot( const  DisplacementVector3D<OtherCoords,Tag>  & v) const {
        return X()*v.X() + Y()*v.Y() + Z()*v.Z();
      }
      /**
          Return the scalar (dot) product of two vectors.
          It is possible to perform the product for any classes
          implementing x(), y() and z() member functions
      */
      template< class OtherVector >
      Scalar Dot( const  OtherVector & v) const {
        return X()*v.x() + Y()*v.y() + Z()*v.z();
      }

      /**
       Return vector (cross) product of two displacement vectors,
       as a vector in the coordinate system of this class.
       It is possible to perform the product for any type of vector coordinates,
       but they must have the same coordinate system tag
      */
      template <class OtherCoords>
      DisplacementVector3D Cross( const DisplacementVector3D<OtherCoords,Tag>  & v) const {
        DisplacementVector3D  result;
        result.SetXYZ (  Y()*v.Z() - v.Y()*Z(),
                         Z()*v.X() - v.Z()*X(),
                         X()*v.Y() - v.X()*Y() );
        return result;
      }
      /**
         Return vector (cross) product of two  vectors,
         as a vector in the coordinate system of this class.
          It is possible to perform the product for any classes
          implementing X(), Y() and Z() member functions
      */
      template <class OtherVector>
      DisplacementVector3D Cross( const OtherVector & v) const {
        DisplacementVector3D  result;
        result.SetXYZ (  Y()*v.z() - v.y()*Z(),
                         Z()*v.x() - v.z()*X(),
                         X()*v.y() - v.x()*Y() );
        return result;
      }



      /**
          Self Addition with a displacement vector.
      */
      template <class OtherCoords>
      DisplacementVector3D & operator+=
                        (const  DisplacementVector3D<OtherCoords,Tag> & v) {
        SetXYZ(  X() + v.X(), Y() + v.Y(), Z() + v.Z() );
        return *this;
      }

      /**
          Self Difference with a displacement vector.
      */
      template <class OtherCoords>
      DisplacementVector3D & operator-=
                        (const  DisplacementVector3D<OtherCoords,Tag> & v) {
        SetXYZ(  x() - v.x(), y() - v.y(), z() - v.z() );
        return *this;
      }


      /**
         multiply this vector by a scalar quantity
      */
      DisplacementVector3D & operator*= (Scalar a) {
        fCoordinates.Scale(a);
        return *this;
      }

      /**
         divide this vector by a scalar quantity
      */
      DisplacementVector3D & operator/= (Scalar a) {
        fCoordinates.Scale(1/a);
        return *this;
      }

      // The following methods (v*a and v/a) could instead be free functions.
      // They were moved into the class to solve a problem on AIX.

      /**
        Multiply a vector by a real number
      */
      DisplacementVector3D operator * ( Scalar a ) const {
        DisplacementVector3D tmp(*this);
        tmp *= a;
        return tmp;
      }

      /**
         Negative of the vector
       */
      DisplacementVector3D operator - ( ) const {
        return operator*( Scalar(-1) );
      }

      /**
         Positive of the vector, return itself
       */
      DisplacementVector3D operator + ( ) const {return *this;}

      /**
         Division of a vector with a real number
       */
      DisplacementVector3D operator/ (Scalar a) const {
        DisplacementVector3D tmp(*this);
        tmp /= a;
        return tmp;
      }


      // Methods providing limited backward name compatibility with CLHEP

      Scalar x()     const { return fCoordinates.X();     }
      Scalar y()     const { return fCoordinates.Y();     }
      Scalar z()     const { return fCoordinates.Z();     }
      Scalar r()     const { return fCoordinates.R();     }
      Scalar theta() const { return fCoordinates.Theta(); }
      Scalar phi()   const { return fCoordinates.Phi();   }
      Scalar eta()   const { return fCoordinates.Eta();   }
      Scalar rho()   const { return fCoordinates.Rho();   }
      Scalar mag2()  const { return fCoordinates.Mag2();  }
      Scalar perp2() const { return fCoordinates.Perp2(); }
      DisplacementVector3D unit() const {return Unit();}


    private:

       CoordSystem fCoordinates;  // internal coordinate system

#ifdef NOT_SURE_THIS_SHOULD_BE_FORBIDDEN
      /**
         Cross product involving a position vector is inappropriate
      */
      template <class T2>
      DisplacementVector3D Cross( const PositionVector3D<T2> & ) const;
#endif

      // the following methods should not compile

      // this should not compile (if from a vector or points with different tag
      template <class OtherCoords, class OtherTag>
      explicit DisplacementVector3D( const DisplacementVector3D<OtherCoords, OtherTag> & ) {}

      template <class OtherCoords, class OtherTag>
      explicit DisplacementVector3D( const PositionVector3D<OtherCoords, OtherTag> & ) {}

      template <class OtherCoords, class OtherTag>
      DisplacementVector3D & operator=( const DisplacementVector3D<OtherCoords, OtherTag> & );


      template <class OtherCoords, class OtherTag>
      DisplacementVector3D & operator=( const PositionVector3D<OtherCoords, OtherTag> & );

      template <class OtherCoords, class OtherTag>
      DisplacementVector3D & operator+=(const  DisplacementVector3D<OtherCoords, OtherTag> & );

      template <class OtherCoords, class OtherTag>
      DisplacementVector3D & operator-=(const  DisplacementVector3D<OtherCoords, OtherTag> & );

      template<class OtherCoords, class OtherTag >
      Scalar Dot( const  DisplacementVector3D<OtherCoords, OtherTag> &  ) const;

      template<class OtherCoords, class OtherTag >
      DisplacementVector3D Cross( const  DisplacementVector3D<OtherCoords, OtherTag> &  ) const;


    };

// ---------- DisplacementVector3D class template ends here ------------
// ---------------------------------------------------------------------



   /**
        Addition of DisplacementVector3D vectors.
        The (coordinate system) type of the returned vector is defined to
        be identical to that of the first vector, which is passed by value
    */
    template <class CoordSystem1, class CoordSystem2, class U>
    inline
    DisplacementVector3D<CoordSystem1,U>
    operator+(       DisplacementVector3D<CoordSystem1,U> v1,
               const DisplacementVector3D<CoordSystem2,U>  & v2) {
      return v1 += v2;
    }

    /**
        Difference between two DisplacementVector3D vectors.
        The (coordinate system) type of the returned vector is defined to
        be identical to that of the first vector.
    */
    template <class CoordSystem1, class CoordSystem2, class U>
    inline
    DisplacementVector3D<CoordSystem1,U>
    operator-( DisplacementVector3D<CoordSystem1,U> v1,
               DisplacementVector3D<CoordSystem2,U> const & v2) {
      return v1 -= v2;
    }

    //#endif // not __CINT__

    /**
       Multiplication of a displacement vector by real number  a*v
    */
    template <class CoordSystem, class U>
    inline
    DisplacementVector3D<CoordSystem,U>
    operator * ( typename DisplacementVector3D<CoordSystem,U>::Scalar a,
                 DisplacementVector3D<CoordSystem,U> v) {
      return v *= a;
      // Note - passing v by value and using operator *= may save one
      // copy relative to passing v by const ref and creating a temporary.
    }


    // v1*v2 notation for Cross product of two vectors is omitted,
    // since it is always confusing as to whether dot product is meant.



    // ------------- I/O to/from streams -------------

    template <class char_t, class traits_t, class T, class U,
              typename std::enable_if<std::is_arithmetic<typename DisplacementVector3D<T, U>::Scalar>::value>::type * =
                 nullptr>
    std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &os,
                                                     DisplacementVector3D<T, U> const &v)
    {
       if (os) {

          typename T::Scalar a, b, c;
          v.GetCoordinates(a, b, c);

          if (detail::get_manip(os, detail::bitforbit)) {
             detail::set_manip(os, detail::bitforbit, '\00');
             typedef GenVector_detail::BitReproducible BR;
             BR::Output(os, a);
             BR::Output(os, b);
             BR::Output(os, c);
          } else {
             os << detail::get_manip(os, detail::open) << a << detail::get_manip(os, detail::sep) << b
                << detail::get_manip(os, detail::sep) << c << detail::get_manip(os, detail::close);
          }
      }
      return os;
    }  // op<< <>()

    template <class char_t, class traits_t, class T, class U,
              typename std::enable_if<!std::is_arithmetic<typename DisplacementVector3D<T, U>::Scalar>::value>::type * =
                 nullptr>
    std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &os,
                                                     DisplacementVector3D<T, U> const &v)
    {
       if (os) {
          os << "{ ";
          for (std::size_t i = 0; i < PositionVector3D<T, U>::Scalar::Size; ++i) {
             os << "(" << v.x()[i] << "," << v.y()[i] << "," << v.z()[i] << ") ";
          }
          os << "}";
       }
       return os;
    } // op<< <>()

    template< class char_t, class traits_t, class T, class U >
      inline
      std::basic_istream<char_t,traits_t> &
      operator >> ( std::basic_istream<char_t,traits_t> & is
                  , DisplacementVector3D<T,U> & v
                  )
    {
      if( !is )  return is;

      typename T::Scalar a, b, c;

      if( detail::get_manip( is, detail::bitforbit ) )  {
        detail::set_manip( is, detail::bitforbit, '\00' );
        typedef GenVector_detail::BitReproducible BR;
        BR::Input(is, a);
        BR::Input(is, b);
        BR::Input(is, c);
      }
      else  {
        detail::require_delim( is, detail::open  );  is >> a;
        detail::require_delim( is, detail::sep   );  is >> b;
        detail::require_delim( is, detail::sep   );  is >> c;
        detail::require_delim( is, detail::close );
      }

      if( is )
        v.SetCoordinates(a, b, c);
      return is;

    }  // op>> <>()



  }  // namespace Math

}  // namespace ROOT


#endif /* ROOT_Math_GenVector_DisplacementVector3D  */

