// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2005 , LCG ROOT MathLib Team  and                    *
  *                      FNAL LCG ROOT MathLib Team                    *
  *                                                                    *
  *                                                                    *
  **********************************************************************/

// Header file for class CylindricalEta3D
//
// Created by: Lorenzo Moneta  at Mon May 30 11:58:46 2005
// Major revamp:  M. Fischler  at Fri Jun 10 2005
//
// Last update: $Id$

//
#ifndef ROOT_Math_GenVector_CylindricalEta3D
#define ROOT_Math_GenVector_CylindricalEta3D  1

#include "Math/Math.h"

#include "Math/GenVector/etaMax.h"


#include <limits>
#include <cmath>


namespace ROOT {

namespace Math {

//__________________________________________________________________________________________
  /**
      Class describing a cylindrical coordinate system based on eta (pseudorapidity) instead of z.
      The base coordinates are rho (transverse component) , eta and phi
      Phi is restricted to be in the range [-PI,PI)

      @ingroup GenVector
  */

template <class T>
class CylindricalEta3D {

public :

  typedef T Scalar;

  /**
     Default constructor with rho=eta=phi=0
   */
  CylindricalEta3D() : fRho(0), fEta(0), fPhi(0) {  }

  /**
     Construct from rho eta and phi values
   */
  CylindricalEta3D(Scalar rho, Scalar eta, Scalar phi) :
             fRho(rho), fEta(eta), fPhi(phi) { Restrict(); }

  /**
     Construct from any Vector or coordinate system implementing
     Rho(), Eta() and Phi()
    */
  template <class CoordSystem >
  explicit CylindricalEta3D( const CoordSystem & v ) :
     fRho(v.Rho() ),  fEta(v.Eta() ),  fPhi(v.Phi() )
  {
     static Scalar bigEta = Scalar(-0.3) * log(std::numeric_limits<Scalar>::epsilon());
     if (std::fabs(fEta) > bigEta) {
        // This gives a small absolute adjustment in rho,
        // which, for large eta, results in a significant
        // improvement in the faithfullness of reproducing z.
        fRho *= v.Z() / Z();
    }
  }

   // for g++  3.2 and 3.4 on 32 bits found that the compiler generated copy ctor and assignment are much slower
   // re-implement them ( there is no no need to have them with g++4)

   /**
      copy constructor
   */
   CylindricalEta3D(const CylindricalEta3D & v) :
      fRho(v.Rho() ),  fEta(v.Eta() ),  fPhi(v.Phi() )  {   }

   /**
      assignment operator
   */
   CylindricalEta3D & operator= (const CylindricalEta3D & v) {
      fRho = v.Rho();
      fEta = v.Eta();
      fPhi = v.Phi();
      return *this;
   }

   /**
      Set internal data based on an array of 3 Scalar numbers
   */
   void SetCoordinates( const Scalar src[] )
   { fRho=src[0]; fEta=src[1]; fPhi=src[2]; Restrict(); }

   /**
      get internal data into an array of 3 Scalar numbers
   */
   void GetCoordinates( Scalar dest[] ) const
   { dest[0] = fRho; dest[1] = fEta; dest[2] = fPhi; }

   /**
      Set internal data based on 3 Scalar numbers
   */
   void SetCoordinates(Scalar  rho, Scalar  eta, Scalar  phi)
   { fRho=rho; fEta=eta; fPhi=phi; Restrict(); }

   /**
      get internal data into 3 Scalar numbers
   */
   void GetCoordinates(Scalar& rho, Scalar& eta, Scalar& phi) const
   {rho=fRho; eta=fEta; phi=fPhi;}

private:
   inline static Scalar pi() { return M_PI; }
   inline void Restrict() {
      if (fPhi <= -pi() || fPhi > pi()) fPhi = fPhi - floor(fPhi / (2 * pi()) + .5) * 2 * pi();
      return;
   }
public:

   // accessors

   T Rho()   const { return fRho; }
   T Eta()   const { return fEta; }
   T Phi()   const { return fPhi; }
   T X() const { return fRho * cos(fPhi); }
   T Y() const { return fRho * sin(fPhi); }
   T Z() const
   {
      return fRho > 0 ? fRho * sinh(fEta) : fEta == 0 ? 0 : fEta > 0 ? fEta - etaMax<T>() : fEta + etaMax<T>();
   }
   T R() const
   {
      return fRho > 0 ? fRho * cosh(fEta)
                      : fEta > etaMax<T>() ? fEta - etaMax<T>() : fEta < -etaMax<T>() ? -fEta - etaMax<T>() : 0;
   }
   T Mag2() const
   {
      const Scalar r = R();
      return r * r;
   }
   T Perp2() const { return fRho*fRho;            }
   T Theta() const { return fRho > 0 ? 2 * atan(exp(-fEta)) : (fEta >= 0 ? 0 : pi()); }

   // setters (only for data members)


   /**
       set the rho coordinate value keeping eta and phi constant
   */
   void SetRho(T rho) {
      fRho = rho;
   }

   /**
       set the eta coordinate value keeping rho and phi constant
   */
   void SetEta(T eta) {
      fEta = eta;
   }

   /**
       set the phi coordinate value keeping rho and eta constant
   */
   void SetPhi(T phi) {
      fPhi = phi;
      Restrict();
   }

   /**
       set all values using cartesian coordinates
   */
   void SetXYZ(Scalar x, Scalar y, Scalar z);


   /**
      scale by a scalar quantity a --
      for cylindrical eta coords, as long as a >= 0, only rho changes!
   */
   void Scale (T a) {
      if (a < 0) {
         Negate();
         a = -a;
      }
      // angles do not change when scaling by a positive quantity
      if (fRho > 0) {
         fRho *= a;
      } else if ( fEta >  etaMax<T>() ) {
         fEta =  ( fEta-etaMax<T>())*a + etaMax<T>();
      } else if ( fEta < -etaMax<T>() ) {
         fEta =  ( fEta+etaMax<T>())*a - etaMax<T>();
      } // when rho==0 and eta is not above etaMax, vector represents 0
      // and remains unchanged
   }

   /**
      negate the vector
   */
   void Negate ( ) {
      fPhi = ( fPhi > 0 ? fPhi - pi() : fPhi + pi()  );
      fEta = -fEta;
   }

   // assignment operators
   /**
      generic assignment operator from any coordinate system
   */
   template <class CoordSystem >
   CylindricalEta3D & operator= ( const CoordSystem & c ) {
      fRho = c.Rho();
      fEta = c.Eta();
      fPhi = c.Phi();
      return *this;
   }

   /**
      Exact component-by-component equality
      Note: Peculiar representaions of the zero vector such as (0,1,0) will
      not test as equal to one another.
   */
   bool operator==(const CylindricalEta3D & rhs) const {
      return fRho == rhs.fRho && fEta == rhs.fEta && fPhi == rhs.fPhi;
   }
   bool operator!= (const CylindricalEta3D & rhs) const
   {return !(operator==(rhs));}


   // ============= Compatibility section ==================

   // The following make this coordinate system look enough like a CLHEP
   // vector that an assignment member template can work with either
   T x() const { return X();}
   T y() const { return Y();}
   T z() const { return Z(); }

   // ============= Specializations for improved speed ==================

   // (none)

#if defined(__MAKECINT__) || defined(G__DICTIONARY)

   // ====== Set member functions for coordinates in other systems =======

   void SetX(Scalar x);

   void SetY(Scalar y);

   void SetZ(Scalar z);

   void SetR(Scalar r);

   void SetTheta(Scalar theta);


#endif


private:
   T fRho;
   T fEta;
   T fPhi;

};

  } // end namespace Math

} // end namespace ROOT


// move implementations here to avoid circle dependencies

#include "Math/GenVector/Cartesian3D.h"

#if defined(__MAKECINT__) || defined(G__DICTIONARY)
#include "Math/GenVector/GenVector_exception.h"
#include "Math/GenVector/Polar3D.h"
#endif

namespace ROOT {

  namespace Math {

template <class T>
void CylindricalEta3D<T>::SetXYZ(Scalar xx, Scalar yy, Scalar zz) {
   *this = Cartesian3D<Scalar>(xx, yy, zz);
}

#if defined(__MAKECINT__) || defined(G__DICTIONARY)


     // ====== Set member functions for coordinates in other systems =======


template <class T>
void CylindricalEta3D<T>::SetX(Scalar xx) {
   GenVector_exception e("CylindricalEta3D::SetX() is not supposed to be called");
   throw e;
   Cartesian3D<Scalar> v(*this); v.SetX(xx);
   *this = CylindricalEta3D<Scalar>(v);
}
template <class T>
void CylindricalEta3D<T>::SetY(Scalar yy) {
   GenVector_exception e("CylindricalEta3D::SetY() is not supposed to be called");
   throw e;
   Cartesian3D<Scalar> v(*this); v.SetY(yy);
   *this = CylindricalEta3D<Scalar>(v);
}
template <class T>
void CylindricalEta3D<T>::SetZ(Scalar zz) {
   GenVector_exception e("CylindricalEta3D::SetZ() is not supposed to be called");
   throw e;
   Cartesian3D<Scalar> v(*this); v.SetZ(zz);
   *this = CylindricalEta3D<Scalar>(v);
}
template <class T>
void CylindricalEta3D<T>::SetR(Scalar r) {
   GenVector_exception e("CylindricalEta3D::SetR() is not supposed to be called");
   throw e;
   Polar3D<Scalar> v(*this); v.SetR(r);
   *this = CylindricalEta3D<Scalar>(v);
}
template <class T>
void CylindricalEta3D<T>::SetTheta(Scalar theta) {
   GenVector_exception e("CylindricalEta3D::SetTheta() is not supposed to be called");
   throw e;
   Polar3D<Scalar> v(*this); v.SetTheta(theta);
   *this = CylindricalEta3D<Scalar>(v);
}

#endif


  } // end namespace Math

} // end namespace ROOT



#endif /* ROOT_Math_GenVector_CylindricalEta3D  */
