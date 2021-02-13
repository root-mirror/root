// @(#)root/mathcore:$Id: 9ef2a4a7bd1b62c1293920c2af2f64791c75bdd8 $
// Authors: W. Brown, M. Fischler, L. Moneta    2005


/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for Vector Utility functions
//
// Created by: moneta  at Tue May 31 21:10:29 2005
//
// Last update: Tue May 31 21:10:29 2005
//
#ifndef ROOT_Math_GenVector_VectorUtil
#define ROOT_Math_GenVector_VectorUtil  1

#include "Math/Math.h"


#include "Math/GenVector/Boost.h"

namespace ROOT {

   namespace Math {


      // utility functions for vector classes



      /**
       Global Helper functions for generic Vector classes. Any Vector classes implementing some defined member functions,
       like  Phi() or Eta() or mag() can use these functions.
       The functions returning a scalar value, returns always double precision number even if the vector are
       based on another precision type

       @ingroup GenVector
       */


      namespace VectorUtil {


         // methods for 3D vectors

         /**
          Find aximutal Angle difference between two generic vectors ( v2.Phi() - v1.Phi() )
          The only requirements on the Vector classes is that they implement the Phi() method
          \param v1  Vector of any type implementing the Phi() operator
          \param v2  Vector of any type implementing the Phi() operator
          \return  Phi difference
          \f[ \Delta \phi = \phi_2 - \phi_1 \f]
          */
         template <class Vector1, class Vector2>
         inline typename Vector1::Scalar DeltaPhi( const Vector1 & v1, const Vector2 & v2) {
            typename Vector1::Scalar dphi = v2.Phi() - v1.Phi();
            if ( dphi > M_PI ) {
               dphi -= 2.0*M_PI;
            } else if ( dphi <= -M_PI ) {
               dphi += 2.0*M_PI;
            }
            return dphi;
         }



         /**
          Find square of the difference in pseudorapidity (Eta) and Phi betwen two generic vectors
          The only requirements on the Vector classes is that they implement the Phi() and Eta() method
          \param v1  Vector 1
          \param v2  Vector 2
          \return   Angle between the two vectors
          \f[ \Delta R2 = ( \Delta \phi )^2 + ( \Delta \eta )^2  \f]
          */
         template <class Vector1, class Vector2>
         inline typename Vector1::Scalar DeltaR2( const Vector1 & v1, const Vector2 & v2) {
            typename Vector1::Scalar dphi = DeltaPhi(v1,v2);
            typename Vector1::Scalar deta = v2.Eta() - v1.Eta();
            return dphi*dphi + deta*deta;
         }

	 /**
	  Find square of the difference in true rapidity (y) and Phi betwen two generic vectors
	  The only requirements on the Vector classes is that they implement the Phi() and Rapidity() method
	  \param v1  Vector 1
	  \param v2  Vector 2
	  \return   Angle between the two vectors
	  \f[ \Delta R2 = ( \Delta \phi )^2 + ( \Delta \y )^2  \f]
	  */
	 template <class Vector1, class Vector2>
	 inline typename Vector1::Scalar DeltaR2RapidityPhi( const Vector1 & v1, const Vector2 & v2) {
	    typename Vector1::Scalar dphi = DeltaPhi(v1,v2);
	    typename Vector1::Scalar drap = v2.Rapidity() - v1.Rapidity();
	    return dphi*dphi + drap*drap;
	 }

         /**
          Find difference in pseudorapidity (Eta) and Phi betwen two generic vectors
          The only requirements on the Vector classes is that they implement the Phi() and Eta() method
          \param v1  Vector 1
          \param v2  Vector 2
          \return   Angle between the two vectors
          \f[ \Delta R = \sqrt{  ( \Delta \phi )^2 + ( \Delta \eta )^2 } \f]
          */
         template <class Vector1, class Vector2>
         inline typename Vector1::Scalar DeltaR( const Vector1 & v1, const Vector2 & v2) {
            return std::sqrt( DeltaR2(v1,v2) );
         }

	/**
          Find difference in Rapidity (y) and Phi betwen two generic vectors
          The only requirements on the Vector classes is that they implement the Phi() and Rapidity() method
          \param v1  Vector 1
          \param v2  Vector 2
          \return   Angle between the two vectors
          \f[ \Delta R = \sqrt{  ( \Delta \phi )^2 + ( \Delta y )^2 } \f],
          */
         template <class Vector1, class Vector2>
         inline typename Vector1::Scalar DeltaRapidityPhi( const Vector1 & v1, const Vector2 & v2) {
            return std::sqrt( DeltaR2RapidityPhi(v1,v2) );
         }

         /**
          Find CosTheta Angle between two generic 3D vectors
          pre-requisite: vectors implement the X(), Y() and Z()
          \param v1  Vector v1
          \param v2  Vector v2
          \return   cosine of Angle between the two vectors
          \f[ \cos \theta = \frac { \vec{v1} \cdot \vec{v2} }{ | \vec{v1} | | \vec{v2} | } \f]
          */
         // this cannot be made all generic since Mag2() for 2, 3 or 4 D is different
         // need to have a specialization for polar Coordinates ??
         template <class Vector1, class Vector2>
         double CosTheta( const Vector1 &  v1, const Vector2  & v2) {
            double arg;
            double v1_r2 = v1.X()*v1.X() + v1.Y()*v1.Y() + v1.Z()*v1.Z();
            double v2_r2 = v2.X()*v2.X() + v2.Y()*v2.Y() + v2.Z()*v2.Z();
            double ptot2 = v1_r2*v2_r2;
            if(ptot2 <= 0) {
               arg = 0.0;
            }else{
               double pdot = v1.X()*v2.X() + v1.Y()*v2.Y() + v1.Z()*v2.Z();
               arg = pdot/std::sqrt(ptot2);
               if(arg >  1.0) arg =  1.0;
               if(arg < -1.0) arg = -1.0;
            }
            return arg;
         }


         /**
          Find Angle between two vectors.
          Use the CosTheta() function
          \param v1  Vector v1
          \param v2  Vector v2
          \return   Angle between the two vectors
          \f[ \theta = \cos ^{-1} \frac { \vec{v1} \cdot \vec{v2} }{ | \vec{v1} | | \vec{v2} | } \f]
          */
         template <class Vector1, class Vector2>
         inline double Angle( const  Vector1 & v1, const Vector2 & v2) {
            return std::acos( CosTheta(v1, v2) );
         }

         /**
          Find the projection of v along the given direction u.
          \param v  Vector v for which the propjection is to be found
          \param u  Vector specifying the direction
          \return   Vector projection (same type of v)
          \f[ \vec{proj} = \frac{ \vec{v}  \cdot \vec{u} }{|\vec{u}|}\vec{u} \f]
          Precondition is that Vector1 implements Dot function and Vector2 implements X(),Y() and Z()
          */
         template <class Vector1, class Vector2>
         Vector1 ProjVector( const  Vector1 & v, const Vector2 & u) {
            double magU2 = u.X()*u.X() + u.Y()*u.Y() + u.Z()*u.Z();
            if (magU2 == 0) return Vector1(0,0,0);
            double d = v.Dot(u)/magU2;
            return Vector1( u.X() * d, u.Y() * d, u.Z() * d);
         }

         /**
          Find the vector component of v perpendicular to the given direction of u
          \param v  Vector v for which the perpendicular component is to be found
          \param u  Vector specifying the direction
          \return   Vector component of v which is perpendicular to u
          \f[ \vec{perp} = \vec{v} -  \frac{ \vec{v}  \cdot \vec{u} }{|\vec{u}|}\vec{u} \f]
          Precondition is that Vector1 implements Dot function and Vector2 implements X(),Y() and Z()
          */
         template <class Vector1, class Vector2>
         inline Vector1 PerpVector( const  Vector1 & v, const Vector2 & u) {
            return v - ProjVector(v,u);
         }

         /**
          Find the magnitude square of the vector component of v perpendicular to the given direction of u
          \param v  Vector v for which the perpendicular component is to be found
          \param u  Vector specifying the direction
          \return   square value of the component of v which is perpendicular to u
          \f[ perp = | \vec{v} -  \frac{ \vec{v}  \cdot \vec{u} }{|\vec{u}|}\vec{u} |^2 \f]
          Precondition is that Vector1 implements Dot function and Vector2 implements X(),Y() and Z()
          */
         template <class Vector1, class Vector2>
         inline double Perp2( const  Vector1 & v, const Vector2 & u) {
            double magU2 = u.X()*u.X() + u.Y()*u.Y() + u.Z()*u.Z();
            double prjvu = v.Dot(u);
            double magV2 = v.Dot(v);
            return magU2 > 0.0 ? magV2-prjvu*prjvu/magU2 : magV2;
         }

         /**
          Find the magnitude of the vector component of v perpendicular to the given direction of u
          \param v  Vector v for which the perpendicular component is to be found
          \param u  Vector specifying the direction
          \return   value of the component of v which is perpendicular to u
          \f[ perp = | \vec{v} -  \frac{ \vec{v}  \cdot \vec{u} }{|\vec{u}|}\vec{u} | \f]
          Precondition is that Vector1 implements Dot function and Vector2 implements X(),Y() and Z()
          */
         template <class Vector1, class Vector2>
         inline double Perp( const  Vector1 & v, const Vector2 & u) {
            return std::sqrt(Perp2(v,u) );
         }



         // Lorentz Vector functions


         /**
          return the invariant mass of two LorentzVector
          The only requirement on the LorentzVector is that they need to implement the
          X() , Y(), Z() and E() methods.
          \param v1 LorenzVector 1
          \param v2 LorenzVector 2
          \return invariant mass M
          \f[ M_{12} = \sqrt{ (\vec{v1} + \vec{v2} ) \cdot (\vec{v1} + \vec{v2} ) } \f]
          */
         template <class Vector1, class Vector2>
         inline typename Vector1::Scalar InvariantMass( const Vector1 & v1, const Vector2 & v2) {
            typedef typename  Vector1::Scalar Scalar;
            Scalar ee = (v1.E() + v2.E() );
            Scalar xx = (v1.X() + v2.X() );
            Scalar yy = (v1.Y() + v2.Y() );
            Scalar zz = (v1.Z() + v2.Z() );
            Scalar mm2 = ee*ee - xx*xx - yy*yy - zz*zz;
            return mm2 < 0.0 ? -std::sqrt(-mm2) : std::sqrt(mm2);
            //  PxPyPzE4D<double> q(xx,yy,zz,ee);
            //  return q.M();
            //return ( v1 + v2).mag();
         }

         template <class Vector1, class Vector2>
         inline typename Vector1::Scalar InvariantMass2( const Vector1 & v1, const Vector2 & v2) {
            typedef typename  Vector1::Scalar Scalar;
            Scalar ee = (v1.E() + v2.E() );
            Scalar xx = (v1.X() + v2.X() );
            Scalar yy = (v1.Y() + v2.Y() );
            Scalar zz = (v1.Z() + v2.Z() );
            Scalar mm2 = ee*ee - xx*xx - yy*yy - zz*zz;
            return mm2 ; // < 0.0 ? -std::sqrt(-mm2) : std::sqrt(mm2);
                         //  PxPyPzE4D<double> q(xx,yy,zz,ee);
                         //  return q.M();
                         //return ( v1 + v2).mag();
         }

         // rotation and transformations


#ifndef __CINT__
         /**
          rotation along X axis for a generic vector by an Angle alpha
          returning a new vector.
          The only pre requisite on the Vector is that it has to implement the X() , Y() and Z()
          and SetXYZ methods.
          */
         template <class Vector>
         Vector RotateX(const Vector & v, double alpha) {
            double sina = sin(alpha);
            double cosa = cos(alpha);
            double y2 = v.Y() * cosa - v.Z()*sina;
            double z2 = v.Z() * cosa + v.Y() * sina;
            Vector vrot;
            vrot.SetXYZ(v.X(), y2, z2);
            return vrot;
         }

         /**
          rotation along Y axis for a generic vector by an Angle alpha
          returning a new vector.
          The only pre requisite on the Vector is that it has to implement the X() , Y() and Z()
          and SetXYZ methods.
          */
         template <class Vector>
         Vector RotateY(const Vector & v, double alpha) {
            double sina = sin(alpha);
            double cosa = cos(alpha);
            double x2 = v.X() * cosa + v.Z() * sina;
            double z2 = v.Z() * cosa - v.X() * sina;
            Vector vrot;
            vrot.SetXYZ(x2, v.Y(), z2);
            return vrot;
         }

         /**
          rotation along Z axis for a generic vector by an Angle alpha
          returning a new vector.
          The only pre requisite on the Vector is that it has to implement the X() , Y() and Z()
          and SetXYZ methods.
          */
         template <class Vector>
         Vector RotateZ(const Vector & v, double alpha) {
            double sina = sin(alpha);
            double cosa = cos(alpha);
            double x2 = v.X() * cosa - v.Y() * sina;
            double y2 = v.Y() * cosa + v.X() * sina;
            Vector vrot;
            vrot.SetXYZ(x2, y2, v.Z());
            return vrot;
         }


         /**
          rotation on a generic vector using a generic rotation class.
          The only requirement on the vector is that implements the
          X(), Y(), Z() and SetXYZ methods.
          The requirement on the rotation matrix is that need to implement the
          (i,j) operator returning the matrix element with R(0,0) = xx element
          */
         template<class Vector, class RotationMatrix>
         Vector Rotate(const Vector &v, const RotationMatrix & rot) {
            double xX = v.X();
            double yY = v.Y();
            double zZ = v.Z();
            double x2 =  rot(0,0)*xX + rot(0,1)*yY + rot(0,2)*zZ;
            double y2 =  rot(1,0)*xX + rot(1,1)*yY + rot(1,2)*zZ;
            double z2 =  rot(2,0)*xX + rot(2,1)*yY + rot(2,2)*zZ;
            Vector vrot;
            vrot.SetXYZ(x2,y2,z2);
            return vrot;
         }

         /**
          Boost a generic Lorentz Vector class using a generic 3D Vector class describing the boost
          The only requirement on the vector is that implements the
          X(), Y(), Z(), T() and SetXYZT methods.
          The requirement on the boost vector is that needs to implement the
          X(), Y() , Z()  retorning the vector elements describing the boost
          The beta of the boost must be <= 1 or a nul Lorentz Vector will be returned
          */
         template <class LVector, class BoostVector>
         LVector boost(const LVector & v, const BoostVector & b) {
            double bx = b.X();
            double by = b.Y();
            double bz = b.Z();
            double b2 = bx*bx + by*by + bz*bz;
            if (b2 >= 1) {
               GenVector::Throw ( "Beta Vector supplied to set Boost represents speed >= c");
               return LVector();
            }
            double gamma = 1.0 / std::sqrt(1.0 - b2);
            double bp = bx*v.X() + by*v.Y() + bz*v.Z();
            double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;
            double x2 = v.X() + gamma2*bp*bx + gamma*bx*v.T();
            double y2 = v.Y() + gamma2*bp*by + gamma*by*v.T();
            double z2 = v.Z() + gamma2*bp*bz + gamma*bz*v.T();
            double t2 = gamma*(v.T() + bp);
            LVector lv;
            lv.SetXYZT(x2,y2,z2,t2);
            return lv;
         }


         /**
          Boost a generic Lorentz Vector class along the X direction with a factor beta
          The only requirement on the vector is that implements the
          X(), Y(), Z(), T()  and SetXYZT methods.
          The beta of the boost must be <= 1 or a nul Lorentz Vector will be returned
          */
         template <class LVector, class T>
         LVector boostX(const LVector & v, T beta) {
            if (beta >= 1) {
               GenVector::Throw ("Beta Vector supplied to set Boost represents speed >= c");
               return LVector();
            }
            T gamma = 1.0/ std::sqrt(1.0 - beta*beta);
            typename LVector::Scalar x2 = gamma * v.X() + gamma * beta * v.T();
            typename LVector::Scalar t2 = gamma * beta * v.X() + gamma * v.T();

            LVector lv;
            lv.SetXYZT(x2,v.Y(),v.Z(),t2);
            return lv;
         }

         /**
          Boost a generic Lorentz Vector class along the Y direction with a factor beta
          The only requirement on the vector is that implements the
          X(), Y(), Z(), T()  methods and be constructed from x,y,z,t values
          The beta of the boost must be <= 1 or a nul Lorentz Vector will be returned
          */
         template <class LVector>
         LVector boostY(const LVector & v, double beta) {
            if (beta >= 1) {
               GenVector::Throw ("Beta Vector supplied to set Boost represents speed >= c");
               return LVector();
            }
            double gamma = 1.0/ std::sqrt(1.0 - beta*beta);
            double y2 = gamma * v.Y() + gamma * beta * v.T();
            double t2 = gamma * beta * v.Y() + gamma * v.T();
            LVector lv;
            lv.SetXYZT(v.X(),y2,v.Z(),t2);
            return lv;
         }

         /**
          Boost a generic Lorentz Vector class along the Z direction with a factor beta
          The only requirement on the vector is that implements the
          X(), Y(), Z(), T()  methods and be constructed from x,y,z,t values
          The beta of the boost must be <= 1 or a nul Lorentz Vector will be returned
          */
         template <class LVector>
         LVector boostZ(const LVector & v, double beta) {
            if (beta >= 1) {
               GenVector::Throw ( "Beta Vector supplied to set Boost represents speed >= c");
               return LVector();
            }
            double gamma = 1.0/ std::sqrt(1.0 - beta*beta);
            double z2 = gamma * v.Z() + gamma * beta * v.T();
            double t2 = gamma * beta * v.Z() + gamma * v.T();
            LVector lv;
            lv.SetXYZT(v.X(),v.Y(),z2,t2);
            return lv;
         }

#endif




         // MATRIX VECTOR MULTIPLICATION
         // cannot define an operator * otherwise conflicts with rotations
         // operations like Rotation3D * vector use Mult

         /**
          Multiplications of a generic matrices with a  DisplacementVector3D of any coordinate system.
          Assume that the matrix implements the operator( i,j) and that it has at least         3 columns and 3 rows. There is no check on the matrix size !!
          */
         template<class Matrix, class CoordSystem, class U>
         inline
         DisplacementVector3D<CoordSystem,U> Mult (const Matrix & m, const DisplacementVector3D<CoordSystem,U> & v) {
            DisplacementVector3D<CoordSystem,U> vret;
            vret.SetXYZ( m(0,0) * v.x() + m(0,1) * v.y() + m(0,2) * v.z() ,
                        m(1,0) * v.x() + m(1,1) * v.y() + m(1,2) * v.z() ,
                        m(2,0) * v.x() + m(2,1) * v.y() + m(2,2) * v.z() );
            return vret;
         }


         /**
          Multiplications of a generic matrices with a generic PositionVector
          Assume that the matrix implements the operator( i,j) and that it has at least         3 columns and 3 rows. There is no check on the matrix size !!
          */
         template<class Matrix, class CoordSystem, class U>
         inline
         PositionVector3D<CoordSystem,U> Mult (const Matrix & m, const PositionVector3D<CoordSystem,U> & p) {
            DisplacementVector3D<CoordSystem,U> pret;
            pret.SetXYZ( m(0,0) * p.x() + m(0,1) * p.y() + m(0,2) * p.z() ,
                        m(1,0) * p.x() + m(1,1) * p.y() + m(1,2) * p.z() ,
                        m(2,0) * p.x() + m(2,1) * p.y() + m(2,2) * p.z() );
            return pret;
         }


         /**
          Multiplications of a generic matrices with a  LorentzVector described
          in any coordinate system.
          Assume that the matrix implements the operator( i,j) and that it has at least         4 columns and 4 rows. There is no check on the matrix size !!
          */
         // this will not be ambigous with operator*(Scalar, LorentzVector) since that one     // Scalar is passed by value
         template<class CoordSystem, class Matrix>
         inline
         LorentzVector<CoordSystem> Mult (const Matrix & m, const LorentzVector<CoordSystem> & v) {
            LorentzVector<CoordSystem> vret;
            vret.SetXYZT( m(0,0)*v.x() + m(0,1)*v.y() + m(0,2)*v.z() + m(0,3)* v.t() ,
                         m(1,0)*v.x() + m(1,1)*v.y() + m(1,2)*v.z() + m(1,3)* v.t() ,
                         m(2,0)*v.x() + m(2,1)*v.y() + m(2,2)*v.z() + m(2,3)* v.t() ,
                         m(3,0)*v.x() + m(3,1)*v.y() + m(3,2)*v.z() + m(3,3)* v.t() );
            return vret;
         }



         // non-template utility functions for all objects


         /**
          Return a phi angle in the interval (0,2*PI]
          */
         double Phi_0_2pi(double phi);
         /**
          Returns phi angle in the interval (-PI,PI]
          */
         double  Phi_mpi_pi(double phi);



      }  // end namespace Vector Util



   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_GenVector_VectorUtil  */
