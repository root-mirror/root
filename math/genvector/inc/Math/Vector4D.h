// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

#ifndef ROOT_Math_Vector4D
#define ROOT_Math_Vector4D

// Defines typedefs to specific vectors and forward declarations.
// Define additional (to Cartesian) coordinate system types.

namespace ROOT {

   namespace Math {

      // forward declarations of Lorentz Vectors and type defs definitions

      template<class CoordSystem> class LorentzVector;

      template<typename T> class PxPyPzE4D;
      template<typename T> class PtEtaPhiE4D;
      template<typename T> class PxPyPzM4D;
      template<typename T> class PtEtaPhiM4D;

      // For LorentzVector have only double classes (define the vector in the global ref frame)

      /// LorentzVector based on x,y,x,t (or px,py,pz,E) coordinates in double precision with metric (-,-,-,+)
      typedef LorentzVector<PxPyPzE4D<double> > XYZTVector;
      // For consistency
      typedef LorentzVector<PxPyPzE4D<double> > PxPyPzEVector;

      /// LorentzVector based on x,y,x,t (or px,py,pz,E) coordinates in float precision with metric (-,-,-,+)
      typedef LorentzVector< PxPyPzE4D <float> > XYZTVectorF;

      /// LorentzVector based on the x, y, z,  and Mass in double precision
      typedef LorentzVector<PxPyPzM4D<double> > PxPyPzMVector;

      /// LorentzVector based on the cylindrical coordinates Pt, eta, phi and E (rho, eta, phi, t) in double precision
      typedef LorentzVector<PtEtaPhiE4D<double> > PtEtaPhiEVector;

      /// LorentzVector based on the cylindrical coordinates pt, eta, phi and Mass in double precision
      typedef LorentzVector<PtEtaPhiM4D<double> > PtEtaPhiMVector;

   } // end namespace Math

} // end namespace ROOT

// Generic LorentzVector class definition.

#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/PtEtaPhiE4D.h"
#include "Math/GenVector/PxPyPzM4D.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

#include "Math/GenVector/LorentzVector.h"

#endif
