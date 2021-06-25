// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2005 , LCG ROOT MathLib Team                         *
  *                                                                    *
  *                                                                    *
  **********************************************************************/

// Header file for class LorentzVectorfwd
//
// Created by: moneta  at Tue May 31 21:06:43 2005
//
// Last update: Tue May 31 21:06:43 2005
//
#ifndef ROOT_Math_Vector4Dfwd
#define ROOT_Math_Vector4Dfwd  1


namespace ROOT {

  namespace Math {

    // Forward declarations of Lorentz Vectors and type defs definitions

    template<class CoordSystem> class LorentzVector;

    template<typename T> class PxPyPzE4D;
    template<typename T> class PtEtaPhiE4D;
    template<typename T> class PxPyPzM4D;
    template<typename T> class PtEtaPhiM4D;

    typedef LorentzVector<PxPyPzE4D<double> > XYZTVector;
    typedef LorentzVector<PxPyPzE4D<double> > PxPyPzEVector;
    typedef LorentzVector< PxPyPzE4D <float> > XYZTVectorF;
    typedef LorentzVector<PxPyPzM4D<double> > PxPyPzMVector;
    typedef LorentzVector<PtEtaPhiE4D<double> > PtEtaPhiEVector;
    typedef LorentzVector<PtEtaPhiM4D<double> > PtEtaPhiMVector;

  } // end namespace Math

} // end namespace ROOT

#endif

