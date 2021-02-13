#ifndef INCLUDE_TRACKMATHCORE_H
#define INCLUDE_TRACKMATHCORE_H

// dummy track class for testing I/o of matric

//#ifdef USE_REFLEX
//#define G__DICTIONARY
//#endif

#include "Math/Point3D.h"
#include "Math/Vector4D.h"
#include "Math/SMatrix.h"
#include "Rtypes.h" // for Double32_t
#include "TError.h"

#include <vector>
#include <string>
#include <iostream>




typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >        Vector4D;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<Double32_t> >    Vector4D32;

typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> >     Vector3D;

typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> >        Point3D;
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<Double32_t> >    Point3D32;

typedef ROOT::Math::SMatrix<double,4,4,ROOT::Math::MatRepStd<double,4,4> >      Matrix4D;
typedef ROOT::Math::SMatrix<Double32_t,4,4,ROOT::Math::MatRepStd<Double32_t,4,4> >      Matrix4D32;

typedef ROOT::Math::SMatrix<double,6,6,ROOT::Math::MatRepSym<double,6> >      SymMatrix6D;
typedef ROOT::Math::SMatrix<Double32_t,6,6,ROOT::Math::MatRepSym<Double32_t,6> >      SymMatrix6D32;;


// track class containing a vector and a point
class  TrackD {

public:

   typedef Matrix4D::const_iterator const_iterator;

   TrackD() {}

   TrackD(double * begin, double * end)  {
      double * itr = begin;
      fPos.SetCoordinates(itr, itr+3); itr +=3;
      fVec.SetCoordinates(itr,itr+4); itr+=4;
      R__ASSERT(itr == end);
   }

   enum {  kSize =        3 + 4    };

   static std::string Type() { return "TrackD"; }

   static bool IsD32() { return false; }

   TrackD & operator += (const TrackD & t) {
      fPos += Vector3D(t.fPos);
      fVec += t.fVec;
      return *this;
   }

   double Sum() const {
      double s = 0;
      double d[4];
      fPos.GetCoordinates(d,d+3);
      for (int i=0; i<3;++i)
         s+= d[i];
      fVec.GetCoordinates(d,d+4);
      for (int i=0; i<4;++i)
         s+= d[i];
      return s;
  }

   void Print() const {
     std::cout << "Point  " << fPos << std::endl;
     std::cout << "Vec    " << fVec << std::endl;
   }

private:

   Point3D    fPos;
   Vector4D   fVec;

};

// track class based on  of Double32


class  TrackD32 {

public:

   typedef Matrix4D::const_iterator const_iterator;

   TrackD32() {}

   enum {  kSize =        3 + 4    };

   static std::string Type() { return "TrackD32"; }

   static bool IsD32() { return true; }

   TrackD32(double * begin, double * end)  {
      double * itr = begin;
      fPos.SetCoordinates(itr, itr+3); itr +=3;
      fVec.SetCoordinates(itr,itr+4); itr+=4;
      R__ASSERT(itr == end);

   }


   TrackD32 & operator += (const TrackD32 & t) {
      fPos += Vector3D(t.fPos);
      fVec += t.fVec;
      return *this;
   }

   double Sum() const {
      double s = 0;
      double d[4];
      fPos.GetCoordinates(d,d+3);
      for (int i=0; i<3;++i)
         s+= d[i];
      fVec.GetCoordinates(d,d+4);
      for (int i=0; i<4;++i)
         s+= d[i];
      return s;
  }

   void Print() const {
     std::cout << "Point  " << fPos << std::endl;
     std::cout << "Vec    " << fVec << std::endl;
   }

private:

   Point3D32     fPos;
   Vector4D32    fVec;


};


// track class  (containing a vector and a point) and matrices

class  TrackErrD {

public:

   typedef Matrix4D::const_iterator const_iterator;

   TrackErrD() {}

   TrackErrD(double * begin, double * end)  {
      double * itr = begin;
      fPos.SetCoordinates(itr, itr+3); itr +=3;
      fVec.SetCoordinates(itr,itr+4); itr+=4;
      fMat = Matrix4D(itr,itr+16); itr += 16;
      fSymMat = SymMatrix6D(itr,itr+21);
      R__ASSERT(itr+21 == end);
   }

   enum {  kSize =
      3 + 4 + Matrix4D32::kSize + SymMatrix6D32::rep_type::kSize
   };

   static std::string Type() { return "TrackErrD"; }

   static bool IsD32() { return false; }

   TrackErrD & operator += (const TrackErrD & t) {
      fPos += Vector3D(t.fPos);
      fVec += t.fVec;
      fMat += t.fMat;
      fSymMat += t.fSymMat;
      return *this;
   }

   double Sum() const {
      double s = 0;
      double d[4];
      fPos.GetCoordinates(d,d+3);
      for (int i=0; i<3;++i)
         s+= d[i];
      fVec.GetCoordinates(d,d+4);
      for (int i=0; i<4;++i)
         s+= d[i];
      for (const_iterator itr = fMat.begin(); itr != fMat.end(); ++itr)
         s+= *itr;
      for (const_iterator itr = fSymMat.begin(); itr != fSymMat.end(); ++itr)
         s+= *itr;
      return s;
  }

   void Print() const {
     std::cout << "Point  " << fPos << std::endl;
     std::cout << "Vec    " << fVec << std::endl;
     std::cout << "Mat    " << fMat << std::endl;
     std::cout << "SymMat " << fSymMat << std::endl;
   }

private:

   Point3D    fPos;
   Vector4D   fVec;
   Matrix4D fMat;
   SymMatrix6D fSymMat;

};

// track class based on  of Double32


class  TrackErrD32 {

public:

   typedef Matrix4D::const_iterator const_iterator;

   TrackErrD32() {}

   enum {  kSize =
      3 + 4 + Matrix4D32::kSize + SymMatrix6D32::rep_type::kSize
   };

   static std::string Type() { return "TrackErrD32"; }

   static bool IsD32() { return true; }

   TrackErrD32(double * begin, double * end)  {
      double * itr = begin;
      fPos.SetCoordinates(itr, itr+3); itr +=3;
      fVec.SetCoordinates(itr,itr+4); itr+=4;
      fMat = Matrix4D32(itr,itr+16); itr += 16;
      fSymMat = SymMatrix6D32(itr,itr+21);
      R__ASSERT(itr+21 == end);

   }


   TrackErrD32 & operator += (const TrackErrD32 & t) {
      fPos += Vector3D(t.fPos);
      fVec += t.fVec;
      fMat += t.fMat;
      fSymMat += t.fSymMat;
      return *this;
   }

   double Sum() const {
      double s = 0;
      double d[4];
      fPos.GetCoordinates(d,d+3);
      for (int i=0; i<3;++i)
         s+= d[i];
      fVec.GetCoordinates(d,d+4);
      for (int i=0; i<4;++i)
         s+= d[i];
      for (const_iterator itr = fMat.begin(); itr != fMat.end(); ++itr)
         s+= *itr;
      for (const_iterator itr = fSymMat.begin(); itr != fSymMat.end(); ++itr)
         s+= *itr;
      return s;
  }

   void Print() const {
     std::cout << "Point  " << fPos << std::endl;
     std::cout << "Vec    " << fVec << std::endl;
     std::cout << "Mat    " << fMat << std::endl;
     std::cout << "SymMat " << fSymMat << std::endl;
   }

private:

   Point3D32     fPos;
   Vector4D32    fVec;
   Matrix4D32    fMat;
   SymMatrix6D32 fSymMat;

};

// class containning a vector of tracks

template<class T>
class VecTrack {

public:

   typedef typename std::vector<T>::const_iterator It;


   VecTrack() {
     // create klen empty trackD
       fTrks.reserve(kLen);
       for (int i = 0; i < kLen; ++i) {
         fTrks.push_back(T() );
       }
   }

   VecTrack(double * ibegin, double * iend)
     {
       fTrks.reserve(kLen);
       double * itr = ibegin;
       for (int i = 0; i < kLen; ++i) {
         fTrks.push_back(T(itr, itr + T::kSize) );
         itr += T::kSize;
       }
       R__ASSERT( itr == iend);
     }

   enum {  kLen = 3, kSize =  kLen*T::kSize };

   static std::string Type() { return "VecTrack<"+T::Type()+">"; }

   static bool IsD32() { return false; }


   VecTrack & operator += (const VecTrack & v) {
    for (unsigned int i = 0; i < fTrks.size() ; ++i)
       fTrks[i] += v.fTrks[i];


    return *this;
   }


  It begin() const { return fTrks.begin(); }
  It end() const  { return fTrks.end(); }


  double Sum() const {
    double s = 0;
    for (unsigned int i = 0; i < fTrks.size() ; ++i)
      s += fTrks[i].Sum();

    return s;
  }


   void Print() const {
     for (unsigned int i = 0; i < fTrks.size() ; ++i) {
       std::cout << "\n======> Track ========<" << i << std::endl;
       fTrks[i].Print();
     }
   }

private:

  std::vector<T>  fTrks;

};

// for instantiating the template VecTrackD class for reflex
struct Dummy {

  VecTrack<TrackD>     v1;
  VecTrack<TrackErrD> v2;
};

#endif
