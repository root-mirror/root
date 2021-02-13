// @(#)root/physics:$Id$
// Author: Pasha Murat , Peter Malzacher  12/02/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TLorentzVector
#define ROOT_TLorentzVector


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TLorentzVector                                                       //
//                                                                      //
// Place holder for real lorentz vector class.                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TVector3.h"
#include "TRotation.h"


class TLorentzRotation;


class TLorentzVector : public TObject {

private:

   TVector3 fP;  // 3 vector component
   Double_t fE;  // time or energy of (x,y,z,t) or (px,py,pz,e)

public:

   typedef Double_t Scalar;   // to be able to use it with the ROOT::Math::VectorUtil functions

   enum { kX=0, kY=1, kZ=2, kT=3, kNUM_COORDINATES=4, kSIZE=kNUM_COORDINATES };
   // Safe indexing of the coordinates when using with matrices, arrays, etc.

   TLorentzVector();

   TLorentzVector(Double_t x, Double_t y, Double_t z, Double_t t);
   // Constructor giving the components x, y, z, t.

   TLorentzVector(const Double_t * carray);
   TLorentzVector(const Float_t * carray);
   // Constructor from an array, not checked!

   TLorentzVector(const TVector3 & vector3, Double_t t);
   // Constructor giving a 3-Vector and a time component.

   TLorentzVector(const TLorentzVector & lorentzvector);
   // Copy constructor.

   virtual ~TLorentzVector(){};
   // Destructor

   // inline operator TVector3 () const;
   // inline operator TVector3 & ();
   // Conversion (cast) to TVector3.

   inline Double_t X() const;
   inline Double_t Y() const;
   inline Double_t Z() const;
   inline Double_t T() const;
   // Get position and time.

   inline void SetX(Double_t a);
   inline void SetY(Double_t a);
   inline void SetZ(Double_t a);
   inline void SetT(Double_t a);
   // Set position and time.

   inline Double_t Px() const;
   inline Double_t Py() const;
   inline Double_t Pz() const;
   inline Double_t P()  const;
   inline Double_t E()  const;
   inline Double_t Energy() const;
   // Get momentum and energy.

   inline void SetPx(Double_t a);
   inline void SetPy(Double_t a);
   inline void SetPz(Double_t a);
   inline void SetE(Double_t a);
   // Set momentum and energy.

   inline TVector3 Vect() const ;
   // Get spatial component.

   inline void SetVect(const TVector3 & vect3);
   // Set spatial component.

   inline Double_t Theta() const;
   inline Double_t CosTheta() const;
   inline Double_t Phi() const; //returns phi from -pi to pi
   inline Double_t Rho() const;
   // Get spatial vector components in spherical coordinate system.

   inline void SetTheta(Double_t theta);
   inline void SetPhi(Double_t phi);
   inline void SetRho(Double_t rho);
   // Set spatial vector components in spherical coordinate system.

   inline void SetPxPyPzE(Double_t px, Double_t py, Double_t pz, Double_t e);
   inline void SetXYZT(Double_t  x, Double_t  y, Double_t  z, Double_t t);
   inline void SetXYZM(Double_t  x, Double_t  y, Double_t  z, Double_t m);
   inline void SetPtEtaPhiM(Double_t pt, Double_t eta, Double_t phi, Double_t m);
   inline void SetPtEtaPhiE(Double_t pt, Double_t eta, Double_t phi, Double_t e);
   // Setters to provide the functionality (but a more meanigful name) of
   // the previous version eg SetV4... PsetV4...

   inline void GetXYZT(Double_t *carray) const;
   inline void GetXYZT(Float_t *carray) const;
   // Getters into an arry
   // no checking!

   Double_t operator () (int i) const;
   inline Double_t operator [] (int i) const;
   // Get components by index.

   Double_t & operator () (int i);
   inline Double_t & operator [] (int i);
   // Set components by index.

   inline TLorentzVector & operator = (const TLorentzVector &);
   // Assignment.

   inline TLorentzVector   operator +  (const TLorentzVector &) const;
   inline TLorentzVector & operator += (const TLorentzVector &);
   // Additions.

   inline TLorentzVector   operator -  (const TLorentzVector &) const;
   inline TLorentzVector & operator -= (const TLorentzVector &);
   // Subtractions.

   inline TLorentzVector operator - () const;
   // Unary minus.

   inline TLorentzVector operator * (Double_t a) const;
   inline TLorentzVector & operator *= (Double_t a);
   // Scaling with real numbers.

   inline Bool_t operator == (const TLorentzVector &) const;
   inline Bool_t operator != (const TLorentzVector &) const;
   // Comparisons.

   inline Double_t Perp2() const;
   // Transverse component of the spatial vector squared.

   inline Double_t Pt() const;
   inline Double_t Perp() const;
   // Transverse component of the spatial vector (R in cylindrical system).

   inline void SetPerp(Double_t);
   // Set the transverse component of the spatial vector.

   inline Double_t Perp2(const TVector3 & v) const;
   // Transverse component of the spatial vector w.r.t. given axis squared.

   inline Double_t Pt(const TVector3 & v) const;
   inline Double_t Perp(const TVector3 & v) const;
   // Transverse component of the spatial vector w.r.t. given axis.

   inline Double_t Et2() const;
   // Transverse energy squared.

   inline Double_t Et() const;
   // Transverse energy.

   inline Double_t Et2(const TVector3 &) const;
   // Transverse energy w.r.t. given axis squared.

   inline Double_t Et(const TVector3 &) const;
   // Transverse energy w.r.t. given axis.

   inline Double_t DeltaPhi(const TLorentzVector &) const;
   inline Double_t DeltaR(const TLorentzVector &, Bool_t useRapidity=kFALSE) const;
   inline Double_t DrEtaPhi(const TLorentzVector &) const;
   inline Double_t DrRapidityPhi(const TLorentzVector &) const;
   inline TVector2 EtaPhiVector();

   inline Double_t Angle(const TVector3 & v) const;
   // Angle wrt. another vector.

   inline Double_t Mag2() const;
   inline Double_t M2() const;
   // Invariant mass squared.

   inline Double_t Mag() const;
   inline Double_t M() const;
   // Invariant mass. If mag2() is negative then -sqrt(-mag2()) is returned.

   inline Double_t Mt2() const;
   // Transverse mass squared.

   inline Double_t Mt() const;
   // Transverse mass.

   inline Double_t Beta() const;
   inline Double_t Gamma() const;

   inline Double_t Dot(const TLorentzVector &) const;
   inline Double_t operator * (const TLorentzVector &) const;
   // Scalar product.

   inline void SetVectMag(const TVector3 & spatial, Double_t magnitude);
   inline void SetVectM(const TVector3 & spatial, Double_t mass);
   // Copy spatial coordinates, and set energy = sqrt(mass^2 + spatial^2)

   inline Double_t Plus() const;
   inline Double_t Minus() const;
   // Returns t +/- z.
   // Related to the positive/negative light-cone component,
   // which some define this way and others define as (t +/- z)/sqrt(2)

   inline TVector3 BoostVector() const ;
   // Returns the spatial components divided by the time component.

   void Boost(Double_t, Double_t, Double_t);
   inline void Boost(const TVector3 &);
   // Lorentz boost.

   Double_t Rapidity() const;
   // Returns the rapidity, i.e. 0.5*ln((E+pz)/(E-pz))

   inline Double_t Eta() const;
   inline Double_t PseudoRapidity() const;
   // Returns the pseudo-rapidity, i.e. -ln(tan(theta/2))

   inline void RotateX(Double_t angle);
   // Rotate the spatial component around the x-axis.

   inline void RotateY(Double_t angle);
   // Rotate the spatial component around the y-axis.

   inline void RotateZ(Double_t angle);
   // Rotate the spatial component around the z-axis.

   inline void RotateUz(const TVector3 & newUzVector);
   // Rotates the reference frame from Uz to newUz (unit vector).

   inline void Rotate(Double_t, const TVector3 &);
   // Rotate the spatial component around specified axis.

   inline TLorentzVector & operator *= (const TRotation &);
   inline TLorentzVector & Transform(const TRotation &);
   // Transformation with HepRotation.

   TLorentzVector & operator *= (const TLorentzRotation &);
   TLorentzVector & Transform(const TLorentzRotation &);
   // Transformation with HepLorenzRotation.

   virtual void        Print(Option_t *option="") const;

   ClassDef(TLorentzVector,4) // A four vector with (-,-,-,+) metric
};


//inline TLorentzVector operator * (const TLorentzVector &, Double_t a);
// moved to TLorentzVector::operator * (Double_t a)
inline TLorentzVector operator * (Double_t a, const TLorentzVector &);
// Scaling LorentzVector with a real number


inline Double_t TLorentzVector::X() const { return fP.X(); }
inline Double_t TLorentzVector::Y() const { return fP.Y(); }
inline Double_t TLorentzVector::Z() const { return fP.Z(); }
inline Double_t TLorentzVector::T() const { return fE; }

inline void TLorentzVector::SetX(Double_t a) { fP.SetX(a); }
inline void TLorentzVector::SetY(Double_t a) { fP.SetY(a); }
inline void TLorentzVector::SetZ(Double_t a) { fP.SetZ(a); }
inline void TLorentzVector::SetT(Double_t a) { fE = a; }

inline Double_t TLorentzVector::Px() const { return X(); }
inline Double_t TLorentzVector::Py() const { return Y(); }
inline Double_t TLorentzVector::Pz() const { return Z(); }
inline Double_t TLorentzVector::P()  const { return fP.Mag(); }
inline Double_t TLorentzVector::E()  const { return T(); }
inline Double_t TLorentzVector::Energy()  const { return T(); }

inline void TLorentzVector::SetPx(Double_t a) { SetX(a); }
inline void TLorentzVector::SetPy(Double_t a) { SetY(a); }
inline void TLorentzVector::SetPz(Double_t a) { SetZ(a); }
inline void TLorentzVector::SetE(Double_t a)  { SetT(a); }

inline TVector3 TLorentzVector::Vect() const { return fP; }

inline void TLorentzVector::SetVect(const TVector3 &p) { fP = p; }

inline Double_t TLorentzVector::Phi() const {
   return fP.Phi();
}

inline Double_t TLorentzVector::Theta() const {
   return fP.Theta();
}

inline Double_t TLorentzVector::CosTheta() const {
   return fP.CosTheta();
}


inline Double_t TLorentzVector::Rho() const {
   return fP.Mag();
}

inline void TLorentzVector::SetTheta(Double_t th) {
   fP.SetTheta(th);
}

inline void TLorentzVector::SetPhi(Double_t phi) {
   fP.SetPhi(phi);
}

inline void TLorentzVector::SetRho(Double_t rho) {
   fP.SetMag(rho);
}

inline void TLorentzVector::SetXYZT(Double_t  x, Double_t  y, Double_t  z, Double_t t) {
   fP.SetXYZ(x, y, z);
   SetT(t);
}

inline void TLorentzVector::SetPxPyPzE(Double_t px, Double_t py, Double_t pz, Double_t e) {
   SetXYZT(px, py, pz, e);
}

inline void TLorentzVector::SetXYZM(Double_t  x, Double_t  y, Double_t  z, Double_t m) {
   if ( m  >= 0 )
      SetXYZT( x, y, z, TMath::Sqrt(x*x+y*y+z*z+m*m) );
   else
      SetXYZT( x, y, z, TMath::Sqrt( TMath::Max((x*x+y*y+z*z-m*m), 0. ) ) );
}

inline void TLorentzVector::SetPtEtaPhiM(Double_t pt, Double_t eta, Double_t phi, Double_t m) {
   pt = TMath::Abs(pt);
   SetXYZM(pt*TMath::Cos(phi), pt*TMath::Sin(phi), pt*sinh(eta) ,m);
}

inline void TLorentzVector::SetPtEtaPhiE(Double_t pt, Double_t eta, Double_t phi, Double_t e) {
   pt = TMath::Abs(pt);
   SetXYZT(pt*TMath::Cos(phi), pt*TMath::Sin(phi), pt*sinh(eta) ,e);
}

inline void TLorentzVector::GetXYZT(Double_t *carray) const {
   fP.GetXYZ(carray);
   carray[3] = fE;
}

inline void TLorentzVector::GetXYZT(Float_t *carray) const{
   fP.GetXYZ(carray);
   carray[3] = fE;
}

inline Double_t & TLorentzVector::operator [] (int i)       { return (*this)(i); }
inline Double_t   TLorentzVector::operator [] (int i) const { return (*this)(i); }

inline TLorentzVector &TLorentzVector::operator = (const TLorentzVector & q) {
   fP = q.Vect();
   fE = q.T();
   return *this;
}

inline TLorentzVector TLorentzVector::operator + (const TLorentzVector & q) const {
   return TLorentzVector(fP+q.Vect(), fE+q.T());
}

inline TLorentzVector &TLorentzVector::operator += (const TLorentzVector & q) {
   fP += q.Vect();
   fE += q.T();
   return *this;
}

inline TLorentzVector TLorentzVector::operator - (const TLorentzVector & q) const {
   return TLorentzVector(fP-q.Vect(), fE-q.T());
}

inline TLorentzVector &TLorentzVector::operator -= (const TLorentzVector & q) {
   fP -= q.Vect();
   fE -= q.T();
   return *this;
}

inline TLorentzVector TLorentzVector::operator - () const {
   return TLorentzVector(-X(), -Y(), -Z(), -T());
}

inline TLorentzVector& TLorentzVector::operator *= (Double_t a) {
   fP *= a;
   fE *= a;
   return *this;
}

inline TLorentzVector TLorentzVector::operator * (Double_t a) const {
   return TLorentzVector(a*X(), a*Y(), a*Z(), a*T());
}

inline Bool_t TLorentzVector::operator == (const TLorentzVector & q) const {
   return (Vect() == q.Vect() && T() == q.T());
}

inline Bool_t TLorentzVector::operator != (const TLorentzVector & q) const {
   return (Vect() != q.Vect() || T() != q.T());
}

inline Double_t TLorentzVector::Perp2() const  { return fP.Perp2(); }

inline Double_t TLorentzVector::Perp()  const  { return fP.Perp(); }

inline Double_t TLorentzVector::Pt() const { return Perp(); }

inline void TLorentzVector::SetPerp(Double_t r) {
   fP.SetPerp(r);
}

inline Double_t TLorentzVector::Perp2(const TVector3 &v) const {
   return fP.Perp2(v);
}

inline Double_t TLorentzVector::Perp(const TVector3 &v) const {
   return fP.Perp(v);
}

inline Double_t TLorentzVector::Pt(const TVector3 &v) const {
   return Perp(v);
}

inline Double_t TLorentzVector::Et2() const {
   Double_t pt2 = fP.Perp2();
   return pt2 == 0 ? 0 : E()*E() * pt2/(pt2+Z()*Z());
}

inline Double_t TLorentzVector::Et() const {
   Double_t etet = Et2();
   return E() < 0.0 ? -sqrt(etet) : sqrt(etet);
}

inline Double_t TLorentzVector::Et2(const TVector3 & v) const {
   Double_t pt2 = fP.Perp2(v);
   Double_t pv = fP.Dot(v.Unit());
   return pt2 == 0 ? 0 : E()*E() * pt2/(pt2+pv*pv);
}

inline Double_t TLorentzVector::Et(const TVector3 & v) const {
   Double_t etet = Et2(v);
   return E() < 0.0 ? -sqrt(etet) : sqrt(etet);
}

inline Double_t TLorentzVector::DeltaPhi(const TLorentzVector & v) const {
   return TVector2::Phi_mpi_pi(Phi()-v.Phi());
}

inline Double_t TLorentzVector::Eta() const {
   return PseudoRapidity();
}

inline Double_t TLorentzVector::DeltaR(const TLorentzVector & v, const Bool_t useRapidity) const {
  if(useRapidity){
     Double_t drap = Rapidity()-v.Rapidity();
     Double_t dphi = TVector2::Phi_mpi_pi(Phi()-v.Phi());
     return TMath::Sqrt( drap*drap+dphi*dphi );
  } else {
    Double_t deta = Eta()-v.Eta();
    Double_t dphi = TVector2::Phi_mpi_pi(Phi()-v.Phi());
    return TMath::Sqrt( deta*deta+dphi*dphi );
  }
}

inline Double_t TLorentzVector::DrEtaPhi(const TLorentzVector & v) const{
   return DeltaR(v);
}

inline Double_t TLorentzVector::DrRapidityPhi(const TLorentzVector & v) const{
   return DeltaR(v, kTRUE);
}

inline TVector2 TLorentzVector::EtaPhiVector() {
   return TVector2 (Eta(),Phi());
}


inline Double_t TLorentzVector::Angle(const TVector3 &v) const {
   return fP.Angle(v);
}

inline Double_t TLorentzVector::Mag2() const {
   return T()*T() - fP.Mag2();
}

inline Double_t TLorentzVector::Mag() const {
   Double_t mm = Mag2();
   return mm < 0.0 ? -TMath::Sqrt(-mm) : TMath::Sqrt(mm);
}

inline Double_t TLorentzVector::M2() const { return Mag2(); }
inline Double_t TLorentzVector::M() const { return Mag(); }

inline Double_t TLorentzVector::Mt2() const {
   return E()*E() - Z()*Z();
}

inline Double_t TLorentzVector::Mt() const {
   Double_t mm = Mt2();
   return mm < 0.0 ? -TMath::Sqrt(-mm) : TMath::Sqrt(mm);
}

inline Double_t TLorentzVector::Beta() const {
   return fP.Mag() / fE;
}

inline Double_t TLorentzVector::Gamma() const {
   Double_t b = Beta();
   return 1.0/TMath::Sqrt(1- b*b);
}

inline void TLorentzVector::SetVectMag(const TVector3 & spatial, Double_t magnitude) {
   SetXYZM(spatial.X(), spatial.Y(), spatial.Z(), magnitude);
}

inline void TLorentzVector::SetVectM(const TVector3 & spatial, Double_t mass) {
   SetVectMag(spatial, mass);
}

inline Double_t TLorentzVector::Dot(const TLorentzVector & q) const {
   return T()*q.T() - Z()*q.Z() - Y()*q.Y() - X()*q.X();
}

inline Double_t TLorentzVector::operator * (const TLorentzVector & q) const {
   return Dot(q);
}

//Member functions Plus() and Minus() return the positive and negative
//light-cone components:
//
//  Double_t pcone = v.Plus();
//  Double_t mcone = v.Minus();
//
//CAVEAT: The values returned are T{+,-}Z. It is known that some authors
//find it easier to define these components as (T{+,-}Z)/sqrt(2). Thus
//check what definition is used in the physics you're working in and adapt
//your code accordingly.

inline Double_t TLorentzVector::Plus() const {
   return T() + Z();
}

inline Double_t TLorentzVector::Minus() const {
   return T() - Z();
}

inline TVector3 TLorentzVector::BoostVector() const {
   return TVector3(X()/T(), Y()/T(), Z()/T());
}

inline void TLorentzVector::Boost(const TVector3 & b) {
   Boost(b.X(), b.Y(), b.Z());
}

inline Double_t TLorentzVector::PseudoRapidity() const {
   return fP.PseudoRapidity();
}

inline void TLorentzVector::RotateX(Double_t angle) {
   fP.RotateX(angle);
}

inline void TLorentzVector::RotateY(Double_t angle) {
   fP.RotateY(angle);
}

inline void TLorentzVector::RotateZ(Double_t angle) {
   fP.RotateZ(angle);
}

inline void TLorentzVector::RotateUz(const TVector3 &newUzVector) {
   fP.RotateUz(newUzVector);
}

inline void TLorentzVector::Rotate(Double_t a, const TVector3 &v) {
   fP.Rotate(a,v);
}

inline TLorentzVector &TLorentzVector::operator *= (const TRotation & m) {
   fP *= m;
   return *this;
}

inline TLorentzVector &TLorentzVector::Transform(const TRotation & m) {
   fP.Transform(m);
   return *this;
}

inline TLorentzVector operator * (Double_t a, const TLorentzVector & p) {
   return TLorentzVector(a*p.X(), a*p.Y(), a*p.Z(), a*p.T());
}

inline TLorentzVector::TLorentzVector()
               : fP(), fE(0.0) {}

inline TLorentzVector::TLorentzVector(Double_t x, Double_t y, Double_t z, Double_t t)
               : fP(x,y,z), fE(t) {}

inline TLorentzVector::TLorentzVector(const Double_t * x0)
               : fP(x0), fE(x0[3]) {}

inline TLorentzVector::TLorentzVector(const Float_t * x0)
               : fP(x0), fE(x0[3]) {}

inline TLorentzVector::TLorentzVector(const TVector3 & p, Double_t e)
               : fP(p), fE(e) {}

inline TLorentzVector::TLorentzVector(const TLorentzVector & p) : TObject(p)
               , fP(p.Vect()), fE(p.T()) {}



inline Double_t TLorentzVector::operator () (int i) const
{
   //dereferencing operator const
   switch(i) {
      case kX:
	 return fP.X();
      case kY:
         return fP.Y();
      case kZ:
         return fP.Z();
      case kT:
         return fE;
      default:
         Error("operator()()", "bad index (%d) returning 0",i);
   }
   return 0.;
}

inline Double_t & TLorentzVector::operator () (int i)
{
   //dereferencing operator
   switch(i) {
      case kX:
         return fP.fX;
      case kY:
         return fP.fY;
      case kZ:
         return fP.fZ;
      case kT:
         return fE;
      default:
         Error("operator()()", "bad index (%d) returning &fE",i);
   }
   return fE;
}

#endif
