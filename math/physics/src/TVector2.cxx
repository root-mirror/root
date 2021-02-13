// @(#)root/physics:$Id$
// Author: Pasha Murat   12/02/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TVector2
    \ingroup Physics

TVector2 is a general two vector class, which can be used for
the description of different vectors in 2D.
*/

#include "TROOT.h"
#include "TVector2.h"
#include "TMath.h"
#include "TBuffer.h"


Double_t const  kPI        = TMath::Pi();
Double_t const  kTWOPI     = 2.*kPI;

ClassImp(TVector2);

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TVector2::TVector2()
{
   fX = 0.;
   fY = 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TVector2::TVector2(Double_t *v)
{
   fX = v[0];
   fY = v[1];
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TVector2::TVector2(Double_t x0, Double_t y0)
{
   fX = x0;
   fY = y0;
}

////////////////////////////////////////////////////////////////////////////////

TVector2::~TVector2()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Return modulo of this vector

Double_t TVector2::Mod() const
{
   return TMath::Sqrt(fX*fX+fY*fY);
}

////////////////////////////////////////////////////////////////////////////////
/// Return module normalized to 1

TVector2 TVector2::Unit() const
{
   return (Mod2()) ? *this/Mod() : TVector2();
}

////////////////////////////////////////////////////////////////////////////////
/// Return vector phi

Double_t TVector2::Phi() const
{
   return TMath::Pi()+TMath::ATan2(-fY,-fX);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns phi angle in the interval [0,2*PI)

Double_t TVector2::Phi_0_2pi(Double_t x) {
   if(TMath::IsNaN(x)){
      gROOT->Error("TVector2::Phi_0_2pi","function called with NaN");
      return x;
   }
   while (x >= kTWOPI) x -= kTWOPI;
   while (x <     0.)  x += kTWOPI;
   return x;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns phi angle in the interval [-PI,PI)

Double_t TVector2::Phi_mpi_pi(Double_t x) {
   if(TMath::IsNaN(x)){
      gROOT->Error("TVector2::Phi_mpi_pi","function called with NaN");
      return x;
   }
   while (x >= kPI) x -= kTWOPI;
   while (x < -kPI) x += kTWOPI;
   return x;
}

////////////////////////////////////////////////////////////////////////////////
/// Rotation by phi

TVector2 TVector2::Rotate (Double_t phi) const
{
   return TVector2( fX*TMath::Cos(phi)-fY*TMath::Sin(phi), fX*TMath::Sin(phi)+fY*TMath::Cos(phi) );
}

////////////////////////////////////////////////////////////////////////////////
/// Set vector using mag and phi

void TVector2::SetMagPhi(Double_t mag, Double_t phi)
{
   Double_t amag = TMath::Abs(mag);
   fX = amag * TMath::Cos(phi);
   fY = amag * TMath::Sin(phi);
}
////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TVector2.

void TVector2::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TVector2::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      if (R__v < 2) TObject::Streamer(R__b);
      R__b >> fX;
      R__b >> fY;
      R__b.CheckByteCount(R__s, R__c, TVector2::IsA());
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TVector2::Class(),this);
   }
}

void TVector2::Print(Option_t*)const
{
   //print vector parameters
   Printf("%s %s (x,y)=(%f,%f) (rho,phi)=(%f,%f)",GetName(),GetTitle(),X(),Y(),
                                          Mod(),Phi()*TMath::RadToDeg());
}
