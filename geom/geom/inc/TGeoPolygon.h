// @(#)root/geom:$Id$
// Author: Mihaela Gheata   05/01/04

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoPolygon
#define ROOT_TGeoPolygon

#include "TObject.h"

class TObjArray;

class TGeoPolygon : public TObject
{
public:
  enum {
      kGeoConvex        = BIT(9),
      kGeoFinishPolygon = BIT(10),
      kGeoACW           = BIT(11)
   };
protected :
// data members
   Int_t               fNvert{0};                  // number of vertices (must be defined clockwise in XY plane)
   Int_t               fNconvex{0};                // number of points of the outscribed convex polygon
   Int_t              *fInd{nullptr};              //[fNvert] list of vertex indices
   Int_t              *fIndc{nullptr};             //[fNconvex] indices of vertices of the outscribed convex polygon
   Double_t           *fX{nullptr};                //! pointer to list of current X coordinates of vertices
   Double_t           *fY{nullptr};                //! pointer to list of current Y coordinates of vertices
   TObjArray          *fDaughters{nullptr};        // list of concave daughters
private:
   void                ConvexCheck(); // force convexity checking
   Bool_t              IsSegConvex(Int_t i1, Int_t i2=-1) const;
   Bool_t              IsRightSided(const Double_t *point, Int_t ind1, Int_t ind2) const;
   void                OutscribedConvex();

   TGeoPolygon(const TGeoPolygon&) = delete;
   TGeoPolygon &operator=(const TGeoPolygon&) = delete;
public:
   // constructors
   TGeoPolygon();
   TGeoPolygon(Int_t nvert);
   // destructor
   virtual ~TGeoPolygon();
   // methods
   Double_t            Area() const;
   Bool_t              Contains(const Double_t *point) const;
   virtual void        Draw(Option_t *option="");
   void                FinishPolygon();
   Int_t               GetNvert() const {return fNvert;}
   Int_t               GetNconvex() const {return fNconvex;}
   Double_t           *GetX() {return fX;}
   Double_t           *GetY() {return fY;}
   void                GetVertices(Double_t *x, Double_t *y) const;
   void                GetConvexVertices(Double_t *x, Double_t *y) const;
   Bool_t              IsClockwise() const {return !TObject::TestBit(kGeoACW);}
   Bool_t              IsConvex() const {return TObject::TestBit(kGeoConvex);}
   Bool_t              IsFinished() const {return TObject::TestBit(kGeoFinishPolygon);}
   Bool_t              IsIllegalCheck() const;
   Double_t            Safety(const Double_t *point, Int_t &isegment) const;
   void                SetConvex(Bool_t flag=kTRUE) {TObject::SetBit(kGeoConvex,flag);}
   void                SetXY(Double_t *x, Double_t *y);
   void                SetNextIndex(Int_t index=-1);

   ClassDef(TGeoPolygon, 2)         // class for handling arbitrary polygons
};

#endif

