// @(#)root/geom:$Id$
// Author: Andrei Gheata   2003/04/10

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoTrack
#define ROOT_TGeoTrack

#include "TVirtualGeoTrack.h"

/////////////////////////////////////////////////////////////////////////////
// TGeoTrack - Tracks attached to a geometry.                              //
//             Tracks are 3D objects made of points and they store a       //
//             pointer to a TParticle. The geometry manager holds a list   //
//             of all tracks that will be deleted on destruction of        //
//             gGeoManager.                                                //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

class TGeoTrack : public TVirtualGeoTrack
{
public:

enum EGeoParticleActions {
   kGeoPDefault      = BIT(7),
   kGeoPOnelevel     = BIT(8),
   kGeoPAllDaughters = BIT(9),
   kGeoPType         = BIT(10),
   kGeoPDrawn        = BIT(11)
};

private :
   Int_t          fPointsSize; // capacity of points array
   Int_t          fNpoints;    // number of stored points
   Double_t      *fPoints;     //[fNpoints] array of points (x,y,z,t) belonging to this track

protected:
   TGeoTrack(const TGeoTrack&) = delete;
   TGeoTrack& operator=(const TGeoTrack&) = delete;

public:
   TGeoTrack();
   TGeoTrack(Int_t id, Int_t pdgcode, TVirtualGeoTrack *parent=0, TObject *particle=0);
   virtual ~TGeoTrack();

   virtual TVirtualGeoTrack *AddDaughter(Int_t id, Int_t pdgcode, TObject *particle=0);
   virtual Int_t       AddDaughter(TVirtualGeoTrack *other);
   virtual void        AddPoint(Double_t x, Double_t y, Double_t z, Double_t t);
   virtual void        AnimateTrack(Double_t tmin=0, Double_t tmax=5E-8, Double_t nframes=200, Option_t *option="/*"); // *MENU*
   void                Browse(TBrowser *b);
   virtual Int_t       DistancetoPrimitive(Int_t px, Int_t py);
   virtual void        Draw(Option_t *option=""); // *MENU*
   virtual void        ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual char       *GetObjectInfo(Int_t px, Int_t py) const;
   virtual Int_t       GetNpoints() const {return (fNpoints>>2);}
   virtual Int_t       GetPoint(Int_t i, Double_t &x, Double_t &y, Double_t &z, Double_t &t) const;
   virtual const Double_t *GetPoint(Int_t i) const;
   Int_t               GetPoint(Double_t tof, Double_t *point, Int_t istart=0) const;
   Bool_t              IsFolder() const {return (GetNdaughters()>0)?kTRUE:kFALSE;}
   virtual void        Paint(Option_t *option="");
   virtual void        PaintCollect(Double_t time, Double_t *box);
   virtual void        PaintCollectTrack(Double_t time, Double_t *box);
   void                PaintMarker(Double_t *point, Option_t *option="");
   virtual void        PaintTrack(Option_t *option="");
   virtual void        Print(Option_t *option="") const; // *MENU*
   virtual void        ResetTrack();
   Int_t               SearchPoint(Double_t time, Int_t istart=0) const;
   void                SetBits(Bool_t is_default=kTRUE, Bool_t is_onelevel=kFALSE,
                               Bool_t is_all=kFALSE, Bool_t is_type=kFALSE);
   Int_t               Size(Int_t &imin, Int_t &imax);
   virtual void        Sizeof3D() const;

   ClassDef(TGeoTrack, 1)              // geometry tracks class
};

#endif




