// @(#)root/geom:$Id$
// Author: Andrei Gheata   2003/04/10

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TVirtualGeoTrack
#define ROOT_TVirtualGeoTrack

#include "TObjArray.h"

#include "TGeoAtt.h"

#include "TAttLine.h"

#include "TAttMarker.h"

class TVirtualGeoTrack : public TObject,
                         public TGeoAtt,
                         public TAttLine,
                         public TAttMarker
{
protected:
   Int_t             fPDG{0};              // track pdg code
   Int_t             fId{-1};              // track id
   TVirtualGeoTrack *fParent{nullptr};     // id of parent
   TObject          *fParticle{nullptr};   // particle for this track
   TObjArray        *fTracks{nullptr};     // daughter tracks

   TVirtualGeoTrack(const TVirtualGeoTrack&) = delete;
   TVirtualGeoTrack& operator=(const TVirtualGeoTrack&) = delete;

public:
   TVirtualGeoTrack();
   TVirtualGeoTrack(Int_t id, Int_t pdgcode, TVirtualGeoTrack *parent=0, TObject *particle=0);
   virtual ~TVirtualGeoTrack();

   virtual TVirtualGeoTrack *AddDaughter(Int_t id, Int_t pdgcode, TObject *particle=0) = 0;
   virtual Int_t       AddDaughter(TVirtualGeoTrack *other) = 0;
   virtual void        AddPoint(Double_t x, Double_t y, Double_t z, Double_t t) = 0;
   virtual TVirtualGeoTrack *FindTrackWithId(Int_t id) const;
   Int_t               GetId() const         {return fId;}
   virtual Int_t       GetDaughterId(Int_t index) const;
   TVirtualGeoTrack   *GetDaughter(Int_t index) const {return (TVirtualGeoTrack*)fTracks->At(index);}
   TVirtualGeoTrack   *GetMother() const {return fParent;}
   TObject            *GetMotherParticle() const {return (fParent)?fParent->GetParticle():0;}
   virtual const char *GetName() const;
   Int_t               GetNdaughters() const {return (fTracks)?fTracks->GetEntriesFast():0;}
   virtual Int_t       GetNpoints() const = 0;
   Int_t               GetParentId() const   {return (fParent)?fParent->GetId():-1;}
   TObject            *GetParticle() const   {return fParticle;}
   Int_t               GetPDG() const        {return fPDG;}
   Int_t               GetLastPoint(Double_t &x, Double_t &y, Double_t &z, Double_t &t) const {return GetPoint(GetNpoints()-1,x,y,z,t);}
   const Double_t     *GetFirstPoint() const {return GetPoint(0);}
   const Double_t     *GetLastPoint() const {return GetPoint(GetNpoints()-1);}
   virtual Int_t       GetPoint(Int_t i, Double_t &x, Double_t &y, Double_t &z, Double_t &t) const = 0;
   virtual const Double_t *GetPoint(Int_t i) const = 0;
   Bool_t              HasPoints() const {return (GetNpoints()==0)?kFALSE:kTRUE;}
   Bool_t              IsInTimeRange() const;
   virtual void        Paint(Option_t *option="") = 0;
   virtual void        PaintCollect(Double_t /*time*/, Double_t * /*box*/) {;}
   virtual void        PaintCollectTrack(Double_t /*time*/, Double_t * /*box*/) {;}
   virtual void        PaintTrack(Option_t *option="") = 0;
   virtual void        ResetTrack() = 0;
   void                SetName(const char *name);
   virtual void        SetParticle(TObject *particle) {fParticle=particle;}
   void                SetParent(TVirtualGeoTrack *parent) {fParent = parent;}
   void                SetId(Int_t id)       {fId = id;}
   virtual void        SetPDG(Int_t pdgcode) {fPDG = pdgcode;}

   ClassDef(TVirtualGeoTrack, 1)              // virtual geometry tracks
};

#endif




