// @(#)root/eve7:$Id$
// Author: Sergey Linev, 27.02.2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveGeoPainter
#define ROOT7_REveGeoPainter

#include "TVirtualGeoPainter.h"

#include <ROOT/REveGeomViewer.hxx>

namespace ROOT {
namespace Experimental {

class REveGeoPainter : public TVirtualGeoPainter {

   TGeoManager *fGeoManager{nullptr};

   std::shared_ptr<REveGeomViewer> fViewer;

public:
   REveGeoPainter(TGeoManager *manager);
   virtual ~REveGeoPainter();

   void       AddSize3D(Int_t, Int_t, Int_t) override {}
   TVirtualGeoTrack *AddTrack(Int_t, Int_t, TObject *) override { return nullptr; }
   void       AddTrackPoint(Double_t *, Double_t *, Bool_t =kFALSE) override {}
   void       BombTranslation(const Double_t *, Double_t *) override {}
   void       CheckPoint(Double_t =0, Double_t =0, Double_t =0, Option_t * ="") override {}
   void       CheckShape(TGeoShape *, Int_t, Int_t, Option_t *) override {}
   void       CheckBoundaryErrors(Int_t =1000000, Double_t =-1.) override {}
   void       CheckBoundaryReference(Int_t =-1) override {}
   void       CheckGeometryFull(Bool_t =kTRUE, Bool_t =kTRUE, Int_t =10000, const Double_t * = nullptr) override {}
   void       CheckGeometry(Int_t, Double_t, Double_t, Double_t) const override {}
   void       CheckOverlaps(const TGeoVolume *, Double_t =0.1, Option_t * ="") const override {}
   Int_t      CountVisibleNodes() override { return 0; }
   void       DefaultAngles() override {}
   void       DefaultColors() override {}
   Int_t      DistanceToPrimitiveVol(TGeoVolume *, Int_t, Int_t) override { return 0; }
   void       Draw(Option_t * ="") override {}
   void       DrawBatemanSol(TGeoBatemanSol *, Option_t * ="") override {}
   void       DrawShape(TGeoShape *, Option_t * ="") override {}
   void       DrawOnly(Option_t * ="") override {}
   void       DrawOverlap(void *, Option_t * ="") override {}
   void       DrawCurrentPoint(Int_t) override {}
   void       DrawPanel() override {}
   void       DrawPath(const char *, Option_t * ="") override {}
   void       DrawPolygon(const TGeoPolygon *) override {}
   void       DrawVolume(TGeoVolume *, Option_t * ="") override;
   void       EditGeometry(Option_t * ="") override {}
   void       EstimateCameraMove(Double_t /*tmin*/, Double_t /*tmax*/, Double_t *, Double_t * ) override {}
   void       ExecuteShapeEvent(TGeoShape *, Int_t, Int_t, Int_t) override {}
   void       ExecuteManagerEvent(TGeoManager *, Int_t, Int_t, Int_t) override {}
   void       ExecuteVolumeEvent(TGeoVolume *, Int_t, Int_t, Int_t) override {}
   Int_t      GetColor(Int_t, Float_t) const override { return 0; }
   Int_t      GetNsegments() const override { return 1; }
   void       GetBombFactors(Double_t &, Double_t &, Double_t &, Double_t &) const override {}
   Int_t      GetBombMode() const override { return 0; }
   const char *GetDrawPath() const override { return ""; }
   TGeoVolume *GetDrawnVolume() const override { return nullptr; }
   TGeoVolume *GetTopVolume() const override { return nullptr; }
   void       GetViewAngles(Double_t &/*longitude*/, Double_t &/*latitude*/, Double_t &/*psi*/) override {}
   Int_t      GetVisLevel() const override { return 0; }
   Int_t      GetVisOption() const override { return 0; }
   const char *GetVolumeInfo(const TGeoVolume *, Int_t, Int_t) const override { return "info"; }
   void       GrabFocus(Int_t =0, Double_t =0, Double_t =0, Double_t =0) override {}
   Double_t  *GetViewBox() override { return nullptr; }
   Bool_t     IsPaintingShape() const override { return kFALSE; }
   Bool_t     IsRaytracing() const override { return kFALSE; }
   Bool_t     IsExplodedView() const override { return kFALSE; }
   TH2F      *LegoPlot(Int_t =60, Double_t =0., Double_t =180.,
                       Int_t =90, Double_t =0., Double_t =360.,
                       Double_t =0., Double_t =9999999, Option_t * ="") override { return nullptr; }
   void       ModifiedPad(Bool_t =kFALSE) const override {}
   void       OpProgress(const char *, Long64_t, Long64_t, TStopwatch * =nullptr, Bool_t =kFALSE, Bool_t =kFALSE, const char * ="") override {}
   void       Paint(Option_t * ="") override {}
   void       PaintNode(TGeoNode *, Option_t * ="", TGeoMatrix * = nullptr) override {}
   void       PaintShape(TGeoShape *, Option_t * ="") override {}
   void       PaintOverlap(void*, Option_t * ="") override {}
   void       PrintOverlaps() const override {}
   void       PaintVolume(TGeoVolume *, Option_t * = "", TGeoMatrix * = nullptr) override {}
   void       RandomPoints(const TGeoVolume *, Int_t, Option_t * = "") override {}
   void       RandomRays(Int_t, Double_t, Double_t, Double_t, const char *, Bool_t) override {}
   void       Raytrace(Option_t* = "") override {}
   TGeoNode  *SamplePoints(Int_t, Double_t &, Double_t, const char*) override { return nullptr; }
   void       SetBombFactors(Double_t =1.3, Double_t =1.3, Double_t =1.3, Double_t =1.3) override {}
   void       SetClippingShape(TGeoShape *) override {}
   void       SetExplodedView(Int_t =0) override {}
   void       SetGeoManager(TGeoManager *) override;
   void       SetIteratorPlugin(TGeoIteratorPlugin *) override {}
   void       SetCheckedNode(TGeoNode *) override {}
   void       SetNsegments(Int_t =20) override {}
   void       SetNmeshPoints(Int_t) override {}
   void       SetRaytracing(Bool_t =kTRUE) override {}
   void       SetTopVisible(Bool_t =kTRUE) override {}
   void       SetTopVolume(TGeoVolume *) override {}
   void       SetVisLevel(Int_t =3) override {}
   void       SetVisOption(Int_t =0) override {}
   Int_t      ShapeDistancetoPrimitive(const TGeoShape *, Int_t, Int_t, Int_t) const override { return 0; }
   void       Test(Int_t, Option_t *) override {}
   void       TestOverlaps(const char *) override {}
   Bool_t     TestVoxels(TGeoVolume *) override { return kFALSE; }
   void       UnbombTranslation(const Double_t *, Double_t *) override {}
   Double_t   Weight(Double_t, Option_t* = "v") override { return 0.; }

   ClassDefOverride(REveGeoPainter,0)  // Web-based geo painter
};

} // namespace Experimental
} // namespace ROOT


#endif
