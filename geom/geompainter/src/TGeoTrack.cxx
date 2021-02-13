// @(#)root/geom:$Id$
// Author: Andrei Gheata  2003/04/10

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TBrowser.h"
#include "TPoint.h"
#include "TVirtualPad.h"
#include "TVirtualX.h"
#include "TView.h"
#include "TGeoManager.h"
#include "TVirtualGeoPainter.h"
#include "TGeoTrack.h"

/** \class TGeoTrack
\ingroup Geometry_classes

Class for user-defined tracks attached to a geometry.
Tracks are 3D objects made of points and they store a
pointer to a TParticle. The geometry manager holds a list
of all tracks that will be deleted on destruction of
gGeoManager.
*/

ClassImp(TGeoTrack);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TGeoTrack::TGeoTrack()
{
   fPointsSize = 0;
   fNpoints    = 0;
   fPoints     = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGeoTrack::TGeoTrack(Int_t id, Int_t pdgcode, TVirtualGeoTrack *parent, TObject *particle)
          :TVirtualGeoTrack(id,pdgcode,parent,particle)
{
   fPointsSize = 0;
   fNpoints    = 0;
   fPoints     = 0;
   if (fParent==0) {
      SetMarkerColor(2);
      SetMarkerStyle(8);
      SetMarkerSize(0.6);
      SetLineColor(2);
      SetLineWidth(2);
   } else {
      SetMarkerColor(4);
      SetMarkerStyle(8);
      SetMarkerSize(0.6);
      SetLineColor(4);
      SetLineWidth(2);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoTrack::~TGeoTrack()
{
   if (fPoints) delete [] fPoints;
//   if (gPad) gPad->GetListOfPrimitives()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a daughter track to this.

TVirtualGeoTrack *TGeoTrack::AddDaughter(Int_t id, Int_t pdgcode, TObject *particle)
{
   if (!fTracks) fTracks = new TObjArray(1);
   Int_t index = fTracks->GetEntriesFast();
   TGeoTrack *daughter = new TGeoTrack(id,pdgcode,this,particle);
   fTracks->AddAtAndExpand(daughter,index);
   return daughter;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a daughter and return its index.

Int_t TGeoTrack::AddDaughter(TVirtualGeoTrack *other)
{
   if (!fTracks) fTracks = new TObjArray(1);
   Int_t index = fTracks->GetEntriesFast();
   fTracks->AddAtAndExpand(other,index);
   other->SetParent(this);
   return index;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw animation of this track

void TGeoTrack::AnimateTrack(Double_t tmin, Double_t tmax, Double_t nframes, Option_t *option)
{
   if (tmin<0 || tmin>=tmax || nframes<1) return;
   gGeoManager->SetAnimateTracks();
   gGeoManager->SetVisLevel(1);
   if (!gPad) {
      gGeoManager->GetMasterVolume()->Draw();
   }
   TList *list = gPad->GetListOfPrimitives();
   TIter next(list);
   TObject *obj;
   while ((obj = next())) {
      if (!strcmp(obj->ClassName(), "TGeoTrack")) list->Remove(obj);
   }
   Double_t dt = (tmax-tmin)/Double_t(nframes);
   Double_t delt = 2E-9;
   Double_t t = tmin;
   Bool_t geomanim = kFALSE;
   Bool_t issave = kFALSE;
   TString fname;

   TString opt(option);
   if (opt.Contains("/G")) geomanim = kTRUE;
   if (opt.Contains("/S")) issave = kTRUE;

   TVirtualGeoPainter *p = gGeoManager->GetGeomPainter();
   Double_t *box = p->GetViewBox();
   box[0] = box[1] = box[2] = 0;
   box[3] = box[4] = box[5] = 100;
   gGeoManager->SetTminTmax(0,0);
   Draw(opt.Data());
   Double_t start[6], end[6];
   Int_t i, j;
   Double_t dlat=0, dlong=0, dpsi=0;
   Double_t dd[6] = {0,0,0,0,0,0};
   if (geomanim) {
      p->EstimateCameraMove(tmin+5*dt, tmin+15*dt, start, end);
      for (i=0; i<3; i++) {
         start[i+3] = 20 + 1.3*start[i+3];
         end[i+3] = 20 + 0.9*end[i+3];
      }
      for (i=0; i<6; i++) {
         dd[i] = (end[i]-start[i])/10.;
      }
      memcpy(box, start, 6*sizeof(Double_t));
      p->GetViewAngles(dlong,dlat,dpsi);
      dlong = (-206-dlong)/Double_t(nframes);
      dlat  = (126-dlat)/Double_t(nframes);
      dpsi  = (75-dpsi)/Double_t(nframes);
      p->GrabFocus();
   }

   for (i=0; i<nframes; i++) {
      if (t-delt<0) gGeoManager->SetTminTmax(0,t);
      else gGeoManager->SetTminTmax(t-delt,t);
      if (geomanim) {
         for (j=0; j<6; j++) box[j]+=dd[j];
         p->GrabFocus(1,dlong,dlat,dpsi);
      } else {
         gPad->Modified();
         gPad->Update();
      }
      if (issave) {
         fname = TString::Format("anim%04d.gif", i);
         gPad->Print(fname);
      }
      t += dt;
   }
   gGeoManager->SetAnimateTracks(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a point on the track.

void TGeoTrack::AddPoint(Double_t x, Double_t y, Double_t z, Double_t t)
{
   if (!fPoints) {
      fPointsSize = 16;
      fPoints = new Double_t[fPointsSize];
   } else {
      if (fNpoints>=fPointsSize) {
         Double_t *temp = new Double_t[2*fPointsSize];
         memcpy(temp, fPoints, fNpoints*sizeof(Double_t));
         fPointsSize *= 2;
         delete [] fPoints;
         fPoints = temp;
      }
   }
   fPoints[fNpoints++] = x;
   fPoints[fNpoints++] = y;
   fPoints[fNpoints++] = z;
   fPoints[fNpoints++] = t;
}

////////////////////////////////////////////////////////////////////////////////
/// How-to-browse for a track.

void TGeoTrack::Browse(TBrowser *b)
{
   if (!b) return;
   Int_t nd = GetNdaughters();
   if (!nd) {
      b->Add(this);
      return;
   }
   for (Int_t i=0; i<nd; i++)
      b->Add(GetDaughter(i));

}

////////////////////////////////////////////////////////////////////////////////
/// Returns distance to track primitive for picking.

Int_t TGeoTrack::DistancetoPrimitive(Int_t px, Int_t py)
{
   const Int_t inaxis = 7;
   const Int_t maxdist = 5;
   Int_t dist = 9999;


   Int_t puxmin = gPad->XtoAbsPixel(gPad->GetUxmin());
   Int_t puymin = gPad->YtoAbsPixel(gPad->GetUymin());
   Int_t puxmax = gPad->XtoAbsPixel(gPad->GetUxmax());
   Int_t puymax = gPad->YtoAbsPixel(gPad->GetUymax());

   // return if point is not in the user area
   if (px < puxmin - inaxis) return dist;
   if (py > puymin + inaxis) return dist;
   if (px > puxmax + inaxis) return dist;
   if (py < puymax - inaxis) return dist;

   TView *view = gPad->GetView();
   if (!view) return dist;
   Int_t imin, imax;
   if (TObject::TestBit(kGeoPDrawn) && Size(imin,imax)>=2) {
      Int_t i, dsegment;
      Double_t x1,y1,x2,y2;
      Double_t xndc[3];
      Int_t np = fNpoints>>2;
      if (imin<0) imin=0;
      if (imax>np-1) imax=np-1;
      for (i=imin;i<imax;i++) {
         view->WCtoNDC(&fPoints[i<<2], xndc);
         x1 = xndc[0];
         y1 = xndc[1];
         view->WCtoNDC(&fPoints[(i+1)<<2], xndc);
         x2 = xndc[0];
         y2 = xndc[1];
         dsegment = DistancetoLine(px,py,x1,y1,x2,y2);
//         printf("%i: dseg=%i\n", i, dsegment);
         if (dsegment < dist) {
            dist = dsegment;
            if (dist<maxdist) {
               gPad->SetSelected(this);
               return 0;
            }
         }
      }
   }
   // check now daughters
   Int_t nd = GetNdaughters();
   if (!nd) return dist;
   TGeoTrack *track;
   for (Int_t id=0; id<nd; id++) {
      track = (TGeoTrack*)GetDaughter(id);
      dist = track->DistancetoPrimitive(px,py);
      if (dist<maxdist) return 0;
   }
   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this track over-imposed on a geometry, according to option.
/// Options (case sensitive):
///  - default : track without daughters
///  - /D      : track and first level descendents only
///  - /*      : track and all descendents
///  - /Ntype  : descendents of this track with particle name matching input type.
///
/// Options can appear only once but can be combined : e.g. Draw("/D /Npion-")
///
/// Time range for visible track segments can be set via TGeoManager::SetTminTmax()

void TGeoTrack::Draw(Option_t *option)
{
   if (!gPad) gGeoManager->GetMasterVolume()->Draw();
   char *opt1 = Compress(option); // we will have to delete this ?
   TString opt(opt1);
   Bool_t is_default = kTRUE;
   Bool_t is_onelevel = kFALSE;
   Bool_t is_all = kFALSE;
   Bool_t is_type = kFALSE;
   if (opt.Contains("/D")) {
      is_onelevel = kTRUE;
      is_default = kFALSE;
   }
   if (opt.Contains("/*")) {
      is_all = kTRUE;
      is_default = kFALSE;
   }
   if (opt.Contains("/N")) {
      is_type = kTRUE;
      Int_t ist = opt.Index("/N")+2;
      Int_t ilast = opt.Index("/",ist);
      if (ilast<0) ilast=opt.Length();
      TString type = opt(ist, ilast-ist);
      gGeoManager->SetParticleName(type.Data());
   }
   SetBits(is_default, is_onelevel, is_all, is_type);
   AppendPad("SAME");
   if (!gGeoManager->IsAnimatingTracks()) {
      gPad->Modified();
      gPad->Update();
   }
   delete [] opt1;
   return;
}

 ///////////////////////////////////////////////////////////////////////////////
 /// Event treatment.

void TGeoTrack::ExecuteEvent(Int_t /*event*/, Int_t /*px*/, Int_t /*py*/)
{
   if (!gPad) return;
   gPad->SetCursor(kHand);
}

////////////////////////////////////////////////////////////////////////////////
/// Get some info about the track.

char *TGeoTrack::GetObjectInfo(Int_t /*px*/, Int_t /*py*/) const
{
   static TString info;
   Double_t x=0,y=0,z=0,t=0;
   GetPoint(0,x,y,z,t);
   info = TString::Format("%s (%g, %g, %g) tof=%g", GetName(),x,y,z,t);
   return (char*)info.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Get coordinates for point I on the track.

Int_t TGeoTrack::GetPoint(Int_t i, Double_t &x, Double_t &y, Double_t &z, Double_t &t) const
{
   Int_t np = fNpoints>>2;
   if (i<0 || i>=np) {
      Error("GetPoint", "no point %i, indmax=%d", i, np-1);
      return -1;
   }
   Int_t icrt = 4*i;
   x = fPoints[icrt];
   y = fPoints[icrt+1];
   z = fPoints[icrt+2];
   t = fPoints[icrt+3];
   return i;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the pointer to the array of points starting with index I.

const Double_t *TGeoTrack::GetPoint(Int_t i) const
{
   if (!fNpoints) return 0;
   return (&fPoints[i<<2]);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the index of point on track having closest TOF smaller than
/// the input value. Output POINT is filled with the interpolated value.

Int_t TGeoTrack::GetPoint(Double_t tof, Double_t *point, Int_t istart) const
{
   Int_t np = fNpoints>>2;
   if (istart>(np-2)) return (np-1);
   Int_t ip = SearchPoint(tof, istart);
   if (ip<0 || ip>(np-2)) return ip;
   // point in segment (ip, ip+1) where 0<=ip<fNpoints-1
   Int_t i;
   Int_t j = ip<<2;
   Int_t k = (ip+1)<<2;
   Double_t dt  = tof-fPoints[j+3];
   Double_t ddt = fPoints[k+3]-fPoints[j+3];
   for (i=0; i<3; i++) point[i] = fPoints[j+i] +(fPoints[k+i]-fPoints[j+i])*dt/ddt;
   return ip;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this track (and descendents) with current attributes.

void TGeoTrack::Paint(Option_t *option)
{
   Bool_t is_default  = TObject::TestBit(kGeoPDefault);
   Bool_t is_onelevel = TObject::TestBit(kGeoPOnelevel);
   Bool_t is_all      = TObject::TestBit(kGeoPAllDaughters);
   Bool_t is_type     = TObject::TestBit(kGeoPType);
   Bool_t match_type  = kTRUE;
   TObject::SetBit(kGeoPDrawn, kFALSE);
   if (is_type) {
      const char *type = gGeoManager->GetParticleName();
      if (strlen(type) && strcmp(type, GetName())) match_type=kFALSE;
   }
   if (match_type) {
      if (is_default || is_onelevel || is_all) PaintTrack(option);
   }
   // paint now daughters
   Int_t nd = GetNdaughters();
   if (!nd || is_default) return;
   TGeoTrack *track;
   for (Int_t i=0; i<nd; i++) {
      track = (TGeoTrack*)GetDaughter(i);
      if (track->IsInTimeRange()) {
         track->SetBits(is_default,kFALSE,is_all,is_type);
         track->Paint(option);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint track and daughters.

void TGeoTrack::PaintCollect(Double_t time, Double_t *box)
{
   Bool_t is_default  = TObject::TestBit(kGeoPDefault);
   Bool_t is_onelevel = TObject::TestBit(kGeoPOnelevel);
   Bool_t is_all      = TObject::TestBit(kGeoPAllDaughters);
   Bool_t is_type     = TObject::TestBit(kGeoPType);
   Bool_t match_type  = kTRUE;
   if (is_type) {
      const char *type = gGeoManager->GetParticleName();
      if (strlen(type) && strcmp(type, GetName())) match_type=kFALSE;
   }
   if (match_type) {
      if (is_default || is_onelevel || is_all) PaintCollectTrack(time, box);
   }
   // loop now daughters
   Int_t nd = GetNdaughters();
   if (!nd || is_default) return;
   TGeoTrack *track;
   for (Int_t i=0; i<nd; i++) {
      track = (TGeoTrack*)GetDaughter(i);
      if (track) track->PaintCollect(time, box);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint just this track.

void TGeoTrack::PaintCollectTrack(Double_t time, Double_t *box)
{
   TVirtualGeoPainter *painter = gGeoManager->GetGeomPainter();
   if (!painter) return;
   Int_t np = fNpoints>>2;
   Double_t point[3], local[3];
   Bool_t convert = (gGeoManager->GetTopVolume() == gGeoManager->GetMasterVolume())?kFALSE:kTRUE;
   Int_t ip = GetPoint(time, point);
   if (ip>=0 && ip<np-1) {
      if (convert) gGeoManager->MasterToTop(point, local);
      else memcpy(local, point, 3*sizeof(Double_t));
      painter->AddTrackPoint(local, box);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint current point of the track as marker.

void TGeoTrack::PaintMarker(Double_t *point, Option_t *)
{
   TPoint p;
   Double_t xndc[3];
   TView *view = gPad->GetView();
   if (!view) return;
   view->WCtoNDC(point, xndc);
   if (xndc[0] < gPad->GetX1() || xndc[0] > gPad->GetX2()) return;
   if (xndc[1] < gPad->GetY1() || xndc[1] > gPad->GetY2()) return;
   p.fX = gPad->XtoPixel(xndc[0]);
   p.fY = gPad->YtoPixel(xndc[1]);
   TAttMarker::Modify();
   gVirtualX->DrawPolyMarker(1, &p);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this track with its current attributes.

void TGeoTrack::PaintTrack(Option_t *option)
{
   // Check whether there is some 3D view class for this TPad
//   TPadView3D *view3D = (TPadView3D*)gPad->GetView3D();
//   if (view3D) view3D->PaintGeoTrack(this,option); // to be implemented

   // Check if option is 'x3d'.      NOTE: This is a simple checking
   //                                      but since there is no other
   //                                      options yet, this works fine.
   TString opt(option);
   opt.ToLower();
   TObject::SetBit(kGeoPDrawn, kFALSE);
   if (opt.Contains("x")) return;
   Int_t np = fNpoints>>2;
   Int_t imin=0;
   Int_t imax=np-1;
   Int_t ip;
   Double_t start[3] = {0.,0.,0.};
   Double_t end[3] = {0.,0.,0.};
   Double_t seg[6] = {0.,0.,0.,0.,0.,0.};
   Bool_t convert = (gGeoManager->GetTopVolume() == gGeoManager->GetMasterVolume())?kFALSE:kTRUE;
   Double_t tmin=0.,tmax=0.;
   Bool_t is_time = gGeoManager->GetTminTmax(tmin,tmax);
   if (is_time) {
      imin = GetPoint(tmin, start);
      if (imin>=0 && imin<np-1) {
      // we have a starting point -> find ending point
         imax = GetPoint(tmax, end, imin);
         if (imax<np-1) {
         // we also have an ending point -> check if on the same segment with imin
            if (imax==imin) {
               // paint the virtual segment between the 2 points
               TAttLine::Modify();
               if (convert) {
                  gGeoManager->MasterToTop(start, &seg[0]);
                  gGeoManager->MasterToTop(end, &seg[3]);
                  gPad->PaintLine3D(&seg[0], &seg[3]);
               } else {
                  gPad->PaintLine3D(start, end);
               }
            } else {
               // paint the starting, ending and connecting segments
               TAttLine::Modify();
               if (convert) {
                  gGeoManager->MasterToTop(start, &seg[0]);
                  gGeoManager->MasterToTop(&fPoints[(imin+1)<<2], &seg[3]);
                  gPad->PaintLine3D(&seg[0], &seg[3]);
                  gGeoManager->MasterToTop(&fPoints[imax<<2], &seg[0]);
                  gGeoManager->MasterToTop(end, &seg[3]);
                  gPad->PaintLine3D(&seg[0], &seg[3]);
                  for (ip=imin+1; ip<imax; ip++) {
                     gGeoManager->MasterToTop(&fPoints[ip<<2], &seg[0]);
                     gGeoManager->MasterToTop(&fPoints[(ip+1)<<2], &seg[3]);
                     gPad->PaintLine3D(&seg[0], &seg[3]);
                  }
               } else {
                  gPad->PaintLine3D(start, &fPoints[(imin+1)<<2]);
                  gPad->PaintLine3D(&fPoints[imax<<2], end);
                  for (ip=imin+1; ip<imax; ip++) {
                     gPad->PaintLine3D(&fPoints[ip<<2], &fPoints[(ip+1)<<2]);
                  }
               }
            }
            if (convert) {
               gGeoManager->MasterToTop(end, &seg[0]);
               PaintMarker(&seg[0]);
            } else {
               PaintMarker(end);
            }
         } else {
            TAttLine::Modify();
            if (convert) {
               gGeoManager->MasterToTop(start, &seg[0]);
               gGeoManager->MasterToTop(&fPoints[(imin+1)<<2], &seg[3]);
               gPad->PaintLine3D(&seg[0], &seg[3]);
               for (ip=imin+1; ip<np-2; ip++) {
                  gGeoManager->MasterToTop(&fPoints[ip<<2], &seg[0]);
                  gGeoManager->MasterToTop(&fPoints[(ip+1)<<2], &seg[3]);
                  gPad->PaintLine3D(&seg[0], &seg[3]);
               }
            } else {
               gPad->PaintLine3D(start, &fPoints[(imin+1)<<2]);
               for (ip=imin+1; ip<np-2; ip++) {
                  gPad->PaintLine3D(&fPoints[ip<<2], &fPoints[(ip+1)<<2]);
               }
            }
         }
      } else {
         imax = GetPoint(tmax, end);
         if (imax<0 || imax>=(np-1)) return;
         // we have to draw just the end of the track
         TAttLine::Modify();
         if (convert) {
            for (ip=0; ip<imax-1; ip++) {
               gGeoManager->MasterToTop(&fPoints[ip<<2], &seg[0]);
               gGeoManager->MasterToTop(&fPoints[(ip+1)<<2], &seg[3]);
               gPad->PaintLine3D(&seg[0], &seg[3]);
            }
         } else {
            for (ip=0; ip<imax-1; ip++) {
               gPad->PaintLine3D(&fPoints[ip<<2], &fPoints[(ip+1)<<2]);
            }
         }
         if (convert) {
            gGeoManager->MasterToTop(&fPoints[imax<<2], &seg[0]);
            gGeoManager->MasterToTop(end, &seg[3]);
            gPad->PaintLine3D(&seg[0], &seg[3]);
            PaintMarker(&seg[3]);
         } else {
            gPad->PaintLine3D(&fPoints[imax<<2], end);
            PaintMarker(end);
         }
      }
      TObject::SetBit(kGeoPDrawn);
      return;
   }

   // paint all segments from track
   TObject::SetBit(kGeoPDrawn);
   TAttLine::Modify();  // change attributes if necessary
   for (ip=imin; ip<imax; ip++) {
      gPad->PaintLine3D(&fPoints[ip<<2], &fPoints[(ip+1)<<2]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print some info about the track.

void TGeoTrack::Print(Option_t * /*option*/) const
{
   Int_t np = fNpoints>>2;
   printf(" TGeoTrack%6i : %s  ===============================\n", fId,GetName());
   printf("   parent =%6i    nd =%3i\n", (fParent)?fParent->GetId():-1, GetNdaughters());
   Double_t x=0,y=0,z=0,t=0;
   GetPoint(0,x,y,z,t);
   printf("   production vertex : (%g, %g, %g) at tof=%g\n", x,y,z,t);
   GetPoint(np-1,x,y,z,t);
   printf("   Npoints =%6i,  last : (%g, %g, %g) at tof=%g\n\n", np,x,y,z,t);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of points within the time interval specified by
/// TGeoManager class and the corresponding indices.

Int_t TGeoTrack::Size(Int_t &imin, Int_t &imax)
{
   Double_t tmin, tmax;
   Int_t np = fNpoints>>2;
   imin = 0;
   imax = np-1;
   Int_t size = np;
   if (!gGeoManager->GetTminTmax(tmin, tmax)) return size;
   imin = SearchPoint(tmin);
   imax = SearchPoint(tmax, imin);
   return (imax-imin+1);
}

////////////////////////////////////////////////////////////////////////////////
/// Search index of track point having the closest time tag smaller than
/// TIME. Optional start index can be provided.

Int_t TGeoTrack::SearchPoint(Double_t time, Int_t istart) const
{
   Int_t nabove, nbelow, middle, midloc;
   Int_t np = fNpoints>>2;
   nabove = np+1;
   nbelow = istart;
   while (nabove-nbelow > 1) {
      middle = (nabove+nbelow)/2;
      midloc = ((middle-1)<<2)+3;
      if (time == fPoints[midloc]) return middle-1;
      if (time < fPoints[midloc])  nabove = middle;
      else                         nbelow = middle;
   }
   return (nbelow-1);
}

////////////////////////////////////////////////////////////////////////////////
/// Set drawing bits for this track

void TGeoTrack::SetBits(Bool_t is_default, Bool_t is_onelevel,
                        Bool_t is_all, Bool_t is_type)
{
   TObject::SetBit(kGeoPDefault, is_default);
   TObject::SetBit(kGeoPOnelevel, is_onelevel);
   TObject::SetBit(kGeoPAllDaughters, is_all);
   TObject::SetBit(kGeoPType, is_type);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns 3D size for the track.

void TGeoTrack::Sizeof3D() const
{
}

////////////////////////////////////////////////////////////////////////////////
/// Reset data for this track.

void TGeoTrack::ResetTrack()
{
   fNpoints    = 0;
   fPointsSize = 0;
   if (fTracks) {fTracks->Delete(); delete fTracks;}
   fTracks = 0;
   if (fPoints) delete [] fPoints;
   fPoints = 0;
}

