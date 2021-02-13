// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveEllipsoid.hxx>
#include <ROOT/REveTrans.hxx>
#include <ROOT/REveProjectionManager.hxx>
#include <ROOT/REveRenderData.hxx>

#include "TMath.h"
#include "TClass.h"

#include <cassert>

using namespace ROOT::Experimental;


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveEllipsoid::REveEllipsoid(const std::string &n , const std::string &t):
   REveStraightLineSet(n, t)
{
   fPhiStep = 0.01f;
}


////////////////////////////////////////////////////////////////////////////////
/// Draw archade as straight line set.

void REveEllipsoid::DrawArch(float phiStart, float phiEnd, float phiStep, REveVector& v0,  REveVector& v1, REveVector& v2)
{
   float phi = phiStart;

   REveVector f =  v1;
   while (phi < phiEnd ) {
      REveVector v = v0 + v1*((float)cos(phi)) + v2*((float)sin(phi));
      AddLine(f, v);
      f=v;
      phi += phiStep;
   }
   REveVector v = v0 + v1*((float)cos(phiEnd)) + v2*((float)sin(phiEnd));
   AddLine(f, v);
}

////////////////////////////////////////////////////////////////////////////////
/// Set size of phi step in archade drawing.

void REveEllipsoid::SetPhiStep(float ps)
{
   fPhiStep = ps;
}

////////////////////////////////////////////////////////////////////////////////
/// Three defining base vectors of ellipse.

void REveEllipsoid::SetBaseVectors(REveVector& v0, REveVector& v1, REveVector& v2)
{
   fV0 = v0;
   fV1 = v1;
   fV2 = v2;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw archade around base vectors.

void REveEllipsoid::Outline()
{
   REveVector v0;
   DrawArch(0, TMath::TwoPi(),fPhiStep, v0, fV0, fV1);
   DrawArch(0, TMath::TwoPi(),fPhiStep, v0, fV0, fV2);
   DrawArch(0, TMath::TwoPi(),fPhiStep, v0, fV1, fV2);
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from REveProjectable, returns REveEllipsoidProjected class.

TClass* REveEllipsoid::ProjectedClass(const REveProjection*) const
{
   return TClass::GetClass<REveEllipsoidProjected>();
}


////////////////////////////////////////////////////////////////////////////////
/// Fill core part of JSON representation.

Int_t REveEllipsoid::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   Int_t ret = REveStraightLineSet::WriteCoreJson(j, rnr_offset);

   j["fSecondarySelect"] = false;
   // printf("REveStraightLineSet::WriteCoreJson %d \n", ret);
   return ret;
}

//==============================================================================
//==============================================================================
//==============================================================================

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveEllipsoidProjected::REveEllipsoidProjected(const std::string& /*n*/, const std::string& /*t*/) :
   REveStraightLineSetProjected()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REveEllipsoidProjected::~REveEllipsoidProjected()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Draw archade around base vectors.

void REveEllipsoidProjected::DrawArchProjected(float phiStart, float phiEnd, float phiStep, REveVector& v0,  REveVector& v1, REveVector& v2)
{
   float phi = phiStart;

   REveVector f =  v1;
   while (phi < phiEnd ) {
      REveVector v = v0 + v1*((float)cos(phi)) + v2*((float)sin(phi));
      fArchPnts.push_back(f);
      fArchPnts.push_back(v);
      f=v;
      phi += phiStep;
   }

   REveVector v = v0 + v1*((float)cos(phiEnd)) + v2*((float)sin(phiEnd));
   fArchPnts.push_back(f);
   fArchPnts.push_back(v);
}

////////////////////////////////////////////////////////////////////////////////
/// Get surface size of projected ellipse

float REveEllipsoidProjected::GetEllipseSurface (const REveVector& v1, const REveVector& v2)
{
   REveEllipsoid&     orig  = * dynamic_cast<REveEllipsoid*>(fProjectable);
   REveTrans *    trans = orig.PtrMainTrans(kFALSE);
   REveProjection&  proj  = * fManager->GetProjection();

   // project center of ellipse
   REveTrans trans0;
   TVector3 v0 = trans->GetPos();
   REveVector p0(v0.x(), v0.y(), v0.z());
   proj.ProjectPointfv(&trans0, p0, p0, fDepth);

   // first axis point
   REveVector p1 = v1;
   proj.ProjectPointfv(trans,v1, p1, fDepth);

   // second axis point
   REveVector p2 = v2;
   proj.ProjectPointfv(trans, v2, p2, fDepth);

   return  (p1-p0).Mag2()+ (p2-p0).Mag2();
}

////////////////////////////////////////////////////////////////////////////////
/// Find longest projection of axes and draw an arch.

void REveEllipsoidProjected::OutlineProjected()
{
   REveEllipsoid&     orig  = * dynamic_cast<REveEllipsoid*>(fProjectable);
   // find ellipse with biggest surface
   float max = 0;
   {
      REveVector v1 = orig.fV0;
      REveVector v2 = orig.fV1;
      float d = GetEllipseSurface(v1, v2);
      if (d > max) {
         fMV0 = v1;
         fMV1 = v2;
         max = d;
      }
   }
   {
      REveVector v1 = orig.fV1;
      REveVector v2 = orig.fV2;
      float d = GetEllipseSurface(v1, v2);
      if (d > max) {
         fMV0 = v1;
         fMV1 = v2;
         max = d;
      }
   }
   {
      REveVector v1 = orig.fV0;
      REveVector v2 = orig.fV2;
      float d = GetEllipseSurface(v1, v2);
      if (d > max) {
         fMV0 = v1;
         fMV1 = v2;
         max = d;
      }
   }
   if (gDebug) {
      printf("REveEllipsoidProjected::OutlineProjected, printing axes %s\n", GetCName());
      fMV0.Dump();
      fMV1.Dump();
   }

   REveVector p0;
   DrawArchProjected(0, TMath::TwoPi(), orig.fPhiStep, p0, fMV0, fMV1);
}

////////////////////////////////////////////////////////////////////////////////
/// Crates 3D point array for rendering.

void  REveEllipsoidProjected::BuildRenderData()
{
   REveStraightLineSetProjected::BuildRenderData();
}



////////////////////////////////////////////////////////////////////////////////
/// This is virtual method from base-class REveProjected.

void REveEllipsoidProjected::SetProjection(REveProjectionManager* mng, REveProjectable* model)
{
   REveProjected::SetProjection(mng, model);
   CopyVizParams(dynamic_cast<REveElement*>(model));
}

////////////////////////////////////////////////////////////////////////////////
/// Callback that actually performs the projection.
/// Called when projection parameters have been updated.

void REveEllipsoidProjected::UpdateProjection()
{
   OutlineProjected();
   REveProjection&      proj = * fManager->GetProjection();
   REveEllipsoid& orig = * dynamic_cast<REveEllipsoid*>(fProjectable);

   REveTrans *trans = orig.PtrMainTrans(kFALSE);

   // Lines
   Int_t num_lines = (int)fArchPnts.size();
   if (proj.HasSeveralSubSpaces())
      num_lines += TMath::Max(1, num_lines/10);
   fLinePlex.Reset(sizeof(Line_t), num_lines);
   REveVector p1, p2;
   for (size_t i = 0; i <fArchPnts.size(); i+=2 )
   {
      proj.ProjectPointfv(trans, fArchPnts[i], p1, fDepth);
      proj.ProjectPointfv(trans, fArchPnts[i+1], p2, fDepth);

      if (proj.AcceptSegment(p1, p2, 0.1f))
      {
         AddLine(p1, p2);
      }
      else
      {
         REveVector bp1(fArchPnts[i]), bp2(fArchPnts[i+1]);
         if (trans) {
            trans->MultiplyIP(bp1);
            trans->MultiplyIP(bp2);
         }
         proj.BisectBreakPoint(bp1, bp2, kTRUE, fDepth);

         AddLine(p1, bp1);
         AddLine(bp2, p2);
      }
   }
   if (proj.HasSeveralSubSpaces())
      fLinePlex.Refit();

   // Markers
   fMarkerPlex.Reset(sizeof(Marker_t), orig.GetMarkerPlex().Size());
   REveChunkManager::iterator mi(orig.GetMarkerPlex());
   REveVector pp;
   while (mi.next())
   {
      Marker_t &m = * (Marker_t*) mi();

      proj.ProjectPointfv(trans, m.fV, pp, fDepth);
      AddMarker(pp, m.fLineId);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill core part of JSON representation.

Int_t REveEllipsoidProjected::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   Int_t ret = REveStraightLineSet::WriteCoreJson(j, rnr_offset);

   j["fSecondarySelect"] = false;
   // printf("REveStraightLineSet::WriteCoreJson %d \n", ret);
   return ret;
}
