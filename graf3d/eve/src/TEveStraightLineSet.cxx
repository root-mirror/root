// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveStraightLineSet.h"

#include "TEveProjectionManager.h"

/** \class TEveStraightLineSet
\ingroup TEve
Set of straight lines with optional markers along the lines.
*/

ClassImp(TEveStraightLineSet);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveStraightLineSet::TEveStraightLineSet(const char* n, const char* t):
   TEveElement (),
   TNamed      (n, t),

   fLinePlex      (sizeof(Line_t), 4),
   fMarkerPlex    (sizeof(Marker_t), 8),
   fOwnLinesIds   (kFALSE),
   fOwnMarkersIds (kFALSE),
   fRnrMarkers    (kTRUE),
   fRnrLines      (kTRUE),
   fDepthTest     (kTRUE),
   fLastLine      (0)
{
   InitMainTrans();
   fPickable = kTRUE;

   fMainColorPtr = &fLineColor;
   fLineColor    = 4;
   fMarkerColor  = 2;
   fMarkerStyle  = 20;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a line.

TEveStraightLineSet::Line_t*
TEveStraightLineSet::AddLine(Float_t x1, Float_t y1, Float_t z1,
                             Float_t x2, Float_t y2, Float_t z2)
{
   fLastLine = new (fLinePlex.NewAtom()) Line_t(x1, y1, z1, x2, y2, z2);
   fLastLine->fId = fLinePlex.Size() - 1;
   return fLastLine;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a line.

TEveStraightLineSet::Line_t*
TEveStraightLineSet::AddLine(const TEveVector& p1, const TEveVector& p2)
{
   return AddLine(p1.fX, p1.fY, p1.fZ, p2.fX, p2.fY, p2.fZ);
}

////////////////////////////////////////////////////////////////////////////////
/// Set line vertices with given index.

void
TEveStraightLineSet::SetLine(int idx,
                             Float_t x1, Float_t y1, Float_t z1,
                             Float_t x2, Float_t y2, Float_t z2)
{
   Line_t* l = (Line_t*) fLinePlex.Atom(idx);

   l->fV1[0] = x1; l->fV1[1] = y1; l->fV1[2] = z1;
   l->fV2[0] = x2; l->fV2[1] = y2; l->fV2[2] = z2;
}

////////////////////////////////////////////////////////////////////////////////
/// Set line vertices with given index.

void
TEveStraightLineSet::SetLine(int idx, const TEveVector& p1, const TEveVector& p2)
{
   SetLine(idx, p1.fX, p1.fY, p1.fZ, p2.fX, p2.fY, p2.fZ);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a marker with given position.

TEveStraightLineSet::Marker_t*
TEveStraightLineSet::AddMarker(Float_t x, Float_t y, Float_t z, Int_t line_id)
{
   Marker_t* marker = new (fMarkerPlex.NewAtom()) Marker_t(x, y, z, line_id);
   return marker;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a marker with given position.

TEveStraightLineSet::Marker_t*
TEveStraightLineSet::AddMarker(const TEveVector& p, Int_t line_id)
{
   return AddMarker(p.fX, p.fY, p.fZ, line_id);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a marker for line with given index on relative position pos.

TEveStraightLineSet::Marker_t*
TEveStraightLineSet::AddMarker(Int_t line_id, Float_t pos)
{
   Line_t& l = * (Line_t*) fLinePlex.Atom(line_id);
   return AddMarker(l.fV1[0] + (l.fV2[0] - l.fV1[0])*pos,
                    l.fV1[1] + (l.fV2[1] - l.fV1[1])*pos,
                    l.fV1[2] + (l.fV2[2] - l.fV1[2])*pos,
                    line_id);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy visualization parameters from element el.

void TEveStraightLineSet::CopyVizParams(const TEveElement* el)
{
   const TEveStraightLineSet* m = dynamic_cast<const TEveStraightLineSet*>(el);
   if (m)
   {
      TAttLine::operator=(*m);
      TAttMarker::operator=(*m);
      fRnrMarkers = m->fRnrMarkers;
      fRnrLines   = m->fRnrLines;
      fDepthTest  = m->fDepthTest;
   }

   TEveElement::CopyVizParams(el);
}

////////////////////////////////////////////////////////////////////////////////
/// Write visualization parameters.

void TEveStraightLineSet::WriteVizParams(std::ostream& out, const TString& var)
{
   TEveElement::WriteVizParams(out, var);

   TString t = "   " + var + "->";
   TAttMarker::SaveMarkerAttributes(out, var);
   TAttLine  ::SaveLineAttributes  (out, var);
   out << t << "SetRnrMarkers(" << ToString(fRnrMarkers) << ");\n";
   out << t << "SetRnrLines("   << ToString(fRnrLines)   << ");\n";
   out << t << "SetDepthTest("  << ToString(fDepthTest)  << ");\n";
}

////////////////////////////////////////////////////////////////////////////////
/// Return class of projected object.
/// Virtual from TEveProjectable.

TClass* TEveStraightLineSet::ProjectedClass(const TEveProjection*) const
{
   return TEveStraightLineSetProjected::Class();
}

////////////////////////////////////////////////////////////////////////////////
/// Compute bounding-box.
/// Virtual from TAttBBox.

void TEveStraightLineSet::ComputeBBox()
{
   if (fLinePlex.Size() == 0 && fMarkerPlex.Size() == 0) {
      BBoxZero();
      return;
   }

   BBoxInit();

   TEveChunkManager::iterator li(fLinePlex);
   while (li.next()) {
      BBoxCheckPoint(((Line_t*)li())->fV1);
      BBoxCheckPoint(((Line_t*)li())->fV2);
   }

   TEveChunkManager::iterator mi(fMarkerPlex);
   while (mi.next())
   {
      BBoxCheckPoint(((Marker_t*)mi())->fV);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the line-set.

void TEveStraightLineSet::Paint(Option_t*)
{
   PaintStandard(this);
}

/** \class TEveStraightLineSetProjected
\ingroup TEve
Projected replica of a TEveStraightLineSet.
*/

ClassImp(TEveStraightLineSetProjected);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveStraightLineSetProjected::TEveStraightLineSetProjected() :
   TEveStraightLineSet(), TEveProjected ()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Set projection manager and model object.

void TEveStraightLineSetProjected::SetProjection(TEveProjectionManager* mng,
                                                 TEveProjectable* model)
{
   TEveProjected::SetProjection(mng, model);

   CopyVizParams(dynamic_cast<TEveElement*>(model));
}

////////////////////////////////////////////////////////////////////////////////
/// Set depth (z-coordinate) of the projected points.

void TEveStraightLineSetProjected::SetDepthLocal(Float_t d)
{
   SetDepthCommon(d, this, fBBox);

   TEveChunkManager::iterator li(fLinePlex);
   while (li.next())
   {
      TEveStraightLineSet::Line_t& l = * (TEveStraightLineSet::Line_t*) li();
      l.fV1[2] = fDepth;
      l.fV2[2] = fDepth;
   }

   TEveChunkManager::iterator mi(fMarkerPlex);
   while (mi.next())
   {
      Marker_t& m = * (Marker_t*) mi();
      m.fV[2] = fDepth;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Callback that actually performs the projection.
/// Called when projection parameters have been updated.

void TEveStraightLineSetProjected::UpdateProjection()
{
   TEveProjection&      proj = * fManager->GetProjection();
   TEveStraightLineSet& orig = * dynamic_cast<TEveStraightLineSet*>(fProjectable);

   TEveTrans *trans = orig.PtrMainTrans(kFALSE);

   BBoxClear();

   // Lines
   Int_t num_lines = orig.GetLinePlex().Size();
   if (proj.HasSeveralSubSpaces())
      num_lines += TMath::Max(1, num_lines/10);
   fLinePlex.Reset(sizeof(Line_t), num_lines);
   TEveVector p1, p2;
   TEveChunkManager::iterator li(orig.GetLinePlex());
   while (li.next())
   {
      Line_t *l = (Line_t*) li();

      proj.ProjectPointfv(trans, l->fV1, p1, fDepth);
      proj.ProjectPointfv(trans, l->fV2, p2, fDepth);

      if (proj.AcceptSegment(p1, p2, 0.1f))
      {
         AddLine(p1, p2)->fId = l->fId;
      }
      else
      {
         TEveVector bp1(l->fV1), bp2(l->fV2);
         if (trans) {
            trans->MultiplyIP(bp1);
            trans->MultiplyIP(bp2);
         }
         proj.BisectBreakPoint(bp1, bp2, kTRUE, fDepth);

         AddLine(p1, bp1)->fId = l->fId;
         AddLine(bp2, p2)->fId = l->fId;
      }
   }
   if (proj.HasSeveralSubSpaces())
      fLinePlex.Refit();

   // Markers
   fMarkerPlex.Reset(sizeof(Marker_t), orig.GetMarkerPlex().Size());
   TEveChunkManager::iterator mi(orig.GetMarkerPlex());
   TEveVector pp;
   while (mi.next())
   {
      Marker_t &m = * (Marker_t*) mi();

      proj.ProjectPointfv(trans, m.fV, pp, fDepth);
      AddMarker(pp, m.fLineId);
   }
}
