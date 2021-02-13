// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveLine.hxx>
#include <ROOT/REveProjectionManager.hxx>
#include <ROOT/REveRenderData.hxx>

#include "TClass.h"

using namespace ROOT::Experimental;

/** \class REveLine
\ingroup REve
An arbitrary polyline with fixed line and marker attributes.
*/

Bool_t REveLine::fgDefaultSmooth = kFALSE;


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveLine::REveLine(const std::string &name, const std::string &title, Int_t n_points) :
   REvePointSet(name, title, n_points),
   fRnrLine   (kTRUE),
   fRnrPoints (kFALSE),
   fSmooth    (fgDefaultSmooth)
{
   fMainColorPtr = &fLineColor;
   fMarkerColor  = kGreen;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

REveLine::REveLine(const REveLine &l) :
   REvePointSet(l),
   TAttLine    (l),
   fRnrLine    (l.fRnrLine),
   fRnrPoints  (l.fRnrPoints),
   fSmooth     (l.fSmooth)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker color. Propagate to projected lines.

void REveLine::SetMarkerColor(Color_t col)
{
   for (auto &pi: fProjectedList)
   {
      REveLine* l = dynamic_cast<REveLine*>(pi);
      if (l && fMarkerColor == l->GetMarkerColor())
      {
         l->SetMarkerColor(col);
         l->StampObjProps();
      }
   }
   TAttMarker::SetMarkerColor(col);
}

////////////////////////////////////////////////////////////////////////////////
/// Set line-style of the line.
/// The style is propagated to projecteds.

void REveLine::SetLineStyle(Style_t lstyle)
{
   for (auto &pi: fProjectedList)
   {
      REveLine* pt = dynamic_cast<REveLine*>(pi);
      if (pt)
      {
         pt->SetLineStyle(lstyle);
         pt->StampObjProps();
      }
   }
   TAttLine::SetLineStyle(lstyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set line-style of the line.
/// The style is propagated to projecteds.

void REveLine::SetLineWidth(Width_t lwidth)
{
   for (auto &pi: fProjectedList)
   {
      REveLine* pt = dynamic_cast<REveLine*>(pi);
      if (pt)
      {
         pt->SetLineWidth(lwidth);
         pt->StampObjProps();
      }
   }
   StampObjProps();
   TAttLine::SetLineWidth(lwidth);
}

////////////////////////////////////////////////////////////////////////////////
/// Set rendering of line. Propagate to projected lines.

void REveLine::SetRnrLine(Bool_t r)
{
   fRnrLine = r;
   for (auto &pi: fProjectedList)
   {
      REveLine* l = dynamic_cast<REveLine*>(pi);
      if (l)
      {
         l->SetRnrLine(r);
         l->StampObjProps();
      }
   }
   StampObjProps();
}

////////////////////////////////////////////////////////////////////////////////
/// Set rendering of points. Propagate to projected lines.

void REveLine::SetRnrPoints(Bool_t r)
{
   fRnrPoints = r;
   for (auto &pi: fProjectedList)
   {
      REveLine *l = dynamic_cast<REveLine*>(pi);
      if (l)
      {
         l->SetRnrPoints(r);
         l->StampObjProps();
      }
   }
   StampObjProps();
}

////////////////////////////////////////////////////////////////////////////////
/// Set smooth rendering. Propagate to projected lines.

void REveLine::SetSmooth(Bool_t r)
{
   fSmooth = r;
   for (auto &pi: fProjectedList)
   {
      REveLine* l = dynamic_cast<REveLine*>(pi);
      if (l)
      {
         l->SetSmooth(r);
         l->StampObjProps();
      }
   }
   StampObjProps();
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that no segment is longer than max.
/// Per point references and integer ids are lost.

void REveLine::ReduceSegmentLengths(Float_t max)
{
   // XXXX rewrite

   const Float_t max2 = max*max;

   Float_t    *p = & fPoints[0].fX;
   Int_t       s = fSize;
   REveVector  a, b, d;

   std::vector<REveVector> q;

   b.Set(p);
   q.push_back(b);
   for (Int_t i = 1; i < s; ++i)
   {
      a = b; b.Set(&p[3*i]); d = b - a;
      Float_t m2 = d.Mag2();
      if (m2 > max2)
      {
         Float_t f = TMath::Sqrt(m2) / max;
         Int_t   n = TMath::FloorNint(f);
         d *= 1.0f / (n + 1);
         for (Int_t j = 0; j < n; ++j)
         {
            a += d;
            q.push_back(a);
         }
      }
      q.push_back(b);
   }

   s = q.size();
   Reset(s);
   for (auto &i: q)
      SetNextPoint(i.fX, i.fY, i.fZ);
}

////////////////////////////////////////////////////////////////////////////////
/// Sum-up lengths of individual segments.

Float_t REveLine::CalculateLineLength() const
{
   Float_t sum = 0;

   for (Int_t i = 1; i < fSize; ++i)
   {
      sum += fPoints[i - 1].Distance(fPoints[i]);
   }

   return sum;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the first point of the line.
/// If there are no points (0,0,0) is returned.

REveVector REveLine::GetLineStart() const
{
   REveVector v;
   if (fSize > 0) v = RefPoint(0);
   return v;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the last point of the line.
/// If there are no points (0,0,0) is returned.

REveVector REveLine::GetLineEnd() const
{
   REveVector v;
   if (fSize > 0) v = RefPoint(fSize - 1);
   return v;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy visualization parameters from element el.

void REveLine::CopyVizParams(const REveElement* el)
{
   const REveLine* m = dynamic_cast<const REveLine*>(el);
   if (m)
   {
      TAttLine::operator=(*m);
      fRnrLine   = m->fRnrLine;
      fRnrPoints = m->fRnrPoints;
      fSmooth    = m->fSmooth;
   }

   REvePointSet::CopyVizParams(el);
}

////////////////////////////////////////////////////////////////////////////////
/// Write visualization parameters.

void REveLine::WriteVizParams(std::ostream& out, const TString& var)
{
   REvePointSet::WriteVizParams(out, var);

   TString t = "   " + var + "->";
   TAttLine::SaveLineAttributes(out, var);
   out << t << "SetRnrLine("   << ToString(fRnrLine)   << ");\n";
   out << t << "SetRnrPoints(" << ToString(fRnrPoints) << ");\n";
   out << t << "SetSmooth("    << ToString(fSmooth)    << ");\n";
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from REveProjectable, returns REvePointSetProjected class.

TClass* REveLine::ProjectedClass(const REveProjection*) const
{
   return TClass::GetClass<REveLineProjected>();
}

//------------------------------------------------------------------------------

Int_t REveLine::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   Int_t ret = REvePointSet::WriteCoreJson(j, rnr_offset);

   j["fLineWidth"] = GetLineWidth();
   j["fLineStyle"] = GetLineStyle();
   j["fLineColor"] = GetLineColor();

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from REveElement. Prepares render data for binary streaming to client

void REveLine::BuildRenderData()
{
   if (fSize > 0)
   {
      fRenderData = std::make_unique<REveRenderData>("makeTrack", 3*fSize);
      fRenderData->PushV(&fPoints[0].fX, 3*fSize);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get default value for smooth-line drawing flag.
/// Static function.

Bool_t REveLine::GetDefaultSmooth()
{
   return fgDefaultSmooth;
}

////////////////////////////////////////////////////////////////////////////////
/// Set default value for smooth-line drawing flag (default kFALSE).
/// Static function.

void REveLine::SetDefaultSmooth(Bool_t r)
{
   fgDefaultSmooth = r;
}

/** \class REveLineProjected
\ingroup REve
Projected copy of a REveLine.
*/


////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

REveLineProjected::REveLineProjected() :
   REveLine      (),
   REveProjected ()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Set projection manager and projection model.
/// Virtual from REveProjected.

void REveLineProjected::SetProjection(REveProjectionManager* mng,
                                      REveProjectable* model)
{
   REveProjected::SetProjection(mng, model);
   CopyVizParams(dynamic_cast<REveElement*>(model));
}

////////////////////////////////////////////////////////////////////////////////
/// Set depth (z-coordinate) of the projected points.

void REveLineProjected::SetDepthLocal(Float_t d)
{
   SetDepthCommon(d, this, fBBox);

   Int_t    n = fSize;
   Float_t *p = & fPoints[0].fZ;
   for (Int_t i = 0; i < n; ++i, p+=3)
      *p = fDepth;
}

////////////////////////////////////////////////////////////////////////////////
/// Re-apply the projection.
/// Virtual from REveProjected.

void REveLineProjected::UpdateProjection()
{
   REveProjection& proj = * fManager->GetProjection();
   REveLine      & als  = * dynamic_cast<REveLine*>(fProjectable);
   REveTrans      *tr   =   als.PtrMainTrans(kFALSE);

   Int_t n = als.GetSize();
   Reset(n);
   fSize = n;
   const Float_t *o = & als.RefPoint(0).fX;
         Float_t *p = & fPoints[0].fX;
   for (Int_t i = 0; i < n; ++i, o+=3, p+=3)
   {
      proj.ProjectPointfv(tr, o, p, fDepth);
   }
}
