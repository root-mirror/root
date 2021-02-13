// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveLine
#define ROOT7_REveLine

#include <ROOT/REvePointSet.hxx>
#include <ROOT/REveVector.hxx>

#include "TAttLine.h"

namespace ROOT {
namespace Experimental {

////////////////////////////////////////////////////////////////////////////////
/// REveLine
/// An arbitrary polyline with fixed line and marker attributes.
////////////////////////////////////////////////////////////////////////////////

class REveLine : public REvePointSet,
                 public TAttLine
{
private:
   REveLine &operator=(const REveLine &) = delete;

protected:
   Bool_t fRnrLine;
   Bool_t fRnrPoints;
   Bool_t fSmooth;

   static Bool_t fgDefaultSmooth;

public:
   REveLine(const std::string &name = "", const std::string &title = "", Int_t n_points = 0);
   REveLine(const REveLine &l);
   virtual ~REveLine() {}

   void SetMarkerColor(Color_t col) override;

   void SetLineColor(Color_t col) override { SetMainColor(col); }
   void SetLineStyle(Style_t lstyle) override;
   void SetLineWidth(Width_t lwidth) override;

   Bool_t GetRnrLine() const { return fRnrLine; }
   Bool_t GetRnrPoints() const { return fRnrPoints; }
   Bool_t GetSmooth() const { return fSmooth; }
   void   SetRnrLine(Bool_t r);
   void   SetRnrPoints(Bool_t r);
   void   SetSmooth(Bool_t r);

   void    ReduceSegmentLengths(Float_t max);
   Float_t CalculateLineLength() const;

   REveVector GetLineStart() const;
   REveVector GetLineEnd() const;

   void CopyVizParams(const REveElement *el) override;
   void WriteVizParams(std::ostream &out, const TString &var) override;

   TClass *ProjectedClass(const REveProjection *p) const override;

   Int_t WriteCoreJson(nlohmann::json &cj, Int_t rnr_offset) override;
   void BuildRenderData() override;

   static Bool_t GetDefaultSmooth();
   static void SetDefaultSmooth(Bool_t r);
};

//------------------------------------------------------------------------------
// REveLineProjected
//------------------------------------------------------------------------------

class REveLineProjected : public REveLine, public REveProjected {
private:
   REveLineProjected(const REveLineProjected &) = delete;
   REveLineProjected &operator=(const REveLineProjected &) = delete;

protected:
   void SetDepthLocal(Float_t d) override;

public:
   REveLineProjected();
   virtual ~REveLineProjected() {}

   void SetProjection(REveProjectionManager *mng, REveProjectable *model) override;
   void UpdateProjection() override;
   REveElement *GetProjectedAsElement() override { return this; }
};

} // namespace Experimental
} // namespace ROOT

#endif
