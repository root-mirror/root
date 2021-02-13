// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_REveStraightLineSet
#define ROOT_REveStraightLineSet

#include "TNamed.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TAttBBox.h"

#include "ROOT/REveElement.hxx"
#include "ROOT/REveProjectionBases.hxx"
#include "ROOT/REveChunkManager.hxx"
#include "ROOT/REveTrans.hxx"
#include "ROOT/REveSecondarySelectable.hxx"

class TRandom;

namespace ROOT {
namespace Experimental {

///////////////////////////////////////////////////////////////////////////////
/// REveStraightLineSet
/// Set of straight lines with optional markers along the lines.
///////////////////////////////////////////////////////////////////////////////

class REveStraightLineSet : public REveElement,
                            public REveProjectable,
                            public TAttLine,
                            public TAttMarker,
                            public TAttBBox,
                            public REveSecondarySelectable
{
private:
   REveStraightLineSet(const REveStraightLineSet&) = delete;
   REveStraightLineSet& operator=(const REveStraightLineSet&) = delete;

public:
   struct Line_t
   {
      Int_t          fId;
      Float_t        fV1[3];
      Float_t        fV2[3];

      Line_t(Float_t x1, Float_t y1, Float_t z1,
             Float_t x2, Float_t y2, Float_t z2) : fId(-1)
      {
         fV1[0] = x1; fV1[1] = y1; fV1[2] = z1;
         fV2[0] = x2; fV2[1] = y2; fV2[2] = z2;
      }
   };

   struct Marker_t
   {
      Float_t      fV[3];
      Int_t        fLineId;

      Marker_t(Float_t x, Float_t y, Float_t z, Int_t line_id) : fLineId(line_id)
      {
         fV[0] = x; fV[1] = y; fV[2] = z;
      }
   };

protected:
   REveChunkManager  fLinePlex;
   REveChunkManager  fMarkerPlex;

   Bool_t            fOwnLinesIds;    // Flag specifying if id-objects are owned by the line-set
   Bool_t            fOwnMarkersIds;  // Flag specifying if id-objects are owned by the line-set

   Bool_t            fRnrMarkers;
   Bool_t            fRnrLines;

   Bool_t            fDepthTest;

   Line_t           *fLastLine{nullptr}; ///<!

public:
   REveStraightLineSet(const std::string &n="StraightLineSet", const std::string &t="");
   virtual ~REveStraightLineSet() {}

   void SetLineColor(Color_t col) override { SetMainColor(col); }

   Line_t   *AddLine(Float_t x1, Float_t y1, Float_t z1, Float_t x2, Float_t y2, Float_t z2);
   Line_t   *AddLine(const REveVector& p1, const REveVector& p2);
   Marker_t *AddMarker(Float_t x, Float_t y, Float_t z, Int_t line_id=-1);
   Marker_t *AddMarker(const REveVector& p, Int_t line_id=-1);
   Marker_t *AddMarker(Int_t line_id, Float_t pos);

   void      SetLine(int idx, Float_t x1, Float_t y1, Float_t z1, Float_t x2, Float_t y2, Float_t z2);
   void      SetLine(int idx, const REveVector& p1, const REveVector& p2);

   REveChunkManager &GetLinePlex()   { return fLinePlex;   }
   REveChunkManager &GetMarkerPlex() { return fMarkerPlex; }

   virtual Bool_t GetRnrMarkers() { return fRnrMarkers; }
   virtual Bool_t GetRnrLines()   { return fRnrLines;   }
   virtual Bool_t GetDepthTest()  { return fDepthTest;   }

   virtual void SetRnrMarkers(Bool_t x) { fRnrMarkers = x; }
   virtual void SetRnrLines(Bool_t x)   { fRnrLines   = x; }
   virtual void SetDepthTest(Bool_t x)  { fDepthTest   = x; }

   void CopyVizParams(const REveElement* el) override;
   void WriteVizParams(std::ostream& out, const TString& var) override;

   TClass* ProjectedClass(const REveProjection* p) const override;

   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset) override;
   void  BuildRenderData() override;

   void ComputeBBox() override;
};


///////////////////////////////////////////////////////////////////////////////
/// REveStraightLineSetProjected
/// Projected copy of a REveStraightLineSet.
///////////////////////////////////////////////////////////////////////////////

class REveStraightLineSetProjected : public REveStraightLineSet,
                                     public REveProjected
{
private:
   REveStraightLineSetProjected(const REveStraightLineSetProjected&) = delete;
   REveStraightLineSetProjected& operator=(const REveStraightLineSetProjected&) = delete;

protected:
   void SetDepthLocal(Float_t d) override;

public:
   REveStraightLineSetProjected();
   virtual ~REveStraightLineSetProjected() {}

   void SetProjection(REveProjectionManager* mng, REveProjectable* model) override;
   void UpdateProjection() override;
   REveElement* GetProjectedAsElement() override { return this; }
};

} // namespace Experimental
} // namespace ROOT

#endif
