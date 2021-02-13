// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveTrackProjected
#define ROOT7_REveTrackProjected

#include <ROOT/REveTrack.hxx>
#include <ROOT/REveProjectionBases.hxx>

namespace ROOT {
namespace Experimental {

////////////////////////////////////////////////////////////////////////////////
/// REveTrackProjected
/// Projected copy of a REveTrack.
////////////////////////////////////////////////////////////////////////////////

class REveTrackProjected : public REveTrack,
                           public REveProjected
{
private:
   REveTrackProjected(const REveTrackProjected &) = delete;
   REveTrackProjected &operator=(const REveTrackProjected &) = delete;

   Int_t GetBreakPointIdx(Int_t start);

   REveVector *fOrigPnts{nullptr}; // original track points

protected:
   std::vector<UInt_t> fBreakPoints; // indices of track break-points

   void SetDepthLocal(Float_t d) override;

public:
   REveTrackProjected() = default;
   virtual ~REveTrackProjected();

   void SetProjection(REveProjectionManager *mng, REveProjectable *model) override;

   void UpdateProjection() override;
   REveElement *GetProjectedAsElement() override { return this; }
   void MakeTrack(Bool_t recurse = kTRUE) override;

   void PrintLineSegments();

   void SecSelected(REveTrack *) override; // marked as signal in REveTrack

   Int_t WriteCoreJson(nlohmann::json &cj, Int_t rnr_offset) override;
   void BuildRenderData() override;
};

////////////////////////////////////////////////////////////////////////////////
/// REveTrackListProjected
/// Specialization of REveTrackList for holding REveTrackProjected objects.
////////////////////////////////////////////////////////////////////////////////

class REveTrackListProjected : public REveTrackList, public REveProjected {
private:
   REveTrackListProjected(const REveTrackListProjected &) = delete;
   REveTrackListProjected &operator=(const REveTrackListProjected &) = delete;

protected:
   void SetDepthLocal(Float_t d) override;

public:
   REveTrackListProjected();
   virtual ~REveTrackListProjected() {}

   void SetProjection(REveProjectionManager *proj, REveProjectable *model) override;
   void UpdateProjection() override {}
   REveElement *GetProjectedAsElement() override { return this; }

   void SetDepth(Float_t d) override;
   virtual void SetDepth(Float_t d, REveElement *el);
};

} // namespace Experimental
} // namespace ROOT

#endif
