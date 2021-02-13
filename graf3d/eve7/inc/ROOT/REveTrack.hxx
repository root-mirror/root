// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveTrack
#define ROOT7_REveTrack

#include <vector>

#include <ROOT/REveVector.hxx>
#include <ROOT/REvePathMark.hxx>
#include <ROOT/REveElement.hxx>
#include <ROOT/REveLine.hxx>
#include <ROOT/REveVSDStructs.hxx>

#include "TPolyMarker3D.h"
#include "TMarker.h"

class TParticle;

namespace ROOT {
namespace Experimental {

class REveTrackPropagator;
class REveTrackList;

////////////////////////////////////////////////////////////////////////////////
/// REveTrack
/// Track with given vertex, momentum and optional referece-points (path-marks) along its path.
////////////////////////////////////////////////////////////////////////////////

class REveTrack : public REveLine
{
   friend class REveTrackPropagator;
   friend class REveTrackList;

private:
   REveTrack &operator=(const REveTrack &) = delete;

public:
   typedef std::vector<REvePathMarkD> vPathMark_t;

   // Deprecated -- to be removed.
   enum EBreakProjectedTracks_e { kBPTDefault, kBPTAlways, kBPTNever };

protected:
   REveVectorD fV;         // Starting vertex
   REveVectorD fP;         // Starting momentum
   REveVectorD fPEnd;      // Momentum at the last point of extrapolation
   Double_t fBeta;         // Relativistic beta factor
   Double_t fDpDs;         // Momentum loss over distance
   Int_t fPdg;             // PDG code
   Int_t fCharge;          // Charge in units of e0
   Int_t fLabel;           // Simulation label
   Int_t fIndex;           // Reconstruction index
   Int_t fStatus;          // Status-word, user-defined.
   Bool_t fLockPoints;     // Lock points that are currently in - do nothing in MakeTrack().
   vPathMark_t fPathMarks; // REveVector of known points along the track
   Int_t fLastPMIdx;       //! Last path-mark index tried in track-propagation.

   REveTrackPropagator *fPropagator{nullptr}; // Pointer to shared render-style

public:
   REveTrack();
   REveTrack(TParticle *t, Int_t label, REveTrackPropagator *prop = nullptr);
   // VSD inputs
   REveTrack(REveMCTrack*  t, REveTrackPropagator* prop=0);
   REveTrack(REveRecTrack* t, REveTrackPropagator* prop=0);
   REveTrack(REveRecTrackD* t, REveTrackPropagator* prop=0);
   REveTrack(const REveTrack &t);
   virtual ~REveTrack();

   void ComputeBBox() override;

   virtual void SetStdTitle();

   virtual void SetTrackParams(const REveTrack &t);
   virtual void SetPathMarks(const REveTrack &t);

   virtual void MakeTrack(Bool_t recurse = kTRUE);

   REveTrackPropagator *GetPropagator() const { return fPropagator; }
   Int_t GetLastPMIdx() const { return fLastPMIdx; }
   void SetPropagator(REveTrackPropagator *prop);
   void SetAttLineAttMarker(REveTrackList *tl);

   const REveVectorD &GetVertex() const { return fV; }
   const REveVectorD &GetMomentum() const { return fP; }
   const REveVectorD &GetEndMomentum() const { return fPEnd; }

   Double_t GetDpDs() const { return fDpDs; }
   void SetDpDs(Double_t dpds) { fDpDs = dpds; }

   Int_t GetPdg() const { return fPdg; }
   void SetPdg(Int_t pdg) { fPdg = pdg; }
   Int_t GetCharge() const { return fCharge; }
   void SetCharge(Int_t chg) { fCharge = chg; }
   Int_t GetLabel() const { return fLabel; }
   void SetLabel(Int_t lbl) { fLabel = lbl; }
   Int_t GetIndex() const { return fIndex; }
   void SetIndex(Int_t idx) { fIndex = idx; }
   Int_t GetStatus() const { return fStatus; }
   void SetStatus(Int_t idx) { fStatus = idx; }

   void AddPathMark(const REvePathMarkD &pm) { fPathMarks.push_back(pm); }
   void AddPathMark(const REvePathMark &pm) { fPathMarks.push_back(pm); }

   void SortPathMarksByTime();
   vPathMark_t &RefPathMarks() { return fPathMarks; }
   const vPathMark_t &RefPathMarks() const { return fPathMarks; }

   void PrintPathMarks(); // *MENU*

   void SetLockPoints(Bool_t l) { fLockPoints = l; }
   Bool_t GetLockPoints() const { return fLockPoints; }

   //-------------------------------------------------------------------

   virtual void SecSelected(REveTrack *); // *SIGNAL*

   void CopyVizParams(const REveElement *el) override;
   void WriteVizParams(std::ostream &out, const TString &var) override;

   TClass *ProjectedClass(const REveProjection *p) const override;

   Int_t WriteCoreJson(nlohmann::json &cj, Int_t rnr_offset) override;
   void BuildRenderData() override;
};

////////////////////////////////////////////////////////////////////////////////
/// REveTrackList
/// A list of tracks supporting change of common attributes and selection based on track parameters.
////////////////////////////////////////////////////////////////////////////////

class REveTrackList : public REveElement,
                      public REveProjectable,
                      public TAttMarker,
                      public TAttLine
{
private:
   REveTrackList(const REveTrackList &) = delete;
   REveTrackList &operator=(const REveTrackList &) = delete;

protected:
   REveTrackPropagator *fPropagator{nullptr}; // Basic track rendering parameters, not enforced to elements.

   Bool_t fRecurse; // Recurse when propagating marker/line/etc attributes to tracks.

   Bool_t fRnrLine;   // Render track as line.
   Bool_t fRnrPoints; // Render track as points.

   Double_t fMinPt; // Minimum track pTfor display selection.
   Double_t fMaxPt; // Maximum track pTfor display selection.
   Double_t fLimPt; // Highest track pT in the container.
   Double_t fMinP;  // Minimum track pfor display selection.
   Double_t fMaxP;  // Maximum track pfor display selection.
   Double_t fLimP;  // Highest track p in the container.

   void FindMomentumLimits(REveElement *el, Bool_t recurse = kTRUE);
   Double_t RoundMomentumLimit(Double_t x);
   void SanitizeMinMaxCuts();

public:
   REveTrackList(REveTrackPropagator *prop = nullptr);
   REveTrackList(const std::string &name, REveTrackPropagator *prop = nullptr);
   virtual ~REveTrackList();

   void MakeTracks(Bool_t recurse = kTRUE);
   void FindMomentumLimits(Bool_t recurse = kTRUE);

   void SetPropagator(REveTrackPropagator *prop);
   REveTrackPropagator *GetPropagator() { return fPropagator; }

   Bool_t GetRecurse() const { return fRecurse; }
   void SetRecurse(Bool_t x) { fRecurse = x; }

   //--------------------------------

   void SetMainColor(Color_t c) override;
   void SetLineColor(Color_t c) override { SetMainColor(c); }
   virtual void SetLineColor(Color_t c, REveElement *el);
   void SetLineWidth(Width_t w) override;
   virtual void SetLineWidth(Width_t w, REveElement *el);
   void SetLineStyle(Style_t s) override;
   virtual void SetLineStyle(Style_t s, REveElement *el);

   void SetMarkerColor(Color_t c) override;
   virtual void SetMarkerColor(Color_t c, REveElement *el);
   void SetMarkerSize(Size_t s) override;
   virtual void SetMarkerSize(Size_t s, REveElement *el);
   void SetMarkerStyle(Style_t s) override;
   virtual void SetMarkerStyle(Style_t s, REveElement *el);

   void SetRnrLine(Bool_t rnr);
   void SetRnrLine(Bool_t rnr, REveElement *el);
   Bool_t GetRnrLine() const { return fRnrLine; }

   void SetRnrPoints(Bool_t r);
   void SetRnrPoints(Bool_t r, REveElement *el);
   Bool_t GetRnrPoints() const { return fRnrPoints; }

   void SelectByPt(Double_t min_pt, Double_t max_pt);
   void SelectByPt(Double_t min_pt, Double_t max_pt, REveElement *el);
   void SelectByP(Double_t min_p, Double_t max_p);
   void SelectByP(Double_t min_p, Double_t max_p, REveElement *el);

   Double_t GetMinPt() const { return fMinPt; }
   Double_t GetMaxPt() const { return fMaxPt; }
   Double_t GetLimPt() const { return fLimPt; }
   Double_t GetMinP() const { return fMinP; }
   Double_t GetMaxP() const { return fMaxP; }
   Double_t GetLimP() const { return fLimP; }

   //-------------------------------------------------------------------

   REveTrack *FindTrackByLabel(Int_t label); // *MENU*
   REveTrack *FindTrackByIndex(Int_t index); // *MENU*

   void CopyVizParams(const REveElement *el) override;
   void WriteVizParams(std::ostream &out, const TString &var) override;

   TClass *ProjectedClass(const REveProjection *p) const override;
};

} // namespace Experimental
} // namespace ROOT

#endif
