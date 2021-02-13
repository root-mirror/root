// @(#)root/eve7:$Id$
// Author: Matevz Tadel, 2010, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveShape
#define ROOT7_REveShape

#include <ROOT/REveElement.hxx>
#include <ROOT/REveVector.hxx>

#include "TAttBBox.h"
#include "TColor.h"

namespace ROOT {
namespace Experimental {

// =========================================================================
// REveShape
// Abstract base-class for 2D/3D shapes.
// =========================================================================

class REveShape : public REveElement,
                  public TAttBBox
{
private:
   REveShape(const REveShape &) = delete;
   REveShape &operator=(const REveShape &) = delete;

public:
   typedef std::vector<REveVector2> vVector2_t;

protected:
   Color_t fFillColor; // fill color of polygons
   Color_t fLineColor; // outline color of polygons
   Float_t fLineWidth; // outline width of polygons

   Bool_t fDrawFrame;      // draw frame
   Bool_t fHighlightFrame; // highlight frame / all shape
   Bool_t fMiniFrame;      // draw minimal frame

public:
   REveShape(const std::string &n = "REveShape", const std::string &t = "");
   virtual ~REveShape();

   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset) override;

   // Rendering parameters.
   void SetMainColor(Color_t color) override;

   virtual Color_t GetFillColor() const { return fFillColor; }
   virtual Color_t GetLineColor() const { return fLineColor; }
   virtual Float_t GetLineWidth() const { return fLineWidth; }
   virtual Bool_t GetDrawFrame() const { return fDrawFrame; }
   virtual Bool_t GetHighlightFrame() const { return fHighlightFrame; }
   virtual Bool_t GetMiniFrame() const { return fMiniFrame; }

   virtual void SetFillColor(Color_t c) { fFillColor = c; }
   virtual void SetLineColor(Color_t c) { fLineColor = c; }
   virtual void SetLineWidth(Float_t lw) { fLineWidth = lw; }
   virtual void SetDrawFrame(Bool_t f) { fDrawFrame = f; }
   virtual void SetHighlightFrame(Bool_t f) { fHighlightFrame = f; }
   virtual void SetMiniFrame(Bool_t r) { fMiniFrame = r; }

   // ----------------------------------------------------------------

   void CopyVizParams(const REveElement *el) override;
   void WriteVizParams(std::ostream &out, const TString &var) override;

   // ----------------------------------------------------------------

   // Abstract function from TAttBBox:
   // virtual void ComputeBBox();

   // ----------------------------------------------------------------

   static Int_t FindConvexHull(const vVector2_t &pin, vVector2_t &pout, REveElement *caller = nullptr);

   static Bool_t IsBoxOrientationConsistentEv(const REveVector box[8]);
   static Bool_t IsBoxOrientationConsistentFv(const Float_t box[8][3]);

   static void CheckAndFixBoxOrientationEv(REveVector box[8]);
   static void CheckAndFixBoxOrientationFv(Float_t box[8][3]);
};

} // namespace Experimental
} // namespace ROOT

#endif
