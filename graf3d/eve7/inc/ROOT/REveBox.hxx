// @(#)root/eve:$Id$
// Author: Matevz Tadel, 2010

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_REveBox
#define ROOT_REveBox

#include "ROOT/REveShape.hxx"

namespace ROOT {
namespace Experimental {

//------------------------------------------------------------------------------
// REveBox
//------------------------------------------------------------------------------

class REveBox : public REveShape,
                public REveProjectable
{
private:
   REveBox(const REveBox&) = delete;
   REveBox& operator=(const REveBox&) = delete;

protected:
   Float_t fVertices[8][3];

public:
   REveBox(const char* n="REveBox", const char* t="");
   virtual ~REveBox();

   void SetVertex(Int_t i, Float_t x, Float_t y, Float_t z);
   void SetVertex(Int_t i, const Float_t* v);
   void SetVertices(const Float_t* vs);

   const Float_t* GetVertex(Int_t i) const { return fVertices[i]; }

   // For TAttBBox:
   virtual void ComputeBBox() override;

   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset) override;
   void BuildRenderData() override;
   // Projectable:
   virtual TClass* ProjectedClass(const REveProjection* p) const override;
};


//------------------------------------------------------------------------------
// REveBoxProjected
//------------------------------------------------------------------------------

class REveBoxProjected : public REveShape,
                         public REveProjected
{
private:
   REveBoxProjected(const REveBoxProjected&) = delete;
   REveBoxProjected& operator=(const REveBoxProjected&) = delete;

protected:
   vVector2_t   fPoints;
   Int_t        fBreakIdx;
   vVector2_t   fDebugPoints;

   Bool_t       fDebugCornerPoints;

   virtual void SetDepthLocal(Float_t d) override;

public:
   REveBoxProjected(const char* n="REveBoxProjected", const char* t="");
   virtual ~REveBoxProjected();

   void BuildRenderData() override;
   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset) override;

   // For TAttBBox:
   void ComputeBBox() override;

   Bool_t GetDebugCornerPoints();
   void   SetDebugCornerPoints(Bool_t d);

   // Projected:
   void SetProjection(REveProjectionManager* mng, REveProjectable* model) override;
   void UpdateProjection() override;

   REveElement* GetProjectedAsElement() override { return this; }
};

} // namespace Experimental
} // namespace ROOT
#endif
