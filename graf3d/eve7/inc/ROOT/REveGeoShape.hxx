// @(#)root/eve7:$Id$
// Author: Matevz Tadel 2007, 1018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveGeoShape
#define ROOT7_REveGeoShape

#include <ROOT/REveShape.hxx>

class TGeoShape;
class TGeoHMatrix;
class TGeoCompositeShape;
class TBuffer3D;

namespace ROOT {
namespace Experimental {

class REveGeoShapeExtract;

// ==========================================================================================
// REveGeoShape
// Wrapper for TGeoShape with absolute positioning and color attributes allowing display of extracted
// TGeoShape's (without an active TGeoManager) and simplified geometries (needed for NLT projections).
// ==========================================================================================

class REveGeoShape : public REveShape,
                     public REveProjectable
{
private:
   REveGeoShape(const REveGeoShape &) = delete;
   REveGeoShape &operator=(const REveGeoShape &) = delete;

protected:
   Int_t fNSegments{0};
   TGeoShape *fShape{nullptr};
   TGeoCompositeShape *fCompositeShape{nullptr}; //! Temporary holder (if passed shape is composite shape).

   static TGeoManager *fgGeoManager;

   static REveGeoShape *SubImportShapeExtract(REveGeoShapeExtract *gse, REveElement *parent);
   REveGeoShapeExtract *DumpShapeTree(REveGeoShape *geon, REveGeoShapeExtract *parent = nullptr);

   TGeoShape *MakePolyShape();

public:
   REveGeoShape(const std::string &name = "REveGeoShape", const std::string &title = "");
   virtual ~REveGeoShape();

   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset) override;
   void BuildRenderData() override;

   Int_t GetNSegments() const { return fNSegments; }
   TGeoShape *GetShape() const { return fShape; }
   void SetNSegments(Int_t s);
   void SetShape(TGeoShape *s);

   void ComputeBBox() override;

   void SaveExtract(const char *file, const char *name);
   void WriteExtract(const char *name);

   static REveGeoShape *ImportShapeExtract(REveGeoShapeExtract *gse, REveElement *parent = nullptr);

   // GeoProjectable
   virtual std::unique_ptr<TBuffer3D> MakeBuffer3D();
   TClass *ProjectedClass(const REveProjection *p) const override;

   static TGeoManager *GetGeoManager();
   static TGeoHMatrix *GetGeoHMatrixIdentity();
};

//------------------------------------------------------------------------------

// ==========================================================================================
// REveGeoShapeProjected
// ==========================================================================================

class REveGeoShapeProjected : public REveShape, public REveProjected {
private:
   REveGeoShapeProjected(const REveGeoShapeProjected &) = delete;
   REveGeoShapeProjected &operator=(const REveGeoShapeProjected &) = delete;

protected:
   std::unique_ptr<TBuffer3D> fBuff;    //! 3d buffer

   void SetDepthLocal(Float_t d) override;

public:
   REveGeoShapeProjected();
   virtual ~REveGeoShapeProjected();

   void SetProjection(REveProjectionManager *proj, REveProjectable *model) override;
   void UpdateProjection() override;
   REveElement *GetProjectedAsElement() override { return this; }

   void ComputeBBox() override;
};

} // namespace Experimental
} // namespace ROOT

#endif
