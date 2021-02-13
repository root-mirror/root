// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_REveBoxSet
#define ROOT_REveBoxSet

#include "ROOT/REveDigitSet.hxx"
#include "ROOT/REveVector.hxx"

class TGeoMatrix;
class TRandom;

namespace ROOT {
namespace Experimental {
class REveBoxSet: public REveDigitSet
{
   friend class REveBoxSetGL;

   REveBoxSet(const REveBoxSet&) = delete;
   REveBoxSet& operator=(const REveBoxSet&) = delete;

public:
   enum EBoxType_e {
      kBT_Undef,           // unknown-ignored
      kBT_FreeBox,         // arbitrary box: specify 8*(x,y,z) box corners
      kBT_AABox,           // axis-aligned box: specify (x,y,z) and (w, h, d)
      kBT_AABoxFixedDim,   // axis-aligned box w/ fixed dimensions: specify (x,y,z)
      kBT_Cone,
      kBT_EllipticCone,
      kBT_Hex
   };

   struct BFreeBox_t       : public DigitBase_t { Float_t fVertices[8][3]; };

   struct BOrigin_t        : public DigitBase_t { Float_t fA, fB, fC; };

   struct BAABox_t         : public BOrigin_t   { Float_t fW, fH, fD; };

   struct BAABoxFixedDim_t : public BOrigin_t   {};

   struct BCone_t          : public DigitBase_t { REveVector fPos, fDir; Float_t fR; };

   struct BEllipticCone_t  : public BCone_t     { Float_t fR2, fAngle; };

   struct BHex_t           : public DigitBase_t { REveVector fPos; Float_t fR, fAngle, fDepth; };

protected:
   EBoxType_e        fBoxType;      // Type of rendered box.

   Float_t           fDefWidth;     // Breadth assigned to first coordinate  (A).
   Float_t           fDefHeight;    // Breadth assigned to second coordinate (B).
   Float_t           fDefDepth;     // Breadth assigned to third coordinate  (C).

   Int_t             fBoxSkip;      // Number of boxes to skip for each drawn box during scene rotation.

   Bool_t            fDrawConeCap;

   static Int_t SizeofAtom(EBoxType_e bt);

public:
   REveBoxSet(const char* n="REveBoxSet", const char* t="");
   virtual ~REveBoxSet() {}

   void Reset(EBoxType_e boxType, Bool_t valIsCol, Int_t chunkSize);
   void Reset();

   void AddBox(const Float_t* verts);
   void AddBox(Float_t a, Float_t b, Float_t c, Float_t w, Float_t h, Float_t d);
   void AddBox(Float_t a, Float_t b, Float_t c);

   void AddCone(const REveVector& pos, const REveVector& dir, Float_t r);
   void AddEllipticCone(const REveVector& pos, const REveVector& dir, Float_t r, Float_t r2, Float_t angle=0);

   void AddHex(const REveVector& pos, Float_t r, Float_t angle, Float_t depth);

   void ComputeBBox() override;

   void Test(Int_t nboxes);

   Float_t GetDefWidth()  const { return fDefWidth;  }
   Float_t GetDefHeight() const { return fDefHeight; }
   Float_t GetDefDepth()  const { return fDefDepth;  }
   Bool_t  GetDrawConeCap() const { return fDrawConeCap;  }

   void SetDefWidth(Float_t v)  { fDefWidth  = v ; }
   void SetDefHeight(Float_t v) { fDefHeight = v ; }
   void SetDefDepth(Float_t v)  { fDefDepth  = v ; }
   void SetDrawConeCap(Bool_t x) { fDrawConeCap=x; StampObjProps(); }

   Int_t GetBoxSkip()   const { return fBoxSkip; }
   void  SetBoxSkip(Int_t bs) { fBoxSkip = bs; }


   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset) override;
   void  BuildRenderData() override;
};

} // namespace Experimental
} // namespace ROOT
#endif
