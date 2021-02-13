// @(#)root/eve:$Id$
// Author: Matevz Tadel, 2010

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveShape
#define ROOT_TEveShape

#include "TEveElement.h"
#include "TEveVector.h"

#include "TAtt3D.h"
#include "TAttBBox.h"
#include "TColor.h"
#include <vector>

class TEveShape : public TEveElementList,
                  public TAtt3D,
                  public TAttBBox
{
   friend class TEveShapeEditor;

private:
   TEveShape(const TEveShape&) = delete;
   TEveShape& operator=(const TEveShape&) = delete;

public:
   typedef std::vector<TEveVector2>           vVector2_t;
   typedef std::vector<TEveVector2>::iterator vVector2_i;

protected:
   Color_t      fFillColor; // fill color of polygons
   Color_t      fLineColor; // outline color of polygons
   Float_t      fLineWidth; // outline width of polygons

   Bool_t       fDrawFrame;      // draw frame
   Bool_t       fHighlightFrame; // highlight frame / all shape
   Bool_t       fMiniFrame;      // draw minimal frame

public:
   TEveShape(const char* n="TEveShape", const char* t="");
   virtual ~TEveShape();

   // Rendering parameters.
   virtual void    SetMainColor(Color_t color);

   virtual Color_t GetFillColor() const { return fFillColor; }
   virtual Color_t GetLineColor() const { return fLineColor; }
   virtual Float_t GetLineWidth() const { return fLineWidth;}
   virtual Bool_t  GetDrawFrame()      const { return fDrawFrame; }
   virtual Bool_t  GetHighlightFrame() const { return fHighlightFrame; }
   virtual Bool_t  GetMiniFrame()      const { return fMiniFrame; }

   virtual void    SetFillColor(Color_t c)  { fFillColor = c; }
   virtual void    SetLineColor(Color_t c)  { fLineColor = c; }
   virtual void    SetLineWidth(Float_t lw) { fLineWidth = lw;}
   virtual void    SetDrawFrame(Bool_t f)      { fDrawFrame = f; }
   virtual void    SetHighlightFrame(Bool_t f) { fHighlightFrame = f; }
   virtual void    SetMiniFrame(Bool_t r)      { fMiniFrame = r; }

   // ----------------------------------------------------------------

   virtual void CopyVizParams(const TEveElement* el);
   virtual void WriteVizParams(std::ostream& out, const TString& var);

   // ----------------------------------------------------------------

   // Virtual from TObject
   virtual void Paint(Option_t* option="");

   // Abstract function from TAttBBox:
   // virtual void ComputeBBox();

   // Abstract from TEveProjectable, overriden in TEveElementList:
   // virtual TClass* ProjectedClass(const TEveProjection* p) const;

   // ----------------------------------------------------------------

   static Int_t  FindConvexHull(const vVector2_t& pin, vVector2_t& pout, TEveElement* caller=0);

   static Bool_t IsBoxOrientationConsistentEv(const TEveVector box[8]);
   static Bool_t IsBoxOrientationConsistentFv(const Float_t    box[8][3]);

   static void   CheckAndFixBoxOrientationEv(TEveVector box[8]);
   static void   CheckAndFixBoxOrientationFv(Float_t    box[8][3]);

   ClassDef(TEveShape, 0); // Abstract base-class for 2D/3D shapes.
};

#endif
