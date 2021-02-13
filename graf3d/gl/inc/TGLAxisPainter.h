// @(#)root/gl:$Id$
// Author: Alja Mrak-Tadel 2009

/*************************************************************************
* Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/

#ifndef ROOT_TGLAxisPainter
#define ROOT_TGLAxisPainter

#include "TGLUtil.h"
#include "TGLFontManager.h"

#include <vector>
#include <utility>

class TAttAxis;
class TAxis;
class TH1;
class TPMERegexp;
class TGLRnrCtx;


//==============================================================================
// TGLAxisPainter
//==============================================================================

class TGLAxisPainter
{
public:
   typedef std::pair  <Float_t, Float_t>    Lab_t; // label <pos, value> pair
   typedef std::vector<Lab_t>               LabVec_t;
   typedef std::pair  <Float_t, Int_t>      TM_t;  // tick-mark <pos, order> pair
   typedef std::vector<TM_t>                TMVec_t; // vector od tick lines

private:
   TGLAxisPainter(const TGLAxisPainter&);            // Not implemented
   TGLAxisPainter& operator=(const TGLAxisPainter&); // Not implemented

   // Print format derived from attributers.
   Int_t fExp;
   Int_t fMaxDigits;
   Int_t fDecimals;
   TString fFormat;

   // Font derived from axis attributes.
   TGLFont fLabelFont;
   TGLFont fTitleFont;

   // Print format.
   void LabelsLimits(const char *label, Int_t &first, Int_t &last) const;
   void FormAxisValue(Double_t x, TString &s) const;

protected:
   TAttAxis        *fAttAxis;    // Model.
   Bool_t           fUseAxisColors; // Use colors from axes or from GL-rnr-ctx.
   TGLFont::EMode   fFontMode;   // To be put into TAttAxis
   LabVec_t         fLabVec;     // List of Labels position-value pairs
   TMVec_t          fTMVec;      // List of tick-mark position-value pairs

   //
   // Additional axis attributes required for GL rendering:

   // Orientation
   TGLVector3 fDir;
   TGLVector3 fTMOff[3];
   Int_t      fTMNDim;

   // Font.
   Int_t    fLabelPixelFontSize;
   Double_t fLabel3DFontSize;
   Int_t    fTitlePixelFontSize;
   Double_t fTitle3DFontSize;

   // Labels options. Allready exist in TAttAxis, but can't be set.
   TGLFont::ETextAlignH_e fLabelAlignH;
   TGLFont::ETextAlignV_e fLabelAlignV;
   TGLVector3  fTitlePos;
   TPMERegexp *fAllZeroesRE;

public:
   TGLAxisPainter();
   virtual ~TGLAxisPainter();

   // GetSets.
   Int_t        GetTMNDim() const { return fTMNDim; }
   void         SetTMNDim(Int_t x) { fTMNDim = x; }

   TGLVector3&  RefDir() { return fDir; }
   TGLVector3&  RefTMOff(Int_t i) { return fTMOff[i]; }

   TGLFont::EMode GetFontMode() const { return fFontMode; }
   void  SetFontMode(TGLFont::EMode m) { fFontMode=m; }

   // this setter not necessary
   void         SetLabelPixelFontSize(Int_t fs) { fLabelPixelFontSize=fs; }
   Int_t        GetLabelPixelFontSize() const { return fLabelPixelFontSize; }
   void         SetTitlePixelFontSize(Int_t fs) { fTitlePixelFontSize=fs; }
   Int_t        GetTitlePixelFontSize() const { return fTitlePixelFontSize; }

   TGLVector3&  RefTitlePos() { return fTitlePos; }


   void         SetLabelAlign(TGLFont::ETextAlignH_e, TGLFont::ETextAlignV_e);

   LabVec_t& RefLabVec() { return fLabVec; }
   TMVec_t&  RefTMVec()  { return fTMVec; }

   void      SetAttAxis(TAttAxis* a) { fAttAxis = a; }
   TAttAxis* GetAttAxis() { return fAttAxis; }

   void   SetUseAxisColors(Bool_t x) { fUseAxisColors = x;    }
   Bool_t GetUseAxisColors() const   { return fUseAxisColors; }

   // Utility.
   void SetLabelFont(TGLRnrCtx &rnrCtx, const char* fontName, Int_t pixelSize = 64, Double_t font3DSize = -1);
   void SetTitleFont(TGLRnrCtx &rnrCtx, const char* fontName, Int_t pixelSize = 64, Double_t font3DSize = -1);

   void SetTextFormat(Double_t min, Double_t max, Double_t binWidth);

   // Renderers.
   void RnrText (const TString &txt, const TGLVector3 &pos, TGLFont::ETextAlignH_e aH, TGLFont::ETextAlignV_e aV, const TGLFont &font) const;
   void RnrTitle(const TString &title, TGLVector3 &pos, TGLFont::ETextAlignH_e aH, TGLFont::ETextAlignV_e aV) const;
   void RnrLabels() const;
   void RnrLines() const;

   void PaintAxis(TGLRnrCtx& ctx, TAxis* ax);

   ClassDef(TGLAxisPainter, 0); // GL axis painter.
};


//==============================================================================
// TGLAxisPainterBox
//==============================================================================

class TGLAxisPainterBox : public TGLAxisPainter
{
protected:
   TGLVector3          fAxisTitlePos[3];
   TAxis*                      fAxis[3];

public:
   TGLAxisPainterBox();
   virtual ~TGLAxisPainterBox();

   void SetAxis3DTitlePos(TGLRnrCtx &rnrCtx);
   void DrawAxis3D(TGLRnrCtx &rnrCtx);

   void PlotStandard(TGLRnrCtx &rnrCtx, TH1* histo, const TGLBoundingBox& bbox);

   ClassDef(TGLAxisPainterBox, 0); // Painter of GL axes for a 3D box.
};

#endif
