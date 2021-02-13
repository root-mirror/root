// @(#)root/gl:$Id$
// Author:  Timur Pocheptsov  14/06/2006

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLLegoPainter
#define ROOT_TGLLegoPainter

#include "TGLPlotPainter.h"
#include "TGLQuadric.h"
#include "TGLUtil.h"
#include "TString.h"

#include <vector>
#include <utility>

class TGLPlotCamera;
class TAxis;
class TH1;

/*
   TGLLegoPainter. The concrete implementation of abstract TGLPlotPainter.
*/

class TGLLegoPainter : public TGLPlotPainter {
private:

   enum ELegoType {
      kColorSimple,
      kColorLevel,
      kCylindricBars
   };
   //Bars, cylinders or textured bars.
   mutable ELegoType     fLegoType;
   Double_t              fMinZ;
   Rgl::Range_t          fMinMaxVal;//For texture coordinates generation.

   std::vector<Rgl::Range_t>  fXEdges;
   std::vector<Rgl::Range_t>  fYEdges;

   typedef std::pair<Double_t, Double_t> CosSin_t;
   std::vector<CosSin_t> fCosSinTableX;
   std::vector<CosSin_t> fCosSinTableY;
   TString               fBinInfo;
   mutable TGLQuadric    fQuadric;
   Bool_t                fDrawErrors;

   mutable TGLLevelPalette       fPalette;
   mutable std::vector<Double_t> fColorLevels;

   TGLLegoPainter(const TGLLegoPainter &);
   TGLLegoPainter &operator = (const TGLLegoPainter &);

public:
   TGLLegoPainter(TH1 *hist, TGLPlotCamera *camera, TGLPlotCoordinates *coord);

   //TGLPlotPainter's final-overriders
   char        *GetPlotInfo(Int_t px, Int_t py);
   Bool_t       InitGeometry();
   void         StartPan(Int_t px, Int_t py);
   void         Pan(Int_t px, Int_t py);
   void         AddOption(const TString &stringOption);
   void         ProcessEvent(Int_t event, Int_t px, Int_t py);

private:
   //Auxilary functions.
   Bool_t       InitGeometryCartesian();
   Bool_t       InitGeometryPolar();
   Bool_t       InitGeometryCylindrical();
   Bool_t       InitGeometrySpherical();
   //Overriders
   void         InitGL()const;
   void         DeInitGL()const;

   void         DrawPlot()const;

   void         DrawLegoCartesian()const;
   void         DrawLegoPolar()const;
   void         DrawLegoCylindrical()const;
   void         DrawLegoSpherical()const;

   void         SetLegoColor()const;

   void         DrawSectionXOZ()const;
   void         DrawSectionYOZ()const;
   void         DrawSectionXOY()const;

   Bool_t       ClampZ(Double_t &zVal)const;
   Bool_t       PreparePalette()const;

   void         DrawPalette()const;
   void         DrawPaletteAxis()const;

   ClassDef(TGLLegoPainter, 0)//Lego painter
};

#endif
