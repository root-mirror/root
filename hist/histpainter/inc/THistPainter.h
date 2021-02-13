// @(#)root/histpainter:$Id$
// Author: Rene Brun, Olivier Couet

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_THistPainter
#define ROOT_THistPainter


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THistPainter                                                         //
//                                                                      //
// helper class to draw histograms                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TVirtualHistPainter.h"
#include "TString.h"

#include <vector>
#include <utility>


class TH1;
class TAxis;
class TCutG;
class TGaxis;
class TPainter3dAlgorithms;
class TGraph2DPainter;
class TPie;
class TF3;

const Int_t kMaxCuts = 16;

struct THistRenderingRegion
{
   std::pair<Int_t, Int_t> fPixelRange;
   std::pair<Int_t, Int_t> fBinRange;
};


class THistPainter : public TVirtualHistPainter {

protected:
   TH1                  *fH;                 //pointer to histogram to paint
   TAxis                *fXaxis;             //pointer to X axis
   TAxis                *fYaxis;             //pointer to Y axis
   TAxis                *fZaxis;             //pointer to Z axis
   TList                *fFunctions;         //pointer to histogram list of functions
   TPainter3dAlgorithms *fLego;              //pointer to a TPainter3dAlgorithms object
   TGraph2DPainter      *fGraph2DPainter;    //pointer to a TGraph2DPainter object
   TPie                 *fPie;               //pointer to a TPie in case of option PIE
   Double_t             *fXbuf;              //X buffer coordinates
   Double_t             *fYbuf;              //Y buffer coordinates
   Int_t                 fNcuts;             //Number of graphical cuts
   Int_t                 fCutsOpt[kMaxCuts]; //sign of each cut
   TCutG                *fCuts[kMaxCuts];    //Pointers to graphical cuts
   TList                *fStack;             //Pointer to stack of histograms (if any)
   Int_t                 fShowProjection;    //True if a projection must be drawn
   TString               fShowOption;        //Option to draw the projection
   Int_t                 fXHighlightBin;     //X highlight bin
   Int_t                 fYHighlightBin;     //Y highlight bin
   TF3                  *fCurrentF3;         //current TF3 function

private:
   mutable TString fObjectInfo;

public:
   THistPainter();
   virtual ~THistPainter();
   virtual void       DefineColorLevels(Int_t ndivz);
   virtual Int_t      DistancetoPrimitive(Int_t px, Int_t py);
   virtual void       DrawPanel();
   virtual void       ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual TList     *GetContourList(Double_t contour) const;
   virtual char      *GetObjectInfo(Int_t px, Int_t py) const;
   virtual TList     *GetStack() const {return fStack;}
   virtual Int_t      GetXHighlightBin() const { return fXHighlightBin; }
   virtual Int_t      GetYHighlightBin() const { return fYHighlightBin; }
   virtual void       HighlightBin(Int_t px, Int_t py);
   virtual Bool_t     IsInside(Int_t x, Int_t y);
   virtual Bool_t     IsInside(Double_t x, Double_t y);
   virtual Int_t      MakeChopt(Option_t *option);
   virtual Int_t      MakeCuts(char *cutsopt);
   virtual void       Paint(Option_t *option="");
   virtual void       PaintArrows(Option_t *option);
   virtual void       PaintAxis(Bool_t drawGridOnly=kFALSE);
   virtual void       PaintBar(Option_t *option);
   virtual void       PaintBarH(Option_t *option);
   virtual void       PaintBoxes(Option_t *option);
   virtual void       PaintCandlePlot(Option_t *option);
   virtual void       PaintColorLevels(Option_t *option);
   virtual void       PaintColorLevelsFast(Option_t *option);
   virtual std::vector<THistRenderingRegion> ComputeRenderingRegions(TAxis *pAxis, Int_t nPixels, bool isLog);

   virtual void       PaintTH2PolyBins(Option_t *option);
   virtual void       PaintTH2PolyColorLevels(Option_t *option);
   virtual void       PaintTH2PolyScatterPlot(Option_t *option);
   virtual void       PaintTH2PolyText(Option_t *option);
   virtual void       PaintContour(Option_t *option);
   virtual Int_t      PaintContourLine(Double_t elev1, Int_t icont1, Double_t x1, Double_t y1,
                                       Double_t elev2, Int_t icont2, Double_t x2, Double_t y2,
                                       Double_t *xarr, Double_t *yarr, Int_t *itarr, Double_t *levels);
   virtual void       PaintErrors(Option_t *option);
   virtual void       Paint2DErrors(Option_t *option);
   virtual void       PaintFrame();
   virtual void       PaintFunction(Option_t *option);
   virtual void       PaintHighlightBin(Option_t *option="");
   virtual void       PaintHist(Option_t *option);
   virtual void       PaintH3(Option_t *option="");
   virtual void       PaintH3Box(Int_t iopt);
   virtual void       PaintH3BoxRaster();
   virtual void       PaintH3Iso();
   virtual Int_t      PaintInit();
   virtual Int_t      PaintInitH();
   virtual void       PaintLego(Option_t *option);
   virtual void       PaintLegoAxis(TGaxis *axis, Double_t ang);
   virtual void       PaintPalette();
   virtual void       PaintScatterPlot(Option_t *option);
   virtual void       PaintStat(Int_t dostat, TF1 *fit);
   virtual void       PaintStat2(Int_t dostat, TF1 *fit);
   virtual void       PaintStat3(Int_t dostat, TF1 *fit);
   virtual void       PaintSurface(Option_t *option);
   virtual void       PaintTriangles(Option_t *option);
   virtual void       PaintTable(Option_t *option);
   virtual void       PaintText(Option_t *option);
   virtual void       PaintTitle();
   virtual void       PaintTF3();
   virtual void       ProcessMessage(const char *mess, const TObject *obj);
   static  Int_t      ProjectAitoff2xy(Double_t l, Double_t b, Double_t &Al, Double_t &Ab);
   static  Int_t      ProjectMercator2xy(Double_t l, Double_t b, Double_t &Al, Double_t &Ab);
   static  Int_t      ProjectSinusoidal2xy(Double_t l, Double_t b, Double_t &Al, Double_t &Ab);
   static  Int_t      ProjectParabolic2xy(Double_t l, Double_t b, Double_t &Al, Double_t &Ab);
   virtual void       RecalculateRange();
   virtual void       RecursiveRemove(TObject *) {;}
   virtual void       SetHighlight();
   virtual void       SetHistogram(TH1 *h);
   virtual void       SetStack(TList *stack) {fStack = stack;}
   virtual void       SetShowProjection(const char *option,Int_t nbins);
   virtual void       ShowProjectionX(Int_t px, Int_t py);
   virtual void       ShowProjectionY(Int_t px, Int_t py);
   virtual void       ShowProjection3(Int_t px, Int_t py);
   virtual Int_t      TableInit();

   static const char *GetBestFormat(Double_t v, Double_t e, const char *f);
   static void        PaintSpecialObjects(const TObject *obj, Option_t *option);

   ClassDef(THistPainter,0)  //Helper class to draw histograms
};

#endif
