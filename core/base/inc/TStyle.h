// @(#)root/base:$Id: 5214d7dc2746ceb518366bc6c4f4068b652d5241 $
// Author: Rene Brun   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TStyle
#define ROOT_TStyle

#include "TNamed.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TAttMarker.h"
#include "TAttText.h"
#include "TAttAxis.h"
#include "TColor.h"

#ifndef R__LESS_INCLUDES
#include "TArrayI.h"
#endif

class TBrowser;

class TStyle : public TNamed, public TAttLine, public TAttFill, public TAttMarker, public TAttText {

private:
   TAttAxis      fXaxis;             ///< X axis attributes
   TAttAxis      fYaxis;             ///< Y axis attributes
   TAttAxis      fZaxis;             ///< Z axis attributes
   Float_t       fBarWidth;          ///< Width of bar for graphs
   Float_t       fBarOffset;         ///< Offset of bar for graphs
   Int_t         fColorModelPS;      ///< PostScript color model: 0 = RGB, 1 = CMYK
   Int_t         fDrawBorder;        ///< Flag to draw border(=1) or not (0)
   Int_t         fOptLogx;           ///< True if log scale in X
   Int_t         fOptLogy;           ///< True if log scale in y
   Int_t         fOptLogz;           ///< True if log scale in z
   Int_t         fOptDate;           ///< True if date option is selected
   Int_t         fOptStat;           ///< True if option Stat is selected
   Int_t         fOptTitle;          ///< True if option Title is selected
   Int_t         fOptFile;           ///< True if option File is selected
   Int_t         fOptFit;            ///< True if option Fit is selected
   Int_t         fShowEventStatus;   ///< Show event status panel
   Int_t         fShowEditor;        ///< Show pad editor
   Int_t         fShowToolBar;       ///< Show toolbar

   Int_t         fNumberContours;    ///< Default number of contours for 2-d plots
   TAttText      fAttDate;           ///< Canvas date attribute
   Float_t       fDateX;             ///< X position of the date in the canvas (in NDC)
   Float_t       fDateY;             ///< Y position of the date in the canvas (in NDC)
   Float_t       fEndErrorSize;      ///< Size of lines at the end of error bars
   Float_t       fErrorX;            ///< Per cent of bin width for errors along X
   Color_t       fFuncColor;         ///< Function color
   Style_t       fFuncStyle;         ///< Function style
   Width_t       fFuncWidth;         ///< Function line width
   Color_t       fGridColor;         ///< Grid line color (if 0 use axis line color)
   Style_t       fGridStyle;         ///< Grid line style
   Width_t       fGridWidth;         ///< Grid line width
   Width_t       fLegendBorderSize;  ///< Legend box border size
   Color_t       fLegendFillColor;   ///< Legend fill color
   Style_t       fLegendFont;        ///< Legend font style
   Double_t      fLegendTextSize;    ///< Legend text size. If 0 the size is computed automatically
   Int_t         fHatchesLineWidth;  ///< Hatches line width for hatch styles > 3100
   Double_t      fHatchesSpacing;    ///< Hatches spacing for hatch styles > 3100
   Color_t       fFrameFillColor;    ///< Pad frame fill color
   Color_t       fFrameLineColor;    ///< Pad frame line color
   Style_t       fFrameFillStyle;    ///< Pad frame fill style
   Style_t       fFrameLineStyle;    ///< Pad frame line style
   Width_t       fFrameLineWidth;    ///< Pad frame line width
   Width_t       fFrameBorderSize;   ///< Pad frame border size
   Int_t         fFrameBorderMode;   ///< Pad frame border mode
   Color_t       fHistFillColor;     ///< Histogram fill color
   Color_t       fHistLineColor;     ///< Histogram line color
   Style_t       fHistFillStyle;     ///< Histogram fill style
   Style_t       fHistLineStyle;     ///< Histogram line style
   Width_t       fHistLineWidth;     ///< Histogram line width
   Bool_t        fHistMinimumZero;   ///< True if default minimum is 0, false if minimum is automatic
   Double_t      fHistTopMargin;     ///< Margin between histogram's top and pad's top
   Bool_t        fCanvasPreferGL;    ///< If true, rendering in canvas is with GL
   Color_t       fCanvasColor;       ///< Canvas color
   Width_t       fCanvasBorderSize;  ///< Canvas border size
   Int_t         fCanvasBorderMode;  ///< Canvas border mode
   Int_t         fCanvasDefH;        ///< Default canvas height
   Int_t         fCanvasDefW;        ///< Default canvas width
   Int_t         fCanvasDefX;        ///< Default canvas top X position
   Int_t         fCanvasDefY;        ///< Default canvas top Y position
   Color_t       fPadColor;          ///< Pad color
   Width_t       fPadBorderSize;     ///< Pad border size
   Int_t         fPadBorderMode;     ///< Pad border mode
   Float_t       fPadBottomMargin;   ///< Pad bottom margin
   Float_t       fPadTopMargin;      ///< Pad top margin
   Float_t       fPadLeftMargin;     ///< Pad left margin
   Float_t       fPadRightMargin;    ///< Pad right margin
   Bool_t        fPadGridX;          ///< True to get the grid along X
   Bool_t        fPadGridY;          ///< True to get the grid along Y
   Int_t         fPadTickX;          ///< True to set special pad ticks along X
   Int_t         fPadTickY;          ///< True to set special pad ticks along Y
   Float_t       fPaperSizeX;        ///< PostScript paper size along X
   Float_t       fPaperSizeY;        ///< PostScript paper size along Y
   Float_t       fScreenFactor;      ///< Multiplication factor for canvas size and position
   Color_t       fStatColor;         ///< Stat fill area color
   Color_t       fStatTextColor;     ///< Stat text color
   Width_t       fStatBorderSize;    ///< Border size of Stats PaveLabel
   Style_t       fStatFont;          ///< Font style of Stats PaveLabel
   Float_t       fStatFontSize;      ///< Font size in pixels for fonts with precision type 3
   Style_t       fStatStyle;         ///< Fill area style of Stats PaveLabel
   TString       fStatFormat;        ///< Printing format for stats
   Float_t       fStatX;             ///< X position of top right corner of stat box
   Float_t       fStatY;             ///< Y position of top right corner of stat box
   Float_t       fStatW;             ///< Width of stat box
   Float_t       fStatH;             ///< Height of stat box
   Bool_t        fStripDecimals;     ///< Strip decimals in axis labels
   Int_t         fTitleAlign;        ///< Title box alignment
   Color_t       fTitleColor;        ///< Title fill area color
   Color_t       fTitleTextColor;    ///< Title text color
   Width_t       fTitleBorderSize;   ///< Border size of Title PavelLabel
   Style_t       fTitleFont;         ///< Font style of Title PaveLabel
   Float_t       fTitleFontSize;     ///< Font size in pixels for fonts with precision type 3
   Style_t       fTitleStyle;        ///< Fill area style of title PaveLabel
   Float_t       fTitleX;            ///< X position of top left corner of title box
   Float_t       fTitleY;            ///< Y position of top left corner of title box
   Float_t       fTitleW;            ///< Width of title box
   Float_t       fTitleH;            ///< Height of title box
   Float_t       fLegoInnerR;        ///< Inner radius for cylindrical legos
   TString       fLineStyle[30];     ///< String describing line style i (for postScript)
   TString       fHeaderPS;          ///< User defined additional Postscript header
   TString       fTitlePS;           ///< User defined Postscript file title
   TString       fFitFormat;         ///< Printing format for fit parameters
   TString       fPaintTextFormat;   ///< Printing format for TH2::PaintText
   Float_t       fLineScalePS;       ///< Line scale factor when drawing lines on Postscript
   Int_t         fJoinLinePS;        ///< Determines the appearance of joining lines on PostScript, PDF and SVG
   Int_t         fCapLinePS;         ///< Determines the appearance of line caps on PostScript, PDF and SVG
   Double_t      fTimeOffset;        ///< Time offset to the beginning of an axis
   Bool_t        fIsReading;         ///<! Set to FALSE when userclass::UseCurrentStyle is called by the style manager
   Float_t       fImageScaling;      ///< Image scaling to produce high definition bitmap images

public:
   enum EPaperSize { kA4, kUSLetter };

   TStyle();
   TStyle(const char *name, const char *title);
   TStyle(const TStyle &style);
   TStyle& operator=(const TStyle& style);
   virtual          ~TStyle();
   inline Int_t     AxisChoice(Option_t *axis) const {
      // Return axis number (1 for X, 2 for Y, 3 for Z)
      UChar_t a = *axis;
      a -= (a >= 'x') ? 'x' : 'X'; // toupper and a-='X'; intentional underflow
      return (a > 2) ? 0 : (Int_t)(a+1);
   };
   virtual void     Browse(TBrowser *b);
   static  void     BuildStyles();
   virtual void     Copy(TObject &style) const;
   virtual void     cd();

   virtual Int_t    DistancetoPrimitive(Int_t px, Int_t py);
   Int_t            GetNdivisions(Option_t *axis="X") const;
   TAttText        *GetAttDate() {return &fAttDate;}
   Color_t          GetAxisColor(Option_t *axis="X") const;
   Color_t          GetLabelColor(Option_t *axis="X") const;
   Style_t          GetLabelFont(Option_t *axis="X") const;
   Float_t          GetLabelOffset(Option_t *axis="X") const;
   Float_t          GetLabelSize(Option_t *axis="X") const;
   Color_t          GetTitleColor(Option_t *axis="X") const;  //return axis title color of pad title color
   Style_t          GetTitleFont(Option_t *axis="X") const;   //return axis title font of pad title font
   Float_t          GetTitleOffset(Option_t *axis="X") const; //return axis title offset
   Float_t          GetTitleSize(Option_t *axis="X") const;   //return axis title size
   Float_t          GetTickLength(Option_t *axis="X") const;

   Float_t          GetBarOffset() const {return fBarOffset;}
   Float_t          GetBarWidth() const {return fBarWidth;}
   Int_t            GetDrawBorder() const {return fDrawBorder;}
   Float_t          GetEndErrorSize() const {return fEndErrorSize;}
   Float_t          GetErrorX() const {return fErrorX;}
   Bool_t           GetCanvasPreferGL() const {return fCanvasPreferGL;}
   Color_t          GetCanvasColor() const {return fCanvasColor;}
   Width_t          GetCanvasBorderSize() const {return fCanvasBorderSize;}
   Int_t            GetCanvasBorderMode() const {return fCanvasBorderMode;}
   Int_t            GetCanvasDefH() const      {return fCanvasDefH;}
   Int_t            GetCanvasDefW() const      {return fCanvasDefW;}
   Int_t            GetCanvasDefX() const      {return fCanvasDefX;}
   Int_t            GetCanvasDefY() const      {return fCanvasDefY;}
   Int_t            GetColorPalette(Int_t i) const;
   Int_t            GetColorModelPS() const    {return fColorModelPS;}
   Float_t          GetDateX()  const          {return fDateX;}
   Float_t          GetDateY() const           {return fDateY;}
   const char      *GetFitFormat()       const {return fFitFormat.Data();}
   Int_t            GetHatchesLineWidth() const {return fHatchesLineWidth;}
   Double_t         GetHatchesSpacing() const  {return fHatchesSpacing;}
   Width_t          GetLegendBorderSize() const   {return fLegendBorderSize;}
   Color_t          GetLegendFillColor() const {return fLegendFillColor;}
   Style_t          GetLegendFont() const {return fLegendFont;}
   Double_t         GetLegendTextSize() const {return fLegendTextSize;}
   Int_t            GetNumberOfColors() const;
   Color_t          GetPadColor() const        {return fPadColor;}
   Width_t          GetPadBorderSize() const   {return fPadBorderSize;}
   Int_t            GetPadBorderMode() const   {return fPadBorderMode;}
   Float_t          GetPadBottomMargin() const {return fPadBottomMargin;}
   Float_t          GetPadTopMargin() const    {return fPadTopMargin;}
   Float_t          GetPadLeftMargin() const   {return fPadLeftMargin;}
   Float_t          GetPadRightMargin() const  {return fPadRightMargin;}
   Bool_t           GetPadGridX() const        {return fPadGridX;}
   Bool_t           GetPadGridY() const        {return fPadGridY;}
   Int_t            GetPadTickX() const        {return fPadTickX;}
   Int_t            GetPadTickY() const        {return fPadTickY;}
   Color_t          GetFuncColor() const       {return fFuncColor;}
   Style_t          GetFuncStyle() const       {return fFuncStyle;}
   Width_t          GetFuncWidth() const       {return fFuncWidth;}
   Color_t          GetGridColor() const       {return fGridColor;}
   Style_t          GetGridStyle() const       {return fGridStyle;}
   Width_t          GetGridWidth() const       {return fGridWidth;}
   Color_t          GetFrameFillColor()  const {return fFrameFillColor;}
   Color_t          GetFrameLineColor()  const {return fFrameLineColor;}
   Style_t          GetFrameFillStyle()  const {return fFrameFillStyle;}
   Style_t          GetFrameLineStyle()  const {return fFrameLineStyle;}
   Width_t          GetFrameLineWidth()  const {return fFrameLineWidth;}
   Width_t          GetFrameBorderSize() const {return fFrameBorderSize;}
   Int_t            GetFrameBorderMode() const {return fFrameBorderMode;}
   Color_t          GetHistFillColor()   const {return fHistFillColor;}
   Color_t          GetHistLineColor()   const {return fHistLineColor;}
   Style_t          GetHistFillStyle()   const {return fHistFillStyle;}
   Style_t          GetHistLineStyle()   const {return fHistLineStyle;}
   Width_t          GetHistLineWidth()   const {return fHistLineWidth;}
   Bool_t           GetHistMinimumZero() const {return fHistMinimumZero;}
   Double_t         GetHistTopMargin()   const {return fHistTopMargin;}
   Float_t          GetImageScaling()    const {return fImageScaling;}
   Float_t          GetLegoInnerR() const {return fLegoInnerR;}
   Int_t            GetNumberContours() const {return fNumberContours;}
   Int_t            GetOptDate() const {return fOptDate;}
   Int_t            GetOptFile() const {return fOptFile;}
   Int_t            GetOptFit() const {return fOptFit;}
   Int_t            GetOptStat() const {return fOptStat;}
   Int_t            GetOptTitle() const {return fOptTitle;}
   Int_t            GetOptLogx() const {return fOptLogx;}
   Int_t            GetOptLogy() const {return fOptLogy;}
   Int_t            GetOptLogz() const {return fOptLogz;}
   const char      *GetPaintTextFormat() const {return fPaintTextFormat.Data();}
   void             GetPaperSize(Float_t &xsize, Float_t &ysize) const;
   Int_t            GetShowEventStatus() const {return fShowEventStatus;}
   Int_t            GetShowEditor() const {return fShowEditor;}
   Int_t            GetShowToolBar() const {return fShowToolBar;}

   Float_t          GetScreenFactor() const {return fScreenFactor;}
   Color_t          GetStatColor() const {return fStatColor;}
   Color_t          GetStatTextColor() const {return fStatTextColor;}
   Width_t          GetStatBorderSize() const {return fStatBorderSize;}
   Style_t          GetStatFont() const  {return fStatFont;}
   Float_t          GetStatFontSize() const  {return fStatFontSize;}
   Style_t          GetStatStyle() const  {return fStatStyle;}
   const char      *GetStatFormat() const {return fStatFormat.Data();}
   Float_t          GetStatX() const     {return fStatX;}
   Float_t          GetStatY() const     {return fStatY;}
   Float_t          GetStatW() const     {return fStatW;}
   Float_t          GetStatH() const     {return fStatH;}
   Int_t            GetStripDecimals() const {return fStripDecimals;}
   Double_t         GetTimeOffset() const {return fTimeOffset;} //return axis time offset
   Int_t            GetTitleAlign() {return fTitleAlign;} // return the histogram title TPaveLabel alignment
   Color_t          GetTitleFillColor() const {return fTitleColor;}  //return histogram title fill area color
   Color_t          GetTitleTextColor() const {return fTitleTextColor;}  //return histogram title text color
   Style_t          GetTitleStyle() const  {return fTitleStyle;}
   Float_t          GetTitleFontSize() const  {return fTitleFontSize;} //return histogram title font size
   Width_t          GetTitleBorderSize() const {return fTitleBorderSize;} //return border size of histogram title TPaveLabel
   Float_t          GetTitleXOffset() const {return GetTitleOffset("X");} //return X axis title offset
   Float_t          GetTitleXSize() const   {return GetTitleSize("X");}   //return X axis title size
   Float_t          GetTitleYOffset() const {return GetTitleOffset("Y");} //return Y axis title offset
   Float_t          GetTitleYSize() const   {return GetTitleSize("Y");}   //return Y axis title size
   Float_t          GetTitleX() const     {return fTitleX;}  //return left X position of histogram title TPavelabel
   Float_t          GetTitleY() const     {return fTitleY;}  //return left bottom position of histogram title TPavelabel
   Float_t          GetTitleW() const     {return fTitleW;}  //return width of histogram title TPaveLabel
   Float_t          GetTitleH() const     {return fTitleH;}  //return height of histogram title TPavelabel
   const char      *GetHeaderPS() const {return fHeaderPS.Data();}
   const char      *GetTitlePS()  const {return fTitlePS.Data();}
   const char      *GetLineStyleString(Int_t i=1) const;
   Int_t            GetJoinLinePS() const {return fJoinLinePS;} ///< Returns the line join method used for PostScript, PDF and SVG output. See `TPostScript::SetLineJoin` for details.
   Int_t            GetCapLinePS()  const {return fCapLinePS;}  ///< Returns the line cap method used for PostScript, PDF and SVG output. See `TPostScript::SetLineCap` for details.
   Float_t          GetLineScalePS() const {return fLineScalePS;}

   Bool_t           IsReading() const {return fIsReading;}
   virtual void     Paint(Option_t *option="");
   virtual void     Reset(Option_t *option="");

   void             SetColorModelPS(Int_t c=0);
   void             SetFitFormat(const char *format="5.4g") {fFitFormat = format;}
   void             SetHeaderPS(const char *header);
   void             SetHatchesLineWidth(Int_t l) {fHatchesLineWidth = l;}
   void             SetHatchesSpacing(Double_t h) {fHatchesSpacing = TMath::Max(0.1,h);}
   void             SetTitlePS(const char *pstitle);
   void             SetJoinLinePS(Int_t joinline=0) {fJoinLinePS=joinline;} ///< Set the line join method used for PostScript, PDF and SVG output. See `TPostScript::SetLineJoin` for details.
   void             SetCapLinePS(Int_t capline=0) {fCapLinePS=capline;}     ///< Set the line cap method used for PostScript, PDF and SVG output. See `TPostScript::SetLineCap` for details.
   void             SetLineScalePS(Float_t scale=3) {fLineScalePS=scale;}
   void             SetLineStyleString(Int_t i, const char *text);
   void             SetNdivisions(Int_t n=510, Option_t *axis="X");
   void             SetAxisColor(Color_t color=1, Option_t *axis="X");
   void             SetLabelColor(Color_t color=1, Option_t *axis="X");
   void             SetLabelFont(Style_t font=62, Option_t *axis="X");
   void             SetLabelOffset(Float_t offset=0.005, Option_t *axis="X");
   void             SetLabelSize(Float_t size=0.04, Option_t *axis="X");
   void             SetLegoInnerR(Float_t rad=0.5) {fLegoInnerR = rad;}
   void             SetScreenFactor(Float_t factor=1) {fScreenFactor = factor;}
   void             SetTickLength(Float_t length=0.03, Option_t *axis="X");
   void             SetTitleColor(Color_t color=1, Option_t *axis="X"); //set axis title color or pad title color
   void             SetTitleFont(Style_t font=62, Option_t *axis="X"); //set axis title font or pad title font
   void             SetTitleOffset(Float_t offset=1, Option_t *axis="X"); //set axis title offset
   void             SetTitleSize(Float_t size=0.02, Option_t *axis="X");  //set axis title size or pad title size
   void             SetNumberContours(Int_t number=20);
   void             SetOptDate(Int_t datefl=1);
   void             SetOptFile(Int_t file=1) {fOptFile = file;}
   void             SetOptFit(Int_t fit=1);
   void             SetOptLogx(Int_t logx=1) {fOptLogx = logx;}
   void             SetOptLogy(Int_t logy=1) {fOptLogy = logy;}
   void             SetOptLogz(Int_t logz=1) {fOptLogz = logz;}
   void             SetOptStat(Int_t stat=1);
   void             SetOptStat(Option_t *stat);
   void             SetOptTitle(Int_t tit=1) {fOptTitle = tit;}
   void             SetBarOffset(Float_t baroff=0.5) {fBarOffset = baroff;}
   void             SetBarWidth(Float_t barwidth=0.5) {fBarWidth = barwidth;}
   void             SetDateX(Float_t x=0.01) {fDateX = x;}
   void             SetDateY(Float_t y=0.01) {fDateY = y;}
   void             SetEndErrorSize(Float_t np=2);
   void             SetErrorX(Float_t errorx=0.5) {fErrorX = errorx;}
   void             SetCanvasPreferGL(Bool_t prefer = kTRUE) {fCanvasPreferGL=prefer;}
   void             SetDrawBorder(Int_t drawborder=1) {fDrawBorder = drawborder;}
   void             SetCanvasColor(Color_t color=19) {fCanvasColor = color;}
   void             SetCanvasBorderSize(Width_t size=1) {fCanvasBorderSize = size;}
   void             SetCanvasBorderMode(Int_t mode=1) {fCanvasBorderMode = mode;}
   void             SetCanvasDefH(Int_t h=500) {fCanvasDefH = h;}
   void             SetCanvasDefW(Int_t w=700) {fCanvasDefW = w;}
   void             SetCanvasDefX(Int_t topx=10) {fCanvasDefX = topx;}
   void             SetCanvasDefY(Int_t topy=10) {fCanvasDefY = topy;}
   void             SetLegendBorderSize(Width_t size=4) {fLegendBorderSize = size;}
   void             SetLegendFillColor(Color_t color=0) {fLegendFillColor = color;}
   void             SetLegendFont(Style_t font=62) {fLegendFont = font;}
   void             SetLegendTextSize(Double_t size=0.) {fLegendTextSize = size;}
   void             SetPadColor(Color_t color=19) {fPadColor = color;}
   void             SetPadBorderSize(Width_t size=1) {fPadBorderSize = size;}
   void             SetPadBorderMode(Int_t mode=1) {fPadBorderMode = mode;}
   void             SetPadBottomMargin(Float_t margin=0.1) {fPadBottomMargin=margin;}
   void             SetPadTopMargin(Float_t margin=0.1)    {fPadTopMargin=margin;}
   void             SetPadLeftMargin(Float_t margin=0.1)   {fPadLeftMargin=margin;}
   void             SetPadRightMargin(Float_t margin=0.1)  {fPadRightMargin=margin;}
   void             SetPadGridX(Bool_t gridx) {fPadGridX = gridx;}
   void             SetPadGridY(Bool_t gridy) {fPadGridY = gridy;}
   void             SetPadTickX(Int_t tickx)  {fPadTickX = tickx;}
   void             SetPadTickY(Int_t ticky)  {fPadTickY = ticky;}
   void             SetFuncStyle(Style_t style=1) {fFuncStyle = style;}
   void             SetFuncColor(Color_t color=1) {fFuncColor = color;}
   void             SetFuncWidth(Width_t width=4) {fFuncWidth = width;}
   void             SetGridStyle(Style_t style=3) {fGridStyle = style;}
   void             SetGridColor(Color_t color=0) {fGridColor = color;}
   void             SetGridWidth(Width_t width=1) {fGridWidth = width;}
   void             SetFrameFillColor(Color_t color=1) {fFrameFillColor = color;}
   void             SetFrameLineColor(Color_t color=1) {fFrameLineColor = color;}
   void             SetFrameFillStyle(Style_t styl=0)  {fFrameFillStyle = styl;}
   void             SetFrameLineStyle(Style_t styl=0)  {fFrameLineStyle = styl;}
   void             SetFrameLineWidth(Width_t width=1) {fFrameLineWidth = width;}
   void             SetFrameBorderSize(Width_t size=1) {fFrameBorderSize = size;}
   void             SetFrameBorderMode(Int_t mode=1) {fFrameBorderMode = mode;}
   void             SetHistFillColor(Color_t color=1) {fHistFillColor = color;}
   void             SetHistLineColor(Color_t color=1) {fHistLineColor = color;}
   void             SetHistFillStyle(Style_t styl=0)  {fHistFillStyle = styl;}
   void             SetHistLineStyle(Style_t styl=0)  {fHistLineStyle = styl;}
   void             SetHistLineWidth(Width_t width=1) {fHistLineWidth = width;}
   void             SetHistMinimumZero(Bool_t zero=kTRUE);
   void             SetHistTopMargin(Double_t hmax=0.05) {fHistTopMargin = hmax;}
   void             SetPaintTextFormat(const char *format="g") {fPaintTextFormat = format;}
   void             SetImageScaling(Float_t s){fImageScaling = s;}
   void             SetPaperSize(EPaperSize size);
   void             SetPaperSize(Float_t xsize=20, Float_t ysize=26);
   void             SetStatColor(Color_t color=19) {fStatColor=color;}
   void             SetStatTextColor(Color_t color=1) {fStatTextColor=color;}
   void             SetStatStyle(Style_t style=1001) {fStatStyle=style;}
   void             SetStatBorderSize(Width_t size=2) {fStatBorderSize=size;}
   void             SetStatFont(Style_t font=62) {fStatFont=font;}
   void             SetStatFontSize(Float_t size=0)  {fStatFontSize=size;}
   void             SetStatFormat(const char *format="6.4g") {fStatFormat = format;}
   void             SetStatX(Float_t x=0)    {fStatX=x;}
   void             SetStatY(Float_t y=0)    {fStatY=y;}
   void             SetStatW(Float_t w=0.19) {fStatW=w;}
   void             SetStatH(Float_t h=0.1)  {fStatH=h;}
   void             SetStripDecimals(Bool_t strip=kTRUE);
   void             SetTimeOffset(Double_t toffset);
   void             SetTitleAlign(Int_t a=13) {fTitleAlign=a;}
   void             SetTitleFillColor(Color_t color=1)   {fTitleColor=color;}
   void             SetTitleTextColor(Color_t color=1)   {fTitleTextColor=color;}
   void             SetTitleStyle(Style_t style=1001)  {fTitleStyle=style;}
   void             SetTitleFontSize(Float_t size=0)   {fTitleFontSize=size;}
   void             SetTitleBorderSize(Width_t size=2) {fTitleBorderSize=size;}
   void             SetTitleXOffset(Float_t offset=1)  {SetTitleOffset(offset,"X");}
   void             SetTitleXSize(Float_t size=0.02)   {SetTitleSize(size,"X");}
   void             SetTitleYOffset(Float_t offset=1)  {SetTitleOffset(offset,"Y");}
   void             SetTitleYSize(Float_t size=0.02)   {SetTitleSize(size,"Y");}
   void             SetTitleX(Float_t x=0)     {fTitleX=x;}
   void             SetTitleY(Float_t y=0.985) {fTitleY=y;}
   void             SetTitleW(Float_t w=0)     {fTitleW=w;}
   void             SetTitleH(Float_t h=0)     {fTitleH=h;}
   void             ToggleEventStatus() { fShowEventStatus = fShowEventStatus ? 0 : 1; }
   void             ToggleEditor() { fShowEditor = fShowEditor ? 0 : 1; }
   void             ToggleToolBar() { fShowToolBar = fShowToolBar ? 0 : 1; }
   void             SetIsReading(Bool_t reading=kTRUE);
   void             SetPalette(Int_t ncolors=kBird, Int_t *colors=0, Float_t alpha=1.);
   void             SavePrimitive(std::ostream &out, Option_t * = "");
   void             SaveSource(const char *filename, Option_t *option=0);

   ClassDef(TStyle, 19);  //A collection of all graphics attributes
};


R__EXTERN TStyle  *gStyle;

#endif
