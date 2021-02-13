// @(#)root/base:$Id$
// Author: Rene Brun   05/12/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TVirtualPad
#define ROOT_TVirtualPad


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualPad                                                          //
//                                                                      //
// Abstract base class for Pads and Canvases                            //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TObject.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TAttPad.h"
#include "TQObject.h"

#include "GuiTypes.h"
#include "TString.h"
#include "Buttons.h"

// forward declarations
class TAxis;
class TObject;
class TObjLink;
class TView;
class TCanvas;
class TCanvasImp;
class TH1F;
class TFrame;
class TLegend;
class TBox;
class TVirtualViewer3D;
class TVirtualPadPainter;

class TVirtualPad : public TObject, public TAttLine, public TAttFill,
                    public TAttPad, public TQObject {

protected:
   Bool_t         fResizing;         //!true when resizing the pad

   void  *GetSender() override { return this; }  //used to set gTQSender

public:
   TVirtualPad();
   TVirtualPad(const char *name, const char *title, Double_t xlow,
               Double_t ylow, Double_t xup, Double_t yup,
               Color_t color=19, Short_t bordersize=4, Short_t bordermode=1);
   virtual ~TVirtualPad();
   virtual void     AbsCoordinates(Bool_t set) = 0;
   virtual Double_t AbsPixeltoX(Int_t px) = 0;
   virtual Double_t AbsPixeltoY(Int_t py) = 0;
   virtual void     AddExec(const char *name, const char *command) = 0;
   virtual TLegend *BuildLegend(Double_t x1=0.3, Double_t y1=0.21, Double_t x2=0.3, Double_t y2=0.21, const char *title="", Option_t *option = "") = 0;
   virtual TVirtualPad* cd(Int_t subpadnumber=0) = 0;
           void     Clear(Option_t *option="") override = 0;
   virtual Int_t    Clip(Double_t *x, Double_t *y, Double_t xclipl, Double_t yclipb, Double_t xclipr, Double_t yclipt) = 0;
   virtual void     Close(Option_t *option="") = 0;
   virtual void     CopyPixmap() = 0;
   virtual void     CopyPixmaps() = 0;
   virtual void     DeleteExec(const char *name) = 0;
   virtual void     Divide(Int_t nx=1, Int_t ny=1, Float_t xmargin=0.01, Float_t ymargin=0.01, Int_t color=0) = 0;
           void     Draw(Option_t *option="") override = 0;
   virtual void     DrawClassObject(const TObject *obj, Option_t *option="") = 0;
   virtual TH1F    *DrawFrame(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, const char *title="") = 0;
   virtual void     ExecuteEventAxis(Int_t event, Int_t px, Int_t py, TAxis *axis) = 0;
   virtual void     UnZoomed() { Emit("UnZoomed()"); } // *SIGNAL*
   virtual Short_t  GetBorderMode() const = 0;
   virtual Short_t  GetBorderSize() const = 0;
   virtual Int_t    GetCanvasID() const = 0;
   virtual TCanvasImp *GetCanvasImp() const = 0;
   virtual TCanvas  *GetCanvas() const = 0;
   virtual TVirtualPad *GetVirtCanvas() const = 0;
   virtual Int_t    GetEvent() const  = 0;
   virtual Int_t    GetEventX() const = 0;
   virtual Int_t    GetEventY() const = 0;
   virtual TFrame   *GetFrame() = 0;
   virtual Color_t  GetHighLightColor() const = 0;
   virtual Int_t    GetNumber() const = 0;
   virtual void     GetRange(Double_t &x1, Double_t &y1, Double_t &x2, Double_t &y2) = 0;
   virtual void     GetRangeAxis(Double_t &xmin, Double_t &ymin, Double_t &xmax, Double_t &ymax) = 0;
   virtual void     GetPadPar(Double_t &xlow, Double_t &ylow, Double_t &xup, Double_t &yup) = 0;
   virtual Double_t GetXlowNDC() const = 0;
   virtual Double_t GetYlowNDC() const = 0;
   virtual Double_t GetWNDC() const = 0;
   virtual Double_t GetHNDC() const = 0;
   virtual UInt_t   GetWw() const = 0;
   virtual UInt_t   GetWh() const = 0;
   virtual Double_t GetAbsXlowNDC() const = 0;
   virtual Double_t GetAbsYlowNDC() const = 0;
   virtual Double_t GetAbsWNDC() const = 0;
   virtual Double_t GetAbsHNDC() const = 0;
   virtual Double_t GetAspectRatio() const = 0;
   virtual Double_t GetPhi() const = 0;
   virtual Double_t GetTheta() const = 0;
   virtual Double_t GetUxmin() const = 0;
   virtual Double_t GetUymin() const = 0;
   virtual Double_t GetUxmax() const = 0;
   virtual Double_t GetUymax() const = 0;
   virtual Bool_t   GetGridx() const = 0;
   virtual Bool_t   GetGridy() const = 0;
   virtual Int_t    GetTickx() const = 0;
   virtual Int_t    GetTicky() const = 0;
   virtual Double_t GetX1() const = 0;
   virtual Double_t GetX2() const = 0;
   virtual Double_t GetY1() const = 0;
   virtual Double_t GetY2() const = 0;
   virtual TList    *GetListOfPrimitives() const = 0;
   virtual TList    *GetListOfExecs() const = 0;
   virtual TObject  *GetPrimitive(const char *name) const = 0;
   virtual TObject  *GetSelected() const = 0;
   virtual TVirtualPad  *GetPad(Int_t subpadnumber) const = 0;
   virtual TObject  *GetPadPointer() const = 0;
   virtual TVirtualPad  *GetPadSave() const = 0;
   virtual TVirtualPad  *GetSelectedPad() const = 0;
   virtual TView    *GetView() const = 0;
   virtual Int_t    GetLogx() const = 0;
   virtual Int_t    GetLogy() const = 0;
   virtual Int_t    GetLogz() const = 0;
   virtual TVirtualPad  *GetMother() const = 0;
           const char *GetName() const override = 0;
           const char *GetTitle() const override = 0;
   virtual Int_t    GetPadPaint() const = 0;
   virtual Int_t    GetPixmapID() const = 0;
   virtual TObject *GetView3D() const = 0;
   virtual Bool_t   HasCrosshair() const = 0;
   virtual void     HighLight(Color_t col=kRed, Bool_t set=kTRUE) = 0;
   virtual Bool_t   HasFixedAspectRatio() const = 0;
   virtual Bool_t   IsBatch() const = 0;
   Bool_t           IsBeingResized() const { return fResizing; }
   virtual Bool_t   IsEditable() const = 0;
   virtual Bool_t   IsModified() const = 0;
   virtual Bool_t   IsRetained() const = 0;
   virtual Bool_t   IsVertical() const = 0;
           void     ls(Option_t *option="") const override = 0;
   virtual void     Modified(Bool_t flag=1) = 0;
   virtual Bool_t   OpaqueMoving() const = 0;
   virtual Bool_t   OpaqueResizing() const = 0;
   virtual Double_t PadtoX(Double_t x) const = 0;
   virtual Double_t PadtoY(Double_t y) const = 0;
           void     Paint(Option_t *option="") override = 0;
   virtual void     PaintBorderPS(Double_t xl,Double_t yl,Double_t xt,Double_t yt,Int_t bmode,Int_t bsize,Int_t dark,Int_t light) = 0;
   virtual void     PaintBox(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="") = 0;
   virtual void     PaintFillArea(Int_t n, Float_t *x, Float_t *y, Option_t *option="") = 0;
   virtual void     PaintFillArea(Int_t n, Double_t *x, Double_t *y, Option_t *option="") = 0;
   virtual void     PaintFillAreaNDC(Int_t n, Double_t *x, Double_t *y, Option_t *option="") = 0;
   virtual void     PaintPadFrame(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax) = 0;
   virtual void     PaintLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2) = 0;
   virtual void     PaintLineNDC(Double_t u1, Double_t v1,Double_t u2, Double_t v2) = 0;
   virtual void     PaintLine3D(Float_t *p1, Float_t *p2) = 0;
   virtual void     PaintLine3D(Double_t *p1, Double_t *p2) = 0;
   virtual void     PaintPolyLine(Int_t n, Float_t *x, Float_t *y, Option_t *option="") = 0;
   virtual void     PaintPolyLine(Int_t n, Double_t *x, Double_t *y, Option_t *option="") = 0;
   virtual void     PaintPolyLine3D(Int_t n, Double_t *p) = 0;
   virtual void     PaintPolyLineNDC(Int_t n, Double_t *x, Double_t *y, Option_t *option="") = 0;
   virtual void     PaintPolyMarker(Int_t n, Float_t *x, Float_t *y, Option_t *option="") = 0;
   virtual void     PaintPolyMarker(Int_t n, Double_t *x, Double_t *y, Option_t *option="") = 0;
   virtual void     PaintModified() = 0;
   virtual void     PaintText(Double_t x, Double_t y, const char *text) = 0;
   virtual void     PaintText(Double_t x, Double_t y, const wchar_t *text) = 0;
   virtual void     PaintTextNDC(Double_t u, Double_t v, const char *text) = 0;
   virtual void     PaintTextNDC(Double_t u, Double_t v, const wchar_t *text) = 0;
   virtual Double_t PixeltoX(Int_t px) = 0;
   virtual Double_t PixeltoY(Int_t py) = 0;
           void     Pop() override = 0;
           void     Print(const char *filename="") const override = 0;
   virtual void     Print(const char *filename, Option_t *option) = 0;
   virtual void     Range(Double_t x1, Double_t y1, Double_t x2, Double_t y2) = 0;
   virtual void     RangeAxis(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax) = 0;
           void     RecursiveRemove(TObject *obj) override = 0;
   virtual void     RedrawAxis(Option_t *option="") = 0;
   virtual void     ResetView3D(TObject *view=0) = 0;
   virtual void     ResizePad(Option_t *option="") = 0;
           void     SaveAs(const char *filename="",Option_t *option="") const override = 0;
   virtual void     SetBatch(Bool_t batch=kTRUE) = 0;
   virtual void     SetBorderMode(Short_t bordermode) = 0;
   virtual void     SetBorderSize(Short_t bordersize) = 0;
   virtual void     SetCanvas(TCanvas *c) = 0;
   virtual void     SetCanvasSize(UInt_t ww, UInt_t wh) = 0;
   virtual void     SetCrosshair(Int_t crhair=1) = 0;
   virtual void     SetCursor(ECursor cursor) = 0;
   virtual void     SetDoubleBuffer(Int_t mode=1) = 0;
   virtual void     SetEditable(Bool_t mode=kTRUE) = 0;
   virtual void     SetFixedAspectRatio(Bool_t fixed = kTRUE) = 0;
   virtual void     SetGrid(Int_t valuex = 1, Int_t valuey = 1) = 0;
   virtual void     SetGridx(Int_t value = 1) = 0;
   virtual void     SetGridy(Int_t value = 1) = 0;
   virtual void     SetLogx(Int_t value = 1) = 0;
   virtual void     SetLogy(Int_t value = 1) = 0;
   virtual void     SetLogz(Int_t value = 1) = 0;
   virtual void     SetPad(const char *name, const char *title,
                           Double_t xlow, Double_t ylow, Double_t xup,
                           Double_t yup, Color_t color=35,
                           Short_t bordersize=5, Short_t bordermode=-1) = 0;
   virtual void     SetPad(Double_t xlow, Double_t ylow, Double_t xup, Double_t yup) = 0;
   virtual void     SetAttFillPS(Color_t color, Style_t style) = 0;
   virtual void     SetAttLinePS(Color_t color, Style_t style, Width_t lwidth) = 0;
   virtual void     SetAttMarkerPS(Color_t color, Style_t style, Size_t msize) = 0;
   virtual void     SetAttTextPS(Int_t align, Float_t angle, Color_t color, Style_t font, Float_t tsize) = 0;
   virtual void     SetName(const char *name) = 0;
   virtual void     SetSelected(TObject *obj) = 0;
   virtual void     SetTicks(Int_t valuex = 1, Int_t valuey = 1) = 0;
   virtual void     SetTickx(Int_t value = 1) = 0;
   virtual void     SetTicky(Int_t value = 1) = 0;
   virtual void     SetTitle(const char *title="") = 0;
   virtual void     SetTheta(Double_t theta=30) = 0;
   virtual void     SetPhi(Double_t phi=30) = 0;
   virtual void     SetToolTipText(const char *text, Long_t delayms = 1000) = 0;
   virtual void     SetVertical(Bool_t vert=kTRUE) = 0;
   virtual void     SetView(TView *view=0) = 0;
   virtual void     SetViewer3D(TVirtualViewer3D * /*viewer3d*/) {}
   virtual void     ShowGuidelines(TObject *object, const Int_t event, const char mode = 'i', const bool cling = true) = 0;
   virtual TObject *WaitPrimitive(const char *pname="", const char *emode="") = 0;
   virtual void     Update() = 0;
   virtual Int_t    UtoAbsPixel(Double_t u) const = 0;
   virtual Int_t    VtoAbsPixel(Double_t v) const = 0;
   virtual Int_t    UtoPixel(Double_t u) const = 0;
   virtual Int_t    VtoPixel(Double_t v) const = 0;
   virtual Int_t    XtoAbsPixel(Double_t x) const = 0;
   virtual Int_t    YtoAbsPixel(Double_t y) const = 0;
   virtual Double_t XtoPad(Double_t x) const = 0;
   virtual Double_t YtoPad(Double_t y) const = 0;
   virtual Int_t    XtoPixel(Double_t x) const = 0;
   virtual Int_t    YtoPixel(Double_t y) const = 0;

   virtual Int_t    IncrementPaletteColor(Int_t i, TString opt) = 0;
   virtual Int_t    NextPaletteColor() = 0;

   virtual Bool_t   PlaceBox(TObject *o, Double_t w, Double_t h, Double_t &xl, Double_t &yb) = 0;

   virtual TObject *CreateToolTip(const TBox *b, const char *text, Long_t delayms) = 0;
   virtual void     DeleteToolTip(TObject *tip) = 0;
   virtual void     ResetToolTip(TObject *tip) = 0;
   virtual void     CloseToolTip(TObject *tip) = 0;

   virtual TVirtualViewer3D *GetViewer3D(Option_t * type = "") = 0;
   virtual Bool_t            HasViewer3D() const = 0;
   virtual void              ReleaseViewer3D(Option_t * type = "")  = 0;

   virtual Int_t               GetGLDevice() = 0;
   virtual void                SetCopyGLDevice(Bool_t copy) = 0;
   virtual TVirtualPadPainter *GetPainter() = 0;

   virtual Bool_t PadInSelectionMode() const;
   virtual Bool_t PadInHighlightMode() const;

   virtual void PushTopLevelSelectable(TObject *top);
   virtual void PushSelectableObject(TObject *obj);
   virtual void PopTopLevelSelectable();

   static TVirtualPad *&Pad();

   ClassDefOverride(TVirtualPad,3)  //Abstract base class for Pads and Canvases
};

//
//Small scope-guard class to add/remove object's into pad's stack of selectable objects.
//Does nothing, unless you implement non-standard picking.
//

class TPickerStackGuard {
public:
   TPickerStackGuard(TObject *obj);
   ~TPickerStackGuard();

private:
   TPickerStackGuard(const TPickerStackGuard &rhs) = delete;
   TPickerStackGuard &operator = (const TPickerStackGuard &rhs) = delete;
};


#ifndef __CINT__
#define gPad (TVirtualPad::Pad())
#endif
R__EXTERN Int_t (*gThreadXAR)(const char *xact, Int_t nb, void **ar, Int_t *iret);

#endif
