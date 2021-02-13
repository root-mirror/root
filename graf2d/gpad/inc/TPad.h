// @(#)root/gpad:$Id$
// Author: Rene Brun   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPad
#define ROOT_TPad

#include "TVirtualPad.h"
#include "TAttBBox2D.h"

class TVirtualViewer3D;
class TVirtualPadPainter;
class TBrowser;
class TBox;
class TLegend;
class TArrow;
class TPoint;

class TPad : public TVirtualPad, public TAttBBox2D {

private:
   TObject      *fTip;              ///<! tool tip associated with box

protected:
   Double_t      fX1;               ///<  X of lower X coordinate
   Double_t      fY1;               ///<  Y of lower Y coordinate
   Double_t      fX2;               ///<  X of upper X coordinate
   Double_t      fY2;               ///<  Y of upper Y coordinate

   Double_t      fXtoAbsPixelk;     ///<  Conversion coefficient for X World to absolute pixel
   Double_t      fXtoPixelk;        ///<  Conversion coefficient for X World to pixel
   Double_t      fXtoPixel;         ///<    xpixel = fXtoPixelk + fXtoPixel*xworld
   Double_t      fYtoAbsPixelk;     ///<  Conversion coefficient for Y World to absolute pixel
   Double_t      fYtoPixelk;        ///<  Conversion coefficient for Y World to pixel
   Double_t      fYtoPixel;         ///<    ypixel = fYtoPixelk + fYtoPixel*yworld

   Double_t      fUtoAbsPixelk;     ///<  Conversion coefficient for U NDC to absolute pixel
   Double_t      fUtoPixelk;        ///<  Conversion coefficient for U NDC to pixel
   Double_t      fUtoPixel;         ///<    xpixel = fUtoPixelk + fUtoPixel*undc
   Double_t      fVtoAbsPixelk;     ///<  Conversion coefficient for V NDC to absolute pixel
   Double_t      fVtoPixelk;        ///<  Conversion coefficient for V NDC to pixel
   Double_t      fVtoPixel;         ///<    ypixel = fVtoPixelk + fVtoPixel*vndc

   Double_t      fAbsPixeltoXk;     ///<  Conversion coefficient for absolute pixel to X World
   Double_t      fPixeltoXk;        ///<  Conversion coefficient for pixel to X World
   Double_t      fPixeltoX;         ///<     xworld = fPixeltoXk + fPixeltoX*xpixel
   Double_t      fAbsPixeltoYk;     ///<  Conversion coefficient for absolute pixel to Y World
   Double_t      fPixeltoYk;        ///<  Conversion coefficient for pixel to Y World
   Double_t      fPixeltoY;         ///<     yworld = fPixeltoYk + fPixeltoY*ypixel

   Double_t      fXlowNDC;          ///<  X bottom left corner of pad in NDC [0,1]
   Double_t      fYlowNDC;          ///<  Y bottom left corner of pad in NDC [0,1]
   Double_t      fXUpNDC;
   Double_t      fYUpNDC;
   Double_t      fWNDC;             ///<  Width of pad along X in Normalized Coordinates (NDC)
   Double_t      fHNDC;             ///<  Height of pad along Y in Normalized Coordinates (NDC)

   Double_t      fAbsXlowNDC;       ///<  Absolute X top left corner of pad in NDC [0,1]
   Double_t      fAbsYlowNDC;       ///<  Absolute Y top left corner of pad in NDC [0,1]
   Double_t      fAbsWNDC;          ///<  Absolute Width of pad along X in NDC
   Double_t      fAbsHNDC;          ///<  Absolute Height of pad along Y in NDC

   Double_t      fUxmin;            ///<  Minimum value on the X axis
   Double_t      fUymin;            ///<  Minimum value on the Y axis
   Double_t      fUxmax;            ///<  Maximum value on the X axis
   Double_t      fUymax;            ///<  Maximum value on the Y axis

   Double_t      fTheta;            ///<  theta angle to view as lego/surface
   Double_t      fPhi;              ///<  phi angle   to view as lego/surface

   Double_t      fAspectRatio;      ///<  ratio of w/h in case of fixed ratio

   Int_t         fPixmapID;         ///<! Off-screen pixmap identifier
   Int_t         fGLDevice;         ///<! OpenGL off-screen pixmap identifier
   Bool_t        fCopyGLDevice;     ///<!
   Bool_t        fEmbeddedGL;       ///<!
   Int_t         fNumber;           ///<  pad number identifier
   Int_t         fTickx;            ///<  Set to 1 if tick marks along X
   Int_t         fTicky;            ///<  Set to 1 if tick marks along Y
   Int_t         fLogx;             ///<  (=0 if X linear scale, =1 if log scale)
   Int_t         fLogy;             ///<  (=0 if Y linear scale, =1 if log scale)
   Int_t         fLogz;             ///<  (=0 if Z linear scale, =1 if log scale)
   Int_t         fPadPaint;         ///<  Set to 1 while painting the pad
   Int_t         fCrosshair;        ///<  Crosshair type (0 if no crosshair requested)
   Int_t         fCrosshairPos;     ///<  Position of crosshair
   Short_t       fBorderSize;       ///<  pad bordersize in pixels
   Short_t       fBorderMode;       ///<  Bordermode (-1=down, 0 = no border, 1=up)
   Bool_t        fModified;         ///<  Set to true when pad is modified
   Bool_t        fGridx;            ///<  Set to true if grid along X
   Bool_t        fGridy;            ///<  Set to true if grid along Y
   Bool_t        fAbsCoord;         ///<  Use absolute coordinates
   Bool_t        fEditable;         ///<  True if canvas is editable
   Bool_t        fFixedAspectRatio; ///<  True if fixed aspect ratio
   TPad         *fMother;           ///<! pointer to mother of the list
   TCanvas      *fCanvas;           ///<! Pointer to mother canvas
   TList        *fPrimitives;       ///<->List of primitives (subpads)
   TList        *fExecs;            ///<  List of commands to be executed when a pad event occurs
   TString       fName;             ///<  Pad name
   TString       fTitle;            ///<  Pad title
   TFrame       *fFrame;            ///<! Pointer to 2-D frame (if one exists)
   TView        *fView;             ///<! Pointer to 3-D view (if one exists)
   TObject      *fPadPointer;       ///<! free pointer
   TObject      *fPadView3D;        ///<! 3D View of this TPad
   static Int_t  fgMaxPickDistance; ///<  Maximum Pick Distance
   Int_t         fNumPaletteColor;  ///<  Number of objects with an automatic color
   Int_t         fNextPaletteColor; ///<  Next automatic color
   Bool_t       *fCollideGrid;      ///<! Grid used to find empty space when adding a box (Legend) in a pad
   Int_t         fCGnx;             ///<! Size of the collide grid along x
   Int_t         fCGny;             ///<! Size of the collide grid along y

   // 3D Viewer support
   TVirtualViewer3D *fViewer3D;     ///<! Current 3D viewer

   void          DestroyExternalViewer3D();
   virtual Int_t DistancetoPrimitive(Int_t px, Int_t py);
   virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual void  HideToolTip(Int_t event);
   void          PaintBorder(Color_t color, Bool_t tops);
   virtual void  PaintBorderPS(Double_t xl,Double_t yl,Double_t xt,Double_t yt,Int_t bmode,Int_t bsize,Int_t dark,Int_t light);
   void          PaintDate();
   virtual void  SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void  SetBatch(Bool_t batch=kTRUE);

private:
   TPad(const TPad &pad);  // cannot copy pads, use TObject::Clone()
   TPad &operator=(const TPad &rhs);  // idem

   void CopyBackgroundPixmap(Int_t x, Int_t y);
   void CopyBackgroundPixmaps(TPad *start, TPad *stop, Int_t x, Int_t y);
   void DrawDist(Rectangle_t aBBox, Rectangle_t bBBox, char mode);

   Bool_t            Collide(Int_t i, Int_t j, Int_t w, Int_t h);
   void              FillCollideGrid(TObject *o);
   void              FillCollideGridTBox(TObject *o);
   void              FillCollideGridTFrame(TObject *o);
   void              FillCollideGridTGraph(TObject *o);
   void              FillCollideGridTH1(TObject *o);
   void              LineNotFree(Int_t x1, Int_t x2, Int_t y1, Int_t y2);

public:
   // TPad status bits
   enum {
      kFraming      = BIT(6),  ///< Frame is requested
      kHori         = BIT(9),  ///< Pad is horizontal
      kClipFrame    = BIT(10), ///< Clip on frame
      kPrintingPS   = BIT(11), ///< PS Printing
      kCannotMove   = BIT(12), ///< Fixed position
      kClearAfterCR = BIT(14)  ///< Clear after CR
   };

   TPad();
   TPad(const char *name, const char *title, Double_t xlow,
        Double_t ylow, Double_t xup, Double_t yup,
        Color_t color=-1, Short_t bordersize=-1, Short_t bordermode=-2);
   virtual ~TPad();
   void              AbsCoordinates(Bool_t set) { fAbsCoord = set; }
   Double_t          AbsPixeltoX(Int_t px) {return fAbsPixeltoXk + px*fPixeltoX;}
   Double_t          AbsPixeltoY(Int_t py) {return fAbsPixeltoYk + py*fPixeltoY;}
   virtual void      AbsPixeltoXY(Int_t xpixel, Int_t ypixel, Double_t &x, Double_t &y);
   virtual void      AddExec(const char *name, const char *command);
   virtual void      AutoExec();
   virtual void      Browse(TBrowser *b);
   virtual TLegend  *BuildLegend(Double_t x1=0.3, Double_t y1=0.21, Double_t x2=0.3, Double_t y2=0.21, const char *title="", Option_t *option = ""); // *MENU*
   TVirtualPad*      cd(Int_t subpadnumber=0); // *MENU*
   void              Clear(Option_t *option="");
   virtual Int_t     Clip(Float_t *x, Float_t *y, Float_t xclipl, Float_t yclipb, Float_t xclipr, Float_t yclipt);
   virtual Int_t     Clip(Double_t *x, Double_t *y, Double_t xclipl, Double_t yclipb, Double_t xclipr, Double_t yclipt);
   virtual Int_t     ClippingCode(Double_t x, Double_t y, Double_t xcl1, Double_t ycl1, Double_t xcl2, Double_t ycl2);
   virtual Int_t     ClipPolygon(Int_t n, Double_t *x, Double_t *y, Int_t nn, Double_t *xc, Double_t *yc, Double_t xclipl, Double_t yclipb, Double_t xclipr, Double_t yclipt);
   virtual void      Close(Option_t *option="");
   virtual void      Closed() { Emit("Closed()"); } // *SIGNAL*
   virtual void      CopyPixmap();
   virtual void      CopyPixmaps();
   virtual void      DeleteExec(const char *name);
   virtual void      Divide(Int_t nx=1, Int_t ny=1, Float_t xmargin=0.01, Float_t ymargin=0.01, Int_t color=0); // *MENU*
   virtual void      DivideSquare(Int_t n, Float_t xmargin=0.01, Float_t ymargin=0.01, Int_t color=0);
   virtual void      Draw(Option_t *option="");
   virtual void      DrawClassObject(const TObject *obj, Option_t *option="");
   static  void      DrawColorTable();
   virtual void      DrawCrosshair();
   TH1F             *DrawFrame(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, const char *title="");
   virtual void      ExecuteEventAxis(Int_t event, Int_t px, Int_t py, TAxis *axis);
   virtual TObject  *FindObject(const char *name) const;
   virtual TObject  *FindObject(const TObject *obj) const;
   virtual void      UseCurrentStyle();  // *MENU*
   virtual Short_t   GetBorderMode() const { return fBorderMode;}
   virtual Short_t   GetBorderSize() const { return fBorderSize;}
   Int_t             GetCrosshair() const;
   virtual Int_t     GetCanvasID() const;
   virtual TCanvasImp *GetCanvasImp() const;
   TFrame           *GetFrame();
   virtual Int_t     GetEvent() const;
   virtual Int_t     GetEventX() const;
   virtual Int_t     GetEventY() const;
   virtual Color_t   GetHighLightColor() const;
   virtual void      GetRange(Double_t &x1, Double_t &y1, Double_t &x2, Double_t &y2);
   virtual void      GetRangeAxis(Double_t &xmin, Double_t &ymin, Double_t &xmax, Double_t &ymax);
   virtual void      GetPadPar(Double_t &xlow, Double_t &ylow, Double_t &xup, Double_t &yup);
   Double_t          GetXlowNDC() const {return fXlowNDC;}
   Double_t          GetYlowNDC() const {return fYlowNDC;}
   /// Get width of pad along X in Normalized Coordinates (NDC)
   Double_t          GetWNDC() const {return fWNDC;}
   /// Get height of pad along Y in Normalized Coordinates (NDC)
   Double_t          GetHNDC() const {return fHNDC;}
   virtual UInt_t    GetWw() const;
   virtual UInt_t    GetWh() const;
   Double_t          GetAbsXlowNDC() const {return fAbsXlowNDC;}
   Double_t          GetAbsYlowNDC() const {return fAbsYlowNDC;}
   Double_t          GetAbsWNDC() const {return fAbsWNDC;}
   Double_t          GetAbsHNDC() const {return fAbsHNDC;}
   Double_t          GetAspectRatio() const { return fAspectRatio; }
   Double_t          GetPhi() const   {return fPhi;}
   Double_t          GetTheta() const {return fTheta;}
   ///Returns the minimum x-coordinate value visible on the pad. If log axis the returned value is in decades.
   Double_t          GetUxmin() const {return fUxmin;}
   ///Returns the minimum y-coordinate value visible on the pad. If log axis the returned value is in decades.
   Double_t          GetUymin() const {return fUymin;}
   ///Returns the maximum x-coordinate value visible on the pad. If log axis the returned value is in decades.
   Double_t          GetUxmax() const {return fUxmax;}
   ///Returns the maximum y-coordinate value visible on the pad. If log axis the returned value is in decades.
   Double_t          GetUymax() const {return fUymax;}
   Bool_t            GetGridx() const {return fGridx;}
   Bool_t            GetGridy() const {return fGridy;}
   Int_t             GetNumber() const {return fNumber;}
   Int_t             GetTickx() const {return fTickx;}
   Int_t             GetTicky() const {return fTicky;}
   Double_t          GetX1() const { return fX1; }
   Double_t          GetX2() const { return fX2; }
   Double_t          GetY1() const { return fY1; }
   Double_t          GetY2() const { return fY2; }
   static Int_t      GetMaxPickDistance();
   TList            *GetListOfPrimitives() const {return fPrimitives;}
   TList            *GetListOfExecs() const {return fExecs;}
   virtual TObject  *GetPrimitive(const char *name) const;  //obsolete, use FindObject instead
   virtual TObject  *GetSelected() const;
   virtual TVirtualPad  *GetPad(Int_t subpadnumber) const;
   virtual TObject  *GetPadPointer() const {return fPadPointer;}
   TVirtualPad      *GetPadSave() const;
   TVirtualPad      *GetSelectedPad() const;
   Int_t             GetGLDevice();
   TView            *GetView() const {return fView;}
   TObject          *GetView3D() const {return fPadView3D;}// Return 3D View of this TPad
   Int_t             GetLogx() const {return fLogx;}
   Int_t             GetLogy() const {return fLogy;}
   Int_t             GetLogz() const {return fLogz;}
   virtual TVirtualPad *GetMother() const {return fMother;}
   const char       *GetName() const {return fName.Data();}
   const char       *GetTitle() const {return fTitle.Data();}
   virtual TCanvas  *GetCanvas() const { return fCanvas; }
   virtual TVirtualPad *GetVirtCanvas() const ;
   virtual TVirtualPadPainter *GetPainter();
   Int_t             GetPadPaint() const {return fPadPaint;}
   Int_t             GetPixmapID() const {return fPixmapID;}
   ULong_t           Hash() const { return fName.Hash(); }
   virtual Bool_t    HasCrosshair() const;
   void              HighLight(Color_t col=kRed, Bool_t set=kTRUE);
   Bool_t            HasFixedAspectRatio() const { return fFixedAspectRatio; }
   virtual Bool_t    IsBatch() const;
   virtual Bool_t    IsEditable() const {return fEditable;}
   Bool_t            IsFolder() const {return kTRUE;}
   Bool_t            IsModified() const {return fModified;}
   virtual Bool_t    IsRetained() const;
   virtual Bool_t    IsVertical() const {return !TestBit(kHori);}
   virtual void      ls(Option_t *option="") const;
   void              Modified(Bool_t flag=1);  // *SIGNAL*
   virtual Bool_t    OpaqueMoving() const;
   virtual Bool_t    OpaqueResizing() const;
   Double_t          PadtoX(Double_t x) const;
   Double_t          PadtoY(Double_t y) const;
   virtual void      Paint(Option_t *option="");
   void              PaintBox(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="");
   void              PaintFillArea(Int_t n, Float_t *x, Float_t *y, Option_t *option=""); // Obsolete
   void              PaintFillArea(Int_t n, Double_t *x, Double_t *y, Option_t *option="");
   void              PaintFillAreaNDC(Int_t n, Double_t *x, Double_t *y, Option_t *option="");
   void              PaintFillAreaHatches(Int_t n, Double_t *x, Double_t *y, Int_t FillStyle);
   void              PaintHatches(Double_t dy, Double_t angle, Int_t nn, Double_t *xx, Double_t *yy);
   void              PaintPadFrame(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax);
   void              PaintLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2);
   void              PaintLineNDC(Double_t u1, Double_t v1,Double_t u2, Double_t v2);
   void              PaintLine3D(Float_t *p1, Float_t *p2);
   void              PaintLine3D(Double_t *p1, Double_t *p2);
   void              PaintPolyLine(Int_t n, Float_t *x, Float_t *y, Option_t *option="");
   void              PaintPolyLine(Int_t n, Double_t *x, Double_t *y, Option_t *option="");
   void              PaintPolyLine3D(Int_t n, Double_t *p);
   void              PaintPolyLineNDC(Int_t n, Double_t *x, Double_t *y, Option_t *option="");
   void              PaintPolyMarker(Int_t n, Float_t *x, Float_t *y, Option_t *option="");
   void              PaintPolyMarker(Int_t n, Double_t *x, Double_t *y, Option_t *option="");
   virtual void      PaintModified();
   void              PaintText(Double_t x, Double_t y, const char *text);
   void              PaintText(Double_t x, Double_t y, const wchar_t *text);
   void              PaintTextNDC(Double_t u, Double_t v, const char *text);
   void              PaintTextNDC(Double_t u, Double_t v, const wchar_t *text);
   virtual TPad     *Pick(Int_t px, Int_t py, TObjLink *&pickobj);
   Double_t          PixeltoX(Int_t px);
   Double_t          PixeltoY(Int_t py);
   virtual void      PixeltoXY(Int_t xpixel, Int_t ypixel, Double_t &x, Double_t &y);
   virtual void      Pop();  // *MENU*
   virtual void      Print(const char *filename="") const;
   virtual void      Print(const char *filename, Option_t *option);
   virtual void      Range(Double_t x1, Double_t y1, Double_t x2, Double_t y2); // *MENU* *ARGS={x1=>fX1,y1=>fY1,x2=>fX2,y2=>fY2}
   virtual void      RangeChanged() { Emit("RangeChanged()"); } // *SIGNAL*
   virtual void      RangeAxis(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax);
   virtual void      RangeAxisChanged() { Emit("RangeAxisChanged()"); } // *SIGNAL*
   virtual void      RecursiveRemove(TObject *obj);
   virtual void      RedrawAxis(Option_t *option="");
   virtual void      ResetView3D(TObject *view=0){fPadView3D=view;}
   virtual void      ResizePad(Option_t *option="");
   virtual void      Resized() { Emit("Resized()"); } // *SIGNAL*
   virtual void      SaveAs(const char *filename="",Option_t *option="") const; // *MENU*
   virtual void      SetBorderMode(Short_t bordermode) {fBorderMode = bordermode; Modified();} // *MENU*
   virtual void      SetBorderSize(Short_t bordersize) {fBorderSize = bordersize; Modified();} // *MENU*
   void              SetCanvas(TCanvas *c) { fCanvas = c; }
   virtual void      SetCanvasSize(UInt_t ww, UInt_t wh);
   virtual void      SetCrosshair(Int_t crhair=1); // *TOGGLE*
   virtual void      SetCursor(ECursor cursor);
   virtual void      SetDoubleBuffer(Int_t mode=1);
   virtual void      SetDrawOption(Option_t *option="");
   virtual void      SetEditable(Bool_t mode=kTRUE); // *TOGGLE*
   virtual void      SetFixedAspectRatio(Bool_t fixed = kTRUE);  // *TOGGLE*
   virtual void      SetGrid(Int_t valuex = 1, Int_t valuey = 1) {fGridx = valuex; fGridy = valuey; Modified();}
   virtual void      SetGridx(Int_t value = 1) {fGridx = value; Modified();} // *TOGGLE*
   virtual void      SetGridy(Int_t value = 1) {fGridy = value; Modified();} // *TOGGLE*
   virtual void      SetFillStyle(Style_t fstyle);
   virtual void      SetLogx(Int_t value = 1); // *TOGGLE*
   virtual void      SetLogy(Int_t value = 1); // *TOGGLE*
   virtual void      SetLogz(Int_t value = 1); // *TOGGLE*
   virtual void      SetNumber(Int_t number) {fNumber = number;}
   virtual void      SetPad(const char *name, const char *title,
                           Double_t xlow, Double_t ylow, Double_t xup,
                           Double_t yup, Color_t color=35,
                           Short_t bordersize=5, Short_t bordermode=-1);
   virtual void      SetPad(Double_t xlow, Double_t ylow, Double_t xup, Double_t yup);
   virtual void      SetAttFillPS(Color_t color, Style_t style);
   virtual void      SetAttLinePS(Color_t color, Style_t style, Width_t lwidth);
   virtual void      SetAttMarkerPS(Color_t color, Style_t style, Size_t msize);
   virtual void      SetAttTextPS(Int_t align, Float_t angle, Color_t color, Style_t font, Float_t tsize);
   static  void      SetMaxPickDistance(Int_t maxPick=5);
   virtual void      SetName(const char *name) {fName = name;} // *MENU*
   virtual void      SetSelected(TObject *obj);
   virtual void      SetTicks(Int_t valuex = 1, Int_t valuey = 1) {fTickx = valuex; fTicky = valuey; Modified();}
   virtual void      SetTickx(Int_t value = 1) {fTickx = value; Modified();} // *TOGGLE*
   virtual void      SetTicky(Int_t value = 1) {fTicky = value; Modified();} // *TOGGLE*
   virtual void      SetTitle(const char *title="") {fTitle = title;}
   virtual void      SetTheta(Double_t theta=30) {fTheta = theta; Modified();}
   virtual void      SetPhi(Double_t phi=30) {fPhi = phi; Modified();}
   virtual void      SetToolTipText(const char *text, Long_t delayms = 1000);
   virtual void      SetVertical(Bool_t vert=kTRUE);
   virtual void      SetView(TView *view = 0);
   virtual void      SetViewer3D(TVirtualViewer3D *viewer3d) {fViewer3D = viewer3d;}

   virtual void      SetGLDevice(Int_t dev) {fGLDevice = dev;}
   virtual void      SetCopyGLDevice(Bool_t copy) {fCopyGLDevice = copy;}

   virtual void      ShowGuidelines(TObject *object, const Int_t event, const char mode = 'i', const bool cling = true);
   virtual void      Update();
   Int_t             UtoAbsPixel(Double_t u) const {return Int_t(fUtoAbsPixelk + u*fUtoPixel);}
   Int_t             VtoAbsPixel(Double_t v) const {return Int_t(fVtoAbsPixelk + v*fVtoPixel);}
   Int_t             UtoPixel(Double_t u) const;
   Int_t             VtoPixel(Double_t v) const;
   virtual TObject  *WaitPrimitive(const char *pname="", const char *emode="");
   Int_t             XtoAbsPixel(Double_t x) const;
   Int_t             YtoAbsPixel(Double_t y) const;
   Double_t          XtoPad(Double_t x) const;
   Double_t          YtoPad(Double_t y) const;
   Int_t             XtoPixel(Double_t x) const;
   Int_t             YtoPixel(Double_t y) const;
   virtual void      XYtoAbsPixel(Double_t x, Double_t y, Int_t &xpixel, Int_t &ypixel) const;
   virtual void      XYtoPixel(Double_t x, Double_t y, Int_t &xpixel, Int_t &ypixel) const;

   virtual TObject  *CreateToolTip(const TBox *b, const char *text, Long_t delayms);
   virtual void      DeleteToolTip(TObject *tip);
   virtual void      ResetToolTip(TObject *tip);
   virtual void      CloseToolTip(TObject *tip);

   Int_t             IncrementPaletteColor(Int_t i, TString opt);
   Int_t             NextPaletteColor();

   void              DrawCollideGrid();
   Bool_t            PlaceBox(TObject *o, Double_t w, Double_t h, Double_t &xl, Double_t &yb);

   virtual void      x3d(Option_t *type=""); // Depreciated

   virtual TVirtualViewer3D *GetViewer3D(Option_t * type = "");
   virtual Bool_t            HasViewer3D() const { return (fViewer3D); }
   virtual void              ReleaseViewer3D(Option_t * type = "");

   virtual Rectangle_t  GetBBox();
   virtual TPoint       GetBBoxCenter();
   virtual void         SetBBoxCenter(const TPoint &p);
   virtual void         SetBBoxCenterX(const Int_t x);
   virtual void         SetBBoxCenterY(const Int_t y);
   virtual void         SetBBoxX1(const Int_t x);
   virtual void         SetBBoxX2(const Int_t x);
   virtual void         SetBBoxY1(const Int_t y);
   virtual void         SetBBoxY2(const Int_t y);

   virtual void      RecordPave(const TObject *obj);              // *SIGNAL*
   virtual void      RecordLatex(const TObject *obj);             // *SIGNAL*
   virtual void      EventPave() { Emit("EventPave()"); }         // *SIGNAL*
   virtual void      StartEditing() { Emit("StartEditing()"); }   // *SIGNAL*

   ClassDef(TPad,13)  //A Graphics pad
};


//______________________________________________________________________________
inline void TPad::Modified(Bool_t flag)
{
   if (!fModified && flag) Emit("Modified()");
   fModified = flag;
}


//______________________________________________________________________________
inline void TPad::AbsPixeltoXY(Int_t xpixel, Int_t ypixel, Double_t &x, Double_t &y)
{
   x = AbsPixeltoX(xpixel);
   y = AbsPixeltoY(ypixel);
}


//______________________________________________________________________________
inline Double_t TPad::PixeltoX(Int_t px)
{
   if (fAbsCoord) return fAbsPixeltoXk + px*fPixeltoX;
   else           return fPixeltoXk    + px*fPixeltoX;
}


//______________________________________________________________________________
inline Double_t TPad::PixeltoY(Int_t py)
{
   if (fAbsCoord) return fAbsPixeltoYk + py*fPixeltoY;
   else           return fPixeltoYk    + py*fPixeltoY;
}


//______________________________________________________________________________
inline void TPad::PixeltoXY(Int_t xpixel, Int_t ypixel, Double_t &x, Double_t &y)
{
   x = PixeltoX(xpixel);
   y = PixeltoY(ypixel);
}


//______________________________________________________________________________
inline Int_t TPad::UtoPixel(Double_t u) const
{
   Double_t val;
   if (fAbsCoord) val = fUtoAbsPixelk + u*fUtoPixel;
   else           val = u*fUtoPixel;
   if (val < -kMaxPixel) return -kMaxPixel;
   if (val >  kMaxPixel) return  kMaxPixel;
   return Int_t(val);
}


//______________________________________________________________________________
inline Int_t TPad::VtoPixel(Double_t v) const
{
   Double_t val;
   if (fAbsCoord) val = fVtoAbsPixelk + v*fVtoPixel;
   else           val = fVtoPixelk    + v*fVtoPixel;
   if (val < -kMaxPixel) return -kMaxPixel;
   if (val >  kMaxPixel) return  kMaxPixel;
   return Int_t(val);
}


//______________________________________________________________________________
inline Int_t TPad::XtoAbsPixel(Double_t x) const
{
   Double_t val = fXtoAbsPixelk + x*fXtoPixel;
   if (val < -kMaxPixel) return -kMaxPixel;
   if (val >  kMaxPixel) return  kMaxPixel;
   return Int_t(val);
}


//______________________________________________________________________________
inline Int_t TPad::XtoPixel(Double_t x) const
{
   Double_t val;
   if (fAbsCoord) val = fXtoAbsPixelk + x*fXtoPixel;
   else           val = fXtoPixelk    + x*fXtoPixel;
   if (val < -kMaxPixel) return -kMaxPixel;
   if (val >  kMaxPixel) return  kMaxPixel;
   return Int_t(val);
}


//______________________________________________________________________________
inline Int_t TPad::YtoAbsPixel(Double_t y) const
{
   Double_t val = fYtoAbsPixelk + y*fYtoPixel;
   if (val < -kMaxPixel) return -kMaxPixel;
   if (val >  kMaxPixel) return  kMaxPixel;
   return Int_t(val);
}


//______________________________________________________________________________
inline Int_t TPad::YtoPixel(Double_t y) const
{
   Double_t val;
   if (fAbsCoord) val = fYtoAbsPixelk + y*fYtoPixel;
   else           val = fYtoPixelk    + y*fYtoPixel;
   if (val < -kMaxPixel) return -kMaxPixel;
   if (val >  kMaxPixel) return  kMaxPixel;
   return Int_t(val);
}


//______________________________________________________________________________
inline void TPad::XYtoAbsPixel(Double_t x, Double_t y, Int_t &xpixel, Int_t &ypixel) const
{
   xpixel = XtoAbsPixel(x);
   ypixel = YtoAbsPixel(y);
}


//______________________________________________________________________________
inline void TPad::XYtoPixel(Double_t x, Double_t y, Int_t &xpixel, Int_t &ypixel) const
{
   xpixel = XtoPixel(x);
   ypixel = YtoPixel(y);
}


//______________________________________________________________________________
inline void TPad::SetDrawOption(Option_t *)
{ }

#endif

