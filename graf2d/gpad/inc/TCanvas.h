// @(#)root/gpad:$Id$
// Author: Rene Brun   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TCanvas
#define ROOT_TCanvas

#include "TPad.h"

#include "TAttCanvas.h"

class TCanvasImp;
class TContextMenu;
class TControlBar;

class TCanvas : public TPad {

friend class TCanvasImp;
friend class TThread;
friend class TInterpreter;

protected:
   TAttCanvas    fCatt;            ///< Canvas attributes
   TString       fDISPLAY;         ///< Name of destination screen
   Size_t        fXsizeUser;       ///< User specified size of canvas along X in CM
   Size_t        fYsizeUser;       ///< User specified size of canvas along Y in CM
   Size_t        fXsizeReal;       ///< Current size of canvas along X in CM
   Size_t        fYsizeReal;       ///< Current size of canvas along Y in CM
   Color_t       fHighLightColor;  ///< Highlight color of active pad
   Int_t         fDoubleBuffer;    ///< Double buffer flag (0=off, 1=on)
   Int_t         fWindowTopX;      ///< Top X position of window (in pixels)
   Int_t         fWindowTopY;      ///< Top Y position of window (in pixels)
   UInt_t        fWindowWidth;     ///< Width of window (including borders, etc.)
   UInt_t        fWindowHeight;    ///< Height of window (including menubar, borders, etc.)
   UInt_t        fCw;              ///< Width of the canvas along X (pixels)
   UInt_t        fCh;              ///< Height of the canvas along Y (pixels)
   Int_t         fEvent;           ///<! Type of current or last handled event
   Int_t         fEventX;          ///<! Last X mouse position in canvas
   Int_t         fEventY;          ///<! Last Y mouse position in canvas
   Int_t         fCanvasID;        ///<! Canvas identifier
   TObject      *fSelected;        ///<! Currently selected object
   TObject      *fClickSelected;   ///<! Currently click-selected object
   Int_t         fSelectedX;       ///<! X of selected object
   Int_t         fSelectedY;       ///<! Y of selected object
   TString       fSelectedOpt;     ///<! Drawing option of selected object
   TPad         *fSelectedPad;     ///<! Pad containing currently selected object
   TPad         *fClickSelectedPad;///<! Pad containing currently click-selected object
   TPad         *fPadSave;         ///<! Pointer to saved pad in HandleInput
   TCanvasImp   *fCanvasImp;       ///<! Window system specific canvas implementation
   TContextMenu *fContextMenu;     ///<! Context menu pointer
   Bool_t        fBatch;           ///<! True when in batchmode
   Bool_t        fUpdating;        ///<! True when Updating the canvas
   Bool_t        fRetained;        ///< Retain structure flag
   Bool_t        fUseGL;           ///<! True when rendering is with GL
   Bool_t        fDrawn;           ///<! Set to True when the Draw method is called
   //
   TVirtualPadPainter *fPainter;   ///<! Canvas (pad) painter.

   static Bool_t fgIsFolder;       ///< Indicates if canvas can be browsed as a folder

private:
   TCanvas(const TCanvas &canvas);  // cannot copy canvas, use TObject::Clone()
   TCanvas &operator=(const TCanvas &rhs);  // idem
   void     Build();
   void     CopyPixmaps();
   void     DrawEventStatus(Int_t event, Int_t x, Int_t y, TObject *selected);
   void     RunAutoExec();

   //Initialize PadPainter.
   void     CreatePainter();

protected:
   virtual void ExecuteEvent(Int_t event, Int_t px, Int_t py);
   //-- used by friend TThread class
   void Init();

public:
   // TCanvas status bits
   enum {
      kShowEventStatus  = BIT(15),
      kAutoExec         = BIT(16),
      kMenuBar          = BIT(17),
      kShowToolBar      = BIT(18),
      kShowEditor       = BIT(19),
      kMoveOpaque       = BIT(20),
      kResizeOpaque     = BIT(21),
      kIsGrayscale      = BIT(22),
      kShowToolTips     = BIT(23)
   };

   TCanvas(Bool_t build=kTRUE);
   TCanvas(const char *name, const char *title="", Int_t form=1);
   TCanvas(const char *name, const char *title, Int_t ww, Int_t wh);
   TCanvas(const char *name, const char *title, Int_t wtopx, Int_t wtopy,
           Int_t ww, Int_t wh);
   TCanvas(const char *name, Int_t ww, Int_t wh, Int_t winid);
   virtual ~TCanvas();

   //-- used by friend TThread class
   void Constructor();
   void Constructor(const char *name, const char *title, Int_t form);
   void Constructor(const char *name, const char *title, Int_t ww, Int_t wh);
   void Constructor(const char *name, const char *title,
           Int_t wtopx, Int_t wtopy, Int_t ww, Int_t wh);
   void Destructor();

   TVirtualPad      *cd(Int_t subpadnumber=0);
   virtual void      Browse(TBrowser *b);
   void              Clear(Option_t *option="");
   void              Close(Option_t *option="");
   virtual void      Delete(Option_t * = "") { MayNotUse("Delete()"); }
   void              DisconnectWidget();  // used by TCanvasImp
   virtual void      Draw(Option_t *option="");
   virtual TObject  *DrawClone(Option_t *option="") const; // *MENU*
   virtual TObject  *DrawClonePad(); // *MENU*
   virtual void      EditorBar();
   void              EmbedInto(Int_t winid, Int_t ww, Int_t wh);
   void              EnterLeave(TPad *prevSelPad, TObject *prevSelObj);
   void              FeedbackMode(Bool_t set);
   void              Flush();
   void              UseCurrentStyle(); // *MENU*
   void              ForceUpdate();
   const char       *GetDISPLAY() const {return fDISPLAY.Data();}
   TContextMenu     *GetContextMenu() const {return fContextMenu;};
   Int_t             GetDoubleBuffer() const {return fDoubleBuffer;}
   Int_t             GetEvent() const { return fEvent; }
   Int_t             GetEventX() const { return fEventX; }
   Int_t             GetEventY() const { return fEventY; }
   Color_t           GetHighLightColor() const { return fHighLightColor; }
   TVirtualPad      *GetPadSave() const { return fPadSave; }
   void              ClearPadSave() { fPadSave = 0; }
   TObject          *GetSelected() const {return fSelected;}
   TObject          *GetClickSelected() const {return fClickSelected;}
   Int_t             GetSelectedX() const {return fSelectedX;}
   Int_t             GetSelectedY() const {return fSelectedY;}
   Option_t         *GetSelectedOpt() const {return fSelectedOpt.Data();}
   TVirtualPad      *GetSelectedPad() const { return fSelectedPad; }
   TVirtualPad      *GetClickSelectedPad() const { return fClickSelectedPad; }
   Bool_t            GetShowEventStatus() const { return TestBit(kShowEventStatus); }
   Bool_t            GetShowToolBar() const { return TestBit(kShowToolBar); }
   Bool_t            GetShowEditor() const { return TestBit(kShowEditor); }
   Bool_t            GetShowToolTips() const { return TestBit(kShowToolTips); }
   Bool_t            GetAutoExec() const { return TestBit(kAutoExec); }
   Size_t            GetXsizeUser() const {return fXsizeUser;}
   Size_t            GetYsizeUser() const {return fYsizeUser;}
   Size_t            GetXsizeReal() const {return fXsizeReal;}
   Size_t            GetYsizeReal() const {return fYsizeReal;}
   Int_t             GetCanvasID() const {return fCanvasID;}
   TCanvasImp       *GetCanvasImp() const {return fCanvasImp;}
   Int_t             GetWindowTopX();
   Int_t             GetWindowTopY();
   UInt_t            GetWindowWidth() const { return fWindowWidth; }
   UInt_t            GetWindowHeight() const { return fWindowHeight; }
   UInt_t            GetWw() const { return fCw; }
   UInt_t            GetWh() const { return fCh; }
   virtual void      GetCanvasPar(Int_t &wtopx, Int_t &wtopy, UInt_t &ww, UInt_t &wh)
                     {wtopx=GetWindowTopX(); wtopy=fWindowTopY; ww=fWindowWidth; wh=fWindowHeight;}
   virtual void      HandleInput(EEventType button, Int_t x, Int_t y);
   Bool_t            HasMenuBar() const { return TestBit(kMenuBar); }
   virtual void      HighlightConnect(const char *slot);
   void              Iconify();
   Bool_t            IsBatch() const { return fBatch; }
   Bool_t            IsDrawn() { return fDrawn; }
   Bool_t            IsFolder() const;
   Bool_t            IsGrayscale();
   Bool_t            IsRetained() const { return fRetained; }
   Bool_t            IsWeb() const;
   virtual void      ls(Option_t *option="") const;
   void              MoveOpaque(Int_t set=1);
   Bool_t            OpaqueMoving() const { return TestBit(kMoveOpaque); }
   Bool_t            OpaqueResizing() const { return TestBit(kResizeOpaque); }
   virtual void      Paint(Option_t *option="");
   virtual TPad     *Pick(Int_t px, Int_t py, TObjLink *&pickobj) { return TPad::Pick(px, py, pickobj); }
   virtual TPad     *Pick(Int_t px, Int_t py, TObject *prevSelObj);
   virtual void      Picked(TPad *selpad, TObject *selected, Int_t event);             // *SIGNAL*
   virtual void      Highlighted(TVirtualPad *pad, TObject *obj, Int_t x, Int_t y);    // *SIGNAL*
   virtual void      ProcessedEvent(Int_t event, Int_t x, Int_t y, TObject *selected); // *SIGNAL*
   virtual void      Selected(TVirtualPad *pad, TObject *obj, Int_t event);            // *SIGNAL*
   virtual void      Cleared(TVirtualPad *pad);                                        // *SIGNAL*
   virtual void      Closed();                                                         // *SIGNAL*
   void              RaiseWindow();
   void              ResetDrawn() { fDrawn=kFALSE; }
   virtual void      Resize(Option_t *option="");
   void              ResizeOpaque(Int_t set=1);
   void              SaveSource(const char *filename="", Option_t *option="");
   void              SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void      SetCursor(ECursor cursor);
   virtual void      SetDoubleBuffer(Int_t mode=1);
   virtual void      SetName(const char *name="");
   virtual void      SetFixedAspectRatio(Bool_t fixed = kTRUE);  // *TOGGLE*
   void              SetGrayscale(Bool_t set = kTRUE); // *TOGGLE* *GETTER=IsGrayscale
   void              SetWindowPosition(Int_t x, Int_t y);
   void              SetWindowSize(UInt_t ww, UInt_t wh);
   void              SetCanvasImp(TCanvasImp *i) { fCanvasImp = i; }
   void              SetCanvasSize(UInt_t ww, UInt_t wh); // *MENU*
   void              SetHighLightColor(Color_t col) { fHighLightColor = col; }
   void              SetSelected(TObject *obj);
   void              SetClickSelected(TObject *obj) { fClickSelected = obj; }
   void              SetSelectedPad(TPad *pad) { fSelectedPad = pad; }
   void              SetClickSelectedPad(TPad *pad) { fClickSelectedPad = pad; }
   void              Show();
   virtual void      Size(Float_t xsizeuser=0, Float_t ysizeuser=0);
   void              SetBatch(Bool_t batch=kTRUE);
   static  void      SetFolder(Bool_t isfolder=kTRUE);
   void              SetPadSave(TPad *pad) {fPadSave = pad;}
   bool              SetRealAspectRatio(const Int_t axis = 1); // *MENU*
   void              SetRetained(Bool_t retained=kTRUE) { fRetained=retained;}
   void              SetTitle(const char *title="");
   virtual void      ToggleEventStatus();
   virtual void      ToggleAutoExec();
   virtual void      ToggleToolBar();
   virtual void      ToggleEditor();
   virtual void      ToggleToolTips();
   virtual void      Update();

   Bool_t              UseGL() const { return fUseGL; }
   void                SetSupportGL(Bool_t support) {fUseGL = support;}
   TVirtualPadPainter *GetCanvasPainter();
   void                DeleteCanvasPainter();

   static TCanvas   *MakeDefCanvas();
   static Bool_t     SupportAlpha();

   ClassDef(TCanvas,8)  //Graphics canvas
};

#endif
