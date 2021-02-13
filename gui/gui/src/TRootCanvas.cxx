// @(#)root/gui:$Id: b4c21444ab4f787f65b2b44199fc0440c3c2ce81 $
// Author: Fons Rademakers   15/01/98

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRootCanvas                                                          //
//                                                                      //
// This class creates a main window with menubar, scrollbars and a      //
// drawing area. The widgets used are the new native ROOT GUI widgets.  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RConfigure.h"

#include "TRootCanvas.h"
#include "TRootApplication.h"
#include "TRootHelpDialog.h"
#include "TGClient.h"
#include "TGCanvas.h"
#include "TGMenu.h"
#include "TGWidget.h"
#include "TGFileDialog.h"
#include "TGStatusBar.h"
#include "TGTextEditDialogs.h"
#include "TROOT.h"
#include "TClass.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "TClassTree.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TColorWheel.h"
#include "TVirtualX.h"
#include "TApplication.h"
#include "TFile.h"
#include "TInterpreter.h"
#include "TEnv.h"
#include "TMath.h"
#include <iostream>
#include "TGDockableFrame.h"

#include "TG3DLine.h"
#include "TGToolBar.h"
#include "TGToolTip.h"
#include "TVirtualPadEditor.h"
#include "TRootControlBar.h"
#include "TGuiBuilder.h"
#include "TImage.h"
#include "TError.h"
#include "TGDNDManager.h"
#include "TBufferFile.h"
#include "TRootBrowser.h"
#include "TGTab.h"
#include "TGedEditor.h"

#include "TPluginManager.h"
#include "TVirtualGL.h"

#ifdef WIN32
#include "TWin32SplashThread.h"
#endif

#include "HelpText.h"


// Canvas menu command ids
enum ERootCanvasCommands {
   kFileNewCanvas,
   kFileOpen,
   kFileSaveAs,
   kFileSaveAsRoot,
   kFileSaveAsC,
   kFileSaveAsPS,
   kFileSaveAsEPS,
   kFileSaveAsPDF,
   kFileSaveAsGIF,
   kFileSaveAsJPG,
   kFileSaveAsPNG,
   kFileSaveAsTEX,
   kFilePrint,
   kFileCloseCanvas,
   kFileQuit,

   kEditStyle,
   kEditCut,
   kEditCopy,
   kEditPaste,
   kEditClearPad,
   kEditClearCanvas,
   kEditUndo,
   kEditRedo,

   kViewEditor,
   kViewToolbar,
   kViewEventStatus,
   kViewToolTips,
   kViewColors,
   kViewFonts,
   kViewMarkers,
   kViewIconify,
   kViewX3D,
   kViewOpenGL,

   kOptionAutoResize,
   kOptionResizeCanvas,
   kOptionMoveOpaque,
   kOptionResizeOpaque,
   kOptionInterrupt,
   kOptionRefresh,
   kOptionAutoExec,
   kOptionStatistics,
   kOptionHistTitle,
   kOptionFitParams,
   kOptionCanEdit,

   kInspectRoot,
   kClassesTree,
   kFitPanel,
   kToolsBrowser,
   kToolsBuilder,
   kToolsRecorder,

   kHelpAbout,
   kHelpOnCanvas,
   kHelpOnMenus,
   kHelpOnGraphicsEd,
   kHelpOnBrowser,
   kHelpOnObjects,
   kHelpOnPS,

   kToolModify,
   kToolArc,
   kToolLine,
   kToolArrow,
   kToolDiamond,
   kToolEllipse,
   kToolPad,
   kToolPave,
   kToolPLabel,
   kToolPText,
   kToolPsText,
   kToolGraph,
   kToolCurlyLine,
   kToolCurlyArc,
   kToolLatex,
   kToolMarker,
   kToolCutG

};

static const char *gOpenTypes[] = { "ROOT files",   "*.root",
                                    "All files",    "*",
                                    0,              0 };

static const char *gSaveAsTypes[] = { "PDF",          "*.pdf",
                                      "PostScript",   "*.ps",
                                      "Encapsulated PostScript", "*.eps",
                                      "SVG",          "*.svg",
                                      "TeX",          "*.tex",
                                      "GIF",          "*.gif",
                                      "ROOT macros",  "*.C",
                                      "ROOT files",   "*.root",
                                      "XML",          "*.xml",
                                      "PNG",          "*.png",
                                      "XPM",          "*.xpm",
                                      "JPEG",         "*.jpg",
                                      "TIFF",         "*.tiff",
                                      "XCF",          "*.xcf",
                                      "All files",    "*",
                                      0,              0 };

static ToolBarData_t gToolBarData[] = {
   // { filename,      tooltip,            staydown,  id,              button}
   { "newcanvas.xpm",  "New",              kFALSE,    kFileNewCanvas,  0 },
   { "open.xpm",       "Open",             kFALSE,    kFileOpen,       0 },
   { "save.xpm",       "Save As",          kFALSE,    kFileSaveAs,     0 },
   { "printer.xpm",    "Print",            kFALSE,    kFilePrint,      0 },
   { "",               "",                 kFALSE,    -1,              0 },
   { "interrupt.xpm",  "Interrupt",        kFALSE,    kOptionInterrupt,0 },
   { "refresh2.xpm",   "Refresh",          kFALSE,    kOptionRefresh,  0 },
   { "",               "",                 kFALSE,    -1,              0 },
   { "inspect.xpm",    "Inspect",          kFALSE,    kInspectRoot,    0 },
   { "browser.xpm",    "Browser",          kFALSE,    kToolsBrowser, 0 },
   { 0,                0,                  kFALSE,    0,               0 }
};

static ToolBarData_t gToolBarData1[] = {
   { "pointer.xpm",    "Modify",           kFALSE,    kToolModify,     0 },
   { "arc.xpm",        "Arc",              kFALSE,    kToolArc,        0 },
   { "line.xpm",       "Line",             kFALSE,    kToolLine,       0 },
   { "arrow.xpm",      "Arrow",            kFALSE,    kToolArrow,      0 },
   { "diamond.xpm",    "Diamond",          kFALSE,    kToolDiamond,    0 },
   { "ellipse.xpm",    "Ellipse",          kFALSE,    kToolEllipse,    0 },
   { "pad.xpm",        "Pad",              kFALSE,    kToolPad,        0 },
   { "pave.xpm",       "Pave",             kFALSE,    kToolPave,       0 },
   { "pavelabel.xpm",  "Pave Label",       kFALSE,    kToolPLabel,     0 },
   { "pavetext.xpm",   "Pave Text",        kFALSE,    kToolPText,      0 },
   { "pavestext.xpm",  "Paves Text",       kFALSE,    kToolPsText,     0 },
   { "graph.xpm",      "Graph",            kFALSE,    kToolGraph,      0 },
   { "curlyline.xpm",  "Curly Line",       kFALSE,    kToolCurlyLine,  0 },
   { "curlyarc.xpm",   "Curly Arc",        kFALSE,    kToolCurlyArc,   0 },
   { "latex.xpm",      "Text/Latex",       kFALSE,    kToolLatex,      0 },
   { "marker.xpm",     "Marker",           kFALSE,    kToolMarker,     0 },
   { "cut.xpm",        "Graphical Cut",    kFALSE,    kToolCutG,       0 },
   { 0,                0,                  kFALSE,    0,               0 }
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRootContainer                                                       //
//                                                                      //
// Utility class used by TRootCanvas. The TRootContainer is the frame   //
// embedded in the TGCanvas widget. The ROOT graphics goes into this    //
// frame. This class is used to enable input events on this graphics    //
// frame and forward the events to the TRootCanvas handlers.            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TRootContainer : public TGCompositeFrame {
private:
   TRootCanvas  *fCanvas;    // pointer back to canvas imp
public:
   TRootContainer(TRootCanvas *c, Window_t id, const TGWindow *parent);

   Bool_t  HandleButton(Event_t *ev);
   Bool_t  HandleDoubleClick(Event_t *ev)
               { return fCanvas->HandleContainerDoubleClick(ev); }
   Bool_t  HandleConfigureNotify(Event_t *ev)
               { TGFrame::HandleConfigureNotify(ev);
                  return fCanvas->HandleContainerConfigure(ev); }
   Bool_t  HandleKey(Event_t *ev)
               { return fCanvas->HandleContainerKey(ev); }
   Bool_t  HandleMotion(Event_t *ev)
               { return fCanvas->HandleContainerMotion(ev); }
   Bool_t  HandleExpose(Event_t *ev)
               { return fCanvas->HandleContainerExpose(ev); }
   Bool_t  HandleCrossing(Event_t *ev)
               { return fCanvas->HandleContainerCrossing(ev); }
   void    SavePrimitive(std::ostream &out, Option_t * = "");
   void    SetEditable(Bool_t) { }
};

////////////////////////////////////////////////////////////////////////////////
/// Create a canvas container.

TRootContainer::TRootContainer(TRootCanvas *c, Window_t id, const TGWindow *p)
   : TGCompositeFrame(gClient, id, p)
{
   fCanvas = c;

   gVirtualX->GrabButton(fId, kAnyButton, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask |
                         kPointerMotionMask, kNone, kNone);

   AddInput(kKeyPressMask | kKeyReleaseMask | kPointerMotionMask |
            kExposureMask | kStructureNotifyMask | kLeaveWindowMask);
   fEditDisabled = kEditDisable;
}

////////////////////////////////////////////////////////////////////////////////
/// Directly handle scroll mouse buttons (4 and 5), only pass buttons
/// 1, 2 and 3 on to the TCanvas.

Bool_t TRootContainer::HandleButton(Event_t *event)
{
   TGViewPort *vp = (TGViewPort*)fParent;
   UInt_t page = vp->GetHeight()/4;
   Int_t newpos;

   gVirtualX->SetInputFocus(GetMainFrame()->GetId());

   if (event->fCode == kButton4) {
      //scroll up
      newpos = fCanvas->fCanvasWindow->GetVsbPosition() - page;
      if (newpos < 0) newpos = 0;
      fCanvas->fCanvasWindow->SetVsbPosition(newpos);
//      return kTRUE;
   }
   if (event->fCode == kButton5) {
      // scroll down
      newpos = fCanvas->fCanvasWindow->GetVsbPosition() + page;
      fCanvas->fCanvasWindow->SetVsbPosition(newpos);
//      return kTRUE;
   }
   return fCanvas->HandleContainerButton(event);
}

ClassImp(TRootCanvas);

////////////////////////////////////////////////////////////////////////////////
/// Create a basic ROOT canvas.

TRootCanvas::TRootCanvas(TCanvas *c, const char *name, UInt_t width, UInt_t height)
   : TGMainFrame(gClient->GetRoot(), width, height), TCanvasImp(c)
{
   CreateCanvas(name);

   ShowToolBar(kFALSE);
   ShowEditor(kFALSE);

   Resize(width, height);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a basic ROOT canvas.

TRootCanvas::TRootCanvas(TCanvas *c, const char *name, Int_t x, Int_t y, UInt_t width, UInt_t height)
   : TGMainFrame(gClient->GetRoot(), width, height), TCanvasImp(c)
{
   CreateCanvas(name);

   ShowToolBar(kFALSE);
   ShowEditor(kFALSE);

   MoveResize(x, y, width, height);
   SetWMPosition(x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Create the actual canvas.

void TRootCanvas::CreateCanvas(const char *name)
{
   fButton    = 0;
   fAutoFit   = kTRUE;   // check also menu entry
   fEditor    = 0;
   fEmbedded  = kFALSE;

   // Create menus
   fFileSaveMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fFileSaveMenu->AddEntry(Form("%s.&ps",  name), kFileSaveAsPS);
   fFileSaveMenu->AddEntry(Form("%s.&eps", name), kFileSaveAsEPS);
   fFileSaveMenu->AddEntry(Form("%s.p&df", name), kFileSaveAsPDF);
   fFileSaveMenu->AddEntry(Form("%s.&tex", name), kFileSaveAsTEX);
   fFileSaveMenu->AddEntry(Form("%s.&gif", name), kFileSaveAsGIF);

   static Int_t img = 0;

   if (!img) {
      Int_t sav = gErrorIgnoreLevel;
      gErrorIgnoreLevel = kFatal;
      TImage* itmp = TImage::Create();
      img = itmp ? 1 : -1;
      if (itmp) {
         delete itmp;
         itmp=NULL;
      }
      gErrorIgnoreLevel = sav;
   }
   if (img > 0) {
      fFileSaveMenu->AddEntry(Form("%s.&jpg",name),  kFileSaveAsJPG);
      fFileSaveMenu->AddEntry(Form("%s.&png",name),  kFileSaveAsPNG);
   }

   fFileSaveMenu->AddEntry(Form("%s.&C",   name), kFileSaveAsC);
   fFileSaveMenu->AddEntry(Form("%s.&root",name), kFileSaveAsRoot);

   fFileMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fFileMenu->AddEntry("&New Canvas",   kFileNewCanvas);
   fFileMenu->AddEntry("&Open...",      kFileOpen);
   fFileMenu->AddEntry("&Close Canvas", kFileCloseCanvas);
   fFileMenu->AddSeparator();
   fFileMenu->AddPopup("&Save",         fFileSaveMenu);
   fFileMenu->AddEntry("Save &As...",   kFileSaveAs);
   fFileMenu->AddSeparator();
   fFileMenu->AddEntry("&Print...",     kFilePrint);
   fFileMenu->AddSeparator();
   fFileMenu->AddEntry("&Quit ROOT",    kFileQuit);

   fEditClearMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fEditClearMenu->AddEntry("&Pad",     kEditClearPad);
   fEditClearMenu->AddEntry("&Canvas",  kEditClearCanvas);

   fEditMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fEditMenu->AddEntry("&Style...",     kEditStyle);
   fEditMenu->AddSeparator();
   fEditMenu->AddEntry("Cu&t",          kEditCut);
   fEditMenu->AddEntry("&Copy",         kEditCopy);
   fEditMenu->AddEntry("&Paste",        kEditPaste);
   fEditMenu->AddSeparator();
   fEditMenu->AddPopup("C&lear",        fEditClearMenu);
   fEditMenu->AddSeparator();
   fEditMenu->AddEntry("&Undo",         kEditUndo);
   fEditMenu->AddEntry("&Redo",         kEditRedo);

   fEditMenu->DisableEntry(kEditCut);
   fEditMenu->DisableEntry(kEditCopy);
   fEditMenu->DisableEntry(kEditPaste);
   fEditMenu->DisableEntry(kEditUndo);
   fEditMenu->DisableEntry(kEditRedo);

   fViewWithMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fViewWithMenu->AddEntry("&X3D",      kViewX3D);
   fViewWithMenu->AddEntry("&OpenGL",   kViewOpenGL);

   fViewMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fViewMenu->AddEntry("&Editor",       kViewEditor);
   fViewMenu->AddEntry("&Toolbar",      kViewToolbar);
   fViewMenu->AddEntry("Event &Statusbar", kViewEventStatus);
   fViewMenu->AddEntry("T&oolTip Info", kViewToolTips);
   fViewMenu->AddSeparator();
   fViewMenu->AddEntry("&Colors",       kViewColors);
   fViewMenu->AddEntry("&Fonts",        kViewFonts);
   fViewMenu->AddEntry("&Markers",      kViewMarkers);
   fViewMenu->AddSeparator();
   fViewMenu->AddEntry("&Iconify",      kViewIconify);
   fViewMenu->AddSeparator();
   fViewMenu->AddPopup("&View With",    fViewWithMenu);

   fViewMenu->DisableEntry(kViewFonts);

   fOptionMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fOptionMenu->AddEntry("&Auto Resize Canvas",  kOptionAutoResize);
   fOptionMenu->AddEntry("&Resize Canvas",       kOptionResizeCanvas);
   fOptionMenu->AddEntry("&Move Opaque",         kOptionMoveOpaque);
   fOptionMenu->AddEntry("Resize &Opaque",       kOptionResizeOpaque);
   fOptionMenu->AddSeparator();
   fOptionMenu->AddEntry("&Interrupt",           kOptionInterrupt);
   fOptionMenu->AddEntry("R&efresh",             kOptionRefresh);
   fOptionMenu->AddSeparator();
   fOptionMenu->AddEntry("&Pad Auto Exec",       kOptionAutoExec);
   fOptionMenu->AddSeparator();
   fOptionMenu->AddEntry("&Statistics",          kOptionStatistics);
   fOptionMenu->AddEntry("Histogram &Title",     kOptionHistTitle);
   fOptionMenu->AddEntry("&Fit Parameters",      kOptionFitParams);
   fOptionMenu->AddEntry("Can Edit &Histograms", kOptionCanEdit);

   // Opaque options initialized in InitWindow()
   fOptionMenu->CheckEntry(kOptionAutoResize);
   if (gStyle->GetOptStat())
      fOptionMenu->CheckEntry(kOptionStatistics);
   if (gStyle->GetOptTitle())
      fOptionMenu->CheckEntry(kOptionHistTitle);
   if (gStyle->GetOptFit())
      fOptionMenu->CheckEntry(kOptionFitParams);
   if (gROOT->GetEditHistograms())
      fOptionMenu->CheckEntry(kOptionCanEdit);

   fToolsMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fToolsMenu->AddEntry("&Inspect ROOT",   kInspectRoot);
   fToolsMenu->AddEntry("&Class Tree",     kClassesTree);
   fToolsMenu->AddEntry("&Fit Panel",      kFitPanel);
   fToolsMenu->AddEntry("&Start Browser",  kToolsBrowser);
   fToolsMenu->AddEntry("&Gui Builder",    kToolsBuilder);
   fToolsMenu->AddEntry("&Event Recorder", kToolsRecorder);

   fHelpMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fHelpMenu->AddLabel("Basic Help On...");
   fHelpMenu->AddSeparator();
   fHelpMenu->AddEntry("&Canvas",          kHelpOnCanvas);
   fHelpMenu->AddEntry("&Menus",           kHelpOnMenus);
   fHelpMenu->AddEntry("&Graphics Editor", kHelpOnGraphicsEd);
   fHelpMenu->AddEntry("&Browser",         kHelpOnBrowser);
   fHelpMenu->AddEntry("&Objects",         kHelpOnObjects);
   fHelpMenu->AddEntry("&PostScript",      kHelpOnPS);
   fHelpMenu->AddSeparator();
   fHelpMenu->AddEntry("&About ROOT...",   kHelpAbout);

   // This main frame will process the menu commands
   fFileMenu->Associate(this);
   fFileSaveMenu->Associate(this);
   fEditMenu->Associate(this);
   fEditClearMenu->Associate(this);
   fViewMenu->Associate(this);
   fViewWithMenu->Associate(this);
   fOptionMenu->Associate(this);
   fToolsMenu->Associate(this);
   fHelpMenu->Associate(this);

   // Create menubar layout hints
   fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 1, 1);
   fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
   fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

   // Create menubar
   fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);
   fMenuBar->AddPopup("&File",    fFileMenu,    fMenuBarItemLayout);
   fMenuBar->AddPopup("&Edit",    fEditMenu,    fMenuBarItemLayout);
   fMenuBar->AddPopup("&View",    fViewMenu,    fMenuBarItemLayout);
   fMenuBar->AddPopup("&Options", fOptionMenu,  fMenuBarItemLayout);
   fMenuBar->AddPopup("&Tools",   fToolsMenu,   fMenuBarItemLayout);
   fMenuBar->AddPopup("&Help",    fHelpMenu,    fMenuBarHelpLayout);

   AddFrame(fMenuBar, fMenuBarLayout);

   fHorizontal1 = new TGHorizontal3DLine(this);
   fHorizontal1Layout = new TGLayoutHints(kLHintsTop | kLHintsExpandX);
   AddFrame(fHorizontal1, fHorizontal1Layout);

   // Create toolbar dock
   fToolDock = new TGDockableFrame(this);
   fToolDock->SetCleanup();
   fToolDock->EnableHide(kFALSE);
   AddFrame(fToolDock, fDockLayout = new TGLayoutHints(kLHintsExpandX));

   // will alocate it later
   fToolBar = 0;
   fVertical1 = 0;
   fVertical2 = 0;
   fVertical1Layout = 0;
   fVertical2Layout = 0;

   fToolBarSep = new TGHorizontal3DLine(this);
   fToolBarLayout = new TGLayoutHints(kLHintsTop |  kLHintsExpandX);
   AddFrame(fToolBarSep, fToolBarLayout);

   fMainFrame = new TGCompositeFrame(this, GetWidth() + 4, GetHeight() + 4,
                                      kHorizontalFrame);
   fMainFrameLayout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);

   // Create editor frame that will host the pad editor
   fEditorFrame = new TGCompositeFrame(fMainFrame, 175, fMainFrame->GetHeight()+4, kFixedWidth);
   fEditorLayout = new TGLayoutHints(kLHintsExpandY | kLHintsLeft);
   fMainFrame->AddFrame(fEditorFrame, fEditorLayout);

   // Create canvas and canvas container that will host the ROOT graphics
   fCanvasWindow = new TGCanvas(fMainFrame, GetWidth()+4, GetHeight()+4,
                                kSunkenFrame | kDoubleBorder);

   fCanvasID = -1;

   if (fCanvas->UseGL()) {
      fCanvas->SetSupportGL(kFALSE);
      //first, initialize GL (if not yet)
      if (!gGLManager) {
         TString x = "win32";
         if (gVirtualX->InheritsFrom("TGX11"))
            x = "x11";
         else if (gVirtualX->InheritsFrom("TGCocoa"))
            x = "osx";

         TPluginHandler *ph = gROOT->GetPluginManager()->FindHandler("TGLManager", x);

         if (ph && ph->LoadPlugin() != -1) {
            if (!ph->ExecPlugin(0))
               Error("CreateCanvas", "GL manager plugin failed");
         }
      }

      if (gGLManager) {
         fCanvasID = gGLManager->InitGLWindow((ULong_t)fCanvasWindow->GetViewPort()->GetId());
         if (fCanvasID != -1) {
            //Create gl context.
            const Int_t glCtx = gGLManager->CreateGLContext(fCanvasID);
            if (glCtx != -1) {
               fCanvas->SetSupportGL(kTRUE);
               fCanvas->SetGLDevice(glCtx);//Now, fCanvas is responsible for context deletion!
            } else
               Error("CreateCanvas", "GL context creation failed.");
         } else
            Error("CreateCanvas", "GL window creation failed\n");
      }
   }

   if (fCanvasID == -1)
      fCanvasID = gVirtualX->InitWindow((ULong_t)fCanvasWindow->GetViewPort()->GetId());

   Window_t win = gVirtualX->GetWindowID(fCanvasID);
   fCanvasContainer = new TRootContainer(this, win, fCanvasWindow->GetViewPort());
   fCanvasWindow->SetContainer(fCanvasContainer);
   fCanvasLayout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY | kLHintsRight);

   fMainFrame->AddFrame(fCanvasWindow, fCanvasLayout);
   AddFrame(fMainFrame, fMainFrameLayout);

   // create the tooltip with a timeout of 250 ms
   fToolTip = new TGToolTip(fClient->GetDefaultRoot(), fCanvasWindow, "", 250);

   fCanvas->Connect("ProcessedEvent(Int_t, Int_t, Int_t, TObject*)",
                    "TRootCanvas", this,
                    "EventInfo(Int_t, Int_t, Int_t, TObject*)");

   // Create status bar
   int parts[] = { 33, 10, 10, 47 };
   fStatusBar = new TGStatusBar(this, 10, 10);
   fStatusBar->SetParts(parts, 4);

   fStatusBarLayout = new TGLayoutHints(kLHintsBottom | kLHintsLeft | kLHintsExpandX, 2, 2, 1, 1);

   AddFrame(fStatusBar, fStatusBarLayout);

   // Misc
   SetWindowName(name);
   SetIconName(name);
   fIconPic = SetIconPixmap("macro_s.xpm");
   SetClassHints("ROOT", "Canvas");

   SetEditDisabled(kEditDisable);
   MapSubwindows();

   // by default status bar, tool bar and pad editor are hidden
   HideFrame(fStatusBar);
   HideFrame(fToolDock);
   HideFrame(fToolBarSep);
   HideFrame(fHorizontal1);

   ShowToolBar(kFALSE);
   ShowEditor(kFALSE);

   // we need to use GetDefaultSize() to initialize the layout algorithm...
   Resize(GetDefaultSize());

   gVirtualX->SetDNDAware(fId, fDNDTypeList);
   SetDNDTarget(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete ROOT basic canvas. Order is significant. Delete in reverse
/// order of creation.

TRootCanvas::~TRootCanvas()
{
   delete fToolTip;
   if (fIconPic) gClient->FreePicture(fIconPic);
   if (fEditor && !fEmbedded) delete fEditor;
   if (fToolBar) {
      Disconnect(fToolDock, "Docked()",   this, "AdjustSize()");
      Disconnect(fToolDock, "Undocked()", this, "AdjustSize()");
      fToolBar->Cleanup();
      delete fToolBar;
   }

   if (!MustCleanup()) {
      delete fStatusBar;
      delete fStatusBarLayout;
      delete fCanvasContainer;
      delete fCanvasWindow;

      delete fEditorFrame;
      delete fEditorLayout;
      delete fMainFrame;
      delete fMainFrameLayout;
      delete fToolBarSep;
      delete fToolDock;
      delete fToolBarLayout;
      delete fHorizontal1;
      delete fHorizontal1Layout;

      delete fMenuBar;
      delete fMenuBarLayout;
      delete fMenuBarItemLayout;
      delete fMenuBarHelpLayout;
      delete fCanvasLayout;
      delete fDockLayout;
   }

   delete fFileMenu;
   delete fFileSaveMenu;
   delete fEditMenu;
   delete fEditClearMenu;
   delete fViewMenu;
   delete fViewWithMenu;
   delete fOptionMenu;
   delete fToolsMenu;
   delete fHelpMenu;
}

////////////////////////////////////////////////////////////////////////////////
/// Called via TCanvasImp interface by TCanvas.

void TRootCanvas::Close()
{
   TVirtualPadEditor* gged = TVirtualPadEditor::GetPadEditor(kFALSE);
   if(gged && gged->GetCanvas() == fCanvas) {
      if (fEmbedded) {
         ((TGedEditor *)gged)->SetModel(0, 0, kButton1Down);
         ((TGedEditor *)gged)->SetCanvas(0);
      }
      else gged->Hide();
   }

   gVirtualX->CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Really delete the canvas and this GUI.

void TRootCanvas::ReallyDelete()
{
   TVirtualPadEditor* gged = TVirtualPadEditor::GetPadEditor(kFALSE);
   if(gged && gged->GetCanvas() == fCanvas) {
      if (fEmbedded) {
         ((TGedEditor *)gged)->SetModel(0, 0, kButton1Down);
         ((TGedEditor *)gged)->SetCanvas(0);
      }
      else gged->Hide();
   }

   fToolTip->Hide();
   Disconnect(fCanvas, "ProcessedEvent(Int_t, Int_t, Int_t, TObject*)",
              this, "EventInfo(Int_t, Int_t, Int_t, TObject*)");

   fCanvas->SetCanvasImp(0);
   fCanvas->Clear();
   fCanvas->SetName("");
   if (gPad && gPad->GetCanvas() == fCanvas)
      gPad = 0;
   delete this;
}

////////////////////////////////////////////////////////////////////////////////
/// In case window is closed via WM we get here.

void TRootCanvas::CloseWindow()
{
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Return width of canvas container.

UInt_t TRootCanvas::GetCwidth() const
{
   return fCanvasContainer->GetWidth();
}

////////////////////////////////////////////////////////////////////////////////
/// Return height of canvas container.

UInt_t TRootCanvas::GetCheight() const
{
   return fCanvasContainer->GetHeight();
}

////////////////////////////////////////////////////////////////////////////////
/// Gets the size and position of the window containing the canvas. This
/// size includes the menubar and borders.

UInt_t TRootCanvas::GetWindowGeometry(Int_t &x, Int_t &y, UInt_t &w, UInt_t &h)
{
   gVirtualX->GetWindowSize(fId, x, y, w, h);

   Window_t childdum;
   gVirtualX->TranslateCoordinates(fId, gClient->GetDefaultRoot()->GetId(),
                                   0, 0, x, y, childdum);
   if (!fCanvas->GetShowEditor()) return 0;
   return fEditorFrame->GetWidth();
}

////////////////////////////////////////////////////////////////////////////////
/// Set text in status bar.

void TRootCanvas::SetStatusText(const char *txt, Int_t partidx)
{
   fStatusBar->SetText(txt, partidx);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle menu and other command generated by the user.

Bool_t TRootCanvas::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   TRootHelpDialog *hd;
   TList *lc;

   switch (GET_MSG(msg)) {

      case kC_COMMAND:

         switch (GET_SUBMSG(msg)) {

            case kCM_BUTTON:
            case kCM_MENU:

               switch (parm1) {
                  // Handle toolbar items...
                  case kToolModify:
                     gROOT->SetEditorMode();
                     break;
                  case kToolArc:
                     gROOT->SetEditorMode("Arc");
                     break;
                  case kToolLine:
                     gROOT->SetEditorMode("Line");
                     break;
                  case kToolArrow:
                     gROOT->SetEditorMode("Arrow");
                     break;
                  case kToolDiamond:
                     gROOT->SetEditorMode("Diamond");
                     break;
                  case kToolEllipse:
                     gROOT->SetEditorMode("Ellipse");
                     break;
                  case kToolPad:
                     gROOT->SetEditorMode("Pad");
                     break;
                  case kToolPave:
                     gROOT->SetEditorMode("Pave");
                     break;
                  case kToolPLabel:
                     gROOT->SetEditorMode("PaveLabel");
                     break;
                  case kToolPText:
                     gROOT->SetEditorMode("PaveText");
                     break;
                  case kToolPsText:
                     gROOT->SetEditorMode("PavesText");
                     break;
                  case kToolGraph:
                     gROOT->SetEditorMode("PolyLine");
                     break;
                  case kToolCurlyLine:
                     gROOT->SetEditorMode("CurlyLine");
                     break;
                  case kToolCurlyArc:
                     gROOT->SetEditorMode("CurlyArc");
                     break;
                  case kToolLatex:
                     gROOT->SetEditorMode("Text");
                     break;
                  case kToolMarker:
                     gROOT->SetEditorMode("Marker");
                     break;
                  case kToolCutG:
                     gROOT->SetEditorMode("CutG");
                     break;

                  // Handle File menu items...
                  case kFileNewCanvas:
                     gROOT->MakeDefCanvas();
                     break;
                  case kFileOpen:
                     {
                        static TString dir(".");
                        TGFileInfo fi;
                        fi.fFileTypes = gOpenTypes;
                        fi.SetIniDir(dir);
                        new TGFileDialog(fClient->GetDefaultRoot(), this, kFDOpen,&fi);
                        if (!fi.fFilename) return kTRUE;
                        dir = fi.fIniDir;
                        new TFile(fi.fFilename, "update");
                     }
                     break;
                  case kFileSaveAs:
                     {
                        TString workdir = gSystem->WorkingDirectory();
                        static TString dir(".");
                        static Int_t typeidx = 0;
                        static Bool_t overwr = kFALSE;
                        TGFileInfo fi;
                        TString defaultType = gEnv->GetValue("Canvas.SaveAsDefaultType", ".pdf");
                        if (typeidx == 0) {
                           for (int i=1;gSaveAsTypes[i];i+=2) {
                              TString ftype = gSaveAsTypes[i];
                              if (ftype.EndsWith(defaultType.Data())) {
                                 typeidx = i-1;
                                 break;
                              }
                           }
                        }
                        fi.fFileTypes   = gSaveAsTypes;
                        fi.SetIniDir(dir);
                        fi.fFileTypeIdx = typeidx;
                        fi.fOverwrite = overwr;
                        new TGFileDialog(fClient->GetDefaultRoot(), this, kFDSave, &fi);
                        gSystem->ChangeDirectory(workdir.Data());
                        if (!fi.fFilename) return kTRUE;
                        Bool_t  appendedType = kFALSE;
                        TString fn = fi.fFilename;
                        TString ft = fi.fFileTypes[fi.fFileTypeIdx+1];
                        dir     = fi.fIniDir;
                        typeidx = fi.fFileTypeIdx;
                        overwr  = fi.fOverwrite;
again:
                        if (fn.EndsWith(".root") ||
                            fn.EndsWith(".ps")   ||
                            fn.EndsWith(".eps")  ||
                            fn.EndsWith(".pdf")  ||
                            fn.EndsWith(".svg")  ||
                            fn.EndsWith(".tex")  ||
                            fn.EndsWith(".gif")  ||
                            fn.EndsWith(".xml")  ||
                            fn.EndsWith(".xpm")  ||
                            fn.EndsWith(".jpg")  ||
                            fn.EndsWith(".png")  ||
                            fn.EndsWith(".xcf")  ||
                            fn.EndsWith(".tiff")) {
                           fCanvas->SaveAs(fn);
                        } else if (fn.EndsWith(".C"))
                           fCanvas->SaveSource(fn);
                        else {
                           if (!appendedType) {
                              if (ft.Index(".") != kNPOS) {
                                 fn += ft(ft.Index("."), ft.Length());
                                 appendedType = kTRUE;
                                 goto again;
                              }
                           }
                           Warning("ProcessMessage", "file %s cannot be saved with this extension", fi.fFilename);
                        }
                        for (int i=1;gSaveAsTypes[i];i+=2) {
                           TString ftype = gSaveAsTypes[i];
                           ftype.ReplaceAll("*.", ".");
                           if (fn.EndsWith(ftype.Data())) {
                              typeidx = i-1;
                              break;
                           }
                        }
                     }
                     break;
                  case kFileSaveAsRoot:
                     fCanvas->SaveAs(".root");
                     break;
                  case kFileSaveAsC:
                     fCanvas->SaveSource();
                     break;
                  case kFileSaveAsPS:
                     fCanvas->SaveAs();
                     break;
                  case kFileSaveAsEPS:
                     fCanvas->SaveAs(".eps");
                     break;
                  case kFileSaveAsPDF:
                     fCanvas->SaveAs(".pdf");
                     break;
                  case kFileSaveAsGIF:
                     fCanvas->SaveAs(".gif");
                     break;
                  case kFileSaveAsJPG:
                     fCanvas->SaveAs(".jpg");
                     break;
                  case kFileSaveAsPNG:
                     fCanvas->SaveAs(".png");
                     break;
                  case kFileSaveAsTEX:
                     fCanvas->SaveAs(".tex");
                     break;
                  case kFilePrint:
                     PrintCanvas();
                     break;
                  case kFileCloseCanvas:
                     SendCloseMessage();
                     break;
                  case kFileQuit:
                     if (!gApplication->ReturnFromRun()) {
                        if ((TVirtualPadEditor::GetPadEditor(kFALSE) != 0))
                           TVirtualPadEditor::Terminate();
                        SendCloseMessage();
                     }
                     if (TVirtualPadEditor::GetPadEditor(kFALSE) != 0)
                        TVirtualPadEditor::Terminate();
                     if (TClass::GetClass("TStyleManager"))
                        gROOT->ProcessLine("TStyleManager::Terminate()");
                     gApplication->Terminate(0);
                     break;

                  // Handle Edit menu items...
                  case kEditStyle:
                     if (!TClass::GetClass("TStyleManager"))
                        gSystem->Load("libGed");
                     gROOT->ProcessLine("TStyleManager::Show()");
                     break;
                  case kEditCut:
                     // still noop
                     break;
                  case kEditCopy:
                     // still noop
                     break;
                  case kEditPaste:
                     // still noop
                     break;
                  case kEditUndo:
                     // noop
                     break;
                  case kEditRedo:
                     // noop
                     break;
                  case kEditClearPad:
                     gPad->Clear();
                     gPad->Modified();
                     gPad->Update();
                     break;
                  case kEditClearCanvas:
                     fCanvas->Clear();
                     fCanvas->Modified();
                     fCanvas->Update();
                     break;

                  // Handle View menu items...
                  case kViewEditor:
                     fCanvas->ToggleEditor();
                     break;
                  case kViewToolbar:
                     fCanvas->ToggleToolBar();
                     break;
                  case kViewEventStatus:
                     fCanvas->ToggleEventStatus();
                     break;
                  case kViewToolTips:
                     fCanvas->ToggleToolTips();
                     break;
                  case kViewColors:
                     {
                        TVirtualPad *padsav = gPad->GetCanvas();
                        //This was the code with the old color table
                        //   TCanvas *m = new TCanvas("colors","Color Table");
                        //   TPad::DrawColorTable();
                        //   m->Update();
                        TColorWheel *wheel = new TColorWheel();
                        wheel->Draw();

                        //tp: with Cocoa, window is visible (and repainted)
                        //before wheel->Draw() was called and you can see "empty"
                        //canvas.
                        gPad->Update();
                        //
                        if (padsav) padsav->cd();
                     }
                     break;
                  case kViewFonts:
                     // noop
                     break;
                  case kViewMarkers:
                     {
                        TVirtualPad *padsav = gPad->GetCanvas();
                        TCanvas *m = new TCanvas("markers","Marker Types",600,200);
                        TMarker::DisplayMarkerTypes();
                        m->Update();
                        if (padsav) padsav->cd();
                     }
                     break;
                  case kViewIconify:
                     Iconify();
                     break;
                  case kViewX3D:
                     gPad->GetViewer3D("x3d");
                     break;
                  case kViewOpenGL:
                     gPad->GetViewer3D("ogl");
                     break;

                  // Handle Option menu items...
                  case kOptionAutoExec:
                     fCanvas->ToggleAutoExec();
                     if (fCanvas->GetAutoExec()) {
                        fOptionMenu->CheckEntry(kOptionAutoExec);
                     } else {
                        fOptionMenu->UnCheckEntry(kOptionAutoExec);
                     }
                     break;
                  case kOptionAutoResize:
                     {
                        fAutoFit = fAutoFit ? kFALSE : kTRUE;
                        int opt = fCanvasContainer->GetOptions();
                        if (fAutoFit) {
                           opt &= ~kFixedSize;
                           fOptionMenu->CheckEntry(kOptionAutoResize);
                        } else {
                           opt |= kFixedSize;
                           fOptionMenu->UnCheckEntry(kOptionAutoResize);
                        }
                        fCanvasContainer->ChangeOptions(opt);
                        // in case of autofit this will generate a configure
                        // event for the container and this will force the
                        // update of the TCanvas
                        //Layout();
                     }
                     Layout();
                     break;
                  case kOptionResizeCanvas:
                     FitCanvas();
                     break;
                  case kOptionMoveOpaque:
                     if (fCanvas->OpaqueMoving()) {
                        fCanvas->MoveOpaque(0);
                        fOptionMenu->UnCheckEntry(kOptionMoveOpaque);
                     } else {
                        fCanvas->MoveOpaque(1);
                        fOptionMenu->CheckEntry(kOptionMoveOpaque);
                     }
                     break;
                  case kOptionResizeOpaque:
                     if (fCanvas->OpaqueResizing()) {
                        fCanvas->ResizeOpaque(0);
                        fOptionMenu->UnCheckEntry(kOptionResizeOpaque);
                     } else {
                        fCanvas->ResizeOpaque(1);
                        fOptionMenu->CheckEntry(kOptionResizeOpaque);
                     }
                     break;
                  case kOptionInterrupt:
                     gROOT->SetInterrupt();
                     break;
                  case kOptionRefresh:
                     fCanvas->Paint();
                     fCanvas->Update();
                     break;
                  case kOptionStatistics:
                     if (gStyle->GetOptStat()) {
                        gStyle->SetOptStat(0);
                        delete gPad->FindObject("stats");
                        fOptionMenu->UnCheckEntry(kOptionStatistics);
                     } else {
                        gStyle->SetOptStat(1);
                        fOptionMenu->CheckEntry(kOptionStatistics);
                     }
                     gPad->Modified();
                     fCanvas->Update();
                     break;
                  case kOptionHistTitle:
                     if (gStyle->GetOptTitle()) {
                        gStyle->SetOptTitle(0);
                        delete gPad->FindObject("title");
                        fOptionMenu->UnCheckEntry(kOptionHistTitle);
                     } else {
                        gStyle->SetOptTitle(1);
                        fOptionMenu->CheckEntry(kOptionHistTitle);
                     }
                     gPad->Modified();
                     fCanvas->Update();
                     break;
                  case kOptionFitParams:
                     if (gStyle->GetOptFit()) {
                        gStyle->SetOptFit(0);
                        fOptionMenu->UnCheckEntry(kOptionFitParams);
                     } else {
                        gStyle->SetOptFit(1);
                        fOptionMenu->CheckEntry(kOptionFitParams);
                     }
                     gPad->Modified();
                     fCanvas->Update();
                     break;
                  case kOptionCanEdit:
                     if (gROOT->GetEditHistograms()) {
                        gROOT->SetEditHistograms(kFALSE);
                        fOptionMenu->UnCheckEntry(kOptionCanEdit);
                     } else {
                        gROOT->SetEditHistograms(kTRUE);
                        fOptionMenu->CheckEntry(kOptionCanEdit);
                     }
                     break;

                  // Handle Tools menu items...
                  case kInspectRoot:
                     fCanvas->cd();
                     gROOT->Inspect();
                     fCanvas->Update();
                     break;
                  case kToolsBrowser:
                     new TBrowser("browser");
                     break;
                  case kToolsBuilder:
                     TGuiBuilder::Instance();
                     break;
                  case kToolsRecorder:
                     gROOT->ProcessLine("new TGRecorder()");
                     break;

                  // Handle Tools menu items...
                  case kClassesTree:
                     {
                        TString cdef;
                        lc = (TList*)gROOT->GetListOfCanvases();
                        if (lc->FindObject("ClassTree")) {
                           cdef = TString::Format("ClassTree_%d", lc->GetSize()+1);
                        } else {
                           cdef = "ClassTree";
                        }
                        new TClassTree(cdef.Data(), "TObject");
                        fCanvas->Update();
                     }
                     break;

               case kFitPanel:
                     {
                        // use plugin manager to create instance of TFitEditor
                        TPluginHandler *handler = gROOT->GetPluginManager()->FindHandler("TFitEditor");
                        if (handler && handler->LoadPlugin() != -1) {
                           if (handler->ExecPlugin(2, fCanvas, 0) == 0)
                              Error("FitPanel", "Unable to crate the FitPanel");
                        }
                        else
                           Error("FitPanel", "Unable to find the FitPanel plug-in");
                     }
                     break;

                  // Handle Help menu items...
                  case kHelpAbout:
                     {
#ifdef R__UNIX
                        TString rootx = TROOT::GetBinDir() + "/root -a &";
                        gSystem->Exec(rootx);
#else
#ifdef WIN32
                        new TWin32SplashThread(kTRUE);
#else

                        char str[32];
                        sprintf(str, "About ROOT %s...", gROOT->GetVersion());
                        hd = new TRootHelpDialog(this, str, 600, 400);
                        hd->SetText(gHelpAbout);
                        hd->Popup();
#endif
#endif
                     }
                     break;
                  case kHelpOnCanvas:
                     hd = new TRootHelpDialog(this, "Help on Canvas...", 600, 400);
                     hd->SetText(gHelpCanvas);
                     hd->Popup();
                     break;
                  case kHelpOnMenus:
                     hd = new TRootHelpDialog(this, "Help on Menus...", 600, 400);
                     hd->SetText(gHelpPullDownMenus);
                     hd->Popup();
                     break;
                  case kHelpOnGraphicsEd:
                     hd = new TRootHelpDialog(this, "Help on Graphics Editor...", 600, 400);
                     hd->SetText(gHelpGraphicsEditor);
                     hd->Popup();
                     break;
                  case kHelpOnBrowser:
                     hd = new TRootHelpDialog(this, "Help on Browser...", 600, 400);
                     hd->SetText(gHelpBrowser);
                     hd->Popup();
                     break;
                  case kHelpOnObjects:
                     hd = new TRootHelpDialog(this, "Help on Objects...", 600, 400);
                     hd->SetText(gHelpObjects);
                     hd->Popup();
                     break;
                  case kHelpOnPS:
                     hd = new TRootHelpDialog(this, "Help on PostScript...", 600, 400);
                     hd->SetText(gHelpPostscript);
                     hd->Popup();
                     break;
               }
            default:
               break;
         }
      default:
         break;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Called by TCanvas ctor to get window indetifier.

Int_t TRootCanvas::InitWindow()
{
   if (fCanvas->OpaqueMoving())
      fOptionMenu->CheckEntry(kOptionMoveOpaque);
   if (fCanvas->OpaqueResizing())
      fOptionMenu->CheckEntry(kOptionResizeOpaque);

   return fCanvasID;
}

////////////////////////////////////////////////////////////////////////////////
/// Set size of canvas container. Units in pixels.

void TRootCanvas::SetCanvasSize(UInt_t w, UInt_t h)
{
   // turn off autofit, we want to stay at the given size
   fAutoFit = kFALSE;
   fOptionMenu->UnCheckEntry(kOptionAutoResize);
   int opt = fCanvasContainer->GetOptions();
   opt |= kFixedSize;    // turn on fixed size mode
   fCanvasContainer->ChangeOptions(opt);
   fCanvasContainer->SetWidth(w);
   fCanvasContainer->SetHeight(h);
   Layout();  // force layout (will update container to given size)
   fCanvas->Resize();
   fCanvas->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Set canvas position (units in pixels).

void TRootCanvas::SetWindowPosition(Int_t x, Int_t y)
{
   Move(x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Set size of canvas (units in pixels).

void TRootCanvas::SetWindowSize(UInt_t w, UInt_t h)
{
   Resize(w, h);

   // Make sure the change of size is really done.
   gVirtualX->Update(1);
   if (!gThreadXAR) {
      gSystem->Sleep(100);
      gSystem->ProcessEvents();
      gSystem->Sleep(10);
      gSystem->ProcessEvents();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Put canvas window on top of the window stack.

void TRootCanvas::RaiseWindow()
{
   gVirtualX->RaiseWindow(GetId());
}

////////////////////////////////////////////////////////////////////////////////
/// Change title on window.

void TRootCanvas::SetWindowTitle(const char *title)
{
   SetWindowName(title);
   SetIconName(title);
   fToolDock->SetWindowName(Form("ToolBar: %s", title));
}

////////////////////////////////////////////////////////////////////////////////
/// Fit canvas container to current window size.

void TRootCanvas::FitCanvas()
{
   if (!fAutoFit) {
      int opt = fCanvasContainer->GetOptions();
      int oopt = opt;
      opt &= ~kFixedSize;   // turn off fixed size mode
      fCanvasContainer->ChangeOptions(opt);
      Layout();  // force layout
      fCanvas->Resize();
      fCanvas->Update();
      fCanvasContainer->ChangeOptions(oopt);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print the canvas.

void TRootCanvas::PrintCanvas()
{
   Int_t ret = 0;
   Bool_t pname = kTRUE;
   char *printer, *printCmd;
   static TString sprinter, sprintCmd;

   if (sprinter == "")
      printer = StrDup(gEnv->GetValue("Print.Printer", ""));
   else
      printer = StrDup(sprinter);
   if (sprintCmd == "")
#ifndef WIN32
      printCmd = StrDup(gEnv->GetValue("Print.Command", ""));
#else
      printCmd = StrDup(gEnv->GetValue("Print.Command", "start AcroRd32.exe /p"));
#endif
   else
      printCmd = StrDup(sprintCmd);

   new TGPrintDialog(fClient->GetDefaultRoot(), this, 400, 150,
                     &printer, &printCmd, &ret);
   if (ret) {
      sprinter  = printer;
      sprintCmd = printCmd;

      if (sprinter == "")
         pname = kFALSE;

      TString fn = "rootprint";
      FILE *f = gSystem->TempFileName(fn, gEnv->GetValue("Print.Directory", gSystem->TempDirectory()));
      if (f) fclose(f);
      fn += TString::Format(".%s",gEnv->GetValue("Print.FileType", "pdf"));
      fCanvas->Print(fn);

      TString cmd = sprintCmd;
      if (cmd.Contains("%p"))
         cmd.ReplaceAll("%p", sprinter);
      else if (pname) {
         cmd += " "; cmd += sprinter; cmd += " ";
      }

      if (cmd.Contains("%f"))
         cmd.ReplaceAll("%f", fn);
      else {
         cmd += " "; cmd += fn; cmd += " ";
      }

      gSystem->Exec(cmd);
#ifndef WIN32
      gSystem->Unlink(fn);
#endif
   }
   delete [] printer;
   delete [] printCmd;
}

////////////////////////////////////////////////////////////////////////////////
/// Display a tooltip with infos about the primitive below the cursor.

void TRootCanvas::EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected)
{
   fToolTip->Hide();
   if (!fCanvas->GetShowToolTips() || selected == 0 ||
       event != kMouseMotion || fButton != 0)
      return;
   TString tipInfo;
   TString objInfo = selected->GetObjectInfo(px, py);
   if (objInfo.BeginsWith("-")) {
      // if the string begins with '-', display only the object info
      objInfo.Remove(TString::kLeading, '-');
      tipInfo = objInfo;
   }
   else {
      const char *title = selected->GetTitle();
      tipInfo += TString::Format("%s::%s", selected->ClassName(),
                                 selected->GetName());
      if (title && strlen(title))
         tipInfo += TString::Format("\n%s", selected->GetTitle());
      tipInfo += TString::Format("\n%d, %d", px, py);
      if (!objInfo.IsNull())
         tipInfo += TString::Format("\n%s", objInfo.Data());
   }
   fToolTip->SetText(tipInfo.Data());
   fToolTip->SetPosition(px+15, py+15);
   fToolTip->Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Show or hide menubar.

void TRootCanvas::ShowMenuBar(Bool_t show)
{
   if (show)  ShowFrame(fMenuBar);
   else       HideFrame(fMenuBar);
}

////////////////////////////////////////////////////////////////////////////////
/// Show or hide statusbar.

void TRootCanvas::ShowStatusBar(Bool_t show)
{
   UInt_t dh = fClient->GetDisplayHeight();
   UInt_t ch = fCanvas->GetWindowHeight();

   UInt_t h = GetHeight();
   UInt_t sh = fStatusBar->GetHeight()+2;

   if (show) {
      ShowFrame(fStatusBar);
      fViewMenu->CheckEntry(kViewEventStatus);
      if (dh - ch >= sh) h = h + sh;
      else h = ch;
   } else {
      HideFrame(fStatusBar);
      fViewMenu->UnCheckEntry(kViewEventStatus);
      if (dh - ch < sh) h = ch;
      else h = h - sh;
   }
   Resize(GetWidth(), h);
}

////////////////////////////////////////////////////////////////////////////////
/// Show or hide side frame.

void TRootCanvas::ShowEditor(Bool_t show)
{
   TVirtualPad *savedPad = 0;
   savedPad = (TVirtualPad *) gPad;
   gPad = Canvas();

   UInt_t w = GetWidth();
   UInt_t e = fEditorFrame->GetWidth();
   UInt_t h = GetHeight();
   UInt_t s = fHorizontal1->GetHeight();

   if (fParent && fParent != fClient->GetDefaultRoot()) {
      TGMainFrame *main = (TGMainFrame *)fParent->GetMainFrame();
      fMainFrame->HideFrame(fEditorFrame);
      if (main && main->InheritsFrom("TRootBrowser")) {
         TRootBrowser *browser = (TRootBrowser *)main;
         if (!fEmbedded)
            browser->GetTabRight()->Connect("Selected(Int_t)", "TRootCanvas",
                                            this, "Activated(Int_t)");
         fEmbedded = kTRUE;
         if (show && (!fEditor || !((TGedEditor *)fEditor)->IsMapped())) {
            if (!browser->GetTabLeft()->GetTabTab("Pad Editor")) {
               if (browser->GetActFrame()) { //already in edit mode
                  TTimer::SingleShot(200, "TRootCanvas", this, "ShowEditor(=kTRUE)");
               } else {
                  browser->StartEmbedding(TRootBrowser::kLeft);
                  if (!fEditor)
                     fEditor = TVirtualPadEditor::GetPadEditor(kTRUE);
                  else {
                     ((TGedEditor *)fEditor)->ReparentWindow(fClient->GetRoot());
                     ((TGedEditor *)fEditor)->MapWindow();
                  }
                  browser->StopEmbedding("Pad Editor");
                  if (fEditor) {
                     fEditor->SetGlobal(kFALSE);
                     gROOT->GetListOfCleanups()->Remove((TGedEditor *)fEditor);
                     ((TGedEditor *)fEditor)->SetCanvas(fCanvas);
                     ((TGedEditor *)fEditor)->SetModel(fCanvas, fCanvas, kButton1Down);
                  }
               }
            }
            else
               fEditor = TVirtualPadEditor::GetPadEditor(kFALSE);
         }
         if (show) browser->GetTabLeft()->SetTab("Pad Editor");
      }
   }
   else {
      if (show) {
         if (!fEditor) CreateEditor();
         TVirtualPadEditor* gged = TVirtualPadEditor::GetPadEditor(kFALSE);
         if(gged && gged->GetCanvas() == fCanvas){
            gged->Hide();
         }
         if (!fViewMenu->IsEntryChecked(kViewToolbar) || fToolDock->IsUndocked()) {
            ShowFrame(fHorizontal1);
            h = h + s;
         }
         fMainFrame->ShowFrame(fEditorFrame);
         fEditor->Show();
         fViewMenu->CheckEntry(kViewEditor);
         w = w + e;
      } else {
         if (!fViewMenu->IsEntryChecked(kViewToolbar) || fToolDock->IsUndocked()) {
            HideFrame(fHorizontal1);
            h = h - s;
         }
         if (fEditor) fEditor->Hide();
         fMainFrame->HideFrame(fEditorFrame);
         fViewMenu->UnCheckEntry(kViewEditor);
         w = w - e;
      }
      Resize(w, h);
   }
   if (savedPad) gPad = savedPad;
}

////////////////////////////////////////////////////////////////////////////////
/// Create embedded editor.

void TRootCanvas::CreateEditor()
{
   fEditorFrame->SetEditDisabled(kEditEnable);
   fEditorFrame->SetEditable();
   gPad = Canvas();
   // next two lines are related to the old editor
   Int_t show = gEnv->GetValue("Canvas.ShowEditor", 0);
   gEnv->SetValue("Canvas.ShowEditor","true");
   fEditor = TVirtualPadEditor::LoadEditor();
   if (fEditor) fEditor->SetGlobal(kFALSE);
   fEditorFrame->SetEditable(kEditDisable);
   fEditorFrame->SetEditable(kFALSE);

   // next line is related to the old editor
   if (show == 0) gEnv->SetValue("Canvas.ShowEditor","false");
}

////////////////////////////////////////////////////////////////////////////////
/// Show or hide toolbar.

void TRootCanvas::ShowToolBar(Bool_t show)
{
   if (show && !fToolBar) {

      fToolBar = new TGToolBar(fToolDock, 60, 20, kHorizontalFrame);
      fToolDock->AddFrame(fToolBar, fHorizontal1Layout);

      Int_t spacing = 6, i;
      for (i = 0; gToolBarData[i].fPixmap; i++) {
         if (strlen(gToolBarData[i].fPixmap) == 0) {
            spacing = 6;
            continue;
         }
         fToolBar->AddButton(this, &gToolBarData[i], spacing);
         spacing = 0;
      }
      fVertical1 = new TGVertical3DLine(fToolBar);
      fVertical2 = new TGVertical3DLine(fToolBar);
      fVertical1Layout = new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 4,2,0,0);
      fVertical2Layout = new TGLayoutHints(kLHintsLeft | kLHintsExpandY);
      fToolBar->AddFrame(fVertical1, fVertical1Layout);
      fToolBar->AddFrame(fVertical2, fVertical2Layout);

      spacing = 6;
      for (i = 0; gToolBarData1[i].fPixmap; i++) {
         if (strlen(gToolBarData1[i].fPixmap) == 0) {
            spacing = 6;
            continue;
         }
         fToolBar->AddButton(this, &gToolBarData1[i], spacing);
         spacing = 0;
      }
      fToolDock->MapSubwindows();
      fToolDock->Layout();
      fToolDock->SetWindowName(Form("ToolBar: %s", GetWindowName()));
      fToolDock->Connect("Docked()", "TRootCanvas", this, "AdjustSize()");
      fToolDock->Connect("Undocked()", "TRootCanvas", this, "AdjustSize()");
   }

   if (!fToolBar) return;

   UInt_t h = GetHeight();
   UInt_t sh = fToolBarSep->GetHeight();
   UInt_t dh = fToolBar->GetHeight();

   if (show) {
      ShowFrame(fToolDock);
      if (!fViewMenu->IsEntryChecked(kViewEditor)) {
         ShowFrame(fHorizontal1);
         h = h + sh;
      }
      ShowFrame(fToolBarSep);
      fViewMenu->CheckEntry(kViewToolbar);
      h = h + dh + sh;
   } else {
      if (fToolDock->IsUndocked()) {
         fToolDock->DockContainer();
         h = h + 2*sh;
      } else h = h - dh;

      HideFrame(fToolDock);
      if (!fViewMenu->IsEntryChecked(kViewEditor)) {
         HideFrame(fHorizontal1);
         h = h - sh;
      }
      HideFrame(fToolBarSep);
      h = h - sh;
      fViewMenu->UnCheckEntry(kViewToolbar);
   }
   Resize(GetWidth(), h);
}

////////////////////////////////////////////////////////////////////////////////
/// Enable or disable tooltip info.

void TRootCanvas::ShowToolTips(Bool_t show)
{
   if (show)
      fViewMenu->CheckEntry(kViewToolTips);
   else
      fViewMenu->UnCheckEntry(kViewToolTips);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if the editor is shown.

Bool_t TRootCanvas::HasEditor() const
{
   return (fEditor) && fViewMenu->IsEntryChecked(kViewEditor);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if the menu bar is shown.

Bool_t TRootCanvas::HasMenuBar() const
{
   return (fMenuBar) && fMenuBar->IsMapped();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if the status bar is shown.

Bool_t TRootCanvas::HasStatusBar() const
{
   return (fStatusBar) && fStatusBar->IsMapped();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if the tool bar is shown.

Bool_t TRootCanvas::HasToolBar() const
{
   return (fToolBar) && fToolBar->IsMapped();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if the tooltips are enabled.

Bool_t TRootCanvas::HasToolTips() const
{
   return (fCanvas) && fCanvas->GetShowToolTips();
}

////////////////////////////////////////////////////////////////////////////////
/// Keep the same canvas size while docking/undocking toolbar.

void TRootCanvas::AdjustSize()
{
   UInt_t h = GetHeight();
   UInt_t dh = fToolBar->GetHeight();
   UInt_t sh = fHorizontal1->GetHeight();

   if (fToolDock->IsUndocked()) {
      if (!fViewMenu->IsEntryChecked(kViewEditor)) {
         HideFrame(fHorizontal1);
         h = h - sh;
      }
      HideFrame(fToolBarSep);
      h = h - dh - sh;
   } else {
      if (!fViewMenu->IsEntryChecked(kViewEditor)) {
         ShowFrame(fHorizontal1);
         h = h + sh;
      }
      ShowFrame(fToolBarSep);
      h = h + dh + sh;
   }
   Resize(GetWidth(), h);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button events in the canvas container.

Bool_t TRootCanvas::HandleContainerButton(Event_t *event)
{
   Int_t button = event->fCode;
   Int_t x = event->fX;
   Int_t y = event->fY;

   if (event->fType == kButtonPress) {
      if (fToolTip && fCanvas->GetShowToolTips()) {
         fToolTip->Hide();
         gVirtualX->UpdateWindow(0);
         gSystem->ProcessEvents();
      }
      fButton = button;
      if (button == kButton1) {
         if (event->fState & kKeyShiftMask)
            fCanvas->HandleInput(kButton1Shift, x, y);
         else
            fCanvas->HandleInput(kButton1Down, x, y);
      }
      if (button == kButton2)
         fCanvas->HandleInput(kButton2Down, x, y);
      if (button == kButton3) {
         fCanvas->HandleInput(kButton3Down, x, y);
         fButton = 0;  // button up is consumed by TContextMenu
      }

   } else if (event->fType == kButtonRelease) {
      if (button == kButton4)
         fCanvas->HandleInput(kWheelUp, x, y);
      if (button == kButton5)
         fCanvas->HandleInput(kWheelDown, x, y);
      if (button == kButton1)
         fCanvas->HandleInput(kButton1Up, x, y);
      if (button == kButton2)
         fCanvas->HandleInput(kButton2Up, x, y);
      if (button == kButton3)
         fCanvas->HandleInput(kButton3Up, x, y);

      fButton = 0;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button double click events in the canvas container.

Bool_t TRootCanvas::HandleContainerDoubleClick(Event_t *event)
{
   Int_t button = event->fCode;
   Int_t x = event->fX;
   Int_t y = event->fY;

   if (button == kButton1)
      fCanvas->HandleInput(kButton1Double, x, y);
   if (button == kButton2)
      fCanvas->HandleInput(kButton2Double, x, y);
   if (button == kButton3)
      fCanvas->HandleInput(kButton3Double, x, y);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle configure (i.e. resize) event.

Bool_t TRootCanvas::HandleContainerConfigure(Event_t *)
{
   if (fAutoFit) {
      fCanvas->Resize();
      fCanvas->Update();
   }

   if (fCanvas->HasFixedAspectRatio()) {
      // get menu height
      static Int_t dh = 0;
      if (!dh)
         dh = GetHeight() - fCanvasContainer->GetHeight();
      UInt_t h = TMath::Nint(fCanvasContainer->GetWidth()/
                             fCanvas->GetAspectRatio()) + dh;
      SetWindowSize(GetWidth(), h);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle keyboard events in the canvas container.

Bool_t TRootCanvas::HandleContainerKey(Event_t *event)
{
   static EGEventType previous_event = kOtherEvent;
   static UInt_t previous_keysym = 0;

   if (event->fType == kGKeyPress) {
      fButton = event->fCode;
      UInt_t keysym;
      char str[2];
      gVirtualX->LookupString(event, str, sizeof(str), keysym);

      if (str[0] == kESC){   // ESC sets the escape flag
         gROOT->SetEscape();
         fCanvas->HandleInput(kButton1Up, 0, 0);
         fCanvas->HandleInput(kMouseMotion, 0, 0);
         gPad->Modified();
         return kTRUE;
      }
      if (str[0] == 3)   // ctrl-c sets the interrupt flag
         gROOT->SetInterrupt();

      // handle arrow keys
      if (keysym > 0x1011 && keysym < 0x1016) {
         Window_t dum1, dum2, wid;
         UInt_t mask = 0;
         Int_t mx, my, tx, ty;
         wid = gVirtualX->GetDefaultRootWindow();
         gVirtualX->QueryPointer(wid, dum1, dum2, mx, my, mx, my, mask);
         gVirtualX->TranslateCoordinates(gClient->GetDefaultRoot()->GetId(),
                                         fCanvasContainer->GetId(),
                                         mx, my, tx, ty, dum1);
         fCanvas->HandleInput(kArrowKeyPress, tx, ty);
         // handle case where we got consecutive same keypressed events coming
         // from auto-repeat on Windows (as it fires only successive keydown events)
         if ((previous_keysym == keysym) && (previous_event == kGKeyPress)) {
            switch (keysym) {
               case 0x1012: // left
                  gVirtualX->Warp(--mx, my, wid); --tx;
                  break;
               case 0x1013: // up
                  gVirtualX->Warp(mx, --my, wid); --ty;
                  break;
               case 0x1014: // right
                  gVirtualX->Warp(++mx, my, wid); ++tx;
                  break;
               case 0x1015: // down
                  gVirtualX->Warp(mx, ++my, wid); ++ty;
                  break;
               default:
                  break;
            }
            fCanvas->HandleInput(kArrowKeyRelease, tx, ty);
         }
         previous_keysym = keysym;
      }
      else {
         fCanvas->HandleInput(kKeyPress, str[0], keysym);
      }
   } else if (event->fType == kKeyRelease) {
      UInt_t keysym;
      char str[2];
      gVirtualX->LookupString(event, str, sizeof(str), keysym);

      if (keysym > 0x1011 && keysym < 0x1016) {
         Window_t dum1, dum2, wid;
         UInt_t mask = 0;
         Int_t mx, my, tx, ty;
         wid = gVirtualX->GetDefaultRootWindow();
         gVirtualX->QueryPointer(wid, dum1, dum2, mx, my, mx, my, mask);
         switch (keysym) {
            case 0x1012: // left
               gVirtualX->Warp(--mx, my, wid);
               break;
            case 0x1013: // up
               gVirtualX->Warp(mx, --my, wid);
               break;
            case 0x1014: // right
               gVirtualX->Warp(++mx, my, wid);
               break;
            case 0x1015: // down
               gVirtualX->Warp(mx, ++my, wid);
               break;
            default:
               break;
         }
         gVirtualX->TranslateCoordinates(gClient->GetDefaultRoot()->GetId(),
                                         fCanvasContainer->GetId(),
                                         mx, my, tx, ty, dum1);
         fCanvas->HandleInput(kArrowKeyRelease, tx, ty);
         previous_keysym = keysym;
      }
      fButton = 0;
   }
   previous_event = event->fType;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion event in the canvas container.

Bool_t TRootCanvas::HandleContainerMotion(Event_t *event)
{
   Int_t x = event->fX;
   Int_t y = event->fY;

   if (fButton == 0)
      fCanvas->HandleInput(kMouseMotion, x, y);
   if (fButton == kButton1) {
      if (event->fState & kKeyShiftMask)
         fCanvas->HandleInput(EEventType(8), x, y);
      else
         fCanvas->HandleInput(kButton1Motion, x, y);
   }
   if (fButton == kButton2)
      fCanvas->HandleInput(kButton2Motion, x, y);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle expose events.

Bool_t TRootCanvas::HandleContainerExpose(Event_t *event)
{
   if (event->fCount == 0) {
      fCanvas->Flush();
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle enter/leave events. Only leave is activated at the moment.

Bool_t TRootCanvas::HandleContainerCrossing(Event_t *event)
{
   Int_t x = event->fX;
   Int_t y = event->fY;

   // pointer grabs create also an enter and leave event but with fCode
   // either kNotifyGrab or kNotifyUngrab, don't propagate these events
   if (event->fType == kLeaveNotify && event->fCode == kNotifyNormal)
      fCanvas->HandleInput(kMouseLeave, x, y);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle drop events.

Bool_t TRootCanvas::HandleDNDDrop(TDNDData *data)
{
   static Atom_t rootObj  = gVirtualX->InternAtom("application/root", kFALSE);
   static Atom_t uriObj  = gVirtualX->InternAtom("text/uri-list", kFALSE);

   if (data->fDataType == rootObj) {
      TBufferFile buf(TBuffer::kRead, data->fDataLength, (void *)data->fData);
      buf.SetReadMode();
      TObject *obj = (TObject *)buf.ReadObjectAny(TObject::Class());
      if (!obj) return kTRUE;
      gPad->Clear();
      if (obj->InheritsFrom("TKey")) {
         TObject *object = (TObject *)gROOT->ProcessLine(Form("((TKey *)0x%lx)->ReadObj();", (ULong_t)obj));
         if (!object) return kTRUE;
         if (object->InheritsFrom("TGraph"))
            object->Draw("ALP");
         else if (object->InheritsFrom("TImage"))
            object->Draw("x");
         else if (object->IsA()->GetMethodAllAny("Draw"))
            object->Draw();
      }
      else if (obj->InheritsFrom("TGraph"))
         obj->Draw("ALP");
      else if (obj->IsA()->GetMethodAllAny("Draw"))
         obj->Draw();
      gPad->Modified();
      gPad->Update();
      return kTRUE;
   }
   else if (data->fDataType == uriObj) {
      TString sfname((char *)data->fData);
      if (sfname.Length() > 7) {
         sfname.ReplaceAll("\r\n", "");
         TUrl uri(sfname.Data());
         if (sfname.EndsWith(".bmp") ||
            sfname.EndsWith(".gif") ||
            sfname.EndsWith(".jpg") ||
            sfname.EndsWith(".png") ||
            sfname.EndsWith(".ps")  ||
            sfname.EndsWith(".eps") ||
            sfname.EndsWith(".pdf") ||
            sfname.EndsWith(".tiff") ||
            sfname.EndsWith(".xpm")) {
            TImage *img = TImage::Open(uri.GetFile());
            if (img) {
               img->Draw("x");
               img->SetEditable(kTRUE);
            }
         }
         gPad->Modified();
         gPad->Update();
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle dragging position events.

Atom_t TRootCanvas::HandleDNDPosition(Int_t x, Int_t y, Atom_t action,
                                      Int_t /*xroot*/, Int_t /*yroot*/)
{
   TPad *pad = fCanvas->Pick(x, y, 0);
   if (pad) {
      pad->cd();
      gROOT->SetSelectedPad(pad);
      // make sure the pad is highlighted (on Windows)
      pad->Update();
   }
   return action;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle drag enter events.

Atom_t TRootCanvas::HandleDNDEnter(Atom_t *typelist)
{
   static Atom_t rootObj  = gVirtualX->InternAtom("application/root", kFALSE);
   static Atom_t uriObj  = gVirtualX->InternAtom("text/uri-list", kFALSE);
   Atom_t ret = kNone;
   for (int i = 0; typelist[i] != kNone; ++i) {
      if (typelist[i] == rootObj)
         ret = rootObj;
      if (typelist[i] == uriObj)
         ret = uriObj;
   }
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle drag leave events.

Bool_t TRootCanvas::HandleDNDLeave()
{
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Slot handling tab switching in the browser, to properly set the canvas
/// and the model to the editor.

void TRootCanvas::Activated(Int_t id)
{
   if (fEmbedded) {
      TGTab *sender = (TGTab *)gTQSender;
      if (sender) {
         TGCompositeFrame *cont = sender->GetTabContainer(id);
         if (cont == fParent) {
            if (!fEditor)
               fEditor = TVirtualPadEditor::GetPadEditor(kFALSE);
            if (fEditor && ((TGedEditor *)fEditor)->IsMapped()) {
               ((TGedEditor *)fEditor)->SetCanvas(fCanvas);
               ((TGedEditor *)fEditor)->SetModel(fCanvas, fCanvas, kButton1Down);
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save a canvas container as a C++ statement(s) on output stream out.

void TRootContainer::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   out << std::endl << "   // canvas container" << std::endl;
   out << "   Int_t canvasID = gVirtualX->InitWindow((ULong_t)"
       << GetParent()->GetParent()->GetName() << "->GetId());" << std::endl;
   out << "   Window_t winC = gVirtualX->GetWindowID(canvasID);" << std::endl;
   out << "   TGCompositeFrame *";
   out << GetName() << " = new TGCompositeFrame(gClient,winC"
       << "," << GetParent()->GetName() << ");" << std::endl;
}
