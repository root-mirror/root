// @(#)root/gl:$Id$
// Author:  Timur Pocheptsov / Richard Maunder

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <memory>

#include "TRootHelpDialog.h"
#include "TPluginManager.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TGCanvas.h"
#include "HelpText.h"
#include "GuiTypes.h"
#include "TG3DLine.h"
#include "TSystem.h"
#include "TGFrame.h"
#include "TGMenu.h"
#include "TGSplitter.h"
#include "TGButton.h"
#include "snprintf.h"

#include "TVirtualPad.h"
#include "TGedEditor.h"
#include "TRootEmbeddedCanvas.h"
#include "TString.h"
#include "TGFileDialog.h"
#include "TVirtualX.h"

#include "TGLOutput.h"
#include "TGLFormat.h"

#include "TGLLogicalShape.h"
#include "TGLPhysicalShape.h"
#include "TGLPShapeObj.h"
#include "TGLClip.h"
#include "TROOT.h"

#ifdef WIN32
#include "TWin32SplashThread.h"
#endif

#include "TGLWidget.h"
#include "TGLSAViewer.h"
#include "TGLSAFrame.h"
#include "TGLEventHandler.h"


const char * TGLSAViewer::fgHelpText1 = "\
DIRECT SCENE INTERACTIONS\n\n\
   Press:\n\
   \tw          --- wireframe mode\n\
   \te          --- switch between dark / light color-set\n\
   \tr          --- filled polygons mode\n\
   \tt          --- outline mode\n\
   \tj          --- ZOOM in\n\
   \tk          --- ZOOM out\n\
   \ta          --- switch on/off arc-ball camera rotation control\n\
   \tArrow Keys --- PAN (TRUCK) across scene\n\
   \tHome       --- reset current camera\n\
   \tCtrl-Home  --- switch external/automatic camera center\n\
\n\
   LEFT mouse button -- ROTATE (ORBIT) the scene by holding the mouse button and moving\n\
   the mouse (perspective camera, needs to be enabled in menu for orthographic cameras).\n\
   By default, the scene will be rotated about its center. To select arbitrary center\n\
   bring up the viewer-editor (e.g., shift-click into empty background) and use\n\
   'Camera center' controls in the 'Guides' tab.\n\
\n\
   MIDDLE mouse button or arrow keys --  PAN (TRUCK) the camera.\n\
\n\
   RIGHT mouse button action depends on camera type:\n\
     orthographic -- zoom,\n\
     perspective  -- move camera forwards / backwards\n\
\n\
   By pressing Ctrl and Shift keys the mouse precision can be changed:\n\
     Shift      -- 10 times less precise\n\
     Ctrl       -- 10 times more precise\n\
     Ctrl Shift -- 100 times more precise\n\
\n\
   Mouse wheel action depends on camera type:\n\
     orthographic -- zoom,\n\
     perspective  -- change field-of-view (focal length)\n\
\n\
   To invert direction of mouse and key actions from scene-centric\n\
   to viewer-centric, set in your .rootrc file:\n\
      OpenGL.EventHandler.ViewerCentricControls: 1\n\
\n\
   Double click will show GUI editor of the viewer (if assigned).\n\
\n\
   RESET the camera via the button in viewer-editor or Home key.\n\
\n\
   SELECT a shape with Shift+Left mouse button click.\n\
\n\
   SELECT the viewer with Shift+Left mouse button click on a free space.\n\
\n\
   MOVE a selected shape using Shift+Mid mouse drag.\n\
\n\
   Invoke the CONTEXT menu with Shift+Right mouse click.\n\n"
   "Secondary selection and direct render object interaction is initiated\n\
   by Alt+Left mouse click (Mod1, actually). Only few classes support this option.\n\
   When 'Alt' is taken by window manager, try Alt-Ctrl-Left.\n\
\n\
CAMERA\n\
\n\
   The \"Camera\" menu is used to select the different projections from \n\
   the 3D world onto the 2D viewport. There are three perspective cameras:\n\
\n\
   \tPerspective (Floor XOZ)\n\
   \tPerspective (Floor YOZ)\n\
   \tPerspective (Floor XOY)\n\
\n\
   In each case the floor plane (defined by two axes) is kept level.\n\
\n\
   There are also four orthographic cameras:\n\
\n\
   \tOrthographic (XOY)\n\
   \tOrthographic (XOZ)\n\
   \tOrthographic (ZOY)\n\
   \tOrthographic (ZOX)\n\
\n\
   In each case the first axis is placed horizontal, the second vertical e.g.\n\
   XOY means X horizontal, Y vertical.\n\n";

const char * TGLSAViewer::fgHelpText2 = "\
SHAPES COLOR AND MATERIAL\n\
\n\
   The selected shape's color can be modified in the Shapes-Color tabs.\n\
   Shape's color is specified by the percentage of red, green, blue light\n\
   it reflects. A surface can reflect DIFFUSE, AMBIENT and SPECULAR light.\n\
   A surface can also emit light. The EMISSIVE parameter allows to define it.\n\
   The surface SHININESS can also be modified.\n\
\n\
SHAPES GEOMETRY\n\
\n\
   The selected shape's location and geometry can be modified in the Shapes-Geom\n\
   tabs by entering desired values in respective number entry controls.\n\
\n\
SCENE CLIPPING\n\
\n\
   In the Scene-Clipping tabs select a 'Clip Type': None, Plane, Box\n\
\n\
   For 'Plane' and 'Box' the lower pane shows the relevant parameters:\n\
\n\
\tPlane: Equation coefficients of form aX + bY + cZ + d = 0\n\
\tBox: Center X/Y/Z and Length X/Y/Z\n\n"
   "For Box checking the 'Show / Edit' checkbox shows the clip box (in light blue)\n\
   in viewer. It also attaches the current manipulator to the box - enabling\n\
   direct editing in viewer.\n\
\n\
MANIPULATORS\n\
\n\
   A widget attached to the selected object - allowing direct manipulation\n\
   of the object with respect to its local axes.\n\
\n\
   There are three modes, toggled with keys while manipulator is active, that is,\n\
   mouse pointer is above it (switches color to yellow):\n\
   \tMode\t\tWidget Component Style\t\tKey\n\
   \t----\t\t----------------------\t\t---\n\
   \tTranslation\tLocal axes with arrows\t\tv\n\
   \tScale\t\tLocal axes with boxes\t\tx\n\
   \tRotate\t\tLocal axes rings\t\tc\n\
\n\
   Each widget has three axis components - red (X), green (Y) and blue (Z).\n\
   The component turns yellow, indicating an active state, when the mouse is moved\n\
   over it. Left click and drag on the active component to adjust the objects\n\
   translation, scale or rotation.\n\
   Some objects do not support all manipulations (e.g. clipping planes cannot be \n\
   scaled). If a manipulation is not permitted the component it drawn in grey and \n\
   cannot be selected/dragged.\n";


/** \class TGLSAViewer
\ingroup opengl
The top level standalone GL-viewer - created via plugin manager.
*/

ClassImp(TGLSAViewer);

Long_t TGLSAViewer::fgMenuHidingTimeout = 400;

const Int_t TGLSAViewer::fgInitX = 0;
const Int_t TGLSAViewer::fgInitY = 0;
const Int_t TGLSAViewer::fgInitW = 780;
const Int_t TGLSAViewer::fgInitH = 670;

// A lot of raw pointers/naked new-expressions - good way to discredit C++ (or C++ programmer
// ROOT has system to cleanup - I'll try to use it

const char *gGLSaveAsTypes[] = {"Encapsulated PostScript", "*.eps",
                                "PDF",                     "*.pdf",
                                "GIF",                     "*.gif",
                                "Animated GIF",            "*.gif+",
                                "JPEG",                    "*.jpg",
                                "PNG",                     "*.png",
                                0, 0};

////////////////////////////////////////////////////////////////////////////////
/// Construct a standalone viewer, bound to supplied 'pad'.

TGLSAViewer::TGLSAViewer(TVirtualPad *pad, TGLFormat* format) :
   TGLViewer(pad, fgInitX, fgInitY, fgInitW, fgInitH),
   fFrame(0),
   fFormat(format),
   fFileMenu(0),
   fFileSaveMenu(0),
   fCameraMenu(0),
   fHelpMenu(0),
   fLeftVerticalFrame(0),
   fRightVerticalFrame(0),
   fDirName("."),
   fTypeIdx(0),
   fOverwrite(kFALSE),
   fMenuBar(0),
   fMenuBut(0),
   fHideMenuBar(kFALSE),
   fMenuHidingTimer(0),
   fMenuHidingShowMenu(kTRUE),
   fDeleteMenuBar(kFALSE)
{
   fFrame = new TGLSAFrame(*this);

   CreateMenus();
   CreateFrames();

   fFrame->SetWindowName("ROOT's GL viewer");
   fFrame->SetClassHints("GLViewer", "GLViewer");
   fFrame->SetMWMHints(kMWMDecorAll, kMWMFuncAll, kMWMInputModeless);
   fFrame->MapSubwindows();
   fFrame->HideFrame(fMenuBut);

   fFrame->Resize(fFrame->GetDefaultSize());
   fFrame->MoveResize(fgInitX, fgInitY, fgInitW, fgInitH);
   fFrame->SetWMPosition(fgInitX, fgInitY);

   // set recursive cleanup, but exclude fGedEditor
   // destructor of fGedEditor has own way of handling child nodes
   TObject* fe = fLeftVerticalFrame->GetList()->First();
   fLeftVerticalFrame->GetList()->Remove(fe);
   fFrame->SetCleanup(kDeepCleanup);
   fLeftVerticalFrame->GetList()->AddFirst(fe);

   Show();
}

////////////////////////////////////////////////////////////////////////////////
/// Construct an embedded standalone viewer, bound to supplied 'pad'.
/// If format is passed, it gets adopted by the viewer as it might
/// need to be reused several times when recreating the GL-widget.
///
/// Modified version of the previous constructor for embedding the
/// viewer into another frame (parent).

TGLSAViewer::TGLSAViewer(const TGWindow *parent, TVirtualPad *pad, TGedEditor *ged,
                         TGLFormat* format) :
   TGLViewer(pad, fgInitX, fgInitY, fgInitW, fgInitH),
   fFrame(0),
   fFormat(format),
   fFileMenu(0),
   fCameraMenu(0),
   fHelpMenu(0),
   fLeftVerticalFrame(0),
   fRightVerticalFrame(0),
   fTypeIdx(0),
   fMenuBar(0),
   fMenuBut(0),
   fHideMenuBar(kFALSE),
   fMenuHidingTimer(0),
   fMenuHidingShowMenu(kTRUE),
   fDeleteMenuBar(kFALSE)
{
   fGedEditor = ged;
   fFrame = new TGLSAFrame(parent, *this);

   CreateMenus();
   CreateFrames();

   fFrame->MapSubwindows();
   fFrame->HideFrame(fMenuBut);
   fFrame->Resize(fFrame->GetDefaultSize());
   fFrame->Resize(fgInitW, fgInitH);

   // set recursive cleanup, but exclude fGedEditor
   // destructor of fGedEditor has own way of handling child nodes
   if (fLeftVerticalFrame)
   {
      TObject* fe = fLeftVerticalFrame->GetList()->First();
      fLeftVerticalFrame->GetList()->Remove(fe);
      fFrame->SetCleanup(kDeepCleanup);
      fLeftVerticalFrame->GetList()->AddFirst(fe);
   }

   Show();
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy standalone viewer object.

TGLSAViewer::~TGLSAViewer()
{
   fGedEditor->DisconnectFromCanvas();

   DisableMenuBarHiding();

   delete fHelpMenu;
   delete fCameraMenu;
   delete fFileSaveMenu;
   delete fFileMenu;
   if(fDeleteMenuBar) {
      delete fMenuBar;
   }
   delete fFormat;
   delete fFrame;
   fGLWidget = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the main-frame.

TGCompositeFrame* TGLSAViewer::GetFrame() const
{
   return fFrame;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a GLwidget, it is an error if it is already created.
/// This is needed for frame-swapping on mac.

void TGLSAViewer::CreateGLWidget()
{
   if (fGLWidget) {
      Error("CreateGLWidget", "Widget already exists.");
      return;
   }

   if (fFormat == 0)
      fFormat = new TGLFormat;

   fGLWidget = TGLWidget::Create(*fFormat, fRightVerticalFrame, kTRUE, kTRUE, 0, 10, 10);
   fGLWidget->SetEventHandler(fEventHandler);

   fRightVerticalFrame->AddFrame(fGLWidget, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   fFrame->Layout();

   fGLWidget->MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy the GLwidget, it is an error if it does not exist.
/// This is needed for frame-swapping on mac.

void TGLSAViewer::DestroyGLWidget()
{
   if (fGLWidget == 0) {
      Error("DestroyGLWidget", "Widget does not exist.");
      return;
   }

   fGLWidget->UnmapWindow();
   fGLWidget->SetEventHandler(0);

   fRightVerticalFrame->RemoveFrame(fGLWidget);
   fGLWidget->DeleteWindow();
   fGLWidget = 0;
}

////////////////////////////////////////////////////////////////////////////////
///File/Camera/Help menus.

void TGLSAViewer::CreateMenus()
{
   fFileMenu = new TGPopupMenu(fFrame->GetClient()->GetDefaultRoot());
   fFileMenu->AddEntry("&Hide Menus", kGLHideMenus);
   fFileMenu->AddEntry("&Edit Object", kGLEditObject);
   fFileMenu->AddSeparator();
   fFileMenu->AddEntry("&Close Viewer", kGLCloseViewer);
   fFileMenu->AddSeparator();
   fFileSaveMenu = new TGPopupMenu(fFrame->GetClient()->GetDefaultRoot());
   fFileSaveMenu->AddEntry("viewer.&eps", kGLSaveEPS);
   fFileSaveMenu->AddEntry("viewer.&pdf", kGLSavePDF);
   fFileSaveMenu->AddEntry("viewer.&gif", kGLSaveGIF);
   fFileSaveMenu->AddEntry("viewer.g&if+", kGLSaveAnimGIF);
   fFileSaveMenu->AddEntry("viewer.&jpg", kGLSaveJPG);
   fFileSaveMenu->AddEntry("viewer.p&ng", kGLSavePNG);
   fFileMenu->AddPopup("&Save", fFileSaveMenu);
   fFileMenu->AddEntry("Save &As...", kGLSaveAS);
   fFileMenu->AddSeparator();
   fFileMenu->AddEntry("&Quit ROOT", kGLQuitROOT);
   fFileMenu->Associate(fFrame);

   fCameraMenu = new TGPopupMenu(fFrame->GetClient()->GetDefaultRoot());
   fCameraMenu->AddEntry("Perspective (Floor XOZ)", kGLPerspXOZ);
   fCameraMenu->AddEntry("Perspective (Floor YOZ)", kGLPerspYOZ);
   fCameraMenu->AddEntry("Perspective (Floor XOY)", kGLPerspXOY);
   fCameraMenu->AddEntry("Orthographic (XOY)", kGLXOY);
   fCameraMenu->AddEntry("Orthographic (XOZ)", kGLXOZ);
   fCameraMenu->AddEntry("Orthographic (ZOY)", kGLZOY);
   fCameraMenu->AddEntry("Orthographic (ZOX)", kGLZOX);
   fCameraMenu->AddEntry("Orthographic (XnOY)", kGLXnOY);
   fCameraMenu->AddEntry("Orthographic (XnOZ)", kGLXnOZ);
   fCameraMenu->AddEntry("Orthographic (ZnOY)", kGLZnOY);
   fCameraMenu->AddEntry("Orthographic (ZnOX)", kGLZnOX);
   fCameraMenu->AddSeparator();
   fCameraMenu->AddEntry("Ortho allow rotate", kGLOrthoRotate);
   fCameraMenu->AddEntry("Ortho allow dolly",  kGLOrthoDolly);
   fCameraMenu->Associate(fFrame);

   fHelpMenu = new TGPopupMenu(fFrame->GetClient()->GetDefaultRoot());
   fHelpMenu->AddEntry("Help on GL Viewer...", kGLHelpViewer);
   fHelpMenu->AddSeparator();
   fHelpMenu->AddEntry("&About ROOT...", kGLHelpAbout);
   fHelpMenu->Associate(fFrame);

   // Create menubar
   fMenuBar = new TGMenuBar(fFrame);
   fMenuBar->AddPopup("&File", fFileMenu, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Camera", fCameraMenu, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Help",   fHelpMenu,   new TGLayoutHints(kLHintsTop | kLHintsRight));
   fFrame->AddFrame(fMenuBar, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 1, 1));
   gVirtualX->SelectInput(fMenuBar->GetId(),
                          kKeyPressMask | kExposureMask | kPointerMotionMask
                          | kStructureNotifyMask | kFocusChangeMask
                          | kEnterWindowMask | kLeaveWindowMask);

   fMenuBut = new TGButton(fFrame);
   fMenuBut->ChangeOptions(kRaisedFrame | kFixedHeight);
   fMenuBut->Resize(20, 4);
   fMenuBut->SetBackgroundColor(0x80A0C0);
   fFrame->AddFrame(fMenuBut, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0, 0, 1, 1));
}

////////////////////////////////////////////////////////////////////////////////
/// Internal frames creation.

void TGLSAViewer::CreateFrames()
{
   TGCompositeFrame* compositeFrame = fFrame;
   if (fGedEditor == 0)
   {
      compositeFrame = new TGCompositeFrame(fFrame, 100, 100, kHorizontalFrame | kRaisedFrame);
      fFrame->AddFrame(compositeFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

      fLeftVerticalFrame = new TGVerticalFrame(compositeFrame, 195, 10, kFixedWidth);
      compositeFrame->AddFrame(fLeftVerticalFrame, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 2, 2, 2));

      const TGWindow* cw =  fFrame->GetClient()->GetRoot();
      fFrame->GetClient()->SetRoot(fLeftVerticalFrame);

      fGedEditor = new TGedEditor();
      fGedEditor->GetTGCanvas()->ChangeOptions(0);
      fLeftVerticalFrame->RemoveFrame(fGedEditor);
      fLeftVerticalFrame->AddFrame(fGedEditor, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 0, 0, 2, 2));
      fLeftVerticalFrame->GetClient()->SetRoot((TGWindow*)cw);
      fLeftVerticalFrame->MapSubwindows();

      TGVSplitter *splitter = new TGVSplitter(compositeFrame);
      splitter->SetFrame(fLeftVerticalFrame, kTRUE);
      compositeFrame->AddFrame(splitter, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 0,1,2,2) );
   }

   // SunkenFrame introduces 1-pixel offset - in TGFrame.cxx:163
   //
   // TGVerticalFrame *rightVerticalFrame = new TGVerticalFrame(compositeFrame, 10, 10, kSunkenFrame);
   // compositeFrame->AddFrame(rightVerticalFrame, new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY,0,2,2,2));
   fRightVerticalFrame = new TGVerticalFrame(compositeFrame, 10, 10);
   compositeFrame->AddFrame(fRightVerticalFrame, new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY));

   fEventHandler = new TGLEventHandler(0, this);
   CreateGLWidget();
}

////////////////////////////////////////////////////////////////////////////////
/// Update GUI components for embedded viewer selection change.
/// Override from TGLViewer.

void TGLSAViewer::SelectionChanged()
{
   TGLPhysicalShape *selected = const_cast<TGLPhysicalShape*>(GetSelected());

   if (selected) {
      fPShapeWrap->fPShape = selected;
      if (fFileMenu->IsEntryChecked(kGLEditObject))
         fGedEditor->SetModel(fPad, selected->GetLogical()->GetExternal(), kButton1Down);
      else
         fGedEditor->SetModel(fPad, fPShapeWrap, kButton1Down);
   } else {
      fPShapeWrap->fPShape = 0;
      fGedEditor->SetModel(fPad, this, kButton1Down);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Show the viewer

void TGLSAViewer::Show()
{
   fFrame->MapRaised();
   fGedEditor->SetModel(fPad, this, kButton1Down);
   RequestDraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Close the viewer - destructed.

void TGLSAViewer::Close()
{
   // Commit suicide when contained GUI is closed.
   delete this;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete the menu bar.

void TGLSAViewer::DeleteMenuBar()
{
   fDeleteMenuBar=kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Deactivate menu entries for closing the GL window and exiting ROOT.

void TGLSAViewer::DisableCloseMenuEntries()
{
   fFileMenu->DeleteEntry(kGLCloseViewer);
   fFileMenu->DeleteEntry(kGLQuitROOT);
}

////////////////////////////////////////////////////////////////////////////////
/// Enable hiding of menu bar.

void TGLSAViewer::EnableMenuBarHiding()
{
   if (fHideMenuBar)
      return;

   fHideMenuBar = kTRUE;

   fMenuBar->Connect("ProcessedEvent(Event_t*)", "TGLSAViewer", this, "HandleMenuBarHiding(Event_t*)");
   fMenuBut->Connect("ProcessedEvent(Event_t*)", "TGLSAViewer", this, "HandleMenuBarHiding(Event_t*)");

   fFrame->HideFrame(fMenuBar);
   fFrame->ShowFrame(fMenuBut);
   fFrame->Layout();

   fMenuHidingTimer = new TTimer;
   fMenuHidingTimer->Connect("Timeout()", "TGLSAViewer", this, "MenuHidingTimeout()");

   fFileMenu->CheckEntry(kGLHideMenus);
}

////////////////////////////////////////////////////////////////////////////////
/// Disable hiding of menu bar.

void TGLSAViewer::DisableMenuBarHiding()
{
   if (!fHideMenuBar)
      return;

   fHideMenuBar = kFALSE;

   fMenuBar->Disconnect("ProcessedEvent(Event_t*)", this, "HandleMenuBarHiding(Event_t*)");
   fMenuBut->Disconnect("ProcessedEvent(Event_t*)", this, "HandleMenuBarHiding(Event_t*)");

   fFrame->ShowFrame(fMenuBar);
   fFrame->HideFrame(fMenuBut);
   fFrame->Layout();

   fMenuHidingTimer->TurnOff();
   delete fMenuHidingTimer;
   fMenuHidingTimer = 0;

   fFileMenu->UnCheckEntry(kGLHideMenus);
}

////////////////////////////////////////////////////////////////////////////////
/// Maybe switch menu-bar / menu-button.

void TGLSAViewer::HandleMenuBarHiding(Event_t* ev)
{
   TGFrame *f = (TGFrame*) gTQSender;

   if (f == fMenuBut)
   {
      if (ev->fType == kEnterNotify)
         ResetMenuHidingTimer(kTRUE);
      else
         fMenuHidingTimer->TurnOff();
   }
   else if (f == fMenuBar)
   {
      if (ev->fType == kLeaveNotify &&
          (ev->fX < 0 || ev->fX >= (Int_t) f->GetWidth() ||
           ev->fY < 0 || ev->fY >= (Int_t) f->GetHeight()))
      {
         if (fMenuBar->GetCurrent() == 0)
            ResetMenuHidingTimer(kFALSE);
         else
            fMenuBar->GetCurrent()->Connect("ProcessedEvent(Event_t*)", "TGLSAViewer", this, "HandleMenuBarHiding(Event_t*)");
      }
      else
      {
         fMenuHidingTimer->TurnOff();
      }
   }
   else
   {
      f->Disconnect("ProcessedEvent(Event_t*)", this);
      ResetMenuHidingTimer(kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the timer for menu-bar hiding.

void TGLSAViewer::ResetMenuHidingTimer(Bool_t show_menu)
{
   // This happens, mysteriously.
   if (fMenuHidingTimer == 0)
      return;

   fMenuHidingTimer->TurnOff();

   fMenuHidingShowMenu = show_menu;

   fMenuHidingTimer->SetTime(fgMenuHidingTimeout);
   fMenuHidingTimer->Reset();
   fMenuHidingTimer->TurnOn();
}

////////////////////////////////////////////////////////////////////////////////
/// Action for menu-hiding timeout.

void TGLSAViewer::MenuHidingTimeout()
{
   fMenuHidingTimer->TurnOff();
   if (fMenuHidingShowMenu) {
      fFrame->HideFrame(fMenuBut);
      fFrame->ShowFrame(fMenuBar);
   } else {
      fFrame->HideFrame(fMenuBar);
      fFrame->ShowFrame(fMenuBut);
   }
   fFrame->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Set global timeout for menu-hiding in mili-seconds.
/// Static function.

void TGLSAViewer::SetMenuHidingTimeout(Long_t timeout)
{
   fgMenuHidingTimeout = timeout;
}

////////////////////////////////////////////////////////////////////////////////
/// Process GUI message capture by the main GUI frame (TGLSAFrame).

Bool_t TGLSAViewer::ProcessFrameMessage(Long_t msg, Long_t parm1, Long_t)
{
   switch (GET_MSG(msg)) {
   case kC_COMMAND:
      switch (GET_SUBMSG(msg)) {
      case kCM_BUTTON:
      case kCM_MENU:
         switch (parm1) {
         case kGLHelpAbout: {
#ifdef R__UNIX
            TString rootx = TROOT::GetBinDir() + "/root -a &";
            gSystem->Exec(rootx);
#else
#ifdef WIN32
            new TWin32SplashThread(kTRUE);
#else
            char str[32];
            snprintf(str,32, "About ROOT %s...", gROOT->GetVersion());
            hd = new TRootHelpDialog(this, str, 600, 400);
            hd->SetText(gHelpAbout);
            hd->Popup();
#endif
#endif
            break;
         }
         case kGLHelpViewer: {
            TRootHelpDialog * hd = new TRootHelpDialog(fFrame, "Help on GL Viewer...", 660, 400);
            hd->AddText(fgHelpText1);
            hd->AddText(fgHelpText2);
            hd->Popup();
            break;
         }
         case kGLPerspYOZ:
            SetCurrentCamera(TGLViewer::kCameraPerspYOZ);
            break;
         case kGLPerspXOZ:
            SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
            break;
         case kGLPerspXOY:
            SetCurrentCamera(TGLViewer::kCameraPerspXOY);
            break;
         case kGLXOY:
            SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
            break;
         case kGLXOZ:
            SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
            break;
         case kGLZOY:
            SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
            break;
         case kGLZOX:
            SetCurrentCamera(TGLViewer::kCameraOrthoZOX);
            break;
         case kGLXnOY:
            SetCurrentCamera(TGLViewer::kCameraOrthoXnOY);
            break;
         case kGLXnOZ:
            SetCurrentCamera(TGLViewer::kCameraOrthoXnOZ);
            break;
         case kGLZnOY:
            SetCurrentCamera(TGLViewer::kCameraOrthoZnOY);
            break;
         case kGLZnOX:
            SetCurrentCamera(TGLViewer::kCameraOrthoZnOX);
            break;
         case kGLOrthoRotate:
            ToggleOrthoRotate();
            break;
         case kGLOrthoDolly:
            ToggleOrthoDolly();
            break;
         case kGLSaveEPS:
            SavePicture("viewer.eps");
            break;
         case kGLSavePDF:
            SavePicture("viewer.pdf");
            break;
         case kGLSaveGIF:
            SavePicture("viewer.gif");
            break;
         case kGLSaveAnimGIF:
            SavePicture("viewer.gif+");
            break;
         case kGLSaveJPG:
            SavePicture("viewer.jpg");
            break;
         case kGLSavePNG:
            SavePicture("viewer.png");
            break;
         case kGLSaveAS:
            {
               TGFileInfo fi;
               fi.fFileTypes   = gGLSaveAsTypes;
               fi.SetIniDir(fDirName);
               fi.fFileTypeIdx = fTypeIdx;
               fi.fOverwrite   = fOverwrite;
               new TGFileDialog(gClient->GetDefaultRoot(), fFrame, kFDSave, &fi);
               if (!fi.fFilename) return kTRUE;
               TString ft(fi.fFileTypes[fi.fFileTypeIdx+1]);
               fDirName   = fi.fIniDir;
               fTypeIdx   = fi.fFileTypeIdx;
               fOverwrite = fi.fOverwrite;

               TString file = fi.fFilename;
               Bool_t  match = kFALSE;
               const char** fin = gGLSaveAsTypes; ++fin;
               while (*fin != 0)
               {
                  if (file.EndsWith(*fin + 1))
                  {
                     match = kTRUE;
                     break;
                  }
                  fin += 2;
               }
               if ( ! match)
               {
                  file += ft(ft.Index("."), ft.Length());
               }
               SavePicture(file);
            }
            break;
         case kGLHideMenus:
            if (fHideMenuBar)
               DisableMenuBarHiding();
            else
               EnableMenuBarHiding();
            break;
         case kGLEditObject:
            ToggleEditObject();
            break;
         case kGLCloseViewer:
            // Exit needs to be delayed to avoid bad drawable X ids - GUI
            // will all be changed in future anyway
            TTimer::SingleShot(50, "TGLSAFrame", fFrame, "SendCloseMessage()");
            break;
         case kGLQuitROOT:
            if (!gApplication->ReturnFromRun())
               delete this;
            gApplication->Terminate(0);
            break;
         default:
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
/// Toggle state of the 'Edit Object' menu entry.

void TGLSAViewer::ToggleEditObject()
{
   if (fFileMenu->IsEntryChecked(kGLEditObject))
      fFileMenu->UnCheckEntry(kGLEditObject);
   else
      fFileMenu->CheckEntry(kGLEditObject);
   SelectionChanged();
}

////////////////////////////////////////////////////////////////////////////////
/// Toggle state of the 'Ortho allow rotate' menu entry.

void TGLSAViewer::ToggleOrthoRotate()
{
   if (fCameraMenu->IsEntryChecked(kGLOrthoRotate))
      fCameraMenu->UnCheckEntry(kGLOrthoRotate);
   else
      fCameraMenu->CheckEntry(kGLOrthoRotate);
   Bool_t state = fCameraMenu->IsEntryChecked(kGLOrthoRotate);
   fOrthoXOYCamera.SetEnableRotate(state);
   fOrthoXOZCamera.SetEnableRotate(state);
   fOrthoZOYCamera.SetEnableRotate(state);
   fOrthoXnOYCamera.SetEnableRotate(state);
   fOrthoXnOZCamera.SetEnableRotate(state);
   fOrthoZnOYCamera.SetEnableRotate(state);
}

////////////////////////////////////////////////////////////////////////////////
/// Toggle state of the 'Ortho allow dolly' menu entry.

void TGLSAViewer::ToggleOrthoDolly()
{
   if (fCameraMenu->IsEntryChecked(kGLOrthoDolly))
      fCameraMenu->UnCheckEntry(kGLOrthoDolly);
   else
      fCameraMenu->CheckEntry(kGLOrthoDolly);
   Bool_t state = ! fCameraMenu->IsEntryChecked(kGLOrthoDolly);
   fOrthoXOYCamera.SetDollyToZoom(state);
   fOrthoXOZCamera.SetDollyToZoom(state);
   fOrthoZOYCamera.SetDollyToZoom(state);
}
