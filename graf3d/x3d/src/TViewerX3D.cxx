// @(#)root/x3d:$Id$
// Author: Rene Brun   05/09/99
/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TViewerX3D                                                           //
//                                                                      //
// C++ interface to the X3D viewer                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TViewerX3D.h"
#include "X3DBuffer.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TMath.h"
#include "TROOT.h"
#include "snprintf.h"

#include "TRootHelpDialog.h"
#include "TGClient.h"
#include "TGCanvas.h"
#include "TGMenu.h"
#include "TGWidget.h"
#include "TGMsgBox.h"
#include "TVirtualX.h"

#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"

#include "HelpText.h"

#include <cassert>

const char gHelpX3DViewer[] = "\
     PRESS \n\
     \tw\t--- wireframe mode\n\
     \te\t--- hidden line mode\n\
     \tr\t--- hidden surface mode\n\
     \tu\t--- move object down\n\
     \ti\t--- move object up\n\
     \to\t--- toggle controls style\n\
     \ts\t--- toggle stereo display\n\
     \td\t--- toggle blue stereo view\n\
     \tf\t--- toggle double buffer\n\
     \th\t--- move object right\n\
     \tj\t--- move object forward\n\
     \tk\t--- move object backward\n\
     \tl\t--- move object left\n\
     \tx a\t--- rotate about x\n\
     \ty b\t--- rotate about y\n\
     \tz c\t--- rotate about z\n\
     \t1 2 3\t--- autorotate about x\n\
     \t4 5 6\t--- autorotate about y\n\
     \t7 8 9\t--- autorotate about z\n\
     \t[ ] { }\t--- adjust focus\n\n\
     HOLD the left mouse button and MOVE mouse to ROTATE object\n\n\
";


extern "C" {
   Window_t x3d_main(Float_t *longitude, Float_t *latitude, Float_t *psi,
                     Option_t *option, Window_t parent);
   void     x3d_set_display(Display_t display);
   int      x3d_dispatch_event(Handle_t event);
   void     x3d_update();
   void     x3d_get_position(Float_t *longitude, Float_t *latitude, Float_t *psi);
   int      x3d_exec_command(int px, int py, char command);
   void     x3d_terminate();
}


// Canvas menu command ids
enum EX3DViewerCommands {
   kFileNewViewer,
   kFileSave,
   kFileSaveAs,
   kFilePrint,
   kFileCloseViewer,

   kHelpAbout,
   kHelpOnViewer
};

Bool_t TViewerX3D::fgCreated = kFALSE;


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TX3DContainer                                                        //
//                                                                      //
// Utility class used by TViewerX3D. The TX3DContainer is the frame     //
// embedded in the TGCanvas widget. The X3D graphics goes into this     //
// frame. This class is used to enable input events on this graphics    //
// frame and forward the events to X3D.                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TX3DContainer : public TGCompositeFrame {
private:
   TViewerX3D  *fViewer;    // pointer back to viewer imp
public:
   TX3DContainer(TViewerX3D *c, Window_t id, const TGWindow *parent);

   Bool_t  HandleButton(Event_t *ev)
                  {  x3d_dispatch_event(gVirtualX->GetNativeEvent());
                     return fViewer->HandleContainerButton(ev); }
   Bool_t  HandleConfigureNotify(Event_t *ev)
                  {  TGFrame::HandleConfigureNotify(ev);
                     return x3d_dispatch_event(gVirtualX->GetNativeEvent()); }
   Bool_t  HandleKey(Event_t *)
                  {  return x3d_dispatch_event(gVirtualX->GetNativeEvent()); }
   Bool_t  HandleMotion(Event_t *)
                  {  return x3d_dispatch_event(gVirtualX->GetNativeEvent()); }
   Bool_t  HandleExpose(Event_t *)
                  {  return x3d_dispatch_event(gVirtualX->GetNativeEvent()); }
   Bool_t  HandleColormapChange(Event_t *)
                  {  return x3d_dispatch_event(gVirtualX->GetNativeEvent()); }
};


////////////////////////////////////////////////////////////////////////////////
/// Create a canvas container.

TX3DContainer::TX3DContainer(TViewerX3D *c, Window_t id, const TGWindow *p)
   : TGCompositeFrame(gClient, id, p)
{
   fViewer = c;
}


ClassImp(TViewerX3D);


////////////////////////////////////////////////////////////////////////////////
/// Create ROOT X3D viewer.

TViewerX3D::TViewerX3D(TVirtualPad *pad)
   : TVirtualViewer3D(),
     fCanvas(0), fContainer(0), fMenuBar(0), fFileMenu(0),
     fHelpMenu(0), fMenuBarLayout(0), fMenuBarItemLayout(0),
     fMenuBarHelpLayout(0), fCanvasLayout(0),
     fPad(pad), fBuildingScene(kFALSE), fPass(kSize)
{
   fMainFrame = new TX3DFrame(*this, gClient->GetRoot(), 800, 600);
   fOption = "x3d";
   fX3DWin = 0;
   fWidth  = 800;
   fHeight = 600;
   fXPos   = 0;
   fYPos   = 0;
   fTitle  = "x3d";
}


////////////////////////////////////////////////////////////////////////////////
/// Create ROOT X3D viewer.

TViewerX3D::TViewerX3D(TVirtualPad *pad, Option_t *option, const char *title,
                       UInt_t width, UInt_t height)
   : TVirtualViewer3D(),
     fCanvas(0), fContainer(0), fMenuBar(0), fFileMenu(0),
     fHelpMenu(0), fMenuBarLayout(0), fMenuBarItemLayout(0),
     fMenuBarHelpLayout(0), fCanvasLayout(0),
     fPad(pad), fBuildingScene(kFALSE), fPass(kSize)
{
   fMainFrame = new TX3DFrame(*this, gClient->GetRoot(), 800, 600);
   fOption = option;
   fX3DWin = 0;
   fWidth  = width;
   fHeight = height;
   fXPos   = 0;
   fYPos   = 0;
   fTitle  = title;
}


////////////////////////////////////////////////////////////////////////////////
/// Create ROOT X3D viewer.

TViewerX3D::TViewerX3D(TVirtualPad *pad, Option_t *option, const char *title,
                       Int_t x, Int_t y, UInt_t width, UInt_t height)
   : TVirtualViewer3D(),
     fCanvas(0), fContainer(0), fMenuBar(0), fFileMenu(0),
     fHelpMenu(0), fMenuBarLayout(0), fMenuBarItemLayout(0),
     fMenuBarHelpLayout(0), fCanvasLayout(0),
     fPad(pad), fBuildingScene(kFALSE), fPass(kSize)
{
   fMainFrame = new TX3DFrame(*this, gClient->GetRoot(), 800, 600);
   fOption = option;
   fX3DWin = 0;
   fWidth  = width;
   fHeight = height;
   fXPos   = x;
   fYPos   = y;
   fTitle  = title;
}


////////////////////////////////////////////////////////////////////////////////
/// Delete ROOT X3D viewer.

TViewerX3D::~TViewerX3D()
{
   if (!fPad) return;

   if (fgCreated) {
      DeleteX3DWindow();
   }
   delete fCanvasLayout;
   delete fMenuBarHelpLayout;
   delete fMenuBarItemLayout;
   delete fMenuBarLayout;
   delete fHelpMenu;
   delete fFileMenu;
   delete fMenuBar;
   delete fContainer;
   delete fCanvas;
   delete fMainFrame;
   fgCreated = kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Close X3D Viewer

void TViewerX3D::Close()
{
   assert(!fBuildingScene);
   fPad->ReleaseViewer3D();
   delete this;
}


////////////////////////////////////////////////////////////////////////////////
/// Create the actual canvas.

void TViewerX3D::CreateViewer(const char *name)
{
   // Create menus
   fFileMenu = new TGPopupMenu(fMainFrame->GetClient()->GetRoot());
   fFileMenu->AddEntry("&New Viewer",         kFileNewViewer);
   fFileMenu->AddSeparator();
   fFileMenu->AddEntry("Save",                kFileSave);
   fFileMenu->AddEntry("Save As...",          kFileSaveAs);
   fFileMenu->AddSeparator();
   fFileMenu->AddEntry("&Print...",           kFilePrint);
   fFileMenu->AddSeparator();
   fFileMenu->AddEntry("&Close Viewer",       kFileCloseViewer);

   //fFileMenu->DefaultEntry(kFileNewViewer);
   fFileMenu->DisableEntry(kFileNewViewer);
   fFileMenu->DisableEntry(kFileSave);
   fFileMenu->DisableEntry(kFileSaveAs);
   fFileMenu->DisableEntry(kFilePrint);

   fHelpMenu = new TGPopupMenu(fMainFrame->GetClient()->GetRoot());
   fHelpMenu->AddEntry("&About ROOT...",           kHelpAbout);
   fHelpMenu->AddSeparator();
   fHelpMenu->AddEntry("Help On X3D Viewer...", kHelpOnViewer);

   // This main frame will process the menu commands
   fFileMenu->Associate(fMainFrame);
   fHelpMenu->Associate(fMainFrame);

   // Create menubar layout hints
   fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 1, 1);
   fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
   fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

   // Create menubar
   fMenuBar = new TGMenuBar(fMainFrame, 1, 1, kHorizontalFrame);
   fMenuBar->AddPopup("&File",    fFileMenu,    fMenuBarItemLayout);
   fMenuBar->AddPopup("&Help",    fHelpMenu,    fMenuBarHelpLayout);

   fMainFrame->AddFrame(fMenuBar, fMenuBarLayout);

   // Create canvas and canvas container that will host the ROOT graphics
   fCanvas = new TGCanvas(fMainFrame, fMainFrame->GetWidth()+4, fMainFrame->GetHeight()+4,
                          kSunkenFrame | kDoubleBorder);
   InitX3DWindow();
   if (!fX3DWin) {
      fContainer    = 0;
      fCanvasLayout = 0;
      return;
   }
   fContainer = new TX3DContainer(this, fX3DWin, fCanvas->GetViewPort());
   fCanvas->SetContainer(fContainer);
   fCanvasLayout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
   fMainFrame->AddFrame(fCanvas, fCanvasLayout);

   // Misc

   fMainFrame->SetWindowName(name);
   fMainFrame->SetIconName(name);
   fMainFrame->SetClassHints("X3DViewer", "X3DViewer");

   fMainFrame->SetMWMHints(kMWMDecorAll, kMWMFuncAll, kMWMInputModeless);

   fMainFrame->MapSubwindows();

   // we need to use GetDefaultSize() to initialize the layout algorithm...
   fMainFrame->Resize(fMainFrame->GetDefaultSize());

   fMainFrame->MoveResize(fXPos, fYPos, fWidth, fHeight);
   fMainFrame->SetWMPosition(fXPos, fYPos);
   fgCreated = kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Setup geometry and initialize X3D.

void TViewerX3D::InitX3DWindow()
{
   TView *view = fPad->GetView();
   if (!view) {
      Error("InitX3DWindow", "view is not set");
      return;
   }

   const Float_t kPI = Float_t (TMath::Pi());

   Float_t longitude_rad = ( 90 + view->GetLongitude()) * kPI/180.0;
   Float_t  latitude_rad = (-90 + view->GetLatitude() ) * kPI/180.0;
   Float_t       psi_rad = (      view->GetPsi()      ) * kPI/180.0;

   // Call 'x3d' package
   x3d_set_display(gVirtualX->GetDisplay());
   fX3DWin = (Window_t) x3d_main(&longitude_rad, &latitude_rad, &psi_rad,
   fOption.Data(), fCanvas->GetViewPort()->GetId());
}


////////////////////////////////////////////////////////////////////////////////
/// The x3d viewer cannot rebuild a scene once created

void TViewerX3D::BeginScene()
{
   if (fgCreated) {
      return;
   }

   fBuildingScene = kTRUE;

   if (fPass == kSize) {
      gSize3D.numPoints = 0;
      gSize3D.numSegs   = 0;
      gSize3D.numPolys  = 0;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// The x3d viewer cannot rebuild a scene once created

void  TViewerX3D::EndScene()
{
   if (fgCreated) {
      return;
   }

   fBuildingScene = kFALSE;

   // Size pass done - and some points actually added
   if (gSize3D.numPoints != 0) {
      if (fPass == kSize) {
         // Allocate the X3D viewer buffer with sizes if any
         if (!AllocateX3DBuffer()) {
            Error("InitX3DWindow", "x3d buffer allocation failure");
            return;
         }

         // Enter draw pass and invoke another paint
         fPass = kDraw;
         fPad->Paint();
         fPass = kSize;
         CreateViewer(fTitle);
         Show();
      }
   } else {
      Int_t retval;
      new TGMsgBox(gClient->GetRoot(), gClient->GetRoot(),
                   "X3D Viewer", "Cannot display this content in the X3D viewer",
                   kMBIconExclamation, kMBOk, &retval);
      Close();
   }
}


////////////////////////////////////////////////////////////////////////////////
/// The x3d viewer cannot rebuild a scene once created

Int_t TViewerX3D::AddObject(const TBuffer3D & buffer, Bool_t * addChildren)
{
   if (fgCreated) {
      if (addChildren) {
         *addChildren = kFALSE;
      }
      return TBuffer3D::kNone;
   }
   else if (addChildren) {
      *addChildren = kTRUE;
   }
   // Ensure we have the required sections
   UInt_t reqSections = TBuffer3D::kCore|TBuffer3D::kRawSizes;

   // Sizing does not require actual raw tesselation information
   if (fPass == kDraw) {
      reqSections |= TBuffer3D::kRaw;
   }

   if (!buffer.SectionsValid(reqSections)) {
      return reqSections;
   }

   if (buffer.Type() == TBuffer3DTypes::kMarker) {
      PaintPolyMarker(buffer);
      return TBuffer3D::kNone;
   }

   switch(fPass) {
      case(kSize): {
         gSize3D.numPoints += buffer.NbPnts();
         gSize3D.numSegs   += buffer.NbSegs();
         gSize3D.numPolys  += buffer.NbPols();
         break;
      }
      case (kDraw): {
         X3DBuffer *x3dBuff = new X3DBuffer;
         x3dBuff->numPoints = buffer.NbPnts();
         x3dBuff->numSegs   = buffer.NbSegs();
         x3dBuff->numPolys  = buffer.NbPols();
         x3dBuff->points    = new Float_t[3*buffer.NbPnts()];
         for (UInt_t i=0; i<3*buffer.NbPnts();i++)
            x3dBuff->points[i] = (Float_t)buffer.fPnts[i];
         x3dBuff->segs      = buffer.fSegs;
         x3dBuff->polys     = buffer.fPols;
         FillX3DBuffer(x3dBuff);
         delete [] x3dBuff->points;
         delete x3dBuff;
         break;
      }
      default: {
         assert(kFALSE);
         break;
      }
   }

   return TBuffer3D::kNone;
}


////////////////////////////////////////////////////////////////////////////////
/// We don't support placed IDs - discard

Int_t TViewerX3D::AddObject(UInt_t /* placedID */, const TBuffer3D & buffer, Bool_t * addChildren)
{
   return AddObject(buffer,addChildren);
}


////////////////////////////////////////////////////////////////////////////////
/// Paint 3D PolyMarker

void TViewerX3D::PaintPolyMarker(const TBuffer3D & buffer) const
{
   if (fgCreated) {
      return;
   }
   UInt_t mode;

   if (buffer.NbPnts() > 10000) mode = 1;     // One line marker    '-'
   else if (buffer.NbPnts() > 3000) mode = 2; // Two lines marker   '+'
   else mode = 3;                           // Three lines marker '*'

   switch(fPass) {
      case(kSize): {
         gSize3D.numPoints += 2*mode*buffer.NbPnts();
         gSize3D.numSegs   += mode*buffer.NbPnts();
         break;
      }
      case (kDraw): {
         X3DBuffer *x3dBuff = new X3DBuffer;
         x3dBuff->numPoints = 2*mode*buffer.NbPnts();
         x3dBuff->numSegs   = mode*buffer.NbPnts();
         x3dBuff->numPolys  = 0;
         x3dBuff->points    = new Float_t[3*x3dBuff->numPoints];
         x3dBuff->segs      = new Int_t[3*x3dBuff->numSegs];
         x3dBuff->polys     = 0;

         Double_t delta = 0.002;

         for (UInt_t i = 0; i < buffer.NbPnts(); i++) {
            for (UInt_t j = 0; j < mode; j++) {
               for (UInt_t k = 0; k < 2; k++) {
                  delta *= -1;
                  for (UInt_t n = 0; n < 3; n++) {
                     x3dBuff->points[mode*6*i+6*j+3*k+n] =
                     buffer.fPnts[3*i+n] * (1 + (j == n ? delta : 0));
                  }
               }
            }
         }

         for (Int_t i=0; i<x3dBuff->numSegs; i++) {
            x3dBuff->segs[3*i  ] = buffer.fSegs[0];
            x3dBuff->segs[3*i+1] = 2*i;
            x3dBuff->segs[3*i+2] = 2*i+1;
         }

         FillX3DBuffer(x3dBuff);
         delete [] x3dBuff->points;
         delete [] x3dBuff->segs;
         delete x3dBuff;
         break;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// This function may be called from a script to animate an X3D picture
/// px, py  mouse position
///command = 0       --- move to px,py
///        = w       --- wireframe mode
///        = e       --- hidden line mode
///        = r       --- hidden surface mode
///        = u       --- move object down
///        = i       --- move object up
///        = o       --- toggle controls style
///        = s       --- toggle stereo display
///        = d       --- toggle blue stereo view
///        = f       --- toggle double buffer
///        = h       --- move object right
///        = j       --- move object forward
///        = k       --- move object backward
///        = l       --- move object left
///        = x a     --- rotate about x
///        = y b     --- rotate about y
///        = z c     --- rotate about z
///        = 1 2 3   --- autorotate about x
///        = 4 5 6   --- autorotate about y
///        = 7 8 9   --- autorotate about z
///        = [ ] { } --- adjust focus
/// Example:

Int_t TViewerX3D::ExecCommand(Int_t px, Int_t py, char command)
{
/*
{
   gSystem->Load("libX3d");
   TCanvas *c1 = new TCanvas("c1");
   TFile *f = new TFile("hsimple.root");
   TTree *ntuple = (TTree*)f->Get("ntuple");
   ntuple->SetMarkerColor(kYellow);
   ntuple->Draw("px:py:pz");
   TViewerX3D *x3d = new TViewerX3D(c1,"");
   for (Int_t i=0;i<500;i++) {
      Int_t px = i%500;
      Int_t py = (2*i)%200;
      x3d->ExecCommand(px,py,0);  //rotate
      if (i%20 >10) x3d->ExecCommand(px,py,'j'); //zoom
      if (i%20 <10) x3d->ExecCommand(px,py,'k'); //unzoom
   }
}
*/

   return x3d_exec_command(px,py,command);
}


////////////////////////////////////////////////////////////////////////////////
/// Get position

void TViewerX3D::GetPosition(Float_t &longitude, Float_t &latitude, Float_t &psi)
{
   x3d_get_position(&longitude, &latitude, &psi);
}

////////////////////////////////////////////////////////////////////////////////
/// Close X3D window.

void TViewerX3D::DeleteX3DWindow()
{
   x3d_terminate();
}


////////////////////////////////////////////////////////////////////////////////
/// Update X3D viewer.

void TViewerX3D::Update()
{
   x3d_update();
}


////////////////////////////////////////////////////////////////////////////////
/// Handle menu and other command generated by the user.

Bool_t TViewerX3D::ProcessFrameMessage(Long_t msg, Long_t parm1, Long_t)
{
   TRootHelpDialog *hd;

   switch (GET_MSG(msg)) {

      case kC_COMMAND:

         switch (GET_SUBMSG(msg)) {

            case kCM_BUTTON:
            case kCM_MENU:

               switch (parm1) {
                  // Handle File menu items...
                  case kFileNewViewer:
                     if (fPad) fPad->GetViewer3D("x3d");
                     break;
                  case kFileSave:
                  case kFileSaveAs:
                  case kFilePrint:
                     break;
                  case kFileCloseViewer:
                     fMainFrame->SendCloseMessage();
                     break;

                  // Handle Help menu items...
                  case kHelpAbout:
                     {
                        char str[32];
                        snprintf(str,32, "About ROOT %s...", gROOT->GetVersion());
                        hd = new TRootHelpDialog(fMainFrame, str, 600, 400);
                        hd->SetText(gHelpAbout);
                        hd->Popup();
                     }
                     break;
                  case kHelpOnViewer:
                     hd = new TRootHelpDialog(fMainFrame, "Help on X3D Viewer...", 600, 400);
                     hd->SetText(gHelpX3DViewer);
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
/// After button release get current position and update associated pad.

Bool_t TViewerX3D::HandleContainerButton(Event_t * /*ev */ )
{
   // Currently disabled as only drawing into one view at a time
   // Re-enalbe when multiple viewer implemented on pad
   /*if (ev->fType == kButtonRelease) {
      Float_t longitude_rad;
      Float_t latitude_rad;
      Float_t psi_rad;
      const Float_t kPI = Float_t (TMath::Pi());

      x3d_get_position(&longitude_rad, &latitude_rad, &psi_rad);

      Int_t irep;

      Float_t longitude_deg = longitude_rad * 180.0/kPI - 90;
      Float_t  latitude_deg = latitude_rad  * 180.0/kPI + 90;
      Float_t       psi_deg = psi_rad       * 180.0/kPI;

      fPad->GetView()->SetView(longitude_deg, latitude_deg, psi_deg, irep);

      fPad->SetPhi(-90 - longitude_deg);
      fPad->SetTheta(90 - latitude_deg);

      fPad->Modified(kTRUE);
      fPad->Update();
}*/
   return kTRUE;
}
