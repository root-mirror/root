// @(#)root/gui:$Id$
// Author: Fons Rademakers   15/07/98

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRootEmbeddedCanvas                                                  //
//                                                                      //
// This class creates a TGCanvas in which a TCanvas is created. Use     //
// GetCanvas() to get a pointer to the TCanvas.                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TPluginManager.h"
#include "TVirtualGL.h"
#include "TGDNDManager.h"
#include "TBufferFile.h"
#include "TImage.h"
#include "TClass.h"
#include "TUrl.h"
#include "TVirtualX.h"

#include <iostream>

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRootEmbeddedContainer                                               //
//                                                                      //
// Utility class used by TRootEmbeddedCanvas. The TRootEmbeddedContainer//
// is the frame embedded in the TGCanvas widget. The ROOT graphics goes //
// into this frame. This class is used to enable input events on this   //
// graphics frame and forward the events to the TRootEmbeddedCanvas     //
// handlers.                                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TRootEmbeddedContainer : public TGCompositeFrame {
private:
   TRootEmbeddedCanvas  *fCanvas;    // pointer back to embedded canvas
public:
   TRootEmbeddedContainer(TRootEmbeddedCanvas *c, Window_t id, const TGWindow *parent);

   Bool_t  HandleButton(Event_t *ev)
               { return fCanvas->HandleContainerButton(ev); }
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
   void    SetEditable(Bool_t) { }
};

////////////////////////////////////////////////////////////////////////////////
/// Create a canvas container.

TRootEmbeddedContainer::TRootEmbeddedContainer(TRootEmbeddedCanvas *c, Window_t id,
   const TGWindow *p) : TGCompositeFrame(gClient, id, p)
{
   fCanvas = c;

   gVirtualX->GrabButton(fId, kAnyButton, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask |
                         kPointerMotionMask, kNone, kNone);

   AddInput(kKeyPressMask | kKeyReleaseMask | kPointerMotionMask |
            kExposureMask | kStructureNotifyMask | kLeaveWindowMask);

   fEditDisabled = kEditDisableGrab;
}




ClassImp(TRootEmbeddedCanvas);

////////////////////////////////////////////////////////////////////////////////
/// Create an TCanvas embedded in a TGFrame. A pointer to the TCanvas can
/// be obtained via the GetCanvas() member function. To embed a canvas
/// derived from a TCanvas do the following:
/// TRootEmbeddedCanvas *embedded = new TRootEmbeddedCanvas(0, p, w, h);
///      [note name must be 0, not null string ""]
/// Int_t wid = embedded->GetCanvasWindowId();
/// TMyCanvas *myc = new TMyCanvas("myname", 10, 10, wid);
/// embedded->AdoptCanvas(myc);
///      [ the MyCanvas is adopted by the embedded canvas and will be
///        destroyed by it ]

TRootEmbeddedCanvas::TRootEmbeddedCanvas(const char *name, const TGWindow *p,
            UInt_t w, UInt_t h, UInt_t options, ULong_t back)
   : TGCanvas(p, w, h, options, back)
{
   fCanvas  = 0;
   fButton  = 0;
   fAutoFit = kTRUE;
   fEditDisabled = kEditDisableLayout;

   fCWinId = -1;

   if (gStyle->GetCanvasPreferGL()) {
      //first, initialize GL (if not yet)
      if (!gGLManager) {
         TString x = "win32";
         if (gVirtualX->InheritsFrom("TGX11"))
            x = "x11";
         else if (gVirtualX->InheritsFrom("TGCocoa"))
            x = "osx";

         TPluginHandler *ph = gROOT->GetPluginManager()->FindHandler("TGLManager", x);

         if (ph && ph->LoadPlugin() != -1) {
            if (!ph->ExecPlugin(0)) {
               Warning("CreateCanvas",
                       "Cannot load GL, will use default canvas imp instead\n");
            }
         }
      }

      if (gGLManager)
         fCWinId = gGLManager->InitGLWindow((ULong_t)GetViewPort()->GetId());
      //Context creation deferred till TCanvas creation (since there is no way to pass it to TCanvas).

      if (!gGLManager || fCWinId == -1)
         gStyle->SetCanvasPreferGL(kFALSE);//TCanvas should not use gl.
   }
   if (fCWinId == -1)
      fCWinId = gVirtualX->InitWindow((ULong_t)GetViewPort()->GetId());

   Window_t win = gVirtualX->GetWindowID(fCWinId);
   fCanvasContainer = new TRootEmbeddedContainer(this, win, GetViewPort());
   SetContainer(fCanvasContainer);

   TString cname;
   if (name) cname = name;
   else cname = TString::Format("%s_canvas", GetName());
   fCanvas = new TCanvas(cname.Data(), w, h, fCWinId);

   // define DND types
   fDNDTypeList = new Atom_t[3];
   fDNDTypeList[0] = gVirtualX->InternAtom("application/root", kFALSE);
   fDNDTypeList[1] = gVirtualX->InternAtom("text/uri-list", kFALSE);
   fDNDTypeList[2] = 0;
   gVirtualX->SetDNDAware(fId, fDNDTypeList);
   SetDNDTarget(kTRUE);

   if (!p) {
      fCanvas->SetBorderMode(0);
      MapSubwindows();
      Resize(100, 100);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete embedded ROOT canvas.

TRootEmbeddedCanvas::~TRootEmbeddedCanvas()
{
   if (!MustCleanup()) {
      delete fCanvas;
      delete fCanvasContainer;
   }
   delete [] fDNDTypeList;
}

////////////////////////////////////////////////////////////////////////////////
/// Canvas c is adopted from this embedded canvas.

void TRootEmbeddedCanvas::AdoptCanvas(TCanvas *c)
{
   if(c == 0) return;
   c->EmbedInto(fCWinId, fWidth, fHeight);
   fCanvas = c;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button events in the canvas container.

Bool_t TRootEmbeddedCanvas::HandleContainerButton(Event_t *event)
{
   if (!fCanvas) return kTRUE;

   Int_t button = event->fCode;
   Int_t x = event->fX;
   Int_t y = event->fY;

   if (event->fType == kButtonPress) {
      fButton = button;
      if (button == kButton1) {
         if (event->fState & kKeyShiftMask)
            fCanvas->HandleInput(EEventType(7), x, y);
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
      if (button == kButton1)
         fCanvas->HandleInput(kButton1Up, x, y);
      if (button == kButton2)
         fCanvas->HandleInput(kButton2Up, x, y);
      if (button == kButton3)
         fCanvas->HandleInput(kButton3Up, x, y);
      if (button == kButton4)
         fCanvas->HandleInput(EEventType(5), x, y);//hack
      if (button == kButton5)
         fCanvas->HandleInput(EEventType(6), x, y);//hack

      fButton = 0;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button double click events in the canvas container.

Bool_t TRootEmbeddedCanvas::HandleContainerDoubleClick(Event_t *event)
{
   if (!fCanvas) return kTRUE;

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

Bool_t TRootEmbeddedCanvas::HandleContainerConfigure(Event_t *)
{
   if (fAutoFit && fCanvas) {
      fCanvas->Resize();
      fCanvas->Update();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle keyboard events in the canvas container.

Bool_t TRootEmbeddedCanvas::HandleContainerKey(Event_t *event)
{
   static EGEventType previous_event = kOtherEvent;
   static UInt_t previous_keysym = 0;

   if (!fCanvas) return kTRUE;

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

Bool_t TRootEmbeddedCanvas::HandleContainerMotion(Event_t *event)
{
   if (!fCanvas) return kTRUE;

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

Bool_t TRootEmbeddedCanvas::HandleContainerExpose(Event_t *event)
{
   if (!fCanvas) return kTRUE;

   if (event->fCount == 0)
      fCanvas->Flush();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle enter/leave events. Only leave is activated at the moment.

Bool_t TRootEmbeddedCanvas::HandleContainerCrossing(Event_t *event)
{
   if (!fCanvas) return kTRUE;

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

Bool_t TRootEmbeddedCanvas::HandleDNDDrop(TDNDData *data)
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

Atom_t TRootEmbeddedCanvas::HandleDNDPosition(Int_t /*x*/, Int_t /*y*/, Atom_t action,
                                              Int_t xroot, Int_t yroot)
{
   Int_t    px = 0, py = 0;
   Window_t wtarget;

   gVirtualX->TranslateCoordinates(gClient->GetDefaultRoot()->GetId(),
                                   gVirtualX->GetWindowID(fCanvas->GetCanvasID()),
                                   xroot, yroot, px, py, wtarget);

   TPad *pad = fCanvas->Pick(px, py, 0);
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

Atom_t TRootEmbeddedCanvas::HandleDNDEnter(Atom_t *typelist)
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

Bool_t TRootEmbeddedCanvas::HandleDNDLeave()
{
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Save an embedded canvas as a C++ statement(s) on output stream out.

void TRootEmbeddedCanvas::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (!GetCanvas()) return;

   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   char quote ='"';

   out << std::endl << "   // embedded canvas" << std::endl;
   out << "   TRootEmbeddedCanvas *";
   out << GetName() << " = new TRootEmbeddedCanvas(0" << "," << fParent->GetName()
       << "," << GetWidth() << "," << GetHeight();

   if (fBackground == GetDefaultFrameBackground()) {
      if (GetOptions() == (kSunkenFrame | kDoubleBorder)) {
         out <<");" << std::endl;
      } else {
         out << "," << GetOptionString() <<");" << std::endl;
      }
   } else {
      out << "," << GetOptionString() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   out << "   Int_t w" << GetName() << " = " << GetName()
       << "->GetCanvasWindowId();" << std::endl;

   static int n = 123;
   TString cname = TString::Format("c%d", n);

   out << "   TCanvas *";
   out <<  cname << " = new TCanvas(";
   out << quote << cname.Data() << quote << ", 10, 10, w"
       << GetName() << ");" << std::endl;
   out << "   " << GetName() << "->AdoptCanvas(" << cname
       << ");" << std::endl;

   n++;
   //Next line is a connection to TCanvas::SavePrimitives()
   //GetCanvas()->SavePrimitive(out,option);
}
