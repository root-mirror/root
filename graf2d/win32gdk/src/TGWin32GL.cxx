// @(#)root/win32gdk:$Id$
// Author: Valeriy Onuchin(TGWin32GL)/ Timur Pocheptsov (TGWin32GLManager)

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TGWin32GL
\ingroup win32

The TGWin32GL is win32gdk implementation of TVirtualGLImp class.
*/

#include <deque>

#include "TGWin32GL.h"
#include "TGWin32VirtualGLProxy.h"
#include "TVirtualViewer3D.h"
#include "TVirtualX.h"
#include "TError.h"
#include "TROOT.h"
#include "TList.h"

#include "Windows4Root.h"
#include "gdk/gdk.h"
#include "gdk/win32/gdkwin32.h"

#include <GL/gl.h>
#include <GL/glu.h>


// Win32 GL Manager's stuff

struct TGWin32GLManager::TGLContext {
   Int_t        fWindowIndex;
   Int_t        fPixmapIndex;
   //
   HDC          fDC;
   HBITMAP      fHBitmap;
   HGLRC        fGLContext;
   //
   UInt_t       fW;
   UInt_t       fH;
   //
   Int_t        fX;
   Int_t        fY;
   //
   Bool_t       fHighColor;
   //
   Bool_t       fDirect;
   //
   UChar_t     *fDIBData;
   //
   TGLContext  *fNextFreeContext;
};

namespace {

   //RAII class for HDC, returned by CreateCompatibleDC
   class CDCGuard {
   private:
      HDC fHDC;

      CDCGuard(const CDCGuard &);
      CDCGuard &operator = (const CDCGuard &);

   public:
      explicit CDCGuard(HDC hDC) : fHDC(hDC)
      {}
      ~CDCGuard()
      {
         if (fHDC)
            DeleteDC(fHDC);
      }
      void Stop()
      {
         fHDC = 0;
      }
   };

   //RAII class for HDC, returned by GetWindowDC
   class WDCGuard {
   private:
      HDC fHDC;
      Window_t fWinID;

      WDCGuard(const WDCGuard &);
      WDCGuard &operator = (const WDCGuard &);

   public:
      WDCGuard(HDC hDC, Window_t winID) : fHDC(hDC), fWinID(winID)
      {}
      ~WDCGuard()
      {
         if (fHDC)
            ReleaseDC((HWND)GDK_DRAWABLE_XID((GdkWindow *)fWinID), fHDC);
      }
      void Stop()
      {
         fHDC = 0;
      }
   };

   //RAII class for HBITMAP
   class BMPGuard {
   private:
      HBITMAP fBMP;

      BMPGuard(const BMPGuard &);
      BMPGuard &operator = (const BMPGuard &);

   public:
      explicit BMPGuard(HBITMAP bmp) : fBMP(bmp)
      {}
      ~BMPGuard()
      {
         if (fBMP)
            DeleteObject(fBMP);
      }
      void Stop()
      {
         fBMP = 0;
      }
   };

   //RAII class for HGLRC
   class WGLGuard {
   private:
      HGLRC fCtx;

      WGLGuard(const WGLGuard &);
      WGLGuard &operator = (const WGLGuard &);

   public:
      explicit WGLGuard(HGLRC glrc) : fCtx(glrc)
      {}
      ~WGLGuard()
      {
         if (fCtx)
            wglDeleteContext(fCtx);
      }
      void Stop()
      {
         fCtx = 0;
      }
   };
}

const PIXELFORMATDESCRIPTOR
doubleBufferDesc = {
   sizeof doubleBufferDesc,        // size of this pfd
   1,                              // version number
   PFD_DRAW_TO_WINDOW |            // support window
   PFD_SUPPORT_OPENGL |            // support OpenGL
   PFD_DOUBLEBUFFER,               // double buffered
   PFD_TYPE_RGBA,                  // RGBA type
   24,                             // 24-bit color depth
   0, 0, 0, 0, 0, 0,               // color bits ignored
   0,                              // no alpha buffer
   0,                              // shift bit ignored
   0,                              // no accumulation buffer
   0, 0, 0, 0,                     // accum bits ignored
   32,                             // 32-bit z-buffer
   8,                              // stencil buffer depth
   0,                              // no auxiliary buffer
   PFD_MAIN_PLANE                  // main layer
};

const PIXELFORMATDESCRIPTOR
singleScreenDesc = {
   sizeof singleScreenDesc,        // size of this pfd
   1,                              // version number
   PFD_DRAW_TO_BITMAP |            // draw into bitmap
   PFD_SUPPORT_OPENGL,             // support OpenGL
   PFD_TYPE_RGBA,                  // RGBA type
   24,                             // 24-bit color depth
   0, 0, 0, 0, 0, 0,               // color bits ignored
   0,                              // no alpha buffer
   0,                              // shift bit ignored
   0,                              // no accumulation buffer
   0, 0, 0, 0,                     // accum bits ignored
   32,                             // 32-bit z-buffer
   8,                              // stencil buffer depth
   0,                              // no auxiliary buffer
   PFD_MAIN_PLANE                  // main layer
};

class TGWin32GLManager::TGWin32GLImpl {
public:
   TGWin32GLImpl() : fNextFreeContext(nullptr)
   {}
   ~TGWin32GLImpl();
   std::deque<TGLContext> fGLContexts;
   TGLContext *fNextFreeContext;
};

TGWin32GLManager::TGWin32GLImpl::~TGWin32GLImpl()
{
   //all devices should be destroyed at this moment
   std::deque<TGLContext>::size_type i = 0;

   for (; i < fGLContexts.size(); ++i) {
      TGLContext &ctx = fGLContexts[i];

      if (ctx.fGLContext) {
         //gl context (+DIB, if exists) must be destroyed from outside, by pad.
         ::Warning("TGWin32GLManager::~TGLWin32GLManager", "You forget to destroy gl-context %d\n", i);
         //destroy hdc and glrc, pixmap will be destroyed by TVirtualX
         if (ctx.fPixmapIndex != -1) {
            gVirtualX->SelectWindow(ctx.fPixmapIndex);
            gVirtualX->ClosePixmap();
         }

         wglDeleteContext(ctx.fGLContext);
         ReleaseDC((HWND)GDK_DRAWABLE_XID((GdkWindow *)gVirtualX->GetWindowID(ctx.fWindowIndex)),
                   ctx.fDC);
      }
   }
}

ClassImp(TGWin32GLManager);

////////////////////////////////////////////////////////////////////////////////

TGWin32GLManager::TGWin32GLManager() : fPimpl(new TGWin32GLImpl)
{
   gPtr2GLManager = &TGWin32GLManagerProxy::ProxyObject;
   gROOT->GetListOfSpecials()->AddLast(this);
   gGLManager = this;
}

////////////////////////////////////////////////////////////////////////////////

TGWin32GLManager::~TGWin32GLManager()
{
   delete fPimpl;
}

////////////////////////////////////////////////////////////////////////////////

Int_t TGWin32GLManager::InitGLWindow(Window_t winID)
{
   return gVirtualX->InitWindow(winID);
}

////////////////////////////////////////////////////////////////////////////////
///winInd is TGWin32 index, returned by previous call gGLManager->InitGLWindow
///returns descripto (index) of gl context or -1 if failed

Int_t TGWin32GLManager::CreateGLContext(Int_t winInd)
{
   Window_t winID = gVirtualX->GetWindowID(winInd);
   HDC hDC = GetWindowDC((HWND)GDK_DRAWABLE_XID((GdkWindow *)winID));

   if (!hDC) {
      Error("CreateGLContext", "GetWindowDC failed\n");
      return -1;
   }

   WDCGuard dcGuard(hDC, winID);

   if (Int_t pixFormat = ChoosePixelFormat(hDC, &doubleBufferDesc)) {
      if (SetPixelFormat(hDC, pixFormat, &doubleBufferDesc)) {
         HGLRC glCtx = wglCreateContext(hDC);

         if (!glCtx) {
            Error("CreateGLContext", "wglCreateContext failed\n");
            return -1;
         }

         TGLContext newDevice = {winInd, -1, hDC, 0, glCtx};
         PIXELFORMATDESCRIPTOR testFormat = {};
         DescribePixelFormat(hDC, pixFormat, sizeof testFormat, &testFormat);
         newDevice.fHighColor = testFormat.cColorBits < 24 ? kTRUE : kFALSE;

         if (TGLContext *ctx = fPimpl->fNextFreeContext) {
            Int_t ind = ctx->fWindowIndex;
            fPimpl->fNextFreeContext = fPimpl->fNextFreeContext->fNextFreeContext;
            *ctx = newDevice;
            dcGuard.Stop();
            return ind;
         } else {
            WGLGuard wglGuard(glCtx);
            fPimpl->fGLContexts.push_back(newDevice);
            wglGuard.Stop();
            dcGuard.Stop();
            return fPimpl->fGLContexts.size() - 1;
         }
      } else
         Error("CreateGLContext", "SetPixelFormat failed\n");
   } else
      Error("CreateGLContext", "ChoosePixelFormat failed\n");

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
///Create DIB section to read GL buffer into

Bool_t TGWin32GLManager::CreateDIB(TGLContext &ctx)const
{
   HDC dibDC = CreateCompatibleDC(0);

   if (!dibDC) {
      Error("CreateDIB", "CreateCompatibleDC failed\n");
      return kFALSE;
   }

   CDCGuard dcGuard(dibDC);

   BITMAPINFOHEADER bmpHeader = {sizeof bmpHeader, (LONG) ctx.fW, (LONG) ctx.fH, 1, 32, BI_RGB};
   void *bmpCnt = nullptr;
   HBITMAP hDIB = CreateDIBSection(dibDC, (BITMAPINFO*)&bmpHeader, DIB_RGB_COLORS, &bmpCnt, 0, 0);

   if (!hDIB) {
      Error("CreateDIB", "CreateDIBSection failed\n");
      return kFALSE;
   }

   BMPGuard bmpGuard(hDIB);

   ctx.fPixmapIndex = gVirtualX->AddPixmap((ULong_t)hDIB, ctx.fW, ctx.fH);
   ctx.fHBitmap = hDIB;
   ctx.fDIBData = static_cast<UChar_t *>(bmpCnt);

   bmpGuard.Stop();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TGWin32GLManager::AttachOffScreenDevice(Int_t ctxInd, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   TGLContext &ctx = fPimpl->fGLContexts[ctxInd];
   TGLContext newCtx = {ctx.fWindowIndex, -1, ctx.fDC, 0, ctx.fGLContext, w, h, x, y, ctx.fHighColor};

   if (CreateDIB(newCtx)) {
      ctx = newCtx;
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
///Create new DIB if needed

Bool_t TGWin32GLManager::ResizeOffScreenDevice(Int_t ctxInd, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   TGLContext &ctx = fPimpl->fGLContexts[ctxInd];

   if (ctx.fPixmapIndex != -1)
      if (TMath::Abs(Int_t(w) - Int_t(ctx.fW)) > 1 || TMath::Abs(Int_t(h) - Int_t(ctx.fH)) > 1) {
         TGLContext newCtx = {ctx.fWindowIndex, -1, ctx.fDC, 0, ctx.fGLContext, w, h, x, y, ctx.fHighColor};
         if (CreateDIB(newCtx)) {
            //new DIB created
            gVirtualX->SelectWindow(ctx.fPixmapIndex);
            gVirtualX->ClosePixmap();
            ctx = newCtx;
         } else {
            Error("ResizeOffScreenDevice", "Error trying to create new DIB\n");
            return kFALSE;
         }
      } else {
         ctx.fX = x;
         ctx.fY = y;
      }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////

void TGWin32GLManager::SelectOffScreenDevice(Int_t ctxInd)
{
   gVirtualX->SelectWindow(fPimpl->fGLContexts[ctxInd].fPixmapIndex);
}

////////////////////////////////////////////////////////////////////////////////

void TGWin32GLManager::MarkForDirectCopy(Int_t pixInd, Bool_t isDirect)
{
   if (fPimpl->fGLContexts[pixInd].fPixmapIndex != -1)
      fPimpl->fGLContexts[pixInd].fDirect = isDirect;
}

////////////////////////////////////////////////////////////////////////////////

void TGWin32GLManager::ReadGLBuffer(Int_t ctxInd)
{
   TGLContext &ctx = fPimpl->fGLContexts[ctxInd];

   if (ctx.fPixmapIndex != -1) {
      glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
      glReadBuffer(GL_BACK);
      glReadPixels(0, 0, ctx.fW, ctx.fH, GL_BGRA_EXT, GL_UNSIGNED_BYTE, ctx.fDIBData);
   }
}

////////////////////////////////////////////////////////////////////////////////

Int_t TGWin32GLManager::GetVirtualXInd(Int_t ctxInd)
{
   return fPimpl->fGLContexts[ctxInd].fPixmapIndex;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TGWin32GLManager::MakeCurrent(Int_t ctxInd)
{
   TGLContext &ctx = fPimpl->fGLContexts[ctxInd];
   return (Bool_t)wglMakeCurrent(ctx.fDC, ctx.fGLContext);
}

////////////////////////////////////////////////////////////////////////////////

void TGWin32GLManager::Flush(Int_t ctxInd)
{
   TGLContext &ctx = fPimpl->fGLContexts[ctxInd];

   if (ctx.fPixmapIndex == -1) {
      //doube-buffered OpenGL
      wglSwapLayerBuffers(ctx.fDC, WGL_SWAP_MAIN_PLANE);
   } else if (ctx.fDirect) {
      //DIB is flushed by viewer directly
      HDC hDC = CreateCompatibleDC(0);

      if (!hDC) {
         Error("Flush", "CreateCompatibleDC failed\n");
         return;
      }

      HBITMAP oldDIB = (HBITMAP)SelectObject(hDC, ctx.fHBitmap);

      if (!BitBlt(ctx.fDC, ctx.fX, ctx.fY, ctx.fW, ctx.fH, hDC, 0, 0, SRCCOPY))
         ctx.fDirect = kFALSE;

      SelectObject(hDC, oldDIB);
      DeleteDC(hDC);
   }
   //do nothing for non-direct off-screen device
}

////////////////////////////////////////////////////////////////////////////////

void TGWin32GLManager::DeleteGLContext(Int_t ctxInd)
{
   TGLContext &ctx = fPimpl->fGLContexts[ctxInd];

   if (ctx.fPixmapIndex != -1) {
      gVirtualX->SelectWindow(ctx.fPixmapIndex);
      gVirtualX->ClosePixmap();
      ctx.fPixmapIndex = -1;
   }

   wglDeleteContext(ctx.fGLContext);
   ctx.fGLContext = 0;
   ReleaseDC((HWND)GDK_DRAWABLE_XID((GdkWindow *)gVirtualX->GetWindowID(ctx.fWindowIndex)),
             ctx.fDC);
   //now, save its own index before putting into list of free devices
   ctx.fWindowIndex = ctxInd;
   ctx.fNextFreeContext = fPimpl->fNextFreeContext;
   fPimpl->fNextFreeContext = &ctx;
}

////////////////////////////////////////////////////////////////////////////////

void TGWin32GLManager::ExtractViewport(Int_t ctxInd, Int_t *viewport)
{
   TGLContext &ctx = fPimpl->fGLContexts[ctxInd];

   if (ctx.fPixmapIndex != -1) {
      viewport[0] = 0;
      viewport[1] = 0;
      viewport[2] = ctx.fW;
      viewport[3] = ctx.fH;
   }
}

////////////////////////////////////////////////////////////////////////////////

void TGWin32GLManager::PaintSingleObject(TVirtualGLPainter *p)
{
   p->Paint();
}

////////////////////////////////////////////////////////////////////////////////

void TGWin32GLManager::PrintViewer(TVirtualViewer3D *vv)
{
   vv->PrintObjects();
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TGWin32GLManager::SelectManip(TVirtualGLManip *manip, const TGLCamera * camera, const TGLRect * rect, const TGLBoundingBox * sceneBox)
{
   return manip->Select(*camera, *rect, *sceneBox);
}

////////////////////////////////////////////////////////////////////////////////

void TGWin32GLManager::PanObject(TVirtualGLPainter *o, Int_t x, Int_t y)
{
   return o->Pan(x, y);
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TGWin32GLManager::PlotSelected(TVirtualGLPainter *plot, Int_t px, Int_t py)
{
   return plot->PlotSelected(px, py);
}

////////////////////////////////////////////////////////////////////////////////

char *TGWin32GLManager::GetPlotInfo(TVirtualGLPainter *plot, Int_t px, Int_t py)
{
    return plot->GetPlotInfo(px, py);
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TGWin32GLManager::HighColorFormat(Int_t ctxInd)
{
   if (ctxInd == -1)
      return kFALSE;

   return fPimpl->fGLContexts[ctxInd].fHighColor;
}
