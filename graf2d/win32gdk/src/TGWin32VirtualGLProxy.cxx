// @(#)root/win32gdk:$Id$
// Author: Valeriy Onuchin   05/08/04

/*************************************************************************
 * Copyright (C) 1995-2004,Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGWin32ProxyDefs.h"
#include "TGWin32VirtualGLProxy.h"
#include "TGWin32.h"
#include "TROOT.h"
#include "TList.h"

static TGLManager *gManager = 0;

////////////////////////////////////////////////////////////////////////////////

TGLManager *TGWin32GLManagerProxy::RealObject()
{
   if (!gManager) {
      gManager = (TGLManager *)gROOT->GetListOfSpecials()->FindObject("gGLManager");
   }

   return gManager;
}

TGWin32GLManagerProxy::TGWin32GLManagerProxy()
{
   fIsVirtualX = kFALSE;
}

RETURN_PROXY_OBJECT(GLManager)
RETURN_METHOD_ARG1(GLManager, Int_t, InitGLWindow, Window_t, winID)
RETURN_METHOD_ARG1(GLManager, Int_t, CreateGLContext, Int_t, winInd)
RETURN_METHOD_ARG5(GLManager, Bool_t, AttachOffScreenDevice, Int_t, winInd, Int_t, x, Int_t, y, UInt_t, w, UInt_t, h)
RETURN_METHOD_ARG5(GLManager, Bool_t, ResizeOffScreenDevice, Int_t, pixInd, Int_t, x, Int_t, y, UInt_t, w, UInt_t, h)
VOID_METHOD_ARG1(GLManager, SelectOffScreenDevice, Int_t, pixInd, 1)
RETURN_METHOD_ARG1(GLManager, Int_t, GetVirtualXInd, Int_t, pixInd)
VOID_METHOD_ARG2(GLManager, MarkForDirectCopy, Int_t, pixInd, Bool_t, direct, 1)
RETURN_METHOD_ARG1(GLManager, Bool_t, MakeCurrent, Int_t, devInd)
VOID_METHOD_ARG1(GLManager, Flush, Int_t, devInd, 1)
VOID_METHOD_ARG1(GLManager, ReadGLBuffer, Int_t, devInd, 1)
VOID_METHOD_ARG1(GLManager, DeleteGLContext, Int_t, devInd, 1)
VOID_METHOD_ARG2(GLManager, ExtractViewport, Int_t, pixInd, Int_t *, vp, 1)
VOID_METHOD_ARG1(GLManager, PrintViewer, TVirtualViewer3D *, vv, 1)
VOID_METHOD_ARG1(GLManager, PaintSingleObject, TVirtualGLPainter *, p, 1)
VOID_METHOD_ARG3(GLManager, PanObject, TVirtualGLPainter *, p, Int_t, x, Int_t, y, 1)
RETURN_METHOD_ARG3(GLManager, Bool_t, PlotSelected, TVirtualGLPainter *, plot, Int_t, x, Int_t, y)
RETURN_METHOD_ARG3(GLManager, char *, GetPlotInfo, TVirtualGLPainter *, plot, Int_t, x, Int_t, y)
RETURN_METHOD_ARG4(GLManager, Bool_t, SelectManip, TVirtualGLManip *, manip, const TGLCamera *, camera, const TGLRect *, rect, const TGLBoundingBox *, box)
RETURN_METHOD_ARG1(GLManager, Bool_t, HighColorFormat, Int_t, ctx)
