// @(#)root/gl:$Id$
// Author: Bertrand Bellenot 23/01/2008

/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGLayout.h"
#include "TGLWidget.h"
#include "TGLSAFrame.h"
#include "TString.h"
#include "TGLPShapeObj.h"
#include "TGLPhysicalShape.h"
#include "TGLLogicalShape.h"
#include "TGLEmbeddedViewer.h"
#include "TGLEventHandler.h"

/** \class TGLEmbeddedViewer
\ingroup opengl
Minimal GL-viewer that can be embedded in a standard ROOT frames.
*/

ClassImp(TGLEmbeddedViewer);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.
/// Argument 'border' specifies how many pixels to pad on each side of the
/// viewer. This area can be used for highlightning of the active viewer.

TGLEmbeddedViewer::TGLEmbeddedViewer(const TGWindow *parent, TVirtualPad *pad,
                                     Int_t border) :
   TGLViewer(pad, 0, 0, 400, 300),
   fFrame(0),
   fBorder(border)
{
   Init(parent);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor allowing to also specify an GED editor to use.
/// Argument 'border' specifies how many pixels to pad on each side of the
/// viewer. This area can be used for highlightning of the active viewer.

TGLEmbeddedViewer::TGLEmbeddedViewer(const TGWindow *parent, TVirtualPad *pad,
                                     TGedEditor *ged, Int_t border) :
   TGLViewer(pad, 0, 0, 400, 300),
   fFrame(0),
   fBorder(border)
{
   fGedEditor = ged;
   Init(parent);
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy standalone viewer object.

TGLEmbeddedViewer::~TGLEmbeddedViewer()
{
   delete fFrame;
   fGLWidget = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Common initialization from all constructors.

void TGLEmbeddedViewer::Init(const TGWindow *parent)
{
   fFrame = new TGCompositeFrame(parent);

   CreateFrames();

   fFrame->MapSubwindows();
   fFrame->Resize(fFrame->GetDefaultSize());
   fFrame->Resize(400, 300);
}

////////////////////////////////////////////////////////////////////////////////
/// Internal frames creation.

void TGLEmbeddedViewer::CreateFrames()
{
   fGLWidget = TGLWidget::Create(fFrame, kTRUE, kTRUE, 0, 10, 10);

   // Direct events from the TGWindow directly to the base viewer
   fEventHandler = new TGLEventHandler(0, this);
   fGLWidget->SetEventHandler(fEventHandler);

   fFrame->AddFrame(fGLWidget, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,
                                                 fBorder, fBorder, fBorder, fBorder));
}

////////////////////////////////////////////////////////////////////////////////
/// Create a GLwidget, it is an error if it is already created.
/// This is needed for frame-swapping on mac.

void TGLEmbeddedViewer::CreateGLWidget()
{
   if (fGLWidget) {
      Error("CreateGLWidget", "Widget already exists.");
      return;
   }

   fGLWidget = TGLWidget::Create(fFrame, kTRUE, kTRUE, 0, 10, 10);
   fGLWidget->SetEventHandler(fEventHandler);

   fFrame->AddFrame(fGLWidget, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,
                                                 fBorder, fBorder, fBorder, fBorder));
   fFrame->Layout();

   fGLWidget->MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy the GLwidget, it is an error if it does not exist.
/// This is needed for frame-swapping on mac.

void TGLEmbeddedViewer::DestroyGLWidget()
{
   if (fGLWidget == 0) {
      Error("DestroyGLWidget", "Widget does not exist.");
      return;
   }

   fGLWidget->UnmapWindow();
   fGLWidget->SetEventHandler(0);

   fFrame->RemoveFrame(fGLWidget);
   fGLWidget->DeleteWindow();
   fGLWidget = 0;
}
