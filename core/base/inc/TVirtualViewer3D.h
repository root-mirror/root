// @(#)root/base:$Id$
// Author: Olivier Couet 05/10/2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TVirtualViewer3D
#define ROOT_TVirtualViewer3D

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualViewer3D                                                     //
//                                                                      //
// Abstract 3D shapes viewer. The concrete implementations are:         //
//                                                                      //
// TViewerX3D   : X3d viewer                                            //
// TViewerOpenGL: OpenGL viewer                                         //
// TViewerPad3D : visualise the 3D scene in the current Pad             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TBuffer3D;
class TVirtualPad;
class TGLRect;

class TVirtualViewer3D : public TObject
{
public:
   virtual ~TVirtualViewer3D();

   // Viewers must always handle master (absolute) positions - and
   // buffer producers must be able to supply them. Some viewers may
   // prefer local frame & translation - and producers can optionally
   // supply them
   virtual Bool_t PreferLocalFrame() const = 0;

   // Viewers can implement their own loop over pad's primitive list.
   virtual Bool_t CanLoopOnPrimitives() const { return kFALSE; }
   // When they can, TPad::Paint() and TPad::PaintModified() simply
   // call the following function:
   virtual void   PadPaint(TVirtualPad*) {}
   virtual void   ObjectPaint(TObject*, Option_t* = "")  {}

   // Addition/removal of objects must occur between Begin/EndUpdate calls
   virtual void   BeginScene() = 0;
   virtual Bool_t BuildingScene() const = 0;
   virtual void   EndScene() = 0;

   // Simple object addition - buffer represents a unique single positioned object
   virtual Int_t  AddObject(const TBuffer3D & buffer, Bool_t * addChildren = 0) = 0;

   // Complex object addition - for adding physical objects which have common logical
   // shapes. In this case buffer describes template shape (aside from kCore).
   virtual Int_t  AddObject(UInt_t physicalID, const TBuffer3D & buffer, Bool_t * addChildren = 0) = 0;

   virtual Bool_t OpenComposite(const TBuffer3D & buffer, Bool_t * addChildren = 0) = 0;
   virtual void   CloseComposite() = 0;
   virtual void   AddCompositeOp(UInt_t operation) = 0;

   virtual TObject *SelectObject(Int_t, Int_t){return 0;}
   virtual void   DrawViewer(){}

   virtual void PrintObjects(){}
   virtual void ResetCameras(){}
   virtual void ResetCamerasAfterNextUpdate(){}

   static  TVirtualViewer3D *Viewer3D(TVirtualPad *pad = 0, Option_t *type = "");

   ClassDef(TVirtualViewer3D,0) // Abstract interface to 3D viewers
};

#endif
