// @(#)root/gl:$Id$
// Author:  Matevz Tadel, Feb 2007

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLManipSet
#define ROOT_TGLManipSet

#include "TGLOverlay.h"
#include "TGLPShapeRef.h"

class TGLManip;

class TGLManipSet : public TGLOverlayElement,
                    public TGLPShapeRef
{
public:
   enum EManip { kTrans, kScale, kRotate, kEndType };

private:
   TGLManipSet(const TGLManipSet&);            // Not implemented
   TGLManipSet& operator=(const TGLManipSet&); // Not implemented

protected:
   TGLManip  *fManip[3]; //! manipulator store
   EManip     fType;     //! current manipulator

   Bool_t     fDrawBBox; //! also draw bounding-box around physical

public:
   TGLManipSet();
   virtual ~TGLManipSet();

   virtual void SetPShape(TGLPhysicalShape* shape);

   virtual Bool_t MouseEnter(TGLOvlSelectRecord& selRec);
   virtual Bool_t Handle(TGLRnrCtx& rnrCtx, TGLOvlSelectRecord& selRec,
                         Event_t* event);
   virtual void   MouseLeave();

   virtual void Render(TGLRnrCtx& rnrCtx);

   TGLManip* GetCurrentManip() const { return fManip[fType]; }

   Int_t  GetManipType()   const { return fType; }
   void   SetManipType(Int_t type);
   Bool_t GetDrawBBox()    const { return fDrawBBox; }
   void   SetDrawBBox(Bool_t bb) { fDrawBBox = bb; }

   ClassDef(TGLManipSet, 0); // A collection of available manipulators.
}; // endclass TGLManipSet


#endif
