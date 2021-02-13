// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveDigitSetGL
#define ROOT_TEveDigitSetGL

#include "TGLObject.h"
#include "TEveDigitSet.h"

#include <set>

class TGLViewer;
class TGLScene;

class TEveDigitSet;

class TEveDigitSetGL : public TGLObject
{
private:
   TEveDigitSetGL(const TEveDigitSetGL&);            // Not implemented
   TEveDigitSetGL& operator=(const TEveDigitSetGL&); // Not implemented

protected:
   mutable const std::set<Int_t> *fHighlightSet;

   Bool_t SetupColor(const TEveDigitSet::DigitBase_t& q) const;
   void   DrawFrameIfNeeded(TGLRnrCtx& rnrCtx) const;

public:
   TEveDigitSetGL();
   virtual ~TEveDigitSetGL() {}

   virtual void   SetBBox();

   virtual void   DrawHighlight(TGLRnrCtx& rnrCtx, const TGLPhysicalShape* pshp, Int_t lvl=-1) const;

   virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
   virtual Bool_t AlwaysSecondarySelect()   const { return ((TEveDigitSet*)fExternalObj)->GetAlwaysSecSelect(); }
   virtual void   ProcessSelection(TGLRnrCtx& rnrCtx, TGLSelectRecord& rec);

   ClassDef(TEveDigitSetGL, 0); // GL renderer class for TEveDigitSet.
};

#endif
