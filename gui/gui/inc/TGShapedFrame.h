// @(#)root/gui:$Id$
// Author: Bertrand Bellenot 23/01/2008

/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGShapedFrame
#define ROOT_TGShapedFrame

#include "TGFrame.h"

#include "TImage.h"

#include "TGPicture.h"

class TGShapedFrame : public TGCompositeFrame {

private:
   TGShapedFrame(const TGShapedFrame&) = delete;
   TGShapedFrame& operator=(const TGShapedFrame&) = delete;

protected:
   const TGPicture      *fBgnd;     // picture used as background/shape
   TImage               *fImage;    // image used as background/shape
   virtual void          DoRedraw() {}

public:
   TGShapedFrame(const char *fname=0, const TGWindow *p=0, UInt_t w=1, UInt_t h=1, UInt_t options=0);
   virtual ~TGShapedFrame();

   const TGPicture   GetPicture() const { return *fBgnd; }
   TImage            GetImage() const { return *fImage; }

   virtual void      SavePrimitive(std::ostream &out, Option_t *option = "");

   ClassDef(TGShapedFrame, 0) // Shaped composite frame
};

#endif
