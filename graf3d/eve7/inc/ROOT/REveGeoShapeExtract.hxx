// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveGeoShapeExtract
#define ROOT7_REveGeoShapeExtract

#include "TNamed.h"

class TList;
class TGeoShape;

namespace ROOT {
namespace Experimental {

class REveGeoShapeExtract : public TNamed
{
   REveGeoShapeExtract(const REveGeoShapeExtract&) = delete;
   REveGeoShapeExtract& operator=(const REveGeoShapeExtract&) = delete;

protected:
   Double_t    fTrans[16];   // Transformation matrix, 4x4 column major.
   Float_t     fRGBA[4];     // RGBA color.
   Float_t     fRGBALine[4]; // RGBA color.
   Bool_t      fRnrSelf;     // Render this object.
   Bool_t      fRnrElements; // Render children of this object.
   Bool_t      fRnrFrame;    // Also draw shape outline.
   Bool_t      fMiniFrame;   // Minimize shape outline when drawing.
   TGeoShape*  fShape;       // Shape to be drawn for this object.
   TList*      fElements;    // Children elements.

public:
   REveGeoShapeExtract(const char *n = "REveGeoShapeExtract", const char *t = nullptr);
   ~REveGeoShapeExtract();

   Bool_t HasElements();
   void   AddElement(REveGeoShapeExtract* gse);

   void SetTrans(const Double_t arr[16]);
   void SetRGBA (const Float_t  arr[4]);
   void SetRGBALine(const Float_t  arr[4]);
   void SetRnrSelf(Bool_t r)     { fRnrSelf = r;     }
   void SetRnrElements(Bool_t r) { fRnrElements = r; }
   void SetRnrFrame(Bool_t r)    { fRnrFrame = r; }
   void SetMiniFrame(Bool_t r)   { fMiniFrame = r; }
   void SetShape(TGeoShape* s)   { fShape = s;       }
   void SetElements(TList* e)    { fElements = e;    }

   Double_t*  GetTrans()       { return fTrans; }
   Float_t*   GetRGBA()        { return fRGBA;  }
   Float_t*   GetRGBALine()    { return fRGBALine; }
   Bool_t     GetRnrSelf()     { return fRnrSelf;     }
   Bool_t     GetRnrElements() { return fRnrElements; }
   Bool_t     GetRnrFrame()    { return fRnrFrame; }
   Bool_t     GetMiniFrame()   { return fMiniFrame; }
   TGeoShape* GetShape()       { return fShape;    }
   TList*     GetElements()    { return fElements; }

   ClassDef(REveGeoShapeExtract, 1); // Globally positioned TGeoShape with rendering attributes and an optional list of daughter shape-extracts.
};

} // namespace Experimental
} // namespace ROOT

#endif
