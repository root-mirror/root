// @(#)root/graf:$Id$
// Author: Rene Brun   12/05/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMarker
#define ROOT_TMarker


#include "TObject.h"
#include "TAttMarker.h"
#include "TAttBBox2D.h"

class TPoint;

class TMarker : public TObject, public TAttMarker, public TAttBBox2D {

protected:
   Double_t     fX;           ///< X position of marker (left,center,etc..)
   Double_t     fY;           ///< Y position of marker (left,center,etc..)

public:
   // TMarker status bits
   enum {
      kMarkerNDC = BIT(14)  ///< Marker position is in NDC
   };

   TMarker();
   TMarker(Double_t x, Double_t y, Int_t marker);
   TMarker(const TMarker &marker);
   virtual ~TMarker();

   void             Copy(TObject &marker) const;
   virtual Int_t    DistancetoPrimitive(Int_t px, Int_t py);
   virtual void     Draw(Option_t *option="");
   virtual void     DrawMarker(Double_t x, Double_t y);
   virtual void     ExecuteEvent(Int_t event, Int_t px, Int_t py);
   Double_t         GetX() const  {return fX;}
   Double_t         GetY() const  {return fY;}
   virtual void     ls(Option_t *option="") const;
   virtual void     Paint(Option_t *option="");
   virtual void     PaintMarker(Double_t x, Double_t y);
   virtual void     PaintMarkerNDC(Double_t u, Double_t v);
   virtual void     Print(Option_t *option="") const;
   virtual void     SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void     SetNDC(Bool_t isNDC=kTRUE);
   virtual void     SetX(Double_t x) { fX = x;} // *MENU*
   virtual void     SetY(Double_t y) { fY = y;} // *MENU*

   virtual Rectangle_t  GetBBox();
   virtual TPoint       GetBBoxCenter();
   virtual void         SetBBoxCenter(const TPoint &p);
   virtual void         SetBBoxCenterX(const Int_t x);
   virtual void         SetBBoxCenterY(const Int_t y);
   virtual void         SetBBoxX1(const Int_t x);
   virtual void         SetBBoxX2(const Int_t x);
   virtual void         SetBBoxY1(const Int_t y);
   virtual void         SetBBoxY2(const Int_t y);

   static  void     DisplayMarkerTypes();
   static  void     DisplayMarkerLineWidths();

   ClassDef(TMarker,3)  //Marker
};

#endif

