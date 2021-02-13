// @(#)root/graf:$Id$
// Author: Rene Brun   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBox
#define ROOT_TBox

#include "TObject.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TAttBBox2D.h"

class TPoint;

class TBox : public TObject, public TAttLine, public TAttFill, public TAttBBox2D {

private:
   TObject     *fTip;          ///<! tool tip associated with box

protected:
   Double_t     fX1;           ///< X of 1st point
   Double_t     fY1;           ///< Y of 1st point
   Double_t     fX2;           ///< X of 2nd point
   Double_t     fY2;           ///< Y of 2nd point
   Bool_t       fResizing;     ///<! True if box is being resized

public:
   // Private bits, clients can only test but not change them
   enum {
      kCannotMove    = BIT(12)  //if set the box cannot be moved/resized
   };
   TBox();
   TBox(Double_t x1, Double_t y1,Double_t x2, Double_t  y2);
   TBox(const TBox &box);
   TBox& operator=(const TBox&);
   virtual ~TBox();
   void Copy(TObject &box) const;
   virtual Int_t DistancetoPrimitive(Int_t px, Int_t py);
   virtual void  Draw(Option_t *option="");
   virtual TBox *DrawBox(Double_t x1, Double_t y1, Double_t x2, Double_t  y2);
   virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);
   Bool_t        IsBeingResized() const { return fResizing; }
   Double_t      GetX1() const { return fX1; }
   Double_t      GetX2() const { return fX2; }
   Double_t      GetY1() const { return fY1; }
   Double_t      GetY2() const { return fY2; }
   virtual void  HideToolTip(Int_t event);
   virtual Int_t IsInside(Double_t x, Double_t y) const;
   virtual void  ls(Option_t *option="") const;
   virtual void  Paint(Option_t *option="");
   virtual void  PaintBox(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="");
   virtual void  Print(Option_t *option="") const;
   virtual void  SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void  SetX1(Double_t x1) {fX1=x1;}
   virtual void  SetX2(Double_t x2) {fX2=x2;}
   virtual void  SetY1(Double_t y1) {fY1=y1;}
   virtual void  SetY2(Double_t y2) {fY2=y2;}
   virtual void  SetToolTipText(const char *text, Long_t delayms = 1000);
   virtual Rectangle_t  GetBBox();
   virtual TPoint       GetBBoxCenter();
   virtual void         SetBBoxCenter(const TPoint &p);
   virtual void         SetBBoxCenterX(const Int_t x);
   virtual void         SetBBoxCenterY(const Int_t y);
   virtual void         SetBBoxX1(const Int_t x);
   virtual void         SetBBoxX2(const Int_t x);
   virtual void         SetBBoxY1(const Int_t y);
   virtual void         SetBBoxY2(const Int_t y);

   ClassDef(TBox,3)  //Box class
};

#endif

