// @(#)root/graf:$Id$
// Author: Rene Brun   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TLine
#define ROOT_TLine


#include "TObject.h"
#include "TAttLine.h"
#include "TAttBBox2D.h"

class TPoint;

class TLine : public TObject, public TAttLine, public TAttBBox2D {

protected:
   Double_t      fX1{0};           ///< X of 1st point
   Double_t      fY1{0};           ///< Y of 1st point
   Double_t      fX2{0};           ///< X of 2nd point
   Double_t      fY2{0};           ///< Y of 2nd point

public:
   // TLine status bits
   enum {
      kLineNDC    = BIT(14), ///< Use NDC coordinates
      kVertical   = BIT(15), ///< Line is vertical
      kHorizontal = BIT(16)  ///< Line is horizontal
   };

   TLine() {}
   TLine(Double_t x1, Double_t y1, Double_t x2, Double_t  y2);
   TLine(const TLine &line);
   virtual ~TLine() = default;

   TLine &operator=(const TLine &src);

   void                 Copy(TObject &line) const;
   virtual Int_t        DistancetoPrimitive(Int_t px, Int_t py);
   virtual TLine       *DrawLine(Double_t x1, Double_t y1,Double_t x2, Double_t y2);
   virtual TLine       *DrawLineNDC(Double_t x1, Double_t y1,Double_t x2, Double_t y2);
   virtual void         ExecuteEvent(Int_t event, Int_t px, Int_t py);
   Double_t             GetX1() const {return fX1;}
   Double_t             GetX2() const {return fX2;}
   Double_t             GetY1() const {return fY1;}
   Double_t             GetY2() const {return fY2;}
   Bool_t               IsHorizontal();
   Bool_t               IsVertical();
   virtual void         ls(Option_t *option="") const;
   virtual void         Paint(Option_t *option="");
   virtual void         PaintLine(Double_t x1, Double_t y1,Double_t x2, Double_t  y2);
   virtual void         PaintLineNDC(Double_t u1, Double_t v1,Double_t u2, Double_t  v2);
   virtual void         Print(Option_t *option="") const;
   virtual void         SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void         SetNDC(Bool_t isNDC=kTRUE);
   void                 SetHorizontal(Bool_t set = kTRUE); // *TOGGLE* *GETTER=IsHorizontal
   void                 SetVertical(Bool_t set = kTRUE); // *TOGGLE* *GETTER=IsVertical
   virtual void         SetX1(Double_t x1) {fX1=x1;}
   virtual void         SetX2(Double_t x2) {fX2=x2;}
   virtual void         SetY1(Double_t y1) {fY1=y1;}
   virtual void         SetY2(Double_t y2) {fY2=y2;}
   virtual Rectangle_t  GetBBox();
   virtual TPoint       GetBBoxCenter();
   virtual void         SetBBoxCenter(const TPoint &p);
   virtual void         SetBBoxCenterX(const Int_t x);
   virtual void         SetBBoxCenterY(const Int_t y);
   virtual void         SetBBoxX1(const Int_t x);
   virtual void         SetBBoxX2(const Int_t x);
   virtual void         SetBBoxY1(const Int_t y);
   virtual void         SetBBoxY2(const Int_t y);

   ClassDef(TLine,3)  //A line segment
};

#endif
