// @(#)root/graf:$Id$
// Author: Rene Brun   16/10/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEllipse
#define ROOT_TEllipse


#include "TObject.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TAttBBox2D.h"

class TPoint;

class TEllipse : public TObject, public TAttLine, public TAttFill, public TAttBBox2D {

protected:
   Double_t    fX1;        ///< X coordinate of centre
   Double_t    fY1;        ///< Y coordinate of centre
   Double_t    fR1;        ///< first radius
   Double_t    fR2;        ///< second radius
   Double_t    fPhimin;    ///< Minimum angle (degrees)
   Double_t    fPhimax;    ///< Maximum angle (degrees)
   Double_t    fTheta;     ///< Rotation angle (degrees)

public:
   // TEllipse status bits
   enum {
      kNoEdges     = BIT(9)   // don't draw lines connecting center to edges
   };
   TEllipse();
   TEllipse(Double_t x1, Double_t y1,Double_t r1,Double_t r2=0,Double_t phimin=0, Double_t phimax=360,Double_t theta=0);
   TEllipse(const TEllipse &ellipse);
   virtual ~TEllipse();
   void   Copy(TObject &ellipse) const;
   virtual Int_t        DistancetoPrimitive(Int_t px, Int_t py);
   virtual void         Draw(Option_t *option="");
   virtual void         DrawEllipse(Double_t x1, Double_t y1, Double_t r1,Double_t r2,Double_t phimin, Double_t phimax,Double_t theta,Option_t *option="");
   virtual void         ExecuteEvent(Int_t event, Int_t px, Int_t py);
   Double_t             GetX1() const {return fX1;}
   Double_t             GetY1() const {return fY1;}
   Double_t             GetR1() const {return fR1;}
   Double_t             GetR2() const {return fR2;}
   Double_t             GetPhimin() const {return fPhimin;}
   Double_t             GetPhimax() const {return fPhimax;}
   Double_t             GetTheta() const  {return fTheta;}
   Bool_t               GetNoEdges() const;
   virtual void         ls(Option_t *option="") const;
   virtual void         Paint(Option_t *option="");
   virtual void         PaintEllipse(Double_t x1, Double_t y1, Double_t r1,Double_t r2,Double_t phimin, Double_t phimax,Double_t theta,Option_t *option="");
   virtual void         Print(Option_t *option="") const;
   virtual void         SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void         SetNoEdges(Bool_t noEdges=kTRUE); // *TOGGLE* *GETTER=GetNoEdges
   virtual void         SetPhimin(Double_t phi=0)   {fPhimin=phi;} // *MENU*
   virtual void         SetPhimax(Double_t phi=360) {fPhimax=phi;} // *MENU*
   virtual void         SetR1(Double_t r1) {fR1=r1;} // *MENU*
   virtual void         SetR2(Double_t r2) {fR2=r2;} // *MENU*
   virtual void         SetTheta(Double_t theta=0) {fTheta=theta;} // *MENU*
   virtual void         SetX1(Double_t x1) {fX1=x1;} // *MENU*
   virtual void         SetY1(Double_t y1) {fY1=y1;} // *MENU*
   virtual Rectangle_t  GetBBox();
   virtual TPoint       GetBBoxCenter();
   virtual void         SetBBoxCenter(const TPoint &p);
   virtual void         SetBBoxCenterX(const Int_t x);
   virtual void         SetBBoxCenterY(const Int_t y);
   virtual void         SetBBoxX1(const Int_t x);
   virtual void         SetBBoxX2(const Int_t x);
   virtual void         SetBBoxY1(const Int_t y);
   virtual void         SetBBoxY2(const Int_t y);


   ClassDef(TEllipse,3)  //An ellipse
};

#endif
