// @(#)root/graf:$Id$
// Author: Nicolas Brun   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TText
#define ROOT_TText


#include "TAttText.h"
#include "TNamed.h"
#include "TAttBBox2D.h"

class TPoint;

class TText : public TNamed, public TAttText, public TAttBBox2D {

protected:
   Double_t     fX{0};                 ///< X position of text (left,center,etc..)
   Double_t     fY{0};                 ///< Y position of text (left,center,etc..)
   void         *fWcsTitle{nullptr};   ///<!Used by TMathText

public:
   // TText status bits
   enum {
      kTextNDC = BIT(14)  ///< The text position is in the NDC space
   };

   TText() {} // NOLINT: not allowed to use = default because of TObject::kIsOnHeap detection, see ROOT-10300
   TText(Double_t x, Double_t y, const char *text);
   TText(Double_t x, Double_t y, const wchar_t *text);
   TText(const TText &text);
   virtual ~TText();

   TText &operator=(const TText &src);

   void             Copy(TObject &text) const;
   virtual Int_t    DistancetoPrimitive(Int_t px, Int_t py);
   virtual TText   *DrawText(Double_t x, Double_t y, const char *text);
   virtual TText   *DrawText(Double_t x, Double_t y, const wchar_t *text);
   virtual TText   *DrawTextNDC(Double_t x, Double_t y, const char *text);
   virtual TText   *DrawTextNDC(Double_t x, Double_t y, const wchar_t *text);
   virtual void     ExecuteEvent(Int_t event, Int_t px, Int_t py);

   virtual void     GetControlBox(Int_t x, Int_t y, Double_t theta,
                                  Int_t cBoxX[4], Int_t cBoxY[4]);
   Double_t         GetX() const  { return fX; }
   virtual void     GetBoundingBox(UInt_t &w, UInt_t &h, Bool_t angle = kFALSE);
   virtual void     GetTextAscentDescent(UInt_t &a, UInt_t &d, const char *text) const;
   virtual void     GetTextAscentDescent(UInt_t &a, UInt_t &d, const wchar_t *text) const;
   virtual void     GetTextExtent(UInt_t &w, UInt_t &h, const char *text) const;
   virtual void     GetTextExtent(UInt_t &w, UInt_t &h, const wchar_t *text) const;
   virtual void     GetTextAdvance(UInt_t &a, const char *text, const Bool_t kern=kTRUE) const;
   const void *     GetWcsTitle(void) const;
   Double_t         GetY() const  { return fY; }

   virtual void     ls(Option_t *option="") const;
   virtual void     Paint(Option_t *option="");
   virtual void     PaintControlBox(Int_t x, Int_t y, Double_t theta);
   virtual void     PaintText(Double_t x, Double_t y, const char *text);
   virtual void     PaintText(Double_t x, Double_t y, const wchar_t *text);
   virtual void     PaintTextNDC(Double_t u, Double_t v, const char *text);
   virtual void     PaintTextNDC(Double_t u, Double_t v, const wchar_t *text);
   virtual void     Print(Option_t *option="") const;
   virtual void     SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void     SetMbTitle(const wchar_t *title=L"");
   virtual void     SetNDC(Bool_t isNDC=kTRUE);
   virtual void     SetText(Double_t x, Double_t y, const char *text) {fX=x; fY=y; SetTitle(text);} // *MENU* *ARGS={x=>fX,y=>fY,text=>fTitle}
   virtual void     SetText(Double_t x, Double_t y, const wchar_t *text) {fX=x; fY=y; SetMbTitle(text);}
   virtual void     SetX(Double_t x) { fX = x; } // *MENU*
   virtual void     SetY(Double_t y) { fY = y; } // *MENU*

   virtual Rectangle_t  GetBBox();
   virtual TPoint       GetBBoxCenter();
   virtual void         SetBBoxCenter(const TPoint &p);
   virtual void         SetBBoxCenterX(const Int_t x);
   virtual void         SetBBoxCenterY(const Int_t y);
   virtual void         SetBBoxX1(const Int_t x); //Not Implemented
   virtual void         SetBBoxX2(const Int_t x); //Not Implemented
   virtual void         SetBBoxY1(const Int_t y); //Not Implemented
   virtual void         SetBBoxY2(const Int_t y); //Not Implemented

   ClassDef(TText,3)  //Text
};

#endif
