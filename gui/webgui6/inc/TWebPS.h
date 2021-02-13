// Author:  Sergey Linev, GSI  23/10/2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TWebPS
#define ROOT_TWebPS

#include "TVirtualPS.h"

#include "TWebPainting.h"

#include <memory>
#include <string>

class TWebPS : public TVirtualPS {

   std::unique_ptr<TWebPainting> fPainting;    ///!< object to store all painting

   enum EAttrKinds { attrLine = 0x1, attrFill = 0x2, attrMarker = 0x4, attrText = 0x8 };

   Float_t *StoreOperation(const std::string &oper, unsigned attrkind, int opersize = 0);

public:
   TWebPS();

   Bool_t IsEmptyPainting() const { return fPainting ? fPainting->IsEmpty() : kTRUE; }
   TWebPainting *GetPainting() { return fPainting.get(); }
   TWebPainting *TakePainting()
   {
      fPainting->FixSize();
      return fPainting.release();
   }
   void CreatePainting();

   /// not yet implemented

   void CellArrayBegin(Int_t, Int_t, Double_t, Double_t, Double_t, Double_t) override  {}
   void CellArrayFill(Int_t, Int_t, Int_t)  override {}
   void CellArrayEnd()  override  {}
   void Close(Option_t * = "")  override  {}
   void DrawFrame(Double_t, Double_t, Double_t, Double_t, Int_t, Int_t, Int_t, Int_t) override {}
   void NewPage() override {}
   void Open(const char *, Int_t = -111) override {}
   void SetColor(Float_t, Float_t, Float_t) override {}

   // overwritten methods
   void DrawBox(Double_t x1, Double_t y1, Double_t x2, Double_t y2) override;
   void DrawPolyMarker(Int_t n, Float_t *x, Float_t *y) override;
   void DrawPolyMarker(Int_t n, Double_t *x, Double_t *y) override;
   void DrawPS(Int_t n, Float_t *xw, Float_t *yw) override;
   void DrawPS(Int_t n, Double_t *xw, Double_t *yw) override;
   void Text(Double_t x, Double_t y, const char *str) override;
   void Text(Double_t x, Double_t y, const wchar_t *str) override;

private:
   //Let's make this clear:
   TWebPS(const TWebPS &rhs) = delete;
   TWebPS(TWebPS && rhs) = delete;
   TWebPS & operator = (const TWebPS &rhs) = delete;
   TWebPS & operator = (TWebPS && rhs) = delete;

   ClassDefOverride(TWebPS, 0) // Redirection of VirtualPS to Web painter
};

#endif
