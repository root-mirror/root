// @(#)root/base:$Id$
// Author: Rene Brun   05/09/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TVirtualPS
#define ROOT_TVirtualPS


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualPS                                                           //
//                                                                      //
// Abstract interface to a PostScript driver.                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TNamed.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TAttText.h"
#include "TAttMarker.h"

class TVirtualPS : public TNamed, public TAttLine, public TAttFill, public TAttMarker, public TAttText {

private:
   TVirtualPS(const TVirtualPS&) = delete;
   TVirtualPS& operator=(const TVirtualPS&) = delete;

protected:
   Int_t        fNByte;           //Number of bytes written in the file (PDF)
   Int_t        fLenBuffer;       //Buffer length
   Int_t        fSizBuffer;       //Buffer size
   Bool_t       fPrinted;         //True when a page must be printed
   std::ofstream    *fStream;          //File stream identifier
   char        *fBuffer;          //File buffer
   const char  *fImplicitCREsc;   //Escape symbol before enforced new line

public:
   TVirtualPS();
   TVirtualPS(const char *filename, Int_t type=-111);
   virtual     ~TVirtualPS();
   virtual void  CellArrayBegin(Int_t W, Int_t H, Double_t x1, Double_t x2, Double_t y1, Double_t y2) = 0;
   virtual void  CellArrayFill(Int_t r, Int_t g, Int_t b) = 0;
   virtual void  CellArrayEnd() = 0;
   virtual void  Close(Option_t *opt="") = 0;
   virtual void  DrawBox(Double_t x1, Double_t y1,Double_t x2, Double_t  y2) = 0;
   virtual void  DrawFrame(Double_t xl, Double_t yl, Double_t xt, Double_t  yt,
                           Int_t mode, Int_t border, Int_t dark, Int_t light) = 0;
   virtual void  DrawPolyMarker(Int_t n, Float_t *x, Float_t *y) = 0;
   virtual void  DrawPolyMarker(Int_t n, Double_t *x, Double_t *y) = 0;
   virtual void  DrawPS(Int_t n, Float_t *xw, Float_t *yw) = 0;
   virtual void  DrawPS(Int_t n, Double_t *xw, Double_t *yw) = 0;
   virtual void  NewPage() = 0;
   virtual void  Open(const char *filename, Int_t type=-111) = 0;
   virtual void  Text(Double_t x, Double_t y, const char *string) = 0;
   virtual void  Text(Double_t x, Double_t y, const wchar_t *string) = 0;
   virtual void  SetColor(Float_t r, Float_t g, Float_t b) = 0;

   virtual void  PrintFast(Int_t nch, const char *string="");
   virtual void  PrintStr(const char *string="");
   virtual void  WriteInteger(Int_t i, Bool_t space=kTRUE);
   virtual void  WriteReal(Float_t r, Bool_t space=kTRUE);
   virtual void  PrintRaw(Int_t len, const char *str);
   virtual void *GetStream() const {  return (void*)fStream; }
   virtual void  SetStream(std::ofstream *os) {  fStream = os; }

   virtual void  SetType(Int_t /*type*/ = -111) { }
   virtual Int_t GetType() const { return 111; }

   ClassDef(TVirtualPS,0)  //Abstract interface to a PostScript driver
};


R__EXTERN TVirtualPS  *gVirtualPS;

#endif
