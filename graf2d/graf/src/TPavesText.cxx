// @(#)root/graf:$Id$
// Author: Rene Brun   19/11/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <iostream>
#include "TROOT.h"
#include "TPavesText.h"
#include "TVirtualPad.h"
#include "TBufferFile.h"
#include "TError.h"

ClassImp(TPavesText);


/** \class TPavesText
\ingroup BasicGraphics

A PaveText (see TPaveText) with several stacked paves.

\image html graf_pavestext.png
*/

////////////////////////////////////////////////////////////////////////////////
/// Pavestext default constructor.

TPavesText::TPavesText(): TPaveText()
{
   fNpaves = 5;
}

////////////////////////////////////////////////////////////////////////////////
/// Pavestext normal constructor.
///
/// The PavesText is by default defined both bordersize=1 and option ="br".
///  - option = "T" Top frame
///  - option = "B" Bottom frame
///  - option = "R" Right frame
///  - option = "L" Left frame
///  - option = "NDC" x1,y1,x2,y2 are given in NDC
///  - option = "ARC" corners are rounded

TPavesText::TPavesText(Double_t x1, Double_t y1,Double_t x2, Double_t  y2, Int_t npaves,Option_t *option)
           :TPaveText(x1,y1,x2,y2,option)
{
   fNpaves = npaves;
   SetBorderSize(1);
}

////////////////////////////////////////////////////////////////////////////////
/// Pavestext default destructor.

TPavesText::~TPavesText()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Pavestext copy constructor.

TPavesText::TPavesText(const TPavesText &pavestext) : TPaveText()
{
   TBufferFile b(TBuffer::kWrite);
   TPavesText *p = (TPavesText*)(&pavestext);
   p->Streamer(b);
   b.SetReadMode();
   b.SetBufferOffset(0);
   Streamer(b);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this pavestext with its current attributes.

void TPavesText::Draw(Option_t *option)
{
   AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this pavestext with its current attributes.

void TPavesText::Paint(Option_t *option)
{
   // Draw the fNpaves-1 stacked paves
   // The spacing between paves is set to 3 times the bordersize
   Int_t bordersize = GetBorderSize();
   const char *opt = GetOption();
   Double_t signx, signy;
   if (strstr(opt,"l")) signx = -1;
   else                 signx =  1;
   if (strstr(opt,"b")) signy = -1;
   else                 signy =  1;
   Double_t dx = 3*signx*(gPad->PixeltoX(bordersize) - gPad->PixeltoX(0));
   Double_t dy = 3*signy*(gPad->PixeltoY(bordersize) - gPad->PixeltoY(0));

   TPave::ConvertNDCtoPad();

   for (Int_t ipave=fNpaves;ipave>1;ipave--) {
      Double_t x1 = fX1 + dx*Double_t(ipave-1);
      Double_t y1 = fY1 - dy*Double_t(ipave-1);
      Double_t x2 = fX2 + dx*Double_t(ipave-1);
      Double_t y2 = fY2 - dy*Double_t(ipave-1);
      TPave::PaintPave(x1,y1,x2,y2,bordersize,option);
   }

   // Draw the top pavetext
   TPaveText::Paint(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TPavesText::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   if (!strcmp(GetName(),"stats")) return;
   if (!strcmp(GetName(),"title")) return;
   Bool_t saved = gROOT->ClassSaved(TPavesText::Class());
   char quote = '"';
   out<<"   "<<std::endl;
   if (saved) {
      out<<"   ";
   } else {
      out<<"   TPavesText *";
   }
   out<<"pst = new TPavesText("<<fX1<<","<<fY1<<","<<fX2<<","<<fY2
      <<","<<fNpaves<<","<<quote<<fOption<<quote<<");"<<std::endl;

   if (strcmp(GetName(),"TPave")) {
      out<<"   pst->SetName("<<quote<<GetName()<<quote<<");"<<std::endl;
   }
   if (fLabel.Length() > 0) {
      out<<"   pst->SetLabel("<<quote<<fLabel<<quote<<");"<<std::endl;
   }
   if (fBorderSize != 4) {
      out<<"   pst->SetBorderSize("<<fBorderSize<<");"<<std::endl;
   }
   SaveFillAttributes(out,"pst",0,1001);
   SaveLineAttributes(out,"pst",1,1,1);
   SaveTextAttributes(out,"pst",22,0,1,62,0);
   TPaveText::SaveLines(out,"pst",saved);
   out<<"   pst->Draw();"<<std::endl;
}
