// @(#)root/graf:$Id$
// Author: Rene Brun   16/05/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TCutG
\ingroup BasicGraphics

Graphical cut class.

A TCutG object is a closed polygon defining a closed region in a x,y plot.
It can be created via the graphics editor option "CutG" or directly by
invoking its constructor. The first and last points should be the same.

To create a TCutG via the graphics editor, use the left button to select the
points building the contour of the cut. Click on the right button to close the
TCutG. When it is created via the graphics editor, the TCutG object is named
"CUTG". It is recommended to immediately change the name by using the context
menu item "SetName". When the graphics editor is used, the names of the
variables X,Y are automatically taken from the current pad title.

Example:

Assume a TTree object T and:
~~~ {.cpp}
    Root > T.Draw("abs(fMomemtum):fEtot")
~~~
the TCutG members fVarX, fVary will be set to:
~~~ {.cpp}
    fVarx = fEtot
    fVary = abs(fMomemtum)
~~~
A graphical cut can be used in a TTree selection expression:
~~~ {.cpp}
    Root > T.Draw("fEtot","cutg1")
~~~
where "cutg1" is the name of an existing graphical cut.

Note that, as shown in the example above, a graphical cut may be used in a
selection expression when drawing TTrees expressions of 1-d, 2-d or
3-dimensions. The expressions used in TTree::Draw can reference the variables in
the fVarX, fVarY of the graphical cut plus other variables.

When the TCutG object is created by TTree::Draw, it is added to the list of special objects in
the main TROOT object pointed by gROOT. To retrieve a pointer to this object
from the code or command line, do:
~~~ {.cpp}
    TCutG *mycutg;
    mycutg = (TCutG*)gROOT->GetListOfSpecials()->FindObject("CUTG")
    mycutg->SetName("mycutg");
~~~
When the TCutG is not created via TTree::Draw, one must set the variable names
corresponding to x,y if one wants to use the cut as input to TTree::Draw,eg
~~~ {.cpp}
    TCutG *cutg = new TCutG("mycut",6);
    cutg->SetVarX("y");
    cutg->SetVarY("x");
    cutg->SetPoint(0,-0.3586207,1.509534);
    cutg->SetPoint(1,-1.894181,-0.529661);
    cutg->SetPoint(2,0.07780173,-1.21822);
    cutg->SetPoint(3,-1.0375,-0.07944915);
    cutg->SetPoint(4,0.756681,0.1853814);
    cutg->SetPoint(5,-0.3586207,1.509534);
~~~
Example of use of a TCutG in TTree::Draw:
~~~ {.cpp}
    tree.Draw("x:y","mycutg && z>0 && sqrt(x)>1")
~~~
A Graphical cut may be drawn via TGraph::Draw. It can be edited like a normal
TGraph. Being a TGraph the drawing options and behavior relatives to graphs apply.
They are listed in the TGraphPainter description.
See in particular "Graphs in logarithmic scale".

A Graphical cut may be saved to a file via TCutG::Write.
*/

#include <cstring>
#include <iostream>

#include "TROOT.h"
#include "TCutG.h"
#include "TBuffer.h"
#include "TVirtualPad.h"
#include "TPaveText.h"
#include "TH2.h"
#include "strlcpy.h"

ClassImp(TCutG);

////////////////////////////////////////////////////////////////////////////////
/// TCutG default constructor.

TCutG::TCutG() : TGraph()
{
   fObjectX  = 0;
   fObjectY  = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// TCutG copy constructor

TCutG::TCutG(const TCutG &cutg)
      :TGraph(cutg)
{
   fVarX    = cutg.fVarX;
   fVarY    = cutg.fVarY;
   fObjectX = cutg.fObjectX ? cutg.fObjectX->Clone() : nullptr;
   fObjectY = cutg.fObjectY ? cutg.fObjectY->Clone() : nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// TCutG normal constructor.

TCutG::TCutG(const char *name, Int_t n)
      :TGraph(n)
{
   fObjectX  = 0;
   fObjectY  = 0;
   SetName(name);
   delete gROOT->GetListOfSpecials()->FindObject(name);
   gROOT->GetListOfSpecials()->Add(this);

   // Take name of cut variables from pad title if title contains ":"
   if (gPad) {
      TPaveText *ptitle = (TPaveText*)gPad->FindObject("title");
      if (!ptitle) return;
      TText *ttitle = ptitle->GetLineWith(":");
      if (!ttitle) ttitle = ptitle->GetLineWith("{");
      if (!ttitle) ttitle = ptitle->GetLine(0);
      if (!ttitle) return;
      const char *title = ttitle->GetTitle();
      Int_t nch = strlen(title);
      char *vars = new char[nch+1];
      strlcpy(vars,title,nch+1);
      char *col = strstr(vars,":");
      if (col) {
         *col = 0;
         col++;
         char *brak = strstr(col," {");
         if (brak) *brak = 0;
         fVarY = vars;
         fVarX = col;
      } else {
         char *brak = strstr(vars," {");
         if (brak) *brak = 0;
         fVarX = vars;
      }
      delete [] vars;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TCutG normal constructor.

TCutG::TCutG(const char *name, Int_t n, const Float_t *x, const Float_t *y)
      :TGraph(n,x,y)
{
   fObjectX  = 0;
   fObjectY  = 0;
   SetName(name);
   delete gROOT->GetListOfSpecials()->FindObject(name);
   gROOT->GetListOfSpecials()->Add(this);

   // Take name of cut variables from pad title if title contains ":"
   if (gPad) {
      TPaveText *ptitle = (TPaveText*)gPad->FindObject("title");
      if (!ptitle) return;
      TText *ttitle = ptitle->GetLineWith(":");
      if (!ttitle) ttitle = ptitle->GetLineWith("{");
      if (!ttitle) ttitle = ptitle->GetLine(0);
      if (!ttitle) return;
      const char *title = ttitle->GetTitle();
      Int_t nch = strlen(title);
      char *vars = new char[nch+1];
      strlcpy(vars,title,nch+1);
      char *col = strstr(vars,":");
      if (col) {
         *col = 0;
         col++;
         char *brak = strstr(col," {");
         if (brak) *brak = 0;
         fVarY = vars;
         fVarX = col;
      } else {
         char *brak = strstr(vars," {");
         if (brak) *brak = 0;
         fVarX = vars;
      }
      delete [] vars;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TCutG normal constructor.

TCutG::TCutG(const char *name, Int_t n, const Double_t *x, const Double_t *y)
      :TGraph(n,x,y)
{
   fObjectX  = 0;
   fObjectY  = 0;
   SetName(name);
   delete gROOT->GetListOfSpecials()->FindObject(name);
   gROOT->GetListOfSpecials()->Add(this);

   // Take name of cut variables from pad title if title contains ":"
   if (gPad) {
      TPaveText *ptitle = (TPaveText*)gPad->FindObject("title");
      if (!ptitle) return;
      TText *ttitle = ptitle->GetLineWith(":");
      if (!ttitle) ttitle = ptitle->GetLineWith("{");
      if (!ttitle) ttitle = ptitle->GetLine(0);
      if (!ttitle) return;
      const char *title = ttitle->GetTitle();
      Int_t nch = strlen(title);
      char *vars = new char[nch+1];
      strlcpy(vars,title,nch+1);
      char *col = strstr(vars,":");
      if (col) {
         *col = 0;
         col++;
         char *brak = strstr(col," {");
         if (brak) *brak = 0;
         fVarY = vars;
         fVarX = col;
      } else {
         char *brak = strstr(vars," {");
         if (brak) *brak = 0;
         fVarX = vars;
      }
      delete [] vars;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TCutG destructor.

TCutG::~TCutG()
{
   delete fObjectX;
   delete fObjectY;
   if ( gROOT && !gROOT->TestBit( TObject::kInvalidObject ) )
      gROOT->GetListOfSpecials()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TCutG &TCutG::operator=(const TCutG &rhs)
{
   if (this != &rhs) {
      TGraph::operator=(rhs);
      delete fObjectX;
      delete fObjectY;
      fObjectX = rhs.fObjectX ? rhs.fObjectX->Clone() : nullptr;
      fObjectY = rhs.fObjectY ? rhs.fObjectY->Clone() : nullptr;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the area inside this TCutG
/// The algorithm uses Stoke's theorem over the border of the closed polygon.
/// Just as a reminder: Stoke's theorem reduces a surface integral
/// to a line integral over the border of the surface integral.

Double_t TCutG::Area() const
{
   Double_t a = 0;
   Int_t n = GetN();
   for (Int_t i=0;i<n-1;i++) {
      a += (fX[i]-fX[i+1])*(fY[i]+fY[i+1]);
   }
   a *= 0.5;
   return a;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the center x,y of this TCutG
/// The algorithm uses Stoke's theorem over the border of the closed polygon.
/// Just as a reminder: Stoke's theorem reduces a surface integral
/// to a line integral over the border of the surface integral.

void TCutG::Center(Double_t &cx, Double_t &cy) const
{
   Int_t n = GetN();
   Double_t a  = 0;
   cx = cy = 0;
   Double_t t;
   for (Int_t i=0;i<n-1;i++) {
      t   = 2*fX[i]*fY[i] + fY[i]*fX[i+1] + fX[i]*fY[i+1] + 2*fX[i+1]*fY[i+1];
      cx += (fX[i]-fX[i+1])*t;
      cy += (-fY[i]+fY[i+1])*t;
      a  += (fX[i]-fX[i+1])*(fY[i]+fY[i+1]);
   }
   a  *= 0.5;
   cx *= 1./(6*a);
   cy *= 1./(6*a);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the integral of 2-d histogram h for all bins inside the cut
/// if option "width" is specified, the integral is the sum of
/// the bin contents multiplied by the bin width in x and in y.

Double_t TCutG::IntegralHist(TH2 *h, Option_t *option) const
{
   if (!h) return 0;
   Int_t n = GetN();
   Double_t xmin= 1e200;
   Double_t xmax = -xmin;
   Double_t ymin = xmin;
   Double_t ymax = xmax;
   for (Int_t i=0;i<n;i++) {
      if (fX[i] < xmin) xmin = fX[i];
      if (fX[i] > xmax) xmax = fX[i];
      if (fY[i] < ymin) ymin = fY[i];
      if (fY[i] > ymax) ymax = fY[i];
   }
   TAxis *xaxis = h->GetXaxis();
   TAxis *yaxis = h->GetYaxis();
   Int_t binx1 = xaxis->FindBin(xmin);
   Int_t binx2 = xaxis->FindBin(xmax);
   Int_t biny1 = yaxis->FindBin(ymin);
   Int_t biny2 = yaxis->FindBin(ymax);
   Int_t nbinsx = h->GetNbinsX();
   Stat_t integral = 0;

   // Loop on bins for which the bin center is in the cut
   TString opt = option;
   opt.ToLower();
   Bool_t width = kFALSE;
   if (opt.Contains("width")) width = kTRUE;
   Int_t bin, binx, biny;
   for (biny=biny1;biny<=biny2;biny++) {
      Double_t y = yaxis->GetBinCenter(biny);
      for (binx=binx1;binx<=binx2;binx++) {
         Double_t x = xaxis->GetBinCenter(binx);
         if (!IsInside(x,y)) continue;
         bin = binx +(nbinsx+2)*biny;
         if (width) integral += h->GetBinContent(bin)*xaxis->GetBinWidth(binx)*yaxis->GetBinWidth(biny);
         else       integral += h->GetBinContent(bin);
      }
   }
   return integral;
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out.

void TCutG::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';
   out<<"   "<<std::endl;
   if (gROOT->ClassSaved(TCutG::Class())) {
      out<<"   ";
   } else {
      out<<"   TCutG *";
   }
   out<<"cutg = new TCutG("<<quote<<GetName()<<quote<<","<<fNpoints<<");"<<std::endl;
   out<<"   cutg->SetVarX("<<quote<<GetVarX()<<quote<<");"<<std::endl;
   out<<"   cutg->SetVarY("<<quote<<GetVarY()<<quote<<");"<<std::endl;
   out<<"   cutg->SetTitle("<<quote<<GetTitle()<<quote<<");"<<std::endl;

   SaveFillAttributes(out,"cutg",0,1001);
   SaveLineAttributes(out,"cutg",1,1,1);
   SaveMarkerAttributes(out,"cutg",1,1,1);

   for (Int_t i=0;i<fNpoints;i++) {
      out<<"   cutg->SetPoint("<<i<<","<<fX[i]<<","<<fY[i]<<");"<<std::endl;
   }
   out<<"   cutg->Draw("
      <<quote<<option<<quote<<");"<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the X object (and delete the previous one if any).

void TCutG::SetObjectX(TObject *obj)
{
   delete fObjectX;
   fObjectX = obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the Y object (and delete the previous one if any).

void TCutG::SetObjectY(TObject *obj)
{
   delete fObjectY;
   fObjectY = obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Set X variable.

void TCutG::SetVarX(const char *varx)
{
   fVarX = varx;
   delete fObjectX;
   fObjectX = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set Y variable.

void TCutG::SetVarY(const char *vary)
{
   fVarY = vary;
   delete fObjectY;
   fObjectY = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TCutG.

void TCutG::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TCutG::Class(), this);
      gROOT->GetListOfSpecials()->Add(this);
   } else {
      R__b.WriteClassBuffer(TCutG::Class(), this);
   }
}
