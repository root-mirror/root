// @(#)root/graf:$Id$
// Author: Matthew.Adam.Dobbs   06/09/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <cstdio>
#include <iostream>

#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TList.h"
#include "TVirtualPad.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLegendEntry.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TH1.h"
#include "THStack.h"

ClassImp(TLegend);

/** \class TLegend
\ingroup BasicGraphics

This class displays a legend box (TPaveText) containing several legend entries.

Each legend entry is made of a reference to a ROOT object, a text label and an
option specifying which graphical attributes (marker/line/fill) should be
displayed.

The following example shows how to create a legend. In this example the legend
contains a histogram, a function and a graph. The histogram is put in the legend
using its reference pointer whereas the graph and the function are added
using their names. Note that, because `TGraph` constructors do not have the
`TGraph` name as parameter, the graph name should be specified using the
`SetName` method.

When an object is added by name, a scan is performed on the list of objects
contained in the current pad (`gPad`) and also in the possible
`TMultiGraph` and `THStack` present in the pad. If a matching
name is found, the corresponding object is added in the legend using its pointer.

Begin_Macro(source)
{
   auto c1 = new TCanvas("c1","c1",600,500);
   gStyle->SetOptStat(0);

   auto h1 = new TH1F("h1","TLegend Example",200,-10,10);
   h1->FillRandom("gaus",30000);
   h1->SetFillColor(kGreen);
   h1->SetFillStyle(3003);
   h1->Draw();

   auto f1=new TF1("f1","1000*TMath::Abs(sin(x)/x)",-10,10);
   f1->SetLineColor(kBlue);
   f1->SetLineWidth(4);
   f1->Draw("same");

   const Int_t n = 20;
   Double_t x[n], y[n], ex[n], ey[n];
   for (Int_t i=0;i<n;i++) {
      x[i]  = i*0.1;
      y[i]  = 1000*sin(x[i]+0.2);
      x[i]  = 17.8*x[i]-8.9;
      ex[i] = 1.0;
      ey[i] = 10.*i;
   }
   auto gr = new TGraphErrors(n,x,y,ex,ey);
   gr->SetName("gr");
   gr->SetLineColor(kRed);
   gr->SetLineWidth(2);
   gr->SetMarkerStyle(21);
   gr->SetMarkerSize(1.3);
   gr->SetMarkerColor(7);
   gr->Draw("P");

   auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend->AddEntry(h1,"Histogram filled with random numbers","f");
   legend->AddEntry("f1","Function abs(#frac{sin(x)}{x})","l");
   legend->AddEntry("gr","Graph with error bars","lep");
   legend->Draw();
}
End_Macro


`TLegend` inherits from `TAttText` therefore changing any
text attributes (text alignment, font, color...) on a legend will changed the
text attributes on each line.

In particular it can be interesting to change the text alignement that way. In
order to have a base-line vertical alignment instead of a centered one simply do:
~~~ {.cpp}
   legend->SetTextAlign(13);
~~~
or
~~~ {.cpp}
   legend->SetTextAlign(11);
~~~
The default value of some `TLegend` attributes can be changed using
`gStyle`. The default settings are:
~~~ {.cpp}
   SetLegendBorderSize(1);
   SetLegendFillColor(0);
   SetLegendFont(42);
   SetLegendTextSize(0.);
~~~
The global attributes change the default values for the next created legends.

Text attributes can be also changed individually on each legend entry:
~~~ {.cpp}
   TLegendEntry *le = leg->AddEntry(h1,"Histogram filled with random numbers","f");
   le->SetTextColor(kBlue);;
~~~

Note that the `TPad` class has a method to build automatically a legend
for all objects in the pad. It is called `TPad::BuildLegend()`.

Each item in the legend is added using the `AddEntry` method. This
method defines the object to be added (by reference or name), the label
associated to this object and an option which a combination of:

  - L: draw line associated with TAttLine if obj inherits from TAttLine
  - P: draw polymarker associated with TAttMarker if obj inherits from TAttMarker
  - F: draw a box with fill associated wit TAttFill if obj inherits TAttFill
  - E: draw vertical error bar

As shown in the following example, passing a NULL pointer as first parameter in
`AddEntry` is also valid. This allows to add text or blank lines in a
legend.

Begin_Macro(source)
{
   auto c2 = new TCanvas("c2","c2",500,300);

   auto* legend = new TLegend(0.2, 0.2, .8, .8);
   auto h = new TH1F("", "", 1, 0, 1);

   legend->AddEntry(h, "Histogram \"h\"", "l");
   legend->AddEntry((TObject*)0, "", "");
   legend->AddEntry((TObject*)0, "Some text", "");
   legend->AddEntry((TObject*)0, "", "");
   legend->AddEntry(h, "Histogram \"h\" again", "l");

   legend->Draw();
}
End_Macro

It is possible to draw the legend entries over several columns using
the method `SetNColumns()` like in the following example.

Begin_Macro(source)
{
   auto c3 = new TCanvas("c2","c2",500,300);

   auto legend = new TLegend(0.2, 0.2, .8, .8);
   auto h = new TH1F("", "", 1, 0, 1);

   legend->SetNColumns(2);

   legend->AddEntry(h, "Column 1 line 1", "l");
   legend->AddEntry(h, "Column 2 line 1", "l");
   legend->AddEntry(h, "Column 1 line 2", "l");
   legend->AddEntry(h, "Column 2 line 2", "l");

   legend->Draw();
}
End_Macro

\since **ROOT version 6.09/03**

The legend can be placed automatically in the current pad in an empty space
found at painting time.

The following example illustrate this facility. Only the width and height of the
legend is specified in percentage of the pad size.

Begin_Macro(source)
../../../tutorials/hist/legendautoplaced.C
End_Macro

*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.
/// This constructor allows to place automatically the legend with a default
/// width(0.3) and a default height (0.15) in normalize coordinates.

TLegend::TLegend(): TPave(0.3,0.15,0.3,0.15,4,"brNDC"),
                    TAttText(12,0,1,gStyle->GetLegendFont(),0)
{
   fPrimitives = 0;
   SetDefaults();
   SetBorderSize(gStyle->GetLegendBorderSize());
   SetFillColor(gStyle->GetLegendFillColor());
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.
///
/// A TLegend is a Pave with several TLegendEntry(s).
///
/// x1,y1,x2,y2 are the coordinates of the Legend in the current pad
/// (in normalised coordinates by default)
///
/// `header` is the title displayed at the top of the legend
/// it is a TLatex string treated like a regular entry. The default
/// is no header (header = 0).
///
/// The options are the same as for TPave.

TLegend::TLegend( Double_t x1, Double_t y1,Double_t x2, Double_t y2,
                  const char *header, Option_t *option)
        :TPave(x1,y1,x2,y2,4,option), TAttText(12,0,1,gStyle->GetLegendFont(),0)
{
   fPrimitives = new TList;
   if ( header && strlen(header) > 0) {
      TLegendEntry *headerEntry = new TLegendEntry( 0, header, "h" );
      headerEntry->SetTextAlign(0);
      headerEntry->SetTextAngle(0);
      headerEntry->SetTextColor(0);
      headerEntry->SetTextFont(gStyle->GetLegendFont());
      headerEntry->SetTextSize(0);
      fPrimitives->AddFirst(headerEntry);
   }
   SetDefaults();
   SetBorderSize(gStyle->GetLegendBorderSize());
   SetFillColor(gStyle->GetLegendFillColor());
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with automatic placement.
///
/// A TLegend is a Pave with several TLegendEntry(s).
///
/// This constructor doesn't define the legend position. `w` and `h` are the
/// width and height of the legend in percentage of the current pad size.
/// The position will be automatically defined at painting time.
///
/// `header` is the title displayed at the top of the legend
/// it is a TLatex string treated like a regular entry. The default
/// is no header (header = 0).
///
/// The options are the same as for TPave.

TLegend::TLegend( Double_t w, Double_t h, const char *header, Option_t *option)
        :TPave(w,h,w,h,4,option), TAttText(12,0,1,gStyle->GetLegendFont(),0)
{
   fPrimitives = new TList;
   if ( header && strlen(header) > 0) {
      TLegendEntry *headerEntry = new TLegendEntry( 0, header, "h" );
      headerEntry->SetTextAlign(0);
      headerEntry->SetTextAngle(0);
      headerEntry->SetTextColor(0);
      headerEntry->SetTextFont(gStyle->GetLegendFont());
      headerEntry->SetTextSize(0);
      fPrimitives->AddFirst(headerEntry);
   }
   SetDefaults();
   SetBorderSize(gStyle->GetLegendBorderSize());
   SetFillColor(gStyle->GetLegendFillColor());
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TLegend::TLegend( const TLegend &legend ) : TPave(legend), TAttText(legend),
                                            fPrimitives(0)
{
  if (legend.fPrimitives) {
      fPrimitives = new TList();
      TListIter it(legend.fPrimitives);
      while (TLegendEntry *e = (TLegendEntry *)it.Next()) {
         TLegendEntry *newentry = new TLegendEntry(*e);
         fPrimitives->Add(newentry);
      }
   }
   ((TLegend&)legend).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TLegend& TLegend::operator=(const TLegend &lg)
{
   if(this!=&lg) {
      TPave::operator=(lg);
      TAttText::operator=(lg);
      fPrimitives=lg.fPrimitives;
      fEntrySeparation=lg.fEntrySeparation;
      fMargin=lg.fMargin;
      fNColumns=lg.fNColumns;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Default destructor.

TLegend::~TLegend()
{
   if (fPrimitives) fPrimitives->Delete();
   delete fPrimitives;
   fPrimitives = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new entry to this legend. "obj" is the object to be represented.
/// "label" is the text you wish to associate with obj in the legend.
/// If "label" is null or empty, the title of the object will be used.
///
/// Options are:
///
///  - L: draw line associated with TAttLine if obj inherits from TAttLine
///  - P: draw polymarker associated with TAttMarker if obj inherits from TAttMarker
///  - F: draw a box with fill associated wit TAttFill if obj inherits TAttFill
///  - E: draw vertical error bar if option "L" is also specified

TLegendEntry *TLegend::AddEntry(const TObject *obj, const char *label, Option_t *option)
{
   const char *lab = label;

   if (obj && (!label || strlen(label)==0)) lab = obj->GetTitle();
   TLegendEntry *newentry = new TLegendEntry( obj, lab, option );
   if ( !fPrimitives ) fPrimitives = new TList;
   fPrimitives->Add(newentry);
   return newentry;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new entry to this legend. "name" is the name of an object in the pad to
/// be represented label is the text you wish to associate with obj in the legend
/// if label is null or empty, the title of the object will be used.
///
/// Options are:
///
///  - L: draw line associated with TAttLine if obj inherits from TAttLine
///  - P: draw polymarker associated with TAttMarker if obj inherits from TAttMarker
///  - F: draw a box with fill associated wit TAttFill if obj inherits TAttFill
///  - E: draw vertical error bar if option "L" is also specified

TLegendEntry *TLegend::AddEntry(const char *name, const char *label, Option_t *option)
{
   if (!gPad) {
      Error("AddEntry", "need to create a canvas first");
      return 0;
   }

   TObject *obj = gPad->FindObject(name);

   // If the object "name" has not been found, the following code tries to
   // find it in TMultiGraph or THStack possibly present in the current pad.
   if (!obj) {
      TList *lop = gPad->GetListOfPrimitives();
      if (lop) {
         TObject *o=0;
         TIter next(lop);
         while( (o=next()) ) {
            if ( o->InheritsFrom(TMultiGraph::Class() ) ) {
               TList * grlist = ((TMultiGraph *)o)->GetListOfGraphs();
               obj = grlist->FindObject(name);
               if (obj) break;
            }
            if ( o->InheritsFrom(THStack::Class() ) ) {
               TList * hlist = ((THStack *)o)->GetHists();
               obj = hlist->FindObject(name);
               if (obj) break;
            }
         }
      }
   }

   return AddEntry( obj, label, option );
}

////////////////////////////////////////////////////////////////////////////////
/// Clear all entries in this legend, including the header.

void TLegend::Clear( Option_t *)
{
   if (!fPrimitives) return;
   fPrimitives->Delete();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this legend into "obj".

void TLegend::Copy( TObject &obj ) const
{
   TPave::Copy(obj);
   TAttText::Copy((TLegend&)obj);
   ((TLegend&)obj).fEntrySeparation = fEntrySeparation;
   ((TLegend&)obj).fMargin = fMargin;
   ((TLegend&)obj).fNColumns = fNColumns;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete entry at the mouse position.

void TLegend::DeleteEntry()
{
   if ( !fPrimitives ) return;
   TLegendEntry* entry = GetEntry();   // get entry pointed by the mouse
   if ( !entry ) return;
   fPrimitives->Remove(entry);
   delete entry;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this legend with its current attributes.

void TLegend::Draw( Option_t *option )
{
   AppendPad(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Edit the fill attributes for the entry pointed by the mouse.

void TLegend::EditEntryAttFill()
{
   TLegendEntry* entry = GetEntry();   // get entry pointed by the mouse
   if ( !entry ) return;
   gROOT->SetSelectedPrimitive( entry );
   entry->SetFillAttributes();
}

////////////////////////////////////////////////////////////////////////////////
/// Edit the line attributes for the entry pointed by the mouse.

void TLegend::EditEntryAttLine()
{
   TLegendEntry* entry = GetEntry();   // get entry pointed by the mouse
   if ( !entry ) return;
   gROOT->SetSelectedPrimitive( entry );
   entry->SetLineAttributes();
}

////////////////////////////////////////////////////////////////////////////////
/// Edit the marker attributes for the entry pointed by the mouse.

void TLegend::EditEntryAttMarker()
{
   TLegendEntry* entry = GetEntry();   // get entry pointed by the mouse
   if ( !entry ) return;
   gROOT->SetSelectedPrimitive( entry );
   entry->SetMarkerAttributes();
}

////////////////////////////////////////////////////////////////////////////////
/// Edit the text attributes for the entry pointed by the mouse.

void TLegend::EditEntryAttText()
{
   TLegendEntry* entry = GetEntry();   // get entry pointed by the mouse
   if ( !entry ) return;
   gROOT->SetSelectedPrimitive( entry );
   entry->SetTextAttributes();
}

////////////////////////////////////////////////////////////////////////////////
/// Get entry pointed to by the mouse.
/// This method is mostly a tool for other methods inside this class.

TLegendEntry *TLegend::GetEntry() const
{
   if (!gPad) {
      Error("GetEntry", "need to create a canvas first");
      return 0;
   }

   Int_t nRows = GetNRows();
   if ( nRows == 0 ) return 0;

   Double_t ymouse = gPad->AbsPixeltoY(gPad->GetEventY())-fY1;
   Double_t yspace = (fY2 - fY1)/nRows;

   Int_t nColumns = GetNColumns();
   Double_t xmouse = gPad->AbsPixeltoX(gPad->GetEventX())-fX1;
   Double_t xspace = 0.;
   if (nColumns > 0) xspace = (fX2 - fX1)/nColumns;

   Int_t ix = 1;
   if (xspace > 0.) ix = (Int_t)(xmouse/xspace)+1;
   if (ix > nColumns) ix = nColumns;
   if (ix < 1)        ix = 1;

   Int_t iy = nRows-(Int_t)(ymouse/yspace);
   if (iy > nRows) iy = nRows;
   if (iy < 1)     iy = 1;

   Int_t nloops = TMath::Min(ix+(nColumns*(iy-1)), fPrimitives->GetSize());

   TIter next(fPrimitives);
   TLegendEntry *entry = 0;

   for (Int_t i=1; i<= nloops; i++) entry = (TLegendEntry *)next();

   return entry;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the header, which is the title that appears at the top
/// of the legend.

const char *TLegend::GetHeader() const
{
   if ( !fPrimitives ) return 0;
      TIter next(fPrimitives);
   TLegendEntry *first;   // header is always the first entry
   if ((  first = (TLegendEntry*)next()  )) {
      TString opt = first->GetOption();
      opt.ToLower();
      if ( opt.Contains("h") ) return first->GetLabel();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new entry before the entry at the mouse position.

void TLegend::InsertEntry( const char* objectName, const char* label, Option_t* option)
{
   if (!gPad) {
      Error("InsertEntry", "need to create a canvas first");
      return;
   }

   TLegendEntry* beforeEntry = GetEntry();   // get entry pointed by the mouse
   TObject *obj = gPad->FindObject( objectName );

   // note either obj OR beforeEntry may be zero at this point

   TLegendEntry *newentry = new TLegendEntry( obj, label, option );

   if ( !fPrimitives ) fPrimitives = new TList;
   if ( beforeEntry ) {
      fPrimitives->AddBefore( (TObject*)beforeEntry, (TObject*)newentry );
   } else {
      fPrimitives->Add((TObject*)newentry);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this legend with its current attributes.

void TLegend::Paint( Option_t* option )
{
   // The legend need to be placed automatically in some empty space
   if (fX1 == fX2 && fY1 == fY2) {
      if (gPad->PlaceBox(this, fX1, fY1, fX1, fY1)) {
         fY2 = fY2+fY1;
         fX2 = fX2+fX1;
      } else {
         Warning("Paint", "Legend too large to be automatically placed; a default position is used");
         fX1 = 0.5;
         fY1 = 0.67;
         fX2 = 0.88;
         fY2 = 0.88;
      }
   }

   // Paint the Legend
   TPave::ConvertNDCtoPad();
   TPave::PaintPave(fX1,fY1,fX2,fY2,GetBorderSize(),option);
   PaintPrimitives();
}

////////////////////////////////////////////////////////////////////////////////
/// Get the number of rows.

Int_t TLegend::GetNRows() const
{
   Int_t nEntries = 0;
   if ( fPrimitives ) nEntries = fPrimitives->GetSize();
   if ( nEntries == 0 ) return 0;

   Int_t nRows;
   if(GetHeader() != NULL) nRows = 1 + (Int_t) TMath::Ceil((Double_t) (nEntries-1)/fNColumns);
   else  nRows = (Int_t) TMath::Ceil((Double_t) nEntries/fNColumns);

   return nRows;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the number of columns for the legend. The header, if set, is given
/// its own row. After that, every nColumns entries are inserted into the
/// same row. For example, if one calls legend.SetNColumns(2), and there
/// is no header, then the first two TObjects added to the legend will be
/// in the first row, the next two will appear in the second row, and so on.

void TLegend::SetNColumns(Int_t nColumns)
{
   if(nColumns < 1) {
      Warning("TLegend::SetNColumns", "illegal value nColumns = %d; keeping fNColumns = %d", nColumns, fNColumns);
      return;
   }
   fNColumns = nColumns;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the entries (list of primitives) for this legend.

void TLegend::PaintPrimitives()
{
   Int_t nRows = GetNRows();
   if ( nRows == 0 ) return;

   // Evaluate text size as a function of the number of entries
   //  taking into account their real size after drawing latex
   // Note: in pixel coords y1 > y2=0, but x2 > x1=0
   //       in NDC          y2 > y1,   and x2 > x1

   Double_t x1 = fX1NDC;
   Double_t y1 = fY1NDC;
   Double_t x2 = fX2NDC;
   Double_t y2 = fY2NDC;
   Double_t margin = fMargin*( x2-x1 )/fNColumns;
   Double_t boxwidth = margin;
   Double_t boxw = boxwidth*0.35;
   Double_t yspace = (y2-y1)/nRows;
   Double_t yspace2 = yspace/2.;
   Double_t textsize = GetTextSize();
   Double_t save_textsize = textsize;
   if (textsize==0.) {
      SetTextSize(gStyle->GetLegendTextSize());
      textsize = GetTextSize();
   }
   Bool_t autosize = kFALSE;
   Double_t* columnWidths = new Double_t[fNColumns];
   memset(columnWidths, 0, fNColumns*sizeof(Double_t));

   if ( textsize == 0 ) {
      autosize = kTRUE;
      textsize = ( 1. - fEntrySeparation ) * yspace;

      // find the max width and height (in pad coords) of one latex entry label
      Double_t maxentrywidth = 0, maxentryheight = 0;
      TIter nextsize(fPrimitives);
      TLegendEntry *entrysize;
      Int_t iColumn = 0;
      while (( entrysize = (TLegendEntry *)nextsize() )) {
         TLatex entrytex( 0, 0, entrysize->GetLabel() );
         entrytex.SetNDC();
         Style_t tfont = entrysize->GetTextFont();
         if (tfont == 0) tfont = GetTextFont();
         if (tfont%10 == 3) --tfont;
         entrytex.SetTextFont(tfont);
         entrytex.SetTextSize(textsize);
         if ( entrytex.GetYsize() > maxentryheight ) {
            maxentryheight = entrytex.GetYsize();
         }
         TString opt = entrysize->GetOption();
         opt.ToLower();
         if ( opt.Contains("h") ) {
            if ( entrytex.GetXsize() > maxentrywidth ) {
               maxentrywidth = entrytex.GetXsize();
            }
         } else {
            if ( entrytex.GetXsize() > columnWidths[iColumn] ) {
               columnWidths[iColumn] = entrytex.GetXsize();
            }
            iColumn++;
            iColumn %= fNColumns;
         }
         Double_t tmpMaxWidth = 0.0;
         for(int i=0; i<fNColumns; i++) tmpMaxWidth += columnWidths[i];
         if ( tmpMaxWidth > maxentrywidth) maxentrywidth = tmpMaxWidth;
      }
      // make sure all labels fit in the allotted space
      Double_t tmpsize_h = maxentryheight /(gPad->GetY2() - gPad->GetY1());
      textsize = TMath::Min( textsize, tmpsize_h );
      Double_t tmpsize_w = textsize*(fX2-fX1)*(1.0-fMargin)/maxentrywidth;
      if(fNColumns > 1) tmpsize_w = textsize*(fX2-fX1)*(1.0-fMargin-fColumnSeparation)/maxentrywidth;
      textsize = TMath::Min( textsize, tmpsize_w );
      SetTextSize( textsize );
   }

   // Update column widths, put into NDC units
   // block off this section of code to make sure all variables are local:
   // don't want to ruin initialisation of these variables later on
   {
      TIter next(fPrimitives);
      TLegendEntry *entry;
      Int_t iColumn = 0;
      memset(columnWidths, 0, fNColumns*sizeof(Double_t));
      while (( entry = (TLegendEntry *)next() )) {
         TLatex entrytex( 0, 0, entry->GetLabel() );
         entrytex.SetNDC();
         Style_t tfont = entry->GetTextFont();
         if (tfont == 0) tfont = GetTextFont();
         if (autosize && tfont%10 == 3) --tfont;
         entrytex.SetTextFont(tfont);
         if(entry->GetTextSize() == 0) entrytex.SetTextSize(textsize);
         TString opt = entry->GetOption();
         opt.ToLower();
         if (!opt.Contains("h")) {
            if ( entrytex.GetXsize() > columnWidths[iColumn] ) {
               columnWidths[iColumn] = entrytex.GetXsize();
            }
            iColumn++;
            iColumn %= fNColumns;
         }
      }
      double totalWidth = 0.0;
      for(int i=0; i<fNColumns; i++) totalWidth += columnWidths[i];
      if(fNColumns > 1) totalWidth /= (1.0-fMargin-fColumnSeparation);
      else totalWidth /= (1.0 - fMargin);
      for(int i=0; i<fNColumns; i++) {
         columnWidths[i] = columnWidths[i]/totalWidth*(x2-x1) + margin;
      }
   }

   Double_t ytext = y2 + 0.5*yspace;  // y-location of 0th entry

   // iterate over and paint all the TLegendEntries
   TIter next(fPrimitives);
   TLegendEntry *entry;
   Int_t iColumn = 0;
   while (( entry = (TLegendEntry *)next() )) {
      if(iColumn == 0) ytext -= yspace;

      // Draw Label in Latexmargin

      Short_t talign = entry->GetTextAlign();
      Float_t tangle = entry->GetTextAngle();
      Color_t tcolor = entry->GetTextColor();
      Style_t tfont  = entry->GetTextFont();
      Size_t  tsize  = entry->GetTextSize();
      // if the user hasn't set a parameter, then set it to the TLegend value
      if (talign == 0) entry->SetTextAlign(GetTextAlign());
      if (tangle == 0) entry->SetTextAngle(GetTextAngle());
      if (tcolor == 0) entry->SetTextColor(GetTextColor());
      if (tfont  == 0) {
         tfont = GetTextFont();
         if (autosize && tfont%10 == 3) --tfont;
         entry->SetTextFont(tfont);
      }
      if (tsize  == 0) entry->SetTextSize(GetTextSize());
      // set x,y according to the requested alignment
      Double_t x=0,y=0;
      Int_t halign = entry->GetTextAlign()/10;
      Double_t entrymargin = margin;
      // for the header the margin is near zero
      TString opt = entry->GetOption();
      opt.ToLower();
      x1 = fX1NDC;
      x2 = fX2NDC;
      if ( opt.Contains("h") ) entrymargin = margin/10.;
      else if (fNColumns > 1) {
         for(int i=0; i<iColumn; i++) x1 += columnWidths[i] + fColumnSeparation*(fX2NDC-fX1NDC)/(fNColumns-1);
         x2 = x1 + columnWidths[iColumn];
         iColumn++;
         iColumn %= fNColumns;
      }
      if (halign == 1) x = x1 + entrymargin;
      if (halign == 2) x = 0.5*( (x1+entrymargin) + x2 );
      if (halign == 3) x = x2 - entrymargin/10.;
      Int_t valign = entry->GetTextAlign()%10;

      if (valign == 1) y = ytext - (1. - fEntrySeparation)* yspace2;
      if (valign == 3) y = ytext + (1. - fEntrySeparation)* yspace2;

      // The vertical alignment "centered" is treated in a special way
      // to ensure a better spacing between lines.
      if (valign == 2) {
         Float_t tsizepad = textsize;
         if (tfont%10 == 3) tsizepad = (gPad->AbsPixeltoY(0) - gPad->AbsPixeltoY(textsize))/(gPad->GetY2() - gPad->GetY1());
         if (yspace2 < tsizepad) {
            entry->SetTextAlign(10*halign+1);
            y = ytext - (1. - fEntrySeparation)* yspace2/2.;
         } else {
            y = ytext;
         }
      }

      TLatex entrytex( x, y, entry->GetLabel() );
      entrytex.SetNDC();
      entry->TAttText::Copy(entrytex);
      entrytex.Paint();

      // reset attributes back to their original values
      entry->SetTextAlign(talign);
      entry->SetTextAngle(tangle);
      entry->SetTextColor(tcolor);
      entry->SetTextFont(tfont);
      entry->SetTextSize(tsize);

      // define x,y as the center of the symbol for this entry
      Double_t xsym = x1 + margin/2.;
      Double_t ysym = ytext;

      TObject *eobj = entry->GetObject();

      // depending on the object drawing option, the endcaps for error
      // bar are drawn differently.
      Int_t endcaps  = 0; // no endcaps.
      if (eobj) { // eobj == nullptr for the legend header
         TString eobjopt = eobj->GetDrawOption();
         eobjopt.ToLower();
         if (eobjopt.Contains("e1") && eobj->InheritsFrom(TH1::Class())) endcaps = 1; // a bar
         if (eobj->InheritsFrom(TGraph::Class())) {
            endcaps = 1; // a bar, default for TGraph
            if (eobjopt.Contains("z"))  endcaps = 0; // no endcaps.
            if (eobjopt.Contains(">"))  endcaps = 2; // empty arrow.
            if (eobjopt.Contains("|>")) endcaps = 3; // filled arrow.
         }
      }

      // Draw fill pattern (in a box)

      if ( opt.Contains("f")) {
         if (eobj && eobj->InheritsFrom(TAttFill::Class())) {
            dynamic_cast<TAttFill*>(eobj)->Copy(*entry);
         }

         // box total height is yspace*0.7
         entry->TAttFill::Modify();
         Double_t xf[4],yf[4];
         xf[0] = xsym - boxw;
         yf[0] = ysym - yspace*0.35;
         xf[1] = xsym + boxw;
         yf[1] = yf[0];
         xf[2] = xf[1];
         yf[2] = ysym + yspace*0.35;
         xf[3] = xf[0];
         yf[3] = yf[2];
         for (Int_t i=0;i<4;i++) {
            xf[i] = gPad->GetX1() + xf[i]*(gPad->GetX2()-gPad->GetX1());
            yf[i] = gPad->GetY1() + yf[i]*(gPad->GetY2()-gPad->GetY1());
         }
         gPad->PaintFillArea(4,xf,yf);
      }

      // Get Polymarker size

      Double_t symbolsize = 0.;
      TMarker entrymarker( xsym, ysym, 0 );

      if ( opt.Contains("p")) {
         if (eobj && eobj->InheritsFrom(TAttMarker::Class())) {
            dynamic_cast<TAttMarker*>(eobj)->Copy(*entry);
         }
         entrymarker.SetNDC();
         entry->TAttMarker::Copy(entrymarker);
         if (entrymarker.GetMarkerStyle() >= 5 ) symbolsize = entrymarker.GetMarkerSize();
      }

      // Draw line

      if ( opt.Contains("l") || opt.Contains("f")) {

         if (eobj && eobj->InheritsFrom(TAttLine::Class())) {
            dynamic_cast<TAttLine*>(eobj)->Copy(*entry);
         }
         // line total length (in x) is margin*0.8
         TLine entryline( xsym - boxw, ysym, xsym + boxw, ysym );
         entryline.SetBit(TLine::kLineNDC);
         entry->TAttLine::Copy(entryline);
         // if the entry is filled, then surround the box with the line instead
         if ( opt.Contains("f") && !opt.Contains("l")) {
            // box total height is yspace*0.7
            boxwidth = yspace*
               (gPad->GetX2()-gPad->GetX1())/(gPad->GetY2()-gPad->GetY1());
            if ( boxwidth > margin ) boxwidth = margin;

            entryline.PaintLineNDC( xsym - boxw, ysym + yspace*0.35,
                                    xsym + boxw, ysym + yspace*0.35);
            entryline.PaintLineNDC( xsym - boxw, ysym - yspace*0.35,
                                    xsym + boxw, ysym - yspace*0.35);
            entryline.PaintLineNDC( xsym + boxw, ysym - yspace*0.35,
                                    xsym + boxw, ysym + yspace*0.35);
            entryline.PaintLineNDC( xsym - boxw, ysym - yspace*0.35,
                                    xsym - boxw, ysym + yspace*0.35);
         } else {
            entryline.Paint();
            if (opt.Contains("e")) {
               if ( !opt.Contains("p")) {
                  entryline.PaintLineNDC( xsym, ysym - yspace*0.30,
                                          xsym, ysym + yspace*0.30);
               } else {
                  Double_t sy  = (fY2NDC-fY1NDC)*((0.5*(gPad->PixeltoY(0) - gPad->PixeltoY(Int_t(symbolsize*8.))))/(fY2-fY1));
                  TLine entryline1(xsym, ysym + sy, xsym, ysym + yspace*0.30);
                  entryline1.SetBit(TLine::kLineNDC);
                  entry->TAttLine::Copy(entryline1);
                  entryline1.Paint();
                  TLine entryline2(xsym, ysym - sy, xsym, ysym - yspace*0.30);
                  entryline2.SetBit(TLine::kLineNDC);
                  entry->TAttLine::Copy(entryline2);
                  entryline2.Paint();
               }
               Double_t barw = boxw*0.1*gStyle->GetEndErrorSize();
               if (endcaps == 1) {
                  TLine entrytop1(xsym-barw, ysym + yspace*0.30, xsym+barw, ysym + yspace*0.30);
                  entrytop1.SetBit(TLine::kLineNDC);
                  entry->TAttLine::Copy(entrytop1);
                  entrytop1.Paint();
                  TLine entrytop2(xsym-barw, ysym - yspace*0.30, xsym+barw, ysym - yspace*0.30);
                  entrytop2.SetBit(TLine::kLineNDC);
                  entry->TAttLine::Copy(entrytop2);
                  entrytop2.Paint();
               } else if (endcaps == 2) {
                  Double_t xe1[3] = {xsym-barw, xsym ,xsym+barw};
                  Double_t ye1[3] = {ysym+yspace*0.20, ysym + yspace*0.30 ,ysym+yspace*0.20};
                  TPolyLine ple1(3,xe1,ye1);
                  ple1.SetBit(TLine::kLineNDC);
                  entry->TAttLine::Copy(ple1);
                  ple1.Paint();
                  Double_t xe2[3] = {xsym-barw, xsym ,xsym+barw};
                  Double_t ye2[3] = {ysym-yspace*0.20, ysym - yspace*0.30 ,ysym-yspace*0.20};
                  TPolyLine ple2(3,xe2,ye2);
                  ple2.SetBit(TLine::kLineNDC);
                  entry->TAttLine::Copy(ple2);
               } else if (endcaps == 3) {
                  Double_t xe1[3] = {xsym-barw, xsym ,xsym+barw};
                  Double_t ye1[3] = {ysym+yspace*0.20, ysym + yspace*0.30 ,ysym+yspace*0.20};
                  Double_t xe2[3] = {xsym-barw, xsym ,xsym+barw};
                  Double_t ye2[3] = {ysym-yspace*0.20, ysym - yspace*0.30 ,ysym-yspace*0.20};
                  for (Int_t i=0;i<3;i++) {
                     xe1[i] = gPad->GetX1() + xe1[i]*(gPad->GetX2()-gPad->GetX1());
                     ye1[i] = gPad->GetY1() + ye1[i]*(gPad->GetY2()-gPad->GetY1());
                     xe2[i] = gPad->GetX1() + xe2[i]*(gPad->GetX2()-gPad->GetX1());
                     ye2[i] = gPad->GetY1() + ye2[i]*(gPad->GetY2()-gPad->GetY1());
                  }
                  TPolyLine ple1(3,xe1,ye1);
                  ple1.SetFillColor(entry->GetLineColor());
                  ple1.SetFillStyle(1001);
                  ple1.Paint("f");
                  TPolyLine ple2(3,xe2,ye2);
                  ple2.SetFillColor(entry->GetLineColor());
                  ple2.SetFillStyle(1001);
                  ple2.Paint("f");
               }
            }
         }
      }

      // Draw error only

      if (opt.Contains("e") && !(opt.Contains("l") || opt.Contains("f"))) {
         if (eobj && eobj->InheritsFrom(TAttLine::Class())) {
            dynamic_cast<TAttLine*>(eobj)->Copy(*entry);
         }
         if ( !opt.Contains("p")) {
            TLine entryline(xsym, ysym - yspace*0.30,
                            xsym, ysym + yspace*0.30);
            entryline.SetBit(TLine::kLineNDC);
            entry->TAttLine::Copy(entryline);
            entryline.Paint();
         } else {
            Double_t sy  = (fY2NDC-fY1NDC)*((0.5*(gPad->PixeltoY(0) - gPad->PixeltoY(Int_t(symbolsize*8.))))/(fY2-fY1));
            TLine entryline1(xsym, ysym + sy, xsym, ysym + yspace*0.30);
            entryline1.SetBit(TLine::kLineNDC);
            entry->TAttLine::Copy(entryline1);
            entryline1.Paint();
            TLine entryline2(xsym, ysym - sy, xsym, ysym - yspace*0.30);
            entryline2.SetBit(TLine::kLineNDC);
            entry->TAttLine::Copy(entryline2);
            entryline2.Paint();
         }
         Double_t barw = boxw*0.1*gStyle->GetEndErrorSize();
         if (endcaps == 1) {
            TLine entrytop1(xsym-barw, ysym + yspace*0.30, xsym+barw, ysym + yspace*0.30);
            entrytop1.SetBit(TLine::kLineNDC);
            entry->TAttLine::Copy(entrytop1);
            entrytop1.Paint();
            TLine entrytop2(xsym-barw, ysym - yspace*0.30, xsym+barw, ysym - yspace*0.30);
            entrytop2.SetBit(TLine::kLineNDC);
            entry->TAttLine::Copy(entrytop2);
            entrytop2.Paint();
         } else if (endcaps == 2) {
            Double_t xe1[3] = {xsym-barw, xsym ,xsym+barw};
            Double_t ye1[3] = {ysym+yspace*0.20, ysym + yspace*0.30 ,ysym+yspace*0.20};
            TPolyLine ple1(3,xe1,ye1);
            ple1.SetBit(TLine::kLineNDC);
            entry->TAttLine::Copy(ple1);
            ple1.Paint();
            Double_t xe2[3] = {xsym-barw, xsym ,xsym+barw};
            Double_t ye2[3] = {ysym-yspace*0.20, ysym - yspace*0.30 ,ysym-yspace*0.20};
            TPolyLine ple2(3,xe2,ye2);
            ple2.SetBit(TLine::kLineNDC);
            entry->TAttLine::Copy(ple2);
            ple2.Paint();
         } else if (endcaps == 3) {
            Double_t xe1[3] = {xsym-barw, xsym ,xsym+barw};
            Double_t ye1[3] = {ysym+yspace*0.20, ysym + yspace*0.30 ,ysym+yspace*0.20};
            Double_t xe2[3] = {xsym-barw, xsym ,xsym+barw};
            Double_t ye2[3] = {ysym-yspace*0.20, ysym - yspace*0.30 ,ysym-yspace*0.20};
            for (Int_t i=0;i<3;i++) {
               xe1[i] = gPad->GetX1() + xe1[i]*(gPad->GetX2()-gPad->GetX1());
               ye1[i] = gPad->GetY1() + ye1[i]*(gPad->GetY2()-gPad->GetY1());
               xe2[i] = gPad->GetX1() + xe2[i]*(gPad->GetX2()-gPad->GetX1());
               ye2[i] = gPad->GetY1() + ye2[i]*(gPad->GetY2()-gPad->GetY1());
            }
            TPolyLine ple1(3,xe1,ye1);
            ple1.SetFillColor(entry->GetLineColor());
            ple1.SetFillStyle(1001);
            ple1.Paint("f");
            TPolyLine ple2(3,xe2,ye2);
            ple2.SetFillColor(entry->GetLineColor());
            ple2.SetFillStyle(1001);
            ple2.Paint("f");
         }
      }

      // Draw Polymarker
      if ( opt.Contains("p"))  entrymarker.Paint();
   }
   SetTextSize(save_textsize);
   delete [] columnWidths;
}

////////////////////////////////////////////////////////////////////////////////
/// Dump this TLegend and its contents.

void TLegend::Print( Option_t* option ) const
{
   TPave::Print( option );
   if (fPrimitives) fPrimitives->Print();
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the legend entries pointing to "obj".

void TLegend::RecursiveRemove(TObject *obj)
{
   TIter next(fPrimitives);
   TLegendEntry *entry;
   while (( entry = (TLegendEntry *)next() )) {
      if (entry->GetObject() == obj) entry->SetObject((TObject*)0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save this legend as C++ statements on output stream out
/// to be used with the SaveAs .C option.

void TLegend::SavePrimitive(std::ostream &out, Option_t* )
{

   out << "   " << std::endl;
   char quote = '"';
   if ( gROOT->ClassSaved( TLegend::Class() ) ) {
      out << "   ";
   } else {
      out << "   TLegend *";
   }
   // note, we can always use NULL header, since its included in primitives
   out << "leg = new TLegend("<<GetX1NDC()<<","<<GetY1NDC()<<","
       <<GetX2NDC()<<","<<GetY2NDC()<<","
       << "NULL" << "," <<quote<< fOption <<quote<<");" << std::endl;
   if (fBorderSize != 4) {
      out<<"   leg->SetBorderSize("<<fBorderSize<<");"<<std::endl;
   }
   SaveTextAttributes(out,"leg",12,0,1,42,0);
   SaveLineAttributes(out,"leg",-1,-1,-1);
   SaveFillAttributes(out,"leg",-1,-1);
   if ( fPrimitives ) {
      TIter next(fPrimitives);
      TLegendEntry *entry;
      while (( entry = (TLegendEntry *)next() )) entry->SaveEntry(out,"leg");
   }
   out << "   leg->Draw();"<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Edit the label of the entry pointed to by the mouse.

void TLegend::SetEntryLabel( const char* label )
{
   TLegendEntry* entry = GetEntry();   // get entry pointed by the mouse
   if ( entry ) entry->SetLabel( label );
}

////////////////////////////////////////////////////////////////////////////////
/// Edit the option of the entry pointed to by the mouse.

void TLegend::SetEntryOption( Option_t* option )
{
   TLegendEntry* entry = GetEntry();   // get entry pointed by the mouse
   if ( entry ) entry->SetOption( option );
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the header, which is the "title" that appears at the top of the legend.
/// If `option` contains `C`, the title will be centered.

void TLegend::SetHeader( const char *header, Option_t* option )
{
   TString opt;

   if ( !fPrimitives ) fPrimitives = new TList;
   TIter next(fPrimitives);
   TLegendEntry *first;   // header is always the first entry
   if ((  first = (TLegendEntry*)next() )) {
      opt = first->GetOption();
      opt.ToLower();
      if ( opt.Contains("h") ) {
         first->SetLabel(header);
         opt = option;
         opt.ToLower();
         if ( opt.Contains("c") ) first->SetTextAlign(22);
         else                     first->SetTextAlign(0);
         return;
      }
   }
   first = new TLegendEntry( 0, header, "h" );
   opt = option;
   opt.ToLower();
   if ( opt.Contains("c") ) first->SetTextAlign(22);
   else                     first->SetTextAlign(0);
   first->SetTextAngle(0);
   first->SetTextColor(0);
   first->SetTextFont(GetTextFont()); // default font is TLegend font for the header
   first->SetTextSize(0);
   fPrimitives->AddFirst((TObject*)first);
}
