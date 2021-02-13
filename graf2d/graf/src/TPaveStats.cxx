// @(#)root/graf:$Id$
// Author: Rene Brun   15/03/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "TROOT.h"
#include "TBuffer.h"
#include "TPaveStats.h"
#include "TPaveLabel.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "TLatex.h"
#include "strlcpy.h"

ClassImp(TPaveStats);

/** \class TPaveStats
\ingroup BasicGraphics

The histogram statistics painter class.

To draw histogram statistics and fit parameters.

- [Statistics Display](#PS01)
- [Fit Statistics](#PS02)
- [Statistics box editing](#PS03)

## <a name="PS01"></a> Statistics Display

The type of information shown in the histogram statistics box can be selected
with:
~~~ {.cpp}
      gStyle->SetOptStat(mode);
~~~

The "`mode`" has up to nine digits that can be set to on (1 or 2), off (0).
~~~ {.cpp}
      mode = ksiourmen  (default = 000001111)
      k = 1;  kurtosis printed
      k = 2;  kurtosis and kurtosis error printed
      s = 1;  skewness printed
      s = 2;  skewness and skewness error printed
      i = 1;  integral of bins printed
      o = 1;  number of overflows printed
      u = 1;  number of underflows printed
      r = 1;  rms printed
      r = 2;  rms and rms error printed
      m = 1;  mean value printed
      m = 2;  mean and mean error values printed
      e = 1;  number of entries printed
      n = 1;  name of histogram is printed
~~~

For example:
~~~ {.cpp}
      gStyle->SetOptStat(11);
~~~
displays only the name of histogram and the number of entries, whereas:
~~~ {.cpp}
      gStyle->SetOptStat(1101);
~~~
displays the name of histogram, mean value and RMS.

<b>WARNING 1:</b> never do:
~~~ {.cpp}
      gStyle->SetOptStat(0001111);
~~~
but instead do:
~~~ {.cpp}
      gStyle->SetOptStat(1111);
~~~
because `0001111` will be taken as an octal number!

<b>WARNING 2:</b> for backward compatibility with older versions
~~~ {.cpp}
      gStyle->SetOptStat(1);
~~~
is taken as:
~~~ {.cpp}
      gStyle->SetOptStat(1111)
~~~
To print only the name of the histogram do:
~~~ {.cpp}
      gStyle->SetOptStat(1000000001);
~~~

<b>NOTE</b> that in case of 2D histograms, when selecting only underflow
(10000) or overflow (100000), the statistics box will show all combinations
of underflow/overflows and not just one single number.

The parameter mode can be any combination of the letters `kKsSiourRmMen`
~~~ {.cpp}
      k :  kurtosis printed
      K :  kurtosis and kurtosis error printed
      s :  skewness printed
      S :  skewness and skewness error printed
      i :  integral of bins printed
      o :  number of overflows printed
      u :  number of underflows printed
      r :  rms printed
      R :  rms and rms error printed
      m :  mean value printed
      M :  mean value mean error values printed
      e :  number of entries printed
      n :  name of histogram is printed
~~~

For example, to print only name of histogram and number of entries do:
~~~ {.cpp}
      gStyle->SetOptStat("ne");
~~~

To print only the name of the histogram do:
~~~ {.cpp}
      gStyle->SetOptStat("n");
~~~

The default value is:
~~~ {.cpp}
      gStyle->SetOptStat("nemr");
~~~

When a histogram is painted, a `TPaveStats` object is created and added
to the list of functions of the histogram. If a `TPaveStats` object
already exists in the histogram list of functions, the existing object is just
updated with the current histogram parameters.

Once a histogram is painted, the statistics box can be accessed using
`h->FindObject("stats")`. In the command line it is enough to do:
~~~ {.cpp}
      Root > h->Draw()
      Root > TPaveStats *st = (TPaveStats*)h->FindObject("stats")
~~~

because after `h->Draw()` the histogram is automatically painted. But
in a script file the painting should be forced using `gPad->Update()`
in order to make sure the statistics box is created:
~~~ {.cpp}
      h->Draw();
      gPad->Update();
      TPaveStats *st = (TPaveStats*)h->FindObject("stats");
~~~

Without `gPad->Update()` the line `h->FindObject("stats")`
returns a null pointer.

When a histogram is drawn with the option "`SAME`", the statistics box
is not drawn. To force the statistics box drawing with the option
"`SAME`", the option "`SAMES`" must be used.
If the new statistics box hides the previous statistics box, one can change
its position with these lines ("`h`" being the pointer to the histogram):
~~~ {.cpp}
      Root > TPaveStats *st = (TPaveStats*)h->FindObject("stats")
      Root > st->SetX1NDC(newx1); //new x start position
      Root > st->SetX2NDC(newx2); //new x end position
~~~

To change the type of information for an histogram with an existing
`TPaveStats` one should do:
~~~ {.cpp}
      st->SetOptStat(mode);
~~~
Where "`mode`" has the same meaning than when calling
`gStyle->SetOptStat(mode)` (see above).

One can delete the statistics box for a histogram `TH1* h` with:
~~~ {.cpp}
      h->SetStats(0)
~~~

and activate it again with:
~~~ {.cpp}
      h->SetStats(1).
~~~

## <a name="PS02"></a> Fit Statistics

The type of information about fit parameters printed in the histogram statistics
box can be selected via the parameter mode. The parameter mode can be
`= pcev`  (default `= 0111`)
~~~ {.cpp}
      p = 1;  print Probability
      c = 1;  print Chisquare/Number of degrees of freedom
      e = 1;  print errors (if e=1, v must be 1)
      v = 1;  print name/values of parameters
~~~
Example:
~~~ {.cpp}
      gStyle->SetOptFit(1011);
~~~
print fit probability, parameter names/values and errors.

  1. When `"v" = 1` is specified, only the non-fixed parameters are
     shown.
  2. When `"v" = 2` all parameters are shown.

Note: `gStyle->SetOptFit(1)` means "default value", so it is equivalent
to `gStyle->SetOptFit(111)`

## <a name="PS03"></a> Statistics box editing

The following example show how to remove and add a line in a statistics box.

Begin_Macro(source)
../../../tutorials/hist/statsEditing.C
End_Macro
*/


const UInt_t kTakeStyle = BIT(17); //see TStyle::SetOptFit/Stat

////////////////////////////////////////////////////////////////////////////////
/// TPaveStats default constructor.

TPaveStats::TPaveStats(): TPaveText()
{
   fParent  = 0;
   fOptFit  = gStyle->GetOptFit();
   fOptStat = gStyle->GetOptStat();
}

////////////////////////////////////////////////////////////////////////////////
/// TPaveStats normal constructor.

TPaveStats::TPaveStats(Double_t x1, Double_t y1,Double_t x2, Double_t  y2, Option_t *option)
           :TPaveText(x1,y1,x2,y2,option)
{
   fParent = 0;
   fOptFit  = gStyle->GetOptFit();
   fOptStat = gStyle->GetOptStat();
   SetFitFormat(gStyle->GetFitFormat());
   SetStatFormat(gStyle->GetStatFormat());
}

////////////////////////////////////////////////////////////////////////////////
/// TPaveStats default destructor.

TPaveStats::~TPaveStats()
{
   if ( fParent && !fParent->TestBit(kInvalidObject)) fParent->RecursiveRemove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the fit option.

Int_t TPaveStats::GetOptFit() const
{
   if (TestBit(kTakeStyle)) return gStyle->GetOptFit();
   return fOptFit;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the stat option.

Int_t TPaveStats::GetOptStat() const
{
   if (TestBit(kTakeStyle)) return gStyle->GetOptStat();
   return fOptStat;
}

////////////////////////////////////////////////////////////////////////////////
/// Save This TPaveStats options in current style.

void TPaveStats::SaveStyle()
{
   gStyle->SetOptFit(fOptFit);
   gStyle->SetOptStat(fOptStat);
   gStyle->SetFitFormat(fFitFormat.Data());
   gStyle->SetStatFormat(fStatFormat.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Change (i.e. set) the format for printing fit parameters in statistics box.

void TPaveStats::SetFitFormat(const char *form)
{
   fFitFormat = form;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the fit option.

void TPaveStats::SetOptFit(Int_t fit)
{
   fOptFit = fit;
   ResetBit(kTakeStyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the stat option.

void TPaveStats::SetOptStat(Int_t stat)
{
   fOptStat = stat;
   ResetBit(kTakeStyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Change (i.e. set) the format for printing statistics.

void TPaveStats::SetStatFormat(const char *form)
{
   fStatFormat = form;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint the pave stat.

void TPaveStats::Paint(Option_t *option)
{
   TPave::ConvertNDCtoPad();
   TPave::PaintPave(fX1,fY1,fX2,fY2,GetBorderSize(),option);

   if (!fLines) return;
   TString typolabel;
   Double_t y2ref = TMath::Max(fY1,fY2);
   Double_t x1ref = TMath::Min(fX1,fX2);
   Double_t x2ref = TMath::Max(fX1,fX2);
   Double_t dx    = TMath::Abs(fX2 - fX1);
   Double_t dy    = TMath::Abs(fY2 - fY1);
   Double_t titlesize=0;
   Double_t textsize = GetTextSize();
   Int_t nlines = GetSize();
   if (nlines == 0) nlines = 5;
   Int_t print_name = fOptStat%10;

   // Evaluate text size as a function of the number of lines
   Double_t y1       = gPad->GetY1();
   Double_t y2       = gPad->GetY2();
   Float_t margin    = fMargin*dx;
   Double_t yspace   = dy/Double_t(nlines);
   Double_t textsave = textsize;
   TObject *line;
   TLatex *latex, *latex_tok;
   TIter next(fLines);
   Double_t longest = 0, titlelength = 0;
   Double_t w, wtok[2];
   char *st, *sl=0;
   if (textsize == 0)  {
      textsize = 0.92*yspace/(y2 - y1);
      titlesize = textsize;
      wtok[0] = 0; wtok[1] = 0;
      while ((line = (TObject*) next())) {
         if (line->IsA() == TLatex::Class()) {
            latex = (TLatex*)line;
            Int_t nchs = strlen(latex->GetTitle());
            sl = new char[nchs+1];
            strlcpy(sl, latex->GetTitle(),nchs+1);
            if (strpbrk(sl, "=") !=0 && print_name == 0) {
               st = strtok(sl, "=");
               Int_t itok = 0;
               while ( st !=0 ) {
                  latex_tok = new TLatex(0.,0.,st);
                  Style_t tfont = latex->GetTextFont();
                  if (tfont == 0) tfont = GetTextFont();
                  latex_tok->SetTextFont(tfont);
                  latex_tok->SetTextSize(textsize);
                  w = latex_tok->GetXsize();
                  if (w > wtok[itok]) wtok[itok] = w;
                  st = strtok(0, "=");
                  ++itok;
                  delete latex_tok;
               }
            } else if (strpbrk(sl, "|") !=0) {
            } else {
               print_name = 0;
               Style_t tfont = latex->GetTextFont();
               if (tfont == 0) latex->SetTextFont(GetTextFont());
               latex->SetTextSize(titlesize);
               titlelength = latex->GetXsize()+2.*margin;
               if (titlelength > 0.98*dx) titlesize *= 0.98*dx/titlelength;
               latex->SetTextFont(tfont);
            }
            delete [] sl; sl = 0;
         }
      }
      longest = wtok[0]+wtok[1]+2.*margin;
      if (longest > 0.98*dx) textsize *= 0.98*dx/longest;
      SetTextSize(textsize);
   } else {
      titlesize = textsize;
   }
   Double_t ytext = y2ref + 0.5*yspace;
   Double_t xtext = 0;
   print_name = fOptStat%10;

   // Iterate over all lines
   // Copy pavetext attributes to line attributes if line attributes not set
   next.Reset();
   while ((line = (TObject*) next())) {
      if (line->IsA() == TLatex::Class()) {
         latex = (TLatex*)line;
         ytext -= yspace;
         Double_t xl    = latex->GetX();
         Double_t yl    = latex->GetY();
         Short_t talign = latex->GetTextAlign();
         Color_t tcolor = latex->GetTextColor();
         Style_t tfont  = latex->GetTextFont();
         Size_t  tsize  = latex->GetTextSize();
         if (tcolor == 0) latex->SetTextColor(GetTextColor());
         if (tfont  == 0) latex->SetTextFont(GetTextFont());
         if (tsize  == 0) latex->SetTextSize(GetTextSize());

         Int_t nchs = strlen(latex->GetTitle());
         sl = new char[nchs+1];
         strlcpy(sl, latex->GetTitle(),nchs+1);
         // Draw all the histogram stats except the 2D under/overflow
         if (strpbrk(sl, "=") !=0 && print_name == 0) {
            st = strtok(sl, "=");
            Int_t halign = 12;
            while ( st !=0 ) {
               typolabel = st;
               latex->SetTextAlign(halign);
               if (halign == 12) xtext = x1ref + margin;
               if (halign == 32) {
                  xtext = x2ref - margin;
                  typolabel = typolabel.Strip();
                  typolabel.ReplaceAll("-","#minus");
               }
               latex->PaintLatex(xtext,ytext,latex->GetTextAngle(),
                                             latex->GetTextSize(),
                                             typolabel.Data());
               st = strtok(0, "=");
               halign = 32;
            }
         // Draw the 2D under/overflow
         } else if (strpbrk(sl, "|") !=0) {
            Double_t yline1 = ytext+yspace/2.;
            Double_t yline2 = ytext-yspace/2.;
            Double_t xline1 = dx/3+x1ref;
            Double_t xline2 = 2*dx/3+x1ref;
            gPad->PaintLine(x1ref,yline1,x2ref,yline1);
            gPad->PaintLine(xline1,yline1,xline1,yline2);
            gPad->PaintLine(xline2,yline1,xline2,yline2);
            st = strtok(sl, "|");
            Int_t theIndex = 0;
            while ( st !=0 ) {
               latex->SetTextAlign(22);
               if (theIndex == 0) xtext = 0.5*(x1ref+xline1);
               if (theIndex == 1) xtext = 0.5*(x1ref+x2ref);
               if (theIndex == 2) xtext = 0.5*(xline2+x2ref);
               typolabel = st;
               typolabel.ReplaceAll("-", "#minus");
               latex->PaintLatex(xtext,ytext,latex->GetTextAngle(),
                                             latex->GetTextSize(),
                                             typolabel.Data());
               theIndex++;
               st = strtok(0, "|");
            }
         // Draw the histogram identifier
         } else {
            print_name = 0;
            latex->SetTextAlign(22);
            xtext = 0.5*(x1ref+x2ref);
            latex->PaintLatex(xtext,ytext,latex->GetTextAngle(),
                                          titlesize,
                                          sl);
            gPad->PaintLine(x1ref,y2ref-yspace,x2ref,y2ref-yspace);
         }
         delete [] sl;

         latex->SetTextAlign(talign);
         latex->SetTextColor(tcolor);
         latex->SetTextFont(tfont);
         latex->SetTextSize(tsize);
         latex->SetX(xl);  //paintlatex modifies fX and fY
         latex->SetY(yl);
      }
   }
   SetTextSize(textsave);

   // if a label create & paint a pavetext title
   if (fLabel.Length() > 0) {
      Double_t x1,x2;
      dy = gPad->GetY2() - gPad->GetY1();
      x1 = x1ref + 0.25*dx;
      x2 = x2ref - 0.25*dx;
      y1 = y2ref - 0.02*dy;
      y2 = y2ref + 0.02*dy;
      TPaveLabel *title = new TPaveLabel(x1,y1,x2,y2,fLabel.Data(),GetDrawOption());
      title->SetFillColor(GetFillColor());
      title->SetTextColor(GetTextColor());
      title->SetTextFont(GetTextFont());
      title->Paint();
      delete title;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out.

void TPaveStats::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   char quote = '"';
   out<<"   "<<std::endl;
   Bool_t saved = gROOT->ClassSaved(TPaveStats::Class());
   if (saved) {
      out<<"   ";
   } else {
      out<<"   "<<ClassName()<<" *";
   }
   if (fOption.Contains("NDC")) {
      out<<"ptstats = new "<<ClassName()<<"("<<fX1NDC<<","<<fY1NDC<<","<<fX2NDC<<","<<fY2NDC
      <<","<<quote<<fOption<<quote<<");"<<std::endl;
   } else {
      out<<"ptstats = new "<<ClassName()<<"("<<fX1<<","<<fY1<<","<<fX2<<","<<fY2
      <<","<<quote<<fOption<<quote<<");"<<std::endl;
   }
   if (strcmp(GetName(),"TPave")) {
      out<<"   ptstats->SetName("<<quote<<GetName()<<quote<<");"<<std::endl;
   }
   if (fBorderSize != 4) {
      out<<"   ptstats->SetBorderSize("<<fBorderSize<<");"<<std::endl;
   }
   SaveFillAttributes(out,"ptstats",19,1001);
   SaveLineAttributes(out,"ptstats",1,1,1);
   SaveTextAttributes(out,"ptstats",22,0,1,62,0);
   SaveLines(out,"ptstats",saved);
   out<<"   ptstats->SetOptStat("<<GetOptStat()<<");"<<std::endl;
   out<<"   ptstats->SetOptFit("<<GetOptFit()<<");"<<std::endl;
   out<<"   ptstats->Draw();"<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TPaveStats.

void TPaveStats::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         R__b.ReadClassBuffer(TPaveStats::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TPaveText::Streamer(R__b);
      R__b >> fOptFit;
      R__b >> fOptStat;
      if (R__v > 1 || R__b.GetVersionOwner() == 22304) {
         fFitFormat.Streamer(R__b);
         fStatFormat.Streamer(R__b);
      } else {
         SetFitFormat();
         SetStatFormat();
      }
      R__b.CheckByteCount(R__s, R__c, TPaveStats::IsA());
      //====end of old versions

   } else {
      R__b.WriteClassBuffer(TPaveStats::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Replace current attributes by current style.

void TPaveStats::UseCurrentStyle()
{
   if (gStyle->IsReading()) {
      SetOptStat(gStyle->GetOptStat());
      SetOptFit(gStyle->GetOptFit());
      SetStatFormat(gStyle->GetStatFormat());
      SetFitFormat(gStyle->GetFitFormat());
      SetBorderSize(gStyle->GetStatBorderSize());
      SetFillColor(gStyle->GetStatColor());
      SetFillStyle(gStyle->GetStatStyle());
      SetTextFont(gStyle->GetStatFont());
      SetTextSize(gStyle->GetStatFontSize());
      SetTextColor(gStyle->GetStatTextColor());
      SetX2NDC(gStyle->GetStatX());
      SetY2NDC(gStyle->GetStatY());
      SetX1NDC(gStyle->GetStatX()-gStyle->GetStatW());
      SetY1NDC(gStyle->GetStatY()-gStyle->GetStatH());
   } else {
      gStyle->SetOptStat(GetOptStat());
      gStyle->SetOptFit(GetOptFit());
      gStyle->SetStatFormat(GetStatFormat());
      gStyle->SetFitFormat(GetFitFormat());
      gStyle->SetStatBorderSize(GetBorderSize());
      gStyle->SetTextColor(GetTextColor());
      gStyle->SetStatColor(GetFillColor());
      gStyle->SetStatStyle(GetFillStyle());
      gStyle->SetTextFont(GetTextFont());
      gStyle->SetStatFontSize(GetTextSize());
      gStyle->SetStatTextColor(GetTextColor());
      gStyle->SetStatX(GetX2NDC());
      gStyle->SetStatY(GetY2NDC());
      gStyle->SetStatW(GetX2NDC()-GetX1NDC());
      gStyle->SetStatH(GetY2NDC()-GetY1NDC());
   }
}
