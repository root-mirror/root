// @(#)root/graf:$Id$
// Author: Rene Brun   15/03/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TPaveStats
#define ROOT_TPaveStats


#include "TPaveText.h"
#include "TVirtualPaveStats.h"

class TPaveStats : public TPaveText, public TVirtualPaveStats {

protected:
   Int_t         fOptFit;            ///< option Fit
   Int_t         fOptStat;           ///< option Stat
   TString       fFitFormat;         ///< Printing format for fit parameters
   TString       fStatFormat;        ///< Printing format for stats
   TObject      *fParent;            ///< owner of this TPaveStats

public:
   TPaveStats();
   TPaveStats(Double_t x1, Double_t y1,Double_t x2 ,Double_t y2, Option_t *option="br");
   virtual ~TPaveStats();
   virtual TBox    *AddBox(Double_t , Double_t , Double_t , Double_t) {return 0;}
   virtual TLine   *AddLine(Double_t , Double_t , Double_t, Double_t) {return 0;}
   virtual void     DeleteText() { }
   virtual void     EditText() { }
   virtual const char  *GetFitFormat()  const {return fFitFormat.Data();}
   virtual const char  *GetStatFormat() const {return fStatFormat.Data();}
   Int_t            GetOptFit() const;
   Int_t            GetOptStat() const;
   virtual TObject *GetParent() const {return fParent;}
   virtual void     Paint(Option_t *option="");
   virtual void     InsertText(const char *) { }
   virtual void     InsertLine() { }
   virtual void     ReadFile(const char *, Option_t *, Int_t, Int_t) { }
   virtual void     SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void     SaveStyle(); // *MENU*
   virtual void     SetAllWith(const char *, Option_t *, Double_t) { }
   virtual void     SetMargin(Float_t) { }
   virtual void     SetFitFormat(const char *format="5.4g");    // *MENU*
   virtual void     SetStatFormat(const char *format="6.4g");   // *MENU*
   void             SetOptFit(Int_t fit=1);                     // *MENU*
   void             SetOptStat(Int_t stat=1);                   // *MENU*
   virtual void     SetParent(TObject*obj) {fParent = obj;}
   virtual void     UseCurrentStyle();

   ClassDef(TPaveStats,5)  //A special TPaveText to draw histogram statistics.
};

#endif
