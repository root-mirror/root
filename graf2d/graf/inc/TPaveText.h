// @(#)root/graf:$Id$
// Author: Rene Brun   20/10/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPaveText
#define ROOT_TPaveText


#include "TPave.h"
#include "TText.h"

class TLine;

class TPaveText : public TPave, public TAttText {

protected:
   TString    fLabel;          ///< Label written at the top of the pavetext
   Int_t      fLongest;        ///< Length of the longest line
   Float_t    fMargin;         ///< Text margin
   TList      *fLines;         ///< List of labels


public:
   TPaveText();
   TPaveText(Double_t x1, Double_t y1,Double_t x2 ,Double_t y2, Option_t *option="br");
   TPaveText(const TPaveText &pavetext);
   virtual ~TPaveText();
   TPaveText& operator=(const TPaveText&);

   virtual TBox    *AddBox(Double_t x1, Double_t y1, Double_t x2, Double_t y2);
   virtual TLine   *AddLine(Double_t x1=0, Double_t y1=0, Double_t x2=0, Double_t y2=0);
   virtual TText   *AddText(Double_t x1, Double_t y1, const char *label);
   virtual TText   *AddText(const char *label);
   virtual void     Clear(Option_t *option="");  // *MENU*
   virtual void     DeleteText(); // *MENU*
   virtual void     Draw(Option_t *option="");
   virtual void     DrawFile(const char *filename, Option_t *option="");
   virtual void     EditText(); // *MENU*
   const char      *GetLabel() const {return fLabel.Data();}
   virtual TText   *GetLine(Int_t number) const;
   virtual TText   *GetLineWith(const char *text) const;
   virtual TList   *GetListOfLines() const {return fLines;}
   Float_t          GetMargin() const {return fMargin;}
   virtual TObject *GetObject(Double_t &ymouse, Double_t &yobj) const;
   virtual Int_t    GetSize() const;
   virtual void     InsertLine(); // *MENU*
   virtual void     InsertText(const char *label); // *MENU*
   virtual void     Paint(Option_t *option="");
   virtual void     PaintPrimitives(Int_t mode);
   virtual void     Print(Option_t *option="") const;
   virtual void     ReadFile(const char *filename, Option_t *option="", Int_t nlines=50, Int_t fromline=0); // *MENU*
   virtual void     SaveLines(std::ostream &out, const char *name, Bool_t saved);
   virtual void     SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void     SetAllWith(const char *text, Option_t *option, Double_t value); // *MENU*
   virtual void     SetLabel(const char *label) {fLabel = label;} // *MENU*
   virtual void     SetMargin(Float_t margin=0.05) {fMargin=margin;} // *MENU*
   virtual void     UseCurrentStyle();

   ClassDef(TPaveText,2)  //PaveText. A Pave with several lines of text.
};

#endif

