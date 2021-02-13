// @(#)root/graf:$Id$
// Author: Rene Brun, Olivier Couet   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGaxis
#define ROOT_TGaxis

#include "TLine.h"
#include "TAttText.h"
#include "TString.h"

class TF1;
class TAxis;
class TLatex;

class TGaxis : public TLine, public TAttText {

protected:

   Double_t   fWmin;                ///< Lowest value on the axis
   Double_t   fWmax;                ///< Highest value on the axis
   Float_t    fGridLength;          ///< Length of the grid in NDC
   Float_t    fTickSize;            ///< Size of primary tick mark in NDC
   Float_t    fLabelOffset;         ///< Offset of label wrt axis
   Float_t    fLabelSize;           ///< Size of labels in NDC
   Float_t    fTitleOffset;         ///< Offset of title wrt axis
   Float_t    fTitleSize;           ///< Size of title in NDC
   Int_t      fNdiv;                ///< Number of divisions
   Int_t      fLabelColor;          ///< Color for labels
   Int_t      fLabelFont;           ///< Font for labels
   Int_t      fNModLabs;            ///< Number of modified labels
   TString    fChopt;               ///< Axis options
   TString    fName;                ///< Axis name
   TString    fTitle;               ///< Axis title
   TString    fTimeFormat;          ///< Time format, ex: 09/12/99 12:34:00
   TString    fFunctionName;        ///< Name of mapping function pointed by fFunction
   TF1       *fFunction;            ///<! Pointer to function computing axis values
   TAxis     *fAxis;                ///<! Pointer to original TAxis axis (if any)
   TList     *fModLabs;             ///<  List of modified labels.

   static Int_t fgMaxDigits;        ///<! Number of digits above which the 10>N notation is used
   static Float_t fXAxisExpXOffset; ///<! Exponent X offset for the X axis
   static Float_t fXAxisExpYOffset; ///<! Exponent Y offset for the X axis
   static Float_t fYAxisExpXOffset; ///<! Exponent X offset for the Y axis
   static Float_t fYAxisExpYOffset; ///<! Exponent Y offset for the Y axis

   TGaxis(const TGaxis&);
   TGaxis& operator=(const TGaxis&);

public:

   TGaxis();
   TGaxis(Double_t xmin,Double_t ymin,Double_t xmax,Double_t ymax,
          Double_t wmin,Double_t wmax,Int_t ndiv=510, Option_t *chopt="",
          Double_t gridlength = 0);
   TGaxis(Double_t xmin,Double_t ymin,Double_t xmax,Double_t ymax,
          const char *funcname, Int_t ndiv=510, Option_t *chopt="",
          Double_t gridlength = 0);
   virtual ~TGaxis();

   virtual void        AdjustBinSize(Double_t A1,  Double_t A2,  Int_t nold
                                    ,Double_t &BinLow, Double_t &BinHigh, Int_t &nbins, Double_t &BinWidth);
   virtual void        CenterLabels(Bool_t center=kTRUE);
   virtual void        CenterTitle(Bool_t center=kTRUE);
   void                ChangeLabelAttributes(Int_t i, Int_t nlabels, TLatex* t, char* c);
   virtual void        DrawAxis(Double_t xmin,Double_t ymin,Double_t xmax,Double_t ymax,
                                Double_t wmin,Double_t wmax,Int_t ndiv=510, Option_t *chopt="",
                                Double_t gridlength = 0);
   Float_t             GetGridLength() const   {return fGridLength;}
   TF1                *GetFunction() const     {return fFunction;}
   Int_t               GetLabelColor() const   {return fLabelColor;}
   Int_t               GetLabelFont() const    {return fLabelFont;}
   Float_t             GetLabelOffset() const  {return fLabelOffset;}
   Float_t             GetLabelSize() const    {return fLabelSize;}
   Float_t             GetTitleOffset() const  {return fTitleOffset;}
   Float_t             GetTitleSize() const    {return fTitleSize;}
   virtual const char *GetName() const  {return fName.Data();}
   virtual const char *GetOption() const {return fChopt.Data();}
   virtual const char *GetTitle() const {return fTitle.Data();}
   static Int_t        GetMaxDigits();
   Int_t               GetNdiv() const         {return fNdiv;}
   Double_t            GetWmin() const         {return fWmin;}
   Double_t            GetWmax()  const        {return fWmax;}
   Float_t             GetTickSize() const     {return fTickSize;}
   virtual void        ImportAxisAttributes(TAxis *axis);
   void                LabelsLimits(const char *label, Int_t &first, Int_t &last);
   virtual void        Paint(Option_t *chopt="");
   virtual void        PaintAxis(Double_t xmin,Double_t ymin,Double_t xmax,Double_t ymax,
                                 Double_t &wmin,Double_t &wmax,Int_t &ndiv, Option_t *chopt="",
                                 Double_t gridlength = 0, Bool_t drawGridOnly = kFALSE);
   virtual void        Rotate(Double_t X,  Double_t Y,  Double_t CFI, Double_t SFI
                             ,Double_t XT, Double_t YT, Double_t &U,   Double_t &V);
   void                ResetLabelAttributes(TLatex* t);
   virtual void        SavePrimitive(std::ostream &out, Option_t *option = "");
   void                SetFunction(const char *funcname="");
   void                SetOption(Option_t *option="");
   void                SetLabelColor(Int_t labelcolor) {fLabelColor = labelcolor;} // *MENU*
   void                SetLabelFont(Int_t labelfont) {fLabelFont = labelfont;} // *MENU*
   void                SetLabelOffset(Float_t labeloffset) {fLabelOffset = labeloffset;} // *MENU*
   void                SetLabelSize(Float_t labelsize) {fLabelSize = labelsize;} // *MENU*
   void                ChangeLabel(Int_t labNum=0, Double_t labAngle = -1.,
                                          Double_t labSize = -1., Int_t labAlign = -1,
                                          Int_t labColor = -1 , Int_t labFont = -1,
                                          TString labText = ""); // *MENU*
   static void         SetMaxDigits(Int_t maxd=5);
   virtual void        SetName(const char *name); // *MENU*
   virtual void        SetNdivisions(Int_t ndiv) {fNdiv = ndiv;} // *MENU*
   virtual void        SetMoreLogLabels(Bool_t more=kTRUE);  // *MENU*
   virtual void        SetNoExponent(Bool_t noExponent=kTRUE);  // *MENU*
   virtual void        SetDecimals(Bool_t dot=kTRUE);  // *MENU*
   void                SetTickSize(Float_t ticksize) {fTickSize = ticksize;} // *MENU*
   void                SetTickLength(Float_t ticklength) {SetTickSize(ticklength);}
   void                SetGridLength(Float_t gridlength) {fGridLength = gridlength;}
   void                SetTimeFormat(const char *tformat);
   void                SetTimeOffset(Double_t toffset, Option_t *option="local");
   virtual void        SetTitle(const char *title=""); // *MENU*
   void                SetTitleOffset(Float_t titleoffset=1) {fTitleOffset = titleoffset;} // *MENU*
   void                SetTitleSize(Float_t titlesize) {fTitleSize = titlesize;} // *MENU*
   void                SetTitleFont(Int_t titlefont) {SetTextFont(titlefont);} // *MENU*
   void                SetTitleColor(Int_t titlecolor) {SetTextColor(titlecolor);} // *MENU*
   void                SetWmin(Double_t wmin) {fWmin = wmin;}
   void                SetWmax(Double_t wmax) {fWmax = wmax;}
   static void         SetExponentOffset(Float_t xoff=0., Float_t yoff=0., Option_t *axis="xy");

   ClassDef(TGaxis,6)  //Graphics axis
};

#endif
