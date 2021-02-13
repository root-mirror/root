// @(#)root/hist:$Id$
// Author: Rene Brun   12/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TAxis
#define ROOT_TAxis


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TAxis                                                                //
//                                                                      //
// Axis class.                                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TAttAxis.h"
#include "TArrayD.h"

class THashList;

class TAxis : public TNamed, public TAttAxis {

private:
   Int_t        fNbins;          //Number of bins
   Double_t     fXmin;           //low edge of first bin
   Double_t     fXmax;           //upper edge of last bin
   TArrayD      fXbins;          //Bin edges array in X
   Int_t        fFirst;          //first bin to display
   Int_t        fLast;           //last bin to display
   UShort_t     fBits2;          //second bit status word
   Bool_t       fTimeDisplay;    //on/off displaying time values instead of numerics
   TString      fTimeFormat;     //Date&time format, ex: 09/12/99 12:34:00
   TObject     *fParent;         //!Object owning this axis
   THashList   *fLabels;         //List of labels
   TList       *fModLabs;        //List of modified labels

   // TAxis extra status bits (stored in fBits2)
   enum {
      kAlphanumeric = BIT(0),   // axis is alphanumeric
      kCanExtend = BIT(1),      // axis can be extended
      kNotAlpha = BIT(2)    // axis is forced to be not alphanumeric
   };

   Bool_t       HasBinWithoutLabel() const;

public:
   // TAxis status bits
   enum EStatusBits {
      kDecimals      = BIT(7),
      kTickPlus      = BIT(9),
      kTickMinus     = BIT(10),
      kAxisRange     = BIT(11),
      kCenterTitle   = BIT(12),
      kCenterLabels  = BIT(14), //bit 13 is used by TObject
      kRotateTitle   = BIT(15),
      kPalette       = BIT(16),
      kNoExponent    = BIT(17),
      kLabelsHori    = BIT(18),
      kLabelsVert    = BIT(19),
      kLabelsDown    = BIT(20),
      kLabelsUp      = BIT(21),
      kIsInteger     = BIT(22),
      kMoreLogLabels = BIT(23)
   };

   TAxis();
   TAxis(Int_t nbins, Double_t xmin, Double_t xmax);
   TAxis(Int_t nbins, const Double_t *xbins);
   TAxis(const TAxis &axis);
   virtual ~TAxis();
   TAxis& operator=(const TAxis&);

   Bool_t     CanExtend() const { return (fBits2 & kCanExtend);  }
   Bool_t     CanBeAlphanumeric() { return !(fBits2 & kNotAlpha); }
   Bool_t     IsAlphanumeric() const { return fBits2 & kAlphanumeric; }
   void       SetAlphanumeric(Bool_t alphanumeric = kTRUE);
   void       SetCanExtend(Bool_t canExtend) { fBits2 = canExtend ? (fBits2 | kCanExtend) : (fBits2 & ~kCanExtend); }
   void       SetNoAlphanumeric(Bool_t noalpha = kTRUE) {
      fBits2 = noalpha ? (fBits2 | kNotAlpha) : (fBits2 & ~kNotAlpha);
      if (IsAlphanumeric() ) {
         SetCanExtend(kFALSE);
         SetAlphanumeric(kFALSE);
      }
   }
   void               CenterLabels(Bool_t center=kTRUE);
   void               CenterTitle(Bool_t center=kTRUE);
   const char        *ChooseTimeFormat(Double_t axislength=0);
   virtual void       Copy(TObject &axis) const;
   virtual void       Delete(Option_t * /*option*/ ="") { }
   virtual Int_t      DistancetoPrimitive(Int_t px, Int_t py);
   virtual TObject   *DrawClone(Option_t * /*option*/ ="") const {return 0;}
   virtual void       ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual Int_t      FindBin(Double_t x);
   virtual Int_t      FindBin(Double_t x) const { return FindFixBin(x); }
   virtual Int_t      FindBin(const char *label);
   virtual Int_t      FindFixBin(Double_t x) const;
   virtual Int_t      FindFixBin(const char *label) const;
   virtual Double_t   GetBinCenter(Int_t bin) const;
   virtual Double_t   GetBinCenterLog(Int_t bin) const;
   const char        *GetBinLabel(Int_t bin) const;
   virtual Double_t   GetBinLowEdge(Int_t bin) const;
   virtual Double_t   GetBinUpEdge(Int_t bin) const;
   virtual Double_t   GetBinWidth(Int_t bin) const;
   virtual void       GetCenter(Double_t *center) const;
           Bool_t     GetCenterLabels() const { return TestBit(kCenterLabels); }
           Bool_t     GetCenterTitle() const { return TestBit(kCenterTitle); }
           Bool_t     GetDecimals() const { return TestBit(kDecimals); }
   THashList         *GetLabels() const { return fLabels; }
   TList             *GetModifiedLabels() const { return fModLabs; }
   virtual void       GetLowEdge(Double_t *edge) const;
           Bool_t     GetMoreLogLabels() const { return TestBit(kMoreLogLabels); }
           Int_t      GetNbins() const { return fNbins; }
           Bool_t     GetNoExponent() const { return TestBit(kNoExponent); }
   virtual TObject   *GetParent() const {return fParent;}
           Bool_t     GetRotateTitle() const { return TestBit(kRotateTitle); }
   virtual const char *GetTicks() const;
   virtual Bool_t     GetTimeDisplay() const {return fTimeDisplay;}
   virtual const char *GetTimeFormat() const {return fTimeFormat.Data();}
   virtual const char *GetTimeFormatOnly() const;
   const char        *GetTitle() const {return fTitle.Data();}
   const TArrayD     *GetXbins() const {return &fXbins;}
           Int_t      GetFirst() const;
           Int_t      GetLast() const;
           Double_t   GetXmin() const {return fXmin;}
           Double_t   GetXmax() const {return fXmax;}
   virtual void       ImportAttributes(const TAxis *axis);
   Bool_t             IsVariableBinSize() const {
                         // true if axis has variable bin sizes, false otherwise
                         return (fXbins.GetSize() != 0);
                      }
   virtual void       LabelsOption(Option_t *option="h");  // *MENU*
           void       RotateTitle(Bool_t rotate=kTRUE); // *TOGGLE* *GETTER=GetRotateTitle
   virtual void       SaveAttributes(std::ostream &out, const char *name, const char *subname);
   virtual void       Set(Int_t nbins, Double_t xmin, Double_t xmax);
   virtual void       Set(Int_t nbins, const Float_t *xbins);
   virtual void       Set(Int_t nbins, const Double_t *xbins);
   virtual void       SetBinLabel(Int_t bin, const char *label);
           void       SetDecimals(Bool_t dot = kTRUE); // *TOGGLE* *GETTER=GetDecimals
   virtual void       SetDefaults();
   virtual void       SetDrawOption(Option_t * /*option*/ ="") { }
   void               ChangeLabel(Int_t labNum=0, Double_t labAngle = -1.,
                                  Double_t labSize = -1., Int_t labAlign = -1,
                                  Int_t labColor = -1 , Int_t labFont = -1,
                                  TString labText = ""); // *MENU*
   virtual void       SetLimits(Double_t xmin, Double_t xmax) { /* set axis limits */ fXmin = xmin; fXmax = xmax; } // *MENU*
           void       SetMoreLogLabels(Bool_t more=kTRUE);  // *TOGGLE* *GETTER=GetMoreLogLabels
           void       SetNoExponent(Bool_t noExponent=kTRUE);  // *TOGGLE* *GETTER=GetNoExponent
   virtual void       SetParent(TObject *obj) {fParent = obj;}
   virtual void       SetRange(Int_t first=0, Int_t last=0);  // *MENU*
   virtual void       SetRangeUser(Double_t ufirst, Double_t ulast);  // *MENU*
   virtual void       SetTicks(Option_t *option="+"); // *MENU*
   virtual void       SetTimeDisplay(Int_t value) {fTimeDisplay = (value != 0);} // *TOGGLE*
   virtual void       SetTimeFormat(const char *format="");  // *MENU*
   virtual void       SetTimeOffset(Double_t toffset, Option_t *option="local");
   virtual void       UnZoom();  // *MENU*
   virtual void       ZoomOut(Double_t factor=0, Double_t offset=0);  // *MENU*

   ClassDef(TAxis,10)  //Axis class
};

////////////////////////////////////////////////////////////////////////////////
/// Center axis labels. If center = kTRUE axis labels will be centered
/// (horizontal axes only) on the bin center. Default is to center on the primary tick marks
/// This option does not make sense if there are more bins than tick marks

inline void TAxis::CenterLabels(Bool_t center)
{
   SetBit(kCenterLabels, center);
}

////////////////////////////////////////////////////////////////////////////////
/// Center axis title. If center = kTRUE axis title will be centered
/// default is right adjusted

inline void TAxis::CenterTitle(Bool_t center)
{
   SetBit(kCenterTitle, center);
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate title by 180 degrees. By default the title is drawn right adjusted.
/// If rotate is TRUE, the title is left adjusted at the end of the axis and rotated by 180 degrees

inline void TAxis::RotateTitle(Bool_t rotate)
{
   SetBit(kRotateTitle, rotate);
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the decimals flag
/// By default, blank characters are stripped, and then the label is correctly aligned.
/// If the dot is the last character of the string, it is also stripped, unless this option is specified.

inline void TAxis::SetDecimals(Bool_t dot) {
   SetBit(kDecimals, dot);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the kMoreLogLabels bit flag
/// When this option is selected more labels are drawn when in log scale and there is a small number
/// of decades  (<3).
/// The flag (in fBits) is passed to the drawing function TGaxis::PaintAxis

inline void TAxis::SetMoreLogLabels(Bool_t more)
{
    SetBit(kMoreLogLabels, more);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the NoExponent flag
/// By default, an exponent of the form 10^N is used when the label value are either all very small or very large.
/// The flag (in fBits) is passed to the drawing function TGaxis::PaintAxis

inline void TAxis::SetNoExponent(Bool_t noExponent)
{
   SetBit(kNoExponent, noExponent);
}


#endif
