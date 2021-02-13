// @(#)root/hist:$Id$
// Author: Rene Brun   26/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TH1
#define ROOT_TH1


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TH1                                                                  //
//                                                                      //
// 1-Dim histogram base class.                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TAxis.h"

#include "TAttLine.h"

#include "TAttFill.h"

#include "TAttMarker.h"

#include "TArrayC.h"
#include "TArrayS.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TArrayD.h"
#include "Foption.h"

#include "TVectorFfwd.h"
#include "TVectorDfwd.h"

#include "TFitResultPtr.h"

#include <cfloat>
#include <string>

class TF1;
class TH1D;
class TBrowser;
class TDirectory;
class TList;
class TCollection;
class TVirtualFFT;
class TVirtualHistPainter;
class TRandom;


class TH1 : public TNamed, public TAttLine, public TAttFill, public TAttMarker {

public:

   // enumeration specifying type of statistics for bin errors
   enum  EBinErrorOpt {
         kNormal = 0,    ///< errors with Normal (Wald) approximation: errorUp=errorLow= sqrt(N)
         kPoisson = 1 ,  ///< errors from Poisson interval at 68.3% (1 sigma)
         kPoisson2 = 2   ///< errors from Poisson interval at 95% CL (~ 2 sigma)
   };

   // enumeration specifying which axes can be extended
   enum {
      kNoAxis  = 0,      ///< NOTE: Must always be 0 !!!
      kXaxis = BIT(0),
      kYaxis = BIT(1),
      kZaxis = BIT(2),
      kAllAxes = kXaxis | kYaxis | kZaxis
   };

   /// Enumeration specifying the way to treat statoverflow
   enum  EStatOverflows {
         kIgnore = 0,   ///< Override global flag ignoring the overflows
         kConsider = 1, ///< Override global flag considering the overflows
         kNeutral = 2,  ///< Adapt to the global flag
   };

   friend class TH1Merger;

protected:
    Int_t         fNcells;          ///< number of bins(1D), cells (2D) +U/Overflows
    TAxis         fXaxis;           ///< X axis descriptor
    TAxis         fYaxis;           ///< Y axis descriptor
    TAxis         fZaxis;           ///< Z axis descriptor
    Short_t       fBarOffset;       ///< (1000*offset) for bar charts or legos
    Short_t       fBarWidth;        ///< (1000*width) for bar charts or legos
    Double_t      fEntries;         ///< Number of entries
    Double_t      fTsumw;           ///< Total Sum of weights
    Double_t      fTsumw2;          ///< Total Sum of squares of weights
    Double_t      fTsumwx;          ///< Total Sum of weight*X
    Double_t      fTsumwx2;         ///< Total Sum of weight*X*X
    Double_t      fMaximum;         ///< Maximum value for plotting
    Double_t      fMinimum;         ///< Minimum value for plotting
    Double_t      fNormFactor;      ///< Normalization factor
    TArrayD       fContour;         ///< Array to display contour levels
    TArrayD       fSumw2;           ///< Array of sum of squares of weights
    TString       fOption;          ///< histogram options
    TList        *fFunctions;       ///<->Pointer to list of functions (fits and user)
    Int_t         fBufferSize;      ///< fBuffer size
    Double_t     *fBuffer;          ///<[fBufferSize] entry buffer
    TDirectory   *fDirectory;       ///<!Pointer to directory holding this histogram
    Int_t         fDimension;       ///<!Histogram dimension (1, 2 or 3 dim)
    Double_t     *fIntegral;        ///<!Integral of bins used by GetRandom
    TVirtualHistPainter *fPainter;  ///<!pointer to histogram painter
    EBinErrorOpt  fBinStatErrOpt;   ///< option for bin statistical errors
    EStatOverflows fStatOverflows;  ///< per object flag to use under/overflows in statistics
    static Int_t  fgBufferSize;     ///<!default buffer size for automatic histograms
    static Bool_t fgAddDirectory;   ///<!flag to add histograms to the directory
    static Bool_t fgStatOverflows;  ///<!flag to use under/overflows in statistics
    static Bool_t fgDefaultSumw2;   ///<!flag to call TH1::Sumw2 automatically at histogram creation time

public:
   static Int_t FitOptionsMake(Option_t *option, Foption_t &Foption);

private:
   void    Build();

   TH1(const TH1&);
   TH1& operator=(const TH1&); // Not implemented


protected:
   TH1();
   TH1(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup);
   TH1(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins);
   TH1(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins);

   Int_t            AxisChoice(Option_t *axis) const;
   virtual Int_t    BufferFill(Double_t x, Double_t w);
   virtual Bool_t   FindNewAxisLimits(const TAxis* axis, const Double_t point, Double_t& newMin, Double_t &newMax);
   virtual void     SavePrimitiveHelp(std::ostream &out, const char *hname, Option_t *option = "");
   static Bool_t    RecomputeAxisLimits(TAxis& destAxis, const TAxis& anAxis);
   static Bool_t    SameLimitsAndNBins(const TAxis& axis1, const TAxis& axis2);
   Bool_t   IsEmpty() const;
   UInt_t GetAxisLabelStatus() const;

   inline static Double_t AutoP2GetPower2(Double_t x, Bool_t next = kTRUE);
   inline static Int_t AutoP2GetBins(Int_t n);
   virtual Int_t AutoP2FindLimits(Double_t min, Double_t max);

   virtual Double_t DoIntegral(Int_t ix1, Int_t ix2, Int_t iy1, Int_t iy2, Int_t iz1, Int_t iz2, Double_t & err,
                               Option_t * opt, Bool_t doerr = kFALSE) const;

   virtual void     DoFillN(Int_t ntimes, const Double_t *x, const Double_t *w, Int_t stride=1);
   Bool_t    GetStatOverflowsBehaviour() const { return EStatOverflows::kNeutral == fStatOverflows ? fgStatOverflows : EStatOverflows::kConsider == fStatOverflows; }

   static bool CheckAxisLimits(const TAxis* a1, const TAxis* a2);
   static bool CheckBinLimits(const TAxis* a1, const TAxis* a2);
   static bool CheckBinLabels(const TAxis* a1, const TAxis* a2);
   static bool CheckEqualAxes(const TAxis* a1, const TAxis* a2);
   static bool CheckConsistentSubAxes(const TAxis *a1, Int_t firstBin1, Int_t lastBin1, const TAxis *a2, Int_t firstBin2=0, Int_t lastBin2=0);
   static bool CheckConsistency(const TH1* h1, const TH1* h2);

public:
   // TH1 status bits
   enum EStatusBits {
      kNoStats     = BIT(9),   ///< don't draw stats box
      kUserContour = BIT(10),  ///< user specified contour levels
      // kCanRebin    = BIT(11), ///< FIXME DEPRECATED - to be removed, replaced by SetCanExtend / CanExtendAllAxes
      kLogX        = BIT(15),  ///< X-axis in log scale
      kIsZoomed   = BIT(16),   ///< bit set when zooming on Y axis
      kNoTitle     = BIT(17),  ///< don't draw the histogram title
      kIsAverage   = BIT(18),  ///< Bin contents are average (used by Add)
      kIsNotW      = BIT(19),  ///< Histogram is forced to be not weighted even when the histogram is filled with weighted
                               /// different than 1.
      kAutoBinPTwo = BIT(20),  ///< Use Power(2)-based algorithm for autobinning
      kIsHighlight = BIT(21)   ///< bit set if histo is highlight
   };
   // size of statistics data (size of  array used in GetStats()/ PutStats )
   // s[0]  = sumw       s[1]  = sumw2
   // s[2]  = sumwx      s[3]  = sumwx2
   // s[4]  = sumwy      s[5]  = sumwy2   s[6]  = sumwxy
   // s[7]  = sumwz      s[8]  = sumwz2   s[9]  = sumwxz   s[10]  = sumwyz
   // s[11] = sumwt      s[12] = sumwt2                 (11 and 12 used only by TProfile3D)
   enum {
      kNstat       = 13  // size of statistics data (up to TProfile3D)
   };


   virtual ~TH1();

   virtual Bool_t   Add(TF1 *h1, Double_t c1=1, Option_t *option="");
   virtual Bool_t   Add(const TH1 *h1, Double_t c1=1);
   virtual Bool_t   Add(const TH1 *h, const TH1 *h2, Double_t c1=1, Double_t c2=1); // *MENU*
   virtual void     AddBinContent(Int_t bin);
   virtual void     AddBinContent(Int_t bin, Double_t w);
   static  void     AddDirectory(Bool_t add=kTRUE);
   static  Bool_t   AddDirectoryStatus();
   virtual void     Browse(TBrowser *b);
   virtual Bool_t   CanExtendAllAxes() const;
   virtual Double_t Chi2Test(const TH1* h2, Option_t *option = "UU", Double_t *res = 0) const;
   virtual Double_t Chi2TestX(const TH1* h2, Double_t &chi2, Int_t &ndf, Int_t &igood,Option_t *option = "UU",  Double_t *res = 0) const;
   virtual Double_t Chisquare(TF1 * f1, Option_t *option = "") const;
   virtual void     ClearUnderflowAndOverflow();
   virtual Double_t ComputeIntegral(Bool_t onlyPositive = false);
   TObject*         Clone(const char* newname=0) const;
   virtual void     Copy(TObject &hnew) const;
   virtual void     DirectoryAutoAdd(TDirectory *);
   virtual Int_t    DistancetoPrimitive(Int_t px, Int_t py);
   virtual Bool_t   Divide(TF1 *f1, Double_t c1=1);
   virtual Bool_t   Divide(const TH1 *h1);
   virtual Bool_t   Divide(const TH1 *h1, const TH1 *h2, Double_t c1=1, Double_t c2=1, Option_t *option=""); // *MENU*
   virtual void     Draw(Option_t *option="");
   virtual TH1     *DrawCopy(Option_t *option="", const char * name_postfix = "_copy") const;
   virtual TH1     *DrawNormalized(Option_t *option="", Double_t norm=1) const;
   virtual void     DrawPanel(); // *MENU*
   virtual Int_t    BufferEmpty(Int_t action=0);
   virtual void     Eval(TF1 *f1, Option_t *option="");
   virtual void     ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual void     ExtendAxis(Double_t x, TAxis *axis);
   virtual TH1     *FFT(TH1* h_output, Option_t *option);
   virtual Int_t    Fill(Double_t x);
   virtual Int_t    Fill(Double_t x, Double_t w);
   virtual Int_t    Fill(const char *name, Double_t w);
   virtual void     FillN(Int_t ntimes, const Double_t *x, const Double_t *w, Int_t stride=1);
   virtual void     FillN(Int_t, const Double_t *, const Double_t *, const Double_t *, Int_t) {;}
   virtual void     FillRandom(const char *fname, Int_t ntimes=5000, TRandom * rng = nullptr);
   virtual void     FillRandom(TH1 *h, Int_t ntimes=5000, TRandom * rng = nullptr);
   virtual Int_t    FindBin(Double_t x, Double_t y=0, Double_t z=0);
   virtual Int_t    FindFixBin(Double_t x, Double_t y=0, Double_t z=0) const;
   virtual Int_t    FindFirstBinAbove(Double_t threshold=0, Int_t axis=1, Int_t firstBin=1, Int_t lastBin=-1) const;
   virtual Int_t    FindLastBinAbove (Double_t threshold=0, Int_t axis=1, Int_t firstBin=1, Int_t lastBin=-1) const;
   virtual TObject *FindObject(const char *name) const;
   virtual TObject *FindObject(const TObject *obj) const;
   virtual TFitResultPtr    Fit(const char *formula ,Option_t *option="" ,Option_t *goption="", Double_t xmin=0, Double_t xmax=0); // *MENU*
   virtual TFitResultPtr    Fit(TF1 *f1 ,Option_t *option="" ,Option_t *goption="", Double_t xmin=0, Double_t xmax=0);
   virtual void     FitPanel(); // *MENU*
   TH1             *GetAsymmetry(TH1* h2, Double_t c2=1, Double_t dc2=0);
   Int_t            GetBufferLength() const {return fBuffer ? (Int_t)fBuffer[0] : 0;}
   Int_t            GetBufferSize  () const {return fBufferSize;}
   const   Double_t *GetBuffer() const {return fBuffer;}
   static  Int_t    GetDefaultBufferSize();
   virtual Double_t *GetIntegral();
   TH1             *GetCumulative(Bool_t forward = kTRUE, const char* suffix = "_cumulative") const;

   TList           *GetListOfFunctions() const { return fFunctions; }

   virtual Int_t    GetNdivisions(Option_t *axis="X") const;
   virtual Color_t  GetAxisColor(Option_t *axis="X") const;
   virtual Color_t  GetLabelColor(Option_t *axis="X") const;
   virtual Style_t  GetLabelFont(Option_t *axis="X") const;
   virtual Float_t  GetLabelOffset(Option_t *axis="X") const;
   virtual Float_t  GetLabelSize(Option_t *axis="X") const;
   virtual Style_t  GetTitleFont(Option_t *axis="X") const;
   virtual Float_t  GetTitleOffset(Option_t *axis="X") const;
   virtual Float_t  GetTitleSize(Option_t *axis="X") const;
   virtual Float_t  GetTickLength(Option_t *axis="X") const;
   virtual Float_t  GetBarOffset() const {return Float_t(0.001*Float_t(fBarOffset));}
   virtual Float_t  GetBarWidth() const  {return Float_t(0.001*Float_t(fBarWidth));}
   virtual Int_t    GetContour(Double_t *levels=0);
   virtual Double_t GetContourLevel(Int_t level) const;
   virtual Double_t GetContourLevelPad(Int_t level) const;

   virtual Int_t    GetBin(Int_t binx, Int_t biny=0, Int_t binz=0) const;
   virtual void     GetBinXYZ(Int_t binglobal, Int_t &binx, Int_t &biny, Int_t &binz) const;
   virtual Double_t GetBinCenter(Int_t bin) const;
   virtual Double_t GetBinContent(Int_t bin) const;
   virtual Double_t GetBinContent(Int_t bin, Int_t) const { return GetBinContent(bin); }
   virtual Double_t GetBinContent(Int_t bin, Int_t, Int_t) const { return GetBinContent(bin); }
   virtual Double_t GetBinError(Int_t bin) const;
   virtual Double_t GetBinError(Int_t binx, Int_t biny) const { return GetBinError(GetBin(binx, biny)); } // for 2D histograms only
   virtual Double_t GetBinError(Int_t binx, Int_t biny, Int_t binz) const { return GetBinError(GetBin(binx, biny, binz)); } // for 3D histograms only
   virtual Double_t GetBinErrorLow(Int_t bin) const;
   virtual Double_t GetBinErrorUp(Int_t bin) const;
   virtual EBinErrorOpt  GetBinErrorOption() const { return fBinStatErrOpt; }
   virtual Double_t GetBinLowEdge(Int_t bin) const;
   virtual Double_t GetBinWidth(Int_t bin) const;
   virtual Double_t GetBinWithContent(Double_t c, Int_t &binx, Int_t firstx=0, Int_t lastx=0,Double_t maxdiff=0) const;
   virtual void     GetCenter(Double_t *center) const;
   static  Bool_t   GetDefaultSumw2();
   TDirectory      *GetDirectory() const {return fDirectory;}
   virtual Double_t GetEntries() const;
   virtual Double_t GetEffectiveEntries() const;
   virtual TF1     *GetFunction(const char *name) const;
   virtual Int_t    GetDimension() const { return fDimension; }
   virtual Double_t GetKurtosis(Int_t axis=1) const;
   virtual void     GetLowEdge(Double_t *edge) const;
   virtual Double_t GetMaximum(Double_t maxval=FLT_MAX) const;
   virtual Int_t    GetMaximumBin() const;
   virtual Int_t    GetMaximumBin(Int_t &locmax, Int_t &locmay, Int_t &locmaz) const;
   virtual Double_t GetMaximumStored() const {return fMaximum;}
   virtual Double_t GetMinimum(Double_t minval=-FLT_MAX) const;
   virtual Int_t    GetMinimumBin() const;
   virtual Int_t    GetMinimumBin(Int_t &locmix, Int_t &locmiy, Int_t &locmiz) const;
   virtual Double_t GetMinimumStored() const {return fMinimum;}
   virtual void     GetMinimumAndMaximum(Double_t& min, Double_t& max) const;
   virtual Double_t GetMean(Int_t axis=1) const;
   virtual Double_t GetMeanError(Int_t axis=1) const;
   virtual Int_t    GetNbinsX() const {return fXaxis.GetNbins();}
   virtual Int_t    GetNbinsY() const {return fYaxis.GetNbins();}
   virtual Int_t    GetNbinsZ() const {return fZaxis.GetNbins();}
   virtual Int_t    GetNcells() const {return fNcells; }
   virtual Double_t GetNormFactor() const {return fNormFactor;}
   virtual char    *GetObjectInfo(Int_t px, Int_t py) const;
   Option_t        *GetOption() const {return fOption.Data();}

   TVirtualHistPainter *GetPainter(Option_t *option="");

   virtual Int_t    GetQuantiles(Int_t nprobSum, Double_t *q, const Double_t *probSum=0);
   virtual Double_t GetRandom(TRandom * rng = nullptr) const;
   virtual void     GetStats(Double_t *stats) const;
   virtual Double_t GetStdDev(Int_t axis=1) const;
   virtual Double_t GetStdDevError(Int_t axis=1) const;
   virtual Double_t GetSumOfWeights() const;
   virtual TArrayD *GetSumw2() {return &fSumw2;}
   virtual const TArrayD *GetSumw2() const {return &fSumw2;}
   virtual Int_t    GetSumw2N() const {return fSumw2.fN;}
           Double_t GetRMS(Int_t axis=1) const { return GetStdDev(axis); }
           Double_t GetRMSError(Int_t axis=1) const { return GetStdDevError(axis); }

   virtual Double_t GetSkewness(Int_t axis=1) const;
           EStatOverflows GetStatOverflows() const {return fStatOverflows; }; ///< Get the behaviour adopted by the object about the statoverflows. See EStatOverflows for more information.
           TAxis*   GetXaxis()  { return &fXaxis; }
           TAxis*   GetYaxis()  { return &fYaxis; }
           TAxis*   GetZaxis()  { return &fZaxis; }
     const TAxis*   GetXaxis() const { return &fXaxis; }
     const TAxis*   GetYaxis() const { return &fYaxis; }
     const TAxis*   GetZaxis() const { return &fZaxis; }
   virtual Double_t Integral(Option_t *option="") const;
   virtual Double_t Integral(Int_t binx1, Int_t binx2, Option_t *option="") const;
   virtual Double_t IntegralAndError(Int_t binx1, Int_t binx2, Double_t & err, Option_t *option="") const;
   virtual Double_t Interpolate(Double_t x) const;
   virtual Double_t Interpolate(Double_t x, Double_t y) const;
   virtual Double_t Interpolate(Double_t x, Double_t y, Double_t z) const;
           Bool_t   IsBinOverflow(Int_t bin, Int_t axis = 0) const;
           Bool_t   IsBinUnderflow(Int_t bin, Int_t axis = 0) const;
   virtual Bool_t   IsHighlight() const { return TestBit(kIsHighlight); }
   virtual Double_t AndersonDarlingTest(const TH1 *h2, Option_t *option="") const;
   virtual Double_t AndersonDarlingTest(const TH1 *h2, Double_t &advalue) const;
   virtual Double_t KolmogorovTest(const TH1 *h2, Option_t *option="") const;
   virtual void     LabelsDeflate(Option_t *axis="X");
   virtual void     LabelsInflate(Option_t *axis="X");
   virtual void     LabelsOption(Option_t *option="h", Option_t *axis="X");
   virtual Long64_t Merge(TCollection *list) { return Merge(list,""); }
           Long64_t Merge(TCollection *list, Option_t * option);
   virtual Bool_t   Multiply(TF1 *f1, Double_t c1=1);
   virtual Bool_t   Multiply(const TH1 *h1);
   virtual Bool_t   Multiply(const TH1 *h1, const TH1 *h2, Double_t c1=1, Double_t c2=1, Option_t *option=""); // *MENU*
   virtual void     Paint(Option_t *option="");
   virtual void     Print(Option_t *option="") const;
   virtual void     PutStats(Double_t *stats);
   virtual TH1     *Rebin(Int_t ngroup=2, const char*newname="", const Double_t *xbins=0);  // *MENU*
   virtual TH1     *RebinX(Int_t ngroup=2, const char*newname="") { return Rebin(ngroup,newname, (Double_t*) 0); }
   virtual void     Rebuild(Option_t *option="");
   virtual void     RecursiveRemove(TObject *obj);
   virtual void     Reset(Option_t *option="");
   virtual void     ResetStats();
   virtual void     SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void     Scale(Double_t c1=1, Option_t *option="");
   virtual void     SetAxisColor(Color_t color=1, Option_t *axis="X");
   virtual void     SetAxisRange(Double_t xmin, Double_t xmax, Option_t *axis="X");
   virtual void     SetBarOffset(Float_t offset=0.25) {fBarOffset = Short_t(1000*offset);}
   virtual void     SetBarWidth(Float_t width=0.5) {fBarWidth = Short_t(1000*width);}
   virtual void     SetBinContent(Int_t bin, Double_t content);
   virtual void     SetBinContent(Int_t bin, Int_t, Double_t content) { SetBinContent(bin, content); }
   virtual void     SetBinContent(Int_t bin, Int_t, Int_t, Double_t content) { SetBinContent(bin, content); }
   virtual void     SetBinError(Int_t bin, Double_t error);
   virtual void     SetBinError(Int_t binx, Int_t biny, Double_t error);
   virtual void     SetBinError(Int_t binx, Int_t biny, Int_t binz, Double_t error);
   virtual void     SetBins(Int_t nx, Double_t xmin, Double_t xmax);
   virtual void     SetBins(Int_t nx, const Double_t *xBins);
   virtual void     SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax);
   virtual void     SetBins(Int_t nx, const Double_t *xBins, Int_t ny, const Double_t *yBins);
   virtual void     SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax,
                            Int_t nz, Double_t zmin, Double_t zmax);
   virtual void     SetBins(Int_t nx, const Double_t *xBins, Int_t ny, const Double_t * yBins, Int_t nz,
                            const Double_t *zBins);
   virtual void     SetBinsLength(Int_t = -1) { } //redefined in derived classes
   virtual void     SetBinErrorOption(EBinErrorOpt type) { fBinStatErrOpt = type; }
   virtual void     SetBuffer(Int_t buffersize, Option_t *option="");
   virtual UInt_t   SetCanExtend(UInt_t extendBitMask);
   virtual void     SetContent(const Double_t *content);
   virtual void     SetContour(Int_t nlevels, const Double_t *levels=0);
   virtual void     SetContourLevel(Int_t level, Double_t value);
   static  void     SetDefaultBufferSize(Int_t buffersize=1000);
   static  void     SetDefaultSumw2(Bool_t sumw2=kTRUE);
   virtual void     SetDirectory(TDirectory *dir);
   virtual void     SetEntries(Double_t n) {fEntries = n;};
   virtual void     SetError(const Double_t *error);
   virtual void     SetHighlight(Bool_t set = kTRUE); // *TOGGLE* *GETTER=IsHighlight
   virtual void     SetLabelColor(Color_t color=1, Option_t *axis="X");
   virtual void     SetLabelFont(Style_t font=62, Option_t *axis="X");
   virtual void     SetLabelOffset(Float_t offset=0.005, Option_t *axis="X");
   virtual void     SetLabelSize(Float_t size=0.02, Option_t *axis="X");

   /*
    * Set the minimum / maximum value for the Y axis (1-D histograms) or Z axis (2-D histograms)
    *   By default the maximum / minimum value used in drawing is the maximum / minimum value of the histogram
    * plus a margin of 10%. If these functions are called, the values are used without any extra margin.
    */
   virtual void     SetMaximum(Double_t maximum = -1111) { fMaximum = maximum; }; // *MENU*
   virtual void     SetMinimum(Double_t minimum = -1111) { fMinimum = minimum; }; // *MENU*

   virtual void     SetName(const char *name); // *MENU*
   virtual void     SetNameTitle(const char *name, const char *title);
   virtual void     SetNdivisions(Int_t n=510, Option_t *axis="X");
   virtual void     SetNormFactor(Double_t factor=1) {fNormFactor = factor;}
   virtual void     SetStats(Bool_t stats=kTRUE); // *MENU*
   virtual void     SetOption(Option_t *option=" ") {fOption = option;}
   virtual void     SetTickLength(Float_t length=0.02, Option_t *axis="X");
   virtual void     SetTitleFont(Style_t font=62, Option_t *axis="X");
   virtual void     SetTitleOffset(Float_t offset=1, Option_t *axis="X");
   virtual void     SetTitleSize(Float_t size=0.02, Option_t *axis="X");
           void     SetStatOverflows(EStatOverflows statOverflows) {fStatOverflows = statOverflows;}; ///< See GetStatOverflows for more information.
   virtual void     SetTitle(const char *title);  // *MENU*
   virtual void     SetXTitle(const char *title) {fXaxis.SetTitle(title);}
   virtual void     SetYTitle(const char *title) {fYaxis.SetTitle(title);}
   virtual void     SetZTitle(const char *title) {fZaxis.SetTitle(title);}
   virtual TH1     *ShowBackground(Int_t niter=20, Option_t *option="same"); // *MENU*
   virtual Int_t    ShowPeaks(Double_t sigma=2, Option_t *option="", Double_t threshold=0.05); // *MENU*
   virtual void     Smooth(Int_t ntimes=1, Option_t *option=""); // *MENU*
   static  void     SmoothArray(Int_t NN, Double_t *XX, Int_t ntimes=1);
   static  void     StatOverflows(Bool_t flag=kTRUE);
   virtual void     Sumw2(Bool_t flag = kTRUE);
   void             UseCurrentStyle();
   static  TH1     *TransformHisto(TVirtualFFT *fft, TH1* h_output,  Option_t *option);


   // TODO: Remove obsolete methods in v6-04
   virtual Double_t GetCellContent(Int_t binx, Int_t biny) const
                        { Obsolete("GetCellContent", "v6-00", "v6-04"); return GetBinContent(GetBin(binx, biny)); }
   virtual Double_t GetCellError(Int_t binx, Int_t biny) const
                        { Obsolete("GetCellError", "v6-00", "v6-04"); return GetBinError(binx, biny); }
   virtual void     RebinAxis(Double_t x, TAxis *axis)
                        { Obsolete("RebinAxis", "v6-00", "v6-04"); ExtendAxis(x, axis); }
   virtual void     SetCellContent(Int_t binx, Int_t biny, Double_t content)
                        { Obsolete("SetCellContent", "v6-00", "v6-04"); SetBinContent(GetBin(binx, biny), content); }
   virtual void     SetCellError(Int_t binx, Int_t biny, Double_t content)
                        { Obsolete("SetCellError", "v6-00", "v6-04"); SetBinError(binx, biny, content); }

   ClassDef(TH1,8)  //1-Dim histogram base class

protected:
   virtual Double_t RetrieveBinContent(Int_t bin) const;
   virtual void     UpdateBinContent(Int_t bin, Double_t content);
   virtual Double_t GetBinErrorSqUnchecked(Int_t bin) const { return fSumw2.fN ? fSumw2.fArray[bin] : RetrieveBinContent(bin); }
};

namespace cling {
  std::string printValue(TH1 *val);
}

//________________________________________________________________________

class TH1C : public TH1, public TArrayC {

public:
   TH1C();
   TH1C(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup);
   TH1C(const char *name,const char *title,Int_t nbinsx,const Float_t  *xbins);
   TH1C(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins);
   TH1C(const TH1C &h1c);
   TH1C& operator=(const TH1C &h1);
   virtual ~TH1C();

   virtual void     AddBinContent(Int_t bin);
   virtual void     AddBinContent(Int_t bin, Double_t w);
   virtual void     Copy(TObject &hnew) const;
   virtual void     Reset(Option_t *option="");
   virtual void     SetBinsLength(Int_t n=-1);

   ClassDef(TH1C,3)  //1-Dim histograms (one char per channel)

   friend  TH1C     operator*(Double_t c1, const TH1C &h1);
   friend  TH1C     operator*(const TH1C &h1, Double_t c1);
   friend  TH1C     operator+(const TH1C &h1, const TH1C &h2);
   friend  TH1C     operator-(const TH1C &h1, const TH1C &h2);
   friend  TH1C     operator*(const TH1C &h1, const TH1C &h2);
   friend  TH1C     operator/(const TH1C &h1, const TH1C &h2);

protected:
   virtual Double_t RetrieveBinContent(Int_t bin) const { return Double_t (fArray[bin]); }
   virtual void     UpdateBinContent(Int_t bin, Double_t content) { fArray[bin] = Char_t (content); }
};

TH1C operator*(Double_t c1, const TH1C &h1);
inline
TH1C operator*(const TH1C &h1, Double_t c1) {return operator*(c1,h1);}
TH1C operator+(const TH1C &h1, const TH1C &h2);
TH1C operator-(const TH1C &h1, const TH1C &h2);
TH1C operator*(const TH1C &h1, const TH1C &h2);
TH1C operator/(const TH1C &h1, const TH1C &h2);

//________________________________________________________________________

class TH1S : public TH1, public TArrayS {

public:
   TH1S();
   TH1S(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup);
   TH1S(const char *name,const char *title,Int_t nbinsx,const Float_t  *xbins);
   TH1S(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins);
   TH1S(const TH1S &h1s);
   TH1S& operator=(const TH1S &h1);
   virtual ~TH1S();

   virtual void     AddBinContent(Int_t bin);
   virtual void     AddBinContent(Int_t bin, Double_t w);
   virtual void     Copy(TObject &hnew) const;
   virtual void     Reset(Option_t *option="");
   virtual void     SetBinsLength(Int_t n=-1);

   ClassDef(TH1S,3)  //1-Dim histograms (one short per channel)

   friend  TH1S     operator*(Double_t c1, const TH1S &h1);
   friend  TH1S     operator*(const TH1S &h1, Double_t c1);
   friend  TH1S     operator+(const TH1S &h1, const TH1S &h2);
   friend  TH1S     operator-(const TH1S &h1, const TH1S &h2);
   friend  TH1S     operator*(const TH1S &h1, const TH1S &h2);
   friend  TH1S     operator/(const TH1S &h1, const TH1S &h2);

protected:
   virtual Double_t RetrieveBinContent(Int_t bin) const { return Double_t (fArray[bin]); }
   virtual void     UpdateBinContent(Int_t bin, Double_t content) { fArray[bin] = Short_t (content); }
};

TH1S operator*(Double_t c1, const TH1S &h1);
inline
TH1S operator*(const TH1S &h1, Double_t c1) {return operator*(c1,h1);}
TH1S operator+(const TH1S &h1, const TH1S &h2);
TH1S operator-(const TH1S &h1, const TH1S &h2);
TH1S operator*(const TH1S &h1, const TH1S &h2);
TH1S operator/(const TH1S &h1, const TH1S &h2);

//________________________________________________________________________

class TH1I: public TH1, public TArrayI {

public:
   TH1I();
   TH1I(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup);
   TH1I(const char *name,const char *title,Int_t nbinsx,const Float_t  *xbins);
   TH1I(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins);
   TH1I(const TH1I &h1i);
   TH1I& operator=(const TH1I &h1);
   virtual ~TH1I();

   virtual void     AddBinContent(Int_t bin);
   virtual void     AddBinContent(Int_t bin, Double_t w);
   virtual void     Copy(TObject &hnew) const;
   virtual void     Reset(Option_t *option="");
   virtual void     SetBinsLength(Int_t n=-1);

   ClassDef(TH1I,3)  //1-Dim histograms (one 32 bits integer per channel)

   friend  TH1I     operator*(Double_t c1, const TH1I &h1);
   friend  TH1I     operator*(const TH1I &h1, Double_t c1);
   friend  TH1I     operator+(const TH1I &h1, const TH1I &h2);
   friend  TH1I     operator-(const TH1I &h1, const TH1I &h2);
   friend  TH1I     operator*(const TH1I &h1, const TH1I &h2);
   friend  TH1I     operator/(const TH1I &h1, const TH1I &h2);

protected:
   virtual Double_t RetrieveBinContent(Int_t bin) const { return Double_t (fArray[bin]); }
   virtual void     UpdateBinContent(Int_t bin, Double_t content) { fArray[bin] = Int_t (content); }
};

TH1I operator*(Double_t c1, const TH1I &h1);
inline
TH1I operator*(const TH1I &h1, Double_t c1) {return operator*(c1,h1);}
TH1I operator+(const TH1I &h1, const TH1I &h2);
TH1I operator-(const TH1I &h1, const TH1I &h2);
TH1I operator*(const TH1I &h1, const TH1I &h2);
TH1I operator/(const TH1I &h1, const TH1I &h2);

//________________________________________________________________________

class TH1F : public TH1, public TArrayF {

public:
   TH1F();
   TH1F(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup);
   TH1F(const char *name,const char *title,Int_t nbinsx,const Float_t  *xbins);
   TH1F(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins);
   explicit TH1F(const TVectorF &v);
   TH1F(const TH1F &h1f);
   TH1F& operator=(const TH1F &h1);
   virtual ~TH1F();

   virtual void     AddBinContent(Int_t bin) {++fArray[bin];}
   virtual void     AddBinContent(Int_t bin, Double_t w)
                                 {fArray[bin] += Float_t (w);}
   virtual void     Copy(TObject &hnew) const;
   virtual void     Reset(Option_t *option="");
   virtual void     SetBinsLength(Int_t n=-1);

   ClassDef(TH1F,3)  //1-Dim histograms (one float per channel)

   friend  TH1F     operator*(Double_t c1, const TH1F &h1);
   friend  TH1F     operator*(const TH1F &h1, Double_t c1);
   friend  TH1F     operator+(const TH1F &h1, const TH1F &h2);
   friend  TH1F     operator-(const TH1F &h1, const TH1F &h2);
   friend  TH1F     operator*(const TH1F &h1, const TH1F &h2);
   friend  TH1F     operator/(const TH1F &h1, const TH1F &h2);

protected:
   virtual Double_t RetrieveBinContent(Int_t bin) const { return Double_t (fArray[bin]); }
   virtual void     UpdateBinContent(Int_t bin, Double_t content) { fArray[bin] = Float_t (content); }
};

TH1F operator*(Double_t c1, const TH1F &h1);
inline
TH1F operator*(const TH1F &h1, Double_t c1) {return operator*(c1,h1);}
TH1F operator+(const TH1F &h1, const TH1F &h2);
TH1F operator-(const TH1F &h1, const TH1F &h2);
TH1F operator*(const TH1F &h1, const TH1F &h2);
TH1F operator/(const TH1F &h1, const TH1F &h2);

//________________________________________________________________________

class TH1D : public TH1, public TArrayD {

public:
   TH1D();
   TH1D(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup);
   TH1D(const char *name,const char *title,Int_t nbinsx,const Float_t  *xbins);
   TH1D(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins);
   explicit TH1D(const TVectorD &v);
   TH1D(const TH1D &h1d);
   TH1D& operator=(const TH1D &h1);
   virtual ~TH1D();

   virtual void     AddBinContent(Int_t bin) {++fArray[bin];}
   virtual void     AddBinContent(Int_t bin, Double_t w)
                                 {fArray[bin] += Double_t (w);}
   virtual void     Copy(TObject &hnew) const;
   virtual void     Reset(Option_t *option="");
   virtual void     SetBinsLength(Int_t n=-1);

   ClassDef(TH1D,3)  //1-Dim histograms (one double per channel)

   friend  TH1D     operator*(Double_t c1, const TH1D &h1);
   friend  TH1D     operator*(const TH1D &h1, Double_t c1);
   friend  TH1D     operator+(const TH1D &h1, const TH1D &h2);
   friend  TH1D     operator-(const TH1D &h1, const TH1D &h2);
   friend  TH1D     operator*(const TH1D &h1, const TH1D &h2);
   friend  TH1D     operator/(const TH1D &h1, const TH1D &h2);

protected:
   virtual Double_t RetrieveBinContent(Int_t bin) const { return fArray[bin]; }
   virtual void     UpdateBinContent(Int_t bin, Double_t content) { fArray[bin] = content; }
};

TH1D operator*(Double_t c1, const TH1D &h1);
inline
TH1D operator*(const TH1D &h1, Double_t c1) {return operator*(c1,h1);}
TH1D operator+(const TH1D &h1, const TH1D &h2);
TH1D operator-(const TH1D &h1, const TH1D &h2);
TH1D operator*(const TH1D &h1, const TH1D &h2);
TH1D operator/(const TH1D &h1, const TH1D &h2);

   extern TH1 *R__H(Int_t hid);
   extern TH1 *R__H(const char *hname);

#endif
