// @(#)root/hist:$Id$
// Author: Rene Brun   26/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cctype>
#include <sstream>
#include <cmath>
#include <iostream>

#include "TROOT.h"
#include "TBuffer.h"
#include "TEnv.h"
#include "TClass.h"
#include "TMath.h"
#include "THashList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF2.h"
#include "TF3.h"
#include "TPluginManager.h"
#include "TVirtualPad.h"
#include "TRandom.h"
#include "TVirtualFitter.h"
#include "THLimitsFinder.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TVectorF.h"
#include "TVectorD.h"
#include "TBrowser.h"
#include "TError.h"
#include "TVirtualHistPainter.h"
#include "TVirtualFFT.h"
#include "TVirtualPaveStats.h"

#include "HFitInterface.h"
#include "Fit/DataRange.h"
#include "Fit/BinData.h"
#include "Math/GoFTest.h"
#include "Math/MinimizerOptions.h"
#include "Math/QuantFuncMathCore.h"

#include "TH1Merger.h"

/** \addtogroup Hist
@{
\class TH1C
\brief 1-D histogram with a byte per channel (see TH1 documentation)
\class TH1S
\brief 1-D histogram with a short per channel (see TH1 documentation)
\class TH1I
\brief 1-D histogram with an int per channel (see TH1 documentation)}
\class TH1F
\brief 1-D histogram with a float per channel (see TH1 documentation)}
\class TH1D
\brief 1-D histogram with a double per channel (see TH1 documentation)}
@}
*/

/** \class TH1
TH1 is the base class of all histogramm classes in ROOT. It provides the common interface for
operations such as binning, filling, drawing, which will be detailed below.

-# [Creating histograms](\ref creating-histograms)
-# [Binning](\ref binning)
-# [Filling histograms](\ref filling-histograms)
-# [Operations on histograms](\ref operations-on-histograms)
-# [Drawing histograms](\ref drawing-histograms)
-# [Saving/reading histograms to/from a ROOT file](\ref saving-histograms)
-# [Miscellaneous operations](\ref misc)

ROOT supports the following histogram types:

  - 1-D histograms:
      - TH1C : histograms with one byte per channel.   Maximum bin content = 127
      - TH1S : histograms with one short per channel.  Maximum bin content = 32767
      - TH1I : histograms with one int per channel.    Maximum bin content = 2147483647
      - TH1F : histograms with one float per channel.  Maximum precision 7 digits
      - TH1D : histograms with one double per channel. Maximum precision 14 digits
  - 2-D histograms:
      - TH2C : histograms with one byte per channel.   Maximum bin content = 127
      - TH2S : histograms with one short per channel.  Maximum bin content = 32767
      - TH2I : histograms with one int per channel.    Maximum bin content = 2147483647
      - TH2F : histograms with one float per channel.  Maximum precision 7 digits
      - TH2D : histograms with one double per channel. Maximum precision 14 digits
  - 3-D histograms:
      - TH3C : histograms with one byte per channel.   Maximum bin content = 127
      - TH3S : histograms with one short per channel.  Maximum bin content = 32767
      - TH3I : histograms with one int per channel.    Maximum bin content = 2147483647
      - TH3F : histograms with one float per channel.  Maximum precision 7 digits
      - TH3D : histograms with one double per channel. Maximum precision 14 digits
  - Profile histograms: See classes  TProfile, TProfile2D and TProfile3D.
      Profile histograms are used to display the mean value of Y and its standard deviation
      for each bin in X. Profile histograms are in many cases an elegant
      replacement of two-dimensional histograms : the inter-relation of two
      measured quantities X and Y can always be visualized by a two-dimensional
      histogram or scatter-plot; If Y is an unknown (but single-valued)
      approximate function of X, this function is displayed by a profile
      histogram with much better precision than by a scatter-plot.

The inheritance hierarchy looks as follows:
~~~
                                TH1
                                 ^
                                 |
                                 |
                                 |
                +----------------+-------+------+------+-----+-----+
                |                |       |      |      |     |     |
                |                |      TH1C   TH1S   TH1I  TH1F  TH1D
                |                |                                 |
                |                |                                 |
                |               TH2                             TProfile
                |                |
                |                |
                |                +-------+------+------+-----+-----+
                |                        |      |      |     |     |
                |                       TH2C   TH2S   TH2I  TH2F  TH2D
                |                                                  |
               TH3                                                 |
                |                                               TProfile2D
                |
                +-------+------+------+------+------+
                        |      |      |      |      |
                       TH3C   TH3S   TH3I   TH3F   TH3D
                                                    |
                                                    |
                                                 TProfile3D

      The TH*C classes also inherit from the array class TArrayC.
      The TH*S classes also inherit from the array class TArrayS.
      The TH*I classes also inherit from the array class TArrayI.
      The TH*F classes also inherit from the array class TArrayF.
      The TH*D classes also inherit from the array class TArrayD.
~~~


\anchor creating-histograms
## Creating histograms

Histograms are created by invoking one of the constructors, e.g.
~~~ {.cpp}
       TH1F *h1 = new TH1F("h1", "h1 title", 100, 0, 4.4);
       TH2F *h2 = new TH2F("h2", "h2 title", 40, 0, 4, 30, -3, 3);
~~~
Histograms may also be created by:

  -  calling the Clone function, see below
  -  making a projection from a 2-D or 3-D histogram, see below
  -  reading an histogram from a file

 When a histogram is created, a reference to it is automatically added
 to the list of in-memory objects for the current file or directory.
 This default behaviour can be changed by:
~~~ {.cpp}
       h->SetDirectory(0);          // for the current histogram h
       TH1::AddDirectory(kFALSE);   // sets a global switch disabling the referencing
~~~
 When the histogram is deleted, the reference to it is removed from
 the list of objects in memory.
 When a file is closed, all histograms in memory associated with this file
 are automatically deleted.

### Labelling axes

 Axis titles can be specified in the title argument of the constructor.
 They must be separated by ";":
~~~ {.cpp}
        TH1F* h=new TH1F("h", "Histogram title;X Axis;Y Axis", 100, 0, 1);
~~~
 The histogram title and the axis titles can be any TLatex string, and
 are persisted if a histogram is written to a file.

 Any title can be omitted:
~~~ {.cpp}
        TH1F* h=new TH1F("h", "Histogram title;;Y Axis", 100, 0, 1);
        TH1F* h=new TH1F("h", ";;Y Axis", 100, 0, 1);
~~~
 The method SetTitle() has the same syntax:
~~~ {.cpp}
        h->SetTitle("Histogram title;Another X title Axis");
~~~
Alternatively, the title of each axis can be set directly:
~~~ {.cpp}
       h->GetXaxis()->SetTitle("X axis title");
       h->GetYaxis()->SetTitle("Y axis title");
~~~
For bin labels see \ref binning.

\anchor binning
## Binning

#### Fix or variable bin size

 All histogram types support either fix or variable bin sizes.
 2-D histograms may have fix size bins along X and variable size bins
 along Y or vice-versa. The functions to fill, manipulate, draw or access
 histograms are identical in both cases.

 Each histogram always contains 3 axis objects of type TAxis: fXaxis, fYaxis and fZaxis.
 To access the axis parameters, use:
~~~ {.cpp}
        TAxis *xaxis = h->GetXaxis(); etc.
        Double_t binCenter = xaxis->GetBinCenter(bin), etc.
~~~
 See class TAxis for a description of all the access functions.
 The axis range is always stored internally in double precision.

#### Convention for numbering bins

 For all histogram types: nbins, xlow, xup
~~~ {.cpp}
        bin = 0;       underflow bin
        bin = 1;       first bin with low-edge xlow INCLUDED
        bin = nbins;   last bin with upper-edge xup EXCLUDED
        bin = nbins+1; overflow bin
~~~
 In case of 2-D or 3-D histograms, a "global bin" number is defined.
 For example, assuming a 3-D histogram with (binx, biny, binz), the function
~~~ {.cpp}
        Int_t gbin = h->GetBin(binx, biny, binz);
~~~
 returns a global/linearized gbin number. This global gbin is useful
 to access the bin content/error information independently of the dimension.
 Note that to access the information other than bin content and errors
 one should use the TAxis object directly with e.g.:
~~~ {.cpp}
         Double_t xcenter = h3->GetZaxis()->GetBinCenter(27);
~~~
 returns the center along z of bin number 27 (not the global bin)
 in the 3-D histogram h3.

#### Alphanumeric Bin Labels

 By default, an histogram axis is drawn with its numeric bin labels.
 One can specify alphanumeric labels instead with:

  - call TAxis::SetBinLabel(bin, label);
    This can always be done before or after filling.
    When the histogram is drawn, bin labels will be automatically drawn.
    See examples labels1.C and labels2.C
  - call to a Fill function with one of the arguments being a string, e.g.
~~~ {.cpp}
           hist1->Fill(somename, weight);
           hist2->Fill(x, somename, weight);
           hist2->Fill(somename, y, weight);
           hist2->Fill(somenamex, somenamey, weight);
~~~
    See examples hlabels1.C and hlabels2.C
  - via TTree::Draw. see for example cernstaff.C
~~~ {.cpp}
           tree.Draw("Nation::Division");
~~~
    where "Nation" and "Division" are two branches of a Tree.

When using the options 2 or 3 above, the labels are automatically
 added to the list (THashList) of labels for a given axis.
 By default, an axis is drawn with the order of bins corresponding
 to the filling sequence. It is possible to reorder the axis

  - alphabetically
  - by increasing or decreasing values

 The reordering can be triggered via the TAxis context menu by selecting
 the menu item "LabelsOption" or by calling directly
 TH1::LabelsOption(option, axis) where

  - axis may be "X", "Y" or "Z"
  - option may be:
    - "a" sort by alphabetic order
    - ">" sort by decreasing values
    - "<" sort by increasing values
    - "h" draw labels horizontal
    - "v" draw labels vertical
    - "u" draw labels up (end of label right adjusted)
    - "d" draw labels down (start of label left adjusted)

 When using the option 2 above, new labels are added by doubling the current
 number of bins in case one label does not exist yet.
 When the Filling is terminated, it is possible to trim the number
 of bins to match the number of active labels by calling
~~~ {.cpp}
           TH1::LabelsDeflate(axis) with axis = "X", "Y" or "Z"
~~~
 This operation is automatic when using TTree::Draw.
 Once bin labels have been created, they become persistent if the histogram
 is written to a file or when generating the C++ code via SavePrimitive.

#### Histograms with automatic bins

 When an histogram is created with an axis lower limit greater or equal
 to its upper limit, the SetBuffer is automatically called with an
 argument fBufferSize equal to fgBufferSize (default value=1000).
 fgBufferSize may be reset via the static function TH1::SetDefaultBufferSize.
 The axis limits will be automatically computed when the buffer will
 be full or when the function BufferEmpty is called.

#### Rebinning

 At any time, an histogram can be rebinned via TH1::Rebin. This function
 returns a new histogram with the rebinned contents.
 If bin errors were stored, they are recomputed during the rebinning.



\anchor filling-histograms
## Filling histograms

 An histogram is typically filled with statements like:
~~~ {.cpp}
       h1->Fill(x);
       h1->Fill(x, w); //fill with weight
       h2->Fill(x, y)
       h2->Fill(x, y, w)
       h3->Fill(x, y, z)
       h3->Fill(x, y, z, w)
~~~
 or via one of the Fill functions accepting names described above.
 The Fill functions compute the bin number corresponding to the given
 x, y or z argument and increment this bin by the given weight.
 The Fill functions return the bin number for 1-D histograms or global
 bin number for 2-D and 3-D histograms.
 If TH1::Sumw2 has been called before filling, the sum of squares of
 weights is also stored.
 One can also increment directly a bin number via TH1::AddBinContent
 or replace the existing content via TH1::SetBinContent.
 To access the bin content of a given bin, do:
~~~ {.cpp}
       Double_t binContent = h->GetBinContent(bin);
~~~

 By default, the bin number is computed using the current axis ranges.
 If the automatic binning option has been set via
~~~ {.cpp}
       h->SetCanExtend(TH1::kAllAxes);
~~~
 then, the Fill Function will automatically extend the axis range to
 accomodate the new value specified in the Fill argument. The method
 used is to double the bin size until the new value fits in the range,
 merging bins two by two. This automatic binning options is extensively
 used by the TTree::Draw function when histogramming Tree variables
 with an unknown range.
 This automatic binning option is supported for 1-D, 2-D and 3-D histograms.

 During filling, some statistics parameters are incremented to compute
 the mean value and Root Mean Square with the maximum precision.

 In case of histograms of type TH1C, TH1S, TH2C, TH2S, TH3C, TH3S
 a check is made that the bin contents do not exceed the maximum positive
 capacity (127 or 32767). Histograms of all types may have positive
 or/and negative bin contents.

### Associated errors
 By default, for each bin, the sum of weights is computed at fill time.
 One can also call TH1::Sumw2 to force the storage and computation
 of the sum of the square of weights per bin.
 If Sumw2 has been called, the error per bin is computed as the
 sqrt(sum of squares of weights), otherwise the error is set equal
 to the sqrt(bin content).
 To return the error for a given bin number, do:
~~~ {.cpp}
        Double_t error = h->GetBinError(bin);
~~~

### Associated functions
 One or more object (typically a TF1*) can be added to the list
 of functions (fFunctions) associated to each histogram.
 When TH1::Fit is invoked, the fitted function is added to this list.
 Given an histogram h, one can retrieve an associated function
 with:
~~~ {.cpp}
        TF1 *myfunc = h->GetFunction("myfunc");
~~~


\anchor operations-on-histograms
## Operations on histograms

 Many types of operations are supported on histograms or between histograms

  -  Addition of an histogram to the current histogram.
  -  Additions of two histograms with coefficients and storage into the current
     histogram.
  -  Multiplications and Divisions are supported in the same way as additions.
  -  The Add, Divide and Multiply functions also exist to add, divide or multiply
     an histogram by a function.

 If an histogram has associated error bars (TH1::Sumw2 has been called),
 the resulting error bars are also computed assuming independent histograms.
 In case of divisions, Binomial errors are also supported.
 One can mark a histogram to be an "average" histogram by setting its bit kIsAverage via
   myhist.SetBit(TH1::kIsAverage);
 When adding (see TH1::Add) average histograms, the histograms are averaged and not summed.

#### Fitting histograms

 Histograms (1-D, 2-D, 3-D and Profiles) can be fitted with a user
 specified function via TH1::Fit. When an histogram is fitted, the
 resulting function with its parameters is added to the list of functions
 of this histogram. If the histogram is made persistent, the list of
 associated functions is also persistent. Given a pointer (see above)
 to an associated function myfunc, one can retrieve the function/fit
 parameters with calls such as:
~~~ {.cpp}
       Double_t chi2 = myfunc->GetChisquare();
       Double_t par0 = myfunc->GetParameter(0); value of 1st parameter
       Double_t err0 = myfunc->GetParError(0);  error on first parameter
~~~

#### Projections of histograms

 One can:

  -  make a 1-D projection of a 2-D histogram or Profile
     see functions TH2::ProjectionX,Y, TH2::ProfileX,Y, TProfile::ProjectionX
  -  make a 1-D, 2-D or profile out of a 3-D histogram
     see functions TH3::ProjectionZ, TH3::Project3D.

 One can fit these projections via:
~~~ {.cpp}
      TH2::FitSlicesX,Y, TH3::FitSlicesZ.
~~~

#### Random Numbers and histograms

 TH1::FillRandom can be used to randomly fill an histogram using
 the contents of an existing TF1 function or another
 TH1 histogram (for all dimensions).
 For example the following two statements create and fill an histogram
 10000 times with a default gaussian distribution of mean 0 and sigma 1:
~~~ {.cpp}
       TH1F h1("h1", "histo from a gaussian", 100, -3, 3);
       h1.FillRandom("gaus", 10000);
~~~
 TH1::GetRandom can be used to return a random number distributed
 according the contents of an histogram.

#### Making a copy of an histogram
 Like for any other ROOT object derived from TObject, one can use
 the Clone() function. This makes an identical copy of the original
 histogram including all associated errors and functions, e.g.:
~~~ {.cpp}
       TH1F *hnew = (TH1F*)h->Clone("hnew");
~~~

#### Normalizing histograms

 One can scale an histogram such that the bins integral is equal to
 the normalization parameter via TH1::Scale(Double_t norm), where norm
 is the desired normalization divided by the integral of the histogram.


\anchor drawing-histograms
## Drawing histograms

 Histograms are drawn via the THistPainter class. Each histogram has
 a pointer to its own painter (to be usable in a multithreaded program).
 Many drawing options are supported.
 See THistPainter::Paint() for more details.

 The same histogram can be drawn with different options in different pads.
 When an histogram drawn in a pad is deleted, the histogram is
 automatically removed from the pad or pads where it was drawn.
 If an histogram is drawn in a pad, then filled again, the new status
 of the histogram will be automatically shown in the pad next time
 the pad is updated. One does not need to redraw the histogram.
 To draw the current version of an histogram in a pad, one can use
~~~ {.cpp}
        h->DrawCopy();
~~~
 This makes a clone (see Clone below) of the histogram. Once the clone
 is drawn, the original histogram may be modified or deleted without
 affecting the aspect of the clone.

 One can use TH1::SetMaximum() and TH1::SetMinimum() to force a particular
 value for the maximum or the minimum scale on the plot. (For 1-D
 histograms this means the y-axis, while for 2-D histograms these
 functions affect the z-axis).

 TH1::UseCurrentStyle() can be used to change all histogram graphics
 attributes to correspond to the current selected style.
 This function must be called for each histogram.
 In case one reads and draws many histograms from a file, one can force
 the histograms to inherit automatically the current graphics style
 by calling before gROOT->ForceStyle().

#### Setting Drawing histogram contour levels (2-D hists only)

 By default contours are automatically generated at equidistant
 intervals. A default value of 20 levels is used. This can be modified
 via TH1::SetContour() or TH1::SetContourLevel().
 the contours level info is used by the drawing options "cont", "surf",
 and "lego".

#### Setting histogram graphics attributes

 The histogram classes inherit from the attribute classes:
 TAttLine, TAttFill, and TAttMarker.
 See the member functions of these classes for the list of options.

#### Customising how axes are drawn

 Use the functions of TAxis, such as
~~~ {.cpp}
 histogram.GetXaxis()->SetTicks("+");
 histogram.GetYaxis()->SetRangeUser(1., 5.);
~~~

\anchor saving-histograms
## Saving/Reading histograms in ROOT files

 The following statements create a ROOT file and store an histogram
 on the file. Because TH1 derives from TNamed, the key identifier on
 the file is the histogram name:
~~~ {.cpp}
        TFile f("histos.root", "new");
        TH1F h1("hgaus", "histo from a gaussian", 100, -3, 3);
        h1.FillRandom("gaus", 10000);
        h1->Write();
~~~
 To read this histogram in another Root session, do:
~~~ {.cpp}
        TFile f("histos.root");
        TH1F *h = (TH1F*)f.Get("hgaus");
~~~
 One can save all histograms in memory to the file by:
~~~ {.cpp}
        file->Write();
~~~


\anchor misc
## Miscellaneous operations

~~~ {.cpp}
        TH1::KolmogorovTest(): statistical test of compatibility in shape
                             between two histograms
        TH1::Smooth() smooths the bin contents of a 1-d histogram
        TH1::Integral() returns the integral of bin contents in a given bin range
        TH1::GetMean(int axis) returns the mean value along axis
        TH1::GetStdDev(int axis)  returns the sigma distribution along axis
        TH1::GetEntries() returns the number of entries
        TH1::Reset() resets the bin contents and errors of an histogram
~~~
 IMPORTANT NOTE: The returned values for GetMean and GetStdDev depend on how the
 histogram statistics are calculated. By default, if no range has been set, the
 returned values are the (unbinned) ones calculated at fill time. If a range has been
 set, however, the values are calculated using the bins in range; THIS IS TRUE EVEN
 IF THE RANGE INCLUDES ALL BINS--use TAxis::SetRange(0, 0) to unset the range.
 To ensure that the returned values are always those of the binned data stored in the
 histogram, call TH1::ResetStats. See TH1::GetStats.
*/

TF1 *gF1=0;  //left for back compatibility (use TVirtualFitter::GetUserFunc instead)

Int_t  TH1::fgBufferSize   = 1000;
Bool_t TH1::fgAddDirectory = kTRUE;
Bool_t TH1::fgDefaultSumw2 = kFALSE;
Bool_t TH1::fgStatOverflows= kFALSE;

extern void H1InitGaus();
extern void H1InitExpo();
extern void H1InitPolynom();
extern void H1LeastSquareFit(Int_t n, Int_t m, Double_t *a);
extern void H1LeastSquareLinearFit(Int_t ndata, Double_t &a0, Double_t &a1, Int_t &ifail);
extern void H1LeastSquareSeqnd(Int_t n, Double_t *a, Int_t idim, Int_t &ifail, Int_t k, Double_t *b);

// Internal exceptions for the CheckConsistency method
class DifferentDimension: public std::exception {};
class DifferentNumberOfBins: public std::exception {};
class DifferentAxisLimits: public std::exception {};
class DifferentBinLimits: public std::exception {};
class DifferentLabels: public std::exception {};

ClassImp(TH1);

////////////////////////////////////////////////////////////////////////////////
/// Histogram default constructor.

TH1::TH1(): TNamed(), TAttLine(), TAttFill(), TAttMarker()
{
   fDirectory     = 0;
   fFunctions     = new TList;
   fNcells        = 0;
   fIntegral      = 0;
   fPainter       = 0;
   fEntries       = 0;
   fNormFactor    = 0;
   fTsumw         = fTsumw2=fTsumwx=fTsumwx2=0;
   fMaximum       = -1111;
   fMinimum       = -1111;
   fBufferSize    = 0;
   fBuffer        = 0;
   fBinStatErrOpt = kNormal;
   fStatOverflows = EStatOverflows::kNeutral;
   fXaxis.SetName("xaxis");
   fYaxis.SetName("yaxis");
   fZaxis.SetName("zaxis");
   fXaxis.SetParent(this);
   fYaxis.SetParent(this);
   fZaxis.SetParent(this);
   UseCurrentStyle();
}

////////////////////////////////////////////////////////////////////////////////
/// Histogram default destructor.

TH1::~TH1()
{
   if (!TestBit(kNotDeleted)) {
      return;
   }
   delete[] fIntegral;
   fIntegral = 0;
   delete[] fBuffer;
   fBuffer = 0;
   if (fFunctions) {
      R__WRITE_LOCKGUARD(ROOT::gCoreMutex);

      fFunctions->SetBit(kInvalidObject);
      TObject* obj = 0;
      //special logic to support the case where the same object is
      //added multiple times in fFunctions.
      //This case happens when the same object is added with different
      //drawing modes
      //In the loop below we must be careful with objects (eg TCutG) that may
      // have been added to the list of functions of several histograms
      //and may have been already deleted.
      while ((obj  = fFunctions->First())) {
         while(fFunctions->Remove(obj)) { }
         if (!obj->TestBit(kNotDeleted)) {
            break;
         }
         delete obj;
         obj = 0;
      }
      delete fFunctions;
      fFunctions = 0;
   }
   if (fDirectory) {
      fDirectory->Remove(this);
      fDirectory = 0;
   }
   delete fPainter;
   fPainter = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for fix bin size histograms.
/// Creates the main histogram structure.
///
/// \param[in] name name of histogram (avoid blanks)
/// \param[in] title histogram title.
///            If title is of the form `stringt;stringx;stringy;stringz`,
///            the histogram title is set to `stringt`,
///            the x axis title to `stringx`, the y axis title to `stringy`, etc.
/// \param[in] nbins number of bins
/// \param[in] xlow low edge of first bin
/// \param[in] xup upper edge of last bin (not included in last bin)
///
/// When an histogram is created, it is automatically added to the list
/// of special objects in the current directory.
/// To find the pointer to this histogram in the current directory
/// by its name, do:
/// ~~~ {.cpp}
///  TH1F *h1 = (TH1F*)gDirectory->FindObject(name);
/// ~~~

TH1::TH1(const char *name,const char *title,Int_t nbins,Double_t xlow,Double_t xup)
    :TNamed(name,title), TAttLine(), TAttFill(), TAttMarker()
{
   Build();
   if (nbins <= 0) {Warning("TH1","nbins is <=0 - set to nbins = 1"); nbins = 1; }
   fXaxis.Set(nbins,xlow,xup);
   fNcells = fXaxis.GetNbins()+2;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size histograms.
/// Creates the main histogram structure.
///
/// \param[in] name name of histogram (avoid blanks)
/// \param[in] title histogram title.
///            If title is of the form `stringt;stringx;stringy;stringz`
///            the histogram title is set to `stringt`,
///            the x axis title to `stringx`, the y axis title to `stringy`, etc.
/// \param[in] nbins number of bins
/// \param[in] xbins array of low-edges for each bin.
///            This is an array of size nbins+1

TH1::TH1(const char *name,const char *title,Int_t nbins,const Float_t *xbins)
    :TNamed(name,title), TAttLine(), TAttFill(), TAttMarker()
{
   Build();
   if (nbins <= 0) {Warning("TH1","nbins is <=0 - set to nbins = 1"); nbins = 1; }
   if (xbins) fXaxis.Set(nbins,xbins);
   else       fXaxis.Set(nbins,0,1);
   fNcells = fXaxis.GetNbins()+2;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for variable bin size histograms.
///
/// \param[in] name name of histogram (avoid blanks)
/// \param[in] title histogram title.
///        If title is of the form `stringt;stringx;stringy;stringz`
///        the histogram title is set to `stringt`,
///        the x axis title to `stringx`, the y axis title to `stringy`, etc.
/// \param[in] nbins number of bins
/// \param[in] xbins array of low-edges for each bin.
///        This is an array of size nbins+1

TH1::TH1(const char *name,const char *title,Int_t nbins,const Double_t *xbins)
    :TNamed(name,title), TAttLine(), TAttFill(), TAttMarker()
{
   Build();
   if (nbins <= 0) {Warning("TH1","nbins is <=0 - set to nbins = 1"); nbins = 1; }
   if (xbins) fXaxis.Set(nbins,xbins);
   else       fXaxis.Set(nbins,0,1);
   fNcells = fXaxis.GetNbins()+2;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.
/// The list of functions is not copied. (Use Clone if needed)

TH1::TH1(const TH1 &h) : TNamed(), TAttLine(), TAttFill(), TAttMarker()
{
   ((TH1&)h).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Static function: cannot be inlined on Windows/NT.

Bool_t TH1::AddDirectoryStatus()
{
   return fgAddDirectory;
}

////////////////////////////////////////////////////////////////////////////////
/// Browse the Histogram object.

void TH1::Browse(TBrowser *b)
{
   Draw(b ? b->GetDrawOption() : "");
   gPad->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates histogram basic data structure.

void TH1::Build()
{
   fDirectory     = 0;
   fPainter       = 0;
   fIntegral      = 0;
   fEntries       = 0;
   fNormFactor    = 0;
   fTsumw         = fTsumw2=fTsumwx=fTsumwx2=0;
   fMaximum       = -1111;
   fMinimum       = -1111;
   fBufferSize    = 0;
   fBuffer        = 0;
   fBinStatErrOpt = kNormal;
   fStatOverflows = EStatOverflows::kNeutral;
   fXaxis.SetName("xaxis");
   fYaxis.SetName("yaxis");
   fZaxis.SetName("zaxis");
   fYaxis.Set(1,0.,1.);
   fZaxis.Set(1,0.,1.);
   fXaxis.SetParent(this);
   fYaxis.SetParent(this);
   fZaxis.SetParent(this);

   SetTitle(fTitle.Data());

   fFunctions = new TList;

   UseCurrentStyle();

   if (TH1::AddDirectoryStatus()) {
      fDirectory = gDirectory;
      if (fDirectory) {
         fFunctions->UseRWLock();
         fDirectory->Append(this,kTRUE);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: `this = this + c1*f1`
/// if errors are defined (see TH1::Sumw2), errors are also recalculated.
///
/// By default, the function is computed at the centre of the bin.
/// if option "I" is specified (1-d histogram only), the integral of the
/// function in each bin is used instead of the value of the function at
/// the centre of the bin.
///
/// Only bins inside the function range are recomputed.
///
/// IMPORTANT NOTE: If you intend to use the errors of this histogram later
/// you should call Sumw2 before making this operation.
/// This is particularly important if you fit the histogram after TH1::Add
///
/// The function return kFALSE if the Add operation failed

Bool_t TH1::Add(TF1 *f1, Double_t c1, Option_t *option)
{
   if (!f1) {
      Error("Add","Attempt to add a non-existing function");
      return kFALSE;
   }

   TString opt = option;
   opt.ToLower();
   Bool_t integral = kFALSE;
   if (opt.Contains("i") && fDimension == 1) integral = kTRUE;

   Int_t ncellsx = GetNbinsX() + 2; // cells = normal bins + underflow bin + overflow bin
   Int_t ncellsy = GetNbinsY() + 2;
   Int_t ncellsz = GetNbinsZ() + 2;
   if (fDimension < 2) ncellsy = 1;
   if (fDimension < 3) ncellsz = 1;

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   //   - Add statistics
   Double_t s1[10];
   for (Int_t i = 0; i < 10; ++i) s1[i] = 0;
   PutStats(s1);
   SetMinimum();
   SetMaximum();

   //   - Loop on bins (including underflows/overflows)
   Int_t bin, binx, biny, binz;
   Double_t cu=0;
   Double_t xx[3];
   Double_t *params = 0;
   f1->InitArgs(xx,params);
   for (binz = 0; binz < ncellsz; ++binz) {
      xx[2] = fZaxis.GetBinCenter(binz);
      for (biny = 0; biny < ncellsy; ++biny) {
         xx[1] = fYaxis.GetBinCenter(biny);
         for (binx = 0; binx < ncellsx; ++binx) {
            xx[0] = fXaxis.GetBinCenter(binx);
            if (!f1->IsInside(xx)) continue;
            TF1::RejectPoint(kFALSE);
            bin = binx + ncellsx * (biny + ncellsy * binz);
            if (integral) {
               cu = c1*f1->Integral(fXaxis.GetBinLowEdge(binx), fXaxis.GetBinUpEdge(binx), 0.) / fXaxis.GetBinWidth(binx);
            } else {
               cu  = c1*f1->EvalPar(xx);
            }
            if (TF1::RejectedPoint()) continue;
            AddBinContent(bin,cu);
         }
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: `this = this + c1*h1`
/// If errors are defined (see TH1::Sumw2), errors are also recalculated.
///
/// Note that if h1 has Sumw2 set, Sumw2 is automatically called for this
/// if not already set.
///
/// Note also that adding histogram with labels is not supported, histogram will be
/// added merging them by bin number independently of the labels.
/// For adding histogram with labels one should use TH1::Merge
///
/// SPECIAL CASE (Average/Efficiency histograms)
/// For histograms representing averages or efficiencies, one should compute the average
/// of the two histograms and not the sum. One can mark a histogram to be an average
/// histogram by setting its bit kIsAverage with
/// myhist.SetBit(TH1::kIsAverage);
/// Note that the two histograms must have their kIsAverage bit set
///
/// IMPORTANT NOTE1: If you intend to use the errors of this histogram later
/// you should call Sumw2 before making this operation.
/// This is particularly important if you fit the histogram after TH1::Add
///
/// IMPORTANT NOTE2: if h1 has a normalisation factor, the normalisation factor
/// is used , ie  this = this + c1*factor*h1
/// Use the other TH1::Add function if you do not want this feature
///
/// IMPORTANT NOTE3: You should be careful about the statistics of the
/// returned histogram, whose statistics may be binned or unbinned,
/// depending on whether c1 is negative, whether TAxis::kAxisRange is true,
/// and whether TH1::ResetStats has been called on either this or h1.
/// See TH1::GetStats.
///
/// The function return kFALSE if the Add operation failed

Bool_t TH1::Add(const TH1 *h1, Double_t c1)
{
   if (!h1) {
      Error("Add","Attempt to add a non-existing histogram");
      return kFALSE;
   }

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   bool useMerge = (c1 == 1. &&  !this->TestBit(kIsAverage) && !h1->TestBit(kIsAverage) );
   try {
      CheckConsistency(this,h1);
      useMerge = kFALSE;
   } catch(DifferentNumberOfBins&) {
      if (useMerge)
         Info("Add","Attempt to add histograms with different number of bins - trying to use TH1::Merge");
      else {
         Error("Add","Attempt to add histograms with different number of bins : nbins h1 = %d , nbins h2 =  %d",GetNbinsX(), h1->GetNbinsX());
         return kFALSE;
      }
   } catch(DifferentAxisLimits&) {
      if (useMerge)
         Info("Add","Attempt to add histograms with different axis limits - trying to use TH1::Merge");
      else
         Warning("Add","Attempt to add histograms with different axis limits");
   } catch(DifferentBinLimits&) {
      if (useMerge)
         Info("Add","Attempt to add histograms with different bin limits - trying to use TH1::Merge");
      else
         Warning("Add","Attempt to add histograms with different bin limits");
   } catch(DifferentLabels&) {
      // in case of different labels -
      if (useMerge)
         Info("Add","Attempt to add histograms with different labels - trying to use TH1::Merge");
      else
         Info("Warning","Attempt to add histograms with different labels");
   }

   if (useMerge) {
      TList l;
      l.Add(const_cast<TH1*>(h1));
      auto iret = Merge(&l);
      return (iret >= 0);
   }

   //    Create Sumw2 if h1 has Sumw2 set
   if (fSumw2.fN == 0 && h1->GetSumw2N() != 0) Sumw2();

   //   - Add statistics
   Double_t entries = TMath::Abs( GetEntries() + c1 * h1->GetEntries() );

   // statistics can be preserved only in case of positive coefficients
   // otherwise with negative c1 (histogram subtraction) one risks to get negative variances
   Bool_t resetStats = (c1 < 0);
   Double_t s1[kNstat] = {0};
   Double_t s2[kNstat] = {0};
   if (!resetStats) {
      // need to initialize to zero s1 and s2 since
      // GetStats fills only used elements depending on dimension and type
      GetStats(s1);
      h1->GetStats(s2);
   }

   SetMinimum();
   SetMaximum();

   //   - Loop on bins (including underflows/overflows)
   Double_t factor = 1;
   if (h1->GetNormFactor() != 0) factor = h1->GetNormFactor()/h1->GetSumOfWeights();;
   Double_t c1sq = c1 * c1;
   Double_t factsq = factor * factor;

   for (Int_t bin = 0; bin < fNcells; ++bin) {
      //special case where histograms have the kIsAverage bit set
      if (this->TestBit(kIsAverage) && h1->TestBit(kIsAverage)) {
         Double_t y1 = h1->RetrieveBinContent(bin);
         Double_t y2 = this->RetrieveBinContent(bin);
         Double_t e1sq = h1->GetBinErrorSqUnchecked(bin);
         Double_t e2sq = this->GetBinErrorSqUnchecked(bin);
         Double_t w1 = 1., w2 = 1.;

         // consider all special cases  when bin errors are zero
         // see http://root-forum.cern.ch/viewtopic.php?f=3&t=13299
         if (e1sq) w1 = 1. / e1sq;
         else if (h1->fSumw2.fN) {
            w1 = 1.E200; // use an arbitrary huge value
            if (y1 == 0) {
               // use an estimated error from the global histogram scale
               double sf = (s2[0] != 0) ? s2[1]/s2[0] : 1;
               w1 = 1./(sf*sf);
            }
         }
         if (e2sq) w2 = 1. / e2sq;
         else if (fSumw2.fN) {
            w2 = 1.E200; // use an arbitrary huge value
            if (y2 == 0) {
               // use an estimated error from the global histogram scale
               double sf = (s1[0] != 0) ? s1[1]/s1[0] : 1;
               w2 = 1./(sf*sf);
            }
         }

         double y =  (w1*y1 + w2*y2)/(w1 + w2);
         UpdateBinContent(bin, y);
         if (fSumw2.fN) {
            double err2 =  1./(w1 + w2);
            if (err2 < 1.E-200) err2 = 0;  // to remove arbitrary value when e1=0 AND e2=0
            fSumw2.fArray[bin] = err2;
         }
      } else { // normal case of addition between histograms
         AddBinContent(bin, c1 * factor * h1->RetrieveBinContent(bin));
         if (fSumw2.fN) fSumw2.fArray[bin] += c1sq * factsq * h1->GetBinErrorSqUnchecked(bin);
      }
   }

   // update statistics (do here to avoid changes by SetBinContent)
   if (resetStats)  {
      // statistics need to be reset in case coefficient are negative
      ResetStats();
   }
   else {
      for (Int_t i=0;i<kNstat;i++) {
         if (i == 1) s1[i] += c1*c1*s2[i];
         else        s1[i] += c1*s2[i];
      }
      PutStats(s1);
      SetEntries(entries);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace contents of this histogram by the addition of h1 and h2.
///
/// `this = c1*h1 + c2*h2`
/// if errors are defined (see TH1::Sumw2), errors are also recalculated
///
/// Note that if h1 or h2 have Sumw2 set, Sumw2 is automatically called for this
/// if not already set.
///
/// Note also that adding histogram with labels is not supported, histogram will be
/// added merging them by bin number independently of the labels.
/// For adding histogram ith labels one should use TH1::Merge
///
/// SPECIAL CASE (Average/Efficiency histograms)
/// For histograms representing averages or efficiencies, one should compute the average
/// of the two histograms and not the sum. One can mark a histogram to be an average
/// histogram by setting its bit kIsAverage with
/// myhist.SetBit(TH1::kIsAverage);
/// Note that the two histograms must have their kIsAverage bit set
///
/// IMPORTANT NOTE: If you intend to use the errors of this histogram later
/// you should call Sumw2 before making this operation.
/// This is particularly important if you fit the histogram after TH1::Add
///
/// IMPORTANT NOTE2: You should be careful about the statistics of the
/// returned histogram, whose statistics may be binned or unbinned,
/// depending on whether c1 is negative, whether TAxis::kAxisRange is true,
/// and whether TH1::ResetStats has been called on either this or h1.
/// See TH1::GetStats.
///
/// ANOTHER SPECIAL CASE : h1 = h2 and c2 < 0
/// do a scaling   this = c1 * h1 / (bin Volume)
///
/// The function returns kFALSE if the Add operation failed

Bool_t TH1::Add(const TH1 *h1, const TH1 *h2, Double_t c1, Double_t c2)
{

   if (!h1 || !h2) {
      Error("Add","Attempt to add a non-existing histogram");
      return kFALSE;
   }

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   Bool_t normWidth = kFALSE;
   if (h1 == h2 && c2 < 0) {c2 = 0; normWidth = kTRUE;}

   if (h1 != h2) {
      bool useMerge = (c1 == 1. && c2 == 1. &&  !this->TestBit(kIsAverage) && !h1->TestBit(kIsAverage) );

      try {
         CheckConsistency(h1,h2);
         CheckConsistency(this,h1);
         useMerge = kFALSE;
      } catch(DifferentNumberOfBins&) {
         if (useMerge)
            Info("Add","Attempt to add histograms with different number of bins - trying to use TH1::Merge");
         else {
            Error("Add","Attempt to add histograms with different number of bins : nbins h1 = %d , nbins h2 =  %d",GetNbinsX(), h1->GetNbinsX());
            return kFALSE;
         }
      } catch(DifferentAxisLimits&) {
         if (useMerge)
            Info("Add","Attempt to add histograms with different axis limits - trying to use TH1::Merge");
         else
            Warning("Add","Attempt to add histograms with different axis limits");
      } catch(DifferentBinLimits&) {
         if (useMerge)
            Info("Add","Attempt to add histograms with different bin limits - trying to use TH1::Merge");
         else
            Warning("Add","Attempt to add histograms with different bin limits");
      } catch(DifferentLabels&) {
         // in case of different labels -
         if (useMerge)
            Info("Add","Attempt to add histograms with different labels - trying to use TH1::Merge");
         else
            Info("Warning","Attempt to add histograms with different labels");
      }

      if (useMerge) {
         TList l;
         // why TList takes non-const pointers ????
         l.Add(const_cast<TH1*>(h1));
         l.Add(const_cast<TH1*>(h2));
         Reset("ICE");
         auto iret = Merge(&l);
         return (iret >= 0);
      }
   }

   //    Create Sumw2 if h1 or h2 have Sumw2 set
   if (fSumw2.fN == 0 && (h1->GetSumw2N() != 0 || h2->GetSumw2N() != 0)) Sumw2();

   //   - Add statistics
   Double_t nEntries = TMath::Abs( c1*h1->GetEntries() + c2*h2->GetEntries() );

   // TODO remove
   // statistics can be preserved only in case of positive coefficients
   // otherwise with negative c1 (histogram subtraction) one risks to get negative variances
   // also in case of scaling with the width we cannot preserve the statistics
   Double_t s1[kNstat] = {0};
   Double_t s2[kNstat] = {0};
   Double_t s3[kNstat];


   Bool_t resetStats = (c1*c2 < 0) || normWidth;
   if (!resetStats) {
      // need to initialize to zero s1 and s2 since
      // GetStats fills only used elements depending on dimension and type
      h1->GetStats(s1);
      h2->GetStats(s2);
      for (Int_t i=0;i<kNstat;i++) {
         if (i == 1) s3[i] = c1*c1*s1[i] + c2*c2*s2[i];
         //else        s3[i] = TMath::Abs(c1)*s1[i] + TMath::Abs(c2)*s2[i];
         else        s3[i] = c1*s1[i] + c2*s2[i];
      }
   }

   SetMinimum();
   SetMaximum();

   if (normWidth) { // DEPRECATED CASE: belongs to fitting / drawing modules

      Int_t nbinsx = GetNbinsX() + 2; // normal bins + underflow, overflow
      Int_t nbinsy = GetNbinsY() + 2;
      Int_t nbinsz = GetNbinsZ() + 2;

      if (fDimension < 2) nbinsy = 1;
      if (fDimension < 3) nbinsz = 1;

      Int_t bin, binx, biny, binz;
      for (binz = 0; binz < nbinsz; ++binz) {
         Double_t wz = h1->GetZaxis()->GetBinWidth(binz);
         for (biny = 0; biny < nbinsy; ++biny) {
            Double_t wy = h1->GetYaxis()->GetBinWidth(biny);
            for (binx = 0; binx < nbinsx; ++binx) {
               Double_t wx = h1->GetXaxis()->GetBinWidth(binx);
               bin = GetBin(binx, biny, binz);
               Double_t w = wx*wy*wz;
               UpdateBinContent(bin, c1 * h1->RetrieveBinContent(bin) / w);
               if (fSumw2.fN) {
                  Double_t e1 = h1->GetBinError(bin)/w;
                  fSumw2.fArray[bin] = c1*c1*e1*e1;
               }
            }
         }
      }
   } else if (h1->TestBit(kIsAverage) && h2->TestBit(kIsAverage)) {
      for (Int_t i = 0; i < fNcells; ++i) { // loop on cells (bins including underflow / overflow)
         // special case where histograms have the kIsAverage bit set
         Double_t y1 = h1->RetrieveBinContent(i);
         Double_t y2 = h2->RetrieveBinContent(i);
         Double_t e1sq = h1->GetBinErrorSqUnchecked(i);
         Double_t e2sq = h2->GetBinErrorSqUnchecked(i);
         Double_t w1 = 1., w2 = 1.;

         // consider all special cases  when bin errors are zero
         // see http://root-forum.cern.ch/viewtopic.php?f=3&t=13299
         if (e1sq) w1 = 1./ e1sq;
         else if (h1->fSumw2.fN) {
            w1 = 1.E200; // use an arbitrary huge value
            if (y1 == 0 ) { // use an estimated error from the global histogram scale
               double sf = (s1[0] != 0) ? s1[1]/s1[0] : 1;
               w1 = 1./(sf*sf);
            }
         }
         if (e2sq) w2 = 1./ e2sq;
         else if (h2->fSumw2.fN) {
            w2 = 1.E200; // use an arbitrary huge value
            if (y2 == 0) { // use an estimated error from the global histogram scale
               double sf = (s2[0] != 0) ? s2[1]/s2[0] : 1;
               w2 = 1./(sf*sf);
            }
         }

         double y =  (w1*y1 + w2*y2)/(w1 + w2);
         UpdateBinContent(i, y);
         if (fSumw2.fN) {
            double err2 =  1./(w1 + w2);
            if (err2 < 1.E-200) err2 = 0;  // to remove arbitrary value when e1=0 AND e2=0
            fSumw2.fArray[i] = err2;
         }
      }
   } else { // case of simple histogram addition
      Double_t c1sq = c1 * c1;
      Double_t c2sq = c2 * c2;
      for (Int_t i = 0; i < fNcells; ++i) { // Loop on cells (bins including underflows/overflows)
         UpdateBinContent(i, c1 * h1->RetrieveBinContent(i) + c2 * h2->RetrieveBinContent(i));
         if (fSumw2.fN) {
            fSumw2.fArray[i] = c1sq * h1->GetBinErrorSqUnchecked(i) + c2sq * h2->GetBinErrorSqUnchecked(i);
         }
      }
   }

   if (resetStats)  {
      // statistics need to be reset in case coefficient are negative
      ResetStats();
   }
   else {
      // update statistics (do here to avoid changes by SetBinContent)  FIXME remove???
      PutStats(s3);
      SetEntries(nEntries);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by 1.

void TH1::AddBinContent(Int_t)
{
   AbstractMethod("AddBinContent");
}

////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by a weight w.

void TH1::AddBinContent(Int_t, Double_t)
{
   AbstractMethod("AddBinContent");
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the flag controlling the automatic add of histograms in memory
///
/// By default (fAddDirectory = kTRUE), histograms are automatically added
/// to the list of objects in memory.
/// Note that one histogram can be removed from its support directory
/// by calling h->SetDirectory(0) or h->SetDirectory(dir) to add it
/// to the list of objects in the directory dir.
///
/// NOTE that this is a static function. To call it, use;
/// TH1::AddDirectory

void TH1::AddDirectory(Bool_t add)
{
   fgAddDirectory = add;
}

////////////////////////////////////////////////////////////////////////////////
/// Auxilliary function to get the power of 2 next (larger) or previous (smaller)
/// a given x
///
///    next = kTRUE  : next larger
///    next = kFALSE : previous smaller
///
/// Used by the autobin power of 2 algorithm

inline Double_t TH1::AutoP2GetPower2(Double_t x, Bool_t next)
{
   Int_t nn;
   Double_t f2 = std::frexp(x, &nn);
   return ((next && x > 0.) || (!next && x <= 0.)) ? std::ldexp(std::copysign(1., f2), nn)
                                                   : std::ldexp(std::copysign(1., f2), --nn);
}

////////////////////////////////////////////////////////////////////////////////
/// Auxilliary function to get the next power of 2 integer value larger then n
///
/// Used by the autobin power of 2 algorithm

inline Int_t TH1::AutoP2GetBins(Int_t n)
{
   Int_t nn;
   Double_t f2 = std::frexp(n, &nn);
   if (TMath::Abs(f2 - .5) > 0.001)
      return (Int_t)std::ldexp(1., nn);
   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Buffer-based estimate of the histogram range using the power of 2 algorithm.
///
/// Used by the autobin power of 2 algorithm.
///
/// Works on arguments (min and max from fBuffer) and internal inputs: fXmin,
/// fXmax, NBinsX (from fXaxis), ...
/// Result save internally in fXaxis.
///
/// Overloaded by TH2 and TH3.
///
/// Return -1 if internal inputs are incosistent, 0 otherwise.
///

Int_t TH1::AutoP2FindLimits(Double_t xmi, Double_t xma)
{
   // We need meaningful raw limits
   if (xmi >= xma)
      return -1;

   THLimitsFinder::GetLimitsFinder()->FindGoodLimits(this, xmi, xma);
   Double_t xhmi = fXaxis.GetXmin();
   Double_t xhma = fXaxis.GetXmax();

   // Now adjust
   if (TMath::Abs(xhma) > TMath::Abs(xhmi)) {
      // Start from the upper limit
      xhma = TH1::AutoP2GetPower2(xhma);
      xhmi = xhma - TH1::AutoP2GetPower2(xhma - xhmi);
   } else {
      // Start from the lower limit
      xhmi = TH1::AutoP2GetPower2(xhmi, kFALSE);
      xhma = xhmi + TH1::AutoP2GetPower2(xhma - xhmi);
   }

   // Round the bins to the next power of 2; take into account the possible inflation
   // of the range
   Double_t rr = (xhma - xhmi) / (xma - xmi);
   Int_t nb = TH1::AutoP2GetBins((Int_t)(rr * GetNbinsX()));

   // Adjust using the same bin width and offsets
   Double_t bw = (xhma - xhmi) / nb;
   // Bins to left free on each side
   Double_t autoside = gEnv->GetValue("Hist.Binning.Auto.Side", 0.05);
   Int_t nbside = (Int_t)(nb * autoside);

   // Side up
   Int_t nbup = (xhma - xma) / bw;
   if (nbup % 2 != 0)
      nbup++; // Must be even
   if (nbup != nbside) {
      // Accounts also for both case: larger or smaller
      xhma -= bw * (nbup - nbside);
      nb -= (nbup - nbside);
   }

   // Side low
   Int_t nblw = (xmi - xhmi) / bw;
   if (nblw % 2 != 0)
      nblw++; // Must be even
   if (nblw != nbside) {
      // Accounts also for both case: larger or smaller
      xhmi += bw * (nblw - nbside);
      nb -= (nblw - nbside);
   }

   // Set everything and project
   SetBins(nb, xhmi, xhma);

   // Done
   return 0;
}

/// Fill histogram with all entries in the buffer.
///
///  - action = -1 histogram is reset and refilled from the buffer (called by THistPainter::Paint)
///  - action =  0 histogram is reset and filled from the buffer. When the histogram is filled from the
///                buffer the value fBuffer[0] is set to a negative number (= - number of entries)
///                When calling with action == 0 the histogram is NOT refilled when fBuffer[0] is < 0
///                While when calling with action = -1 the histogram is reset and ALWAYS refilled independently if
///                the histogram was filled before. This is needed when drawing the histogram
///  - action =  1 histogram is filled and buffer is deleted
///                The buffer is automatically deleted when filling the histogram and the entries is
///                larger than the buffer size

Int_t TH1::BufferEmpty(Int_t action)
{
   // do we need to compute the bin size?
   if (!fBuffer) return 0;
   Int_t nbentries = (Int_t)fBuffer[0];

   // nbentries correspond to the number of entries of histogram

   if (nbentries == 0) {
      // if action is 1 we delete the buffer
      // this will avoid infinite recursion
      if (action > 0) {
         delete [] fBuffer;
         fBuffer = 0;
         fBufferSize = 0;
      }
      return 0;
   }
   if (nbentries < 0 && action == 0) return 0;    // case histogram has been already filled from the buffer

   Double_t *buffer = fBuffer;
   if (nbentries < 0) {
      nbentries  = -nbentries;
      //  a reset might call BufferEmpty() giving an infinite recursion
      // Protect it by setting fBuffer = 0
      fBuffer=0;
       //do not reset the list of functions
      Reset("ICES");
      fBuffer = buffer;
   }
   if (CanExtendAllAxes() || (fXaxis.GetXmax() <= fXaxis.GetXmin())) {
      //find min, max of entries in buffer
      Double_t xmin = fBuffer[2];
      Double_t xmax = xmin;
      for (Int_t i=1;i<nbentries;i++) {
         Double_t x = fBuffer[2*i+2];
         if (x < xmin) xmin = x;
         if (x > xmax) xmax = x;
      }
      if (fXaxis.GetXmax() <= fXaxis.GetXmin()) {
         Int_t rc = -1;
         if (TestBit(TH1::kAutoBinPTwo)) {
            if ((rc = AutoP2FindLimits(xmin, xmax)) < 0)
               Warning("BufferEmpty",
                       "incosistency found by power-of-2 autobin algorithm: fallback to standard method");
         }
         if (rc < 0)
            THLimitsFinder::GetLimitsFinder()->FindGoodLimits(this, xmin, xmax);
      } else {
         fBuffer = 0;
         Int_t keep = fBufferSize; fBufferSize = 0;
         if (xmin <  fXaxis.GetXmin()) ExtendAxis(xmin, &fXaxis);
         if (xmax >= fXaxis.GetXmax()) ExtendAxis(xmax, &fXaxis);
         fBuffer = buffer;
         fBufferSize = keep;
      }
   }

   // call DoFillN which will not put entries in the buffer as FillN does
   // set fBuffer to zero to avoid re-emptying the buffer from functions called
   // by DoFillN (e.g Sumw2)
   buffer = fBuffer; fBuffer = 0;
   DoFillN(nbentries,&buffer[2],&buffer[1],2);
   fBuffer = buffer;

   // if action == 1 - delete the buffer
   if (action > 0) {
      delete [] fBuffer;
      fBuffer = 0;
      fBufferSize = 0;
   } else {
      // if number of entries is consistent with buffer - set it negative to avoid
      // refilling the histogram every time BufferEmpty(0) is called
      // In case it is not consistent, by setting fBuffer[0]=0 is like resetting the buffer
      // (it will not be used anymore the next time BufferEmpty is called)
      if (nbentries == (Int_t)fEntries)
         fBuffer[0] = -nbentries;
      else
         fBuffer[0] = 0;
   }
   return nbentries;
}

////////////////////////////////////////////////////////////////////////////////
/// accumulate arguments in buffer. When buffer is full, empty the buffer
///
///  - `fBuffer[0]` = number of entries in buffer
///  - `fBuffer[1]` = w of first entry
///  - `fBuffer[2]` = x of first entry

Int_t TH1::BufferFill(Double_t x, Double_t w)
{
   if (!fBuffer) return -2;
   Int_t nbentries = (Int_t)fBuffer[0];


   if (nbentries < 0) {
      // reset nbentries to a positive value so next time BufferEmpty()  is called
      // the histogram will be refilled
      nbentries  = -nbentries;
      fBuffer[0] =  nbentries;
      if (fEntries > 0) {
         // set fBuffer to zero to avoid calling BufferEmpty in Reset
         Double_t *buffer = fBuffer; fBuffer=0;
         Reset("ICES");  // do not reset list of functions
         fBuffer = buffer;
      }
   }
   if (2*nbentries+2 >= fBufferSize) {
      BufferEmpty(1);
      if (!fBuffer)
         // to avoid infinite recursion Fill->BufferFill->Fill
         return Fill(x,w);
      // this cannot happen
      R__ASSERT(0);
   }
   fBuffer[2*nbentries+1] = w;
   fBuffer[2*nbentries+2] = x;
   fBuffer[0] += 1;
   return -2;
}

////////////////////////////////////////////////////////////////////////////////
/// Check bin limits.

bool TH1::CheckBinLimits(const TAxis* a1, const TAxis * a2)
{
   const TArrayD * h1Array = a1->GetXbins();
   const TArrayD * h2Array = a2->GetXbins();
   Int_t fN = h1Array->fN;
   if ( fN != 0 ) {
      if ( h2Array->fN != fN ) {
         throw DifferentBinLimits();
         return false;
      }
      else {
         for ( int i = 0; i < fN; ++i ) {
            // for i==fN (nbin+1) a->GetBinWidth() returns last bin width
            // we do not need to exclude that case
            double binWidth = a1->GetBinWidth(i);
            if ( ! TMath::AreEqualAbs( h1Array->GetAt(i), h2Array->GetAt(i), binWidth*1E-10 ) ) {
               throw DifferentBinLimits();
               return false;
            }
         }
      }
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Check that axis have same labels.

bool TH1::CheckBinLabels(const TAxis* a1, const TAxis * a2)
{
   THashList *l1 = a1->GetLabels();
   THashList *l2 = a2->GetLabels();

   if (!l1 && !l2 )
      return true;
   if (!l1 ||  !l2 ) {
      throw DifferentLabels();
      return false;
   }
   // check now labels sizes  are the same
   if (l1->GetSize() != l2->GetSize() ) {
      throw DifferentLabels();
      return false;
   }
   for (int i = 1; i <= a1->GetNbins(); ++i) {
      TString label1 = a1->GetBinLabel(i);
      TString label2 = a2->GetBinLabel(i);
      if (label1 != label2) {
         throw DifferentLabels();
         return false;
      }
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Check that the axis limits of the histograms are the same.
/// If a first and last bin is passed the axis is compared between the given range

bool TH1::CheckAxisLimits(const TAxis *a1, const TAxis *a2 )
{
   double firstBin = a1->GetBinWidth(1);
   double lastBin = a1->GetBinWidth( a1->GetNbins() );
   if ( ! TMath::AreEqualAbs(a1->GetXmin(), a2->GetXmin(), firstBin* 1.E-10) ||
        ! TMath::AreEqualAbs(a1->GetXmax(), a2->GetXmax(), lastBin*1.E-10) ) {
      throw DifferentAxisLimits();
      return false;
   }
   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Check that the axis are the same

bool TH1::CheckEqualAxes(const TAxis *a1, const TAxis *a2 )
{
   if (a1->GetNbins() != a2->GetNbins() ) {
      //throw DifferentNumberOfBins();
      ::Info("CheckEqualAxes","Axes have different number of bins : nbin1 = %d nbin2 = %d",a1->GetNbins(),a2->GetNbins() );
      return false;
   }
   try {
      CheckAxisLimits(a1,a2);
   } catch (DifferentAxisLimits&) {
      ::Info("CheckEqualAxes","Axes have different limits");
      return false;
   }
   try {
      CheckBinLimits(a1,a2);
   } catch (DifferentBinLimits&) {
      ::Info("CheckEqualAxes","Axes have different bin limits");
      return false;
   }

   // check labels
   try {
      CheckBinLabels(a1,a2);
   } catch (DifferentLabels&) {
      ::Info("CheckEqualAxes","Axes have different labels");
      return false;
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Check that two sub axis are the same.
/// The limits are defined by first bin and last bin
/// N.B. no check is done in this case for variable bins

bool TH1::CheckConsistentSubAxes(const TAxis *a1, Int_t firstBin1, Int_t lastBin1, const TAxis * a2, Int_t firstBin2, Int_t lastBin2 )
{
   // By default is assumed that no bins are given for the second axis
   Int_t nbins1   = lastBin1-firstBin1 + 1;
   Double_t xmin1 = a1->GetBinLowEdge(firstBin1);
   Double_t xmax1 = a1->GetBinUpEdge(lastBin1);

   Int_t nbins2 = a2->GetNbins();
   Double_t xmin2 = a2->GetXmin();
   Double_t xmax2 = a2->GetXmax();

   if (firstBin2 <  lastBin2) {
      // in this case assume no bins are given for the second axis
      nbins2   = lastBin1-firstBin1 + 1;
      xmin2 = a1->GetBinLowEdge(firstBin1);
      xmax2 = a1->GetBinUpEdge(lastBin1);
   }

   if (nbins1 != nbins2 ) {
      ::Info("CheckConsistentSubAxes","Axes have different number of bins");
      return false;
   }

   Double_t firstBin = a1->GetBinWidth(firstBin1);
   Double_t lastBin = a1->GetBinWidth(lastBin1);
   if ( ! TMath::AreEqualAbs(xmin1,xmin2,1.E-10 * firstBin) ||
        ! TMath::AreEqualAbs(xmax1,xmax2,1.E-10 * lastBin) ) {
      ::Info("CheckConsistentSubAxes","Axes have different limits");
      return false;
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Check histogram compatibility.

bool TH1::CheckConsistency(const TH1* h1, const TH1* h2)
{
   if (h1 == h2) return true;

   if (h1->GetDimension() != h2->GetDimension() ) {
      throw DifferentDimension();
      return false;
   }
   Int_t dim = h1->GetDimension();

   // returns kTRUE if number of bins and bin limits are identical
   Int_t nbinsx = h1->GetNbinsX();
   Int_t nbinsy = h1->GetNbinsY();
   Int_t nbinsz = h1->GetNbinsZ();

   // Check whether the histograms have the same number of bins.
   if (nbinsx != h2->GetNbinsX() ||
       (dim > 1 && nbinsy != h2->GetNbinsY())  ||
       (dim > 2 && nbinsz != h2->GetNbinsZ()) ) {
      throw DifferentNumberOfBins();
      return false;
   }

   bool ret = true;

   // check axis limits
   ret &= CheckAxisLimits(h1->GetXaxis(), h2->GetXaxis());
   if (dim > 1) ret &= CheckAxisLimits(h1->GetYaxis(), h2->GetYaxis());
   if (dim > 2) ret &= CheckAxisLimits(h1->GetZaxis(), h2->GetZaxis());

   // check bin limits
   ret &= CheckBinLimits(h1->GetXaxis(), h2->GetXaxis());
   if (dim > 1) ret &= CheckBinLimits(h1->GetYaxis(), h2->GetYaxis());
   if (dim > 2) ret &= CheckBinLimits(h1->GetZaxis(), h2->GetZaxis());

   // check labels if histograms are both not empty
   if ( !h1->IsEmpty() && !h2->IsEmpty() ) {
      ret &= CheckBinLabels(h1->GetXaxis(), h2->GetXaxis());
      if (dim > 1) ret &= CheckBinLabels(h1->GetYaxis(), h2->GetYaxis());
      if (dim > 2) ret &= CheckBinLabels(h1->GetZaxis(), h2->GetZaxis());
   }

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// \f$ \chi^{2} \f$ test for comparing weighted and unweighted histograms
///
/// Function: Returns p-value. Other return values are specified by the 3rd parameter
///
/// \param[in] h2 the second histogram
/// \param[in] option
///   - "UU" = experiment experiment comparison (unweighted-unweighted)
///   - "UW" = experiment MC comparison (unweighted-weighted). Note that
///      the first histogram should be unweighted
///   - "WW" = MC MC comparison (weighted-weighted)
///   - "NORM" = to be used when one or both of the histograms is scaled
///              but the histogram originally was unweighted
///   - by default underflows and overflows are not included:
///      * "OF" = overflows included
///      * "UF" = underflows included
///   - "P" = print chi2, ndf, p_value, igood
///   - "CHI2" = returns chi2 instead of p-value
///   - "CHI2/NDF" = returns \f$ \chi^{2} \f$/ndf
/// \param[in] res not empty - computes normalized residuals and returns them in this array
///
/// The current implementation is based on the papers \f$ \chi^{2} \f$ test for comparison
/// of weighted and unweighted histograms" in Proceedings of PHYSTAT05 and
/// "Comparison weighted and unweighted histograms", arXiv:physics/0605123
/// by N.Gagunashvili. This function has been implemented by Daniel Haertl in August 2006.
///
/// #### Introduction:
///
/// A frequently used technique in data analysis is the comparison of
/// histograms. First suggested by Pearson [1] the \f$ \chi^{2} \f$  test of
/// homogeneity is used widely for comparing usual (unweighted) histograms.
/// This paper describes the implementation modified \f$ \chi^{2} \f$ tests
/// for comparison of weighted and unweighted  histograms and two weighted
/// histograms [2] as well as usual Pearson's \f$ \chi^{2} \f$ test for
/// comparison two usual (unweighted) histograms.
///
/// #### Overview:
///
/// Comparison of two histograms expect hypotheses that two histograms
/// represent identical distributions. To make a decision p-value should
/// be calculated. The hypotheses of identity is rejected if the p-value is
/// lower then some significance level. Traditionally significance levels
/// 0.1, 0.05 and 0.01 are used. The comparison procedure should include an
/// analysis of the residuals which is often helpful in identifying the
/// bins of histograms responsible for a significant overall \f$ \chi^{2} \f$ value.
/// Residuals are the difference between bin contents and expected bin
/// contents. Most convenient for analysis are the normalized residuals. If
/// hypotheses of identity are valid then normalized residuals are
/// approximately independent and identically distributed random variables
/// having N(0,1) distribution. Analysis of residuals expect test of above
/// mentioned properties of residuals. Notice that indirectly the analysis
/// of residuals increase the power of \f$ \chi^{2} \f$ test.
///
/// #### Methods of comparison:
///
/// \f$ \chi^{2} \f$ test for comparison two (unweighted) histograms:
/// Let us consider two  histograms with the  same binning and the  number
/// of bins equal to r. Let us denote the number of events in the ith bin
/// in the first histogram as ni and as mi in the second one. The total
/// number of events in the first histogram is equal to:
/// \f[
///  N = \sum_{i=1}^{r} n_{i}
/// \f]
/// and
/// \f[
///  M = \sum_{i=1}^{r} m_{i}
/// \f]
/// in the second histogram. The hypothesis of identity (homogeneity) [3]
/// is that the two histograms represent random values with identical
/// distributions. It is equivalent that there exist r constants p1,...,pr,
/// such that
/// \f[
///\sum_{i=1}^{r} p_{i}=1
/// \f]
/// and the probability of belonging to the ith bin for some measured value
/// in both experiments is equal to pi. The number of events in the ith
/// bin is a random variable with a distribution approximated by a Poisson
/// probability distribution
/// \f[
///\frac{e^{-Np_{i}}(Np_{i})^{n_{i}}}{n_{i}!}
/// \f]
///for the first histogram and with distribution
/// \f[
///\frac{e^{-Mp_{i}}(Mp_{i})^{m_{i}}}{m_{i}!}
/// \f]
/// for the second histogram. If the hypothesis of homogeneity is valid,
/// then the  maximum likelihood estimator of pi, i=1,...,r, is
/// \f[
///\hat{p}_{i}= \frac{n_{i}+m_{i}}{N+M}
/// \f]
/// and then
/// \f[
///  X^{2} = \sum_{i=1}^{r}\frac{(n_{i}-N\hat{p}_{i})^{2}}{N\hat{p}_{i}} + \sum_{i=1}^{r}\frac{(m_{i}-M\hat{p}_{i})^{2}}{M\hat{p}_{i}} =\frac{1}{MN} \sum_{i=1}^{r}\frac{(Mn_{i}-Nm_{i})^{2}}{n_{i}+m_{i}}
/// \f]
/// has approximately a \f$ \chi^{2}_{(r-1)} \f$ distribution [3].
/// The comparison procedure can include an analysis of the residuals which
/// is often helpful in identifying the bins of histograms responsible for
/// a significant overall \f$ \chi^{2} \f$ value. Most convenient for
/// analysis are the adjusted (normalized) residuals [4]
/// \f[
///  r_{i} = \frac{n_{i}-N\hat{p}_{i}}{\sqrt{N\hat{p}_{i}}\sqrt{(1-N/(N+M))(1-(n_{i}+m_{i})/(N+M))}}
/// \f]
/// If hypotheses of  homogeneity are valid then residuals ri are
/// approximately independent and identically distributed random variables
/// having N(0,1) distribution. The application of the \f$ \chi^{2} \f$ test has
/// restrictions related to the value of the expected frequencies Npi,
/// Mpi, i=1,...,r. A conservative rule formulated in [5] is that all the
/// expectations must be 1 or greater for both histograms. In practical
/// cases when expected frequencies are not known the estimated expected
/// frequencies \f$ M\hat{p}_{i}, N\hat{p}_{i}, i=1,...,r \f$ can be used.
///
/// #### Unweighted and weighted histograms comparison:
///
/// A simple modification of the ideas described above can be used for the
/// comparison of the usual (unweighted) and weighted histograms. Let us
/// denote the number of events in the ith bin in the unweighted
/// histogram as ni and the common weight of events in the ith bin of the
/// weighted histogram as wi. The total number of events in the
/// unweighted histogram is equal to
///\f[
///  N = \sum_{i=1}^{r} n_{i}
///\f]
/// and the total weight of events in the weighted histogram is equal to
///\f[
///  W = \sum_{i=1}^{r} w_{i}
///\f]
/// Let us formulate the hypothesis of identity of an unweighted histogram
/// to a weighted histogram so that there exist r constants p1,...,pr, such
/// that
///\f[
///  \sum_{i=1}^{r} p_{i} = 1
///\f]
/// for the unweighted histogram. The weight wi is a random variable with a
/// distribution approximated by the normal probability distribution
/// \f$ N(Wp_{i},\sigma_{i}^{2}) \f$ where \f$ \sigma_{i}^{2} \f$ is the variance of the weight wi.
/// If we replace the variance \f$ \sigma_{i}^{2} \f$
/// with estimate \f$ s_{i}^{2} \f$ (sum of squares of weights of
/// events in the ith bin) and the hypothesis of identity is valid, then the
/// maximum likelihood estimator of  pi,i=1,...,r, is
///\f[
///  \hat{p}_{i} = \frac{Ww_{i}-Ns_{i}^{2}+\sqrt{(Ww_{i}-Ns_{i}^{2})^{2}+4W^{2}s_{i}^{2}n_{i}}}{2W^{2}}
///\f]
/// We may then use the test statistic
///\f[
///  X^{2} = \sum_{i=1}^{r} \frac{(n_{i}-N\hat{p}_{i})^{2}}{N\hat{p}_{i}} + \sum_{i=1}^{r} \frac{(w_{i}-W\hat{p}_{i})^{2}}{s_{i}^{2}}
///\f]
/// and it has approximately a \f$ \sigma^{2}_{(r-1)} \f$ distribution [2]. This test, as well
/// as the original one [3], has a restriction on the expected frequencies. The
/// expected frequencies recommended for the weighted histogram is more than 25.
/// The value of the minimal expected frequency can be decreased down to 10 for
/// the case when the weights of the events are close to constant. In the case
/// of a weighted histogram if the number of events is unknown, then we can
/// apply this recommendation for the equivalent number of events as
///\f[
///  n_{i}^{equiv} = \frac{ w_{i}^{2} }{ s_{i}^{2} }
///\f]
/// The minimal expected frequency for an unweighted histogram must be 1. Notice
/// that any usual (unweighted) histogram can be considered as a weighted
/// histogram with events that have constant weights equal to 1.
/// The variance \f$ z_{i}^{2} \f$ of the difference between the weight wi
/// and the estimated expectation value of the weight is approximately equal to:
///\f[
///  z_{i}^{2} = Var(w_{i}-W\hat{p}_{i}) = N\hat{p}_{i}(1-N\hat{p}_{i})\left(\frac{Ws_{i}^{2}}{\sqrt{(Ns_{i}^{2}-w_{i}W)^{2}+4W^{2}s_{i}^{2}n_{i}}}\right)^{2}+\frac{s_{i}^{2}}{4}\left(1+\frac{Ns_{i}^{2}-w_{i}W}{\sqrt{(Ns_{i}^{2}-w_{i}W)^{2}+4W^{2}s_{i}^{2}n_{i}}}\right)^{2}
///\f]
/// The  residuals
///\f[
///  r_{i} = \frac{w_{i}-W\hat{p}_{i}}{z_{i}}
///\f]
/// have approximately a normal distribution with mean equal to 0 and standard
/// deviation  equal to 1.
///
/// #### Two weighted histograms comparison:
///
/// Let us denote the common  weight of events of the ith bin in the first
/// histogram as w1i and as w2i in the second one. The total weight of events
/// in the first histogram is equal to
///\f[
///  W_{1} = \sum_{i=1}^{r} w_{1i}
///\f]
/// and
///\f[
///  W_{2} = \sum_{i=1}^{r} w_{2i}
///\f]
/// in the second histogram. Let us formulate the hypothesis of identity of
/// weighted histograms so that there exist r constants p1,...,pr, such that
///\f[
///  \sum_{i=1}^{r} p_{i} = 1
///\f]
/// and also expectation value of weight w1i equal to W1pi and expectation value
/// of weight w2i equal to W2pi. Weights in both the histograms are random
/// variables with distributions which can be approximated by a normal
/// probability distribution \f$ N(W_{1}p_{i},\sigma_{1i}^{2}) \f$ for the first histogram
/// and by a distribution \f$ N(W_{2}p_{i},\sigma_{2i}^{2}) \f$ for the second.
/// Here \f$ \sigma_{1i}^{2} \f$ and \f$ \sigma_{2i}^{2} \f$ are the variances
/// of w1i and w2i with estimators \f$ s_{1i}^{2} \f$ and \f$ s_{2i}^{2} \f$ respectively.
/// If the hypothesis of identity is valid, then the maximum likelihood and
/// Least Square Method estimator of pi,i=1,...,r, is
///\f[
///  \hat{p}_{i} = \frac{w_{1i}W_{1}/s_{1i}^{2}+w_{2i}W_{2} /s_{2i}^{2}}{W_{1}^{2}/s_{1i}^{2}+W_{2}^{2}/s_{2i}^{2}}
///\f]
/// We may then use the test statistic
///\f[
/// X^{2} = \sum_{i=1}^{r} \frac{(w_{1i}-W_{1}\hat{p}_{i})^{2}}{s_{1i}^{2}} + \sum_{i=1}^{r} \frac{(w_{2i}-W_{2}\hat{p}_{i})^{2}}{s_{2i}^{2}} = \sum_{i=1}^{r} \frac{(W_{1}w_{2i}-W_{2}w_{1i})^{2}}{W_{1}^{2}s_{2i}^{2}+W_{2}^{2}s_{1i}^{2}}
///\f]
/// and it has approximately a \f$ \chi^{2}_{(r-1)} \f$ distribution [2].
/// The normalized or studentised residuals [6]
///\f[
///  r_{i} = \frac{w_{1i}-W_{1}\hat{p}_{i}}{s_{1i}\sqrt{1 - \frac{1}{(1+W_{2}^{2}s_{1i}^{2}/W_{1}^{2}s_{2i}^{2})}}}
///\f]
/// have approximately a normal distribution with mean equal to 0 and standard
/// deviation 1. A recommended minimal expected frequency is equal to 10 for
/// the proposed test.
///
/// #### Numerical examples:
///
/// The method described herein is now illustrated with an example.
/// We take a distribution
///\f[
/// \phi(x) = \frac{2}{(x-10)^{2}+1} + \frac{1}{(x-14)^{2}+1}       (1)
///\f]
/// defined on the interval [4,16]. Events distributed according to the formula
/// (1) are simulated to create the unweighted histogram. Uniformly distributed
/// events are simulated for the weighted histogram with weights calculated by
/// formula (1). Each histogram has the same number of bins: 20. Fig.1 shows
/// the result of comparison of the unweighted histogram with 200 events
/// (minimal expected frequency equal to one) and the weighted histogram with
/// 500 events (minimal expected frequency equal to 25)
/// Begin_Macro
/// ../../../tutorials/math/chi2test.C
/// End_Macro
/// Fig 1. An example of comparison of the unweighted histogram with 200 events
/// and the weighted histogram with 500 events:
///   1. unweighted histogram;
///   2. weighted histogram;
///   3. normalized residuals plot;
///   4. normal Q-Q plot of residuals.
///
/// The value of the test statistic \f$ \chi^{2} \f$ is equal to
/// 21.09 with p-value equal to 0.33, therefore the hypothesis of identity of
/// the two histograms can be accepted for 0.05 significant level. The behavior
/// of the normalized residuals plot (see Fig. 1c) and the normal Q-Q plot
/// (see Fig. 1d) of residuals are regular and we cannot identify the outliers
/// or bins with a big influence on \f$ \chi^{2} \f$.
///
/// The second example presents the same two histograms but 17 events was added
/// to content of bin number 15 in unweighted histogram. Fig.2 shows the result
/// of comparison of the unweighted histogram with 217 events (minimal expected
/// frequency equal to one) and the weighted histogram with 500 events (minimal
/// expected frequency equal to 25)
/// Begin_Macro
/// ../../../tutorials/math/chi2test.C(17)
/// End_Macro
/// Fig 2. An example of comparison of the unweighted histogram with 217 events
/// and the weighted histogram with 500 events:
///   1. unweighted histogram;
///   2. weighted histogram;
///   3. normalized residuals plot;
///   4. normal Q-Q plot of residuals.
///
/// The value of the test statistic \f$ \chi^{2} \f$ is equal to
/// 32.33 with p-value equal to 0.029, therefore the hypothesis of identity of
/// the two histograms is rejected for 0.05 significant level. The behavior of
/// the normalized residuals plot (see Fig. 2c) and the normal Q-Q plot (see
/// Fig. 2d) of residuals are not regular and we can identify the outlier or
/// bin with a big influence on \f$ \chi^{2} \f$.
///
/// #### References:
///
///  - [1] Pearson, K., 1904. On the Theory of Contingency and Its Relation to
///    Association and Normal Correlation. Drapers' Co. Memoirs, Biometric
///    Series No. 1, London.
///  - [2] Gagunashvili, N., 2006. \f$ \sigma^{2} \f$ test for comparison
///    of weighted and unweighted histograms. Statistical Problems in Particle
///    Physics, Astrophysics and Cosmology, Proceedings of PHYSTAT05,
///    Oxford, UK, 12-15 September 2005, Imperial College Press, London, 43-44.
///    Gagunashvili,N., Comparison of weighted and unweighted histograms,
///    arXiv:physics/0605123, 2006.
///  - [3] Cramer, H., 1946. Mathematical methods of statistics.
///    Princeton University Press, Princeton.
///  - [4] Haberman, S.J., 1973. The analysis of residuals in cross-classified tables.
///    Biometrics 29, 205-220.
///  - [5] Lewontin, R.C. and Felsenstein, J., 1965. The robustness of homogeneity
///    test in 2xN tables. Biometrics 21, 19-33.
///  - [6] Seber, G.A.F., Lee, A.J., 2003, Linear Regression Analysis.
///    John Wiley & Sons Inc., New York.

Double_t TH1::Chi2Test(const TH1* h2, Option_t *option, Double_t *res) const
{
   Double_t chi2 = 0;
   Int_t ndf = 0, igood = 0;

   TString opt = option;
   opt.ToUpper();

   Double_t prob = Chi2TestX(h2,chi2,ndf,igood,option,res);

   if(opt.Contains("P")) {
      printf("Chi2 = %f, Prob = %g, NDF = %d, igood = %d\n", chi2,prob,ndf,igood);
   }
   if(opt.Contains("CHI2/NDF")) {
      if (ndf == 0) return 0;
      return chi2/ndf;
   }
   if(opt.Contains("CHI2")) {
      return chi2;
   }

   return prob;
}

////////////////////////////////////////////////////////////////////////////////
/// The computation routine of the Chisquare test. For the method description,
/// see Chi2Test() function.
///
/// \return p-value
/// \param[in] h2 the second histogram
/// \param[in] option
///  - "UU" = experiment experiment comparison (unweighted-unweighted)
///  - "UW" = experiment MC comparison (unweighted-weighted). Note that the first
///        histogram should be unweighted
///  - "WW" = MC MC comparison (weighted-weighted)
///  - "NORM" = if one or both histograms is scaled
///  - "OF" = overflows included
///  - "UF" = underflows included
///      by default underflows and overflows are not included
/// \param[out] igood test output
///    - igood=0 - no problems
///    - For unweighted unweighted  comparison
///      - igood=1'There is a bin in the 1st histogram with less than 1 event'
///      - igood=2'There is a bin in the 2nd histogram with less than 1 event'
///      - igood=3'when the conditions for igood=1 and igood=2 are satisfied'
///    - For  unweighted weighted  comparison
///      - igood=1'There is a bin in the 1st histogram with less then 1 event'
///      - igood=2'There is a bin in the 2nd histogram with less then 10 effective number of events'
///      - igood=3'when the conditions for igood=1 and igood=2 are satisfied'
///    - For  weighted weighted  comparison
///      - igood=1'There is a bin in the 1st  histogram with less then 10 effective
///        number of events'
///      - igood=2'There is a bin in the 2nd  histogram with less then 10 effective
///        number of events'
///      - igood=3'when the conditions for igood=1 and igood=2 are satisfied'
/// \param[out] chi2 chisquare of the test
/// \param[out] ndf number of degrees of freedom (important, when both histograms have the same empty bins)
/// \param[out] res normalized residuals for further analysis

Double_t TH1::Chi2TestX(const TH1* h2,  Double_t &chi2, Int_t &ndf, Int_t &igood, Option_t *option,  Double_t *res) const
{

   Int_t i_start, i_end;
   Int_t j_start, j_end;
   Int_t k_start, k_end;

   Double_t sum1 = 0.0, sumw1 = 0.0;
   Double_t sum2 = 0.0, sumw2 = 0.0;

   chi2 = 0.0;
   ndf = 0;

   TString opt = option;
   opt.ToUpper();

   if (fBuffer) const_cast<TH1*>(this)->BufferEmpty();

   const TAxis *xaxis1 = GetXaxis();
   const TAxis *xaxis2 = h2->GetXaxis();
   const TAxis *yaxis1 = GetYaxis();
   const TAxis *yaxis2 = h2->GetYaxis();
   const TAxis *zaxis1 = GetZaxis();
   const TAxis *zaxis2 = h2->GetZaxis();

   Int_t nbinx1 = xaxis1->GetNbins();
   Int_t nbinx2 = xaxis2->GetNbins();
   Int_t nbiny1 = yaxis1->GetNbins();
   Int_t nbiny2 = yaxis2->GetNbins();
   Int_t nbinz1 = zaxis1->GetNbins();
   Int_t nbinz2 = zaxis2->GetNbins();

   //check dimensions
   if (this->GetDimension() != h2->GetDimension() ){
      Error("Chi2TestX","Histograms have different dimensions.");
      return 0.0;
   }

   //check number of channels
   if (nbinx1 != nbinx2) {
      Error("Chi2TestX","different number of x channels");
   }
   if (nbiny1 != nbiny2) {
      Error("Chi2TestX","different number of y channels");
   }
   if (nbinz1 != nbinz2) {
      Error("Chi2TestX","different number of z channels");
   }

   //check for ranges
   i_start = j_start = k_start = 1;
   i_end = nbinx1;
   j_end = nbiny1;
   k_end = nbinz1;

   if (xaxis1->TestBit(TAxis::kAxisRange)) {
      i_start = xaxis1->GetFirst();
      i_end   = xaxis1->GetLast();
   }
   if (yaxis1->TestBit(TAxis::kAxisRange)) {
      j_start = yaxis1->GetFirst();
      j_end   = yaxis1->GetLast();
   }
   if (zaxis1->TestBit(TAxis::kAxisRange)) {
      k_start = zaxis1->GetFirst();
      k_end   = zaxis1->GetLast();
   }


   if (opt.Contains("OF")) {
      if (GetDimension() == 3) k_end = ++nbinz1;
      if (GetDimension() >= 2) j_end = ++nbiny1;
      if (GetDimension() >= 1) i_end = ++nbinx1;
   }

   if (opt.Contains("UF")) {
      if (GetDimension() == 3) k_start = 0;
      if (GetDimension() >= 2) j_start = 0;
      if (GetDimension() >= 1) i_start = 0;
   }

   ndf = (i_end - i_start + 1) * (j_end - j_start + 1) * (k_end - k_start + 1) - 1;

   Bool_t comparisonUU = opt.Contains("UU");
   Bool_t comparisonUW = opt.Contains("UW");
   Bool_t comparisonWW = opt.Contains("WW");
   Bool_t scaledHistogram  = opt.Contains("NORM");

   if (scaledHistogram && !comparisonUU) {
      Info("Chi2TestX", "NORM option should be used together with UU option. It is ignored");
   }

   // look at histo global bin content and effective entries
   Stat_t s[kNstat];
   GetStats(s);// s[1] sum of squares of weights, s[0] sum of weights
   Double_t sumBinContent1 = s[0];
   Double_t effEntries1 = (s[1] ? s[0] * s[0] / s[1] : 0.0);

   h2->GetStats(s);// s[1] sum of squares of weights, s[0] sum of weights
   Double_t sumBinContent2 = s[0];
   Double_t effEntries2 = (s[1] ? s[0] * s[0] / s[1] : 0.0);

   if (!comparisonUU && !comparisonUW && !comparisonWW ) {
      // deduce automatically from type of histogram
      if (TMath::Abs(sumBinContent1 - effEntries1) < 1) {
         if ( TMath::Abs(sumBinContent2 - effEntries2) < 1) comparisonUU = true;
         else comparisonUW = true;
      }
      else comparisonWW = true;
   }
   // check unweighted histogram
   if (comparisonUW) {
      if (TMath::Abs(sumBinContent1 - effEntries1) >= 1) {
         Warning("Chi2TestX","First histogram is not unweighted and option UW has been requested");
      }
   }
   if ( (!scaledHistogram && comparisonUU)   ) {
      if ( ( TMath::Abs(sumBinContent1 - effEntries1) >= 1) || (TMath::Abs(sumBinContent2 - effEntries2) >= 1) ) {
         Warning("Chi2TestX","Both histograms are not unweighted and option UU has been requested");
      }
   }


   //get number of events in histogram
   if (comparisonUU && scaledHistogram) {
      for (Int_t i = i_start; i <= i_end; ++i) {
         for (Int_t j = j_start; j <= j_end; ++j) {
            for (Int_t k = k_start; k <= k_end; ++k) {

               Int_t bin = GetBin(i, j, k);

               Double_t cnt1 = RetrieveBinContent(bin);
               Double_t cnt2 = h2->RetrieveBinContent(bin);
               Double_t e1sq = GetBinErrorSqUnchecked(bin);
               Double_t e2sq = h2->GetBinErrorSqUnchecked(bin);

               if (e1sq > 0.0) cnt1 = TMath::Floor(cnt1 * cnt1 / e1sq + 0.5); // avoid rounding errors
               else cnt1 = 0.0;

               if (e2sq > 0.0) cnt2 = TMath::Floor(cnt2 * cnt2 / e2sq + 0.5); // avoid rounding errors
               else cnt2 = 0.0;

               // sum contents
               sum1 += cnt1;
               sum2 += cnt2;
               sumw1 += e1sq;
               sumw2 += e2sq;
            }
         }
      }
      if (sumw1 <= 0.0 || sumw2 <= 0.0) {
         Error("Chi2TestX", "Cannot use option NORM when one histogram has all zero errors");
         return 0.0;
      }

   } else {
      for (Int_t i = i_start; i <= i_end; ++i) {
         for (Int_t j = j_start; j <= j_end; ++j) {
            for (Int_t k = k_start; k <= k_end; ++k) {

               Int_t bin = GetBin(i, j, k);

               sum1 += RetrieveBinContent(bin);
               sum2 += h2->RetrieveBinContent(bin);

               if ( comparisonWW ) sumw1 += GetBinErrorSqUnchecked(bin);
               if ( comparisonUW || comparisonWW ) sumw2 += h2->GetBinErrorSqUnchecked(bin);
            }
         }
      }
   }
   //checks that the histograms are not empty
   if (sum1 == 0.0 || sum2 == 0.0) {
      Error("Chi2TestX","one histogram is empty");
      return 0.0;
   }

   if ( comparisonWW  && ( sumw1 <= 0.0 && sumw2 <= 0.0 ) ){
      Error("Chi2TestX","Hist1 and Hist2 have both all zero errors\n");
      return 0.0;
   }

   //THE TEST
   Int_t m = 0, n = 0;

   //Experiment - experiment comparison
   if (comparisonUU) {
      Double_t sum = sum1 + sum2;
      for (Int_t i = i_start; i <= i_end; ++i) {
         for (Int_t j = j_start; j <= j_end; ++j) {
            for (Int_t k = k_start; k <= k_end; ++k) {

               Int_t bin = GetBin(i, j, k);

               Double_t cnt1 = RetrieveBinContent(bin);
               Double_t cnt2 = h2->RetrieveBinContent(bin);

               if (scaledHistogram) {
                  // scale bin value to effective bin entries
                  Double_t e1sq = GetBinErrorSqUnchecked(bin);
                  Double_t e2sq = h2->GetBinErrorSqUnchecked(bin);

                  if (e1sq > 0) cnt1 = TMath::Floor(cnt1 * cnt1 / e1sq + 0.5); // avoid rounding errors
                  else cnt1 = 0;

                  if (e2sq > 0) cnt2 = TMath::Floor(cnt2 * cnt2 / e2sq + 0.5); // avoid rounding errors
                  else cnt2 = 0;
               }

               if (Int_t(cnt1) == 0 && Int_t(cnt2) == 0) --ndf;  // no data means one degree of freedom less
               else {

                  Double_t cntsum = cnt1 + cnt2;
                  Double_t nexp1 = cntsum * sum1 / sum;
                  //Double_t nexp2 = binsum*sum2/sum;

                  if (res) res[i - i_start] = (cnt1 - nexp1) / TMath::Sqrt(nexp1);

                  if (cnt1 < 1) ++m;
                  if (cnt2 < 1) ++n;

                  //Habermann correction for residuals
                  Double_t correc = (1. - sum1 / sum) * (1. - cntsum / sum);
                  if (res) res[i - i_start] /= TMath::Sqrt(correc);

                  Double_t delta = sum2 * cnt1 - sum1 * cnt2;
                  chi2 += delta * delta / cntsum;
               }
            }
         }
      }
      chi2 /= sum1 * sum2;

      // flag error only when of the two histogram is zero
      if (m) {
         igood += 1;
         Info("Chi2TestX","There is a bin in h1 with less than 1 event.\n");
      }
      if (n) {
         igood += 2;
         Info("Chi2TestX","There is a bin in h2 with less than 1 event.\n");
      }

      Double_t prob = TMath::Prob(chi2,ndf);
      return prob;

   }

   // unweighted - weighted  comparison
   // case of error = 0 and content not zero is treated without problems by excluding second chi2 sum
   // and can be considered as a data-theory comparison
   if ( comparisonUW ) {
      for (Int_t i = i_start; i <= i_end; ++i) {
         for (Int_t j = j_start; j <= j_end; ++j) {
            for (Int_t k = k_start; k <= k_end; ++k) {

               Int_t bin = GetBin(i, j, k);

               Double_t cnt1 = RetrieveBinContent(bin);
               Double_t cnt2 = h2->RetrieveBinContent(bin);
               Double_t e2sq = h2->GetBinErrorSqUnchecked(bin);

               // case both histogram have zero bin contents
               if (cnt1 * cnt1 == 0 && cnt2 * cnt2 == 0) {
                  --ndf;  //no data means one degree of freedom less
                  continue;
               }

               // case weighted histogram has zero bin content and error
               if (cnt2 * cnt2 == 0 && e2sq == 0) {
                  if (sumw2 > 0) {
                     // use as approximated  error as 1 scaled by a scaling ratio
                     // estimated from the total sum weight and sum weight squared
                     e2sq = sumw2 / sum2;
                  }
                  else {
                     // return error because infinite discrepancy here:
                     // bin1 != 0 and bin2 =0 in a histogram with all errors zero
                     Error("Chi2TestX","Hist2 has in bin (%d,%d,%d) zero content and zero errors\n", i, j, k);
                     chi2 = 0; return 0;
                  }
               }

               if (cnt1 < 1) m++;
               if (e2sq > 0 && cnt2 * cnt2 / e2sq < 10) n++;

               Double_t var1 = sum2 * cnt2 - sum1 * e2sq;
               Double_t var2 = var1 * var1 + 4. * sum2 * sum2 * cnt1 * e2sq;

               // if cnt1 is zero and cnt2 = 1 and sum1 = sum2 var1 = 0 && var2 == 0
               // approximate by incrementing cnt1
               // LM (this need to be fixed for numerical errors)
               while (var1 * var1 + cnt1 == 0 || var1 + var2 == 0) {
                  sum1++;
                  cnt1++;
                  var1 = sum2 * cnt2 - sum1 * e2sq;
                  var2 = var1 * var1 + 4. * sum2 * sum2 * cnt1 * e2sq;
               }
               var2 = TMath::Sqrt(var2);

               while (var1 + var2 == 0) {
                  sum1++;
                  cnt1++;
                  var1 = sum2 * cnt2 - sum1 * e2sq;
                  var2 = var1 * var1 + 4. * sum2 * sum2 * cnt1 * e2sq;
                  while (var1 * var1 + cnt1 == 0 || var1 + var2 == 0) {
                     sum1++;
                     cnt1++;
                     var1 = sum2 * cnt2 - sum1 * e2sq;
                     var2 = var1 * var1 + 4. * sum2 * sum2 * cnt1 * e2sq;
                  }
                  var2 = TMath::Sqrt(var2);
               }

               Double_t probb = (var1 + var2) / (2. * sum2 * sum2);

               Double_t nexp1 = probb * sum1;
               Double_t nexp2 = probb * sum2;

               Double_t delta1 = cnt1 - nexp1;
               Double_t delta2 = cnt2 - nexp2;

               chi2 += delta1 * delta1 / nexp1;

               if (e2sq > 0) {
                  chi2 += delta2 * delta2 / e2sq;
               }

               if (res) {
                  if (e2sq > 0) {
                     Double_t temp1 = sum2 * e2sq / var2;
                     Double_t temp2 = 1.0 + (sum1 * e2sq - sum2 * cnt2) / var2;
                     temp2 = temp1 * temp1 * sum1 * probb * (1.0 - probb) + temp2 * temp2 * e2sq / 4.0;
                     // invert sign here
                     res[i - i_start] = - delta2 / TMath::Sqrt(temp2);
                  }
                  else
                     res[i - i_start] = delta1 / TMath::Sqrt(nexp1);
               }
            }
         }
      }

      if (m) {
         igood += 1;
         Info("Chi2TestX","There is a bin in h1 with less than 1 event.\n");
      }
      if (n) {
         igood += 2;
         Info("Chi2TestX","There is a bin in h2 with less than 10 effective events.\n");
      }

      Double_t prob = TMath::Prob(chi2, ndf);

      return prob;
   }

   // weighted - weighted  comparison
   if (comparisonWW) {
      for (Int_t i = i_start; i <= i_end; ++i) {
         for (Int_t j = j_start; j <= j_end; ++j) {
            for (Int_t k = k_start; k <= k_end; ++k) {

               Int_t bin = GetBin(i, j, k);
               Double_t cnt1 = RetrieveBinContent(bin);
               Double_t cnt2 = h2->RetrieveBinContent(bin);
               Double_t e1sq = GetBinErrorSqUnchecked(bin);
               Double_t e2sq = h2->GetBinErrorSqUnchecked(bin);

               // case both histogram have zero bin contents
               // (use square of content to avoid numerical errors)
                if (cnt1 * cnt1 == 0 && cnt2 * cnt2 == 0) {
                   --ndf;  //no data means one degree of freedom less
                   continue;
                }

                if (e1sq == 0 && e2sq == 0) {
                   // cannot treat case of booth histogram have zero zero errors
                  Error("Chi2TestX","h1 and h2 both have bin %d,%d,%d with all zero errors\n", i,j,k);
                  chi2 = 0; return 0;
               }

               Double_t sigma = sum1 * sum1 * e2sq + sum2 * sum2 * e1sq;
               Double_t delta = sum2 * cnt1 - sum1 * cnt2;
               chi2 += delta * delta / sigma;

               if (res) {
                  Double_t temp = cnt1 * sum1 * e2sq + cnt2 * sum2 * e1sq;
                  Double_t probb = temp / sigma;
                  Double_t z = 0;
                  if (e1sq > e2sq) {
                     Double_t d1 = cnt1 - sum1 * probb;
                     Double_t s1 = e1sq * ( 1. - e2sq * sum1 * sum1 / sigma );
                     z = d1 / TMath::Sqrt(s1);
                  }
                  else {
                     Double_t d2 = cnt2 - sum2 * probb;
                     Double_t s2 = e2sq * ( 1. - e1sq * sum2 * sum2 / sigma );
                     z = -d2 / TMath::Sqrt(s2);
                  }
                  res[i - i_start] = z;
               }

               if (e1sq > 0 && cnt1 * cnt1 / e1sq < 10) m++;
               if (e2sq > 0 && cnt2 * cnt2 / e2sq < 10) n++;
            }
         }
      }
      if (m) {
         igood += 1;
         Info("Chi2TestX","There is a bin in h1 with less than 10 effective events.\n");
      }
      if (n) {
         igood += 2;
         Info("Chi2TestX","There is a bin in h2 with less than 10 effective events.\n");
      }
      Double_t prob = TMath::Prob(chi2, ndf);
      return prob;
   }
   return 0;
}
////////////////////////////////////////////////////////////////////////////////
/// Compute and return the chisquare of this histogram with respect to a function
/// The chisquare is computed by weighting each histogram point by the bin error
/// By default the full range of the histogram is used.
/// Use option "R" for restricting the chisquare calculation to the given range of the function
/// Use option "L" for using the chisquare based on the poisson likelihood (Baker-Cousins Chisquare)

Double_t TH1::Chisquare(TF1 * func, Option_t *option) const
{
   if (!func) {
      Error("Chisquare","Function pointer is Null - return -1");
      return -1;
   }

   TString opt(option); opt.ToUpper();
   bool useRange = opt.Contains("R");
   bool usePL = opt.Contains("L");

   return ROOT::Fit::Chisquare(*this, *func, useRange, usePL);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all the content from the underflow and overflow bins, without changing the number of entries
/// After calling this method, every undeflow and overflow bins will have content 0.0
/// The Sumw2 is also cleared, since there is no more content in the bins

void TH1::ClearUnderflowAndOverflow()
{
   for (Int_t bin = 0; bin < fNcells; ++bin)
      if (IsBinUnderflow(bin) || IsBinOverflow(bin)) {
         UpdateBinContent(bin, 0.0);
         if (fSumw2.fN) fSumw2.fArray[bin] = 0.0;
      }
}

////////////////////////////////////////////////////////////////////////////////
///  Compute integral (cumulative sum of bins)
///  The result stored in fIntegral is used by the GetRandom functions.
///  This function is automatically called by GetRandom when the fIntegral
///  array does not exist or when the number of entries in the histogram
///  has changed since the previous call to GetRandom.
///  The resulting integral is normalized to 1
///  If the routine is called with the onlyPositive flag set an error will
///  be produced in case of negative bin content and a NaN value returned

Double_t TH1::ComputeIntegral(Bool_t onlyPositive)
{
   if (fBuffer) BufferEmpty();

   // delete previously computed integral (if any)
   if (fIntegral) delete [] fIntegral;

   //   - Allocate space to store the integral and compute integral
   Int_t nbinsx = GetNbinsX();
   Int_t nbinsy = GetNbinsY();
   Int_t nbinsz = GetNbinsZ();
   Int_t nbins  = nbinsx * nbinsy * nbinsz;

   fIntegral = new Double_t[nbins + 2];
   Int_t ibin = 0; fIntegral[ibin] = 0;

   for (Int_t binz=1; binz <= nbinsz; ++binz) {
      for (Int_t biny=1; biny <= nbinsy; ++biny) {
         for (Int_t binx=1; binx <= nbinsx; ++binx) {
            ++ibin;
            Double_t y = RetrieveBinContent(GetBin(binx, biny, binz));
            if (onlyPositive && y < 0) {
                 Error("ComputeIntegral","Bin content is negative - return a NaN value");
                 fIntegral[nbins] = TMath::QuietNaN();
                 break;
             }
            fIntegral[ibin] = fIntegral[ibin - 1] + y;
         }
      }
   }

   //   - Normalize integral to 1
   if (fIntegral[nbins] == 0 ) {
      Error("ComputeIntegral", "Integral = zero"); return 0;
   }
   for (Int_t bin=1; bin <= nbins; ++bin)  fIntegral[bin] /= fIntegral[nbins];
   fIntegral[nbins+1] = fEntries;
   return fIntegral[nbins];
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the array of bins integral.
/// if the pointer fIntegral is null, TH1::ComputeIntegral is called
/// The array dimension is the number of bins in the histograms
/// including underflow and overflow (fNCells)
/// the last value integral[fNCells] is set to the number of entries of
/// the histogram

Double_t *TH1::GetIntegral()
{
   if (!fIntegral) ComputeIntegral();
   return fIntegral;
}

////////////////////////////////////////////////////////////////////////////////
///  Return a pointer to an histogram containing the cumulative content.
///  The cumulative can be computed both in the forward (default) or backward
///  direction; the name of the new histogram is constructed from
///  the name of this histogram with the suffix "suffix" appended provided
///  by the user. If not provided a default suffix="_cumulative" is used.
///
/// The cumulative distribution is formed by filling each bin of the
/// resulting histogram with the sum of that bin and all previous
/// (forward == kTRUE) or following (forward = kFALSE) bins.
///
/// Note: while cumulative distributions make sense in one dimension, you
/// may not be getting what you expect in more than 1D because the concept
/// of a cumulative distribution is much trickier to define; make sure you
/// understand the order of summation before you use this method with
/// histograms of dimension >= 2.
///
/// Note 2: By default the cumulative is computed from bin 1 to Nbins
/// If an axis range is set, values between the minimum and maximum of the range
/// are set.
/// Setting an axis range can also be used for including underflow and overflow in
/// the cumulative (e.g. by setting h->GetXaxis()->SetRange(0, h->GetNbinsX()+1); )
///

TH1 *TH1::GetCumulative(Bool_t forward, const char* suffix) const
{
   const Int_t firstX = fXaxis.GetFirst();
   const Int_t lastX  = fXaxis.GetLast();
   const Int_t firstY = (fDimension > 1) ? fYaxis.GetFirst() : 1;
   const Int_t lastY = (fDimension > 1) ? fYaxis.GetLast() : 1;
   const Int_t firstZ = (fDimension > 1) ? fZaxis.GetFirst() : 1;
   const Int_t lastZ = (fDimension > 1) ? fZaxis.GetLast() : 1;

   TH1* hintegrated = (TH1*) Clone(fName + suffix);
   hintegrated->Reset();
   Double_t sum = 0.;
   Double_t esum = 0;
   if (forward) { // Forward computation
      for (Int_t binz = firstZ; binz <= lastZ; ++binz) {
         for (Int_t biny = firstY; biny <= lastY; ++biny) {
            for (Int_t binx = firstX; binx <= lastX; ++binx) {
               const Int_t bin = hintegrated->GetBin(binx, biny, binz);
               sum += RetrieveBinContent(bin);
               hintegrated->AddBinContent(bin, sum);
               if (fSumw2.fN) {
                  esum += GetBinErrorSqUnchecked(bin);
                  fSumw2.fArray[bin] = esum;
               }
            }
         }
      }
   } else { // Backward computation
      for (Int_t binz = lastZ; binz >= firstZ; --binz) {
         for (Int_t biny = lastY; biny >= firstY; --biny) {
            for (Int_t binx = lastX; binx >= firstX; --binx) {
               const Int_t bin = hintegrated->GetBin(binx, biny, binz);
               sum += RetrieveBinContent(bin);
               hintegrated->AddBinContent(bin, sum);
               if (fSumw2.fN) {
                  esum += GetBinErrorSqUnchecked(bin);
                  fSumw2.fArray[bin] = esum;
               }
            }
         }
      }
   }
   return hintegrated;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this histogram structure to newth1.
///
/// Note that this function does not copy the list of associated functions.
/// Use TObject::Clone to make a full copy of an histogram.
///
/// Note also that the histogram it will be created in gDirectory (if AddDirectoryStatus()=true)
/// or will not be added to any directory if  AddDirectoryStatus()=false
/// independently of the current directory stored in the original histogram

void TH1::Copy(TObject &obj) const
{
   if (((TH1&)obj).fDirectory) {
      // We are likely to change the hash value of this object
      // with TNamed::Copy, to keep things correct, we need to
      // clean up its existing entries.
      ((TH1&)obj).fDirectory->Remove(&obj);
      ((TH1&)obj).fDirectory = 0;
   }
   TNamed::Copy(obj);
   ((TH1&)obj).fDimension = fDimension;
   ((TH1&)obj).fNormFactor= fNormFactor;
   ((TH1&)obj).fNcells    = fNcells;
   ((TH1&)obj).fBarOffset = fBarOffset;
   ((TH1&)obj).fBarWidth  = fBarWidth;
   ((TH1&)obj).fOption    = fOption;
   ((TH1&)obj).fBinStatErrOpt = fBinStatErrOpt;
   ((TH1&)obj).fBufferSize= fBufferSize;
   // copy the Buffer
   // delete first a previously existing buffer
   if (((TH1&)obj).fBuffer != 0)  {
      delete []  ((TH1&)obj).fBuffer;
      ((TH1&)obj).fBuffer = 0;
   }
   if (fBuffer) {
      Double_t *buf = new Double_t[fBufferSize];
      for (Int_t i=0;i<fBufferSize;i++) buf[i] = fBuffer[i];
      // obj.fBuffer has been deleted before
      ((TH1&)obj).fBuffer    = buf;
   }


   TArray* a = dynamic_cast<TArray*>(&obj);
   if (a) a->Set(fNcells);
   for (Int_t i = 0; i < fNcells; i++) ((TH1&)obj).UpdateBinContent(i, RetrieveBinContent(i));

   ((TH1&)obj).fEntries   = fEntries;

   // which will call BufferEmpty(0) and set fBuffer[0] to a  Maybe one should call
   // assignment operator on the TArrayD

   ((TH1&)obj).fTsumw     = fTsumw;
   ((TH1&)obj).fTsumw2    = fTsumw2;
   ((TH1&)obj).fTsumwx    = fTsumwx;
   ((TH1&)obj).fTsumwx2   = fTsumwx2;
   ((TH1&)obj).fMaximum   = fMaximum;
   ((TH1&)obj).fMinimum   = fMinimum;

   TAttLine::Copy(((TH1&)obj));
   TAttFill::Copy(((TH1&)obj));
   TAttMarker::Copy(((TH1&)obj));
   fXaxis.Copy(((TH1&)obj).fXaxis);
   fYaxis.Copy(((TH1&)obj).fYaxis);
   fZaxis.Copy(((TH1&)obj).fZaxis);
   ((TH1&)obj).fXaxis.SetParent(&obj);
   ((TH1&)obj).fYaxis.SetParent(&obj);
   ((TH1&)obj).fZaxis.SetParent(&obj);
   fContour.Copy(((TH1&)obj).fContour);
   fSumw2.Copy(((TH1&)obj).fSumw2);
   //   fFunctions->Copy(((TH1&)obj).fFunctions);
   // when copying an histogram if the AddDirectoryStatus() is true it
   // will be added to gDirectory independently of the fDirectory stored.
   // and if the AddDirectoryStatus() is false it will not be added to
   // any directory (fDirectory = 0)
   if (fgAddDirectory && gDirectory) {
      gDirectory->Append(&obj);
      ((TH1&)obj).fFunctions->UseRWLock();
      ((TH1&)obj).fDirectory = gDirectory;
   } else
      ((TH1&)obj).fDirectory = 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Make a complete copy of the underlying object.  If 'newname' is set,
/// the copy's name will be set to that name.

TObject* TH1::Clone(const char* newname) const
{
   TH1* obj = (TH1*)IsA()->GetNew()(0);
   Copy(*obj);

   // Now handle the parts that Copy doesn't do
   if(fFunctions) {
      // The Copy above might have published 'obj' to the ListOfCleanups.
      // Clone can call RecursiveRemove, for example via TCheckHashRecursiveRemoveConsistency
      // when dictionary information is initialized, so we need to
      // keep obj->fFunction valid during its execution and
      // protect the update with the write lock.

      // Reset stats parent - else cloning the stats will clone this histogram, too.
      auto oldstats = dynamic_cast<TVirtualPaveStats*>(fFunctions->FindObject("stats"));
      TObject *oldparent = nullptr;
      if (oldstats) {
         oldparent = oldstats->GetParent();
         oldstats->SetParent(nullptr);
      }

      auto newlist = (TList*)fFunctions->Clone();

      if (oldstats)
         oldstats->SetParent(oldparent);
      auto newstats = dynamic_cast<TVirtualPaveStats*>(obj->fFunctions->FindObject("stats"));
      if (newstats)
         newstats->SetParent(obj);

      auto oldlist = obj->fFunctions;
      {
         R__WRITE_LOCKGUARD(ROOT::gCoreMutex);
         obj->fFunctions = newlist;
      }
      delete oldlist;
   }
   if(newname && strlen(newname) ) {
      obj->SetName(newname);
   }
   return obj;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform the automatic addition of the histogram to the given directory
///
/// Note this function is called in place when the semantic requires
/// this object to be added to a directory (I.e. when being read from
/// a TKey or being Cloned)

void TH1::DirectoryAutoAdd(TDirectory *dir)
{
   Bool_t addStatus = TH1::AddDirectoryStatus();
   if (addStatus) {
      SetDirectory(dir);
      if (dir) {
         ResetBit(kCanDelete);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to a line.
///
///  Compute the closest distance of approach from point px,py to elements
///  of an histogram.
///  The distance is computed in pixels units.
///
///  #### Algorithm:
///  Currently, this simple model computes the distance from the mouse
///  to the histogram contour only.

Int_t TH1::DistancetoPrimitive(Int_t px, Int_t py)
{
   if (!fPainter) return 9999;
   return fPainter->DistancetoPrimitive(px,py);
}

////////////////////////////////////////////////////////////////////////////////
/// Performs the operation: `this = this/(c1*f1)`
/// if errors are defined (see TH1::Sumw2), errors are also recalculated.
///
/// Only bins inside the function range are recomputed.
/// IMPORTANT NOTE: If you intend to use the errors of this histogram later
/// you should call Sumw2 before making this operation.
/// This is particularly important if you fit the histogram after TH1::Divide
///
/// The function return kFALSE if the divide operation failed

Bool_t TH1::Divide(TF1 *f1, Double_t c1)
{
   if (!f1) {
      Error("Divide","Attempt to divide by a non-existing function");
      return kFALSE;
   }

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   Int_t nx = GetNbinsX() + 2; // normal bins + uf / of
   Int_t ny = GetNbinsY() + 2;
   Int_t nz = GetNbinsZ() + 2;
   if (fDimension < 2) ny = 1;
   if (fDimension < 3) nz = 1;


   SetMinimum();
   SetMaximum();

   //   - Loop on bins (including underflows/overflows)
   Int_t bin, binx, biny, binz;
   Double_t cu, w;
   Double_t xx[3];
   Double_t *params = 0;
   f1->InitArgs(xx,params);
   for (binz = 0; binz < nz; ++binz) {
      xx[2] = fZaxis.GetBinCenter(binz);
      for (biny = 0; biny < ny; ++biny) {
         xx[1] = fYaxis.GetBinCenter(biny);
         for (binx = 0; binx < nx; ++binx) {
            xx[0] = fXaxis.GetBinCenter(binx);
            if (!f1->IsInside(xx)) continue;
            TF1::RejectPoint(kFALSE);
            bin = binx + nx * (biny + ny * binz);
            cu  = c1 * f1->EvalPar(xx);
            if (TF1::RejectedPoint()) continue;
            if (cu) w = RetrieveBinContent(bin) / cu;
            else    w = 0;
            UpdateBinContent(bin, w);
            if (fSumw2.fN) {
               if (cu != 0) fSumw2.fArray[bin] = GetBinErrorSqUnchecked(bin) / (cu * cu);
               else         fSumw2.fArray[bin] = 0;
            }
         }
      }
   }
   ResetStats();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this histogram by h1.
///
/// `this = this/h1`
/// if errors are defined (see TH1::Sumw2), errors are also recalculated.
/// Note that if h1 has Sumw2 set, Sumw2 is automatically called for this
/// if not already set.
/// The resulting errors are calculated assuming uncorrelated histograms.
/// See the other TH1::Divide that gives the possibility to optionally
/// compute binomial errors.
///
/// IMPORTANT NOTE: If you intend to use the errors of this histogram later
/// you should call Sumw2 before making this operation.
/// This is particularly important if you fit the histogram after TH1::Scale
///
/// The function return kFALSE if the divide operation failed

Bool_t TH1::Divide(const TH1 *h1)
{
   if (!h1) {
      Error("Divide", "Input histogram passed does not exist (NULL).");
      return kFALSE;
   }

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   try {
      CheckConsistency(this,h1);
   } catch(DifferentNumberOfBins&) {
      Error("Divide","Cannot divide histograms with different number of bins");
      return kFALSE;
   } catch(DifferentAxisLimits&) {
      Warning("Divide","Dividing histograms with different axis limits");
   } catch(DifferentBinLimits&) {
      Warning("Divide","Dividing histograms with different bin limits");
   } catch(DifferentLabels&) {
      Warning("Divide","Dividing histograms with different labels");
   }

   //    Create Sumw2 if h1 has Sumw2 set
   if (fSumw2.fN == 0 && h1->GetSumw2N() != 0) Sumw2();

   //   - Loop on bins (including underflows/overflows)
   for (Int_t i = 0; i < fNcells; ++i) {
      Double_t c0 = RetrieveBinContent(i);
      Double_t c1 = h1->RetrieveBinContent(i);
      if (c1) UpdateBinContent(i, c0 / c1);
      else UpdateBinContent(i, 0);

      if(fSumw2.fN) {
         if (c1 == 0) { fSumw2.fArray[i] = 0; continue; }
         Double_t c1sq = c1 * c1;
         fSumw2.fArray[i] = (GetBinErrorSqUnchecked(i) * c1sq + h1->GetBinErrorSqUnchecked(i) * c0 * c0) / (c1sq * c1sq);
      }
   }
   ResetStats();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace contents of this histogram by the division of h1 by h2.
///
/// `this = c1*h1/(c2*h2)`
///
/// If errors are defined (see TH1::Sumw2), errors are also recalculated
/// Note that if h1 or h2 have Sumw2 set, Sumw2 is automatically called for this
/// if not already set.
/// The resulting errors are calculated assuming uncorrelated histograms.
/// However, if option ="B" is specified, Binomial errors are computed.
/// In this case c1 and c2 do not make real sense and they are ignored.
///
/// IMPORTANT NOTE: If you intend to use the errors of this histogram later
/// you should call Sumw2 before making this operation.
/// This is particularly important if you fit the histogram after TH1::Divide
///
///  Please note also that in the binomial case errors are calculated using standard
///  binomial statistics, which means when b1 = b2, the error is zero.
///  If you prefer to have efficiency errors not going to zero when the efficiency is 1, you must
///  use the function TGraphAsymmErrors::BayesDivide, which will return an asymmetric and non-zero lower
///  error for the case b1=b2.
///
/// The function return kFALSE if the divide operation failed

Bool_t TH1::Divide(const TH1 *h1, const TH1 *h2, Double_t c1, Double_t c2, Option_t *option)
{

   TString opt = option;
   opt.ToLower();
   Bool_t binomial = kFALSE;
   if (opt.Contains("b")) binomial = kTRUE;
   if (!h1 || !h2) {
      Error("Divide", "At least one of the input histograms passed does not exist (NULL).");
      return kFALSE;
   }

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   try {
      CheckConsistency(h1,h2);
      CheckConsistency(this,h1);
   } catch(DifferentNumberOfBins&) {
      Error("Divide","Cannot divide histograms with different number of bins");
      return kFALSE;
   } catch(DifferentAxisLimits&) {
      Warning("Divide","Dividing histograms with different axis limits");
   } catch(DifferentBinLimits&) {
      Warning("Divide","Dividing histograms with different bin limits");
   }  catch(DifferentLabels&) {
      Warning("Divide","Dividing histograms with different labels");
   }


   if (!c2) {
      Error("Divide","Coefficient of dividing histogram cannot be zero");
      return kFALSE;
   }

   //    Create Sumw2 if h1 or h2 have Sumw2 set, or if binomial errors are explicitly requested
   if (fSumw2.fN == 0 && (h1->GetSumw2N() != 0 || h2->GetSumw2N() != 0 || binomial)) Sumw2();

   SetMinimum();
   SetMaximum();

   //   - Loop on bins (including underflows/overflows)
   for (Int_t i = 0; i < fNcells; ++i) {
      Double_t b1 = h1->RetrieveBinContent(i);
      Double_t b2 = h2->RetrieveBinContent(i);
      if (b2) UpdateBinContent(i, c1 * b1 / (c2 * b2));
      else UpdateBinContent(i, 0);

      if (fSumw2.fN) {
         if (b2 == 0) { fSumw2.fArray[i] = 0; continue; }
         Double_t b1sq = b1 * b1; Double_t b2sq = b2 * b2;
         Double_t c1sq = c1 * c1; Double_t c2sq = c2 * c2;
         Double_t e1sq = h1->GetBinErrorSqUnchecked(i);
         Double_t e2sq = h2->GetBinErrorSqUnchecked(i);
         if (binomial) {
            if (b1 != b2) {
               // in the case of binomial statistics c1 and c2 must be 1 otherwise it does not make sense
               // c1 and c2 are ignored
               //fSumw2.fArray[bin] = TMath::Abs(w*(1-w)/(c2*b2));//this is the formula in Hbook/Hoper1
               //fSumw2.fArray[bin] = TMath::Abs(w*(1-w)/b2);     // old formula from G. Flucke
               // formula which works also for weighted histogram (see http://root-forum.cern.ch/viewtopic.php?t=3753 )
               fSumw2.fArray[i] = TMath::Abs( ( (1. - 2.* b1 / b2) * e1sq  + b1sq * e2sq / b2sq ) / b2sq );
            } else {
               //in case b1=b2 error is zero
               //use  TGraphAsymmErrors::BayesDivide for getting the asymmetric error not equal to zero
               fSumw2.fArray[i] = 0;
            }
         } else {
            fSumw2.fArray[i] = c1sq * c2sq * (e1sq * b2sq + e2sq * b1sq) / (c2sq * c2sq * b2sq * b2sq);
         }
      }
   }
   ResetStats();
   if (binomial)
      // in case of binomial division use denominator for number of entries
      SetEntries ( h2->GetEntries() );

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw this histogram with options.
///
/// Histograms are drawn via the THistPainter class. Each histogram has
/// a pointer to its own painter (to be usable in a multithreaded program).
/// The same histogram can be drawn with different options in different pads.
/// When an histogram drawn in a pad is deleted, the histogram is
/// automatically removed from the pad or pads where it was drawn.
/// If an histogram is drawn in a pad, then filled again, the new status
/// of the histogram will be automatically shown in the pad next time
/// the pad is updated. One does not need to redraw the histogram.
/// To draw the current version of an histogram in a pad, one can use
/// `h->DrawCopy();`
/// This makes a clone of the histogram. Once the clone is drawn, the original
/// histogram may be modified or deleted without affecting the aspect of the
/// clone.
/// By default, TH1::Draw clears the current pad.
///
/// One can use TH1::SetMaximum and TH1::SetMinimum to force a particular
/// value for the maximum or the minimum scale on the plot.
///
/// TH1::UseCurrentStyle can be used to change all histogram graphics
/// attributes to correspond to the current selected style.
/// This function must be called for each histogram.
/// In case one reads and draws many histograms from a file, one can force
/// the histograms to inherit automatically the current graphics style
/// by calling before gROOT->ForceStyle();
///
/// See the THistPainter class for a description of all the drawing options.

void TH1::Draw(Option_t *option)
{
   TString opt1 = option; opt1.ToLower();
   TString opt2 = option;
   Int_t index  = opt1.Index("same");

   // Check if the string "same" is part of a TCutg name.
   if (index>=0) {
      Int_t indb = opt1.Index("[");
      if (indb>=0) {
         Int_t indk = opt1.Index("]");
         if (index>indb && index<indk) index = -1;
      }
   }

   // If there is no pad or an empty pad the "same" option is ignored.
   if (gPad) {
      if (!gPad->IsEditable()) gROOT->MakeDefCanvas();
      if (index>=0) {
         if (gPad->GetX1() == 0   && gPad->GetX2() == 1 &&
             gPad->GetY1() == 0   && gPad->GetY2() == 1 &&
             gPad->GetListOfPrimitives()->GetSize()==0) opt2.Remove(index,4);
      } else {
         //the following statement is necessary in case one attempts to draw
         //a temporary histogram already in the current pad
         if (TestBit(kCanDelete)) gPad->GetListOfPrimitives()->Remove(this);
         gPad->Clear();
      }
      gPad->IncrementPaletteColor(1, opt1);
   } else {
      if (index>=0) opt2.Remove(index,4);
   }

   AppendPad(opt2.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this histogram and Draw in the current pad.
///
/// Once the histogram is drawn into the pad, any further modification
/// using graphics input will be made on the copy of the histogram,
/// and not to the original object.
/// By default a postfix "_copy" is added to the histogram name. Pass an empty postfix in case
/// you want to draw an histogram with the same name
///
/// See Draw for the list of options

TH1 *TH1::DrawCopy(Option_t *option, const char * name_postfix) const
{
   TString opt = option;
   opt.ToLower();
   if (gPad && !opt.Contains("same")) gPad->Clear();
   TString newName = (name_postfix) ?  TString::Format("%s%s",GetName(),name_postfix) : "";
   TH1 *newth1 = (TH1 *)Clone(newName);
   newth1->SetDirectory(0);
   newth1->SetBit(kCanDelete);
   if (gPad) gPad->IncrementPaletteColor(1, opt);

   newth1->AppendPad(option);
   return newth1;
}

////////////////////////////////////////////////////////////////////////////////
///  Draw a normalized copy of this histogram.
///
///  A clone of this histogram is normalized to norm and drawn with option.
///  A pointer to the normalized histogram is returned.
///  The contents of the histogram copy are scaled such that the new
///  sum of weights (excluding under and overflow) is equal to norm.
///  Note that the returned normalized histogram is not added to the list
///  of histograms in the current directory in memory.
///  It is the user's responsibility to delete this histogram.
///  The kCanDelete bit is set for the returned object. If a pad containing
///  this copy is cleared, the histogram will be automatically deleted.
///
///  See Draw for the list of options

TH1 *TH1::DrawNormalized(Option_t *option, Double_t norm) const
{
   Double_t sum = GetSumOfWeights();
   if (sum == 0) {
      Error("DrawNormalized","Sum of weights is null. Cannot normalize histogram: %s",GetName());
      return 0;
   }
   Bool_t addStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   TH1 *h = (TH1*)Clone();
   h->SetBit(kCanDelete);
   // in case of drawing with error options - scale correctly the error
   TString opt(option); opt.ToUpper();
   if (fSumw2.fN == 0) {
      h->Sumw2();
      // do not use in this case the "Error option " for drawing which is enabled by default since the normalized histogram has now errors
      if (opt.IsNull() || opt == "SAME") opt += "HIST";
   }
   h->Scale(norm/sum);
   if (TMath::Abs(fMaximum+1111) > 1e-3) h->SetMaximum(fMaximum*norm/sum);
   if (TMath::Abs(fMinimum+1111) > 1e-3) h->SetMinimum(fMinimum*norm/sum);
   h->Draw(opt);
   TH1::AddDirectory(addStatus);
   return h;
}

////////////////////////////////////////////////////////////////////////////////
/// Display a panel with all histogram drawing options.
///
///   See class TDrawPanelHist for example

void TH1::DrawPanel()
{
   if (!fPainter) {Draw(); if (gPad) gPad->Update();}
   if (fPainter) fPainter->DrawPanel();
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate function f1 at the center of bins of this histogram.
///
///  - If option "R" is specified, the function is evaluated only
///    for the bins included in the function range.
///  - If option "A" is specified, the value of the function is added to the
///    existing bin contents
///  - If option "S" is specified, the value of the function is used to
///    generate a value, distributed according to the Poisson
///    distribution, with f1 as the mean.

void TH1::Eval(TF1 *f1, Option_t *option)
{
   Double_t x[3];
   Int_t range, stat, add;
   if (!f1) return;

   TString opt = option;
   opt.ToLower();
   if (opt.Contains("a")) add   = 1;
   else                   add   = 0;
   if (opt.Contains("s")) stat  = 1;
   else                   stat  = 0;
   if (opt.Contains("r")) range = 1;
   else                   range = 0;

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   Int_t nbinsx  = fXaxis.GetNbins();
   Int_t nbinsy  = fYaxis.GetNbins();
   Int_t nbinsz  = fZaxis.GetNbins();
   if (!add) Reset();

   for (Int_t binz = 1; binz <= nbinsz; ++binz) {
      x[2]  = fZaxis.GetBinCenter(binz);
      for (Int_t biny = 1; biny <= nbinsy; ++biny) {
         x[1]  = fYaxis.GetBinCenter(biny);
         for (Int_t binx = 1; binx <= nbinsx; ++binx) {
            Int_t bin = GetBin(binx,biny,binz);
            x[0]  = fXaxis.GetBinCenter(binx);
            if (range && !f1->IsInside(x)) continue;
            Double_t fu = f1->Eval(x[0], x[1], x[2]);
            if (stat) fu = gRandom->PoissonD(fu);
            AddBinContent(bin, fu);
            if (fSumw2.fN) fSumw2.fArray[bin] += TMath::Abs(fu);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Execute action corresponding to one event.
///
/// This member function is called when a histogram is clicked with the locator
///
/// If Left button clicked on the bin top value, then the content of this bin
/// is modified according to the new position of the mouse when it is released.

void TH1::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   if (fPainter) fPainter->ExecuteEvent(event, px, py);
}

////////////////////////////////////////////////////////////////////////////////
/// This function allows to do discrete Fourier transforms of TH1 and TH2.
/// Available transform types and flags are described below.
///
/// To extract more information about the transform, use the function
/// TVirtualFFT::GetCurrentTransform() to get a pointer to the current
/// transform object.
///
/// \param[out] h_output histogram for the output. If a null pointer is passed, a new histogram is created
///          and returned, otherwise, the provided histogram is used and should be big enough
/// \param[in] option option parameters consists of 3 parts:
/// - option on what to return
///   - "RE" - returns a histogram of the real part of the output
///   - "IM" - returns a histogram of the imaginary part of the output
///   - "MAG"- returns a histogram of the magnitude of the output
///   - "PH" - returns a histogram of the phase of the output
/// - option of transform type
///   - "R2C"  - real to complex transforms - default
///   - "R2HC" - real to halfcomplex (special format of storing output data,
///     results the same as for R2C)
///   - "DHT" - discrete Hartley transform
///     real to real transforms (sine and cosine):
///   - "R2R_0", "R2R_1", "R2R_2", "R2R_3" - discrete cosine transforms of types I-IV
///   - "R2R_4", "R2R_5", "R2R_6", "R2R_7" - discrete sine transforms of types I-IV
///     To specify the type of each dimension of a 2-dimensional real to real
///     transform, use options of form "R2R_XX", for example, "R2R_02" for a transform,
///     which is of type "R2R_0" in 1st dimension and  "R2R_2" in the 2nd.
/// - option of transform flag
///   - "ES" (from "estimate") - no time in preparing the transform, but probably sub-optimal
///     performance
///   - "M" (from "measure")   - some time spend in finding the optimal way to do the transform
///   - "P" (from "patient")   - more time spend in finding the optimal way to do the transform
///   - "EX" (from "exhaustive") - the most optimal way is found
///     This option should be chosen depending on how many transforms of the same size and
///     type are going to be done. Planning is only done once, for the first transform of this
///     size and type. Default is "ES".
///
/// Examples of valid options: "Mag R2C M" "Re R2R_11" "Im R2C ES" "PH R2HC EX"

TH1* TH1::FFT(TH1* h_output, Option_t *option)
{

   Int_t ndim[3];
   ndim[0] = this->GetNbinsX();
   ndim[1] = this->GetNbinsY();
   ndim[2] = this->GetNbinsZ();

   TVirtualFFT *fft;
   TString opt = option;
   opt.ToUpper();
   if (!opt.Contains("2R")){
      if (!opt.Contains("2C") && !opt.Contains("2HC") && !opt.Contains("DHT")) {
         //no type specified, "R2C" by default
         opt.Append("R2C");
      }
      fft = TVirtualFFT::FFT(this->GetDimension(), ndim, opt.Data());
   }
   else {
      //find the kind of transform
      Int_t ind = opt.Index("R2R", 3);
      Int_t *kind = new Int_t[2];
      char t;
      t = opt[ind+4];
      kind[0] = atoi(&t);
      if (h_output->GetDimension()>1) {
         t = opt[ind+5];
         kind[1] = atoi(&t);
      }
      fft = TVirtualFFT::SineCosine(this->GetDimension(), ndim, kind, option);
      delete [] kind;
   }

   if (!fft) return 0;
   Int_t in=0;
   for (Int_t binx = 1; binx<=ndim[0]; binx++) {
      for (Int_t biny=1; biny<=ndim[1]; biny++) {
         for (Int_t binz=1; binz<=ndim[2]; binz++) {
            fft->SetPoint(in, this->GetBinContent(binx, biny, binz));
            in++;
         }
      }
   }
   fft->Transform();
   h_output = TransformHisto(fft, h_output, option);
   return h_output;
}

////////////////////////////////////////////////////////////////////////////////
/// Increment bin with abscissa X by 1.
///
/// if x is less than the low-edge of the first bin, the Underflow bin is incremented
/// if x is equal to or greater than the upper edge of last bin, the Overflow bin is incremented
///
/// If the storage of the sum of squares of weights has been triggered,
/// via the function Sumw2, then the sum of the squares of weights is incremented
/// by 1 in the bin corresponding to x.
///
/// The function returns the corresponding bin number which has its content incremented by 1

Int_t TH1::Fill(Double_t x)
{
   if (fBuffer)  return BufferFill(x,1);

   Int_t bin;
   fEntries++;
   bin =fXaxis.FindBin(x);
   if (bin <0) return -1;
   AddBinContent(bin);
   if (fSumw2.fN) ++fSumw2.fArray[bin];
   if (bin == 0 || bin > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   ++fTsumw;
   ++fTsumw2;
   fTsumwx  += x;
   fTsumwx2 += x*x;
   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Increment bin with abscissa X with a weight w.
///
/// if x is less than the low-edge of the first bin, the Underflow bin is incremented
/// if x is equal to or greater than the upper edge of last bin, the Overflow bin is incremented
///
/// If the weight is not equal to 1, the storage of the sum of squares of
/// weights is automatically triggered and the sum of the squares of weights is incremented
/// by \f$ w^2 \f$ in the bin corresponding to x.
///
/// The function returns the corresponding bin number which has its content incremented by w

Int_t TH1::Fill(Double_t x, Double_t w)
{

   if (fBuffer) return BufferFill(x,w);

   Int_t bin;
   fEntries++;
   bin =fXaxis.FindBin(x);
   if (bin <0) return -1;
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW) )  Sumw2();   // must be called before AddBinContent
   if (fSumw2.fN)  fSumw2.fArray[bin] += w*w;
   AddBinContent(bin, w);
   if (bin == 0 || bin > fXaxis.GetNbins()) {
      if (!GetStatOverflowsBehaviour()) return -1;
   }
   Double_t z= w;
   fTsumw   += z;
   fTsumw2  += z*z;
   fTsumwx  += z*x;
   fTsumwx2 += z*x*x;
   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Increment bin with namex with a weight w
///
/// if x is less than the low-edge of the first bin, the Underflow bin is incremented
/// if x is equal to or greater than the upper edge of last bin, the Overflow bin is incremented
///
/// If the weight is not equal to 1, the storage of the sum of squares of
/// weights is automatically triggered and the sum of the squares of weights is incremented
/// by \f$ w^2 \f$ in the bin corresponding to x.
///
/// The function returns the corresponding bin number which has its content
/// incremented by w.

Int_t TH1::Fill(const char *namex, Double_t w)
{
   Int_t bin;
   fEntries++;
   bin =fXaxis.FindBin(namex);
   if (bin <0) return -1;
   if (!fSumw2.fN && w != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();
   if (fSumw2.fN) fSumw2.fArray[bin] += w*w;
   AddBinContent(bin, w);
   if (bin == 0 || bin > fXaxis.GetNbins()) return -1;
   Double_t z= w;
   fTsumw   += z;
   fTsumw2  += z*z;
   // this make sense if the histogram is not expanding (the x axis cannot be extended)
   if (!fXaxis.CanExtend() || !fXaxis.IsAlphanumeric()) {
      Double_t x = fXaxis.GetBinCenter(bin);
      fTsumwx  += z*x;
      fTsumwx2 += z*x*x;
   }
   return bin;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill this histogram with an array x and weights w.
///
/// \param[in] ntimes number of entries in arrays x and w (array size must be ntimes*stride)
/// \param[in] x array of values to be histogrammed
/// \param[in] w array of weighs
/// \param[in] stride step size through arrays x and w
///
/// If the weight is not equal to 1, the storage of the sum of squares of
/// weights is automatically triggered and the sum of the squares of weights is incremented
/// by \f$ w^2 \f$ in the bin corresponding to x.
/// if w is NULL each entry is assumed a weight=1

void TH1::FillN(Int_t ntimes, const Double_t *x, const Double_t *w, Int_t stride)
{
   //If a buffer is activated, fill buffer
   if (fBuffer) {
      ntimes *= stride;
      Int_t i = 0;
      for (i=0;i<ntimes;i+=stride) {
         if (!fBuffer) break;   // buffer can be deleted in BufferFill when is empty
         if (w) BufferFill(x[i],w[i]);
         else BufferFill(x[i], 1.);
      }
      // fill the remaining entries if the buffer has been deleted
      if (i < ntimes && fBuffer==0) {
         auto weights = w ? &w[i] : nullptr;
         DoFillN((ntimes-i)/stride,&x[i],weights,stride);
      }
      return;
   }
   // call internal method
   DoFillN(ntimes, x, w, stride);
}

////////////////////////////////////////////////////////////////////////////////
/// Internal method to fill histogram content from a vector
/// called directly by TH1::BufferEmpty

void TH1::DoFillN(Int_t ntimes, const Double_t *x, const Double_t *w, Int_t stride)
{
   Int_t bin,i;

   fEntries += ntimes;
   Double_t ww = 1;
   Int_t nbins   = fXaxis.GetNbins();
   ntimes *= stride;
   for (i=0;i<ntimes;i+=stride) {
      bin =fXaxis.FindBin(x[i]);
      if (bin <0) continue;
      if (w) ww = w[i];
      if (!fSumw2.fN && ww != 1.0 && !TestBit(TH1::kIsNotW))  Sumw2();
      if (fSumw2.fN) fSumw2.fArray[bin] += ww*ww;
      AddBinContent(bin, ww);
      if (bin == 0 || bin > nbins) {
         if (!GetStatOverflowsBehaviour()) continue;
      }
      Double_t z= ww;
      fTsumw   += z;
      fTsumw2  += z*z;
      fTsumwx  += z*x[i];
      fTsumwx2 += z*x[i]*x[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill histogram following distribution in function fname.
///
///  @param fname  : Function name used for filling the historam
///  @param ntimes : number of times the histogram is filled
///  @param rng    : (optional) Random number generator used to sample
///
///
/// The distribution contained in the function fname (TF1) is integrated
/// over the channel contents for the bin range of this histogram.
/// It is normalized to 1.
///
/// Getting one random number implies:
///  - Generating a random number between 0 and 1 (say r1)
///  - Look in which bin in the normalized integral r1 corresponds to
///  - Fill histogram channel
///    ntimes random numbers are generated
///
/// One can also call TF1::GetRandom to get a random variate from a function.

void TH1::FillRandom(const char *fname, Int_t ntimes, TRandom * rng)
{
   Int_t bin, binx, ibin, loop;
   Double_t r1, x;
   //   - Search for fname in the list of ROOT defined functions
   TF1 *f1 = (TF1*)gROOT->GetFunction(fname);
   if (!f1) { Error("FillRandom", "Unknown function: %s",fname); return; }

   //   - Allocate temporary space to store the integral and compute integral

   TAxis * xAxis = &fXaxis;

   // in case axis of histogram is not defined use the function axis
   if (fXaxis.GetXmax() <= fXaxis.GetXmin()) {
      Double_t xmin,xmax;
      f1->GetRange(xmin,xmax);
      Info("FillRandom","Using function axis and range [%g,%g]",xmin, xmax);
      xAxis = f1->GetHistogram()->GetXaxis();
   }

   Int_t first  = xAxis->GetFirst();
   Int_t last   = xAxis->GetLast();
   Int_t nbinsx = last-first+1;

   Double_t *integral = new Double_t[nbinsx+1];
   integral[0] = 0;
   for (binx=1;binx<=nbinsx;binx++) {
      Double_t fint = f1->Integral(xAxis->GetBinLowEdge(binx+first-1),xAxis->GetBinUpEdge(binx+first-1), 0.);
      integral[binx] = integral[binx-1] + fint;
   }

   //   - Normalize integral to 1
   if (integral[nbinsx] == 0 ) {
      delete [] integral;
      Error("FillRandom", "Integral = zero"); return;
   }
   for (bin=1;bin<=nbinsx;bin++)  integral[bin] /= integral[nbinsx];

   //   --------------Start main loop ntimes
   for (loop=0;loop<ntimes;loop++) {
      r1 = (rng) ? rng->Rndm() : gRandom->Rndm();
      ibin = TMath::BinarySearch(nbinsx,&integral[0],r1);
      //binx = 1 + ibin;
      //x    = xAxis->GetBinCenter(binx); //this is not OK when SetBuffer is used
      x    = xAxis->GetBinLowEdge(ibin+first)
             +xAxis->GetBinWidth(ibin+first)*(r1-integral[ibin])/(integral[ibin+1] - integral[ibin]);
      Fill(x);
   }
   delete [] integral;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill histogram following distribution in histogram h.
///
///  @param h      : Histogram  pointer used for smpling random number
///  @param ntimes : number of times the histogram is filled
///  @param rng    : (optional) Random number generator used for sampling
///
/// The distribution contained in the histogram h (TH1) is integrated
/// over the channel contents for the bin range of this histogram.
/// It is normalized to 1.
///
/// Getting one random number implies:
///  - Generating a random number between 0 and 1 (say r1)
///  - Look in which bin in the normalized integral r1 corresponds to
///  - Fill histogram channel ntimes random numbers are generated
///
/// SPECIAL CASE when the target histogram has the same binning as the source.
/// in this case we simply use a poisson distribution where
/// the mean value per bin = bincontent/integral.

void TH1::FillRandom(TH1 *h, Int_t ntimes, TRandom * rng)
{
   if (!h) { Error("FillRandom", "Null histogram"); return; }
   if (fDimension != h->GetDimension()) {
      Error("FillRandom", "Histograms with different dimensions"); return;
   }
   if (std::isnan(h->ComputeIntegral(true))) {
      Error("FillRandom", "Histograms contains negative bins, does not represent probabilities");
      return;
   }

   //in case the target histogram has the same binning and ntimes much greater
   //than the number of bins we can use a fast method
   Int_t first  = fXaxis.GetFirst();
   Int_t last   = fXaxis.GetLast();
   Int_t nbins = last-first+1;
   if (ntimes > 10*nbins) {
      try {
         CheckConsistency(this,h);
         Double_t sumw = h->Integral(first,last);
         if (sumw == 0) return;
         Double_t sumgen = 0;
         for (Int_t bin=first;bin<=last;bin++) {
            Double_t mean = h->RetrieveBinContent(bin)*ntimes/sumw;
            Double_t cont = (rng) ? rng->Poisson(mean) : gRandom->Poisson(mean);
            sumgen += cont;
            AddBinContent(bin,cont);
            if (fSumw2.fN) fSumw2.fArray[bin] += cont;
         }

         // fix for the fluctuations in the total number n
         // since we use Poisson instead of multinomial
         // add a correction to have ntimes as generated entries
         Int_t i;
         if (sumgen < ntimes) {
            // add missing entries
            for (i = Int_t(sumgen+0.5); i < ntimes; ++i)
            {
               Double_t x = h->GetRandom();
               Fill(x);
            }
         }
         else if (sumgen > ntimes) {
            // remove extra entries
            i =  Int_t(sumgen+0.5);
            while( i > ntimes) {
               Double_t x = h->GetRandom(rng);
               Int_t ibin = fXaxis.FindBin(x);
               Double_t y = RetrieveBinContent(ibin);
               // skip in case bin is empty
               if (y > 0) {
                  SetBinContent(ibin, y-1.);
                  i--;
               }
            }
         }

         ResetStats();
         return;
      }
      catch(std::exception&) {}  // do nothing
   }
   // case of different axis and not too large ntimes

   if (h->ComputeIntegral() ==0) return;
   Int_t loop;
   Double_t x;
   for (loop=0;loop<ntimes;loop++) {
      x = h->GetRandom();
      Fill(x);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return Global bin number corresponding to x,y,z
///
/// 2-D and 3-D histograms are represented with a one dimensional
/// structure. This has the advantage that all existing functions, such as
/// GetBinContent, GetBinError, GetBinFunction work for all dimensions.
/// This function tries to extend the axis if the given point belongs to an
/// under-/overflow bin AND if CanExtendAllAxes() is true.
///
/// See also TH1::GetBin, TAxis::FindBin and TAxis::FindFixBin

Int_t TH1::FindBin(Double_t x, Double_t y, Double_t z)
{
   if (GetDimension() < 2) {
      return fXaxis.FindBin(x);
   }
   if (GetDimension() < 3) {
      Int_t nx   = fXaxis.GetNbins()+2;
      Int_t binx = fXaxis.FindBin(x);
      Int_t biny = fYaxis.FindBin(y);
      return  binx + nx*biny;
   }
   if (GetDimension() < 4) {
      Int_t nx   = fXaxis.GetNbins()+2;
      Int_t ny   = fYaxis.GetNbins()+2;
      Int_t binx = fXaxis.FindBin(x);
      Int_t biny = fYaxis.FindBin(y);
      Int_t binz = fZaxis.FindBin(z);
      return  binx + nx*(biny +ny*binz);
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return Global bin number corresponding to x,y,z.
///
/// 2-D and 3-D histograms are represented with a one dimensional
/// structure. This has the advantage that all existing functions, such as
/// GetBinContent, GetBinError, GetBinFunction work for all dimensions.
/// This function DOES NOT try to extend the axis if the given point belongs
/// to an under-/overflow bin.
///
/// See also TH1::GetBin, TAxis::FindBin and TAxis::FindFixBin

Int_t TH1::FindFixBin(Double_t x, Double_t y, Double_t z) const
{
   if (GetDimension() < 2) {
      return fXaxis.FindFixBin(x);
   }
   if (GetDimension() < 3) {
      Int_t nx   = fXaxis.GetNbins()+2;
      Int_t binx = fXaxis.FindFixBin(x);
      Int_t biny = fYaxis.FindFixBin(y);
      return  binx + nx*biny;
   }
   if (GetDimension() < 4) {
      Int_t nx   = fXaxis.GetNbins()+2;
      Int_t ny   = fYaxis.GetNbins()+2;
      Int_t binx = fXaxis.FindFixBin(x);
      Int_t biny = fYaxis.FindFixBin(y);
      Int_t binz = fZaxis.FindFixBin(z);
      return  binx + nx*(biny +ny*binz);
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Find first bin with content > threshold for axis (1=x, 2=y, 3=z)
/// if no bins with content > threshold is found the function returns -1.
/// The search will occur between the specified first and last bin. Specifying
/// the value of the last bin to search to less than zero will search until the
/// last defined bin.

Int_t TH1::FindFirstBinAbove(Double_t threshold, Int_t axis, Int_t firstBin, Int_t lastBin) const
{
   if (fBuffer) ((TH1*)this)->BufferEmpty();

   if (axis < 1 || (axis > 1 && GetDimension() == 1 ) ||
       ( axis > 2 && GetDimension() == 2 ) || ( axis > 3 && GetDimension() > 3 ) ) {
      Warning("FindFirstBinAbove","Invalid axis number : %d, axis x assumed\n",axis);
      axis = 1;
   }
   if (firstBin < 1) {
      firstBin = 1;
   }
   Int_t nbinsx = fXaxis.GetNbins();
   Int_t nbinsy = (GetDimension() > 1 ) ? fYaxis.GetNbins() : 1;
   Int_t nbinsz = (GetDimension() > 2 ) ? fZaxis.GetNbins() : 1;

   if (axis == 1) {
      if (lastBin < 0 || lastBin > fXaxis.GetNbins()) {
         lastBin = fXaxis.GetNbins();
      }
      for (Int_t binx = firstBin; binx <= lastBin; binx++) {
         for (Int_t biny = 1; biny <= nbinsy; biny++) {
            for (Int_t binz = 1; binz <= nbinsz; binz++) {
               if (RetrieveBinContent(GetBin(binx,biny,binz)) > threshold) return binx;
            }
         }
      }
   }
   else if (axis == 2) {
      if (lastBin < 0 || lastBin > fYaxis.GetNbins()) {
         lastBin = fYaxis.GetNbins();
      }
      for (Int_t biny = firstBin; biny <= lastBin; biny++) {
         for (Int_t binx = 1; binx <= nbinsx; binx++) {
            for (Int_t binz = 1; binz <= nbinsz; binz++) {
               if (RetrieveBinContent(GetBin(binx,biny,binz)) > threshold) return biny;
           }
         }
      }
   }
   else if (axis == 3) {
      if (lastBin < 0 || lastBin > fZaxis.GetNbins()) {
         lastBin = fZaxis.GetNbins();
      }
      for (Int_t binz = firstBin; binz <= lastBin; binz++) {
         for (Int_t binx = 1; binx <= nbinsx; binx++) {
            for (Int_t biny = 1; biny <= nbinsy; biny++) {
               if (RetrieveBinContent(GetBin(binx,biny,binz)) > threshold) return binz;
            }
         }
      }
   }

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Find last bin with content > threshold for axis (1=x, 2=y, 3=z)
/// if no bins with content > threshold is found the function returns -1.
/// The search will occur between the specified first and last bin. Specifying
/// the value of the last bin to search to less than zero will search until the
/// last defined bin.

Int_t TH1::FindLastBinAbove(Double_t threshold, Int_t axis, Int_t firstBin, Int_t lastBin) const
{
   if (fBuffer) ((TH1*)this)->BufferEmpty();


   if (axis < 1 || ( axis > 1 && GetDimension() == 1 ) ||
       ( axis > 2 && GetDimension() == 2 ) || ( axis > 3 && GetDimension() > 3) ) {
      Warning("FindFirstBinAbove","Invalid axis number : %d, axis x assumed\n",axis);
      axis = 1;
   }
   if (firstBin < 1) {
      firstBin = 1;
   }
   Int_t nbinsx = fXaxis.GetNbins();
   Int_t nbinsy = (GetDimension() > 1 ) ? fYaxis.GetNbins() : 1;
   Int_t nbinsz = (GetDimension() > 2 ) ? fZaxis.GetNbins() : 1;

   if (axis == 1) {
      if (lastBin < 0 || lastBin > fXaxis.GetNbins()) {
         lastBin = fXaxis.GetNbins();
      }
      for (Int_t binx = lastBin; binx >= firstBin; binx--) {
         for (Int_t biny = 1; biny <= nbinsy; biny++) {
            for (Int_t binz = 1; binz <= nbinsz; binz++) {
               if (RetrieveBinContent(GetBin(binx, biny, binz)) > threshold) return binx;
            }
         }
      }
   }
   else if (axis == 2) {
      if (lastBin < 0 || lastBin > fYaxis.GetNbins()) {
         lastBin = fYaxis.GetNbins();
      }
      for (Int_t biny = lastBin; biny >= firstBin; biny--) {
         for (Int_t binx = 1; binx <= nbinsx; binx++) {
            for (Int_t binz = 1; binz <= nbinsz; binz++) {
               if (RetrieveBinContent(GetBin(binx, biny, binz)) > threshold) return biny;
            }
         }
      }
   }
   else if (axis == 3) {
      if (lastBin < 0 || lastBin > fZaxis.GetNbins()) {
         lastBin = fZaxis.GetNbins();
      }
      for (Int_t binz = lastBin; binz >= firstBin; binz--) {
         for (Int_t binx = 1; binx <= nbinsx; binx++) {
            for (Int_t biny = 1; biny <= nbinsy; biny++) {
               if (RetrieveBinContent(GetBin(binx, biny, binz)) > threshold) return binz;
            }
         }
      }
   }

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Search object named name in the list of functions.

TObject *TH1::FindObject(const char *name) const
{
   if (fFunctions) return fFunctions->FindObject(name);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Search object obj in the list of functions.

TObject *TH1::FindObject(const TObject *obj) const
{
   if (fFunctions) return fFunctions->FindObject(obj);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Fit histogram with function fname.
///
/// fname is the name of an already predefined function created by TF1 or TF2
/// Predefined functions such as gaus, expo and poln are automatically
/// created by ROOT.
/// fname can also be a formula, accepted by the linear fitter (linear parts divided
/// by "++" sign), for example "x++sin(x)" for fitting "[0]*x+[1]*sin(x)"
///
///  This function finds a pointer to the TF1 object with name fname
///  and calls TH1::Fit(TF1 *f1,...)

TFitResultPtr TH1::Fit(const char *fname ,Option_t *option ,Option_t *goption, Double_t xxmin, Double_t xxmax)
{
   char *linear;
   linear= (char*)strstr(fname, "++");
   Int_t ndim=GetDimension();
   if (linear){
      if (ndim<2){
         TF1 f1(fname, fname, xxmin, xxmax);
         return Fit(&f1,option,goption,xxmin,xxmax);
      }
      else if (ndim<3){
         TF2 f2(fname, fname);
         return Fit(&f2,option,goption,xxmin,xxmax);
      }
      else{
         TF3 f3(fname, fname);
         return Fit(&f3,option,goption,xxmin,xxmax);
      }
   }
   else{
      TF1 * f1 = (TF1*)gROOT->GetFunction(fname);
      if (!f1) { Printf("Unknown function: %s",fname); return -1; }
      return Fit(f1,option,goption,xxmin,xxmax);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fit histogram with function f1.
///
/// \param[in] option fit options is given in parameter option.
///        - "W"  Ignore the bin uncertainties when fitting using the default least square (chi2) method but skip empty bins
///        - "WW" Ignore the bin uncertainties when fitting using the default least square (chi2) method and include also the empty bins
///        - "I"  Use integral of function in bin, normalized by the bin volume,
///          instead of value at bin center
///        -  "L"  Use Loglikelihood method (default is chisquare method)
///        - "WL" Use Loglikelihood method and bin contents are not integer,
///          i.e. histogram is weighted (must have Sumw2() set)
///        -"MULTI" Use Loglikelihood method based on multi-nomial distribution.
///              In this case function must be normalized and one fits only the function shape (a not extended binned
///              likelihood fit)
///        - "P"  Use Pearson chi2 (using expected errors instead of observed errors)
///        - "U"  Use a User specified fitting algorithm (via SetFCN)
///        - "Q"  Quiet mode (minimum printing)
///        - "V"  Verbose mode (default is between Q and V)
///        - "E"  Perform better Errors estimation using Minos technique
///        - "B"  User defined parameter settings are used for predefined functions
///          like "gaus", "expo", "poln", "landau".
///          Use this option when you want to fix one or more parameters for these functions.
///        - "M"  More. Improve fit results.
///          It uses the IMPROVE command of TMinuit (see TMinuit::mnimpr).
///          This algorithm attempts to improve the found local minimum by searching for a
///          better one.
///        - "R"  Use the Range specified in the function range
///        - "N"  Do not store the graphics function, do not draw
///        - "0"  Do not plot the result of the fit. By default the fitted function
///          is drawn unless the option"N" above is specified.
///        - "+"  Add this new fitted function to the list of fitted functions
///          (by default, any previous function is deleted)
///        - "C"  In case of linear fitting, don't calculate the chisquare
///          (saves time)
///        - "F"  If fitting a polN, switch to minuit fitter
///        - "S"  The result of the fit is returned in the TFitResultPtr
///          (see below Access to the Fit Result)
/// \param[in] goption specify a list of graphics options. See TH1::Draw for a complete list of these options.
/// \param[in] xxmin range
/// \param[in] xxmax range
///
/// In order to use the Range option, one must first create a function
/// with the expression to be fitted. For example, if your histogram
/// has a defined range between -4 and 4 and you want to fit a gaussian
/// only in the interval 1 to 3, you can do:
///
/// ~~~ {.cpp}
///      TF1 *f1 = new TF1("f1", "gaus", 1, 3);
///      histo->Fit("f1", "R");
/// ~~~
///
/// ## Setting initial conditions
/// Parameters must be initialized before invoking the Fit function.
/// The setting of the parameter initial values is automatic for the
/// predefined functions : poln, expo, gaus, landau. One can however disable
/// this automatic computation by specifying the option "B".
/// Note that if a predefined function is defined with an argument,
/// eg, gaus(0), expo(1), you must specify the initial values for
/// the parameters.
/// You can specify boundary limits for some or all parameters via
///
/// ~~~ {.cpp}
///      f1->SetParLimits(p_number, parmin, parmax);
/// ~~~
///
/// if parmin>=parmax, the parameter is fixed
/// Note that you are not forced to fix the limits for all parameters.
/// For example, if you fit a function with 6 parameters, you can do:
///
/// ~~~ {.cpp}
///      func->SetParameters(0, 3.1, 1.e-6, -8, 0, 100);
///      func->SetParLimits(3, -10, -4);
///      func->FixParameter(4, 0);
///      func->SetParLimits(5, 1, 1);
/// ~~~
///
/// With this setup, parameters 0->2 can vary freely
/// Parameter 3 has boundaries [-10,-4] with initial value -8
/// Parameter 4 is fixed to 0
/// Parameter 5 is fixed to 100.
/// When the lower limit and upper limit are equal, the parameter is fixed.
/// However to fix a parameter to 0, one must call the FixParameter function.
///
///
/// #### Changing the fitting objective function
///
/// By default a chi square function is used for fitting. When option "L" (or "LL") is used
/// a Poisson likelihood function (see note below) is used.
/// Using option "MULTI" a multinomial likelihood fit is used. In this case the function normalization is not fitted
/// but only the function shape. Therefore the provided function must be normalized.
/// The functions are defined in the header Fit/Chi2Func.h or Fit/PoissonLikelihoodFCN and they
/// are implemented using the routines FitUtil::EvaluateChi2 or FitUtil::EvaluatePoissonLogL in
/// the file math/mathcore/src/FitUtil.cxx.
/// To specify a User defined fitting function, specify option "U" and
/// call the following functions:
///
/// ~~~ {.cpp}
///      TVirtualFitter::Fitter(myhist)->SetFCN(MyFittingFunction)
/// ~~~
///
/// where MyFittingFunction is of type:
///
/// ~~~ {.cpp}
///      extern void MyFittingFunction(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag);
/// ~~~
///
/// #### Chi2 Fits
///
/// By default a chi2 (least-square) fit is performed on the histogram. The so-called modified least-square method
/// is used where the residual for each bin is computed using as error the observed value (the bin error)
///
/// \f[
///      Chi2 = \sum{ \left(\frac{y(i) - f(x(i) | p )}{e(i)} \right)^2 }
/// \f]
///
/// where y(i) is the bin content for each bin i, x(i) is the bin center and e(i) is the bin error (sqrt(y(i) for
/// an un-weighted histogram. Bins with zero errors are excluded from the fit. See also later the note on the treatment
/// of empty bins. When using option "I" the residual is computed not using the function value at the bin center, f
/// (x(i) | p), but the integral of the function in the bin,   Integral{ f(x|p)dx } divided by the bin volume
///
/// #### Likelihood Fits
///
/// When using option "L" a likelihood fit is used instead of the default chi2 square fit.
/// The likelihood is built assuming a Poisson probability density function for each bin.
/// The negative log-likelihood to be minimized is
///
/// \f[
///       NLL = \sum{ log Poisson ( y(i) | f(x(i) | p ) ) }
/// \f]
///
/// The exact likelihood used is the Poisson likelihood described in this paper:
/// S. Baker and R. D. Cousins, “Clarification of the use of chi-square and likelihood functions in fits to histograms,”
/// Nucl. Instrum. Meth. 221 (1984) 437.
///
/// This method can then be used only when the bin content represents counts (i.e. errors are sqrt(N) ).
/// The likelihood method has the advantage of treating correctly bins with low statistics. In case of high
/// statistics/bin the distribution of the bin content becomes a normal distribution and the likelihood and chi2 fit
/// give the same result.
///
/// The likelihood method, although a bit slower, it is therefore the recommended method in case of low
/// bin statistics, where the chi2 method may give incorrect results, in particular when there are
/// several empty bins (see also below).
/// In case of a weighted histogram, it is possible to perform a likelihood fit by using the
/// option "WL". Note a weighted histogram is an histogram which has been filled with weights and it
/// contains the sum of the weight square ( TH1::Sumw2() has been called). The bin error for a weighted
/// histogram is the square root of the sum of the weight square.
///
/// #### Treatment of Empty Bins
///
/// Empty bins, which have the content equal to zero AND error equal to zero,
/// are excluded by default from the chisquare fit, but they are considered in the likelihood fit.
/// since they affect the likelihood if the function value in these bins is not negligible.
/// When using option "WW" these bins will be considered in the chi2 fit with an error of 1.
/// Note that if the histogram is having bins with zero content and non zero-errors they are considered as
/// any other bins in the fit. Instead bins with zero error and non-zero content are excluded in the chi2 fit.
/// A likelihood fit should also not be performed on such an histogram, since we are assuming a wrong pdf for each bin.
/// In general, one should not fit an histogram with non-empty bins and zero errors, apart if all the bins have zero
/// errors. In this case one could use the option "w", which gives a weight=1 for each bin (unweighted least-square
/// fit).
/// Note that in case of histogram with no errors (chi2 fit with option W or W1) the resulting fitted parameter errors
/// are corrected by the obtained chi2 value using this  expression:  errorp *= sqrt(chisquare/(ndf-1))
///
/// #### Fitting a histogram of dimension N with a function of dimension N-1
///
/// It is possible to fit a TH2 with a TF1 or a TH3 with a TF2.
/// In this case the option "Integral" is not allowed and each cell has
/// equal weight. Also in this case th eobtained parameter error are corrected as in the case when the
/// option "W" is used (see above)
///
/// #### Associated functions
///
/// One or more object (typically a TF1*) can be added to the list
/// of functions (fFunctions) associated to each histogram.
/// When TH1::Fit is invoked, the fitted function is added to this list.
/// Given an histogram h, one can retrieve an associated function
/// with:
///
/// ~~~ {.cpp}
///  TF1 *myfunc = h->GetFunction("myfunc");
/// ~~~
///
/// #### Access to the fit result
///
/// The function returns a TFitResultPtr which can hold a  pointer to a TFitResult object.
/// By default the TFitResultPtr contains only the status of the fit which is return by an
/// automatic conversion of the TFitResultPtr to an integer. One can write in this case directly:
///
/// ~~~ {.cpp}
///     Int_t fitStatus =  h->Fit(myFunc)
/// ~~~
///
/// If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart
/// pointer to it. For example one can do:
///
/// ~~~ {.cpp}
///     TFitResultPtr r = h->Fit(myFunc,"S");
///     TMatrixDSym cov = r->GetCovarianceMatrix();  //  to access the covariance matrix
///     Double_t chi2   = r->Chi2(); // to retrieve the fit chi2
///     Double_t par0   = r->Parameter(0); // retrieve the value for the parameter 0
///     Double_t err0   = r->ParError(0); // retrieve the error for the parameter 0
///     r->Print("V");     // print full information of fit including covariance matrix
///     r->Write();        // store the result in a file
/// ~~~
///
/// The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also
/// from the fitted function.
/// If the histogram is made persistent, the list of
/// associated functions is also persistent. Given a pointer (see above)
/// to an associated function myfunc, one can retrieve the function/fit
/// parameters with calls such as:
///
/// ~~~ {.cpp}
///     Double_t chi2 = myfunc->GetChisquare();
///     Double_t par0 = myfunc->GetParameter(0); //value of 1st parameter
///     Double_t err0 = myfunc->GetParError(0);  //error on first parameter
/// ~~~
///
/// #### Access to the fit status
///
/// The status of the fit can be obtained converting the TFitResultPtr to an integer
/// independently if the fit option "S" is used or not:
///
/// ~~~ {.cpp}
///     TFitResultPtr r = h->Fit(myFunc,opt);
///     Int_t fitStatus = r;
/// ~~~
///
/// The fitStatus is 0 if the fit is OK (i.e no error occurred).
/// The value of the fit status code is negative in case of an error not connected with the
/// minimization procedure, for example  when a wrong function is used.
/// Otherwise the return value is the one returned from the minimization procedure.
/// When TMinuit (default case) or Minuit2 are used as minimizer the status returned is :
/// `fitStatus =  migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult`.
/// TMinuit will return 0 (for migrad, minos, hesse or improve) in case of success and 4 in
/// case of error (see the documentation of TMinuit::mnexcm). So for example, for an error
/// only in Minos but not in Migrad a fitStatus of 40 will be returned.
/// Minuit2 will return also 0 in case of success and different values in migrad minos or
/// hesse depending on the error. See in this case the documentation of
/// Minuit2Minimizer::Minimize for the migradResult, Minuit2Minimizer::GetMinosError for the
/// minosResult and Minuit2Minimizer::Hesse for the hesseResult.
/// If other minimizers are used see their specific documentation for the status code returned.
/// For example in the case of Fumili, for the status returned see TFumili::Minimize.
///
/// #### Excluding points
///
/// Use TF1::RejectPoint inside your fitting function to exclude points
/// within a certain range from the fit. Example:
///
/// ~~~ {.cpp}
///     Double_t fline(Double_t *x, Double_t *par)
///     {
///        if (x[0] > 2.5 && x[0] < 3.5) {
///           TF1::RejectPoint();
///           return 0;
///        }
///        return par[0] + par[1]*x[0];
///     }
///
///     void exclude() {
///        TF1 *f1 = new TF1("f1", "[0] +[1]*x +gaus(2)", 0, 5);
///        f1->SetParameters(6, -1,5, 3, 0.2);
///        TH1F *h = new TH1F("h", "background + signal", 100, 0, 5);
///        h->FillRandom("f1", 2000);
///        TF1 *fline = new TF1("fline", fline, 0, 5, 2);
///        fline->SetParameters(2, -1);
///        h->Fit("fline", "l");
///     }
/// ~~~
///
/// #### Warning when using the option "0"
///
/// When selecting the option "0", the fitted function is added to
/// the list of functions of the histogram, but it is not drawn.
/// You can undo what you disabled in the following way:
///
/// ~~~ {.cpp}
///     h.Fit("myFunction", "0"); // fit, store function but do not draw
///     h.Draw(); function is not drawn
///     const Int_t kNotDraw = 1<<9;
///     h.GetFunction("myFunction")->ResetBit(kNotDraw);
///     h.Draw();  // function is visible again
/// ~~~
///
/// #### Access to the Minimizer information during fitting
///
/// This function calls, the ROOT::Fit::FitObject function implemented in HFitImpl.cxx
/// which uses the ROOT::Fit::Fitter class. The Fitter class creates the objective function
/// (e.g. chi2 or likelihood) and uses an implementation of the  Minimizer interface for minimizing
/// the function.
/// The default minimizer is Minuit (class TMinuitMinimizer which calls TMinuit).
/// The default  can be set in the resource file in etc/system.rootrc. For example
///
/// ~~~ {.cpp}
///     Root.Fitter:      Minuit2
/// ~~~
///
/// A different fitter can also be set via ROOT::Math::MinimizerOptions::SetDefaultMinimizer
/// (or TVirtualFitter::SetDefaultFitter).
/// For example ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSLMultiMin","BFGS");
/// will set the usage of the BFGS algorithm of the GSL multi-dimensional minimization
/// (implemented in libMathMore). ROOT::Math::MinimizerOptions can be used also to set other
/// default options, like maximum number of function calls, minimization tolerance or print
/// level. See the documentation of this class.
///
/// For fitting linear functions (containing the "++" sign" and polN functions,
/// the linear fitter is automatically initialized.

TFitResultPtr TH1::Fit(TF1 *f1 ,Option_t *option ,Option_t *goption, Double_t xxmin, Double_t xxmax)
{
   // implementation of Fit method is in file hist/src/HFitImpl.cxx
   Foption_t fitOption;
   ROOT::Fit::FitOptionsMake(ROOT::Fit::kHistogram,option,fitOption);

   // create range and minimizer options with default values
   ROOT::Fit::DataRange range(xxmin,xxmax);
   ROOT::Math::MinimizerOptions minOption;

   // need to empty the buffer before
   // (t.b.d. do a ML unbinned fit with buffer data)
   if (fBuffer) BufferEmpty();

   return ROOT::Fit::FitObject(this, f1 , fitOption , minOption, goption, range);
}

////////////////////////////////////////////////////////////////////////////////
/// Display a panel with all histogram fit options.
///
/// See class TFitPanel for example

void TH1::FitPanel()
{
   if (!gPad)
      gROOT->MakeDefCanvas();

   if (!gPad) {
      Error("FitPanel", "Unable to create a default canvas");
      return;
   }


   // use plugin manager to create instance of TFitEditor
   TPluginHandler *handler = gROOT->GetPluginManager()->FindHandler("TFitEditor");
   if (handler && handler->LoadPlugin() != -1) {
      if (handler->ExecPlugin(2, gPad, this) == 0)
         Error("FitPanel", "Unable to create the FitPanel");
   }
   else
         Error("FitPanel", "Unable to find the FitPanel plug-in");
}

////////////////////////////////////////////////////////////////////////////////
/// Return an histogram containing the asymmetry of this histogram with h2,
/// where the asymmetry is defined as:
///
/// ~~~ {.cpp}
/// Asymmetry = (h1 - h2)/(h1 + h2)  where h1 = this
/// ~~~
///
/// works for 1D, 2D, etc. histograms
/// c2 is an optional argument that gives a relative weight between the two
/// histograms, and dc2 is the error on this weight.  This is useful, for example,
/// when forming an asymmetry between two histograms from 2 different data sets that
/// need to be normalized to each other in some way.  The function calculates
/// the errors assuming Poisson statistics on h1 and h2 (that is, dh = sqrt(h)).
///
/// example:  assuming 'h1' and 'h2' are already filled
///
/// ~~~ {.cpp}
/// h3 = h1->GetAsymmetry(h2)
/// ~~~
///
/// then 'h3' is created and filled with the asymmetry between 'h1' and 'h2';
/// h1 and h2 are left intact.
///
/// Note that it is the user's responsibility to manage the created histogram.
/// The name of the returned histogram will be `Asymmetry_nameOfh1-nameOfh2`
///
/// code proposed by Jason Seely (seely@mit.edu) and adapted by R.Brun
///
/// clone the histograms so top and bottom will have the
/// correct dimensions:
/// Sumw2 just makes sure the errors will be computed properly
/// when we form sums and ratios below.

TH1 *TH1::GetAsymmetry(TH1* h2, Double_t c2, Double_t dc2)
{
   TH1 *h1 = this;
   TString name =  TString::Format("Asymmetry_%s-%s",h1->GetName(),h2->GetName() );
   TH1 *asym   = (TH1*)Clone(name);

   // set also the title
   TString title = TString::Format("(%s - %s)/(%s+%s)",h1->GetName(),h2->GetName(),h1->GetName(),h2->GetName() );
   asym->SetTitle(title);

   asym->Sumw2();
   Bool_t addStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   TH1 *top    = (TH1*)asym->Clone();
   TH1 *bottom = (TH1*)asym->Clone();
   TH1::AddDirectory(addStatus);

   // form the top and bottom of the asymmetry, and then divide:
   top->Add(h1,h2,1,-c2);
   bottom->Add(h1,h2,1,c2);
   asym->Divide(top,bottom);

   Int_t   xmax = asym->GetNbinsX();
   Int_t   ymax = asym->GetNbinsY();
   Int_t   zmax = asym->GetNbinsZ();

   if (h1->fBuffer) h1->BufferEmpty(1);
   if (h2->fBuffer) h2->BufferEmpty(1);
   if (bottom->fBuffer) bottom->BufferEmpty(1);

   // now loop over bins to calculate the correct errors
   // the reason this error calculation looks complex is because of c2
   for(Int_t i=1; i<= xmax; i++){
      for(Int_t j=1; j<= ymax; j++){
         for(Int_t k=1; k<= zmax; k++){
            Int_t bin = GetBin(i, j, k);
            // here some bin contents are written into variables to make the error
            // calculation a little more legible:
            Double_t a   = h1->RetrieveBinContent(bin);
            Double_t b   = h2->RetrieveBinContent(bin);
            Double_t bot = bottom->RetrieveBinContent(bin);

            // make sure there are some events, if not, then the errors are set = 0
            // automatically.
            //if(bot < 1){} was changed to the next line from recommendation of Jason Seely (28 Nov 2005)
            if(bot < 1e-6){}
            else{
               // computation of errors by Christos Leonidopoulos
               Double_t dasq  = h1->GetBinErrorSqUnchecked(bin);
               Double_t dbsq  = h2->GetBinErrorSqUnchecked(bin);
               Double_t error = 2*TMath::Sqrt(a*a*c2*c2*dbsq + c2*c2*b*b*dasq+a*a*b*b*dc2*dc2)/(bot*bot);
               asym->SetBinError(i,j,k,error);
            }
         }
      }
   }
   delete top;
   delete bottom;

   return asym;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function
/// return the default buffer size for automatic histograms
/// the parameter fgBufferSize may be changed via SetDefaultBufferSize

Int_t TH1::GetDefaultBufferSize()
{
   return fgBufferSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if TH1::Sumw2 must be called when creating new histograms.
/// see TH1::SetDefaultSumw2.

Bool_t TH1::GetDefaultSumw2()
{
   return fgDefaultSumw2;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current number of entries.

Double_t TH1::GetEntries() const
{
   if (fBuffer) {
      Int_t nentries = (Int_t) fBuffer[0];
      if (nentries > 0) return nentries;
   }

   return fEntries;
}

////////////////////////////////////////////////////////////////////////////////
/// Number of effective entries of the histogram.
///
/// \f[
/// neff = \frac{(\sum Weights )^2}{(\sum Weight^2 )}
/// \f]
///
/// In case of an unweighted histogram this number is equivalent to the
/// number of entries of the histogram.
/// For a weighted histogram, this number corresponds to the hypothetical number of unweighted entries
/// a histogram would need to have the same statistical power as this weighted histogram.
/// Note: The underflow/overflow are included if one has set the TH1::StatOverFlows flag
/// and if the statistics has been computed at filling time.
/// If a range is set in the histogram the number is computed from the given range.

Double_t TH1::GetEffectiveEntries() const
{
   Stat_t s[kNstat];
   this->GetStats(s);// s[1] sum of squares of weights, s[0] sum of weights
   return (s[1] ? s[0]*s[0]/s[1] : TMath::Abs(s[0]) );
}

////////////////////////////////////////////////////////////////////////////////
/// Set highlight (enable/disable) mode for the histogram
/// by default highlight mode is disable

void TH1::SetHighlight(Bool_t set)
{
   if (IsHighlight() == set) return;
   if (fDimension > 2) {
      Info("SetHighlight", "Supported only 1-D or 2-D histograms");
      return;
   }

   if (!fPainter) {
      Info("SetHighlight", "Need to draw histogram first");
      return;
   }
   SetBit(kIsHighlight, set);
   fPainter->SetHighlight();
}

////////////////////////////////////////////////////////////////////////////////
/// Redefines TObject::GetObjectInfo.
/// Displays the histogram info (bin number, contents, integral up to bin
/// corresponding to cursor position px,py

char *TH1::GetObjectInfo(Int_t px, Int_t py) const
{
   return ((TH1*)this)->GetPainter()->GetObjectInfo(px,py);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to painter.
/// If painter does not exist, it is created

TVirtualHistPainter *TH1::GetPainter(Option_t *option)
{
   if (!fPainter) {
      TString opt = option;
      opt.ToLower();
      if (opt.Contains("gl") || gStyle->GetCanvasPreferGL()) {
         //try to create TGLHistPainter
         TPluginHandler *handler = gROOT->GetPluginManager()->FindHandler("TGLHistPainter");

         if (handler && handler->LoadPlugin() != -1)
            fPainter = reinterpret_cast<TVirtualHistPainter *>(handler->ExecPlugin(1, this));
      }
   }

   if (!fPainter) fPainter = TVirtualHistPainter::HistPainter(this);

   return fPainter;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute Quantiles for this histogram
/// Quantile x_q of a probability distribution Function F is defined as
///
/// ~~~ {.cpp}
///      F(x_q) = q with 0 <= q <= 1.
/// ~~~
///
/// For instance the median x_0.5 of a distribution is defined as that value
/// of the random variable for which the distribution function equals 0.5:
///
/// ~~~ {.cpp}
///      F(x_0.5) = Probability(x < x_0.5) = 0.5
/// ~~~
///
/// code from Eddy Offermann, Renaissance
///
/// \param[in] nprobSum maximum size of array q and size of array probSum (if given)
/// \param[in] probSum array of positions where quantiles will be computed.
///   - if probSum is null, probSum will be computed internally and will
///     have a size = number of bins + 1 in h. it will correspond to the
///     quantiles calculated at the lowest edge of the histogram (quantile=0) and
///     all the upper edges of the bins.
///   - if probSum is not null, it is assumed to contain at least nprobSum values.
/// \param[out] q array q filled with nq quantiles
/// \return value nq (<=nprobSum) with the number of quantiles computed
///
/// Note that the Integral of the histogram is automatically recomputed
/// if the number of entries is different of the number of entries when
/// the integral was computed last time. In case you do not use the Fill
/// functions to fill your histogram, but SetBinContent, you must call
/// TH1::ComputeIntegral before calling this function.
///
/// Getting quantiles q from two histograms and storing results in a TGraph,
/// a so-called QQ-plot
///
/// ~~~ {.cpp}
/// TGraph *gr = new TGraph(nprob);
/// h1->GetQuantiles(nprob,gr->GetX());
/// h2->GetQuantiles(nprob,gr->GetY());
/// gr->Draw("alp");
/// ~~~
///
/// Example:
///
/// ~~~ {.cpp}
/// void quantiles() {
///    // demo for quantiles
///    const Int_t nq = 20;
///    TH1F *h = new TH1F("h","demo quantiles",100,-3,3);
///    h->FillRandom("gaus",5000);
///
///    Double_t xq[nq];  // position where to compute the quantiles in [0,1]
///    Double_t yq[nq];  // array to contain the quantiles
///    for (Int_t i=0;i<nq;i++) xq[i] = Float_t(i+1)/nq;
///    h->GetQuantiles(nq,yq,xq);
///
///    //show the original histogram in the top pad
///    TCanvas *c1 = new TCanvas("c1","demo quantiles",10,10,700,900);
///    c1->Divide(1,2);
///    c1->cd(1);
///    h->Draw();
///
///    // show the quantiles in the bottom pad
///    c1->cd(2);
///    gPad->SetGrid();
///    TGraph *gr = new TGraph(nq,xq,yq);
///    gr->SetMarkerStyle(21);
///    gr->Draw("alp");
/// }
/// ~~~

Int_t TH1::GetQuantiles(Int_t nprobSum, Double_t *q, const Double_t *probSum)
{
   if (GetDimension() > 1) {
      Error("GetQuantiles","Only available for 1-d histograms");
      return 0;
   }

   const Int_t nbins = GetXaxis()->GetNbins();
   if (!fIntegral) ComputeIntegral();
   if (fIntegral[nbins+1] != fEntries) ComputeIntegral();

   Int_t i, ibin;
   Double_t *prob = (Double_t*)probSum;
   Int_t nq = nprobSum;
   if (probSum == 0) {
      nq = nbins+1;
      prob = new Double_t[nq];
      prob[0] = 0;
      for (i=1;i<nq;i++) {
         prob[i] = fIntegral[i]/fIntegral[nbins];
      }
   }

   for (i = 0; i < nq; i++) {
      ibin = TMath::BinarySearch(nbins,fIntegral,prob[i]);
      while (ibin < nbins-1 && fIntegral[ibin+1] == prob[i]) {
         if (fIntegral[ibin+2] == prob[i]) ibin++;
         else break;
      }
      q[i] = GetBinLowEdge(ibin+1);
      const Double_t dint = fIntegral[ibin+1]-fIntegral[ibin];
      if (dint > 0) q[i] += GetBinWidth(ibin+1)*(prob[i]-fIntegral[ibin])/dint;
   }

   if (!probSum) delete [] prob;
   return nq;
}

////////////////////////////////////////////////////////////////////////////////
/// Decode string choptin and fill fitOption structure.

Int_t TH1::FitOptionsMake(Option_t *choptin, Foption_t &fitOption)
{
   ROOT::Fit::FitOptionsMake(ROOT::Fit::kHistogram, choptin,fitOption);
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute Initial values of parameters for a gaussian.

void H1InitGaus()
{
   Double_t allcha, sumx, sumx2, x, val, stddev, mean;
   Int_t bin;
   const Double_t sqrtpi = 2.506628;

   //   - Compute mean value and StdDev of the histogram in the given range
   TVirtualFitter *hFitter = TVirtualFitter::GetFitter();
   TH1 *curHist = (TH1*)hFitter->GetObjectFit();
   Int_t hxfirst = hFitter->GetXfirst();
   Int_t hxlast  = hFitter->GetXlast();
   Double_t valmax  = curHist->GetBinContent(hxfirst);
   Double_t binwidx = curHist->GetBinWidth(hxfirst);
   allcha = sumx = sumx2 = 0;
   for (bin=hxfirst;bin<=hxlast;bin++) {
      x       = curHist->GetBinCenter(bin);
      val     = TMath::Abs(curHist->GetBinContent(bin));
      if (val > valmax) valmax = val;
      sumx   += val*x;
      sumx2  += val*x*x;
      allcha += val;
   }
   if (allcha == 0) return;
   mean = sumx/allcha;
   stddev  = sumx2/allcha - mean*mean;
   if (stddev > 0) stddev  = TMath::Sqrt(stddev);
   else         stddev  = 0;
   if (stddev == 0) stddev = binwidx*(hxlast-hxfirst+1)/4;
   //if the distribution is really gaussian, the best approximation
   //is binwidx*allcha/(sqrtpi*stddev)
   //However, in case of non-gaussian tails, this underestimates
   //the normalisation constant. In this case the maximum value
   //is a better approximation.
   //We take the average of both quantities
   Double_t constant = 0.5*(valmax+binwidx*allcha/(sqrtpi*stddev));

   //In case the mean value is outside the histo limits and
   //the StdDev is bigger than the range, we take
   //  mean = center of bins
   //  stddev  = half range
   Double_t xmin = curHist->GetXaxis()->GetXmin();
   Double_t xmax = curHist->GetXaxis()->GetXmax();
   if ((mean < xmin || mean > xmax) && stddev > (xmax-xmin)) {
      mean = 0.5*(xmax+xmin);
      stddev  = 0.5*(xmax-xmin);
   }
   TF1 *f1 = (TF1*)hFitter->GetUserFunc();
   f1->SetParameter(0,constant);
   f1->SetParameter(1,mean);
   f1->SetParameter(2,stddev);
   f1->SetParLimits(2,0,10*stddev);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute Initial values of parameters for an exponential.

void H1InitExpo()
{
   Double_t constant, slope;
   Int_t ifail;
   TVirtualFitter *hFitter = TVirtualFitter::GetFitter();
   Int_t hxfirst = hFitter->GetXfirst();
   Int_t hxlast  = hFitter->GetXlast();
   Int_t nchanx  = hxlast - hxfirst + 1;

   H1LeastSquareLinearFit(-nchanx, constant, slope, ifail);

   TF1 *f1 = (TF1*)hFitter->GetUserFunc();
   f1->SetParameter(0,constant);
   f1->SetParameter(1,slope);

}

////////////////////////////////////////////////////////////////////////////////
/// Compute Initial values of parameters for a polynom.

void H1InitPolynom()
{
   Double_t fitpar[25];

   TVirtualFitter *hFitter = TVirtualFitter::GetFitter();
   TF1 *f1 = (TF1*)hFitter->GetUserFunc();
   Int_t hxfirst = hFitter->GetXfirst();
   Int_t hxlast  = hFitter->GetXlast();
   Int_t nchanx  = hxlast - hxfirst + 1;
   Int_t npar    = f1->GetNpar();

   if (nchanx <=1 || npar == 1) {
      TH1 *curHist = (TH1*)hFitter->GetObjectFit();
      fitpar[0] = curHist->GetSumOfWeights()/Double_t(nchanx);
   } else {
      H1LeastSquareFit( nchanx, npar, fitpar);
   }
   for (Int_t i=0;i<npar;i++) f1->SetParameter(i, fitpar[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Least squares lpolynomial fitting without weights.
///
/// \param[in] n number of points to fit
/// \param[in] m number of parameters
/// \param[in] a array of parameters
///
///  based on CERNLIB routine LSQ: Translated to C++ by Rene Brun
///  (E.Keil.  revised by B.Schorr, 23.10.1981.)

void H1LeastSquareFit(Int_t n, Int_t m, Double_t *a)
{
   const Double_t zero = 0.;
   const Double_t one = 1.;
   const Int_t idim = 20;

   Double_t  b[400]        /* was [20][20] */;
   Int_t i, k, l, ifail;
   Double_t power;
   Double_t da[20], xk, yk;

   if (m <= 2) {
      H1LeastSquareLinearFit(n, a[0], a[1], ifail);
      return;
   }
   if (m > idim || m > n) return;
   b[0]  = Double_t(n);
   da[0] = zero;
   for (l = 2; l <= m; ++l) {
      b[l-1]           = zero;
      b[m + l*20 - 21] = zero;
      da[l-1]          = zero;
   }
   TVirtualFitter *hFitter = TVirtualFitter::GetFitter();
   TH1 *curHist  = (TH1*)hFitter->GetObjectFit();
   Int_t hxfirst = hFitter->GetXfirst();
   Int_t hxlast  = hFitter->GetXlast();
   for (k = hxfirst; k <= hxlast; ++k) {
      xk     = curHist->GetBinCenter(k);
      yk     = curHist->GetBinContent(k);
      power  = one;
      da[0] += yk;
      for (l = 2; l <= m; ++l) {
         power   *= xk;
         b[l-1]  += power;
         da[l-1] += power*yk;
      }
      for (l = 2; l <= m; ++l) {
         power            *= xk;
         b[m + l*20 - 21] += power;
      }
   }
   for (i = 3; i <= m; ++i) {
      for (k = i; k <= m; ++k) {
         b[k - 1 + (i-1)*20 - 21] = b[k + (i-2)*20 - 21];
      }
   }
   H1LeastSquareSeqnd(m, b, idim, ifail, 1, da);

   for (i=0; i<m; ++i) a[i] = da[i];

}

////////////////////////////////////////////////////////////////////////////////
/// Least square linear fit without weights.
///
///  extracted from CERNLIB LLSQ: Translated to C++ by Rene Brun
///  (added to LSQ by B. Schorr, 15.02.1982.)

void H1LeastSquareLinearFit(Int_t ndata, Double_t &a0, Double_t &a1, Int_t &ifail)
{
   Double_t xbar, ybar, x2bar;
   Int_t i, n;
   Double_t xybar;
   Double_t fn, xk, yk;
   Double_t det;

   n     = TMath::Abs(ndata);
   ifail = -2;
   xbar  = ybar  = x2bar = xybar = 0;
   TVirtualFitter *hFitter = TVirtualFitter::GetFitter();
   TH1 *curHist  = (TH1*)hFitter->GetObjectFit();
   Int_t hxfirst = hFitter->GetXfirst();
   Int_t hxlast  = hFitter->GetXlast();
   for (i = hxfirst; i <= hxlast; ++i) {
      xk = curHist->GetBinCenter(i);
      yk = curHist->GetBinContent(i);
      if (ndata < 0) {
         if (yk <= 0) yk = 1e-9;
         yk = TMath::Log(yk);
      }
      xbar  += xk;
      ybar  += yk;
      x2bar += xk*xk;
      xybar += xk*yk;
   }
   fn    = Double_t(n);
   det   = fn*x2bar - xbar*xbar;
   ifail = -1;
   if (det <= 0) {
      a0 = ybar/fn;
      a1 = 0;
      return;
   }
   ifail = 0;
   a0 = (x2bar*ybar - xbar*xybar) / det;
   a1 = (fn*xybar - xbar*ybar) / det;

}

////////////////////////////////////////////////////////////////////////////////
/// Extracted from CERN Program library routine DSEQN.
///
/// Translated to C++ by Rene Brun

void H1LeastSquareSeqnd(Int_t n, Double_t *a, Int_t idim, Int_t &ifail, Int_t k, Double_t *b)
{
   Int_t a_dim1, a_offset, b_dim1, b_offset;
   Int_t nmjp1, i, j, l;
   Int_t im1, jp1, nm1, nmi;
   Double_t s1, s21, s22;
   const Double_t one = 1.;

   /* Parameter adjustments */
   b_dim1 = idim;
   b_offset = b_dim1 + 1;
   b -= b_offset;
   a_dim1 = idim;
   a_offset = a_dim1 + 1;
   a -= a_offset;

   if (idim < n) return;

   ifail = 0;
   for (j = 1; j <= n; ++j) {
      if (a[j + j*a_dim1] <= 0) { ifail = -1; return; }
      a[j + j*a_dim1] = one / a[j + j*a_dim1];
      if (j == n) continue;
      jp1 = j + 1;
      for (l = jp1; l <= n; ++l) {
         a[j + l*a_dim1] = a[j + j*a_dim1] * a[l + j*a_dim1];
         s1 = -a[l + (j+1)*a_dim1];
         for (i = 1; i <= j; ++i) { s1 = a[l + i*a_dim1] * a[i + (j+1)*a_dim1] + s1; }
         a[l + (j+1)*a_dim1] = -s1;
      }
   }
   if (k <= 0) return;

   for (l = 1; l <= k; ++l) {
      b[l*b_dim1 + 1] = a[a_dim1 + 1]*b[l*b_dim1 + 1];
   }
   if (n == 1) return;
   for (l = 1; l <= k; ++l) {
      for (i = 2; i <= n; ++i) {
         im1 = i - 1;
         s21 = -b[i + l*b_dim1];
         for (j = 1; j <= im1; ++j) {
            s21 = a[i + j*a_dim1]*b[j + l*b_dim1] + s21;
         }
         b[i + l*b_dim1] = -a[i + i*a_dim1]*s21;
      }
      nm1 = n - 1;
      for (i = 1; i <= nm1; ++i) {
         nmi = n - i;
         s22 = -b[nmi + l*b_dim1];
         for (j = 1; j <= i; ++j) {
            nmjp1 = n - j + 1;
            s22 = a[nmi + nmjp1*a_dim1]*b[nmjp1 + l*b_dim1] + s22;
         }
         b[nmi + l*b_dim1] = -s22;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return Global bin number corresponding to binx,y,z.
///
/// 2-D and 3-D histograms are represented with a one dimensional
/// structure.
/// This has the advantage that all existing functions, such as
/// GetBinContent, GetBinError, GetBinFunction work for all dimensions.
///
/// In case of a TH1x, returns binx directly.
/// see TH1::GetBinXYZ for the inverse transformation.
///
/// Convention for numbering bins
///
/// For all histogram types: nbins, xlow, xup
///
///  - bin = 0;       underflow bin
///  - bin = 1;       first bin with low-edge xlow INCLUDED
///  - bin = nbins;   last bin with upper-edge xup EXCLUDED
///  - bin = nbins+1; overflow bin
///
/// In case of 2-D or 3-D histograms, a "global bin" number is defined.
/// For example, assuming a 3-D histogram with binx,biny,binz, the function
///
/// ~~~ {.cpp}
///     Int_t bin = h->GetBin(binx,biny,binz);
/// ~~~
///
/// returns a global/linearized bin number. This global bin is useful
/// to access the bin information independently of the dimension.

Int_t TH1::GetBin(Int_t binx, Int_t, Int_t) const
{
   Int_t ofx = fXaxis.GetNbins() + 1; // overflow bin
   if (binx < 0) binx = 0;
   if (binx > ofx) binx = ofx;

   return binx;
}

////////////////////////////////////////////////////////////////////////////////
/// Return binx, biny, binz corresponding to the global bin number globalbin
/// see TH1::GetBin function above

void TH1::GetBinXYZ(Int_t binglobal, Int_t &binx, Int_t &biny, Int_t &binz) const
{
   Int_t nx  = fXaxis.GetNbins()+2;
   Int_t ny  = fYaxis.GetNbins()+2;

   if (GetDimension() == 1) {
      binx = binglobal%nx;
      biny = 0;
      binz = 0;
      return;
   }
   if (GetDimension() == 2) {
      binx = binglobal%nx;
      biny = ((binglobal-binx)/nx)%ny;
      binz = 0;
      return;
   }
   if (GetDimension() == 3) {
      binx = binglobal%nx;
      biny = ((binglobal-binx)/nx)%ny;
      binz = ((binglobal-binx)/nx -biny)/ny;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return a random number distributed according the histogram bin contents.
/// This function checks if the bins integral exists. If not, the integral
/// is evaluated, normalized to one.
///
/// @param rng (optional) Random number generator pointer used (default is gRandom)
///
/// The integral is automatically recomputed if the number of entries
/// is not the same then when the integral was computed.
/// NB Only valid for 1-d histograms. Use GetRandom2 or 3 otherwise.
/// If the histogram has a bin with negative content a NaN is returned

Double_t TH1::GetRandom(TRandom * rng) const
{
   if (fDimension > 1) {
      Error("GetRandom","Function only valid for 1-d histograms");
      return 0;
   }
   Int_t nbinsx = GetNbinsX();
   Double_t integral = 0;
   // compute integral checking that all bins have positive content (see ROOT-5894)
   if (fIntegral) {
      if (fIntegral[nbinsx+1] != fEntries) integral = ((TH1*)this)->ComputeIntegral(true);
      else  integral = fIntegral[nbinsx];
   } else {
      integral = ((TH1*)this)->ComputeIntegral(true);
   }
   if (integral == 0) return 0;
   // return a NaN in case some bins have negative content
   if (integral == TMath::QuietNaN() ) return TMath::QuietNaN();

   Double_t r1 = (rng) ? rng->Rndm() : gRandom->Rndm();
   Int_t ibin = TMath::BinarySearch(nbinsx,fIntegral,r1);
   Double_t x = GetBinLowEdge(ibin+1);
   if (r1 > fIntegral[ibin]) x +=
      GetBinWidth(ibin+1)*(r1-fIntegral[ibin])/(fIntegral[ibin+1] - fIntegral[ibin]);
   return x;
}

////////////////////////////////////////////////////////////////////////////////
/// Return content of bin number bin.
///
/// Implemented in TH1C,S,F,D
///
///  Convention for numbering bins
///
///  For all histogram types: nbins, xlow, xup
///
///  - bin = 0;       underflow bin
///  - bin = 1;       first bin with low-edge xlow INCLUDED
///  - bin = nbins;   last bin with upper-edge xup EXCLUDED
///  - bin = nbins+1; overflow bin
///
///  In case of 2-D or 3-D histograms, a "global bin" number is defined.
///  For example, assuming a 3-D histogram with binx,biny,binz, the function
///
/// ~~~ {.cpp}
///    Int_t bin = h->GetBin(binx,biny,binz);
/// ~~~
///
///  returns a global/linearized bin number. This global bin is useful
///  to access the bin information independently of the dimension.

Double_t TH1::GetBinContent(Int_t bin) const
{
   if (fBuffer) const_cast<TH1*>(this)->BufferEmpty();
   if (bin < 0) bin = 0;
   if (bin >= fNcells) bin = fNcells-1;

   return RetrieveBinContent(bin);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute first binx in the range [firstx,lastx] for which
/// diff = abs(bin_content-c) <= maxdiff
///
/// In case several bins in the specified range with diff=0 are found
/// the first bin found is returned in binx.
/// In case several bins in the specified range satisfy diff <=maxdiff
/// the bin with the smallest difference is returned in binx.
/// In all cases the function returns the smallest difference.
///
/// NOTE1: if firstx <= 0, firstx is set to bin 1
///    if (lastx < firstx then firstx is set to the number of bins
///    ie if firstx=0 and lastx=0 (default) the search is on all bins.
///
/// NOTE2: if maxdiff=0 (default), the first bin with content=c is returned.

Double_t TH1::GetBinWithContent(Double_t c, Int_t &binx, Int_t firstx, Int_t lastx,Double_t maxdiff) const
{
   if (fDimension > 1) {
      binx = 0;
      Error("GetBinWithContent","function is only valid for 1-D histograms");
      return 0;
   }

   if (fBuffer) ((TH1*)this)->BufferEmpty();

   if (firstx <= 0) firstx = 1;
   if (lastx < firstx) lastx = fXaxis.GetNbins();
   Int_t binminx = 0;
   Double_t diff, curmax = 1.e240;
   for (Int_t i=firstx;i<=lastx;i++) {
      diff = TMath::Abs(RetrieveBinContent(i)-c);
      if (diff <= 0) {binx = i; return diff;}
      if (diff < curmax && diff <= maxdiff) {curmax = diff, binminx=i;}
   }
   binx = binminx;
   return curmax;
}

////////////////////////////////////////////////////////////////////////////////
/// Given a point x, approximates the value via linear interpolation
/// based on the two nearest bin centers
///
/// Andy Mastbaum 10/21/08

Double_t TH1::Interpolate(Double_t x) const
{
   if (fBuffer) ((TH1*)this)->BufferEmpty();

   Int_t xbin = fXaxis.FindFixBin(x);
   Double_t x0,x1,y0,y1;

   if(x<=GetBinCenter(1)) {
      return RetrieveBinContent(1);
   } else if(x>=GetBinCenter(GetNbinsX())) {
      return RetrieveBinContent(GetNbinsX());
   } else {
      if(x<=GetBinCenter(xbin)) {
         y0 = RetrieveBinContent(xbin-1);
         x0 = GetBinCenter(xbin-1);
         y1 = RetrieveBinContent(xbin);
         x1 = GetBinCenter(xbin);
      } else {
         y0 = RetrieveBinContent(xbin);
         x0 = GetBinCenter(xbin);
         y1 = RetrieveBinContent(xbin+1);
         x1 = GetBinCenter(xbin+1);
      }
      return y0 + (x-x0)*((y1-y0)/(x1-x0));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// 2d Interpolation. Not yet implemented.

Double_t TH1::Interpolate(Double_t, Double_t) const
{
   Error("Interpolate","This function must be called with 1 argument for a TH1");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// 3d Interpolation. Not yet implemented.

Double_t TH1::Interpolate(Double_t, Double_t, Double_t) const
{
   Error("Interpolate","This function must be called with 1 argument for a TH1");
   return 0;
}

///////////////////////////////////////////////////////////////////////////////
/// Check if an histogram is empty
///  (this a protected method used mainly by TH1Merger )

Bool_t TH1::IsEmpty() const
{
   // if fTsumw or fentries are not zero histogram is not empty
   // need to use GetEntries() instead of fEntries in case of bugger histograms
   // so we will flash the buffer
   if (fTsumw != 0) return kFALSE;
   if (GetEntries() != 0) return kFALSE;
   // case fTSumw == 0 amd entries are also zero
   // this should not really happening, but if one sets content by hand
   // it can happen. a call to ResetStats() should be done in such cases
   double sumw = 0;
   for (int i = 0; i< GetNcells(); ++i) sumw += RetrieveBinContent(i);
   return (sumw != 0) ? kFALSE : kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the bin is overflow.

Bool_t TH1::IsBinOverflow(Int_t bin, Int_t iaxis) const
{
   Int_t binx, biny, binz;
   GetBinXYZ(bin, binx, biny, binz);

   if (iaxis == 0) {
      if ( fDimension == 1 )
         return binx >= GetNbinsX() + 1;
      if ( fDimension == 2 )
         return (binx >= GetNbinsX() + 1) ||
            (biny >= GetNbinsY() + 1);
      if ( fDimension == 3 )
         return (binx >= GetNbinsX() + 1) ||
            (biny >= GetNbinsY() + 1) ||
            (binz >= GetNbinsZ() + 1);
      return kFALSE;
   }
   if (iaxis == 1)
      return binx >= GetNbinsX() + 1;
   if (iaxis == 2)
      return biny >= GetNbinsY() + 1;
   if (iaxis == 3)
      return binz >= GetNbinsZ() + 1;

   Error("IsBinOverflow","Invalid axis value");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the bin is underflow.
/// If iaxis = 0  make OR with all axes otherwise check only for the given axis

Bool_t TH1::IsBinUnderflow(Int_t bin, Int_t iaxis) const
{
   Int_t binx, biny, binz;
   GetBinXYZ(bin, binx, biny, binz);

   if (iaxis == 0) {
      if ( fDimension == 1 )
         return (binx <= 0);
      else if ( fDimension == 2 )
         return (binx <= 0 || biny <= 0);
      else if ( fDimension == 3 )
         return (binx <= 0 || biny <= 0 || binz <= 0);
      else
         return kFALSE;
   }
   if (iaxis == 1)
       return (binx <= 0);
   if (iaxis == 2)
      return (biny <= 0);
   if (iaxis == 3)
      return (binz <= 0);

   Error("IsBinUnderflow","Invalid axis value");
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Reduce the number of bins for the axis passed in the option to the number of bins having a label.
/// The method will remove only the extra bins existing after the last "labeled" bin.
/// Note that if there are "un-labeled" bins present between "labeled" bins they will not be removed

void TH1::LabelsDeflate(Option_t *ax)
{
   Int_t iaxis = AxisChoice(ax);
   TAxis *axis = 0;
   if (iaxis == 1) axis = GetXaxis();
   if (iaxis == 2) axis = GetYaxis();
   if (iaxis == 3) axis = GetZaxis();
   if (!axis) {
      Error("LabelsDeflate","Invalid axis option %s",ax);
      return;
   }
   if (!axis->GetLabels()) return;

   // find bin with last labels
   // bin number is object ID in list of labels
   // therefore max bin number is number of bins of the deflated histograms
   TIter next(axis->GetLabels());
   TObject *obj;
   Int_t nbins = 0;
   while ((obj = next())) {
      Int_t ibin = obj->GetUniqueID();
      if (ibin > nbins) nbins = ibin;
   }
   if (nbins < 1) nbins = 1;

   // Do nothing in case it was the last bin
   if (nbins==axis->GetNbins()) return;

   TH1 *hold = (TH1*)IsA()->New();
   R__ASSERT(hold);
   hold->SetDirectory(0);
   Copy(*hold);

   Bool_t timedisp = axis->GetTimeDisplay();
   Double_t xmin = axis->GetXmin();
   Double_t xmax = axis->GetBinUpEdge(nbins);
   if (xmax <= xmin) xmax = xmin +nbins;
   axis->SetRange(0,0);
   axis->Set(nbins,xmin,xmax);
   SetBinsLength(-1);  // reset the number of cells
   Int_t errors = fSumw2.fN;
   if (errors) fSumw2.Set(fNcells);
   axis->SetTimeDisplay(timedisp);
   // reset histogram content
   Reset("ICE");

   //now loop on all bins and refill
   // NOTE that if the bins without labels have content
   // it will be put in the underflow/overflow.
   // For this reason we use AddBinContent method
   Double_t oldEntries = fEntries;
   Int_t bin,binx,biny,binz;
   for (bin=0; bin < hold->fNcells; ++bin) {
      hold->GetBinXYZ(bin,binx,biny,binz);
      Int_t ibin = GetBin(binx,biny,binz);
      Double_t cu = hold->RetrieveBinContent(bin);
      AddBinContent(ibin,cu);
      if (errors) {
         fSumw2.fArray[ibin] += hold->fSumw2.fArray[bin];
      }
   }
   fEntries = oldEntries;
   delete hold;
}

////////////////////////////////////////////////////////////////////////////////
/// Double the number of bins for axis.
/// Refill histogram
/// This function is called by TAxis::FindBin(const char *label)

void TH1::LabelsInflate(Option_t *ax)
{
   Int_t iaxis = AxisChoice(ax);
   TAxis *axis = 0;
   if (iaxis == 1) axis = GetXaxis();
   if (iaxis == 2) axis = GetYaxis();
   if (iaxis == 3) axis = GetZaxis();
   if (!axis) return;

   TH1 *hold = (TH1*)IsA()->New();;
   hold->SetDirectory(0);
   Copy(*hold);

   Bool_t timedisp = axis->GetTimeDisplay();
   Int_t nbins   = axis->GetNbins();
   Double_t xmin = axis->GetXmin();
   Double_t xmax = axis->GetXmax();
   xmax = xmin + 2*(xmax-xmin);
   axis->SetRange(0,0);
   // double the bins and recompute ncells
   axis->Set(2*nbins,xmin,xmax);
   SetBinsLength(-1);
   Int_t errors = fSumw2.fN;
   if (errors) fSumw2.Set(fNcells);
   axis->SetTimeDisplay(timedisp);

   Reset("ICE");  // reset content and error

   //now loop on all bins and refill
   Double_t oldEntries = fEntries;
   Int_t bin,ibin,binx,biny,binz;
   for (ibin =0; ibin < hold->fNcells; ibin++) {
      // get the binx,y,z values . The x-y-z (axis) bin values will stay the same between new-old after the expanding
      hold->GetBinXYZ(ibin,binx,biny,binz);
      bin = GetBin(binx,biny,binz);

      // underflow and overflow will be cleaned up because their meaning has been altered
      if (hold->IsBinUnderflow(ibin,iaxis) || hold->IsBinOverflow(ibin,iaxis)) {
         continue;
      }
      else {
         AddBinContent(bin, hold->RetrieveBinContent(ibin));
         if (errors) fSumw2.fArray[bin] += hold->fSumw2.fArray[ibin];
      }
   }
   fEntries = oldEntries;
   delete hold;
}

////////////////////////////////////////////////////////////////////////////////
/// Sort bins with labels or set option(s) to draw axis with labels
/// \param[in] option
///     - "a" sort by alphabetic order
///     - ">" sort by decreasing values
///     - "<" sort by increasing values
///     - "h" draw labels horizontal
///     - "v" draw labels vertical
///     - "u" draw labels up (end of label right adjusted)
///     - "d" draw labels down (start of label left adjusted)
///
/// In case not all bins have labels sorting will work only in the case
/// the first `n` consecutive bins have all labels and sorting will be performed on
/// those label bins.
///
/// \param[in] ax axis

void TH1::LabelsOption(Option_t *option, Option_t *ax)
{
   Int_t iaxis = AxisChoice(ax);
   TAxis *axis = 0;
   if (iaxis == 1)
      axis = GetXaxis();
   if (iaxis == 2)
      axis = GetYaxis();
   if (iaxis == 3)
      axis = GetZaxis();
   if (!axis)
      return;
   THashList *labels = axis->GetLabels();
   if (!labels) {
      Warning("LabelsOption", "Axis %s has no labels!",axis->GetName());
      return;
   }
   TString opt = option;
   opt.ToLower();
   Int_t iopt = -1;
   if (opt.Contains("h")) {
      axis->SetBit(TAxis::kLabelsHori);
      axis->ResetBit(TAxis::kLabelsVert);
      axis->ResetBit(TAxis::kLabelsDown);
      axis->ResetBit(TAxis::kLabelsUp);
      iopt = 0;
   }
   if (opt.Contains("v")) {
      axis->SetBit(TAxis::kLabelsVert);
      axis->ResetBit(TAxis::kLabelsHori);
      axis->ResetBit(TAxis::kLabelsDown);
      axis->ResetBit(TAxis::kLabelsUp);
      iopt = 1;
   }
   if (opt.Contains("u")) {
      axis->SetBit(TAxis::kLabelsUp);
      axis->ResetBit(TAxis::kLabelsVert);
      axis->ResetBit(TAxis::kLabelsDown);
      axis->ResetBit(TAxis::kLabelsHori);
      iopt = 2;
   }
   if (opt.Contains("d")) {
      axis->SetBit(TAxis::kLabelsDown);
      axis->ResetBit(TAxis::kLabelsVert);
      axis->ResetBit(TAxis::kLabelsHori);
      axis->ResetBit(TAxis::kLabelsUp);
      iopt = 3;
   }
   Int_t sort = -1;
   if (opt.Contains("a"))
      sort = 0;
   if (opt.Contains(">"))
      sort = 1;
   if (opt.Contains("<"))
      sort = 2;
   if (sort < 0) {
      if (iopt < 0)
         Error("LabelsOption", "%s is an invalid label placement option!",opt.Data());
      return;
   }

   // Code works only if first n bins have labels if we uncomment following line
   // but we don't want to support this special case
   // Int_t n = TMath::Min(axis->GetNbins(), labels->GetSize());

   // support only cases where each bin has a labels (should be when axis is alphanumeric)
   Int_t n = labels->GetSize();
   if (n != axis->GetNbins()) {
      // check if labels are all consecutive and starts from the first bin
      // in that case the current code will work fine
      Int_t firstLabelBin = axis->GetNbins()+1;
      Int_t lastLabelBin = -1;
      for (Int_t i = 0; i < n; ++i) {
         Int_t bin  = labels->At(i)->GetUniqueID();
         if (bin < firstLabelBin) firstLabelBin = bin;
         if (bin > lastLabelBin) lastLabelBin = bin;
      }
      if (firstLabelBin != 1 || lastLabelBin-firstLabelBin +1 != n) {
         Error("LabelsOption", "%s of Histogram %s contains bins without labels. Sorting will not work correctly - return",
            axis->GetName(), GetName());
         return;
      }
      // case where label bins are consecutive starting from first bin will work
      // calling before a TH1::LabelsDeflate() will avoid this error message
      Warning("LabelsOption", "axis %s of Histogram %s has extra following bins without labels. Sorting will work only for first label bins",
            axis->GetName(), GetName());
   }
   std::vector<Int_t> a(n);
   std::vector<Int_t> b(n);


   Int_t i, j, k;
   std::vector<Double_t> cont;
   std::vector<Double_t> errors2;
   THashList *labold = new THashList(labels->GetSize(), 1);
   TIter nextold(labels);
   TObject *obj = nullptr;
   labold->AddAll(labels);
   labels->Clear();

   // delete buffer if it is there since bins will be reordered.
   if (fBuffer)
      BufferEmpty(1);

   if (sort > 0) {
      //---sort by values of bins
      if (GetDimension() == 1) {
         cont.resize(n);
         if (fSumw2.fN)
            errors2.resize(n);
         for (i = 0; i < n; i++) {
            cont[i] = RetrieveBinContent(i + 1);
            if (!errors2.empty())
               errors2[i] = GetBinErrorSqUnchecked(i + 1);
            b[i] = labold->At(i)->GetUniqueID(); // this is the bin corresponding to the label
            a[i] = i;
         }
         if (sort == 1)
            TMath::Sort(n, cont.data(), a.data(), kTRUE); // sort by decreasing values
         else
            TMath::Sort(n, cont.data(), a.data(), kFALSE); // sort by increasing values
         for (i = 0; i < n; i++) {
            // use UpdateBinCOntent to not screw up histogram entries
            UpdateBinContent(i + 1, cont[b[a[i]] - 1]); // b[a[i]] returns bin number. .we need to subtract 1
            if (gDebug)
               Info("LabelsOption","setting bin %d value %f from bin %d label %s at pos %d ",
                         i+1,cont[b[a[i]] - 1],b[a[i]],labold->At(a[i])->GetName(),a[i]);
            if (!errors2.empty())
               fSumw2.fArray[i + 1] =  errors2[b[a[i]] - 1];
         }
         for (i = 0; i < n; i++) {
            obj = labold->At(a[i]);
            labels->Add(obj);
            obj->SetUniqueID(i + 1);
         }
      } else if (GetDimension() == 2) {
         std::vector<Double_t> pcont(n + 2);
         Int_t nx = fXaxis.GetNbins() + 2;
         Int_t ny = fYaxis.GetNbins() + 2;
         cont.resize((nx + 2) * (ny + 2));
         if (fSumw2.fN)
            errors2.resize((nx + 2) * (ny + 2));
         for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
               Int_t bin = GetBin(i,j);
               cont[i + nx * j] = RetrieveBinContent(bin);
               if (!errors2.empty())
                  errors2[i + nx * j] = GetBinErrorSqUnchecked(bin);
               if (axis == GetXaxis())
                  k = i - 1;
               else
                  k = j - 1;
               if (k >= 0 && k < n) { // we consider underflow/overflows in y for ordering the bins
                  pcont[k] += cont[i + nx * j];
                  a[k] = k;
               }
            }
         }
         if (sort == 1)
            TMath::Sort(n, pcont.data(), a.data(), kTRUE); // sort by decreasing values
         else
            TMath::Sort(n, pcont.data(), a.data(), kFALSE); // sort by increasing values
         for (i = 0; i < n; i++) {
            // iterate on old label  list to find corresponding bin match
            TIter next(labold);
            UInt_t bin = a[i] + 1;
            while ((obj = next())) {
               if (obj->GetUniqueID() == (UInt_t)bin)
                  break;
               else
                  obj = nullptr;
            }
            if (!obj) {
               // this should not really happen
               R__ASSERT("LabelsOption - No corresponding bin found when ordering labels");
               return;
            }

            labels->Add(obj);
            if (gDebug)
               std::cout << " set label " << obj->GetName() << " to bin " << i + 1 << " from order " << a[i] << " bin "
                         << b[a[i]] << "content " << pcont[a[i]] << std::endl;
         }
         // need to set here new ordered labels - otherwise loop before does not work since labold and labels list
         // contain same objects
         for (i = 0; i < n; i++) {
            labels->At(i)->SetUniqueID(i + 1);
         }
         // set now the bin contents
         if (axis == GetXaxis()) {
            for (i = 0; i < n; i++) {
               Int_t ix = a[i] + 1;
               for (j = 0; j < ny; j++) {
                  Int_t bin = GetBin(i + 1, j);
                  UpdateBinContent(bin, cont[ix + nx * j]);
                  if (!errors2.empty())
                     fSumw2.fArray[bin] = errors2[ix + nx * j];
               }
            }
         } else {
            // using y axis
            for (i = 0; i < nx; i++) {
               for (j = 0; j < n; j++) {
                  Int_t iy = a[j] + 1;
                  Int_t bin = GetBin(i, j + 1);
                  UpdateBinContent(bin, cont[i + nx * iy]);
                  if (!errors2.empty())
                     fSumw2.fArray[bin] = errors2[i + nx * iy];
               }
            }
         }
      } else {
         // sorting histograms: 3D case
         std::vector<Double_t> pcont(n + 2);
         Int_t nx = fXaxis.GetNbins() + 2;
         Int_t ny = fYaxis.GetNbins() + 2;
         Int_t nz = fZaxis.GetNbins() + 2;
         Int_t l = 0;
         cont.resize((nx + 2) * (ny + 2) * (nz + 2));
         if (fSumw2.fN)
            errors2.resize((nx + 2) * (ny + 2) * (nz + 2));
         for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
               for (k = 0; k < nz; k++) {
                  Int_t bin  = GetBin(i,j,k);
                  Double_t c  = RetrieveBinContent(bin);
                  if (axis == GetXaxis())
                     l = i - 1;
                  else if (axis == GetYaxis())
                     l = j - 1;
                  else
                     l = k - 1;
                  if (l >= 0 && l < n) { // we consider underflow/overflows in y for ordering the bins
                     pcont[l] += c;
                     a[l] = l;
                  }
                  cont[i + nx * (j + ny * k)] = c;
                  if (!errors2.empty())
                     errors2[i + nx * (j + ny * k)] = GetBinErrorSqUnchecked(bin);
               }
            }
         }
         if (sort == 1)
            TMath::Sort(n, pcont.data(), a.data(), kTRUE); // sort by decreasing values
         else
            TMath::Sort(n, pcont.data(), a.data(), kFALSE); // sort by increasing values
         for (i = 0; i < n; i++) {
            // iterate on the old label  list to find corresponding bin match
            TIter next(labold);
            UInt_t bin = a[i] + 1;
            obj = nullptr;
            while ((obj = next())) {
               if (obj->GetUniqueID() == (UInt_t)bin) {
                  break;
               }
               else
                  obj = nullptr;
            }
            if (!obj) {
               R__ASSERT("LabelsOption - No corresponding bin found when ordering labels");
               return;
            }
            labels->Add(obj);
            if (gDebug)
               std::cout << " set label " << obj->GetName() << " to bin " << i + 1 << " from bin " << a[i] << "content "
                         << pcont[a[i]] << std::endl;
         }

         // need to set here new ordered labels - otherwise loop before does not work since labold and llabels list
         // contain same objects
         for (i = 0; i < n; i++) {
            labels->At(i)->SetUniqueID(i + 1);
         }
         // set now the bin contents
         if (axis == GetXaxis()) {
            for (i = 0; i < n; i++) {
               Int_t ix = a[i] + 1;
               for (j = 0; j < ny; j++) {
                  for (k = 0; k < nz; k++) {
                     Int_t bin = GetBin(i + 1, j, k);
                     UpdateBinContent(bin, cont[ix + nx * (j + ny * k)]);
                     if (!errors2.empty())
                        fSumw2.fArray[bin] = errors2[ix + nx * (j + ny * k)];
                  }
               }
            }
         } else if (axis == GetYaxis()) {
            // using y axis
            for (i = 0; i < nx; i++) {
               for (j = 0; j < n; j++) {
                  Int_t iy = a[j] + 1;
                  for (k = 0; k < nz; k++) {
                     Int_t bin = GetBin(i, j + 1, k);
                     UpdateBinContent(bin, cont[i + nx * (iy + ny * k)]);
                     if (!errors2.empty())
                       fSumw2.fArray[bin] = errors2[i + nx * (iy + ny * k)];
                  }
               }
            }
         } else {
            // using z axis
            for (i = 0; i < nx; i++) {
               for (j = 0; j < ny; j++) {
                  for (k = 0; k < n; k++) {
                     Int_t iz = a[k] + 1;
                     Int_t bin = GetBin(i, j , k +1);
                     UpdateBinContent(bin, cont[i + nx * (j + ny * iz)]);
                     if (!errors2.empty())
                         fSumw2.fArray[bin] = errors2[i + nx * (j + ny * iz)];
                  }
               }
            }
         }
      }
   } else {
      //---alphabetic sort
      // sort labels using vector of strings and TMath::Sort
      // I need to array because labels order in list is not necessary that of the bins
      std::vector<std::string> vecLabels(n);
      for (i = 0; i < n; i++) {
         vecLabels[i] = labold->At(i)->GetName();
         b[i] = labold->At(i)->GetUniqueID(); // this is the bin corresponding to the label
         a[i] = i;
      }
      // sort in ascending order for strings
      TMath::Sort(n, vecLabels.data(), a.data(), kFALSE);
      // set the new labels
      for (i = 0; i < n; i++) {
         TObject *labelObj = labold->At(a[i]);
         labels->Add(labold->At(a[i]));
         // set the corresponding bin. NB bin starts from 1
         labelObj->SetUniqueID(i + 1);
         if (gDebug)
            std::cout << "bin " << i + 1 << " setting new labels for axis " << labold->At(a[i])->GetName() << " from "
                      << b[a[i]] << std::endl;
      }

      if (GetDimension() == 1) {
         cont.resize(n + 2);
         if (fSumw2.fN)
            errors2.resize(n + 2);
         for (i = 0; i < n; i++) {
            cont[i] = RetrieveBinContent(b[a[i]]);
            if (!errors2.empty())
               errors2[i] = GetBinErrorSqUnchecked(b[a[i]]);
         }
         for (i = 0; i < n; i++) {
            UpdateBinContent(i + 1, cont[i]);
            if (!errors2.empty())
               fSumw2.fArray[i+1] = errors2[i];
         }
      } else if (GetDimension() == 2) {
         Int_t nx = fXaxis.GetNbins() + 2;
         Int_t ny = fYaxis.GetNbins() + 2;
         cont.resize(nx * ny);
         if (fSumw2.fN)
            errors2.resize(nx * ny);
         // copy old bin contents and then set to new ordered bins
         // N.B. bin in histograms starts from 1, but in y we consider under/overflows
         for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) { // ny is nbins+2
               Int_t bin = GetBin(i, j);
               cont[i + nx * j] = RetrieveBinContent(bin);
               if (!errors2.empty())
                  errors2[i + nx * j] = GetBinErrorSqUnchecked(bin);
            }
         }
         if (axis == GetXaxis()) {
            for (i = 0; i < n; i++) {
               for (j = 0; j < ny; j++) {
                  Int_t bin = GetBin(i + 1 , j);
                  UpdateBinContent(bin, cont[b[a[i]] + nx * j]);
                  if (!errors2.empty())
                     fSumw2.fArray[bin] = errors2[b[a[i]] + nx * j];
               }
            }
         } else {
            for (i = 0; i < nx; i++) {
               for (j = 0; j < n; j++) {
                  Int_t bin = GetBin(i, j + 1);
                  UpdateBinContent(bin, cont[i + nx * b[a[j]]]);
                  if (!errors2.empty())
                     fSumw2.fArray[bin] = errors2[i + nx * b[a[j]]];
               }
            }
         }
      } else {
         // case of 3D (needs to be tested)
         Int_t nx = fXaxis.GetNbins() + 2;
         Int_t ny = fYaxis.GetNbins() + 2;
         Int_t nz = fZaxis.GetNbins() + 2;
         cont.resize(nx * ny * nz);
         if (fSumw2.fN)
            errors2.resize(nx * ny * nz);
         for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
               for (k = 0; k < nz; k++) {
                  Int_t bin = GetBin(i, j, k);
                  cont[i + nx * (j + ny * k)] = RetrieveBinContent(bin);
                  if (!errors2.empty())
                     errors2[i + nx * (j + ny * k)] = GetBinErrorSqUnchecked(bin);
               }
            }
         }
         if (axis == GetXaxis()) {
            // labels on x axis
            for (i = 0; i < n; i++) { // for x we lool only on bins with the labels
               for (j = 0; j < ny; j++) {
                  for (k = 0; k < nz; k++) {
                     Int_t bin = GetBin(i + 1, j, k);
                     UpdateBinContent(bin, cont[b[a[i]] + nx * (j + ny * k)]);
                     if (!errors2.empty())
                        fSumw2.fArray[bin] = errors2[b[a[i]] + nx * (j + ny * k)];
                  }
               }
            }
         } else if (axis == GetYaxis()) {
            // labels on y axis
            for (i = 0; i < nx; i++) {
               for (j = 0; j < n; j++) {
                  for (k = 0; k < nz; k++) {
                     Int_t bin = GetBin(i, j+1, k);
                     UpdateBinContent(bin, cont[i + nx * (b[a[j]] + ny * k)]);
                     if (!errors2.empty())
                        fSumw2.fArray[bin] = errors2[i + nx * (b[a[j]] + ny * k)];
                  }
               }
            }
         } else {
            // labels on z axis
            for (i = 0; i < nx; i++) {
               for (j = 0; j < ny; j++) {
                  for (k = 0; k < n; k++) {
                     Int_t bin = GetBin(i, j, k+1);
                     UpdateBinContent(bin, cont[i + nx * (j + ny * b[a[k]])]);
                     if (!errors2.empty())
                        fSumw2.fArray[bin] = errors2[i + nx * (j + ny * b[a[k]])];
                  }
               }
            }
         }
      }
   }
   // need to set to zero the statistics if axis has been sorted
   // see for example TH3::PutStats for definition of s vector
   bool labelsAreSorted = kFALSE;
   for (i = 0; i < n; ++i) {
      if (a[i] != i) {
         labelsAreSorted = kTRUE;
         break;
      }
   }
   if (labelsAreSorted) {
      double s[TH1::kNstat];
      GetStats(s);
      if (iaxis == 1) {
         s[2] = 0; // fTsumwx
         s[3] = 0; // fTsumwx2
         s[6] = 0; // fTsumwxy
         s[9] = 0; // fTsumwxz
      } else if (iaxis == 2) {
         s[4] = 0;  // fTsumwy
         s[5] = 0;  // fTsumwy2
         s[6] = 0;  // fTsumwxy
         s[10] = 0; // fTsumwyz
      } else if (iaxis == 3) {
         s[7] = 0;  // fTsumwz
         s[8] = 0;  // fTsumwz2
         s[9] = 0;  // fTsumwxz
         s[10] = 0; // fTsumwyz
      }
      PutStats(s);
   }
   delete labold;
}

////////////////////////////////////////////////////////////////////////////////
/// Test if two double are almost equal.

static inline Bool_t AlmostEqual(Double_t a, Double_t b, Double_t epsilon = 0.00000001)
{
   return TMath::Abs(a - b) < epsilon;
}

////////////////////////////////////////////////////////////////////////////////
/// Test if a double is almost an integer.

static inline Bool_t AlmostInteger(Double_t a, Double_t epsilon = 0.00000001)
{
   return AlmostEqual(a - TMath::Floor(a), 0, epsilon) ||
      AlmostEqual(a - TMath::Floor(a), 1, epsilon);
}

////////////////////////////////////////////////////////////////////////////////
/// Test if the binning is equidistant.

static inline bool IsEquidistantBinning(const TAxis& axis)
{
   // check if axis bin are equals
   if (!axis.GetXbins()->fN) return true;  //
   // not able to check if there is only one axis entry
   bool isEquidistant = true;
   const Double_t firstBinWidth = axis.GetBinWidth(1);
   for (int i = 1; i < axis.GetNbins(); ++i) {
      const Double_t binWidth = axis.GetBinWidth(i);
      const bool match = TMath::AreEqualRel(firstBinWidth, binWidth, 1.E-10);
      isEquidistant &= match;
      if (!match)
         break;
   }
   return isEquidistant;
}

////////////////////////////////////////////////////////////////////////////////
/// Same limits and bins.

Bool_t TH1::SameLimitsAndNBins(const TAxis &axis1, const TAxis &axis2){
   return axis1.GetNbins() == axis2.GetNbins() &&
          TMath::AreEqualAbs(axis1.GetXmin(), axis2.GetXmin(), axis1.GetBinWidth(axis1.GetNbins()) * 1.E-10) &&
          TMath::AreEqualAbs(axis1.GetXmax(), axis2.GetXmax(), axis1.GetBinWidth(axis1.GetNbins()) * 1.E-10);
}

////////////////////////////////////////////////////////////////////////////////
/// Finds new limits for the axis for the Merge function.
/// returns false if the limits are incompatible

Bool_t TH1::RecomputeAxisLimits(TAxis &destAxis, const TAxis &anAxis)
{
   if (SameLimitsAndNBins(destAxis, anAxis))
      return kTRUE;

   if (!IsEquidistantBinning(destAxis) || !IsEquidistantBinning(anAxis))
      return kFALSE;       // not equidistant user binning not supported

   Double_t width1 = destAxis.GetBinWidth(0);
   Double_t width2 = anAxis.GetBinWidth(0);
   if (width1 == 0 || width2 == 0)
      return kFALSE;       // no binning not supported

   Double_t xmin = TMath::Min(destAxis.GetXmin(), anAxis.GetXmin());
   Double_t xmax = TMath::Max(destAxis.GetXmax(), anAxis.GetXmax());
   Double_t width = TMath::Max(width1, width2);

   // check the bin size
   if (!AlmostInteger(width/width1) || !AlmostInteger(width/width2))
      return kFALSE;

   // std::cout << "Find new limit using given axis " << anAxis.GetXmin() << " , " <<  anAxis.GetXmax() << " bin width " << width2 << std::endl;
   // std::cout << "           and destination axis " << destAxis.GetXmin() << " , " <<  destAxis.GetXmax() << " bin width " << width1 << std::endl;


   // check the limits
   Double_t delta;
   delta = (destAxis.GetXmin() - xmin)/width1;
   if (!AlmostInteger(delta))
      xmin -= (TMath::Ceil(delta) - delta)*width1;

   delta = (anAxis.GetXmin() - xmin)/width2;
   if (!AlmostInteger(delta))
      xmin -= (TMath::Ceil(delta) - delta)*width2;


   delta = (destAxis.GetXmin() - xmin)/width1;
   if (!AlmostInteger(delta))
      return kFALSE;


   delta = (xmax - destAxis.GetXmax())/width1;
   if (!AlmostInteger(delta))
      xmax += (TMath::Ceil(delta) - delta)*width1;


   delta = (xmax - anAxis.GetXmax())/width2;
   if (!AlmostInteger(delta))
      xmax += (TMath::Ceil(delta) - delta)*width2;


   delta = (xmax - destAxis.GetXmax())/width1;
   if (!AlmostInteger(delta))
      return kFALSE;
#ifdef DEBUG
   if (!AlmostInteger((xmax - xmin) / width)) {   // unnecessary check
      printf("TH1::RecomputeAxisLimits - Impossible\n");
      return kFALSE;
   }
#endif


   destAxis.Set(TMath::Nint((xmax - xmin)/width), xmin, xmax);

   //std::cout << "New re-computed axis : [ " << xmin << " , " << xmax << " ] width = " << width << " nbins " << destAxis.GetNbins() << std::endl;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add all histograms in the collection to this histogram.
/// This function computes the min/max for the x axis,
/// compute a new number of bins, if necessary,
/// add bin contents, errors and statistics.
/// If all histograms have bin labels, bins with identical labels
/// will be merged, no matter what their order is.
/// If overflows are present and limits are different the function will fail.
/// The function returns the total number of entries in the result histogram
/// if the merge is successful, -1 otherwise.
///
/// Possible option:
///   -NOL : the merger will ignore the labels and merge the histograms bin by bin using bin center values to match bins
///   -NOCHECK:  the histogram will not perform a check for duplicate labels in case of axes with labels. The check
///              (enabled by default) slows down the merging
///
/// IMPORTANT remark. The axis x may have different number
/// of bins and different limits, BUT the largest bin width must be
/// a multiple of the smallest bin width and the upper limit must also
/// be a multiple of the bin width.
/// Example:
///
/// ~~~ {.cpp}
/// void atest() {
///   TH1F *h1 = new TH1F("h1","h1",110,-110,0);
///   TH1F *h2 = new TH1F("h2","h2",220,0,110);
///   TH1F *h3 = new TH1F("h3","h3",330,-55,55);
///   TRandom r;
///   for (Int_t i=0;i<10000;i++) {
///      h1->Fill(r.Gaus(-55,10));
///      h2->Fill(r.Gaus(55,10));
///      h3->Fill(r.Gaus(0,10));
///   }
///
///   TList *list = new TList;
///   list->Add(h1);
///   list->Add(h2);
///   list->Add(h3);
///   TH1F *h = (TH1F*)h1->Clone("h");
///   h->Reset();
///   h->Merge(list);
///   h->Draw();
/// }
/// ~~~

Long64_t TH1::Merge(TCollection *li,Option_t * opt)
{
    if (!li) return 0;
    if (li->IsEmpty()) return (Long64_t) GetEntries();

    // use TH1Merger class
    TH1Merger merger(*this,*li,opt);
    Bool_t ret =  merger();

    return (ret) ? GetEntries() : -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Performs the operation:
///
/// `this = this*c1*f1`
///
/// If errors are defined (see TH1::Sumw2), errors are also recalculated.
///
/// Only bins inside the function range are recomputed.
/// IMPORTANT NOTE: If you intend to use the errors of this histogram later
/// you should call Sumw2 before making this operation.
/// This is particularly important if you fit the histogram after TH1::Multiply
///
/// The function return kFALSE if the Multiply operation failed

Bool_t TH1::Multiply(TF1 *f1, Double_t c1)
{
   if (!f1) {
      Error("Multiply","Attempt to multiply by a non-existing function");
      return kFALSE;
   }

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   Int_t nx = GetNbinsX() + 2; // normal bins + uf / of (cells)
   Int_t ny = GetNbinsY() + 2;
   Int_t nz = GetNbinsZ() + 2;
   if (fDimension < 2) ny = 1;
   if (fDimension < 3) nz = 1;

   // reset min-maximum
   SetMinimum();
   SetMaximum();

   //   - Loop on bins (including underflows/overflows)
   Double_t xx[3];
   Double_t *params = 0;
   f1->InitArgs(xx,params);

   for (Int_t binz = 0; binz < nz; ++binz) {
      xx[2] = fZaxis.GetBinCenter(binz);
      for (Int_t biny = 0; biny < ny; ++biny) {
         xx[1] = fYaxis.GetBinCenter(biny);
         for (Int_t binx = 0; binx < nx; ++binx) {
            xx[0] = fXaxis.GetBinCenter(binx);
            if (!f1->IsInside(xx)) continue;
            TF1::RejectPoint(kFALSE);
            Int_t bin = binx + nx * (biny + ny *binz);
            Double_t cu  = c1*f1->EvalPar(xx);
            if (TF1::RejectedPoint()) continue;
            UpdateBinContent(bin, RetrieveBinContent(bin) * cu);
            if (fSumw2.fN) {
               fSumw2.fArray[bin] = cu * cu * GetBinErrorSqUnchecked(bin);
            }
         }
      }
   }
   ResetStats();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply this histogram by h1.
///
/// `this = this*h1`
///
/// If errors of this are available (TH1::Sumw2), errors are recalculated.
/// Note that if h1 has Sumw2 set, Sumw2 is automatically called for this
/// if not already set.
///
/// IMPORTANT NOTE: If you intend to use the errors of this histogram later
/// you should call Sumw2 before making this operation.
/// This is particularly important if you fit the histogram after TH1::Multiply
///
/// The function return kFALSE if the Multiply operation failed

Bool_t TH1::Multiply(const TH1 *h1)
{
   if (!h1) {
      Error("Multiply","Attempt to multiply by a non-existing histogram");
      return kFALSE;
   }

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   try {
      CheckConsistency(this,h1);
   } catch(DifferentNumberOfBins&) {
      Error("Multiply","Attempt to multiply histograms with different number of bins");
      return kFALSE;
   } catch(DifferentAxisLimits&) {
      Warning("Multiply","Attempt to multiply histograms with different axis limits");
   } catch(DifferentBinLimits&) {
      Warning("Multiply","Attempt to multiply histograms with different bin limits");
   } catch(DifferentLabels&) {
      Warning("Multiply","Attempt to multiply histograms with different labels");
   }

   //    Create Sumw2 if h1 has Sumw2 set
   if (fSumw2.fN == 0 && h1->GetSumw2N() != 0) Sumw2();

   //   - Reset min-  maximum
   SetMinimum();
   SetMaximum();

   //   - Loop on bins (including underflows/overflows)
   for (Int_t i = 0; i < fNcells; ++i) {
      Double_t c0 = RetrieveBinContent(i);
      Double_t c1 = h1->RetrieveBinContent(i);
      UpdateBinContent(i, c0 * c1);
      if (fSumw2.fN) {
         fSumw2.fArray[i] = GetBinErrorSqUnchecked(i) * c1 * c1 + h1->GetBinErrorSqUnchecked(i) * c0 * c0;
      }
   }
   ResetStats();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace contents of this histogram by multiplication of h1 by h2.
///
/// `this = (c1*h1)*(c2*h2)`
///
/// If errors of this are available (TH1::Sumw2), errors are recalculated.
/// Note that if h1 or h2 have Sumw2 set, Sumw2 is automatically called for this
/// if not already set.
///
/// IMPORTANT NOTE: If you intend to use the errors of this histogram later
/// you should call Sumw2 before making this operation.
/// This is particularly important if you fit the histogram after TH1::Multiply
///
/// The function return kFALSE if the Multiply operation failed

Bool_t TH1::Multiply(const TH1 *h1, const TH1 *h2, Double_t c1, Double_t c2, Option_t *option)
{
   TString opt = option;
   opt.ToLower();
   //   Bool_t binomial = kFALSE;
   //   if (opt.Contains("b")) binomial = kTRUE;
   if (!h1 || !h2) {
      Error("Multiply","Attempt to multiply by a non-existing histogram");
      return kFALSE;
   }

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   try {
      CheckConsistency(h1,h2);
      CheckConsistency(this,h1);
   } catch(DifferentNumberOfBins&) {
      Error("Multiply","Attempt to multiply histograms with different number of bins");
      return kFALSE;
   } catch(DifferentAxisLimits&) {
      Warning("Multiply","Attempt to multiply histograms with different axis limits");
   } catch(DifferentBinLimits&) {
      Warning("Multiply","Attempt to multiply histograms with different bin limits");
   } catch(DifferentLabels&) {
      Warning("Multiply","Attempt to multiply histograms with different labels");
   }

   //    Create Sumw2 if h1 or h2 have Sumw2 set
   if (fSumw2.fN == 0 && (h1->GetSumw2N() != 0 || h2->GetSumw2N() != 0)) Sumw2();

   //   - Reset min - maximum
   SetMinimum();
   SetMaximum();

   //   - Loop on bins (including underflows/overflows)
   Double_t c1sq = c1 * c1; Double_t c2sq = c2 * c2;
   for (Int_t i = 0; i < fNcells; ++i) {
      Double_t b1 = h1->RetrieveBinContent(i);
      Double_t b2 = h2->RetrieveBinContent(i);
      UpdateBinContent(i, c1 * b1 * c2 * b2);
      if (fSumw2.fN) {
         fSumw2.fArray[i] = c1sq * c2sq * (h1->GetBinErrorSqUnchecked(i) * b2 * b2 + h2->GetBinErrorSqUnchecked(i) * b1 * b1);
      }
   }
   ResetStats();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Control routine to paint any kind of histograms.
///
/// This function is automatically called by TCanvas::Update.
/// (see TH1::Draw for the list of options)

void TH1::Paint(Option_t *option)
{
   GetPainter(option);

   if (fPainter) {
      if (strlen(option) > 0) fPainter->Paint(option);
      else                    fPainter->Paint(fOption.Data());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Rebin this histogram
///
/// #### case 1  xbins=0
///
/// If newname is blank (default), the current histogram is modified and
/// a pointer to it is returned.
///
/// If newname is not blank, the current histogram is not modified, and a
/// new histogram is returned which is a Clone of the current histogram
/// with its name set to newname.
///
/// The parameter ngroup indicates how many bins of this have to be merged
/// into one bin of the result.
///
/// If the original histogram has errors stored (via Sumw2), the resulting
/// histograms has new errors correctly calculated.
///
/// examples: if h1 is an existing TH1F histogram with 100 bins
///
/// ~~~ {.cpp}
///     h1->Rebin();  //merges two bins in one in h1: previous contents of h1 are lost
///     h1->Rebin(5); //merges five bins in one in h1
///     TH1F *hnew = dynamic_cast<TH1F*>(h1->Rebin(5,"hnew")); // creates a new histogram hnew
///                                                            // merging 5 bins of h1 in one bin
/// ~~~
///
/// NOTE:  If ngroup is not an exact divider of the number of bins,
/// the top limit of the rebinned histogram is reduced
/// to the upper edge of the last bin that can make a complete
/// group. The remaining bins are added to the overflow bin.
/// Statistics will be recomputed from the new bin contents.
///
/// #### case 2  xbins!=0
///
/// A new histogram is created (you should specify newname).
/// The parameter ngroup is the number of variable size bins in the created histogram.
/// The array xbins must contain ngroup+1 elements that represent the low-edges
/// of the bins.
/// If the original histogram has errors stored (via Sumw2), the resulting
/// histograms has new errors correctly calculated.
///
/// NOTE:  The bin edges specified in xbins should correspond to bin edges
/// in the original histogram. If a bin edge in the new histogram is
/// in the middle of a bin in the original histogram, all entries in
/// the split bin in the original histogram will be transfered to the
/// lower of the two possible bins in the new histogram. This is
/// probably not what you want. A warning message is emitted in this
/// case
///
/// examples: if h1 is an existing TH1F histogram with 100 bins
///
/// ~~~ {.cpp}
///     Double_t xbins[25] = {...} array of low-edges (xbins[25] is the upper edge of last bin
///     h1->Rebin(24,"hnew",xbins);  //creates a new variable bin size histogram hnew
/// ~~~

TH1 *TH1::Rebin(Int_t ngroup, const char*newname, const Double_t *xbins)
{
   Int_t nbins    = fXaxis.GetNbins();
   Double_t xmin  = fXaxis.GetXmin();
   Double_t xmax  = fXaxis.GetXmax();
   if ((ngroup <= 0) || (ngroup > nbins)) {
      Error("Rebin", "Illegal value of ngroup=%d",ngroup);
      return 0;
   }

   if (fDimension > 1 || InheritsFrom(TProfile::Class())) {
      Error("Rebin", "Operation valid on 1-D histograms only");
      return 0;
   }
   if (!newname && xbins) {
      Error("Rebin","if xbins is specified, newname must be given");
      return 0;
   }

   Int_t newbins = nbins/ngroup;
   if (!xbins) {
      Int_t nbg = nbins/ngroup;
      if (nbg*ngroup != nbins) {
         Warning("Rebin", "ngroup=%d is not an exact divider of nbins=%d.",ngroup,nbins);
      }
   }
   else {
   // in the case that xbins is given (rebinning in variable bins), ngroup is
   // the new number of bins and number of grouped bins is not constant.
   // when looping for setting the contents for the new histogram we
   // need to loop on all bins of original histogram.  Then set ngroup=nbins
      newbins = ngroup;
      ngroup = nbins;
   }

   // Save old bin contents into a new array
   Double_t entries = fEntries;
   Double_t *oldBins = new Double_t[nbins+2];
   Int_t bin, i;
   for (bin=0;bin<nbins+2;bin++) oldBins[bin] = RetrieveBinContent(bin);
   Double_t *oldErrors = 0;
   if (fSumw2.fN != 0) {
      oldErrors = new Double_t[nbins+2];
      for (bin=0;bin<nbins+2;bin++) oldErrors[bin] = GetBinError(bin);
   }
   // rebin will not include underflow/overflow if new axis range is larger than old axis range
   if (xbins) {
      if (xbins[0] < fXaxis.GetXmin() && oldBins[0] != 0 )
         Warning("Rebin","underflow entries will not be used when rebinning");
      if (xbins[newbins] > fXaxis.GetXmax() && oldBins[nbins+1] != 0 )
         Warning("Rebin","overflow entries will not be used when rebinning");
   }


   // create a clone of the old histogram if newname is specified
   TH1 *hnew = this;
   if ((newname && strlen(newname) > 0) || xbins) {
      hnew = (TH1*)Clone(newname);
   }

   //reset can extend bit to avoid an axis extension in SetBinContent
   UInt_t oldExtendBitMask = hnew->SetCanExtend(kNoAxis);

   // save original statistics
   Double_t stat[kNstat];
   GetStats(stat);
   bool resetStat = false;
   // change axis specs and rebuild bin contents array::RebinAx
   if(!xbins && (newbins*ngroup != nbins)) {
      xmax = fXaxis.GetBinUpEdge(newbins*ngroup);
      resetStat = true; //stats must be reset because top bins will be moved to overflow bin
   }
   // save the TAttAxis members (reset by SetBins)
   Int_t    nDivisions  = fXaxis.GetNdivisions();
   Color_t  axisColor   = fXaxis.GetAxisColor();
   Color_t  labelColor  = fXaxis.GetLabelColor();
   Style_t  labelFont   = fXaxis.GetLabelFont();
   Float_t  labelOffset = fXaxis.GetLabelOffset();
   Float_t  labelSize   = fXaxis.GetLabelSize();
   Float_t  tickLength  = fXaxis.GetTickLength();
   Float_t  titleOffset = fXaxis.GetTitleOffset();
   Float_t  titleSize   = fXaxis.GetTitleSize();
   Color_t  titleColor  = fXaxis.GetTitleColor();
   Style_t  titleFont   = fXaxis.GetTitleFont();

   if(!xbins && (fXaxis.GetXbins()->GetSize() > 0)){ // variable bin sizes
      Double_t *bins = new Double_t[newbins+1];
      for(i = 0; i <= newbins; ++i) bins[i] = fXaxis.GetBinLowEdge(1+i*ngroup);
      hnew->SetBins(newbins,bins); //this also changes errors array (if any)
      delete [] bins;
   } else if (xbins) {
      hnew->SetBins(newbins,xbins);
   } else {
      hnew->SetBins(newbins,xmin,xmax);
   }

   // Restore axis attributes
   fXaxis.SetNdivisions(nDivisions);
   fXaxis.SetAxisColor(axisColor);
   fXaxis.SetLabelColor(labelColor);
   fXaxis.SetLabelFont(labelFont);
   fXaxis.SetLabelOffset(labelOffset);
   fXaxis.SetLabelSize(labelSize);
   fXaxis.SetTickLength(tickLength);
   fXaxis.SetTitleOffset(titleOffset);
   fXaxis.SetTitleSize(titleSize);
   fXaxis.SetTitleColor(titleColor);
   fXaxis.SetTitleFont(titleFont);

   // copy merged bin contents (ignore under/overflows)
   // Start merging only once the new lowest edge is reached
   Int_t startbin = 1;
   const Double_t newxmin = hnew->GetXaxis()->GetBinLowEdge(1);
   while( fXaxis.GetBinCenter(startbin) < newxmin && startbin <= nbins ) {
      startbin++;
   }
   Int_t oldbin = startbin;
   Double_t binContent, binError;
   for (bin = 1;bin<=newbins;bin++) {
      binContent = 0;
      binError   = 0;
      Int_t imax = ngroup;
      Double_t xbinmax = hnew->GetXaxis()->GetBinUpEdge(bin);
      // check bin edges for the cases when we provide an array of bins
      // be careful in case bins can have zero width
      if (xbins && !TMath::AreEqualAbs(fXaxis.GetBinLowEdge(oldbin),
                                       hnew->GetXaxis()->GetBinLowEdge(bin),
                                       TMath::Max(1.E-8 * fXaxis.GetBinWidth(oldbin), 1.E-16 )) )
      {
         Warning("Rebin","Bin edge %d of rebinned histogram does not match any bin edges of the old histogram. Result can be inconsistent",bin);
      }
      for (i=0;i<ngroup;i++) {
         if( (oldbin+i > nbins) ||
             ( hnew != this && (fXaxis.GetBinCenter(oldbin+i) > xbinmax)) ) {
            imax = i;
            break;
         }
         binContent += oldBins[oldbin+i];
         if (oldErrors) binError += oldErrors[oldbin+i]*oldErrors[oldbin+i];
      }
      hnew->SetBinContent(bin,binContent);
      if (oldErrors) hnew->SetBinError(bin,TMath::Sqrt(binError));
      oldbin += imax;
   }

   // sum underflow and overflow contents until startbin
   binContent = 0;
   binError = 0;
   for (i = 0; i < startbin; ++i)  {
      binContent += oldBins[i];
      if (oldErrors) binError += oldErrors[i]*oldErrors[i];
   }
   hnew->SetBinContent(0,binContent);
   if (oldErrors) hnew->SetBinError(0,TMath::Sqrt(binError));
   // sum overflow
   binContent = 0;
   binError = 0;
   for (i = oldbin; i <= nbins+1; ++i)  {
      binContent += oldBins[i];
      if (oldErrors) binError += oldErrors[i]*oldErrors[i];
   }
   hnew->SetBinContent(newbins+1,binContent);
   if (oldErrors) hnew->SetBinError(newbins+1,TMath::Sqrt(binError));

   hnew->SetCanExtend(oldExtendBitMask); // restore previous state

   // restore statistics and entries modified by SetBinContent
   hnew->SetEntries(entries);
   if (!resetStat) hnew->PutStats(stat);
   delete [] oldBins;
   if (oldErrors) delete [] oldErrors;
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// finds new limits for the axis so that *point* is within the range and
/// the limits are compatible with the previous ones (see TH1::Merge).
/// new limits are put into *newMin* and *newMax* variables.
/// axis - axis whose limits are to be recomputed
/// point - point that should fit within the new axis limits
/// newMin - new minimum will be stored here
/// newMax - new maximum will be stored here.
/// false if failed (e.g. if the initial axis limits are wrong
/// or the new range is more than \f$ 2^{64} \f$ times the old one).

Bool_t TH1::FindNewAxisLimits(const TAxis* axis, const Double_t point, Double_t& newMin, Double_t &newMax)
{
   Double_t xmin = axis->GetXmin();
   Double_t xmax = axis->GetXmax();
   if (xmin >= xmax) return kFALSE;
   Double_t range = xmax-xmin;
   Double_t binsize = range / axis->GetNbins();

   //recompute new axis limits by doubling the current range
   Int_t ntimes = 0;
   while (point < xmin) {
      if (ntimes++ > 64)
         return kFALSE;
      xmin = xmin - range;
      range *= 2;
      binsize *= 2;
      // // make sure that the merging will be correct
      // if ( xmin / binsize - TMath::Floor(xmin / binsize) >= 0.5) {
      //    xmin += 0.5 * binsize;
      //    xmax += 0.5 * binsize;  // won't work with a histogram with only one bin, but I don't care
      // }
   }
   while (point >= xmax) {
      if (ntimes++ > 64)
         return kFALSE;
      xmax = xmax + range;
      range *= 2;
      binsize *= 2;
      // // make sure that the merging will be correct
      // if ( xmin / binsize - TMath::Floor(xmin / binsize) >= 0.5) {
      //    xmin -= 0.5 * binsize;
      //    xmax -= 0.5 * binsize;  // won't work with a histogram with only one bin, but I don't care
      // }
   }
   newMin = xmin;
   newMax = xmax;
   //   Info("FindNewAxisLimits", "OldAxis: (%lf, %lf), new: (%lf, %lf), point: %lf",
   //      axis->GetXmin(), axis->GetXmax(), xmin, xmax, point);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Histogram is resized along axis such that x is in the axis range.
/// The new axis limits are recomputed by doubling iteratively
/// the current axis range until the specified value x is within the limits.
/// The algorithm makes a copy of the histogram, then loops on all bins
/// of the old histogram to fill the extended histogram.
/// Takes into account errors (Sumw2) if any.
/// The algorithm works for 1-d, 2-D and 3-D histograms.
/// The axis must be extendable before invoking this function.
/// Ex:
///
/// ~~~ {.cpp}
/// h->GetXaxis()->SetCanExtend(kTRUE);
/// ~~~

void TH1::ExtendAxis(Double_t x, TAxis *axis)
{
   if (!axis->CanExtend()) return;
   if (TMath::IsNaN(x)) {         // x may be a NaN
      SetCanExtend(kNoAxis);
      return;
   }

   if (axis->GetXmin() >= axis->GetXmax()) return;
   if (axis->GetNbins() <= 0) return;

   Double_t xmin, xmax;
   if (!FindNewAxisLimits(axis, x, xmin, xmax))
      return;

   //save a copy of this histogram
   TH1 *hold = (TH1*)IsA()->New();
   hold->SetDirectory(0);
   Copy(*hold);
   //set new axis limits
   axis->SetLimits(xmin,xmax);


   //now loop on all bins and refill
   Int_t errors = GetSumw2N();

   Reset("ICE"); //reset only Integral, contents and Errors

   int iaxis = 0;
   if (axis == &fXaxis) iaxis = 1;
   if (axis == &fYaxis) iaxis = 2;
   if (axis == &fZaxis) iaxis = 3;
   bool firstw = kTRUE;
   Int_t binx,biny, binz = 0;
   Int_t ix = 0,iy = 0,iz = 0;
   Double_t bx,by,bz;
   Int_t ncells = hold->GetNcells();
   for (Int_t bin = 0; bin < ncells; ++bin) {
      hold->GetBinXYZ(bin,binx,biny,binz);
      bx = hold->GetXaxis()->GetBinCenter(binx);
      ix  = fXaxis.FindFixBin(bx);
      if (fDimension > 1) {
         by  = hold->GetYaxis()->GetBinCenter(biny);
         iy  = fYaxis.FindFixBin(by);
         if (fDimension > 2) {
            bz  = hold->GetZaxis()->GetBinCenter(binz);
            iz  = fZaxis.FindFixBin(bz);
         }
      }
      // exclude underflow/overflow
      double content = hold->RetrieveBinContent(bin);
      if (content == 0) continue;
      if (IsBinUnderflow(bin,iaxis) || IsBinOverflow(bin,iaxis) ) {
         if (firstw) {
            Warning("ExtendAxis","Histogram %s has underflow or overflow in the axis that is extendable"
                    " their content will be lost",GetName() );
            firstw= kFALSE;
         }
         continue;
      }
      Int_t ibin= GetBin(ix,iy,iz);
      AddBinContent(ibin, content);
      if (errors) {
         fSumw2.fArray[ibin] += hold->GetBinErrorSqUnchecked(bin);
      }
   }
   delete hold;
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively remove object from the list of functions

void TH1::RecursiveRemove(TObject *obj)
{
   // Rely on TROOT::RecursiveRemove to take the readlock.

   if (fFunctions) {
      if (!fFunctions->TestBit(kInvalidObject)) fFunctions->RecursiveRemove(obj);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply this histogram by a constant c1.
///
/// `this = c1*this`
///
/// Note that both contents and errors (if any) are scaled.
/// This function uses the services of TH1::Add
///
/// IMPORTANT NOTE: Sumw2() is called automatically when scaling.
/// If you are not interested in the histogram statistics you can call
/// Sumw2(kFALSE) or use the option "nosw2"
///
/// One can scale an histogram such that the bins integral is equal to
/// the normalization parameter via TH1::Scale(Double_t norm), where norm
/// is the desired normalization divided by the integral of the histogram.
///
/// If option contains "width" the bin contents and errors are divided
/// by the bin width.

void TH1::Scale(Double_t c1, Option_t *option)
{

   TString opt = option; opt.ToLower();
   // store bin errors when scaling since cannot anymore be computed as sqrt(N)
   if (!opt.Contains("nosw2") && GetSumw2N() == 0) Sumw2();
   if (opt.Contains("width")) Add(this, this, c1, -1);
   else {
      if (fBuffer) BufferEmpty(1);
      for(Int_t i = 0; i < fNcells; ++i) UpdateBinContent(i, c1 * RetrieveBinContent(i));
      if (fSumw2.fN) for(Int_t i = 0; i < fNcells; ++i) fSumw2.fArray[i] *= (c1 * c1); // update errors
      // update global histograms statistics
      Double_t s[kNstat] = {0};
      GetStats(s);
      for (Int_t i=0 ; i < kNstat; i++) {
         if (i == 1)   s[i] = c1*c1*s[i];
         else          s[i] = c1*s[i];
      }
      PutStats(s);
      SetMinimum(); SetMaximum(); // minimum and maximum value will be recalculated the next time
   }

   // if contours set, must also scale contours
   Int_t ncontours = GetContour();
   if (ncontours == 0) return;
   Double_t* levels = fContour.GetArray();
   for (Int_t i = 0; i < ncontours; ++i) levels[i] *= c1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if all axes are extendable.

Bool_t TH1::CanExtendAllAxes() const
{
   Bool_t canExtend = fXaxis.CanExtend();
   if (GetDimension() > 1) canExtend &= fYaxis.CanExtend();
   if (GetDimension() > 2) canExtend &= fZaxis.CanExtend();

   return canExtend;
}

////////////////////////////////////////////////////////////////////////////////
/// Make the histogram axes extendable / not extendable according to the bit mask
/// returns the previous bit mask specifying which axes are extendable

UInt_t TH1::SetCanExtend(UInt_t extendBitMask)
{
   UInt_t oldExtendBitMask = kNoAxis;

   if (fXaxis.CanExtend()) oldExtendBitMask |= kXaxis;
   if (extendBitMask & kXaxis) fXaxis.SetCanExtend(kTRUE);
   else fXaxis.SetCanExtend(kFALSE);

   if (GetDimension() > 1) {
      if (fYaxis.CanExtend()) oldExtendBitMask |= kYaxis;
      if (extendBitMask & kYaxis) fYaxis.SetCanExtend(kTRUE);
      else fYaxis.SetCanExtend(kFALSE);
   }

   if (GetDimension() > 2) {
      if (fZaxis.CanExtend()) oldExtendBitMask |= kZaxis;
      if (extendBitMask & kZaxis) fZaxis.SetCanExtend(kTRUE);
      else fZaxis.SetCanExtend(kFALSE);
   }

   return oldExtendBitMask;
}

///////////////////////////////////////////////////////////////////////////////
/// Internal function used in TH1::Fill to see which axis is full alphanumeric
/// i.e. can be extended and is alphanumeric
UInt_t TH1::GetAxisLabelStatus() const
{
   UInt_t bitMask = kNoAxis;
   if (fXaxis.CanExtend() && fXaxis.IsAlphanumeric() ) bitMask |= kXaxis;
   if (GetDimension() > 1 && fYaxis.CanExtend() && fYaxis.IsAlphanumeric())
      bitMask |= kYaxis;
   if (GetDimension() > 2 && fZaxis.CanExtend() && fZaxis.IsAlphanumeric())
      bitMask |= kZaxis;

   return bitMask;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function to set the default buffer size for automatic histograms.
/// When an histogram is created with one of its axis lower limit greater
/// or equal to its upper limit, the function SetBuffer is automatically
/// called with the default buffer size.

void TH1::SetDefaultBufferSize(Int_t buffersize)
{
   fgBufferSize = buffersize > 0 ? buffersize : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// When this static function is called with `sumw2=kTRUE`, all new
/// histograms will automatically activate the storage
/// of the sum of squares of errors, ie TH1::Sumw2 is automatically called.

void TH1::SetDefaultSumw2(Bool_t sumw2)
{
   fgDefaultSumw2 = sumw2;
}

////////////////////////////////////////////////////////////////////////////////
/// Change (i.e. set) the title
///
/// if title is in the form `stringt;stringx;stringy;stringz`
/// the histogram title is set to `stringt`, the x axis title to `stringx`,
/// the y axis title to `stringy`, and the z axis title to `stringz`.
///
/// To insert the character `;` in one of the titles, one should use `#;`
/// or `#semicolon`.

void TH1::SetTitle(const char *title)
{
   fTitle = title;
   fTitle.ReplaceAll("#;",2,"#semicolon",10);

   // Decode fTitle. It may contain X, Y and Z titles
   TString str1 = fTitle, str2;
   Int_t isc = str1.Index(";");
   Int_t lns = str1.Length();

   if (isc >=0 ) {
      fTitle = str1(0,isc);
      str1   = str1(isc+1, lns);
      isc    = str1.Index(";");
      if (isc >=0 ) {
         str2 = str1(0,isc);
         str2.ReplaceAll("#semicolon",10,";",1);
         fXaxis.SetTitle(str2.Data());
         lns  = str1.Length();
         str1 = str1(isc+1, lns);
         isc  = str1.Index(";");
         if (isc >=0 ) {
            str2 = str1(0,isc);
            str2.ReplaceAll("#semicolon",10,";",1);
            fYaxis.SetTitle(str2.Data());
            lns  = str1.Length();
            str1 = str1(isc+1, lns);
            str1.ReplaceAll("#semicolon",10,";",1);
            fZaxis.SetTitle(str1.Data());
         } else {
            str1.ReplaceAll("#semicolon",10,";",1);
            fYaxis.SetTitle(str1.Data());
         }
      } else {
         str1.ReplaceAll("#semicolon",10,";",1);
         fXaxis.SetTitle(str1.Data());
      }
   }

   fTitle.ReplaceAll("#semicolon",10,";",1);

   if (gPad && TestBit(kMustCleanup)) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// Smooth array xx, translation of Hbook routine hsmoof.F
/// based on algorithm 353QH twice presented by J. Friedman
/// in Proc.of the 1974 CERN School of Computing, Norway, 11-24 August, 1974.

void  TH1::SmoothArray(Int_t nn, Double_t *xx, Int_t ntimes)
{
   if (nn < 3 ) {
      ::Error("SmoothArray","Need at least 3 points for smoothing: n = %d",nn);
      return;
   }

   Int_t ii;
   Double_t hh[6] = {0,0,0,0,0,0};

   std::vector<double> yy(nn);
   std::vector<double> zz(nn);
   std::vector<double> rr(nn);

   for (Int_t pass=0;pass<ntimes;pass++) {
      // first copy original data into temp array
      std::copy(xx, xx+nn, zz.begin() );

      for (int noent = 0; noent < 2; ++noent) { // run algorithm two times

         //  do 353 i.e. running median 3, 5, and 3 in a single loop
         for  (int kk = 0; kk < 3; kk++)  {
            std::copy(zz.begin(), zz.end(), yy.begin());
            int medianType = (kk != 1)  ?  3 : 5;
            int ifirst      = (kk != 1 ) ?  1 : 2;
            int ilast       = (kk != 1 ) ? nn-1 : nn -2;
            //nn2 = nn - ik - 1;
            // do all elements beside the first and last point for median 3
            //  and first two and last 2 for median 5
            for  ( ii = ifirst; ii < ilast; ii++)  {
               assert(ii - ifirst >= 0);
               for  (int jj = 0; jj < medianType; jj++)   {
                  hh[jj] = yy[ii - ifirst + jj ];
               }
               zz[ii] = TMath::Median(medianType, hh);
            }

            if  (kk == 0)  {   // first median 3
               // first point
               hh[0] = zz[1];
               hh[1] = zz[0];
               hh[2] = 3*zz[1] - 2*zz[2];
               zz[0] = TMath::Median(3, hh);
               // last point
               hh[0] = zz[nn - 2];
               hh[1] = zz[nn - 1];
               hh[2] = 3*zz[nn - 2] - 2*zz[nn - 3];
               zz[nn - 1] = TMath::Median(3, hh);
            }

            if  (kk == 1)  {   //  median 5
               for  (ii = 0; ii < 3; ii++) {
                  hh[ii] = yy[ii];
               }
               zz[1] = TMath::Median(3, hh);
               // last two points
               for  (ii = 0; ii < 3; ii++) {
                  hh[ii] = yy[nn - 3 + ii];
               }
               zz[nn - 2] = TMath::Median(3, hh);
            }

         }

         std::copy ( zz.begin(), zz.end(), yy.begin() );

         // quadratic interpolation for flat segments
         for (ii = 2; ii < (nn - 2); ii++) {
            if  (zz[ii - 1] != zz[ii]) continue;
            if  (zz[ii] != zz[ii + 1]) continue;
            hh[0] = zz[ii - 2] - zz[ii];
            hh[1] = zz[ii + 2] - zz[ii];
            if  (hh[0] * hh[1] <= 0) continue;
            int jk = 1;
            if  ( TMath::Abs(hh[1]) > TMath::Abs(hh[0]) ) jk = -1;
            yy[ii] = -0.5*zz[ii - 2*jk] + zz[ii]/0.75 + zz[ii + 2*jk] /6.;
            yy[ii + jk] = 0.5*(zz[ii + 2*jk] - zz[ii - 2*jk]) + zz[ii];
         }

         // running means
         //std::copy(zz.begin(), zz.end(), yy.begin());
         for  (ii = 1; ii < nn - 1; ii++) {
            zz[ii] = 0.25*yy[ii - 1] + 0.5*yy[ii] + 0.25*yy[ii + 1];
         }
         zz[0] = yy[0];
         zz[nn - 1] = yy[nn - 1];

         if (noent == 0) {

            // save computed values
            std::copy(zz.begin(), zz.end(), rr.begin());

            // COMPUTE  residuals
            for  (ii = 0; ii < nn; ii++)  {
               zz[ii] = xx[ii] - zz[ii];
            }
         }

      }  // end loop on noent


      double xmin = TMath::MinElement(nn,xx);
      for  (ii = 0; ii < nn; ii++) {
         if (xmin < 0) xx[ii] = rr[ii] + zz[ii];
         // make smoothing defined positive - not better using 0 ?
         else  xx[ii] = TMath::Max((rr[ii] + zz[ii]),0.0 );
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Smooth bin contents of this histogram.
/// if option contains "R" smoothing is applied only to the bins
/// defined in the X axis range (default is to smooth all bins)
/// Bin contents are replaced by their smooth values.
/// Errors (if any) are not modified.
/// the smoothing procedure is repeated ntimes (default=1)

void  TH1::Smooth(Int_t ntimes, Option_t *option)
{
   if (fDimension != 1) {
      Error("Smooth","Smooth only supported for 1-d histograms");
      return;
   }
   Int_t nbins = fXaxis.GetNbins();
   if (nbins < 3) {
      Error("Smooth","Smooth only supported for histograms with >= 3 bins. Nbins = %d",nbins);
      return;
   }

   // delete buffer if it is there since it will become invalid
   if (fBuffer) BufferEmpty(1);

   Int_t firstbin = 1, lastbin = nbins;
   TString opt = option;
   opt.ToLower();
   if (opt.Contains("r")) {
      firstbin= fXaxis.GetFirst();
      lastbin  = fXaxis.GetLast();
   }
   nbins = lastbin - firstbin + 1;
   Double_t *xx = new Double_t[nbins];
   Double_t nent = fEntries;
   Int_t i;
   for (i=0;i<nbins;i++) {
      xx[i] = RetrieveBinContent(i+firstbin);
   }

   TH1::SmoothArray(nbins,xx,ntimes);

   for (i=0;i<nbins;i++) {
      UpdateBinContent(i+firstbin,xx[i]);
   }
   fEntries = nent;
   delete [] xx;

   if (gPad) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
/// if flag=kTRUE, underflows and overflows are used by the Fill functions
/// in the computation of statistics (mean value, StdDev).
/// By default, underflows or overflows are not used.

void  TH1::StatOverflows(Bool_t flag)
{
   fgStatOverflows = flag;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.

void TH1::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (fDirectory) fDirectory->Remove(this);
      fDirectory = 0;
      if (R__v > 2) {
         b.ReadClassBuffer(TH1::Class(), this, R__v, R__s, R__c);

         ResetBit(kMustCleanup);
         fXaxis.SetParent(this);
         fYaxis.SetParent(this);
         fZaxis.SetParent(this);
         TIter next(fFunctions);
         TObject *obj;
         while ((obj=next())) {
            if (obj->InheritsFrom(TF1::Class())) ((TF1*)obj)->SetParent(this);
         }
         return;
      }
      //process old versions before automatic schema evolution
      TNamed::Streamer(b);
      TAttLine::Streamer(b);
      TAttFill::Streamer(b);
      TAttMarker::Streamer(b);
      b >> fNcells;
      fXaxis.Streamer(b);
      fYaxis.Streamer(b);
      fZaxis.Streamer(b);
      fXaxis.SetParent(this);
      fYaxis.SetParent(this);
      fZaxis.SetParent(this);
      b >> fBarOffset;
      b >> fBarWidth;
      b >> fEntries;
      b >> fTsumw;
      b >> fTsumw2;
      b >> fTsumwx;
      b >> fTsumwx2;
      if (R__v < 2) {
         Float_t maximum, minimum, norm;
         Float_t *contour=0;
         b >> maximum; fMaximum = maximum;
         b >> minimum; fMinimum = minimum;
         b >> norm;    fNormFactor = norm;
         Int_t n = b.ReadArray(contour);
         fContour.Set(n);
         for (Int_t i=0;i<n;i++) fContour.fArray[i] = contour[i];
         delete [] contour;
      } else {
         b >> fMaximum;
         b >> fMinimum;
         b >> fNormFactor;
         fContour.Streamer(b);
      }
      fSumw2.Streamer(b);
      fOption.Streamer(b);
      fFunctions->Delete();
      fFunctions->Streamer(b);
      b.CheckByteCount(R__s, R__c, TH1::IsA());

   } else {
      b.WriteClassBuffer(TH1::Class(),this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print some global quantities for this histogram.
/// \param[in] option
///   - "base" is given, number of bins and ranges are also printed
///   - "range" is given, bin contents and errors are also printed
///     for all bins in the current range (default 1-->nbins)
///   - "all" is given, bin contents and errors are also printed
///     for all bins including under and overflows.

void TH1::Print(Option_t *option) const
{
   if (fBuffer) const_cast<TH1*>(this)->BufferEmpty();
   printf( "TH1.Print Name  = %s, Entries= %d, Total sum= %g\n",GetName(),Int_t(fEntries),GetSumOfWeights());
   TString opt = option;
   opt.ToLower();
   Int_t all;
   if      (opt.Contains("all"))   all = 0;
   else if (opt.Contains("range")) all = 1;
   else if (opt.Contains("base"))  all = 2;
   else                            return;

   Int_t bin, binx, biny, binz;
   Int_t firstx=0,lastx=0,firsty=0,lasty=0,firstz=0,lastz=0;
   if (all == 0) {
      lastx  = fXaxis.GetNbins()+1;
      if (fDimension > 1) lasty = fYaxis.GetNbins()+1;
      if (fDimension > 2) lastz = fZaxis.GetNbins()+1;
   } else {
      firstx = fXaxis.GetFirst(); lastx  = fXaxis.GetLast();
      if (fDimension > 1) {firsty = fYaxis.GetFirst(); lasty = fYaxis.GetLast();}
      if (fDimension > 2) {firstz = fZaxis.GetFirst(); lastz = fZaxis.GetLast();}
   }

   if (all== 2) {
      printf("          Title = %s\n", GetTitle());
      printf("          NbinsX= %d, xmin= %g, xmax=%g", fXaxis.GetNbins(), fXaxis.GetXmin(), fXaxis.GetXmax());
      if( fDimension > 1) printf(", NbinsY= %d, ymin= %g, ymax=%g", fYaxis.GetNbins(), fYaxis.GetXmin(), fYaxis.GetXmax());
      if( fDimension > 2) printf(", NbinsZ= %d, zmin= %g, zmax=%g", fZaxis.GetNbins(), fZaxis.GetXmin(), fZaxis.GetXmax());
      printf("\n");
      return;
   }

   Double_t w,e;
   Double_t x,y,z;
   if (fDimension == 1) {
      for (binx=firstx;binx<=lastx;binx++) {
         x = fXaxis.GetBinCenter(binx);
         w = RetrieveBinContent(binx);
         e = GetBinError(binx);
         if(fSumw2.fN) printf(" fSumw[%d]=%g, x=%g, error=%g\n",binx,w,x,e);
         else          printf(" fSumw[%d]=%g, x=%g\n",binx,w,x);
      }
   }
   if (fDimension == 2) {
      for (biny=firsty;biny<=lasty;biny++) {
         y = fYaxis.GetBinCenter(biny);
         for (binx=firstx;binx<=lastx;binx++) {
            bin = GetBin(binx,biny);
            x = fXaxis.GetBinCenter(binx);
            w = RetrieveBinContent(bin);
            e = GetBinError(bin);
            if(fSumw2.fN) printf(" fSumw[%d][%d]=%g, x=%g, y=%g, error=%g\n",binx,biny,w,x,y,e);
            else          printf(" fSumw[%d][%d]=%g, x=%g, y=%g\n",binx,biny,w,x,y);
         }
      }
   }
   if (fDimension == 3) {
      for (binz=firstz;binz<=lastz;binz++) {
         z = fZaxis.GetBinCenter(binz);
         for (biny=firsty;biny<=lasty;biny++) {
            y = fYaxis.GetBinCenter(biny);
            for (binx=firstx;binx<=lastx;binx++) {
               bin = GetBin(binx,biny,binz);
               x = fXaxis.GetBinCenter(binx);
               w = RetrieveBinContent(bin);
               e = GetBinError(bin);
               if(fSumw2.fN) printf(" fSumw[%d][%d][%d]=%g, x=%g, y=%g, z=%g, error=%g\n",binx,biny,binz,w,x,y,z,e);
               else          printf(" fSumw[%d][%d][%d]=%g, x=%g, y=%g, z=%g\n",binx,biny,binz,w,x,y,z);
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Using the current bin info, recompute the arrays for contents and errors

void TH1::Rebuild(Option_t *)
{
   SetBinsLength();
   if (fSumw2.fN) {
      fSumw2.Set(fNcells);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset this histogram: contents, errors, etc.
/// \param[in] option
///   - if "ICE" is specified, resets only Integral, Contents and Errors.
///   - if "ICES" is specified, resets only Integral, Contents, Errors and Statistics
///     This option is used
///   - if "M" is specified, resets also Minimum and Maximum

void TH1::Reset(Option_t *option)
{
   // The option "ICE" is used when extending the histogram (in ExtendAxis, LabelInflate, etc..)
   // The option "ICES is used in combination with the buffer (see BufferEmpty and BufferFill)

   TString opt = option;
   opt.ToUpper();
   fSumw2.Reset();
   if (fIntegral) {delete [] fIntegral; fIntegral = 0;}

   if (opt.Contains("M")) {
      SetMinimum();
      SetMaximum();
   }

   if (opt.Contains("ICE") && !opt.Contains("S")) return;

   // Setting fBuffer[0] = 0 is like resetting the buffer but not deleting it
   // But what is the sense of calling BufferEmpty() ? For making the axes ?
   // BufferEmpty will update contents that later will be
   // reset in calling TH1D::Reset. For this we need to reset the stats afterwards
   // It may be needed for computing the axis limits....
   if (fBuffer) {BufferEmpty(); fBuffer[0] = 0;}

   // need to reset also the statistics
   // (needs to be done after calling BufferEmpty() )
   fTsumw       = 0;
   fTsumw2      = 0;
   fTsumwx      = 0;
   fTsumwx2     = 0;
   fEntries     = 0;

   if (opt == "ICES") return;


   TObject *stats = fFunctions->FindObject("stats");
   fFunctions->Remove(stats);
   //special logic to support the case where the same object is
   //added multiple times in fFunctions.
   //This case happens when the same object is added with different
   //drawing modes
   TObject *obj;
   while ((obj  = fFunctions->First())) {
      while(fFunctions->Remove(obj)) { }
      delete obj;
   }
   if(stats) fFunctions->Add(stats);
   fContour.Set(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TH1::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   // empty the buffer before if it exists
   if (fBuffer) BufferEmpty();

   Bool_t nonEqiX = kFALSE;
   Bool_t nonEqiY = kFALSE;
   Bool_t nonEqiZ = kFALSE;
   Int_t i;
   static Int_t nxaxis = 0;
   static Int_t nyaxis = 0;
   static Int_t nzaxis = 0;
   TString sxaxis="xAxis",syaxis="yAxis",szaxis="zAxis";

   // Check if the histogram has equidistant X bins or not.  If not, we
   // create an array holding the bins.
   if (GetXaxis()->GetXbins()->fN && GetXaxis()->GetXbins()->fArray) {
      nonEqiX = kTRUE;
      nxaxis++;
      sxaxis += nxaxis;
      out << "   Double_t "<<sxaxis<<"[" << GetXaxis()->GetXbins()->fN
         << "] = {";
      for (i = 0; i < GetXaxis()->GetXbins()->fN; i++) {
         if (i != 0) out << ", ";
         out << GetXaxis()->GetXbins()->fArray[i];
      }
      out << "}; " << std::endl;
   }
   // If the histogram is 2 or 3 dimensional, check if the histogram
   // has equidistant Y bins or not.  If not, we create an array
   // holding the bins.
   if (fDimension > 1 && GetYaxis()->GetXbins()->fN &&
      GetYaxis()->GetXbins()->fArray) {
         nonEqiY = kTRUE;
         nyaxis++;
         syaxis += nyaxis;
         out << "   Double_t "<<syaxis<<"[" << GetYaxis()->GetXbins()->fN
            << "] = {";
         for (i = 0; i < GetYaxis()->GetXbins()->fN; i++) {
            if (i != 0) out << ", ";
            out << GetYaxis()->GetXbins()->fArray[i];
         }
         out << "}; " << std::endl;
   }
   // IF the histogram is 3 dimensional, check if the histogram
   // has equidistant Z bins or not.  If not, we create an array
   // holding the bins.
   if (fDimension > 2 && GetZaxis()->GetXbins()->fN &&
      GetZaxis()->GetXbins()->fArray) {
         nonEqiZ = kTRUE;
         nzaxis++;
         szaxis += nzaxis;
         out << "   Double_t "<<szaxis<<"[" << GetZaxis()->GetXbins()->fN
            << "] = {";
         for (i = 0; i < GetZaxis()->GetXbins()->fN; i++) {
            if (i != 0) out << ", ";
            out << GetZaxis()->GetXbins()->fArray[i];
         }
         out << "}; " << std::endl;
   }

   char quote = '"';
   out <<"   "<<std::endl;
   out <<"   "<< ClassName() <<" *";

   // Histogram pointer has by default the histogram name with an incremental suffix.
   // If the histogram belongs to a graph or a stack the suffix is not added because
   // the graph and stack objects are not aware of this new name. Same thing if
   // the histogram is drawn with the option COLZ because the TPaletteAxis drawn
   // when this option is selected, does not know this new name either.
   TString opt = option;
   opt.ToLower();
   static Int_t hcounter = 0;
   TString histName = GetName();
   if (     !histName.Contains("Graph")
         && !histName.Contains("_stack_")
         && !opt.Contains("colz")) {
      hcounter++;
      histName += "__";
      histName += hcounter;
   }
   histName = gInterpreter-> MapCppName(histName);
   const char *hname = histName.Data();
   if (!strlen(hname)) hname = "unnamed";
   TString savedName = GetName();
   this->SetName(hname);
   TString t(GetTitle());
   t.ReplaceAll("\\","\\\\");
   t.ReplaceAll("\"","\\\"");
   out << hname << " = new " << ClassName() << "(" << quote
       << hname << quote << "," << quote<< t.Data() << quote
       << "," << GetXaxis()->GetNbins();
   if (nonEqiX)
      out << ", "<<sxaxis;
   else
      out << "," << GetXaxis()->GetXmin()
      << "," << GetXaxis()->GetXmax();
   if (fDimension > 1) {
      out << "," << GetYaxis()->GetNbins();
      if (nonEqiY)
         out << ", "<<syaxis;
      else
         out << "," << GetYaxis()->GetXmin()
         << "," << GetYaxis()->GetXmax();
   }
   if (fDimension > 2) {
      out << "," << GetZaxis()->GetNbins();
      if (nonEqiZ)
         out << ", "<<szaxis;
      else
         out << "," << GetZaxis()->GetXmin()
         << "," << GetZaxis()->GetXmax();
   }
   out << ");" << std::endl;

   // save bin contents
   Int_t bin;
   for (bin=0;bin<fNcells;bin++) {
      Double_t bc = RetrieveBinContent(bin);
      if (bc) {
         out<<"   "<<hname<<"->SetBinContent("<<bin<<","<<bc<<");"<<std::endl;
      }
   }

   // save bin errors
   if (fSumw2.fN) {
      for (bin=0;bin<fNcells;bin++) {
         Double_t be = GetBinError(bin);
         if (be) {
            out<<"   "<<hname<<"->SetBinError("<<bin<<","<<be<<");"<<std::endl;
         }
      }
   }

   TH1::SavePrimitiveHelp(out, hname, option);
   this->SetName(savedName.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Helper function for the SavePrimitive functions from TH1
/// or classes derived from TH1, eg TProfile, TProfile2D.

void TH1::SavePrimitiveHelp(std::ostream &out, const char *hname, Option_t *option /*= ""*/)
{
   char quote = '"';
   if (TMath::Abs(GetBarOffset()) > 1e-5) {
      out<<"   "<<hname<<"->SetBarOffset("<<GetBarOffset()<<");"<<std::endl;
   }
   if (TMath::Abs(GetBarWidth()-1) > 1e-5) {
      out<<"   "<<hname<<"->SetBarWidth("<<GetBarWidth()<<");"<<std::endl;
   }
   if (fMinimum != -1111) {
      out<<"   "<<hname<<"->SetMinimum("<<fMinimum<<");"<<std::endl;
   }
   if (fMaximum != -1111) {
      out<<"   "<<hname<<"->SetMaximum("<<fMaximum<<");"<<std::endl;
   }
   if (fNormFactor != 0) {
      out<<"   "<<hname<<"->SetNormFactor("<<fNormFactor<<");"<<std::endl;
   }
   if (fEntries != 0) {
      out<<"   "<<hname<<"->SetEntries("<<fEntries<<");"<<std::endl;
   }
   if (fDirectory == 0) {
      out<<"   "<<hname<<"->SetDirectory(0);"<<std::endl;
   }
   if (TestBit(kNoStats)) {
      out<<"   "<<hname<<"->SetStats(0);"<<std::endl;
   }
   if (fOption.Length() != 0) {
      out<<"   "<<hname<<"->SetOption("<<quote<<fOption.Data()<<quote<<");"<<std::endl;
   }

   // save contour levels
   Int_t ncontours = GetContour();
   if (ncontours > 0) {
      out<<"   "<<hname<<"->SetContour("<<ncontours<<");"<<std::endl;
      Double_t zlevel;
      for (Int_t bin=0;bin<ncontours;bin++) {
         if (gPad->GetLogz()) {
            zlevel = TMath::Power(10,GetContourLevel(bin));
         } else {
            zlevel = GetContourLevel(bin);
         }
         out<<"   "<<hname<<"->SetContourLevel("<<bin<<","<<zlevel<<");"<<std::endl;
      }
   }

   // save list of functions
   TObjOptLink *lnk = (TObjOptLink*)fFunctions->FirstLink();
   TObject *obj;
   static Int_t funcNumber = 0;
   while (lnk) {
      obj = lnk->GetObject();
      obj->SavePrimitive(out,Form("nodraw #%d\n",++funcNumber));
      if (obj->InheritsFrom(TF1::Class())) {
         TString fname;
         fname.Form("%s%d",obj->GetName(),funcNumber);
         out << "   " << fname << "->SetParent(" << hname << ");\n";
         out<<"   "<<hname<<"->GetListOfFunctions()->Add("
            << fname <<");"<<std::endl;
      } else if (obj->InheritsFrom("TPaveStats")) {
         out<<"   "<<hname<<"->GetListOfFunctions()->Add(ptstats);"<<std::endl;
         out<<"   ptstats->SetParent("<<hname<<");"<<std::endl;
      } else if (obj->InheritsFrom("TPolyMarker")) {
         out<<"   "<<hname<<"->GetListOfFunctions()->Add("
            <<"pmarker ,"<<quote<<lnk->GetOption()<<quote<<");"<<std::endl;
      } else {
         out<<"   "<<hname<<"->GetListOfFunctions()->Add("
            <<obj->GetName()
            <<","<<quote<<lnk->GetOption()<<quote<<");"<<std::endl;
      }
      lnk = (TObjOptLink*)lnk->Next();
   }

   // save attributes
   SaveFillAttributes(out,hname,0,1001);
   SaveLineAttributes(out,hname,1,1,1);
   SaveMarkerAttributes(out,hname,1,1,1);
   fXaxis.SaveAttributes(out,hname,"->GetXaxis()");
   fYaxis.SaveAttributes(out,hname,"->GetYaxis()");
   fZaxis.SaveAttributes(out,hname,"->GetZaxis()");
   TString opt = option;
   opt.ToLower();
   if (!opt.Contains("nodraw")) {
      out<<"   "<<hname<<"->Draw("
         <<quote<<option<<quote<<");"<<std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy current attributes from/to current style

void TH1::UseCurrentStyle()
{
   if (!gStyle) return;
   if (gStyle->IsReading()) {
      fXaxis.ResetAttAxis("X");
      fYaxis.ResetAttAxis("Y");
      fZaxis.ResetAttAxis("Z");
      SetBarOffset(gStyle->GetBarOffset());
      SetBarWidth(gStyle->GetBarWidth());
      SetFillColor(gStyle->GetHistFillColor());
      SetFillStyle(gStyle->GetHistFillStyle());
      SetLineColor(gStyle->GetHistLineColor());
      SetLineStyle(gStyle->GetHistLineStyle());
      SetLineWidth(gStyle->GetHistLineWidth());
      SetMarkerColor(gStyle->GetMarkerColor());
      SetMarkerStyle(gStyle->GetMarkerStyle());
      SetMarkerSize(gStyle->GetMarkerSize());
      Int_t dostat = gStyle->GetOptStat();
      if (gStyle->GetOptFit() && !dostat) dostat = 1000000001;
      SetStats(dostat);
   } else {
      gStyle->SetBarOffset(fBarOffset);
      gStyle->SetBarWidth(fBarWidth);
      gStyle->SetHistFillColor(GetFillColor());
      gStyle->SetHistFillStyle(GetFillStyle());
      gStyle->SetHistLineColor(GetLineColor());
      gStyle->SetHistLineStyle(GetLineStyle());
      gStyle->SetHistLineWidth(GetLineWidth());
      gStyle->SetMarkerColor(GetMarkerColor());
      gStyle->SetMarkerStyle(GetMarkerStyle());
      gStyle->SetMarkerSize(GetMarkerSize());
      gStyle->SetOptStat(TestBit(kNoStats));
   }
   TIter next(GetListOfFunctions());
   TObject *obj;

   while ((obj = next())) {
      obj->UseCurrentStyle();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// For axis = 1,2 or 3 returns the mean value of the histogram along
/// X,Y or Z axis.
///
/// For axis = 11, 12, 13 returns the standard error of the mean value
/// of the histogram along X, Y or Z axis
///
/// Note that the mean value/StdDev is computed using the bins in the currently
/// defined range (see TAxis::SetRange). By default the range includes
/// all bins from 1 to nbins included, excluding underflows and overflows.
/// To force the underflows and overflows in the computation, one must
/// call the static function TH1::StatOverflows(kTRUE) before filling
/// the histogram.
///
/// IMPORTANT NOTE: The returned value depends on how the histogram statistics
/// are calculated. By default, if no range has been set, the returned mean is
/// the (unbinned) one calculated at fill time. If a range has been set, however,
/// the mean is calculated using the bins in range, as described above; THIS
/// IS TRUE EVEN IF THE RANGE INCLUDES ALL BINS--use TAxis::SetRange(0, 0) to unset
/// the range. To ensure that the returned mean (and all other statistics) is
/// always that of the binned data stored in the histogram, call TH1::ResetStats.
/// See TH1::GetStats.
///
/// Return mean value of this histogram along the X axis.
///
/// Note that the mean value/StdDev is computed using the bins in the currently
/// defined range (see TAxis::SetRange). By default the range includes
/// all bins from 1 to nbins included, excluding underflows and overflows.
/// To force the underflows and overflows in the computation, one must
/// call the static function TH1::StatOverflows(kTRUE) before filling
/// the histogram.

Double_t TH1::GetMean(Int_t axis) const
{
   if (axis<1 || (axis>3 && axis<11) || axis>13) return 0;
   Double_t stats[kNstat];
   for (Int_t i=4;i<kNstat;i++) stats[i] = 0;
   GetStats(stats);
   if (stats[0] == 0) return 0;
   if (axis<4){
      Int_t ax[3] = {2,4,7};
      return stats[ax[axis-1]]/stats[0];
   } else {
      // mean error = StdDev / sqrt( Neff )
      Double_t stddev = GetStdDev(axis-10);
      Double_t neff = GetEffectiveEntries();
      return ( neff > 0 ? stddev/TMath::Sqrt(neff) : 0. );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return standard error of mean of this histogram along the X axis.
///
/// Note that the mean value/StdDev is computed using the bins in the currently
/// defined range (see TAxis::SetRange). By default the range includes
/// all bins from 1 to nbins included, excluding underflows and overflows.
/// To force the underflows and overflows in the computation, one must
/// call the static function TH1::StatOverflows(kTRUE) before filling
/// the histogram.
///
/// Also note, that although the definition of standard error doesn't include the
/// assumption of normality, many uses of this feature implicitly assume it.
///
/// IMPORTANT NOTE: The returned value depends on how the histogram statistics
/// are calculated. By default, if no range has been set, the returned value is
/// the (unbinned) one calculated at fill time. If a range has been set, however,
/// the value is calculated using the bins in range, as described above; THIS
/// IS TRUE EVEN IF THE RANGE INCLUDES ALL BINS--use TAxis::SetRange(0, 0) to unset
/// the range. To ensure that the returned value (and all other statistics) is
/// always that of the binned data stored in the histogram, call TH1::ResetStats.
/// See TH1::GetStats.

Double_t TH1::GetMeanError(Int_t axis) const
{
   return GetMean(axis+10);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the Standard Deviation (Sigma).
/// The Sigma estimate is computed as
/// \f[
/// \sqrt{\frac{1}{N}(\sum(x_i-x_{mean})^2)}
/// \f]
/// For axis = 1,2 or 3 returns the Sigma value of the histogram along
/// X, Y or Z axis
/// For axis = 11, 12 or 13 returns the error of StdDev estimation along
/// X, Y or Z axis for Normal distribution
///
/// Note that the mean value/sigma is computed using the bins in the currently
/// defined range (see TAxis::SetRange). By default the range includes
/// all bins from 1 to nbins included, excluding underflows and overflows.
/// To force the underflows and overflows in the computation, one must
/// call the static function TH1::StatOverflows(kTRUE) before filling
/// the histogram.
///
/// IMPORTANT NOTE: The returned value depends on how the histogram statistics
/// are calculated. By default, if no range has been set, the returned standard
/// deviation is the (unbinned) one calculated at fill time. If a range has been
/// set, however, the standard deviation is calculated using the bins in range,
/// as described above; THIS IS TRUE EVEN IF THE RANGE INCLUDES ALL BINS--use
/// TAxis::SetRange(0, 0) to unset the range. To ensure that the returned standard
/// deviation (and all other statistics) is always that of the binned data stored
/// in the histogram, call TH1::ResetStats. See TH1::GetStats.

Double_t TH1::GetStdDev(Int_t axis) const
{
   if (axis<1 || (axis>3 && axis<11) || axis>13) return 0;

   Double_t x, stddev2, stats[kNstat];
   for (Int_t i=4;i<kNstat;i++) stats[i] = 0;
   GetStats(stats);
   if (stats[0] == 0) return 0;
   Int_t ax[3] = {2,4,7};
   Int_t axm = ax[axis%10 - 1];
   x    = stats[axm]/stats[0];
   // for negative stddev (e.g. when having negative weights) - return stdev=0
   stddev2 = TMath::Max( stats[axm+1]/stats[0] -x*x, 0.0 );
   if (axis<10)
      return TMath::Sqrt(stddev2);
   else {
      // The right formula for StdDev error depends on 4th momentum (see Kendall-Stuart Vol 1 pag 243)
      // formula valid for only gaussian distribution ( 4-th momentum =  3 * sigma^4 )
      Double_t neff = GetEffectiveEntries();
      return ( neff > 0 ? TMath::Sqrt(stddev2/(2*neff) ) : 0. );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return error of standard deviation estimation for Normal distribution
///
/// Note that the mean value/StdDev is computed using the bins in the currently
/// defined range (see TAxis::SetRange). By default the range includes
/// all bins from 1 to nbins included, excluding underflows and overflows.
/// To force the underflows and overflows in the computation, one must
/// call the static function TH1::StatOverflows(kTRUE) before filling
/// the histogram.
///
/// Value returned is standard deviation of sample standard deviation.
/// Note that it is an approximated value which is valid only in the case that the
/// original data distribution is Normal. The correct one would require
/// the 4-th momentum value, which cannot be accurately estimated from an histogram since
/// the x-information for all entries is not kept.
///
/// IMPORTANT NOTE: The returned value depends on how the histogram statistics
/// are calculated. By default, if no range has been set, the returned value is
/// the (unbinned) one calculated at fill time. If a range has been set, however,
/// the value is calculated using the bins in range, as described above; THIS
/// IS TRUE EVEN IF THE RANGE INCLUDES ALL BINS--use TAxis::SetRange(0, 0) to unset
/// the range. To ensure that the returned value (and all other statistics) is
/// always that of the binned data stored in the histogram, call TH1::ResetStats.
/// See TH1::GetStats.

Double_t TH1::GetStdDevError(Int_t axis) const
{
   return GetStdDev(axis+10);
}

////////////////////////////////////////////////////////////////////////////////
///  - For axis = 1, 2 or 3 returns skewness of the histogram along x, y or z axis.
///  - For axis = 11, 12 or 13 returns the approximate standard error of skewness
///    of the histogram along x, y or z axis
///
///Note, that since third and fourth moment are not calculated
///at the fill time, skewness and its standard error are computed bin by bin
///
/// IMPORTANT NOTE: The returned value depends on how the histogram statistics
/// are calculated. See TH1::GetMean and TH1::GetStdDev.

Double_t TH1::GetSkewness(Int_t axis) const
{

   if (axis > 0 && axis <= 3){

      Double_t mean = GetMean(axis);
      Double_t stddev = GetStdDev(axis);
      Double_t stddev3 = stddev*stddev*stddev;

      Int_t firstBinX = fXaxis.GetFirst();
      Int_t lastBinX  = fXaxis.GetLast();
      Int_t firstBinY = fYaxis.GetFirst();
      Int_t lastBinY  = fYaxis.GetLast();
      Int_t firstBinZ = fZaxis.GetFirst();
      Int_t lastBinZ  = fZaxis.GetLast();
      // include underflow/overflow if TH1::StatOverflows(kTRUE) in case no range is set on the axis
      if (GetStatOverflowsBehaviour()) {
        if ( !fXaxis.TestBit(TAxis::kAxisRange) ) {
            if (firstBinX == 1) firstBinX = 0;
            if (lastBinX ==  fXaxis.GetNbins() ) lastBinX += 1;
         }
         if ( !fYaxis.TestBit(TAxis::kAxisRange) ) {
            if (firstBinY == 1) firstBinY = 0;
            if (lastBinY ==  fYaxis.GetNbins() ) lastBinY += 1;
         }
         if ( !fZaxis.TestBit(TAxis::kAxisRange) ) {
            if (firstBinZ == 1) firstBinZ = 0;
            if (lastBinZ ==  fZaxis.GetNbins() ) lastBinZ += 1;
         }
      }

      Double_t x = 0;
      Double_t sum=0;
      Double_t np=0;
      for (Int_t  binx = firstBinX; binx <= lastBinX; binx++) {
         for (Int_t biny = firstBinY; biny <= lastBinY; biny++) {
            for (Int_t binz = firstBinZ; binz <= lastBinZ; binz++) {
               if (axis==1 ) x = fXaxis.GetBinCenter(binx);
               else if (axis==2 ) x = fYaxis.GetBinCenter(biny);
               else if (axis==3 ) x = fZaxis.GetBinCenter(binz);
               Double_t w = GetBinContent(binx,biny,binz);
               np+=w;
               sum+=w*(x-mean)*(x-mean)*(x-mean);
            }
         }
      }
      sum/=np*stddev3;
      return sum;
   }
   else if (axis > 10 && axis <= 13) {
      //compute standard error of skewness
      // assume parent normal distribution use formula from  Kendall-Stuart, Vol 1 pag 243, second edition
      Double_t neff = GetEffectiveEntries();
      return ( neff > 0 ? TMath::Sqrt(6./neff ) : 0. );
   }
   else {
      Error("GetSkewness", "illegal value of parameter");
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
///  - For axis =1, 2 or 3 returns kurtosis of the histogram along x, y or z axis.
///    Kurtosis(gaussian(0, 1)) = 0.
///  - For axis =11, 12 or 13 returns the approximate standard error of kurtosis
///    of the histogram along x, y or z axis
////
/// Note, that since third and fourth moment are not calculated
/// at the fill time, kurtosis and its standard error are computed bin by bin
///
/// IMPORTANT NOTE: The returned value depends on how the histogram statistics
/// are calculated. See TH1::GetMean and TH1::GetStdDev.

Double_t TH1::GetKurtosis(Int_t axis) const
{
   if (axis > 0 && axis <= 3){

      Double_t mean = GetMean(axis);
      Double_t stddev = GetStdDev(axis);
      Double_t stddev4 = stddev*stddev*stddev*stddev;

      Int_t firstBinX = fXaxis.GetFirst();
      Int_t lastBinX  = fXaxis.GetLast();
      Int_t firstBinY = fYaxis.GetFirst();
      Int_t lastBinY  = fYaxis.GetLast();
      Int_t firstBinZ = fZaxis.GetFirst();
      Int_t lastBinZ  = fZaxis.GetLast();
      // include underflow/overflow if TH1::StatOverflows(kTRUE) in case no range is set on the axis
      if (GetStatOverflowsBehaviour()) {
        if ( !fXaxis.TestBit(TAxis::kAxisRange) ) {
            if (firstBinX == 1) firstBinX = 0;
            if (lastBinX ==  fXaxis.GetNbins() ) lastBinX += 1;
         }
         if ( !fYaxis.TestBit(TAxis::kAxisRange) ) {
            if (firstBinY == 1) firstBinY = 0;
            if (lastBinY ==  fYaxis.GetNbins() ) lastBinY += 1;
         }
         if ( !fZaxis.TestBit(TAxis::kAxisRange) ) {
            if (firstBinZ == 1) firstBinZ = 0;
            if (lastBinZ ==  fZaxis.GetNbins() ) lastBinZ += 1;
         }
      }

      Double_t x = 0;
      Double_t sum=0;
      Double_t np=0;
      for (Int_t binx = firstBinX; binx <= lastBinX; binx++) {
         for (Int_t biny = firstBinY; biny <= lastBinY; biny++) {
            for (Int_t binz = firstBinZ; binz <= lastBinZ; binz++) {
               if (axis==1 ) x = fXaxis.GetBinCenter(binx);
               else if (axis==2 ) x = fYaxis.GetBinCenter(biny);
               else if (axis==3 ) x = fZaxis.GetBinCenter(binz);
               Double_t w = GetBinContent(binx,biny,binz);
               np+=w;
               sum+=w*(x-mean)*(x-mean)*(x-mean)*(x-mean);
            }
         }
      }
      sum/=(np*stddev4);
      return sum-3;

   } else if (axis > 10 && axis <= 13) {
      //compute standard error of skewness
      // assume parent normal distribution use formula from  Kendall-Stuart, Vol 1 pag 243, second edition
      Double_t neff = GetEffectiveEntries();
      return ( neff > 0 ? TMath::Sqrt(24./neff ) : 0. );
   }
   else {
      Error("GetKurtosis", "illegal value of parameter");
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// fill the array stats from the contents of this histogram
/// The array stats must be correctly dimensioned in the calling program.
///
/// ~~~ {.cpp}
///      stats[0] = sumw
///      stats[1] = sumw2
///      stats[2] = sumwx
///      stats[3] = sumwx2
/// ~~~
///
/// If no axis-subrange is specified (via TAxis::SetRange), the array stats
/// is simply a copy of the statistics quantities computed at filling time.
/// If a sub-range is specified, the function recomputes these quantities
/// from the bin contents in the current axis range.
///
/// IMPORTANT NOTE: This means that the returned statistics are context-dependent.
/// If TAxis::kAxisRange, the returned statistics are dependent on the binning;
/// otherwise, they are a copy of the histogram statistics computed at fill time,
/// which are unbinned by default (calling TH1::ResetStats forces them to use
/// binned statistics). You can reset TAxis::kAxisRange using TAxis::SetRange(0, 0).
///
/// Note that the mean value/StdDev is computed using the bins in the currently
/// defined range (see TAxis::SetRange). By default the range includes
/// all bins from 1 to nbins included, excluding underflows and overflows.
/// To force the underflows and overflows in the computation, one must
/// call the static function TH1::StatOverflows(kTRUE) before filling
/// the histogram.

void TH1::GetStats(Double_t *stats) const
{
   if (fBuffer) ((TH1*)this)->BufferEmpty();

   // Loop on bins (possibly including underflows/overflows)
   Int_t bin, binx;
   Double_t w,err;
   Double_t x;
   // identify the case of labels with extension of axis range
   // in this case the statistics in x does not make any sense
   Bool_t labelHist =  ((const_cast<TAxis&>(fXaxis)).GetLabels() && fXaxis.CanExtend() );
   // fTsumw == 0 && fEntries > 0 is a special case when uses SetBinContent or calls ResetStats before
   if ( (fTsumw == 0 && fEntries > 0) || fXaxis.TestBit(TAxis::kAxisRange) ) {
      for (bin=0;bin<4;bin++) stats[bin] = 0;

      Int_t firstBinX = fXaxis.GetFirst();
      Int_t lastBinX  = fXaxis.GetLast();
      // include underflow/overflow if TH1::StatOverflows(kTRUE) in case no range is set on the axis
      if (GetStatOverflowsBehaviour() && !fXaxis.TestBit(TAxis::kAxisRange)) {
         if (firstBinX == 1) firstBinX = 0;
         if (lastBinX ==  fXaxis.GetNbins() ) lastBinX += 1;
      }
      for (binx = firstBinX; binx <= lastBinX; binx++) {
         x   = fXaxis.GetBinCenter(binx);
         //w   = TMath::Abs(RetrieveBinContent(binx));
         // not sure what to do here if w < 0
         w   = RetrieveBinContent(binx);
         err = TMath::Abs(GetBinError(binx));
         stats[0] += w;
         stats[1] += err*err;
         // statistics in x makes sense only for not labels histograms
         if (!labelHist)  {
            stats[2] += w*x;
            stats[3] += w*x*x;
         }
      }
      // if (stats[0] < 0) {
      //    // in case total is negative do something ??
      //    stats[0] = 0;
      // }
   } else {
      stats[0] = fTsumw;
      stats[1] = fTsumw2;
      stats[2] = fTsumwx;
      stats[3] = fTsumwx2;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Replace current statistics with the values in array stats

void TH1::PutStats(Double_t *stats)
{
   fTsumw   = stats[0];
   fTsumw2  = stats[1];
   fTsumwx  = stats[2];
   fTsumwx2 = stats[3];
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the statistics including the number of entries
/// and replace with values calculated from bin content
///
/// The number of entries is set to the total bin content or (in case of weighted histogram)
/// to number of effective entries
///
/// Note that, by default, before calling this function, statistics are those
/// computed at fill time, which are unbinned. See TH1::GetStats.

void TH1::ResetStats()
{
   Double_t stats[kNstat] = {0};
   fTsumw = 0;
   fEntries = 1; // to force re-calculation of the statistics in TH1::GetStats
   GetStats(stats);
   PutStats(stats);
   fEntries = TMath::Abs(fTsumw);
   // use effective entries for weighted histograms:  (sum_w) ^2 / sum_w2
   if (fSumw2.fN > 0 && fTsumw > 0 && stats[1] > 0 ) fEntries = stats[0]*stats[0]/ stats[1];
}

////////////////////////////////////////////////////////////////////////////////
/// Return the sum of weights excluding under/overflows.

Double_t TH1::GetSumOfWeights() const
{
   if (fBuffer) const_cast<TH1*>(this)->BufferEmpty();

   Int_t bin,binx,biny,binz;
   Double_t sum =0;
   for(binz=1; binz<=fZaxis.GetNbins(); binz++) {
      for(biny=1; biny<=fYaxis.GetNbins(); biny++) {
         for(binx=1; binx<=fXaxis.GetNbins(); binx++) {
            bin = GetBin(binx,biny,binz);
            sum += RetrieveBinContent(bin);
         }
      }
   }
   return sum;
}

////////////////////////////////////////////////////////////////////////////////
///Return integral of bin contents. Only bins in the bins range are considered.
///
/// By default the integral is computed as the sum of bin contents in the range.
/// if option "width" is specified, the integral is the sum of
/// the bin contents multiplied by the bin width in x.

Double_t TH1::Integral(Option_t *option) const
{
   return Integral(fXaxis.GetFirst(),fXaxis.GetLast(),option);
}

////////////////////////////////////////////////////////////////////////////////
/// Return integral of bin contents in range [binx1,binx2].
///
/// By default the integral is computed as the sum of bin contents in the range.
/// if option "width" is specified, the integral is the sum of
/// the bin contents multiplied by the bin width in x.

Double_t TH1::Integral(Int_t binx1, Int_t binx2, Option_t *option) const
{
   double err = 0;
   return DoIntegral(binx1,binx2,0,-1,0,-1,err,option);
}

////////////////////////////////////////////////////////////////////////////////
/// Return integral of bin contents in range [binx1,binx2] and its error.
///
/// By default the integral is computed as the sum of bin contents in the range.
/// if option "width" is specified, the integral is the sum of
/// the bin contents multiplied by the bin width in x.
/// the error is computed using error propagation from the bin errors assuming that
/// all the bins are uncorrelated

Double_t TH1::IntegralAndError(Int_t binx1, Int_t binx2, Double_t & error, Option_t *option) const
{
   return DoIntegral(binx1,binx2,0,-1,0,-1,error,option,kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Internal function compute integral and optionally the error  between the limits
/// specified by the bin number values working for all histograms (1D, 2D and 3D)

Double_t TH1::DoIntegral(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Double_t & error ,
                          Option_t *option, Bool_t doError) const
{
   if (fBuffer) ((TH1*)this)->BufferEmpty();

   Int_t nx = GetNbinsX() + 2;
   if (binx1 < 0) binx1 = 0;
   if (binx2 >= nx || binx2 < binx1) binx2 = nx - 1;

   if (GetDimension() > 1) {
      Int_t ny = GetNbinsY() + 2;
      if (biny1 < 0) biny1 = 0;
      if (biny2 >= ny || biny2 < biny1) biny2 = ny - 1;
   } else {
      biny1 = 0; biny2 = 0;
   }

   if (GetDimension() > 2) {
      Int_t nz = GetNbinsZ() + 2;
      if (binz1 < 0) binz1 = 0;
      if (binz2 >= nz || binz2 < binz1) binz2 = nz - 1;
   } else {
      binz1 = 0; binz2 = 0;
   }

   //   - Loop on bins in specified range
   TString opt = option;
   opt.ToLower();
   Bool_t width   = kFALSE;
   if (opt.Contains("width")) width = kTRUE;


   Double_t dx = 1., dy = .1, dz =.1;
   Double_t integral = 0;
   Double_t igerr2 = 0;
   for (Int_t binx = binx1; binx <= binx2; ++binx) {
      if (width) dx = fXaxis.GetBinWidth(binx);
      for (Int_t biny = biny1; biny <= biny2; ++biny) {
         if (width) dy = fYaxis.GetBinWidth(biny);
         for (Int_t binz = binz1; binz <= binz2; ++binz) {
            Int_t bin = GetBin(binx, biny, binz);
            Double_t dv = 0.0;
            if (width) {
               dz = fZaxis.GetBinWidth(binz);
               dv = dx * dy * dz;
               integral += RetrieveBinContent(bin) * dv;
            } else {
              integral += RetrieveBinContent(bin);
            }
            if (doError) {
               if (width)  igerr2 += GetBinErrorSqUnchecked(bin) * dv * dv;
               else        igerr2 += GetBinErrorSqUnchecked(bin);
            }
         }
      }
   }

   if (doError) error = TMath::Sqrt(igerr2);
   return integral;
}

////////////////////////////////////////////////////////////////////////////////
/// Statistical test of compatibility in shape between
/// this histogram and h2, using the Anderson-Darling 2 sample test.
///
/// The AD 2 sample test formula are derived from the paper
/// F.W Scholz, M.A. Stephens "k-Sample Anderson-Darling Test".
///
/// The test is implemented in root in the ROOT::Math::GoFTest class
/// It is the same formula ( (6) in the paper), and also shown in
/// [this preprint](http://arxiv.org/pdf/0804.0380v1.pdf)
///
/// Binned data are considered as un-binned data
/// with identical observation happening in the bin center.
///
/// \param[in] option is a character string to specify options
///    - "D" Put out a line of "Debug" printout
///    - "T" Return the normalized A-D test statistic
///
///  - Note1: Underflow and overflow are not considered in the test
///  - Note2:  The test works only for un-weighted histogram (i.e. representing counts)
///  - Note3:  The histograms are not required to have the same X axis
///  - Note4:  The test works only for 1-dimensional histograms

Double_t TH1::AndersonDarlingTest(const TH1 *h2, Option_t *option) const
{
   Double_t advalue = 0;
   Double_t pvalue = AndersonDarlingTest(h2, advalue);

   TString opt = option;
   opt.ToUpper();
   if (opt.Contains("D") ) {
      printf(" AndersonDarlingTest Prob     = %g, AD TestStatistic  = %g\n",pvalue,advalue);
   }
   if (opt.Contains("T") ) return advalue;

   return pvalue;
}

////////////////////////////////////////////////////////////////////////////////
/// Same function as above but returning also the test statistic value

Double_t TH1::AndersonDarlingTest(const TH1 *h2, Double_t & advalue) const
{
   if (GetDimension() != 1 || h2->GetDimension() != 1) {
      Error("AndersonDarlingTest","Histograms must be 1-D");
      return -1;
   }

   // empty the buffer. Probably we could add as an unbinned test
   if (fBuffer) ((TH1*)this)->BufferEmpty();

   // use the BinData class
   ROOT::Fit::BinData data1;
   ROOT::Fit::BinData data2;

   ROOT::Fit::FillData(data1, this, 0);
   ROOT::Fit::FillData(data2, h2, 0);

   double pvalue;
   ROOT::Math::GoFTest::AndersonDarling2SamplesTest(data1,data2, pvalue,advalue);

   return pvalue;
}

////////////////////////////////////////////////////////////////////////////////
/// Statistical test of compatibility in shape between
/// this histogram and h2, using Kolmogorov test.
/// Note that the KolmogorovTest (KS) test should in theory be used only for unbinned data
/// and not for binned data as in the case of the histogram (see NOTE 3 below).
/// So, before using this method blindly, read the NOTE 3.
///
/// Default: Ignore under- and overflow bins in comparison
///
/// \param[in] h2 histogram
/// \param[in] option is a character string to specify options
///    - "U" include Underflows in test  (also for 2-dim)
///    -  "O" include Overflows     (also valid for 2-dim)
///    -  "N" include comparison of normalizations
///    -  "D" Put out a line of "Debug" printout
///    -  "M" Return the Maximum Kolmogorov distance instead of prob
///    -  "X" Run the pseudo experiments post-processor with the following procedure:
///       make pseudoexperiments based on random values from the parent distribution,
///       compare the KS distance of the pseudoexperiment to the parent
///       distribution, and count all the KS values above the value
///       obtained from the original data to Monte Carlo distribution.
///       The number of pseudo-experiments nEXPT is currently fixed at 1000.
///       The function returns the probability.
///       (thanks to Ben Kilminster to submit this procedure). Note that
///       this option "X" is much slower.
///
/// The returned function value is the probability of test
///   (much less than one means NOT compatible)
///
/// Code adapted by Rene Brun from original HBOOK routine HDIFF
///
/// NOTE1
/// A good description of the Kolmogorov test can be seen at:
/// http://www.itl.nist.gov/div898/handbook/eda/section3/eda35g.htm
///
/// NOTE2
/// see also alternative function TH1::Chi2Test
/// The Kolmogorov test is assumed to give better results than Chi2Test
/// in case of histograms with low statistics.
///
/// NOTE3 (Jan Conrad, Fred James)
/// "The returned value PROB is calculated such that it will be
/// uniformly distributed between zero and one for compatible histograms,
/// provided the data are not binned (or the number of bins is very large
/// compared with the number of events). Users who have access to unbinned
/// data and wish exact confidence levels should therefore not put their data
/// into histograms, but should call directly TMath::KolmogorovTest. On
/// the other hand, since TH1 is a convenient way of collecting data and
/// saving space, this function has been provided. However, the values of
/// PROB for binned data will be shifted slightly higher than expected,
/// depending on the effects of the binning. For example, when comparing two
/// uniform distributions of 500 events in 100 bins, the values of PROB,
/// instead of being exactly uniformly distributed between zero and one, have
/// a mean value of about 0.56. We can apply a useful
/// rule: As long as the bin width is small compared with any significant
/// physical effect (for example the experimental resolution) then the binning
/// cannot have an important effect. Therefore, we believe that for all
/// practical purposes, the probability value PROB is calculated correctly
/// provided the user is aware that:
///
///  1. The value of PROB should not be expected to have exactly the correct
///     distribution for binned data.
///  2. The user is responsible for seeing to it that the bin widths are
///     small compared with any physical phenomena of interest.
///  3. The effect of binning (if any) is always to make the value of PROB
///     slightly too big. That is, setting an acceptance criterion of (PROB>0.05
///     will assure that at most 5% of truly compatible histograms are rejected,
///     and usually somewhat less."
///
/// Note also that for GoF test of unbinned data ROOT provides also the class
/// ROOT::Math::GoFTest. The class has also method for doing one sample tests
/// (i.e. comparing the data with a given distribution).

Double_t TH1::KolmogorovTest(const TH1 *h2, Option_t *option) const
{
   TString opt = option;
   opt.ToUpper();

   Double_t prob = 0;
   TH1 *h1 = (TH1*)this;
   if (h2 == 0) return 0;
   const TAxis *axis1 = h1->GetXaxis();
   const TAxis *axis2 = h2->GetXaxis();
   Int_t ncx1   = axis1->GetNbins();
   Int_t ncx2   = axis2->GetNbins();

   // Check consistency of dimensions
   if (h1->GetDimension() != 1 || h2->GetDimension() != 1) {
      Error("KolmogorovTest","Histograms must be 1-D\n");
      return 0;
   }

   // Check consistency in number of channels
   if (ncx1 != ncx2) {
      Error("KolmogorovTest","Histograms have different number of bins, %d and %d\n",ncx1,ncx2);
      return 0;
   }

   // empty the buffer. Probably we could add as an unbinned test
   if (fBuffer) ((TH1*)this)->BufferEmpty();

   // Check consistency in bin edges
   for(Int_t i = 1; i <= axis1->GetNbins() + 1; ++i) {
      if(!TMath::AreEqualRel(axis1->GetBinLowEdge(i), axis2->GetBinLowEdge(i), 1.E-15)) {
         Error("KolmogorovTest","Histograms are not consistent: they have different bin edges");
         return 0;
      }
   }

   Bool_t afunc1 = kFALSE;
   Bool_t afunc2 = kFALSE;
   Double_t sum1 = 0, sum2 = 0;
   Double_t ew1, ew2, w1 = 0, w2 = 0;
   Int_t bin;
   Int_t ifirst = 1;
   Int_t ilast = ncx1;
   // integral of all bins (use underflow/overflow if option)
   if (opt.Contains("U")) ifirst = 0;
   if (opt.Contains("O")) ilast = ncx1 +1;
   for (bin = ifirst; bin <= ilast; bin++) {
      sum1 += h1->RetrieveBinContent(bin);
      sum2 += h2->RetrieveBinContent(bin);
      ew1   = h1->GetBinError(bin);
      ew2   = h2->GetBinError(bin);
      w1   += ew1*ew1;
      w2   += ew2*ew2;
   }
   if (sum1 == 0) {
      Error("KolmogorovTest","Histogram1 %s integral is zero\n",h1->GetName());
      return 0;
   }
   if (sum2 == 0) {
      Error("KolmogorovTest","Histogram2 %s integral is zero\n",h2->GetName());
      return 0;
   }

   // calculate the effective entries.
   // the case when errors are zero (w1 == 0 or w2 ==0) are equivalent to
   // compare to a function. In that case the rescaling is done only on sqrt(esum2) or sqrt(esum1)
   Double_t esum1 = 0, esum2 = 0;
   if (w1 > 0)
      esum1 = sum1 * sum1 / w1;
   else
      afunc1 = kTRUE;    // use later for calculating z

   if (w2 > 0)
      esum2 = sum2 * sum2 / w2;
   else
      afunc2 = kTRUE;    // use later for calculating z

   if (afunc2 && afunc1) {
      Error("KolmogorovTest","Errors are zero for both histograms\n");
      return 0;
   }


   Double_t s1 = 1/sum1;
   Double_t s2 = 1/sum2;

   // Find largest difference for Kolmogorov Test
   Double_t dfmax =0, rsum1 = 0, rsum2 = 0;

   for (bin=ifirst;bin<=ilast;bin++) {
      rsum1 += s1*h1->RetrieveBinContent(bin);
      rsum2 += s2*h2->RetrieveBinContent(bin);
      dfmax = TMath::Max(dfmax,TMath::Abs(rsum1-rsum2));
   }

   // Get Kolmogorov probability
   Double_t z, prb1=0, prb2=0, prb3=0;

   // case h1 is exact (has zero errors)
  if  (afunc1)
      z = dfmax*TMath::Sqrt(esum2);
  // case h2 has zero errors
  else if (afunc2)
      z = dfmax*TMath::Sqrt(esum1);
  else
     // for comparison between two data sets
     z = dfmax*TMath::Sqrt(esum1*esum2/(esum1+esum2));

   prob = TMath::KolmogorovProb(z);

   // option N to combine normalization makes sense if both afunc1 and afunc2 are false
   if (opt.Contains("N") && !(afunc1 || afunc2 ) ) {
      // Combine probabilities for shape and normalization,
      prb1 = prob;
      Double_t d12    = esum1-esum2;
      Double_t chi2   = d12*d12/(esum1+esum2);
      prb2 = TMath::Prob(chi2,1);
      // see Eadie et al., section 11.6.2
      if (prob > 0 && prb2 > 0) prob *= prb2*(1-TMath::Log(prob*prb2));
      else                      prob = 0;
   }
   // X option. Pseudo-experiments post-processor to determine KS probability
   const Int_t nEXPT = 1000;
   if (opt.Contains("X") && !(afunc1 || afunc2 ) ) {
      Double_t dSEXPT;
      TH1 *h1_cpy =  (TH1 *)(gDirectory ? gDirectory->CloneObject(this, kFALSE) : gROOT->CloneObject(this, kFALSE));
      TH1 *h1Expt = (TH1*)(gDirectory ? gDirectory->CloneObject(this,kFALSE) : gROOT->CloneObject(this,kFALSE));
      TH1 *h2Expt = (TH1*)(gDirectory ? gDirectory->CloneObject(this,kFALSE) : gROOT->CloneObject(this,kFALSE));

      if (GetMinimum() < 0.0) {
         // we need to create a new histogram
         // With negative bins we can't draw random samples in a meaningful way.
         Warning("KolmogorovTest", "Detected bins with negative weights, these have been ignored and output might be "
                                   "skewed. Reduce number of bins for histogram?");
         while (h1_cpy->GetMinimum() < 0.0) {
            Int_t idx = h1_cpy->GetMinimumBin();
            h1_cpy->SetBinContent(idx, 0.0);
         }
      }

      // make nEXPT experiments (this should be a parameter)
      prb3 = 0;
      for (Int_t i=0; i < nEXPT; i++) {
         h1Expt->Reset();
         h2Expt->Reset();
         h1Expt->FillRandom(h1_cpy, (Int_t)esum1);
         h2Expt->FillRandom(h1_cpy, (Int_t)esum2);
         dSEXPT = h1Expt->KolmogorovTest(h2Expt,"M");
         if (dSEXPT>dfmax) prb3 += 1.0;
      }
      prb3 /= (Double_t)nEXPT;
      delete h1_cpy;
      delete h1Expt;
      delete h2Expt;
   }

   // debug printout
   if (opt.Contains("D")) {
      printf(" Kolmo Prob  h1 = %s, sum bin content =%g  effective entries =%g\n",h1->GetName(),sum1,esum1);
      printf(" Kolmo Prob  h2 = %s, sum bin content =%g  effective entries =%g\n",h2->GetName(),sum2,esum2);
      printf(" Kolmo Prob     = %g, Max Dist = %g\n",prob,dfmax);
      if (opt.Contains("N"))
         printf(" Kolmo Prob     = %f for shape alone, =%f for normalisation alone\n",prb1,prb2);
      if (opt.Contains("X"))
         printf(" Kolmo Prob     = %f with %d pseudo-experiments\n",prb3,nEXPT);
   }
   // This numerical error condition should never occur:
   if (TMath::Abs(rsum1-1) > 0.002) Warning("KolmogorovTest","Numerical problems with h1=%s\n",h1->GetName());
   if (TMath::Abs(rsum2-1) > 0.002) Warning("KolmogorovTest","Numerical problems with h2=%s\n",h2->GetName());

   if(opt.Contains("M"))      return dfmax;
   else if(opt.Contains("X")) return prb3;
   else                       return prob;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace bin contents by the contents of array content

void TH1::SetContent(const Double_t *content)
{
   fEntries = fNcells;
   fTsumw = 0;
   for (Int_t i = 0; i < fNcells; ++i) UpdateBinContent(i, content[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Return contour values into array levels if pointer levels is non zero.
///
/// The function returns the number of contour levels.
/// see GetContourLevel to return one contour only

Int_t TH1::GetContour(Double_t *levels)
{
   Int_t nlevels = fContour.fN;
   if (levels) {
      if (nlevels == 0) {
         nlevels = 20;
         SetContour(nlevels);
      } else {
         if (TestBit(kUserContour) == 0) SetContour(nlevels);
      }
      for (Int_t level=0; level<nlevels; level++) levels[level] = fContour.fArray[level];
   }
   return nlevels;
}

////////////////////////////////////////////////////////////////////////////////
/// Return value of contour number level.
/// Use GetContour to return the array of all contour levels

Double_t TH1::GetContourLevel(Int_t level) const
{
   return (level >= 0 && level < fContour.fN) ? fContour.fArray[level] : 0.0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the value of contour number "level" in Pad coordinates.
/// ie: if the Pad is in log scale along Z it returns le log of the contour level
/// value. See GetContour to return the array of all contour levels

Double_t TH1::GetContourLevelPad(Int_t level) const
{
   if (level <0 || level >= fContour.fN) return 0;
   Double_t zlevel = fContour.fArray[level];

   // In case of user defined contours and Pad in log scale along Z,
   // fContour.fArray doesn't contain the log of the contour whereas it does
   // in case of equidistant contours.
   if (gPad && gPad->GetLogz() && TestBit(kUserContour)) {
      if (zlevel <= 0) return 0;
      zlevel = TMath::Log10(zlevel);
   }
   return zlevel;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the maximum number of entries to be kept in the buffer.

void TH1::SetBuffer(Int_t buffersize, Option_t * /*option*/)
{
   if (fBuffer) {
      BufferEmpty();
      delete [] fBuffer;
      fBuffer = 0;
   }
   if (buffersize <= 0) {
      fBufferSize = 0;
      return;
   }
   if (buffersize < 100) buffersize = 100;
   fBufferSize = 1 + buffersize*(fDimension+1);
   fBuffer = new Double_t[fBufferSize];
   memset(fBuffer,0,sizeof(Double_t)*fBufferSize);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the number and values of contour levels.
///
/// By default the number of contour levels is set to 20. The contours values
/// in the array "levels" should be specified in increasing order.
///
/// if argument levels = 0 or missing, equidistant contours are computed

void TH1::SetContour(Int_t  nlevels, const Double_t *levels)
{
   Int_t level;
   ResetBit(kUserContour);
   if (nlevels <=0 ) {
      fContour.Set(0);
      return;
   }
   fContour.Set(nlevels);

   //   -  Contour levels are specified
   if (levels) {
      SetBit(kUserContour);
      for (level=0; level<nlevels; level++) fContour.fArray[level] = levels[level];
   } else {
      //   - contour levels are computed automatically as equidistant contours
      Double_t zmin = GetMinimum();
      Double_t zmax = GetMaximum();
      if ((zmin == zmax) && (zmin != 0)) {
         zmax += 0.01*TMath::Abs(zmax);
         zmin -= 0.01*TMath::Abs(zmin);
      }
      Double_t dz   = (zmax-zmin)/Double_t(nlevels);
      if (gPad && gPad->GetLogz()) {
         if (zmax <= 0) return;
         if (zmin <= 0) zmin = 0.001*zmax;
         zmin = TMath::Log10(zmin);
         zmax = TMath::Log10(zmax);
         dz   = (zmax-zmin)/Double_t(nlevels);
      }
      for (level=0; level<nlevels; level++) {
         fContour.fArray[level] = zmin + dz*Double_t(level);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set value for one contour level.

void TH1::SetContourLevel(Int_t level, Double_t value)
{
   if (level < 0 || level >= fContour.fN) return;
   SetBit(kUserContour);
   fContour.fArray[level] = value;
}

////////////////////////////////////////////////////////////////////////////////
/// Return maximum value smaller than maxval of bins in the range,
/// unless the value has been overridden by TH1::SetMaximum,
/// in which case it returns that value. (This happens, for example,
/// when the histogram is drawn and the y or z axis limits are changed
///
/// To get the maximum value of bins in the histogram regardless of
/// whether the value has been overridden, use
///
/// ~~~ {.cpp}
///  h->GetBinContent(h->GetMaximumBin())
/// ~~~

Double_t TH1::GetMaximum(Double_t maxval) const
{
   if (fMaximum != -1111) return fMaximum;

   // empty the buffer
   if (fBuffer) ((TH1*)this)->BufferEmpty();

   Int_t bin, binx, biny, binz;
   Int_t xfirst  = fXaxis.GetFirst();
   Int_t xlast   = fXaxis.GetLast();
   Int_t yfirst  = fYaxis.GetFirst();
   Int_t ylast   = fYaxis.GetLast();
   Int_t zfirst  = fZaxis.GetFirst();
   Int_t zlast   = fZaxis.GetLast();
   Double_t maximum = -FLT_MAX, value;
   for (binz=zfirst;binz<=zlast;binz++) {
      for (biny=yfirst;biny<=ylast;biny++) {
         for (binx=xfirst;binx<=xlast;binx++) {
            bin = GetBin(binx,biny,binz);
            value = RetrieveBinContent(bin);
            if (value > maximum && value < maxval) maximum = value;
         }
      }
   }
   return maximum;
}

////////////////////////////////////////////////////////////////////////////////
/// Return location of bin with maximum value in the range.

Int_t TH1::GetMaximumBin() const
{
   Int_t locmax, locmay, locmaz;
   return GetMaximumBin(locmax, locmay, locmaz);
}

////////////////////////////////////////////////////////////////////////////////
/// Return location of bin with maximum value in the range.

Int_t TH1::GetMaximumBin(Int_t &locmax, Int_t &locmay, Int_t &locmaz) const
{
      // empty the buffer
   if (fBuffer) ((TH1*)this)->BufferEmpty();

   Int_t bin, binx, biny, binz;
   Int_t locm;
   Int_t xfirst  = fXaxis.GetFirst();
   Int_t xlast   = fXaxis.GetLast();
   Int_t yfirst  = fYaxis.GetFirst();
   Int_t ylast   = fYaxis.GetLast();
   Int_t zfirst  = fZaxis.GetFirst();
   Int_t zlast   = fZaxis.GetLast();
   Double_t maximum = -FLT_MAX, value;
   locm = locmax = locmay = locmaz = 0;
   for (binz=zfirst;binz<=zlast;binz++) {
      for (biny=yfirst;biny<=ylast;biny++) {
         for (binx=xfirst;binx<=xlast;binx++) {
            bin = GetBin(binx,biny,binz);
            value = RetrieveBinContent(bin);
            if (value > maximum) {
               maximum = value;
               locm    = bin;
               locmax  = binx;
               locmay  = biny;
               locmaz  = binz;
            }
         }
      }
   }
   return locm;
}

////////////////////////////////////////////////////////////////////////////////
/// Return minimum value larger than minval of bins in the range,
/// unless the value has been overridden by TH1::SetMinimum,
/// in which case it returns that value. (This happens, for example,
/// when the histogram is drawn and the y or z axis limits are changed
///
/// To get the minimum value of bins in the histogram regardless of
/// whether the value has been overridden, use
///
/// ~~~ {.cpp}
/// h->GetBinContent(h->GetMinimumBin())
/// ~~~

Double_t TH1::GetMinimum(Double_t minval) const
{
   if (fMinimum != -1111) return fMinimum;

   // empty the buffer
   if (fBuffer) ((TH1*)this)->BufferEmpty();

   Int_t bin, binx, biny, binz;
   Int_t xfirst  = fXaxis.GetFirst();
   Int_t xlast   = fXaxis.GetLast();
   Int_t yfirst  = fYaxis.GetFirst();
   Int_t ylast   = fYaxis.GetLast();
   Int_t zfirst  = fZaxis.GetFirst();
   Int_t zlast   = fZaxis.GetLast();
   Double_t minimum=FLT_MAX, value;
   for (binz=zfirst;binz<=zlast;binz++) {
      for (biny=yfirst;biny<=ylast;biny++) {
         for (binx=xfirst;binx<=xlast;binx++) {
            bin = GetBin(binx,biny,binz);
            value = RetrieveBinContent(bin);
            if (value < minimum && value > minval) minimum = value;
         }
      }
   }
   return minimum;
}

////////////////////////////////////////////////////////////////////////////////
/// Return location of bin with minimum value in the range.

Int_t TH1::GetMinimumBin() const
{
   Int_t locmix, locmiy, locmiz;
   return GetMinimumBin(locmix, locmiy, locmiz);
}

////////////////////////////////////////////////////////////////////////////////
/// Return location of bin with minimum value in the range.

Int_t TH1::GetMinimumBin(Int_t &locmix, Int_t &locmiy, Int_t &locmiz) const
{
      // empty the buffer
   if (fBuffer) ((TH1*)this)->BufferEmpty();

   Int_t bin, binx, biny, binz;
   Int_t locm;
   Int_t xfirst  = fXaxis.GetFirst();
   Int_t xlast   = fXaxis.GetLast();
   Int_t yfirst  = fYaxis.GetFirst();
   Int_t ylast   = fYaxis.GetLast();
   Int_t zfirst  = fZaxis.GetFirst();
   Int_t zlast   = fZaxis.GetLast();
   Double_t minimum = FLT_MAX, value;
   locm = locmix = locmiy = locmiz = 0;
   for (binz=zfirst;binz<=zlast;binz++) {
      for (biny=yfirst;biny<=ylast;biny++) {
         for (binx=xfirst;binx<=xlast;binx++) {
            bin = GetBin(binx,biny,binz);
            value = RetrieveBinContent(bin);
            if (value < minimum) {
               minimum = value;
               locm    = bin;
               locmix  = binx;
               locmiy  = biny;
               locmiz  = binz;
            }
         }
      }
   }
   return locm;
}

///////////////////////////////////////////////////////////////////////////////
/// Retrieve the minimum and maximum values in the histogram
///
/// This will not return a cached value and will always search the
/// histogram for the min and max values. The user can condition whether
/// or not to call this with the GetMinimumStored() and GetMaximumStored()
/// methods. If the cache is empty, then the value will be -1111. Users
/// can then use the SetMinimum() or SetMaximum() methods to cache the results.
/// For example, the following recipe will make efficient use of this method
/// and the cached minimum and maximum values.
//
/// \code{.cpp}
/// Double_t currentMin = pHist->GetMinimumStored();
/// Double_t currentMax = pHist->GetMaximumStored();
/// if ((currentMin == -1111) || (currentMax == -1111)) {
///    pHist->GetMinimumAndMaximum(currentMin, currentMax);
///    pHist->SetMinimum(currentMin);
///    pHist->SetMaximum(currentMax);
/// }
/// \endcode
///
/// \param min     reference to variable that will hold found minimum value
/// \param max     reference to variable that will hold found maximum value

void TH1::GetMinimumAndMaximum(Double_t& min, Double_t& max) const
{
   // empty the buffer
   if (fBuffer) ((TH1*)this)->BufferEmpty();

   Int_t bin, binx, biny, binz;
   Int_t xfirst  = fXaxis.GetFirst();
   Int_t xlast   = fXaxis.GetLast();
   Int_t yfirst  = fYaxis.GetFirst();
   Int_t ylast   = fYaxis.GetLast();
   Int_t zfirst  = fZaxis.GetFirst();
   Int_t zlast   = fZaxis.GetLast();
   min=TMath::Infinity();
   max=-TMath::Infinity();
   Double_t value;
   for (binz=zfirst;binz<=zlast;binz++) {
      for (biny=yfirst;biny<=ylast;biny++) {
         for (binx=xfirst;binx<=xlast;binx++) {
            bin = GetBin(binx,biny,binz);
            value = RetrieveBinContent(bin);
            if (value < min) min = value;
            if (value > max) max = value;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Redefine  x axis parameters.
///
/// The X axis parameters are modified.
/// The bins content array is resized
/// if errors (Sumw2) the errors array is resized
/// The previous bin contents are lost
/// To change only the axis limits, see TAxis::SetRange

void TH1::SetBins(Int_t nx, Double_t xmin, Double_t xmax)
{
   if (GetDimension() != 1) {
      Error("SetBins","Operation only valid for 1-d histograms");
      return;
   }
   fXaxis.SetRange(0,0);
   fXaxis.Set(nx,xmin,xmax);
   fYaxis.Set(1,0,1);
   fZaxis.Set(1,0,1);
   fNcells = nx+2;
   SetBinsLength(fNcells);
   if (fSumw2.fN) {
      fSumw2.Set(fNcells);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Redefine  x axis parameters with variable bin sizes.
///
/// The X axis parameters are modified.
/// The bins content array is resized
/// if errors (Sumw2) the errors array is resized
/// The previous bin contents are lost
/// To change only the axis limits, see TAxis::SetRange
/// xBins is supposed to be of length nx+1

void TH1::SetBins(Int_t nx, const Double_t *xBins)
{
   if (GetDimension() != 1) {
      Error("SetBins","Operation only valid for 1-d histograms");
      return;
   }
   fXaxis.SetRange(0,0);
   fXaxis.Set(nx,xBins);
   fYaxis.Set(1,0,1);
   fZaxis.Set(1,0,1);
   fNcells = nx+2;
   SetBinsLength(fNcells);
   if (fSumw2.fN) {
      fSumw2.Set(fNcells);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Redefine  x and y axis parameters.
///
/// The X and Y axis parameters are modified.
/// The bins content array is resized
/// if errors (Sumw2) the errors array is resized
/// The previous bin contents are lost
/// To change only the axis limits, see TAxis::SetRange

void TH1::SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax)
{
   if (GetDimension() != 2) {
      Error("SetBins","Operation only valid for 2-D histograms");
      return;
   }
   fXaxis.SetRange(0,0);
   fYaxis.SetRange(0,0);
   fXaxis.Set(nx,xmin,xmax);
   fYaxis.Set(ny,ymin,ymax);
   fZaxis.Set(1,0,1);
   fNcells = (nx+2)*(ny+2);
   SetBinsLength(fNcells);
   if (fSumw2.fN) {
      fSumw2.Set(fNcells);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Redefine  x and y axis parameters with variable bin sizes.
///
/// The X and Y axis parameters are modified.
/// The bins content array is resized
/// if errors (Sumw2) the errors array is resized
/// The previous bin contents are lost
/// To change only the axis limits, see TAxis::SetRange
/// xBins is supposed to be of length nx+1, yBins is supposed to be of length ny+1

void TH1::SetBins(Int_t nx, const Double_t *xBins, Int_t ny, const Double_t *yBins)
{
   if (GetDimension() != 2) {
      Error("SetBins","Operation only valid for 2-D histograms");
      return;
   }
   fXaxis.SetRange(0,0);
   fYaxis.SetRange(0,0);
   fXaxis.Set(nx,xBins);
   fYaxis.Set(ny,yBins);
   fZaxis.Set(1,0,1);
   fNcells = (nx+2)*(ny+2);
   SetBinsLength(fNcells);
   if (fSumw2.fN) {
      fSumw2.Set(fNcells);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Redefine  x, y and z axis parameters.
///
/// The X, Y and Z axis parameters are modified.
/// The bins content array is resized
/// if errors (Sumw2) the errors array is resized
/// The previous bin contents are lost
/// To change only the axis limits, see TAxis::SetRange

void TH1::SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax, Int_t nz, Double_t zmin, Double_t zmax)
{
   if (GetDimension() != 3) {
      Error("SetBins","Operation only valid for 3-D histograms");
      return;
   }
   fXaxis.SetRange(0,0);
   fYaxis.SetRange(0,0);
   fZaxis.SetRange(0,0);
   fXaxis.Set(nx,xmin,xmax);
   fYaxis.Set(ny,ymin,ymax);
   fZaxis.Set(nz,zmin,zmax);
   fNcells = (nx+2)*(ny+2)*(nz+2);
   SetBinsLength(fNcells);
   if (fSumw2.fN) {
      fSumw2.Set(fNcells);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Redefine  x, y and z axis parameters with variable bin sizes.
///
/// The X, Y and Z axis parameters are modified.
/// The bins content array is resized
/// if errors (Sumw2) the errors array is resized
/// The previous bin contents are lost
/// To change only the axis limits, see TAxis::SetRange
/// xBins is supposed to be of length nx+1, yBins is supposed to be of length ny+1,
/// zBins is supposed to be of length nz+1

void TH1::SetBins(Int_t nx, const Double_t *xBins, Int_t ny, const Double_t *yBins, Int_t nz, const Double_t *zBins)
{
   if (GetDimension() != 3) {
      Error("SetBins","Operation only valid for 3-D histograms");
      return;
   }
   fXaxis.SetRange(0,0);
   fYaxis.SetRange(0,0);
   fZaxis.SetRange(0,0);
   fXaxis.Set(nx,xBins);
   fYaxis.Set(ny,yBins);
   fZaxis.Set(nz,zBins);
   fNcells = (nx+2)*(ny+2)*(nz+2);
   SetBinsLength(fNcells);
   if (fSumw2.fN) {
      fSumw2.Set(fNcells);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// By default when an histogram is created, it is added to the list
/// of histogram objects in the current directory in memory.
/// Remove reference to this histogram from current directory and add
/// reference to new directory dir. dir can be 0 in which case the
/// histogram does not belong to any directory.
///
/// Note that the directory is not a real property of the histogram and
/// it will not be copied when the histogram is copied or cloned.
/// If the user wants to have the copied (cloned) histogram in the same
/// directory, he needs to set again the directory using SetDirectory to the
/// copied histograms

void TH1::SetDirectory(TDirectory *dir)
{
   if (fDirectory == dir) return;
   if (fDirectory) fDirectory->Remove(this);
   fDirectory = dir;
   if (fDirectory) {
      fFunctions->UseRWLock();
      fDirectory->Append(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Replace bin errors by values in array error.

void TH1::SetError(const Double_t *error)
{
   for (Int_t i = 0; i < fNcells; ++i) SetBinError(i, error[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Change the name of this histogram
///

void TH1::SetName(const char *name)
{
   //  Histograms are named objects in a THashList.
   //  We must update the hashlist if we change the name
   //  We protect this operation
   R__LOCKGUARD(gROOTMutex);
   if (fDirectory) fDirectory->Remove(this);
   fName = name;
   if (fDirectory) fDirectory->Append(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Change the name and title of this histogram

void TH1::SetNameTitle(const char *name, const char *title)
{
   //  Histograms are named objects in a THashList.
   //  We must update the hashlist if we change the name
   SetName(name);
   SetTitle(title);
}

////////////////////////////////////////////////////////////////////////////////
/// Set statistics option on/off.
///
/// By default, the statistics box is drawn.
/// The paint options can be selected via gStyle->SetOptStats.
/// This function sets/resets the kNoStats bit in the histogram object.
/// It has priority over the Style option.

void TH1::SetStats(Bool_t stats)
{
   ResetBit(kNoStats);
   if (!stats) {
      SetBit(kNoStats);
      //remove the "stats" object from the list of functions
      if (fFunctions) {
         TObject *obj = fFunctions->FindObject("stats");
         if (obj) {
            fFunctions->Remove(obj);
            delete obj;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create structure to store sum of squares of weights.
///
/// if histogram is already filled, the sum of squares of weights
/// is filled with the existing bin contents
///
/// The error per bin will be computed as sqrt(sum of squares of weight)
/// for each bin.
///
/// This function is automatically called when the histogram is created
/// if the static function TH1::SetDefaultSumw2 has been called before.
/// If flag = false the structure containing the sum of the square of weights
/// is rest and it will be empty, but it is not deleted (i.e. GetSumw2()->fN = 0)

void TH1::Sumw2(Bool_t flag)
{
   if (!flag) {
      // clear the array if existing - do nothing otherwise
      if (fSumw2.fN > 0 ) fSumw2.Set(0);
      return;
   }

   if (fSumw2.fN == fNcells) {
      if (!fgDefaultSumw2 )
         Warning("Sumw2","Sum of squares of weights structure already created");
      return;
   }

   fSumw2.Set(fNcells);

   // empty the buffer
   if (fBuffer) BufferEmpty();

   if (fEntries > 0)
      for (Int_t i = 0; i < fNcells; ++i)
         fSumw2.fArray[i] = TMath::Abs(RetrieveBinContent(i));
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to function with name.
///
///
/// Functions such as TH1::Fit store the fitted function in the list of
/// functions of this histogram.

TF1 *TH1::GetFunction(const char *name) const
{
   return (TF1*)fFunctions->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Return value of error associated to bin number bin.
///
/// if the sum of squares of weights has been defined (via Sumw2),
/// this function returns the sqrt(sum of w2).
/// otherwise it returns the sqrt(contents) for this bin.

Double_t TH1::GetBinError(Int_t bin) const
{
   if (bin < 0) bin = 0;
   if (bin >= fNcells) bin = fNcells-1;
   if (fBuffer) ((TH1*)this)->BufferEmpty();
   if (fSumw2.fN) return TMath::Sqrt(fSumw2.fArray[bin]);

   return TMath::Sqrt(TMath::Abs(RetrieveBinContent(bin)));
}

////////////////////////////////////////////////////////////////////////////////
/// Return lower error associated to bin number bin.
///
/// The error will depend on the statistic option used will return
/// the binContent - lower interval value

Double_t TH1::GetBinErrorLow(Int_t bin) const
{
   if (fBinStatErrOpt == kNormal) return GetBinError(bin);
   // in case of weighted histogram check if it is really weighted
   if (fSumw2.fN && fTsumw != fTsumw2) return GetBinError(bin);

   if (bin < 0) bin = 0;
   if (bin >= fNcells) bin = fNcells-1;
   if (fBuffer) ((TH1*)this)->BufferEmpty();

   Double_t alpha = 1.- 0.682689492;
   if (fBinStatErrOpt == kPoisson2) alpha = 0.05;

   Double_t c = RetrieveBinContent(bin);
   Int_t n = int(c);
   if (n < 0) {
      Warning("GetBinErrorLow","Histogram has negative bin content-force usage to normal errors");
      ((TH1*)this)->fBinStatErrOpt = kNormal;
      return GetBinError(bin);
   }

   if (n == 0) return 0;
   return c - ROOT::Math::gamma_quantile( alpha/2, n, 1.);
}

////////////////////////////////////////////////////////////////////////////////
/// Return upper error associated to bin number bin.
///
/// The error will depend on the statistic option used will return
/// the binContent - upper interval value

Double_t TH1::GetBinErrorUp(Int_t bin) const
{
   if (fBinStatErrOpt == kNormal) return GetBinError(bin);
   // in case of weighted histogram check if it is really weighted
   if (fSumw2.fN && fTsumw != fTsumw2) return GetBinError(bin);
   if (bin < 0) bin = 0;
   if (bin >= fNcells) bin = fNcells-1;
   if (fBuffer) ((TH1*)this)->BufferEmpty();

   Double_t alpha = 1.- 0.682689492;
   if (fBinStatErrOpt == kPoisson2) alpha = 0.05;

   Double_t c = RetrieveBinContent(bin);
   Int_t n = int(c);
   if (n < 0) {
      Warning("GetBinErrorUp","Histogram has negative bin content-force usage to normal errors");
      ((TH1*)this)->fBinStatErrOpt = kNormal;
      return GetBinError(bin);
   }

   // for N==0 return an upper limit at 0.68 or (1-alpha)/2 ?
   // decide to return always (1-alpha)/2 upper interval
   //if (n == 0) return ROOT::Math::gamma_quantile_c(alpha,n+1,1);
   return ROOT::Math::gamma_quantile_c( alpha/2, n+1, 1) - c;
}

//L.M. These following getters are useless and should be probably deprecated
////////////////////////////////////////////////////////////////////////////////
/// Return bin center for 1D histogram.
/// Better to use h1.GetXaxis()->GetBinCenter(bin)

Double_t TH1::GetBinCenter(Int_t bin) const
{
   if (fDimension == 1) return  fXaxis.GetBinCenter(bin);
   Error("GetBinCenter","Invalid method for a %d-d histogram - return a NaN",fDimension);
   return TMath::QuietNaN();
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin lower edge for 1D histogram.
/// Better to use h1.GetXaxis()->GetBinLowEdge(bin)

Double_t TH1::GetBinLowEdge(Int_t bin) const
{
   if (fDimension == 1) return  fXaxis.GetBinLowEdge(bin);
   Error("GetBinLowEdge","Invalid method for a %d-d histogram - return a NaN",fDimension);
   return TMath::QuietNaN();
}

////////////////////////////////////////////////////////////////////////////////
/// Return bin width for 1D histogram.
/// Better to use h1.GetXaxis()->GetBinWidth(bin)

Double_t TH1::GetBinWidth(Int_t bin) const
{
   if (fDimension == 1) return  fXaxis.GetBinWidth(bin);
   Error("GetBinWidth","Invalid method for a %d-d histogram - return a NaN",fDimension);
   return TMath::QuietNaN();
}

////////////////////////////////////////////////////////////////////////////////
/// Fill array with center of bins for 1D histogram
/// Better to use h1.GetXaxis()->GetCenter(center)

void TH1::GetCenter(Double_t *center) const
{
   if (fDimension == 1) {
      fXaxis.GetCenter(center);
      return;
   }
   Error("GetCenter","Invalid method for a %d-d histogram ",fDimension);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill array with low edge of bins for 1D histogram
/// Better to use h1.GetXaxis()->GetLowEdge(edge)

void TH1::GetLowEdge(Double_t *edge) const
{
   if (fDimension == 1) {
      fXaxis.GetLowEdge(edge);
      return;
   }
   Error("GetLowEdge","Invalid method for a %d-d histogram ",fDimension);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the bin Error
/// Note that this resets the bin eror option to be of Normal Type and for the
/// non-empty bin the bin error is set by default to the square root of their content.
/// Note that in case the user sets after calling SetBinError explicitly a new bin content (e.g. using SetBinContent)
/// he needs then to provide also the corresponding bin error (using SetBinError) since the bin error
/// will not recalcualated after setting the content and a default error = 0 will be used for those bins.
///
/// See convention for numbering bins in TH1::GetBin

void TH1::SetBinError(Int_t bin, Double_t error)
{
   if (bin < 0 || bin>= fNcells) return;
   if (!fSumw2.fN) Sumw2();
   fSumw2.fArray[bin] = error * error;
   // reset the bin error option
   SetBinErrorOption(kNormal);
}

////////////////////////////////////////////////////////////////////////////////
/// Set bin content
/// see convention for numbering bins in TH1::GetBin
/// In case the bin number is greater than the number of bins and
/// the timedisplay option is set or CanExtendAllAxes(),
/// the number of bins is automatically doubled to accommodate the new bin

void TH1::SetBinContent(Int_t bin, Double_t content)
{
   fEntries++;
   fTsumw = 0;
   if (bin < 0) return;
   if (bin >= fNcells-1) {
      if (fXaxis.GetTimeDisplay() || CanExtendAllAxes() ) {
         while (bin >=  fNcells-1)  LabelsInflate();
      } else {
         if (bin == fNcells-1) UpdateBinContent(bin, content);
         return;
      }
   }
   UpdateBinContent(bin, content);
}

////////////////////////////////////////////////////////////////////////////////
/// See convention for numbering bins in TH1::GetBin

void TH1::SetBinError(Int_t binx, Int_t biny, Double_t error)
{
   if (binx < 0 || binx > fXaxis.GetNbins() + 1) return;
   if (biny < 0 || biny > fYaxis.GetNbins() + 1) return;
   SetBinError(GetBin(binx, biny), error);
}

////////////////////////////////////////////////////////////////////////////////
/// See convention for numbering bins in TH1::GetBin

void TH1::SetBinError(Int_t binx, Int_t biny, Int_t binz, Double_t error)
{
   if (binx < 0 || binx > fXaxis.GetNbins() + 1) return;
   if (biny < 0 || biny > fYaxis.GetNbins() + 1) return;
   if (binz < 0 || binz > fZaxis.GetNbins() + 1) return;
   SetBinError(GetBin(binx, biny, binz), error);
}

////////////////////////////////////////////////////////////////////////////////
/// This function calculates the background spectrum in this histogram.
/// The background is returned as a histogram.
///
/// \param[in] niter number of iterations (default value = 2)
///     Increasing niter make the result smoother and lower.
/// \param[in] option may contain one of the following options
///  - to set the direction parameter
///    "BackDecreasingWindow". By default the direction is BackIncreasingWindow
///  - filterOrder-order of clipping filter (default "BackOrder2")
///    possible values= "BackOrder4" "BackOrder6" "BackOrder8"
///  - "nosmoothing" - if selected, the background is not smoothed
///    By default the background is smoothed.
///  - smoothWindow - width of smoothing window, (default is "BackSmoothing3")
///    possible values= "BackSmoothing5" "BackSmoothing7" "BackSmoothing9"
///    "BackSmoothing11" "BackSmoothing13" "BackSmoothing15"
///  - "nocompton" - if selected the estimation of Compton edge
///    will be not be included   (by default the compton estimation is set)
///  - "same" if this option is specified, the resulting background
///    histogram is superimposed on the picture in the current pad.
///    This option is given by default.
///
/// NOTE that the background is only evaluated in the current range of this histogram.
/// i.e., if this has a bin range (set via h->GetXaxis()->SetRange(binmin, binmax),
/// the returned histogram will be created with the same number of bins
/// as this input histogram, but only bins from binmin to binmax will be filled
/// with the estimated background.

TH1 *TH1::ShowBackground(Int_t niter, Option_t *option)
{

   return (TH1*)gROOT->ProcessLineFast(Form("TSpectrum::StaticBackground((TH1*)0x%lx,%d,\"%s\")",
                                            (ULong_t)this, niter, option));
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to TSpectrum::Search.
/// The function finds peaks in this histogram where the width is > sigma
/// and the peak maximum greater than threshold*maximum bin content of this.
/// For more details see TSpectrum::Search.
/// Note the difference in the default value for option compared to TSpectrum::Search
/// option="" by default (instead of "goff").

Int_t TH1::ShowPeaks(Double_t sigma, Option_t *option, Double_t threshold)
{
   return (Int_t)gROOT->ProcessLineFast(Form("TSpectrum::StaticSearch((TH1*)0x%lx,%g,\"%s\",%g)",
                                             (ULong_t)this, sigma, option, threshold));
}

////////////////////////////////////////////////////////////////////////////////
/// For a given transform (first parameter), fills the histogram (second parameter)
/// with the transform output data, specified in the third parameter
/// If the 2nd parameter h_output is empty, a new histogram (TH1D or TH2D) is created
/// and the user is responsible for deleting it.
///
/// Available options:
///  - "RE" - real part of the output
///  - "IM" - imaginary part of the output
///  - "MAG" - magnitude of the output
///  - "PH"  - phase of the output

TH1* TH1::TransformHisto(TVirtualFFT *fft, TH1* h_output,  Option_t *option)
{
   if (!fft ||  !fft->GetN() ) {
      ::Error("TransformHisto","Invalid FFT transform class");
      return 0;
   }

   if (fft->GetNdim()>2){
      ::Error("TransformHisto","Only 1d and 2D transform are supported");
      return 0;
   }
   Int_t binx,biny;
   TString opt = option;
   opt.ToUpper();
   Int_t *n = fft->GetN();
   TH1 *hout=0;
   if (h_output) {
      hout = h_output;
   }
   else {
      TString name = TString::Format("out_%s", opt.Data());
      if (fft->GetNdim()==1)
         hout = new TH1D(name, name,n[0], 0, n[0]);
      else if (fft->GetNdim()==2)
         hout = new TH2D(name, name, n[0], 0, n[0], n[1], 0, n[1]);
   }
   R__ASSERT(hout != 0);
   TString type=fft->GetType();
   Int_t ind[2];
   if (opt.Contains("RE")){
      if (type.Contains("2C") || type.Contains("2HC")) {
         Double_t re, im;
         for (binx = 1; binx<=hout->GetNbinsX(); binx++) {
            for (biny=1; biny<=hout->GetNbinsY(); biny++) {
               ind[0] = binx-1; ind[1] = biny-1;
               fft->GetPointComplex(ind, re, im);
               hout->SetBinContent(binx, biny, re);
            }
         }
      } else {
         for (binx = 1; binx<=hout->GetNbinsX(); binx++) {
            for (biny=1; biny<=hout->GetNbinsY(); biny++) {
               ind[0] = binx-1; ind[1] = biny-1;
               hout->SetBinContent(binx, biny, fft->GetPointReal(ind));
            }
         }
      }
   }
   if (opt.Contains("IM")) {
      if (type.Contains("2C") || type.Contains("2HC")) {
         Double_t re, im;
         for (binx = 1; binx<=hout->GetNbinsX(); binx++) {
            for (biny=1; biny<=hout->GetNbinsY(); biny++) {
               ind[0] = binx-1; ind[1] = biny-1;
               fft->GetPointComplex(ind, re, im);
               hout->SetBinContent(binx, biny, im);
            }
         }
      } else {
         ::Error("TransformHisto","No complex numbers in the output");
         return 0;
      }
   }
   if (opt.Contains("MA")) {
      if (type.Contains("2C") || type.Contains("2HC")) {
         Double_t re, im;
         for (binx = 1; binx<=hout->GetNbinsX(); binx++) {
            for (biny=1; biny<=hout->GetNbinsY(); biny++) {
               ind[0] = binx-1; ind[1] = biny-1;
               fft->GetPointComplex(ind, re, im);
               hout->SetBinContent(binx, biny, TMath::Sqrt(re*re + im*im));
            }
         }
      } else {
         for (binx = 1; binx<=hout->GetNbinsX(); binx++) {
            for (biny=1; biny<=hout->GetNbinsY(); biny++) {
               ind[0] = binx-1; ind[1] = biny-1;
               hout->SetBinContent(binx, biny, TMath::Abs(fft->GetPointReal(ind)));
            }
         }
      }
   }
   if (opt.Contains("PH")) {
      if (type.Contains("2C") || type.Contains("2HC")){
         Double_t re, im, ph;
         for (binx = 1; binx<=hout->GetNbinsX(); binx++){
            for (biny=1; biny<=hout->GetNbinsY(); biny++){
               ind[0] = binx-1; ind[1] = biny-1;
               fft->GetPointComplex(ind, re, im);
               if (TMath::Abs(re) > 1e-13){
                  ph = TMath::ATan(im/re);
                  //find the correct quadrant
                  if (re<0 && im<0)
                     ph -= TMath::Pi();
                  if (re<0 && im>=0)
                     ph += TMath::Pi();
               } else {
                  if (TMath::Abs(im) < 1e-13)
                     ph = 0;
                  else if (im>0)
                     ph = TMath::Pi()*0.5;
                  else
                     ph = -TMath::Pi()*0.5;
               }
               hout->SetBinContent(binx, biny, ph);
            }
         }
      } else {
         printf("Pure real output, no phase");
         return 0;
      }
   }

   return hout;
}

////////////////////////////////////////////////////////////////////////////////
/// Raw retrieval of bin content on internal data structure
/// see convention for numbering bins in TH1::GetBin

Double_t TH1::RetrieveBinContent(Int_t) const
{
   AbstractMethod("RetrieveBinContent");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Raw update of bin content on internal data structure
/// see convention for numbering bins in TH1::GetBin

void TH1::UpdateBinContent(Int_t, Double_t)
{
   AbstractMethod("UpdateBinContent");
}

////////////////////////////////////////////////////////////////////////////////
/// Print value overload

std::string cling::printValue(TH1 *val) {
  std::ostringstream strm;
  strm << cling::printValue((TObject*)val) << " NbinsX: " << val->GetNbinsX();
  return strm.str();
}

//______________________________________________________________________________
//                     TH1C methods
// TH1C : histograms with one byte per channel.   Maximum bin content = 127
//______________________________________________________________________________

ClassImp(TH1C);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH1C::TH1C(): TH1(), TArrayC()
{
   fDimension = 1;
   SetBinsLength(3);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with fix bins of type char (one byte per channel)
/// (see TH1::TH1 for explanation of parameters)

TH1C::TH1C(const char *name,const char *title,Int_t nbins,Double_t xlow,Double_t xup)
: TH1(name,title,nbins,xlow,xup)
{
   fDimension = 1;
   TArrayC::Set(fNcells);

   if (xlow >= xup) SetBuffer(fgBufferSize);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with variable bins of type char (one byte per channel)
/// (see TH1::TH1 for explanation of parameters)

TH1C::TH1C(const char *name,const char *title,Int_t nbins,const Float_t *xbins)
: TH1(name,title,nbins,xbins)
{
   fDimension = 1;
   TArrayC::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with variable bins of type char (one byte per channel)
/// (see TH1::TH1 for explanation of parameters)

TH1C::TH1C(const char *name,const char *title,Int_t nbins,const Double_t *xbins)
: TH1(name,title,nbins,xbins)
{
   fDimension = 1;
   TArrayC::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH1C::~TH1C()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TH1C::TH1C(const TH1C &h1c) : TH1(), TArrayC()
{
   ((TH1C&)h1c).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by 1.

void TH1C::AddBinContent(Int_t bin)
{
   if (fArray[bin] < 127) fArray[bin]++;
}

////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by w.

void TH1C::AddBinContent(Int_t bin, Double_t w)
{
   Int_t newval = fArray[bin] + Int_t(w);
   if (newval > -128 && newval < 128) {fArray[bin] = Char_t(newval); return;}
   if (newval < -127) fArray[bin] = -127;
   if (newval >  127) fArray[bin] =  127;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this to newth1

void TH1C::Copy(TObject &newth1) const
{
   TH1::Copy(newth1);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset.

void TH1C::Reset(Option_t *option)
{
   TH1::Reset(option);
   TArrayC::Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH1C::SetBinsLength(Int_t n)
{
   if (n < 0) n = fXaxis.GetNbins() + 2;
   fNcells = n;
   TArrayC::Set(n);
}

////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH1C& TH1C::operator=(const TH1C &h1)
{
   if (this != &h1)  ((TH1C&)h1).Copy(*this);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH1C operator*(Double_t c1, const TH1C &h1)
{
   TH1C hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH1C operator+(const TH1C &h1, const TH1C &h2)
{
   TH1C hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH1C operator-(const TH1C &h1, const TH1C &h2)
{
   TH1C hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH1C operator*(const TH1C &h1, const TH1C &h2)
{
   TH1C hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH1C operator/(const TH1C &h1, const TH1C &h2)
{
   TH1C hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}

//______________________________________________________________________________
//                     TH1S methods
// TH1S : histograms with one short per channel.  Maximum bin content = 32767
//______________________________________________________________________________

ClassImp(TH1S);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH1S::TH1S(): TH1(), TArrayS()
{
   fDimension = 1;
   SetBinsLength(3);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with fix bins of type short
/// (see TH1::TH1 for explanation of parameters)

TH1S::TH1S(const char *name,const char *title,Int_t nbins,Double_t xlow,Double_t xup)
: TH1(name,title,nbins,xlow,xup)
{
   fDimension = 1;
   TArrayS::Set(fNcells);

   if (xlow >= xup) SetBuffer(fgBufferSize);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with variable bins of type short
/// (see TH1::TH1 for explanation of parameters)

TH1S::TH1S(const char *name,const char *title,Int_t nbins,const Float_t *xbins)
: TH1(name,title,nbins,xbins)
{
   fDimension = 1;
   TArrayS::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with variable bins of type short
/// (see TH1::TH1 for explanation of parameters)

TH1S::TH1S(const char *name,const char *title,Int_t nbins,const Double_t *xbins)
: TH1(name,title,nbins,xbins)
{
   fDimension = 1;
   TArrayS::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH1S::~TH1S()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TH1S::TH1S(const TH1S &h1s) : TH1(), TArrayS()
{
   ((TH1S&)h1s).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by 1.

void TH1S::AddBinContent(Int_t bin)
{
   if (fArray[bin] < 32767) fArray[bin]++;
}

////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by w

void TH1S::AddBinContent(Int_t bin, Double_t w)
{
   Int_t newval = fArray[bin] + Int_t(w);
   if (newval > -32768 && newval < 32768) {fArray[bin] = Short_t(newval); return;}
   if (newval < -32767) fArray[bin] = -32767;
   if (newval >  32767) fArray[bin] =  32767;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this to newth1

void TH1S::Copy(TObject &newth1) const
{
   TH1::Copy(newth1);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset.

void TH1S::Reset(Option_t *option)
{
   TH1::Reset(option);
   TArrayS::Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH1S::SetBinsLength(Int_t n)
{
   if (n < 0) n = fXaxis.GetNbins() + 2;
   fNcells = n;
   TArrayS::Set(n);
}

////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH1S& TH1S::operator=(const TH1S &h1)
{
   if (this != &h1)  ((TH1S&)h1).Copy(*this);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH1S operator*(Double_t c1, const TH1S &h1)
{
   TH1S hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH1S operator+(const TH1S &h1, const TH1S &h2)
{
   TH1S hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH1S operator-(const TH1S &h1, const TH1S &h2)
{
   TH1S hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH1S operator*(const TH1S &h1, const TH1S &h2)
{
   TH1S hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH1S operator/(const TH1S &h1, const TH1S &h2)
{
   TH1S hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}

//______________________________________________________________________________
//                     TH1I methods
// TH1I : histograms with one int per channel.    Maximum bin content = 2147483647
//______________________________________________________________________________

ClassImp(TH1I);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH1I::TH1I(): TH1(), TArrayI()
{
   fDimension = 1;
   SetBinsLength(3);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with fix bins of type integer
/// (see TH1::TH1 for explanation of parameters)

TH1I::TH1I(const char *name,const char *title,Int_t nbins,Double_t xlow,Double_t xup)
: TH1(name,title,nbins,xlow,xup)
{
   fDimension = 1;
   TArrayI::Set(fNcells);

   if (xlow >= xup) SetBuffer(fgBufferSize);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with variable bins of type integer
/// (see TH1::TH1 for explanation of parameters)

TH1I::TH1I(const char *name,const char *title,Int_t nbins,const Float_t *xbins)
: TH1(name,title,nbins,xbins)
{
   fDimension = 1;
   TArrayI::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with variable bins of type integer
/// (see TH1::TH1 for explanation of parameters)

TH1I::TH1I(const char *name,const char *title,Int_t nbins,const Double_t *xbins)
: TH1(name,title,nbins,xbins)
{
   fDimension = 1;
   TArrayI::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH1I::~TH1I()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TH1I::TH1I(const TH1I &h1i) : TH1(), TArrayI()
{
   ((TH1I&)h1i).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by 1.

void TH1I::AddBinContent(Int_t bin)
{
   if (fArray[bin] < 2147483647) fArray[bin]++;
}

////////////////////////////////////////////////////////////////////////////////
/// Increment bin content by w

void TH1I::AddBinContent(Int_t bin, Double_t w)
{
   Long64_t newval = fArray[bin] + Long64_t(w);
   if (newval > -2147483647 && newval < 2147483647) {fArray[bin] = Int_t(newval); return;}
   if (newval < -2147483647) fArray[bin] = -2147483647;
   if (newval >  2147483647) fArray[bin] =  2147483647;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this to newth1

void TH1I::Copy(TObject &newth1) const
{
   TH1::Copy(newth1);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset.

void TH1I::Reset(Option_t *option)
{
   TH1::Reset(option);
   TArrayI::Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH1I::SetBinsLength(Int_t n)
{
   if (n < 0) n = fXaxis.GetNbins() + 2;
   fNcells = n;
   TArrayI::Set(n);
}

////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH1I& TH1I::operator=(const TH1I &h1)
{
   if (this != &h1)  ((TH1I&)h1).Copy(*this);
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH1I operator*(Double_t c1, const TH1I &h1)
{
   TH1I hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH1I operator+(const TH1I &h1, const TH1I &h2)
{
   TH1I hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH1I operator-(const TH1I &h1, const TH1I &h2)
{
   TH1I hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH1I operator*(const TH1I &h1, const TH1I &h2)
{
   TH1I hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH1I operator/(const TH1I &h1, const TH1I &h2)
{
   TH1I hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}

//______________________________________________________________________________
//                     TH1F methods
// TH1F : histograms with one float per channel.  Maximum precision 7 digits
//______________________________________________________________________________

ClassImp(TH1F);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH1F::TH1F(): TH1(), TArrayF()
{
   fDimension = 1;
   SetBinsLength(3);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with fix bins of type float
/// (see TH1::TH1 for explanation of parameters)

TH1F::TH1F(const char *name,const char *title,Int_t nbins,Double_t xlow,Double_t xup)
: TH1(name,title,nbins,xlow,xup)
{
   fDimension = 1;
   TArrayF::Set(fNcells);

   if (xlow >= xup) SetBuffer(fgBufferSize);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with variable bins of type float
/// (see TH1::TH1 for explanation of parameters)

TH1F::TH1F(const char *name,const char *title,Int_t nbins,const Float_t *xbins)
: TH1(name,title,nbins,xbins)
{
   fDimension = 1;
   TArrayF::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with variable bins of type float
/// (see TH1::TH1 for explanation of parameters)

TH1F::TH1F(const char *name,const char *title,Int_t nbins,const Double_t *xbins)
: TH1(name,title,nbins,xbins)
{
   fDimension = 1;
   TArrayF::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a histogram from a TVectorF
/// by default the histogram name is "TVectorF" and title = ""

TH1F::TH1F(const TVectorF &v)
: TH1("TVectorF","",v.GetNrows(),0,v.GetNrows())
{
   TArrayF::Set(fNcells);
   fDimension = 1;
   Int_t ivlow  = v.GetLwb();
   for (Int_t i=0;i<fNcells-2;i++) {
      SetBinContent(i+1,v(i+ivlow));
   }
   TArrayF::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy Constructor.

TH1F::TH1F(const TH1F &h) : TH1(), TArrayF()
{
   ((TH1F&)h).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH1F::~TH1F()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this to newth1.

void TH1F::Copy(TObject &newth1) const
{
   TH1::Copy(newth1);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset.

void TH1F::Reset(Option_t *option)
{
   TH1::Reset(option);
   TArrayF::Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH1F::SetBinsLength(Int_t n)
{
   if (n < 0) n = fXaxis.GetNbins() + 2;
   fNcells = n;
   TArrayF::Set(n);
}

////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH1F& TH1F::operator=(const TH1F &h1)
{
   if (this != &h1)  ((TH1F&)h1).Copy(*this);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH1F operator*(Double_t c1, const TH1F &h1)
{
   TH1F hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH1F operator+(const TH1F &h1, const TH1F &h2)
{
   TH1F hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH1F operator-(const TH1F &h1, const TH1F &h2)
{
   TH1F hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH1F operator*(const TH1F &h1, const TH1F &h2)
{
   TH1F hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH1F operator/(const TH1F &h1, const TH1F &h2)
{
   TH1F hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}

//______________________________________________________________________________
//                     TH1D methods
// TH1D : histograms with one double per channel. Maximum precision 14 digits
//______________________________________________________________________________

ClassImp(TH1D);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH1D::TH1D(): TH1(), TArrayD()
{
   fDimension = 1;
   SetBinsLength(3);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with fix bins of type double
/// (see TH1::TH1 for explanation of parameters)

TH1D::TH1D(const char *name,const char *title,Int_t nbins,Double_t xlow,Double_t xup)
: TH1(name,title,nbins,xlow,xup)
{
   fDimension = 1;
   TArrayD::Set(fNcells);

   if (xlow >= xup) SetBuffer(fgBufferSize);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with variable bins of type double
/// (see TH1::TH1 for explanation of parameters)

TH1D::TH1D(const char *name,const char *title,Int_t nbins,const Float_t *xbins)
: TH1(name,title,nbins,xbins)
{
   fDimension = 1;
   TArrayD::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a 1-Dim histogram with variable bins of type double
/// (see TH1::TH1 for explanation of parameters)

TH1D::TH1D(const char *name,const char *title,Int_t nbins,const Double_t *xbins)
: TH1(name,title,nbins,xbins)
{
   fDimension = 1;
   TArrayD::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a histogram from a TVectorD
/// by default the histogram name is "TVectorD" and title = ""

TH1D::TH1D(const TVectorD &v)
: TH1("TVectorD","",v.GetNrows(),0,v.GetNrows())
{
   TArrayD::Set(fNcells);
   fDimension = 1;
   Int_t ivlow  = v.GetLwb();
   for (Int_t i=0;i<fNcells-2;i++) {
      SetBinContent(i+1,v(i+ivlow));
   }
   TArrayD::Set(fNcells);
   if (fgDefaultSumw2) Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH1D::~TH1D()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH1D::TH1D(const TH1D &h1d) : TH1(), TArrayD()
{
   ((TH1D&)h1d).Copy(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this to newth1

void TH1D::Copy(TObject &newth1) const
{
   TH1::Copy(newth1);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset.

void TH1D::Reset(Option_t *option)
{
   TH1::Reset(option);
   TArrayD::Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Set total number of bins including under/overflow
/// Reallocate bin contents array

void TH1D::SetBinsLength(Int_t n)
{
   if (n < 0) n = fXaxis.GetNbins() + 2;
   fNcells = n;
   TArrayD::Set(n);
}

////////////////////////////////////////////////////////////////////////////////
/// Operator =

TH1D& TH1D::operator=(const TH1D &h1)
{
   if (this != &h1)  ((TH1D&)h1).Copy(*this);
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH1D operator*(Double_t c1, const TH1D &h1)
{
   TH1D hnew = h1;
   hnew.Scale(c1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator +

TH1D operator+(const TH1D &h1, const TH1D &h2)
{
   TH1D hnew = h1;
   hnew.Add(&h2,1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator -

TH1D operator-(const TH1D &h1, const TH1D &h2)
{
   TH1D hnew = h1;
   hnew.Add(&h2,-1);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator *

TH1D operator*(const TH1D &h1, const TH1D &h2)
{
   TH1D hnew = h1;
   hnew.Multiply(&h2);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator /

TH1D operator/(const TH1D &h1, const TH1D &h2)
{
   TH1D hnew = h1;
   hnew.Divide(&h2);
   hnew.SetDirectory(0);
   return hnew;
}

////////////////////////////////////////////////////////////////////////////////
///return pointer to histogram with name
///hid if id >=0
///h_id if id <0

TH1 *R__H(Int_t hid)
{
   TString hname;
   if(hid >= 0) hname.Form("h%d",hid);
   else         hname.Form("h_%d",hid);
   return (TH1*)gDirectory->Get(hname);
}

////////////////////////////////////////////////////////////////////////////////
///return pointer to histogram with name hname

TH1 *R__H(const char * hname)
{
   return (TH1*)gDirectory->Get(hname);
}


/// \fn void TH1::SetBarOffset(Float_t offset)
/// Set the bar offset as fraction of the bin width for drawing mode "B".
/// This shifts bars to the right on the x axis, and helps to draw bars next to each other.
/// \see THistPainter, SetBarWidth()

/// \fn void TH1::SetBarWidth(Float_t width)
/// Set the width of bars as fraction of the bin width for drawing mode "B".
/// This allows for making bars narrower than the bin width. With SetBarOffset(), this helps to draw multiple bars next to each other.
/// \see THistPainter, SetBarOffset()
