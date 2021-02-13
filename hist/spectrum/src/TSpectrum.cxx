// @(#)root/spectrum:$Id$
// Author: Miroslav Morhac   27/05/99

#include "TSpectrum.h"
#include "TPolyMarker.h"
#include "TVirtualPad.h"
#include "TList.h"
#include "TH1.h"
#include "TMath.h"
#include "snprintf.h"

/** \class TSpectrum
    \ingroup Spectrum
    \brief Advanced Spectra Processing
    \author Miroslav Morhac

 This class contains advanced spectra processing functions for:

 -   One-dimensional background estimation
 -   One-dimensional smoothing
 -   One-dimensional deconvolution
 -   One-dimensional peak search

 The algorithms in this class have been published in the following references:

 1.  M.Morhac et al.: Background elimination methods for multidimensional coincidence gamma-ray spectra. Nuclear Instruments and Methods in Physics Research A 401 (1997) 113-132.
 2.  M.Morhac et al.: Efficient one- and two-dimensional Gold deconvolution and its application to gamma-ray spectra decomposition. Nuclear Instruments and Methods in Physics Research A 401 (1997) 385-408.
 3.  M.Morhac et al.: Identification of peaks in multidimensional coincidence gamma-ray spectra. Nuclear Instruments and Methods in Research Physics A 443(2000), 108-125.

 These NIM papers are also available as doc or ps files from:

 - [Spectrum.doc](https://root.cern.ch/download/Spectrum.doc)
 - [SpectrumDec.ps.gz](https://root.cern.ch/download/SpectrumDec.ps.gz)
 - [SpectrumSrc.ps.gz](https://root.cern.ch/download/SpectrumSrc.ps.gz)
 - [SpectrumBck.ps.gz](https://root.cern.ch/download/SpectrumBck.ps.gz)

 See also the
 [online documentation](https://root.cern.ch/guides/tspectrum-manual) and
 [tutorials](https://root.cern.ch/doc/master/group__tutorial__spectrum.html).
*/

Int_t TSpectrum::fgIterations    = 3;
Int_t TSpectrum::fgAverageWindow = 3;

#define PEAK_WINDOW 1024
ClassImp(TSpectrum);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TSpectrum::TSpectrum() :TNamed("Spectrum", "Miroslav Morhac peak finder")
{
   Int_t n = 100;
   fMaxPeaks  = n;
   fPosition   = new Double_t[n];
   fPositionX  = new Double_t[n];
   fPositionY  = new Double_t[n];
   fResolution = 1;
   fHistogram  = 0;
   fNPeaks     = 0;
}

////////////////////////////////////////////////////////////////////////////////
///   - maxpositions: maximum number of peaks
///   - resolution: *NOT USED* determines resolution of the neighbouring peaks
///                   default value is 1 correspond to 3 sigma distance
///                   between peaks. Higher values allow higher resolution
///                   (smaller distance between peaks.
///                   May be set later through SetResolution.

TSpectrum::TSpectrum(Int_t maxpositions, Double_t resolution)
          :TNamed("Spectrum", "Miroslav Morhac peak finder")
{
   Int_t n = maxpositions;
   if (n <= 0) n = 1;
   fMaxPeaks  = n;
   fPosition  = new Double_t[n];
   fPositionX = new Double_t[n];
   fPositionY = new Double_t[n];
   fHistogram = 0;
   fNPeaks    = 0;
   SetResolution(resolution);
}

////////////////////////////////////////////////////////////////////////////////
///   Destructor.

TSpectrum::~TSpectrum()
{
   delete [] fPosition;
   delete [] fPositionX;
   delete [] fPositionY;
   delete    fHistogram;
}

////////////////////////////////////////////////////////////////////////////////
///   Static function: Set average window of searched peaks
///   (see TSpectrum::SearchHighRes).

void TSpectrum::SetAverageWindow(Int_t w)
{
   fgAverageWindow = w;
}

////////////////////////////////////////////////////////////////////////////////
///   Static function: Set max number of decon iterations in deconvolution
///   operation (see TSpectrum::SearchHighRes).

void TSpectrum::SetDeconIterations(Int_t n)
{
   fgIterations = n;
}

////////////////////////////////////////////////////////////////////////////////
///   One-dimensional background estimation function.
///
///   This function calculates the background spectrum in the input histogram h.
///   The background is returned as a histogram.
///
/// #### Parameters:
///
///   - h: input 1-d histogram
///   - numberIterations, (default value = 20)
///      Increasing numberIterations make the result smoother and lower.
///   - option: may contain one of the following options:
///
///      - to set the direction parameter
///      "BackIncreasingWindow". By default the direction is BackDecreasingWindow
///      - filterOrder-order of clipping filter,  (default "BackOrder2")
///                  -possible values= "BackOrder4"
///                                    "BackOrder6"
///                                    "BackOrder8"
///      - "nosmoothing"- if selected, the background is not smoothed
///           By default the background is smoothed.
///      - smoothWindow-width of smoothing window, (default is "BackSmoothing3")
///                  -possible values= "BackSmoothing5"
///                                    "BackSmoothing7"
///                                    "BackSmoothing9"
///                                    "BackSmoothing11"
///                                    "BackSmoothing13"
///                                    "BackSmoothing15"
///      - "Compton" if selected the estimation of Compton edge
///                  will be included.
///      - "same" : if this option is specified, the resulting background
///                 histogram is superimposed on the picture in the current pad.
///
///   NOTE that the background is only evaluated in the current range of h.
///   ie, if h has a bin range (set via `h->GetXaxis()->SetRange(binmin,binmax)`,
///   the returned histogram will be created with the same number of bins
///   as the input histogram h, but only bins from `binmin` to `binmax` will be filled
///   with the estimated background.

TH1 *TSpectrum::Background(const TH1 * h, int numberIterations,
                           Option_t * option)
{
   if (h == 0) return 0;
   Int_t dimension = h->GetDimension();
   if (dimension > 1) {
      Error("Search", "Only implemented for 1-d histograms");
      return 0;
   }
   TString opt = option;
   opt.ToLower();

   //set options
   Int_t direction = kBackDecreasingWindow;
   if (opt.Contains("backincreasingwindow")) direction = kBackIncreasingWindow;
   Int_t filterOrder = kBackOrder2;
   if (opt.Contains("backorder4")) filterOrder = kBackOrder4;
   if (opt.Contains("backorder6")) filterOrder = kBackOrder6;
   if (opt.Contains("backorder8")) filterOrder = kBackOrder8;
   Bool_t smoothing = kTRUE;
   if (opt.Contains("nosmoothing")) smoothing = kFALSE;
   Int_t smoothWindow = kBackSmoothing3;
   if (opt.Contains("backsmoothing5"))  smoothWindow = kBackSmoothing5;
   if (opt.Contains("backsmoothing7"))  smoothWindow = kBackSmoothing7;
   if (opt.Contains("backsmoothing9"))  smoothWindow = kBackSmoothing9;
   if (opt.Contains("backsmoothing11")) smoothWindow = kBackSmoothing11;
   if (opt.Contains("backsmoothing13")) smoothWindow = kBackSmoothing13;
   if (opt.Contains("backsmoothing15")) smoothWindow = kBackSmoothing15;
   Bool_t compton = kFALSE;
   if (opt.Contains("compton")) compton = kTRUE;

   Int_t first = h->GetXaxis()->GetFirst();
   Int_t last  = h->GetXaxis()->GetLast();
   Int_t size = last-first+1;
   Int_t i;
   Double_t * source = new Double_t[size];
   for (i = 0; i < size; i++) source[i] = h->GetBinContent(i + first);

   //find background (source is input and in output contains the background
   Background(source,size,numberIterations, direction, filterOrder,smoothing,
              smoothWindow,compton);

   //create output histogram containing background
   //only bins in the range of the input histogram are filled
   Int_t nch = strlen(h->GetName());
   char *hbname = new char[nch+20];
   snprintf(hbname,nch+20,"%s_background",h->GetName());
   TH1 *hb = (TH1*)h->Clone(hbname);
   hb->Reset();
   hb->GetListOfFunctions()->Delete();
   hb->SetLineColor(2);
   for (i=0; i< size; i++) hb->SetBinContent(i+first,source[i]);
   hb->SetEntries(size);

   //if option "same is specified, draw the result in the pad
   if (opt.Contains("same")) {
      if (gPad) delete gPad->GetPrimitive(hbname);
      hb->Draw("same");
   }
   delete [] source;
   delete [] hbname;
   return hb;
}

////////////////////////////////////////////////////////////////////////////////
///   Print the array of positions.

void TSpectrum::Print(Option_t *) const
{
   printf("\nNumber of positions = %d\n",fNPeaks);
   for (Int_t i=0;i<fNPeaks;i++) {
      printf(" x[%d] = %g, y[%d] = %g\n",i,fPositionX[i],i,fPositionY[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
///   One-dimensional peak search function
///
///   This function searches for peaks in source spectrum in hin
///   The number of found peaks and their positions are written into
///   the members fNpeaks and fPositionX.
///   The search is performed in the current histogram range.
///
/// #### Parameters:
///
///   - hin:       pointer to the histogram of source spectrum
///   - sigma:   sigma of searched peaks, for details we refer to manual
///   - threshold: (default=0.05)  peaks with amplitude less than
///       threshold*highest_peak are discarded.  0<threshold<1
///
///   By default, the background is removed before deconvolution.
///   Specify the option "nobackground" to not remove the background.
///
///   By default the "Markov" chain algorithm is used.
///   Specify the option "noMarkov" to disable this algorithm
///   Note that by default the source spectrum is replaced by a new spectrum
///
///   By default a polymarker object is created and added to the list of
///   functions of the histogram. The histogram is drawn with the specified
///   option and the polymarker object drawn on top of the histogram.
///   The polymarker coordinates correspond to the npeaks peaks found in
///   the histogram.
///
///   A pointer to the polymarker object can be retrieved later via:
/// ~~~ {.cpp}
///    TList *functions = hin->GetListOfFunctions();
///    TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
/// ~~~
///   Specify the option "goff" to disable the storage and drawing of the
///   polymarker.
///
///   To disable the final drawing of the histogram with the search results (in case
///   you want to draw it yourself) specify "nodraw" in the options parameter.

Int_t TSpectrum::Search(const TH1 * hin, Double_t sigma, Option_t * option,
                        Double_t threshold)
{
   if (hin == 0) return 0;
   Int_t dimension = hin->GetDimension();
   if (dimension > 2) {
      Error("Search", "Only implemented for 1-d and 2-d histograms");
      return 0;
   }
   if (threshold <=0 || threshold >= 1) {
      Warning("Search","threshold must 0<threshold<1, threshold=0.05 assumed");
      threshold = 0.05;
   }
   TString opt = option;
   opt.ToLower();
   Bool_t background = kTRUE;
   if (opt.Contains("nobackground")) {
      background = kFALSE;
      opt.ReplaceAll("nobackground","");
   }
   Bool_t markov = kTRUE;
   if (opt.Contains("nomarkov")) {
      markov = kFALSE;
      opt.ReplaceAll("nomarkov","");
   }
   Bool_t draw = kTRUE;
   if (opt.Contains("nodraw")) {
      draw = kFALSE;
      opt.ReplaceAll("nodraw","");
   }
   if (dimension == 1) {
      Int_t first = hin->GetXaxis()->GetFirst();
      Int_t last  = hin->GetXaxis()->GetLast();
      Int_t size = last-first+1;
      Int_t i, bin, npeaks;
      Double_t * source = new Double_t[size];
      Double_t * dest   = new Double_t[size];
      for (i = 0; i < size; i++) source[i] = hin->GetBinContent(i + first);
      if (sigma < 1) {
         sigma = size/fMaxPeaks;
         if (sigma < 1) sigma = 1;
         if (sigma > 8) sigma = 8;
      }
      npeaks = SearchHighRes(source, dest, size, sigma, 100*threshold,
                             background, fgIterations, markov, fgAverageWindow);

      for (i = 0; i < npeaks; i++) {
         bin = first + Int_t(fPositionX[i] + 0.5);
         fPositionX[i] = hin->GetBinCenter(bin);
         fPositionY[i] = hin->GetBinContent(bin);
      }
      delete [] source;
      delete [] dest;

      if (opt.Contains("goff"))
         return npeaks;
      if (!npeaks) return 0;
      TPolyMarker * pm =
         (TPolyMarker*)hin->GetListOfFunctions()->FindObject("TPolyMarker");
      if (pm) {
         hin->GetListOfFunctions()->Remove(pm);
         delete pm;
      }
      pm = new TPolyMarker(npeaks, fPositionX, fPositionY);
      hin->GetListOfFunctions()->Add(pm);
      pm->SetMarkerStyle(23);
      pm->SetMarkerColor(kRed);
      pm->SetMarkerSize(1.3);
      opt.ReplaceAll(" ","");
      opt.ReplaceAll(",","");
      if (draw)
         ((TH1*)hin)->Draw(opt.Data());
      return npeaks;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// *NOT USED*
///  resolution: determines resolution of the neighbouring peaks
///              default value is 1 correspond to 3 sigma distance
///              between peaks. Higher values allow higher resolution
///              (smaller distance between peaks.
///              May be set later through SetResolution.

void TSpectrum::SetResolution(Double_t resolution)
{
   if (resolution > 1)
      fResolution = resolution;
   else
      fResolution = 1;
}

////////////////////////////////////////////////////////////////////////////////
/// This function calculates background spectrum from source spectrum.
/// The result is placed in the vector pointed by spe1945ctrum pointer.
/// The goal is to separate the useful information (peaks) from useless
/// information (background).
///
/// - method is based on Sensitive Nonlinear Iterative Peak (SNIP) clipping
///      algorithm.
/// - new value in the channel "i" is calculated
///
/// \f[
/// v_p(i) = min \left\{ v_{p-1}(i)^{\frac{\left[v_{p-1}(i+p)+v_{p-1}(i-p)\right]}{2}}   \right\}
/// \f]
///
/// where p = 1, 2, ..., numberIterations. In fact it represents second order
/// difference filter (-1,2,-1).
///
/// One can also change the
/// direction of the change of the clipping window, the order of the clipping
/// filter, to include smoothing, to set width of smoothing window and to include
/// the estimation of Compton edges. On successful completion it returns 0. On
/// error it returns pointer to the string describing error.
///
/// #### Parameters:
///
/// - spectrum: pointer to the vector of source spectrum
/// - ssize: length of the spectrum vector
/// - numberIterations: maximal width of clipping window,
/// - direction:  direction of change of clipping window.
///      Possible values: kBackIncreasingWindow, kBackDecreasingWindow
/// - filterOrder: order of clipping filter.
///      Possible values: kBackOrder2, kBackOrder4, kBackOrder6, kBackOrder8
/// - smoothing: logical variable whether the smoothing operation in the
///      estimation of background will be included.
///      Possible values: kFALSE, kTRUE
/// - smoothWindow: width of smoothing window.
///      Possible values: kBackSmoothing3, kBackSmoothing5, kBackSmoothing7,
///      kBackSmoothing9, kBackSmoothing11, kBackSmoothing13, kBackSmoothing15.
/// - compton: logical variable whether the estimation of Compton edge will be
///      included. Possible values: kFALSE, kTRUE.
///
/// #### References:
///
///   1. C. G Ryan et al.: SNIP, a statistics-sensitive background treatment for the
/// quantitative analysis of PIXE spectra in geoscience applications. NIM, B34
/// (1988), 396-402.
///
///   2. M. Morhac;, J. Kliman, V. Matouoek, M. Veselsky, I. Turzo:
/// Background elimination methods for multidimensional gamma-ray spectra. NIM,
/// A401 (1997) 113-132.
///
///   3. D. D. Burgess, R. J. Tervo: Background estimation for gamma-ray
/// spectroscopy. NIM 214 (1983), 431-434.
///
/// ### Example 1 script Background_incr.C:
///
/// Example of the estimation of background for number of iterations=6.
/// Original spectrum is shown in black color, estimated background in red color.
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/Background_incr.C
/// End_Macro
///
/// ### Example 2 script Background_decr.C:
///
/// In Example 1. one can notice that at the edges of the peaks the estimated
/// background goes under the peaks. An alternative approach is to decrease the
/// clipping window from a given value numberIterations to the value of one, which
/// is presented in this example.
///
/// Example of the estimation of background for numberIterations=6 using
/// decreasing clipping window algorithm. Original spectrum is shown in black
/// color, estimated background in red color.
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/Background_decr.C
/// End_Macro
///
/// ### Example 3 script Background_width.C:
///
/// The question is how to choose the width of the clipping window, i.e.,
/// numberIterations parameter. The influence of this parameter on the estimated
/// background is illustrated in Example 3.
///
/// Example of the influence of clipping window width on the estimated background
/// for numberIterations=4 (red line), 6 (orange line) 8 (green line) using decreasing
/// clipping window algorithm.
///
/// in general one should set this parameter so that the value
/// 2*numberIterations+1 was greater than the widths of preserved objects (peaks).
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/Background_width.C
/// End_Macro
///
/// ### Example 4 script Background_width2.C:
///
/// another example for very complex spectrum is given here.
///
/// Example of the influence of clipping window width on the estimated background
/// for numberIterations=10 (red line), 20 (blue line), 30 (green line) and
/// 40 (magenta line) using decreasing clipping window algorithm.
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/Background_width2.C
/// End_Macro
///
/// ### Example 5 script Background_order.C:
///
/// Second order difference filter removes linear (quasi-linear) background and
/// preserves symmetrical peaks. However if the shape of the background is more
/// complex one can employ higher-order clipping filters.
///
/// Example of the influence of clipping filter difference order on the estimated
/// background for fNnumberIterations=40, 2-nd order red line, 4-th order blue line,
/// 6-th order green line and 8-th order magenta line, and using decreasing
/// clipping window algorithm.
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/Background_order.C
/// End_Macro
///
/// ### Example 6 script Background_smooth.C:
///
/// The estimate of the background can be influenced by noise present in the
/// spectrum.  We proposed the algorithm of the background estimate with
/// simultaneous smoothing.  In the original algorithm without smoothing, the
/// estimated background snatches the lower spikes in the noise. Consequently,
/// the areas of peaks are biased by this error.
///
/// \image html TSpectrum_Background_smooth1.jpg Principle of background estimation algorithm with simultaneous smoothing.
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/Background_smooth.C
/// End_Macro
///
/// ### Example 8 script Background_compton.C:
///
/// Sometimes it is necessary to include also the Compton edges into the estimate of
/// the background. This example presents the synthetic spectrum
/// with Compton edges. The background was estimated using the 8-th order filter
/// with the estimation of the Compton edges using decreasing
/// clipping window algorithm (numberIterations=10) with smoothing
/// (smoothingWindow=5).
///
/// Example of the estimate of the background with Compton edges (red line) for
/// numberIterations=10, 8-th order difference filter, using decreasing clipping
/// window algorithm and smoothing (smoothingWindow=5).
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/Background_compton.C
/// End_Macro

const char *TSpectrum::Background(Double_t *spectrum, int ssize,
                                          int numberIterations,
                                          int direction, int filterOrder,
                                          bool smoothing,int smoothWindow,
                                          bool compton)
{
   int i, j, w, bw, b1, b2, priz;
   Double_t a, b, c, d, e, yb1, yb2, ai, av, men, b4, c4, d4, e4, b6, c6, d6, e6, f6, g6, b8, c8, d8, e8, f8, g8, h8, i8;
   if (ssize <= 0)
      return "Wrong Parameters";
   if (numberIterations < 1)
      return "Width of Clipping Window Must Be Positive";
   if (ssize < 2 * numberIterations + 1)
      return "Too Large Clipping Window";
   if (smoothing == kTRUE && smoothWindow != kBackSmoothing3 && smoothWindow != kBackSmoothing5 && smoothWindow != kBackSmoothing7 && smoothWindow != kBackSmoothing9 && smoothWindow != kBackSmoothing11 && smoothWindow != kBackSmoothing13 && smoothWindow != kBackSmoothing15)
      return "Incorrect width of smoothing window";
   Double_t *working_space = new Double_t[2 * ssize];
   for (i = 0; i < ssize; i++){
      working_space[i] = spectrum[i];
      working_space[i + ssize] = spectrum[i];
   }
   bw=(smoothWindow-1)/2;
   if (direction == kBackIncreasingWindow)
      i = 1;
   else if(direction == kBackDecreasingWindow)
      i = numberIterations;
   if (filterOrder == kBackOrder2) {
      do{
         for (j = i; j < ssize - i; j++) {
            if (smoothing == kFALSE){
               a = working_space[ssize + j];
               b = (working_space[ssize + j - i] + working_space[ssize + j + i]) / 2.0;
               if (b < a)
                  a = b;
               working_space[j] = a;
            }

            else if (smoothing == kTRUE){
               a = working_space[ssize + j];
               av = 0;
               men = 0;
               for (w = j - bw; w <= j + bw; w++){
                  if ( w >= 0 && w < ssize){
                     av += working_space[ssize + w];
                     men +=1;
                  }
               }
               av = av / men;
               b = 0;
               men = 0;
               for (w = j - i - bw; w <= j - i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     b += working_space[ssize + w];
                     men +=1;
                  }
               }
               b = b / men;
               c = 0;
               men = 0;
               for (w = j + i - bw; w <= j + i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     c += working_space[ssize + w];
                     men +=1;
                  }
               }
               c = c / men;
               b = (b + c) / 2;
               if (b < a)
                  av = b;
               working_space[j]=av;
            }
         }
         for (j = i; j < ssize - i; j++)
            working_space[ssize + j] = working_space[j];
         if (direction == kBackIncreasingWindow)
            i+=1;
         else if(direction == kBackDecreasingWindow)
            i-=1;
      }while((direction == kBackIncreasingWindow && i <= numberIterations) || (direction == kBackDecreasingWindow && i >= 1));
   }

   else if (filterOrder == kBackOrder4) {
      do{
         for (j = i; j < ssize - i; j++) {
            if (smoothing == kFALSE){
               a = working_space[ssize + j];
               b = (working_space[ssize + j - i] + working_space[ssize + j + i]) / 2.0;
               c = 0;
               ai = i / 2;
               c -= working_space[ssize + j - (Int_t) (2 * ai)] / 6;
               c += 4 * working_space[ssize + j - (Int_t) ai] / 6;
               c += 4 * working_space[ssize + j + (Int_t) ai] / 6;
               c -= working_space[ssize + j + (Int_t) (2 * ai)] / 6;
               if (b < c)
                  b = c;
               if (b < a)
                  a = b;
               working_space[j] = a;
            }

            else if (smoothing == kTRUE){
               a = working_space[ssize + j];
               av = 0;
               men = 0;
               for (w = j - bw; w <= j + bw; w++){
                  if ( w >= 0 && w < ssize){
                     av += working_space[ssize + w];
                     men +=1;
                  }
               }
               av = av / men;
               b = 0;
               men = 0;
               for (w = j - i - bw; w <= j - i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     b += working_space[ssize + w];
                     men +=1;
                  }
               }
               b = b / men;
               c = 0;
               men = 0;
               for (w = j + i - bw; w <= j + i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     c += working_space[ssize + w];
                     men +=1;
                  }
               }
               c = c / men;
               b = (b + c) / 2;
               ai = i / 2;
               b4 = 0, men = 0;
               for (w = j - (Int_t)(2 * ai) - bw; w <= j - (Int_t)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     b4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               b4 = b4 / men;
               c4 = 0, men = 0;
               for (w = j - (Int_t)ai - bw; w <= j - (Int_t)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     c4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               c4 = c4 / men;
               d4 = 0, men = 0;
               for (w = j + (Int_t)ai - bw; w <= j + (Int_t)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     d4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               d4 = d4 / men;
               e4 = 0, men = 0;
               for (w = j + (Int_t)(2 * ai) - bw; w <= j + (Int_t)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     e4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               e4 = e4 / men;
               b4 = (-b4 + 4 * c4 + 4 * d4 - e4) / 6;
               if (b < b4)
                  b = b4;
               if (b < a)
                  av = b;
               working_space[j]=av;
            }
         }
         for (j = i; j < ssize - i; j++)
            working_space[ssize + j] = working_space[j];
         if (direction == kBackIncreasingWindow)
            i+=1;
         else if(direction == kBackDecreasingWindow)
            i-=1;
      }while((direction == kBackIncreasingWindow && i <= numberIterations) || (direction == kBackDecreasingWindow && i >= 1));
   }

   else if (filterOrder == kBackOrder6) {
      do{
         for (j = i; j < ssize - i; j++) {
            if (smoothing == kFALSE){
               a = working_space[ssize + j];
               b = (working_space[ssize + j - i] + working_space[ssize + j + i]) / 2.0;
               c = 0;
               ai = i / 2;
               c -= working_space[ssize + j - (Int_t) (2 * ai)] / 6;
               c += 4 * working_space[ssize + j - (Int_t) ai] / 6;
               c += 4 * working_space[ssize + j + (Int_t) ai] / 6;
               c -= working_space[ssize + j + (Int_t) (2 * ai)] / 6;
               d = 0;
               ai = i / 3;
               d += working_space[ssize + j - (Int_t) (3 * ai)] / 20;
               d -= 6 * working_space[ssize + j - (Int_t) (2 * ai)] / 20;
               d += 15 * working_space[ssize + j - (Int_t) ai] / 20;
               d += 15 * working_space[ssize + j + (Int_t) ai] / 20;
               d -= 6 * working_space[ssize + j + (Int_t) (2 * ai)] / 20;
               d += working_space[ssize + j + (Int_t) (3 * ai)] / 20;
               if (b < d)
                  b = d;
               if (b < c)
                  b = c;
               if (b < a)
                  a = b;
               working_space[j] = a;
            }

            else if (smoothing == kTRUE){
               a = working_space[ssize + j];
               av = 0;
               men = 0;
               for (w = j - bw; w <= j + bw; w++){
                  if ( w >= 0 && w < ssize){
                     av += working_space[ssize + w];
                     men +=1;
                  }
               }
               av = av / men;
               b = 0;
               men = 0;
               for (w = j - i - bw; w <= j - i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     b += working_space[ssize + w];
                     men +=1;
                  }
               }
               b = b / men;
               c = 0;
               men = 0;
               for (w = j + i - bw; w <= j + i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     c += working_space[ssize + w];
                     men +=1;
                  }
               }
               c = c / men;
               b = (b + c) / 2;
               ai = i / 2;
               b4 = 0, men = 0;
               for (w = j - (Int_t)(2 * ai) - bw; w <= j - (Int_t)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     b4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               b4 = b4 / men;
               c4 = 0, men = 0;
               for (w = j - (Int_t)ai - bw; w <= j - (Int_t)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     c4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               c4 = c4 / men;
               d4 = 0, men = 0;
               for (w = j + (Int_t)ai - bw; w <= j + (Int_t)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     d4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               d4 = d4 / men;
               e4 = 0, men = 0;
               for (w = j + (Int_t)(2 * ai) - bw; w <= j + (Int_t)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     e4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               e4 = e4 / men;
               b4 = (-b4 + 4 * c4 + 4 * d4 - e4) / 6;
               ai = i / 3;
               b6 = 0, men = 0;
               for (w = j - (Int_t)(3 * ai) - bw; w <= j - (Int_t)(3 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     b6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               b6 = b6 / men;
               c6 = 0, men = 0;
               for (w = j - (Int_t)(2 * ai) - bw; w <= j - (Int_t)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     c6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               c6 = c6 / men;
               d6 = 0, men = 0;
               for (w = j - (Int_t)ai - bw; w <= j - (Int_t)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     d6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               d6 = d6 / men;
               e6 = 0, men = 0;
               for (w = j + (Int_t)ai - bw; w <= j + (Int_t)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     e6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               e6 = e6 / men;
               f6 = 0, men = 0;
               for (w = j + (Int_t)(2 * ai) - bw; w <= j + (Int_t)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     f6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               f6 = f6 / men;
               g6 = 0, men = 0;
               for (w = j + (Int_t)(3 * ai) - bw; w <= j + (Int_t)(3 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     g6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               g6 = g6 / men;
               b6 = (b6 - 6 * c6 + 15 * d6 + 15 * e6 - 6 * f6 + g6) / 20;
               if (b < b6)
                  b = b6;
               if (b < b4)
                  b = b4;
               if (b < a)
                  av = b;
               working_space[j]=av;
            }
         }
         for (j = i; j < ssize - i; j++)
            working_space[ssize + j] = working_space[j];
         if (direction == kBackIncreasingWindow)
            i+=1;
         else if(direction == kBackDecreasingWindow)
            i-=1;
      }while((direction == kBackIncreasingWindow && i <= numberIterations) || (direction == kBackDecreasingWindow && i >= 1));
   }

   else if (filterOrder == kBackOrder8) {
      do{
         for (j = i; j < ssize - i; j++) {
            if (smoothing == kFALSE){
               a = working_space[ssize + j];
               b = (working_space[ssize + j - i] + working_space[ssize + j + i]) / 2.0;
               c = 0;
               ai = i / 2;
               c -= working_space[ssize + j - (Int_t) (2 * ai)] / 6;
               c += 4 * working_space[ssize + j - (Int_t) ai] / 6;
               c += 4 * working_space[ssize + j + (Int_t) ai] / 6;
               c -= working_space[ssize + j + (Int_t) (2 * ai)] / 6;
               d = 0;
               ai = i / 3;
               d += working_space[ssize + j - (Int_t) (3 * ai)] / 20;
               d -= 6 * working_space[ssize + j - (Int_t) (2 * ai)] / 20;
               d += 15 * working_space[ssize + j - (Int_t) ai] / 20;
               d += 15 * working_space[ssize + j + (Int_t) ai] / 20;
               d -= 6 * working_space[ssize + j + (Int_t) (2 * ai)] / 20;
               d += working_space[ssize + j + (Int_t) (3 * ai)] / 20;
               e = 0;
               ai = i / 4;
               e -= working_space[ssize + j - (Int_t) (4 * ai)] / 70;
               e += 8 * working_space[ssize + j - (Int_t) (3 * ai)] / 70;
               e -= 28 * working_space[ssize + j - (Int_t) (2 * ai)] / 70;
               e += 56 * working_space[ssize + j - (Int_t) ai] / 70;
               e += 56 * working_space[ssize + j + (Int_t) ai] / 70;
               e -= 28 * working_space[ssize + j + (Int_t) (2 * ai)] / 70;
               e += 8 * working_space[ssize + j + (Int_t) (3 * ai)] / 70;
               e -= working_space[ssize + j + (Int_t) (4 * ai)] / 70;
               if (b < e)
                  b = e;
               if (b < d)
                  b = d;
               if (b < c)
                  b = c;
               if (b < a)
                  a = b;
               working_space[j] = a;
            }

            else if (smoothing == kTRUE){
               a = working_space[ssize + j];
               av = 0;
               men = 0;
               for (w = j - bw; w <= j + bw; w++){
                  if ( w >= 0 && w < ssize){
                     av += working_space[ssize + w];
                     men +=1;
                  }
               }
               av = av / men;
               b = 0;
               men = 0;
               for (w = j - i - bw; w <= j - i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     b += working_space[ssize + w];
                     men +=1;
                  }
               }
               b = b / men;
               c = 0;
               men = 0;
               for (w = j + i - bw; w <= j + i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     c += working_space[ssize + w];
                     men +=1;
                  }
               }
               c = c / men;
               b = (b + c) / 2;
               ai = i / 2;
               b4 = 0, men = 0;
               for (w = j - (Int_t)(2 * ai) - bw; w <= j - (Int_t)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     b4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               b4 = b4 / men;
               c4 = 0, men = 0;
               for (w = j - (Int_t)ai - bw; w <= j - (Int_t)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     c4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               c4 = c4 / men;
               d4 = 0, men = 0;
               for (w = j + (Int_t)ai - bw; w <= j + (Int_t)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     d4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               d4 = d4 / men;
               e4 = 0, men = 0;
               for (w = j + (Int_t)(2 * ai) - bw; w <= j + (Int_t)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     e4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               e4 = e4 / men;
               b4 = (-b4 + 4 * c4 + 4 * d4 - e4) / 6;
               ai = i / 3;
               b6 = 0, men = 0;
               for (w = j - (Int_t)(3 * ai) - bw; w <= j - (Int_t)(3 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     b6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               b6 = b6 / men;
               c6 = 0, men = 0;
               for (w = j - (Int_t)(2 * ai) - bw; w <= j - (Int_t)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     c6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               c6 = c6 / men;
               d6 = 0, men = 0;
               for (w = j - (Int_t)ai - bw; w <= j - (Int_t)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     d6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               d6 = d6 / men;
               e6 = 0, men = 0;
               for (w = j + (Int_t)ai - bw; w <= j + (Int_t)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     e6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               e6 = e6 / men;
               f6 = 0, men = 0;
               for (w = j + (Int_t)(2 * ai) - bw; w <= j + (Int_t)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     f6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               f6 = f6 / men;
               g6 = 0, men = 0;
               for (w = j + (Int_t)(3 * ai) - bw; w <= j + (Int_t)(3 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     g6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               g6 = g6 / men;
               b6 = (b6 - 6 * c6 + 15 * d6 + 15 * e6 - 6 * f6 + g6) / 20;
               ai = i / 4;
               b8 = 0, men = 0;
               for (w = j - (Int_t)(4 * ai) - bw; w <= j - (Int_t)(4 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     b8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               b8 = b8 / men;
               c8 = 0, men = 0;
               for (w = j - (Int_t)(3 * ai) - bw; w <= j - (Int_t)(3 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     c8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               c8 = c8 / men;
               d8 = 0, men = 0;
               for (w = j - (Int_t)(2 * ai) - bw; w <= j - (Int_t)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     d8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               d8 = d8 / men;
               e8 = 0, men = 0;
               for (w = j - (Int_t)ai - bw; w <= j - (Int_t)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     e8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               e8 = e8 / men;
               f8 = 0, men = 0;
               for (w = j + (Int_t)ai - bw; w <= j + (Int_t)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     f8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               f8 = f8 / men;
               g8 = 0, men = 0;
               for (w = j + (Int_t)(2 * ai) - bw; w <= j + (Int_t)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     g8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               g8 = g8 / men;
               h8 = 0, men = 0;
               for (w = j + (Int_t)(3 * ai) - bw; w <= j + (Int_t)(3 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     h8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               h8 = h8 / men;
               i8 = 0, men = 0;
               for (w = j + (Int_t)(4 * ai) - bw; w <= j + (Int_t)(4 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     i8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               i8 = i8 / men;
               b8 = ( -b8 + 8 * c8 - 28 * d8 + 56 * e8 - 56 * f8 - 28 * g8 + 8 * h8 - i8)/70;
               if (b < b8)
                  b = b8;
               if (b < b6)
                  b = b6;
               if (b < b4)
                  b = b4;
               if (b < a)
                  av = b;
               working_space[j]=av;
            }
         }
         for (j = i; j < ssize - i; j++)
            working_space[ssize + j] = working_space[j];
         if (direction == kBackIncreasingWindow)
            i += 1;
         else if(direction == kBackDecreasingWindow)
            i -= 1;
      }while((direction == kBackIncreasingWindow && i <= numberIterations) || (direction == kBackDecreasingWindow && i >= 1));
   }

   if (compton == kTRUE) {
      for (i = 0, b2 = 0; i < ssize; i++){
         b1 = b2;
         a = working_space[i], b = spectrum[i];
         j = i;
         if (TMath::Abs(a - b) >= 1) {
            b1 = i - 1;
            if (b1 < 0)
               b1 = 0;
            yb1 = working_space[b1];
            for (b2 = b1 + 1, c = 0, priz = 0; priz == 0 && b2 < ssize; b2++){
               a = working_space[b2], b = spectrum[b2];
               c = c + b - yb1;
               if (TMath::Abs(a - b) < 1) {
                  priz = 1;
                  yb2 = b;
               }
            }
            if (b2 == ssize)
               b2 -= 1;
            yb2 = working_space[b2];
            if (yb1 <= yb2){
               for (j = b1, c = 0; j <= b2; j++){
                  b = spectrum[j];
                  c = c + b - yb1;
               }
               if (c > 1){
                  c = (yb2 - yb1) / c;
                  for (j = b1, d = 0; j <= b2 && j < ssize; j++){
                     b = spectrum[j];
                     d = d + b - yb1;
                     a = c * d + yb1;
                     working_space[ssize + j] = a;
                  }
               }
            }

            else{
               for (j = b2, c = 0; j >= b1; j--){
                  b = spectrum[j];
                  c = c + b - yb2;
               }
               if (c > 1){
                  c = (yb1 - yb2) / c;
                  for (j = b2, d = 0;j >= b1 && j >= 0; j--){
                     b = spectrum[j];
                     d = d + b - yb2;
                     a = c * d + yb2;
                     working_space[ssize + j] = a;
                  }
               }
            }
            i=b2;
         }
      }
   }

   for (j = 0; j < ssize; j++)
      spectrum[j] = working_space[ssize + j];
   delete[]working_space;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// One-dimensional markov spectrum smoothing function
///
/// This function calculates smoothed spectrum from source spectrum based on
/// Markov chain method. The result is placed in the array pointed by source
/// pointer. On successful completion it returns 0. On error it returns pointer
/// to the string describing error.
///
/// #### Parameters:
///
///  -  source: pointer to the array of source spectrum
///  -  ssize: length of source array
///  -  averWindow: width of averaging smoothing window
///
/// The goal of this function is the suppression of the statistical fluctuations.
/// The algorithm is based on discrete Markov chain, which has very simple
/// invariant distribution:
///
/// \f[
/// U_2 = \frac{p_{1,2}}{p_{2,1}}U_1, U_3 = \frac{p_{2,3}}{p_{3,2}}U_2U_1, ... , U_n = \frac{p_{n-1,n}}{p_{n,n-1}}U_{n-1}...U_2U_1
/// \f]
/// \f$ U_1\f$ being defined from the normalization condition
/// \f$ \sum_{i=1}^{n} U_i=1\f$. \f$ n \f$ is the length of the smoothed spectrum and
/// \f[
/// p_{i,i\pm 1} = A_i\sum_{k=1}^{m} exp\left[ \frac{y(i\pm k)-y(i)}{y(i\pm k)+y(i)}\right]
/// \f]
///
/// #### Reference:
///
///  1. Z.K. Silagadze, A new algorithm for automatic photopeak searches.
/// NIM A 376 (1996), 451.
///
/// ### Example 14 - script Smoothing.C
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/Smoothing.C
/// End_Macro

const char* TSpectrum::SmoothMarkov(Double_t *source, int ssize, int averWindow)
{
   int xmin, xmax, i, l;
   Double_t a, b, maxch;
   Double_t nom, nip, nim, sp, sm, area = 0;
   if(averWindow <= 0)
      return "Averaging Window must be positive";
   Double_t *working_space = new Double_t[ssize];
   xmin = 0,xmax = ssize - 1;
   for(i = 0, maxch = 0; i < ssize; i++){
      working_space[i]=0;
      if(maxch < source[i])
         maxch = source[i];

      area += source[i];
   }
   if(maxch == 0) {
      delete [] working_space;
      return 0 ;
   }

   nom = 1;
   working_space[xmin] = 1;
   for(i = xmin; i < xmax; i++){
      nip = source[i] / maxch;
      nim = source[i + 1] / maxch;
      sp = 0,sm = 0;
      for(l = 1; l <= averWindow; l++){
         if((i + l) > xmax)
            a = source[xmax] / maxch;

         else
            a = source[i + l] / maxch;
         b = a - nip;
         if(a + nip <= 0)
            a = 1;

         else
            a = TMath::Sqrt(a + nip);
         b = b / a;
         b = TMath::Exp(b);
         sp = sp + b;
         if((i - l + 1) < xmin)
            a = source[xmin] / maxch;

         else
            a = source[i - l + 1] / maxch;
         b = a - nim;
         if(a + nim <= 0)
            a = 1;
         else
            a = TMath::Sqrt(a + nim);
         b = b / a;
         b = TMath::Exp(b);
         sm = sm + b;
      }
      a = sp / sm;
      a = working_space[i + 1] = working_space[i] * a;
      nom = nom + a;
   }
   for(i = xmin; i <= xmax; i++){
      working_space[i] = working_space[i] / nom;
   }
   for(i = 0; i < ssize; i++)
      source[i] = working_space[i] * area;
   delete[]working_space;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// One-dimensional deconvolution function
///
/// This function calculates deconvolution from source spectrum according to
/// response spectrum using Gold deconvolution algorithm. The result is placed
/// in the vector pointed by source pointer. On successful completion it
/// returns 0. On error it returns pointer to the string describing error. If
/// desired after every numberIterations one can apply boosting operation
/// (exponential function with exponent given by boost coefficient) and repeat
/// it numberRepetitions times.
///
/// #### Parameters:
///
///  - source:  pointer to the vector of source spectrum
///  - response:     pointer to the vector of response spectrum
///  - ssize:    length of source and response spectra
///  - numberIterations, for details we refer to the reference given below
///  - numberRepetitions, for repeated boosted deconvolution
///  - boost, boosting coefficient
///
/// The goal of this function is the improvement of the resolution in spectra,
/// decomposition of multiplets. The mathematical formulation of
/// the convolution system is:
///
/// \f[
/// y(i) = \sum_{k=0}^{N-1} h(i-k)x(k), i=0,1,2,...,N-1
/// \f]
///
/// where h(i) is the impulse response function, x, y are input and output
/// vectors, respectively, N is the length of x and h vectors. In matrix form
/// we have:
/**
 \f[
 \begin{bmatrix}
     y(0)    \\
     y(1)    \\
     \dots   \\
     y(2N-2) \\
     y(2N-1)
 \end{bmatrix}
 =
 \begin{bmatrix}
     h(0)    & 0      & 0      & \dots & 0         \\
     h(1)    & h(0)   & 0      & \dots & \dots     \\
     \dots   & h(1)   & h(0)   & \dots & \dots     \\
     \dots   & \dots  & h(1)   & \dots & \dots     \\
     \dots   & \dots  & \dots  & \dots & \dots     \\
     h(N-1)  & \dots  & \dots  &\dots  & 0         \\
     0       & h(N-1) & \dots   & \dots & h(0)     \\
     0       & 0      & h(N-1)  & \dots & h(1)     \\
     \dots   & \dots  & \dots   & \dots & \dots    \\
     0       & 0      & 0       & \dots & h(N-1)
 \end{bmatrix}
 \begin{bmatrix}
     x(0)  \\
     x(1)  \\
     \dots \\
     x(N-1)
 \end{bmatrix}
 \f]
*/
/// Let us assume that we know the response and the output vector (spectrum) of
/// the above given system. The deconvolution represents solution of the
/// overdetermined system of linear equations, i.e., the calculation of the
/// vector x. From numerical stability point of view the operation of
/// deconvolution is extremely critical (ill-posed problem) as well as time
/// consuming operation. The Gold deconvolution algorithm proves to work very
/// well, other methods (Fourier, VanCittert etc) oscillate. It is suitable to
/// process positive definite data (e.g. histograms).
///
/// #### Gold deconvolution algorithm:
/**
 \f[
 y = Hx \\
 H^{T}y = H^{T}Hx \\
 y^{'} = H^{'}x \\
 x_{i}^{(k+1)} = \frac{y_{i}^{'}}{\sum_{m=0}^{N-1}H_{im}^{'}x_{m}{(k)}}x_{i}{(k)},  i=0,1,2,...,N-1 \\
 k = 1,2,3,...,L
 x^{0} = [1,1, ..., 1]^T
 \f]
*/
/// Where L is given number of iterations (numberIterations parameter).
///
/// #### Boosted deconvolution:
///
///  1. Set the initial solution:
///       \f$ x^{(0)} = [1,1,...,1]^{T} \f$
///  2. Set required number of repetitions R and iterations L.
///  3. Set r = 1.
///  4. Using Gold deconvolution algorithm for k=1,2,...,L find
///      \f$ x^{(L)} \f$
///  5. If r = R stop calculation, else
///
///    1. Apply boosting operation, i.e., set
///        \f$ x^{(0)}(i) = [x^{(L)}(i)]^{p} \f$
///       i=0,1,...N-1 and p is boosting coefficient >0.
///    2. r = r + 1
///    3. continue in 4.
///
/// #### References:
///
///  1. Gold R., ANL-6984, Argonne National Laboratories, Argonne Ill, 1964.
///  2. Coote G.E., Iterative smoothing and deconvolution of one- and two-dimensional
///     elemental distribution data, NIM B 130 (1997) 118.
///  3. M. Morhac;, J. Kliman, V.  Matouoek, M. Veselsky,
///     I. Turzo: Efficient one- and two-dimensional Gold deconvolution and
///     its application to gamma-ray spectra decomposition. NIM, A401 (1997) 385-408.
///  4. Morhac; M., Matouoek V., Kliman J., Efficient algorithm of multidimensional
///     deconvolution and its application to nuclear data processing, Digital Signal
///     Processing 13 (2003) 144.
///
/// ### Example 8 - script Deconvolution.C :
///
/// response function (usually peak) should be shifted left to the first
/// non-zero channel (bin).
///
/// \image html TSpectrum_Deconvolution2.jpg Principle how the response matrix is composed inside of the Deconvolution function.
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/Deconvolution.C
/// End_Macro
///
/// ### Examples of Gold deconvolution method:
///
/// First let us study the influence of the number of iterations on the
/// deconvolved spectrum (Fig. 12).
///
/// \image html TSpectrum_Deconvolution_wide1.jpg Fig. 12 Study of Gold deconvolution algorithm.The original source spectrum is drawn with black color, spectrum after 100 iterations with red color, spectrum after 1000 iterations with blue color, spectrum after 10000 iterations with green color and spectrum after 100000 iterations with magenta color.
///
/// For relatively narrow peaks in the above given example the Gold
/// deconvolution method is able to decompose overlapping peaks practically to
/// delta - functions. In the next example we have chosen a synthetic data
/// (spectrum, 256 channels) consisting of 5 very closely positioned, relatively
/// wide peaks (sigma =5), with added noise (Fig. 13). Thin lines represent
/// pure Gaussians (see Table 1); thick line is a resulting spectrum with
/// additive noise (10% of the amplitude of small peaks).
///
/// \image html TSpectrum_Deconvolution_wide2.jpg Fig. 13 Testing example of synthetic spectrum composed of 5 Gaussians with added noise.
///
///   | Peak #   | Position | Height | Area   |
///   |----------|----------|--------|--------|
///   | 1        | 50       | 500    | 10159  |
///   | 2        | 70       | 3000   | 60957  |
///   | 3        | 80       | 1000   | 20319  |
///   | 4        | 100      | 5000   | 101596 |
///   | 5        | 110      | 500    | 10159  |
///
/// Table 1 Positions, heights and areas of peaks in the spectrum shown in Fig. 13.
///
/// In ideal case, we should obtain the result given in Fig. 14. The areas of
/// the Gaussian components of the spectrum are concentrated completely to
/// delta-functions. When solving the overdetermined system of linear equations
/// with data from Fig. 13 in the sense of minimum least squares criterion
/// without any regularisation we obtain the result with large oscillations
/// (Fig. 15). From mathematical point of view, it is the optimal solution in
/// the unconstrained space of independent variables. From physical point of
/// view we are interested only in a meaningful solution. Therefore, we have to
/// employ regularisation techniques (e.g. Gold deconvolution) and/or to
/// confine the space of allowed solutions to subspace of positive solutions.
///
/// \image html TSpectrum_Deconvolution_wide3.jpg Fig. 14 The same spectrum like in Fig. 13, outlined bars show the contents of present components (peaks).
/// \image html TSpectrum_Deconvolution_wide4.jpg Fig. 15 Least squares solution of the system of linear equations without regularisation.
///
/// ### Example 9 - script Deconvolution_wide.C
///
/// When we employ Gold deconvolution algorithm we obtain the result given in
/// Fig. 16. One can observe that the resulting spectrum is smooth. On the
/// other hand the method is not able to decompose completely the peaks in the
/// spectrum.
///
/// Example of Gold deconvolution for closely positioned wide peaks. The original
/// source spectrum is drawn with black color, the spectrum after the deconvolution
/// (10000 iterations) with red color.
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/Deconvolution_wide.C
/// End_Macro
///
/// ### Example 10 - script Deconvolution_wide_boost.C :
///
/// Further let us employ boosting operation into deconvolution (Fig. 17).
///
/// The original source spectrum is drawn with black color, the spectrum after
/// the deconvolution with red color. Number of iterations = 200, number of
/// repetitions = 50 and boosting coefficient = 1.2.
///
/// One can observe that peaks are decomposed practically to delta functions.
/// Number of peaks is correct, positions of big peaks as well as their areas
/// are relatively well estimated. However there is a considerable error in
/// the estimation of the position of small right hand peak.
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/Deconvolution_wide_boost.C
/// End_Macro

const char *TSpectrum::Deconvolution(Double_t *source, const Double_t *response,
                                      int ssize, int numberIterations,
                                      int numberRepetitions, Double_t boost )
{
   if (ssize <= 0)
      return "Wrong Parameters";

   if (numberRepetitions <= 0)
      return "Wrong Parameters";

       //   working_space-pointer to the working vector
       //   (its size must be 4*ssize of source spectrum)
   Double_t *working_space = new Double_t[4 * ssize];
   int i, j, k, lindex, posit, lh_gold, l, repet;
   Double_t lda, ldb, ldc, area, maximum;
   area = 0;
   lh_gold = -1;
   posit = 0;
   maximum = 0;

//read response vector
   for (i = 0; i < ssize; i++) {
      lda = response[i];
      if (lda != 0)
         lh_gold = i + 1;
      working_space[i] = lda;
      area += lda;
      if (lda > maximum) {
         maximum = lda;
         posit = i;
      }
   }
   if (lh_gold == -1) {
      delete [] working_space;
      return "ZERO RESPONSE VECTOR";
   }

//read source vector
   for (i = 0; i < ssize; i++)
      working_space[2 * ssize + i] = source[i];

// create matrix at*a and vector at*y
   for (i = 0; i < ssize; i++){
      lda = 0;
      for (j = 0; j < ssize; j++){
         ldb = working_space[j];
         k = i + j;
         if (k < ssize){
            ldc = working_space[k];
            lda = lda + ldb * ldc;
         }
      }
      working_space[ssize + i] = lda;
      lda = 0;
      for (k = 0; k < ssize; k++){
         l = k - i;
         if (l >= 0){
            ldb = working_space[l];
            ldc = working_space[2 * ssize + k];
            lda = lda + ldb * ldc;
         }
      }
      working_space[3 * ssize + i]=lda;
   }

// move vector at*y
   for (i = 0; i < ssize; i++){
      working_space[2 * ssize + i] = working_space[3 * ssize + i];
   }

//initialization of resulting vector
   for (i = 0; i < ssize; i++)
      working_space[i] = 1;

       //**START OF ITERATIONS**
   for (repet = 0; repet < numberRepetitions; repet++) {
      if (repet != 0) {
         for (i = 0; i < ssize; i++)
            working_space[i] = TMath::Power(working_space[i], boost);
      }
      for (lindex = 0; lindex < numberIterations; lindex++) {
         for (i = 0; i < ssize; i++) {
            if (working_space[2 * ssize + i] > 0.000001
                 && working_space[i] > 0.000001) {
               lda = 0;
               for (j = 0; j < lh_gold; j++) {
                  ldb = working_space[j + ssize];
                  if (j != 0){
                     k = i + j;
                     ldc = 0;
                     if (k < ssize)
                        ldc = working_space[k];
                     k = i - j;
                     if (k >= 0)
                        ldc += working_space[k];
                  }

                  else
                     ldc = working_space[i];
                  lda = lda + ldb * ldc;
               }
               ldb = working_space[2 * ssize + i];
               if (lda != 0)
                  lda = ldb / lda;

               else
                  lda = 0;
               ldb = working_space[i];
               lda = lda * ldb;
               working_space[3 * ssize + i] = lda;
            }
         }
         for (i = 0; i < ssize; i++)
            working_space[i] = working_space[3 * ssize + i];
      }
   }

//shift resulting spectrum
   for (i = 0; i < ssize; i++) {
      lda = working_space[i];
      j = i + posit;
      j = j % ssize;
      working_space[ssize + j] = lda;
   }

//write back resulting spectrum
   for (i = 0; i < ssize; i++)
      source[i] = area * working_space[ssize + i];
   delete[]working_space;
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// One-dimensional deconvolution function.
///
/// This function calculates deconvolution from source spectrum according to
/// response spectrum using Richardson-Lucy deconvolution algorithm. The result
/// is placed in the vector pointed by source pointer. On successful completion
/// it returns 0. On error it returns pointer to the string describing error.
/// If desired after every numberIterations one can apply boosting operation
/// (exponential function with exponent given by boost coefficient) and repeat
/// it numberRepetitions times (see Gold deconvolution).
///
/// #### Parameters:
///
///  - source:  pointer to the vector of source spectrum
///  - response:     pointer to the vector of response spectrum
///  - ssize:    length of source and response spectra
///  - numberIterations, for details we refer to the reference given above
///  - numberRepetitions, for repeated boosted deconvolution
///  - boost, boosting coefficient
///
/// ### Richardson-Lucy deconvolution algorithm:
///
/// For discrete systems it has the form:
/**
 \f[
 x^{(n)}(i) = x^{(n-1)}(i) \sum_{j=0}^{N-1}h(i,j)\frac{y(j)}{\sum_{k=0}^{M-1}h(j,k)x^{(n-1)}(k)} \\
 i \in \left<0,M-1\right>
 \f]
*/
/// for positive input data and response matrix this iterative method forces
/// the deconvoluted spectra to be non-negative. The Richardson-Lucy
/// iteration converges to the maximum likelihood solution for Poisson statistics
/// in the data.
///
/// #### References:
///
///  1. Abreu M.C. et al., A four-dimensional deconvolution method to correct NA38
///     experimental data, NIM A 405 (1998) 139.
///  2. Lucy L.B., A.J. 79 (1974) 745.
///  3. Richardson W.H., J. Opt. Soc. Am. 62 (1972) 55.
///
/// ### Examples of Richardson-Lucy deconvolution method:
///
/// ### Example 11 - script DeconvolutionRL_wide.C :
///
/// When we employ Richardson-Lucy deconvolution algorithm to our data from
/// Fig. 13 we obtain the following result. One can observe improvements
/// as compared to the result achieved by Gold deconvolution. Nevertheless it is
/// unable to decompose the multiplet.
///
/// Example of Richardson-Lucy deconvolution for closely positioned wide peaks.
/// The original source spectrum is drawn with black color, the spectrum after
/// the deconvolution (10000 iterations) with red color.
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/DeconvolutionRL_wide.C
/// End_Macro
///
/// ### Example 12 - script DeconvolutionRL_wide_boost.C :
///
/// Further let us employ boosting operation into deconvolution.
///
/// The original source spectrum is drawn with black color, the spectrum after
/// the deconvolution with red color. Number of iterations = 200, number of
/// repetitions = 50 and boosting coefficient = 1.2.
///
/// One can observe improvements in the estimation of peak positions as compared
/// to the results achieved by Gold deconvolution.
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/DeconvolutionRL_wide_boost.C
/// End_Macro

const char *TSpectrum::DeconvolutionRL(Double_t *source, const Double_t *response,
                                      int ssize, int numberIterations,
                                      int numberRepetitions, Double_t boost )
{
   if (ssize <= 0)
      return "Wrong Parameters";

   if (numberRepetitions <= 0)
      return "Wrong Parameters";

       //   working_space-pointer to the working vector
       //   (its size must be 4*ssize of source spectrum)
   Double_t *working_space = new Double_t[4 * ssize];
   int i, j, k, lindex, posit, lh_gold, repet, kmin, kmax;
   Double_t lda, ldb, ldc, maximum;
   lh_gold = -1;
   posit = 0;
   maximum = 0;

//read response vector
   for (i = 0; i < ssize; i++) {
      lda = response[i];
      if (lda != 0)
         lh_gold = i + 1;
      working_space[ssize + i] = lda;
      if (lda > maximum) {
         maximum = lda;
         posit = i;
      }
   }
   if (lh_gold == -1) {
      delete [] working_space;
      return "ZERO RESPONSE VECTOR";
   }

//read source vector
   for (i = 0; i < ssize; i++)
      working_space[2 * ssize + i] = source[i];

//initialization of resulting vector
   for (i = 0; i < ssize; i++){
      if (i <= ssize - lh_gold)
         working_space[i] = 1;

      else
         working_space[i] = 0;

   }
       //**START OF ITERATIONS**
   for (repet = 0; repet < numberRepetitions; repet++) {
      if (repet != 0) {
         for (i = 0; i < ssize; i++)
            working_space[i] = TMath::Power(working_space[i], boost);
      }
      for (lindex = 0; lindex < numberIterations; lindex++) {
         for (i = 0; i <= ssize - lh_gold; i++){
            lda = 0;
            if (working_space[i] > 0){//x[i]
               for (j = i; j < i + lh_gold; j++){
                  ldb = working_space[2 * ssize + j];//y[j]
                  if (j < ssize){
                     if (ldb > 0){//y[j]
                        kmax = j;
                        if (kmax > lh_gold - 1)
                           kmax = lh_gold - 1;
                        kmin = j + lh_gold - ssize;
                        if (kmin < 0)
                           kmin = 0;
                        ldc = 0;
                        for (k = kmax; k >= kmin; k--){
                           ldc += working_space[ssize + k] * working_space[j - k];//h[k]*x[j-k]
                        }
                        if (ldc > 0)
                           ldb = ldb / ldc;

                        else
                           ldb = 0;
                     }
                     ldb = ldb * working_space[ssize + j - i];//y[j]*h[j-i]/suma(h[j][k]x[k])
                  }
                  lda += ldb;
               }
               lda = lda * working_space[i];
            }
            working_space[3 * ssize + i] = lda;
         }
         for (i = 0; i < ssize; i++)
            working_space[i] = working_space[3 * ssize + i];
      }
   }

//shift resulting spectrum
   for (i = 0; i < ssize; i++) {
      lda = working_space[i];
      j = i + posit;
      j = j % ssize;
      working_space[ssize + j] = lda;
   }

//write back resulting spectrum
   for (i = 0; i < ssize; i++)
      source[i] = working_space[ssize + i];
   delete[]working_space;
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// One-dimensional unfolding function
///
/// This function unfolds source spectrum according to response matrix columns.
/// The result is placed in the vector pointed by source pointer.
/// The coefficients of the resulting vector represent contents of the columns
/// (weights) in the input vector. On successful completion it returns 0. On
/// error it returns pointer to the string describing error. If desired after
/// every numberIterations one can apply boosting operation (exponential
/// function with exponent given by boost coefficient) and repeat it
/// numberRepetitions times. For details we refer to [1].
///
/// #### Parameters:
///
///  -  source: pointer to the vector of source spectrum
///  -  respMatrix: pointer to the matrix of response spectra
///  -  ssizex: length of source spectrum and # of rows of the response
///      matrix. ssizex must be >= ssizey.
///  -  ssizey: length of destination coefficients and # of columns of the response
///      matrix.
///  -  numberIterations: number of iterations
///  -  numberRepetitions: number of repetitions for boosted deconvolution.
///      It must be greater or equal to one.
///  -  boost: boosting coefficient, applies only if numberRepetitions is
///      greater than one.
///
/// ### Unfolding:
///
/// The goal is the decomposition of spectrum to a given set of component spectra.
///
/// The mathematical formulation of the discrete linear system is:
///
/// \f[
/// y(i) = \sum_{k=0}^{N_y-1} h(i,k)x(k), i = 0,1,2,...,N_x-1
/// \f]
/**
 \f[
 \begin{bmatrix}
     y(0) \\
     y(1) \\
     \dots \\
     y(N_x-1)
 \end{bmatrix}
 =
 \begin{bmatrix}
     h(0,0)       & h(0,1) & \dots & h(0,N_y-1) \\
     h(1,0)       & h(1,1) & \dots & h(1,N_y-1) \\
 \dots \\
     h(N_x-1,0)   & h(N_x-1,1) & \dots & h(N_x-1,N_y-1)
 \end{bmatrix}
 \begin{bmatrix}
     x(0) \\
     x(1) \\
     \dots \\
     x(N_y-1)
 \end{bmatrix}
 \f]
*/
///
/// #### References:
///
///  1. Jandel M., Morhac; M., Kliman J., Krupa L., Matouoek
///     V., Hamilton J. H., Ramaya A. V.:
///     Decomposition of continuum gamma-ray spectra using synthesised response matrix.
///     NIM A 516 (2004), 172-183.
///
/// ### Example of unfolding:
///
/// ### Example 13 - script Unfolding.C:
///
/// \image html TSpectrum_Unfolding3.gif Fig. 20 Response matrix composed of neutron spectra of pure chemical elements.
/// \image html TSpectrum_Unfolding2.jpg Fig. 21 Source neutron spectrum to be decomposed
/// \image html TSpectrum_Unfolding3.jpg Fig. 22 Spectrum after decomposition, contains 10 coefficients, which correspond to contents of chemical components (dominant 8-th and 10-th components, i.e. O, Si)
///
/// #### Script:
///
/// ~~~ {.cpp}
///  // Example to illustrate unfolding function (class TSpectrum).
///  // To execute this example, do
///  // root > .x Unfolding.C
///
/// void Unfolding() {
///    Int_t i, j;
///    Int_t nbinsx = 2048;
///    Int_t nbinsy = 10;
///    double xmin  = 0;
///    double xmax  = nbinsx;
///    double ymin  = 0;
///    double ymax  = nbinsy;
///    double *source    = new double[nbinsx];
///    double **response = new double *[nbinsy];
///    for (i=0;i<nbinsy;i++) response[i] = new double[nbinsx];
///    TH1F *h = new TH1F("h","",nbinsx,xmin,xmax);
///    TH1F *d = new TH1F("d","Decomposition - unfolding",nbinsx,xmin,xmax);
///    TH2F *decon_unf_resp = new TH2F("decon_unf_resp","Root File",nbinsy,ymin,ymax,nbinsx,xmin,xmax);
///    TFile *f = new TFile("TSpectrum.root");
///    h = (TH1F*) f->Get("decon_unf_in;1");
///    TFile *fr = new TFile("TSpectrum.root");
///    decon_unf_resp = (TH2F*) fr->Get("decon_unf_resp;1");
///    for (i = 0; i < nbinsx; i++) source[i] = h->GetBinContent(i + 1);
///    for (i = 0; i < nbinsy; i++){
///       for (j = 0; j< nbinsx; j++){
///          response[i][j] = decon_unf_resp->GetBinContent(i + 1, j + 1);
///       }
///    }
///    TCanvas *Decon1 = (TCanvas *)gROOT->GetListOfCanvases()->FindObject("Decon1");
///    if (!Decon1) Decon1 = new TCanvas("Decon1","Decon1",10,10,1000,700);
///    h->Draw("L");
///    TSpectrum *s = new TSpectrum();
///    s->Unfolding(source,(const double**)response,nbinsx,nbinsy,1000,1,1);
///    for (i = 0; i < nbinsy; i++) d->SetBinContent(i + 1,source[i]);
///    d->SetLineColor(kRed);
///    d->SetAxisRange(0,nbinsy);
///    d->Draw("");
/// }
/// ~~~

const char *TSpectrum::Unfolding(Double_t *source,
                                 const Double_t **respMatrix,
                                 int ssizex, int ssizey,
                                 int numberIterations,
                                 int numberRepetitions, Double_t boost)
{
   int i, j, k, lindex, lhx = 0, repet;
   Double_t lda, ldb, ldc, area;
   if (ssizex <= 0 || ssizey <= 0)
      return "Wrong Parameters";
   if (ssizex < ssizey)
      return "Sizex must be greater than sizey)";
   if (numberIterations <= 0)
      return "Number of iterations must be positive";
   Double_t *working_space =
       new Double_t[ssizex * ssizey + 2 * ssizey * ssizey + 4 * ssizex];

/*read response matrix*/
   for (j = 0; j < ssizey && lhx != -1; j++) {
      area = 0;
      lhx = -1;
      for (i = 0; i < ssizex; i++) {
         lda = respMatrix[j][i];
         if (lda != 0) {
            lhx = i + 1;
         }
         working_space[j * ssizex + i] = lda;
         area = area + lda;
      }
      if (lhx != -1) {
         for (i = 0; i < ssizex; i++)
            working_space[j * ssizex + i] /= area;
      }
   }
   if (lhx == -1) {
      delete [] working_space;
      return ("ZERO COLUMN IN RESPONSE MATRIX");
   }

/*read source vector*/
   for (i = 0; i < ssizex; i++)
      working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex + i] =
          source[i];

/*create matrix at*a + at*y */
   for (i = 0; i < ssizey; i++) {
      for (j = 0; j < ssizey; j++) {
         lda = 0;
         for (k = 0; k < ssizex; k++) {
            ldb = working_space[ssizex * i + k];
            ldc = working_space[ssizex * j + k];
            lda = lda + ldb * ldc;
         }
         working_space[ssizex * ssizey + ssizey * i + j] = lda;
      }
      lda = 0;
      for (k = 0; k < ssizex; k++) {
         ldb = working_space[ssizex * i + k];
         ldc =
             working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex +
                           k];
         lda = lda + ldb * ldc;
      }
      working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i] =
          lda;
   }

/*move vector at*y*/
   for (i = 0; i < ssizey; i++)
      working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex + i] =
          working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i];

/*create matrix at*a*at*a + vector at*a*at*y */
   for (i = 0; i < ssizey; i++) {
      for (j = 0; j < ssizey; j++) {
         lda = 0;
         for (k = 0; k < ssizey; k++) {
            ldb = working_space[ssizex * ssizey + ssizey * i + k];
            ldc = working_space[ssizex * ssizey + ssizey * j + k];
            lda = lda + ldb * ldc;
         }
         working_space[ssizex * ssizey + ssizey * ssizey + ssizey * i + j] =
             lda;
      }
      lda = 0;
      for (k = 0; k < ssizey; k++) {
         ldb = working_space[ssizex * ssizey + ssizey * i + k];
         ldc =
             working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex +
                           k];
         lda = lda + ldb * ldc;
      }
      working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i] =
          lda;
   }

/*move at*a*at*y*/
   for (i = 0; i < ssizey; i++)
      working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex + i] =
          working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i];

/*initialization in resulting vector */
   for (i = 0; i < ssizey; i++)
      working_space[ssizex * ssizey + 2 * ssizey * ssizey + i] = 1;

        /***START OF ITERATIONS***/
   for (repet = 0; repet < numberRepetitions; repet++) {
      if (repet != 0) {
         for (i = 0; i < ssizey; i++)
            working_space[ssizex * ssizey + 2 * ssizey * ssizey + i] = TMath::Power(working_space[ssizex * ssizey + 2 * ssizey * ssizey + i], boost);
      }
      for (lindex = 0; lindex < numberIterations; lindex++) {
         for (i = 0; i < ssizey; i++) {
            lda = 0;
            for (j = 0; j < ssizey; j++) {
               ldb =
                   working_space[ssizex * ssizey + ssizey * ssizey + ssizey * i + j];
               ldc = working_space[ssizex * ssizey + 2 * ssizey * ssizey + j];
               lda = lda + ldb * ldc;
            }
            ldb =
                working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex + i];
            if (lda != 0) {
               lda = ldb / lda;
            }

            else
               lda = 0;
            ldb = working_space[ssizex * ssizey + 2 * ssizey * ssizey + i];
            lda = lda * ldb;
            working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i] = lda;
         }
         for (i = 0; i < ssizey; i++)
            working_space[ssizex * ssizey + 2 * ssizey * ssizey + i] =
                working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i];
      }
   }

/*write back resulting spectrum*/
   for (i = 0; i < ssizex; i++) {
      if (i < ssizey)
         source[i] = working_space[ssizex * ssizey + 2 * ssizey * ssizey + i];

      else
         source[i] = 0;
   }
   delete[]working_space;
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// One-dimensional high-resolution peak search function
///
/// This function searches for peaks in source spectrum. It is based on
/// deconvolution method. First the background is removed (if desired), then
/// Markov smoothed spectrum is calculated (if desired), then the response
/// function is generated according to given sigma and deconvolution is
/// carried out. The order of peaks is arranged according to their heights in
/// the spectrum after background elimination. The highest peak is the first in
/// the list. On success it returns number of found peaks.
///
/// #### Parameters:
///
///  -  source: pointer to the vector of source spectrum.
///  -  destVector: pointer to the vector of resulting deconvolved spectrum.
///  -  ssize: length of source spectrum.
///  -  sigma: sigma of searched peaks, for details we refer to manual.
///  -  threshold: threshold value in % for selected peaks, peaks with
///     amplitude less than threshold*highest_peak/100
///     are ignored, see manual.
///  -  backgroundRemove: logical variable, set if the removal of
///     background before deconvolution is desired.
///  -  deconIterations-number of iterations in deconvolution operation.
///  -  markov: logical variable, if it is true, first the source spectrum
///     is replaced by new spectrum calculated using Markov
///     chains method.
///  -  averWindow: averaging window of searched peaks, for details
///     we refer to manual (applies only for Markov method).
///
/// ### Peaks searching:
///
/// The goal of this function is to identify automatically the peaks in spectrum
/// with the presence of the continuous background and statistical
/// fluctuations - noise.
///
/// The common problems connected with correct peak identification are:
///
///  - non-sensitivity to noise, i.e., only statistically
///    relevant peaks should be identified.
///  - non-sensitivity of the algorithm to continuous
///    background.
///  - ability to identify peaks close to the edges of the
///    spectrum region. Usually peak finders fail to detect them.
///  - resolution, decomposition of Double_tts and multiplets.
///    The algorithm should be able to recognise close positioned peaks.
///  - ability to identify peaks with different sigma.
///
/// \image html TSpectrum_Searching1.jpg Fig. 27 An example of one-dimensional synthetic spectrum with found peaks denoted by markers.
///
/// #### References:
///
///  1. M.A. Mariscotti: A method for identification of peaks in the presence of
///     background and its application to spectrum analysis. NIM 50 (1967),
///     309-320.
///  2. M. Morhac;, J. Kliman, V.  Matouoek, M. Veselsky,
///     I. Turzo.:Identification of peaks in
///     multidimensional coincidence gamma-ray spectra. NIM, A443 (2000) 108-125.
///  3. Z.K. Silagadze, A new algorithm for automatic photopeak searches. NIM
///     A 376 (1996), 451.
///
/// Examples of peak searching method:
///
/// The SearchHighRes function provides users with the possibility to vary the
/// input parameters and with the access to the output deconvolved data in the
/// destination spectrum. Based on the output data one can tune the parameters.
///
/// ### Example 15 - script SearchHR1.C:
///
/// One-dimensional spectrum with found peaks denoted by markers, 3 iterations
/// steps in the deconvolution.
///
/// #### Script:
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/SearchHR1.C
/// End_Macro
///
/// ### Example 16 - script SearchHR3.C:
///
/// Influence of number of iterations (3-red, 10-blue, 100- green, 1000-magenta),
/// sigma=8, smoothing width=3.
///
/// Begin_Macro(source)
/// ../../../tutorials/spectrum/SearchHR3.C
/// End_Macro

Int_t TSpectrum::SearchHighRes(Double_t *source,Double_t *destVector, int ssize,
                                     Double_t sigma, Double_t threshold,
                                     bool backgroundRemove,int deconIterations,
                                     bool markov, int averWindow)
{
   int i, j, numberIterations = (Int_t)(7 * sigma + 0.5);
   Double_t a, b, c;
   int k, lindex, posit, imin, imax, jmin, jmax, lh_gold, priz;
   Double_t lda, ldb, ldc, area, maximum, maximum_decon;
   int xmin, xmax, l, peak_index = 0, size_ext = ssize + 2 * numberIterations, shift = numberIterations, bw = 2, w;
   Double_t maxch;
   Double_t nom, nip, nim, sp, sm, plocha = 0;
   Double_t m0low=0,m1low=0,m2low=0,l0low=0,l1low=0,detlow,av,men;
   if (sigma < 1) {
      Error("SearchHighRes", "Invalid sigma, must be greater than or equal to 1");
      return 0;
   }

   if(threshold<=0 || threshold>=100){
      Error("SearchHighRes", "Invalid threshold, must be positive and less than 100");
      return 0;
   }

   j = (Int_t) (5.0 * sigma + 0.5);
   if (j >= PEAK_WINDOW / 2) {
      Error("SearchHighRes", "Too large sigma");
      return 0;
   }

   if (markov == true) {
      if (averWindow <= 0) {
         Error("SearchHighRes", "Averaging window must be positive");
         return 0;
      }
   }

   if(backgroundRemove == true){
      if(ssize < 2 * numberIterations + 1){
         Error("SearchHighRes", "Too large clipping window");
         return 0;
      }
   }

   k = int(2 * sigma+0.5);
   if(k >= 2){
      for(i = 0;i < k;i++){
         a = i,b = source[i];
         m0low += 1,m1low += a,m2low += a * a,l0low += b,l1low += a * b;
      }
      detlow = m0low * m2low - m1low * m1low;
      if(detlow != 0)
         l1low = (-l0low * m1low + l1low * m0low) / detlow;

      else
         l1low = 0;
      if(l1low > 0)
         l1low=0;
   }

   else{
      l1low = 0;
   }

   i = (Int_t)(7 * sigma + 0.5);
   i = 2 * i;
   Double_t *working_space = new Double_t [7 * (ssize + i)];
   for (j=0;j<7 * (ssize + i);j++) working_space[j] = 0;
   for(i = 0; i < size_ext; i++){
      if(i < shift){
         a = i - shift;
         working_space[i + size_ext] = source[0] + l1low * a;
         if(working_space[i + size_ext] < 0)
            working_space[i + size_ext]=0;
      }

      else if(i >= ssize + shift){
         a = i - (ssize - 1 + shift);
         working_space[i + size_ext] = source[ssize - 1];
         if(working_space[i + size_ext] < 0)
            working_space[i + size_ext]=0;
      }

      else
         working_space[i + size_ext] = source[i - shift];
   }

   if(backgroundRemove == true){
      for(i = 1; i <= numberIterations; i++){
         for(j = i; j < size_ext - i; j++){
            if(markov == false){
               a = working_space[size_ext + j];
               b = (working_space[size_ext + j - i] + working_space[size_ext + j + i]) / 2.0;
               if(b < a)
                  a = b;

               working_space[j]=a;
            }

            else{
               a = working_space[size_ext + j];
               av = 0;
               men = 0;
               for (w = j - bw; w <= j + bw; w++){
                  if ( w >= 0 && w < size_ext){
                     av += working_space[size_ext + w];
                     men +=1;
                  }
               }
               av = av / men;
               b = 0;
               men = 0;
               for (w = j - i - bw; w <= j - i + bw; w++){
                  if ( w >= 0 && w < size_ext){
                     b += working_space[size_ext + w];
                     men +=1;
                  }
               }
               b = b / men;
               c = 0;
               men = 0;
               for (w = j + i - bw; w <= j + i + bw; w++){
                  if ( w >= 0 && w < size_ext){
                     c += working_space[size_ext + w];
                     men +=1;
                  }
               }
               c = c / men;
               b = (b + c) / 2;
               if (b < a)
                  av = b;
               working_space[j]=av;
            }
         }
         for(j = i; j < size_ext - i; j++)
            working_space[size_ext + j] = working_space[j];
      }
      for(j = 0;j < size_ext; j++){
         if(j < shift){
                  a = j - shift;
                  b = source[0] + l1low * a;
                  if (b < 0) b = 0;
            working_space[size_ext + j] = b - working_space[size_ext + j];
         }

         else if(j >= ssize + shift){
                  a = j - (ssize - 1 + shift);
                  b = source[ssize - 1];
                  if (b < 0) b = 0;
            working_space[size_ext + j] = b - working_space[size_ext + j];
         }

         else{
            working_space[size_ext + j] = source[j - shift] - working_space[size_ext + j];
         }
      }
      for(j = 0;j < size_ext; j++){
         if(working_space[size_ext + j] < 0) working_space[size_ext + j] = 0;
      }
   }

   for(i = 0; i < size_ext; i++){
      working_space[i + 6*size_ext] = working_space[i + size_ext];
   }

   if(markov == true){
      for(j = 0; j < size_ext; j++)
         working_space[2 * size_ext + j] = working_space[size_ext + j];
      xmin = 0,xmax = size_ext - 1;
      for(i = 0, maxch = 0; i < size_ext; i++){
         working_space[i] = 0;
         if(maxch < working_space[2 * size_ext + i])
            maxch = working_space[2 * size_ext + i];
         plocha += working_space[2 * size_ext + i];
      }
      if(maxch == 0) {
         delete [] working_space;
         return 0;
      }

      nom = 1;
      working_space[xmin] = 1;
      for(i = xmin; i < xmax; i++){
         nip = working_space[2 * size_ext + i] / maxch;
         nim = working_space[2 * size_ext + i + 1] / maxch;
         sp = 0,sm = 0;
         for(l = 1; l <= averWindow; l++){
            if((i + l) > xmax)
               a = working_space[2 * size_ext + xmax] / maxch;

            else
               a = working_space[2 * size_ext + i + l] / maxch;

            b = a - nip;
            if(a + nip <= 0)
               a=1;

            else
               a = TMath::Sqrt(a + nip);

            b = b / a;
            b = TMath::Exp(b);
            sp = sp + b;
            if((i - l + 1) < xmin)
               a = working_space[2 * size_ext + xmin] / maxch;

            else
               a = working_space[2 * size_ext + i - l + 1] / maxch;

            b = a - nim;
            if(a + nim <= 0)
               a = 1;

            else
               a = TMath::Sqrt(a + nim);

            b = b / a;
            b = TMath::Exp(b);
            sm = sm + b;
         }
         a = sp / sm;
         a = working_space[i + 1] = working_space[i] * a;
         nom = nom + a;
      }
      for(i = xmin; i <= xmax; i++){
         working_space[i] = working_space[i] / nom;
      }
      for(j = 0; j < size_ext; j++)
         working_space[size_ext + j] = working_space[j] * plocha;
      for(j = 0; j < size_ext; j++){
         working_space[2 * size_ext + j] = working_space[size_ext + j];
      }
      if(backgroundRemove == true){
         for(i = 1; i <= numberIterations; i++){
            for(j = i; j < size_ext - i; j++){
               a = working_space[size_ext + j];
               b = (working_space[size_ext + j - i] + working_space[size_ext + j + i]) / 2.0;
               if(b < a)
                  a = b;
               working_space[j] = a;
            }
            for(j = i; j < size_ext - i; j++)
               working_space[size_ext + j] = working_space[j];
         }
         for(j = 0; j < size_ext; j++){
            working_space[size_ext + j] = working_space[2 * size_ext + j] - working_space[size_ext + j];
         }
      }
   }
//deconvolution starts
   area = 0;
   lh_gold = -1;
   posit = 0;
   maximum = 0;
//generate response vector
   for(i = 0; i < size_ext; i++){
      lda = (Double_t)i - 3 * sigma;
      lda = lda * lda / (2 * sigma * sigma);
      j = (Int_t)(1000 * TMath::Exp(-lda));
      lda = j;
      if(lda != 0)
         lh_gold = i + 1;

      working_space[i] = lda;
      area = area + lda;
      if(lda > maximum){
         maximum = lda;
         posit = i;
      }
   }
//read source vector
   for(i = 0; i < size_ext; i++)
      working_space[2 * size_ext + i] = TMath::Abs(working_space[size_ext + i]);
//create matrix at*a(vector b)
   i = lh_gold - 1;
   if(i > size_ext)
      i = size_ext;

   imin = -i,imax = i;
   for(i = imin; i <= imax; i++){
      lda = 0;
      jmin = 0;
      if(i < 0)
         jmin = -i;
      jmax = lh_gold - 1 - i;
      if(jmax > (lh_gold - 1))
         jmax = lh_gold - 1;

      for(j = jmin;j <= jmax; j++){
         ldb = working_space[j];
         ldc = working_space[i + j];
         lda = lda + ldb * ldc;
      }
      working_space[size_ext + i - imin] = lda;
   }
//create vector p
   i = lh_gold - 1;
   imin = -i,imax = size_ext + i - 1;
   for(i = imin; i <= imax; i++){
      lda = 0;
      for(j = 0; j <= (lh_gold - 1); j++){
         ldb = working_space[j];
         k = i + j;
         if(k >= 0 && k < size_ext){
            ldc = working_space[2 * size_ext + k];
            lda = lda + ldb * ldc;
         }

      }
      working_space[4 * size_ext + i - imin] = lda;
   }
//move vector p
   for(i = imin; i <= imax; i++)
      working_space[2 * size_ext + i - imin] = working_space[4 * size_ext + i - imin];
//initialization of resulting vector
   for(i = 0; i < size_ext; i++)
      working_space[i] = 1;
//START OF ITERATIONS
   for(lindex = 0; lindex < deconIterations; lindex++){
      for(i = 0; i < size_ext; i++){
         if(TMath::Abs(working_space[2 * size_ext + i]) > 0.00001 && TMath::Abs(working_space[i]) > 0.00001){
            lda=0;
            jmin = lh_gold - 1;
            if(jmin > i)
               jmin = i;

            jmin = -jmin;
            jmax = lh_gold - 1;
            if(jmax > (size_ext - 1 - i))
               jmax=size_ext-1-i;

            for(j = jmin; j <= jmax; j++){
               ldb = working_space[j + lh_gold - 1 + size_ext];
               ldc = working_space[i + j];
               lda = lda + ldb * ldc;
            }
            ldb = working_space[2 * size_ext + i];
            if(lda != 0)
               lda = ldb / lda;

            else
               lda = 0;

            ldb = working_space[i];
            lda = lda * ldb;
            working_space[3 * size_ext + i] = lda;
         }
      }
      for(i = 0; i < size_ext; i++){
         working_space[i] = working_space[3 * size_ext + i];
      }
   }
//shift resulting spectrum
   for(i=0;i<size_ext;i++){
      lda = working_space[i];
      j = i + posit;
      j = j % size_ext;
      working_space[size_ext + j] = lda;
   }
//write back resulting spectrum
   maximum = 0, maximum_decon = 0;
   j = lh_gold - 1;
   for(i = 0; i < size_ext - j; i++){
      if(i >= shift && i < ssize + shift){
         working_space[i] = area * working_space[size_ext + i + j];
         if(maximum_decon < working_space[i])
            maximum_decon = working_space[i];
         if(maximum < working_space[6 * size_ext + i])
            maximum = working_space[6 * size_ext + i];
      }

      else
         working_space[i] = 0;
   }
   lda=1;
   if(lda>threshold)
      lda=threshold;
   lda=lda/100;

//searching for peaks in deconvolved spectrum
   for(i = 1; i < size_ext - 1; i++){
      if(working_space[i] > working_space[i - 1] && working_space[i] > working_space[i + 1]){
         if(i >= shift && i < ssize + shift){
            if(working_space[i] > lda*maximum_decon && working_space[6 * size_ext + i] > threshold * maximum / 100.0){
               for(j = i - 1, a = 0, b = 0; j <= i + 1; j++){
                  a += (Double_t)(j - shift) * working_space[j];
                  b += working_space[j];
               }
               a = a / b;
               if(a < 0)
                  a = 0;

               if(a >= ssize)
                  a = ssize - 1;
               if(peak_index == 0){
                  fPositionX[0] = a;
                  peak_index = 1;
               }

               else{
                  for(j = 0, priz = 0; j < peak_index && priz == 0; j++){
                     if(working_space[6 * size_ext + shift + (Int_t)a] > working_space[6 * size_ext + shift + (Int_t)fPositionX[j]])
                        priz = 1;
                  }
                  if(priz == 0){
                     if(j < fMaxPeaks){
                        fPositionX[j] = a;
                     }
                  }

                  else{
                     for(k = peak_index; k >= j; k--){
                        if(k < fMaxPeaks){
                           fPositionX[k] = fPositionX[k - 1];
                        }
                     }
                     fPositionX[j - 1] = a;
                  }
                  if(peak_index < fMaxPeaks)
                     peak_index += 1;
               }
            }
         }
      }
   }

   for(i = 0; i < ssize; i++) destVector[i] = working_space[i + shift];
   delete [] working_space;
   fNPeaks = peak_index;
   if(peak_index == fMaxPeaks)
      Warning("SearchHighRes", "Peak buffer full");
   return fNPeaks;
}

////////////////////////////////////////////////////////////////////////////////
/// Old name of SearcHighRes introduced for back compatibility.
/// This function will be removed after the June 2006 release

Int_t TSpectrum::Search1HighRes(Double_t *source,Double_t *destVector, int ssize,
                                     Double_t sigma, Double_t threshold,
                                     bool backgroundRemove,int deconIterations,
                                     bool markov, int averWindow)
{


   return SearchHighRes(source,destVector,ssize,sigma,threshold,backgroundRemove,
                        deconIterations,markov,averWindow);
}

////////////////////////////////////////////////////////////////////////////////
/// Static function, interface to TSpectrum::Search.

Int_t TSpectrum::StaticSearch(const TH1 *hist, Double_t sigma, Option_t *option, Double_t threshold)
{

   TSpectrum s;
   return s.Search(hist,sigma,option,threshold);
}

////////////////////////////////////////////////////////////////////////////////
/// Static function, interface to TSpectrum::Background.

TH1 *TSpectrum::StaticBackground(const TH1 *hist,Int_t niter, Option_t *option)
{
   TSpectrum s;
   return s.Background(hist,niter,option);
}
