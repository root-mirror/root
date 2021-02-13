// @(#)root/base:$Id$
// Author: Anna Kreshuk  10/04/2006

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TVirtualFFT
\ingroup Base

TVirtualFFT is an interface class for Fast Fourier Transforms.

The default FFT library is FFTW. To use it, FFTW3 library should already
be installed, and ROOT should be have fftw3 module enabled, with the directories
of fftw3 include file and library specified (see installation instructions).
Function SetDefaultFFT() allows to change the default library.

## Available transform types:
FFT:
  - "C2CFORWARD" - a complex input/output discrete Fourier transform (DFT)
                   in one or more dimensions, -1 in the exponent
  - "C2CBACKWARD"- a complex input/output discrete Fourier transform (DFT)
                   in one or more dimensions, +1 in the exponent
  - "R2C"        - a real-input/complex-output discrete Fourier transform (DFT)
                   in one or more dimensions,
  - "C2R"        - inverse transforms to "R2C", taking complex input
                   (storing the non-redundant half of a logically Hermitian array)
                   to real output
  - "R2HC"       - a real-input DFT with output in ¡Èhalfcomplex¡É format,
                   i.e. real and imaginary parts for a transform of size n stored as
                   r0, r1, r2, ..., rn/2, i(n+1)/2-1, ..., i2, i1
  - "HC2R"       - computes the reverse of FFTW_R2HC, above
  - "DHT"        - computes a discrete Hartley transform

## Sine/cosine transforms:
Different types of transforms are specified by parameter kind of the SineCosine() static
function. 4 different kinds of sine and cosine transforms are available

  - DCT-I  (REDFT00 in FFTW3 notation)- kind=0
  - DCT-II (REDFT01 in FFTW3 notation)- kind=1
  - DCT-III(REDFT10 in FFTW3 notation)- kind=2
  - DCT-IV (REDFT11 in FFTW3 notation)- kind=3
  - DST-I  (RODFT00 in FFTW3 notation)- kind=4
  - DST-II (RODFT01 in FFTW3 notation)- kind=5
  - DST-III(RODFT10 in FFTW3 notation)- kind=6
  - DST-IV (RODFT11 in FFTW3 notation)- kind=7

Formulas and detailed descriptions can be found in the chapter
"What FFTW really computes" of the FFTW manual

NOTE: FFTW computes unnormalized transforms, so doing a transform, followed by its
      inverse will give the original array, multiplied by normalization constant
      (transform size(N) for FFT, 2*(N-1) for DCT-I, 2*(N+1) for DST-I, 2*N for
      other sine/cosine transforms)

## How to use it:
Call to the static function FFT returns a pointer to a fast Fourier transform
with requested parameters. Call to the static function SineCosine returns a
pointer to a sine or cosine transform with requested parameters. Example:
~~~ {.cpp}
{
   Int_t N = 10; Double_t *in = new Double_t[N];
   TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &N, "R2C");
   fftr2c->SetPoints(in);
   fftr2c->Transform();
   Double_t re, im;
   for (Int_t i=0; i<N; i++)
      fftr2c->GetPointComplex(i, re, im);
   ...
   fftr2c->SetPoints(in2);
   ...
   fftr2c->SetPoints(in3);
   ...
}
~~~
Different options are explained in the function comments
*/

#include "TROOT.h"
#include "TVirtualFFT.h"
#include "TPluginManager.h"
#include "TError.h"

TVirtualFFT *TVirtualFFT::fgFFT    = 0;
TString      TVirtualFFT::fgDefault   = "";

ClassImp(TVirtualFFT);

////////////////////////////////////////////////////////////////////////////////
///destructor

TVirtualFFT::~TVirtualFFT()
{
   if (this==fgFFT)
      fgFFT = 0;
}

////////////////////////////////////////////////////////////////////////////////
///Returns a pointer to the FFT of requested size and type.
///
/// \param[in] ndim      number of transform dimensions
/// \param[in] n         sizes of each dimension (an array at least ndim long)
/// \param [in] option   consists of 3 parts - flag option and an option to create a new TVirtualFFT
///                      1. transform type option:
///                         Available transform types are:
///                         C2CForward, C2CBackward, C2R, R2C, R2HC, HC2R, DHT
///                         see class description for details
///                      2. flag option: choosing how much time should be spent in planning the transform:
///                         Possible options:
///                         - "ES" (from "estimate") - no time in preparing the transform,
///                            but probably sub-optimal  performance
///                         - "M"  (from "measure")  - some time spend in finding the optimal way
///                            to do the transform
///                         - "P" (from "patient")   - more time spend in finding the optimal way
///                            to do the transform
///                         - "EX" (from "exhaustive") - the most optimal way is found
///                            This option should be chosen depending on how many transforms of the
///                            same size and type are going to be done.
///                         Planning is only done once, for the first transform of this size and type.
///                      3. option allowing to choose between the global fgFFT and a new TVirtualFFT object
///                         ""  - default, changes and returns the global fgFFT variable
///                         "K" (from "keep")- without touching the global fgFFT,
///                         creates and returns a new TVirtualFFT*. User is then responsible for deleting it.
///
/// Examples of valid options: "R2C ES K", "C2CF M", "DHT P K", etc.

TVirtualFFT* TVirtualFFT::FFT(Int_t ndim, Int_t *n, Option_t *option)
{

   Int_t inputtype=0, currenttype=0;
   TString opt = option;
   opt.ToUpper();
   //find the tranform flag
   Option_t *flag;
   flag = "ES";
   if (opt.Contains("ES")) flag = "ES";
   if (opt.Contains("M"))  flag = "M";
   if (opt.Contains("P"))  flag = "P";
   if (opt.Contains("EX")) flag = "EX";

   Int_t ndiff = 0;

   if (!opt.Contains("K")) {
      if (fgFFT){
         //if the global transform exists, check if it should be changed
         if (fgFFT->GetNdim()!=ndim)
            ndiff++;
         else {
            Int_t *ncurrent = fgFFT->GetN();
            for (Int_t i=0; i<ndim; i++){
               if (n[i]!=ncurrent[i])
                  ndiff++;
            }
         }
         Option_t *t = fgFFT->GetType();
         if (!opt.Contains(t)) {
            if (opt.Contains("HC") || opt.Contains("DHT"))
               inputtype = 1;
            if (strcmp(t,"R2HC")==0 || strcmp(t,"HC2R")==0 || strcmp(t,"DHT")==0)
               currenttype=1;

            if (!(inputtype==1 && currenttype==1))
               ndiff++;
         }
         if (ndiff>0){
            delete fgFFT;
            fgFFT = 0;
         }
      }
   }

   Int_t sign = 0;
   if (opt.Contains("C2CB") || opt.Contains("C2R"))
      sign = 1;
   if (opt.Contains("C2CF") || opt.Contains("R2C"))
      sign = -1;

   TVirtualFFT *fft = 0;
   if (opt.Contains("K") || !fgFFT) {

      R__LOCKGUARD(gROOTMutex);

      TPluginHandler *h;
      TString pluginname;
      if (fgDefault.Length()==0) fgDefault="fftw";
      if (strcmp(fgDefault.Data(),"fftw")==0) {
         if (opt.Contains("C2C")) pluginname = "fftwc2c";
         if (opt.Contains("C2R")) pluginname = "fftwc2r";
         if (opt.Contains("R2C")) pluginname = "fftwr2c";
         if (opt.Contains("HC") || opt.Contains("DHT")) pluginname = "fftwr2r";
         if ((h=gROOT->GetPluginManager()->FindHandler("TVirtualFFT", pluginname))) {
            if (h->LoadPlugin()==-1) {
               ::Error("TVirtualFFT::FFT", "handler not found");
               return 0;
            }
            fft = (TVirtualFFT*)h->ExecPlugin(3, ndim, n, kFALSE);
            if (!fft) {
               ::Error("TVirtualFFT::FFT", "plugin failed to create TVirtualFFT object");
               return 0;
            }
            Int_t *kind = new Int_t[1];
            if (pluginname=="fftwr2r") {
               if (opt.Contains("R2HC")) kind[0] = 10;
               if (opt.Contains("HC2R")) kind[0] = 11;
               if (opt.Contains("DHT")) kind[0] = 12;
            }
            fft->Init(flag, sign, kind);
            if (!opt.Contains("K")) {
               fgFFT = fft;
            }
            delete [] kind;
            return fft;
         }
         else {
            ::Error("TVirtualFFT::FFT", "plugin not found");
            return 0;
         }
      }
   } else {

      R__LOCKGUARD(gROOTMutex);

      //if the global transform already exists and just needs to be reinitialised
      //with different parameters
      if (fgFFT->GetSign()!=sign || !opt.Contains(fgFFT->GetTransformFlag()) || !opt.Contains(fgFFT->GetType())) {
         Int_t *kind = new Int_t[1];
         if (inputtype==1) {
            if (opt.Contains("R2HC")) kind[0] = 10;
            if (opt.Contains("HC2R")) kind[0] = 11;
            if (opt.Contains("DHT")) kind[0] = 12;
         }
         fgFFT->Init(flag, sign, kind);
         delete [] kind;
      }
   }
   return fgFFT;
}

////////////////////////////////////////////////////////////////////////////////
///Returns a pointer to a sine or cosine transform of requested size and kind
///
///Parameters:
/// \param [in] ndim      number of transform dimensions
/// \param [in] n         sizes of each dimension (an array at least ndim long)
/// \param [in] r2rkind   transform kind for each dimension
///                       4 different kinds of sine and cosine transforms are available
///                       - DCT-I    - kind=0
///                       - DCT-II   - kind=1
///                       - DCT-III  - kind=2
///                       - DCT-IV   - kind=3
///                       - DST-I    - kind=4
///                       - DST-II   - kind=5
///                       - DST-III  - kind=6
///                       - DST-IV   - kind=7
/// \param [in] option : consists of 2 parts
///         - flag option and an option to create a new TVirtualFFT
///         - flag option: choosing how much time should be spent in planning the transform:
///           Possible options:
///           - "ES" (from "estimate") - no time in preparing the transform,
///                                    but probably sub-optimal  performance
///           - "M"  (from "measure")  - some time spend in finding the optimal way
///                                    to do the transform
///           - "P" (from "patient")   - more time spend in finding the optimal way
///                                    to do the transform
///           - "EX" (from "exhaustive") - the most optimal way is found
///              This option should be chosen depending on how many transforms of the
///              same size and type are going to be done.
///              Planning is only done once, for the first transform of this size and type.
///              - option allowing to choose between the global fgFFT and a new TVirtualFFT object
///           - ""  - default, changes and returns the global fgFFT variable
///           - "K" (from "keep")- without touching the global fgFFT,
///             creates and returns a new TVirtualFFT*. User is then responsible for deleting it.
/// Examples of valid options: "ES K", "EX", etc

TVirtualFFT* TVirtualFFT::SineCosine(Int_t ndim, Int_t *n, Int_t *r2rkind, Option_t *option)
{
   TString opt = option;
   //find the tranform flag
   Option_t *flag;
   flag = "ES";
   if (opt.Contains("ES")) flag = "ES";
   if (opt.Contains("M"))  flag = "M";
   if (opt.Contains("P"))  flag = "P";
   if (opt.Contains("EX")) flag = "EX";

   if (!opt.Contains("K")) {
      if (fgFFT){
         Int_t ndiff = 0;
         if (fgFFT->GetNdim()!=ndim || strcmp(fgFFT->GetType(),"R2R")!=0)
            ndiff++;
         else {
            Int_t *ncurrent = fgFFT->GetN();
            for (Int_t i=0; i<ndim; i++) {
               if (n[i] != ncurrent[i])
                  ndiff++;
            }

         }
         if (ndiff>0) {
            delete fgFFT;
            fgFFT = 0;
         }
      }
   }
   TVirtualFFT *fft = 0;

   R__LOCKGUARD(gROOTMutex);

   if (!fgFFT || opt.Contains("K")) {
      TPluginHandler *h;
      TString pluginname;
      if (fgDefault.Length()==0) fgDefault="fftw";
      if (strcmp(fgDefault.Data(),"fftw")==0) {
         pluginname = "fftwr2r";
         if ((h=gROOT->GetPluginManager()->FindHandler("TVirtualFFT", pluginname))) {
            if (h->LoadPlugin()==-1){
               ::Error("TVirtualFFT::SineCosine", "handler not found");
               return 0;
            }
            fft = (TVirtualFFT*)h->ExecPlugin(3, ndim, n, kFALSE);
            if (!fft) {
               ::Error("TVirtualFFT::SineCosine", "plugin failed to create TVirtualFFT object");
               return 0;
            }
            fft->Init(flag, 0, r2rkind);
            if (!opt.Contains("K"))
               fgFFT = fft;
            return fft;
         } else {
            ::Error("TVirtualFFT::SineCosine", "handler not found");
            return 0;
         }
      }
   }

   //if (fgFFT->GetTransformFlag()!=flag)
   fgFFT->Init(flag,0, r2rkind);
   return fgFFT;
}

////////////////////////////////////////////////////////////////////////////////
/// static: return current fgFFT

TVirtualFFT* TVirtualFFT::GetCurrentTransform()
{
   if (fgFFT)
      return fgFFT;
   else{
      ::Warning("TVirtualFFT::GetCurrentTransform", "fgFFT is not defined yet");
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// static: set the current transfrom to parameter

void TVirtualFFT::SetTransform(TVirtualFFT* fft)
{
   fgFFT = fft;
}

////////////////////////////////////////////////////////////////////////////////
/// static: return the name of the default fft

const char *TVirtualFFT::GetDefaultFFT()
{
   return fgDefault.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// static: set name of default fft

void TVirtualFFT::SetDefaultFFT(const char *name)
{
   if (fgDefault == name) return;
   delete fgFFT;
   fgFFT = 0;
   fgDefault = name;
}
