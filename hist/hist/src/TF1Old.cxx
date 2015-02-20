// @(#)root/hist:$Id$
// Author: Rene Brun   18/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TF1Old.h"

#include "TH1.h"

ClassImp(TF1Old)

//______________________________________________________________________________
TF1Old::TF1Old(): TFormulaOld(), TAttLine(), TAttFill(), TAttMarker()
{
   // F1 default constructor.

   fXmin      = 0;
   fXmax      = 0;
   fNpx       = 100;
   fType      = 0;
   fNpfits    = 0;
   fNDF       = 0;
   fNsave     = 0;
   fChisquare = 0;
   fParErrors = 0;
   fParMin    = 0;
   fParMax    = 0;
   fSave      = 0;
   fMinimum   = -1111;
   fMaximum   = -1111;
   SetFillStyle(0);
}

//______________________________________________________________________________
TF1Old::~TF1Old()
{
   // TF1 default destructor.

   if (fParMin)    delete [] fParMin;
   if (fParMax)    delete [] fParMax;
   if (fParErrors) delete [] fParErrors;
   if (fSave)      delete [] fSave;
}

//______________________________________________________________________________
void TF1Old::Streamer(TBuffer &b)
{
   // Stream a class object.

   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t v = b.ReadVersion(&R__s, &R__c);
      if (v > 4) {
         b.ReadClassBuffer(TF1Old::Class(), this, v, R__s, R__c);
         if (v == 5 && fNsave > 0) {
            //correct badly saved fSave in 3.00/06
            Int_t np = fNsave - 3;
            fSave[np]   = fSave[np-1];
            fSave[np+1] = fXmin;
            fSave[np+2] = fXmax;
         }
         return;
      }
      //====process old versions before automatic schema evolution
      TFormulaOld::Streamer(b);
      TAttLine::Streamer(b);
      TAttFill::Streamer(b);
      TAttMarker::Streamer(b);
      if (v < 4) {
         Float_t xmin,xmax;
         b >> xmin; fXmin = xmin;
         b >> xmax; fXmax = xmax;
      } else {
         b >> fXmin;
         b >> fXmax;
      }
      b >> fNpx;
      b >> fType;
      b >> fChisquare;
      b.ReadArray(fParErrors);
      if (v > 1) {
         b.ReadArray(fParMin);
         b.ReadArray(fParMax);
      } else {
         fParMin = new Double_t[fNpar+1];
         fParMax = new Double_t[fNpar+1];
      }
      b >> fNpfits;
      if (v == 1) {
         TH1 * histogram;
         b >> histogram;
         delete histogram; //fHistogram = 0;
      }
      if (v > 1) {
         if (v < 4) {
            Float_t minimum,maximum;
            b >> minimum; fMinimum =minimum;
            b >> maximum; fMaximum =maximum;
         } else {
            b >> fMinimum;
            b >> fMaximum;
         }
      }
      if (v > 2) {
         b >> fNsave;
         if (fNsave > 0) {
            fSave = new Double_t[fNsave+10];
            b.ReadArray(fSave);
            //correct fSave limits to match new version
            fSave[fNsave]   = fSave[fNsave-1];
            fSave[fNsave+1] = fSave[fNsave+2];
            fSave[fNsave+2] = fSave[fNsave+3];
            fNsave += 3;
         } else fSave = 0;
      }
      b.CheckByteCount(R__s, R__c, TF1Old::IsA());
      //====end of old versions

   } else {
      // this will be needed if we want to write in old format
      //Int_t saved = 0;
      // if (fType > 0 && fNsave <= 0) { saved = 1; Save(fXmin,fXmax,0,0,0,0);}

      b.WriteClassBuffer(TF1Old::Class(),this);

      //if (saved) {delete [] fSave; fSave = 0; fNsave = 0;}
   }
}
