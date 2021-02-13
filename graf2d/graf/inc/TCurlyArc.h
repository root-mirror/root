// @(#)root/graf:$Id$
// Author: Otto Schaile   20/11/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TCurlyArc
#define ROOT_TCurlyArc

#include "TCurlyLine.h"

class TPoint;

class TCurlyArc : public TCurlyLine {

private:
   Double_t fR1;                  ///< Radius of arc
   Double_t fPhimin;              ///< start phi (degrees)
   Double_t fPhimax;              ///< end phi (degrees)
   Double_t fTheta;               ///< used internally

   static Double_t fgDefaultWaveLength;   ///< default wavelength
   static Double_t fgDefaultAmplitude;    ///< default amplitude
   static Bool_t   fgDefaultIsCurly;      ///< default curly type

public:
   TCurlyArc();
   TCurlyArc(Double_t x1, Double_t y1, Double_t rad,
             Double_t phimin, Double_t phimax,
             Double_t wl = .02, Double_t amp = .01);
   virtual     ~TCurlyArc(){;}
   virtual void Build();
   Int_t        DistancetoPrimitive(Int_t px, Int_t py);
   void         ExecuteEvent(Int_t event, Int_t px, Int_t py);
   Double_t     GetRadius() const {return fR1;}
   Double_t     GetPhimin() const {return fPhimin;}
   Double_t     GetPhimax() const {return fPhimax;}
   virtual void SetCenter(Double_t x1, Double_t y1); // *MENU* *ARGS={x1=>fX1,y1=>fY1}
   virtual void SetRadius(Double_t radius);          // *MENU* *ARGS={radius=>fR1}
   virtual void SetPhimin(Double_t phimin);          // *MENU* *ARGS={phimin=>fPhimin}
   virtual void SetPhimax(Double_t phimax);          // *MENU* *ARGS={phimax=>fPhimax}
   virtual void SavePrimitive(std::ostream &out, Option_t * = "");

   static void          SetDefaultWaveLength(Double_t WaveLength);
   static void          SetDefaultAmplitude (Double_t Amplitude );
   static void          SetDefaultIsCurly   (Bool_t   IsCurly   );
   static Double_t      GetDefaultWaveLength();
   static Double_t      GetDefaultAmplitude ();
   static Bool_t        GetDefaultIsCurly   ();
   virtual Rectangle_t  GetBBox();
   virtual TPoint       GetBBoxCenter();
   virtual void         SetBBoxCenter(const TPoint &p);
   virtual void         SetBBoxCenterX(const Int_t x);
   virtual void         SetBBoxCenterY(const Int_t y);
   virtual void         SetBBoxX1(const Int_t x);
   virtual void         SetBBoxX2(const Int_t x);
   virtual void         SetBBoxY1(const Int_t y);
   virtual void         SetBBoxY2(const Int_t y);

   ClassDef(TCurlyArc,3) // A curly arc
};

#endif

