// @(#)root/graf:$Id$
// Author: Otto Schaile   20/11/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TCurlyLine
#define ROOT_TCurlyLine

#include "TPolyLine.h"
#include "TAttBBox2D.h"

class TPoint;

class TCurlyLine : public TPolyLine, public TAttBBox2D {

protected:
   Double_t fX1;             ///< start x, center for arc
   Double_t fY1;             ///< start y, center for arc
   Double_t fX2;             ///< end x
   Double_t fY2;             ///< end y
   Double_t fWaveLength;     ///< wavelength of sinusoid in percent of pad height
   Double_t fAmplitude;      ///< amplitude of sinusoid in percent of pad height
   Int_t    fNsteps;         ///< used internally (controls precision)
   Bool_t   fIsCurly;        ///< true: Gluon, false: Gamma

   static Double_t fgDefaultWaveLength;   ///< default wavelength
   static Double_t fgDefaultAmplitude;    ///< default amplitude
   static Bool_t   fgDefaultIsCurly;      ///< default curly type

public:
   TCurlyLine();
   TCurlyLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2,
              Double_t wl = .02,
              Double_t amp = .01);
   virtual ~TCurlyLine(){;}
   virtual void Build();
   Int_t        DistancetoPrimitive(Int_t px, Int_t py);
   void         ExecuteEvent(Int_t event, Int_t px, Int_t py);
   Bool_t       GetCurly() const     {return fIsCurly;}
   Double_t     GetWaveLength() const{return fWaveLength;}
   Double_t     GetAmplitude() const {return fAmplitude;}
   Double_t     GetStartX() const    {return fX1;}
   Double_t     GetEndX() const      {return fX2;}
   Double_t     GetStartY() const    {return fY1;}
   Double_t     GetEndY() const      {return fY2;}
   virtual void SetCurly();                             // *MENU*
   virtual void SetWavy();                              // *MENU*
   virtual void SetWaveLength(Double_t WaveLength);     // *MENU* *ARGS={WaveLength=>fWaveLength}
   virtual void SetAmplitude(Double_t x);               // *MENU* *ARGS={x=>fAmplitude}
   virtual void SetStartPoint(Double_t x1, Double_t y1);
   virtual void SetEndPoint  (Double_t x2, Double_t y2);
   virtual void SavePrimitive(std::ostream &out, Option_t * = "");

   static void     SetDefaultWaveLength(Double_t WaveLength);
   static void     SetDefaultAmplitude (Double_t Amplitude );
   static void     SetDefaultIsCurly   (Bool_t   IsCurly   );
   static Double_t GetDefaultWaveLength();
   static Double_t GetDefaultAmplitude ();
   static Bool_t   GetDefaultIsCurly   ();
   virtual Rectangle_t  GetBBox();
   virtual TPoint       GetBBoxCenter();
   virtual void         SetBBoxCenter(const TPoint &p);
   virtual void         SetBBoxCenterX(const Int_t x);
   virtual void         SetBBoxCenterY(const Int_t y);
   virtual void         SetBBoxX1(const Int_t x);
   virtual void         SetBBoxX2(const Int_t x);
   virtual void         SetBBoxY1(const Int_t y);
   virtual void         SetBBoxY2(const Int_t y);

   ClassDef(TCurlyLine,3) // A curly polyline
};

#endif
