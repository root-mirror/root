// @(#)root/geom:$Id$
// Author: Andrei Gheata   09/02/03

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoOverlap
#define ROOT_TGeoOverlap

#include "TNamed.h"

#include "TAttLine.h"

#include "TAttFill.h"

#include "TAtt3D.h"

#include "TGeoMatrix.h"

///////////////////////////////////////////////////////////////////////////
//                                                                       //
// TGeoOverlap - base class describing geometry overlaps. Overlaps apply //
//   to the nodes contained inside a volume. These should not overlap to //
//   each other nor extrude the shape of their mother volume.            //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

class TGeoVolume;
class TPolyMarker3D;
class TBrowser;

class TGeoOverlap : public TNamed,
                    public TAttLine,
                    public TAttFill,
                    public TAtt3D
{
public:
enum EOverlapType {
   kGeoOverlap    = BIT(14),
   kGeoExtrusion  = BIT(15)
};

private:
   TGeoOverlap(const TGeoOverlap&) = delete;
   TGeoOverlap& operator=(const TGeoOverlap&) = delete;

protected:
   Double_t         fOverlap;     // overlap distance
   TGeoVolume      *fVolume1;     // first volume
   TGeoVolume      *fVolume2;     // second volume
   TGeoHMatrix     *fMatrix1;     // positioning matrix for first volume
   TGeoHMatrix     *fMatrix2;     // positioning matrix for second volume
   TPolyMarker3D   *fMarker;      // points in the overlapping region

public:
   TGeoOverlap();
   TGeoOverlap(const char *name, TGeoVolume *vol1, TGeoVolume *vol2,
               const TGeoMatrix *matrix1, const TGeoMatrix *matrix2,
               Bool_t isovlp=kTRUE,  Double_t ovlp=0.01);
   virtual           ~TGeoOverlap();

   void              Browse(TBrowser *b);
   virtual Int_t     Compare(const TObject *obj) const;
   virtual Int_t     DistancetoPrimitive(Int_t px, Int_t py);
   virtual void      Draw(Option_t *option=""); // *MENU*
   virtual void      ExecuteEvent(Int_t event, Int_t px, Int_t py);
   TPolyMarker3D    *GetPolyMarker() const {return fMarker;}
   TGeoVolume       *GetFirstVolume() const {return fVolume1;}
   TGeoVolume       *GetSecondVolume() const {return fVolume2;}
   TGeoHMatrix      *GetFirstMatrix() const {return fMatrix1;}
   TGeoHMatrix      *GetSecondMatrix() const {return fMatrix2;}
   Double_t          GetOverlap() const {return fOverlap;}
   Bool_t            IsExtrusion() const {return TObject::TestBit(kGeoExtrusion);}
   Bool_t            IsOverlap() const {return TObject::TestBit(kGeoOverlap);}
   Bool_t            IsFolder() const {return kFALSE;}
   virtual Bool_t    IsSortable() const {return kTRUE;}
   virtual void      Paint(Option_t *option="");
   virtual void      Print(Option_t *option="") const; // *MENU*
   virtual void      PrintInfo() const;
   virtual void      Sizeof3D() const;
   void              SampleOverlap(Int_t npoints=1000000); // *MENU*
   void              SetIsExtrusion(Bool_t flag=kTRUE) {TObject::SetBit(kGeoExtrusion,flag); TObject::SetBit(kGeoOverlap,!flag);}
   void              SetIsOverlap(Bool_t flag=kTRUE) {TObject::SetBit(kGeoOverlap,flag); TObject::SetBit(kGeoExtrusion,!flag);}
   void              SetNextPoint(Double_t x, Double_t y, Double_t z);
   void              SetFirstVolume(TGeoVolume *vol) {fVolume1=vol;}
   void              SetSecondVolume(TGeoVolume *vol) {fVolume2=vol;}
   void              SetFirstMatrix(TGeoMatrix *matrix) {*fMatrix1 = matrix;}
   void              SetSecondMatrix(TGeoMatrix *matrix) {*fMatrix2 = matrix;}
   void              SetOverlap(Double_t ovlp)  {fOverlap=ovlp;}
   void              Validate() const; // *MENU*

   ClassDef(TGeoOverlap, 2)         // base class for geometical overlaps
};

#endif

