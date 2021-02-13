// @(#)root/graf:$Id$
// Author: Guido Volpi, Olivier Couet  03/11/2006

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPie
#define ROOT_TPie

#include <TNamed.h>
#include <TAttText.h>
#include <TString.h>

class TH1;
class TPieSlice;
class TLegend;

class TPie : public TNamed , public TAttText {
private:
   void Init(Int_t np, Double_t ao, Double_t x, Double_t y, Double_t r);
   void DrawGhost();

   Float_t  fSum;             ///<!Sum for the slice values
   Float_t *fSlices;          ///<!Subdivisions of the slices
   TLegend *fLegend;          ///<!Legend for this piechart

protected:
   Double_t    fX;              ///< X coordinate of the pie centre
   Double_t    fY;              ///< Y coordinate of the pie centre
   Double_t    fRadius;         ///< Radius Pie radius
   Double_t    fAngularOffset;  ///< Offset angular offset for the first slice
   Float_t     fLabelsOffset;   ///< LabelsOffset offset of label
   TString     fLabelFormat;    ///< Format format of the slices' label
   TString     fValueFormat;    ///< Vform numeric format for the value
   TString     fFractionFormat; ///< Rform numeric format for the fraction of a slice
   TString     fPercentFormat;  ///< Pfrom numeric format for the percent of a slice
   Int_t       fNvals;          ///< Number of elements
   TPieSlice **fPieSlices;      ///<[fNvals] Slice array of this pie-chart
   Bool_t      fIs3D;           ///<! true if the pseudo-3d is enabled
   Double_t    fHeight;         ///< Height of the slice in pixel
   Float_t     fAngle3D;        ///< The angle of the pseudo-3d view

public:
   TPie();
   TPie(const char *,const char *, Int_t);
   TPie(const char *,const char *, Int_t, Double_t *,Int_t *cols=0, const char *lbls[]=0);
   TPie(const char *,const char *, Int_t, Float_t *,Int_t *cols=0, const char *lbls[]=0);
   TPie(const TH1 *h);
   TPie(const TPie&);
   ~TPie();

   virtual Int_t  DistancetoPrimitive(Int_t px, Int_t py);
   Int_t          DistancetoSlice(Int_t,Int_t);
   virtual void   Draw(Option_t *option="l"); // *MENU*
   virtual void   ExecuteEvent(Int_t,Int_t,Int_t);
   Float_t        GetAngle3D() { return fAngle3D; }
   Double_t       GetAngularOffset() { return fAngularOffset; }
   Int_t          GetEntryFillColor(Int_t);
   Int_t          GetEntryFillStyle(Int_t);
   const char*    GetEntryLabel(Int_t);
   Int_t          GetEntryLineColor(Int_t);
   Int_t          GetEntryLineStyle(Int_t);
   Int_t          GetEntryLineWidth(Int_t);
   Double_t       GetEntryRadiusOffset(Int_t);
   Double_t       GetEntryVal(Int_t);
   const char    *GetFractionFormat() { return fFractionFormat.Data(); }
   Double_t       GetHeight() { return fHeight; }
   const char    *GetLabelFormat() { return fLabelFormat.Data(); }
   Float_t        GetLabelsOffset() { return fLabelsOffset; }
   TLegend       *GetLegend();
   Int_t          GetEntries() { return fNvals; }
   const char    *GetPercentFormat() { return fPercentFormat.Data(); }
   Double_t       GetRadius() { return fRadius;}
   TPieSlice     *GetSlice(Int_t i);
   const char    *GetValueFormat() { return fValueFormat.Data(); }
   Double_t       GetX() { return fX; }
   Double_t       GetY() { return fY; }
   TLegend       *MakeLegend(Double_t x1=.65,Double_t y1=.65,Double_t x2=.95, Double_t y2=.95, const char *leg_header="");
   void           MakeSlices(Bool_t force=kFALSE);
   virtual void   Paint(Option_t *);
   void           SavePrimitive(std::ostream &out, Option_t *opts="");
   void           SetAngle3D(Float_t val = 30.); // *MENU*
   void           SetAngularOffset(Double_t);
   void           SetCircle(Double_t x=.5, Double_t y=.5, Double_t rad=.4);
   void           SetEntryLabel(Int_t, const char *text="Slice");
   void           SetEntryLineColor(Int_t, Int_t);
   void           SetEntryLineStyle(Int_t, Int_t);
   void           SetEntryLineWidth(Int_t, Int_t);
   void           SetEntryFillColor(Int_t, Int_t);
   void           SetEntryFillStyle(Int_t, Int_t);
   void           SetEntryRadiusOffset(Int_t, Double_t);
   void           SetEntryVal(Int_t, Double_t);
   void           SetFillColors(Int_t*);
   void           SetFractionFormat(const char*); // *MENU*
   void           SetHeight(Double_t val=.08); // *MENU*
   void           SetLabelFormat(const char *); // *MENU*
   void           SetLabels(const char *[]);
   void           SetLabelsOffset(Float_t); // *MENU*
   void           SetPercentFormat(const char *); // *MENU*
   void           SetRadius(Double_t); // *MENU*
   void           SetValueFormat(const char *); // *MENU*
   void           SetX(Double_t); // *MENU*
   void           SetY(Double_t); // *MENU*
   void           SortSlices(Bool_t amode=kTRUE,Float_t merge_thresold=.0);

   ClassDef(TPie,1) //Pie chart graphics class
};

#endif // ROOT_TPie
