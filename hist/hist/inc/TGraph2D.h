// @(#)root/hist:$Id: TGraph2D.h,v 1.00
// Author: Olivier Couet

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGraph2D
#define ROOT_TGraph2D


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGraph2D                                                             //
//                                                                      //
// Graph 2D graphics class.                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TAttMarker.h"

class TAxis;
class TList;
class TF2;
class TH1;
class TH2;
class TH2D;
class TView;
class TDirectory;
class TVirtualHistPainter;

#include "TFitResultPtr.h"

class TGraph2D : public TNamed, public TAttLine, public TAttFill, public TAttMarker {

protected:

   Int_t       fNpoints;          ///< Number of points in the data set
   Int_t       fNpx;              ///< Number of bins along X in fHistogram
   Int_t       fNpy;              ///< Number of bins along Y in fHistogram
   Int_t       fMaxIter;          ///< Maximum number of iterations to find Delaunay triangles
   Int_t       fSize;             ///<!Real size of fX, fY and fZ
   Double_t   *fX;                ///<[fNpoints]
   Double_t   *fY;                ///<[fNpoints] Data set to be plotted
   Double_t   *fZ;                ///<[fNpoints]
   Double_t    fMinimum;          ///< Minimum value for plotting along z
   Double_t    fMaximum;          ///< Maximum value for plotting along z
   Double_t    fMargin;           ///< Extra space (in %) around interpolated area for fHistogram
   Double_t    fZout;             ///< fHistogram bin height for points lying outside the interpolated area
   TList      *fFunctions;        ///< Pointer to list of functions (fits and user)
   TH2D       *fHistogram;        ///<!2D histogram of z values linearly interpolated on the triangles
   TObject    *fDelaunay;         ///<! Pointer to Delaunay interpolator object
   TDirectory *fDirectory;        ///<!Pointer to directory holding this 2D graph
   TVirtualHistPainter *fPainter; ///<!Pointer to histogram painter

   void     Build(Int_t n);

private:

   Bool_t      fUserHisto;   // True when SetHistogram has been called

   enum EStatusBits {
      kOldInterpolation =  BIT(15)
   };

   void CreateInterpolator(Bool_t oldInterp);

protected:

public:

   TGraph2D();
   TGraph2D(Int_t n);
   TGraph2D(Int_t n, Int_t *x, Int_t *y, Int_t *z);
   TGraph2D(Int_t n, Float_t *x, Float_t *y, Float_t *z);
   TGraph2D(Int_t n, Double_t *x, Double_t *y, Double_t *z);
   TGraph2D(TH2 *h2);
   TGraph2D(const char *name, const char *title, Int_t n, Double_t *x, Double_t *y, Double_t *z);
   TGraph2D(const char *filename, const char *format="%lg %lg %lg", Option_t *option="");
   TGraph2D(const TGraph2D &);

   virtual ~TGraph2D();

   TGraph2D& operator=(const TGraph2D &);

   virtual void          AddPoint(Double_t x, Double_t y, Double_t z) { SetPoint(fNpoints, x, y, z); } ///< Append a new point to the graph.
   virtual void          Browse(TBrowser *);
   virtual void          Clear(Option_t *option="");
   virtual void          DirectoryAutoAdd(TDirectory *);
   Int_t                 DistancetoPrimitive(Int_t px, Int_t py);
   virtual void          Draw(Option_t *option="P0");
   void                  ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual TObject      *FindObject(const char *name) const;
   virtual TObject      *FindObject(const TObject *obj) const;
   virtual TFitResultPtr Fit(const char *formula ,Option_t *option="" ,Option_t *goption=""); // *MENU*
   virtual TFitResultPtr Fit(TF2 *f2 ,Option_t *option="" ,Option_t *goption=""); // *MENU*
   virtual void          FitPanel(); // *MENU*
   TList                *GetContourList(Double_t contour);
   TDirectory           *GetDirectory() const {return fDirectory;}
   Int_t                 GetNpx() const {return fNpx;}
   Int_t                 GetNpy() const {return fNpy;}
   TH2D                 *GetHistogram(Option_t *option="");
   TList                *GetListOfFunctions() const { return fFunctions; }
   virtual Double_t      GetErrorX(Int_t bin) const;
   virtual Double_t      GetErrorY(Int_t bin) const;
   virtual Double_t      GetErrorZ(Int_t bin) const;
   Double_t              GetMargin() const {return fMargin;}
   Double_t              GetMaximum() const {return fMaximum;};
   Double_t              GetMinimum() const {return fMinimum;};
   TAxis                *GetXaxis() const ;
   TAxis                *GetYaxis() const ;
   TAxis                *GetZaxis() const ;
   Int_t                 GetN() const {return fNpoints;}
   Double_t             *GetX() const {return fX;}
   Double_t             *GetY() const {return fY;}
   Double_t             *GetZ() const {return fZ;}
   virtual Double_t     *GetEX() const {return 0;}
   virtual Double_t     *GetEY() const {return 0;}
   virtual Double_t     *GetEZ() const {return 0;}
   Double_t              GetXmax() const;
   Double_t              GetXmin() const;
   Double_t              GetYmax() const;
   Double_t              GetYmin() const;
   Double_t              GetZmax() const;
   Double_t              GetZmin() const;
   virtual Double_t      GetXmaxE() const {return GetXmax();};
   virtual Double_t      GetXminE() const {return GetXmin();};
   virtual Double_t      GetYmaxE() const {return GetYmax();};
   virtual Double_t      GetYminE() const {return GetYmin();};
   virtual Double_t      GetZmaxE() const {return GetZmax();};
   virtual Double_t      GetZminE() const {return GetZmin();};
   virtual Int_t         GetPoint(Int_t i, Double_t &x, Double_t &y, Double_t &z) const;
   Double_t              Interpolate(Double_t x, Double_t y);
   void                  Paint(Option_t *option="");
   virtual void          Print(Option_t *chopt="") const;
   TH1                  *Project(Option_t *option="x") const; // *MENU*
   Int_t                 RemovePoint(Int_t ipoint); // *MENU*
   virtual void          SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void          Scale(Double_t c1=1., Option_t *option="z"); // *MENU*
   virtual void          Set(Int_t n);
   virtual void          SetDirectory(TDirectory *dir);
   virtual void          SetHistogram(TH2 *h);
   void                  SetMargin(Double_t m=0.1); // *MENU*
   void                  SetMarginBinsContent(Double_t z=0.); // *MENU*
   void                  SetMaximum(Double_t maximum=-1111); // *MENU*
   void                  SetMinimum(Double_t minimum=-1111); // *MENU*
   void                  SetMaxIter(Int_t n=100000) {fMaxIter = n;} // *MENU*
   virtual void          SetName(const char *name); // *MENU*
   virtual void          SetNameTitle(const char *name, const char *title);
   void                  SetNpx(Int_t npx=40); // *MENU*
   void                  SetNpy(Int_t npx=40); // *MENU*
   virtual void          SetPoint(Int_t point, Double_t x, Double_t y, Double_t z); // *MENU*
   virtual void          SetTitle(const char *title=""); // *MENU*


   ClassDef(TGraph2D,1)  //Set of n x[n],y[n],z[n] points with 3-d graphics including Delaunay triangulation
};

#endif
