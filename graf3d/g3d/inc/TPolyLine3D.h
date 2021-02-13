// @(#)root/g3d:$Id$
// Author: Nenad Buncic   17/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPolyLine3D
#define ROOT_TPolyLine3D


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPolyLine3D                                                          //
//                                                                      //
// A 3-D polyline.                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TString.h"
#include "TAttLine.h"
#include "TAtt3D.h"

class TList;


class TPolyLine3D : public TObject, public TAttLine, public TAtt3D
{
protected:
   Int_t        fN;            ///< Number of points
   Float_t     *fP;            ///< [3*fN] Array of 3-D coordinates  (x,y,z)
   TString      fOption;       ///< options
   Int_t        fLastPoint;    ///< The index of the last filled point

public:
   TPolyLine3D();
   TPolyLine3D(Int_t n, Option_t *option="");
   TPolyLine3D(Int_t n, Float_t const* p, Option_t *option="");
   TPolyLine3D(Int_t n, Double_t const* p, Option_t *option="");
   TPolyLine3D(Int_t n, Float_t const* x, Float_t const* y, Float_t const* z, Option_t *option="");
   TPolyLine3D(Int_t n, Double_t const* x, Double_t const* y, Double_t const* z, Option_t *option="");
   TPolyLine3D(const TPolyLine3D &polylin);
   TPolyLine3D& operator=(const TPolyLine3D &polylin);
   virtual ~TPolyLine3D();

   virtual void      Copy(TObject &polyline) const;
   virtual Int_t     DistancetoPrimitive(Int_t px, Int_t py);
   virtual void      Draw(Option_t *option="");
   virtual void      DrawPolyLine(Int_t n, Float_t *p, Option_t *option="");
   virtual void      ExecuteEvent(Int_t event, Int_t px, Int_t py);
   Int_t             GetLastPoint() const {return fLastPoint;}
   Int_t             GetN() const {return fN;}
   Float_t          *GetP() const {return fP;}
   Option_t         *GetOption() const {return fOption.Data();}
   virtual void      ls(Option_t *option="") const;
   virtual Int_t     Merge(TCollection *list);
   virtual void      Paint(Option_t *option="");
   virtual void      Print(Option_t *option="") const;
   virtual void      SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual Int_t     SetNextPoint(Double_t x, Double_t y, Double_t z); // *MENU*
   virtual void      SetOption(Option_t *option="") {fOption = option;}
   virtual void      SetPoint(Int_t point, Double_t x, Double_t y, Double_t z); // *MENU*
   virtual void      SetPolyLine(Int_t n, Option_t *option="");
   virtual void      SetPolyLine(Int_t n, Float_t *p, Option_t *option="");
   virtual void      SetPolyLine(Int_t n, Double_t *p, Option_t *option="");
   virtual Int_t     Size() const { return fLastPoint+1;}

   static  void      DrawOutlineCube(TList *outline, Double_t *rmin, Double_t *rmax);

   ClassDef(TPolyLine3D,1)  //A 3-D polyline
};

#endif
