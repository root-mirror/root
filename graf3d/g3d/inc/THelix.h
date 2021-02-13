// @(#)root/g3d:$Id$
// Author: Ping Yeh   19/12/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_THelix
#define ROOT_THelix

#include "TPolyLine3D.h"
#include "TRotMatrix.h"

enum EHelixRangeType {
   kHelixT, kHelixX, kHelixY, kHelixZ, kLabX, kLabY, kLabZ, kUnchanged
};


class THelix : public TPolyLine3D {

protected:
   Double_t     fX0;       //Pivot's x position (see parametrization in class doc)
   Double_t     fY0;       //Pivot's y position (see parametrization in class doc)
   Double_t     fZ0;       //Pivot's z position (see parametrization in class doc)
   Double_t     fVt;       //Transverse velocity (constant of motion)
   Double_t     fPhi0;     //Initial phase, so vx0 = fVt*cos(fPhi0)
   Double_t     fVz;       //Z velocity (constant of motion)
   Double_t     fW;        //Angular frequency
   Double_t     fAxis[3];  //Direction unit vector of the helix axis
   TRotMatrix  *fRotMat;   //Rotation matrix: axis // z  -->  axis // fAxis
   Double_t     fRange[2]; //Range of helix parameter t

   THelix& operator=(const THelix&);

   void         SetRotMatrix();    //Set rotation matrix
   Double_t     FindClosestPhase(Double_t phi0,  Double_t cosine);

   static Int_t fgMinNSeg;   //minimal number of segments in polyline

public:
   THelix();
   THelix(Double_t x,  Double_t y,  Double_t z,
          Double_t vx, Double_t vy, Double_t vz,
          Double_t w);
   THelix(Double_t const* xyz, Double_t const* v, Double_t w,
          Double_t const* range=0, EHelixRangeType rtype=kHelixZ,
          Double_t const* axis=0);
   THelix(const THelix &helix);
   virtual ~THelix();

   virtual void    Copy(TObject &helix) const;
   virtual void    Draw(Option_t *option="");
   Option_t       *GetOption() const {return fOption.Data();}
   virtual void    Print(Option_t *option="") const;
   virtual void    SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void    SetOption(Option_t *option="") {fOption = option;}
   virtual void    SetAxis(Double_t const* axis);       //Define new axis
   virtual void    SetAxis(Double_t x, Double_t y, Double_t z);
   virtual void    SetRange(Double_t * range, EHelixRangeType rtype=kHelixZ);
   virtual void    SetRange(Double_t r1, Double_t r2, EHelixRangeType rtype=kHelixZ);
   void            SetHelix(Double_t const* xyz,  Double_t const* v, Double_t w,
                            Double_t const* range=0, EHelixRangeType type=kUnchanged,
                            Double_t const* axis=0);

   ClassDef(THelix,2)  //A Helix drawn as a PolyLine3D
};

#endif
