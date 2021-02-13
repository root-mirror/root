// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveCalo2DGL
#define ROOT_TEveCalo2DGL

#include "TGLObject.h"
#include "TEveCaloData.h"
#include "TEveCalo.h"
#include <vector>

class TGLViewer;
class TGLScene;

class TEveCalo2D;
class TEveProjection;

class TEveCalo2DGL : public TGLObject
{
private:
   TEveCalo2DGL(const TEveCalo2DGL&) = delete;
   TEveCalo2DGL& operator=(const TEveCalo2DGL&) = delete;

protected:
   TEveCalo2D   *fM;  // Model object.

   void      MakeRhoZCell(Float_t thetaMin, Float_t thetaMax, Float_t& offset, Bool_t isBarrel, Bool_t phiPlus, Float_t towerH) const;

   void      MakeRPhiCell(Float_t phiMin, Float_t phiMax, Float_t towerH, Float_t offset) const;

   void      DrawRPhi(TGLRnrCtx & rnrCtx, TEveCalo2D::vBinCells_t&) const;
   void      DrawRPhiHighlighted(std::vector<TEveCaloData::vCellId_t*>& cellLists) const;
   void      DrawRhoZ(TGLRnrCtx & rnrCtx, TEveCalo2D::vBinCells_t&) const;
   void      DrawRhoZHighlighted(std::vector<TEveCaloData::vCellId_t*>& cellLists) const;

   Bool_t    IsRPhi() const;

public:
   TEveCalo2DGL();
   virtual ~TEveCalo2DGL() {}

   virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
   virtual void   SetBBox();

   virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;
   virtual void DrawHighlight(TGLRnrCtx& rnrCtx, const TGLPhysicalShape* ps, Int_t lvl=-1) const;

   // To support two-level selection
   virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
   virtual Bool_t AlwaysSecondarySelect()   const { return kTRUE; }
   virtual void ProcessSelection(TGLRnrCtx & rnrCtx, TGLSelectRecord & rec);

   ClassDef(TEveCalo2DGL, 0); // GL renderer class for TEveCalo2D.
};

#endif
