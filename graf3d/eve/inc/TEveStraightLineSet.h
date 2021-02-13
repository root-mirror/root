// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveStraightLineSet
#define ROOT_TEveStraightLineSet

#include "TEveUtil.h"

#include "TNamed.h"
#include "TQObject.h"
#include "TAtt3D.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TAttBBox.h"

#include "TEveElement.h"
#include "TEveProjectionBases.h"
#include "TEveChunkManager.h"
#include "TEveTrans.h"

class TRandom;

class TEveStraightLineSet : public TEveElement,
                            public TEveProjectable,
                            public TNamed,
                            public TQObject,
                            public TAtt3D,
                            public TAttLine,
                            public TAttMarker,
                            public TAttBBox
{
private:
   TEveStraightLineSet(const TEveStraightLineSet&);            // Not implemented
   TEveStraightLineSet& operator=(const TEveStraightLineSet&); // Not implemented

public:
   struct Line_t
   {
      Int_t          fId;
      Float_t        fV1[3];
      Float_t        fV2[3];
      TRef           fRef;

      Line_t(Float_t x1, Float_t y1, Float_t z1,
             Float_t x2, Float_t y2, Float_t z2) : fId(-1), fRef()
      {
         fV1[0] = x1, fV1[1] = y1, fV1[2] = z1;
         fV2[0] = x2, fV2[1] = y2, fV2[2] = z2;
      }
   };

   struct Marker_t
   {
      Float_t      fV[3];
      Int_t        fLineId;
      TRef         fRef;

      Marker_t(Float_t x, Float_t y, Float_t z, Int_t line_id) : fLineId(line_id), fRef()
      {
         fV[0] = x, fV[1] = y, fV[2] = z;
      }
   };

protected:
   TEveChunkManager  fLinePlex;
   TEveChunkManager  fMarkerPlex;

   Bool_t            fOwnLinesIds;    // Flag specifying if id-objects are owned by the line-set
   Bool_t            fOwnMarkersIds;  // Flag specifying if id-objects are owned by the line-set

   Bool_t            fRnrMarkers;
   Bool_t            fRnrLines;

   Bool_t            fDepthTest;

   Line_t*           fLastLine; //!

public:
   TEveStraightLineSet(const char* n="StraightLineSet", const char* t="");
   virtual ~TEveStraightLineSet() {}

   virtual void SetLineColor(Color_t col) { SetMainColor(col); }

   Line_t*   AddLine(Float_t x1, Float_t y1, Float_t z1, Float_t x2, Float_t y2, Float_t z2);
   Line_t*   AddLine(const TEveVector& p1, const TEveVector& p2);
   Marker_t* AddMarker(Float_t x, Float_t y, Float_t z, Int_t line_id=-1);
   Marker_t* AddMarker(const TEveVector& p, Int_t line_id=-1);
   Marker_t* AddMarker(Int_t line_id, Float_t pos);

   void      SetLine(int idx, Float_t x1, Float_t y1, Float_t z1, Float_t x2, Float_t y2, Float_t z2);
   void      SetLine(int idx, const TEveVector& p1, const TEveVector& p2);

   TEveChunkManager& GetLinePlex()   { return fLinePlex;   }
   TEveChunkManager& GetMarkerPlex() { return fMarkerPlex; }

   virtual Bool_t GetRnrMarkers() { return fRnrMarkers; }
   virtual Bool_t GetRnrLines()   { return fRnrLines;   }
   virtual Bool_t GetDepthTest()  { return fDepthTest;   }

   virtual void SetRnrMarkers(Bool_t x) { fRnrMarkers = x; }
   virtual void SetRnrLines(Bool_t x)   { fRnrLines   = x; }
   virtual void SetDepthTest(Bool_t x)  { fDepthTest   = x; }

   virtual void CopyVizParams(const TEveElement* el);
   virtual void WriteVizParams(std::ostream& out, const TString& var);

   virtual TClass* ProjectedClass(const TEveProjection* p) const;

   virtual void ComputeBBox();
   virtual void Paint(Option_t* option="");

   ClassDef(TEveStraightLineSet, 0); // Set of straight lines with optional markers along the lines.
};


/******************************************************************************/

class TEveStraightLineSetProjected : public TEveStraightLineSet,
                                     public TEveProjected
{
private:
   TEveStraightLineSetProjected(const TEveStraightLineSetProjected&);            // Not implemented
   TEveStraightLineSetProjected& operator=(const TEveStraightLineSetProjected&); // Not implemented

protected:
   virtual void SetDepthLocal(Float_t d);

public:
   TEveStraightLineSetProjected();
   virtual ~TEveStraightLineSetProjected() {}

   virtual void SetProjection(TEveProjectionManager* mng, TEveProjectable* model);
   virtual void UpdateProjection();
   virtual TEveElement* GetProjectedAsElement() { return this; }

   ClassDef(TEveStraightLineSetProjected, 0); // Projected copy of a TEveStraightLineSet.
};

#endif
