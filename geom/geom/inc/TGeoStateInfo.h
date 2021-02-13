// @(#):$Id$
// Author: Andrei Gheata   07/02/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoStateInfo
#define ROOT_TGeoStateInfo

#include "TGeoMatrix.h"

class TGeoNode;
class TGeoPolygon;
struct TGeoStateInfo;

struct TGeoStateInfo {
   TGeoNode            *fNode;           // Node to which applies
   // Assembly data
   Int_t                fAsmCurrent;     // Index for current entered node (assemblies)
   Int_t                fAsmNext;        // Index for next entered node (assemblies)
   // Divisions data
   Int_t                fDivCurrent;     // Index for the current division node
   Int_t                fDivNext;        // Index for the next division node
   TGeoTranslation      fDivTrans;       // Translation used by current division node
   TGeoRotation         fDivRot;         // Rotation used by current division node
   TGeoCombiTrans       fDivCombi;       // Combi transformation used by current division
   // Voxels data
   Int_t                fVoxNcandidates; // Number of candidates
   Int_t                fVoxCurrent;     // Index of current voxel in sorted list
   Int_t               *fVoxCheckList;   // List of candidates
   UChar_t             *fVoxBits1;       // Bits used for list intersection
   Int_t                fVoxSlices[3];   // Slice indices for current voxel
   Int_t                fVoxInc[3];      // Slice index increment
   Double_t             fVoxInvdir[3];   // 1/current director cosines
   Double_t             fVoxLimits[3];   // Limits on X,Y,Z
   // Composite shape data
   Int_t                fBoolSelected;   // Selected Boolean node
   // Xtru shape data
   Int_t                fXtruSeg;        // current segment [0,fNvert-1]
   Int_t                fXtruIz;         // current z plane [0,fNz-1]
   Double_t            *fXtruXc;         // [fNvert] current X positions for polygon vertices
   Double_t            *fXtruYc;         // [fNvert] current Y positions for polygon vertices
   TGeoPolygon         *fXtruPoly;       // polygon defining section shape

   TGeoStateInfo(Int_t maxdaughters=0);
   TGeoStateInfo(const TGeoStateInfo &other) = delete;
   TGeoStateInfo & operator=(const TGeoStateInfo &other) = delete;
   virtual ~TGeoStateInfo();

   ClassDef(TGeoStateInfo, 0)  // No I/O for this structure
};

#endif
