// @(#):$Id$
// Author: Andrei Gheata   07/02/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGeoStateInfo.h"

#include "Rtypes.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoPolygon.h"

/** \class TGeoStateInfo
\ingroup Geometry_classes
Statefull info for the current geometry level.
*/

ClassImp(TGeoStateInfo);

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TGeoStateInfo::TGeoStateInfo(Int_t maxdaughters)
              :fNode(0),
               fAsmCurrent(0),
               fAsmNext(0),
               fDivCurrent(0),
               fDivNext(0),
               fDivTrans(),
               fDivRot(),
               fDivCombi(),
               fVoxNcandidates(0),
               fVoxCurrent(0),
               fVoxCheckList(0),
               fVoxBits1(0),
               fBoolSelected(0),
               fXtruSeg(0),
               fXtruIz(0),
               fXtruXc(0),
               fXtruYc(0),
               fXtruPoly(0)
{
   Int_t maxDaughters = (maxdaughters>0) ? maxdaughters : TGeoManager::GetMaxDaughters();
   Int_t maxXtruVert  = TGeoManager::GetMaxXtruVert();
   fVoxCheckList = new Int_t[maxDaughters];
   fVoxBits1 = new UChar_t[2 + ((maxDaughters-1)>>3)];
   fXtruXc = new Double_t[maxXtruVert];
   fXtruYc = new Double_t[maxXtruVert];
   fVoxSlices[0] = fVoxSlices[1] = fVoxSlices[2] = -1;
   fVoxInc[0] = fVoxInc[1] = fVoxInc[2] = 0;
   fVoxInvdir[0] = fVoxInvdir[1] = fVoxInvdir[2] = 0;
   fVoxLimits[0] = fVoxLimits[1] = fVoxLimits[2] = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoStateInfo::~TGeoStateInfo()
{
   delete [] fVoxCheckList;
   delete [] fVoxBits1;
   delete [] fXtruXc;
   delete [] fXtruYc;
}
