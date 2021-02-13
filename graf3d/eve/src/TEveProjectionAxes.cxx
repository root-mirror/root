// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
* Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/

#include "TEveProjectionAxes.h"
#include "TEveProjectionManager.h"


/** \class TEveProjectionAxes
\ingroup TEve
Axes for non-linear projections. Show scale of TEveProjectionManager
children. With different step mode tick-marks can positioned
equidistant or placed with value monotonically increasing from lower left corner
of bounding box.
*/

ClassImp(TEveProjectionAxes);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveProjectionAxes::TEveProjectionAxes(TEveProjectionManager* m, Bool_t useCS) :
   TEveElement(),
   TNamed("TEveProjectionAxes", ""),
   fManager(m),

   fUseColorSet(useCS),

   fLabMode(kValue),
   fAxesMode(kAll),

   fDrawCenter(kFALSE),
   fDrawOrigin(kFALSE)
{
   fCanEditMainTrans = kFALSE;
   fManager->AddDependent(this);

   // Axis attributes.
   fNdivisions = 1010;
   fLabelSize = 0.015;
   fLabelColor = kGray+1;
   fAxisColor = kGray+1;
   fTickLength = 0.015;
   fLabelOffset = 0.01;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TEveProjectionAxes::~TEveProjectionAxes()
{
   fManager->RemoveDependent(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this object. Only direct rendering is supported.

void TEveProjectionAxes::Paint(Option_t*)
{
   PaintStandard(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from TAttBBox; fill bounding-box information.

void TEveProjectionAxes::ComputeBBox()
{
   static const TEveException eH("TEveProjectionManager::ComputeBBox ");

   BBoxZero();
   if(fManager == 0 || fManager->GetBBox() == 0)
      return;

   for (Int_t i=0; i<6; ++i)
      fBBox[i] = fManager->GetBBox()[i];

   AssertBBoxExtents(0.1);
}

////////////////////////////////////////////////////////////////////////////////
/// Return TEveProjectionAxes icon.

const TGPicture* TEveProjectionAxes::GetListTreeIcon(Bool_t)
{
   return TEveElement::fgListTreeIcons[6];
}
