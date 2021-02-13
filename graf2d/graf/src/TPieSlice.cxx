/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <TPieSlice.h>

#include <TError.h>
#include <TVirtualPad.h>
#include <TPie.h>

ClassImp(TPieSlice);

/** \class TPieSlice
\ingroup BasicGraphics

A slice of a piechart, see the TPie class.

This class describe the property of single
*/

////////////////////////////////////////////////////////////////////////////////
/// This is the default constructor, used to create the standard.

TPieSlice::TPieSlice() : TNamed(), TAttFill(), TAttLine()
{
   fPie = 0;
   fValue = 1;
   fRadiusOffset = 0;
   fIsActive = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// This constructor create a slice with a particular value.

TPieSlice::TPieSlice(const char *name, const char *title,
                     TPie *pie, Double_t val) :
                     TNamed(name, title), TAttFill(), TAttLine()
{
   fPie = pie;
   fValue = val;
   fRadiusOffset = 0;
   fIsActive = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Eval if the mouse is over the area associated with this slice.

Int_t TPieSlice::DistancetoPrimitive(Int_t /*px*/, Int_t /*py*/)
{
   Int_t dist = 9999;

   if (fIsActive) {
      dist = 0;
      fIsActive = kFALSE;
      gPad->SetCursor(kHand);
   }

   return dist;
}

////////////////////////////////////////////////////////////////////////////////
/// return the value of the offset in radial direction for this slice.

Double_t TPieSlice::GetRadiusOffset()
{
   return fRadiusOffset;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the value of this slice.

Double_t TPieSlice::GetValue()
{
   return fValue;
}

////////////////////////////////////////////////////////////////////////////////
/// Do nothing.

void TPieSlice::SavePrimitive(std::ostream &/*out*/, Option_t * /*opts*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Set the radial offset of this slice.

void TPieSlice::SetRadiusOffset(Double_t val)
{
   fRadiusOffset = val;
   if (fRadiusOffset<.0) fRadiusOffset = .0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the value for this slice.
/// Negative values are changed with its absolute value.

void TPieSlice::SetValue(Double_t val)
{
   fValue = val;
   if (fValue<.0) {
      Warning("SetValue","Invalid negative value. Absolute value taken");
      fValue *= -1;
   }

   fPie->MakeSlices(kTRUE);
}
