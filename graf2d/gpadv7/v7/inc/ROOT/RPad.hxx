/// \file ROOT/RPad.hxx
/// \ingroup Gpad ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \author Sergey Linev <s.linev@gsi.de>
/// \date 2017-07-06
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RPad
#define ROOT7_RPad

#include "ROOT/RPadBase.hxx"

namespace ROOT {
namespace Experimental {


/** \class ROOT::Experimental::RPad
  Graphic container for `RDrawable`-s.
  */

class RPad: public RPadBase {

   /// Pad containing this pad as a sub-pad.
   RPadBase *fParent{nullptr};             /// The parent pad, if this pad has one.

   RPadPos fPos;                           ///<! pad position
   RPadExtent fSize;                       ///<! pad size

   RAttrLine fLineAttr{this, "border_"};   ///<! border attributes

public:
   /// Create a topmost, non-paintable pad.
   RPad() = default;

   /// Create a child pad.
   RPad(RPadBase *parent, const RPadPos &pos, const RPadExtent &size): fParent(parent) { fPos = pos; fSize = size; }

   /// Destructor to have a vtable.
   virtual ~RPad();

   /// Access to the parent pad (const version).
   const RPadBase *GetParent() const { return fParent; }

   /// Access to the parent pad (non-const version).
   RPadBase *GetParent() { return fParent; }

   /// Access to the top-most canvas (const version).
   const RCanvas *GetCanvas() const override { return fParent ? fParent->GetCanvas() : nullptr; }

   /// Access to the top-most canvas (non-const version).
   RCanvas *GetCanvas() override { return fParent ? fParent->GetCanvas() : nullptr; }

   /// Get the position of the pad in parent (!) coordinates.
   const RPadPos &GetPos() const { return fPos; }

   /// Get the size of the pad in parent (!) coordinates.
   const RPadExtent &GetSize() const { return fSize; }

   /// Set the size of the pad in parent (!) coordinates.
   void SetSize(const RPadExtent &sz) { fSize = sz; }

   /// Set position
   void SetPos(const RPadPos &p) { fPos = p; }

   RAttrLine &AttrLine() { return fLineAttr; }
   const RAttrLine &AttrLine() const { return fLineAttr; }

   void Paint(Internal::RPadPainter &) final;

   /// Convert a `Pixel` position to Canvas-normalized positions.
   std::array<RPadLength::Normal, 2> PixelsToNormal(const std::array<RPadLength::Pixel, 2> &pos) const override
   {
      std::array<RPadLength::Normal, 2> posInParentNormal = fParent->PixelsToNormal(pos);
      std::array<RPadLength::Normal, 2> myPixelInNormal =
         fParent->PixelsToNormal({{fSize.Horiz().GetPixel(), fSize.Vert().GetPixel()}});
      std::array<RPadLength::Normal, 2> myUserInNormal =
         fParent->UserToNormal({{fSize.Horiz().GetUser(), fSize.Vert().GetUser()}});
      // If the parent says pos is at 0.6 in normal coords, and our size converted to normal is 0.2, then pos in our
      // coord system is 3.0!
      return {{posInParentNormal[0] / (fSize.Horiz().GetNormal() + myPixelInNormal[0] + myUserInNormal[0]),
               posInParentNormal[1] / (fSize.Vert().GetNormal() + myPixelInNormal[1] + myUserInNormal[1])}};
   }

   /// Convert a RPadPos to [x, y] of normalized coordinates.
   std::array<RPadLength::Normal, 2> ToNormal(const RPadPos &pos) const
   {
      std::array<RPadLength::Normal, 2> pixelsInNormal = PixelsToNormal({{pos.Horiz().GetPixel(), pos.Vert().GetPixel()}});
      std::array<RPadLength::Normal, 2> userInNormal = UserToNormal({{pos.Horiz().GetUser(), pos.Vert().GetUser()}});
      return {{pos.Horiz().GetNormal() + pixelsInNormal[0] + userInNormal[0],
               pos.Vert().GetNormal() + pixelsInNormal[1] + userInNormal[1]}};
   }


};

} // namespace Experimental
} // namespace ROOT

#endif
