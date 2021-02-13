/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAxisDrawable
#define ROOT7_RAxisDrawable

#include <ROOT/RDrawable.hxx>
#include <ROOT/RAttrAxis.hxx>
#include <ROOT/RPadPos.hxx>
#include <ROOT/RPadLength.hxx>
#include <ROOT/RLogger.hxx>

#include <vector>
#include <string>

namespace ROOT {
namespace Experimental {

/** \class RAxisDrawableBase
\ingroup GrafROOT7
\brief Axis base drawing - only attributes and position.
\author Sergey Linev <S.Linev@gsi.de>
\date 2020-11-03
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RAxisDrawableBase : public RDrawable {

   RPadPos fPos;                          ///< axis start point
   bool fVertical{false};                 ///< is vertical axis
   RPadLength fLength;                    ///< axis length
   RAttrAxis fAttrAxis{this, "axis"};     ///<! axis attributes

public:

   RAxisDrawableBase() : RDrawable("axis") {}

   RAxisDrawableBase(const RPadPos &pos, bool vertical, const RPadLength &len) : RAxisDrawableBase()
   {
      SetPos(pos);
      SetVertical(vertical);
      SetLength(len);
   }

   RAxisDrawableBase &SetPos(const RPadPos &pos) { fPos = pos; return *this; }
   RAxisDrawableBase &SetVertical(bool vertical = true) { fVertical = vertical; return *this; }
   RAxisDrawableBase &SetLength(const RPadLength &len) { fLength = len; return *this; }

   bool IsVertical() const { return fVertical; }
   const RPadPos& GetPos() const { return fPos; }
   const RPadLength& GetLength() const { return fLength; }

   const RAttrAxis &GetAttrAxis() const { return fAttrAxis; }
   RAxisDrawableBase &SetAttrAxis(const RAttrAxis &attr) { fAttrAxis = attr; return *this; }
   RAttrAxis &AttrAxis() { return fAttrAxis; }
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////

/** \class RAxisDrawable
\ingroup GrafROOT7
\brief Plain axis drawing.
\author Sergey Linev <S.Linev@gsi.de>
\date 2020-11-03
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RAxisDrawable : public RAxisDrawableBase {

   double fMin{0}, fMax{100};      ///< axis minimum and maximum

public:

   RAxisDrawable() = default;

   RAxisDrawable(const RPadPos &pos, bool vertical, const RPadLength &len) : RAxisDrawableBase(pos, vertical, len) {}

   RAxisDrawable &SetMinMax(double min, double max) { fMin = min; fMax = max; return *this; }
   double GetMin() const { return fMin; }
   double GetMax() const { return fMax; }
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////

/** \class RAxisLabelsDrawable
\ingroup GrafROOT7
\brief Labels axis drawing.
\author Sergey Linev <S.Linev@gsi.de>
\date 2020-11-03
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RAxisLabelsDrawable : public RAxisDrawableBase {

   std::vector<std::string> fLabels;      ///< axis labels

public:

   RAxisLabelsDrawable() = default;

   RAxisLabelsDrawable(const RPadPos &pos, bool vertical, const RPadLength &len) : RAxisDrawableBase(pos, vertical, len) {}

   RAxisLabelsDrawable &SetLabels(const std::vector<std::string> &lbls) { fLabels = lbls; return *this; }
   const std::vector<std::string> &GetLabels() const { return fLabels; }
};


} // namespace Experimental
} // namespace ROOT

#endif
