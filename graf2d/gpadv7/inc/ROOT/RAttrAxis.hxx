/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAttrAxis
#define ROOT7_RAttrAxis

#include <ROOT/RAttrBase.hxx>
#include <ROOT/RAttrLine.hxx>
#include <ROOT/RAttrText.hxx>
#include <ROOT/RAttrValue.hxx>
#include <ROOT/RPadLength.hxx>

namespace ROOT {
namespace Experimental {

/** \class RAttrAxis
\ingroup GpadROOT7
\author Sergey Linev <s.linev@gsi.de>
\date 2020-02-20
\brief All supported axes attributes for: line, ticks, labels, title, min/max, log, reverse, ...
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RAttrAxis : public RAttrBase {

   RAttrValue<double> fZoomMin{this, "zoommin", 0.};                     ///<! axis zoom min
   RAttrValue<double> fZoomMax{this, "zoommax", 0.};                     ///<! axis zoom max
   RAttrValue<double> fLog{this, "log", 0};                              ///<! log scale, <1 off, 1 - base10, 2 - base 2, 2.71 - exp, 3, 4, ...
   RAttrValue<bool> fReverse{this, "reverse", false};                    ///<! reverse scale
   RAttrValue<bool> fTimeDisplay{this, "time", false};                   ///<! time display
   RAttrValue<double> fTimeOffset{this, "time_offset", 0};               ///<! time offset to display
   RAttrValue<std::string> fTimeFormat{this, "time_format", ""};         ///<! time format
   RAttrLine fAttrLine{this, "line"};                                    ///<! line attributes
   RAttrValue<std::string> fEndingStyle{this, "ending_style", ""};       ///<! axis ending style - none, arrow, circle
   RAttrValue<RPadLength> fEndingSize{this, "ending_size", 0.02_normal}; ///<! axis ending size
   RAttrValue<std::string> fTicksSide{this, "ticks_side", "normal"};     ///<! ticks position - normal, invert, both
   RAttrValue<RPadLength> fTicksSize{this, "ticks_size", 0.02_normal};   ///<! ticks size
   RAttrColor fTicksColor{this, "ticks_color"};                          ///<! ticks color
   RAttrText fLabelsAttr{this, "labels"};                                ///<! text attributes for labels
   RAttrValue<RPadLength> fLabelsOffset{this, "labels_offset", {}};      ///<! axis labels offset - relative
   RAttrValue<bool> fLabelsCenter{this, "labels_center", false};         ///<! center labels
   RAttrText fTitleAttr{this, "title"};                                  ///<! axis title text attributes
   RAttrValue<std::string> fTitle{this, "title", ""};                    ///<! axis title
   RAttrValue<std::string> fTitlePos{this, "title_position", "right"};   ///<! axis title position - left, right, center
   RAttrValue<RPadLength> fTitleOffset{this, "title_offset", {}};        ///<! axis title offset - relative

   R__ATTR_CLASS(RAttrAxis, "axis");

   RAttrAxis &SetZoomMin(double min) { fZoomMin = min; return *this; }
   RAttrAxis &SetZoomMax(double max) { fZoomMax = max; return *this; }
   double GetZoomMin() const { return fZoomMin; }
   double GetZoomMax() const { return fZoomMax; }
   bool HasZoomMin() const { return fZoomMin.Has(); }
   bool HasZoomMax() const { return fZoomMax.Has(); }

   RAttrAxis &SetZoom(double min, double max) { SetZoomMin(min); SetZoomMax(max); return *this; }
   void ClearZoom() { fZoomMin.Clear(); fZoomMax.Clear(); }

   RAttrAxis &SetLog(double base = 10) { fLog = (base < 1) ? 0 : base; return *this; }
   double GetLog() const { return fLog; }

   RAttrAxis &SetReverse(bool on = true) { fReverse = on; return *this; }
   bool GetReverse() const { return fReverse; }

   RAttrAxis &SetTimeDisplay(const std::string &fmt = "", double offset = -1)
   {
      fTimeDisplay = true;
      if (!fmt.empty()) fTimeFormat = fmt;
      if (offset >= 0) fTimeOffset = offset;
      return *this;
   }
   bool GetTimeDisplay() const { return fTimeDisplay; }
   void ClearTimeDisplay()
   {
      fTimeDisplay.Clear();
      fTimeOffset.Clear();
      fTimeFormat.Clear();
   }

   const RAttrLine &GetAttrLine() const { return fAttrLine; }
   RAttrAxis &SetAttrLine(const RAttrLine &line) { fAttrLine = line; return *this; }
   RAttrLine &AttrLine() { return fAttrLine; }

   RAttrAxis &SetEndingSize(const RPadLength &sz) { fEndingSize = sz; return *this; }
   RPadLength GetEndingSize() const { return fEndingSize; }

   RAttrAxis &SetEndingStyle(const std::string &st) { fEndingStyle = st; return *this; }
   RAttrAxis &SetEndingArrow() { return SetEndingStyle("arrow"); }
   RAttrAxis &SetEndingCircle() { return SetEndingStyle("cicrle"); }
   std::string GetEndingStyle() const { return fEndingStyle; }

   RAttrAxis &SetTicksSize(const RPadLength &sz) { fTicksSize = sz; return *this; }
   RPadLength GetTicksSize() const { return fTicksSize; }

   RAttrAxis &SetTicksSide(const std::string &side) { fTicksSide = side; return *this; }
   RAttrAxis &SetTicksNormal() { return SetTicksSide("normal"); }
   RAttrAxis &SetTicksInvert() { return SetTicksSide("invert"); }
   RAttrAxis &SetTicksBoth() { return SetTicksSide("both"); }
   std::string GetTicksSide() const { return fTicksSide; }

   RAttrAxis &SetTicksColor(const RColor &color) { fTicksColor = color; return *this; }
   RColor GetTicksColor() const { return fTicksColor.GetColor(); }
   RAttrColor &AttrTiksColor() { return fTicksColor; }

   RAttrAxis &SetLabelsOffset(const RPadLength &len) { fLabelsOffset = len; return *this; }
   RPadLength GetLabelsOffset() const { return fLabelsOffset; }

   const RAttrText &GetLabelsAttr() const { return fLabelsAttr; }
   RAttrAxis &SetLabelsAttr(const RAttrText &attr) { fLabelsAttr = attr; return *this; }
   RAttrText &LabelsAttr() { return fLabelsAttr; }

   RAttrAxis &SetLabelsCenter(bool on = true) { fLabelsCenter = on; return *this; }
   bool GetLabelsCenter() const { return fLabelsCenter; }

   const RAttrText &GetTitleAttr() const { return fTitleAttr; }
   RAttrAxis &SetTitleAttr(const RAttrText &attr) { fTitleAttr = attr; return *this; }
   RAttrText &TitleAttr() { return fTitleAttr; }

   RAttrAxis &SetTitle(const std::string &title) { fTitle = title; return *this; }
   std::string GetTitle() const { return fTitle; }

   RAttrAxis &SetTitlePos(const std::string &pos) { fTitlePos = pos; return *this; }
   RAttrAxis &SetTitleLeft() { return SetTitlePos("left"); }
   RAttrAxis &SetTitleCenter() { return SetTitlePos("center"); }
   RAttrAxis &SetTitleRight() { return SetTitlePos("right"); }
   std::string GetTitlePos() const { return fTitlePos; }

   RAttrAxis &SetTitleOffset(const RPadLength &len) { fTitleOffset = len; return *this; }
   RPadLength GetTitleOffset() const { return fTitleOffset; }
};

} // namespace Experimental
} // namespace ROOT

#endif
