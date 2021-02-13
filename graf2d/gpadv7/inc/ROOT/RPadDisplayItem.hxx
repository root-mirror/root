/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RPadDisplayItem
#define ROOT7_RPadDisplayItem

#include <ROOT/RDisplayItem.hxx>
#include <ROOT/RPad.hxx>
#include "ROOT/RStyle.hxx"

namespace ROOT {
namespace Experimental {


/** class RPadBaseDisplayItem
\ingroup GpadROOT7
\brief Display item for the RPadBase class, includes primitives, attributes and frame
\author Sergey Linev
\date 2017-05-31
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RPadBaseDisplayItem : public RDisplayItem {
public:
   // list of snapshot for primitives in pad
   using PadPrimitives_t = std::vector<std::unique_ptr<RDisplayItem>>;

protected:
   const RAttrMap *fAttr{nullptr};      ///< temporary pointer on attributes
   PadPrimitives_t fPrimitives;         ///< display items for all primitives in the pad
   std::vector<std::shared_ptr<RStyle>> fStyles; ///<! locked styles of the objects and pad until streaming is performed
public:
   RPadBaseDisplayItem() = default;
   virtual ~RPadBaseDisplayItem() = default;
   void SetAttributes(const RAttrMap *f) { fAttr = f; }
   /// Add display item and style which should be used for it
   void Add(std::unique_ptr<RDisplayItem> &&item, std::shared_ptr<RStyle> &&style)
   {
      if (style) {
         item->SetStyle(style.get());
         fStyles.emplace_back(std::move(style));
      }
      fPrimitives.push_back(std::move(item));
   }
   /// Assign style for the pad
   void SetPadStyle(std::shared_ptr<RStyle> &&style)
   {
      if (style) {
         SetStyle(style.get());
         fStyles.emplace_back(std::move(style));
      }
   }
};

/** class RPadDisplayItem
\ingroup GpadROOT7
\brief Display item for the RPad class, add pad position and size
\author Sergey Linev
\date 2017-05-31
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RPadDisplayItem : public RPadBaseDisplayItem {

protected:
   const RPadPos *fPos{nullptr};        ///< pad position
   const RPadExtent *fSize{nullptr};    ///< pad size
public:
   RPadDisplayItem() = default;
   virtual ~RPadDisplayItem() {}
   void SetPadPosSize(const RPadPos *pos, const RPadExtent *size) { fPos = pos; fSize = size; }

   void BuildFullId(const std::string &prefix) override
   {
      RDisplayItem::BuildFullId(prefix);
      std::string subprefix = prefix + std::to_string(GetIndex()) + "_";
      for (auto &item : fPrimitives)
         item->BuildFullId(subprefix);
   }
};


/** class RCanvasDisplayItem
\ingroup GpadROOT7
\brief Display item for the RCanvas class, add canvas title and size
\author Sergey Linev
\date 2017-05-31
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RCanvasDisplayItem : public RPadBaseDisplayItem {

protected:
   std::string fTitle;                  ///< title of the canvas
   std::array<int, 2> fWinSize;         ///< canvas window size
public:
   RCanvasDisplayItem() = default;
   virtual ~RCanvasDisplayItem() = default;
   void SetTitle(const std::string &title) { fTitle = title; }
   void SetWindowSize(const std::array<RPadLength::Pixel, 2> &win) { fWinSize[0] = (int) win[0].fVal; fWinSize[1] = (int) win[1].fVal; }

   void BuildFullId(const std::string &prefix) override
   {
      for (auto &item : fPrimitives)
         item->BuildFullId(prefix);
   }
};



} // Experimental
} // ROOT

#endif
