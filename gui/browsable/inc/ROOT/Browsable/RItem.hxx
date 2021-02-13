/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_Browsable_RItem
#define ROOT7_Browsable_RItem

#include <string>

namespace ROOT {
namespace Experimental {
namespace Browsable {

/** \class RItem
\ingroup rbrowser
\brief Representation of single item in the browser
*/

class RItem {
protected:
   std::string name;     ///< item name
   int nchilds{0};       ///< number of childs
   std::string icon;     ///< icon associated with item
   std::string title;    ///< item title
   bool checked{false};  ///< is checked, not used yet
   bool expanded{false}; ///< is expanded, not used yet
public:

   RItem() = default;
   RItem(const std::string &_name, int _nchilds = 0, const std::string &_icon = "") : name(_name), nchilds(_nchilds), icon(_icon) {}
   // must be here, one needs virtual table for correct streaming of sub-classes
   virtual ~RItem() = default;

   const std::string &GetName() const { return name; }
   const std::string &GetIcon() const { return icon; }
   const std::string &GetTitle() const { return title; }
   virtual bool IsFolder() const { return false; }
   virtual bool IsHidden() const { return false; }

   void SetChecked(bool on = true) { checked = on; }
   void SetExpanded(bool on = true) { expanded = on; }
   void SetIcon(const std::string &_icon) { icon = _icon; }
   void SetTitle(const std::string &_title) { title = _title; }

   virtual bool Compare(const RItem *b, const std::string &) const
   {
      if (IsFolder() != b->IsFolder())
         return IsFolder();
      return GetName() < b->GetName();
   }
};


} // namespace Browsable
} // namespace Experimental
} // namespace ROOT

#endif


