/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAttrValue
#define ROOT7_RAttrValue

#include <ROOT/RAttrBase.hxx>

namespace ROOT {
namespace Experimental {

/** \class RAttrValue
\ingroup GpadROOT7
\author Sergey Linev <s.linev@gsi.de>
\date 2020-06-24
\brief Template class to access single value from drawable or other attributes
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/


template<typename T>
class RAttrValue : public RAttrBase {
protected:

   RAttrMap  fDefaults;    ///<!    map with default values

   const RAttrMap &GetDefaults() const override { return fDefaults; }

   bool IsValue() const override { return true; }

public:

   RAttrValue() = default;

   RAttrValue(RDrawable *drawable, const std::string &name, const T &dflt = T())
   {
      fDefaults.AddValue("", dflt);
      AssignDrawable(drawable, name);
   }

   RAttrValue(RAttrBase *parent, const std::string &name, const T &dflt = T())
   {
      fDefaults.AddValue("", dflt);
      AssignParent(parent, name);
   }

   void Set(const T &v) { SetValue("", v); }
   T Get() const { return GetValue<T>(""); }
   void Clear() { ClearValue(""); }
   bool Has() const { return HasValue<T>(""); }

   RAttrValue &operator=(const T &v) { Set(v); return *this; }

   operator T() const { return Get(); }
};

} // namespace Experimental
} // namespace ROOT

#endif // ROOT7_RAttrValue
