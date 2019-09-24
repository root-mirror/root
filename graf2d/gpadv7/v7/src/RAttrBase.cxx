/// \file ROOT/RAttrBase.cxx
/// \ingroup Gpad ROOT7
/// \author Sergey Linev <s.linev@gsi.de>
/// \date 2019-09-17
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RAttrBase.hxx>

#include <ROOT/RLogger.hxx>


///////////////////////////////////////////////////////////////////////////////
/// Return default values for attributes, empty for base class

const ROOT::Experimental::RAttrMap &ROOT::Experimental::RAttrBase::GetDefaults() const
{
   static RAttrMap empty;
   return empty;
}

///////////////////////////////////////////////////////////////////////////////
/// Copy attributes from other object

bool ROOT::Experimental::RAttrBase::CopyValue(const std::string &name, const RAttrMap::Value_t &value, bool check_type)
{
   if (check_type) {
      const auto *dvalue = GetDefaults().Find(name);
      if (!dvalue || !dvalue->Compatible(value.Kind()))
         return false;
   }

   if (auto access = EnsureAttr(name)) {
      access.attr->Add(access.fullname, value.Copy());
      return true;
   }

   return false;
}

///////////////////////////////////////////////////////////////////////////////
/// Copy attributes into target object

bool ROOT::Experimental::RAttrBase::IsValueEqual(const std::string &name, const RAttrMap::Value_t &value, bool use_style) const
{
   if (auto v = AccessValue(name, use_style))
      return v.value->IsEqual(value);

   return false;
}

///////////////////////////////////////////////////////////////////////////////
/// Copy attributes into target object

void ROOT::Experimental::RAttrBase::CopyTo(RAttrBase &tgt, bool use_style) const
{
   for (const auto &entry : GetDefaults()) {
      if (auto v = AccessValue(entry.first, use_style))
         tgt.CopyValue(entry.first, *v.value);
   }
}

///////////////////////////////////////////////////////////////////////////////
/// Check if all values which are evaluated in this object are exactly the same as in tgt object

bool ROOT::Experimental::RAttrBase::IsSame(const RAttrBase &tgt, bool use_style) const
{
   for (const auto &entry : GetDefaults()) {
      if (auto v = AccessValue(entry.first, use_style))
         if (!tgt.IsValueEqual(entry.first, *v.value, use_style)) return false;
   }
   return true;
}

///////////////////////////////////////////////////////////////////////////////
/// Return value from attributes container - no style or defaults are used

void ROOT::Experimental::RAttrBase::AssignDrawable(RDrawable *drawable, const std::string &prefix)
{
   fDrawable = drawable;
   fOwnAttr.reset();
   fPrefix = prefix;
   fParent = nullptr;
}

void ROOT::Experimental::RAttrBase::AssignParent(const RAttrBase *parent, const std::string &prefix)
{
   fDrawable = nullptr;
   fOwnAttr.reset();
   fPrefix = prefix;
   fParent = parent;
}

void ROOT::Experimental::RAttrBase::ClearValue(const std::string &name)
{
   if (auto access = AccessAttr(name))
       const_cast<RAttrMap*>(access.attr)->Clear(access.fullname);
}

void ROOT::Experimental::RAttrBase::SetValue(const std::string &name, int value)
{
   if (auto access = EnsureAttr(name))
      access.attr->AddInt(access.fullname, value);
}

void ROOT::Experimental::RAttrBase::SetValue(const std::string &name, double value)
{
   if (auto access = EnsureAttr(name))
      access.attr->AddDouble(access.fullname, value);
}

void ROOT::Experimental::RAttrBase::SetValue(const std::string &name, const std::string &value)
{
   if (auto access = EnsureAttr(name))
      access.attr->AddString(access.fullname, value);
}

/** Clear all respective values from drawable. Only defaults can be used */
void ROOT::Experimental::RAttrBase::Clear()
{
   for (const auto &entry : GetDefaults())
      ClearValue(entry.first);
}
