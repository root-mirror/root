/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RAttrMap.hxx"

#include "ROOT/RAttrBase.hxx"
#include "ROOT/RLogger.hxx"

#include <string>
#include <algorithm>
#include <limits>

using namespace ROOT::Experimental;

using namespace std::string_literals;

template<> bool RAttrMap::Value_t::Get<bool>() const { return GetBool(); }
template<> int RAttrMap::Value_t::Get<int>() const { return GetInt(); }
template<> double RAttrMap::Value_t::Get<double>() const { return GetDouble(); }
template<> std::string RAttrMap::Value_t::Get<std::string>() const { return GetString(); }
template<> RPadLength RAttrMap::Value_t::Get<RPadLength>() const { return GetString(); }

template<> bool RAttrMap::Value_t::GetValue<bool,void>(const Value_t *rec) { return rec ? rec->GetBool() : false; }
template<> int RAttrMap::Value_t::GetValue<int,void>(const Value_t *rec) { return rec ? rec->GetInt() : 0; }
template<> double RAttrMap::Value_t::GetValue<double,void>(const Value_t *rec) { return rec ? rec->GetDouble() : 0.; }
template<> std::string RAttrMap::Value_t::GetValue<std::string,void>(const Value_t *rec) { return rec ? rec->GetString() : ""s; }
template<> RPadLength RAttrMap::Value_t::GetValue<RPadLength,void>(const Value_t *rec) { return rec ? rec->GetString() : ""s; }

template<> const RAttrMap::Value_t *RAttrMap::Value_t::GetValue<const RAttrMap::Value_t *,void>(const Value_t *rec) { return rec; }
template<> const RAttrMap::Value_t *RAttrMap::Value_t::GetValue<const RAttrMap::Value_t *,bool>(const Value_t *rec) { return rec && rec->CanConvertTo(RAttrMap::kBool) ? rec : nullptr; }
template<> const RAttrMap::Value_t *RAttrMap::Value_t::GetValue<const RAttrMap::Value_t *,int>(const Value_t *rec) { return rec && rec->CanConvertTo(RAttrMap::kInt) ? rec : nullptr; }
template<> const RAttrMap::Value_t *RAttrMap::Value_t::GetValue<const RAttrMap::Value_t *,double>(const Value_t *rec) { return rec && rec->CanConvertTo(RAttrMap::kDouble) ? rec : nullptr; }
template<> const RAttrMap::Value_t *RAttrMap::Value_t::GetValue<const RAttrMap::Value_t *,std::string>(const Value_t *rec) { return rec && rec->CanConvertTo(RAttrMap::kString) ? rec : nullptr; }
template<> const RAttrMap::Value_t *RAttrMap::Value_t::GetValue<const RAttrMap::Value_t *,RPadLength>(const Value_t *rec) { return rec && rec->CanConvertTo(RAttrMap::kString) ? rec : nullptr; }


/////////////////////////////////////////////////////////////////////////////////////////////////
/// Add defaults values form sub attribute

RAttrMap &RAttrMap::AddDefaults(const RAttrBase &vis)
{
   auto prefix = vis.GetPrefix();

   for (const auto &entry : vis.GetDefaults())
      m[prefix+entry.first] = entry.second->Copy();

   return *this;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Add attribute, converting to best possible type
/// Tested boolean, int, double. If none works - store as a string

void RAttrMap::AddBestMatch(const std::string &name, const std::string &value)
{
   if ((value == "none"s) || (value == "null"s) || value.empty()) {
      AddNoValue(name);
      return;
   }

   if (value == "true"s) {
      AddBool(name, true);
      return;
   }

   if (value == "false"s) {
      AddBool(name, false);
      return;
   }

   auto beg = value.begin();
   int base = 10;
   bool int_conversion_fails = false;

   if (*beg == '-') {
      ++beg;
   } else if ((value.length() > 2) && (*beg == '0') && (value[1] == 'x')) {
      beg += 2;
      base = 16;
   }

   // check if only digits are present
   if (std::find_if(beg, value.end(), [](unsigned char c) { return !std::isdigit(c); }) == value.end()) {

      try {

         auto ivalue = std::stoll(base == 16 ? value.substr(2) : value, nullptr, base);

         if ((ivalue >= std::numeric_limits<int>::min()) && (ivalue <= std::numeric_limits<int>::max()))
            AddInt(name, ivalue);
         else
            AddDouble(name, ivalue);

         return;
      } catch (...) {
         // int conversion fails
         int_conversion_fails = true;
      }
   }

   // check if characters for double is present
   if (!int_conversion_fails && std::find_if(beg, value.end(), [](unsigned char c) {
                                   return !std::isdigit(c) && (c != '.') && (c != '-') && (c != '+') && (c != 'e');
                                }) == value.end()) {
      try {
         double dvalue = std::stod(value);
         AddDouble(name, dvalue);
         return;
      } catch (...) {
         // do nothing
      }
   }

   AddString(name, value);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Change attribute using string value and kind
/// Used to change attributes from JS side
/// Returns true if value was really changed

bool RAttrMap::Change(const std::string &name, Value_t *value)
{
   auto entry = m.find(name);
   if ((entry == m.end()) || (entry->second->Kind() == kNoValue)) {
      if (!value) return false;
      m[name] = value->Copy();
      return true;
   }

   // specify nullptr means clear attribute
   if (!value) {
      m.erase(entry);
      return true;
   }

   // error situation - conversion cannot be performed
   if(!value->CanConvertTo(entry->second->Kind())) {
      R__LOG_ERROR(GPadLog()) << "Wrong data type provided for attribute " << name;
      return false;
   }

   // no need to change something
   if (entry->second->IsEqual(*value))
      return false;

   switch (entry->second->Kind()) {
      case kNoValue: break; // just to avoid compiler warnings
      case kBool: AddBool(name, value->GetBool()); break;
      case kInt: AddInt(name, value->GetInt()); break;
      case kDouble: AddDouble(name, value->GetDouble()); break;
      case kString: AddString(name, value->GetString()); break;
   }

   return true;
}
