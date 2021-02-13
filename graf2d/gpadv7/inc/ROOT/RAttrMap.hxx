/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAttrMap
#define ROOT7_RAttrMap

#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>

#include <ROOT/RPadLength.hxx>
#include <ROOT/RMakeUnique.hxx>

namespace ROOT {
namespace Experimental {

class RAttrBase;
class RStyle;

/** \class RAttrMap
\ingroup GpadROOT7
\authors Axel Naumann <axel@cern.ch> Sergey Linev <s.linev@gsi.de>
\date 2017-09-26
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RAttrMap {

   friend class RAttrBase;
   friend class RStyle;

public:

   enum EValuesKind { kNoValue, kBool, kInt, kDouble, kString };

   class Value_t {
   public:
      virtual ~Value_t() = default;
      virtual EValuesKind Kind() const = 0;
      virtual bool CanConvertFrom(EValuesKind kind) const { return kind == Kind(); }
      virtual bool CanConvertTo(EValuesKind kind) const { return kind == Kind(); }
      virtual bool GetBool() const { return false; }
      virtual int GetInt() const { return 0; }
      virtual double GetDouble() const { return 0; }
      virtual std::string GetString() const { return ""; }
      virtual bool IsEqual(const Value_t &) const { return false; }
      virtual std::unique_ptr<Value_t> Copy() const = 0;

      template<typename T> T Get() const;

      template <typename RET_TYPE, typename MATCH_TYPE = void>
      static RET_TYPE GetValue(const Value_t *rec);
   };

   class NoValue_t : public Value_t {
   public:
      explicit NoValue_t() {}
      EValuesKind Kind() const final { return kNoValue; }
      std::unique_ptr<Value_t> Copy() const final { return std::make_unique<NoValue_t>(); }
      bool IsEqual(const Value_t &tgt) const final { return (tgt.Kind() == kNoValue); }
   };

   class BoolValue_t : public Value_t {
      bool v{false}; ///< integer value
   public:
      explicit BoolValue_t(bool _v = false) : v(_v) {}
      EValuesKind Kind() const final { return kBool; }
      bool CanConvertFrom(EValuesKind kind) const final { return (kind == kDouble) || (kind == kInt) || (kind == kBool) || (kind == kString); }
      bool CanConvertTo(EValuesKind kind) const final { return (kind == kDouble) || (kind == kInt) || (kind == kBool) || (kind == kString); }
      bool GetBool() const final { return v; }
      int GetInt() const final { return v ? 1 : 0; }
      double GetDouble() const final { return v ? 1 : 0; }
      std::string GetString() const final { return v ? "true" : "false"; }
      std::unique_ptr<Value_t> Copy() const final { return std::make_unique<BoolValue_t>(v); }
      bool IsEqual(const Value_t &tgt) const final { return tgt.GetBool() == v; }
   };

   class IntValue_t : public Value_t {
      int v{0}; ///< integer value
   public:
      IntValue_t(int _v = 0) : v(_v) {}
      EValuesKind Kind() const final { return kInt; }
      bool CanConvertFrom(EValuesKind kind) const final { return (kind == kInt) || (kind == kBool); }
      bool CanConvertTo(EValuesKind kind) const final { return (kind == kDouble) || (kind == kInt) || (kind == kBool) || (kind == kString); }
      bool GetBool() const final { return v ? true : false; }
      int GetInt() const final { return v; }
      double GetDouble() const final { return v; }
      std::string GetString() const final { return std::to_string(v); }
      std::unique_ptr<Value_t> Copy() const final { return std::make_unique<IntValue_t>(v); }
      bool IsEqual(const Value_t &tgt) const final { return tgt.GetInt() == v; }
   };

   class DoubleValue_t : public Value_t {
      double v{0}; ///< double value
   public:
      DoubleValue_t(double _v = 0) : v(_v) {}
      EValuesKind Kind() const final { return kDouble; }
      bool CanConvertFrom(EValuesKind kind) const final { return (kind == kDouble) || (kind == kInt) || (kind == kBool); }
      bool CanConvertTo(EValuesKind kind) const final { return (kind == kDouble) || (kind == kBool) || (kind == kString); }
      bool GetBool() const final { return v ? true : false; }
      int GetInt() const final { return (int) v; }
      double GetDouble() const final { return v; }
      std::string GetString() const final { return std::to_string(v); }
      std::unique_ptr<Value_t> Copy() const final { return std::make_unique<DoubleValue_t>(v); }
      bool IsEqual(const Value_t &tgt) const final { return tgt.GetDouble() == v; }
   };

   class StringValue_t : public Value_t {
      std::string v; ///< string value
   public:
      StringValue_t(const std::string _v = "") : v(_v) {}
      EValuesKind Kind() const final { return kString; }
      // all values can be converted into the string
      bool CanConvertFrom(EValuesKind) const final { return true; }
      // can convert into string and boolean
      bool CanConvertTo(EValuesKind kind) const final { return kind == kString; }
      bool GetBool() const final { return v.compare("true") == 0; }
      std::string GetString() const final { return v; }
      bool IsEqual(const Value_t &tgt) const final { return tgt.GetString() == v; }
      std::unique_ptr<Value_t> Copy() const final { return std::make_unique<StringValue_t>(v); }
   };

private:

   // FIXME: due to ROOT-10306 only data member of such kind can be correctly stored by ROOT I/O
   // Once problem fixed, one could make this container a base class
   std::unordered_map<std::string, std::unique_ptr<Value_t>> m; ///< JSON_object

   void AddBestMatch(const std::string &name, const std::string &value);

public:

   RAttrMap() = default; ///< JSON_asbase - store as map object

   RAttrMap &Add(const std::string &name, std::unique_ptr<Value_t> &&value) { m[name] = std::move(value); return *this; }
   RAttrMap &AddNoValue(const std::string &name) { m[name] = std::make_unique<NoValue_t>(); return *this; }
   RAttrMap &AddBool(const std::string &name, bool value) { m[name] = std::make_unique<BoolValue_t>(value); return *this; }
   RAttrMap &AddInt(const std::string &name, int value) { m[name] = std::make_unique<IntValue_t>(value); return *this; }
   RAttrMap &AddDouble(const std::string &name, double value) { m[name] = std::make_unique<DoubleValue_t>(value); return *this; }
   RAttrMap &AddString(const std::string &name, const std::string &value) { m[name] = std::make_unique<StringValue_t>(value); return *this; }
   RAttrMap &AddPadLength(const std::string &name, const RPadLength &value)
   {
      if (value.Empty())
         Clear(name);
      else
         m[name] = std::make_unique<StringValue_t>(value.AsString());
      return *this;
   }
   RAttrMap &AddDefaults(const RAttrBase &vis);

   RAttrMap &AddValue(const std::string &name, bool value) { return AddBool(name, value); }
   RAttrMap &AddValue(const std::string &name, int value) { return AddInt(name, value); }
   RAttrMap &AddValue(const std::string &name, double value) { return AddDouble(name, value); }
   RAttrMap &AddValue(const std::string &name, const std::string &value) { return AddString(name, value); }
   RAttrMap &AddValue(const std::string &name, const RPadLength &value) { return AddPadLength(name, value); }

   RAttrMap(const RAttrMap &src)
   {
      for (const auto &pair : src.m)
         m[pair.first] = pair.second->Copy();
   }

   RAttrMap &operator=(const RAttrMap &src)
   {
      m.clear();
      for (const auto &pair : src.m)
         m[pair.first] = pair.second->Copy();
      return *this;
   }

   const Value_t *Find(const std::string &name) const
   {
      auto entry = m.find(name);
      return (entry != m.end()) ? entry->second.get() : nullptr;
   }

   /** Clear specified attribute */
   void Clear(const std::string &name)
   {
      auto entry = m.find(name);
      if (entry != m.end())
         m.erase(entry);
   }

   bool Change(const std::string &name, Value_t *value = nullptr);

   auto begin() const { return m.begin(); }
   auto end() const { return m.end(); }
};

template<> bool RAttrMap::Value_t::Get<bool>() const;
template<> int RAttrMap::Value_t::Get<int>() const;
template<> double RAttrMap::Value_t::Get<double>() const;
template<> std::string RAttrMap::Value_t::Get<std::string>() const;
template<> RPadLength RAttrMap::Value_t::Get<RPadLength>() const;

template<> bool RAttrMap::Value_t::GetValue<bool,void>(const Value_t *rec);
template<> int RAttrMap::Value_t::GetValue<int,void>(const Value_t *rec);
template<> double RAttrMap::Value_t::GetValue<double,void>(const Value_t *rec);
template<> std::string RAttrMap::Value_t::GetValue<std::string,void>(const Value_t *rec);
template<> RPadLength RAttrMap::Value_t::GetValue<RPadLength,void>(const Value_t *rec);

template<> const RAttrMap::Value_t *RAttrMap::Value_t::GetValue<const RAttrMap::Value_t *,void>(const Value_t *rec);
template<> const RAttrMap::Value_t *RAttrMap::Value_t::GetValue<const RAttrMap::Value_t *,bool>(const Value_t *rec);
template<> const RAttrMap::Value_t *RAttrMap::Value_t::GetValue<const RAttrMap::Value_t *,int>(const Value_t *rec);
template<> const RAttrMap::Value_t *RAttrMap::Value_t::GetValue<const RAttrMap::Value_t *,double>(const Value_t *rec);
template<> const RAttrMap::Value_t *RAttrMap::Value_t::GetValue<const RAttrMap::Value_t *,std::string>(const Value_t *rec);
template<> const RAttrMap::Value_t *RAttrMap::Value_t::GetValue<const RAttrMap::Value_t *,RPadLength>(const Value_t *rec);

} // namespace Experimental
} // namespace ROOT

#endif // ROOT7_RAttrMap
