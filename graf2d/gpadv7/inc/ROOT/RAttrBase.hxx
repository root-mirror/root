/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAttrBase
#define ROOT7_RAttrBase

#include <ROOT/RAttrMap.hxx>
#include <ROOT/RStyle.hxx>
#include <ROOT/RDrawable.hxx>

namespace ROOT {
namespace Experimental {

class RLogChannel;
/// Log channel for GPad diagnostics.
RLogChannel &GPadLog();

/** \class RAttrBase
\ingroup GpadROOT7
\author Sergey Linev <s.linev@gsi.de>
\date 2019-09-17
\brief Base class for all attributes, used with RDrawable
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RAttrBase {

   friend class RAttrMap;

   enum {kDrawable, kParent, kOwnAttr} fKind{kDrawable}; ///<!  kind of data

   union {
      RDrawable *drawable;  // either drawable to which attributes belongs to
      RAttrBase *parent;    // or aggregation of attributes
      RAttrMap  *ownattr;   // or just own container with values
   } fD{nullptr};  ///<!  data

   std::string fPrefix;                ///<! name prefix for all attributes values

   void ClearData();
   RAttrMap *CreateOwnAttr();

protected:

   RDrawable *GetDrawable() const { return fKind == kDrawable ? fD.drawable : nullptr; }
   RAttrBase *GetParent() const { return fKind == kParent ? fD.parent : nullptr; }
   RAttrMap *GetOwnAttr() const { return fKind == kOwnAttr ? fD.ownattr : nullptr; }

   void AssignDrawable(RDrawable *drawable, const std::string &prefix);

   void AssignParent(RAttrBase *parent, const std::string &prefix);


   virtual const RAttrMap &GetDefaults() const;

   bool CopyValue(const std::string &name, const RAttrMap::Value_t &value, bool check_type = true);

   bool IsValueEqual(const std::string &name, const RAttrMap::Value_t &value, bool use_style = false) const;

   ///////////////////////////////////////////////////////////////////////////////

   struct Rec_t {
      RAttrMap *attr{nullptr};
      std::string fullname;
      RDrawable *drawable{nullptr};
      operator bool() const { return !!attr; }
   };

   /// Find attributes container and full-qualified name for value
   const Rec_t AccessAttr(const std::string &name) const
   {
      const RAttrBase *prnt = this;
      std::string fullname = name;
      while (prnt) {
         fullname.insert(0, prnt->fPrefix); // fullname = prnt->fPrefix + fullname
         if (auto dr = prnt->GetDrawable())
            return { &dr->fAttr, fullname, dr };
         if (auto attr = prnt->GetOwnAttr())
            return { attr, fullname, nullptr };
         prnt = prnt->GetParent();
      }
      return {nullptr, fullname, nullptr};
   }

   struct Val_t {
      const RAttrMap::Value_t *value{nullptr};
      std::shared_ptr<RStyle> style;
      operator bool() const { return !!value; }
   };

   /** Search value with given name in attributes */
   const Val_t AccessValue(const std::string &name, bool use_style = true) const
   {
      if (auto access = AccessAttr(name)) {
         if (auto rec = access.attr->Find(access.fullname))
            return {rec, nullptr};
         if (access.drawable && use_style)
            if (auto observe = access.drawable->fStyle.lock()) {
               if (auto rec = observe->Eval(access.fullname, *access.drawable))
                  return {rec, observe};
            }
      }

      return {nullptr, nullptr};
   }

   /// Ensure attribute with give name exists - creates container for attributes if required

   Rec_t EnsureAttr(const std::string &name)
   {
      auto prnt = this;
      std::string fullname = name;
      while (prnt) {
         fullname.insert(0, prnt->fPrefix); // fullname = prnt->fPrefix + fullname
         if (auto dr = prnt->GetDrawable())
            return { &dr->fAttr, fullname, dr };
         if (prnt->fKind != kParent)
            return { prnt->CreateOwnAttr(), fullname, nullptr };
         prnt = prnt->GetParent();
      }
      return {nullptr, fullname, nullptr};
   }

   /// Evaluate attribute value

   template <typename RET_TYPE,typename MATCH_TYPE = void>
   auto Eval(const std::string &name, bool use_dflts = true) const
   {
      if (auto v = AccessValue(name, true))
         return RAttrMap::Value_t::GetValue<RET_TYPE,MATCH_TYPE>(v.value);

      const RAttrMap::Value_t *rec = nullptr;

      if (use_dflts)
         rec = GetDefaults().Find(name);

      return RAttrMap::Value_t::GetValue<RET_TYPE,MATCH_TYPE>(rec);
   }

   void CopyTo(RAttrBase &tgt, bool use_style = true) const;

   void MoveTo(RAttrBase &tgt);

   bool IsSame(const RAttrBase &src, bool use_style = true) const;

   RAttrBase(RDrawable *drawable, const std::string &prefix) { AssignDrawable(drawable, prefix); }

   RAttrBase(RAttrBase *parent, const std::string &prefix) { AssignParent(parent, prefix); }

   RAttrBase(const RAttrBase &src) { src.CopyTo(*this); }

   RAttrBase &operator=(const RAttrBase &src)
   {
      Clear();
      src.CopyTo(*this);
      return *this;
   }

   void SetNoValue(const std::string &name);
   void SetValue(const std::string &name, bool value);
   void SetValue(const std::string &name, double value);
   void SetValue(const std::string &name, int value);
   void SetValue(const std::string &name, const std::string &value);
   void SetValue(const std::string &name, const RPadLength &value);
   void SetValue(const std::string &name, const RColor &value);

   const std::string &GetPrefix() const { return fPrefix; }

   void ClearValue(const std::string &name);

   void Clear();

   template <typename T = void>
   bool HasValue(const std::string &name, bool check_defaults = false) const
   {
      auto res = Eval<const RAttrMap::Value_t *, T>(name, check_defaults);
      return res ? (res->Kind() != RAttrMap::kNoValue) : false;
   }

   template <typename T>
   T GetValue(const std::string &name) const
   {
      return Eval<T>(name);
   }

   virtual RAttrMap CollectDefaults() const;

   virtual bool IsValue() const { return false; }

public:
   RAttrBase() = default;

   virtual ~RAttrBase() { ClearData(); }

   friend bool operator==(const RAttrBase& lhs, const RAttrBase& rhs) { return lhs.IsSame(rhs) && rhs.IsSame(lhs); }
   friend bool operator!=(const RAttrBase& lhs, const RAttrBase& rhs) { return !lhs.IsSame(rhs) || !rhs.IsSame(lhs); }

};


} // namespace Experimental
} // namespace ROOT

#define R__ATTR_CLASS(ClassName,dflt_prefix) \
protected: \
const RAttrMap &GetDefaults() const override \
{ \
   static auto dflts = CollectDefaults(); \
   return dflts; \
} \
public: \
   ClassName() = default; \
   ClassName(RDrawable *drawable, const std::string &prefix = dflt_prefix) { AssignDrawable(drawable, prefix); } \
   ClassName(RAttrBase *parent, const std::string &prefix = dflt_prefix) { AssignParent(parent, prefix); } \
   ClassName(const ClassName &src) : ClassName() { src.CopyTo(*this); } \
   ClassName(ClassName &&src) : ClassName() { src.MoveTo(*this); } \
   ClassName &operator=(ClassName &&src) { src.MoveTo(*this); return *this; } \
   ClassName &operator=(const ClassName &src) { Clear(); src.CopyTo(*this); return *this; } \

#endif
