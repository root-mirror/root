// Author: Brian Bockelman UNL 09/2017

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/TIOFeatures.hxx"
#include "TBasket.h"
#include "TBranch.h"
#include "TEnum.h"
#include "TEnumConstant.h"
#include "TError.h"
#include "TTree.h"

#include <bitset>
#include <sstream>

using namespace ROOT;

/**
 * \class ROOT::TIOFeatures
 * \ingroup tree
 *
 * `TIOFeatures` provides the end-user with the ability to change the IO behavior
 * of data written via a `TTree`.  This class allows access to experimental and non-default
 * features.
 *
 * When one of these features are activated, forward compatibility breaks may occur.
 * That is, older versions of ROOT may not be able to read files written by this version
 * of ROOT that have enabled these non-default features.
 *
 * To utilize `TIOFeatures`, create the object, set the desired feature flags, then attach
 * it to a `TTree`.  All subsequently created branches (and their baskets) will be serialized
 * using those particular features.
 *
 * Example usage:
 * ~~~{.cpp}
 * ROOT::TIOFeatures features;
 * features.Set(ROOT::Experimental::EIOFeatures::kGenerateOffsetMap);
 * ttree_ref.SetIOFeatures(features);
 * ~~~
 *
 * The method `TTree::SetIOFeatures` creates a copy of the feature set; subsequent changes
 * to the `TIOFeatures` object do not propogate to the `TTree`.
 */


////////////////////////////////////////////////////////////////////////////
/// \brief Clear a specific IO feature from this set.
/// \param[in] input_bits The specific feature to disable.
///
/// Removes a feature from the `TIOFeatures` object; emits an Error message if
/// the IO feature is not supported by this version of ROOT.
void TIOFeatures::Clear(Experimental::EIOFeatures input_bits)
{
   Clear(static_cast<EIOFeatures>(input_bits));
}

////////////////////////////////////////////////////////////////////////////
/// \brief Clear a specific IO feature from this set.
/// \param[in] input_bits The specific feature to disable.
///
/// Removes a feature from the `TIOFeatures` object; emits an Error message if
/// the IO feature is not supported by this version of ROOT.
void TIOFeatures::Clear(Experimental::EIOUnsupportedFeatures input_bits)
{
   Clear(static_cast<EIOFeatures>(input_bits));
}

////////////////////////////////////////////////////////////////////////////
/// \brief Clear a specific IO feature from this set.
/// \param[in] input_bits The specific feature to disable.
///
/// Removes a feature from the `TIOFeatures` object; emits an Error message if
/// the IO feature is not supported by this version of ROOT.
void TIOFeatures::Clear(EIOFeatures input_bits)
{
   TBasket::EIOBits enum_bits = static_cast<TBasket::EIOBits>(input_bits);
   auto bits = static_cast<UChar_t>(enum_bits);
   if (R__unlikely((bits & static_cast<UChar_t>(TBasket::EIOBits::kSupported)) != bits)) {
      Error("TestFeature", "A feature is being cleared that is not supported.");
      return;
   }
   fIOBits &= ~bits;
}

static std::string GetUnsupportedName(TBasket::EUnsupportedIOBits enum_flag)
{
   UChar_t flag = static_cast<UChar_t>(enum_flag);

   std::string retval = "unknown";

   TClass *cl = TBasket::Class();
   if (cl == nullptr) {
      return retval;
   }

   TEnum *eUnsupportedIOBits = (TEnum *)cl->GetListOfEnums()->FindObject("EUnsupportedIOBits");
   if (eUnsupportedIOBits == nullptr) {
      return retval;
   }

   for (auto constant : ROOT::Detail::TRangeStaticCast<TEnumConstant>(eUnsupportedIOBits->GetConstants())) {
      if (constant->GetValue() == flag) {
         retval = constant->GetName();
         break;
      }
   }
   return retval;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Set a specific IO feature.
/// \param[in] input_bits The specific feature to enable.
///
/// Sets a feature in the `TIOFeatures` object; emits an Error message if
/// the IO feature is not supported by this version of ROOT.
///
/// If the feature is supported by ROOT, this function returns kTRUE; otherwise,
/// it returns kFALSE.
bool TIOFeatures::Set(Experimental::EIOFeatures input_bits)
{
   return Set(static_cast<EIOFeatures>(input_bits));
}

////////////////////////////////////////////////////////////////////////////
/// \brief Set a specific IO feature.
/// \param[in] input_bits The specific feature to enable.
///
/// Sets a feature in the `TIOFeatures` object; emits an Error message if
/// the IO feature is not supported by this version of ROOT.
///
/// If the feature is supported by ROOT, this function returns kTRUE; otherwise,
/// it returns kFALSE.
bool TIOFeatures::Set(EIOFeatures input_bits)
{
   TBasket::EIOBits enum_bits = static_cast<TBasket::EIOBits>(input_bits);
   auto bits = static_cast<UChar_t>(enum_bits);
   if (R__unlikely((bits & static_cast<UChar_t>(TBasket::EIOBits::kSupported)) != bits)) {
      UChar_t unsupported = bits & static_cast<UChar_t>(TBasket::EUnsupportedIOBits::kUnsupported);
      if (unsupported) {
         Error("SetFeature", "A feature was request (%s) but this feature is no longer supported.",
               GetUnsupportedName(static_cast<TBasket::EUnsupportedIOBits>(unsupported)).c_str());
      } else {
         Error("SetFeature", "An unknown feature was requested (flag=%s); cannot enable it.",
               std::bitset<32>(unsupported).to_string().c_str());
      }
      return kFALSE;
   }
   fIOBits |= bits;
   return kTRUE;
}


/////////////////////////////////////////////////////////////////////////////
/// \brief Given a IO feature string, set the corresponding feature
/// \param [in] value Feature name to test.
///
/// This allows one to set a feature given a specific string from the
/// TBasket::EIOBits enum.
///
/// *NOTE* this function is quite slow and users are strongly encouraged to
/// use the type-safe `Set` version instead.  This has been added for better
/// CLI interfaces.
///
/// Returns kTRUE only if a new feature was set; otherwise emits an error message
/// and returns kFALSE.
bool TIOFeatures::Set(const std::string &value)
{
   TClass *cl = TBasket::Class();
   if (cl == nullptr) {
      Error("Set", "Could not retrieve TBasket's class");
      return kFALSE;
   }
   TEnum *eIOBits = static_cast<TEnum*>(cl->GetListOfEnums()->FindObject("EIOBits"));
   if (eIOBits == nullptr) {
      Error("Set", "Could not locate TBasket::EIOBits enum");
      return kFALSE;
   }
   for (auto constant : ROOT::Detail::TRangeStaticCast<TEnumConstant>(eIOBits->GetConstants())) {
      if (!strcmp(constant->GetName(), value.c_str())) {
         return Set(static_cast<EIOFeatures>(constant->GetValue()));
      }
   }
   Error("Set", "Could not locate %s in TBasket::EIOBits", value.c_str());
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////
/// \brief Print a human-readable representation of the TIOFeatures to stdout
///
/// Prints a string with the names of all enabled IO features.
void TIOFeatures::Print() const
{
   TClass *cl = TBasket::Class();
   if (cl == nullptr) {
      Error("Print", "Could not retrieve TBasket's class");
      return;
   }
   TEnum *eIOBits = static_cast<TEnum *>(cl->GetListOfEnums()->FindObject("EIOBits"));
   if (eIOBits == nullptr) {
      Error("Print", "Could not locate TBasket::EIOBits enum");
      return;
   }
   std::stringstream ss;
   bool hasFeatures = false;
   ss << "TIOFeatures{";
   for (auto constant : ROOT::Detail::TRangeStaticCast<TEnumConstant>(eIOBits->GetConstants())) {
      if ((constant->GetValue() & fIOBits) == constant->GetValue()) {
         ss << (hasFeatures ? ", " : "") << constant->GetName();
         hasFeatures = true;
      }
   }
   ss << "}";
   Printf("%s", ss.str().c_str());
}

////////////////////////////////////////////////////////////////////////////
/// \brief Test to see if a given feature is set
/// \param[in] input_bits The specific feature to test.
///
/// Returns kTRUE if the feature is enables in this object and supported by
/// this version of ROOT.
bool TIOFeatures::Test(Experimental::EIOFeatures input_bits) const
{
   return Test(static_cast<EIOFeatures>(input_bits));
}

////////////////////////////////////////////////////////////////////////////
/// \brief Test to see if a given feature is set
/// \param[in] input_bits The specific feature to test.
///
/// Returns kTRUE if the feature is enables in this object and supported by
/// this version of ROOT.
bool TIOFeatures::Test(EIOFeatures input_bits) const
{
   TBasket::EIOBits enum_bits = static_cast<TBasket::EIOBits>(input_bits);
   auto bits = static_cast<UChar_t>(enum_bits);
   if (R__unlikely((bits & static_cast<UChar_t>(TBasket::EIOBits::kSupported)) != bits)) {
      Error("TestFeature", "A feature is being tested for that is not supported or known.");
      return kFALSE;
   }
   return (fIOBits & bits) == bits;
}

UChar_t TIOFeatures::GetFeatures() const
{
   return fIOBits;
}
