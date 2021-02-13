// Author: Philippe Canal, 2017

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TStatusBitsChecker

   TStatusBitsChecker::Check and TStatusBitsChecker::CheckAllClasses will
   determine if the set of "status bit" declared in the class and its
   base classes presents any overlap.  The status bit are declared in
   a given class by declaring an enum type named EStatusBits.
   If some of the duplication is intentional, those duplication can
   be registered in an enum type named EStatusBitsDupExceptions.

   ~~~ {.cpp}
   // TStreamerElement status bits
   enum EStatusBits {
      kHasRange     = BIT(6),
      kCache        = BIT(9),
      kRepeat       = BIT(10),
      kRead         = BIT(11),
      kWrite        = BIT(12),
      kDoNotDelete  = BIT(13),
      kWholeObject  = BIT(14)
   };

   enum class EStatusBitsDupExceptions {
      // This bit duplicates TObject::kInvalidObject. As the semantic of kDoNotDelete is a persistent,
      // we can not change its value without breaking forward compatibility.
      // Furthermore, TObject::kInvalidObject and its semantic is not (and should not be)
      // used in TStreamerElement
      kDoNotDelete  = TStreamerElement::kDoNotDelete,

      // This bit duplicates TObject::kCannotPick. As the semantic of kHasRange is a persistent,
      // we can not change its value without breaking forward compatibility.
      // Furthermore, TObject::kCannotPick and its semantic is not (and should not be)
      // used in TStreamerElement
      kHasRange = TStreamerElement::kHasRange
   };
   ~~~ {.cpp}

  Without the EStatusBitsDupExceptions enum you would see

  ~~~ {.cpp}
TStatusBitsChecker::Check("TStreamerElement");

Error in <TStatusBitsChecker>: In TStreamerElement class hierarchy, there are duplicates bits:
Error in <TStatusBitsChecker>:    Bit   6 used in TStreamerElement as kHasRange
Error in <TStatusBitsChecker>:    Bit   6 used in TObject as kCannotPick
Error in <TStatusBitsChecker>:    Bit  13 used in TStreamerElement as kDoNotDelete
Error in <TStatusBitsChecker>:    Bit  13 used in TObject as kInvalidObject
  ~~~ {.cpp}

*/

#include "TStatusBitsChecker.h"

#include "TBaseClass.h"
#include "TClass.h"
#include "TClassTable.h"
#include "TEnum.h"
#include "TEnumConstant.h"
#include "TError.h"

#include <cmath>
#include <set>

namespace ROOT {
namespace Detail {

struct TStatusBitsChecker::Registry::Info {
   Info() = default;
   Info(const Info &) = default;
   Info(Info &&) = default;

   Info(TClass &o, std::string &&n, bool intentionalDup) : fOwner(&o), fConstantName(n), fIntentionalDup(intentionalDup)
   {
   }

   ~Info() = default;

   TClass *fOwner;
   std::string fConstantName;
   bool fIntentionalDup = false;
};

/// Default constructor. Implemented in source file to allow hiding of the Info struct.
TStatusBitsChecker::Registry::Registry() = default;

/// Default destructor. Implemented in source file to allow hiding of the Info struct.
TStatusBitsChecker::Registry::~Registry() = default;

/// Figure out which bit the constant has been set from/to.
/// Return 255 if the constant is not an integer or out of range.
UChar_t TStatusBitsChecker::ConvertToBit(Long64_t constant, TClass &classRef, const char *constantName)
{

   if (constant <= 0) {
      Error("TStatusBitsChecker::ConvertBit", "In %s the value of %s is %lld which was not produced by BIT macro.",
            classRef.GetName(), constantName, constant);
      return 255;
   }

   int backshift;
   double fraction = std::frexp(constant, &backshift);
   // frexp doc is:
   //    if no errors occur,
   //    returns the value x in the range (-1;-0.5], [0.5; 1)
   //    and stores an integer value in *exp such that x×2^(*exp) = arg
   --backshift; // frexp is such that BIT(0) == 1 == 0.5*2^(*exp) with *exp == 1

   if (backshift < 0 || std::abs(0.5 - fraction) > 0.00001f) {
      Error("TStatusBitsChecker::ConvertBit", "In %s the value of %s is %lld which was not produced by BIT macro.",
            classRef.GetName(), constantName, constant);
      return 255;
   }

   if (backshift > 24) {
      Error("TStatusBitsChecker::ConvertBit", "In %s the value of %s is %lld (>23) which is ignored by SetBit.",
            classRef.GetName(), constantName, constant);

      if (backshift > 255) // really we could snip it sooner.
         backshift = 255;
   }

   return backshift;
}

/// @brief Add to fRegister the Info about the bits in this class and its base
/// classes.
void TStatusBitsChecker::Registry::RegisterBits(TClass &classRef /* = false */)
{
   TEnum *eStatusBits = (TEnum *)classRef.GetListOfEnums()->FindObject("EStatusBits");
   TEnum *exceptionBits = (TEnum *)classRef.GetListOfEnums()->FindObject("EStatusBitsDupExceptions");

   if (eStatusBits) {

      for (auto constant : TRangeStaticCast<TEnumConstant>(*eStatusBits->GetConstants())) {

         // Ignore the known/intentional duplication.
         bool intentionalDup = exceptionBits && exceptionBits->GetConstant(constant->GetName());

         auto value = constant->GetValue();
         auto bit = ConvertToBit(value, classRef, constant->GetName());

         if (bit < 24) {
            bool need = true;
            for (auto reg : fRegister[bit]) {
               if (reg.fOwner == &classRef) {
                  // We have a duplicate declared in the same class
                  // let's accept this as an alias.
                  need = false;
               }
            }

            if (need)
               fRegister[bit].emplace_back(classRef, std::string(constant->GetName()), intentionalDup);
         }
      }
   }

   TList *lb = classRef.GetListOfBases();
   if (lb) {
      for (auto base : TRangeStaticCast<TBaseClass>(*lb)) {
         TClass *bcl = base->GetClassPointer();
         if (bcl)
            RegisterBits(*bcl);
      }
   }
}

/// @brief Return false and print error messages if there is any unexpected
/// duplicates BIT constant in the class hierarchy or any of the bits
/// already registered.
/// If verbose is true, also print all the bit declare in this class
/// and all its bases.
bool TStatusBitsChecker::Registry::Check(TClass &classRef, bool verbose /* = false */)
{
   RegisterBits(classRef);

   if (verbose) {
      for (auto cursor : fRegister) {
         for (auto constant : cursor.second) {
            Printf("Bit %3d declared in %s as %s", cursor.first, constant.fOwner->GetName(),
                   constant.fConstantName.c_str());
         }
      }
   }

   bool issuedHeader = false;
   bool result = true;
   for (auto cursor : fRegister) {
      unsigned int nDuplicate = 0;
      for (auto constant : cursor.second) {
         if (!constant.fIntentionalDup)
            ++nDuplicate;
      }
      if (nDuplicate > 1) {
         if (!issuedHeader) {
            Error("TStatusBitsChecker", "In %s class hierarchy, there are duplicates bits:", classRef.GetName());
            issuedHeader = true;
         }
         for (auto constant : cursor.second) {
            if (!constant.fIntentionalDup) {
               Error("TStatusBitsChecker", "   Bit %3d used in %s as %s", cursor.first, constant.fOwner->GetName(),
                     constant.fConstantName.c_str());
               result = false;
            }
         }
      }
   }

   return result;
}

/// Return false and print error messages if there is any unexpected
/// duplicates BIT constant in the class hierarchy.
/// If verbose is true, also print all the bit declare in this class
/// and all its bases.
bool TStatusBitsChecker::Check(TClass &classRef, bool verbose /* = false */)
{
   return Registry().Check(classRef, verbose);
}

/// Return false and print error messages if there is any unexpected
/// duplicates BIT constant in the class hierarchy.
/// If verbose is true, also print all the bit declare in this class
/// and all its bases.
bool TStatusBitsChecker::Check(const char *classname, bool verbose)
{
   TClass *cl = TClass::GetClass(classname);
   if (cl)
      return Check(*cl, verbose);
   else
      return true;
}

/// Return false and print error messages if there is any unexpected
/// duplicates BIT constant in any of the class hierarchy knows
/// to TClassTable.
/// If verbose is true, also print all the bit declare in eacho of the classes
/// and all their bases.
bool TStatusBitsChecker::CheckAllClasses(bool verbosity /* = false */)
{
   bool result = true;

   // Start from beginning
   gClassTable->Init();

   std::set<std::string> rootLibs;
   TList classesDeclFileNotFound;
   TList classesImplFileNotFound;

   Int_t totalNumberOfClasses = gClassTable->Classes();
   for (Int_t i = 0; i < totalNumberOfClasses; i++) {

      // get class name
      const char *cname = gClassTable->Next();
      if (!cname)
         continue;

      // get class & filename - use TROOT::GetClass, as we also
      // want those classes without decl file name!
      TClass *classPtr = TClass::GetClass((const char *)cname, kTRUE);
      if (!classPtr)
         continue;

      result = Check(*classPtr, verbosity) && result;
   }

   return result;
}

} // namespace Detail
} // namespace ROOT
