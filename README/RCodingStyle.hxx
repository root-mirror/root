/// \file RCodingStyle.hxx
/// \author Axel Naumann <axel@cern.ch>
/// \author Other Author <other@author>
/// \date 2018-07-24 or \date July, 2018
// The above entries are mostly for giving some context.
// The "author" field gives a hint whom to contact in case of questions, also
// from within the team. The date shows whether this is ancient or only part
// of the latest release. It's the date of creation / last massive rewrite.

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// This file demonstrates the coding styles to be used for new ROOT files.
// See the accompanying source file for context.

// File names correspond to the name of the main class that they contain.
// Our users will expect that a file called "RCodingStyle.hxx" contains a class
// named `RCodingStyle`.

// All files start with a Doxygen header, followed by the copyright statement.
// The "-year" part is kept up to date by a script.

// The file must contain no trailing whitespace including no indentation-only lines.
// Except for extreme cases, line breaks should happen between 80 and 120 characters.

// The first non-comment, non-empty line must be the #include guard.
// Acceptable forms are:
//   ROOT_RCodingStyle, ROOT_RCodingStyle_hxx, ROOT_RCodingStyle_Hxx,
//   ROOT_RCODINGSTYLE, ROOT_RCODINGSTYLE_HXX
// The first is preferred.
#ifndef ROOT_RCodingStyle

// Especially within header files, using forward declarations instead of
// including the relevant headers wherever possible.
// EXCEPTION: standard library classes must be #included, not forward declared.
// All files must be stand-alone: all types must be available through a
// forward declaration or an #include; it is not acceptable to rely on an
// included header to itself include another header needed by the file.

// Include ROOT's headers first, in alphabetical order:
//   Rationale: make sure that ROOT headers are standalone.
// Do not use #include guards at the point of `#include`.
//   Rationale: compilers optimize this themselves, this is interfering with their optimization.
// Use <>, not "".
//   Rationale: "" looks at ./ first; as ROOT headers are included by user code
//   that adds context fragility.
#include <ROOT/TSeq.hxx>
// Include v6 ROOT headers last.
// You're welcome to state what's used from a header, i.e. why it's included.
#include <Rtypes.h> // for ULong64_t

// Now include other C++ headers hat are not from the stdlib:
// e.g. #include <cppyy/CPpYy.h>

// Include stdlib headers next, in alphabetical order.
#include <string>

// Include C headers last.
// Rationale: non-ROOT, non-stdlib headers often #define like mad. Reduce interference
// with ROOT or stdlib headers.
#include <curses.h>

// Say NO to
// - `using namespace` in at global scope in headers
// - variables, types, functions etc at global scope
// - preprocessor `#define`s except for #include guards (and are cases where
//   `#include`d headers require them)

namespace ROOT {
// All new classes are within the namespace ROOT.

// Forward declare before class definitions.
namespace Experimental {
class RAxisConfig;
}

///\class RExampleClass
/// Classes and structs start with 'R'. And they are documented.
class RExampleClass {
public:
   // FIRST section: public types.

   /// Nested types are documented, too. Do not use `typedef` but `using`.
   /// Type aliases end with "_t".
   using SeqULL_t = TSeq<ULong64_t>;

   ///\class Nested
   /// This is a nested class. Only use `struct`s for simple aggregates without
   /// methods and all public members.
   struct Nested {
      std::string fNestedName; ///< The nested name.
   };

private:
   // SECOND section: private types.

   // Any type used by any public / protected interface must also be public / protected!
   // This is not enforced by C++ but required by this coding style.
   using Nesteds_t = std::vector<Nested>;

   // THIRD section: data members.

   // Data member names start with `f`.
   // Documentation comments can be in front of data members or right behind them, see below.
   // For the latter case, comments must be column-aligned:

   int fIntMember = -2;                  ///< Use inline initialization at least for fundamental types.
   std::vector<Nested> fManyStrings;     ///<! This marks a member that does not get serialized (I/O comment).
   static constexpr const int kSeq = -7; ///< A static constexpr (const) variable starts with `k`.

   /// Don't explicitly default initialize (`{}`)
   std::vector<ROOT::Experimental::RAxisConfig> fAxes; //! For pre-comment cases, ROOT I/O directives go here.

   /// A static variable starts with `fg`. Avoid all non-constexpr static variables! Long data
   /// member documentation is just fine, but then put it in front of the data member.
   static std::string fgNameBadPleaseUseStaticFunctionStaticVariable;

protected:
   /// FOURTH section: non-public functions.

   /// Instead of static data members (whether private or public), use outlined static
   /// functions. See the implementation in the source file.
   static const std::string &AccessStaticVar();

public:
   // FIFTH section: methods, starting with constructors, then alphabetical order, possibly
   // grouped by functionality (\ingroup or \{ \}).

   /// Use `= default` inside the class whenever possible
   RExampleClass() = default;

   /// Classes with virtual functions must have a virtual destructor. If all other functions are
   /// pure virtual or inlined, the destructor should be outlined to pin the class's vtable to an
   /// object file. I.e. say `= default` in the source, not here.
   virtual ~RExampleClass();

   /// For trivial functions, use `const` and possibly `constexpr`, but not `noexcept`.
   int GetAsInt() const { return fIntMember; }

   /// Static functions line up alphabetically with the others.
   static int GetSeq() { return kSeq; }

   /// It is fine to have one-line function declarations.
   void SetInt(int i) { fIntMember = i; }
   
   /// \brief Replace the many strings
   /// \param[in] input A collection of strings to use.
   /// \return Always return true but if it could fail would return false in that case.
   ///
   /// All functions must have at least one line of comments, parseable/useable by doxygen.
   bool SetManyStrings(const std::vector<Nested>& input) { fManyStrings = input; return true; }

   /// Virtual functions are still fine (though we try to avoid them as they come with a call
   /// and optimization penalty). See `RCodingStyle` for overriding virtual functions.
   virtual void VirtualFunction() = 0;
};

///\class RCodingStyle
/// This is exhaustive documentation of the class template, explaining the constraints
/// on the template parameters, e.g.
///   `IDX` must be positive, `T` must be copyable.
/// Disabled if `T` is of reference type.
template <int IDX, // Template parameters are all-capital letters
          class T, // Use `typename` for T being a builtin (int, float,...) else `class`.
                   // template parameter name `T` is fine for generic class names
                   // we use `enable_if` through unnamed template parameters and provide a `assert`-style message
                   // to make the diagnostics more understandable.
          class = std::enable_if<"Must only be used with non-reference types" && !std::is_reference<T>::value>>
class RCodingStyle : public RExampleClass {
private:
   T fMember; ///< The wrapped value.

public:
   // Use static_assert excessively: it makes diagnostics much more readable,
   // and helps readers of your code.
   static_assert(IDX > 0, "IDX must be >0.");

   /// Defaulted default constructor.
   RCodingStyle() = default;

   // Open a group:
   /// \{ \name copy pperations

   // Use `= delete` instead of private, unimplemented.
   // To avoid copy and move construction and assignment, `delete` the copy
   // constructor and copy-assignment operator, see https://stackoverflow.com/a/15181645/6182509
   RCodingStyle(const RCodingStyle &) = delete;
   RCodingStyle &operator=(const RCodingStyle &) = delete;
   /// \}

   /// Virtual destructors of derived classes may be decorated with `virtual` (even though
   /// they "inherit" the virtuality of the virtual base class destructor).
   virtual ~RCodingStyle() {}

   /// Overridden virtual functions do not specify `virtual` but either `override` or `final`.
   void VirtualFunction() override {}

   /// "Free" (i.e. non-class member) operators should be declared as "hidden friends".
   /// This makes them available only when needed, instead of adding to the long list of
   /// "no matching overload for a+b, did you mean ...(hundreds of unrelated types)..."
   friend RCodingStyle operator+(const RCodingStyle &a, const RCodingStyle &b)
   {
      RCodingStyle ret{a};
      ret.SetInt(ret.GetAsInt() + b.GetAsInt());
      return ret;
   }

   /// Heterogenous operators should be implemented as hidden-friend functions:
   friend RCodingStyle operator+(const RCodingStyle &a, int i)
   {
      RCodingStyle ret{a};
      ret.SetInt(ret.GetAsInt() + i);
      return ret;
   }

   /// ...as the symmetrical / inverse version can only be specified as non-member:
   friend RCodingStyle operator+(int i, const RCodingStyle &a) { return a + i; }
};

/// Functions that access only public members of ROOT classes shall be in here.
/// \param a - first RExampleClass to add.
/// \param b - second RExampleClass to add.
double Add(const RExampleClass &a, const RExampleClass &b);

// Rules on sub-namespaces.
// All namespaces are CamelCase, i.e. starting with an upper case character.
namespace CodingStyle {
// Parts of ROOT might want to use sub-namespaces for *internal* code.
// Such code might include the following two namespaces:
namespace Detail {
// Contains code that is not expected to be seen by the everage user.
// It might provide customization points for advanced users, or
// interfaces that are useful for performance-critical code, at the
// expense of robustness.
}

namespace Internal {
// Contains code that is meant to be used only by ROOT itself. We do not guarantee
// interface stability for any type inside any sub-namespace of `ROOT::` called
// `Internal`, i.e. `ROOT::Internal` or `ROOT::...::Internal`.
}
}

} // namespace ROOT
#endif // ROOT_RCodingStyle
// The file must end on a trailing new-line.
