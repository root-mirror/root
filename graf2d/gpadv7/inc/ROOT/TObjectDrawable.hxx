/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_TObjectDrawable
#define ROOT7_TObjectDrawable

#include <ROOT/RDrawable.hxx>

class TObject;
class TColor;

namespace ROOT {
namespace Experimental {

class RPadBase;

/** \class TObjectDrawable
\ingroup GpadROOT7
\brief Provides v7 drawing facilities for TObject types (TGraph etc).
\author Sergey Linev
\date 2017-05-31
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class TObjectDrawable final : public RDrawable {
private:

   enum {
      kNone = 0,     ///< empty container
      kObject = 1,   ///< plain object
   };

   int fKind{kNone};                   ///< object kind
   Internal::RIOShared<TObject> fObj;  ///< The object to be painted
   std::string fOpts;                  ///< drawing options

   const char *GetColorCode(TColor *col);

   std::unique_ptr<TObject> CreateSpecials(int kind);

protected:

   void CollectShared(Internal::RIOSharedVector_t &vect) final { vect.emplace_back(&fObj); }

   std::unique_ptr<RDisplayItem> Display(const RDisplayContext &) override;

   void PopulateMenu(RMenuItems &) final;

   void Execute(const std::string &) final;

public:
   // special kinds, see TWebSnapshot enums
   enum EKind {
      kColors = 4,   ///< list of ROOT colors
      kStyle = 5,    ///< instance of TStyle object
      kPalette = 6   ///< list of colors from palette
   };

   TObjectDrawable() : RDrawable("tobject") {}

   TObjectDrawable(const std::shared_ptr<TObject> &obj, const std::string &opt = "") : RDrawable("tobject"), fKind(kObject), fObj(obj), fOpts(opt) {}

   TObjectDrawable(EKind kind, bool persistent = false);

   virtual ~TObjectDrawable() = default;
};

} // namespace Experimental
} // namespace ROOT


#endif
