/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RLegend
#define ROOT7_RLegend

#include <ROOT/RPave.hxx>
#include <ROOT/RAttrLine.hxx>
#include <ROOT/RAttrFill.hxx>
#include <ROOT/RAttrMarker.hxx>
#include <ROOT/RDisplayItem.hxx>

#include <initializer_list>
#include <memory>
#include <vector>

namespace ROOT {
namespace Experimental {

/** \class RLegend
\ingroup GrafROOT7
\brief A legend for several drawables
\author Sergey Linev <S.Linev@gsi.de>
\date 2019-10-09
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RLegend : public RPave {

public:
   /** \class REntry
   \ingroup GrafROOT7
   \brief An entry in RLegend, references RDrawable and its attributes
   \author Sergey Linev <S.Linev@gsi.de>
   \date 2019-10-09
   \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
   is welcome!
   */

   class REntry {

      friend class RLegend;

      std::string fLabel; ///< label shown for the entry

      bool fLine{true}, fFill{false}, fMarker{false}; ///< enable line, fill, marker showing

      Internal::RIOShared<RDrawable> fDrawable; ///< reference to RDrawable

      std::string fDrawableId; ///< drawable id, used only when display item

      bool IsCustomDrawable() const { return fDrawableId == "custom"; }

      bool EnsureCustomDrawable()
      {
         if (!IsCustomDrawable())
            return false;

         if (!fDrawable)
            fDrawable = std::make_shared<RDrawable>("lentry");

         return true;
      }

   public:
      REntry() = default;

      /** Create entry without reference to existing drawable object, can assign attributes */
      REntry(const std::string &lbl)
      {
         fLabel = lbl;
         fDrawableId = "custom";
      }

      /** Create entry with reference to existing drawable object */
      REntry(std::shared_ptr<RDrawable> drawable, const std::string &lbl)
      {
         fDrawable = drawable;
         fLabel = lbl;
      }

      REntry &SetLabel(const std::string &lbl)
      {
         fLabel = lbl;
         return *this;
      }
      const std::string &GetLabel() const { return fLabel; }

      REntry &SetLine(bool on = true)
      {
         fLine = on;
         return *this;
      }
      bool GetLine() const { return fLine; }

      REntry &SetAttrLine(const RAttrLine &attr)
      {
         if (EnsureCustomDrawable()) {
            RAttrLine(fDrawable.get()) = attr;
            SetLine(true);
         }
         return *this;
      }

      RAttrLine GetAttrLine() const
      {
         if (IsCustomDrawable() && fDrawable)
            return RAttrLine(const_cast<RDrawable *>(fDrawable.get()));
         return {};
      }

      REntry &SetFill(bool on = true)
      {
         fFill = on;
         return *this;
      }
      bool GetFill() const { return fFill; }

      REntry &SetAttrFill(const RAttrFill &attr)
      {
         if (EnsureCustomDrawable()) {
            RAttrFill(fDrawable.get()) = attr;
            SetFill(true);
         }
         return *this;
      }

      RAttrFill GetAttrFill() const
      {
         if (IsCustomDrawable() && fDrawable)
            return RAttrFill(const_cast<RDrawable *>(fDrawable.get()));
         return {};
      }

      REntry &SetMarker(bool on = true)
      {
         fMarker = on;
         return *this;
      }
      bool GetMarker() const { return fMarker; }

      REntry &SetAttrMarker(const RAttrMarker &attr)
      {
         if (EnsureCustomDrawable()) {
            RAttrMarker(fDrawable.get()) = attr;
            SetMarker(true);
         }
         return *this;
      }

      RAttrMarker GetAttrMarker() const
      {
         if (IsCustomDrawable() && fDrawable)
            return RAttrMarker(const_cast<RDrawable *>(fDrawable.get()));
         return {};
      }
   };

private:
   std::string fTitle; ///< legend title

   std::vector<REntry> fEntries; ///< list of entries which should be displayed

protected:
   void CollectShared(Internal::RIOSharedVector_t &vect) override
   {
      for (auto &entry : fEntries) {
         vect.emplace_back(&entry.fDrawable);
         if (entry.fDrawable)
            entry.fDrawable->CollectShared(vect);
      }
   }

   /** hide I/O pointers when creating display item */
   std::unique_ptr<RDisplayItem> Display(const RDisplayContext &) override
   {
      for (auto &entry : fEntries) {
         if (!entry.IsCustomDrawable()) {
            entry.fDrawableId = RDisplayItem::ObjectIDFromPtr(entry.fDrawable.get());
            entry.fDrawable.reset_io();
         }
      }

      return std::make_unique<RDrawableDisplayItem>(*this);
   }

   /** when display item destroyed - restore I/O pointers */
   void OnDisplayItemDestroyed(RDisplayItem *) const override
   {
      for (auto &centry : fEntries) {
         if (!centry.IsCustomDrawable()) {
            auto entry = const_cast<REntry *>(&centry);
            entry->fDrawable.restore_io();
            entry->fDrawableId.clear();
         }
      }
   }

public:
   RLegend() : RPave("legend") {}

   RLegend(const std::string &title) : RLegend() { SetTitle(title); }

   RLegend &SetTitle(const std::string &title)
   {
      fTitle = title;
      return *this;
   }
   const std::string &GetTitle() const { return fTitle; }

   REntry &AddEntry(const std::string &lbl)
   {
      fEntries.emplace_back(lbl);
      return fEntries.back();
   }

   REntry &AddEntry(std::shared_ptr<RDrawable> drawable, const std::string &lbl)
   {
      fEntries.emplace_back(drawable, lbl);
      return fEntries.back();
   }

   auto NumEntries() const { return fEntries.size(); }

   auto &GetEntry(int n) { return fEntries[n]; }
};

} // namespace Experimental
} // namespace ROOT

#endif
