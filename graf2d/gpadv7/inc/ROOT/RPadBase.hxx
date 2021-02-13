/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RPadBase
#define ROOT7_RPadBase

#include "ROOT/RDrawable.hxx"
#include "ROOT/RFrame.hxx"
#include "ROOT/RPadExtent.hxx"
#include "ROOT/RPadPos.hxx"
#include "ROOT/RPadUserAxis.hxx"
#include "ROOT/TypeTraits.hxx"

#include <memory>
#include <vector>

namespace ROOT {
namespace Experimental {

class RPad;
class RCanvas;
class RPadBaseDisplayItem;

/** \class ROOT::Experimental::RPadBase
\ingroup GpadROOT7
\brief Base class for graphic containers for `RDrawable`-s.
\authors Axel Naumann <axel@cern.ch> Sergey Linev <s.linev@gsi.de>
\date 2019-10-02
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RPadBase : public RDrawable {

private:

   using Primitive_t = Internal::RIOShared<RDrawable>;

   /// Content of the pad.

   std::vector<Primitive_t> fPrimitives;

   /// Disable copy construction.
   RPadBase(const RPadBase &) = delete;

   /// Disable assignment.
   RPadBase &operator=(const RPadBase &) = delete;

   void TestIfFrameRequired(const RDrawable *drawable)
   {
      if (drawable->IsFrameRequired())
         GetOrCreateFrame();
   }

protected:
   /// Allow derived classes to default construct a RPadBase.
   RPadBase() : RDrawable("pad") {}

   void CollectShared(Internal::RIOSharedVector_t &) override;

   void DisplayPrimitives(RPadBaseDisplayItem &paditem, RDisplayContext &ctxt);

   void SetDrawableVersion(Version_t vers) override;

public:

   using Primitives_t = std::vector<std::shared_ptr<RDrawable>>;

   virtual ~RPadBase();

   void UseStyle(const std::shared_ptr<RStyle> &style) override;

   /// Divide this pad into a grid of subpads with padding in between.
   /// \param nHoriz Number of horizontal pads.
   /// \param nVert Number of vertical pads.
   /// \param padding Padding between pads.
   /// \returns vector of vector (ret[x][y]) of created pads.
   std::vector<std::vector<std::shared_ptr<RPad>>> Divide(int nHoriz, int nVert, const RPadExtent &padding = {});

   /// Create drawable of specified class T
   template<class T, class... ARGS>
   auto Draw(ARGS... args)
   {
      auto drawable = std::make_shared<T>(args...);

      TestIfFrameRequired(drawable.get());

      fPrimitives.emplace_back(drawable);

      return drawable;
   }

   /// Add existing drawable instance to canvas
   auto Draw(std::shared_ptr<RDrawable> &&drawable)
   {
      TestIfFrameRequired(drawable.get());

      fPrimitives.emplace_back(std::move(drawable));

      return fPrimitives.back().get_shared();
   }

   /// Add object to be painted.
   /// Correspondent drawable will be created via GetDrawable() function which should be defined and be accessed at calling time.
   /// If required, extra arguments for GetDrawable() function can be provided.
   template <class T, class... ARGS>
   auto Draw(const std::shared_ptr<T> &what, ARGS... args)
   {
      // Requires GetDrawable(what) to be known!
      auto drawable = GetDrawable(what, args...);

      TestIfFrameRequired(drawable.get());

      fPrimitives.emplace_back(drawable);

      return drawable;
   }

   /// returns number of primitives in the pad
   unsigned NumPrimitives() const { return fPrimitives.size(); }

   /// returns primitive of given number
   std::shared_ptr<RDrawable> GetPrimitive(unsigned num) const
   {
      if (num >= fPrimitives.size()) return nullptr;
      return fPrimitives[num].get_shared();
   }

   std::shared_ptr<RDrawable> FindPrimitive(const std::string &id) const;

   std::shared_ptr<RDrawable> FindPrimitiveByDisplayId(const std::string &display_id) const;

   const RPadBase *FindPadForPrimitiveWithDisplayId(const std::string &display_id) const;

   /// Get all primitives contained in the pad.
   auto GetPrimitives() const
   {
      Primitives_t res;
      for (auto &entry : fPrimitives)
         res.emplace_back(entry.get_shared());
      return res;
   }

   /// Remove an object from the list of primitives.
   bool Remove(const std::string &id)
   {
      auto iter = std::find_if(fPrimitives.begin(), fPrimitives.end(),
         [&id](const Internal::RIOShared<RDrawable>& dr) { return dr->GetId() == id; });
      if (iter == fPrimitives.end())
         return false;
      iter->reset();
      fPrimitives.erase(iter);
      return true;
   }

   /// Remove drawable from list of primitives
   bool Remove(const std::shared_ptr<RDrawable> &drawable)
   {
      auto iter = std::find_if(fPrimitives.begin(), fPrimitives.end(),
         [&drawable](const Internal::RIOShared<RDrawable>& dr) { return drawable.get() == dr.get(); });
      if (iter == fPrimitives.end())
         return false;
      iter->reset();
      fPrimitives.erase(iter);
      return true;
   }

   /// Remove drawable at specified position
   bool RemoveAt(unsigned indx)
   {
      if (indx >= fPrimitives.size()) return false;
      fPrimitives[indx].reset();
      fPrimitives.erase(fPrimitives.begin() + indx);
      return true;
   }

   /// Wipe the pad by clearing the list of primitives.
   void Wipe() { fPrimitives.clear(); }

   std::shared_ptr<RFrame> GetOrCreateFrame();
   std::shared_ptr<RFrame> GetFrame();
   const std::shared_ptr<RFrame> GetFrame() const;

   RPadUserAxisBase* GetOrCreateAxis(size_t dimension);
   RPadUserAxisBase* GetAxis(size_t dimension) const;

   void SetAxisBounds(int dimension, double begin, double end);
   void SetAxisBound(int dimension, RPadUserAxisBase::EAxisBoundsKind boundsKind, double bound);
   void SetAxisAutoBounds(int dimension);

   void SetAllAxisBounds(const std::vector<std::array<double, 2>> &vecBeginAndEnd);

   /// Simple struct representing an axis bound.
   struct BoundKindAndValue {
      RPadUserAxisBase::EAxisBoundsKind fKind = RPadUserAxisBase::kAxisBoundsAuto;
      double fBound = 0.;
   };
   void SetAllAxisBound(const std::vector<BoundKindAndValue> &vecBoundAndKind);
   void SetAllAxisAutoBounds();

   /// Convert a `Pixel` position to Canvas-normalized positions.
   virtual std::array<RPadLength::Normal, 2> PixelsToNormal(const std::array<RPadLength::Pixel, 2> &pos) const = 0;

   /// Access to the top-most canvas, if any (const version).
   virtual const RCanvas *GetCanvas() const = 0;

   /// Access to the top-most canvas, if any (non-const version).
   virtual RCanvas *GetCanvas() = 0;

   /// Convert user coordinates to normal coordinates.
   std::array<RPadLength::Normal, 2> UserToNormal(const std::array<RPadLength::User, 2> &pos) const;

   void AssignAutoColors();
};

} // namespace Experimental
} // namespace ROOT

#endif
