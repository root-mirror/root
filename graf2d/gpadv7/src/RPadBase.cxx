/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RPadBase.hxx"

#include "ROOT/RLogger.hxx"
#include "ROOT/RPadExtent.hxx"
#include "ROOT/RPadPos.hxx"
#include <ROOT/RPad.hxx>
#include <ROOT/RCanvas.hxx>
#include <ROOT/RPadDisplayItem.hxx>

#include <cassert>
#include <limits>

using namespace std::string_literals;

ROOT::Experimental::RPadBase::~RPadBase() = default;

///////////////////////////////////////////////////////////////////////////
/// Use provided style for pad and all primitives inside

void ROOT::Experimental::RPadBase::UseStyle(const std::shared_ptr<RStyle> &style)
{
   RDrawable::UseStyle(style);
   for (auto &drawable : fPrimitives)
      drawable->UseStyle(style);
}

///////////////////////////////////////////////////////////////////////////
/// Find primitive with specified id

std::shared_ptr<ROOT::Experimental::RDrawable> ROOT::Experimental::RPadBase::FindPrimitive(const std::string &id) const
{
   for (auto &drawable : fPrimitives) {

      if (drawable->GetId() == id)
         return drawable.get_shared();

      const RPadBase *pad_draw = dynamic_cast<const RPadBase *> (drawable.get());

      if (pad_draw) {
         auto subelem = pad_draw->FindPrimitive(id);

         if (subelem)
            return subelem;
      }
   }

   return nullptr;
}

///////////////////////////////////////////////////////////////////////////
/// Find primitive with unique id, produce for RDisplayItem
/// Such id used for client-server identification of objects

std::shared_ptr<ROOT::Experimental::RDrawable> ROOT::Experimental::RPadBase::FindPrimitiveByDisplayId(const std::string &id) const
{
   auto p = id.find("_");
   if (p == std::string::npos)
      return nullptr;

   auto prim = GetPrimitive(std::stoul(id.substr(0,p)));
   if (!prim)
      return nullptr;

   auto subid = id.substr(p+1);

   if (RDisplayItem::ObjectIDFromPtr(prim.get()) == subid)
      return prim;

   auto subpad = std::dynamic_pointer_cast<RPadBase>(prim);

   return subpad ? subpad->FindPrimitiveByDisplayId(subid) : nullptr;
}

///////////////////////////////////////////////////////////////////////////
/// Find subpad which contains primitive with given display id

const ROOT::Experimental::RPadBase *ROOT::Experimental::RPadBase::FindPadForPrimitiveWithDisplayId(const std::string &id) const
{
   auto p = id.find("_");
   if (p == std::string::npos)
      return nullptr;

   auto prim = GetPrimitive(std::stoul(id.substr(0,p)));
   if (!prim)
      return nullptr;

   auto subid = id.substr(p+1);

   if (RDisplayItem::ObjectIDFromPtr(prim.get()) == subid)
      return this;

   auto subpad = std::dynamic_pointer_cast<RPadBase>(prim);

   return subpad ? subpad->FindPadForPrimitiveWithDisplayId(subid) : nullptr;
}

///////////////////////////////////////////////////////////////////////////
/// Method collect existing colors and assign new values if required

void ROOT::Experimental::RPadBase::AssignAutoColors()
{
   int cnt = 0;
   RColor col;

   for (auto &drawable : fPrimitives) {
      for (auto &attr: drawable->fAttr) {
         // only boolean attribute can return true
         if (!attr.second->GetBool()) continue;
         auto pos = attr.first.rfind("_color_auto");
         if ((pos > 0) && (pos == attr.first.length() - 11)) {
            // FIXME: dummy code to assign autocolors, later should use RPalette
            switch (cnt++ % 3) {
              case 0: col = RColor::kRed; break;
              case 1: col = RColor::kGreen; break;
              case 2: col = RColor::kBlue; break;
            }
            drawable->fAttr.AddString(attr.first.substr(0,pos) + "_color_rgb", col.AsHex());
         }
      }
   }
}

///////////////////////////////////////////////////////////////////////////
/// Create display items for all primitives in the pad
/// Each display item gets its special id, which used later for client-server communication
/// Second parameter is version id which already delivered to the client

void ROOT::Experimental::RPadBase::DisplayPrimitives(RPadBaseDisplayItem &paditem, RDisplayContext &ctxt)
{
   paditem.SetAttributes(&GetAttrMap());
   paditem.SetPadStyle(fStyle.lock());

   unsigned indx = 0;

   for (auto &drawable : fPrimitives) {

      ctxt.SetDrawable(drawable.get(), indx++);

      auto item = drawable->Display(ctxt);

      if (!item)
         item = std::make_unique<RDisplayItem>(true);

      item->SetObjectIDAsPtr(drawable.get());
      item->SetIndex(ctxt.GetIndex());
      // add object with the style
      paditem.Add(std::move(item), drawable->fStyle.lock());
   }
}

///////////////////////////////////////////////////////////////////////////
/// Divide pad on nHoriz X nVert subpads
/// Return array of array of pads

std::vector<std::vector<std::shared_ptr<ROOT::Experimental::RPad>>>
ROOT::Experimental::RPadBase::Divide(int nHoriz, int nVert, const RPadExtent &padding)
{
   std::vector<std::vector<std::shared_ptr<RPad>>> ret;
   if (!nHoriz)
      R__LOG_ERROR(GPadLog()) << "Cannot divide into 0 horizontal sub-pads!";
   if (!nVert)
      R__LOG_ERROR(GPadLog()) << "Cannot divide into 0 vertical sub-pads!";
   if (!nHoriz || !nVert)
      return ret;

   // Start with the whole (sub-)pad:
   RPadExtent offset{1._normal, 1._normal};
   /// We need n Pads plus n-1 padding. Thus each `(subPadSize + padding)` is `(parentPadSize + padding) / n`.
   offset = (offset + padding);
   offset *= {1. / nHoriz, 1. / nVert};
   const RPadExtent size = offset - padding;

   for (int iHoriz = 0; iHoriz < nHoriz; ++iHoriz) {
      ret.emplace_back();
      for (int iVert = 0; iVert < nVert; ++iVert) {
         RPadPos subPos = offset;
         subPos *= {1. * iHoriz, 1. * iVert};

         auto subpad = Draw<RPad>(this, subPos, size);

         ret.back().emplace_back(subpad);
         // printf("Create subpad pos %5.2f %5.2f\n", subPos.fHoriz.fNormal.fVal, subPos.fVert.fNormal.fVal);
      }
   }
   return ret;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Get a frame object for the pad.
/// If frame not exists - creates and add to the end of primitives list


std::shared_ptr<ROOT::Experimental::RFrame> ROOT::Experimental::RPadBase::GetOrCreateFrame()
{
   auto frame = GetFrame();
   if (!frame) {
      frame.reset(new RFrame());
      fPrimitives.emplace_back(frame);
   }

   return frame;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Get a frame object if exists

const std::shared_ptr<ROOT::Experimental::RFrame> ROOT::Experimental::RPadBase::GetFrame() const
{
   for (auto &drawable : fPrimitives) {
      if (drawable->GetCssType() == "frame") {
         const std::shared_ptr<RFrame> frame = std::dynamic_pointer_cast<RFrame>(drawable.get_shared());
         if (frame) return frame;
      }
   }
   return nullptr;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Get a frame object if exists

std::shared_ptr<ROOT::Experimental::RFrame> ROOT::Experimental::RPadBase::GetFrame()
{
   for (auto &drawable : fPrimitives) {
      if (drawable->GetCssType() == "frame") {
         std::shared_ptr<RFrame> frame = std::dynamic_pointer_cast<RFrame>(drawable.get_shared());
         if (frame) return frame;
      }
   }
   return nullptr;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Get a pad axis from the RFrame.
/// \param dimension - Index of the dimension of the RFrame user coordinate system.

ROOT::Experimental::RPadUserAxisBase* ROOT::Experimental::RPadBase::GetAxis(size_t dimension) const
{
   auto frame = GetFrame();

   if (frame && dimension < frame->GetNDimensions())
      return &frame->GetUserAxis(dimension);
   return nullptr;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Get a pad axis from the RFrame.
/// \param dimension - Index of the dimension of the RFrame user coordinate system.

ROOT::Experimental::RPadUserAxisBase* ROOT::Experimental::RPadBase::GetOrCreateAxis(size_t dimension)
{
   auto frame = GetOrCreateFrame();
   frame->GrowToDimensions(dimension);
   return &frame->GetUserAxis(dimension);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Set the range of an axis as begin, end.

void ROOT::Experimental::RPadBase::SetAxisBounds(int dimension, double begin, double end)
{
   GetOrCreateFrame()->GrowToDimensions(dimension);
   GetAxis(dimension)->SetBounds(begin, end);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Set the range of an axis as bound kind and bound (up or down).

void ROOT::Experimental::RPadBase::SetAxisBound(int dimension, RPadUserAxisBase::EAxisBoundsKind boundsKind, double bound)
{
   GetOrCreateFrame()->GrowToDimensions(dimension);
   GetAxis(dimension)->SetBound(boundsKind, bound);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Set the range of an axis as bound kind and bound (up or down).

void ROOT::Experimental::RPadBase::SetAxisAutoBounds(int dimension)
{
   GetOrCreateFrame()->GrowToDimensions(dimension);
   GetAxis(dimension)->SetAutoBounds();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Set the range of an axis as bound kind and bound (up or down).

void ROOT::Experimental::RPadBase::SetAllAxisBounds(const std::vector<std::array<double, 2>> &vecBeginAndEnd)
{
   auto frame = GetOrCreateFrame();

   frame->GrowToDimensions(vecBeginAndEnd.size());
   if (vecBeginAndEnd.size() != frame->GetNDimensions()) {
      R__LOG_ERROR(GPadLog())
         << "Array of axis bound has wrong size " <<  vecBeginAndEnd.size()
         << " versus numer of axes in frame " << frame->GetNDimensions();
      return;
   }

   for (size_t i = 0, n = frame->GetNDimensions(); i < n; ++i)
      frame->GetUserAxis(i).SetBounds(vecBeginAndEnd[i][0], vecBeginAndEnd[i][1]);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Set the range of an axis as bound kind and bound (up or down).

void ROOT::Experimental::RPadBase::SetAllAxisBound(const std::vector<BoundKindAndValue> &vecBoundAndKind)
{
   auto frame = GetOrCreateFrame();

   frame->GrowToDimensions(vecBoundAndKind.size());
   if (vecBoundAndKind.size() != frame->GetNDimensions()) {
      R__LOG_ERROR(GPadLog())
         << "Array of axis bound has wrong size " << vecBoundAndKind.size()
         << " versus numer of axes in frame " << frame->GetNDimensions();
      return;
   }

   for (size_t i = 0, n = frame->GetNDimensions(); i < n; ++i)
      frame->GetUserAxis(i).SetBound(vecBoundAndKind[i].fKind, vecBoundAndKind[i].fBound);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Collect all shared items to resolve shared_ptr after IO

void ROOT::Experimental::RPadBase::CollectShared(Internal::RIOSharedVector_t &vect)
{
   for (auto &handle : fPrimitives) {
      vect.emplace_back(&handle);
      auto drawable = handle.get();
      if (drawable) drawable->CollectShared(vect);
   }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Set the range of an axis as bound kind and bound (up or down).

void ROOT::Experimental::RPadBase::SetAllAxisAutoBounds()
{
   auto frame = GetOrCreateFrame();

   for (size_t i = 0, n = frame->GetNDimensions(); i < n; ++i)
      frame->GetUserAxis(i).SetAutoBounds();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Convert user coordinates to normal coordinates.

std::array<ROOT::Experimental::RPadLength::Normal, 2> ROOT::Experimental::RPadBase::UserToNormal(const std::array<RPadLength::User, 2> &pos) const
{
   auto frame = GetFrame();
   if (!frame) return {};

   return frame->UserToNormal(pos);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Assign drawable version - for pad itself and all primitives

void ROOT::Experimental::RPadBase::SetDrawableVersion(Version_t vers)
{
   RDrawable::SetDrawableVersion(vers);

   for (auto &drawable : fPrimitives)
      drawable->SetDrawableVersion(vers);
}
