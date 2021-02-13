/// \file ROOT/RHistDrawable.hxx
/// \ingroup HistDraw ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-07-09
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RHistDrawable
#define ROOT7_RHistDrawable

#include <ROOT/RDrawable.hxx>
#include <ROOT/RAttrLine.hxx>
#include <ROOT/RAttrText.hxx>
#include <ROOT/RAttrMarker.hxx>
#include <ROOT/RAttrFill.hxx>
#include <ROOT/RAttrValue.hxx>
#include <ROOT/RHist.hxx>
#include <ROOT/RHistImpl.hxx>

// TODO: move to separate file
#include <ROOT/RDrawableRequest.hxx>
#include <ROOT/RDisplayItem.hxx>

#include <memory>

namespace ROOT {
namespace Experimental {

class RHistDrawableBase : public RDrawable {
   RAttrValue<std::string>  fKind{this, "kind", ""};        ///<! hist draw kind
   RAttrValue<int>          fSub{this, "sub", -1};          ///<! hist draw sub kind
   RAttrLine                fAttrLine{this, "line"};        ///<! hist line attributes
   RAttrFill                fAttrFill{this, "fill"};        ///<! hist fill attributes
   RAttrText                fAttrText{this, "text"};        ///<! hist text attributes
   RAttrMarker              fMarkerAttr{this, "marker"};    ///<! hist marker attributes
   RAttrValue<bool>         fOptimize{this, "optimize", false}; ///<! optimize drawing

protected:

   bool IsFrameRequired() const final { return true; }

   void PopulateMenu(RMenuItems &) override { }

   void SetDrawKind(const std::string &kind, int sub = -1)
   {
      fKind = kind;
      if (sub >= 0)
         fSub = sub;
      else
         fSub.Clear();
   }

   std::string GetDrawKind() const { return fKind; }

   virtual std::unique_ptr<RDisplayItem> CreateHistDisplay(const RDisplayContext &) = 0;

   virtual bool Is3D() const { return false; }

   std::unique_ptr<RDisplayItem> Display(const RDisplayContext &ctxt) override
   {
      if (fOptimize)
         return CreateHistDisplay(ctxt);

      return RDrawable::Display(ctxt);
   }

public:

   class RReply : public RDrawableReply {
   public:
      std::unique_ptr<RDisplayItem> item;
   };

   class RRequest : public RDrawableRequest {
      std::unique_ptr<RDrawableReply> Process() override
      {
         auto hdraw = dynamic_cast<RHistDrawableBase *>(GetContext().GetDrawable());

         auto reply = std::make_unique<RReply>();
         if (hdraw)
            reply->item = hdraw->CreateHistDisplay(GetContext());
         return reply;
      }
   };

   friend class RRequest;

   RHistDrawableBase() : RDrawable("hist") {}

   const RAttrLine &GetAttrLine() const { return fAttrLine; }
   RHistDrawableBase &SetAttrLine(const RAttrLine &attr) { fAttrLine = attr; return *this; }
   RAttrLine &AttrLine() { return fAttrLine; }

   const RAttrFill &GetAttrFill() const { return fAttrFill; }
   RHistDrawableBase &SetAttrFill(const RAttrFill &fill) { fAttrFill = fill; return *this; }
   RAttrFill &AttrFill() { return fAttrFill; }

   const RAttrText &GetAttrText() const { return fAttrText; }
   RHistDrawableBase &SetAttrText(const RAttrText &attr) { fAttrText = attr; return *this; }
   RAttrText &AttrText() { return fAttrText; }

   const RAttrMarker &GetAttrMarker() const { return fMarkerAttr; }
   RHistDrawableBase &SetAttrMarker(const RAttrMarker &attr) { fMarkerAttr = attr; return *this; }
   RAttrMarker &AttrMarker() { return fMarkerAttr; }

   RHistDrawableBase &Optimize(bool on = true) { fOptimize = on; return *this; }
};


template <int DIMENSIONS>
class RHistDrawable : public RHistDrawableBase {
public:
   using HistImpl_t = Detail::RHistImplPrecisionAgnosticBase<DIMENSIONS>;

protected:

   Internal::RIOShared<HistImpl_t> fHistImpl;             ///< I/O capable reference on histogram

   void CollectShared(Internal::RIOSharedVector_t &vect) override { vect.emplace_back(&fHistImpl); }

public:
   RHistDrawable() = default;
   virtual ~RHistDrawable() = default;

   template <class HIST>
   RHistDrawable(const std::shared_ptr<HIST> &hist) : RHistDrawableBase()
   {
      fHistImpl = std::shared_ptr<HistImpl_t>(hist, hist->GetImpl());
   }

   std::shared_ptr<HistImpl_t> GetHist() const { return fHistImpl.get_shared(); }
};


class RHist1Drawable final : public RHistDrawable<1> {
   RAttrValue<double> fBarOffset{this, "bar_offset", 0.}; ///<!  bar offset
   RAttrValue<double> fBarWidth{this, "bar_width", 1.};   ///<!  bar width
   RAttrValue<bool> fText{this, "text", false};           ///<! draw text

protected:
   std::unique_ptr<RDisplayItem> CreateHistDisplay(const RDisplayContext &) override;

   bool Is3D() const final { return GetDrawKind() == "lego"; }

public:
   RHist1Drawable() = default;

   template <class HIST>
   RHist1Drawable(const std::shared_ptr<HIST> &hist) : RHistDrawable<1>(hist) {}

   RHist1Drawable &Bar() { SetDrawKind("bar", 0); fBarOffset.Clear(); fBarWidth.Clear(); return *this; }
   RHist1Drawable &Bar(double offset, double width, bool mode3d = false) { SetDrawKind("bar", mode3d ? 1 : 0); fBarOffset = offset; fBarWidth = width; return *this; }
   RHist1Drawable &Error(int kind = 0) { SetDrawKind("err", kind); return *this; }
   RHist1Drawable &Marker() { SetDrawKind("p"); return *this; }
   RHist1Drawable &Star() { AttrMarker().SetStyle(3); return Marker(); }
   RHist1Drawable &Hist() { SetDrawKind("hist"); return *this; }
   RHist1Drawable &Line() { SetDrawKind("l"); return *this; }
   RHist1Drawable &Lego(int kind = 0) { SetDrawKind("lego", kind); return *this; }
   RHist1Drawable &Text(bool on = true) { fText = on; return *this; }

   double GetBarOffset() const { return fBarOffset; }
   double GetBarWidth() const { return fBarWidth; }
};


class RHist2Drawable final : public RHistDrawable<2> {
   RAttrValue<bool> fText{this, "text", false};               ///<! draw text

protected:

   std::unique_ptr<RDisplayItem> CreateHistDisplay(const RDisplayContext &) override;

   bool Is3D() const final { return (GetDrawKind() == "lego") || (GetDrawKind() == "surf") || (GetDrawKind() == "err"); }

public:
   RHist2Drawable() = default;

   template <class HIST>
   RHist2Drawable(const std::shared_ptr<HIST> &hist) : RHistDrawable<2>(hist) {}

   RHist2Drawable &Color() { SetDrawKind("col"); return *this; }
   RHist2Drawable &Box(int kind = 0) { SetDrawKind("box", kind); return *this; }
   RHist2Drawable &Lego(int kind = 0) { SetDrawKind("lego", kind); return *this; }
   RHist2Drawable &Surf(int kind = 0) { SetDrawKind("surf", kind); return *this; }
   RHist2Drawable &Error() { SetDrawKind("err"); return *this; }
   RHist2Drawable &Contour(int kind = 0) { SetDrawKind("cont", kind); return *this; }
   RHist2Drawable &Scatter() { SetDrawKind("scat"); return *this; }
   RHist2Drawable &Arrow() { SetDrawKind("arr"); return *this; }
   RHist2Drawable &Text(bool on = true) { fText = on; return *this; }
};


class RHist3Drawable final : public RHistDrawable<3> {
protected:
   std::unique_ptr<RDisplayItem> CreateHistDisplay(const RDisplayContext &) override;

   bool Is3D() const final { return true; }

public:
   RHist3Drawable() = default;

   template <class HIST>
   RHist3Drawable(const std::shared_ptr<HIST> &hist) : RHistDrawable<3>(hist) {}

   RHist3Drawable &Color() { SetDrawKind("col"); return *this; }
   RHist3Drawable &Box(int kind = 0) { SetDrawKind("box", kind); return *this; }
   RHist3Drawable &Sphere(int kind = 0) { SetDrawKind("sphere", kind); return *this; }
   RHist3Drawable &Scatter() { SetDrawKind("scat"); return *this; }
};


inline auto GetDrawable(const std::shared_ptr<RH1D> &histimpl)
{
   return std::make_shared<RHist1Drawable>(histimpl);
}

inline auto GetDrawable(const std::shared_ptr<RH1I> &histimpl)
{
   return std::make_shared<RHist1Drawable>(histimpl);
}

inline auto GetDrawable(const std::shared_ptr<RH1C> &histimpl)
{
   return std::make_shared<RHist1Drawable>(histimpl);
}

inline auto GetDrawable(const std::shared_ptr<RH1F> &histimpl)
{
   return std::make_shared<RHist1Drawable>(histimpl);
}

inline auto GetDrawable(const std::shared_ptr<RH2D> &histimpl)
{
   return std::make_shared<RHist2Drawable>(histimpl);
}

inline auto GetDrawable(const std::shared_ptr<RH2I> &histimpl)
{
   return std::make_shared<RHist2Drawable>(histimpl);
}

inline auto GetDrawable(const std::shared_ptr<RH2C> &histimpl)
{
   return std::make_shared<RHist2Drawable>(histimpl);
}

inline auto GetDrawable(const std::shared_ptr<RH2F> &histimpl)
{
   return std::make_shared<RHist2Drawable>(histimpl);
}

inline auto GetDrawable(const std::shared_ptr<RH3D> &histimpl)
{
   return std::make_shared<RHist3Drawable>(histimpl);
}

inline auto GetDrawable(const std::shared_ptr<RH3I> &histimpl)
{
   return std::make_shared<RHist3Drawable>(histimpl);
}

inline auto GetDrawable(const std::shared_ptr<RH3C> &histimpl)
{
   return std::make_shared<RHist3Drawable>(histimpl);
}

inline auto GetDrawable(const std::shared_ptr<RH3F> &histimpl)
{
   return std::make_shared<RHist3Drawable>(histimpl);
}

} // namespace Experimental
} // namespace ROOT

#endif
