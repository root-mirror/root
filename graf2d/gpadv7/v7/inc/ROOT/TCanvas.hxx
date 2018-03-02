/// \file ROOT/TCanvas.hxx
/// \ingroup Gpad ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-07-08
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_TCanvas
#define ROOT7_TCanvas

#include "ROOT/TPad.hxx"
#include "ROOT/TVirtualCanvasPainter.hxx"

#include <memory>
#include <string>
#include <vector>

namespace ROOT {
namespace Experimental {

/** \class ROOT::Experimental::TCanvas
  A window's topmost `TPad`.
  */

class TCanvas: public TPadBase {
private:
   /// Title of the canvas.
   std::string fTitle;

   /// Size of the canvas in pixels,
   std::array<TPadLength::Pixel, 2> fSize;

   /// Modify counter, incremented every time canvas is changed
   uint64_t fModified; ///<!

   /// The painter of this canvas, bootstrapping the graphics connection.
   /// Unmapped canvases (those that never had `Draw()` invoked) might not have
   /// a painter.
   std::unique_ptr<Internal::TVirtualCanvasPainter> fPainter; ///<!

   /// Disable copy construction for now.
   TCanvas(const TCanvas &) = delete;

   /// Disable assignment for now.
   TCanvas &operator=(const TCanvas &) = delete;

public:
   static std::shared_ptr<TCanvas> Create(const std::string &title);

   /// Create a temporary TCanvas; for long-lived ones please use Create().
   TCanvas() = default;

   ~TCanvas() = default;

   const TCanvas *GetCanvas() const override { return this; }

   /// Access to the top-most canvas, if any (non-const version).
   TCanvas *GetCanvas() override { return this; }

   /// Return canvas pixel size as array with two elements - width and height
   const std::array<TPadLength::Pixel, 2> &GetSize() const { return fSize; }

   /// Set canvas pixel size as array with two elements - width and height
   void SetSize(const std::array<TPadLength::Pixel, 2> &sz) { fSize = sz; }

   /// Set canvas pixel size - width and height
   void SetSize(const TPadLength::Pixel &width, const TPadLength::Pixel &height)
   {
      fSize[0] = width;
      fSize[1] = height;
   }

   /// Display the canvas.
   void Show(const std::string &where = "");

   /// Close all canvas displays
   void Hide();

   /// Insert panel into the canvas, canvas should be shown at this moment
   template <class PANEL>
   bool AddPanel(std::shared_ptr<PANEL> &panel) {
      if (!fPainter) return false;
      return fPainter->AddPanel(panel->GetWindow());
   }

   // Indicates that primitives list was changed or any primitive was modified
   void Modified() { fModified++; }

   // Return if canvas was modified and not yet updated
   bool IsModified() const;

   /// update drawing
   void Update(bool async = false, CanvasCallback_t callback = nullptr);

   /// Save canvas in image file
   void SaveAs(const std::string &filename, bool async = false, CanvasCallback_t callback = nullptr);

   /// Get the canvas's title.
   const std::string &GetTitle() const { return fTitle; }

   /// Set the canvas's title.
   void SetTitle(const std::string &title) { fTitle = title; }

   /// Convert a `Pixel` position to Canvas-normalized positions.
   std::array<TPadLength::Normal, 2> PixelsToNormal(const std::array<TPadLength::Pixel, 2> &pos) const final
   {
      return {{pos[0] / fSize[0], pos[1] / fSize[1]}};
   }

   static const std::vector<std::shared_ptr<TCanvas>> &GetCanvases();
};

} // namespace Experimental
} // namespace ROOT

#endif
