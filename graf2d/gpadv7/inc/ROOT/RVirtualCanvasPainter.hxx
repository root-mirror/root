/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RVirtualCanvasPainter
#define ROOT7_RVirtualCanvasPainter

#include <memory>
#include <functional>
#include <string>

namespace ROOT {
namespace Experimental {

using CanvasCallback_t = std::function<void(bool)>;

class RCanvas;
class RWebWindow;

namespace Internal {

/** \class ROOT::Experimental::Internal::RVirtualCanvasPainter
\ingroup GpadROOT7
\brief Abstract interface for painting a canvas.
\author Axel Naumann <axel@cern.ch>
\date 2017-05-31
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RVirtualCanvasPainter {
protected:
   class Generator {
   public:
      /// Abstract interface to create a RVirtualCanvasPainter implementation.
      virtual std::unique_ptr<RVirtualCanvasPainter> Create(RCanvas &canv) const = 0;
      /// Default destructor.
      virtual ~Generator();
   };

   /// generator getter
   static std::unique_ptr<Generator> &GetGenerator();

public:
   /// Default destructor.
   virtual ~RVirtualCanvasPainter();

   /// indicate that canvas changed, provides current version of the canvas
   virtual void CanvasUpdated(uint64_t, bool, CanvasCallback_t) = 0;

   /// return true if canvas modified since last painting
   virtual bool IsCanvasModified(uint64_t) const = 0;

   /// perform special action when drawing is ready
   virtual void DoWhenReady(const std::string &, const std::string &, bool, CanvasCallback_t) = 0;

   /// produce file output in batch mode like png, jpeg, svg or pdf
   virtual bool ProduceBatchOutput(const std::string &, int, int) = 0;

   virtual void NewDisplay(const std::string &where) = 0;

   virtual int NumDisplays() const = 0;

   virtual std::string GetWindowAddr() const = 0;

   /// run canvas functionality in caller thread, not needed when main thread is used
   virtual void Run(double tm = 0.) = 0;

   virtual bool AddPanel(std::shared_ptr<RWebWindow>) { return false; }

   /// Loads the plugin that implements this class.
   static std::unique_ptr<RVirtualCanvasPainter> Create(RCanvas &canv);
};
} // namespace Internal
} // namespace Experimental
} // namespace ROOT

#endif // ROOT7_RVirtualCanvasPainter
