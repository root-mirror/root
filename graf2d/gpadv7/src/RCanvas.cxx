/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RCanvas.hxx"

#include "ROOT/RLogger.hxx"

#include <algorithm>
#include <memory>
#include <mutex>
#include <thread>
#include <chrono>
#include <stdio.h>
#include <string.h>

#include "TList.h"
#include "TROOT.h"

namespace {

static std::mutex &GetHeldCanvasesMutex()
{
   static std::mutex sMutex;
   return sMutex;
}

static std::vector<std::shared_ptr<ROOT::Experimental::RCanvas>> &GetHeldCanvases()
{
   static std::vector<std::shared_ptr<ROOT::Experimental::RCanvas>> sCanvases;
   return sCanvases;
}


} // namespace

///////////////////////////////////////////////////////////////////////////////////////
/// Returns list of created canvases

const std::vector<std::shared_ptr<ROOT::Experimental::RCanvas>> ROOT::Experimental::RCanvas::GetCanvases()
{
   std::lock_guard<std::mutex> grd(GetHeldCanvasesMutex());

   return GetHeldCanvases();
}

///////////////////////////////////////////////////////////////////////////////////////
/// Release list of held canvases pointers
/// If no other shared pointers exists on the canvas, object will be destroyed

void ROOT::Experimental::RCanvas::ReleaseHeldCanvases()
{
   std::vector<std::shared_ptr<ROOT::Experimental::RCanvas>> vect;

   {
      std::lock_guard<std::mutex> grd(GetHeldCanvasesMutex());

      std::swap(vect, GetHeldCanvases());
   }
}

///////////////////////////////////////////////////////////////////////////////////////
/// Returns true is canvas was modified since last painting

bool ROOT::Experimental::RCanvas::IsModified() const
{
   return fPainter ? fPainter->IsCanvasModified(fModified) : fModified;
}

///////////////////////////////////////////////////////////////////////////////////////
/// Update canvas

void ROOT::Experimental::RCanvas::Update(bool async, CanvasCallback_t callback)
{
   if (fPainter)
      fPainter->CanvasUpdated(fModified, async, callback);
}

class RCanvasCleanup : public TObject {
public:

   static RCanvasCleanup *gInstance;

   RCanvasCleanup() : TObject() { gInstance = this; }

   virtual ~RCanvasCleanup()
   {
      gInstance = nullptr;
      ROOT::Experimental::RCanvas::ReleaseHeldCanvases();
   }
};

RCanvasCleanup *RCanvasCleanup::gInstance = nullptr;


///////////////////////////////////////////////////////////////////////////////////////
/// Create new canvas instance

std::shared_ptr<ROOT::Experimental::RCanvas> ROOT::Experimental::RCanvas::Create(const std::string &title)
{
   auto pCanvas = std::make_shared<RCanvas>();
   pCanvas->SetTitle(title);
   {
      std::lock_guard<std::mutex> grd(GetHeldCanvasesMutex());
      GetHeldCanvases().emplace_back(pCanvas);
   }

   if (!RCanvasCleanup::gInstance) {
      auto cleanup = new RCanvasCleanup();
      TDirectory *dummydir = new TDirectory("rcanvas_cleanup_dummydir","title");
      dummydir->GetList()->Add(cleanup);
      gROOT->GetListOfClosedObjects()->Add(dummydir);
   }

   return pCanvas;
}

//////////////////////////////////////////////////////////////////////////
/// Create new display for the canvas
/// The parameter `where` specifies which program could be used for display creation
/// Possible values:
///
///  - `cef` Chromium Embeded Framework, local display, local communication
///  - `qt5` Qt5 WebEngine (when running via rootqt5), local display, local communication
///  - `browser` default system web-browser, communication via random http port from range 8800 - 9800
///  - `<prog>` any program name which will be started instead of default browser, like firefox or /usr/bin/opera
///     one could also specify $url in program name, which will be replaced with canvas URL
///  - `native` either any available local display or default browser
///
///  Canvas can be displayed in several different places

void ROOT::Experimental::RCanvas::Show(const std::string &where)
{
   if (fPainter) {
      bool isany = (fPainter->NumDisplays() > 0);

      if (!where.empty())
         fPainter->NewDisplay(where);

      if (isany) return;
   }

   if (!fModified)
      fModified = 1; // 0 is special value, means no changes and no drawings

   if (!fPainter)
      fPainter = Internal::RVirtualCanvasPainter::Create(*this);

   if (fPainter) {
      fPainter->NewDisplay(where);
      fPainter->CanvasUpdated(fModified, true, nullptr); // trigger async display
   }
}

//////////////////////////////////////////////////////////////////////////
/// Returns window name for canvas

std::string ROOT::Experimental::RCanvas::GetWindowAddr() const
{
   if (fPainter)
      return fPainter->GetWindowAddr();

   return "";
}


//////////////////////////////////////////////////////////////////////////
/// Hide all canvas displays

void ROOT::Experimental::RCanvas::Hide()
{
   if (fPainter)
      delete fPainter.release();
}

//////////////////////////////////////////////////////////////////////////
/// Create image file for the canvas
/// Supported SVG (extension .svg), JPEG (extension .jpg or .jpeg) and PNG (extension .png)

bool ROOT::Experimental::RCanvas::SaveAs(const std::string &filename)
{
   if (!fPainter)
      fPainter = Internal::RVirtualCanvasPainter::Create(*this);

   if (!fPainter)
      return false;

   auto width = fSize[0].fVal;
   auto height = fSize[1].fVal;

   return fPainter->ProduceBatchOutput(filename, width > 1 ? (int) width : 800, height > 1 ? (int) height : 600);
}

//////////////////////////////////////////////////////////////////////////
/// Remove canvas from global canvas lists, will be destroyed once last shared_ptr is disappear

void ROOT::Experimental::RCanvas::Remove()
{
   std::lock_guard<std::mutex> grd(GetHeldCanvasesMutex());
   auto &held = GetHeldCanvases();
   auto indx = held.size();
   while (indx-- > 0) {
      if (held[indx].get() == this)
         held.erase(held.begin() + indx);
   }
}

//////////////////////////////////////////////////////////////////////////
/// Run canvas functionality for the given time (in seconds)
/// Used to process canvas-related actions in the appropriate thread context.
/// Must be regularly called when canvas created and used in extra thread.
/// Time parameter specifies minimal execution time in seconds - if default value 0 is used,
/// just all pending actions will be performed.
/// When canvas is not yet displayed - just performs sleep for given time interval.
///
/// Example of usage:
///
/// ~~~ {.cpp}
/// void draw_canvas(bool &run_loop, std::make_shared<RH1D> hist)
/// {
///   auto canvas = RCanvas::Create("Canvas title");
///   canvas->Draw(hist)->SetLineColor(RColor::kBlue);
///   canvas->Show();
///   while (run_loop) {
///      pHist->Fill(1);
///      canvas->Modified();
///      canvas->Update();
///      canvas->Run(0.1); // process canvas events
///   }
///
///   canvas->Remove();
/// }
///
/// int main()
/// {
///    RAxisConfig xaxis(100, -10., 10.);
///    auto pHist = std::make_shared<RH1D>(xaxis);
///    bool run_loop = true;
///
///    std::thread thrd(draw_canvas, run_loop, pHist);
///    std::this_thread::sleep_for(std::chrono::seconds(100));
///    run_loop = false;
///    thrd.join();
///    return 0;
/// }
/// ~~~

void ROOT::Experimental::RCanvas::Run(double tm)
{
   if (fPainter) {
      fPainter->Run(tm);
   } else if (tm>0) {
      std::this_thread::sleep_for(std::chrono::milliseconds(int(tm*1000)));
   }
}

//////////////////////////////////////////////////////////////////////////
/// To resolve problem with storing of shared pointers
/// Call this method when reading canvas from the file
/// Can be called many times - after reinitialization of shared pointers no changes will be performed

void ROOT::Experimental::RCanvas::ResolveSharedPtrs()
{
   Internal::RIOSharedVector_t vect;

   CollectShared(vect);

   for (unsigned n = 0; n < vect.size(); ++n) {
      if (vect[n]->HasShared() || !vect[n]->GetIOPtr()) continue;

      auto shrd_ptr = vect[n]->MakeShared();

      for (auto n2 = n+1; n2 < vect.size(); ++n2) {
         if (vect[n2]->GetIOPtr() == vect[n]->GetIOPtr()) {
            if (vect[n2]->HasShared())
               R__LOG_ERROR(GPadLog()) << "FATAL Shared pointer for same IO ptr already exists";
            else
               vect[n2]->SetShared(shrd_ptr);
         }
      }

   }
}


/////////////////////////////////////////////////////////////////////////////////////////////////
/// Apply attributes changes to the drawable
/// Return mask with actions which were really applied

std::unique_ptr<ROOT::Experimental::RDrawableReply> ROOT::Experimental::RChangeAttrRequest::Process()
{
   // suppress all changes coming from non-main connection
   if (!GetContext().IsMainConn())
      return nullptr;

   auto canv = const_cast<ROOT::Experimental::RCanvas *>(GetContext().GetCanvas());
   if (!canv) return nullptr;

   if ((ids.size() != names.size()) || (ids.size() != values.size())) {
      R__LOG_ERROR(GPadLog()) << "Mismatch of arrays size in RChangeAttrRequest";
      return nullptr;
   }

   Version_t vers = 0;

   for(int indx = 0; indx < (int) ids.size(); indx++) {
      if (ids[indx] == "canvas") {
         if (canv->GetAttrMap().Change(names[indx], values[indx].get())) {
            if (!vers) vers = canv->IncModified();
            canv->SetDrawableVersion(vers);
         }
      } else {
         auto drawable = canv->FindPrimitiveByDisplayId(ids[indx]);
         if (drawable && drawable->GetAttrMap().Change(names[indx], values[indx].get())) {
            if (!vers) vers = canv->IncModified();
            drawable->SetDrawableVersion(vers);
         }
      }
   }

   fNeedUpdate = (vers > 0) && update;

   return nullptr; // no need for any reply
}
