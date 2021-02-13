// @(#)root/graf2d:$Id$
// Author: Timur Pocheptsov   29/02/2012

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//#define NDEBUG

#include <stdexcept>
#include <cstring>
#include <cassert>
#include <memory>

#include <Cocoa/Cocoa.h>

#include "ROOTOpenGLView.h"
#include "CocoaPrivate.h"
#include "QuartzWindow.h"
#include "QuartzPixmap.h"
#include "QuartzUtils.h"
#include "X11Drawable.h"
#include "X11Buffer.h"
#include "TGWindow.h"
#include "TGClient.h"
#include "TGCocoa.h"
#include "TError.h"

namespace ROOT {
namespace MacOSX {
namespace X11 {

//______________________________________________________________________________
Command::Command(Drawable_t wid, const GCValues_t &gc)
            : fID(wid),
              fGC(gc)
{
}

//______________________________________________________________________________
Command::Command(Drawable_t wid)
            : fID(wid),
              fGC()
{
}

//______________________________________________________________________________
Command::~Command()
{
}

//______________________________________________________________________________
void Command::Execute(CGContextRef /*ctx*/)const
{
}

//______________________________________________________________________________
bool Command::HasOperand(Drawable_t wid)const
{
   return wid == fID;
}

//______________________________________________________________________________
bool Command::IsGraphicsCommand()const
{
   return false;
}

//______________________________________________________________________________
DrawLine::DrawLine(Drawable_t wid, const GCValues_t &gc, const Point &p1, const Point &p2)
            : Command(wid, gc),
              fP1(p1),
              fP2(p2)
{
}

//______________________________________________________________________________
void DrawLine::Execute()const
{
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "Execute, gVirtualX is either null or not of TGCocoa type");
   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->DrawLineAux(fID, fGC, fP1.fX, fP1.fY, fP2.fX, fP2.fY);
}

//______________________________________________________________________________
DrawSegments::DrawSegments(Drawable_t wid, const GCValues_t &gc,
                           const Segment_t *segments, Int_t nSegments)
                 : Command(wid, gc)
{
   assert(segments != 0 && "DrawSegments, segments parameter is null");
   assert(nSegments > 0 && "DrawSegments, nSegments <= 0");

   fSegments.assign(segments, segments + nSegments);
}

//______________________________________________________________________________
void DrawSegments::Execute()const
{
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "Execute, gVirtualX is either null or not of TGCocoa type");
   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->DrawSegmentsAux(fID, fGC, &fSegments[0], (Int_t)fSegments.size());
}

//______________________________________________________________________________
ClearArea::ClearArea(Window_t wid, const Rectangle_t &area)
             : Command(wid),
               fArea(area)
{
}

//______________________________________________________________________________
void ClearArea::Execute()const
{
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "Execute, gVirtualX is either null or not of TGCocoa type");
   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->ClearAreaAux(fID, fArea.fX, fArea.fY, fArea.fWidth, fArea.fHeight);
}

//______________________________________________________________________________
CopyArea::CopyArea(Drawable_t src, Drawable_t dst, const GCValues_t &gc,
                   const Rectangle_t &area, const Point &dstPoint)
               : Command(dst, gc),
                 fSrc(src),
                 fArea(area),
                 fDstPoint(dstPoint)
{
}

//______________________________________________________________________________
bool CopyArea::HasOperand(Drawable_t drawable)const
{
   return fID == drawable || fSrc == drawable || fGC.fClipMask == drawable;
}

//______________________________________________________________________________
void CopyArea::Execute()const
{
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "Execute, gVirtualX is either null or not of TGCocoa type");

   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->CopyAreaAux(fSrc, fID, fGC, fArea.fX, fArea.fY, fArea.fWidth,
                   fArea.fHeight, fDstPoint.fX, fDstPoint.fY);
}

//______________________________________________________________________________
DrawString::DrawString(Drawable_t wid, const GCValues_t &gc, const Point &point,
                       const std::string &text)
               : Command(wid, gc),
                 fPoint(point),
                 fText(text)
{
}

//______________________________________________________________________________
void DrawString::Execute()const
{
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "Execute, gVirtualX is either null or not of TGCocoa type");

   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->DrawStringAux(fID, fGC, fPoint.fX, fPoint.fY, fText.c_str(), fText.length());
}

//______________________________________________________________________________
FillRectangle::FillRectangle(Drawable_t wid, const GCValues_t &gc,
                             const Rectangle_t &rectangle)
                  : Command(wid, gc),
                    fRectangle(rectangle)
{
}

//______________________________________________________________________________
void FillRectangle::Execute()const
{
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "Execute, gVirtualX is either null or not of TGCocoa type");

   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->FillRectangleAux(fID, fGC, fRectangle.fX, fRectangle.fY,
                        fRectangle.fWidth, fRectangle.fHeight);
}

//______________________________________________________________________________
FillPolygon::FillPolygon(Drawable_t wid, const GCValues_t &gc,
                         const Point_t *points, Int_t nPoints)
                : Command(wid, gc)
{
   assert(points != 0 && "FillPolygon, points parameter is null");
   assert(nPoints > 0 && "FillPolygon, nPoints <= 0");

   fPolygon.assign(points, points + nPoints);
}

//______________________________________________________________________________
void FillPolygon::Execute()const
{
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "Execute, gVirtualX is either null or not of TGCocoa type");

   ((TGCocoa *)gVirtualX)->FillPolygonAux(fID, fGC, &fPolygon[0], (Int_t)fPolygon.size());
}

//______________________________________________________________________________
DrawRectangle::DrawRectangle(Drawable_t wid, const GCValues_t &gc,
                             const Rectangle_t &rectangle)
                 : Command(wid, gc),
                   fRectangle(rectangle)
{
}

//______________________________________________________________________________
void DrawRectangle::Execute()const
{
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "Execute, gVirtualX is either null or not of TGCocoa type");

   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->DrawRectangleAux(fID, fGC, fRectangle.fX, fRectangle.fY,
                        fRectangle.fWidth, fRectangle.fHeight);
}

//______________________________________________________________________________
UpdateWindow::UpdateWindow(QuartzView *view)
                : Command(view.fID),
                  fView(view)
{
   assert(view != nil && "UpdateWindow, view parameter is nil");//view.fID will be also 0.
}

//______________________________________________________________________________
void UpdateWindow::Execute()const
{
   assert(fView.fContext != 0 && "Execute, view.fContext is null");

   if (QuartzPixmap *pixmap = fView.fBackBuffer)
      [fView copy : pixmap area : Rectangle(0, 0, pixmap.fWidth, pixmap.fHeight)
         withMask : nil clipOrigin : Point() toPoint : Point()];
}

//______________________________________________________________________________
DeletePixmap::DeletePixmap(Pixmap_t pixmap)
                : Command(pixmap, GCValues_t())
{
}

//______________________________________________________________________________
void DeletePixmap::Execute()const
{
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "Execute, gVirtualX is either null or not of TGCocoa type");

   ((TGCocoa *)gVirtualX)->DeletePixmapAux(fID);
}

//______________________________________________________________________________
DrawBoxXor::DrawBoxXor(Window_t windowID, const Point &p1, const Point &p2)
               : Command(windowID, GCValues_t()),
                 fP1(p1),
                 fP2(p2)
{
   if (fP1.fX > fP2.fX)
      std::swap(fP1.fX, fP2.fX);
   if (fP1.fY > fP2.fY)
      std::swap(fP1.fY, fP2.fY);
}

//______________________________________________________________________________
void DrawBoxXor::Execute()const
{
   //Noop.
}

//______________________________________________________________________________
void DrawBoxXor::Execute(CGContextRef ctx)const
{
   //
   assert(ctx != 0 && "Execute, ctx parameter is null");

   CGContextSetRGBStrokeColor(ctx, 0., 0., 0., 1.);
   CGContextSetLineWidth(ctx, 1.);

   CGContextStrokeRect(ctx, CGRectMake(fP1.fX, fP1.fY, fP2.fX - fP1.fX, fP2.fY - fP1.fY));
}

//______________________________________________________________________________
DrawLineXor::DrawLineXor(Window_t windowID, const Point &p1, const Point &p2)
               : Command(windowID, GCValues_t()),
                 fP1(p1),
                 fP2(p2)
{
}

//______________________________________________________________________________
void DrawLineXor::Execute()const
{
   //Noop.
}

//______________________________________________________________________________
void DrawLineXor::Execute(CGContextRef)const
{
   //Noop.
}

//______________________________________________________________________________
CommandBuffer::CommandBuffer()
{
}

//______________________________________________________________________________
CommandBuffer::~CommandBuffer()
{
   ClearCommands();
   ClearXOROperations();
}

//______________________________________________________________________________
void CommandBuffer::AddDrawLine(Drawable_t wid, const GCValues_t &gc, Int_t x1,
                                Int_t y1, Int_t x2, Int_t y2)
{
   try {
      //if this throws, I do not care.
      std::unique_ptr<DrawLine> cmd(new DrawLine(wid, gc, Point(x1, y1), Point(x2, y2)));
      fCommands.push_back(cmd.get());//this can throw.
      cmd.release();
   } catch (const std::exception &) {
      throw;
   }
}

//______________________________________________________________________________
void CommandBuffer::AddDrawSegments(Drawable_t wid, const GCValues_t &gc,
                                    const Segment_t *segments, Int_t nSegments)
{
   assert(segments != 0 && "AddDrawSegments, segments parameter is null");
   assert(nSegments > 0 && "AddDrawSegments, nSegments <= 0");

   try {
      std::unique_ptr<DrawSegments> cmd(new DrawSegments(wid, gc, segments, nSegments));
      fCommands.push_back(cmd.get());
      cmd.release();
   } catch (const std::exception &) {
      throw;
   }
}

//______________________________________________________________________________
void CommandBuffer::AddClearArea(Window_t wid, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   try {
      Rectangle_t r = {};//To be replaced with X11::Rectangle.
      r.fX = x;
      r.fY = y;
      r.fWidth = (UShort_t)w;
      r.fHeight = (UShort_t)h;
      std::unique_ptr<ClearArea> cmd(new ClearArea(wid, r));//Can throw, nothing leaks.
      fCommands.push_back(cmd.get());//this can throw.
      cmd.release();
   } catch (const std::exception &) {
      throw;
   }
}

//______________________________________________________________________________
void CommandBuffer::AddCopyArea(Drawable_t src, Drawable_t dst, const GCValues_t &gc,
                                Int_t srcX, Int_t srcY, UInt_t width, UInt_t height,
                                Int_t dstX, Int_t dstY)
{
   try {
      Rectangle_t area = {};
      area.fX = srcX;
      area.fY = srcY;
      area.fWidth = (UShort_t)width;
      area.fHeight = (UShort_t)height;
      //Can throw, nothing leaks.
      std::unique_ptr<CopyArea> cmd(new CopyArea(src, dst, gc, area, Point(dstX, dstY)));
      fCommands.push_back(cmd.get());//this can throw.
      cmd.release();
   } catch (const std::exception &) {
      throw;
   }
}

//______________________________________________________________________________
void CommandBuffer::AddDrawString(Drawable_t wid, const GCValues_t &gc, Int_t x, Int_t y,
                                  const char *text, Int_t len)
{
   try {
      if (len < 0)//Negative length can come from caller.
         len = std::strlen(text);
      const std::string substr(text, len);//Can throw.
      std::unique_ptr<DrawString> cmd(new DrawString(wid, gc, Point(x, y), substr));//Can throw.
      fCommands.push_back(cmd.get());//can throw.
      cmd.release();
   } catch (const std::exception &) {
      throw;
   }
}

//______________________________________________________________________________
void CommandBuffer::AddFillRectangle(Drawable_t wid, const GCValues_t &gc,
                                     Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   try {
      Rectangle_t r = {};
      r.fX = x;
      r.fY = y;
      r.fWidth = (UShort_t)w;
      r.fHeight = (UShort_t)h;
      std::unique_ptr<FillRectangle> cmd(new FillRectangle(wid, gc, r));
      fCommands.push_back(cmd.get());
      cmd.release();
   } catch (const std::exception &) {
      throw;
   }
}

//______________________________________________________________________________
void CommandBuffer::AddDrawRectangle(Drawable_t wid, const GCValues_t &gc,
                                     Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   try {
      Rectangle_t r = {};
      r.fX = x;
      r.fY = y;
      r.fWidth = (UShort_t)w;
      r.fHeight = (UShort_t)h;
      std::unique_ptr<DrawRectangle> cmd(new DrawRectangle(wid, gc, r));
      fCommands.push_back(cmd.get());
      cmd.release();
   } catch (const std::exception &) {
      throw;
   }
}

//______________________________________________________________________________
void CommandBuffer::AddFillPolygon(Drawable_t wid, const GCValues_t &gc,
                                   const Point_t *polygon, Int_t nPoints)
{
   assert(polygon != 0 && "AddFillPolygon, polygon parameter is null");
   assert(nPoints > 0 && "AddFillPolygon, nPoints <= 0");

   try {
      std::unique_ptr<FillPolygon> cmd(new FillPolygon(wid, gc, polygon, nPoints));
      fCommands.push_back(cmd.get());
      cmd.release();
   } catch (const std::exception &) {
      throw;
   }
}

//______________________________________________________________________________
void CommandBuffer::AddUpdateWindow(QuartzView *view)
{
   assert(view != nil && "AddUpdateWindow, view parameter is nil");

   try {
      std::unique_ptr<UpdateWindow> cmd(new UpdateWindow(view));
      fCommands.push_back(cmd.get());
      cmd.release();
   } catch (const std::exception &) {
      throw;
   }
}

//______________________________________________________________________________
void CommandBuffer::AddDeletePixmap(Pixmap_t pixmapID)
{
   try {
      std::unique_ptr<DeletePixmap> cmd(new DeletePixmap(pixmapID));
      fCommands.push_back(cmd.get());
      cmd.release();
   } catch (const std::exception &) {
      throw;
   }
}

//______________________________________________________________________________
void CommandBuffer::AddDrawBoxXor(Window_t windowID, Int_t x1, Int_t y1, Int_t x2, Int_t y2)
{
   try {
      std::unique_ptr<DrawBoxXor> cmd(new DrawBoxXor(windowID, Point(x1, y1), Point(x2, y2)));
      fXorOps.push_back(cmd.get());
      cmd.release();
   } catch (const std::exception &) {
      throw;
   }
}

//______________________________________________________________________________
void CommandBuffer::AddDrawLineXor(Window_t windowID, Int_t x1, Int_t y1, Int_t x2, Int_t y2)
{
   try {
      std::unique_ptr<DrawLineXor> cmd(new DrawLineXor(windowID, Point(x1, y1), Point(x2, y2)));
      fXorOps.push_back(cmd.get());
      cmd.release();
   } catch (const std::exception &) {
      throw;
   }
}

//______________________________________________________________________________
void CommandBuffer::Flush(Details::CocoaPrivate *impl)
{
   assert(impl != 0 && "Flush, impl parameter is null");

   //Basic es-guarantee: state is unknown, but valid, no
   //resource leaks, no locked focus.

   //All magic is here.
   CGContextRef prevContext = 0;
   CGContextRef currContext = 0;
   QuartzView *prevView = nil;

   for (size_type i = 0, e = fCommands.size(); i < e; ++i) {
      const Command *cmd = fCommands[i];
      if (!cmd)//Command was deleted by RemoveOperation/RemoveGraphicsOperation.
         continue;

      NSObject<X11Drawable> *drawable = impl->GetDrawable(cmd->fID);
      if (drawable.fIsPixmap) {
         cmd->Execute();//Can throw, ok.
         continue;
      }

      QuartzView *view = (QuartzView *)impl->GetWindow(cmd->fID).fContentView;

      if (prevView != view)
         ClipOverlaps(view);//Can throw, ok.

      prevView = view;

      try {
         if ([view lockFocusIfCanDraw]) {
            NSGraphicsContext *nsContext = [NSGraphicsContext currentContext];
            assert(nsContext != nil && "Flush, currentContext is nil");
            currContext = (CGContextRef)[nsContext graphicsPort];
            assert(currContext != 0 && "Flush, graphicsPort is null");//remove this assert?

            view.fContext = currContext;
            if (prevContext && prevContext != currContext)
               CGContextFlush(prevContext);
            prevContext = currContext;

            const Quartz::CGStateGuard ctxGuard(currContext);

            //Clip regions first.
            if (fClippedRegion.size())
               CGContextClipToRects(currContext, &fClippedRegion[0], fClippedRegion.size());

            //Now add also shape combine mask.
            if (view.fQuartzWindow.fShapeCombineMask)
               ClipToShapeMask(view, currContext);

            cmd->Execute();//This can throw, we should restore as much as we can here.

            if (view.fBackBuffer) {
               //Very "special" window.
               const Rectangle copyArea(0, 0, view.fBackBuffer.fWidth, view.fBackBuffer.fHeight);
               [view copy : view.fBackBuffer area : copyArea
                 withMask : nil clipOrigin : Point() toPoint : Point()];
            }

            [view unlockFocus];

            view.fContext = 0;
         }
      } catch (const std::exception &) {
         //Focus was locked, roll-back:
         [view unlockFocus];
         //View's context was modified, roll-back:
         view.fContext = 0;
         //Re-throw, something really bad happened (std::bad_alloc).
         throw;
      }
   }

   if (currContext)
      CGContextFlush(currContext);

   ClearCommands();
}

//______________________________________________________________________________
void CommandBuffer::FlushXOROps(Details::CocoaPrivate *impl)
{
   // The only XOR operations we ever had to support was drawing
   // lines of a crosshair in a TCanvas. We were using a deprecated
   // (since 10.14) trick with locking a focus on a view, drawing,
   // flushing CGContext and then unlocking. This is not working
   // starting from 10.15. So now the only thing we do - we draw
   // a crosshair into the special transparent window which is
   // attached on top of the canvas.
   assert(impl != 0 && "FlushXOROps, impl parameter is null");

   if (fXorOps.size() < 2) {
       ClearXOROperations();
       return;
   }

   NSObject<X11Drawable> * const drawable = impl->GetDrawable(fXorOps.back()->fID);
   assert([drawable isKindOfClass : [QuartzView class]] &&
          "FlushXOROps, drawable must be of type QuartzView");
   QuartzView * const view = (QuartzView *)drawable;
   QuartzWindow * const window = view.fQuartzWindow;
   auto crosshairWindow = [window findCrosshairWindow];
   if (!crosshairWindow) {
      ::Warning("FlushXOROps", "No CrosshairWindow found to draw into");
       ClearXOROperations();
       return;
   }

   for (auto *candidateOp : fXorOps) {
       if (!dynamic_cast<DrawLineXor *>(candidateOp)) {
           NSLog(@"XOR mode supports DrawLine only, unexpected operation was found");
           ClearXOROperations();
           return;
       }
   }

   // If size is > 2, we take the last two operations.
   DrawLineXor *line1 = static_cast<DrawLineXor *>(*(fXorOps.end() - 2));
   DrawLineXor *line2 = static_cast<DrawLineXor *>(fXorOps.back());

   auto rootToNs = [](Point rp) {
      return NSPoint{CGFloat(rp.fX), CGFloat(rp.fY)};
   };

   NSPoint cross[] = {rootToNs(line1->start()), rootToNs(line1->end()),
                      rootToNs(line2->start()), rootToNs(line2->end())};

   for (auto &point : cross) {
      // From the view's coordinates to window's:
      point = [view convertPoint : point toView : nil];
   }

   CrosshairView *cv = (CrosshairView *)crosshairWindow.contentView;
   cv.start1 = cross[0];
   cv.end1 = cross[1];
   cv.start2 = cross[2];
   cv.end2 = cross[3];
   [cv setNeedsDisplay : YES];

   ClearXOROperations();
}

//______________________________________________________________________________
void CommandBuffer::RemoveOperationsForDrawable(Drawable_t drawable)
{
   for (size_type i = 0; i < fCommands.size(); ++i) {
      if (fCommands[i] && fCommands[i]->HasOperand(drawable)) {
         delete fCommands[i];
         fCommands[i] = 0;
      }
   }

   for (size_type i = 0; i < fXorOps.size(); ++i) {
      if (fXorOps[i] && fXorOps[i]->HasOperand(drawable)) {
         delete fXorOps[i];
         fXorOps[i] = 0;
      }
   }
}

//______________________________________________________________________________
void CommandBuffer::RemoveGraphicsOperationsForWindow(Window_t wid)
{
   for (size_type i = 0; i < fCommands.size(); ++i) {
      if (fCommands[i] && fCommands[i]->HasOperand(wid) &&
          fCommands[i]->IsGraphicsCommand())
      {
         delete fCommands[i];
         fCommands[i] = 0;
      }
   }
}

//______________________________________________________________________________
void CommandBuffer::RemoveXORGraphicsOperationsForWindow(Window_t wid)
{
   for (size_type i = 0; i < fCommands.size(); ++i) {
      if (fXorOps[i] && fXorOps[i]->HasOperand(wid)) {
         delete fXorOps[i];
         fXorOps[i] = 0;
      }
   }
}

//______________________________________________________________________________
void CommandBuffer::ClearCommands()
{
   for (size_type i = 0, e = fCommands.size(); i < e; ++i)
      delete fCommands[i];

   fCommands.clear();
}

//______________________________________________________________________________
void CommandBuffer::ClearXOROperations()
{
   for (size_type i = 0, e = fXorOps.size(); i < e; ++i)
      delete fXorOps[i];

   fXorOps.clear();
}

//Clipping machinery.

namespace {

//________________________________________________________________________________________
bool RectsOverlap(const NSRect &r1, const NSRect &r2)
{
   if (r2.origin.x >= r1.origin.x + r1.size.width)
      return false;
   if (r2.origin.x + r2.size.width <= r1.origin.x)
      return false;
   if (r2.origin.y >= r1.origin.y + r1.size.height)
      return false;
   if (r2.origin.y + r2.size.height <= r1.origin.y)
      return false;

   return true;
}

}

//______________________________________________________________________________
void CommandBuffer::ClipOverlaps(QuartzView *view)
{
   //QuartzViews do not have backing store.
   //But ROOT calls gClient->NeedRedraw ignoring
   //children or overlapping siblings. This leads
   //to obvious problems, for example, parent
   //erasing every child inside while repainting itself.
   //To fix this and emulate window with backing store
   //without real backing store, I'm calculating the
   //area of a view this is visible and not overlapped.

   //Who can overlap our view?
   //1. Its own siblings and, probably, siblings of its ancestors.
   //2. Children views.

   assert(view != nil && "ClipOverlaps, view parameter is nil");

   typedef std::vector<QuartzView *>::reverse_iterator reverse_iterator;
   typedef std::vector<CGRect>::iterator rect_iterator;

   fRectsToClip.clear();
   fClippedRegion.clear();

   //Check siblings and ancestors' siblings:

   //1. Remember the whole branch starting from our view
   //up to a top-level window.
   fViewBranch.clear();
   for (QuartzView *v = view; v; v = v.fParentView)
      fViewBranch.push_back(v);

   //We do not need content view, since it does not have any siblings.
   if (fViewBranch.size())
      fViewBranch.pop_back();

   //For every fViewBranch[i] in our branch, we're looking for overlapping siblings.
   //Calculations are in view.fParentView's coordinate system.

   WidgetRect clipRect;
   NSRect frame1 = {};

   const NSRect frame2 = view.frame;

   for (reverse_iterator it = fViewBranch.rbegin(), eIt = fViewBranch.rend(); it != eIt; ++it) {
      QuartzView *ancestorView = *it;//This is either one of ancestors, or a view itself.
      bool doCheck = false;
      for (QuartzView *sibling in [ancestorView.fParentView subviews]) {
         if (ancestorView == sibling) {
            //View has its children in an array, and for every subviews[i] in this array,
            //only views with index > i can overlap subviews[i].
            doCheck = true;//all views after this must be checked.
            continue;
         } else if (!doCheck || sibling.fMapState != kIsViewable) {
            continue;
         }

         frame1 = sibling.frame;

         if (!frame1.size.width || !frame1.size.height)
            continue;

         frame1.origin = [sibling.fParentView convertPoint : frame1.origin
                          toView : view.fParentView];

         //Check if two rects intersect.
         if (RectsOverlap(frame2, frame1)) {
            //Substruct frame1 from our view's rect.
            clipRect.fX1 = frame1.origin.x;
            clipRect.fX2 = clipRect.fX1 + frame1.size.width;
            clipRect.fY1 = frame1.origin.y;
            clipRect.fY2 = clipRect.fY1 + frame1.size.height;
            fRectsToClip.push_back(clipRect);
         }
      }
   }

   //Substruct children.

   for (QuartzView *child in [view subviews]) {
      if (child.fMapState != kIsViewable)
         continue;

      frame1 = child.frame;

      if (!frame1.size.width || !frame1.size.height)
         continue;

      if (view.fParentView)//view can also be a content view.
         frame1.origin = [view convertPoint : frame1.origin toView : view.fParentView];

      if (RectsOverlap(frame2, frame1)) {
         clipRect.fX1 = frame1.origin.x;
         clipRect.fX2 = clipRect.fX1 + frame1.size.width;
         clipRect.fY1 = frame1.origin.y;
         clipRect.fY2 = clipRect.fY1 + frame1.size.height;
         fRectsToClip.push_back(clipRect);
      }
   }

   if (fRectsToClip.size()) {
      //Now, if we have any rectanges to substruct them from our view's frame,
      //we are building a set of rectangles, which represents visible part of view.

      WidgetRect rect(frame2.origin.x, frame2.origin.y, frame2.origin.x + frame2.size.width,
                      frame2.origin.y + frame2.size.height);

      BuildClipRegion(rect);

      if (view.fParentView) {
         //To able to use this set of rectangles with CGContextClipToRects,
         //convert them (if needed) into view's own coordinate system.
         rect_iterator recIt = fClippedRegion.begin(), eIt = fClippedRegion.end();
         for (; recIt != eIt; ++recIt) {
            if (!recIt->size.width && !recIt->size.height) {
               //This is a special 'empty' rectangle, which means our view is completely hidden.
               assert(fClippedRegion.size() == 1 && "ClipOverlaps, internal logic error");
               break;
            }
            recIt->origin = NSPointToCGPoint([view.fParentView convertPoint :
                                              NSPointFromCGPoint(recIt->origin) toView : view]);
         }
      }
   }
}

namespace {

typedef std::vector<int>::iterator int_iterator;

//_____________________________________________________________________________________________________
int_iterator BinarySearchLeft(int_iterator first, int_iterator last, int value)
{
   if (first == last)
      return last;

   const int_iterator it = std::lower_bound(first, last, value);
   assert(it != last && (it == first || *it == value) && "internal logic error");

   //If value < *first, return last (not found).
   return it == first && *it != value ? last : it;
}

//_____________________________________________________________________________________________________
int_iterator BinarySearchRight(int_iterator first, int_iterator last, int value)
{
   if (first == last)
      return last;

   const int_iterator it = std::lower_bound(first, last, value);
   assert((it == last || *it == value) && "internal logic error");

   return it;
}

}//unnamed namespace.

//_____________________________________________________________________________________________________
void CommandBuffer::BuildClipRegion(const WidgetRect &rect)
{
   //Input requirements:
   // 1) all rects are valid (non-empty and x1 < x2, y1 < y2);
   // 2) all rects intersect with widget's rect.
   //I do not check these conditions here, this is done when filling rectsToClip.

   //I did not find any reasonable algorithm (have to search better?),
   //code in gdk and pixman has to many dependencies and is lib-specific +
   //they require input to be quite special:
   // a) no overlaps (in my case I have overlaps)
   // b) sorted in a special way.
   //To convert my input into such a format
   //means to implement everything myself (for example, to work out overlaps).

   //Also, my case is more simple: gdk and pixman substract region (== set of rectangles)
   //from another region, I have to substract region from _one_ rectangle.

   //This is quite straightforward implementation - I'm calculation rectangles, which are part of
   //a widget's rect, not hidden by any of fRectsToClip.

   typedef std::vector<WidgetRect>::const_iterator rect_const_iterator;

   assert(fRectsToClip.size() != 0 && "BuildClipRegion, nothing to clip");

   fClippedRegion.clear();
   fXBounds.clear();
   fYBounds.clear();

   //[First, we "cut" the original rect into stripes.
   rect_const_iterator recIt = fRectsToClip.begin(), endIt = fRectsToClip.end();
   for (; recIt != endIt; ++recIt) {
      if (recIt->fX1 <= rect.fX1 && recIt->fX2 >= rect.fX2 &&
          recIt->fY1 <= rect.fY1 && recIt->fY2 >= rect.fY2) {
         //this rect completely overlaps our view, not need to calculate anything at all.
         fClippedRegion.push_back(CGRectMake(0., 0., 0., 0.));
         return;
      }

      if (recIt->fX1 > rect.fX1)//recIt->x1 is always < rect.x2 (input validation).
         fXBounds.push_back(recIt->fX1);

      if (recIt->fX2 < rect.fX2)//recIt->x2 is always > rect.x1 (input validation).
         fXBounds.push_back(recIt->fX2);

      if (recIt->fY1 > rect.fY1)
         fYBounds.push_back(recIt->fY1);

      if (recIt->fY2 < rect.fY2)
         fYBounds.push_back(recIt->fY2);
   }

   std::sort(fXBounds.begin(), fXBounds.end());
   std::sort(fYBounds.begin(), fYBounds.end());

   //We do not need duplicates.
   const int_iterator xBoundsEnd = std::unique(fXBounds.begin(), fXBounds.end());
   const int_iterator yBoundsEnd = std::unique(fYBounds.begin(), fYBounds.end());
   //Rectangle is now "cut into pieces"].

   const size_type nXBands = size_type(xBoundsEnd - fXBounds.begin()) + 1;
   const size_type nYBands = size_type(yBoundsEnd - fYBounds.begin()) + 1;

   fGrid.assign(nXBands * nYBands, false);

   //Mark the overlapped parts.
   recIt = fRectsToClip.begin(), endIt = fRectsToClip.end();
   for (; recIt != endIt; ++recIt) {
      const int_iterator left = BinarySearchLeft(fXBounds.begin(), xBoundsEnd, recIt->fX1);
      const size_type firstXBand = left == xBoundsEnd ? 0 : left - fXBounds.begin() + 1;

      const int_iterator right = BinarySearchRight(fXBounds.begin(), xBoundsEnd, recIt->fX2);
      const size_type lastXBand = right - fXBounds.begin() + 1;

      const int_iterator bottom = BinarySearchLeft(fYBounds.begin(), yBoundsEnd, recIt->fY1);
      const size_type firstYBand = bottom == yBoundsEnd ? 0 : bottom - fYBounds.begin() + 1;

      const int_iterator top = BinarySearchRight(fYBounds.begin(), yBoundsEnd, recIt->fY2);
      const size_type lastYBand = top - fYBounds.begin() + 1;

      for (size_type i = firstYBand; i < lastYBand; ++i) {
         const size_type baseIndex = i * nXBands;
         for (size_type j = firstXBand; j < lastXBand; ++j)
            fGrid[baseIndex + j] = true;
      }
   }

   //I do not merge rectangles.
   //Search for non-overlapped parts and create rectangles for them.
   CGRect newRect = {};

   for (size_type i = 0; i < nYBands; ++i) {
      const size_type baseIndex = i * nXBands;
      for (size_type j = 0; j < nXBands; ++j) {
         if (!fGrid[baseIndex + j]) {
            newRect.origin.x = j ? fXBounds[j - 1] : rect.fX1;
            newRect.origin.y = i ? fYBounds[i - 1] : rect.fY1;

            newRect.size.width = (j == nXBands - 1 ? rect.fX2 : fXBounds[j]) - newRect.origin.x;
            newRect.size.height = (i == nYBands - 1 ? rect.fY2 : fYBounds[i]) - newRect.origin.y;

            fClippedRegion.push_back(newRect);
         }
      }
   }

   if (!fClippedRegion.size())//Completely hidden
      fClippedRegion.push_back(CGRectMake(0., 0., 0., 0.));
}

}//X11
}//MacOSX
}//ROOT
