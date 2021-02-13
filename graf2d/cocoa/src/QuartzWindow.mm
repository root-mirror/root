// @(#)root/graf2d:$Id$
// Author: Timur Pocheptsov   16/02/2012

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//#define DEBUG_ROOT_COCOA

//#define NDEBUG

#ifdef DEBUG_ROOT_COCOA
#include <iostream>
#include <fstream>

#include "TClass.h"
#endif

#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <vector>

#include <Availability.h>

#include "ROOTOpenGLView.h"
#include "CocoaConstants.h"
#include "QuartzWindow.h"
#include "QuartzPixmap.h"
#include "QuartzUtils.h"
#include "CocoaUtils.h"
#include "RConfigure.h"
#include "X11Colors.h"
#include "X11Buffer.h"
#include "TGWindow.h"
#include "TGClient.h"
#include "TSystem.h"
#include "TGCocoa.h"
#include "TROOT.h"
#include "TGTextView.h"
#include "TGView.h"
#include "TGCanvas.h"


namespace ROOT {
namespace MacOSX {
namespace X11 {

#pragma mark - Create a window or a view.

//______________________________________________________________________________
QuartzWindow *CreateTopLevelWindow(Int_t x, Int_t y, UInt_t w, UInt_t h, UInt_t /*border*/, Int_t depth,
                                   UInt_t clss, void */*visual*/, SetWindowAttributes_t *attr, UInt_t)
{
   using namespace Details;

   NSRect winRect = {};
   winRect.origin.x = GlobalXROOTToCocoa(x);
   winRect.origin.y = GlobalYROOTToCocoa(y + h);
   winRect.size.width = w;
   winRect.size.height = h;

   const NSUInteger styleMask = kTitledWindowMask | kClosableWindowMask |
                                kMiniaturizableWindowMask | kResizableWindowMask;

   QuartzWindow * const newWindow = [[QuartzWindow alloc] initWithContentRect : winRect
                                                                    styleMask : styleMask
                                                                      backing : NSBackingStoreBuffered
                                                                        defer : YES
                                                             windowAttributes : attr];
   if (!newWindow)
      throw std::runtime_error("CreateTopLevelWindow failed");

   newWindow.fDepth = depth;
   newWindow.fClass = clss;

   return newWindow;
}

//______________________________________________________________________________
QuartzView *CreateChildView(QuartzView * /*parent*/, Int_t x, Int_t y, UInt_t w, UInt_t h, UInt_t /*border*/, Int_t /*depth*/,
                            UInt_t /*clss*/, void * /*visual*/, SetWindowAttributes_t *attr, UInt_t /*wtype*/)
{
   NSRect viewRect = {};
   viewRect.origin.x = x;
   viewRect.origin.y = y;
   viewRect.size.width = w;
   viewRect.size.height = h;

   QuartzView * const view = [[QuartzView alloc] initWithFrame : viewRect windowAttributes : attr];
   if (!view)
      throw std::runtime_error("CreateChildView failed");

   return view;
}

#pragma mark - root window (does not really exist, it's our desktop built of all screens).

//______________________________________________________________________________
void GetRootWindowAttributes(WindowAttributes_t *attr)
{
   //'root' window does not exist, but we can request its attributes.
   assert(attr != 0 && "GetRootWindowAttributes, parameter 'attr' is null");


   NSArray * const screens = [NSScreen screens];
   assert(screens != nil && "screens array is nil");
   NSScreen * const mainScreen = [screens objectAtIndex : 0];
   assert(mainScreen != nil && "screen with index 0 is nil");

   *attr = WindowAttributes_t();

   assert(dynamic_cast<TGCocoa *>(gVirtualX) &&
          "GetRootWindowAttributes, gVirtualX is either null or has a wrong type");

   TGCocoa * const gCocoa = static_cast<TGCocoa *>(gVirtualX);

   const Rectangle &frame = gCocoa->GetDisplayGeometry();

   attr->fX = 0;
   attr->fY = 0;
   attr->fWidth = frame.fWidth;
   attr->fHeight = frame.fHeight;
   attr->fBorderWidth = 0;
   attr->fYourEventMask = 0;
   attr->fAllEventMasks = 0;//???

   attr->fDepth = NSBitsPerPixelFromDepth([mainScreen depth]);
   attr->fVisual = 0;
   attr->fRoot = 0;
}


#pragma mark - Coordinate conversions.

//______________________________________________________________________________
NSPoint ConvertPointFromBaseToScreen(NSWindow *window, NSPoint windowPoint)
{
   assert(window != nil && "ConvertPointFromBaseToScreen, parameter 'window' is nil");

   //I have no idea why apple deprecated function for a point conversion and requires rect conversion,
   //point conversion seems to produce wrong results with HiDPI.

   NSRect tmpRect = {};
   tmpRect.origin = windowPoint;
   tmpRect.size = NSMakeSize(1., 1.);//This is strange size :) But if they require rect, 0,0 - will not work?
   tmpRect = [window convertRectToScreen : tmpRect];

   return tmpRect.origin;
}

//______________________________________________________________________________
NSPoint ConvertPointFromScreenToBase(NSPoint screenPoint, NSWindow *window)
{
   assert(window != nil && "ConvertPointFromScreenToBase, parameter 'window' is nil");

   //I have no idea why apple deprecated function for a point conversion and requires rect conversion,
   //point conversion seems to produce wrong results with HiDPI.

   NSRect tmpRect = {};
   tmpRect.origin = screenPoint;
   tmpRect.size = NSMakeSize(1., 1.);
   tmpRect = [window convertRectFromScreen : tmpRect];

   return tmpRect.origin;
}

//______________________________________________________________________________
int GlobalYCocoaToROOT(CGFloat yCocoa)
{
   //We can have several physical displays and thus several NSScreens in some arbitrary order.
   //With Cocoa, some screens can have negative coordinates - to the left ro down to the primary
   //screen (whatever it means). With X11 (XQuartz) though it's always 0,0.

   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "GlobalYCocoaToROOT, gVirtualX is either nul or has a wrong type");

   const Rectangle frame = ((TGCocoa *)gVirtualX)->GetDisplayGeometry();

   return int(frame.fHeight - (yCocoa - frame.fY));
}

//______________________________________________________________________________
int GlobalXCocoaToROOT(CGFloat xCocoa)
{
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "GlobalXCocoaToROOT, gVirtualX is either nul or has a wrong type");
   const Rectangle frame = ((TGCocoa *)gVirtualX)->GetDisplayGeometry();
   //With X11 coordinate space always starts from 0, 0
   return int(xCocoa - frame.fX);
}

//______________________________________________________________________________
int GlobalYROOTToCocoa(CGFloat yROOT)
{
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "GlobalYROOTToCocoa, gVirtualX is either nul or has a wrong type");
   const Rectangle frame = ((TGCocoa *)gVirtualX)->GetDisplayGeometry();

   return int(frame.fY + (frame.fHeight - yROOT));
}

//______________________________________________________________________________
int GlobalXROOTToCocoa(CGFloat xROOT)
{
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "GlobalXROOTToCocoa, gVirtualX is either nul or has a wrong type");
   const Rectangle frame = ((TGCocoa *)gVirtualX)->GetDisplayGeometry();
   //With X11 coordinate space always starts from 0, 0
   return int(frame.fX + xROOT);
}

//______________________________________________________________________________
int LocalYCocoaToROOT(NSView<X11Window> *parentView, CGFloat yCocoa)
{
   assert(parentView != nil && "LocalYCocoaToROOT, parent view is nil");

   return int(parentView.frame.size.height - yCocoa);
}

//______________________________________________________________________________
int LocalYROOTToCocoa(NSView<X11Window> *parentView, CGFloat yROOT)
{
   //:)
   assert(parentView != nil && "LocalYROOTToCocoa, parent view is nil");

   return int(parentView.frame.size.height - yROOT);
}


//______________________________________________________________________________
int LocalYROOTToCocoa(NSObject<X11Drawable> *drawable, CGFloat yROOT)
{
   //:)
   assert(drawable != nil && "LocalYROOTToCocoa, drawable is nil");

   return int(drawable.fHeight - yROOT);
}

//______________________________________________________________________________
NSPoint TranslateToScreen(NSView<X11Window> *from, NSPoint point)
{
   assert(from != nil && "TranslateToScreen, parameter 'from' is nil");

   const NSPoint winPoint = [from convertPoint : point toView : nil];

   NSPoint screenPoint = ConvertPointFromBaseToScreen([from window], winPoint);
   screenPoint.x = GlobalXCocoaToROOT(screenPoint.x);
   screenPoint.y = GlobalYCocoaToROOT(screenPoint.y);

   return screenPoint;
}

//______________________________________________________________________________
NSPoint TranslateFromScreen(NSPoint point, NSView<X11Window> *to)
{
   assert(to != nil && "TranslateFromScreen, parameter 'to' is nil");

   point.x = GlobalXROOTToCocoa(point.x);
   point.y = GlobalYROOTToCocoa(point.y);
   point = ConvertPointFromScreenToBase(point, [to window]);

   return [to convertPoint : point fromView : nil];
}

//______________________________________________________________________________
NSPoint TranslateCoordinates(NSView<X11Window> *from, NSView<X11Window> *to, NSPoint sourcePoint)
{
   //Both views are valid.
   assert(from != nil && "TranslateCoordinates, parameter 'from' is nil");
   assert(to != nil && "TranslateCoordinates, parameter 'to' is nil");

   if ([from window] == [to window]) {
      //Both views are in the same window.
      return [to convertPoint : sourcePoint fromView : from];
   } else {
      //May be, I can do it in one call, but it's not obvious for me
      //what is 'pixel aligned backing store coordinates' and
      //if they are the same as screen coordinates.

      //Many thanks to Apple for deprecated functions!!!

      const NSPoint win1Point = [from convertPoint : sourcePoint toView : nil];
      const NSPoint screenPoint = ConvertPointFromBaseToScreen([from window], win1Point);
      const NSPoint win2Point = ConvertPointFromScreenToBase(screenPoint, [to window]);

      return [to convertPoint : win2Point fromView : nil];
   }
}

//______________________________________________________________________________
bool ScreenPointIsInView(NSView<X11Window> *view, Int_t x, Int_t y)
{
   assert(view != nil && "ScreenPointIsInView, parameter 'view' is nil");

   NSPoint point = {};
   point.x = x, point.y = y;
   point = TranslateFromScreen(point, view);
   const NSRect viewFrame = view.frame;

   if (point.x < 0 || point.x > viewFrame.size.width)
      return false;
   if (point.y < 0 || point.y > viewFrame.size.height)
      return false;

   return true;
}

#pragma mark - Different FindView/Window functions iterating on the ROOT's windows/views.

//______________________________________________________________________________
QuartzWindow *FindWindowInPoint(Int_t x, Int_t y)
{
   //array's counter is increased, all object in array are also retained.
   const Util::AutoreleasePool pool;

   NSArray * const orderedWindows = [NSApp orderedWindows];
   for (NSWindow *window in orderedWindows) {
      if (![window isKindOfClass : [QuartzWindow class]])
         continue;
      QuartzWindow * const qw = (QuartzWindow *)window;
      if (qw.fIsDeleted)//Because of reference counting this can happen.
         continue;
      //Check if point is inside.
      if (ScreenPointIsInView(qw.fContentView, x, y))
         return qw;
   }

   return nil;
}

//______________________________________________________________________________
NSView<X11Window> *FindDNDAwareViewInPoint(NSArray *children, Window_t dragWinID,
                                           Window_t inputWinID, Int_t x, Int_t y, Int_t maxDepth)
{
   assert(children != nil && "FindDNDAwareViewInPoint, parameter 'children' is nil");

   if (maxDepth <= 0)
      return nil;

   NSEnumerator * const reverseEnumerator = [children reverseObjectEnumerator];
   for (NSView<X11Window> *child in reverseEnumerator) {
      if (!ScreenPointIsInView(child, x, y))
         continue;
      if (child.fIsDNDAware && child.fID != dragWinID && child.fID != inputWinID)
         return child;//got it!

      NSView<X11Window> * const testView = FindDNDAwareViewInPoint([child subviews], dragWinID,
                                                                   inputWinID, x, y, maxDepth - 1);
      if (testView)
         return testView;
   }

   return nil;
}

//______________________________________________________________________________
NSView<X11Window> *FindDNDAwareViewInPoint(NSView *parentView, Window_t dragWinID, Window_t inputWinID,
                                           Int_t x, Int_t y, Int_t maxDepth)
{
   //X and Y are ROOT's screen coordinates (Y is inverted).
   if (maxDepth <= 0)
      return nil;

   const Util::AutoreleasePool pool;

   if (!parentView) {//Start from the screen as a 'root' window.
      NSArray * const orderedWindows = [NSApp orderedWindows];
      for (NSWindow *window in orderedWindows) {
         if (![window isKindOfClass : [QuartzWindow class]])
            continue;
         QuartzWindow * const qw = (QuartzWindow *)window;

         if (qw.fIsDeleted)//Because of reference counting this can happen.
            continue;

         if (qw.fMapState != kIsViewable)
            continue;

         //First, check this view itself, my be we found what we need already.
         NSView<X11Window> *testView = qw.fContentView;
         if (!ScreenPointIsInView(testView, x, y))
            continue;

         if (testView.fIsDNDAware && testView.fID != dragWinID && testView.fID != inputWinID)
            return testView;

         //Recursive part, check children.
         NSArray * const children = [testView subviews];
         testView = FindDNDAwareViewInPoint(children, dragWinID, inputWinID, x, y, maxDepth - 1);
         if (testView)
            return testView;
      }

      //We did not find anything for 'root' window as parent.
      return nil;
   } else {
      //Parent view is tested already (or should not be tested at all, check children.
      return FindDNDAwareViewInPoint([parentView subviews], dragWinID, inputWinID, x, y, maxDepth);
   }
}

//______________________________________________________________________________
QuartzWindow *FindWindowUnderPointer()
{
   const Util::AutoreleasePool pool;

   NSArray * const orderedWindows = [NSApp orderedWindows];
   for (NSWindow *nsWindow in orderedWindows) {
      if (![nsWindow isKindOfClass : [QuartzWindow class]])
         continue;

      QuartzWindow * const qWindow = (QuartzWindow *)nsWindow;

      if (qWindow.fIsDeleted)//Because of reference counting this can happen.
         continue;

      if (qWindow.fMapState != kIsViewable)//Can it be false and still in this array???
         continue;

      const NSPoint mousePosition = [qWindow mouseLocationOutsideOfEventStream];
      const NSSize windowSize = qWindow.frame.size;
      if (mousePosition.x >= 0 && mousePosition.x <= windowSize.width &&
          mousePosition.y >= 0 && mousePosition.y <= windowSize.height)
         return qWindow;
   }

   return nil;
}

//______________________________________________________________________________
NSView<X11Window> *FindViewUnderPointer()
{
   const Util::AutoreleasePool pool;

   if (QuartzWindow *topLevel = FindWindowUnderPointer()) {
      const NSPoint mousePosition = [topLevel mouseLocationOutsideOfEventStream];
      return (NSView<X11Window> *)[[topLevel contentView] hitTest : mousePosition];
   }

   return nil;
}

//______________________________________________________________________________
QuartzWindow *FindWindowForPointerEvent(NSEvent *pointerEvent)
{
   //FindWindowForPointerEvent is required because due to grabs
   //the receiver of the event can be different from the actual
   //window under cursor.

   assert(pointerEvent != nil &&
          "FindWindowForPointerEvent, parameter 'pointerEvent' is nil");

   const Util::AutoreleasePool pool;

   NSArray * const orderedWindows = [NSApp orderedWindows];
   for (NSWindow *nsWindow in orderedWindows) {
      if (![nsWindow isKindOfClass : [QuartzWindow class]])
         continue;

      QuartzWindow * const qWindow = (QuartzWindow *)nsWindow;

      if (qWindow.fIsDeleted)//Because of reference counting this can happen.
         continue;

      //Can it be false and still in this array???
      if (qWindow.fMapState != kIsViewable)
         continue;

      NSPoint mousePosition = [pointerEvent locationInWindow];
      //The event has a window, so position is in this window's coordinate system,
      //convert it into screen point first.
      if ([pointerEvent window]) {
         //convertBaseToScreen is deprecated.
         //mousePosition = [[pointerEvent window] convertBaseToScreen : mousePosition];
         mousePosition = ConvertPointFromBaseToScreen([pointerEvent window], mousePosition);
      }

      //convertScreenToBase is deprecated.
      //mousePosition = [qWindow convertScreenToBase : mousePosition];
      mousePosition = ConvertPointFromScreenToBase(mousePosition, qWindow);

      const NSSize windowSize = qWindow.frame.size;
      if (mousePosition.x >= 0 && mousePosition.x <= windowSize.width &&
          mousePosition.y >= 0 && mousePosition.y <= windowSize.height)
         return qWindow;
   }

   return nil;
}

//______________________________________________________________________________
NSView<X11Window> *FindViewForPointerEvent(NSEvent *pointerEvent)
{
   //FindViewForPointerEvent is required because of grabs - the receiver of the
   //event can be different from the actual window under cursor.

   assert(pointerEvent != nil &&
          "FindViewForPointerEvent, parameter 'pointerEvent' is nil");

   const Util::AutoreleasePool pool;

   if (QuartzWindow *topLevel = FindWindowForPointerEvent(pointerEvent)) {
      NSPoint mousePosition = [pointerEvent locationInWindow];
      if ([pointerEvent window])
         mousePosition = ConvertPointFromBaseToScreen([pointerEvent window], mousePosition);

      //convertScreenToBase is deprecated.
      //mousePosition = [topLevel convertScreenToBase : mousePosition];
      mousePosition = ConvertPointFromScreenToBase(mousePosition, topLevel);

      return (NSView<X11Window> *)[[topLevel contentView] hitTest : mousePosition];
   }

   return nil;
}

#pragma mark - Downscale image ("reading color bits" on retina macs).

//Hoho, we support C++11?? Let's return by value then!!!
std::vector<unsigned char> DownscaledImageData(unsigned w, unsigned h, CGImageRef image)
{
   assert(w != 0 && h != 0 && "DownscaledImageData, invalid geometry");
   assert(image != nullptr && "DonwscaledImageData, invalid parameter 'image'");

   std::vector<unsigned char> result;
   try {
      result.resize(w * h * 4);
   } catch (const std::bad_alloc &) {
      NSLog(@"DownscaledImageData, memory allocation failed");
      return result;
   }

   const Util::CFScopeGuard<CGColorSpaceRef> colorSpace(CGColorSpaceCreateDeviceRGB());//[1]
   if (!colorSpace.Get()) {
      NSLog(@"DownscaledImageData, CGColorSpaceCreateDeviceRGB failed");
      return result;
   }

   Util::CFScopeGuard<CGContextRef> ctx(CGBitmapContextCreateWithData(&result[0], w, h, 8,
                                                                      w * 4, colorSpace.Get(),
                                                                      kCGImageAlphaPremultipliedLast, NULL, 0));
   if (!ctx.Get()) {
      NSLog(@"DownscaledImageData, CGBitmapContextCreateWithData failed");
      return result;
   }

   CGContextDrawImage(ctx.Get(), CGRectMake(0, 0, w, h), image);

   return result;
}

#pragma mark - "Focus management" - just make another window key window.

//______________________________________________________________________________
void WindowLostFocus(Window_t winID)
{
   //XQuartz (and other X11
   if (![NSApp isActive])
      return;

   const Util::AutoreleasePool pool;

   NSArray * const orderedWindows = [NSApp orderedWindows];
   for (NSWindow *nsWindow in orderedWindows) {
      if (![nsWindow isKindOfClass : [QuartzWindow class]])
         continue;

      QuartzWindow * const qWindow = (QuartzWindow *)nsWindow;

      if (qWindow.fIsDeleted || qWindow.fMapState != kIsViewable || qWindow.fID == winID)
         continue;

      if (qWindow.fContentView.fOverrideRedirect)
         continue;

      [qWindow makeKeyAndOrderFront : qWindow];
      break;
   }
}

#pragma mark - 'shape mask' - to create a window with arbitrary (probably non-rectangle) shape.

//______________________________________________________________________________
void ClipToShapeMask(NSView<X11Window> *view, CGContextRef ctx)
{
   assert(view != nil && "ClipToShapeMask, parameter 'view' is nil");
   assert(ctx != 0 && "ClipToShapeMask, parameter 'ctx' is null");

   QuartzWindow * const topLevelParent = view.fQuartzWindow;
   assert(topLevelParent.fShapeCombineMask != nil &&
          "ClipToShapeMask, fShapeCombineMask is nil on a top-level window");
   assert(topLevelParent.fShapeCombineMask.fImage != 0 &&
          "ClipToShapeMask, shape mask is null");

   //Important: shape mask should have the same width and height as
   //a top-level window. In ROOT it does not :( Say hello to visual artifacts.

   //Attach clip mask to the context.
   if (!view.fParentView) {
      //'view' is a top-level view.
      const CGRect clipRect = CGRectMake(0, 0, topLevelParent.fShapeCombineMask.fWidth,
                                               topLevelParent.fShapeCombineMask.fHeight);
      CGContextClipToMask(ctx, clipRect, topLevelParent.fShapeCombineMask.fImage);
   } else {
      NSRect clipRect = view.frame;
      //More complex case: 'self' is a child view, we have to create a subimage from shape mask.
      clipRect.origin = [view.fParentView convertPoint : clipRect.origin
                         toView : [view window].contentView];
      clipRect.origin.y = X11::LocalYROOTToCocoa((NSView<X11Window> *)[view window].contentView,
                                                  clipRect.origin.y + clipRect.size.height);

      if (AdjustCropArea(topLevelParent.fShapeCombineMask, clipRect)) {
         const Util::CFScopeGuard<CGImageRef>
            clipImageGuard(CGImageCreateWithImageInRect(topLevelParent.fShapeCombineMask.fImage,
                                                        NSRectToCGRect(clipRect)));
         clipRect.origin = NSPoint();
         CGContextClipToMask(ctx, NSRectToCGRect(clipRect), clipImageGuard.Get());
      } else {
         //View is invisible.
         CGRect rect = {};
         CGContextClipToRect(ctx, rect);
      }
   }
}

#pragma mark - Window's geometry and attributes.

//______________________________________________________________________________
void SetWindowAttributes(const SetWindowAttributes_t *attr, NSObject<X11Window> *window)
{
   assert(attr != 0 && "SetWindowAttributes, parameter 'attr' is null");
   assert(window != nil && "SetWindowAttributes, parameter 'window' is nil");

   const Mask_t mask = attr->fMask;

   if (mask & kWABackPixel)
      window.fBackgroundPixel = attr->fBackgroundPixel;

   if (mask & kWAEventMask)
      window.fEventMask = attr->fEventMask;

   if (mask & kWABitGravity)
      window.fBitGravity = attr->fBitGravity;

   if (mask & kWAWinGravity)
      window.fWinGravity = attr->fWinGravity;

   if (mask & kWAOverrideRedirect) {
      if ([(NSObject *)window isKindOfClass : [QuartzWindow class]]) {
         QuartzWindow * const qw = (QuartzWindow *)window;
         [qw setStyleMask : Details::kBorderlessWindowMask];
         [qw setAlphaValue : 0.95];
      }

      window.fOverrideRedirect = YES;
   }
}

//______________________________________________________________________________
void GetWindowGeometry(NSObject<X11Window> *win, WindowAttributes_t *dst)
{
   assert(win != nil && "GetWindowGeometry, parameter 'win' is nil");
   assert(dst != 0 && "GetWindowGeometry, parameter 'dst' is null");

   dst->fX = win.fX;
   dst->fY = win.fY;

   dst->fWidth = win.fWidth;
   dst->fHeight = win.fHeight;
}

//______________________________________________________________________________
void GetWindowAttributes(NSObject<X11Window> *window, WindowAttributes_t *dst)
{
   assert(window != nil && "GetWindowAttributes, parameter 'window' is nil");
   assert(dst != 0 && "GetWindowAttributes, parameter 'attr' is null");

   *dst = WindowAttributes_t();

   //fX, fY, fWidth, fHeight.
   GetWindowGeometry(window, dst);

   //Actually, most of them are not used by GUI.
   dst->fBorderWidth = 0;
   dst->fDepth = window.fDepth;
   //Dummy value.
   dst->fVisual = 0;
   //Dummy value.
   dst->fRoot = 0;
   dst->fClass = window.fClass;
   dst->fBitGravity = window.fBitGravity;
   dst->fWinGravity = window.fWinGravity;
   //Dummy value.
   dst->fBackingStore = kAlways;//??? CHECK
   dst->fBackingPlanes = 0;

   //Dummy value.
   dst->fBackingPixel = 0;

   dst->fSaveUnder = 0;

   //Dummy value.
   dst->fColormap = 0;
   //Dummy value.
   dst->fMapInstalled = kTRUE;

   dst->fMapState = window.fMapState;

   dst->fAllEventMasks = window.fEventMask;
   dst->fYourEventMask = window.fEventMask;

   //Not used by GUI.
   //dst->fDoNotPropagateMask

   dst->fOverrideRedirect = window.fOverrideRedirect;
   //Dummy value.
   dst->fScreen = 0;
}

//With Apple's poor man's objective C/C++ + "brilliant" Cocoa you never know, what should be
//the linkage of callback functions, API + language dialects == MESS. I declare/define this comparators here
//as having "C++" linkage. If one good day clang will start to complane, I'll have to change this.

#pragma mark - Comparators (I need them when changing a window's z-order).

//______________________________________________________________________________
// SDK 10.11 and above ...
#ifdef MAC_OS_X_VERSION_10_11
NSComparisonResult CompareViewsToLower(__kindof NSView *view1, __kindof NSView *view2, void *context)
#else
NSComparisonResult CompareViewsToLower(id view1, id view2, void *context)
#endif
{
    id topView = (id)context;
    if (view1 == topView)
        return NSOrderedAscending;
    if (view2 == topView)
        return NSOrderedDescending;

    return NSOrderedSame;
}

//______________________________________________________________________________
// SDK 10.11 and above ...
#ifdef MAC_OS_X_VERSION_10_11
NSComparisonResult CompareViewsToRaise(__kindof NSView *view1, __kindof NSView *view2, void *context)
#else
NSComparisonResult CompareViewsToRaise(id view1, id view2, void *context)
#endif
{
   id topView = (id)context;
   if (view1 == topView)
      return NSOrderedDescending;
   if (view2 == topView)
      return NSOrderedAscending;

   return NSOrderedSame;
}

#pragma mark - Cursor's area.

//______________________________________________________________________________
NSPoint GetCursorHotStop(NSImage *image, ECursor cursor)
{
   assert(image != nil && "CursroHotSpot, parameter 'image' is nil");

   const NSSize imageSize = image.size;

   if (cursor == kArrowRight)
      return NSMakePoint(imageSize.width, imageSize.height / 2);

   return NSMakePoint(imageSize.width / 2, imageSize.height / 2);
}

//______________________________________________________________________________
NSCursor *CreateCustomCursor(ECursor currentCursor)
{
   // Returns auto-released cursor object.
   const char *pngFileName = 0;

   switch (currentCursor) {
   case kMove:
      pngFileName = "move_cursor.png";
      break;
   case kArrowHor:
      pngFileName = "hor_arrow_cursor.png";
      break;
   case kArrowVer:
      pngFileName = "ver_arrow_cursor.png";
      break;
   case kArrowRight:
      pngFileName = "right_arrow_cursor.png";
      break;
   case kRotate:
      pngFileName = "rotate.png";
      break;
   case kBottomLeft:
   case kTopRight:
      pngFileName = "top_right_cursor.png";
      break;
   case kTopLeft:
   case kBottomRight:
      pngFileName = "top_left_cursor.png";
      break;
   default:;
   }

   if (pngFileName) {
      const char * const path = gSystem->Which(TROOT::GetIconPath(), pngFileName, kReadPermission);
      const Util::ScopedArray<const char> arrayGuard(path);

      if (!path || path[0] == 0) {
         //File was not found.
         return nil;
      }

      NSString *nsPath = [NSString stringWithFormat : @"%s", path];//in autorelease pool.
      NSImage * const cursorImage = [[NSImage alloc] initWithContentsOfFile : nsPath];

      if (!cursorImage)
         return nil;

      const NSPoint hotSpot(X11::GetCursorHotStop(cursorImage, currentCursor));
      NSCursor * const customCursor = [[[NSCursor alloc] initWithImage : cursorImage
                                                         hotSpot : hotSpot] autorelease];

      [cursorImage release];

      return customCursor;
   }

   return nil;
}

//______________________________________________________________________________
NSCursor *CreateCursor(ECursor currentCursor)
{
   // Returns auto-released cursor object.

   //Cursors from TVirtaulX:
   // kBottomLeft, kBottomRight, kTopLeft,  kTopRight,
   // kBottomSide, kLeftSide,    kTopSide,  kRightSide,
   // kMove,       kCross,       kArrowHor, kArrowVer,
   // kHand,       kRotate,      kPointer,  kArrowRight,
   // kCaret,      kWatch

   NSCursor *cursor = nil;
   switch (currentCursor) {
   case kCross:
      cursor = [NSCursor crosshairCursor];
      break;
   case kPointer:
      cursor = [NSCursor arrowCursor];
      break;
   case kHand:
      cursor = [NSCursor openHandCursor];
      break;
   case kLeftSide:
      cursor = [NSCursor resizeLeftCursor];
      break;
   case kRightSide:
      cursor = [NSCursor resizeRightCursor];
      break;
   case kTopSide:
      cursor = [NSCursor resizeUpCursor];
      break;
   case kBottomSide:
      cursor = [NSCursor resizeDownCursor];
      break;
   case kCaret:
      cursor = [NSCursor IBeamCursor];
      break;
   case kRotate:
   case kWatch:
   default:
      cursor = CreateCustomCursor(currentCursor);
   }

   return cursor;
}

//TGTextView/TGHtml is a very special window: it's a TGCompositeFrame,
//which has TGCompositeFrame inside (TGViewFrame). This TGViewFrame
//delegates Expose events to its parent, and parent tries to draw
//inside a TGViewFrame. This does not work with default
//QuartzView -drawRect/TGCocoa. So I need a trick to identify
//this special window.

#pragma mark - Workarounds for a text view and its descendants.

//______________________________________________________________________________
bool ViewIsTextView(unsigned viewID)
{
   const TGWindow * const window = gClient->GetWindowById(viewID);
   if (!window)
      return false;
   // This code used to use TObject::InheritsFrom, however since this is
   // run under the AppKit, we can not call core/meta functions, otherwise
   // we will run into deadlocks.
   return dynamic_cast<const TGTextView*>(window);
}

//______________________________________________________________________________
bool ViewIsTextView(NSView<X11Window> *view)
{
   assert(view != nil && "ViewIsTextView, parameter 'view' is nil");

   return ViewIsTextView(view.fID);
}

//______________________________________________________________________________
bool ViewIsTextViewFrame(NSView<X11Window> *view, bool checkParent)
{
   assert(view != nil && "ViewIsTextViewFrame, parameter 'view' is nil");

   const TGWindow * const window = gClient->GetWindowById(view.fID);
   if (!window)
      return false;

   // This code used to use TObject::InheritsFrom, however since this is
   // run under the AppKit, we can not call core/meta functions, otherwise
   // we will run into deadlocks.
   if (!dynamic_cast<const TGViewFrame*>(window))
      return false;

   if (!checkParent)
      return true;

   if (!view.fParentView)
      return false;

   return ViewIsTextView(view.fParentView);
}

//______________________________________________________________________________
bool ViewIsHtmlView(unsigned viewID)
{
   const TGWindow * const window = gClient->GetWindowById(viewID);
   if (!window)
      return false;
   // This code used to use TObject::InheritsFrom, however since this is
   // run under the AppKit, we can not call core/meta functions, otherwise
   // we will run into deadlocks.
   return window->TestBit(TGWindow::kIsHtmlView);
}

//______________________________________________________________________________
bool ViewIsHtmlView(NSView<X11Window> *view)
{
   assert(view != nil && "ViewIsHtmlView, parameter 'view' is nil");

   return ViewIsHtmlView(view.fID);
}

//______________________________________________________________________________
bool ViewIsHtmlViewFrame(NSView<X11Window> *view, bool checkParent)
{
   //
   assert(view != nil && "ViewIsHtmlViewFrame, parameter 'view' is nil");

   const TGWindow * const window = gClient->GetWindowById(view.fID);
   if (!window)
      return false;

   // This code used to use TObject::InheritsFrom, however since this is
   // run under the AppKit, we can not call core/meta functions, otherwise
   // we will run into deadlocks.
   if (!dynamic_cast<const TGViewFrame*>(window))
      return false;

   if (!checkParent)
      return true;

   if (!view.fParentView)
      return false;

   return ViewIsHtmlView(view.fParentView);
}

//______________________________________________________________________________
NSView<X11Window> *FrameForTextView(NSView<X11Window> *textView)
{
   assert(textView != nil && "FrameForTextView, parameter 'textView' is nil");

   for (NSView<X11Window> *child in [textView subviews]) {
      if (ViewIsTextViewFrame(child, false))
         return child;
   }

   return nil;
}

//______________________________________________________________________________
NSView<X11Window> *FrameForHtmlView(NSView<X11Window> *htmlView)
{
   assert(htmlView != nil && "FrameForHtmlView, parameter 'htmlView' is nil");

   for (NSView<X11Window> *child in [htmlView subviews]) {
      if (ViewIsHtmlViewFrame(child, false))
         return child;
   }

   return nil;
}

#pragma mark - Workarounds for 'paint out of paint events'.

//______________________________________________________________________________
bool LockFocus(NSView<X11Window> *view)
{
   assert(view != nil && "LockFocus, parameter 'view' is nil");
   assert([view isKindOfClass : [QuartzView class]] &&
          "LockFocus, QuartzView is expected");

   if ([view lockFocusIfCanDraw]) {
      NSGraphicsContext *nsContext = [NSGraphicsContext currentContext];
      assert(nsContext != nil && "LockFocus, currentContext is nil");
      CGContextRef currContext = (CGContextRef)[nsContext graphicsPort];
      assert(currContext != 0 && "LockFocus, graphicsPort is null");//remove this assert?

      ((QuartzView *)view).fContext = currContext;

      return true;
   }

   return false;
}

//______________________________________________________________________________
void UnlockFocus(NSView<X11Window> *view)
{
   assert(view != nil && "UnlockFocus, parameter 'view' is nil");
   assert([view isKindOfClass : [QuartzView class]] &&
          "UnlockFocus, QuartzView is expected");

   [view unlockFocus];
   ((QuartzView *)view).fContext = 0;
}

}//X11
}//MacOSX
}//ROOT

namespace Quartz = ROOT::Quartz;
namespace Util = ROOT::MacOSX::Util;
namespace X11 = ROOT::MacOSX::X11;
namespace Details = ROOT::MacOSX::Details;

#ifdef DEBUG_ROOT_COCOA

#pragma mark - 'loggers'.

namespace {

//______________________________________________________________________________
void log_attributes(const SetWindowAttributes_t *attr, unsigned winID)
{
   //This function is loggin requests, at the moment I can not set all
   //of these attributes, so I first have to check, what is actually
   //requested by ROOT.
   static std::ofstream logfile("win_attr.txt");

   const Mask_t mask = attr->fMask;
   if (mask & kWABackPixmap)
      logfile<<"win "<<winID<<": BackPixmap\n";
   if (mask & kWABackPixel)
      logfile<<"win "<<winID<<": BackPixel\n";
   if (mask & kWABorderPixmap)
      logfile<<"win "<<winID<<": BorderPixmap\n";
   if (mask & kWABorderPixel)
      logfile<<"win "<<winID<<": BorderPixel\n";
   if (mask & kWABorderWidth)
      logfile<<"win "<<winID<<": BorderWidth\n";
   if (mask & kWABitGravity)
      logfile<<"win "<<winID<<": BitGravity\n";
   if (mask & kWAWinGravity)
      logfile<<"win "<<winID<<": WinGravity\n";
   if (mask & kWABackingStore)
      logfile<<"win "<<winID<<": BackingStore\n";
   if (mask & kWABackingPlanes)
      logfile<<"win "<<winID<<": BackingPlanes\n";
   if (mask & kWABackingPixel)
      logfile<<"win "<<winID<<": BackingPixel\n";
   if (mask & kWAOverrideRedirect)
      logfile<<"win "<<winID<<": OverrideRedirect\n";
   if (mask & kWASaveUnder)
      logfile<<"win "<<winID<<": SaveUnder\n";
   if (mask & kWAEventMask)
      logfile<<"win "<<winID<<": EventMask\n";
   if (mask & kWADontPropagate)
      logfile<<"win "<<winID<<": DontPropagate\n";
   if (mask & kWAColormap)
      logfile<<"win "<<winID<<": Colormap\n";
   if (mask & kWACursor)
      logfile<<"win "<<winID<<": Cursor\n";
}

//______________________________________________________________________________
void print_mask_info(ULong_t mask)
{
   if (mask & kButtonPressMask)
      NSLog(@"button press mask");
   if (mask & kButtonReleaseMask)
      NSLog(@"button release mask");
   if (mask & kExposureMask)
      NSLog(@"exposure mask");
   if (mask & kPointerMotionMask)
      NSLog(@"pointer motion mask");
   if (mask & kButtonMotionMask)
      NSLog(@"button motion mask");
   if (mask & kEnterWindowMask)
      NSLog(@"enter notify mask");
   if (mask & kLeaveWindowMask)
      NSLog(@"leave notify mask");
}

}
#endif

@implementation CrosshairView

- (void) drawRect : (NSRect) dirtyRect
{
    [super drawRect:dirtyRect];

    NSGraphicsContext *nsContext = [NSGraphicsContext currentContext];
    if (!nsContext)
        return;

    CGContextRef cgContext = nsContext.CGContext;
    if (!cgContext)
        return;

    const Quartz::CGStateGuard ctxGuard(cgContext);

    CGContextSetRGBStrokeColor(cgContext, 0., 0., 0., 1.);
    CGContextSetLineWidth(cgContext, 1.);
    // Horizontal line:
    CGContextBeginPath(cgContext);
    CGContextMoveToPoint(cgContext, self.start1.x, self.start1.y);
    CGContextAddLineToPoint(cgContext, self.end1.x, self.end1.y);
    CGContextStrokePath(cgContext);
    // Vertical line:
    CGContextBeginPath(cgContext);
    CGContextMoveToPoint(cgContext, self.start2.x, self.start2.y);
    CGContextAddLineToPoint(cgContext, self.end2.x, self.end2.y);
    CGContextStrokePath(cgContext);
}

@end

@implementation CrosshairWindow

- (instancetype) init
{
    if (self = [super init])
    {
        self.styleMask = NSWindowStyleMaskBorderless; // No titlebar, buttons, etc.
        self.opaque = NO;
        self.hasShadow = NO;
        self.backgroundColor = NSColor.clearColor; // No background.
        self.ignoresMouseEvents = YES; // Lets mouse events pass through.
        self.contentView = [[CrosshairView alloc] init];
    }
    return self;
}

#pragma mark - suppress the normal window behavior.
- (BOOL)canBecomeKeyWindow
{
    return NO;
}
- (BOOL)canBecomeMainWindow
{
    return NO;
}

@end


@implementation QuartzWindow

@synthesize fMainWindow;
@synthesize fHasFocus;

#pragma mark - QuartzWindow's life cycle.

//______________________________________________________________________________
- (id) initWithContentRect : (NSRect) contentRect styleMask : (NSUInteger) windowStyle
       backing : (NSBackingStoreType) bufferingType defer : (BOOL) deferCreation
       windowAttributes : (const SetWindowAttributes_t *) attr
{
   self = [super initWithContentRect : contentRect styleMask : windowStyle
                             backing : bufferingType defer : deferCreation];

   if (self) {
      //ROOT's not able to draw GUI concurrently, thanks to global variables and gVirtualX itself.
      [self setAllowsConcurrentViewDrawing : NO];

      self.delegate = self;
      //create content view here.
      NSRect contentViewRect = contentRect;
      contentViewRect.origin.x = 0.f;
      contentViewRect.origin.y = 0.f;

      fContentView = [[QuartzView alloc] initWithFrame : contentViewRect windowAttributes : 0];

      [self setContentView : fContentView];

      [fContentView release];
      fDelayedTransient = NO;

      if (attr)
         X11::SetWindowAttributes(attr, self);

      fIsDeleted = NO;
      fHasFocus = NO;
   }

   return self;
}

//______________________________________________________________________________
- (id) initWithGLView : (ROOTOpenGLView *) glView
{
   using namespace Details;

   assert(glView != nil && "-initWithGLView, parameter 'glView' is nil");

   const NSUInteger styleMask = kTitledWindowMask | kClosableWindowMask |
                                kMiniaturizableWindowMask | kResizableWindowMask;

   NSRect contentRect = glView.frame;
   contentRect.origin = NSPoint();

   self = [super initWithContentRect : contentRect styleMask : styleMask
                             backing : NSBackingStoreBuffered defer : NO];

   if (self) {
      //ROOT's not able to draw GUI concurrently, thanks to global variables and gVirtualX itself.
      [self setAllowsConcurrentViewDrawing : NO];
      self.delegate = self;
      fContentView = glView;
      [self setContentView : fContentView];
      fDelayedTransient = NO;
      fIsDeleted = NO;
      fHasFocus = NO;
   }

   return self;
}

//______________________________________________________________________________
- (void) dealloc
{
   [fShapeCombineMask release];
   [super dealloc];
}

//______________________________________________________________________________
- (BOOL) fIsDeleted
{
   return fIsDeleted;
}

//______________________________________________________________________________
- (void) setContentView:(NSView *)cv
{
    [super setContentView:cv];
    if ([cv isKindOfClass:[QuartzView class]])
        fContentView = (QuartzView *)cv;
    else
        fContentView = nil;
}

//______________________________________________________________________________
- (void) setFIsDeleted : (BOOL) deleted
{
   fIsDeleted = deleted;
}

#pragma mark - Forwaring: I want to forward a lot of property setters/getters to the content view.

//______________________________________________________________________________
- (void) forwardInvocation : (NSInvocation *) anInvocation
{
   if (!fContentView)
      return;

   if ([fContentView respondsToSelector : [anInvocation selector]]) {
      [anInvocation invokeWithTarget : fContentView];
   } else {
      [super forwardInvocation : anInvocation];
   }
}

//______________________________________________________________________________
- (NSMethodSignature*) methodSignatureForSelector : (SEL) selector
{
   NSMethodSignature *signature = [super methodSignatureForSelector : selector];

   if (!signature) {
      if (fContentView)
         signature = [fContentView methodSignatureForSelector : selector];
   }

   return signature;
}

//______________________________________________________________________________
- (void) addTransientWindow : (QuartzWindow *) window
{
   //Transient window: all the popups (menus, dialogs, popups, comboboxes, etc.)
   //always should be on the top of its 'parent' window.
   //To enforce this ordering, I have to connect such windows with parent/child
   //relation (it's not the same as a view hierarchy - both child and parent
   //windows are top-level windows).

   assert(window != nil && "-addTransientWindow:, parameter 'window' is nil");

   window.fMainWindow = self;

   if (window.fMapState != kIsViewable) {
      //If I add it as child, it'll immediately make a window visible
      //and this thing sucks.
      window.fDelayedTransient = YES;
   } else {
      [self addChildWindow : window ordered : NSWindowAbove];
      window.fDelayedTransient = NO;
   }
}

//______________________________________________________________________________
- (void) makeKeyAndOrderFront : (id) sender
{
#pragma unused(sender)

   //The more I know Cocoa, the less I like it.
   //Window behavior between spaces is a total mess.
   //Set the window to join all spaces.
#ifdef MAC_OS_X_VERSION_10_9
   [self setCollectionBehavior : NSWindowCollectionBehaviorMoveToActiveSpace];
#else
   [self setCollectionBehavior : NSWindowCollectionBehaviorCanJoinAllSpaces];
#endif
   //now bring it to the front, it will appear on the active space.
   [super makeKeyAndOrderFront : self];
   //then reset the collection behavior to default, so the window
   [self setCollectionBehavior : NSWindowCollectionBehaviorDefault];
}

//______________________________________________________________________________
- (void) setFDelayedTransient : (BOOL) d
{
   fDelayedTransient = d;
}

//______________________________________________________________________________
- (QuartzImage *) fShapeCombineMask
{
   return fShapeCombineMask;
}

//______________________________________________________________________________
- (void) setFShapeCombineMask : (QuartzImage *) mask
{
   if (mask != fShapeCombineMask) {
      [fShapeCombineMask release];
      if (mask) {
         fShapeCombineMask = [mask retain];
         //TODO: Check window's shadow???
      }
   }
}

#pragma mark - X11Drawable's protocol.

//______________________________________________________________________________
- (BOOL) fIsPixmap
{
   //Never.
   return NO;
}

//______________________________________________________________________________
- (BOOL) fIsOpenGLWidget
{
   //Never.
   return NO;
}

//______________________________________________________________________________
- (CGFloat) fScaleFactor
{
   if (!self.screen)
      return 1.;
   return self.screen.backingScaleFactor;
}

//______________________________________________________________________________
- (int) fX
{
   return X11::GlobalXCocoaToROOT(self.frame.origin.x);
}

//______________________________________________________________________________
- (int) fY
{
   return X11::GlobalYCocoaToROOT(self.frame.origin.y + self.frame.size.height);
}

//______________________________________________________________________________
- (unsigned) fWidth
{
   return self.frame.size.width;
}

//______________________________________________________________________________
- (unsigned) fHeight
{
   //NSWindow's frame (height component) also includes title-bar.
   //So I have to use content view's height.
   //Obviously, there is a "hole" == 22 pixels.
   assert(fContentView != nil && "-fHeight:, content view is nil");

   return fContentView.frame.size.height;
}

//______________________________________________________________________________
- (void) setDrawableSize : (NSSize) newSize
{
   //Can not simply do self.frame.size = newSize.
   assert(!(newSize.width < 0) && "-setDrawableSize:, width is negative");
   assert(!(newSize.height < 0) && "-setDrawableSize:, height is negative");

   NSRect frame = self.frame;
   //dY is potentially a titlebar height.
   const CGFloat dY = fContentView ? frame.size.height - fContentView.frame.size.height : 0.;
   //Adjust the frame.
   frame.origin.y = frame.origin.y + frame.size.height - newSize.height - dY;
   frame.size = newSize;
   frame.size.height += dY;
   [self setFrame : frame display : YES];
}

//______________________________________________________________________________
- (void) setX : (int) x Y : (int) y width : (unsigned) w height : (unsigned) h
{
   NSSize newSize = {};
   newSize.width = w;
   newSize.height = h;
   [self setContentSize : newSize];

   //Check how this is affected by title bar's height.
   NSPoint topLeft = {};
   topLeft.x = X11::GlobalXROOTToCocoa(x);
   topLeft.y = X11::GlobalYROOTToCocoa(y);

   [self setFrameTopLeftPoint : topLeft];
}

//______________________________________________________________________________
- (void) setX : (int) x Y : (int) y
{
   NSPoint topLeft = {};
   topLeft.x = X11::GlobalXROOTToCocoa(x);
   topLeft.y = X11::GlobalYROOTToCocoa(y);

   [self setFrameTopLeftPoint : topLeft];
}

//______________________________________________________________________________
- (void) copy : (NSObject<X11Drawable> *) src area : (X11::Rectangle) area withMask : (QuartzImage *) mask
         clipOrigin : (X11::Point) clipXY toPoint : (X11::Point) dstPoint
{
   if (!fContentView)
      return;

   [fContentView copy : src area : area withMask : mask clipOrigin : clipXY toPoint : dstPoint];
}

//______________________________________________________________________________
- (unsigned char *) readColorBits : (X11::Rectangle) area
{
   if (!fContentView)
      return nullptr;

   return [fContentView readColorBits : area];
}

#pragma mark - CrosshairWindow/View

//______________________________________________________________________________
- (void) addCrosshairWindow
{
    if ([self findCrosshairWindow])
        return;

    CrosshairWindow *special = [[CrosshairWindow alloc] init];
    [self adjustCrosshairWindowGeometry:special];
    [self addChildWindow : special ordered : NSWindowAbove];
    [special release];
}

//______________________________________________________________________________
- (void) adjustCrosshairWindowGeometry
{
    if (auto win = [self findCrosshairWindow])
        [self adjustCrosshairWindowGeometry:win];
}

//______________________________________________________________________________
- (void) adjustCrosshairWindowGeometry : (CrosshairWindow *) win
{
    assert(win && "invalid (nil) parameter 'win'");
    auto frame = self.contentView.frame;
    frame = [self convertRectToScreen:frame];
    [win setFrame:frame display:NO];
}

//______________________________________________________________________________
- (void) removeCrosshairWindow
{
    if (auto win = [self findCrosshairWindow]) {
        // For some reason, without ordeing out, the crosshair window's content stays
        // in the parent's window. Thus we first have to order out the crosshair window.
        [win orderOut:nil];
        [self removeChildWindow : win];
    }
}

//______________________________________________________________________________
- (CrosshairWindow *) findCrosshairWindow
{
    auto children = [self childWindows];
    for (NSWindow *child in children) {
        if ([child isKindOfClass : CrosshairWindow.class])
            return (CrosshairWindow *)child;
    }
    return nil;
}


#pragma mark - X11Window protocol's implementation.

//______________________________________________________________________________
- (QuartzView *) fParentView
{
   return nil;
}

//______________________________________________________________________________
- (void) setFParentView : (QuartzView *) parent
{
#pragma unused(parent)
}

//______________________________________________________________________________
- (NSView<X11Window> *) fContentView
{
   return fContentView;
}

//______________________________________________________________________________
- (QuartzWindow *) fQuartzWindow
{
   return self;
}

//... many forwards to fContentView.

//______________________________________________________________________________
- (void) setFBackgroundPixel : (unsigned long) backgroundColor
{
   if (!fContentView)
      return;

   if (!fShapeCombineMask) {
      CGFloat rgba[] = {0., 0., 0., 1.};
      X11::PixelToRGB(backgroundColor, rgba);

      [self setBackgroundColor : [NSColor colorWithColorSpace : [NSColorSpace deviceRGBColorSpace] components : rgba count : 4]];
   }

   fContentView.fBackgroundPixel = backgroundColor;
}

//______________________________________________________________________________
- (unsigned long) fBackgroundPixel
{
   if (!fContentView)
       return 0;

   return fContentView.fBackgroundPixel;
}

//______________________________________________________________________________
- (int) fMapState
{
   //Top-level window can be only kIsViewable or kIsUnmapped (not unviewable).
   if (!fContentView)
      return kIsUnmapped;

   if ([fContentView isHidden])
      return kIsUnmapped;

   return kIsViewable;
}

//______________________________________________________________________________
- (void) addChild : (NSView<X11Window> *) child
{
   assert(child != nil && "-addChild:, parameter 'child' is nil");

   if (!fContentView) {
      //This can happen only in case of re-parent operation.
      assert([child isKindOfClass : [QuartzView class]] &&
             "-addChild: gl view in a top-level window as content view is not supported");

      fContentView = (QuartzView *)child;
      [self setContentView : child];
      fContentView.fParentView = nil;
   } else
      [fContentView addChild : child];
}

//______________________________________________________________________________
- (void) getAttributes : (WindowAttributes_t *) attr
{
   if (!fContentView)
      return;

   assert(attr && "-getAttributes:, parameter 'attr' is nil");

   X11::GetWindowAttributes(self, attr);
}

//______________________________________________________________________________
- (void) setAttributes : (const SetWindowAttributes_t *) attr
{
   assert(attr != 0 && "-setAttributes:, parameter 'attr' is null");

#ifdef DEBUG_ROOT_COCOA
   log_attributes(attr, self.fID);
#endif

   X11::SetWindowAttributes(attr, self);
}

//______________________________________________________________________________
- (void) mapRaised
{
   if (!fContentView)
      return;

   const Util::AutoreleasePool pool;

   [fContentView setHidden : NO];
   [self makeKeyAndOrderFront : self];
   [fContentView configureNotifyTree];

   if (fDelayedTransient) {
      fDelayedTransient = NO;
      [fMainWindow addChildWindow : self ordered : NSWindowAbove];
   }
}

//______________________________________________________________________________
- (void) mapWindow
{
   if (!fContentView)
      return;

   const Util::AutoreleasePool pool;

   [fContentView setHidden : NO];
   [self makeKeyAndOrderFront : self];
   [fContentView configureNotifyTree];

   if (fDelayedTransient) {
      fDelayedTransient = NO;
      [fMainWindow addChildWindow : self ordered : NSWindowAbove];
   }
}

//______________________________________________________________________________
- (void) mapSubwindows
{
   if (!fContentView)
      return;

   const Util::AutoreleasePool pool;

   [fContentView mapSubwindows];
   [fContentView configureNotifyTree];
}

//______________________________________________________________________________
- (void) unmapWindow
{
   if (!fContentView)
      return;

   [fContentView setHidden : YES];
   [self orderOut : self];

   if (fMainWindow && !fDelayedTransient) {
      [fMainWindow removeChildWindow : self];
      fMainWindow = nil;
   }
}

#pragma mark - Events.

//______________________________________________________________________________
- (void) sendEvent : (NSEvent *) theEvent
{
   //With XQuartz, if you open a menu and try to move a window without closing this menu,
   //window does not move, menu closes, and after that you can start draggin a window again.
   //With Cocoa I can not do such a thing (window WILL move), but still can report button release event
   //to close a menu.
   if (!fContentView)
      return;

   if (theEvent.type == Details::kLeftMouseDown || theEvent.type == Details::kRightMouseDown) {
      bool generateFakeRelease = false;

      const NSPoint windowPoint = [theEvent locationInWindow];

      if (windowPoint.x <= 4 || windowPoint.x >= self.fWidth - 4)
         generateFakeRelease = true;

      if (windowPoint.y <= 4 || windowPoint.y >= self.fHeight - 4)
         generateFakeRelease = true;

      const NSPoint viewPoint =  [fContentView convertPoint : windowPoint fromView : nil];

      if (viewPoint.y <= 0 && windowPoint.y >= 0)
         generateFakeRelease = true;

      assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
             "-sendEvent:, gVirtualX is either null or not of TGCocoa type");

      TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
      if (vx->GetEventTranslator()->HasPointerGrab() && generateFakeRelease) {
          vx->GetEventTranslator()->GenerateButtonReleaseEvent(fContentView, theEvent,
                                                               theEvent.type == Details::kLeftMouseDown ?
                                                               kButton1 : kButton3);
         //Yes, ignore this event completely (this means, you are not able to immediately start
         //resizing a window, if some popup is open. Actually, this is more or less
         //the same as with XQuartz and X11 version.
         return;
      }
   }

   [super sendEvent : theEvent];
}

#pragma mark - NSWindowDelegate's methods.

//______________________________________________________________________________
- (BOOL) windowShouldClose : (id) sender
{
#pragma unused(sender)
   if (!fContentView)
      return NO;

   if ([[self childWindows] count])
      return NO;

   //Prepare client message for a window.
   Event_t closeEvent = {};
   closeEvent.fWindow = fContentView.fID;
   closeEvent.fType = kClientMessage;
   closeEvent.fFormat = 32;//Taken from GUI classes.
   closeEvent.fHandle = TGCocoa::fgDeleteWindowAtom;
   closeEvent.fUser[0] = TGCocoa::fgDeleteWindowAtom;
   //Place it into the queue.
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "-windowShouldClose:, gVirtualX is either null or has a type different from TGCocoa");
   ((TGCocoa *)gVirtualX)->SendEvent(fContentView.fID, &closeEvent);

   //Do not let AppKit to close a window,
   //ROOT will do.
   return NO;
}

//______________________________________________________________________________
- (void) windowDidBecomeKey : (NSNotification *) aNotification
{
#pragma unused(aNotification)

   if (!fContentView)
      return;

   if (!fContentView.fOverrideRedirect) {
      fHasFocus = YES;
      //
      assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
             "-windowDidBecomeKey:, gVirtualX is null or not of TGCocoa type");
      TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
      vx->GetEventTranslator()->GenerateFocusChangeEvent(self.fContentView);
   }
}


//______________________________________________________________________________
- (void) windowDidResignKey : (NSNotification *) aNotification
{
#pragma unused(aNotification)
   fHasFocus = NO;
}

@end

#pragma mark - Passive key grab info.

@implementation PassiveKeyGrab

//______________________________________________________________________________
- (id) initWithKey : (unichar) keyCode modifiers : (NSUInteger) modifiers
{
   if (self = [super init]) {
      fKeyCode = keyCode;
      fModifiers = modifiers;
   }

   return self;
}

//______________________________________________________________________________
- (BOOL) matchKey : (unichar) keyCode modifiers : (NSUInteger) modifiers
{
   return keyCode == fKeyCode && modifiers == fModifiers;
}

//______________________________________________________________________________
- (BOOL) matchKey : (unichar) keyCode
{
   return keyCode == fKeyCode;
}

//______________________________________________________________________________
- (unichar) fKeyCode
{
   return fKeyCode;
}

//______________________________________________________________________________
- (NSUInteger) fModifiers
{
   return fModifiers;
}

@end

#pragma mark - X11 property emulation.

@interface QuartzWindowProperty : NSObject {
   NSData *fPropertyData;
   Atom_t fType;
   unsigned fFormat;
}

@property (nonatomic, readonly) Atom_t fType;

@end

@implementation QuartzWindowProperty

@synthesize fType;

//______________________________________________________________________________
- (id) initWithData : (unsigned char *) data size : (unsigned) dataSize type : (Atom_t) type format : (unsigned) format
{
   if (self = [super init]) {
      //Memory is zero-initialized, but just to make it explicit:
      fPropertyData = nil;
      fType = 0;
      fFormat = 0;

      [self resetPropertyData : data size : dataSize type : type format : format];
   }

   return self;
}

//______________________________________________________________________________
- (void) dealloc
{
   [fPropertyData release];

   [super dealloc];
}

//______________________________________________________________________________
- (void) resetPropertyData : (unsigned char *) data size : (unsigned) dataSize
         type : (Atom_t) type format : (unsigned) format
{
   [fPropertyData release];

   fFormat = format;
   if (format == 16)
      dataSize *= 2;
   else if (format == 32)
      dataSize *= 4;

   fPropertyData = [[NSData dataWithBytes : data length : dataSize] retain];

   fType = type;
}

//______________________________________________________________________________
- (NSData *) fPropertyData
{
   return fPropertyData;
}

//______________________________________________________________________________
- (unsigned) fFormat
{
   return fFormat;
}

@end

#pragma mark - QuartzView.

//
//QuartzView is a children view (also is a content view for a top-level QuartzWindow).
//

@implementation QuartzView

@synthesize fID;
@synthesize fContext;
/////////////////////
//SetWindowAttributes_t/WindowAttributes_t
@synthesize fEventMask;
@synthesize fClass;
@synthesize fDepth;
@synthesize fBitGravity;
@synthesize fWinGravity;
@synthesize fBackgroundPixel;
@synthesize fOverrideRedirect;
//SetWindowAttributes_t/WindowAttributes_t
/////////////////////
@synthesize fHasFocus;
@synthesize fParentView;

@synthesize fPassiveGrabButton;
@synthesize fPassiveGrabEventMask;
@synthesize fPassiveGrabKeyModifiers;
@synthesize fActiveGrabEventMask;
@synthesize fPassiveGrabOwnerEvents;
@synthesize fSnapshotDraw;
@synthesize fCurrentCursor;
@synthesize fIsDNDAware;

#pragma mark - Lifetime.

//______________________________________________________________________________
- (id) initWithFrame : (NSRect) frame windowAttributes : (const SetWindowAttributes_t *)attr
{
   if (self = [super initWithFrame : frame]) {
      //Make this explicit (though memory is zero initialized).
      fBackBuffer = nil;
      fID = 0;

      //Passive grab parameters.
      fPassiveGrabButton = -1;//0 is kAnyButton.
      fPassiveGrabEventMask = 0;
      fPassiveGrabOwnerEvents = NO;

      fPassiveKeyGrabs = [[NSMutableArray alloc] init];

      [self setCanDrawConcurrently : NO];

      [self setHidden : YES];
      //Actually, check if view need this.
      //
      if (attr)
         X11::SetWindowAttributes(attr, self);

      fCurrentCursor = kPointer;
      fX11Properties = [[NSMutableDictionary alloc] init];

      fCurrentGrabType = X11::kPGNoGrab;
      fActiveGrabEventMask = 0;
      fActiveGrabOwnerEvents = YES;
   }

   return self;
}

//______________________________________________________________________________
- (void) dealloc
{
   [fBackBuffer release];
   [fPassiveKeyGrabs release];
   [fX11Properties release];
   [fBackgroundPixmap release];
   [super dealloc];
}

#pragma mark - Tracking area.

//Tracking area is required to ... track mouse motion events inside a view.

//______________________________________________________________________________
- (void) updateTrackingAreas
{
   [super updateTrackingAreas];

   if (!fID)
      return;

   const Util::AutoreleasePool pool;

   if (NSArray *trackingArray = [self trackingAreas]) {
      const NSUInteger size = [trackingArray count];
      for (NSUInteger i = 0; i < size; ++i) {
         NSTrackingArea * const t = [trackingArray objectAtIndex : i];
         [self removeTrackingArea : t];
      }
   }

   const NSUInteger trackerOptions = NSTrackingMouseMoved | NSTrackingMouseEnteredAndExited |
                                     NSTrackingActiveInActiveApp | NSTrackingInVisibleRect |
                                     NSTrackingEnabledDuringMouseDrag | NSTrackingCursorUpdate;

   NSRect frame = {};
   frame.size.width = self.fWidth;
   frame.size.height = self.fHeight;

   NSTrackingArea * const tracker = [[NSTrackingArea alloc] initWithRect : frame
                                     options : trackerOptions owner : self userInfo : nil];
   [self addTrackingArea : tracker];
   [tracker release];
}

//______________________________________________________________________________
- (void) updateTrackingAreasAfterRaise
{
   [self updateTrackingAreas];

   for (QuartzView *childView in [self subviews])
      [childView updateTrackingAreasAfterRaise];
}

#pragma mark - X11Drawable protocol.

//______________________________________________________________________________
- (BOOL) fIsPixmap
{
   return NO;
}

//______________________________________________________________________________
- (BOOL) fIsOpenGLWidget
{
   return NO;
}

//______________________________________________________________________________
- (CGFloat) fScaleFactor
{
   return self.fQuartzWindow.fScaleFactor;
}

//______________________________________________________________________________
- (int) fX
{
   return self.frame.origin.x;
}

//______________________________________________________________________________
- (int) fY
{
   return self.frame.origin.y;
}

//______________________________________________________________________________
- (unsigned) fWidth
{
   return self.frame.size.width;
}

//______________________________________________________________________________
- (unsigned) fHeight
{
   return self.frame.size.height;
}

//______________________________________________________________________________
- (void) setDrawableSize : (NSSize) newSize
{
   assert(!(newSize.width < 0) && "-setDrawableSize, width is negative");
   assert(!(newSize.height < 0) && "-setDrawableSize, height is negative");

   //This will cause redraw(?)

   //In X11, resize changes the size, but upper-left corner is not changed.
   //In Cocoa, bottom-left is fixed.
   NSRect frame = self.frame;
   frame.size = newSize;

   self.frame = frame;
}

//______________________________________________________________________________
- (void) setX : (int) x Y : (int) y width : (unsigned) w height : (unsigned) h
{
   NSRect newFrame = {};
   newFrame.origin.x = x;
   newFrame.origin.y = y;
   newFrame.size.width = w;
   newFrame.size.height = h;

   self.frame = newFrame;
}

//______________________________________________________________________________
- (void) setX : (int) x Y : (int) y
{
   NSRect newFrame = self.frame;
   newFrame.origin.x = x;
   newFrame.origin.y = y;

   self.frame = newFrame;
}

//______________________________________________________________________________
- (void) copyImage : (QuartzImage *) srcImage area : (X11::Rectangle) area
          withMask : (QuartzImage *) mask clipOrigin : (X11::Point) clipXY
          toPoint : (X11::Point) dstPoint
{
   //Check parameters.
   assert(srcImage != nil &&
          "-copyImage:area:withMask:clipOrigin:toPoint:, parameter 'srcImage' is nil");
   assert(srcImage.fImage != nil &&
          "-copyImage:area:withMask:clipOrigin:toPoint:, srcImage.fImage is nil");

   //Check self.
   assert(self.fContext != 0 &&
          "-copyImage:area:withMask:clipOrigin:toPoint:, self.fContext is null");

   if (!X11::AdjustCropArea(srcImage, area)) {
      NSLog(@"QuartzView: -copyImage:area:withMask:clipOrigin:toPoint:,"
             " srcRect and copyRect do not intersect");
      return;
   }

   //No RAII for subImage, since it can be really subimage or image itself and
   //in these cases there is no need to release image.
   CGImageRef subImage = 0;
   bool needSubImage = false;
   if (area.fX || area.fY || area.fWidth != srcImage.fWidth || area.fHeight != srcImage.fHeight) {
      needSubImage = true;
      subImage = X11::CreateSubImage(srcImage, area);
      if (!subImage) {
         NSLog(@"QuartzView: -copyImage:area:withMask:clipOrigin:toPoint:,"
                " subimage creation failed");
         return;
      }
   } else
      subImage = srcImage.fImage;

   //Save context state.
   const Quartz::CGStateGuard ctxGuard(self.fContext);

   //Scale and translate to undo isFlipped.
   CGContextTranslateCTM(self.fContext, 0., self.fHeight);
   CGContextScaleCTM(self.fContext, 1., -1.);
   //Set clip mask on a context.

   if (mask) {
      assert(mask.fImage != nil &&
             "-copyImage:area:withMask:clipOrigin:toPoint:, mask.fImage is nil");
      assert(CGImageIsMask(mask.fImage) == true &&
             "-copyImage:area:withMask:clipOrigin:toPoint:, mask.fImage is not a mask");
      //clipXY.fY = X11::LocalYROOTToCocoa(self, clipXY.fY + mask.fHeight);
      const CGFloat clipY = X11::LocalYROOTToCocoa(self, CGFloat(clipXY.fY) + mask.fHeight);
      //const CGRect clipRect = CGRectMake(clipXY.fX, clipXY.fY, mask.fWidth, mask.fHeight);
      const CGRect clipRect = CGRectMake(clipXY.fX, clipY, mask.fWidth, mask.fHeight);
      CGContextClipToMask(self.fContext, clipRect, mask.fImage);
   }

   //Convert from X11 to Cocoa (as soon as we scaled y * -1).
   //dstPoint.fY = X11::LocalYROOTToCocoa(self, dstPoint.fY + area.fHeight);
   const CGFloat dstY = X11::LocalYROOTToCocoa(self, CGFloat(dstPoint.fY) + area.fHeight);
   //const CGRect imageRect = CGRectMake(dstPoint.fX, dstPoint.fY, area.fWidth, area.fHeight);
   const CGRect imageRect = CGRectMake(dstPoint.fX, dstY, area.fWidth, area.fHeight);
   CGContextDrawImage(self.fContext, imageRect, subImage);

   if (needSubImage)
      CGImageRelease(subImage);
}

//______________________________________________________________________________
- (void) copyView : (QuartzView *) srcView area : (X11::Rectangle) area toPoint : (X11::Point) dstPoint
{
   //To copy one "window" to another "window", I have to ask source QuartzView to draw intself into
   //bitmap, and copy this bitmap into the destination view.

   assert(srcView != nil && "-copyView:area:toPoint:, parameter 'srcView' is nil");

   const NSRect frame = [srcView frame];
   //imageRep is in autorelease pool now.
   NSBitmapImageRep * const imageRep = [srcView bitmapImageRepForCachingDisplayInRect : frame];
   if (!imageRep) {
      NSLog(@"QuartzView: -copyView:area:toPoint failed");
      return;
   }

   assert(srcView != nil && "-copyView:area:toPoint:, parameter 'srcView' is nil");
   assert(self.fContext != 0 && "-copyView:area:toPoint, self.fContext is null");

   //It can happen, that src and self are the same.
   //cacheDisplayInRect calls drawRect with bitmap context
   //(and this will reset self.fContext: I have to save/restore it.
   CGContextRef ctx = srcView.fContext;
   srcView.fSnapshotDraw = YES;
   [srcView cacheDisplayInRect : frame toBitmapImageRep : imageRep];
   srcView.fSnapshotDraw = NO;
   srcView.fContext = ctx;

   const CGRect subImageRect = CGRectMake(area.fX, area.fY, area.fWidth, area.fHeight);
   const Util::CFScopeGuard<CGImageRef> subImage(CGImageCreateWithImageInRect(imageRep.CGImage, subImageRect));

   if (!subImage.Get()) {
      NSLog(@"QuartzView: -copyView:area:toPoint, CGImageCreateWithImageInRect failed");
      return;
   }

   const Quartz::CGStateGuard ctxGuard(self.fContext);
   const CGRect imageRect = CGRectMake(dstPoint.fX,
                                       [self visibleRect].size.height - (CGFloat(dstPoint.fY) + area.fHeight),
                                       area.fWidth, area.fHeight);

   CGContextTranslateCTM(self.fContext, 0., [self visibleRect].size.height);
   CGContextScaleCTM(self.fContext, 1., -1.);

   CGContextDrawImage(self.fContext, imageRect, subImage.Get());
}

//______________________________________________________________________________
- (void) copyPixmap : (QuartzPixmap *) srcPixmap area : (X11::Rectangle) area
         withMask : (QuartzImage *) mask clipOrigin : (X11::Point) clipXY toPoint : (X11::Point) dstPoint
{
   //Check parameters.
   assert(srcPixmap != nil && "-copyPixmap:area:withMask:clipOrigin:toPoint:, parameter 'srcPixmap' is nil");

   if (!X11::AdjustCropArea(srcPixmap, area)) {
      NSLog(@"QuartzView: -copyPixmap:area:withMask:clipOrigin:toPoint,"
             " no intersection between pixmap rectangle and cropArea");
      return;
   }

   //Check self.
   assert(self.fContext != 0 &&
          "-copyPixmap:area:withMask:clipOrigin:toPoint:, self.fContext is null");

   //Save context state.
   const Quartz::CGStateGuard ctxGuard(self.fContext);

   CGContextTranslateCTM(self.fContext, 0., self.frame.size.height);//???
   CGContextScaleCTM(self.fContext, 1., -1.);

   const Util::CFScopeGuard<CGImageRef> imageFromPixmap([srcPixmap createImageFromPixmap]);
   assert(imageFromPixmap.Get() != 0 &&
          "-copyPixmap:area:withMask:clipOrigin:toPoint:, createImageFromPixmap failed");

   CGImageRef subImage = 0;
   bool needSubImage = false;
   if (area.fX || area.fY || area.fWidth != srcPixmap.fWidth || area.fHeight != srcPixmap.fHeight) {
      needSubImage = true;
      const CGRect subImageRect = CGRectMake(area.fX, area.fY, area.fHeight, area.fWidth);
      subImage = CGImageCreateWithImageInRect(imageFromPixmap.Get(), subImageRect);
      if (!subImage) {
         NSLog(@"QuartzView: -copyImage:area:withMask:clipOrigin:toPoint:,"
                " subimage creation failed");
         return;
      }
   } else
      subImage = imageFromPixmap.Get();

   if (mask) {
      assert(mask.fImage != nil &&
             "-copyPixmap:area:withMask:clipOrigin:toPoint:, mask.fImage is nil");
      assert(CGImageIsMask(mask.fImage) == true &&
             "-copyPixmap:area:withMask:clipOrigin:toPoint:, mask.fImage is not a mask");

      //clipXY.fY = X11::LocalYROOTToCocoa(self, clipXY.fY + mask.fHeight);
      const CGFloat clipY = X11::LocalYROOTToCocoa(self, CGFloat(clipXY.fY) + mask.fHeight);
      //const CGRect clipRect = CGRectMake(clipXY.fX, clipXY.fY, mask.fWidth, mask.fHeight);
      const CGRect clipRect = CGRectMake(clipXY.fX, clipY, mask.fWidth, mask.fHeight);
      CGContextClipToMask(self.fContext, clipRect, mask.fImage);
   }

   //dstPoint.fY = X11::LocalYCocoaToROOT(self, dstPoint.fY + area.fHeight);
   const CGFloat dstY = X11::LocalYCocoaToROOT(self, CGFloat(dstPoint.fY) + area.fHeight);
   const CGRect imageRect = CGRectMake(dstPoint.fX, dstY, area.fWidth, area.fHeight);
   CGContextDrawImage(self.fContext, imageRect, imageFromPixmap.Get());

   if (needSubImage)
      CGImageRelease(subImage);
}


//______________________________________________________________________________
- (void) copyImage : (QuartzImage *) srcImage area : (X11::Rectangle) area
           toPoint : (X11::Point) dstPoint
{
   assert(srcImage != nil && "-copyImage:area:toPoint:, parameter 'srcImage' is nil");
   assert(srcImage.fImage != nil && "-copyImage:area:toPoint:, srcImage.fImage is nil");
   assert(self.fContext != 0 && "-copyImage:area:toPoint:, fContext is null");

   if (!X11::AdjustCropArea(srcImage, area)) {
      NSLog(@"QuartzView: -copyImage:area:toPoint, image and copy area do not intersect");
      return;
   }

   CGImageRef subImage = 0;
   bool needSubImage = false;
   if (area.fX || area.fY || area.fWidth != srcImage.fWidth || area.fHeight != srcImage.fHeight) {
      needSubImage = true;
      subImage = X11::CreateSubImage(srcImage, area);
      if (!subImage) {
         NSLog(@"QuartzView: -copyImage:area:toPoint:, subimage creation failed");
         return;
      }
   } else
      subImage = srcImage.fImage;

   const Quartz::CGStateGuard ctxGuard(self.fContext);

   CGContextTranslateCTM(self.fContext, 0., self.fHeight);
   CGContextScaleCTM(self.fContext, 1., -1.);

   //dstPoint.fY = X11::LocalYCocoaToROOT(self, dstPoint.fY + area.fHeight);
   const CGFloat dstY = X11::LocalYCocoaToROOT(self, CGFloat(dstPoint.fY) + area.fHeight);
   //const CGRect imageRect = CGRectMake(dstPoint.fX, dstPoint.fY, area.fWidth, area.fHeight);
   const CGRect imageRect = CGRectMake(dstPoint.fX, dstY, area.fWidth, area.fHeight);
   CGContextDrawImage(self.fContext, imageRect, subImage);

   if (needSubImage)
      CGImageRelease(subImage);
}

//______________________________________________________________________________
- (void) copy : (NSObject<X11Drawable> *) src area : (X11::Rectangle) area
     withMask : (QuartzImage *)mask clipOrigin : (X11::Point) clipXY toPoint : (X11::Point) dstPoint
{
   assert(src != nil && "-copy:area:withMask:clipOrigin:toPoint:, parameter 'src' is nil");
   assert(area.fWidth && area.fHeight && "-copy:area:withMask:clipOrigin:toPoint:, area to copy is empty");

   if ([src isKindOfClass : [QuartzWindow class]]) {
      //Forget about mask (can I have it???)
      QuartzWindow * const qw = (QuartzWindow *)src;
      //Will not work with OpenGL.
      [self copyView : (QuartzView *)qw.fContentView area : area toPoint : dstPoint];
   } else if ([src isKindOfClass : [QuartzView class]]) {
      //Forget about mask (can I have it???)
      [self copyView : (QuartzView *)src area : area toPoint : dstPoint];
   } else if ([src isKindOfClass : [QuartzPixmap class]]) {
      [self copyPixmap : (QuartzPixmap *)src area : area withMask : mask clipOrigin : clipXY toPoint : dstPoint];
   } else if ([src isKindOfClass : [QuartzImage class]]) {
      [self copyImage : (QuartzImage *)src area : area withMask : mask clipOrigin : clipXY toPoint : dstPoint];
   } else {
      assert(0 && "-copy:area:withMask:clipOrigin:toPoint:, src is of unknown type");
   }
}

//______________________________________________________________________________
- (unsigned char *) readColorBits : (X11::Rectangle) area
{
   //This is quite a bad idea - to read pixels back from a view,
   //but our GUI does exactly this. In case of Cocoa it's expensive
   //and not guaranteed to work.

   assert(area.fWidth && area.fHeight && "-readColorBits:, area to copy is empty");

   //int, not unsigned or something - to keep it simple.
   const NSRect visRect = [self visibleRect];
   const X11::Rectangle srcRect(int(visRect.origin.x), int(visRect.origin.y),
                                unsigned(visRect.size.width), unsigned(visRect.size.height));

   if (!X11::AdjustCropArea(srcRect, area)) {
      NSLog(@"QuartzView: -readColorBits:, visible rect of view and copy area do not intersect");
      return nullptr;
   }

   //imageRep is autoreleased.
   NSBitmapImageRep * const imageRep = [self bitmapImageRepForCachingDisplayInRect : visRect];
   if (!imageRep) {
      NSLog(@"QuartzView: -readColorBits:, bitmapImageRepForCachingDisplayInRect failed");
      return nullptr;
   }

   CGContextRef ctx = self.fContext; //Save old context if any.
   [self cacheDisplayInRect : visRect toBitmapImageRep : imageRep];
   self.fContext = ctx; //Restore old context.
   //
   const NSInteger bitsPerPixel = [imageRep bitsPerPixel];

   assert(bitsPerPixel == 32 && "-readColorBits:, no alpha channel???");
   const NSInteger bytesPerRow = [imageRep bytesPerRow];
   unsigned dataWidth = bytesPerRow / (bitsPerPixel / 8);//assume an octet :(

   unsigned char *srcData = nullptr;
   std::vector<unsigned char> downscaled;
   if ([self.window.screen backingScaleFactor] > 1 && imageRep.CGImage) {
      downscaled = X11::DownscaledImageData(area.fWidth, area.fHeight, imageRep.CGImage);
      if (downscaled.size())
         srcData = &downscaled[0];
      dataWidth = area.fWidth;
   } else
      srcData = [imageRep bitmapData];

   if (!srcData) {
      NSLog(@"QuartzView: -readColorBits:, failed to obtain backing store contents");
      return nullptr;
   }

   //We have a source data now. Let's allocate buffer for ROOT's GUI and convert source data.
   unsigned char *data = nullptr;

   try {
      data = new unsigned char[area.fWidth * area.fHeight * 4];//bgra?
   } catch (const std::bad_alloc &) {
      NSLog(@"QuartzView: -readColorBits:, memory allocation failed");
      return nullptr;
   }

   unsigned char *dstPixel = data;
   const unsigned char *line = srcData + area.fY * dataWidth * 4;
   const unsigned char *srcPixel = line + area.fX * 4;

   for (unsigned i = 0; i < area.fHeight; ++i) {
      for (unsigned j = 0; j < area.fWidth; ++j, srcPixel += 4, dstPixel += 4) {
         dstPixel[0] = srcPixel[2];
         dstPixel[1] = srcPixel[1];
         dstPixel[2] = srcPixel[0];
         dstPixel[3] = srcPixel[3];
      }

      line += dataWidth * 4;
      srcPixel = line + area.fX * 4;
   }

   return data;
}

//______________________________________________________________________________
- (void) setFBackgroundPixmap : (QuartzImage *) pixmap
{
   if (fBackgroundPixmap != pixmap) {
      [fBackgroundPixmap release];
      if (pixmap)
         fBackgroundPixmap = [pixmap retain];
      else
         fBackgroundPixmap = nil;
   }
}

//______________________________________________________________________________
- (QuartzImage *) fBackgroundPixmap
{
   //I do not autorelease, screw this idiom!

   return fBackgroundPixmap;
}

//______________________________________________________________________________
- (int) fMapState
{
   if ([self isHidden])
      return kIsUnmapped;

   for (QuartzView *parent = fParentView; parent; parent = parent.fParentView) {
      if ([parent isHidden])
         return kIsUnviewable;
   }

   return kIsViewable;
}

//______________________________________________________________________________
- (BOOL) fHasFocus
{
   //With the latest update clang became a bit more stupid.
   //Let's write a stupid useless cargo cult code
   //to make IT SHUT THE F... UP.
   (void)fHasFocus;
   return NO;
}

//______________________________________________________________________________
- (void) setFHasFocus : (BOOL) focus
{
#pragma unused(focus)
   //With the latest update clang became a bit more stupid.
   //Let's write a stupid useless cargo cult code
   //to make IT SHUT THE F... UP.
   (void)fHasFocus;
}

//______________________________________________________________________________
- (QuartzPixmap *) fBackBuffer
{
   return fBackBuffer;//No autorelease, I know the object's lifetime myself.
}

//______________________________________________________________________________
- (void) setFBackBuffer : (QuartzPixmap *) backBuffer
{
   if (fBackBuffer != backBuffer) {
      [fBackBuffer release];

      if (backBuffer)
         fBackBuffer = [backBuffer retain];
      else
         fBackBuffer = nil;
   }
}

//______________________________________________________________________________
- (NSView<X11Window> *) fContentView
{
   return self;
}

//______________________________________________________________________________
- (QuartzWindow *) fQuartzWindow
{
   return (QuartzWindow *)[self window];
}

//______________________________________________________________________________
- (void) activatePassiveGrab
{
   fCurrentGrabType = X11::kPGPassiveGrab;
}

//______________________________________________________________________________
- (void) activateImplicitGrab
{
   fCurrentGrabType = X11::kPGImplicitGrab;
}

//______________________________________________________________________________
- (void) activateGrab : (unsigned) eventMask ownerEvents : (BOOL) ownerEvents
{
   fCurrentGrabType = X11::kPGActiveGrab;
   fActiveGrabEventMask = eventMask;
   fActiveGrabOwnerEvents = ownerEvents;
}

//______________________________________________________________________________
- (void) cancelGrab
{
   fCurrentGrabType = X11::kPGNoGrab;
   fActiveGrabEventMask = 0;
   fActiveGrabOwnerEvents = YES;
}

//______________________________________________________________________________
- (BOOL) acceptsCrossingEvents : (unsigned) eventMask
{
   bool accepts = fEventMask & eventMask;

   //In ROOT passive grabs are always with owner_events == true.
   if (fCurrentGrabType == X11::kPGPassiveGrab)
      accepts = accepts || (fPassiveGrabEventMask & eventMask);

   if (fCurrentGrabType == X11::kPGActiveGrab) {
      if (fActiveGrabOwnerEvents)
         accepts = accepts || (fActiveGrabOwnerEvents & eventMask);
      else
         accepts = fActiveGrabOwnerEvents & eventMask;
   }

   return accepts;
}

//______________________________________________________________________________
- (void) addChild : (NSView<X11Window> *) child
{
   assert(child != nil && "-addChild:, parameter 'child' is nil");

   [self addSubview : child];
   child.fParentView = self;
}

//______________________________________________________________________________
- (void) getAttributes : (WindowAttributes_t *) attr
{
   assert(attr != 0 && "-getAttributes:, parameter 'attr' is null");

   X11::GetWindowAttributes(self, attr);
}

//______________________________________________________________________________
- (void) setAttributes : (const SetWindowAttributes_t *)attr
{
   assert(attr != 0 && "-setAttributes:, parameter 'attr' is null");

#ifdef DEBUG_ROOT_COCOA
   log_attributes(attr, fID);
#endif

   X11::SetWindowAttributes(attr, self);
}

//______________________________________________________________________________
- (void) mapRaised
{
   //Move view to the top of subviews.
   QuartzView * const parent = fParentView;
   [self removeFromSuperview];
   [parent addSubview : self];
   [self setHidden : NO];
}

//______________________________________________________________________________
- (void) mapWindow
{
   [self setHidden : NO];
}

//______________________________________________________________________________
- (void) mapSubwindows
{
   for (QuartzView * v in [self subviews])
      [v setHidden : NO];
}

//______________________________________________________________________________
- (void) unmapWindow
{
   [self setHidden : YES];
}

//______________________________________________________________________________
- (BOOL) fIsOverlapped
{
   return fIsOverlapped;
}

//______________________________________________________________________________
- (void) setOverlapped : (BOOL) overlap
{
   fIsOverlapped = overlap;
   for (NSView<X11Window> *child in [self subviews])
      [child setOverlapped : overlap];
}

//______________________________________________________________________________
- (void) raiseWindow
{
   //Now, I can not remove window and add it ...
   //For example, if you click on a tab, this:
   //1. Creates (potentially) a passive button grab
   //2. Raises this tab - changes the window order.
   //3. On a button release - grab is release.
   //The tough problem is, if I remove a view from subviews
   //and add it ... it will never receve the
   //release event thus a grab will 'hang' on
   //view leading to bugs and artifacts.
   //So instead I have to ... SORT!!!!!

   using namespace X11;//Comparators.

   for (QuartzView *sibling in [fParentView subviews]) {
      if (self == sibling)
         continue;
      if ([sibling isHidden])
         continue;

      if (NSEqualRects(sibling.frame, self.frame)) {
         [sibling setOverlapped : YES];
         [sibling setHidden : YES];
      }
   }

   [self setOverlapped : NO];
   //
   [self setHidden : NO];
   //
   [fParentView sortSubviewsUsingFunction : CompareViewsToRaise context : (void *)self];
   //
   [self updateTrackingAreasAfterRaise];
   //
   [self setNeedsDisplay : YES];
}

//______________________________________________________________________________
- (void) lowerWindow
{
   //See comment about sorting in -raiseWindow.

   using namespace X11;

   NSEnumerator * const reverseEnumerator = [[fParentView subviews] reverseObjectEnumerator];
   for (QuartzView *sibling in reverseEnumerator) {
      if (sibling == self)
         continue;

      if (NSEqualRects(sibling.frame, self.frame)) {
         [sibling setOverlapped : NO];
         //
         [sibling setHidden : NO];
         //
         [sibling setNeedsDisplay : YES];
         [self setOverlapped : YES];
         //
         [self setHidden : YES];
         //
         break;
      }
   }

   [fParentView sortSubviewsUsingFunction : CompareViewsToLower context : (void*)self];
}

//______________________________________________________________________________
- (BOOL) isFlipped
{
   //Now view's placement, geometry, moving and resizing can be
   //done with ROOT's (X11) coordinates without conversion - we're are 'flipped'.
   return YES;
}

//______________________________________________________________________________
- (void) configureNotifyTree
{
   if (self.fMapState == kIsViewable || fIsOverlapped == YES) {
      if (fEventMask & kStructureNotifyMask) {
         assert(dynamic_cast<TGCocoa *>(gVirtualX) &&
                "-configureNotifyTree, gVirtualX is either null or has type different from TGCocoa");
         TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
         vx->GetEventTranslator()->GenerateConfigureNotifyEvent(self, self.frame);
      }

      for (NSView<X11Window> *v in [self subviews])
         [v configureNotifyTree];
   }
}

#pragma mark - Key grabs.

//______________________________________________________________________________
- (void) addPassiveKeyGrab : (unichar) keyCode modifiers : (NSUInteger) modifiers
{
   [self removePassiveKeyGrab : keyCode modifiers : modifiers];
   PassiveKeyGrab * const newGrab = [[PassiveKeyGrab alloc] initWithKey : keyCode
                                                              modifiers : modifiers];
   [fPassiveKeyGrabs addObject : newGrab];
   [newGrab release];
}

//______________________________________________________________________________
- (void) removePassiveKeyGrab : (unichar) keyCode modifiers : (NSUInteger) modifiers
{
   const NSUInteger count = [fPassiveKeyGrabs count];
   for (NSUInteger i = 0; i < count; ++i) {
      PassiveKeyGrab *grab = [fPassiveKeyGrabs objectAtIndex : i];
      if ([grab matchKey : keyCode modifiers : modifiers]) {
         [fPassiveKeyGrabs removeObjectAtIndex : i];
         break;
      }
   }
}

//______________________________________________________________________________
- (PassiveKeyGrab *) findPassiveKeyGrab : (unichar) keyCode modifiers : (NSUInteger) modifiers
{
   NSEnumerator * const enumerator = [fPassiveKeyGrabs objectEnumerator];
   while (PassiveKeyGrab *grab = (PassiveKeyGrab *)[enumerator nextObject]) {
      if ([grab matchKey : keyCode modifiers : modifiers])
         return grab;
   }

   return nil;
}

//______________________________________________________________________________
- (PassiveKeyGrab *) findPassiveKeyGrab : (unichar) keyCode
{
   //Do not check modifiers.
   NSEnumerator * const enumerator = [fPassiveKeyGrabs objectEnumerator];
   while (PassiveKeyGrab *grab = (PassiveKeyGrab *)[enumerator nextObject]) {
      if ([grab matchKey : keyCode])
         return grab;
   }

   return nil;
}

#pragma mark - Painting mechanics.

//______________________________________________________________________________
- (void) drawRect : (NSRect) dirtyRect
{
#pragma unused(dirtyRect)

   using namespace X11;

   if (fID) {
      if (TGWindow * const window = gClient->GetWindowById(fID)) {
         //It's never painted, parent renders child. true == check the parent also.
         if (ViewIsTextViewFrame(self, true) ||ViewIsHtmlViewFrame(self, true))
            return;

         NSGraphicsContext * const nsContext = [NSGraphicsContext currentContext];
         assert(nsContext != nil && "-drawRect:, currentContext returned nil");

         TGCocoa * const vx = (TGCocoa *)gVirtualX;
         vx->CocoaDrawON();

         fContext = (CGContextRef)[nsContext graphicsPort];
         assert(fContext != 0 && "-drawRect:, graphicsPort returned null");

         const Quartz::CGStateGuard ctxGuard(fContext);

         //Non-rectangular windows.
         if (self.fQuartzWindow.fShapeCombineMask)
            X11::ClipToShapeMask(self, fContext);

         // This code used to use TObject::InheritsFrom, however since this is
         // run under the AppKit, we can not call core/meta functions, otherwise
         // we will run into deadlocks.
         if (dynamic_cast<const TGContainer*>(window))//It always has an ExposureMask.
            vx->GetEventTranslator()->GenerateExposeEvent(self, [self visibleRect]);

         if (fEventMask & kExposureMask) {
            if (ViewIsTextView(self)) {
               //Send Expose event, using child view (this is how it's done in GUI :( ).
               [NSColor.whiteColor setFill];
               NSRectFill(dirtyRect);
               NSView<X11Window> * const viewFrame = FrameForTextView(self);
               if (viewFrame)//Now we set fExposedRegion for TGView.
                  vx->GetEventTranslator()->GenerateExposeEvent(viewFrame, viewFrame.visibleRect);
            }

            if (ViewIsHtmlView(self)) {
               NSView<X11Window> *const viewFrame = FrameForHtmlView(self);
               if (viewFrame)
                  vx->GetEventTranslator()->GenerateExposeEvent(viewFrame, viewFrame.visibleRect);
            }

            //Ask ROOT's widget/window to draw itself.
            gClient->NeedRedraw(window, kTRUE);

            if (!fSnapshotDraw) {
               //If Cocoa repaints widget, cancel all ROOT's "outside of paint event"
               //rendering into this widget.
               gClient->CancelRedraw(window);
               vx->GetCommandBuffer()->RemoveGraphicsOperationsForWindow(fID);
            }
         }

         if (fBackBuffer) {
            //Very "special" window.
            const X11::Rectangle copyArea(0, 0, fBackBuffer.fWidth, fBackBuffer.fHeight);
            [self copy : fBackBuffer area : copyArea withMask : nil
                  clipOrigin : X11::Point() toPoint : X11::Point()];
         }

         vx->CocoaDrawOFF();
#ifdef DEBUG_ROOT_COCOA
         CGContextSetRGBStrokeColor(fContext, 1., 0., 0., 1.);
         CGContextStrokeRect(fContext, dirtyRect);
#endif

         fContext = 0;
      } else {
#ifdef DEBUG_ROOT_COCOA
         NSLog(@"QuartzView: -drawRect: method, no window for id %u was found", fID);
#endif
      }
   }
}

#pragma mark - Geometry.

//______________________________________________________________________________
- (void) setFrame : (NSRect) newFrame
{
   //In case of TBrowser, setFrame started infinite recursion:
   //HandleConfigure for embedded main frame emits signal, slot
   //calls layout, layout calls setFrame -> HandleConfigure and etc. etc.
   if (NSEqualRects(newFrame, self.frame))
      return;

   [super setFrame : newFrame];
}

//______________________________________________________________________________
- (void) setFrameSize : (NSSize) newSize
{
   //Check, if setFrameSize calls setFrame.

   [super setFrameSize : newSize];

   if ((fEventMask & kStructureNotifyMask) && (self.fMapState == kIsViewable || fIsOverlapped == YES)) {
      assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
             "setFrameSize:, gVirtualX is either null or has a type, different from TGCocoa");
      TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
      vx->GetEventTranslator()->GenerateConfigureNotifyEvent(self, self.frame);
   }

   [self setNeedsDisplay : YES];//?
}

#pragma mark - Event handling.

//______________________________________________________________________________
- (void) mouseDown : (NSEvent *) theEvent
{
   assert(fID != 0 && "-mouseDown:, fID is 0");

   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "-mouseDown:, gVirtualX is either null or has a type, different from TGCocoa");
   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->GetEventTranslator()->GenerateButtonPressEvent(self, theEvent, kButton1);
}

//______________________________________________________________________________
- (void) scrollWheel : (NSEvent*) theEvent
{
   assert(fID != 0 && "-scrollWheel:, fID is 0");


   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "-scrollWheel:, gVirtualX is either null or has a type, different from TGCocoa");

   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   const CGFloat deltaY = [theEvent deltaY];
   if (deltaY < 0) {
      vx->GetEventTranslator()->GenerateButtonPressEvent(self, theEvent, kButton5);
      vx->GetEventTranslator()->GenerateButtonReleaseEvent(self, theEvent, kButton5);
   } else if (deltaY > 0) {
      vx->GetEventTranslator()->GenerateButtonPressEvent(self, theEvent, kButton4);
      vx->GetEventTranslator()->GenerateButtonReleaseEvent(self, theEvent, kButton4);
   }
}

#ifdef DEBUG_ROOT_COCOA
//______________________________________________________________________________
- (void) printViewInformation
{
   assert(fID != 0 && "-printWindowInformation, fID is 0");
   const TGWindow * const window = gClient->GetWindowById(fID);
   assert(window != 0 && "printWindowInformation, window not found");

   NSLog(@"-----------------View %u info:---------------------", fID);
   NSLog(@"ROOT's window class is %s", window->IsA()->GetName());
   NSLog(@"event mask is:");
   print_mask_info(fEventMask);
   NSLog(@"grab mask is:");
   print_mask_info(fPassiveGrabEventMask);
   NSLog(@"view's geometry: x == %g, y == %g, w == %g, h == %g", self.frame.origin.x,
         self.frame.origin.y, self.frame.size.width, self.frame.size.height);
   NSLog(@"----------------End of view info------------------");
}
#endif

//______________________________________________________________________________
- (void) rightMouseDown : (NSEvent *) theEvent
{
   assert(fID != 0 && "-rightMouseDown:, fID is 0");

#ifdef DEBUG_ROOT_COCOA
   [self printViewInformation];
#endif

   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "-rightMouseDown:, gVirtualX is either null or has type different from TGCocoa");
   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->GetEventTranslator()->GenerateButtonPressEvent(self, theEvent, kButton3);
}

//______________________________________________________________________________
- (void) otherMouseDown : (NSEvent *) theEvent
{
   assert(fID != 0 && "-otherMouseDown:, fID is 0");

   //Funny enough, [theEvent buttonNumber] is not the same thing as button masked in [NSEvent pressedMouseButtons],
   //button number actually is a kind of right operand for bitshift for pressedMouseButtons.
   if ([theEvent buttonNumber] == 2) {//this '2' will correspond to '4' in pressedMouseButtons.
      //I do not care about mouse buttons after left/right/wheel - ROOT does not have
      //any code for this.
      assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
             "-otherMouseDown:, gVirtualX is either null or has type different from TGCocoa");
      TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
      vx->GetEventTranslator()->GenerateButtonPressEvent(self, theEvent, kButton2);
   }
}

//______________________________________________________________________________
- (void) mouseUp : (NSEvent *) theEvent
{
   assert(fID != 0 && "-mouseUp:, fID is 0");

   assert(dynamic_cast<TGCocoa *>(gVirtualX) &&
          "-mouseUp:, gVirtualX is either null or has type different from TGCocoa");
   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->GetEventTranslator()->GenerateButtonReleaseEvent(self, theEvent, kButton1);
}

//______________________________________________________________________________
- (void) rightMouseUp : (NSEvent *) theEvent
{

   assert(fID != 0 && "-rightMouseUp:, fID is 0");

   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "-rightMouseUp:, gVirtualX is either null or has type different from TGCocoa");

   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->GetEventTranslator()->GenerateButtonReleaseEvent(self, theEvent, kButton3);
}

//______________________________________________________________________________
- (void) otherMouseUp : (NSEvent *) theEvent
{
   assert(fID != 0 && "-otherMouseUp:, fID is 0");

   //Here I assume it's always kButton2.
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "-otherMouseUp:, gVirtualX is either null or has type different from TGCocoa");
   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->GetEventTranslator()->GenerateButtonReleaseEvent(self, theEvent, kButton2);
}

//______________________________________________________________________________
- (void) mouseEntered : (NSEvent *) theEvent
{
   assert(fID != 0 && "-mouseEntered:, fID is 0");
   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "-mouseEntered:, gVirtualX is null or not of TGCocoa type");

   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->GetEventTranslator()->GenerateCrossingEvent(theEvent);
}

//______________________________________________________________________________
- (void) mouseExited : (NSEvent *) theEvent
{
   assert(fID != 0 && "-mouseExited:, fID is 0");

   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "-mouseExited:, gVirtualX is null or not of TGCocoa type");

   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->GetEventTranslator()->GenerateCrossingEvent(theEvent);
}

//______________________________________________________________________________
- (void) mouseMoved : (NSEvent *) theEvent
{
   assert(fID != 0 && "-mouseMoved:, fID is 0");

   if (fParentView)//Suppress events in all views, except the top-level one.
      return;

   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "-mouseMoved:, gVirtualX is null or not of TGCocoa type");

   TGCocoa *vx = static_cast<TGCocoa *>(gVirtualX);
   vx->GetEventTranslator()->GeneratePointerMotionEvent(theEvent);
}

//______________________________________________________________________________
- (void) mouseDragged : (NSEvent *) theEvent
{
   assert(fID != 0 && "-mouseDragged:, fID is 0");

   TGCocoa * const vx = dynamic_cast<TGCocoa *>(gVirtualX);
   assert(vx != 0 && "-mouseDragged:, gVirtualX is null or not of TGCocoa type");

   vx->GetEventTranslator()->GeneratePointerMotionEvent(theEvent);
}

//______________________________________________________________________________
- (void) rightMouseDragged : (NSEvent *) theEvent
{
   assert(fID != 0 && "-rightMouseDragged:, fID is 0");

   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "-rightMouseDragged:, gVirtualX is null or not of TGCocoa type");

   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->GetEventTranslator()->GeneratePointerMotionEvent(theEvent);
}

//______________________________________________________________________________
- (void) otherMouseDragged : (NSEvent *) theEvent
{
   assert(fID != 0 && "-otherMouseDragged:, fID is 0");

   if ([theEvent buttonNumber] == 2) {
      assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
             "-otherMouseDragged:, gVirtualX is null or not of TGCocoa type");
      TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
      vx->GetEventTranslator()->GeneratePointerMotionEvent(theEvent);
   }
}

//______________________________________________________________________________
- (void) keyDown : (NSEvent *) theEvent
{
   assert(fID != 0 && "-keyDown:, fID is 0");

   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "-keyDown:, gVirtualX is null or not of TGCocoa type");

   NSView<X11Window> *eventView = self;
   if (NSView<X11Window> *pointerView = X11::FindViewUnderPointer())
      eventView = pointerView;

   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   vx->GetEventTranslator()->GenerateKeyPressEvent(eventView, theEvent);
}

//______________________________________________________________________________
- (void) keyUp : (NSEvent *) theEvent
{
   assert(fID != 0 && "-keyUp:, fID is 0");

   assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 &&
          "-keyUp:, gVirtualX is null or not of TGCocoa type");

   TGCocoa * const vx = static_cast<TGCocoa *>(gVirtualX);
   NSView<X11Window> *eventView = self;
   if (NSView<X11Window> *pointerView = X11::FindViewUnderPointer())
      eventView = pointerView;

   vx->GetEventTranslator()->GenerateKeyReleaseEvent(eventView, theEvent);
}

#pragma mark - First responder stuff.

//______________________________________________________________________________
- (BOOL) acceptsFirstMouse : (NSEvent *) theEvent
{
#pragma unused(theEvent)
   return YES;
}

//______________________________________________________________________________
- (BOOL) acceptsFirstResponder
{
   return YES;
}

#pragma mark - Cursors.

//______________________________________________________________________________
- (void) setFCurrentCursor : (ECursor) cursor
{
   if (cursor != fCurrentCursor) {
      fCurrentCursor = cursor;
      [self.fQuartzWindow invalidateCursorRectsForView : self];
   }
}

//______________________________________________________________________________
- (NSCursor *) createCustomCursor
{
   const char *pngFileName = 0;

   switch (fCurrentCursor) {
   case kMove:
      pngFileName = "move_cursor.png";
      break;
   case kArrowHor:
      pngFileName = "hor_arrow_cursor.png";
      break;
   case kArrowVer:
      pngFileName = "ver_arrow_cursor.png";
      break;
   case kArrowRight:
      pngFileName = "right_arrow_cursor.png";
      break;
   case kRotate:
      pngFileName = "rotate.png";
      break;
   case kBottomLeft:
   case kTopRight:
      pngFileName = "top_right_cursor.png";
      break;
   case kTopLeft:
   case kBottomRight:
      pngFileName = "top_left_cursor.png";
      break;
   default:;
   }

   if (pngFileName) {
      const char * const path = gSystem->Which(TROOT::GetIconPath(), pngFileName, kReadPermission);
      const Util::ScopedArray<const char> arrayGuard(path);

      if (!path || path[0] == 0) {
         //File was not found.
         return nil;
      }

      NSString *nsPath = [NSString stringWithFormat : @"%s", path];//in autorelease pool.
      NSImage * const cursorImage = [[NSImage alloc] initWithContentsOfFile : nsPath];

      if (!cursorImage)
         return nil;

      NSPoint hotSpot = X11::GetCursorHotStop(cursorImage, fCurrentCursor);
      NSCursor * const customCursor = [[[NSCursor alloc] initWithImage : cursorImage
                                                         hotSpot : hotSpot] autorelease];

      [cursorImage release];

      return customCursor;
   }

   return nil;
}

//______________________________________________________________________________
- (void) resetCursorRects
{
   if (NSCursor * const cursor = X11::CreateCursor(fCurrentCursor))
      [self addCursorRect : self.visibleRect cursor : cursor];
}

//______________________________________________________________________________
- (void) cursorUpdate
{
   if (NSCursor * const cursor = X11::CreateCursor(fCurrentCursor)) {
      // NB: [window invalidateCursorRectsForView] called here has the
      // same problem as commented below in -cursorUpdate:.
      [cursor set];
   }
}

//______________________________________________________________________________
- (void) cursorUpdate : (NSEvent *) event
{
#pragma unused(event)
   // It looks like [NSCursor set] method does not work properly when called from
   // cursorUpdate:, having, say, a parent frame with 'arrow' cursor and a child (completely
   // filling its parent's area) with 'cross', it happens the 'cross' cursor is not always
   // set correctly, for example:
   // if we have a TCanvas and resize it, cursor is 'arrow' inside this canvas,
   // though it must be 'cross'. This all, as it always happesn with "thinking different"
   // Apple is somehow related to run loop or something. As always, it's not documented,
   // so Apple can continue to think different. The idea with performSelector comes from:
   // http://stackoverflow.com/questions/8430236/nscursor-set-method-has-no-effect
   // Or may be it's just a bug:
   // http://stackoverflow.com/questions/13901232/nscursor-set-not-working-on-unfocused-window
   [self performSelector : @selector(cursorUpdate) withObject : nil afterDelay : 0.05f];
}

#pragma mark - Emulated X11 properties.

//______________________________________________________________________________
- (void) setProperty : (const char *) propName data : (unsigned char *) propData
         size : (unsigned) dataSize forType : (Atom_t) dataType format : (unsigned) format
{
   assert(propName != 0 && "-setProperty:data:size:forType:, parameter 'propName' is null");
   assert(propData != 0 && "-setProperty:data:size:forType:, parameter 'propData' is null");
   assert(dataSize != 0 && "-setProperty:data:size:forType:, parameter 'dataSize' is 0");

   NSString * const key = [NSString stringWithCString : propName encoding : NSASCIIStringEncoding];
   QuartzWindowProperty * property = (QuartzWindowProperty *)[fX11Properties valueForKey : key];

   //At the moment (and I think this will never change) TGX11 always calls XChangeProperty with PropModeReplace.
   if (property)
      [property resetPropertyData : propData size : dataSize type : dataType format : format];
   else {
      //No property found, add a new one.
      property = [[QuartzWindowProperty alloc] initWithData : propData size : dataSize
                                               type : dataType format : format];
      [fX11Properties setObject : property forKey : key];
      [property release];
   }
}

//______________________________________________________________________________
- (BOOL) hasProperty : (const char *) propName
{
   assert(propName != 0 && "-hasProperty:, propName parameter is null");

   NSString * const key = [NSString stringWithCString : propName encoding : NSASCIIStringEncoding];
   QuartzWindowProperty * const property = (QuartzWindowProperty *)[fX11Properties valueForKey : key];

   return property != nil;
}

//______________________________________________________________________________
- (unsigned char *) getProperty : (const char *) propName returnType : (Atom_t *) type
   returnFormat : (unsigned *) format nElements : (unsigned *) nElements
{
   assert(propName != 0 &&
          "-getProperty:returnType:returnFormat:nElements:, parameter 'propName' is null");
   assert(type != 0 &&
          "-getProperty:returnType:returnFormat:nElements:, parameter 'type' is null");
   assert(format != 0 &&
          "-getProperty:returnType:returnFormat:nElements:, parameter 'format' is null");
   assert(nElements != 0 &&
          "-getProperty:returnType:returnFormat:nElements:, parameter 'nElements' is null");

   NSString * const key = [NSString stringWithCString : propName encoding : NSASCIIStringEncoding];
   QuartzWindowProperty * const property = (QuartzWindowProperty *)[fX11Properties valueForKey : key];
   assert(property != 0 &&
          "-getProperty:returnType:returnFormat:nElements, property not found");

   NSData * const propData = property.fPropertyData;

   const NSUInteger dataSize = [propData length];
   unsigned char *buff = 0;
   try {
      buff = new unsigned char[dataSize]();
   } catch (const std::bad_alloc &) {
      //Hmm, can I log, if new failed? :)
      NSLog(@"QuartzWindow: -getProperty:returnType:returnFormat:nElements:,"
            " memory allocation failed");
      return 0;
   }

   [propData getBytes : buff length : dataSize];
   *format = property.fFormat;

   *nElements = dataSize;

   if (*format == 16)
      *nElements= dataSize / 2;
   else if (*format == 32)
      *nElements = dataSize / 4;

   *type = property.fType;

   return buff;
}

//______________________________________________________________________________
- (void) removeProperty : (const char *) propName
{
   assert(propName != 0 && "-removeProperty:, parameter 'propName' is null");

   NSString * const key = [NSString stringWithCString : propName
                           encoding : NSASCIIStringEncoding];
   [fX11Properties removeObjectForKey : key];
}

//DND
//______________________________________________________________________________
- (NSDragOperation) draggingEntered : (id<NSDraggingInfo>) sender
{
   NSPasteboard * const pasteBoard = [sender draggingPasteboard];
   const NSDragOperation sourceDragMask = [sender draggingSourceOperationMask];

   if ([[pasteBoard types] containsObject : NSFilenamesPboardType] && (sourceDragMask & NSDragOperationCopy))
      return NSDragOperationCopy;

   return NSDragOperationNone;
}

//______________________________________________________________________________
- (BOOL) performDragOperation : (id<NSDraggingInfo>) sender
{
   //We can drag some files (images, pdfs, source code files) from
   //finder to ROOT's window (mainly TCanvas or text editor).
   //The logic is totally screwed here :((( - ROOT will try to
   //read a property of some window (not 'self', unfortunately) -
   //this works since on Window all data is in a global clipboard
   //(on X11 it simply does not work at all).
   //I'm attaching the file name as a property for the top level window,
   //there is no other way to make this data accessible for ROOT.

   NSPasteboard * const pasteBoard = [sender draggingPasteboard];
   const NSDragOperation sourceDragMask = [sender draggingSourceOperationMask];

   if ([[pasteBoard types] containsObject : NSFilenamesPboardType] && (sourceDragMask & NSDragOperationCopy)) {

      //Here I try to put string ("file://....") into window's property to make
      //it accesible from ROOT's GUI.
      const Atom_t textUriAtom = gVirtualX->InternAtom("text/uri-list", kFALSE);

      NSArray * const files = [pasteBoard propertyListForType : NSFilenamesPboardType];
      for (NSString *path in files) {
         //ROOT can not process several files, use the first one.
         NSString * const item = [@"file://" stringByAppendingString : path];
         //Yes, only ASCII encoding, but after all, ROOT's not able to work with NON-ASCII strings.
         const NSUInteger len = [item lengthOfBytesUsingEncoding : NSASCIIStringEncoding] + 1;
         try {
            std::vector<unsigned char> propertyData(len);
            [item getCString : (char *)&propertyData[0] maxLength : propertyData.size()
             encoding : NSASCIIStringEncoding];
            //There is no any guarantee, that this will ever work, logic in TGDNDManager is totally crazy.
            NSView<X11Window> * const targetView = self.fQuartzWindow.fContentView;
            [targetView setProperty : "_XC_DND_DATA" data : &propertyData[0]
                        size : propertyData.size() forType : textUriAtom format : 8];
         } catch (const std::bad_alloc &) {
            //Hehe, can I log something in case of bad_alloc??? ;)
            NSLog(@"QuartzView: -performDragOperation:, memory allocation failed");
            return NO;
         }

         break;
      }

      //Property is attached now.

      //Gdk on windows creates three events on file drop (WM_DROPFILES): XdndEnter, XdndPosition, XdndDrop.
      //1. Dnd enter.
      Event_t event1 = {};
      event1.fType = kClientMessage;
      event1.fWindow = fID;
      event1.fHandle = gVirtualX->InternAtom("XdndEnter", kFALSE);
      event1.fUser[0] = long(fID);
      event1.fUser[2] = textUriAtom;//gVirtualX->InternAtom("text/uri-list", kFALSE);
      //
      gVirtualX->SendEvent(fID, &event1);

      //2. Dnd position.
      Event_t event2 = {};
      event2.fType = kClientMessage;
      event2.fWindow = fID;
      event2.fHandle = gVirtualX->InternAtom("XdndPosition", kFALSE);
      event2.fUser[0] = long(fID);
      event2.fUser[2] = 0;//Here I have to pack x and y for drop coordinates, shifting by 16 bits.
      NSPoint dropPoint = [sender draggingLocation];
      //convertPointFromBase is deprecated.
      //dropPoint = [self convertPointFromBase : dropPoint];
      dropPoint = [self convertPoint : dropPoint fromView : nil];
      //
      dropPoint = X11::TranslateToScreen(self, dropPoint);
      event2.fUser[2] = UShort_t(dropPoint.y) | (UShort_t(dropPoint.x) << 16);

      gVirtualX->SendEvent(fID, &event2);

      Event_t event3 = {};
      event3.fType = kClientMessage;
      event3.fWindow = fID;
      event3.fHandle = gVirtualX->InternAtom("XdndDrop", kFALSE);

      gVirtualX->SendEvent(fID, &event3);
   }

   return YES;//Always ok, even if file type is not supported - no need in "animation".
}

@end
