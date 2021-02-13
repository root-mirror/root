// @(#)root/graf2d:$Id$
// Author: Timur Pocheptsov   22/11/2011

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//#define NDEBUG

#include "TGCocoa.h"

// We want to pickup ROOT's glew and not the system OpenGL coming from:
// ROOTOpenGLView.h ->QuartzWindow.h->Cocoa.h
// Allowing TU's which include the system GL and then glew (from TGLIncludes)
// leads to gltypes.h redefinition errors.
#include "TGLIncludes.h"

#include "ROOTOpenGLView.h"
#include "CocoaConstants.h"
#include "TMacOSXSystem.h"
#include "CocoaPrivate.h"
#include "QuartzWindow.h"
#include "QuartzPixmap.h"
#include "QuartzUtils.h"
#include "X11Drawable.h"
#include "QuartzText.h"
#include "CocoaUtils.h"
#include "MenuLoader.h"
#include "TVirtualGL.h"
#include "X11Events.h"
#include "X11Buffer.h"
#include "TGClient.h"
#include "TGWindow.h"
#include "TSystem.h"
#include "TGFrame.h"
#include "TGLIncludes.h"
#include "TError.h"
#include "TColor.h"
#include "TROOT.h"
#include "TEnv.h"
#include "TVirtualMutex.h"

#include <ApplicationServices/ApplicationServices.h>
#include <OpenGL/OpenGL.h>
#include <Cocoa/Cocoa.h>

#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <cstring>
#include <cstddef>
#include <limits>


//Style notes: I'm using a lot of asserts to check pre-conditions - mainly function parameters.
//In asserts, expression always looks like 'p != 0' for "C++ pointer" (either object of built-in type
//or C++ class), and 'p != nil' for object from Objective-C. There is no difference, this is to make
//asserts more explicit. In conditional statement, it'll always be 'if (p)'  or 'if (!p)' for both
//C++ and Objective-C pointers/code.

//I never use const qualifier for pointers to Objective-C objects since they are useless:
//there are no cv-qualified methods (member-functions in C++) in Objective-C, and I do not use
//'->' operator to access instance variables (data-members in C++) of Objective-C's object.
//I also declare pointer as a const, if it's const:
//NSWindow * const topLevelWindow = ... (and note, not pointer to const - no use with Obj-C).

//Asserts on drawables ids usually only check, that it's not a 'root' window id (unless operation
//is permitted on a 'root' window):
//a) assert(!fPimpl->IsRootWindow(windowID)) and later I also check that windowID != 0 (kNone).
//b) assert(drawableID > fPimpl->GetRootWindowID()) so drawableID can not be kNone and
//   can not be a 'root' window.

//ROOT window has id 1. So if id > 1 (id > fPimpl->GetRootWindowID())
//id is considered as valid (if it's out of range and > maximum valid id, this will be
//caught by CocoaPrivate.

namespace Details = ROOT::MacOSX::Details;
namespace Util = ROOT::MacOSX::Util;
namespace X11 = ROOT::MacOSX::X11;
namespace Quartz = ROOT::Quartz;
namespace OpenGL = ROOT::MacOSX::OpenGL;

namespace {

#pragma mark - Display configuration management.

//______________________________________________________________________________
void DisplayReconfigurationCallback(CGDirectDisplayID /*display*/, CGDisplayChangeSummaryFlags flags, void * /*userInfo*/)
{
   if (flags & kCGDisplayBeginConfigurationFlag)
      return;

   if (flags & kCGDisplayDesktopShapeChangedFlag) {
      assert(dynamic_cast<TGCocoa *>(gVirtualX) != 0 && "DisplayReconfigurationCallback, gVirtualX"
                                                        " is either null or has a wrong type");
      TGCocoa * const gCocoa = static_cast<TGCocoa *>(gVirtualX);
      gCocoa->ReconfigureDisplay();
   }
}

#pragma mark - Aux. functions called from GUI-rendering part.

//______________________________________________________________________________
void SetStrokeForegroundColorFromX11Context(CGContextRef ctx, const GCValues_t &gcVals)
{
   assert(ctx != 0 && "SetStrokeForegroundColorFromX11Context, parameter 'ctx' is null");

   CGFloat rgb[3] = {};
   if (gcVals.fMask & kGCForeground)
      X11::PixelToRGB(gcVals.fForeground, rgb);
   else
      ::Warning("SetStrokeForegroundColorFromX11Context",
                "x11 context does not have line color information");

   CGContextSetRGBStrokeColor(ctx, rgb[0], rgb[1], rgb[2], 1.);
}

//______________________________________________________________________________
void SetStrokeDashFromX11Context(CGContextRef ctx, const GCValues_t &gcVals)
{
   //Set line dash pattern (X11's LineOnOffDash line style).
   assert(ctx != 0 && "SetStrokeDashFromX11Context, ctx parameter is null");

   SetStrokeForegroundColorFromX11Context(ctx, gcVals);

   static const std::size_t maxLength = sizeof gcVals.fDashes / sizeof gcVals.fDashes[0];
   assert(maxLength >= std::size_t(gcVals.fDashLen) &&
          "SetStrokeDashFromX11Context, x11 context has bad dash length > sizeof(fDashes)");

   CGFloat dashes[maxLength] = {};
   for (Int_t i = 0; i < gcVals.fDashLen; ++i)
      dashes[i] = gcVals.fDashes[i];

   CGContextSetLineDash(ctx, gcVals.fDashOffset, dashes, gcVals.fDashLen);
}

//______________________________________________________________________________
void SetStrokeDoubleDashFromX11Context(CGContextRef /*ctx*/, const GCValues_t & /*gcVals*/)
{
   //assert(ctx != 0 && "SetStrokeDoubleDashFromX11Context, ctx parameter is null");
   ::Warning("SetStrokeDoubleDashFromX11Context", "Not implemented yet, kick tpochep!");
}

//______________________________________________________________________________
void SetStrokeParametersFromX11Context(CGContextRef ctx, const GCValues_t &gcVals)
{
   //Set line width and color from GCValues_t object.
   //(GUI rendering).
   assert(ctx != 0 && "SetStrokeParametersFromX11Context, parameter 'ctx' is null");

   const Mask_t mask = gcVals.fMask;
   if ((mask & kGCLineWidth) && gcVals.fLineWidth > 1)
      CGContextSetLineWidth(ctx, gcVals.fLineWidth);
   else
      CGContextSetLineWidth(ctx, 1.);

   CGContextSetLineDash(ctx, 0., 0, 0);

   if (mask & kGCLineStyle) {
      if (gcVals.fLineStyle == kLineSolid)
         SetStrokeForegroundColorFromX11Context(ctx, gcVals);
      else if (gcVals.fLineStyle == kLineOnOffDash)
         SetStrokeDashFromX11Context(ctx, gcVals);
      else if (gcVals.fLineStyle == kLineDoubleDash)
         SetStrokeDoubleDashFromX11Context(ctx ,gcVals);
      else {
         ::Warning("SetStrokeParametersFromX11Context", "line style bit is set,"
                                                        " but line style is unknown");
         SetStrokeForegroundColorFromX11Context(ctx, gcVals);
      }
   } else
      SetStrokeForegroundColorFromX11Context(ctx, gcVals);
}

//______________________________________________________________________________
void SetFilledAreaColorFromX11Context(CGContextRef ctx, const GCValues_t &gcVals)
{
   //Set fill color from "foreground" pixel color.
   //(GUI rendering).
   assert(ctx != 0 && "SetFilledAreaColorFromX11Context, parameter 'ctx' is null");

   CGFloat rgb[3] = {};
   if (gcVals.fMask & kGCForeground)
      X11::PixelToRGB(gcVals.fForeground, rgb);
   else
      ::Warning("SetFilledAreaColorFromX11Context", "no fill color found in x11 context");

   CGContextSetRGBFillColor(ctx, rgb[0], rgb[1], rgb[2], 1.);
}

struct PatternContext {
   Mask_t                 fMask;
   Int_t                  fFillStyle;
   ULong_t                fForeground;
   ULong_t                fBackground;
   NSObject<X11Drawable> *fImage;//Either stipple or tile image.
   CGSize                 fPhase;
};


//______________________________________________________________________________
bool HasFillTiledStyle(Mask_t mask, Int_t fillStyle)
{
   return (mask & kGCFillStyle) && (fillStyle == kFillTiled);
}

//______________________________________________________________________________
bool HasFillTiledStyle(const GCValues_t &gcVals)
{
   return HasFillTiledStyle(gcVals.fMask, gcVals.fFillStyle);
}

//______________________________________________________________________________
bool HasFillStippledStyle(Mask_t mask, Int_t fillStyle)
{
   return (mask & kGCFillStyle) && (fillStyle == kFillStippled);
}

//______________________________________________________________________________
bool HasFillStippledStyle(const GCValues_t &gcVals)
{
   return HasFillStippledStyle(gcVals.fMask, gcVals.fFillStyle);
}

//______________________________________________________________________________
bool HasFillOpaqueStippledStyle(Mask_t mask, Int_t fillStyle)
{
   return (mask & kGCFillStyle) && (fillStyle == kFillOpaqueStippled);
}

//______________________________________________________________________________
bool HasFillOpaqueStippledStyle(const GCValues_t &gcVals)
{
   return HasFillOpaqueStippledStyle(gcVals.fMask, gcVals.fFillStyle);
}

//______________________________________________________________________________
void DrawTile(NSObject<X11Drawable> *patternImage, CGContextRef ctx)
{
   assert(patternImage != nil && "DrawTile, parameter 'patternImage' is nil");
   assert(ctx != 0 && "DrawTile, ctx parameter is null");

   const CGRect patternRect = CGRectMake(0, 0, patternImage.fWidth, patternImage.fHeight);
   if ([patternImage isKindOfClass : [QuartzImage class]]) {
      CGContextDrawImage(ctx, patternRect, ((QuartzImage *)patternImage).fImage);
   } else if ([patternImage isKindOfClass : [QuartzPixmap class]]){
      const Util::CFScopeGuard<CGImageRef> imageFromPixmap([((QuartzPixmap *)patternImage) createImageFromPixmap]);
      assert(imageFromPixmap.Get() != 0 && "DrawTile, createImageFromPixmap failed");
      CGContextDrawImage(ctx, patternRect, imageFromPixmap.Get());
   } else
      assert(0 && "DrawTile, pattern is neither a QuartzImage, nor a QuartzPixmap");
}

//______________________________________________________________________________
void DrawPattern(void *info, CGContextRef ctx)
{
   //Pattern callback, either use foreground (and background, if any)
   //color and stipple mask to draw a pattern, or use pixmap
   //as a pattern image.
   //(GUI rendering).
   assert(info != 0 && "DrawPattern, parameter 'info' is null");
   assert(ctx != 0 && "DrawPattern, parameter 'ctx' is null");

   const PatternContext * const patternContext = (PatternContext *)info;
   const Mask_t mask = patternContext->fMask;
   const Int_t fillStyle = patternContext->fFillStyle;

   NSObject<X11Drawable> * const patternImage = patternContext->fImage;
   assert(patternImage != nil && "DrawPattern, pattern (stipple) image is nil");
   const CGRect patternRect = CGRectMake(0, 0, patternImage.fWidth, patternImage.fHeight);

   if (HasFillTiledStyle(mask, fillStyle)) {
      DrawTile(patternImage, ctx);
   } else if (HasFillStippledStyle(mask, fillStyle) || HasFillOpaqueStippledStyle(mask, fillStyle)) {
      assert([patternImage isKindOfClass : [QuartzImage class]] &&
             "DrawPattern, stipple must be a QuartzImage object");
      QuartzImage * const image = (QuartzImage *)patternImage;
      assert(image.fIsStippleMask == YES && "DrawPattern, image is not a stipple mask");

      CGFloat rgb[3] = {};

      if (HasFillOpaqueStippledStyle(mask,fillStyle)) {
         //Fill background first.
         assert((mask & kGCBackground) &&
                "DrawPattern, fill style is FillOpaqueStippled, but background color is not set in a context");
         X11::PixelToRGB(patternContext->fBackground, rgb);
         CGContextSetRGBFillColor(ctx, rgb[0], rgb[1], rgb[2], 1.);
         CGContextFillRect(ctx, patternRect);
      }

      //Fill rectangle with foreground colour, using stipple mask.
      assert((mask & kGCForeground) && "DrawPattern, foreground color is not set");
      X11::PixelToRGB(patternContext->fForeground, rgb);
      CGContextSetRGBFillColor(ctx, rgb[0], rgb[1], rgb[2], 1.);
      CGContextClipToMask(ctx, patternRect, image.fImage);
      CGContextFillRect(ctx, patternRect);
   } else {
      //This can be a window background pixmap
      DrawTile(patternImage, ctx);
   }
}

//______________________________________________________________________________
void SetFillPattern(CGContextRef ctx, const PatternContext *patternContext)
{
   //Create CGPatternRef to fill GUI elements with pattern.
   //Pattern is a QuartzImage object, it can be either a mask,
   //or pattern image itself.
   //(GUI-rendering).
   assert(ctx != 0 && "SetFillPattern, parameter 'ctx' is null");
   assert(patternContext != 0 && "SetFillPattern, parameter 'patternContext' is null");
   assert(patternContext->fImage != nil && "SetFillPattern, pattern image is nil");

   const Util::CFScopeGuard<CGColorSpaceRef> patternColorSpace(CGColorSpaceCreatePattern(0));
   CGContextSetFillColorSpace(ctx, patternColorSpace.Get());

   CGPatternCallbacks callbacks = {};
   callbacks.drawPattern = DrawPattern;
   const CGRect patternRect = CGRectMake(0, 0, patternContext->fImage.fWidth, patternContext->fImage.fHeight);
   const Util::CFScopeGuard<CGPatternRef> pattern(CGPatternCreate((void *)patternContext, patternRect, CGAffineTransformIdentity,
                                                                  patternContext->fImage.fWidth, patternContext->fImage.fHeight,
                                                                  kCGPatternTilingNoDistortion, true, &callbacks));
   const CGFloat alpha = 1.;
   CGContextSetFillPattern(ctx, pattern.Get(), &alpha);
   CGContextSetPatternPhase(ctx, patternContext->fPhase);
}

//______________________________________________________________________________
bool ParentRendersToChild(NSView<X11Window> *child)
{
   assert(child != nil && "ParentRendersToChild, parameter 'child' is nil");

   //Adovo poluchaetsia, tashhem-ta! ;)
   return (X11::ViewIsTextViewFrame(child, true) || X11::ViewIsHtmlViewFrame(child, true)) && !child.fContext &&
           child.fMapState == kIsViewable && child.fParentView.fContext &&
           !child.fIsOverlapped;
}

class ViewFixer final {
public:
   ViewFixer(QuartzView *&viewToFix, Drawable_t &widToFix)
   {
      if (ParentRendersToChild(viewToFix) && [viewToFix.fParentView isKindOfClass:[QuartzView class]]) {
         const auto origin = viewToFix.frame.origin;
         viewToFix = viewToFix.fParentView;
         widToFix = viewToFix.fID;
         if ((context = viewToFix.fContext)) {
            CGContextSaveGState(context);
            CGContextTranslateCTM(context, origin.x, origin.y);
         }
      }
   }
   ~ViewFixer()
   {
       if (context)
          CGContextRestoreGState(context);
   }
   ViewFixer(const ViewFixer &rhs) = delete;
   ViewFixer &operator = (const ViewFixer &) = delete;

private:
   CGContextRef context = nullptr;
};

//______________________________________________________________________________
bool IsNonPrintableAsciiCharacter(UniChar c)
{
   if (c == 9 || (c >= 32 && c < 127))
      return false;

   return true;
}

//______________________________________________________________________________
void FixAscii(std::vector<UniChar> &text)
{
   //GUI text is essentially ASCII. Our GUI
   //calculates text metrix 'per-symbol', this means,
   //it never asks about 'Text' metrics, but 'T', 'e', 'x', 't'.
   //Obviously, text does not fit any widget because of
   //this and I have to place all glyphs manually.
   //And here I have another problem from our GUI - it
   //can easily feed TGCocoa with non-printable symbols
   //(this is a bug). Obviously, I do not have glyphs for, say, form feed
   //or 'data link escape'. So I have to fix ascii text before
   //manual glyph rendering: DLE symbol - replaced by space (this
   //is done in TGText, but due to a bug it fails to replace them all)
   //Other non-printable symbols simply removed (and thus ignored).

   //Replace remaining ^P symbols with whitespaces, I have not idea why
   //TGTextView replaces only part of them and not all of them.
   std::replace(text.begin(), text.end(), UniChar(16), UniChar(' '));

   //Now, remove remaining non-printable characters (no glyphs exist for them).
   text.erase(std::remove_if(text.begin(), text.end(), IsNonPrintableAsciiCharacter), text.end());
}

}

ClassImp(TGCocoa)

Atom_t TGCocoa::fgDeleteWindowAtom = 0;

//______________________________________________________________________________
TGCocoa::TGCocoa()
            : fSelectedDrawable(0),
              fCocoaDraw(0),
              fDrawMode(kCopy),
              fDirectDraw(false),
              fForegroundProcess(false),
              fSetApp(true),
              fDisplayShapeChanged(true)
{
   assert(dynamic_cast<TMacOSXSystem *>(gSystem) != nullptr &&
          "TGCocoa, gSystem is eihter null or has a wrong type");
   TMacOSXSystem * const system = (TMacOSXSystem *)gSystem;

   if (!system->CocoaInitialized())
      system->InitializeCocoa();

   fPimpl.reset(new Details::CocoaPrivate);

   X11::InitWithPredefinedAtoms(fNameToAtom, fAtomToName);
   fgDeleteWindowAtom = FindAtom("WM_DELETE_WINDOW", true);

   CGDisplayRegisterReconfigurationCallback (DisplayReconfigurationCallback, 0);
}

//______________________________________________________________________________
TGCocoa::TGCocoa(const char *name, const char *title)
            : TVirtualX(name, title),
              fSelectedDrawable(0),
              fCocoaDraw(0),
              fDrawMode(kCopy),
              fDirectDraw(false),
              fForegroundProcess(false),
              fSetApp(true),
              fDisplayShapeChanged(true)
{
   assert(dynamic_cast<TMacOSXSystem *>(gSystem) != nullptr &&
          "TGCocoa, gSystem is eihter null or has a wrong type");
   TMacOSXSystem * const system = (TMacOSXSystem *)gSystem;

   if (!system->CocoaInitialized())
      system->InitializeCocoa();

   fPimpl.reset(new Details::CocoaPrivate);

   X11::InitWithPredefinedAtoms(fNameToAtom, fAtomToName);
   fgDeleteWindowAtom = FindAtom("WM_DELETE_WINDOW", true);

   CGDisplayRegisterReconfigurationCallback (DisplayReconfigurationCallback, 0);
}

//______________________________________________________________________________
TGCocoa::~TGCocoa()
{
   //
   CGDisplayRemoveReconfigurationCallback (DisplayReconfigurationCallback, 0);
}

//General part (empty, since it's not an X server.

//______________________________________________________________________________
Bool_t TGCocoa::Init(void * /*display*/)
{
   //Nothing to initialize here, return true to make
   //a caller happy.
   return kTRUE;
}


//______________________________________________________________________________
Int_t TGCocoa::OpenDisplay(const char * /*dpyName*/)
{
   //Noop.
   return 0;
}

//______________________________________________________________________________
const char *TGCocoa::DisplayName(const char *)
{
   //Noop.
   return "dummy";
}

//______________________________________________________________________________
Int_t TGCocoa::SupportsExtension(const char *) const
{
   //No, thank you, I'm not supporting any of X11 extensions!
   return -1;
}

//______________________________________________________________________________
void TGCocoa::CloseDisplay()
{
   //Noop.
}

//______________________________________________________________________________
Display_t TGCocoa::GetDisplay() const
{
   //Noop.
   return 0;
}

//______________________________________________________________________________
Visual_t TGCocoa::GetVisual() const
{
   //Noop.
   return 0;
}

//______________________________________________________________________________
Int_t TGCocoa::GetScreen() const
{
   //Noop.
   return 0;
}

//______________________________________________________________________________
UInt_t TGCocoa::ScreenWidthMM() const
{
   //Comment from TVirtualX:
   // Returns the width of the screen in millimeters.
   //End of comment.

   return CGDisplayScreenSize(CGMainDisplayID()).width;
}

//______________________________________________________________________________
Int_t TGCocoa::GetDepth() const
{
   //Comment from TVirtualX:
   // Returns depth of screen (number of bit planes).
   // Equivalent to GetPlanes().
   //End of comment.

   NSArray * const screens = [NSScreen screens];
   assert(screens != nil && "screens array is nil");

   NSScreen * const mainScreen = [screens objectAtIndex : 0];
   assert(mainScreen != nil && "screen with index 0 is nil");

   return NSBitsPerPixelFromDepth([mainScreen depth]);
}

//______________________________________________________________________________
void TGCocoa::Update(Int_t mode)
{
   R__LOCKGUARD(gROOTMutex);

   if (mode == 2) {
      assert(gClient != 0 && "Update, gClient is null");
      gClient->DoRedraw();//Call DoRedraw for all widgets, who need to be updated.
   } else if (mode > 0) {
      //Execute buffered commands.
      fPimpl->fX11CommandBuffer.Flush(fPimpl.get());
   }

   if (fDirectDraw && mode != 2)
      fPimpl->fX11CommandBuffer.FlushXOROps(fPimpl.get());
}

//______________________________________________________________________________
void TGCocoa::ReconfigureDisplay()
{
   fDisplayShapeChanged = true;
}

//______________________________________________________________________________
X11::Rectangle TGCocoa::GetDisplayGeometry()const
{
   if (fDisplayShapeChanged) {
      NSArray * const screens = [NSScreen screens];
      assert(screens != nil && screens.count != 0 && "GetDisplayGeometry, no screens found");

      NSRect frame = [(NSScreen *)[screens objectAtIndex : 0] frame];
      CGFloat xMin = frame.origin.x, xMax = xMin + frame.size.width;
      CGFloat yMin = frame.origin.y, yMax = yMin + frame.size.height;

      for (NSUInteger i = 1, e = screens.count; i < e; ++i) {
         frame = [(NSScreen *)[screens objectAtIndex : i] frame];
         xMin = std::min(xMin, frame.origin.x);
         xMax = std::max(xMax, frame.origin.x + frame.size.width);
         yMin = std::min(yMin, frame.origin.y);
         yMax = std::max(yMax, frame.origin.y + frame.size.height);
      }

      fDisplayRect.fX = int(xMin);
      fDisplayRect.fY = int(yMin);
      fDisplayRect.fWidth = unsigned(xMax - xMin);
      fDisplayRect.fHeight = unsigned(yMax - yMin);

      fDisplayShapeChanged = false;
   }

   return fDisplayRect;
}

#pragma mark - Window management part.

//______________________________________________________________________________
Window_t TGCocoa::GetDefaultRootWindow() const
{
   //Index, fixed and used only by 'root' window.
   return fPimpl->GetRootWindowID();
}

//______________________________________________________________________________
Int_t TGCocoa::InitWindow(ULong_t parentID)
{
   //InitWindow is a bad name, since this function
   //creates a window, but this name comes from the TVirtualX interface.
   //Actually, there is no special need in this function,
   //it's a kind of simplified CreateWindow (with only
   //one parameter). This function is called by TRootCanvas,
   //to create a special window inside TGCanvas (thus parentID must be a valid window ID).
   //TGX11/TGWin32 have internal array of such special windows,
   //they return index into this array, instead of drawable's ids.
   //I simply re-use CreateWindow and return a drawable's id.

   assert(parentID != 0 && "InitWindow, parameter 'parentID' is 0");

   //Use parent's attributes (as it's done in TGX11).
   WindowAttributes_t attr = {};
   if (fPimpl->IsRootWindow(parentID))
      ROOT::MacOSX::X11::GetRootWindowAttributes(&attr);
   else
      [fPimpl->GetWindow(parentID) getAttributes : &attr];

   return CreateWindow(parentID, 0, 0, attr.fWidth, attr.fHeight, 0, attr.fDepth, attr.fClass, 0, 0, 0);
}

//______________________________________________________________________________
Window_t TGCocoa::GetWindowID(Int_t windowID)
{
   //In case of TGX11/TGWin32, there is a mixture of
   //casted X11 ids (Window_t) and indices in some internal array, which
   //contains such an id. On Mac I always have indices. Yes, I'm smart.
   return windowID;
}

//______________________________________________________________________________
void TGCocoa::SelectWindow(Int_t windowID)
{
   //This function can be called from pad/canvas, both for window and for pixmap.
   fSelectedDrawable = windowID;
}

//______________________________________________________________________________
void TGCocoa::ClearWindow()
{
   //Clear the selected drawable OR pixmap (the name - from TVirtualX interface - is bad).
   assert(fSelectedDrawable > fPimpl->GetRootWindowID() &&
          "ClearWindow, fSelectedDrawable is invalid");

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(fSelectedDrawable);
   if (drawable.fIsPixmap) {
      //Pixmaps are white by default.
      //This is bad - we can not have transparent sub-pads (in TCanvas)
      //because of this. But there is no way how gVirtualX can
      //obtain real pad's color and check for its transparency.
      CGContextRef pixmapCtx = drawable.fContext;
      assert(pixmapCtx != 0 && "ClearWindow, pixmap's context is null");
      //const Quartz::CGStateGuard ctxGuard(pixmapCtx);
      //CGContextSetRGBFillColor(pixmapCtx, 1., 1., 1., 1.);
      //CGContextFillRect(pixmapCtx, CGRectMake(0, 0, drawable.fWidth, drawable.fHeight));
      //Now we really clear!
      CGContextClearRect(pixmapCtx, CGRectMake(0, 0, drawable.fWidth, drawable.fHeight));
   } else {
      //For a window ClearArea with w == 0 and h == 0 means the whole window.
      ClearArea(fSelectedDrawable, 0, 0, 0, 0);
   }
}

//______________________________________________________________________________
void TGCocoa::GetGeometry(Int_t windowID, Int_t & x, Int_t &y, UInt_t &w, UInt_t &h)
{
   //In TGX11, GetGeometry works with special windows, created by InitWindow
   //(thus this function is called from TCanvas/TGCanvas/TRootCanvas).

   //IMPORTANT: this function also translates x and y
   //from parent's coordinates into screen coordinates - so, again, name "GetGeometry"
   //from the TVirtualX interface is bad and misleading.

   if (windowID < 0 || fPimpl->IsRootWindow(windowID)) {
      //Comment in TVirtualX suggests, that wid can be < 0.
      //This will be a screen's geometry.
      WindowAttributes_t attr = {};
      ROOT::MacOSX::X11::GetRootWindowAttributes(&attr);
      x = attr.fX;
      y = attr.fY;
      w = attr.fWidth;
      h = attr.fHeight;
   } else {
      NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(windowID);
      x = drawable.fX;
      y = drawable.fY;
      w = drawable.fWidth;
      h = drawable.fHeight;

      if (!drawable.fIsPixmap) {
         NSObject<X11Window> * const window = (NSObject<X11Window> *)drawable;
         NSPoint srcPoint = {};
         srcPoint.x = x;
         srcPoint.y = y;
         NSView<X11Window> * const view = window.fContentView.fParentView ? window.fContentView.fParentView : window.fContentView;
         //View parameter for TranslateToScreen call must
         //be parent view, since x and y are in parent's
         //coordinate system.
         const NSPoint dstPoint = X11::TranslateToScreen(view, srcPoint);
         x = dstPoint.x;
         y = dstPoint.y;
      }
   }
}

//______________________________________________________________________________
void TGCocoa::MoveWindow(Int_t windowID, Int_t x, Int_t y)
{
   //windowID is either kNone or a valid window id.
   //x and y are coordinates of a top-left corner relative to the parent's coordinate system.

   assert(!fPimpl->IsRootWindow(windowID) && "MoveWindow, called for root window");

   if (!windowID)//From TGX11.
      return;

   [fPimpl->GetWindow(windowID) setX : x Y : y];
}

//______________________________________________________________________________
void TGCocoa::RescaleWindow(Int_t /*wid*/, UInt_t /*w*/, UInt_t /*h*/)
{
   //This function is for TRootCanvas and related stuff, never gets
   //called/used from/by any our GUI class.
   //Noop.
}

//______________________________________________________________________________
void TGCocoa::ResizeWindow(Int_t windowID)
{
   //This function does not resize window (it was done already by layout management?),
   //it resizes "back buffer" if any.

   if (!windowID)//From TGX11.
      return;

   assert(!fPimpl->IsRootWindow(windowID) &&
          "ResizeWindow, parameter 'windowID' is a root window's id");

   const Util::AutoreleasePool pool;

   NSObject<X11Window> * const window = fPimpl->GetWindow(windowID);
   if (window.fBackBuffer) {
      const Drawable_t currentDrawable = fSelectedDrawable;
      fSelectedDrawable = windowID;
      SetDoubleBufferON();
      fSelectedDrawable = currentDrawable;
   }
}

//______________________________________________________________________________
void TGCocoa::UpdateWindow(Int_t /*mode*/)
{
   //This function is used by TCanvas/TPad:
   //draw "back buffer" image into the view.
   //fContentView (destination) MUST be a QuartzView.

   //Basic es-guarantee: X11Buffer::AddUpdateWindow modifies vector with commands,
   //if the following call to TGCocoa::Update will produce an exception dusing X11Buffer::Flush,
   //initial state of X11Buffer can not be restored, but it still must be in some valid state.

   assert(fSelectedDrawable > fPimpl->GetRootWindowID() &&
          "UpdateWindow, fSelectedDrawable is not a valid window id");

   //Have no idea, why this can happen with ROOT - done by TGDNDManager :(
   if (fPimpl->GetDrawable(fSelectedDrawable).fIsPixmap == YES)
      return;

   NSObject<X11Window> * const window = fPimpl->GetWindow(fSelectedDrawable);

   if (QuartzPixmap * const pixmap = window.fBackBuffer) {
      assert([window.fContentView isKindOfClass : [QuartzView class]] && "UpdateWindow, content view is not a QuartzView");
      QuartzView *dstView = (QuartzView *)window.fContentView;

      if (dstView.fIsOverlapped)
         return;

      if (dstView.fContext) {
         //We can draw directly.
         const X11::Rectangle copyArea(0, 0, pixmap.fWidth, pixmap.fHeight);
         [dstView copy : pixmap area : copyArea withMask : nil clipOrigin : X11::Point() toPoint : X11::Point()];
      } else {
         //Have to wait.
         fPimpl->fX11CommandBuffer.AddUpdateWindow(dstView);
         Update(1);
      }
   }
}

//______________________________________________________________________________
Window_t TGCocoa::GetCurrentWindow() const
{
   //Window selected by SelectWindow.
   return fSelectedDrawable;
}

//______________________________________________________________________________
void TGCocoa::CloseWindow()
{
   //Deletes selected window.
}

//______________________________________________________________________________
Int_t TGCocoa::AddWindow(ULong_t /*qwid*/, UInt_t /*w*/, UInt_t /*h*/)
{
   //Should register a window created by Qt as a ROOT window,
   //but since Qt-ROOT does not work on Mac and will never work,
   //especially with version 4.8 - this implementation will always
   //be empty.
   return 0;
}

//______________________________________________________________________________
void TGCocoa::RemoveWindow(ULong_t /*qwid*/)
{
   //Remove window, created by Qt.
}

//______________________________________________________________________________
Window_t TGCocoa::CreateWindow(Window_t parentID, Int_t x, Int_t y, UInt_t w, UInt_t h, UInt_t border, Int_t depth,
                               UInt_t clss, void *visual, SetWindowAttributes_t *attr, UInt_t wtype)
{
   //Create new window (top-level == QuartzWindow + QuartzView, or child == QuartzView)

   //Strong es-guarantee - exception can be only during registration, class state will remain
   //unchanged, no leaks (scope guards).

   const Util::AutoreleasePool pool;

   if (fPimpl->IsRootWindow(parentID)) {//parent == root window.
      //Can throw:
      QuartzWindow * const newWindow = X11::CreateTopLevelWindow(x, y, w, h, border,
                                                                 depth, clss, visual, attr, wtype);
      //Something like unique_ptr would perfectly solve the problem with raw pointer + a separate
      //guard for this pointer, but it requires move semantics.
      const Util::NSScopeGuard<QuartzWindow> winGuard(newWindow);
      const Window_t result = fPimpl->RegisterDrawable(newWindow);//Can throw.
      newWindow.fID = result;
      [newWindow setAcceptsMouseMovedEvents : YES];

      return result;
   } else {
      NSObject<X11Window> * const parentWin = fPimpl->GetWindow(parentID);
      //OpenGL view can not have children.
      assert([parentWin.fContentView isKindOfClass : [QuartzView class]] &&
             "CreateWindow, parent view must be QuartzView");

      //Can throw:
      QuartzView * const childView = X11::CreateChildView((QuartzView *)parentWin.fContentView,
                                                           x, y, w, h, border, depth, clss, visual, attr, wtype);
      const Util::NSScopeGuard<QuartzView> viewGuard(childView);
      const Window_t result = fPimpl->RegisterDrawable(childView);//Can throw.
      childView.fID = result;
      [parentWin addChild : childView];

      return result;
   }
}

//______________________________________________________________________________
void TGCocoa::DestroyWindow(Window_t wid)
{
   //The XDestroyWindow function destroys the specified window as well as all of its subwindows
   //and causes the X server to generate a DestroyNotify event for each window.  The window
   //should never be referenced again.  If the window specified by the w argument is mapped,
   //it is unmapped automatically.  The ordering of the
   //DestroyNotify events is such that for any given window being destroyed, DestroyNotify is generated
   //on any inferiors of the window before being generated on the window itself.  The ordering
   //among siblings and across subhierarchies is not otherwise constrained.
   //If the window you specified is a root window, no windows are destroyed. Destroying a mapped window
   //will generate Expose events on other windows that were obscured by the window being destroyed.

   //No-throw guarantee???

   //I have NO idea why ROOT's GUI calls DestroyWindow with illegal
   //window id, but it does.

   if (!wid)
      return;

   if (fPimpl->IsRootWindow(wid))
      return;

   BOOL needFocusChange = NO;

   {//Block to force autoreleasepool to drain.
   const Util::AutoreleasePool pool;

   fPimpl->fX11EventTranslator.CheckUnmappedView(wid);

   assert(fPimpl->GetDrawable(wid).fIsPixmap == NO &&
          "DestroyWindow, can not be called for QuartzPixmap or QuartzImage object");

   NSObject<X11Window> * const window = fPimpl->GetWindow(wid);
   if (fPimpl->fX11CommandBuffer.BufferSize())
      fPimpl->fX11CommandBuffer.RemoveOperationsForDrawable(wid);

   //TEST: "fix" a keyboard focus.
   if ((needFocusChange = window == window.fQuartzWindow && window.fQuartzWindow.fHasFocus))
      window.fHasFocus = NO;//If any.

   DestroySubwindows(wid);
   if (window.fEventMask & kStructureNotifyMask)
      fPimpl->fX11EventTranslator.GenerateDestroyNotify(wid);

   //Interrupt modal loop (TGClient::WaitFor).
   if (gClient->GetWaitForEvent() == kDestroyNotify && wid == gClient->GetWaitForWindow())
      gClient->SetWaitForWindow(kNone);

   fPimpl->DeleteDrawable(wid);
   }

   //"Fix" a keyboard focus.
   if (needFocusChange)
      X11::WindowLostFocus(wid);
}

//______________________________________________________________________________
void TGCocoa::DestroySubwindows(Window_t wid)
{
   // The DestroySubwindows function destroys all inferior windows of the
   // specified window, in bottom-to-top stacking order.

   //No-throw guarantee??

   //From TGX11:
   if (!wid)
      return;

   if (fPimpl->IsRootWindow(wid))
      return;

   const Util::AutoreleasePool pool;

   assert(fPimpl->GetDrawable(wid).fIsPixmap == NO &&
          "DestroySubwindows, can not be called for QuartzPixmap or QuartzImage object");

   NSObject<X11Window> *window = fPimpl->GetWindow(wid);

   //I can not iterate on subviews array directly, since it'll be modified
   //during this iteration - create a copy (and it'll also increase references,
   //which will be decreased by guard's dtor).
   const Util::NSScopeGuard<NSArray> children([[window.fContentView subviews] copy]);

   for (NSView<X11Window> *child in children.Get())
      DestroyWindow(child.fID);
}

//______________________________________________________________________________
void TGCocoa::GetWindowAttributes(Window_t wid, WindowAttributes_t &attr)
{
   //No-throw guarantee.

   if (!wid)//X11's None?
      return;

   if (fPimpl->IsRootWindow(wid))
      ROOT::MacOSX::X11::GetRootWindowAttributes(&attr);
   else
      [fPimpl->GetWindow(wid) getAttributes : &attr];
}

//______________________________________________________________________________
void TGCocoa::ChangeWindowAttributes(Window_t wid, SetWindowAttributes_t *attr)
{
   //No-throw guarantee.

   if (!wid)//From TGX11
      return;

   const Util::AutoreleasePool pool;

   assert(!fPimpl->IsRootWindow(wid) && "ChangeWindowAttributes, called for root window");
   assert(attr != 0 && "ChangeWindowAttributes, parameter 'attr' is null");

   [fPimpl->GetWindow(wid) setAttributes : attr];
}

//______________________________________________________________________________
void TGCocoa::SelectInput(Window_t windowID, UInt_t eventMask)
{
   //No-throw guarantee.

   // Defines which input events the window is interested in. By default
   // events are propageted up the window stack. This mask can also be
   // set at window creation time via the SetWindowAttributes_t::fEventMask
   // attribute.

   //TGuiBldDragManager selects input on a 'root' window.
   //TGWin32 has a check on windowID == 0.
   if (windowID <= fPimpl->GetRootWindowID())
      return;

   NSObject<X11Window> * const window = fPimpl->GetWindow(windowID);
   //XSelectInput overrides a previous mask.
   window.fEventMask = eventMask;
}

//______________________________________________________________________________
void TGCocoa::ReparentChild(Window_t wid, Window_t pid, Int_t x, Int_t y)
{
   //Reparent view.
   using namespace Details;

   assert(!fPimpl->IsRootWindow(wid) && "ReparentChild, can not re-parent root window");

   const Util::AutoreleasePool pool;

   NSView<X11Window> * const view = fPimpl->GetWindow(wid).fContentView;
   if (fPimpl->IsRootWindow(pid)) {
      //Make a top-level view from a child view.
      [view retain];
      [view removeFromSuperview];
      view.fParentView = nil;

      NSRect frame = view.frame;
      frame.origin = NSPoint();

      NSUInteger styleMask = kClosableWindowMask | kMiniaturizableWindowMask | kResizableWindowMask;
      if (!view.fOverrideRedirect)
         styleMask |= kTitledWindowMask;

      QuartzWindow * const newTopLevel = [[QuartzWindow alloc] initWithContentRect : frame
                                                                         styleMask : styleMask
                                                                           backing : NSBackingStoreBuffered
                                                                             defer : NO];
      [view setX : x Y : y];
      [newTopLevel addChild : view];

      fPimpl->ReplaceDrawable(wid, newTopLevel);

      [view release];
      [newTopLevel release];
   } else {
      [view retain];
      [view removeFromSuperview];
      //
      NSObject<X11Window> * const newParent = fPimpl->GetWindow(pid);
      assert(newParent.fIsPixmap == NO && "ReparentChild, pixmap can not be a new parent");
      [view setX : x Y : y];
      [newParent addChild : view];//It'll also update view's level, no need to call updateLevel.
      [view release];
   }
}

//______________________________________________________________________________
void TGCocoa::ReparentTopLevel(Window_t wid, Window_t pid, Int_t x, Int_t y)
{
   //Reparent top-level window.
   //I have to delete QuartzWindow here and place in its slot content view +
   //reparent this view into pid.
   if (fPimpl->IsRootWindow(pid))//Nothing to do, wid is already a top-level window.
      return;

   const Util::AutoreleasePool pool;

   NSView<X11Window> * const contentView = fPimpl->GetWindow(wid).fContentView;
   QuartzWindow * const topLevel = (QuartzWindow *)[contentView window];
   [contentView retain];
   [contentView removeFromSuperview];
   [topLevel setContentView : nil];
   fPimpl->ReplaceDrawable(wid, contentView);
   [contentView setX : x Y : y];
   [fPimpl->GetWindow(pid) addChild : contentView];//Will also replace view's level.
   [contentView release];
}

//______________________________________________________________________________
void TGCocoa::ReparentWindow(Window_t wid, Window_t pid, Int_t x, Int_t y)
{
   //Change window's parent (possibly creating new top-level window or destroying top-level window).

   if (!wid) //From TGX11.
      return;

   assert(!fPimpl->IsRootWindow(wid) && "ReparentWindow, can not re-parent root window");

   NSView<X11Window> * const view = fPimpl->GetWindow(wid).fContentView;
   if (view.fParentView)
      ReparentChild(wid, pid, x, y);
   else
      //wid is a top-level window (or content view of such a window).
      ReparentTopLevel(wid, pid, x, y);
}

//______________________________________________________________________________
void TGCocoa::MapWindow(Window_t wid)
{
   // Maps the window "wid" and all of its subwindows that have had map
   // requests. This function has no effect if the window is already mapped.

   assert(!fPimpl->IsRootWindow(wid) && "MapWindow, called for root window");

   const Util::AutoreleasePool pool;

   if (MakeProcessForeground())
      [fPimpl->GetWindow(wid) mapWindow];

   if (fSetApp) {
      SetApplicationIcon();
      Details::PopulateMainMenu();
      fSetApp = false;
   }
}

//______________________________________________________________________________
void TGCocoa::MapSubwindows(Window_t wid)
{
   // Maps all subwindows for the specified window "wid" in top-to-bottom
   // stacking order.

   assert(!fPimpl->IsRootWindow(wid) && "MapSubwindows, called for 'root' window");

   const Util::AutoreleasePool pool;

   if (MakeProcessForeground())
      [fPimpl->GetWindow(wid) mapSubwindows];
}

//______________________________________________________________________________
void TGCocoa::MapRaised(Window_t wid)
{
   // Maps the window "wid" and all of its subwindows that have had map
   // requests on the screen and put this window on the top of of the
   // stack of all windows.

   assert(!fPimpl->IsRootWindow(wid) && "MapRaised, called for root window");

   const Util::AutoreleasePool pool;

   if (MakeProcessForeground())
      [fPimpl->GetWindow(wid) mapRaised];

   if (fSetApp) {
      SetApplicationIcon();
      Details::PopulateMainMenu();
      fSetApp = false;
   }
}

//______________________________________________________________________________
void TGCocoa::UnmapWindow(Window_t wid)
{
   // Unmaps the specified window "wid". If the specified window is already
   // unmapped, this function has no effect. Any child window will no longer
   // be visible (but they are still mapped) until another map call is made
   // on the parent.
   assert(!fPimpl->IsRootWindow(wid) && "UnmapWindow, called for root window");

   const Util::AutoreleasePool pool;

   //If this window is a grab window or a parent of a grab window.
   fPimpl->fX11EventTranslator.CheckUnmappedView(wid);

   NSObject<X11Window> * const win = fPimpl->GetWindow(wid);
   [win unmapWindow];

   if (win == win.fQuartzWindow && win.fQuartzWindow.fHasFocus)
      X11::WindowLostFocus(win.fID);

   win.fHasFocus = NO;

   //if (window.fEventMask & kStructureNotifyMask)
   //   fPimpl->fX11EventTranslator.GenerateUnmapNotify(wid);

   //Interrupt modal loop (TGClient::WaitForUnmap).
   if (gClient->GetWaitForEvent() == kUnmapNotify && gClient->GetWaitForWindow() == wid)
      gClient->SetWaitForWindow(kNone);
}

//______________________________________________________________________________
void TGCocoa::RaiseWindow(Window_t wid)
{
   // Raises the specified window to the top of the stack so that no
   // sibling window obscures it.

   if (!wid)//From TGX11.
      return;

   assert(!fPimpl->IsRootWindow(wid) && "RaiseWindow, called for root window");

   if (!fPimpl->GetWindow(wid).fParentView)
      return;

   [fPimpl->GetWindow(wid) raiseWindow];
}

//______________________________________________________________________________
void TGCocoa::LowerWindow(Window_t wid)
{
   // Lowers the specified window "wid" to the bottom of the stack so
   // that it does not obscure any sibling windows.

   if (!wid)//From TGX11.
      return;

   assert(!fPimpl->IsRootWindow(wid) && "LowerWindow, called for root window");

   if (!fPimpl->GetWindow(wid).fParentView)
      return;

   [fPimpl->GetWindow(wid) lowerWindow];
}

//______________________________________________________________________________
void TGCocoa::MoveWindow(Window_t wid, Int_t x, Int_t y)
{
   // Moves the specified window to the specified x and y coordinates.
   // It does not change the window's size, raise the window, or change
   // the mapping state of the window.
   //
   // x, y - coordinates, which define the new position of the window
   //        relative to its parent.

   if (!wid)//From TGX11.
      return;

   assert(!fPimpl->IsRootWindow(wid) && "MoveWindow, called for root window");
   const Util::AutoreleasePool pool;
   [fPimpl->GetWindow(wid) setX : x Y : y];
}

//______________________________________________________________________________
void TGCocoa::MoveResizeWindow(Window_t wid, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   // Changes the size and location of the specified window "wid" without
   // raising it.
   //
   // x, y - coordinates, which define the new position of the window
   //        relative to its parent.
   // w, h - the width and height, which define the interior size of
   //        the window

   if (!wid)//From TGX11.
      return;

   assert(!fPimpl->IsRootWindow(wid) && "MoveResizeWindow, called for 'root' window");

   const Util::AutoreleasePool pool;
   [fPimpl->GetWindow(wid) setX : x Y : y width : w height : h];
}

//______________________________________________________________________________
void TGCocoa::ResizeWindow(Window_t wid, UInt_t w, UInt_t h)
{
   if (!wid)//From TGX11.
      return;

   assert(!fPimpl->IsRootWindow(wid) && "ResizeWindow, called for 'root' window");

   const Util::AutoreleasePool pool;

   //We can have this unfortunately.
   const UInt_t siMax = std::numeric_limits<Int_t>::max();
   if (w > siMax || h > siMax)
      return;

   NSSize newSize = {};
   newSize.width = w;
   newSize.height = h;

   [fPimpl->GetWindow(wid) setDrawableSize : newSize];
}

//______________________________________________________________________________
void TGCocoa::IconifyWindow(Window_t wid)
{
   // Iconifies the window "wid".
   if (!wid)
      return;

   assert(!fPimpl->IsRootWindow(wid) && "IconifyWindow, can not iconify the root window");
   assert(fPimpl->GetWindow(wid).fIsPixmap == NO && "IconifyWindow, invalid window id");

   NSObject<X11Window> * const win = fPimpl->GetWindow(wid);
   assert(win.fQuartzWindow == win && "IconifyWindow, can be called only for a top level window");

   fPimpl->fX11EventTranslator.CheckUnmappedView(wid);

   NSObject<X11Window> * const window = fPimpl->GetWindow(wid);
   if (fPimpl->fX11CommandBuffer.BufferSize())
      fPimpl->fX11CommandBuffer.RemoveOperationsForDrawable(wid);

   if (window.fQuartzWindow.fHasFocus) {
      X11::WindowLostFocus(win.fID);
      window.fQuartzWindow.fHasFocus = NO;
   }

   [win.fQuartzWindow miniaturize : win.fQuartzWindow];
}

//______________________________________________________________________________
void TGCocoa::TranslateCoordinates(Window_t srcWin, Window_t dstWin, Int_t srcX, Int_t srcY, Int_t &dstX, Int_t &dstY, Window_t &child)
{
   // Translates coordinates in one window to the coordinate space of another
   // window. It takes the "src_x" and "src_y" coordinates relative to the
   // source window's origin and returns these coordinates to "dest_x" and
   // "dest_y" relative to the destination window's origin.

   // child          - returns the child of "dest" if the coordinates
   //                  are contained in a mapped child of the destination
   //                  window; otherwise, child is set to 0
   child = 0;
   if (!srcWin || !dstWin)//This is from TGX11, looks like this can happen.
      return;

   const bool srcIsRoot = fPimpl->IsRootWindow(srcWin);
   const bool dstIsRoot = fPimpl->IsRootWindow(dstWin);

   if (srcIsRoot && dstIsRoot) {
      //This can happen with ROOT's GUI. Set dstX/Y equal to srcX/Y.
      //From man for XTranslateCoordinates it's not clear, what should be in child.
      dstX = srcX;
      dstY = srcY;

      if (QuartzWindow * const qw = X11::FindWindowInPoint(srcX, srcY))
         child = qw.fID;

      return;
   }

   NSPoint srcPoint = {};
   srcPoint.x = srcX;
   srcPoint.y = srcY;

   NSPoint dstPoint = {};


   if (dstIsRoot) {
      NSView<X11Window> * const srcView = fPimpl->GetWindow(srcWin).fContentView;
      dstPoint = X11::TranslateToScreen(srcView, srcPoint);
   } else if (srcIsRoot) {
      NSView<X11Window> * const dstView = fPimpl->GetWindow(dstWin).fContentView;
      dstPoint = X11::TranslateFromScreen(srcPoint, dstView);

      if ([dstView superview]) {
         //hitTest requires a point in a superview's coordinate system.
         //Even contentView of QuartzWindow has a superview (NSThemeFrame),
         //so this should always work.
         dstPoint = [[dstView superview] convertPoint : dstPoint fromView : dstView];
         if (NSView<X11Window> * const view = (NSView<X11Window> *)[dstView hitTest : dstPoint]) {
            if (view != dstView && view.fMapState == kIsViewable)
               child = view.fID;
         }
      }
   } else {
      NSView<X11Window> * const srcView = fPimpl->GetWindow(srcWin).fContentView;
      NSView<X11Window> * const dstView = fPimpl->GetWindow(dstWin).fContentView;

      dstPoint = X11::TranslateCoordinates(srcView, dstView, srcPoint);
      if ([dstView superview]) {
         //hitTest requires a point in a view's superview coordinate system.
         //Even contentView of QuartzWindow has a superview (NSThemeFrame),
         //so this should always work.
         const NSPoint pt = [[dstView superview] convertPoint : dstPoint fromView : dstView];
         if (NSView<X11Window> * const view = (NSView<X11Window> *)[dstView hitTest : pt]) {
            if (view != dstView && view.fMapState == kIsViewable)
               child = view.fID;
         }
      }
   }

   dstX = dstPoint.x;
   dstY = dstPoint.y;
}

//______________________________________________________________________________
void TGCocoa::GetWindowSize(Drawable_t wid, Int_t &x, Int_t &y, UInt_t &w, UInt_t &h)
{
   // Returns the location and the size of window "wid"
   //
   // x, y - coordinates of the upper-left outer corner relative to the
   //        parent window's origin
   // w, h - the size of the window, not including the border.

   //From GX11Gui.cxx:
   if (!wid)
      return;

   if (fPimpl->IsRootWindow(wid)) {
      WindowAttributes_t attr = {};
      ROOT::MacOSX::X11::GetRootWindowAttributes(&attr);
      x = attr.fX;
      y = attr.fY;
      w = attr.fWidth;
      h = attr.fHeight;
   } else {
      NSObject<X11Drawable> *window = fPimpl->GetDrawable(wid);
      //ROOT can ask window size for ... non-window drawable.
      if (!window.fIsPixmap) {
         x = window.fX;
         y = window.fY;
      } else {
         x = 0;
         y = 0;
      }

      w = window.fWidth;
      h = window.fHeight;
   }
}

//______________________________________________________________________________
void TGCocoa::SetWindowBackground(Window_t wid, ULong_t color)
{
   //From TGX11:
   if (!wid)
      return;

   assert(!fPimpl->IsRootWindow(wid) && "SetWindowBackground, can not set color for root window");

   fPimpl->GetWindow(wid).fBackgroundPixel = color;
}

//______________________________________________________________________________
void TGCocoa::SetWindowBackgroundPixmap(Window_t windowID, Pixmap_t pixmapID)
{
   // Sets the background pixmap of the window "wid" to the specified
   // pixmap "pxm".

   //From TGX11/TGWin32:
   if (!windowID)
      return;

   assert(!fPimpl->IsRootWindow(windowID) &&
          "SetWindowBackgroundPixmap, can not set background for a root window");
   assert(fPimpl->GetDrawable(windowID).fIsPixmap == NO &&
          "SetWindowBackgroundPixmap, invalid window id");

   NSObject<X11Window> * const window = fPimpl->GetWindow(windowID);
   if (pixmapID == kNone) {
      window.fBackgroundPixmap = nil;
      return;
   }

   assert(pixmapID > fPimpl->GetRootWindowID() &&
          "SetWindowBackgroundPixmap, parameter 'pixmapID' is not a valid pixmap id");
   assert(fPimpl->GetDrawable(pixmapID).fIsPixmap == YES &&
          "SetWindowBackgroundPixmap, bad drawable");

   NSObject<X11Drawable> * const pixmapOrImage = fPimpl->GetDrawable(pixmapID);
   //X11 doc says, that pixmap can be freed immediately after call
   //XSetWindowBackgroundPixmap, so I have to copy a pixmap.
   Util::NSScopeGuard<QuartzImage> backgroundImage;

   if ([pixmapOrImage isKindOfClass : [QuartzPixmap class]]) {
      backgroundImage.Reset([[QuartzImage alloc] initFromPixmap : (QuartzPixmap *)pixmapOrImage]);
      if (backgroundImage.Get())
         window.fBackgroundPixmap = backgroundImage.Get();//the window is retaining the image.
   } else {
      backgroundImage.Reset([[QuartzImage alloc] initFromImage : (QuartzImage *)pixmapOrImage]);
      if (backgroundImage.Get())
         window.fBackgroundPixmap = backgroundImage.Get();//the window is retaining the image.
   }

   if (!backgroundImage.Get())
      //Detailed error message was issued by QuartzImage at this point.
      Error("SetWindowBackgroundPixmap", "QuartzImage initialization failed");
}

//______________________________________________________________________________
Window_t TGCocoa::GetParent(Window_t windowID) const
{
   // Returns the parent of the window "windowID".

   //0 or root (checked in TGX11):
   if (windowID <= fPimpl->GetRootWindowID())
      return windowID;

   NSView<X11Window> *view = fPimpl->GetWindow(windowID).fContentView;
   return view.fParentView ? view.fParentView.fID : fPimpl->GetRootWindowID();
}

//______________________________________________________________________________
void TGCocoa::SetWindowName(Window_t wid, char *name)
{
   if (!wid || !name)//From TGX11.
      return;

   const Util::AutoreleasePool pool;

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);

   if ([(NSObject *)drawable isKindOfClass : [NSWindow class]]) {
      NSString * const windowTitle = [NSString stringWithCString : name encoding : NSASCIIStringEncoding];
      [(NSWindow *)drawable setTitle : windowTitle];
   }
}

//______________________________________________________________________________
void TGCocoa::SetIconName(Window_t /*wid*/, char * /*name*/)
{
   //Noop.
}

//______________________________________________________________________________
void TGCocoa::SetIconPixmap(Window_t /*wid*/, Pixmap_t /*pix*/)
{
   //Noop.
}

//______________________________________________________________________________
void TGCocoa::SetClassHints(Window_t /*wid*/, char * /*className*/, char * /*resourceName*/)
{
   //Noop.
}

//______________________________________________________________________________
void TGCocoa::ShapeCombineMask(Window_t windowID, Int_t /*x*/, Int_t /*y*/, Pixmap_t pixmapID)
{
   //Comment from TVirtualX:
   // The Nonrectangular Window Shape Extension adds nonrectangular
   // windows to the System.
   // This allows for making shaped (partially transparent) windows

   assert(!fPimpl->IsRootWindow(windowID) &&
          "ShapeCombineMask, windowID parameter is a 'root' window");
   assert(fPimpl->GetDrawable(windowID).fIsPixmap == NO &&
          "ShapeCombineMask, windowID parameter is a bad window id");
   assert([fPimpl->GetDrawable(pixmapID) isKindOfClass : [QuartzImage class]] &&
          "ShapeCombineMask, pixmapID parameter must point to QuartzImage object");

   if (fPimpl->GetWindow(windowID).fContentView.fParentView)
      return;

   QuartzImage * const srcImage = (QuartzImage *)fPimpl->GetDrawable(pixmapID);
   assert(srcImage.fIsStippleMask == YES && "ShapeCombineMask, source image is not a stipple mask");

   // There is some kind of problems with shape masks and
   // flipped views, I have to do an image flip here.
   const Util::NSScopeGuard<QuartzImage> image([[QuartzImage alloc] initFromImageFlipped : srcImage]);
   if (image.Get()) {
      QuartzWindow * const qw = fPimpl->GetWindow(windowID).fQuartzWindow;
      qw.fShapeCombineMask = image.Get();
      [qw setOpaque : NO];
      [qw setBackgroundColor : [NSColor clearColor]];
   }
}

#pragma mark - "Window manager hints" set of functions.

//______________________________________________________________________________
void TGCocoa::SetMWMHints(Window_t wid, UInt_t value, UInt_t funcs, UInt_t /*input*/)
{
   // Sets decoration style.
   using namespace Details;

   assert(!fPimpl->IsRootWindow(wid) && "SetMWMHints, called for 'root' window");

   QuartzWindow * const qw = fPimpl->GetWindow(wid).fQuartzWindow;
   NSUInteger newMask = 0;

   if ([qw styleMask] & kTitledWindowMask) {//Do not modify this.
      newMask |= kTitledWindowMask;
      newMask |= kClosableWindowMask;
   }

   if (value & kMWMFuncAll) {
      newMask |= kMiniaturizableWindowMask | kResizableWindowMask;
   } else {
      if (value & kMWMDecorMinimize)
         newMask |= kMiniaturizableWindowMask;
      if (funcs & kMWMFuncResize)
         newMask |= kResizableWindowMask;
   }

   [qw setStyleMask : newMask];

   if (funcs & kMWMDecorAll) {
      if (!qw.fMainWindow) {//Do not touch buttons for transient window.
         [[qw standardWindowButton : NSWindowZoomButton] setEnabled : YES];
         [[qw standardWindowButton : NSWindowMiniaturizeButton] setEnabled : YES];
      }
   } else {
      if (!qw.fMainWindow) {//Do not touch transient window's titlebar.
         [[qw standardWindowButton : NSWindowZoomButton] setEnabled : funcs & kMWMDecorMaximize];
         [[qw standardWindowButton : NSWindowMiniaturizeButton] setEnabled : funcs & kMWMDecorMinimize];
      }
   }
}

//______________________________________________________________________________
void TGCocoa::SetWMPosition(Window_t /*wid*/, Int_t /*x*/, Int_t /*y*/)
{
   //Noop.
}

//______________________________________________________________________________
void TGCocoa::SetWMSize(Window_t /*wid*/, UInt_t /*w*/, UInt_t /*h*/)
{
   //Noop.
}

//______________________________________________________________________________
void TGCocoa::SetWMSizeHints(Window_t wid, UInt_t wMin, UInt_t hMin, UInt_t wMax, UInt_t hMax, UInt_t /*wInc*/, UInt_t /*hInc*/)
{
   using namespace Details;

   assert(!fPimpl->IsRootWindow(wid) && "SetWMSizeHints, called for root window");

   const NSUInteger styleMask = kTitledWindowMask | kClosableWindowMask | kMiniaturizableWindowMask | kResizableWindowMask;
   const NSRect minRect = [NSWindow frameRectForContentRect : NSMakeRect(0., 0., wMin, hMin) styleMask : styleMask];
   const NSRect maxRect = [NSWindow frameRectForContentRect : NSMakeRect(0., 0., wMax, hMax) styleMask : styleMask];

   QuartzWindow * const qw = fPimpl->GetWindow(wid).fQuartzWindow;
   [qw setMinSize : minRect.size];
   [qw setMaxSize : maxRect.size];
}

//______________________________________________________________________________
void TGCocoa::SetWMState(Window_t /*wid*/, EInitialState /*state*/)
{
   //Noop.
}

//______________________________________________________________________________
void TGCocoa::SetWMTransientHint(Window_t wid, Window_t mainWid)
{
   //Comment from TVirtualX:
   // Tells window manager that the window "wid" is a transient window
   // of the window "main_id". A window manager may decide not to decorate
   // a transient window or may treat it differently in other ways.
   //End of TVirtualX's comment.

   //TGTransientFrame uses this hint to attach a window to some "main" window,
   //so that transient window is alway above the main window. This is used for
   //dialogs and dockable panels.
   assert(wid > fPimpl->GetRootWindowID() && "SetWMTransientHint, wid parameter is not a valid window id");

   if (fPimpl->IsRootWindow(mainWid))
      return;

   QuartzWindow * const mainWindow = fPimpl->GetWindow(mainWid).fQuartzWindow;

   if (![mainWindow isVisible])
      return;

   QuartzWindow * const transientWindow = fPimpl->GetWindow(wid).fQuartzWindow;

   if (mainWindow != transientWindow) {
      if (transientWindow.fMainWindow) {
         if (transientWindow.fMainWindow != mainWindow)
            Error("SetWMTransientHint", "window is already transient for other window");
      } else {
         [[transientWindow standardWindowButton : NSWindowZoomButton] setEnabled : NO];
         [mainWindow addTransientWindow : transientWindow];
      }
   } else
      Warning("SetWMTransientHint", "transient and main windows are the same window");
}

#pragma mark - GUI-rendering part.

//______________________________________________________________________________
void TGCocoa::DrawLineAux(Drawable_t wid, const GCValues_t &gcVals, Int_t x1, Int_t y1, Int_t x2, Int_t y2)
{
   //Can be called directly of when flushing command buffer.
   assert(!fPimpl->IsRootWindow(wid) && "DrawLineAux, called for root window");

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);
   CGContextRef ctx = drawable.fContext;
   assert(ctx != 0 && "DrawLineAux, context is null");

   const Quartz::CGStateGuard ctxGuard(ctx);//Will restore state back.
   //Draw a line.
   //This draw line is a special GUI method, it's used not by ROOT's graphics, but
   //widgets. The problem is:
   //-I have to switch off anti-aliasing, since if anti-aliasing is on,
   //the line is thick and has different color.
   //-As soon as I switch-off anti-aliasing, and line is precise, I can not
   //draw a line [0, 0, -> w, 0].
   //I use a small translation, after all, this is ONLY gui method and it
   //will not affect anything except GUI.

   CGContextSetAllowsAntialiasing(ctx, false);//Smoothed line is of wrong color and in a wrong position - this is bad for GUI.

   if (!drawable.fIsPixmap)
      CGContextTranslateCTM(ctx, 0.5, 0.5);
   else {
      //Pixmap uses native Cocoa's left-low-corner system.
      y1 = Int_t(X11::LocalYROOTToCocoa(drawable, y1));
      y2 = Int_t(X11::LocalYROOTToCocoa(drawable, y2));
   }

   SetStrokeParametersFromX11Context(ctx, gcVals);
   CGContextBeginPath(ctx);
   CGContextMoveToPoint(ctx, x1, y1);
   CGContextAddLineToPoint(ctx, x2, y2);
   CGContextStrokePath(ctx);

   CGContextSetAllowsAntialiasing(ctx, true);//Somehow, it's not saved/restored, this affects ... window's titlebar.
}

//______________________________________________________________________________
void TGCocoa::DrawLine(Drawable_t wid, GContext_t gc, Int_t x1, Int_t y1, Int_t x2, Int_t y2)
{
   //This function can be called:
   //a)'normal' way - from view's drawRect method.
   //b) for 'direct rendering' - operation was initiated by ROOT's GUI, not by
   //   drawRect.

   //From TGX11:
   if (!wid)
      return;

   assert(!fPimpl->IsRootWindow(wid) && "DrawLine, called for root window");
   assert(gc > 0 && gc <= fX11Contexts.size() && "DrawLine, invalid context index");

   const GCValues_t &gcVals = fX11Contexts[gc - 1];

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);
   if (!drawable.fIsPixmap) {
      NSObject<X11Window> * const window = (NSObject<X11Window> *)drawable;
      QuartzView *view = (QuartzView *)window.fContentView;
      const ViewFixer fixer(view, wid);
      if (!view.fIsOverlapped && view.fMapState == kIsViewable) {
         if (!view.fContext)
            fPimpl->fX11CommandBuffer.AddDrawLine(wid, gcVals, x1, y1, x2, y2);
         else
            DrawLineAux(wid, gcVals, x1, y1, x2, y2);
      }
   } else {
      if (!IsCocoaDraw()) {
         fPimpl->fX11CommandBuffer.AddDrawLine(wid, gcVals, x1, y1, x2, y2);
      } else {
         DrawLineAux(wid, gcVals, x1, y1, x2, y2);
      }
   }
}

//______________________________________________________________________________
void TGCocoa::DrawSegmentsAux(Drawable_t wid, const GCValues_t &gcVals, const Segment_t *segments, Int_t nSegments)
{
   assert(!fPimpl->IsRootWindow(wid) && "DrawSegmentsAux, called for root window");
   assert(segments != 0 && "DrawSegmentsAux, segments parameter is null");
   assert(nSegments > 0 && "DrawSegmentsAux, nSegments <= 0");

   for (Int_t i = 0; i < nSegments; ++i)
      DrawLineAux(wid, gcVals, segments[i].fX1, segments[i].fY1 - 3, segments[i].fX2, segments[i].fY2 - 3);
}

//______________________________________________________________________________
void TGCocoa::DrawSegments(Drawable_t wid, GContext_t gc, Segment_t *segments, Int_t nSegments)
{
   //Draw multiple line segments. Each line is specified by a pair of points.

   //From TGX11:
   if (!wid)
      return;

   assert(!fPimpl->IsRootWindow(wid) && "DrawSegments, called for root window");
   assert(gc > 0 && gc <= fX11Contexts.size() && "DrawSegments, invalid context index");
   assert(segments != 0 && "DrawSegments, parameter 'segments' is null");
   assert(nSegments > 0 && "DrawSegments, number of segments <= 0");

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);
   const GCValues_t &gcVals = fX11Contexts[gc - 1];

   if (!drawable.fIsPixmap) {
      QuartzView *view = (QuartzView *)fPimpl->GetWindow(wid).fContentView;
      const ViewFixer fixer(view, wid);

      if (!view.fIsOverlapped && view.fMapState == kIsViewable) {
         if (!view.fContext)
            fPimpl->fX11CommandBuffer.AddDrawSegments(wid, gcVals, segments, nSegments);
         else
            DrawSegmentsAux(wid, gcVals, segments, nSegments);
      }
   } else {
      if (!IsCocoaDraw())
         fPimpl->fX11CommandBuffer.AddDrawSegments(wid, gcVals, segments, nSegments);
      else
         DrawSegmentsAux(wid, gcVals, segments, nSegments);
   }
}

//______________________________________________________________________________
void TGCocoa::DrawRectangleAux(Drawable_t wid, const GCValues_t &gcVals, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   //Can be called directly or during flushing command buffer.
   assert(!fPimpl->IsRootWindow(wid) && "DrawRectangleAux, called for root window");

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);

   if (!drawable.fIsPixmap) {
      //I can not draw a line at y == 0, shift the rectangle to 1 pixel (and reduce its height).
      if (!y) {
         y = 1;
         if (h)
            h -= 1;
      }
   } else {
      //Pixmap has native Cocoa's low-left-corner system.
      y = Int_t(X11::LocalYROOTToCocoa(drawable, y + h));
   }

   CGContextRef ctx = fPimpl->GetDrawable(wid).fContext;
   assert(ctx && "DrawRectangleAux, context is null");
   const Quartz::CGStateGuard ctxGuard(ctx);//Will restore context state.

   CGContextSetAllowsAntialiasing(ctx, false);
   //Line color from X11 context.
   SetStrokeParametersFromX11Context(ctx, gcVals);

   const CGRect rect = CGRectMake(x, y, w, h);
   CGContextStrokeRect(ctx, rect);

   CGContextSetAllowsAntialiasing(ctx, true);
}

//______________________________________________________________________________
void TGCocoa::DrawRectangle(Drawable_t wid, GContext_t gc, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   //Can be called in a 'normal way' - from drawRect method (QuartzView)
   //or directly by ROOT.

   if (!wid)//From TGX11.
      return;

   assert(!fPimpl->IsRootWindow(wid) && "DrawRectangle, called for root window");
   assert(gc > 0 && gc <= fX11Contexts.size() && "DrawRectangle, invalid context index");

   const GCValues_t &gcVals = fX11Contexts[gc - 1];

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);

   if (!drawable.fIsPixmap) {
      NSObject<X11Window> * const window = (NSObject<X11Window> *)drawable;
      QuartzView *view = (QuartzView *)window.fContentView;
      const ViewFixer fixer(view, wid);

      if (!view.fIsOverlapped && view.fMapState == kIsViewable) {
         if (!view.fContext)
            fPimpl->fX11CommandBuffer.AddDrawRectangle(wid, gcVals, x, y, w, h);
         else
            DrawRectangleAux(wid, gcVals, x, y, w, h);
      }
   } else {
      if (!IsCocoaDraw())
         fPimpl->fX11CommandBuffer.AddDrawRectangle(wid, gcVals, x, y, w, h);
      else
         DrawRectangleAux(wid, gcVals, x, y, w, h);
   }
}

//______________________________________________________________________________
void TGCocoa::FillRectangleAux(Drawable_t wid, const GCValues_t &gcVals, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   //Can be called directly or when flushing command buffer.
  //Can be called directly or when flushing command buffer.

   //From TGX11:
   if (!wid)
      return;

   assert(!fPimpl->IsRootWindow(wid) && "FillRectangleAux, called for root window");

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);
   CGContextRef ctx = drawable.fContext;
   CGSize patternPhase = {};

   if (drawable.fIsPixmap) {
      //Pixmap has low-left-corner based system.
      y = Int_t(X11::LocalYROOTToCocoa(drawable, y + h));
   }

   const CGRect fillRect = CGRectMake(x, y, w, h);

   if (!drawable.fIsPixmap) {
      QuartzView * const view = (QuartzView *)fPimpl->GetWindow(wid).fContentView;
      if (view.fParentView) {
         const NSPoint origin = [view.fParentView convertPoint : view.frame.origin toView : nil];
         patternPhase.width = origin.x;
         patternPhase.height = origin.y;
      }
   }

   const Quartz::CGStateGuard ctxGuard(ctx);//Will restore context state.

   if (HasFillStippledStyle(gcVals) || HasFillOpaqueStippledStyle(gcVals) ||  HasFillTiledStyle(gcVals)) {
      PatternContext patternContext = {gcVals.fMask, gcVals.fFillStyle, 0, 0, nil, patternPhase};

      if (HasFillStippledStyle(gcVals) || HasFillOpaqueStippledStyle(gcVals)) {
         assert(gcVals.fStipple != kNone &&
                "FillRectangleAux, fill_style is FillStippled/FillOpaqueStippled,"
                " but no stipple is set in a context");

         patternContext.fForeground = gcVals.fForeground;
         patternContext.fImage = fPimpl->GetDrawable(gcVals.fStipple);

         if (HasFillOpaqueStippledStyle(gcVals))
            patternContext.fBackground = gcVals.fBackground;
      } else {
         assert(gcVals.fTile != kNone &&
                "FillRectangleAux, fill_style is FillTiled, but not tile is set in a context");

         patternContext.fImage = fPimpl->GetDrawable(gcVals.fTile);
      }

      SetFillPattern(ctx, &patternContext);
      CGContextFillRect(ctx, fillRect);

      return;
   }

   SetFilledAreaColorFromX11Context(ctx, gcVals);
   CGContextFillRect(ctx, fillRect);
}

//______________________________________________________________________________
void TGCocoa::FillRectangle(Drawable_t wid, GContext_t gc, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   //Can be called in a 'normal way' - from drawRect method (QuartzView)
   //or directly by ROOT.

   //From TGX11:
   if (!wid)
      return;

   assert(!fPimpl->IsRootWindow(wid) && "FillRectangle, called for root window");
   assert(gc > 0 && gc <= fX11Contexts.size() && "FillRectangle, invalid context index");

   const GCValues_t &gcVals = fX11Contexts[gc - 1];
   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);

   if (!drawable.fIsPixmap) {
      NSObject<X11Window> * const window = (NSObject<X11Window> *)drawable;
      QuartzView *view = (QuartzView *)window.fContentView;
      const ViewFixer fixer(view, wid);
      if (!view.fIsOverlapped && view.fMapState == kIsViewable) {
         if (!view.fContext)
            fPimpl->fX11CommandBuffer.AddFillRectangle(wid, gcVals, x, y, w, h);
         else
            FillRectangleAux(wid, gcVals, x, y, w, h);
      }
   } else
      FillRectangleAux(wid, gcVals, x, y, w, h);
}

//______________________________________________________________________________
void TGCocoa::FillPolygonAux(Window_t wid, const GCValues_t &gcVals, const Point_t *polygon, Int_t nPoints)
{
   //Can be called directly or when flushing command buffer.

   //From TGX11:
   if (!wid)
      return;

   assert(!fPimpl->IsRootWindow(wid) && "FillPolygonAux, called for root window");
   assert(polygon != 0 && "FillPolygonAux, parameter 'polygon' is null");
   assert(nPoints > 0 && "FillPolygonAux, number of points must be positive");

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);
   CGContextRef ctx = drawable.fContext;

   CGSize patternPhase = {};

   if (!drawable.fIsPixmap) {
      QuartzView * const view = (QuartzView *)fPimpl->GetWindow(wid).fContentView;
      const NSPoint origin = [view convertPoint : view.frame.origin toView : nil];
      patternPhase.width = origin.x;
      patternPhase.height = origin.y;
   }

   const Quartz::CGStateGuard ctxGuard(ctx);//Will restore context state.

   CGContextSetAllowsAntialiasing(ctx, false);

   if (HasFillStippledStyle(gcVals) || HasFillOpaqueStippledStyle(gcVals) ||  HasFillTiledStyle(gcVals)) {
      PatternContext patternContext = {gcVals.fMask, gcVals.fFillStyle, 0, 0, nil, patternPhase};

      if (HasFillStippledStyle(gcVals) || HasFillOpaqueStippledStyle(gcVals)) {
         assert(gcVals.fStipple != kNone &&
                "FillRectangleAux, fill style is FillStippled/FillOpaqueStippled,"
                " but no stipple is set in a context");

         patternContext.fForeground = gcVals.fForeground;
         patternContext.fImage = fPimpl->GetDrawable(gcVals.fStipple);

         if (HasFillOpaqueStippledStyle(gcVals))
            patternContext.fBackground = gcVals.fBackground;
      } else {
         assert(gcVals.fTile != kNone &&
                "FillRectangleAux, fill_style is FillTiled, but not tile is set in a context");

         patternContext.fImage = fPimpl->GetDrawable(gcVals.fTile);
      }

      SetFillPattern(ctx, &patternContext);
   } else
      SetFilledAreaColorFromX11Context(ctx, gcVals);

   //This +2 -2 shit is the result of ROOT's GUI producing strange coordinates out of ....
   // - first noticed on checkmarks in a menu - they were all shifted.

   CGContextBeginPath(ctx);
   if (!drawable.fIsPixmap) {
      CGContextMoveToPoint(ctx, polygon[0].fX, polygon[0].fY - 2);
      for (Int_t i = 1; i < nPoints; ++i)
         CGContextAddLineToPoint(ctx, polygon[i].fX, polygon[i].fY - 2);
   } else {
      CGContextMoveToPoint(ctx, polygon[0].fX, X11::LocalYROOTToCocoa(drawable, polygon[0].fY + 2));
      for (Int_t i = 1; i < nPoints; ++i)
         CGContextAddLineToPoint(ctx, polygon[i].fX, X11::LocalYROOTToCocoa(drawable, polygon[i].fY + 2));
   }

   CGContextFillPath(ctx);
   CGContextSetAllowsAntialiasing(ctx, true);
}

//______________________________________________________________________________
void TGCocoa::FillPolygon(Window_t wid, GContext_t gc, Point_t *polygon, Int_t nPoints)
{
   // Fills the region closed by the specified path. The path is closed
   // automatically if the last point in the list does not coincide with the
   // first point.
   //
   // Point_t *points - specifies an array of points
   // Int_t npnt      - specifies the number of points in the array
   //
   // GC components in use: function, plane-mask, fill-style, fill-rule,
   // subwindow-mode, clip-x-origin, clip-y-origin, and clip-mask.  GC
   // mode-dependent components: foreground, background, tile, stipple,
   // tile-stipple-x-origin, and tile-stipple-y-origin.
   // (see also the GCValues_t structure)

   //From TGX11:
   if (!wid)
      return;

   assert(polygon != 0 && "FillPolygon, parameter 'polygon' is null");
   assert(nPoints > 0 && "FillPolygon, number of points must be positive");
   assert(gc > 0 && gc <= fX11Contexts.size() && "FillPolygon, invalid context index");

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);
   const GCValues_t &gcVals = fX11Contexts[gc - 1];

   if (!drawable.fIsPixmap) {
      QuartzView *view = (QuartzView *)fPimpl->GetWindow(wid).fContentView;
      const ViewFixer fixer(view, wid);

      if (!view.fIsOverlapped && view.fMapState == kIsViewable) {
         if (!view.fContext)
            fPimpl->fX11CommandBuffer.AddFillPolygon(wid, gcVals, polygon, nPoints);
         else
            FillPolygonAux(wid, gcVals, polygon, nPoints);
      }
   } else {
      if (!IsCocoaDraw())
         fPimpl->fX11CommandBuffer.AddFillPolygon(wid, gcVals, polygon, nPoints);
      else
         FillPolygonAux(wid, gcVals, polygon, nPoints);
   }
}

//______________________________________________________________________________
void TGCocoa::CopyAreaAux(Drawable_t src, Drawable_t dst, const GCValues_t &gcVals, Int_t srcX, Int_t srcY,
                          UInt_t width, UInt_t height, Int_t dstX, Int_t dstY)
{
   //Called directly or when flushing command buffer.
   if (!src || !dst)//Can this happen? From TGX11.
      return;

   assert(!fPimpl->IsRootWindow(src) && "CopyAreaAux, src parameter is root window");
   assert(!fPimpl->IsRootWindow(dst) && "CopyAreaAux, dst parameter is root window");

   //Some copy operations create autoreleased cocoa objects,
   //I do not want them to wait till run loop's iteration end to die.
   const Util::AutoreleasePool pool;

   NSObject<X11Drawable> * const srcDrawable = fPimpl->GetDrawable(src);
   NSObject<X11Drawable> * const dstDrawable = fPimpl->GetDrawable(dst);

   const X11::Point dstPoint(dstX, dstY);
   const X11::Rectangle copyArea(srcX, srcY, width, height);

   QuartzImage *mask = nil;
   if ((gcVals.fMask & kGCClipMask) && gcVals.fClipMask) {
      assert(fPimpl->GetDrawable(gcVals.fClipMask).fIsPixmap == YES &&
             "CopyArea, mask is not a pixmap");
      mask = (QuartzImage *)fPimpl->GetDrawable(gcVals.fClipMask);
   }

   X11::Point clipOrigin;
   if (gcVals.fMask & kGCClipXOrigin)
      clipOrigin.fX = gcVals.fClipXOrigin;
   if (gcVals.fMask & kGCClipYOrigin)
      clipOrigin.fY = gcVals.fClipYOrigin;

   [dstDrawable copy : srcDrawable area : copyArea withMask : mask clipOrigin : clipOrigin toPoint : dstPoint];
}

//______________________________________________________________________________
void TGCocoa::CopyArea(Drawable_t src, Drawable_t dst, GContext_t gc, Int_t srcX, Int_t srcY,
                       UInt_t width, UInt_t height, Int_t dstX, Int_t dstY)
{
   if (!src || !dst)//Can this happen? From TGX11.
      return;

   assert(!fPimpl->IsRootWindow(src) && "CopyArea, src parameter is root window");
   assert(!fPimpl->IsRootWindow(dst) && "CopyArea, dst parameter is root window");
   assert(gc > 0 && gc <= fX11Contexts.size() && "CopyArea, invalid context index");

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(dst);
   const GCValues_t &gcVals = fX11Contexts[gc - 1];

   if (!drawable.fIsPixmap) {
      QuartzView *view = (QuartzView *)fPimpl->GetWindow(dst).fContentView;
      const ViewFixer fixer(view, dst);

      if (!view.fIsOverlapped && view.fMapState == kIsViewable) {
         if (!view.fContext)
            fPimpl->fX11CommandBuffer.AddCopyArea(src, dst, gcVals, srcX, srcY, width, height, dstX, dstY);
         else
            CopyAreaAux(src, dst, gcVals, srcX, srcY, width, height, dstX, dstY);
      }
   } else {
      if (fPimpl->GetDrawable(src).fIsPixmap) {
         //Both are pixmaps, nothing is buffered for src (???).
         CopyAreaAux(src, dst, gcVals, srcX, srcY, width, height, dstX, dstY);
      } else {
         if (!IsCocoaDraw())
            fPimpl->fX11CommandBuffer.AddCopyArea(src, dst, gcVals, srcX, srcY, width, height, dstX, dstY);
         else
            CopyAreaAux(src, dst, gcVals, srcX, srcY, width, height, dstX, dstY);
      }
   }
}

//______________________________________________________________________________
void TGCocoa::DrawStringAux(Drawable_t wid, const GCValues_t &gcVals, Int_t x, Int_t y, const char *text, Int_t len)
{
   //Can be called by ROOT directly, or indirectly by AppKit.
   assert(!fPimpl->IsRootWindow(wid) && "DrawStringAux, called for root window");

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);
   CGContextRef ctx = drawable.fContext;
   assert(ctx != 0 && "DrawStringAux, context is null");

   const Quartz::CGStateGuard ctxGuard(ctx);//Will reset parameters back.

   CGContextSetTextMatrix(ctx, CGAffineTransformIdentity);

   //View is flipped, I have to transform for text to work.
   if (!drawable.fIsPixmap) {
      CGContextTranslateCTM(ctx, 0., drawable.fHeight);
      CGContextScaleCTM(ctx, 1., -1.);
   }

   //Text must be antialiased
   CGContextSetAllowsAntialiasing(ctx, true);

   assert(gcVals.fMask & kGCFont && "DrawString, font is not set in a context");

   if (len < 0)//Negative length can come from caller.
      len = std::strlen(text);
   //Text can be not black, for example, highlighted label.
   CGFloat textColor[4] = {0., 0., 0., 1.};//black by default.
   //I do not check the results here, it's ok to have a black text.
   if (gcVals.fMask & kGCForeground)
      X11::PixelToRGB(gcVals.fForeground, textColor);

   CGContextSetRGBFillColor(ctx, textColor[0], textColor[1], textColor[2], textColor[3]);

   //Do a simple text layout using CGGlyphs.
   //GUI uses non-ascii symbols, and does not care about signed/unsigned - just dump everything
   //into a char and be happy. I'm not.
   std::vector<UniChar> unichars((unsigned char *)text, (unsigned char *)text + len);
   FixAscii(unichars);

   Quartz::DrawTextLineNoKerning(ctx, (CTFontRef)gcVals.fFont, unichars, x,  X11::LocalYROOTToCocoa(drawable, y));
}

//______________________________________________________________________________
void TGCocoa::DrawString(Drawable_t wid, GContext_t gc, Int_t x, Int_t y, const char *text, Int_t len)
{
   //Can be called by ROOT directly, or indirectly by AppKit.
   if (!wid)//from TGX11.
      return;

   assert(!fPimpl->IsRootWindow(wid) && "DrawString, called for root window");
   assert(gc > 0 && gc <= fX11Contexts.size() && "DrawString, invalid context index");

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);
   const GCValues_t &gcVals = fX11Contexts[gc - 1];
   assert(gcVals.fMask & kGCFont && "DrawString, font is not set in a context");

   if (!drawable.fIsPixmap) {
      QuartzView *view = (QuartzView *)fPimpl->GetWindow(wid).fContentView;
      const ViewFixer fixer(view, wid);

      if (!view.fIsOverlapped && view.fMapState == kIsViewable) {
         if (!view.fContext)
            fPimpl->fX11CommandBuffer.AddDrawString(wid, gcVals, x, y, text, len);
         else
            DrawStringAux(wid, gcVals, x, y, text, len);
      }

   } else {
      if (!IsCocoaDraw())
         fPimpl->fX11CommandBuffer.AddDrawString(wid, gcVals, x, y, text, len);
      else
         DrawStringAux(wid, gcVals, x, y, text, len);
   }
}

//______________________________________________________________________________
void TGCocoa::ClearAreaAux(Window_t windowID, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   assert(!fPimpl->IsRootWindow(windowID) && "ClearAreaAux, called for root window");

   QuartzView * const view = (QuartzView *)fPimpl->GetWindow(windowID).fContentView;
   assert(view.fContext != 0 && "ClearAreaAux, view.fContext is null");

   //w and h can be 0 (comment from TGX11) - clear the entire window.
   if (!w)
      w = view.fWidth;
   if (!h)
      h = view.fHeight;

   if (!view.fBackgroundPixmap) {
      //Simple solid fill.
      CGFloat rgb[3] = {};
      X11::PixelToRGB(view.fBackgroundPixel, rgb);

      const Quartz::CGStateGuard ctxGuard(view.fContext);
      CGContextSetRGBFillColor(view.fContext, rgb[0], rgb[1], rgb[2], 1.);//alpha can be also used.
      CGContextFillRect(view.fContext, CGRectMake(x, y, w, h));
   } else {
      const CGRect fillRect = CGRectMake(x, y, w, h);

      CGSize patternPhase = {};
      if (view.fParentView) {
         const NSPoint origin = [view.fParentView convertPoint : view.frame.origin toView : nil];
         patternPhase.width = origin.x;
         patternPhase.height = origin.y;
      }
      const Quartz::CGStateGuard ctxGuard(view.fContext);//Will restore context state.

      PatternContext patternContext = {Mask_t(), 0, 0, 0, view.fBackgroundPixmap, patternPhase};
      SetFillPattern(view.fContext, &patternContext);
      CGContextFillRect(view.fContext, fillRect);
   }
}

//______________________________________________________________________________
void TGCocoa::ClearArea(Window_t wid, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   //Can be called from drawRect method and also by ROOT's GUI directly.
   //Should not be called for pixmap?

   //From TGX11:
   if (!wid)
      return;

   assert(!fPimpl->IsRootWindow(wid) && "ClearArea, called for root window");

   //If wid is pixmap or image, this will crush.
   QuartzView *view = (QuartzView *)fPimpl->GetWindow(wid).fContentView;
   if (ParentRendersToChild(view))
       return;

   if (!view.fIsOverlapped && view.fMapState == kIsViewable) {
      if (!view.fContext)
         fPimpl->fX11CommandBuffer.AddClearArea(wid, x, y, w, h);
      else
         ClearAreaAux(wid, x, y, w, h);
   }
}

//______________________________________________________________________________
void TGCocoa::ClearWindow(Window_t wid)
{
   //Clears the entire area in the specified window (comment from TGX11).

   //From TGX11:
   if (!wid)
      return;

   ClearArea(wid, 0, 0, 0, 0);
}

#pragma mark - Pixmap management.

//______________________________________________________________________________
Int_t TGCocoa::OpenPixmap(UInt_t w, UInt_t h)
{
   //Two stage creation.
   NSSize newSize = {};
   newSize.width = w;
   newSize.height = h;

   Util::NSScopeGuard<QuartzPixmap> pixmap([[QuartzPixmap alloc] initWithW : w H : h
                                           scaleFactor : [[NSScreen mainScreen] backingScaleFactor]]);
   if (pixmap.Get()) {
      pixmap.Get().fID = fPimpl->RegisterDrawable(pixmap.Get());//Can throw.
      return (Int_t)pixmap.Get().fID;
   } else {
      //Detailed error message was issued by QuartzPixmap by this point:
      Error("OpenPixmap", "QuartzPixmap initialization failed");
      return -1;
   }
}

//______________________________________________________________________________
Int_t TGCocoa::ResizePixmap(Int_t wid, UInt_t w, UInt_t h)
{
   assert(!fPimpl->IsRootWindow(wid) && "ResizePixmap, called for root window");

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);
   assert(drawable.fIsPixmap == YES && "ResizePixmap, invalid drawable");

   QuartzPixmap *pixmap = (QuartzPixmap *)drawable;
   if (w == pixmap.fWidth && h == pixmap.fHeight)
      return 1;

   if ([pixmap resizeW : w H : h scaleFactor : [[NSScreen mainScreen] backingScaleFactor]])
      return 1;

   return -1;
}

//______________________________________________________________________________
void TGCocoa::SelectPixmap(Int_t pixmapID)
{
   assert(pixmapID > (Int_t)fPimpl->GetRootWindowID() &&
          "SelectPixmap, parameter 'pixmapID' is not a valid id");

   fSelectedDrawable = pixmapID;
}

//______________________________________________________________________________
void TGCocoa::CopyPixmap(Int_t pixmapID, Int_t x, Int_t y)
{
   assert(pixmapID > (Int_t)fPimpl->GetRootWindowID() &&
          "CopyPixmap, parameter 'pixmapID' is not a valid id");
   assert(fSelectedDrawable > fPimpl->GetRootWindowID() &&
          "CopyPixmap, fSelectedDrawable is not a valid window id");

   NSObject<X11Drawable> * const source = fPimpl->GetDrawable(pixmapID);
   assert([source isKindOfClass : [QuartzPixmap class]] &&
          "CopyPixmap, source is not a pixmap");
   QuartzPixmap * const pixmap = (QuartzPixmap *)source;

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(fSelectedDrawable);
   NSObject<X11Drawable> * destination = nil;

   if (drawable.fIsPixmap) {
      destination = drawable;
   } else {
      NSObject<X11Window> * const window = fPimpl->GetWindow(fSelectedDrawable);
      if (window.fBackBuffer) {
         destination = window.fBackBuffer;
      } else {
         Warning("CopyPixmap", "Operation skipped, since destination"
                               " window is not double buffered");
         return;
      }
   }

   const X11::Rectangle copyArea(0, 0, pixmap.fWidth, pixmap.fHeight);
   const X11::Point dstPoint(x, y);

   [destination copy : pixmap area : copyArea withMask : nil clipOrigin : X11::Point() toPoint : dstPoint];
}

//______________________________________________________________________________
void TGCocoa::ClosePixmap()
{
   // Deletes current pixmap.
   assert(fSelectedDrawable > fPimpl->GetRootWindowID() && "ClosePixmap, no drawable selected");
   assert(fPimpl->GetDrawable(fSelectedDrawable).fIsPixmap == YES && "ClosePixmap, selected drawable is not a pixmap");

   DeletePixmap(fSelectedDrawable);
   fSelectedDrawable = 0;
}

#pragma mark - Different functions to create pixmap from different data sources. Used by GUI.
#pragma mark - These functions implement TVirtualX interface, some of them dupilcate others.

//______________________________________________________________________________
Pixmap_t TGCocoa::CreatePixmap(Drawable_t /*wid*/, UInt_t w, UInt_t h)
{
   //
   return OpenPixmap(w, h);
}

//______________________________________________________________________________
Pixmap_t TGCocoa::CreatePixmap(Drawable_t /*wid*/, const char *bitmap, UInt_t width, UInt_t height,
                               ULong_t foregroundPixel, ULong_t backgroundPixel, Int_t depth)
{
   //Create QuartzImage, using bitmap and foregroundPixel/backgroundPixel,
   //if depth is one - create an image mask instead.

   assert(bitmap != 0 && "CreatePixmap, parameter 'bitmap' is null");
   assert(width > 0 && "CreatePixmap, parameter 'width' is 0");
   assert(height > 0 && "CreatePixmap, parameter 'height' is 0");

   std::vector<unsigned char> imageData (depth > 1 ? width * height * 4 : width * height);

   X11::FillPixmapBuffer((unsigned char*)bitmap, width, height, foregroundPixel,
                          backgroundPixel, depth, &imageData[0]);

   //Now we can create CGImageRef.
   Util::NSScopeGuard<QuartzImage> image;

   if (depth > 1)
      image.Reset([[QuartzImage alloc] initWithW : width H : height data: &imageData[0]]);
   else
      image.Reset([[QuartzImage alloc] initMaskWithW : width H : height bitmapMask : &imageData[0]]);

   if (!image.Get()) {
      Error("CreatePixmap", "QuartzImage initialization failed");//More concrete message was issued by QuartzImage.
      return kNone;
   }

   image.Get().fID = fPimpl->RegisterDrawable(image.Get());//This can throw.
   return image.Get().fID;
}

//______________________________________________________________________________
Pixmap_t TGCocoa::CreatePixmapFromData(unsigned char *bits, UInt_t width, UInt_t height)
{
   //Create QuartzImage, using "bits" (data in bgra format).
   assert(bits != 0 && "CreatePixmapFromData, data parameter is null");
   assert(width != 0 && "CreatePixmapFromData, width parameter is 0");
   assert(height != 0 && "CreatePixmapFromData, height parameter is 0");

   //I'm not using vector here, since I have to pass this pointer to Obj-C code
   //(and Obj-C object will own this memory later).
   std::vector<unsigned char> imageData(bits, bits + width * height * 4);

   //Convert bgra to rgba.
   unsigned char *p = &imageData[0];
   for (unsigned i = 0, e = width * height; i < e; ++i, p += 4)
      std::swap(p[0], p[2]);

   //Now we can create CGImageRef.
   Util::NSScopeGuard<QuartzImage> image([[QuartzImage alloc] initWithW : width
                                          H : height data : &imageData[0]]);

   if (!image.Get()) {
      //Detailed error message was issued by QuartzImage.
      Error("CreatePixmapFromData", "QuartzImage initialziation failed");
      return kNone;
   }

   image.Get().fID = fPimpl->RegisterDrawable(image.Get());//This can throw.
   return image.Get().fID;
}

//______________________________________________________________________________
Pixmap_t TGCocoa::CreateBitmap(Drawable_t /*wid*/, const char *bitmap, UInt_t width, UInt_t height)
{
   //Create QuartzImage with image mask.
   assert(std::numeric_limits<unsigned char>::digits == 8 && "CreateBitmap, ASImage requires octets");

   //I'm not using vector here, since I have to pass this pointer to Obj-C code
   //(and Obj-C object will own this memory later).

   //TASImage has a bug, it calculates size in pixels (making a with to multiple-of eight and
   //allocates memory as each bit occupies one byte, and later packs bits into bytes.

   std::vector<unsigned char> imageData(width * height);

   //TASImage assumes 8-bit bytes and packs mask bits.
   for (unsigned i = 0, j = 0, e = width / 8 * height; i < e; ++i) {
      for(unsigned bit = 0; bit < 8; ++bit, ++j) {
         if (bitmap[i] & (1 << bit))
            imageData[j] = 0;//Opaque.
         else
            imageData[j] = 255;//Masked out bit.
      }
   }

   //Now we can create CGImageRef.
   Util::NSScopeGuard<QuartzImage> image([[QuartzImage alloc] initMaskWithW : width
                                         H : height bitmapMask : &imageData[0]]);
   if (!image.Get()) {
      //Detailed error message was issued by QuartzImage.
      Error("CreateBitmap", "QuartzImage initialization failed");
      return kNone;
   }

   image.Get().fID = fPimpl->RegisterDrawable(image.Get());//This can throw.
   return image.Get().fID;
}

//______________________________________________________________________________
void TGCocoa::DeletePixmapAux(Pixmap_t pixmapID)
{
   fPimpl->DeleteDrawable(pixmapID);
}

//______________________________________________________________________________
void TGCocoa::DeletePixmap(Pixmap_t pixmapID)
{
   // Explicitely deletes the pixmap resource "pmap".
   assert(fPimpl->GetDrawable(pixmapID).fIsPixmap == YES && "DeletePixmap, object is not a pixmap");
   fPimpl->fX11CommandBuffer.AddDeletePixmap(pixmapID);
}

//______________________________________________________________________________
Int_t TGCocoa::AddPixmap(ULong_t /*pixind*/, UInt_t /*w*/, UInt_t /*h*/)
{
   // Registers a pixmap created by TGLManager as a ROOT pixmap
   //
   // w, h - the width and height, which define the pixmap size
   return 0;
}

//______________________________________________________________________________
unsigned char *TGCocoa::GetColorBits(Drawable_t wid, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   //Can be also in a window management part, since window is also drawable.
   if (fPimpl->IsRootWindow(wid)) {
      Warning("GetColorBits", "Called for root window");
   } else {
      assert(x >= 0 && "GetColorBits, parameter 'x' is negative");
      assert(y >= 0 && "GetColorBits, parameter 'y' is negative");
      assert(w != 0 && "GetColorBits, parameter 'w' is 0");
      assert(h != 0 && "GetColorBits, parameter 'h' is 0");

      const X11::Rectangle area(x, y, w, h);
      return [fPimpl->GetDrawable(wid) readColorBits : area];//readColorBits can throw std::bad_alloc, no resource will leak.
   }

   return 0;
}

#pragma mark - XImage emulation.

//______________________________________________________________________________
Drawable_t TGCocoa::CreateImage(UInt_t width, UInt_t height)
{
   // Allocates the memory needed for a drawable.
   //
   // width  - the width of the image, in pixels
   // height - the height of the image, in pixels
   return OpenPixmap(width, height);
}

//______________________________________________________________________________
void TGCocoa::GetImageSize(Drawable_t wid, UInt_t &width, UInt_t &height)
{
   // Returns the width and height of the image wid
   assert(wid > fPimpl->GetRootWindowID() && "GetImageSize, parameter 'wid' is invalid");

   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(wid);
   width = drawable.fWidth;
   height = drawable.fHeight;
}

//______________________________________________________________________________
void TGCocoa::PutPixel(Drawable_t imageID, Int_t x, Int_t y, ULong_t pixel)
{
   // Overwrites the pixel in the image with the specified pixel value.
   // The image must contain the x and y coordinates.
   //
   // imageID - specifies the image
   // x, y  - coordinates
   // pixel - the new pixel value

   assert([fPimpl->GetDrawable(imageID) isKindOfClass : [QuartzPixmap class]] &&
          "PutPixel, parameter 'imageID' is a bad pixmap id");
   assert(x >= 0 && "PutPixel, parameter 'x' is negative");
   assert(y >= 0 && "PutPixel, parameter 'y' is negative");

   QuartzPixmap * const pixmap = (QuartzPixmap *)fPimpl->GetDrawable(imageID);

   unsigned char rgb[3] = {};
   X11::PixelToRGB(pixel, rgb);
   [pixmap putPixel : rgb X : x Y : y];
}

//______________________________________________________________________________
void TGCocoa::PutImage(Drawable_t drawableID, GContext_t gc, Drawable_t imageID, Int_t dstX, Int_t dstY,
                       Int_t srcX, Int_t srcY, UInt_t width, UInt_t height)
{
   //TGX11 uses ZPixmap in CreateImage ... so background/foreground
   //in gc can NEVER be used (and the depth is ALWAYS > 1).
   //This means .... I can call CopyArea!

   CopyArea(imageID, drawableID, gc, srcX, srcY, width, height, dstX, dstY);
}

//______________________________________________________________________________
void TGCocoa::DeleteImage(Drawable_t imageID)
{
   // Deallocates the memory associated with the image img
   assert([fPimpl->GetDrawable(imageID) isKindOfClass : [QuartzPixmap class]] &&
          "DeleteImage, imageID parameter is not a valid image id");
   DeletePixmap(imageID);
}

#pragma mark - Mouse related code.

//______________________________________________________________________________
void TGCocoa::GrabButton(Window_t wid, EMouseButton button, UInt_t keyModifiers, UInt_t eventMask,
                         Window_t /*confine*/, Cursor_t /*cursor*/, Bool_t grab)
{
   //Emulate "passive grab" feature of X11 (similar to "implicit grab" in Cocoa
   //and implicit grab on X11, the difference is that "implicit grab" works as
   //if owner_events parameter for XGrabButton was False, but in ROOT
   //owner_events for XGrabButton is _always_ True.
   //Confine will never be used - no such feature on MacOSX and
   //I'm not going to emulate it..
   //This function also does ungrab.

   //From TGWin32:
   if (!wid)
      return;

   assert(!fPimpl->IsRootWindow(wid) && "GrabButton, called for 'root' window");

   NSObject<X11Window> * const widget = fPimpl->GetWindow(wid);

   if (grab) {
      widget.fPassiveGrabOwnerEvents = YES;   //This is how TGX11 works.
      widget.fPassiveGrabButton = button;
      widget.fPassiveGrabEventMask = eventMask;
      widget.fPassiveGrabKeyModifiers = keyModifiers;
      //Set the cursor.
   } else {
      widget.fPassiveGrabOwnerEvents = NO;
      widget.fPassiveGrabButton = -1;//0 is kAnyButton.
      widget.fPassiveGrabEventMask = 0;
      widget.fPassiveGrabKeyModifiers = 0;
   }
}

//______________________________________________________________________________
void TGCocoa::GrabPointer(Window_t wid, UInt_t eventMask, Window_t /*confine*/, Cursor_t /*cursor*/, Bool_t grab, Bool_t ownerEvents)
{
   //Emulate pointer grab from X11.
   //Confine will never be used - no such feature on MacOSX and
   //I'm not going to emulate it..
   //This function also does ungrab.

   if (grab) {
      NSView<X11Window> * const view = fPimpl->GetWindow(wid).fContentView;
      assert(!fPimpl->IsRootWindow(wid) && "GrabPointer, called for 'root' window");
      //set the cursor.
      //set active grab.
      fPimpl->fX11EventTranslator.SetPointerGrab(view, eventMask, ownerEvents);
   } else {
      //unset cursor?
      //cancel grab.
      fPimpl->fX11EventTranslator.CancelPointerGrab();
   }
}

//______________________________________________________________________________
void TGCocoa::ChangeActivePointerGrab(Window_t, UInt_t, Cursor_t)
{
   // Changes the specified dynamic parameters if the pointer is actively
   // grabbed by the client and if the specified time is no earlier than the
   // last-pointer-grab time and no later than the current X server time.
   //Noop.
}

//______________________________________________________________________________
void TGCocoa::SetKeyAutoRepeat(Bool_t /*on*/)
{
   // Turns key auto repeat on (kTRUE) or off (kFALSE).
   //Noop.
}

//______________________________________________________________________________
void TGCocoa::GrabKey(Window_t wid, Int_t keyCode, UInt_t rootKeyModifiers, Bool_t grab)
{
   //Comment from TVirtualX:
   // Establishes a passive grab on the keyboard. In the future, the
   // keyboard is actively grabbed, the last-keyboard-grab time is set
   // to the time at which the key was pressed (as transmitted in the
   // KeyPress event), and the KeyPress event is reported if all of the
   // following conditions are true:
   //    - the keyboard is not grabbed and the specified key (which can
   //      itself be a modifier key) is logically pressed when the
   //      specified modifier keys are logically down, and no other
   //      modifier keys are logically down;
   //    - either the grab window "id" is an ancestor of (or is) the focus
   //      window, or "id" is a descendant of the focus window and contains
   //      the pointer;
   //    - a passive grab on the same key combination does not exist on any
   //      ancestor of grab_window
   //
   // id       - window id
   // keycode  - specifies the KeyCode or AnyKey
   // modifier - specifies the set of keymasks or AnyModifier; the mask is
   //            the bitwise inclusive OR of the valid keymask bits
   // grab     - a switch between grab/ungrab key
   //            grab = kTRUE  grab the key and modifier
   //            grab = kFALSE ungrab the key and modifier
   //End of comment.


   //Key code already must be Cocoa's key code, this is done by GUI classes,
   //they call KeySymToKeyCode.
   assert(!fPimpl->IsRootWindow(wid) && "GrabKey, called for root window");

   NSView<X11Window> * const view = fPimpl->GetWindow(wid).fContentView;
   const NSUInteger cocoaKeyModifiers = X11::GetCocoaKeyModifiersFromROOTKeyModifiers(rootKeyModifiers);

   if (grab)
      [view addPassiveKeyGrab : keyCode modifiers : cocoaKeyModifiers];
   else
      [view removePassiveKeyGrab : keyCode modifiers : cocoaKeyModifiers];
}

//______________________________________________________________________________
Int_t TGCocoa::KeysymToKeycode(UInt_t keySym)
{
   // Converts the "keysym" to the appropriate keycode. For example,
   // keysym is a letter and keycode is the matching keyboard key (which
   // is dependend on the current keyboard mapping). If the specified
   // "keysym" is not defined for any keycode, returns zero.

   return X11::MapKeySymToKeyCode(keySym);
}

//______________________________________________________________________________
Window_t TGCocoa::GetInputFocus()
{
   // Returns the window id of the window having the input focus.

   return fPimpl->fX11EventTranslator.GetInputFocus();
}

//______________________________________________________________________________
void TGCocoa::SetInputFocus(Window_t wid)
{
   // Changes the input focus to specified window "wid".
   assert(!fPimpl->IsRootWindow(wid) && "SetInputFocus, called for root window");

   if (wid == kNone)
      fPimpl->fX11EventTranslator.SetInputFocus(nil);
   else
      fPimpl->fX11EventTranslator.SetInputFocus(fPimpl->GetWindow(wid).fContentView);
}

//______________________________________________________________________________
void TGCocoa::LookupString(Event_t *event, char *buf, Int_t length, UInt_t &keysym)
{
   // Converts the keycode from the event structure to a key symbol (according
   // to the modifiers specified in the event structure and the current
   // keyboard mapping). In "buf" a null terminated ASCII string is returned
   // representing the string that is currently mapped to the key code.
   //
   // event  - specifies the event structure to be used
   // buf    - returns the translated characters
   // buflen - the length of the buffer
   // keysym - returns the "keysym" computed from the event
   //          if this argument is not NULL
   assert(buf != 0 && "LookupString, parameter 'buf' is null");
   assert(length >= 2 && "LookupString, parameter 'length' - not enough memory to return null-terminated ASCII string");

   X11::MapUnicharToKeySym(event->fCode, buf, length, keysym);
}

#pragma mark - Font management.

//______________________________________________________________________________
FontStruct_t TGCocoa::LoadQueryFont(const char *fontName)
{
   //fontName is in XLFD format:
   //-foundry-family- ..... etc., some components can be omitted and replaced by *.
   assert(fontName != 0 && "LoadQueryFont, fontName is null");

   X11::XLFDName xlfd;
   if (ParseXLFDName(fontName, xlfd)) {
      //Make names more flexible: fFamilyName can be empty or '*'.
      if (!xlfd.fFamilyName.length() || xlfd.fFamilyName == "*")
         xlfd.fFamilyName = "Courier";//Up to me, right?
      if (!xlfd.fPixelSize)
         xlfd.fPixelSize = 11;//Again, up to me.
      return fPimpl->fFontManager.LoadFont(xlfd);
   }

   return FontStruct_t();
}

//______________________________________________________________________________
FontH_t TGCocoa::GetFontHandle(FontStruct_t fs)
{
   return (FontH_t)fs;
}

//______________________________________________________________________________
void TGCocoa::DeleteFont(FontStruct_t fs)
{
   fPimpl->fFontManager.UnloadFont(fs);
}

//______________________________________________________________________________
Bool_t TGCocoa::HasTTFonts() const
{
   // Returns True when TrueType fonts are used
   //No, we use Core Text and do not want TTF to calculate metrics.
   return kFALSE;
}

//______________________________________________________________________________
Int_t TGCocoa::TextWidth(FontStruct_t font, const char *s, Int_t len)
{
   // Return lenght of the string "s" in pixels. Size depends on font.
   return fPimpl->fFontManager.GetTextWidth(font, s, len);
}

//______________________________________________________________________________
void TGCocoa::GetFontProperties(FontStruct_t font, Int_t &maxAscent, Int_t &maxDescent)
{
   // Returns the font properties.
   fPimpl->fFontManager.GetFontProperties(font, maxAscent, maxDescent);
}

//______________________________________________________________________________
FontStruct_t TGCocoa::GetFontStruct(FontH_t fh)
{
   // Retrieves the associated font structure of the font specified font
   // handle "fh".
   //
   // Free returned FontStruct_t using FreeFontStruct().

   return (FontStruct_t)fh;
}

//______________________________________________________________________________
void TGCocoa::FreeFontStruct(FontStruct_t /*fs*/)
{
   // Frees the font structure "fs". The font itself will be freed when
   // no other resource references it.
   //Noop.
}

//______________________________________________________________________________
char **TGCocoa::ListFonts(const char *fontName, Int_t maxNames, Int_t &count)
{
   count = 0;

   if (fontName && fontName[0]) {
      X11::XLFDName xlfd;
      if (X11::ParseXLFDName(fontName, xlfd))
         return fPimpl->fFontManager.ListFonts(xlfd, maxNames, count);
   }

   return 0;
}

//______________________________________________________________________________
void TGCocoa::FreeFontNames(char **fontList)
{
   // Frees the specified the array of strings "fontlist".
   if (!fontList)
      return;

   fPimpl->fFontManager.FreeFontNames(fontList);
}

#pragma mark - Color management.

//______________________________________________________________________________
Bool_t TGCocoa::ParseColor(Colormap_t /*cmap*/, const char *colorName, ColorStruct_t &color)
{
   //"Color" passed as colorName, can be one of the names, defined in X11/rgb.txt,
   //or rgb triplet, which looks like: #rgb #rrggbb #rrrgggbbb #rrrrggggbbbb,
   //where r, g, and b - are hex digits.
   return fPimpl->fX11ColorParser.ParseColor(colorName, color);
}

//______________________________________________________________________________
Bool_t TGCocoa::AllocColor(Colormap_t /*cmap*/, ColorStruct_t &color)
{
   const unsigned red = unsigned(double(color.fRed) / 0xFFFF * 0xFF);
   const unsigned green = unsigned(double(color.fGreen) / 0xFFFF * 0xFF);
   const unsigned blue = unsigned(double(color.fBlue) / 0xFFFF * 0xFF);
   color.fPixel = red << 16 | green << 8 | blue;
   return kTRUE;
}

//______________________________________________________________________________
void TGCocoa::QueryColor(Colormap_t /*cmap*/, ColorStruct_t & color)
{
   // Returns the current RGB value for the pixel in the "color" structure
   color.fRed = (color.fPixel >> 16 & 0xFF) * 0xFFFF / 0xFF;
   color.fGreen = (color.fPixel >> 8 & 0xFF) * 0xFFFF / 0xFF;
   color.fBlue = (color.fPixel & 0xFF) * 0xFFFF / 0xFF;
}

//______________________________________________________________________________
void TGCocoa::FreeColor(Colormap_t /*cmap*/, ULong_t /*pixel*/)
{
   // Frees color cell with specified pixel value.
}

//______________________________________________________________________________
ULong_t TGCocoa::GetPixel(Color_t rootColorIndex)
{
   ULong_t pixel = 0;
   if (const TColor * const color = gROOT->GetColor(rootColorIndex)) {
      Float_t red = 0.f, green = 0.f, blue = 0.f;
      color->GetRGB(red, green, blue);
      pixel = unsigned(red * 255) << 16;
      pixel |= unsigned(green * 255) << 8;
      pixel |= unsigned(blue * 255);
   }

   return pixel;
}

//______________________________________________________________________________
void TGCocoa::GetPlanes(Int_t &nPlanes)
{
   //Implemented as NSBitsPerPixelFromDepth([mainScreen depth]);
   nPlanes = GetDepth();
}

//______________________________________________________________________________
void TGCocoa::GetRGB(Int_t /*index*/, Float_t &/*r*/, Float_t &/*g*/, Float_t &/*b*/)
{
   // Returns RGB values for color "index".
}

//______________________________________________________________________________
void TGCocoa::SetRGB(Int_t /*cindex*/, Float_t /*r*/, Float_t /*g*/, Float_t /*b*/)
{
   // Sets color intensities the specified color index "cindex".
   //
   // cindex  - color index
   // r, g, b - the red, green, blue intensities between 0.0 and 1.0
}

//______________________________________________________________________________
Colormap_t TGCocoa::GetColormap() const
{
   return Colormap_t();
}

#pragma mark - Graphical context management.

//______________________________________________________________________________
GContext_t TGCocoa::CreateGC(Drawable_t /*wid*/, GCValues_t *gval)
{
   //Here I have to imitate graphics context that exists in X11.
   fX11Contexts.push_back(*gval);
   return fX11Contexts.size();
}

//______________________________________________________________________________
void TGCocoa::SetForeground(GContext_t gc, ULong_t foreground)
{
   // Sets the foreground color for the specified GC (shortcut for ChangeGC
   // with only foreground mask set).
   //
   // gc         - specifies the GC
   // foreground - the foreground you want to set
   // (see also the GCValues_t structure)

   assert(gc <= fX11Contexts.size() && gc > 0 && "ChangeGC, invalid context id");

   GCValues_t &x11Context = fX11Contexts[gc - 1];
   x11Context.fMask |= kGCForeground;
   x11Context.fForeground = foreground;
}

//______________________________________________________________________________
void TGCocoa::ChangeGC(GContext_t gc, GCValues_t *gval)
{
   //
   assert(gc <= fX11Contexts.size() && gc > 0 && "ChangeGC, invalid context id");
   assert(gval != 0 && "ChangeGC, gval parameter is null");

   GCValues_t &x11Context = fX11Contexts[gc - 1];
   const Mask_t &mask = gval->fMask;
   x11Context.fMask |= mask;

   //Not all of GCValues_t members are used, but
   //all can be copied/set without any problem.

   if (mask & kGCFunction)
      x11Context.fFunction = gval->fFunction;
   if (mask & kGCPlaneMask)
      x11Context.fPlaneMask = gval->fPlaneMask;
   if (mask & kGCForeground)
      x11Context.fForeground = gval->fForeground;
   if (mask & kGCBackground)
      x11Context.fBackground = gval->fBackground;
   if (mask & kGCLineWidth)
      x11Context.fLineWidth = gval->fLineWidth;
   if (mask & kGCLineStyle)
      x11Context.fLineStyle = gval->fLineStyle;
   if (mask & kGCCapStyle)//nobody uses
      x11Context.fCapStyle = gval->fCapStyle;
   if (mask & kGCJoinStyle)//nobody uses
      x11Context.fJoinStyle = gval->fJoinStyle;
   if (mask & kGCFillRule)//nobody uses
      x11Context.fFillRule = gval->fFillRule;
   if (mask & kGCArcMode)//nobody uses
      x11Context.fArcMode = gval->fArcMode;
   if (mask & kGCFillStyle)
      x11Context.fFillStyle = gval->fFillStyle;
   if (mask & kGCTile)
      x11Context.fTile = gval->fTile;
   if (mask & kGCStipple)
      x11Context.fStipple = gval->fStipple;
   if (mask & kGCTileStipXOrigin)
      x11Context.fTsXOrigin = gval->fTsXOrigin;
   if (mask & kGCTileStipYOrigin)
      x11Context.fTsYOrigin = gval->fTsYOrigin;
   if (mask & kGCFont)
      x11Context.fFont = gval->fFont;
   if (mask & kGCSubwindowMode)
      x11Context.fSubwindowMode = gval->fSubwindowMode;
   if (mask & kGCGraphicsExposures)
      x11Context.fGraphicsExposures = gval->fGraphicsExposures;
   if (mask & kGCClipXOrigin)
      x11Context.fClipXOrigin = gval->fClipXOrigin;
   if (mask & kGCClipYOrigin)
      x11Context.fClipYOrigin = gval->fClipYOrigin;
   if (mask & kGCClipMask)
      x11Context.fClipMask = gval->fClipMask;
   if (mask & kGCDashOffset)
      x11Context.fDashOffset = gval->fDashOffset;
   if (mask & kGCDashList) {
      const unsigned nDashes = sizeof x11Context.fDashes / sizeof x11Context.fDashes[0];
      for (unsigned i = 0; i < nDashes; ++i)
         x11Context.fDashes[i] = gval->fDashes[i];
      x11Context.fDashLen = gval->fDashLen;
   }
}

//______________________________________________________________________________
void TGCocoa::CopyGC(GContext_t src, GContext_t dst, Mask_t mask)
{
   assert(src <= fX11Contexts.size() && src > 0 && "CopyGC, bad source context");
   assert(dst <= fX11Contexts.size() && dst > 0 && "CopyGC, bad destination context");

   GCValues_t srcContext = fX11Contexts[src - 1];
   srcContext.fMask = mask;

   ChangeGC(dst, &srcContext);
}

//______________________________________________________________________________
void TGCocoa::GetGCValues(GContext_t gc, GCValues_t &gval)
{
   // Returns the components specified by the mask in "gval" for the
   // specified GC "gc" (see also the GCValues_t structure)
   const GCValues_t &gcVal = fX11Contexts[gc - 1];
   gval = gcVal;
}

//______________________________________________________________________________
void TGCocoa::DeleteGC(GContext_t /*gc*/)
{
   // Deletes the specified GC "gc".
}

#pragma mark - Cursor management.

//______________________________________________________________________________
Cursor_t TGCocoa::CreateCursor(ECursor cursor)
{
   // Creates the specified cursor. (just return cursor from cursor pool).
   // The cursor can be:
   //
   // kBottomLeft, kBottomRight, kTopLeft,  kTopRight,
   // kBottomSide, kLeftSide,    kTopSide,  kRightSide,
   // kMove,       kCross,       kArrowHor, kArrowVer,
   // kHand,       kRotate,      kPointer,  kArrowRight,
   // kCaret,      kWatch

   return Cursor_t(cursor + 1);//HAHAHAHAHA!!! CREATED!!!
}

//______________________________________________________________________________
void TGCocoa::SetCursor(Int_t wid, ECursor cursor)
{
   // The cursor "cursor" will be used when the pointer is in the
   // window "wid".
   assert(!fPimpl->IsRootWindow(wid) && "SetCursor, called for root window");

   NSView<X11Window> * const view = fPimpl->GetWindow(wid).fContentView;
   view.fCurrentCursor = cursor;
}

//______________________________________________________________________________
void TGCocoa::SetCursor(Window_t wid, Cursor_t cursorID)
{
   // Sets the cursor "curid" to be used when the pointer is in the
   // window "wid".
   if (cursorID > 0)
      SetCursor(Int_t(wid), ECursor(cursorID - 1));
   else
      SetCursor(Int_t(wid), kPointer);
}

//______________________________________________________________________________
void TGCocoa::QueryPointer(Int_t &x, Int_t &y)
{
   // Returns the pointer position.

   //I ignore fSelectedDrawable here. If you have any problems with this, hehe, you can ask me :)
   const NSPoint screenPoint = [NSEvent mouseLocation];
   x = X11::GlobalXCocoaToROOT(screenPoint.x);
   y = X11::GlobalYCocoaToROOT(screenPoint.y);
}

//______________________________________________________________________________
void TGCocoa::QueryPointer(Window_t winID, Window_t &rootWinID, Window_t &childWinID,
                           Int_t &rootX, Int_t &rootY, Int_t &winX, Int_t &winY, UInt_t &mask)
{
   //Emulate XQueryPointer.

   //From TGX11/TGWin32:
   if (!winID)
      return;//Neither TGX11, nor TGWin32 set any of out parameters.

   //We have only one root window.
   rootWinID = fPimpl->GetRootWindowID();
   //Find cursor position (screen coordinates).
   NSPoint screenPoint = [NSEvent mouseLocation];
   screenPoint.x = X11::GlobalXCocoaToROOT(screenPoint.x);
   screenPoint.y = X11::GlobalYCocoaToROOT(screenPoint.y);
   rootX = screenPoint.x;
   rootY = screenPoint.y;

   //Convert a screen point to winID's coordinate system.
   if (winID > fPimpl->GetRootWindowID()) {
      NSObject<X11Window> * const window = fPimpl->GetWindow(winID);
      const NSPoint winPoint = X11::TranslateFromScreen(screenPoint, window.fContentView);
      winX = winPoint.x;
      winY = winPoint.y;
   } else {
      winX = screenPoint.x;
      winY = screenPoint.y;
   }

   //Find child window in these coordinates (?).
   if (QuartzWindow * const childWin = X11::FindWindowInPoint(screenPoint.x, screenPoint.y)) {
      childWinID = childWin.fID;
      mask = X11::GetModifiers();
   } else {
      childWinID = 0;
      mask = 0;
   }
}

#pragma mark - OpenGL management.

//______________________________________________________________________________
Double_t TGCocoa::GetOpenGLScalingFactor()
{
   //Scaling factor to let our OpenGL code know, that we probably
   //work on a retina display.

   return [[NSScreen mainScreen] backingScaleFactor];
}

//______________________________________________________________________________
Window_t TGCocoa::CreateOpenGLWindow(Window_t parentID, UInt_t width, UInt_t height,
                                     const std::vector<std::pair<UInt_t, Int_t> > &formatComponents)
{
   //ROOT never creates GL widgets with 'root' as a parent (so not top-level gl-windows).
   //If this change, assert must be deleted.
   typedef std::pair<UInt_t, Int_t> component_type;
   typedef std::vector<component_type>::size_type size_type;

   //Convert pairs into Cocoa's GL attributes.
   std::vector<NSOpenGLPixelFormatAttribute> attribs;
   for (size_type i = 0, e = formatComponents.size(); i < e; ++i) {
      const component_type &comp = formatComponents[i];

      if (comp.first == Rgl::kDoubleBuffer) {
         attribs.push_back(NSOpenGLPFADoubleBuffer);
      } else if (comp.first == Rgl::kDepth) {
         attribs.push_back(NSOpenGLPFADepthSize);
         attribs.push_back(comp.second > 0 ? comp.second : 32);
      } else if (comp.first == Rgl::kAccum) {
         attribs.push_back(NSOpenGLPFAAccumSize);
         attribs.push_back(comp.second > 0 ? comp.second : 1);
      } else if (comp.first == Rgl::kStencil) {
         attribs.push_back(NSOpenGLPFAStencilSize);
         attribs.push_back(comp.second > 0 ? comp.second : 8);
      } else if (comp.first == Rgl::kMultiSample) {
         attribs.push_back(NSOpenGLPFAMultisample);
         attribs.push_back(NSOpenGLPFASampleBuffers);
         attribs.push_back(1);
         attribs.push_back(NSOpenGLPFASamples);
         attribs.push_back(comp.second ? comp.second : 8);
      }
   }

   attribs.push_back(0);

   NSOpenGLPixelFormat * const pixelFormat = [[NSOpenGLPixelFormat alloc] initWithAttributes : &attribs[0]];
   const Util::NSScopeGuard<NSOpenGLPixelFormat> formatGuard(pixelFormat);

   NSView<X11Window> *parentView = nil;
   if (!fPimpl->IsRootWindow(parentID)) {
      parentView = fPimpl->GetWindow(parentID).fContentView;
      assert([parentView isKindOfClass : [QuartzView class]] &&
             "CreateOpenGLWindow, parent view must be QuartzView");
   }

   NSRect viewFrame = {};
   viewFrame.size.width = width;
   viewFrame.size.height = height;

   ROOTOpenGLView * const glView = [[ROOTOpenGLView alloc] initWithFrame : viewFrame pixelFormat : pixelFormat];
   const Util::NSScopeGuard<ROOTOpenGLView> viewGuard(glView);

   Window_t glID = kNone;

   if (parentView) {
      [parentView addChild : glView];
      glID = fPimpl->RegisterDrawable(glView);
      glView.fID = glID;
   } else {
      //"top-level glview".
      //Create a window to be parent of this gl-view.
      QuartzWindow *parent = [[QuartzWindow alloc] initWithGLView : glView];
      const Util::NSScopeGuard<QuartzWindow> winGuard(parent);


      if (!parent) {
         Error("CreateOpenGLWindow", "QuartzWindow allocation/initialization"
                                     " failed for a top-level GL widget");
         return kNone;
      }

      glID = fPimpl->RegisterDrawable(parent);
      parent.fID = glID;
   }

   return glID;
}

//______________________________________________________________________________
Handle_t TGCocoa::CreateOpenGLContext(Window_t windowID, Handle_t sharedID)
{
   assert(!fPimpl->IsRootWindow(windowID) &&
          "CreateOpenGLContext, parameter 'windowID' is a root window");
   assert([fPimpl->GetWindow(windowID).fContentView isKindOfClass : [ROOTOpenGLView class]] &&
          "CreateOpenGLContext, view is not an OpenGL view");

   NSOpenGLContext * const sharedContext = fPimpl->GetGLContextForHandle(sharedID);
   ROOTOpenGLView * const glView = (ROOTOpenGLView *)fPimpl->GetWindow(windowID);

   const Util::NSScopeGuard<NSOpenGLContext>
      newContext([[NSOpenGLContext alloc] initWithFormat : glView.pixelFormat shareContext : sharedContext]);
   glView.fOpenGLContext = newContext.Get();
   const Handle_t ctxID = fPimpl->RegisterGLContext(newContext.Get());

   return ctxID;
}

//______________________________________________________________________________
void TGCocoa::CreateOpenGLContext(Int_t /*wid*/)
{
   // Creates OpenGL context for window "wid"
}

//______________________________________________________________________________
Bool_t TGCocoa::MakeOpenGLContextCurrent(Handle_t ctxID, Window_t windowID)
{
   using namespace Details;

   assert(ctxID > 0 && "MakeOpenGLContextCurrent, invalid context id");

   NSOpenGLContext * const glContext = fPimpl->GetGLContextForHandle(ctxID);
   if (!glContext) {
      Error("MakeOpenGLContextCurrent", "No OpenGL context found for id %d", int(ctxID));

      return kFALSE;
   }

   ROOTOpenGLView * const glView = (ROOTOpenGLView *)fPimpl->GetWindow(windowID).fContentView;

   if (OpenGL::GLViewIsValidDrawable(glView)) {
      if ([glContext view] != glView)
         [glContext setView : glView];

      if (glView.fUpdateContext) {
         [glContext update];
         glView.fUpdateContext = NO;
      }

      glView.fOpenGLContext = glContext;
      [glContext makeCurrentContext];

      return kTRUE;
   } else {
      //Oh, here's the real black magic.
      //Our brilliant GL code is sure that MakeCurrent always succeeds.
      //But it does not: if view is not visible, context can not be attached,
      //gl operations will fail.
      //Funny enough, but if you have invisible window with visible view,
      //this trick works.

      NSView *fakeView = nil;
      QuartzWindow *fakeWindow = fPimpl->GetFakeGLWindow();

      if (!fakeWindow) {
         //We did not find any window. Create a new one.
         SetWindowAttributes_t attr = {};
         //100 - is just a stupid hardcoded value:
         const UInt_t width = std::max(glView.frame.size.width, CGFloat(100));
         const UInt_t height = std::max(glView.frame.size.height, CGFloat(100));

         NSRect viewFrame = {};
         viewFrame.size.width = width;
         viewFrame.size.height = height;

         const NSUInteger styleMask = kTitledWindowMask | kClosableWindowMask |
                                      kMiniaturizableWindowMask | kResizableWindowMask;

         //NOTE: defer parameter is 'NO', otherwise this trick will not help.
         fakeWindow = [[QuartzWindow alloc] initWithContentRect : viewFrame styleMask : styleMask
                        backing : NSBackingStoreBuffered defer : NO windowAttributes : &attr];
         Util::NSScopeGuard<QuartzWindow> winGuard(fakeWindow);

         fakeView = fakeWindow.fContentView;
         [fakeView setHidden : NO];//!

         fPimpl->SetFakeGLWindow(fakeWindow);//Can throw.
         winGuard.Release();
      } else {
         fakeView = fakeWindow.fContentView;
         [fakeView setHidden : NO];
      }

      glView.fOpenGLContext = nil;
      [glContext setView : fakeView];
      [glContext makeCurrentContext];
   }

   return kTRUE;
}

//______________________________________________________________________________
Handle_t TGCocoa::GetCurrentOpenGLContext()
{
   NSOpenGLContext * const currentContext = [NSOpenGLContext currentContext];
   if (!currentContext) {
      Error("GetCurrentOpenGLContext", "The current OpenGL context is null");
      return kNone;
   }

   const Handle_t contextID = fPimpl->GetHandleForGLContext(currentContext);
   if (!contextID)
      Error("GetCurrentOpenGLContext", "The current OpenGL context was"
                                       " not created/registered by TGCocoa");

   return contextID;
}

//______________________________________________________________________________
void TGCocoa::FlushOpenGLBuffer(Handle_t ctxID)
{
   assert(ctxID > 0 && "FlushOpenGLBuffer, invalid context id");

   NSOpenGLContext * const glContext = fPimpl->GetGLContextForHandle(ctxID);
   assert(glContext != nil && "FlushOpenGLBuffer, bad context id");

   if (glContext != [NSOpenGLContext currentContext])//???
      return;

   glFlush();//???
   [glContext flushBuffer];
}

//______________________________________________________________________________
void TGCocoa::DeleteOpenGLContext(Int_t ctxID)
{
   //Historically, DeleteOpenGLContext was accepting window id,
   //now it's a context id. DeleteOpenGLContext is not used in ROOT,
   //only in TGLContext for Cocoa.
   NSOpenGLContext * const glContext = fPimpl->GetGLContextForHandle(ctxID);
   if (NSView * const v = [glContext view]) {
      if ([v isKindOfClass : [ROOTOpenGLView class]])
         ((ROOTOpenGLView *)v).fOpenGLContext = nil;

      [glContext clearDrawable];
   }

   if (glContext == [NSOpenGLContext currentContext])
      [NSOpenGLContext clearCurrentContext];

   fPimpl->DeleteGLContext(ctxID);
}

#pragma mark - Off-screen rendering for TPad/TCanvas.

//______________________________________________________________________________
void TGCocoa::SetDoubleBuffer(Int_t windowID, Int_t mode)
{
   //In ROOT, canvas has a "double buffer" - pixmap attached to 'wid'.
   assert(windowID > (Int_t)fPimpl->GetRootWindowID() && "SetDoubleBuffer called for root window");

   if (windowID == 999) {//Comment in TVirtaulX suggests, that 999 means all windows.
      Warning("SetDoubleBuffer", "called with wid == 999");
      //Window with id 999 can not exists - this is checked in CocoaPrivate.
   } else {
      fSelectedDrawable = windowID;
      mode ? SetDoubleBufferON() : SetDoubleBufferOFF();
   }
}

//______________________________________________________________________________
void TGCocoa::SetDoubleBufferOFF()
{
   fDirectDraw = true;
}

//______________________________________________________________________________
void TGCocoa::SetDoubleBufferON()
{
   //Attach pixmap to the selected window (view).
   fDirectDraw = false;

   assert(fSelectedDrawable > fPimpl->GetRootWindowID() &&
          "SetDoubleBufferON, called, but no correct window was selected before");

   NSObject<X11Window> * const window = fPimpl->GetWindow(fSelectedDrawable);

   if (!window) return;

   assert(window.fIsPixmap == NO &&
          "SetDoubleBufferON, selected drawable is a pixmap, can not attach pixmap to pixmap");

   const unsigned currW = window.fWidth;
   const unsigned currH = window.fHeight;

   if (QuartzPixmap *const currentPixmap = window.fBackBuffer) {
      if (currH == currentPixmap.fHeight && currW == currentPixmap.fWidth)
         return;
   }

   Util::NSScopeGuard<QuartzPixmap> pixmap([[QuartzPixmap alloc] initWithW : currW
                                            H : currH scaleFactor : [[NSScreen mainScreen] backingScaleFactor]]);
   if (pixmap.Get())
      window.fBackBuffer = pixmap.Get();
   else
      //Detailed error message was issued by QuartzPixmap.
      Error("SetDoubleBufferON", "QuartzPixmap initialization failed");
}

//______________________________________________________________________________
void TGCocoa::SetDrawMode(EDrawMode mode)
{
   // Sets the drawing mode.
   //
   //EDrawMode{kCopy, kXor, kInvert};
   if (fDrawMode == kInvert && mode != kInvert) {
       // Remove previously added CrosshairWindow.
       auto windows = NSApplication.sharedApplication.windows;
       for (NSWindow *candidate : windows) {
           if ([candidate isKindOfClass:QuartzWindow.class])
               [(QuartzWindow *)candidate removeCrosshairWindow];
       }
       fPimpl->fX11CommandBuffer.ClearXOROperations();
   }

   fDrawMode = mode;
}

#pragma mark - Event management part.

//______________________________________________________________________________
void TGCocoa::SendEvent(Window_t wid, Event_t *event)
{
   if (fPimpl->IsRootWindow(wid))//ROOT's GUI can send events to root window.
      return;

   //From TGX11:
   if (!wid || !event)
      return;

   Event_t newEvent = *event;
   newEvent.fWindow = wid;
   fPimpl->fX11EventTranslator.fEventQueue.push_back(newEvent);
}

//______________________________________________________________________________
void TGCocoa::NextEvent(Event_t &event)
{
   assert(fPimpl->fX11EventTranslator.fEventQueue.size() > 0 && "NextEvent, event queue is empty");

   event = fPimpl->fX11EventTranslator.fEventQueue.front();
   fPimpl->fX11EventTranslator.fEventQueue.pop_front();
}

//______________________________________________________________________________
Int_t TGCocoa::EventsPending()
{
   return (Int_t)fPimpl->fX11EventTranslator.fEventQueue.size();
}


//______________________________________________________________________________
Bool_t TGCocoa::CheckEvent(Window_t windowID, EGEventType type, Event_t &event)
{
   typedef X11::EventQueue_t::iterator iterator_type;

   iterator_type it = fPimpl->fX11EventTranslator.fEventQueue.begin();
   iterator_type eIt = fPimpl->fX11EventTranslator.fEventQueue.end();

   for (; it != eIt; ++it) {
      const Event_t &queuedEvent = *it;
      if (queuedEvent.fWindow == windowID && queuedEvent.fType == type) {
         event = queuedEvent;
         fPimpl->fX11EventTranslator.fEventQueue.erase(it);
         return kTRUE;
      }
   }

   return kFALSE;
}

//______________________________________________________________________________
Handle_t TGCocoa::GetNativeEvent() const
{
   //I can not give an access to the native event,
   //it even, probably, does not exist already.
   return kNone;
}

#pragma mark - "Drag and drop", "Copy and paste", X11 properties.

//______________________________________________________________________________
Atom_t  TGCocoa::InternAtom(const char *name, Bool_t onlyIfExist)
{
   //X11 properties emulation.

   assert(name != 0 && "InternAtom, parameter 'name' is null");
   return FindAtom(name, !onlyIfExist);
}

//______________________________________________________________________________
void TGCocoa::SetPrimarySelectionOwner(Window_t windowID)
{
   //Comment from TVirtualX:
   // Makes the window "wid" the current owner of the primary selection.
   // That is the window in which, for example some text is selected.
   //End of comment.

   //It's not clear, why SetPrimarySelectionOwner and SetSelectionOwner have different return types.

   if (!windowID)//From TGWin32.
      return;

   assert(!fPimpl->IsRootWindow(windowID) &&
          "SetPrimarySelectionOwner, windowID parameter is a 'root' window");
   assert(fPimpl->GetDrawable(windowID).fIsPixmap == NO &&
          "SetPrimarySelectionOwner, windowID parameter is not a valid window");

   const Atom_t primarySelectionAtom = FindAtom("XA_PRIMARY", false);
   assert(primarySelectionAtom != kNone &&
          "SetPrimarySelectionOwner, predefined XA_PRIMARY atom was not found");

   fSelectionOwners[primarySelectionAtom] = windowID;
   //No events will be send - I do not have different clients, so nobody to send SelectionClear.
}

//______________________________________________________________________________
Bool_t TGCocoa::SetSelectionOwner(Window_t windowID, Atom_t &selection)
{
   //Comment from TVirtualX:
   // Changes the owner and last-change time for the specified selection.
   //End of comment.

   //It's not clear, why SetPrimarySelectionOwner and SetSelectionOwner have different return types.

   if (!windowID)
      return kFALSE;

   assert(!fPimpl->IsRootWindow(windowID) &&
          "SetSelectionOwner, windowID parameter is a 'root' window'");
   assert(fPimpl->GetDrawable(windowID).fIsPixmap == NO &&
          "SetSelectionOwner, windowID parameter is not a valid window");

   fSelectionOwners[selection] = windowID;
   //No messages, since I do not have different clients.

   return kTRUE;
}

//______________________________________________________________________________
Window_t TGCocoa::GetPrimarySelectionOwner()
{
   //Comment from TVirtualX:
   // Returns the window id of the current owner of the primary selection.
   // That is the window in which, for example some text is selected.
   //End of comment.
   const Atom_t primarySelectionAtom = FindAtom("XA_PRIMARY", false);
   assert(primarySelectionAtom != kNone &&
          "GetPrimarySelectionOwner, predefined XA_PRIMARY atom was not found");

   return fSelectionOwners[primarySelectionAtom];
}

//______________________________________________________________________________
void TGCocoa::ConvertPrimarySelection(Window_t windowID, Atom_t clipboard, Time_t when)
{
   //Comment from TVirtualX:
   // Causes a SelectionRequest event to be sent to the current primary
   // selection owner. This event specifies the selection property
   // (primary selection), the format into which to convert that data before
   // storing it (target = XA_STRING), the property in which the owner will
   // place the information (sel_property), the window that wants the
   // information (id), and the time of the conversion request (when).
   // The selection owner responds by sending a SelectionNotify event, which
   // confirms the selected atom and type.
   //End of comment.

   //From TGWin32:
   if (!windowID)
      return;

   assert(!fPimpl->IsRootWindow(windowID) &&
          "ConvertPrimarySelection, parameter 'windowID' is root window");
   assert(fPimpl->GetDrawable(windowID).fIsPixmap == NO &&
          "ConvertPrimarySelection, parameter windowID parameter is not a window id");

   Atom_t primarySelectionAtom = FindAtom("XA_PRIMARY", false);
   assert(primarySelectionAtom != kNone &&
          "ConvertPrimarySelection, XA_PRIMARY predefined atom not found");

   Atom_t stringAtom = FindAtom("XA_STRING", false);
   assert(stringAtom != kNone &&
          "ConvertPrimarySelection, XA_STRING predefined atom not found");

   ConvertSelection(windowID, primarySelectionAtom, stringAtom, clipboard, when);
}

//______________________________________________________________________________
void TGCocoa::ConvertSelection(Window_t windowID, Atom_t &selection, Atom_t &target,
                               Atom_t &property, Time_t &/*timeStamp*/)
{
   // Requests that the specified selection be converted to the specified
   // target type.

   // Requests that the specified selection be converted to the specified
   // target type.

   if (!windowID)
      return;

   assert(!fPimpl->IsRootWindow(windowID) &&
          "ConvertSelection, parameter 'windowID' is root window'");
   assert(fPimpl->GetDrawable(windowID).fIsPixmap == NO &&
          "ConvertSelection, parameter 'windowID' is not a window id");

   Event_t newEvent = {};
   selection_iterator selIter = fSelectionOwners.find(selection);

   if (selIter != fSelectionOwners.end())
      newEvent.fType = kSelectionRequest;
   else
      newEvent.fType = kSelectionNotify;

   newEvent.fWindow = windowID;
   newEvent.fUser[0] = windowID;//requestor
   newEvent.fUser[1] = selection;
   newEvent.fUser[2] = target;
   newEvent.fUser[3] = property;

   SendEvent(windowID, &newEvent);
}

//______________________________________________________________________________
Int_t TGCocoa::GetProperty(Window_t windowID, Atom_t propertyID, Long_t, Long_t, Bool_t, Atom_t,
                           Atom_t *actualType, Int_t *actualFormat, ULong_t *nItems,
                           ULong_t *bytesAfterReturn, unsigned char **propertyReturn)
{
   //Comment from TVirtualX:
   // Returns the actual type of the property; the actual format of the property;
   // the number of 8-bit, 16-bit, or 32-bit items transferred; the number of
   // bytes remaining to be read in the property; and a pointer to the data
   // actually returned.
   //End of comment.

   if (fPimpl->IsRootWindow(windowID))
      return 0;

   assert(fPimpl->GetDrawable(windowID).fIsPixmap == NO &&
          "GetProperty, parameter 'windowID' is not a valid window id");
   assert(propertyID > 0 && propertyID <= fAtomToName.size() &&
          "GetProperty, parameter 'propertyID' is not a valid atom");
   assert(actualType != 0 && "GetProperty, parameter 'actualType' is null");
   assert(actualFormat != 0 && "GetProperty, parameter 'actualFormat' is null");
   assert(bytesAfterReturn != 0 && "GetProperty, parameter 'bytesAfterReturn' is null");
   assert(propertyReturn != 0 && "GetProperty, parameter 'propertyReturn' is null");

   const Util::AutoreleasePool pool;

   *bytesAfterReturn = 0;//In TGWin32 the value set to .. nItems?
   *propertyReturn = 0;
   *nItems = 0;

   const std::string &atomName = fAtomToName[propertyID - 1];
   NSObject<X11Window> *window = fPimpl->GetWindow(windowID);

   if (![window hasProperty : atomName.c_str()]) {
      Error("GetProperty", "Unknown property %s requested", atomName.c_str());
      return 0;//actually, 0 is ... Success (X11)?
   }

   unsigned tmpFormat = 0, tmpElements = 0;
   *propertyReturn = [window getProperty : atomName.c_str() returnType : actualType
                      returnFormat : &tmpFormat nElements : &tmpElements];
   *actualFormat = (Int_t)tmpFormat;
   *nItems = tmpElements;

   return *nItems;//Success (X11) is 0?
}

//______________________________________________________________________________
void TGCocoa::GetPasteBuffer(Window_t windowID, Atom_t propertyID, TString &text,
                             Int_t &nChars, Bool_t clearBuffer)
{
   //Comment from TVirtualX:
   // Gets contents of the paste buffer "atom" into the string "text".
   // (nchar = number of characters) If "del" is true deletes the paste
   // buffer afterwards.
   //End of comment.

   //From TGX11:
   if (!windowID)
      return;

   assert(!fPimpl->IsRootWindow(windowID) &&
          "GetPasteBuffer, parameter 'windowID' is root window");
   assert(fPimpl->GetDrawable(windowID).fIsPixmap == NO &&
          "GetPasteBuffer, parameter 'windowID' is not a valid window");
   assert(propertyID && propertyID <= fAtomToName.size() &&
          "GetPasteBuffer, parameter 'propertyID' is not a valid atom");

   const Util::AutoreleasePool pool;

   const std::string &atomString = fAtomToName[propertyID - 1];
   NSObject<X11Window> *window = fPimpl->GetWindow(windowID);

   if (![window hasProperty : atomString.c_str()]) {
      Error("GetPasteBuffer", "No property %s on a window", atomString.c_str());
      return;
   }

   Atom_t tmpType = 0;
   unsigned tmpFormat = 0, nElements = 0;

   const Util::ScopedArray<char>
      propertyData((char *)[window getProperty : atomString.c_str()
                            returnType : &tmpType returnFormat : &tmpFormat
                            nElements : &nElements]);

   assert(tmpFormat == 8 && "GetPasteBuffer, property has wrong format");

   text.Insert(0, propertyData.Get(), nElements);
   nChars = (Int_t)nElements;

   if (clearBuffer) {
      //For the moment - just remove the property
      //(anyway, ChangeProperty/ChangeProperties will re-create it).
      [window removeProperty : atomString.c_str()];
   }
}

//______________________________________________________________________________
void TGCocoa::ChangeProperty(Window_t windowID, Atom_t propertyID, Atom_t type,
                             UChar_t *data, Int_t len)
{
   //Comment from TVirtualX:
   // Alters the property for the specified window and causes the X server
   // to generate a PropertyNotify event on that window.
   //
   // wid       - the window whose property you want to change
   // property - specifies the property name
   // type     - the type of the property; the X server does not
   //            interpret the type but simply passes it back to
   //            an application that might ask about the window
   //            properties
   // data     - the property data
   // len      - the length of the specified data format
   //End of comment.

   //TGX11 always calls XChangeProperty with PropModeReplace.
   //I simply reset the property (or create a new one).

   if (!windowID) //From TGWin32.
      return;

   if (!data || !len) //From TGWin32.
      return;

   assert(!fPimpl->IsRootWindow(windowID) &&
          "ChangeProperty, parameter 'windowID' is root window");
   assert(fPimpl->GetDrawable(windowID).fIsPixmap == NO &&
          "ChangeProperty, parameter 'windowID' is not a valid window id");
   assert(propertyID && propertyID <= fAtomToName.size() &&
          "ChangeProperty, parameter 'propertyID' is not a valid atom");

   const Util::AutoreleasePool pool;

   const std::string &atomString = fAtomToName[propertyID - 1];

   NSObject<X11Window> * const window = fPimpl->GetWindow(windowID);
   [window setProperty : atomString.c_str() data : data size : len forType : type format : 8];
   //ROOT ignores PropertyNotify events.
}

//______________________________________________________________________________
void TGCocoa::ChangeProperties(Window_t windowID, Atom_t propertyID, Atom_t type,
                               Int_t format, UChar_t *data, Int_t len)
{
   //Comment from TVirtualX:
   // Alters the property for the specified window and causes the X server
   // to generate a PropertyNotify event on that window.
   //End of comment.

   //TGX11 always calls XChangeProperty with PropModeReplace.
   //I simply reset the property (or create a new one).

   if (!windowID)//From TGWin32.
      return;

   if (!data || !len)//From TGWin32.
      return;

   assert(!fPimpl->IsRootWindow(windowID) &&
          "ChangeProperties, parameter 'windowID' is root window");
   assert(fPimpl->GetDrawable(windowID).fIsPixmap == NO &&
          "ChangeProperties, parameter 'windowID' is not a valid window id");
   assert(propertyID && propertyID <= fAtomToName.size() &&
          "ChangeProperties, parameter 'propertyID' is not a valid atom");

   const Util::AutoreleasePool pool;

   const std::string &atomName = fAtomToName[propertyID - 1];

   NSObject<X11Window> * const window = fPimpl->GetWindow(windowID);
   [window setProperty : atomName.c_str() data : data
                  size : len forType : type format : format];
   //No property notify, ROOT does not know about this.
}

//______________________________________________________________________________
void TGCocoa::DeleteProperty(Window_t windowID, Atom_t &propertyID)
{
   //Comment from TVirtualX:
   // Deletes the specified property only if the property was defined on the
   // specified window and causes the X server to generate a PropertyNotify
   // event on the window unless the property does not exist.
   //End of comment.

   if (!windowID)//Can this happen?
      return;

   //Strange signature - why propertyID is a reference?
   assert(!fPimpl->IsRootWindow(windowID) &&
          "DeleteProperty, parameter 'windowID' is root window");
   assert(fPimpl->GetDrawable(windowID).fIsPixmap == NO &&
          "DeleteProperty, parameter 'windowID' is not a valid window");
   assert(propertyID && propertyID <= fAtomToName.size() &&
          "DeleteProperty, parameter 'propertyID' is not a valid atom");

   const std::string &atomString = fAtomToName[propertyID - 1];
   [fPimpl->GetWindow(windowID) removeProperty : atomString.c_str()];
}

//______________________________________________________________________________
void TGCocoa::SetDNDAware(Window_t windowID, Atom_t *typeList)
{
   //Comment from TVirtaulX:
   // Add XdndAware property and the list of drag and drop types to the
   // Window win.
   //End of comment.


   //TGX11 first replaces XdndAware property for a windowID, and then appends atoms from a typelist.
   //I simply put all data for a property into a vector and set the property (either creating
   //a new property or replacing the existing).

   assert(windowID > fPimpl->GetRootWindowID() &&
          "SetDNDAware, parameter 'windowID' is not a valid window id");
   assert(fPimpl->GetDrawable(windowID).fIsPixmap == NO &&
          "SetDNDAware, parameter 'windowID' is not a window");

   const Util::AutoreleasePool pool;

   QuartzView * const view = (QuartzView *)fPimpl->GetWindow(windowID).fContentView;
   NSArray * const supportedTypes = [NSArray arrayWithObjects : NSFilenamesPboardType, nil];//In a pool.

   //Do this for Cocoa - to make it possible to drag something to a
   //ROOT's window (also this will change cursor shape while dragging).
   [view registerForDraggedTypes : supportedTypes];
   //Declared property - for convenience, not to check atoms/shmatoms or X11 properties.
   view.fIsDNDAware = YES;

   FindAtom("XdndAware", true);//Add it, if not yet.
   const Atom_t xaAtomAtom = FindAtom("XA_ATOM", false);

   assert(xaAtomAtom == 4 && "SetDNDAware, XA_ATOM is not defined");//This is a predefined atom.

   //ROOT's GUI uses Atom_t, which is unsigned long, and it's 64-bit.
   //While calling XChangeProperty, it passes the address of this typelist
   //and format is ... 32. I have to pack data into unsigned and force the size:
   assert(sizeof(unsigned) == 4 && "SetDNDAware, sizeof(unsigned) must be 4");

   std::vector<unsigned> propertyData;
   propertyData.push_back(4);//This '4' is from TGX11 (is it XA_ATOM???)

   if (typeList) {
      for (unsigned i = 0; typeList[i]; ++i)
         propertyData.push_back(unsigned(typeList[i]));//hehe.
   }

   [view setProperty : "XdndAware" data : (unsigned char *)&propertyData[0]
                size : propertyData.size() forType : xaAtomAtom format : 32];
}

//______________________________________________________________________________
Bool_t TGCocoa::IsDNDAware(Window_t windowID, Atom_t * /*typeList*/)
{
   //Checks if the Window is DND aware. typeList is ignored.

   if (windowID <= fPimpl->GetRootWindowID())//kNone or root.
      return kFALSE;

   assert(fPimpl->GetDrawable(windowID).fIsPixmap == NO &&
          "IsDNDAware, windowID parameter is not a window");

   QuartzView * const view = (QuartzView *)fPimpl->GetWindow(windowID).fContentView;
   return view.fIsDNDAware;
}

//______________________________________________________________________________
void TGCocoa::SetTypeList(Window_t, Atom_t, Atom_t *)
{
   // Add the list of drag and drop types to the Window win.
   //It's never called from GUI.
   ::Warning("SetTypeList", "Not implemented");
}

//______________________________________________________________________________
Window_t TGCocoa::FindRWindow(Window_t winID, Window_t dragWinID, Window_t inputWinID, int x, int y, int maxDepth)
{
   //Comment from TVirtualX:

   // Recursively search in the children of Window for a Window which is at
   // location x, y and is DND aware, with a maximum depth of maxd.
   // Ignore dragwin and input (???)
   //End of comment from TVirtualX.


   //Now my comments. The name of this function, as usually, says nothing about what it does.
   //It's searching for some window, probably child of winID, or may be winID itself(?) and
   //window must be DND aware. So the name should be FindDNDAwareWindowRecursively or something like this.

   //This function is not documented, comments suck as soon as they are simply wrong - the
   //first return statement in X11 version contradicts with comments
   //about child. Since X11 version is more readable, I'm reproducing X11 version here,
   //and ... my code can't be wrong, since there is nothing right about this function.

   NSView<X11Window> * const testView = X11::FindDNDAwareViewInPoint(
                                                      fPimpl->IsRootWindow(winID) ? nil : fPimpl->GetWindow(winID).fContentView,
                                                      dragWinID, inputWinID, x, y, maxDepth);
   if (testView)
      return testView.fID;

   return kNone;
}

#pragma mark - Noops.

//______________________________________________________________________________
UInt_t TGCocoa::ExecCommand(TGWin32Command * /*code*/)
{
   // Executes the command "code" coming from the other threads (Win32)
   return 0;
}

//______________________________________________________________________________
Int_t TGCocoa::GetDoubleBuffer(Int_t /*wid*/)
{
   // Queries the double buffer value for the window "wid".
   return 0;
}

//______________________________________________________________________________
void TGCocoa::GetCharacterUp(Float_t &chupx, Float_t &chupy)
{
   // Returns character up vector.
   chupx = chupy = 0.f;
}

//______________________________________________________________________________
Pixmap_t TGCocoa::ReadGIF(Int_t /*x0*/, Int_t /*y0*/, const char * /*file*/, Window_t /*id*/)
{
   // If id is NULL - loads the specified gif file at position [x0,y0] in the
   // current window. Otherwise creates pixmap from gif file

   return kNone;
}

//______________________________________________________________________________
Int_t TGCocoa::RequestLocator(Int_t /*mode*/, Int_t /*ctyp*/, Int_t &/*x*/, Int_t &/*y*/)
{
   // Requests Locator position.
   // x,y  - cursor position at moment of button press (output)
   // ctyp - cursor type (input)
   //        ctyp = 1 tracking cross
   //        ctyp = 2 cross-hair
   //        ctyp = 3 rubber circle
   //        ctyp = 4 rubber band
   //        ctyp = 5 rubber rectangle
   //
   // mode - input mode
   //        mode = 0 request
   //        mode = 1 sample
   //
   // The returned value is:
   //        in request mode:
   //                       1 = left is pressed
   //                       2 = middle is pressed
   //                       3 = right is pressed
   //        in sample mode:
   //                       11 = left is released
   //                       12 = middle is released
   //                       13 = right is released
   //                       -1 = nothing is pressed or released
   //                       -2 = leave the window
   //                     else = keycode (keyboard is pressed)

   return 0;
}

//______________________________________________________________________________
Int_t TGCocoa::RequestString(Int_t /*x*/, Int_t /*y*/, char * /*text*/)
{
   // Requests string: text is displayed and can be edited with Emacs-like
   // keybinding. Returns termination code (0 for ESC, 1 for RETURN)
   //
   // x,y  - position where text is displayed
   // text - displayed text (as input), edited text (as output)
   return 0;
}

//______________________________________________________________________________
void TGCocoa::SetCharacterUp(Float_t /*chupx*/, Float_t /*chupy*/)
{
   // Sets character up vector.
}

//______________________________________________________________________________
void TGCocoa::SetClipOFF(Int_t /*wid*/)
{
   // Turns off the clipping for the window "wid".
}

//______________________________________________________________________________
void TGCocoa::SetClipRegion(Int_t /*wid*/, Int_t /*x*/, Int_t /*y*/, UInt_t /*w*/, UInt_t /*h*/)
{
   // Sets clipping region for the window "wid".
   //
   // wid  - window indentifier
   // x, y - origin of clipping rectangle
   // w, h - the clipping rectangle dimensions

}

//______________________________________________________________________________
void TGCocoa::SetTextMagnitude(Float_t /*mgn*/)
{
   // Sets the current text magnification factor to "mgn"
}

//______________________________________________________________________________
void TGCocoa::Sync(Int_t /*mode*/)
{
   // Set synchronisation on or off.
   // mode : synchronisation on/off
   //    mode=1  on
   //    mode<>0 off
}

//______________________________________________________________________________
void TGCocoa::Warp(Int_t ix, Int_t iy, Window_t winID)
{
   // Sets the pointer position.
   // ix - new X coordinate of pointer
   // iy - new Y coordinate of pointer
   // Coordinates are relative to the origin of the window id
   // or to the origin of the current window if id == 0.

   if (!winID)
      return;

   NSPoint newCursorPosition = {};
   newCursorPosition.x = ix;
   newCursorPosition.y = iy;

   if (fPimpl->GetRootWindowID() == winID) {
      //Suddenly .... top-left - based!
      newCursorPosition.x = X11::GlobalXROOTToCocoa(newCursorPosition.x);
   } else {
      assert(fPimpl->GetDrawable(winID).fIsPixmap == NO &&
             "Warp, drawable is not a window");
      newCursorPosition = X11::TranslateToScreen(fPimpl->GetWindow(winID).fContentView,
                                                 newCursorPosition);
   }

   CGWarpMouseCursorPosition(NSPointToCGPoint(newCursorPosition));
}

//______________________________________________________________________________
Int_t TGCocoa::WriteGIF(char * /*name*/)
{
   // Writes the current window into GIF file.
   // Returns 1 in case of success, 0 otherwise.

   return 0;
}

//______________________________________________________________________________
void TGCocoa::WritePixmap(Int_t /*wid*/, UInt_t /*w*/, UInt_t /*h*/, char * /*pxname*/)
{
   // Writes the pixmap "wid" in the bitmap file "pxname".
   //
   // wid    - the pixmap address
   // w, h   - the width and height of the pixmap.
   // pxname - the file name
}

//______________________________________________________________________________
Bool_t TGCocoa::NeedRedraw(ULong_t /*tgwindow*/, Bool_t /*force*/)
{
   // Notify the low level GUI layer ROOT requires "tgwindow" to be
   // updated
   //
   // Returns kTRUE if the notification was desirable and it was sent
   //
   // At the moment only Qt4 layer needs that
   //
   // One needs explicitly cast the first parameter to TGWindow to make
   // it working in the implementation.
   //
   // One needs to process the notification to confine
   // all paint operations within "expose" / "paint" like low level event
   // or equivalent

   return kFALSE;
}

//______________________________________________________________________________
Bool_t TGCocoa::CreatePictureFromFile(Drawable_t /*wid*/,
                                      const char * /*filename*/,
                                      Pixmap_t &/*pict*/,
                                      Pixmap_t &/*pict_mask*/,
                                      PictureAttributes_t &/*attr*/)
{
   // Creates a picture pict from data in file "filename". The picture
   // attributes "attr" are used for input and output. Returns kTRUE in
   // case of success, kFALSE otherwise. If the mask "pict_mask" does not
   // exist it is set to kNone.

   return kFALSE;
}

//______________________________________________________________________________
Bool_t TGCocoa::CreatePictureFromData(Drawable_t /*wid*/, char ** /*data*/,
                                      Pixmap_t &/*pict*/,
                                      Pixmap_t &/*pict_mask*/,
                                      PictureAttributes_t & /*attr*/)
{
   // Creates a picture pict from data in bitmap format. The picture
   // attributes "attr" are used for input and output. Returns kTRUE in
   // case of success, kFALSE otherwise. If the mask "pict_mask" does not
   // exist it is set to kNone.

   return kFALSE;
}
//______________________________________________________________________________
Bool_t TGCocoa::ReadPictureDataFromFile(const char * /*filename*/, char *** /*ret_data*/)
{
   // Reads picture data from file "filename" and store it in "ret_data".
   // Returns kTRUE in case of success, kFALSE otherwise.

   return kFALSE;
}

//______________________________________________________________________________
void TGCocoa::DeletePictureData(void * /*data*/)
{
   // Delete picture data created by the function ReadPictureDataFromFile.
}

//______________________________________________________________________________
void TGCocoa::SetDashes(GContext_t /*gc*/, Int_t /*offset*/, const char * /*dash_list*/, Int_t /*n*/)
{
   // Sets the dash-offset and dash-list attributes for dashed line styles
   // in the specified GC. There must be at least one element in the
   // specified dash_list. The initial and alternating elements (second,
   // fourth, and so on) of the dash_list are the even dashes, and the
   // others are the odd dashes. Each element in the "dash_list" array
   // specifies the length (in pixels) of a segment of the pattern.
   //
   // gc        - specifies the GC (see GCValues_t structure)
   // offset    - the phase of the pattern for the dashed line-style you
   //             want to set for the specified GC.
   // dash_list - the dash-list for the dashed line-style you want to set
   //             for the specified GC
   // n         - the number of elements in dash_list
   // (see also the GCValues_t structure)
}

//______________________________________________________________________________
void TGCocoa::Bell(Int_t /*percent*/)
{
   // Sets the sound bell. Percent is loudness from -100% .. 100%.
}

//______________________________________________________________________________
void TGCocoa::WMDeleteNotify(Window_t /*wid*/)
{
   // Tells WM to send message when window is closed via WM.
}

//______________________________________________________________________________
void TGCocoa::SetClipRectangles(GContext_t /*gc*/, Int_t /*x*/, Int_t /*y*/,
                                Rectangle_t * /*recs*/, Int_t /*n*/)
{
   // Sets clipping rectangles in graphics context. [x,y] specify the origin
   // of the rectangles. "recs" specifies an array of rectangles that define
   // the clipping mask and "n" is the number of rectangles.
   // (see also the GCValues_t structure)
}

//______________________________________________________________________________
Region_t TGCocoa::CreateRegion()
{
   // Creates a new empty region.

   return 0;
}

//______________________________________________________________________________
void TGCocoa::DestroyRegion(Region_t /*reg*/)
{
   // Destroys the region "reg".
}

//______________________________________________________________________________
void TGCocoa::UnionRectWithRegion(Rectangle_t * /*rect*/, Region_t /*src*/, Region_t /*dest*/)
{
   // Updates the destination region from a union of the specified rectangle
   // and the specified source region.
   //
   // rect - specifies the rectangle
   // src  - specifies the source region to be used
   // dest - returns the destination region
}

//______________________________________________________________________________
Region_t TGCocoa::PolygonRegion(Point_t * /*points*/, Int_t /*np*/, Bool_t /*winding*/)
{
   // Returns a region for the polygon defined by the points array.
   //
   // points  - specifies an array of points
   // np      - specifies the number of points in the polygon
   // winding - specifies the winding-rule is set (kTRUE) or not(kFALSE)

   return 0;
}

//______________________________________________________________________________
void TGCocoa::UnionRegion(Region_t /*rega*/, Region_t /*regb*/, Region_t /*result*/)
{
   // Computes the union of two regions.
   //
   // rega, regb - specify the two regions with which you want to perform
   //              the computation
   // result     - returns the result of the computation

}

//______________________________________________________________________________
void TGCocoa::IntersectRegion(Region_t /*rega*/, Region_t /*regb*/, Region_t /*result*/)
{
   // Computes the intersection of two regions.
   //
   // rega, regb - specify the two regions with which you want to perform
   //              the computation
   // result     - returns the result of the computation
}

//______________________________________________________________________________
void TGCocoa::SubtractRegion(Region_t /*rega*/, Region_t /*regb*/, Region_t /*result*/)
{
   // Subtracts regb from rega and stores the results in result.
}

//______________________________________________________________________________
void TGCocoa::XorRegion(Region_t /*rega*/, Region_t /*regb*/, Region_t /*result*/)
{
   // Calculates the difference between the union and intersection of
   // two regions.
   //
   // rega, regb - specify the two regions with which you want to perform
   //              the computation
   // result     - returns the result of the computation

}

//______________________________________________________________________________
Bool_t  TGCocoa::EmptyRegion(Region_t /*reg*/)
{
   // Returns kTRUE if the region reg is empty.

   return kFALSE;
}

//______________________________________________________________________________
Bool_t  TGCocoa::PointInRegion(Int_t /*x*/, Int_t /*y*/, Region_t /*reg*/)
{
   // Returns kTRUE if the point [x, y] is contained in the region reg.

   return kFALSE;
}

//______________________________________________________________________________
Bool_t  TGCocoa::EqualRegion(Region_t /*rega*/, Region_t /*regb*/)
{
   // Returns kTRUE if the two regions have the same offset, size, and shape.

   return kFALSE;
}

//______________________________________________________________________________
void TGCocoa::GetRegionBox(Region_t /*reg*/, Rectangle_t * /*rect*/)
{
   // Returns smallest enclosing rectangle.
}

#pragma mark - Details and aux. functions.

//______________________________________________________________________________
ROOT::MacOSX::X11::EventTranslator *TGCocoa::GetEventTranslator()const
{
   return &fPimpl->fX11EventTranslator;
}

//______________________________________________________________________________
ROOT::MacOSX::X11::CommandBuffer *TGCocoa::GetCommandBuffer()const
{
   return &fPimpl->fX11CommandBuffer;
}

//______________________________________________________________________________
void TGCocoa::CocoaDrawON()
{
   ++fCocoaDraw;
}

//______________________________________________________________________________
void TGCocoa::CocoaDrawOFF()
{
   assert(fCocoaDraw > 0 && "CocoaDrawOFF, was already off");
   --fCocoaDraw;
}

//______________________________________________________________________________
bool TGCocoa::IsCocoaDraw()const
{
   return bool(fCocoaDraw);
}

//______________________________________________________________________________
void *TGCocoa::GetCurrentContext()
{
   NSObject<X11Drawable> * const drawable = fPimpl->GetDrawable(fSelectedDrawable);
   if (!drawable.fIsPixmap) {
      Error("GetCurrentContext", "TCanvas/TPad's internal error,"
                                 " selected drawable is not a pixmap!");
      return 0;
   }

   return drawable.fContext;
}

//______________________________________________________________________________
bool TGCocoa::MakeProcessForeground()
{
   //We start ROOT in a terminal window, so it's considered as a
   //background process. Background process has a lot of problems
   //if it tries to create and manage windows.
   //So, first time we convert process to foreground, next time
   //we make it front.

   if (!fForegroundProcess) {
      ProcessSerialNumber psn = {0, kCurrentProcess};

      const OSStatus res1 = TransformProcessType(&psn, kProcessTransformToForegroundApplication);

      //When TGCocoa's functions are called from the python (Apple's system version),
      //TransformProcessType fails with paramErr (looks like process is _already_ foreground),
      //why is it a paramErr - I've no idea.
      if (res1 != noErr && res1 != paramErr) {
         Error("MakeProcessForeground", "TransformProcessType failed with code %d", int(res1));
         return false;
      }
#ifdef MAC_OS_X_VERSION_10_9
      //Instead of quite transparent Carbon calls we now have another black-box function.
      [[NSApplication sharedApplication] activateIgnoringOtherApps : YES];
#else
      const OSErr res2 = SetFrontProcess(&psn);
      if (res2 != noErr) {
         Error("MakeProcessForeground", "SetFrontProcess failed with code %d", res2);
         return false;
      }
#endif

      fForegroundProcess = true;
   } else {
#ifdef MAC_OS_X_VERSION_10_9
      //Instead of quite transparent Carbon calls we now have another black-box function.
      [[NSApplication sharedApplication] activateIgnoringOtherApps : YES];
#else
      ProcessSerialNumber psn = {};

      OSErr res = GetCurrentProcess(&psn);
      if (res != noErr) {
         Error("MakeProcessForeground", "GetCurrentProcess failed with code %d", res);
         return false;
      }

      res = SetFrontProcess(&psn);
      if (res != noErr) {
         Error("MapProcessForeground", "SetFrontProcess failed with code %d", res);
         return false;
      }
#endif
   }

   return true;
}

//______________________________________________________________________________
Atom_t TGCocoa::FindAtom(const std::string &atomName, bool addIfNotFound)
{
   const std::map<std::string, Atom_t>::const_iterator it = fNameToAtom.find(atomName);

   if (it != fNameToAtom.end())
      return it->second;
   else if (addIfNotFound) {
      //Create a new atom.
      fAtomToName.push_back(atomName);
      fNameToAtom[atomName] = Atom_t(fAtomToName.size());

      return Atom_t(fAtomToName.size());
   }

   return kNone;
}

//______________________________________________________________________________
void TGCocoa::SetApplicationIcon()
{
   if (gEnv) {
      const char * const iconDirectoryPath = gEnv->GetValue("Gui.IconPath",TROOT::GetIconPath());
      if (iconDirectoryPath) {
         const Util::ScopedArray<char> fileName(gSystem->Which(iconDirectoryPath, "Root6Icon.png", kReadPermission));
         if (fileName.Get()) {
            const Util::AutoreleasePool pool;
            //Aha, ASCII ;) do not install ROOT in ...
            NSString *cocoaStr = [NSString stringWithCString : fileName.Get() encoding : NSASCIIStringEncoding];
            NSImage *image = [[[NSImage alloc] initWithContentsOfFile : cocoaStr] autorelease];
            [NSApp setApplicationIconImage : image];
         }
      }
   }
}
