// @(#)root/graf2d:$Id$
// Author: Timur Pocheptsov   16/02/2012

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_QuartzWindow
#define ROOT_QuartzWindow

#include <Cocoa/Cocoa.h>

#include "CocoaGuiTypes.h"
#include "X11Drawable.h"
#include "X11Events.h"
#include "GuiTypes.h"

/////////////////////////////////////////////////////////////////////
//                                                                 //
// CrosshairView comprises the content view of a CrosshairWindow.  //
// It's purpose is to render two lines ... a crosshair.            //
//                                                                 //
/////////////////////////////////////////////////////////////////////
@interface CrosshairView: NSView
// Line1:
@property (nonatomic, assign)NSPoint start1;
@property (nonatomic, assign)NSPoint end1;
// Line2:
@property (nonatomic, assign)NSPoint start2;
@property (nonatomic, assign)NSPoint end2;
@end

// CrosshairWindow is a special window: a transparent
// transient child window that we attach to a canvas
// to draw a crosshair on top of the pad's contents.
// It's transparent to all mouse events and can never
// be main or a key window. It has a transparent
// background.
@interface CrosshairWindow : NSWindow
- (instancetype) init;
@end

////////////////////////////////////////////////
//                                            //
// QuartzWindow class : top-level window.     //
//                                            //
////////////////////////////////////////////////

@class ROOTOpenGLView;
@class QuartzImage;

@interface QuartzWindow : NSWindow<X11Window, NSWindowDelegate> {
@private
   QuartzWindow *fMainWindow;
   BOOL fHasFocus;
   
   QuartzView *fContentView;
   BOOL fDelayedTransient;
   QuartzImage *fShapeCombineMask;
   BOOL fIsDeleted;
}

- (id) initWithContentRect : (NSRect) contentRect styleMask : (NSUInteger) windowStyle
       backing : (NSBackingStoreType) bufferingType defer : (BOOL) deferCreation
       windowAttributes : (const SetWindowAttributes_t *) attr;

- (id) initWithGLView : (ROOTOpenGLView *) glView;

- (void) dealloc;

//With reference counting and autorelease pools, it's possible that
//TGCocoa::DestroyWindow was called and window was correctly deleted,
//but it's still on screen and if used in some functions (like FindWindowForPointerEvent)
//and this ends in a segmentation fault.
//fIsDeleted property is here to solve this problem.
- (BOOL) fIsDeleted;
- (void) setFIsDeleted : (BOOL) deleted;

//Many properties in QuartzWindow just forwards to fContentView.
- (void) forwardInvocation : (NSInvocation *) anInvocation;
- (NSMethodSignature*) methodSignatureForSelector : (SEL) selector;

//This is to emulate "transient" window/main window relationship:
@property (nonatomic, assign) QuartzWindow *fMainWindow;
- (void) addTransientWindow : (QuartzWindow *) window;

//Shape mask - non-rectangular window.
@property (nonatomic, assign) QuartzImage *fShapeCombineMask;
//@property (nonatomic, assign) NSPoint fShapeMaskShift;

//1. X11Drawable protocol.

- (BOOL) fIsPixmap;
- (BOOL) fIsOpenGLWidget;
- (CGFloat) fScaleFactor;

//Geometry.
- (int) fX;
- (int) fY;

- (unsigned) fWidth;
- (unsigned) fHeight;

- (void) setDrawableSize : (NSSize) newSize;
- (void) setX : (int) x Y : (int) y width : (unsigned) w height : (unsigned) h;
- (void) setX : (int) x Y : (int) y;

//
- (void) copy : (NSObject<X11Drawable> *) src area : (ROOT::MacOSX::X11::Rectangle) area withMask : (QuartzImage *) mask
         clipOrigin : (ROOT::MacOSX::X11::Point) origin toPoint : (ROOT::MacOSX::X11::Point) dstPoint;

- (unsigned char *) readColorBits : (ROOT::MacOSX::X11::Rectangle) area;

// Trick for crosshair drawing in TCanvas ("pseudo-XOR")
- (void) addCrosshairWindow;
- (void) adjustCrosshairWindowGeometry;
- (void) adjustCrosshairWindowGeometry : (CrosshairWindow *)win;
- (void) removeCrosshairWindow;
- (CrosshairWindow *) findCrosshairWindow;

//X11Window protocol.

/////////////////////////////////////////////////////////////////
//SetWindowAttributes_t/WindowAttributes_t
@property (nonatomic, assign) unsigned long fBackgroundPixel;
@property (nonatomic, readonly) int         fMapState;

//End of SetWindowAttributes_t/WindowAttributes_t
/////////////////////////////////////////////////////////////////

@property (nonatomic, assign) BOOL          fHasFocus;

//"Back buffer" is a bitmap, attached to a window by TCanvas.
@property (nonatomic, assign) QuartzView            *fParentView;
@property (nonatomic, readonly) NSView<X11Window>   *fContentView;
@property (nonatomic, readonly) QuartzWindow        *fQuartzWindow;

//Children subviews.
- (void) addChild : (NSView<X11Window> *) child;

//X11/ROOT GUI's attributes.
- (void) getAttributes : (WindowAttributes_t *) attr;
- (void) setAttributes : (const SetWindowAttributes_t *) attr;

//X11's XMapWindow etc.
- (void) mapRaised;
- (void) mapWindow;
- (void) mapSubwindows;
- (void) unmapWindow;

@end

//////////////////////////////////////////////////////////////
//                                                          //
// I have to attach passive key grabs to a view.            //
//                                                          //
//////////////////////////////////////////////////////////////

@interface PassiveKeyGrab : NSObject {
@private
   unichar fKeyCode;
   NSUInteger fModifiers;
}
- (unichar)    fKeyCode;
- (NSUInteger) fModifiers;
- (id) initWithKey : (unichar) keyCode modifiers : (NSUInteger) modifiers;
- (BOOL) matchKey : (unichar) keyCode modifiers : (NSUInteger) modifiers;
- (BOOL) matchKey : (unichar) keyCode;
@end

////////////////////////////////////////
//                                    //
// QuartzView class - child window.   //
//                                    //
////////////////////////////////////////

@class QuartzImage;

@interface QuartzView : NSView<X11Window> {
@protected
   unsigned fID;
   CGContextRef fContext;
   long fEventMask;
   int fClass;
   int fDepth;
   int fBitGravity;
   int fWinGravity;
   unsigned long fBackgroundPixel;
   BOOL fOverrideRedirect;

   BOOL fHasFocus;
   QuartzView *fParentView;

   int fPassiveGrabButton;
   unsigned fPassiveGrabEventMask;
   unsigned fPassiveGrabKeyModifiers;
   unsigned fActiveGrabEventMask;
   BOOL fPassiveGrabOwnerEvents;
   BOOL fSnapshotDraw;
   ECursor fCurrentCursor;
   BOOL fIsDNDAware;

   QuartzPixmap   *fBackBuffer;
   NSMutableArray *fPassiveKeyGrabs;
   BOOL            fIsOverlapped;

   NSMutableDictionary   *fX11Properties;
   QuartzImage           *fBackgroundPixmap;

   ROOT::MacOSX::X11::PointerGrab fCurrentGrabType;

   BOOL             fActiveGrabOwnerEvents;
}

//Life-cycle.
- (id) initWithFrame : (NSRect) frame windowAttributes : (const SetWindowAttributes_t *) attr;

//X11Drawable protocol.

@property (nonatomic, assign) unsigned fID;

- (BOOL) fIsPixmap;
- (BOOL) fIsOpenGLWidget;
- (CGFloat) fScaleFactor;

@property (nonatomic, assign) CGContextRef fContext;

//Geometry.
- (int)      fX;
- (int)      fY;
- (unsigned) fWidth;
- (unsigned) fHeight;
- (void)     setDrawableSize : (NSSize) newSize;
- (void)     setX : (int) x Y : (int) y width : (unsigned) w height : (unsigned) h;
- (void)     setX : (int) x Y : (int) y;

- (void)     copy : (NSObject<X11Drawable> *) src area : (ROOT::MacOSX::X11::Rectangle) area withMask : (QuartzImage *)mask
             clipOrigin : (ROOT::MacOSX::X11::Point) origin toPoint : (ROOT::MacOSX::X11::Point) dstPoint;
- (unsigned char *) readColorBits : (ROOT::MacOSX::X11::Rectangle) area;

//X11Window protocol.

/////////////////////////////////////////////////////////////////
//SetWindowAttributes_t/WindowAttributes_t

@property (nonatomic, assign) long          fEventMask;
@property (nonatomic, assign) int           fClass;
@property (nonatomic, assign) int           fDepth;
@property (nonatomic, assign) int           fBitGravity;
@property (nonatomic, assign) int           fWinGravity;
@property (nonatomic, assign) unsigned long fBackgroundPixel;
@property (nonatomic, retain) QuartzImage  *fBackgroundPixmap;
@property (nonatomic, readonly) int         fMapState;
@property (nonatomic, assign) BOOL          fOverrideRedirect;

//End of SetWindowAttributes_t/WindowAttributes_t
/////////////////////////////////////////////////////////////////

@property (nonatomic, assign) BOOL          fHasFocus;


@property (nonatomic, retain) QuartzPixmap        *fBackBuffer;
@property (nonatomic, assign) QuartzView          *fParentView;
@property (nonatomic, readonly) NSView<X11Window> *fContentView;
@property (nonatomic, readonly) QuartzWindow      *fQuartzWindow;

@property (nonatomic, assign) int      fPassiveGrabButton;
@property (nonatomic, assign) unsigned fPassiveGrabEventMask;
@property (nonatomic, assign) unsigned fPassiveGrabKeyModifiers;

@property (nonatomic, assign) BOOL     fPassiveGrabOwnerEvents;

- (void) activatePassiveGrab;
- (void) activateImplicitGrab;
- (void) activateGrab : (unsigned) eventMask ownerEvents : (BOOL) ownerEvents;
- (void) cancelGrab;

- (BOOL) acceptsCrossingEvents : (unsigned) eventMask;

//Children subviews.
- (void) addChild : (NSView<X11Window> *)child;

//X11/ROOT GUI's attributes.
- (void) getAttributes : (WindowAttributes_t *) attr;
- (void) setAttributes : (const SetWindowAttributes_t *) attr;

- (void) mapRaised;
- (void) mapWindow;
- (void) mapSubwindows;

- (void) unmapWindow;

- (void) raiseWindow;
- (void) lowerWindow;

- (BOOL) fIsOverlapped;
- (void) setOverlapped : (BOOL) overlap;
- (void) configureNotifyTree;

//Additional methods and properties.

@property (nonatomic, assign) BOOL   fSnapshotDraw;
- (BOOL) isFlipped;//override method from NSView.

//Keyboard:
- (void) addPassiveKeyGrab : (unichar) keyCode modifiers : (NSUInteger) modifiers;
- (void) removePassiveKeyGrab : (unichar) keyCode modifiers : (NSUInteger) modifiers;
- (PassiveKeyGrab *) findPassiveKeyGrab : (unichar) keyCode modifiers : (NSUInteger) modifiers;
- (PassiveKeyGrab *) findPassiveKeyGrab : (unichar) keyCode;

//Cursors.
@property (nonatomic, assign) ECursor fCurrentCursor;

//X11 "properties".
- (void) setProperty : (const char *) propName data : (unsigned char *) propData size : (unsigned) dataSize
         forType : (Atom_t) dataType format : (unsigned) format;
- (BOOL) hasProperty : (const char *) propName;
- (unsigned char *) getProperty : (const char *) propName returnType : (Atom_t *) type
   returnFormat : (unsigned *) format nElements : (unsigned *) nElements;
- (void) removeProperty : (const char *) propName;

//DND
@property (nonatomic, assign) BOOL fIsDNDAware;

- (NSDragOperation) draggingEntered : (id<NSDraggingInfo>) sender;
- (BOOL) performDragOperation : (id<NSDraggingInfo>) sender;

@end

namespace ROOT {
namespace MacOSX {
namespace X11 {

QuartzWindow *CreateTopLevelWindow(Int_t x, Int_t y, UInt_t w, UInt_t h, UInt_t border, Int_t depth,
                                   UInt_t clss, void *visual, SetWindowAttributes_t *attr, UInt_t);
QuartzView *CreateChildView(QuartzView *parent, Int_t x, Int_t y, UInt_t w, UInt_t h, UInt_t border, Int_t depth,
                            UInt_t clss, void *visual, SetWindowAttributes_t *attr, UInt_t wtype);

void GetRootWindowAttributes(WindowAttributes_t *attr);
void GetWindowAttributes(NSObject<X11Window> *window, WindowAttributes_t *dst);

//Coordinate conversion.

//This two functions operate with Cocoa's coordinate system (so, 'to screen' will return Cocoa's
//point, and 'from screen' expects Cocoa's point (not ROOT)).
NSPoint ConvertPointFromBaseToScreen(NSWindow *window, NSPoint windowPoint);
NSPoint ConvertPointFromScreenToBase(NSPoint screenPoint, NSWindow *window);

int GlobalXCocoaToROOT(CGFloat xCocoa);
int GlobalYCocoaToROOT(CGFloat yCocoa);
int GlobalXROOTToCocoa(CGFloat xROOT);
int GlobalYROOTToCocoa(CGFloat yROOT);

int LocalYCocoaToROOT(NSView<X11Window> *parentView, CGFloat yCocoa);
int LocalYROOTToCocoa(NSView<X11Window> *parentView, CGFloat yROOT);
int LocalYROOTToCocoa(NSObject<X11Drawable> *parentView, CGFloat yROOT);

NSPoint TranslateToScreen(NSView<X11Window> *from, NSPoint point);
NSPoint TranslateFromScreen(NSPoint point, NSView<X11Window> *to);
NSPoint TranslateCoordinates(NSView<X11Window> *fromView, NSView<X11Window> *toView, NSPoint sourcePoint);

bool ViewIsTextViewFrame(NSView<X11Window> *view, bool checkParent);
bool ViewIsHtmlViewFrame(NSView<X11Window> *view, bool checkParent);
bool LockFocus(NSView<X11Window> *view);
void UnlockFocus(NSView<X11Window> *view);

bool ScreenPointIsInView(NSView<X11Window> *view, Int_t x, Int_t y);
QuartzWindow *FindWindowInPoint(Int_t x, Int_t y);
NSView<X11Window> *FindDNDAwareViewInPoint(NSView *parentView, Window_t dragWinID, Window_t inputWinID, Int_t x, Int_t y, Int_t maxDepth);

//Pointer == cursor in X11's terms.

//These two functions use "mouse location outside of event stream" - simply
//asks for the current cursor location
//("regardless of the current event being handled or of any events pending").
QuartzWindow *FindWindowUnderPointer();
NSView<X11Window> *FindViewUnderPointer();

//These two functions use coordinates from the event to find a window/view.
QuartzWindow *FindWindowForPointerEvent(NSEvent *pointerEvent);
NSView<X11Window> *FindViewForPointerEvent(NSEvent *pointerEvent);
void WindowLostFocus(Window_t winID);

//Add shape mask to context.
void ClipToShapeMask(NSView<X11Window> *view, CGContextRef ctx);

}//X11
}//MacOSX
}//ROOT

#endif
