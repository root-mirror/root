/* @(#)root/base:$Id$ */

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_GuiTypes
#define ROOT_GuiTypes

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// GuiTypes                                                             //
//                                                                      //
// Types used by the GUI classes.                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"

// Basic GUI types
typedef ULong_t            Handle_t;     //Generic resource handle
typedef Handle_t           Display_t;    //Display handle
typedef Handle_t           Visual_t;     //Visual handle
typedef Handle_t           Window_t;     //Window handle
typedef Handle_t           Pixmap_t;     //Pixmap handle
typedef Handle_t           Drawable_t;   //Drawable handle
typedef Handle_t           Region_t;     //Region handle
typedef Handle_t           Colormap_t;   //Colormap handle
typedef Handle_t           Cursor_t;     //Cursor handle
typedef Handle_t           FontH_t;      //Font handle (as opposed to Font_t which is an index)
typedef Handle_t           KeySym_t;     //Key symbol handle
typedef Handle_t           Atom_t;       //WM token
typedef Handle_t           GContext_t;   //Graphics context handle
typedef Handle_t           FontStruct_t; //Pointer to font structure
typedef ULong_t            Pixel_t;      //Pixel value
typedef UInt_t             Mask_t;       //Structure mask type
typedef ULong_t            Time_t;       //Event time

enum EGuiConstants {
   kNotUseful = 0, kWhenMapped = 1, kAlways = 2,
   kIsUnmapped = 0, kIsUnviewable = 1, kIsViewable = 2,
   kInputOutput = 1, kInputOnly = 2,
   kLineSolid = 0, kLineOnOffDash = 1, kLineDoubleDash = 2,
   kCapNotLast = 0, kCapButt = 1, kCapRound = 2, kCapProjecting = 3,
   kJoinMiter = 0, kJoinRound = 1, kJoinBevel = 2,
   kFillSolid = 0, kFillTiled = 1, kFillStippled = 2, kFillOpaqueStippled = 3,
   kEvenOddRule = 0, kWindingRule = 1,
   kClipByChildren = 0, kIncludeInferiors = 1,
   kArcChord = 0, kArcPieSlice = 1
};

// GUI event types. Later merge with EEventType in Button.h and rename to
// EEventTypes. Also rename in that case kGKeyPress to kKeyPress.
enum EGEventType {
   kGKeyPress, kKeyRelease, kButtonPress, kButtonRelease,
   kMotionNotify, kEnterNotify, kLeaveNotify, kFocusIn, kFocusOut,
   kExpose, kConfigureNotify, kMapNotify, kUnmapNotify, kDestroyNotify,
   kClientMessage, kSelectionClear, kSelectionRequest, kSelectionNotify,
   kColormapNotify, kButtonDoubleClick, kOtherEvent
};

enum EGraphicsFunction {
   kGXclear = 0,               // 0
   kGXand,                     // src AND dst
   kGXandReverse,              // src AND NOT dst
   kGXcopy,                    // src
   kGXandInverted,             // NOT src AND dst
   kGXnoop,                    // dst
   kGXxor,                     // src XOR dst
   kGXor,                      // src OR dst
   kGXnor,                     // NOT src AND NOT dst
   kGXequiv,                   // NOT src XOR dst
   kGXinvert,                  // NOT dst
   kGXorReverse,               // src OR NOT dst
   kGXcopyInverted,            // NOT src
   kGXorInverted,              // NOT src OR dst
   kGXnand,                    // NOT src OR NOT dst
   kGXset                      // 1
};

enum { kDefaultScrollBarWidth = 16 };

const Handle_t kNone = 0;
const Handle_t kCopyFromParent = 0;
const Handle_t kParentRelative = 1;

// Attributes that can be used when creating or changing a window
struct SetWindowAttributes_t {
   Pixmap_t   fBackgroundPixmap;     // background or kNone or kParentRelative
   ULong_t    fBackgroundPixel;      // background pixel
   Pixmap_t   fBorderPixmap;         // border of the window
   ULong_t    fBorderPixel;          // border pixel value
   UInt_t     fBorderWidth;          // border width in pixels
   Int_t      fBitGravity;           // one of bit gravity values
   Int_t      fWinGravity;           // one of the window gravity values
   Int_t      fBackingStore;         // kNotUseful, kWhenMapped, kAlways
   ULong_t    fBackingPlanes;        // planes to be preseved if possible
   ULong_t    fBackingPixel;         // value to use in restoring planes
   Bool_t     fSaveUnder;            // should bits under be saved (popups)?
   Long_t     fEventMask;            // set of events that should be saved
   Long_t     fDoNotPropagateMask;   // set of events that should not propagate
   Bool_t     fOverrideRedirect;     // boolean value for override-redirect
   Colormap_t fColormap;             // color map to be associated with window
   Cursor_t   fCursor;               // cursor to be displayed (or kNone)
   Mask_t     fMask;                 // bit mask specifying which fields are valid
};

// Window attributes that can be inquired
struct WindowAttributes_t {
   Int_t      fX, fY;                 // location of window
   Int_t      fWidth, fHeight;        // width and height of window
   Int_t      fBorderWidth;           // border width of window
   Int_t      fDepth;                 // depth of window
   void      *fVisual;                // the associated visual structure
   Window_t   fRoot;                  // root of screen containing window
   Int_t      fClass;                 // kInputOutput, kInputOnly
   Int_t      fBitGravity;            // one of bit gravity values
   Int_t      fWinGravity;            // one of the window gravity values
   Int_t      fBackingStore;          // kNotUseful, kWhenMapped, kAlways
   ULong_t    fBackingPlanes;         // planes to be preserved if possible
   ULong_t    fBackingPixel;          // value to be used when restoring planes
   Bool_t     fSaveUnder;             // boolean, should bits under be saved?
   Colormap_t fColormap;              // color map to be associated with window
   Bool_t     fMapInstalled;          // boolean, is color map currently installed
   Int_t      fMapState;              // kIsUnmapped, kIsUnviewable, kIsViewable
   Long_t     fAllEventMasks;         // set of events all people have interest in
   Long_t     fYourEventMask;         // my event mask
   Long_t     fDoNotPropagateMask;    // set of events that should not propagate
   Bool_t     fOverrideRedirect;      // boolean value for override-redirect
   void      *fScreen;                // back pointer to correct screen
};

// Bits telling which SetWindowAttributes_t fields are valid
const Mask_t kWABackPixmap       = BIT(0);
const Mask_t kWABackPixel        = BIT(1);
const Mask_t kWABorderPixmap     = BIT(2);
const Mask_t kWABorderPixel      = BIT(3);
const Mask_t kWABorderWidth      = BIT(4);
const Mask_t kWABitGravity       = BIT(5);
const Mask_t kWAWinGravity       = BIT(6);
const Mask_t kWABackingStore     = BIT(7);
const Mask_t kWABackingPlanes    = BIT(8);
const Mask_t kWABackingPixel     = BIT(9);
const Mask_t kWAOverrideRedirect = BIT(10);
const Mask_t kWASaveUnder        = BIT(11);
const Mask_t kWAEventMask        = BIT(12);
const Mask_t kWADontPropagate    = BIT(13);
const Mask_t kWAColormap         = BIT(14);
const Mask_t kWACursor           = BIT(15);

// Input event masks, used to set SetWindowAttributes_t::fEventMask
// and to be passed to TVirtualX::SelectInput()
const Mask_t kNoEventMask         = 0;
const Mask_t kKeyPressMask        = BIT(0);
const Mask_t kKeyReleaseMask      = BIT(1);
const Mask_t kButtonPressMask     = BIT(2);
const Mask_t kButtonReleaseMask   = BIT(3);
const Mask_t kPointerMotionMask   = BIT(4);
const Mask_t kButtonMotionMask    = BIT(5);
const Mask_t kExposureMask        = BIT(6);
const Mask_t kStructureNotifyMask = BIT(7);
const Mask_t kEnterWindowMask     = BIT(8);
const Mask_t kLeaveWindowMask     = BIT(9);
const Mask_t kFocusChangeMask     = BIT(10);
const Mask_t kOwnerGrabButtonMask = BIT(11);
const Mask_t kColormapChangeMask  = BIT(12);

// Event structure
struct Event_t {
   EGEventType fType;              // of event (see EGEventType)
   Window_t    fWindow;            // window reported event is relative to
   Time_t      fTime;              // time event event occured in ms
   Int_t       fX, fY;             // pointer x, y coordinates in event window
   Int_t       fXRoot, fYRoot;     // coordinates relative to root
   UInt_t      fCode;              // key or button code
   UInt_t      fState;             // key or button mask
   UInt_t      fWidth, fHeight;    // width and height of exposed area
   Int_t       fCount;             // if non-zero, at least this many more exposes
   Bool_t      fSendEvent;         // true if event came from SendEvent
   Handle_t    fHandle;            // general resource handle (used for atoms or windows)
   Int_t       fFormat;            // Next fields only used by kClientMessageEvent
   Long_t      fUser[5];           // 5 longs can be used by client message events
                                   // NOTE: only [0], [1] and [2] may be used.
                                   // [1] and [2] may contain >32 bit quantities
                                   // (i.e. pointers on 64 bit machines)
};

// Key masks, used as modifiers to GrabButton and GrabKey and
// in Event_t::fState in various key-, mouse-, and button-related events
const Mask_t kKeyShiftMask   = BIT(0);
const Mask_t kKeyLockMask    = BIT(1);
const Mask_t kKeyControlMask = BIT(2);
const Mask_t kKeyMod1Mask    = BIT(3);   // typically the Alt key
const Mask_t kKeyMod2Mask    = BIT(4);   // typically mod on numeric keys
const Mask_t kKeyMod3Mask    = BIT(5);
const Mask_t kKeyMod4Mask    = BIT(6);
const Mask_t kKeyMod5Mask    = BIT(7);
const Mask_t kButton1Mask    = BIT(8);
const Mask_t kButton2Mask    = BIT(9);
const Mask_t kButton3Mask    = BIT(10);
const Mask_t kButton4Mask    = BIT(11);
const Mask_t kButton5Mask    = BIT(12);
const Mask_t kButton6Mask    = BIT(13);
const Mask_t kButton7Mask    = BIT(14);
const Mask_t kAnyModifier    = BIT(15);

// Button names. Used as arguments to GrabButton and as Event_t::fCode
// for button events. Maps to the X11 values.
enum EMouseButton { kAnyButton, kButton1, kButton2, kButton3,
                    kButton4, kButton5, kButton6, kButton7 };

// Some magic X notify modes used in TGTextEntry widget.
// Values must match the ones in /usr/include/X11/X.h. Check when porting.
enum EXMagic { kNotifyNormal = 0, kNotifyGrab = 1, kNotifyUngrab = 2,
               kNotifyPointer = 5, kColormapUninstalled = 0,
               kColormapInstalled = 1 };

// Graphics context structure
struct GCValues_t {
   EGraphicsFunction fFunction;  // logical operation
   ULong_t  fPlaneMask;          // plane mask
   ULong_t  fForeground;         // foreground pixel
   ULong_t  fBackground;         // background pixel
   Int_t    fLineWidth;          // line width
   Int_t    fLineStyle;          // kLineSolid, kLineOnOffDash, kLineDoubleDash
   Int_t    fCapStyle;           // kCapNotLast, kCapButt,
                                 // kCapRound, kCapProjecting
   Int_t    fJoinStyle;          // kJoinMiter, kJoinRound, kJoinBevel
   Int_t    fFillStyle;          // kFillSolid, kFillTiled,
                                 // kFillStippled, kFillOpaeueStippled
   Int_t    fFillRule;           // kEvenOddRule, kWindingRule
   Int_t    fArcMode;            // kArcChord, kArcPieSlice
   Pixmap_t fTile;               // tile pixmap for tiling operations
   Pixmap_t fStipple;            // stipple 1 plane pixmap for stipping
   Int_t    fTsXOrigin;          // offset for tile or stipple operations
   Int_t    fTsYOrigin;
   FontH_t  fFont;               // default text font for text operations
   Int_t    fSubwindowMode;      // kClipByChildren, kIncludeInferiors
   Bool_t   fGraphicsExposures;  // boolean, should exposures be generated
   Int_t    fClipXOrigin;        // origin for clipping
   Int_t    fClipYOrigin;
   Pixmap_t fClipMask;           // bitmap clipping; other calls for rects
   Int_t    fDashOffset;         // patterned/dashed line information
   Char_t   fDashes[8];          // dash pattern list (dash length per byte)
   Int_t    fDashLen;            // number of dashes in fDashes
   Mask_t   fMask;               // bit mask specifying which fields are valid

   GCValues_t() : // default constructor
      fFunction  (kGXcopy),
      fPlaneMask  (0),
      fForeground  (0),
      fBackground  (1),
      fLineWidth  (0),
      fLineStyle  (kLineSolid),
      fCapStyle  (kCapButt),
      fJoinStyle  (kJoinMiter),
      fFillStyle  (kFillSolid),
      fFillRule  (kEvenOddRule),
      fArcMode  (kArcPieSlice),
      fTile  (0),
      fStipple  (0),
      fTsXOrigin  (0),
      fTsYOrigin  (0),
      fFont  (0),
      fSubwindowMode  (kClipByChildren),
      fGraphicsExposures  (kTRUE),
      fClipXOrigin  (0),
      fClipYOrigin  (0),
      fClipMask  (0),
      fDashOffset  (0),
      fDashLen  (2),
      fMask  (0)
   {
      for (int i = 2; i < 8; i++) fDashes[i] = 0;
      fDashes[0] = 5; // dashed
      fDashes[1] = 5;
   }
};

// Bits telling which GCValues_t fields are valid
const Mask_t kGCFunction          = BIT(0);
const Mask_t kGCPlaneMask         = BIT(1);
const Mask_t kGCForeground        = BIT(2);
const Mask_t kGCBackground        = BIT(3);
const Mask_t kGCLineWidth         = BIT(4);
const Mask_t kGCLineStyle         = BIT(5);
const Mask_t kGCCapStyle          = BIT(6);
const Mask_t kGCJoinStyle         = BIT(7);
const Mask_t kGCFillStyle         = BIT(8);
const Mask_t kGCFillRule          = BIT(9);
const Mask_t kGCTile              = BIT(10);
const Mask_t kGCStipple           = BIT(11);
const Mask_t kGCTileStipXOrigin   = BIT(12);
const Mask_t kGCTileStipYOrigin   = BIT(13);
const Mask_t kGCFont              = BIT(14);
const Mask_t kGCSubwindowMode     = BIT(15);
const Mask_t kGCGraphicsExposures = BIT(16);
const Mask_t kGCClipXOrigin       = BIT(17);
const Mask_t kGCClipYOrigin       = BIT(18);
const Mask_t kGCClipMask          = BIT(19);
const Mask_t kGCDashOffset        = BIT(20);
const Mask_t kGCDashList          = BIT(21);
const Mask_t kGCArcMode           = BIT(22);

struct ColorStruct_t {
   ULong_t   fPixel;    // color pixel value (index in color table)
   UShort_t  fRed;      // red component (0..65535)
   UShort_t  fGreen;    // green component (0..65535)
   UShort_t  fBlue;     // blue component (0..65535)
   UShort_t  fMask;     // mask telling which color components are valid
};

// Bits telling which ColorStruct_t fields are valid
const Mask_t kDoRed   = BIT(0);
const Mask_t kDoGreen = BIT(1);
const Mask_t kDoBlue  = BIT(2);

struct PictureAttributes_t {
   Colormap_t   fColormap;   // colormap to use
   Int_t        fDepth;      // depth of window
   UInt_t       fWidth;      // width of picture
   UInt_t       fHeight;     // height of picture
   UInt_t       fXHotspot;   // picture x hotspot coordinate
   UInt_t       fYHotspot;   // picture y hotspot coordinate
   ULong_t     *fPixels;     // list of used color pixels (if set use delete[])
   UInt_t       fNpixels;    // number of used color pixels
   UInt_t       fCloseness;  // allowable RGB deviation
   Mask_t       fMask;       // mask specifying which attributes are defined
};

// PictureAttributes_t masks bits
const Mask_t kPAColormap     = BIT(0);
const Mask_t kPADepth        = BIT(1);
const Mask_t kPASize         = BIT(2);   // width and height
const Mask_t kPAHotspot      = BIT(3);   // x and y hotspot
const Mask_t kPAReturnPixels = BIT(4);
const Mask_t kPACloseness    = BIT(5);

// Initial window mapping state
enum EInitialState {
   kNormalState = BIT(0),
   kIconicState = BIT(1)
};

// Used for drawing line segments (maps to the X11 XSegments structure)
struct Segment_t {
   Short_t fX1, fY1, fX2, fY2;
};

// Point structure (maps to the X11 XPoint structure)
struct Point_t {
   Short_t fX, fY;
};

// Rectangle structure (maps to the X11 XRectangle structure)
struct Rectangle_t {
   Short_t  fX, fY;
   UShort_t fWidth, fHeight;
};

// Atoms used for text cut and paste between windows
const Atom_t kPrimarySelection = 1;  // magic values, must match the ones
const Atom_t kCutBuffer        = 9;  // in /usr/include/X11/Xatom.h
const Int_t  kMaxPixel         = 32000;

const int kNumCursors = 19;
enum ECursor { kBottomLeft, kBottomRight, kTopLeft, kTopRight,
               kBottomSide, kLeftSide, kTopSide, kRightSide,
               kMove, kCross, kArrowHor, kArrowVer, kHand, kRotate,
               kPointer, kArrowRight, kCaret, kWatch, kNoDrop };

//---- types of frames (and borders), moved from TGFrame.h

enum EFrameType {
   kChildFrame      = 0,
   kMainFrame       = BIT(0),
   kVerticalFrame   = BIT(1),
   kHorizontalFrame = BIT(2),
   kSunkenFrame     = BIT(3),
   kRaisedFrame     = BIT(4),
   kDoubleBorder    = BIT(5),
   kFitWidth        = BIT(6),
   kFixedWidth      = BIT(7),
   kFitHeight       = BIT(8),
   kFixedHeight     = BIT(9),
   kFixedSize       = (kFixedWidth | kFixedHeight),
   kOwnBackground   = BIT(10),
   kTransientFrame  = BIT(11),
   kTempFrame       = BIT(12),
   kMdiMainFrame    = BIT(13),
   kMdiFrame        = BIT(14)
};


#endif
