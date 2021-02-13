// @(#)root/x11:$Id$
// Author: Rene Brun, Olivier Couet, Fons Rademakers   28/11/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/// \defgroup x11 X11 backend
/// \brief Interface to X11 graphics.
/// \ingroup GraphicsBackends

/** \class TGX11
\ingroup x11
This class is the basic interface to the X11 (Xlib) graphics system.
It is an implementation of the abstract TVirtualX class.

This class gives access to basic X11 graphics, pixmap, text and font handling
routines.

The companion class for Win32 is TGWin32.

The file G11Gui.cxx contains the implementation of the GUI methods of the
TGX11 class. Most of the methods are used by the machine independent
GUI classes (libGUI.so).

This code was initially developed in the context of HIGZ and PAW
by Olivier Couet (package X11INT).
*/

#include "TROOT.h"
#include "TColor.h"
#include "TGX11.h"
#include "TPoint.h"
#include "TMath.h"
#include "TStorage.h"
#include "TStyle.h"
#include "TExMap.h"
#include "TEnv.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "RStipples.h"
#include "strlcpy.h"

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/cursorfont.h>
#include <X11/keysym.h>
#include <X11/xpm.h>

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <unistd.h>
#ifdef R__AIX
#   include <sys/socket.h>
#endif

extern float   XRotVersion(char*, int);
extern void    XRotSetMagnification(float);
extern void    XRotSetBoundingBoxPad(int);
extern int     XRotDrawString(Display*, XFontStruct*, float,
                              Drawable, GC, int, int, char*);
extern int     XRotDrawImageString(Display*, XFontStruct*, float,
                                   Drawable, GC, int, int, char*);
extern int     XRotDrawAlignedString(Display*, XFontStruct*, float,
                                     Drawable, GC, int, int, char*, int);
extern int     XRotDrawAlignedImageString(Display*, XFontStruct*, float,
                                          Drawable, GC, int, int, char*, int);
extern XPoint *XRotTextExtents(Display*, XFontStruct*, float,
                               int, int, char*, int);

//---- globals

static XWindow_t *gCws;      // gCws: pointer to the current window
static XWindow_t *gTws;      // gTws: temporary pointer

const Int_t kBIGGEST_RGB_VALUE = 65535;

//
// Primitives Graphic Contexts global for all windows
//
const int kMAXGC = 7;
static GC gGClist[kMAXGC];
static GC *gGCline = &gGClist[0];  // PolyLines
static GC *gGCmark = &gGClist[1];  // PolyMarker
static GC *gGCfill = &gGClist[2];  // Fill areas
static GC *gGCtext = &gGClist[3];  // Text
static GC *gGCinvt = &gGClist[4];  // Inverse text
static GC *gGCdash = &gGClist[5];  // Dashed lines
static GC *gGCpxmp = &gGClist[6];  // Pixmap management

static GC gGCecho;                 // Input echo

static Int_t  gFillHollow;         // Flag if fill style is hollow
static Pixmap gFillPattern = 0;    // Fill pattern

//
// Text management
//
const Int_t kMAXFONT = 4;
static struct {
   XFontStruct *id;
   char         name[80];                    // Font name
} gFont[kMAXFONT];                          // List of fonts loaded

static XFontStruct *gTextFont;              // Current font
static Int_t        gCurrentFontNumber = 0; // Current font number in gFont[]

//
// Markers
//
const Int_t kMAXMK = 100;
static struct {
   int    type;
   int    n;
   XPoint xy[kMAXMK];
} gMarker;                        // Point list to draw marker
static int  gMarkerLineWidth = 0;
static int  gMarkerLineStyle = LineSolid;
static int  gMarkerCapStyle  = CapRound;
static int  gMarkerJoinStyle = JoinRound;

//
// Keep style values for line GC
//
static int  gLineWidth = 0;
static int  gLineStyle = LineSolid;
static int  gCapStyle  = CapButt;
static int  gJoinStyle = JoinMiter;
static char gDashList[10];
static int  gDashLength = 0;
static int  gDashOffset = 0;
static int  gDashSize   = 0;

//
// Event masks
//
static ULong_t gMouseMask =   ButtonPressMask   | ButtonReleaseMask |
                              EnterWindowMask   | LeaveWindowMask   |
                              PointerMotionMask | KeyPressMask      |
                              KeyReleaseMask;
static ULong_t gKeybdMask =   ButtonPressMask | KeyPressMask |
                              EnterWindowMask | LeaveWindowMask;

//
// Data to create an invisible cursor
//
const char null_cursor_bits[] = {
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
static Cursor gNullCursor = 0;

struct RXGCValues:XGCValues{};
struct RXColor:XColor{};
struct RXImage:XImage{};
struct RXPoint:XPoint{};
struct RXVisualInfo:XVisualInfo{};
struct RVisual:Visual{};

ClassImp(TGX11);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TGX11::TGX11()
{
   int i;
   fDisplay            = 0;
   fScreenNumber       = 0;
   fVisual             = 0;
   fRootWin            = 0;
   fVisRootWin         = 0;
   fColormap           = 0;
   fBlackPixel         = 0;
   fWhitePixel         = 0;
   fWindows            = 0;
   fColors             = 0;
   fXEvent             = new XEvent;
   fRedDiv             = -1;
   fGreenDiv           = -1;
   fBlueDiv            = -1;
   fRedShift           = -1;
   fGreenShift         = -1;
   fBlueShift          = -1;
   fCharacterUpX       = 1;
   fCharacterUpY       = 1;
   fDepth              = 0;
   fHasTTFonts         = kFALSE;
   fHasXft             = kFALSE;
   fMaxNumberOfWindows = 10;
   fTextAlignH         = 1;
   fTextAlignV         = 1;
   fTextAlign          = 7;
   fTextMagnitude      = 1;
   for (i = 0; i < kNumCursors; i++) fCursors[i] = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal Constructor.

TGX11::TGX11(const char *name, const char *title) : TVirtualX(name, title)
{
   int i;
   fDisplay            = 0;
   fScreenNumber       = 0;
   fVisual             = 0;
   fRootWin            = 0;
   fVisRootWin         = 0;
   fColormap           = 0;
   fBlackPixel         = 0;
   fWhitePixel         = 0;
   fDrawMode           = kCopy;
   fXEvent             = new XEvent;
   fRedDiv             = -1;
   fGreenDiv           = -1;
   fBlueDiv            = -1;
   fRedShift           = -1;
   fGreenShift         = -1;
   fBlueShift          = -1;
   fCharacterUpX       = 1;
   fCharacterUpY       = 1;
   fDepth              = 0;
   fHasTTFonts         = kFALSE;
   fHasXft             = kFALSE;
   fMaxNumberOfWindows = 10;
   fTextAlignH         = 1;
   fTextAlignV         = 1;
   fTextAlign          = 7;
   fTextMagnitude      = 1;
   for (i = 0; i < kNumCursors; i++) fCursors[i] = 0;

   //fWindows = new XWindow_t[fMaxNumberOfWindows];
   fWindows = (XWindow_t*) TStorage::Alloc(fMaxNumberOfWindows*sizeof(XWindow_t));
   for (i = 0; i < fMaxNumberOfWindows; i++)
      fWindows[i].fOpen = 0;

   fColors = new TExMap;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor. Currently only used by TGX11TTF.

TGX11::TGX11(const TGX11 &org) : TVirtualX(org)
{
   int i;

   fDisplay         = org.fDisplay;
   fScreenNumber    = org.fScreenNumber;
   fVisual          = org.fVisual;
   fRootWin         = org.fRootWin;
   fVisRootWin      = org.fVisRootWin;
   fColormap        = org.fColormap;
   fBlackPixel      = org.fBlackPixel;
   fWhitePixel      = org.fWhitePixel;
   fHasTTFonts      = org.fHasTTFonts;
   fHasXft          = org.fHasXft;
   fTextAlignH      = org.fTextAlignH;
   fTextAlignV      = org.fTextAlignV;
   fTextAlign       = org.fTextAlign;
   fTextMagnitude   = org.fTextMagnitude;
   fCharacterUpX    = org.fCharacterUpX;
   fCharacterUpY    = org.fCharacterUpY;
   fDepth           = org.fDepth;
   fRedDiv          = org.fRedDiv;
   fGreenDiv        = org.fGreenDiv;
   fBlueDiv         = org.fBlueDiv;
   fRedShift        = org.fRedShift;
   fGreenShift      = org.fGreenShift;
   fBlueShift       = org.fBlueShift;
   fDrawMode        = org.fDrawMode;
   fXEvent          = new XEvent;

   fMaxNumberOfWindows = org.fMaxNumberOfWindows;
   //fWindows = new XWindow_t[fMaxNumberOfWindows];
   fWindows = (XWindow_t*) TStorage::Alloc(fMaxNumberOfWindows*sizeof(XWindow_t));
   for (i = 0; i < fMaxNumberOfWindows; i++) {
      fWindows[i].fOpen         = org.fWindows[i].fOpen;
      fWindows[i].fDoubleBuffer = org.fWindows[i].fDoubleBuffer;
      fWindows[i].fIsPixmap     = org.fWindows[i].fIsPixmap;
      fWindows[i].fDrawing      = org.fWindows[i].fDrawing;
      fWindows[i].fWindow       = org.fWindows[i].fWindow;
      fWindows[i].fBuffer       = org.fWindows[i].fBuffer;
      fWindows[i].fWidth        = org.fWindows[i].fWidth;
      fWindows[i].fHeight       = org.fWindows[i].fHeight;
      fWindows[i].fClip         = org.fWindows[i].fClip;
      fWindows[i].fXclip        = org.fWindows[i].fXclip;
      fWindows[i].fYclip        = org.fWindows[i].fYclip;
      fWindows[i].fWclip        = org.fWindows[i].fWclip;
      fWindows[i].fHclip        = org.fWindows[i].fHclip;
      fWindows[i].fNewColors    = org.fWindows[i].fNewColors;
      fWindows[i].fNcolors      = org.fWindows[i].fNcolors;
      fWindows[i].fShared       = org.fWindows[i].fShared;
   }

   for (i = 0; i < kNumCursors; i++)
      fCursors[i] = org.fCursors[i];

   fColors = new TExMap;
   Long64_t key, value;
   TExMapIter it(org.fColors);
   while (it.Next(key, value)) {
      XColor_t *colo = (XColor_t *) (Long_t)value;
      XColor_t *col  = new XColor_t;
      col->fPixel   = colo->fPixel;
      col->fRed     = colo->fRed;
      col->fGreen   = colo->fGreen;
      col->fBlue    = colo->fBlue;
      col->fDefined = colo->fDefined;
      fColors->Add(key, (Long_t) col);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGX11::~TGX11()
{
   delete (XEvent*)fXEvent;
   if (fWindows) TStorage::Dealloc(fWindows);

   if (!fColors) return;
   Long64_t key, value;
   TExMapIter it(fColors);
   while (it.Next(key, value)) {
      XColor_t *col = (XColor_t *) (Long_t)value;
      delete col;
   }
   delete fColors;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize X11 system. Returns kFALSE in case of failure.

Bool_t TGX11::Init(void *display)
{
   if (OpenDisplay((Display *) display) == -1) return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Allocate color in colormap. If we are on an <= 8 plane machine
/// we will use XAllocColor. If we are on a >= 15 (15, 16 or 24) plane
/// true color machine we will calculate the pixel value using:
/// for 15 and 16 bit true colors have 6 bits precision per color however
/// only the 5 most significant bits are used in the color index.
/// Except for 16 bits where green uses all 6 bits. I.e.:
/// ~~~ {.cpp}
///   15 bits = rrrrrgggggbbbbb
///   16 bits = rrrrrggggggbbbbb
/// ~~~
/// for 24 bits each r, g and b are represented by 8 bits.
///
/// Since all colors are set with a max of 65535 (16 bits) per r, g, b
/// we just right shift them by 10, 11 and 10 bits for 16 planes, and
/// (10, 10, 10 for 15 planes) and by 8 bits for 24 planes.
/// Returns kFALSE in case color allocation failed.

Bool_t TGX11::AllocColor(Colormap cmap, RXColor *color)
{
   if (fRedDiv == -1) {
      if (XAllocColor((Display*)fDisplay, cmap, color))
         return kTRUE;
   } else {
      color->pixel = (color->red   >> fRedDiv)   << fRedShift |
                     (color->green >> fGreenDiv) << fGreenShift |
                     (color->blue  >> fBlueDiv)  << fBlueShift;
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the current RGB value for the pixel in the XColor structure.

void TGX11::QueryColors(Colormap cmap, RXColor *color, Int_t ncolors)
{
   if (fRedDiv == -1) {
      XQueryColors((Display*)fDisplay, cmap, color, ncolors);
   } else {
      ULong_t r, g, b;
      for (Int_t i = 0; i < ncolors; i++) {
         r = (color[i].pixel & fVisual->red_mask) >> fRedShift;
         color[i].red = UShort_t(r*kBIGGEST_RGB_VALUE/(fVisual->red_mask >> fRedShift));

         g = (color[i].pixel & fVisual->green_mask) >> fGreenShift;
         color[i].green = UShort_t(g*kBIGGEST_RGB_VALUE/(fVisual->green_mask >> fGreenShift));

         b = (color[i].pixel & fVisual->blue_mask) >> fBlueShift;
         color[i].blue = UShort_t(b*kBIGGEST_RGB_VALUE/(fVisual->blue_mask >> fBlueShift));

         color[i].flags = DoRed | DoGreen | DoBlue;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the pixmap pix.

void TGX11::ClearPixmap(Drawable *pix)
{
   Window root;
   int xx, yy;
   unsigned int ww, hh, border, depth;
   XGetGeometry((Display*)fDisplay, *pix, &root, &xx, &yy, &ww, &hh, &border, &depth);
   SetColor(gGCpxmp, 0);
   XFillRectangle((Display*)fDisplay, *pix, *gGCpxmp, 0 ,0 ,ww ,hh);
   SetColor(gGCpxmp, 1);
   XFlush((Display*)fDisplay);
}

////////////////////////////////////////////////////////////////////////////////
/// Clear current window.

void TGX11::ClearWindow()
{
   if (!gCws->fIsPixmap && !gCws->fDoubleBuffer) {
      XSetWindowBackground((Display*)fDisplay, gCws->fDrawing, GetColor(0).fPixel);
      XClearWindow((Display*)fDisplay, gCws->fDrawing);
      XFlush((Display*)fDisplay);
   } else {
      SetColor(gGCpxmp, 0);
      XFillRectangle((Display*)fDisplay, gCws->fDrawing, *gGCpxmp,
                     0, 0, gCws->fWidth, gCws->fHeight);
      SetColor(gGCpxmp, 1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete current pixmap.

void TGX11::ClosePixmap()
{
   CloseWindow1();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete current window.

void TGX11::CloseWindow()
{
   if (gCws->fShared)
      gCws->fOpen = 0;
   else
      CloseWindow1();

   // Never close connection. TApplication takes care of that
   //   if (!gCws) Close();    // close X when no open window left
}

////////////////////////////////////////////////////////////////////////////////
/// Delete current window.

void TGX11::CloseWindow1()
{
   int wid;

   if (gCws->fIsPixmap)
      XFreePixmap((Display*)fDisplay, gCws->fWindow);
   else
      XDestroyWindow((Display*)fDisplay, gCws->fWindow);

   if (gCws->fBuffer) XFreePixmap((Display*)fDisplay, gCws->fBuffer);

   if (gCws->fNewColors) {
      if (fRedDiv == -1)
         XFreeColors((Display*)fDisplay, fColormap, gCws->fNewColors, gCws->fNcolors, 0);
      delete [] gCws->fNewColors;
      gCws->fNewColors = 0;
   }

   XFlush((Display*)fDisplay);

   gCws->fOpen = 0;

   // make first window in list the current window
   for (wid = 0; wid < fMaxNumberOfWindows; wid++)
      if (fWindows[wid].fOpen) {
         gCws = &fWindows[wid];
         return;
      }

   gCws = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the pixmap wid at the position xpos, ypos in the current window.

void TGX11::CopyPixmap(int wid, int xpos, int ypos)
{
   gTws = &fWindows[wid];

   XCopyArea((Display*)fDisplay, gTws->fDrawing, gCws->fDrawing, *gGCpxmp, 0, 0, gTws->fWidth,
             gTws->fHeight, xpos, ypos);
   XFlush((Display*)fDisplay);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy area of current window in the pixmap pix.

void TGX11::CopyWindowtoPixmap(Drawable *pix, int xpos, int ypos )
{
   Window root;
   int xx, yy;
   unsigned int ww, hh, border, depth;

   XGetGeometry((Display*)fDisplay, *pix, &root, &xx, &yy, &ww, &hh, &border, &depth);
   XCopyArea((Display*)fDisplay, gCws->fDrawing, *pix, *gGCpxmp, xpos, ypos, ww, hh, 0, 0);
   XFlush((Display*)fDisplay);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a box.
///
///  - mode=0 hollow  (kHollow)
///  - mode=1 solid   (kSolid)

void TGX11::DrawBox(int x1, int y1, int x2, int y2, EBoxMode mode)
{
   Int_t x = TMath::Min(x1, x2);
   Int_t y = TMath::Min(y1, y2);
   Int_t w = TMath::Abs(x2 - x1);
   Int_t h = TMath::Abs(y2 - y1);

   switch (mode) {

      case kHollow:
         XDrawRectangle((Display*)fDisplay, gCws->fDrawing, *gGCline, x, y, w, h);
         break;

      case kFilled:
         XFillRectangle((Display*)fDisplay, gCws->fDrawing, *gGCfill, x, y, w, h);
         break;

      default:
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a cell array.
//
///  \param [in] x1,y1        : left down corner
///  \param [in] x2,y2        : right up corner
///  \param [in] nx,ny        : array size
///  \param [in] ic           : array
///
/// Draw a cell array. The drawing is done with the pixel precision
/// if (X2-X1)/NX (or Y) is not a exact pixel number the position of
/// the top right corner may be wrong.

void TGX11::DrawCellArray(int x1, int y1, int x2, int y2, int nx, int ny, int *ic)
{
   int i, j, icol, ix, iy, w, h, current_icol;

   current_icol = -1;
   w            = TMath::Max((x2-x1)/(nx),1);
   h            = TMath::Max((y1-y2)/(ny),1);
   ix           = x1;

   for (i = 0; i < nx; i++) {
      iy = y1-h;
      for (j = 0; j < ny; j++) {
         icol = ic[i+(nx*j)];
         if (icol != current_icol) {
            XSetForeground((Display*)fDisplay, *gGCfill, GetColor(icol).fPixel);
            current_icol = icol;
         }
         XFillRectangle((Display*)fDisplay, gCws->fDrawing, *gGCfill, ix, iy, w, h);
         iy = iy-h;
      }
      ix = ix+w;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill area described by polygon.
///
///  \param [in] n     number of points
///  \param [in] xyt   list of points

void TGX11::DrawFillArea(int n, TPoint *xyt)
{
   XPoint *xy = (XPoint*)xyt;

   if (gFillHollow)
      XDrawLines((Display*)fDisplay, gCws->fDrawing, *gGCfill, xy, n, CoordModeOrigin);

   else {
      XFillPolygon((Display*)fDisplay, gCws->fDrawing, *gGCfill,
                   xy, n, Nonconvex, CoordModeOrigin);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a line.
///
///  \param [in] x1,y1        : begin of line
///  \param [in] x2,y2        : end of line

void TGX11::DrawLine(int x1, int y1, int x2, int y2)
{
   if (gLineStyle == LineSolid)
      XDrawLine((Display*)fDisplay, gCws->fDrawing, *gGCline, x1, y1, x2, y2);
   else {
      XSetDashes((Display*)fDisplay, *gGCdash, gDashOffset, gDashList, gDashSize);
      XDrawLine((Display*)fDisplay, gCws->fDrawing, *gGCdash, x1, y1, x2, y2);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a line through all points.
///
///  \param [in] n     number of points
///  \param [in] xyt   list of points

void TGX11::DrawPolyLine(int n, TPoint *xyt)
{
   XPoint *xy = (XPoint*)xyt;

   const Int_t kMaxPoints = 1000001;

   if (n > kMaxPoints) {
      int ibeg = 0;
      int iend = kMaxPoints - 1;
      while (iend < n) {
         DrawPolyLine( kMaxPoints, &xyt[ibeg] );
         ibeg = iend;
         iend += kMaxPoints - 1;
      }
      if (ibeg < n) {
         int npt = n - ibeg;
         DrawPolyLine( npt, &xyt[ibeg] );
      }
   } else if (n > 1) {
      if (gLineStyle == LineSolid)
         XDrawLines((Display*)fDisplay, gCws->fDrawing, *gGCline, xy, n, CoordModeOrigin);
      else {
         int i;
         XSetDashes((Display*)fDisplay, *gGCdash,
                    gDashOffset, gDashList, gDashSize);
         XDrawLines((Display*)fDisplay, gCws->fDrawing, *gGCdash, xy, n, CoordModeOrigin);

         // calculate length of line to update dash offset
         for (i = 1; i < n; i++) {
            int dx = xy[i].x - xy[i-1].x;
            int dy = xy[i].y - xy[i-1].y;
            if (dx < 0) dx = - dx;
            if (dy < 0) dy = - dy;
            gDashOffset += dx > dy ? dx : dy;
         }
         gDashOffset %= gDashLength;
      }
   } else {
      int px,py;
      px=xy[0].x;
      py=xy[0].y;
      XDrawPoint((Display*)fDisplay, gCws->fDrawing,
                 gLineStyle == LineSolid ? *gGCline : *gGCdash, px, py);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw n markers with the current attributes at position x, y.
///
///  \param [in] n     number of markers to draw
///  \param [in] xyt   x,y coordinates of markers

void TGX11::DrawPolyMarker(int n, TPoint *xyt)
{
   XPoint *xy = (XPoint*)xyt;

   if (gMarker.n <= 0) {
      const int kNMAX = 1000000;
      int nt = n/kNMAX;
      for (int it=0;it<=nt;it++) {
         if (it < nt) {
            XDrawPoints((Display*)fDisplay, gCws->fDrawing, *gGCmark, &xy[it*kNMAX], kNMAX, CoordModeOrigin);
         } else {
            XDrawPoints((Display*)fDisplay, gCws->fDrawing, *gGCmark, &xy[it*kNMAX], n-it*kNMAX, CoordModeOrigin);
         }
      }
   } else {
      int r = gMarker.n / 2;
      int m;

      for (m = 0; m < n; m++) {
         int hollow = 0;

         switch (gMarker.type) {
            int i;

            case 0:        // hollow circle
               XDrawArc((Display*)fDisplay, gCws->fDrawing, *gGCmark,
                        xy[m].x - r, xy[m].y - r, gMarker.n, gMarker.n, 0, 360*64);
               break;

            case 1:        // filled circle
               XFillArc((Display*)fDisplay, gCws->fDrawing, *gGCmark,
                        xy[m].x - r, xy[m].y - r, gMarker.n, gMarker.n, 0, 360*64);
               break;

            case 2:        // hollow polygon
               hollow = 1;
            case 3:        // filled polygon
               for (i = 0; i < gMarker.n; i++) {
                  gMarker.xy[i].x += xy[m].x;
                  gMarker.xy[i].y += xy[m].y;
               }
               if (hollow)
                  XDrawLines((Display*)fDisplay, gCws->fDrawing, *gGCmark,
                             gMarker.xy, gMarker.n, CoordModeOrigin);
               else
                  XFillPolygon((Display*)fDisplay, gCws->fDrawing, *gGCmark,
                               gMarker.xy, gMarker.n, Nonconvex, CoordModeOrigin);
               for (i = 0; i < gMarker.n; i++) {
                  gMarker.xy[i].x -= xy[m].x;
                  gMarker.xy[i].y -= xy[m].y;
               }
               break;

            case 4:        // segmented line
               for (i = 0; i < gMarker.n; i += 2)
                  XDrawLine((Display*)fDisplay, gCws->fDrawing, *gGCmark,
                            xy[m].x + gMarker.xy[i].x, xy[m].y + gMarker.xy[i].y,
                            xy[m].x + gMarker.xy[i+1].x, xy[m].y + gMarker.xy[i+1].y);
               break;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a text string using current font.
///
///  \param [in] mode       : drawing mode
///              - mode=0     : the background is not drawn (kClear)
///              - mode=1     : the background is drawn (kOpaque)
///  \param [in] x,y        : text position
///  \param [in] angle      : text angle
///  \param [in] mgn        : magnification factor
///  \param [in] text       : text string

void TGX11::DrawText(int x, int y, float angle, float mgn,
                     const char *text, ETextMode mode)
{
   XRotSetMagnification(mgn);

   if (!text) return;

   switch (mode) {

      case kClear:
         XRotDrawAlignedString((Display*)fDisplay, gTextFont, angle,
                      gCws->fDrawing, *gGCtext, x, y, (char*)text, fTextAlign);
         break;

      case kOpaque:
         XRotDrawAlignedImageString((Display*)fDisplay, gTextFont, angle,
                      gCws->fDrawing, *gGCtext, x, y, (char*)text, fTextAlign);
         break;

      default:
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find best visual, i.e. the one with the most planes and TrueColor or
/// DirectColor. Sets fVisual, fDepth, fRootWin, fColormap, fBlackPixel
/// and fWhitePixel.

void TGX11::FindBestVisual()
{
   Int_t findvis = gEnv->GetValue("X11.FindBestVisual", 1);

   Visual *vis = DefaultVisual((Display*)fDisplay, fScreenNumber);
   if (((vis->c_class != TrueColor && vis->c_class != DirectColor) ||
       DefaultDepth((Display*)fDisplay, fScreenNumber) < 15) && findvis) {

      // try to find better visual
      static XVisualInfo templates[] = {
         // Visual, visualid, screen, depth, class      , red_mask, green_mask, blue_mask, colormap_size, bits_per_rgb
         { 0     , 0       , 0     , 24   , TrueColor  , 0       , 0         , 0        , 0            , 0 },
         { 0     , 0       , 0     , 32   , TrueColor  , 0       , 0         , 0        , 0            , 0 },
         { 0     , 0       , 0     , 16   , TrueColor  , 0       , 0         , 0        , 0            , 0 },
         { 0     , 0       , 0     , 15   , TrueColor  , 0       , 0         , 0        , 0            , 0 },
         // no suitable TrueColorMode found - now do the same thing to DirectColor
         { 0     , 0       , 0     , 24   , DirectColor, 0       , 0         , 0        , 0            , 0 },
         { 0     , 0       , 0     , 32   , DirectColor, 0       , 0         , 0        , 0            , 0 },
         { 0     , 0       , 0     , 16   , DirectColor, 0       , 0         , 0        , 0            , 0 },
         { 0     , 0       , 0     , 15   , DirectColor, 0       , 0         , 0        , 0            , 0 },
         { 0     , 0       , 0     , 0    , 0          , 0       , 0         , 0        , 0            , 0 },
      };

      Int_t nitems = 0;
      XVisualInfo *vlist = 0;
      for (Int_t i = 0; templates[i].depth != 0; i++) {
         Int_t mask = VisualScreenMask|VisualDepthMask|VisualClassMask;
         templates[i].screen = fScreenNumber;
         if ((vlist = XGetVisualInfo((Display*)fDisplay, mask, &(templates[i]), &nitems))) {
            FindUsableVisual((RXVisualInfo*)vlist, nitems);
            XFree(vlist);
            vlist = 0;
            if (fVisual)
               break;
         }
      }
   }

   fRootWin = RootWindow((Display*)fDisplay, fScreenNumber);

   if (!fVisual) {
      fDepth      = DefaultDepth((Display*)fDisplay, fScreenNumber);
      fVisual     = (RVisual*)DefaultVisual((Display*)fDisplay, fScreenNumber);
      fVisRootWin = fRootWin;
      if (fDepth > 1)
         fColormap = DefaultColormap((Display*)fDisplay, fScreenNumber);
      fBlackPixel = BlackPixel((Display*)fDisplay, fScreenNumber);
      fWhitePixel = WhitePixel((Display*)fDisplay, fScreenNumber);
   }
   if (gDebug > 1)
      Printf("Selected visual 0x%lx: depth %d, class %d, colormap: %s",
             fVisual->visualid, fDepth, fVisual->c_class,
             fColormap == DefaultColormap((Display*)fDisplay, fScreenNumber) ? "default" :
             "custom");
}

////////////////////////////////////////////////////////////////////////////////
/// Dummy error handler for X11. Used by FindUsableVisual().

static Int_t DummyX11ErrorHandler(Display *, XErrorEvent *)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if visual is usable, if so set fVisual, fDepth, fColormap,
/// fBlackPixel and fWhitePixel.

void TGX11::FindUsableVisual(RXVisualInfo *vlist, Int_t nitems)
{
   Int_t (*oldErrorHandler)(Display *, XErrorEvent *) =
       XSetErrorHandler(DummyX11ErrorHandler);

   XSetWindowAttributes attr;
   memset(&attr, 0, sizeof(attr));

   Window root = RootWindow((Display*)fDisplay, fScreenNumber);

   for (Int_t i = 0; i < nitems; i++) {
      Window w = None, wjunk;
      UInt_t width, height, ujunk;
      Int_t  junk;

      // try and use default colormap when possible
      if (vlist[i].visual == DefaultVisual((Display*)fDisplay, fScreenNumber)) {
         attr.colormap = DefaultColormap((Display*)fDisplay, fScreenNumber);
      } else {
         attr.colormap = XCreateColormap((Display*)fDisplay, root, vlist[i].visual, AllocNone);
      }

      static XColor black_xcol = { 0, 0x0000, 0x0000, 0x0000, DoRed|DoGreen|DoBlue, 0 };
      static XColor white_xcol = { 0, 0xFFFF, 0xFFFF, 0xFFFF, DoRed|DoGreen|DoBlue, 0 };
      XAllocColor((Display*)fDisplay, attr.colormap, &black_xcol);
      XAllocColor((Display*)fDisplay, attr.colormap, &white_xcol);
      attr.border_pixel = black_xcol.pixel;
      attr.override_redirect = True;

      w = XCreateWindow((Display*)fDisplay, root, -20, -20, 10, 10, 0, vlist[i].depth,
                        CopyFromParent, vlist[i].visual,
                        CWColormap|CWBorderPixel|CWOverrideRedirect, &attr);
      if (w != None && XGetGeometry((Display*)fDisplay, w, &wjunk, &junk, &junk,
                                    &width, &height, &ujunk, &ujunk)) {
         fVisual     = (RVisual*)vlist[i].visual;
         fDepth      = vlist[i].depth;
         fColormap   = attr.colormap;
         fBlackPixel = black_xcol.pixel;
         fWhitePixel = white_xcol.pixel;
         fVisRootWin = w;
         break;
      }
      if (attr.colormap != DefaultColormap((Display*)fDisplay, fScreenNumber))
         XFreeColormap((Display*)fDisplay, attr.colormap);
   }
   XSetErrorHandler(oldErrorHandler);
}

////////////////////////////////////////////////////////////////////////////////
/// Return character up vector.

void TGX11::GetCharacterUp(Float_t &chupx, Float_t &chupy)
{
   chupx = fCharacterUpX;
   chupy = fCharacterUpY;
}

////////////////////////////////////////////////////////////////////////////////
/// Return reference to internal color structure associated
/// to color index cid.

XColor_t &TGX11::GetColor(Int_t cid)
{
   XColor_t *col = (XColor_t*) (Long_t)fColors->GetValue(cid);
   if (!col) {
      col = new XColor_t;
      fColors->Add(cid, (Long_t) col);
   }
   return *col;
}

////////////////////////////////////////////////////////////////////////////////
/// Return current window pointer. Protected method used by TGX11TTF.

Window_t TGX11::GetCurrentWindow() const
{
   return (Window_t)(gCws ? gCws->fDrawing : 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Return desired Graphics Context ("which" maps directly on gGCList[]).
/// Protected method used by TGX11TTF.

void *TGX11::GetGC(Int_t which) const
{
   if (which >= kMAXGC || which < 0) {
      Error("GetGC", "trying to get illegal GC (which = %d)", which);
      return 0;
   }
   return &gGClist[which];
}

////////////////////////////////////////////////////////////////////////////////
/// Query the double buffer value for the window wid.

Int_t TGX11::GetDoubleBuffer(int wid)
{
   gTws = &fWindows[wid];
   if (!gTws->fOpen)
      return -1;
   else
      return gTws->fDoubleBuffer;
}

////////////////////////////////////////////////////////////////////////////////
/// Return position and size of window wid.
///
///  \param [in] wid        : window identifier
///  \param [in] x,y        : window position (output)
///  \param [in] w,h        : window size (output)
///
/// if wid < 0 the size of the display is returned

void TGX11::GetGeometry(int wid, int &x, int &y, unsigned int &w, unsigned int &h)
{
   Window junkwin=0;

   if (wid < 0) {
      x = 0;
      y = 0;
      w = DisplayWidth((Display*)fDisplay,fScreenNumber);
      h = DisplayHeight((Display*)fDisplay,fScreenNumber);
   } else {
      Window root;
      unsigned int border, depth;
      unsigned int width, height;

      gTws = &fWindows[wid];
      XGetGeometry((Display*)fDisplay, gTws->fWindow, &root, &x, &y,
                   &width, &height, &border, &depth);
      XTranslateCoordinates((Display*)fDisplay, gTws->fWindow, fRootWin,
                            0, 0, &x, &y, &junkwin);
      if (width >= 65535)
         width = 1;
      if (height >= 65535)
         height = 1;
      if (width > 0 && height > 0) {
         gTws->fWidth  = width;
         gTws->fHeight = height;
      }
      w = gTws->fWidth;
      h = gTws->fHeight;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return hostname on which the display is opened.

const char *TGX11::DisplayName(const char *dpyName)
{
   return XDisplayName(dpyName);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pixel value associated to specified ROOT color number.

ULong_t TGX11::GetPixel(Color_t ci)
{
   TColor *color = gROOT->GetColor(ci);
   if (color)
      SetRGB(ci, color->GetRed(), color->GetGreen(), color->GetBlue());
//   else
//      Warning("GetPixel", "color with index %d not defined", ci);

   XColor_t &col = GetColor(ci);
   return col.fPixel;
}

////////////////////////////////////////////////////////////////////////////////
/// Get maximum number of planes.

void TGX11::GetPlanes(int &nplanes)
{
   nplanes = fDepth;
}

////////////////////////////////////////////////////////////////////////////////
/// Get rgb values for color "index".

void TGX11::GetRGB(int index, float &r, float &g, float &b)
{
   if (index == 0) {
      r = g = b = 1.0;
   } else if (index == 1) {
      r = g = b = 0.0;
   } else {
      XColor_t &col = GetColor(index);
      r = ((float) col.fRed) / ((float) kBIGGEST_RGB_VALUE);
      g = ((float) col.fGreen) / ((float) kBIGGEST_RGB_VALUE);
      b = ((float) col.fBlue) / ((float) kBIGGEST_RGB_VALUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the size of a character string.
///
///  \param [in] w          : text width
///  \param [in] h          : text height
///  \param [in] mess        : message

void TGX11::GetTextExtent(unsigned int &w, unsigned int &h, char *mess)
{
   w=0; h=0;
   if (strlen(mess)==0) return;

   XPoint *cBox;
   XRotSetMagnification(fTextMagnitude);
   cBox = XRotTextExtents((Display*)fDisplay, gTextFont, 0., 0, 0, mess, 0);
   if (cBox) {
      w    = cBox[2].x;
      h    = -cBox[2].y;
      free((char *)cBox);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the X11 window identifier.
///
///  \param [in] wid      : Workstation identifier (input)

Window_t TGX11::GetWindowID(int wid)
{
   return (Window_t) fWindows[wid].fWindow;
}

////////////////////////////////////////////////////////////////////////////////
/// Move the window wid.
///
///  \param [in] wid  : Window identifier.
///  \param [in] x    : x new window position
///  \param [in] y    : y new window position

void TGX11::MoveWindow(int wid, int x, int y)
{
   gTws = &fWindows[wid];
   if (!gTws->fOpen) return;

   XMoveWindow((Display*)fDisplay, gTws->fWindow, x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Open the display. Return -1 if the opening fails, 0 when ok.

Int_t TGX11::OpenDisplay(void *disp)
{
   Pixmap  pixmp1, pixmp2;
   XColor  fore, back;
   char  **fontlist;
   int     fontcount = 0;
   int     i;

   if (fDisplay) return 0;

   fDisplay      = disp;
   fScreenNumber = DefaultScreen((Display*)fDisplay);

   FindBestVisual();

   GetColor(1).fDefined = kTRUE; // default foreground
   GetColor(1).fPixel = fBlackPixel;
   GetColor(0).fDefined = kTRUE; // default background
   GetColor(0).fPixel = fWhitePixel;

   // Inquire the the XServer Vendor
   char vendor[132];
   strlcpy(vendor, XServerVendor((Display*)fDisplay),132);

   // Create primitives graphic contexts
   for (i = 0; i < kMAXGC; i++)
      gGClist[i] = XCreateGC((Display*)fDisplay, fVisRootWin, 0, 0);

   XGCValues values;
   if (XGetGCValues((Display*)fDisplay, *gGCtext, GCForeground|GCBackground, &values)) {
      XSetForeground((Display*)fDisplay, *gGCinvt, values.background);
      XSetBackground((Display*)fDisplay, *gGCinvt, values.foreground);
   } else {
      Error("OpenDisplay", "cannot get GC values");
   }

   // Turn-off GraphicsExpose and NoExpose event reporting for the pixmap
   // manipulation GC, this to prevent these events from being stacked up
   // without ever being processed and thereby wasting a lot of memory.
   XSetGraphicsExposures((Display*)fDisplay, *gGCpxmp, False);

   // Create input echo graphic context
   XGCValues echov;
   echov.foreground = fBlackPixel;
   echov.background = fWhitePixel;
   if (strstr(vendor,"Hewlett"))
      echov.function   = GXxor;
   else
      echov.function   = GXinvert;

   gGCecho = XCreateGC((Display*)fDisplay, fVisRootWin,
                       GCForeground | GCBackground | GCFunction,
                       &echov);

   // Load a default Font
   static int isdisp = 0;
   if (!isdisp) {
      for (i = 0; i < kMAXFONT; i++) {
         gFont[i].id = 0;
         strcpy(gFont[i].name, " ");
      }
      fontlist = XListFonts((Display*)fDisplay, "*courier*", 1, &fontcount);
      if (fontlist && fontcount != 0) {
         gFont[gCurrentFontNumber].id = XLoadQueryFont((Display*)fDisplay, fontlist[0]);
         gTextFont = gFont[gCurrentFontNumber].id;
         strcpy(gFont[gCurrentFontNumber].name, "*courier*");
         gCurrentFontNumber++;
         XFreeFontNames(fontlist);
      } else {
         // emergency: try fixed font
         fontlist = XListFonts((Display*)fDisplay, "fixed", 1, &fontcount);
         if (fontlist && fontcount != 0) {
            gFont[gCurrentFontNumber].id = XLoadQueryFont((Display*)fDisplay, fontlist[0]);
            gTextFont = gFont[gCurrentFontNumber].id;
            strcpy(gFont[gCurrentFontNumber].name, "fixed");
            gCurrentFontNumber++;
            XFreeFontNames(fontlist);
         } else {
            Warning("OpenDisplay", "no default font loaded");
         }
      }
      isdisp = 1;
   }

   // Create a null cursor
   pixmp1 = XCreateBitmapFromData((Display*)fDisplay, fRootWin,
                                  null_cursor_bits, 16, 16);
   pixmp2 = XCreateBitmapFromData((Display*)fDisplay, fRootWin,
                                  null_cursor_bits, 16, 16);
   gNullCursor = XCreatePixmapCursor((Display*)fDisplay,pixmp1,pixmp2,&fore,&back,0,0);

   // Create cursors
   fCursors[kBottomLeft]  = XCreateFontCursor((Display*)fDisplay, XC_bottom_left_corner);
   fCursors[kBottomRight] = XCreateFontCursor((Display*)fDisplay, XC_bottom_right_corner);
   fCursors[kTopLeft]     = XCreateFontCursor((Display*)fDisplay, XC_top_left_corner);
   fCursors[kTopRight]    = XCreateFontCursor((Display*)fDisplay, XC_top_right_corner);
   fCursors[kBottomSide]  = XCreateFontCursor((Display*)fDisplay, XC_bottom_side);
   fCursors[kLeftSide]    = XCreateFontCursor((Display*)fDisplay, XC_left_side);
   fCursors[kTopSide]     = XCreateFontCursor((Display*)fDisplay, XC_top_side);
   fCursors[kRightSide]   = XCreateFontCursor((Display*)fDisplay, XC_right_side);
   fCursors[kMove]        = XCreateFontCursor((Display*)fDisplay, XC_fleur);
   fCursors[kCross]       = XCreateFontCursor((Display*)fDisplay, XC_tcross);
   fCursors[kArrowHor]    = XCreateFontCursor((Display*)fDisplay, XC_sb_h_double_arrow);
   fCursors[kArrowVer]    = XCreateFontCursor((Display*)fDisplay, XC_sb_v_double_arrow);
   fCursors[kHand]        = XCreateFontCursor((Display*)fDisplay, XC_hand2);
   fCursors[kRotate]      = XCreateFontCursor((Display*)fDisplay, XC_exchange);
   fCursors[kPointer]     = XCreateFontCursor((Display*)fDisplay, XC_left_ptr);
   fCursors[kArrowRight]  = XCreateFontCursor((Display*)fDisplay, XC_arrow);
   fCursors[kCaret]       = XCreateFontCursor((Display*)fDisplay, XC_xterm);
   fCursors[kWatch]       = XCreateFontCursor((Display*)fDisplay, XC_watch);
   fCursors[kNoDrop]      = XCreateFontCursor((Display*)fDisplay, XC_pirate);

   // Setup color information
   fRedDiv = fGreenDiv = fBlueDiv = fRedShift = fGreenShift = fBlueShift = -1;

   if (fVisual->c_class == TrueColor) {
      for (i = 0; i < int(sizeof(fVisual->blue_mask)*kBitsPerByte); i++) {
         if (fBlueShift == -1 && ((fVisual->blue_mask >> i) & 1))
            fBlueShift = i;
         if ((fVisual->blue_mask >> i) == 1) {
            fBlueDiv = sizeof(UShort_t)*kBitsPerByte - i - 1 + fBlueShift;
            break;
         }
      }
      for (i = 0; i < int(sizeof(fVisual->green_mask)*kBitsPerByte); i++) {
         if (fGreenShift == -1 && ((fVisual->green_mask >> i) & 1))
            fGreenShift = i;
         if ((fVisual->green_mask >> i) == 1) {
            fGreenDiv = sizeof(UShort_t)*kBitsPerByte - i - 1 + fGreenShift;
            break;
         }
      }
      for (i = 0; i < int(sizeof(fVisual->red_mask)*kBitsPerByte); i++) {
         if (fRedShift == -1 && ((fVisual->red_mask >> i) & 1))
            fRedShift = i;
         if ((fVisual->red_mask >> i) == 1) {
            fRedDiv = sizeof(UShort_t)*kBitsPerByte - i - 1 + fRedShift;
            break;
         }
      }
      //printf("fRedDiv = %d, fGreenDiv = %d, fBlueDiv = %d, fRedShift = %d, fGreenShift = %d, fBlueShift = %d\n",
      //       fRedDiv, fGreenDiv, fBlueDiv, fRedShift, fGreenShift, fBlueShift);
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a new pixmap.
///
///  \param [in] w,h : Width and height of the pixmap.

Int_t TGX11::OpenPixmap(unsigned int w, unsigned int h)
{
   Window root;
   unsigned int wval, hval;
   int xx, yy, i, wid;
   unsigned int ww, hh, border, depth;
   wval = w;
   hval = h;

   // Select next free window number

again:
   for (wid = 0; wid < fMaxNumberOfWindows; wid++)
      if (!fWindows[wid].fOpen) {
         fWindows[wid].fOpen = 1;
         gCws = &fWindows[wid];
         break;
      }

   if (wid == fMaxNumberOfWindows) {
      int newsize = fMaxNumberOfWindows + 10;
      fWindows = (XWindow_t*) TStorage::ReAlloc(fWindows, newsize*sizeof(XWindow_t),
                                                fMaxNumberOfWindows*sizeof(XWindow_t));
      for (i = fMaxNumberOfWindows; i < newsize; i++)
         fWindows[i].fOpen = 0;
      fMaxNumberOfWindows = newsize;
      goto again;
   }

   gCws->fWindow = XCreatePixmap((Display*)fDisplay, fRootWin, wval, hval, fDepth);
   XGetGeometry((Display*)fDisplay, gCws->fWindow, &root, &xx, &yy, &ww, &hh, &border, &depth);

   for (i = 0; i < kMAXGC; i++)
      XSetClipMask((Display*)fDisplay, gGClist[i], None);

   SetColor(gGCpxmp, 0);
   XFillRectangle((Display*)fDisplay, gCws->fWindow, *gGCpxmp, 0, 0, ww, hh);
   SetColor(gGCpxmp, 1);

   // Initialise the window structure
   gCws->fDrawing       = gCws->fWindow;
   gCws->fBuffer        = 0;
   gCws->fDoubleBuffer  = 0;
   gCws->fIsPixmap      = 1;
   gCws->fClip          = 0;
   gCws->fWidth         = wval;
   gCws->fHeight        = hval;
   gCws->fNewColors     = 0;
   gCws->fShared        = kFALSE;

   return wid;
}

////////////////////////////////////////////////////////////////////////////////
/// Open window and return window number.
///
/// \return -1 if window initialization fails.

Int_t TGX11::InitWindow(ULong_t win)
{
   XSetWindowAttributes attributes;
   ULong_t attr_mask = 0;
   int wid;
   int xval, yval;
   unsigned int wval, hval, border, depth;
   Window root;

   Window wind = (Window) win;

   XGetGeometry((Display*)fDisplay, wind, &root, &xval, &yval, &wval, &hval, &border, &depth);

   // Select next free window number

again:
   for (wid = 0; wid < fMaxNumberOfWindows; wid++)
      if (!fWindows[wid].fOpen) {
         fWindows[wid].fOpen = 1;
         fWindows[wid].fDoubleBuffer = 0;
         gCws = &fWindows[wid];
         break;
      }

   if (wid == fMaxNumberOfWindows) {
      int newsize = fMaxNumberOfWindows + 10;
      fWindows = (XWindow_t*) TStorage::ReAlloc(fWindows, newsize*sizeof(XWindow_t),
                                                fMaxNumberOfWindows*sizeof(XWindow_t));
      for (int i = fMaxNumberOfWindows; i < newsize; i++)
         fWindows[i].fOpen = 0;
      fMaxNumberOfWindows = newsize;
      goto again;
   }

   // Create window

   attributes.background_pixel = GetColor(0).fPixel;
   attr_mask |= CWBackPixel;
   attributes.border_pixel = GetColor(1).fPixel;
   attr_mask |= CWBorderPixel;
   attributes.event_mask = NoEventMask;
   attr_mask |= CWEventMask;
   attributes.backing_store = Always;
   attr_mask |= CWBackingStore;
   attributes.bit_gravity = NorthWestGravity;
   attr_mask |= CWBitGravity;
   if (fColormap) {
      attributes.colormap = fColormap;
      attr_mask |= CWColormap;
   }

   gCws->fWindow = XCreateWindow((Display*)fDisplay, wind,
                                 xval, yval, wval, hval, 0, fDepth,
                                 InputOutput, fVisual,
                                 attr_mask, &attributes);

   XMapWindow((Display*)fDisplay, gCws->fWindow);
   XFlush((Display*)fDisplay);

   // Initialise the window structure

   gCws->fDrawing      = gCws->fWindow;
   gCws->fBuffer       = 0;
   gCws->fDoubleBuffer = 0;
   gCws->fIsPixmap     = 0;
   gCws->fClip         = 0;
   gCws->fWidth        = wval;
   gCws->fHeight       = hval;
   gCws->fNewColors    = 0;
   gCws->fShared       = kFALSE;

   return wid;
}

////////////////////////////////////////////////////////////////////////////////
/// Register a window created by Qt as a ROOT window (like InitWindow()).

Int_t TGX11::AddWindow(ULong_t qwid, UInt_t w, UInt_t h)
{
   Int_t wid;

   // Select next free window number

again:
   for (wid = 0; wid < fMaxNumberOfWindows; wid++)
      if (!fWindows[wid].fOpen) {
         fWindows[wid].fOpen = 1;
         fWindows[wid].fDoubleBuffer = 0;
         gCws = &fWindows[wid];
         break;
      }

   if (wid == fMaxNumberOfWindows) {
      int newsize = fMaxNumberOfWindows + 10;
      fWindows = (XWindow_t*) TStorage::ReAlloc(fWindows, newsize*sizeof(XWindow_t),
                                                fMaxNumberOfWindows*sizeof(XWindow_t));
      for (int i = fMaxNumberOfWindows; i < newsize; i++)
         fWindows[i].fOpen = 0;
      fMaxNumberOfWindows = newsize;
      goto again;
   }

   gCws->fWindow = qwid;

   //init Xwindow_t struct
   gCws->fDrawing       = gCws->fWindow;
   gCws->fBuffer        = 0;
   gCws->fDoubleBuffer  = 0;
   gCws->fIsPixmap      = 0;
   gCws->fClip          = 0;
   gCws->fWidth         = w;
   gCws->fHeight        = h;
   gCws->fNewColors     = 0;
   gCws->fShared        = kTRUE;

   return wid;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a window created by Qt (like CloseWindow1()).

void TGX11::RemoveWindow(ULong_t qwid)
{
   SelectWindow((int)qwid);

   if (gCws->fBuffer) XFreePixmap((Display*)fDisplay, gCws->fBuffer);

   if (gCws->fNewColors) {
      if (fRedDiv == -1)
         XFreeColors((Display*)fDisplay, fColormap, gCws->fNewColors, gCws->fNcolors, 0);
      delete [] gCws->fNewColors;
      gCws->fNewColors = 0;
   }

   gCws->fOpen = 0;

   // make first window in list the current window
   for (Int_t wid = 0; wid < fMaxNumberOfWindows; wid++)
      if (fWindows[wid].fOpen) {
         gCws = &fWindows[wid];
         return;
      }

   gCws = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Query pointer position.
///
///  \param [in] ix       : X coordinate of pointer
///  \param [in] iy       : Y coordinate of pointer
/// (both coordinates are relative to the origin of the root window)

void TGX11::QueryPointer(int &ix, int &iy)
{
   Window    root_return, child_return;
   int       win_x_return, win_y_return;
   int       root_x_return, root_y_return;
   unsigned int mask_return;

   XQueryPointer((Display*)fDisplay,gCws->fWindow, &root_return,
                 &child_return, &root_x_return, &root_y_return, &win_x_return,
                 &win_y_return, &mask_return);

   ix = root_x_return;
   iy = root_y_return;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the pixmap pix.

void  TGX11::RemovePixmap(Drawable *pix)
{
   XFreePixmap((Display*)fDisplay,*pix);
}

////////////////////////////////////////////////////////////////////////////////
/// Request Locator position.
///
///  \param [in] x,y       : cursor position at moment of button press (output)
///  \param [in] ctyp      : cursor type (input)
///              - ctyp=1 tracking cross
///              - ctyp=2 cross-hair
///              - ctyp=3 rubber circle
///              - ctyp=4 rubber band
///              - ctyp=5 rubber rectangle
///
///  \param [in] mode      : input mode
///              - mode=0 request
///              - mode=1 sample
///
/// Request locator:
/// return button number:
///                     - 1 = left is pressed
///                     - 2 = middle is pressed
///                     - 3 = right is pressed
///        in sample mode:
///                     - 11 = left is released
///                     - 12 = middle is released
///                     - 13 = right is released
///                     - -1 = nothing is pressed or released
///                     - -2 = leave the window
///                     - else = keycode (keyboard is pressed)

Int_t TGX11::RequestLocator(int mode, int ctyp, int &x, int &y)
{
   static int xloc  = 0;
   static int yloc  = 0;
   static int xlocp = 0;
   static int ylocp = 0;
   static Cursor cursor = 0;

   XEvent event;
   int button_press;
   int radius;

   // Change the cursor shape
   if (cursor == 0) {
      if (ctyp > 1) {
         XDefineCursor((Display*)fDisplay, gCws->fWindow, gNullCursor);
         XSetForeground((Display*)fDisplay, gGCecho, GetColor(0).fPixel);
      } else {
         cursor = XCreateFontCursor((Display*)fDisplay, XC_crosshair);
         XDefineCursor((Display*)fDisplay, gCws->fWindow, cursor);
      }
   }

   // Event loop

   button_press = 0;

   while (button_press == 0) {

      switch (ctyp) {

         case 1 :
            break;

         case 2 :
            XDrawLine((Display*)fDisplay, gCws->fWindow, gGCecho,
                      xloc, 0, xloc, gCws->fHeight);
            XDrawLine((Display*)fDisplay, gCws->fWindow, gGCecho,
                      0, yloc, gCws->fWidth, yloc);
            break;

         case 3 :
            radius = (int) TMath::Sqrt((double)((xloc-xlocp)*(xloc-xlocp) +
                                       (yloc-ylocp)*(yloc-ylocp)));
            XDrawArc((Display*)fDisplay, gCws->fWindow, gGCecho,
                     xlocp-radius, ylocp-radius,
                     2*radius, 2*radius, 0, 23040);
            break;

         case 4 :
            XDrawLine((Display*)fDisplay, gCws->fWindow, gGCecho,
                      xlocp, ylocp, xloc, yloc);
            break;

         case 5 :
            XDrawRectangle((Display*)fDisplay, gCws->fWindow, gGCecho,
                           TMath::Min(xlocp,xloc), TMath::Min(ylocp,yloc),
                           TMath::Abs(xloc-xlocp), TMath::Abs(yloc-ylocp));
            break;

         default:
            break;
      }

      while (XEventsQueued( (Display*)fDisplay, QueuedAlready) > 1) {
         XNextEvent((Display*)fDisplay, &event);
      }
      XWindowEvent((Display*)fDisplay, gCws->fWindow, gMouseMask, &event);

      switch (ctyp) {

         case 1 :
            break;

         case 2 :
            XDrawLine((Display*)fDisplay, gCws->fWindow, gGCecho,
                      xloc, 0, xloc, gCws->fHeight);
            XDrawLine((Display*)fDisplay, gCws->fWindow, gGCecho,
                      0, yloc, gCws->fWidth, yloc);
            break;

         case 3 :
            radius = (int) TMath::Sqrt((double)((xloc-xlocp)*(xloc-xlocp) +
                                           (yloc-ylocp)*(yloc-ylocp)));
            XDrawArc((Display*)fDisplay, gCws->fWindow, gGCecho,
                     xlocp-radius, ylocp-radius,
                     2*radius, 2*radius, 0, 23040);
            break;

         case 4 :
            XDrawLine((Display*)fDisplay, gCws->fWindow, gGCecho,
                      xlocp, ylocp, xloc, yloc);
            break;

         case 5 :
            XDrawRectangle((Display*)fDisplay, gCws->fWindow, gGCecho,
                           TMath::Min(xlocp,xloc), TMath::Min(ylocp,yloc),
                           TMath::Abs(xloc-xlocp), TMath::Abs(yloc-ylocp));
            break;

         default:
            break;
      }

      xloc = event.xbutton.x;
      yloc = event.xbutton.y;

      switch (event.type) {

         case LeaveNotify :
            if (mode == 0) {
               while (1) {
                  XNextEvent((Display*)fDisplay, &event);
                  if (event.type == EnterNotify) break;
               }
            } else {
               button_press = -2;
            }
            break;

         case ButtonPress :
            button_press = event.xbutton.button ;
            xlocp = event.xbutton.x;
            ylocp = event.xbutton.y;
            XUndefineCursor( (Display*)fDisplay, gCws->fWindow );
            cursor = 0;
            break;

         case ButtonRelease :
            if (mode == 1) {
               button_press = 10+event.xbutton.button ;
               xlocp = event.xbutton.x;
               ylocp = event.xbutton.y;
            }
            break;

         case KeyPress :
            if (mode == 1) {
               button_press = event.xkey.keycode;
               xlocp = event.xbutton.x;
               ylocp = event.xbutton.y;
            }
            break;

         case KeyRelease :
            if (mode == 1) {
               button_press = -event.xkey.keycode;
               xlocp = event.xbutton.x;
               ylocp = event.xbutton.y;
            }
            break;

         default :
            break;
      }

      if (mode == 1) {
         if (button_press == 0)
            button_press = -1;
         break;
      }
   }
   x = event.xbutton.x;
   y = event.xbutton.y;

   return button_press;
}

////////////////////////////////////////////////////////////////////////////////
/// Request a string.
///
///  \param [in] x,y         : position where text is displayed
///  \param [in] text        : text displayed (input), edited text (output)
///
/// Request string:
/// text is displayed and can be edited with Emacs-like keybinding
/// return termination code (0 for ESC, 1 for RETURN)

Int_t TGX11::RequestString(int x, int y, char *text)
{
   static Cursor cursor = 0;
   static int percent = 0;  // bell volume
   Window focuswindow;
   int focusrevert;
   XEvent event;
   KeySym keysym;
   int key = -1;
   int len_text = strlen(text);
   int nt;         // defined length of text
   int pt;         // cursor position in text

   // change the cursor shape
   if (cursor == 0) {
      XKeyboardState kbstate;
      cursor = XCreateFontCursor((Display*)fDisplay, XC_question_arrow);
      XGetKeyboardControl((Display*)fDisplay, &kbstate);
      percent = kbstate.bell_percent;
   }
   if (cursor != 0)
      XDefineCursor((Display*)fDisplay, gCws->fWindow, cursor);
   for (nt = len_text; nt > 0 && text[nt-1] == ' '; nt--) { }
      pt = nt;
   XGetInputFocus((Display*)fDisplay, &focuswindow, &focusrevert);
   XSetInputFocus((Display*)fDisplay, gCws->fWindow, focusrevert, CurrentTime);
   while (key < 0) {
      char keybuf[8];
      char nbytes;
      int dx;
      int i;
      XDrawImageString((Display*)fDisplay, gCws->fWindow, *gGCtext, x, y, text, nt);
      dx = XTextWidth(gTextFont, text, nt);
      XDrawImageString((Display*)fDisplay, gCws->fWindow, *gGCtext, x + dx, y, " ", 1);
      dx = pt == 0 ? 0 : XTextWidth(gTextFont, text, pt);
      XDrawImageString((Display*)fDisplay, gCws->fWindow, *gGCinvt,
                       x + dx, y, pt < len_text ? &text[pt] : " ", 1);
      XWindowEvent((Display*)fDisplay, gCws->fWindow, gKeybdMask, &event);
      switch (event.type) {
         case ButtonPress:
         case EnterNotify:
            XSetInputFocus((Display*)fDisplay, gCws->fWindow, focusrevert, CurrentTime);
            break;
         case LeaveNotify:
            XSetInputFocus((Display*)fDisplay, focuswindow, focusrevert, CurrentTime);
            break;
         case KeyPress:
            nbytes = XLookupString(&event.xkey, keybuf, sizeof(keybuf),
                                   &keysym, 0);
            switch (keysym) {      // map cursor keys
               case XK_Left:
                  keybuf[0] = '\002';  // Control-B
                  nbytes = 1;
                  break;
               case XK_Right:
                  keybuf[0] = '\006';  // Control-F
                  nbytes = 1;
                  break;
            }
            if (nbytes == 1) {
            if (isascii(keybuf[0]) && isprint(keybuf[0])) {
               // insert character
               if (nt < len_text)
                  nt++;
               for (i = nt - 1; i > pt; i--)
                  text[i] = text[i-1];
               if (pt < len_text) {
                  text[pt] = keybuf[0];
                  pt++;
               }
            } else
               switch (keybuf[0]) {
                  // Emacs-like editing keys

                  case '\010':    // backspace
                  case '\177':    // delete
                     // delete backward
                     if (pt > 0) {
                        for (i = pt; i < nt; i++)
                           text[i-1] = text[i];
                        text[nt-1] = ' ';
                        nt--;
                        pt--;
                     }
                     break;
                  case '\001':    // ^A
                     // beginning of line
                     pt = 0;
                     break;
                  case '\002':    // ^B
                     // move backward
                     if (pt > 0)
                        pt--;
                     break;
                  case '\004':    // ^D
                     // delete forward
                     if (pt > 0) {
                        for (i = pt; i < nt; i++)
                           text[i-1] = text[i];
                        text[nt-1] = ' ';
                        pt--;
                     }
                     break;
                  case '\005':    // ^E
                     // end of line
                     pt = nt;
                     break;

                  case '\006':    // ^F
                     // move forward
                     if (pt < nt)
                        pt++;
                     break;
                  case '\013':    // ^K
                     // delete to end of line
                     for (i = pt; i < nt; i++)
                        text[i] = ' ';
                     nt = pt;
                     break;
                  case '\024':    // ^T
                     // transpose
                     if (pt > 0) {
                        char c = text[pt];
                        text[pt] = text[pt-1];
                        text[pt-1] = c;
                     }
                     break;
                  case '\012':    // newline
                  case '\015':    // return
                     key = 1;
                     break;
                  case '\033':    // escape
                     key = 0;
                     break;

                  default:
                     XBell((Display*)fDisplay, percent);
               }
            }
      }
   }
   XSetInputFocus((Display*)fDisplay, focuswindow, focusrevert, CurrentTime);

   if (cursor != 0) {
      XUndefineCursor((Display*)fDisplay, gCws->fWindow);
      cursor = 0;
   }

   return key;
}

////////////////////////////////////////////////////////////////////////////////
/// Rescale the window wid.
///
///  \param [in] wid  : Window identifier
///  \param [in] w    : Width
///  \param [in] h    : Height

void TGX11::RescaleWindow(int wid, unsigned int w, unsigned int h)
{
   int i;

   gTws = &fWindows[wid];
   if (!gTws->fOpen) return;

   // don't do anything when size did not change
   if (gTws->fWidth == w && gTws->fHeight == h) return;

   XResizeWindow((Display*)fDisplay, gTws->fWindow, w, h);

   if (gTws->fBuffer) {
      // don't free and recreate pixmap when new pixmap is smaller
      if (gTws->fWidth < w || gTws->fHeight < h) {
         XFreePixmap((Display*)fDisplay,gTws->fBuffer);
         gTws->fBuffer = XCreatePixmap((Display*)fDisplay, fRootWin, w, h, fDepth);
      }
      for (i = 0; i < kMAXGC; i++) XSetClipMask((Display*)fDisplay, gGClist[i], None);
      SetColor(gGCpxmp, 0);
      XFillRectangle( (Display*)fDisplay, gTws->fBuffer, *gGCpxmp, 0, 0, w, h);
      SetColor(gGCpxmp, 1);
      if (gTws->fDoubleBuffer) gTws->fDrawing = gTws->fBuffer;
   }
   gTws->fWidth  = w;
   gTws->fHeight = h;
}

////////////////////////////////////////////////////////////////////////////////
/// Resize a pixmap.
///
///  \param [in] wid : pixmap to be resized
///  \param [in] w,h : Width and height of the pixmap

int TGX11::ResizePixmap(int wid, unsigned int w, unsigned int h)
{
   Window root;
   unsigned int wval, hval;
   int xx, yy, i;
   unsigned int ww, hh, border, depth;
   wval = w;
   hval = h;

   gTws = &fWindows[wid];

   // don't do anything when size did not change
   //  if (gTws->fWidth == wval && gTws->fHeight == hval) return 0;

   // due to round-off errors in TPad::Resize() we might get +/- 1 pixel
   // change, in those cases don't resize pixmap
   if (gTws->fWidth  >= wval-1 && gTws->fWidth  <= wval+1 &&
       gTws->fHeight >= hval-1 && gTws->fHeight <= hval+1) return 0;

   // don't free and recreate pixmap when new pixmap is smaller
   if (gTws->fWidth < wval || gTws->fHeight < hval) {
      XFreePixmap((Display*)fDisplay, gTws->fWindow);
      gTws->fWindow = XCreatePixmap((Display*)fDisplay, fRootWin, wval, hval, fDepth);
   }
   XGetGeometry((Display*)fDisplay, gTws->fWindow, &root, &xx, &yy, &ww, &hh, &border, &depth);

   for (i = 0; i < kMAXGC; i++)
      XSetClipMask((Display*)fDisplay, gGClist[i], None);

   SetColor(gGCpxmp, 0);
   XFillRectangle((Display*)fDisplay, gTws->fWindow, *gGCpxmp, 0, 0, ww, hh);
   SetColor(gGCpxmp, 1);

   // Initialise the window structure
   gTws->fDrawing = gTws->fWindow;
   gTws->fWidth   = wval;
   gTws->fHeight  = hval;

   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Resize the current window if necessary.

void TGX11::ResizeWindow(int wid)
{
   int i;
   int xval=0, yval=0;
   Window win, root=0;
   unsigned int wval=0, hval=0, border=0, depth=0;

   gTws = &fWindows[wid];

   win = gTws->fWindow;

   XGetGeometry((Display*)fDisplay, win, &root,
                &xval, &yval, &wval, &hval, &border, &depth);
   if (wval >= 65500) wval = 1;
   if (hval >= 65500) hval = 1;

   // don't do anything when size did not change
   if (gTws->fWidth == wval && gTws->fHeight == hval) return;

   XResizeWindow((Display*)fDisplay, gTws->fWindow, wval, hval);

   if (gTws->fBuffer) {
      if (gTws->fWidth < wval || gTws->fHeight < hval) {
         XFreePixmap((Display*)fDisplay,gTws->fBuffer);
         gTws->fBuffer = XCreatePixmap((Display*)fDisplay, fRootWin, wval, hval, fDepth);
      }
      for (i = 0; i < kMAXGC; i++) XSetClipMask((Display*)fDisplay, gGClist[i], None);
      SetColor(gGCpxmp, 0);
      XFillRectangle((Display*)fDisplay, gTws->fBuffer, *gGCpxmp, 0, 0, wval, hval);
      SetColor(gGCpxmp, 1);
      if (gTws->fDoubleBuffer) gTws->fDrawing = gTws->fBuffer;
   }
   gTws->fWidth  = wval;
   gTws->fHeight = hval;
}

////////////////////////////////////////////////////////////////////////////////
/// Select window to which subsequent output is directed.

void TGX11::SelectWindow(int wid)
{
   XRectangle region;
   int i;

   if (wid < 0 || wid >= fMaxNumberOfWindows || !fWindows[wid].fOpen) return;

   gCws = &fWindows[wid];

   if (gCws->fClip && !gCws->fIsPixmap && !gCws->fDoubleBuffer) {
      region.x      = gCws->fXclip;
      region.y      = gCws->fYclip;
      region.width  = gCws->fWclip;
      region.height = gCws->fHclip;
      for (i = 0; i < kMAXGC; i++)
         XSetClipRectangles((Display*)fDisplay, gGClist[i], 0, 0, &region, 1, YXBanded);
   } else {
      for (i = 0; i < kMAXGC; i++)
         XSetClipMask((Display*)fDisplay, gGClist[i], None);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set character up vector.

void TGX11::SetCharacterUp(Float_t chupx, Float_t chupy)
{
   if (chupx == fCharacterUpX  && chupy == fCharacterUpY) return;

   if      (chupx == 0  && chupy == 0)  fTextAngle = 0;
   else if (chupx == 0  && chupy == 1)  fTextAngle = 0;
   else if (chupx == -1 && chupy == 0)  fTextAngle = 90;
   else if (chupx == 0  && chupy == -1) fTextAngle = 180;
   else if (chupx == 1  && chupy ==  0) fTextAngle = 270;
   else {
      fTextAngle = ((TMath::ACos(chupx/TMath::Sqrt(chupx*chupx +chupy*chupy))*180.)/TMath::Pi())-90;
      if (chupy < 0) fTextAngle = 180 - fTextAngle;
      if (TMath::Abs(fTextAngle) <= 0.01) fTextAngle = 0;
   }
   fCharacterUpX = chupx;
   fCharacterUpY = chupy;
}

////////////////////////////////////////////////////////////////////////////////
/// Turn off the clipping for the window wid.

void TGX11::SetClipOFF(int wid)
{
   gTws       = &fWindows[wid];
   gTws->fClip = 0;

   for (int i = 0; i < kMAXGC; i++)
      XSetClipMask( (Display*)fDisplay, gGClist[i], None );
}

////////////////////////////////////////////////////////////////////////////////
/// Set clipping region for the window wid.
///
///  \param [in] wid        : Window identifier
///  \param [in] x,y        : origin of clipping rectangle
///  \param [in] w,h        : size of clipping rectangle;

void TGX11::SetClipRegion(int wid, int x, int y, unsigned int w, unsigned int h)
{

   gTws = &fWindows[wid];
   gTws->fXclip = x;
   gTws->fYclip = y;
   gTws->fWclip = w;
   gTws->fHclip = h;
   gTws->fClip  = 1;
   if (gTws->fClip && !gTws->fIsPixmap && !gTws->fDoubleBuffer) {
      XRectangle region;
      region.x      = gTws->fXclip;
      region.y      = gTws->fYclip;
      region.width  = gTws->fWclip;
      region.height = gTws->fHclip;
      for (int i = 0; i < kMAXGC; i++)
         XSetClipRectangles((Display*)fDisplay, gGClist[i], 0, 0, &region, 1, YXBanded);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the foreground color in GC.

void  TGX11::SetColor(void *gci, int ci)
{
   GC gc = *(GC *)gci;

   TColor *color = gROOT->GetColor(ci);
   if (color)
      SetRGB(ci, color->GetRed(), color->GetGreen(), color->GetBlue());

   XColor_t &col = GetColor(ci);
   if (fColormap && !col.fDefined) {
      col = GetColor(0);
   } else if (!fColormap && (ci < 0 || ci > 1)) {
      col = GetColor(0);
   }

   if (fDrawMode == kXor) {
      XGCValues values;
      XGetGCValues((Display*)fDisplay, gc, GCBackground, &values);
      XSetForeground((Display*)fDisplay, gc, col.fPixel ^ values.background);
   } else {
      XSetForeground((Display*)fDisplay, gc, col.fPixel);

      // make sure that foreground and background are different
      XGCValues values;
      XGetGCValues((Display*)fDisplay, gc, GCForeground | GCBackground, &values);
      if (values.foreground == values.background)
         XSetBackground((Display*)fDisplay, gc, GetColor(!ci).fPixel);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the cursor.

void  TGX11::SetCursor(int wid, ECursor cursor)
{
   gTws = &fWindows[wid];
   XDefineCursor((Display*)fDisplay, gTws->fWindow, fCursors[cursor]);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the double buffer on/off on window wid.
///
///  \param [in] wid  : Window identifier.
///       - 999 means all the opened windows.
///  \param [in] mode :
///       - 1 double buffer is on
///       - 0 double buffer is off

void TGX11::SetDoubleBuffer(int wid, int mode)
{
   if (wid == 999) {
      for (int i = 0; i < fMaxNumberOfWindows; i++) {
         gTws = &fWindows[i];
         if (gTws->fOpen) {
            switch (mode) {
               case 1 :
                  SetDoubleBufferON();
                  break;
               default:
                  SetDoubleBufferOFF();
                  break;
            }
         }
      }
   } else {
      gTws = &fWindows[wid];
      if (!gTws->fOpen) return;
      switch (mode) {
         case 1 :
            SetDoubleBufferON();
            return;
         default:
            SetDoubleBufferOFF();
            return;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Turn double buffer mode off.

void TGX11::SetDoubleBufferOFF()
{
   if (!gTws->fDoubleBuffer) return;
   gTws->fDoubleBuffer = 0;
   gTws->fDrawing      = gTws->fWindow;
}

////////////////////////////////////////////////////////////////////////////////
/// Turn double buffer mode on.

void TGX11::SetDoubleBufferON()
{
   if (gTws->fDoubleBuffer || gTws->fIsPixmap) return;
   if (!gTws->fBuffer) {
      gTws->fBuffer = XCreatePixmap((Display*)fDisplay, fRootWin,
                                   gTws->fWidth, gTws->fHeight, fDepth);
      SetColor(gGCpxmp, 0);
      XFillRectangle((Display*)fDisplay, gTws->fBuffer, *gGCpxmp, 0, 0, gTws->fWidth, gTws->fHeight);
      SetColor(gGCpxmp, 1);
   }
   for (int i = 0; i < kMAXGC; i++) XSetClipMask((Display*)fDisplay, gGClist[i], None);
   gTws->fDoubleBuffer  = 1;
   gTws->fDrawing       = gTws->fBuffer;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the drawing mode.
///
///  \param [in] mode : drawing mode
///            - mode=1 copy
///            - mode=2 xor
///            - mode=3 invert
///            - mode=4 set the suitable mode for cursor echo according to
///                     the vendor

void TGX11::SetDrawMode(EDrawMode mode)
{
   int i;
   if (fDisplay) {
      switch (mode) {
         case kCopy:
            for (i = 0; i < kMAXGC; i++) XSetFunction((Display*)fDisplay, gGClist[i], GXcopy);
            break;

         case kXor:
            for (i = 0; i < kMAXGC; i++) XSetFunction((Display*)fDisplay, gGClist[i], GXxor);
            break;

         case kInvert:
            for (i = 0; i < kMAXGC; i++) XSetFunction((Display*)fDisplay, gGClist[i], GXinvert);
            break;
      }
   }
   fDrawMode = mode;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for fill areas.

void TGX11::SetFillColor(Color_t cindex)
{
   if (!gStyle->GetFillColor() && cindex > 1) cindex = 0;
   if (cindex >= 0) SetColor(gGCfill, Int_t(cindex));
   fFillColor = cindex;

   // invalidate fill pattern
   if (gFillPattern != 0) {
      XFreePixmap((Display*)fDisplay, gFillPattern);
      gFillPattern = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set fill area style.
///
///  \param [in] fstyle   : compound fill area interior style
///            - fstyle = 1000*interiorstyle + styleindex

void TGX11::SetFillStyle(Style_t fstyle)
{
   if (fFillStyle == fstyle) return;
   fFillStyle = fstyle;
   Int_t style = fstyle/1000;
   Int_t fasi  = fstyle%1000;
   SetFillStyleIndex(style,fasi);
}

////////////////////////////////////////////////////////////////////////////////
/// Set fill area style index.

void TGX11::SetFillStyleIndex(Int_t style, Int_t fasi)
{
   static int current_fasi = 0;

   fFillStyle = 1000*style + fasi;

   switch (style) {

      case 1:         // solid
         gFillHollow = 0;
         XSetFillStyle((Display*)fDisplay, *gGCfill, FillSolid);
         break;

      case 2:         // pattern
         gFillHollow = 1;
         break;

      case 3:         // hatch
         gFillHollow = 0;
         XSetFillStyle((Display*)fDisplay, *gGCfill, FillStippled);
         if (fasi != current_fasi) {
            if (gFillPattern != 0) {
               XFreePixmap((Display*)fDisplay, gFillPattern);
               gFillPattern = 0;
            }
            int stn = (fasi >= 1 && fasi <=25) ? fasi : 2;

            gFillPattern = XCreateBitmapFromData((Display*)fDisplay, fRootWin,
                                                 (const char*)gStipples[stn], 16, 16);

            XSetStipple( (Display*)fDisplay, *gGCfill, gFillPattern );
            current_fasi = fasi;
         }
         break;

      default:
         gFillHollow = 1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set input on or off.

void TGX11::SetInput(int inp)
{
   XSetWindowAttributes attributes;
   ULong_t attr_mask;

   if (inp == 1) {
      attributes.event_mask = gMouseMask | gKeybdMask;
      attr_mask = CWEventMask;
      XChangeWindowAttributes((Display*)fDisplay, gCws->fWindow, attr_mask, &attributes);
   } else {
      attributes.event_mask = NoEventMask;
      attr_mask = CWEventMask;
      XChangeWindowAttributes((Display*)fDisplay, gCws->fWindow, attr_mask, &attributes);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for lines.

void TGX11::SetLineColor(Color_t cindex)
{
   if (cindex < 0) return;

   TAttLine::SetLineColor(cindex);

   SetColor(gGCline, Int_t(cindex));
   SetColor(gGCdash, Int_t(cindex));
}

////////////////////////////////////////////////////////////////////////////////
/// Set line type.
///
///  \param [in] n         : length of dash list
///  \param [in] dash(n)   : dash segment lengths
///
///  - if n <= 0 use solid lines
///  - if n >  0 use dashed lines described by DASH(N)
///       e.g. N=4,DASH=(6,3,1,3) gives a dashed-dotted line with dash length 6
///       and a gap of 7 between dashes

void TGX11::SetLineType(int n, int *dash)
{
   if (n <= 0) {
      gLineStyle = LineSolid;
      XSetLineAttributes((Display*)fDisplay, *gGCline, gLineWidth,
                         gLineStyle, gCapStyle, gJoinStyle);
   } else {
      gDashSize = TMath::Min((int)sizeof(gDashList),n);
      gDashLength = 0;
      for (int i = 0; i < gDashSize; i++ ) {
         gDashList[i] = dash[i];
         gDashLength += gDashList[i];
      }
      gDashOffset = 0;
      gLineStyle = LineOnOffDash;
      if (gLineWidth == 0) gLineWidth =1;
      XSetLineAttributes((Display*)fDisplay, *gGCline, gLineWidth,
                         gLineStyle, gCapStyle, gJoinStyle);
      XSetLineAttributes((Display*)fDisplay, *gGCdash, gLineWidth,
                         gLineStyle, gCapStyle, gJoinStyle);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set line style.

void TGX11::SetLineStyle(Style_t lstyle)
{
   static Int_t dashed[2] = {3,3};
   static Int_t dotted[2] = {1,2};
   static Int_t dasheddotted[4] = {3,4,1,4};

   if (fLineStyle != lstyle) { //set style index only if different
      fLineStyle = lstyle;
      if (lstyle <= 1 ) {
         SetLineType(0,0);
      } else if (lstyle == 2 ) {
         SetLineType(2,dashed);
      } else if (lstyle == 3 ) {
         SetLineType(2,dotted);
      } else if (lstyle == 4 ) {
         SetLineType(4,dasheddotted);
      } else {
         TString st = (TString)gStyle->GetLineStyleString(lstyle);
         TObjArray *tokens = st.Tokenize(" ");
         Int_t nt;
         nt = tokens->GetEntries();
         Int_t *linestyle = new Int_t[nt];
         for (Int_t j = 0; j<nt; j++) {
            Int_t it;
            sscanf(((TObjString*)tokens->At(j))->GetName(), "%d", &it);
            linestyle[j] = (Int_t)(it/4);
         }
         SetLineType(nt,linestyle);
         delete [] linestyle;
         delete tokens;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set line width.
///
///  \param [in] width   : line width in pixels

void TGX11::SetLineWidth(Width_t width )
{
   if (fLineWidth == width) return;
   fLineWidth = width;

   if (width == 1 && gLineStyle == LineSolid) gLineWidth = 0;
   else                                       gLineWidth = width;

   if (gLineWidth < 0) return;

   XSetLineAttributes((Display*)fDisplay, *gGCline, gLineWidth,
                      gLineStyle, gCapStyle, gJoinStyle);
   XSetLineAttributes((Display*)fDisplay, *gGCdash, gLineWidth,
              gLineStyle, gCapStyle, gJoinStyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for markers.

void TGX11::SetMarkerColor(Color_t cindex)
{
   if (cindex < 0) return;

   TAttMarker::SetMarkerColor(cindex);

   SetColor(gGCmark, Int_t(cindex));
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker size index.
///
///  \param [in] msize  : marker scale factor

void TGX11::SetMarkerSize(Float_t msize)
{
   if (msize == fMarkerSize) return;

   fMarkerSize = msize;
   if (msize < 0) return;

   SetMarkerStyle(-fMarkerStyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker type.
///
///  \param [in] type      : marker type
///  \param [in] n         : length of marker description
///  \param [in] xy        : list of points describing marker shape
///
///  - if n == 0 marker is a single point
///  - if TYPE == 0 marker is hollow circle of diameter N
///  - if TYPE == 1 marker is filled circle of diameter N
///  - if TYPE == 2 marker is a hollow polygon describe by line XY
///  - if TYPE == 3 marker is a filled polygon describe by line XY
///  - if TYPE == 4 marker is described by segmented line XY
///     e.g. TYPE=4,N=4,XY=(-3,0,3,0,0,-3,0,3) sets a plus shape of 7x7 pixels

void TGX11::SetMarkerType(int type, int n, RXPoint *xy)
{
   gMarker.type = type;
   gMarker.n = n < kMAXMK ? n : kMAXMK;
   if (gMarker.type >= 2) {
      for (int i = 0; i < gMarker.n; i++) {
         gMarker.xy[i].x = xy[i].x;
         gMarker.xy[i].y = xy[i].y;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker style.

void TGX11::SetMarkerStyle(Style_t markerstyle)
{
   if (fMarkerStyle == markerstyle) return;
   static RXPoint shape[30];
   fMarkerStyle = TMath::Abs(markerstyle);
   markerstyle = TAttMarker::GetMarkerStyleBase(fMarkerStyle);
   gMarkerLineWidth = TAttMarker::GetMarkerLineWidth(fMarkerStyle);

   // The fast pixel markers need to be treated separately
   if (markerstyle == 1 || markerstyle == 6 || markerstyle == 7) {
       XSetLineAttributes((Display*)fDisplay, *gGCmark, 0, LineSolid, CapButt, JoinMiter);
   } else {
       XSetLineAttributes((Display*)fDisplay, *gGCmark, gMarkerLineWidth,
                          gMarkerLineStyle, gMarkerCapStyle, gMarkerJoinStyle);
   }

   Float_t MarkerSizeReduced = fMarkerSize - TMath::Floor(gMarkerLineWidth/2.)/4.;
   Int_t im = Int_t(4*MarkerSizeReduced + 0.5);
   if (markerstyle == 2) {
      // + shaped marker
      shape[0].x = -im;  shape[0].y = 0;
      shape[1].x =  im;  shape[1].y = 0;
      shape[2].x = 0  ;  shape[2].y = -im;
      shape[3].x = 0  ;  shape[3].y = im;
      SetMarkerType(4,4,shape);
   } else if (markerstyle == 3 || markerstyle == 31) {
      // * shaped marker
      shape[0].x = -im;  shape[0].y = 0;
      shape[1].x =  im;  shape[1].y = 0;
      shape[2].x = 0  ;  shape[2].y = -im;
      shape[3].x = 0  ;  shape[3].y = im;
      im = Int_t(0.707*Float_t(im) + 0.5);
      shape[4].x = -im;  shape[4].y = -im;
      shape[5].x =  im;  shape[5].y = im;
      shape[6].x = -im;  shape[6].y = im;
      shape[7].x =  im;  shape[7].y = -im;
      SetMarkerType(4,8,shape);
   } else if (markerstyle == 4 || markerstyle == 24) {
      // O shaped marker
      SetMarkerType(0,im*2,shape);
   } else if (markerstyle == 5) {
      // X shaped marker
      im = Int_t(0.707*Float_t(im) + 0.5);
      shape[0].x = -im;  shape[0].y = -im;
      shape[1].x =  im;  shape[1].y = im;
      shape[2].x = -im;  shape[2].y = im;
      shape[3].x =  im;  shape[3].y = -im;
      SetMarkerType(4,4,shape);
   } else if (markerstyle == 6) {
      // + shaped marker (with 1 pixel)
      shape[0].x = -1 ;  shape[0].y = 0;
      shape[1].x =  1 ;  shape[1].y = 0;
      shape[2].x =  0 ;  shape[2].y = -1;
      shape[3].x =  0 ;  shape[3].y = 1;
      SetMarkerType(4,4,shape);
   } else if (markerstyle == 7) {
      // . shaped marker (with 9 pixel)
      shape[0].x = -1 ;  shape[0].y = 1;
      shape[1].x =  1 ;  shape[1].y = 1;
      shape[2].x = -1 ;  shape[2].y = 0;
      shape[3].x =  1 ;  shape[3].y = 0;
      shape[4].x = -1 ;  shape[4].y = -1;
      shape[5].x =  1 ;  shape[5].y = -1;
      SetMarkerType(4,6,shape);
   } else if (markerstyle == 8 || markerstyle == 20) {
      // O shaped marker (filled)
      SetMarkerType(1,im*2,shape);
   } else if (markerstyle == 21) {
      // full square
      shape[0].x = -im;  shape[0].y = -im;
      shape[1].x =  im;  shape[1].y = -im;
      shape[2].x =  im;  shape[2].y = im;
      shape[3].x = -im;  shape[3].y = im;
      shape[4].x = -im;  shape[4].y = -im;
      SetMarkerType(3,5,shape);
   } else if (markerstyle == 22) {
      // full triangle up
      shape[0].x = -im;  shape[0].y = im;
      shape[1].x =  im;  shape[1].y = im;
      shape[2].x =   0;  shape[2].y = -im;
      shape[3].x = -im;  shape[3].y = im;
      SetMarkerType(3,4,shape);
   } else if (markerstyle == 23) {
      // full triangle down
      shape[0].x =   0;  shape[0].y = im;
      shape[1].x =  im;  shape[1].y = -im;
      shape[2].x = -im;  shape[2].y = -im;
      shape[3].x =   0;  shape[3].y = im;
      SetMarkerType(3,4,shape);
   } else if (markerstyle == 25) {
      // open square
      shape[0].x = -im;  shape[0].y = -im;
      shape[1].x =  im;  shape[1].y = -im;
      shape[2].x =  im;  shape[2].y = im;
      shape[3].x = -im;  shape[3].y = im;
      shape[4].x = -im;  shape[4].y = -im;
      SetMarkerType(2,5,shape);
   } else if (markerstyle == 26) {
      // open triangle up
      shape[0].x = -im;  shape[0].y = im;
      shape[1].x =  im;  shape[1].y = im;
      shape[2].x =   0;  shape[2].y = -im;
      shape[3].x = -im;  shape[3].y = im;
      SetMarkerType(2,4,shape);
   } else if (markerstyle == 27) {
      // open losange
      Int_t imx = Int_t(2.66*MarkerSizeReduced + 0.5);
      shape[0].x =-imx;  shape[0].y = 0;
      shape[1].x =   0;  shape[1].y = -im;
      shape[2].x = imx;  shape[2].y = 0;
      shape[3].x =   0;  shape[3].y = im;
      shape[4].x =-imx;  shape[4].y = 0;
      SetMarkerType(2,5,shape);
   } else if (markerstyle == 28) {
      // open cross
      Int_t imx = Int_t(1.33*MarkerSizeReduced + 0.5);
      shape[0].x = -im;  shape[0].y =-imx;
      shape[1].x =-imx;  shape[1].y =-imx;
      shape[2].x =-imx;  shape[2].y = -im;
      shape[3].x = imx;  shape[3].y = -im;
      shape[4].x = imx;  shape[4].y =-imx;
      shape[5].x =  im;  shape[5].y =-imx;
      shape[6].x =  im;  shape[6].y = imx;
      shape[7].x = imx;  shape[7].y = imx;
      shape[8].x = imx;  shape[8].y = im;
      shape[9].x =-imx;  shape[9].y = im;
      shape[10].x=-imx;  shape[10].y= imx;
      shape[11].x= -im;  shape[11].y= imx;
      shape[12].x= -im;  shape[12].y=-imx;
      SetMarkerType(2,13,shape);
   } else if (markerstyle == 29) {
      // full star pentagone
      Int_t im1 = Int_t(0.66*MarkerSizeReduced + 0.5);
      Int_t im2 = Int_t(2.00*MarkerSizeReduced + 0.5);
      Int_t im3 = Int_t(2.66*MarkerSizeReduced + 0.5);
      Int_t im4 = Int_t(1.33*MarkerSizeReduced + 0.5);
      shape[0].x = -im;  shape[0].y = im4;
      shape[1].x =-im2;  shape[1].y =-im1;
      shape[2].x =-im3;  shape[2].y = -im;
      shape[3].x =   0;  shape[3].y =-im2;
      shape[4].x = im3;  shape[4].y = -im;
      shape[5].x = im2;  shape[5].y =-im1;
      shape[6].x =  im;  shape[6].y = im4;
      shape[7].x = im4;  shape[7].y = im4;
      shape[8].x =   0;  shape[8].y = im;
      shape[9].x =-im4;  shape[9].y = im4;
      shape[10].x= -im;  shape[10].y= im4;
      SetMarkerType(3,11,shape);
   } else if (markerstyle == 30) {
      // open star pentagone
      Int_t im1 = Int_t(0.66*MarkerSizeReduced + 0.5);
      Int_t im2 = Int_t(2.00*MarkerSizeReduced + 0.5);
      Int_t im3 = Int_t(2.66*MarkerSizeReduced + 0.5);
      Int_t im4 = Int_t(1.33*MarkerSizeReduced + 0.5);
      shape[0].x = -im;  shape[0].y = im4;
      shape[1].x =-im2;  shape[1].y =-im1;
      shape[2].x =-im3;  shape[2].y = -im;
      shape[3].x =   0;  shape[3].y =-im2;
      shape[4].x = im3;  shape[4].y = -im;
      shape[5].x = im2;  shape[5].y =-im1;
      shape[6].x =  im;  shape[6].y = im4;
      shape[7].x = im4;  shape[7].y = im4;
      shape[8].x =   0;  shape[8].y = im;
      shape[9].x =-im4;  shape[9].y = im4;
      shape[10].x= -im;  shape[10].y= im4;
      SetMarkerType(2,11,shape);
   } else if (markerstyle == 32) {
      // open triangle down
      shape[0].x =   0;  shape[0].y = im;
      shape[1].x =  im;  shape[1].y = -im;
      shape[2].x = -im;  shape[2].y = -im;
      shape[3].x =   0;  shape[3].y = im;
      SetMarkerType(2,4,shape);
   } else if (markerstyle == 33) {
      // full losange
      Int_t imx = Int_t(2.66*MarkerSizeReduced + 0.5);
      shape[0].x =-imx;  shape[0].y = 0;
      shape[1].x =   0;  shape[1].y = -im;
      shape[2].x = imx;  shape[2].y = 0;
      shape[3].x =   0;  shape[3].y = im;
      shape[4].x =-imx;  shape[4].y = 0;
      SetMarkerType(3,5,shape);
   } else if (markerstyle == 34) {
      // full cross
      Int_t imx = Int_t(1.33*MarkerSizeReduced + 0.5);
      shape[0].x = -im;  shape[0].y =-imx;
      shape[1].x =-imx;  shape[1].y =-imx;
      shape[2].x =-imx;  shape[2].y = -im;
      shape[3].x = imx;  shape[3].y = -im;
      shape[4].x = imx;  shape[4].y =-imx;
      shape[5].x =  im;  shape[5].y =-imx;
      shape[6].x =  im;  shape[6].y = imx;
      shape[7].x = imx;  shape[7].y = imx;
      shape[8].x = imx;  shape[8].y = im;
      shape[9].x =-imx;  shape[9].y = im;
      shape[10].x=-imx;  shape[10].y= imx;
      shape[11].x= -im;  shape[11].y= imx;
      shape[12].x= -im;  shape[12].y=-imx;
      SetMarkerType(3,13,shape);
   } else if (markerstyle == 35) {
      // diamond with cross
      shape[0].x =-im;  shape[0].y = 0;
      shape[1].x =  0;  shape[1].y = -im;
      shape[2].x = im;  shape[2].y = 0;
      shape[3].x =  0;  shape[3].y = im;
      shape[4].x =-im;  shape[4].y = 0;
      shape[5].x = im;  shape[5].y = 0;
      shape[6].x =  0;  shape[6].y = im;
      shape[7].x =  0;  shape[7].y =-im;
      SetMarkerType(2,8,shape);
   } else if (markerstyle == 36) {
      // square with diagonal cross
      shape[0].x = -im;  shape[0].y = -im;
      shape[1].x =  im;  shape[1].y = -im;
      shape[2].x =  im;  shape[2].y = im;
      shape[3].x = -im;  shape[3].y = im;
      shape[4].x = -im;  shape[4].y = -im;
      shape[5].x =  im;  shape[5].y = im;
      shape[6].x = -im;  shape[6].y = im;
      shape[7].x =  im;  shape[7].y = -im;
      SetMarkerType(2,8,shape);
   } else if (markerstyle == 37) {
      // open three triangles
      Int_t im2 = Int_t(2.0*MarkerSizeReduced + 0.5);
      shape[0].x =   0;  shape[0].y =   0;
      shape[1].x =-im2;  shape[1].y =  im;
      shape[2].x = im2;  shape[2].y =  im;
      shape[3].x =   0;  shape[3].y =   0;
      shape[4].x =-im2;  shape[4].y = -im;
      shape[5].x = -im;  shape[5].y =   0;
      shape[6].x =   0;  shape[6].y =   0;
      shape[7].x =  im;  shape[7].y =   0;
      shape[8].x = im2;  shape[8].y =  -im;
      shape[9].x =   0;  shape[9].y =   0;
      SetMarkerType(2, 10,shape);
   } else if (markerstyle == 38) {
      // + shaped marker with octagon
      Int_t im2 = Int_t(2.0*MarkerSizeReduced + 0.5);
      shape[0].x = -im;  shape[0].y = 0;
      shape[1].x = -im;  shape[1].y =-im2;
      shape[2].x =-im2;  shape[2].y = -im;
      shape[3].x = im2;  shape[3].y = -im;
      shape[4].x =  im;  shape[4].y =-im2;
      shape[5].x =  im;  shape[5].y = im2;
      shape[6].x = im2;  shape[6].y = im;
      shape[7].x =-im2;  shape[7].y = im;
      shape[8].x = -im;  shape[8].y = im2;
      shape[9].x = -im;  shape[9].y = 0;
      shape[10].x = im;  shape[10].y = 0;
      shape[11].x =  0;  shape[11].y = 0;
      shape[12].x =  0;  shape[12].y = -im;
      shape[13].x =  0;  shape[13].y = im;
      shape[14].x =  0;  shape[14].y = 0;
      SetMarkerType(2,15,shape);
   } else if (markerstyle == 39) {
      // filled three triangles
      Int_t im2 = Int_t(2.0*MarkerSizeReduced + 0.5);
      shape[0].x =   0;  shape[0].y =   0;
      shape[1].x =-im2;  shape[1].y =  im;
      shape[2].x = im2;  shape[2].y =  im;
      shape[3].x =   0;  shape[3].y =   0;
      shape[4].x =-im2;  shape[4].y = -im;
      shape[5].x = -im;  shape[5].y =   0;
      shape[6].x =   0;  shape[6].y =   0;
      shape[7].x =  im;  shape[7].y =   0;
      shape[8].x = im2;  shape[8].y =  -im;
      SetMarkerType(3,9,shape);
   } else if (markerstyle == 40) {
      // four open triangles X
      Int_t im2 = Int_t(2.0*MarkerSizeReduced + 0.5);
      shape[0].x =     0;  shape[0].y =    0;
      shape[1].x =   im2;  shape[1].y =   im;
      shape[2].x =    im;  shape[2].y =  im2;
      shape[3].x =     0;  shape[3].y =    0;
      shape[4].x =    im;  shape[4].y = -im2;
      shape[5].x =   im2;  shape[5].y =  -im;
      shape[6].x =     0;  shape[6].y =    0;
      shape[7].x =  -im2;  shape[7].y =  -im;
      shape[8].x =   -im;  shape[8].y = -im2;
      shape[9].x =     0;  shape[9].y =    0;
      shape[10].x =   -im;  shape[10].y =  im2;
      shape[11].x =  -im2;  shape[11].y =   im;
      shape[12].x =     0;  shape[12].y =  0;
      SetMarkerType(2,13,shape);
   } else if (markerstyle == 41) {
      // four filled triangles X
      Int_t im2 = Int_t(2.0*MarkerSizeReduced + 0.5);
      shape[0].x =     0;  shape[0].y =    0;
      shape[1].x =   im2;  shape[1].y =   im;
      shape[2].x =    im;  shape[2].y =  im2;
      shape[3].x =     0;  shape[3].y =    0;
      shape[4].x =    im;  shape[4].y = -im2;
      shape[5].x =   im2;  shape[5].y =  -im;
      shape[6].x =     0;  shape[6].y =    0;
      shape[7].x =  -im2;  shape[7].y =  -im;
      shape[8].x =   -im;  shape[8].y = -im2;
      shape[9].x =     0;  shape[9].y =    0;
      shape[10].x =   -im;  shape[10].y =  im2;
      shape[11].x =  -im2;  shape[11].y =   im;
      shape[12].x =     0;  shape[12].y =  0;
      SetMarkerType(3,13,shape);
   } else if (markerstyle == 42) {
      // open double diamonds
      Int_t imx = Int_t(MarkerSizeReduced + 0.5);
      shape[0].x=     0;   shape[0].y= im;
      shape[1].x=  -imx;   shape[1].y= imx;
      shape[2].x  = -im;   shape[2].y = 0;
      shape[3].x = -imx;   shape[3].y = -imx;
      shape[4].x =    0;   shape[4].y = -im;
      shape[5].x =  imx;   shape[5].y = -imx;
      shape[6].x =   im;   shape[6].y = 0;
      shape[7].x=   imx;   shape[7].y= imx;
      shape[8].x=     0;   shape[8].y= im;
      SetMarkerType(2,9,shape);
   } else if (markerstyle == 43) {
      // filled double diamonds
      Int_t imx = Int_t(MarkerSizeReduced + 0.5);
      shape[0].x =    0;   shape[0].y =   im;
      shape[1].x = -imx;   shape[1].y =  imx;
      shape[2].x =  -im;   shape[2].y =    0;
      shape[3].x = -imx;   shape[3].y = -imx;
      shape[4].x =    0;   shape[4].y =  -im;
      shape[5].x =  imx;   shape[5].y = -imx;
      shape[6].x =   im;   shape[6].y =    0;
      shape[7].x =  imx;   shape[7].y =  imx;
      shape[8].x =    0;   shape[8].y =   im;
      SetMarkerType(3,9,shape);
   } else if (markerstyle == 44) {
      // open four triangles plus
      Int_t im2 = Int_t(2.0*MarkerSizeReduced + 0.5);
      shape[0].x =    0;  shape[0].y =    0;
      shape[1].x =  im2;  shape[1].y =   im;
      shape[2].x = -im2;  shape[2].y =   im;
      shape[3].x =  im2;  shape[3].y =  -im;
      shape[4].x = -im2;  shape[4].y =  -im;
      shape[5].x =    0;  shape[5].y =    0;
      shape[6].x =   im;  shape[6].y =  im2;
      shape[7].x =   im;  shape[7].y = -im2;
      shape[8].x =  -im;  shape[8].y =  im2;
      shape[9].x =  -im;  shape[9].y = -im2;
      shape[10].x =    0;  shape[10].y =    0;
      SetMarkerType(2,11,shape);
   } else if (markerstyle == 45) {
      // filled four triangles plus
      Int_t im0 = Int_t(0.4*MarkerSizeReduced + 0.5);
      Int_t im2 = Int_t(2.0*MarkerSizeReduced + 0.5);
      shape[0].x =  im0;  shape[0].y =  im0;
      shape[1].x =  im2;  shape[1].y =   im;
      shape[2].x = -im2;  shape[2].y =   im;
      shape[3].x = -im0;  shape[3].y =  im0;
      shape[4].x =  -im;  shape[4].y =  im2;
      shape[5].x =  -im;  shape[5].y = -im2;
      shape[6].x = -im0;  shape[6].y = -im0;
      shape[7].x = -im2;  shape[7].y =  -im;
      shape[8].x =  im2;  shape[8].y =  -im;
      shape[9].x =  im0;  shape[9].y = -im0;
      shape[10].x =   im;  shape[10].y = -im2;
      shape[11].x =   im;  shape[11].y =  im2;
      shape[12].x =  im0;  shape[12].y =  im0;
      SetMarkerType(3,13,shape);
   } else if (markerstyle == 46) {
      // open four triangles X
      Int_t im2 = Int_t(2.0*MarkerSizeReduced + 0.5);
      shape[0].x =    0;  shape[0].y =  im2;
      shape[1].x = -im2;  shape[1].y =   im;
      shape[2].x =  -im;  shape[2].y =  im2;
      shape[3].x = -im2;  shape[3].y =    0;
      shape[4].x =  -im;  shape[4].y = -im2;
      shape[5].x = -im2;  shape[5].y =  -im;
      shape[6].x =    0;  shape[6].y = -im2;
      shape[7].x =  im2;  shape[7].y =  -im;
      shape[8].x =   im;  shape[8].y = -im2;
      shape[9].x =  im2;  shape[9].y =    0;
      shape[10].x =  im;  shape[10].y = im2;
      shape[11].x = im2;  shape[11].y =  im;
      shape[12].x =   0;  shape[12].y = im2;
      SetMarkerType(2,13,shape);
   } else if (markerstyle == 47) {
      // filled four triangles X
      Int_t im2 = Int_t(2.0*MarkerSizeReduced + 0.5);
      shape[0].x =    0;  shape[0].y =  im2;
      shape[1].x = -im2;  shape[1].y =   im;
      shape[2].x =  -im;  shape[2].y =  im2;
      shape[3].x = -im2;  shape[3].y =    0;
      shape[4].x =  -im;  shape[4].y = -im2;
      shape[5].x = -im2;  shape[5].y =  -im;
      shape[6].x =    0;  shape[6].y = -im2;
      shape[7].x =  im2;  shape[7].y =  -im;
      shape[8].x =   im;  shape[8].y = -im2;
      shape[9].x =  im2;  shape[9].y =    0;
      shape[10].x =  im;  shape[10].y = im2;
      shape[11].x = im2;  shape[11].y =  im;
      shape[12].x =   0;  shape[12].y = im2;
      SetMarkerType(3,13,shape);
   } else if (markerstyle == 48) {
      // four filled squares X
      Int_t im2 = Int_t(2.0*MarkerSizeReduced + 0.5);
      shape[0].x =    0;  shape[0].y =  im2*1.005;
      shape[1].x = -im2;  shape[1].y =   im;
      shape[2].x =  -im;  shape[2].y =  im2;
      shape[3].x = -im2;  shape[3].y =    0;
      shape[4].x =  -im;  shape[4].y = -im2;
      shape[5].x = -im2;  shape[5].y =  -im;
      shape[6].x =    0;  shape[6].y = -im2;
      shape[7].x =  im2;  shape[7].y =  -im;
      shape[8].x =   im;  shape[8].y = -im2;
      shape[9].x =  im2;  shape[9].y =    0;
      shape[10].x =  im;  shape[10].y = im2;
      shape[11].x = im2;  shape[11].y =  im;
      shape[12].x =   0;  shape[12].y = im2*0.995;
      shape[13].x =  im2*0.995;  shape[13].y =    0;
      shape[14].x =    0;  shape[14].y = -im2*0.995;
      shape[15].x = -im2*0.995;  shape[15].y =    0;
      shape[16].x =    0;  shape[16].y =  im2*0.995;
      SetMarkerType(3,16,shape);
   } else if (markerstyle == 49) {
      // four filled squares plus
      Int_t imx = Int_t(1.33*MarkerSizeReduced + 0.5);
      shape[0].x =-imx;  shape[0].y =-imx*1.005;
      shape[1].x =-imx;  shape[1].y = -im;
      shape[2].x = imx;  shape[2].y = -im;
      shape[3].x = imx;  shape[3].y =-imx;
      shape[4].x =  im;  shape[4].y =-imx;
      shape[5].x =  im;  shape[5].y = imx;
      shape[6].x = imx;  shape[6].y = imx;
      shape[7].x = imx;  shape[7].y = im;
      shape[8].x =-imx;  shape[8].y = im;
      shape[9].x =-imx;  shape[9].y = imx;
      shape[10].x = -im;  shape[10].y = imx;
      shape[11].x = -im;  shape[11].y =-imx;
      shape[12].x =-imx;  shape[12].y =-imx*0.995;
      shape[13].x =-imx;  shape[13].y = imx;
      shape[14].x = imx;  shape[14].y = imx;
      shape[15].x = imx;  shape[15].y =-imx;
      shape[16].x =-imx;  shape[16].y =-imx*1.005;
      SetMarkerType(3,17,shape);
   } else {
      // single dot
      SetMarkerType(0,0,shape);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set opacity of a window. This image manipulation routine works
/// by adding to a percent amount of neutral to each pixels RGB.
/// Since it requires quite some additional color map entries is it
/// only supported on displays with more than > 8 color planes (> 256
/// colors).

void TGX11::SetOpacity(Int_t percent)
{
   if (fDepth <= 8) return;
   if (percent == 0) return;
   // if 100 percent then just make white

   ULong_t *orgcolors = 0, *tmpc = 0;
   Int_t    maxcolors = 0, ncolors = 0, ntmpc = 0;

   // save previous allocated colors, delete at end when not used anymore
   if (gCws->fNewColors) {
      tmpc = gCws->fNewColors;
      ntmpc = gCws->fNcolors;
   }

   // get pixmap from server as image
   XImage *image = XGetImage((Display*)fDisplay, gCws->fDrawing, 0, 0, gCws->fWidth,
                             gCws->fHeight, AllPlanes, ZPixmap);
   if (!image) return;
   // collect different image colors
   int x, y;
   for (y = 0; y < (int) gCws->fHeight; y++) {
      for (x = 0; x < (int) gCws->fWidth; x++) {
         ULong_t pixel = XGetPixel(image, x, y);
         CollectImageColors(pixel, orgcolors, ncolors, maxcolors);
      }
   }
   if (ncolors == 0) {
      XDestroyImage(image);
      ::operator delete(orgcolors);
      return;
   }

   // create opaque counter parts
   MakeOpaqueColors(percent, orgcolors, ncolors);

   if (gCws->fNewColors) {
      // put opaque colors in image
      for (y = 0; y < (int) gCws->fHeight; y++) {
         for (x = 0; x < (int) gCws->fWidth; x++) {
            ULong_t pixel = XGetPixel(image, x, y);
            Int_t idx = FindColor(pixel, orgcolors, ncolors);
            XPutPixel(image, x, y, gCws->fNewColors[idx]);
         }
      }
   }

   // put image back in pixmap on server
   XPutImage((Display*)fDisplay, gCws->fDrawing, *gGCpxmp, image, 0, 0, 0, 0,
             gCws->fWidth, gCws->fHeight);
   XFlush((Display*)fDisplay);

   // clean up
   if (tmpc) {
      if (fRedDiv == -1)
         XFreeColors((Display*)fDisplay, fColormap, tmpc, ntmpc, 0);
      delete [] tmpc;
   }
   XDestroyImage(image);
   ::operator delete(orgcolors);
}

////////////////////////////////////////////////////////////////////////////////
/// Collect in orgcolors all different original image colors.

void TGX11::CollectImageColors(ULong_t pixel, ULong_t *&orgcolors, Int_t &ncolors,
                               Int_t &maxcolors)
{
   if (maxcolors == 0) {
      ncolors   = 0;
      maxcolors = 100;
      orgcolors = (ULong_t*) ::operator new(maxcolors*sizeof(ULong_t));
   }

   for (int i = 0; i < ncolors; i++)
      if (pixel == orgcolors[i]) return;

   if (ncolors >= maxcolors) {
      orgcolors = (ULong_t*) TStorage::ReAlloc(orgcolors,
          maxcolors*2*sizeof(ULong_t), maxcolors*sizeof(ULong_t));
      maxcolors *= 2;
   }

   orgcolors[ncolors++] = pixel;
}

////////////////////////////////////////////////////////////////////////////////
/// Get RGB values for orgcolors, add percent neutral to the RGB and
/// allocate fNewColors.

void TGX11::MakeOpaqueColors(Int_t percent, ULong_t *orgcolors, Int_t ncolors)
{
   if (ncolors == 0) return;

   RXColor *xcol = new RXColor[ncolors];

   int i;
   for (i = 0; i < ncolors; i++) {
      xcol[i].pixel = orgcolors[i];
      xcol[i].red   = xcol[i].green = xcol[i].blue = 0;
      xcol[i].flags = DoRed | DoGreen | DoBlue;
   }
   QueryColors(fColormap, xcol, ncolors);

   UShort_t add = percent * kBIGGEST_RGB_VALUE / 100;

   Int_t val;
   for (i = 0; i < ncolors; i++) {
      val = xcol[i].red + add;
      if (val > kBIGGEST_RGB_VALUE) val = kBIGGEST_RGB_VALUE;
      xcol[i].red = (UShort_t) val;
      val = xcol[i].green + add;
      if (val > kBIGGEST_RGB_VALUE) val = kBIGGEST_RGB_VALUE;
      xcol[i].green = (UShort_t) val;
      val = xcol[i].blue + add;
      if (val > kBIGGEST_RGB_VALUE) val = kBIGGEST_RGB_VALUE;
      xcol[i].blue = (UShort_t) val;
      if (!AllocColor(fColormap, &xcol[i]))
         Warning("MakeOpaqueColors", "failed to allocate color %hd, %hd, %hd",
                 xcol[i].red, xcol[i].green, xcol[i].blue);
      // assumes that in case of failure xcol[i].pixel is not changed
   }

   gCws->fNewColors = new ULong_t[ncolors];
   gCws->fNcolors   = ncolors;

   for (i = 0; i < ncolors; i++)
      gCws->fNewColors[i] = xcol[i].pixel;

   delete [] xcol;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns index in orgcolors (and fNewColors) for pixel.

Int_t TGX11::FindColor(ULong_t pixel, ULong_t *orgcolors, Int_t ncolors)
{
   for (int i = 0; i < ncolors; i++)
      if (pixel == orgcolors[i]) return i;

   Error("FindColor", "did not find color, should never happen!");

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color intensities for given color index.
///
///  \param [in] cindex     : color index
///  \param [in] r,g,b      : red, green, blue intensities between 0.0 and 1.0

void TGX11::SetRGB(int cindex, float r, float g, float b)
{
   if (fColormap) {
      RXColor xcol;
      xcol.red   = (UShort_t)(r * kBIGGEST_RGB_VALUE);
      xcol.green = (UShort_t)(g * kBIGGEST_RGB_VALUE);
      xcol.blue  = (UShort_t)(b * kBIGGEST_RGB_VALUE);
      xcol.flags = DoRed | DoGreen | DoBlue;
      XColor_t &col = GetColor(cindex);
      if (col.fDefined) {
         // if color is already defined with same rgb just return
         if (col.fRed  == xcol.red && col.fGreen == xcol.green &&
             col.fBlue == xcol.blue)
            return;
         col.fDefined = kFALSE;
         if (fRedDiv == -1)
            XFreeColors((Display*)fDisplay, fColormap, &col.fPixel, 1, 0);
      }
      if (AllocColor(fColormap, &xcol)) {
         col.fDefined = kTRUE;
         col.fPixel   = xcol.pixel;
         col.fRed     = xcol.red;
         col.fGreen   = xcol.green;
         col.fBlue    = xcol.blue;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set text alignment.
///
///  \param [in] talign   text alignment

void TGX11::SetTextAlign(Short_t talign)
{
   Int_t txalh = talign/10;
   Int_t txalv = talign%10;
   fTextAlignH = txalh;
   fTextAlignV = txalv;

   switch (txalh) {

      case 0 :
      case 1 :
         switch (txalv) {  //left
            case 1 :
               fTextAlign = 7;   //bottom
               break;
            case 2 :
               fTextAlign = 4;   //center
               break;
            case 3 :
               fTextAlign = 1;   //top
               break;
         }
         break;
      case 2 :
         switch (txalv) { //center
            case 1 :
               fTextAlign = 8;   //bottom
               break;
            case 2 :
               fTextAlign = 5;   //center
               break;
            case 3 :
               fTextAlign = 2;   //top
               break;
         }
         break;
      case 3 :
         switch (txalv) {  //right
            case 1 :
               fTextAlign = 9;   //bottom
               break;
            case 2 :
               fTextAlign = 6;   //center
               break;
            case 3 :
               fTextAlign = 3;   //top
               break;
         }
         break;
   }

   TAttText::SetTextAlign(fTextAlign);
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for text.

void TGX11::SetTextColor(Color_t cindex)
{
   if (cindex < 0) return;

   TAttText::SetTextColor(cindex);

   SetColor(gGCtext, Int_t(cindex));

   XGCValues values;
   if (XGetGCValues((Display*)fDisplay, *gGCtext, GCForeground | GCBackground, &values)) {
      XSetForeground( (Display*)fDisplay, *gGCinvt, values.background );
      XSetBackground( (Display*)fDisplay, *gGCinvt, values.foreground );
   } else {
      Error("SetTextColor", "cannot get GC values");
   }
   XSetBackground((Display*)fDisplay, *gGCtext, GetColor(0).fPixel);
}

////////////////////////////////////////////////////////////////////////////////
/// Set text font to specified name.
///
///  \param [in] fontname   font name
///  \param [in] mode       loading flag
///            - mode=0     search if the font exist (kCheck)
///            - mode=1     search the font and load it if it exists (kLoad)
///
/// Set text font to specified name. This function returns 0 if
/// the specified font is found, 1 if not.

Int_t TGX11::SetTextFont(char *fontname, ETextSetMode mode)
{
   char **fontlist;
   int fontcount;
   int i;

   if (mode == kLoad) {
      for (i = 0; i < kMAXFONT; i++) {
         if (strcmp(fontname, gFont[i].name) == 0) {
            gTextFont = gFont[i].id;
            XSetFont((Display*)fDisplay, *gGCtext, gTextFont->fid);
            XSetFont((Display*)fDisplay, *gGCinvt, gTextFont->fid);
            return 0;
         }
      }
   }

   fontlist = XListFonts((Display*)fDisplay, fontname, 1, &fontcount);

   if (fontlist && fontcount != 0) {
      if (mode == kLoad) {
         if (gFont[gCurrentFontNumber].id)
            XFreeFont((Display*)fDisplay, gFont[gCurrentFontNumber].id);
         gTextFont = XLoadQueryFont((Display*)fDisplay, fontlist[0]);
         XSetFont((Display*)fDisplay, *gGCtext, gTextFont->fid);
         XSetFont((Display*)fDisplay, *gGCinvt, gTextFont->fid);
         gFont[gCurrentFontNumber].id = gTextFont;
         strlcpy(gFont[gCurrentFontNumber].name,fontname,80);
         gCurrentFontNumber++;
         if (gCurrentFontNumber == kMAXFONT) gCurrentFontNumber = 0;
      }
      XFreeFontNames(fontlist);
      return 0;
   } else {
      return 1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set current text font number.

void TGX11::SetTextFont(Font_t fontnumber)
{
   fTextFont = fontnumber;
}

////////////////////////////////////////////////////////////////////////////////
/// Set current text size.

void TGX11::SetTextSize(Float_t textsize)
{
   fTextSize = textsize;
}

////////////////////////////////////////////////////////////////////////////////
/// Set synchronisation on or off.
///
///  \param [in] mode : synchronisation on/off
///            - mode=1  on
///            - mode<>0 off

void TGX11::Sync(int mode)
{
   switch (mode) {

      case 1 :
         XSynchronize((Display*)fDisplay,1);
         break;

      default:
         XSynchronize((Display*)fDisplay,0);
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update display.
///
///  \param [in] mode : (1) update (0) sync
///
/// Synchronise client and server once (not permanent).
/// Copy the pixmap gCws->fDrawing on the window gCws->fWindow
/// if the double buffer is on.

void TGX11::UpdateWindow(int mode)
{
   if (gCws->fDoubleBuffer) {
      XCopyArea((Display*)fDisplay, gCws->fDrawing, gCws->fWindow,
                *gGCpxmp, 0, 0, gCws->fWidth, gCws->fHeight, 0, 0);
   }
   if (mode == 1) {
      XFlush((Display*)fDisplay);
   } else {
      XSync((Display*)fDisplay, False);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set pointer position.
///
/// \param [in] ix   New X coordinate of pointer
/// \param [in] iy   New Y coordinate of pointer
/// \param [in] id   Window identifier
///
/// Coordinates are relative to the origin of the window id
/// or to the origin of the current window if id == 0.

void TGX11::Warp(Int_t ix, Int_t iy, Window_t id)
{
   if (!id) {
      // Causes problems when calling ProcessEvents()... BadWindow
      //XWarpPointer((Display*)fDisplay, None, gCws->fWindow, 0, 0, 0, 0, ix, iy);
   } else {
      XWarpPointer((Display*)fDisplay, None, (Window) id, 0, 0, 0, 0, ix, iy);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write the pixmap wid in the bitmap file pxname.
///
///  \param [in] wid         : Pixmap address
///  \param [in] w,h         : Width and height of the pixmap.
///  \param [in] pxname      : pixmap name

void TGX11::WritePixmap(int wid, unsigned int w, unsigned int h, char *pxname)
{
   unsigned int wval, hval;
   wval = w;
   hval = h;

   gTws = &fWindows[wid];
   XWriteBitmapFile((Display*)fDisplay, pxname, gTws->fDrawing, wval, hval, -1, -1);
}


//
// Functions for GIFencode()
//

static FILE *gOut;                      // output unit used WriteGIF and PutByte
static XImage *gXimage = 0;             // image used in WriteGIF and GetPixel

extern "C" {
   int GIFquantize(UInt_t width, UInt_t height, Int_t *ncol, Byte_t *red, Byte_t *green,
                   Byte_t *blue, Byte_t *outputBuf, Byte_t *outputCmap);
   long GIFencode(int Width, int Height, Int_t Ncol, Byte_t R[], Byte_t G[], Byte_t B[], Byte_t ScLine[],
                  void (*get_scline) (int, int, Byte_t *), void (*pb)(Byte_t));
   int GIFdecode(Byte_t *gifArr, Byte_t *pixArr, int *Width, int *Height, int *Ncols, Byte_t *R, Byte_t *G, Byte_t *B);
   int GIFinfo(Byte_t *gifArr, int *Width, int *Height, int *Ncols);
}

////////////////////////////////////////////////////////////////////////////////
/// Get pixels in line y and put in array scline.

static void GetPixel(int y, int width, Byte_t *scline)
{
   for (int i = 0; i < width; i++)
      scline[i] = Byte_t(XGetPixel(gXimage, i, y));
}

////////////////////////////////////////////////////////////////////////////////
/// Put byte b in output stream.

static void PutByte(Byte_t b)
{
   if (ferror(gOut) == 0) fputc(b, gOut);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns in R G B the ncol colors of the palette used by the image.
/// The image pixels are changed to index values in these R G B arrays.
/// This produces a colormap with only the used colors (so even on displays
/// with more than 8 planes we will be able to create GIF's when the image
/// contains no more than 256 different colors). If it does contain more
/// colors we will have to use GIFquantize to reduce the number of colors.
/// The R G B arrays must be deleted by the caller.

void TGX11::ImgPickPalette(RXImage *image, Int_t &ncol, Int_t *&R, Int_t *&G, Int_t *&B)
{
   ULong_t *orgcolors = 0;
   Int_t    maxcolors = 0, ncolors = 0;

   // collect different image colors
   int x, y;
   for (x = 0; x < (int) gCws->fWidth; x++) {
      for (y = 0; y < (int) gCws->fHeight; y++) {
         ULong_t pixel = XGetPixel(image, x, y);
         CollectImageColors(pixel, orgcolors, ncolors, maxcolors);
      }
   }

   // get RGB values belonging to pixels
   RXColor *xcol = new RXColor[ncolors];

   int i;
   for (i = 0; i < ncolors; i++) {
      xcol[i].pixel = orgcolors[i];
      xcol[i].red   = xcol[i].green = xcol[i].blue = 0;
      xcol[i].flags = DoRed | DoGreen | DoBlue;
   }
   QueryColors(fColormap, xcol, ncolors);

   // create RGB arrays and store RGB's for each color and set number of colors
   // (space must be delete by caller)
   R = new Int_t[ncolors];
   G = new Int_t[ncolors];
   B = new Int_t[ncolors];

   for (i = 0; i < ncolors; i++) {
      R[i] = xcol[i].red;
      G[i] = xcol[i].green;
      B[i] = xcol[i].blue;
   }
   ncol = ncolors;

   // update image with indices (pixels) into the new RGB colormap
   for (x = 0; x < (int) gCws->fWidth; x++) {
      for (y = 0; y < (int) gCws->fHeight; y++) {
         ULong_t pixel = XGetPixel(image, x, y);
         Int_t idx = FindColor(pixel, orgcolors, ncolors);
         XPutPixel(image, x, y, idx);
      }
   }

   // cleanup
   delete [] xcol;
   ::operator delete(orgcolors);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes the current window into GIF file. Returns 1 in case of success,
/// 0 otherwise.

Int_t TGX11::WriteGIF(char *name)
{
   Byte_t    scline[2000], r[256], b[256], g[256];
   Int_t    *red, *green, *blue;
   Int_t     ncol, maxcol, i;

   if (gXimage) {
      XDestroyImage(gXimage);
      gXimage = 0;
   }

   gXimage = XGetImage((Display*)fDisplay, gCws->fDrawing, 0, 0,
                       gCws->fWidth, gCws->fHeight,
                       AllPlanes, ZPixmap);

   ImgPickPalette((RXImage*)gXimage, ncol, red, green, blue);

   if (ncol > 256) {
      //GIFquantize(...);
      Error("WriteGIF", "Cannot create GIF of image containing more than 256 colors. Try in batch mode.");
      delete [] red;
      delete [] green;
      delete [] blue;
      return 0;
   }

   maxcol = 0;
   for (i = 0; i < ncol; i++) {
      if (maxcol < red[i] )   maxcol = red[i];
      if (maxcol < green[i] ) maxcol = green[i];
      if (maxcol < blue[i] )  maxcol = blue[i];
      r[i] = 0;
      g[i] = 0;
      b[i] = 0;
   }
   if (maxcol != 0) {
      for (i = 0; i < ncol; i++) {
         r[i] = red[i] * 255/maxcol;
         g[i] = green[i] * 255/maxcol;
         b[i] = blue[i] * 255/maxcol;
      }
   }

   gOut = fopen(name, "w+");

   if (gOut) {
      GIFencode(gCws->fWidth, gCws->fHeight,
             ncol, r, g, b, scline, ::GetPixel, PutByte);
      fclose(gOut);
      i = 1;
   } else {
      Error("WriteGIF","cannot write file: %s",name);
      i = 0;
   }
   delete [] red;
   delete [] green;
   delete [] blue;
   return i;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw image.

void TGX11::PutImage(int offset,int itran,int x0,int y0,int nx,int ny,int xmin,
                     int ymin,int xmax,int ymax, unsigned char *image,Drawable_t wid)
{
   const int maxSegment = 20;
   int           i, n, x, y, xcur, x1, x2, y1, y2;
   unsigned char *jimg, *jbase, icol;
   int           nlines[256];
   XSegment      lines[256][maxSegment];
   Drawable_t    id;

   if (wid) {
      id = wid;
   } else {
      id = gCws->fDrawing;
   }

   for (i = 0; i < 256; i++) nlines[i] = 0;

   x1 = x0 + xmin; y1 = y0 + ny - ymax - 1;
   x2 = x0 + xmax; y2 = y0 + ny - ymin - 1;
   jbase = image + (ymin-1)*nx + xmin;

   for (y = y2; y >= y1; y--) {
      xcur = x1; jbase += nx;
      for (jimg = jbase, icol = *jimg++, x = x1+1; x <= x2; jimg++, x++) {
         if (icol != *jimg) {
            if (icol != itran) {
               n = nlines[icol]++;
               lines[icol][n].x1 = xcur; lines[icol][n].y1 = y;
               lines[icol][n].x2 = x-1;  lines[icol][n].y2 = y;
               if (nlines[icol] == maxSegment) {
                  SetColor(gGCline,(int)icol+offset);
                  XDrawSegments((Display*)fDisplay,id,*gGCline,&lines[icol][0],
                                maxSegment);
                  nlines[icol] = 0;
               }
            }
            icol = *jimg; xcur = x;
         }
      }
      if (icol != itran) {
         n = nlines[icol]++;
         lines[icol][n].x1 = xcur; lines[icol][n].y1 = y;
         lines[icol][n].x2 = x-1;  lines[icol][n].y2 = y;
         if (nlines[icol] == maxSegment) {
            SetColor(gGCline,(int)icol+offset);
            XDrawSegments((Display*)fDisplay,id,*gGCline,&lines[icol][0],
                          maxSegment);
            nlines[icol] = 0;
         }
      }
   }

   for (i = 0; i < 256; i++) {
      if (nlines[i] != 0) {
         SetColor(gGCline,i+offset);
         XDrawSegments((Display*)fDisplay,id,*gGCline,&lines[i][0],nlines[i]);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If id is NULL - loads the specified gif file at position [x0,y0] in the
/// current window. Otherwise creates pixmap from gif file

Pixmap_t TGX11::ReadGIF(int x0, int y0, const char *file, Window_t id)
{
   FILE  *fd;
   Seek_t filesize = 0;
   unsigned char *gifArr, *pixArr, red[256], green[256], blue[256], *j1, *j2, icol;
   int   i, j, k, width, height, ncolor, irep, offset;
   float rr, gg, bb;
   Pixmap_t pic = 0;

   fd = fopen(file, "r");
   if (!fd) {
      Error("ReadGIF", "unable to open GIF file");
      return pic;
   }

   fseek(fd, 0L, 2);
   long ft = ftell(fd);
   if (ft <=0) {
      Error("ReadGIF", "unable to open GIF file");
      fclose(fd);
      return pic;
   } else {
      filesize = Seek_t(ft);
   }
   fseek(fd, 0L, 0);

   if (!(gifArr = (unsigned char *) calloc(filesize+256,1))) {
      Error("ReadGIF", "unable to allocate array for gif");
      fclose(fd);
      return pic;
   }

   if (fread(gifArr, filesize, 1, fd) != 1) {
      Error("ReadGIF", "GIF file read failed");
      free(gifArr);
      fclose(fd);
      return pic;
   }
   fclose(fd);

   irep = GIFinfo(gifArr, &width, &height, &ncolor);
   if (irep != 0) {
      free(gifArr);
      return pic;
   }

   if (!(pixArr = (unsigned char *) calloc((width*height),1))) {
      Error("ReadGIF", "unable to allocate array for image");
      free(gifArr);
      return pic;
   }

   irep = GIFdecode(gifArr, pixArr, &width, &height, &ncolor, red, green, blue);
   if (irep != 0) {
      free(gifArr);
      free(pixArr);
      return pic;
   }

   // S E T   P A L E T T E

   offset = 8;

   for (i = 0; i < ncolor; i++) {
      rr = red[i]/255.;
      gg = green[i]/255.;
      bb = blue[i]/255.;
      j = i+offset;
      SetRGB(j,rr,gg,bb);
   }

   // O U T P U T   I M A G E

   for (i = 1; i <= height/2; i++) {
      j1 = pixArr + (i-1)*width;
      j2 = pixArr + (height-i)*width;
      for (k = 0; k < width; k++) {
         icol = *j1; *j1++ = *j2; *j2++ = icol;
      }
   }
   if (id) pic = CreatePixmap(id, width, height);
   PutImage(offset,-1,x0,y0,width,height,0,0,width-1,height-1,pixArr,pic);

   free(gifArr);
   free(pixArr);

   if (pic)
      return pic;
   else if (gCws->fDrawing)
      return (Pixmap_t)gCws->fDrawing;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns an array of pixels created from a part of drawable
/// (defined by x, y, w, h) in format:
/// `b1, g1, r1, 0,  b2, g2, r2, 0,  ...,  bn, gn, rn, 0`.
///
/// Pixels are numbered from left to right and from top to bottom.
/// By default all pixels from the whole drawable are returned.
///
/// Note that return array is 32-bit aligned

unsigned char *TGX11::GetColorBits(Drawable_t /*wid*/, Int_t /*x*/, Int_t /*y*/,
                                       UInt_t /*w*/, UInt_t /*h*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// create pixmap from RGB data. RGB data is in format :
/// b1, g1, r1, 0,  b2, g2, r2, 0 ... bn, gn, rn, 0 ..
///
/// Pixels are numbered from left to right and from top to bottom.
/// Note that data must be 32-bit aligned

Pixmap_t TGX11::CreatePixmapFromData(unsigned char * /*bits*/, UInt_t /*width*/,
                                       UInt_t /*height*/)
{
   return (Pixmap_t)0;
}

////////////////////////////////////////////////////////////////////////////////
/// Register pixmap created by gVirtualGL
///
/// \param [in] pixid   Pixmap identifier
/// \param [in] w,h     Width and height of the pixmap
///
/// register new pixmap

Int_t TGX11::AddPixmap(ULong_t pixid, UInt_t w, UInt_t h)
{
   Int_t wid = 0;

   // Select next free window number
   for (; wid < fMaxNumberOfWindows; ++wid)
      if (!fWindows[wid].fOpen)
         break;

   if (wid == fMaxNumberOfWindows) {
      Int_t newsize = fMaxNumberOfWindows + 10;
      fWindows = (XWindow_t*) TStorage::ReAlloc(
                                                fWindows, newsize * sizeof(XWindow_t),
                                                fMaxNumberOfWindows*sizeof(XWindow_t)
                                               );

      for (Int_t i = fMaxNumberOfWindows; i < newsize; ++i)
         fWindows[i].fOpen = 0;

      fMaxNumberOfWindows = newsize;
   }

   fWindows[wid].fOpen = 1;
   gCws = fWindows + wid;
   gCws->fWindow = pixid;
   gCws->fDrawing = gCws->fWindow;
   gCws->fBuffer = 0;
   gCws->fDoubleBuffer = 0;
   gCws->fIsPixmap = 1;
   gCws->fClip = 0;
   gCws->fWidth = w;
   gCws->fHeight = h;
   gCws->fNewColors = 0;
   gCws->fShared = kFALSE;

   return wid;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns 1 if window system server supports extension given by the
/// argument, returns 0 in case extension is not supported and returns -1
/// in case of error (like server not initialized).
/// Examples:
///  - "Apple-WM" - does server run on MacOS X;
///  - "XINERAMA" - does server support Xinerama.
/// See also the output of xdpyinfo.

Int_t TGX11::SupportsExtension(const char *ext) const
{
   Int_t major_opcode, first_event, first_error;
   if (!(Display*)fDisplay)
      return -1;
   return XQueryExtension((Display*)fDisplay, ext, &major_opcode, &first_event, &first_error);
}
