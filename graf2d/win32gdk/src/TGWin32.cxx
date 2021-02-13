// @(#)root/win32gdk:$Id$
// Author: Rene Brun, Olivier Couet, Fons Rademakers, Valeri Onuchin, Bertrand Bellenot 27/11/01

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/// \defgroup win32 Win32 backend
/// \brief Interface to Windows graphics.
/// \ingroup GraphicsBackends

/** \class TGWin32
\ingroup win32
This class is the basic interface to the Win32 graphics system.
It is  an implementation of the abstract TVirtualX class.

This code was initially developed in the context of HIGZ and PAW
by Olivier Couet (package X11INT).
*/

#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_GLYPH_H
#include "TGWin32.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>
#include <process.h>
#include <wchar.h>
#include "gdk/gdkkeysyms.h"
#include "xatom.h"
#include <winuser.h>

#include "TROOT.h"
#include "TApplication.h"
#include "TColor.h"
#include "TPoint.h"
#include "TMath.h"
#include "TStorage.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TError.h"
#include "TException.h"
#include "TClassTable.h"
#include "KeySymbols.h"
#include "TWinNTSystem.h"
#include "TGWin32VirtualXProxy.h"
#include "TWin32SplashThread.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TExMap.h"
#include "TEnv.h"
#include "RStipples.h"
#include "GuiTypes.h"

// DND protocol version
#define XDND_PROTOCOL_VERSION   5
#ifndef IDC_HAND
#define IDC_HAND  MAKEINTRESOURCE(32649)
#endif

extern "C" {
void gdk_win32_draw_rectangle (GdkDrawable    *drawable,
                  GdkGC          *gc,
                  gint            filled,
                  gint            x,
                  gint            y,
                  gint            width,
                  gint            height);
void gdk_win32_draw_arc       (GdkDrawable    *drawable,
                  GdkGC          *gc,
                  gint            filled,
                  gint            x,
                  gint            y,
                  gint            width,
                  gint            height,
                  gint            angle1,
                  gint            angle2);
void gdk_win32_draw_polygon   (GdkDrawable    *drawable,
                  GdkGC          *gc,
                  gint            filled,
                  GdkPoint       *points,
                  gint            npoints);
void gdk_win32_draw_text      (GdkDrawable    *drawable,
                  GdkFont        *font,
                  GdkGC          *gc,
                  gint            x,
                  gint            y,
                  const gchar    *text,
                  gint            text_length);
void gdk_win32_draw_points    (GdkDrawable    *drawable,
                  GdkGC          *gc,
                  GdkPoint       *points,
                  gint            npoints);
void gdk_win32_draw_segments  (GdkDrawable    *drawable,
                  GdkGC          *gc,
                  GdkSegment     *segs,
                  gint            nsegs);
void gdk_win32_draw_lines     (GdkDrawable    *drawable,
                  GdkGC          *gc,
                  GdkPoint       *points,
                  gint            npoints);

};

//////////// internal classes & structures (very private) ////////////////

struct XWindow_t {
   Int_t    open;                 // 1 if the window is open, 0 if not
   Int_t    double_buffer;        // 1 if the double buffer is on, 0 if not
   Int_t    ispixmap;             // 1 if pixmap, 0 if not
   GdkDrawable *drawing;          // drawing area, equal to window or buffer
   GdkDrawable *window;           // win32 window
   GdkDrawable *buffer;           // pixmap used for double buffer
   UInt_t   width;                // width of the window
   UInt_t   height;               // height of the window
   Int_t    clip;                 // 1 if the clipping is on
   Int_t    xclip;                // x coordinate of the clipping rectangle
   Int_t    yclip;                // y coordinate of the clipping rectangle
   UInt_t   wclip;                // width of the clipping rectangle
   UInt_t   hclip;                // height of the clipping rectangle
   ULong_t *new_colors;           // new image colors (after processing)
   Int_t    ncolors;              // number of different colors
};


/////////////////////////////////// globals //////////////////////////////////
int gdk_debug_level;

namespace {
/////////////////////////////////// globals //////////////////////////////////

GdkAtom gClipboardAtom = GDK_NONE;
static XWindow_t *gCws;         // gCws: pointer to the current window
static XWindow_t *gTws;         // gTws: temporary pointer

//
// gColors[0]           : background also used for b/w screen
// gColors[1]           : foreground also used for b/w screen
// gColors[2..kMAXCOL-1]: colors which can be set by SetColor
//
const Int_t kBIGGEST_RGB_VALUE = 65535;
//const Int_t kMAXCOL = 1000;
//static struct {
//   Int_t defined;
//   GdkColor color;
//} gColors[kMAXCOL];

//
// Primitives Graphic Contexts global for all windows
//
const int kMAXGC = 7;
static GdkGC *gGClist[kMAXGC];
static GdkGC *gGCline;          // = gGClist[0];  // PolyLines
static GdkGC *gGCmark;          // = gGClist[1];  // PolyMarker
static GdkGC *gGCfill;          // = gGClist[2];  // Fill areas
static GdkGC *gGCtext;          // = gGClist[3];  // Text
static GdkGC *gGCinvt;          // = gGClist[4];  // Inverse text
static GdkGC *gGCdash;          // = gGClist[5];  // Dashed lines
static GdkGC *gGCpxmp;          // = gGClist[6];  // Pixmap management

static GdkGC *gGCecho;          // Input echo

static Int_t gFillHollow;       // Flag if fill style is hollow
static GdkPixmap *gFillPattern; // Fill pattern

//
// Text management
//
static const char *gTextFont = "arial.ttf";      // Current font

//
// Markers
//
const Int_t kMAXMK = 100;
static struct {
   int type;
   int n;
   GdkPoint xy[kMAXMK];
} gMarker;                      // Point list to draw marker
static int  gMarkerLineWidth = 0;
static int  gMarkerLineStyle = GDK_LINE_SOLID;
static int  gMarkerCapStyle  = GDK_CAP_ROUND;
static int  gMarkerJoinStyle = GDK_JOIN_ROUND;

//
// Keep style values for line GdkGC
//
static int  gLineWidth = 0;
static int  gLineStyle = GDK_LINE_SOLID;
static int  gCapStyle  = GDK_CAP_BUTT;
static int  gJoinStyle = GDK_JOIN_MITER;
static char gDashList[10];
static int  gDashLength = 0;
static int  gDashOffset = 0;
static int  gDashSize   = 0;

//
// Event masks
//
static ULong_t gMouseMask =
    GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK | GDK_ENTER_NOTIFY_MASK
    | GDK_LEAVE_NOTIFY_MASK | GDK_POINTER_MOTION_MASK | GDK_KEY_PRESS_MASK
    | GDK_KEY_RELEASE_MASK;
static ULong_t gKeybdMask =
    GDK_BUTTON_PRESS_MASK | GDK_KEY_PRESS_MASK | GDK_ENTER_NOTIFY_MASK |
    GDK_LEAVE_NOTIFY_MASK;

//
// Data to create an invisible cursor
//
const char null_cursor_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};
static GdkCursor *gNullCursor;

static bool gdk_initialized = false;

//---- MWM Hints stuff

struct MWMHintsProperty_t {
   Handle_t fFlags;
   Handle_t fFunctions;
   Handle_t fDecorations;
   Int_t fInputMode;
};

//---- hints

const ULong_t kMWMHintsFunctions = BIT(0);
const ULong_t kMWMHintsDecorations = BIT(1);
const ULong_t kMWMHintsInputMode = BIT(2);

const Int_t kPropMotifWMHintsElements = 4;
const Int_t kPropMWMHintElements = kPropMotifWMHintsElements;


//---- Key symbol mapping

struct KeySymbolMap_t {
   KeySym fXKeySym;
   EKeySym fKeySym;
};

static const char *keyCodeToString[] = {
   "",                          /* 0x000 */
   "",                          /* 0x001, VK_LBUTTON */
   "",                          /* 0x002, VK_RBUTTON */
   "",                          /* 0x003, VK_CANCEL */
   "",                          /* 0x004, VK_MBUTTON */
   "",                          /* 0x005 */
   "",                          /* 0x006 */
   "",                          /* 0x007 */
   "\015",                      /* 0x008, VK_BACK */
   "\t",                        /* 0x009, VK_TAB */
   "",                          /* 0x00A */
   "",                          /* 0x00B */
   "",                          /* 0x00C, VK_CLEAR */
   "\r",                        /* 0x00D, VK_RETURN */
   "",                          /* 0x00E */
   "",                          /* 0x00F */
   "",                          /* 0x010, VK_SHIFT */
   "",                          /* 0x011, VK_CONTROL */
   "",                          /* 0x012, VK_MENU */
   "",                          /* 0x013, VK_PAUSE */
   "",                          /* 0x014, VK_CAPITAL */
   "",                          /* 0x015, VK_KANA */
   "",                          /* 0x016 */
   "",                          /* 0x017 */
   "",                          /* 0x018 */
   "",                          /* 0x019, VK_KANJI */
   "",                          /* 0x01A */
   "",                          /* 0x01B, VK_ESCAPE */
   "",                          /* 0x01C, VK_CONVERT */
   "",                          /* 0x01D, VK_NONCONVERT */
   "",                          /* 0x01E */
   "",                          /* 0x01F */
   " ",                         /* 0x020, VK_SPACE */
   "",                          /* 0x021, VK_PRIOR */
   "",                          /* 0x022, VK_NEXT */
   "",                          /* 0x023, VK_END */
   "",                          /* 0x024, VK_HOME */
   "",                          /* 0x025, VK_LEFT */
   "",                          /* 0x026, VK_UP */
   "",                          /* 0x027, VK_RIGHT */
   "",                          /* 0x028, VK_DOWN */
   "",                          /* 0x029, VK_SELECT */
   "",                          /* 0x02A, VK_PRINT */
   "",                          /* 0x02B, VK_EXECUTE */
   "",                          /* 0x02C, VK_SNAPSHOT */
   "",                          /* 0x02D, VK_INSERT */
   "\037",                      /* 0x02E, VK_DELETE */
   "",                          /* 0x02F, VK_HELP */
};

//---- Mapping table of all non-trivial mappings (the ASCII keys map
//---- one to one so are not included)

static KeySymbolMap_t gKeyMap[] = {
   {GDK_Escape, kKey_Escape},
   {GDK_Tab, kKey_Tab},
#ifndef GDK_ISO_Left_Tab
   {0xFE20, kKey_Backtab},
#else
   {GDK_ISO_Left_Tab, kKey_Backtab},
#endif
   {GDK_BackSpace, kKey_Backspace},
   {GDK_Return, kKey_Return},
   {GDK_Insert, kKey_Insert},
   {GDK_Delete, kKey_Delete},
   {GDK_Clear, kKey_Delete},
   {GDK_Pause, kKey_Pause},
   {GDK_Print, kKey_Print},
   {0x1005FF60, kKey_SysReq},   // hardcoded Sun SysReq
   {0x1007ff00, kKey_SysReq},   // hardcoded X386 SysReq
   {GDK_Home, kKey_Home},       // cursor movement
   {GDK_End, kKey_End},
   {GDK_Left, kKey_Left},
   {GDK_Up, kKey_Up},
   {GDK_Right, kKey_Right},
   {GDK_Down, kKey_Down},
   {GDK_Prior, kKey_Prior},
   {GDK_Next, kKey_Next},
   {GDK_Shift_L, kKey_Shift},   // modifiers
   {GDK_Shift_R, kKey_Shift},
   {GDK_Shift_Lock, kKey_Shift},
   {GDK_Control_L, kKey_Control},
   {GDK_Control_R, kKey_Control},
   {GDK_Meta_L, kKey_Meta},
   {GDK_Meta_R, kKey_Meta},
   {GDK_Alt_L, kKey_Alt},
   {GDK_Alt_R, kKey_Alt},
   {GDK_Caps_Lock, kKey_CapsLock},
   {GDK_Num_Lock, kKey_NumLock},
   {GDK_Scroll_Lock, kKey_ScrollLock},
   {GDK_KP_Space, kKey_Space},  // numeric keypad
   {GDK_KP_Tab, kKey_Tab},
   {GDK_KP_Enter, kKey_Enter},
   {GDK_KP_Equal, kKey_Equal},
   {GDK_KP_Multiply, kKey_Asterisk},
   {GDK_KP_Add, kKey_Plus},
   {GDK_KP_Separator, kKey_Comma},
   {GDK_KP_Subtract, kKey_Minus},
   {GDK_KP_Decimal, kKey_Period},
   {GDK_KP_Divide, kKey_Slash},
   {0, (EKeySym) 0}
};


/////////////////////static auxilary functions /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static Int_t _lookup_string(Event_t * event, char *buf, Int_t buflen)
{
   int i;
   int n = event->fUser[1];
   if (n > 0) {
      for (i = 0; i < n; i++) {
         buf[i] = event->fUser[2 + i];
      }
      buf[n] = 0;
   } else {
      buf[0] = 0;
   }
   if (event->fCode <= 0x20) {
      strncpy(buf, keyCodeToString[event->fCode], buflen - 1);
   }
   return n;
}

////////////////////////////////////////////////////////////////////////////////

inline void SplitLong(Long_t ll, Long_t & i1, Long_t & i2)
{
   union {
      Long_t l;
      Int_t i[2];
   } conv;

   conv.l    = 0L;
   conv.i[0] = 0;
   conv.i[1] = 0;

   conv.l = ll;
   i1 = conv.i[0];
   i2 = conv.i[1];
}

////////////////////////////////////////////////////////////////////////////////

inline void AsmLong(Long_t i1, Long_t i2, Long_t & ll)
{
   union {
      Long_t l;
      Int_t i[2];
   } conv;

   conv.i[0] = (Int_t) i1;
   conv.i[1] = (Int_t) i2;
   ll = conv.l;
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure the child window is visible.

static BOOL CALLBACK EnumChildProc(HWND hwndChild, LPARAM lParam)
{
   ::ShowWindow(hwndChild, SW_SHOWNORMAL);
   GdkWindow *child = gdk_window_lookup(hwndChild);
   if (child)
      ((GdkWindowPrivate *) child)->mapped = TRUE;
   return TRUE;
}

////////////////////////////////////////////////////////////////////////////////

static void _ChangeProperty(HWND w, char *np, char *dp, int n, Atom_t type)
{
   HGLOBAL hMem;
   char *p;

   hMem = ::GetProp(w, np);
   if (hMem != NULL) {
      ::GlobalFree(hMem);
   }
   hMem = ::GlobalAlloc(GHND, n + sizeof(Atom_t));
   p = (char *) ::GlobalLock(hMem);
   memcpy(p, &type, sizeof(Atom_t));
   memcpy(p + sizeof(Atom_t), dp, n);
   ::GlobalUnlock(hMem);
   ::SetProp(w, np, hMem);
   ::GlobalFree(hMem);
}

////////////////////////////////////////////////////////////////////////////////
///

static void W32ChangeProperty(HWND w, Atom_t property, Atom_t type,
                       int format, int mode, const unsigned char *data,
                       int nelements)
{
   char *atomName;
   char buffer[256];
   int len;
   char propName[32];

   if (mode == GDK_PROP_MODE_REPLACE || mode == GDK_PROP_MODE_PREPEND) {
      len = (int) ::GlobalGetAtomName(property, buffer, sizeof(buffer));
      if ((atomName = (char *) malloc(len + 1)) == NULL) {
         return;
      } else {
         strcpy(atomName, buffer);
      }
      sprintf(propName, "#0x%0.4x", (unsigned) atomName);
      _ChangeProperty(w, propName, (char *) data, nelements, type);
      free(atomName);
   }
}

////////////////////////////////////////////////////////////////////////////////
///

static int _GetWindowProperty(GdkWindow * id, Atom_t property, Long_t long_offset,
                       Long_t long_length, Bool_t delete_it, Atom_t req_type,
                       Atom_t * actual_type_return,
                       Int_t * actual_format_return, ULong_t * nitems_return,
                       ULong_t * bytes_after_return, UChar_t ** prop_return)
{
   if (!id) return 0;

   char *data, *destPtr;
   char propName[32];
   HGLOBAL handle;
   HWND w;

   w = (HWND) GDK_DRAWABLE_XID(id);

   if (::IsClipboardFormatAvailable(CF_TEXT) && ::OpenClipboard(NULL)) {
      handle = ::GetClipboardData(CF_TEXT);
      if (handle != NULL) {
         data = (char *) ::GlobalLock(handle);
         *nitems_return = strlen(data);
         *prop_return = (UChar_t *) malloc(*nitems_return + 1);
         destPtr = (char *) *prop_return;
         while (*data != '\0') {
            if (*data != '\r') {
               *destPtr = *data;
               destPtr++;
            }
            data++;
         }
         *destPtr = '\0';
         ::GlobalUnlock(handle);
         *actual_type_return = XA_STRING;
         *bytes_after_return = 0;
      }
      ::CloseClipboard();
      return 1;
   }
   if (delete_it) {
      ::RemoveProp(w, propName);
   }
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
///

static ULong_t GetPixelImage(Drawable_t id, Int_t x, Int_t y)
{
   if (!id) return 0;

   GdkImage *image = (GdkImage *)id;
   ULong_t pixel;

   if (image->depth == 1) {
      pixel = (((char *) image->mem)[y * image->bpl + (x >> 3)] & (1 << (7 - (x & 0x7)))) != 0;
   } else {
      UChar_t *pixelp = (UChar_t *) image->mem + y * image->bpl + x * image->bpp;
      switch (image->bpp) {
         case 1:
            pixel = *pixelp;
            break;
         // Windows is always LSB, no need to check image->byte_order.
         case 2:
            pixel = pixelp[0] | (pixelp[1] << 8);
            break;
         case 3:
            pixel = pixelp[0] | (pixelp[1] << 8) | (pixelp[2] << 16);
            break;
         case 4:
            pixel = pixelp[0] | (pixelp[1] << 8) | (pixelp[2] << 16);
            break;
      }
   }
   return pixel;
}

////////////////////////////////////////////////////////////////////////////////
/// Collect in orgcolors all different original image colors.

static void CollectImageColors(ULong_t pixel, ULong_t * &orgcolors,
                                 Int_t & ncolors, Int_t & maxcolors)
{
   if (maxcolors == 0) {
      ncolors = 0;
      maxcolors = 100;
      orgcolors = (ULong_t*) ::operator new(maxcolors*sizeof(ULong_t));
   }

   for (int i = 0; i < ncolors; i++) {
      if (pixel == orgcolors[i]) return;
   }
   if (ncolors >= maxcolors) {
      orgcolors = (ULong_t *) TStorage::ReAlloc(orgcolors,
                                                maxcolors * 2 *
                                                sizeof(ULong_t),
                                                maxcolors *
                                                sizeof(ULong_t));
      maxcolors *= 2;
   }
   orgcolors[ncolors++] = pixel;
}

////////////////////////////////////////////////////////////////////////////////
/// debug function for printing event mask

static char *EventMask2String(UInt_t evmask)
{
   static char bfr[500];
   char *p = bfr;

   *p = '\0';
#define BITmask(x) \
  if (evmask & k##x##Mask) \
    p += sprintf (p, "%s" #x, (p > bfr ? " " : ""))
   BITmask(Exposure);
   BITmask(PointerMotion);
   BITmask(ButtonMotion);
   BITmask(ButtonPress);
   BITmask(ButtonRelease);
   BITmask(KeyPress);
   BITmask(KeyRelease);
   BITmask(EnterWindow);
   BITmask(LeaveWindow);
   BITmask(FocusChange);
   BITmask(StructureNotify);
#undef BITmask

   return bfr;
}

///////////////////////////////////////////////////////////////////////////////
class TGWin32MainThread {

public:
   void     *fHandle;                     // handle of server (aka command) thread
   DWORD    fId;                          // id of server (aka command) thread
   static LPCRITICAL_SECTION  fCritSec;      // general mutex
   static LPCRITICAL_SECTION  fMessageMutex; // message queue mutex

   TGWin32MainThread();
   ~TGWin32MainThread();
   static void LockMSG();
   static void UnlockMSG();
};

TGWin32MainThread *gMainThread = 0;
LPCRITICAL_SECTION TGWin32MainThread::fCritSec = 0;
LPCRITICAL_SECTION TGWin32MainThread::fMessageMutex = 0;


////////////////////////////////////////////////////////////////////////////////
/// dtor

TGWin32MainThread::~TGWin32MainThread()
{
   if (fCritSec) {
      ::LeaveCriticalSection(fCritSec);
      ::DeleteCriticalSection(fCritSec);
      delete fCritSec;
   }
   fCritSec = 0;

   if (fMessageMutex) {
      ::LeaveCriticalSection(fMessageMutex);
      ::DeleteCriticalSection(fMessageMutex);
      delete fMessageMutex;
   }
   fMessageMutex = 0;

   if(fHandle) {
      ::PostThreadMessage(fId, WM_QUIT, 0, 0);
      ::CloseHandle(fHandle);
   }
   fHandle = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// lock message queue

void TGWin32MainThread::LockMSG()
{
   if (fMessageMutex) ::EnterCriticalSection(fMessageMutex);
}

////////////////////////////////////////////////////////////////////////////////
/// unlock message queue

void TGWin32MainThread::UnlockMSG()
{
   if (fMessageMutex) ::LeaveCriticalSection(fMessageMutex);
}


////////////////////////////////////////////////////////////////////////////////
/// Windows timer handling events while moving/resizing windows

VOID CALLBACK MyTimerProc(HWND hwnd, UINT message, UINT idTimer, DWORD dwTime)
{
   gSystem->ProcessEvents();
   //gVirtualX->UpdateWindow(1); // cause problems with OpenGL in pad...
}

////////////////////////////////////////////////////////////////////////////////
/// Message processing function for the GUI thread.
/// Kicks in once TGWin32 becomes active, and "replaces" the dummy one
/// in TWinNTSystem; see TWinNTSystem.cxx's GUIThreadMessageProcessingLoop().

Bool_t GUIThreadMessageFunc(MSG *msg)
{
   Bool_t ret = kFALSE;
   static Int_t m_timer = 0;

   if ( (msg->message == WM_NCLBUTTONDOWN) ) {
      if (m_timer == 0)
         m_timer = SetTimer(NULL, 1, 20, (TIMERPROC) MyTimerProc);
   }
   else if (msg->message == WM_NCMOUSELEAVE ) {
      if (m_timer) {
         KillTimer(NULL, m_timer);
      }
      m_timer = 0;
   }

   if (msg->message == TGWin32ProxyBase::fgPostMessageId) {
      if (msg->wParam) {
         TGWin32ProxyBase *proxy = (TGWin32ProxyBase*)msg->wParam;
         proxy->ExecuteCallBack(kTRUE);
      } else {
         ret = kTRUE;
      }
   } else if (msg->message == TGWin32ProxyBase::fgPingMessageId) {
      TGWin32ProxyBase::GlobalUnlock();
   } else {
      //if ( (msg->message >= WM_NCMOUSEMOVE) &&
      //     (msg->message <= WM_NCMBUTTONDBLCLK) ) {
      //   TGWin32ProxyBase::GlobalLock();
      //}
      TGWin32MainThread::LockMSG();
      TranslateMessage(msg);
      DispatchMessage(msg);
      TGWin32MainThread::UnlockMSG();
   }
   return ret;
}


///////////////////////////////////////////////////////////////////////////////
class TGWin32RefreshTimer : public TTimer {

public:
   TGWin32RefreshTimer() : TTimer(10, kTRUE) { if (gSystem) gSystem->AddTimer(this); }
   ~TGWin32RefreshTimer() { if (gSystem) gSystem->RemoveTimer(this); }
   Bool_t Notify()
   {
      Reset();
      MSG msg;

      while (::PeekMessage(&msg, NULL, NULL, NULL, PM_NOREMOVE)) {
         ::PeekMessage(&msg, NULL, NULL, NULL, PM_REMOVE);
         if (!gVirtualX)
            Sleep(200); // avoid start-up race
         if (gVirtualX)
            GUIThreadMessageFunc(&msg);
      }
      return kFALSE;
   }
};
/*
////////////////////////////////////////////////////////////////////////////////
/// thread for processing windows messages (aka Main/Server thread)

static DWORD WINAPI MessageProcessingLoop(void *p)
{
   MSG msg;
   Int_t erret;
   Bool_t endLoop = kFALSE;
   TGWin32RefreshTimer *refersh = 0;

   // force to create message queue
   ::PeekMessage(&msg, NULL, WM_USER, WM_USER, PM_NOREMOVE);

   // periodically we refresh windows
   // Don't create refresh timer if the application has been created inside PVSS
   if (gApplication) {
      TString arg = gSystem->BaseName(gApplication->Argv(0));
      if (!arg.Contains("PVSS"))
         refersh = new TGWin32RefreshTimer();
   }

   while (!endLoop) {
      erret = ::GetMessage(&msg, NULL, NULL, NULL);
      if (erret <= 0) endLoop = kTRUE;
      endLoop = MessageProcessingFunc(&msg);
   }

   TGWin32::Instance()->CloseDisplay();
   if (refersh)
      delete refersh;

   // exit thread
   if (erret == -1) {
      erret = ::GetLastError();
      Error("MsgLoop", "Error in GetMessage");
      ::ExitThread(-1);
   } else {
      ::ExitThread(0);
   }
   return 0;
}
*/

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGWin32MainThread::TGWin32MainThread()
{
   fCritSec = new CRITICAL_SECTION;
   ::InitializeCriticalSection(fCritSec);
   fMessageMutex = new CRITICAL_SECTION;
   ::InitializeCriticalSection(fMessageMutex);
   fHandle = ((TWinNTSystem*)gSystem)->GetGUIThreadHandle();
   fId = ((TWinNTSystem*)gSystem)->GetGUIThreadId();
   ((TWinNTSystem*)gSystem)->SetGUIThreadMsgHandler(GUIThreadMessageFunc);
}

} // unnamed namespace

///////////////////////// TGWin32 implementation ///////////////////////////////
ClassImp(TGWin32);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TGWin32::TGWin32(): fRefreshTimer(0)
{
   fScreenNumber = 0;
   fWindows      = 0;
   fColors       = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal Constructor.

TGWin32::TGWin32(const char *name, const char *title) : TVirtualX(name,title), fRefreshTimer(0)
{
   fScreenNumber = 0;
   fHasTTFonts = kFALSE;
   fUseSysPointers = kFALSE;
   fTextAlignH = 1;
   fTextAlignV = 1;
   fTextAlign = 7;
   fTextMagnitude = 1;
   fCharacterUpX = 1;
   fCharacterUpY = 1;
   fDrawMode = kCopy;
   fWindows = 0;
   fMaxNumberOfWindows = 10;
   fXEvent = 0;
   fFillColorModified = kFALSE;
   fFillStyleModified = kFALSE;
   fLineColorModified = kFALSE;
   fPenModified = kFALSE;
   fMarkerStyleModified = kFALSE;
   fMarkerColorModified = kFALSE;

   fWindows = (XWindow_t*) TStorage::Alloc(fMaxNumberOfWindows*sizeof(XWindow_t));
   for (int i = 0; i < fMaxNumberOfWindows; i++) fWindows[i].open = 0;

   fColors = new TExMap;

   if (gApplication) {
      TString arg = gSystem->BaseName(gApplication->Argv(0));
      if (!arg.Contains("PVSS"))
         fRefreshTimer = new TGWin32RefreshTimer();
   } else {
      fRefreshTimer = new TGWin32RefreshTimer();
   }

   // initialize GUI thread and proxy objects
   if (!gROOT->IsBatch() && !gMainThread) {
      gMainThread = new TGWin32MainThread();
      TGWin32ProxyBase::fgMainThreadId = ::GetCurrentThreadId(); // gMainThread->fId;
      TGWin32VirtualXProxy::fgRealObject = this;
      gPtr2VirtualX = &TGWin32VirtualXProxy::ProxyObject;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// destructor.

TGWin32::~TGWin32()
{
   CloseDisplay();
   if (fRefreshTimer)
      delete fRefreshTimer;
   if (!fColors) return;
   Long64_t key, value;
   TExMapIter it(fColors);
   while (it.Next(key, value)) {
      XColor_t *col = (XColor_t *) value;
      delete col;
   }
   delete fColors;
}

////////////////////////////////////////////////////////////////////////////////
/// returns kTRUE if we are inside cmd/server thread

Bool_t TGWin32::IsCmdThread() const
{
#ifdef OLD_THREAD_IMPLEMENTATION
   return ((::GetCurrentThreadId() == TGWin32ProxyBase::fgMainThreadId) ||
           (::GetCurrentThreadId() == TGWin32ProxyBase::fgUserThreadId));
#else
   return kTRUE;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// close display (terminate server/gMainThread thread)

void TGWin32::CloseDisplay()
{
   // disable any processing while exiting
   TGWin32ProxyBase::GlobalLock();

   // terminate server thread
   gPtr2VirtualX = 0;
   gVirtualX = TGWin32VirtualXProxy::RealObject();

   // The lock above does not work, so at least
   // minimize the risk
   TGWin32MainThread *delThread = gMainThread;
   if (gMainThread) {
      gMainThread = 0;
      delete delThread;
   }

   TGWin32ProxyBase::fgMainThreadId = 0;

   // terminate ROOT logo splash thread
   TWin32SplashThread *delSplash = gSplash;
   if (gSplash) {
      gSplash = 0;
      delete delSplash;
   }

   if (fWindows) TStorage::Dealloc(fWindows);
   fWindows = 0;

   if (fXEvent) gdk_event_free((GdkEvent*)fXEvent);

   TGWin32ProxyBase::GlobalUnlock();

   gROOT->SetBatch(kTRUE); // no GUI is possible
}

////////////////////////////////////////////////////////////////////////////////
///

void  TGWin32::Lock()
{
   if (gMainThread && gMainThread->fCritSec) ::EnterCriticalSection(gMainThread->fCritSec);
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::Unlock()
{
   if (gMainThread && gMainThread->fCritSec) ::LeaveCriticalSection(gMainThread->fCritSec);
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize Win32 system. Returns kFALSE in case of failure.

Bool_t TGWin32::Init(void *display)
{
   if (!gdk_initialized) {
      if (!gdk_init_check(NULL, NULL)) return kFALSE;
      gdk_initialized = true;
   }

   if (!gClipboardAtom) {
      gClipboardAtom = gdk_atom_intern("CLIPBOARD", kFALSE);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Open the display. Return -1 if the opening fails, 0 when ok.

Int_t TGWin32::OpenDisplay(const char *dpyName)
{
   GdkPixmap *pixmp1, *pixmp2;
   GdkColor fore, back;
   GdkColor color;
   GdkGCValues gcvals;
   int i;

   if (!Init((void*)dpyName)) {
      return -1;
   }

   if (gDebug <= 4) {
      gdk_debug_level = gDebug;
   } else {
      gdk_debug_level = 0;
   }

   fore.red = fore.green = fore.blue = 0;
   back.red = back.green = back.blue = 0;
   color.red = color.green = color.blue = 0;

   fScreenNumber = 0;           //DefaultScreen(fDisplay);
   fVisual = gdk_visual_get_best();
   fColormap = gdk_colormap_get_system();
   fDepth = gdk_visual_get_best_depth();

   GetColor(1).fDefined = kTRUE; // default foreground
   gdk_color_black((GdkColormap *)fColormap, &GetColor(1).color);

   GetColor(0).fDefined = kTRUE; // default background
   gdk_color_white((GdkColormap *)fColormap, &GetColor(0).color);

   // Create primitives graphic contexts
   for (i = 0; i < kMAXGC; i++) {
      gGClist[i]  = gdk_gc_new(GDK_ROOT_PARENT());
      gdk_gc_set_foreground(gGClist[i], &GetColor(1).color);
      gdk_gc_set_background(gGClist[i], &GetColor(0).color);
   }

   gGCline = gGClist[0];        // PolyLines
   gGCmark = gGClist[1];        // PolyMarker
   gGCfill = gGClist[2];        // Fill areas
   gGCtext = gGClist[3];        // Text
   gGCinvt = gGClist[4];        // Inverse text
   gGCdash = gGClist[5];        // Dashed lines
   gGCpxmp = gGClist[6];        // Pixmap management

   gdk_gc_get_values(gGCtext, &gcvals);
   gdk_gc_set_foreground(gGCinvt, &gcvals.background);
   gdk_gc_set_background(gGCinvt, &gcvals.foreground);

   // Create input echo graphic context
   GdkGCValues echov;
   gdk_color_black(fColormap, &echov.foreground); // = BlackPixel(fDisplay, fScreenNumber);
   gdk_color_white(fColormap, &echov.background); // = WhitePixel(fDisplay, fScreenNumber);
   echov.function = GDK_INVERT;
   echov.subwindow_mode = GDK_CLIP_BY_CHILDREN;
   gGCecho =
       gdk_gc_new_with_values((GdkWindow *) GDK_ROOT_PARENT(), &echov,
                              (GdkGCValuesMask) (GDK_GC_FOREGROUND |
                                                 GDK_GC_BACKGROUND |
                                                 GDK_GC_FUNCTION |
                                                 GDK_GC_SUBWINDOW));
   // Create a null cursor
   pixmp1 = gdk_bitmap_create_from_data(GDK_ROOT_PARENT(),
                                       (const char *)null_cursor_bits, 16,16);

   pixmp2 = gdk_bitmap_create_from_data(GDK_ROOT_PARENT(),
                                       (const char *)null_cursor_bits, 16, 16);

   gNullCursor = gdk_cursor_new_from_pixmap((GdkDrawable *)pixmp1, (GdkDrawable *)pixmp2,
                                             &fore, &back, 0, 0);
   // Create cursors
   if (gEnv->GetValue("Win32.UseSysPointers", 0)) {
      fUseSysPointers = kTRUE;
      fCursors[kBottomLeft] = gdk_syscursor_new((ULong_t)IDC_SIZENESW);
      fCursors[kBottomRight] = gdk_syscursor_new((ULong_t)IDC_SIZENWSE);
      fCursors[kTopLeft] = gdk_syscursor_new((ULong_t)IDC_SIZENWSE);
      fCursors[kTopRight] = gdk_syscursor_new((ULong_t)IDC_SIZENESW);
      fCursors[kBottomSide] =  gdk_syscursor_new((ULong_t)IDC_SIZENS);
      fCursors[kLeftSide] = gdk_syscursor_new((ULong_t)IDC_SIZEWE);
      fCursors[kTopSide] = gdk_syscursor_new((ULong_t)IDC_SIZENS);
      fCursors[kRightSide] = gdk_syscursor_new((ULong_t)IDC_SIZEWE);
      fCursors[kMove] = gdk_syscursor_new((ULong_t)IDC_SIZEALL);
      fCursors[kCross] =gdk_syscursor_new((ULong_t)IDC_CROSS);
      fCursors[kArrowHor] = gdk_syscursor_new((ULong_t)IDC_SIZEWE);
      fCursors[kArrowVer] = gdk_syscursor_new((ULong_t)IDC_SIZENS);
      fCursors[kHand] = gdk_syscursor_new((ULong_t)IDC_HAND);
      fCursors[kPointer] = gdk_syscursor_new((ULong_t)IDC_ARROW);
      fCursors[kCaret] =  gdk_syscursor_new((ULong_t)IDC_IBEAM);
      fCursors[kWatch] = gdk_syscursor_new((ULong_t)IDC_WAIT);
      fCursors[kNoDrop] = gdk_syscursor_new((ULong_t)IDC_NO);
   }
   else {
      fUseSysPointers = kFALSE;
      fCursors[kBottomLeft] = gdk_cursor_new(GDK_BOTTOM_LEFT_CORNER);
      fCursors[kBottomRight] = gdk_cursor_new(GDK_BOTTOM_RIGHT_CORNER);
      fCursors[kTopLeft] = gdk_cursor_new(GDK_TOP_LEFT_CORNER);
      fCursors[kTopRight] = gdk_cursor_new(GDK_TOP_RIGHT_CORNER);
      fCursors[kBottomSide] =  gdk_cursor_new(GDK_BOTTOM_SIDE);
      fCursors[kLeftSide] = gdk_cursor_new(GDK_LEFT_SIDE);
      fCursors[kTopSide] = gdk_cursor_new(GDK_TOP_SIDE);
      fCursors[kRightSide] = gdk_cursor_new(GDK_RIGHT_SIDE);
      fCursors[kMove] = gdk_cursor_new(GDK_FLEUR);
      fCursors[kCross] =gdk_cursor_new(GDK_CROSSHAIR);
      fCursors[kArrowHor] = gdk_cursor_new(GDK_SB_H_DOUBLE_ARROW);
      fCursors[kArrowVer] = gdk_cursor_new(GDK_SB_V_DOUBLE_ARROW);
      fCursors[kHand] = gdk_cursor_new(GDK_HAND2);
      fCursors[kPointer] = gdk_cursor_new(GDK_LEFT_PTR);
      fCursors[kCaret] =  gdk_cursor_new(GDK_XTERM);
      //fCursors[kWatch] = gdk_cursor_new(GDK_WATCH);
      fCursors[kWatch] = gdk_cursor_new(GDK_BUSY);
      fCursors[kNoDrop] = gdk_cursor_new(GDK_PIRATE);
   }
   fCursors[kRotate] = gdk_cursor_new(GDK_EXCHANGE);
   fCursors[kArrowRight] = gdk_cursor_new(GDK_ARROW);

   // Setup color information
   fRedDiv = fGreenDiv = fBlueDiv = fRedShift = fGreenShift = fBlueShift = -1;

   if ( gdk_visual_get_best_type() == GDK_VISUAL_TRUE_COLOR) {
      int i;
      for (i = 0; i < int(sizeof(fVisual->blue_mask)*kBitsPerByte); i++) {
         if (fBlueShift == -1 && ((fVisual->blue_mask >> i) & 1)) {
            fBlueShift = i;
         }
         if ((fVisual->blue_mask >> i) == 1) {
            fBlueDiv = sizeof(UShort_t)*kBitsPerByte - i - 1 + fBlueShift;
            break;
         }
      }
      for (i = 0; i < int(sizeof(fVisual->green_mask)*kBitsPerByte); i++) {
         if (fGreenShift == -1 && ((fVisual->green_mask >> i) & 1)) {
            fGreenShift = i;
         }
         if ((fVisual->green_mask >> i) == 1) {
            fGreenDiv = sizeof(UShort_t)*kBitsPerByte - i - 1 + fGreenShift;
            break;
         }
      }
      for (i = 0; i < int(sizeof(fVisual->red_mask)*kBitsPerByte); i++) {
         if (fRedShift == -1 && ((fVisual->red_mask >> i) & 1)) {
            fRedShift = i;
         }
         if ((fVisual->red_mask >> i) == 1) {
            fRedDiv = sizeof(UShort_t)*kBitsPerByte - i - 1 + fRedShift;
            break;
         }
      }
   }

   SetName("Win32TTF");
   SetTitle("ROOT interface to Win32 with TrueType fonts");

   if (!TTF::IsInitialized()) TTF::Init();

   if (fDepth > 8) {
      TTF::SetSmoothing(kTRUE);
   } else {
      TTF::SetSmoothing(kFALSE);
   }

   TGWin32VirtualXProxy::fMaxResponseTime = 1000;
   fHasTTFonts = kTRUE;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Allocate color in colormap. If we are on an <= 8 plane machine
/// we will use XAllocColor. If we are on a >= 15 (15, 16 or 24) plane
/// true color machine we will calculate the pixel value using:
/// for 15 and 16 bit true colors have 6 bits precision per color however
/// only the 5 most significant bits are used in the color index.
/// Except for 16 bits where green uses all 6 bits. I.e.:
///   15 bits = rrrrrgggggbbbbb
///   16 bits = rrrrrggggggbbbbb
/// for 24 bits each r, g and b are represented by 8 bits.
///
/// Since all colors are set with a max of 65535 (16 bits) per r, g, b
/// we just right shift them by 10, 11 and 10 bits for 16 planes, and
/// (10, 10, 10 for 15 planes) and by 8 bits for 24 planes.
/// Returns kFALSE in case color allocation failed.

Bool_t TGWin32::AllocColor(GdkColormap *cmap, GdkColor *color)
{
   if (fRedDiv == -1) {
      if ( gdk_color_alloc((GdkColormap *)cmap, (GdkColor *)color) ) return kTRUE;
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

void TGWin32::QueryColors(GdkColormap *cmap, GdkColor *color, Int_t ncolors)
{
   ULong_t r, g, b;

   if (fRedDiv == -1) {
      GdkColorContext *cc = gdk_color_context_new(gdk_visual_get_system(), cmap);
      gdk_color_context_query_colors(cc, color, ncolors);
      gdk_color_context_free(cc);
   } else {
      for (Int_t i = 0; i < ncolors; i++) {
         r = (color[i].pixel & fVisual->red_mask) >> fRedShift;
         color[i].red = UShort_t(r*kBIGGEST_RGB_VALUE/(fVisual->red_mask >> fRedShift));

         g = (color[i].pixel & fVisual->green_mask) >> fGreenShift;
         color[i].green = UShort_t(g*kBIGGEST_RGB_VALUE/(fVisual->green_mask >> fGreenShift));

         b = (color[i].pixel & fVisual->blue_mask) >> fBlueShift;
         color[i].blue = UShort_t(b*kBIGGEST_RGB_VALUE/(fVisual->blue_mask >> fBlueShift));
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute alignment variables. The alignment is done on the horizontal string
/// then the rotation is applied on the alignment variables.
/// SetRotation and LayoutGlyphs should have been called before.

void TGWin32::Align(void)
{
   EAlign align = (EAlign) fTextAlign;

   // vertical alignment
   if (align == kTLeft || align == kTCenter || align == kTRight) {
      fAlign.y = TTF::GetAscent();
   } else if (align == kMLeft || align == kMCenter || align == kMRight) {
      fAlign.y = TTF::GetAscent()/2;
   } else {
      fAlign.y = 0;
   }
   // horizontal alignment
   if (align == kTRight || align == kMRight || align == kBRight) {
      fAlign.x = TTF::GetWidth();
   } else if (align == kTCenter || align == kMCenter || align == kBCenter) {
      fAlign.x = TTF::GetWidth()/2;
   } else {
      fAlign.x = 0;
   }

   FT_Vector_Transform(&fAlign, TTF::GetRotMatrix());
   fAlign.x = fAlign.x >> 6;
   fAlign.y = fAlign.y >> 6;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw FT_Bitmap bitmap to xim image at position bx,by using specified
/// foreground color.

void TGWin32::DrawImage(FT_Bitmap *source, ULong_t fore, ULong_t back,
                         GdkImage *xim, Int_t bx, Int_t by)
{
   UChar_t d = 0, *s = source->buffer;

   if (TTF::GetSmoothing()) {

      static GdkColor col[5];
      GdkColor *bcol = 0, *bc;
      Int_t    x, y;

      // background kClear, i.e. transparent, we take as background color
      // the average of the rgb values of all pixels covered by this character
      if (back == (ULong_t) -1 && (UInt_t)source->width) {
         ULong_t r, g, b;
         Int_t   dots, dotcnt;
         const Int_t maxdots = 50000;

         dots = Int_t(source->width * source->rows);
         dots = dots > maxdots ? maxdots : dots;
         bcol = new GdkColor[dots];
         if (!bcol) return;

         bc = bcol;
         dotcnt = 0;
         for (y = 0; y < (int) source->rows; y++) {
            for (x = 0; x < (int) source->width; x++, bc++) {
               bc->pixel = GetPixelImage((Drawable_t)xim, bx + x, by + y);
               if (++dotcnt >= maxdots) break;
            }
         }
         QueryColors(fColormap, bcol, dots);
         r = g = b = 0;
         bc = bcol;
         dotcnt = 0;
         for (y = 0; y < (int) source->rows; y++) {
            for (x = 0; x < (int) source->width; x++, bc++) {
               r += bc->red;
               g += bc->green;
               b += bc->blue;
               if (++dotcnt >= maxdots) break;
            }
         }
         if (dots != 0) {
            r /= dots;
            g /= dots;
            b /= dots;
         }
         bc = &col[0];
         if (bc->red == r && bc->green == g && bc->blue == b) {
            bc->pixel = back;
         } else {
            bc->pixel = ~back;
            bc->red   = (UShort_t) r;
            bc->green = (UShort_t) g;
            bc->blue  = (UShort_t) b;
         }
      }
      delete [] bcol;

      // if fore or background have changed from previous character
      // recalculate the 3 smooting colors (interpolation between fore-
      // and background colors)
      if (fore != col[4].pixel || back != col[0].pixel) {
         col[4].pixel = fore;
         if (back != (ULong_t) -1) {
            col[3].pixel = back;
            QueryColors(fColormap, &col[3], 2);
            col[0] = col[3];
         } else {
            QueryColors(fColormap, &col[4], 1);
         }

         // interpolate between fore and backgound colors
         for (x = 3; x > 0; x--) {
            col[x].red   = (col[4].red  *x + col[0].red  *(4-x)) /4;
            col[x].green = (col[4].green*x + col[0].green*(4-x)) /4;
            col[x].blue  = (col[4].blue *x + col[0].blue *(4-x)) /4;
            if (!AllocColor(fColormap, &col[x])) {
               Warning("DrawImage", "cannot allocate smoothing color");
               col[x].pixel = col[x+1].pixel;
            }
         }
      }

      // put smoothed character, character pixmap values are an index
      // into the 5 colors used for aliasing (4 = foreground, 0 = background)
      for (y = 0; y < (int) source->rows; y++) {
         for (x = 0; x < (int) source->width; x++) {
            d = *s++ & 0xff;
            d = ((d + 10) * 5) / 256;
            if (d > 4) d = 4;
            if (d && x < (int) source->width) {
               ULong_t p = col[d].pixel;
               PutPixel((Drawable_t)xim, bx + x, by + y, p);
            }
         }
      }
   } else {
      // no smoothing, just put character using foreground color
      UChar_t* row=s;
      for (int y = 0; y < (int) source->rows; y++) {
         int n = 0;
         s = row;
         for (int x = 0; x < (int) source->width; x++) {
            if (n == 0) d = *s++;
            if (TESTBIT(d,7-n)) {
               PutPixel((Drawable_t)xim, bx + x, by + y, fore);
            }
            if (++n == (int) kBitsPerByte) n = 0;
         }
         row += source->pitch;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw text using TrueType fonts. If TrueType fonts are not available the
/// text is drawn with TGWin32::DrawText.

void TGWin32::DrawText(Int_t x, Int_t y, Float_t angle, Float_t mgn,
                       const char *text, ETextMode mode)
{
   if (!TTF::IsInitialized()) TTF::Init();
   TTF::SetRotationMatrix(angle);
   TTF::PrepareString(text);
   TTF::LayoutGlyphs();
   Align();
   RenderString(x, y, mode);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw text using TrueType fonts. If TrueType fonts are not available the
/// text is drawn with TGWin32::DrawText.

void TGWin32::DrawText(Int_t x, Int_t y, Float_t angle, Float_t mgn,
                       const wchar_t *text, ETextMode mode)
{
   if (!TTF::IsInitialized()) TTF::Init();
   TTF::SetRotationMatrix(angle);
   TTF::PrepareString(text);
   TTF::LayoutGlyphs();
   Align();
   RenderString(x, y, mode);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the background of the current window in an XImage.

GdkImage *TGWin32::GetBackground(Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   Window_t cws = GetCurrentWindow();
   UInt_t width;
   UInt_t height;
   Int_t xy;
   gVirtualX->GetWindowSize(cws, xy, xy, width, height);

   if (x < 0) {
      w += x;
      x  = 0;
   }
   if (y < 0) {
      h += y;
      y  = 0;
   }

   if (x+w > width)  w = width - x;
   if (y+h > height) h = height - y;

   return gdk_image_get((GdkDrawable*)cws, x, y, w, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Test if there is really something to render

Bool_t TGWin32::IsVisible(Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   Window_t cws = GetCurrentWindow();
   UInt_t width;
   UInt_t height;
   Int_t xy;
   gVirtualX->GetWindowSize(cws, xy, xy, width, height);

   // If w or h is 0, very likely the string is only blank characters
   if ((int)w == 0 || (int)h == 0)  return kFALSE;

   // If string falls outside window, there is probably no need to draw it.
   if (x + (int)w <= 0 || x >= (int)width)  return kFALSE;
   if (y + (int)h <= 0 || y >= (int)height) return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform the string rendering in the pad.
/// LayoutGlyphs should have been called before.

void TGWin32::RenderString(Int_t x, Int_t y, ETextMode mode)
{
   TTF::TTGlyph* glyph = TTF::GetGlyphs();
   GdkGCValues gcvals;

   // compute the size and position of the XImage that will contain the text
   Int_t Xoff = 0; if (TTF::GetBox().xMin < 0) Xoff = -TTF::GetBox().xMin;
   Int_t Yoff = 0; if (TTF::GetBox().yMin < 0) Yoff = -TTF::GetBox().yMin;
   Int_t w    = TTF::GetBox().xMax + Xoff;
   Int_t h    = TTF::GetBox().yMax + Yoff;
   Int_t x1   = x-Xoff-fAlign.x;
   Int_t y1   = y+Yoff+fAlign.y-h;

   if (!IsVisible(x1, y1, w, h)) {
       return;
   }

   // create the XImage that will contain the text
   UInt_t depth = fDepth;
   GdkImage *xim  = gdk_image_new(GDK_IMAGE_SHARED, gdk_visual_get_best(), w, h);

   // use malloc since Xlib will use free() in XDestroyImage
//   xim->data = (char *) malloc(xim->bytes_per_line * h);
//   memset(xim->data, 0, xim->bytes_per_line * h);

   ULong_t   pixel;
   ULong_t   bg;

   gdk_gc_get_values((GdkGC*)GetGC(3), &gcvals);

   // get the background
   if (mode == kClear) {
      // if mode == kClear we need to get an image of the background
      GdkImage *bim = GetBackground(x1, y1, w, h);
      if (!bim) {
         Error("DrawText", "error getting background image");
         return;
      }

      // and copy it into the text image
      Int_t xo = 0, yo = 0;
      if (x1 < 0) xo = -x1;
      if (y1 < 0) yo = -y1;

      for (int yp = 0; yp < (int) bim->height; yp++) {
         for (int xp = 0; xp < (int) bim->width; xp++) {
            pixel = GetPixelImage((Drawable_t)bim, xp, yp);
            PutPixel((Drawable_t)xim, xo+xp, yo+yp, pixel);
         }
      }

      gdk_image_unref((GdkImage *)bim);

      bg = (ULong_t) -1;
   } else {
      // if mode == kOpaque its simple, we just draw the background

      GdkImage *bim = GetBackground(x1, y1, w, h);
      if (!bim) {
         pixel = gcvals.background.pixel;
      } else {
         pixel = GetPixelImage((Drawable_t)bim, 0, 0);
      }
      Int_t xo = 0, yo = 0;
      if (x1 < 0) xo = -x1;
      if (y1 < 0) yo = -y1;

      for (int yp = 0; yp < h; yp++) {
         for (int xp = 0; xp < (int) w; xp++) {
            PutPixel((Drawable_t)xim, xo+xp, yo+yp, pixel);
         }
      }
      if (bim) {
         gdk_image_unref((GdkImage *)bim);
         bg = (ULong_t) -1;
      } else {
         bg = pixel;
      }
   }

   // paint the glyphs in the XImage
   glyph = TTF::GetGlyphs();
   for (int n = 0; n < TTF::GetNumGlyphs(); n++, glyph++) {
      if (FT_Glyph_To_Bitmap(&glyph->fImage,
                             TTF::GetSmoothing() ? ft_render_mode_normal
                                              : ft_render_mode_mono,
                             0, 1 )) continue;
      FT_BitmapGlyph bitmap = (FT_BitmapGlyph)glyph->fImage;
      FT_Bitmap*     source = &bitmap->bitmap;
      Int_t          bx, by;

      bx = bitmap->left+Xoff;
      by = h - bitmap->top-Yoff;
      DrawImage(source, gcvals.foreground.pixel, bg, xim, bx, by);
   }

   // put the Ximage on the screen
   Window_t cws = GetCurrentWindow();
   gdk_draw_image((GdkDrawable *)cws, GetGC(6), xim, 0, 0, x1, y1, w, h);

   gdk_image_unref(xim);
}

////////////////////////////////////////////////////////////////////////////////
/// Set specified font.

void TGWin32::SetTextFont(Font_t fontnumber)
{
   fTextFont = fontnumber;
   TTF::SetTextFont(fontnumber);
}

////////////////////////////////////////////////////////////////////////////////
/// Set text font to specified name.
/// mode       : loading flag
/// mode=0     : search if the font exist (kCheck)
/// mode=1     : search the font and load it if it exists (kLoad)
/// font       : font name
///
/// Set text font to specified name. This function returns 0 if
/// the specified font is found, 1 if not.

Int_t TGWin32::SetTextFont(char *fontname, ETextSetMode mode)
{
   return TTF::SetTextFont(fontname);
}

////////////////////////////////////////////////////////////////////////////////
/// Set current text size.

void TGWin32::SetTextSize(Float_t textsize)
{
   fTextSize = textsize;
   TTF::SetTextSize(textsize);
}

////////////////////////////////////////////////////////////////////////////////
/// Clear current window.

void TGWin32::ClearWindow()
{
   if (!fWindows) return;

   if (!gCws->ispixmap && !gCws->double_buffer) {
      gdk_window_set_background(gCws->drawing, (GdkColor *) & GetColor(0).color);
      gdk_window_clear(gCws->drawing);
      GdiFlush();
   } else {
      SetColor(gGCpxmp, 0);
      gdk_win32_draw_rectangle(gCws->drawing, gGCpxmp, 1,
                         0, 0, gCws->width, gCws->height);
      SetColor(gGCpxmp, 1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete current pixmap.

void TGWin32::ClosePixmap()
{
   CloseWindow1();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete current window.

void TGWin32::CloseWindow()
{
   CloseWindow1();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete current window.

void TGWin32::CloseWindow1()
{
   int wid;

   if (gCws->ispixmap) {
      gdk_pixmap_unref(gCws->window);
   } else {
      gdk_window_destroy(gCws->window, kTRUE);
   }

   if (gCws->buffer) {
      gdk_pixmap_unref(gCws->buffer);
   }
   if (gCws->new_colors) {
      gdk_colormap_free_colors((GdkColormap *) fColormap,
                               (GdkColor *)gCws->new_colors, gCws->ncolors);

      delete [] gCws->new_colors;
      gCws->new_colors = 0;
   }

   GdiFlush();
   gCws->open = 0;

   if (!fWindows) return;

   // make first window in list the current window
   for (wid = 0; wid < fMaxNumberOfWindows; wid++) {
      if (fWindows[wid].open) {
         gCws = &fWindows[wid];
         return;
      }
   }
   gCws = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy the pixmap wid at the position xpos, ypos in the current window.

void TGWin32::CopyPixmap(int wid, int xpos, int ypos)
{
   if (!fWindows) return;

   gTws = &fWindows[wid];
   gdk_window_copy_area(gCws->drawing, gGCpxmp, xpos, ypos, gTws->drawing,
                        0, 0, gTws->width, gTws->height);
   GdiFlush();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a box.
/// mode=0 hollow  (kHollow)
/// mode=1 solid   (kSolid)

void TGWin32::DrawBox(int x1, int y1, int x2, int y2, EBoxMode mode)
{
   if (!fWindows) return;

   Int_t x = TMath::Min(x1, x2);
   Int_t y = TMath::Min(y1, y2);
   Int_t w = TMath::Abs(x2 - x1);
   Int_t h = TMath::Abs(y2 - y1);

   switch (mode) {

   case kHollow:
      if (fLineColorModified) UpdateLineColor();
      if (fPenModified) UpdateLineStyle();
      gdk_win32_draw_rectangle(gCws->drawing, gGCline, 0, x, y, w, h);
      break;

   case kFilled:
      if (fFillStyleModified) UpdateFillStyle();
      if (fFillColorModified) UpdateFillColor();
      gdk_win32_draw_rectangle(gCws->drawing, gGCfill, 1, x, y, w, h);
      break;

   default:
      break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a cell array.
/// x1,y1        : left down corner
/// x2,y2        : right up corner
/// nx,ny        : array size
/// ic           : array
///
/// Draw a cell array. The drawing is done with the pixel presicion
/// if (X2-X1)/NX (or Y) is not a exact pixel number the position of
/// the top rigth corner may be wrong.

void TGWin32::DrawCellArray(Int_t x1, Int_t y1, Int_t x2, Int_t y2,
                            Int_t nx, Int_t ny, Int_t *ic)
{
   int i, j, icol, ix, iy, w, h, current_icol;

   if (!fWindows) return;

   current_icol = -1;
   w = TMath::Max((x2 - x1) / (nx), 1);
   h = TMath::Max((y1 - y2) / (ny), 1);
   ix = x1;

   if (fFillStyleModified) UpdateFillStyle();
   if (fFillColorModified) UpdateFillColor();

   for (i = 0; i < nx; i++) {
      iy = y1 - h;
      for (j = 0; j < ny; j++) {
         icol = ic[i + (nx * j)];
         if (icol != current_icol) {
            gdk_gc_set_foreground(gGCfill, (GdkColor *) & GetColor(icol).color);
            current_icol = icol;
         }

         gdk_win32_draw_rectangle(gCws->drawing, gGCfill, kTRUE, ix,  iy, w, h);
         iy = iy - h;
      }
      ix = ix + w;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill area described by polygon.
/// n         : number of points
/// xy(2,n)   : list of points

void TGWin32::DrawFillArea(int n, TPoint *xyt)
{
   int i;
   static int lastn = 0;
   static GdkPoint *xy = 0;

   if (!fWindows) return;

   if (fFillStyleModified) UpdateFillStyle();
   if (fFillColorModified) UpdateFillColor();

   if (lastn!=n) {
      delete [] (GdkPoint *)xy;
      xy = new GdkPoint[n];
      lastn = n;
   }
   for (i = 0; i < n; i++) {
      xy[i].x = xyt[i].fX;
      xy[i].y = xyt[i].fY;
   }

   if (gFillHollow) {
      gdk_win32_draw_lines(gCws->drawing, gGCfill, xy, n);
   } else {
      gdk_win32_draw_polygon(gCws->drawing, gGCfill, 1, xy, n);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a line.
/// x1,y1        : begin of line
/// x2,y2        : end of line

void TGWin32::DrawLine(int x1, int y1, int x2, int y2)
{
   if (!fWindows) return;

   if (fLineColorModified) UpdateLineColor();
   if (fPenModified) UpdateLineStyle();

   if (gLineStyle == GDK_LINE_SOLID) {
      gdk_draw_line(gCws->drawing, gGCline, x1, y1, x2, y2);
   } else {
      int i;
      gint8 dashes[32];
      for (i = 0; i < gDashSize; i++) {
         dashes[i] = (gint8) gDashList[i];
      }
      for (i = gDashSize; i < 32; i++) {
         dashes[i] = (gint8) 0;
      }
      gdk_gc_set_dashes(gGCdash, gDashOffset, dashes, gDashSize);
      gdk_draw_line(gCws->drawing, gGCdash, x1, y1, x2, y2);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a line through all points.
/// n         : number of points
/// xy        : list of points

void TGWin32::DrawPolyLine(int n, TPoint * xyt)
{
   int i;

   if (!fWindows) return;

   Point_t *xy = new Point_t[n];

   for (i = 0; i < n; i++) {
      xy[i].fX = xyt[i].fX;
      xy[i].fY = xyt[i].fY;
   }

   if (fLineColorModified) UpdateLineColor();
   if (fPenModified) UpdateLineStyle();

   if (n > 1) {
      if (gLineStyle == GDK_LINE_SOLID) {
         gdk_win32_draw_lines(gCws->drawing, gGCline, (GdkPoint *)xy, n);
      } else {
         int i;
         gint8 dashes[32];

         for (i = 0; i < gDashSize; i++) {
            dashes[i] = (gint8) gDashList[i];
         }
         for (i = gDashSize; i < 32; i++) {
            dashes[i] = (gint8) 0;
         }

         gdk_gc_set_dashes(gGCdash, gDashOffset, dashes, gDashSize);
         gdk_win32_draw_lines(gCws->drawing, (GdkGC*)gGCdash, (GdkPoint *)xy, n);

         // calculate length of line to update dash offset
         for (i = 1; i < n; i++) {
            int dx = xy[i].fX - xy[i - 1].fX;
            int dy = xy[i].fY - xy[i - 1].fY;

            if (dx < 0) dx = -dx;
            if (dy < 0) dy = -dy;
            gDashOffset += dx > dy ? dx : dy;
         }
         gDashOffset %= gDashLength;
      }
   } else {
      gdk_win32_draw_points( gCws->drawing, gLineStyle == GDK_LINE_SOLID ?
                              gGCline : gGCdash, (GdkPoint *)xy,1);
   }
   delete [] xy;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw n markers with the current attributes at position x, y.
/// n    : number of markers to draw
/// xy   : x,y coordinates of markers

void TGWin32::DrawPolyMarker(int n, TPoint *xyt)
{
   int i;
   static int lastn = 0;
   static GdkPoint *xy = 0;

   if (!fWindows) return;

   if (fMarkerStyleModified) UpdateMarkerStyle();
   if (fMarkerColorModified) UpdateMarkerColor();

   if (lastn!=n) {
      delete [] (GdkPoint *)xy;
      xy = new GdkPoint[n];
      lastn = n;
   }

   for (i = 0; i < n; i++) {
      xy[i].x = xyt[i].fX;
      xy[i].y = xyt[i].fY;
   }

   if (gMarker.n <= 0) {
       gdk_win32_draw_points(gCws->drawing, gGCmark, xy, n);
   } else {
      int r = gMarker.n / 2;
      int m;

      for (m = 0; m < n; m++) {
         int hollow = 0;
         switch (gMarker.type) {
            int i;

         case 0:               // hollow circle
            gdk_win32_draw_arc(gCws->drawing, gGCmark, kFALSE, xy[m].x-r, xy[m].y-r,
                              gMarker.n, gMarker.n, 0, 23040);
            break;

         case 1:               // filled circle
            gdk_win32_draw_arc(gCws->drawing, gGCmark, kTRUE, xy[m].x-r, xy[m].y-r,
                              gMarker.n, gMarker.n, 0, 23040);
            break;

         case 2:               // hollow polygon
            hollow = 1;
         case 3:               // filled polygon
            for (i = 0; i < gMarker.n; i++) {
               gMarker.xy[i].x += xy[m].x;
               gMarker.xy[i].y += xy[m].y;
            }
            if (hollow) {
               gdk_win32_draw_lines(gCws->drawing, gGCmark, (GdkPoint *)gMarker.xy, gMarker.n);
            } else {
               gdk_win32_draw_polygon(gCws->drawing, gGCmark, 1, (GdkPoint *)gMarker.xy, gMarker.n);
            }
            for (i = 0; i < gMarker.n; i++) {
               gMarker.xy[i].x -= xy[m].x;
               gMarker.xy[i].y -= xy[m].y;
            }
            break;

         case 4:               // segmented line
            for (i = 0; i < gMarker.n; i += 2) {
               gdk_draw_line(gCws->drawing, gGCmark,
                             xy[m].x + gMarker.xy[i].x,
                             xy[m].y + gMarker.xy[i].y,
                             xy[m].x + gMarker.xy[i + 1].x,
                             xy[m].y + gMarker.xy[i + 1].y);
            }
            break;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return character up vector.

void TGWin32::GetCharacterUp(Float_t & chupx, Float_t & chupy)
{
   chupx = fCharacterUpX;
   chupy = fCharacterUpY;
}

////////////////////////////////////////////////////////////////////////////////
/// Return reference to internal color structure associated
/// to color index cid.

XColor_t &TGWin32::GetColor(Int_t cid)
{
   XColor_t *col = (XColor_t*) fColors->GetValue(cid);
   if (!col) {
      col = new XColor_t;
      fColors->Add(cid, (Long_t) col);
   }
   return *col;
}

////////////////////////////////////////////////////////////////////////////////
/// Return current window pointer. Protected method used by TGWin32TTF.

Window_t TGWin32::GetCurrentWindow() const
{
   return (Window_t)(gCws ? gCws->drawing : 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Return desired Graphics Context ("which" maps directly on gGCList[]).
/// Protected method used by TGWin32TTF.

GdkGC *TGWin32::GetGC(Int_t which) const
{
   if (which >= kMAXGC || which < 0) {
      Error("GetGC", "trying to get illegal GdkGC (which = %d)", which);
      return 0;
   }

   return gGClist[which];
}

////////////////////////////////////////////////////////////////////////////////
/// Query the double buffer value for the window wid.

Int_t TGWin32::GetDoubleBuffer(int wid)
{
   if (!fWindows) return 0;

   gTws = &fWindows[wid];

   if (!gTws->open) {
      return -1;
   } else {
      return gTws->double_buffer;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return position and size of window wid.
/// wid        : window identifier
/// x,y        : window position (output)
/// w,h        : window size (output)
/// if wid < 0 the size of the display is returned

void TGWin32::GetGeometry(int wid, int &x, int &y, unsigned int &w,
                          unsigned int &h)
{
   if (!fWindows) return;

   if (wid < 0) {
      x = 0;
      y = 0;

      w = gdk_screen_width();
      h = gdk_screen_height();
   } else {
      int depth;
      int width, height;

      gTws = &fWindows[wid];
      gdk_window_get_geometry((GdkDrawable *) gTws->window, &x, &y,
                              &width, &height, &depth);

      gdk_window_get_deskrelative_origin((GdkDrawable *) gTws->window, &x, &y);

      if (width > 0 && height > 0) {
         gTws->width = width;
         gTws->height = height;
      }
      w = gTws->width;
      h = gTws->height;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return hostname on which the display is opened.

const char *TGWin32::DisplayName(const char *dpyName)
{
   return "localhost";          //return gdk_get_display();
}

////////////////////////////////////////////////////////////////////////////////
/// Get maximum number of planes.

void TGWin32::GetPlanes(int &nplanes)
{
   nplanes = gdk_visual_get_best_depth();
}

////////////////////////////////////////////////////////////////////////////////
/// Get rgb values for color "index".

void TGWin32::GetRGB(int index, float &r, float &g, float &b)
{
   if (index == 0) {
      r = g = b = 1.0;
   } else if (index == 1) {
      r = g = b = 0.0;
   } else {
      XColor_t &col = GetColor(index);
      r = ((float) col.color.red) / ((float) kBIGGEST_RGB_VALUE);
      g = ((float) col.color.green) / ((float) kBIGGEST_RGB_VALUE);
      b = ((float) col.color.blue) / ((float) kBIGGEST_RGB_VALUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the size of a character string.
/// iw          : text width
/// ih          : text height
/// mess        : message

void TGWin32::GetTextExtent(unsigned int &w, unsigned int &h, char *mess)
{
   TTF::SetTextFont(gTextFont);
   TTF::SetTextSize(fTextSize);
   TTF::GetTextExtent(w, h, mess);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the X11 window identifier.
/// wid      : Workstation identifier (input)

Window_t TGWin32::GetWindowID(int wid)
{
   if (!fWindows) return 0;
   return (Window_t) fWindows[wid].window;
}

////////////////////////////////////////////////////////////////////////////////
/// Move the window wid.
/// wid  : GdkWindow identifier.
/// x    : x new window position
/// y    : y new window position

void TGWin32::MoveWindow(int wid, int x, int y)
{
   if (!fWindows) return;

   gTws = &fWindows[wid];
   if (!gTws->open) return;

   gdk_window_move((GdkDrawable *) gTws->window, x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Open a new pixmap.
/// w,h : Width and height of the pixmap.

Int_t TGWin32::OpenPixmap(unsigned int w, unsigned int h)
{
   int wval, hval;
   int i, wid;
   int ww, hh, depth;
   wval = w;
   hval = h;

   // Select next free window number
 again:
   for (wid = 0; wid < fMaxNumberOfWindows; wid++) {
      if (!fWindows[wid].open) {
         fWindows[wid].open = 1;
         gCws = &fWindows[wid];
         break;
      }
   }
   if (wid == fMaxNumberOfWindows) {
      int newsize = fMaxNumberOfWindows + 10;
      fWindows = (XWindow_t *) TStorage::ReAlloc(fWindows,
                                                 newsize * sizeof(XWindow_t),
                                                 fMaxNumberOfWindows *
                                                 sizeof(XWindow_t));

      for (i = fMaxNumberOfWindows; i < newsize; i++) fWindows[i].open = 0;
      fMaxNumberOfWindows = newsize;
      goto again;
   }

   depth =gdk_visual_get_best_depth();
   gCws->window = (GdkPixmap *) gdk_pixmap_new(GDK_ROOT_PARENT(),wval,hval,depth);
   gdk_drawable_get_size((GdkDrawable *) gCws->window, &ww, &hh);

   for (i = 0; i < kMAXGC; i++) {
      gdk_gc_set_clip_mask((GdkGC *) gGClist[i], (GdkDrawable *)None);
   }

   SetColor(gGCpxmp, 0);
   gdk_win32_draw_rectangle(gCws->window,(GdkGC *)gGCpxmp, kTRUE,
                           0, 0, ww, hh);
   SetColor(gGCpxmp, 1);

   // Initialise the window structure
   gCws->drawing = gCws->window;
   gCws->buffer = 0;
   gCws->double_buffer = 0;
   gCws->ispixmap = 1;
   gCws->clip = 0;
   gCws->width = wval;
   gCws->height = hval;
   gCws->new_colors = 0;

   return wid;
}

////////////////////////////////////////////////////////////////////////////////
/// Open window and return window number.
/// Return -1 if window initialization fails.

Int_t TGWin32::InitWindow(ULong_t win)
{
   GdkWindowAttr attributes;
   unsigned long attr_mask = 0;
   int wid;
   int xval, yval;
   int wval, hval, depth;

   GdkWindow *wind = (GdkWindow *) win;

   gdk_window_get_geometry(wind, &xval, &yval, &wval, &hval, &depth);

   // Select next free window number

 again:
   for (wid = 0; wid < fMaxNumberOfWindows; wid++) {
      if (!fWindows[wid].open) {
         fWindows[wid].open = 1;
         fWindows[wid].double_buffer = 0;
         gCws = &fWindows[wid];
         break;
      }
   }

   if (wid == fMaxNumberOfWindows) {
      int newsize = fMaxNumberOfWindows + 10;
      fWindows =
          (XWindow_t *) TStorage::ReAlloc(fWindows,
                                          newsize * sizeof(XWindow_t),
                                          fMaxNumberOfWindows *
                                          sizeof(XWindow_t));

      for (int i = fMaxNumberOfWindows; i < newsize; i++) {
         fWindows[i].open = 0;
      }

      fMaxNumberOfWindows = newsize;
      goto again;
   }
   // Create window
   attributes.wclass = GDK_INPUT_OUTPUT;
   attributes.event_mask = 0L;  //GDK_ALL_EVENTS_MASK;
   attributes.event_mask |= GDK_EXPOSURE_MASK | GDK_STRUCTURE_MASK |
       GDK_PROPERTY_CHANGE_MASK;
//                            GDK_ENTER_NOTIFY_MASK | GDK_LEAVE_NOTIFY_MASK;
   if (xval >= 0) {
      attributes.x = xval;
   } else {
      attributes.x = -1.0 * xval;
   }

   if (yval >= 0) {
      attributes.y = yval;
   } else {
      attributes.y = -1.0 * yval;
   }
   attributes.width = wval;
   attributes.height = hval;
   attributes.colormap = gdk_colormap_get_system();
   attributes.visual = gdk_window_get_visual(wind);
   attributes.override_redirect = TRUE;

   if ((attributes.y > 0) && (attributes.x > 0)) {
      attr_mask = GDK_WA_X | GDK_WA_Y | GDK_WA_COLORMAP |
          GDK_WA_WMCLASS | GDK_WA_NOREDIR;
   } else {
      attr_mask = GDK_WA_COLORMAP | GDK_WA_WMCLASS | GDK_WA_NOREDIR;
   }

   if (attributes.visual != NULL) {
      attr_mask |= GDK_WA_VISUAL;
   }
   attributes.window_type = GDK_WINDOW_CHILD;
   gCws->window = gdk_window_new(wind, &attributes, attr_mask);
   HWND window = (HWND)GDK_DRAWABLE_XID((GdkWindow *)gCws->window);
   ::ShowWindow(window, SW_SHOWNORMAL);
   ::ShowWindow(window, SW_RESTORE);
   ::BringWindowToTop(window);

   if (!fUseSysPointers) {
      ::SetClassLong(window, GCL_HCURSOR,
                    (LONG)GDK_CURSOR_XID(fCursors[kPointer]));
   }

   // Initialise the window structure

   gCws->drawing = gCws->window;
   gCws->buffer = 0;
   gCws->double_buffer = 0;
   gCws->ispixmap = 0;
   gCws->clip = 0;
   gCws->width = wval;
   gCws->height = hval;
   gCws->new_colors = 0;

   return wid;
}

////////////////////////////////////////////////////////////////////////////////
/// Query pointer position.
/// ix       : X coordinate of pointer
/// iy       : Y coordinate of pointer
/// (both coordinates are relative to the origin of the root window)

void TGWin32::QueryPointer(int &ix, int &iy)
{
   //GdkModifierType mask;
   //GdkWindow *retw = gdk_window_get_pointer((GdkWindow *) gCws->window,
   //                                          &ix, &iy, &mask);
   POINT cpt;
   GetCursorPos(&cpt);
   ix = cpt.x;
   iy = cpt.y;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the pixmap pix.

void TGWin32::RemovePixmap(GdkDrawable *pix)
{
   gdk_pixmap_unref((GdkPixmap *)pix);
}

////////////////////////////////////////////////////////////////////////////////
/// Request Locator position.
/// x,y       : cursor position at moment of button press (output)
/// ctyp      : cursor type (input)
///   ctyp=1 tracking cross
///   ctyp=2 cross-hair
///   ctyp=3 rubber circle
///   ctyp=4 rubber band
///   ctyp=5 rubber rectangle
///
/// mode      : input mode
///   mode=0 request
///   mode=1 sample
///
/// Request locator:
/// return button number  1 = left is pressed
///                       2 = middle is pressed
///                       3 = right is pressed
///        in sample mode:
///                      11 = left is released
///                      12 = middle is released
///                      13 = right is released
///                      -1 = nothing is pressed or released
///                      -2 = leave the window
///                    else = keycode (keyboard is pressed)

Int_t TGWin32::RequestLocator(Int_t mode, Int_t ctyp, Int_t & x, Int_t & y)
{
   static int xloc = 0;
   static int yloc = 0;
   static int xlocp = 0;
   static int ylocp = 0;
   static GdkCursor *cursor = NULL;
   Int_t  xtmp, ytmp;

   GdkEvent *event;
   int button_press;
   int radius;

   // Change the cursor shape
   if (cursor == NULL) {
      if (ctyp > 1) {
         gdk_window_set_cursor((GdkWindow *)gCws->window, (GdkCursor *)gNullCursor);
         gdk_gc_set_foreground((GdkGC *) gGCecho, &GetColor(0).color);
      } else {
         if (fUseSysPointers)
            cursor = gdk_syscursor_new((ULong_t)IDC_CROSS);
         else
            cursor = gdk_cursor_new((GdkCursorType)GDK_CROSSHAIR);
         gdk_window_set_cursor((GdkWindow *)gCws->window, (GdkCursor *)cursor);
      }
   }

   // Event loop
   button_press = 0;

   // Set max response time to 2 minutes to avoid timeout
   // in TGWin32ProxyBase::ForwardCallBack during RequestLocator
   TGWin32VirtualXProxy::fMaxResponseTime = 120000;
   while (button_press == 0) {
      event = gdk_event_get();

      switch (ctyp) {

      case 1:
         break;

      case 2:
         gdk_draw_line(gCws->window, gGCecho, xloc, 0, xloc, gCws->height);
         gdk_draw_line(gCws->window, gGCecho, 0, yloc, gCws->width, yloc);
         break;

      case 3:
         radius = (int) TMath::Sqrt((double)((xloc - xlocp) * (xloc - xlocp) +
                                             (yloc - ylocp) * (yloc - ylocp)));

         gdk_win32_draw_arc(gCws->window, gGCecho, kFALSE,
                            xlocp - radius, ylocp - radius,
                            2 * radius, 2 * radius, 0, 23040);
         break;

      case 4:
         gdk_draw_line(gCws->window, gGCecho, xlocp, ylocp, xloc, yloc);
         break;

      case 5:
         gdk_win32_draw_rectangle( gCws->window, gGCecho, kFALSE,
                                 TMath::Min(xlocp, xloc), TMath::Min(ylocp, yloc),
                                 TMath::Abs(xloc - xlocp), TMath::Abs(yloc - ylocp));
         break;

      default:
         break;
      }

      xloc = event->button.x;
      yloc = event->button.y;

      switch (event->type) {

      case GDK_LEAVE_NOTIFY:
         if (mode == 0) {
            while (1) {
               event = gdk_event_get();

               if (event->type == GDK_ENTER_NOTIFY) {
                  gdk_event_free(event);
                  break;
               }
               gdk_event_free(event);
            }
         } else {
            button_press = -2;
         }
         break;

      case GDK_BUTTON_PRESS:
         button_press = event->button.button;
         xlocp = event->button.x;
         ylocp = event->button.y;
         gdk_cursor_unref(cursor);
         cursor = 0;
         break;

      case GDK_BUTTON_RELEASE:
         if (mode == 1) {
            button_press = 10 + event->button.button;
            xlocp = event->button.x;
            ylocp = event->button.y;
         }
         break;

      case GDK_KEY_PRESS:
         if (mode == 1) {
            button_press = event->key.keyval;
            xlocp = event->button.x;
            ylocp = event->button.y;
         }
         break;

      case GDK_KEY_RELEASE:
         if (mode == 1) {
            button_press = -1 * (int)(event->key.keyval);
            xlocp = event->button.x;
            ylocp = event->button.y;
         }
         break;

      default:
         break;
      }

      xtmp = event->button.x;
      ytmp = event->button.y;

      gdk_event_free(event);

      if (mode == 1) {
         if (button_press == 0) {
            button_press = -1;
         }
         break;
      }
   }
   TGWin32VirtualXProxy::fMaxResponseTime = 1000;

   x = xtmp;
   y = ytmp;

   return button_press;
}

////////////////////////////////////////////////////////////////////////////////
/// Request a string.
/// x,y         : position where text is displayed
/// text        : text displayed (input), edited text (output)
///
/// Request string:
/// text is displayed and can be edited with Emacs-like keybinding
/// return termination code (0 for ESC, 1 for RETURN)

Int_t TGWin32::RequestString(int x, int y, char *text)
{
   static GdkCursor *cursor = NULL;
   static int percent = 0;      // bell volume
   static GdkWindow *CurWnd;
   HWND focuswindow;
   GdkEvent *event;
   KeySym keysym;
   int key = -1;
   int len_text = strlen(text);
   int nt;                      // defined length of text
   int pt;                      // cursor position in text

   CurWnd = (GdkWindow *)gCws->window;
   // change the cursor shape
   if (cursor == NULL) {
      if (fUseSysPointers)
         cursor = gdk_syscursor_new((ULong_t)IDC_HELP);
      else
         cursor = gdk_cursor_new((GdkCursorType)GDK_QUESTION_ARROW);
   }
   if (cursor != 0) {
      gdk_window_set_cursor(CurWnd, cursor);
   }
   for (nt = len_text; nt > 0 && text[nt - 1] == ' '; nt--);

   pt = nt;
   focuswindow = ::SetFocus((HWND)GDK_DRAWABLE_XID(CurWnd));

   // Set max response time to 2 minutes to avoid timeout
   // in TGWin32ProxyBase::ForwardCallBack during RequestString
   TGWin32VirtualXProxy::fMaxResponseTime = 120000;
   TTF::SetTextFont(gTextFont);
   TTF::SetTextSize(fTextSize);
   do {
      char tmp[2];
      char keybuf[8];
      char nbytes;
      UInt_t dx, ddx, h;
      int i;

      if (EventsPending()) {
         event = gdk_event_get();
      } else {
         gSystem->ProcessEvents();
         ::SleepEx(10, kTRUE);
         continue;
      }

      DrawText(x, y, 0.0, 1.0, text, kOpaque);
      TTF::GetTextExtent(dx, h, text);
      DrawText(x+dx, y, 0.0, 1.0, " ", kOpaque);

      if (pt == 0) {
         dx = 0;
      } else {
         char *stmp = new char[pt+1];
         strncpy(stmp, text, pt);
         stmp[pt] = '\0';
         TTF::GetTextExtent(ddx, h, stmp);
         dx = ddx;
         delete[] stmp;
      }

      if (pt < len_text) {
         tmp[0] = text[pt];
         tmp[1] = '\0';
         DrawText(x+dx, y, 0.0, 1.0, tmp, kOpaque);
      } else {
         DrawText(x+dx, y, 0.0, 1.0, " ", kOpaque);
      }

      if (event != NULL) {
         switch (event->type) {
         case GDK_BUTTON_PRESS:
         case GDK_ENTER_NOTIFY:
            focuswindow = ::SetFocus((HWND)GDK_DRAWABLE_XID(CurWnd));
            break;

         case GDK_LEAVE_NOTIFY:
            ::SetFocus(focuswindow);
            break;
         case GDK_KEY_PRESS:
            nbytes = event->key.length;
            for (i = 0; i < nbytes; i++) {
               keybuf[i] = event->key.string[i];
            }
            keysym = event->key.keyval;
            switch (keysym) {   // map cursor keys
            case GDK_BackSpace:
               keybuf[0] = 0x08; // backspace
               nbytes = 1;
               break;
            case GDK_Return:
               keybuf[0] = 0x0d; // return
               nbytes = 1;
               break;
            case GDK_Delete:
               keybuf[0] = 0x7f; // del
               nbytes = 1;
               break;
            case GDK_Escape:
               keybuf[0] = 0x1b; // esc
               nbytes = 1;
               break;
            case GDK_Home:
               keybuf[0] = 0x01; // home
               nbytes = 1;
               break;
            case GDK_Left:
               keybuf[0] = 0x02; // backward
               nbytes = 1;
               break;
            case GDK_Right:
               keybuf[0] = 0x06; // forward
               nbytes = 1;
               break;
            case GDK_End:
               keybuf[0] = 0x05; // end
               nbytes = 1;
               break;
            }
            if (nbytes == 1) {
               if (isascii(keybuf[0]) && isprint(keybuf[0])) {
                  // insert character
                  if (nt < len_text) {
                     nt++;
                  }
                  for (i = nt - 1; i > pt; i--) {
                     text[i] = text[i - 1];
                  }
                  if (pt < len_text) {
                     text[pt] = keybuf[0];
                     pt++;
                  }
               } else {
                  switch (keybuf[0]) {
                     // Emacs-like editing keys

                  case 0x08:   //'\010':    // backspace
                  case 0x7f:   //'\177':    // delete
                     // delete backward
                     if (pt > 0) {
                        for (i = pt; i < nt; i++) {
                           text[i - 1] = text[i];
                        }
                        text[nt - 1] = ' ';
                        nt--;
                        pt--;
                     }
                     break;
                  case 0x01:   //'\001':    // ^A
                     // beginning of line
                     pt = 0;
                     break;
                  case 0x02:   //'\002':    // ^B
                     // move backward
                     if (pt > 0) {
                        pt--;
                     }
                     break;
                  case 0x04:   //'\004':    // ^D
                     // delete forward
                     if (pt > 0) {
                        for (i = pt; i < nt; i++) {
                           text[i - 1] = text[i];
                        }
                        text[nt - 1] = ' ';
                        pt--;
                     }
                     break;
                  case 0x05:   //'\005':    // ^E
                     // end of line
                     pt = nt;
                     break;

                  case 0x06:   //'\006':    // ^F
                     // move forward
                     if (pt < nt) {
                        pt++;
                     }
                     break;
                  case 0x0b:   //'\013':    // ^K
                     // delete to end of line
                     for (i = pt; i < nt; i++)
                        text[i] = ' ';
                     nt = pt;
                     break;
                  case 0x14:   //'\024':    // ^T
                     // transpose
                     if (pt > 0) {
                        char c = text[pt];
                        text[pt] = text[pt - 1];
                        text[pt - 1] = c;
                     }
                     break;
                  case 0x0A:   //'\012':    // newline
                  case 0x0D:   //'\015':    // return
                     key = 1;
                     break;
                  case 0x1B:   //'\033':    // escape
                     key = 0;
                     break;

                  default:
                     gSystem->Beep();
                     break;
                  }
               }
            }
         default:
            SetInputFocus((Window_t)gCws->window);
            break;
         }
         gdk_event_free(event);
      }
   } while (key < 0);
   TGWin32VirtualXProxy::fMaxResponseTime = 1000;
   ::SetFocus(focuswindow);
   SetInputFocus((Window_t)CurWnd);

   gdk_window_set_cursor(CurWnd, (GdkCursor *)fCursors[kPointer]);
   if (cursor != 0) {
      gdk_cursor_unref(cursor);
      cursor = 0;
   }

   return key;
}

////////////////////////////////////////////////////////////////////////////////
/// Rescale the window wid.
/// wid  : GdkWindow identifier
/// w    : Width
/// h    : Heigth

void TGWin32::RescaleWindow(int wid, unsigned int w, unsigned int h)
{
    int i;

   if (!fWindows) return;

   gTws = &fWindows[wid];
   if (!gTws->open)
      return;

   // don't do anything when size did not change
   if (gTws->width == w && gTws->height == h)
      return;

   gdk_window_resize((GdkWindow *) gTws->window, w, h);

   if (gTws->buffer) {
      // don't free and recreate pixmap when new pixmap is smaller
      if (gTws->width < w || gTws->height < h) {
         gdk_pixmap_unref(gTws->buffer);
         gTws->buffer = gdk_pixmap_new(GDK_ROOT_PARENT(), // NULL,
                                       w, h, gdk_visual_get_best_depth());
      }
      for (i = 0; i < kMAXGC; i++) {
         gdk_gc_set_clip_mask(gGClist[i], None);
      }
      SetColor(gGCpxmp, 0);
      gdk_win32_draw_rectangle(gTws->buffer, gGCpxmp, 1, 0, 0, w, h);
      SetColor(gGCpxmp, 1);

      if (gTws->double_buffer) gTws->drawing = gTws->buffer;
   }
   gTws->width = w;
   gTws->height = h;
}

////////////////////////////////////////////////////////////////////////////////
/// Resize a pixmap.
/// wid : pixmap to be resized
/// w,h : Width and height of the pixmap

int TGWin32::ResizePixmap(int wid, unsigned int w, unsigned int h)
{
   int wval, hval;
   int i;
   int ww, hh, depth;
   wval = w;
   hval = h;

   if (!fWindows) return 0;

   gTws = &fWindows[wid];

   // don't do anything when size did not change
   //  if (gTws->width == wval && gTws->height == hval) return 0;

   // due to round-off errors in TPad::Resize() we might get +/- 1 pixel
   // change, in those cases don't resize pixmap
   if (gTws->width >= wval - 1 && gTws->width <= wval + 1 &&
       gTws->height >= hval - 1 && gTws->height <= hval + 1)
      return 0;

   // don't free and recreate pixmap when new pixmap is smaller
   if (gTws->width < wval || gTws->height < hval) {
      gdk_pixmap_unref((GdkPixmap *)gTws->window);
      depth = gdk_visual_get_best_depth();
      gTws->window = gdk_pixmap_new(GDK_ROOT_PARENT(), wval, hval, depth);
   }

   gdk_drawable_get_size(gTws->window, &ww, &hh);

   for (i = 0; i < kMAXGC; i++) {
      gdk_gc_set_clip_mask((GdkGC *) gGClist[i], (GdkDrawable *)None);
   }

   SetColor(gGCpxmp, 0);
   gdk_win32_draw_rectangle(gTws->window,(GdkGC *)gGCpxmp, kTRUE, 0, 0, ww, hh);
   SetColor(gGCpxmp, 1);

   // Initialise the window structure
   gTws->drawing = gTws->window;
   gTws->width = wval;
   gTws->height = hval;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Resize the current window if necessary.

void TGWin32::ResizeWindow(int wid)
{
   int i;
   int xval = 0, yval = 0;
   GdkWindow *win, *root = NULL;
   int wval = 0, hval = 0, depth = 0;

   if (!fWindows) return;

   gTws = &fWindows[wid];

   win = (GdkWindow *) gTws->window;
   gdk_window_get_geometry(win, &xval, &yval,
                           &wval, &hval, &depth);

   // don't do anything when size did not change
   if (gTws->width == wval && gTws->height == hval) {
      return;
   }

   gdk_window_resize((GdkWindow *) gTws->window, wval, hval);

   if (gTws->buffer) {
      if (gTws->width < wval || gTws->height < hval) {
         gdk_pixmap_unref((GdkPixmap *)gTws->buffer);
         depth = gdk_visual_get_best_depth();
         gTws->buffer = (GdkPixmap *) gdk_pixmap_new(GDK_ROOT_PARENT(),
                                                     wval, hval, depth);
      }

      for (i = 0; i < kMAXGC; i++) {
         gdk_gc_set_clip_mask((GdkGC *) gGClist[i], (GdkDrawable *)None);
      }

      SetColor(gGCpxmp, 0);
      gdk_win32_draw_rectangle(gTws->buffer,(GdkGC *)gGCpxmp, kTRUE, 0, 0, wval, hval);

      SetColor(gGCpxmp, 1);

      if (gTws->double_buffer) gTws->drawing = gTws->buffer;
   }

   gTws->width = wval;
   gTws->height = hval;
}

////////////////////////////////////////////////////////////////////////////////
/// Select window to which subsequent output is directed.

void TGWin32::SelectWindow(int wid)
{
   int i;
   GdkRectangle rect;

   if (!fWindows || wid < 0 || wid >= fMaxNumberOfWindows || !fWindows[wid].open) {
      return;
   }

   gCws = &fWindows[wid];

   if (gCws->clip && !gCws->ispixmap && !gCws->double_buffer) {
      rect.x = gCws->xclip;
      rect.y = gCws->yclip;
      rect.width = gCws->wclip;
      rect.height = gCws->hclip;

      for (i = 0; i < kMAXGC; i++) {
         gdk_gc_set_clip_rectangle((GdkGC *) gGClist[i], &rect);
      }
   } else {
      for (i = 0; i < kMAXGC; i++) {
         gdk_gc_set_clip_mask((GdkGC *) gGClist[i], (GdkDrawable *)None);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set character up vector.

void TGWin32::SetCharacterUp(Float_t chupx, Float_t chupy)
{
   if (chupx == fCharacterUpX && chupy == fCharacterUpY) return;

   if (chupx == 0 && chupy == 0) {
      fTextAngle = 0;
   } else if (chupx == 0 && chupy == 1) {
      fTextAngle = 0;
   } else if (chupx == -1 && chupy == 0) {
      fTextAngle = 90;
   } else if (chupx == 0 && chupy == -1) {
      fTextAngle = 180;
   } else if (chupx == 1 && chupy == 0) {
      fTextAngle = 270;
   } else {
      fTextAngle =
          ((TMath::
            ACos(chupx / TMath::Sqrt(chupx * chupx + chupy * chupy)) *
            180.) / 3.14159) - 90;
      if (chupy < 0) fTextAngle = 180 - fTextAngle;
      if (TMath::Abs(fTextAngle) <= 0.01) fTextAngle = 0;
   }
   fCharacterUpX = chupx;
   fCharacterUpY = chupy;
}

////////////////////////////////////////////////////////////////////////////////
/// Turn off the clipping for the window wid.

void TGWin32::SetClipOFF(int wid)
{
   if (!fWindows) return;

   gTws = &fWindows[wid];
   gTws->clip = 0;

   for (int i = 0; i < kMAXGC; i++) {
      gdk_gc_set_clip_mask((GdkGC *) gGClist[i], (GdkDrawable *)None);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set clipping region for the window wid.
/// wid        : GdkWindow indentifier
/// x,y        : origin of clipping rectangle
/// w,h        : size of clipping rectangle;

void TGWin32::SetClipRegion(int wid, int x, int y, unsigned int w,
                            unsigned int h)
{
   if (!fWindows) return;

   gTws = &fWindows[wid];
   gTws->xclip = x;
   gTws->yclip = y;
   gTws->wclip = w;
   gTws->hclip = h;
   gTws->clip = 1;
   GdkRectangle rect;

   if (gTws->clip && !gTws->ispixmap && !gTws->double_buffer) {
      rect.x = gTws->xclip;
      rect.y = gTws->yclip;
      rect.width = gTws->wclip;
      rect.height = gTws->hclip;

      for (int i = 0; i < kMAXGC; i++) {
         gdk_gc_set_clip_rectangle((GdkGC *)gGClist[i], &rect);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return pixel value associated to specified ROOT color number.

ULong_t TGWin32::GetPixel(Color_t ci)
{
   TColor *color = gROOT->GetColor(ci);
   if (color)
      SetRGB(ci, color->GetRed(), color->GetGreen(), color->GetBlue());
   XColor_t &col = GetColor(ci);
   return col.color.pixel;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the foreground color in GdkGC.

void TGWin32::SetColor(GdkGC *gc, int ci)
{
   GdkGCValues gcvals;
   GdkColor color;

   if (ci<=0) ci = 10; //white

   TColor *clr = gROOT->GetColor(ci);
   if (clr)
      SetRGB(ci, clr->GetRed(), clr->GetGreen(), clr->GetBlue());

   XColor_t &col = GetColor(ci);
   if (fColormap && !col.fDefined) {
      col = GetColor(0);
   } else if (!fColormap && (ci < 0 || ci > 1)) {
      col = GetColor(0);
   }

   if (fDrawMode == kXor) {
      gdk_gc_get_values(gc, &gcvals);

      color.pixel = col.color.pixel ^ gcvals.background.pixel;
      color.red = GetRValue(color.pixel);
      color.green = GetGValue(color.pixel);
      color.blue = GetBValue(color.pixel);
      gdk_gc_set_foreground(gc, &color);

   } else {
      gdk_gc_set_foreground(gc, &col.color);

      // make sure that foreground and background are different
      gdk_gc_get_values(gc, &gcvals);

      if (gcvals.foreground.pixel != gcvals.background.pixel) {
         gdk_gc_set_background(gc, &GetColor(!ci).color);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the cursor.

void TGWin32::SetCursor(int wid, ECursor cursor)
{
   if (!fWindows) return;

   gTws = &fWindows[wid];
   gdk_window_set_cursor((GdkWindow *)gTws->window, (GdkCursor *)fCursors[cursor]);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the specified cursor.

void TGWin32::SetCursor(Window_t id, Cursor_t curid)
{
   if (!id) return;

   static GdkWindow *lid = 0;
   static GdkCursor *lcur = 0;

   if ((lid == (GdkWindow *)id) && (lcur==(GdkCursor *)curid)) return;
   lid = (GdkWindow *)id;
   lcur = (GdkCursor *)curid;

   gdk_window_set_cursor((GdkWindow *) id, (GdkCursor *)curid);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the double buffer on/off on window wid.
/// wid  : GdkWindow identifier.
///        999 means all the opened windows.
/// mode : 1 double buffer is on
///        0 double buffer is off

void TGWin32::SetDoubleBuffer(int wid, int mode)
{
   if (!fWindows) return;

   if (wid == 999) {
      for (int i = 0; i < fMaxNumberOfWindows; i++) {
         gTws = &fWindows[i];
         if (gTws->open) {
            switch (mode) {
            case 1:
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
      if (!gTws->open) return;

      switch (mode) {
      case 1:
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

void TGWin32::SetDoubleBufferOFF()
{
   if (!gTws->double_buffer) return;
   gTws->double_buffer = 0;
   gTws->drawing = gTws->window;
}

////////////////////////////////////////////////////////////////////////////////
/// Turn double buffer mode on.

void TGWin32::SetDoubleBufferON()
{
   if (!fWindows || gTws->double_buffer || gTws->ispixmap) return;

   if (!gTws->buffer) {
      gTws->buffer = gdk_pixmap_new(GDK_ROOT_PARENT(), //NULL,
                                    gTws->width, gTws->height,
                                    gdk_visual_get_best_depth());
      SetColor(gGCpxmp, 0);
      gdk_win32_draw_rectangle(gTws->buffer, gGCpxmp, 1, 0, 0, gTws->width,
                         gTws->height);
      SetColor(gGCpxmp, 1);
   }
   for (int i = 0; i < kMAXGC; i++) {
      gdk_gc_set_clip_mask(gGClist[i], None);
   }
   gTws->double_buffer = 1;
   gTws->drawing = gTws->buffer;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the drawing mode.
/// mode : drawing mode
///   mode=1 copy
///   mode=2 xor
///   mode=3 invert
///   mode=4 set the suitable mode for cursor echo according to
///          the vendor

void TGWin32::SetDrawMode(EDrawMode mode)
{
   int i;

   switch (mode) {
   case kCopy:
      for (i = 0; i < kMAXGC; i++) {
         gdk_gc_set_function(gGClist[i], GDK_COPY);
      }
      break;

   case kXor:
      for (i = 0; i < kMAXGC; i++) {
         gdk_gc_set_function(gGClist[i], GDK_XOR);
      }
      break;

   case kInvert:
      for (i = 0; i < kMAXGC; i++) {
         gdk_gc_set_function(gGClist[i], GDK_INVERT);
      }
      break;
   }
   fDrawMode = mode;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for fill areas.

void TGWin32::SetFillColor(Color_t cindex)
{
   Int_t indx = Int_t(cindex);

   if (!gStyle->GetFillColor() && cindex > 1) {
      indx = 0;
   }

   fFillColor = indx;
   fFillColorModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::UpdateFillColor()
{
   if (fFillColor >= 0) {
      SetColor(gGCfill, fFillColor);
   }

   // invalidate fill pattern
   if (gFillPattern != NULL) {
      gdk_pixmap_unref(gFillPattern);
      gFillPattern = NULL;
   }
   fFillColorModified = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set fill area style.
/// fstyle   : compound fill area interior style
///    fstyle = 1000*interiorstyle + styleindex

void TGWin32::SetFillStyle(Style_t fstyle)
{
   if (fFillStyle==fstyle) return;

   fFillStyle = fstyle;
   fFillStyleModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set fill area style index.

void TGWin32::UpdateFillStyle()
{
   static int current_fasi = 0;

   Int_t style = fFillStyle / 1000;
   Int_t fasi = fFillStyle % 1000;

   switch (style) {

   case 1:                     // solid
      gFillHollow = 0;
      gdk_gc_set_fill(gGCfill, GDK_SOLID);
      break;

   case 2:                     // pattern
      gFillHollow = 1;
      break;

   case 3:                     // hatch
      gFillHollow = 0;
      gdk_gc_set_fill(gGCfill, GDK_STIPPLED);

      if (fasi != current_fasi) {
         if (gFillPattern != NULL) {
            gdk_pixmap_unref(gFillPattern);
            gFillPattern = NULL;
         }
         int stn = (fasi >= 1 && fasi <=25) ? fasi : 2;
         char pattern[32];
         for (int i=0;i<32;++i)
            pattern[i] = ~gStipples[stn][i];
         gFillPattern = gdk_bitmap_create_from_data(GDK_ROOT_PARENT(),
                                                    (const char *)&pattern, 16, 16);
         gdk_gc_set_stipple(gGCfill, gFillPattern);
         current_fasi = fasi;
      }
      break;

   default:
      gFillHollow = 1;
   }

   fFillStyleModified = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set input on or off.

void TGWin32::SetInput(int inp)
{
   EnableWindow((HWND) GDK_DRAWABLE_XID(gCws->window), inp);
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for lines.

void TGWin32::SetLineColor(Color_t cindex)
{
   if ((cindex < 0) || (cindex==fLineColor)) return;

   fLineColor =  cindex;
   fLineColorModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::UpdateLineColor()
{
   SetColor(gGCline, Int_t(fLineColor));
   SetColor(gGCdash, Int_t(fLineColor));
   fLineColorModified = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set line type.
/// n         : length of dash list
/// dash(n)   : dash segment lengths
///
/// if n <= 0 use solid lines
/// if n >  0 use dashed lines described by DASH(N)
///    e.g. N=4,DASH=(6,3,1,3) gives a dashed-dotted line with dash length 6
///    and a gap of 7 between dashes

void TGWin32::SetLineType(int n, int *dash)
{
   if (n <= 0) {
      gLineStyle = GDK_LINE_SOLID;
      gdk_gc_set_line_attributes(gGCline, gLineWidth,
                                 (GdkLineStyle)gLineStyle,
                                 (GdkCapStyle) gCapStyle,
                                 (GdkJoinStyle) gJoinStyle);
   } else {
      int i;
      gDashSize = TMath::Min((int)sizeof(gDashList),n);
      gDashLength = 0;
      for (i = 0; i < gDashSize; i++) {
         gDashList[i] = dash[i];
         gDashLength += gDashList[i];
      }
      gDashOffset = 0;
      gLineStyle = GDK_LINE_ON_OFF_DASH;
      if (gLineWidth == 0) gLineWidth =1;
      gdk_gc_set_line_attributes(gGCdash, gLineWidth,
                                 (GdkLineStyle) gLineStyle,
                                 (GdkCapStyle) gCapStyle,
                                 (GdkJoinStyle) gJoinStyle);
   }
   fPenModified = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set line style.

void TGWin32::SetLineStyle(Style_t lstyle)
{
   if (fLineStyle == lstyle) return;

   fLineStyle = lstyle;
   fPenModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Update line style

void TGWin32::UpdateLineStyle()
{
   static Int_t dashed[2] = { 3, 3 };
   static Int_t dotted[2] = { 1, 2 };
   static Int_t dasheddotted[4] = { 3, 4, 1, 4 };

   if (fLineStyle <= 1) {
      SetLineType(0, 0);
   } else if (fLineStyle == 2) {
      SetLineType(2, dashed);
   } else if (fLineStyle == 3) {
      SetLineType(2, dotted);
   } else if (fLineStyle == 4) {
      SetLineType(4, dasheddotted);
   } else {
      TString st = (TString)gStyle->GetLineStyleString(fLineStyle);
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
   fPenModified = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set line width.
/// width   : line width in pixels

void TGWin32::SetLineWidth(Width_t width)
{
   if (fLineWidth == width) return;
   fLineWidth = width;

   if (width == 1 && gLineStyle == GDK_LINE_SOLID) gLineWidth = 0;
   else gLineWidth = width;

   fPenModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for markers.

void TGWin32::SetMarkerColor(Color_t cindex)
{
   if ((cindex<0) || (cindex==fMarkerColor)) return;
   fMarkerColor = cindex;
   fMarkerColorModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::UpdateMarkerColor()
{
   SetColor(gGCmark, Int_t(fMarkerColor));
   fMarkerColorModified = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker size index.
/// msize  : marker scale factor

void TGWin32::SetMarkerSize(Float_t msize)
{
   if ((msize==fMarkerSize) || (msize<0)) return;

   fMarkerSize = msize;
   SetMarkerStyle(-fMarkerStyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker type.
/// type      : marker type
/// n         : length of marker description
/// xy        : list of points describing marker shape
///
/// if n == 0 marker is a single point
/// if TYPE == 0 marker is hollow circle of diameter N
/// if TYPE == 1 marker is filled circle of diameter N
/// if TYPE == 2 marker is a hollow polygon describe by line XY
/// if TYPE == 3 marker is a filled polygon describe by line XY
/// if TYPE == 4 marker is described by segmented line XY
///   e.g. TYPE=4,N=4,XY=(-3,0,3,0,0,-3,0,3) sets a plus shape of 7x7 pixels

void TGWin32::SetMarkerType(int type, int n, GdkPoint * xy)
{
   gMarker.type = type;
   gMarker.n = n < kMAXMK ? n : kMAXMK;
   if (gMarker.type >= 2) {
      for (int i = 0; i < gMarker.n; i++) {
         gMarker.xy[i] = xy[i];
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker style.

void TGWin32::SetMarkerStyle(Style_t markerstyle)
{
   if (fMarkerStyle == markerstyle) return;
   fMarkerStyle = TMath::Abs(markerstyle);
   fMarkerStyleModified = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::UpdateMarkerStyle()
{
   Style_t markerstyle = TAttMarker::GetMarkerStyleBase(fMarkerStyle);
   gMarkerLineWidth = TAttMarker::GetMarkerLineWidth(fMarkerStyle);

   // The fast pixel markers need to be treated separately
   if (markerstyle == 1 || markerstyle == 6 || markerstyle == 7) {
       gdk_gc_set_line_attributes(gGCmark, 0, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
   } else {
       gdk_gc_set_line_attributes(gGCmark, gMarkerLineWidth,
                                  (GdkLineStyle) gMarkerLineStyle,
                                  (GdkCapStyle)  gMarkerCapStyle,
                                  (GdkJoinStyle) gMarkerJoinStyle);
   }

   static GdkPoint shape[30];

   Float_t MarkerSizeReduced = fMarkerSize - TMath::Floor(gMarkerLineWidth/2.)/4.;
   Int_t im = Int_t(4 * MarkerSizeReduced + 0.5);

   if (markerstyle == 2) {
      // + shaped marker
      shape[0].x = -im;
      shape[0].y = 0;
      shape[1].x = im;
      shape[1].y = 0;
      shape[2].x = 0;
      shape[2].y = -im;
      shape[3].x = 0;
      shape[3].y = im;
      SetMarkerType(4, 4, shape);
   } else if (markerstyle == 3 || markerstyle == 31) {
      // * shaped marker
      shape[0].x = -im;
      shape[0].y = 0;
      shape[1].x = im;
      shape[1].y = 0;
      shape[2].x = 0;
      shape[2].y = -im;
      shape[3].x = 0;
      shape[3].y = im;
      im = Int_t(0.707 * Float_t(im) + 0.5);
      shape[4].x = -im;
      shape[4].y = -im;
      shape[5].x = im;
      shape[5].y = im;
      shape[6].x = -im;
      shape[6].y = im;
      shape[7].x = im;
      shape[7].y = -im;
      SetMarkerType(4, 8, shape);
   } else if (markerstyle == 4 || markerstyle == 24) {
      // O shaped marker
      SetMarkerType(0, im * 2, shape);
   } else if (markerstyle == 5) {
      // X shaped marker
      im = Int_t(0.707 * Float_t(im) + 0.5);
      shape[0].x = -im;
      shape[0].y = -im;
      shape[1].x = im;
      shape[1].y = im;
      shape[2].x = -im;
      shape[2].y = im;
      shape[3].x = im;
      shape[3].y = -im;
      SetMarkerType(4, 4, shape);
   } else if (markerstyle == 6) {
      // + shaped marker (with 1 pixel)
      shape[0].x = -1;
      shape[0].y = 0;
      shape[1].x = 1;
      shape[1].y = 0;
      shape[2].x = 0;
      shape[2].y = -1;
      shape[3].x = 0;
      shape[3].y = 1;
      SetMarkerType(4, 4, shape);
   } else if (markerstyle == 7) {
      // . shaped marker (with 9 pixel)
      shape[0].x = -1;
      shape[0].y = 1;
      shape[1].x = 1;
      shape[1].y = 1;
      shape[2].x = -1;
      shape[2].y = 0;
      shape[3].x = 1;
      shape[3].y = 0;
      shape[4].x = -1;
      shape[4].y = -1;
      shape[5].x = 1;
      shape[5].y = -1;
      SetMarkerType(4, 6, shape);
   } else if (markerstyle == 8 || markerstyle == 20) {
      // O shaped marker (filled)
      SetMarkerType(1, im * 2, shape);
   } else if (markerstyle == 21) {
      // full square
      shape[0].x = -im;
      shape[0].y = -im;
      shape[1].x = im;
      shape[1].y = -im;
      shape[2].x = im;
      shape[2].y = im;
      shape[3].x = -im;
      shape[3].y = im;
      shape[4].x = -im;
      shape[4].y = -im;
      SetMarkerType(3, 5, shape);
   } else if (markerstyle == 22) {
      // full triangle up
      shape[0].x = -im;
      shape[0].y = im;
      shape[1].x = im;
      shape[1].y = im;
      shape[2].x = 0;
      shape[2].y = -im;
      shape[3].x = -im;
      shape[3].y = im;
      SetMarkerType(3, 4, shape);
   } else if (markerstyle == 23) {
      // full triangle down
      shape[0].x = 0;
      shape[0].y = im;
      shape[1].x = im;
      shape[1].y = -im;
      shape[2].x = -im;
      shape[2].y = -im;
      shape[3].x = 0;
      shape[3].y = im;
      SetMarkerType(3, 4, shape);
   } else if (markerstyle == 25) {
      // open square
      shape[0].x = -im;
      shape[0].y = -im;
      shape[1].x = im;
      shape[1].y = -im;
      shape[2].x = im;
      shape[2].y = im;
      shape[3].x = -im;
      shape[3].y = im;
      shape[4].x = -im;
      shape[4].y = -im;
      SetMarkerType(2, 5, shape);
   } else if (markerstyle == 26) {
      // open triangle up
      shape[0].x = -im;
      shape[0].y = im;
      shape[1].x = im;
      shape[1].y = im;
      shape[2].x = 0;
      shape[2].y = -im;
      shape[3].x = -im;
      shape[3].y = im;
      SetMarkerType(2, 4, shape);
   } else if (markerstyle == 27) {
      // open losange
      Int_t imx = Int_t(2.66 * MarkerSizeReduced + 0.5);
      shape[0].x = -imx;
      shape[0].y = 0;
      shape[1].x = 0;
      shape[1].y = -im;
      shape[2].x = imx;
      shape[2].y = 0;
      shape[3].x = 0;
      shape[3].y = im;
      shape[4].x = -imx;
      shape[4].y = 0;
      SetMarkerType(2, 5, shape);
   } else if (markerstyle == 28) {
      // open cross
      Int_t imx = Int_t(1.33 * MarkerSizeReduced + 0.5);
      shape[0].x = -im;
      shape[0].y = -imx;
      shape[1].x = -imx;
      shape[1].y = -imx;
      shape[2].x = -imx;
      shape[2].y = -im;
      shape[3].x = imx;
      shape[3].y = -im;
      shape[4].x = imx;
      shape[4].y = -imx;
      shape[5].x = im;
      shape[5].y = -imx;
      shape[6].x = im;
      shape[6].y = imx;
      shape[7].x = imx;
      shape[7].y = imx;
      shape[8].x = imx;
      shape[8].y = im;
      shape[9].x = -imx;
      shape[9].y = im;
      shape[10].x = -imx;
      shape[10].y = imx;
      shape[11].x = -im;
      shape[11].y = imx;
      shape[12].x = -im;
      shape[12].y = -imx;
      SetMarkerType(2, 13, shape);
   } else if (markerstyle == 29) {
      // full star pentagone
      Int_t im1 = Int_t(0.66 * MarkerSizeReduced + 0.5);
      Int_t im2 = Int_t(2.00 * MarkerSizeReduced + 0.5);
      Int_t im3 = Int_t(2.66 * MarkerSizeReduced + 0.5);
      Int_t im4 = Int_t(1.33 * MarkerSizeReduced + 0.5);
      shape[0].x = -im;
      shape[0].y = im4;
      shape[1].x = -im2;
      shape[1].y = -im1;
      shape[2].x = -im3;
      shape[2].y = -im;
      shape[3].x = 0;
      shape[3].y = -im2;
      shape[4].x = im3;
      shape[4].y = -im;
      shape[5].x = im2;
      shape[5].y = -im1;
      shape[6].x = im;
      shape[6].y = im4;
      shape[7].x = im4;
      shape[7].y = im4;
      shape[8].x = 0;
      shape[8].y = im;
      shape[9].x = -im4;
      shape[9].y = im4;
      shape[10].x = -im;
      shape[10].y = im4;
      SetMarkerType(3, 11, shape);
   } else if (markerstyle == 30) {
      // open star pentagone
      Int_t im1 = Int_t(0.66 * MarkerSizeReduced + 0.5);
      Int_t im2 = Int_t(2.00 * MarkerSizeReduced + 0.5);
      Int_t im3 = Int_t(2.66 * MarkerSizeReduced + 0.5);
      Int_t im4 = Int_t(1.33 * MarkerSizeReduced + 0.5);
      shape[0].x = -im;
      shape[0].y = im4;
      shape[1].x = -im2;
      shape[1].y = -im1;
      shape[2].x = -im3;
      shape[2].y = -im;
      shape[3].x = 0;
      shape[3].y = -im2;
      shape[4].x = im3;
      shape[4].y = -im;
      shape[5].x = im2;
      shape[5].y = -im1;
      shape[6].x = im;
      shape[6].y = im4;
      shape[7].x = im4;
      shape[7].y = im4;
      shape[8].x = 0;
      shape[8].y = im;
      shape[9].x = -im4;
      shape[9].y = im4;
      shape[10].x = -im;
      shape[10].y = im4;
      SetMarkerType(2, 11, shape);
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
   } else if (markerstyle == 36) {
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
   } else if (markerstyle == 37) {
      // open three triangles
      Int_t im2 = Int_t(2.0*MarkerSizeReduced + 0.5);
      shape[0].x =   0;  shape[0].y =   0;
      shape[1].x =-im2;  shape[1].y =  im;
      shape[2].x = -im;  shape[2].y =   0;
      shape[3].x =   0;  shape[3].y =   0;
      shape[4].x =-im2;  shape[4].y = -im;
      shape[5].x = im2;  shape[5].y = -im;
      shape[6].x =   0;  shape[6].y =   0;
      shape[7].x =  im;  shape[7].y =   0;
      shape[8].x = im2;  shape[8].y =  im;
      shape[9].x =   0;  shape[9].y =   0;
      SetMarkerType(2,10,shape);
   } else if (markerstyle == 38) {
      // + shaped marker with octagon
      Int_t im2 = Int_t(2.0*MarkerSizeReduced + 0.5);
      shape[0].x = -im;  shape[0].y = 0;
      shape[1].x = -im;  shape[1].y =-im2;
      shape[2].x =-im2;  shape[2].y =-im;
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
      shape[2].x = -im;  shape[2].y =   0;
      shape[3].x =   0;  shape[3].y =   0;
      shape[4].x =-im2;  shape[4].y = -im;
      shape[5].x = im2;  shape[5].y = -im;
      shape[6].x =   0;  shape[6].y =   0;
      shape[7].x =  im;  shape[7].y =   0;
      shape[8].x = im2;  shape[8].y =  im;
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
      SetMarkerType(0, 0, shape);
   }
   fMarkerStyleModified = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set opacity of a window. This image manipulation routine works
/// by adding to a percent amount of neutral to each pixels RGB.
/// Since it requires quite some additional color map entries is it
/// only supported on displays with more than > 8 color planes (> 256
/// colors)

void TGWin32::SetOpacity(Int_t percent)
{
   Int_t depth = gdk_visual_get_best_depth();

   if (depth <= 8) return;
   if (percent == 0) return;

   // if 100 percent then just make white
   ULong_t *orgcolors = 0, *tmpc = 0;
   Int_t maxcolors = 0, ncolors, ntmpc = 0;

   // save previous allocated colors, delete at end when not used anymore
   if (gCws->new_colors) {
      tmpc = gCws->new_colors;
      ntmpc = gCws->ncolors;
   }
   // get pixmap from server as image
   GdkImage *image = gdk_image_get((GdkDrawable*)gCws->drawing, 0, 0,
                                   gCws->width, gCws->height);

   // collect different image colors
   int x, y;
   for (y = 0; y < (int) gCws->height; y++) {
      for (x = 0; x < (int) gCws->width; x++) {
         ULong_t pixel = GetPixelImage((Drawable_t)image, x, y);
         CollectImageColors(pixel, orgcolors, ncolors, maxcolors);
      }
   }
   if (ncolors == 0) {
      gdk_image_unref(image);
      ::operator delete(orgcolors);
      return;
   }
   // create opaque counter parts
   MakeOpaqueColors(percent, orgcolors, ncolors);

   // put opaque colors in image
   for (y = 0; y < (int) gCws->height; y++) {
      for (x = 0; x < (int) gCws->width; x++) {
         ULong_t pixel = GetPixelImage((Drawable_t)image, x, y);
         Int_t idx = FindColor(pixel, orgcolors, ncolors);
         PutPixel((Drawable_t)image, x, y, gCws->new_colors[idx]);
      }
   }

   // put image back in pixmap on server
   gdk_draw_image(gCws->drawing, gGCpxmp, (GdkImage *)image,
                  0, 0, 0, 0, gCws->width, gCws->height);
   GdiFlush();

   // clean up
   if (tmpc) {
      gdk_colors_free((GdkColormap *)fColormap, tmpc, ntmpc, 0);
      delete[]tmpc;
   }
   gdk_image_unref(image);
   ::operator delete(orgcolors);
}

////////////////////////////////////////////////////////////////////////////////
/// Get RGB values for orgcolors, add percent neutral to the RGB and
/// allocate new_colors.

void TGWin32::MakeOpaqueColors(Int_t percent, ULong_t *orgcolors, Int_t ncolors)
{
   Int_t ret;
   if (ncolors <= 0) return;
   GdkColor *xcol = new GdkColor[ncolors];

   int i;
   for (i = 0; i < ncolors; i++) {
      xcol[i].pixel = orgcolors[i];
      xcol[i].red = xcol[i].green = xcol[i].blue = 0;
   }

   GdkColorContext *cc;
   cc = gdk_color_context_new(gdk_visual_get_system(), (GdkColormap *)fColormap);
   gdk_color_context_query_colors(cc, xcol, ncolors);
   gdk_color_context_free(cc);

   UShort_t add = percent * kBIGGEST_RGB_VALUE / 100;

   Int_t val;
   for (i = 0; i < ncolors; i++) {
      val = xcol[i].red + add;
      if (val > kBIGGEST_RGB_VALUE) {
         val = kBIGGEST_RGB_VALUE;
      }
      xcol[i].red = (UShort_t) val;
      val = xcol[i].green + add;
      if (val > kBIGGEST_RGB_VALUE) {
         val = kBIGGEST_RGB_VALUE;
      }
      xcol[i].green = (UShort_t) val;
      val = xcol[i].blue + add;
      if (val > kBIGGEST_RGB_VALUE) {
         val = kBIGGEST_RGB_VALUE;
      }
      xcol[i].blue = (UShort_t) val;

      ret = gdk_color_alloc((GdkColormap *)fColormap, &xcol[i]);

      if (!ret) {
         Warning("MakeOpaqueColors",
                 "failed to allocate color %hd, %hd, %hd", xcol[i].red,
                 xcol[i].green, xcol[i].blue);
      // assumes that in case of failure xcol[i].pixel is not changed
      }
   }

   gCws->new_colors = new ULong_t[ncolors];
   gCws->ncolors = ncolors;

   for (i = 0; i < ncolors; i++) {
      gCws->new_colors[i] = xcol[i].pixel;
   }

   delete []xcol;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns index in orgcolors (and new_colors) for pixel.

Int_t TGWin32::FindColor(ULong_t pixel, ULong_t * orgcolors, Int_t ncolors)
{
   for (int i = 0; i < ncolors; i++) {
      if (pixel == orgcolors[i]) return i;
   }
   Error("FindColor", "did not find color, should never happen!");

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color intensities for given color index.
/// cindex     : color index
/// r,g,b      : red, green, blue intensities between 0.0 and 1.0

void TGWin32::SetRGB(int cindex, float r, float g, float b)
{
   GdkColor xcol;

   if (fColormap && cindex >= 0) {
      xcol.red = (unsigned short) (r * kBIGGEST_RGB_VALUE);
      xcol.green = (unsigned short) (g * kBIGGEST_RGB_VALUE);
      xcol.blue = (unsigned short) (b * kBIGGEST_RGB_VALUE);
      xcol.pixel = RGB(xcol.red, xcol.green, xcol.blue);

      XColor_t &col = GetColor(cindex);
      if (col.fDefined) {
         // if color is already defined with same rgb just return
         if (col.color.red  == xcol.red && col.color.green == xcol.green &&
             col.color.blue == xcol.blue)
            return;
         col.fDefined = kFALSE;
         gdk_colormap_free_colors((GdkColormap *) fColormap,
                                  (GdkColor *)&col, 1);
      }

      Int_t ret = gdk_colormap_alloc_color(fColormap, &xcol, 1, 1);
      if (ret != 0) {
         col.fDefined = kTRUE;
         col.color.pixel   = xcol.pixel;
         col.color.red     = xcol.red;
         col.color.green   = xcol.green;
         col.color.blue    = xcol.blue;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set text alignment.
/// txalh   : horizontal text alignment
/// txalv   : vertical text alignment

void TGWin32::SetTextAlign(Short_t talign)
{
   static Short_t current = 0;
   if (talign==current) return;
   current = talign;

   Int_t txalh = talign / 10;
   Int_t txalv = talign % 10;
   fTextAlignH = txalh;
   fTextAlignV = txalv;

   switch (txalh) {

   case 0:
   case 1:
      switch (txalv) {          //left
      case 1:
         fTextAlign = 7;        //bottom
         break;
      case 2:
         fTextAlign = 4;        //center
         break;
      case 3:
         fTextAlign = 1;        //top
         break;
      }
      break;
   case 2:
      switch (txalv) {          //center
      case 1:
         fTextAlign = 8;        //bottom
         break;
      case 2:
         fTextAlign = 5;        //center
         break;
      case 3:
         fTextAlign = 2;        //top
         break;
      }
      break;
   case 3:
      switch (txalv) {          //right
      case 1:
         fTextAlign = 9;        //bottom
         break;
      case 2:
         fTextAlign = 6;        //center
         break;
      case 3:
         fTextAlign = 3;        //top
         break;
      }
      break;
   }
   TAttText::SetTextAlign(fTextAlign);
}

////////////////////////////////////////////////////////////////////////////////
/// Set color index for text.

void TGWin32::SetTextColor(Color_t cindex)
{
   static Int_t current = 0;
   GdkGCValues values;
   if ((cindex < 0) || (Int_t(cindex)==current)) return;

   TAttText::SetTextColor(cindex);

   SetColor(gGCtext, Int_t(cindex));
   gdk_gc_get_values(gGCtext, &values);
   gdk_gc_set_foreground(gGCinvt, &values.background);
   gdk_gc_set_background(gGCinvt, &values.foreground);
   gdk_gc_set_background(gGCtext, (GdkColor *) & GetColor(0).color);
   current = Int_t(cindex);
}

////////////////////////////////////////////////////////////////////////////////

void TGWin32::Sync(int mode)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Update display.
/// mode : (1) update
///        (0) sync
///
/// Synchronise client and server once (not permanent).
/// Copy the pixmap gCws->drawing on the window gCws->window
/// if the double buffer is on.

void TGWin32::UpdateWindow(int mode)
{
   if (gCws && gCws->double_buffer) {
      gdk_window_copy_area(gCws->window, gGCpxmp, 0, 0,
                           gCws->drawing, 0, 0, gCws->width, gCws->height);
   }
   Update(mode);
}

////////////////////////////////////////////////////////////////////////////////
/// Set pointer position.
/// ix       : New X coordinate of pointer
/// iy       : New Y coordinate of pointer
/// Coordinates are relative to the origin of the window id
/// or to the origin of the current window if id == 0.

void TGWin32::Warp(int ix, int iy, Window_t id)
{
   if (!id) return;

   POINT cpt, tmp;
   HWND dw;
   if (!id)
      dw = (HWND) GDK_DRAWABLE_XID((GdkWindow *)gCws->window);
   else
      dw = (HWND) GDK_DRAWABLE_XID((GdkWindow *)id);
   GetCursorPos(&cpt);
   tmp.x = ix > 0 ? ix : cpt.x;
   tmp.y = iy > 0 ? iy : cpt.y;
   ClientToScreen(dw, &tmp);
   SetCursorPos(tmp.x, tmp.y);
}

////////////////////////////////////////////////////////////////////////////////
/// Write the pixmap wid in the bitmap file pxname.
/// wid         : Pixmap address
/// w,h         : Width and height of the pixmap.
/// lenname     : pixmap name length
/// pxname      : pixmap name

void TGWin32::WritePixmap(int wid, unsigned int w, unsigned int h,
                          char *pxname)
{
   int wval, hval;
   wval = w;
   hval = h;

   if (!fWindows) return;
   gTws = &fWindows[wid];
//   XWriteBitmapFile(fDisplay,pxname,(Pixmap)gTws->drawing,wval,hval,-1,-1);
}


//
// Functions for GIFencode()
//

static FILE *gGifFile;           // output unit used WriteGIF and PutByte
static GdkImage *gGifImage = 0;  // image used in WriteGIF and GetPixel

extern "C" {
   int GIFquantize(UInt_t width, UInt_t height, Int_t * ncol, Byte_t * red,
                   Byte_t * green, Byte_t * blue, Byte_t * outputBuf,
                   Byte_t * outputCmap);
   long GIFencode(int Width, int Height, Int_t Ncol, Byte_t R[],
                  Byte_t G[], Byte_t B[], Byte_t ScLine[],
                  void (*get_scline) (int, int, Byte_t *),
                  void (*pb) (Byte_t));
   int GIFdecode(Byte_t * GIFarr, Byte_t * PIXarr, int *Width, int *Height,
                 int *Ncols, Byte_t * R, Byte_t * G, Byte_t * B);
   int GIFinfo(Byte_t * GIFarr, int *Width, int *Height, int *Ncols);
}


////////////////////////////////////////////////////////////////////////////////
/// Get pixels in line y and put in array scline.

static void GetPixel(int y, int width, Byte_t * scline)
{
   for (int i = 0; i < width; i++) {
       scline[i] = Byte_t(GetPixelImage((Drawable_t)gGifImage, i, y));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Put byte b in output stream.

static void PutByte(Byte_t b)
{
   if (ferror(gGifFile) == 0) fputc(b, gGifFile);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns in R G B the ncol colors of the palette used by the image.
/// The image pixels are changed to index values in these R G B arrays.
/// This produces a colormap with only the used colors (so even on displays
/// with more than 8 planes we will be able to create GIF's when the image
/// contains no more than 256 different colors). If it does contain more
/// colors we will have to use GIFquantize to reduce the number of colors.
/// The R G B arrays must be deleted by the caller.

void TGWin32::ImgPickPalette(GdkImage * image, Int_t & ncol, Int_t * &R,
                             Int_t * &G, Int_t * &B)
{
   ULong_t *orgcolors = 0;
   Int_t maxcolors = 0, ncolors;

   // collect different image colors
   int x, y;
   for (x = 0; x < (int) gCws->width; x++) {
      for (y = 0; y < (int) gCws->height; y++) {
         ULong_t pixel = GetPixelImage((Drawable_t)image, x, y);
         CollectImageColors(pixel, orgcolors, ncolors, maxcolors);
      }
   }

   // get RGB values belonging to pixels
   GdkColor *xcol = new GdkColor[ncolors];

   int i;
   for (i = 0; i < ncolors; i++) {
      xcol[i].pixel = orgcolors[i];
//      xcol[i].red   = xcol[i].green = xcol[i].blue = 0;
      xcol[i].red = GetRValue(xcol[i].pixel);
      xcol[i].green = GetGValue(xcol[i].pixel);
      xcol[i].blue = GetBValue(xcol[i].pixel);
   }

   GdkColorContext *cc;
   cc =  gdk_color_context_new(gdk_visual_get_system(), (GdkColormap *)fColormap);
   gdk_color_context_query_colors(cc, xcol, ncolors);
   gdk_color_context_free(cc);

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
   for (x = 0; x < (int) gCws->width; x++) {
      for (y = 0; y < (int) gCws->height; y++) {
         ULong_t pixel = GetPixelImage((Drawable_t)image, x, y);
         Int_t idx = FindColor(pixel, orgcolors, ncolors);
         PutPixel((Drawable_t)image, x, y, idx);
      }
   }

   // cleanup
   delete[]xcol;
   ::operator delete(orgcolors);
}

////////////////////////////////////////////////////////////////////////////////
/// Writes the current window into GIF file.

Int_t TGWin32::WriteGIF(char *name)
{
   Byte_t scline[2000], r[256], b[256], g[256];
   Int_t *R, *G, *B;
   Int_t ncol, maxcol, i;

   if (gGifImage) {
      gdk_image_unref((GdkImage *)gGifImage);
   }

   gGifImage = gdk_image_get((GdkDrawable*)gCws->drawing, 0, 0,
                             gCws->width, gCws->height);

   ImgPickPalette(gGifImage, ncol, R, G, B);

   if (ncol > 256) {
      //GIFquantize(...);
      Error("WriteGIF",
            "can not create GIF of image containing more than 256 colors");
      delete[]R;
      delete[]G;
      delete[]B;
      return 0;
   }

   maxcol = 0;
   for (i = 0; i < ncol; i++) {
      if (maxcol < R[i]) maxcol = R[i];
      if (maxcol < G[i]) maxcol = G[i];
      if (maxcol < B[i]) maxcol = B[i];
      r[i] = 0;
      g[i] = 0;
      b[i] = 0;
   }
   if (maxcol != 0) {
      for (i = 0; i < ncol; i++) {
         r[i] = R[i] * 255 / maxcol;
         g[i] = G[i] * 255 / maxcol;
         b[i] = B[i] * 255 / maxcol;
      }
   }

   gGifFile = fopen(name, "wb");

   if (gGifFile) {
      GIFencode(gCws->width, gCws->height,
          ncol, r, g, b, scline, ::GetPixel, PutByte);
      fclose(gGifFile);
      i = 1;
    } else {
      Error("WriteGIF","cannot write file: %s",name);
      i = 0;
   }
   delete[]R;
   delete[]G;
   delete[]B;

   return i;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw image.

void TGWin32::PutImage(int offset, int itran, int x0, int y0, int nx,
                       int ny, int xmin, int ymin, int xmax, int ymax,
                       unsigned char *image, Drawable_t wid)
{
   const int MAX_SEGMENT = 20;
   int i, n, x, y, xcur, x1, x2, y1, y2;
   unsigned char *jimg, *jbase, icol;
   int nlines[256];
   GdkSegment lines[256][MAX_SEGMENT];
   GdkDrawable *id;

   if (wid) {
      id = (GdkDrawable*)wid;
   } else {
      id = gCws->drawing;
   }

   for (i = 0; i < 256; i++) nlines[i] = 0;

   x1 = x0 + xmin;
   y1 = y0 + ny - ymax - 1;
   x2 = x0 + xmax;
   y2 = y0 + ny - ymin - 1;
   jbase = image + (ymin - 1) * nx + xmin;

   for (y = y2; y >= y1; y--) {
      xcur = x1;
      jbase += nx;
      for (jimg = jbase, icol = *jimg++, x = x1 + 1; x <= x2; jimg++, x++) {
         if (icol != *jimg) {
            if (icol != itran) {
               n = nlines[icol]++;
               lines[icol][n].x1 = xcur;
               lines[icol][n].y1 = y;
               lines[icol][n].x2 = x - 1;
               lines[icol][n].y2 = y;
               if (nlines[icol] == MAX_SEGMENT) {
                  SetColor(gGCline, (int) icol + offset);
                  gdk_win32_draw_segments(id, (GdkGC *) gGCline,
                                       (GdkSegment *) &lines[icol][0], MAX_SEGMENT);
                  nlines[icol] = 0;
               }
            }
            icol = *jimg;
            xcur = x;
         }
      }
      if (icol != itran) {
         n = nlines[icol]++;
         lines[icol][n].x1 = xcur;
         lines[icol][n].y1 = y;
         lines[icol][n].x2 = x - 1;
         lines[icol][n].y2 = y;
         if (nlines[icol] == MAX_SEGMENT) {
            SetColor(gGCline, (int) icol + offset);
            gdk_win32_draw_segments(id, (GdkGC *) gGCline,
                              (GdkSegment *)&lines[icol][0], MAX_SEGMENT);
            nlines[icol] = 0;
         }
      }
   }

   for (i = 0; i < 256; i++) {
      if (nlines[i] != 0) {
         SetColor(gGCline, i + offset);
         gdk_win32_draw_segments(id, (GdkGC *) gGCline,
                           (GdkSegment *)&lines[icol][0], nlines[i]);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If id is NULL - loads the specified gif file at position [x0,y0] in the
/// current window. Otherwise creates pixmap from gif file

Pixmap_t TGWin32::ReadGIF(int x0, int y0, const char *file, Window_t id)
{
   FILE *fd;
   Seek_t filesize;
   unsigned char *GIFarr, *PIXarr, R[256], G[256], B[256], *j1, *j2, icol;
   int i, j, k, width, height, ncolor, irep, offset;
   float rr, gg, bb;
   Pixmap_t pic = 0;

   fd = fopen(file, "r+b");
   if (!fd) {
      Error("ReadGIF", "unable to open GIF file");
      return pic;
   }

   fseek(fd, 0L, 2);
   filesize = Seek_t(ftell(fd));
   fseek(fd, 0L, 0);

   if (!(GIFarr = (unsigned char *) calloc(filesize + 256, 1))) {
      fclose(fd);
      Error("ReadGIF", "unable to allocate array for gif");
      return pic;
   }

   if (fread(GIFarr, filesize, 1, fd) != 1) {
      fclose(fd);
      Error("ReadGIF", "GIF file read failed");
      free(GIFarr);
      return pic;
   }
   fclose(fd);

   irep = GIFinfo(GIFarr, &width, &height, &ncolor);
   if (irep != 0) {
      return pic;
   }

   if (!(PIXarr = (unsigned char *) calloc((width * height), 1))) {
      Error("ReadGIF", "unable to allocate array for image");
      return pic;
   }

   irep = GIFdecode(GIFarr, PIXarr, &width, &height, &ncolor, R, G, B);
   if (irep != 0) {
      return pic;
   }
   // S E T   P A L E T T E

   offset = 8;

   for (i = 0; i < ncolor; i++) {
      rr = R[i] / 255.;
      gg = G[i] / 255.;
      bb = B[i] / 255.;
      j = i + offset;
      SetRGB(j, rr, gg, bb);
   }

   // O U T P U T   I M A G E

   for (i = 1; i <= height / 2; i++) {
      j1 = PIXarr + (i - 1) * width;
      j2 = PIXarr + (height - i) * width;
      for (k = 0; k < width; k++) {
         icol = *j1;
         *j1++ = *j2;
         *j2++ = icol;
      }
   }

   if (id) pic = CreatePixmap(id, width, height);
   PutImage(offset, -1, x0, y0, width, height, 0, 0, width-1, height-1, PIXarr, pic);

   if (pic) return pic;
   else if (gCws->drawing) return  (Pixmap_t)gCws->drawing;
   else return 0;
}

//////////////////////////// GWin32Gui //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/// Map window on screen.

void TGWin32::MapWindow(Window_t id)
{
   if (!id) return;

   gdk_window_show((GdkWindow *)id);
   if ((GDK_DRAWABLE_TYPE((GdkWindow *)id) != GDK_WINDOW_TEMP) &&
       (GetParent(id) == GetDefaultRootWindow())) {
      HWND window = (HWND)GDK_DRAWABLE_XID((GdkWindow *)id);
      ::SetForegroundWindow(window);
   }
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::MapSubwindows(Window_t id)
{
   if (!id) return;

   EnumChildWindows((HWND)GDK_DRAWABLE_XID((GdkWindow *)id),
                    EnumChildProc, (LPARAM) NULL);
}

////////////////////////////////////////////////////////////////////////////////
/// Map window on screen and put on top of all windows.

void TGWin32::MapRaised(Window_t id)
{
   if (!id) return;

   HWND hwnd = ::GetForegroundWindow();
   HWND window = (HWND)GDK_DRAWABLE_XID((GdkWindow *)id);
   gdk_window_show((GdkWindow *)id);
   if (GDK_DRAWABLE_TYPE((GdkWindow *)id) != GDK_WINDOW_TEMP) {
      ::BringWindowToTop(window);
      if (GDK_DRAWABLE_TYPE((GdkWindow *)id) != GDK_WINDOW_CHILD)
         ::SetForegroundWindow(window);
   }

   if (gConsoleWindow && (hwnd == (HWND)gConsoleWindow)) {
      RECT r1, r2, r3;
      ::GetWindowRect((HWND)gConsoleWindow, &r1);
      HWND fore = ::GetForegroundWindow();
      ::GetWindowRect(fore, &r2);
      if (!::IntersectRect(&r3, &r2, &r1)) {
         ::SetForegroundWindow((HWND)gConsoleWindow);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Unmap window from screen.

void TGWin32::UnmapWindow(Window_t id)
{
   if (!id) return;

   gdk_window_hide((GdkWindow *) id);
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy window.

void TGWin32::DestroyWindow(Window_t id)
{
   if (!id) return;

   // we need to unmap the window before to destroy it, in order to properly
   // receive kUnmapNotify needed by gClient->WaitForUnmap()...
   gdk_window_hide((GdkWindow *) id);
   gdk_window_destroy((GdkDrawable *) id, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy all internal subwindows

void TGWin32::DestroySubwindows(Window_t id)
{
   if (!id) return;

   gdk_window_destroy((GdkDrawable *) id, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Put window on top of window stack.

void TGWin32::RaiseWindow(Window_t id)
{
   if (!id) return;

   HWND window = (HWND)GDK_DRAWABLE_XID((GdkWindow *)id);
   if (GDK_DRAWABLE_TYPE((GdkWindow *)id) == GDK_WINDOW_TEMP) {
       ::SetWindowPos(window, HWND_TOPMOST,  0, 0, 0, 0,
                      SWP_NOACTIVATE | SWP_NOMOVE | SWP_NOSIZE);
   }
   else {
      ::BringWindowToTop(window);
      if (GDK_DRAWABLE_TYPE((GdkWindow *)id) != GDK_WINDOW_CHILD)
         ::SetForegroundWindow(window);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Lower window so it lays below all its siblings.

void TGWin32::LowerWindow(Window_t id)
{
   if (!id) return;

   HWND window = (HWND)GDK_DRAWABLE_XID((GdkWindow *)id);
   ::SetWindowPos(window, HWND_BOTTOM, 0, 0, 0, 0,
                  SWP_NOACTIVATE | SWP_NOMOVE | SWP_NOSIZE);
}

////////////////////////////////////////////////////////////////////////////////
/// Move a window.

void TGWin32::MoveWindow(Window_t id, Int_t x, Int_t y)
{
   if (!id) return;

   gdk_window_move((GdkDrawable *) id, x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Move and resize a window.

void TGWin32::MoveResizeWindow(Window_t id, Int_t x, Int_t y, UInt_t w,
                               UInt_t h)
{
   if (!id) return;

   gdk_window_move_resize((GdkWindow *) id, x, y, w, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Resize the window.

void TGWin32::ResizeWindow(Window_t id, UInt_t w, UInt_t h)
{
   if (!id) return;

   // protect against potential negative values
   if (w >= (UInt_t)INT_MAX || h >= (UInt_t)INT_MAX)
      return;
   gdk_window_resize((GdkWindow *) id, w, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Iconify the window.

void TGWin32::IconifyWindow(Window_t id)
{
   if (!id) return;

   gdk_window_lower((GdkWindow *) id);
   ::CloseWindow((HWND)GDK_DRAWABLE_XID((GdkWindow *)id));
}

////////////////////////////////////////////////////////////////////////////////
/// Reparent window, make pid the new parent and position the window at
/// position (x,y) in new parent.

void TGWin32::ReparentWindow(Window_t id, Window_t pid, Int_t x, Int_t y)
{
   if (!id) return;

   gdk_window_reparent((GdkWindow *)id, (GdkWindow *)pid, x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the window background color.

void TGWin32::SetWindowBackground(Window_t id, ULong_t color)
{
   if (!id) return;

   GdkColor back;
   back.pixel = color;
   back.red = GetRValue(color);
   back.green = GetGValue(color);
   back.blue = GetBValue(color);

   gdk_window_set_background((GdkWindow *) id, &back);
}

////////////////////////////////////////////////////////////////////////////////
/// Set pixmap as window background.

void TGWin32::SetWindowBackgroundPixmap(Window_t id, Pixmap_t pxm)
{
   if (!id) return;

   gdk_window_set_back_pixmap((GdkWindow *) id, (GdkPixmap *) pxm, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Return handle to newly created gdk window.

Window_t TGWin32::CreateWindow(Window_t parent, Int_t x, Int_t y,
                               UInt_t w, UInt_t h, UInt_t border,
                               Int_t depth, UInt_t clss,
                               void *visual, SetWindowAttributes_t * attr,
                               UInt_t wtype)
{
   GdkWindowAttr xattr;
   GdkWindow *newWin;
   GdkColor background_color;
   ULong_t xmask = 0;

   if (attr) {
      MapSetWindowAttributes(attr, xmask, xattr);
      xattr.window_type = GDK_WINDOW_CHILD;
      if (wtype & kMainFrame) {
         xattr.window_type = GDK_WINDOW_TOPLEVEL;
      }
      if (wtype & kTransientFrame) {
         xattr.window_type = GDK_WINDOW_DIALOG;
      }
      if (wtype & kTempFrame) {
         xattr.window_type = GDK_WINDOW_TEMP;
      }
      newWin = gdk_window_new((GdkWindow *) parent, &xattr, xmask);
   } else {
      xattr.width = w;
      xattr.height = h;
      xattr.wclass = GDK_INPUT_OUTPUT;
      xattr.event_mask = 0L;    //GDK_ALL_EVENTS_MASK;
      xattr.event_mask |= GDK_EXPOSURE_MASK | GDK_STRUCTURE_MASK |
          GDK_PROPERTY_CHANGE_MASK;
//                          GDK_ENTER_NOTIFY_MASK | GDK_LEAVE_NOTIFY_MASK;
      if (x >= 0) {
         xattr.x = x;
      } else {
         xattr.x = -1.0 * x;
      }
      if (y >= 0) {
         xattr.y = y;
      } else {
         xattr.y = -1.0 * y;
      }
      xattr.colormap = gdk_colormap_get_system();
      xattr.cursor = NULL;
      xattr.override_redirect = TRUE;
      if ((xattr.y > 0) && (xattr.x > 0)) {
         xmask = GDK_WA_X | GDK_WA_Y | GDK_WA_COLORMAP |
             GDK_WA_WMCLASS | GDK_WA_NOREDIR;
      } else {
         xmask = GDK_WA_COLORMAP | GDK_WA_WMCLASS | GDK_WA_NOREDIR;
      }
      if (visual != NULL) {
         xattr.visual = (GdkVisual *) visual;
         xmask |= GDK_WA_VISUAL;
      } else {
         xattr.visual = gdk_visual_get_system();
         xmask |= GDK_WA_VISUAL;
      }
      xattr.window_type = GDK_WINDOW_CHILD;
      if (wtype & kMainFrame) {
         xattr.window_type = GDK_WINDOW_TOPLEVEL;
      }
      if (wtype & kTransientFrame) {
         xattr.window_type = GDK_WINDOW_DIALOG;
      }
      if (wtype & kTempFrame) {
         xattr.window_type = GDK_WINDOW_TEMP;
      }
      newWin = gdk_window_new((GdkWindow *) parent, &xattr, xmask);
      gdk_window_set_events(newWin, (GdkEventMask) 0L);
   }
   if (border > 0) {
      gdk_window_set_decorations(newWin,
                                 (GdkWMDecoration) GDK_DECOR_BORDER);
   }
   if (attr) {
      if ((attr->fMask & kWABackPixmap)) {
         if (attr->fBackgroundPixmap == kNone) {
            gdk_window_set_back_pixmap(newWin, (GdkPixmap *) GDK_NONE, 0);
         } else if (attr->fBackgroundPixmap == kParentRelative) {
            gdk_window_set_back_pixmap(newWin, (GdkPixmap *) GDK_NONE, 1);
         } else {
            gdk_window_set_back_pixmap(newWin,
                                       (GdkPixmap *) attr->
                                       fBackgroundPixmap, 0);
         }
      }
      if ((attr->fMask & kWABackPixel)) {
         background_color.pixel = attr->fBackgroundPixel;
         background_color.red = GetRValue(attr->fBackgroundPixel);
         background_color.green = GetGValue(attr->fBackgroundPixel);
         background_color.blue = GetBValue(attr->fBackgroundPixel);
         gdk_window_set_background(newWin, &background_color);
      }
   }
   if (!fUseSysPointers) {
      ::SetClassLong((HWND)GDK_DRAWABLE_XID(newWin), GCL_HCURSOR,
                     (LONG)GDK_CURSOR_XID(fCursors[kPointer]));
   }
   return (Window_t) newWin;
}

////////////////////////////////////////////////////////////////////////////////
/// Map event mask to or from gdk.

void TGWin32::MapEventMask(UInt_t & emask, UInt_t & xemask, Bool_t tox)
{
   if (tox) {
      Long_t lxemask = 0L;
      if ((emask & kKeyPressMask)) {
         lxemask |= GDK_KEY_PRESS_MASK;
      }
      if ((emask & kKeyReleaseMask)) {
         lxemask |= GDK_KEY_RELEASE_MASK;
      }
      if ((emask & kButtonPressMask)) {
         lxemask |= GDK_BUTTON_PRESS_MASK;
      }
      if ((emask & kButtonReleaseMask)) {
         lxemask |= GDK_BUTTON_RELEASE_MASK;
      }
      if ((emask & kPointerMotionMask)) {
         lxemask |= GDK_POINTER_MOTION_MASK;
      }
      if ((emask & kButtonMotionMask)) {
         lxemask |= GDK_BUTTON_MOTION_MASK;
      }
      if ((emask & kExposureMask)) {
         lxemask |= GDK_EXPOSURE_MASK;
      }
      if ((emask & kStructureNotifyMask)) {
         lxemask |= GDK_STRUCTURE_MASK;
      }
      if ((emask & kEnterWindowMask)) {
         lxemask |= GDK_ENTER_NOTIFY_MASK;
      }
      if ((emask & kLeaveWindowMask)) {
         lxemask |= GDK_LEAVE_NOTIFY_MASK;
      }
      if ((emask & kFocusChangeMask)) {
         lxemask |= GDK_FOCUS_CHANGE_MASK;
      }
      xemask = (UInt_t) lxemask;
   } else {
      emask = 0;
      if ((xemask & GDK_KEY_PRESS_MASK)) {
         emask |= kKeyPressMask;
      }
      if ((xemask & GDK_KEY_RELEASE_MASK)) {
         emask |= kKeyReleaseMask;
      }
      if ((xemask & GDK_BUTTON_PRESS_MASK)) {
         emask |= kButtonPressMask;
      }
      if ((xemask & GDK_BUTTON_RELEASE_MASK)) {
         emask |= kButtonReleaseMask;
      }
      if ((xemask & GDK_POINTER_MOTION_MASK)) {
         emask |= kPointerMotionMask;
      }
      if ((xemask & GDK_BUTTON_MOTION_MASK)) {
         emask |= kButtonMotionMask;
      }
      if ((xemask & GDK_EXPOSURE_MASK)) {
         emask |= kExposureMask;
      }
      if ((xemask & GDK_STRUCTURE_MASK)) {
         emask |= kStructureNotifyMask;
      }
      if ((xemask & GDK_ENTER_NOTIFY_MASK)) {
         emask |= kEnterWindowMask;
      }
      if ((xemask & GDK_LEAVE_NOTIFY_MASK)) {
         emask |= kLeaveWindowMask;
      }
      if ((xemask & GDK_FOCUS_CHANGE_MASK)) {
         emask |= kFocusChangeMask;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Map a SetWindowAttributes_t to a GdkWindowAttr structure.

void TGWin32::MapSetWindowAttributes(SetWindowAttributes_t * attr,
                                     ULong_t & xmask,
                                     GdkWindowAttr & xattr)
{
   Mask_t mask = attr->fMask;
   xmask = 0;

   if ((mask & kWAOverrideRedirect)) {
      xmask |= GDK_WA_NOREDIR;
      xattr.override_redirect = attr->fOverrideRedirect;
   }
   if ((mask & kWAEventMask)) {
      UInt_t xmsk, msk = (UInt_t) attr->fEventMask;
      MapEventMask(msk, xmsk, kTRUE);
      xattr.event_mask = xmsk;
   }
   if ((mask & kWAColormap)) {
      xmask |= GDK_WA_COLORMAP;
      xattr.colormap = (GdkColormap *) attr->fColormap;
   }
   if ((mask & kWACursor)) {
      xmask |= GDK_WA_CURSOR;
      if (attr->fCursor != kNone) {
         xattr.cursor = (GdkCursor *) attr->fCursor;
      }
   }
   xattr.wclass = GDK_INPUT_OUTPUT;
}

////////////////////////////////////////////////////////////////////////////////
/// Map a GCValues_t to a XCGValues structure if tox is true. Map
/// the other way in case tox is false.

void TGWin32::MapGCValues(GCValues_t & gval,
                          ULong_t & xmask, GdkGCValues & xgval, Bool_t tox)
{
   if (tox) {
      // map GCValues_t to XGCValues
      Mask_t mask = gval.fMask;
      xmask = 0;

      if ((mask & kGCFunction)) {
         xmask |= GDK_GC_FUNCTION;
         switch (gval.fFunction) {
         case kGXclear:
            xgval.function = GDK_CLEAR;
            break;
         case kGXand:
            xgval.function = GDK_AND;
            break;
         case kGXandReverse:
            xgval.function = GDK_AND_REVERSE;
            break;
         case kGXcopy:
            xgval.function = GDK_COPY;
            break;
         case kGXandInverted:
            xgval.function = GDK_AND_INVERT;
            break;
         case kGXnoop:
            xgval.function = GDK_NOOP;
            break;
         case kGXxor:
            xgval.function = GDK_XOR;
            break;
         case kGXor:
            xgval.function = GDK_OR;
            break;
         case kGXequiv:
            xgval.function = GDK_EQUIV;
            break;
         case kGXinvert:
            xgval.function = GDK_INVERT;
            break;
         case kGXorReverse:
            xgval.function = GDK_OR_REVERSE;
            break;
         case kGXcopyInverted:
            xgval.function = GDK_COPY_INVERT;
            break;
         case kGXorInverted:
            xgval.function = GDK_OR_INVERT;
            break;
         case kGXnand:
            xgval.function = GDK_NAND;
            break;
         case kGXset:
            xgval.function = GDK_SET;
            break;
         }
      }
      if (mask & kGCSubwindowMode) {
         xmask |= GDK_GC_SUBWINDOW;
         if (gval.fSubwindowMode == kIncludeInferiors) {
            xgval.subwindow_mode = GDK_INCLUDE_INFERIORS;
         } else {
            xgval.subwindow_mode = GDK_CLIP_BY_CHILDREN;
         }
      }
      if (mask & kGCForeground) {
         xmask |= GDK_GC_FOREGROUND;
         xgval.foreground.pixel = gval.fForeground;
         xgval.foreground.red = GetRValue(gval.fForeground);
         xgval.foreground.green = GetGValue(gval.fForeground);
         xgval.foreground.blue = GetBValue(gval.fForeground);
      }
      if (mask & kGCBackground) {
         xmask |= GDK_GC_BACKGROUND;
         xgval.background.pixel = gval.fBackground;
         xgval.background.red = GetRValue(gval.fBackground);
         xgval.background.green = GetGValue(gval.fBackground);
         xgval.background.blue = GetBValue(gval.fBackground);
      }
      if (mask & kGCLineWidth) {
         xmask |= GDK_GC_LINE_WIDTH;
         xgval.line_width = gval.fLineWidth;
      }
      if (mask & kGCLineStyle) {
         xmask |= GDK_GC_LINE_STYLE;
         xgval.line_style = (GdkLineStyle) gval.fLineStyle; // ident mapping
      }
      if (mask & kGCCapStyle) {
         xmask |= GDK_GC_CAP_STYLE;
         xgval.cap_style = (GdkCapStyle) gval.fCapStyle; // ident mapping
      }
      if (mask & kGCJoinStyle) {
         xmask |= GDK_GC_JOIN_STYLE;
         xgval.join_style = (GdkJoinStyle) gval.fJoinStyle; // ident mapping
      }
      if ((mask & kGCFillStyle)) {
         xmask |= GDK_GC_FILL;
         xgval.fill = (GdkFill) gval.fFillStyle;   // ident mapping
      }
      if ((mask & kGCTile)) {
         xmask |= GDK_GC_TILE;
         xgval.tile = (GdkPixmap *) gval.fTile;
      }
      if ((mask & kGCStipple)) {
         xmask |= GDK_GC_STIPPLE;
         xgval.stipple = (GdkPixmap *) gval.fStipple;
      }
      if ((mask & kGCTileStipXOrigin)) {
         xmask |= GDK_GC_TS_X_ORIGIN;
         xgval.ts_x_origin = gval.fTsXOrigin;
      }
      if ((mask & kGCTileStipYOrigin)) {
         xmask |= GDK_GC_TS_Y_ORIGIN;
         xgval.ts_y_origin = gval.fTsYOrigin;
      }
      if ((mask & kGCFont)) {
         xmask |= GDK_GC_FONT;
         xgval.font = (GdkFont *) gval.fFont;
      }
      if ((mask & kGCGraphicsExposures)) {
         xmask |= GDK_GC_EXPOSURES;
         xgval.graphics_exposures = gval.fGraphicsExposures;
      }
      if ((mask & kGCClipXOrigin)) {
         xmask |= GDK_GC_CLIP_X_ORIGIN;
         xgval.clip_x_origin = gval.fClipXOrigin;
      }
      if ((mask & kGCClipYOrigin)) {
         xmask |= GDK_GC_CLIP_Y_ORIGIN;
         xgval.clip_y_origin = gval.fClipYOrigin;
      }
      if ((mask & kGCClipMask)) {
         xmask |= GDK_GC_CLIP_MASK;
         xgval.clip_mask = (GdkPixmap *) gval.fClipMask;
      }
   } else {
      // map XValues to GCValues_t
      Mask_t mask = 0;

      if ((xmask & GDK_GC_FUNCTION)) {
         mask |= kGCFunction;
         gval.fFunction = (EGraphicsFunction) xgval.function; // ident mapping
         switch (xgval.function) {
         case GDK_CLEAR:
            gval.fFunction = kGXclear;
            break;
         case GDK_AND:
            gval.fFunction = kGXand;
            break;
         case GDK_AND_REVERSE:
            gval.fFunction = kGXandReverse;
            break;
         case GDK_COPY:
            gval.fFunction = kGXcopy;
            break;
         case GDK_AND_INVERT:
            gval.fFunction = kGXandInverted;
            break;
         case GDK_NOOP:
            gval.fFunction = kGXnoop;
            break;
         case GDK_XOR:
            gval.fFunction = kGXxor;
            break;
         case GDK_OR:
            gval.fFunction = kGXor;
            break;
         case GDK_EQUIV:
            gval.fFunction = kGXequiv;
            break;
         case GDK_INVERT:
            gval.fFunction = kGXinvert;
            break;
         case GDK_OR_REVERSE:
            gval.fFunction = kGXorReverse;
            break;
         case GDK_COPY_INVERT:
            gval.fFunction = kGXcopyInverted;
            break;
         case GDK_OR_INVERT:
            gval.fFunction = kGXorInverted;
            break;
         case GDK_NAND:
            gval.fFunction = kGXnand;
            break;
         case GDK_SET:
            gval.fFunction = kGXset;
            break;
         }
      }
      if (xmask & GDK_GC_SUBWINDOW) {
         mask |= kGCSubwindowMode;
         if (xgval.subwindow_mode == GDK_INCLUDE_INFERIORS)
            gval.fSubwindowMode = kIncludeInferiors;
         else
            gval.fSubwindowMode = kClipByChildren;
      }
      if ((xmask & GDK_GC_FOREGROUND)) {
         mask |= kGCForeground;
         gval.fForeground = xgval.foreground.pixel;
      }
      if ((xmask & GDK_GC_BACKGROUND)) {
         mask |= kGCBackground;
         gval.fBackground = xgval.background.pixel;
      }
      if ((xmask & GDK_GC_LINE_WIDTH)) {
         mask |= kGCLineWidth;
         gval.fLineWidth = xgval.line_width;
      }
      if ((xmask & GDK_GC_LINE_STYLE)) {
         mask |= kGCLineStyle;
         gval.fLineStyle = xgval.line_style; // ident mapping
      }
      if ((xmask & GDK_GC_CAP_STYLE)) {
         mask |= kGCCapStyle;
         gval.fCapStyle = xgval.cap_style;   // ident mapping
      }
      if ((xmask & GDK_GC_JOIN_STYLE)) {
         mask |= kGCJoinStyle;
         gval.fJoinStyle = xgval.join_style; // ident mapping
      }
      if ((xmask & GDK_GC_FILL)) {
         mask |= kGCFillStyle;
         gval.fFillStyle = xgval.fill;       // ident mapping
      }
      if ((xmask & GDK_GC_TILE)) {
         mask |= kGCTile;
         gval.fTile = (Pixmap_t) xgval.tile;
      }
      if ((xmask & GDK_GC_STIPPLE)) {
         mask |= kGCStipple;
         gval.fStipple = (Pixmap_t) xgval.stipple;
      }
      if ((xmask & GDK_GC_TS_X_ORIGIN)) {
         mask |= kGCTileStipXOrigin;
         gval.fTsXOrigin = xgval.ts_x_origin;
      }
      if ((xmask & GDK_GC_TS_Y_ORIGIN)) {
         mask |= kGCTileStipYOrigin;
         gval.fTsYOrigin = xgval.ts_y_origin;
      }
      if ((xmask & GDK_GC_FONT)) {
         mask |= kGCFont;
         gval.fFont = (FontH_t) xgval.font;
      }
      if ((xmask & GDK_GC_EXPOSURES)) {
         mask |= kGCGraphicsExposures;
         gval.fGraphicsExposures = (Bool_t) xgval.graphics_exposures;
      }
      if ((xmask & GDK_GC_CLIP_X_ORIGIN)) {
         mask |= kGCClipXOrigin;
         gval.fClipXOrigin = xgval.clip_x_origin;
      }
      if ((xmask & GDK_GC_CLIP_Y_ORIGIN)) {
         mask |= kGCClipYOrigin;
         gval.fClipYOrigin = xgval.clip_y_origin;
      }
      if ((xmask & GDK_GC_CLIP_MASK)) {
         mask |= kGCClipMask;
         gval.fClipMask = (Pixmap_t) xgval.clip_mask;
      }
      gval.fMask = mask;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get window attributes and return filled in attributes structure.

void TGWin32::GetWindowAttributes(Window_t id, WindowAttributes_t & attr)
{
   if (!id) return;

   RECT rcClient, rcWind;
   ::GetClientRect((HWND)GDK_DRAWABLE_XID((GdkWindow *) id), &rcClient);
   ::GetWindowRect((HWND)GDK_DRAWABLE_XID((GdkWindow *) id), &rcWind);

   gdk_window_get_geometry((GdkWindow *) id, &attr.fX, &attr.fY,
                           &attr.fWidth, &attr.fHeight, &attr.fDepth);
   attr.fX = ((rcWind.right - rcWind.left) - rcClient.right) / 2;
   attr.fY = ((rcWind.bottom - rcWind.top) - rcClient.bottom) - attr.fX;

   attr.fRoot = (Window_t) GDK_ROOT_PARENT();
   attr.fColormap = (Colormap_t) gdk_window_get_colormap((GdkWindow *) id);
   attr.fBorderWidth = 0;
   attr.fVisual = gdk_window_get_visual((GdkWindow *) id);
   attr.fClass = kInputOutput;
   attr.fBackingStore = kNotUseful;
   attr.fSaveUnder = kFALSE;
   attr.fMapInstalled = kTRUE;
   attr.fOverrideRedirect = kFALSE;   // boolean value for override-redirect

   if (!gdk_window_is_visible((GdkWindow *) id)) {
      attr.fMapState = kIsUnmapped;
   } else if (!gdk_window_is_viewable((GdkWindow *) id)) {
      attr.fMapState = kIsUnviewable;
   } else {
      attr.fMapState = kIsViewable;
   }

   UInt_t tmp_mask = (UInt_t)gdk_window_get_events((GdkWindow *) id);
   UInt_t evmask;
   MapEventMask(evmask, tmp_mask, kFALSE);

   attr.fYourEventMask = evmask;
}

////////////////////////////////////////////////////////////////////////////////
///

Display_t TGWin32::GetDisplay() const
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get maximum number of planes.

Int_t TGWin32::GetDepth() const
{
   return gdk_visual_get_best_depth();
}

////////////////////////////////////////////////////////////////////////////////
/// Return atom handle for atom_name. If it does not exist
/// create it if only_if_exist is false. Atoms are used to communicate
/// between different programs (i.e. window manager) via the X server.

Atom_t TGWin32::InternAtom(const char *atom_name, Bool_t only_if_exist)
{
   GdkAtom a = gdk_atom_intern((const gchar *) atom_name, only_if_exist);

   if (a == None) return kNone;
   return (Atom_t) a;
}

////////////////////////////////////////////////////////////////////////////////
/// Return handle to the default root window created when calling
/// XOpenDisplay().

Window_t TGWin32::GetDefaultRootWindow() const
{
   return (Window_t) GDK_ROOT_PARENT();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the parent of the window.

Window_t TGWin32::GetParent(Window_t id) const
{
   if (!id) return (Window_t)0;

   return (Window_t)gdk_window_get_parent((GdkWindow *) id);
}

////////////////////////////////////////////////////////////////////////////////
/// Load font and query font. If font is not found 0 is returned,
/// otherwise an opaque pointer to the FontStruct_t.
/// Free the loaded font using DeleteFont().

FontStruct_t TGWin32::LoadQueryFont(const char *font_name)
{
   char  family[100], weight[32], slant[32], fontname[256];
   Int_t n1, pixel, numfields;

   numfields = sscanf(font_name, "%s -%d%n", family, &pixel, &n1);
   if (numfields == 2) {
      sprintf(weight,"medium");
      if (strstr(font_name, "bold"))
         sprintf(weight,"bold");
      sprintf(slant,"r");
      if (strstr(font_name, "italic"))
         sprintf(slant,"i");
      sprintf(fontname, "-*-%s-%s-%s-*-*-%d-*-*-*-*-*-iso8859-1",
              family, weight, slant, pixel);
   }
   else
      sprintf(fontname, "%s", font_name);
   return (FontStruct_t) gdk_font_load(fontname);
}

////////////////////////////////////////////////////////////////////////////////
/// Return handle to font described by font structure.

FontH_t TGWin32::GetFontHandle(FontStruct_t fs)
{
   if (fs) {
      return (FontH_t)gdk_font_ref((GdkFont *) fs);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Explicitely delete font structure obtained with LoadQueryFont().

void TGWin32::DeleteFont(FontStruct_t fs)
{
   gdk_font_unref((GdkFont *) fs);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a graphics context using the values set in gval (but only for
/// those entries that are in the mask).

GContext_t TGWin32::CreateGC(Drawable_t id, GCValues_t *gval)
{
   if (!id) return (GContext_t)0;

   GdkGCValues xgval;
   ULong_t xmask = 0;

   if (gval) MapGCValues(*gval, xmask, xgval, kTRUE);

   xgval.subwindow_mode = GDK_CLIP_BY_CHILDREN; // GDK_INCLUDE_INFERIORS;

   GdkGC *gc = gdk_gc_new_with_values((GdkDrawable *) id,
                                      &xgval, (GdkGCValuesMask)xmask);
   return (GContext_t) gc;
}

////////////////////////////////////////////////////////////////////////////////
/// Change entries in an existing graphics context, gc, by values from gval.

void TGWin32::ChangeGC(GContext_t gc, GCValues_t * gval)
{
   GdkGCValues xgval;
   ULong_t xmask = 0;
   Mask_t mask = 0;

   if (gval) {
      mask = gval->fMask;
      MapGCValues(*gval, xmask, xgval, kTRUE);
   }
   if (mask & kGCForeground) {
      gdk_gc_set_foreground((GdkGC *) gc, &xgval.foreground);
   }
   if (mask & kGCBackground) {
      gdk_gc_set_background((GdkGC *) gc, &xgval.background);
   }
   if (mask & kGCFont) {
      gdk_gc_set_font((GdkGC *) gc, xgval.font);
   }
   if (mask & kGCFunction) {
      gdk_gc_set_function((GdkGC *) gc, xgval.function);
   }
   if (mask & kGCFillStyle) {
      gdk_gc_set_fill((GdkGC *) gc, xgval.fill);
   }
   if (mask & kGCTile) {
      gdk_gc_set_tile((GdkGC *) gc, xgval.tile);
   }
   if (mask & kGCStipple) {
      gdk_gc_set_stipple((GdkGC *) gc, xgval.stipple);
   }
   if ((mask & kGCTileStipXOrigin) || (mask & kGCTileStipYOrigin)) {
      gdk_gc_set_ts_origin((GdkGC *) gc, xgval.ts_x_origin,
                           xgval.ts_y_origin);
   }
   if ((mask & kGCClipXOrigin) || (mask & kGCClipYOrigin)) {
      gdk_gc_set_clip_origin((GdkGC *) gc, xgval.clip_x_origin,
                             xgval.clip_y_origin);
   }
   if (mask & kGCClipMask) {
      gdk_gc_set_clip_mask((GdkGC *) gc, xgval.clip_mask);
   }
   if (mask & kGCGraphicsExposures) {
      gdk_gc_set_exposures((GdkGC *) gc, xgval.graphics_exposures);
   }
   if (mask & kGCLineWidth) {
      gdk_gc_set_values((GdkGC *) gc, &xgval, GDK_GC_LINE_WIDTH);
   }
   if (mask & kGCLineStyle) {
      gdk_gc_set_values((GdkGC *) gc, &xgval, GDK_GC_LINE_STYLE);
   }
   if (mask & kGCCapStyle) {
      gdk_gc_set_values((GdkGC *) gc, &xgval, GDK_GC_CAP_STYLE);
   }
   if (mask & kGCJoinStyle) {
      gdk_gc_set_values((GdkGC *) gc, &xgval, GDK_GC_JOIN_STYLE);
   }
   if (mask & kGCSubwindowMode) {
      gdk_gc_set_subwindow((GdkGC *) gc, xgval.subwindow_mode);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copies graphics context from org to dest. Only the values specified
/// in mask are copied. Both org and dest must exist.

void TGWin32::CopyGC(GContext_t org, GContext_t dest, Mask_t mask)
{
   GCValues_t gval;
   GdkGCValues xgval;
   ULong_t xmask;

   if (!mask) {
      // in this case copy all fields
      mask = (Mask_t) - 1;
   }

   gval.fMask = mask;           // only set fMask used to convert to xmask
   MapGCValues(gval, xmask, xgval, kTRUE);

   gdk_gc_copy((GdkGC *) dest, (GdkGC *) org);
}

////////////////////////////////////////////////////////////////////////////////
/// Explicitely delete a graphics context.

void TGWin32::DeleteGC(GContext_t gc)
{
   gdk_gc_unref((GdkGC *) gc);
}

////////////////////////////////////////////////////////////////////////////////
/// Create cursor handle (just return cursor from cursor pool fCursors).

Cursor_t TGWin32::CreateCursor(ECursor cursor)
{
   return (Cursor_t) fCursors[cursor];
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a pixmap of the width and height you specified
/// and returns a pixmap ID that identifies it.

Pixmap_t TGWin32::CreatePixmap(Drawable_t id, UInt_t w, UInt_t h)
{
   GdkWindow *wid = (GdkWindow *)id;
   if (!id) wid =  GDK_ROOT_PARENT();

   return (Pixmap_t) gdk_pixmap_new(wid, w, h, gdk_visual_get_best_depth());
}

////////////////////////////////////////////////////////////////////////////////
/// Create a pixmap from bitmap data. Ones will get foreground color and
/// zeroes background color.

Pixmap_t TGWin32::CreatePixmap(Drawable_t id, const char *bitmap,
                               UInt_t width, UInt_t height,
                               ULong_t forecolor, ULong_t backcolor,
                               Int_t depth)
{
   GdkColor fore, back;
   fore.pixel = forecolor;
   fore.red = GetRValue(forecolor);
   fore.green = GetGValue(forecolor);
   fore.blue = GetBValue(forecolor);

   back.pixel = backcolor;
   back.red = GetRValue(backcolor);
   back.green = GetGValue(backcolor);
   back.blue = GetBValue(backcolor);

   GdkWindow *wid = (GdkWindow *)id;
   if (!id) wid =  GDK_ROOT_PARENT();

   return (Pixmap_t) gdk_pixmap_create_from_data(wid, (char *) bitmap, width,
                                                 height, depth, &fore, &back);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a bitmap (i.e. pixmap with depth 1) from the bitmap data.

Pixmap_t TGWin32::CreateBitmap(Drawable_t id, const char *bitmap,
                               UInt_t width, UInt_t height)
{
   GdkWindow *wid = (GdkWindow *)id;
   if (!id) wid =  GDK_ROOT_PARENT();

   Pixmap_t ret = (Pixmap_t) gdk_bitmap_create_from_data(wid,
                                                 (char *)bitmap, width, height);
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Explicitely delete pixmap resource.

void TGWin32::DeletePixmap(Pixmap_t pmap)
{
   gdk_pixmap_unref((GdkPixmap *) pmap);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a picture pixmap from data on file. The picture attributes
/// are used for input and output. Returns kTRUE in case of success,
/// kFALSE otherwise. If mask does not exist it is set to kNone.

Bool_t TGWin32::CreatePictureFromFile(Drawable_t id, const char *filename,
                                      Pixmap_t & pict,
                                      Pixmap_t & pict_mask,
                                      PictureAttributes_t & attr)
{
   GdkBitmap *gdk_pixmap_mask;
   if (strstr(filename, ".xpm") || strstr(filename, ".XPM")) {
      GdkWindow *wid = (GdkWindow *)id;
      if (!id) wid =  GDK_ROOT_PARENT();

      pict = (Pixmap_t) gdk_pixmap_create_from_xpm(wid, &gdk_pixmap_mask, 0,
                                                filename);
      pict_mask = (Pixmap_t) gdk_pixmap_mask;
   } else if (strstr(filename, ".gif") || strstr(filename, ".GIF")) {
      pict = ReadGIF(0, 0, filename, id);
      pict_mask = kNone;
   }

   gdk_drawable_get_size((GdkPixmap *) pict, (int *) &attr.fWidth,
                         (int *) &attr.fHeight);
   if (pict) {
      return kTRUE;
   }
   if (pict_mask) {
      pict_mask = kNone;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a pixture pixmap from data. The picture attributes
/// are used for input and output. Returns kTRUE in case of success,
/// kFALSE otherwise. If mask does not exist it is set to kNone.

Bool_t TGWin32::CreatePictureFromData(Drawable_t id, char **data,
                                      Pixmap_t & pict,
                                      Pixmap_t & pict_mask,
                                      PictureAttributes_t & attr)
{
   GdkBitmap *gdk_pixmap_mask;
   GdkWindow *wid = (GdkWindow *)id;
   if (!id) wid =  GDK_ROOT_PARENT();

   pict = (Pixmap_t) gdk_pixmap_create_from_xpm_d(wid, &gdk_pixmap_mask, 0,
                                                  data);
   pict_mask = (Pixmap_t) gdk_pixmap_mask;

   if (pict) {
      return kTRUE;
   }
   if (pict_mask) {
      pict_mask = kNone;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Read picture data from file and store in ret_data. Returns kTRUE in
/// case of success, kFALSE otherwise.

Bool_t TGWin32::ReadPictureDataFromFile(const char *filename, char ***ret_data)
{
   Bool_t ret = kFALSE;
   GdkPixmap *pxm = gdk_pixmap_create_from_xpm(NULL, NULL, NULL, filename);
   ret_data = 0;

   if (pxm==NULL) return kFALSE;

   HBITMAP hbm = (HBITMAP)GDK_DRAWABLE_XID(pxm);
   BITMAP bitmap;

   ret = ::GetObject(hbm, sizeof(HBITMAP), (LPVOID)&bitmap);
   ret_data = (char ***)&bitmap.bmBits;
   gdk_pixmap_unref(pxm);
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete picture data created by the function ReadPictureDataFromFile.

void TGWin32::DeletePictureData(void *data)
{
   free(data);
}

////////////////////////////////////////////////////////////////////////////////
/// Specify a dash pattertn. Offset defines the phase of the pattern.
/// Each element in the dash_list array specifies the length (in pixels)
/// of a segment of the pattern. N defines the length of the list.

void TGWin32::SetDashes(GContext_t gc, Int_t offset, const char *dash_list,
                        Int_t n)
{
   int i;
   gint8 dashes[32];
   for (i = 0; i < n; i++) {
      dashes[i] = (gint8) dash_list[i];
   }
   for (i = n; i < 32; i++) {
      dashes[i] = (gint8) 0;
   }

   gdk_gc_set_dashes((GdkGC *) gc, offset, dashes, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Map a ColorStruct_t to a XColor structure.

void TGWin32::MapColorStruct(ColorStruct_t * color, GdkColor & xcolor)
{
   xcolor.pixel = color->fPixel;
   xcolor.red = color->fRed;
   xcolor.green = color->fGreen;
   xcolor.blue = color->fBlue;
}

////////////////////////////////////////////////////////////////////////////////
/// Parse string cname containing color name, like "green" or "#00FF00".
/// It returns a filled in ColorStruct_t. Returns kFALSE in case parsing
/// failed, kTRUE in case of success. On success, the ColorStruct_t
/// fRed, fGreen and fBlue fields are all filled in and the mask is set
/// for all three colors, but fPixel is not set.

Bool_t TGWin32::ParseColor(Colormap_t cmap, const char *cname,
                           ColorStruct_t & color)
{
   GdkColor xc;

   if (gdk_color_parse((char *)cname, &xc)) {
      color.fPixel = xc.pixel = RGB(xc.red, xc.green, xc.blue);
      color.fRed = xc.red;
      color.fGreen = xc.green;
      color.fBlue = xc.blue;
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Find and allocate a color cell according to the color values specified
/// in the ColorStruct_t. If no cell could be allocated it returns kFALSE,
/// otherwise kTRUE.

Bool_t TGWin32::AllocColor(Colormap_t cmap, ColorStruct_t & color)
{
   int status;
   GdkColor xc;

   xc.red = color.fRed;
   xc.green = color.fGreen;
   xc.blue = color.fBlue;

   status = gdk_colormap_alloc_color((GdkColormap *) cmap, &xc, FALSE, TRUE);
   color.fPixel = xc.pixel;

   return kTRUE;                // status != 0 ? kTRUE : kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill in the primary color components for a specific pixel value.
/// On input fPixel should be set on return the fRed, fGreen and
/// fBlue components will be set.

void TGWin32::QueryColor(Colormap_t cmap, ColorStruct_t & color)
{
   GdkColor xc;
   xc.pixel = color.fPixel;

   GdkColorContext *cc = gdk_color_context_new(gdk_visual_get_system(), fColormap);
   gdk_color_context_query_color(cc, &xc);
   gdk_color_context_free(cc);

   color.fPixel = xc.pixel;
   color.fRed = xc.red;
   color.fGreen = xc.green;
   color.fBlue = xc.blue;
}

////////////////////////////////////////////////////////////////////////////////
/// Free color cell with specified pixel value.

void TGWin32::FreeColor(Colormap_t cmap, ULong_t pixel)
{
   // FIXME: to be implemented.
}

////////////////////////////////////////////////////////////////////////////////
/// Check if there is for window "id" an event of type "type". If there
/// is fill in the event structure and return true. If no such event
/// return false.

Bool_t TGWin32::CheckEvent(Window_t id, EGEventType type, Event_t & ev)
{
   if (!id) return kFALSE;

   Event_t tev;
   GdkEvent xev;

   tev.fType = type;
   tev.fWindow = (Window_t) id;
   tev.fTime = 0;
   tev.fX = tev.fY = 0;
   tev.fXRoot = tev.fYRoot = 0;
   tev.fCode = 0;
   tev.fState = 0;
   tev.fWidth = tev.fHeight = 0;
   tev.fCount = 0;
   tev.fSendEvent = kFALSE;
   tev.fHandle = 0;
   tev.fFormat = 0;
   tev.fUser[0] = tev.fUser[1] = tev.fUser[2] = tev.fUser[3] = tev.fUser[4] = 0L;

   TGWin32MainThread::LockMSG();
   MapEvent(tev, xev, kTRUE);
   Bool_t r = gdk_check_typed_window_event((GdkWindow *) id, xev.type, &xev);

   if (r) MapEvent(ev, xev, kFALSE);
   TGWin32MainThread::UnlockMSG();

   return r ? kTRUE : kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Send event ev to window id.

void TGWin32::SendEvent(Window_t id, Event_t * ev)
{
   if (!ev || !id) return;

   TGWin32MainThread::LockMSG();
   GdkEvent xev;
   MapEvent(*ev, xev, kTRUE);
   gdk_event_put(&xev);
   TGWin32MainThread::UnlockMSG();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns number of pending events.

Int_t TGWin32::EventsPending()
{
   Int_t ret;

   TGWin32MainThread::LockMSG();
   ret = (Int_t)gdk_event_queue_find_first();
   TGWin32MainThread::UnlockMSG();

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Copies first pending event from event queue to Event_t structure
/// and removes event from queue. Not all of the event fields are valid
/// for each event type, except fType and fWindow.

void TGWin32::NextEvent(Event_t & event)
{
   TGWin32MainThread::LockMSG();
   GdkEvent *xev = gdk_event_unqueue();

   // fill in Event_t
   event.fType = kOtherEvent;   // bb add
   if (xev == NULL) {
      TGWin32MainThread::UnlockMSG();
      return;
   }
   MapEvent(event, *xev, kFALSE);
   gdk_event_free (xev);
   TGWin32MainThread::UnlockMSG();
}

////////////////////////////////////////////////////////////////////////////////
/// Map modifier key state to or from X.

void TGWin32::MapModifierState(UInt_t & state, UInt_t & xstate, Bool_t tox)
{
   if (tox) {
      xstate = state;
      if (state & kAnyModifier) {
         xstate = GDK_MODIFIER_MASK;
      }
   } else {
      state = xstate;
   }
}

static void _set_event_time(GdkEvent &event, UInt_t time)
{
   // set gdk event time

   switch (event.type) {
      case GDK_MOTION_NOTIFY:
         event.motion.time = time;
      case GDK_BUTTON_PRESS:
      case GDK_2BUTTON_PRESS:
      case GDK_3BUTTON_PRESS:
      case GDK_BUTTON_RELEASE:
      case GDK_SCROLL:
         event.button.time = time;
      case GDK_KEY_PRESS:
      case GDK_KEY_RELEASE:
         event.key.time = time;
      case GDK_ENTER_NOTIFY:
      case GDK_LEAVE_NOTIFY:
         event.crossing.time = time;
      case GDK_PROPERTY_NOTIFY:
         event.property.time = time;
      case GDK_SELECTION_CLEAR:
      case GDK_SELECTION_REQUEST:
      case GDK_SELECTION_NOTIFY:
         event.selection.time = time;
      case GDK_PROXIMITY_IN:
      case GDK_PROXIMITY_OUT:
         event.proximity.time = time;
      case GDK_DRAG_ENTER:
      case GDK_DRAG_LEAVE:
      case GDK_DRAG_MOTION:
      case GDK_DRAG_STATUS:
      case GDK_DROP_START:
      case GDK_DROP_FINISHED:
         event.dnd.time = time;
      default:                 /* use current time */
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Map Event_t structure to gdk_event structure. If tox is false
/// map the other way.

void TGWin32::MapEvent(Event_t & ev, GdkEvent & xev, Bool_t tox)
{
   if (tox) {
      // map from Event_t to gdk_event
      xev.type = GDK_NOTHING;
      if (ev.fType == kGKeyPress)
         xev.type = GDK_KEY_PRESS;
      if (ev.fType == kKeyRelease)
         xev.type = GDK_KEY_RELEASE;
      if (ev.fType == kButtonPress)
         xev.type = GDK_BUTTON_PRESS;
      if (ev.fType == kButtonRelease)
         xev.type = GDK_BUTTON_RELEASE;
      if (ev.fType == kMotionNotify)
         xev.type = GDK_MOTION_NOTIFY;
      if (ev.fType == kEnterNotify)
         xev.type = GDK_ENTER_NOTIFY;
      if (ev.fType == kLeaveNotify)
         xev.type = GDK_LEAVE_NOTIFY;
      if (ev.fType == kExpose)
         xev.type = GDK_EXPOSE;
      if (ev.fType == kConfigureNotify)
         xev.type = GDK_CONFIGURE;
      if (ev.fType == kMapNotify)
         xev.type = GDK_MAP;
      if (ev.fType == kUnmapNotify)
         xev.type = GDK_UNMAP;
      if (ev.fType == kDestroyNotify)
         xev.type = GDK_DESTROY;
      if (ev.fType == kClientMessage)
         xev.type = GDK_CLIENT_EVENT;
      if (ev.fType == kSelectionClear)
         xev.type = GDK_SELECTION_CLEAR;
      if (ev.fType == kSelectionRequest)
         xev.type = GDK_SELECTION_REQUEST;
      if (ev.fType == kSelectionNotify)
         xev.type = GDK_SELECTION_NOTIFY;

      xev.any.type = xev.type;
      xev.any.send_event = ev.fSendEvent;
      if (ev.fType == kDestroyNotify) {
         xev.any.window = (GdkWindow *) ev.fWindow;
      }
      if (ev.fType == kFocusIn) {
         xev.type = GDK_FOCUS_CHANGE;
         xev.focus_change.type = xev.type;
         xev.focus_change.window = (GdkWindow *) ev.fWindow;
         xev.focus_change.in = TRUE;
      }
      if (ev.fType == kFocusOut) {
         xev.type = GDK_FOCUS_CHANGE;
         xev.focus_change.type = xev.type;
         xev.focus_change.window = (GdkWindow *) ev.fWindow;
         xev.focus_change.in = FALSE;
      }
      if (ev.fType == kGKeyPress || ev.fType == kKeyRelease) {
         xev.key.window = (GdkWindow *) ev.fWindow;
         xev.key.type = xev.type;
         MapModifierState(ev.fState, xev.key.state, kTRUE); // key mask
         xev.key.keyval = ev.fCode; // key code
      }
      if (ev.fType == kButtonPress || ev.fType == kButtonRelease) {
         xev.button.window = (GdkWindow *) ev.fWindow;
         xev.button.type = xev.type;
         xev.button.x = ev.fX;
         xev.button.y = ev.fY;
         xev.button.x_root = ev.fXRoot;
         xev.button.y_root = ev.fYRoot;
         MapModifierState(ev.fState, xev.button.state, kTRUE); // button mask
         xev.button.button = ev.fCode; // button code
      }
      if (ev.fType == kSelectionNotify) {
         xev.selection.window = (GdkWindow *) ev.fUser[0];
         xev.selection.requestor = (guint32) ev.fUser[0];
         xev.selection.selection = (GdkAtom) ev.fUser[1];
         xev.selection.target = (GdkAtom) ev.fUser[2];
         xev.selection.property = (GdkAtom) ev.fUser[3];
         xev.selection.type = xev.type;
      }
      if (ev.fType == kClientMessage) {
         if ((ev.fFormat == 32) && (ev.fHandle == gWM_DELETE_WINDOW)) {
            xev.type = GDK_DELETE;
            xev.any.type = xev.type;
            xev.any.window = (GdkWindow *) ev.fWindow;
         } else {
            xev.client.window = (GdkWindow *) ev.fWindow;
            xev.client.type = xev.type;
            xev.client.message_type = (GdkAtom) ev.fHandle;
            xev.client.data_format = ev.fFormat;
            xev.client.data.l[0] = ev.fUser[0];
            if (sizeof(ev.fUser[0]) > 4) {
               SplitLong(ev.fUser[1], xev.client.data.l[1],
                         xev.client.data.l[3]);
               SplitLong(ev.fUser[2], xev.client.data.l[2],
                         xev.client.data.l[4]);
            } else {
               xev.client.data.l[1] = ev.fUser[1];
               xev.client.data.l[2] = ev.fUser[2];
               xev.client.data.l[3] = ev.fUser[3];
               xev.client.data.l[4] = ev.fUser[4];
            }
         }
      }
      if (ev.fType == kMotionNotify) {
         xev.motion.window = (GdkWindow *) ev.fWindow;
         xev.motion.type = xev.type;
         xev.motion.x = ev.fX;
         xev.motion.y = ev.fY;
         xev.motion.x_root = ev.fXRoot;
         xev.motion.y_root = ev.fYRoot;
      }
      if ((ev.fType == kEnterNotify) || (ev.fType == kLeaveNotify)) {
         xev.crossing.window = (GdkWindow *) ev.fWindow;
         xev.crossing.type = xev.type;
         xev.crossing.x = ev.fX;
         xev.crossing.y = ev.fY;
         xev.crossing.x_root = ev.fXRoot;
         xev.crossing.y_root = ev.fYRoot;
         xev.crossing.mode = (GdkCrossingMode) ev.fCode; // NotifyNormal, NotifyGrab, NotifyUngrab
         MapModifierState(ev.fState, xev.crossing.state, kTRUE);  // key or button mask
      }
      if (ev.fType == kExpose) {
         xev.expose.window = (GdkWindow *) ev.fWindow;
         xev.expose.type = xev.type;
         xev.expose.area.x = ev.fX;
         xev.expose.area.y = ev.fY;
         xev.expose.area.width = ev.fWidth;  // width and
         xev.expose.area.height = ev.fHeight;   // height of exposed area
         xev.expose.count = ev.fCount; // number of expose events still to come
      }
      if (ev.fType == kConfigureNotify) {
         xev.configure.window = (GdkWindow *) ev.fWindow;
         xev.configure.type = xev.type;
         xev.configure.x = ev.fX;
         xev.configure.y = ev.fY;
         xev.configure.width = ev.fWidth;
         xev.configure.height = ev.fHeight;
      }
      if (ev.fType == kSelectionClear) {
         xev.selection.window = (GdkWindow *) ev.fWindow;
         xev.selection.type = xev.type;
         xev.selection.selection = ev.fUser[0];
      }
      if (ev.fType == kSelectionRequest) {
         xev.selection.window = (GdkWindow *) ev.fUser[0];
         xev.selection.type = xev.type;
         xev.selection.selection = ev.fUser[1];
         xev.selection.target = ev.fUser[2];
         xev.selection.property = ev.fUser[3];
      }
      if ((ev.fType == kMapNotify) || (ev.fType == kUnmapNotify)) {
         xev.any.window = (GdkWindow *) ev.fWindow;
      }
      if (xev.type != GDK_CLIENT_EVENT)
         _set_event_time(xev, ev.fTime);
   } else {
      // map from gdk_event to Event_t
      ev.fType = kOtherEvent;
      if (xev.type == GDK_KEY_PRESS)
         ev.fType = kGKeyPress;
      if (xev.type == GDK_KEY_RELEASE)
         ev.fType = kKeyRelease;
      if (xev.type == GDK_BUTTON_PRESS)
         ev.fType = kButtonPress;
      if (xev.type == GDK_BUTTON_RELEASE)
         ev.fType = kButtonRelease;
      if (xev.type == GDK_MOTION_NOTIFY)
         ev.fType = kMotionNotify;
      if (xev.type == GDK_ENTER_NOTIFY)
         ev.fType = kEnterNotify;
      if (xev.type == GDK_LEAVE_NOTIFY)
         ev.fType = kLeaveNotify;
      if (xev.type == GDK_EXPOSE)
         ev.fType = kExpose;
      if (xev.type == GDK_CONFIGURE)
         ev.fType = kConfigureNotify;
      if (xev.type == GDK_MAP)
         ev.fType = kMapNotify;
      if (xev.type == GDK_UNMAP)
         ev.fType = kUnmapNotify;
      if (xev.type == GDK_DESTROY)
         ev.fType = kDestroyNotify;
      if (xev.type == GDK_SELECTION_CLEAR)
         ev.fType = kSelectionClear;
      if (xev.type == GDK_SELECTION_REQUEST)
         ev.fType = kSelectionRequest;
      if (xev.type == GDK_SELECTION_NOTIFY)
         ev.fType = kSelectionNotify;

      ev.fSendEvent = kFALSE; //xev.any.send_event ? kTRUE : kFALSE;
      ev.fTime = gdk_event_get_time((GdkEvent *)&xev);
      ev.fWindow = (Window_t) xev.any.window;

      if ((xev.type == GDK_MAP) || (xev.type == GDK_UNMAP)) {
         ev.fWindow = (Window_t) xev.any.window;
      }
      if (xev.type == GDK_DELETE) {
         ev.fWindow = (Window_t) xev.any.window;
         ev.fType = kClientMessage;
         ev.fFormat = 32;
         ev.fHandle = gWM_DELETE_WINDOW;
         ev.fUser[0] = (Long_t) gWM_DELETE_WINDOW;
         if (sizeof(ev.fUser[0]) > 4) {
            AsmLong(xev.client.data.l[1], xev.client.data.l[3],
                    ev.fUser[1]);
            AsmLong(xev.client.data.l[2], xev.client.data.l[4],
                    ev.fUser[2]);
         } else {
            ev.fUser[1] = 0; //xev.client.data.l[1];
            ev.fUser[2] = 0; //xev.client.data.l[2];
            ev.fUser[3] = 0; //xev.client.data.l[3];
            ev.fUser[4] = 0; //xev.client.data.l[4];
         }
      }
      if (xev.type == GDK_DESTROY) {
         ev.fType = kDestroyNotify;
         ev.fHandle = (Window_t) xev.any.window;   // window to be destroyed
         ev.fWindow = (Window_t) xev.any.window;
      }
      if (xev.type == GDK_FOCUS_CHANGE) {
         ev.fWindow = (Window_t) xev.focus_change.window;
         ev.fCode = kNotifyNormal;
         ev.fState = 0;
         if (xev.focus_change.in == TRUE) {
            ev.fType = kFocusIn;
         } else {
            ev.fType = kFocusOut;
         }
      }
      if (ev.fType == kGKeyPress || ev.fType == kKeyRelease) {
         ev.fWindow = (Window_t) xev.key.window;
         MapModifierState(ev.fState, xev.key.state, kFALSE);   // key mask
         ev.fCode = xev.key.keyval;  // key code
         ev.fUser[1] = xev.key.length;
         if (xev.key.length > 0) ev.fUser[2] = xev.key.string[0];
         if (xev.key.length > 1) ev.fUser[3] = xev.key.string[1];
         if (xev.key.length > 2) ev.fUser[4] = xev.key.string[2];
         HWND tmpwin = (HWND) GetWindow((HWND) GDK_DRAWABLE_XID((GdkWindow *)xev.key.window), GW_CHILD);
         if (tmpwin) {
            ev.fUser[0] = (ULong_t) gdk_xid_table_lookup((HANDLE)tmpwin);
         } else {
            ev.fUser[0] = (ULong_t) xev.key.window;
         }
      }
      if (ev.fType == kButtonPress || ev.fType == kButtonRelease) {
         ev.fWindow = (Window_t) xev.button.window;
         ev.fX = xev.button.x;
         ev.fY = xev.button.y;
         ev.fXRoot = xev.button.x_root;
         ev.fYRoot = xev.button.y_root;
         MapModifierState(ev.fState, xev.button.state, kFALSE);   // button mask
         ev.fCode = xev.button.button; // button code
         POINT tpoint;
         tpoint.x = xev.button.x;
         tpoint.y = xev.button.y;
         HWND tmpwin = ChildWindowFromPoint((HWND) GDK_DRAWABLE_XID((GdkWindow *)xev.button.window), tpoint);
         if (tmpwin) {
             ev.fUser[0] = (ULong_t) gdk_xid_table_lookup((HANDLE)tmpwin);
         } else {
            ev.fUser[0] = (ULong_t) 0;
         }
      }
      if (ev.fType == kMotionNotify) {
         ev.fWindow = (Window_t) xev.motion.window;
         ev.fX = xev.motion.x;
         ev.fY = xev.motion.y;
         ev.fXRoot = xev.motion.x_root;
         ev.fYRoot = xev.motion.y_root;
         MapModifierState(ev.fState, xev.motion.state, kFALSE);   // key or button mask

         POINT tpoint;
         tpoint.x = xev.button.x;
         tpoint.y = xev.button.y;
         HWND tmpwin = ChildWindowFromPoint((HWND) GDK_DRAWABLE_XID((GdkWindow *)xev.motion.window), tpoint);
         if (tmpwin) {
             ev.fUser[0] = (ULong_t)gdk_xid_table_lookup((HANDLE)tmpwin);
         } else {
            ev.fUser[0] = (ULong_t) xev.motion.window;
         }
      }
      if (ev.fType == kEnterNotify || ev.fType == kLeaveNotify) {
         ev.fWindow = (Window_t) xev.crossing.window;
         ev.fX = xev.crossing.x;
         ev.fY = xev.crossing.y;
         ev.fXRoot = xev.crossing.x_root;
         ev.fYRoot = xev.crossing.y_root;
         ev.fCode = xev.crossing.mode; // NotifyNormal, NotifyGrab, NotifyUngrab
         MapModifierState(ev.fState, xev.crossing.state, kFALSE); // key or button mask
      }
      if (ev.fType == kExpose) {
         ev.fWindow = (Window_t) xev.expose.window;
         ev.fX = xev.expose.area.x;
         ev.fY = xev.expose.area.y;
         ev.fWidth = xev.expose.area.width;  // width and
         ev.fHeight = xev.expose.area.height;   // height of exposed area
         ev.fCount = xev.expose.count; // number of expose events still to come
      }
      if (ev.fType == kConfigureNotify) {
         ev.fWindow = (Window_t) xev.configure.window;
         ev.fX = xev.configure.x;
         ev.fY = xev.configure.y;
         ev.fWidth = xev.configure.width;
         ev.fHeight = xev.configure.height;
      }
      if (xev.type == GDK_CLIENT_EVENT) {
         ev.fWindow = (Window_t) xev.client.window;
         ev.fType = kClientMessage;
         ev.fHandle = xev.client.message_type;
         ev.fFormat = xev.client.data_format;
         ev.fUser[0] = xev.client.data.l[0];
         if (sizeof(ev.fUser[0]) > 4) {
            AsmLong(xev.client.data.l[1], xev.client.data.l[3],
                    ev.fUser[1]);
            AsmLong(xev.client.data.l[2], xev.client.data.l[4],
                    ev.fUser[2]);
         } else {
            ev.fUser[1] = xev.client.data.l[1];
            ev.fUser[2] = xev.client.data.l[2];
            ev.fUser[3] = xev.client.data.l[3];
            ev.fUser[4] = xev.client.data.l[4];
         }
      }
      if (ev.fType == kSelectionClear) {
         ev.fWindow = (Window_t) xev.selection.window;
         ev.fUser[0] = xev.selection.selection;
      }
      if (ev.fType == kSelectionRequest) {
         ev.fWindow = (Window_t) xev.selection.window;
         ev.fUser[0] = (ULong_t) xev.selection.window;
         ev.fUser[1] = xev.selection.selection;
         ev.fUser[2] = xev.selection.target;
         ev.fUser[3] = xev.selection.property;
      }
      if (ev.fType == kSelectionNotify) {
         ev.fWindow = (Window_t) xev.selection.window;
         ev.fUser[0] = (ULong_t) xev.selection.window;
         ev.fUser[1] = xev.selection.selection;
         ev.fUser[2] = xev.selection.target;
         ev.fUser[3] = xev.selection.property;
      }
      if (xev.type == GDK_SCROLL) {
         ev.fType = kButtonRelease;
         if (xev.scroll.direction == GDK_SCROLL_UP) {
            ev.fCode = kButton4;
         } else if (xev.scroll.direction == GDK_SCROLL_DOWN) {
            ev.fCode = kButton5;
         }
         ev.fWindow = (Window_t) xev.scroll.window;
         ev.fX = xev.scroll.x;
         ev.fY = xev.scroll.y;
         ev.fXRoot = xev.scroll.x_root;
         ev.fYRoot = xev.scroll.y_root;
         POINT tpoint;
         tpoint.x = xev.scroll.x;
         tpoint.y = xev.scroll.y;
         HWND tmpwin = ChildWindowFromPoint((HWND) GDK_DRAWABLE_XID((GdkWindow *)xev.scroll.window), tpoint);
         if (tmpwin) {
             ev.fUser[0] = (ULong_t)gdk_xid_table_lookup((HANDLE)tmpwin);
         } else {
            ev.fUser[0] = (ULong_t) 0;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::Bell(Int_t percent)
{
   gSystem->Beep();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a drawable (i.e. pixmap) to another drawable (pixmap, window).
/// The graphics context gc will be used and the source will be copied
/// from src_x,src_y,src_x+width,src_y+height to dest_x,dest_y.

void TGWin32::CopyArea(Drawable_t src, Drawable_t dest, GContext_t gc,
                       Int_t src_x, Int_t src_y, UInt_t width,
                       UInt_t height, Int_t dest_x, Int_t dest_y)
{
   if (!src || !dest) return;

   gdk_window_copy_area((GdkDrawable *) dest, (GdkGC *) gc, dest_x, dest_y,
                        (GdkDrawable *) src, src_x, src_y, width, height);
}

////////////////////////////////////////////////////////////////////////////////
/// Change window attributes.

void TGWin32::ChangeWindowAttributes(Window_t id, SetWindowAttributes_t * attr)
{
   if (!id) return;

   GdkColor color;
   UInt_t xevmask;
   Mask_t evmask;

   if (attr && (attr->fMask & kWAEventMask)) {
      evmask = (Mask_t) attr->fEventMask;
      MapEventMask(evmask, xevmask);
      gdk_window_set_events((GdkWindow *) id, (GdkEventMask) xevmask);
   }
   if (attr && (attr->fMask & kWABackPixel)) {
      color.pixel = attr->fBackgroundPixel;
      color.red = GetRValue(attr->fBackgroundPixel);
      color.green = GetGValue(attr->fBackgroundPixel);
      color.blue = GetBValue(attr->fBackgroundPixel);
      gdk_window_set_background((GdkWindow *) id, &color);
   }
//   if (attr && (attr->fMask & kWAOverrideRedirect))
//      gdk_window_set_override_redirect ((GdkWindow *) id, attr->fOverrideRedirect);
   if (attr && (attr->fMask & kWABackPixmap)) {
      gdk_window_set_back_pixmap((GdkWindow *) id,
                                 (GdkPixmap *) attr->fBackgroundPixmap, 0);
   }
   if (attr && (attr->fMask & kWACursor)) {
      gdk_window_set_cursor((GdkWindow *) id, (GdkCursor *) attr->fCursor);
   }
   if (attr && (attr->fMask & kWAColormap)) {
      gdk_window_set_colormap((GdkWindow *) id,(GdkColormap *) attr->fColormap);
   }
   if (attr && (attr->fMask & kWABorderWidth)) {
      if (attr->fBorderWidth > 0) {
         gdk_window_set_decorations((GdkWindow *) id,
                                    (GdkWMDecoration) GDK_DECOR_BORDER);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// This function alters the property for the specified window and
/// causes the X server to generate a PropertyNotify event on that
/// window.

void TGWin32::ChangeProperty(Window_t id, Atom_t property, Atom_t type,
                             UChar_t * data, Int_t len)
{
   if (!id) return;

   gdk_property_change((GdkWindow *) id, (GdkAtom) property,
                       (GdkAtom) type, 8, GDK_PROP_MODE_REPLACE, data,len);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a line.

void TGWin32::DrawLine(Drawable_t id, GContext_t gc, Int_t x1, Int_t y1,
                       Int_t x2, Int_t y2)
{
   if (!id) return;

   gdk_draw_line((GdkDrawable *) id, (GdkGC *) gc, x1, y1, x2, y2);
}

////////////////////////////////////////////////////////////////////////////////
/// Clear a window area to the bakcground color.

void TGWin32::ClearArea(Window_t id, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   if (!id) return;

   gdk_window_clear_area((GdkWindow *) id, x, y, w, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Tell WM to send message when window is closed via WM.

void TGWin32::WMDeleteNotify(Window_t id)
{
   if (!id) return;

   Atom prop;
   prop = (Atom_t) gdk_atom_intern("WM_DELETE_WINDOW", FALSE);

   W32ChangeProperty((HWND) GDK_DRAWABLE_XID((GdkWindow *) id),
                     prop, XA_ATOM, 32, GDK_PROP_MODE_REPLACE,
                     (unsigned char *) &gWM_DELETE_WINDOW, 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Turn key auto repeat on or off.

void TGWin32::SetKeyAutoRepeat(Bool_t on)
{
   if (on) {
      gdk_key_repeat_restore();
    } else {
      gdk_key_repeat_disable();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Establish passive grab on a certain key. That is, when a certain key
/// keycode is hit while certain modifier's (Shift, Control, Meta, Alt)
/// are active then the keyboard will be grabed for window id.
/// When grab is false, ungrab the keyboard for this key and modifier.

void TGWin32::GrabKey(Window_t id, Int_t keycode, UInt_t modifier, Bool_t grab)
{
   UInt_t xmod;

   MapModifierState(modifier, xmod);

   if (grab) {
      gdk_key_grab(keycode, (GdkEventMask)xmod, (GdkWindow *)id);
   } else {
      gdk_key_ungrab(keycode, (GdkEventMask)xmod, (GdkWindow *)id);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Establish passive grab on a certain mouse button. That is, when a
/// certain mouse button is hit while certain modifier's (Shift, Control,
/// Meta, Alt) are active then the mouse will be grabed for window id.
/// When grab is false, ungrab the mouse button for this button and modifier.

void TGWin32::GrabButton(Window_t id, EMouseButton button, UInt_t modifier,
                         UInt_t evmask, Window_t confine, Cursor_t cursor,
                         Bool_t grab)
{
   UInt_t xevmask;
   UInt_t xmod;

   if (!id) return;

   MapModifierState(modifier, xmod);

   if (grab) {
      MapEventMask(evmask, xevmask);
      gdk_button_grab(button, xmod, ( GdkWindow *)id, 1,  (GdkEventMask)xevmask,
                      (GdkWindow*)confine,  (GdkCursor*)cursor);
   } else {
      gdk_button_ungrab(button, xmod, ( GdkWindow *)id);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Establish an active pointer grab. While an active pointer grab is in
/// effect, further pointer events are only reported to the grabbing
/// client window.

void TGWin32::GrabPointer(Window_t id, UInt_t evmask, Window_t confine,
                          Cursor_t cursor, Bool_t grab, Bool_t owner_events)
{
   UInt_t xevmask;
   MapEventMask(evmask, xevmask);

   if (grab) {
      if(!::IsWindowVisible((HWND)GDK_DRAWABLE_XID(id))) return;
      gdk_pointer_grab((GdkWindow *) id, owner_events, (GdkEventMask) xevmask,
                       (GdkWindow *) confine, (GdkCursor *) cursor,
                       GDK_CURRENT_TIME);
   } else {
      gdk_pointer_ungrab(GDK_CURRENT_TIME);
      ::SetCursor((HCURSOR)GDK_CURSOR_XID(fCursors[kPointer]));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set window name.

void TGWin32::SetWindowName(Window_t id, char *name)
{
   if (!id) return;

   gdk_window_set_title((GdkWindow *) id, name);
}

////////////////////////////////////////////////////////////////////////////////
/// Set window icon name.

void TGWin32::SetIconName(Window_t id, char *name)
{
   if (!id) return;

   gdk_window_set_icon_name((GdkWindow *) id, name);
}

////////////////////////////////////////////////////////////////////////////////
/// Set pixmap the WM can use when the window is iconized.

void TGWin32::SetIconPixmap(Window_t id, Pixmap_t pic)
{
   if (!id) return;

   gdk_window_set_icon((GdkWindow *)id, NULL, (GdkPixmap *)pic, (GdkPixmap *)pic);
}

#define safestrlen(s) ((s) ? strlen(s) : 0)

////////////////////////////////////////////////////////////////////////////////
/// Set the windows class and resource name.

void TGWin32::SetClassHints(Window_t id, char *className, char *resourceName)
{
   if (!id) return;

   char *class_string;
   char *s;
   int len_nm, len_cl;
   GdkAtom prop;

   prop = gdk_atom_intern("WM_CLASS", kFALSE);

   len_nm = safestrlen(resourceName);
   len_cl = safestrlen(className);

   if ((class_string = s =
        (char *) malloc((unsigned) (len_nm + len_cl + 2)))) {
      if (len_nm) {
         strcpy(s, resourceName);
         s += len_nm + 1;
      } else
         *s++ = '\0';
      if (len_cl) {
         strcpy(s, className);
      } else {
         *s = '\0';
      }

      W32ChangeProperty((HWND) GDK_DRAWABLE_XID((GdkWindow *) id),
                            (Atom) XA_WM_CLASS, (Atom) XA_WM_CLASS, 8,
                            GDK_PROP_MODE_REPLACE,
                            (unsigned char *) class_string,
                            len_nm + len_cl + 2);
      free(class_string);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set decoration style for MWM-compatible wm (mwm, ncdwm, fvwm?).

void TGWin32::SetMWMHints(Window_t id, UInt_t value, UInt_t funcs,
                          UInt_t input)
{
   if (!id) return;

   gdk_window_set_decorations((GdkDrawable *) id, (GdkWMDecoration) value);
   gdk_window_set_functions((GdkDrawable *) id, (GdkWMFunction) funcs);
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::SetWMPosition(Window_t id, Int_t x, Int_t y)
{
   if (!id) return;

   gdk_window_move((GdkDrawable *) id, x, y);
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::SetWMSize(Window_t id, UInt_t w, UInt_t h)
{
   if (!id) return;

   gdk_window_resize((GdkWindow *) id, w, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Give the window manager minimum and maximum size hints. Also
/// specify via winc and hinc the resize increments.

void TGWin32::SetWMSizeHints(Window_t id, UInt_t wmin, UInt_t hmin,
                             UInt_t wmax, UInt_t hmax,
                             UInt_t winc, UInt_t hinc)
{
   if (!id) return;

   GdkGeometry hints;
   GdkWindowHints flags;

   flags = (GdkWindowHints) (GDK_HINT_MIN_SIZE | GDK_HINT_MAX_SIZE |
                             GDK_HINT_RESIZE_INC);
   hints.min_width = (Int_t) wmin;
   hints.max_width = (Int_t) wmax;
   hints.min_height = (Int_t) hmin;
   hints.max_height = (Int_t) hmax;
   hints.width_inc = (Int_t) winc;
   hints.height_inc = (Int_t) hinc;

   gdk_window_set_geometry_hints((GdkWindow *) id, (GdkGeometry *) &hints,
                                 (GdkWindowHints) flags);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the initial state of the window. Either kNormalState or kIconicState.

void TGWin32::SetWMState(Window_t id, EInitialState state)
{
   if (!id) return;

#if 0
   XWMHints hints;
   Int_t xstate = NormalState;

   if (state == kNormalState)
      xstate = NormalState;
   if (state == kIconicState)
      xstate = IconicState;

   hints.flags = StateHint;
   hints.initial_state = xstate;

   XSetWMHints((GdkWindow *) id, &hints);
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Tell window manager that window is a transient window of gdk_parent_root.

void TGWin32::SetWMTransientHint(Window_t id, Window_t main_id)
{
   if (!id) return;

   gdk_window_set_transient_for((GdkWindow *) id, (GdkWindow *) main_id);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a string using a specific graphics context in position (x,y).

void TGWin32::DrawString(Drawable_t id, GContext_t gc, Int_t x, Int_t y,
                         const char *s, Int_t len)
{
   if (!id) return;

   GdkGCValues values;
   gdk_gc_get_values((GdkGC *) gc, &values);
   gdk_win32_draw_text((GdkDrawable *) id, (GdkFont *) values.font,
                 (GdkGC *) gc, x, y, (const gchar *)s, len);
}

////////////////////////////////////////////////////////////////////////////////
/// Return length of string in pixels. Size depends on font.

Int_t TGWin32::TextWidth(FontStruct_t font, const char *s, Int_t len)
{
   return gdk_text_width((GdkFont *)font, s, len);
}

////////////////////////////////////////////////////////////////////////////////
/// Return some font properties.

void TGWin32::GetFontProperties(FontStruct_t font, Int_t & max_ascent,
                                Int_t & max_descent)
{
   GdkFont *f = (GdkFont *) font;
   max_ascent = f->ascent;
   max_descent = f->descent;
}

////////////////////////////////////////////////////////////////////////////////
/// Get current values from graphics context gc. Which values of the
/// context to get is encoded in the GCValues::fMask member.

void TGWin32::GetGCValues(GContext_t gc, GCValues_t & gval)
{
   GdkGCValues xgval;
   ULong_t xmask;

   MapGCValues(gval, xmask, xgval, kTRUE);
   gdk_gc_get_values((GdkGC *) gc, &xgval);
   MapGCValues(gval, xmask, xgval, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve associated font structure once we have the font handle.
/// Free returned FontStruct_t using FreeFontStruct().

FontStruct_t TGWin32::GetFontStruct(FontH_t fh)
{
   return (FontStruct_t) gdk_font_ref((GdkFont *) fh);
}

////////////////////////////////////////////////////////////////////////////////
/// Free font structure returned by GetFontStruct().

void TGWin32::FreeFontStruct(FontStruct_t fs)
{
   gdk_font_unref((GdkFont *) fs);
}

////////////////////////////////////////////////////////////////////////////////
/// Clear window.

void TGWin32::ClearWindow(Window_t id)
{
   if (!id) return;

   gdk_window_clear((GdkDrawable *) id);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert a keysym to the appropriate keycode. For example keysym is
/// a letter and keycode is the matching keyboard key (which is dependend
/// on the current keyboard mapping).

Int_t TGWin32::KeysymToKeycode(UInt_t keysym)
{
   UInt_t xkeysym;
   MapKeySym(keysym, xkeysym);
   return xkeysym;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a filled rectangle. Filling is done according to the gc.

void TGWin32::FillRectangle(Drawable_t id, GContext_t gc, Int_t x, Int_t y,
                            UInt_t w, UInt_t h)
{
   if (!id) return;

   gdk_win32_draw_rectangle((GdkDrawable *) id, (GdkGC *) gc, kTRUE, x, y, w, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a rectangle outline.

void TGWin32::DrawRectangle(Drawable_t id, GContext_t gc, Int_t x, Int_t y,
                            UInt_t w, UInt_t h)
{
   if (!id) return;

   gdk_win32_draw_rectangle((GdkDrawable *) id, (GdkGC *) gc, kFALSE, x, y, w, h);
}

////////////////////////////////////////////////////////////////////////////////
/// Draws multiple line segments. Each line is specified by a pair of points.

void TGWin32::DrawSegments(Drawable_t id, GContext_t gc, Segment_t * seg,
                           Int_t nseg)
{
   if (!id) return;

   gdk_win32_draw_segments((GdkDrawable *) id, (GdkGC *) gc, (GdkSegment *)seg, nseg);
}

////////////////////////////////////////////////////////////////////////////////
/// Defines which input events the window is interested in. By default
/// events are propageted up the window stack. This mask can also be
/// set at window creation time via the SetWindowAttributes_t::fEventMask
/// attribute.

void TGWin32::SelectInput(Window_t id, UInt_t evmask)
{
   if (!id) return;

   UInt_t xevmask;
   MapEventMask(evmask, xevmask, kTRUE);
   gdk_window_set_events((GdkWindow *) id, (GdkEventMask)xevmask);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the window id of the window having the input focus.

Window_t TGWin32::GetInputFocus()
{
   HWND hwnd = ::GetFocus();
   return (Window_t) gdk_xid_table_lookup(hwnd);
}

////////////////////////////////////////////////////////////////////////////////
/// Set keyboard input focus to window id.

void TGWin32::SetInputFocus(Window_t id)
{
   if (!id) return;

   HWND hwnd = (HWND)GDK_DRAWABLE_XID((GdkWindow *)id);
   ::SetFocus(hwnd);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the window id of the current owner of the primary selection.
/// That is the window in which, for example some text is selected.

Window_t TGWin32::GetPrimarySelectionOwner()
{
   return (Window_t)gdk_selection_owner_get(gClipboardAtom);
}

////////////////////////////////////////////////////////////////////////////////
/// Makes the window id the current owner of the primary selection.
/// That is the window in which, for example some text is selected.

void TGWin32::SetPrimarySelectionOwner(Window_t id)
{
   if (!id) return;

   gdk_selection_owner_set((GdkWindow *) id, gClipboardAtom, GDK_CURRENT_TIME, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// XConvertSelection() causes a SelectionRequest event to be sent to the
/// current primary selection owner. This event specifies the selection
/// property (primary selection), the format into which to convert that
/// data before storing it (target = XA_STRING), the property in which
/// the owner will place the information (sel_property), the window that
/// wants the information (id), and the time of the conversion request
/// (when).
/// The selection owner responds by sending a SelectionNotify event, which
/// confirms the selected atom and type.

void TGWin32::ConvertPrimarySelection(Window_t id, Atom_t clipboard, Time_t when)
{
   if (!id) return;

   gdk_selection_convert((GdkWindow *) id, clipboard,
                         gdk_atom_intern("GDK_TARGET_STRING", 0), when);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert the keycode from the event structure to a key symbol (according
/// to the modifiers specified in the event structure and the current
/// keyboard mapping). In buf a null terminated ASCII string is returned
/// representing the string that is currently mapped to the key code.

void TGWin32::LookupString(Event_t * event, char *buf, Int_t buflen,
                           UInt_t & keysym)
{
   _lookup_string(event, buf, buflen);
   UInt_t ks, xks = (UInt_t) event->fCode;
   MapKeySym(ks, xks, kFALSE);
   keysym = (Int_t) ks;
}

////////////////////////////////////////////////////////////////////////////////
/// Map to and from X key symbols. Keysym are the values returned by
/// XLookUpString.

void TGWin32::MapKeySym(UInt_t & keysym, UInt_t & xkeysym, Bool_t tox)
{
   if (tox) {
      xkeysym = GDK_VoidSymbol;
      if (keysym < 127) {
         xkeysym = keysym;
      } else if (keysym >= kKey_F1 && keysym <= kKey_F35) {
         xkeysym = GDK_F1 + (keysym - (UInt_t) kKey_F1); // function keys
      } else {
         for (int i = 0; gKeyMap[i].fKeySym; i++) {      // any other keys
            if (keysym == (UInt_t) gKeyMap[i].fKeySym) {
               xkeysym = (UInt_t) gKeyMap[i].fXKeySym;
               break;
            }
         }
      }
   } else {
      keysym = kKey_Unknown;
      // commentary in X11/keysymdef says that X codes match ASCII
      if (xkeysym < 127) {
         keysym = xkeysym;
      } else if (xkeysym >= GDK_F1 && xkeysym <= GDK_F35) {
         keysym = kKey_F1 + (xkeysym - GDK_F1);    // function keys
      } else if (xkeysym >= GDK_KP_0 && xkeysym <= GDK_KP_9) {
         keysym = kKey_0 + (xkeysym - GDK_KP_0);   // numeric keypad keys
      } else {
         for (int i = 0; gKeyMap[i].fXKeySym; i++) { // any other keys
            if (xkeysym == gKeyMap[i].fXKeySym) {
               keysym = (UInt_t) gKeyMap[i].fKeySym;
               break;
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get contents of paste buffer atom into string. If del is true delete
/// the paste buffer afterwards.

void TGWin32::GetPasteBuffer(Window_t id, Atom_t atom, TString & text,
                             Int_t & nchar, Bool_t del)
{
   if (!id) return;

   char *data;
   int nread, actual_format;

   nread = gdk_selection_property_get((GdkWindow *) id,
                                      (unsigned char **) &data,
                                      (GdkAtom *) & atom, &actual_format);

   if ((nread == 0) || (data == NULL)) {
      nchar = 0;
      return;
   }

   text.Insert(0, (const char *) data);
   nchar = 1;                   //strlen(data);
   g_free(data);

   // if (del)
   gdk_property_delete((GdkWindow *) id,
                       gdk_atom_intern("GDK_SELECTION", FALSE));
}

////////////////////////////////////////////////////////////////////////////////
/// TranslateCoordinates translates coordinates from the frame of
/// reference of one window to another. If the point is contained
/// in a mapped child of the destination, the id of that child is
/// returned as well.

void TGWin32::TranslateCoordinates(Window_t src, Window_t dest,
                                   Int_t src_x, Int_t src_y,
                                   Int_t &dest_x, Int_t &dest_y,
                                   Window_t &child)
{
   if (!src || !dest) return;

   HWND sw, dw, ch = NULL;
   POINT point;
   sw = (HWND)GDK_DRAWABLE_XID((GdkWindow *)src);
   dw = (HWND)GDK_DRAWABLE_XID((GdkWindow *)dest);
   point.x = src_x;
   point.y = src_y;
   ::MapWindowPoints(sw,        // handle of window to be mapped from
                   dw,          // handle to window to be mapped to
                   &point,      // pointer to array with points to map
                   1);          // number of structures in array
   ch = ::ChildWindowFromPointEx(dw, point, CWP_SKIPDISABLED | CWP_SKIPINVISIBLE);
   child = (Window_t)gdk_xid_table_lookup(ch);

   if (child == src) {
      child = (Window_t) 0;
   }
   dest_x = point.x;
   dest_y = point.y;
}

////////////////////////////////////////////////////////////////////////////////
/// Return geometry of window (should be called GetGeometry but signature
/// already used).

void TGWin32::GetWindowSize(Drawable_t id, Int_t & x, Int_t & y,
                            UInt_t & w, UInt_t & h)
{
   if (!id) return;

   Int_t ddum;
   if (GDK_DRAWABLE_TYPE(id) == GDK_DRAWABLE_PIXMAP) {
      x = y = 0;
      gdk_drawable_get_size((GdkDrawable *)id, (int*)&w, (int*)&h);
   }
   else {
      gdk_window_get_geometry((GdkDrawable *) id, &x, &y, (int*)&w,
                              (int*)&h, &ddum);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// FillPolygon fills the region closed by the specified path.
/// The path is closed automatically if the last point in the list does
/// not coincide with the first point. All point coordinates are
/// treated as relative to the origin. For every pair of points
/// inside the polygon, the line segment connecting them does not
/// intersect the path.

void TGWin32::FillPolygon(Window_t id, GContext_t gc, Point_t * points,
                          Int_t npnt)
{
   if (!id) return;

   gdk_win32_draw_polygon((GdkWindow *) id, (GdkGC *) gc, 1, (GdkPoint *) points, npnt);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the root window the pointer is logically on and the pointer
/// coordinates relative to the root window's origin.
/// The pointer coordinates returned to win_x and win_y are relative to
/// the origin of the specified window. In this case, QueryPointer returns
/// the child that contains the pointer, if any, or else kNone to
/// childw. QueryPointer returns the current logical state of the
/// keyboard buttons and the modifier keys in mask.

void TGWin32::QueryPointer(Window_t id, Window_t &rootw,
                           Window_t &childw, Int_t &root_x,
                           Int_t &root_y, Int_t &win_x, Int_t &win_y,
                           UInt_t &mask)
{
   if (!id) return;

   POINT currPt;
   HWND chw, window;
   UInt_t umask = 0;
   BYTE kbd[256];

   window = (HWND)GDK_DRAWABLE_XID((GdkWindow *)id);
   rootw = (Window_t)GDK_ROOT_PARENT();
   ::GetCursorPos(&currPt);
   chw = ::WindowFromPoint(currPt);
   childw = (Window_t)gdk_xid_table_lookup(chw);
   root_x = currPt.x;
   root_y = currPt.y;

   ::ScreenToClient(window, &currPt);
   win_x = currPt.x;
   win_y = currPt.y;

   ::GetKeyboardState (kbd);

   if (kbd[VK_SHIFT] & 0x80) {
      umask |= GDK_SHIFT_MASK;
   }
   if (kbd[VK_CAPITAL] & 0x80) {
      umask |= GDK_LOCK_MASK;
   }
   if (kbd[VK_CONTROL] & 0x80) {
      umask |= GDK_CONTROL_MASK;
   }
   if (kbd[VK_MENU] & 0x80) {
      umask |= GDK_MOD1_MASK;
   }
   if (kbd[VK_LBUTTON] & 0x80) {
      umask |= GDK_BUTTON1_MASK;
   }
   if (kbd[VK_MBUTTON] & 0x80) {
      umask |= GDK_BUTTON2_MASK;
   }
   if (kbd[VK_RBUTTON] & 0x80) {
      umask |= GDK_BUTTON3_MASK;
   }

   MapModifierState(mask, umask, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set foreground color in graphics context (shortcut for ChangeGC with
/// only foreground mask set).

void TGWin32::SetForeground(GContext_t gc, ULong_t foreground)
{
   GdkColor fore;
   fore.pixel = foreground;
   fore.red = GetRValue(foreground);
   fore.green = GetGValue(foreground);
   fore.blue = GetBValue(foreground);
   gdk_gc_set_foreground((GdkGC *) gc, &fore);
}

////////////////////////////////////////////////////////////////////////////////
/// Set clipping rectangles in graphics context. X, Y specify the origin
/// of the rectangles. Recs specifies an array of rectangles that define
/// the clipping mask and n is the number of rectangles.

void TGWin32::SetClipRectangles(GContext_t gc, Int_t x, Int_t y,
                                Rectangle_t * recs, Int_t n)
{
   Int_t i;
   GdkRectangle *grects = new GdkRectangle[n];

   for (i = 0; i < n; i++) {
      grects[i].x = x+recs[i].fX;
      grects[i].y = y+recs[i].fY;
      grects[i].width = recs[i].fWidth;
      grects[i].height = recs[i].fHeight;
   }

   for (i = 0; i < n; i++) {
      gdk_gc_set_clip_rectangle((GdkGC *)gc, (GdkRectangle*)recs);
   }
   delete [] grects;
}

////////////////////////////////////////////////////////////////////////////////
/// Flush (mode = 0, default) or synchronize (mode = 1) X output buffer.
/// Flush flushes output buffer. Sync flushes buffer and waits till all
/// requests have been processed by X server.

void TGWin32::Update(Int_t mode)
{
   GdiFlush();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new empty region.

Region_t TGWin32::CreateRegion()
{
   return (Region_t) gdk_region_new();
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy region.

void TGWin32::DestroyRegion(Region_t reg)
{
   gdk_region_destroy((GdkRegion *) reg);
}

////////////////////////////////////////////////////////////////////////////////
/// Union of rectangle with a region.

void TGWin32::UnionRectWithRegion(Rectangle_t * rect, Region_t src, Region_t dest)
{
   GdkRectangle r;
   r.x = rect->fX;
   r.y = rect->fY;
   r.width = rect->fWidth;
   r.height = rect->fHeight;
   dest = (Region_t) gdk_region_union_with_rect((GdkRegion *) src, &r);
}

////////////////////////////////////////////////////////////////////////////////
/// Create region for the polygon defined by the points array.
/// If winding is true use WindingRule else EvenOddRule as fill rule.

Region_t TGWin32::PolygonRegion(Point_t * points, Int_t np, Bool_t winding)
{
   return (Region_t) gdk_region_polygon((GdkPoint*)points, np,
                                 winding ? GDK_WINDING_RULE : GDK_EVEN_ODD_RULE);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the union of rega and regb and return result region.
/// The output region may be the same result region.

void TGWin32::UnionRegion(Region_t rega, Region_t regb, Region_t result)
{
   result = (Region_t) gdk_regions_union((GdkRegion *) rega, (GdkRegion *) regb);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the intersection of rega and regb and return result region.
/// The output region may be the same as the result region.

void TGWin32::IntersectRegion(Region_t rega, Region_t regb,
                              Region_t result)
{
   result = (Region_t) gdk_regions_intersect((GdkRegion *) rega,(GdkRegion *) regb);
}

////////////////////////////////////////////////////////////////////////////////
/// Subtract rega from regb.

void TGWin32::SubtractRegion(Region_t rega, Region_t regb, Region_t result)
{
   result = (Region_t)gdk_regions_subtract((GdkRegion *) rega,(GdkRegion *) regb);
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the difference between the union and intersection of
/// two regions.

void TGWin32::XorRegion(Region_t rega, Region_t regb, Region_t result)
{
   result = (Region_t) gdk_regions_xor((GdkRegion *) rega, (GdkRegion *) regb);
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the region is empty.

Bool_t TGWin32::EmptyRegion(Region_t reg)
{
   return (Bool_t) gdk_region_empty((GdkRegion *) reg);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if the point x,y is in the region.

Bool_t TGWin32::PointInRegion(Int_t x, Int_t y, Region_t reg)
{
   return (Bool_t) gdk_region_point_in((GdkRegion *) reg, x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if two regions are equal.

Bool_t TGWin32::EqualRegion(Region_t rega, Region_t regb)
{
   return (Bool_t) gdk_region_equal((GdkRegion *) rega, (GdkRegion *) regb);
}

////////////////////////////////////////////////////////////////////////////////
/// Return smallest enclosing rectangle.

void TGWin32::GetRegionBox(Region_t reg, Rectangle_t * rect)
{
   GdkRectangle r;
   gdk_region_get_clipbox((GdkRegion *) reg, &r);
   rect->fX = r.x;
   rect->fY = r.y;
   rect->fWidth = r.width;
   rect->fHeight = r.height;
}

////////////////////////////////////////////////////////////////////////////////
/// Return list of font names matching "fontname".

char **TGWin32::ListFonts(const char *fontname, Int_t /*max*/, Int_t &count)
{
   char  foundry[32], family[100], weight[32], slant[32], font_name[256];
   char  **fontlist;
   Int_t n1, fontcount = 0;

   sscanf(fontname, "-%30[^-]-%100[^-]-%30[^-]-%30[^-]-%n",
          foundry, family, weight, slant, &n1);
   // replace "medium" by "normal"
   if(!stricmp(weight,"medium")) {
      sprintf(weight,"normal");
   }
   // since all sizes are allowed with TTF, just forget it...
   sprintf(font_name, "-%s-%s-%s-%s-*", foundry, family, weight, slant);
   fontlist = gdk_font_list_new(font_name, &fontcount);
   count = fontcount;

   if (fontcount > 0) return fontlist;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::FreeFontNames(char **fontlist)
{
   gdk_font_list_free(fontlist);
}

////////////////////////////////////////////////////////////////////////////////
///

Drawable_t TGWin32::CreateImage(UInt_t width, UInt_t height)
{
   return (Drawable_t) gdk_image_new(GDK_IMAGE_SHARED, gdk_visual_get_best(),
                                     width, height);
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::GetImageSize(Drawable_t id, UInt_t &width, UInt_t &height)
{
   width  = ((GdkImage*)id)->width;
   height = ((GdkImage*)id)->height;
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::PutPixel(Drawable_t id, Int_t x, Int_t y, ULong_t pixel)
{
   if (!id) return;

   GdkImage *image = (GdkImage *)id;
   if (image->depth == 1) {
      if (pixel & 1) {
         ((UChar_t *) image->mem)[y * image->bpl + (x >> 3)] |= (1 << (7 - (x & 0x7)));
      } else {
         ((UChar_t *) image->mem)[y * image->bpl + (x >> 3)] &= ~(1 << (7 - (x & 0x7)));
      }
   } else {
      UChar_t *pixelp = (UChar_t *) image->mem + y * image->bpl + x * image->bpp;
      // Windows is always LSB, no need to check image->byte_order.
      switch (image->bpp) {
         case 4:
            pixelp[3] = 0;
         case 3:
            pixelp[2] = ((pixel >> 16) & 0xFF);
         case 2:
            pixelp[1] = ((pixel >> 8) & 0xFF);
         case 1:
            pixelp[0] = (pixel & 0xFF);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::PutImage(Drawable_t id, GContext_t gc, Drawable_t img, Int_t dx,
                       Int_t dy, Int_t x, Int_t y, UInt_t w, UInt_t h)
{
   if (!id) return;

   gdk_draw_image((GdkDrawable *) id, (GdkGC *)gc, (GdkImage *)img,
                  x, y, dx, dy, w, h);
   ::GdiFlush();
}

////////////////////////////////////////////////////////////////////////////////
///

void TGWin32::DeleteImage(Drawable_t img)
{
   gdk_image_unref((GdkImage *)img);
}

////////////////////////////////////////////////////////////////////////////////
/// Gets DIB bits
/// x, y, width, height - position of bitmap
/// returns a pointer on bitmap bits array
/// in format:
/// b1, g1, r1, 0,  b2, g2, r2, 0 ... bn, gn, rn, 0 ..
///
/// Pixels are numbered from left to right and from top to bottom.
/// By default all pixels from the whole drawable are returned.

unsigned char *TGWin32::GetColorBits(Drawable_t wid,  Int_t x, Int_t y,
                                     UInt_t width, UInt_t height)
{
   HDC hdc, memdc;
   BITMAPINFO bmi;
   HGDIOBJ oldbitmap1, oldbitmap2;
   BITMAP bm;
   HBITMAP ximage = 0;
   VOID  *bmbits = 0;
   unsigned char *ret = 0;

   if (GDK_DRAWABLE_TYPE(wid) == GDK_DRAWABLE_PIXMAP) {
      hdc = ::CreateCompatibleDC(NULL);
      oldbitmap1 = ::SelectObject(hdc, GDK_DRAWABLE_XID(wid));
      ::GetObject(GDK_DRAWABLE_XID(wid), sizeof(BITMAP), &bm);
   } else {
      hdc = ::GetDC((HWND)GDK_DRAWABLE_XID(wid));
   }
   memdc = ::CreateCompatibleDC(hdc);

   bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
   bmi.bmiHeader.biWidth = width;
   bmi.bmiHeader.biHeight = -1 * (int)(height);
   bmi.bmiHeader.biPlanes = 1;
   bmi.bmiHeader.biBitCount = 32;
   bmi.bmiHeader.biCompression = BI_RGB;
   bmi.bmiHeader.biSizeImage = 0;
   bmi.bmiHeader.biXPelsPerMeter = bmi.bmiHeader.biYPelsPerMeter = 0;
   bmi.bmiHeader.biClrUsed = 0;
   bmi.bmiHeader.biClrImportant = 0;

   ximage = ::CreateDIBSection(hdc, (BITMAPINFO *) &bmi, DIB_RGB_COLORS, &bmbits, NULL, 0);

   if (ximage && bmbits) {
      oldbitmap2 = ::SelectObject(memdc, ximage);
      ::BitBlt(memdc, x, y, width, height, hdc, 0, 0, SRCCOPY);
      ::SelectObject(memdc, oldbitmap2);
   }
   ::DeleteDC(memdc);
   if (GDK_DRAWABLE_TYPE(wid) == GDK_DRAWABLE_PIXMAP) {
      ::SelectObject(hdc, oldbitmap1);
      ::DeleteDC(hdc);
   } else {
      ::ReleaseDC((HWND)GDK_DRAWABLE_XID(wid), hdc);
   }
   if (ximage && bmbits) {
      ULong_t sz = width*height*4;
      ret = new unsigned char[sz];
      memcpy(ret, bmbits, sz);
      ::DeleteObject(ximage);
   }
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// create an image from RGB data. RGB data is in format :
/// b1, g1, r1, 0,  b2, g2, r2, 0 ... bn, gn, rn, 0 ..
///
/// Pixels are numbered from left to right and from top to bottom.
/// Note that data must be 32-bit aligned

Pixmap_t TGWin32::CreatePixmapFromData(unsigned char *bits, UInt_t width, UInt_t height)
{
   BITMAPINFO bmp_info;
   bmp_info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
   bmp_info.bmiHeader.biWidth = width;
   bmp_info.bmiHeader.biHeight = -1 * (int)(height);
   bmp_info.bmiHeader.biPlanes = 1;
   bmp_info.bmiHeader.biBitCount = 32;
   bmp_info.bmiHeader.biCompression = BI_RGB;
   bmp_info.bmiHeader.biSizeImage = 0;
   bmp_info.bmiHeader.biClrUsed = 0;
   bmp_info.bmiHeader.biXPelsPerMeter = 0L;
   bmp_info.bmiHeader.biYPelsPerMeter = 0L;
   bmp_info.bmiHeader.biClrImportant = 0;
   bmp_info.bmiColors[0].rgbRed = 0;
   bmp_info.bmiColors[0].rgbGreen = 0;
   bmp_info.bmiColors[0].rgbBlue = 0;
   bmp_info.bmiColors[0].rgbReserved = 0;

   HDC hdc = ::GetDC(NULL);
   HBITMAP hbitmap = ::CreateDIBitmap(hdc, &bmp_info.bmiHeader, CBM_INIT,
                                      (void *)bits, &bmp_info, DIB_RGB_COLORS);
   ::ReleaseDC(NULL, hdc);

   SIZE size;
   // For an obscure reason, we have to set the size of the
   // bitmap this way before to call gdk_pixmap_foreign_new
   // otherwise, it fails...
   ::SetBitmapDimensionEx(hbitmap,width, height, &size);

   return (Pixmap_t)gdk_pixmap_foreign_new((guint32)hbitmap);
}

////////////////////////////////////////////////////////////////////////////////
///register pixmap created by TGWin32GLManager

Int_t TGWin32::AddPixmap(ULong_t pix, UInt_t w, UInt_t h)
{
   HBITMAP hBmp = reinterpret_cast<HBITMAP>(pix);
   SIZE sz = SIZE();

   SetBitmapDimensionEx(hBmp, w, h, &sz);
   GdkPixmap *newPix = gdk_pixmap_foreign_new(reinterpret_cast<guint32>(hBmp));

   Int_t wid = 0;
   for(; wid < fMaxNumberOfWindows; ++wid)
      if (!fWindows[wid].open)
         break;

   if (wid == fMaxNumberOfWindows) {
      Int_t newSize = fMaxNumberOfWindows + 10;

      fWindows = (XWindow_t *)TStorage::ReAlloc(fWindows, newSize * sizeof(XWindow_t),
                                                fMaxNumberOfWindows * sizeof(XWindow_t));

      for (Int_t i = fMaxNumberOfWindows; i < newSize; ++i)
         fWindows[i].open = 0;

      fMaxNumberOfWindows = newSize;
   }

   fWindows[wid].open = 1;
   gCws = fWindows + wid;
   gCws->window = newPix;
   gCws->drawing = gCws->window;
   gCws->buffer = 0;
   gCws->double_buffer = 0;
   gCws->ispixmap = 1;
   gCws->clip = 0;
   gCws->width = w;
   gCws->height = h;
   gCws->new_colors = 0;

   return wid;
}

////////////////////////////////////////////////////////////////////////////////
/// Register a window created by Qt as a ROOT window (like InitWindow()).

Int_t TGWin32::AddWindow(ULong_t qwid, UInt_t w, UInt_t h)
{
   Int_t wid;
   // Select next free window number

 again:
   for (wid = 0; wid < fMaxNumberOfWindows; wid++) {
      if (!fWindows[wid].open) {
         fWindows[wid].open = 1;
         fWindows[wid].double_buffer = 0;
         gCws = &fWindows[wid];
         break;
      }
   }

   if (wid == fMaxNumberOfWindows) {
      int newsize = fMaxNumberOfWindows + 10;
      fWindows =
          (XWindow_t *) TStorage::ReAlloc(fWindows,
                                          newsize * sizeof(XWindow_t),
                                          fMaxNumberOfWindows *
                                          sizeof(XWindow_t));

      for (int i = fMaxNumberOfWindows; i < newsize; i++) {
         fWindows[i].open = 0;
      }

      fMaxNumberOfWindows = newsize;
      goto again;
   }

   gCws->window = gdk_window_foreign_new((guint32)qwid);

   gCws->drawing       = gCws->window;
   gCws->buffer        = 0;
   gCws->double_buffer = 0;
   gCws->ispixmap      = 0;
   gCws->clip          = 0;
   gCws->width         = w;
   gCws->height        = h;
   gCws->new_colors    = 0;

   return wid;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a window created by Qt (like CloseWindow1()).

void TGWin32::RemoveWindow(ULong_t qwid)
{
   int wid;

   SelectWindow((int)qwid);

   if (gCws->buffer) {
      gdk_pixmap_unref(gCws->buffer);
   }
   if (gCws->new_colors) {
      gdk_colormap_free_colors((GdkColormap *) fColormap,
                               (GdkColor *)gCws->new_colors, gCws->ncolors);

      delete [] gCws->new_colors;
      gCws->new_colors = 0;
   }

   GdiFlush();
   gCws->open = 0;

   if (!fWindows) return;

   // make first window in list the current window
   for (wid = 0; wid < fMaxNumberOfWindows; wid++) {
      if (fWindows[wid].open) {
         gCws = &fWindows[wid];
         return;
      }
   }
   gCws = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// The Nonrectangular Window Shape Extension adds nonrectangular
/// windows to the System.
/// This allows for making shaped (partially transparent) windows

void TGWin32::ShapeCombineMask(Window_t id, Int_t x, Int_t y, Pixmap_t mask)
{
   gdk_window_shape_combine_mask((GdkWindow *)id, (GdkBitmap *) mask, x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the width of the screen in millimeters.

UInt_t TGWin32::ScreenWidthMM() const
{
   return (UInt_t)gdk_screen_width_mm();
}

//------------------------------ Drag and Drop ---------------------------------

////////////////////////////////////////////////////////////////////////////////
/// Deletes the specified property on the specified window.

void TGWin32::DeleteProperty(Window_t win, Atom_t& prop)
{
   HWND hWnd = (HWND)GDK_DRAWABLE_XID((GdkWindow *)win);
   Atom_t atom = (Atom_t)GetProp(hWnd,(LPCTSTR)MAKELONG(prop,0));
   if (atom != 0) {
      GlobalDeleteAtom(atom);
   }
   RemoveProp(hWnd,(LPCTSTR)MAKELONG(prop,0));
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the actual type of the property, the actual format of the property,
/// and a pointer to the data actually returned.

Int_t TGWin32::GetProperty(Window_t win, Atom_t prop, Long_t offset, Long_t len,
                         Bool_t del, Atom_t req_type, Atom_t *act_type,
                         Int_t *act_format, ULong_t *nitems, ULong_t *bytes,
                         unsigned char **prop_list)
{
   HGLOBAL hdata;
   UChar_t *ptr, *data;
   UInt_t i, n, length;

   HWND hWnd = (HWND)GDK_DRAWABLE_XID((GdkWindow *)win);
   if (hWnd == NULL)
      return 0;

   Atom_t dndproxy = InternAtom("XdndProxy", kFALSE);
   Atom_t dndtypelist = InternAtom("XdndTypeList", kFALSE);

   if (prop == dndproxy)
      return 0;
   if (prop == dndtypelist) {
      *act_type = XA_ATOM;
      *prop_list = (unsigned char *)GetProp(hWnd, (LPCTSTR)MAKELONG(prop,0));
      for (n = 0; prop_list[n]; n++);
      *nitems = n;
      return n;
   }
   else {
      if (!OpenClipboard((HWND)GDK_DRAWABLE_XID((GdkWindow *)win))) {
         return 0;
      }
      hdata = GetClipboardData(CF_PRIVATEFIRST);
      ptr = (UChar_t *)GlobalLock(hdata);
      length = GlobalSize(hdata);
      data = (UChar_t *)malloc(length + 1);
      for (i = 0; i < length; i++) {
         data[i] = ptr[i];
      }
      GlobalUnlock(hdata);
      CloseClipboard();
      *prop_list = data;
      *bytes = *nitems = length;
      return length;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Changes the active cursor of the specified window.

void TGWin32::ChangeActivePointerGrab(Window_t win, UInt_t mask, Cursor_t cur)
{
   UInt_t xevmask;
   MapEventMask(mask, xevmask);
   if (cur == kNone)
      gdk_window_set_cursor((GdkWindow *) win, fCursors[kHand]);
   else
      gdk_window_set_cursor((GdkWindow *) win, (GdkCursor *)cur);
}

////////////////////////////////////////////////////////////////////////////////
/// Get Clipboard data.

void TGWin32::ConvertSelection(Window_t win, Atom_t &sel, Atom_t &target,
                             Atom_t &prop, Time_t &stamp)
{
   HGLOBAL hdata;

   static UINT gdk_selection_notify_msg =
      RegisterWindowMessage("gdk-selection-notify");
   HWND hWnd = (HWND)GDK_DRAWABLE_XID((GdkWindow *)win);
   if (!OpenClipboard((HWND)GDK_DRAWABLE_XID((GdkWindow *)win))) {
      return;
   }
   hdata = GetClipboardData(CF_PRIVATEFIRST);
   CloseClipboard();
   if (hdata == 0)
      return;
   /* Send ourselves an ersatz selection notify message so that we actually
    * fetch the data.
    */
   PostMessage(hWnd, gdk_selection_notify_msg, sel, target);
}

////////////////////////////////////////////////////////////////////////////////
/// Assigns owner of Clipboard.

Bool_t TGWin32::SetSelectionOwner(Window_t owner, Atom_t &sel)
{
   static UINT gdk_selection_request_msg =
      RegisterWindowMessage("gdk-selection-request");
   HWND hWnd = (HWND)GDK_DRAWABLE_XID((GdkWindow *)owner);
   OpenClipboard(hWnd);
   EmptyClipboard();
   CloseClipboard();
   if (owner) {
      ::PostMessage(hWnd, gdk_selection_request_msg, sel, 0);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Put data into Clipboard.

void TGWin32::ChangeProperties(Window_t id, Atom_t property, Atom_t type,
                               Int_t format, UChar_t *data, Int_t len)
{
   HGLOBAL hdata;
   Int_t i;
   UChar_t *ptr;

   if (data == 0 || len == 0)
      return;
   if (!OpenClipboard((HWND)GDK_DRAWABLE_XID((GdkWindow *)id))) {
      return;
   }
   hdata = GlobalAlloc(GMEM_MOVEABLE | GMEM_DDESHARE, len + 1);
   ptr = (UChar_t *)GlobalLock(hdata);
   for (i = 0; i < len; i++) {
      *ptr++ = *data++;
   }
   GlobalUnlock(hdata);
   SetClipboardData(CF_PRIVATEFIRST, hdata);
   CloseClipboard();
}

////////////////////////////////////////////////////////////////////////////////
/// Add the list of drag and drop types to the Window win.

void TGWin32::SetTypeList(Window_t win, Atom_t prop, Atom_t *typelist)
{
   SetProp((HWND)GDK_DRAWABLE_XID((GdkWindow *)win),
           (LPCTSTR)MAKELONG(prop,0),
           (HANDLE)typelist);
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively search in the children of Window for a Window which is at
/// location x, y and is DND aware, with a maximum depth of maxd.
/// Possibility to exclude dragwin and input.

Window_t TGWin32::FindRWindow(Window_t root, Window_t dragwin, Window_t input,
                              int x, int y, int maxd)
{
   POINT point;
   POINT cpt;
   RECT  rect;
   HWND hwnd, hwndc;
   HWND hwndt;
   Window_t win, retwin = kNone;
   Atom_t version = 0;
   Atom_t dndaware = InternAtom("XdndAware", kFALSE);

   cpt.x = x;
   cpt.y = y;
   hwnd = ::ChildWindowFromPointEx((HWND)GDK_DRAWABLE_XID((GdkWindow *)root),
                                    cpt, CWP_ALL);
   while (hwnd) {
      GetWindowRect(hwnd, &rect);
      if (PtInRect(&rect, cpt)) {
         if (GetProp(hwnd,(LPCTSTR)MAKELONG(dndaware,0))) {
            win = (Window_t) gdk_xid_table_lookup(hwnd);
            if (win && win != dragwin && win != input)
               return win;
         }
         Bool_t done = kFALSE;
         hwndt = hwnd;
         while (!done) {
            point = cpt;
            ::MapWindowPoints(NULL, hwndt, &point, 1);
            hwndc = ChildWindowFromPoint (hwndt, point);
            if (GetProp(hwnd,(LPCTSTR)MAKELONG(dndaware,0))) {
               win = (Window_t) gdk_xid_table_lookup(hwndc);
               if (win && win != dragwin && win != input)
                  return win;
            }
            if (hwndc == NULL)
               done = TRUE;
            else if (hwndc == hwndt)
               done = TRUE;
            else
               hwndt = hwndc;
            if (GetProp(hwndt,(LPCTSTR)MAKELONG(dndaware,0))) {
               win = (Window_t) gdk_xid_table_lookup(hwndt);
               if (win && win != dragwin && win != input)
                  return win;
            }
         }
      }
      hwnd = GetNextWindow(hwnd, GW_HWNDNEXT);
   }
   return kNone;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if Window win is DND aware, and knows any of the DND formats
/// passed in argument.

Bool_t TGWin32::IsDNDAware(Window_t win, Atom_t *typelist)
{
   if (!win) return kFALSE;

   Atom_t version = 0;
   Atom_t dndaware = InternAtom("XdndAware", kFALSE);
   HWND window = (HWND)GDK_DRAWABLE_XID((GdkWindow *)win);
   while (window) {
      version = (Atom_t)GetProp(window,(LPCTSTR)MAKELONG(dndaware,0));
      if (version) return kTRUE;
      window = ::GetParent(window);
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add XdndAware property and the list of drag and drop types to the
/// Window win.

void TGWin32::SetDNDAware(Window_t id, Atom_t *typelist)
{
   int n;
   if (!id) return;

   DWORD dwStyle = GetWindowLong((HWND)GDK_DRAWABLE_XID((GdkWindow *)id),
                                 GWL_EXSTYLE);
   SetWindowLong((HWND)GDK_DRAWABLE_XID((GdkWindow *)id), GWL_EXSTYLE,
                 dwStyle | WS_EX_ACCEPTFILES);
   Atom_t dndaware = InternAtom("XdndAware", kFALSE);
   SetProp((HWND)GDK_DRAWABLE_XID((GdkWindow *)id),
           (LPCTSTR)MAKELONG(dndaware,0),
           (HANDLE)XDND_PROTOCOL_VERSION);

   if (typelist == 0)
      return;
   for (n = 0; typelist[n]; n++);
   Atom_t dndtypelist = InternAtom("XdndTypeList", kFALSE);
   SetProp((HWND)GDK_DRAWABLE_XID((GdkWindow *)id),
           (LPCTSTR)MAKELONG(dndtypelist,0),
           (HANDLE)typelist);

}

////////////////////////////////////////////////////////////////////////////////
/// Set user thread id. This is used when an extra thread is created
/// to process events.

void TGWin32::SetUserThreadId(ULong_t id)
{
   if (id == 0) {
      TGWin32ProxyBase::fgMainThreadId = ((TWinNTSystem*)gSystem)->GetGUIThreadId();
   }
   else {
      TGWin32ProxyBase::fgUserThreadId = id;
   }
}

