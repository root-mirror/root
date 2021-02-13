// @(#)root/rootx:$Id$
// Author: Fons Rademakers   19/02/98
// Re-written for ROOT 6 by Timur Pocheptsov 11/03/2014.

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Rootxx                                                               //
//                                                                      //
// X11 based routines used to display the splash screen for rootx,      //
// the root front-end program.                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RConfigure.h"

#include <cassert>
#include <csignal>

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cstdio>

#include <unistd.h>
#include <pwd.h>
#include <sys/types.h>
#include <X11/Xlib.h>
#include <X11/xpm.h>
#include <X11/extensions/shape.h>

#include "Rtypes.h"
#include "snprintf.h"

#include "rootcoreteam.h"

#if defined(R__AIX) || defined(R__SOLARIS)
#   include <sys/select.h>
#endif
#include <time.h>
#include <sys/time.h>

namespace ROOT {
namespace ROOTX {

//Globals! A bunch of functions
//messing around with these globals :)
//But it's ok for our humble lil program :)
Display     *gDisplay       = 0;
Window       gLogoWindow    = 0;
Pixmap       gLogoPixmap    = 0;

//Non-rect window:
Pixmap       gShapeMask     = 0;
//
Pixmap       gCreditsPixmap = 0;
GC           gGC            = 0;
XFontStruct *gFont          = 0;
bool         gDone          = false;

//Can be accessed from a signal handler:
volatile sig_atomic_t gMayPopdown = 0;

bool         gAbout         = false;

//Non-rect window:
bool         gHasShapeExt   = false;
//
Colormap     gColormap      = Colormap();
Pixel        gBackground    = Pixel();
Pixel        gTextColor     = Pixel();

unsigned int gWidth         = 0;
unsigned int gHeight        = 0;
int          gStayUp        = 4000;   // 4 seconds

//gCreditsRect: x and y to be set at runtime (depending on shape extension and
//images).
XRectangle   gCreditsRect   = {0, 0, 455, 80}; // clip rect in logo

unsigned int gCreditsWidth  = gCreditsRect.width; // credits pixmap size
unsigned int gCreditsHeight = 0;

struct timeval gPopupTime;

const char *gConception[] = {
   "Rene Brun",
   "Fons Rademakers",
   0
};

char **gContributors = 0;


//
//Our "private API" - different aux. functions.
//
bool CreateSplashscreenWindow();
bool CreateSplashscreenImageAndShapeMask();
bool CreateFont();
bool CreateGC();
//
bool CreateCustomColors();
void FreeCustomColors();
//
void CreateTextPixmap();
//
void SetBackgroundPixmapAndMask();
void SetSplashscreenPosition();

//Event-handling and rendering:
bool StayUp(int milliSec);
void Sleep(int milliSec);
void ScrollCredits(int ypos);
void DrawVersion();
int DrawCreditItem(const char *creditItem, const char **members, int y, bool draw);
int DrawCredits(bool draw, bool extended);

//Non gui functions.
void ReadContributors();

//The final cleanup at the end.
void Cleanup();

}//namespace ROOTX
}//namespace ROOT

//The "public interface" - PopupLogo/WaitLogo.
////////////////////////////////////////////////////////////////////////////////
///Popup a splashscreen window.

void PopupLogo(bool about)
{
   using namespace ROOT::ROOTX;

   //Initialize and check what we can do:
   gDisplay = XOpenDisplay("");
   if (!gDisplay) {
      printf("PopupLogo, XOpenDisplay failed\n");
      return;
   }

   //Create a window.
   if (!CreateSplashscreenWindow()) {
      printf("PopupLogo, CreateSplashscreenWindow failed\n");
      XCloseDisplay(gDisplay);
      gDisplay = 0;
      return;
   }

   if (!CreateSplashscreenImageAndShapeMask()) {
      //If we failed to load a background image,
      //we have nothing to show at all.
      printf("PopupLogo, failed to create a background pixmap\n");
      XDestroyWindow(gDisplay, gLogoWindow);
      XCloseDisplay(gDisplay);
      gLogoWindow = 0;
      gDisplay = 0;
      return;
   }

   gAbout = about;

   //Background and (probably) shape combine mask, if we can.
   SetBackgroundPixmapAndMask();
   //Center splashscreen.
   SetSplashscreenPosition();

   if (CreateFont()) {//Request a custom font
      //Create a custom context
      if (CreateGC()) {
         //Docs say nothing about XSetFont's return values :(
         XSetFont(gDisplay, gGC, gFont->fid);
         CreateCustomColors();//Try to allocate special colors, result is not important here.
      } else {
         //GC creation failed, no need in a custom font anymore.
         XFreeFont(gDisplay, gFont);
         gFont = 0;
      }
   }

   if (gFont) {
      //We have both context and custom font,
      //so we can render credits.
      if (gAbout)
         ReadContributors();

      CreateTextPixmap();

      if (gCreditsPixmap)
         ScrollCredits(0);
      else {
         //Error while creating a pixmap, we
         //do not need our context and font anymore.
         XFreeFont(gDisplay, gFont);
         gFont = 0;
         //
         FreeCustomColors();
         //
         XFreeGC(gDisplay, gGC);
         gGC = 0;
         //We still can show an empty splashscreen with
         //our nice logo! ;) - so this error is not fatal.
      }
   }

   XSelectInput(gDisplay, gLogoWindow, ButtonPressMask | ExposureMask);
   XMapRaised(gDisplay, gLogoWindow);

   gettimeofday(&gPopupTime, 0);
}

////////////////////////////////////////////////////////////////////////////////
///From original version.
/// Main event loop waiting till time arrives to pop down logo
/// or when forced by button press event.
///From me: this even loop seems to be quite twisted and ugly.
///The original code does not work now though - event queue is
///growing and at some point you already not able to extract a
///button press event (for example) and close a splashscreen window.
///With my first version I also had problems - somehow I was missing
///the first expose event (from time to time).
///Now I empty the event queue on every iteration (selecting
///interesting events only).
///Why original version ignore this - I have no idea.
///We have at least NoExpose events in a queue (generated by XCopyArea)
///and somebody obviously have to remove them eventually.

void WaitLogo()
{
   using namespace ROOT::ROOTX;

   if (!gDisplay || !gLogoWindow)
      return;

   int ypos = 0;
   bool stopScroll = false;

   ScrollCredits(ypos);
   DrawVersion();
   XFlush(gDisplay);

   while (!gDone) {
      XEvent event = XEvent();//And here g++ does not complain??
      bool foundExpose = false;

      const int nEvents = XPending(gDisplay);
      for (int i = 0; i < nEvents; ++i) {
         XNextEvent(gDisplay, &event);
         if (event.type == ButtonPress) {
            if (gAbout && event.xbutton.button == 3) {
               stopScroll = stopScroll ? false : true;
            } else {
               gDone = true;
               break;
            }
         } else if (event.type == Expose)
            foundExpose = true;
         //We also have other events in a queue, get rid of them!
      }

      if (gDone)//Early exit.
         break;

      if (foundExpose) {
         //So we ... kind of ... paint it twice????
         ScrollCredits(ypos);
         DrawVersion();
      }

      Sleep(100);

      if (!gAbout && !StayUp(gStayUp) && gMayPopdown)
         gDone = true;

      if (gAbout && !stopScroll) {
         if (ypos == 0) Sleep(2000);
         ypos++;
         if (ypos > (int) (gCreditsHeight - gCreditsRect.height - 50))
            ypos = -int(gCreditsRect.height);
         ScrollCredits(ypos);
         XFlush(gDisplay);
      }
   }

   //Free/destroy windows, pixmaps, fonts, colors etc..
   Cleanup();

}

////////////////////////////////////////////////////////////////////////////////
///ROOT is ready to run, may pop down the logo if stay up time expires.

void PopdownLogo()
{
   ROOT::ROOTX::gMayPopdown = 1;
}

////////////////////////////////////////////////////////////////////////////////
///Close connection to X server (called by child).

void CloseDisplay()
{
   using ROOT::ROOTX::gDisplay;

   if (gDisplay)
      close(ConnectionNumber(gDisplay));
}


//Our "private API":

namespace ROOT {
namespace ROOTX {

////////////////////////////////////////////////////////////////////////////////

bool CreateSplashscreenWindow()
{
   assert(gDisplay != 0 && "CreateSplashscreenWindow, gDisplay is None");
   assert(gLogoWindow == 0 && "CreateSplashscreenWindow, gLogoWindow is already initialized");

   //TODO: check the screen???
   const int screen = DefaultScreen(gDisplay);

   gBackground = WhitePixel(gDisplay, screen);
   const Pixel foreground = BlackPixel(gDisplay, screen);

   gLogoWindow = XCreateSimpleWindow(gDisplay, DefaultRootWindow(gDisplay),
                                     -100, -100, 50, 50, 0,
                                     foreground, gBackground);

   return gLogoWindow;
}

////////////////////////////////////////////////////////////////////////////////
///1. Test if X11 supports shape combine mask.
///2.a if no - go to 3.
///2.b If yes - try to read both background image
///    and the mask. If any of operations failed - go to 3.
///    If both succeeded - return true.
///3. Try to read image without transparency (mask not needed anymore).

bool CreateSplashscreenImageAndShapeMask()
{
   bool LoadROOTSplashscreenPixmap(const char *imageFileName, bool needMask);

   assert(gDisplay != 0 && "CreateSplashscreenImageAndShapeMask, gDisplay is None");
   assert(gLogoPixmap == 0 &&
          "CreateSplashscreenImageAndShapeMask, gLogoPixmap is initialized already");
   assert(gShapeMask == 0 &&
          "CreateSplashscreenImageAndShapeMask, gShapeMask is initialized already");

   int eventBase = 0, errorBase = 0;
   gHasShapeExt = XShapeQueryExtension(gDisplay, &eventBase, &errorBase);

   if (gHasShapeExt) {
      if (!LoadROOTSplashscreenPixmap("Root6SplashEXT.xpm", true)) {//true - mask is needed.
         //We do not have a mask (or image not found)
         //and we can not call shape combine.
         gHasShapeExt = false;
         //"fallback".
         LoadROOTSplashscreenPixmap("Root6Splash.xpm", false);
      }
   } else
      LoadROOTSplashscreenPixmap("Root6Splash.xpm", false);

   return gLogoPixmap;
}

////////////////////////////////////////////////////////////////////////////////
///Splashscreen background image and (probably) a mask (if we want to use
///shape combine - non-rect window).

bool LoadROOTSplashscreenPixmap(const char *imageFileName, bool needMask)
{
   assert(imageFileName != 0 && "LoadROOTSplashscreenPixmap, parameter 'imageFileName' is null");
   //'None' instead of '0'?
   assert(gDisplay != 0 && "LoadROOTSplashscreenPixmap, gDisplay is null");
   assert(gLogoWindow != 0 && "LoadROOTSplashscreenPixmap, gLogoWindow is None");

   Screen * const screen = XDefaultScreenOfDisplay(gDisplay);
   if (!screen)
      return false;

   //TODO: Check the result?
   const int depth = PlanesOfScreen(screen);

   XWindowAttributes winAttr = XWindowAttributes();//DIE GCC!
   XGetWindowAttributes(gDisplay, gLogoWindow, &winAttr);

   XpmAttributes xpmAttr = XpmAttributes();//I have gcc!!!
   xpmAttr.valuemask    = XpmVisual | XpmColormap | XpmDepth;
   xpmAttr.visual       = winAttr.visual;
   xpmAttr.colormap     = winAttr.colormap;
   xpmAttr.depth        = winAttr.depth;

#ifdef XpmColorKey              // Not available in XPM 3.2 and earlier
   xpmAttr.valuemask |= XpmColorKey;
   if (depth > 4)
      xpmAttr.color_key = XPM_COLOR;
   else if (depth > 2)
      xpmAttr.color_key = XPM_GRAY4;
   else if (depth > 1)
      xpmAttr.color_key = XPM_GRAY;
   else if (depth == 1)
      xpmAttr.color_key = XPM_MONO;
   else
      xpmAttr.valuemask &= ~XpmColorKey;
#endif // defined(XpmColorKey)

   if (!strlen(imageFileName) || !imageFileName[0])//at least some check.
      return false;

   std::string path(100, ' ');

#ifdef ROOTICONPATH
   assert(strlen(ROOTICONPATH) != 0 &&
          "LoadROOTSplashscreenPixmap, invalid 'ROOTICONPATH'");

   path = ROOTICONPATH;
   path += "/";
   path += imageFileName;
#else
   assert(strlen(getenv("ROOTSYS")) != 0 &&
          "LoadROOTSplashscreenPixmap, the $ROOTSYS string is too long");
   path = getenv("ROOTSYS");
   path += "/icons/";
   path += imageFileName;
#endif

   Pixmap logo = None, mask = None;

   //Bertrand! Many thanks for this simple but ... smart and not so obvious (??) idea
   //with a mask :) Without you I'll have two separate xpms :)
   const int ret = XpmReadFileToPixmap(gDisplay, gLogoWindow, (char *)path.c_str(), &logo,
                                       gHasShapeExt ? &mask : 0, &xpmAttr);
   XpmFreeAttributes(&xpmAttr);

   if ((ret == XpmSuccess || ret == XpmColorError) && logo) {
      if (needMask) {
         if (mask) {
            gLogoPixmap = logo;
            gShapeMask = mask;
            return true;
         } //We need a mask, but
           //it's creation failed.
           //It's an error.
      } else {
         gLogoPixmap = logo;
         return true;
      }
   }

   printf("LoadROOTSplashscreenPixmap, failed to load a splashscreen image\n");

   if (logo)
      XFreePixmap(gDisplay, logo);
   if (mask)
      XFreePixmap(gDisplay, mask);

   return false;
}

////////////////////////////////////////////////////////////////////////////////

bool CreateFont()
{
   assert(gDisplay != 0 && "CreateFont, gDisplay is null");
   assert(gFont == 0 && "CreateFont, gFont exists already");

   gFont = XLoadQueryFont(gDisplay, "-adobe-helvetica-medium-r-*-*-10-*-*-*-*-*-iso8859-1");
   if (!gFont) {
      gFont = XLoadQueryFont(gDisplay, "fixed");
      if (!gFont)
         printf("Font creation failed\n");
   }

   return gFont;
}

////////////////////////////////////////////////////////////////////////////////

bool CreateGC()
{
   assert(gDisplay != 0 && "CreateGC, gDisplay is None");
   assert(gLogoWindow != 0 && "CreateGC, gLogoWindow is None");
   //Call it only once.
   assert(gGC == 0 && "CreateGC, gGC exists already");

   if (!(gGC = XCreateGC(gDisplay, gLogoWindow, 0, 0))) {
      printf("rootx - XCreateGC failed\n");
      return false;
   }

   return gGC;
}

////////////////////////////////////////////////////////////////////////////////

bool CreateCustomColors()
{
   assert(gDisplay != 0 && "CreateCustomColors, gDisplay is None");
   //Call it only once.
   assert(gTextColor == 0 && "CreateCustomColors, gTextColor exists already");

   const int screen = DefaultScreen(gDisplay);
   //TODO: check the result?
   gColormap = DefaultColormap(gDisplay, screen);
   if (!gColormap) {
      printf("CreateCustomColors, failed to aquire a default colormap");
      return false;
   }

   //Blue-ish color for the main text body.
   XColor textColor = XColor();
   textColor.red    = 39976;
   textColor.green  = 49151;
   textColor.blue   = 58981;
   textColor.flags  = DoRed | DoGreen | DoBlue;

   if (XAllocColor(gDisplay, gColormap, &textColor))
      gTextColor = textColor.pixel;

   return gTextColor;
}

////////////////////////////////////////////////////////////////////////////////

void FreeCustomColors()
{
   if (gTextColor) {
      assert(gDisplay != 0 && "FreeCustomColors, gDisplay is None");
      assert(gColormap != 0 && "FreeCustomColors, gColormap is None");
      XFreeColors(gDisplay, gColormap, &gTextColor, 1, 0);
      gTextColor = 0;
      gColormap = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////

void CreateTextPixmap()
{
   assert(gDisplay != 0 && "CreateTextPixmap, gDisplay is None");
   assert(gLogoWindow != 0 && "CreateTextPixmap, gLogoWindow is None");

   if (!gGC)//Something is wrong and we don't need pixmap anymore.
      return;

   Window rootWindow = Window();
   int x = 0, y = 0;
   unsigned int borderWidth = 0, depth = 0;
   unsigned int width = 0, height = 0;
   XGetGeometry(gDisplay, gLogoWindow, &rootWindow, &x, &y, &width, &height, &borderWidth, &depth);

   gCreditsHeight = DrawCredits(false, gAbout) + gCreditsRect.height + 50;
   gCreditsPixmap = XCreatePixmap(gDisplay, gLogoWindow, gCreditsWidth, gCreditsHeight, depth);

   if (gHasShapeExt)
      gCreditsRect.x = 115;
   else
      gCreditsRect.x = 15;

   assert(gHeight > 105 && "CreateTextPixmap, internal error - unexpected geometry");
   gCreditsRect.y = gHeight - 105;
}

////////////////////////////////////////////////////////////////////////////////

void SetBackgroundPixmapAndMask()
{
   assert(gDisplay != 0 && "SetBackgroundPixmapAndMask, gDisplay is None");
   assert(gLogoWindow != 0 && "SetBackgroundPixmapAndMask, gLogoWindow is None");
   assert(gLogoPixmap != 0 && "SetBackgroundPixmapAndMask, gLogoPixmap is None");

   unsigned long mask = CWBackPixmap | CWOverrideRedirect;
   XSetWindowAttributes winAttr = XSetWindowAttributes();//I had {} but f...g g++ is so STUPID!
   winAttr.background_pixmap = gLogoPixmap;
   winAttr.override_redirect = True;
   XChangeWindowAttributes(gDisplay, gLogoWindow, mask, &winAttr);

   if (gHasShapeExt) {
      assert(gShapeMask != 0 && "SetBackgroundPixmapAndMask, gShapeMask is None");
      XShapeCombineMask(gDisplay, gLogoWindow, ShapeBounding, 0, 0, gShapeMask, ShapeSet);
   }
}

////////////////////////////////////////////////////////////////////////////////

void SetSplashscreenPosition()
{
   assert(gDisplay != 0 && "SetSplashscreenPosition, gDisplay is None");
   assert(gLogoWindow != 0 && "SetSplashscreenPosition, gLogoWindow is None");

   Window rootWindow = Window();
   int x = 0, y = 0;
   unsigned int borderWidth = 0, depth = 0;
   XGetGeometry(gDisplay, gLogoPixmap, &rootWindow, &x, &y,
                &gWidth, &gHeight, &borderWidth, &depth);

   //TODO: this is wrong with multi-head display setup!
   if (Screen * const screen = XDefaultScreenOfDisplay(gDisplay)) {
      x = (WidthOfScreen(screen) - gWidth) / 2;
      y = (HeightOfScreen(screen) - gHeight) / 2;
   } else {
      //Some stupid numbers.
      x = 100;
      y = 100;
   }

   XMoveResizeWindow(gDisplay, gLogoWindow, x, y, gWidth, gHeight);
   XSync(gDisplay, False);// make sure move & resize is done before mapping
}

////////////////////////////////////////////////////////////////////////////////
///Returns false if milliSec milliseconds have passed since logo
///was popped up, true otherwise.

bool StayUp(int milliSec)
{
   struct timeval ctv, dtv, tv, ptv = gPopupTime;

   tv.tv_sec  = milliSec / 1000;
   tv.tv_usec = (milliSec % 1000) * 1000;

   gettimeofday(&ctv, 0);
   if ((dtv.tv_usec = ctv.tv_usec - ptv.tv_usec) < 0) {
      dtv.tv_usec += 1000000;
      ptv.tv_sec++;
   }
   dtv.tv_sec = ctv.tv_sec - ptv.tv_sec;

   if ((ctv.tv_usec = tv.tv_usec - dtv.tv_usec) < 0) {
      ctv.tv_usec += 1000000;
      dtv.tv_sec++;
   }
   ctv.tv_sec = tv.tv_sec - dtv.tv_sec;

   if (ctv.tv_sec < 0) return false;

   return true;
}

////////////////////////////////////////////////////////////////////////////////
///Sleep for specified amount of milliseconds.

void Sleep(int milliSec)
{
   struct timeval tv;

   tv.tv_sec  = milliSec / 1000;
   tv.tv_usec = (milliSec % 1000) * 1000;

   select(0, 0, 0, 0, &tv);
}

////////////////////////////////////////////////////////////////////////////////

void ScrollCredits(int ypos)
{
   assert(gDisplay != 0 && "ScrollCredits, gDisplay is None");

   if (!gGC || !gLogoPixmap || !gCreditsPixmap)
      return;

   XCopyArea(gDisplay, gLogoPixmap, gCreditsPixmap, gGC, gCreditsRect.x, gCreditsRect.y,
             gCreditsRect.width, gCreditsRect.height, 0, ypos);

   DrawCredits(true, gAbout);

   XRectangle crect[1];
   crect[0] = gCreditsRect;
   XSetClipRectangles(gDisplay, gGC, 0, 0, crect, 1, Unsorted);

   XCopyArea(gDisplay, gCreditsPixmap, gLogoWindow, gGC,
             0, ypos, gCreditsWidth, gCreditsHeight, gCreditsRect.x, gCreditsRect.y);

   XSetClipMask(gDisplay, gGC, None);
}

////////////////////////////////////////////////////////////////////////////////
///Draw version string.

void DrawVersion()
{
#ifndef ROOT_RELEASE
   assert(0 && "DrawVersion, 'ROOT_RELEASE' is not defined");
   return;
#endif

   assert(gDisplay != 0 && "DrawVersion, gDisplay is None");
   assert(gLogoWindow != 0 && "DrawVersion, gLogoWindow is None");
   assert(gGC != 0 && "DrawVersion, gGC is None");

   std::string version("Version ");
   version += ROOT_RELEASE;

   XDrawString(gDisplay, gLogoWindow, gGC, gWidth - 90, gHeight - 15, version.data(),
               version.length());
}

////////////////////////////////////////////////////////////////////////////////
///Draw credit item.

int DrawCreditItem(const char *creditItem, const char **members, int y, bool draw)
{
   assert(creditItem != 0 && "DrawCreditItem, parameter 'creditItem' is null");

   assert(gFont != 0 && "DrawCreditItem, gFont is None");
   assert(gDisplay != 0 && "DrawCreditItem, gDisplay is None");

   const int lineSpacing = gFont->max_bounds.ascent + gFont->max_bounds.descent;
   assert(lineSpacing > 0 && "DrawCreditItem, lineSpacing must be positive");

   std::string credit(creditItem);

   for (unsigned i = 0; members && members[i]; ++i) {
      if (i)
         credit += ", ";

      if (XTextWidth(gFont, credit.data(), credit.length()) +
          XTextWidth(gFont, members[i], strlen(members[i])) > (int) gCreditsWidth)
      {
         if (draw) {
            XDrawString(gDisplay, gCreditsPixmap, gGC, 0, y,
                        credit.data(), credit.length());
         }

         y += lineSpacing;
         credit = "   ";
      }

      credit += members[i];
   }

   if (draw) {
      XDrawString(gDisplay, gCreditsPixmap, gGC, 0, y,
                  credit.data(), credit.length());
   }

   return y;
}

////////////////////////////////////////////////////////////////////////////////
///Draw credits. If draw is true draw credits,
///otherwise just return size of all credit text.
///If extended is true draw or returns size for extended full
///credits list.

int DrawCredits(bool draw, bool extended)
{
   if (!gFont || !gGC)
      return 150;  //Size not important no text will be drawn anyway

   if (gTextColor)
      XSetForeground(gDisplay, gGC, gTextColor);//Text is white!
   else
      XSetForeground(gDisplay, gGC, gBackground);//Text is white!

   assert((draw == false || gCreditsPixmap != 0) &&
          "DrawCredits, parameter 'draw' is true, but destination pixmap is None");
   const int lineSpacing = gFont->max_bounds.ascent + gFont->max_bounds.descent;
   assert(lineSpacing > 0 && "DrawCredits, lineSpacing must be positive");

   int y = lineSpacing;
   y = DrawCreditItem("Conception: ", gConception, y, draw);
   y += 2 * lineSpacing;

   y = DrawCreditItem("Core Engineering: ", ROOT::ROOTX::gROOTCoreTeam, y, draw);

   if (extended && gContributors) {
      y += 2 * lineSpacing;
      y = DrawCreditItem("Contributors: ", (const char **)gContributors, y, draw);

      y += 2 * lineSpacing;
      y = DrawCreditItem("Our sincere thanks and apologies to anyone who deserves", 0, y, draw);
      y += lineSpacing;
      y = DrawCreditItem("credit but fails to appear in this list.", 0, y, draw);

      struct passwd *pwd = getpwuid(getuid());
      if (pwd) {
         char *name = new char [strlen(pwd->pw_gecos)+1];
         strcpy(name, pwd->pw_gecos);
         char *s = strchr(name, ',');
         if (s) *s = 0;
         char line[1024];
         if (strlen(name))
            snprintf(line, sizeof(line), "Extra special thanks go to %s,", name);
         else
            snprintf(line, sizeof(line), "Extra special thanks go to %s,", pwd->pw_name);
         delete [] name;
         y += 2*lineSpacing;
         y = DrawCreditItem(line, 0, y, draw);
         y += lineSpacing;
         y = DrawCreditItem("one of our favorite users.", 0, y, draw);
      }
   }

   if (gTextColor)
      XSetForeground(gDisplay, gGC, gBackground);//Text is white!

   return y;
}

////////////////////////////////////////////////////////////////////////////////
///Read the file $ROOTSYS/README/CREDITS for the names of the
///contributors.

void ReadContributors()
{
   //TODO: re-write without ... good old C :)
   char buf[2048];
#ifdef ROOTDOCDIR
   snprintf(buf, sizeof(buf), "%s/CREDITS", ROOTDOCDIR);
#else
   snprintf(buf, sizeof(buf), "%s/README/CREDITS", getenv("ROOTSYS"));
#endif

   gContributors = 0;

   FILE *f = fopen(buf, "r");
   if (!f) return;

   int cnt = 0;
   while (fgets(buf, sizeof(buf), f)) {
      if (!strncmp(buf, "N: ", 3)) {
         cnt++;
      }
   }
   gContributors = new char*[cnt+1];

   cnt = 0;
   rewind(f);
   while (fgets(buf, sizeof(buf), f)) {
      if (!strncmp(buf, "N: ", 3)) {
         int len = strlen(buf);
         buf[len-1] = 0;    // remove \n
         len -= 3;          // remove "N: "
         gContributors[cnt] = new char[len];
         strncpy(gContributors[cnt], buf+3, len);
         cnt++;
      }
   }

   gContributors[cnt] = 0;

   fclose(f);
}

////////////////////////////////////////////////////////////////////////////////

void Cleanup()
{
   if (gLogoWindow) {
      assert(gDisplay != 0 && "Cleanup, gDisplay is None");
      XUnmapWindow(gDisplay, gLogoWindow);
      XDestroyWindow(gDisplay, gLogoWindow);
      gLogoWindow = 0;
   }

   if (gLogoPixmap) {
      assert(gDisplay != 0 && "Cleanup, gDisplay is None");
      XFreePixmap(gDisplay, gLogoPixmap);
      gLogoPixmap = 0;
   }

   if (gShapeMask) {
      assert(gDisplay != 0 && "Cleanup, gDisplay is None");
      XFreePixmap(gDisplay, gShapeMask);
      gShapeMask = 0;
   }

   if (gCreditsPixmap) {
      assert(gDisplay != 0 && "Cleanup, gDisplay is None");
      XFreePixmap(gDisplay, gCreditsPixmap);
      gCreditsPixmap = 0;
   }

   if (gFont) {
      assert(gDisplay != 0 && "Cleanup, gDisplay is None");
      XFreeFont(gDisplay, gFont);
      gFont = 0;
   }

   //If any.
   FreeCustomColors();

   if (gGC) {
      assert(gDisplay != 0 && "Cleanup, gDisplay is None");
      XFreeGC(gDisplay, gGC);
      gGC = 0;
   }

   if (gDisplay) {
      XSync(gDisplay, False);
      XCloseDisplay(gDisplay);
      gDisplay = 0;
   }
}

}//namespace ROOTX
}//namespace ROOT
