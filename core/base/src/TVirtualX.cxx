// @(#)root/base:$Id$
// Author: Fons Rademakers   3/12/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**
\defgroup GraphicsBackends Graphics' Backends
\ingroup Graphics
Graphics' Backends interface classes.
Graphics classes interfacing ROOT graphics with the low level
native graphics backends(s) like X11, Cocoa, Win32 etc...
These classes are not meant to be used directly by ROOT users.
*/

/** \class TVirtualX
\ingroup GraphicsBackends
\ingroup Base
Semi-Abstract base class defining a generic interface to the underlying, low
level, native graphics backend (X11, Win32, MacOS, OpenGL...).
An instance of TVirtualX itself defines a batch interface to the graphics system.
*/

#include "TVirtualX.h"
#include "TString.h"


Atom_t    gWM_DELETE_WINDOW;
Atom_t    gMOTIF_WM_HINTS;
Atom_t    gROOT_MESSAGE;


TVirtualX     *gGXBatch;  //Global pointer to batch graphics interface
TVirtualX*   (*gPtr2VirtualX)() = 0; // returns pointer to global object


ClassImp(TVirtualX);


////////////////////////////////////////////////////////////////////////////////
/// Ctor of ABC

TVirtualX::TVirtualX(const char *name, const char *title) : TNamed(name, title),
      TAttLine(1,1,1),TAttFill(1,1),TAttText(11,0,1,62,0.01), TAttMarker(1,1,1),
      fDrawMode()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns gVirtualX global

TVirtualX *&TVirtualX::Instance()
{
   static TVirtualX *instance = nullptr;
   if (gPtr2VirtualX) instance = gPtr2VirtualX();
   return instance;
}

////////////////////////////////////////////////////////////////////////////////
/// The WindowAttributes_t structure is set to default.

void TVirtualX::GetWindowAttributes(Window_t /*id*/, WindowAttributes_t &attr)
{
   attr.fX = attr.fY = 0;
   attr.fWidth = attr.fHeight = 0;
   attr.fVisual   = 0;
   attr.fMapState = kIsUnmapped;
   attr.fScreen   = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Looks up the string name of a color "cname" with respect to the screen
/// associated with the specified colormap. It returns the exact color value.
/// If the color name is not in the Host Portable Character Encoding,
/// the result is implementation dependent.
///
/// \param [in] cmap    the colormap
/// \param [in] cname   the color name string; use of uppercase or lowercase
///            does not matter
/// \param [in] color   returns the exact color value for later use
///
/// The ColorStruct_t structure is set to default. Let system think we
/// could parse color.

Bool_t TVirtualX::ParseColor(Colormap_t /*cmap*/, const char * /*cname*/,
                             ColorStruct_t &color)
{
   color.fPixel = 0;
   color.fRed   = 0;
   color.fGreen = 0;
   color.fBlue  = 0;
   color.fMask  = kDoRed | kDoGreen | kDoBlue;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Allocates a read-only colormap entry corresponding to the closest RGB
/// value supported by the hardware. If no cell could be allocated it
/// returns kFALSE, otherwise kTRUE.
///
/// The pixel value is set to default. Let system think we could allocate
/// color.
///
/// \param [in] cmap    the colormap
/// \param [in] color   specifies and returns the values actually used in the cmap

Bool_t TVirtualX::AllocColor(Colormap_t /*cmap*/, ColorStruct_t &color)
{
   color.fPixel = 0;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the current RGB value for the pixel in the "color" structure
///
/// The color components are set to default.
///
/// \param [in] cmap    the colormap
/// \param [in] color   specifies and returns the RGB values for the pixel specified
///         in the structure

void TVirtualX::QueryColor(Colormap_t /*cmap*/, ColorStruct_t &color)
{
   color.fRed = color.fGreen = color.fBlue = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// The "event" is set to default event.
/// This method however, should never be called.

void TVirtualX::NextEvent(Event_t &event)
{
   event.fType   = kButtonPress;
   event.fWindow = 0;
   event.fTime   = 0;
   event.fX      = 0;
   event.fY      = 0;
   event.fXRoot  = 0;
   event.fYRoot  = 0;
   event.fState  = 0;
   event.fCode   = 0;
   event.fWidth  = 0;
   event.fHeight = 0;
   event.fCount  = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Gets contents of the paste buffer "atom" into the string "text".
/// (nchar = number of characters) If "del" is true deletes the paste
/// buffer afterwards.

void TVirtualX::GetPasteBuffer(Window_t /*id*/, Atom_t /*atom*/, TString &text,
                               Int_t &nchar, Bool_t /*del*/)
{
   text = "";
   nchar = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Initializes the X system. Returns kFALSE in case of failure.
/// It is implementation dependent.

Bool_t TVirtualX::Init(void * /*display*/)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Clears the entire area of the current window.

void TVirtualX::ClearWindow()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Deletes current window.

void TVirtualX::CloseWindow()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Deletes current pixmap.

void TVirtualX::ClosePixmap()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copies the pixmap "wid" at the position [xpos,ypos] in the current window.

void TVirtualX::CopyPixmap(Int_t /*wid*/, Int_t /*xpos*/, Int_t /*ypos*/)
{
}

////////////////////////////////////////////////////////////////////////////////
///On a HiDPI resolution it can be > 1., this means glViewport should use
///scaled width and height.

Double_t TVirtualX::GetOpenGLScalingFactor()
{
   return 1.;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates OpenGL context for window "wid"

void TVirtualX::CreateOpenGLContext(Int_t /*wid*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Deletes OpenGL context for window "wid"

void TVirtualX::DeleteOpenGLContext(Int_t /*wid*/)
{
}

////////////////////////////////////////////////////////////////////////////////
///Create window with special pixel format. Noop everywhere except Cocoa.

Window_t TVirtualX::CreateOpenGLWindow(Window_t /*parentID*/, UInt_t /*width*/, UInt_t /*height*/, const std::vector<std::pair<UInt_t, Int_t> > &/*format*/)
{
   return Window_t();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates OpenGL context for window "windowID".

Handle_t TVirtualX::CreateOpenGLContext(Window_t /*windowID*/, Handle_t /*shareWith*/)
{
   return Handle_t();
}

////////////////////////////////////////////////////////////////////////////////
/// Makes context ctx current OpenGL context.

Bool_t TVirtualX::MakeOpenGLContextCurrent(Handle_t /*ctx*/, Window_t /*windowID*/)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Asks OpenGL subsystem about the current OpenGL context.

Handle_t TVirtualX::GetCurrentOpenGLContext()
{
   return Handle_t();
}

////////////////////////////////////////////////////////////////////////////////
/// Flushes OpenGL buffer.

void TVirtualX::FlushOpenGLBuffer(Handle_t /*ctx*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Draws a box between [x1,y1] and [x2,y2] according to the "mode".
///
/// \param [in] x1,y1   left down corner
/// \param [in] x2,y2   right up corner
/// \param [in] mode    drawing mode:
///             - mode = 0 hollow  (kHollow)
///             - mode = 1 solid   (kSolid)

void TVirtualX::DrawBox(Int_t /*x1*/, Int_t /*y1*/, Int_t /*x2*/, Int_t /*y2*/,
                        EBoxMode /*mode*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Draws a cell array. The drawing is done with the pixel precision
/// if (x2-x1)/nx (or y) is not a exact pixel number the position of
/// the top right corner may be wrong.
///
/// \param [in] x1,y1   left down corner
/// \param [in] x2,y2   right up corner
/// \param [in] nx,ny   array size
/// \param [in] ic      array

void TVirtualX::DrawCellArray(Int_t /*x1*/, Int_t /*y1*/,
                              Int_t /*x2*/, Int_t /*y2*/,
                              Int_t /*nx*/, Int_t /*ny*/, Int_t * /*ic*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Fills area described by the polygon.
///
/// \param [in] n    number of points
/// \param [in] xy   list of points. xy(2,n)

void TVirtualX::DrawFillArea(Int_t /*n*/, TPoint * /*xy*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Draws a line.
///
/// \param [in] x1,y1   begin of line
/// \param [in] x2,y2   end of line

void TVirtualX::DrawLine(Int_t /*x1*/, Int_t /*y1*/, Int_t /*x2*/, Int_t /*y2*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Draws a line through all points in the list.
///
/// \param [in] n    number of points
/// \param [in] xy   list of points

void TVirtualX::DrawPolyLine(Int_t /*n*/, TPoint * /*xy*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Draws "n" markers with the current attributes at position [x,y].
///
/// \param [in] n    number of markers to draw
/// \param [in] xy   an array of x,y marker coordinates

void TVirtualX::DrawPolyMarker(Int_t /*n*/, TPoint * /*xy*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Draws a text string using current font.
///
/// \param [in] x,y     text position
/// \param [in] angle   text angle
/// \param [in] mgn     magnification factor
/// \param [in] text    text string
/// \param [in] mode    drawing mode:
///           - mode = 0 the background is not drawn (kClear)
///           - mode = 1 the background is drawn (kOpaque)

void TVirtualX::DrawText(Int_t /*x*/, Int_t /*y*/, Float_t /*angle*/,
                         Float_t /*mgn*/, const char * /*text*/,
                         ETextMode /*mode*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Draws a text string using current font.
///
/// \param [in] x,y     text position
/// \param [in] angle   text angle
/// \param [in] mgn     magnification factor
/// \param [in] text    text string
/// \param [in] mode    drawing mode:
///           - mode = 0 the background is not drawn (kClear)
///           - mode = 1 the background is drawn (kOpaque)

void TVirtualX::DrawText(Int_t /*x*/, Int_t /*y*/, Float_t /*angle*/,
                         Float_t /*mgn*/, const wchar_t * /*text*/,
                         ETextMode /*mode*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Executes the command "code" coming from the other threads (Win32)

UInt_t TVirtualX::ExecCommand(TGWin32Command * /*code*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Queries the double buffer value for the window "wid".

Int_t TVirtualX::GetDoubleBuffer(Int_t /*wid*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns character up vector.

void TVirtualX::GetCharacterUp(Float_t &chupx, Float_t &chupy)
{
   chupx = chupy = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns position and size of window "wid".
///
/// \param [in] wid    window identifier
///                    if wid < 0 the size of the display is returned
/// \param [in] x, y   returned window position
/// \param [in] w, h   returned window size

void TVirtualX::GetGeometry(Int_t /*wid*/, Int_t &x, Int_t &y,
                            UInt_t &w, UInt_t &h)
{
   x = y = 0;
   w = h = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns hostname on which the display is opened.

const char *TVirtualX::DisplayName(const char *)
{
   return "batch";
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the current native event handle.

Handle_t  TVirtualX::GetNativeEvent() const
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns pixel value associated to specified ROOT color number "cindex".

ULong_t TVirtualX::GetPixel(Color_t /*cindex*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the maximum number of planes.

void TVirtualX::GetPlanes(Int_t &nplanes)
{
   nplanes = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns RGB values for color "index".

void TVirtualX::GetRGB(Int_t /*index*/, Float_t &r, Float_t &g, Float_t &b)
{
   r = g = b = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the size of the specified character string "mess".
///
/// \param [in] w      the text width
/// \param [in] h      the text height
/// \param [in] mess   the string

void TVirtualX::GetTextExtent(UInt_t &w, UInt_t &h, char * /*mess*/)
{
   w = h = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the size of the specified character string "mess".
///
/// \param [in] w      the text width
/// \param [in] h      the text height
/// \param [in] mess   the string

void TVirtualX::GetTextExtent(UInt_t &w, UInt_t &h, wchar_t * /*mess*/)
{
   w = h = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the ascent of the current font (in pixels).
/// The ascent of a font is the distance from the baseline
/// to the highest position characters extend to

Int_t   TVirtualX::GetFontAscent() const
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default version is noop, but in principle what
/// ROOT understands as ascent is text related.

Int_t   TVirtualX::GetFontAscent(const char * /*mess*/) const
{
   return GetFontAscent();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the descent of the current font (in pixels.
/// The descent is the distance from the base line
/// to the lowest point characters extend to.

Int_t   TVirtualX::GetFontDescent() const
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default version is noop, but in principle what
/// ROOT understands as descent requires a certain text.

Int_t   TVirtualX::GetFontDescent(const char * /*mess*/) const
{
   return GetFontDescent();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the current font magnification factor

Float_t TVirtualX::GetTextMagnitude()
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns True when TrueType fonts are used

Bool_t TVirtualX::HasTTFonts() const
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the X11 window identifier.
///
/// \param [in] wid   workstation identifier (input)

Window_t TVirtualX::GetWindowID(Int_t /*wid*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a new window and return window number.
/// Returns -1 if window initialization fails.

Int_t TVirtualX::InitWindow(ULong_t /*window*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Registers a window created by Qt as a ROOT window
///
/// \param [in] qwid   window identifier
/// \param [in] w, h   the width and height, which define the window size

Int_t TVirtualX::AddWindow(ULong_t /*qwid*/, UInt_t /*w*/, UInt_t /*h*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Registers a pixmap created by TGLManager as a ROOT pixmap
///
/// \param [in] pixid  pixmap identifier
/// \param [in] w, h   the width and height, which define the pixmap size

Int_t TVirtualX::AddPixmap(ULong_t /*pixid*/, UInt_t /*w*/, UInt_t /*h*/)
{
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Removes the created by Qt window "qwid".

void TVirtualX::RemoveWindow(ULong_t /*qwid*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Moves the window "wid" to the specified x and y coordinates.
/// It does not change the window's size, raise the window, or change
/// the mapping state of the window.
///
/// \param [in] wid    window identifier
/// \param [in] x, y   coordinates, which define the new position of the window
///                    relative to its parent.

void TVirtualX::MoveWindow(Int_t /*wid*/, Int_t /*x*/, Int_t /*y*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a pixmap of the width "w" and height "h" you specified.

Int_t TVirtualX::OpenPixmap(UInt_t /*w*/, UInt_t /*h*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the pointer position.

void TVirtualX::QueryPointer(Int_t &ix, Int_t &iy)
{
   ix = iy = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// If id is NULL - loads the specified gif file at position [x0,y0] in the
/// current window. Otherwise creates pixmap from gif file

Pixmap_t TVirtualX::ReadGIF(Int_t /*x0*/, Int_t /*y0*/, const char * /*file*/,
                            Window_t /*id*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Requests Locator position.
///
/// \param [in] x,y    cursor position at moment of button press (output)
/// \param [in] ctyp   cursor type (input)
///       - ctyp = 1 tracking cross
///       - ctyp = 2 cross-hair
///       - ctyp = 3 rubber circle
///       - ctyp = 4 rubber band
///       - ctyp = 5 rubber rectangle
///
/// \param [in] mode   input mode
///       - mode = 0 request
///       - mode = 1 sample
///
/// \return
///       - in request mode:
///                     -  1 = left is pressed
///                     -  2 = middle is pressed
///                     -  3 = right is pressed
///       - in sample mode:
///                     -  11 = left is released
///                     -  12 = middle is released
///                     -  13 = right is released
///                     -  -1 = nothing is pressed or released
///                     -  -2 = leave the window
///                     - else = keycode (keyboard is pressed)

Int_t TVirtualX::RequestLocator(Int_t /*mode*/, Int_t /*ctyp*/,
                                Int_t &x, Int_t &y)
{
   x = y = 0;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Requests string: text is displayed and can be edited with Emacs-like
/// keybinding. Returns termination code (0 for ESC, 1 for RETURN)
///
/// \param [in] x,y    position where text is displayed
/// \param [in] text   displayed text (as input), edited text (as output)

Int_t TVirtualX::RequestString(Int_t /*x*/, Int_t /*y*/, char *text)
{
   if (text) *text = 0;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Rescales the window "wid".
///
/// \param [in] wid   window identifier
/// \param [in] w     the width
/// \param [in] h     the height

void TVirtualX::RescaleWindow(Int_t /*wid*/, UInt_t /*w*/, UInt_t /*h*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Resizes the specified pixmap "wid".
///
/// \param [in] wid    window identifier
/// \param [in] w, h   the width and height which define the pixmap dimensions

Int_t TVirtualX::ResizePixmap(Int_t /*wid*/, UInt_t /*w*/, UInt_t /*h*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Resizes the window "wid" if necessary.

void TVirtualX::ResizeWindow(Int_t /*wid*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Selects the window "wid" to which subsequent output is directed.

void TVirtualX::SelectWindow(Int_t /*wid*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Selects the pixmap "qpixid".

void TVirtualX::SelectPixmap(Int_t /*qpixid*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets character up vector.

void TVirtualX::SetCharacterUp(Float_t /*chupx*/, Float_t /*chupy*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Turns off the clipping for the window "wid".

void TVirtualX::SetClipOFF(Int_t /*wid*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets clipping region for the window "wid".
///
/// \param [in] wid    window identifier
/// \param [in] x, y   origin of clipping rectangle
/// \param [in] w, h   the clipping rectangle dimensions

void TVirtualX::SetClipRegion(Int_t /*wid*/, Int_t /*x*/, Int_t /*y*/,
                              UInt_t /*w*/, UInt_t /*h*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// The cursor "cursor" will be used when the pointer is in the
/// window "wid".

void TVirtualX::SetCursor(Int_t /*win*/, ECursor /*cursor*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the double buffer on/off on the window "wid".
///
/// \param [in] wid    window identifier.
///       - 999 means all opened windows.
/// \param [in] mode   the on/off switch
///       - mode = 1 double buffer is on
///       - mode = 0 double buffer is off

void TVirtualX::SetDoubleBuffer(Int_t /*wid*/, Int_t /*mode*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Turns double buffer mode off.

void TVirtualX::SetDoubleBufferOFF()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Turns double buffer mode on.

void TVirtualX::SetDoubleBufferON()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the drawing mode.
///
/// \param [in] mode    drawing mode.
///       - mode = 1 copy
///       - mode = 2 xor
///       - mode = 3 invert
///       - mode = 4 set the suitable mode for cursor echo according to the vendor

void TVirtualX::SetDrawMode(EDrawMode /*mode*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets color index "cindex" for fill areas.

void TVirtualX::SetFillColor(Color_t /*cindex*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets fill area style.
///
/// \param [in] style   compound fill area interior style
///        - style = 1000 * interiorstyle + styleindex

void TVirtualX::SetFillStyle(Style_t /*style*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets color index "cindex" for drawing lines.

void TVirtualX::SetLineColor(Color_t /*cindex*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the line type.
///
/// \param [in] n         length of the dash list
///          - n <= 0 use solid lines
///          - n >  0 use dashed lines described by dash(n)
///                 e.g. n = 4,dash = (6,3,1,3) gives a dashed-dotted line
///                 with dash length 6 and a gap of 7 between dashes
/// \param [in] dash(n)   dash segment lengths

void TVirtualX::SetLineType(Int_t /*n*/, Int_t * /*dash*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the line style.
///
/// \param [in] linestyle   line style.
///        - linestyle <= 1 solid
///        - linestyle  = 2 dashed
///        - linestyle  = 3 dotted
///        - linestyle  = 4 dashed-dotted

void TVirtualX::SetLineStyle(Style_t /*linestyle*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the line width.
///
/// \param [in] width   the line width in pixels

void TVirtualX::SetLineWidth(Width_t /*width*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets color index "cindex" for markers.

void TVirtualX::SetMarkerColor(Color_t /*cindex*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets marker size index.
///
/// \param [in] markersize   the marker scale factor

void TVirtualX::SetMarkerSize(Float_t /*markersize*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets marker style.

void TVirtualX::SetMarkerStyle(Style_t /*markerstyle*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets opacity of the current window. This image manipulation routine
/// works by adding to a percent amount of neutral to each pixels RGB.
/// Since it requires quite some additional color map entries is it
/// only supported on displays with more than > 8 color planes (> 256
/// colors).

void TVirtualX::SetOpacity(Int_t /*percent*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets color intensities the specified color index "cindex".
///
/// \param [in] cindex    color index
/// \param [in] r, g, b   the red, green, blue intensities between 0.0 and 1.0

void TVirtualX::SetRGB(Int_t /*cindex*/, Float_t /*r*/, Float_t /*g*/,
                       Float_t /*b*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the text alignment.
///
/// \param [in] talign   text alignment.
///        - talign = txalh horizontal text alignment
///        - talign = txalv vertical text alignment

void TVirtualX::SetTextAlign(Short_t /*talign*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the color index "cindex" for text.

void TVirtualX::SetTextColor(Color_t /*cindex*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets text font to specified name "fontname".This function returns 0 if
/// the specified font is found, 1 if it is not.
///
/// \param [in] fontname   font name
/// \param [in] mode       loading flag
///           - mode = 0 search if the font exist (kCheck)
///           - mode = 1 search the font and load it if it exists (kLoad)

Int_t TVirtualX::SetTextFont(char * /*fontname*/, ETextSetMode /*mode*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the current text font number.

void TVirtualX::SetTextFont(Font_t /*fontnumber*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the current text magnification factor to "mgn"

void TVirtualX::SetTextMagnitude(Float_t /*mgn*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the current text size to "textsize"

void TVirtualX::SetTextSize(Float_t /*textsize*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Set synchronisation on or off.
///
/// \param [in] mode   synchronisation on/off
///    - mode=1  on
///    - mode<>0 off

void TVirtualX::Sync(Int_t /*mode*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Updates or synchronises client and server once (not permanent).
/// according to "mode".
///
/// \param [in] mode   update mode.
///        - mode = 1 update
///        - mode = 0 sync

void TVirtualX::UpdateWindow(Int_t /*mode*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the pointer position.
///
/// \param [in] ix   new X coordinate of pointer
/// \param [in] iy   new Y coordinate of pointer
/// \param [in] id   window identifier
///
/// Coordinates are relative to the origin of the window id
/// or to the origin of the current window if id == 0.

void TVirtualX::Warp(Int_t /*ix*/, Int_t /*iy*/, Window_t /*id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Writes the current window into GIF file.
/// Returns 1 in case of success, 0 otherwise.

Int_t TVirtualX::WriteGIF(char * /*name*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Writes the pixmap "wid" in the bitmap file "pxname".
///
/// \param [in] wid      the pixmap address
/// \param [in] w, h     the width and height of the pixmap.
/// \param [in] pxname   the file name

void TVirtualX::WritePixmap(Int_t /*wid*/, UInt_t /*w*/, UInt_t /*h*/,
                            char * /*pxname*/)
{
}


//---- Methods used for GUI -----
////////////////////////////////////////////////////////////////////////////////
/// Maps the window "id" and all of its subwindows that have had map
/// requests. This function has no effect if the window is already mapped.

void TVirtualX::MapWindow(Window_t /*id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Maps all subwindows for the specified window "id" in top-to-bottom
/// stacking order.

void TVirtualX::MapSubwindows(Window_t /*id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Maps the window "id" and all of its subwindows that have had map
/// requests on the screen and put this window on the top of of the
/// stack of all windows.

void TVirtualX::MapRaised(Window_t /*id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Unmaps the specified window "id". If the specified window is already
/// unmapped, this function has no effect. Any child window will no longer
/// be visible (but they are still mapped) until another map call is made
/// on the parent.

void TVirtualX::UnmapWindow(Window_t /*id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destroys the window "id" as well as all of its subwindows.
/// The window should never be referenced again. If the window specified
/// by the "id" argument is mapped, it is unmapped automatically.

void TVirtualX::DestroyWindow(Window_t /*id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// The DestroySubwindows function destroys all inferior windows of the
/// specified window, in bottom-to-top stacking order.

void TVirtualX::DestroySubwindows(Window_t /*id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Raises the specified window to the top of the stack so that no
/// sibling window obscures it.

void TVirtualX::RaiseWindow(Window_t /*id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Lowers the specified window "id" to the bottom of the stack so
/// that it does not obscure any sibling windows.

void TVirtualX::LowerWindow(Window_t /*id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Moves the specified window to the specified x and y coordinates.
/// It does not change the window's size, raise the window, or change
/// the mapping state of the window.
///
/// \param [in] id     window identifier
/// \param [in] x, y   coordinates, which define the new position of the window
///                    relative to its parent.

void TVirtualX::MoveWindow(Window_t /*id*/, Int_t /*x*/, Int_t /*y*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Changes the size and location of the specified window "id" without
/// raising it.
///
/// \param [in] id     window identifier
/// \param [in] x, y   coordinates, which define the new position of the window
///                    relative to its parent.
/// \param [in] w, h   the width and height, which define the interior size of
///                    the window

void TVirtualX::MoveResizeWindow(Window_t /*id*/, Int_t /*x*/, Int_t /*y*/,
                                   UInt_t /*w*/, UInt_t /*h*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Changes the width and height of the specified window "id", not
/// including its borders. This function does not change the window's
/// upper-left coordinate.
///
/// \param [in] id     window identifier
/// \param [in] w, h   the width and height, which are the interior dimensions of
///                    the window after the call completes.

void TVirtualX::ResizeWindow(Window_t /*id*/, UInt_t /*w*/, UInt_t /*h*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Iconifies the window "id".

void TVirtualX::IconifyWindow(Window_t /*id*/)
{
}
////////////////////////////////////////////////////////////////////////////////
/// Notify the low level GUI layer ROOT requires "tgwindow" to be
/// updated
///
/// Returns kTRUE if the notification was desirable and it was sent
///
/// At the moment only Qt4 layer needs that
///
/// One needs explicitly cast the first parameter to TGWindow to make
/// it working in the implementation.
///
/// One needs to process the notification to confine
/// all paint operations within "expose" / "paint" like low level event
/// or equivalent

Bool_t TVirtualX::NeedRedraw(ULong_t /*tgwindow*/, Bool_t /*force*/)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// If the specified window is mapped, ReparentWindow automatically
/// performs an UnmapWindow request on it, removes it from its current
/// position in the hierarchy, and inserts it as the child of the specified
/// parent. The window is placed in the stacking order on top with respect
/// to sibling windows.

void TVirtualX::ReparentWindow(Window_t /*id*/, Window_t /*pid*/,
                               Int_t /*x*/, Int_t /*y*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the background of the window "id" to the specified color value
/// "color". Changing the background does not cause the window contents
/// to be changed.

void TVirtualX::SetWindowBackground(Window_t /*id*/, ULong_t /*color*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the background pixmap of the window "id" to the specified
/// pixmap "pxm".

void TVirtualX::SetWindowBackgroundPixmap(Window_t /*id*/, Pixmap_t /*pxm*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Creates an unmapped subwindow for a specified parent window and returns
/// the created window. The created window is placed on top in the stacking
/// order with respect to siblings. The coordinate system has the X axis
/// horizontal and the Y axis vertical with the origin [0,0] at the
/// upper-left corner. Each window and pixmap has its own coordinate system.
///
/// \param [in] parent   the parent window
/// \param [in] x, y     coordinates, the top-left outside corner of the window's
///             borders; relative to the inside of the parent window's borders
/// \param [in] w, h     width and height of the created window; do not include the
///             created window's borders
/// \param [in] border   the border pixel value of the window
/// \param [in] depth    the window's depth
/// \param [in] clss     the created window's class; can be InputOutput, InputOnly, or
///             CopyFromParent
/// \param [in] visual   the visual type
/// \param [in] attr     the structure from which the values are to be taken.
/// \param [in] wtype    the window type

Window_t TVirtualX::CreateWindow(Window_t /*parent*/, Int_t /*x*/, Int_t /*y*/,
                                 UInt_t /*w*/, UInt_t /*h*/,
                                 UInt_t /*border*/, Int_t /*depth*/,
                                 UInt_t /*clss*/, void * /*visual*/,
                                 SetWindowAttributes_t * /*attr*/,
                                 UInt_t /*wtype*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Opens connection to display server (if such a thing exist on the
/// current platform). The encoding and interpretation of the display
/// name.
///
/// On X11 this method returns on success the X display socket descriptor
/// >0, 0 in case of batch mode, and <0 in case of failure (cannot connect
/// to display dpyName).

Int_t TVirtualX::OpenDisplay(const char * /*dpyName*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Closes connection to display server and destroys all windows.

void TVirtualX::CloseDisplay()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns handle to display (might be useful in some cases where
/// direct X11 manipulation outside of TVirtualX is needed, e.g. GL
/// interface).

Display_t TVirtualX::GetDisplay() const
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns handle to visual.
///
/// Might be useful in some cases where direct X11 manipulation outside
/// of TVirtualX is needed, e.g. GL interface.

Visual_t TVirtualX::GetVisual() const
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns screen number.
///
/// Might be useful in some cases where direct X11 manipulation outside
/// of TVirtualX is needed, e.g. GL interface.

Int_t TVirtualX::GetScreen() const
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns depth of screen (number of bit planes).
/// Equivalent to GetPlanes().

Int_t TVirtualX::GetDepth() const
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns handle to colormap.
///
/// Might be useful in some cases where direct X11 manipulation outside
/// of TVirtualX is needed, e.g. GL interface.

Colormap_t TVirtualX::GetColormap() const
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns handle to the default root window created when calling
/// XOpenDisplay().

Window_t TVirtualX::GetDefaultRootWindow() const
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the atom identifier associated with the specified "atom_name"
/// string. If "only_if_exists" is False, the atom is created if it does
/// not exist. If the atom name is not in the Host Portable Character
/// Encoding, the result is implementation dependent. Uppercase and
/// lowercase matter; the strings "thing", "Thing", and "thinG" all
/// designate different atoms.

Atom_t  TVirtualX::InternAtom(const char * /*atom_name*/,
                              Bool_t /*only_if_exist*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the parent of the window "id".

Window_t TVirtualX::GetParent(Window_t /*id*/) const
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Provides the most common way for accessing a font: opens (loads) the
/// specified font and returns a pointer to the appropriate FontStruct_t
/// structure. If the font does not exist, it returns NULL.

FontStruct_t TVirtualX::LoadQueryFont(const char * /*font_name*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the font handle of the specified font structure "fs".

FontH_t TVirtualX::GetFontHandle(FontStruct_t /*fs*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Explicitly deletes the font structure "fs" obtained via LoadQueryFont().

void TVirtualX::DeleteFont(FontStruct_t /*fs*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a graphics context using the provided GCValues_t *gval structure.
/// The mask data member of gval specifies which components in the GC are
/// to be set using the information in the specified values structure.
/// It returns a graphics context handle GContext_t that can be used with any
/// destination drawable or O if the creation falls.

GContext_t TVirtualX::CreateGC(Drawable_t /*id*/, GCValues_t * /*gval*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Changes the components specified by the mask in gval for the specified GC.
///
/// \param [in] gc     specifies the GC to be changed
/// \param [in] gval   specifies the mask and the values to be set
///
/// (see also the GCValues_t structure)

void TVirtualX::ChangeGC(GContext_t /*gc*/, GCValues_t * /*gval*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copies the specified components from the source GC "org" to the
/// destination GC "dest". The "mask" defines which component to copy
/// and it is a data member of GCValues_t.

void TVirtualX::CopyGC(GContext_t /*org*/, GContext_t /*dest*/, Mask_t /*mask*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Deletes the specified GC "gc".

void TVirtualX::DeleteGC(GContext_t /*gc*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Creates the specified cursor. (just return cursor from cursor pool).
/// The cursor can be:
/// ~~~ {.cpp}
/// kBottomLeft, kBottomRight, kTopLeft,  kTopRight,
/// kBottomSide, kLeftSide,    kTopSide,  kRightSide,
/// kMove,       kCross,       kArrowHor, kArrowVer,
/// kHand,       kRotate,      kPointer,  kArrowRight,
/// kCaret,      kWatch
/// ~~~

Cursor_t TVirtualX::CreateCursor(ECursor /*cursor*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the cursor "curid" to be used when the pointer is in the
/// window "id".

void TVirtualX::SetCursor(Window_t /*id*/, Cursor_t /*curid*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a pixmap of the specified width and height and returns
/// a pixmap ID that identifies it.

Pixmap_t TVirtualX::CreatePixmap(Drawable_t /*id*/, UInt_t /*w*/, UInt_t /*h*/)
{
   return kNone;
}
////////////////////////////////////////////////////////////////////////////////
/// Creates a pixmap from bitmap data of the width, height, and depth you
/// specified and returns a pixmap that identifies it. The width and height
/// arguments must be nonzero. The depth argument must be one of the depths
/// supported by the screen of the specified drawable.
///
/// \param [in] id              specifies which screen the pixmap is created on
/// \param [in] bitmap          the data in bitmap format
/// \param [in] width, height   define the dimensions of the pixmap
/// \param [in] forecolor       the foreground pixel values to use
/// \param [in] backcolor       the background pixel values to use
/// \param [in] depth           the depth of the pixmap

Pixmap_t TVirtualX::CreatePixmap(Drawable_t /*id*/, const char * /*bitmap*/,
                                 UInt_t /*width*/, UInt_t /*height*/,
                                 ULong_t /*forecolor*/, ULong_t /*backcolor*/,
                                 Int_t /*depth*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a bitmap (i.e. pixmap with depth 1) from the bitmap data.
///
/// \param [in] id              specifies which screen the pixmap is created on
/// \param [in] bitmap          the data in bitmap format
/// \param [in] width, height   define the dimensions of the pixmap

Pixmap_t TVirtualX::CreateBitmap(Drawable_t /*id*/, const char * /*bitmap*/,
                                 UInt_t /*width*/, UInt_t /*height*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Explicitly deletes the pixmap resource "pmap".

void TVirtualX::DeletePixmap(Pixmap_t /*pmap*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a picture pict from data in file "filename". The picture
/// attributes "attr" are used for input and output. Returns kTRUE in
/// case of success, kFALSE otherwise. If the mask "pict_mask" does not
/// exist it is set to kNone.

Bool_t TVirtualX::CreatePictureFromFile(Drawable_t /*id*/,
                                        const char * /*filename*/,
                                        Pixmap_t &/*pict*/,
                                        Pixmap_t &/*pict_mask*/,
                                        PictureAttributes_t &/*attr*/)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a picture pict from data in bitmap format. The picture
/// attributes "attr" are used for input and output. Returns kTRUE in
/// case of success, kFALSE otherwise. If the mask "pict_mask" does not
/// exist it is set to kNone.

Bool_t TVirtualX::CreatePictureFromData(Drawable_t /*id*/, char ** /*data*/,
                                        Pixmap_t &/*pict*/,
                                        Pixmap_t &/*pict_mask*/,
                                        PictureAttributes_t & /*attr*/)
{
   return kFALSE;
}
////////////////////////////////////////////////////////////////////////////////
/// Reads picture data from file "filename" and store it in "ret_data".
/// Returns kTRUE in case of success, kFALSE otherwise.

Bool_t TVirtualX::ReadPictureDataFromFile(const char * /*filename*/,
                                          char *** /*ret_data*/)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete picture data created by the function ReadPictureDataFromFile.

void TVirtualX::DeletePictureData(void * /*data*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the dash-offset and dash-list attributes for dashed line styles
/// in the specified GC. There must be at least one element in the
/// specified dash_list. The initial and alternating elements (second,
/// fourth, and so on) of the dash_list are the even dashes, and the
/// others are the odd dashes. Each element in the "dash_list" array
/// specifies the length (in pixels) of a segment of the pattern.
///
/// \param [in] gc          specifies the GC (see GCValues_t structure)
/// \param [in] offset      the phase of the pattern for the dashed line-style you
///                want to set for the specified GC.
/// \param [in] dash_list   the dash-list for the dashed line-style you want to set
///                for the specified GC
/// \param [in] n           the number of elements in dash_list
/// (see also the GCValues_t structure)

void TVirtualX::SetDashes(GContext_t /*gc*/, Int_t /*offset*/,
                          const char * /*dash_list*/, Int_t /*n*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Frees color cell with specified pixel value.

void TVirtualX::FreeColor(Colormap_t /*cmap*/, ULong_t /*pixel*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the number of events that have been received from the X server
/// but have not been removed from the event queue.

Int_t TVirtualX::EventsPending()
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the sound bell. Percent is loudness from -100% to 100%.

void TVirtualX::Bell(Int_t /*percent*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Combines the specified rectangle of "src" with the specified rectangle
/// of "dest" according to the "gc".
///
/// \param [in] src              source rectangle
/// \param [in] dest             destination rectangle
/// \param [in] gc               graphics context
/// \param [in] src_x, src_y     specify the x and y coordinates, which are relative
///                              to the origin of the source rectangle and specify
///                              upper-left corner.
/// \param [in] width, height    the width and height, which are the dimensions of both
///                              the source and destination rectangles
/// \param [in] dest_x, dest_y   specify the upper-left corner of the destination
///                              rectangle
///
/// GC components in use: function, plane-mask, subwindow-mode,
/// graphics-exposure, clip-x-origin, clip-y-origin, and clip-mask.
/// (see also the GCValues_t structure)

void TVirtualX::CopyArea(Drawable_t /*src*/, Drawable_t /*dest*/,
                         GContext_t /*gc*/, Int_t /*src_x*/, Int_t /*src_y*/,
                         UInt_t /*width*/, UInt_t /*height*/,
                         Int_t /*dest_x*/, Int_t /*dest_y*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Changes the attributes of the specified window "id" according the
/// values provided in "attr". The mask data member of "attr" specifies
/// which window attributes are defined in the attributes argument.
/// This mask is the bitwise inclusive OR of the valid attribute mask
/// bits; if it is zero, the attributes are ignored.

void TVirtualX::ChangeWindowAttributes(Window_t /*id*/,
                                       SetWindowAttributes_t * /*attr*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Alters the property for the specified window and causes the X server
/// to generate a PropertyNotify event on that window.
///
/// \param [in] id         the window whose property you want to change
/// \param [in] property   specifies the property name
/// \param [in] type       the type of the property; the X server does not
///               interpret the type but simply passes it back to
///               an application that might ask about the window
///               properties
/// \param [in] data       the property data
/// \param [in] len        the length of the specified data format

void TVirtualX::ChangeProperty(Window_t /*id*/, Atom_t /*property*/,
                               Atom_t /*type*/, UChar_t * /*data*/,
                               Int_t /*len*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Uses the components of the specified GC to draw a line between the
/// specified set of points (x1, y1) and (x2, y2).
///
/// GC components in use: function, plane-mask, line-width, line-style,
/// cap-style, fill-style, subwindow-mode, clip-x-origin, clip-y-origin,
/// and clip-mask.
///
/// GC mode-dependent components: foreground, background, tile, stipple,
/// tile-stipple-x-origin, tile-stipple-y-origin, dash-offset, dash-list.
/// (see also the GCValues_t structure)

void TVirtualX::DrawLine(Drawable_t /*id*/, GContext_t /*gc*/,
                         Int_t /*x1*/, Int_t /*y1*/, Int_t /*x2*/, Int_t /*y2*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Paints a rectangular area in the specified window "id" according to
/// the specified dimensions with the window's background pixel or pixmap.
///
/// \param [in] id   specifies the window
/// \param [in] x, y   coordinates, which are relative to the origin
/// \param [in] w, h   the width and height which define the rectangle dimensions

void TVirtualX::ClearArea(Window_t /*id*/, Int_t /*x*/, Int_t /*y*/,
                          UInt_t /*w*/, UInt_t /*h*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Check if there is for window "id" an event of type "type". If there
/// is it fills in the event structure and return true. If no such event
/// return false.

Bool_t TVirtualX::CheckEvent(Window_t /*id*/, EGEventType /*type*/,
                             Event_t &/*ev*/)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Specifies the event "ev" is to be sent to the window "id".
/// This function requires you to pass an event mask.

void TVirtualX::SendEvent(Window_t /*id*/, Event_t * /*ev*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Force processing of event, sent by SendEvent before.

void TVirtualX::DispatchClientMessage(UInt_t /*messageID*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Tells WM to send message when window is closed via WM.

void TVirtualX::WMDeleteNotify(Window_t /*id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Turns key auto repeat on (kTRUE) or off (kFALSE).

void TVirtualX::SetKeyAutoRepeat(Bool_t /*on = kTRUE*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Establishes a passive grab on the keyboard. In the future, the
/// keyboard is actively grabbed, the last-keyboard-grab time is set
/// to the time at which the key was pressed (as transmitted in the
/// KeyPress event), and the KeyPress event is reported if all of the
/// following conditions are true:
///
///  - the keyboard is not grabbed and the specified key (which can
///    itself be a modifier key) is logically pressed when the
///    specified modifier keys are logically down, and no other
///    modifier keys are logically down;
///  - either the grab window "id" is an ancestor of (or is) the focus
///    window, or "id" is a descendant of the focus window and contains
///    the pointer;
///  - a passive grab on the same key combination does not exist on any
///    ancestor of grab_window
///
/// \param [in] id         window id
/// \param [in] keycode    specifies the KeyCode or AnyKey
/// \param [in] modifier   specifies the set of keymasks or AnyModifier; the mask is
///               the bitwise inclusive OR of the valid keymask bits
/// \param [in] grab       a switch between grab/ungrab key
///               grab = kTRUE  grab the key and modifier
///               grab = kFALSE ungrab the key and modifier

void TVirtualX::GrabKey(Window_t /*id*/, Int_t /*keycode*/, UInt_t /*modifier*/,
                        Bool_t /*grab = kTRUE*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Establishes a passive grab on a certain mouse button. That is, when a
/// certain mouse button is hit while certain modifier's (Shift, Control,
/// Meta, Alt) are active then the mouse will be grabbed for window id.
/// When grab is false, ungrab the mouse button for this button and modifier.

void TVirtualX::GrabButton(Window_t /*id*/, EMouseButton /*button*/,
                           UInt_t /*modifier*/, UInt_t /*evmask*/,
                           Window_t /*confine*/, Cursor_t /*cursor*/,
                           Bool_t /*grab = kTRUE*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Establishes an active pointer grab. While an active pointer grab is in
/// effect, further pointer events are only reported to the grabbing
/// client window.

void TVirtualX::GrabPointer(Window_t /*id*/, UInt_t /*evmask*/,
                            Window_t /*confine*/, Cursor_t /*cursor*/,
                            Bool_t /*grab = kTRUE*/,
                            Bool_t /*owner_events = kTRUE*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the window name.

void TVirtualX::SetWindowName(Window_t /*id*/, char * /*name*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the window icon name.

void TVirtualX::SetIconName(Window_t /*id*/, char * /*name*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the icon name pixmap.

void TVirtualX::SetIconPixmap(Window_t /*id*/, Pixmap_t /*pix*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the windows class and resource name.

void TVirtualX::SetClassHints(Window_t /*id*/, char * /*className*/,
                              char * /*resourceName*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets decoration style.

void TVirtualX::SetMWMHints(Window_t /*id*/, UInt_t /*value*/, UInt_t /*funcs*/,
                            UInt_t /*input*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Tells the window manager the desired position [x,y] of window "id".

void TVirtualX::SetWMPosition(Window_t /*id*/, Int_t /*x*/, Int_t /*y*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Tells window manager the desired size of window "id".
///
/// \param [in] id   window identifier
/// \param [in] w    the width
/// \param [in] h    the height

void TVirtualX::SetWMSize(Window_t /*id*/, UInt_t /*w*/, UInt_t /*h*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Gives the window manager minimum and maximum size hints of the window
/// "id". Also specify via "winc" and "hinc" the resize increments.
///
/// \param [in] id           window identifier
/// \param [in] wmin, hmin   specify the minimum window size
/// \param [in] wmax, hmax   specify the maximum window size
/// \param [in] winc, hinc   define an arithmetic progression of sizes into which
///                          the window to be resized (minimum to maximum)

void TVirtualX::SetWMSizeHints(Window_t /*id*/, UInt_t /*wmin*/, UInt_t /*hmin*/,
                               UInt_t /*wmax*/, UInt_t /*hmax*/,
                               UInt_t /*winc*/, UInt_t /*hinc*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the initial state of the window "id": either kNormalState
/// or kIconicState.

void TVirtualX::SetWMState(Window_t /*id*/, EInitialState /*state*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Tells window manager that the window "id" is a transient window
/// of the window "main_id". A window manager may decide not to decorate
/// a transient window or may treat it differently in other ways.

void TVirtualX::SetWMTransientHint(Window_t /*id*/, Window_t /*main_id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Each character image, as defined by the font in the GC, is treated as an
/// additional mask for a fill operation on the drawable.
///
/// \param [in] id     the drawable
/// \param [in] gc     the GC
/// \param [in] x, y   coordinates, which are relative to the origin of the specified
///           drawable and define the origin of the first character
/// \param [in] s      the character string
/// \param [in] len    the number of characters in the string argument
///
/// GC components in use: function, plane-mask, fill-style, font,
/// subwindow-mode, clip-x-origin, clip-y-origin, and clip-mask.
/// GC mode-dependent components: foreground, background, tile, stipple,
/// tile-stipple-x-origin, and tile-stipple-y-origin.
/// (see also the GCValues_t structure)

void TVirtualX::DrawString(Drawable_t /*id*/, GContext_t /*gc*/, Int_t /*x*/,
                           Int_t /*y*/, const char * /*s*/, Int_t /*len*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Return length of the string "s" in pixels. Size depends on font.

Int_t TVirtualX::TextWidth(FontStruct_t /*font*/, const char * /*s*/,
                             Int_t /*len*/)
{
   return 5;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the font properties.

void TVirtualX::GetFontProperties(FontStruct_t /*font*/, Int_t &max_ascent,
                                  Int_t &max_descent)
{
   max_ascent = 5;
   max_descent = 5;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the components specified by the mask in "gval" for the
/// specified GC "gc" (see also the GCValues_t structure)

void TVirtualX::GetGCValues(GContext_t /*gc*/, GCValues_t &gval)
{
   gval.fMask = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the font associated with the graphics context gc

FontStruct_t TVirtualX::GetGCFont(GContext_t /*gc*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieves the associated font structure of the font specified font
/// handle "fh".
///
/// Free returned FontStruct_t using FreeFontStruct().

FontStruct_t TVirtualX::GetFontStruct(FontH_t /*fh*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Frees the font structure "fs". The font itself will be freed when
/// no other resource references it.

void TVirtualX::FreeFontStruct(FontStruct_t /*fs*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Clears the entire area in the specified window and it is equivalent to
/// ClearArea(id, 0, 0, 0, 0)

void TVirtualX::ClearWindow(Window_t /*id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Converts the "keysym" to the appropriate keycode. For example,
/// keysym is a letter and keycode is the matching keyboard key (which
/// is dependent on the current keyboard mapping). If the specified
/// "keysym" is not defined for any keycode, returns zero.

Int_t TVirtualX::KeysymToKeycode(UInt_t /*keysym*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills the specified rectangle defined by [x,y] [x+w,y] [x+w,y+h] [x,y+h].
/// using the GC you specify.
///
/// GC components in use are: function, plane-mask, fill-style,
/// subwindow-mode, clip-x-origin, clip-y-origin, clip-mask.
/// GC mode-dependent components: foreground, background, tile, stipple,
/// tile-stipple-x-origin, and tile-stipple-y-origin.
/// (see also the GCValues_t structure)

void TVirtualX::FillRectangle(Drawable_t /*id*/, GContext_t /*gc*/,
                              Int_t /*x*/, Int_t /*y*/,
                              UInt_t /*w*/, UInt_t /*h*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Draws rectangle outlines of [x,y] [x+w,y] [x+w,y+h] [x,y+h]
///
/// GC components in use: function, plane-mask, line-width, line-style,
/// cap-style, join-style, fill-style, subwindow-mode, clip-x-origin,
/// clip-y-origin, clip-mask.
/// GC mode-dependent components: foreground, background, tile, stipple,
/// tile-stipple-x-origin, tile-stipple-y-origin, dash-offset, dash-list.
/// (see also the GCValues_t structure)

void TVirtualX::DrawRectangle(Drawable_t /*id*/, GContext_t /*gc*/,
                              Int_t /*x*/, Int_t /*y*/,
                              UInt_t /*w*/, UInt_t /*h*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Draws multiple line segments. Each line is specified by a pair of points.
///
/// \param [in] id     Drawable identifier
/// \param [in] gc     graphics context
/// \param [in] *seg   specifies an array of segments
/// \param [in] nseg   specifies the number of segments in the array
///
/// GC components in use: function, plane-mask, line-width, line-style,
/// cap-style, join-style, fill-style, subwindow-mode, clip-x-origin,
/// clip-y-origin, clip-mask.
///
/// GC mode-dependent components: foreground, background, tile, stipple,
/// tile-stipple-x-origin, tile-stipple-y-origin, dash-offset, and dash-list.
/// (see also the GCValues_t structure)

void TVirtualX::DrawSegments(Drawable_t /*id*/, GContext_t /*gc*/,
                             Segment_t * /*seg*/, Int_t /*nseg*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Defines which input events the window is interested in. By default
/// events are propagated up the window stack. This mask can also be
/// set at window creation time via the SetWindowAttributes_t::fEventMask
/// attribute.

void TVirtualX::SelectInput(Window_t /*id*/, UInt_t /*evmask*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the window id of the window having the input focus.

Window_t TVirtualX::GetInputFocus()
{
   return kNone;
}

////////////////////////////////////////////////////////////////////////////////
/// Changes the input focus to specified window "id".

void TVirtualX::SetInputFocus(Window_t /*id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the window id of the current owner of the primary selection.
/// That is the window in which, for example some text is selected.

Window_t TVirtualX::GetPrimarySelectionOwner()
{
   return kNone;
}

////////////////////////////////////////////////////////////////////////////////
/// Makes the window "id" the current owner of the primary selection.
/// That is the window in which, for example some text is selected.

void TVirtualX::SetPrimarySelectionOwner(Window_t /*id*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Causes a SelectionRequest event to be sent to the current primary
/// selection owner. This event specifies the selection property
/// (primary selection), the format into which to convert that data before
/// storing it (target = XA_STRING), the property in which the owner will
/// place the information (sel_property), the window that wants the
/// information (id), and the time of the conversion request (when).
/// The selection owner responds by sending a SelectionNotify event, which
/// confirms the selected atom and type.

void TVirtualX::ConvertPrimarySelection(Window_t /*id*/, Atom_t /*clipboard*/,
                                        Time_t /*when*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Converts the keycode from the event structure to a key symbol (according
/// to the modifiers specified in the event structure and the current
/// keyboard mapping). In "buf" a null terminated ASCII string is returned
/// representing the string that is currently mapped to the key code.
///
/// \param [in] event    specifies the event structure to be used
/// \param [in] buf      returns the translated characters
/// \param [in] buflen   the length of the buffer
/// \param [in] keysym   returns the "keysym" computed from the event
///             if this argument is not NULL

void TVirtualX::LookupString(Event_t * /*event*/, char * /*buf*/,
                             Int_t /*buflen*/, UInt_t &keysym)
{
   keysym = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Translates coordinates in one window to the coordinate space of another
/// window. It takes the "src_x" and "src_y" coordinates relative to the
/// source window's origin and returns these coordinates to "dest_x" and
/// "dest_y" relative to the destination window's origin.
///
/// \param [in] src              the source window
/// \param [in] dest             the destination window
/// \param [in] src_x, src_y     coordinates within the source window
/// \param [in] dest_x, dest_y   coordinates within the destination window
/// \param [in] child            returns the child of "dest" if the coordinates
///                     are contained in a mapped child of the destination
///                     window; otherwise, child is set to 0

void TVirtualX::TranslateCoordinates(Window_t /*src*/, Window_t /*dest*/,
                                     Int_t /*src_x*/, Int_t /*src_y*/,
                                     Int_t &dest_x, Int_t &dest_y,
                                     Window_t &child)
{
   dest_x = dest_y = 0;
   child = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the location and the size of window "id"
///
/// \param [in] id     drawable identifier
/// \param [in] x, y   coordinates of the upper-left outer corner relative to the
///                    parent window's origin
/// \param [in] w, h   the inside size of the window, not including the border

void TVirtualX::GetWindowSize(Drawable_t /*id*/, Int_t &x, Int_t &y,
                              UInt_t &w, UInt_t &h)
{
   x = y = 0;
   w = h = 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills the region closed by the specified path. The path is closed
/// automatically if the last point in the list does not coincide with the
/// first point.
///
/// \param [in] id        window identifier
/// \param [in] gc        graphics context
/// \param [in] *points   specifies an array of points
/// \param [in] npnt      specifies the number of points in the array
///
/// GC components in use: function, plane-mask, fill-style, fill-rule,
/// subwindow-mode, clip-x-origin, clip-y-origin, and clip-mask.  GC
/// mode-dependent components: foreground, background, tile, stipple,
/// tile-stipple-x-origin, and tile-stipple-y-origin.
/// (see also the GCValues_t structure)

void TVirtualX::FillPolygon(Window_t /*id*/, GContext_t /*gc*/, Point_t *
                            /*points*/, Int_t /*npnt*/) {
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the root window the pointer is logically on and the pointer
/// coordinates relative to the root window's origin.
///
/// \param [in] id               specifies the window
/// \param [in] rootw            the root window that the pointer is in
/// \param [in] childw           the child window that the pointer is located in, if any
/// \param [in] root_x, root_y   the pointer coordinates relative to the root window's
///                              origin
/// \param [in] win_x, win_y     the pointer coordinates relative to the specified
///                              window "id"
/// \param [in] mask             the current state of the modifier keys and pointer
///                              buttons

void TVirtualX::QueryPointer(Window_t /*id*/, Window_t &rootw, Window_t &childw,
                             Int_t &root_x, Int_t &root_y, Int_t &win_x,
                             Int_t &win_y, UInt_t &mask)
{
   rootw = childw = kNone;
   root_x = root_y = win_x = win_y = 0;
   mask = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the foreground color for the specified GC (shortcut for ChangeGC
/// with only foreground mask set).
///
/// \param [in] gc           specifies the GC
/// \param [in] foreground   the foreground you want to set
///
/// (see also the GCValues_t structure)

void TVirtualX::SetForeground(GContext_t /*gc*/, ULong_t /*foreground*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Sets clipping rectangles in graphics context. [x,y] specify the origin
/// of the rectangles. "recs" specifies an array of rectangles that define
/// the clipping mask and "n" is the number of rectangles.
/// (see also the GCValues_t structure)

void TVirtualX::SetClipRectangles(GContext_t /*gc*/, Int_t /*x*/, Int_t /*y*/,
                                  Rectangle_t * /*recs*/, Int_t /*n*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Flushes (mode = 0, default) or synchronizes (mode = 1) X output buffer.
/// Flush flushes output buffer. Sync flushes buffer and waits till all
/// requests have been processed by X server.

void TVirtualX::Update(Int_t /*mode = 0*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a new empty region.

Region_t TVirtualX::CreateRegion()
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destroys the region "reg".

void TVirtualX::DestroyRegion(Region_t /*reg*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Updates the destination region from a union of the specified rectangle
/// and the specified source region.
///
/// \param [in] rect   specifies the rectangle
/// \param [in] src    specifies the source region to be used
/// \param [in] dest   returns the destination region

void TVirtualX::UnionRectWithRegion(Rectangle_t * /*rect*/, Region_t /*src*/,
                                    Region_t /*dest*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a region for the polygon defined by the points array.
///
/// \param [in] points    specifies an array of points
/// \param [in] np        specifies the number of points in the polygon
/// \param [in] winding   specifies the winding-rule is set (kTRUE) or not(kFALSE)

Region_t TVirtualX::PolygonRegion(Point_t * /*points*/, Int_t /*np*/,
                                  Bool_t /*winding*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Computes the union of two regions.
///
/// \param [in] rega, regb   specify the two regions with which you want to perform
///                 the computation
/// \param [in] result       returns the result of the computation

void TVirtualX::UnionRegion(Region_t /*rega*/, Region_t /*regb*/,
                            Region_t /*result*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Computes the intersection of two regions.
///
/// \param [in] rega, regb   specify the two regions with which you want to perform
///                 the computation
/// \param [in] result       returns the result of the computation

void TVirtualX::IntersectRegion(Region_t /*rega*/, Region_t /*regb*/,
                                Region_t /*result*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Subtracts regb from rega and stores the results in result.

void TVirtualX::SubtractRegion(Region_t /*rega*/, Region_t /*regb*/,
                               Region_t /*result*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates the difference between the union and intersection of
/// two regions.
///
/// \param [in] rega, regb   specify the two regions with which you want to perform
///                 the computation
/// \param [in] result       returns the result of the computation

void TVirtualX::XorRegion(Region_t /*rega*/, Region_t /*regb*/,
                          Region_t /*result*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if the region reg is empty.

Bool_t  TVirtualX::EmptyRegion(Region_t /*reg*/)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if the point [x, y] is contained in the region reg.

Bool_t  TVirtualX::PointInRegion(Int_t /*x*/, Int_t /*y*/, Region_t /*reg*/)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if the two regions have the same offset, size, and shape.

Bool_t  TVirtualX::EqualRegion(Region_t /*rega*/, Region_t /*regb*/)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns smallest enclosing rectangle.

void TVirtualX::GetRegionBox(Region_t /*reg*/, Rectangle_t * /*rect*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns list of font names matching fontname regexp, like "-*-times-*".
/// The pattern string can contain any characters, but each asterisk (*)
/// is a wildcard for any number of characters, and each question mark (?)
/// is a wildcard for a single character. If the pattern string is not in
/// the Host Portable Character Encoding, the result is implementation
/// dependent. Use of uppercase or lowercase does not matter. Each returned
/// string is null-terminated.
///
/// \param [in] fontname   specifies the null-terminated pattern string that can
///               contain wildcard characters
/// \param [in] max        specifies the maximum number of names to be returned
/// \param [in] count      returns the actual number of font names

char **TVirtualX::ListFonts(const char * /*fontname*/, Int_t /*max*/, Int_t & count)
{
   count=0;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Frees the specified the array of strings "fontlist".

void TVirtualX::FreeFontNames(char ** /*fontlist*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Allocates the memory needed for an drawable.
///
/// \param [in] width    the width of the image, in pixels
/// \param [in] height   the height of the image, in pixels

Drawable_t TVirtualX::CreateImage(UInt_t /*width*/, UInt_t /*height*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the width and height of the image id

void TVirtualX::GetImageSize(Drawable_t /*id*/, UInt_t &/*width*/,
                             UInt_t &/*height*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Overwrites the pixel in the image with the specified pixel value.
/// The image must contain the x and y coordinates.
///
/// \param [in] id      specifies the image
/// \param [in] x, y    coordinates
/// \param [in] pixel   the new pixel value

void TVirtualX::PutPixel(Drawable_t /*id*/, Int_t /*x*/, Int_t /*y*/,
                         ULong_t /*pixel*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Combines an image with a rectangle of the specified drawable. The
/// section of the image defined by the x, y, width, and height arguments
/// is drawn on the specified part of the drawable.
///
/// \param [in] id     the drawable
/// \param [in] gc     the GC
/// \param [in] img    the image you want combined with the rectangle
/// \param [in] dx     the offset in X from the left edge of the image
/// \param [in] dy     the offset in Y from the top edge of the image
/// \param [in] x, y   coordinates, which are relative to the origin of the
///           drawable and are the coordinates of the subimage
/// \param [in] w, h   the width and height of the subimage, which define the
///           rectangle dimensions
///
/// GC components in use: function, plane-mask, subwindow-mode,
/// clip-x-origin, clip-y-origin, and clip-mask.
/// GC mode-dependent components: foreground and background.
/// (see also the GCValues_t structure)

void TVirtualX::PutImage(Drawable_t /*id*/, GContext_t /*gc*/,
                         Drawable_t /*img*/, Int_t /*dx*/, Int_t /*dy*/,
                         Int_t /*x*/, Int_t /*y*/, UInt_t /*w*/, UInt_t /*h*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Deallocates the memory associated with the image img

void TVirtualX::DeleteImage(Drawable_t /*img*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// pointer to the current internal window used in canvas graphics

Window_t TVirtualX::GetCurrentWindow() const
{
   return (Window_t)0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns an array of pixels created from a part of drawable (defined by x, y, w, h)
/// in format:
///
/// ~~~ {.cpp}
/// b1, g1, r1, 0,  b2, g2, r2, 0 ... bn, gn, rn, 0 ..
/// ~~~
///
/// Pixels are numbered from left to right and from top to bottom.
/// By default all pixels from the whole drawable are returned.
///
/// Note that return array is 32-bit aligned

unsigned char *TVirtualX::GetColorBits(Drawable_t /*wid*/, Int_t /*x*/, Int_t /*y*/,
                                       UInt_t /*w*/, UInt_t /*h*/)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// create pixmap from RGB data. RGB data is in format:
///
/// ~~~ {.cpp}
/// b1, g1, r1, 0,  b2, g2, r2, 0 ... bn, gn, rn, 0 ..
/// ~~~
///
/// Pixels are numbered from left to right and from top to bottom.
/// Note that data must be 32-bit aligned

Pixmap_t TVirtualX::CreatePixmapFromData(unsigned char * /*bits*/, UInt_t /*width*/,
                                       UInt_t /*height*/)
{
   return (Pixmap_t)0;
}

////////////////////////////////////////////////////////////////////////////////
/// The Non-rectangular Window Shape Extension adds non-rectangular
/// windows to the System.
/// This allows for making shaped (partially transparent) windows

void TVirtualX::ShapeCombineMask(Window_t, Int_t, Int_t, Pixmap_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the width of the screen in millimeters.

UInt_t TVirtualX::ScreenWidthMM() const
{
   return 400;
}

////////////////////////////////////////////////////////////////////////////////
/// Deletes the specified property only if the property was defined on the
/// specified window and causes the X server to generate a PropertyNotify
/// event on the window unless the property does not exist.

void TVirtualX::DeleteProperty(Window_t, Atom_t&)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the actual type of the property; the actual format of the property;
/// the number of 8-bit, 16-bit, or 32-bit items transferred; the number of
/// bytes remaining to be read in the property; and a pointer to the data
/// actually returned.

Int_t TVirtualX::GetProperty(Window_t, Atom_t, Long_t, Long_t, Bool_t, Atom_t,
                             Atom_t*, Int_t*, ULong_t*, ULong_t*, unsigned char**)
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Changes the specified dynamic parameters if the pointer is actively
/// grabbed by the client and if the specified time is no earlier than the
/// last-pointer-grab time and no later than the current X server time.

void TVirtualX::ChangeActivePointerGrab(Window_t, UInt_t, Cursor_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Requests that the specified selection be converted to the specified
/// target type.

void TVirtualX::ConvertSelection(Window_t, Atom_t&, Atom_t&, Atom_t&, Time_t&)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Changes the owner and last-change time for the specified selection.

Bool_t TVirtualX::SetSelectionOwner(Window_t, Atom_t&)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Alters the property for the specified window and causes the X server
/// to generate a PropertyNotify event on that window.

void TVirtualX::ChangeProperties(Window_t, Atom_t, Atom_t, Int_t, UChar_t *, Int_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Add XdndAware property and the list of drag and drop types to the
/// Window win.

void TVirtualX::SetDNDAware(Window_t, Atom_t *)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Add the list of drag and drop types to the Window win.

void TVirtualX::SetTypeList(Window_t, Atom_t, Atom_t *)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively search in the children of Window for a Window which is at
/// location x, y and is DND aware, with a maximum depth of maxd.

Window_t TVirtualX::FindRWindow(Window_t, Window_t, Window_t, int, int, int)
{
   return kNone;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the Window is DND aware, and knows any of the DND formats
/// passed in argument.

Bool_t TVirtualX::IsDNDAware(Window_t, Atom_t *)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Start a modal session for a dialog window.

void TVirtualX::BeginModalSessionFor(Window_t)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns 1 if window system server supports extension given by the
/// argument, returns 0 in case extension is not supported and returns -1
/// in case of error (like server not initialized).

Int_t TVirtualX::SupportsExtension(const char *) const
{
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Map the XftFont with the Graphics Context using it.

void TVirtualX::MapGCFont(GContext_t, FontStruct_t)
{
}


