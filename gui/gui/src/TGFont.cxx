// @(#)root/gui:$Id$
// Author: Fons Rademakers   20/5/2003

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
/**************************************************************************

    This source is based on Xclass95, a Win95-looking GUI toolkit.
    Copyright (C) 1996, 1997 David Barth, Ricky Ralston, Hector Peraza.

    Xclass95 is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

**************************************************************************/


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGFont and TGFontPool                                                //
//                                                                      //
// Encapsulate fonts used in the GUI system.                            //
// TGFontPool provides a pool of fonts.                                 //
// TGTextLayout is used to keep track of string  measurement            //
// information when  using the text layout facilities.                  //
// It can be displayed with respect to any origin.                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGFont.h"
#include "TGClient.h"
#include "THashTable.h"
#include "TVirtualX.h"
#include "TObjString.h"
#include "TGWidget.h"
#include "TROOT.h"
#include "TError.h"
#include "TMath.h"
#include "strlcpy.h"
#include "snprintf.h"

#include <cerrno>
#include <cstdlib>
#include <climits>
#include <iostream>


ClassImp(TGFont);
ClassImp(TGFontPool);
ClassImp(TGTextLayout);

#define FONT_FAMILY     0
#define FONT_SIZE       1
#define FONT_WEIGHT     2
#define FONT_SLANT      3
#define FONT_UNDERLINE  4
#define FONT_OVERSTRIKE 5
#define FONT_NUMFIELDS  6

// The following defines specify the meaning of the fields in a fully
// qualified XLFD.

#define XLFD_FOUNDRY        0
#define XLFD_FAMILY         1
#define XLFD_WEIGHT         2
#define XLFD_SLANT          3
#define XLFD_SETWIDTH       4
#define XLFD_ADD_STYLE      5
#define XLFD_PIXEL_SIZE     6
#define XLFD_POINT_SIZE     7
#define XLFD_RESOLUTION_X   8
#define XLFD_RESOLUTION_Y   9
#define XLFD_SPACING        10
#define XLFD_AVERAGE_WIDTH  11
#define XLFD_REGISTRY       12
#define XLFD_ENCODING       13
#define XLFD_NUMFIELDS      14   // Number of fields in XLFD.


// A LayoutChunk_t represents a contiguous range of text that can be measured
// and displayed by low-level text calls. In general, chunks will be
// delimited by newlines and tabs. Low-level, platform-specific things
// like kerning and non-integer character widths may occur between the
// characters in a single chunk, but not between characters in different
// chunks.

struct LayoutChunk_t {

   const char *fStart;     // Pointer to simple string to be displayed.
                           // This is a pointer into the TGTextLayout's
                           // string.
   Int_t fNumChars;        // The number of characters in this chunk.
   Int_t fNumDisplayChars; // The number of characters to display when
                           // this chunk is displayed. Can be less than
                           // numChars if extra space characters were
                           // absorbed by the end of the chunk. This
                           // will be < 0 if this is a chunk that is
                           // holding a tab or newline.
   Int_t fX;               // The x origin and
   Int_t fY;               // the y origin of the first character in this
                           // chunk with respect to the upper-left hand
                           // corner of the TGTextLayout.
   Int_t fTotalWidth;      // Width in pixels of this chunk. Used
                           // when hit testing the invisible spaces at
                           // the end of a chunk.
   Int_t fDisplayWidth;    // Width in pixels of the displayable
                           // characters in this chunk. Can be less than
                           // width if extra space characters were
                           // absorbed by the end of the chunk.
};


// The following structure is used to return attributes when parsing an
// XLFD. The extra information is used to find the closest matching font.

struct XLFDAttributes_t {
   FontAttributes_t fFA; // Standard set of font attributes.
   const char *fFoundry; // The foundry of the font.
   Int_t fSlant;         // The tristate value for the slant
   Int_t fSetwidth;      // The proportionate width
   Int_t fCharset;       // The character set encoding (the glyph family).
   Int_t fEncoding;      // Variations within a charset for the glyphs above character 127.

   XLFDAttributes_t() :  // default constructor
      fFA(),
      fFoundry(0),
      fSlant(0),
      fSetwidth(0),
      fCharset(0),
      fEncoding(0) { }
};


// The following data structure is used to keep track of the font attributes
// for each named font that has been defined. The named font is only deleted
// when the last reference to it goes away.

class TNamedFont : public TObjString, public TRefCnt {
public:
   Int_t            fDeletePending; // Non-zero if font should be deleted when last reference goes away.
   FontAttributes_t fFA;            // Desired attributes for named font.
};

// enums
enum EFontSpacing { kFontProportional = 0,
                    kFontFixed = 1,
                    kFontMono = 1,
                    kFontCharcell = 2 };

enum EFontSetWidth { kFontSWNormal = 0,
                     kFontSWCondence = 1,
                     kFontSWExpand = 2,
                     kFontSWUnknown = 3 };

enum EFontCharset { kFontCSNormal = 0,
                    kFontCSSymbol = 1,
                    kFontCSOther = 2 };


// Possible values for entries in the "types" field in a TGFont structure,
// which classifies the types of all characters in the given font. This
// information is used when measuring and displaying characters.
//
// kCharNormal:         Standard character.
// kCharReplace:        This character doesn't print: instead of displaying
//                      character, display a replacement sequence like "\n"
//                      (for those characters where ANSI C defines such a
//                      sequence) or a sequence of the form "\xdd" where dd
//                      is the hex equivalent of the character.
// kCharSkip:           Don't display anything for this character. This is
//                      only used where the font doesn't contain all the
//                      characters needed to generate replacement sequences.
enum ECharType { kCharNormal, kCharReplace, kCharSkip };


// The following structures are used as two-way maps between the values for
// the fields in the FontAttributes_t structure and the strings used when
// parsing both option-value format and style-list format font name strings.

struct FontStateMap_t { Int_t fNumKey; const char *fStrKey; };

static const FontStateMap_t gWeightMap[] = {
   { kFontWeightNormal,  "normal" },
   { kFontWeightBold,    "bold"   },
   { kFontWeightUnknown, 0        }
};

static const FontStateMap_t gSlantMap[] = {
   { kFontSlantRoman,   "roman"  },
   { kFontSlantItalic,  "italic" },
   { kFontSlantUnknown, 0        }
};

static const FontStateMap_t gUnderlineMap[] = {
   { 1, "underline" },
   { 0, 0           }
};

static const FontStateMap_t gOverstrikeMap[] = {
   { 1, "overstrike" },
   { 0, 0            }
};

// The following structures are used when parsing XLFD's into a set of
// FontAttributes_t.

static const FontStateMap_t gXlfdgWeightMap[] = {
   { kFontWeightNormal, "normal"   },
   { kFontWeightNormal, "medium"   },
   { kFontWeightNormal, "book"     },
   { kFontWeightNormal, "light"    },
   { kFontWeightBold,   "bold"     },
   { kFontWeightBold,   "demi"     },
   { kFontWeightBold,   "demibold" },
   { kFontWeightNormal,  0         }  // Assume anything else is "normal".
};

static const FontStateMap_t gXlfdSlantMap[] = {
   { kFontSlantRoman,   "r"  },
   { kFontSlantItalic,  "i"  },
   { kFontSlantOblique, "o"  },
   { kFontSlantRoman,   0    }  // Assume anything else is "roman".
};

static const FontStateMap_t gXlfdSetwidthMap[] = {
   { kFontSWNormal,   "normal"        },
   { kFontSWCondence, "narrow"        },
   { kFontSWCondence, "semicondensed" },
   { kFontSWCondence, "condensed"     },
   { kFontSWUnknown,  0               }
};

static const FontStateMap_t gXlfdCharsetMap[] = {
   { kFontCSNormal, "iso8859" },
   { kFontCSSymbol, "adobe"   },
   { kFontCSSymbol, "sun"     },
   { kFontCSOther,  0         }
};


// Characters used when displaying control sequences.

static char gHexChars[] = "0123456789abcdefxtnvr\\";


// The following table maps some control characters to sequences like '\n'
// rather than '\x10'. A zero entry in the table means no such mapping
// exists, and the table only maps characters less than 0x10.

static char gMapChars[] = {
   0, 0, 0, 0, 0, 0, 0, 'a', 'b', 't', 'n', 'v', 'f', 'r', 0
};

static int GetControlCharSubst(int c, char buf[4]);


////////////////////////////////////////////////////////////////////////////////
/// Delete font.

TGFont::~TGFont()
{
   if (fFontStruct) {
      gVirtualX->DeleteFont(fFontStruct);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get font metrics.

void TGFont::GetFontMetrics(FontMetrics_t *m) const
{
   if (!m) {
      Error("GetFontMetrics", "argument may not be 0");
      return;
   }

   *m = fFM;
   m->fLinespace = fFM.fAscent + fFM.fDescent;
}

////////////////////////////////////////////////////////////////////////////////
/// Not inline due to a bug in g++ 2.96 20000731 (Red Hat Linux 7.0)

FontStruct_t TGFont::operator()() const
{
   return fFontStruct;
}

////////////////////////////////////////////////////////////////////////////////
/// Print font info.

void TGFont::Print(Option_t *option) const
{
   TString opt = option;

   if ((opt == "full") && fNamedHash) {
      Printf("TGFont: %s, %s, ref cnt = %u",
              fNamedHash->GetName(),
              fFM.fFixed ? "fixed" : "prop", References());
   } else {
      Printf("TGFont: %s, %s, ref cnt = %u", fName.Data(),
              fFM.fFixed ? "fixed" : "prop", References());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the name of the corresponding Postscript font for this TGFont.
///
/// The return value is the pointsize of the TGFont. The name of the
/// Postscript font is appended to ds.
///
/// If the font does not exist on the printer, the print job will fail at
/// print time. Given a "reasonable" Postscript printer, the following
/// TGFont font families should print correctly:
///
///     Avant Garde, Arial, Bookman, Courier, Courier New, Geneva,
///     Helvetica, Monaco, New Century Schoolbook, New York,
///     Palatino, Symbol, Times, Times New Roman, Zapf Chancery,
///     and Zapf Dingbats.
///
/// Any other TGFont font families may not print correctly because the
/// computed Postscript font name may be incorrect.
///
/// dst -- Pointer to an initialized TString object to which the name of the
///        Postscript font that corresponds to the font will be appended.

Int_t TGFont::PostscriptFontName(TString *dst) const
{
   const char *family;
   TString weightString;
   TString slantString;
   char *src, *dest;
   Int_t upper, len;

   len = dst->Length();

   // Convert the case-insensitive TGFont family name to the
   // case-sensitive Postscript family name. Take out any spaces and
   // capitalize the first letter of each word.

   family = fFA.fFamily;
   if (strncasecmp(family, "itc ", 4) == 0) {
      family = family + 4;
   }
   if ((strcasecmp(family, "Arial") == 0)
       || (strcasecmp(family, "Geneva") == 0)) {
      family = "Helvetica";
   } else if ((strcasecmp(family, "Times New Roman") == 0)
              || (strcasecmp(family, "New York") == 0)) {
      family = "Times";
   } else if ((strcasecmp(family, "Courier New") == 0)
              || (strcasecmp(family, "Monaco") == 0)) {
      family = "Courier";
   } else if (strcasecmp(family, "AvantGarde") == 0) {
      family = "AvantGarde";
   } else if (strcasecmp(family, "ZapfChancery") == 0) {
      family = "ZapfChancery";
   } else if (strcasecmp(family, "ZapfDingbats") == 0) {
      family = "ZapfDingbats";
   } else {

      // Inline, capitalize the first letter of each word, lowercase the
      // rest of the letters in each word, and then take out the spaces
      // between the words. This may make the TString shorter, which is
      // safe to do.

      dst->Append(family);

      src = dest = (char*)dst->Data() + len;
      upper = 1;
      for (; *src != '\0'; src++, dest++) {
         while (isspace(UChar_t(*src))) {
            src++;
            upper = 1;
         }
         *dest = *src;
         if ((upper != 0) && (islower(UChar_t(*src)))) {
            *dest = toupper(UChar_t(*src));
         }
         upper = 0;
      }
      *dest = '\0';
      //dst->SetLength(dest - dst->GetString()); // dst->ResetLength(); may be better
      family = (char *) dst->Data() + len;
   }
   if (family != (char *) dst->Data() + len) {
      dst->Append(family);
      family = (char *) dst->Data() + len;
   }
   if (strcasecmp(family, "NewCenturySchoolbook") == 0) {
//      dst->SetLength(len);
      dst->Append("NewCenturySchlbk");
      family = (char *) dst->Data() + len;
   }

   // Get the string to use for the weight.

   weightString = "";
   if (fFA.fWeight == kFontWeightNormal) {
      if (strcmp(family, "Bookman") == 0) {
         weightString = "Light";
      } else if (strcmp(family, "AvantGarde") == 0) {
         weightString = "Book";
      } else if (strcmp(family, "ZapfChancery") == 0) {
         weightString = "Medium";
      }
   } else {
      if ((strcmp(family, "Bookman") == 0)
           || (strcmp(family, "AvantGarde") == 0)) {
         weightString = "Demi";
      } else {
         weightString = "Bold";
      }
   }

   // Get the string to use for the slant.

   slantString = "";
   if (fFA.fSlant == kFontSlantRoman) {
      ;
   } else {
      if ((strcmp(family, "Helvetica") == 0)
           || (strcmp(family, "Courier") == 0)
           || (strcmp(family, "AvantGarde") == 0)) {
         slantString = "Oblique";
      } else {
         slantString = "Italic";
      }
   }

   // The string "Roman" needs to be added to some fonts that are not bold
   // and not italic.

   if ((slantString.IsNull()) && (weightString.IsNull())) {
      if ((strcmp(family, "Times") == 0)
           || (strcmp(family, "NewCenturySchlbk") == 0)
           || (strcmp(family, "Palatino") == 0)) {
         dst->Append("-Roman");
      }
   } else {
      dst->Append("-");
      if (!weightString.IsNull()) dst->Append(weightString);
      if (!slantString.IsNull()) dst->Append(slantString);
   }

   return fFA.fPointsize;
}

////////////////////////////////////////////////////////////////////////////////
/// Determine the number of characters from the string that will fit in the
/// given horizontal span. The measurement is done under the assumption that
/// DrawChars() will be used to actually display the characters.
///
/// The return value is the number of characters from source that fit into
/// the span that extends from 0 to maxLength. *length is filled with the
/// x-coordinate of the right edge of the last character that did fit.
///
/// source    -- Characters to be displayed. Need not be '\0' terminated.
/// numChars  -- Maximum number of characters to consider from source string.
/// maxLength -- If > 0, maxLength specifies the longest permissible line
///              length; don't consider any character that would cross this
///              x-position. If <= 0, then line length is unbounded and the
///              flags argument is ignored.
/// flags     -- Various flag bits OR-ed together:
///              TEXT_PARTIAL_OK means include the last char which only
///              partially fit on this line.
///              TEXT_WHOLE_WORDS means stop on a word boundary, if possible.
///              TEXT_AT_LEAST_ONE means return at least one character even
///              if no characters fit.
/// *length   -- Filled with x-location just after the terminating character.

Int_t TGFont::MeasureChars(const char *source, Int_t numChars, Int_t maxLength,
                          Int_t flags, Int_t *length) const
{
   const char *p;    // Current character.
   const char *term; // Pointer to most recent character that may legally be a terminating character.
   Int_t termX;      // X-position just after term.
   Int_t curX;       // X-position corresponding to p.
   Int_t newX;       // X-position corresponding to p+1.
   Int_t c, sawNonSpace;

   if (!numChars) {
      *length = 0;
      return 0;
   }
   if (maxLength <= 0) {
      maxLength = INT_MAX;
   }
   newX = curX = termX = 0;
   p = term = source;
   sawNonSpace = !isspace(UChar_t(*p));

   // Scan the input string one character at a time, calculating width.

   for (c = UChar_t(*p);;) {
      newX += fWidths[c];
      if (newX > maxLength) {
         break;
      }
      curX = newX;
      numChars--;
      p++;
      if (!numChars) {
         term = p;
         termX = curX;
         break;
      }
      c = UChar_t(*p);
      if (isspace(c)) {
         if (sawNonSpace) {
            term = p;
            termX = curX;
            sawNonSpace = 0;
         }
      } else {
         sawNonSpace = 1;
      }
   }

   // P points to the first character that doesn't fit in the desired
   // span. Use the flags to figure out what to return.

   if ((flags & kTextPartialOK) && (numChars > 0) && (curX < maxLength)) {

      // Include the first character that didn't quite fit in the desired
      // span. The width returned will include the width of that extra
      // character.

      numChars--;
      curX = newX;
      p++;
   }
   if ((flags & kTextAtLeastOne) && (term == source) && (numChars > 0)) {
      term = p;
      termX = curX;
      if (term == source) {
         term++;
         termX = newX;
      }
   } else if ((numChars == 0) || !(flags & kTextWholeWords)) {
      term = p;
      termX = curX;
   }
   *length = termX;

   return term - source;
}

////////////////////////////////////////////////////////////////////////////////
/// A wrapper function for the more complicated interface of MeasureChars.
/// Computes how much space the given simple string needs.
///
/// The return value is the width (in pixels) of the given string.
///
/// string   -- String whose width will be computed.
/// numChars -- Number of characters to consider from string, or < 0 for
///             strlen().

Int_t TGFont::TextWidth(const char *string, Int_t numChars) const
{
   Int_t width;

   if (numChars < 0) {
      numChars = strlen(string);
   }
   MeasureChars(string, numChars, 0, 0, &width);

   return width;
}

////////////////////////////////////////////////////////////////////////////////
/// Return text widht in pixels

Int_t TGFont::XTextWidth(const char *string, Int_t numChars) const
{
   int width;

   if (numChars < 0) {
      numChars = strlen(string);
   }
   width = gVirtualX->TextWidth(fFontStruct, string, numChars);

   return width;
}

////////////////////////////////////////////////////////////////////////////////
/// This procedure draws an underline for a given range of characters in a
/// given string. It doesn't draw the characters (which are assumed to have
/// been displayed previously); it just draws the underline. This procedure
/// would mainly be used to quickly underline a few characters without having
/// to construct an underlined font. To produce properly underlined text, the
/// appropriate underlined font should be constructed and used.
///
/// dst       -- Window or pixmap in which to draw.
/// gc        -- Graphics context for actually drawing line.
/// string    -- String containing characters to be underlined or overstruck.
/// x, y      -- Coordinates at which first character of string is drawn.
/// firstChar -- Index of first character.
/// lastChar  -- Index of one after the last character.

void TGFont::UnderlineChars(Drawable_t dst, GContext_t gc,
                            const char *string, Int_t x, Int_t y,
                            Int_t firstChar, Int_t lastChar) const
{
   Int_t startX, endX;

   MeasureChars(string, firstChar, 0, 0, &startX);
   MeasureChars(string, lastChar, 0, 0, &endX);

   gVirtualX->FillRectangle(dst, gc, x + startX, y + fUnderlinePos,
                            (UInt_t) (endX - startX),
                            (UInt_t) fUnderlineHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes the amount of screen space needed to display a multi-line,
/// justified string of text. Records all the measurements that were done
/// to determine to size and positioning of the individual lines of text;
/// this information can be used by the TGTextLayout::DrawText() procedure
/// to display the text quickly (without remeasuring it).
///
/// This procedure is useful for simple widgets that want to display
/// single-font, multi-line text and want TGFont to handle the details.
///
/// The return value is a TGTextLayout token that holds the measurement
/// information for the given string. The token is only valid for the given
/// string. If the string is freed, the token is no longer valid and must
/// also be deleted.
///
/// The dimensions of the screen area needed to display the text are stored
/// in *width and *height.
///
/// string     -- String whose dimensions are to be computed.
/// numChars   -- Number of characters to consider from string, or < 0 for
///               strlen().
/// wrapLength -- Longest permissible line length, in pixels. <= 0 means no
///               automatic wrapping: just let lines get as long as needed.
/// justify    -- How to justify lines.
/// flags      -- Flag bits OR-ed together. kTextIgnoreTabs means that tab
///               characters should not be expanded. kTextIgnoreNewlines
///               means that newline characters should not cause a line break.
/// width      -- Filled with width of string.
/// height     -- Filled with height of string.

TGTextLayout *TGFont::ComputeTextLayout(const char *string, Int_t numChars,
                                        Int_t wrapLength, Int_t justify, Int_t flags,
                                        UInt_t *width, UInt_t *height) const
{
   const char *start, *end, *special;
   Int_t n, y=0, charsThisChunk, maxChunks;
   Int_t baseline, h, curX, newX, maxWidth;
   TGTextLayout *layout;
   LayoutChunk_t *chunk;

#define MAX_LINES 50
   Int_t staticLineLengths[MAX_LINES];
   Int_t *lineLengths;
   Int_t maxLines, curLine, layoutHeight;

   lineLengths = staticLineLengths;
   maxLines = MAX_LINES;

   h = fFM.fAscent + fFM.fDescent;

   if (numChars < 0) {
      numChars = strlen(string);
   }
   maxChunks = 0;

   layout = new TGTextLayout;
   layout->fFont = this;
   layout->fString = string;
   layout->fNumChunks = 0;
   layout->fChunks = 0;

   baseline = fFM.fAscent;
   maxWidth = 0;

   // Divide the string up into simple strings and measure each string.

   curX = 0;

   end = string + numChars;
   special = string;

   flags &= kTextIgnoreTabs | kTextIgnoreNewlines;
   flags |= kTextWholeWords | kTextAtLeastOne;
   curLine = 0;

   for (start = string; start < end;) {
      if (start >= special) {
         // Find the next special character in the string.

         for (special = start; special < end; special++) {
            if (!(flags & kTextIgnoreNewlines)) {
               if ((*special == '\n') || (*special == '\r')) {
                  break;
               }
            }
            if (!(flags & kTextIgnoreTabs)) {
               if (*special == '\t') {
                  break;
               }
            }
         }
      }

      // Special points at the next special character (or the end of the
      // string). Process characters between start and special.

      chunk = 0;
      if (start < special) {
         charsThisChunk = MeasureChars(start, special - start,
                                       wrapLength - curX, flags, &newX);
         newX += curX;
         flags &= ~kTextAtLeastOne;
         if (charsThisChunk > 0) {
            chunk = NewChunk(layout, &maxChunks, start,
                             charsThisChunk, curX, newX, baseline);

            start += charsThisChunk;
            curX = newX;
         }
      }
      if ((start == special) && (special < end)) {
         // Handle the special character.
         LayoutChunk_t *newchunk = 0;

         chunk = 0;
         if (*special == '\t') {
            newX = curX + fTabWidth;
            newX -= newX % fTabWidth;
            newchunk = NewChunk(layout, &maxChunks, start, 1, curX, newX, baseline);
            if (newchunk) newchunk->fNumDisplayChars = -1;
            start++;
            if ((start < end) && ((wrapLength <= 0) || (newX <= wrapLength))) {

               // More chars can still fit on this line.

               curX = newX;
               flags &= ~kTextAtLeastOne;
               continue;
            }
         } else {
            newchunk = NewChunk(layout, &maxChunks, start, 1, curX, 1000000000, baseline);
            if (newchunk) newchunk->fNumDisplayChars = -1;
            start++;
            goto wrapLine;
         }
      }

      // No more characters are going to go on this line, either because
      // no more characters can fit or there are no more characters left.
      // Consume all extra spaces at end of line.

      while ((start < end) && isspace(UChar_t(*start))) {
         if (!(flags & kTextIgnoreNewlines)) {
            if ((*start == '\n') || (*start == '\r')) {
               break;
            }
         }
         if (!(flags & kTextIgnoreTabs)) {
            if (*start == '\t') {
               break;
            }
         }
         start++;
      }
      if (chunk) {
         // Append all the extra spaces on this line to the end of the
         // last text chunk.

         charsThisChunk = start - (chunk->fStart + chunk->fNumChars);
         if (charsThisChunk > 0) {
            chunk->fNumChars += MeasureChars(chunk->fStart + chunk->fNumChars,
                                             charsThisChunk, 0, 0, &chunk->fTotalWidth);
            chunk->fTotalWidth += curX;
         }
      }
wrapLine:
      flags |= kTextAtLeastOne;

      // Save current line length, then move current position to start of
      // next line.

      if (curX > maxWidth) {
         maxWidth = curX;
      }

      // Remember width of this line, so that all chunks on this line
      // can be centered or right justified, if necessary.

      if (curLine >= maxLines) {
         int *newLengths;

         newLengths = new int[2 * maxLines];
         memcpy((void *) newLengths, lineLengths, maxLines * sizeof (int));

         if (lineLengths != staticLineLengths) {
            delete[] lineLengths;
         }
         lineLengths = newLengths;
         maxLines *= 2;
      }
      lineLengths[curLine] = curX;
      curLine++;

      curX = 0;
      baseline += h;
   }

   // If last line ends with a newline, then we need to make a 0 width
   // chunk on the next line. Otherwise "Hello" and "Hello\n" are the
   // same height.

   if ((layout->fNumChunks > 0) && ((flags & kTextIgnoreNewlines) == 0)) {
      if (layout->fChunks[layout->fNumChunks - 1].fStart[0] == '\n') {
         chunk = NewChunk(layout, &maxChunks, start, 0, curX, 1000000000, baseline);
         chunk->fNumDisplayChars = -1;
         baseline += h;
      }
   }

   // Using maximum line length, shift all the chunks so that the lines are
   // all justified correctly.

   curLine = 0;
   chunk = layout->fChunks;
   if (chunk) y = chunk->fY;
   for (n = 0; n < layout->fNumChunks; n++) {
      int extra;

      if (chunk->fY != y) {
         curLine++;
         y = chunk->fY;
      }
      extra = maxWidth - lineLengths[curLine];
      if (justify == kTextCenterX) {
         chunk->fX += extra / 2;
      } else if (justify == kTextRight) {
         chunk->fX += extra;
      }
      ++chunk;
   }

   layout->fWidth = maxWidth;
   layoutHeight = baseline - fFM.fAscent;
   if (layout->fNumChunks == 0) {
      layoutHeight = h;

      // This fake chunk is used by the other procedures so that they can
      // pretend that there is a chunk with no chars in it, which makes
      // the coding simpler.

      layout->fNumChunks = 1;
      layout->fChunks = new LayoutChunk_t[1];
      layout->fChunks[0].fStart = string;
      layout->fChunks[0].fNumChars = 0;
      layout->fChunks[0].fNumDisplayChars = -1;
      layout->fChunks[0].fX = 0;
      layout->fChunks[0].fY = fFM.fAscent;
      layout->fChunks[0].fTotalWidth = 0;
      layout->fChunks[0].fDisplayWidth = 0;
   }
   if (width) {
      *width = layout->fWidth;
   }
   if (height) {
      *height = layoutHeight;
   }
   if (lineLengths != staticLineLengths) {
      delete[] lineLengths;
   }

   return layout;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGTextLayout::~TGTextLayout()
{
   if (fChunks) {
      delete[] fChunks;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Use the information in the TGTextLayout object to display a multi-line,
/// justified string of text.
///
/// This procedure is useful for simple widgets that need to display
/// single-font, multi-line text and want TGFont to handle the details.
///
/// dst       -- Window or pixmap in which to draw.
/// gc        -- Graphics context to use for drawing text.
/// x, y      -- Upper-left hand corner of rectangle in which to draw
///              (pixels).
/// firstChar -- The index of the first character to draw from the given
///              text item. 0 specfies the beginning.
/// lastChar  -- The index just after the last character to draw from the
///              given text item. A number < 0 means to draw all characters.

void TGTextLayout::DrawText(Drawable_t dst, GContext_t gc,
                            Int_t x, Int_t y, Int_t firstChar, Int_t lastChar) const
{
   Int_t i, numDisplayChars, drawX;
   LayoutChunk_t *chunk;

   if (lastChar < 0) lastChar = 100000000;
   chunk = fChunks;

   for (i = 0; i < fNumChunks; i++) {
      numDisplayChars = chunk->fNumDisplayChars;
      if ((numDisplayChars > 0) && (firstChar < numDisplayChars)) {
         if (firstChar <= 0) {
            drawX = 0;
            firstChar = 0;
         } else {
            fFont->MeasureChars(chunk->fStart, firstChar, 0, 0, &drawX);
         }
         if (lastChar < numDisplayChars) numDisplayChars = lastChar;
         fFont->DrawChars(dst, gc, chunk->fStart + firstChar, numDisplayChars - firstChar,
                           x + chunk->fX + drawX, y + chunk->fY);
      }
      firstChar -= chunk->fNumChars;
      lastChar -= chunk->fNumChars;

      if (lastChar <= 0) break;
      chunk++;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Use the information in the TGTextLayout object to display an underline
/// below an individual character. This procedure does not draw the text,
/// just the underline.
///
/// This procedure is useful for simple widgets that need to display
/// single-font, multi-line text with an individual character underlined
/// and want TGFont to handle the details. To display larger amounts of
/// underlined text, construct and use an underlined font.
///
/// dst       -- Window or pixmap in which to draw.
/// gc        -- Graphics context to use for drawing text.
/// x, y      -- Upper-left hand corner of rectangle in which to draw
///              (pixels).
/// underline -- Index of the single character to underline, or -1 for
///              no underline.

void TGTextLayout::UnderlineChar(Drawable_t dst, GContext_t gc,
                                Int_t x, Int_t y, Int_t underline) const
{
   int xx, yy, width, height;

   if ((CharBbox(underline, &xx, &yy, &width, &height) != 0)
      && (width != 0)) {
      gVirtualX->FillRectangle(dst, gc, x + xx,
                               y + yy + fFont->fFM.fAscent + fFont->fUnderlinePos,
                               (UInt_t) width, (UInt_t) fFont->fUnderlineHeight);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Use the information in the TGTextLayout token to determine the character
/// closest to the given point. The point must be specified with respect to
/// the upper-left hand corner of the text layout, which is considered to be
/// located at (0, 0).
///
/// Any point whose y-value is less that 0 will be considered closest to the
/// first character in the text layout; any point whose y-value is greater
/// than the height of the text layout will be considered closest to the last
/// character in the text layout.
///
/// Any point whose x-value is less than 0 will be considered closest to the
/// first character on that line; any point whose x-value is greater than the
/// width of the text layout will be considered closest to the last character
/// on that line.
///
/// The return value is the index of the character that was closest to the
/// point. Given a text layout with no characters, the value 0 will always
/// be returned, referring to a hypothetical zero-width placeholder character.

Int_t TGTextLayout::PointToChar(Int_t x, Int_t y) const
{
   LayoutChunk_t *chunk, *last;
   Int_t i, n, dummy, baseline, pos;

   if (y < 0) {
      // Point lies above any line in this layout. Return the index of
      // the first char.

      return 0;
   }

   // Find which line contains the point.

   last = chunk = fChunks;
   for (i = 0; i < fNumChunks; i++) {
      baseline = chunk->fY;
      if (y < baseline + fFont->fFM.fDescent) {
         if (x < chunk->fX) {
            // Point is to the left of all chunks on this line. Return
            // the index of the first character on this line.

            return (chunk->fStart - fString);
         }
         if (x >= fWidth) {

         // If point lies off right side of the text layout, return
         // the last char in the last chunk on this line. Without
         // this, it might return the index of the first char that
         // was located outside of the text layout.

            x = INT_MAX;
         }

         // Examine all chunks on this line to see which one contains
         // the specified point.

         last = chunk;
         while ((i < fNumChunks) && (chunk->fY == baseline)) {
            if (x < chunk->fX + chunk->fTotalWidth) {

               // Point falls on one of the characters in this chunk.

               if (chunk->fNumDisplayChars < 0) {

                  // This is a special chunk that encapsulates a single
                  // tab or newline char.

                  return (chunk->fStart - fString);
               }
               n = fFont->MeasureChars(chunk->fStart, chunk->fNumChars,
                                       x + 1 - chunk->fX, kTextPartialOK, &dummy);
               return ((chunk->fStart + n - 1) - fString);
            }
            last = chunk;
            chunk++;
            i++;
         }

         // Point is to the right of all chars in all the chunks on this
         // line. Return the index just past the last char in the last
         // chunk on this line.

         pos = (last->fStart + last->fNumChars) - fString;
         if (i < fNumChunks) pos--;
         return pos;
      }
      last = chunk;
      chunk++;
   }

   // Point lies below any line in this text layout. Return the index
   // just past the last char.

   return ((last->fStart + last->fNumChars) - fString);
}

////////////////////////////////////////////////////////////////////////////////
/// Use the information in the TGTextLayout token to return the bounding box
/// for the character specified by index.
///
/// The width of the bounding box is the advance width of the character, and
/// does not include and left- or right-bearing. Any character that extends
/// partially outside of the text layout is considered to be truncated at the
/// edge. Any character which is located completely outside of the text
/// layout is considered to be zero-width and pegged against the edge.
///
/// The height of the bounding box is the line height for this font,
/// extending from the top of the ascent to the bottom of the descent.
/// Information about the actual height of the individual letter is not
/// available.
///
/// A text layout that contains no characters is considered to contain a
/// single zero-width placeholder character.
///
/// The return value is 0 if the index did not specify a character in the
/// text layout, or non-zero otherwise. In that case, *bbox is filled with
/// the bounding box of the character.
///
/// layout -- Layout information, from a previous call to ComputeTextLayout().
/// index  -- The index of the character whose bbox is desired.
/// x, y   -- Filled with the upper-left hand corner, in pixels, of the
///           bounding box for the character specified by index, if non-NULL.
/// w, h   -- Filled with the width and height of the bounding box for the
///           character specified by index, if non-NULL.

Int_t TGTextLayout::CharBbox(Int_t index, Int_t *x, Int_t *y, Int_t *w, Int_t *h) const
{
   LayoutChunk_t *chunk;
   Int_t i, xx = 0, ww = 0;

   if (index < 0) {
      return 0;
   }

   chunk = fChunks;

   for (i = 0; i < fNumChunks; i++) {
      if (chunk->fNumDisplayChars < 0) {
         if (!index) {
            xx = chunk->fX;
            ww = chunk->fTotalWidth;
            goto check;
         }
      } else if (index < chunk->fNumChars) {
         if (x) {
            fFont->MeasureChars(chunk->fStart, index, 0, 0, &xx);
            xx += chunk->fX;
         }
         if (w) {
            fFont->MeasureChars(chunk->fStart + index, 1, 0, 0, &ww);
         }
         goto check;
      }
      index -= chunk->fNumChars;
      chunk++;
   }
   if (!index) {

      // Special case to get location just past last char in layout.

      chunk--;
      xx = chunk->fX + chunk->fTotalWidth;
      ww = 0;
   } else {
      return 0;
   }

   // Ensure that the bbox lies within the text layout. This forces all
   // chars that extend off the right edge of the text layout to have
   // truncated widths, and all chars that are completely off the right
   // edge of the text layout to peg to the edge and have 0 width.

check:
   if (y) {
      *y = chunk->fY - fFont->fFM.fAscent;
   }
   if (h) {
      *h = fFont->fFM.fAscent + fFont->fFM.fDescent;
   }
   if (xx > fWidth) {
      xx = fWidth;
   }
   if (x) {
      *x = xx;
   }
   if (w) {
      if (xx + ww > fWidth) {
         ww = fWidth - xx;
      }
      *w = ww;
   }
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Computes the distance in pixels from the given point to the given
/// text layout. Non-displaying space characters that occur at the end of
/// individual lines in the text layout are ignored for hit detection
/// purposes.
///
/// The return value is 0 if the point (x, y) is inside the text layout.
/// If the point isn't inside the text layout then the return value is the
/// distance in pixels from the point to the text item.
///
/// x, y -- Coordinates of point to check, with respect to the upper-left
///         corner of the text layout (in pixels).

Int_t TGTextLayout::DistanceToText(Int_t x, Int_t y) const
{
   Int_t i, x1, x2, y1, y2, xDiff, yDiff, dist, minDist, ascent, descent;
   LayoutChunk_t *chunk;

   ascent = fFont->fFM.fAscent;
   descent = fFont->fFM.fDescent;

   minDist = 0;
   chunk = fChunks;
   for (i = 0; i < fNumChunks; i++) {
      if (chunk->fStart[0] == '\n') {

         // Newline characters are not counted when computing distance
         // (but tab characters would still be considered).

         chunk++;
         continue;
      }
      x1 = chunk->fX;
      y1 = chunk->fY - ascent;
      x2 = chunk->fX + chunk->fDisplayWidth;
      y2 = chunk->fY + descent;

      if (x < x1) {
         xDiff = x1 - x;
      } else if (x >= x2) {
         xDiff = x - x2 + 1;
      } else {
         xDiff = 0;
      }

      if (y < y1) {
         yDiff = y1 - y;
      } else if (y >= y2) {
         yDiff = y - y2 + 1;
      } else {
         yDiff = 0;
      }
      if ((xDiff == 0) && (yDiff == 0)) {
         return 0;
      }
      dist = (int) TMath::Hypot((Double_t) xDiff, (Double_t) yDiff);
      if ((dist < minDist) || !minDist) {
         minDist = dist;
      }
      chunk++;
   }
   return minDist;
}

////////////////////////////////////////////////////////////////////////////////
/// Determines whether a text layout lies entirely inside, entirely outside,
/// or overlaps a given rectangle. Non-displaying space characters that occur
/// at the end of individual lines in the text layout are ignored for
/// intersection calculations.
///
/// The return value is -1 if the text layout is entirely outside of the
/// rectangle, 0 if it overlaps, and 1 if it is entirely inside of the
/// rectangle.
///
/// x, y -- Upper-left hand corner, in pixels, of rectangular area to compare
///         with text layout. Coordinates are with respect to the upper-left
///         hand corner of the text layout itself.
/// w, h -- The width and height of the above rectangular area, in pixels.

Int_t TGTextLayout::IntersectText(Int_t x, Int_t y, Int_t w, Int_t h) const
{
   Int_t result, i, x1, y1, x2, y2;
   LayoutChunk_t *chunk;
   Int_t left, top, right, bottom;

   // Scan the chunks one at a time, seeing whether each is entirely in,
   // entirely out, or overlapping the rectangle.  If an overlap is
   // detected, return immediately; otherwise wait until all chunks have
   // been processed and see if they were all inside or all outside.

   chunk = fChunks;

   left = x;
   top = y;
   right = x + w;
   bottom = y + h;

   result = 0;
   for (i = 0; i < fNumChunks; i++) {
      if (chunk->fStart[0] == '\n') {

         // Newline characters are not counted when computing area
         // intersection (but tab characters would still be considered).

         chunk++;
         continue;
      }
      x1 = chunk->fX;
      y1 = chunk->fY - fFont->fFM.fAscent;
      x2 = chunk->fX + chunk->fDisplayWidth;
      y2 = chunk->fY + fFont->fFM.fDescent;

      if ((right < x1) || (left >= x2) || (bottom < y1) || (top >= y2)) {
         if (result == 1) {
            return 0;
         }
         result = -1;
      } else if ((x1 < left) || (x2 >= right) || (y1 < top) || (y2 >= bottom)) {
         return 0;
      } else if (result == -1) {
         return 0;
      } else {
         result = 1;
      }
      chunk++;
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Outputs the contents of a text layout in Postscript format. The set of
/// lines in the text layout will be rendered by the user supplied Postscript
/// function. The function should be of the form:
///
///     justify x y string  function  --
///
/// Justify is -1, 0, or 1, depending on whether the following string should
/// be left, center, or right justified, x and y is the location for the
/// origin of the string, string is the sequence of characters to be printed,
/// and function is the name of the caller-provided function; the function
/// should leave nothing on the stack.
///
/// The meaning of the origin of the string (x and y) depends on the
/// justification. For left justification, x is where the left edge of the
/// string should appear. For center justification, x is where the center of
/// the string should appear. And for right justification, x is where the
/// right edge of the string should appear. This behavior is necessary
/// because, for example, right justified text on the screen is justified
/// with screen metrics. The same string needs to be justified with printer
/// metrics on the printer to appear in the correct place with respect to
/// other similarly justified strings. In all circumstances, y is the
/// location of the baseline for the string.
///
/// result is modified to hold the Postscript code that will render the text
/// layout.

void TGTextLayout::ToPostscript(TString *result) const
{
#define MAXUSE 128
   char buf[MAXUSE + 10];
   LayoutChunk_t *chunk;
   Int_t i, j, used, c, baseline;

   chunk = fChunks;
   baseline = chunk->fY;
   used = 0;
   buf[used++] = '(';

   for (i = 0; i < fNumChunks; i++) {
      if (baseline != chunk->fY) {
         buf[used++] = ')';
         buf[used++] = '\n';
         buf[used++] = '(';
         baseline = chunk->fY;
      }
      if (chunk->fNumDisplayChars <= 0) {
         if (chunk->fStart[0] == '\t') {
            buf[used++] = '\\';
            buf[used++] = 't';
         }
      } else {
         for (j = 0; j < chunk->fNumDisplayChars; j++) {
            c = UChar_t(chunk->fStart[j]);
            if ((c == '(') || (c == ')') || (c == '\\') || (c < 0x20) || (c >= UChar_t(0x7f))) {

            // Tricky point: the "03" is necessary in the sprintf
            // below, so that a full three digits of octal are
            // always generated. Without the "03", a number
            // following this sequence could be interpreted by
            // Postscript as part of this sequence.

               // coverity[secure_coding]
               sprintf(buf + used, "\\%03o", c);
               used += 4;
            } else {
               buf[used++] = c;
            }
            if (used >= MAXUSE) {
               buf[used] = '\0';
               result->Append(buf);
               used = 0;
            }
         }
      }
      if (used >= MAXUSE) {
         // If there are a whole bunch of returns or tabs in a row,
         // then buf[] could get filled up.

         buf[used] = '\0';
         result->Append(buf);
         used = 0;
      }
      chunk++;
   }
   buf[used++] = ')';
   buf[used++] = '\n';
   buf[used]   = '\0';

   result->Append(buf);
}

////////////////////////////////////////////////////////////////////////////////
/// Helper function for ComputeTextLayout(). Encapsulates a measured set of
/// characters in a chunk that can be quickly drawn.
///
/// Returns a pointer to the new chunk in the text layout. The text layout is
/// reallocated to hold more chunks as necessary.
///
/// Currently, ComputeTextLayout() stores contiguous ranges of "normal"
/// characters in a chunk, along with individual tab and newline chars in
/// their own chunks. All characters in the text layout are accounted for.

LayoutChunk_t *TGFont::NewChunk(TGTextLayout *layout, Int_t *maxPtr,
                                const char *start, Int_t numChars,
                                Int_t curX, Int_t newX, Int_t y) const
{
   LayoutChunk_t *chunk;
   Int_t i, maxChunks;

   maxChunks = *maxPtr;
   if (layout->fNumChunks == maxChunks) {
      if (maxChunks == 0) {
         maxChunks = 1;
      } else {
         maxChunks *= 2;
      }
      chunk = new LayoutChunk_t[maxChunks];

      if (layout->fNumChunks > 0) {
         for (i=0; i<layout->fNumChunks; ++i) chunk[i] = layout->fChunks[i];
         delete[] layout->fChunks;
      }
      layout->fChunks = chunk;
      *maxPtr = maxChunks;
   }

   chunk = &layout->fChunks[layout->fNumChunks];
   chunk->fStart = start;
   chunk->fNumChars = numChars;
   chunk->fNumDisplayChars = numChars;
   chunk->fX = curX;
   chunk->fY = y;
   chunk->fTotalWidth = newX - curX;
   chunk->fDisplayWidth = newX - curX;
   layout->fNumChunks++;

   return chunk;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw a string of characters on the screen. DrawCharsExp() expands
/// control characters that occur in the string to \X or \xXX sequences.
/// DrawChars() just draws the strings.
///
/// dst      -- Window or pixmap in which to draw.
/// gc       -- Graphics context for drawing characters.
/// source   -- Characters to be displayed. Need not be'\0' terminated.
///             For DrawChars(), all meta-characters (tabs, control
///             characters, and newlines) should be stripped out of the
///             string that is passed to this function. If they are not
///             stripped out, they will be displayed as regular printing
///             characters.
/// numChars -- Number of characters in string.
/// x, y     -- Coordinates at which to place origin of string when drawing.

void TGFont::DrawCharsExp(Drawable_t dst, GContext_t gc,
                          const char *source, Int_t numChars,
                          Int_t x, Int_t y) const
{
   const char *p;
   Int_t i, type;
   char buf[4];

   p = source;
   for (i = 0; i < numChars; i++) {
      type = fTypes[UChar_t(*p)];
      if (type != kCharNormal) {
         DrawChars(dst, gc, source, p - source, x, y);
         x += gVirtualX->TextWidth(fFontStruct, source, p - source);
         if (type == kCharReplace) {
            DrawChars(dst, gc, buf, GetControlCharSubst(UChar_t(*p), buf), x, y);
            x += fWidths[UChar_t(*p)];
         }
         source = p + 1;
      }
      p++;
   }

   DrawChars(dst, gc, source, p - source, x, y);
}

////////////////////////////////////////////////////////////////////////////////
/// Perform a quick sanity check to ensure we won't overflow the X
/// coordinate space.

void TGFont::DrawChars(Drawable_t dst, GContext_t gc,
                       const char *source, Int_t numChars,
                       Int_t x, Int_t y) const
{
   Int_t max_width =  gVirtualX->TextWidth(fFontStruct, "@", 1);

   if ((x + (max_width * numChars) > 0x7fff)) {
      int length;

      // The string we are being asked to draw is too big and would overflow
      // the X coordinate space. Unfortunatley X servers aren't too bright
      // and so they won't deal with this case cleanly. We need to truncate
      // the string before sending it to X.

      numChars = MeasureChars(source, numChars, 0x7fff - x, 0, &length);
   }

   gVirtualX->DrawString(dst, gc, x, y, source, numChars);

   if (fFA.fUnderline != 0) {
      gVirtualX->FillRectangle(dst, gc, x,  y + fUnderlinePos,
                          (UInt_t) gVirtualX->TextWidth(fFontStruct, source, numChars),
                          (UInt_t) fBarHeight);
   }
   if (fFA.fOverstrike != 0) {
      y -= fFM.fDescent + fFM.fAscent / 10;
      gVirtualX->FillRectangle(dst, gc, x, y,
                          (UInt_t) gVirtualX->TextWidth(fFontStruct, source, numChars),
                          (UInt_t) fBarHeight);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a font pool.

TGFontPool::TGFontPool(TGClient *client)
{
   fClient = client;
   fList   = new THashTable(50);
   fList->SetOwner();

   fNamedTable = new THashTable(50);
   fNamedTable->SetOwner();

   fUidTable = new THashTable(50);
   fUidTable->SetOwner();
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup font pool.

TGFontPool::~TGFontPool()
{
   delete fList;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the specified font.
/// The font can be one of the following forms:
///        XLFD (see X documentation)
///        "Family [size [style] [style ...]]"
/// Returns 0 if error or no font can be found.
/// If fixedDefault is false the "fixed" font will not be substituted
/// as fallback when the asked for font does not exist.

TGFont *TGFontPool::GetFont(const char *font, Bool_t fixedDefault)
{
   if (!font || !*font) {
      Error("GetFont", "argument may not be 0 or empty");
      return 0;
   }

   TGFont *f = (TGFont*)fList->FindObject(font);

   if (f) {
      f->AddReference();
      return f;
   }

   TNamedFont *nf = (TNamedFont*)fNamedTable->FindObject(font);

   if (nf) {
      // Construct a font based on a named font.
      nf->AddReference();
      f = GetFontFromAttributes(&nf->fFA, 0);

   } else {

      // Native font (aka string in XLFD format)?
      Int_t errsav = gErrorIgnoreLevel;
      gErrorIgnoreLevel = kFatal;

      f = GetNativeFont(font, fixedDefault);
      gErrorIgnoreLevel = errsav;

      if (!f) {
         FontAttributes_t fa;

         if (!ParseFontName(font, &fa)) {
            //fontCache.DeleteHashEntry(cacheHash);

            return 0;
         }

         // String contained the attributes inline.
         f = GetFontFromAttributes(&fa, 0);
      }
   }

   if (!f) return 0;

   fList->Add(f);

   f->SetRefCount(1);
   //f->cacheHash = cacheHash;
   f->fNamedHash = nf;

   f->MeasureChars("0", 1, 0, 0, &f->fTabWidth);

   if (!f->fTabWidth) {
      f->fTabWidth = f->fFM.fMaxWidth;
   }
   f->fTabWidth *= 8;

   // Make sure the tab width isn't zero (some fonts may not have enough
   // information to set a reasonable tab width).

   if (!f->fTabWidth) {
      f->fTabWidth = 1;
   }

   // Get information used for drawing underlines in generic code on a
   // non-underlined font.

   Int_t descent = f->fFM.fDescent;
   f->fUnderlinePos = descent/2;  // ==!== could be set by MakeFont()
   f->fUnderlineHeight = f->fFA.fPointsize/10;

   if (!f->fUnderlineHeight) {
      f->fUnderlineHeight = 1;
   }
   if (f->fUnderlinePos + f->fUnderlineHeight > descent) {

      // If this set of values would cause the bottom of the underline
      // bar to stick below the descent of the font, jack the underline
      // up a bit higher.

      f->fUnderlineHeight = descent - f->fUnderlinePos;

      if (!f->fUnderlineHeight) {
         f->fUnderlinePos--;
         f->fUnderlineHeight = 1;
      }
   }

   return f;
}

////////////////////////////////////////////////////////////////////////////////
/// Use font, i.e. increases ref count of specified font. Returns 0
/// if font is not found.

TGFont *TGFontPool::GetFont(const TGFont *font)
{
   TGFont *f = (TGFont*)fList->FindObject(font);

   if (f) {
      f->AddReference();
      return f;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Use font, i.e. increases ref count of specified font.

TGFont *TGFontPool::GetFont(FontStruct_t fs)
{
   TGFont *f = FindFont(fs);

   if (f) {
      f->AddReference();
      return f;
   }

   static int i = 0;
   f = MakeFont(0, fs, TString::Format("unknown-%d", i));
   fList->Add(f);
   i++;

   return f;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns font specified bay family, pixel/point size, weight and slant
///  negative value of ptsize means size in pixels
///  positive value of ptsize means size in points
///
///  For example:
///    TGFont *font = fpool->GetFont("helvetica", -9, kFontWeightNormal, kFontSlantRoman);
///    font->Print();

TGFont *TGFontPool::GetFont(const char *family, Int_t ptsize, Int_t weight, Int_t slant)
{
   const char *s;
   TString tmp;

   tmp = TString::Format("%s %d", family, ptsize);
   s = FindStateString(gWeightMap, weight);
   if (s) {
      tmp += " ";
      tmp + s;
   }
   s = FindStateString(gSlantMap, slant);
   if (s) {
      tmp += " ";
      tmp += s;
   }
   return GetFont(tmp.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Free font. If ref count is 0 delete font.

void TGFontPool::FreeFont(const TGFont *font)
{
   TGFont *f = (TGFont*) fList->FindObject(font);
   if (f) {
      if (f->RemoveReference() == 0) {
         if (font->fNamedHash) {

            // The font is being deleted. Determine if the associated named
            // font definition should and/or can be deleted too.

            TNamedFont *nf = (TNamedFont *) font->fNamedHash;

            if ((nf->RemoveReference() == 0) && (nf->fDeletePending != 0)) {
               fNamedTable->Remove(nf);
               delete nf;
            }
         }
         fList->Remove(f);
         delete font;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find font based on its font struct. Returns 0 if font is not found.

TGFont *TGFontPool::FindFont(FontStruct_t font) const
{
   TIter next(fList);
   TGFont *f = 0;

   while ((f = (TGFont*) next())) {
      if (f->fFontStruct == font) {
         return f;
      }
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Find font based on its font handle. Returns 0 if font is not found.

TGFont *TGFontPool::FindFontByHandle(FontH_t font) const
{
   TIter next(fList);
   TGFont *f = 0;

   while ((f = (TGFont*) next())) {
      if (f->fFontH == font) {
         return f;
      }
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Given a string, this procedure returns a unique identifier for the string.
///
/// This procedure returns a pointer to a new char string corresponding to
/// the "string" argument. The new string has a value identical to string
/// (strcmp will return 0), but it's guaranteed that any other calls to this
/// procedure with a string equal to "string" will return exactly the same
/// result (i.e. can compare pointer *values* directly, without having to
/// call strcmp on what they point to).

const char *TGFontPool::GetUid(const char *string)
{
   TObjString *obj = 0;
   obj = (TObjString*)fUidTable->FindObject(string);

   if (!obj) {
      obj = new TObjString(string);
      fUidTable->Add(obj);
   }

   return (const char *)obj->GetName();
}

////////////////////////////////////////////////////////////////////////////////
/// Return information about the font attributes as an array of strings.
///
/// An array of FONT_NUMFIELDS strings is returned holding the value of the
/// font attributes in the following order:
/// family size weight slant underline overstrike

char **TGFontPool::GetAttributeInfo(const FontAttributes_t *fa)
{
   Int_t i, num;
   const char *str = 0;

   char **result = new char*[FONT_NUMFIELDS];

   for (i = 0; i < FONT_NUMFIELDS; ++i) {
      str = 0;
      num = 0;

      switch (i) {
      case FONT_FAMILY:
         str = fa->fFamily;
         if (!str) str = "";
         break;

      case FONT_SIZE:
         num = fa->fPointsize;
         break;

      case FONT_WEIGHT:
         str = FindStateString(gWeightMap, fa->fWeight);
         break;

      case FONT_SLANT:
         str = FindStateString(gSlantMap, fa->fSlant);
         break;

      case FONT_UNDERLINE:
         num = fa->fUnderline;
         break;

      case FONT_OVERSTRIKE:
         num = fa->fOverstrike;
         break;
      }

      if (str) {
         int len = strlen(str)+1;
         result[i] = new char[len];
         strlcpy(result[i], str, len);
      } else {
         result[i] = new char[20];
         snprintf(result[i], 20, "%d", num);
      }
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Free attributes info.

void TGFontPool::FreeAttributeInfo(char **info)
{
   Int_t i;

   if (info) {
      for (i = 0; i < FONT_NUMFIELDS; ++i) {
         if (info[i]) {
            delete[] info[i];
         }
      }
      delete[] info;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// List all fonts in the pool.

void TGFontPool::Print(Option_t *opt) const
{
   fList->Print(opt);
}

////////////////////////////////////////////////////////////////////////////////
/// Save the used font as a C++ statement(s) on output stream out.

void TGFont::SavePrimitive(std::ostream &out, Option_t * /*= ""*/)
{
   char quote = '"';

   if (gROOT->ClassSaved(TGFont::Class())) {
      out << std::endl;
   } else {
      //  declare a font object to reflect required user changes
      out << std::endl;
      out << "   TGFont *ufont;         // will reflect user font changes" << std::endl;
   }
   out << "   ufont = gClient->GetFont(" << quote << GetName() << quote << ");" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

static char *GetToken(char *str)
{
   static char *p = 0;
   char *retp;

   if (str) p = str;

   if (!p) {
      return 0;
   }
   if (!*p) {
      return 0;
   }

   while (*p && ((*p == ' ') || (*p == '\t'))) {   // skip spaces
      ++p;
   }

   if (!*p) {
      return 0;
   }

   if (*p == '"') {  // quoted string
      retp = ++p;

      if (!*p) {
         return 0;
      }

      while (*p && (*p != '"')) {
         ++p;
      }

      if (*p == '"') {
         *p++ = '\0';
      }
   } else {
      retp = p;
      while (*p && (*p != ' ') && (*p != '\t')) {
         ++p;
      }
      if (*p) {
         *p++ = '\0';
      }
   }

   return retp;
}

////////////////////////////////////////////////////////////////////////////////
/// Converts a string into a set of font attributes that can be used to
/// construct a font.
///
/// The string can be one of the following forms:
///        XLFD (see X documentation)
///        "Family [size [style] [style ...]]"
///
/// The return value is kFALSE if the object was syntactically
/// invalid. Otherwise, fills the font attribute buffer with the values
/// parsed from the string and returns kTRUE. The structure must already be
/// properly initialized.

Bool_t TGFontPool::ParseFontName(const char *string, FontAttributes_t *fa)
{
   char *s;
   int n, result;

   XLFDAttributes_t xa;

   int len = strlen(string)+1;
   char *str = new char[len];
   strlcpy(str, string, len);

   if (*str == '-' || *str == '*') {

    // This appears to be an XLFD.

      xa.fFA = *fa;
      result = ParseXLFD(str, &xa);
      if (result) {
         *fa = xa.fFA;
         delete[] str;
         return kTRUE;  //OK
      }
   }

   // Wasn't an XLFD or "-option value" string. Try it as a
   // "font size style" list.

   s = GetToken(str);
   if (!s) {
      delete[] str;
      return kFALSE;
   }
   fa->fFamily = GetUid(s);

   s = GetToken(0);

   if (s) {
      char *end;

      fa->fPointsize = strtol(s, &end, 0);
      if ((errno == ERANGE) || (end == s)) {
         delete[] str;
         return kFALSE;
      }
   }

   while ((s = GetToken(0))) {
      n = FindStateNum(gWeightMap, s);
      if ((EFontWeight)n != kFontWeightUnknown) {
         fa->fWeight = n;
         continue;
      }
      n = FindStateNum(gSlantMap, s);
      // tell coverity that n is an integer value, and not an enum, even if
      // we compare it with an enum value (which is -1 in both case anyway)
      // coverity[mixed_enums]
      if ((EFontSlant)n != kFontSlantUnknown) {
         fa->fSlant = n;
         continue;
      }
      n = FindStateNum(gUnderlineMap, s);
      if (n) {
         fa->fUnderline = n;
         continue;
      }
      n = FindStateNum(gOverstrikeMap, s);
      if (n) {
         fa->fOverstrike = n;
         continue;
      }

      // Unknown style.

      delete[] str;
      return kFALSE;
   }

   delete[] str;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Break up a fully specified XLFD into a set of font attributes.
///
/// Return value is kFALSE if string was not a fully specified XLFD.
/// Otherwise, fills font attribute buffer with the values parsed from
/// the XLFD and returns kTRUE.
///
/// string -- Parseable font description string.
/// xa     -- XLFD attributes structure whose fields are to be modified.
///           Structure must already be properly initialized.

Bool_t TGFontPool::ParseXLFD(const char *string, XLFDAttributes_t *xa)
{
   char *src;
   const char *str;
   int i, j;
   char *field[XLFD_NUMFIELDS + 2];
   TString ds("");

   memset(field, '\0', sizeof (field));

   str = string;
   if (*str == '-') str++;

   ds.Append((char *) str);
   src = (char*)ds.Data();

   field[0] = src;
   for (i = 0; *src != '\0'; src++) {
      if (isupper(UChar_t(*src))) {
         *src = tolower(UChar_t(*src));
      }
      if (*src == '-') {
         i++;
         if (i > XLFD_NUMFIELDS) {
           break;
         }
         *src = '\0';
         field[i] = src + 1;
      }
   }

   // An XLFD of the form -adobe-times-medium-r-*-12-*-* is pretty common,
   // but it is (strictly) malformed, because the first * is eliding both
   // the Setwidth and the Addstyle fields. If the Addstyle field is a
   // number, then assume the above incorrect form was used and shift all
   // the rest of the fields up by one, so the number gets interpreted
   // as a pixelsize.

   if ((i > XLFD_ADD_STYLE) && (FieldSpecified(field[XLFD_ADD_STYLE]))) {
      if (atoi(field[XLFD_ADD_STYLE]) != 0) {
         for (j = XLFD_NUMFIELDS - 1; j >= XLFD_ADD_STYLE; j--) {
            field[j + 1] = field[j];
         }
         field[XLFD_ADD_STYLE] = 0;
         i++;
      }
   }

   // Bail if we don't have enough of the fields (up to pointsize).

   if (i < XLFD_FAMILY) {
      return kFALSE;
   }
   if (FieldSpecified(field[XLFD_FOUNDRY])) {
      xa->fFoundry = GetUid(field[XLFD_FOUNDRY]);
   }
   if (FieldSpecified(field[XLFD_FAMILY])) {
      xa->fFA.fFamily = GetUid(field[XLFD_FAMILY]);
   }
   if (FieldSpecified(field[XLFD_WEIGHT])) {
      xa->fFA.fWeight = FindStateNum(gXlfdgWeightMap, field[XLFD_WEIGHT]);
   }
   if (FieldSpecified(field[XLFD_SLANT])) {
      xa->fSlant = FindStateNum(gXlfdSlantMap, field[XLFD_SLANT]);
      if (xa->fSlant == kFontSlantRoman) {
         xa->fFA.fSlant = kFontSlantRoman;
      } else {
         xa->fFA.fSlant = kFontSlantItalic;
      }
   }
   if (FieldSpecified(field[XLFD_SETWIDTH])) {
      xa->fSetwidth = FindStateNum(gXlfdSetwidthMap, field[XLFD_SETWIDTH]);
   }
   // XLFD_ADD_STYLE ignored.

   // Pointsize in tenths of a point, but treat it as tenths of a pixel.

   if (FieldSpecified(field[XLFD_POINT_SIZE])) {
      if (field[XLFD_POINT_SIZE][0] == '[') {

         // Some X fonts have the point size specified as follows:
         //
         //      [ N1 N2 N3 N4 ]
         //
         // where N1 is the point size (in points, not decipoints!), and
         // N2, N3, and N4 are some additional numbers that I don't know
         // the purpose of, so I ignore them.

         xa->fFA.fPointsize = atoi(field[XLFD_POINT_SIZE] + 1);
      } else {
         char *end;

         xa->fFA.fPointsize = strtol(field[XLFD_POINT_SIZE], &end, 0);
         if (errno == ERANGE || end == field[XLFD_POINT_SIZE]) {
            return kFALSE;
         }
         xa->fFA.fPointsize /= 10;
      }
   }

   // Pixel height of font. If specified, overrides pointsize.

   if (FieldSpecified(field[XLFD_PIXEL_SIZE])) {
      if (field[XLFD_PIXEL_SIZE][0] == '[') {

         // Some X fonts have the pixel size specified as follows:
         //
         //      [ N1 N2 N3 N4 ]
         //
         // where N1 is the pixel size, and where N2, N3, and N4
         // are some additional numbers that I don't know
         // the purpose of, so I ignore them.

         xa->fFA.fPointsize = atoi(field[XLFD_PIXEL_SIZE] + 1);
      } else {
         char *end;

         xa->fFA.fPointsize = strtol(field[XLFD_PIXEL_SIZE], &end, 0);
         if (errno == ERANGE || end == field[XLFD_PIXEL_SIZE]) {
            return kFALSE;
         }
      }
   }
   xa->fFA.fPointsize = -xa->fFA.fPointsize;

   // XLFD_RESOLUTION_X ignored.

   // XLFD_RESOLUTION_Y ignored.

   // XLFD_SPACING ignored.

   // XLFD_AVERAGE_WIDTH ignored.

   if (FieldSpecified(field[XLFD_REGISTRY])) {
      xa->fCharset = FindStateNum(gXlfdCharsetMap, field[XLFD_REGISTRY]);
   }
   if (FieldSpecified(field[XLFD_ENCODING])) {
      xa->fEncoding = atoi(field[XLFD_ENCODING]);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Given a lookup table, map a string to a number in the table.
///
/// If strKey was equal to the string keys of one of the elements in the
/// table, returns the numeric key of that element. Returns the numKey
/// associated with the last element (the NULL string one) in the table
/// if strKey was not equal to any of the string keys in the table.

Int_t TGFontPool::FindStateNum(const FontStateMap_t *map, const char *strKey)
{
   const FontStateMap_t *m;

   if (!map->fStrKey) {
      return 0;
   }

   for (m = map; m->fStrKey != 0; m++) {
      if (strcasecmp(strKey, m->fStrKey) == 0) {
         return m->fNumKey;
      }
   }
   return m->fNumKey;
}

////////////////////////////////////////////////////////////////////////////////
/// Given a lookup table, map a number to a string in the table.
///
/// If numKey was equal to the numeric key of one of the elements in the
/// table, returns the string key of that element. Returns NULL if numKey
/// was not equal to any of the numeric keys in the table

const char *TGFontPool::FindStateString(const FontStateMap_t *map, Int_t numKey)
{
   for ( ; map->fStrKey != 0; map++) {
      if (numKey == map->fNumKey) return map->fStrKey;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper function for ParseXLFD(). Determines if a field in the XLFD was
/// set to a non-null, non-don't-care value.
///
/// The return value is kFALSE if the field in the XLFD was not set and
/// should be ignored, kTRUE otherwise.
///
/// field -- The field of the XLFD to check. Strictly speaking, only when
///          the string is "*" does it mean don't-care. However, an
///          unspecified or question mark is also interpreted as don't-care.

Bool_t TGFontPool::FieldSpecified(const char *field)
{
   char ch;

   if (!field) {
      return kFALSE;
   }
   ch = field[0];

   return (ch != '*' && ch != '?');
}

////////////////////////////////////////////////////////////////////////////////
/// Given a font, return a textual string identifying it.

const char *TGFontPool::NameOfFont(TGFont *font)
{
   return font->GetName();
}

////////////////////////////////////////////////////////////////////////////////
/// Return information about the font families that are available on the
/// current display.
///
/// An array of strings is returned holding a list of all the available font
/// families. The array is terminated with a NULL pointer.

char **TGFontPool::GetFontFamilies()
{
   Int_t i, numNames;
   char *family, *end, *p;

   THashTable familyTable(100);
   familyTable.SetOwner();

   char **nameList;
   char **dst;

   // coverity[returned_null]
   // coverity[dereference]
   nameList = gVirtualX->ListFonts("*", 10000, numNames);

   for (i = 0; i < numNames; i++) {
      if (nameList[i][0] != '-') {
         continue;
      }
      family = strchr(nameList[i] + 1, '-');
      if (!family) {
         continue;
      }
      family++;
      end = strchr(family, '-');
      if (!end) {
         continue;
      }
      *end = '\0';
      for (p = family; *p != '\0'; p++) {
         if (isupper(UChar_t(*p))) {
            *p = tolower(UChar_t(*p));
         }
      }
      if (!familyTable.FindObject(family)) {
         familyTable.Add(new TObjString(family));
      }
   }

   UInt_t entries = familyTable.GetEntries();
   dst = new char*[entries+1];

   TIter next(&familyTable);
   i = 0;
   TObject *obj;

   while ((obj = next())) {
      dst[i] = StrDup(obj->GetName());
      i++;
   }
   dst[i] = 0;

   gVirtualX->FreeFontNames(nameList);
   return dst;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete an array of families allocated GetFontFamilies() method

void TGFontPool::FreeFontFamilies(char **f)
{
   Int_t i;

   if (!f) return;

   for (i = 0; f[i] != 0; ++i) {
      delete[] f[i];
   }
   delete[] f;
}

////////////////////////////////////////////////////////////////////////////////
/// Given a desired set of attributes for a font, find a font with the
/// closest matching attributes and create a new TGFont object.
/// The return value is a pointer to a TGFont object that represents the
/// font with the desired attributes. If a font with the desired attributes
/// could not be constructed, some other font will be substituted
/// automatically.
///
/// Every call to this procedure returns a new TGFont object, even if the
/// specified attributes have already been seen before.

TGFont *TGFontPool::GetFontFromAttributes(FontAttributes_t *fa, TGFont *fontPtr)
{
   Int_t numNames, score, i, scaleable, pixelsize, xaPixelsize;
   Int_t bestIdx, bestScore, bestScaleableIdx, bestScaleableScore;
   XLFDAttributes_t xa;
   TString buf;
   char **nameList;
   TGFont *font;
   FontStruct_t fontStruct;
   const char *fmt, *family;

   family = fa->fFamily;
   if (!family) {
      family = "*";
   }
   pixelsize = -fa->fPointsize;

   if (pixelsize < 0) {
      double d;
      d = -pixelsize * 25.4/72;
      Int_t xx; Int_t yy; UInt_t ww; UInt_t hh;
      gVirtualX->GetWindowSize(gVirtualX->GetDefaultRootWindow(), xx, yy, ww, hh);
      d *= ww;

      d /= gVirtualX->ScreenWidthMM();
      d += 0.5;
      pixelsize = (int) d;
   }

   fontStruct = 0;

   // Couldn't find exact match. Now fall back to other available physical fonts.

   fmt = "-*-%.240s-*-*-*-*-*-*-*-*-*-*-*-*";
   buf = TString::Format(fmt, family);
   nameList = gVirtualX->ListFonts(buf.Data(), 32768, numNames);
   if (!numNames) {
      // Try getting some system font.

      buf = TString::Format(fmt, "fixed");
      // coverity[returned_null]
      // coverity[dereference]
      nameList = gVirtualX->ListFonts(buf.Data(), 32768, numNames);

      if (!numNames) {

getsystem:
         fontStruct = gVirtualX->LoadQueryFont("fixed");

         if (!fontStruct) {
            fontStruct = gVirtualX->LoadQueryFont("*");
            if (!fontStruct) {
               return 0;
            }
         }
         goto end;
      }
   }

   // Inspect each of the XLFDs and pick the one that most closely
   // matches the desired attributes.

   bestIdx = 0;
   bestScore = kMaxInt;
   bestScaleableIdx = 0;
   bestScaleableScore = kMaxInt;

   for (i = 0; i < numNames; i++) {
      score = 0;
      scaleable = 0;
      if (!ParseXLFD(nameList[i], &xa)) {
         continue;
      }
      xaPixelsize = -xa.fFA.fPointsize;

      // Since most people used to use -adobe-* in their XLFDs,
      // preserve the preference for "adobe" foundry. Otherwise
      // some applications looks may change slightly if another foundry
      // is chosen.

      if (xa.fFoundry && (strcasecmp(xa.fFoundry, "adobe") != 0)) {
         score += 3000;
      }
      if (!xa.fFA.fPointsize) {

         // A scaleable font is almost always acceptable, but the
         // corresponding bitmapped font would be better.

         score += 10;
         scaleable = 1;
      } else {

         // A font that is too small is better than one that is too big.

         if (xaPixelsize > pixelsize) {
            score += (xaPixelsize - pixelsize) * 120;
         } else {
            score += (pixelsize - xaPixelsize) * 100;
         }
      }

      score += TMath::Abs(xa.fFA.fWeight - fa->fWeight) * 30;
      score += TMath::Abs(xa.fFA.fSlant - fa->fSlant) * 25;

      if (xa.fSlant == kFontSlantOblique) {

         // Italic fonts are preferred over oblique.

         //score += 4;
      }
      if (xa.fSetwidth != kFontSWNormal) {

         // The normal setwidth is highly preferred.

         score += 2000;
      }
      if (xa.fCharset == kFontCSOther) {

         // The standard character set is highly preferred over
         // foreign languages charsets (because we don't support
         // other languages yet).

         score += 11000;
      }
      if ((xa.fCharset == kFontCSNormal) && (xa.fEncoding != 1)) {

         // The '1' encoding for the characters above 0x7f is highly
         // preferred over the other encodings.

         score += 8000;
      }
      if (scaleable) {
         if (score < bestScaleableScore) {
            bestScaleableIdx = i;
            bestScaleableScore = score;
         }
      } else {
         if (score < bestScore) {
            bestIdx = i;
            bestScore = score;
         }
      }
      if (!score) {
         break;
      }
   }

   // Now we know which is the closest matching scaleable font and the
   // closest matching bitmapped font. If the scaleable font was a
   // better match, try getting the scaleable font; however, if the
   // scalable font was not actually available in the desired pointsize,
   // fall back to the closest bitmapped font.

   fontStruct = 0;

   if (bestScaleableScore < bestScore) {
      char *str, *rest;

      // Fill in the desired pointsize info for this font.

tryscale:
      str = nameList[bestScaleableIdx];
      for (i = 0; i < XLFD_PIXEL_SIZE - 1; i++) {
         str = strchr(str + 1, '-');
      }
      rest = str;
      for (i = XLFD_PIXEL_SIZE - 1; i < XLFD_REGISTRY; i++) {
         rest = strchr(rest + 1, '-');
      }
      *str = '\0';
      buf = TString::Format("%.240s-*-%d-*-*-*-*-*%s", nameList[bestScaleableIdx], pixelsize, rest);
      *str = '-';
      fontStruct = gVirtualX->LoadQueryFont(buf.Data());
      bestScaleableScore = kMaxInt;
   }
   if (!fontStruct) {
      buf = nameList[bestIdx];
      fontStruct = gVirtualX->LoadQueryFont(buf.Data());

      if (!fontStruct) {

         // This shouldn't happen because the font name is one of the
         // names that X gave us to use, but it does anyhow.

         if (bestScaleableScore < kMaxInt) {
            goto tryscale;
         } else {
            gVirtualX->FreeFontNames(nameList);
            goto getsystem;
         }
      }
   }
   gVirtualX->FreeFontNames(nameList);

end:
   font = MakeFont(fontPtr, fontStruct, buf);
   font->fFA.fUnderline = fa->fUnderline;
   font->fFA.fOverstrike = fa->fOverstrike;

   return font;
}

////////////////////////////////////////////////////////////////////////////////
/// The return value is a pointer to an TGFont object that represents the
/// native font. If a native font by the given name could not be found,
/// the return value is NULL.
///
/// Every call to this procedure returns a new TGFont object, even if the
/// name has already been seen before. The caller should call FreeFont
/// when the font is no longer needed.

TGFont *TGFontPool::GetNativeFont(const char *name, Bool_t fixedDefault)
{
   FontStruct_t fontStruct;
   fixedDefault = fixedDefault && ((*name == '-') || (*name == '*'));
   fontStruct = fClient->GetFontByName(name, fixedDefault);

   if (!fontStruct) {
      return 0;
   }

   return MakeFont(0, fontStruct, name);
}

////////////////////////////////////////////////////////////////////////////////
/// Helper for GetNativeFont() and GetFontFromAttributes(). Creates and
/// intializes a new TGFont object.
///
/// font       -- If non-NULL, store the information in this existing TGFont
///               object, rather than creating a new one; the existing
///               contents of the font will be released. If NULL, a new
///               TGFont object is created.
/// fontStruct -- information about font.
/// fontName   -- The string passed to TVirtualX::LoadQueryFont() to construct the
///               fontStruct.

TGFont *TGFontPool::MakeFont(TGFont *font, FontStruct_t fontStruct,
                             const char *fontName)
{
   TGFont *newFont;

   Int_t i, width, firstChar, lastChar, n, replaceOK;
   char *p;
   char buf[4];
   XLFDAttributes_t xa;

   if (font) {
      gVirtualX->FreeFontStruct(font->fFontStruct);
      newFont = font;
   } else {
      newFont = new TGFont(fontName);
   }

   if (!ParseXLFD(fontName, &xa)) {
      newFont->fFA.fFamily = GetUid(fontName);
   } else {
      newFont->fFA = xa.fFA;
   }

   if (newFont->fFA.fPointsize < 0) {
      double d;
      Int_t xx; Int_t yy; UInt_t ww; UInt_t hh;
      gVirtualX->GetWindowSize(gVirtualX->GetDefaultRootWindow(), xx, yy, ww, hh);
      d = -newFont->fFA.fPointsize * 72/25.4;
      d *= gVirtualX->ScreenWidthMM();
      d /= ww;
      d += 0.5;
      newFont->fFA.fPointsize = (int) d;
   }

   Int_t ascent;
   Int_t descent;
   gVirtualX->GetFontProperties(fontStruct, ascent, descent);

   newFont->fFM.fAscent = ascent;
   newFont->fFM.fDescent = descent;
   newFont->fFM.fLinespace = ascent + descent;
   newFont->fFM.fMaxWidth = gVirtualX->TextWidth(fontStruct, "@", 1);
   newFont->fFM.fFixed = kTRUE;
   newFont->fFontStruct = fontStruct;
   newFont->fFontH      = gVirtualX->GetFontHandle(fontStruct);

   // Classify the characters.

   firstChar = 0x20; //fontStruct->min_char_or_byte2;
   lastChar = 0xff; //fontStruct->max_char_or_byte2;

   for (i = 0; i < 256; i++) {
      if ((i == 160) || (i == 173) || (i == 177) ||
          (i < firstChar) || (i > lastChar)) {
         newFont->fTypes[i] = kCharReplace;
      } else {
         newFont->fTypes[i] = kCharNormal;
      }
   }

   // Compute the widths for all the normal characters. Any other
   // characters are given an initial width of 0. Also, this determines
   // if this is a fixed or variable width font, by comparing the widths
   // of all the normal characters.

   char ch[2] = {0, 0};
   width = 0;
   for (i = 0; i < 256; i++) {
      if (newFont->fTypes[i] != kCharNormal) {
         n = 0;
      } else {
         ch[0] = i;
         n = gVirtualX->TextWidth(fontStruct, ch, 1);
      }
      newFont->fWidths[i] = n;
      if (n) {
         if (!width) {
            width = n;
         } else if (width != n) {
            newFont->fFM.fFixed = kFALSE;
         }
      }
   }

   // Compute the widths of the characters that should be replaced with
   // control character expansions. If the appropriate chars are not
   // available in this font, then control character expansions will not
   // be used; control chars will be invisible & zero-width.

   replaceOK = kTRUE;
   for (p = gHexChars; *p != '\0'; p++) {
      if ((UChar_t(*p) < firstChar) || (UChar_t(*p) > lastChar)) {
         replaceOK = kFALSE;
         break;
      }
   }
   for (i = 0; i < 256; i++) {
      if (newFont->fTypes[i] == kCharReplace) {
         if (replaceOK) {
            n = GetControlCharSubst(i, buf);
            for (; --n >= 0;) {
               newFont->fWidths[i] += newFont->fWidths[UChar_t(buf[n])];
            }
         } else {
            newFont->fTypes[i] = kCharSkip;
         }
      }
   }

   newFont->fUnderlinePos = descent >> 1;
   newFont->fBarHeight = newFont->fWidths[(int)'I']/3;

   if (newFont->fBarHeight == 0) {
      newFont->fBarHeight = 1;
   }

   if (newFont->fUnderlinePos + newFont->fBarHeight > descent) {

      // If this set of cobbled together values would cause the bottom of
      // the underline bar to stick below the descent of the font, jack
      // the underline up a bit higher.

      newFont->fBarHeight = descent - newFont->fUnderlinePos;

      if (!newFont->fBarHeight) {
         newFont->fUnderlinePos--;
         newFont->fBarHeight = 1;
      }
   }

   return newFont;
}

////////////////////////////////////////////////////////////////////////////////
/// When displaying text in a widget, a backslashed escape sequence is
/// substituted for control characters that occur in the text. Given a
/// control character, fill in a buffer with the replacement string that
/// should be displayed.
///
/// The return value is the length of the substitute string, buf is
/// filled with the substitute string; it is not '\0' terminated.
///
/// c   -- The control character to be replaced.
/// buf -- Buffer that gets replacement string. It only needs to be
///        4 characters long.

static Int_t GetControlCharSubst(Int_t c, char buf[4])
{
   buf[0] = '\\';

   if (((UInt_t)c < sizeof(gMapChars)) && (gMapChars[c] != 0)) {
      buf[1] = gMapChars[c];
      return 2;
   } else {
      buf[1] = 'x';
      buf[2] = gHexChars[(c >> 4) & 0xf];
      buf[3] = gHexChars[c & 0xf];
      return 4;
   }
}
