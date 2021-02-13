// $Id: TGHtmlParse.cxx,v 1.1 2007/05/04 17:07:01 brun Exp $
// Author:  Valeriy Onuchin   03/05/2007

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun, Fons Rademakers and Reiner Rohlfs *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**************************************************************************

    HTML widget for xclass. Based on tkhtml 1.28
    Copyright (C) 1997-2000 D. Richard Hipp <drh@acm.org>
    Copyright (C) 2002-2003 Hector Peraza.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library; if not, write to the Free
    Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

**************************************************************************/

// A tokenizer that converts raw HTML into a linked list of HTML elements.

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cctype>

#include "TGHtml.h"
#include "TGHtmlTokens.h"
#include "strlcpy.h"
#include "snprintf.h"

//----------------------------------------------------------------------

extern SHtmlTokenMap_t HtmlMarkupMap[];


/****************** Begin Escape Sequence Translator *************/

// The next section of code implements routines used to translate
// the '&' escape sequences of SGML to individual characters.
// Examples:
//
//         &amp;          &
//         &lt;           <
//         &gt;           >
//         &nbsp;         nonbreakable space
//

// Each escape sequence is recorded as an instance of the following
// structure

struct SgEsc_t {
   const char *fZName;      // The name of this escape sequence.  ex:  "amp"
   char  fValue[8];         // The value for this sequence.       ex:  "&"
   SgEsc_t *fPNext;         // Next sequence with the same hash on zName
};

// The following is a table of all escape sequences.  Add new sequences
// by adding entries to this table.

static struct SgEsc_t gEscSequences[] = {
   { "quot",      "\"",    0 },
   { "amp",       "&",     0 },
   { "lt",        "<",     0 },
   { "gt",        ">",     0 },
   { "nbsp",      " ",     0 },
   { "iexcl",     "\241",  0 },
   { "cent",      "\242",  0 },
   { "pound",     "\243",  0 },
   { "curren",    "\244",  0 },
   { "yen",       "\245",  0 },
   { "brvbar",    "\246",  0 },
   { "sect",      "\247",  0 },
   { "uml",       "\250",  0 },
   { "copy",      "\251",  0 },
   { "ordf",      "\252",  0 },
   { "laquo",     "\253",  0 },
   { "not",       "\254",  0 },
   { "shy",       "\255",  0 },
   { "reg",       "\256",  0 },
   { "macr",      "\257",  0 },
   { "deg",       "\260",  0 },
   { "plusmn",    "\261",  0 },
   { "sup2",      "\262",  0 },
   { "sup3",      "\263",  0 },
   { "acute",     "\264",  0 },
   { "micro",     "\265",  0 },
   { "para",      "\266",  0 },
   { "middot",    "\267",  0 },
   { "cedil",     "\270",  0 },
   { "sup1",      "\271",  0 },
   { "ordm",      "\272",  0 },
   { "raquo",     "\273",  0 },
   { "frac14",    "\274",  0 },
   { "frac12",    "\275",  0 },
   { "frac34",    "\276",  0 },
   { "iquest",    "\277",  0 },
   { "Agrave",    "\300",  0 },
   { "Aacute",    "\301",  0 },
   { "Acirc",     "\302",  0 },
   { "Atilde",    "\303",  0 },
   { "Auml",      "\304",  0 },
   { "Aring",     "\305",  0 },
   { "AElig",     "\306",  0 },
   { "Ccedil",    "\307",  0 },
   { "Egrave",    "\310",  0 },
   { "Eacute",    "\311",  0 },
   { "Ecirc",     "\312",  0 },
   { "Euml",      "\313",  0 },
   { "Igrave",    "\314",  0 },
   { "Iacute",    "\315",  0 },
   { "Icirc",     "\316",  0 },
   { "Iuml",      "\317",  0 },
   { "ETH",       "\320",  0 },
   { "Ntilde",    "\321",  0 },
   { "Ograve",    "\322",  0 },
   { "Oacute",    "\323",  0 },
   { "Ocirc",     "\324",  0 },
   { "Otilde",    "\325",  0 },
   { "Ouml",      "\326",  0 },
   { "times",     "\327",  0 },
   { "Oslash",    "\330",  0 },
   { "Ugrave",    "\331",  0 },
   { "Uacute",    "\332",  0 },
   { "Ucirc",     "\333",  0 },
   { "Uuml",      "\334",  0 },
   { "Yacute",    "\335",  0 },
   { "THORN",     "\336",  0 },
   { "szlig",     "\337",  0 },
   { "agrave",    "\340",  0 },
   { "aacute",    "\341",  0 },
   { "acirc",     "\342",  0 },
   { "atilde",    "\343",  0 },
   { "auml",      "\344",  0 },
   { "aring",     "\345",  0 },
   { "aelig",     "\346",  0 },
   { "ccedil",    "\347",  0 },
   { "egrave",    "\350",  0 },
   { "eacute",    "\351",  0 },
   { "ecirc",     "\352",  0 },
   { "euml",      "\353",  0 },
   { "igrave",    "\354",  0 },
   { "iacute",    "\355",  0 },
   { "icirc",     "\356",  0 },
   { "iuml",      "\357",  0 },
   { "eth",       "\360",  0 },
   { "ntilde",    "\361",  0 },
   { "ograve",    "\362",  0 },
   { "oacute",    "\363",  0 },
   { "ocirc",     "\364",  0 },
   { "otilde",    "\365",  0 },
   { "ouml",      "\366",  0 },
   { "divide",    "\367",  0 },
   { "oslash",    "\370",  0 },
   { "ugrave",    "\371",  0 },
   { "uacute",    "\372",  0 },
   { "ucirc",     "\373",  0 },
   { "uuml",      "\374",  0 },
   { "yacute",    "\375",  0 },
   { "thorn",     "\376",  0 },
   { "yuml",      "\377",  0 },
};


// The size of the handler hash table.  For best results this should
// be a prime number which is about the same size as the number of
// escape sequences known to the system.

#define ESC_HASH_SIZE (sizeof(gEscSequences)/sizeof(gEscSequences[0])+7)


// The hash table
//
// If the name of an escape sequence hashes to the value H, then
// gApEscHash[H] will point to a linked list of Esc structures, one of
// which will be the Esc structure for that escape sequence.

static struct SgEsc_t *gApEscHash[ESC_HASH_SIZE];


// Hash a escape sequence name. The value returned is an integer
// between 0 and ESC_HASH_SIZE-1, inclusive.

static int EscHash(const char *zName) {
   int h = 0;      // The hash value to be returned
   char c;         // The next character in the name being hashed

   while ((c = *zName) != 0) {
      h = h<<5 ^ h ^ c;
      zName++;
   }
   if (h < 0) h = -h;

   return h % ESC_HASH_SIZE;
}

#ifdef TEST
// Compute the longest and average collision chain length for the
// escape sequence hash table

static void EscHashStats()
{
   int i;
   int sum = 0;
   int max = 0;
   int cnt;
   int notempty = 0;
   struct SgEsc_t *p;

   for (i = 0; i < sizeof(gEscSequences) / sizeof(gEscSequences[0]); i++) {
      cnt = 0;
      p = gApEscHash[i];
      if (p) notempty++;
      while (p) {
         ++cnt;
         p = p->fPNext;
      }
      sum += cnt;
      if (cnt > max) max = cnt;
   }
   printf("Longest chain=%d  avg=%g  slots=%d  empty=%d (%g%%)\n",
          max, (double)sum/(double)notempty, i, i-notempty,
          100.0*(i-notempty)/(double)i);
}
#endif

// Initialize the escape sequence hash table

static void EscInit() {
   int i;  /* For looping thru the list of escape sequences */
   int h;  /* The hash on a sequence */

   for (i = 0; i < int(sizeof(gEscSequences) / sizeof(gEscSequences[i])); i++) {
/* #ifdef XCLASS_UTF_MAX */
#if 0
      int c = gEscSequences[i].value[0];
      xclass::UniCharToUtf(c, gEscSequences[i].value);
   }
#endif
      h = EscHash(gEscSequences[i].fZName);
      gEscSequences[i].fPNext = gApEscHash[h];
      gApEscHash[h] = &gEscSequences[i];
   }
#ifdef TEST
   EscHashStats();
#endif
}


// This table translates the non-standard microsoft characters between 0x80
// and 0x9f into plain ASCII so that the characters will be visible on Unix
// systems. Care is taken to translate the characters into values less than
// 0x80, to avoid UTF-8 problems.

static char gAcMsChar[] = {
   /* 0x80 */ 'C',
   /* 0x81 */ ' ',
   /* 0x82 */ ',',
   /* 0x83 */ 'f',
   /* 0x84 */ '"',
   /* 0x85 */ '.',
   /* 0x86 */ '*',
   /* 0x87 */ '*',
   /* 0x88 */ '^',
    /* 0x89 */ '%',
   /* 0x8a */ 'S',
   /* 0x8b */ '<',
   /* 0x8c */ 'O',
   /* 0x8d */ ' ',
   /* 0x8e */ 'Z',
   /* 0x8f */ ' ',
   /* 0x90 */ ' ',
   /* 0x91 */ '\'',
   /* 0x92 */ '\'',
   /* 0x93 */ '"',
   /* 0x94 */ '"',
   /* 0x95 */ '*',
   /* 0x96 */ '-',
   /* 0x97 */ '-',
   /* 0x98 */ '~',
   /* 0x99 */ '@',
   /* 0x9a */ 's',
   /* 0x9b */ '>',
   /* 0x9c */ 'o',
   /* 0x9d */ ' ',
   /* 0x9e */ 'z',
   /* 0x9f */ 'Y',
};


////////////////////////////////////////////////////////////////////////////////
/// Translate escape sequences in the string "z".  "z" is overwritten
/// with the translated sequence.
///
/// Unrecognized escape sequences are unaltered.
///
/// Example:
///
///      input  = "AT&amp;T &gt MCI"
///      output = "AT&T > MCI"

void HtmlTranslateEscapes(char *z)
{
   int from;   // Read characters from this position in z[]
   int to;     // Write characters into this position in z[]
   int h;      // A hash on the escape sequence
   struct SgEsc_t *p;  // For looping down the escape sequence collision chain
   static int isInit = 0;   // True after initialization

   from = to = 0;
   if (!isInit) {
      EscInit();
      isInit = 1;
   }
   while (z[from]) {
      if (z[from] == '&') {
         if (z[from+1] == '#') {
            int i = from + 2;
            int v = 0;
            while (isdigit(z[i])) {
               v = v*10 + z[i] - '0';
               i++;
            }
            if (z[i] == ';') { i++; }

            // Translate the non-standard microsoft characters in the range of
            // 0x80 to 0x9f into something we can see.

            if (v >= 0x80 && v < 0xa0) {
               v = gAcMsChar[v & 0x1f];
            }

            // Put the character in the output stream in place of the "&#000;".
            // How we do this depends on whether or not we are using UTF-8.

            z[to++] = v;
            from = i;
         } else {
            int i = from+1;
            int c;
            while (z[i] && isalnum(z[i])) ++i;
            c = z[i];
            z[i] = 0;
            h = EscHash(&z[from+1]);
            p = gApEscHash[h];
            while (p && strcmp(p->fZName, &z[from+1]) != 0) p = p->fPNext;
            z[i] = c;
            if (p) {
               int j;
               for (j = 0; p->fValue[j]; ++j) z[to++] = p->fValue[j];
               from = i;
               if (c == ';') from++;
            } else {
               z[to++] = z[from++];
            }
         }

         // Look for the non-standard microsoft characters between 0x80 and 0x9f
         // and translate them into printable ASCII codes. Separate algorithms
         // are required to do this for plain ascii and for utf-8.

      } else if (((unsigned char) z[from]) >= 0x80 &&
                 ((unsigned char) z[from]) < 0xa0) {
         z[to++] = gAcMsChar[z[from++] & 0x1f];
      } else {
         z[to++] = z[from++];
      }
   }
   z[to] = 0;
}

/******************* End Escape Sequence Translator ***************/

/******************* Begin HTML tokenizer code *******************/

// The following variable becomes TRUE when the markup hash table
// (stored in HtmlMarkupMap[]) is initialized.

static int gIsInit = 0;

// The hash table for HTML markup names.
//
// If an HTML markup name hashes to H, then gApMap[H] will point to
// a linked list of sgMap structure, one of which will describe the
// the particular markup (if it exists.)

static SHtmlTokenMap_t *gApMap[HTML_MARKUP_HASH_SIZE];

// Hash a markup name
//
// HTML markup is case insensitive, so this function will give the
// same hash regardless of the case of the markup name.
//
// The value returned is an integer between 0 and HTML_MARKUP_HASH_SIZE-1,
// inclusive.

static int HtmlHash(const char *zName) {
  int h = 0;
  char c;

  while ((c = *zName) != 0) {
    if (isupper(c)) {   // do we have to check for this??????
      c = tolower(c);
    }
    h = h<<5 ^ h ^ c;
    zName++;
  }
  if (h < 0) {
    h = -h;
  }

  return h % HTML_MARKUP_HASH_SIZE;
}


#ifdef TEST
// Compute the longest and average collision chain length for the
// markup hash table

static void HtmlHashStats() {
  int i;
  int sum = 0;
  int max = 0;
  int cnt;
  int notempty = 0;
  struct sgMap *p;

   for (i = 0; i < HTML_MARKUP_COUNT; i++) {
      cnt = 0;
      p = gApMap[i];
      if (p) notempty++;
      while (p) {
         cnt++;
         p = p->fPCollide;
      }
      sum += cnt;
      if (cnt > max) max = cnt;
  }

  printf("longest chain=%d  avg=%g  slots=%d  empty=%d (%g%%)\n",
         max, (double)sum/(double)notempty, i, i-notempty,
         100.0*(i-notempty)/(double)i);
}
#endif


// Initialize the escape sequence hash table

static void HtmlHashInit(void){
  int i;
  int h;  // The hash on a markup name

  for (i = 0; i < HTML_MARKUP_COUNT; i++) {
    h = HtmlHash(HtmlMarkupMap[i].fZName);
    HtmlMarkupMap[i].fPCollide = gApMap[h];
    gApMap[h] = &HtmlMarkupMap[i];
  }
#ifdef TEST
  HtmlHashStats();
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Append the given TGHtmlElement to the tokenizers list of elements

void TGHtml::AppendElement(TGHtmlElement *pElem)
{
   pElem->fPNext = 0;
   pElem->fPPrev = fPLast;
   if (fPFirst == 0) {
      fPFirst = pElem;
   } else {
      fPLast->fPNext = pElem;
   }
   fPLast = pElem;
   fNToken++;
}

////////////////////////////////////////////////////////////////////////////////
/// Insert token pNew before token p

void TGHtml::AppToken(TGHtmlElement *pNew, TGHtmlElement *p, int offs)
{
   if (offs < 0) {
      if (p) {
         offs = p->fOffs;
      } else {
         offs = fNText;
      }
   }

////if (p) { pNew->fStyle = p->fStyle; pNew->fFlags = p->fFlags; }

//  pNew->fCount = 0;
   pNew->fOffs = offs;
   pNew->fPNext = p;
   if (p) {
      pNew->fElId = p->fElId;
      p->fElId = ++fIdind;
      pNew->fPPrev = p->fPPrev;
      if (p->fPPrev) p->fPPrev->fPNext = pNew;
      if (fPFirst == p) fPFirst = pNew;
      p->fPPrev = pNew;
   } else {
      pNew->fElId = ++fIdind;
      AppendElement(pNew);
   }
   fNToken++;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the new column index following the given character.

static int NextColumn(int iCol, char c)
{
   switch (c) {
      case '\n': return 0;
      case '\t': return (iCol | 7) + 1;
      default:   return iCol+1;
   }
   /* NOT REACHED */
}

////////////////////////////////////////////////////////////////////////////////
/// Convert a string to all lower-case letters.

void ToLower(char *z)
{
   while (*z) {
      if (isupper(*z)) *z = tolower(*z);
      z++;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Process as much of the input HTML as possible. Construct new
/// TGHtmlElement objects and appended them to the list. Return
/// the number of characters actually processed.
///
/// This routine may invoke a callback procedure which could delete
/// the HTML widget.
///
/// This routine is not reentrant for the same HTML widget.  To
/// prevent reentrancy (during a callback), the p->fICol field is
/// set to a negative number. This is a flag to future invocations
/// not to reentry this routine. The p->fICol field is restored
/// before exiting, of course.

int TGHtml::Tokenize()
{
   char *z;             // The input HTML text
   int c;               // The next character of input
   int n;               // Number of characters processed so far
   int inpCol;          // Column of input
   int i, j;            // Loop counters
   int h;               // Result from HtmlHash()
   TGHtmlElement *pElem;// A new HTML element
   int selfClose;       // True for content free elements. Ex: <br/>
   int argc;            // The number of arguments on a markup
   SHtmlTokenMap_t *pMap; // For searching the markup name hash table
# define mxARG 200      // Maximum number of parameters in a single markup
   char *argv[mxARG];   // Pointers to each markup argument.
   int arglen[mxARG];   // Length of each markup argument
   //int rl, ol;
#ifdef pIsInMeachnism
   int pIsInScript = 0;
   int pIsInNoScript = 0;
   int pIsInNoFrames = 0;
#endif
   int sawdot = 0;
   int inLi = 0;

   static char null[1] = { "" };

   inpCol = fICol;
   n = fNComplete;
   z = fZText;
   if (inpCol < 0) return n;   // Prevents recursion
   fICol = -1;
   pElem = 0;

   while ((c = z[n]) != 0) {

      sawdot--;
      if (c == -64 && z[n+1] == -128) {
         n += 2;
         continue;
      }

      if (fPScript) {

         // We are in the middle of <SCRIPT>...</SCRIPT>.  Just look for
         // the </SCRIPT> markup.  (later:)  Treat <STYLE>...</STYLE> the
         // same way.

         TGHtmlScript *pScr = fPScript;
         const char *zEnd;
         int nEnd;
         //int curline, curch, curlast = n;
         int sqcnt;
         if (pScr->fType == Html_SCRIPT) {
            zEnd = "</script>";
            nEnd = 9;
         } else if (pScr->fType == Html_NOSCRIPT) {
            zEnd = "</noscript>";
            nEnd = 11;
         } else if (pScr->fType == Html_NOFRAMES) {
            zEnd = "</noframes>";
            nEnd = 11;
         } else {
            zEnd = "</style>";
            nEnd = 8;
         }
         if (pScr->fNStart < 0) {
            pScr->fNStart = n;
            pScr->fNScript = 0;
         }
         sqcnt = 0;
         for (i = n /*pScr->fNStart + pScr->fNScript*/; z[i]; i++) {
            if (z[i] == '\'' || z[i] == '"') {
               sqcnt++; // Skip if odd # quotes
            } else if (z[i] == '\n') {
               sqcnt = 0;
            }
            if (z[i] == '<' && z[i+1] == '/' &&
               strncasecmp(&z[i], zEnd, nEnd) == 0) {
               if (zEnd[3] == 'c' && ((sqcnt % 2) == 1)) continue;
               pScr->fNScript = i - n;
               fPScript = 0;
               n = i + nEnd;
               break;
            }
         }
         if (z[i] == 0) goto incomplete;
         if (fPScript) {
            pScr->fNScript = i - n;
            n = i;
         }
         else {
#ifdef pIsInMeachnism
            // If there is a script, execute it now and insert any output
            // to the html stream for parsing as html. (ie. client side scripting)

            if (pIsInScript && !pIsInNoScript && !pIsInNoFrames) {

               //for (curch = 0, curline = 1; curch <= curlast; curch++)
               //  if (z[curch] == '\n') curline++;

               // arglist in pElem and text pointers in pScr?
               // Inline scripts can contain unmatched brackets :-)
               //char varind[50];
               //sprintf(varind, "HtmlScrVar%d", p->varind++);
               //char savech = fZText[pScr->fNStart + pScr->fNScript];
               //fZText[pScr->fNStart + pScr->fNScript] = 0;
               //char *scriptBody = StrDup(fZText[pScr->fNStart]);
               //fZText[pScr->fNStart + pScr->fNScript] = savech;
               AdvanceLayout(p);
               inParse++;
               char *result = ProcessScript((TGHtmlScript *) pElem);  // pElem or pScr??
               inParse--;
               if (result) {
                  ol = fNAlloc;
                  rl = strlen(result);
                  fNAlloc += rl;
                  z = fZText = HtmlRealloc(z, ol+rl);
                  memmove(z + n + rl, z+n, ol - n);
                  memmove(z + n, result, rl);
               }
            }
            pIsInScript = 0;
            pIsInNoScript = 0;
            pIsInNoFrames = 0;
#endif
         }
         //continue;

      }
      else if (isspace((unsigned char)c)) {

         // White space
         for (i = 0;
             (c = z[n+i]) != 0 && isspace((unsigned char)c) && c != '\n' && c != '\r';
              i++) { }
         if (c == '\r' && z[n+i+1] == '\n') ++i;
#if 0  // this is certainly NOT OK, since it alters pre-formatted text
         if (sawdot == 1) {
            pElem = new TGHtmlTextElement(2);
            strcpy(((TGHtmlTextElement *)pElem)->fZText, " ");
            pElem->fElId = ++fIdind;
            pElem->fOffs = n;
            pElem->fCount = 1;
            AppendElement(pElem);
         }
#endif
         pElem = new TGHtmlSpaceElement;
         if (pElem == 0) goto incomplete;
         ((TGHtmlSpaceElement *)pElem)->fW = 0;
         pElem->fOffs = n;
         pElem->fElId = ++fIdind;
         if (c == '\n' || c == '\r') {
            pElem->fFlags = HTML_NewLine;
            pElem->fCount = 1;
            i++;
            inpCol = 0;
         } else {
            int iColStart = inpCol;
            pElem->fFlags = 0;
            for (j = 0; j < i; j++) {
               inpCol = NextColumn(inpCol, z[n+j]);
            }
            pElem->fCount = inpCol - iColStart;
         }
         AppendElement(pElem);
         n += i;

      }
      else if (c != '<' || fIPlaintext != 0 ||
              (!isalpha(z[n+1]) && z[n+1] != '/' && z[n+1] != '!' && z[n+1] != '?')) {

         // Ordinary text
         for (i = 1; (c = z[n+i]) != 0 && !isspace((unsigned char)c) && c != '<'; i++) {}
         if (z[n+i-1] == '.' || z[n+i-1] == '!' || z[n+i-1] == '?') sawdot = 2;
         if (c == 0) goto incomplete;
         if (fIPlaintext != 0 && z[n] == '<') {
            switch (fIPlaintext) {
               case Html_LISTING:
                  if (i >= 10 && strncasecmp(&z[n], "</listing>", 10) == 0) {
                     fIPlaintext = 0;
                     goto doMarkup;
                  }
                  break;

               case Html_XMP:
                  if (i >= 6 && strncasecmp(&z[n], "</xmp>", 6) == 0) {
                     fIPlaintext = 0;
                     goto doMarkup;
                  }
                  break;

               case Html_TEXTAREA:
                  if (i >= 11 && strncasecmp(&z[n], "</textarea>", 11) == 0) {
                     fIPlaintext = 0;
                     goto doMarkup;
                  }
                  break;

               default:
                  break;
            }
         }
         pElem = new TGHtmlTextElement(i);
         if (pElem == 0) goto incomplete;
         TGHtmlTextElement *tpElem = (TGHtmlTextElement *) pElem;
         tpElem->fElId = ++fIdind;
         tpElem->fOffs = n;
         strncpy(tpElem->fZText, &z[n], i);
         tpElem->fZText[i] = 0;
         AppendElement(pElem);
         if (fIPlaintext == 0 || fIPlaintext == Html_TEXTAREA) {
            HtmlTranslateEscapes(tpElem->fZText);
         }
         pElem->fCount = (Html_16_t) strlen(tpElem->fZText);
         n += i;
         inpCol += i;

      } else if (strncmp(&z[n], "<!--", 4) == 0) {

         // An HTML comment. Just skip it.
         for (i = 4; z[n+i]; i++) {
            if (z[n+i] == '-' && strncmp(&z[n+i], "-->", 3) == 0) break;
         }
         if (z[n+i] == 0) goto incomplete;

         pElem = new TGHtmlTextElement(i);
         if (pElem == 0) goto incomplete;
         TGHtmlTextElement *tpElem = (TGHtmlTextElement *) pElem;
         tpElem->fType = Html_COMMENT;
         tpElem->fElId = ++fIdind;
         tpElem->fOffs = n;
         strncpy(tpElem->fZText, &z[n+4], i-4);
         tpElem->fZText[i-4] = 0;
         tpElem->fCount = 0;
         AppendElement(pElem);

         pElem = new TGHtmlElement(Html_EndCOMMENT);
         AppToken(pElem, 0, n+4);

         for (j = 0; j < i+3; j++) {
           inpCol = NextColumn(inpCol, z[n+j]);
         }
         n += i + 3;

      }
      else {

         // Markup.
         //
         // First get the name of the markup
doMarkup:
         argc = 1;
         argv[0] = &z[n+1];
         for (i = 1;
             (c = z[n+i]) != 0 && !isspace((unsigned char)c) && c != '>' && (i < 2 || c != '/');
             i++) {}
         arglen[0] = i - 1;
         if (c == 0) goto incomplete;

         // Now parse up the arguments

         while (isspace((unsigned char)z[n+i])) ++i;
         while ((c = z[n+i]) != 0 && c != '>' && (c != '/' || z[n+i+1] != '>')) {
            if (argc > mxARG - 3) argc = mxARG - 3;
            argv[argc] = &z[n+i];
            j = 0;
            while ((c = z[n+i+j]) != 0 && !isspace((unsigned char)c) && c != '>' &&
                    c != '=' && (c != '/' || z[n+i+j+1] != '>')) ++j;
            arglen[argc] = j;
            if (c == 0) goto incomplete;
            i += j;
            while (isspace((unsigned char)c)) {
               i++;
               c = z[n+i];
            }
            if (c == 0) goto incomplete;
            argc++;
            if (c != '=') {
               argv[argc] = null;
               arglen[argc] = 0;
               argc++;
               continue;
            }
            i++;
            c = z[n+i];
            while (isspace((unsigned char)c)) {
               i++;
               c = z[n+i];
            }
            if (c == 0) goto incomplete;
            if (c == '\'' || c == '"') {
               int cQuote = c;
               i++;
               argv[argc] = &z[n+i];
               for (j = 0; (c = z[n+i+j]) != 0 && c != cQuote; j++) {}
               if (c == 0) goto incomplete;
               arglen[argc] = j;
               i += j+1;
            } else {
               argv[argc] = &z[n+i];
               for (j = 0; (c = z[n+i+j]) != 0 && !isspace((unsigned char)c) && c != '>'; j++) {}
               if (c == 0) goto incomplete;
               arglen[argc] = j;
               i += j;
            }
            argc++;
            while (isspace(z[n+i])) ++i;
         }
         if (c == '/') {
            i++;
            c = z[n+i];
            selfClose = 1;
         } else {
            selfClose = 0;
         }
         if (c == 0) goto incomplete;
         for (j = 0; j < i+1; j++) {
            inpCol = NextColumn(inpCol, z[n+j]);
         }
         n += i + 1;

         // Lookup the markup name in the hash table

         if (!gIsInit) {
            HtmlHashInit();
            gIsInit = 1;
         }
         c = argv[0][arglen[0]];
         argv[0][arglen[0]] = 0;
         h = HtmlHash(argv[0]);
         for (pMap = gApMap[h]; pMap; pMap = pMap->fPCollide) {
            if (strcasecmp(pMap->fZName, argv[0]) == 0) break;
         }
         argv[0][arglen[0]] = c;
         if (pMap == 0) continue;  // Ignore unknown markup

makeMarkupEntry:
         // Construct a TGHtmlMarkupElement object for this markup.

         pElem = MakeMarkupEntry(pMap->fObjType, pMap->fType, argc, arglen, argv);
         if (pElem == 0) goto incomplete;

         pElem->fElId = ++fIdind;
         pElem->fOffs = n;

         AddFormInfo(pElem);

         // The new markup has now been constructed in pElem. But before
         // appending it to the list, check to see if there is a special
         // handler for this markup type.

         if (ProcessToken(pElem, pMap->fZName, pMap->fType)) {
            // delete pElem;

            // Tricky, tricky. The user function might have caused the p->fZText
            // pointer to change, so renew our copy of that pointer.

            z = fZText;
            if (z == 0) {
               n = 0;
               inpCol = 0;
               goto incomplete;
            }
            continue;
         }

         // No special handler for this markup. Just append it to the
         // list of all tokens.

         AppendElement(pElem);
         switch (pMap->fType) {
            case Html_TABLE:
               break;

            case Html_PLAINTEXT:
            case Html_LISTING:
            case Html_XMP:
            case Html_TEXTAREA:
               fIPlaintext = pMap->fType;
               break;

            case Html_NOFRAMES:
               if (!fHasFrames) break;
#ifdef pIsInMeachnism
               pIsInNoFrames = 1;
#endif
            case Html_NOSCRIPT:
               break;
               // coverity[unreachable]
               if (!fHasScript) break;
#ifdef pIsInMeachnism
               pIsInNoScript = 1;
#endif
            case Html_SCRIPT:
#ifdef pIsInMeachnism
               pIsInScript = 1;
#endif
               // fallthrough
            case Html_STYLE:
               fPScript = (TGHtmlScript *) pElem;
               break;

            case Html_LI:
               if (!fAddEndTags) break;
               if (inLi) {
                  TGHtmlElement *e = new TGHtmlMarkupElement(Html_EndLI, 1, 0, 0);
                  AppToken(e, pElem, n);
               } else {
                  inLi = 1;
               }
               break;

            case Html_EndLI:
               inLi=0;
               break;

            case Html_EndOL:
            case Html_EndUL:
               if (!fAddEndTags) break;
               if (inLi) {
                  TGHtmlElement *e = new TGHtmlMarkupElement(Html_EndLI, 1, 0, 0);
                  AppToken(e, pElem, n);
               } else {
                  inLi = 0;
               }
               break;

            default:
               break;
         }

         // If this is self-closing markup (ex: <br/> or <img/>) then
         // synthesize a closing token.

         if (selfClose && argv[0][0] != '/' &&
             strcmp(&pMap[1].fZName[1], pMap->fZName) == 0) {
            selfClose = 0;
            pMap++;
            argc = 1;
            goto makeMarkupEntry;
         }
      }
   }

incomplete:
   fICol = inpCol;
   ////fPScript = 0;

   return n;
}

/************************** End HTML Tokenizer Code ***************************/

////////////////////////////////////////////////////////////////////////////////
/// Make one markup entry.

TGHtmlMarkupElement *TGHtml::MakeMarkupEntry(int objType, int type, int argc,
                                             int arglen[], char *argv[])
{
   TGHtmlMarkupElement *e;

   switch (objType) {
    case O_HtmlCell:
      e = new TGHtmlCell(type, argc, arglen, argv);
      break;

    case O_HtmlTable:
      e = new TGHtmlTable(type, argc, arglen, argv);
      break;

    case O_HtmlRef:
      e = new TGHtmlRef(type, argc, arglen, argv);
      break;

    case O_HtmlLi:
      e = new TGHtmlLi(type, argc, arglen, argv);
      break;

    case O_HtmlListStart:
      e = new TGHtmlListStart(type, argc, arglen, argv);
      break;

    case O_HtmlImageMarkup:
      e = new TGHtmlImageMarkup(type, argc, arglen, argv);
      break;

    case O_HtmlInput:
      e = new TGHtmlInput(type, argc, arglen, argv);
      break;

    case O_HtmlForm:
      e = new TGHtmlForm(type, argc, arglen, argv);
      break;

    case O_HtmlHr:
      e = new TGHtmlHr(type, argc, arglen, argv);
      break;

    case O_HtmlAnchor:
      e = new TGHtmlAnchor(type, argc, arglen, argv);
      break;

    case O_HtmlScript:
      e = new TGHtmlScript(type, argc, arglen, argv);
      break;

    case O_HtmlMapArea:
      e = new TGHtmlMapArea(type, argc, arglen, argv);
      break;

    default:
      e = new TGHtmlMarkupElement(type, argc, arglen, argv);
      break;
   }

   return e;
}

////////////////////////////////////////////////////////////////////////////////
/// Append text to the tokenizer engine.

void TGHtml::TokenizerAppend(const char *text)
{
   int len = strlen(text);

   if (fNText == 0) {
      fNAlloc = len + 100;
      fZText = new char [fNAlloc];
   } else if (fNText + len >= fNAlloc) {
      fNAlloc += len + 100;
      char *tmp = new char[fNAlloc];
      // coverity[secure_coding]
      strcpy(tmp, fZText);
      delete[] fZText;
      fZText = tmp;
   }

   if (fZText == 0) {
      fNText = 0;
      UNTESTED;
      return;
   }

   // coverity[secure_coding]
   strcpy(&fZText[fNText], text);
   fNText += len;
   fNComplete = Tokenize();
}

////////////////////////////////////////////////////////////////////////////////
/// This routine takes a text representation of a token, converts it into an
/// TGHtmlElement object and inserts it immediately prior to pToken. If pToken
/// is 0, then the newly created TGHtmlElement is appended.
///
/// This routine does nothing to resize, restyle, relayout or redisplay
/// the HTML. That is the calling routines responsibility.
///
/// Return the new TGHtmlElement object if successful. Return zero if
/// zType is not a known markup name.
///
///  pToken  - Insert before this. Append if pToken == 0
///  zType   - Type of markup. Ex: "/a" or "table"
///  zArgs   - List of arguments
///  offs    - Calculate offset, and insert changed text into fZText!

TGHtmlElement *TGHtml::InsertToken(TGHtmlElement *pToken,
                                  char *zType, char *zArgs, int offs)
{
   SHtmlTokenMap_t *pMap;     // For searching the markup name hash table
   int h;                   // The hash on zType
   TGHtmlElement *pElem;     // The new element
   //int nByte;               // How many bytes to allocate
   //int i;                   // Loop counter

   if (!gIsInit) {
      HtmlHashInit();
      gIsInit = 1;
   }

   if (strcmp(zType, "Text") == 0) {
      pElem = new TGHtmlTextElement(zArgs ? strlen(zArgs) : 0);
      if (pElem == 0) return 0;
      if (zArgs) {
         // coverity[secure_coding]
         strcpy (((TGHtmlTextElement *)pElem)->fZText, zArgs);
         pElem->fCount = (Html_16_t) strlen(zArgs);
      }
   } else if (!strcmp(zType, "Space")) {
      pElem = new TGHtmlSpaceElement();
      if (pElem == 0) return 0;
   } else {
      h = HtmlHash(zType);
      for (pMap = gApMap[h]; pMap; pMap = pMap->fPCollide) {
         if (strcasecmp(pMap->fZName, zType) == 0) break;
      }
      if (pMap == 0) return 0;
      if (zArgs == 0 || *zArgs == 0) {
         // Special case of no arguments. This is a lot easier...
         // well... now its the same thing!
         pElem = MakeMarkupEntry(pMap->fObjType, pMap->fType, 1, 0, 0);
         if (pElem == 0) return 0;
      } else {
         // The general case. There are arguments that need to be parsed
         // up. This is slower, but we gotta do it.
         //int  argc;
         //char **argv;
         //char *zBuf;

#if 0
      if (!SplitList(zArgs, &argc, &argv)) return 0;

      // shall we insert a dummy argv[0]?

      pElem = MakeMarkupEntry(pMap->fObjType, pMap->fType, argc/*+1??*/, 0, argv);
      if (pElem == 0) return 1;

      while (--argc >= 0) if (argv[argc]) delete[] argv[argc];
      delete[] argv;
#else
         return 0;
#endif
      }
   }

   pElem->fElId = ++fIdind;

   AppToken(pElem, pToken, offs);

   return pElem;
}

////////////////////////////////////////////////////////////////////////////////
/// Insert text into text token, or break token into two text tokens.
/// Also, handle backspace char by deleting text.
/// Should also handle newline char by splitting text.

int TGHtml::TextInsertCmd(int /*argc*/, char ** /*argv*/)
{
#if 0
  TGHtmlElement *p, *pElem;
  int i, l, n = 0;
  int idx = 0;
  int ptyp = Html_Unknown;
  int istxt = 0;
  char *cp = 0, c, *cp2;

  if (GetIndex(argv[3], &p, &i) != 0) {
    // sprintf(tmp, "malformed index: \"%s\"", argv[3]);
    return 0;
  }
  if (p) {
    ptyp = p->fType;
    if ((istxt = (ptyp == Html_Text))) {
      l = p->fCount;
      cp = ((TGHtmlTextElement *)p)->fZText;
    }
  }
  if (argv[2][0] == 'b') {  // Break text token into two.
    if (!istxt) return 1;
    if (i == 0 || i == l) return 1;
    pElem = InsertToken(p->fPNext, "Text", cp + i, -1);
    cp[i] = 0;
    p->fCount = i;
    return 1;
  }
  c = argv[4][0];
  if (!c) return 1;
  if (c == '\b') {
    if ((!istxt) || (!l) || (!i)) {
      if (!p) return 1;
      if (p->fType == Html_BR)
        RemoveElements(p, p);
      return 1;
    }
    if (p && l == 1) {
      RemoveElements(p, p);
      return 1;
    }
    if (i == l)
      cp[p->fCount] = 0;
    else
      memcpy(cp+i-1, cp+i, l-i+1);

    cp[--p->fCount] = 0;
    if (ins.i-- <= 0) ins.i = 0;
    ins.p = p;
    return 1;
  }
  if (c == '\n' || c == '\r') {
  }
  if (istxt) {
    char *cp;
    int t, j, alen = strlen(argv[4]);
    n = alen + l;

    TGHtmlTextElement *text = (TGHtmlTextElement *) p;

    if (text->fZText == (char*) ((&text->fZText)+1)) {
      cp = new char[n+1];
      strcpy(cp, text->fZText);
    } else {
      cp = new char[n+1];
      strcpy(cp, text->fZText);
    }
    cp2 = new char[alen+1];
    memcpy(cp2, argv[4], alen+1);
    HtmlTranslateEscapes(cp2);
    alen = strlen(cp2);
    memmove(cp+alen+i, cp+i, l-i+1);
    for (j = 0; j < alen; j++) cp[i+j] = cp2[j];
    delete[] cp2;
    delete[] text->fZText;
    text->fZText = cp;
    p->fCount = strlen(cp);
    ins.p = p;
    ins.i = i+alen;
  } else {
    p = InsertToken(p ? p->fPNext : 0, "Text", argv[4], -1);
    AddStyle(p);
    i = 0;
    ins.p = p;
    ins.i = 1;
  }
  if (p) {
    idx = p->base.id;
    AddStrOffset(p, argv[4], i);
  }
#endif
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns token map matching zType name.

SHtmlTokenMap_t *TGHtml::NameToPmap(char *zType)
{
   SHtmlTokenMap_t *pMap;     // For searching the markup name hash table
   int h;                   // The hash on zType

   if (!gIsInit) {
      HtmlHashInit();
      gIsInit = 1;
   }
   h = HtmlHash(zType);
   for (pMap = gApMap[h]; pMap; pMap = pMap->fPCollide) {
      if (strcasecmp(pMap->fZName, zType) == 0) break;
   }

   return pMap;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert a markup name into a type integer

int TGHtml::NameToType(char *zType)
{
   SHtmlTokenMap_t *pMap = NameToPmap(zType);
   return pMap ? pMap->fType : (int)Html_Unknown;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert a type into a symbolic name

const char *TGHtml::TypeToName(int type)
{
   if (type >= Html_A && type <= Html_EndXMP) {
      SHtmlTokenMap_t *pMap = gApMap[type - Html_A];
      return pMap->fZName;
   } else {
      return "???";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// For debugging purposes, print information about a token

char *TGHtml::DumpToken(TGHtmlElement *p)
{
//#ifdef DEBUG
   static char zBuf[200];
   int j;
   const char *zName;

   if (p == 0) {
      snprintf(zBuf, 200, "NULL");
      return zBuf;
   }
   switch (p->fType) {
      case Html_Text:
         snprintf(zBuf, 200, "text: \"%.*s\"", p->fCount, ((TGHtmlTextElement *)p)->fZText);
         break;

      case Html_Space:
         if (p->fFlags & HTML_NewLine) {
            snprintf(zBuf, 200, "space: \"\\n\"");
         } else {
            snprintf(zBuf, 200, "space: \" \"");
         }
         break;

      case Html_Block: {
         TGHtmlBlock *block = (TGHtmlBlock *) p;
         if (block->fN > 0) {
            int n = block->fN;
            if (n > 150) n = 150;
               snprintf(zBuf, 200, "<Block z=\"%.*s\">", n, block->fZ);
            } else {
               snprintf(zBuf, 200, "<Block>");
            }
            break;
      }

      default:
         if (p->fType >= HtmlMarkupMap[0].fType
             && p->fType <= HtmlMarkupMap[HTML_MARKUP_COUNT-1].fType) {
            zName = HtmlMarkupMap[p->fType - HtmlMarkupMap[0].fType].fZName;
         } else {
            zName = "Unknown";
         }
         snprintf(zBuf, 200, "markup (%d) <%s", p->fType, zName);
         for (j = 1 ; j < p->fCount; j += 2) {
            snprintf(&zBuf[strlen(zBuf)], 200-strlen(zBuf), " %s=\"%s\"",
                    ((TGHtmlMarkupElement *)p)->fArgv[j-1],
                    ((TGHtmlMarkupElement *)p)->fArgv[j]);
         }
         // coverity[secure_coding]
         strcat(zBuf, ">");
         break;
   }
   return zBuf;
//#else
//  return 0;
//#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Append all the arguments of the given markup to the given TGString.
///
/// Example:  If the markup is <IMG SRC=image.gif ALT="hello!">
/// then the following text is appended to the TGString:
///
///       "src image.gif alt hello!"
///
/// Notice how all attribute names are converted to lower case.
/// This conversion happens in the parser.

void TGHtml::AppendArglist(TGString *str, TGHtmlMarkupElement *pElem)
{
   int i;

   for (i = 0; i + 1 < pElem->fCount; i += 2) {
      str->Append(pElem->fArgv[i]);
      str->Append("=");
      str->Append(pElem->fArgv[i+1]);
      str->Append(" ");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns token name of html element p.

char *TGHtml::GetTokenName(TGHtmlElement *p)
{
   static char zBuf[200];
   //int j;
   const char *zName;

   zBuf[0] = 0;
   if (p == 0) {
      // coverity[secure_coding]: zBuf is large enough
      strcpy(zBuf, "NULL");
      return zBuf;
   }
   switch (p->fType) {
      case Html_Text:
      case Html_Space:
         break;

      case Html_Block:
         break;

      default:
         if (p->fType >= HtmlMarkupMap[0].fType &&
             p->fType <= HtmlMarkupMap[HTML_MARKUP_COUNT-1].fType) {
            zName = HtmlMarkupMap[p->fType - HtmlMarkupMap[0].fType].fZName;
         } else {
            zName = "Unknown";
         }
         strlcpy(zBuf, zName, sizeof(zBuf));
         break;
   }

   return zBuf;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns token map at location n.

SHtmlTokenMap_t* TGHtml::GetMarkupMap(int n)
{
   return HtmlMarkupMap+n;
}

////////////////////////////////////////////////////////////////////////////////
/// Return all tokens between the two elements as a string list.

TGString *TGHtml::ListTokens(TGHtmlElement *p, TGHtmlElement *pEnd)
{
   TGString *str;
   int i;
   const char *zName;
   char zLine[100];

   str = new TGString("");
   while (p && p != pEnd) {
      switch (p->fType) {
         case Html_Block:
            break;

         case Html_Text:
            str->Append("{ Text \"");
            str->Append(((TGHtmlTextElement *)p)->fZText);
            str->Append("\" } ");
            break;

         case Html_Space:
            snprintf(zLine, 100, "Space %d %d ",
                     p->fCount, (p->fFlags & HTML_NewLine) != 0);
            str->Append(zLine);
            break;

         case Html_Unknown:
            str->Append("Unknown ");
            break;

         default:
            str->Append("{ Markup ");
            if (p->fType >= HtmlMarkupMap[0].fType &&
                p->fType <= HtmlMarkupMap[HTML_MARKUP_COUNT-1].fType) {
               zName = HtmlMarkupMap[p->fType - HtmlMarkupMap[0].fType].fZName;
            } else {
               zName = "Unknown";
            }
            str->Append(zName);
            str->Append(" ");
            for (i = 0; i < p->fCount; ++i) {
               str->Append(((TGHtmlMarkupElement *)p)->fArgv[i]);
               str->Append(" ");
            }
            str->Append("} ");
            break;
      }
      p = p->fPNext;
   }

   return str;
}

////////////////////////////////////////////////////////////////////////////////
/// Print a list of tokens

void TGHtml::PrintList(TGHtmlElement *first, TGHtmlElement *last)
{
   TGHtmlElement *p;

   for (p = first; p != last; p = p->fPNext) {
      if (p->fType == Html_Block) {
         TGHtmlBlock *block = (TGHtmlBlock *) p;
         const char *z = block->fZ;
         int n = block->fN;
         if (n == 0 || z == 0) {
            n = 1;
            z = "";
         }
         printf("Block flags=%02x cnt=%d x=%d..%d y=%d..%d z=\"%.*s\"\n",
                p->fFlags, p->fCount, block->fLeft, block->fRight,
                block->fTop, block->fBottom, n, z);
      } else {
         printf("Token font=%2d color=%2d align=%d flags=0x%04x name=%s\n",
                p->fStyle.fFont, p->fStyle.fColor,
                p->fStyle.fAlign, p->fStyle.fFlags, DumpToken(p));
      }
   }
}
