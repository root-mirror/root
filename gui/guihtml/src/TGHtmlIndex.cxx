// $Id: TGHtmlIndex.cxx,v 1.2 2007/05/20 23:22:13 rdm Exp $
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

// Routines that deal with indexes

#include <cctype>
#include <cstring>
#include <cstdlib>

#include "TGHtml.h"
#include "strlcpy.h"

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the Nth TGHtmlElement in the list. If there
/// is no Nth element, return 0 if flag==0 and return either the first
/// or last element (whichever is closest) if flag!=0

TGHtmlElement *TGHtml::TokenByIndex(int N, int /*flag*/)
{
   TGHtmlElement *p;
   int n;

   if (N == 0) return fPFirst;

   if (N > fNToken / 2) {
      // Start at the end and work back toward the beginning
      for (p = fPLast, n = fNToken; p; p = p->fPPrev) {
         if (p->fType != Html_Block) {
            if (p->fElId == N) break;
            --n;
         }
      }
   } else {
      // Start at the beginning and work forward
      for (p = fPFirst; p; p = p->fPNext) {
         if (p->fType != Html_Block) {
            --N;
            if (N == p->fElId) break;
         }
      }
   }

   return p;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the token number for the given TGHtmlElement

int TGHtml::TokenNumber(TGHtmlElement *p)
{
   //int n = 0;

   if (!p) return -1;
   return p->fElId;

///  while (p) {
///    if (p->fType != Html_Block) ++n;
///    p = p->fPPrev;
///  }
///
///  return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Find the maximum index for the given token

void TGHtml::MaxIndex(TGHtmlElement *p, int *pIndex, int isLast)
{
   if (p == 0) {
      *pIndex = 0;
   } else {
      switch (p->fType) {
         case Html_Text:
            *pIndex = p->fCount - isLast;
            break;

         case Html_Space:
            if (p->fStyle.fFlags & STY_Preformatted) {
               *pIndex = p->fCount - isLast;
            } else {
               *pIndex = 0;
            }
            break;

         default:
            *pIndex = 0;
            break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Given a Block and an x coordinate, find the Index of the character
/// that is closest to the given x coordinate.
///
/// The x-coordinate might specify a point to the left of the block,
/// in which case the procedure returns the first token and a character
/// index of 0.  Or the x-coordinate might specify a point to the right
/// of the block, in which case the last token is returned with an index
/// equal to its last character.

void TGHtml::FindIndexInBlock(TGHtmlBlock *pBlock, int x,
                              TGHtmlElement **ppToken, int *pIndex)
{
   TGHtmlElement *p;
   TGFont *font;
   int len;
   int n;

   p = pBlock->fPNext;
   font = GetFont(p->fStyle.fFont);
   if (x <= pBlock->fLeft) {
      *ppToken = p;
      *pIndex = 0;
      return;
   } else if (x >= pBlock->fRight) {
      *ppToken = p;
      *pIndex = 0;
      while (p && p->fType != Html_Block) {
         *ppToken = p;
         p = p->fPNext;
      }
      p = *ppToken;
      if (p && p->fType == Html_Text) {
         *pIndex = p->fCount - 1;
      }
      return;
   }
   if (pBlock->fN == 0) {
      *ppToken = p;
      *pIndex = 0;
   }
   n = font->MeasureChars(pBlock->fZ, pBlock->fN, x - pBlock->fLeft, 0, &len);
   *pIndex = 0;
   *ppToken = 0;
   while (p && n >= 0) {
      switch (p->fType) {
         case Html_Text:
            if (n < p->fCount) {
               *pIndex = n;
            } else {
               *pIndex = p->fCount - 1;
            }
            *ppToken = p;
            n -= p->fCount;
            break;

         case Html_Space:
            if (p->fStyle.fFlags & STY_Preformatted) {
               if (n < p->fCount) {
                  *pIndex = n;
               } else {
                  *pIndex = p->fCount - 1;
               }
               *ppToken = p;
               n -= p->fCount;
            } else {
               *pIndex = 0;
               *ppToken = p;
               --n;
            }
            break;

         default:
            break;
      }
      if (p) p = p->fPNext;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Convert an Element-based index into a Block-based index.
///
/// In other words, given a pointer to an element and an index
/// of a particular character within that element, compute a
/// pointer to the TGHtmlBlock used to display that character and
/// the index in the TGHtmlBlock of the character.

void TGHtml::IndexToBlockIndex(SHtmlIndex_t sIndex,
                               TGHtmlBlock **ppBlock, int *piIndex)
{
   int n = sIndex.fI;
   TGHtmlElement *p;

   if (sIndex.fP == 0) {
      *ppBlock = 0;
      *piIndex = 0;
      return;
   }
   p = sIndex.fP->fPPrev;
   while (p && p->fType != Html_Block) {
      switch (p->fType) {
         case Html_Text:
            n += p->fCount;
            break;
         case Html_Space:
            if (p->fStyle.fFlags & STY_Preformatted) {
               n += p->fCount;
            } else {
               n++;
            }
            break;
         default:
            break;
      }
      p = p->fPPrev;
   }
   if (p) {
      *ppBlock = (TGHtmlBlock *) p;
      *piIndex = n;
      return;
   }
   for (p = sIndex.fP; p && p->fType != Html_Block; p = p->fPNext) {}
   *ppBlock = (TGHtmlBlock *) p;
   *piIndex = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Modify an index for both pointer and char +/-/=N

int TGHtml::IndexMod(TGHtmlElement **pp, int *ip, char *cp)
{
   char nbuf[50];
   int i, x, cnt, ccnt[2], cflag[2];

   if (pp == 0 || !*pp) return -1;
   ccnt[0] = ccnt[1] = cflag[0] = cflag[1] = 0;
   x = 0;
   while (*cp && x < 2) {
      cnt = 0;
      i = 1;
      while (i < 45 && isdigit(cp[i])) {
         nbuf[i-1] = cp[i];
         i++;
      }
      if (i > 1) {
         nbuf[i-1] = 0;
         cnt = atoi(nbuf);
         if (cnt < 0) return -1;
      }
      switch (*cp) {
         case '+': if (i == 1) ccnt[x] = 1; else ccnt[x] = cnt; break;
         case '-': if (i == 1) ccnt[x] = -1; else ccnt[x] = -cnt; break;
         case '=': ccnt[x] = 0; cflag[x] = 1; break;
         default: return -1;
      }
      cp += i;
      ++x;
   }
   if (ccnt[0] > 0) {
      for (i = 0; i < ccnt[0] && (*pp)->fPNext; ++i) {
         *pp = (*pp)->fPNext;
         while ((*pp)->fType == Html_Block && (*pp)->fPNext) {
            *pp = (*pp)->fPNext;
         }
      }
   } else if (ccnt[0] < 0) {
      for (i = 0; ccnt[0] < i && (*pp)->fPPrev; --i) {
         //printf("i=%d, cnt=%d\n", i, ccnt[0]);
         *pp = (*pp)->fPPrev;
         while ((*pp)->fType == Html_Block && (*pp)->fPPrev) {
            *pp = (*pp)->fPPrev;
         }
      }
   }
   if (ccnt[1] > 0) {
      for (i = 0; i < ccnt[1]; ++i) (*ip)++;
   } else if (ccnt[1] < 0) {
      for (i = 0; i > ccnt[1]; --i) (*ip)--;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Given a base index name (without any modifiers) return a pointer
/// to the token described, and the character within that token.
///
/// Valid input forms include:
///
///       N.M          Token number N (with numbering starting at 1) and
///                    character number M (with numbering starting at 0).
///
///       M.X          Like above, but token is markup and X is an attribute.
///
///       begin        The start of all text
///
///       end          The end of all text
///
///       N.last       Last character of token number N.
///
///       N.end        One past last character of token number N.
///
///       sel.first    First character of the selection.
///
///       sel.last     Last character of the selection.
///
///       sel.end      On past last character of the selection.
///
///       insert       The character holding the insertion cursor.
///
///       @X,Y         The character a location X,Y of the clipping window.
///
///       &DOM         The DOM Address of a token.
///
/// Zero is returned if we are successful and non-zero if there is
/// any kind of error.
///
/// If the given token doesn't exist (for example if there are only 10
/// tokens and 11.5 is requested) then *ppToken is left pointing to NULL.
/// But the function still returns 0 for success.

int TGHtml::DecodeBaseIndex(const char *baseIx,
                            TGHtmlElement **ppToken, int *pIndex)
{
   int i, n, x, y;
   TGHtmlElement *p = 0;
   TGHtmlBlock *pBlock;
   TGHtmlBlock *pNearby;
   int dist = 1000000;
   int rc = 0;
   char buf[200], *base = buf, *suffix, *ep;

   strlcpy(buf, baseIx, sizeof(buf));

   while (isspace((unsigned char)*base)) base++;
   ep = base;
   while (*ep && !isspace((unsigned char)*ep)) ep++;
   *ep = 0;

   if ((suffix = strchr(base, ':'))) *suffix = 0;

   switch (*base) {
      case '1': case '2': case '3': case '4': case '5':
      case '6': case '7': case '8': case '9': case '0':
         // coverity[secure_coding]
         n = sscanf(base, "%d.%d", &x, &y);
         if (n > 0) {
            p = *ppToken = TokenByIndex(x, 0);
         }
         if (n == 2) {
            *pIndex = y;
         } else {
            for (i = 1; isdigit(base[i]); ++i) {}
            if (base[i] == 0) {
               *pIndex = 0;
            } else if (strcmp(&base[i], ".last") == 0) {
               MaxIndex(p, pIndex, 1);
            } else if (strcmp(&base[i], ".end") == 0) {
               MaxIndex(p, pIndex, 0);
               (*pIndex)++;
            } else {
               if (n == 1 && p && p->IsMarkup() && base[i] == '.' &&
                   p->MarkupArg(fZBase + i + 1, 0)) {
                  *pIndex = 0;
               } else {
                  rc = 1;
               }
            }
         }
         break;

      case 'b':
         if (strcmp(base, "begin") == 0) {
            p = *ppToken = fPFirst;
            *pIndex = 0;
         } else {
            rc = 1;
         }
         break;

      case 'e':
         if (strcmp(base, "end") == 0) {
            p = *ppToken = fPLast;
            MaxIndex(p, pIndex, 0);
         } else {
            rc = 1;
         }
         break;

      case 'l':
         if (strcmp(base, "last") == 0) {
            p = *ppToken = fPLast;
            MaxIndex(p, pIndex, 1);
         } else {
            rc = 1;
         }
         break;

      case 's':
         if (strcmp(base, "sel.first") == 0) {
            *ppToken = fSelBegin.fP;
            *pIndex = fSelBegin.fI;
         } else if (strcmp(base, "sel.last") == 0) {
            *ppToken = fSelEnd.fP;
            *pIndex = fSelEnd.fI;
         } else if (strcmp(base, "sel.end") == 0) {
            *ppToken = fSelEnd.fP;
            *pIndex = fSelEnd.fI + 1;
         } else {
            rc = 1;
         }
         break;

      case 'i':
         if (strcmp(baseIx, "insert") == 0) {
            *ppToken = fIns.fP;
            *pIndex = fIns.fI;
         } else {
            rc = 1;
         }
         break;

#if 0
      case '&':
         *pIndex = 0;
         if (DomIdLookup("id", base + 1, ppToken)) rc = 1;
         break;
#endif

      case '@':
         n = sscanf(base, "@%d,%d", &x, &y);
         if (n != 2) {
            rc = 1;
            break;
         }
         x += fVisible.fX;
         y += fVisible.fY;
         pNearby = 0;
         *ppToken = fPLast;
         *pIndex = 0;
         for (pBlock = fFirstBlock; pBlock; pBlock = pBlock->fBNext) {
            int dotest;
            if (pBlock->fN == 0) {
               switch (pBlock->fPNext->fType) {
                  case Html_LI:
                  case Html_IMG:
                  case Html_INPUT:
                  case Html_TEXTAREA:
                  case Html_SELECT:
                     dotest = 1;
                     break;
                  default:
                     dotest = 0;
                     break;
               }
            } else {
               dotest = 1;
            }
            if (dotest) {
               if (pBlock->fTop <= y && pBlock->fBottom >= y) {
                  if (pBlock->fLeft > x) {
                     if (pBlock->fLeft - x < dist) {
                        dist = pBlock->fLeft - x;
                        pNearby = pBlock;
                     }
                  } else if (pBlock->fRight < x) {
                     if (x - pBlock->fRight < dist) {
                        dist = x - pBlock->fRight;
                        pNearby = pBlock;
                     }
                  } else {
                     FindIndexInBlock(pBlock, x, ppToken, pIndex);
                     break;
                  }
               } else {
                  int distY;
                  int distX;

                  if (pBlock->fBottom < y) {
                     distY = y - pBlock->fBottom;
                  } else {
                     distY = pBlock->fTop - y;
                  }
                  if (pBlock->fLeft > x) {
                     distX = pBlock->fLeft - x;
                  } else if (pBlock->fRight < x) {
                     distX = x - pBlock->fRight;
                  } else {
                     distX = 0;
                  }
                  if (distX + 4*distY < dist) {
                     dist = distX + 4*distY;
                     pNearby = pBlock;
                  }
               }
            }
         }
         if (pBlock == 0) {
            if (pNearby) {
               FindIndexInBlock(pNearby, x, ppToken, pIndex);
            }
         }
         break;

      default:
         rc = 1;
         break;
   }
   if (suffix) IndexMod(ppToken, pIndex, suffix + 1);

   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// This routine decodes a complete index specification. A complete
/// index consists of the base specification followed by modifiers.

int TGHtml::GetIndex(const char *zIndex,
                     TGHtmlElement **ppToken, int *pIndex)
{
   return DecodeBaseIndex(zIndex, ppToken, pIndex);
}
