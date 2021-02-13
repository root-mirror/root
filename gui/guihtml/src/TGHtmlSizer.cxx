// $Id: TGHtmlSizer.cxx,v 1.2 2007/05/04 20:33:16 rdm Exp $
// Author:  Valeriy Onuchin   03/05/2007

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

// Routines used to compute the style and size of individual elements.

#include <cstring>
#include <cstdlib>
#include <cctype>

#include "TGHtml.h"
#include "TImage.h"
#include "TVirtualX.h"
#include "snprintf.h"

////////////////////////////////////////////////////////////////////////////////
/// Get the current rendering style. In other words, get the style
/// that is currently on the top of the style stack.

SHtmlStyle_t TGHtml::GetCurrentStyle()
{
   SHtmlStyle_t style;

   if (fStyleStack) {
      style = fStyleStack->fStyle;
   } else {
      style.fFont = NormalFont(2);
      style.fColor = COLOR_Normal;
      style.fBgcolor = COLOR_Background;
      style.fSubscript = 0;
      style.fAlign = ALIGN_Left;
      style.fFlags = 0;
      style.fExpbg = 0;
   }

   return style;
}

////////////////////////////////////////////////////////////////////////////////
/// Push a new rendering style onto the stack.
///
///  tag   - Tag for this style. Normally the end-tag such as \</h3\> or \</em\>.
///  style - The style to push

void TGHtml::PushStyleStack(int tag, SHtmlStyle_t style)
{
   SHtmlStyleStack_t *p;

   p = new SHtmlStyleStack_t;
   p->fPNext = fStyleStack;
   p->fType = tag;
   p->fStyle = style;
   fStyleStack = p;
}

////////////////////////////////////////////////////////////////////////////////
/// Pop a rendering style off of the stack.
///
/// The top-most style on the stack should have a tag equal to "tag".
/// If not, then we have an HTML coding error.  Perhaps something
/// like this:  "Some text <em>Enphasized</i> more text".  It is an
/// interesting problem to figure out how to respond sanely to this
/// kind of error.  Our solution is to keep popping the stack until
/// we find the correct tag, or until the stack is empty.

SHtmlStyle_t TGHtml::PopStyleStack(int tag)
{
   int i, type;
   SHtmlStyleStack_t *p;
   static Html_u8_t priority[Html_TypeCount+1];

   if (priority[Html_TABLE] == 0) {
      for (i = 0; i <= Html_TypeCount; i++) priority[i] = 1;
      priority[Html_TD] = 2;
      priority[Html_EndTD] = 2;
      priority[Html_TH] = 2;
      priority[Html_EndTH] = 2;
      priority[Html_TR] = 3;
      priority[Html_EndTR] = 3;
      priority[Html_TABLE] = 4;
      priority[Html_EndTABLE] = 4;
   }
   if (tag <= 0 || tag > Html_TypeCount) {
      CANT_HAPPEN;
      return GetCurrentStyle();
   }
   while ((p = fStyleStack) != 0) {
      type = p->fType;
      if (type <= 0 || type > Html_TypeCount) {
         CANT_HAPPEN;
         return GetCurrentStyle();
      }
      if (type != tag && priority[type] > priority[tag]) {
         return GetCurrentStyle();
      }
      fStyleStack = p->fPNext;
      delete p;
      if (type == tag) break;
   }

   return GetCurrentStyle();
}

////////////////////////////////////////////////////////////////////////////////
/// Change the font size on the given style by the delta-amount given

static void ScaleFont(SHtmlStyle_t *pStyle, int delta)
{
   int size = FontSize(pStyle->fFont) + delta;

   if (size < 0) {
      delta -= size;
   } else if (size > 6) {
      delta -= size-6;
   }

   pStyle->fFont += delta;
}

////////////////////////////////////////////////////////////////////////////////
/// Add the STY_Invisible style to every token between p_first and p_last.

void TGHtml::MakeInvisible(TGHtmlElement *p_first, TGHtmlElement *p_last)
{
   if (p_first == 0) return;
   p_first = p_first->fPNext;
   while (p_first && p_first != p_last) {
      p_first->fStyle.fFlags |= STY_Invisible;
      p_first = p_first->fPNext;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// For the markup \<a href=XXX\>, find out if the URL has been visited
/// before or not.  Return COLOR_Visited or COLOR_Unvisited, as
/// appropriate.

int TGHtml::GetLinkColor(const char *zURL)
{
   return IsVisited(zURL) ? COLOR_Visited : COLOR_Unvisited;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns coordinates of string str.

static int *GetCoords(const char *str, int *nptr)
{
   const char *cp = str;
   char *ncp;
   int  *cr, i, n = 0, sz = 4;

   cr = new int[sz];
   while (cp) {
      while (*cp && (!isdigit(*cp))) cp++;
      if ((!*cp) || (!isdigit(*cp))) break;
      cr[n] = (int) strtol(cp, &ncp, 10);
      if (cp == ncp) break;
      cp = ncp;
      n++;
      if (n >= sz) {
         int *tmp = new int[sz+4];
         for (i = 0; i < sz; ++i) tmp[i] = cr[i];
         delete[] cr;
         cr = tmp;
         sz += 4;
      }
   }
   *nptr = n;

   return cr;
}

////////////////////////////////////////////////////////////////////////////////
/// This routine adds information to the input texts that doesn't change
/// when the display is resized or when new fonts are selected, etc.
/// Mostly this means adding style attributes.  But other constant
/// information (such as numbering on <li> and images used for <IMG>)
/// is also obtained.  The key is that this routine is only called
/// once, where the sizer and layout routines can be called many times.
///
/// This routine is called whenever the list of elements grows.  The
/// style stack is stored as part of the HTML widget so that we can
/// always continue where we left off the last time.
///
/// In addition to adding style, this routine will invoke methods
/// needed to acquire information about a markup. The IsVisitied()
/// method is called for each <a> and the GetImage() is called
/// for each <IMG> or for each <LI> that has a SRC= field.
///
/// When a markup is inserted or deleted from the token list, the
/// style routine must be completely rerun from the beginning.  So
/// what we said above, that this routine is only run once, is not
/// strictly true.

void TGHtml::AddStyle(TGHtmlElement *p)
{
   SHtmlStyle_t style;       // Current style
   int size;                 // A new font size
   int i;                    // Loop counter
   int paraAlign;            // Current paragraph alignment
   int rowAlign;             // Current table row alignment
   SHtmlStyle_t nextStyle;   // Style for next token if useNextStyle==1
   int useNextStyle = 0;     // True if nextStyle is valid
   const char *z;            // A tag parameter's value

   // The size of header fonts relative to the current font size
   static int header_sizes[] = { +2, +1, 1, 1, -1, -1 };

   // Don't allow recursion
   if (fFlags & STYLER_RUNNING) return;
   fFlags |= STYLER_RUNNING;

   // Load the style state out of the TGHtml object and into local
   // variables. This is purely a matter of convenience...

   style = GetCurrentStyle();
   nextStyle = style;   //ia: nextStyle was not initialized
   paraAlign = fParaAlignment;
   rowAlign = fRowAlignment;

   // Loop over tokens
   while (fPFirst && p) {
      switch (p->fType) {
         case Html_A:
            if (fAnchorStart) {
               style = PopStyleStack(Html_EndA);
               fAnchorStart = 0;
               fAnchorFlags = 0;
            }
            z = p->MarkupArg("href", 0);
            if (z) {
               style.fColor = GetLinkColor(z);
               if (fUnderlineLinks) style.fFlags |= STY_Underline;
               fAnchorFlags |= STY_Anchor;
               PushStyleStack(Html_EndA, style);
               fAnchorStart = (TGHtmlAnchor *) p;
            }
            break;

         case Html_EndA:
            if (fAnchorStart) {
               ((TGHtmlRef *)p)->fPOther = fAnchorStart;
               style = PopStyleStack(Html_EndA);
               fAnchorStart = 0;
               fAnchorFlags = 0;
            }
            break;

         case Html_MAP:
         break;

         case Html_EndMAP:
         break;

         case Html_AREA: {
            TGHtmlMapArea *area = (TGHtmlMapArea *) p;
            z = p->MarkupArg("shape", 0);
            area->fMType = HTML_MAP_RECT;
            if (z) {
               if (strcasecmp(z, "circle") == 0) {
                  area->fMType = HTML_MAP_CIRCLE;
               } else if (strcasecmp(z,"poly") == 0) {
                  area->fMType = HTML_MAP_POLY;
               }
            }
            z = p->MarkupArg("coords", 0);
            if (z) {
               area->fCoords = GetCoords(z, &area->fNum);
            }
            break;
         }

         case Html_ADDRESS:
         case Html_EndADDRESS:
         case Html_BLOCKQUOTE:
         case Html_EndBLOCKQUOTE:
            paraAlign = ALIGN_None;
            break;

         case Html_APPLET:
            if (0 /* has ProcessApplet() */) {
               nextStyle = style;
               nextStyle.fFlags |= STY_Invisible;
               PushStyleStack(Html_EndAPPLET, nextStyle);
               useNextStyle = 1;
            } else {
               PushStyleStack(Html_EndAPPLET, style);
            }
            break;

         case Html_B:
            style.fFont = BoldFont(style.fFont);
            PushStyleStack(Html_EndB, style);
            break;

         case Html_BODY:
            z = p->MarkupArg("text", 0);
            if (z) {
               //FreeColor(fApColor[COLOR_Normal]);
               fApColor[COLOR_Normal] = AllocColor(z);
            }
            z = p->MarkupArg("bgcolor", 0);
            if (z) {
               //FreeColor(fApColor[COLOR_Background]);
               fApColor[COLOR_Background] = AllocColor(z);
               SetBackgroundColor(fApColor[COLOR_Background]->fPixel);
               SetBackgroundPixmap(0);
            }
            z = p->MarkupArg("link", 0);
            if (z) {
               //FreeColor(fApColor[COLOR_Unvisited]);
               fApColor[COLOR_Unvisited] = AllocColor(z);
            }
            z = p->MarkupArg("vlink", 0);
            if (z) {
               //FreeColor(fApColor[COLOR_Visited]);
               fApColor[COLOR_Visited] = AllocColor(z);
            }
            z = p->MarkupArg("alink", 0);
            if (z) {
            }
            z = p->MarkupArg("background", 0);
            if (z) {
               z = ResolveUri(z);
               if (z) {
                  TImage *img = LoadImage(z, 0, 0);
                  if (img) {
#if 0
                     SetupBackgroundPic(img->GetPicture());
#else
                     GCValues_t gcv;
                     // unsigned int mask;

                     // mask = kGCTile | kGCFillStyle | kGCGraphicsExposures;
                     gcv.fTile = img->GetPixmap();
                     gcv.fFillStyle = kFillTiled;
                     gcv.fGraphicsExposures = kTRUE;
                     fCanvas->SetBackgroundPixmap(img->GetPixmap());

                     gVirtualX->ChangeGC(fWhiteGC.GetGC(), &gcv);

                     //NeedRedraw(TGRectangle(fVisible, fCanvas->GetSize()));
#endif
                     fBgImage = img;//delete img;
                  }
                  delete [] z;
               }
            }
            break;

         case Html_EndBODY:
            break;

         case Html_EndAPPLET:
         case Html_EndB:
         case Html_EndBIG:
         case Html_EndCENTER:
         case Html_EndCITE:
         case Html_EndCODE:
         case Html_EndCOMMENT:
         case Html_EndDFN:
         case Html_EndEM:
         case Html_EndFONT:
         case Html_EndI:
         case Html_EndKBD:
         case Html_EndMARQUEE:
         case Html_EndNOBR:
         case Html_EndNOFRAMES:
         case Html_EndNOSCRIPT:
         case Html_EndNOEMBED:
         case Html_EndS:
         case Html_EndSAMP:
         case Html_EndSMALL:
         case Html_EndSTRIKE:
         case Html_EndSTRONG:
         case Html_EndSUB:
         case Html_EndSUP:
         case Html_EndTITLE:
         case Html_EndTT:
         case Html_EndU:
         case Html_EndVAR:
            style = PopStyleStack(p->fType);
            break;

         case Html_BASE:
            z = p->MarkupArg("href", 0);
            if (z) {
               char *z1 = ResolveUri(z);
               if (z1 != 0) {
                  if (fZBaseHref) delete[] fZBaseHref;
                  fZBaseHref = z1;
               }
            }
            break;

         case Html_EndDIV:
            paraAlign = ALIGN_None;
            style = PopStyleStack(p->fType);
            break;

         case Html_EndBASEFONT:
            style = PopStyleStack(Html_EndBASEFONT);
            style.fFont = FontFamily(style.fFont) + 2;
            break;

         case Html_BIG:
            ScaleFont(&style, 1);
            PushStyleStack(Html_EndBIG, style);
            break;

         case Html_CAPTION:
            paraAlign = p->GetAlignment(paraAlign);
            break;

         case Html_EndCAPTION:
            paraAlign = ALIGN_None;
            break;

         case Html_CENTER:
            paraAlign = ALIGN_None;
            style.fAlign = ALIGN_Center;
            PushStyleStack(Html_EndCENTER, style);
            break;

         case Html_CITE:
            style.fFont = ItalicFont(style.fFont);
            PushStyleStack(Html_EndCITE, style);
            break;

         case Html_CODE:
            style.fFont = CWFont(style.fFont);
            PushStyleStack(Html_EndCODE, style);
            break;

         case Html_COMMENT:
            style.fFlags |= STY_Invisible;
            PushStyleStack(Html_EndCOMMENT, style);
            break;

         case Html_DD:
            if (fInnerList && fInnerList->fType == Html_DL) {
               ((TGHtmlRef *)p)->fPOther = fInnerList;
            } else {
               ((TGHtmlRef *)p)->fPOther = 0;
            }
            fInDt = 0;
            break;

         case Html_DFN:
            style.fFont = ItalicFont(style.fFont);
            PushStyleStack(Html_EndDFN, style);
            break;

         case Html_DIR:
         case Html_MENU:
         case Html_UL: {
            TGHtmlListStart *list = (TGHtmlListStart *) p;
            list->fLPrev = fInnerList;
            list->fCnt = 0;
            fInnerList = list;
            if (list->fLPrev == 0) {
               list->fLtype = LI_TYPE_Bullet1;
               list->fCompact = (list->MarkupArg("compact", 0) != 0);
            } else if (list->fLPrev->fLPrev == 0) {
               list->fLtype = LI_TYPE_Bullet2;
               list->fCompact = 1;
            } else {
               list->fLtype = LI_TYPE_Bullet3;
               list->fCompact = 1;
            }
            list->fLtype = list->GetUnorderedListType(list->fLtype);
            break;
         }

         case Html_EndDL:
            fInDt = 0;
            /* Fall thru into the next case */
         case Html_EndDIR:
         case Html_EndMENU:
         case Html_EndOL:
         case Html_EndUL:
            ((TGHtmlRef *)p)->fPOther = fInnerList;
            if (fInnerList) fInnerList = fInnerList->fLPrev;
            break;

         case Html_DIV:
            paraAlign = ALIGN_None;
            style.fAlign = p->GetAlignment(style.fAlign);
            PushStyleStack(Html_EndDIV, style);
            break;

         case Html_DT:
            if (fInnerList && fInnerList->fType == Html_DL) {
               ((TGHtmlRef *)p)->fPOther = fInnerList;
            } else {
               ((TGHtmlRef *)p)->fPOther = 0;
            }
            fInDt = STY_DT;
            break;

         case Html_EndDD:
         case Html_EndDT:
            fInDt = 0;
            break;

         case Html_DL: {
            TGHtmlListStart *list = (TGHtmlListStart *) p;
            list->fLPrev = fInnerList;
            list->fCnt = 0;
            fInnerList = list;
            list->fCompact = (list->MarkupArg("compact", 0) != 0);
            fInDt = 0;
            break;
         }

         case Html_EM:
            style.fFont = ItalicFont(style.fFont);
            PushStyleStack(Html_EndEM, style);
            break;

         case Html_EMBED:
            break;

         case Html_BASEFONT:
         case Html_FONT:
            z = p->MarkupArg("size", 0);
            if (z && !fOverrideFonts) {
               if (*z == '-') {
                  size = FontSize(style.fFont) - atoi(&z[1]) +1;
               } else if (*z == '+') {
                  size = FontSize(style.fFont) + atoi(&z[1]) +1;
               } else {
                  size = atoi(z);
               }
               if (size <= 0) size = 1;
               if (size >= N_FONT_SIZE) size = N_FONT_SIZE - 1;
               style.fFont = FontFamily(style.fFont) + size - 1;
            }
            z = p->MarkupArg("color", 0);
            if (z && *z && !fOverrideColors) style.fColor = GetColorByName(z);
            PushStyleStack(p->fType == Html_FONT ?
                           Html_EndFONT : Html_EndBASEFONT, style);
            break;

         case Html_FORM: {
            TGHtmlForm *form = (TGHtmlForm *) p;

            const char *zUrl;
            const char *zMethod;
            TGString cmd("");
            // int result;
            char zToken[50];

            fFormStart = 0;
            //form->fFormId = 0;

            zUrl = p->MarkupArg("action", 0);
            if (zUrl == 0) zUrl = fZBase;
            zUrl = ResolveUri(zUrl);
            if (zUrl == 0) zUrl = StrDup("");
            zMethod = p->MarkupArg("method", "GET");
            snprintf(zToken, 50, " %d form ", form->fFormId);
            cmd.Append("Form:");
            cmd.Append(zToken);
            cmd.Append(zUrl);
            cmd.Append(" ");
            cmd.Append(zMethod);
            cmd.Append(" { ");
            AppendArglist(&cmd, (TGHtmlMarkupElement *) p);
            cmd.Append("} ");
            /* result = */ FormCreate(form, zUrl, cmd.GetString());
            delete[] zUrl;

            /*if (result)*/ fFormStart = form;

            break;
         }

         case Html_EndFORM:
            ((TGHtmlRef *)p)->fPOther = fFormStart;
            if (fFormStart) fFormStart->fPEnd = p;
            fFormStart = 0;
            break;

         case Html_H1:
         case Html_H2:
         case Html_H3:
         case Html_H4:
         case Html_H5:
         case Html_H6:
            if (!fInTr) paraAlign = ALIGN_None;
            i = (p->fType - Html_H1) / 2 + 1;
            if (i >= 1 && i <= 6) {
               ScaleFont(&style, header_sizes[i-1]);
            }
            style.fFont = BoldFont(style.fFont);
            style.fAlign = p->GetAlignment(style.fAlign);
            PushStyleStack(Html_EndH1, style);
            break;

         case Html_EndH1:
         case Html_EndH2:
         case Html_EndH3:
         case Html_EndH4:
         case Html_EndH5:
         case Html_EndH6:
            paraAlign = ALIGN_None;
            style = PopStyleStack(Html_EndH1);
            break;

         case Html_HR:
            nextStyle = style;
            style.fAlign = p->GetAlignment(ALIGN_None);
            useNextStyle = 1;
            break;

         case Html_I:
            style.fFont = ItalicFont(style.fFont);
            PushStyleStack(Html_EndI, style);
            break;

         case Html_IMG:
            if (style.fFlags & STY_Invisible) break;
            ((TGHtmlImageMarkup *)p)->fPImage = GetImage((TGHtmlImageMarkup *) p);
            break;

         case Html_OPTION:
            break;

         case Html_INPUT:
            ((TGHtmlInput *)p)->fPForm = fFormStart;
            ////ControlSize((TGHtmlInput *) p);
            break;

         case Html_KBD:
            style.fFont = CWFont(style.fFont);
            PushStyleStack(Html_EndKBD, style);
            break;

         case Html_LI:
            if (fInnerList) {
               TGHtmlLi *li = (TGHtmlLi *) p;
               li->fLtype = fInnerList->fLtype;
               if (fInnerList->fType == Html_OL) {
                  z = li->MarkupArg("value", 0);
                  if (z) {
                     int n = atoi(z);
                     if (n > 0) {
                        li->fCnt = n;
                        fInnerList->fCnt = n+1;
                     }
                  } else {
                     li->fCnt = fInnerList->fCnt++;
                  }
                  li->fLtype = li->GetOrderedListType(li->fLtype);
               } else {
                  li->fLtype = li->GetUnorderedListType(li->fLtype);
               }
            } else {
               p->fFlags &= ~HTML_Visible;
            }
            break;

         case Html_MARQUEE:
            style.fFlags |= STY_Invisible;
            PushStyleStack(Html_EndMARQUEE, style);
            break;

         case Html_NOBR:
            style.fFlags |= STY_NoBreak;
            PushStyleStack(Html_EndNOBR, style);
            break;

         case Html_NOFRAMES:
            if (0 /* has ProcessFrame()*/) {
               nextStyle = style;
               nextStyle.fFlags |= STY_Invisible;
               PushStyleStack(Html_EndNOFRAMES, nextStyle);
               useNextStyle = 1;
            } else {
               PushStyleStack(Html_EndNOFRAMES, style);
            }
            break;

         case Html_NOEMBED:
            if (0 /* has ProcessScript() && HasScript */) {
               nextStyle = style;
               nextStyle.fFlags |= STY_Invisible;
               PushStyleStack(Html_EndNOEMBED, nextStyle);
               useNextStyle = 1;
            } else {
               PushStyleStack(Html_EndNOEMBED, style);
            }
            break;

         case Html_NOSCRIPT:
            if (0 /* has ProcessScript() && HasScript */) {
               nextStyle = style;
               nextStyle.fFlags |= STY_Invisible;
               PushStyleStack(Html_EndNOSCRIPT, nextStyle);
               useNextStyle = 1;
            } else {
               PushStyleStack(Html_EndNOSCRIPT, style);
            }
            break;

         case Html_OL: {
            TGHtmlListStart *list = (TGHtmlListStart *) p;
            list->fLPrev = fInnerList;
            list->fLtype = list->GetOrderedListType(LI_TYPE_Enum_1);
            list->fCnt = 1;
            z = list->MarkupArg("start", 0);
            if (z) {
               int n = atoi(z);
               if (n > 0) list->fCnt = n;
            }
            list->fCompact = (fInnerList != 0 || list->MarkupArg("compact", 0) != 0);
            fInnerList = list;
            break;
         }

         case Html_P:
            paraAlign = p->GetAlignment(ALIGN_None);
            break;

         case Html_EndP:
            paraAlign = ALIGN_None;
            break;

         case Html_PRE:
         case Html_LISTING:
         case Html_XMP:
         case Html_PLAINTEXT:
            paraAlign = ALIGN_None;
            style.fFont = CWFont(style.fFont);
            style.fFlags |= STY_Preformatted;
            PushStyleStack(Html_EndPRE, style);
            break;

         case Html_EndPRE:
         case Html_EndLISTING:
         case Html_EndXMP:
            style = PopStyleStack(Html_EndPRE);
            break;

         case Html_S:
            style.fFlags |= STY_StrikeThru;
            PushStyleStack(Html_EndS, style);
            break;

         case Html_SCRIPT: {
            char *result;
            result = ProcessScript((TGHtmlScript *) p);   // fZText[script->nStart .. script->nScript]
            if (result) {
               TGHtmlElement *b2 = p->fPNext, *b3, *e1 = p, *e2 = b2, *e3;
               if (e2) while (e2->fPNext) e2 = e2->fPNext;
               TokenizerAppend(result);
               if (e2 && e2 != p && ((e3 = b3 = e2->fPNext))) {
                  while (e3->fPNext) e3 = e3->fPNext;
                  e1->fPNext = b3;
                  e2->fPNext = 0;   b2->fPPrev = e3;
                  e3->fPNext = b2;  b3->fPPrev = e1;
               }
               delete[] result;
            }
            nextStyle = style;
            style.fFlags |= STY_Invisible;
            useNextStyle = 1;
            break;
         }

         case Html_SELECT:
            ((TGHtmlInput *)p)->fPForm = fFormStart;
            nextStyle.fFlags |= STY_Invisible;
            useNextStyle = 1;
            PushStyleStack(Html_EndSELECT, style);
            fFormElemStart = (TGHtmlInput *) p;
            break;

         case Html_EndSELECT:
            style = PopStyleStack(Html_EndSELECT);
            if (fFormElemStart && fFormElemStart->fType == Html_SELECT) {
                ((TGHtmlRef *)p)->fPOther = fFormElemStart;
               MakeInvisible(((TGHtmlRef *)p)->fPOther, p);
            } else {
               ((TGHtmlRef *)p)->fPOther = 0;
            }
            fFormElemStart = 0;
            break;

         case Html_STRIKE:
            style.fFlags |= STY_StrikeThru;
            PushStyleStack(Html_EndSTRIKE, style);
            break;

         case Html_STYLE:
            // Ignore style sheets
            break;

         case Html_SAMP:
            style.fFont = CWFont(style.fFont);
            PushStyleStack(Html_EndSAMP, style);
            break;

         case Html_SMALL:
            ScaleFont(&style, -1);
            PushStyleStack(Html_EndSMALL, style);
            break;

         case Html_STRONG:
            style.fFont = BoldFont(style.fFont);
            PushStyleStack(Html_EndSTRONG, style);
            break;

         case Html_SUB:
            ScaleFont(&style, -1);
            if (style.fSubscript > -6 ) style.fSubscript--;
            PushStyleStack(Html_EndSUB, style);
            break;

         case Html_SUP:
            ScaleFont(&style, -1);
            if (style.fSubscript < 6) style.fSubscript++;
            PushStyleStack(Html_EndSUP, style);
            break;

         case Html_TABLE:
            paraAlign = ALIGN_None;
            nextStyle = style;
            if (style.fFlags & STY_Preformatted) {
               nextStyle.fFlags &= ~STY_Preformatted;
               style.fFlags |= STY_Preformatted;
            }
            nextStyle.fAlign = ALIGN_Left;
            z = p->MarkupArg("bgcolor", 0);
            if (z && *z && !fOverrideColors) {
               style.fBgcolor = nextStyle.fBgcolor = GetColorByName(z);
               style.fExpbg = 1;
//            } else {
//               nextStyle.fBgcolor = COLOR_Background;
            }
            TableBgndImage(p);
            PushStyleStack(Html_EndTABLE, nextStyle);
            useNextStyle = 1;
            fInTd = 0;
            fInTr = 0;
            break;

         case Html_EndTABLE:
            paraAlign = ALIGN_None;
            if (fInTd) {
               style = PopStyleStack(Html_EndTD);
               fInTd = 0;
            }
            if (fInTr) {
               style = PopStyleStack(Html_EndTR);
               fInTr = 0;
            }
            style = PopStyleStack(p->fType);
            break;

         case Html_TD:
            if (fInTd) style = PopStyleStack(Html_EndTD);
            fInTd = 1;
            paraAlign = p->GetAlignment(rowAlign);
            z = p->MarkupArg("bgcolor", 0);
            if (z && *z && !fOverrideColors) {
               style.fBgcolor = GetColorByName(z);
               style.fExpbg = 1;
            }
            TableBgndImage(p);
            PushStyleStack(Html_EndTD, style);
            break;

         case Html_TEXTAREA:
            ((TGHtmlInput *)p)->fPForm = fFormStart;
            nextStyle = style;
            nextStyle.fFlags |= STY_Invisible;
            PushStyleStack(Html_EndTEXTAREA, nextStyle);
            fFormElemStart = (TGHtmlInput *) p;
            useNextStyle = 1;
            break;

         case Html_EndTEXTAREA:
            style = PopStyleStack(Html_EndTEXTAREA);
            if (fFormElemStart && fFormElemStart->fType == Html_TEXTAREA) {
               ((TGHtmlRef *)p)->fPOther = fFormElemStart;
            } else {
               ((TGHtmlRef *)p)->fPOther = 0;
            }
            fFormElemStart = 0;
            break;

         case Html_TH:
            //paraAlign = p->GetAlignment(rowAlign);
            if (fInTd) style = PopStyleStack(Html_EndTD);
            paraAlign = p->GetAlignment(ALIGN_Center);
            style.fFont = BoldFont(style.fFont);
            z = p->MarkupArg("bgcolor", 0);
            if (z && *z && !fOverrideColors) {
               style.fBgcolor = GetColorByName(z);
               style.fExpbg = 1;
            }
            PushStyleStack(Html_EndTD, style);
            fInTd = 1;
            break;

         case Html_TR:
            if (fInTd) {
               style = PopStyleStack(Html_EndTD);
               fInTd = 0;
            }
            if (fInTr) {
               style = PopStyleStack(Html_EndTR);
            }
            rowAlign = p->GetAlignment(ALIGN_None);
            z = p->MarkupArg("bgcolor", 0);
            if (z && *z && !fOverrideColors) {
               style.fBgcolor = GetColorByName(z);
               style.fExpbg = 1;
            }
            TableBgndImage(p);
            PushStyleStack(Html_EndTR, style);
            fInTr = 1;
            break;

         case Html_EndTR:
            if (fInTd) {
               style = PopStyleStack(Html_EndTD);
               fInTd = 0;
            }
            style = PopStyleStack(Html_EndTR);
            fInTr = 0;
            paraAlign = ALIGN_None;
            rowAlign = ALIGN_None;
            break;

         case Html_EndTD:
         case Html_EndTH:
            style = PopStyleStack(Html_EndTD);
            fInTd = 0;
            paraAlign = ALIGN_None;
            //rowAlign = ALIGN_None;
            break;

         case Html_TITLE:
            style.fFlags |= STY_Invisible;
            PushStyleStack(Html_EndTITLE, style);
            break;

         case Html_TT:
            style.fFont = CWFont(style.fFont);
            PushStyleStack(Html_EndTT, style);
            break;

         case Html_U:
            style.fFlags |= STY_Underline;
            PushStyleStack(Html_EndU, style);
            break;

         case Html_VAR:
            style.fFont = ItalicFont(style.fFont);
            PushStyleStack(Html_EndVAR, style);
            break;

         default:
            break;
      }

      p->fStyle = style;
      p->fStyle.fFlags |= fAnchorFlags | fInDt;
      if (paraAlign != ALIGN_None) {
         p->fStyle.fAlign = paraAlign;
      }
      if (useNextStyle) {
         style = nextStyle;
         style.fExpbg = 0;
         useNextStyle = 0;
      }

      TRACE(HtmlTrace_Style,
          ("Style font=%02d color=%02d bg=%02d "
           "align=%d flags=0x%04x token=%s\n",
           p->fStyle.fFont, p->fStyle.fColor, p->fStyle.fBgcolor,
           p->fStyle.fAlign, p->fStyle.fFlags, DumpToken(p)));

      p = p->fPNext;
   }

   // Copy state information back into the TGHtml object for safe keeping.

   fParaAlignment = paraAlign;
   fRowAlignment = rowAlign;

   fFlags &= ~STYLER_RUNNING;
}

////////////////////////////////////////////////////////////////////////////////
/// Set background picture of a html table.

void TGHtml::TableBgndImage(TGHtmlElement *p)
{
   const char *z;

   z = p->MarkupArg("background", 0);
   if (!z) return;

   char *z1 = ResolveUri(z);
   TImage *img = LoadImage(z1, 0, 0);
   delete [] z1;

   switch (p->fType) {
      case Html_TABLE: {
         TGHtmlTable *table = (TGHtmlTable *) p;
         if (table->fBgImage) delete table->fBgImage;
         table->fBgImage = img;
         break;
      }
      case Html_TR: {
         TGHtmlRef *ref = (TGHtmlRef *) p;
         if (ref->fBgImage) delete ref->fBgImage;
            ref->fBgImage = img;
            break;
      }
      case Html_TH:
      case Html_TD: {
         TGHtmlCell *cell = (TGHtmlCell *) p;
         if (cell->fBgImage) delete cell->fBgImage;
            cell->fBgImage = img;
            break;
         }
      default:
         if (img) delete img;
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the size of all elements in the widget. Assume that a style has
/// already been assigned to all elements.
///
/// Some of the elements might have already been sized. Refer to the
/// fLastSized and only compute sizes for elements that follow this one. If
/// fLastSized is 0, then size everything.
///
/// This routine only computes the sizes of individual elements. The size of
/// aggregate elements (like tables) are computed separately.
///
/// The HTML_Visible flag is also set on every element that results in ink on
/// the page.
///
/// This routine may invoke a callback procedure which could delete the HTML
/// widget.

void TGHtml::Sizer()
{
   TGHtmlElement *p;
   int iFont = -1;
   TGFont *font=0;
   int spaceWidth = 0;
   FontMetrics_t fontMetrics;
   const char *z;
   int stop = 0;

   if (fPFirst == 0) return;

   if (fLastSized == 0) {
      p = fPFirst;
   } else {
      p = fLastSized->fPNext;
   }

   // coverity[dead_error_line]
   for (; !stop && p; p = p ? p->fPNext : 0) {
      if (p->fStyle.fFlags & STY_Invisible) {
         p->fFlags &= ~HTML_Visible;
         continue;
      }
      if (iFont != (int)p->fStyle.fFont) {
         iFont = p->fStyle.fFont;
         font = GetFont(iFont);
         font->GetFontMetrics(&fontMetrics);
         spaceWidth = 0;
      }
      switch (p->fType) {
         case Html_Text: {
            TGHtmlTextElement *text = (TGHtmlTextElement *) p;
            text->fW = font->TextWidth(text->fZText, p->fCount);
            p->fFlags |= HTML_Visible;
            text->fDescent = fontMetrics.fDescent;
            text->fAscent = fontMetrics.fAscent;
            if (spaceWidth == 0) spaceWidth = font->TextWidth(" ", 1);
            text->fSpaceWidth = spaceWidth;
            break;
         }

         case Html_Space: {
            TGHtmlSpaceElement *space = (TGHtmlSpaceElement *) p;
            if (spaceWidth == 0) spaceWidth = font->TextWidth(" ", 1);
            space->fW = spaceWidth;
            space->fDescent = fontMetrics.fDescent;
            space->fAscent = fontMetrics.fAscent;
            p->fFlags &= ~HTML_Visible;
            break;
         }

         case Html_TD:
         case Html_TH: {
            TGHtmlCell *cell = (TGHtmlCell *) p;
            z = p->MarkupArg("rowspan", "1");
            cell->fRowspan = z ? atoi(z) : 1;
            z = p->MarkupArg("colspan", "1");
            cell->fColspan = z ? atoi(z) : 1;
            p->fFlags |= HTML_Visible;
            break;
         }

         case Html_LI: {
            TGHtmlLi *li = (TGHtmlLi *) p;
            li->fDescent = fontMetrics.fDescent;
            li->fAscent = fontMetrics.fAscent;
            p->fFlags |= HTML_Visible;
            break;
         }

         case Html_IMG: {
            TGHtmlImageMarkup *image = (TGHtmlImageMarkup *) p;
            z = p->MarkupArg("usemap", 0);
            if (z && *z == '#') {
               image->fPMap = GetMap(z+1);
            } else {
               image->fPMap = 0;
            }
            p->fFlags |= HTML_Visible;
            image->fRedrawNeeded = 0;
            image->fTextAscent = fontMetrics.fAscent;
            image->fTextDescent = fontMetrics.fDescent;
            image->fAlign = GetImageAlignment(p);
            if (image->fPImage == 0) {
               image->fAscent = fontMetrics.fAscent;
               image->fDescent = fontMetrics.fDescent;
               image->fZAlt = p->MarkupArg("alt", "<image>");
               if (image->fZAlt == 0) image->fZAlt = "<image>";
               image->fW = font->TextWidth(image->fZAlt, strlen(image->fZAlt));
            } else {
               int w, h;
               image->fINext = image->fPImage->fPList;
               image->fPImage->fPList = image;
               w = image->fPImage->fImage->GetWidth();
               h = image->fPImage->fImage->GetHeight();
               image->fH = h;
               image->fW = w;
               image->fAscent = h / 2;
               image->fDescent = h - image->fAscent;
            }
            if ((z = p->MarkupArg("width", 0)) != 0) {
               int w = atoi(z);
               if (z[strlen(z)-1] == '%') w = 0; //// -- HP
               if (w > 0) image->fW = w;
            }
            if ((z = p->MarkupArg("height", 0)) != 0) {
               int h = atoi(z);
               if (h > 0) image->fH = h;
            }

#if 1  // --HP
            if (image->fPImage == 0 && !*image->fZAlt) {
               image->fAscent = image->fH / 2;
               image->fDescent = image->fH - image->fAscent;
            }
#endif
            break;
         }

         case Html_TABLE:
            p->fFlags |= HTML_Visible;
            break;

         case Html_HR:
            p->fFlags |= HTML_Visible;
            break;

         case Html_APPLET:
         case Html_EMBED:
         case Html_INPUT: {
            TGHtmlInput *input = (TGHtmlInput *) p;
            input->fTextAscent = fontMetrics.fAscent;
            input->fTextDescent = fontMetrics.fDescent;
            stop = ControlSize(input);
            break;
         }

         case Html_SELECT:
         case Html_TEXTAREA: {
            TGHtmlInput *input = (TGHtmlInput *) p;
            input->fTextAscent = fontMetrics.fAscent;
            input->fTextDescent = fontMetrics.fDescent;
            break;
         }

         case Html_EndSELECT:
         case Html_EndTEXTAREA: {
            TGHtmlRef *ref = (TGHtmlRef *) p;
            if (ref->fPOther) {
               ((TGHtmlInput *)ref->fPOther)->fPEnd = p;
               stop = ControlSize((TGHtmlInput *) ref->fPOther);
            }
            break;
         }

         default:
            p->fFlags &= ~HTML_Visible;
            break;
      }
   }

   if (p) {
      fLastSized = p;
   } else {
      fLastSized = fPLast;
   }
}
