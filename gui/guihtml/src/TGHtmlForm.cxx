// $Id: TGHtmlForm.cxx,v 1.3 2007/05/18 16:00:28 brun Exp $
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

// Routines used for processing HTML makeup for forms.

#include <cstring>
#include <cstdlib>
#include <cstdarg>

#include "TGHtml.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGListBox.h"
#include "TGTextEdit.h"
#include "TGComboBox.h"
#include "snprintf.h"

////////////////////////////////////////////////////////////////////////////////
/// Unmap any input control that is currently mapped.

void TGHtml::UnmapControls()
{
   TGHtmlInput *p;

   for (p = fFirstInput; p; p = p->fINext) {
      if (p->fFrame != 0 /*&& p->fFrame->IsMapped()*/) {
         p->fFrame->UnmapWindow();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Map any control that should be visible according to the
/// current scroll position. At the same time, if any controls that
/// should not be visible are mapped, unmap them. After this routine
/// finishes, all <INPUT> controls should be in their proper places
/// regardless of where they might have been before.
///
/// Return the number of controls that are currently visible.

int TGHtml::MapControls()
{
   TGHtmlInput *p;     // For looping over all controls
   int x, y, w, h;    // Part of the virtual canvas that is visible
   int cnt = 0;       // Number of visible controls

   x = fVisible.fX;
   y = fVisible.fY;
   w = fCanvas->GetWidth();
   h = fCanvas->GetHeight();
   for (p = fFirstInput; p; p = p->fINext) {
      if (p->fFrame == 0) continue;
      if (p->fY < y + h && p->fY + p->fH > y &&
          p->fX < x + w && p->fX + p->fW > x) {
         // The control should be visible. Make is so if it isn't already
         p->fFrame->MoveResize(p->fX - x, p->fY + fFormPadding/2 - y,
                             p->fW, p->fH - fFormPadding);
         /*if (!p->fFrame->IsMapped())*/ p->fFrame->MapWindow();
         ++cnt;
      } else {
         // This control should not be visible. Unmap it.
         /*if (p->fFrame->IsMapped())*/ p->fFrame->UnmapWindow();
      }
   }

   return cnt;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete all input controls. This happens when the TGHtml widget
/// is cleared.

void TGHtml::DeleteControls()
{
   TGHtmlInput *p;        // For looping over all controls

   p = fFirstInput;
   fFirstInput = 0;
   fLastInput = 0;
   fNInput = 0;

   if (p == 0) return;

   for (; p; p = p->fINext) {
      if (p->fPForm && ((TGHtmlForm *)p->fPForm)->fHasctl) {
         ((TGHtmlForm *)p->fPForm)->fHasctl = 0;
      }
      if (p->fFrame) {
         if (!fExiting) p->fFrame->DestroyWindow();
         delete p->fFrame;
         p->fFrame = 0;
      }
      p->fSized = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return an appropriate type value for the given <INPUT> markup.

static int InputType(TGHtmlElement *p)
{
   int type = INPUT_TYPE_Unknown;
   const char *z;
   int i;
   static struct {
      const char *zName;
      int type;
   } types[] = {
      { "checkbox",  INPUT_TYPE_Checkbox },
      { "file",      INPUT_TYPE_File     },
      { "hidden",    INPUT_TYPE_Hidden   },
      { "image",     INPUT_TYPE_Image    },
      { "password",  INPUT_TYPE_Password },
      { "radio",     INPUT_TYPE_Radio    },
      { "reset",     INPUT_TYPE_Reset    },
      { "submit",    INPUT_TYPE_Submit   },
      { "text",      INPUT_TYPE_Text     },
      { "name",      INPUT_TYPE_Text     },
      { "textfield", INPUT_TYPE_Text     },
      { "button",    INPUT_TYPE_Button   },
      { "name",      INPUT_TYPE_Text     },
   };

   switch (p->fType) {
      case Html_INPUT:
         z = p->MarkupArg("type", "text");
         if (z == 0) break;
         for (i = 0; i < int(sizeof(types) / sizeof(types[0])); i++) {
            if (strcasecmp(types[i].zName, z) == 0) {
               type = types[i].type;
               break;
            }
         }
         break;

      case Html_SELECT:
         type = INPUT_TYPE_Select;
         break;

      case Html_TEXTAREA:
         type = INPUT_TYPE_TextArea;
         break;

      case Html_APPLET:
      case Html_IFRAME:
      case Html_EMBED:
         type = INPUT_TYPE_Applet;
         break;

      default:
         CANT_HAPPEN;
         break;
   }
   return type;
}

////////////////////////////////////////////////////////////////////////////////
/// 'frame' is the child widget that is used to implement an input
/// element. Query the widget for its size and put that information
/// in the pElem structure that represents the input.

void TGHtml::SizeAndLink(TGFrame *frame, TGHtmlInput *pElem)
{

   pElem->fFrame = frame;
   if (pElem->fFrame == 0) {
      pElem->Empty();
   } else if (pElem->fItype == INPUT_TYPE_Hidden) {
      pElem->fW = 0;
      pElem->fH = 0;
      pElem->fFlags &= ~HTML_Visible;
      pElem->fStyle.fFlags |= STY_Invisible;
   } else {
      pElem->fW = frame->GetDefaultWidth();
      pElem->fH = frame->GetDefaultHeight() + fFormPadding;
      pElem->fFlags |= HTML_Visible;
      pElem->fHtml = this;
   }
   pElem->fINext = 0;
   if (fFirstInput == 0) {
      fFirstInput = pElem;
   } else {
      fLastInput->fINext = pElem;
   }
   fLastInput = pElem;
   pElem->fSized = 1;

#if 0
   if (pElem->fFrame) {
      pElem->fFrame->ChangeOptions(pElem->fFrame->GetOptions() | kOwnBackground);
      pElem->fFrame->SetBackgroundColor(_defaultFrameBackground);
   }
#else
   if (pElem->fFrame) {
      int bg = pElem->fStyle.fBgcolor;
      //int fg = pElem->fStyle.color;
      ColorStruct_t *cbg = fApColor[bg];
      //ColorStruct_t *cfg = fApColor[fg];
      pElem->fFrame->ChangeOptions(pElem->fFrame->GetOptions() | kOwnBackground);
      pElem->fFrame->SetBackgroundColor(cbg->fPixel);
   }
#endif

   if (pElem->fFrame) {
      // the following is needed by some embedded widgets like
      // TGListBox and TGTextEdit
      pElem->fFrame->MapSubwindows();
      pElem->fFrame->Layout();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Append all text and space tokens between pStart and pEnd to
/// the given TString.  [ TGTextEdit ]

void TGHtml::AppendText(TGString *str, TGHtmlElement *pFirs,
                        TGHtmlElement *pEnd)
{
   while (pFirs && pFirs != pEnd) {
      switch (pFirs->fType) {
         case Html_Text:
            str->Append(((TGHtmlTextElement *)pFirs)->fZText);
            break;

         case Html_Space:
            if (pFirs->fFlags & HTML_NewLine) {
               str->Append("\n");
            } else {
               int cnt;
               static char zSpaces[] = "                             ";
               cnt = pFirs->fCount;
               while (cnt > (int)sizeof(zSpaces) - 1) {
                  str->Append(zSpaces, sizeof(zSpaces) - 1);
                  cnt -= sizeof(zSpaces) - 1;
               }
               if (cnt > 0) {
                  str->Append(zSpaces, cnt);
               }
            }
            break;

         default:
            break;
      }
      pFirs = pFirs->fPNext;
   }
}


class TGHtmlLBEntry : public TGTextLBEntry {
public:
   TGHtmlLBEntry(const TGWindow *p, TGString *s, TGString *val, int ID) :
      TGTextLBEntry(p, s, ID) { fVal = val; }
   virtual ~TGHtmlLBEntry() { if (fVal) delete fVal; }

   const char *GetValue() const { return fVal ? fVal->GetString() : 0; }

protected:
   TGString *fVal;
};


////////////////////////////////////////////////////////////////////////////////
/// The "p" argument points to a <select>.  This routine scans all
/// subsequent elements (up to the next </select>) looking for
/// <option> tags.  For each option tag, it appends the corresponding
/// entry to the "lb" listbox element.
///
/// lb   -- An TGListBox object
/// p    -- The <SELECT> markup
/// pEnd -- The </SELECT> markup

void TGHtml::AddSelectOptions(TGListBox *lb, TGHtmlElement *p,
                              TGHtmlElement *pEnd)
{
   int id = 0;

   while (p && p != pEnd && p->fType != Html_EndSELECT) {
      if (p->fType == Html_OPTION) {
         TGString *str;
         int selected = -1;

         const char *zValue = p->MarkupArg("value", "");
         const char *sel = p->MarkupArg("selected", "");
         if (sel && !strcmp(sel, "selected"))
            selected = id;

         p = p->fPNext;

         str = new TGString("");
         while (p && p != pEnd &&
                p->fType != Html_EndOPTION &&
                p->fType != Html_OPTION &&
                p->fType != Html_EndSELECT) {
            if (p->fType == Html_Text) {
               str->Append(((TGHtmlTextElement *)p)->fZText);
            } else if (p->fType == Html_Space) {
               str->Append(" ");
            }
            p = p->fPNext;
         }
         lb->AddEntry(new TGHtmlLBEntry(lb->GetContainer(), str,
                      new TGString(zValue), id),
                      new TGLayoutHints(kLHintsTop | kLHintsExpandX));
         //if (p->MarkupArg("selected", 0) != 0) lb->Select(id);
         if (selected >= 0)
            lb->Select(selected);
         ++id;
      } else {
         p = p->fPNext;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// This routine implements the Sizer() function for <INPUT>,
/// <SELECT> and <TEXTAREA> markup.
///
/// A side effect of sizing these markups is that widgets are
/// created to represent the corresponding input controls.
///
/// The function normally returns 0.  But if it is dealing with
/// a <SELECT> or <TEXTAREA> that is incomplete, 1 is returned.
/// In that case, the sizer will be called again at some point in
/// the future when more information is available.

int TGHtml::ControlSize(TGHtmlInput *pElem)
{
   int incomplete = 0;    // kTRUE if data is incomplete

   if (pElem->fSized) return 0;

   pElem->fItype = InputType(pElem);   //// pElem->InputType();
                                     //// or done in the constructor!

//   if (pElem->fPForm == 0) {
//      pElem->Empty();
//      return incomplete;
//   }

   switch (pElem->fItype) {
      case INPUT_TYPE_File:
      case INPUT_TYPE_Hidden:
      case INPUT_TYPE_Image:
         pElem->Empty();
         SizeAndLink(0, pElem);
         break;

      case INPUT_TYPE_Checkbox: {
         pElem->fCnt = ++fNInput;
         TGCheckButton *f = new TGCheckButton(fCanvas, "", pElem->fCnt);
         if (pElem->MarkupArg("checked", 0))
            ((TGCheckButton *)f)->SetState(kButtonDown);
         f->Associate(this);
         f->Resize(f->GetDefaultSize());
         SizeAndLink(f, pElem);
         break;
      }

      case INPUT_TYPE_Radio: {
         pElem->fCnt = ++fNInput;
         TGRadioButton *f = new TGRadioButton(fCanvas, "", pElem->fCnt);
         if (pElem->MarkupArg("checked", 0))
            ((TGRadioButton *)f)->SetState(kButtonDown);
         f->Associate(this);
         f->Resize(f->GetDefaultSize());
         SizeAndLink(f, pElem);
         break;
      }

      case INPUT_TYPE_Reset: {
         pElem->fCnt = ++fNInput;
         const char *z = pElem->MarkupArg("value", 0);
         if (!z) z = "Reset";
         TGTextButton *f = new TGTextButton(fCanvas, new TGHotString(z), pElem->fCnt);
         f->RequestFocus();
         f->Associate(this);
         f->Resize(f->GetDefaultSize());
         SizeAndLink(f, pElem);
         break;
      }

      case INPUT_TYPE_Button:
      case INPUT_TYPE_Submit: {
         pElem->fCnt = ++fNInput;
         const char *z = pElem->MarkupArg("value", 0);
         if (!z) z = "Submit";
         TGTextButton *f = new TGTextButton(fCanvas, new TGHotString(z), pElem->fCnt);
         f->RequestFocus();
         f->Associate(this);
         // TODO: bg color!
         f->Resize(f->GetDefaultSize());
         SizeAndLink(f, pElem);
         break;
      }

      case INPUT_TYPE_Text: {
         pElem->fCnt = ++fNInput;
         const char *z = pElem->MarkupArg("maxlength", 0);
         int maxlen = z ? atoi(z) : 256;
         if (maxlen < 2) maxlen = 2;
         z = pElem->MarkupArg("size", 0);
         int size = z ? atoi(z) * 5 : 150;
         TGTextEntry *f = new TGTextEntry(fCanvas, new TGTextBuffer(maxlen),
                                          pElem->fCnt);
         z = pElem->MarkupArg("value", 0);
         if (z) f->AppendText(z);
         f->Resize(size, f->GetDefaultHeight());
         SizeAndLink(f, pElem);
         break;
      }

      case INPUT_TYPE_Password: {
         pElem->fCnt = ++fNInput;
         const char *z = pElem->MarkupArg("maxlength", 0);
         int maxlen = z ? atoi(z) : 256;
         if (maxlen < 2) maxlen = 2;
         z = pElem->MarkupArg("size", 0);
         int size = z ? atoi(z) * 5 : 150;
         TGTextEntry *f = new TGTextEntry(fCanvas, new TGTextBuffer(maxlen),
                                          pElem->fCnt);
         f->SetEchoMode(TGTextEntry::kPassword);
         z = pElem->MarkupArg("value", 0);
         if (z) f->AppendText(z);
         f->Resize(size, f->GetDefaultHeight());
         SizeAndLink(f, pElem);
         break;
      }

      case INPUT_TYPE_Select: {  // listbox or dd-listbox?
         pElem->fCnt = ++fNInput;
         const char *z = pElem->MarkupArg("size", 0);
         int size = z ? atoi(z) : 1;
         UInt_t width = 0, height = 0;
         if (size == 1) {
            TGComboBox *cb = new TGComboBox(fCanvas, pElem->fCnt);
            TGListBox *lb = cb->GetListBox();
            AddSelectOptions(lb, pElem, pElem->fPEnd);
            TGLBEntry *e = lb->GetSelectedEntry();
            if (e) lb->Select(e->EntryId(), kFALSE);
            lb->MapSubwindows();
            lb->Layout();
            for (int i=0;i<lb->GetNumberOfEntries();++i) {
               TGHtmlLBEntry *te = (TGHtmlLBEntry *)lb->GetEntry(i);
               if (te && te->GetText())
                  width = TMath::Max(width, te->GetDefaultWidth());
            }
            height = lb->GetItemVsize() ? lb->GetItemVsize()+4 : 22;
            cb->Resize(width > 0 ? width+30 : 200,
                       height > 22 ? height : 22);
            if (e) cb->Select(e->EntryId(), kFALSE);
            SizeAndLink(cb, pElem);
         } else {
            TGListBox *lb = new TGListBox(fCanvas, pElem->fCnt);
            z = pElem->MarkupArg("multiple", 0);
            if (z) lb->SetMultipleSelections(kTRUE);
            AddSelectOptions(lb, pElem, pElem->fPEnd);
            for (int i=0;i<lb->GetNumberOfEntries();++i) {
               TGHtmlLBEntry *te = (TGHtmlLBEntry *)lb->GetEntry(i);
               if (te && te->GetText())
                  width = TMath::Max(width, te->GetDefaultWidth());
            }
            height = lb->GetItemVsize() ? lb->GetItemVsize() : 22;
            lb->Resize(width > 0 ? width+30 : 200, height * size);
            lb->Associate(this);
            SizeAndLink(lb, pElem);
         }
         break;
      }

      case INPUT_TYPE_TextArea: {
         pElem->fCnt = ++fNInput;
         // const char *z = pElem->MarkupArg("rows", 0);
         //int rows = z ? atoi(z) : 10;
         // coverity[returned_pointer]
         // z = pElem->MarkupArg("cols", 0);
         //int cols = z ? atoi(z) : 10;
         TGTextEdit *f = new TGTextEdit(fCanvas, 300, 200, pElem->fCnt);
         TGString str("");
         AppendText(&str, pElem, pElem->fPEnd);
         //f->InsertText(&str);
         SizeAndLink(f, pElem);
         break;
      }

      case INPUT_TYPE_Applet: {
         //int result;

         TGFrame *f = ProcessApplet(pElem);
         if (!f) {
            pElem->Empty();
            break;
         }
         pElem->fCnt = ++fNInput;
         SizeAndLink(f, pElem);
         break;
      }

      default: {
         CANT_HAPPEN;
         pElem->fFlags &= ~HTML_Visible;
         pElem->fStyle.fFlags |= STY_Invisible;
         pElem->fFrame = 0;
         break;
      }
   }
   return incomplete;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of elments of type p in a form.

int TGHtml::FormCount(TGHtmlInput *p, int radio)
{
   TGHtmlElement *q = p;

   switch (p->fType) {
      case Html_SELECT:
         return p->fSubId;
      case Html_TEXTAREA:
      case Html_INPUT:
         if (radio && p->fType == INPUT_TYPE_Radio)
            return p->fSubId;
         return ((TGHtmlForm *)p->fPForm)->fElements;
      case Html_OPTION:
         while ((q = q->fPPrev))
            if (q->fType == Html_SELECT) return ((TGHtmlInput *)q)->fSubId;
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Add the DOM control information for form elements.

void TGHtml::AddFormInfo(TGHtmlElement *p)
{
   TGHtmlElement *q;
   TGHtmlForm *f;
   const char *name, *z;
   int t;

   switch (p->fType) {
      case Html_SELECT:
      case Html_TEXTAREA:
      case Html_INPUT: {
         TGHtmlInput *input = (TGHtmlInput *) p;
         if (!(f = fFormStart)) return;
         input->fPForm = fFormStart;
         if (!f->fPFirst)
            f->fPFirst = p;
         if (fFormElemLast)
            fFormElemLast->fINext = input;
         fFormElemLast = input;
         input->fInpId = fInputIdx++;
         t = input->fItype = InputType(input);
         if (t == INPUT_TYPE_Radio) {
            if ((name = p->MarkupArg("name", 0))) {
               for (q = f->fPFirst; q; q = ((TGHtmlInput *)q)->fINext) {
                  if ((z = q->MarkupArg("name", 0)) && !strcmp(z, name)) {
                     input->fSubId = fRadioIdx++;
                     break;
                  }
               }
               if (!q) input->fSubId = fRadioIdx = 0;
            }
         }
         break;
      }

      case Html_FORM:
         fFormStart = (TGHtmlForm *) p;
         ((TGHtmlForm *)p)->fFormId = fNForm++;
         break;

      case Html_EndTEXTAREA:
      case Html_EndSELECT:
      case Html_EndFORM:
         fFormStart = 0;
         fInputIdx = 0;
         fRadioIdx = 0;
         fFormElemLast = 0;
         break;

      case Html_OPTION:
         if (fFormElemLast && fFormElemLast->fType == Html_SELECT)
            fFormElemLast->fSubId++;
         break;

      default:
         break;
   }
}

// The following array determines which characters can be put directly
// in a query string and which must be escaped.

static char gNeedEscape[] = {
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1,
};
#define NeedToEscape(C) ((C)>0 && (C)<127 && gNeedEscape[(int)(C)])

////////////////////////////////////////////////////////////////////////////////
/// Append to the given TString an encoded version of the given text.

void TGHtml::EncodeText(TGString *str, const char *z)
{
   int i;

   while (*z) {
      for (i = 0; z[i] && !NeedToEscape(z[i]); ++i) {}
      if (i > 0) str->Append(z, i);
      z += i;
      while (*z && NeedToEscape(*z)) {
         if (*z == ' ') {
            str->Append("+", 1);
         } else if (*z == '\n') {
            str->Append("%0D%0A", 6);
         } else if (*z == '\r') {
            // Ignore it...
         } else {
            char zBuf[10];
            snprintf(zBuf, 10, "%%%02X", 0xff & *z);
            str->Append(zBuf, 3);
         }
         z++;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Process messages (GUI events) in the html widget.

Bool_t TGHtml::ProcessMessage(Long_t msg, Long_t p1, Long_t p2)
{
/*
  OWidgetMessage *wmsg = (OWidgetMessage *) msg;
  TGHtmlInput *p;

  switch (msg->fType) {
    case MSG_BUTTON:
    case MSG_RADIOBUTTON:
    case MSG_CHECKBUTTON:
    case MSG_LISTBOX:
    case MSG_DDLISTBOX:
      for (p = fFirstInput; p; p = p->fINext) {
        if (p->fCnt == wmsg->id) {
          switch (p->fItype) {
            case INPUT_TYPE_Button:
            case INPUT_TYPE_Submit:
              if (p->fPForm) {
                FormAction(p->fPForm, wmsg->id);
              } else {
                printf("action, but no form!\n");
              }
              break;

            case INPUT_TYPE_Reset: {
              //ClearForm(p->fPForm);
              TGHtmlInput *pr;
              for (pr = fFirstInput; pr; pr = pr->fINext) {
                if (pr->fPForm == p->fPForm) {
                  switch (pr->fItype) {
                    case INPUT_TYPE_Radio: {
                      TGRadioButton *rb = (TGRadioButton *) pr->fFrame;
                      if (pr->MarkupArg("checked", 0))
                        rb->SetState(kButtonDown);
                      else
                        rb->SetState(kButtonUp);
                      break;
                    }

                    case INPUT_TYPE_Checkbox: {
                      TGCheckButton *cb = (TGCheckButton *) pr->fFrame;
                      if (pr->MarkupArg("checked", 0))
                        cb->SetState(kButtonDown);
                      else
                        cb->SetState(kButtonUp);
                      break;
                    }

                    case INPUT_TYPE_Text:
                    case INPUT_TYPE_Password: {
                      TGTextEntry *te = (TGTextEntry *) pr->fFrame;
                      te->Clear();
                      const char *z = pr->MarkupArg("value", 0);
                      if (z) te->AddText(0, z);
                      break;
                    }

                    case INPUT_TYPE_Select: {
                      break;
                    }

                    default:
                      break;
                  }
                }
              }
              break;
            }

            case INPUT_TYPE_Radio: {
              TGHtmlInput *pr;
              for (pr = fFirstInput; pr; pr = pr->fINext) {
                if ((pr->fPForm == p->fPForm) &&
                    (pr->fItype == INPUT_TYPE_Radio)) {
                  if (pr != p) {
                    if (strcmp(pr->MarkupArg("name", ""),
                               p->MarkupArg("name", "")) == 0)
                      ((TGRadioButton *)pr->fFrame)->SetState(kButtonUp);
                  }
                }
              }
              break;
            }

            case INPUT_TYPE_Select: {
              break;
            }

            default:
              break;
          }
          return kTRUE;
        }
      }
      break;

    default:
      break;
  }
*/
   return TGView::ProcessMessage(msg, p1, p2);
}
