// @(#)root/gui:$Id$
// Author: Valeriy Onuchin & Fons Rademakers   16/10/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// The TGButtonGroup widget organizes TGButton widgets in a group.      //
//                                                                      //
// A button group widget makes it easier to deal with groups of buttons.//
// A button in a button group is associated with a unique identifier.   //
// The button group emits a Clicked() signal with this identifier when  //
// the button is clicked. Thus, a button group is an ideal solution     //
// when you have several similar buttons and want to connect all their  //
// Clicked() signals, for example, to one slot.                         //
//                                                                      //
// An exclusive button group switches off all toggle buttons except     //
// the one that was clicked. A button group is by default non-exclusive.//
// All radio buttons that are inserted, will be mutually exclusive even //
// if the button group is non-exclusive.                                //
//                                                                      //
//                                                                      //
// There are two ways of using a button group:                          //
//                                                                      //
//    The button group is a parent widget of a number of buttons,       //
//    i.e. the button group is the parent argument in the button        //
//    constructor. The buttons are assigned identifiers 1, 2, 3 etc.    //
//    in the order they are created or you can specify button id in     //
//    the button constructor. A TGButtonGroup can display a frame and   //
//    a title because it inherits from TGGroupFrame.                    //
//                                                                      //
// Example:                                                             //
//                                                                      //
//    // vertical frame without border and title                        //
//    TGVButtonGroup *bg = new TGVButtonGroup(main_frame);              //
//                                                                      //
//    // create text button with id=1                                   //
//    TGTextButton *button1 = new TGTextButton(bg,"some text");         //
//                                                                      //
//    // create another text button with id=2                           //
//    TGTextButton *button2 = new TGTextButton(bg,"another text");      //
//                                                                      //
//    // map all buttons                                                //
//    bg->Show();                                                       //
//                                                                      //
// NOTE: there is no need to call AddFrame() since the buttons are      //
// automatically added with a default layout hint to their parent,      //
// i.e. the buttongroup. To override the default layout hints use the   //
// SetLayoutHints() method.                                             //
//                                                                      //
//  ButtonGroup Signals:                                                //
//                                                                      //
//    Pressed(Int_t id)  -->  is emitted when a button in the group is  //
//                            pressed down. The id argument is the      //
//                            button's identifier.                      //
//    Released(Int_t id) -->  is emitted when a button in the group is  //
//                            released. The id argument is the button's //
//                            identifier.                               //
//    Clicked(Int_t id)  -->  is emitted when a button in the group is  //
//                            clicked. The id argument is the button's  //
//                            identifier.                               //
//                                                                      //
//                                                                      //
// The TGHButtonGroup widget organizes TGButton widgets in a group      //
// with one horizontal row. TGHButtonGroup is a convenience class that  //
// offers a thin layer on top of TGButtonGroup. It inherits from        //
// TGButtonGroup.                                                       //
//                                                                      //
// The TGVButtonGroup widget organizes TGButton widgets in a group      //
// with one vertical column. TGVButtonGroup is a convenience class that //
// offers a thin layer on top of TGButtonGroup. It inherits from        //
// TGButtonGroup.                                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGButtonGroup.h"
#include "TGButton.h"
#include "TClass.h"
#include "TGLayout.h"
#include "TList.h"
#include "TGResourcePool.h"
#include "TVirtualX.h"

#include <iostream>


ClassImp(TGButtonGroup);
ClassImp(TGHButtonGroup);
ClassImp(TGVButtonGroup);

////////////////////////////////////////////////////////////////////////////////
/// Constructor. Layout 1 row or 1 column.

TGButtonGroup::TGButtonGroup(const TGWindow *parent,
                             const TString &title,
                             UInt_t options,
                             GContext_t norm,
                             FontStruct_t font,
                             ULong_t back) :
   TGGroupFrame(parent, new TGString(title), options, norm, font, back)
{
   Init();
   if (options & kVerticalFrame) {
      SetLayoutManager(new TGVerticalLayout(this));
   } else {
      SetLayoutManager(new TGHorizontalLayout(this));
   }

   fDrawBorder = !title.IsNull();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor. Layout defined by TGMatrixLayout:
///    r = number of rows
///    c = number of columns
///    s = interval between frames
///    h = layout hints

TGButtonGroup::TGButtonGroup(const TGWindow *parent,
                             UInt_t r, UInt_t c,
                             Int_t s, Int_t h,
                             const TString &title,
                             GContext_t norm ,
                             FontStruct_t font ,
                             ULong_t back) :
   TGGroupFrame(parent, new TGString(title), 0, norm, font, back)
{
   Init();
   fDrawBorder = !title.IsNull();
   SetLayoutManager(new TGMatrixLayout(this,r,c,s,h));
}

////////////////////////////////////////////////////////////////////////////////
/// Default init.

void TGButtonGroup::Init()
{
   fState        = kTRUE;
   fMapOfButtons = new TMap();  // map of button/id pairs
   fExclGroup    = kFALSE;
   fRadioExcl    = kFALSE;
   fDrawBorder   = kTRUE;

   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor, we do not delete the buttons.

TGButtonGroup::~TGButtonGroup()
{
   TIter next(fMapOfButtons);
   TGButton *item = 0;

   while ((item = (TGButton*)next())) {
      item->SetGroup(0);
   }

   SafeDelete(fMapOfButtons);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw the group frame. Need special DoRedraw() since we need to
/// redraw with fBorderWidth=0.

void TGButtonGroup::DoRedraw()
{
   gVirtualX->ClearArea(fId, 0, 0, fWidth, fHeight);

   DrawBorder();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw border of around the group frame.
///
/// if frame is kRaisedFrame  - a frame border is of "wall style",
/// otherwise of "groove style".

void TGButtonGroup::DrawBorder()
{
   if (!fDrawBorder) return;

   Int_t x, y, l, t, r, b, gl, gr, sep, max_ascent, max_descent;

   UInt_t tw = gVirtualX->TextWidth(fFontStruct, fText->GetString(), fText->GetLength());
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);

   l = 0;
   t = (max_ascent + max_descent + 2) >> 1;
   r = fWidth - 1;
   // next three lines are for backward compatibility in case of horizontal layout
   TGLayoutManager * lm = GetLayoutManager();
   // coverity[returned_null]
   // coverity[dereference]
   if ((lm->InheritsFrom(TGHorizontalLayout::Class())) ||
       (lm->InheritsFrom(TGMatrixLayout::Class())))
      b = fHeight - 1;
   else
      b = fHeight - t;

   sep = 3;
   UInt_t rr = 5 + (sep << 1) + tw;

   switch (fTitlePos) {
      case kRight:
         gl = fWidth>rr ? Int_t(fWidth - rr) : 5 + sep;
         break;
      case kCenter:
         gl = fWidth>tw ? Int_t((fWidth - tw)>>1) - sep : 5 + sep;
         break;
      case kLeft:
      default:
         gl = 5 + sep;
   }
   gr = gl + tw + (sep << 1);

   switch (fOptions & (kSunkenFrame | kRaisedFrame)) {
      case kRaisedFrame:
         gVirtualX->DrawLine(fId, GetHilightGC()(),  l,   t,   gl,  t);
         gVirtualX->DrawLine(fId, GetShadowGC()(), l+1, t+1, gl,  t+1);

         gVirtualX->DrawLine(fId, GetHilightGC()(),  gr,  t,   r-1, t);
         gVirtualX->DrawLine(fId, GetShadowGC()(), gr,  t+1, r-2, t+1);

         gVirtualX->DrawLine(fId, GetHilightGC()(),  r-1, t,   r-1, b-1);
         gVirtualX->DrawLine(fId, GetShadowGC()(), r,   t,   r,   b);

         gVirtualX->DrawLine(fId, GetHilightGC()(),  r-1, b-1, l,   b-1);
         gVirtualX->DrawLine(fId, GetShadowGC()(), r,   b,   l,   b);

         gVirtualX->DrawLine(fId, GetHilightGC()(),  l,   b-1, l,   t);
         gVirtualX->DrawLine(fId, GetShadowGC()(), l+1, b-2, l+1, t+1);
         break;
      case kSunkenFrame:
      default:
         gVirtualX->DrawLine(fId, GetShadowGC()(),  l,   t,   gl,  t);
         gVirtualX->DrawLine(fId, GetHilightGC()(), l+1, t+1, gl,  t+1);

         gVirtualX->DrawLine(fId, GetShadowGC()(),  gr,  t,   r-1, t);
         gVirtualX->DrawLine(fId, GetHilightGC()(), gr,  t+1, r-2, t+1);

         gVirtualX->DrawLine(fId, GetShadowGC()(),  r-1, t,   r-1, b-1);
         gVirtualX->DrawLine(fId, GetHilightGC()(), r,   t,   r,   b);

         gVirtualX->DrawLine(fId, GetShadowGC()(),  r-1, b-1, l,   b-1);
         gVirtualX->DrawLine(fId, GetHilightGC()(), r,   b,   l,   b);

         gVirtualX->DrawLine(fId, GetShadowGC()(),  l,   b-1, l,   t);
         gVirtualX->DrawLine(fId, GetHilightGC()(), l+1, b-2, l+1, t+1);
         break;
   }

   x = gl + sep;
   y = 1;

   if (fState) {
      fText->Draw(fId, fNormGC, x, y + max_ascent);
   } else {
      fText->Draw(fId, GetHilightGC()(), x, y + 1 + max_ascent);
      fText->Draw(fId, GetShadowGC()(), x, y + max_ascent);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Makes border to be visible/invisible.

void TGButtonGroup::SetBorderDrawn(Bool_t enable)
{
   if (enable != IsBorderDrawn()) {
      fDrawBorder = enable;
      ChangedBy("SetBorderDrawn");        // emit signal
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the button group to be exclusive if enable is kTRUE,
/// or to be non-exclusive if enable is kFALSE.
/// An exclusive button group switches off all other toggle buttons when
/// one is switched on. This is ideal for groups of radio-buttons
/// A non-exclusive group allow many buttons to be switched on at the same
/// time. The default setting is kFALSE.

void TGButtonGroup::SetExclusive(Bool_t enable)
{
   if (enable != IsExclusive()) {
      fExclGroup = enable;
      ChangedBy("SetExclusive");  // emit signal
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If enable is kTRUE, this button group will treat radio buttons as
/// mutually exclusive, and other buttons according to IsExclusive().
/// This function is called automatically whenever a TGRadioButton
/// is inserted, so you should normally never have to call it.

void TGButtonGroup::SetRadioButtonExclusive(Bool_t enable)
{
   if (enable != IsRadioButtonExclusive()) {
      fRadioExcl = enable;
      ChangedBy("SetRadioButtonExclusive"); // emit signal
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the state of all the buttons in the group to enable or disable.

void TGButtonGroup::SetState(Bool_t state)
{
   fState = state;

   TIter next(fMapOfButtons);
   TGButton *item = 0;

   while ((item = (TGButton*)next())) {    // loop over all buttons
      if (state) {
         item->SetState(kButtonUp);
      } else {
         item->SetState(kButtonDisabled);
      }
   }
   DoRedraw();
}
////////////////////////////////////////////////////////////////////////////////
/// Sets the button with id to be on/down, and if this is an
/// exclusive group, all other button in the group to be off/up.

void TGButtonGroup::SetButton(Int_t id, Bool_t down)
{
   TGButton *b = Find(id);

   if (b && (b->IsDown() != down)) {
      b->SetState(kButtonDown, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Inserts a button with the identifier id into the button group.
/// Returns the button identifier.
///
/// It is not necessary to manually insert buttons that have this button
/// group as their parent widget. An exception is when you want custom
/// identifiers instead of the default 1, 2, 3 etc.
///
/// The button is assigned the identifier id or an automatically
/// generated identifier.  It works as follows: If id > 0, this
/// identifier is assigned.  If id == -1 (default), the identifier is
/// equal to the number of buttons in the group+1.  If id is any other
/// negative integer, for instance -2, a unique identifier (negative
/// integer <= -2) is generated.
///
/// Inserting several buttons with id = -1 assigns the identifiers 1,
/// 2, 3, etc.

Int_t TGButtonGroup::Insert(TGButton *button, Int_t id)
{
   if (button->fGroup && button->fGroup != this)
      button->fGroup->Remove(button);

   if (button->fGroup == this) {
      if (id == -1)
         return GetId(button);    // the button is already in group
      else
         button->fGroup->Remove(button);  // want to set a new id
   }

   button->fGroup = this;
   button->Associate(this);

   static Int_t seq_no = -2;
   Long_t bid;

   if (id < -1)       bid = seq_no--;
   else if (id == -1) bid = GetCount()+1;
   else               bid = id;

   fMapOfButtons->Add(button, (TObject*)bid);
   AddFrame(button);

   // coverity[returned_null]
   // coverity[dereference]
   SetRadioButtonExclusive(button->IsA()->InheritsFrom(TGRadioButton::Class()));

   Connect(button, "Clicked()" , "TGButtonGroup", this, "ReleaseButtons()");
   Connect(button, "Pressed()" , "TGButtonGroup", this, "ButtonPressed()");
   Connect(button, "Released()", "TGButtonGroup", this, "ButtonReleased()");
   Connect(button, "Clicked()" , "TGButtonGroup", this, "ButtonClicked()");

   return (Int_t) bid;
}

////////////////////////////////////////////////////////////////////////////////
/// Removes a button from the button group.

void TGButtonGroup::Remove(TGButton *button)
{
   TGButton *item = (TGButton*) fMapOfButtons->Remove(button);
   if (item) {
      button->SetGroup(0);
      button->Disconnect(this);
      button->DestroyWindow();
   }

   RemoveFrame(button);
}

////////////////////////////////////////////////////////////////////////////////
/// Finds and returns a pointer to the button with the specified
/// identifier id. Returns null if the button was not found.

TGButton *TGButtonGroup::Find(Int_t id) const
{
   TIter next(fMapOfButtons);
   TGButton *item = 0;

   while ((item = (TGButton*)next())) {
      if ((Long_t)fMapOfButtons->GetValue(item) == id) break;   // found
   }

   return item;
}

////////////////////////////////////////////////////////////////////////////////
/// Finds and returns the id of the button.
/// Returns -1 if the button is not a member of this group.

Int_t TGButtonGroup::GetId(TGButton *button) const
{
   TPair *a = (TPair*) fMapOfButtons->FindObject(button);
   if (a)
      return (Int_t)Long_t(a->Value());
   else
      return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// This slot is activated when one of the buttons in the group emits the
/// Pressed() signal.

void TGButtonGroup::ButtonPressed()
{
#if 0
   // Is here for historical purposes and example. Now this is not needed
   // anymore since TGButton has has its own GetSender() method returning
   // the TGButton proper.

   // This is needed since gTQSender points to TQObject part of TGButton
   TGButton *btn = dynamic_cast<TGButton*>((TQObject*)gTQSender);

   if (!btn) {
      Error("ButtonPressed", "gTQSender not a TGButton");
      return;
   }
#else
      TGButton *btn = (TGButton*)gTQSender;
#endif

   TPair *a = (TPair*) fMapOfButtons->FindObject(btn);
   if (a) {
      Int_t id = (Int_t)Long_t(a->Value());
      Pressed(id);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// This slot is activated when one of the buttons in the group emits the
/// Released() signal.

void TGButtonGroup::ButtonReleased()
{
   TGButton *btn = (TGButton*)gTQSender;

   TPair *a = (TPair*) fMapOfButtons->FindObject(btn);
   if (a) {
      Int_t id = (Int_t)Long_t(a->Value());
      Released(id);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// This slot is activated when one of the buttons in the group emits the
/// Clicked() signal.

void TGButtonGroup::ButtonClicked()
{
   TGButton *btn = (TGButton*)gTQSender;

   TPair *a = (TPair*) fMapOfButtons->FindObject(btn);
   if (a) {
      Int_t id = (Int_t)Long_t(a->Value());
      Clicked(id);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// This slot is activated when one of the buttons in the
/// exclusive group emits the Pressed() signal.

void TGButtonGroup::ReleaseButtons()
{
   if (!fExclGroup && !fRadioExcl) return;

   TGButton *btn = (TGButton*)gTQSender;

   if (!fExclGroup && !btn)
      return;

   TIter next(fMapOfButtons);
   TGButton *item = 0;

   while ((item = (TGButton*)next())) {    // loop over all buttons
      // coverity[returned_null]
      // coverity[dereference]
      if (btn != item && item->IsToggleButton() && item->IsOn() &&
          (fExclGroup || (item->IsA()->InheritsFrom(TGRadioButton::Class())
                          && btn->IsA()->InheritsFrom(TGRadioButton::Class())))) {
         item->SetOn(kFALSE);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Show group of buttons.

void TGButtonGroup::Show()
{
   MapSubwindows();
   Resize();
   MapRaised();
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Hide group of buttons.

void TGButtonGroup::Hide()
{
   UnmapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Set or change title.

void TGButtonGroup::SetTitle(TGString *title)
{
   if (!title) {
      Error("SetTitle", "title cannot be 0, try \"\"");
      return;
   }

   if (strcmp(fText->GetString(), title->GetString())) {
      SetBorderDrawn(title->GetLength() ? kTRUE : kFALSE);
      TGGroupFrame::SetTitle(title);
      ChangedBy("SetTitle");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set or change title.

void TGButtonGroup::SetTitle(const char *title)
{
   if (!title) {
      Error("SetTitle", "title cannot be 0, try \"\"");
      return;
   }

   if (strcmp(fText->GetString(), title)) {
      SetBorderDrawn(title && strlen(title));
      TGGroupFrame::SetTitle(title);
      ChangedBy("SetTitle");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set layout hints for the specified button or if button=0 for all
/// buttons.

void TGButtonGroup::SetLayoutHints(TGLayoutHints *l, TGButton *button)
{
   TGFrameElement *el;
   TIter next(fList);

   while ((el = (TGFrameElement *)next())) {
      if ((el->fFrame==(TGFrame*)button) || !button) {
         el->fLayout = l ? l : fgDefaultHints;
      }
   }
   Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Save a button group widget as a C++ statement(s) on output stream out.

void TGButtonGroup::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote ='"';

   // font + GC
   option = GetName()+5;         // unique digit id of the name
   TString parGC, parFont;
   // coverity[returned_null]
   // coverity[dereference]
   parFont.Form("%s::GetDefaultFontStruct()",IsA()->GetName());
   // coverity[returned_null]
   // coverity[dereference]
   parGC.Form("%s::GetDefaultGC()()",IsA()->GetName());

   if ((GetDefaultFontStruct() != fFontStruct) || (GetDefaultGC()() != fNormGC)) {
      TGFont *ufont = gClient->GetResourcePool()->GetFontPool()->FindFont(fFontStruct);
      if (ufont) {
         ufont->SavePrimitive(out, option);
         parFont.Form("ufont->GetFontStruct()");
      }

      TGGC *userGC = gClient->GetResourcePool()->GetGCPool()->FindGC(fNormGC);
      if (userGC) {
         userGC->SavePrimitive(out, option);
         parGC.Form("uGC->GetGC()");
      }
   }

   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   out << std::endl << "   // buttongroup frame" << std::endl;

   out << "   TGButtonGroup *";
   out << GetName() << " = new TGButtonGroup(" << fParent->GetName()
       << ","<< quote << fText->GetString() << quote;

   if (fBackground == GetDefaultFrameBackground()) {
      if (fFontStruct == GetDefaultFontStruct()) {
         if (fNormGC == GetDefaultGC()()) {
            if (!GetOptions()) {
               out <<");" << std::endl;
            } else {
               out << "," << GetOptionString() <<");" << std::endl;
            }
         } else {
            out << "," << GetOptionString() << "," << parGC.Data() <<");" << std::endl;
         }
      } else {
         out << "," << GetOptionString() << "," << parGC.Data() << "," << parFont.Data() <<");" << std::endl;
      }
   } else {
      out << "," << GetOptionString() << "," << parGC.Data() << "," << parFont.Data() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   // setting layout manager
   out << "   " << GetName() <<"->SetLayoutManager(";
   // coverity[returned_null]
   // coverity[dereference]
   GetLayoutManager()->SavePrimitive(out, option);
   out << ");"<< std::endl;

   TGFrameElement *f;
   TIter next(GetList());
   while ((f = (TGFrameElement *)next())) {
      f->fFrame->SavePrimitive(out,option);
      if (f->fFrame->InheritsFrom("TGButton")) continue;
      else {
         out << "   " << GetName() << "->AddFrame(" << f->fFrame->GetName();
         f->fLayout->SavePrimitive(out, option);
         out << ");"<< std::endl;
      }
   }

   if (IsExclusive())
      out << "   " << GetName() <<"->SetExclusive(kTRUE);" << std::endl;

   if (IsRadioButtonExclusive())
      out << "   " << GetName() <<"->SetRadioButtonExclusive(kTRUE);" << std::endl;

   if (!IsBorderDrawn())
      out << "   " << GetName() <<"->SetBorderDrawn(kFALSE);" << std::endl;


   out << "   " << GetName() << "->Resize(" << GetWidth()
       << "," << GetHeight() << ");" << std::endl;

   if (!IsEnabled())
      out << "   " << GetName() <<"->SetState(kFALSE);" << std::endl;

   out << "   " << GetName() << "->Show();" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a button group widget as a C++ statement(s) on output stream out.

void TGHButtonGroup::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote ='"';

   // font + GC
   option = GetName()+5;         // unique digit id of the name
   TString parGC, parFont;
   parFont.Form("%s::GetDefaultFontStruct()",IsA()->GetName());
   parGC.Form("%s::GetDefaultGC()()",IsA()->GetName());

   if ((GetDefaultFontStruct() != fFontStruct) || (GetDefaultGC()() != fNormGC)) {
      TGFont *ufont = gClient->GetResourcePool()->GetFontPool()->FindFont(fFontStruct);
      if (ufont) {
         ufont->SavePrimitive(out, option);
         parFont.Form("ufont->GetFontStruct()");
      }

      TGGC *userGC = gClient->GetResourcePool()->GetGCPool()->FindGC(fNormGC);
      if (userGC) {
         userGC->SavePrimitive(out, option);
         parGC.Form("uGC->GetGC()");
      }
   }

   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   out << std::endl << "   // horizontal buttongroup frame" << std::endl;

   out << "   TGHButtonGroup *";
   out << GetName() << " = new TGHButtonGroup(" << fParent->GetName()
       << "," << quote << fText->GetString() << quote;
   if (fBackground == GetDefaultFrameBackground()) {

      if (fFontStruct == GetDefaultFontStruct()) {

         if (fNormGC == GetDefaultGC()()) {
            out << ");" << std::endl;
         } else {
            out << "," << parGC.Data() <<");" << std::endl;
         }
      } else {
         out << "," << parGC.Data() << "," << parFont.Data() <<");" << std::endl;
      }
   } else {
      out << "," << parGC.Data() << "," << parFont.Data() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   TGFrameElement *f;
   TIter next(GetList());
   while ((f = (TGFrameElement *)next())) {
      f->fFrame->SavePrimitive(out,option);
      if (f->fFrame->InheritsFrom("TGButton")){
         out << "   " << GetName() << "->SetLayoutHints(";
         f->fLayout->SavePrimitive(out, "nocoma");
         out << "," << f->fFrame->GetName();
         out << ");"<< std::endl;
      }
      else {
         out << "   " << GetName() << "->AddFrame(" << f->fFrame->GetName();
         f->fLayout->SavePrimitive(out, option);
         out << ");"<< std::endl;
      }
   }

   if (!IsEnabled())
      out << "   " << GetName() <<"->SetState(kFALSE);" << std::endl;

   if (IsExclusive())
      out << "   " << GetName() <<"->SetExclusive(kTRUE);" << std::endl;

   if (IsRadioButtonExclusive())
      out << "   " << GetName() <<"->SetRadioButtonExclusive(kTRUE);" << std::endl;

   if (!IsBorderDrawn())
      out << "   " << GetName() <<"->SetBorderDrawn(kFALSE);" << std::endl;

   out << "   " << GetName() <<"->Resize(" << GetWidth() << ","
       << GetHeight() << ");" << std::endl;

   out << "   " << GetName() << "->Show();" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a button group widget as a C++ statement(s) on output stream out.

void TGVButtonGroup::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote ='"';

   // font + GC
   option = GetName()+5;         // unique digit id of the name
   TString parGC, parFont;
   parFont.Form("%s::GetDefaultFontStruct()",IsA()->GetName());
   parGC.Form("%s::GetDefaultGC()()",IsA()->GetName());

   if ((GetDefaultFontStruct() != fFontStruct) || (GetDefaultGC()() != fNormGC)) {
      TGFont *ufont = gClient->GetResourcePool()->GetFontPool()->FindFont(fFontStruct);
      if (ufont) {
         ufont->SavePrimitive(out, option);
         parFont.Form("ufont->GetFontStruct()");
      }

      TGGC *userGC = gClient->GetResourcePool()->GetGCPool()->FindGC(fNormGC);
      if (userGC) {
         userGC->SavePrimitive(out, option);
         parGC.Form("uGC->GetGC()");
      }
   }

   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   out << std::endl << "   // vertical buttongroup frame" << std::endl;

   out << "   TGVButtonGroup *";
   out << GetName() << " = new TGVButtonGroup(" << fParent->GetName()
       << "," << quote << fText->GetString() << quote;

   if (fBackground == GetDefaultFrameBackground()) {
      if (fFontStruct == GetDefaultFontStruct()) {
         if (fNormGC == GetDefaultGC()()) {
            out <<");" << std::endl;
         } else {
            out << "," << parGC.Data() <<");" << std::endl;
         }
      } else {
         out << "," << parGC.Data() << "," << parFont.Data() <<");" << std::endl;
      }
   } else {
      out << "," << parGC.Data() << "," << parFont.Data() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   TGFrameElement *f;
   TIter next(GetList());
   while ((f = (TGFrameElement *)next())) {
      f->fFrame->SavePrimitive(out,option);
      if (f->fFrame->InheritsFrom("TGButton")) continue;
      else {
         out << "   " << GetName() << "->AddFrame(" << f->fFrame->GetName();
         f->fLayout->SavePrimitive(out, option);
         out << ");"<< std::endl;
      }
   }

   if (!IsEnabled())
      out << "   " << GetName() <<"->SetState(kFALSE);" << std::endl;

   if (IsExclusive())
      out << "   " << GetName() <<"->SetExclusive(kTRUE);" << std::endl;

   if (IsRadioButtonExclusive())
      out << "   " << GetName() <<"->SetRadioButtonExclusive(kTRUE);" << std::endl;

   if (!IsBorderDrawn())
      out << "   " << GetName() <<"->SetBorderDrawn(kFALSE);" << std::endl;

   out << "   " << GetName() << "->Resize(" << GetWidth()
       << "," << GetHeight() << ");"<< std::endl;

   out << "   " << GetName() << "->Show();" << std::endl;
}
