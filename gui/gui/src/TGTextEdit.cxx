// @(#)root/gui:$Id$
// Author: Fons Rademakers   3/7/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
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
// TGTextEdit                                                           //
//                                                                      //
// A TGTextEdit is a specialization of TGTextView. It provides the      //
// text edit functionality to the static text viewing widget.           //
// For the messages supported by this widget see the TGView class.      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGTextEdit.h"
#include "TGTextEditDialogs.h"
#include "TGResourcePool.h"
#include "TSystem.h"
#include "TTimer.h"
#include "TGMenu.h"
#include "TGMsgBox.h"
#include "TGFileDialog.h"
#include "TGScrollBar.h"
#include "KeySymbols.h"
#include "RConfigure.h"
#include "TVirtualX.h"

#include <iostream>

static const char *gFiletypes[] = { "All files",     "*",
                                    "Text files",    "*.txt",
                                    "ROOT macros",   "*.C",
                                    0,               0 };
static char *gPrinter      = nullptr;
static char *gPrintCommand = nullptr;


TGGC *TGTextEdit::fgCursor0GC = nullptr;
TGGC *TGTextEdit::fgCursor1GC = nullptr;


///////////////////////////////////////////////////////////////////////////////
class TGTextEditHist : public TList {

public:
   TGTextEditHist() {}
   virtual ~TGTextEditHist() { Delete(); }

   Bool_t Notify() { //
      TObject *obj = Last();
      if (!obj) return kFALSE;

      obj->Notify(); // execute undo action
      RemoveLast();
      delete obj;
      return kTRUE;
   }
};

///////////////////////////////////////////////////////////////////////////////
class TGTextEditCommand : public TObject {
protected:
   TGTextEdit     *fEdit;
   TGLongPosition  fPos;

public:
   TGTextEditCommand(TGTextEdit *te) : fEdit(te) {
      fPos = fEdit->GetCurrentPos();
      fEdit->GetHistory()->Add(this);
   }
   void SetPos(TGLongPosition pos) { fPos = pos; }
};

///////////////////////////////////////////////////////////////////////////////
class TInsCharCom : public TGTextEditCommand {

public:
   TInsCharCom(TGTextEdit *te, char ch) : TGTextEditCommand(te) {
      fEdit->InsChar(ch);
   }
   Bool_t Notify() { //
      fEdit->SetCurrent(fPos);
      fEdit->NextChar();
      fEdit->DelChar();
      return kTRUE;
   }
};

///////////////////////////////////////////////////////////////////////////////
class TDelCharCom : public TGTextEditCommand {

private:
   char fChar;

public:
   TDelCharCom(TGTextEdit *te) : TGTextEditCommand(te)  {
      fPos.fX--;
      fChar = fEdit->GetText()->GetChar(fPos);
      fEdit->DelChar();
   }
   Bool_t Notify() { //
      if (fChar > 0) {
         fEdit->SetCurrent(fPos);
         fEdit->InsChar(fChar);
      } else {
         fPos.fY--;
         fEdit->BreakLine();
      }
      return kTRUE;
   }
};

///////////////////////////////////////////////////////////////////////////////
class TBreakLineCom : public TGTextEditCommand {

public:
   TBreakLineCom(TGTextEdit *te) : TGTextEditCommand(te)  {
      fEdit->BreakLine();
      fPos.fX = 0;
      fPos.fY++;
   }

   Bool_t Notify() { //
      fEdit->SetCurrent(fPos);
      fEdit->DelChar();
      return kTRUE;
   }
};

///////////////////////////////////////////////////////////////////////////////
class TInsTextCom : public TGTextEditCommand {
private:
   TGLongPosition  fEndPos;

public:
   char            fChar;

   TInsTextCom(TGTextEdit *te) : TGTextEditCommand(te), fChar(0)  {
   }

   void SetEndPos(TGLongPosition end) {
      fEndPos = end;
   }

   Bool_t Notify() { //
      fEdit->GetText()->DelText(fPos, fEndPos);

      if (fChar > 0) {
         fEdit->GetText()->InsChar(fPos, fChar);
      } else if (fPos.fY != fEndPos.fY) {
         fEdit->GetText()->BreakLine(fPos);
      }
      fEdit->SetCurrent(fPos);
      fEdit->Update();
      return kTRUE;
   }
};

///////////////////////////////////////////////////////////////////////////////
class TDelTextCom : public TGTextEditCommand {

private:
   TGText         *fText;
   TGLongPosition  fEndPos;
   Bool_t          fBreakLine;

public:
   TDelTextCom(TGTextEdit *te, TGText *txt) : TGTextEditCommand(te)  {
      fText = new TGText(txt);
      fBreakLine = kFALSE;
   }
   TDelTextCom(const TDelTextCom &dtc) : TGTextEditCommand(dtc) {
      fText = new TGText(dtc.fText);
      fBreakLine = dtc.fBreakLine;
   }
   virtual ~TDelTextCom() { delete fText; }

   TDelTextCom &operator=(const TDelTextCom &dtc) {
      if (this != &dtc) {
         if (fText) delete fText;
         fText = new TGText(dtc.fText);
         fBreakLine = dtc.fBreakLine;
      }
      return *this;
   }

   void SetEndPos(TGLongPosition end) {
      fEndPos = end;
   }

   void SetBreakLine(Bool_t on) { fBreakLine = on; }

   Bool_t Notify() { //
      TGLongPosition start_src, end_src;
      start_src.fX = start_src.fY = 0;
      end_src.fY   = fText->RowCount() - 1;
      end_src.fX   = fText->GetLineLength(end_src.fY) - 1;

      fEdit->GetText()->InsText(fPos, fText, start_src, end_src);

      if (fBreakLine) {
         fEndPos.fY++;
         fEdit->GetText()->BreakLine(fEndPos);
         fEndPos.fX = fEdit->GetText()->GetLineLength(fEndPos.fY);
      } else {
         fEndPos.fX++;
      }

      fEdit->SetCurrent(fEndPos);
      fEdit->Update();
      return kTRUE;
   }
};


ClassImp(TGTextEdit);


////////////////////////////////////////////////////////////////////////////////
/// Create a text edit widget.

TGTextEdit::TGTextEdit(const TGWindow *parent, UInt_t w, UInt_t h, Int_t id,
                       UInt_t sboptions, ULong_t back) :
     TGTextView(parent, w, h, id, sboptions, back)
{
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a text edit widget. Initialize it with the specified text buffer.

TGTextEdit::TGTextEdit(const TGWindow *parent, UInt_t w, UInt_t h, TGText *text,
                       Int_t id, UInt_t sboptions, ULong_t back) :
     TGTextView(parent, w, h, text, id, sboptions, back)
{
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a text edit widget. Initialize it with the specified string.

TGTextEdit::TGTextEdit(const TGWindow *parent, UInt_t w, UInt_t h,
                       const char *string, Int_t id, UInt_t sboptions,
                       ULong_t back) :
     TGTextView(parent, w, h, string, id, sboptions, back)
{
   Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup text edit widget.

TGTextEdit::~TGTextEdit()
{
   if (TGSearchDialog::SearchDialog()) {
      TQObject::Disconnect(TGSearchDialog::SearchDialog(), 0, this);
   }
   delete fCurBlink;
   delete fMenu;
   delete fHistory;
}

////////////////////////////////////////////////////////////////////////////////
/// Initiliaze a text edit widget.

void TGTextEdit::Init()
{
   fCursor0GC   = GetCursor0GC()();
   fCursor1GC   = GetCursor1GC()();
   fCursorState = 1;
   fCurrent.fY  = fCurrent.fX = 0;
   fInsertMode  = kInsert;
   fCurBlink    = 0;
   fSearch      = 0;
   fEnableMenu  = kTRUE;
   fEnableCursorWithoutFocus = kTRUE;

   gVirtualX->SetCursor(fCanvas->GetId(), fClient->GetResourcePool()->GetTextCursor());

   // create popup menu with default editor actions
   fMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fMenu->AddEntry("New", kM_FILE_NEW);
   fMenu->AddEntry("Open...", kM_FILE_OPEN);
   fMenu->AddSeparator();
   fMenu->AddEntry("Close", kM_FILE_CLOSE);
   fMenu->AddEntry("Save", kM_FILE_SAVE);
   fMenu->AddEntry("Save As...", kM_FILE_SAVEAS);
   fMenu->AddSeparator();
   fMenu->AddEntry("Print...", kM_FILE_PRINT);
   fMenu->AddSeparator();
   fMenu->AddEntry("Cut", kM_EDIT_CUT);
   fMenu->AddEntry("Copy", kM_EDIT_COPY);
   fMenu->AddEntry("Paste", kM_EDIT_PASTE);
   fMenu->AddEntry("Select All", kM_EDIT_SELECTALL);
   fMenu->AddSeparator();
   fMenu->AddEntry("Find...", kM_SEARCH_FIND);
   fMenu->AddEntry("Find Again", kM_SEARCH_FINDAGAIN);
   fMenu->AddEntry("Goto...", kM_SEARCH_GOTO);

   fMenu->Associate(this);

   fHistory = new TGTextEditHist();
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/disable menu items in function of what is possible.

void TGTextEdit::SetMenuState()
{
   if (fText->RowCount() == 1 && fText->GetLineLength(0) <= 0) {
      fMenu->DisableEntry(kM_FILE_CLOSE);
      fMenu->DisableEntry(kM_FILE_SAVE);
      fMenu->DisableEntry(kM_FILE_SAVEAS);
      fMenu->DisableEntry(kM_FILE_PRINT);
      fMenu->DisableEntry(kM_EDIT_SELECTALL);
      fMenu->DisableEntry(kM_SEARCH_FIND);
      fMenu->DisableEntry(kM_SEARCH_FINDAGAIN);
      fMenu->DisableEntry(kM_SEARCH_GOTO);
   } else {
      fMenu->EnableEntry(kM_FILE_CLOSE);
      fMenu->EnableEntry(kM_FILE_SAVE);
      fMenu->EnableEntry(kM_FILE_SAVEAS);
      fMenu->EnableEntry(kM_FILE_PRINT);
      fMenu->EnableEntry(kM_EDIT_SELECTALL);
      fMenu->EnableEntry(kM_SEARCH_FIND);
      fMenu->EnableEntry(kM_SEARCH_FINDAGAIN);
      fMenu->EnableEntry(kM_SEARCH_GOTO);
   }

   if (IsSaved())
      fMenu->DisableEntry(kM_FILE_SAVE);
   else
      fMenu->EnableEntry(kM_FILE_SAVE);

   if (fIsMarked) {
      fMenu->EnableEntry(kM_EDIT_CUT);
      fMenu->EnableEntry(kM_EDIT_COPY);
   } else {
      fMenu->DisableEntry(kM_EDIT_CUT);
      fMenu->DisableEntry(kM_EDIT_COPY);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return width of longest line in widget.

Long_t TGTextEdit::ReturnLongestLineWidth()
{
   Long_t linewidth = TGTextView::ReturnLongestLineWidth();
   linewidth += 3*fScrollVal.fX;
   return linewidth;
}

////////////////////////////////////////////////////////////////////////////////
/// Clear text edit widget.

void TGTextEdit::Clear(Option_t *)
{
   fCursorState = 1;
   fCurrent.fY = fCurrent.fX = 0;
   TGTextView::Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Save file. If filename==0 ask user via dialog for a filename, if in
/// addition saveas==kTRUE always ask for new filename. Returns
/// kTRUE if file was correctly saved, kFALSE otherwise.

Bool_t TGTextEdit::SaveFile(const char *filename, Bool_t saveas)
{
   if (!filename) {
      Bool_t untitled = !strlen(fText->GetFileName()) ? kTRUE : kFALSE;
      if (untitled || saveas) {
         static TString dir(".");
         static Bool_t overwr = kFALSE;
         TGFileInfo fi;
         fi.fFileTypes = gFiletypes;
         fi.SetIniDir(dir);
         fi.fOverwrite = overwr;
         new TGFileDialog(fClient->GetDefaultRoot(), this, kFDSave, &fi);
         overwr = fi.fOverwrite;
         if (fi.fFilename && strlen(fi.fFilename)) {
            dir = fi.fIniDir;
            return fText->Save(fi.fFilename);
         }
         return kFALSE;
      }
      return fText->Save(fText->GetFileName());
   }

   return fText->Save(filename);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy text.

Bool_t TGTextEdit::Copy()
{
   if (!fIsMarked || ((fMarkedStart.fX == fMarkedEnd.fX) &&
       (fMarkedStart.fY == fMarkedEnd.fY))) {
      return kFALSE;
   }

   TGTextView::Copy();

   Bool_t del = !fCurrent.fX && (fCurrent.fY == fMarkedEnd.fY) && !fMarkedEnd.fX;
   del = del || (!fMarkedEnd.fX && (fCurrent.fY != fMarkedEnd.fY));
   del = del && fClipText->AsString().Length() > 0;

   if (del) {
      TGLongPosition pos;
      pos.fY = fClipText->RowCount();
      pos.fX = 0;
      fClipText->InsText(pos, 0);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Cut text.

Bool_t TGTextEdit::Cut()
{
   if (!Copy()) {
      return kFALSE;
   }
   Delete();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Paste text into widget.

Bool_t TGTextEdit::Paste()
{
   if (fReadOnly) {
      return kFALSE;
   }

   if (fIsMarked) {
      TString sav = fClipText->AsString();
      TGTextView::Copy();
      Delete();
      fClipText->Clear();
      fClipText->LoadBuffer(sav.Data());
   }

   gVirtualX->ConvertPrimarySelection(fId, fClipboard, 0);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Send current buffer to printer.

void TGTextEdit::Print(Option_t *) const
{
   TString msg;

   msg.Form("%s -P%s\n", gPrintCommand, gPrinter);
   FILE *p = gSystem->OpenPipe(msg.Data(), "w");
   if (p) {
      char   *buf1, *buf2;
      Long_t  len;
      ULong_t i = 0;
      TGLongPosition pos;

      pos.fX = pos.fY = 0;
      while (pos.fY < fText->RowCount()) {
         len = fText->GetLineLength(pos.fY);
         if (len < 0) len = 0;
         buf1 = fText->GetLine(pos, len);
         buf2 = new char[len + 2];
         strncpy(buf2, buf1, (UInt_t)len);
         buf2[len]   = '\n';
         buf2[len+1] = '\0';
         while (buf2[i] != '\0') {
            if (buf2[i] == '\t') {
               ULong_t j = i+1;
               while (buf2[j] == 16)
                  j++;
               // coverity[secure_coding]
               strcpy(buf2+i+1, buf2+j);
            }
            i++;
         }
         fwrite(buf2, sizeof(char), strlen(buf2)+1, p);

         delete [] buf1;
         delete [] buf2;
         pos.fY++;
      }
      gSystem->ClosePipe(p);

      Bool_t untitled = !strlen(fText->GetFileName()) ? kTRUE : kFALSE;
      msg.Form("Printed: %s\nLines: %ld\nUsing: %s -P%s",
               untitled ? "Untitled" : fText->GetFileName(),
               fText->RowCount() - 1, gPrintCommand, gPrinter);
      new TGMsgBox(fClient->GetDefaultRoot(), this, "Editor", msg.Data(),
                   kMBIconAsterisk, kMBOk, 0);
   } else {
      msg.Form("Could not execute: %s -P%s\n", gPrintCommand, gPrinter);
      new TGMsgBox(fClient->GetDefaultRoot(), this, "Editor", msg.Data(),
                   kMBIconExclamation, kMBOk, 0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete selection.

void TGTextEdit::Delete(Option_t *)
{
   if (!fIsMarked || fReadOnly) {
      return;
   }

   if (fMarkedStart.fX == fMarkedEnd.fX &&
       fMarkedStart.fY == fMarkedEnd.fY) {
      Long_t len = fText->GetLineLength(fCurrent.fY);

      if (fCurrent.fY == fText->RowCount()-1 && fCurrent.fX == len) {
         gVirtualX->Bell(0);
         return;
      }

      new TDelCharCom(this);
      return;
   }

   TGLongPosition pos, endPos;
   Bool_t delast = kFALSE;

   endPos.fX = fMarkedEnd.fX - 1;
   endPos.fY = fMarkedEnd.fY;

   if (endPos.fX == -1) {
      pos = fCurrent;
      if (endPos.fY > 0) {
         SetCurrent(endPos);
         DelChar();
         endPos.fY--;
         SetCurrent(pos);
      }
      endPos.fX = fText->GetLineLength(endPos.fY);
      if (endPos.fX < 0) {
         endPos.fX = 0;
      }
      delast = kTRUE;
   }

   // delete command for undo
   TDelTextCom *dcom = new TDelTextCom(this, fClipText);
   dcom->SetPos(fMarkedStart);
   dcom->SetEndPos(endPos);

   if (delast || ((fText->GetLineLength(endPos.fY) == endPos.fX+1) &&
       (fClipText->RowCount() > 1))) {
      TGLongPosition p = endPos;

      p.fY--;
      if (!delast) p.fX++;
      dcom->SetEndPos(p);
      dcom->SetBreakLine(kTRUE);
   }

   fText->DelText(fMarkedStart, endPos);

   pos.fY = ToObjYCoord(fVisible.fY);

   if (fMarkedStart.fY < pos.fY) {
      pos.fY = fMarkedStart.fY;
   }
   pos.fX = ToObjXCoord(fVisible.fX, pos.fY);
   if (fMarkedStart.fX < pos.fX) {
      pos.fX = fMarkedStart.fX;
   }

   Int_t th = (Int_t)ToScrYCoord(fText->RowCount());
   Int_t ys = (Int_t)ToScrYCoord(fMarkedStart.fY);
   th = th < 0 ? 0 : th;
   ys = ys < 0 ? 0 : ys;

   // clear
   if ((th < 0) || (th < (Int_t)fCanvas->GetHeight())) {
      gVirtualX->ClearArea(fCanvas->GetId(), 0, ys,
                               fCanvas->GetWidth(), fCanvas->GetHeight() - ys);
   }

   UpdateRegion(0, ys, fCanvas->GetWidth(), UInt_t(fCanvas->GetHeight() - ys));

   SetVsbPosition((ToScrYCoord(pos.fY) + fVisible.fY)/fScrollVal.fY);
   SetHsbPosition((ToScrXCoord(pos.fX, pos.fY) + fVisible.fX)/fScrollVal.fX);
   SetSBRange(kHorizontal);
   SetSBRange(kVertical);
   SetCurrent(fMarkedStart);

   SendMessage(fMsgWindow, MK_MSG(kC_TEXTVIEW, kTXT_ISMARKED), fWidgetId, kFALSE);
   UnMark();

   // only to make sure that IsSaved() returns true in case everything has
   // been deleted
   if (fText->RowCount() == 1 && fText->GetLineLength(0) == 0) {
      delete fText;
      fText = new TGText();
      fText->Clear();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Search for string in the specified direction. If direction is true
/// the search will be in forward direction.

Bool_t TGTextEdit::Search(const char *string, Bool_t direction,
                          Bool_t caseSensitive)
{
   if (!IsMapped()) return kFALSE;

   if (gTQSender && (gTQSender == TGSearchDialog::SearchDialog())) {
      caseSensitive = TGSearchDialog::SearchDialog()->GetType()->fCaseSensitive;
      direction = TGSearchDialog::SearchDialog()->GetType()->fDirection;
      fSearch = TGSearchDialog::SearchDialog()->GetType();
   }

   TGLongPosition pos;
   if (!fText->Search(&pos, fCurrent, string, direction, caseSensitive)) {
      fCurrent.fX = 1;
      fCurrent.fY = 1;

      if (!fText->Search(&pos, fCurrent, string, direction, caseSensitive)) { //try again
         TString msg;
         msg.Form("Couldn't find \"%s\"", string);
         gVirtualX->Bell(20);
         new TGMsgBox(fClient->GetDefaultRoot(), fCanvas, "TextEdit",
                      msg.Data(), kMBIconExclamation, kMBOk, 0);
         return kFALSE;
      }
      return kTRUE;
   }
   UnMark();
   fIsMarked = kTRUE;
   fMarkedStart.fY = fMarkedEnd.fY = pos.fY;
   fMarkedStart.fX = pos.fX;
   fMarkedEnd.fX = fMarkedStart.fX + strlen(string);

   if (direction) {
      SetCurrent(fMarkedEnd);
   } else {
      SetCurrent(fMarkedStart);
   }

   pos.fY = ToObjYCoord(fVisible.fY);
   if (fCurrent.fY < pos.fY ||
       ToScrYCoord(fCurrent.fY) >= (Int_t)fCanvas->GetHeight()) {
      pos.fY = fMarkedStart.fY;
   }
   pos.fX = ToObjXCoord(fVisible.fX, pos.fY);

   if (fCurrent.fX < pos.fX ||
       ToScrXCoord(fCurrent.fX, pos.fY) >= (Int_t)fCanvas->GetWidth()) {
      pos.fX = fMarkedStart.fX;
   }

   SetVsbPosition((ToScrYCoord(pos.fY)+fVisible.fY)/fScrollVal.fY);
   SetHsbPosition((ToScrXCoord(pos.fX, pos.fY)+fVisible.fX)/fScrollVal.fX);

   UpdateRegion(0, (Int_t)ToScrYCoord(fMarkedStart.fY), fCanvas->GetWidth(),
                UInt_t(ToScrYCoord(fMarkedEnd.fY+1)-ToScrYCoord(fMarkedEnd.fY)));

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Replace text starting at textPos.

Bool_t TGTextEdit::Replace(TGLongPosition textPos, const char *oldText,
                           const char *newText, Bool_t direction, Bool_t caseSensitive)
{
   TGLongPosition pos;
   if (!fText->Replace(textPos, oldText, newText, direction, caseSensitive)) {
      return kFALSE;
   }
   UnMark();
   fIsMarked = kTRUE;
   fMarkedStart.fY = fMarkedEnd.fY = textPos.fY;
   fMarkedStart.fX = textPos.fX;
   fMarkedEnd.fX = fMarkedStart.fX + strlen(newText);

   if (direction) {
      SetCurrent(fMarkedEnd);
   } else {
      SetCurrent(fMarkedStart);
   }

   pos.fY = ToObjYCoord(fVisible.fY);
   if (fCurrent.fY < pos.fY ||
       ToScrYCoord(fCurrent.fY) >= (Int_t)fCanvas->GetHeight()) {
      pos.fY = fMarkedStart.fY;
   }
   pos.fX = ToObjXCoord(fVisible.fX, pos.fY);
   if (fCurrent.fX < pos.fX ||
       ToScrXCoord(fCurrent.fX, pos.fY) >= (Int_t)fCanvas->GetWidth()) {
      pos.fX = fMarkedStart.fX;
   }

   SetVsbPosition((ToScrYCoord(pos.fY)+fVisible.fY)/fScrollVal.fY);
   SetHsbPosition((ToScrXCoord(pos.fX, pos.fY)+fVisible.fX)/fScrollVal.fX);

   UpdateRegion(0, (Int_t)ToScrYCoord(fMarkedStart.fY), fCanvas->GetWidth(),
                UInt_t(ToScrYCoord(fMarkedEnd.fY+1)-ToScrYCoord(fMarkedEnd.fY)));

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Goto the specified line.

Bool_t TGTextEdit::Goto(Long_t line, Long_t column)
{
   if (line < 0)
      line = 0;
   if (line >= fText->RowCount())
      line = fText->RowCount() - 1;
   if (column < 0)
      column = 0;
   if (column > fText->GetLineLength(line))
      column = fText->GetLineLength(line);

   TGLongPosition gotopos, pos;
   gotopos.fY = line;
   gotopos.fX = column;
   SetCurrent(gotopos);

   pos.fY = ToObjYCoord(fVisible.fY);
   if (fCurrent.fY < pos.fY ||
       ToScrYCoord(fCurrent.fY) >= (Int_t)fCanvas->GetHeight())
      pos.fY = gotopos.fY;

   SetVsbPosition((ToScrYCoord(pos.fY)+fVisible.fY)/fScrollVal.fY);
   SetHsbPosition(0);

   UnMark();
   fIsMarked = kTRUE;
   fMarkedStart.fY = fMarkedEnd.fY = line;
   fMarkedStart.fX = 0;
   fMarkedEnd.fX = fCanvas->GetWidth();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the mode how characters are entered.

void TGTextEdit::SetInsertMode(EInsertMode mode)
{
   if (fInsertMode == mode) return;

   fInsertMode = mode;
}

////////////////////////////////////////////////////////////////////////////////
/// If cursor if on, turn it off.

void TGTextEdit::CursorOff()
{
   if (fCursorState == 1) {
      DrawCursor(2);
   }
   fCursorState = 2;
}

////////////////////////////////////////////////////////////////////////////////
/// Turn cursor on.

void TGTextEdit::CursorOn()
{
   DrawCursor(1);
   fCursorState = 1;

   if (fCurBlink) {
      fCurBlink->Reset();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Make the specified position the current position.

void TGTextEdit::SetCurrent(TGLongPosition new_coord)
{
   CursorOff();

   fCurrent.fY = new_coord.fY;
   fCurrent.fX = new_coord.fX;

   CursorOn();

   SendMessage(fMsgWindow, MK_MSG(kC_TEXTVIEW, kTXT_DATACHANGE), fWidgetId, 0);
   DataChanged();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw cursor. If mode = 1 draw cursor, if mode = 2 erase cursor.

void TGTextEdit::DrawCursor(Int_t mode)
{
   char count = -1;
   char cursor = ' ';
   if (fCurrent.fY >= fText->RowCount() || fCurrent.fX > fText->GetLineLength(fCurrent.fY) || fReadOnly) {
      return;
   }

   if (fCurrent.fY >= ToObjYCoord(fVisible.fY) &&
       fCurrent.fY <= ToObjYCoord(fVisible.fY+fCanvas->GetHeight()) &&
       fCurrent.fX >= ToObjXCoord(fVisible.fX, fCurrent.fY) &&
       fCurrent.fX <= ToObjXCoord(fVisible.fX+fCanvas->GetWidth(),fCurrent.fY)) {
      if (fCurrent.fY < fText->RowCount()) {
         count = fText->GetChar(fCurrent);
      }
      if (count == -1 || count == '\t') {
         cursor = ' ';
      } else {
         cursor = count;
      }

      if (mode == 2) {
         if (fIsMarked && count != -1) {
            if ((fCurrent.fY > fMarkedStart.fY && fCurrent.fY < fMarkedEnd.fY) ||
                (fCurrent.fY == fMarkedStart.fY && fCurrent.fX >= fMarkedStart.fX &&
                 fCurrent.fY < fMarkedEnd.fY) ||
                (fCurrent.fY == fMarkedEnd.fY && fCurrent.fX < fMarkedEnd.fX &&
                 fCurrent.fY > fMarkedStart.fY) ||
                (fCurrent.fY == fMarkedStart.fY && fCurrent.fY == fMarkedEnd.fY &&
                 fCurrent.fX >= fMarkedStart.fX && fCurrent.fX < fMarkedEnd.fX &&
                 fMarkedStart.fX != fMarkedEnd.fX)) {
               // back ground fillrectangle
               gVirtualX->FillRectangle(fCanvas->GetId(), fSelbackGC(),
                                     Int_t(ToScrXCoord(fCurrent.fX, fCurrent.fY)),
                                     Int_t(ToScrYCoord(fCurrent.fY)),
                                     UInt_t(ToScrXCoord(fCurrent.fX+1, fCurrent.fY) -
                                     ToScrXCoord(fCurrent.fX, fCurrent.fY)),
                                     UInt_t(ToScrYCoord(fCurrent.fY+1)-ToScrYCoord(fCurrent.fY)));
               if (count != -1)
                  gVirtualX->DrawString(fCanvas->GetId(), fSelGC(), (Int_t)ToScrXCoord(fCurrent.fX,fCurrent.fY),
                       Int_t(ToScrYCoord(fCurrent.fY+1) - fMaxDescent), &cursor, 1);
            } else {
               gVirtualX->ClearArea(fCanvas->GetId(),
                                    Int_t(ToScrXCoord(fCurrent.fX, fCurrent.fY)),
                                    Int_t(ToScrYCoord(fCurrent.fY)),
                                    UInt_t(ToScrXCoord(fCurrent.fX+1, fCurrent.fY) -
                                    ToScrXCoord(fCurrent.fX, fCurrent.fY)),
                                    UInt_t(ToScrYCoord(fCurrent.fY+1)-ToScrYCoord(fCurrent.fY)));
               if (count != -1)
                  gVirtualX->DrawString(fCanvas->GetId(), fNormGC(), (Int_t)ToScrXCoord(fCurrent.fX,fCurrent.fY),
                       Int_t(ToScrYCoord(fCurrent.fY+1) - fMaxDescent), &cursor, 1);
            }
         } else {
            gVirtualX->ClearArea(fCanvas->GetId(),
                                 Int_t(ToScrXCoord(fCurrent.fX, fCurrent.fY)),
                                 Int_t(ToScrYCoord(fCurrent.fY)),
                                 UInt_t(ToScrXCoord(fCurrent.fX+1, fCurrent.fY) -
                                 ToScrXCoord(fCurrent.fX, fCurrent.fY)),
                                 UInt_t(ToScrYCoord(fCurrent.fY+1)-ToScrYCoord(fCurrent.fY)));
            gVirtualX->DrawString(fCanvas->GetId(), fNormGC(), (Int_t)ToScrXCoord(fCurrent.fX,fCurrent.fY),
                       Int_t(ToScrYCoord(fCurrent.fY+1) - fMaxDescent), &cursor, 1);
         }
      } else {
         if (mode == 1) {
            gVirtualX->FillRectangle(fCanvas->GetId(), fCursor1GC,
                                     Int_t(ToScrXCoord(fCurrent.fX, fCurrent.fY)),
                                     Int_t(ToScrYCoord(fCurrent.fY)),
                                     2,
                                     UInt_t(ToScrYCoord(fCurrent.fY+1)-ToScrYCoord(fCurrent.fY)));
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Adjust current position.

void TGTextEdit::AdjustPos()
{
   TGLongPosition pos;
   pos.fY = fCurrent.fY;
   pos.fX = fCurrent.fX;

   if (pos.fY < ToObjYCoord(fVisible.fY)) {
      pos.fY = ToObjYCoord(fVisible.fY);
   } else if (ToScrYCoord(pos.fY+1) >= (Int_t) fCanvas->GetHeight()) {
      pos.fY = ToObjYCoord(fVisible.fY + fCanvas->GetHeight())-1;
   }
   if (pos.fX < ToObjXCoord(fVisible.fX, pos.fY)) {
      pos.fX = ToObjXCoord(fVisible.fX, pos.fY);
   } else if (ToScrXCoord(pos.fX, pos.fY) >= (Int_t) fCanvas->GetWidth()) {
      pos.fX = ToObjXCoord(fVisible.fX + fCanvas->GetWidth(), pos.fY)-1;
   }
   if (pos.fY != fCurrent.fY || pos.fX != fCurrent.fX) {
      SetCurrent(pos);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle timer cursor blink timer.

Bool_t TGTextEdit::HandleTimer(TTimer *t)
{
   if (t != fCurBlink) {
      TGTextView::HandleTimer(t);
      return kTRUE;
   }

   if (fCursorState == 1) {
      fCursorState = 2;
   } else {
      fCursorState = 1;
   }

   DrawCursor(fCursorState);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle selection notify event.

Bool_t TGTextEdit::HandleSelection(Event_t *event)
{
   TString data;
   Int_t   nchar;

   gVirtualX->GetPasteBuffer((Window_t)event->fUser[0], (Atom_t)event->fUser[3],
                             data, nchar, kFALSE);

   if (!nchar) return kTRUE;

   delete fClipText;

   fClipText = new TGText;
   fClipText->LoadBuffer(data.Data());

   TGLongPosition start_src, end_src, pos;

   pos.fX = pos.fY = 0;
   start_src.fY = start_src.fX = 0;
   end_src.fY = fClipText->RowCount()-1;
   end_src.fX = fClipText->GetLineLength(end_src.fY)-1;

   if (end_src.fX < 0) {
      end_src.fX = 0;
   }

   // undo command
   TInsTextCom *icom = new TInsTextCom(this);
   icom->fChar = fText->GetChar(fCurrent);
   fText->InsText(fCurrent, fClipText, start_src, end_src);

   fIsMarked = kFALSE;

   fExposedRegion.fX = 0;
   fExposedRegion.fY = ToScrYCoord(fCurrent.fY);

   pos.fY = fCurrent.fY + fClipText->RowCount()-1;
   pos.fX = fClipText->GetLineLength(fClipText->RowCount()-1);

   if (start_src.fY == end_src.fY) {
      pos.fX = pos.fX + fCurrent.fX;
   }

   icom->SetEndPos(pos);

   // calculate exposed region
   fExposedRegion.fW = fCanvas->GetWidth();
   fExposedRegion.fH = fCanvas->GetHeight() - fExposedRegion.fY;

   SetCurrent(pos);

   if (ToScrYCoord(pos.fY) >= (Int_t)fCanvas->GetHeight()) {
      pos.fY = ToScrYCoord(pos.fY) + fVisible.fY - fCanvas->GetHeight()/2;
      fExposedRegion.fX = fExposedRegion.fY = 0;
      fExposedRegion.fH = fCanvas->GetHeight();
   } else {
      pos.fY = fVisible.fY;
   }
   if (ToScrXCoord(pos.fX, fCurrent.fY) >= (Int_t) fCanvas->GetWidth()) {
      pos.fX = ToScrXCoord(pos.fX, fCurrent.fY) + fVisible.fX + fCanvas->GetWidth()/2;
   } else if (ToScrXCoord(pos.fX, fCurrent.fY < 0) && pos.fX != 0) {
      if (fVisible.fX - (Int_t)fCanvas->GetWidth()/2 > 0) {
         pos.fX = fVisible.fX - fCanvas->GetWidth()/2;
      } else {
         pos.fX = 0;
      }
   } else {
      pos.fX = fVisible.fX;
   }

   SetSBRange(kHorizontal);
   SetSBRange(kVertical);
   SetVsbPosition(pos.fY/fScrollVal.fY);
   SetHsbPosition(pos.fX/fScrollVal.fX);

   fClient->NeedRedraw(this);

   return kTRUE;
}

static Bool_t gDbl_clk = kFALSE;
static Bool_t gTrpl_clk = kFALSE;

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse button event in text edit widget.

Bool_t TGTextEdit::HandleButton(Event_t *event)
{
   if (event->fWindow != fCanvas->GetId()) {
      return kFALSE;
   }

   TGLongPosition pos;

   TGTextView::HandleButton(event);

   if (event->fType == kButtonPress) {
      SetFocus();
      //Update();

      if (event->fCode == kButton1 || event->fCode == kButton2) {
         pos.fY = ToObjYCoord(fVisible.fY + event->fY);
         if (pos.fY >= fText->RowCount()) {
            pos.fY = fText->RowCount()-1;
         }
         pos.fX = ToObjXCoord(fVisible.fX+event->fX, pos.fY);
         if (pos.fX >= fText->GetLineLength(pos.fY)) {
            pos.fX = fText->GetLineLength(pos.fY);
         }
         while (fText->GetChar(pos) == 16) {
            pos.fX++;
         }

         SetCurrent(pos);

         TGTextLine *line = fText->GetCurrentLine();
         char *word = line->GetWord(pos.fX);
         Clicked((const char*)word);   // emit signal
         delete [] word;
      }
      if (event->fCode == kButton2) {
         if (gVirtualX->GetPrimarySelectionOwner() != kNone) {
            gVirtualX->ConvertPrimarySelection(fId, fClipboard, event->fTime);
            Update();
            return kTRUE;
         }
      }
      if (event->fCode == kButton3) {
         // do not handle during guibuilding
         if (fClient->IsEditable() || !fEnableMenu) {
            return kTRUE;
         }
         SetMenuState();
         fMenu->PlaceMenu(event->fXRoot, event->fYRoot, kTRUE, kTRUE);
      }
      gDbl_clk = kFALSE;
      gTrpl_clk = kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle double click event.

Bool_t TGTextEdit::HandleDoubleClick(Event_t *event)
{
   if (event->fWindow != fCanvas->GetId()) {
      return kFALSE;
   }

   if (event->fCode != kButton1) {
      return kFALSE;
   }
   if (!fText->GetCurrentLine()->GetText()) {// empty line
      return kFALSE;
   }

   SetFocus();
   TGLongPosition pos;
   pos.fY = ToObjYCoord(fVisible.fY + event->fY);

   if (gDbl_clk && (event->fTime - fgLastClick < 350)) { // triple click
      fgLastClick  = event->fTime;
      gDbl_clk = kFALSE;
      gTrpl_clk = kTRUE;
      fMarkedStart.fY = fMarkedEnd.fY = pos.fY;
      fIsMarked = kTRUE;
      fMarkedStart.fX = 0;
      fMarkedEnd.fX = strlen(fText->GetCurrentLine()->GetText());
      Marked(kTRUE);
      UpdateRegion(0, (Int_t)ToScrYCoord(fMarkedStart.fY), fCanvas->GetWidth(),
                 UInt_t(ToScrYCoord(fMarkedEnd.fY + 1) - ToScrYCoord(fMarkedStart.fY)));
      return kTRUE;
   }

   if (gTrpl_clk && (event->fTime - fgLastClick < 350)) { // 4 click
      fgLastClick  = event->fTime;
      gTrpl_clk = kFALSE;
      fIsMarked = kTRUE;
      fMarkedStart.fY = 0;
      fMarkedStart.fX = 0;
      fMarkedEnd.fY = fText->RowCount()-1;
      fMarkedEnd.fX = fText->GetLineLength(fMarkedEnd.fY);
      if (fMarkedEnd.fX < 0) {
         fMarkedEnd.fX = 0;
      }
      UpdateRegion(0, 0, fCanvas->GetWidth(), fCanvas->GetHeight());
      return kTRUE;
   }

   gDbl_clk = kTRUE;
   gTrpl_clk = kFALSE;

   if (pos.fY >= fText->RowCount()) {
      pos.fY = fText->RowCount() - 1;
   }
   pos.fX = ToObjXCoord(fVisible.fX + event->fX, pos.fY);

   if (pos.fX >= fText->GetLineLength(pos.fY)) {
      pos.fX = fText->GetLineLength(pos.fY);
   }
   while (fText->GetChar(pos) == 16) {
      pos.fX++;
   }

   SetCurrent(pos);

   fMarkedStart.fY = fMarkedEnd.fY = pos.fY;
   char *line = fText->GetCurrentLine()->GetText();
   UInt_t len = (UInt_t)fText->GetCurrentLine()->GetLineLength();
   Int_t start = pos.fX;
   Int_t end = pos.fX;
   Int_t i = pos.fX;

   if (line[i] == ' ' || line[i] == '\t') {
      while (start >= 0) {
         if (line[start] == ' ' || line[start] == '\t') --start;
         else break;
      }
      ++start;
      while (end < (Int_t)len) {
         if (line[end] == ' ' || line[end] == '\t') ++end;
         else break;
      }
   } else if (isalnum(line[i])) {
      while (start >= 0) {
         if (isalnum(line[start])) --start;
         else break;
      }
      ++start;
      while (end < (Int_t)len) {
         if (isalnum(line[end])) ++end;
         else break;
      }
   } else {
      while (start >= 0) {
         if (isalnum(line[start]) || line[start] == ' ' || line[start] == '\t') {
            break;
         } else {
            --start;
         }
      }
      ++start;
      while (end < (Int_t)len) {
         if (isalnum(line[end]) || line[end] == ' ' || line[end] == '\t') {
            break;
         } else {
            ++end;
         }
      }
   }

   fMarkedStart.fX = start;
   fIsMarked = kTRUE;
   fMarkedEnd.fX = end;
   Marked(kTRUE);

   len = end - start; //length
   char *word = new char[len + 1];
   word[len] = '\0';
   strncpy(word, line+start, (UInt_t)len);
   DoubleClicked((const char *)word);  // emit signal

   delete [] word;
//   delete [] line;

   UpdateRegion(0, (Int_t)ToScrYCoord(fMarkedStart.fY), fCanvas->GetWidth(),
                UInt_t(ToScrYCoord(fMarkedEnd.fY + 1) - ToScrYCoord(fMarkedStart.fY)));

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion event in text edit widget.

Bool_t TGTextEdit::HandleMotion(Event_t *event)
{
   TGLongPosition pos;
   if (event->fWindow != fCanvas->GetId()) {
      return kTRUE;
   }

   if (fScrolling == -1) {
      pos.fY = ToObjYCoord(fVisible.fY+event->fY);
      if (pos.fY >= fText->RowCount()) {
         pos.fY = fText->RowCount()-1;
      }
      pos.fX = ToObjXCoord(fVisible.fX+event->fX, pos.fY);
      if (pos.fX > fText->GetLineLength(pos.fY)) {
         pos.fX = fText->GetLineLength(pos.fY);
      }
      if (fText->GetChar(pos) == 16) {
         if (pos.fX < fCurrent.fX) {
            pos.fX = fCurrent.fX;
         }
         if (pos.fX > fCurrent.fX) {
            do {
               pos.fX++;
            } while (fText->GetChar(pos) == 16);
         }
      }
      event->fY = (Int_t)ToScrYCoord(pos.fY);
      event->fX = (Int_t)ToScrXCoord(pos.fX, pos.fY);
      if (pos.fY != fCurrent.fY || pos.fX != fCurrent.fX) {
         TGTextView::HandleMotion(event);
         SetCurrent(pos);
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// The key press event handler converts a key press to some line editor
/// action.

Bool_t TGTextEdit::HandleKey(Event_t *event)
{
   Bool_t mark_ok = kFALSE;
   char   input[10];
   Int_t  n;
   UInt_t keysym;

   if (event->fType == kGKeyPress) {
      gVirtualX->LookupString(event, input, sizeof(input), keysym);
      n = strlen(input);

      AdjustPos();

      switch ((EKeySym)keysym) {   // ignore these keys
         case kKey_Shift:
         case kKey_Control:
         case kKey_Meta:
         case kKey_Alt:
         case kKey_CapsLock:
         case kKey_NumLock:
         case kKey_ScrollLock:
            return kTRUE;
         default:
            break;
      }
      if (event->fState & kKeyControlMask) {   // Cntrl key modifier pressed
         switch((EKeySym)keysym & ~0x20) {   // treat upper and lower the same
            case kKey_A:
               SelectAll();
               return kTRUE;
            case kKey_B:
               mark_ok = kTRUE;
               PrevChar();
               break;
            case kKey_C:
               Copy();
               return kTRUE;
            case kKey_D:
               if (fIsMarked) {
                  Cut();
               } else {
                  Long_t len = fText->GetLineLength(fCurrent.fY);
                  if (fCurrent.fY == fText->RowCount()-1 && fCurrent.fX == len) {
                     gVirtualX->Bell(0);
                     return kTRUE;
                  }
                  NextChar();
                  new TDelCharCom(this);
               }
               break;
            case kKey_E:
               mark_ok = kTRUE;
               End();
               break;
            case kKey_H:
               if (fCurrent.fX || fCurrent.fY) new TDelCharCom(this);
               else gVirtualX->Bell(0);
               break;
            case kKey_K:
               End();
               fIsMarked = kTRUE;
               Mark(fCurrent.fX, fCurrent.fY);
               Cut();
               break;
            case kKey_U:
               Home();
               UnMark();
               fMarkedStart.fY = fMarkedEnd.fY = fCurrent.fY;
               fMarkedStart.fX = fMarkedEnd.fX = fCurrent.fX;
               End();
               fIsMarked = kTRUE;
               Mark(fCurrent.fX, fCurrent.fY);
               Cut();
               break;
            case kKey_V:
            case kKey_Y:
               Paste();
               return kTRUE;
            case kKey_X:
               Cut();
               return kTRUE;
            case kKey_Z:
               fHistory->Notify();  // undo action
               return kTRUE;
            case kKey_F:
               Search(kFALSE);
               return kTRUE;
            case kKey_L:
            {
               Long_t ret = fCurrent.fY+1;
               new TGGotoDialog(fClient->GetDefaultRoot(), this, 400, 150, &ret);
               if (ret > -1) {
                  ret--;   // user specifies lines starting at 1
                  Goto(ret);
               }
               return kTRUE;
            }
            case kKey_Home:
               {
                  TGLongPosition pos;
                  pos.fY = 0;
                  pos.fX = 0;
                  SetHsbPosition(0);
                  SetVsbPosition(0);
                  SetCurrent(pos);
               }
               break;
            case kKey_End:
               {
                  TGLongPosition pos;
                  pos.fY = fText->RowCount()-1;
                  pos.fX = fText->GetLineLength(pos.fY);
                  if (fVsb && fVsb->IsMapped())
                     SetVsbPosition((ToScrYCoord(pos.fY)+fVisible.fY)/fScrollVal.fY);
                  SetCurrent(pos);
               }
               break;
            default:
               return kTRUE;
         }
      }
      if (n && keysym >= 32 && keysym < 127 &&     // printable keys
          !(event->fState & kKeyControlMask) &&
          (EKeySym)keysym != kKey_Delete &&
          (EKeySym)keysym != kKey_Backspace) {

         if (fIsMarked) {
            Cut();
         }
         new TInsCharCom(this, input[0]);

      } else {

         switch ((EKeySym)keysym) {
            case kKey_F3:
               // typically FindAgain action
               SendMessage(fMsgWindow, MK_MSG(kC_TEXTVIEW, kTXT_F3), fWidgetId,
                           kTRUE);
               SetMenuState();
               if (fMenu->IsEntryEnabled(kM_SEARCH_FINDAGAIN)) {
                  SendMessage(this, MK_MSG(kC_COMMAND, kCM_MENU),
                              kM_SEARCH_FINDAGAIN, 0);
                  FindAgain();
               }
               break;
            case kKey_Delete:
               if (fIsMarked) {
                  Cut();
               } else {
                  Long_t len = fText->GetLineLength(fCurrent.fY);
                  if (fCurrent.fY == fText->RowCount()-1 && fCurrent.fX == len) {
                     gVirtualX->Bell(0);
                     return kTRUE;
                  }
                  NextChar();
                  new TDelCharCom(this);
               }
               break;
            case kKey_Return:
            case kKey_Enter:
               new TBreakLineCom(this);
               break;
            case kKey_Tab:
               new TInsCharCom(this, '\t');
               break;
            case kKey_Backspace:
               if (fIsMarked) {
                  Cut();
               } else {
                  if (fCurrent.fX || fCurrent.fY) {
                     new TDelCharCom(this);
                  } else {
                     gVirtualX->Bell(0);
                  }
               }
               break;
            case kKey_Left:
               mark_ok = kTRUE;
               PrevChar();
               break;
            case kKey_Right:
               mark_ok = kTRUE;
               NextChar();
               break;
            case kKey_Up:
               mark_ok = kTRUE;
               LineUp();
               break;
            case kKey_Down:
               mark_ok = kTRUE;
               LineDown();
               break;
            case kKey_PageUp:
               mark_ok = kTRUE;
               ScreenUp();
               break;
            case kKey_PageDown:
               mark_ok = kTRUE;
               ScreenDown();
               break;
            case kKey_Home:
               mark_ok = kTRUE;
               Home();
               break;
            case kKey_End:
               mark_ok = kTRUE;
               End();
               break;
            case kKey_Insert:           // switch on/off insert mode
               SetInsertMode(GetInsertMode() == kInsert ? kReplace : kInsert);
               break;
            default:
               break;
         }
      }
      if ((event->fState & kKeyShiftMask) && mark_ok) {
         fIsMarked = kTRUE;
         Mark(fCurrent.fX, fCurrent.fY);
         Copy();
         SendMessage(fMsgWindow, MK_MSG(kC_TEXTVIEW, kTXT_ISMARKED), fWidgetId,
                     kTRUE);
         Marked(kTRUE);
      } else {
         UnMark();
         SendMessage(fMsgWindow, MK_MSG(kC_TEXTVIEW, kTXT_ISMARKED),
                        fWidgetId, kFALSE);
         fMarkedStart.fY = fMarkedEnd.fY = fCurrent.fY;
         fMarkedStart.fX = fMarkedEnd.fX = fCurrent.fX;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse crossing event.

Bool_t TGTextEdit::HandleCrossing(Event_t *event)
{
   if (event->fWindow != fCanvas->GetId()) {
      return kTRUE;
   }
   if (gVirtualX->GetInputFocus() != fCanvas->GetId()) {
      if (event->fType == kEnterNotify) {
         if (!fCurBlink) {
            fCurBlink = new TViewTimer(this, 500);
         }
         fCurBlink->Reset();
         gSystem->AddTimer(fCurBlink);
      } else {
         if (fCurBlink) fCurBlink->Remove();
         if (!fEnableCursorWithoutFocus && (fCursorState == 1)) {
            DrawCursor(2);
            fCursorState = 2;
         } else if (fCursorState == 2) {
            DrawCursor(1);
            fCursorState = 1;
         }
      }
   }

   TGTextView::HandleCrossing(event);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle focus change event in text edit widget.

Bool_t TGTextEdit::HandleFocusChange(Event_t *event)
{
   if (event->fWindow != fCanvas->GetId()) {
      return kTRUE;
   }

   // check this when porting to Win32
   if ((event->fCode == kNotifyNormal) && (event->fState != kNotifyPointer)) {
      if (event->fType == kFocusIn) {
         if (!fCurBlink) {
            fCurBlink = new TViewTimer(this, 500);
         }
         fCurBlink->Reset();
         gSystem->AddTimer(fCurBlink);
      } else {
         if (fCurBlink) fCurBlink->Remove();
         if (fCursorState == 2) {
            DrawCursor(1);
            fCursorState = 1;
         }
      }
      fClient->NeedRedraw(this);
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Invokes search dialog.

void TGTextEdit::Search(Bool_t close)
{
   static TGSearchType *srch = 0;
   Int_t ret = 0;

   if (!srch) srch = new TGSearchType;
   srch->fClose = close;

   if (!close) {
      if (!TGSearchDialog::SearchDialog()) {
         TGSearchDialog::SearchDialog() = new TGSearchDialog(fClient->GetDefaultRoot(),
                                                        fCanvas, 400, 150, srch, &ret);
      }
      TGSearchDialog::SearchDialog()->Connect("TextEntered(char *)", "TGTextEdit",
                                          this, "Search(char *,Bool_t,Bool_t)");
      TGSearchDialog::SearchDialog()->MapRaised();
   } else {
      new TGSearchDialog(fClient->GetDefaultRoot(), fCanvas, 400, 150, srch, &ret);
      if (ret) {
         Search(srch->fBuffer);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Process context menu messages.

Bool_t TGTextEdit::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
   TString msg2;
   TGTextView::ProcessMessage(msg, parm1, parm2);

   switch(GET_MSG(msg)) {
      case kC_COMMAND:
         switch(GET_SUBMSG(msg)) {
            case kCM_MENU:
               switch (parm1) {
                  case kM_FILE_NEW:
                  case kM_FILE_CLOSE:
                  case kM_FILE_OPEN:
                     if (!IsSaved()) {
                        Int_t retval;
                        Bool_t untitled = !strlen(fText->GetFileName()) ? kTRUE : kFALSE;

                        msg2.Form("Save \"%s\"?",
                                  untitled ? "Untitled" : fText->GetFileName());
                        new TGMsgBox(fClient->GetDefaultRoot(), this, "Editor",
                                     msg2.Data(), kMBIconExclamation,
                                     kMBYes | kMBNo | kMBCancel, &retval);

                        if (retval == kMBCancel)
                           return kTRUE;
                        if (retval == kMBYes)
                           if (!SaveFile(0))
                              return kTRUE;
                     }
                     Clear();
                     if (parm1 == kM_FILE_CLOSE) {
                        SendMessage(fMsgWindow, MK_MSG(kC_TEXTVIEW, kTXT_CLOSE),
                                    fWidgetId, 0);
                        Closed();
                     }
                     if (parm1 == kM_FILE_OPEN) {
                        TGFileInfo fi;
                        fi.fFileTypes = gFiletypes;
                        new TGFileDialog(fClient->GetDefaultRoot(), this, kFDOpen, &fi);
                        if (fi.fFilename && strlen(fi.fFilename)) {
                           LoadFile(fi.fFilename);
                           SendMessage(fMsgWindow, MK_MSG(kC_TEXTVIEW, kTXT_OPEN),
                                       fWidgetId, 0);
                           Opened();
                        }
                     }
                     break;
                  case kM_FILE_SAVE:
                     if (SaveFile(0)) {
                        SendMessage(fMsgWindow, MK_MSG(kC_TEXTVIEW, kTXT_SAVE),
                                    fWidgetId, 0);
                        Saved();
                     }
                     break;
                  case kM_FILE_SAVEAS:
                     if (SaveFile(0, kTRUE)) {
                        SendMessage(fMsgWindow, MK_MSG(kC_TEXTVIEW, kTXT_SAVE),
                                    fWidgetId, 0);
                        SavedAs();
                     }
                     break;
                  case kM_FILE_PRINT:
                     {
                        Int_t ret = 0;
                        if (!gPrinter) {
                           gPrinter = StrDup("892_2_cor"); // use gEnv
                           gPrintCommand = StrDup("xprint");
                        }
                        new TGPrintDialog(fClient->GetDefaultRoot(), this, 400, 150,
                                          &gPrinter, &gPrintCommand, &ret);
                        if (ret)
                           Print();
                     }
                     break;
                  case kM_EDIT_CUT:
                     Cut();
                     break;
                  case kM_EDIT_COPY:
                     Copy();
                     break;
                  case kM_EDIT_PASTE:
                     Paste();
                     break;
                  case kM_EDIT_SELECTALL:
                     SelectAll();
                     break;
                  case kM_SEARCH_FIND:
                     {
                        Search(kFALSE);
                     }
                     break;
                  case kM_SEARCH_FINDAGAIN:
                     if (!fSearch) {
                        SendMessage(this, MK_MSG(kC_COMMAND, kCM_MENU),
                                    kM_SEARCH_FIND, 0);
                        return kTRUE;
                     }
                     if (!Search(fSearch->fBuffer, fSearch->fDirection,
                                 fSearch->fCaseSensitive)) {
                        msg2.Form("Couldn't find \"%s\"", fSearch->fBuffer);
                        new TGMsgBox(fClient->GetDefaultRoot(), this, "Editor",
                                     msg2.Data(), kMBIconExclamation, kMBOk, 0);
                     }
                     break;
                  case kM_SEARCH_GOTO:
                     {
                        Long_t ret = fCurrent.fY+1;
                        new TGGotoDialog(fClient->GetDefaultRoot(), this, 400, 150, &ret);
                        if (ret > -1) {
                           ret--;   // user specifies lines starting at 1
                           Goto(ret);
                        }
                     }
                     break;
                  default:
                     printf("No action implemented for menu id %ld\n", parm1);
                     break;
               }
            default:
               break;
         }
         break;

      default:
         break;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Insert a character in the text edit widget.

void TGTextEdit::InsChar(char character)
{
   if (fReadOnly) return;

   char *charstring = 0;
   TGLongPosition pos;

   if (character == '\t') {
      pos.fX = fCurrent.fX;
      pos.fY = fCurrent.fY;
      fText->InsChar(pos, '\t');
      pos.fX++;
      while (pos.fX & 0x7) {
         pos.fX++;
      }
      fText->ReTab(pos.fY);
      UpdateRegion(0, (Int_t)ToScrYCoord(pos.fY), fCanvas->GetWidth(),
                   UInt_t(ToScrYCoord(pos.fY+1) - ToScrYCoord(pos.fY)));
      SetSBRange(kHorizontal);
      if (ToScrXCoord(pos.fX, pos.fY) >= (Int_t)fCanvas->GetWidth()) {
         if (pos.fX != fText->GetLineLength(fCurrent.fY)) {
            SetHsbPosition((fVisible.fX+fCanvas->GetWidth()/2)/fScrollVal.fX);
         } else {
            SetHsbPosition(fVisible.fX/fScrollVal.fX);
         }
      }
      SetCurrent(pos);
      return;
   } else {
      if (fInsertMode == kReplace) {
         fCurrent.fX++;
         new TDelCharCom(this);
      }
      fText->InsChar(fCurrent, character);
      pos.fX = fCurrent.fX + 1;
      pos.fY = fCurrent.fY;
      charstring = new char[2];
      charstring[1] = '\0';
      charstring[0] = character;
   }
   SetSBRange(kHorizontal);
   if (ToScrXCoord(pos.fX, pos.fY) >= (Int_t)fCanvas->GetWidth()) {
      if (pos.fX != fText->GetLineLength(fCurrent.fY)) {
         SetHsbPosition((fVisible.fX+fCanvas->GetWidth()/2)/fScrollVal.fX);
      } else {
         SetHsbPosition(fVisible.fX/fScrollVal.fX+strlen(charstring));
      }
      if (!fHsb)
         gVirtualX->DrawString(fCanvas->GetId(), fNormGC(),
                               (Int_t)ToScrXCoord(fCurrent.fX, fCurrent.fY),
                               Int_t(ToScrYCoord(fCurrent.fY+1) - fMaxDescent),
                               charstring, strlen(charstring));
   } else {
#ifdef R__HAS_COCOA
      //I would use const, but some members of TGTextLine are non-const.
      if (TGTextLine *currentLine = fText->GetCurrentLine()) {
         const ULong_t lineStart = ToObjXCoord(fVisible.fX, fCurrent.fY);
         if (lineStart < currentLine->GetLineLength()) {
            const char *textToRender = currentLine->GetText(lineStart, currentLine->GetLineLength() - lineStart);
            //The next two lines can throw and textToRender will leak, but ROOT does not care about such things. :(
            gVirtualX->ClearArea(fCanvas->GetId(), Int_t(ToScrXCoord(0, fCurrent.fY)),
                                 Int_t(ToScrYCoord(fCurrent.fY)), UInt_t(ToScrXCoord(currentLine->GetLineLength(), fCurrent.fY)),
                                 UInt_t(ToScrYCoord(fCurrent.fY+1)-ToScrYCoord(fCurrent.fY)));
            gVirtualX->DrawString(fCanvas->GetId(), fNormGC(), Int_t(ToScrXCoord(0, fCurrent.fY)),
                                  Int_t(ToScrYCoord(fCurrent.fY + 1) - fMaxDescent),
                                  textToRender, -1);
            delete [] textToRender;
         }
      }
#else
      gVirtualX->CopyArea(fCanvas->GetId(), fCanvas->GetId(), fNormGC(),
                          (Int_t)ToScrXCoord(fCurrent.fX, fCurrent.fY),
                          (Int_t)ToScrYCoord(fCurrent.fY), fCanvas->GetWidth(),
                          UInt_t(ToScrYCoord(fCurrent.fY+1)-ToScrYCoord(fCurrent.fY)),
                          (Int_t)ToScrXCoord(pos.fX, fCurrent.fY),
                          (Int_t)ToScrYCoord(fCurrent.fY));
      gVirtualX->ClearArea(fCanvas->GetId(),
                           Int_t(ToScrXCoord(fCurrent.fX, fCurrent.fY)),
                           Int_t(ToScrYCoord(fCurrent.fY)),
                           UInt_t(ToScrXCoord(fCurrent.fX+strlen(charstring), fCurrent.fY) -
                           ToScrXCoord(fCurrent.fX, fCurrent.fY)),
                           UInt_t(ToScrYCoord(fCurrent.fY+1)-ToScrYCoord(fCurrent.fY)));
      gVirtualX->DrawString(fCanvas->GetId(), fNormGC(),
                            Int_t(ToScrXCoord(fCurrent.fX, fCurrent.fY)),
                            Int_t(ToScrYCoord(fCurrent.fY+1) - fMaxDescent),
                            charstring, strlen(charstring));
      fCursorState = 2;  // the ClearArea effectively turned off the cursor
#endif
   }
   delete [] charstring;
   SetCurrent(pos);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a character from the text edit widget.

void TGTextEdit::DelChar()
{
   if (fReadOnly) {
      return;
   }

   char *buffer;
   TGLongPosition pos, pos2;
   Long_t len;

   pos.fY = fCurrent.fY;
   pos.fX = fCurrent.fX;
   UInt_t h = 0;

   if (fCurrent.fX > 0) {
      Int_t y = (Int_t)ToScrYCoord(pos.fY);
      h = UInt_t(ToScrYCoord(pos.fY+2) - y);
      if (!y) h = h << 1;

      pos.fX--;
      if (fText->GetChar(pos) == 16) {
         do {
            pos.fX++;
            fText->DelChar(pos);
            pos.fX -= 2;
         } while (fText->GetChar(pos) != '\t');

         pos.fX++;
         fText->DelChar(pos);
         pos.fX--;
         fText->ReTab(pos.fY);
         UpdateRegion(0, y, fCanvas->GetWidth(), h);
      } else {
         pos.fX = fCurrent.fX;
         fText->DelChar(pos);
         pos.fX = fCurrent.fX - 1;
      }
      if (ToScrXCoord(fCurrent.fX-1, fCurrent.fY) < 0) {
         SetHsbPosition((fVisible.fX-fCanvas->GetWidth()/2)/fScrollVal.fX);
      }
      SetSBRange(kHorizontal);
      UpdateRegion(0, y, fCanvas->GetWidth(), h);
   } else {
      if (fCurrent.fY > 0) {
         len = fText->GetLineLength(fCurrent.fY);
         if (len > 0) {
            buffer = fText->GetLine(fCurrent, len);
            pos.fY--;
            pos.fX = fText->GetLineLength(fCurrent.fY-1);
            fText->InsText(pos, buffer);
            pos.fY++;
            delete [] buffer;
         } else {
            pos.fX = fText->GetLineLength(fCurrent.fY-1);
         }

         pos2.fY = ToScrYCoord(fCurrent.fY+1);
         pos.fY = fCurrent.fY - 1;
         fText->DelLine(fCurrent.fY);
         len = fText->GetLineLength(fCurrent.fY-1);

         if (ToScrXCoord(pos.fX, fCurrent.fY-1) >= (Int_t)fCanvas->GetWidth()) {
            SetHsbPosition((ToScrXCoord(pos.fX, pos.fY)+fVisible.fX-fCanvas->GetWidth()/2)/fScrollVal.fX);
         }

#ifdef R__HAS_COCOA
         UpdateRegion(0, 0, fCanvas->GetWidth(), fCanvas->GetHeight());
#else
         h = UInt_t(fCanvas->GetHeight() - ToScrYCoord(fCurrent.fY));
         gVirtualX->CopyArea(fCanvas->GetId(), fCanvas->GetId(), fNormGC(), 0,
                             Int_t(pos2.fY), fWidth, h, 0, (Int_t)ToScrYCoord(fCurrent.fY));
         if (ToScrYCoord(pos.fY) < 0) {
            SetVsbPosition(fVisible.fY/fScrollVal.fY-1);
         }
         UpdateRegion(0, (Int_t)ToScrYCoord(pos.fY), fCanvas->GetWidth(), h);
#endif
         SetSBRange(kVertical);
         SetSBRange(kHorizontal);
      }
   }

   SetCurrent(pos);
}

////////////////////////////////////////////////////////////////////////////////
/// Break a line.

void TGTextEdit::BreakLine()
{
   if (fReadOnly) return;

   TGLongPosition pos;
   fText->BreakLine(fCurrent);
   if (ToScrYCoord(fCurrent.fY+2) <= (Int_t)fCanvas->GetHeight()) {
#ifdef R__HAS_COCOA
      UpdateRegion(0, (Int_t)ToScrYCoord(fCurrent.fY + 1), fCanvas->GetWidth(), fCanvas->GetHeight());
#else
      gVirtualX->CopyArea(fCanvas->GetId(), fCanvas->GetId(), fNormGC(), 0,
                          (Int_t)ToScrYCoord(fCurrent.fY+1), fCanvas->GetWidth(),
                          UInt_t(fCanvas->GetHeight()-(ToScrYCoord(fCurrent.fY+2)-
                          ToScrYCoord(fCurrent.fY))),
                          0, (Int_t)ToScrYCoord(fCurrent.fY+2));
      UpdateRegion(0, (Int_t)ToScrYCoord(fCurrent.fY), fCanvas->GetWidth(),
                  UInt_t(ToScrYCoord(fCurrent.fY+2) - ToScrYCoord(fCurrent.fY)));
#endif
      if (fVisible.fX != 0) {
         SetHsbPosition(0);
      }
      SetSBRange(kHorizontal);
      SetSBRange(kVertical);
   } else {
      SetSBRange(kHorizontal);
      SetSBRange(kVertical);
      SetVsbPosition(fVisible.fY/fScrollVal.fY + 1);
      UpdateRegion(0, (Int_t)ToScrYCoord(fCurrent.fY), fCanvas->GetWidth(),
                   UInt_t(ToScrYCoord(fCurrent.fY+1) - ToScrYCoord(fCurrent.fY)));
   }
   pos.fY = fCurrent.fY+1;
   pos.fX = 0;
   SetCurrent(pos);
}

////////////////////////////////////////////////////////////////////////////////
/// Scroll the canvas to new_top in the kVertical or kHorizontal direction.

void TGTextEdit::ScrollCanvas(Int_t new_top, Int_t direction)
{
   CursorOff();

   TGTextView::ScrollCanvas(new_top, direction);

   CursorOn();
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw the text edit widget.

void TGTextEdit::DrawRegion(Int_t x, Int_t y, UInt_t width, UInt_t height)
{
   CursorOff();

   TGTextView::DrawRegion(x, y, width, height);

   CursorOn();
}

////////////////////////////////////////////////////////////////////////////////
/// Go to the previous character.

void TGTextEdit::PrevChar()
{
   if (fCurrent.fY == 0 && fCurrent.fX == 0) {
      gVirtualX->Bell(0);
      return;
   }

   TGLongPosition pos;
   Long_t len;

   pos.fY = fCurrent.fY;
   pos.fX = fCurrent.fX;
   if (fCurrent.fX > 0) {
      pos.fX--;
      while (fText->GetChar(pos) == 16) {
         pos.fX--;
      }

      if (ToScrXCoord(pos.fX, pos.fY) < 0) {
         if (fVisible.fX-(Int_t)fCanvas->GetWidth()/2 >= 0) {
            SetHsbPosition((fVisible.fX-fCanvas->GetWidth()/2)/fScrollVal.fX);
         } else {
            SetHsbPosition(0);
         }
      }
   } else {
      if (fCurrent.fY > 0) {
         pos.fY = fCurrent.fY - 1;
         len = fText->GetLineLength(pos.fY);
         if (ToScrYCoord(fCurrent.fY) <= 0) {
            SetVsbPosition(fVisible.fY/fScrollVal.fY-1);
         }
         if (ToScrXCoord(len, pos.fY) >= (Int_t)fCanvas->GetWidth()) {
            SetHsbPosition((ToScrXCoord(len, pos.fY)+fVisible.fX -
                            fCanvas->GetWidth()/2)/fScrollVal.fX);
         }
         pos.fX = len;
      }
   }
   SetCurrent(pos);
}

////////////////////////////////////////////////////////////////////////////////
/// Go to next character.

void TGTextEdit::NextChar()
{
   Long_t len = fText->GetLineLength(fCurrent.fY);

   if (fCurrent.fY == fText->RowCount()-1 && fCurrent.fX == len) {
      gVirtualX->Bell(0);
      return;
   }

   TGLongPosition pos;
   pos.fY = fCurrent.fY;
   if (fCurrent.fX < len) {
      if (fText->GetChar(fCurrent) == '\t') {
         pos.fX = fCurrent.fX + 8 - (fCurrent.fX & 0x7);
      } else {
         pos.fX = fCurrent.fX + 1;
      }

      if (ToScrXCoord(pos.fX, pos.fY) >= (Int_t)fCanvas->GetWidth()) {
         SetHsbPosition(fVisible.fX/fScrollVal.fX+(fCanvas->GetWidth()/2)/fScrollVal.fX);
      }
   } else {
      if (fCurrent.fY < fText->RowCount()-1) {
         pos.fY = fCurrent.fY + 1;
         if (ToScrYCoord(pos.fY+1) >= (Int_t)fCanvas->GetHeight()) {
            SetVsbPosition(fVisible.fY/fScrollVal.fY+1);
         }
         SetHsbPosition(0);
         pos.fX = 0;
      }
   }
   SetCurrent(pos);
}

////////////////////////////////////////////////////////////////////////////////
/// Make current position first line in window by scrolling up.

void TGTextEdit::LineUp()
{
   TGLongPosition pos;
   Long_t len;
   if (fCurrent.fY > 0) {
      pos.fY = fCurrent.fY - 1;
      if (ToScrYCoord(fCurrent.fY) <= 0) {
         SetVsbPosition(fVisible.fY/fScrollVal.fY-1);
      }
      len = fText->GetLineLength(fCurrent.fY-1);
      if (fCurrent.fX > len) {
         if (ToScrXCoord(len, pos.fY) <= 0) {
            if (ToScrXCoord(len, pos.fY) < 0) {
               SetHsbPosition(ToScrXCoord(len, pos.fY)+
                            (fVisible.fX-fCanvas->GetWidth()/2)/fScrollVal.fX);
            } else {
               SetHsbPosition(0);
            }
         }
         pos.fX = len;
      } else {
         pos.fX = ToObjXCoord(ToScrXCoord(fCurrent.fX, fCurrent.fY)+fVisible.fX, pos.fY);
      }

      while (fText->GetChar(pos) == 16) {
         pos.fX++;
      }
      SetCurrent(pos);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Move one line down.

void TGTextEdit::LineDown()
{
   TGLongPosition pos;
   Long_t len;
   if (fCurrent.fY < fText->RowCount()-1) {
      len = fText->GetLineLength(fCurrent.fY+1);
      pos.fY = fCurrent.fY + 1;
      if (ToScrYCoord(pos.fY+1) > (Int_t)fCanvas->GetHeight()) {
         SetVsbPosition(fVisible.fY/fScrollVal.fY+1);
      }
      if (fCurrent.fX > len) {
         if (ToScrXCoord(len, pos.fY) <= 0) {
            if (ToScrXCoord(len, pos.fY) < 0) {
               SetHsbPosition((ToScrXCoord(len, pos.fY)+fVisible.fX-fCanvas->GetWidth()/2)/fScrollVal.fX);
            } else {
               SetHsbPosition(0);
            }
         }
         pos.fX = len;
      } else {
         pos.fX = ToObjXCoord(ToScrXCoord(fCurrent.fX, fCurrent.fY)+fVisible.fX, pos.fY);
      }

      while (fText->GetChar(pos) == 16) {
         pos.fX++;
      }
      SetCurrent(pos);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Move one screen up.

void TGTextEdit::ScreenUp()
{
   TGLongPosition pos;
   pos.fX = fCurrent.fX;
   pos.fY = fCurrent.fY - (ToObjYCoord(fCanvas->GetHeight())-ToObjYCoord(0))-1;
   if (fVisible.fY - (Int_t)fCanvas->GetHeight() >= 0) { // +1
      SetVsbPosition((fVisible.fY - fCanvas->GetHeight())/fScrollVal.fY);
   } else {
      pos.fY = 0;
      SetVsbPosition(0);
   }
   while (fText->GetChar(pos) == 16) {
      pos.fX++;
   }
   SetCurrent(pos);
}

////////////////////////////////////////////////////////////////////////////////
/// Move one screen down.

void TGTextEdit::ScreenDown()
{
   TGLongPosition pos;
   pos.fX = fCurrent.fX;
   pos.fY = fCurrent.fY + (ToObjYCoord(fCanvas->GetHeight()) - ToObjYCoord(0));
   Long_t count = fText->RowCount()-1;
   if ((Int_t)fCanvas->GetHeight() < ToScrYCoord(count)) {
      SetVsbPosition((fVisible.fY+fCanvas->GetHeight())/fScrollVal.fY);
   } else {
      pos.fY = count;
   }
   while (fText->GetChar(pos) == 16) {
      pos.fX++;
   }
   SetCurrent(pos);
}

////////////////////////////////////////////////////////////////////////////////
/// Move to beginning of line.

void TGTextEdit::Home()
{
   TGLongPosition pos;
   pos.fY = fCurrent.fY;
   pos.fX = 0;
   SetHsbPosition(0);
   SetCurrent(pos);
}

////////////////////////////////////////////////////////////////////////////////
/// Move to end of line.

void TGTextEdit::End()
{
   TGLongPosition pos;
   pos.fY = fCurrent.fY;
   pos.fX = fText->GetLineLength(pos.fY);
   if (ToScrXCoord(pos.fX, pos.fY) >= (Int_t)fCanvas->GetWidth()) {
      SetHsbPosition((ToScrXCoord(pos.fX, pos.fY) + fVisible.fX - fCanvas->GetWidth()/2)/fScrollVal.fX);
   }
   SetCurrent(pos);
}

////////////////////////////////////////////////////////////////////////////////
/// Return selection graphics context for text cursor.

const TGGC &TGTextEdit::GetCursor0GC()
{
   if (!fgCursor0GC) {
      fgCursor0GC = new TGGC(GetDefaultSelectedGC());
      fgCursor0GC->SetFunction(kGXxor);
   }
   return *fgCursor0GC;
}

////////////////////////////////////////////////////////////////////////////////
/// Return default graphics context for text cursor.

const TGGC &TGTextEdit::GetCursor1GC()
{
   if (!fgCursor1GC) {
      fgCursor1GC = new TGGC(GetDefaultGC());
      fgCursor1GC->SetFunction(kGXand);
   }
   return *fgCursor1GC;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a text edit widget as a C++ statement(s) on output stream out

void TGTextEdit::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';
   out << "   TGTextEdit *";
   out << GetName() << " = new TGTextEdit(" << fParent->GetName()
       << "," << GetWidth() << "," << GetHeight()
       << ");"<< std::endl;
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   if (IsReadOnly()) {
      out << "   " << GetName() << "->SetReadOnly(kTRUE);" << std::endl;
   }

   if (!IsMenuEnabled()) {
      out << "   " << GetName() << "->EnableMenu(kFALSE);" << std::endl;
   }

   if (fCanvas->GetBackground() != TGFrame::fgWhitePixel) {
      out << "   " << GetName() << "->ChangeBackground(" << fCanvas->GetBackground() << ");" << std::endl;
   }

   TGText *txt = GetText();
   Bool_t fromfile = strlen(txt->GetFileName()) ? kTRUE : kFALSE;
   TString fn;

   if (fromfile) {
      const char *filename = txt->GetFileName();
      fn = gSystem->UnixPathName(filename);
      gSystem->ExpandPathName(fn);
   } else {
      fn = TString::Format("Txt%s", GetName()+5);
      txt->Save(fn.Data());
   }
   out << "   " << GetName() << "->LoadFile(" << quote << fn.Data() << quote << ");" << std::endl;
}
