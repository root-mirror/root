// @(#)root/gui:$Id: ca18aae7f02ae8574b0c5414c78ede465425855d $
// Author: Bertrand Bellenot   20/06/06

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
//  TGTextEditor                                                        //
//                                                                      //
//  A simple text editor that uses the TGTextEdit widget.               //
//  It provides all functionalities of TGTextEdit as copy, paste, cut,  //
//  search, go to a given line number. In addition, it provides the     //
//  possibilities for compiling, executing or interrupting a running    //
//  macro.                                                              //
//                                                                      //
//  This class can be used in following ways:                           //
//  - with file name as argument:                                       //
//    new TGTextEditor("hsimple.C");                                    //
//  - with a TMacro* as argument:                                       //
//    TMacro *macro = new TMacro("hsimple.C");                          //
//    new TGTextEditor(macro);                                          //
//                                                                      //
//  Basic Features:                                                     //
//                                                                      //
//  New Document                                                        //
//                                                                      //
//  To create a new blank document, select File menu / New, or click    //
//  the New toolbar button. It will create a new instance of            //
//  TGTextEditor.                                                       //
//                                                                      //
//  Open/Save File                                                      //
//                                                                      //
//  To open a file, select File menu / Open or click on the Open        //
//  toolbar button. This will bring up the standard File Dialog for     //
//  opening files.                                                      //
//  If the current document has not been saved yet, you will be asked   //
//  either to save or abandon the changes.                              //
//  To save the file using the same name, select File menu / Save or    //
//  the toolbar Save button. To change the file name use File menu /    //
//  Save As... or corresponding SaveAs button on the toolbar.           //
//                                                                      //
//  Text Selection                                                      //
//                                                                      //
//  You can move the cursor by simply clicking on the desired location  //
//  with the left mouse button. To highlight some text, press the mouse //
//  and drag the mouse while holding the left button pressed.           //
//  To select a word, double-click on it;                               //
//  to select the text line - triple-click on it;                       //
//  to select all  do quadruple-click.                                  //
//                                                                      //
//  Cut, Copy, Paste                                                    //
//                                                                      //
//  After selecting some text, you can cut or copy it to the clipboard. //
//  A subsequent paste operation will insert the contents of the        //
//  clipboard at the current cursor location.                           //
//                                                                      //
//  Text Search                                                         //
//                                                                      //
//  The editor uses a standard Search dialog. You can specify a forward //
//  or backward search direction starting from the current cursor       //
//  location according to the selection made of a case sensitive mode   //
//  or not. The last search can be repeated by pressing F3.             //
//                                                                      //
//  Text Font                                                           //
//                                                                      //
//  You can change the text font by selecting Edit menu / Set Font.     //
//  The Font Dialog pops up and shows the Name, Style, and Size of any  //
//  available font. The selected font sample is shown in the preview    //
//  area.                                                               //
//                                                                      //
//  Executing Macros                                                    //
//                                                                      //
//  You can execute the currently loaded macro in the editor by         //
//  selecting Tools menu / Execute Macro; by clicking on the            //
//  corresponding toolbar button, or by using Ctrl+F5 accelerator keys. //
//  This is identical to the command ".x macro.C" in the root prompt    //
//  command line.                                                       //
//                                                                      //
//  Compiling Macros                                                    //
//                                                                      //
//  The currently loaded macro can be compiled with ACLiC if you select //
//  Tools menu / Compile Macro; by clicking on the corresponding        //
//  toolbar button, or by using Ctrl+F7 accelerator keys.               //
//  This is identical to the command ".L macro.C++" in the root prompt  //
//  command line.                                                       //
//                                                                      //
//  Interrupting a Running Macro                                        //
//                                                                      //
//  You can interrupt a running macro by selecting the Tools menu /     //
//  Interrupt; by clicking on the corresponding toolbar button, or by   //
//  using Shift+F5 accelerator keys.                                    //
//                                                                      //
//  Interface to CINT Interpreter                                       //
//                                                                      //
//  Any command entered in the Command combo box will be passed to      //
//  the CINT interpreter. This combo box will keep the commands history //
//  and will allow you to re-execute the same commands during an editor //
//  session.                                                            //
//                                                                      //
//  Keyboard Bindings                                                   //
//                                                                      //
//  The following table lists the keyboard shortcuts and accelerator    //
//  keys.                                                               //
//                                                                      //
//  Key:              Action:                                           //
//  ====              =======                                           //
//                                                                      //
//  Up                Move cursor up.                                   //
//  Shift+Up          Move cursor up and extend selection.              //
//  Down              Move cursor down.                                 //
//  Shift+Down        Move cursor down and extend selection.            //
//  Left              Move cursor left.                                 //
//  Shift+Left        Move cursor left and extend selection.            //
//  Right             Move cursor right.                                //
//  Shift+Right       Move cursor right and extend selection.           //
//  Home              Move cursor to begin of line.                     //
//  Shift+Home        Move cursor to begin of line and extend selection.//
//  Ctrl+Home         Move cursor to top of page.                       //
//  End               Move cursor to end of line.                       //
//  Shift+End         Move cursor to end of line and extend selection.  //
//  Ctrl+End          Move cursor to end of page.                       //
//  PgUp              Move cursor up one page.                          //
//  Shift+PgUp        Move cursor up one page and extend selection.     //
//  PgDn              Move cursor down one page.                        //
//  Shift+PgDn        Move cursor down one page and extend selection.   //
//  Delete            Delete character after cursor, or text selection. //
//  BackSpace         Delete character before cursor, or text selection.//
//  Ctrl+B            Move cursor left.                                 //
//  Ctrl+D            Delete character after cursor, or text selection. //
//  Ctrl+E            Move cursor to end of line.                       //
//  Ctrl+H            Delete character before cursor, or text selection.//
//  Ctrl+K            Delete characters from current position to the    //
//                    end of line.                                      //
//  Ctrl+U            Delete current line.                              //
//                                                                      //
//Begin_Html
/*
<img src="gif/TGTextEditor.gif">
*/
//End_Html
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TMacro.h"
#include "TInterpreter.h"
#include "TGMsgBox.h"
#include "TGFileDialog.h"
#include "TGFontDialog.h"
#include "TGTextEdit.h"
#include "TGMenu.h"
#include "TGStatusBar.h"
#include "KeySymbols.h"
#include "TGToolBar.h"
#include "TG3DLine.h"
#include "TGLabel.h"
#include "TGTextEntry.h"
#include "TGTextEditDialogs.h"
#include "TGTextEditor.h"
#include "TGComboBox.h"
#include "TObjString.h"
#include "TRootHelpDialog.h"
#include "HelpText.h"
#include "TVirtualX.h"
#ifdef WIN32
#include "TWin32SplashThread.h"
#endif

const char *ed_filetypes[] = {
   "ROOT Macros",  "*.C",
   "Source files", "*.cxx",
   "Text files",   "*.txt",
   "All files",    "*",
   0, 0
};

enum ETextEditorCommands {
   kM_FILE_NEW, kM_FILE_OPEN, kM_FILE_SAVE, kM_FILE_SAVEAS, kM_FILE_CLOSE,
   kM_FILE_PRINT, kM_FILE_EXIT, kM_EDIT_CUT, kM_EDIT_COPY, kM_EDIT_PASTE,
   kM_EDIT_DELETE, kM_EDIT_SELECTALL, kM_SEARCH_FIND, kM_SEARCH_FINDNEXT,
   kM_SEARCH_GOTO, kM_TOOLS_COMPILE, kM_TOOLS_EXECUTE, kM_TOOLS_INTERRUPT,
   kM_HELP_CONTENTS, kM_HELP_ABOUT, kM_EDIT_SELFONT
};

ToolBarData_t fTbData[] = {
  { "ed_new.png",       "New File",         kFALSE, kM_FILE_NEW,         0 },
  { "ed_open.png",      "Open File",        kFALSE, kM_FILE_OPEN,        0 },
  { "ed_save.png",      "Save File",        kFALSE, kM_FILE_SAVE,        0 },
  { "ed_saveas.png",    "Save File As...",  kFALSE, kM_FILE_SAVEAS,      0 },
  { "",                 0,                  0,      -1,                  0 },
  { "ed_print.png",     "Print",            kFALSE, kM_FILE_PRINT,       0 },
  { "",                 0,                  0,      -1,                  0 },
  { "ed_cut.png",       "Cut selection",    kFALSE, kM_EDIT_CUT,         0 },
  { "ed_copy.png",      "Copy selection",   kFALSE, kM_EDIT_COPY,        0 },
  { "ed_paste.png",     "Paste selection",  kFALSE, kM_EDIT_PASTE,       0 },
  { "ed_delete.png",    "Delete selection", kFALSE, kM_EDIT_DELETE,      0 },
  { "",                 0,                  0,      -1,                  0 },
  { "ed_find.png",      "Find...",          kFALSE, kM_SEARCH_FIND,      0 },
  { "ed_findnext.png",  "Find next",        kFALSE, kM_SEARCH_FINDNEXT,  0 },
  { "ed_goto.png",      "Goto...",          kFALSE, kM_SEARCH_GOTO,      0 },
  { "",                 0,                  0,      -1,                  0 },
  { "ed_compile.png",   "Compile Macro",    kFALSE, kM_TOOLS_COMPILE,    0 },
  { "ed_execute.png",   "Execute Macro",    kFALSE, kM_TOOLS_EXECUTE,    0 },
  { "ed_interrupt.png", "Interrupt",        kFALSE, kM_TOOLS_INTERRUPT,  0 },
  { "",                 0,                  0,      -1,                  0 },
  { "ed_help.png",      "Help Contents",    kFALSE, kM_HELP_CONTENTS,    0 },
  { "",                 0,                  0,      -1,                  0 },
  { "ed_quit.png",      "Close Editor",     kFALSE, kM_FILE_EXIT,        0 },
  {  0,                 0,                  0,      0,                   0 }
};

static char *gEPrinter      = 0;
static char *gEPrintCommand = 0;

ClassImp(TGTextEditor);

////////////////////////////////////////////////////////////////////////////////
/// TGTextEditor constructor with file name as first argument.

TGTextEditor::TGTextEditor(const char *filename, const TGWindow *p, UInt_t w,
                           UInt_t h) : TGMainFrame(p, w, h)
{
   Build();
   if (p && p != gClient->GetDefaultRoot()) {
      // special case for TRootBrowser
      // remove the command line combo box and its associated label
      fComboCmd->UnmapWindow();
      fToolBar->RemoveFrame(fComboCmd);
      fLabel->UnmapWindow();
      fToolBar->RemoveFrame(fLabel);
      fToolBar->GetButton(kM_FILE_EXIT)->SetState(kButtonDisabled);
      fToolBar->Layout();
   }
   if (filename) {
      LoadFile(filename);
   }
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// TGTextEditor constructor with pointer to a TMacro as first argument.

TGTextEditor::TGTextEditor(TMacro *macro, const TGWindow *p, UInt_t w, UInt_t h) :
              TGMainFrame(p, w, h)
{
   TString tmp;
   Build();
   if (p && p != gClient->GetDefaultRoot()) {
      // special case for TRootBrowser
      // remove the command line combo box and its associated label
      fComboCmd->UnmapWindow();
      fLabel->UnmapWindow();
      if (fToolBar->GetButton(kM_FILE_EXIT))
         fToolBar->GetButton(kM_FILE_EXIT)->SetState(kButtonDisabled);
      fToolBar->Layout();
   }
   if (macro) {
      fMacro = macro;
      TIter next(macro->GetListOfLines());
      TObjString *obj;
      while ((obj = (TObjString*) next())) {
         fTextEdit->AddLine(obj->GetName());
      }
      tmp.Form("TMacro : %s: %ld lines read.",
               macro->GetName(), fTextEdit->ReturnLineCount());
      fStatusBar->SetText(tmp.Data(), 0);
      fFilename = macro->GetName();
      fFilename += ".C";
      tmp.Form("TMacro : %s - TGTextEditor", macro->GetName());
      SetWindowName(tmp.Data());
   }
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// TGTextEditor destructor.

TGTextEditor::~TGTextEditor()
{
   gApplication->Disconnect("Terminate(Int_t)");
   if (fTimer) delete fTimer;
   if (fMenuFile) delete fMenuFile;
   if (fMenuEdit) delete fMenuEdit;
   if (fMenuSearch) delete fMenuSearch;
   if (fMenuTools) delete fMenuTools;
   if (fMenuHelp) delete fMenuHelp;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete TGTextEditor Window.

void TGTextEditor::DeleteWindow()
{
   gApplication->Disconnect("Terminate(Int_t)");
   delete fTimer; fTimer = 0;
   delete fMenuFile; fMenuFile = 0;
   delete fMenuEdit; fMenuEdit = 0;
   delete fMenuSearch; fMenuSearch = 0;
   delete fMenuTools; fMenuTools = 0;
   delete fMenuHelp; fMenuHelp = 0;
   Cleanup();
   TGMainFrame::DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Build TGTextEditor widget.

void TGTextEditor::Build()
{
   SetCleanup(kDeepCleanup);
   fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 1, 1);
   fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);

   fMenuFile = new TGPopupMenu(fClient->GetDefaultRoot());
   fMenuFile->AddEntry("&New", kM_FILE_NEW);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("&Open...", kM_FILE_OPEN);
   fMenuFile->AddEntry("&Close", kM_FILE_CLOSE);
   fMenuFile->AddEntry("&Save", kM_FILE_SAVE);
   fMenuFile->AddEntry("Save &As...", kM_FILE_SAVEAS);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("&Print...", kM_FILE_PRINT);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("E&xit", kM_FILE_EXIT);

   fMenuEdit = new TGPopupMenu(fClient->GetDefaultRoot());
   fMenuEdit->AddEntry("Cu&t\tCtrl+X", kM_EDIT_CUT);
   fMenuEdit->AddEntry("&Copy\tCtrl+C", kM_EDIT_COPY);
   fMenuEdit->AddEntry("&Paste\tCtrl+V", kM_EDIT_PASTE);
   fMenuEdit->AddEntry("De&lete\tDel", kM_EDIT_DELETE);
   fMenuEdit->AddSeparator();
   fMenuEdit->AddEntry("Select &All\tCtrl+A", kM_EDIT_SELECTALL);
   fMenuEdit->AddSeparator();
   fMenuEdit->AddEntry("Set &Font", kM_EDIT_SELFONT);

   fMenuTools = new TGPopupMenu(fClient->GetDefaultRoot());
   fMenuTools->AddEntry("&Compile Macro\tCtrl+F7", kM_TOOLS_COMPILE);
   fMenuTools->AddEntry("&Execute Macro\tCtrl+F5", kM_TOOLS_EXECUTE);
   fMenuTools->AddEntry("&Interrupt\tShift+F5", kM_TOOLS_INTERRUPT);

   fMenuEdit->DisableEntry(kM_EDIT_CUT);
   fMenuEdit->DisableEntry(kM_EDIT_COPY);
   fMenuEdit->DisableEntry(kM_EDIT_DELETE);
   fMenuEdit->DisableEntry(kM_EDIT_PASTE);

   fMenuSearch = new TGPopupMenu(fClient->GetDefaultRoot());
   fMenuSearch->AddEntry("&Find...\tCtrl+F", kM_SEARCH_FIND);
   fMenuSearch->AddEntry("Find &Next\tF3", kM_SEARCH_FINDNEXT);
   fMenuSearch->AddSeparator();
   fMenuSearch->AddEntry("&Goto Line...\tCtrl+L", kM_SEARCH_GOTO);

   fMenuHelp = new TGPopupMenu(fClient->GetDefaultRoot());
   fMenuHelp->AddEntry("&Help Topics\tF1", kM_HELP_CONTENTS);
   fMenuHelp->AddSeparator();
   fMenuHelp->AddEntry("&About...", kM_HELP_ABOUT);

   fMenuFile->Associate(this);
   fMenuEdit->Associate(this);
   fMenuSearch->Associate(this);
   fMenuTools->Associate(this);
   fMenuHelp->Associate(this);

   fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);
   fMenuBar->SetCleanup(kDeepCleanup);
   fMenuBar->AddPopup("&File", fMenuFile, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Edit", fMenuEdit, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Search", fMenuSearch, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Tools", fMenuTools, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Help", fMenuHelp, new TGLayoutHints(kLHintsTop |
                      kLHintsRight));
   AddFrame(fMenuBar, fMenuBarLayout);

   //---- toolbar

   AddFrame(new TGHorizontal3DLine(this),
            new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0,0,2,2));
   Int_t i,spacing = 8;
   fToolBar = new TGToolBar(this, 60, 20, kHorizontalFrame);
   fToolBar->SetCleanup(kDeepCleanup);
   for (i = 0; fTbData[i].fPixmap; i++) {
      if (strlen(fTbData[i].fPixmap) == 0) {
         spacing = 8;
         continue;
      }
      fToolBar->AddButton(this, &fTbData[i], spacing);
      spacing = 0;
   }
   fComboCmd   = new TGComboBox(fToolBar, "");
   fCommand    = fComboCmd->GetTextEntry();
   fCommandBuf = fCommand->GetBuffer();
   fCommand->Associate(this);
   fComboCmd->Resize(200, fCommand->GetDefaultHeight());
   fToolBar->AddFrame(fComboCmd, new TGLayoutHints(kLHintsCenterY |
            kLHintsRight, 5, 5, 1, 1));

   fToolBar->AddFrame(fLabel = new TGLabel(fToolBar, "Command :"),
            new TGLayoutHints(kLHintsCenterY | kLHintsRight, 5, 5, 1, 1));
   AddFrame(fToolBar, new TGLayoutHints(kLHintsTop | kLHintsExpandX,
            0, 0, 0, 0));
   AddFrame(new TGHorizontal3DLine(this),
            new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0,0,2,2));

   fToolBar->GetButton(kM_EDIT_CUT)->SetState(kButtonDisabled);
   fToolBar->GetButton(kM_EDIT_COPY)->SetState(kButtonDisabled);
   fToolBar->GetButton(kM_EDIT_DELETE)->SetState(kButtonDisabled);
   fToolBar->GetButton(kM_EDIT_PASTE)->SetState(kButtonDisabled);

   fTextEdit = new TGTextEdit(this, 10, 10, 1);
   if (gClient->GetStyle() < 2) {
      Pixel_t pxl;
      gClient->GetColorByName("#3399ff", pxl);
      fTextEdit->SetSelectBack(pxl);
      fTextEdit->SetSelectFore(TGFrame::GetWhitePixel());
   }
   fTextEdit->Associate(this);
   AddFrame(fTextEdit, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   Int_t parts[] = { 75, 25 };
   fStatusBar = new TGStatusBar(this);
   fStatusBar->SetCleanup(kDeepCleanup);
   fStatusBar->SetParts(parts, 2);
   AddFrame(fStatusBar, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 0, 0, 3, 0));

   SetClassHints("ROOT", "TGTextEditor");
   SetWindowName("Untitled - TGTextEditor");

   fMacro = 0;
   fFilename = "Untitled";
   fStatusBar->SetText(fFilename.Data(), 0);

   fTextEdit->SetFocus();
   fTextEdit->GetMenu()->DisableEntry(TGTextEdit::kM_FILE_NEW);
   fTextEdit->GetMenu()->DisableEntry(TGTextEdit::kM_FILE_OPEN);
   fTextEdit->Connect("DataChanged()", "TGTextEditor", this, "DataChanged()");
   fTextEdit->Connect("Closed()", "TGTextEditor", this, "ClearText()");
   fTextEdit->Connect("Opened()", "TGTextEditor", this, "ClearText()");
   fTextEdit->Connect("DataDropped(char *)", "TGTextEditor", this, "DataDropped(char *)");
   fTextEdit->MapWindow();

   MapSubwindows();
   Resize(GetDefaultWidth() + 50, GetDefaultHeight() > 500 ? GetDefaultHeight() : 500);
   Layout();

   gApplication->Connect("Terminate(Int_t)", "TGTextEditor", this, "ClearText()");
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_F3), 0, kTRUE);

   AddInput(kKeyPressMask | kEnterWindowMask | kLeaveWindowMask |
            kFocusChangeMask | kStructureNotifyMask);

   fTimer = new TTimer(this, 250);
   fTimer->Reset();
   fTimer->TurnOn();

   fExiting = kFALSE;
   fTextChanged = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Update file information when receiving the signal
/// DataDropped from TGTextEdit widget.

void TGTextEditor::DataDropped(char *fname)
{
   TString tmp;
   fFilename = fname;
   tmp.Form("%s: %ld lines read.", fname, fTextEdit->ReturnLineCount());
   fStatusBar->SetText(tmp, 0);
   tmp.Form("%s - TGTextEditor", fname);
   SetWindowName(tmp.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Slot setting the fTextChanged flag to true when the text has been
/// modified in fTextEdit

void TGTextEditor::DataChanged()
{
   fTextChanged = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Load a file into the editor. If fname is 0, a TGFileDialog will popup.

void TGTextEditor::LoadFile(const char *fname)
{
   TString tmp;
   TGFileInfo fi;
   fi.fFileTypes = ed_filetypes;
   switch (IsSaved()) {
      case kMBCancel:
         return;
      case kMBYes:
         if (!fFilename.CompareTo("Untitled"))
            SaveFileAs();
         else
            SaveFile(fFilename.Data());
         if (fTextChanged) {
            return;
         }
         break;
      case kMBNo:
         break;
      default:
         return;
   }
   if (fname == 0) {
      new TGFileDialog(fClient->GetDefaultRoot(), this, kFDOpen, &fi);
      if (fi.fFilename && strlen(fi.fFilename)) {
         fname = fi.fFilename;
      }
   }
   if (fname) {
      if (!fTextEdit->LoadFile(fname)) {
         tmp.Form("Error opening file \"%s\"", fname);
         new TGMsgBox(fClient->GetRoot(), this, "TGTextEditor",
                      tmp.Data(), kMBIconExclamation, kMBOk);
      } else {
         fFilename = fname;
         tmp.Form("%s: %ld lines read.", fname, fTextEdit->ReturnLineCount());
         fStatusBar->SetText(tmp.Data(), 0);
         tmp.Form("%s - TGTextEditor", fname);
         SetWindowName(tmp.Data());
         fTextChanged = kFALSE;
      }
   }
   fTextEdit->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Save the edited text in the file "fname".

void TGTextEditor::SaveFile(const char *fname)
{
   char *p;
   TString tmp;

   if (!fTextEdit->SaveFile(fname)) {
      tmp.Form("Error saving file \"%s\"", fname);
      new TGMsgBox(fClient->GetRoot(), this, "TGTextEditor",
                   tmp.Data(), kMBIconExclamation, kMBOk);
      return;
   }
   if ((p = (char *)strrchr(fname, '/')) == 0) {
      p = (char *)fname;
   } else {
      ++p;
   }
   tmp.Form("%s: %ld lines written.", p, fTextEdit->ReturnLineCount());
   fStatusBar->SetText(tmp.Data(), 0);

   tmp.Form("%s - TGTextEditor", p);
   SetWindowName(tmp.Data());
   fTextChanged = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Save the edited text in a file selected with TGFileDialog.
/// Shouldn't we create a backup file?

Bool_t TGTextEditor::SaveFileAs()
{
   TString workdir = gSystem->WorkingDirectory();
   static TString dir(".");
   static Bool_t overwr = kFALSE;
   TGFileInfo fi;
   fi.fFileTypes = ed_filetypes;
   fi.SetIniDir(dir);
   fi.fOverwrite = overwr;
   new TGFileDialog(fClient->GetDefaultRoot(), this, kFDSave, &fi);
   gSystem->ChangeDirectory(workdir.Data());
   overwr = fi.fOverwrite;
   if (fi.fFilename && strlen(fi.fFilename)) {
      SaveFile(fi.fFilename);
      fFilename = fi.fFilename;
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if file has to be saved in case of modifications.

Int_t TGTextEditor::IsSaved()
{
   Int_t ret;
   TString tmp;
   Int_t opt = (kMBYes | kMBNo);

   tmp.Form("The text has been modified. Do you want to save the changes?");

   if (!fTextChanged) {
      return kMBNo;
   } else {
      if (fParent == gClient->GetDefaultRoot())
         opt |= kMBCancel;
      new TGMsgBox(fClient->GetRoot(), this, "TGTextEditor",
                   tmp.Data(), kMBIconExclamation, opt, &ret);
      return ret;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Open the print dialog and send current buffer to printer.

void TGTextEditor::PrintText()
{
   TString tmp;
   Int_t ret = 0;
   if (!gEPrinter) {
      gEPrinter = StrDup("892_2_cor"); // use gEnv
      gEPrintCommand = StrDup("xprint");
   }
   new TGPrintDialog(fClient->GetDefaultRoot(), this, 400, 150,
                     &gEPrinter, &gEPrintCommand, &ret);
   if (ret) {
      fTextEdit->Print();
      tmp.Form("Printed: %s", fFilename.Data());
      fStatusBar->SetText(tmp.Data(), 0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Close TGTextEditor window.

void TGTextEditor::CloseWindow()
{
   Bool_t reallyClose = kFALSE;
   if (fExiting) {
      return;
   }
   fExiting = kTRUE;
   switch (IsSaved()) {
      case kMBYes:
         if (!fFilename.CompareTo("Untitled"))
            fTextChanged = !SaveFileAs();
         else
            SaveFile(fFilename.Data());
         // close the editor only if the text has been saved and the
         // TGTextEditor is not embedded in the ROOT browser
         reallyClose = (!fTextChanged || (fParent != gClient->GetDefaultRoot()));
         break;
      case kMBCancel:
         reallyClose = kFALSE;
         break;
      case kMBNo:
         reallyClose = kTRUE;
         break;
   }
   if (reallyClose) {
      gApplication->Disconnect("Terminate(Int_t)");
      TGMainFrame::CloseWindow();
   }
   fExiting = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Keyboard event handler.

Bool_t TGTextEditor::HandleKey(Event_t *event)
{
   char   input[10];
   UInt_t keysym;

   if (event->fType == kGKeyPress) {
      gVirtualX->LookupString(event, input, sizeof(input), keysym);

      switch ((EKeySym)keysym) {   // ignore these keys
         case kKey_Shift:
         case kKey_Control:
         case kKey_Meta:
         case kKey_Alt:
         case kKey_CapsLock:
         case kKey_NumLock:
         case kKey_ScrollLock:
            return kTRUE;
         case kKey_F1:
            SendMessage(this, MK_MSG(kC_COMMAND, kCM_MENU),
                        kM_HELP_CONTENTS, 0);
            return kTRUE;
         case kKey_F3:
            Search(kTRUE);
            return kTRUE;
         default:
            break;
      }
      if (event->fState & kKeyControlMask) {   // Ctrl key modifier pressed
         switch((EKeySym)keysym) {
            case kKey_F5:
               ExecuteMacro();
               return kTRUE;
            case kKey_F7:
               CompileMacro();
               return kTRUE;
            default:
               break;
         }
      }
      if (event->fState & kKeyShiftMask) {   // Shift key modifier pressed
         switch((EKeySym)keysym) {
            case kKey_F5:
               InterruptMacro();
               return kTRUE;
            default:
               break;
         }
      }
   }
   return TGMainFrame::HandleKey(event);
}

////////////////////////////////////////////////////////////////////////////////
/// Clear text edit widget.

void TGTextEditor::ClearText()
{
   fTextEdit->Clear();
   fMacro = 0;
   fFilename = "Untitled";
   SetWindowName("Untitled - TGTextEditor");
   fStatusBar->SetText("New File", 0);
   fTextChanged = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Invokes search dialog, or just search previous string if again is true.

void TGTextEditor::Search(Bool_t again)
{
   if (again) {
      SendMessage(fTextEdit, MK_MSG(kC_COMMAND, kCM_MENU),
                  TGTextEdit::kM_SEARCH_FINDAGAIN, 0);
   }
   else {
      fTextEdit->Search(kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Invokes goto dialog, and go to the specified line.

void TGTextEditor::Goto()
{
   Long_t ret;

   new TGGotoDialog(fClient->GetDefaultRoot(), this, 400, 150, &ret);

   if (ret >= 0)
      fTextEdit->Goto(ret-1);
}

////////////////////////////////////////////////////////////////////////////////
/// Save the edited text in a temporary macro, then compile it.

void TGTextEditor::CompileMacro()
{
   if (fTextEdit->ReturnLineCount() < 3)
      return;
   if ((fMacro) || (!fFilename.CompareTo("Untitled"))) {
      if (!SaveFileAs())
         return;
   }
   char *tmpfile = gSystem->ConcatFileName(gSystem->TempDirectory(),
                                gSystem->BaseName(fFilename.Data()));
   fTextEdit->SaveFile(tmpfile, kFALSE);
   gSystem->CompileMacro(tmpfile);
   gSystem->Unlink(tmpfile);
   delete [] tmpfile;
}

////////////////////////////////////////////////////////////////////////////////
/// Save the edited text in a temporary macro, execute it, and then delete
/// the temporary file.

void TGTextEditor::ExecuteMacro()
{
   if (fTextEdit->ReturnLineCount() < 3)
      return;
   if (fMacro) {
      fMacro->Exec();
      return;
   }
   if (fTextChanged) {
      Int_t ret;
      new TGMsgBox(fClient->GetRoot(), this, "TGTextEditor",
            "The text has been modified. Do you want to save the changes?",
            kMBIconExclamation, kMBYes | kMBNo | kMBCancel, &ret);
      if (ret == kMBYes) {
         if (!fFilename.CompareTo("Untitled"))
            SaveFileAs();
         else
            SaveFile(fFilename.Data());
         fTextChanged = kFALSE;
      }
      if (ret == kMBCancel)
         return;
   }
   if (!fFilename.CompareTo("Untitled")) {
      //if (!SaveFileAs())
      //   return;
      fFilename += ".C";
   }
   gInterpreter->SaveContext();
   TString savdir = gSystem->WorkingDirectory();
   TString tmpfile = gSystem->BaseName(fFilename.Data());
   tmpfile += "_exec";
   gSystem->ChangeDirectory(gSystem->GetDirName(fFilename.Data()).Data());
   fTextEdit->SaveFile(tmpfile.Data(), kFALSE);
   gROOT->SetExecutingMacro(kTRUE);
   gROOT->Macro(tmpfile.Data());
   gROOT->SetExecutingMacro(kFALSE);
   if (gInterpreter->IsLoaded(tmpfile.Data()))
      gInterpreter->UnloadFile(tmpfile.Data());
   gSystem->Unlink(tmpfile.Data());
   gSystem->ChangeDirectory(savdir.Data());
   gInterpreter->Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Interrupt execution of a macro.

void TGTextEditor::InterruptMacro()
{
   gROOT->SetInterrupt(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Display ROOT splash screen.

void TGTextEditor::About()
{
#ifdef R__UNIX
   TString rootx = TROOT::GetBinDir() + "/root -a &";
   gSystem->Exec(rootx);
#else
#ifdef WIN32
   new TWin32SplashThread(kTRUE);
#else
   char str[32];
   sprintf(str, "About ROOT %s...", gROOT->GetVersion());
   TRootHelpDialog *hd = new TRootHelpDialog(this, str, 600, 400);
   hd->SetText(gHelpAbout);
   hd->Popup();
#endif
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Handle timer event.

Bool_t TGTextEditor::HandleTimer(TTimer *t)
{
   TString tmp;
   if (t != fTimer) return kTRUE;
   // check if some text is available in the clipboard
   if ((gVirtualX->InheritsFrom("TGX11")) &&
      (gVirtualX->GetPrimarySelectionOwner() == kNone)) {
      fMenuEdit->DisableEntry(kM_EDIT_PASTE);
      fToolBar->GetButton(kM_EDIT_PASTE)->SetState(kButtonDisabled);
   }
   else {
      fMenuEdit->EnableEntry(kM_EDIT_PASTE);
      if (fToolBar->GetButton(kM_EDIT_PASTE)->GetState() == kButtonDisabled)
         fToolBar->GetButton(kM_EDIT_PASTE)->SetState(kButtonUp);
   }
   // check if text is selected in the editor
   if (fTextEdit->IsMarked()) {
      fMenuEdit->EnableEntry(kM_EDIT_CUT);
      fMenuEdit->EnableEntry(kM_EDIT_COPY);
      fMenuEdit->EnableEntry(kM_EDIT_DELETE);
      if (fToolBar->GetButton(kM_EDIT_CUT)->GetState() == kButtonDisabled) {
         fToolBar->GetButton(kM_EDIT_CUT)->SetState(kButtonUp);
         fToolBar->GetButton(kM_EDIT_COPY)->SetState(kButtonUp);
         fToolBar->GetButton(kM_EDIT_DELETE)->SetState(kButtonUp);
      }
   }
   else {
      fMenuEdit->DisableEntry(kM_EDIT_CUT);
      fMenuEdit->DisableEntry(kM_EDIT_COPY);
      fMenuEdit->DisableEntry(kM_EDIT_DELETE);
      if (fToolBar->GetButton(kM_EDIT_CUT)->GetState() == kButtonUp) {
         fToolBar->GetButton(kM_EDIT_CUT)->SetState(kButtonDisabled);
         fToolBar->GetButton(kM_EDIT_COPY)->SetState(kButtonDisabled);
         fToolBar->GetButton(kM_EDIT_DELETE)->SetState(kButtonDisabled);
      }
   }
   // get cursor position
   TGLongPosition pos = fTextEdit->GetCurrentPos();
   tmp.Form("Ln %ld, Ch %ld", pos.fY, pos.fX);
   fStatusBar->SetText(tmp.Data(), 1);
   fTimer->Reset();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle menu and other command generated by the user.

Bool_t TGTextEditor::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   TRootHelpDialog *hd;

   switch(GET_MSG(msg)) {
      case kC_COMMAND:
         switch(GET_SUBMSG(msg)) {
            case kCM_BUTTON:
            case kCM_MENU:
               switch (parm1) {
                  // "File" menu related events
                  case kM_FILE_NEW:
                     new TGTextEditor();
                     break;
                  case kM_FILE_OPEN:
                     LoadFile();
                     break;
                  case kM_FILE_CLOSE:
                     ClearText();
                     break;
                  case kM_FILE_SAVE:
                     if (!fFilename.CompareTo("Untitled"))
                        SaveFileAs();
                     else
                        SaveFile(fFilename.Data());
                     break;
                  case kM_FILE_SAVEAS:
                     SaveFileAs();
                     break;
                  case kM_FILE_PRINT:
                     PrintText();
                     break;
                  case kM_FILE_EXIT:
                     CloseWindow();
                     break;

                  // "Edit" menu related events
                  case kM_EDIT_CUT:
                     fTextEdit->Cut();
                     break;
                  case kM_EDIT_COPY:
                     fTextEdit->Copy();
                     break;
                  case kM_EDIT_PASTE:
                     fTextEdit->Paste();
                     break;
                  case kM_EDIT_DELETE:
                     fTextEdit->Delete();
                     break;
                  case kM_EDIT_SELECTALL:
                     fTextEdit->SelectAll();
                     if (fTextEdit->IsMarked()) {
                        fMenuEdit->EnableEntry(kM_EDIT_CUT);
                        fMenuEdit->EnableEntry(kM_EDIT_COPY);
                        fMenuEdit->EnableEntry(kM_EDIT_DELETE);
                        if (fToolBar->GetButton(kM_EDIT_CUT)->GetState() == kButtonDisabled) {
                           fToolBar->GetButton(kM_EDIT_CUT)->SetState(kButtonUp);
                           fToolBar->GetButton(kM_EDIT_COPY)->SetState(kButtonUp);
                           fToolBar->GetButton(kM_EDIT_DELETE)->SetState(kButtonUp);
                        }
                     }
                     break;
                  case kM_EDIT_SELFONT:
                     {
                        Int_t count;
                        TString fontname;
                        TGFontDialog::FontProp_t prop;
                        new TGFontDialog(fClient->GetRoot(), this, &prop);
                        if (prop.fName != "") {
                           fontname.Form("-*-%s-%s-%c-*-*-%d-*-*-*-*-*-*-*",
                                         prop.fName.Data(),
                                         prop.fBold ? "bold" : "medium",
                                         prop.fItalic ? 'i' : 'r',
                                         prop.fSize);
                           if (!gVirtualX->ListFonts(fontname, 10, count)) {
                              fontname.Form("-*-%s-%s-%c-*-*-%d-*-*-*-*-*-*-*",
                                            prop.fName.Data(),
                                            prop.fBold ? "bold" : "medium",
                                            prop.fItalic ? 'o' : 'r',
                                            prop.fSize);
                           }
                           TGFont *font = fClient->GetFont(fontname);
                           if (font) {
                              FontStruct_t editorfont = font->GetFontStruct();
                              fTextEdit->SetFont(editorfont);
                              fTextEdit->Update();
                           }
                        }
                     }
                     break;

                  // "Tools" menu related events
                  case kM_TOOLS_COMPILE:
                     CompileMacro();
                     break;
                  case kM_TOOLS_EXECUTE:
                     ExecuteMacro();
                     break;
                  case kM_TOOLS_INTERRUPT:
                     InterruptMacro();
                     break;

                  // "Search" menu related events
                  case kM_SEARCH_FIND:
                     Search(kFALSE);
                     break;
                  case kM_SEARCH_FINDNEXT:
                     Search(kTRUE);
                     break;
                  case kM_SEARCH_GOTO:
                     Goto();
                     break;

                  // "Help" menu related events
                  case kM_HELP_CONTENTS:
                     hd = new TRootHelpDialog(this, "Help on Editor...", 600, 400);
                     hd->SetText(gHelpTextEditor);
                     hd->Popup();
                     break;
                  case kM_HELP_ABOUT:
                     About();
                     break;
               }
               break;
         }
         break;
      case kC_TEXTENTRY:
         switch (GET_SUBMSG(msg)) {
            case kTE_ENTER:
               {
                  // here copy the string from text buffer to return variable
                  const char *string = fCommandBuf->GetString();
                  if(strlen(string) > 1) {
                     gROOT->ProcessLine(string);
                     fComboCmd->ReturnPressed();
                  }
               }
               break;
            default:
               break;
         }
         break;

      default:
         break;
   }
   return kTRUE;
}
