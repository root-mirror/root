// @(#)root/test/RootIDE/:$Id$
// Author: Bertrand Bellenot   20/04/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TGRootIDE                                                           //
//                                                                      //
//  A simple IDE editor that uses the TGTextEdit widget.                //
//  It provides all functionalities of TGTextEdit as copy, paste, cut,  //
//  search, go to a given line number. In addition, it provides the     //
//  possibilities for compiling, executing or interrupting a running    //
//  macro.                                                              //
//                                                                      //
//  This class can be used in following ways:                           //
//  - with file name as argument:                                       //
//    new TGRootIDE("hsimple.C");                                       //
//  - with a TMacro* as argument:                                       //
//    TMacro *macro = new TMacro("hsimple.C");                          //
//    new TGRootIDE(macro);                                             //
//                                                                      //
//  Basic Features:                                                     //
//                                                                      //
//  New Document                                                        //
//                                                                      //
//  To create a new blank document, select File menu / New, or click    //
//  the New toolbar button. It will create a new instance of            //
//  TGRootIDE.                                                          //
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
<img src="gif/TGRootIDE.gif">
*/
//End_Html
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TMacro.h"
#include "TGMsgBox.h"
#include "TGFileDialog.h"
#include "TGFontDialog.h"
#include "TGTextEdit.h"
#include "TGMenu.h"
#include "TGButton.h"
#include "TGStatusBar.h"
#include "KeySymbols.h"
#include "TGToolBar.h"
#include "TG3DLine.h"
#include "TGLabel.h"
#include "TGTextEntry.h"
#include "TGTextEditDialogs.h"
#include "TGRootIDE.h"
#include "TGComboBox.h"
#include "TGTab.h"
#include "TGFSContainer.h"
#include "TGListView.h"
#include "TBrowser.h"
#include "TFile.h"
#include "TKey.h"
#include "TObjString.h"
#include "TRootHelpDialog.h"
#include "TGSplitter.h"
#include "TObjArray.h"
#include "HelpText.h"
#include "TGHtml.h"
#include "TUrl.h"
#include "TSocket.h"
#include "TImage.h"
#include "THtml.h"
#include "TRint.h"
#include "TProcessID.h"
#include "Getline.h"
#ifdef WIN32
#include "TWin32SplashThread.h"
#endif
#include <string>

const char *ed_filetypes[] = {
   "ROOT Macros",  "*.C",
   "Source files", "*.cxx",
   "Text files",   "*.txt",
   "All files",    "*",
   0, 0
};

const char *filters[] = {
   "",
   "*.*",
   "*.[C|c|h]*",
   "*.txt"
};

const char *HtmlError[] = {
"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3c.org/TR/1999/REC-html401-19991224/loose.dtd\"> ",
"<HTML><HEAD><TITLE>RHTML cannot display the webpage</TITLE> ",
"<META http-equiv=Content-Type content=\"text/html; charset=UTF-8\"></HEAD> ",
"<BODY> ",
"<TABLE cellSpacing=0 cellPadding=0 width=730 border=0> ",
"  <TBODY> ",
"  <TR> ",
"    <TD id=infoIconAlign vAlign=top align=left width=60 rowSpan=2> ",
"    <IMG src=\"info.gif\"> ",
"    </TD> ",
"    <TD id=mainTitleAlign vAlign=center align=left width=*> ",
"      <H1 id=mainTitle>RHTML cannot display the webpage</H1></TD></TR> ",
"  <TR> ",
"    <TD class=errorCodeAndDivider id=errorCodeAlign align=right>&nbsp;  ",
"      <DIV class=divider></DIV></TD></TR> ",
"  <TR> ",
"      <UL> ",
"      </UL> ",
"    <TD>&nbsp; </TD> ",
"    <TD id=MostLikelyAlign vAlign=top align=left> ",
"      <H3 id=likelyCauses>Most likely causes:</H3> ",
"      <UL> ",
"        <LI id=causeNotConnected>You are not connected to the Internet.  ",
"        <LI id=causeSiteProblem>The website is encountering problems.  ",
"        <LI id=causeErrorInAddress>There might be a typing error in the address.  ",
"        <LI id=causeOtherError>  ",
"        </LI></UL></TD></TR> ",
"  <TR> ",
"    <TD id=infoBlockAlign vAlign=top align=right>&nbsp; </TD> ",
"    <TD id=moreInformationAlign vAlign=center align=left> ",
"      <H4> ",
"      <TABLE> ",
"        <TBODY> ",
"        <TR> ",
"          <TD vAlign=top><SPAN id=moreInfoContainer></SPAN><ID  ",
"            id=moreInformation>More information</ID> ",
"      </TD></TR></TBODY></TABLE></H4> ",
"      <DIV class=infoBlock id=infoBlockID> ",
"      <P><ID id=errorExpl1>This problem can be caused by a variety of issues,  ",
"      including:</ID>  ",
"      <UL> ",
"        <LI id=errorExpl2>Internet connectivity has been lost.  ",
"        <LI id=errorExpl3>The website is temporarily unavailable.  ",
"        <LI id=errorExpl4>The Domain Name Server (DNS) is not reachable.  ",
"        <LI id=errorExpl5>The Domain Name Server (DNS) does not have a listing  ",
"        for the website's domain.  ",
"      <P></P> ",
"      <P></P></DIV></TD></TR></TBODY></TABLE></BODY></HTML> ",
0
};

enum ETextEditorCommands {
   kM_FILE_NEW, kM_FILE_CLOSE, kM_FILE_OPEN, kM_FILE_SAVE, kM_FILE_SAVEAS, kM_FILE_PRINT,
   kM_FILE_EXIT, kM_EDIT_CUT, kM_EDIT_COPY, kM_EDIT_PASTE, kM_EDIT_DELETE,
   kM_EDIT_SELECTALL, kM_SEARCH_FIND, kM_SEARCH_FINDNEXT, kM_SEARCH_GOTO,
   kM_TOOLS_COMPILE, kM_TOOLS_EXECUTE, kM_TOOLS_INTERRUPT, kM_TOOLS_BROWSER,
   kM_TOOLS_CLEAN_LOG, kM_HELP_CONTENTS, kM_HELP_ABOUT, kM_EDIT_SELFONT
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

ClassImp(TGRootIDE);


////////////////////////////////////////////////////////////////////////////////

TGDocument::TGDocument(const char *fname, const char *title, Int_t tabid,
                     TGTab *tab, TGTabElement *tabel, TGTextEdit *edit,
                     TObjArray *doclist) :  TNamed(fname, title)
{
   fModified = kFALSE;
   fTabId    = tabid;
   fEditor   = edit;
   fTab      = tab;
   fTabEl    = tabel;
   fDocList  = doclist;
   Open(fname);
}

////////////////////////////////////////////////////////////////////////////////
/// Close the current active document.

Bool_t TGDocument::Close()
{
   Int_t ret;
   if (fModified) {
      // the current active document has been modified
      // then ask the user if they want to save it
      TString sfname(GetName());
      new TGMsgBox(gClient->GetRoot(), fTab, "TGRootIDE",
                   Form("%s has been modified. Do you want to save the changes?",
                        sfname.Data()),
                   kMBIconExclamation, kMBYes | kMBNo | kMBCancel, &ret);
      if (ret == kMBYes) {
         Save();
      }
      if (ret != kMBCancel) {
         // always keep the first two tabs
         if (fTab->GetNumberOfTabs() > 2) {
            fTab->RemoveTab(fTab->GetCurrent());
            fTab->Layout();
         }
         else {
            fEditor->Clear();
            fTabEl->SetText(new TGString("Untitled"));
            fTab->MapSubwindows();
            fTab->Layout();
         }
         fDocList->Remove((TObject *)this);
         ((TGRootIDE *)fTab->GetMainFrame())->DoTab(fTab->GetCurrent());
         delete this;
         return kTRUE;
      }
   }
   else {
      // no changes, so just close it
      if (fTab->GetNumberOfTabs() > 2) {
         // always keep the first two tabs
         fTab->RemoveTab(fTab->GetCurrent());
         fTab->Layout();
         fDocList->Remove((TObject *)this);
         ((TGRootIDE *)fTab->GetMainFrame())->DoTab(fTab->GetCurrent());
         delete this;
         return kTRUE;
      }
      else {
         fEditor->Clear();
         fTabEl->SetText(new TGString("Untitled"));
         fTab->MapSubwindows();
         fTab->Layout();
      }
      ((TGRootIDE *)fTab->GetMainFrame())->DoTab(fTab->GetCurrent());
      return kTRUE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Open file and create new document. If fname is NULL, create new Untitled
/// document. Add a new tab element for this documment.

Bool_t TGDocument::Open(const char *fname)
{
   TGFileInfo fi;
   fi.fFileTypes = ed_filetypes;
   if (fname == 0) {
      // no filename provided --> empty (untitled) document
      TGCompositeFrame *tf = fTab->AddTab("Untitled");
      tf->SetLayoutManager(new TGHorizontalLayout(tf));
      TGTextEdit *textEdit = new TGTextEdit(tf, 10, 10, 1);
      Pixel_t pxl;
      gClient->GetColorByName("#ccccff", pxl);
      textEdit->SetSelectBack(pxl);
      textEdit->SetSelectFore(TGFrame::GetBlackPixel());
      textEdit->Connect("DataChanged()", "TGDocument", this, "DataChanged()");
      textEdit->Connect("DataDropped(char *)", "TGDocument", this, "DataDropped(char *)");
      tf->AddFrame(textEdit, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
      fEditor = textEdit;
      fEditor->Associate(fTab->GetMainFrame());
      fEditor->MapWindow();
      fTab->SetTab(fTab->GetNumberOfTabs()-1, kFALSE);
      fTab->MapSubwindows();
      fTab->Layout();
      fTabEl = fTab->GetTabTab(fTab->GetCurrent());
      fTabEl->ShowClose();
   }
   if (fname) {
      if ((fEditor == 0) || (fTabEl == 0) || (fEditor &&
           fEditor->GetText()->RowCount() > 1 &&
           fEditor->GetText()->ColCount() > 1)) {
         // if no current editor, or if current text editor already has
         // text, add a new tab
         TGCompositeFrame *tf = fTab->AddTab(gSystem->BaseName(fname));
         tf->SetLayoutManager(new TGHorizontalLayout(tf));
         TGTextEdit *textEdit = new TGTextEdit(tf, 10, 10, 1);
         Pixel_t pxl;
         gClient->GetColorByName("#ccccff", pxl);
         textEdit->SetSelectBack(pxl);
         textEdit->SetSelectFore(TGFrame::GetBlackPixel());
         textEdit->Connect("DataChanged()", "TGDocument", this, "DataChanged()");
         textEdit->Connect("DataDropped(char *)", "TGDocument", this, "DataDropped(char *)");
         tf->AddFrame(textEdit, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
         fEditor = textEdit;
         fEditor->Associate(fTab->GetMainFrame());
         fEditor->MapWindow();
         fTab->SetTab(fTab->GetNumberOfTabs()-1, kFALSE);
         fTab->MapSubwindows();
         fTab->Layout();
         fTabEl = fTab->GetTabTab(fTab->GetCurrent());
         fTabEl->ShowClose();
      }
      if (!fname[0]) {
         // no filename provided --> empty (untitled) document
         SetName("Untitled");
         SetTitle("Untitled");
         fEditor->Layout();
         fEditor->Connect("DataChanged()", "TGDocument", this, "DataChanged()");
         fEditor->Connect("DataDropped(char *)", "TGDocument", this, "DataDropped(char *)");
         fTabEl->SetText(new TGString("Untitled"));
         fTab->MapSubwindows();
         fTab->Layout();
      }
      else {
         // current tab is empty --> open document inside
         if (!fEditor->LoadFile(fname)) {
            new TGMsgBox(gClient->GetRoot(), fTab, "TGRootIDE",
                         Form("Error opening file \"%s\"", fname),
                         kMBIconExclamation, kMBOk);
         } else {
            fEditor->Layout();
            fTabEl->SetText(new TGString(gSystem->BaseName(fname)));
            fTab->MapSubwindows();
            fTab->Layout();
         }
      }
   }
   fEditor->Layout();
   fModified = kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Save current document.

Bool_t TGDocument::Save(const char *fname)
{
   TString sname;
   if (fname && strlen(fname) > 3) {
      sname = fname;
   }
   else sname = GetName();
   sname.Remove(TString::kLeading, '*');
   if (!fEditor->SaveFile(sname.Data())) {
      new TGMsgBox(gClient->GetRoot(), fTab, "TGRootIDE",
                   Form("Error saving file \"%s\"", sname.Data()),
                   kMBIconExclamation, kMBOk);
      return kFALSE;
   }
   fTabEl->SetText(new TGString(gSystem->BaseName(sname.Data())));
   fTabEl->Layout();
   fModified = kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle drop event.

void TGDocument::DataDropped(char *fname)
{
   if (strstr(GetName(),"Untitled")) {
      if (!fModified) {
         fTabEl->SetText(new TGString(Form("*%s", fTabEl->GetString())));
         fTabEl->Layout();
      }
      fModified = kTRUE;
   }
   SetName(fname);
   SetTitle(gSystem->BaseName(fname));
   fTabEl->SetText(new TGString(gSystem->BaseName(fname)));
   fTab->MapSubwindows();
   fTab->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if current document has been modified, and add or remove
/// a mark in front of tab name to indicate status.

void TGDocument::DataChanged()
{
   TList *hist = fEditor->GetHistory();
   if (hist->GetSize()) {
      if (!fModified) {
         fTabEl->SetText(new TGString(Form("*%s", fTabEl->GetString())));
         fTab->Layout();
      }
      fModified = kTRUE;
   }
   else {
      if (fModified) {
         TString sname(fTabEl->GetString());
         sname.Remove(TString::kLeading, '*');
         fTabEl->SetText(new TGString(Form("%s", sname.Data())));
         fTab->Layout();
      }
      fModified = kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TGRootIDE constructor with file name as first argument.

TGRootIDE::TGRootIDE(const char *filename, const TGWindow *p, UInt_t w,
                           UInt_t h) : TGMainFrame(p, w, h)
{
   Build();
   if (filename) {
      fTab->SetTab(fTab->GetNumberOfTabs()-1, kFALSE);
      fTab->MapSubwindows();
      fTab->Layout();
      LoadFile((char *)filename);
   }
   if (w > 0 && h > 0) {
      Resize(w, h > 500 ? h : 500);
      Layout();
   }
   MapWindow();
   fContents->DisplayDirectory();
   fContents->AddFile("..");        // up level directory
   fContents->SetViewMode(kLVDetails);
   fContents->Sort(kSortByType);
}

////////////////////////////////////////////////////////////////////////////////
/// TGRootIDE constructor with pointer to a TMacro as first argument.

TGRootIDE::TGRootIDE(TMacro *macro, const TGWindow *p, UInt_t w, UInt_t h) :
              TGMainFrame(p, w, h)
{
   Build();
   if (macro) {
      fTab->SetTab(fTab->GetNumberOfTabs()-1, kFALSE);
      fTab->MapSubwindows();
      fTab->Layout();
      fMacro = macro;
      TIter next(macro->GetListOfLines());
      TObjString *obj;
      while ((obj = (TObjString*) next())) {
         fTextEdit->AddLine(obj->GetName());
      }
      fStatusBar->SetText(Form("TMacro : %s: %ld lines read.",
                          macro->GetName(), fTextEdit->ReturnLineCount()), 0);
      fFilename = macro->GetName();
      fFilename += ".C";
      SetWindowName(Form("TMacro : %s - TGRootIDE", macro->GetName()));
   }
   if (w > 0 && h > 0) {
      Resize(w, h > 500 ? h : 500);
      Layout();
   }
   MapWindow();
   fContents->DisplayDirectory();
   fContents->AddFile("..");        // up level directory
   fContents->SetViewMode(kLVDetails);
   fContents->Sort(kSortByType);
}

////////////////////////////////////////////////////////////////////////////////
/// TGRootIDE destructor.

TGRootIDE::~TGRootIDE()
{
   fDocList->Delete();
   delete fDocList;
   delete fTimer;
}

////////////////////////////////////////////////////////////////////////////////
/// Build TGRootIDE widget.

void TGRootIDE::Build()
{
   fDocList = new TObjArray(100);
   fCurrent = 0;
   fNbDoc   = 0;
   fCurrentDoc = 0;
   fPid = gSystem->GetPid(); //TProcessID::GetSessionProcessID();
   fHtml = new THtml();

   SetCleanup(kDeepCleanup);
   fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 1, 1);
   fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);

   fMenuFile = new TGPopupMenu(fClient->GetRoot());
   fMenuFile->AddEntry(" &New", kM_FILE_NEW, 0,
                       gClient->GetPicture("ed_new.png"));
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry(" &Open...", kM_FILE_OPEN, 0,
                       gClient->GetPicture("ed_open.png"));
   fMenuFile->AddEntry(" &Close", kM_FILE_CLOSE);
   fMenuFile->AddEntry(" &Save", kM_FILE_SAVE, 0,
                       gClient->GetPicture("ed_save.png"));
   fMenuFile->AddEntry(" Save &As...", kM_FILE_SAVEAS, 0,
                       gClient->GetPicture("ed_saveas.png"));
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry(" &Print...", kM_FILE_PRINT, 0,
                       gClient->GetPicture("ed_print.png"));
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry(" E&xit", kM_FILE_EXIT, 0,
                       gClient->GetPicture("bld_exit.png"));

   fMenuEdit = new TGPopupMenu(fClient->GetRoot());
   fMenuEdit->AddEntry(" Cu&t             Ctrl+X", kM_EDIT_CUT, 0,
                       gClient->GetPicture("ed_cut.png"));
   fMenuEdit->AddEntry(" &Copy          Ctrl+C", kM_EDIT_COPY, 0,
                       gClient->GetPicture("ed_copy.png"));
   fMenuEdit->AddEntry(" &Paste         Ctrl+V", kM_EDIT_PASTE, 0,
                       gClient->GetPicture("ed_paste.png"));
   fMenuEdit->AddEntry(" De&lete        Del", kM_EDIT_DELETE, 0,
                       gClient->GetPicture("ed_delete.png"));
   fMenuEdit->AddSeparator();
   fMenuEdit->AddEntry(" Select &All   Ctrl+A", kM_EDIT_SELECTALL);
   fMenuEdit->AddSeparator();
   fMenuEdit->AddEntry(" Set &Font", kM_EDIT_SELFONT);

   fMenuTools = new TGPopupMenu(fClient->GetRoot());
   fMenuTools->AddEntry(" &Compile Macro  Ctrl+F7", kM_TOOLS_COMPILE, 0,
                       gClient->GetPicture("ed_compile.png"));
   fMenuTools->AddEntry(" &Execute Macro   Ctrl+F5", kM_TOOLS_EXECUTE, 0,
                       gClient->GetPicture("ed_execute.png"));
   fMenuTools->AddEntry(" &Interrupt              Shift+F5", kM_TOOLS_INTERRUPT, 0,
                       gClient->GetPicture("ed_interrupt.png"));
   fMenuTools->AddSeparator();
   fMenuTools->AddEntry(" Start &Browser", kM_TOOLS_BROWSER);
   fMenuTools->AddSeparator();
   fMenuTools->AddEntry(" Cleanup &Log Files", kM_TOOLS_CLEAN_LOG);

   fMenuEdit->DisableEntry(kM_EDIT_CUT);
   fMenuEdit->DisableEntry(kM_EDIT_COPY);
   fMenuEdit->DisableEntry(kM_EDIT_DELETE);
   fMenuEdit->DisableEntry(kM_EDIT_PASTE);

   fMenuSearch = new TGPopupMenu(fClient->GetRoot());
   fMenuSearch->AddEntry(" &Find...         Ctrl+F", kM_SEARCH_FIND, 0,
                       gClient->GetPicture("ed_find.png"));
   fMenuSearch->AddEntry(" Find &Next    F3", kM_SEARCH_FINDNEXT, 0,
                       gClient->GetPicture("ed_findnext.png"));
   fMenuSearch->AddSeparator();
   fMenuSearch->AddEntry(" &Goto Line... Ctrl+L", kM_SEARCH_GOTO, 0,
                       gClient->GetPicture("ed_goto.png"));

   fMenuHelp = new TGPopupMenu(fClient->GetRoot());
   fMenuHelp->AddEntry(" &Help Topics    F1", kM_HELP_CONTENTS, 0,
                       gClient->GetPicture("ed_help.png"));
   fMenuHelp->AddSeparator();
   fMenuHelp->AddEntry(" &About...", kM_HELP_ABOUT, 0,
                       gClient->GetPicture("about.xpm"));

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
   fCommandBuf = new TGTextBuffer(256);
   fComboCmd   = new TGComboBox(fToolBar, "", 1);
   fCommand    = fComboCmd->GetTextEntry();
   fCommandBuf = fCommand->GetBuffer();
   fCommand->Associate(this);
   fComboCmd->Resize(200, fCommand->GetDefaultHeight());
   fToolBar->AddFrame(fComboCmd, new TGLayoutHints(kLHintsCenterY |
            kLHintsRight | kLHintsExpandX, 5, 5, 1, 1));

   TString defhist(Form("%s/.root_hist", gSystem->UnixPathName(gSystem->HomeDirectory())));
   FILE *lunin = fopen(defhist.Data(), "rt");
   if (lunin) {
      char histline[256];
      while (fgets(histline, 256, lunin)) {
         histline[strlen(histline)-1] = 0; // remove trailing "\n"
         fComboCmd->InsertEntry(histline, 0, -1);
      }
      fclose(lunin);
   }

   fToolBar->AddFrame(fLabel = new TGLabel(fToolBar, "Command (local):"),
            new TGLayoutHints(kLHintsCenterY | kLHintsRight, 5, 5, 1, 1));
   AddFrame(fToolBar, new TGLayoutHints(kLHintsTop | kLHintsExpandX,
            0, 0, 0, 0));
   AddFrame(new TGHorizontal3DLine(this),
            new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0,0,2,2));

   fToolBar->GetButton(kM_EDIT_CUT)->SetState(kButtonDisabled);
   fToolBar->GetButton(kM_EDIT_COPY)->SetState(kButtonDisabled);
   fToolBar->GetButton(kM_EDIT_DELETE)->SetState(kButtonDisabled);
   fToolBar->GetButton(kM_EDIT_PASTE)->SetState(kButtonDisabled);

   TGHorizontalFrame *hf = new TGHorizontalFrame(this, 100, 100);

   TGVerticalFrame *vf3 = new TGVerticalFrame(hf, 160, 100, kSunkenFrame | kFixedWidth);

   fDirBuf   = new TGTextBuffer(256);
   fDirCombo = new TGComboBox(vf3, "");
   fDir      = fDirCombo->GetTextEntry();
   fDirBuf   = fDir->GetBuffer();
   fDirCombo->Resize(200, fDir->GetDefaultHeight());
   fDir->Connect("ReturnPressed()", "TGRootIDE", this, "DirChanged()");

   fDirCombo->AddEntry(gSystem->WorkingDirectory(), 1);
   gSystem->ChangeDirectory(gSystem->WorkingDirectory());
   fDir->SetText(gSystem->WorkingDirectory());

   fDirCombo->Select(0);
   fDirCombo->Connect("Selected(char *)", "TGRootIDE", this, "DirSelected(char *)");

   vf3->AddFrame(fDirCombo, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,2,2));

   TGListView* lv = new TGListView(vf3, 200, 100);
   vf3->AddFrame(lv,new TGLayoutHints(kLHintsExpandY | kLHintsExpandX));
   Pixel_t white;
   gClient->GetColorByName("white",white);
   fContents = new TGFileContainer(lv, kSunkenFrame, white);
   fContents->Associate(this);

   fFileType = new TGComboBox(vf3, " All Files (*.*)");
   Int_t dropt = 1;
   fFileType->AddEntry(" All Files (*.*)", dropt++);
   fFileType->AddEntry(" C/C++ Files (*.c;*.cxx;*.h;...)", dropt++);
   fFileType->AddEntry(" Text Files (*.txt)", dropt++);
   fFilter = fFileType->GetTextEntry();
   fFileType->Resize(200, 20);
   vf3->AddFrame(fFileType, new TGLayoutHints(kLHintsExpandX, 1, 1, 1, 1));
   fFileType->Connect("Selected(Int_t)", "TGRootIDE", this, "ApplyFilter(Int_t)");

   hf->AddFrame(vf3, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));

   TGVSplitter *splitter = new TGVSplitter(hf, 4);
   splitter->SetFrame(vf3, kTRUE);
   hf->AddFrame(splitter, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));

   TGVerticalFrame *vf2 = new TGVerticalFrame(hf, 100, 100, kSunkenFrame);

   fTab = new TGTab(vf2, 300, 300);
   TGCompositeFrame *tf = fTab->AddTab("HTML");
   tf->SetLayoutManager(new TGHorizontalLayout(tf));

   // vertical frame
   fVerticalFrame = new TGVerticalFrame(tf,727,600,kVerticalFrame);

   fHorizontalFrame = new TGHorizontalFrame(fVerticalFrame,727,600);

   fBack = new TGPictureButton(fHorizontalFrame,
            gClient->GetPicture("GoBack.gif"));
   fBack->SetToolTipText("Go Back");
   fHorizontalFrame->AddFrame(fBack, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));
   fBack->Connect("Clicked()", "TGRootIDE", this, "Back()");

   fForward = new TGPictureButton(fHorizontalFrame,
            gClient->GetPicture("GoForward.gif"));
   fForward->SetToolTipText("Go Forward");
   fHorizontalFrame->AddFrame(fForward, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));
   fForward->Connect("Clicked()", "TGRootIDE", this, "Forward()");

   fReload = new TGPictureButton(fHorizontalFrame,
            gClient->GetPicture("ReloadPage.gif"));
   fReload->SetToolTipText("Reload Page");
   fHorizontalFrame->AddFrame(fReload, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));
   fReload->Connect("Clicked()", "TGRootIDE", this, "Reload()");

   fStop = new TGPictureButton(fHorizontalFrame,
            gClient->GetPicture("StopLoading.gif"));
   fStop->SetToolTipText("Stop Loading");
   fHorizontalFrame->AddFrame(fStop, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));
   fStop->Connect("Clicked()", "TGRootIDE", this, "Stop()");

   fHome = new TGPictureButton(fHorizontalFrame,
           gClient->GetPicture("GoHome.gif"));
   fHome->SetToolTipText("Go to ROOT HomePage\n  (http://root.cern)");
   fHorizontalFrame->AddFrame(fHome, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY,2,2,2,2));
   fHome->Connect("Clicked()", "TGRootIDE", this, "Selected(=\"http://root.cern\")");

   // combo box
   fURLBuf   = new TGTextBuffer(256);
   fComboBox = new TGComboBox(fHorizontalFrame, "");
   fURL      = fComboBox->GetTextEntry();
   fURLBuf   = fURL->GetBuffer();
   fComboBox->Resize(200, fURL->GetDefaultHeight());
   fURL->Connect("ReturnPressed()", "TGRootIDE", this, "URLChanged()");

   fComboBox->AddEntry("http://root.cern", 1);
   fComboBox->AddEntry("http://root.cern/root/htmldoc/ClassIndex.html", 2);
   fURL->SetText("http://root.cern/root/htmldoc/ClassIndex.html");

   fComboBox->Select(0);
   fComboBox->Connect("Selected(char *)", "TGRootIDE", this, "Selected(char *)");

   fHorizontalFrame->AddFrame(fComboBox, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

   fVerticalFrame->AddFrame(fHorizontalFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,2,2));

   // embedded canvas
   fGuiHtml = new TGHtml(fVerticalFrame, 10, 10, -1);
   fVerticalFrame->AddFrame(fGuiHtml, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));

   tf->AddFrame(fVerticalFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));

   fGuiHtml->Connect("MouseOver(char *)", "TGRootIDE", this, "MouseOver(char *)");
   fGuiHtml->Connect("MouseDown(char *)", "TGRootIDE", this, "MouseDown(char *)");
   Selected("http://root.cern/root/htmldoc/ClassIndex.html");
   fGuiHtml->Layout();

   tf = fTab->AddTab("Untitled");
   tf->SetLayoutManager(new TGHorizontalLayout(tf));
   fTextEdit = new TGTextEdit(tf, 10, 10, 1);
   fTextEdit->Associate(this);
   TGTabElement *tabel = fTab->GetTabTab(1);
   tabel->ShowClose();

   // set selected text colors
   Pixel_t pxl;
   gClient->GetColorByName("#ccccff", pxl);
   fTextEdit->SetSelectBack(pxl);
   fTextEdit->SetSelectFore(TGFrame::GetBlackPixel());
   fCurrentDoc = new TGDocument("", "", 1, fTab,
                  fTab->GetTabTab(1),
                  fTextEdit, fDocList);
   fDocList->Add((TObject *)fCurrentDoc);
   tf->AddFrame(fTextEdit, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   vf2->AddFrame(fTab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   fTab->Connect("Selected(Int_t)", "TGRootIDE", this, "DoTab(Int_t)");
   fTab->Connect("CloseTab(Int_t)", "TGRootIDE", this, "CloseTab(Int_t)");

   hf->AddFrame(vf2, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   AddFrame(hf, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   gClient->GetColorByName("#ccccff", pxl);
   fTextView = new TGTextView(this, 10, 100, 1);
   fTextView->SetSelectBack(pxl);
   fTextView->SetSelectFore(TGFrame::GetBlackPixel());
   fTextView->Associate(this);
   fTextView->ChangeOptions(fTextView->GetOptions() | kFixedHeight);

   TGHSplitter *hsplitter = new TGHSplitter(this, 4);
   hsplitter->SetFrame(fTextView, kFALSE);
   AddFrame(hsplitter, new TGLayoutHints(kLHintsLeft | kLHintsExpandX));

   AddFrame(fTextView, new TGLayoutHints(kLHintsExpandX));

   Int_t parts[] = { 75, 25 };
   fStatusBar = new TGStatusBar(this);
   fStatusBar->SetCleanup(kDeepCleanup);
   fStatusBar->SetParts(parts, 2);
   AddFrame(fStatusBar, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 0, 0, 3, 0));

   SetClassHints("TGRootIDE", "TGRootIDE");
   SetWindowName("Untitled - TGRootIDE");

   fMacro = 0;
   fFilename = "Untitled";
   fStatusBar->SetText(fFilename.Data(), 0);

   MapSubwindows();
   Resize(GetDefaultWidth() + 50, GetDefaultHeight() > 500 ? GetDefaultHeight() : 500);
   Layout();

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
/// Load a file into the editor. If fname is 0, a TGFileDialog will popup.

void TGRootIDE::LoadFile(char *fname)
{
   TGFileInfo fi;
   fi.fFileTypes = ed_filetypes;
   if (fname == 0) {
      new TGFileDialog(fClient->GetDefaultRoot(), this, kFDOpen, &fi);
      if (fi.fFilename && strlen(fi.fFilename)) {
         fname = fi.fFilename;
      }
   }
   if (fname) {
      const char *p = fTab->GetTabTab(fTab->GetCurrent())->GetString();
      if (!strcmp(p, "HTML")) {
         TString filename(fname);
         if (filename.EndsWith(".htm") ||
             filename.EndsWith(".html")) {
            Selected(Form("file://%s", filename.Data()));
         }
         else {
            TString pathtmp = Form("%s/%s.html",
               gSystem->UnixPathName(gSystem->TempDirectory()),
               gSystem->BaseName(fname));
            fHtml->Convert(fname, fname,
               gSystem->UnixPathName(gSystem->TempDirectory()),
               gSystem->UnixPathName(gSystem->TempDirectory()));
            Selected(Form("file://%s", pathtmp.Data()));
         }
         //gSystem->Unlink(pathtmp.Data());
      }
      else {
         TGDocument *doc = new TGDocument(fname, gSystem->BaseName(fname),
                                        fTab->GetNumberOfTabs()+1, fTab,
                                        0, 0, fDocList);
         fDocList->Add((TObject *)doc);
         fCurrent = fTab->GetCurrent();
         fCurrentDoc = doc;
         fFilename = fCurrentDoc->GetName();
         fTextEdit = fCurrentDoc->GetTextEdit();
         fTextEdit->SetFocus();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save the edited text in the file "fname".

void TGRootIDE::SaveFile(const char *fname)
{
   char *p;
   if (!fCurrentDoc) return;
   fCurrentDoc->Save(fname);
   if ((p = (char *)strrchr(fname, '/')) == 0) {
      p = (char *)fname;
   } else {
      ++p;
   }
   fStatusBar->SetText(Form("%s: %ld lines written.", p,
                       fTextEdit->ReturnLineCount()), 0);
   SetWindowName(Form("%s - TGRootIDE", p));
   fTextChanged = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Save the edited text in a file selected with TGFileDialog.
/// Shouldn't we create a backup file?

Bool_t TGRootIDE::SaveFileAs()
{
   if (!fCurrentDoc) return kFALSE;
   static TString dir(".");
   static Bool_t overwr = kFALSE;
   TGFileInfo fi;
   fi.fFileTypes = ed_filetypes;
   fi.SetIniDir(dir);
   fi.fOverwrite = overwr;
   new TGFileDialog(fClient->GetDefaultRoot(), this, kFDSave, &fi);
   overwr = fi.fOverwrite;
   if (fi.fFilename && strlen(fi.fFilename)) {
      fCurrentDoc->Save(fi.fFilename);
      fFilename = fi.fFilename;
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if file has to be saved in case of modifications.

Int_t TGRootIDE::IsSaved()
{
   Int_t ret;
   TGDocument *doc = 0;
   TIter next(fDocList);
   while ((doc = (TGDocument *)next())) {
      if (doc->IsModified()) {
         TString sfname(doc->GetName());
         new TGMsgBox(fClient->GetRoot(), this, "TGRootIDE",
                      Form("%s has been modified. Do you want to save the changes?",
                           sfname.Data()),
                      kMBIconExclamation, kMBYes | kMBNo | kMBCancel, &ret);
         if (ret == kMBYes) {
            doc->Save();
         }
      }
   }
   return kMBNo; //ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Open the print dialog and send current buffer to printer.

void TGRootIDE::PrintText()
{
   Int_t ret = 0;
   if (!gEPrinter) {
      gEPrinter = StrDup("892_2_cor"); // use gEnv
      gEPrintCommand = StrDup("xprint");
   }
   new TGPrintDialog(fClient->GetDefaultRoot(), this, 400, 150,
                     &gEPrinter, &gEPrintCommand, &ret);
   if (ret) {
      fTextEdit->Print();
      fStatusBar->SetText(Form("Printed: %s", fFilename.Data()), 0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Close TGRootIDE window.

void TGRootIDE::CloseWindow()
{
   if (fExiting) {
      return;
   }
   fExiting = kTRUE;
   if (IsSaved() == kMBCancel) {
      fExiting = kFALSE;
      return;
   }
   fExiting = kFALSE;
   Cleanup();
#ifdef WIN32
   gSystem->Exec(Form("del %s\\*.html", gSystem->TempDirectory()));
   gSystem->Exec(Form("del %s\\*.C", gSystem->TempDirectory()));
#else
   gSystem->Exec(Form("rm -f %s/*.html", gSystem->TempDirectory()));
   gSystem->Exec(Form("rm -f %s/*.C", gSystem->TempDirectory()));
#endif
   delete this;
   gApplication->Terminate(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Keyboard event handler.

Bool_t TGRootIDE::HandleKey(Event_t *event)
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
            case kKey_F4:
            case kKey_W:
               SendMessage(this, MK_MSG(kC_COMMAND, kCM_MENU),
                           kM_FILE_CLOSE, 0);
               return kTRUE;
            case kKey_F5:
               ExecuteMacro();
               return kTRUE;
            case kKey_F7:
               CompileMacro();
               return kTRUE;
            case kKey_Tab:
               if (fTab->GetCurrent() == fTab->GetNumberOfTabs()-1)
                  fTab->SetTab(0);
               else
                  fTab->SetTab(fTab->GetCurrent()+1);
               break;
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

void TGRootIDE::ClearText()
{
   fTextEdit->Clear();
   fMacro = 0;
   fFilename = "Untitled";
   SetWindowName("Untitled - TGRootIDE");
   fStatusBar->SetText("New File", 0);
   fTextChanged = kFALSE;
   fTab->GetTabTab(fTab->GetCurrent())->SetText(new TGString("Untitled"));
   fTab->MapSubwindows();
   fTab->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Invokes search dialog, or just search previous string if again is true.

void TGRootIDE::Search(Bool_t again)
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

void TGRootIDE::Goto()
{
   Long_t ret;

   new TGGotoDialog(fClient->GetDefaultRoot(), this, 400, 150, &ret);

   if (ret >= 0)
      fTextEdit->Goto(ret-1);
}

////////////////////////////////////////////////////////////////////////////////
/// Save the edited text in a temporary macro, then compile it.

void TGRootIDE::CompileMacro()
{
   if (fTextEdit->ReturnLineCount() < 3)
      return;
   if ((fMacro) || (fFilename.Contains("Untitled"))) {
      if (!SaveFileAs())
         return;
   }
   TString pathtmp = Form("%s/ride.%d.log", gSystem->TempDirectory(), fPid);
   gSystem->RedirectOutput(pathtmp.Data(), "a");

   TString tmpfile = gSystem->BaseName(fFilename.Data());
   gSystem->PrependPathName(gSystem->TempDirectory(), tmpfile);
   fTextEdit->SaveFile(tmpfile, kFALSE);
   gSystem->CompileMacro(tmpfile);
   gSystem->Unlink(tmpfile);

   gSystem->RedirectOutput(0);
   fTextView->LoadFile(pathtmp.Data());
   if (fTextView->ReturnLineCount() > 7)
      fTextView->SetVsbPosition(fTextView->ReturnLineCount());
}

////////////////////////////////////////////////////////////////////////////////
/// Save the edited text in a temporary macro, execute it, and then delete
/// the temporary file.

void TGRootIDE::ExecuteMacro()
{
   if (fTextEdit->ReturnLineCount() < 3)
      return;
   if (fMacro) {
      fMacro->Exec();
      return;
   }
   if (fTextChanged) {
      Int_t ret;
      new TGMsgBox(fClient->GetRoot(), this, "TGRootIDE",
            "The text has been modified. Do you want to save the changes?",
            kMBIconExclamation, kMBYes | kMBNo | kMBCancel, &ret);
      if (ret == kMBYes) {
         if (fFilename.Contains("Untitled"))
            SaveFileAs();
         else
            SaveFile(fFilename.Data());
         fTextChanged = kFALSE;
      }
      if (ret == kMBCancel)
         return;
   }
   if (fFilename.Contains("Untitled")) {
      if (!SaveFileAs())
         return;
   }
   TString pathtmp = Form("%s/ride.%d.log", gSystem->TempDirectory(), fPid);
   gSystem->RedirectOutput(pathtmp.Data(), "a");

   TString tmpfile = gSystem->BaseName(fFilename.Data());
   gSystem->PrependPathName(gSystem->TempDirectory(), tmpfile);
   gROOT->SetExecutingMacro(kTRUE);
   fTextEdit->SaveFile(tmpfile, kFALSE);
   gROOT->Macro(tmpfile);
   gSystem->Unlink(tmpfile);
   gROOT->SetExecutingMacro(kFALSE);
   gSystem->RedirectOutput(0);
   fTextView->LoadFile(pathtmp.Data());
   if (fTextView->ReturnLineCount() > 7)
      fTextView->SetVsbPosition(fTextView->ReturnLineCount());
}

////////////////////////////////////////////////////////////////////////////////
/// Interrupt execution of a macro.

void TGRootIDE::InterruptMacro()
{
   gROOT->SetInterrupt(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Display ROOT splash screen.

void TGRootIDE::About()
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

Bool_t TGRootIDE::HandleTimer(TTimer *t)
{
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
   if (fTextEdit && fTextEdit->IsMarked()) {
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
   if (fTextEdit) {
      TGLongPosition pos = fTextEdit->GetCurrentPos();
      fStatusBar->SetText(Form("Ln %ld, Ch %ld", pos.fY, pos.fX), 1);
   }
   fTimer->Reset();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle menu and other command generated by the user.

Bool_t TGRootIDE::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   TRootHelpDialog *hd;

   switch(GET_MSG(msg)) {
      case kC_CONTAINER:
         switch (GET_SUBMSG(msg)) {
            case kCT_ITEMDBLCLICK:
               if (parm1==kButton1) OnDoubleClick((TGLVEntry*)fContents->GetLastActive(), parm1);
               break;
         }
         break;
      case kC_COMMAND:
         switch(GET_SUBMSG(msg)) {
            case kCM_BUTTON:
            case kCM_MENU:
               switch (parm1) {
                  // "File" menu related events
                  case kM_FILE_NEW:
                     {
                        fDocList->Add((TObject *)new TGDocument(0, 0,
                           fTab->GetNumberOfTabs(), fTab, 0, 0, fDocList));
                     }
                     break;
                  case kM_FILE_CLOSE:
                     if (fCurrentDoc) {
                        fCurrentDoc->Close();
                     }
                     fTab->Layout();
                     break;
                  case kM_FILE_OPEN:
                     LoadFile();
                     break;
                  case kM_FILE_SAVE:
                     if (fFilename.Contains("Untitled"))
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
                        char fontname[256];
                        TGFontDialog::FontProp_t prop;
                        new TGFontDialog(fClient->GetRoot(), this, &prop);
                        if (prop.fName != "") {
                           sprintf(fontname,"-*-%s-%s-%c-*-*-%d-*-*-*-*-*-*-*",
                                   prop.fName.Data(), prop.fBold ? "bold" : "medium",
                                   prop.fItalic ? 'i' : 'r',
                                   prop.fSize);
                           if (!gVirtualX->ListFonts(fontname, 10, count)) {
                              sprintf(fontname,"-*-%s-%s-%c-*-*-%d-*-*-*-*-*-*-*",
                                      prop.fName.Data(), prop.fBold ? "bold" : "medium",
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
                  case kM_TOOLS_BROWSER:
                     new TBrowser();
                     break;
                  case kM_TOOLS_CLEAN_LOG:
#ifdef WIN32
                     gSystem->Exec(Form("del %s\\ride.*.log", gSystem->TempDirectory()));
#else
                     gSystem->Exec(Form("rm -f %s/ride.*.log", gSystem->TempDirectory()));
#endif
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
               if (parm1 == 1) {
                  // here copy the string from text buffer to return variable
                  const char *string = fCommandBuf->GetString();
                  if (strlen(string) > 1) {
                     // form temporary file path
                     TString pathtmp = Form("%s/ride.%d.log", gSystem->TempDirectory(),
                                             fPid);
                     TString sPrompt = ((TRint*)gROOT->GetApplication())->GetPrompt();
                     FILE *lunout = fopen(pathtmp.Data(), "a+t");
                     if (lunout) {
                        fputs(Form("%s%s\n",sPrompt.Data(), string), lunout);
                        fclose(lunout);
                     }
                     gSystem->RedirectOutput(pathtmp.Data(), "a");
                     gApplication->SetBit(TApplication::kProcessRemotely);
                     gROOT->ProcessLine(string);
                     //fComboCmd->ReturnPressed();
                     fComboCmd->InsertEntry(string, 0, -1);
                     Gl_histadd((char *)string);
                     gSystem->RedirectOutput(0);
                     fTextView->LoadFile(pathtmp.Data());
                     if (fTextView->ReturnLineCount() > 7)
                        fTextView->SetVsbPosition(fTextView->ReturnLineCount());
                     CheckRemote(string);
                     fCommand->Clear();
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

////////////////////////////////////////////////////////////////////////////////
/// Display content of ROOT file.

void TGRootIDE::DisplayFile(const TString &fname)
{
   TFile file(fname);
   fContents->RemoveAll();
   fContents->AddFile(gSystem->WorkingDirectory());
   fContents->SetPagePosition(0,0);
   fContents->SetColHeaders("Name","Title");

   TIter next(file.GetListOfKeys());
   TKey *key;

   while ((key=(TKey*)next())) {
      TString cname = key->GetClassName();
      TString name = key->GetName();
      TGLVEntry *entry = new TGLVEntry(fContents,name,cname);
      entry->SetSubnames(key->GetTitle());
      fContents->AddItem(entry);

      // user data is a filename
      entry->SetUserData((void*)StrDup(fname));
   }
   fContents->Sort(kSortByType);
   Resize();
}

////////////////////////////////////////////////////////////////////////////////
/// Display content of directory.

void TGRootIDE::DisplayDirectory(const TString &fname)
{
   fContents->SetDefaultHeaders();
   gSystem->ChangeDirectory(fname);
   fContents->ChangeDirectory(fname);
   fContents->DisplayDirectory();
   fContents->AddFile("..");  // up level directory
   fContents->Sort(kSortByType);
   Resize();
   fDir->SetText(gSystem->WorkingDirectory());
   if (!fDirCombo->FindEntry(gSystem->WorkingDirectory()))
      fDirCombo->AddEntry(gSystem->WorkingDirectory(),
                          fDirCombo->GetNumberOfEntries()+1);
}

////////////////////////////////////////////////////////////////////////////////
/// Display object located in file.

void TGRootIDE::DisplayObject(const TString& fname, const TString& name)
{
   TDirectory *sav = gDirectory;

   static TFile *file = 0;
   if (file) delete file;     // close
   file = new TFile(fname);   // reopen

   TObject* obj = file->Get(name);
   if (obj) {
      if (!obj->IsFolder()) {
         obj->Browse(0);
      } else obj->Print();
   }
   gDirectory = sav;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if given a text file
/// Uses the specification given on p86 of the Camel book
/// - Text files have no NULLs in the first block
/// - and less than 30% of characters with high bit set

static Bool_t IsTextFile(const char *candidate)
{
   Int_t i;
   Int_t nchars;
   Int_t weirdcount = 0;
   char buffer[512];
   FILE *infile;
   FileStat_t buf;

   gSystem->GetPathInfo(candidate, buf);
   if (!(buf.fMode & kS_IFREG))
      return kFALSE;

   infile = fopen(candidate, "r");
   if (infile) {
      // Read a block
      nchars = fread(buffer, 1, 512, infile);
      fclose (infile);
      // Examine the block
      for (i = 0; i < nchars; i++) {
         if (buffer[i] & 128)
            weirdcount++;
         if (buffer[i] == '\0')
            // No NULLs in text files
            return kFALSE;
      }
      if ((nchars > 0) && ((weirdcount * 100 / nchars) > 30))
         return kFALSE;
   } else {
      // Couldn't open it. Not a text file then
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle double click in TGListView.

void TGRootIDE::OnDoubleClick(TGLVEntry* f, Int_t btn)
{
   if (btn!=kButton1) return;
   gVirtualX->SetCursor(fContents->GetId(),gVirtualX->CreateCursor(kWatch));

   TString name(f->GetTitle());
   const char* fname = (const char*)f->GetUserData();

   if (IsTextFile(name.Data())) {
      LoadFile((char *)name.Data());
   }
   else if (fname) {
      DisplayObject(fname,name);
   } else if (name.EndsWith(".root")) {
      DisplayFile(name);
   } else {
      DisplayDirectory(name);
   }
   fContents->Sort(kSortByType);
   gVirtualX->SetCursor(fContents->GetId(),gVirtualX->CreateCursor(kPointer));
}

////////////////////////////////////////////////////////////////////////////////
/// Handle Tab navigation.

void TGRootIDE::DoTab(Int_t id)
{
   fCurrentDoc = 0;
   fFilename = "";
   //fTextEdit = 0;
   TGDocument *doc = 0;
   TIter next(fDocList);
   while ((doc = (TGDocument *) next())) {
      if (doc->GetTabEl() == fTab->GetTabTab(id)) {
         fCurrentDoc = doc;
         fFilename = fCurrentDoc->GetName();
         fTextEdit = fCurrentDoc->GetTextEdit();
         break;
      }
   }
   const char *p = fTab->GetTabTab(id)->GetString();
   fFilename = p;
   SetWindowName(Form("%s - TGRootIDE", p));
   fCurrent = id;
}

////////////////////////////////////////////////////////////////////////////////
/// Close tab "id".

void TGRootIDE::CloseTab(Int_t id)
{
   if (fCurrentDoc) {
      fCurrentDoc->Close();
   }
   else if (fTab->GetNumberOfTabs() > 2) {
      fTab->RemoveTab(id);
   }
   fTab->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Apply filter selected in combo box to the file list view.

void TGRootIDE::ApplyFilter(Int_t id)
{
   fContents->SetFilter(filters[id]);
   fContents->DisplayDirectory();
   fContents->AddFile("..");        // up level directory
   fContents->Sort(kSortByType);
}

////////////////////////////////////////////////////////////////////////////////
/// A directory has been selected in the navigation history.

void TGRootIDE::DirSelected(const char *uri)
{
   fDir->SetText(uri);
   if (!fDirCombo->FindEntry(uri))
      fDirCombo->AddEntry(uri, fDirCombo->GetNumberOfEntries()+1);
   DisplayDirectory(uri);
}

////////////////////////////////////////////////////////////////////////////////
/// A directory has been typed in the text entry of the navigation history.

void TGRootIDE::DirChanged()
{
   const char *string = fDir->GetText();
   if (string) {
      TString buf = string;
      DirSelected(buf.Data());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Temporary function to read remote pictures

static char *ReadRemote(const char *url)
{
   static char *buf = 0;
   TUrl fUrl(url);

   TString msg = "GET ";
   msg += fUrl.GetProtocol();
   msg += "://";
   msg += fUrl.GetHost();
   msg += ":";
   msg += fUrl.GetPort();
   msg += "/";
   msg += fUrl.GetFile();
   msg += "\r\n";

   TString uri(url);
   if (!uri.BeginsWith("http://"))
      return 0;
   TSocket s(fUrl.GetHost(), fUrl.GetPort());
   if (!s.IsValid())
      return 0;
   if (s.SendRaw(msg.Data(), msg.Length()) == -1)
      return 0;
   Int_t size = 1024*1024;
   buf = (char *)calloc(size, sizeof(char));
   if (s.RecvRaw(buf, size) == -1) {
      free(buf);
      return 0;
   }
   return buf;
}

////////////////////////////////////////////////////////////////////////////////
/// A URL has been selected, either by a click on a link or by the
/// navigation buttons, or by history combobox / text entry.

void TGRootIDE::Selected(const char *uri)
{
   char *buf = 0;
   FILE *f;

   gVirtualX->SetCursor(fGuiHtml->GetId(), gVirtualX->CreateCursor(kWatch));
   TString surl(gSystem->UnixPathName(uri));
   // if url does not contains "http://", prepend "file://" (local navigation)
   if (!surl.BeginsWith("http://") && !surl.BeginsWith("file://"))
      surl.Prepend("file://");
   if (surl.EndsWith(".root")) {
      // Open Root files directly and open a Root browser.
      TFile *f = TFile::Open(surl.Data());
      if (f && !f->IsZombie()) {
         f->Browse(new TBrowser());
      }
      gVirtualX->SetCursor(fGuiHtml->GetId(), gVirtualX->CreateCursor(kPointer));
      return;
   }
   TUrl url(surl.Data());
   if ((!strcmp(url.GetProtocol(), "http"))) {
      // web file...
      buf = ReadRemote(url.GetUrl());
      if (buf) {
         // display html page
         fGuiHtml->Clear();
         fGuiHtml->Layout();
         fGuiHtml->SetBaseUri(url.GetUrl());
         fGuiHtml->ParseText(buf);
         free(buf);
         fURL->SetText(surl.Data());
         if (!fComboBox->FindEntry(surl.Data()))
            fComboBox->AddEntry(surl.Data(), fComboBox->GetNumberOfEntries()+1);
      }
      else {
         // something went wrong --> display error
         fGuiHtml->Clear();
         fGuiHtml->Layout();
         fGuiHtml->SetBaseUri("");
         for (int i=0; HtmlError[i]; i++) {
            fGuiHtml->ParseText((char *)HtmlError[i]);
         }
      }
   }
   else {
      // local file...
      f = fopen(url.GetFile(), "r");
      if (f) {
         // file is opened (and valid)
         fGuiHtml->Clear();
         fGuiHtml->Layout();
         fGuiHtml->SetBaseUri("");
         buf = (char *)calloc(4096, sizeof(char));
         while (fgets(buf, 4096, f)) {
            fGuiHtml->ParseText(buf);
         }
         free(buf);
         fclose(f);
         fURL->SetText(surl.Data());
         if (!fComboBox->FindEntry(surl.Data()))
            fComboBox->AddEntry(surl.Data(), fComboBox->GetNumberOfEntries()+1);
      }
      else {
         // something went wrong --> display error
         fGuiHtml->Clear();
         fGuiHtml->Layout();
         fGuiHtml->SetBaseUri("");
         for (int i=0; HtmlError[i]; i++) {
            fGuiHtml->ParseText((char *)HtmlError[i]);
         }
      }
   }
   gVirtualX->SetCursor(fGuiHtml->GetId(), gVirtualX->CreateCursor(kPointer));
   fGuiHtml->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// The text entry of navigation history has changed.

void TGRootIDE::URLChanged()
{
   const char *string = fURL->GetText();
   if (string) {
      TString buf = gSystem->UnixPathName(string);
      Selected(buf.Data());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle "Back" navigation button.

void TGRootIDE::Back()
{
   Int_t index = 0;
   const char *string = fURL->GetText();
   TGLBEntry * lbe1 = fComboBox->FindEntry(string);
   if (lbe1)
      index = lbe1->EntryId();
   if (index > 0) {
      fComboBox->Select(index - 1, kTRUE);
      TGTextLBEntry *entry = (TGTextLBEntry *)fComboBox->GetSelectedEntry();
      if (entry) {
         const char *string = entry->GetTitle();
         if (string)
            Selected(string);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle "Forward" navigation button.

void TGRootIDE::Forward()
{
   Int_t index = 0;
   const char *string = fURL->GetText();
   TGLBEntry * lbe1 = fComboBox->FindEntry(string);
   if (lbe1)
      index = lbe1->EntryId();
   if (index < fComboBox->GetNumberOfEntries()) {
      fComboBox->Select(index + 1, kTRUE);
      TGTextLBEntry *entry = (TGTextLBEntry *)fComboBox->GetSelectedEntry();
      if (entry) {
         const char *string = entry->GetTitle();
         if (string)
            Selected(string);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle "Reload" navigation button.

void TGRootIDE::Reload()
{
   const char *string = fURL->GetText();
   if (string)
      Selected(string);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle "Stop Loading" navigation button.
/// Not active for the time being.

void TGRootIDE::Stop()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Handle MouseOver signal from TGHtml widget.

void TGRootIDE::MouseOver(char *url)
{
   fStatusBar->SetText(url, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle MouseDown signal from TGHtml widget.

void TGRootIDE::MouseDown(char *url)
{
   Selected(url);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if actual ROOT session is a remote one or a local one.

void TGRootIDE::CheckRemote(const char * /*str*/)
{
   Pixel_t pxl;
   TString sPrompt = ((TRint*)gROOT->GetApplication())->GetPrompt();
   Int_t end = sPrompt.Index(":root [", 0);
   if (end > 0 && end != kNPOS) {
      // remote session
      sPrompt.Remove(end);
      gClient->GetColorByName("#ff0000", pxl);
      fLabel->SetTextColor(pxl);
      fLabel->SetText(Form("Command (%s):", sPrompt.Data()));
   }
   else {
      // local session
      gClient->GetColorByName("#000000", pxl);
      fLabel->SetTextColor(pxl);
      fLabel->SetText("Command (local):");
   }
   fToolBar->Layout();
}


