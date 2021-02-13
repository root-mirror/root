// @(#)root/gui:$Id: f3cd439bd51d763ffd53693e89c42b2eaa27c520 $
// Author: Fons Rademakers   20/01/98

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGFileDialog                                                         //
//                                                                      //
// This class creates a file selection dialog. It contains a combo box  //
// to select the desired directory. A listview with the different       //
// files in the current directory and a combo box with which you can    //
// select a filter (on file extensions).                                //
// When creating a file dialog one passes a pointer to a TGFileInfo     //
// object. In this object you can set the fFileTypes and fIniDir to     //
// specify the list of file types for the filter combo box and the      //
// initial directory. When the TGFileDialog ctor returns the selected   //
// file name can be found in the TGFileInfo::fFilename field and the    //
// selected directory in TGFileInfo::fIniDir. The fFilename and         //
// fIniDir are deleted by the TGFileInfo dtor.                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGFileDialog.h"
#include "TGLabel.h"
#include "TGTextEntry.h"
#include "TGComboBox.h"
#include "TGListView.h"
#include "TGFSContainer.h"
#include "TGFSComboBox.h"
#include "TGMsgBox.h"
#include "TGInputDialog.h"
#include "TSystem.h"
#include "TObjString.h"
#include "strlcpy.h"

#include <sys/stat.h>

enum EFileFialog {
   kIDF_CDUP,
   kIDF_NEW_FOLDER,
   kIDF_LIST,
   kIDF_DETAILS,
   kIDF_CHECKB,
   kIDF_FSLB,
   kIDF_FTYPESLB,
   kIDF_OK,
   kIDF_CANCEL
};

static const char *gDefTypes[] = { "All files",     "*",
                                   "ROOT files",    "*.root",
                                   "ROOT macros",   "*.C",
                                    nullptr,         nullptr };

static TGFileInfo gInfo;


ClassImp(TGFileDialog);

////////////////////////////////////////////////////////////////////////////////
/// TGFileInfo Destructor.

TGFileInfo::~TGFileInfo()
{
   delete [] fFilename;
   delete [] fIniDir;
   DeleteFileNamesList();
}


////////////////////////////////////////////////////////////////////////////////
/// Delete file names list

void TGFileInfo::DeleteFileNamesList()
{
   if (fFileNamesList) {
      fFileNamesList->Delete();
      delete fFileNamesList;
      fFileNamesList = nullptr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Turn on/off multiple selection.

void TGFileInfo::SetMultipleSelection(Bool_t option)
{
   if ( fMultipleSelection != option ) {
      fMultipleSelection = option;
      DeleteFileNamesList();
      if (fMultipleSelection)
         fFileNamesList = new TList();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set file name

void TGFileInfo::SetFilename(const char *fname)
{
   delete [] fFilename;
   fFilename = fname ? StrDup(fname) : nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Set directory name

void TGFileInfo::SetIniDir(const char *inidir)
{
   delete [] fIniDir;
   fIniDir = inidir ? StrDup(inidir) : nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a file selection dialog. Depending on the dlg_type it can be
/// used for opening or saving a file.
/// About the first two arguments, p is the parent Window, usually the
/// desktop (root) window, and main is the main (TGMainFrame) application
/// window (the one opening the dialog), onto which the dialog is
/// usually centered, and which is waiting for it to close.

TGFileDialog::TGFileDialog(const TGWindow *p, const TGWindow *main,
                           EFileDialogMode dlg_type, TGFileInfo *file_info) :
   TGTransientFrame(p, main, 10, 10, kVerticalFrame), fTbfname(0), fName(0),
   fTypes(0), fTreeLB(0), fCdup(0), fNewf(0), fList(0), fDetails(0), fCheckB(0),
   fPcdup(0), fPnewf(0), fPlist(0), fPdetails(0), fOk(0), fCancel(0), fFv(0),
   fFc(0), fFileInfo(0)
{
   SetCleanup(kDeepCleanup);
   Connect("CloseWindow()", "TGFileDialog", this, "CloseWindow()");
   DontCallClose();

   int i;

   if (!p && !main) {
      MakeZombie();
      return;
   }
   if (!file_info) {
      Error("TGFileDialog", "file_info argument not set");
      fFileInfo = &gInfo;
      fFileInfo->SetIniDir(nullptr);
      fFileInfo->SetFilename(nullptr);
      fFileInfo->fFileTypeIdx = 0;
   } else
      fFileInfo = file_info;

   if (!fFileInfo->fFileTypes)
      fFileInfo->fFileTypes = gDefTypes;

   if (!fFileInfo->fIniDir)
      fFileInfo->SetIniDir(".");

   TGHorizontalFrame *fHtop = new TGHorizontalFrame(this, 10, 10);

   //--- top toolbar elements
   TGLabel *fLookin = new TGLabel(fHtop, new TGHotString((dlg_type == kFDSave)
                                                  ? "S&ave in:" : "&Look in:"));
   fTreeLB = new TGFSComboBox(fHtop, kIDF_FSLB);
   fTreeLB->Associate(this);

   fPcdup = fClient->GetPicture("tb_uplevel.xpm");
   fPnewf = fClient->GetPicture("tb_newfolder.xpm");
   fPlist = fClient->GetPicture("tb_list.xpm");
   fPdetails = fClient->GetPicture("tb_details.xpm");

   if (!(fPcdup && fPnewf && fPlist && fPdetails))
      Error("TGFileDialog", "missing toolbar pixmap(s).\n");

   fCdup    = new TGPictureButton(fHtop, fPcdup, kIDF_CDUP);
   fNewf    = new TGPictureButton(fHtop, fPnewf, kIDF_NEW_FOLDER);
   fList    = new TGPictureButton(fHtop, fPlist, kIDF_LIST);
   fDetails = new TGPictureButton(fHtop, fPdetails, kIDF_DETAILS);

   fCdup->SetStyle(gClient->GetStyle());
   fNewf->SetStyle(gClient->GetStyle());
   fList->SetStyle(gClient->GetStyle());
   fDetails->SetStyle(gClient->GetStyle());

   fCdup->SetToolTipText("Up One Level");
   fNewf->SetToolTipText("Create New Folder");
   fList->SetToolTipText("List");
   fDetails->SetToolTipText("Details");

   fCdup->Associate(this);
   fNewf->Associate(this);
   fList->Associate(this);
   fDetails->Associate(this);

   fList->AllowStayDown(kTRUE);
   fDetails->AllowStayDown(kTRUE);

   fTreeLB->Resize(200, fTreeLB->GetDefaultHeight());

   fHtop->AddFrame(fLookin, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 2, 2));
   fHtop->AddFrame(fTreeLB, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 3, 0, 2, 2));
   fHtop->AddFrame(fCdup, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 0, 2, 2));
   fHtop->AddFrame(fNewf, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 0, 2, 2));
   fHtop->AddFrame(fList, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 0, 2, 2));
   fHtop->AddFrame(fDetails, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 8, 2, 2));

   if (dlg_type == kFDSave) {
      fCheckB = new TGCheckButton(fHtop, "&Overwrite", kIDF_CHECKB);
      fCheckB->SetToolTipText("Overwrite a file without displaying a message if selected");
   } else {
      fCheckB = new TGCheckButton(fHtop, "&Multiple files", kIDF_CHECKB);
      fCheckB->SetToolTipText("Allows multiple file selection when SHIFT is pressed");
      fCheckB->Connect("Toggled(Bool_t)","TGFileInfo",fFileInfo,"SetMultipleSelection(Bool_t)");
   }
   fHtop->AddFrame(fCheckB, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
   fCheckB->SetOn(fFileInfo->fMultipleSelection);
   AddFrame(fHtop, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 4, 4, 3, 1));

   //--- file view

   fFv = new TGListView(this, 400, 161);

   fFc = new TGFileContainer(fFv->GetViewPort(),
                             10, 10, kHorizontalFrame, fgWhitePixel);
   fFc->Associate(this);

   fFv->GetViewPort()->SetBackgroundColor(fgWhitePixel);
   fFv->SetContainer(fFc);
   fFv->SetViewMode(kLVList);
   fFv->SetIncrements(1, 19); // set vertical scroll one line height at a time

   TGTextButton** buttons = fFv->GetHeaderButtons();
   if (buttons) {
      buttons[0]->Connect("Clicked()", "TGFileContainer", fFc, "Sort(=kSortByName)");
      buttons[1]->Connect("Clicked()", "TGFileContainer", fFc, "Sort(=kSortByType)");
      buttons[2]->Connect("Clicked()", "TGFileContainer", fFc, "Sort(=kSortBySize)");
      buttons[3]->Connect("Clicked()", "TGFileContainer", fFc, "Sort(=kSortByOwner)");
      buttons[4]->Connect("Clicked()", "TGFileContainer", fFc, "Sort(=kSortByGroup)");
      buttons[5]->Connect("Clicked()", "TGFileContainer", fFc, "Sort(=kSortByDate)");
   }

   fFc->SetFilter(fFileInfo->fFileTypes[fFileInfo->fFileTypeIdx+1]);
   fFc->Sort(kSortByName);
   fFc->ChangeDirectory(fFileInfo->fIniDir);
   fFc->SetMultipleSelection(fFileInfo->fMultipleSelection);
   fTreeLB->Update(fFc->GetDirectory());

   fList->SetState(kButtonEngaged);

   AddFrame(fFv, new TGLayoutHints(kLHintsTop | kLHintsExpandX | kLHintsExpandY, 4, 4, 3, 1));

   if (dlg_type == kFDOpen) {
      fCheckB->Connect("Toggled(Bool_t)","TGFileContainer",fFc,"SetMultipleSelection(Bool_t)");
      fCheckB->Connect("Toggled(Bool_t)","TGFileContainer",fFc,"UnSelectAll()");
   }

   //--- file name and types

   TGHorizontalFrame *fHf = new TGHorizontalFrame(this, 10, 10);

   TGVerticalFrame *fVf = new TGVerticalFrame(fHf, 10, 10);

   TGHorizontalFrame *fHfname = new TGHorizontalFrame(fVf, 10, 10);

   TGLabel *fLfname = new TGLabel(fHfname, new TGHotString("File &name:"));
   fTbfname = new TGTextBuffer(1034);
   fName = new TGTextEntry(fHfname, fTbfname);
   fName->Resize(230, fName->GetDefaultHeight());
   fName->Associate(this);

   fHfname->AddFrame(fLfname, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 2, 2));
   fHfname->AddFrame(fName, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 0, 20, 2, 2));

   fVf->AddFrame(fHfname, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX));

   TGHorizontalFrame *fHftype = new TGHorizontalFrame(fVf, 10, 10);

   TGLabel *fLftypes = new TGLabel(fHftype, new TGHotString("Files of &type:"));
   fTypes = new TGComboBox(fHftype, kIDF_FTYPESLB);
   fTypes->Associate(this);
   fTypes->Resize(230, fName->GetDefaultHeight());

   TString s;
   for (i = 0; fFileInfo->fFileTypes[i] != 0; i += 2) {
      s.Form("%s (%s)", fFileInfo->fFileTypes[i], fFileInfo->fFileTypes[i+1]);
      fTypes->AddEntry(s.Data(), i);
   }
   fTypes->Select(fFileInfo->fFileTypeIdx);

   // Show all items in combobox without scrollbar
   //TGDimension fw = fTypes->GetListBox()->GetContainer()->GetDefaultSize();
   //fTypes->GetListBox()->Resize(fw.fWidth, fw.fHeight);

   if (fFileInfo->fFilename && fFileInfo->fFilename[0]) {
      fTbfname->AddText(0, fFileInfo->fFilename);
   } else {
      fTbfname->Clear();
      if (dlg_type == kFDSave) {
         fTbfname->AddText(0, "unnamed");
         fName->SelectAll();
         if (fFileInfo->fFileTypes[fFileInfo->fFileTypeIdx+1] &&
             strstr(fFileInfo->fFileTypes[fFileInfo->fFileTypeIdx+1], "*.")) {
            TString ext = fFileInfo->fFileTypes[fFileInfo->fFileTypeIdx+1];
            ext.ReplaceAll("*.", ".");
            fTbfname->AddText(7, ext.Data());
         }
         fName->SetFocus();
      }
   }

   fTypes->GetListBox()->Resize(230, 120);
   fHftype->AddFrame(fLftypes, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 2, 2));
   fHftype->AddFrame(fTypes, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 0, 20, 2, 2));

   fVf->AddFrame(fHftype, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX));

   fHf->AddFrame(fVf, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX));

   //--- Open/Save and Cancel buttons

   TGVerticalFrame *fVbf = new TGVerticalFrame(fHf, 10, 10, kFixedWidth);

   fOk = new TGTextButton(fVbf, new TGHotString((dlg_type == kFDSave)
                                                 ? "&Save" : "&Open"), kIDF_OK);
   fCancel = new TGTextButton(fVbf, new TGHotString("Cancel"), kIDF_CANCEL);

   fOk->Associate(this);
   fCancel->Associate(this);

   fVbf->AddFrame(fOk, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 2, 2));
   fVbf->AddFrame(fCancel, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 2, 2));

   UInt_t width = TMath::Max(fOk->GetDefaultWidth(), fCancel->GetDefaultWidth()) + 20;
   fVbf->Resize(width + 20, fVbf->GetDefaultHeight());

   fHf->AddFrame(fVbf, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

   AddFrame(fHf, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 4, 4, 3, 1));
   SetEditDisabled(kEditDisable);

   MapSubwindows();

   TGDimension size = GetDefaultSize();

   Resize(size);

   //---- position relative to the parent's window

   CenterOnParent();

   //---- make the message box non-resizable

   SetWMSize(size.fWidth, size.fHeight);
   SetWMSizeHints(size.fWidth, size.fHeight, 10000, 10000, 1, 1);

   const char *wname = (dlg_type == kFDSave) ? "Save As..." : "Open";
   SetWindowName(wname);
   SetIconName(wname);
   SetClassHints("ROOT", "FileDialog");

   SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
                              kMWMDecorMinimize | kMWMDecorMenu,
               kMWMFuncAll |  kMWMFuncResize    | kMWMFuncMaximize |
                              kMWMFuncMinimize,
               kMWMInputModeless);

   MapWindow();
   fFc->DisplayDirectory();
   if (dlg_type == kFDSave)
      fName->SetFocus();
   fClient->WaitFor(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete file dialog.

TGFileDialog::~TGFileDialog()
{
   if (IsZombie()) return;
   TString str = fCheckB->GetString();
   if (str.Contains("Multiple"))
      fCheckB->Disconnect("Toggled(Bool_t)");
   fClient->FreePicture(fPcdup);
   fClient->FreePicture(fPnewf);
   fClient->FreePicture(fPlist);
   fClient->FreePicture(fPdetails);
   delete fFc;
}

////////////////////////////////////////////////////////////////////////////////
/// Close file dialog.

void TGFileDialog::CloseWindow()
{
   fFileInfo->SetFilename(nullptr);
   fFileInfo->DeleteFileNamesList();
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Small function used to prevent memory leaks with TSystem::ExpandPathName,
/// which returns a string created by StrDup, that has to be deleted

namespace {
   static inline void pExpandUnixPathName(TGFileInfo &file_info) {
      char *tmpPath = gSystem->ExpandPathName(file_info.fFilename);
      file_info.SetFilename(gSystem->UnixPathName(tmpPath));
      delete[] tmpPath;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Process messages generated by the user input in the file dialog.

Bool_t TGFileDialog::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   if (!fFc->GetDisplayStat()) return kTRUE;  // Cancel button was pressed

   TGTreeLBEntry *e;
   TGTextLBEntry *te;
   TGFileItem *f;
   void *p = 0;
   TString txt;
   TString sdir = gSystem->WorkingDirectory();

   switch (GET_MSG(msg)) {
      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
               switch (parm1) {
                  case kIDF_OK:
                     // same code as under kTE_ENTER
                     if (fTbfname->GetTextLength() == 0) {
                        txt = "Please provide file name or use \"Cancel\"";
                        new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                                     "Missing File Name", txt, kMBIconExclamation,
                                     kMBOk);
                        return kTRUE;
                     } else if (!gSystem->AccessPathName(fTbfname->GetString(), kFileExists) &&
                                !strcmp(fOk->GetTitle(), "Save") &&
                                (!(fCheckB->GetState() == kButtonDown))) {
                        Int_t ret;
                        txt = TString::Format("File name %s already exists, OK to overwrite it?",
                                              fTbfname->GetString());
                        new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                                     "File Name Exist", txt.Data(), kMBIconExclamation,
                                     kMBYes | kMBNo, &ret);
                        if (ret == kMBNo)
                           return kTRUE;
                     }
                     if (fFileInfo->fMultipleSelection) {
                        fFileInfo->SetFilename(nullptr);
                     } else {
                        fFileInfo->SetFilename(nullptr);
                        // FIXME: once appropriate gSystem method exists, use SetFilename here
                        if (gSystem->IsAbsoluteFileName(fTbfname->GetString()))
                           fFileInfo->SetFilename(fTbfname->GetString());
                        else
                           fFileInfo->fFilename = gSystem->ConcatFileName(fFc->GetDirectory(),
                                                                          fTbfname->GetString());
                        pExpandUnixPathName(*fFileInfo);
                     }
                     if (fCheckB && (fCheckB->GetState() == kButtonDown))
                        fFileInfo->fOverwrite = kTRUE;
                     else
                        fFileInfo->fOverwrite = kFALSE;
                     DeleteWindow();
                     break;

                  case kIDF_CANCEL:
                     fFileInfo->SetFilename(nullptr);
                     if (fFc->GetDisplayStat())
                        fFc->SetDisplayStat(kFALSE);
                     fFileInfo->DeleteFileNamesList();
                     DeleteWindow();
                     return kTRUE;   //no need to redraw fFc
                     break;

                  case kIDF_CDUP:
                     fFc->ChangeDirectory("..");
                     fTreeLB->Update(fFc->GetDirectory());
                     fFileInfo->SetIniDir(fFc->GetDirectory());
                     if (strcmp(gSystem->WorkingDirectory(),fFc->GetDirectory())) {
                        gSystem->cd(fFc->GetDirectory());
                     }
                     break;

                  case kIDF_NEW_FOLDER: {
                     char answer[128];
                     strlcpy(answer, "(empty)", sizeof(answer));
                     new TGInputDialog(gClient->GetRoot(), GetMainFrame(),
                                       "Enter directory name:",
                                       answer/*"(empty)"*/, answer);

                     while ( strcmp(answer, "(empty)") == 0 ) {
                        new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error",
                                     "Please enter a valid directory name.",
                                     kMBIconStop, kMBOk);
                        new TGInputDialog(gClient->GetRoot(), GetMainFrame(),
                                          "Enter directory name:",
                                          answer, answer);
                     }
                     if ( strcmp(answer, "") == 0 )  // Cancel button was pressed
                        break;

                     if (strcmp(gSystem->WorkingDirectory(),fFc->GetDirectory())) {
                        gSystem->cd(fFc->GetDirectory());
                     }
                     if ( gSystem->MakeDirectory(answer) != 0 )
                        new TGMsgBox(gClient->GetRoot(), GetMainFrame(), "Error",
                                     TString::Format("Directory name \'%s\' already exists!", answer),
                                     kMBIconStop, kMBOk);
                     else {
                        fFc->DisplayDirectory();
                     }
                     gSystem->ChangeDirectory(sdir.Data());
                     break;
                  }

                  case kIDF_LIST:
                     fFv->SetViewMode(kLVList);
                     fDetails->SetState(kButtonUp);
                     break;

                  case kIDF_DETAILS:
                     fFv->SetViewMode(kLVDetails);
                     fList->SetState(kButtonUp);
                     break;
               }
               break;

            case kCM_COMBOBOX:
               switch (parm1) {
                  case kIDF_FSLB:
                     e = (TGTreeLBEntry *) fTreeLB->GetSelectedEntry();
                     if (e) {
                        fFc->ChangeDirectory(e->GetPath()->GetString());
                        fTreeLB->Update(fFc->GetDirectory());
                        fFileInfo->SetIniDir(fFc->GetDirectory());
                        if (strcmp(gSystem->WorkingDirectory(),fFc->GetDirectory())) {
                           gSystem->cd(fFc->GetDirectory());
                        }
                     }
                     break;

                  case kIDF_FTYPESLB:
                     te = (TGTextLBEntry *) fTypes->GetSelectedEntry();
                     if (te) {
                        //fTbfname->Clear();
                        //fTbfname->AddText(0, fFileInfo->fFileTypes[te->EntryId()+1]);
                        fFileInfo->fFileTypeIdx = te->EntryId();
                        fFc->SetFilter(fFileInfo->fFileTypes[fFileInfo->fFileTypeIdx+1]);
                        fFc->DisplayDirectory();
                        fClient->NeedRedraw(fName);
                     }
                     break;
               }
               break;

            default:
               break;
         } // switch(GET_SUBMSG(msg))
         break;

      case kC_CONTAINER:
         switch (GET_SUBMSG(msg)) {
            case kCT_ITEMCLICK:
               if (parm1 == kButton1) {
                  if (fFc->NumSelected() > 0) {
                     if ( fFileInfo->fMultipleSelection == kFALSE ) {
                        TGLVEntry *e2 = (TGLVEntry *) fFc->GetNextSelected(&p);
                        if ((e2) && !R_ISDIR(((TGFileItem *)e2)->GetType())) {
                           fTbfname->Clear();
                           if (e2->GetItemName())
                              fTbfname->AddText(0, e2->GetItemName()->GetString());
                           fClient->NeedRedraw(fName);
                        }
                     }
                     else {
                        TString tmpString;
                        TList *tmp = fFc->GetSelectedItems();
                        TObjString *el;
                        TIter next(tmp);
                        if (fFileInfo->fFileNamesList) {
                           fFileInfo->fFileNamesList->Delete();
                        } else {
                           fFileInfo->fFileNamesList = new TList();
                        }
                        while ((el = (TObjString *) next())) {
                           char *s = gSystem->ConcatFileName(fFc->GetDirectory(),
                                                             el->GetString());
                           tmpString += "\"" + el->GetString() + "\" ";
                           fFileInfo->fFileNamesList->Add(new TObjString(s));
                           delete [] s;
                        }
                        tmp->Delete();
                        delete tmp;
                        fTbfname->Clear();
                        fTbfname->AddText(0, tmpString);
                        fClient->NeedRedraw(fName);
                     }
                  }
               }
               break;

            case kCT_ITEMDBLCLICK:

               if (parm1 == kButton1) {
                  if (fFc->NumSelected() == 1) {
                     f = (TGFileItem *) fFc->GetNextSelected(&p);
                     if (f && R_ISDIR(f->GetType())) {
                        fFc->ChangeDirectory(f->GetItemName()->GetString());
                        fTreeLB->Update(fFc->GetDirectory());
                        fFileInfo->SetIniDir(fFc->GetDirectory());
                        if (strcmp(gSystem->WorkingDirectory(),fFc->GetDirectory())) {
                           gSystem->cd(fFc->GetDirectory());
                        }
                     } else {
                        if (!strcmp(fOk->GetTitle(), "Save") &&
                            (!(fCheckB->GetState() == kButtonDown))) {

                           Int_t ret;
                           txt = TString::Format("File name %s already exists, OK to overwrite it?",
                                                 fTbfname->GetString());
                           new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                                        "File Name Exist", txt.Data(), kMBIconExclamation,
                                        kMBYes | kMBNo, &ret);
                           if (ret == kMBNo)
                              return kTRUE;
                        }
                        fFileInfo->SetFilename(nullptr);
                        // FIXME: once appropriate gSystem method exists, use SetFilename here
                        if (gSystem->IsAbsoluteFileName(fTbfname->GetString()))
                           fFileInfo->SetFilename(fTbfname->GetString());
                        else
                           fFileInfo->fFilename = gSystem->ConcatFileName(fFc->GetDirectory(),
                                                                          fTbfname->GetString());
                        pExpandUnixPathName(*fFileInfo);
                        if (fCheckB && (fCheckB->GetState() == kButtonDown))
                           fFileInfo->fOverwrite = kTRUE;
                        else
                           fFileInfo->fOverwrite = kFALSE;

                        DeleteWindow();
                     }
                  }
               }

               break;

            default:
               break;

         } // switch(GET_SUBMSG(msg))
         break;

      case kC_TEXTENTRY:
         switch (GET_SUBMSG(msg)) {
            case kTE_ENTER:
               // same code as under kIDF_OK
               if (fTbfname->GetTextLength() == 0) {
                  const char *txt2 = "Please provide file name or use \"Cancel\"";
                  new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                               "Missing File Name", txt2, kMBIconExclamation,
                               kMBOk);
                  return kTRUE;
               } else if (!gSystem->AccessPathName(fTbfname->GetString(), kFileExists)) {
                  FileStat_t buf;
                  if (!gSystem->GetPathInfo(fTbfname->GetString(), buf) &&
                      R_ISDIR(buf.fMode)) {
                     fFc->ChangeDirectory(fTbfname->GetString());
                     fTreeLB->Update(fFc->GetDirectory());
                     if (strcmp(gSystem->WorkingDirectory(), fFc->GetDirectory())) {
                        gSystem->cd(fFc->GetDirectory());
                     }
                     fName->SetText("", kFALSE);
                     return kTRUE;
                  }
                  else if (!strcmp(fOk->GetTitle(), "Save") &&
                          (!(fCheckB->GetState() == kButtonDown))) {
                     Int_t ret;
                     txt = TString::Format("File name %s already exists, OK to overwrite it?",
                                           fTbfname->GetString());
                     new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                                  "File Name Exist", txt.Data(), kMBIconExclamation,
                                  kMBYes | kMBNo, &ret);
                     if (ret == kMBNo)
                        return kTRUE;
                  }
               }
               fFileInfo->SetFilename(nullptr);
               // FIXME: once appropriate gSystem method exists, use SetFilename here
               fFileInfo->fFilename = gSystem->ConcatFileName(fFc->GetDirectory(),
                                                              fTbfname->GetString());
               pExpandUnixPathName(*fFileInfo);
               if (fCheckB && (fCheckB->GetState() == kButtonDown))
                  fFileInfo->fOverwrite = kTRUE;
               else
                  fFileInfo->fOverwrite = kFALSE;
               DeleteWindow();
               break;

            default:
               break;
         }
         break;

      default:
         break;

   } // switch(GET_MSG(msg))

   fClient->NeedRedraw(fFc);
   return kTRUE;
}
