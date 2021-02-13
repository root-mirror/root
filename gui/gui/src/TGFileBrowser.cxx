// @(#)root/gui:$Id$
// Author: Bertrand Bellenot   26/09/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TGListTree.h"
#include "TGLayout.h"
#include "TGComboBox.h"
#include "TContextMenu.h"
#include "TGTextEntry.h"
#include "TGTab.h"
#include "TGLabel.h"
#include "TSystemDirectory.h"
#include "TGMimeTypes.h"
#include "TClass.h"
#include "TQClass.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TInterpreter.h"
#include "TRegexp.h"
#include "TEnv.h"
#include "TImage.h"
#include "TBrowser.h"
#include "TRemoteObject.h"
#include "TKey.h"
#include "TKeyMapFile.h"
#include "TVirtualPad.h"
#include "Getline.h"
#include "snprintf.h"
#include <ctime>
#include <cstring>
#include <cstdlib>

#include "TGFileBrowser.h"
#include "TRootBrowser.h"
#include "TGInputDialog.h"

#include "TVirtualPadEditor.h"
#include "TGedEditor.h"
#include "TBaseClass.h"
#include "TVirtualX.h"

#include "RConfigure.h"

#ifdef WIN32
const char rootdir[] = "\\";
#else
const char rootdir[] = "/";
#endif

const char *filters[] = {
   "",
   "*.*",
   "*.[C|c|h]*",
   "*.root",
   "*.txt"
};

//_____________________________________________________________________________
//
// TCursorSwitcher
//
// Helper class used to change the cursor in a method and restore the
// original one when going out of the method scope.
//_____________________________________________________________________________

///////////////////////////////////////////////////////////////////////////////
class TCursorSwitcher {
private:
   TGWindow *fW1;
   TGWindow *fW2;
public:
   TCursorSwitcher(TGWindow *w1, TGWindow *w2) : fW1(w1), fW2(w2) {
      if (w1) gVirtualX->SetCursor(w1->GetId(), gVirtualX->CreateCursor(kWatch));
      if (w2) gVirtualX->SetCursor(w2->GetId(), gVirtualX->CreateCursor(kWatch));
   }
   ~TCursorSwitcher() {
      if (fW1) gVirtualX->SetCursor(fW1->GetId(), gVirtualX->CreateCursor(kPointer));
      if (fW2) gVirtualX->SetCursor(fW2->GetId(), gVirtualX->CreateCursor(kPointer));
   }
};

//_____________________________________________________________________________
//
// TGFileBrowser
//
// System file browser, used as TRootBrowser plug-in.
// This class is the real core of the ROOT browser.
//_____________________________________________________________________________

ClassImp(TGFileBrowser);

////////////////////////////////////////////////////////////////////////////////
/// TGFileBrowser constructor.

TGFileBrowser::TGFileBrowser(const TGWindow *p, TBrowser* b, UInt_t w, UInt_t h)
   : TGMainFrame(p, w, h), TBrowserImp(b), fNewBrowser(0)
{
   if (p && p != gClient->GetDefaultRoot())
      fNewBrowser = (TRootBrowser *)p->GetMainFrame();
   if (fNewBrowser)
      fNewBrowser->SetActBrowser(this);
   CreateBrowser();
   Resize(w, h);
   if (fBrowser) Show();
}

////////////////////////////////////////////////////////////////////////////////
/// Create the actual file browser.

void TGFileBrowser::CreateBrowser()
{
   fCachedPic  = 0;
   SetCleanup(kDeepCleanup);

   fTopFrame = new TGHorizontalFrame(this, 100, 30);
   fDrawOption = new TGComboBox(fTopFrame, "");
   TGTextEntry *dropt_entry = fDrawOption->GetTextEntry();
   dropt_entry->SetToolTipText("Object Draw Option", 300);
   fDrawOption->Resize(80, 20);
   TGListBox *lb = fDrawOption->GetListBox();
   lb->Resize(lb->GetWidth(), 120);
   Int_t dropt = 1;
   fDrawOption->AddEntry("", dropt++);
   fDrawOption->AddEntry("box", dropt++);
   fDrawOption->AddEntry("colz", dropt++);
   fDrawOption->AddEntry("lego", dropt++);
   fDrawOption->AddEntry("lego1", dropt++);
   fDrawOption->AddEntry("lego2", dropt++);
   fDrawOption->AddEntry("same", dropt++);
   fDrawOption->AddEntry("surf", dropt++);
   fDrawOption->AddEntry("surf1", dropt++);
   fDrawOption->AddEntry("surf2", dropt++);
   fDrawOption->AddEntry("surf3", dropt++);
   fDrawOption->AddEntry("surf4", dropt++);
   fDrawOption->AddEntry("surf5", dropt++);
   fDrawOption->AddEntry("text", dropt++);
   fTopFrame->AddFrame(fDrawOption, new TGLayoutHints(kLHintsCenterY |
                       kLHintsRight, 2, 2, 2, 2));
   fTopFrame->AddFrame(new TGLabel(fTopFrame, "Draw Option:"),
                       new TGLayoutHints(kLHintsCenterY | kLHintsRight,
                       5, 2, 2, 2));

   fSortButton = new TGPictureButton(fTopFrame, "bld_sortup.png");
   fSortButton->SetStyle(gClient->GetStyle());
   fSortButton->SetToolTipText("Sort Alphabetically\n(Current folder only)");
   fTopFrame->AddFrame(fSortButton, new TGLayoutHints(kLHintsCenterY |
                       kLHintsLeft, 2, 2, 2, 2));
   fSortButton->Connect("Clicked()", "TGFileBrowser", this, "ToggleSort()");

   fFilterButton = new TGPictureButton(fTopFrame, "filter.png");
   fFilterButton->SetStyle(gClient->GetStyle());
   fFilterButton->SetToolTipText("Filter Content");
   fTopFrame->AddFrame(fFilterButton, new TGLayoutHints(kLHintsCenterY |
                       kLHintsLeft, 2, 2, 2, 2));
   fFilterButton->Connect("Clicked()", "TGFileBrowser", this, "RequestFilter()");

   fRefreshButton = new TGPictureButton(fTopFrame, "refresh.png");
   fRefreshButton->SetStyle(gClient->GetStyle());
   fRefreshButton->SetToolTipText("Refresh Current Folder");
   fTopFrame->AddFrame(fRefreshButton, new TGLayoutHints(kLHintsCenterY |
                       kLHintsLeft, 2, 5, 2, 2));
   fRefreshButton->Connect("Clicked()", "TGFileBrowser", this, "Refresh()");

   AddFrame(fTopFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop |
            kLHintsExpandX, 2, 2, 2, 2));
   fCanvas   = new TGCanvas(this, 100, 100);
   fListTree = new TGListTree(fCanvas, kHorizontalFrame);
   AddFrame(fCanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop |
                kLHintsExpandX | kLHintsExpandY));
   fListTree->Connect("DoubleClicked(TGListTreeItem *, Int_t)",
      "TGFileBrowser", this, "DoubleClicked(TGListTreeItem *, Int_t)");
   fListTree->Connect("Clicked(TGListTreeItem *, Int_t, Int_t, Int_t)",
      "TGFileBrowser", this, "Clicked(TGListTreeItem *, Int_t, Int_t, Int_t)");
   fListTree->Connect("Checked(TObject*, Bool_t)", "TGFileBrowser",
      this, "Checked(TObject*, Bool_t)");

   fRootIcon = gClient->GetPicture("rootdb_t.xpm");
   fFileIcon = gClient->GetPicture("doc_t.xpm");

   fBotFrame = new TGHorizontalFrame(this, 100, 30);
   fBotFrame->AddFrame(new TGLabel(fBotFrame, "Filter: "),
                       new TGLayoutHints(kLHintsCenterY | kLHintsLeft,
                       2, 2, 2, 2));
   fFileType = new TGComboBox(fBotFrame, " All Files (*.*)");
   Int_t ftype = 1;
   fFileType->AddEntry(" All Files (*.*)", ftype++);
   fFileType->AddEntry(" C/C++ Files (*.c;*.cxx;*.h;...)", ftype++);
   fFileType->AddEntry(" ROOT Files (*.root)", ftype++);
   fFileType->AddEntry(" Text Files (*.txt)", ftype++);
   fFileType->Resize(200, 20);
   fBotFrame->AddFrame(fFileType, new TGLayoutHints(kLHintsLeft | kLHintsTop |
                kLHintsExpandX, 2, 2, 2, 2));
   fFileType->Connect("Selected(Int_t)", "TGFileBrowser", this,
                      "ApplyFilter(Int_t)");
   fFileType->GetTextEntry()->Connect("ReturnPressed()", "TGFileBrowser",
                                      this, "ApplyFilter(Int_t = -1)");
   AddFrame(fBotFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop |
            kLHintsExpandX, 2, 2, 2, 2));

   fContextMenu = new TContextMenu("FileBrowserContextMenu");
   fFilter      = 0;
   fGroupSize   = 1000;
   fListLevel   = 0;
   fCurrentDir  = 0;
   fRootDir     = 0;
   fDir         = 0;
   fFile        = 0;
   fNKeys       = 0;
   fCnt         = 0;
   fFilterStr   = "*";

   TString gv = gEnv->GetValue("Browser.GroupView", "1000");
   Int_t igv = atoi(gv.Data());
   if (igv > 10)
      fGroupSize = igv;

   if (gEnv->GetValue("Browser.ShowHidden", 0))
      fShowHidden = kTRUE;
   else
      fShowHidden = kFALSE;

   fDblClick = kFALSE;

   if (TClass::GetClass("TGHtmlBrowser"))
      TQObject::Connect("TGHtmlBrowser", "Clicked(char*)",
                        "TGFileBrowser", this, "Selected(char*)");

   TQObject::Connect("TPad", "Modified()",
                     "TGFileBrowser", this, "PadModified()");

   fListLevel = 0;
   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGFileBrowser::~TGFileBrowser()
{
   if (TClass::GetClass("TGHtmlBrowser"))
      TQObject::Disconnect("TGHtmlBrowser", "Clicked(char*)");
   TQObject::Disconnect("TPad", "Modified()");

   delete fContextMenu;
   delete fListTree;
   if (fRootIcon) fClient->FreePicture(fRootIcon);
   if (fCachedPic && (fCachedPic != fFileIcon))
      fClient->FreePicture(fCachedPic);
   if (fFileIcon) fClient->FreePicture(fFileIcon);
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Helper function checking if a class has a graphic properties editor.

static Bool_t IsObjectEditable(TClass *cl)
{
   TBaseClass *base;
   TList* bcl = cl->GetListOfBases();
   TIter next(bcl);
   while ((base = (TBaseClass*) next())) {
      cl = base->GetClassPointer();
      if (cl && TClass::GetClass(Form("%sEditor", cl->GetName())))
         return kTRUE;
      if (IsObjectEditable(cl))
         return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Format the tooltip information, based on the object passed in argument.

static const char *FormatToolTip(TObject *obj, Int_t maxlen=0)
{
   static TString infos;
   if (!obj) {
      infos.Clear();
      return 0;
   }
   infos = obj->GetName();
   if (obj->GetTitle()) {
      infos += "\n";
      infos += obj->GetTitle();
   }
   if (maxlen > 0 && infos.Length() > maxlen) {
      infos.Remove(maxlen - 3);
      infos += "...";
   }
   TString objinfo = obj->GetObjectInfo(1, 1);
   if (!objinfo.IsNull() && !objinfo.BeginsWith("x=")) {
      Long64_t bsize, fsize, objsize;
      objsize = objinfo.Atoll();
      if (objsize > 0) {
         infos += "\n";
         bsize = fsize = objsize;
         if (fsize > 1024) {
            fsize /= 1024;
            if (fsize > 1024) {
               // 3.7MB is more informative than just 3MB
               infos += TString::Format("Size: %lld.%lldM", fsize/1024,
                                        (fsize%1024)/103);
            } else {
               infos += TString::Format("Size: %lld.%lldK", bsize/1024,
                                        (bsize%1024)/103);
            }
         } else {
            infos += TString::Format("Size: %lld bytes", bsize);
         }
      }
   }
   return infos.Data();
}

/**************************************************************************/
// TBrowserImp virtuals
/**************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Add items to the browser. This function has to be called
/// by the Browse() member function of objects when they are
/// called by a browser. If check < 0 (default) no check box is drawn,
/// if 0 then unchecked checkbox is added, if 1 checked checkbox is added.

void TGFileBrowser::Add(TObject *obj, const char *name, Int_t check)
{
   if (fListLevel && !strcmp(fListLevel->GetText(), "Classes") &&
      fListLevel->GetParent() &&
      !strcmp(fListLevel->GetParent()->GetText(), "root")) {
      // Browsing list of root classes...
   }
   else {
      if (obj && obj->InheritsFrom("TApplication"))
         fListLevel = 0;
      if (obj && obj->InheritsFrom("TSystemDirectory"))
         return;
   }
   if (fListLevel) {
      TString oname = "";
      if (name)
         oname = name;
      else if (obj)
         oname = obj->GetName();
      // check if the current item is filtered
      auto it = fFilteredItems.find(fListLevel);
      if  (it != fFilteredItems.end()) {
         // check if the item (object) name match the filter
         const char *filter = it->second.c_str();
         TRegexp re(filter, kTRUE);
         // if not, simply return, so no item will be added
         if (oname.Index(re) == kNPOS) return;
      }
   }
   const TGPicture *pic = nullptr;
   if (obj && obj->InheritsFrom("TKey") && (obj->IsA() != TClass::Class()))
      AddKey(fListLevel, obj, name);
   else if (obj) {
      GetObjPicture(&pic, obj);
      if (!name) name = obj->GetName();
      if (check > -1) {
         if (!fListTree->FindChildByName(fListLevel, name)) {
            TGListTreeItem *item = fListTree->AddItem(fListLevel, name, obj,
                                                      pic, pic, kTRUE);
            if ((pic != fFileIcon) && (pic != fCachedPic))
               fClient->FreePicture(pic);
            if (item) fListTree->CheckItem(item, (Bool_t)check);
            fListTree->SetToolTipItem(item, FormatToolTip(obj, 32));
         }
      }
      else {
         // special case for remote object
         Bool_t isRemote = kFALSE;
         if (obj->InheritsFrom("TRemoteObject"))
            isRemote = kTRUE;
         else if (fListLevel) {
            // check also if one of its parents is a remote object
            TGListTreeItem *top = fListLevel;
            while (top->GetParent()) {
               TObject *tobj = (TObject *) top->GetUserData();
               if (tobj && (tobj->InheritsFrom("TRemoteObject") ||
                  tobj->InheritsFrom("TApplicationRemote"))) {
                  isRemote = kTRUE;
                  break;
               }
               top = top->GetParent();
            }
         }
         if (isRemote) {
            TRemoteObject *robj = (TRemoteObject *)obj;
            if (!strcmp(robj->GetClassName(), "TKey")) {
               AddKey(fListLevel, obj, name);
            }
            else {
               TString fname = name;
               // add the remote object only if not already in the list
               if (!fShowHidden && fname.BeginsWith("."))
                  return;
               AddRemoteFile(obj);
            }
         }
         else {
            if (!fListTree->FindChildByName(fListLevel, name)) {
               TGListTreeItem *item = fListTree->AddItem(fListLevel, name, obj, pic, pic);
               if ((pic != fFileIcon) && (pic != fCachedPic))
                  fClient->FreePicture(pic);
               if (item && obj && obj->InheritsFrom("TObject"))
                  item->SetDNDSource(kTRUE);
               fListTree->SetToolTipItem(item, FormatToolTip(obj, 32));
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add remote file in list tree.

void TGFileBrowser::AddRemoteFile(TObject *obj)
{
   Bool_t      is_link;
   Int_t       type;
   TString     filename;
   const TGPicture *spic;
   TGPicture *pic;

   FileStat_t sbuf;

   type    = 0;
   is_link = kFALSE;

   TRemoteObject *robj = (TRemoteObject *)obj;

   robj->GetFileStat(&sbuf);
   is_link = sbuf.fIsLink;
   type    = sbuf.fMode;
   filename = robj->GetName();
   if (R_ISDIR(type) || fFilter == 0 ||
       (fFilter && filename.Index(*fFilter) != kNPOS)) {

      GetFilePictures(&spic, type, is_link, filename);

      pic = (TGPicture*)spic; pic->AddReference();

      if ((!fListTree->FindChildByName(fListLevel, filename)) &&
         (!fListTree->FindChildByData(fListLevel, obj)))
         fListTree->AddItem(fListLevel, filename, obj, pic, pic);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Browse object. This, in turn, will trigger the calling of
/// TBrowser::Add() which will fill the IconBox and the tree.
/// Emits signal "BrowseObj(TObject*)".

void TGFileBrowser::BrowseObj(TObject *obj)
{
   if (fNewBrowser)
      fNewBrowser->SetActBrowser(this);
   if (obj != gROOT) {
      if (!fListTree->FindItemByObj(fListTree->GetFirstItem(), obj)) {
         fListLevel = 0;
         Add(obj);
         fListLevel = fListTree->FindItemByObj(fListTree->GetFirstItem(), obj);
         fListTree->HighlightItem(fListLevel);
         if (obj->IsFolder())
            fListTree->OpenItem(fListLevel);
         fListTree->ClearViewPort();
         fListTree->AdjustPosition(fListLevel);
      }
   }
   if (fBrowser) obj->Browse(fBrowser);
   if (obj == gROOT) {
      TList *volumes = gSystem->GetVolumes("all");
      TList *curvol  = gSystem->GetVolumes("cur");
      if (volumes && curvol) {
         const char *curdrive;
         TNamed *named = (TNamed *)curvol->At(0);
         if (named)
            curdrive = named->GetName();
         else
            curdrive = "C:";
         TIter next(volumes);
         TNamed *drive;
         while ((drive = (TNamed *)next())) {
            AddFSDirectory(TString::Format("%s\\", drive->GetName()), drive->GetTitle(),
                           (strcmp(drive->GetName(), curdrive) == 0) ?
                           "SetRootDir" : "Add");
         }
         delete volumes;
         delete curvol;
      }
      else {
         AddFSDirectory("/");
      }
      GotoDir(gSystem->WorkingDirectory());
      if (gROOT->GetListOfFiles() && !gROOT->GetListOfFiles()->IsEmpty())
         Selected(0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Emits signal when double clicking on icon.

void TGFileBrowser::Checked(TObject *obj, Bool_t checked)
{
   if (fNewBrowser)
      fNewBrowser->Checked(obj, checked);
}

////////////////////////////////////////////////////////////////////////////////
/// returns drawing option

Option_t *TGFileBrowser::GetDrawOption() const
{
   return fDrawOption->GetTextEntry()->GetText();
}

////////////////////////////////////////////////////////////////////////////////
/// Determine the file picture for the given file type.

void TGFileBrowser::GetFilePictures(const TGPicture **pic, Int_t file_type,
                                    Bool_t is_link, const char *name)
{
   static TString cached_ext;
   static const TGPicture *cached_spic = 0;
   const char *ext = name ? strrchr(name, '.') : 0;
   TString sname = name ? name : " ";
   *pic = 0;

   if (ext && cached_spic && (cached_ext == ext)) {
      *pic = cached_spic;
      return;
   }

   if (R_ISREG(file_type)) {
      *pic = gClient->GetMimeTypeList()->GetIcon(name, kTRUE);

      if (*pic) {
         if (ext) {
            cached_ext = ext;
            cached_spic = *pic;
            return;
         }
      }
   } else {
      *pic = 0;
   }

   if (*pic == 0) {
      *pic = gClient->GetPicture("doc_t.xpm");

      if (R_ISREG(file_type) && (file_type) & kS_IXUSR) {
         *pic = gClient->GetPicture("app_t.xpm");
      }
      if (R_ISDIR(file_type)) {
         *pic = gClient->GetPicture("folder_t.xpm");
      }
      if(sname.EndsWith(".root")) {
         *pic = gClient->GetPicture("rootdb_t.xpm");
      }

   }
   if (is_link) {
      *pic = gClient->GetPicture("slink_t.xpm");
   }

   cached_spic = 0;
   cached_ext = "";
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively remove object.

void TGFileBrowser::RecursiveRemove(TObject *obj)
{
   TGListTreeItem *itm = nullptr, *item = nullptr;
   if (obj->InheritsFrom("TFile")) {
      itm = fListTree->FindChildByData(0, gROOT->GetListOfFiles());
      if (itm)
         item = fListTree->FindChildByData(itm, obj);
      if (item) {
         // if the item to be deleted has a filter,
         // delete its entry in the map
         if (CheckFiltered(item))
            fFilteredItems.erase(item);
         fListTree->DeleteItem(item);
      }
      itm = fRootDir ? fRootDir->GetFirstChild() : nullptr;
      while (itm) {
         item = fListTree->FindItemByObj(itm, obj);
         if (item) {
            fListTree->DeleteChildren(item);
            item->SetUserData(0);
         }
         itm = itm->GetNextSibling();
      }
   }
   if (!obj->InheritsFrom("TFile") && fRootDir) {
      item = fListTree->FindItemByObj(fRootDir, obj);
      // if the item to be deleted has a filter, delete its entry in the map
      if (item && CheckFiltered(item))
         fFilteredItems.erase(item);
      fListTree->RecursiveDeleteItem(fRootDir, obj);
   }
   //fListTree->ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Refresh content of the list tree.

void TGFileBrowser::Refresh(Bool_t /*force*/)
{
   TTimer::SingleShot(200, "TGFileBrowser", this, "Update()");
   return; // disable refresh for the time being...
   // coverity[unreachable]
   TCursorSwitcher cursorSwitcher(this, fListTree);
   static UInt_t prev = 0;
   UInt_t curr =  gROOT->GetListOfBrowsables()->GetSize();
   if (!prev) prev = curr;

   if (prev != curr) { // refresh gROOT
      TGListTreeItem *sav = fListLevel;
      fListLevel = 0;
      BrowseObj(gROOT);
      fListLevel = sav;
      prev = curr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update content of the list tree.

void TGFileBrowser::Update()
{
   Long64_t size = 0;
   Long_t id = 0, flags = 0, modtime = 0;
   char path[1024];
   TGListTreeItem *item = fCurrentDir;
   if (!item) item = fRootDir;
   if (!item) return;
   //fListTree->DeleteChildren(item);
   TGListTreeItem *curr = fListTree->GetSelected(); // GetCurrent() ??
   if (curr) {
      TObject *obj = (TObject *) curr->GetUserData();
      if (obj && !obj->TestBit(kNotDeleted)) {
         // if the item to be deleted has a filter,
         // delete its entry in the map
         if (CheckFiltered(curr))
            fFilteredItems.erase(curr);
         fListTree->DeleteItem(curr);
         curr = 0;
         obj = 0;
      }
      else if (obj && obj->TestBit(kNotDeleted) &&
               obj->InheritsFrom("TObjString") && curr->GetParent()) {
         fListTree->GetPathnameFromItem(curr->GetParent(), path);
         if (strlen(path) > 1) {
            TString dirpath = FullPathName(curr->GetParent());
            Int_t res = gSystem->GetPathInfo(dirpath.Data(), &id, &size,
                                             &flags, &modtime);
            if ((res == 0) && (flags & 2)) {
               TString fullpath = FullPathName(curr);
               if (gSystem->AccessPathName(fullpath.Data())) {
                  // if the item to be deleted has a filter,
                  // delete its entry in the map
                  if (CheckFiltered(curr))
                     fFilteredItems.erase(curr);
                  fListTree->DeleteItem(curr);
                  curr = 0;
                  obj = 0;
               }
            }
         }
      }
   }
   TString actpath = FullPathName(item);
   flags = id = size = modtime = 0;
   if (gSystem->GetPathInfo(actpath.Data(), &id, &size, &flags, &modtime) == 0) {
      Int_t isdir = (Int_t)flags & 2;

      TString savdir = gSystem->WorkingDirectory();
      if (isdir) {
         TGListTreeItem *del = 0, *itm = item->GetFirstChild();
         while (itm) {
            fListTree->GetPathnameFromItem(itm, path);
            if (strlen(path) > 1) {
               TString recpath = FullPathName(itm);
               if (gSystem->AccessPathName(recpath.Data())) {
                  del = itm;
                  itm = itm->GetNextSibling();
                  // if the item to be deleted has a filter,
                  // delete its entry in the map
                  if (CheckFiltered(del))
                     fFilteredItems.erase(del);
                  fListTree->DeleteItem(del);
               }
            }
            if (del)
               del = 0;
            else
               itm = itm->GetNextSibling();
         }
      }
   }
   TGListTreeItem *sav = fListLevel;
   DoubleClicked(item, 1);
   fListLevel = sav;
   CheckFiltered(fListLevel, kTRUE);
}

/**************************************************************************/
// Other
/**************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Add file system directory in the list tree.

void TGFileBrowser::AddFSDirectory(const char *entry, const char *path,
                                   Option_t *opt)
{
   TGListTreeItem *item = 0;
   if ((opt == 0) || (!opt[0])) {
      if (fRootDir == 0 && !fListTree->FindChildByName(0, rootdir))
         item = fRootDir = fListTree->AddItem(0, rootdir);
      return;
   }
   if (strstr(opt, "SetRootDir")) {
      if (!fListTree->FindChildByName(0, entry))
         item = fRootDir = fListTree->AddItem(0, entry);
   }
   else if (strstr(opt, "Add")) {
      // MT: i give up! wanted to place entries for selected
      // directories like home, pwd, alice-macros.
      // TGListTreeItem *lti = fListTree->AddItem(0, entry);
      //
      if (!fListTree->FindChildByName(0, entry))
         item = fListTree->AddItem(0, entry);
   }
   if (item && path) {
      TString infos = path;
      item->SetTipText(path);
      TGPicture *pic = 0;
      if (infos.Contains("Removable"))
         pic = (TGPicture *)gClient->GetPicture("fdisk_t.xpm");
      else if (infos.Contains("Local"))
         pic = (TGPicture *)gClient->GetPicture("hdisk_t.xpm");
      else if (infos.Contains("CD"))
         pic = (TGPicture *)gClient->GetPicture("cdrom_t.xpm");
      else if (infos.Contains("Network"))
         pic = (TGPicture *)gClient->GetPicture("netdisk_t.xpm");
      if (pic)
         item->SetPictures(pic, pic);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// display content of ROOT file

void TGFileBrowser::AddKey(TGListTreeItem *itm, TObject *obj, const char *name)
{
   // Int_t from, to;
   TGListTreeItem *where;
   static TGListTreeItem *olditem = itm;
   static TGListTreeItem *item = itm;
   const TGPicture *pic;

   if (itm == 0) return;

   if ((fCnt == 0) || (olditem != itm)) {
      olditem = item = itm;
   }
   if (!name) name = obj->GetName();
   if (fNKeys > fGroupSize) {
      where = itm->GetFirstChild();
      while (where) {
         if (fListTree->FindItemByObj(where, obj))
            return;
         where = where->GetNextSibling();
      }
   }
   if ((fNKeys > fGroupSize) && (fCnt % fGroupSize == 0)) {
      if (item != itm) {
         TString newname = TString::Format("%s-%s", item->GetText(), name);
         item->Rename(newname.Data());
      }
      item = fListTree->AddItem(itm, name);
      item->SetDNDSource(kTRUE);
   }
   if ((fCnt > fGroupSize) && (fCnt >= fNKeys-1)) {
      TString newname = TString::Format("%s-%s", item->GetText(), name);
      item->Rename(newname.Data());
   }
   GetObjPicture(&pic, obj);
   if (!fListTree->FindChildByName(item, name)) {
      TGListTreeItem *it = fListTree->AddItem(item, name, obj, pic, pic);
      if (pic && (pic != fFileIcon) && (pic != fCachedPic))
         fClient->FreePicture(pic);
      it->SetDNDSource(kTRUE);
      it->SetTipText(FormatToolTip(obj, 32));
   }
   fCnt++;
}

////////////////////////////////////////////////////////////////////////////////
/// Apply filter selected in combo box to the file tree view.

void TGFileBrowser::ApplyFilter(Int_t id)
{
   // Long64_t size;
   // Long_t fid, flags, modtime;

   if (fFilter) delete fFilter;
   fFilter = 0;
   if ((id > 1) && (id < 5))
      fFilter = new TRegexp(filters[id], kTRUE);
   else if ((id < 0) || (id > 4)) {
      TGTextLBEntry *lbe = (TGTextLBEntry *)fFileType->GetSelectedEntry();
      if (lbe) {
         const char *text = lbe->GetTitle();
         fFilter = new TRegexp(text, kTRUE);
      }
   }
   TGListTreeItem *item = fCurrentDir;
   if (!item)
      item = fRootDir;
   if (!item) return;
   fListTree->DeleteChildren(item);
   DoubleClicked(item, 1);
   //fListTree->AdjustPosition(item);
   fListTree->ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Make object associated with item the current directory.

void TGFileBrowser::Chdir(TGListTreeItem *item)
{
   if (item) {
      TGListTreeItem *i = item;
      while (i) {
         TObject *obj = (TObject*) i->GetUserData();
         if ((obj) && obj->InheritsFrom("TDirectory")) {
            ((TDirectory *)obj)->cd();
            break;
         }
         i = i->GetParent();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the current list tree item points to a remote object.

void TGFileBrowser::CheckRemote(TGListTreeItem *item)
{
   if (!item) return;
   TObject *obj = (TObject *) item->GetUserData();
   if (obj) {
      if (obj->InheritsFrom("TApplicationRemote")) {
         if (!gApplication->GetAppRemote()) {
            gROOT->ProcessLine(TString::Format(".R %s", item->GetText()));
            if (gApplication->GetAppRemote()) {
               Getlinem(kInit, TString::Format("\n%s:root [0]",
                        gApplication->GetAppRemote()->ApplicationName()));
            }
         }
      }
      if (item->GetParent() && item->GetParent()->GetUserData() &&
         ((TObject *)item->GetParent()->GetUserData())->InheritsFrom("TApplicationRemote")) {
         // switch to remote session
         if (!gApplication->GetAppRemote()) {
            gROOT->ProcessLine(TString::Format(".R %s", item->GetParent()->GetText()));
            if (gApplication->GetAppRemote()) {
               Getlinem(kInit, TString::Format("\n%s:root [0]",
                        gApplication->GetAppRemote()->ApplicationName()));
            }
         }
         else if (!strcmp(item->GetText(), "ROOT Files")) {
            // update list of files opened in the remote session
            gApplication->SetBit(TApplication::kProcessRemotely);
            gApplication->ProcessLine("((TApplicationServer *)gApplication)->BrowseFile(0);");
         }
      }
      else {
         // check if the listtree item is from a local session or
         // from a remote session, then switch to the session it belongs to
         TGListTreeItem *top = item;
         while (top->GetParent()) {
            top = top->GetParent();
         }
         TObject *topobj = (TObject *) top->GetUserData();
         if (topobj && topobj->InheritsFrom("TApplicationRemote")) {
            // it belongs to a remote session
            if (!gApplication->GetAppRemote()) {
               // switch to remote session if not already in
               gROOT->ProcessLine(TString::Format(".R %s", top->GetText()));
               if (gApplication->GetAppRemote()) {
                  Getlinem(kInit, TString::Format("\n%s:root [0]",
                           gApplication->GetAppRemote()->ApplicationName()));
               }
            }
         }
         else if (gApplication->GetAppRemote()) {
            // switch back to local session if not already in
            gApplication->ProcessLine(".R");
            Getlinem(kInit, "\nroot [0]");
         }
      }
   }
   else if (gApplication->GetAppRemote()) {
      // switch back to local session if not already in
      gApplication->ProcessLine(".R");
      Getlinem(kInit, "\nroot [0]");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check if there is a filter active on the children of the list tree item.
/// If the but argument is true, the "filter" button state is set accordingly,
/// and its tooltip will show the filter used.

Bool_t TGFileBrowser::CheckFiltered(TGListTreeItem *item, Bool_t but)
{
   Bool_t found = kFALSE;
   TString filter;
   // if there is no filter (the map is empty) then just return
   if (fFilteredItems.empty())
      return kFALSE;
   auto it = fFilteredItems.find(item);
   if  (it != fFilteredItems.end()) {
      // if the item is in the map, take the filter regexp string
      filter = it->second.c_str();
      fFilterStr = filter;
      found = kTRUE;
   }
   if (but) {
      // if the but argument is true, change the button state
      // to reflect the filtering state
      fFilterButton->SetState(found ? kButtonEngaged : kButtonUp);
      if (found) {
         // format the tooltip to display the regexp used as filter
         filter.Prepend("Showing only \'");
         filter += "\'";
         fFilterButton->SetToolTipText(filter.Data());
      }
      else {
         // reset the tooltip text
         fFilterButton->SetToolTipText("Filter Content...");
      }
   }
   return found;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the list tree item children are alphabetically sorted.
/// If the but argument is true, the "sort" button state is set accordingly.

Bool_t TGFileBrowser::CheckSorted(TGListTreeItem *item, Bool_t but)
{
   Bool_t found = kFALSE;
   TGListTreeItem *i, *itm;
   if (item->GetFirstChild())
      itm = item;
   else
      itm = item->GetParent();
   for (sLTI_i p=fSortedItems.begin(); p!=fSortedItems.end(); ++p) {
      i = (TGListTreeItem *)(*p);
      if (itm == i) {
         found = kTRUE;
         break;
      }
   }
   if (but) fSortButton->SetState(found ? kButtonEngaged : kButtonUp);
   return found;
}

////////////////////////////////////////////////////////////////////////////////
/// Process mouse clicks in TGListTree.

void TGFileBrowser::Clicked(TGListTreeItem *item, Int_t btn, Int_t x, Int_t y)
{
   char path[1024];
   Long64_t size = 0;
   Long_t id = 0, flags = 0, modtime = 0;
   fListLevel = item;
   if (!item) return;
   CheckSorted(item, kTRUE);
   CheckFiltered(item, kTRUE);
   CheckRemote(item);
   TObject *selected = 0;
   TString fullpath = FullPathName(item);
   TObject *obj = (TObject *) item->GetUserData();
   if (obj && (!obj->InheritsFrom("TObjString") ||
       gSystem->AccessPathName(fullpath.Data()))) {
      if (obj->InheritsFrom("TKey") && (obj->IsA() != TClass::Class())) {
         Chdir(item);
         const char *clname = ((TKey *)obj)->GetClassName();
         if (clname && strcmp(clname, "TGeoManager")) {
            TClass *cl = TClass::GetClass(clname);
            TString name = ((TKey *)obj)->GetName();
            name += ";";
            name += ((TKey *)obj)->GetCycle();
            void *add = gDirectory->FindObjectAny((char *) name.Data());
            if (add && cl->IsTObject()) {
               obj = (TObject*)add;
               // don't change the user data, to avoid deletion of the
               // list tree item by RecursiveRemove()
               // it is better to read the object each time anyway,
               // as it may have changed in the file
               if (obj->InheritsFrom("TDirectory") || obj->InheritsFrom("TList"))
                  item->SetUserData(obj);
            }
         }
      }
      if (obj->InheritsFrom("TLeaf") ||
          obj->InheritsFrom("TBranch")) {
         Chdir(item);
      }
      if (btn == kButton3)
        fContextMenu->Popup(x, y, obj, fBrowser);
      selected = obj;
   }
   else {
      fListTree->GetPathnameFromItem(item, path);
      if (strlen(path) > 3) {
         if (gSystem->GetPathInfo(fullpath.Data(), &id, &size, &flags, &modtime) == 0) {
            if (flags & 2) {
               fCurrentDir = item;
               if (btn == kButton3) {
                  if (fDir) delete fDir;
                  fDir = new TSystemDirectory(item->GetText(), fullpath.Data());
                  fContextMenu->Popup(x, y, fDir, fBrowser);
               }
            }
            else {
               fCurrentDir = item->GetParent();
               if (btn == kButton3) {
                  if (fFile) delete fFile;
                  fFile = new TSystemFile(item->GetText(), fullpath.Data());
                  fContextMenu->Popup(x, y, fFile, fBrowser);
               }
            }
         }
      }
   }
   fListTree->ClearViewPort();
   if (selected && (selected->IsA() != TClass::Class())) {
      if (selected->InheritsFrom("TLeaf"))
         selected = (TObject *)gROOT->ProcessLine(TString::Format("((TLeaf *)0x%lx)->GetBranch()->GetTree();", (ULong_t)selected));
      if (selected->InheritsFrom("TBranch"))
         selected = (TObject *)gROOT->ProcessLine(TString::Format("((TBranch *)0x%lx)->GetTree();", (ULong_t)selected));
      if (selected->InheritsFrom("TTree")) {
         // if a tree not attached to any directory (e.g. in a TFolder)
         // then attach it to the current directory (gDirectory)
         TDirectory *tdir = (TDirectory *)gROOT->ProcessLine(TString::Format("((TTree *)0x%lx)->GetDirectory();", (ULong_t)selected));
         if (!tdir) {
            gROOT->ProcessLine(TString::Format("((TTree *)0x%lx)->SetDirectory(gDirectory);", (ULong_t)selected));
         }
      }
   }
   if (selected && gPad && IsObjectEditable(selected->IsA())) {
      TVirtualPadEditor *ved = TVirtualPadEditor::GetPadEditor(kFALSE);
      if (ved) {
         TGedEditor *ged = (TGedEditor *)ved;
         ged->SetModel(gPad, selected, kButton1Down);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// returns an absolute path

TString TGFileBrowser::FullPathName(TGListTreeItem* item)
{
   TGListTreeItem *parent, *itm = item;
   TString dirname = itm->GetText();

   while ((parent=itm->GetParent())) {
      char *s = gSystem->ConcatFileName(parent->GetText(), dirname);
      dirname = s;
      delete [] s;
      itm = parent;
   }
   gSystem->ExpandPathName(dirname);
#ifdef R__WIN32
   // only handle .lnk files on Windows
   while (dirname.Contains(".lnk")) {
      Ssiz_t idx = dirname.Index(".lnk") + 4;
      TString resolved = dirname;
      resolved.Remove(idx);
      gSystem->ExpandPathName(resolved);
      dirname = resolved.Append(dirname.Remove(0, idx));
   }
#endif
   return dirname;
}

////////////////////////////////////////////////////////////////////////////////
/// returns the directory path

TString TGFileBrowser::DirName(TGListTreeItem* item)
{
   TString dirname;
   TString fullpath = FullPathName(item);

#ifdef WIN32
   char   winDrive[256];
   char   winDir[256];
   char   winName[256];
   char   winExt[256];
   _splitpath(fullpath.Data(), winDrive, winDir, winName, winExt);
   dirname = TString::Format("%s%s", winDrive, winDir);
#else
   dirname = gSystem->GetDirName(fullpath);
#endif
   return dirname;
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

   if (gSystem->GetPathInfo(candidate, buf) || !(buf.fMode & kS_IFREG))
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
/// Create a symlink (shortcut on Windows) icon by merging the picture
/// passed as argument and the slink_t.xpm icon (small arrow)

static const TGPicture *MakeLinkPic(const TGPicture *pic)
{
   const TGPicture *merged;
   TImage *img1, *img2;
   if (pic) {
      img1 = TImage::Create();
      if (img1 == 0) return pic;
      img1->SetImage(((const TGPicture *)pic)->GetPicture(),
                     ((const TGPicture *)pic)->GetMask());
      img2 = TImage::Open("slink_t.xpm");
      if (img2) img1->Merge(img2);
      TString lnk_name = ((const TGPicture *)pic)->GetName();
      lnk_name.Prepend("lnk_");
      merged = gClient->GetPicturePool()->GetPicture(lnk_name.Data(),
                                          img1->GetPixmap(), img1->GetMask());
      if (img2) delete img2;
      delete img1;
      return merged;
   }
   return pic;
}

////////////////////////////////////////////////////////////////////////////////
/// Process double clicks in TGListTree.

void TGFileBrowser::DoubleClicked(TGListTreeItem *item, Int_t /*btn*/)
{
   const TGPicture *pic=0;
   TString dirname = DirName(item);
   TString fullpath = FullPathName(item);
   TGListTreeItem *itm;
   FileStat_t sbuf;
   Long64_t size;
   Long_t id, flags, modtime;
   char action[512];
   TString act;
   Bool_t is_link = kFALSE;
   if (!gSystem->GetPathInfo(item->GetText(), sbuf) && sbuf.fIsLink) {
      is_link = kTRUE;
      fullpath = item->GetText();
      gSystem->ExpandPathName(fullpath);
   }

   if (fNewBrowser)
      fNewBrowser->SetActBrowser(this);
   TCursorSwitcher switcher(this, fListTree);
   fListLevel = item;
   CheckSorted(item, kTRUE);
   CheckFiltered(item, kTRUE);
   CheckRemote(item);
   TGListTreeItem *pitem = item->GetParent();
   TObject *obj = (TObject *) item->GetUserData();
   if (obj && !obj->InheritsFrom("TSystemFile")) {
      TString ext = obj->GetName();
      if (obj->InheritsFrom("TDirectory") && (obj->IsA() != TClass::Class())) {
         if (((TDirectory *)obj)->GetListOfKeys())
            fNKeys = ((TDirectory *)obj)->GetListOfKeys()->GetEntries();
         else
            fNKeys = 0;
      }
      else if (obj->InheritsFrom("TKey") && (obj->IsA() != TClass::Class())) {
         Chdir(item);
         const char *clname = ((TKey *)obj)->GetClassName();
         if (clname) {
            TClass *cl = TClass::GetClass(clname);
            TString name = ((TKey *)obj)->GetName();
            name += ";";
            name += ((TKey *)obj)->GetCycle();
            void *add = gDirectory->FindObjectAny((char *) name.Data());
            if (add && cl->IsTObject()) {
               obj = (TObject*)add;
               // don't change the user data, to avoid deletion of the
               // list tree item by RecursiveRemove()
               // it is better to read the object each time anyway,
               // as it may have changed in the file
               if (obj->InheritsFrom("TDirectory") || obj->InheritsFrom("TList"))
                  item->SetUserData(obj);
            }
         }
      }
      else if (obj->InheritsFrom("TLeaf") || obj->InheritsFrom("TBranch")) {
         Chdir(item);
      }
      else if (obj->InheritsFrom("TRemoteObject")) {
         // the real object is a TKey
         TRemoteObject *robj = (TRemoteObject *)obj;
         if (!strcmp(robj->GetClassName(), "TKey")) {
            TGListTreeItem *parent = item;
            TRemoteObject *probj = (TRemoteObject *)parent->GetUserData();
            // find the TFile remote object containing the TKey
            while ( probj && strcmp(probj->GetClassName(), "TFile")) {
               parent = parent->GetParent();
               probj = (TRemoteObject *)parent->GetUserData();
            }
            if (probj && !strcmp(probj->GetClassName(), "TFile")) {
               // remotely browse file (remotely call TFile::cd())
               gApplication->SetBit(TApplication::kProcessRemotely);
               gApplication->ProcessLine(
                  TString::Format("((TApplicationServer *)gApplication)->BrowseFile(\"%s\");",
                       probj->GetName()));
               gSystem->Sleep(250);
            }
         }
         if (gClient->GetMimeTypeList()->GetAction(obj->GetName(), action)) {
            act = action;
            act.ReplaceAll("%s", obj->GetName());
            if ((act[0] != '!') && (strcmp(pitem->GetText(), "ROOT Files"))) {
               // special case for remote object: remote process
               gApplication->SetBit(TApplication::kProcessRemotely);
               gApplication->ProcessLine(act.Data());
            }
         }
         if ((ext.EndsWith(".root")) && (strcmp(pitem->GetText(), "ROOT Files"))) {
            gApplication->SetBit(TApplication::kProcessRemotely);
            gApplication->ProcessLine("((TApplicationServer *)gApplication)->BrowseFile(0);");
         }
      }
      if (!obj->InheritsFrom("TObjString") ||
          gSystem->AccessPathName(fullpath.Data())) {
         if (fBrowser) fBrowser->SetDrawOption(GetDrawOption());
         fDblClick = kTRUE;
         if (gClient->GetMimeTypeList()->GetAction(obj->IsA()->GetName(), action)) {
            act = action;
            if (fBrowser && act.Contains("->Browse()")) obj->Browse(fBrowser);
            else if (act.Contains("->Draw()")) obj->Draw(GetDrawOption());
            else {
               if (act.Contains("%s")) act.ReplaceAll("%s", obj->GetName());
               else act.Prepend(obj->GetName());
               gInterpreter->SaveGlobalsContext();
               if (act[0] == '!') {
                  act.Remove(0, 1);
                  gSystem->Exec(act.Data());
               } else {
                  // special case for remote object: remote process
                  if (obj->InheritsFrom("TRemoteObject"))
                     gApplication->SetBit(TApplication::kProcessRemotely);
                  gApplication->ProcessLine(act.Data());
               }
            }
         }
         else if (obj->InheritsFrom("TCanvas") && fNewBrowser &&
                  fNewBrowser->GetTabRight() &&
                  fNewBrowser->GetTabRight()->GetTabTab(obj->GetName())) {
            // avoid potential crash when drawing a canvas with the same name
            // than a canvas already embedded in one of the browser's tab
            obj->DrawClone();
         }
         else if (fBrowser && !obj->InheritsFrom("TFormula") &&
                  !obj->InheritsFrom("TMethodBrowsable"))
            obj->Browse(fBrowser);
         fDblClick = kFALSE;
         fNKeys = 0;
         fCnt = 0;
         fListTree->ClearViewPort();
         if (gPad) gPad->Update();
         return;
      }
   }
   flags = id = size = modtime = 0;
   if (gSystem->GetPathInfo(fullpath.Data(), &id, &size, &flags, &modtime) != 0)
      return;
   Int_t isdir = (Int_t)flags & 2;

   TString savdir = gSystem->WorkingDirectory();
   if (isdir) {
      fCurrentDir = item;
      //fListTree->DeleteChildren(item);
      TSystemDirectory dir(item->GetText(),FullPathName(item));
      TList *files = dir.GetListOfFiles();
      if (files) {
         files->Sort();
         TIter next(files);
         TSystemFile *file;
         TString fname, pname;
         // directories first
         //fListTree->DeleteChildren(item);
         while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            if (file->IsDirectory()) {
               if (!fShowHidden && fname.BeginsWith("."))
                  continue;
               if ((fname!="..") && (fname!=".")) { // skip it
                  if (!fListTree->FindChildByName(item, fname)) {
                     itm = fListTree->AddItem(item, fname);
                     if (!gSystem->GetPathInfo(fname, sbuf) &&
                         sbuf.fIsLink) {
                        // change the pictures if it is a symlink
                        // (shortcut on Windows)
                        const TGPicture *opened = 0, *l_opened = 0;
                        const TGPicture *closed = 0, *l_closed = 0;
                        opened = fClient->GetPicture("ofolder_t.xpm");
                        if (opened) l_opened = MakeLinkPic(opened);
                        closed = fClient->GetPicture("folder_t.xpm");
                        if (closed) l_closed = MakeLinkPic(closed);
                        if (l_opened && l_closed)
                           itm->SetPictures(l_opened, l_closed);
                        if (opened) fClient->FreePicture(opened);
                        if (closed) fClient->FreePicture(closed);
                        if (l_opened) fClient->FreePicture(l_opened);
                        if (l_closed) fClient->FreePicture(l_closed);
                     }
                     // uncomment line below to set directories as
                     // DND targets
                     //itm->SetDNDTarget(kTRUE);
                     itm->SetUserData(0);
                  }
               }
            }
         }
         // then files...
         TIter nextf(files);
         while ((file=(TSystemFile*)nextf())) {
            fname = pname = file->GetName();
            if (!file->IsDirectory() && (fFilter == 0 ||
               (fFilter && fname.Index(*fFilter) != kNPOS))) {
               if (!fShowHidden && fname.BeginsWith("."))
                  continue;
               size = modtime = 0;
               if (gSystem->GetPathInfo(fname, sbuf) == 0) {
                  size    = sbuf.fSize;
                  modtime = sbuf.fMtime;
               }
               if (sbuf.fIsLink && pname.EndsWith(".lnk"))
                  pname.Remove(pname.Length()-4);
               pic = gClient->GetMimeTypeList()->GetIcon(pname, kTRUE);
               if (!pic)
                  pic = fFileIcon;
               if (sbuf.fIsLink)
                  pic = MakeLinkPic(pic);
               if (!fListTree->FindChildByName(item, fname)) {
                  itm = fListTree->AddItem(item, fname, pic, pic);
                  if (pic != fFileIcon)
                     fClient->FreePicture(pic);
                  if (sbuf.fIsLink) {
                     TString fullname = file->GetName();
                     gSystem->ExpandPathName(fullname);
                     itm->SetUserData(new TObjString(TString::Format("file://%s\r\n",fullname.Data())), kTRUE);
                  } else {
                     itm->SetUserData(new TObjString(TString::Format("file://%s/%s\r\n",
                                      gSystem->UnixPathName(file->GetTitle()),
                                      file->GetName())), kTRUE);
                  }
                  itm->SetDNDSource(kTRUE);
                  if (size && modtime) {
                     char *tiptext = FormatFileInfo(fname.Data(), size, modtime);
                     itm->SetTipText(tiptext);
                     delete [] tiptext;
                  }
               }
            }
         }
         files->Delete();
         delete files;
      }
   }
   else {
      TString lnkname = item->GetText();
      if (is_link && lnkname.EndsWith(".lnk"))
         lnkname.Remove(lnkname.Length()-4);
      fCurrentDir = item->GetParent();
      TSystemFile f(lnkname.Data(), fullpath.Data());
      TString fname = f.GetName();
      if (fname.EndsWith(".root")) {
         TDirectory *rfile = 0;
         gSystem->ChangeDirectory(dirname.Data());
         rfile = (TDirectory *)gROOT->GetListOfFiles()->FindObject(obj);
         if (!rfile) {
            rfile = (TDirectory *)gROOT->ProcessLine(TString::Format("new TFile(\"%s\")",fname.Data()));
         }
         if (rfile) {
            // replace actual user data (TObjString) by the TDirectory...
            if (item->GetUserData()) {
               // first delete the data to avoid memory leaks
               TObject *obj2 = static_cast<TObject *>(item->GetUserData());
               // only delete TObjString as they are the only objects
               // created who have to be deleted
               TObjString *ostr = dynamic_cast<TObjString *>(obj2);
               if (ostr) delete ostr;
            }
            item->SetUserData(rfile);
            fNKeys = rfile->GetListOfKeys()->GetEntries();
            fCnt = 0;
            if (fBrowser) rfile->Browse(fBrowser);
            fNKeys = 0;
            fCnt = 0;
         }
      }
      else if (fname.EndsWith(".png")) {
         gSystem->ChangeDirectory(dirname.Data());
         XXExecuteDefaultAction(&f);
         gSystem->ChangeDirectory(savdir.Data());
      }
      else if (IsTextFile(fullpath.Data())) {
         gSystem->ChangeDirectory(dirname.Data());
         if (fNewBrowser) {
            TGFrameElement *fe = 0;
            TGTab *tabRight = fNewBrowser->GetTabRight();
            TGCompositeFrame *frame = tabRight->GetCurrentContainer();
            if (frame)
               fe = (TGFrameElement *)frame->GetList()->First();
            if (fe) {
               TGCompositeFrame *embed = (TGCompositeFrame *)fe->fFrame;
               TString fullname = f.GetTitle();
               fullname.ReplaceAll("\\", "\\\\");
               if (embed->InheritsFrom("TGTextEditor")) {
                  gROOT->ProcessLine(TString::Format("((TGTextEditor *)0x%lx)->LoadFile(\"%s\");",
                                     (ULong_t)embed, fullname.Data()));
               }
               else if (embed->InheritsFrom("TGTextEdit")) {
                  gROOT->ProcessLine(TString::Format("((TGTextEdit *)0x%lx)->LoadFile(\"%s\");",
                                     (ULong_t)embed, fullname.Data()));
               }
               else {
                  XXExecuteDefaultAction(&f);
               }
            }
            else {
               XXExecuteDefaultAction(&f);
            }
         }
         gSystem->ChangeDirectory(savdir.Data());
      }
      else {
         gSystem->ChangeDirectory(dirname.Data());
         XXExecuteDefaultAction(&f);
         gSystem->ChangeDirectory(savdir.Data());
      }
   }
   //gSystem->ChangeDirectory(savdir.Data());
   fListTree->ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Execute default action for selected object (action is specified
/// in the $HOME/.root.mimes or $ROOTSYS/etc/root.mimes file.

Long_t TGFileBrowser::XXExecuteDefaultAction(TObject *obj)
{
   char action[512];
   TString act;
   TString ext = obj->GetName();
   if (fBrowser) fBrowser->SetDrawOption(GetDrawOption());

   if (gClient->GetMimeTypeList()->GetAction(obj->GetName(), action)) {
      act = action;
      act.ReplaceAll("%s", obj->GetName());
      gInterpreter->SaveGlobalsContext();

      if (act[0] == '!') {
         act.Remove(0, 1);
         gSystem->Exec(act.Data());
         return 0;
      } else {
         // special case for remote object: remote process
         if (obj->InheritsFrom("TRemoteObject"))
            gApplication->SetBit(TApplication::kProcessRemotely);

         const Long_t res = gApplication->ProcessLine(act.Data());
#ifdef R__HAS_COCOA
         if (act.Contains(".x") || act.Contains(".X")) {
            if (gPad) gPad->Update();
         }
#endif
         return res;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Format file information to be displayed in the tooltip.

char *TGFileBrowser::FormatFileInfo(const char *fname, Long64_t size, Long_t modtime)
{
   Long64_t fsize, bsize;
   TString infos = fname;
   infos += "\n";

   fsize = bsize = size;
   if (fsize > 1024) {
      fsize /= 1024;
      if (fsize > 1024) {
         // 3.7MB is more informative than just 3MB
         infos += TString::Format("Size: %lld.%lldM", fsize/1024, (fsize%1024)/103);
      } else {
         infos += TString::Format("Size: %lld.%lldK", bsize/1024, (bsize%1024)/103);
      }
   } else {
      infos += TString::Format("Size: %lld", bsize);
   }
   struct tm *newtime;
   time_t loctime = (time_t) modtime;
   newtime = localtime(&loctime);
   if (newtime) {
      infos += "\n";
      infos += TString::Format("%d-%02d-%02d %02d:%02d",
                               newtime->tm_year + 1900,
                               newtime->tm_mon+1, newtime->tm_mday,
                               newtime->tm_hour, newtime->tm_min);
   }
   return StrDup(infos.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve icons associated with class "name". Association is made
/// via the user's ~/.root.mimes file or via $ROOTSYS/etc/root.mimes.

void TGFileBrowser::GetObjPicture(const TGPicture **pic, TObject *obj)
{
   const char *clname = 0;
   TClass *objClass = 0;
   static TImage *im = 0;
   if (!im) {
      im = TImage::Create();
   }

   if (obj->IsA() == TClass::Class()) {
      objClass = obj->IsA();
      if (objClass)
         clname = objClass->GetName();
   }
   else if (obj->InheritsFrom("TKey")) {
      clname = ((TKey *)obj)->GetClassName();
   }
   else if (obj->InheritsFrom("TKeyMapFile")) {
      clname = ((TKeyMapFile *)obj)->GetTitle();
   }
   else if (obj->InheritsFrom("TRemoteObject")) {
      // special case for remote object: get real object class
      TRemoteObject *robj = (TRemoteObject *)obj;
      if (!strcmp(robj->GetClassName(), "TKey"))
         clname = robj->GetKeyClassName();
      else
         clname = robj->GetClassName();
   }
   else {
      objClass = obj->IsA();
      if (objClass)
         clname = objClass->GetName();
   }
   if (!clname) {
      clname = "Unknown";
   }
   const char *name = obj->GetIconName() ? obj->GetIconName() : clname;
   TString xpm_magic(name, 3);
   Bool_t xpm = xpm_magic == "/* ";
   const char *iconname = xpm ? obj->GetName() : name;

   if (obj->IsA()->InheritsFrom("TGeoVolume")) {
      iconname = obj->GetIconName() ? obj->GetIconName() : obj->IsA()->GetName();
   }

   if (fCachedPicName == iconname) {
      *pic = fCachedPic;
      return;
   }
   *pic = gClient->GetMimeTypeList()->GetIcon(iconname, kTRUE);
   if (!(*pic) && xpm) {
      if (im && im->SetImageBuffer((char**)&name, TImage::kXpm)) {
         im->Scale(im->GetWidth()/4, im->GetHeight()/4);
         *pic = gClient->GetPicturePool()->GetPicture(iconname, im->GetPixmap(),
                                                      im->GetMask());
      }
      gClient->GetMimeTypeList()->AddType("[thumbnail]", iconname, iconname, iconname, "->Browse()");
      return;
   }
   if (fCachedPic && (fCachedPic != fFileIcon))
      fClient->FreePicture(fCachedPic);
   if (*pic == 0) {
      if (!obj->IsFolder())
         *pic = fFileIcon;
   }
   fCachedPic = *pic;
   fCachedPicName = iconname;
}

////////////////////////////////////////////////////////////////////////////////
/// Go to the directory "path" and open all the parent list tree items.

void TGFileBrowser::GotoDir(const char *path)
{
   TGListTreeItem *item, *itm;
   ULong_t id;
   Long_t bsize, blocks, bfree;
   Bool_t expand = kTRUE;
   TString sPath(gSystem->UnixPathName(path));
   item = fRootDir;
   if (item == 0) return;
   fListTree->OpenItem(item);
   TObjArray *tokens = sPath.Tokenize("/");
   if (tokens->IsEmpty()) {
      fListTree->HighlightItem(item);
      DoubleClicked(item, 1);
      delete tokens;
      fListTree->ClearViewPort();
      fListTree->AdjustPosition(item);
      return;
   }
   // if the Browser.ExpandDirectories option is set to "no", then don't
   // expand the parent directory tree (for example on nfs)
   TString str = gEnv->GetValue("Browser.ExpandDirectories", "yes");
   str.ToLower();
   expand = (str == "yes") ? kTRUE : kFALSE;
   TString first = ((TObjString*)tokens->At(0))->GetName();
   // always prevent expanding the parent directory tree on afs
   if (first == "afs")
      expand = kFALSE;
   // check also AFS_SUPER_MAGIC, NFS_SUPER_MAGIC, FUSE_SUPER_MAGIC,
   // CIFS_MAGIC_NUMBER and SMB_SUPER_MAGIC
   if (!gSystem->GetFsInfo(path, (Long_t *)&id, &bsize, &blocks, &bfree))
      if (id == 0x5346414f || id == 0x6969 || id == 0x65735546 || id == 0xff534d42 || id == 0x517b)
         expand = kFALSE;
   if (first.Length() == 2 && first.EndsWith(":")) {
      TList *curvol  = gSystem->GetVolumes("cur");
      if (curvol) {
         TNamed *drive = (TNamed *)curvol->At(0);
         if (first == drive->GetName()) {
            TString infos = drive->GetTitle();
            if (infos.Contains("Network"))
               expand = kFALSE;
         }
         delete curvol;
      }
   }
   for (Int_t i = 0; i < tokens->GetEntriesFast(); ++i) {
      TString token = ((TObjString*)tokens->At(i))->GetName();
      if (token.Length() == 2 && token.EndsWith(":")) {
         token.Append("\\");
         itm = fListTree->FindChildByName(0, token);
         if (itm) {
            item = itm;
            fListTree->OpenItem(item);
            if (expand)
               DoubleClicked(item, 1);
         }
         continue;
      }
      itm = fListTree->FindChildByName(item, token);
      if (itm) {
         item = itm;
         fListTree->OpenItem(item);
         if (expand)
            DoubleClicked(item, 1);
      }
      else {
         itm = fListTree->AddItem(item, token);
         item = itm;
         fListTree->OpenItem(item);
         if (expand)
            DoubleClicked(item, 1);
      }
   }
   fListTree->HighlightItem(item);
   DoubleClicked(item, 1);
   delete tokens;
   fListTree->ClearViewPort();
   fListTree->AdjustPosition(item);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot used to switch to the tab containing the current pad/canvas (gPad)
/// used e.g. when drawing a histogram by double-clicking on its list tree
/// item in a root file.

void TGFileBrowser::PadModified()
{
   if (fDblClick && fNewBrowser) {
      Int_t i;
      TGTab *tabRight = fNewBrowser->GetTabRight();
      for (i=0;i<tabRight->GetNumberOfTabs();++i) {
         TGFrameElement *fe = 0;
         TGCompositeFrame *embed = 0;
         TGCompositeFrame *frame = tabRight->GetTabContainer(i);
         if (frame)
            fe = (TGFrameElement *)frame->GetList()->First();
         if (fe)
            embed = (TGCompositeFrame *)fe->fFrame;
         if (embed && embed->InheritsFrom("TRootCanvas")) {
            ULong_t canvas = gROOT->ProcessLine(TString::Format("((TRootCanvas *)0x%lx)->Canvas();",
                                                (ULong_t)embed));
            if ((canvas) && (canvas == (ULong_t)gPad ||
                canvas == (ULong_t)gPad->GetCanvas())) {
               tabRight->SetTab(i, kTRUE);
               break;
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Open a dialog box asking for a string to be used as filter (regexp), and
/// add an entry in the map of filtered entries. Entering "*" or empty string
/// ("") will disable filtering on the current list tree item.

void TGFileBrowser::RequestFilter()
{
   char filter[1024];
   if (!fListLevel)
      return;
   // initialize with previous (active) filter string
   snprintf(filter, sizeof(filter), "%s", fFilterStr.Data());
   new TGInputDialog(gClient->GetRoot(), this,
                     "Enter filter expression:\n(empty string \"\" or \"*\" to remove filter)",
                      filter, filter);
   // if user pressed cancel, update the status of the current list tree
   // item and return
   if ((filter[0] == 0) && (filter[1] == 0)) {
      CheckFiltered(fListLevel, kTRUE);
      return;
   }
   else if (((filter[0] == 0) && (filter[1] == 1)) || !strcmp(filter, "*")) {
      // if user entered "*" or "", just disable filtering for the current
      // list tree item
      fFilterButton->SetState(kButtonUp);
      fFilteredItems.erase(fListLevel);
   } else {
      // if user entered a string different from "*", use it to create an
      // entry in the filter map
      fFilterStr = filter;
      fFilterButton->SetState(kButtonEngaged);
      // if there is already a filter on this item, delete it
      if (CheckFiltered(fListLevel))
         fFilteredItems.erase(fListLevel);
      // insert a new entry for the current list tree item
      fFilteredItems.emplace(fListLevel, filter);
   }
   // finally update the list tree
   fListTree->DeleteChildren(fListLevel);
   DoubleClicked(fListLevel, 1);
   fListTree->ClearViewPort();
   fListTree->AdjustPosition(fListLevel);
}

////////////////////////////////////////////////////////////////////////////////
/// A ROOT File has been selected in TGHtmlBrowser.

void TGFileBrowser::Selected(char *)
{
   TGListTreeItem *itm = fListTree->FindChildByData(0, gROOT->GetListOfFiles());
   if (itm) {
      fListTree->ClearHighlighted();
      fListLevel = itm;
      fListTree->HighlightItem(fListLevel);
      fListTree->OpenItem(fListLevel);
      BrowseObj(gROOT->GetListOfFiles());
      fListTree->ClearViewPort();
      fListTree->AdjustPosition(fListLevel);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Toggle the sort mode and set the "sort button" state accordingly.

void TGFileBrowser::ToggleSort()
{
   if (!fListLevel) return;
   TString itemname;
   TGListTreeItem *item = fListLevel;
   if (!fListLevel->GetFirstChild()) {
      item = fListLevel->GetParent();
      itemname = fListLevel->GetText();
   }
   if (!item)
      return;
   Bool_t is_sorted = CheckSorted(item);
   if (!is_sorted) {
      //alphabetical sorting
      fListTree->SortChildren(item);
      fSortedItems.push_back(item);
      fSortButton->SetState(kButtonEngaged);
   } else {
      fListTree->DeleteChildren(item);
      DoubleClicked(item, 1);
      fSortedItems.remove(item);
      fSortButton->SetState(kButtonUp);
      gClient->NeedRedraw(fListTree, kTRUE);
      gClient->HandleInput();
      if (itemname.Length() > 0) {
         TGListTreeItem *itm = fListTree->FindChildByName(item, itemname);
         if (itm) {
            fListTree->ClearHighlighted();
            Clicked(itm, 1, 0, 0);
            itm->SetActive(kTRUE);
            fListTree->SetSelected(itm);
            fListTree->HighlightItem(itm, kTRUE, kTRUE);
         }
      }
   }
   fListTree->ClearViewPort();
   fListTree->AdjustPosition(fListLevel);
}


