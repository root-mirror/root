// @(#)root/gui:$Id$
// Author: Fons Rademakers   27/02/98

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRootBrowserLite                                                     //
//                                                                      //
// This class creates a ROOT object browser (looking like Windows       //
// Explorer). The widgets used are the new native ROOT GUI widgets.     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RConfigure.h"

#include "TRootBrowserLite.h"
#include "TRootApplication.h"
#include "TGCanvas.h"
#include "TGMenu.h"
#include "TGFileDialog.h"
#include "TGStatusBar.h"
#include "TGFSComboBox.h"
#include "TGLabel.h"
#include "TGButton.h"
#include "TGListView.h"
#include "TGListTree.h"
#include "TGToolBar.h"
#include "TGSplitter.h"
#include "TG3DLine.h"
#include "TGFSContainer.h"
#include "TGMimeTypes.h"
#include "TRootHelpDialog.h"
#include "TGTextEntry.h"
#include "TGTextEdit.h"
#include "TGTextEditDialogs.h"

#include "TROOT.h"
#include "TEnv.h"
#include "TBrowser.h"
#include "TApplication.h"
#include "TFile.h"
#include "TKey.h"
#include "TKeyMapFile.h"
#include "TClass.h"
#include "TContextMenu.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TRemoteObject.h"
#include "TInterpreter.h"
#include "TGuiBuilder.h"
#include "TImage.h"
#include "TVirtualPad.h"
#include "KeySymbols.h"
#include "THashTable.h"
#include "TColor.h"
#include "TObjString.h"
#include "TGDNDManager.h"
#include "TBufferFile.h"
#include "TFolder.h"
#include "Getline.h"
#include "snprintf.h"
#include "TVirtualX.h"

#include "HelpText.h"

#ifdef WIN32
#include "TWin32SplashThread.h"
#endif

// Browser menu command ids
enum ERootBrowserCommands {
   kFileNewBrowserLite,
   kFileNewBrowser,
   kFileNewCanvas,
   kFileNewBuilder,
   kFileOpen,
   kFileSave,
   kFileSaveAs,
   kFilePrint,
   kFileCloseBrowser,
   kFileQuit,

   kViewToolBar,
   kViewStatusBar,
   kViewLargeIcons,
   kViewSmallIcons,
   kViewList,
   kViewDetails,
   kViewLineUp,
   kViewHidden,
   kViewRefresh,
   kViewFind,
   kViewExec,
   kViewInterrupt,
   kViewSave,

   kViewArrangeByName,     // Arrange submenu
   kViewArrangeByType,
   kViewArrangeBySize,
   kViewArrangeByDate,
   kViewArrangeAuto,
   kViewGroupLV,

   kHistoryBack,
   kHistoryForw,

   kOptionShowCycles,
   kOptionAutoThumbnail,

   kOneLevelUp,            // One level up toolbar button
   kFSComboBox,            // File system combobox in toolbar

   kHelpAbout,
   kHelpOnBrowser,
   kHelpOnCanvas,
   kHelpOnMenus,
   kHelpOnGraphicsEd,
   kHelpOnObjects,
   kHelpOnPS,
   kHelpOnRemote
};


//----- Struct for default icons

struct DefaultIcon_t {
   const char      *fPicnamePrefix;
   const TGPicture *fIcon[2];
};

#if 0
static DefaultIcon_t gDefaultIcon[] = {
   { "folder",  { 0, 0 } },
   { "app",     { 0, 0 } },
   { "doc",     { 0, 0 } },
   { "slink",   { 0, 0 } },
   { "histo",   { 0, 0 } },
   { "object",  { 0, 0 } }
};
#endif


//----- Toolbar stuff...

static ToolBarData_t gToolBarData[] = {
   { "tb_uplevel.xpm",   "Up One Level",   kFALSE, kOneLevelUp, 0 },
   { "",                 "",               kFALSE, -1, 0 },
   { "tb_bigicons.xpm",  "Large Icons",    kTRUE,  kViewLargeIcons, 0 },
   { "tb_smicons.xpm",   "Small Icons",    kTRUE,  kViewSmallIcons, 0 },
   { "tb_list.xpm",      "List",           kTRUE,  kViewList, 0 },
   { "tb_details.xpm",   "Details",        kTRUE,  kViewDetails, 0 },
   { "",                 "",               kFALSE, -1, 0 },
   { "tb_back.xpm",      "Back",           kFALSE, kHistoryBack, 0 },
   { "tb_forw.xpm",      "Forward",        kFALSE, kHistoryForw, 0 },
   { "tb_refresh.xpm",   "Refresh (F5)",   kFALSE, kViewRefresh, 0 },
   { "",                 "",               kFALSE, -1, 0 },
   { "tb_find.xpm",      "Find (Ctrl-F)",  kFALSE, kViewFind, 0 },
   { "",                 "",               kFALSE, -1, 0 },
   { "macro_t.xpm",      "Execute Macro",  kFALSE, kViewExec, 0 },
   { "interrupt.xpm",    "Interrupt Macro",kFALSE, kViewInterrupt, 0 },
   { "filesaveas.xpm",   "Save Macro",     kFALSE, kViewSave, 0 },
   { 0,                  0,                kFALSE, 0, 0 }
};


//----- TGFileDialog file types

static const char *gOpenTypes[] = { "ROOT files",   "*.root",
                                    "All files",    "*",
                                    0,              0 };

////////////////////////////////////////////////////////////////////////////////////
class TRootBrowserHistoryCursor : public TObject {
public:
   TGListTreeItem *fItem;

   TRootBrowserHistoryCursor(TGListTreeItem *item) : fItem(item) {}
   void Print(Option_t *) const {  if (fItem) printf("%s\n", fItem->GetText()); }
};


////////////////////////////////////////////////////////////////////////////////////
class TRootBrowserHistory : public TList {
public:
   void RecursiveRemove(TObject *obj) {
      TRootBrowserHistoryCursor *cur;
      TIter next(this);

      while ((cur = (TRootBrowserHistoryCursor*)next())) {
         if (cur->fItem->GetUserData() == obj) {
            Remove(cur);
            delete cur;
         }
      }
   }

   void DeleteItem(TGListTreeItem *item) {
      TRootBrowserHistoryCursor *cur;
      TIter next(this);

      while ((cur = (TRootBrowserHistoryCursor*)next())) {
         if (cur->fItem == item) {
            Remove(cur);
            delete cur;
         }
      }
   }
};


////////////////////////////////////////////////////////////////////////////////////
class TRootBrowserCursorSwitcher {
private:
   TGWindow *fW1;
   TGWindow *fW2;
public:
   TRootBrowserCursorSwitcher(TGWindow *w1, TGWindow *w2) : fW1(w1), fW2(w2) {
      if (w1) gVirtualX->SetCursor(w1->GetId(), gVirtualX->CreateCursor(kWatch));
      if (w2) gVirtualX->SetCursor(w2->GetId(), gVirtualX->CreateCursor(kWatch));
   }
   ~TRootBrowserCursorSwitcher() {
      if (fW1) gVirtualX->SetCursor(fW1->GetId(), gVirtualX->CreateCursor(kPointer));
      if (fW2) gVirtualX->SetCursor(fW2->GetId(), gVirtualX->CreateCursor(kPointer));
   }
};

////////////////////////////////////////////////////////////////////////////////////
class TIconBoxThumb : public TObject {
public:
   TString fName;
   const TGPicture *fSmall;
   const TGPicture *fLarge;

   TIconBoxThumb(const char *name, const TGPicture *spic, const TGPicture *pic) {
      fName = name;
      fSmall = spic;
      fLarge = pic;
   }
   ULong_t Hash() const { return fName.Hash(); }
   const char *GetName() const { return fName.Data(); }
};



//----- Special ROOT object item (this are items in the icon box, see
//----- TRootIconBox)
////////////////////////////////////////////////////////////////////////////////////
class TRootObjItem : public TGFileItem {
public:
   TRootObjItem(const TGWindow *p, const TGPicture *bpic,
                const TGPicture *spic, TGString *name,
                TObject *obj, TClass *cl, EListViewMode viewMode = kLVSmallIcons);

   virtual TDNDData *GetDNDData(Atom_t) {
      TObject *object = 0;
      if (fObj->IsA() == TKey::Class())
         object = ((TKey *)fObj)->ReadObj();
      else
         object = fObj;
      if (object) {
         if (!fBuf) fBuf = new TBufferFile(TBuffer::kWrite);
         fBuf->WriteObject(object);
         fDNDData.fData = fBuf->Buffer();
         fDNDData.fDataLength = fBuf->Length();
      }
      fDNDData.fDataType = gVirtualX->InternAtom("application/root", kFALSE);
      return &fDNDData;
   }

   virtual Bool_t HandleDNDFinished() {
      if (GetParent())
         return ((TGFrame *)GetParent())->HandleDNDFinished();
      return kFALSE;
   }

protected:
   TObject     *fObj;
   TDNDData     fDNDData;
};

////////////////////////////////////////////////////////////////////////////////
/// Create an icon box item.

TRootObjItem::TRootObjItem(const TGWindow *p, const TGPicture *bpic,
                           const TGPicture *spic, TGString *name,
                           TObject *obj, TClass *, EListViewMode viewMode) :
   TGFileItem(p, bpic, 0, spic, 0, name, 0, 0, 0, 0, 0, viewMode)
{
   fObj = obj;
   fDNDData.fData = 0;
   fDNDData.fDataLength = 0;

   if (fSubnames) {
      for (Int_t i = 0; fSubnames[i] != 0; ++i) delete fSubnames[i];
   }
   delete [] fSubnames;
   fSubnames = new TGString* [2];

   fSubnames[0] = new TGString(obj->GetTitle());

   fSubnames[1] = 0;

   if (obj->IsA()->HasDefaultConstructor()) {
      SetDNDSource(kTRUE);
   }
   if ((obj->IsA() == TFolder::Class()) ||
       (obj->IsA() == TClass::Class())) {
      SetDNDSource(kFALSE);
   }

   int i;
   for (i = 0; fSubnames[i] != 0; ++i)
      ;
   fCtw = new int[i];
   for (i = 0; fSubnames[i] != 0; ++i)
      fCtw[i] = gVirtualX->TextWidth(fFontStruct, fSubnames[i]->GetString(),
                                     fSubnames[i]->GetLength());
}

class TRootIconBox;
////////////////////////////////////////////////////////////////////////////////////
class TRootIconList : public TList {

private:
   TRootIconBox    *fIconBox; // iconbox to which list belongs
   const TGPicture *fPic;     // list view icon

public:
   TRootIconList(TRootIconBox* box = 0);
   virtual ~TRootIconList();
   void              UpdateName();
   const char       *GetTitle() const { return "ListView Container"; }
   Bool_t            IsFolder() const { return kFALSE; }
   void              Browse(TBrowser *b);
   const TGPicture  *GetPicture() const { return fPic; }
};

////////////////////////////////////////////////////////////////////////////////
/// constructor

TRootIconList::TRootIconList(TRootIconBox* box)
{
   fPic = gClient->GetPicture("listview.xpm");
   fIconBox = box;
   fName = "empty";
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TRootIconList::~TRootIconList()
{
   gClient->FreePicture(fPic);
}

////////////////////////////////////////////////////////////////////////////////
/// composite name

void TRootIconList::UpdateName()
{
   if (!First()) return;

   if (fSize==1) {
      fName = First()->GetName();
      return;
   }

   fName = First()->GetName();
   fName += "-";
   fName += Last()->GetName();
}

//----- Special ROOT object container (this is the icon box on the
//----- right side of the browser)
////////////////////////////////////////////////////////////////////////////////////
class TRootIconBox : public TGFileContainer {
friend class TRootIconList;
friend class TRootBrowserLite;

private:
   Bool_t           fCheckHeaders;   // if true check headers
   TRootIconList   *fCurrentList;    //
   TRootObjItem    *fCurrentItem;    //
   Bool_t           fGrouped;        //
   TString          fCachedPicName;  //
   TList           *fGarbage;        // garbage for  TRootIconList's
   Int_t            fGroupSize;      // the total number of items when icon box switched to "global view" mode
   TGString        *fCurrentName;    //
   const TGPicture *fLargeCachedPic; //
   const TGPicture *fSmallCachedPic; //
   Bool_t           fWasGrouped;
   TObject         *fActiveObject;   //
   Bool_t           fIsEmpty;
   THashTable      *fThumbnails;     // hash table with thumbnailed pictures
   Bool_t           fAutoThumbnail;  //
   TRootBrowserLite    *fBrowser;

   void  *FindItem(const TString& name,
                   Bool_t direction = kTRUE,
                   Bool_t caseSensitive = kTRUE,
                   Bool_t beginWith = kFALSE);
   void RemoveGarbage();

public:
   TRootIconBox(TRootBrowserLite *browser, TGListView *lv,
                UInt_t options = kSunkenFrame,
                ULong_t back = GetDefaultFrameBackground());

   virtual ~TRootIconBox();

   void   AddObjItem(const char *name, TObject *obj, TClass *cl);
   void   GetObjPictures(const TGPicture **pic, const TGPicture **spic,
                         TObject *obj, const char *name);
   void   SetObjHeaders();
   void   Refresh();
   void   RemoveAll();
   void   SetGroupSize(Int_t siz) { fGroupSize = siz; }
   Int_t  GetGroupSize() const { return fGroupSize; }
   TGFrameElement *FindFrame(Int_t x, Int_t y, Bool_t exclude=kTRUE) { return TGContainer::FindFrame(x,y,exclude); }
   Bool_t WasGrouped() const { return fWasGrouped; }
};

////////////////////////////////////////////////////////////////////////////////
/// Create iconbox containing ROOT objects in browser.

TRootIconBox::TRootIconBox(TRootBrowserLite *browser, TGListView *lv, UInt_t options,
                           ULong_t back) : TGFileContainer(lv, options, back)
{
   fListView = lv;
   fBrowser = browser;

   fCheckHeaders = kTRUE;
   fTotal = 0;
   fGarbage = new TList();
   fCurrentList = 0;
   fCurrentItem = 0;
   fGrouped = kFALSE;
   fGroupSize = 1000;
   fCurrentName = 0;
   fWasGrouped = kFALSE;
   fActiveObject = 0;
   fIsEmpty = kTRUE;
   fLargeCachedPic = 0;
   fSmallCachedPic = 0;

   // Don't use timer HERE (timer is set in TBrowser).
   StopRefreshTimer();
   fRefresh = 0;
   fThumbnails = new THashTable(50);
   fAutoThumbnail = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TRootIconBox::~TRootIconBox()
{
   RemoveAll();
   RemoveGarbage();
   delete fGarbage;
   delete fThumbnails;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve icons associated with class "name". Association is made
/// via the user's ~/.root.mimes file or via $ROOTSYS/etc/root.mimes.

void TRootIconBox::GetObjPictures(const TGPicture **pic, const TGPicture **spic,
                                  TObject *obj, const char *name)
{
   static TImage *im = 0;
   if (!im) {
      im = TImage::Create();
   }

   TString xpm_magic(name, 3);
   Bool_t xpm = xpm_magic == "/* ";
   const char *iconname = xpm ? obj->GetName() : name;

   if (obj->IsA()->InheritsFrom("TGeoVolume")) {
      iconname = obj->GetIconName() ? obj->GetIconName() : obj->IsA()->GetName();
   }

   if (fCachedPicName == iconname) {
      *pic = fLargeCachedPic;
      *spic = fSmallCachedPic;
      return;
   }

   *pic = fClient->GetMimeTypeList()->GetIcon(iconname, kFALSE);

   if (!(*pic) && xpm) {
      if (im && im->SetImageBuffer((char**)&name, TImage::kXpm)) {
         *pic = fClient->GetPicturePool()->GetPicture(iconname, im->GetPixmap(),
                                                      im->GetMask());
         im->Scale(im->GetWidth()/2, im->GetHeight()/2);
         *spic = fClient->GetPicturePool()->GetPicture(iconname, im->GetPixmap(),
                                                      im->GetMask());
      }

      fClient->GetMimeTypeList()->AddType("[thumbnail]", iconname, iconname, iconname, "->Browse()");
      return;
   }

   if (*pic == 0) {
      if (obj->IsFolder()) {
         *pic = fFolder_s;
      } else {
         *pic = fDoc_s;
      }
   }
   fLargeCachedPic = *pic;

   *spic = fClient->GetMimeTypeList()->GetIcon(iconname, kTRUE);

   if (*spic == 0) {
      if (obj->IsFolder())
         *spic = fFolder_t;
      else
         *spic = fDoc_t;
   }
   fSmallCachedPic = *spic;
   fCachedPicName = iconname;
}

////////////////////////////////////////////////////////////////////////////////
/// delete all TRootIconLists from garbage

void TRootIconBox::RemoveGarbage()
{
   TIter next(fGarbage);
   TList *li;

   while ((li=(TList *)next())) {
      li->Clear("nodelete");
   }
   fGarbage->Delete();
}

////////////////////////////////////////////////////////////////////////////////
/// Add object to iconbox. Class is used to get the associated icons
/// via the mime file (see GetObjPictures()).

void TRootIconBox::AddObjItem(const char *name, TObject *obj, TClass *cl)
{
   if (!cl) return;

   Bool_t isSystemFile = kFALSE;
   TGFileItem *fi;
   fWasGrouped = kFALSE;
   const TGPicture *pic = 0;
   const TGPicture *spic = 0;

   if (obj->InheritsFrom("TRemoteObject")) {
      // check if the real remote object is a system file or directory
      TRemoteObject *robj = (TRemoteObject *)obj;
      if ((TString(robj->GetClassName()) == "TSystemFile") ||
          (TString(robj->GetClassName()) == "TSystemDirectory"))
         isSystemFile = kTRUE;
   }

   if (isSystemFile || obj->IsA() == TSystemFile::Class() ||
       obj->IsA() == TSystemDirectory::Class()) {
      if (fCheckHeaders) {
         if (strcmp(fListView->GetHeader(1), "Attributes")) {
            fListView->SetDefaultHeaders();
            TGTextButton** buttons = fListView->GetHeaderButtons();
            if (buttons) {
               buttons[0]->Connect("Clicked()", "TRootBrowserLite", fBrowser,
                                   TString::Format("SetSortMode(=%d)", kViewArrangeByName));
               buttons[1]->Connect("Clicked()", "TRootBrowserLite", fBrowser,
                                   TString::Format("SetSortMode(=%d)", kViewArrangeByType));
               buttons[2]->Connect("Clicked()", "TRootBrowserLite", fBrowser,
                                   TString::Format("SetSortMode(=%d)", kViewArrangeBySize));
               buttons[5]->Connect("Clicked()", "TRootBrowserLite", fBrowser,
                                   TString::Format("SetSortMode(=%d)", kViewArrangeByDate));
            }
         }
         fCheckHeaders = kFALSE;
      }

      TIconBoxThumb *thumb = 0;
      char *thumbname = gSystem->ConcatFileName(gSystem->WorkingDirectory(), name);
      thumb = (TIconBoxThumb *)fThumbnails->FindObject(gSystem->IsAbsoluteFileName(name) ? name :
                                                       thumbname);
      delete []thumbname;

      if (thumb) {
         spic = thumb->fSmall;
         pic =  thumb->fLarge;
      }

      if (obj->InheritsFrom("TRemoteObject"))
         // special case for remote object
         fi = AddRemoteFile(obj, spic, pic);
      else
         fi = AddFile(name, spic, pic);
      if (fi) {
         fi->SetUserData(obj);
         if (obj->IsA() == TSystemFile::Class()) {
            TString str;
            TDNDData data;
            str = TString::Format("file://%s/%s\r\n",
                                  gSystem->UnixPathName(obj->GetTitle()),
                                  obj->GetName());
            data.fData = (void *)str.Data();
            data.fDataLength = str.Length()+1;
            data.fDataType = gVirtualX->InternAtom("text/uri-list", kFALSE);
            fi->SetDNDData(&data);
            fi->SetDNDSource(kTRUE);
         }
      }

      fIsEmpty = kFALSE;
      return;
   }

   if (!fCurrentList) {
      fCurrentList = new TRootIconList(this);
      fGarbage->Add(fCurrentList);
   }

   fCurrentList->Add(obj);
   fCurrentList->UpdateName();
   fIsEmpty = kFALSE;

   TGFrameElement *el;
   TIter next(fList);
   while ((el = (TGFrameElement *) next())) {
      TGLVEntry *f = (TGLVEntry *) el->fFrame;
      if (f->GetUserData() == obj) {
         return;
      }
   }

   if (fGrouped && fCurrentItem && (fCurrentList->GetSize()>1)) {
      fCurrentName->SetString(fCurrentList->GetName());
   }

   EListViewMode view = fListView->GetViewMode();

   if ((fCurrentList->GetSize() < fGroupSize) && !fGrouped) {
      GetObjPictures(&pic, &spic, obj, obj->GetIconName() ?
                     obj->GetIconName() : cl->GetName());

      if (fCheckHeaders) {
         if (strcmp(fListView->GetHeader(1), "Title")) {
            SetObjHeaders();
         }
         fCheckHeaders = kFALSE;
      }

      fi = new TRootObjItem(this, pic, spic, new TGString(name), obj, cl, view);

      fi->SetUserData(obj);
      AddItem(fi);
      return;
   }

   if (fGrouped && (fCurrentList->GetSize()==1)) {
      fCurrentName = new TGString(fCurrentList->GetName());
      fCurrentItem = new TRootObjItem(this, fCurrentList->GetPicture(), fCurrentList->GetPicture(),
                                      fCurrentName,fCurrentList, TList::Class(), view);
      fCurrentItem->SetUserData(fCurrentList);
      AddItem(fCurrentItem);
      fTotal = fList->GetSize();
      return;
   }

   if ((fCurrentList->GetSize()==fGroupSize) && !fGrouped) {
      fGrouped = kTRUE;

      // clear fList
      TGFrameElement *el2;
      TIter nextl(fList);

      while ((el2 = (TGFrameElement *) nextl())) {
         el2->fFrame->DestroyWindow();
         delete el2->fFrame;
         fList->Remove(el2);
         delete el2;
      }

      fCurrentName = new TGString(fCurrentList->GetName());
      fi = new TRootObjItem(this, fCurrentList->GetPicture(), fCurrentList->GetPicture(),
                            fCurrentName, fCurrentList, TList::Class(), view);
      fi->SetUserData(fCurrentList);
      AddItem(fi);

      fCurrentList = new TRootIconList(this);
      fGarbage->Add(fCurrentList);
      fTotal = 1;
      return;
   }

   if ((fCurrentList->GetSize()==fGroupSize) && fGrouped) {
      fCurrentList = new TRootIconList(this);
      fGarbage->Add(fCurrentList);
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// browse icon list

void TRootIconList::Browse(TBrowser *)
{
   if (!fIconBox) return;

   TObject *obj;
   TGFileItem *fi;
   const TGPicture *pic = 0;
   const TGPicture *spic = 0;
   TClass *cl;
   TString name;
   TKey *key = 0;

   fIconBox->RemoveAll();
   TObjLink *lnk = FirstLink();

   while (lnk) {
      obj = lnk->GetObject();
      lnk = lnk->Next();

      if (obj->IsA() == TKey::Class()) {
         cl = TClass::GetClass(((TKey *)obj)->GetClassName());
         key = (TKey *)obj;
      } else if (obj->IsA() == TKeyMapFile::Class()) {
         cl = TClass::GetClass(((TKeyMapFile *)obj)->GetTitle());
      } else if (obj->InheritsFrom("TRemoteObject")) {
         // special case for remote object: get real object class
         TRemoteObject *robj = (TRemoteObject *)obj;
         cl = TClass::GetClass(robj->GetClassName());
      } else {
         cl = obj->IsA();
      }

      name = obj->GetName();

      if (key && obj->IsA() == TKey::Class()) {
         name += ";";
         name +=  key->GetCycle();
      }

      fIconBox->GetObjPictures(&pic, &spic, obj, obj->GetIconName() ?
                               obj->GetIconName() : cl->GetName());

      fi = new TRootObjItem((const TGWindow*)fIconBox, pic, spic, new TGString(name.Data()),
                             obj, cl, (EListViewMode)fIconBox->GetViewMode());
      fi->SetUserData(obj);
      fIconBox->AddItem(fi);
      fIconBox->fTotal++;

      if (obj==fIconBox->fActiveObject) {
         fIconBox->ActivateItem((TGFrameElement*)fIconBox->fList->Last());
      }
   }

   fIconBox->fGarbage->Remove(this);
   fIconBox->RemoveGarbage();
   fIconBox->fGarbage->Add(this); // delete this later

   fIconBox->Refresh();
   fIconBox->AdjustPosition();

   fIconBox->fWasGrouped = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Find a frame which assosiated object has a name containing a "name" string.

void *TRootIconBox::FindItem(const TString& name, Bool_t direction,
                             Bool_t caseSensitive,Bool_t beginWith)
{
   if (!fGrouped) {
      return TGContainer::FindItem(name, direction, caseSensitive, beginWith);
   }

   if (name.IsNull()) return 0;
   int idx = kNPOS;

   TGFrameElement* el = 0;
   TString str;
   TString::ECaseCompare cmp = caseSensitive ? TString::kExact : TString::kIgnoreCase;

   fLastDir = direction;
   fLastCase = caseSensitive;
   fLastName = name;

   if (fLastActiveEl) {
      el = fLastActiveEl;

      if (direction) {
         el = (TGFrameElement *)fList->After(el);
      } else {
         el = (TGFrameElement *)fList->Before(el);
      }
   } else {
      if (direction) el = (TGFrameElement *)fList->First();
      else el  = (TGFrameElement *)fList->Last();
   }

   TGLVEntry* lv = 0;
   TObject* obj = 0;
   TList* li = 0;

   while (el) {
      lv = (TGLVEntry*)el->fFrame;
      li = (TList*)lv->GetUserData();

      TIter next(li);

      while ((obj=next())) {
         str = obj->GetName();
         idx = str.Index(name,0,cmp);

         if (idx!=kNPOS) {
            if (beginWith) {
               if (idx==0) {
                  fActiveObject = obj;
                  return el;
               }
            } else {
               fActiveObject = obj;
               return el;
            }
         }
      }
      if (direction) {
         el = (TGFrameElement *)fList->After(el);
      } else {
         el = (TGFrameElement *)fList->Before(el);
      }
   }
   fActiveObject = 0;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set list box headers used to display detailed object iformation.
/// Currently this is only "Name" and "Title".

void TRootIconBox::SetObjHeaders()
{
   fListView->SetHeaders(2);
   fListView->SetHeader("Name",  kTextLeft, kTextLeft, 0);
   fListView->SetHeader("Title", kTextLeft, kTextLeft, 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Sort icons, and send message to browser with number of objects
/// in box.

void TRootIconBox::Refresh()
{
   // This automatically calls layout
   Sort(fSortType);

   // Make TRootBrowserLite display total objects in status bar
   SendMessage(fMsgWindow, MK_MSG(kC_CONTAINER, kCT_SELCHANGED), fTotal, fSelected);

   MapSubwindows();
   fListView->AdjustHeaders();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all items from icon box

void TRootIconBox::RemoveAll()
{
   if (fIsEmpty) return;

   fCheckHeaders = kTRUE;
   TGFileContainer::RemoveAll();
   fGrouped = kFALSE;
   fCurrentItem = 0;
   fCurrentList = 0;
   fIsEmpty = kTRUE;
}


//_____________________________________________________________________________
//
// TRootBrowserLite
//
// ROOT object browser (looking like Windows Explorer).
//_____________________________________________________________________________

ClassImp(TRootBrowserLite);

////////////////////////////////////////////////////////////////////////////////
/// Create browser with a specified width and height.

TRootBrowserLite::TRootBrowserLite(TBrowser *b, const char *name, UInt_t width, UInt_t height)
   : TGMainFrame(gClient->GetDefaultRoot(), width, height), TBrowserImp(b)
{
   CreateBrowser(name);

   Resize(width, height);
   if (b) Show();
}

////////////////////////////////////////////////////////////////////////////////
/// Create browser with a specified width and height and at position x, y.

TRootBrowserLite::TRootBrowserLite(TBrowser *b, const char *name, Int_t x, Int_t y,
                           UInt_t width, UInt_t height)
   : TGMainFrame(gClient->GetDefaultRoot(), width, height), TBrowserImp(b)
{
   CreateBrowser(name);

   MoveResize(x, y, width, height);
   SetWMPosition(x, y);
   if (b) Show();
}

////////////////////////////////////////////////////////////////////////////////
/// Browser destructor.

TRootBrowserLite::~TRootBrowserLite()
{
   if (fIconPic) gClient->FreePicture(fIconPic);

   delete fToolBarSep;

   fToolBar->Cleanup();
   delete fToolBar;
   delete fStatusBar;
   delete fV1;
   delete fV2;
   delete fLbl1;
   delete fLbl2;
   delete fHf;
   delete fTreeHdr;
   delete fListHdr;
   delete fIconBox;
   delete fListView;
   delete fLt;
   delete fTreeView;

   delete fMenuBar;
   delete fFileMenu;
   delete fViewMenu;
   delete fOptionMenu;
   delete fHelpMenu;
   delete fSortMenu;

   delete fMenuBarLayout;
   delete fMenuBarItemLayout;
   delete fMenuBarHelpLayout;
   delete fBarLayout;

   delete fTextEdit;

   if (fWidgets) fWidgets->Delete();
   delete fWidgets;

   fHistory->Delete();
   delete fHistory;
}

////////////////////////////////////////////////////////////////////////////////
/// Create the actual browser.

void TRootBrowserLite::CreateBrowser(const char *name)
{
   fWidgets = new TList;
   fHistory = new TRootBrowserHistory;
   fHistoryCursor = 0;
   fBrowseTextFile = kFALSE;

   // Create menus
   fFileMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fFileMenu->AddEntry("&New Browser",        kFileNewBrowser);
   fFileMenu->AddEntry("New Browser &Lite",   kFileNewBrowserLite);
   fFileMenu->AddEntry("New Canvas",          kFileNewCanvas);
   fFileMenu->AddEntry("&Gui Builder",        kFileNewBuilder);
   fFileMenu->AddEntry("&Open...",            kFileOpen);
   fFileMenu->AddSeparator();
   fFileMenu->AddEntry("&Save",               kFileSave);
   fFileMenu->AddEntry("Save &As...",         kFileSaveAs);
   fFileMenu->AddSeparator();
   fFileMenu->AddEntry("&Print...",           kFilePrint);
   fFileMenu->AddSeparator();
   fFileMenu->AddEntry("&Close Browser",      kFileCloseBrowser);
   fFileMenu->AddSeparator();
   fFileMenu->AddEntry("&Quit ROOT",          kFileQuit);

   //fFileMenu->DefaultEntry(kFileNewCanvas);
   fFileMenu->DisableEntry(kFileSave);
   fFileMenu->DisableEntry(kFileSaveAs);
   fFileMenu->DisableEntry(kFilePrint);

   fSortMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fSortMenu->AddEntry("By &Name",            kViewArrangeByName);
   fSortMenu->AddEntry("By &Type",            kViewArrangeByType);
   fSortMenu->AddEntry("By &Size",            kViewArrangeBySize);
   fSortMenu->AddEntry("By &Date",            kViewArrangeByDate);
   fSortMenu->AddSeparator();
   fSortMenu->AddEntry("&Auto Arrange",       kViewArrangeAuto);

   fSortMenu->CheckEntry(kViewArrangeAuto);

   fViewMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fViewMenu->AddEntry("&Toolbar",            kViewToolBar);
   fViewMenu->AddEntry("Status &Bar",         kViewStatusBar);
   fViewMenu->AddSeparator();
   fViewMenu->AddEntry("Lar&ge Icons",        kViewLargeIcons);
   fViewMenu->AddEntry("S&mall Icons",        kViewSmallIcons);
   fViewMenu->AddEntry("&List",               kViewList);
   fViewMenu->AddEntry("&Details",            kViewDetails);
   fViewMenu->AddSeparator();
   fViewMenu->AddEntry("Show &Hidden",        kViewHidden);
   fViewMenu->AddPopup("Arrange &Icons",      fSortMenu);
   fViewMenu->AddEntry("Lin&e up Icons",      kViewLineUp);
   fViewMenu->AddEntry("&Group Icons",        kViewGroupLV);

   fViewMenu->AddSeparator();
   fViewMenu->AddEntry("&Refresh (F5)",       kViewRefresh);

   fViewMenu->CheckEntry(kViewToolBar);
   fViewMenu->CheckEntry(kViewStatusBar);

   if (fBrowser) {
      if (gEnv->GetValue("Browser.ShowHidden", 0)) {
         fViewMenu->CheckEntry(kViewHidden);
         fBrowser->SetBit(TBrowser::kNoHidden, kFALSE);
      } else {
         fViewMenu->UnCheckEntry(kViewHidden);
         fBrowser->SetBit(TBrowser::kNoHidden, kTRUE);
      }
   }

   fOptionMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fOptionMenu->AddEntry("&Show Cycles",        kOptionShowCycles);
   fOptionMenu->AddEntry("&AutoThumbnail",      kOptionAutoThumbnail);

   fHelpMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fHelpMenu->AddEntry("&About ROOT...",        kHelpAbout);
   fHelpMenu->AddSeparator();
   fHelpMenu->AddEntry("Help On Browser...",    kHelpOnBrowser);
   fHelpMenu->AddEntry("Help On Canvas...",     kHelpOnCanvas);
   fHelpMenu->AddEntry("Help On Menus...",      kHelpOnMenus);
   fHelpMenu->AddEntry("Help On Graphics Editor...", kHelpOnGraphicsEd);
   fHelpMenu->AddEntry("Help On Objects...",    kHelpOnObjects);
   fHelpMenu->AddEntry("Help On PostScript...", kHelpOnPS);
   fHelpMenu->AddEntry("Help On Remote Session...", kHelpOnRemote);

   // This main frame will process the menu commands
   fFileMenu->Associate(this);
   fViewMenu->Associate(this);
   fSortMenu->Associate(this);
   fOptionMenu->Associate(this);
   fHelpMenu->Associate(this);

   // Create menubar layout hints
   fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 1, 1);
   fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
   fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

   // Create menubar
   fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);
   fMenuBar->AddPopup("&File",    fFileMenu,    fMenuBarItemLayout);
   fMenuBar->AddPopup("&View",    fViewMenu,    fMenuBarItemLayout);
   fMenuBar->AddPopup("&Options", fOptionMenu,  fMenuBarItemLayout);
   fMenuBar->AddPopup("&Help",    fHelpMenu,    fMenuBarHelpLayout);

   AddFrame(fMenuBar, fMenuBarLayout);

   // Create toolbar and separator

   fToolBarSep = new TGHorizontal3DLine(this);
   fToolBar = new TGToolBar(this, 60, 20, kHorizontalFrame);
   fFSComboBox = new TGFSComboBox(fToolBar, kFSComboBox);

   fComboLayout = new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 0, 0, 2, 2);
   fToolBar->AddFrame(fFSComboBox, fComboLayout);
   fFSComboBox->Resize(190, fFSComboBox->GetDefaultHeight());
   fFSComboBox->Associate(this);

   int spacing = 8;

   for (int i = 0; gToolBarData[i].fPixmap; i++) {
      if (strlen(gToolBarData[i].fPixmap) == 0) {
         spacing = 8;
         continue;
      }
      fToolBar->AddButton(this, &gToolBarData[i], spacing);
      spacing = 0;
   }

   fDrawOption = new TGComboBox(fToolBar, "");
   TGTextEntry *dropt_entry = fDrawOption->GetTextEntry();
   dropt_entry->SetToolTipText("Object Draw Option", 300);
   fDrawOption->Resize(80, 10);
   TGListBox *lb = fDrawOption->GetListBox();
   lb->Resize(lb->GetWidth(), 120);
   Int_t dropt = 1;
   fDrawOption->AddEntry("", dropt++);
   fDrawOption->AddEntry("same", dropt++);
   fDrawOption->AddEntry("box", dropt++);
   fDrawOption->AddEntry("lego", dropt++);
   fDrawOption->AddEntry("colz", dropt++);
   fDrawOption->AddEntry("alp", dropt++);
   fDrawOption->AddEntry("text", dropt++);

   fToolBar->AddFrame(fDrawOption, new TGLayoutHints(kLHintsCenterY | kLHintsRight | kLHintsExpandY,2,2,2,0));
   fToolBar->AddFrame(new TGLabel(fToolBar,"Option"),
                      new TGLayoutHints(kLHintsCenterY | kLHintsRight, 2,2,2,0));

   fBarLayout = new TGLayoutHints(kLHintsTop | kLHintsExpandX);
   AddFrame(fToolBarSep, fBarLayout);
   AddFrame(fToolBar, fBarLayout);

   // Create panes

   fHf = new TGHorizontalFrame(this, 10, 10);

   fV1 = new TGVerticalFrame(fHf, 10, 10, kFixedWidth);
   fV2 = new TGVerticalFrame(fHf, 10, 10);
   fTreeHdr = new TGCompositeFrame(fV1, 10, 10, kSunkenFrame);
   fListHdr = new TGCompositeFrame(fV2, 10, 10, kSunkenFrame);

   fLbl1 = new TGLabel(fTreeHdr, "All Folders");
   fLbl2 = new TGLabel(fListHdr, "Contents of \".\"");

   TGLayoutHints *lo;

   lo = new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 0, 0, 0);
   fWidgets->Add(lo);
   fTreeHdr->AddFrame(fLbl1, lo);
   fListHdr->AddFrame(fLbl2, lo);

   lo = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 1, 2);
   fWidgets->Add(lo);
   fV1->AddFrame(fTreeHdr, lo);
   fV2->AddFrame(fListHdr, lo);

   fV1->Resize(fTreeHdr->GetDefaultWidth()+100, fV1->GetDefaultHeight());

   lo = new TGLayoutHints(kLHintsLeft | kLHintsExpandY);
   fWidgets->Add(lo);
   fHf->AddFrame(fV1, lo);

   TGVSplitter *splitter = new TGVSplitter(fHf);
   splitter->SetFrame(fV1, kTRUE);
   lo = new TGLayoutHints(kLHintsLeft | kLHintsExpandY);
   fWidgets->Add(splitter);
   fWidgets->Add(lo);
   fHf->AddFrame(splitter, lo);

   lo = new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY);
   fWidgets->Add(lo);
   fHf->AddFrame(fV2, lo);

   // Create tree
   fTreeView = new TGCanvas(fV1, 10, 10, kSunkenFrame | kDoubleBorder); // canvas
   fLt = new TGListTree(fTreeView, kHorizontalFrame,fgWhitePixel); // container
   fLt->Associate(this);
   fLt->SetAutoTips();

   fExpandLayout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
   fWidgets->Add(fExpandLayout);
   fV1->AddFrame(fTreeView, fExpandLayout);

   // Create list view (icon box)
   fListView = new TGListView(fV2, 520, 250); // canvas
   // container
   fIconBox = new TRootIconBox(this, fListView, kHorizontalFrame, fgWhitePixel);
   fIconBox->Associate(this);
   fListView->SetIncrements(1, 19); // set vertical scroll one line height at a time
   fViewMode = fListView->GetViewMode();

   TString str = gEnv->GetValue("Browser.AutoThumbnail", "yes");
   str.ToLower();
   fIconBox->fAutoThumbnail = (str == "yes") || atoi(str.Data());
   fIconBox->fAutoThumbnail ? fOptionMenu->CheckEntry(kOptionAutoThumbnail) :
                              fOptionMenu->UnCheckEntry(kOptionAutoThumbnail);

   str = gEnv->GetValue("Browser.GroupView", "10000");
   Int_t igv = atoi(str.Data());

   if (igv>10) {
      fViewMenu->CheckEntry(kViewGroupLV);
      fIconBox->SetGroupSize(igv);
   }

   // reuse lo from "create tree"
   fV2->AddFrame(fListView, fExpandLayout);

   AddFrame(fHf, lo);

   // Statusbar

   int parts[] = { 26, 74 };
   fStatusBar = new TGStatusBar(this, 60, 10);
   fStatusBar->SetParts(parts, 2);
   lo = new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 0, 0, 3, 0);
   AddFrame(fStatusBar, lo);

   fTextEdit = 0;

   // Misc
   TString bname(name);
   bname.Prepend("Old ");
   SetWindowName(bname.Data());
   SetIconName(bname.Data());
   fIconPic = SetIconPixmap("rootdb_s.xpm");
   SetClassHints("ROOT", "Browser");

   SetWMSizeHints(600, 350, 10000, 10000, 2, 2);

   fListLevel = 0;
   fTreeLock  = kFALSE;

   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Escape), 0, kTRUE);
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_F5), 0, kTRUE);
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Right), kKeyMod1Mask, kTRUE);
   gVirtualX->GrabKey(fId, gVirtualX->KeysymToKeycode(kKey_Left), kKeyMod1Mask, kTRUE);
   ClearHistory();
   SetEditDisabled(kEditDisable);

   gVirtualX->SetDNDAware(fId, fDNDTypeList);
   MapSubwindows();
   SetDefaults();
   Resize();
   ShowMacroButtons(kFALSE);

   printf("\n You are using the old ROOT browser! A new version is available. To use it:\n");
   printf(" Select the \"New Browser\" entry from the \"File\" menu in the browser, or change\n");
   printf(" \"Browser.Name:\" from \"TRootBrowserLite\" to \"TRootBrowser\" in system.rootrc\n\n");

   Connect(fLt, "Checked(TObject*, Bool_t)", "TRootBrowserLite",
           this, "Checked(TObject *,Bool_t)");
}

////////////////////////////////////////////////////////////////////////////////
/// handle keys

Bool_t TRootBrowserLite::HandleKey(Event_t *event)
{
   if (event->fType == kGKeyPress) {
      UInt_t keysym;
      char input[10];
      gVirtualX->LookupString(event, input, sizeof(input), keysym);

      if (!event->fState && (EKeySym)keysym == kKey_F5) {
         Refresh(kTRUE);
         return kTRUE;
      }
      if (!event->fState && (EKeySym)keysym == kKey_Escape) {
         if (gDNDManager->IsDragging()) gDNDManager->EndDrag();
      }

      if (event->fState & kKeyMod1Mask) {
         switch ((EKeySym)keysym & ~0x20) {
            case kKey_Right:
               HistoryForward();
               return kTRUE;
            case kKey_Left:
               HistoryBackward();
               return kTRUE;
            default:
               break;
         }
      }
   }
   return TGMainFrame::HandleKey(event);
}

////////////////////////////////////////////////////////////////////////////////
/// Add items to the browser. This function has to be called
/// by the Browse() member function of objects when they are
/// called by a browser. If check < 0 (default) no check box is drawn,
/// if 0 then unchecked checkbox is added, if 1 checked checkbox is added.

void TRootBrowserLite::Add(TObject *obj, const char *name, Int_t check)
{
   if (!obj)
      return;
   if (obj->InheritsFrom("TObjectSpy"))
      return;
   if (!name) name = obj->GetName();

   AddToBox(obj, name);
   if (check > -1) {
      TGFrameElement *el;
      TIter next(fIconBox->fList);
      if (!obj->IsFolder()) {
         while ((el = (TGFrameElement *) next())) {
            TGLVEntry *f = (TGLVEntry *) el->fFrame;
            if (f->GetUserData() == obj) {
               f->SetCheckedEntry(check);
            }
         }
      }
   }

   // Don't show current dir and up dir links in the tree
   if (name[0] == '.' && ((name[1] == '\0') || (name[1] == '.' && name[2] == '\0')))
      return;

   if (obj->IsFolder())
      AddToTree(obj, name, check);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a checkbox in the TGListTreeItem corresponding to obj
/// and a checkmark on TGLVEntry if check = kTRUE.

void TRootBrowserLite::AddCheckBox(TObject *obj, Bool_t check)
{
   if (obj) {
      TGListTreeItem *item = fLt->FindItemByObj(fLt->GetFirstItem(), obj);
      while (item) {
         fLt->SetCheckBox(item, kTRUE);
         fLt->CheckItem(item, check);
         item = fLt->FindItemByObj(item->GetNextSibling(), obj);
      }
      TGFrameElement *el;
      TIter next(fIconBox->fList);
      while ((el = (TGFrameElement *) next())) {
         TGLVEntry *f = (TGLVEntry *) el->fFrame;
         if (f->GetUserData() == obj) {
            f->SetCheckedEntry(check);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check / uncheck the TGListTreeItem corresponding to this
/// object and add a checkmark on TGLVEntry if check = kTRUE.

void TRootBrowserLite::CheckObjectItem(TObject *obj, Bool_t check)
{
   if (obj) {
      TGListTreeItem *item = fLt->FindItemByObj(fLt->GetFirstItem(), obj);
      while (item) {
         fLt->CheckItem(item, check);
         item = fLt->FindItemByObj(item->GetNextSibling(), obj);
         TGFrameElement *el;
         TIter next(fIconBox->fList);
         if (!obj->IsFolder()) {
            while ((el = (TGFrameElement *) next())) {
               TGLVEntry *f = (TGLVEntry *) el->fFrame;
               if (f->GetUserData() == obj) {
                  f->SetCheckedEntry(check);
               }
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove checkbox from TGListTree and checkmark from TGListView.

void TRootBrowserLite::RemoveCheckBox(TObject *obj)
{
   if (obj) {
      TGListTreeItem *item = fLt->FindItemByObj(fLt->GetFirstItem(), obj);
      while (item) {
         fLt->SetCheckBox(item, kFALSE);
         item = fLt->FindItemByObj(item->GetNextSibling(), obj);
         TGFrameElement *el;
         TIter next(fIconBox->fList);
         if (!obj->IsFolder()) {
            while ((el = (TGFrameElement *) next())) {
               TGLVEntry *f = (TGLVEntry *) el->fFrame;
               if (f->GetUserData() == obj) {
                  f->SetCheckedEntry(kFALSE);
               }
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add items to the iconbox of the browser.

void TRootBrowserLite::AddToBox(TObject *obj, const char *name)
{
   if (obj) {
      if (!name) name = obj->GetName() ? obj->GetName() : "NoName";
      //const char *titlePtr = obj->GetTitle() ? obj->GetTitle() : " ";

      TClass *objClass = 0;

      if (obj->IsA() == TKey::Class())
         objClass = TClass::GetClass(((TKey *)obj)->GetClassName());
      else if (obj->IsA() == TKeyMapFile::Class())
         objClass = TClass::GetClass(((TKeyMapFile *)obj)->GetTitle());
      else if (obj->InheritsFrom("TRemoteObject")) {
         // special case for remote object: get real object class
         TRemoteObject *robj = (TRemoteObject *)obj;
         if (!strcmp(robj->GetClassName(), "TKey"))
            objClass = TClass::GetClass(robj->GetKeyClassName());
         else
            objClass = TClass::GetClass(robj->GetClassName());
      }
      else
         objClass = obj->IsA();

      fIconBox->AddObjItem(name, obj, objClass);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add items to the current TGListTree of the browser.

void TRootBrowserLite::AddToTree(TObject *obj, const char *name, Int_t check)
{
   if (!obj)
      return;
   if (obj->InheritsFrom("TApplication"))
      fListLevel = 0;
   if (!fTreeLock) {
      if (!name) name = obj->GetName();
      if (name[0] == '.' && name[1] == '.')
         Info("AddToTree", "up one level %s", name);
      if(check > -1) {
         TGListTreeItem *item = fLt->AddItem(fListLevel, name, obj, 0, 0, kTRUE);
         if (item) fLt->CheckItem(item, (Bool_t)check);
         TString tip(obj->ClassName());
         if (obj->GetTitle()) {
            tip += " ";
            tip += obj->GetTitle();
         }
         fLt->SetToolTipItem(item, tip.Data());
      } else {
         // special case for remote object
         if (obj->InheritsFrom("TRemoteObject")) {
            // Nothing to do
         } else if (fListLevel) {
            // check also if one of its parents is a remote object
            TGListTreeItem *top = fListLevel;
            while (top->GetParent()) {
               TObject *tobj = (TObject *) top->GetUserData();
               if (tobj && (tobj->InheritsFrom("TRemoteObject") ||
                  tobj->InheritsFrom("TApplicationRemote"))) {
                  break;
               }
               top = top->GetParent();
            }
         }
         // add the object only if not already in the list
         if ((!fLt->FindChildByName(fListLevel, name)) &&
             (!fLt->FindChildByData(fListLevel, obj))) {
            TGListTreeItem *it = fLt->AddItem(fListLevel, name, obj);
            Long64_t bsize, fsize, objsize = 0;
            TString objinfo = obj->GetObjectInfo(1, 1);
            TString infos = obj->GetName();
            infos += "\n";
            infos += obj->GetTitle();
            if (!objinfo.IsNull() && !objinfo.BeginsWith("x=")) {
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
            if (it)
               it->SetTipText(infos.Data());
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Browse object. This, in turn, will trigger the calling of
/// TRootBrowserLite::Add() which will fill the IconBox and the tree.
/// Emits signal "BrowseObj(TObject*)".

void TRootBrowserLite::BrowseObj(TObject *obj)
{
   TGPosition pos = fIconBox->GetPagePosition();
   Emit("BrowseObj(TObject*)", (Long_t)obj);

   if (obj != gROOT) {
      if (!fLt->FindItemByObj(fLt->GetFirstItem(), obj)) {
         fListLevel = 0;
         Add(obj);
         fListLevel = fLt->FindItemByObj(fLt->GetFirstItem(), obj);
         fLt->HighlightItem(fListLevel);
         if (obj->IsFolder())
            fLt->OpenItem(fListLevel);
         fLt->ClearViewPort();
         fLt->AdjustPosition(fListLevel);
      }
   }

   if (obj->IsFolder()) fIconBox->RemoveAll();
   obj->Browse(fBrowser);
   if ((fListLevel && obj->IsFolder()) || (!fListLevel && (obj == gROOT))) {
      fIconBox->Refresh();
   }

   if (fBrowser) {
      fBrowser->SetRefreshFlag(kFALSE);
   }
   UpdateDrawOption();

   fIconBox->SetHsbPosition(pos.fX);
   fIconBox->SetVsbPosition(pos.fY);
}

////////////////////////////////////////////////////////////////////////////////
/// add new draw option to the "history"

void TRootBrowserLite::UpdateDrawOption()
{
   TString opt = GetDrawOption();
   TGListBox *lb = fDrawOption->GetListBox();
   TGLBContainer *lbc = (TGLBContainer *)lb->GetContainer();

   TIter next(lbc->GetList());
   TGFrameElement *el;

   while ((el = (TGFrameElement *)next())) {
      TGTextLBEntry *lbe = (TGTextLBEntry *)el->fFrame;
      if (lbe->GetText()->GetString() == opt) {
         return;
      }
   }

   Int_t nn = fDrawOption->GetNumberOfEntries() + 1;
   fDrawOption->AddEntry(opt.Data(), nn);
   fDrawOption->Select(nn);
}

////////////////////////////////////////////////////////////////////////////////
/// returns pointer to fIconBox object

TGFileContainer *TRootBrowserLite::GetIconBox() const
{
   return (TGFileContainer*)fIconBox;
}

////////////////////////////////////////////////////////////////////////////////
/// Really delete the browser and the this GUI.

void TRootBrowserLite::ReallyDelete()
{
   fBrowser->SetBrowserImp(0);
   delete this;
}

////////////////////////////////////////////////////////////////////////////////
/// In case window is closed via WM we get here.

void TRootBrowserLite::CloseWindow()
{
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Display in statusbar total number of objects and number of
/// selected objects in IconBox.

void TRootBrowserLite::DisplayTotal(Int_t total, Int_t selected)
{
   char tmp[64];
   const char *fmt;

   if (selected)
      fmt = "%d Object%s, %d selected.";
   else
      fmt = "%d Object%s.";

   snprintf(tmp, 64, fmt, total, (total == 1) ? "" : "s", selected);
   fStatusBar->SetText(tmp, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Display current directory in second label, fLbl2.

void TRootBrowserLite::DisplayDirectory()
{
   char *p, path[1024];

   fLt->GetPathnameFromItem(fListLevel, path, 12);
   p = path;
   while (*p && *(p+1) == '/') ++p;
   if (!p[0])
      fLbl2->SetText(new TGString("Contents of \".\""));
   else
      fLbl2->SetText(new TGString(Form("Contents of \"%s\"", p)));
   fListHdr->Layout();

   // Get full pathname for FS combobox (previously truncated to 12 levels deep)
   fLt->GetPathnameFromItem(fListLevel, path);
   p = path;
   while (*p && *(p+1) == '/') ++p;
   fFSComboBox->Update(p);

   if (fListLevel) {
      // disable/enable up level navigation
      TGButton *btn = fToolBar->GetButton(kOneLevelUp);
      TString dirname = gSystem->GetDirName(p);

      TObject *obj = (TObject*)fListLevel->GetUserData();
      Bool_t disableUp = dirname == "/";

      // normal file directory
      if (disableUp && (obj) && (obj->IsA() == TSystemDirectory::Class())) {
         disableUp = strlen(p) == 1;
      }
      btn->SetState(disableUp ? kButtonDisabled : kButtonUp);
      AddToHistory(fListLevel);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Execute default action for selected object (action is specified
/// in the $HOME/.root.mimes or $ROOTSYS/etc/root.mimes file.
/// Emits signal "ExecuteDefaultAction(TObject*)".

void TRootBrowserLite::ExecuteDefaultAction(TObject *obj)
{
   TRootBrowserCursorSwitcher cursorSwitcher(fIconBox, fLt);
   char action[512];
   fBrowser->SetDrawOption(GetDrawOption());
   TVirtualPad *wasp = gPad ? (TVirtualPad*)gPad->GetCanvas() : 0;
   TFile *wasf = gFile;

   // Special case for file system objects...
   if (obj->IsA() == TSystemFile::Class() ||
       obj->InheritsFrom("TRemoteObject")) {
      TString act;
      TString ext = obj->GetName();

      if (fClient->GetMimeTypeList()->GetAction(obj->GetName(), action)) {
         act = action;
         act.ReplaceAll("%s", obj->GetName());
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
         Emit("ExecuteDefaultAction(TObject*)", (Long_t)obj);
      }

      // special case for remote object: browse real object
      if (obj->InheritsFrom("TRemoteObject") && ext.EndsWith(".root")) {
         TRootBrowserCursorSwitcher cursorSwitcher2(fIconBox, fLt);
         gApplication->SetBit(TApplication::kProcessRemotely);
         gApplication->ProcessLine("((TApplicationServer *)gApplication)->BrowseFile(0);");
         Refresh();
      }
      ////////// new TFile was opened. Add it to the browser /////
      if (gFile && (wasf != gFile) && ext.EndsWith(".root")) {
         TGListTreeItem *itm = fLt->FindChildByData(0, gROOT->GetListOfFiles());

         if (itm) {
            fLt->ClearHighlighted();
            fListLevel = itm;
            ListTreeHighlight(fListLevel);
            fLt->OpenItem(fListLevel);
            itm = fLt->AddItem(fListLevel, gFile->GetName());
            itm->SetUserData(gFile);
            fClient->NeedRedraw(fLt, kTRUE);
            return;
         }
      }

      // only valid for local text files
      if (!obj->InheritsFrom("TRemoteObject"))
         BrowseTextFile(obj->GetName());

      /////////////// cache and change file's icon ///////////////////////
      TVirtualPad *nowp = gPad ? (TVirtualPad*)gPad->GetCanvas() : 0;

      if (fIconBox->fAutoThumbnail && nowp && (nowp != wasp)) {
         TSystemFile *sf = (TSystemFile*)obj;
         const TGPicture *pic, *spic;

         TString path;
         if (gSystem->IsAbsoluteFileName(sf->GetName())) {
            path = sf->GetName();
         } else {
            char *buf = gSystem->ConcatFileName(gSystem->WorkingDirectory(), sf->GetName());
            path = buf;
            delete [] buf;
         }

         TIconBoxThumb *thumb = (TIconBoxThumb*)fIconBox->fThumbnails->FindObject(path.Data());

         if (thumb) {
            spic = thumb->fSmall;
            pic = thumb->fLarge;
         } else {
            TImage *img = TImage::Create();
            nowp->Modified();
            nowp->Update();
            img->FromPad(nowp);

            if (!img->IsValid()) {
               return;
            }

            static const UInt_t sz = 72;
            UInt_t w = sz;
            UInt_t h = sz;

            auto imgw = img->GetWidth();
            auto imgh = img->GetHeight();

            if (imgw > imgh) {
               h = (imgh*sz) / (imgw > 0 ? imgw : 1);
            } else {
               w = (imgw*sz) / (imgh > 0 ? imgh : 1);
            }

            w = (w < 54) ? 54 : w;
            h = (h < 54) ? 54 : h;

            img->Scale(w, h);
            img->Merge(img, "tint");   // contrasting
            img->DrawBox(0, 0, w, h, "#ffff00", 1); // yellow frame

            pic = fClient->GetPicturePool()->GetPicture(path.Data(), img->GetPixmap(), 0);
            img->Scale(w/3, h/3);
            spic = fClient->GetPicturePool()->GetPicture(path.Data(), img->GetPixmap(), 0);

            thumb = new TIconBoxThumb(path.Data(), spic, pic);
            fIconBox->fThumbnails->Add(thumb);
            delete img;
         }
      }
      return;
   }

   // For other objects the default action is still hard coded in
   // their Browse() member function.
}

////////////////////////////////////////////////////////////////////////////////
/// Handle menu and other command generated by the user.

Bool_t TRootBrowserLite::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
   TRootHelpDialog *hd;
   TRootBrowserCursorSwitcher *cursorSwitcher = 0;
   TDirectory *tdir = 0;
   TString cmd;

   if (GET_SUBMSG(msg) != kCT_SELCHANGED) {
      cursorSwitcher = new TRootBrowserCursorSwitcher(fIconBox, fLt);
   }

   TObject *obj;
   TGListTreeItem *item = 0;

   gVirtualX->Update();

   switch (GET_MSG(msg)) {

      case kC_COMMAND:

         switch (GET_SUBMSG(msg)) {

            case kCM_BUTTON:
               // fallthrough
            case kCM_MENU:

               switch ((ERootBrowserCommands)parm1) {
                  // Handle File menu items...
                  case kFileNewBrowserLite:
                     new TBrowser("Browser", "ROOT Object Browser");
                     break;
                  case kFileNewBrowser:
                     gEnv->SetValue("Browser.Name", "TRootBrowser");
                     new TBrowser();
                     gEnv->SetValue("Browser.Name", "TRootBrowserLite");
                     break;
                  case kFileNewCanvas:
                     gROOT->MakeDefCanvas();
                     break;
                  case kFileNewBuilder:
                     TGuiBuilder::Instance();
                     break;
                  case kFileOpen:
                     {
                        static TString dir(".");
                        TGFileInfo fi;
                        fi.fFileTypes = gOpenTypes;
                        fi.SetIniDir(dir);
                        new TGFileDialog(fClient->GetDefaultRoot(), this,
                                         kFDOpen,&fi);
                        dir = fi.fIniDir;
                        if (fi.fMultipleSelection && fi.fFileNamesList) {
                           TObjString *el;
                           TIter next(fi.fFileNamesList);
                           while ((el = (TObjString *) next())) {
                              new TFile(el->GetString(), "update");
                           }
                        }
                        else if (fi.fFilename) {
                           new TFile(fi.fFilename, "update");
                        }
                     }
                     break;
                  case kFileSave:
                  case kFileSaveAs:
                     break;
                  case kFilePrint:
                     break;
                  case kFileCloseBrowser:
                     SendCloseMessage();
                     break;
                  case kFileQuit:
                     gApplication->Terminate(0);
                     break;

                  // Handle View menu items...
                  case kViewToolBar:
                     if (fViewMenu->IsEntryChecked(kViewToolBar))
                        ShowToolBar(kFALSE);
                     else
                        ShowToolBar();
                     break;
                  case kViewStatusBar:
                     if (fViewMenu->IsEntryChecked(kViewStatusBar))
                        ShowStatusBar(kFALSE);
                     else
                        ShowStatusBar();
                     break;
                  case kViewLargeIcons:
                  case kViewSmallIcons:
                  case kViewList:
                  case kViewDetails:
                     SetViewMode((Int_t)parm1);
                     break;
                  case kViewHidden:
                     if (fBrowser->TestBit(TBrowser::kNoHidden)) {
                        fViewMenu->CheckEntry(kViewHidden);
                        fBrowser->SetBit(TBrowser::kNoHidden, kFALSE);
                     } else {
                        fViewMenu->UnCheckEntry(kViewHidden);
                        fBrowser->SetBit(TBrowser::kNoHidden, kTRUE);
                     }
                     Refresh(kTRUE);
                     break;
                  case kViewArrangeByName:
                  case kViewArrangeByType:
                  case kViewArrangeBySize:
                  case kViewArrangeByDate:
                     SetSortMode((Int_t)parm1);
                     break;
                  case kViewLineUp:
                     break;
                  case kViewRefresh:
                     Refresh(kTRUE);
                     break;
                  case kViewGroupLV:
                     if (!fViewMenu->IsEntryChecked(kViewGroupLV)) {
                        fViewMenu->CheckEntry(kViewGroupLV);
                        TString gv = gEnv->GetValue("Browser.GroupView", "10000");
                        Int_t igv = atoi(gv.Data());

                        if (igv > 10) {
                           fIconBox->SetGroupSize(igv);
                        }
                     } else {
                        fViewMenu->UnCheckEntry(kViewGroupLV);
                        fIconBox->SetGroupSize(10000000); // very large value
                     }
                     break;

                  // Handle Option menu items...
                  case kOptionShowCycles:
                     printf("Currently the browser always shows all cycles\n");
                     break;

                  case kOptionAutoThumbnail:
                     if (fOptionMenu->IsEntryChecked(kOptionAutoThumbnail)) {
                        fOptionMenu->UnCheckEntry(kOptionAutoThumbnail);
                        fIconBox->fThumbnails->Delete();
                        fIconBox->fAutoThumbnail = kFALSE;
                        Refresh(kTRUE);
                     } else {
                        fOptionMenu->CheckEntry(kOptionAutoThumbnail);
                        fIconBox->fAutoThumbnail = kTRUE;
                     }
                     break;

                  // Handle toolbar button...
                  case kOneLevelUp:
                  {
                     if (fBrowseTextFile) {
                        HideTextEdit();
                        break;
                     }
                     if (!fListLevel || !fListLevel->IsActive()) break;

                     if (fListLevel && fIconBox->WasGrouped()) {
                        if (fListLevel) {
                           item = fListLevel->GetParent();
                           if (item) fListLevel = item;
                           obj = (TObject *) fListLevel->GetUserData();
                           HighlightListLevel();
                           if (obj) BrowseObj(obj);
                        }

                        fClient->NeedRedraw(fLt, kTRUE);
                        break;
                     }
                     if (fListLevel) item = fListLevel->GetParent();


                     if (item) {
                        fListLevel = item;
                        obj = (TObject *)fListLevel->GetUserData();
                        HighlightListLevel();
                        DisplayDirectory();
                        if (obj) BrowseObj(obj);
                        fClient->NeedRedraw(fLt, kTRUE);
                     } else {
                        obj = (TObject *)fListLevel->GetUserData();
                        if (obj) ToSystemDirectory(gSystem->GetDirName(obj->GetTitle()));
                     }
                     break;
                  }

                  // toolbar buttons
                  case kHistoryBack:
                     HistoryBackward();
                     break;
                  case kHistoryForw:
                     HistoryForward();
                     break;

                  case kViewFind:
                     Search();
                     break;

                  // Handle Help menu items...
                  case kHelpAbout:
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
                        hd = new TRootHelpDialog(this, str, 600, 400);
                        hd->SetText(gHelpAbout);
                        hd->Popup();
#endif
#endif
                     }
                     break;
                  case kHelpOnCanvas:
                     hd = new TRootHelpDialog(this, "Help on Canvas...", 600, 400);
                     hd->SetText(gHelpCanvas);
                     hd->Popup();
                     break;
                  case kHelpOnMenus:
                     hd = new TRootHelpDialog(this, "Help on Menus...", 600, 400);
                     hd->SetText(gHelpPullDownMenus);
                     hd->Popup();
                     break;
                  case kHelpOnGraphicsEd:
                     hd = new TRootHelpDialog(this, "Help on Graphics Editor...", 600, 400);
                     hd->SetText(gHelpGraphicsEditor);
                     hd->Popup();
                     break;
                  case kHelpOnBrowser:
                     hd = new TRootHelpDialog(this, "Help on Browser...", 600, 400);
                     hd->SetText(gHelpBrowserLite);
                     hd->Popup();
                     break;
                  case kHelpOnObjects:
                     hd = new TRootHelpDialog(this, "Help on Objects...", 600, 400);
                     hd->SetText(gHelpObjects);
                     hd->Popup();
                     break;
                  case kHelpOnPS:
                     hd = new TRootHelpDialog(this, "Help on PostScript...", 600, 400);
                     hd->SetText(gHelpPostscript);
                     hd->Popup();
                     break;
                  case kHelpOnRemote:
                     hd = new TRootHelpDialog(this, "Help on Browser...", 600, 400);
                     hd->SetText(gHelpRemote);
                     hd->Popup();
                     break;
                  default:
                     break;
               }
               break;
            case kCM_COMBOBOX:
               if (parm1 == kFSComboBox) {
                  TGTreeLBEntry *e = (TGTreeLBEntry *) fFSComboBox->GetSelectedEntry();
                  if (e) {
                     const char *dirname = e->GetPath()->GetString();
                     item = fLt->FindItemByPathname(dirname);
                     if (item) {
                        fListLevel = item;
                        HighlightListLevel();
                        DisplayDirectory();
                        fClient->NeedRedraw(fLt, kTRUE);
                     } else {
                        ToSystemDirectory(dirname);
                     }
                  }
               }
               break;
            default:
               break;
         }

         break;

      case kC_LISTTREE:
         switch (GET_SUBMSG(msg)) {

            case kCT_ITEMCLICK:
               // tell coverity that parm1 is a Long_t, and not an enum (even
               // if we compare it with an enum value) and the meaning of
               // parm1 depends on GET_MSG(msg) and GET_SUBMSG(msg)
               // coverity[mixed_enums]
               if (((EMouseButton)parm1 == kButton1) ||
                   ((EMouseButton)parm1 == kButton3)) {
                  HideTextEdit();
                  TGListTreeItem *item2;
                  TObject *obj2 = 0;
                  if ((item2 = fLt->GetSelected()) != 0 ) {
                     ListTreeHighlight(item2);
                     fStatusBar->SetText("", 1);   // clear
                  }
                  if (item2 && parm1 == kButton3) {
                     Int_t x = (Int_t)(parm2 & 0xffff);
                     Int_t y = (Int_t)((parm2 >> 16) & 0xffff);
                     obj2 = (TObject *) item2->GetUserData();
                     if (obj2) {
                        if (obj2->InheritsFrom("TTree")) {
                           // if a tree not attached to any directory (e.g. in a TFolder)
                           // then attach it to the current directory (gDirectory)
                           cmd = TString::Format("((TTree *)0x%lx)->GetDirectory();",
                                                 (ULong_t)obj2);
                           tdir = (TDirectory *)gROOT->ProcessLine(cmd.Data());
                           if (!tdir) {
                              cmd = TString::Format("((TTree *)0x%lx)->SetDirectory(gDirectory);",
                                                    (ULong_t)obj2);
                              gROOT->ProcessLine(cmd.Data());
                           }
                        }
                        fBrowser->GetContextMenu()->Popup(x, y, obj2, fBrowser);
                     }
                  }
                  fClient->NeedRedraw(fLt);
                  fListView->LayoutHeader(0);
                  fListView->Layout();
               }
               break;

            case kCT_ITEMDBLCLICK:
               if (parm1 == kButton1) {
                  if (fBrowseTextFile) {
                     HideTextEdit();
                  }
                  if (fListLevel && fIconBox->WasGrouped()) {
                     TObject *obj2;
                     TGListTreeItem *item2;

                     if (fListLevel) {
                        item2 = fListLevel->GetParent();
                        if (item2) fListLevel = item2;

                        obj2 = (TObject *) fListLevel->GetUserData();
                        HighlightListLevel();
                        if (obj2) {
                           BrowseObj(obj2);
                        }
                     }
                     break;
                  }
               }

            default:
               break;
         }
         break;

      case kC_CONTAINER:
         switch (GET_SUBMSG(msg)) {

            case kCT_ITEMCLICK:
               if (fIconBox->NumSelected() == 1) {
                  // display title of selected object
                  TGFileItem *item2;
                  void *p = 0;
                  if ((item2 = (TGFileItem *)fIconBox->GetNextSelected(&p)) != 0) {
                     TObject *obj2 = (TObject *)item2->GetUserData();

                     TGListTreeItem *itm = 0;
                     if (!fListLevel) itm = fLt->GetFirstItem();
                     else itm = fListLevel->GetFirstChild();
                     //Bool_t found = kFALSE;

                     while (itm) {
                        if (itm->GetUserData() == obj2) break;
                        itm = itm->GetNextSibling();
                     }

                     if (itm) {
                        if ((fListLevel && fListLevel->IsOpen()) || !fListLevel) {
                           fLt->ClearHighlighted();
                           fLt->HighlightItem(itm);
                           fClient->NeedRedraw(fLt, kTRUE);
                        }
                     }

                     if (obj2) fStatusBar->SetText(obj2->GetName(), 1);
                  }
               }
               if (parm1 == kButton3) {
                  // show context menu for selected object
                  if (fIconBox->NumSelected() == 1) {
                     void *p = 0;
                     TGFileItem *item2;
                     if ((item2 = (TGFileItem *) fIconBox->GetNextSelected(&p)) != 0) {
                        Int_t x = (Int_t)(parm2 & 0xffff);
                        Int_t y = (Int_t)((parm2 >> 16) & 0xffff);
                        TObject *obj2 = (TObject *)item2->GetUserData();
                        if (obj2) {
                           if (obj2->IsA() == TKey::Class()) {
                              TKey *key = (TKey*)obj2;
                              TClass *cl = TClass::GetClass(key->GetClassName());
                              TString name = key->GetName();
                              name += ";";
                              name += key->GetCycle();
                              //void *add = gROOT->FindObject((char *) name.Data());//key->GetName());
                              void *add = gDirectory->FindObjectAny((char *) name.Data());
                              if (cl->IsTObject()) {
                                 obj2 = (TObject*)add; // cl->DynamicCast(TObject::Class(),startadd);
                                 item2->SetUserData(obj2);
                              } else {
                                 Error("ProcessMessage","do not support non TObject (like %s) yet",
                                       cl->GetName());
                                 break;
                              }
                           }
                           if (obj2 && obj2->InheritsFrom("TTree")) {
                              // if a tree not attached to any directory (e.g. in a TFolder)
                              // then attach it to the current directory (gDirectory)
                              cmd = TString::Format("((TTree *)0x%lx)->GetDirectory();",
                                                    (ULong_t)obj2);
                              tdir = (TDirectory *)gROOT->ProcessLine(cmd.Data());
                              if (!tdir) {
                                 cmd = TString::Format("((TTree *)0x%lx)->SetDirectory(gDirectory);",
                                                       (ULong_t)obj2);
                                 gROOT->ProcessLine(cmd.Data());
                              }
                           }
                           fBrowser->GetContextMenu()->Popup(x, y, obj2, fBrowser);
                        }
                     }
                  }
               }
               break;
            case kCT_ITEMDBLCLICK:
               if (parm1 == kButton1) {
                  if (fIconBox->NumSelected() == 1) {
                     void *p = 0;
                     TGFileItem *item2;
                     if ((item2 = (TGFileItem *) fIconBox->GetNextSelected(&p)) != 0) {
                        TObject *obj2 = (TObject *)item2->GetUserData();
                        if (obj2) {
                           DoubleClicked(obj2);
                           IconBoxAction(obj2);
                        }
                        delete cursorSwitcher;
                        return kTRUE; //
                     }
                  }
               }
               break;
            case kCT_SELCHANGED:
               DisplayTotal((Int_t)parm1, (Int_t)parm2);
               break;
            default:
               break;
         }

         break;

      default:
         break;
   }

   delete cursorSwitcher;

   fClient->NeedRedraw(fIconBox);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Make object associated with item the current directory.

void TRootBrowserLite::Chdir(TGListTreeItem *item)
{
   if (item) {
      TGListTreeItem *i = item;
      TString dir;
      while (i) {
         TObject *obj = (TObject*) i->GetUserData();
         if (obj) {
            if (obj->IsA() == TDirectoryFile::Class()) {
               dir = "/" + dir;
               dir = obj->GetName() + dir;
            }
            if (obj->IsA() == TFile::Class()) {
               dir = ":/" + dir;
               dir = obj->GetName() + dir;
            }
            if (obj->IsA() == TKey::Class()) {
               if (strcmp(((TKey*)obj)->GetClassName(), "TDirectoryFile") == 0) {
                  dir = "/" + dir;
                  dir = obj->GetName() + dir;
               }
            }
         }
         i = i->GetParent();
      }

      if (gDirectory && dir.Length()) gDirectory->cd(dir.Data());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// helper method  to track history

void TRootBrowserLite::HighlightListLevel()
{
   if (!fListLevel) return;

   fLt->ClearHighlighted();
   fLt->HighlightItem(fListLevel);
}

////////////////////////////////////////////////////////////////////////////////
/// helper method to track history

void TRootBrowserLite::AddToHistory(TGListTreeItem *item)
{
   TGButton *btn = fToolBar->GetButton(kHistoryBack);

   if (!item || (fHistoryCursor &&
       (item == ((TRootBrowserHistoryCursor*)fHistoryCursor)->fItem))) return;

   TRootBrowserHistoryCursor *cur = (TRootBrowserHistoryCursor*)fHistoryCursor;

   while ((cur = (TRootBrowserHistoryCursor*)fHistory->After(fHistoryCursor))) {
      fHistory->Remove(cur);
      delete cur;
   }

   cur = new TRootBrowserHistoryCursor(item);
   fHistory->Add(cur);
   fHistoryCursor = cur;
   btn->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// clear navigation history

void TRootBrowserLite::ClearHistory()
{
   fHistory->Delete();
   TGButton *btn = fToolBar->GetButton(kHistoryBack);
   TGButton *btn2 = fToolBar->GetButton(kHistoryForw);
   btn->SetState(kButtonDisabled);
   btn2->SetState(kButtonDisabled);
}

////////////////////////////////////////////////////////////////////////////////
/// go to the past

Bool_t TRootBrowserLite::HistoryBackward()
{
   if (fBrowseTextFile) {
      HideTextEdit();
      return kFALSE;
   }
   TRootBrowserHistoryCursor *cur = (TRootBrowserHistoryCursor*)fHistory->Before(fHistoryCursor);
   TGButton *btn = fToolBar->GetButton(kHistoryBack);
   TGButton *btn2 = fToolBar->GetButton(kHistoryForw);

   if (!cur) {
      btn->SetState(kButtonDisabled);
      return kFALSE;
   }

   fLt->ClearHighlighted();
   fHistoryCursor = cur;
   fListLevel = cur->fItem;
   ListTreeHighlight(fListLevel);
   fLt->AdjustPosition();
   fClient->NeedRedraw(fLt, kTRUE);

   btn2->SetState(kButtonUp);
   cur = (TRootBrowserHistoryCursor*)fHistory->Before(fHistoryCursor);
   if (!cur) {
      btn->SetState(kButtonDisabled);
      return kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///  go to the future

Bool_t TRootBrowserLite::HistoryForward()
{
   if (fBrowseTextFile) {
      HideTextEdit();
      return kFALSE;
   }

   TRootBrowserHistoryCursor *cur = (TRootBrowserHistoryCursor*)fHistory->After(fHistoryCursor);
   TGButton *btn = fToolBar->GetButton(kHistoryForw);
   TGButton *btn2 = fToolBar->GetButton(kHistoryBack);

   if (!cur) {
      btn->SetState(kButtonDisabled);
      return kFALSE;
   }

   fLt->ClearHighlighted();
   fHistoryCursor = cur;
   fListLevel = cur->fItem;
   ListTreeHighlight(fListLevel);
   fLt->AdjustPosition();
   fClient->NeedRedraw(fLt, kTRUE);

   btn2->SetState(kButtonUp);

   cur = (TRootBrowserHistoryCursor*)fHistory->After(fHistoryCursor);
   if (!cur) {
      btn->SetState(kButtonDisabled);
      return kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// delete list tree item, remove it from history

void TRootBrowserLite::DeleteListTreeItem(TGListTreeItem *item)
{
   ((TRootBrowserHistory*)fHistory)->DeleteItem(item);
   fLt->DeleteItem(item);
}

////////////////////////////////////////////////////////////////////////////////
/// Open tree item and list in iconbox its contents.

void TRootBrowserLite::ListTreeHighlight(TGListTreeItem *item)
{
   if (item) {
      TObject *obj = (TObject *) item->GetUserData();

      if (obj) {
         if (obj->IsA() == TKey::Class()) {

            TKey *key = (TKey *)obj;
            TString name = obj->GetName();
            name += ";";
            name += key->GetCycle();
            Chdir(item->GetParent());
            //TObject *k_obj = gROOT->FindObject(name);
            TObject *k_obj = gDirectory->FindObjectAny(name);

            if (k_obj) {
               item->SetUserData(k_obj);
               obj = k_obj;
            }
         } else if (obj->InheritsFrom(TDirectoryFile::Class())) {
            Chdir(item->GetParent());
         }
         else if (obj->InheritsFrom("TApplicationRemote")) {
            if (!gApplication->GetAppRemote()) {
               gROOT->ProcessLine(Form(".R %s", item->GetText()));
               if (gApplication->GetAppRemote()) {
                  Getlinem(kInit, TString::Format("\n%s:root [0]",
                           gApplication->GetAppRemote()->ApplicationName()));
               }
            }
         }
         else if (obj->InheritsFrom("TRemoteObject")) {
            // special case for remote object
            TRemoteObject *robj = (TRemoteObject *)obj;
            // the real object is a TKey
            if (!strcmp(robj->GetClassName(), "TKey")) {
               TGListTreeItem *parent = item;
               TRemoteObject *probj = (TRemoteObject *)parent->GetUserData();
               // find the TFile remote object containing the TKey
               while ( probj && strcmp(probj->GetClassName(), "TFile")) {
                  parent = parent->GetParent();
                  probj = (TRemoteObject *)parent->GetUserData();
               }
               if (probj) {
                  // remotely browse file (remotely call TFile::cd())
                  gApplication->SetBit(TApplication::kProcessRemotely);
                  gApplication->ProcessLine(
                     TString::Format("((TApplicationServer *)gApplication)->BrowseFile(\"%s\");",
                                     probj->GetName()));
               }
            }
         }
         if (item->GetParent() && item->GetParent()->GetUserData() &&
            ((TObject *)item->GetParent()->GetUserData())->InheritsFrom("TApplicationRemote")) {
            // switch to remote session
            if (!gApplication->GetAppRemote()) {
               gROOT->ProcessLine(Form(".R %s", item->GetParent()->GetText()));
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
                  gROOT->ProcessLine(Form(".R %s", top->GetText()));
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

         if (!fListLevel || !fListLevel->IsActive()) {
            fListLevel = item;
            BrowseObj(obj);
            fLt->HighlightItem(fListLevel);
         }
      }
      DisplayDirectory();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// display  directory

void TRootBrowserLite::ToSystemDirectory(const char *dirname)
{
   TString dir = dirname;

   if (fListLevel) {
      TObject* obj = (TObject*)fListLevel->GetUserData();

      if (obj && (obj->IsA() == TSystemDirectory::Class())) {
         TObject* old = obj;
         fListLevel->Rename(dir.Data());
         obj = new TSystemDirectory(dir.Data(), dir.Data());
         while (fListLevel->GetFirstChild())
            fLt->RecursiveDeleteItem(fListLevel->GetFirstChild(),
                                     fListLevel->GetFirstChild()->GetUserData());

         fListLevel->SetUserData(obj);
         gROOT->GetListOfBrowsables()->Remove(old);
         delete old;
         gROOT->GetListOfBrowsables()->Add(obj);
         fTreeLock = kTRUE;
         BrowseObj(obj);
         fTreeLock = kFALSE;

         fClient->NeedRedraw(fLt, kTRUE);
         fClient->NeedRedraw(fIconBox);
         DisplayDirectory();
         //gSystem->ChangeDirectory(dir.Data());
         fStatusBar->SetText(dir.Data(), 1);
         ClearHistory();   // clear browsing history
      }
   }
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// sets drawing option

void TRootBrowserLite::SetDrawOption(Option_t *option)
{
   fDrawOption->GetTextEntry()->SetText(option);
}

////////////////////////////////////////////////////////////////////////////////
/// returns drawing option

Option_t *TRootBrowserLite::GetDrawOption() const
{
   return fDrawOption->GetTextEntry()->GetText();
}
////////////////////////////////////////////////////////////////////////////////
/// Emits signal when double clicking on icon.

void TRootBrowserLite::DoubleClicked(TObject *obj)
{
   Emit("DoubleClicked(TObject*)", (Long_t)obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Emits signal when double clicking on icon.

void TRootBrowserLite::Checked(TObject *obj, Bool_t checked)
{
   Long_t args[2];

   args[0] = (Long_t)obj;
   args[1] = checked;

   Emit("Checked(TObject*,Bool_t)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Default action when double clicking on icon.

void TRootBrowserLite::IconBoxAction(TObject *obj)
{
   Bool_t browsable = kFALSE;
   if (obj) {

      TRootBrowserCursorSwitcher cursorSwitcher(fIconBox, fLt);

      Bool_t useLock = kTRUE;

      if (obj->IsA()->GetMethodWithPrototype("Browse", "TBrowser*"))
         browsable = kTRUE;

      if (obj->InheritsFrom("TLeaf")) {
         TObject *dir = (TObject *)gROOT->ProcessLine(Form("((%s *)0x%lx)->GetBranch()->GetDirectory();",
                                                      obj->ClassName(), (ULong_t)obj));
         if (!dir) {
            browsable = kFALSE;
         }
      }
      if (obj->InheritsFrom("TBranchElement")) {
         TObject *dir = (TObject *)gROOT->ProcessLine(Form("((%s *)0x%lx)->GetDirectory();",
                                                      obj->ClassName(), (ULong_t)obj));
         if (!dir) {
            browsable = kFALSE;
         }
      }

      if (obj->InheritsFrom("TKey")) {
         TKey *key = dynamic_cast<TKey*>(obj);
         if (key && key->GetClassName() && (!strcmp(key->GetClassName(), "TFormula")))
            browsable = kFALSE;
      }

      if (obj->IsA() == TSystemDirectory::Class()) {
         useLock = kFALSE;

         TString t(obj->GetName());
         if (t == ".") goto out;
         if (t == "..") {
            if (fListLevel && fListLevel->GetParent()) {
               fListLevel = fListLevel->GetParent();
               obj = (TObject*)fListLevel->GetUserData();
               if (fListLevel->GetParent()) {
                  fListLevel = fListLevel->GetParent();
               } else  {
                  obj = (TObject*)fListLevel->GetUserData();
                  fListLevel = 0;
               }
            } else {
               TString dirname = gSystem->GetDirName(gSystem->pwd());
               ToSystemDirectory(dirname);
               return;
            }
         }
      }

      if (obj && obj->IsFolder()) {
         fIconBox->RemoveAll();
         TGListTreeItem *itm = nullptr;

         if (fListLevel) {
            fLt->OpenItem(fListLevel);
            itm = fListLevel->GetFirstChild();
         } else {
            itm = fLt->GetFirstItem();
         }

         while (itm && (itm->GetUserData() != obj)) {
            itm = itm->GetNextSibling();
         }

         if (!itm && fListLevel) {
            // special case for remote objects
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
               // add the remote object only if not already in the list
               if ((!fLt->FindChildByName(fListLevel, obj->GetName())) &&
                   (!fLt->FindChildByData(fListLevel, obj))) {
                  itm = fLt->AddItem(fListLevel, obj->GetName());
                  if (itm) itm->SetUserData(obj);
               }
               else {
                  // set the current item to the one found in the list
                  itm = fLt->FindChildByData(fListLevel, obj) ?
                     fLt->FindChildByData(fListLevel, obj) :
                     fLt->FindChildByName(fListLevel, obj->GetName());
               }
            }
            else {
               itm = fLt->AddItem(fListLevel, obj->GetName());
               if (itm) itm->SetUserData(obj);
            }
         }

         if (itm) {
            fListLevel = itm;
            DisplayDirectory();
            TObject *kobj = (TObject *)itm->GetUserData();

            if (kobj && kobj->IsA() == TKey::Class()) {
               Chdir(fListLevel->GetParent());
               //kobj = gROOT->FindObject(kobj->GetName());
               kobj = gDirectory->FindObjectAny(kobj->GetName());

               if (kobj) {
                  TGListTreeItem *parent = fListLevel->GetParent();
                  DeleteListTreeItem(fListLevel);
                  TGListTreeItem *kitem = fLt->AddItem(parent, kobj->GetName(), kobj);
                  if (kitem) {
                     obj = kobj;
                     useLock = kFALSE;
                     kitem->SetUserData(kobj);
                     fListLevel = kitem;
                  } else
                     fListLevel = parent;
               }
            }
            HighlightListLevel();
         }
      }

      if (browsable) {
         if (useLock) fTreeLock = kTRUE;
         Emit("BrowseObj(TObject*)", (Long_t)obj);
         if (obj) obj->Browse(fBrowser);
         if (useLock) fTreeLock = kFALSE;
      }

out:
      if (obj && obj->IsA() != TSystemFile::Class()) {
         if (obj->IsFolder()) {
            fIconBox->Refresh();
         }

         if (fBrowser) {
            fBrowser->SetRefreshFlag(kFALSE);
         }

         fClient->NeedRedraw(fIconBox);
         fClient->NeedRedraw(fLt, kTRUE);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively remove object from browser.

void TRootBrowserLite::RecursiveRemove(TObject *obj)
{
   // don't delete fIconBox items here (it's status will be updated
   // via TBrowser::Refresh() which should be called once all objects have
   // been removed.

   TGListTreeItem *item = fLt->FindItemByObj(fLt->GetFirstItem(), obj);
   if (item == 0)
      return;
   if (fListLevel && (item == fListLevel)) {
      TGListTreeItem *parent = item->GetParent();
      if (parent) {
         fListLevel = parent;
         fLt->ClearHighlighted();
         fLt->HighlightItem(fListLevel);
         fLt->OpenItem(fListLevel);
      }
      else
         fListLevel = 0;
   }
   DeleteListTreeItem(item);
}

////////////////////////////////////////////////////////////////////////////////
/// Refresh the browser contents.

void TRootBrowserLite::Refresh(Bool_t force)
{
   Bool_t refresh = fBrowser && fBrowser->GetRefreshFlag();

   if (fTextEdit && !gROOT->IsExecutingMacro() && force) {
      fTextEdit->LoadFile(fTextFileName.Data());
      fClient->NeedRedraw(fTextEdit);
      return;
   }

   if ( (refresh || force) && !fIconBox->WasGrouped()
      && fIconBox->NumItems()<fIconBox->GetGroupSize() ) {

      TRootBrowserCursorSwitcher cursorSwitcher(fIconBox, fLt);
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

      // Refresh the IconBox
      if (fListLevel) {
         TObject *obj = (TObject *)fListLevel->GetUserData();
         if (obj) {
            fTreeLock = kTRUE;
            BrowseObj(obj);
            fTreeLock = kFALSE;
         }
      }
   }
   fClient->NeedRedraw(fLt, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Show or hide toolbar.

void TRootBrowserLite::ShowToolBar(Bool_t show)
{
   if (show) {
      ShowFrame(fToolBar);
      ShowFrame(fToolBarSep);
      fViewMenu->CheckEntry(kViewToolBar);
   } else {
      HideFrame(fToolBar);
      HideFrame(fToolBarSep);
      fViewMenu->UnCheckEntry(kViewToolBar);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Show or hide statusbar.

void TRootBrowserLite::ShowStatusBar(Bool_t show)
{
   if (show) {
      ShowFrame(fStatusBar);
      fViewMenu->CheckEntry(kViewStatusBar);
   } else {
      HideFrame(fStatusBar);
      fViewMenu->UnCheckEntry(kViewStatusBar);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set defaults depending on settings in the user's .rootrc.

void TRootBrowserLite::SetDefaults(const char *iconStyle, const char *sortBy)
{
   const char *opt;

   // IconStyle: big, small, list, details
   if (iconStyle)
      opt = iconStyle;
   else
      opt = gEnv->GetValue("Browser.IconStyle", "small");
   if (!strcasecmp(opt, "big"))
      SetViewMode(kViewLargeIcons, kTRUE);
   else if (!strcasecmp(opt, "small"))
      SetViewMode(kViewSmallIcons, kTRUE);
   else if (!strcasecmp(opt, "list"))
      SetViewMode(kViewList, kTRUE);
   else if (!strcasecmp(opt, "details"))
      SetViewMode(kViewDetails, kTRUE);
   else
      SetViewMode(kViewSmallIcons, kTRUE);

   // SortBy: name, type, size, date
   if (sortBy)
      opt = sortBy;
   else
      opt = gEnv->GetValue("Browser.SortBy", "name");
   if (!strcasecmp(opt, "name"))
      SetSortMode(kViewArrangeByName);
   else if (!strcasecmp(opt, "type"))
      SetSortMode(kViewArrangeByType);
   else if (!strcasecmp(opt, "size"))
      SetSortMode(kViewArrangeBySize);
   else if (!strcasecmp(opt, "date"))
      SetSortMode(kViewArrangeByDate);
   else
      SetSortMode(kViewArrangeByName);

   fIconBox->Refresh();
}

////////////////////////////////////////////////////////////////////////////////
/// Set iconbox's view mode and update menu and toolbar buttons accordingly.

void TRootBrowserLite::SetViewMode(Int_t new_mode, Bool_t force)
{
   int i, bnum;
   EListViewMode lv;

   if (force || (fViewMode != new_mode)) {

      switch (new_mode) {
         default:
            if (!force)
               return;
            else
               new_mode = kViewLargeIcons;
            // intentionally no break
         case kViewLargeIcons:
            bnum = 2;
            lv = kLVLargeIcons;
            break;
         case kViewSmallIcons:
            bnum = 3;
            lv = kLVSmallIcons;
            break;
         case kViewList:
            bnum = 4;
            lv = kLVList;
            break;
         case kViewDetails:
            bnum = 5;
            lv = kLVDetails;
            break;
      }

      fViewMode = new_mode;
      fViewMenu->RCheckEntry(fViewMode, kViewLargeIcons, kViewDetails);

      for (i = 2; i <= 5; ++i)
         gToolBarData[i].fButton->SetState((i == bnum) ? kButtonEngaged : kButtonUp);

      fListView->SetViewMode(lv);
      TGTextButton** buttons = fListView->GetHeaderButtons();
      if ((lv == kLVDetails) && (buttons)) {
         if (!strcmp(fListView->GetHeader(1), "Attributes")) {
            buttons[0]->Connect("Clicked()", "TRootBrowserLite", this,
                                TString::Format("SetSortMode(=%d)", kViewArrangeByName));
            buttons[1]->Connect("Clicked()", "TRootBrowserLite", this,
                                TString::Format("SetSortMode(=%d)", kViewArrangeByType));
            buttons[2]->Connect("Clicked()", "TRootBrowserLite", this,
                                TString::Format("SetSortMode(=%d)", kViewArrangeBySize));
            buttons[5]->Connect("Clicked()", "TRootBrowserLite", this,
                                TString::Format("SetSortMode(=%d)", kViewArrangeByDate));
         }
      }
      fIconBox->AdjustPosition();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set iconbox's sort mode and update menu radio buttons accordingly.

void TRootBrowserLite::SetSortMode(Int_t new_mode)
{
   EFSSortMode smode;

   switch (new_mode) {
      default:
         new_mode = kViewArrangeByName;
         // intentionally no break
      case kViewArrangeByName:
         smode = kSortByName;
         break;
      case kViewArrangeByType:
         smode = kSortByType;
         break;
      case kViewArrangeBySize:
         smode = kSortBySize;
         break;
      case kViewArrangeByDate:
         smode = kSortByDate;
         break;
   }

   fSortMode = new_mode;
   fSortMenu->RCheckEntry(fSortMode, kViewArrangeByName, kViewArrangeByDate);

   fIconBox->Sort(smode);
}

////////////////////////////////////////////////////////////////////////////////
/// starts serach dialog

void TRootBrowserLite::Search()
{
   if (!fTextEdit) {
      fIconBox->Search(kFALSE);
   } else {
      fTextEdit->Search(kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// test

static Bool_t isBinary(const char *str, int len)
{
   for (int i = 0; i < len; i++) {
      char c = str[i];
      if (((c < 32) || (c > 126)) && (c != '\t') && (c != '\r') && (c != '\n')) {
         return kTRUE;
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// hide text edit

void TRootBrowserLite::HideTextEdit()
{
   if (!fTextEdit) return;

   ShowMacroButtons(kFALSE);
   fTextEdit->UnmapWindow();
   fV2->RemoveFrame(fTextEdit);
   fV2->AddFrame(fListView, fExpandLayout);
   TGButton *savbtn = fToolBar->GetButton(kViewSave);
   savbtn->Disconnect();
   fTextEdit->DestroyWindow();
   delete fTextEdit;
   fTextEdit = 0;
   fListView->Resize(fV2->GetWidth(), fV2->GetHeight());
   fV2->MapSubwindows();
   fV2->Layout();
   fBrowseTextFile = kFALSE;
   fTextFileName = "";
}

////////////////////////////////////////////////////////////////////////////////
/// browse text file

void TRootBrowserLite::BrowseTextFile(const char *file)
{
   Bool_t loaded = (fTextEdit != 0);
   if (gSystem->AccessPathName(file, kReadPermission)) {
      if (loaded) {
         HistoryBackward();
      }
      return;
   }
   const int bufferSize = 1024;
   char buffer[bufferSize];

   FILE *fd = fopen(file, "rb");
   if (fd == 0) {
      if (loaded) {
         HistoryBackward();
      }
      return;
   }
   int sz = fread(buffer, 1, bufferSize, fd);
   fclose(fd);

   if ((sz > 0) && isBinary(buffer, sz)) {
      if (loaded) {
         HistoryBackward();
      }
      return;
   }

   if (!fTextEdit) {
      fTextEdit = new TGTextEdit(fV2, fV2->GetWidth(), fV2->GetHeight(),
                                 kSunkenFrame | kDoubleBorder);
      TColor *col = gROOT->GetColor(19);
      if (col)
         fTextEdit->SetBackgroundColor(col->GetPixel());
      if (TGSearchDialog::SearchDialog()) {
         TGSearchDialog::SearchDialog()->Connect("TextEntered(char *)", "TGTextEdit",
                                                 fTextEdit, "Search(char *,Bool_t,Bool_t)");
      }
      fV2->AddFrame(fTextEdit, fExpandLayout);
      TGButton *savbtn = fToolBar->GetButton(kViewSave);
      savbtn->Connect("Released()", "TGTextEdit", fTextEdit, "SaveFile(=0,kTRUE)");
   }
   fTextFileName = file;
   fTextEdit->LoadFile(file);
   if (loaded) return;

   if (fTextFileName.EndsWith(".C")) {
      ShowMacroButtons();
   } else {
      fTextEdit->SetReadOnly();
   }
   fListView->UnmapWindow();
   fV2->RemoveFrame(fListView);
   fTextEdit->MapWindow();
   fV2->MapSubwindows();
   fV2->Layout();
   fBrowseTextFile = kTRUE;

   if (fListLevel) {
      AddToHistory(fListLevel);
   }
   TGButton *btn = fToolBar->GetButton(kHistoryForw);

   if (btn) {
      btn->SetState(kButtonDisabled);
   }

   TGButton *btn2 = fToolBar->GetButton(kHistoryBack);

   if (btn2) {
      btn2->SetState(kButtonUp);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// executed browsed text macro

void TRootBrowserLite::ExecMacro()
{
   char *tmpfile = gSystem->ConcatFileName(gSystem->TempDirectory(),
                                           fTextFileName.Data());

   gROOT->SetExecutingMacro(kTRUE);
   fTextEdit->SaveFile(tmpfile, kFALSE);
   gROOT->Macro(tmpfile);
   gSystem->Unlink(tmpfile);
   delete [] tmpfile;
   gROOT->SetExecutingMacro(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// interrupt browsed macro execution

void TRootBrowserLite::InterruptMacro()
{
   gROOT->SetInterrupt(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// show/hide macro buttons

void TRootBrowserLite::ShowMacroButtons(Bool_t show)
{
   TGButton *bt1 = fToolBar->GetButton(kViewExec);
   TGButton *bt2 = fToolBar->GetButton(kViewInterrupt);
   TGButton *bt3 = fToolBar->GetButton(kViewSave);

   static Bool_t connected = kFALSE;

   if (!show) {
      bt1->UnmapWindow();
      bt2->UnmapWindow();
      bt3->UnmapWindow();
   } else {
      bt1->MapWindow();
      bt2->MapWindow();
      bt3->MapWindow();

      if (!connected && fTextEdit) {
         bt1->Connect("Pressed()", "TRootBrowserLite", this, "ExecMacro()");
         bt2->Connect("Pressed()", "TRootBrowserLite", this, "InterruptMacro()");
         connected = kTRUE;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set text in column col in status bar.

void TRootBrowserLite::SetStatusText(const char *txt, Int_t col)
{
   ShowStatusBar(kTRUE);
   TGStatusBar* status = GetStatusBar();
   if (status!=0) {
      status->SetText(txt, col);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Interface method to the old browser.

TBrowserImp *TRootBrowserLite::NewBrowser(TBrowser *b, const char *title,
                                      UInt_t width, UInt_t height,
                                      Option_t * /*opt*/)
{
   TRootBrowserLite *browser = new TRootBrowserLite(b, title, width, height);
   return (TBrowserImp *)browser;
}

////////////////////////////////////////////////////////////////////////////////
/// Interface method to the old browser.

TBrowserImp *TRootBrowserLite::NewBrowser(TBrowser *b, const char *title, Int_t x,
                                      Int_t y, UInt_t width, UInt_t height,
                                      Option_t * /*opt*/)
{
   TRootBrowserLite *browser = new TRootBrowserLite(b, title, x, y, width, height);
   return (TBrowserImp *)browser;
}
