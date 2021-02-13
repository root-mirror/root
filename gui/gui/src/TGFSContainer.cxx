// @(#)root/gui:$Id$
// Author: Fons Rademakers   19/01/98

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
// TGFileIcon, TGFileEntry, TGFSContainer                               //
//                                                                      //
// Utility classes used by the file selection dialog (TGFSDialog).      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGFSContainer.h"
#include "TGIcon.h"
#include "TGMsgBox.h"
#include "TGMimeTypes.h"
#include "TRegexp.h"
#include "TList.h"
#include "TSystem.h"
#include "TVirtualX.h"
#include "TGDNDManager.h"
#include "TRemoteObject.h"
#include "TBufferFile.h"
#include "TImage.h"
#include "snprintf.h"

#include <ctime>
#include <cstdlib>
#include <iostream>

ClassImp(TGFileItem);
ClassImp(TGFileContainer);

class TViewUpdateTimer : public TTimer {

private:
   TGFileContainer   *fContainer;

public:
   TViewUpdateTimer(TGFileContainer *t, Long_t ms) : TTimer(ms, kTRUE) { fContainer = t; }
   Bool_t Notify();
};



class TGFileIcon : public TGIcon {

protected:
   const TGPicture *fLpic;   // icon picture

   virtual void DoRedraw();

public:
   TGFileIcon(const TGWindow *p, const TGPicture *pic, const TGPicture *lpic,
              UInt_t options = kChildFrame, Pixel_t back = GetWhitePixel()) :
      TGIcon(p, pic, 0, 0, options, back) { fLpic = lpic; }
};



////////////////////////////////////////////////////////////////////////////////

class TGFSFrameElement : public TGFrameElement {
public:
   TGFileContainer  *fContainer;   // file container

   Bool_t IsSortable() const { return kTRUE; }
   Int_t  Compare(const TObject *obj) const;
};

////////////////////////////////////////////////////////////////////////////////
/// Sort frame elements in file selection list view container.

Int_t TGFSFrameElement::Compare(const TObject *obj) const
{
   Int_t type1, type2;

   TGFileItem *f1 = (TGFileItem *) fFrame;
   TGFileItem *f2 = (TGFileItem *) ((TGFrameElement *) obj)->fFrame;

   switch (fContainer->fSortType) {
      default:
      case kSortByName:
         //--- this is not exactly what I want...
         type1 = f1->GetType();
         type2 = f2->GetType();

         //--- use posix macros
         if (R_ISDIR(type1)) type1 = 1;
         else                type1 = 6;

         if (R_ISDIR(type2)) type2 = 1;
         else                type2 = 6;

         if (type1 < type2)  return -1;
         if (type1 > type2)  return  1;
         return strcmp(f1->GetItemName()->GetString(),
                       f2->GetItemName()->GetString());

      case kSortByOwner:
         if (f1->GetUid() != f2->GetUid()) {
            if (f1->GetUid() < f2->GetUid())
               return -1;
            else
               return +1;
         }

         // else sort by name
         type1 = f1->GetType();
         type2 = f2->GetType();

         //--- use posix macros
         if (R_ISDIR(type1)) type1 = 1;
         else                type1 = 6;

         if (R_ISDIR(type2)) type2 = 1;
         else                type2 = 6;

         if (type1 < type2)  return -1;
         if (type1 > type2)  return  1;
         return strcmp(f1->GetItemName()->GetString(),
                       f2->GetItemName()->GetString());

      case kSortByGroup:
         if (f1->GetGid() != f2->GetGid()) {
            if (f1->GetGid() < f2->GetGid())
               return -1;
            else
               return +1;
         }

         // else sort by name
         type1 = f1->GetType();
         type2 = f2->GetType();

         //--- use posix macros
         if (R_ISDIR(type1)) type1 = 1;
         else                type1 = 6;

         if (R_ISDIR(type2)) type2 = 1;
         else                type2 = 6;

         if (type1 < type2)  return -1;
         if (type1 > type2)  return  1;
         return strcmp(f1->GetItemName()->GetString(),
                       f2->GetItemName()->GetString());

      case kSortByType:
         //--- this is not exactly what I want...

         type1 = f1->GetType();
         type2 = f2->GetType();

         //--- use posix macros

         if (R_ISDIR(type1))         type1 = 1;
         else if (R_ISLNK(type1))    type1 = 2;
         else if (R_ISSOCK(type1))   type1 = 3;
         else if (R_ISFIFO(type1))   type1 = 4;
         else if (R_ISREG(type1) && (type1 & kS_IXUSR)) type1 = 5;
         else                        type1 = 6;

         if (R_ISDIR(type2))         type2 = 1;
         else if (R_ISLNK(type2))    type2 = 2;
         else if (R_ISSOCK(type2))   type2 = 3;
         else if (R_ISFIFO(type2))   type2 = 4;
         else if (R_ISREG(type2) && (type2 & kS_IXUSR)) type2 = 5;
         else                        type2 = 6;

         if (type1 < type2) return -1;
         if (type1 > type2) return 1;
         return strcmp(f1->GetItemName()->GetString(),
                       f2->GetItemName()->GetString());

      case kSortBySize:
         if (f1->GetSize() < f2->GetSize()) return -1;
         if (f1->GetSize() > f2->GetSize()) return 1;
         return strcmp(f1->GetItemName()->GetString(),
                       f2->GetItemName()->GetString());

      case kSortByDate:
         time_t loctimeF1 = (time_t) f1->GetModTime();
         // coverity[returned_null]
         struct tm tmF1 = *localtime(&loctimeF1);

         time_t loctimeF2 = (time_t) f2->GetModTime();
         // coverity[returned_null]
         struct tm tmF2 = *localtime(&loctimeF2);

         if ( tmF1.tm_year != tmF2.tm_year )
            return (tmF1.tm_year < tmF2.tm_year) ? +1 : -1;
         else if ( tmF1.tm_mon != tmF2.tm_mon )
            return (tmF1.tm_mon < tmF2.tm_mon) ? +1 : -1;
         else if ( tmF1.tm_mday != tmF2.tm_mday )
            return (tmF1.tm_mday < tmF2.tm_mday) ? +1 : -1;
         else if ( tmF1.tm_hour != tmF2.tm_hour )
            return (tmF1.tm_hour < tmF2.tm_hour) ? +1 : -1;
         else if ( tmF1.tm_min != tmF2.tm_min )
            return (tmF1.tm_min < tmF2.tm_min) ? +1 : -1;
         else if ( tmF1.tm_sec != tmF2.tm_sec )
            return (tmF1.tm_sec < tmF2.tm_sec) ? +1 : -1;
         else
            return 0;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Reset the timer.

Bool_t TViewUpdateTimer::Notify()
{
   fContainer->HandleTimer(0);
   Reset();
   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// Draw icon.

void TGFileIcon::DoRedraw()
{
   TGIcon::DoRedraw();
   if (fLpic) fLpic->Draw(fId, GetBckgndGC()(), 0, 0);
}


////////////////////////////////////////////////////////////////////////////////
/// Create a list view item.

TGFileItem::TGFileItem(const TGWindow *p,
                       const TGPicture *bpic, const TGPicture *blpic,
                       const TGPicture *spic, const TGPicture *slpic,
                       TGString *name, Int_t type, Long64_t size, Int_t uid,
                       Int_t gid, Long_t modtime, EListViewMode viewMode,
                       UInt_t options, ULong_t back) :
   TGLVEntry(p, bpic, spic, name, 0, viewMode, options, back)
{
   FileStat_t buf;

   buf.fMode   = type;
   buf.fSize   = size;
   buf.fUid    = uid;
   buf.fGid    = gid;
   buf.fMtime  = modtime;
   buf.fIsLink = (blpic != 0);  // FIXME: hack...

   Init(blpic, slpic, buf, viewMode);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a list view item.

TGFileItem::TGFileItem(const TGWindow *p,
                       const TGPicture *bpic, const TGPicture *blpic,
                       const TGPicture *spic, const TGPicture *slpic,
                       TGString *name, FileStat_t &stat, EListViewMode viewMode,
                       UInt_t options, ULong_t back) :
   TGLVEntry(p, bpic, spic, name, 0, viewMode, options, back)
{
   Init(blpic, slpic, stat, viewMode);
}

////////////////////////////////////////////////////////////////////////////////
/// Common initializer for file list view item.

void TGFileItem::Init(const TGPicture *blpic, const TGPicture *slpic,
                      FileStat_t &stat, EListViewMode viewMode)
{
   char tmp[256];
   Long64_t fsize, bsize;

   fBuf = 0;
   fDNDData.fData = 0;
   fDNDData.fDataLength = 0;
   fDNDData.fDataType = 0;
   fLcurrent =
   fBlpic = blpic;
   fSlpic = slpic;

   fViewMode = (EListViewMode) -1;
   SetViewMode(viewMode);

   fType    = stat.fMode;
   fSize    = stat.fSize;
   fUid     = stat.fUid;
   fGid     = stat.fGid;
   fModTime = stat.fMtime;
   fIsLink  = stat.fIsLink;

   fSubnames = new TGString* [6];

   // file type
   snprintf(tmp, sizeof(tmp), "%c%c%c%c%c%c%c%c%c%c",
            (fIsLink ?
             'l' :
             R_ISREG(fType) ?
             '-' :
             (R_ISDIR(fType) ?
              'd' :
              (R_ISCHR(fType) ?
               'c' :
               (R_ISBLK(fType) ?
                'b' :
                (R_ISFIFO(fType) ?
                 'p' :
                 (R_ISSOCK(fType) ?
                  's' : '?' )))))),
            ((fType & kS_IRUSR) ? 'r' : '-'),
            ((fType & kS_IWUSR) ? 'w' : '-'),
            ((fType & kS_ISUID) ? 's' : ((fType & kS_IXUSR) ? 'x' : '-')),
            ((fType & kS_IRGRP) ? 'r' : '-'),
            ((fType & kS_IWGRP) ? 'w' : '-'),
            ((fType & kS_ISGID) ? 's' : ((fType & kS_IXGRP) ? 'x' : '-')),
            ((fType & kS_IROTH) ? 'r' : '-'),
            ((fType & kS_IWOTH) ? 'w' : '-'),
            ((fType & kS_ISVTX) ? 't' : ((fType & kS_IXOTH) ? 'x' : '-')));
   fSubnames[0] = new TGString(tmp);

   // file size
   fsize = bsize = fSize;
   if (fsize > 1024) {
      fsize /= 1024;
      if (fsize > 1024) {
         // 3.7MB is more informative than just 3MB
         snprintf(tmp, sizeof(tmp), "%lld.%lldM", fsize/1024, (fsize%1024)/103);
      } else {
         snprintf(tmp, sizeof(tmp), "%lld.%lldK", bsize/1024, (bsize%1024)/103);
      }
   } else {
      snprintf(tmp, sizeof(tmp), "%lld", bsize);
   }
   fSubnames[1] = new TGString(tmp);

   {
      struct UserGroup_t *user_group;

      user_group = gSystem->GetUserInfo(fUid);
      if (user_group) {
         fSubnames[2] = new TGString(user_group->fUser);
         fSubnames[3] = new TGString(user_group->fGroup);
         delete user_group;
      } else {
         fSubnames[2] = new TGString(TString::Format("%d", fUid));
         fSubnames[3] = new TGString(TString::Format("%d", fGid));
      }
   }

   struct tm *newtime;
   time_t loctime = (time_t) fModTime;
   newtime = localtime(&loctime);
   if (newtime) {
      snprintf(tmp, sizeof(tmp), "%d-%02d-%02d %02d:%02d", newtime->tm_year + 1900,
               newtime->tm_mon+1, newtime->tm_mday, newtime->tm_hour,
               newtime->tm_min);
      fSubnames[4] = new TGString(tmp);
   }
   else
      fSubnames[4] = new TGString("1901-01-01 00:00");

   fSubnames[5] = 0;

   int i;
   for (i = 0; fSubnames[i] != 0; ++i)
      ;
   fCtw = new int[i+1];
   fCtw[i] = 0;
   for (i = 0; fSubnames[i] != 0; ++i)
      fCtw[i] = gVirtualX->TextWidth(fFontStruct, fSubnames[i]->GetString(),
                                     fSubnames[i]->GetLength());

   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGFileItem::~TGFileItem()
{
   delete fBuf;
}

////////////////////////////////////////////////////////////////////////////////
/// Set container item view mode.

void TGFileItem::SetViewMode(EListViewMode viewMode)
{
   TGLVEntry::SetViewMode(viewMode);

   if (viewMode == kLVLargeIcons)
      fLcurrent = fBlpic;
   else
      fLcurrent = fSlpic;

   if (fClient) fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw list view container item.

void TGFileItem::DoRedraw()
{
   int ix, iy;

   TGLVEntry::DoRedraw();
   if (!fLcurrent) return;

   if (fViewMode == kLVLargeIcons) {
      ix = (fWidth - fLcurrent->GetWidth()) >> 1;
      iy = 0;
   } else {
      ix = 0;
      iy = (fHeight - fLcurrent->GetHeight()) >> 1;
   }

   fLcurrent->Draw(fId, fNormGC, ix, iy);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle drag and drop enter

Atom_t TGFileItem::HandleDNDEnter(Atom_t *)
{
   if (!IsDNDTarget()) return kNone;
   return gVirtualX->InternAtom("application/root", kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set drag and drop data

void TGFileItem::SetDNDData(TDNDData *data)
{
   if (fDNDData.fDataLength > 0)
      free(fDNDData.fData);
   fDNDData.fData = calloc(sizeof(unsigned char), data->fDataLength);
   if (fDNDData.fData)
      memcpy(fDNDData.fData, data->fData, data->fDataLength);
   fDNDData.fDataLength = data->fDataLength;
   fDNDData.fDataType = data->fDataType;
}

////////////////////////////////////////////////////////////////////////////////
/// Set drag and drop object

void TGFileItem::SetDNDObject(TObject *obj)
{
   if (fDNDData.fDataLength)
      free(fDNDData.fData);
   fBuf->WriteObject(obj);
   fDNDData.fData = fBuf->Buffer();
   fDNDData.fDataLength = fBuf->Length();
   fDNDData.fDataType = gVirtualX->InternAtom("application/root", kFALSE);
}



////////////////////////////////////////////////////////////////////////////////
/// Create a list view container which will hold the contents of
/// the current directory.

TGFileContainer::TGFileContainer(const TGWindow *p, UInt_t w, UInt_t h,
                                 UInt_t options, ULong_t back) :
   TGLVContainer(p, w, h, options, back)
{
   fSortType  = kSortByName;
   fFilter    = 0;
   fMtime     = 0;
   fDirectory = gSystem->WorkingDirectory();
   fRefresh   = new TViewUpdateTimer(this, 1000);
   gSystem->AddTimer(fRefresh);
   fCachePictures = kTRUE;
   fDisplayStat   = kTRUE;
   fCleanups  = new TList;

   fFolder_s = fClient->GetPicture("folder_s.xpm");
   fFolder_t = fClient->GetPicture("folder_t.xpm");
   fApp_s    = fClient->GetPicture("app_s.xpm");
   fApp_t    = fClient->GetPicture("app_t.xpm");
   fDoc_s    = fClient->GetPicture("doc_s.xpm");
   fDoc_t    = fClient->GetPicture("doc_t.xpm");
   fSlink_s  = fClient->GetPicture("slink_s.xpm");
   fSlink_t  = fClient->GetPicture("slink_t.xpm");

   if (!fFolder_s || !fFolder_t ||
       !fApp_s    || !fApp_t    ||
       !fDoc_s    || !fDoc_t    ||
       !fSlink_s  || !fSlink_t)
      Error("TGFileContainer", "required pixmap(s) missing\n");

   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a list view container which will hold the contents of
/// the current directory.

TGFileContainer::TGFileContainer(TGCanvas *p, UInt_t options, ULong_t back) :
   TGLVContainer(p,options, back)
{
   fSortType  = kSortByName;
   fFilter    = 0;
   fMtime     = 0;
   fDirectory = gSystem->WorkingDirectory();
   fRefresh   = new TViewUpdateTimer(this, 1000);
   gSystem->AddTimer(fRefresh);
   fCachePictures = kTRUE;
   fDisplayStat   = kTRUE;
   fCleanups  = new TList;

   fFolder_s = fClient->GetPicture("folder_s.xpm");
   fFolder_t = fClient->GetPicture("folder_t.xpm");
   fApp_s    = fClient->GetPicture("app_s.xpm");
   fApp_t    = fClient->GetPicture("app_t.xpm");
   fDoc_s    = fClient->GetPicture("doc_s.xpm");
   fDoc_t    = fClient->GetPicture("doc_t.xpm");
   fSlink_s  = fClient->GetPicture("slink_s.xpm");
   fSlink_t  = fClient->GetPicture("slink_t.xpm");

   if (!fFolder_s || !fFolder_t ||
       !fApp_s    || !fApp_t    ||
       !fDoc_s    || !fDoc_t    ||
       !fSlink_s  || !fSlink_t)
      Error("TGFileContainer", "required pixmap(s) missing\n");

   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete list view file container.

TGFileContainer::~TGFileContainer()
{
   if (fRefresh) delete fRefresh;
   if (fFilter)  delete fFilter;
   fClient->FreePicture(fFolder_s);
   fClient->FreePicture(fFolder_t);
   fClient->FreePicture(fApp_s);
   fClient->FreePicture(fApp_t);
   fClient->FreePicture(fDoc_s);
   fClient->FreePicture(fDoc_t);
   fClient->FreePicture(fSlink_s);
   fClient->FreePicture(fSlink_t);
   if (fCleanups) {
      TGPicture *pic;
      TIter nextp(fCleanups);
      while ((pic = (TGPicture *)nextp())) {
         fClient->GetPicturePool()->FreePicture(pic);
      }
      fCleanups->Clear();
      delete fCleanups;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add frame to the composite frame.

void TGFileContainer::AddFrame(TGFrame *f, TGLayoutHints *l)
{
   TGFSFrameElement *nw;

   nw = new TGFSFrameElement;
   nw->fFrame     = f;
   nw->fLayout    = l ? l : fgDefaultHints;
   nw->fState     = 1;
   nw->fContainer = this;
   fList->Add(nw);
}

////////////////////////////////////////////////////////////////////////////////
/// Refresh container contents. Check every 5 seconds to see if the
/// directory modification date has changed.

Bool_t TGFileContainer::HandleTimer(TTimer *)
{
   FileStat_t sbuf;

   if (gSystem->GetPathInfo(fDirectory, sbuf) == 0)
      if (fMtime != (ULong_t)sbuf.fMtime) DisplayDirectory();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set file selection filter.

void TGFileContainer::SetFilter(const char *filter)
{
   if (fFilter) delete fFilter;
   fFilter = new TRegexp(filter, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Sort file system list view container according to sortType.

void TGFileContainer::Sort(EFSSortMode sortType)
{
   fSortType = sortType;

   fList->Sort();

   TGCanvas *canvas = (TGCanvas *) this->GetParent()->GetParent();
   canvas->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Determine the file picture for the given file type.

void TGFileContainer::GetFilePictures(const TGPicture **pic,
             const TGPicture **lpic, Int_t file_type, Bool_t is_link,
             const char *name, Bool_t /*small*/)
{
   static TString cached_ext;
   static const TGPicture *cached_spic = 0;
   static const TGPicture *cached_lpic = 0;
   const char *ext = name ? strrchr(name, '.') : 0;
   *pic = 0;
   *lpic = 0;

   if (fCachePictures && ext && cached_spic && cached_lpic && (cached_ext == ext)) {
      *pic = cached_spic;
      *lpic = cached_lpic;
      if (!is_link) return;
   }

   if (R_ISREG(file_type)) {
      TString fname(name);
      if (is_link && fname.EndsWith(".lnk")) {
         fname.Remove(fname.Length()-4);
      }
      *pic = fClient->GetMimeTypeList()->GetIcon(fname.Data(), kTRUE);
      *lpic = fClient->GetMimeTypeList()->GetIcon(fname.Data(), kFALSE);

      if (*pic) {
         if (!*lpic) *lpic = *pic;
         if (ext) {
            cached_ext = ext;
            cached_spic = *pic;
            cached_lpic = *lpic;
            if (!is_link) return;
         }
      }
   } else {
      *pic = 0;
   }

   if (*pic == 0) {
      *pic = fDoc_t;
      *lpic = fDoc_s;

      if (R_ISREG(file_type) && (file_type) & kS_IXUSR) {
         *pic = fApp_t;
         *lpic = fApp_s;
      }
      if (R_ISDIR(file_type)) {
         *pic = fFolder_t;
         *lpic = fFolder_s;
      }
   }
   if (is_link) {
      TImage *img1, *img2;
      if (*pic && *lpic) {
         TString lnk_name;
         img1 = TImage::Create();
         if (img1) {
            img1->SetImage(((const TGPicture *)*pic)->GetPicture(),
                           ((const TGPicture *)*pic)->GetMask());
            img2 = TImage::Open("slink_t.xpm");
            if (img2) img1->Merge(img2);
            lnk_name = ((const TGPicture *)*pic)->GetName();
            lnk_name.Prepend("lnk_");
            *pic = fClient->GetPicturePool()->GetPicture(lnk_name.Data(),
                                 img1->GetPixmap(), img1->GetMask());
            fCleanups->Add(((TObject *)*pic));
            if (img2) delete img2;
            delete img1;
         }
         img1 = TImage::Create();
         if (img1) {
            img1->SetImage(((const TGPicture *)*lpic)->GetPicture(),
                           ((const TGPicture *)*lpic)->GetMask());
            img2 = TImage::Open("slink_s.xpm");
            if (img2) img1->Merge(img2);
            lnk_name = ((const TGPicture *)*lpic)->GetName();
            lnk_name.Prepend("lnk_");
            *lpic = fClient->GetPicturePool()->GetPicture(lnk_name.Data(),
                                 img1->GetPixmap(), img1->GetMask());
            fCleanups->Add(((TObject *)*lpic));
            if (img2) delete img2;
            delete img1;
         }
      }
      else {
         *pic = fSlink_t;
         *lpic = fSlink_s;
      }
   }

   cached_lpic = 0;
   cached_spic = 0;
   cached_ext = "";
}

////////////////////////////////////////////////////////////////////////////////
/// Change current directory.

void TGFileContainer::ChangeDirectory(const char *path)
{
   TString savdir = gSystem->WorkingDirectory();
   gSystem->ChangeDirectory(fDirectory.Data());   // so path of ".." will work
   char *exppath = gSystem->ExpandPathName(path);
   if (gSystem->ChangeDirectory(exppath)) {
      fDirectory = gSystem->WorkingDirectory();
      gSystem->ChangeDirectory(savdir.Data());
      DisplayDirectory();
   }
   delete[] exppath;
}

////////////////////////////////////////////////////////////////////////////////
/// Display the contents of the current directory in the container.
/// This can be used to refresh the contents of the window.

void TGFileContainer::DisplayDirectory()
{
   RemoveAll();
   CreateFileList();

   // This automatically calls layout
   Sort(fSortType);

   // Make TGExplorerMainFrame display total objects in status bar
   SendMessage(fMsgWindow, MK_MSG(kC_CONTAINER, kCT_SELCHANGED),
               fTotal, fSelected);

   MapSubwindows();
}

////////////////////////////////////////////////////////////////////////////////
/// This function creates the file list from current dir.

void TGFileContainer::CreateFileList()
{
   TString savdir = gSystem->WorkingDirectory();
   if (!gSystem->ChangeDirectory(fDirectory.Data())) return;

   FileStat_t sbuf;
   if (gSystem->GetPathInfo(".", sbuf) == 0)
      fMtime = sbuf.fMtime;

   void *dirp;
   if ((dirp = gSystem->OpenDirectory(".")) == 0) {
      gSystem->ChangeDirectory(savdir.Data());
      return;
   }

   const char *name;
   while ((name = gSystem->GetDirEntry(dirp)) != 0 && fDisplayStat) {
      if (strcmp(name, ".") && strcmp(name, ".."))
         AddFile(name);
      gSystem->ProcessEvents();
   }
   gSystem->FreeDirectory(dirp);

   gSystem->ChangeDirectory(savdir.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Add file in container.

TGFileItem *TGFileContainer::AddFile(const char *name,  const TGPicture *ipic,
                                     const TGPicture *ilpic)
{
   TString     filename;
   TGFileItem *item = 0;
   const TGPicture *spic, *slpic;
   TGPicture *pic, *lpic;

   FileStat_t sbuf;

   if (gSystem->GetPathInfo(name, sbuf)) {
      if (sbuf.fIsLink) {
         Info("AddFile", "Broken symlink of %s.", name);
      } else {
         TString msg;
         msg.Form("Can't read file attributes of \"%s\": %s.",
                  name, gSystem->GetError());
         new TGMsgBox(fClient->GetDefaultRoot(), GetMainFrame(),
                      "Error", msg.Data(), kMBIconStop, kMBOk);
      }
      return item;
   }

   filename = name;
   if (R_ISDIR(sbuf.fMode) || fFilter == 0 ||
       (fFilter && filename.Index(*fFilter) != kNPOS)) {

      if (ipic && ilpic) { // dynamic icons
         spic = ipic;
         slpic = ilpic;
      } else {
         GetFilePictures(&spic, &slpic, sbuf.fMode, sbuf.fIsLink, name, kTRUE);
      }

      pic = (TGPicture*)spic; pic->AddReference();
      lpic = (TGPicture*)slpic; lpic->AddReference();

      item = new TGFileItem(this, lpic, slpic, spic, pic,
                            new TGString(gSystem->BaseName(name)),
                            sbuf, fViewMode);
      AddItem(item);
   }

   return item;
}

////////////////////////////////////////////////////////////////////////////////
/// Add remote file in container.

TGFileItem *TGFileContainer::AddRemoteFile(TObject *obj, const TGPicture *ipic,
                                           const TGPicture *ilpic)
{
   TString     filename;
   TGFileItem *item = 0;
   const TGPicture *spic, *slpic;
   TGPicture *pic, *lpic;

   FileStat_t sbuf;

   TRemoteObject *robj = (TRemoteObject *)obj;

   robj->GetFileStat(&sbuf);
   filename = robj->GetName();

   if (R_ISDIR(sbuf.fMode) || fFilter == 0 ||
       (fFilter && filename.Index(*fFilter) != kNPOS)) {

      if (ipic && ilpic) { // dynamic icons
         spic = ipic;
         slpic = ilpic;
      } else {
         GetFilePictures(&spic, &slpic, sbuf.fMode, sbuf.fIsLink, filename, kTRUE);
      }

      pic = (TGPicture*)spic; pic->AddReference();
      lpic = (TGPicture*)slpic; lpic->AddReference();

      item = new TGFileItem(this, lpic, slpic, spic, pic, new TGString(filename),
                            sbuf, fViewMode);
      AddItem(item);
   }
   return item;
}

////////////////////////////////////////////////////////////////////////////////
/// stop refresh  timer

void TGFileContainer::StopRefreshTimer()
{
   if (fRefresh) delete fRefresh;
   fRefresh = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// start refreshing

void TGFileContainer::StartRefreshTimer(ULong_t msec)
{
   fRefresh = new TViewUpdateTimer(this, msec);
   gSystem->AddTimer(fRefresh);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a file container widget as a C++ statement(s) on output stream out.

void TGFileContainer::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fBackground != GetDefaultFrameBackground()) SaveUserColor(out, option);

   char quote = '"';
   out << std::endl << "   // container frame" << std::endl;
   out << "   TGFileContainer *";

   if ((fParent->GetParent())->InheritsFrom(TGCanvas::Class())) {
      out << GetName() << " = new TGFileContainer(" << GetCanvas()->GetName();
   } else {
      out << GetName() << " = new TGFileContainer(" << fParent->GetName();
      out << "," << GetWidth() << "," << GetHeight();
   }

   if (fBackground == GetDefaultFrameBackground()) {
      if (GetOptions() == kSunkenFrame) {
         out <<");" << std::endl;
      } else {
         out << "," << GetOptionString() <<");" << std::endl;
      }
   } else {
      out << "," << GetOptionString() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;
   out << "   " << GetCanvas()->GetName() << "->SetContainer("
                << GetName() << ");" << std::endl;
   out << "   " << GetName() << "->DisplayDirectory();" << std::endl;
   out << "   " << GetName() << "->AddFile("<< quote << ".." << quote << ");" << std::endl;
   out << "   " << GetName() << "->StopRefreshTimer();" << std::endl;
}
