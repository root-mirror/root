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
// TGFSComboBox, TGTreeLBEntry                                          //
//                                                                      //
// This is a combo box that is used in the File Selection dialog box.   //
// It will allow the file path selection.                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RConfigure.h"

#include "TGFSComboBox.h"
#include "TGResourcePool.h"
#include "TGPicture.h"
#include "TSystem.h"
#include "TVirtualX.h"
#include "strlcpy.h"

#include <cstdlib>
#include <iostream>

const TGFont *TGTreeLBEntry::fgDefaultFont = nullptr;
TGGC         *TGTreeLBEntry::fgDefaultGC = nullptr;

//--- this is temp here...

struct Lbc_t {
   const char *fName;   // root prefix name
   const char *fPath;   // path
   const char *fPixmap; // picture file
   Int_t       fId;     // widget id
   Int_t       fIndent; // identification level
   Int_t       fFlags;  // flag
};

static struct Lbc_t gLbc[32];

ClassImp(TGTreeLBEntry);
ClassImp(TGFSComboBox);

////////////////////////////////////////////////////////////////////////////////
/// Create a tree (i.e. entry can be indented) listbox entry.
/// The strings text and path are adopted by the entry.

TGTreeLBEntry::TGTreeLBEntry(const TGWindow *p, TGString *text,
                             const TGPicture *pic, Int_t id, TGString *path,
                             GContext_t norm, FontStruct_t font, UInt_t options,
                             ULong_t back) :
   TGLBEntry(p, id, options, back)
{
   if (text && !pic)
      Error("TGTreeLBEntry", "icon not found for entry %s", text->GetString());

   fPic    = pic;
   fSelPic = 0;
   fTWidth = 0;
   fText   = text;
   fPath   = path;

   fNormGC     = norm;
   fFontStruct = font;

   fActive = kFALSE;

   int max_ascent, max_descent;
   if (fText)
      fTWidth = gVirtualX->TextWidth(fFontStruct, fText->GetString(), fText->GetLength());
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);
   fTHeight = max_ascent + max_descent;
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete tree listbox entry.

TGTreeLBEntry::~TGTreeLBEntry()
{
   delete fText;
   delete fPath;
   delete fSelPic;
}

////////////////////////////////////////////////////////////////////////////////
/// Make entry active (highlight picture).

void TGTreeLBEntry::Activate(Bool_t a)
{
   if (fActive == a) return;
   fActive = a;

   if (fActive) {
      fSelPic = new TGSelectedPicture(fClient, fPic);
   } else {
      if (fSelPic) delete fSelPic;
      fSelPic = 0;
   }
   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw the tree listbox entry on pixmap/window.

void TGTreeLBEntry::DrawCopy(Handle_t id, Int_t x, Int_t y)
{
   int ix, iy, lx, ly;

   ix = x;
   iy = y + ((fHeight - fPic->GetHeight()) >> 1);
   lx = x + (int)(fPic->GetWidth() + 4);
   ly = y + (int)((fHeight - (fTHeight+1)) >> 1);

   if (fActive) {
      if (fSelPic) fSelPic->Draw(id, fNormGC, ix, iy);
      gVirtualX->SetForeground(fNormGC, fgDefaultSelectedBackground);
      gVirtualX->FillRectangle(id, fNormGC, lx-2, ly, fWidth-(lx-x), fTHeight+1);
      gVirtualX->SetForeground(fNormGC, fClient->GetResourcePool()->GetSelectedFgndColor());
   } else {
      fPic->Draw(id, fNormGC, ix, iy);
      gVirtualX->SetForeground(fNormGC, fgWhitePixel);
      gVirtualX->FillRectangle(id, fNormGC, lx-2, ly, fWidth-(lx-x), fTHeight+1);
      gVirtualX->SetForeground(fNormGC, fgBlackPixel);
   }

   int max_ascent, max_descent;
   gVirtualX->GetFontProperties(fFontStruct, max_ascent, max_descent);

   fText->Draw(id, fNormGC, lx, ly + max_ascent);
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw the tree listbox entry.

void TGTreeLBEntry::DoRedraw()
{
   DrawCopy(fId, 0, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Return default size of tree listbox entry.

TGDimension TGTreeLBEntry::GetDefaultSize() const
{
   TGDimension isize(fPic->GetWidth(), fPic->GetHeight());
   TGDimension lsize(fTWidth, fTHeight+1);

   return TGDimension(isize.fWidth + lsize.fWidth + 4,
                      TMath::Max(isize.fHeight, lsize.fHeight) + 2);
}

////////////////////////////////////////////////////////////////////////////////
/// Update text and picture of a listbox entry.

void TGTreeLBEntry::Update(TGLBEntry *e)
{
   TGTreeLBEntry *te = (TGTreeLBEntry *) e;

   if (fText) delete fText;
   fText = new TGString(te->GetText());
   fPic = te->GetPicture();
   fTWidth = gVirtualX->TextWidth(fFontStruct, fText->GetString(), fText->GetLength());
   gVirtualX->ClearWindow(fId);
   fClient->NeedRedraw(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Return default font structure.

FontStruct_t TGTreeLBEntry::GetDefaultFontStruct()
{
   if (!fgDefaultFont)
      fgDefaultFont = gClient->GetResourcePool()->GetDefaultFont();
   return fgDefaultFont->GetFontStruct();
}

////////////////////////////////////////////////////////////////////////////////
/// Return default graphics context.

const TGGC &TGTreeLBEntry::GetDefaultGC()
{
   if (!fgDefaultGC)
      fgDefaultGC = new TGGC(*gClient->GetResourcePool()->GetFrameGC());
   return *fgDefaultGC;
}


////////////////////////////////////////////////////////////////////////////////
/// Create a file system combobox showing system directories.

TGFSComboBox::TGFSComboBox(const TGWindow *parent, Int_t id, UInt_t options,
                           ULong_t back) :
   TGComboBox(parent, id, options | kOwnBackground, back)
{
   int   i, indent;
   const TGPicture *pic;
   char  *p;

   SetTopEntry(new TGTreeLBEntry(this, new TGString("Current dir"),
               fClient->GetPicture("folder_t.xpm"), 0),
               new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 4, 0, 0, 0));

   fListBox->GetContainer()->AddInput(kButtonPressMask | kButtonReleaseMask |
                                      kPointerMotionMask);

   //--- first check for the existence of some directories...

   const char *homeDir = gSystem->HomeDirectory();
#ifndef ROOTPREFIX
   const char *rootSys = gSystem->Getenv("ROOTSYS");
#else
   // const char *rootSys = ROOTPREFIX;
#endif

   Int_t idx = 0;
   TList *volumes = gSystem->GetVolumes("all");
   TList *curvol  = gSystem->GetVolumes("cur");
   TString infos;
   const char *curdrive = "";
   if (volumes && curvol) {
      TNamed *named = (TNamed *)curvol->At(0);
      if (named) {
         curdrive = named->GetName();
         infos = named->GetTitle();
         gLbc[idx].fName = StrDup(infos.Data());
         gLbc[idx].fPath = StrDup(Form("%s\\", curdrive));
         if (infos.Contains("Removable"))
            gLbc[idx].fPixmap = StrDup("fdisk_t.xpm");
         else if (infos.Contains("Local"))
            gLbc[idx].fPixmap = StrDup("hdisk_t.xpm");
         else if (infos.Contains("CD"))
            gLbc[idx].fPixmap = StrDup("cdrom_t.xpm");
         else if (infos.Contains("Network"))
            gLbc[idx].fPixmap = StrDup("netdisk_t.xpm");
         else
            gLbc[idx].fPixmap = StrDup("hdisk_t.xpm");
         gLbc[idx].fId     = 1000;
         gLbc[idx].fIndent = 0;
         gLbc[idx].fFlags  = 0;
         ++idx;
      }
      else {
         gLbc[idx].fName = StrDup("Root");
         gLbc[idx].fPath = StrDup("/");
         gLbc[idx].fPixmap = StrDup("hdisk_t.xpm");
         gLbc[idx].fId     = 1000;
         gLbc[idx].fIndent = 1;
         gLbc[idx].fFlags  = 0;
         ++idx;
      }
   }
   else {
      gLbc[idx].fName = StrDup("Root");
      gLbc[idx].fPath = StrDup("/");
      gLbc[idx].fPixmap = StrDup("hdisk_t.xpm");
      gLbc[idx].fId     = 1000;
      gLbc[idx].fIndent = 1;
      gLbc[idx].fFlags  = 0;
      ++idx;
      gLbc[idx].fName = StrDup("Floppy");
      gLbc[idx].fPath = StrDup("/floppy");
      gLbc[idx].fPixmap = StrDup("fdisk_t.xpm");
      gLbc[idx].fId     = 2000;
      gLbc[idx].fIndent = 1;
      gLbc[idx].fFlags  = 0;
      ++idx;
      gLbc[idx].fName = StrDup("CD-ROM");
      gLbc[idx].fPath = StrDup("/cdrom");
      gLbc[idx].fPixmap = StrDup("cdrom_t.xpm");
      gLbc[idx].fId     = 3000;
      gLbc[idx].fIndent = 1;
      gLbc[idx].fFlags  = 0;
      ++idx;
   }
   gLbc[idx].fName   = StrDup("Home");
   gLbc[idx].fPath   = StrDup("$HOME");
   gLbc[idx].fPixmap = StrDup("home_t.xpm");
   gLbc[idx].fId     = (idx+1) * 1000;
   gLbc[idx].fIndent = 1;
   gLbc[idx].fFlags  = 0;
   ++idx;
#ifndef ROOTPREFIX
   gLbc[idx].fName   = StrDup("RootSys");
   gLbc[idx].fPath   = StrDup("$ROOTSYS");
#else
   gLbc[idx].fName   = StrDup(ROOTPREFIX);
   gLbc[idx].fPath   = StrDup(ROOTPREFIX);
#endif
   gLbc[idx].fPixmap = StrDup("root_t.xpm");
   gLbc[idx].fId     = (idx+1) * 1000;
   gLbc[idx].fIndent = 1;
   gLbc[idx].fFlags  = 0;
   ++idx;

   if (volumes && curvol) {
      TIter next(volumes);
      TNamed *drive;
      while ((drive = (TNamed *)next())) {
         if (!strcmp(drive->GetName(), curdrive))
            continue;
         infos = drive->GetTitle();
         gLbc[idx].fName   = StrDup(drive->GetTitle());
         gLbc[idx].fPath   = StrDup(Form("%s\\", drive->GetName()));
         if (infos.Contains("Removable"))
            gLbc[idx].fPixmap = StrDup("fdisk_t.xpm");
         else if (infos.Contains("Local"))
            gLbc[idx].fPixmap = StrDup("hdisk_t.xpm");
         else if (infos.Contains("CD"))
            gLbc[idx].fPixmap = StrDup("cdrom_t.xpm");
         else if (infos.Contains("Network"))
            gLbc[idx].fPixmap = StrDup("netdisk_t.xpm");
         else
            gLbc[idx].fPixmap = StrDup("hdisk_t.xpm");
         gLbc[idx].fId     = (idx+1) * 1000;
         gLbc[idx].fIndent = 0;
         gLbc[idx].fFlags  = 0;
         ++idx;
      }
      delete volumes;
      delete curvol;
   }
   gLbc[idx].fName   = 0;
   gLbc[idx].fPath   = 0;
   gLbc[idx].fPixmap = 0;
   gLbc[idx].fId     = (idx+1) * 1000;
   gLbc[idx].fIndent = 0;
   gLbc[idx].fFlags  = 0;

   for (i = 0; gLbc[i].fPath != 0; ++i) {
      if (strstr(gLbc[i].fPath, "$HOME") != 0) {
         if (homeDir) {
            int hlen = strlen(homeDir);
            int blen = hlen + strlen(gLbc[i].fPath) - 3;
            p = new char[blen];
            strlcpy(p, homeDir, blen);
            strlcat(p, &gLbc[i].fPath[5], blen-strlen(&gLbc[i].fPath[5]));
            delete [] gLbc[i].fPath;
            gLbc[i].fPath = p;
         } else {
            gLbc[i].fFlags = 0;
         }
      }
#ifndef ROOTPREFIX
      // Below should _only_ be called if the prefix isn't set at build
      // time. The code below expands the occurance of `$ROOTSYS' in
      // the table above.  However, in the case of prefix being set at
      // build time, we do not need to expand the prefix, as it is
      // already known, so the entries in the table above are actually
      // fully expanded.
      if (strstr(gLbc[i].fPath, "$ROOTSYS") != 0) {
         // Get the size of the prefix template
         const int plen = 8;
         if (rootSys) {
            int hlen = strlen(rootSys);
            // Allocate enough memory to hold prefix (hlen), and
            // what's in the path (strlen(gLbc[i].fPath)) minus the
            // prefix template size, and one character for terminating
            // null.
            int blen = hlen + strlen(gLbc[i].fPath) - plen + 1;
            p = new char[blen];
            strlcpy(p, rootSys, blen);
            strlcat(p, &(gLbc[i].fPath[plen]), blen-strlen(&gLbc[i].fPath[plen]));
            // Figure out where to put the terminating NULL
            int npos = hlen + strlen(&(gLbc[i].fPath[plen]));
            p[npos] = '\0';
            delete [] gLbc[i].fPath;
            gLbc[i].fPath = p;
         } else {
            gLbc[i].fFlags = 0;
         }
      }
#endif
      if (gSystem->AccessPathName(gLbc[i].fPath, kFileExists) == 0)
         gLbc[i].fFlags = 1;
   }

   //--- then init the contents...

   for (i = 0; gLbc[i].fName != 0; ++i) {
      if (gLbc[i].fFlags) {
         indent = 4 + (gLbc[i].fIndent * 10);
         pic = fClient->GetPicture(gLbc[i].fPixmap);
         if (!pic) Error("TGFSComboBox", "pixmap not found: %s", gLbc[i].fPixmap);
         AddEntry(new TGTreeLBEntry(fListBox->GetContainer(),
                  new TGString(gLbc[i].fName), pic, gLbc[i].fId,
                  new TGString(gLbc[i].fPath)),
                  new TGLayoutHints(kLHintsExpandX | kLHintsTop, indent, 0, 0, 0));
      }
   }
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Update file system combo box.

void TGFSComboBox::Update(const char *path)
{
   char dirname[1024], mpath[1024];
   const char *tailpath = 0;
   int  i, indent_lvl = 0, afterID = -1, sel = -1;

   if (!path) return;

   for (i = 0; gLbc[i].fPath != 0; ++i)
      RemoveEntries(gLbc[i].fId+1, gLbc[i+1].fId-1);

   int len = 0;
   for (i = 0; gLbc[i].fName != 0; ++i) {
      if (gLbc[i].fFlags) {
         int slen = strlen(gLbc[i].fPath);
         if (strncmp(path, gLbc[i].fPath, slen) == 0) {
            if (slen > len) {
               sel = afterID = gLbc[i].fId;
               indent_lvl = gLbc[i].fIndent + 1;
               if ((len > 0) && ((path[slen] == '\\') || (path[slen] == '/') ||
                   (path[slen] == 0)))
                  tailpath = path + slen;
               strlcpy(mpath, gLbc[i].fPath, 1024);
               len = slen;
            }
         }
      }
   }

   if (tailpath && *tailpath) {
      if ((*tailpath == '/') || (*tailpath == '\\')) ++tailpath;
      if (*tailpath)
         while (1) {
            const char *picname;
            const char *semi = strchr(tailpath, '/');
            if (semi == 0) semi = strchr(tailpath, '\\');
            if (semi == 0) {
               strlcpy(dirname, tailpath, 1024);
               picname = "ofolder_t.xpm";
            } else {
               strlcpy(dirname, tailpath, (semi-tailpath)+1);
               picname = "folder_t.xpm";
            }
            if ((mpath[strlen(mpath)-1] != '/') &&
                (mpath[strlen(mpath)-1] != '\\')) {
               strlcat(mpath, "/", 1024-strlen(mpath));
            }
            strlcat(mpath, dirname, 1024-strlen(mpath));
            int indent = 4 + (indent_lvl * 10);
            const TGPicture *pic = fClient->GetPicture(picname);
            if (!pic) Error("Update", "pixmap not found: %s", picname);
            InsertEntry(new TGTreeLBEntry(fListBox->GetContainer(),
                        new TGString(dirname), pic, afterID+1,
                        new TGString(mpath)),
                        new TGLayoutHints(kLHintsExpandX | kLHintsTop,
                                          indent, 0, 0, 0),
                        afterID);
            sel = ++afterID;
            ++indent_lvl;
            if (semi == 0) break;
            tailpath = ++semi;
         }
   }
   if (sel > 0) Select(sel);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a file system combo box as a C++ statement(s) on output stream out.

void TGFSComboBox::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   if (fBackground != GetWhitePixel()) SaveUserColor(out, option);

   out << std::endl << "   // file system combo box" << std::endl;
   out << "   TGFSComboBox *";
   out << GetName() << " = new TGFSComboBox(" << fParent->GetName()
                                          << "," << fWidgetId;
   if (fBackground == GetWhitePixel()) {
      if (GetOptions() == (kHorizontalFrame | kSunkenFrame | kDoubleBorder)) {
         out <<");" << std::endl;
      } else {
         out << "," << GetOptionString() <<");" << std::endl;
      }
   } else {
      out << "," << GetOptionString() << ",ucolor);" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;

   out << "   " << GetName() << "->Resize(" << GetWidth()  << ","
       << GetHeight() << ");" << std::endl;
   out << "   " << GetName() << "->Select(" << GetSelected() << ");" << std::endl;

}
