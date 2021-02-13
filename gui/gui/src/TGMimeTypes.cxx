// @(#)root/gui:$Id: c44ba45210ec143ec5bf9aa1708855c60088e954 $
// Author: Fons Rademakers   18/01/98

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
// TGMimeTypes and TGMime                                               //
//                                                                      //
// This class handles mime types, used by browsers to map file types    //
// to applications and icons. TGMime is internally used by TGMimeType.  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGMimeTypes.h"
#include "TOrdCollection.h"
#include "TSystem.h"
#include "TDatime.h"
#include "TRegexp.h"
#include "strlcpy.h"

ClassImp(TGMimeTypes);

////////////////////////////////////////////////////////////////////////////////
/// Create a mime type cache. Read the mime types file "filename" and
/// built a list of mime types.

TGMimeTypes::TGMimeTypes(TGClient *client, const char *filename)
{
   char     line[1024]   = {0};
   char     mime[1024]   = {0};
   char     pattern[256] = {0};
   char     icon[256]    = {0};
   char     sicon[256]   = {0};
   char     action[256]  = {0};
   char    *s;

   fClient   = client;
   fFilename = filename;
   fChanged  = kFALSE;
   fList     = new TOrdCollection(50);

   FILE *mfp;
   mfp = fopen(filename, "r");
   if (!mfp) {
      Warning("TGMimeTypes", "error opening mime type file %s", filename);
      return;
   }

   int cnt = 0;
   while (fgets(line, 1024, mfp)) {
      s = line;
      s[strlen(line)-1] = 0;       // strip off trailing \n
      while (*s == ' ') s++;       // strip leading blanks
      if (*s == '#') continue;     // skip comments
      if (!s[0]) continue;    // skip empty lines

      if (*s == '[') {
         strlcpy(mime, line, 1024);
         cnt = 0;
         continue;
      }
      if (!strncmp(s, "pattern", 7)) {
         if (!(s = strchr(line, '='))) {
            Error("TGMimeTypes", "malformed pattern line, = missing");
            pattern[0] = 0;
         } else {
            s++;
            s = Strip(s);
            strlcpy(pattern, s, 256);
            delete [] s;
         }
         cnt++;
      } else if (!strncmp(s, "icon", 4)) {
         if (!(s = strchr(line, '='))) {
            Error("TGMimeTypes", "malformed icon line, = missing");
            icon[0] = 0;
         } else {
            s++;
            s = Strip(s);
            char *s2;
            if ((s2 = strchr(s, ' '))) {
               *s2 = 0;
               strlcpy(icon, s, 256);
               s2++;
               s2 = Strip(s2);
               strlcpy(sicon, s2, 256);
               delete [] s2;
            } else {
               strlcpy(icon, s, 256);
               strlcpy(sicon, s, 256);
            }
            delete [] s;
         }
         cnt++;
      } else if (!strncmp(s, "action", 6)) {
         if (!(s = strchr(line, '='))) {
            Error("TGMimeTypes", "malformed action line, = missing");
            action[0] = 0;
         } else {
            s++;
            s = Strip(s);
            strlcpy(action, s, 256);
            delete [] s;
         }
         cnt++;
      }

      if (cnt == 3) {
         if (strchr(pattern, ' ')) {
            char *tmppattern = strtok(pattern, " ");
            while (tmppattern && (*tmppattern != ' ')) {
               AddType(mime, tmppattern, icon, sicon, action);
               tmppattern = strtok(0, " ");
            }
         } else {
            AddType(mime, pattern, icon, sicon, action);
         }
      }
   }

   fclose(mfp);

   fChanged = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete mime type pool.

TGMimeTypes::~TGMimeTypes()
{
   if (fChanged) SaveMimes();
   fList->Delete();
   delete fList;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

TGMimeTypes::TGMimeTypes(const TGMimeTypes& gmt) :
   TObject(gmt),
   fClient(gmt.fClient),
   fFilename(gmt.fFilename),
   fChanged(gmt.fChanged),
   fList(gmt.fList)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator

TGMimeTypes& TGMimeTypes::operator=(const TGMimeTypes& gmt)
{
   if(this!=&gmt) {
      TObject::operator=(gmt);
      fClient=gmt.fClient;
      fFilename=gmt.fFilename;
      fChanged=gmt.fChanged;
      fList=gmt.fList;
   }
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// Given a filename find the matching mime type object.

TGMime *TGMimeTypes::Find(const char *filename)
{
   if (!filename) return 0;

   TString fn = filename;

   TGMime  *mime;
   TIter    next(fList);

   while ((mime = (TGMime *) next()))
      if (fn.Index(*(mime->fReg)) != kNPOS) return mime;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return icon belonging to mime type of filename.

const TGPicture *TGMimeTypes::GetIcon(const char *filename, Bool_t small_icon)
{
   TGMime *mime;
   const TGPicture *mypic = 0;

   if ((mime = Find(filename))) {
      Bool_t thumb = (mime->fType == "[thumbnail]");
      if (small_icon) {
         mypic = thumb ? fClient->GetPicture(mime->fSIcon.Data(), 32, 32) :
                         fClient->GetPicture(mime->fSIcon.Data(), 16, 16);
      } else {
         mypic = thumb ? fClient->GetPicture(mime->fIcon.Data(), 64, 64) :
                         fClient->GetPicture(mime->fIcon.Data(), 32, 32);
      }
      return mypic;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return in action the mime action string belonging to filename.

Bool_t TGMimeTypes::GetAction(const char *filename, char *action)
{
   TGMime *mime;

   action[0] = 0;
   if ((mime = Find(filename))) {
      strlcpy(action, mime->fAction.Data(), 512);
      return (strlen(action) > 0);
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return in type the mime type belonging to filename.

Bool_t TGMimeTypes::GetType(const char *filename, char *type)
{
   TGMime *mime;

   memset(type, 0, strlen(type));
   if ((mime = Find(filename))) {
      strlcpy(type, mime->fType.Data(), 256);
      return (strlen(type) > 0);
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Print list of mime types.

void TGMimeTypes::Print(Option_t *) const
{
   TGMime *m;
   TIter next(fList);

   while ((m = (TGMime *) next())) {
      printf("Type:    %s\n", m->fType.Data());
      printf("Pattern: %s\n", m->fPattern.Data());
      if (m->fIcon != m->fSIcon)
         printf("Icon:    %s %s\n", m->fIcon.Data(), m->fSIcon.Data());
      else
         printf("Icon:    %s\n",    m->fIcon.Data());
      printf("Action:  %s\n", m->fAction.Data());
      printf("------------\n\n");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save mime types in user's mime type file.

void TGMimeTypes::SaveMimes()
{
   TString filename;
#ifdef WIN32
   filename.Form("%s\\.root.mimes", gSystem->HomeDirectory());
#else
   filename.Form("%s/.root.mimes", gSystem->HomeDirectory());
#endif

   FILE *fp = fopen(filename.Data(), "wt");

   if (!fp) {
      Error("SaveMimes", "can not open %s to store mime types", filename.Data());
      return;
   }

   TDatime dt;
   fprintf(fp, "# %s written on %s\n\n", filename.Data(), dt.AsString());

   TGMime *m;
   TIter next(fList);

   while ((m = (TGMime *) next())) {
      fprintf(fp, "%s\n",            m->fType.Data());
      fprintf(fp, "pattern = %s\n",  m->fPattern.Data());
      if (m->fIcon != m->fSIcon)
         fprintf(fp, "icon = %s %s\n", m->fIcon.Data(), m->fSIcon.Data());
      else
         fprintf(fp, "icon = %s\n",    m->fIcon.Data());
      fprintf(fp, "action = %s\n\n", m->fAction.Data());
   }

   fclose(fp);

   fChanged = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a mime type to the list of mime types.

void TGMimeTypes::AddType(const char *type, const char *pattern, const char *icon,
                          const char *sicon, const char *action)
{
   TGMime *mime = new TGMime;

   mime->fType    = type;
   mime->fPattern = pattern;
   mime->fIcon    = icon;
   mime->fSIcon   = sicon;
   mime->fAction  = action;

   mime->fReg = new TRegexp(pattern, kTRUE);

   fList->AddFirst(mime);

   fChanged = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete mime object.

TGMime::~TGMime()
{
   delete fReg;
}
