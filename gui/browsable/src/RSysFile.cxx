/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/// \file ROOT/RSysFile.cxx
/// \ingroup rbrowser
/// \author Sergey Linev <S.Linev@gsi.de>
/// \date 2019-10-15
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!


#include <ROOT/Browsable/RSysFile.hxx>

#include <ROOT/Browsable/RSysFileItem.hxx>
#include <ROOT/Browsable/RWrapper.hxx>
#include <ROOT/Browsable/RLevelIter.hxx>
#include <ROOT/Browsable/RProvider.hxx>

#include "ROOT/RLogger.hxx"

#include "TROOT.h"
#include "TList.h"
#include "TBase64.h"
#include "snprintf.h"

#include <sstream>
#include <fstream>
#include <algorithm>

#ifdef _MSC_VER
#include <windows.h>
#include <tchar.h>
#endif

using namespace std::string_literals;

using namespace ROOT::Experimental::Browsable;

namespace ROOT {
namespace Experimental {
namespace Browsable {


#ifdef _MSC_VER
bool IsWindowsLink(const std::string &path)
{
   return (path.length() > 4) && (path.rfind(".lnk") == path.length() - 4);
}
#endif


/** \class RSysDirLevelIter
\ingroup rbrowser

Iterator over files in in sub-directory
*/


class RSysDirLevelIter : public RLevelIter {
   std::string fPath;        ///<! fully qualified path without final slash
   void *fDir{nullptr};      ///<! current directory handle
   std::string fCurrentName; ///<! current file name
   std::string fItemName;    ///<! current item name
   FileStat_t fCurrentStat;  ///<! stat for current file name

   /** Open directory for listing */
   bool OpenDir()
   {
      if (fDir)
         CloseDir();

#ifdef _MSC_VER
      // on Windows path can be redirected via .lnk therefore get real path name before OpenDirectory,
      // otherwise such realname will not be known for us
      if (IsWindowsLink(fPath)) {
         char *realWinPath = gSystem->ExpandPathName(fPath.c_str());
         if (realWinPath) fPath = realWinPath;
         delete [] realWinPath;
      }
#endif

      fDir = gSystem->OpenDirectory(fPath.c_str());

#ifdef _MSC_VER
    // Directory can be an soft link (not as .lnk) and should be tried as well
    if (!fDir) {

      auto hFile = CreateFile(fPath.c_str(),  // file to open
                              0,                // open for reading
                              0,                // share for reading
                              0,                // default security
                              OPEN_EXISTING,    // existing file only
                              FILE_FLAG_BACKUP_SEMANTICS, // flag to work with dirs
                              NULL);                 // no attr. template

      if( hFile != INVALID_HANDLE_VALUE) {
         const int BUFSIZE = 2048;
         TCHAR path[BUFSIZE];
         auto dwRet = GetFinalPathNameByHandle( hFile, path, BUFSIZE, VOLUME_NAME_DOS );
         // produced file name may include \\? symbols, which are indicating long file name
         if ((dwRet > 0) && (dwRet < BUFSIZE))
           if ((path[0] == '\\') && (path[1] == '\\') && (path[2] == '?') && (path[3] == '\\')) {
              R__LOG_DEBUG(0, BrowsableLog()) << "Try to open directory " << (path+4) << " instead of " << fPath;
              fDir = gSystem->OpenDirectory(path + 4);
              if (fDir) fPath = path + 4;
           }
      }

      CloseHandle(hFile);
   }

#endif

      if (!fDir) {
         R__LOG_ERROR(BrowsableLog()) << "Fail to open directory " << fPath;
         return false;
      }

      return true;
   }

   /** Close directory for listing */
   void CloseDir()
   {
      if (fDir)
         gSystem->FreeDirectory(fDir);
      fDir = nullptr;
      fCurrentName.clear();
      fItemName.clear();
   }

   /** Return full dir name with appropriate slash at the end */
   std::string FullDirName() const
   {
      std::string path = fPath;
#ifdef _MSC_VER
      const char *slash = "\\";
#else
      const char *slash = "/";
#endif
      if (path.rfind(slash) != path.length() - 1)
         path.append(slash);
      return path;
   }

   /** Check if entry of that name exists */
   bool TestDirEntry(const std::string &name)
   {
      auto testname = name;

      auto path = FullDirName() + testname;

      auto pathinfores = gSystem->GetPathInfo(path.c_str(), fCurrentStat);

#ifdef _MSC_VER
      if (pathinfores && !IsWindowsLink(path)) {
         std::string lpath = path + ".lnk";
         pathinfores = gSystem->GetPathInfo(lpath.c_str(), fCurrentStat);
         if (!pathinfores) testname.append(".lnk");
      }
#endif

      if (pathinfores) {

         if (fCurrentStat.fIsLink) {
            R__LOG_ERROR(BrowsableLog()) << "Broken symlink of " << path;
         } else {
            R__LOG_ERROR(BrowsableLog()) << "Can't read file attributes of \"" <<  path << "\" err:" << gSystem->GetError();
         }
         return false;
      }

      fItemName = fCurrentName = testname;
#ifdef _MSC_VER
      if (IsWindowsLink(fItemName))
         fItemName.resize(fItemName.length() - 4);
#endif
      return true;
   }

   /** Trying to produce next entry */
   bool NextDirEntry()
   {
      fCurrentName.clear();
      fItemName.clear();

      if (!fDir)
         return false;

      while (fCurrentName.empty()) {

         // one have to use const char* to correctly check for nullptr
         const char *name = gSystem->GetDirEntry(fDir);

         if (!name) {
            CloseDir();
            return false;
         }

         std::string sname = name;

         if ((sname == ".") || (sname == ".."))
            continue;

         TestDirEntry(sname);
      }


      return true;
   }

   std::string GetFileExtension(const std::string &fname) const
   {
      auto pos = fname.rfind(".");
      if ((pos != std::string::npos) && (pos < fname.length() - 1) && (pos > 0))
         return fname.substr(pos+1);

      return ""s;
   }

public:
   explicit RSysDirLevelIter(const std::string &path = "") : fPath(path) { OpenDir(); }

   virtual ~RSysDirLevelIter() { CloseDir(); }

   bool Next() override { return NextDirEntry(); }

   bool Find(const std::string &name, int = -1) override
   {
      // ignore index, it is not possible to have duplicated file names

      if (!fDir && !OpenDir())
         return false;

      return TestDirEntry(name);
   }

   std::string GetItemName() const override { return fItemName; }

   /** Returns true if directory or is file format supported */
   bool CanItemHaveChilds() const override
   {
      if (R_ISDIR(fCurrentStat.fMode))
         return true;

      if (RProvider::IsFileFormatSupported(GetFileExtension(fCurrentName)))
         return true;

      return false;
   }

   std::unique_ptr<RItem> CreateItem() override
   {
      auto item = std::make_unique<RSysFileItem>(GetItemName(), CanItemHaveChilds() ? -1 : 0);

      // this is construction of current item
      char tmp[256];

      item->type     = fCurrentStat.fMode;
      item->size     = fCurrentStat.fSize;
      item->uid      = fCurrentStat.fUid;
      item->gid      = fCurrentStat.fGid;
      item->modtime  = fCurrentStat.fMtime;
      item->islink   = fCurrentStat.fIsLink;
      item->isdir    = R_ISDIR(fCurrentStat.fMode);

      if (item->isdir)
         item->SetIcon("sap-icon://folder-blank"s);
      else
         item->SetIcon(RSysFile::GetFileIcon(GetItemName()));

      // file size
      Long64_t _fsize = item->size, bsize = item->size;
      if (_fsize > 1024) {
         _fsize /= 1024;
         if (_fsize > 1024) {
            // 3.7MB is more informative than just 3MB
            snprintf(tmp, sizeof(tmp), "%lld.%lldM", _fsize/1024, (_fsize%1024)/103);
         } else {
            snprintf(tmp, sizeof(tmp), "%lld.%lldK", bsize/1024, (bsize%1024)/103);
         }
      } else {
         snprintf(tmp, sizeof(tmp), "%lld", bsize);
      }
      item->fsize = tmp;

      // modification time
      time_t loctime = (time_t) item->modtime;
      struct tm *newtime = localtime(&loctime);
      if (newtime) {
         snprintf(tmp, sizeof(tmp), "%d-%02d-%02d %02d:%02d", newtime->tm_year + 1900,
                  newtime->tm_mon+1, newtime->tm_mday, newtime->tm_hour,
                  newtime->tm_min);
         item->mtime = tmp;
      } else {
         item->mtime = "1901-01-01 00:00";
      }

      // file type
      snprintf(tmp, sizeof(tmp), "%c%c%c%c%c%c%c%c%c%c",
               (item->islink ?
                'l' :
                R_ISREG(item->type) ?
                '-' :
                (R_ISDIR(item->type) ?
                 'd' :
                 (R_ISCHR(item->type) ?
                  'c' :
                  (R_ISBLK(item->type) ?
                   'b' :
                   (R_ISFIFO(item->type) ?
                    'p' :
                    (R_ISSOCK(item->type) ?
                     's' : '?' )))))),
               ((item->type & kS_IRUSR) ? 'r' : '-'),
               ((item->type & kS_IWUSR) ? 'w' : '-'),
               ((item->type & kS_ISUID) ? 's' : ((item->type & kS_IXUSR) ? 'x' : '-')),
               ((item->type & kS_IRGRP) ? 'r' : '-'),
               ((item->type & kS_IWGRP) ? 'w' : '-'),
               ((item->type & kS_ISGID) ? 's' : ((item->type & kS_IXGRP) ? 'x' : '-')),
               ((item->type & kS_IROTH) ? 'r' : '-'),
               ((item->type & kS_IWOTH) ? 'w' : '-'),
               ((item->type & kS_ISVTX) ? 't' : ((item->type & kS_IXOTH) ? 'x' : '-')));
      item->ftype = tmp;

      struct UserGroup_t *user_group = gSystem->GetUserInfo(item->uid);
      if (user_group) {
         item->fuid = user_group->fUser;
         item->fgid = user_group->fGroup;
         delete user_group;
      } else {
         item->fuid = std::to_string(item->uid);
         item->fgid = std::to_string(item->gid);
      }

      return item;
   }

   /** Returns full information for current element */
   std::shared_ptr<RElement> GetElement() override
   {
      if (!R_ISDIR(fCurrentStat.fMode)) {
         auto extension = GetFileExtension(fCurrentName);

         if (RProvider::IsFileFormatSupported(extension)) {
            auto elem = RProvider::OpenFile(extension, FullDirName() + fCurrentName);
            if (elem) return elem;
         }
      }

      return std::make_shared<RSysFile>(fCurrentStat, FullDirName(), fCurrentName);
   }

};


} // namespace Browsable
} // namespace Experimental
} // namespace ROOT


/////////////////////////////////////////////////////////////////////////////////
/// Get icon for the type of given file name

std::string RSysFile::GetFileIcon(const std::string &fname)
{
    std::string name = fname;
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);

   auto EndsWith = [name](const std::string &suffix) {
      return (name.length() > suffix.length()) ? (0 == name.compare (name.length() - suffix.length(), suffix.length(), suffix)) : false;
   };

   if ((EndsWith(".c")) ||
       (EndsWith(".cpp")) ||
       (EndsWith(".cxx")) ||
       (EndsWith(".c++")) ||
       (EndsWith(".cxx")) ||
       (EndsWith(".h")) ||
       (EndsWith(".hpp")) ||
       (EndsWith(".hxx")) ||
       (EndsWith(".h++")) ||
       (EndsWith(".py")) ||
       (EndsWith(".txt")) ||
       (EndsWith(".cmake")) ||
       (EndsWith(".dat")) ||
       (EndsWith(".log")) ||
       (EndsWith(".xml")) ||
       (EndsWith(".htm")) ||
       (EndsWith(".html")) ||
       (EndsWith(".json")) ||
       (EndsWith(".sh")) ||
       (EndsWith(".md")) ||
       (EndsWith(".css")) ||
       (EndsWith(".js")))
      return "sap-icon://document-text"s;
   if ((EndsWith(".bmp")) ||
       (EndsWith(".gif")) ||
       (EndsWith(".jpeg")) ||
       (EndsWith(".jpg")) ||
       (EndsWith(".png")) ||
       (EndsWith(".svg")))
      return "sap-icon://picture"s;
  if (EndsWith(".root"))
      return "sap-icon://org-chart"s;

   return "sap-icon://document"s;
}


/////////////////////////////////////////////////////////////////////////////////
/// Create file element

RSysFile::RSysFile(const std::string &filename) : fFileName(filename)
{
   if (gSystem->GetPathInfo(fFileName.c_str(), fStat)) {
      if (fStat.fIsLink) {
         R__LOG_ERROR(BrowsableLog()) << "Broken symlink of " << fFileName;
      } else {
         R__LOG_ERROR(BrowsableLog()) << "Can't read file attributes of \"" << fFileName
                                    << "\" err:" << gSystem->GetError();
      }
   }

   auto pos = fFileName.find_last_of("\\/");
   if ((pos != std::string::npos) && (pos < fFileName.length() - 1)) {
      fDirName = fFileName.substr(0, pos+1);
      fFileName.erase(0, pos+1);
   }
}

/////////////////////////////////////////////////////////////////////////////////
/// Create file element with already provided stats information

RSysFile::RSysFile(const FileStat_t &stat, const std::string &dirname, const std::string &filename)
   : fStat(stat), fDirName(dirname), fFileName(filename)
{
}

/////////////////////////////////////////////////////////////////////////////////
/// return file name

std::string RSysFile::GetName() const
{
   return fFileName;
}

/////////////////////////////////////////////////////////////////////////////////
/// Check if file name the same, ignore case on Windows

bool RSysFile::MatchName(const std::string &name) const
{
   auto ownname = GetName();

#ifdef _MSC_VER

   return std::equal(name.begin(), name.end(),
                     ownname.begin(), ownname.end(),
                     [](char a, char b) {
                         return tolower(a) == tolower(b);
                      });
#else

   return ownname == name;

#endif
}

/////////////////////////////////////////////////////////////////////////////////
/// Get default action for the file
/// Either start text editor or image viewer or just do file browsing

RElement::EActionKind RSysFile::GetDefaultAction() const
{
   if (R_ISDIR(fStat.fMode)) return kActBrowse;

   auto icon = GetFileIcon(GetName());
   if (icon == "sap-icon://document-text"s) return kActEdit;
   if (icon == "sap-icon://picture"s) return kActImage;
   if (icon == "sap-icon://org-chart"s) return kActBrowse;
   return kActNone;
}

/////////////////////////////////////////////////////////////////////////////////
/// Returns full file name - including fully qualified path

std::string RSysFile::GetFullName() const
{
   return fDirName + fFileName;
}

/////////////////////////////////////////////////////////////////////////////////
/// Returns iterator for files in directory

std::unique_ptr<RLevelIter> RSysFile::GetChildsIter()
{
   if (!R_ISDIR(fStat.fMode))
      return nullptr;

   return std::make_unique<RSysDirLevelIter>(GetFullName());
}

/////////////////////////////////////////////////////////////////////////////////
/// Returns file content of requested kind

std::string RSysFile::GetContent(const std::string &kind)
{
   if ((GetContentKind(kind) == kText) && (GetFileIcon(GetName()) == "sap-icon://document-text"s)) {
      std::ifstream t(GetFullName());
      return std::string(std::istreambuf_iterator<char>(t), std::istreambuf_iterator<char>());
   }

   if ((GetContentKind(kind) == kImage) && (GetFileIcon(GetName()) == "sap-icon://picture"s)) {
      std::ifstream t(GetFullName(), std::ios::binary);
      std::string content = std::string(std::istreambuf_iterator<char>(t), std::istreambuf_iterator<char>());

      auto encode = TBase64::Encode(content.data(), content.length());

      auto pos = GetName().rfind(".");

      return "data:image/"s  + GetName().substr(pos+1) + ";base64,"s + encode.Data();
   }

   if (GetContentKind(kind) == kFileName) {
      return GetFullName();
   }

   return ""s;
}


/////////////////////////////////////////////////////////////////////////////////
/// Provide top entries for file system
/// On windows it is list of existing drivers, on Linux it is "Files system" and "Home"

RElementPath_t RSysFile::ProvideTopEntries(std::shared_ptr<RGroup> &comp, const std::string &workdir)
{
   std::string seldir = workdir;

   if (seldir.empty())
      seldir = gSystem->WorkingDirectory();

   seldir = gSystem->UnixPathName(seldir.c_str());

   auto volumes = gSystem->GetVolumes("all");
   if (volumes) {
      // this is Windows
      TIter iter(volumes);
      TObject *obj;
      while ((obj = iter()) != nullptr) {
         std::string name = obj->GetName();
         std::string dir = name + "\\"s;
         comp->Add(std::make_shared<Browsable::RWrapper>(name, std::make_unique<RSysFile>(dir)));
      }
      delete volumes;

   } else {
      comp->Add(std::make_shared<Browsable::RWrapper>("Files system", std::make_unique<RSysFile>("/")));

      seldir = "/Files system"s + seldir;

      std::string homedir = gSystem->UnixPathName(gSystem->HomeDirectory());

      if (!homedir.empty())
         comp->Add(std::make_shared<Browsable::RWrapper>("Home", std::make_unique<RSysFile>(homedir)));
   }

   return RElement::ParsePath(seldir);
}
