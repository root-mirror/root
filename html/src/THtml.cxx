// @(#)root/html:$Id$
// Author: Nenad Buncic (18/10/95), Axel Naumann (09/28/01)

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "THtml.h"
#include "RConfigure.h"
#include "Riostream.h"
#include "TBaseClass.h"
#include "TClass.h"
#include "TClassDocOutput.h"
#include "TClassEdit.h"
#include "TClassTable.h"
#include "TDataType.h"
#include "TDocInfo.h"
#include "TDocOutput.h"
#include "TEnv.h"
#include "TInterpreter.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TPRegexp.h"
#include "TRegexp.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TThread.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <set>
#include <fstream>

THtml *gHtml = 0;

//______________________________________________________________________________
//______________________________________________________________________________
namespace {
   class THtmlThreadInfo {
   public:
      THtmlThreadInfo(THtml* html, bool force): fHtml(html), fForce(force) {}
      Bool_t GetForce() const {return fForce;}
      THtml* GetHtml() const {return fHtml;}

   private:
      THtml* fHtml;
      Bool_t fForce;
   };
};


////////////////////////////////////////////////////////////////////////////////
/// Helper's destructor.
/// Check that no THtml object is attached to the helper - it might still need it!

THtml::THelperBase::~THelperBase()
{
   if (fHtml) {
      fHtml->HelperDeleted(this);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Set the THtml object owning this object; if it's already set to
/// a different THtml object than issue an error message and signal to
/// the currently set object that we are not belonging to it anymore.

void THtml::THelperBase::SetOwner(THtml* html) {
   if (fHtml && html && html != fHtml) {
      Error("SetOwner()", "Object already owned by an THtml instance!");
      fHtml->HelperDeleted(this);
   }
   fHtml = html;
}


////////////////////////////////////////////////////////////////////////////////
/// Set out_modulename to cl's module name; return true if it's valid.
/// If applicable, the module contains super modules separated by "/".
///
/// ROOT takes the directory part of cl's implementation file name
/// (or declaration file name, if the implementation file name is empty),
/// removes the last subdirectory if it is "src/" or "inc/", and interprets
/// the remaining path as the module hierarchy, converting it to upper case.
/// hist/histpainter/src/THistPainter.cxx thus becomes the module
/// HIST/HISTPAINTER. (Node: some ROOT packages get special treatment.)
/// If the file cannot be mapped into this scheme, the class's library
/// name (without directories, leading "lib" prefix or file extensions)
/// ius taken as the module name. If the module cannot be determined it is
/// set to "USER" and false is returned.
///
/// If your software cannot be mapped into this scheme then derive your
/// own class from TModuleDefinition and pass it to THtml::SetModuleDefinition().
///
/// The fse parameter is used to determine the relevant part of the path, i.e.
/// to not include parent directories of a TFileSysRoot.

bool THtml::TModuleDefinition::GetModule(TClass* cl, TFileSysEntry* fse,
                                         TString& out_modulename) const
{
   out_modulename = "USER";
   if (!cl) return false;

   // Filename: impl or decl?
   TString filename;
   if (fse) fse->GetFullName(filename, kFALSE);
   else {
      if (!GetOwner()->GetImplFileName(cl, kFALSE, filename))
         if (!GetOwner()->GetDeclFileName(cl, kFALSE, filename))
            return false;
   }
   TString inputdir = GetOwner()->GetInputPath();
   TString tok;
   Ssiz_t start = 0;
   // For -Idir/sub and A.h in dir/sub/A.h, use sub as module name if
   // it would eb empty otehrwise.
   TString trailingInclude;
   while (inputdir.Tokenize(tok, start, THtml::GetDirDelimiter())) {
      if (filename.BeginsWith(tok)) {
         if (tok.EndsWith("/") || tok.EndsWith("\\"))
            tok.Remove(tok.Length() - 1);
         trailingInclude = gSystem->BaseName(tok);
         filename.Remove(0, tok.Length());
         break;
      }
   }

   // take the directory name without "/" or leading "."
   out_modulename = gSystem->GetDirName(filename);

   while (out_modulename[0] == '.')
      out_modulename.Remove(0, 1);
   out_modulename.ReplaceAll("\\", "/");
   while (out_modulename[0] == '/')
      out_modulename.Remove(0, 1);
   while (out_modulename.EndsWith("/"))
      out_modulename.Remove(out_modulename.Length() - 1);

   if (!out_modulename[0])
      out_modulename = trailingInclude;

   if (!out_modulename[0])
      out_modulename = trailingInclude;

   // remove "/src", "/inc"
   if (out_modulename.EndsWith("/src")
      || out_modulename.EndsWith("/inc"))
      out_modulename.Remove(out_modulename.Length() - 4, 4);
   else {
   // remove "/src/whatever", "/inc/whatever"
      Ssiz_t pos = out_modulename.Index("/src/");
      if (pos == kNPOS)
         pos = out_modulename.Index("/inc/");
      if (pos != kNPOS)
         out_modulename.Remove(pos);
   }

   while (out_modulename.EndsWith("/"))
      out_modulename.Remove(out_modulename.Length() - 1);

   // special treatment:
   if (out_modulename == "MATH/GENVECTOR")
      out_modulename = "MATHCORE";
   else if (out_modulename == "MATH/MATRIX")
      out_modulename = "SMATRIX";
   else if (!out_modulename.Length()) {
      const char* cname= cl->GetName();
      if (strstr(cname, "::SMatrix<") || strstr(cname, "::SVector<"))
         out_modulename = "SMATRIX";
      else if (strstr(cname, "::TArrayProxy<") || strstr(cname, "::TClaArrayProxy<")
               || strstr(cname, "::TImpProxy<") || strstr(cname, "::TClaImpProxy<"))
         out_modulename = "TREEPLAYER";
      else {
         // determine the module name from the library name:
         out_modulename = cl->GetSharedLibs();
         Ssiz_t pos = out_modulename.Index(' ');
         if (pos != kNPOS)
            out_modulename.Remove(pos, out_modulename.Length());
         if (out_modulename.BeginsWith("lib"))
            out_modulename.Remove(0,3);
         pos = out_modulename.Index('.');
         if (pos != kNPOS)
            out_modulename.Remove(pos, out_modulename.Length());

         if (!out_modulename.Length()) {
            out_modulename = "USER";
            return false;
         }
      }
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Create all permutations of path and THtml's input path:
/// path being PP/ and THtml's input being .:include/:src/ gives
/// .:./PP/:include:include/PP/:src/:src/PP

void THtml::TFileDefinition::ExpandSearchPath(TString& path) const
{
   THtml* owner = GetOwner();
   if (!owner) return;

   TString pathext;
   TString inputdir = owner->GetInputPath();
   TString tok;
   Ssiz_t start = 0;
   while (inputdir.Tokenize(tok, start, THtml::GetDirDelimiter())) {
      if (pathext.Length())
         pathext += GetDirDelimiter();
      if (tok.EndsWith("\\"))
         tok.Remove(tok.Length() - 1);
      pathext += tok;
      if (path.BeginsWith(tok))
         pathext += GetDirDelimiter() + path;
      else
         pathext += GetDirDelimiter() + tok + "/" + path;
   }
   path = pathext;

}

////////////////////////////////////////////////////////////////////////////////
/// Given a class name with a scope, split the class name into directory part
/// and file name: A::B::C becomes module B, filename C.

void THtml::TFileDefinition::SplitClassIntoDirFile(const TString& clname, TString& dir,
                                                   TString& filename) const
{
   TString token;
   Ssiz_t from = 0;
   filename = "";
   dir = "";
   while (clname.Tokenize(token, from, "::") ) {
      dir = filename;
      filename = token;
   }

   // convert from Scope, class to module, filename.h
   dir.ToLower();
}


////////////////////////////////////////////////////////////////////////////////
/// Determine cl's declaration file name. Usually it's just
/// cl->GetDeclFileName(), but sometimes conversions need to be done
/// like include/ to abc/cde/inc/. If no declaration file name is
/// available, look for b/inc/C.h for class A::B::C. out_fsys will contain
/// the file system's (i.e. local machine's) full path name to the file.
/// The function returns false if the class's header file cannot be found.
///
/// If your software cannot be mapped into this scheme then derive your
/// own class from TFileDefinition and pass it to THtml::SetFileDefinition().

bool THtml::TFileDefinition::GetDeclFileName(const TClass* cl, TString& out_filename,
                                             TString& out_fsys, TFileSysEntry** fse) const
{
   return GetFileName(cl, true, out_filename, out_fsys, fse);
}

////////////////////////////////////////////////////////////////////////////////
/// Determine cl's implementation file name. Usually it's just
/// cl->GetImplFileName(), but sometimes conversions need to be done.
/// If no implementation file name is available look for b/src/C.cxx for
/// class A::B::C. out_fsys will contain the file system's (i.e. local
/// machine's) full path name to the file.
/// The function returns false if the class's source file cannot be found.
///
/// If your software cannot be mapped into this scheme then derive your
/// own class from TFileDefinition and pass it to THtml::SetFileDefinition().

bool THtml::TFileDefinition::GetImplFileName(const TClass* cl, TString& out_filename,
                                             TString& out_fsys, TFileSysEntry** fse) const
{
   return GetFileName(cl, false, out_filename, out_fsys, fse);
}


////////////////////////////////////////////////////////////////////////////////
/// Remove "/./" and collapse "/subdir/../" to "/"

void THtml::TFileDefinition::NormalizePath(TString& filename) const
{
   static const char* delim[] = {"/", "\\\\"};
   for (int i = 0; i < 2; ++i) {
      const char* d = delim[i];
      filename = filename.ReplaceAll(TString::Format("%c.%c", d[0], d[0]), TString(d[0]));
      TPRegexp reg(TString::Format("%s[^%s]+%s\\.\\.%s", d, d, d, d));
      while (reg.Substitute(filename, TString(d[0]), "", 0, 1)) {}
   }
   if (filename.BeginsWith("./") || filename.BeginsWith(".\\"))
      filename.Remove(0,2);
}


////////////////////////////////////////////////////////////////////////////////
/// Find filename in the list of system files; return the system file name
/// and change filename to the file name as included.
/// filename must be normalized (no "/./" etc) before calling.

TString THtml::TFileDefinition::MatchFileSysName(TString& filename, TFileSysEntry** fse) const
{
   const TList* bucket = GetOwner()->GetLocalFiles()->GetEntries().GetListForObject(gSystem->BaseName(filename));
   TString filesysname;
   if (bucket) {
      TIter iFS(bucket);
      TFileSysEntry* fsentry = 0;
      while ((fsentry = (TFileSysEntry*) iFS())) {
         if (!filename.EndsWith(fsentry->GetName()))
            continue;
         fsentry->GetFullName(filesysname, kTRUE); // get the short version
         filename = filesysname;
         if (!filename.EndsWith(filesysname)) {
            // It's something - let's see whether we find something better
            // else leave it as plan B. This helps finding Reflex sources.
            //filesysname = "";
            continue;
         }
         fsentry->GetFullName(filesysname, kFALSE); // get the long version
         if (fse) *fse = fsentry;
         break;
      }
   }
   return filesysname;
}


////////////////////////////////////////////////////////////////////////////////
/// Common implementation for GetDeclFileName(), GetImplFileName()

bool THtml::TFileDefinition::GetFileName(const TClass* cl, bool decl,
                                         TString& out_filename, TString& out_fsys,
                                         TFileSysEntry** fse) const
{
   out_fsys = "";

   if (!cl) {
      out_filename = "";
      return false;
   }

   TString possibleFileName;
   TString possiblePath;
   TString filesysname;

   TString clfile = decl ? cl->GetDeclFileName() : cl->GetImplFileName();
   NormalizePath(clfile);

   out_filename = clfile;
   if (clfile.Length()) {
      // check that clfile doesn't start with one of the include paths;
      // that's not what we want (include/TObject.h), we want the actual file
      // if it exists (core/base/inc/TObject.h)

      // special case for TMath namespace:
      if (clfile == "include/TMathBase.h") {
         clfile = "math/mathcore/inc/TMath.h";
         out_filename = clfile;
      }

      TString inclDir;
      TString inclPath(GetOwner()->GetPathInfo().fIncludePath);
      Ssiz_t pos = 0;
      Ssiz_t longestMatch = kNPOS;
      while (inclPath.Tokenize(inclDir, pos, GetOwner()->GetDirDelimiter())) {
         if (clfile.BeginsWith(inclDir) && (longestMatch == kNPOS || inclDir.Length() > longestMatch))
            longestMatch = inclDir.Length();
      }
      if (longestMatch != kNPOS) {
         clfile.Remove(0, longestMatch);
         if (clfile.BeginsWith("/") || clfile.BeginsWith("\\"))
            clfile.Remove(0, 1);
         TString asincl(clfile);
         GetOwner()->GetPathDefinition().GetFileNameFromInclude(asincl, clfile);
         out_filename = clfile;
      } else {
         // header file without a -Iinclude-dir prefix
         filesysname = MatchFileSysName(out_filename, fse);
         if (filesysname[0]) {
            clfile = out_filename;
         }
      }
   } else {
      // check for a file named like the class:
      filesysname = cl->GetName();
      int templateLevel = 0;
      Ssiz_t end = filesysname.Length();
      Ssiz_t start = end - 1;
      for (; start >= 0 && (templateLevel || filesysname[start] != ':'); --start) {
         if (filesysname[start] == '>')
            ++templateLevel;
         else if (filesysname[start] == '<') {
            --templateLevel;
            if (!templateLevel)
               end = start;
         }
      }
      filesysname = filesysname(start + 1, end - start - 1);
      if (decl)
         filesysname += ".h";
      else
         filesysname += ".cxx";
      out_filename = filesysname;
      filesysname = MatchFileSysName(out_filename, fse);
      if (filesysname[0]) {
         clfile = out_filename;
      }
   }

   if (!decl && !clfile.Length()) {
      // determine possible impl file name from the decl file name,
      // replacing ".whatever" by ".cxx", and looking for it in the known
      // file names
      TString declSysFileName;
      if (GetFileName(cl, true, filesysname, declSysFileName)) {
         filesysname = gSystem->BaseName(filesysname);
         Ssiz_t posExt = filesysname.Last('.');
         if (posExt != kNPOS)
            filesysname.Remove(posExt);
         filesysname += ".cxx";
         out_filename = filesysname;
         filesysname = MatchFileSysName(out_filename, fse);
         if (filesysname[0]) {
            clfile = out_filename;
         }
      }
   }

   if (clfile.Length() && !decl) {
      // Do not return the source file for these packages, even though we can find them.
      // THtml needs to have the class description in the source file if it finds the
      // source file, and these classes have their class descriptions in the header files.
      // THtml needs to be improved to collect all of a class' documentation before writing
      // it out, so it can take the class doc from the header even though a source exists.
      static const char* vetoClasses[] = {"math/mathcore/", "math/mathmore/", "math/genvector/",
                                          "math/minuit2/", "math/smatrix/"};
      for (unsigned int i = 0; i < sizeof(vetoClasses) / sizeof(char*); ++i) {
         if (clfile.Contains(vetoClasses[i])) {
            // of course there are exceptions from the exceptions:
            // TComplex and TRandom, TRandom1,...
            if (strcmp(cl->GetName(), "TComplex")
                && strcmp(cl->GetName(), "TMath")
                && strncmp(cl->GetName(), "TKDTree", 7)
                && strcmp(cl->GetName(), "TVirtualFitter")
                && strncmp(cl->GetName(), "TRandom", 7)) {
               out_filename = "";
               return false;
            } else break;
         }
      }
   }


   if (!clfile.Length()) {
      // determine possible decl file name from class + scope name:
      // A::B::C::myclass will result in possible file name myclass.h
      // in directory C/inc/
      out_filename = cl->GetName();
      if (!out_filename.Contains("::")) {
         out_filename = "";
         return false;
      }
      SplitClassIntoDirFile(out_filename, possiblePath, possibleFileName);

      // convert from Scope, class to module, filename.h
      if (possibleFileName.Length()) {
         if (decl)
            possibleFileName += ".h";
         else
            possibleFileName += ".cxx";
      }
      if (possiblePath.Length())
         possiblePath += "/";
      if (decl)
         possiblePath += "inc/";
      else
         possiblePath += "src/";
      out_filename = possiblePath + "/" + possibleFileName;
   } else {
      possiblePath = gSystem->GetDirName(clfile);
      possibleFileName = gSystem->BaseName(clfile);
   }

   if (possiblePath.Length())
      ExpandSearchPath(possiblePath);
   else possiblePath=".";

   out_fsys = gSystem->FindFile(possiblePath, possibleFileName, kReadPermission);
   if (out_fsys.Length()) {
      NormalizePath(out_fsys);
      return true;
   }
   out_filename = "";
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Determine the path to look for macros (see TDocMacroDirective) for
/// classes from a given module. If the path was sucessfully determined return true.
/// For ROOT, this directory is the "doc/macros" subdirectory of the module
/// directory; the path returned is GetDocDir(module) + "/macros".
///
/// If your software cannot be mapped into this scheme then derive your
/// own class from TPathDefinition and pass it to THtml::SetPathDefinition().

bool THtml::TPathDefinition::GetMacroPath(const TString& module, TString& out_dir) const
{
   TString moduledoc;
   if (!GetDocDir(module, moduledoc))
      return false;
   if (moduledoc.EndsWith("\\"))
      moduledoc.Remove(moduledoc.Length() - 1);

   TString macropath(GetOwner()->GetMacroPath());
   TString macrodirpart;
   out_dir = "";
   Ssiz_t pos = 0;
   while (macropath.Tokenize(macrodirpart, pos, ":")) {
      out_dir += moduledoc + "/" + macrodirpart + ":";
   }
   return true;
}


////////////////////////////////////////////////////////////////////////////////
/// Determine the module's documentation directory. If module is empty,
/// set doc_dir to the product's documentation directory.
/// If the path was sucessfuly determined return true.
/// For ROOT, this directory is the subdir "doc/" in the
/// module's path; the directory returned is module + "/doc".
///
/// If your software cannot be mapped into this scheme then derive your
/// own class from TPathDefinition and pass it to THtml::SetPathDefinition().

bool THtml::TPathDefinition::GetDocDir(const TString& module, TString& doc_dir) const
{
   doc_dir = "";
   if (GetOwner()->GetProductName() == "ROOT") {
      doc_dir = "$ROOTSYS";
      gSystem->ExpandPathName(doc_dir);
      doc_dir += "/";
   }

   if (module.Length())
      doc_dir += module + "/";
   doc_dir += GetOwner()->GetPathInfo().fDocPath;
   return true;
}


////////////////////////////////////////////////////////////////////////////////
/// Determine the path and filename used in an include statement for the
/// header file of the given class. E.g. the class ROOT::Math::Boost is
/// meant to be included as "Math/Genvector/Boost.h" - which is what
/// out_dir is set to. GetIncludeAs() returns whether the include
/// statement's path was successfully determined.
///
/// Any leading directory part that is part of fIncludePath (see SetIncludePath)
/// will be removed. For ROOT, leading "include/" is removed; everything after
/// is the include path.
///
/// If your software cannot be mapped into this scheme then derive your
/// own class from TPathDefinition and pass it to THtml::SetPathDefinition().

bool THtml::TPathDefinition::GetIncludeAs(TClass* cl, TString& out_dir) const
{
   out_dir = "";
   if (!cl || !GetOwner()) return false;

   TString hdr;
   if (!GetOwner()->GetDeclFileName(cl, kFALSE, hdr))
      return false;

   out_dir = hdr;
   bool includePathMatches = false;
   TString tok;
   Ssiz_t pos = 0;
   while (!includePathMatches && GetOwner()->GetPathInfo().fIncludePath.Tokenize(tok, pos, THtml::GetDirDelimiter()))
      if (out_dir.BeginsWith(tok)) {
         out_dir = hdr(tok.Length(), hdr.Length());
         if (out_dir[0] == '/' || out_dir[0] == '\\')
            out_dir.Remove(0, 1);
         includePathMatches = true;
      }

   if (!includePathMatches) {
      // We probably have a file super/module/inc/optional/filename.h.
      // That gets translated into optional/filename.h.
      // Assume that only one occurrence of "/inc/" exists in hdr.
      // If /inc/ is not part of the include file name then
      // just return the full path.
      // If we have matched any include path then this ROOT-only
      // algorithm is skipped!
      Ssiz_t posInc = hdr.Index("/inc/");
      if (posInc == kNPOS) return true;
      hdr.Remove(0, posInc + 5);
      out_dir = hdr;
   }

   return (out_dir.Length());
}


////////////////////////////////////////////////////////////////////////////////
/// Set out_fsname to the full pathname corresponding to a file
/// included as "included". Return false if this file cannot be determined
/// or found. For ROOT, out_fsname corresponds to included prepended with
/// "include"; only THtml prefers to work on the original files, e.g.
/// core/base/inc/TObject.h instead of include/TObject.h, so the
/// default implementation searches the TFileSysDB for an entry with
/// basename(included) and with matching directory part, setting out_fsname
/// to the TFileSysEntry's path.

bool THtml::TPathDefinition::GetFileNameFromInclude(const char* included, TString& out_fsname) const
{
   if (!included) return false;

   out_fsname = included;

   TString incBase(gSystem->BaseName(included));
   const TList* bucket = GetOwner()->GetLocalFiles()->GetEntries().GetListForObject(incBase);
   if (!bucket) return false;

   TString alldir = gSystem->GetDirName(included);
   TObjArray* arrSubDirs = alldir.Tokenize("/");
   TIter iEntry(bucket);
   TFileSysEntry* entry = 0;
   while ((entry = (TFileSysEntry*) iEntry())) {
      if (incBase != entry->GetName()) continue;
      // find entry with matching enclosing directory
      THtml::TFileSysDir* parent = entry->GetParent();
      for (int i = arrSubDirs->GetEntries() - 1; parent && i >= 0; --i) {
         const TString& subdir(((TObjString*)(*arrSubDirs)[i])->String());
         if (!subdir.Length() || subdir == ".")
            continue;
         if (subdir == parent->GetName())
            parent = parent->GetParent();
         else parent = 0;
      }
      if (parent) {
         // entry found!
         entry->GetFullName(out_fsname, kFALSE);
         delete arrSubDirs;
         return true;
      }
   }
   delete arrSubDirs;
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively fill entries by parsing the contents of path.

void THtml::TFileSysDir::Recurse(TFileSysDB* db, const char* path)
{
   TString dir(path);
   if (gDebug > 0 || GetLevel() < 2)
      Info("Recurse", "scanning %s...", path);
   TPMERegexp regexp(db->GetIgnore());
   dir += "/";
   void* hDir = gSystem->OpenDirectory(dir);
   const char* direntry = 0;
   while ((direntry = gSystem->GetDirEntry(hDir))) {
      if (!direntry[0] || direntry[0] == '.' || regexp.Match(direntry)) continue;
      TString entryPath(dir + direntry);
      if (gSystem->AccessPathName(entryPath, kReadPermission))
         continue;
      FileStat_t buf;
      if (!gSystem->GetPathInfo(entryPath, buf)) {
         if (R_ISDIR(buf.fMode)) {
            // skip if we would nest too deeply,  and skip soft links:
            if (GetLevel() > db->GetMaxLevel()
#ifndef R__WIN32
                || db->GetMapIno().GetValue(buf.fIno)
#endif
                ) continue;
            TFileSysDir* subdir = new TFileSysDir(direntry, this);
            fDirs.Add(subdir);
#ifndef R__WIN32
            db->GetMapIno().Add(buf.fIno, (Long_t)subdir);
#endif
            subdir->Recurse(db, entryPath);
         } else {
            int delen = strlen(direntry);
            // only .cxx and .h, .hxx are taken
            if (strcmp(direntry + delen - 4, ".cxx")
                && strcmp(direntry + delen - 2, ".h")
                && strcmp(direntry + delen - 4, ".hxx"))
               continue;
            TFileSysEntry* entry = new TFileSysEntry(direntry, this);
            db->GetEntries().Add(entry);
            fFiles.Add(entry);
         }
      } // if !gSystem->GetPathInfo()
   } // while dir entry
   gSystem->FreeDirectory(hDir);
}


////////////////////////////////////////////////////////////////////////////////
/// Recursively fill entries by parsing the path specified in GetName();
/// can be a THtml::GetDirDelimiter() delimited list of paths.

void THtml::TFileSysDB::Fill()
{
   TString dir;
   Ssiz_t posPath = 0;
   while (fName.Tokenize(dir, posPath, THtml::GetDirDelimiter())) {
      gSystem->ExpandPathName(dir);
      if (gSystem->AccessPathName(dir, kReadPermission)) {
         Warning("Fill", "Cannot read InputPath \"%s\"!", dir.Data());
         continue;
      }
      FileStat_t buf;
      if (!gSystem->GetPathInfo(dir, buf) && R_ISDIR(buf.fMode)) {
#ifndef R__WIN32
         TFileSysRoot* prevroot = (TFileSysRoot*) (Long_t)GetMapIno().GetValue(buf.fIno);
         if (prevroot != 0) {
            Warning("Fill", "InputPath \"%s\" already present as \"%s\"!", dir.Data(), prevroot->GetName());
            continue;
         }
#endif
         TFileSysRoot* root = new TFileSysRoot(dir, this);
         fDirs.Add(root);
#ifndef R__WIN32
         GetMapIno().Add(buf.fIno, (Long_t)root);
#endif
         root->Recurse(this, dir);
      } else {
         Warning("Fill", "Cannot read InputPath \"%s\"!", dir.Data());
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/* BEGIN_HTML
<p>The THtml class is designed to easily document
classes, code, and code related text files (like change logs). It generates HTML
pages conforming to the XHTML 1.0 transitional specifications; an example of
these pages is ROOT's own <a href="http://root.cern.ch/root/html/ClassIndex.html">
reference guide</a>. This page was verified to be valid XHTML 1.0 transitional,
which proves that all pages generated by THtml can be valid, as long as the user
provided XHTML (documentation, header, etc) is valid. You can check the current
THtml by clicking this icon:
<a href="http://validator.w3.org/check?uri=referer"><img
        src="http://www.w3.org/Icons/valid-xhtml10"
        alt="Valid XHTML 1.0 Transitional" height="31" width="88" style="border: none;"/></a></p>
Overview:
<ol style="list-style-type: upper-roman;">
  <li><a href="#usage">Usage</a></li>
  <li><a href="#conf">Configuration</a>
  <ol><li><a href="#conf:input">Input files</a></li>
  <li><a href="#conf:output">Output directory</a></li>
  <li><a href="#conf:liblink">Linking other documentation</a></li>
  <li><a href="#conf:classdoc">Recognizing class documentation</a></li>
  <li><a href="#conf:tags">Author, copyright, etc.</a></li>
  <li><a href="#conf:header">Header and footer</a></li>
  <li><a href="#conf:search">Links to searches, home page, ViewVC</a></li>
  <li><a href="#conf:charset">HTML Charset</a></li>
  </ol></li>
  <li><a href="#syntax">Documentation syntax</a>
  <ol><li><a href="#syntax:classdesc">Class description</a></li>
  <li><a href="#syntax:classidx">Class index</a></li>
  <li><a href="#syntax:meth">Method documentation</a></li>
  <li><a href="#syntax:datamem">Data member documentation</a></li>
  </ol></li>
  <li><a href="#directive">Documentation directives</a>
  <ol><li><a href="#directive:html"><tt>BEGIN<!-- -->_HTML</tt> <tt>END<!-- -->_HTML</tt>: include 'raw' HTML</a></li>
  <li><a href="#directive:macro"><tt>BEGIN<!-- -->_MACRO</tt> <tt>END<!-- -->_MACRO</tt>: include a picture generated by a macro</a></li>
  <li><a href="#directive:latex"><tt>BEGIN<!-- -->_LATEX</tt> <tt>END<!-- -->_LATEX</tt>: include a latex picture</a></li>
  </ol></li>
  <li><a href="#index">Product and module index</a></li>
  <li><a href="#aux">Auxiliary files: style sheet, JavaScript, help page</a></li>
  <li><a href="#charts">Class Charts</a></li>
  <li><a href="#confvar">Configuration variables</a></li>
  <li><a href="#how">Behind the scenes</a></li>
</ol>


<h3><a name="usage">I. Usage</a></h3>
These are typical things people do with THtml:
<pre>
    root[] <a href="http://root.cern.ch/root/html/THtml.html">THtml</a> html;                // create a <a href="http://root.cern.ch/root/html/THtml.html">THtml</a> object
    root[] html.LoadAllLibs();         // Load all rootmap'ed libraries
    root[] html.MakeAll();             // generate documentation for all changed classes
</pre>
or to run on just a few classes:
<pre>
    root[] <a href="http://root.cern.ch/root/html/THtml.html">THtml</a> html;                // create a <a href="http://root.cern.ch/root/html/THtml.html">THtml</a> object
    root[] html.MakeIndex();           // create auxiliary files (style sheet etc) and indices
    root[] html.MakeClass("TMyClass"); // create documentation for TMyClass only
</pre>
To "beautify" (i.e. create links to documentation for class names etc) some text
file or macro, use:
<pre>
    root[] html.Convert( "hsimple.C", "Histogram example" )
</pre>


<h3><a name="conf">II. Configuration</a></h3>
Most configuration options can be set as a call to THtml, or as a TEnv variable,
which you can set in your .rootrc.

<h4><a name="conf:input">II.1 Input files</a></h4>

<p>In your .rootrc, define Root.Html.SourceDir to point to directories containing
.cxx and .h files (see: <a href="http://root.cern.ch/root/html/TEnv.html">TEnv</a>)
of the classes you want to document, or call THtml::SetInputDir()</p>

<p>Example:</p><pre>
  Root.Html.SourceDir:  .:src:include
  Root.Html.Root:       http://root.cern.ch/root/html</pre>


<h4><a name="conf:output">II.2 Output directory</a></h4>

<p>The output directory can be specified using the Root.Html.OutputDir
configuration variable (default value: "htmldoc"). If that directory
doesn't exist <a href="http://root.cern.ch/root/html/THtml.html">THtml</a>
will create it.</p>

<p>Example:</p><pre>
  Root.Html.OutputDir:         htmldoc</pre>

<h4><a name="conf:liblink">II.3 Linking other documentation</a></h4>

<p>When trying to document a class, THtml searches for a source file in
the directories set via SetInputDir(). If it cannot find it, it assumes
that this class must have been documented before. Based on the library
this class is defined in, it checks the configuration variable
<tt>Root.Html.LibName</tt>, and creates a link using its value.
Alternatively, you can set these URLs via THtml::SetLibURL().</p>

<p>Example:<br/>
If a class MyClass is defined in class mylibs/libMyLib.so, and .rootrc
contains</p><pre>
  Root.Html.MyLib: ../mylib/</pre>
<p>THtml will create a link to "../mylib/MyClass.html".</p>

<p>The library name association can be set up using the rootmap facility.
For the library in the example above, which contains a dictionary
generated from the linkdef MyLinkdef.h, the command to generate the
rootmap file is</p>
<pre>  $ rlibmap -f -r rootmap -l mylib/libMyLib.so -d libCore.so -c MyLinkdef.h</pre>
<p>Here, <tt>-r</tt> specifies that the entries for libMyLib should be updated,
<tt>-l</tt> specifies the library we're dealing with, <tt>-d</tt> its
dependencies, and <tt>-c</tt> its linkdef. The rootmap file must be within
one of the <tt>LD_LIBRARY_PATH</tt> (or <tt>PATH</tt> for Windows) directories
when ROOT is started, otherwise ROOT will not use it.</p>

<h4><a name="conf:classdoc">II.4 Recognizing class documentation</a></h4>

<p>The class documentation has to appear in the header file containing the
class, right in front of its declaration. It is introduced by a string
defined by Root.Html.Description or SetClassDocTag(). See the section on
<a href="#syntax">documentation syntax</a> for further details.</p>

<p>Example:</p><pre>
  Root.Html.Description:       //____________________</pre>

<p>The class documentation will show which include statement is to be used
and which library needs to be linked to access it.
The include file name is determined via
<a href="http://root.cern.ch/root/html/TClass.html#TClass:GetDeclFileName">
TClass::GetDeclFileName()</a>;
leading parts are removed if they match any of the ':' separated entries in
THtml::GetIncludePath().</p>

<h4><a name="conf:tags">II.5 Author, copyright, etc.</a></h4>

<p>During the conversion,
<a href="http://root.cern.ch/root/html/THtml.html">THtml</a> will look for
some strings ("tags") in the source file, which have to appear right in
front of e.g. the author's name, copyright notice, etc. These tags can be
defined with the following environment variables: Root.Html.Author,
Root.Html.LastUpdate and Root.Html.Copyright, or with
SetAuthorTag(), SetLastUpdateTag(), SetCopyrightTag().</p>

<p>If the LastUpdate tag is not found, the current date and time are used.
This is useful when using
<a href="http://root.cern.ch/root/html/THtml.html#THtml:MakeAll">THtml::MakeAll()</a>'s
default option force=kFALSE, in which case
<a href="http://root.cern.ch/root/html/THtml.html">THtml</a> generates
documentation only for changed classes.</p>

Authors can be a comma separated list of author entries. Each entry has
one of the following two formats
<ul><li><tt>Name (non-alpha)</tt>.
<p><a href="http://root.cern.ch/root/html/THtml.html">THtml</a> will generate an
HTML link for <tt>Name</tt>, taking the Root.Html.XWho configuration
variable (defaults to "http://consult.cern.ch/xwho/people?") and adding
all parts of the name with spaces replaces by '+'. Non-alphanumerical
characters are printed out behind <tt>Name</tt>.</p>

<p>Example:</p>
<tt>// Author: Enrico Fermi</tt> appears in the source file.
<a href="http://root.cern.ch/root/html/THtml.html">THtml</a> will generate the link
<tt>http://consult.cern.ch/xwho/people?Enrico+Fermi</tt>. This works well for
people at CERN.</li>

<li><tt>Name &lt;link&gt; Info</tt>.
<p><a href="http://root.cern.ch/root/html/THtml.html">THtml</a> will generate
an HTML link for <tt>Name</tt> as specified by <tt>link</tt> and print
<tt>Info</tt> behind <tt>Name</tt>.</p>

<p>Example:</p>
<tt>// Author: Enrico Fermi &lt;http://www.enricos-home.it&gt;</tt> or<br/>
<tt>// Author: Enrico Fermi &lt;mailto:enrico@fnal.gov&gt;</tt> in the
source file. That's world compatible.</li>
</ul>

<p>Example (with defaults given):</p><pre>
      Root.Html.Author:     // Author:
      Root.Html.LastUpdate: // @(#)
      Root.Html.Copyright:  * Copyright
      Root.Html.XWho:       http://consult.cern.ch/xwho/people?</pre>


<h4><a name="conf:header">II.6 Header and footer</a></h4>

<p><a href="http://root.cern.ch/root/html/THtml.html">THtml</a> generates
a default header and footer for all pages. You can
specify your own versions with the configuration variables Root.Html.Header
and Root.Html.Footer, or by calling SetHeader(), SetFooter().
Both variables default to "", using the standard Root
versions. If it has a "+" appended, <a href="http://root.cern.ch/root/html/THtml.html">THtml</a> will
write both versions (user and root) to a file, for the header in the order
1st root, 2nd user, and for the footer 1st user, 2nd root (the root
versions containing "&lt;html&gt;" and &lt;/html&gt; tags, resp).</p>

<p>If you want to replace root's header you have to write a file containing
all HTML elements necessary starting with the &lt;doctype&gt; tag and ending with
(and including) the &lt;body&gt; tag. If you add your header it will be added
directly after Root's &lt;body&gt; tag. Any occurrence of the string <tt>%TITLE%</tt>
in the user's header file will be replaced by
a sensible, automatically generated title. If the header is generated for a
class, occurrences of <tt>%CLASS%</tt> will be replaced by the current class's name,
<tt>%SRCFILE%</tt> and <tt>%INCFILE%</tt> by the name of the source and header file, resp.
(as given by <a href="http://root.cern.ch/root/html/TClass.html#TClass:GetImplFileLine">TClass::GetImplFileName()</a>,
<a href="http://root.cern.ch/root/html/TClass.html#TClass:GetImplFileLine">TClass::GetDeclFileName()</a>).
If the header is not generated for a class, they will be replaced by "".</p>

<p>Root's footer starts with the tag &lt;!--SIGNATURE--&gt;. It includes the
author(s), last update, copyright, the links to the Root home page, to the
user home page, to the index file (ClassIndex.html), to the top of the page
and <tt>this page is automatically generated</tt> infomation. It ends with the
tags <tt>&lt;/body&gt;&lt;/html&gt;</tt>. If you want to replace it,
<a href="http://root.cern.ch/root/html/THtml.html">THtml</a> will search for some
tags in your footer: Occurrences of the strings <tt>%AUTHOR%</tt>, <tt>%UPDATE%</tt>, and
<tt>%COPYRIGHT%</tt> are replaced by their
corresponding values before writing the html file. The <tt>%AUTHOR%</tt> tag will be
replaced by the exact string that follows Root.Html.Author, no link
generation will occur.</p>


<h4><a name="conf:search">II.7 Links to searches, home page, ViewVC</a></h4>

<p>Additional parameters can be set by Root.Html.Homepage (address of the
user's home page), Root.Html.SearchEngine (search engine for the class
documentation), Root.Html.Search (search URL, where %u is replaced by the
referer and %s by the escaped search expression), and a ViewVC base URL
Root.Html.ViewCVS. For the latter, the file name is appended or, if
the URL contains %f, %f is replaced by the file name.
All values default to "".</p>

<p>Examples:</p><pre>
      Root.Html.Homepage:     http://www.enricos-home.it
      Root.Html.SearchEngine: http://root.cern.ch/root/Search.phtml
      Root.Html.Search:       http://www.google.com/search?q=%s+site%3A%u</pre>


<h4><a name="conf:charset">II.8 HTML Charset</a></h4>

<p>XHTML 1.0 transitional recommends the specification of the charset in the
content type meta tag, see e.g. <a href="http://www.w3.org/TR/2002/REC-xhtml1-20020801/">http://www.w3.org/TR/2002/REC-xhtml1-20020801/</a>
<a href="http://root.cern.ch/root/html/THtml.html">THtml</a> generates it for the HTML output files. It defaults to ISO-8859-1, and
can be changed using Root.Html.Charset.</p>

<p>Example:</p><pre>
      Root.Html.Charset:      EUC-JP</pre>

<h3><a name="syntax">III. Documentation syntax</a></h3>
<h4><a name="syntax:classdesc">III.1 Class description</a></h4>

<p>A class description block, which must be placed before the first
member function, has a following form:</p>
<pre>
////////////////////////////////////////////////////////////////
//                                                            //
// TMyClass                                                   //
//                                                            //
// This is the description block.                             //
//                                                            //
////////////////////////////////////////////////////////////////
</pre>
<p>The environment variable Root.Html.Description
(see: <a href="http://root.cern.ch/root/html/TEnv.html">TEnv</a>) contains
the delimiter string (default value: <tt>//_________________</tt>). It means
that you can also write your class description block like this:</p>
<pre>
   //_____________________________________________________________
   // A description of the class starts with the line above, and
   // will take place here !
   //
</pre>
<p>Note that <b><i>everything</i></b> until the first non-commented line is considered
as a valid class description block.</p>

<h4><a name="syntax:classidx">III.2 Class index</a></h4>

<p>All classes to be documented will have an entry in the ClassIndex.html,
showing their name with a link to their documentation page and a miniature
description. This discription for e.g. the class MyClass has to be given
in MyClass's header as a comment right after ClassDef(MyClass, n).</p>

<h4><a name="syntax:meth">III.3 Method documentation</a></h4>
<p>A member function description block starts immediately after '{'
and looks like this:</p>
<pre>
   void TWorld::HelloWorldFunc(string *text)
   {
      // This is an example of description for the
      // TWorld member function

      helloWorld.Print( text );
   }
</pre>
Like in a class description block, <b><i>everything</i></b> until the first
non-commented line is considered as a valid member function
description block.

If the rootrc variable <tt>Root.Html.DescriptionStyle</tt> is set to
<tt>Doc++</tt> THtml will also look for method documentation in front of
the function implementation. This feature is not recommended; source code
making use of this does not comply to the ROOT documentation standards, which
means future versions of THtml might not support it anymore.

<h4><a name="syntax:datamem">III.4 Data member documentation</a></h4>

<p>Data members are documented by putting a C++ comment behind their
declaration in the header file, e.g.</p>
<pre>
   int fIAmADataMember; // this is a data member
</pre>


<h3><a name="directive">IV. Documentation directives</a></h3>
<em>NOTE that THtml does not yet support nested directives
(i.e. latex inside html etc)!</em>

<h4><a name="directive:html">IV.1 <tt>BEGIN<!-- -->_HTML</tt> <tt>END<!-- -->_HTML</tt>: include 'raw' HTML</a></h4>

<p>You can insert pure html code into your documentation comments. During the
generation of the documentation, this code will be inserted as is
into the html file.</p>
<p>Pure html code must be surrounded by the keywords
<tt>BEGIN<!-- -->_HTML</tt> and <tt>END<!-- -->_HTML</tt>, where the
case is ignored.
An example of pure html code is this class description you are reading right now.
THtml uses a
<a href="http://root.cern.ch/root/html/TDocHtmlDirective.html">TDocHtmlDirective</a>
object to process this directive.</p>

<h4><a name="directive:macro">IV.2 <tt>BEGIN<!-- -->_MACRO</tt> <tt>END<!-- -->_MACRO</tt>: include a picture generated by a macro</a></h4>

<p>THtml can create images from scripts. You can either call an external
script by surrounding it by "begin_macro"/"end_macro", or include an unnamed
macro within these keywords. The macro should return a pointer to an object;
this object will then be saved as a GIF file.</p>
<p>Objects deriving from
<a href="http://root.cern.ch/root/html/TGObject.html">TGObject</a> (GUI elements)
will need to run in graphics mode (non-batch). You must specify this as a parameter:
"Begin_macro(GUI)...".
To create a second tab that displays the source of the macro you can specify
the argument "Begin_macro(source)...".
Of course you can combine them,
e.g. as "Begin_macro(source,gui)...".
THtml uses a
<a href="http://root.cern.ch/root/html/TDocMacroDirective.html">TDocMacroDirective</a>
object to process this directive.</p>
<p>This is an example:</p> END_HTML
BEGIN_MACRO(source)
{
  TCanvas* macro_example_canvas = new TCanvas("macro_example_canvas", "", 150, 150);
  macro_example_canvas->SetBorderSize(0);
  macro_example_canvas->SetFillStyle(1001);
  macro_example_canvas->SetFillColor(kWhite);
  macro_example_canvas->cd();
  TArc* macro_example_arc = new TArc(0.5,0.32,0.11,180,360);
  macro_example_arc->Draw();
  TEllipse* macro_example_ellipsis = new TEllipse(0.42,0.58,0.014,0.014,0,360,0);
  macro_example_ellipsis->SetFillStyle(0);
  macro_example_ellipsis->Draw();
  macro_example_ellipsis = new TEllipse(0.58,0.58,0.014,0.014,0,360,0);
  macro_example_ellipsis->SetFillStyle(0);
  macro_example_ellipsis->Draw();
  macro_example_ellipsis = new TEllipse(0.50,0.48,0.22,0.32,0,360,0);
  macro_example_ellipsis->SetFillStyle(0);
  macro_example_ellipsis->Draw();
  TLine* macro_example_line = new TLine(0.48,0.53,0.52,0.41);
  macro_example_line->Draw();
  return macro_example_canvas;
}
END_MACRO

BEGIN_HTML
<h4><a name="directive:latex">IV.3 <tt>BEGIN<!-- -->_LATEX</tt> <tt>END<!-- -->_LATEX</tt>: include a latex picture</a></h4>

<p>You can specify <a href="http://root.cern.ch/root/html/TLatex.html">TLatex</a>
style text and let THtml convert it into an image by surrounding it by "Begin_Latex", "End_Latex".
You can have multiple lines, and e.g. align each line at the '=' sign by passing
the argument <tt>separator='='</tt>. You can also specify how to align these parts;
if you want the part left of the separator to be right aligned, and the right part
to be left aligned, you could specify <tt>align='rl'</tt>.
THtml uses a <a href="http://root.cern.ch/root/html/TDocLatexDirective.html">TDocLatexDirective</a>
object to process the directive.
This is an example output with arguments <tt>separator='=', align='rl'</tt>:</p>
END_HTML BEGIN_LATEX(separator='=', align='rl')#kappa(x)^{2}=sin(x)^{x}
x=#chi^{2} END_LATEX

BEGIN_HTML

<h3><a name="index">V. Product and module index</a></h3>

<p><a href="#THtml:MakeIndex">THtml::MakeIndex()</a> will generate index files for classes
and types, all modules, and the product which you can set by
<a href="#THtml:SetProductName">THtml::SetProductName()</a>.
THtml will make use of external documentation in the module and product index,
either by linking it or by including it.
The files for modules are searched based on the source file directory of the
module's classes.</p>

<p>A filename starting with "index." will be included in the index page;
all other files will be linked.
Only files ending on <tt>.html</tt> or <tt>.txt</tt> will be taken into account;
the text files will first be run through
<a href="#THtml:Convert">THtml::Convert()</a>.
You can see an example <a href="http://root.cern.ch/root/html/HIST_Index.html">here</a>;
the part between "Index of HIST classes" and "Jump to" is created by parsing
the module's doc directory.</p>

<h3><a name="aux">VI. Auxiliary files: style sheet, JavaScript, help page</a></h3>

<p>The documentation pages share a common set of javascript and CSS files. They
are generated automatically when running <a href="#THtml:MakeAll">MakeAll()</a>;
they can be generated on
demand by calling <a href="#THtml:CreateAuxiliaryFiles">CreateAuxiliaryFiles()</a>.</p>


<h3><a name="charts">VII. Class Charts</a></h3>
THtml can generate a number of graphical representations for a class, which
are displayed as a tabbed set of imaged ontop of the class description.
It can show the inheritance, inherited and hidden members, directly and
indirectly included files, and library dependencies.

These graphs are generated using the <a href="http://www.graphviz.org/">Graphviz</a>
package. You can install it from <a href="http://www.graphviz.org">http://www.graphviz.org</a>.
You can either put it into your $PATH, or tell THtml where to find it by calling
<a href="#THtml:SetDotDir">SetDotDir()</a>.


<h3><a name="confvar">VIII. Configuration variables</a></h3>

<p>Here is a list of all configuration variables that are known to THtml.
You can set them in your .rootrc file, see
<a href="http://root.cern.ch/root/html/TEnv.html">TEnv</a>.</p>

<pre>
  Root.Html.OutputDir    (default: htmldoc)
  Root.Html.SourceDir    (default: .:src/:include/)
  Root.Html.Author       (default: // Author:) - start tag for authors
  Root.Html.LastUpdate   (default: // @(#)) - start tag for last update
  Root.Html.Copyright    (default:  * Copyright) - start tag for copyright notice
  Root.Html.Description  (default: //____________________ ) - start tag for class descr
  Root.Html.HomePage     (default: ) - URL to the user defined home page
  Root.Html.Header       (default: ) - location of user defined header
  Root.Html.Footer       (default: ) - location of user defined footer
  Root.Html.Root         (default: ) - URL of Root's class documentation
  Root.Html.SearchEngine (default: ) - link to the search engine
  Root.Html.Search       (defualt: ) - link to search by replacing "%s" with user input
  Root.Html.ViewCVS      (default: ) - URL of ViewCVS base
  Root.Html.XWho         (default: http://consult.cern.ch/xwho/people?) - URL of CERN's xWho
  Root.Html.Charset      (default: ISO-8859-1) - HTML character set
</pre>

<h3><a name="how">IX. Behind the scene</a></h3>

<p>Internally, THtml is just an API class that sets up the list of known
classes, and forwards API invocations to the "work horses".
<a href="http://root.cern.ch/root/html/TDocOutput.html">TDocOutput</a>
generates the output by letting a
<a href="http://root.cern.ch/root/html/TDocParser.html">TDocParser</a>
object parse the sources, which in turn invokes objects deriving from
<a href="http://root.cern.ch/root/html/TDocDirective.html">TDocDirective</a>
to process directives.</p>

END_HTML */
////////////////////////////////////////////////////////////////////////////////

ClassImp(THtml);
////////////////////////////////////////////////////////////////////////////////
/// Create a THtml object.
/// In case output directory does not exist an error
/// will be printed and gHtml stays 0 also zombie bit will be set.

THtml::THtml():
   fCounterFormat("%12s %5s %s"),
   fProductName("(UNKNOWN PRODUCT)"),
   fThreadedClassIter(0), fThreadedClassCount(0), fMakeClassMutex(0),
   fGClient(0), fPathDef(0), fModuleDef(0), fFileDef(0),
   fLocalFiles(0), fBatch(kFALSE)
{
   // check for source directory
   fPathInfo.fInputPath = gEnv->GetValue("Root.Html.SourceDir", "./:src/:include/");

   // check for output directory
   SetOutputDir(gEnv->GetValue("Root.Html.OutputDir", "htmldoc"));

   fLinkInfo.fXwho = gEnv->GetValue("Root.Html.XWho", "http://consult.cern.ch/xwho/people?");
   fLinkInfo.fROOTURL = gEnv->GetValue("Root.Html.Root", "http://root.cern.ch/root/html");
   fDocSyntax.fClassDocTag = gEnv->GetValue("Root.Html.Description", "//____________________");
   fDocSyntax.fAuthorTag = gEnv->GetValue("Root.Html.Author", "// Author:");
   fDocSyntax.fLastUpdateTag = gEnv->GetValue("Root.Html.LastUpdate", "// @(#)");
   fDocSyntax.fCopyrightTag = gEnv->GetValue("Root.Html.Copyright", "* Copyright");
   fOutputStyle.fHeader = gEnv->GetValue("Root.Html.Header", "");
   fOutputStyle.fFooter = gEnv->GetValue("Root.Html.Footer", "");
   fLinkInfo.fHomepage = gEnv->GetValue("Root.Html.Homepage", "");
   fLinkInfo.fSearchStemURL = gEnv->GetValue("Root.Html.Search", "");
   fLinkInfo.fSearchEngine = gEnv->GetValue("Root.Html.SearchEngine", "");
   fLinkInfo.fViewCVS = gEnv->GetValue("Root.Html.ViewCVS", "");
   fOutputStyle.fCharset = gEnv->GetValue("Root.Html.Charset", "ISO-8859-1");
   fDocSyntax.fDocStyle = gEnv->GetValue("Root.Html.DescriptionStyle", "");

   fDocEntityInfo.fClasses.SetOwner();
   fDocEntityInfo.fModules.SetOwner();
   // insert html object in the list of special ROOT objects
   if (!gHtml) {
      gHtml = this;
      gROOT->GetListOfSpecials()->Add(gHtml);
   }

}


////////////////////////////////////////////////////////////////////////////////
/// Default destructor

THtml::~THtml()
{
   fDocEntityInfo.fClasses.Clear();
   fDocEntityInfo.fModules.Clear();
   if (gHtml == this) {
      gROOT->GetListOfSpecials()->Remove(gHtml);
      gHtml = 0;
   }
   delete fPathDef;
   delete fModuleDef;
   delete fFileDef;
   delete fLocalFiles;
}

////////////////////////////////////////////////////////////////////////////////
/// Add path to the directories to be searched for macro files
/// that are to be executed via the TDocMacroDirective
/// ("Begin_Macro"/"End_Macro"); relative to the source file
/// that the directive is run on.

void THtml::AddMacroPath(const char* path)
{
   const char pathDelimiter =
#ifdef R__WIN32
      ';';
#else
      ':';
#endif
   fPathInfo.fMacroPath += pathDelimiter;
   fPathInfo.fMacroPath += path;
}


////////////////////////////////////////////////////////////////////////////////
/// copy CSS, javascript file, etc to the output dir

void THtml::CreateAuxiliaryFiles() const
{
   CreateJavascript();
   CreateStyleSheet();
   CopyFileFromEtcDir("HELP.html");
}

////////////////////////////////////////////////////////////////////////////////
/// Return the TModuleDefinition (or derived) object as set by
/// SetModuleDefinition(); create and return a TModuleDefinition object
/// if none was set.

const THtml::TModuleDefinition& THtml::GetModuleDefinition() const
{
   if (!fModuleDef) {
      fModuleDef = new TModuleDefinition();
      fModuleDef->SetOwner(const_cast<THtml*>(this));
   }
   return *fModuleDef;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the TFileDefinition (or derived) object as set by
/// SetFileDefinition(); create and return a TFileDefinition object
/// if none was set.

const THtml::TFileDefinition& THtml::GetFileDefinition() const
{
   if (!fFileDef) {
      fFileDef = new TFileDefinition();
      fFileDef->SetOwner(const_cast<THtml*>(this));
   }
   return *fFileDef;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the TModuleDefinition (or derived) object as set by
/// SetModuleDefinition(); create and return a TModuleDefinition object
/// if none was set.

const THtml::TPathDefinition& THtml::GetPathDefinition() const
{
   if (!fPathDef) {
      fPathDef = new TPathDefinition();
      fPathDef->SetOwner(const_cast<THtml*>(this));
   }
   return *fPathDef;
}


////////////////////////////////////////////////////////////////////////////////
/// Get the directory containing THtml's auxiliary files ($ROOTSYS/etc/html)

const char* THtml::GetEtcDir() const
{
   if (fPathInfo.fEtcDir.Length())
      return fPathInfo.fEtcDir;

   R__LOCKGUARD(GetMakeClassMutex());

   fPathInfo.fEtcDir = "html";
   gSystem->PrependPathName(TROOT::GetEtcDir(), fPathInfo.fEtcDir);

   return fPathInfo.fEtcDir;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the next class to be generated for MakeClassThreaded.

TClassDocInfo *THtml::GetNextClass()
{
   if (!fThreadedClassIter) return 0;

   R__LOCKGUARD(GetMakeClassMutex());

   TClassDocInfo* classinfo = 0;
   while ((classinfo = (TClassDocInfo*)(*fThreadedClassIter)())
          && !classinfo->IsSelected()) { }

   if (!classinfo) {
      delete fThreadedClassIter;
      fThreadedClassIter = 0;
   }

   fCounter.Form("%5d", fDocEntityInfo.fClasses.GetSize() - fThreadedClassCount++);

   return classinfo;
}


////////////////////////////////////////////////////////////////////////////////
/// Get the documentation URL for library lib.
/// If lib == 0 or no documentation URL has been set for lib, return the ROOT
/// documentation URL. The return value is always != 0.

const char* THtml::GetURL(const char* lib /*=0*/) const
{
   R__LOCKGUARD(GetMakeClassMutex());

   if (lib && strlen(lib)) {
      std::map<std::string, TString>::const_iterator iUrl = fLinkInfo.fLibURLs.find(lib);
      if (iUrl != fLinkInfo.fLibURLs.end()) return iUrl->second;
      return gEnv->GetValue(TString("Root.Html.") + lib, fLinkInfo.fROOTURL);
   }
   return fLinkInfo.fROOTURL;
}

////////////////////////////////////////////////////////////////////////////////
/// Check whether dot is available in $PATH or in the directory set
/// by SetDotPath()

Bool_t THtml::HaveDot()
{
   if (fPathInfo.fFoundDot != PathInfo_t::kDotUnknown)
      return (fPathInfo.fFoundDot == PathInfo_t::kDotFound);

   R__LOCKGUARD(GetMakeClassMutex());

   Info("HaveDot", "Checking for Graphviz (dot)...");
   TString runDot("dot");
   if (fPathInfo.fDotDir.Length())
      gSystem->PrependPathName(fPathInfo.fDotDir, runDot);
   runDot += " -V";
   if (gDebug > 3)
      Info("HaveDot", "Running: %s", runDot.Data());
   if (gSystem->Exec(runDot)) {
      fPathInfo.fFoundDot = PathInfo_t::kDotNotFound;
      return kFALSE;
   }
   fPathInfo.fFoundDot = PathInfo_t::kDotFound;
   return kTRUE;

}

////////////////////////////////////////////////////////////////////////////////
/// Inform the THtml object that one of its helper objects was deleted.
/// Called by THtml::HelperBase::~HelperBase().

void THtml::HelperDeleted(THtml::THelperBase* who)
{
   THelperBase* helpers[3] = {fPathDef, fModuleDef, fFileDef};
   for (int i = 0; who && i < 3; ++i)
      if (who == helpers[i])
         helpers[i] = who = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// It converts a single text file to HTML
///
///
/// Input: filename - name of the file to convert
///        title    - title which will be placed at the top of the HTML file
///        dirname  - optional parameter, if it's not specified, output will
///                   be placed in htmldoc/examples directory.
///        relpath  - optional parameter pointing to the THtml generated doc
///                   on the server, relative to the current page.
///        includeOutput - if != kNoOutput, run the script passed as filename and
///                   store all created canvases in PNG files that are
///                   shown next to the converted source. Bitwise-ORing with
///                   kForceOutput re-runs the script even if output PNGs exist
///                   that are newer than the script. If kCompiledOutput is
///                   passed, the script is run through ACLiC (.x filename+)
///        context  - line shown verbatim at the top of the page; e.g. for links.
///                   If context is non-empty it is expected to also print the
///                   title.
///
///  NOTE: Output file name is the same as filename, but with extension .html
///

void THtml::Convert(const char *filename, const char *title,
                    const char *dirname /*= ""*/, const char *relpath /*= "../"*/,
                    Int_t includeOutput /* = kNoOutput */,
                    const char* context /* = "" */)
{
   gROOT->GetListOfGlobals(kTRUE);        // force update of this list
   CreateListOfClasses("*");

   const char *dir;
   TString dfltdir;

   // if it's not defined, make the "examples" as a default directory
   if (!*dirname) {
      gSystem->ExpandPathName(fPathInfo.fOutputDir);
      char *tmp0 = gSystem->ConcatFileName(fPathInfo.fOutputDir, "examples");
      dfltdir = tmp0;
      delete [] tmp0;
      dir = dfltdir.Data();
   } else {
      dir = dirname;
   }

   // create directory if necessary
   if (gSystem->AccessPathName(dir))
      gSystem->MakeDirectory(dir);

   // find a file
   char *cRealFilename =
       gSystem->Which(fPathInfo.fInputPath, filename, kReadPermission);

   if (!cRealFilename) {
      Error("Convert", "Can't find file '%s' !", filename);
      return;
   }

   TString realFilename = cRealFilename;
   delete[] cRealFilename;

   // open source file
   std::ifstream sourceFile;
   sourceFile.open(realFilename, std::ios::in);

   if (!sourceFile.good()) {
      Error("Convert", "Can't open file '%s' !", realFilename.Data());
      return;
   }

   if (gSystem->AccessPathName(dir)) {
      Error("Convert",
            "Directory '%s' doesn't exist, or it's write protected !", dir);
      return;
   }
   char *tmp1 =
       gSystem->ConcatFileName(dir, gSystem->BaseName(filename));

   TDocOutput output(*this);
   if (!fGClient)
      gROOT->ProcessLine(TString::Format("*((TGClient**)0x%lx) = gClient;",
                                         (ULong_t)&fGClient));
   if (includeOutput && !fGClient)
      Warning("Convert", "Output requested but cannot initialize graphics: GUI  and GL windows not be available");
   output.Convert(sourceFile, realFilename, tmp1, title, relpath, includeOutput, context, fGClient);

   delete [] tmp1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the module name for a given class.
/// Use the cached information from fDocEntityInfo.fClasses.

void  THtml::GetModuleNameForClass(TString& module, TClass* cl) const
{
   module = "(UNKNOWN)";
   if (!cl) return;

   TClassDocInfo* cdi = (TClassDocInfo*)fDocEntityInfo.fClasses.FindObject(cl->GetName());
   if (!cdi || !cdi->GetModule())
      return;
   module = cdi->GetModule()->GetName();
}


////////////////////////////////////////////////////////////////////////////////
/// Create the list of all known classes

void THtml::CreateListOfClasses(const char* filter)
{
   if (fDocEntityInfo.fClasses.GetSize() && fDocEntityInfo.fClassFilter == filter)
      return;

   Info("CreateListOfClasses", "Initializing - this might take a while...");
   // get total number of classes
   Int_t totalNumberOfClasses = gClassTable->Classes();

   // allocate memory
   fDocEntityInfo.fClasses.Clear();
   fDocEntityInfo.fModules.Clear();

   fDocEntityInfo.fClassFilter = filter;

   // start from beginning
   gClassTable->Init();
   if (filter && (!filter[0] || !strcmp(filter, "*")))
      filter = ".*";
   TString reg = filter;
   TPMERegexp re(reg);

   bool skipROOTClasses = false;
   std::set<std::string> rootLibs;
   TList classesDeclFileNotFound;
   TList classesImplFileNotFound;

   // pre-run TObject at i == -1
   for (Int_t i = -1; i < totalNumberOfClasses; i++) {

      // get class name
      const char *cname = 0;
      if (i < 0) cname = "TObject";
      else cname = gClassTable->Next();
      if (!cname)
         continue;

      if (i >= 0 && !strcmp(cname, "TObject")) {
         // skip the second iteration on TObject
         continue;
      }

      // This is a hack for until after Cint and Reflex are one.
      if (strstr(cname, "__gnu_cxx::")) continue;
      // Work around ROOT-6016
      if (!strcmp(cname, "timespec")) continue;
      // "tuple"s are synthetic in the interpreter
      if (!strncmp(cname, "tuple<", 6)) continue;

      // get class & filename - use TROOT::GetClass, as we also
      // want those classes without decl file name!
      TClass *classPtr = TClass::GetClass((const char *) cname, kTRUE);
      if (!classPtr) continue;

      std::string shortName(ShortType(cname));
      cname = shortName.c_str();

      TString s = cname;
      Bool_t matchesSelection = re.Match(s);


      TString hdr;
      TString hdrFS;
      TString src;
      TString srcFS;
      TString htmlfilename;
      TFileSysEntry* fse = 0;

      TClassDocInfo* cdi = (TClassDocInfo*) fDocEntityInfo.fClasses.FindObject(cname);
      if (cdi) {
         hdr = cdi->GetDeclFileName();
         hdrFS = cdi->GetDeclFileSysName();
         src = cdi->GetImplFileName();
         srcFS = cdi->GetImplFileSysName();
         htmlfilename = cdi->GetHtmlFileName();
      }

      if (!hdrFS.Length()) {
         if (!GetFileDefinition().GetDeclFileName(classPtr, hdr, hdrFS, &fse)) {
            // we don't even know where the class is defined;
            // just skip. Silence if it doesn't match the selection anyway
            if (i == -1 ) {
               skipROOTClasses = true;
               Info("CreateListOfClasses", "Cannot find header file for TObject at %s given the input path %s.",
                  classPtr->GetDeclFileName(), GetInputPath().Data());
               Info("CreateListOfClasses", "Assuming documentation is not for ROOT classes, or you need to pass "
                  "the proper directory to THtml::SetInputDir() so I can find %s.", classPtr->GetDeclFileName());
               continue;
            }
            // ignore STL
            if (classPtr->GetClassInfo() &&
                (gInterpreter->ClassInfo_Property(classPtr->GetClassInfo()) & kIsDefinedInStd))
              continue;
            if (classPtr->GetDeclFileName() && (!strncmp(classPtr->GetDeclFileName(), "prec_stl/", 9) ||
                                                strstr(classPtr->GetDeclFileName(), "include/c++/") ||
                                                !strncmp(classPtr->GetDeclFileName(), "/usr/include",12)))
               continue;
            if (classPtr->GetDeclFileName() && (
                 !strcmp(classPtr->GetDeclFileName(), "vector") ||
                 !strcmp(classPtr->GetDeclFileName(), "string") ||
                 !strcmp(classPtr->GetDeclFileName(), "list") ||
                 !strcmp(classPtr->GetDeclFileName(), "deque") ||
                 !strcmp(classPtr->GetDeclFileName(), "map") ||
                 !strcmp(classPtr->GetDeclFileName(), "valarray") ||
                 !strcmp(classPtr->GetDeclFileName(), "set") ||
                 !strcmp(classPtr->GetDeclFileName(), "typeinfo") ||
                 !strcmp(classPtr->GetDeclFileName(), "stdlib.h") ) )
            {
               // Those are STL header, just ignore.
               continue;
            }
            if (skipROOTClasses) {
               if (classPtr->GetSharedLibs() && classPtr->GetSharedLibs()[0]) {
                  std::string lib(classPtr->GetSharedLibs());
                  size_t posSpace = lib.find(' ');
                  if (posSpace != std::string::npos)
                     lib.erase(posSpace);
                  if (rootLibs.find(lib) == rootLibs.end()) {
                     TString rootlibdir = TROOT::GetLibDir();
                     TString sLib(lib);
                     if (sLib.Index('.') == -1) {
                        sLib += ".";
                        sLib += gSystem->GetSoExt();
                     }
                     gSystem->PrependPathName(rootlibdir, sLib);
                     if (gSystem->AccessPathName(sLib))
                        // the library doesn't exist in $ROOTSYS/lib, so it's not
                        // a root lib and we need to tell the user.
                        classesDeclFileNotFound.AddLast(classPtr);
                     else rootLibs.insert(lib);
                  } // end "if rootLibs does not contain lib"
               } else {
                  // lib name unknown
                  static const char* rootClassesToIgnore[] =
                  { "ColorStruct_t", "CpuInfo_t", "Event_t", "FileStat_t", "GCValues_t", "MemInfo_t",
                     "PictureAttributes_t", "Point_t", "ProcInfo_t", "ROOT", "ROOT::Fit",
                     "Rectangle_t", "RedirectHandle_t", "Segment_t", "SetWindowAttributes_t",
                     "SysInfo_t", "TCint", "UserGroup_t", "WindowAttributes_t", "timespec", 0};
                  static const char* rootClassStemsToIgnore[] =
                  { "ROOT::Math", "TKDTree", "TMatrixT", "TParameter", "vector", 0 };
                  static size_t rootClassStemsToIgnoreLen[] = {0, 0, 0, 0, 0};
                  static std::set<std::string> setRootClassesToIgnore;
                  if (setRootClassesToIgnore.empty()) {
                     for (int ii = 0; rootClassesToIgnore[ii]; ++ii)
                        setRootClassesToIgnore.insert(rootClassesToIgnore[ii]);
                     for (int ii = 0; rootClassStemsToIgnore[ii]; ++ii)
                        rootClassStemsToIgnoreLen[ii] = strlen(rootClassStemsToIgnore[ii]);
                  }
                  // only complain about this class if it should not be ignored:
                  if (setRootClassesToIgnore.find(cname) == setRootClassesToIgnore.end()) {
                     bool matched = false;
                     for (int ii = 0; !matched && rootClassStemsToIgnore[ii]; ++ii)
                        matched = !strncmp(cname, rootClassStemsToIgnore[ii], rootClassStemsToIgnoreLen[ii]);
                     if (!matched)
                        classesDeclFileNotFound.AddLast(classPtr);
                  }
               } // lib name known
               continue;
            } else {
               if (matchesSelection && (!classPtr->GetDeclFileName() ||
                                        !strstr(classPtr->GetDeclFileName(),"prec_stl/") ||
                                        !strstr(classPtr->GetDeclFileName(), "include/c++/") ||
                                        strncmp(classPtr->GetDeclFileName(), "/usr/include",12)))
                  classesDeclFileNotFound.AddLast(classPtr);
               continue;
            }
         }
      }

      Bool_t haveSource = (srcFS.Length());
      if (!haveSource)
         haveSource = GetFileDefinition().GetImplFileName(classPtr, src, srcFS, fse ? 0 : &fse);

      if (!haveSource) {
         classesImplFileNotFound.AddLast(classPtr);
      }

      if (!htmlfilename.Length())
         GetHtmlFileName(classPtr, htmlfilename);

      if (!cdi) {
         cdi = new TClassDocInfo(classPtr, htmlfilename, hdrFS, srcFS, hdr, src);
         fDocEntityInfo.fClasses.Add(cdi);
      } else {
         cdi->SetDeclFileName(hdr);
         cdi->SetImplFileName(src);
         cdi->SetDeclFileSysName(hdrFS);
         cdi->SetImplFileSysName(srcFS);
         cdi->SetHtmlFileName(htmlfilename);
      }

      cdi->SetSelected(matchesSelection);

      TString modulename;
      GetModuleDefinition().GetModule(classPtr, fse, modulename);
      if (!modulename.Length() || modulename == "USER")
         GetModuleNameForClass(modulename, classPtr);

      TModuleDocInfo* module = (TModuleDocInfo*) fDocEntityInfo.fModules.FindObject(modulename);
      if (!module) {
         bool moduleSelected = cdi->IsSelected();

         TString parentModuleName = gSystem->GetDirName(modulename);
         TModuleDocInfo* super = nullptr;
         if (parentModuleName.Length() && parentModuleName != ".") {
            super = (TModuleDocInfo*) fDocEntityInfo.fModules.FindObject(parentModuleName);
            if (!super) {
               // create parents:
               TString token;
               Ssiz_t pos = 0;
               while (parentModuleName.Tokenize(token, pos, "/")) {
                  if (!token.Length() || token == ".") continue;
                  super = new TModuleDocInfo(token, super);
                  super->SetSelected(moduleSelected);
                  fDocEntityInfo.fModules.Add(super);
               }
            }
         }
         module = new TModuleDocInfo(modulename, super);
         module->SetSelected(moduleSelected);
         fDocEntityInfo.fModules.Add(module);
      }

      if (module) {
         module->AddClass(cdi);
         cdi->SetModule(module);
         if (cdi->HaveSource() && cdi->IsSelected())
            module->SetSelected();
      }

      // clear the typedefs; we fill them later
      cdi->GetListOfTypedefs().Clear();

      if (gDebug > 0)
         Info("CreateListOfClasses", "Adding class %s, module %s (%sselected)",
              cdi->GetName(), module ? module->GetName() : "[UNKNOWN]",
              cdi->IsSelected() ? "" : "not ");
   }



   bool cannotFind = false;
   if (!classesDeclFileNotFound.IsEmpty()) {
      Warning("CreateListOfClasses",
         "Cannot find the header for the following classes [reason]:");
      TIter iClassesDeclFileNotFound(&classesDeclFileNotFound);
      TClass* iClass = 0;
      while ((iClass = (TClass*)iClassesDeclFileNotFound())) {
         if (iClass->GetDeclFileName() && iClass->GetDeclFileName()[0]) {
            Warning("CreateListOfClasses", "   %s [header %s not found]", iClass->GetName(), iClass->GetDeclFileName());
            cannotFind = true;
         } else
            Warning("CreateListOfClasses", "   %s [header file is unknown]", iClass->GetName());
      }
   }

   if (!classesImplFileNotFound.IsEmpty() && gDebug > 3) {
      Warning("CreateListOfClasses",
         "Cannot find the source file for the following classes [reason]:");
      TIter iClassesDeclFileNotFound(&classesImplFileNotFound);
      TClass* iClass = 0;
      while ((iClass = (TClass*)iClassesDeclFileNotFound())) {
         if (iClass->GetDeclFileName() && iClass->GetDeclFileName()[0]) {
            Info("CreateListOfClasses", "   %s [source %s not found]", iClass->GetName(), iClass->GetImplFileName());
            cannotFind = true;
         } else
            Info("CreateListOfClasses", "   %s [source file is unknown, add \"ClassImpl(%s)\" to source file if it exists]",
               iClass->GetName(), iClass->GetName());
      }
   }
   if (cannotFind) {
      Warning("CreateListOfClasses", "THtml cannot find all headers and sources. ");
      Warning("CreateListOfClasses",
         "You might need to adjust the input path (currently %s) by calling THtml::SetInputDir()",
         GetInputPath().Data());
   }

   // fill typedefs
   TIter iTypedef(gROOT->GetListOfTypes());
   TDataType* dt = 0;
   TDocOutput output(*this);
   while ((dt = (TDataType*) iTypedef())) {
      if (dt->GetType() != -1) continue;
      TClassDocInfo* cdi = (TClassDocInfo*) fDocEntityInfo.fClasses.FindObject(dt->GetFullTypeName());
      if (cdi) {
         cdi->GetListOfTypedefs().Add(dt);
         if (gDebug > 1)
            Info("CreateListOfClasses", "Adding typedef %s to class %s",
                 dt->GetName(), cdi->GetName());

         bool inNamespace = true;
         TString surroundingNamespace(dt->GetName());
         Ssiz_t posTemplate = surroundingNamespace.Last('>');
         inNamespace = inNamespace && (posTemplate == kNPOS);
         if (inNamespace) {
            Ssiz_t posColumn = surroundingNamespace.Last(':');
            if (posColumn != kNPOS) {
               surroundingNamespace.Remove(posColumn - 1);
               TClass* clSurrounding = GetClass(surroundingNamespace);
               inNamespace = inNamespace && (!clSurrounding || IsNamespace(clSurrounding));
            }
         }
         if (inNamespace && cdi->GetModule()) {
            TString htmlfilename(dt->GetName());
            output.NameSpace2FileName(htmlfilename);
            htmlfilename += ".html";
            TClassDocInfo* cdiTD = new TClassDocInfo(dt, htmlfilename);
            cdiTD->SetModule(cdi->GetModule());
            cdiTD->SetSelected(cdi->IsSelected());
            cdi->GetModule()->AddClass(cdiTD);
         }
      }
   }

   fDocEntityInfo.fClasses.Sort();
   fDocEntityInfo.fModules.Sort();
   TIter iterModule(&fDocEntityInfo.fModules);
   TModuleDocInfo* mdi = 0;
   while ((mdi = (TModuleDocInfo*) iterModule()))
      mdi->GetClasses()->Sort();

   if (fProductName == "(UNKNOWN PRODUCT)"
      && fDocEntityInfo.fModules.FindObject("core/base")
      && fDocEntityInfo.fModules.FindObject("core/cont")
      && fDocEntityInfo.fModules.FindObject("core/rint")
      && gProgName && strstr(gProgName, "root"))
      // if we have these modules we're probably building the root doc
      fProductName = "ROOT";

   if (fProductName == "(UNKNOWN PRODUCT)") {
      Warning("CreateListOfClasses", "Product not set. You should call gHtml->SetProductName(\"MyProductName\");");
   } else if (fProductName != "ROOT") {
      if (GetViewCVS().Contains("http://root.cern.ch/"))
         SetViewCVS("");
   }

   if (fDocEntityInfo.fModules.GetEntries() == 1
      && fDocEntityInfo.fModules.At(0)->GetName()
      && !strcmp(fDocEntityInfo.fModules.At(0)->GetName(), "(UNKNOWN)"))
      // Only one module, and its name is not known.
      // Let's call it "MAIN":
      ((TModuleDocInfo*) fDocEntityInfo.fModules.At(0))->SetName("MAIN");

   Info("CreateListOfClasses", "Initializing - DONE.");
}


////////////////////////////////////////////////////////////////////////////////
/// Create index of all data types and a page for each typedef-to-class

void THtml::CreateListOfTypes()
{
   TDocOutput output(*this);
   output.CreateTypeIndex();
   output.CreateClassTypeDefs();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a file from $ROOTSYS/etc/html into GetOutputDir()

Bool_t THtml::CopyFileFromEtcDir(const char* filename) const {
   R__LOCKGUARD(GetMakeClassMutex());

   TString outFile(filename);

   TString inFile(outFile);
   gSystem->PrependPathName(GetEtcDir(), inFile);

   gSystem->PrependPathName(GetOutputDir(), outFile);

   if (gSystem->CopyFile(inFile, outFile, kTRUE) != 0) {
      Warning("CopyFileFromEtcDir", "Could not copy %s to %s", inFile.Data(), outFile.Data());
      return kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create the inheritance hierarchy diagram for all classes

void THtml::CreateHierarchy()
{
   TDocOutput output(*this);
   output.CreateHierarchy();
}

////////////////////////////////////////////////////////////////////////////////
/// Write the default ROOT style sheet.

void THtml::CreateJavascript() const {
   CopyFileFromEtcDir("ROOT.js");
}

////////////////////////////////////////////////////////////////////////////////
/// Write the default ROOT style sheet.

void THtml::CreateStyleSheet() const {
   CopyFileFromEtcDir("ROOT.css");
   CopyFileFromEtcDir("shadowAlpha.png");
   CopyFileFromEtcDir("shadow.gif");
}



////////////////////////////////////////////////////////////////////////////////
/// fill derived with all classes inheriting from cl and their inheritance
/// distance to cl

void THtml::GetDerivedClasses(TClass* cl, std::map<TClass*, Int_t>& derived) const
{
   TIter iClass(&fDocEntityInfo.fClasses);
   TClassDocInfo* cdi = 0;
   while ((cdi = (TClassDocInfo*) iClass())) {
      TClass* candidate = dynamic_cast<TClass*>(cdi->GetClass());
      if (!candidate) continue;
      if (candidate != cl && candidate->InheritsFrom(cl)) {
         Int_t level = 0;
         TClass* currentBaseOfCandidate = candidate;
         while (currentBaseOfCandidate != cl) {
            TList* bases = currentBaseOfCandidate->GetListOfBases();
            if (!bases) continue;
            TIter iBase(bases);
            TBaseClass* base = 0;
            while ((base = (TBaseClass*) iBase())) {
               TClass* clBase = base->GetClassPointer();
               if (clBase && clBase->InheritsFrom(cl)) {
                  ++level;
                  currentBaseOfCandidate = clBase;
               }
            }
         }
         derived[candidate] = level;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return real HTML filename
///
///
///  Input: classPtr - pointer to a class
///         filename - string containing a full name
///         of the corresponding HTML file after the function returns.
///

void THtml::GetHtmlFileName(TClass * classPtr, TString& filename) const
{
   filename.Remove(0);
   if (!classPtr) return;

   TString cFilename;
   if (!GetImplFileName(classPtr, kFALSE, cFilename))
      GetDeclFileName(classPtr, kFALSE, cFilename);

   // classes without Impl/DeclFileName don't have docs,
   // and classes without docs don't have output file names
   if (!cFilename.Length())
      return;

   TString libName;
   const char *colon = strchr(cFilename, ':');
   if (colon)
      // old version, where source file name is prepended by "TAG:"
      libName = TString(cFilename, colon - cFilename);
   else
      // New version, check class's libname.
      // If libname is dir/libMyLib.so, check Root.Html.MyLib
      // If libname is myOtherLib.so.2.3, check Root.Html.myOtherLib
      // (i.e. remove directories, "lib" prefix, and any "extension")
      if (classPtr->GetSharedLibs()) {
         // first one is the class's lib
         TString libname(classPtr->GetSharedLibs());
         Ssiz_t posSpace = libname.First(' ');
         if (posSpace != kNPOS)
            libname.Remove(posSpace, libname.Length());
         TString libnameBase = gSystem->BaseName(libname);
         if (libnameBase.BeginsWith("lib"))
            libnameBase.Remove(0, 3);
         Ssiz_t posExt = libnameBase.First('.');
         if (posExt != '.')
            libnameBase.Remove(posExt, libnameBase.Length());
         if (libnameBase.Length())
            libName = libnameBase;
      }

   filename = cFilename;
   TString htmlFileName;
   if (!filename.Length() ||
      !gSystem->FindFile(fPathInfo.fInputPath, filename, kReadPermission)) {
      htmlFileName = GetURL(libName);
   } else
      htmlFileName = "./";

   if (htmlFileName.Length()) {
      filename = htmlFileName;
      TString className(classPtr->GetName());
      TDocOutput output(*const_cast<THtml*>(this));
      output.NameSpace2FileName(className);
      gSystem->PrependPathName(filename, className);
      filename = className;
      filename.ReplaceAll("\\", "/");
      filename += ".html";
   } else filename.Remove(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the html file name for a class named classname.
/// Returns 0 if the class is not documented.

const char* THtml::GetHtmlFileName(const char* classname) const
{
   TClassDocInfo* cdi = (TClassDocInfo*) fDocEntityInfo.fClasses.FindObject(classname);
   if (cdi)
      return cdi->GetHtmlFileName();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to class with name.

TClass *THtml::GetClass(const char *name1) const
{
   if(!name1 || !name1[0]) return 0;
   // no doc for internal classes
   if (strstr(name1,"ROOT::")==name1) {
      Bool_t ret = kTRUE;
      if (!strncmp(name1 + 6,"Math", 4))   ret = kFALSE;
      if (ret) return 0;
   }

   TClassDocInfo* cdi = (TClassDocInfo*)fDocEntityInfo.fClasses.FindObject(name1);
   if (!cdi) return 0;
   TClass *cl = dynamic_cast<TClass*>(cdi->GetClass());
   // hack to get rid of prec_stl types
   // TClassEdit checks are far too slow...
   /*
   if (cl && GetDeclFileName(cl) &&
       (strstr(GetDeclFileName(cl),"prec_stl/") || !strstr(classPtr->GetDeclFileName(), "include/c++/") )
      cl = 0;
   */
   TString declFileName;
   if (cl && GetDeclFileName(cl, kFALSE, declFileName))
      return cl;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return declaration file name; return the full path if filesys is true.

bool THtml::GetDeclFileName(TClass * cl, Bool_t filesys, TString& out_name) const
{
   return GetDeclImplFileName(cl, filesys, true, out_name);
}

////////////////////////////////////////////////////////////////////////////////
/// Return implementation file name

bool THtml::GetImplFileName(TClass * cl, Bool_t filesys, TString& out_name) const
{
   return GetDeclImplFileName(cl, filesys, false, out_name);
}

////////////////////////////////////////////////////////////////////////////////
/// Combined implementation for GetDeclFileName(), GetImplFileName():
/// Return declaration / implementation file name (depending on decl);
/// return the full path if filesys is true.

bool THtml::GetDeclImplFileName(TClass * cl, bool filesys, bool decl, TString& out_name) const
{
   out_name = "";

   R__LOCKGUARD(GetMakeClassMutex());
   TClassDocInfo* cdi = (TClassDocInfo*) fDocEntityInfo.fClasses.FindObject(cl->GetName());
   // whether we need to determine the fil name
   bool determine = (!cdi); // no cdi
   if (!determine) determine |=  decl &&  filesys && !cdi->GetDeclFileSysName()[0];
   if (!determine) determine |=  decl && !filesys && !cdi->GetDeclFileName()[0];
   if (!determine) determine |= !decl &&  filesys && !cdi->GetImplFileSysName()[0];
   if (!determine) determine |= !decl && !filesys && !cdi->GetImplFileName()[0];
   if (determine) {
      TString name;
      TString sysname;
      if (decl) {
         if (!GetFileDefinition().GetDeclFileName(cl, name, sysname))
            return false;
      } else {
         if (!GetFileDefinition().GetImplFileName(cl, name, sysname))
            return false;
      }
      if (cdi) {
         if (decl) {
            if (!cdi->GetDeclFileName() || !cdi->GetDeclFileName()[0])
               cdi->SetDeclFileName(name);
            if (!cdi->GetDeclFileSysName() || !cdi->GetDeclFileSysName()[0])
               cdi->SetDeclFileSysName(sysname);
         } else {
            if (!cdi->GetImplFileName() || !cdi->GetImplFileName()[0])
               cdi->SetImplFileName(name);
            if (!cdi->GetImplFileSysName() || !cdi->GetImplFileSysName()[0])
               cdi->SetImplFileSysName(sysname);
         }
      }

      if (filesys) out_name = sysname;
      else         out_name = name;
      return true;
   }
   if (filesys) {
      if (decl) out_name = cdi->GetDeclFileSysName();
      else      out_name = cdi->GetImplFileSysName();
   } else {
      if (decl) out_name = cdi->GetDeclFileName();
      else      out_name = cdi->GetImplFileName();
   }
   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the output directory as set by SetOutputDir().
/// Create it if it doesn't exist and if createDir is kTRUE.

const TString& THtml::GetOutputDir(Bool_t createDir /*= kTRUE*/) const
{
   if (createDir) {
      R__LOCKGUARD(GetMakeClassMutex());

      gSystem->ExpandPathName(const_cast<THtml*>(this)->fPathInfo.fOutputDir);
      Long64_t sSize;
      Long_t sId, sFlags, sModtime;
      if (fPathInfo.fOutputDir.EndsWith("/") || fPathInfo.fOutputDir.EndsWith("\\"))
         fPathInfo.fOutputDir.Remove(fPathInfo.fOutputDir.Length() - 1);
      Int_t st = gSystem->GetPathInfo(fPathInfo.fOutputDir, &sId, &sSize, &sFlags, &sModtime);
      if (st || !(sFlags & 2)) {
         if (st == 0)
            Error("GetOutputDir", "output directory %s is an existing file",
                  fPathInfo.fOutputDir.Data());
         else if (gSystem->MakeDirectory(fPathInfo.fOutputDir) == -1)
            Error("GetOutputDir", "output directory %s does not exist and can't create it", fPathInfo.fOutputDir.Data());
      }
   }
   return fPathInfo.fOutputDir;
}

////////////////////////////////////////////////////////////////////////////////
/// Check whether cl is a namespace

Bool_t THtml::IsNamespace(const TClass*cl)
{
   return (cl->Property() & kIsNamespace);
}

////////////////////////////////////////////////////////////////////////////////
/// Load all libraries known to ROOT via the rootmap system.

void THtml::LoadAllLibs()
{
   gSystem->LoadAllLibraries();
}


////////////////////////////////////////////////////////////////////////////////
/// Produce documentation for all the classes specified in the filter (by default "*")
/// To process all classes having a name starting with XX, do:
///        html.MakeAll(kFALSE,"XX*");
/// If force=kFALSE (default), only the classes that have been modified since
/// the previous call to this function will be generated.
/// If force=kTRUE, all classes passing the filter will be processed.
/// If numthreads is != -1, use numthreads threads, else decide automatically
/// based on the number of CPUs.

void THtml::MakeAll(Bool_t force, const char *filter, int numthreads /*= -1*/)
{
   MakeIndex(filter);

   if (numthreads == 1) {
      // CreateListOfClasses(filter); already done by MakeIndex
      TClassDocInfo* classinfo = 0;
      TIter iClassInfo(&fDocEntityInfo.fClasses);
      UInt_t count = 0;

      while ((classinfo = (TClassDocInfo*)iClassInfo())) {
         if (!classinfo->IsSelected())
            continue;
         fCounter.Form("%5d", fDocEntityInfo.fClasses.GetSize() - count++);
         MakeClass(classinfo, force);
      }
   } else {
      if (numthreads == -1) {
         SysInfo_t sysinfo;
         gSystem->GetSysInfo(&sysinfo);
         numthreads = sysinfo.fCpus;
         if (numthreads < 1)
            numthreads = 2;
      }
      fThreadedClassCount = 0;
      fThreadedClassIter = new TIter(&fDocEntityInfo.fClasses);
      THtmlThreadInfo hti(this, force);
      if (!fMakeClassMutex && gGlobalMutex) {
         gGlobalMutex->Lock();
         fMakeClassMutex = gGlobalMutex->Factory(kTRUE);
         gGlobalMutex->UnLock();
      }

      TList threads;
      gSystem->Load("libThread");
      while (--numthreads >= 0) {
         TThread* thread = new TThread(MakeClassThreaded, &hti);
         thread->Run();
         threads.Add(thread);
      }

      TIter iThread(&threads);
      TThread* thread = 0;
      Bool_t wait = kTRUE;
      while (wait) {
         while (wait && (thread = (TThread*) iThread()))
            wait &= (thread->GetState() == TThread::kRunningState);
         gSystem->ProcessEvents();
         gSystem->Sleep(500);
      }

      iThread.Reset();
      while ((thread = (TThread*) iThread()))
         thread->Join();
   }
   fCounter.Remove(0);
}


////////////////////////////////////////////////////////////////////////////////
/// Make HTML files for a single class
///
///
/// Input: className - name of the class to process
///

void THtml::MakeClass(const char *className, Bool_t force)
{
   CreateListOfClasses("*");

   TClassDocInfo* cdi = (TClassDocInfo*)fDocEntityInfo.fClasses.FindObject(className);
   if (!cdi) {
      if (!TClassEdit::IsStdClass(className)) // stl classes won't be available, so no warning
         Error("MakeClass", "Unknown class '%s'!", className);
      return;
   }

   MakeClass(cdi, force);
}

////////////////////////////////////////////////////////////////////////////////
/// Make HTML files for a single class
///
///
/// Input: cdi - doc info for class to process
///

void THtml::MakeClass(void *cdi_void, Bool_t force)
{
   if (!fDocEntityInfo.fClasses.GetSize())
      CreateListOfClasses("*");

   TClassDocInfo* cdi = (TClassDocInfo*) cdi_void;
   TClass* currentClass = dynamic_cast<TClass*>(cdi->GetClass());

   if (!currentClass) {
      if (!cdi->GetClass() &&
          !TClassEdit::IsStdClass(cdi->GetName())) // stl classes won't be available, so no warning
         Error("MakeClass", "Class '%s' is known, but I cannot find its TClass object!", cdi->GetName());
      return;
   }
   TString htmlFile(cdi->GetHtmlFileName());
   if (htmlFile.Length()
       && (htmlFile.BeginsWith("http://")
           || htmlFile.BeginsWith("https://")
           || gSystem->IsAbsoluteFileName(htmlFile))
       ) {
      htmlFile.Remove(0);
   }
   if (htmlFile.Length()) {
      TClassDocOutput cdo(*this, currentClass, &cdi->GetListOfTypedefs());
      cdo.Class2Html(force);
      cdo.MakeTree(force);
   } else {
      TString what(cdi->GetName());
      what += " (sources not found)";
      Printf(fCounterFormat.Data(), "-skipped-", fCounter.Data(), what.Data());
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Entry point of worker threads for multi-threaded MakeAll().
/// info points to an (internal) THtmlThreadInfo object containing the current
/// THtml object, and whether "force" was passed to MakeAll().
/// The thread will poll GetNextClass() until no further class is available.

void* THtml::MakeClassThreaded(void* info) {
   const THtmlThreadInfo* hti = (const THtmlThreadInfo*)info;
   if (!hti) return 0;
   TClassDocInfo* classinfo = 0;
   while ((classinfo = hti->GetHtml()->GetNextClass()))
      hti->GetHtml()->MakeClass(classinfo, hti->GetForce());

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create the index files for the product, modules, all types, etc.
/// By default all classes are indexed (if filter="*");
/// to generate an index for all classes starting with "XX", do
///    html.MakeIndex("XX*");

void THtml::MakeIndex(const char *filter)
{
   CreateListOfClasses(filter);

   TDocOutput output(*this);
   // create indices
   output.CreateTypeIndex();
   output.CreateClassTypeDefs();
   output.CreateModuleIndex();
   output.CreateClassIndex();
   output.CreateProductIndex();

   // create a class hierarchy
   output.CreateHierarchy();
}


////////////////////////////////////////////////////////////////////////////////
/// Make an inheritance tree
///
///
/// Input: className - name of the class to process
///

void THtml::MakeTree(const char *className, Bool_t force)
{
   // create canvas & set fill color
   TClass *classPtr = GetClass(className);

   if (!classPtr) {
      Error("MakeTree", "Unknown class '%s' !", className);
      return;
   }

   TClassDocOutput cdo(*this, classPtr, 0);
   cdo.MakeTree(force);
}

////////////////////////////////////////////////////////////////////////////////
/// Set whether "dot" (a GraphViz utility) is available

void THtml::SetFoundDot(Bool_t found) {
   if (found) fPathInfo.fFoundDot = PathInfo_t::kDotFound;
   else fPathInfo.fFoundDot = PathInfo_t::kDotNotFound;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the files available in the file system below fPathInfo.fInputPath

void THtml::SetLocalFiles() const
{
   if (fLocalFiles) delete fLocalFiles;
   fLocalFiles = new TFileSysDB(fPathInfo.fInputPath, fPathInfo.fIgnorePath + "|(\\b" + GetOutputDir(kFALSE) + "\\b)" , 6);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the module defining object to be used; can also be a user derived
/// object (a la traits).

void THtml::SetModuleDefinition(const TModuleDefinition& md)
{
   delete fModuleDef;
   fModuleDef = (TModuleDefinition*) md.Clone();
   fModuleDef->SetOwner(const_cast<THtml*>(this));
}


////////////////////////////////////////////////////////////////////////////////
/// Set the file defining object to be used; can also be a user derived
/// object (a la traits).

void THtml::SetFileDefinition(const TFileDefinition& md)
{
   delete fFileDef;
   fFileDef = (TFileDefinition*) md.Clone();
   fFileDef->SetOwner(const_cast<THtml*>(this));
}


////////////////////////////////////////////////////////////////////////////////
/// Set the path defining object to be used; can also be a user derived
/// object (a la traits).

void THtml::SetPathDefinition(const TPathDefinition& md)
{
   delete fPathDef;
   fPathDef = (TPathDefinition*) md.Clone();
   fPathDef->SetOwner(const_cast<THtml*>(this));
}


////////////////////////////////////////////////////////////////////////////////
/// Set the directory containing the source files.
/// The source file for a class MyClass will be searched
/// by prepending dir to the value of
/// MyClass::Class()->GetImplFileName() - which can contain
/// directory information!
/// Also resets the class structure, in case new files can
/// be found after this call.

void THtml::SetInputDir(const char *dir)
{
   fPathInfo.fInputPath = dir;
   gSystem->ExpandPathName(fPathInfo.fInputPath);

   // reset class table
   fDocEntityInfo.fClasses.Clear();
   fDocEntityInfo.fModules.Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the directory where the HTML pages shuold be written to.
/// If the directory does not exist it will be created when needed.

void THtml::SetOutputDir(const char *dir)
{
   fPathInfo.fOutputDir = dir;
#ifdef R__WIN32
   fPathInfo.fOutputDir.ReplaceAll("/","\\");
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Explicitly set a decl file name for TClass cl.

void THtml::SetDeclFileName(TClass* cl, const char* filename)
{
   TClassDocInfo* cdi = (TClassDocInfo*) fDocEntityInfo.fClasses.FindObject(cl->GetName());
   if (!cdi) {
      cdi = new TClassDocInfo(cl, "" /*html*/, "" /*fsdecl*/, "" /*fsimpl*/, filename);
      fDocEntityInfo.fClasses.Add(cdi);
   } else
      cdi->SetDeclFileName(filename);
}

////////////////////////////////////////////////////////////////////////////////
/// Explicitly set a impl file name for TClass cl.

void THtml::SetImplFileName(TClass* cl, const char* filename)
{
   TClassDocInfo* cdi = (TClassDocInfo*) fDocEntityInfo.fClasses.FindObject(cl->GetName());
   if (!cdi) {
      cdi = new TClassDocInfo(cl, "" /*html*/, "" /*fsdecl*/, "" /*fsimpl*/, 0 /*decl*/, filename);
      fDocEntityInfo.fClasses.Add(cdi);
   } else
      cdi->SetImplFileName(filename);
}

////////////////////////////////////////////////////////////////////////////////
/// Get short type name, i.e. with default templates removed.

const char* THtml::ShortType(const char* name) const
{
   const char* tmplt = strchr(name, '<');
   if (!tmplt) return name;
   tmplt = strrchr(tmplt, ':');
   if (tmplt > name && tmplt[-1] == ':') {
      // work-around for CINT bug: template instantiation can produce bogus
      // typedefs e.g. in namespace ROOT::Math::ROOT::Math instead of ROOT::Math.
      TString namesp(name, tmplt - name - 1);
      // is the enclosing namespace known?
      if (!GetClass(namesp)) return name;
   }
   TObject* scn = fDocEntityInfo.fShortClassNames.FindObject(name);
   if (!scn) {
      scn = new TNamed(name, TClassEdit::ShortType(name, 1<<7));
      fDocEntityInfo.fShortClassNames.Add(scn);
   }
   return scn->GetTitle();
}
