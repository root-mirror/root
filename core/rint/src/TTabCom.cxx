// @(#)root/rint:$Id$
// Author: Christian Lacunza <lacunza@cdfsg6.lbl.gov>   27/04/99

// Modified by Artur Szostak <artur@alice.phy.uct.ac.za> : 1 June 2003
//   Added support for namespaces.

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TTabCom                                                                //
//                                                                        //
// This class performs basic tab completion.                              //
// You should be able to hit [TAB] to complete a partially typed:         //
//                                                                        //
//   username                                                             //
//   environment variable                                                 //
//   preprocessor directive                                               //
//   pragma                                                               //
//   filename (with a context-sensitive path)                             //
//   public member function or data member (including base classes)       //
//   global variable, function, or class name                             //
//                                                                        //
// Also, something like                                                   //
//                                                                        //
//   someObject->Func([TAB]                                               //
//   someObject.Func([TAB]                                                //
//   someClass::Func([TAB]                                                //
//   someClass var([TAB]                                                  //
//   new someClass([TAB]                                                  //
//                                                                        //
// will print a list of prototypes for the indicated                      //
// method or constructor.                                                 //
//                                                                        //
// Current limitations and bugs:                                          //
//                                                                        //
//  1. you can only use one member access operator at a time.             //
//     eg, this will work: gROOT->GetListOfG[TAB]                         //
//     but this will not:  gROOT->GetListOfGlobals()->Conta[TAB]          //
//                                                                        //
//  2. nothing is guaranteed to work on windows                           //
//     (for one thing, /bin/env and /etc/passwd are hardcoded)            //
//                                                                        //
//  3. CINT shortcut #2 is deliberately not supported.                    //
//     (using "operator.()" instead of "operator->()")                    //
//                                                                        //
//  4. most identifiers (including C++ identifiers, usernames,            //
//     environment variables, etc)                                        //
//     are restriceted to this character set: [_a-zA-Z0-9]                //
//     therefore, you won't be able to complete things like               //
//                                                                        //
//          operator new                                                  //
//          operator+                                                     //
//          etc                                                           //
//                                                                        //
//  5. ~whatever[TAB] always tries to complete a username.                //
//     use whitespace (~ whatever[TAB]) if you want to complete a global  //
//     identifier.                                                        //
//                                                                        //
//  6. CINT shortcut #3 is not supported when trying to complete          //
//     the name of a global object.  (it is supported when trying to      //
//     complete a member of a global object)                              //
//                                                                        //
//  7. the list of #pragma's is hardcoded                                 //
//     (ie not obtained from the interpreter at runtime)                  //
//     ==> user-defined #pragma's will not be recognized                  //
//                                                                        //
//  8. the system include directories are also hardcoded                  //
//     because i don't know how to get them from the interpreter.         //
//     fons, maybe they should be #ifdef'd for the different sytems?      //
//                                                                        //
//  9. the TabCom.FileIgnore resource is always applied, even if you      //
//     are not trying to complete a filename.                             //
//                                                                        //
// 10. anything in quotes is assumed to be a filename                     //
//     so (among other things) you can't complete a quoted class name:    //
//     eg, TClass class1( "TDict[TAB]                                     //
//     this won't work... looks for a file in pwd starting with TDict     //
//                                                                        //
// 11. the prototypes tend to omit the word "const" a lot.                //
//     this is a problem with ROOT or CINT.                               //
//                                                                        //
// 12. when listing ambiguous matches, only one column is used,           //
//     even if there are many completions.                                //
//                                                                        //
// 13. anonymous objects are not currently identified                     //
//     so, for example,                                                   //
//                                                                        //
//          root> printf( TString([TAB                                    //
//                                                                        //
//     gives an error message instead of listing TString's constructors.  //
//     (this could be fixed)                                              //
//                                                                        //
// 14. the routine that adds the "appendage" isn't smart enough to know   //
//     if it's already there:                                             //
//                                                                        //
//          root> TCanvas::Update()                                       //
//              press [TAB] here ^                                        //
//          root> TCanvas::Update()()                                     //
//     (this could be fixed)                                              //
//                                                                        //
// 15. the appendage is only applied if there is exactly 1 match.         //
//     eg, this                                                           //
//                                                                        //
//          root> G__at[TAB]                                              //
//          root> G__ateval                                               //
//                                                                        //
//     happens instead of this                                            //
//                                                                        //
//          root> G__at[TAB]                                              //
//          root> G__ateval(                                              //
//                                                                        //
//     because there are several overloaded versions of G__ateval().      //
//     (this could be fixed)                                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cassert>
#include <set>

#include "RConfigure.h"
#include "TTabCom.h"
#include "TClass.h"
#include "TClassTable.h"
#include "TDataMember.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMethod.h"
#include "TEnv.h"
#include "TBenchmark.h"
#include "TError.h"
#include "TGlobal.h"
#include "TList.h"
#include "THashList.h"
#include "TObjString.h"
#include "Getline.h"
#include "TFunction.h"
#include "TMethodArg.h"
#include "TInterpreter.h"
#include "Riostream.h"
#include "Rstrstream.h"
#include "strlcpy.h"

#define BUF_SIZE 1024 // must be smaller than/equal to fgLineBufSize in Getline.cxx and
                      // lineBufSize in cppcompleter.py
#define IfDebug(x)  if(gDebug==TTabCom::kDebug) x

#ifdef R__WIN32
const char kDelim = ';';
#else
const char kDelim = ':';
#endif


ClassImp(TTabCom);
// ----------------------------------------------------------------------------
//
//             global/file scope variables
//
TTabCom *gTabCom = 0;

// ----------------------------------------------------------------------------
//
//              constructors
//

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TTabCom::TTabCom():
   fpClasses(0),
   fPrevInterpMarker(0),
   fpDirectives(0),
   fpEnvVars(0),
   fpFiles(0),
   fpGlobals(0),
   fpPragmas(0),
   fpSysIncFiles(0),
   fpUsers(0),
   fBuf(0),
   fpLoc(0),
   fVarIsPointer(kFALSE),
   fLastIter(0)
{
   InitPatterns();
}

//
//              constructors
//
// ----------------------------------------------------------------------------

TTabCom::~TTabCom()
{
   // Destructor.

   ClearAll();
   ClearSysIncFiles(); // this one stays cached
   ClearUsers();       // this one stays cached
}

// ----------------------------------------------------------------------------
//
//              public member functions
//


////////////////////////////////////////////////////////////////////////////////
/// Clear classes and namespace collections.

void TTabCom::ClearClasses()
{
   if (fpClasses) {
      delete fpClasses;
      fpClasses = 0;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Forget all Cpp directives seen so far.

void TTabCom::ClearCppDirectives()
{
   if (!fpDirectives)
      return;
   fpDirectives->Delete(0);
   delete fpDirectives;
   fpDirectives = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Forget all environment variables seen so far.

void TTabCom::ClearEnvVars()
{
   if (!fpEnvVars)
      return;
   fpEnvVars->Delete(0);
   delete fpEnvVars;
   fpEnvVars = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Close all files.

void TTabCom::ClearFiles()
{
   if (!fpFiles)
      return;
   fpFiles->Delete(0);
   delete fpFiles;
   fpFiles = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Forget all global functions seen so far.
/// Not needed anymore. Use gROOT->GetListOfGlobalFunctions()

void TTabCom::ClearGlobalFunctions()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Forget all global variables seen so far.
/// With teh new implamentation the list is gROOT->GetListOfGlobals(true).

void TTabCom::ClearGlobals()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Forget all pragmas seen so far.

void TTabCom::ClearPragmas()
{
   if (!fpPragmas)
      return;
   fpPragmas->Delete(0);
   delete fpPragmas;
   fpPragmas = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Close system files.

void TTabCom::ClearSysIncFiles()
{
   if (!fpSysIncFiles)
      return;
   fpSysIncFiles->Delete(0);
   delete fpSysIncFiles;
   fpSysIncFiles = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Forget all user seen so far.

void TTabCom::ClearUsers()
{
   if (!fpUsers)
      return;
   fpUsers->Delete(0);
   delete fpUsers;
   fpUsers = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// clears all lists
/// except for user names and system include files.

void TTabCom::ClearAll()
{
   ClearClasses();
   ClearCppDirectives();
   ClearEnvVars();
   ClearFiles();
   ClearGlobalFunctions();
   ClearGlobals();
   ClearPragmas();
//   ClearSysIncFiles(); <-- this one stays cached
//   ClearUsers();       <-- this one stays cached
}

////////////////////////////////////////////////////////////////////////////////
/// Do the class rehash.

void TTabCom::RehashClasses()
{
   ClearClasses();
   GetListOfClasses();
}

////////////////////////////////////////////////////////////////////////////////
/// Cpp rehashing.

void TTabCom::RehashCppDirectives()
{
   ClearCppDirectives();
   GetListOfCppDirectives();
}

////////////////////////////////////////////////////////////////////////////////
/// Environemnt variables rehashing.

void TTabCom::RehashEnvVars()
{
   ClearEnvVars();
   GetListOfEnvVars();
}

////////////////////////////////////////////////////////////////////////////////
/// Close files.

void TTabCom::RehashFiles()
{
   ClearFiles();                /* path unknown */
}                               // think about this

////////////////////////////////////////////////////////////////////////////////
/// Reload global functions.

void TTabCom::RehashGlobalFunctions()
{
   GetListOfGlobalFunctions();
}

////////////////////////////////////////////////////////////////////////////////
/// Reload globals.

void TTabCom::RehashGlobals()
{
   ClearGlobals();
   GetListOfGlobals();
}

////////////////////////////////////////////////////////////////////////////////
/// Reload pragmas.

void TTabCom::RehashPragmas()
{
   ClearPragmas();
   GetListOfPragmas();
}

////////////////////////////////////////////////////////////////////////////////
/// Reload system include files.

void TTabCom::RehashSysIncFiles()
{
   ClearSysIncFiles();
   GetListOfSysIncFiles();
}

////////////////////////////////////////////////////////////////////////////////
/// Reload users.

void TTabCom::RehashUsers()
{
   ClearUsers();
   GetListOfUsers();
}

////////////////////////////////////////////////////////////////////////////////
/// clears and then rebuilds all lists
/// except for user names and system include files.

void TTabCom::RehashAll()
{
   RehashClasses();
   RehashCppDirectives();
   RehashEnvVars();
   RehashFiles();
   RehashGlobalFunctions();
   RehashGlobals();
   RehashPragmas();
//   RehashSysIncFiles(); <-- this one stays cached
//   RehashUsers();       <-- this one stays cached
}

////////////////////////////////////////////////////////////////////////////////
/// Return the list of classes.

const TSeqCollection *TTabCom::GetListOfClasses()
{
   if (!fpClasses) {
      fpClasses = new THashList;
      // Iterate over the table from the map file.
      THashList* entries = gInterpreter->GetMapfile()->GetTable();
      TIter next(entries);
      while (const auto key = next()) {
         // This is not needed with the new rootmap format
         const char* className = key->GetName();
         if (!strncmp(className, "Library.", 8))
            className += 8;

         if (!strstr(className, ".h"))
            fpClasses->Add(new TObjString(className));
      }

      // We might have autoload entries that don't have a rootmap entry
      // (libCore) and no interpreter info (not yet loaded).
      TClassTable::Init(); // reset counter
      while (const char* className = TClassTable::Next()) {
         if (!fpClasses->FindObject(className)) {
            fpClasses->Add(new TObjString(className));
         }
      }

      // Add possible identifiers coming from the global module index.
      // FIXME: Consolidate all this -- we have ultimate source of information
      // clang's identifier table and the ASTIndentifiers.
      gInterpreter->AddAvailableIndentifiers(*fpClasses);
   }

   if (fPrevInterpMarker != gInterpreter->GetInterpreterStateMarker()) {
      ClassInfo_t* ci = gInterpreter->ClassInfo_Factory(kFALSE /*all*/);
      while (gInterpreter->ClassInfo_Next(ci)) {
         const char* className = gInterpreter->ClassInfo_FullName(ci);
         if (strstr(className, "(anonymous)"))
            continue;
         if (!fpClasses->FindObject(className)) {
            fpClasses->Add(new TObjString(className));
         }
      }
      gInterpreter->ClassInfo_Delete(ci);
   }

   return fpClasses;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the list of CPP directives.

const TSeqCollection *TTabCom::GetListOfCppDirectives()
{
   if (!fpDirectives) {
      fpDirectives = new TContainer;

      fpDirectives->Add(new TObjString("if"));
      fpDirectives->Add(new TObjString("ifdef"));
      fpDirectives->Add(new TObjString("ifndef"));
      fpDirectives->Add(new TObjString("elif"));
      fpDirectives->Add(new TObjString("else"));
      fpDirectives->Add(new TObjString("endif"));
      fpDirectives->Add(new TObjString("include"));
      fpDirectives->Add(new TObjString("define"));
      fpDirectives->Add(new TObjString("undef"));
      fpDirectives->Add(new TObjString("line"));
      fpDirectives->Add(new TObjString("error"));
      fpDirectives->Add(new TObjString("pragma"));
   }

   return fpDirectives;
}

////////////////////////////////////////////////////////////////////////////////
/// "path" should be initialized with a colon separated list of
/// system directories

const TSeqCollection *TTabCom::GetListOfFilesInPath(const char path[])
{
   static TString previousPath;

   if (path && fpFiles && strcmp(path, previousPath) == 0) {
      return fpFiles;
   } else {
      ClearFiles();

      fpFiles = NewListOfFilesInPath(path);
      previousPath = path;
   }

   return fpFiles;
}

////////////////////////////////////////////////////////////////////////////////
/// Uses "env" (Unix) or "set" (Windows) to get list of environment variables.

const TSeqCollection *TTabCom::GetListOfEnvVars()
{
   if (!fpEnvVars) {
      TString outf = ".TTabCom-";
      FILE *fout = gSystem->TempFileName(outf);
      if (!fout) return 0;
      fclose(fout);
      TString cmd;

#ifndef WIN32
      char *env = gSystem->Which(gSystem->Getenv("PATH"), "env", kExecutePermission);
      if (!env)
         return 0;
      cmd = env;
      cmd += " > ";
      delete [] env;
#else
      cmd = "set > ";
#endif
      cmd += outf;
      cmd += "\n";
      gSystem->Exec(cmd.Data());

      // open the file
      std::ifstream file1(outf);
      if (!file1) {
         Error("TTabCom::GetListOfEnvVars", "could not open file \"%s\"",
               outf.Data());
         gSystem->Unlink(outf);
         return 0;
      }
      // parse, add
      fpEnvVars = new TContainer;
      TString line;
      while (file1)             // i think this loop goes one time extra which
         // results in an empty string in the list, but i don't think it causes any
         // problems.
      {
         line.ReadToDelim(file1, '=');
         file1.ignore(32000, '\n');
         fpEnvVars->Add(new TObjString(line.Data()));
      }

      file1.close();
      gSystem->Unlink(outf);
   }

   return fpEnvVars;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the list of globals.

const TSeqCollection *TTabCom::GetListOfGlobals()
{
   return (TSeqCollection*)gROOT->GetListOfGlobals(true);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the list of global functions.

TCollection *TTabCom::GetListOfGlobalFunctions()
{
   return gROOT->GetListOfGlobalFunctions(true);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the list of pragmas

const TSeqCollection *TTabCom::GetListOfPragmas()
{
   if (!fpPragmas) {
      fpPragmas = new TContainer;

      fpPragmas->Add(new TObjString("ANSI "));
      fpPragmas->Add(new TObjString("autocompile "));
      fpPragmas->Add(new TObjString("bytecode "));
      fpPragmas->Add(new TObjString("compile "));
      fpPragmas->Add(new TObjString("endbytecode "));
      fpPragmas->Add(new TObjString("endcompile "));
      fpPragmas->Add(new TObjString("include "));
      fpPragmas->Add(new TObjString("includepath "));
      fpPragmas->Add(new TObjString("K&R "));
      fpPragmas->Add(new TObjString("link "));
      fpPragmas->Add(new TObjString("preprocess "));
      fpPragmas->Add(new TObjString("preprocessor "));
      fpPragmas->Add(new TObjString("security level"));
      // "setertti "  omitted. Ordinary user should not use this statement
      // "setstdio "  omitted. Ordinary user should not use this statement
      // "setstream " omitted. Ordinary user should not use this statement
      // "stub"       omitted. Ordinary user should not use this statement

   }

   return fpPragmas;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the list of system include files.

const TSeqCollection *TTabCom::GetListOfSysIncFiles()
{
   if (!fpSysIncFiles) {
      fpSysIncFiles = NewListOfFilesInPath(GetSysIncludePath());
   }

   return fpSysIncFiles;
}

////////////////////////////////////////////////////////////////////////////////
/// reads from "/etc/passwd"

const TSeqCollection *TTabCom::GetListOfUsers()
{
   if (!fpUsers) {
      fpUsers = new TContainer;

      std::ifstream passwd;
      TString user;

      passwd.open("/etc/passwd");
      while (passwd) {
         user.ReadToDelim(passwd, ':');
         fpUsers->Add(new TObjString(user));
         passwd.ignore(32000, '\n');
      }
      passwd.close();
   }

   return fpUsers;
}

//
//              public member functions
//
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
//
//                           static utility functions
//

////////////////////////////////////////////////////////////////////////////////
///[static utility function]///////////////////////////////////////////
///
///  if all the strings in "*pList" have the same ith character,
///  that character is returned.
///  otherwise 0 is returned.
///
///  any string "s" for which "ExcludedByFignore(s)" is true
///  will be ignored unless All the strings in "*pList"
///  are "ExcludedByFignore()"
///
///  in addition, the number of strings which were not
///  "ExcludedByFignore()" is returned in "nGoodStrings".
///
//////////////////////////////////////////////////////////////////////////

Char_t TTabCom::AllAgreeOnChar(int i, const TSeqCollection * pList,
                               Int_t & nGoodStrings)
{
   assert(pList != 0);

   TIter next(pList);
   TObject *pObj;
   const char *s = "";
   char ch0;
   Bool_t isGood;
   Bool_t atLeast1GoodString;

   // init
   nGoodStrings = 0;
   atLeast1GoodString = kFALSE;

   // first look for a good string
   do {
      if ((pObj = next())) {
         s = pObj->GetName();
         isGood = !ExcludedByFignore(s);
         if (isGood) {
            atLeast1GoodString = kTRUE;
            nGoodStrings += 1;
         }
      } else {
         // reached end of list without finding a single good string.
         // just use the first one.
         next.Reset();
         pObj = next();
         if (pObj) s = pObj->GetName();
         break;
      }
   }
   while (!isGood);

   // found a good string...
   ch0 = s[i];

   // all subsequent good strings must have the same ith char
   do {
      if ((pObj = next())) {
         s = pObj->GetName();
         isGood = !ExcludedByFignore(s);
         if (isGood)
            nGoodStrings += 1;
      } else
         return ch0;
   }
   while (((int) strlen(s) >= i && s[i] == ch0) ||
          (atLeast1GoodString && !isGood));

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///[static utility function]/////////////////////////////
///
///  adds a TObjString to "*pList"
///  for each entry found in the system directory "dirName"
///
///  directories that do not exist are silently ignored.
///
///////////////////////////////////////////////////////////

void TTabCom::AppendListOfFilesInDirectory(const char dirName[],
                                           TSeqCollection * pList)
{
   assert(dirName != 0);
   assert(pList != 0);

   // open the directory
   void *dir = gSystem->OpenDirectory(dirName);

   // it is normal for this function to receive names of directories that do not exist.
   // they should be ignored and should not generate any error messages.
   if (!dir)
      return;

   // put each filename in the list
   const char *tmp_ptr;         // gSystem->GetDirEntry() returns 0 when no more files.
   TString fileName;

   while ((tmp_ptr = gSystem->GetDirEntry(dir))) {
      fileName = tmp_ptr;

      // skip "." and ".."
      if (fileName == "." || fileName == "..")
         continue;

      // add to list
      pList->Add(new TObjString(dirName + fileName.Prepend("/")));
   }
   // NOTE:
   // with a path like "/usr/include:/usr/include/CC:$ROOTDIR/include:$ROOTDIR/cint/include:..."
   // the above loop could get traversed 700 times or more.
   // ==> keep it minimal or it could cost whole seconds on slower machines.
   // also: TClonesArray doesn't help.

   // close the directory
   gSystem->FreeDirectory(dir);
}

// -----\/-------- homemade RTTI ---------------\/------------------------
////////////////////////////////////////////////////////////////////////////////
///[static utility function]/////////////////////////////
///
///  returns empty string on failure.
///  otherwise returns something like this: "TROOT*".
///  fails for non-class types (ie, int, char, etc).
///  fails for pointers to functions.
///
////////////////////////////////////

TString TTabCom::DetermineClass(const char varName[])
{

   ///////////////////////////////////
   //
   //  note that because of the strange way this function works,
   //  CINT will print
   //
   //     Error: No symbol asdf in current scope  FILE:/var/tmp/gaaa001HR LINE:1
   //
   //  if "varName" is not defined. (in this case, varName=="asdf")
   //  i don't know how to suppress this.
   //
   ///////////////////////////////////

   assert(varName != 0);
   IfDebug(std::cerr << "DetermineClass(\"" << varName << "\");" << std::endl);

   TString outf = ".TTabCom-";
   FILE *fout = gSystem->TempFileName(outf);
   if (!fout) return "";
   fclose(fout);

   TString cmd(".> ");
   cmd += outf;
   //redirect
   gROOT->ProcessLineSync(cmd.Data());

   cmd = "gROOT->ProcessLine(\"";
   cmd += varName;
   cmd += "\");";
   cmd += "\n";
   gROOT->ProcessLineSync(cmd.Data());
   // the type of the variable whose name is "varName"
   // should now be stored on disk in the file "tmpfile"

   gROOT->ProcessLineSync(".>");

   TString type = "";
   int c;

   // open the file
   std::ifstream file1(outf);
   if (!file1) {
      Error("TTabCom::DetermineClass", "could not open file \"%s\"",
            outf.Data());
      goto cleanup;
   }
   // first char should be '(', which we can ignore.
   c = file1.get();
   if (!file1 || c <= 0 || c != '(') {
      Error("TTabCom::DetermineClass", "variable \"%s\" not defined?",
            varName);
      goto cleanup;
   }
   IfDebug(std::cerr << (char) c << std::flush);

   // in case of success, "class TClassName*)0x12345" remains or TClassName 0x[0-9][a-f]+
   // in case of objects
   // since the opening '(' was removed.

   if (type == "const" || type == "class") {
      file1 >> type;
      // ignore ' '
      c = file1.get();
      IfDebug(std::cerr << (char) c << std::flush);
   }

   // this is what we want
   type.ReadToDelim(file1, ')');
   IfDebug(std::cerr << type << std::endl);

   // new version of CINT returns: "class TClassName*const)0x12345"
   // so we have to strip off "const"
   if (type.EndsWith("const"))
      type.Remove(type.Length() - 5);

cleanup:
   // done reading from file
   file1.close();
   gSystem->Unlink(outf);

   return type;
}

////////////////////////////////////////////////////////////////////////////////
///[static utility function]/////////////////////////////
///
///  returns true iff "s" ends with one of
///  the strings listed in the "TabCom.FileIgnore" resource.
///
//////////////////////////////////////////////////////////////

Bool_t TTabCom::ExcludedByFignore(TString s)
{
   const char *fignore = gEnv->GetValue("TabCom.FileIgnore", (char *) 0);

   if (!fignore) {
      return kFALSE;
   } else {
#ifdef R__SSTREAM
      std::istringstream endings((char *) fignore);
#else
      std::istrstream endings((char *) fignore);  // do i need to make a copy first?
#endif
      TString ending;

      ending.ReadToDelim(endings, kDelim);

      while (!ending.IsNull()) {
         if (s.EndsWith(ending))
            return kTRUE;
         else
            ending.ReadToDelim(endings, kDelim);  // next
      }
      return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
///[static utility function]/////////////////////////////
///
///  returns a colon-separated string of directories
///  that CINT will search when you call #include<...>
///
///  returns empty string on failure.
///
////////////////////////////////////////////////////////////

TString TTabCom::GetSysIncludePath()
{
   // >i noticed that .include doesn't list the standard directories like
   // >/usr/include or /usr/include/CC.
   // >
   // >how can i get a list of all the directories the interpreter will
   // >search through when the user does a #include<...> ?
   //
   // Right now, there is no easy command to tell you about it.  Instead, I can
   // describe it here.
   //
   // 1) CINT first searches current working directory for #include "xxx"
   //   (#include <xxx> does not)
   //
   // 2) CINT searches include path directories given by -I option
   //
   // 3) CINT searches following standard include directories.
   //    $CINTSYSDIR/include
   //    $CINTSYSDIR/stl
   //    $CINTSYSDIR/msdev/include   if VC++4.0
   //    $CINTSYSDIR/sc/include      if Symantec C++
   //    /usr/include
   //    /usr/include/g++            if gcc,g++
   //    /usr/include/CC             if HP-UX
   //    /usr/include/codelibs       if HP-UX
   //
   // .include command only displays 2).
   //
   // Thank you
   // Masaharu Goto

   // 1) current dir
   // ----------------------------------------------
   // N/A


   // 2) -I option (and #pragma includepath)
   // ----------------------------------------------

   // get this part of the include path from the interpreter
   // and stick it in a tmp file.
   TString outf = ".TTabCom-";
   FILE *fout = gSystem->TempFileName(outf);
   if (!fout) return "";
   gCling->DisplayIncludePath(fout);
   fclose(fout);

   // open the tmp file
   std::ifstream file1(outf);
   if (!file1) {                // error
      Error("TTabCom::GetSysIncludePath", "could not open file \"%s\"",
            outf.Data());
      gSystem->Unlink(outf);
      return "";
   }
   // parse it.
   TString token;               // input buffer
   TString path;                // all directories so far (colon-separated)
   file1 >> token;              // skip "include"
   file1 >> token;              // skip "path:"
   while (file1) {
      file1 >> token;
      if (!token.IsNull()) {
         if (path.Length() > 0)
            path.Append(":");
         path.Append(token.Data() + 2);  // +2 skips "-I"
      }
   }

   // done with the tmp file
   file1.close();
   gSystem->Unlink(outf);

   // 3) standard directories
   // ----------------------------------------------

#ifndef CINTINCDIR
   TString sCINTSYSDIR("$ROOTSYS/cint");
#else
   TString sCINTSYSDIR(CINTINCDIR);
#endif
   path.Append(":" + sCINTSYSDIR + "/include");
//   path.Append(":"+CINTSYSDIR+"/stl");
//   path.Append(":"+CINTSYSDIR+"/msdev/include");
//   path.Append(":"+CINTSYSDIR+"/sc/include");
   path.Append(":/usr/include");
//   path.Append(":/usr/include/g++");
//   path.Append(":/usr/include/CC");
//   path.Append(":/usr/include/codelibs");

   return path;
}

////////////////////////////////////////////////////////////////////////////////
///[static utility function]/////////////////////////////
///
///  calls TSystem::GetPathInfo() to see if "fileName"
///  is a system directory.
///
////////////////////////////////////////////////////////

Bool_t TTabCom::IsDirectory(const char fileName[])
{
   FileStat_t stat;
   if (!gSystem->GetPathInfo(fileName, stat))
      return R_ISDIR(stat.fMode);
   else
      return false;
}

////////////////////////////////////////////////////////////////////////////////
///[static utility function]/////////////////////////////
///
///  creates a list containing the full path name for each file
///  in the (colon separated) string "path1"
///
///  memory is allocated with "new", so
///  whoever calls this function takes responsibility for deleting it.
///
///////////////////////////////////////////////////////////////////////

TSeqCollection *TTabCom::NewListOfFilesInPath(const char path1[])
{
   assert(path1 != 0);
   if (!path1[0]) path1 = ".";

   TContainer *pList = new TContainer;  // maybe use RTTI here? (since its a static function)
#ifdef R__SSTREAM
   std::istringstream path((char *) path1);
#else
   std::istrstream path((char *) path1);
#endif

   while (path.good())
   {
      TString dirName;
      dirName.ReadToDelim(path, kDelim);
      if (dirName.IsNull())
         continue;

      IfDebug(std::cerr << "NewListOfFilesInPath(): dirName = " << dirName <<
              std::endl);

      AppendListOfFilesInDirectory(dirName, pList);
   }

   return pList;
}

////////////////////////////////////////////////////////////////////////////////
///[static utility function]/////////////////////////////
///
///  true if "fileName"
///  1. is an absolute path ("/tmp/a")
///  2. is a relative path  ("../whatever", "./test")
///  3. starts with user name ("~/mail")
///  4. starts with an environment variable ("$ROOTSYS/bin")
///
///////////////////////////////////////////////////////////////////////////

Bool_t TTabCom::PathIsSpecifiedInFileName(const TString & fileName)
{
   char c1 = (fileName.Length() > 0) ? fileName[0] : 0;
   return c1 == '/' || c1 == '~' || c1 == '$' || fileName.BeginsWith("./")
       || fileName.BeginsWith("../");
}

////////////////////////////////////////////////////////////////////////////////
///[static utility function]/////////////////////////////
///
///  calling "NoMsg( errorLevel )",
///  sets "gErrorIgnoreLevel" to "errorLevel+1" so that
///  all errors with "level < errorLevel" will be ignored.
///
///  calling the function with a negative argument
///  (e.g., "NoMsg( -1 )")
///  resets gErrorIgnoreLevel to its previous value.
///
///////////////////////////////////////////////////////////////////

void TTabCom::NoMsg(Int_t errorLevel)
{
   ////////////////////////////////////////////////////////////////
   //
   // if you call the function twice with a non-negative argument
   // (without an intervening call with a negative argument)
   // it will complain because it is almost certainly an error
   // that will cause the function to loose track of the previous
   // value of gErrorIgnoreLevel.
   //
   // most common causes: 1. suspiciously placed "return;" statement
   //                     2. calling a function that calls "NoMsg()"
   //
   //////////////////////////////////////////////////////////////////

   const Int_t kNotDefined = -2;
   static Int_t old_level = kNotDefined;

   if (errorLevel < 0)          // reset
   {
      if (old_level == kNotDefined) {
         std::cerr << "NoMsg(): ERROR 1. old_level==" << old_level << std::endl;
         return;
      }

      gErrorIgnoreLevel = old_level;  // restore
      old_level = kNotDefined;
   } else                       // set
   {
      if (old_level != kNotDefined) {
         std::cerr << "NoMsg(): ERROR 2. old_level==" << old_level << std::endl;
         return;
      }

      old_level = gErrorIgnoreLevel;
      if (gErrorIgnoreLevel <= errorLevel)
         gErrorIgnoreLevel = errorLevel + 1;
   }
}

//
//                           static utility functions
//
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
//
//                       private member functions
//
//

////////////////////////////////////////////////////////////////////////////////
/// [private]

Int_t TTabCom::Complete(const TRegexp & re,
                        const TSeqCollection * pListOfCandidates,
                        const char appendage[],
                        std::ostream& out,
                        TString::ECaseCompare cmp)
{
   // returns position of first change in buffer
   // ------------------------------------------
   // -2 ==> new line altogether (whole thing needs to be redrawn, including prompt)
   // -1 ==> no changes
   //  0 ==> beginning of line
   //  1 ==> after 1st char
   //  n ==> after nth char

   IfDebug(std::cerr << "TTabCom::Complete() ..." << std::endl);
   assert(fpLoc != 0);
   assert(pListOfCandidates != 0);

   Int_t pos = 0;               // position of first change
   const int loc = *fpLoc;      // location where TAB was pressed

   // -----------------------------------------
   //
   // 1. get the substring we need to complete
   //
   // NOTES:
   // s1 = original buffer
   // s2 = sub-buffer from 0 to wherever the user hit TAB
   // s3 = the actual text that needs completing
   //
   // -----------------------------------------
   TString s1(fBuf);
   TString s2 = s1(0, loc);
   TString s3 = s2(re);

   int start = s2.Index(re);

   IfDebug(std::cerr << "   s1: " << s1 << std::endl);
   IfDebug(std::cerr << "   s2: " << s2 << std::endl);
   IfDebug(std::cerr << "   s3: " << s3 << std::endl);
   IfDebug(std::cerr << "start: " << start << std::endl);
   IfDebug(std::cerr << std::endl);

   // -----------------------------------------
   // 2. go through each possible completion,
   //    keeping track of the number of matches
   // -----------------------------------------
   TList listOfMatches;         // list of matches (local filenames only) (insertion order must agree across these 3 lists)
   TList listOfFullPaths;       // list of matches (full filenames)       (insertion order must agree across these 3 lists)
   listOfMatches.SetOwner();
   listOfFullPaths.SetOwner();

   int nMatches = 0;            // number of matches
   TObject *pObj;               // pointer returned by iterator
   TIter next_candidate(pListOfCandidates);
   TIter next_match(&listOfMatches);
   TIter next_fullpath(&listOfFullPaths);

   // stick all matches into "listOfMatches"
   while ((pObj = next_candidate())) {
      // get the full filename
      const char *s4 = pObj->GetName();

      assert(s4 != 0);

      // pick off tail
      const char *s5 = strrchr(s4, '/');
      if (!s5)
         s5 = s4;               // no '/' found
      else
         s5 += 1;               // advance past '/'

      // if case sensitive (normal behaviour), check for match
      // if case insensitive, convert to TString and compare case insensitively
      if ((cmp == TString::kExact) && (strstr(s5, s3) == s5)) {
         nMatches += 1;
         listOfMatches.Add(new TObjString(s5));
         listOfFullPaths.Add(new TObjString(s4));
         IfDebug(std::cerr << "adding " << s5 << '\t' << s4 << std::endl);
      } else if (cmp == TString::kIgnoreCase) {
         TString ts5(s5);
         if (ts5.BeginsWith(s3, cmp))
         {
            nMatches += 1;
            listOfMatches.Add(new TObjString(s5));
            listOfFullPaths.Add(new TObjString(s4));
            IfDebug(std::cerr << "adding " << s5 << '\t' << s4 << std::endl);
         }
      } else {
//rdm         IfDebug(std::cerr << "considered " << s5 << '\t' << s4 << std::endl);
      }

   }

   // -----------------------------------------
   // 3. beep, list, or complete
   //    depending on how many matches were found
   // -----------------------------------------

   // 3a. no matches ==> bell
   TString partialMatch = "";

   if (nMatches == 0) {
      // Ring a bell!
      gSystem->Beep();
      pos = -1;
      goto done;                //* RETURN *//
   }
   // 3b. one or more matches.
   char match[1024];

   if (nMatches == 1) {
      // get the (lone) match
      const char *short_name = next_match()->GetName();
      const char *full_name = next_fullpath()->GetName();

      pObj = pListOfCandidates->FindObject(short_name);
      if (pObj) {
         IfDebug(std::cerr << std::endl << "class: " << pObj->ClassName() << std::endl);
         TString className = pObj->ClassName();
         if (0);
         else if (className == "TMethod" || className == "TFunction") {
            TFunction *pFunc = (TFunction *) pObj;
            if (0 == pFunc->GetNargs())
               appendage = "()";  // no args
            else
               appendage = "("; // user needs to supply some args
         } else if (className == "TDataMember") {
            appendage = " ";
         }
      }

      CopyMatch(match, short_name, appendage, full_name);
   } else {
      // multiple matches ==> complete as far as possible
      Char_t ch;
      Int_t nGoodStrings;

      for (int i = 0;
           (ch = AllAgreeOnChar(i, &listOfMatches, nGoodStrings));
           i += 1) {
         IfDebug(std::cerr << " i=" << i << " ch=" << ch << std::endl);
         partialMatch.Append(ch);
      }

      const char *s;
      const char *s0;

      // multiple matches, but maybe only 1 of them is any good.
      if (nGoodStrings == 1) {

         // find the 1 good match
         do {
            s = next_match()->GetName();
            s0 = next_fullpath()->GetName();
         }
         while (ExcludedByFignore(s));

         // and use it.
         CopyMatch(match, s, appendage, s0);
      } else {
         IfDebug(std::cerr << "more than 1 GoodString" << std::endl);

         if (partialMatch.Length() > s3.Length())
            // this partial match is our (partial) completion.
         {
            CopyMatch(match, partialMatch.Data());
         } else
            // couldn't do any completing at all,
            // print a list of all the ambiguous matches
            // (except for those excluded by "FileIgnore")
         {
            IfDebug(std::cerr << "printing ambiguous matches" << std::endl);
            std::set<std::string> alreadyPrinted;
            while ((pObj = next_match())) {
               s = pObj->GetName();
               if (alreadyPrinted.insert(s).second) {
                  // Nothing of that name has been printed yet.
                  s0 = next_fullpath()->GetName();
                  if (!ExcludedByFignore(s) || nGoodStrings == 0) {
                     if (IsDirectory(s0))
                        out << s << "/" << std::endl;
                     else
                        out << s << std::endl;
                  }
               }
            }
            pos = -2;
            if (cmp == TString::kExact || partialMatch.Length() < s3.Length()) {
               goto done;          //* RETURN *//
            } // else:
            // update the matching part, will have changed
            // capitalization because only cmp == TString::kIgnoreCase
            // matches.
            CopyMatch(match, partialMatch.Data());
         }
      }
   }


   // ---------------------------------------
   // 4. finally write text into the buffer.
   // ---------------------------------------
   {
      int i = strlen(fBuf);     // old EOL position is i
      int l = strlen(match) - (loc - start);  // new EOL position will be i+L

      // first check for overflow
      if (strlen(fBuf) + strlen(match) + 1 > BUF_SIZE) {
         Error("TTabCom::Complete", "buffer overflow");
         pos = -2;
         goto done;             /* RETURN */
      }
      // debugging output
      IfDebug(std::cerr << "  i=" << i << std::endl);
      IfDebug(std::cerr << "  L=" << l << std::endl);
      IfDebug(std::cerr << "loc=" << loc << std::endl);

      // slide everything (including the null terminator) over to make space
      for (; i >= loc; i -= 1) {
         fBuf[i + l] = fBuf[i];
      }

      // insert match
      strlcpy(fBuf + start, match, BUF_SIZE - start);

      // the "get"->"Get" case of TString::kIgnore sets pos to -2
      // and falls through to update the buffer; we need to return
      // -2 in that case, so check here:
      if (pos != -2) {
         pos = loc;                // position of first change in "fBuf"
         if (cmp == TString::kIgnoreCase && pos < 0) {
            // We might have changed somthing before loc, due to differences in
            // capitalization. So return start:
            pos = start;
         }
      }
      *fpLoc = loc + l;         // new cursor position
   }

done:                         // <----- goto label
   // un-init
   fpLoc = 0;
   fBuf = 0;

   return pos;
}

////////////////////////////////////////////////////////////////////////////////
/// [private]

void TTabCom::CopyMatch(char dest[], const char localName[],
                        const char appendage[],
                        const char fullName[]) const
{
   // if "appendage" is 0, no appendage is applied.
   //
   // if "appendage" is of the form "filenameXXX" then,
   // "filename" is ignored and "XXX" is taken to be the appendage,
   // but it will only be applied if the file is not a directory...
   // if the file is a directory, a "/" will be used for the appendage instead.
   //
   // if "appendage" is of the form "XXX" then "XXX" will be appended to the match.

   assert(dest != 0);
   assert(localName != 0);

   // potential buffer overflow.
   strcpy(dest, localName);

   const char *key = "filename";
   const int key_len = strlen(key);

   IfDebug(std::cerr << "CopyMatch()." << std::endl);
   IfDebug(std::cerr << "localName: " << (localName ? localName : "0") <<
           std::endl);
   IfDebug(std::cerr << "appendage: " << (appendage ? appendage : "0") <<
           std::endl);
   IfDebug(std::cerr << " fullName: " << (fullName ? fullName : "0") <<
           std::endl);


   // check to see if "appendage" starts with "key"
   if (appendage && strncmp(appendage, key, key_len) == 0) {
      // filenames get special treatment
      appendage += key_len;
      IfDebug(std::cerr << "new appendage: " << appendage << std::endl);
      if (IsDirectory(fullName)) {
         if (fullName)
            strcpy(dest + strlen(localName), "/");
      } else {
         if (appendage)
            strcpy(dest + strlen(localName), appendage);
      }
   } else {
      if (appendage)
         strcpy(dest + strlen(localName), appendage);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// [private]

TTabCom::EContext_t TTabCom::DetermineContext() const
{
   assert(fBuf != 0);

   const char *pStart;          // start of match
   const char *pEnd;            // end of match

   for (int context = 0; context < kNUM_PAT; ++context) {
      pEnd = Matchs(fBuf, *fpLoc, fPat[context], &pStart);
      if (pEnd) {
         IfDebug(std::cerr << std::endl
                 << "context=" << context << " "
                 << "RegExp=" << fRegExp[context]
                 << std::endl);
         return EContext_t(context);  //* RETURN *//
      }
   }

   return kUNKNOWN_CONTEXT;     //* RETURN *//
}

////////////////////////////////////////////////////////////////////////////////
/// [private]

TString TTabCom::DeterminePath(const TString & fileName,
                               const char defaultPath[]) const
{
   if (PathIsSpecifiedInFileName(fileName)) {
      TString path = fileName;
      gSystem->ExpandPathName(path);
      Int_t end = path.Length()-1;
      if (end>0 && path[end]!='/' && path[end]!='\\') {
         path = gSystem->GetDirName(path);
      }
      return path;
   } else {
      TString newBase, extendedPath;
      if (fileName.Contains("/")) {
         Int_t end = fileName.Length()-1;
         if (fileName[end] != '/' && fileName[end] != '\\') {
            newBase = gSystem->GetDirName(fileName);
         } else {
            newBase = fileName;
         }
         extendedPath = ExtendPath(defaultPath, newBase);
      } else {
         newBase = "";
         extendedPath = defaultPath;
      }
      IfDebug(std::cerr << std::endl);
      IfDebug(std::cerr << "    fileName: " << fileName << std::endl);
      IfDebug(std::cerr << "    pathBase: " << newBase << std::endl);
      if (defaultPath) {
         IfDebug(std::cerr << " defaultPath: " << defaultPath << std::endl);
      } else {
         IfDebug(std::cerr << " defaultPath: " << std::endl);
      }
      IfDebug(std::cerr << "extendedPath: " << extendedPath << std::endl);
      IfDebug(std::cerr << std::endl);

      return extendedPath;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// [private]

TString TTabCom::ExtendPath(const char originalPath[], TString newBase) const
{
   if (newBase.BeginsWith("/"))
      newBase.Remove(TString::kLeading, '/');
#ifdef R__SSTREAM
   std::stringstream str;
#else
   std::strstream str;
#endif
   TString dir;
   TString newPath;
   if (originalPath) str << originalPath;

   while (str.good())
   {
      dir = "";
      dir.ReadToDelim(str, kDelim);
      if (dir.IsNull())
         continue;              // ignore blank entries
      newPath.Append(dir);
      if (!newPath.EndsWith("/"))
         newPath.Append("/");
      newPath.Append(newBase);
      newPath.Append(kDelim);
   }

   return newPath.Strip(TString::kTrailing, kDelim);
}

////////////////////////////////////////////////////////////////////////////////
/// [private]

Int_t TTabCom::Hook(char *buf, int *pLoc, std::ostream& out)
{
   // initialize
   fBuf = buf;
   fpLoc = pLoc;

   // frodo: iteration counter for recursive MakeClassFromVarName
   fLastIter = 0;

   // default
   Int_t pos = -2;  // position of the first character that was changed in the buffer (needed for redrawing)

   // get the context this tab was triggered in.
   EContext_t context = DetermineContext();

   // get the substring that triggered this tab (as defined by "SetPattern()")
   const char dummy[] = ".";
   TRegexp re1(context == kUNKNOWN_CONTEXT ? dummy : fRegExp[context]);
   TString s1(fBuf);
   TString s2 = s1(0, *fpLoc);
   TString s3 = s2(re1);

   switch (context) {
   case kUNKNOWN_CONTEXT:
      std::cerr << std::endl << "tab completion not implemented for this context" <<
          std::endl;
      pos = -2;
      break;

   case kSYS_UserName:
      {
         const TSeqCollection *pListOfUsers = GetListOfUsers();

         pos = Complete("[^~]*$", pListOfUsers, "/", out);
      }
      break;
   case kSYS_EnvVar:
      {
         const TSeqCollection *pEnv = GetListOfEnvVars();

         pos = Complete("[^$]*$", pEnv, "", out);
      }
      break;

   case kCINT_stdout:
   case kCINT_stderr:
   case kCINT_stdin:
      {
         const TString fileName = s3("[^ ><]*$");
         const TString filePath = DeterminePath(fileName,0);
         const TSeqCollection *pListOfFiles =
             GetListOfFilesInPath(filePath.Data());

//             pos = Complete( "[^ /]*$", pListOfFiles, " ", out );
         pos = Complete("[^ /]*$", pListOfFiles, "filename ", out);
      }
      break;

   case kCINT_Edit:
   case kCINT_Load:
   case kCINT_Exec:
   case kCINT_EXec:
      {
         const TString fileName = s3("[^ ]*$");
         const TString macroPath =
             DeterminePath(fileName, TROOT::GetMacroPath());
         const TSeqCollection *pListOfFiles =
             GetListOfFilesInPath(macroPath.Data());

//             pos = Complete( "[^ /]*$", pListOfFiles, " ", out);
         pos = Complete("[^ /]*$", pListOfFiles, "filename ", out);
      }
      break;

   case kCINT_pragma:
      {
         pos = Complete("[^ ]*$", GetListOfPragmas(), "", out);
      }
      break;
   case kCINT_includeSYS:
      {
         TString fileName = s3("[^<]*$");
         if (PathIsSpecifiedInFileName(fileName) || fileName.Contains("/")) {
            TString includePath =
                DeterminePath(fileName, GetSysIncludePath());

//                  pos = Complete( "[^</]*$", GetListOfFilesInPath( includePath ), "> ", out);
            pos =
                Complete("[^</]*$", GetListOfFilesInPath(includePath),
                         "filename> ", out);
         } else {
//                  pos = Complete( "[^</]*$", GetListOfSysIncFiles(), "> ", out);
            pos =
                Complete("[^</]*$", GetListOfSysIncFiles(), "filename> ", out);
         }
      }
      break;
   case kCINT_includePWD:
      {
         const TString fileName = s3("[^\"]*$");
         const TString includePath = DeterminePath(fileName, ".");
         const TSeqCollection *pListOfFiles =
             GetListOfFilesInPath(includePath.Data());

//             pos = Complete( "[^\"/]*$", pListOfFiles, "\" ", out);
         pos = Complete("[^\"/]*$", pListOfFiles, "filename\" ", out);
      }
      break;

   case kCINT_cpp:
      {
         pos = Complete("[^# ]*$", GetListOfCppDirectives(), " ", out);
      }
      break;

   case kROOT_Load:
      {
         const TString fileName = s3("[^\"]*$");
         const TString dynamicPath = DeterminePath(fileName, gSystem->GetDynamicPath());
         const TSeqCollection *pListOfFiles = GetListOfFilesInPath(dynamicPath);

//             pos = Complete( "[^\"/]*$", pListOfFiles, "\");", out);
         pos = Complete("[^\"/]*$", pListOfFiles, "filename\");", out);
      }
      break;

   case kSYS_FileName:
      {
         const TString fileName = s3("[^ \"]*$");
         const TString filePath = DeterminePath(fileName,".");
         const TSeqCollection *pListOfFiles = GetListOfFilesInPath(filePath.Data());

         pos = Complete("[^\" /]*$", pListOfFiles, "filename\"", out);
      }
      break;

   case kCXX_ScopeMember:
      {
         const EContext_t original_context = context;  // save this for later

         TClass *pClass;
         // may be a namespace, class, object, or pointer
         TString name = s3("^[_a-zA-Z][_a-zA-Z0-9]*");

         IfDebug(std::cerr << std::endl);
         IfDebug(std::cerr << "name: " << '"' << name << '"' << std::endl);

         // We need to decompose s3 a little more:
         // The part name is the partial symbol at the end of ::
         // eg. given s3 = "foo::bar::part" ,  partname = "part"
         TString partname = s3("[_a-zA-Z][_a-zA-Z0-9]*$");

         // The prefix, considering the s3 = "foo::bar::part" example would be
         // prefix = "foo::bar::". prefix equals the empty string if there is only one
         // or no set of colons in s3.
         // Note: we reconstruct the fully qualified name with a while loop because
         // it does not seem that TRegexp can handle something like "([_a-zA-Z][_a-zA-Z0-9]*::)+$"
         TString prefix = "";
         TString str = s2;
         str.Remove(str.Length() - partname.Length(), partname.Length());
         while (1) {
            TString sym = str("[_a-zA-Z][_a-zA-Z0-9]*::$");
            if (sym.Length() == 0)
               break;
            str.Remove(str.Length() - sym.Length(), sym.Length());
            prefix = sym + prefix;
         }

         // Not the preprefix would be = "foo::" from our previous example or the empty
         // string, "" if there is only one or no set of colons in prefix, eg. prefix = "bar::"
         TString preprefix = prefix;
         TString sym = prefix("[_a-zA-Z][_a-zA-Z0-9]*::$");
         preprefix.Remove(preprefix.Length() - sym.Length(), sym.Length());

         IfDebug(std::cerr << "prefix: " << '"' << prefix << '"' << std::endl);
         IfDebug(std::cerr << "preprefix: " << '"' << preprefix << '"' << std::endl);

         TString namesp = prefix;
         if (namesp.Length() >= 2)
            namesp.Remove(namesp.Length() - 2, 2);  // Remove the '::' at the end of the string.
         IfDebug(std::cerr << "namesp: " << '"' << namesp << '"' << std::endl);

         // Make sure autoloading happens (if it can).
         TryMakeClassFromClassName(namesp);

         TContainer *pList = new TContainer;
         // Add all classes to pList that contain the prefix, i.e. are in the
         // specified namespace.
         const TSeqCollection *tmp = GetListOfClasses();
         if (!tmp) break;

         Int_t i;
         for (i = 0; i < tmp->GetSize(); i++) {
            TString astr = ((TObjString *) tmp->At(i))->String();
            TString rxp = "^";
            rxp += prefix;
            if (astr.Contains(TRegexp(rxp))) {
               astr.Remove(0, prefix.Length());
               TString s = astr("^[^: ]*");
               TObjString *ostr = new TObjString(s);
               if (!pList->Contains(ostr))
                  pList->Add(ostr);
               else
                  delete ostr;
            }
         }

         // If a class with the same name as the Namespace name exists then
         // add it to the pList. (I don't think the C++ spec allows for this
         // but do this anyway, cant harm).
         pClass = TryMakeClassFromClassName(preprefix + name);
         if (pClass) {
            pList->AddAll(pClass->GetListOfAllPublicMethods(true));
            pList->AddAll(pClass->GetListOfAllPublicDataMembers(true));
         }

         pos = Complete("[^: ]*$", pList, "", out);

         delete pList;

         if (context != original_context)
            pos = -2;
      }
      break;

   case kCXX_DirectMember:
   case kCXX_IndirectMember:
      {
         const EContext_t original_context = context;  // save this for later

         TClass *pClass;

         // frodo: Instead of just passing the last portion of the string to
         //        MakeClassFromVarName(), we now pass the all string and let
         //        it decide how to handle it... I know it's not the best way
         //        because of the context handling, but I wanted to "minimize"
         //        the changes to the current code and this seemed the best way
         //        to do it
         TString name = s1("[_a-zA-Z][-_a-zA-Z0-9<>():.]*$");

         IfDebug(std::cerr << std::endl);
         IfDebug(std::cerr << "name: " << '"' << name << '"' << std::endl);

         switch (context) {
         case kCXX_DirectMember:
            pClass = MakeClassFromVarName(name, context);
            break;
         case kCXX_IndirectMember:
            pClass = MakeClassFromVarName(name, context);
            break;
         default:
            Fatal("TTabCom::Hook","Conext case %d not handled",context);
            return 0; // Avoid warning about uninitialized pClass.
         }
         if (!pClass) {
            pos = -2;
            break;
         }

         TContainer *pList = new TContainer;

         pList->AddAll(pClass->GetListOfAllPublicMethods());
         pList->AddAll(pClass->GetListOfAllPublicDataMembers());

         switch (context) {
         case kCXX_DirectMember:
            {
               int* store_fpLoc = fpLoc;
               char* store_fBuf = fBuf;
               pos = Complete("[^. ]*$", pList, "(", out);
               if (pos == -1) {
                  fpLoc = store_fpLoc;
                  fBuf = store_fBuf;
                  pos = Complete("[^. ]*$", pList, "(", out, TString::kIgnoreCase);
               }
               break;
            }
         case kCXX_IndirectMember:
            pos = Complete("[^> ]*$", pList, "(", out);
            break;
         default:
            Fatal("TTabCom::Hook","Conext case %d not handled",context);
            return 0; // Avoid warning about uninitialized pClass.
         }

         delete pList;

         if (context != original_context)
            pos = -2;
      }
      break;

   case kCXX_ScopeProto:
      {
         const EContext_t original_context = context;  // save this for later

         // get class
         TClass *pClass;
         TString name = s3("^[_a-zA-Z][_a-zA-Z0-9]*");
         // "name" may now be the name of a class, object, or pointer

         IfDebug(std::cerr << std::endl);
         IfDebug(std::cerr << "name: " << '"' << name << '"' << std::endl);

         // We need to decompose s3 a little more:
         // The partname is the method symbol and a bracket at the end of ::
         // eg. given s3 = "foo::bar::part(" ,  partname = "part("
         TString partname = s3("[_a-zA-Z][_a-zA-Z0-9]* *($");

         // The prefix, considering the s3 = "foo::bar::part" example would be
         // prefix = "foo::bar::". prefix equals the empty string if there is only one
         // or no set of colons in s3.
         // Note: we reconstruct the fully qualified name with a while loop because
         // it does not seem that TRegexp can handle something like "([_a-zA-Z][_a-zA-Z0-9]*::)+$"
         TString prefix = "";
         TString str = s2;
         str.Remove(str.Length() - partname.Length(), partname.Length());
         while (1) {
            TString sym = str("[_a-zA-Z][_a-zA-Z0-9]*::$");
            if (sym.Length() == 0)
               break;
            str.Remove(str.Length() - sym.Length(), sym.Length());
            prefix = sym + prefix;
         }

         // Not the preprefix would be = "foo::" from our previous example or the empty
         // string, "" if there is only one or no set of colons in prefix, eg. prefix = "bar::"
         TString preprefix = prefix;
         TString sym = prefix("[_a-zA-Z][_a-zA-Z0-9]*::$");
         preprefix.Remove(preprefix.Length() - sym.Length(), sym.Length());

         IfDebug(std::cerr << "prefix: " << '"' << prefix << '"' << std::endl);
         IfDebug(std::cerr << "preprefix: " << '"' << preprefix << '"' << std::endl);

         pClass = MakeClassFromClassName(preprefix + name);
         if (!pClass) {
            pos = -2;
            break;
         }
         // get method name
         TString methodName;

         // (normal member function)
         methodName = s3("[^:>\\.(]*($");
         methodName.Chop();
         methodName.Remove(TString::kTrailing, ' ');

         IfDebug(std::cerr << methodName << std::endl);

         // get methods
         TContainer *pList = new TContainer;
         pList->AddAll(pClass->GetListOfAllPublicMethods());

         // print prototypes
         Bool_t foundOne = kFALSE;
         TIter nextMethod(pList);
         TMethod *pMethod;
         while ((pMethod = (TMethod *) nextMethod())) {
            if (methodName == pMethod->GetName()) {
               foundOne = kTRUE;
               out << pMethod->GetReturnTypeName()
                   << " " << pMethod->GetName()
                   << pMethod->GetSignature();
               const char *comment = pMethod->GetCommentString();
               if (comment && comment[0] != '\0') {
                  out << " \t// " << comment;
               }
               out << std::endl;
            }
         }

         // done
         if (foundOne) {
            pos = -2;
         } else {
            gSystem->Beep();
            pos = -1;
         }

         // cleanup
         delete pList;

         if (context != original_context)
            pos = -2;
      }
      break;

   case kCXX_DirectProto:
   case kCXX_IndirectProto:
   case kCXX_NewProto:
   case kCXX_ConstructorProto:
      {
         const EContext_t original_context = context;  // save this for later

         // get class
         TClass *pClass;
         TString name;
         if (context == kCXX_NewProto) {
            name = s3("[_a-zA-Z][_a-zA-Z0-9:]* *($", 3);
            name.Chop();
            name.Remove(TString::kTrailing, ' ');
            // "name" should now be the name of a class
         } else {
            name = s3("^[_a-zA-Z][_a-zA-Z0-9:]*");
            // "name" may now be the name of a class, object, or pointer
         }
         IfDebug(std::cerr << std::endl);
         IfDebug(std::cerr << "name: " << '"' << name << '"' << std::endl);

         // frodo: Again, passing the all string
         TString namerec = s1;

         switch (context) {
         case kCXX_ScopeProto:
            pClass = MakeClassFromClassName(name);
            break;
         case kCXX_DirectProto:
            pClass = MakeClassFromVarName(namerec, context); // frodo
            break;
         case kCXX_IndirectProto:
            pClass = MakeClassFromVarName(namerec, context); // frodo
            break;
         case kCXX_NewProto:
            pClass = MakeClassFromClassName(name);
            break;
         case kCXX_ConstructorProto:
            pClass = MakeClassFromClassName(name);
            break;
         default:
            Fatal("TTabCom::Hook","Conext case %d not handled",context);
            return 0; // Avoid warning about uninitialized pClass.
         }
         if (!pClass) {
            pos = -2;
            break;
         }
         // get method name
         TString methodName;
         if (context == kCXX_ConstructorProto || context == kCXX_NewProto) {
            // (constructor)
            methodName = name("[_a-zA-Z][_a-zA-Z0-9]*$");
         } else {
            // (normal member function)
            methodName = s3("[^:>\\.(]*($");
            methodName.Chop();
            methodName.Remove(TString::kTrailing, ' ');
         }
         IfDebug(std::cerr << methodName << std::endl);

         // get methods
         TContainer *pList = new TContainer;
         pList->AddAll(pClass->GetListOfAllPublicMethods());

         // print prototypes
         Bool_t foundOne = kFALSE;
         TIter nextMethod(pList);
         TMethod *pMethod;
         while ((pMethod = (TMethod *) nextMethod())) {
            if (methodName == pMethod->GetName()) {
               foundOne = kTRUE;
               out << pMethod->GetReturnTypeName()
                   << " " << pMethod->GetName()
                   << pMethod->GetSignature();
               const char *comment = pMethod->GetCommentString();
               if (comment && comment[0] != '\0') {
                  out << " \t// " << comment;
               }
               out << std::endl;
            }
         }

         // done
         if (foundOne) {
            pos = -2;
         } else {
            gSystem->Beep();
            pos = -1;
         }

         // cleanup
         delete pList;

         if (context != original_context)
            pos = -2;
      }
      break;

   case kCXX_Global:
      {
         // first need to veto a few possibilities.
         int l2 = s2.Length(), l3 = s3.Length();

         // "abc().whatever[TAB]"
         if (l2 > l3 && s2[l2 - l3 - 1] == '.') {
            std::cerr << std::endl <<
                "tab completion not implemented for this context" << std::endl;
            break;              // veto
         }
         // "abc()->whatever[TAB]"
         if (l2 > l3 + 1 && s2(l2 - l3 - 2, 2) == "->") {
            std::cerr << std::endl <<
                "tab completion not implemented for this context" << std::endl;
            break;              // veto
         }

         TContainer *pList = new TContainer;

         const TSeqCollection *pL2 = GetListOfClasses();
         if (pL2) pList->AddAll(pL2);

         //
         const TSeqCollection *pC1 = GetListOfGlobals();
         if (pC1) pList->AddAll(pC1);
         //
         TCollection *pC3 = GetListOfGlobalFunctions();
         if (pC3) pList->AddAll(pC3);

         pos = Complete("[_a-zA-Z][_a-zA-Z0-9]*$", pList, "", out);

         delete pList;
      }
      break;

   case kCXX_GlobalProto:
      {
         // get function name
         TString functionName = s3("[_a-zA-Z][_a-zA-Z0-9]*");
         IfDebug(std::cerr << functionName << std::endl);

         TContainer listOfMatchingGlobalFuncs;
         TIter nextGlobalFunc(GetListOfGlobalFunctions());
         TObject *pObj;
         while ((pObj = nextGlobalFunc())) {
            if (strcmp(pObj->GetName(), functionName) == 0) {
               listOfMatchingGlobalFuncs.Add(pObj);
            }
         }

         if (listOfMatchingGlobalFuncs.IsEmpty()) {
            std::cerr << std::endl << "no such function: " << dblquote(functionName)
                << std::endl;
         } else {
            TIter next(&listOfMatchingGlobalFuncs);
            TFunction *pFunction;
            while ((pFunction = (TFunction *) next())) {
               out << pFunction->GetReturnTypeName()
                   << " " << pFunction->GetName()
                   << pFunction->GetSignature()
                   << std::endl;
            }
         }

         pos = -2;
      }
      break;

      /******************************************************************/
      /*                                                                */
      /* default: should never happen                                   */
      /*                                                                */
      /******************************************************************/
   default:
      Fatal("TTabCom::Hook","Conext case %d not handled",context);
      return 0; // Avoid warning about uninitialized pClass.
   }

   return pos;
}

////////////////////////////////////////////////////////////////////////////////
/// [private]

void TTabCom::InitPatterns()
{
   // add more patterns somewhere below.
   // add corresponding enum to "EContext_t"
   //
   // note:
   // 1. in some cases order is important ...
   //
   //    the order of the "case" statements in "switch( context )" in "TTabCom::Hook()" is Not important.
   //
   //    the order of the "SetPattern()" function calls below is Not important.
   //
   //    the order of the initializers in the "EContext_t" enumeration Is important
   //    because DetermineContext() goes through the array in order, and returns at the first match.
   //
   // 2. below, "$" will match cursor position

   SetPattern(kSYS_UserName, "~[_a-zA-Z0-9]*$");
   SetPattern(kSYS_EnvVar, "$[_a-zA-Z0-9]*$");

   SetPattern(kCINT_stdout, "; *>>?.*$");  // stdout
   SetPattern(kCINT_stderr, "; *2>>?.*$"); // stderr
   SetPattern(kCINT_stdin, "; *<.*$");     // stdin

   SetPattern(kCINT_Edit, "^ *\\.E .*$");
   SetPattern(kCINT_Load, "^ *\\.L .*$");
   SetPattern(kCINT_Exec, "^ *\\.x +[-0-9_a-zA-Z~$./]*$");
   SetPattern(kCINT_EXec, "^ *\\.X +[-0-9_a-zA-Z~$./]*$");

   SetPattern(kCINT_pragma, "^# *pragma +[_a-zA-Z0-9]*$");
   SetPattern(kCINT_includeSYS, "^# *include *<[^>]*$");   // system files
   SetPattern(kCINT_includePWD, "^# *include *\"[^\"]*$"); // local files

   SetPattern(kCINT_cpp, "^# *[_a-zA-Z0-9]*$");

   SetPattern(kROOT_Load, "gSystem *-> *Load *( *\"[^\"]*$");

   SetPattern(kCXX_NewProto, "new +[_a-zA-Z][_a-zA-Z0-9:]* *($");
   SetPattern(kCXX_ConstructorProto,
              "[_a-zA-Z][_a-zA-Z0-9:]* +[_a-zA-Z][_a-zA-Z0-9]* *($");
   SetPattern(kCXX_ScopeProto,
              "[_a-zA-Z][_a-zA-Z0-9]* *:: *[_a-zA-Z0-9]* *($");
   SetPattern(kCXX_DirectProto,
              "[_a-zA-Z][_a-zA-Z0-9()]* *\\. *[_a-zA-Z0-9]* *($");
   SetPattern(kCXX_IndirectProto,
              "[_a-zA-Z][_a-zA-Z0-9()]* *-> *[_a-zA-Z0-9]* *($");

   SetPattern(kCXX_ScopeMember,
              "[_a-zA-Z][_a-zA-Z0-9]* *:: *[_a-zA-Z0-9]*$");
   SetPattern(kCXX_DirectMember,
              "[_a-zA-Z][_a-zA-Z0-9()]* *\\. *[_a-zA-Z0-9()]*$");  //

   SetPattern(kCXX_IndirectMember,
              "[_a-zA-Z][_a-zA-Z0-9()]* *-> *[_a-zA-Z0-9()]*$");    // frodo

   SetPattern(kSYS_FileName, "\"[-0-9_a-zA-Z~$./]*$");
   SetPattern(kCXX_Global, "[_a-zA-Z][_a-zA-Z0-9]*$");
   SetPattern(kCXX_GlobalProto, "[_a-zA-Z][_a-zA-Z0-9]* *($");
}

////////////////////////////////////////////////////////////////////////////////
/// [private]
///   (does some specific error handling that makes the function unsuitable for general use.)
///   returns a new'd TClass given the name of a class.
///   user must delete.
///   returns 0 in case of error.

TClass *TTabCom::MakeClassFromClassName(const char className[]) const
{
   // the TClass constructor will print a Warning message for classes that don't exist
   // so, ignore warnings temporarily.
   NoMsg(kWarning);
   TClass *pClass = TClass::GetClass(className);
   NoMsg(-1);

   if (!pClass){
      Error("TTabCom::MakeClassFromClassName", "Unknown class \"%s\"", className);
      return nullptr;
   }

   // make sure "className" exists
   // if (pClass->Size() == 0) {   //namespace has 0 size
   if (pClass->GetListOfAllPublicMethods()->GetSize() == 0 &&
       pClass->GetListOfAllPublicDataMembers()->GetSize() == 0) {
      // i'm assuming this happens iff there was some error.
      // (misspelled the class name, for example)
      Error("TTabCom::MakeClassFromClassName", "class \"%s\" is not defined.", className);
      return 0;
   }

   return pClass;
}

////////////////////////////////////////////////////////////////////////////////
/// Same as above but does not print the error message.

TClass *TTabCom::TryMakeClassFromClassName(const char className[]) const
{
   // the TClass constructor will print a Warning message for classes that don't exist
   // so, ignore warnings temporarily.
   NoMsg(kWarning);
   TClass *pClass = TClass::GetClass(className);
   NoMsg(-1);

   return pClass;
}

////////////////////////////////////////////////////////////////////////////////
/// [private]
///   (does some specific error handling that makes the function unsuitable for general use.)
///   returns a new'd TClass given the name of a variable.
///   user must delete.
///   returns 0 in case of error.
///   if user has operator.() or operator->() backwards, will modify: context, *fpLoc and fBuf.
///   context sensitive behavior.

TClass *TTabCom::MakeClassFromVarName(const char varName[],
                                      EContext_t & context, int iter)
{
   // frodo:
   // Because of the Member and Proto recursion, this has become a bit
   // complicated, so here is how it works:
   //
   // root [1] var.a.b.c[TAB]
   //
   // will generate the sucessive calls:
   // MakeClassFromVarName("var.a.b.c", context, 0) returns the class of "c"
   // MakeClassFromVarName("var.a.b", context, 1)   returns the class of "b"
   // MakeClassFromVarName("var.a", context, 2)     returns the class of "a"
   // MakeClassFromVarName("var", context, 3)

   // need to make sure "varName" exists
   // because "DetermineClass()" prints clumsy error message otherwise.
   Bool_t varName_exists = GetListOfGlobals()->Contains(varName) || // check in list of globals first.
       (gROOT->FindObject(varName) != 0);  // then check CINT "shortcut #3"


   //
   // frodo: Member and Proto recursion code
   //
   if (0) printf("varName is [%s] with iteration [%i]\n", varName, iter);

   // ParseReverse will return 0 if there are no "." or "->" in the varName
   Int_t cut = ParseReverse(varName, strlen(varName));

   // If it's not a "simple" variable and if there is at least one "." or "->"
   if (!varName_exists && cut != 0)
   {
      TString parentName = varName;
      TString memberName = varName;

      // Check to see if this is the last call (last member/method)
      if (iter > fLastIter) fLastIter = iter;

      parentName[cut] = 0;
      if (0) printf("Parent string is [%s]\n", parentName.Data());

      // We are treating here cases like h->SetXTitle(gROOT->Get<TAB>
      // i.e. when the parentName has an unbalanced number of paranthesis.
      if (cut>2) {
         UInt_t level = 0;
         for(Int_t i = cut-1; i>=0; --i) {
            switch (parentName[i]) {
               case '(':
                  if (level) --level;
                  else {
                     parentName = parentName(i+1,cut-i-1);
                     i = 0;
                  }
                  break;
               case ')':
                  ++level; break;
            }
         }
      }

      TClass *pclass;
      // Can be "." or "->"
      if (varName[cut] == '.') {
         memberName = varName+cut+1;
         if (0) printf("Member/method is [%s]\n", memberName.Data());
         EContext_t subcontext = kCXX_DirectMember;
         pclass = MakeClassFromVarName(parentName.Data(), subcontext, iter+1);
      } else {
         memberName = varName+cut+2;
         if (0) printf("Member/method is [%s]\n", memberName.Data());
         EContext_t subcontext = kCXX_IndirectMember;
         pclass = MakeClassFromVarName(parentName.Data(), subcontext, iter+1);
      }

      if (0) printf("I got [%s] from MakeClassFromVarName()\n", pclass->GetName());

      if (pclass)
      {
         if (0) printf("Variable [%s] exists!\n", parentName.Data());

         // If it's back in the first call of the function, return immediatly
         if (iter == 0) return pclass;

         if (0) printf("Trying data member [%s] of class [%s] ...\n",
            memberName.Data(), pclass->GetName());

         // Check if it's a member
         TDataMember *dmptr = 0; //pclass->GetDataMember(memberName.Data());
         TIter   next(pclass->GetListOfAllPublicDataMembers());
         while ((dmptr = (TDataMember *) next())) {
            if (memberName == dmptr->GetName()) break;
         }
         if (dmptr)
         {
            if (0) printf("It's a member!\n");

            TString returnName = dmptr->GetTypeName();
            //              if (returnName[returnName.Length()-1] == '*')
            //                  printf("It's a pointer!\n");

            TClass *mclass = TClass::GetClass(returnName.Data());
            return mclass;
         }


         // Check if it's a proto: must have ()
         // This might not be too safe to use   :(
         char *parentesis_ptr = (char*)strrchr(memberName.Data(), '(');
         if (parentesis_ptr) *parentesis_ptr = 0;


         if (0) printf("Trying method [%s] of class [%s] ...\n",
            memberName.Data(), pclass->GetName());

         // Check if it's a method
         TMethod *mptr = 0; // pclass->GetMethodAny(memberName.Data());
         const TList  *mlist = pclass->GetListOfAllPublicMethods();
         next = mlist;
         while ((mptr = (TMethod *) next())) {
            if (strcmp(memberName.Data(),mptr->GetName())==0) break;
         }
         if (mptr)
         {
            TString returnName = mptr->GetReturnTypeName();

            if (0) printf("It's a method called [%s] with return type [%s]\n",
               memberName.Data(), returnName.Data());

            // This will handle the methods that returns a pointer to a class
            if (returnName[returnName.Length()-1] == '*')
            {
               returnName[returnName.Length()-1] = 0;
               fVarIsPointer = kTRUE;
            }
            else
            {
               fVarIsPointer = kFALSE;
            }

            TClass *mclass = TClass::GetClass(returnName.Data());
            return mclass;
         }
      }
   }

   //
   // frodo: End of Member and Proto recursion code
   //


   // not found...
   if (!varName_exists) {
      std::cerr << std::endl << "variable " << dblquote(varName) << " not defined."
         << std::endl;
      return 0;                 //* RETURN *//
   }

   /*****************************************************************************************/
   /*                                                                                       */
   /*  this section is really ugly.                                                         */
   /*  and slow.                                                                            */
   /*  it could be made a lot better if there was some way to tell whether or not a given   */
   /*  variable is a pointer or a pointer to a pointer.                                     */
   /*                                                                                       */
   /*****************************************************************************************/

   TString className = DetermineClass(varName);

   if (className.IsNull() || className == "*") {
      // this will happen if "varName" is a fundamental type (as opposed to class type).
      // or a pointer to a pointer.
      // or a function pointer.
      std::cerr << std::endl << "problem determining class of " << dblquote(varName)
         << std::endl;
      return 0;                 //* RETURN *//
   }

   fVarIsPointer = className[className.Length() - 1] == '*';

   // frodo: I shouldn't have to do this, but for some reason now I have to
   //        otherwise the varptr->[TAB] won't work    :(
   if (fVarIsPointer || className[className.Length() - 1] == '&')
      className[className.Length()-1] = 0;

   // frodo case '.' dealt with in the previous if statement?!
   // frodo: I wasn't able to put the automatic "." to "->" replacing working
   //        so I just commented out.
   //


   //   Bool_t varIsPointer = className[className.Length() - 1] == '*';

   //printf("Context is %i, fContext is %i, pointer is %i\n", context, fContext, fVarIsPointer);

   if (fVarIsPointer &&
      (context == kCXX_DirectMember || context == kCXX_DirectProto)) {
         // user is using operator.() instead of operator->()
         // ==>
         //      1. we are in wrong context.
         //      2. user is lazy
         //      3. or maybe confused

         // 1. fix the context
         switch (context) {
      case kCXX_DirectMember:
         context = kCXX_IndirectMember;
         break;
      case kCXX_DirectProto:
         context = kCXX_IndirectProto;
         break;
      default:
         Fatal("TTabCom::MakeClassFromVarName","Conext case %d not handled",context);
         return 0; // Avoid warning about uninitialized pClass.
         }

         // 2. fix the operator.
         int i;
         for (i = *fpLoc; fBuf[i] != '.'; i -= 1) {
         }
         int loc = i;
         for (i = strlen(fBuf); i >= loc; i -= 1) {
            fBuf[i + 1] = fBuf[i];
         }
         fBuf[loc] = '-';
         fBuf[loc + 1] = '>';
         *fpLoc += 1;

         // 3. inform the user.
         std::cerr << std::endl << dblquote(varName) <<
            " is of pointer type. Use this operator: ->" << std::endl;
         return 0;
   }

   if (context == kCXX_IndirectMember || context == kCXX_IndirectProto) {
      if (fVarIsPointer) {
         // frodo: This part no longer makes sense...
         className.Chop();      // remove the '*'

         if (className[className.Length() - 1] == '*') {
            std::cerr << std::endl << "can't handle pointers to pointers." << std::endl;
            return 0;           // RETURN
         }
      } else {
         // user is using operator->() instead of operator.()
         // ==>
         //      1. we are in wrong context.
         //      2. user is lazy
         //      3. or maybe confused

         // 1. fix the context
         switch (context) {
         case kCXX_IndirectMember:
            context = kCXX_DirectMember;
            break;
         case kCXX_IndirectProto:
            context = kCXX_DirectProto;
            break;
         default:
            Fatal("TTabCom::MakeClassFromVarName","Conext case %d not handled",context);
            return 0; // Avoid warning about uninitialized pClass.
         }

         // 2. fix the operator.
         int i;
         for (i = *fpLoc; fBuf[i - 1] != '-' && fBuf[i] != '>'; i -= 1) {
         }
         fBuf[i - 1] = '.';
         int len = strlen(fBuf);
         for (; i < len; i += 1) {
            fBuf[i] = fBuf[i + 1];
         }
         *fpLoc -= 1;

         // 3. inform the user.
         std::cerr << std::endl << dblquote(varName) <<
             " is not of pointer type. Use this operator: ." << std::endl;
         return 0;
      }
   }

   className.Strip(TString::kBoth);

   return TClass::GetClass(className);
}

////////////////////////////////////////////////////////////////////////////////
/// [private]

void TTabCom::SetPattern(EContext_t handle, const char regexp[])
{
   // prevent overflow
   if (handle >= kNUM_PAT) {
      std::cerr << std::endl
          << "ERROR: handle="
          << (int) handle << " >= kNUM_PAT=" << (int) kNUM_PAT << std::endl;
      return;
   }

   fRegExp[handle] = regexp;
   Makepat(regexp, fPat[handle], MAX_LEN_PAT);
}



////////////////////////////////////////////////////////////////////////////////
///
/// Returns the place in the string where to put the \0, starting the search
/// from "start"
///

int TTabCom::ParseReverse(const char *var_str, int start)
{
   int end = 0;
   if (start > (int)strlen(var_str)) start = strlen(var_str);

   for (int i = start; i > 0; i--)
   {
      if (var_str[i] == '.') return i;
      if (var_str[i] == '>' && i > 0 && var_str[i-1] == '-') return i-1;
   }

   return end;
}
