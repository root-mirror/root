// @(#)root/base:$Id$
// Author: Fons Rademakers   26/1/2002

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPluginManager
#define ROOT_TPluginManager


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPluginManager                                                       //
//                                                                      //
// This class implements a plugin library manager. It keeps track of    //
// a list of plugin handlers. A plugin handler knows which plugin       //
// library to load to get a specific class that is used to extend the   //
// functionality of a specific base class and how to create an object   //
// of this class. For example, to extend the base class TFile to be     //
// able to read SQLite files one needs to load the plugin library       //
// libRSQLite.so which defines the TSQLiteServer class. This loading    //
// should be triggered when a given URI contains a regular expression   //
// defined by the handler.                                              //
// Plugin handlers can be defined via macros in a list of plugin        //
// directories. With $ROOTSYS/etc/plugins the default top plugin        //
// directory specified in $ROOTSYS/etc/system.rootrc. Additional        //
// directories can be specified by adding them to the end of the list.  //
// Macros for identical plugin handlers in later directories will       //
// override previous ones (the inverse of normal search path behavior). //
// The macros must have names like <BaseClass>/PX0_<PluginClass>.C,     //
// e.g.:                                                                //
//    TSQLServer/P20_TMySQLServer.C, etc.                               //
// to allow easy sorting and grouping. If the BaseClass is in a         //
// namespace the directory must have the name NameSpace@@BaseClass as   //
// : is a reserved pathname character on some operating systems.        //
// Macros not beginning with 'P' and ending with ".C" are ignored.      //
// These macros typically look like:                                    //
//                                                                      //
//   void P10_TDCacheFile()                                             //
//   {                                                                  //
//       gPluginMgr->AddHandler("TFile", "^dcache", "TDCacheFile",      //
//          "DCache", "TDCacheFile(const char*,Option_t*)");            //
//   }                                                                  //
//                                                                      //
// Plugin handlers can also be defined via resources in the .rootrc     //
// file. Although now deprecated this method still works for backward   //
// compatibility, e.g.:                                                 //
//                                                                      //
//   Plugin.TSQLServer:  ^mysql:  TMySQLServer MySQL  "<constructor>"   //
//   +Plugin.TSQLServer: ^pgsql:  TPgSQLServer PgSQL  "<constructor>"   //
//   Plugin.TVirtualFitter: *     TFitter      Minuit "TFitter(Int_t)"  //
//                                                                      //
// Where the + in front of Plugin.TSQLServer says that it extends the   //
// existing definition of TSQLServer, useful when there is more than    //
// one plugin that can extend the same base class. The "<constructor>"  //
// should be the constructor or a static method that generates an       //
// instance of the specified class. Global methods should start with    //
// "::" in their name, like "::CreateFitter()".                         //
// Instead of being a shared library a plugin can also be a CINT        //
// script, so instead of libDialog.so one can have Dialog.C.            //
// The * is a placeholder in case there is no need for a URI to         //
// differentiate between different plugins for the same base class.     //
// For the default plugins see $ROOTSYS/etc/system.rootrc.              //
//                                                                      //
// Plugin handlers can also be registered at run time, e.g.:            //
//                                                                      //
//   gPluginMgr->AddHandler("TSQLServer", "^sqlite:",                   //
//                          "TSQLiteServer", "RSQLite",                 //
//             "TSQLiteServer(const char*,const char*,const char*)");   //
//                                                                      //
// A list of currently defined handlers can be printed using:           //
//                                                                      //
//   gPluginMgr->Print(); // use option="a" to see ctors                //
//                                                                      //
// The use of the plugin library manager removes all textual references //
// to hard-coded class and library names and the resulting dependencies //
// in the base classes. The plugin manager is used to extend a.o.       //
// TFile, TSQLServer, TGrid, etc. functionality.                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TString.h"
#include "TMethodCall.h"
#include "TVirtualMutex.h"
#include "TInterpreter.h"

class TEnv;
class TList;
class THashTable;
class TFunction;
class TPluginManager;

#include <atomic>

class TPluginHandler : public TObject {

friend class TPluginManager;

private:
   using AtomicInt_t = std::atomic<Int_t>;

   TString      fBase;      // base class which will be extended by plugin
   TString      fRegexp;    // regular expression which must be matched in URI
   TString      fClass;     // class to be loaded from plugin library
   TString      fPlugin;    // plugin library which should contain fClass
   TString      fCtor;      // ctor used to instantiate object of fClass
   TString      fOrigin;    // origin of plugin handler definition
   TMethodCall *fCallEnv;   //!ctor method call environment
   TFunction   *fMethod;    //!ctor method or global function
   AtomicInt_t  fCanCall;   //!if 1 fCallEnv is ok, -1 fCallEnv is not ok, 0 fCallEnv not setup yet.
   Bool_t       fIsMacro;   // plugin is a macro and not a library
   Bool_t       fIsGlobal;  // plugin ctor is a global function

   TPluginHandler() :
      fBase(), fRegexp(), fClass(), fPlugin(), fCtor(), fOrigin(),
      fCallEnv(nullptr), fMethod(nullptr), fCanCall(0), fIsMacro(kTRUE), fIsGlobal(kTRUE) { }
   TPluginHandler(const char *base, const char *regexp,
                  const char *className, const char *pluginName,
                  const char *ctor, const char *origin);
   TPluginHandler(const TPluginHandler &) = delete;
   TPluginHandler& operator=(const TPluginHandler &) = delete;

   ~TPluginHandler();

   const char *GetBase() const { return fBase; }
   const char *GetRegexp() const { return fRegexp; }
   const char *GetPlugin() const { return fPlugin; }
   const char *GetCtor() const { return fCtor; }
   const char *GetOrigin() const { return fOrigin; }

   Bool_t CanHandle(const char *base, const char *uri);
   void   SetupCallEnv();

   Bool_t CheckForExecPlugin(Int_t nargs);

public:
   const char *GetClass() const { return fClass; }
   Int_t       CheckPlugin() const;
   Int_t       LoadPlugin();

   template <typename... T> Long_t ExecPluginImpl(const T&... params)
   {
      auto nargs = sizeof...(params);
      if (!CheckForExecPlugin(nargs)) return 0;

      // The fCallEnv object is shared, since the PluginHandler is a global
      // resource ... and both SetParams and Execute ends up taking the lock
      // individually anyway ...

      R__LOCKGUARD(gInterpreterMutex);
      fCallEnv->SetParams(params...);

      Long_t ret;
      fCallEnv->Execute(ret);

      return ret;
   }

   template <typename... T> Long_t ExecPlugin(int nargs, const T&... params)
   {
      // For backward compatibility.
      if ((gDebug > 1) && (nargs != (int)sizeof...(params))) {
         Warning("ExecPlugin","Announced number of args different from the real number of argument passed %d vs %lu",
                 nargs, (unsigned long)sizeof...(params) );
      }
      return ExecPluginImpl(params...);
   }

   void        Print(Option_t *opt = "") const;

   ClassDef(TPluginHandler,3)  // Handler for plugin libraries
};


class TPluginManager : public TObject {

private:
   TList      *fHandlers;     // list of plugin handlers
   THashTable *fBasesLoaded;  //! table of base classes already checked or loaded
   Bool_t      fReadingDirs;  //! true if we are running LoadHandlersFromPluginDirs

   TPluginManager(const TPluginManager &) = delete;
   TPluginManager& operator=(const TPluginManager &) = delete;
   void   LoadHandlerMacros(const char *path);

public:
   TPluginManager() : fHandlers(0), fBasesLoaded(0), fReadingDirs(kFALSE) { }
   ~TPluginManager();

   void   LoadHandlersFromEnv(TEnv *env);
   void   LoadHandlersFromPluginDirs(const char *base = 0);
   void   AddHandler(const char *base, const char *regexp,
                     const char *className, const char *pluginName,
                     const char *ctor = 0, const char *origin = 0);
   void   RemoveHandler(const char *base, const char *regexp = 0);

   TPluginHandler *FindHandler(const char *base, const char *uri = 0);

   void   Print(Option_t *opt = "") const;
   Int_t  WritePluginMacros(const char *dir, const char *plugin = 0) const;
   Int_t  WritePluginRecords(const char *envFile, const char *plugin = 0) const;

   ClassDef(TPluginManager,1)  // Manager for plugin handlers
};

R__EXTERN TPluginManager *gPluginMgr;

#endif
