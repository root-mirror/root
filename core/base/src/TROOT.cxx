// @(#)root/base:$Id$
// Author: Rene Brun   08/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TROOT
\ingroup Base

ROOT top level object description.

The TROOT object is the entry point to the ROOT system.
The single instance of TROOT is accessible via the global gROOT.
Using the gROOT pointer one has access to basically every object
created in a ROOT based program. The TROOT object is essentially a
container of several lists pointing to the main ROOT objects.

The following lists are accessible from gROOT object:

~~~ {.cpp}
      gROOT->GetListOfClasses
      gROOT->GetListOfColors
      gROOT->GetListOfTypes
      gROOT->GetListOfGlobals
      gROOT->GetListOfGlobalFunctions
      gROOT->GetListOfFiles
      gROOT->GetListOfMappedFiles
      gROOT->GetListOfSockets
      gROOT->GetListOfSecContexts
      gROOT->GetListOfCanvases
      gROOT->GetListOfStyles
      gROOT->GetListOfFunctions
      gROOT->GetListOfSpecials (for example graphical cuts)
      gROOT->GetListOfGeometries
      gROOT->GetListOfBrowsers
      gROOT->GetListOfCleanups
      gROOT->GetListOfMessageHandlers
~~~

The TROOT class provides also many useful services:
  - Get pointer to an object in any of the lists above
  - Time utilities TROOT::Time

The ROOT object must be created as a static object. An example
of a main program creating an interactive version is shown below:

### Example of a main program

~~~ {.cpp}
      #include "TRint.h"

      int main(int argc, char **argv)
      {
         TRint *theApp = new TRint("ROOT example", &argc, argv);

         // Init Intrinsics, build all windows, and enter event loop
         theApp->Run();

         return(0);
      }
~~~
*/

#include <ROOT/RConfig.hxx>
#include <ROOT/TErrorDefaultHandler.hxx>
#include "RConfigure.h"
#include "RConfigOptions.h"
#include "RVersion.h"
#include "RGitCommit.h"
#include <string>
#include <map>
#include <cstdlib>
#ifdef WIN32
#include <io.h>
#include "Windows4Root.h"
#include <Psapi.h>
#define RTLD_DEFAULT ((void *)::GetModuleHandle(NULL))
//#define dlsym(library, function_name) ::GetProcAddress((HMODULE)library, function_name)
#define dlopen(library_name, flags) ::LoadLibrary(library_name)
#define dlclose(library) ::FreeLibrary((HMODULE)library)
char *dlerror() {
   static char Msg[1000];
   FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM, NULL, GetLastError(),
                 MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), Msg,
                 sizeof(Msg), NULL);
   return Msg;
}
FARPROC dlsym(void *library, const char *function_name)
{
   HMODULE hMods[1024];
   DWORD cbNeeded;
   FARPROC address = NULL;
   unsigned int i;
   if (library == RTLD_DEFAULT) {
      if (EnumProcessModules(::GetCurrentProcess(), hMods, sizeof(hMods), &cbNeeded)) {
         for (i = 0; i < (cbNeeded / sizeof(HMODULE)); i++) {
            address = ::GetProcAddress((HMODULE)hMods[i], function_name);
            if (address)
               return address;
         }
      }
      return address;
   } else {
      return ::GetProcAddress((HMODULE)library, function_name);
   }
}
#else
#include <dlfcn.h>
#endif

#include <iostream>
#include "ROOT/FoundationUtils.hxx"
#include "TROOT.h"
#include "TClass.h"
#include "TClassEdit.h"
#include "TClassGenerator.h"
#include "TDataType.h"
#include "TStyle.h"
#include "TObjectTable.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "THashList.h"
#include "TObjArray.h"
#include "TEnv.h"
#include "TError.h"
#include "TColor.h"
#include "TGlobal.h"
#include "TFunction.h"
#include "TVirtualPad.h"
#include "TBrowser.h"
#include "TSystemDirectory.h"
#include "TApplication.h"
#include "TInterpreter.h"
#include "TGuiFactory.h"
#include "TMessageHandler.h"
#include "TFolder.h"
#include "TQObject.h"
#include "TProcessUUID.h"
#include "TPluginManager.h"
#include "TVirtualMutex.h"
#include "TListOfTypes.h"
#include "TListOfDataMembers.h"
#include "TListOfEnumsWithLock.h"
#include "TListOfFunctions.h"
#include "TListOfFunctionTemplates.h"
#include "TFunctionTemplate.h"
#include "ThreadLocalStorage.h"
#include "TVirtualRWMutex.h"
#include "TVirtualX.h"

#if defined(R__UNIX)
#if defined(R__HAS_COCOA)
#include "TMacOSXSystem.h"
#include "TUrl.h"
#else
#include "TUnixSystem.h"
#endif
#elif defined(R__WIN32)
#include "TWinNTSystem.h"
#endif

extern "C" void R__SetZipMode(int);

static DestroyInterpreter_t *gDestroyInterpreter = nullptr;
static void *gInterpreterLib = nullptr;

// Mutex for protection of concurrent gROOT access
TVirtualMutex* gROOTMutex = nullptr;
ROOT::TVirtualRWMutex *ROOT::gCoreMutex = nullptr;

// For accessing TThread::Tsd indirectly.
void **(*gThreadTsd)(void*,Int_t) = 0;

//-------- Names of next three routines are a small homage to CMZ --------------
////////////////////////////////////////////////////////////////////////////////
/// Return version id as an integer, i.e. "2.22/04" -> 22204.

static Int_t IVERSQ()
{
   Int_t maj, min, cycle;
   sscanf(ROOT_RELEASE, "%d.%d/%d", &maj, &min, &cycle);
   return 10000*maj + 100*min + cycle;
}

////////////////////////////////////////////////////////////////////////////////
/// Return built date as integer, i.e. "Apr 28 2000" -> 20000428.

static Int_t IDATQQ(const char *date)
{
   static const char *months[] = {"Jan","Feb","Mar","Apr","May",
                                  "Jun","Jul","Aug","Sep","Oct",
                                  "Nov","Dec"};

   char  sm[12];
   Int_t yy, mm=0, dd;
   sscanf(date, "%s %d %d", sm, &dd, &yy);
   for (int i = 0; i < 12; i++)
      if (!strncmp(sm, months[i], 3)) {
         mm = i+1;
         break;
      }
   return 10000*yy + 100*mm + dd;
}

////////////////////////////////////////////////////////////////////////////////
/// Return built time as integer (with min precision), i.e.
/// "17:32:37" -> 1732.

static Int_t ITIMQQ(const char *time)
{
   Int_t hh, mm, ss;
   sscanf(time, "%d:%d:%d", &hh, &mm, &ss);
   return 100*hh + mm;
}

////////////////////////////////////////////////////////////////////////////////
/// Clean up at program termination before global objects go out of scope.

static void CleanUpROOTAtExit()
{
   if (gROOT) {
      R__LOCKGUARD(gROOTMutex);

      if (gROOT->GetListOfFiles())
         gROOT->GetListOfFiles()->Delete("slow");
      if (gROOT->GetListOfSockets())
         gROOT->GetListOfSockets()->Delete();
      if (gROOT->GetListOfMappedFiles())
         gROOT->GetListOfMappedFiles()->Delete("slow");
      if (gROOT->GetListOfClosedObjects())
         gROOT->GetListOfClosedObjects()->Delete("slow");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// A module and its headers. Intentionally not a copy:
/// If these strings end up in this struct they are
/// long lived by definition because they get passed in
/// before initialization of TCling.

namespace {
   struct ModuleHeaderInfo_t {
      ModuleHeaderInfo_t(const char* moduleName,
                         const char** headers,
                         const char** includePaths,
                         const char* payloadCode,
                         const char* fwdDeclCode,
                         void (*triggerFunc)(),
                         const TROOT::FwdDeclArgsToKeepCollection_t& fwdDeclsArgToSkip,
                         const char **classesHeaders,
                         bool hasCxxModule):
                           fModuleName(moduleName),
                           fHeaders(headers),
                           fPayloadCode(payloadCode),
                           fFwdDeclCode(fwdDeclCode),
                           fIncludePaths(includePaths),
                           fTriggerFunc(triggerFunc),
                           fClassesHeaders(classesHeaders),
                           fFwdNargsToKeepColl(fwdDeclsArgToSkip),
                           fHasCxxModule(hasCxxModule) {}

      const char* fModuleName; // module name
      const char** fHeaders; // 0-terminated array of header files
      const char* fPayloadCode; // Additional code to be given to cling at library load
      const char* fFwdDeclCode; // Additional code to let cling know about selected classes and functions
      const char** fIncludePaths; // 0-terminated array of header files
      void (*fTriggerFunc)(); // Pointer to the dict initialization used to find the library name
      const char** fClassesHeaders; // 0-terminated list of classes and related header files
      const TROOT::FwdDeclArgsToKeepCollection_t fFwdNargsToKeepColl; // Collection of
                                                                      // pairs of template fwd decls and number of
      bool fHasCxxModule; // Whether this module has a C++ module alongside it.
   };

   std::vector<ModuleHeaderInfo_t>& GetModuleHeaderInfoBuffer() {
      static std::vector<ModuleHeaderInfo_t> moduleHeaderInfoBuffer;
      return moduleHeaderInfoBuffer;
   }
}

Int_t  TROOT::fgDirLevel = 0;
Bool_t TROOT::fgRootInit = kFALSE;
Bool_t TROOT::fgMemCheck = kFALSE;

static void at_exit_of_TROOT() {
   if (ROOT::Internal::gROOTLocal)
      ROOT::Internal::gROOTLocal->~TROOT();
}

// This local static object initializes the ROOT system
namespace ROOT {
namespace Internal {
   class TROOTAllocator {
      // Simple wrapper to separate, time-wise, the call to the
      // TROOT destructor and the actual free-ing of the memory.
      //
      // Since the interpreter implementation (currently TCling) is
      // loaded via dlopen by libCore, the destruction of its global
      // variable (i.e. in particular clang's) is scheduled before
      // those in libCore so we need to schedule the call to the TROOT
      // destructor before that *but* we want to make sure the memory
      // stay around until libCore itself is unloaded so that code
      // using gROOT can 'properly' check for validity.
      //
      // The order of loading for is:
      //    libCore.so
      //    libRint.so
      //    ... anything other library hard linked to the executable ...
      //    ... for example libEvent
      //    libCling.so
      //    ... other libraries like libTree for example ....
      // and the destruction order is (of course) the reverse.
      // By default the unloading of the dictionary, does use
      // the service of the interpreter ... which of course
      // fails if libCling is already unloaded by that information
      // has not been registered per se.
      //
      // To solve this problem, we now schedule the destruction
      // of the TROOT object to happen _just_ before the
      // unloading/destruction of libCling so that we can
      // maximize the amount of clean-up we can do correctly
      // and we can still allocate the TROOT object's memory
      // statically.
      //
      union {
         TROOT fObj;
         char fHolder[sizeof(TROOT)];
      };
   public:
      TROOTAllocator(): fObj("root", "The ROOT of EVERYTHING")
      {}

      ~TROOTAllocator() {
         if (gROOTLocal) {
            gROOTLocal->~TROOT();
         }
      }
   };

   // The global gROOT is defined to be a function (ROOT::GetROOT())
   // which itself is dereferencing a function pointer.

   // Initially this function pointer's value is & GetROOT1 whose role is to
   // create and initialize the TROOT object itself.
   // At the very end of the TROOT constructor the value of the function pointer
   // is switch to & GetROOT2 whose role is to initialize the interpreter.

   // This mechanism was primarily intended to fix the issues with order in which
   // global TROOT and LLVM globals are initialized. TROOT was initializing
   // Cling, but Cling could not be used yet due to LLVM globals not being
   // Initialized yet.  The solution is to delay initializing the interpreter in
   // TROOT till after main() when all LLVM globals are initialized.

   // Technically, the mechanism used actually delay the interpreter
   // initialization until the first use of gROOT *after* the end of the
   // TROOT constructor.

   // So to delay until after the start of main, we also made sure that none
   // of the ROOT code (mostly the dictionary code) used during library loading
   // is using gROOT (directly or indirectly).

   // In practice, the initialization of the interpreter is now delayed until
   // the first use gROOT (or gInterpreter) after the start of main (but user
   // could easily break this by using gROOT in their library initialization
   // code).

   extern TROOT *gROOTLocal;

   TROOT *GetROOT1() {
      if (gROOTLocal)
         return gROOTLocal;
      static TROOTAllocator alloc;
      return gROOTLocal;
   }

   TROOT *GetROOT2() {
      static Bool_t initInterpreter = kFALSE;
      if (!initInterpreter) {
         initInterpreter = kTRUE;
         gROOTLocal->InitInterpreter();
         // Load and init threads library
         gROOTLocal->InitThreads();
      }
      return gROOTLocal;
   }
   typedef TROOT *(*GetROOTFun_t)();

   static GetROOTFun_t gGetROOT = &GetROOT1;

   static Func_t GetSymInLibImt(const char *funcname)
   {
      const static bool loadSuccess = dlsym(RTLD_DEFAULT, "usedToIdentifyRootClingByDlSym")? false : 0 <= gSystem->Load("libImt");
      if (loadSuccess) {
         if (auto sym = gSystem->DynFindSymbol(nullptr, funcname)) {
            return sym;
         } else {
            Error("GetSymInLibImt", "Cannot get symbol %s.", funcname);
         }
      }
      return nullptr;
   }

   //////////////////////////////////////////////////////////////////////////////
   /// Globally enables the parallel branch processing, which is a case of
   /// implicit multi-threading (IMT) in ROOT, activating the required locks.
   /// This IMT use case, implemented in TTree::GetEntry, spawns a task for
   /// each branch of the tree. Therefore, a task takes care of the reading,
   /// decompression and deserialisation of a given branch.
   void EnableParBranchProcessing()
   {
#ifdef R__USE_IMT
      static void (*sym)() = (void(*)())Internal::GetSymInLibImt("ROOT_TImplicitMT_EnableParBranchProcessing");
      if (sym)
         sym();
#else
      ::Warning("EnableParBranchProcessing", "Cannot enable parallel branch processing, please build ROOT with -Dimt=ON");
#endif
   }

   //////////////////////////////////////////////////////////////////////////////
   /// Globally disables the IMT use case of parallel branch processing,
   /// deactivating the corresponding locks.
   void DisableParBranchProcessing()
   {
#ifdef R__USE_IMT
      static void (*sym)() = (void(*)())Internal::GetSymInLibImt("ROOT_TImplicitMT_DisableParBranchProcessing");
      if (sym)
         sym();
#else
      ::Warning("DisableParBranchProcessing", "Cannot disable parallel branch processing, please build ROOT with -Dimt=ON");
#endif
   }

   //////////////////////////////////////////////////////////////////////////////
   /// Returns true if parallel branch processing is enabled.
   Bool_t IsParBranchProcessingEnabled()
   {
#ifdef R__USE_IMT
      static Bool_t (*sym)() = (Bool_t(*)())Internal::GetSymInLibImt("ROOT_TImplicitMT_IsParBranchProcessingEnabled");
      if (sym)
         return sym();
      else
         return kFALSE;
#else
      return kFALSE;
#endif
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Keeps track of the status of ImplicitMT w/o resorting to the load of
   /// libImt
   static Bool_t &IsImplicitMTEnabledImpl()
   {
      static Bool_t isImplicitMTEnabled = kFALSE;
      return isImplicitMTEnabled;
   }

} // end of Internal sub namespace
// back to ROOT namespace

   TROOT *GetROOT() {
      return (*Internal::gGetROOT)();
   }

   TString &GetMacroPath() {
      static TString macroPath;
      return macroPath;
   }

   // clang-format off
   ////////////////////////////////////////////////////////////////////////////////
   /// Enables the global mutex to make ROOT thread safe/aware.
   ///
   /// The following becomes safe:
   /// - concurrent construction and destruction of TObjects, including the ones registered in ROOT's global lists (e.g. gROOT->GetListOfCleanups(), gROOT->GetListOfFiles())
   /// - concurrent usage of _different_ ROOT objects from different threads, including ones with global state (e.g. TFile, TTree, TChain) with the exception of graphics classes (e.g. TCanvas)
   /// - concurrent calls to ROOT's type system classes, e.g. TClass and TEnum
   /// - concurrent calls to the interpreter through gInterpreter
   /// - concurrent loading of ROOT plug-ins
   ///
   /// In addition, gDirectory, gFile and gPad become a thread-local variable.
   /// In all threads, gDirectory defaults to gROOT, a singleton which supports thread-safe insertion and deletion of contents.
   /// gFile and gPad default to nullptr, as it is for single-thread programs.
   ///
   /// The ROOT graphics subsystem is not made thread-safe by this method. In particular drawing or printing different
   /// canvases from different threads (and analogous operations such as invoking `Draw` on a `TObject`) is not thread-safe.
   ///
   /// Note that there is no `DisableThreadSafety()`. ROOT's thread-safety features cannot be disabled once activated.
   // clang-format on
   void EnableThreadSafety()
   {
      static void (*sym)() = (void(*)())Internal::GetSymInLibImt("ROOT_TThread_Initialize");
      if (sym)
         sym();
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// @param[in] numthreads Number of threads to use. If not specified or
   ///                       set to zero, the number of threads is automatically
   ///                       decided by the implementation. Any other value is
   ///                       used as a hint.
   ///
   /// ROOT must be built with the compilation flag `imt=ON` for this feature to be available.
   /// The following objects and methods automatically take advantage of
   /// multi-threading if a call to `EnableImplicitMT` has been made before usage:
   ///
   ///  - RDataFrame internally runs the event-loop by parallelizing over clusters of entries
   ///  - TTree::GetEntry reads multiple branches in parallel
   ///  - TTree::FlushBaskets writes multiple baskets to disk in parallel
   ///  - TTreeCacheUnzip decompresses the baskets contained in a TTreeCache in parallel
   ///  - THx::Fit performs in parallel the evaluation of the objective function over the data
   ///  - TMVA::DNN trains the deep neural networks in parallel
   ///  - TMVA::BDT trains the classifier in parallel and multiclass BDTs are evaluated in parallel
   ///
   /// EnableImplicitMT calls in turn EnableThreadSafety.
   /// The 'numthreads' parameter allows to control the number of threads to
   /// be used by the implicit multi-threading. However, this parameter is just
   /// a hint for ROOT: it will try to satisfy the request if the execution
   /// scenario allows it. For example, if ROOT is configured to use an external
   /// scheduler, setting a value for 'numthreads' might not have any effect.
   void EnableImplicitMT(UInt_t numthreads)
   {
#ifdef R__USE_IMT
      if (ROOT::Internal::IsImplicitMTEnabledImpl())
         return;
      EnableThreadSafety();
      static void (*sym)(UInt_t) = (void(*)(UInt_t))Internal::GetSymInLibImt("ROOT_TImplicitMT_EnableImplicitMT");
      if (sym)
         sym(numthreads);
      ROOT::Internal::IsImplicitMTEnabledImpl() = true;
#else
      ::Warning("EnableImplicitMT", "Cannot enable implicit multi-threading with %d threads, please build ROOT with -Dimt=ON", numthreads);
#endif
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Disables the implicit multi-threading in ROOT (see EnableImplicitMT).
   void DisableImplicitMT()
   {
#ifdef R__USE_IMT
      static void (*sym)() = (void(*)())Internal::GetSymInLibImt("ROOT_TImplicitMT_DisableImplicitMT");
      if (sym)
         sym();
      ROOT::Internal::IsImplicitMTEnabledImpl() = kFALSE;
#else
      ::Warning("DisableImplicitMT", "Cannot disable implicit multi-threading, please build ROOT with -Dimt=ON");
#endif
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Returns true if the implicit multi-threading in ROOT is enabled.
   Bool_t IsImplicitMTEnabled()
   {
      return ROOT::Internal::IsImplicitMTEnabledImpl();
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Returns the size of ROOT's thread pool
   UInt_t GetThreadPoolSize()
   {
#ifdef R__USE_IMT
      static UInt_t (*sym)() = (UInt_t(*)())Internal::GetSymInLibImt("ROOT_MT_GetThreadPoolSize");
      if (sym)
         return sym();
      else
         return 0;
#else
      return 0;
#endif
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Returns the size of the pool used for implicit multi-threading.
   UInt_t GetImplicitMTPoolSize()
   {
      return GetThreadPoolSize();
   }
} // end of ROOT namespace

TROOT *ROOT::Internal::gROOTLocal = ROOT::GetROOT();

// Global debug flag (set to > 0 to get debug output).
// Can be set either via the interpreter (gDebug is exported to CINT),
// via the rootrc resource "Root.Debug", via the shell environment variable
// ROOTDEBUG, or via the debugger.
Int_t gDebug;


ClassImp(TROOT);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TROOT::TROOT() : TDirectory(),
     fLineIsProcessing(0), fVersion(0), fVersionInt(0), fVersionCode(0),
     fVersionDate(0), fVersionTime(0), fBuiltDate(0), fBuiltTime(0),
     fTimer(0), fApplication(0), fInterpreter(0), fBatch(kTRUE),
     fIsWebDisplay(kFALSE), fIsWebDisplayBatch(kFALSE), fEditHistograms(kTRUE),
     fFromPopUp(kTRUE),fMustClean(kTRUE),fForceStyle(kFALSE),
     fInterrupt(kFALSE),fEscape(kFALSE),fExecutingMacro(kFALSE),fEditorMode(0),
     fPrimitive(0),fSelectPad(0),fClasses(0),fTypes(0),fGlobals(0),fGlobalFunctions(0),
     fClosedObjects(0),fFiles(0),fMappedFiles(0),fSockets(0),fCanvases(0),fStyles(0),fFunctions(0),
     fTasks(0),fColors(0),fGeometries(0),fBrowsers(0),fSpecials(0),fCleanups(0),
     fMessageHandlers(0),fStreamerInfo(0),fClassGenerators(0),fSecContexts(0),
     fProofs(0),fClipboard(0),fDataSets(0),fUUIDs(0),fRootFolder(0),fBrowsables(0),
     fPluginManager(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the ROOT system. The creation of the TROOT object initializes
/// the ROOT system. It must be the first ROOT related action that is
/// performed by a program. The TROOT object must be created on the stack
/// (can not be called via new since "operator new" is protected). The
/// TROOT object is either created as a global object (outside the main()
/// program), or it is one of the first objects created in main().
/// Make sure that the TROOT object stays in scope for as long as ROOT
/// related actions are performed. TROOT is a so called singleton so
/// only one instance of it can be created. The single TROOT object can
/// always be accessed via the global pointer gROOT.
/// The name and title arguments can be used to identify the running
/// application. The initfunc argument can contain an array of
/// function pointers (last element must be 0). These functions are
/// executed at the end of the constructor. This way one can easily
/// extend the ROOT system without adding permanent dependencies
/// (e.g. the graphics system is initialized via such a function).

TROOT::TROOT(const char *name, const char *title, VoidFuncPtr_t *initfunc)
   : TDirectory(), fLineIsProcessing(0), fVersion(0), fVersionInt(0), fVersionCode(0),
     fVersionDate(0), fVersionTime(0), fBuiltDate(0), fBuiltTime(0),
     fTimer(0), fApplication(0), fInterpreter(0), fBatch(kTRUE),
     fIsWebDisplay(kFALSE), fIsWebDisplayBatch(kFALSE), fEditHistograms(kTRUE),
     fFromPopUp(kTRUE),fMustClean(kTRUE),fForceStyle(kFALSE),
     fInterrupt(kFALSE),fEscape(kFALSE),fExecutingMacro(kFALSE),fEditorMode(0),
     fPrimitive(0),fSelectPad(0),fClasses(0),fTypes(0),fGlobals(0),fGlobalFunctions(0),
     fClosedObjects(0),fFiles(0),fMappedFiles(0),fSockets(0),fCanvases(0),fStyles(0),fFunctions(0),
     fTasks(0),fColors(0),fGeometries(0),fBrowsers(0),fSpecials(0),fCleanups(0),
     fMessageHandlers(0),fStreamerInfo(0),fClassGenerators(0),fSecContexts(0),
     fProofs(0),fClipboard(0),fDataSets(0),fUUIDs(0),fRootFolder(0),fBrowsables(0),
     fPluginManager(0)
{
   if (fgRootInit || ROOT::Internal::gROOTLocal) {
      //Warning("TROOT", "only one instance of TROOT allowed");
      return;
   }

   R__LOCKGUARD(gROOTMutex);

   ROOT::Internal::gROOTLocal = this;
   gDirectory = 0;

   SetName(name);
   SetTitle(title);

   // will be used by global "operator delete" so make sure it is set
   // before anything is deleted
   fMappedFiles = 0;

   // create already here, but only initialize it after gEnv has been created
   gPluginMgr = fPluginManager = new TPluginManager;

   // Initialize Operating System interface
   InitSystem();

   // Initialize static directory functions
   GetRootSys();
   GetBinDir();
   GetLibDir();
   GetIncludeDir();
   GetEtcDir();
   GetDataDir();
   GetDocDir();
   GetMacroDir();
   GetTutorialDir();
   GetSourceDir();
   GetIconPath();
   GetTTFFontDir();

   gRootDir = GetRootSys().Data();

   TDirectory::BuildDirectory(nullptr, nullptr);

   // Initialize interface to CINT C++ interpreter
   fVersionInt      = 0;  // check in TROOT dtor in case TCling fails
   fClasses         = 0;  // might be checked via TCling ctor
   fEnums           = 0;

   fConfigOptions   = R__CONFIGUREOPTIONS;
   fConfigFeatures  = R__CONFIGUREFEATURES;
   fVersion         = ROOT_RELEASE;
   fVersionCode     = ROOT_VERSION_CODE;
   fVersionInt      = IVERSQ();
   fVersionDate     = IDATQQ(ROOT_RELEASE_DATE);
   fVersionTime     = ITIMQQ(ROOT_RELEASE_TIME);
   fBuiltDate       = IDATQQ(__DATE__);
   fBuiltTime       = ITIMQQ(__TIME__);

   ReadGitInfo();

   fClasses         = new THashTable(800,3); fClasses->UseRWLock();
   //fIdMap           = new IdMap_t;
   fStreamerInfo    = new TObjArray(100); fStreamerInfo->UseRWLock();
   fClassGenerators = new TList;

   // usedToIdentifyRootClingByDlSym is available when TROOT is part of
   // rootcling.
   if (!dlsym(RTLD_DEFAULT, "usedToIdentifyRootClingByDlSym")) {
      // initialize plugin manager early
      fPluginManager->LoadHandlersFromEnv(gEnv);
   }

   TSystemDirectory *workdir = new TSystemDirectory("workdir", gSystem->WorkingDirectory());

   auto setNameLocked = [](TSeqCollection *l, const char *collection_name) {
      l->SetName(collection_name);
      l->UseRWLock();
      return l;
   };

   fTimer       = 0;
   fApplication = 0;
   fColors      = setNameLocked(new TObjArray(1000), "ListOfColors");
   fTypes       = 0;
   fGlobals     = 0;
   fGlobalFunctions = 0;
   // fList was created in TDirectory::Build but with different sizing.
   delete fList;
   fList        = new THashList(1000,3); fList->UseRWLock();
   fClosedObjects = setNameLocked(new TList, "ClosedFiles");
   fFiles       = setNameLocked(new TList, "Files");
   fMappedFiles = setNameLocked(new TList, "MappedFiles");
   fSockets     = setNameLocked(new TList, "Sockets");
   fCanvases    = setNameLocked(new TList, "Canvases");
   fStyles      = setNameLocked(new TList, "Styles");
   fFunctions   = setNameLocked(new TList, "Functions");
   fTasks       = setNameLocked(new TList, "Tasks");
   fGeometries  = setNameLocked(new TList, "Geometries");
   fBrowsers    = setNameLocked(new TList, "Browsers");
   fSpecials    = setNameLocked(new TList, "Specials");
   fBrowsables  = (TList*)setNameLocked(new TList, "Browsables");
   fCleanups    = setNameLocked(new THashList, "Cleanups");
   fMessageHandlers = setNameLocked(new TList, "MessageHandlers");
   fSecContexts = setNameLocked(new TList, "SecContexts");
   fProofs      = setNameLocked(new TList, "Proofs");
   fClipboard   = setNameLocked(new TList, "Clipboard");
   fDataSets    = setNameLocked(new TList, "DataSets");
   fTypes       = new TListOfTypes; fTypes->UseRWLock();

   TProcessID::AddProcessID();
   fUUIDs = new TProcessUUID();

   fRootFolder = new TFolder();
   fRootFolder->SetName("root");
   fRootFolder->SetTitle("root of all folders");
   fRootFolder->AddFolder("Classes",   "List of Active Classes",fClasses);
   fRootFolder->AddFolder("Colors",    "List of Active Colors",fColors);
   fRootFolder->AddFolder("MapFiles",  "List of MapFiles",fMappedFiles);
   fRootFolder->AddFolder("Sockets",   "List of Socket Connections",fSockets);
   fRootFolder->AddFolder("Canvases",  "List of Canvases",fCanvases);
   fRootFolder->AddFolder("Styles",    "List of Styles",fStyles);
   fRootFolder->AddFolder("Functions", "List of Functions",fFunctions);
   fRootFolder->AddFolder("Tasks",     "List of Tasks",fTasks);
   fRootFolder->AddFolder("Geometries","List of Geometries",fGeometries);
   fRootFolder->AddFolder("Browsers",  "List of Browsers",fBrowsers);
   fRootFolder->AddFolder("Specials",  "List of Special Objects",fSpecials);
   fRootFolder->AddFolder("Handlers",  "List of Message Handlers",fMessageHandlers);
   fRootFolder->AddFolder("Cleanups",  "List of RecursiveRemove Collections",fCleanups);
   fRootFolder->AddFolder("StreamerInfo","List of Active StreamerInfo Classes",fStreamerInfo);
   fRootFolder->AddFolder("SecContexts","List of Security Contexts",fSecContexts);
   fRootFolder->AddFolder("PROOF Sessions", "List of PROOF sessions",fProofs);
   fRootFolder->AddFolder("ROOT Memory","List of Objects in the gROOT Directory",fList);
   fRootFolder->AddFolder("ROOT Files","List of Connected ROOT Files",fFiles);

   // by default, add the list of files, tasks, canvases and browsers in the Cleanups list
   fCleanups->Add(fCanvases); fCanvases->SetBit(kMustCleanup);
   fCleanups->Add(fBrowsers); fBrowsers->SetBit(kMustCleanup);
   fCleanups->Add(fTasks);    fTasks->SetBit(kMustCleanup);
   fCleanups->Add(fFiles);    fFiles->SetBit(kMustCleanup);
   fCleanups->Add(fClosedObjects); fClosedObjects->SetBit(kMustCleanup);
   // And add TROOT's TDirectory personality
   fCleanups->Add(fList);

   fExecutingMacro= kFALSE;
   fForceStyle    = kFALSE;
   fFromPopUp     = kFALSE;
   fInterrupt     = kFALSE;
   fEscape        = kFALSE;
   fMustClean     = kTRUE;
   fPrimitive     = 0;
   fSelectPad     = 0;
   fEditorMode    = 0;
   fDefCanvasName = "c1";
   fEditHistograms= kFALSE;
   fLineIsProcessing = 1;   // This prevents WIN32 "Windows" thread to pick ROOT objects with mouse
   gDirectory     = this;
   gPad           = 0;

   //set name of graphical cut class for the graphics editor
   //cannot call SetCutClassName at this point because the TClass of TCutG
   //is not yet build
   fCutClassName = "TCutG";

   // Create a default MessageHandler
   new TMessageHandler((TClass*)0);

   // Create some styles
   gStyle = 0;
   TStyle::BuildStyles();
   SetStyle(gEnv->GetValue("Canvas.Style", "Modern"));

   // Setup default (batch) graphics and GUI environment
   gBatchGuiFactory = new TGuiFactory;
   gGuiFactory      = gBatchGuiFactory;
   gGXBatch         = new TVirtualX("Batch", "ROOT Interface to batch graphics");
   gVirtualX        = gGXBatch;

#if defined(R__WIN32)
   fBatch = kFALSE;
#elif defined(R__HAS_COCOA)
   fBatch = kFALSE;
#else
   if (gSystem->Getenv("DISPLAY"))
      fBatch = kFALSE;
   else
      fBatch = kTRUE;
#endif

   int i = 0;
   while (initfunc && initfunc[i]) {
      (initfunc[i])();
      fBatch = kFALSE;  // put system in graphics mode (backward compatible)
      i++;
   }

   // Set initial/default list of browsable objects
   fBrowsables->Add(fRootFolder, "root");
   fBrowsables->Add(fProofs, "PROOF Sessions");
   fBrowsables->Add(workdir, gSystem->WorkingDirectory());
   fBrowsables->Add(fFiles, "ROOT Files");

   atexit(CleanUpROOTAtExit);

   ROOT::Internal::gGetROOT = &ROOT::Internal::GetROOT2;
}

////////////////////////////////////////////////////////////////////////////////
/// Clean up and free resources used by ROOT (files, network sockets,
/// shared memory segments, etc.).

TROOT::~TROOT()
{
   using namespace ROOT::Internal;

   if (gROOTLocal == this) {

      // If the interpreter has not yet been initialized, don't bother
      gGetROOT = &GetROOT1;

      // Mark the object as invalid, so that we can veto some actions
      // (like autoloading) while we are in the destructor.
      SetBit(TObject::kInvalidObject);

      // Turn-off the global mutex to avoid recreating mutexes that have
      // already been deleted during the destruction phase
      gGlobalMutex = 0;

      // Return when error occurred in TCling, i.e. when setup file(s) are
      // out of date
      if (!fVersionInt) return;

      // ATTENTION!!! Order is important!

      SafeDelete(fBrowsables);

      // FIXME: Causes rootcling to deadlock, debug and uncomment
      // SafeDelete(fRootFolder);

#ifdef R__COMPLETE_MEM_TERMINATION
      fSpecials->Delete();   SafeDelete(fSpecials);    // delete special objects : PostScript, Minuit, Html
#endif

      fClosedObjects->Delete("slow"); // and closed files
      fFiles->Delete("slow");       // and files
      SafeDelete(fFiles);
      fSecContexts->Delete("slow"); SafeDelete(fSecContexts); // and security contexts
      fSockets->Delete();           SafeDelete(fSockets);     // and sockets
      fMappedFiles->Delete("slow");                     // and mapped files
      TSeqCollection *tl = fMappedFiles; fMappedFiles = 0; delete tl;

      SafeDelete(fClosedObjects);

      delete fUUIDs;
      TProcessID::Cleanup();                            // and list of ProcessIDs

      fFunctions->Delete();  SafeDelete(fFunctions);   // etc..
      fGeometries->Delete(); SafeDelete(fGeometries);
      fBrowsers->Delete();   SafeDelete(fBrowsers);
      SafeDelete(fCanvases);
      fColors->Delete();     SafeDelete(fColors);
      fStyles->Delete();     SafeDelete(fStyles);

#ifdef R__COMPLETE_MEM_TERMINATION
      if (gGuiFactory != gBatchGuiFactory) SafeDelete(gGuiFactory);
      SafeDelete(gBatchGuiFactory);
      if (gGXBatch != gVirtualX) SafeDelete(gGXBatch);
      SafeDelete(gVirtualX);
#endif

      // Stop emitting signals
      TQObject::BlockAllSignals(kTRUE);

      fMessageHandlers->Delete(); SafeDelete(fMessageHandlers);

#ifdef R__COMPLETE_MEM_TERMINATION
      SafeDelete(fCanvases);
      SafeDelete(fTasks);
      SafeDelete(fProofs);
      SafeDelete(fDataSets);
      SafeDelete(fClipboard);

      fCleanups->Clear();
      delete fPluginManager; gPluginMgr = fPluginManager = 0;
      delete gClassTable;  gClassTable = 0;
      delete gEnv; gEnv = 0;

      if (fTypes) fTypes->Delete();
      SafeDelete(fTypes);
      if (fGlobals) fGlobals->Delete();
      SafeDelete(fGlobals);
      if (fGlobalFunctions) fGlobalFunctions->Delete();
      SafeDelete(fGlobalFunctions);
      fEnums.load()->Delete();

      // FIXME: Causes segfault in rootcling, debug and uncomment
      // fClasses->Delete();    SafeDelete(fClasses);     // TClass'es must be deleted last
#endif

      // Remove shared libraries produced by the TSystem::CompileMacro() call
      gSystem->CleanCompiledMacros();

      // Cleanup system class
      ROOT::Internal::SetErrorSystemMsgHandler(ROOT::Internal::ErrorSystemMsgHandlerFunc_t());
      SetErrorHandler(ROOT::Internal::MinimalErrorHandler);
      ROOT::Internal::ReleaseDefaultErrorHandler();
      delete gSystem;

      // ROOT-6022:
      //   if (gInterpreterLib) dlclose(gInterpreterLib);
#ifdef R__COMPLETE_MEM_TERMINATION
      // On some 'newer' platform (Fedora Core 17+, Ubuntu 12), the
      // initialization order is (by default?) is 'wrong' and so we can't
      // delete the interpreter now .. because any of the static in the
      // interpreter's library have already been deleted.
      // On the link line, we must list the most dependent .o file
      // and end with the least dependent (LLVM libraries), unfortunately,
      // Fedora Core 17+ or Ubuntu 12 will also execute the initialization
      // in the same order (hence doing libCore's before LLVM's and
      // vice et versa for both the destructor.  We worked around the
      // initialization order by delay the TROOT creation until first use.
      // We can not do the same for destruction as we have no way of knowing
      // the last access ...
      // So for now, let's avoid delete TCling except in the special build
      // checking the completeness of the termination deletion.

      // TODO: Should we do more cleanup here than just call delete?
      // Segfaults rootcling in some cases, debug and uncomment:
      //
      //    delete fInterpreter;

      // We cannot delete fCleanups because of the logic in atexit which needs it.
      SafeDelete(fCleanups);
#endif

#ifdef _MSC_VER
      // usedToIdentifyRootClingByDlSym is available when TROOT is part of rootcling.
      if (dlsym(RTLD_DEFAULT, "usedToIdentifyRootClingByDlSym")) {
         // deleting the interpreter makes things crash at exit in some cases
         delete fInterpreter;
      }
#else
      // deleting the interpreter makes things crash at exit in some cases
      delete fInterpreter;
#endif

      // Prints memory stats
      TStorage::PrintStatistics();

      gROOTLocal = 0;
      fgRootInit = kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add a class to the list and map of classes.
/// This routine is deprecated, use TClass::AddClass directly.

void TROOT::AddClass(TClass *cl)
{
   TClass::AddClass(cl);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a class generator.  This generator will be called by TClass::GetClass
/// in case its does not find a loaded rootcint dictionary to request the
/// creation of a TClass object.

void TROOT::AddClassGenerator(TClassGenerator *generator)
{
   if (!generator) return;
   fClassGenerators->Add(generator);
}

////////////////////////////////////////////////////////////////////////////////
/// Append object to this directory.
///
/// If replace is true:
///   remove any existing objects with the same same (if the name is not "")

void TROOT::Append(TObject *obj, Bool_t replace /* = kFALSE */)
{
   R__LOCKGUARD(gROOTMutex);
   TDirectory::Append(obj,replace);
}

////////////////////////////////////////////////////////////////////////////////
/// Add browsable objects to TBrowser.

void TROOT::Browse(TBrowser *b)
{
   TObject *obj;
   TIter next(fBrowsables);

   while ((obj = (TObject *) next())) {
      const char *opt = next.GetOption();
      if (opt && strlen(opt))
         b->Add(obj, opt);
      else
         b->Add(obj, obj->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// return class status bit kClassSaved for class cl
/// This function is called by the SavePrimitive functions writing
/// the C++ code for an object.

Bool_t TROOT::ClassSaved(TClass *cl)
{
   if (cl == 0) return kFALSE;
   if (cl->TestBit(TClass::kClassSaved)) return kTRUE;
   cl->SetBit(TClass::kClassSaved);
   return kFALSE;
}

namespace {
   static void R__ListSlowClose(TList *files)
   {
      // Routine to close a list of files using the 'slow' techniques
      // that also for the deletion ot update the list itself.

      static TObject harmless;
      TObjLink *cursor = files->FirstLink();
      while (cursor) {
         TDirectory *dir = static_cast<TDirectory*>( cursor->GetObject() );
         if (dir) {
            // In order for the iterator to stay valid, we must
            // prevent the removal of the object (dir) from the list
            // (which is done in TFile::Close).   We can also can not
            // just move to the next iterator since the Close might
            // also (indirectly) remove that file.
            // So we SetObject to a harmless value, so that 'dir'
            // is not seen as part of the list.
            // We will later, remove all the object (see files->Clear()
            cursor->SetObject(&harmless); // this must not be zero otherwise things go wrong.
            // See related comment at the files->Clear("nodelete");
            dir->Close("nodelete");
            // Put it back
            cursor->SetObject(dir);
         }
         cursor = cursor->Next();
      };
      // Now were done, clear the list but do not delete the objects as
      // they have been moved to the list of closed objects and must be
      // deleted from there in order to avoid a double delete from a
      // use objects (on the interpreter stack).
      files->Clear("nodelete");
   }

   static void R__ListSlowDeleteContent(TList *files)
   {
      // Routine to delete the content of list of files using the 'slow' techniques

      static TObject harmless;
      TObjLink *cursor = files->FirstLink();
      while (cursor) {
         TDirectory *dir = dynamic_cast<TDirectory*>( cursor->GetObject() );
         if (dir) {
            // In order for the iterator to stay valid, we must
            // prevent the removal of the object (dir) from the list
            // (which is done in TFile::Close).   We can also can not
            // just move to the next iterator since the Close might
            // also (indirectly) remove that file.
            // So we SetObject to a harmless value, so that 'dir'
            // is not seen as part of the list.
            // We will later, remove all the object (see files->Clear()
            cursor->SetObject(&harmless); // this must not be zero otherwise things go wrong.
            // See related comment at the files->Clear("nodelete");
            dir->GetList()->Delete("slow");
            // Put it back
            cursor->SetObject(dir);
         }
         cursor = cursor->Next();
      };
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Close any files and sockets that gROOT knows about.
/// This can be used to insures that the files and sockets are closed before any library is unloaded!

void TROOT::CloseFiles()
{
   if (fFiles && fFiles->First()) {
      R__ListSlowClose(static_cast<TList*>(fFiles));
   }
   // and Close TROOT itself.
   Close("slow");
   // Now sockets.
   if (fSockets && fSockets->First()) {
      if (0==fCleanups->FindObject(fSockets) ) {
         fCleanups->Add(fSockets);
         fSockets->SetBit(kMustCleanup);
      }
      CallFunc_t *socketCloser = gInterpreter->CallFunc_Factory();
      Long_t offset = 0;
      TClass *socketClass = TClass::GetClass("TSocket");
      gInterpreter->CallFunc_SetFuncProto(socketCloser, socketClass->GetClassInfo(), "Close", "", &offset);
      if (gInterpreter->CallFunc_IsValid(socketCloser)) {
         static TObject harmless;
         TObjLink *cursor = static_cast<TList*>(fSockets)->FirstLink();
         TList notclosed;
         while (cursor) {
            TObject *socket = cursor->GetObject();
            // In order for the iterator to stay valid, we must
            // prevent the removal of the object (dir) from the list
            // (which is done in TFile::Close).   We can also can not
            // just move to the next iterator since the Close might
            // also (indirectly) remove that file.
            // So we SetObject to a harmless value, so that 'dir'
            // is not seen as part of the list.
            // We will later, remove all the object (see files->Clear()
            cursor->SetObject(&harmless); // this must not be zero otherwise things go wrong.

            if (socket->IsA()->InheritsFrom(socketClass)) {
               gInterpreter->CallFunc_Exec(socketCloser, ((char*)socket)+offset);
               // Put the object in the closed list for later deletion.
               socket->SetBit(kMustCleanup);
               fClosedObjects->AddLast(socket);
            } else {
               // Crap ... this is not a socket, likely Proof or something, let's try to find a Close
               Long_t other_offset;
               CallFunc_t *otherCloser = gInterpreter->CallFunc_Factory();
               gInterpreter->CallFunc_SetFuncProto(otherCloser, socket->IsA()->GetClassInfo(), "Close", "", &other_offset);
               if (gInterpreter->CallFunc_IsValid(otherCloser)) {
                  gInterpreter->CallFunc_Exec(otherCloser, ((char*)socket)+other_offset);
                  // Put the object in the closed list for later deletion.
                  socket->SetBit(kMustCleanup);
                  fClosedObjects->AddLast(socket);
               } else {
                  notclosed.AddLast(socket);
               }
               gInterpreter->CallFunc_Delete(otherCloser);
               // Put it back
               cursor->SetObject(socket);
            }
            cursor = cursor->Next();
         }
         // Now were done, clear the list
         fSockets->Clear();
         // Read the one we did not close
         cursor = notclosed.FirstLink();
         while (cursor) {
            static_cast<TList*>(fSockets)->AddLast(cursor->GetObject());
            cursor = cursor->Next();
         }
      }
      gInterpreter->CallFunc_Delete(socketCloser);
   }
   if (fMappedFiles && fMappedFiles->First()) {
      R__ListSlowClose(static_cast<TList*>(fMappedFiles));
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Execute the cleanups necessary at the end of the process, in particular
/// those that must be executed before the library start being unloaded.

void TROOT::EndOfProcessCleanups()
{
   // This will not delete the objects 'held' by the TFiles so that
   // they can still be 'reacheable' when ResetGlobals is run.
   CloseFiles();

   if (gInterpreter) {
      gInterpreter->ResetGlobals();
   }

   // Now delete the objects 'held' by the TFiles so that it
   // is done before the tear down of the libraries.
   if (fClosedObjects && fClosedObjects->First()) {
      R__ListSlowDeleteContent(static_cast<TList*>(fClosedObjects));
   }

   // Now a set of simpler things to delete.  See the same ordering in
   // TROOT::~TROOT
   fFunctions->Delete();
   fGeometries->Delete();
   fBrowsers->Delete();
   fCanvases->Delete("slow");
   fColors->Delete();
   fStyles->Delete();

   TQObject::BlockAllSignals(kTRUE);

   if (gInterpreter) {
      gInterpreter->ShutDown();
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Find an object in one Root folder

TObject *TROOT::FindObject(const TObject *) const
{
   Error("FindObject","Not yet implemented");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns address of a ROOT object if it exists
///
/// If name contains at least one "/" the function calls FindObjectany
/// else
/// This function looks in the following order in the ROOT lists:
///     - List of files
///     - List of memory mapped files
///     - List of functions
///     - List of geometries
///     - List of canvases
///     - List of styles
///     - List of specials
///     - List of materials in current geometry
///     - List of shapes in current geometry
///     - List of matrices in current geometry
///     - List of Nodes in current geometry
///     - Current Directory in memory
///     - Current Directory on file

TObject *TROOT::FindObject(const char *name) const
{
   if (name && strstr(name,"/")) return FindObjectAny(name);

   TObject *temp = 0;

   temp   = fFiles->FindObject(name);       if (temp) return temp;
   temp   = fMappedFiles->FindObject(name); if (temp) return temp;
   {
      R__LOCKGUARD(gROOTMutex);
      temp   = fFunctions->FindObject(name);if (temp) return temp;
   }
   temp   = fGeometries->FindObject(name);  if (temp) return temp;
   temp   = fCanvases->FindObject(name);    if (temp) return temp;
   temp   = fStyles->FindObject(name);      if (temp) return temp;
   {
      R__LOCKGUARD(gROOTMutex);
      temp = fSpecials->FindObject(name);   if (temp) return temp;
   }
   TIter next(fGeometries);
   TObject *obj;
   while ((obj=next())) {
      temp = obj->FindObject(name);         if (temp) return temp;
   }
   if (gDirectory) temp = gDirectory->Get(name);
   if (temp) return temp;
   if (gPad) {
      TVirtualPad *canvas = gPad->GetVirtCanvas();
      if (fCanvases->FindObject(canvas)) {  //this check in case call from TCanvas ctor
         temp = canvas->FindObject(name);
         if (!temp && canvas != gPad) temp  = gPad->FindObject(name);
      }
   }
   return temp;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns address and folder of a ROOT object if it exists
///
/// This function looks in the following order in the ROOT lists:
///     - List of files
///     - List of memory mapped files
///     - List of functions
///     - List of geometries
///     - List of canvases
///     - List of styles
///     - List of specials
///     - List of materials in current geometry
///     - List of shapes in current geometry
///     - List of matrices in current geometry
///     - List of Nodes in current geometry
///     - Current Directory in memory
///     - Current Directory on file

TObject *TROOT::FindSpecialObject(const char *name, void *&where)
{
   TObject *temp = 0;
   where = 0;

   if (!temp) {
      temp  = fFiles->FindObject(name);
      where = fFiles;
   }
   if (!temp) {
      temp  = fMappedFiles->FindObject(name);
      where = fMappedFiles;
   }
   if (!temp) {
      R__LOCKGUARD(gROOTMutex);
      temp  = fFunctions->FindObject(name);
      where = fFunctions;
   }
   if (!temp) {
      temp  = fCanvases->FindObject(name);
      where = fCanvases;
   }
   if (!temp) {
      temp  = fStyles->FindObject(name);
      where = fStyles;
   }
   if (!temp) {
      temp  = fSpecials->FindObject(name);
      where = fSpecials;
   }
   if (!temp) {
      TObject *glast = fGeometries->Last();
      if (glast) {where = glast; temp = glast->FindObject(name);}
   }
   if (!temp && gDirectory) {
      temp  = gDirectory->Get(name);
      where = gDirectory;
   }
   if (!temp && gPad) {
      TVirtualPad *canvas = gPad->GetVirtCanvas();
      if (fCanvases->FindObject(canvas)) {  //this check in case call from TCanvas ctor
         temp  = canvas->FindObject(name);
         where = canvas;
         if (!temp && canvas != gPad) {
            temp  = gPad->FindObject(name);
            where = gPad;
         }
      }
   }
   if (!temp) return 0;
   if (temp->TestBit(kNotDeleted)) return temp;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the first object with name starting at //root.
/// This function scans the list of all folders.
/// if no object found in folders, it scans the memory list of all files.

TObject *TROOT::FindObjectAny(const char *name) const
{
   TObject *obj = fRootFolder->FindObjectAny(name);
   if (obj) return obj;
   return gDirectory->FindObjectAnyFile(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Scan the memory lists of all files for an object with name

TObject *TROOT::FindObjectAnyFile(const char *name) const
{
   R__LOCKGUARD(gROOTMutex);
   TDirectory *d;
   TIter next(GetListOfFiles());
   while ((d = (TDirectory*)next())) {
      // Call explicitly TDirectory::FindObject to restrict the search to the
      // already in memory object.
      TObject *obj = d->TDirectory::FindObject(name);
      if (obj) return obj;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns class name of a ROOT object including CINT globals.

const char *TROOT::FindObjectClassName(const char *name) const
{
   // Search first in the list of "standard" objects
   TObject *obj = FindObject(name);
   if (obj) return obj->ClassName();

   // Is it a global variable?
   TGlobal *g = GetGlobal(name);
   if (g) return g->GetTypeName();

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return path name of obj somewhere in the //root/... path.
/// The function returns the first occurence of the object in the list
/// of folders. The returned string points to a static char array in TROOT.
/// If this function is called in a loop or recursively, it is the
/// user's responsibility to copy this string in their area.

const char *TROOT::FindObjectPathName(const TObject *) const
{
   Error("FindObjectPathName","Not yet implemented");
   return "??";
}

////////////////////////////////////////////////////////////////////////////////
/// return a TClass object corresponding to 'name' assuming it is an STL container.
/// In particular we looking for possible alternative name (default template
/// parameter, typedefs template arguments, typedefed name).

TClass *TROOT::FindSTLClass(const char *name, Bool_t load, Bool_t silent) const
{
   // Example of inputs are
   //   vector<int>  (*)
   //   vector<Int_t>
   //   vector<long long>
   //   vector<Long_64_t> (*)
   //   vector<int, allocator<int> >
   //   vector<Int_t, allocator<int> >
   //
   //   One of the possibly expensive operation is the resolving of the typedef
   //   which can provoke the parsing of the header files (and/or the loading
   //   of clang pcms information).

   R__LOCKGUARD(gInterpreterMutex);

   // Remove std::, allocator, typedef, add Long64_t, etc. in just one call.
   std::string normalized;
   TClassEdit::GetNormalizedName(normalized, name);

   TClass *cl = 0;
   if (normalized != name) cl = TClass::GetClass(normalized.c_str(),load,silent);

   if (load && cl==0) {
      // Create an Emulated class for this container.
      cl = gInterpreter->GenerateTClass(normalized.c_str(), kTRUE, silent);
   }

   return cl;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to class with name. Obsolete, use TClass::GetClass directly

TClass *TROOT::GetClass(const char *name, Bool_t load, Bool_t silent) const
{
   return TClass::GetClass(name,load,silent);
}


////////////////////////////////////////////////////////////////////////////////
/// Return pointer to class from its name. Obsolete, use TClass::GetClass directly
/// See TClass::GetClass

TClass *TROOT::GetClass(const std::type_info& typeinfo, Bool_t load, Bool_t silent) const
{
   return TClass::GetClass(typeinfo,load,silent);
}

////////////////////////////////////////////////////////////////////////////////
/// Return address of color with index color.

TColor *TROOT::GetColor(Int_t color) const
{
   TColor::InitializeColors();
   TObjArray *lcolors = (TObjArray*) GetListOfColors();
   if (!lcolors) return 0;
   if (color < 0 || color >= lcolors->GetSize()) return 0;
   TColor *col = (TColor*)lcolors->At(color);
   if (col && col->GetNumber() == color) return col;
   TIter   next(lcolors);
   while ((col = (TColor *) next()))
      if (col->GetNumber() == color) return col;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a default canvas.

TCanvas *TROOT::MakeDefCanvas() const
{
   return (TCanvas*)gROOT->ProcessLine("TCanvas::MakeDefCanvas();");
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to type with name.

TDataType *TROOT::GetType(const char *name, Bool_t /* load */) const
{
   return (TDataType*)gROOT->GetListOfTypes()->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to file with name.

TFile *TROOT::GetFile(const char *name) const
{
   R__LOCKGUARD(gROOTMutex);
   return (TFile*)GetListOfFiles()->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to style with name

TStyle *TROOT::GetStyle(const char *name) const
{
   return (TStyle*)GetListOfStyles()->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to function with name.

TObject *TROOT::GetFunction(const char *name) const
{
   if (name == 0 || name[0] == 0) {
      return 0;
   }

   {
      R__LOCKGUARD(gROOTMutex);
      TObject *f1 = fFunctions->FindObject(name);
      if (f1) return f1;
   }

   gROOT->ProcessLine("TF1::InitStandardFunctions();");

   R__LOCKGUARD(gROOTMutex);
   return fFunctions->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////

TFunctionTemplate *TROOT::GetFunctionTemplate(const char *name)
{
   if (!gInterpreter) return 0;

   if (!fFuncTemplate) fFuncTemplate = new TListOfFunctionTemplates(0);

   return (TFunctionTemplate*)fFuncTemplate->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to global variable by name. If load is true force
/// reading of all currently defined globals from CINT (more expensive).

TGlobal *TROOT::GetGlobal(const char *name, Bool_t load) const
{
   return (TGlobal *)gROOT->GetListOfGlobals(load)->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to global variable with address addr.

TGlobal *TROOT::GetGlobal(const TObject *addr, Bool_t /* load */) const
{
   if (addr == 0 || ((Long_t)addr) == -1) return 0;

   TInterpreter::DeclId_t decl = gInterpreter->GetDataMemberAtAddr(addr);
   if (decl) {
      TListOfDataMembers *globals = ((TListOfDataMembers*)(gROOT->GetListOfGlobals(kFALSE)));
      return (TGlobal*)globals->Get(decl);
   }
   // If we are actually looking for a global that is held by a global
   // pointer (for example gRandom), we need to find a pointer with the
   // correct value.
   decl = gInterpreter->GetDataMemberWithValue(addr);
   if (decl) {
      TListOfDataMembers *globals = ((TListOfDataMembers*)(gROOT->GetListOfGlobals(kFALSE)));
      return (TGlobal*)globals->Get(decl);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal routine returning, and creating if necessary, the list
/// of global function.

TListOfFunctions *TROOT::GetGlobalFunctions()
{
   if (!fGlobalFunctions) fGlobalFunctions = new TListOfFunctions(0);
   return fGlobalFunctions;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the collection of functions named "name".

TCollection *TROOT::GetListOfFunctionOverloads(const char* name) const
{
   return ((TListOfFunctions*)fGlobalFunctions)->GetListForObject(name);
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to global function by name.
/// If params != 0 it will also resolve overloading other it returns the first
/// name match.
/// If params == 0 and load is true force reading of all currently defined
/// global functions from Cling.
/// The param string must be of the form: "3189,\"aap\",1.3".

TFunction *TROOT::GetGlobalFunction(const char *function, const char *params,
                                    Bool_t load)
{
   if (!params) {
      R__LOCKGUARD(gROOTMutex);
      return (TFunction *)GetListOfGlobalFunctions(load)->FindObject(function);
   } else {
      if (!fInterpreter)
         Fatal("GetGlobalFunction", "fInterpreter not initialized");

      R__LOCKGUARD(gROOTMutex);
      TInterpreter::DeclId_t decl = gInterpreter->GetFunctionWithValues(0,
                                                                 function, params,
                                                                 false);

      if (!decl) return 0;

      TFunction *f = GetGlobalFunctions()->Get(decl);
      if (f) return f;

      Error("GetGlobalFunction",
            "\nDid not find matching TFunction <%s> with \"%s\".",
            function,params);
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to global function by name. If proto != 0
/// it will also resolve overloading. If load is true force reading
/// of all currently defined global functions from CINT (more expensive).
/// The proto string must be of the form: "int, char*, float".

TFunction *TROOT::GetGlobalFunctionWithPrototype(const char *function,
                                               const char *proto, Bool_t load)
{
   if (!proto) {
      R__LOCKGUARD(gROOTMutex);
      return (TFunction *)GetListOfGlobalFunctions(load)->FindObject(function);
   } else {
      if (!fInterpreter)
         Fatal("GetGlobalFunctionWithPrototype", "fInterpreter not initialized");

      R__LOCKGUARD(gROOTMutex);
      TInterpreter::DeclId_t decl = gInterpreter->GetFunctionWithPrototype(0,
                                                                           function, proto);

      if (!decl) return 0;

      TFunction *f = GetGlobalFunctions()->Get(decl);
      if (f) return f;

      Error("GetGlobalFunctionWithPrototype",
            "\nDid not find matching TFunction <%s> with \"%s\".",
            function,proto);
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to Geometry with name

TObject *TROOT::GetGeometry(const char *name) const
{
   return GetListOfGeometries()->FindObject(name);
}

////////////////////////////////////////////////////////////////////////////////

TCollection *TROOT::GetListOfEnums(Bool_t load /* = kTRUE */)
{
   if(!fEnums.load()) {
      R__LOCKGUARD(gROOTMutex);
      // Test again just in case, another thread did the work while we were
      // waiting.
      if (!fEnums.load()) fEnums = new TListOfEnumsWithLock(0);
   }
   if (load) {
      R__LOCKGUARD(gROOTMutex);
      (*fEnums).Load(); // Refresh the list of enums.
   }
   return fEnums.load();
}

////////////////////////////////////////////////////////////////////////////////

TCollection *TROOT::GetListOfFunctionTemplates()
{
   R__LOCKGUARD(gROOTMutex);
   if(!fFuncTemplate) {
      fFuncTemplate = new TListOfFunctionTemplates(0);
   }
   return fFuncTemplate;
}

////////////////////////////////////////////////////////////////////////////////
/// Return list containing the TGlobals currently defined.
/// Since globals are created and deleted during execution of the
/// program, we need to update the list of globals every time we
/// execute this method. However, when calling this function in
/// a (tight) loop where no interpreter symbols will be created
/// you can set load=kFALSE (default).

TCollection *TROOT::GetListOfGlobals(Bool_t load)
{
   if (!fGlobals) {
      fGlobals = new TListOfDataMembers(nullptr, TDictionary::EMemberSelection::kAlsoUsingDecls);
      // We add to the list the "funcky-fake" globals.

      // provide special functor for gROOT, while ROOT::GetROOT() does not return reference
      TGlobalMappedFunction::MakeFunctor("gROOT", "TROOT*", ROOT::GetROOT, [] {
         ROOT::GetROOT();
         return (void *)&ROOT::Internal::gROOTLocal;
      });

      TGlobalMappedFunction::MakeFunctor("gPad", "TVirtualPad*", TVirtualPad::Pad);
      TGlobalMappedFunction::MakeFunctor("gVirtualX", "TVirtualX*", TVirtualX::Instance);
      TGlobalMappedFunction::MakeFunctor("gDirectory", "TDirectory*", TDirectory::CurrentDirectory);

      // Don't let TGlobalMappedFunction delete our globals, now that we take them.
      fGlobals->AddAll(&TGlobalMappedFunction::GetEarlyRegisteredGlobals());
      TGlobalMappedFunction::GetEarlyRegisteredGlobals().SetOwner(kFALSE);
      TGlobalMappedFunction::GetEarlyRegisteredGlobals().Clear();
   }

   if (!fInterpreter)
      Fatal("GetListOfGlobals", "fInterpreter not initialized");

   if (load) fGlobals->Load();

   return fGlobals;
}

////////////////////////////////////////////////////////////////////////////////
/// Return list containing the TFunctions currently defined.
/// Since functions are created and deleted during execution of the
/// program, we need to update the list of functions every time we
/// execute this method. However, when calling this function in
/// a (tight) loop where no interpreter symbols will be created
/// you can set load=kFALSE (default).

TCollection *TROOT::GetListOfGlobalFunctions(Bool_t load)
{
   R__LOCKGUARD(gROOTMutex);

   if (!fGlobalFunctions) {
      fGlobalFunctions = new TListOfFunctions(0);
   }

   if (!fInterpreter)
      Fatal("GetListOfGlobalFunctions", "fInterpreter not initialized");

   // A thread that calls with load==true and a thread that calls with load==false
   // will conflict here (the load==true will be updating the list while the
   // other is reading it).  To solve the problem, we could use a read-write lock
   // inside the list itself.
   if (load) fGlobalFunctions->Load();

   return fGlobalFunctions;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a dynamic list giving access to all TDataTypes (typedefs)
/// currently defined.
///
/// The list is populated on demand.  Calling
/// ~~~ {.cpp}
///    gROOT->GetListOfTypes()->FindObject(nameoftype);
/// ~~~
/// will return the TDataType corresponding to 'nameoftype'.  If the
/// TDataType is not already in the list itself and the type does exist,
/// a new TDataType will be created and added to the list.
///
/// Calling
/// ~~~ {.cpp}
///    gROOT->GetListOfTypes()->ls(); // or Print()
/// ~~~
/// list only the typedefs that have been previously accessed through the
/// list (plus the builtins types).

TCollection *TROOT::GetListOfTypes(Bool_t /* load */)
{
   if (!fInterpreter)
      Fatal("GetListOfTypes", "fInterpreter not initialized");

   return fTypes;
}

////////////////////////////////////////////////////////////////////////////////
/// Get number of classes.

Int_t TROOT::GetNclasses() const
{
   return fClasses->GetSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Get number of types.

Int_t TROOT::GetNtypes() const
{
   return fTypes->GetSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Execute command when system has been idle for idleTimeInSec seconds.

void TROOT::Idle(UInt_t idleTimeInSec, const char *command)
{
   if (!fApplication.load())
      TApplication::CreateApplication();

   if (idleTimeInSec <= 0)
      (*fApplication).RemoveIdleTimer();
   else
      (*fApplication).SetIdleTimer(idleTimeInSec, command);
}

////////////////////////////////////////////////////////////////////////////////
/// Check whether className is a known class, and only autoload
/// if we can. Helper function for TROOT::IgnoreInclude().

static TClass* R__GetClassIfKnown(const char* className)
{
   // Check whether the class is available for auto-loading first:
   const char* libsToLoad = gInterpreter->GetClassSharedLibs(className);
   TClass* cla = 0;
   if (libsToLoad) {
      // trigger autoload, and only create TClass in this case.
      return TClass::GetClass(className);
   } else if (gROOT->GetListOfClasses()
              && (cla = (TClass*)gROOT->GetListOfClasses()->FindObject(className))) {
      // cla assigned in if statement
   } else if (gClassTable->FindObject(className)) {
      return TClass::GetClass(className);
   }
   return cla;
}

////////////////////////////////////////////////////////////////////////////////
/// Return 1 if the name of the given include file corresponds to a class that
///  is known to ROOT, e.g. "TLorentzVector.h" versus TLorentzVector.

Int_t TROOT::IgnoreInclude(const char *fname, const char * /*expandedfname*/)
{
   if (fname == 0) return 0;

   TString stem(fname);
   // Remove extension if any, ignore files with extension not being .h*
   Int_t where = stem.Last('.');
   if (where != kNPOS) {
      if (stem.EndsWith(".so") || stem.EndsWith(".sl") ||
          stem.EndsWith(".dl") || stem.EndsWith(".a")  ||
          stem.EndsWith(".dll", TString::kIgnoreCase))
         return 0;
      stem.Remove(where);
   }

   TString className = gSystem->BaseName(stem);
   TClass* cla = R__GetClassIfKnown(className);
   if (!cla) {
      // Try again with modifications to the file name:
      className = stem;
      className.ReplaceAll("/", "::");
      className.ReplaceAll("\\", "::");
      if (className.Contains(":::")) {
         // "C:\dir" becomes "C:::dir".
         // fname corresponds to whatever is stated after #include and
         // a full path name usually means that it's not a regular #include
         // but e.g. a ".L", so we can assume that this is not a header of
         // a class in a namespace (a global-namespace class would have been
         // detected already before).
         return 0;
      }
      cla = R__GetClassIfKnown(className);
   }

   if (!cla) {
      return 0;
   }

   // cla is valid, check wether it's actually in the header of the same name:
   if (cla->GetDeclFileLine() <= 0) return 0; // to a void an error with VisualC++
   TString decfile = gSystem->BaseName(cla->GetDeclFileName());
   if (decfile != gSystem->BaseName(fname)) {
      return 0;
   }
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize operating system interface.

void TROOT::InitSystem()
{
   if (gSystem == 0) {
#if defined(R__UNIX)
#if defined(R__HAS_COCOA)
      gSystem = new TMacOSXSystem;
#else
      gSystem = new TUnixSystem;
#endif
#elif defined(R__WIN32)
      gSystem = new TWinNTSystem;
#else
      gSystem = new TSystem;
#endif

      if (gSystem->Init())
         fprintf(stderr, "Fatal in <TROOT::InitSystem>: can't init operating system layer\n");

      if (!gSystem->HomeDirectory()) {
         fprintf(stderr, "Fatal in <TROOT::InitSystem>: HOME directory not set\n");
         fprintf(stderr, "Fix this by defining the HOME shell variable\n");
      }

      // read default files
      gEnv = new TEnv(".rootrc");

      ROOT::Internal::SetErrorSystemMsgHandler([](){ return gSystem->GetError(); });
      SetErrorHandler(DefaultErrorHandler);

      gDebug = gEnv->GetValue("Root.Debug", 0);

      if (!gEnv->GetValue("Root.ErrorHandlers", 1))
         gSystem->ResetSignals();

      // The old "Root.ZipMode" had a discrepancy between documentation vs actual meaning.
      // Also, a value with the meaning "default" wasn't available. To solved this,
      // "Root.ZipMode" was replaced by "Root.CompressionAlgorithm". Warn about usage of
      // the old value, if it's set to 0, but silently translate the setting to
      // "Root.CompressionAlgorithm" for values > 1.
      Int_t oldzipmode = gEnv->GetValue("Root.ZipMode", -1);
      if (oldzipmode == 0) {
         fprintf(stderr, "Warning in <TROOT::InitSystem>: ignoring old rootrc entry \"Root.ZipMode = 0\"!\n");
      } else {
         if (oldzipmode == -1 || oldzipmode == 1) {
            // Not set or default value, use "default" for "Root.CompressionAlgorithm":
            oldzipmode = 0;
         }
         // else keep the old zipmode (e.g. "3") as "Root.CompressionAlgorithm"
         // if "Root.CompressionAlgorithm" isn't set; see below.
      }

      Int_t zipmode = gEnv->GetValue("Root.CompressionAlgorithm", oldzipmode);
      if (zipmode != 0) R__SetZipMode(zipmode);

      const char *sdeb;
      if ((sdeb = gSystem->Getenv("ROOTDEBUG")))
         gDebug = atoi(sdeb);

      if (gDebug > 0 && isatty(2))
         fprintf(stderr, "Info in <TROOT::InitSystem>: running with gDebug = %d\n", gDebug);

      if (gEnv->GetValue("Root.MemStat", 0))
         TStorage::EnableStatistics();
      int msize = gEnv->GetValue("Root.MemStat.size", -1);
      int mcnt  = gEnv->GetValue("Root.MemStat.cnt", -1);
      if (msize != -1 || mcnt != -1)
         TStorage::EnableStatistics(msize, mcnt);

      fgMemCheck = gEnv->GetValue("Root.MemCheck", 0);

#if defined(R__HAS_COCOA)
      // create and delete a dummy TUrl so that TObjectStat table does not contain
      // objects that are deleted after recording is turned-off (in next line),
      // like the TUrl::fgSpecialProtocols list entries which are created in the
      // TMacOSXSystem ctor.
      { TUrl dummy("/dummy"); }
#endif
      TObject::SetObjectStat(gEnv->GetValue("Root.ObjectStat", 0));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Load and initialize thread library.

void TROOT::InitThreads()
{
   if (gEnv->GetValue("Root.UseThreads", 0) || gEnv->GetValue("Root.EnableThreadSafety", 0)) {
      ROOT::EnableThreadSafety();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the interpreter. Should be called only after main(),
/// to make sure LLVM/Clang is fully initialized.

void TROOT::InitInterpreter()
{
   // usedToIdentifyRootClingByDlSym is available when TROOT is part of
   // rootcling.
   if (!dlsym(RTLD_DEFAULT, "usedToIdentifyRootClingByDlSym")
       && !dlsym(RTLD_DEFAULT, "usedToIdentifyStaticRoot")) {
      char *libRIO = gSystem->DynamicPathName("libRIO");
      void *libRIOHandle = dlopen(libRIO, RTLD_NOW|RTLD_GLOBAL);
      delete [] libRIO;
      if (!libRIOHandle) {
         TString err = dlerror();
         fprintf(stderr, "Fatal in <TROOT::InitInterpreter>: cannot load library %s\n", err.Data());
         exit(1);
      }

      char *libcling = gSystem->DynamicPathName("libCling");
      gInterpreterLib = dlopen(libcling, RTLD_LAZY|RTLD_LOCAL);
      delete [] libcling;

      if (!gInterpreterLib) {
         TString err = dlerror();
         fprintf(stderr, "Fatal in <TROOT::InitInterpreter>: cannot load library %s\n", err.Data());
         exit(1);
      }
      dlerror();   // reset error message
   } else {
      gInterpreterLib = RTLD_DEFAULT;
   }
   CreateInterpreter_t *CreateInterpreter = (CreateInterpreter_t*) dlsym(gInterpreterLib, "CreateInterpreter");
   if (!CreateInterpreter) {
      TString err = dlerror();
      fprintf(stderr, "Fatal in <TROOT::InitInterpreter>: cannot load symbol %s\n", err.Data());
      exit(1);
   }
   // Schedule the destruction of TROOT.
   atexit(at_exit_of_TROOT);

   gDestroyInterpreter = (DestroyInterpreter_t*) dlsym(gInterpreterLib, "DestroyInterpreter");
   if (!gDestroyInterpreter) {
      TString err = dlerror();
      fprintf(stderr, "Fatal in <TROOT::InitInterpreter>: cannot load symbol %s\n", err.Data());
      exit(1);
   }

   const char *interpArgs[] = {
#ifdef NDEBUG
      "-DNDEBUG",
#else
      "-UNDEBUG",
#endif
#ifdef DEBUG
      "-DDEBUG",
#else
      "-UDEBUG",
#endif
#ifdef _DEBUG
      "-D_DEBUG",
#else
      "-U_DEBUG",
#endif
      nullptr};

   fInterpreter = CreateInterpreter(gInterpreterLib, interpArgs);

   fCleanups->Add(fInterpreter);
   fInterpreter->SetBit(kMustCleanup);

   fgRootInit = kTRUE;

   // initialize gClassTable is not already done
   if (!gClassTable)
      new TClassTable;

   // Initialize all registered dictionaries.
   for (std::vector<ModuleHeaderInfo_t>::const_iterator
           li = GetModuleHeaderInfoBuffer().begin(),
           le = GetModuleHeaderInfoBuffer().end(); li != le; ++li) {
         // process buffered module registrations
      fInterpreter->RegisterModule(li->fModuleName,
                                   li->fHeaders,
                                   li->fIncludePaths,
                                   li->fPayloadCode,
                                   li->fFwdDeclCode,
                                   li->fTriggerFunc,
                                   li->fFwdNargsToKeepColl,
                                   li->fClassesHeaders,
                                   kTRUE /*lateRegistration*/,
                                   li->fHasCxxModule);
   }
   GetModuleHeaderInfoBuffer().clear();

   fInterpreter->Initialize();
}

////////////////////////////////////////////////////////////////////////////////
/// Helper function used by TClass::GetClass().
/// This function attempts to load the dictionary for 'classname'
/// either from the TClassTable or from the list of generator.
/// If silent is 'true', do not warn about missing dictionary for the class.
/// (typically used for class that are used only for transient members)
///
/// The 'requestedname' is expected to be already normalized.

TClass *TROOT::LoadClass(const char *requestedname, Bool_t silent) const
{
   return TClass::LoadClass(requestedname, silent);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if class "classname" is known to the interpreter (in fact,
/// this check is not needed anymore, so classname is ignored). If
/// not it will load library "libname". If the library couldn't be found with original
/// libname and if the name was not prefixed with lib, try to prefix with "lib" and search again.
/// If DynamicPathName still couldn't find the library, return -1.
/// If check is true it will only check if libname exists and is
/// readable.
/// Returns 0 on successful loading, -1 in case libname does not
/// exist or in case of error and -2 in case of version mismatch.

Int_t TROOT::LoadClass(const char * /*classname*/, const char *libname,
                       Bool_t check)
{
   TString lib(libname);

   // Check if libname exists in path or not
   if (char *path = gSystem->DynamicPathName(lib, kTRUE)) {
      // If check == true, only check if it exists and if it's readable
      if (check) {
         delete [] path;
         return 0;
      }

      // If check == false, try to load the library
      else {
         int err = gSystem->Load(path, 0, kTRUE);
         delete [] path;

         // TSystem::Load returns 1 when the library was already loaded, return success in this case.
         if (err == 1)
            err = 0;
         return err;
      }
   } else {
      // This is the branch where libname didn't exist
      if (check) {
         FileStat_t stat;
         if (!gSystem->GetPathInfo(libname, stat) && (R_ISREG(stat.fMode) &&
             !gSystem->AccessPathName(libname, kReadPermission)))
            return 0;
      }

      // Take care of user who didn't write the whole name
      if (!lib.BeginsWith("lib")) {
         lib = "lib" + lib;
         return LoadClass("", lib.Data(), check);
      }
   }

   // Execution reaches here when library was prefixed with lib, check is false and couldn't find
   // the library name.
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the file is local and is (likely) to be a ROOT file

Bool_t TROOT::IsRootFile(const char *filename) const
{
   Bool_t result = kFALSE;
   FILE *mayberootfile = fopen(filename,"rb");
   if (mayberootfile) {
      char header[5];
      if (fgets(header,5,mayberootfile)) {
         result = strncmp(header,"root",4)==0;
      }
      fclose(mayberootfile);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// To list all objects of the application.
/// Loop on all objects created in the ROOT linked lists.
/// Objects may be files and windows or any other object directly
/// attached to the ROOT linked list.

void TROOT::ls(Option_t *option) const
{
//   TObject::SetDirLevel();
//   GetList()->R__FOR_EACH(TObject,ls)(option);
   TDirectory::ls(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Load a macro in the interpreter's memory. Equivalent to the command line
/// command ".L filename". If the filename has "+" or "++" appended
/// the macro will be compiled by ACLiC. The filename must have the format:
/// [path/]macro.C[+|++[g|O]].
/// The possible error codes are defined by TInterpreter::EErrorCode.
/// If check is true it will only check if filename exists and is
/// readable.
/// Returns 0 on successful loading and -1 in case filename does not
/// exist or in case of error.

Int_t TROOT::LoadMacro(const char *filename, int *error, Bool_t check)
{
   Int_t err = -1;
   Int_t lerr, *terr;
   if (error)
      terr = error;
   else
      terr = &lerr;

   if (fInterpreter) {
      TString aclicMode;
      TString arguments;
      TString io;
      TString fname = gSystem->SplitAclicMode(filename, aclicMode, arguments, io);

      if (arguments.Length()) {
         Warning("LoadMacro", "argument(%s) ignored in %s", arguments.Data(), GetMacroPath());
      }
      char *mac = gSystem->Which(GetMacroPath(), fname, kReadPermission);
      if (!mac) {
         if (!check)
            Error("LoadMacro", "macro %s not found in path %s", fname.Data(), GetMacroPath());
         *terr = TInterpreter::kFatal;
      } else {
         err = 0;
         if (!check) {
            fname = mac;
            fname += aclicMode;
            fname += io;
            gInterpreter->LoadMacro(fname.Data(), (TInterpreter::EErrorCode*)terr);
            if (*terr)
               err = -1;
         }
      }
      delete [] mac;
   }
   return err;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute a macro in the interpreter. Equivalent to the command line
/// command ".x filename". If the filename has "+" or "++" appended
/// the macro will be compiled by ACLiC. The filename must have the format:
/// [path/]macro.C[+|++[g|O]][(args)].
/// The possible error codes are defined by TInterpreter::EErrorCode.
/// If padUpdate is true (default) update the current pad.
/// Returns the macro return value.

Long_t TROOT::Macro(const char *filename, Int_t *error, Bool_t padUpdate)
{
   Long_t result = 0;

   if (fInterpreter) {
      TString aclicMode;
      TString arguments;
      TString io;
      TString fname = gSystem->SplitAclicMode(filename, aclicMode, arguments, io);

      char *mac = gSystem->Which(GetMacroPath(), fname, kReadPermission);
      if (!mac) {
         Error("Macro", "macro %s not found in path %s", fname.Data(), GetMacroPath());
         if (error)
            *error = TInterpreter::kFatal;
      } else {
         fname = mac;
         fname += aclicMode;
         fname += arguments;
         fname += io;
         result = gInterpreter->ExecuteMacro(fname, (TInterpreter::EErrorCode*)error);
      }
      delete [] mac;

      if (padUpdate && gPad)
         gPad->Update();
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Process message id called by obj.

void  TROOT::Message(Int_t id, const TObject *obj)
{
   TIter next(fMessageHandlers);
   TMessageHandler *mh;
   while ((mh = (TMessageHandler*)next())) {
      mh->HandleMessage(id,obj);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Process interpreter command via TApplication::ProcessLine().
/// On Win32 the line will be processed asynchronously by sending
/// it to the CINT interpreter thread. For explicit synchronous processing
/// use ProcessLineSync(). On non-Win32 platforms there is no difference
/// between ProcessLine() and ProcessLineSync().
/// The possible error codes are defined by TInterpreter::EErrorCode. In
/// particular, error will equal to TInterpreter::kProcessing until the
/// CINT interpreted thread has finished executing the line.
/// Returns the result of the command, cast to a Long_t.

Long_t TROOT::ProcessLine(const char *line, Int_t *error)
{
   TString sline = line;
   sline = sline.Strip(TString::kBoth);

   if (!fApplication.load())
      TApplication::CreateApplication();

   return (*fApplication).ProcessLine(sline, kFALSE, error);
}

////////////////////////////////////////////////////////////////////////////////
/// Process interpreter command via TApplication::ProcessLine().
/// On Win32 the line will be processed synchronously (i.e. it will
/// only return when the CINT interpreter thread has finished executing
/// the line). On non-Win32 platforms there is no difference between
/// ProcessLine() and ProcessLineSync().
/// The possible error codes are defined by TInterpreter::EErrorCode.
/// Returns the result of the command, cast to a Long_t.

Long_t TROOT::ProcessLineSync(const char *line, Int_t *error)
{
   TString sline = line;
   sline = sline.Strip(TString::kBoth);

   if (!fApplication.load())
      TApplication::CreateApplication();

   return (*fApplication).ProcessLine(sline, kTRUE, error);
}

////////////////////////////////////////////////////////////////////////////////
/// Process interpreter command directly via CINT interpreter.
/// Only executable statements are allowed (no variable declarations),
/// In all other cases use TROOT::ProcessLine().
/// The possible error codes are defined by TInterpreter::EErrorCode.

Long_t TROOT::ProcessLineFast(const char *line, Int_t *error)
{
   TString sline = line;
   sline = sline.Strip(TString::kBoth);

   if (!fApplication.load())
      TApplication::CreateApplication();

   Long_t result = 0;

   if (fInterpreter) {
      TInterpreter::EErrorCode *code = (TInterpreter::EErrorCode*)error;
      result = gInterpreter->Calc(sline, code);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Read Git commit information and branch name from the
/// etc/gitinfo.txt file.

void TROOT::ReadGitInfo()
{
#ifdef ROOT_GIT_COMMIT
   fGitCommit = ROOT_GIT_COMMIT;
#endif
#ifdef ROOT_GIT_BRANCH
   fGitBranch = ROOT_GIT_BRANCH;
#endif

   TString gitinfo = "gitinfo.txt";
   char *filename = gSystem->ConcatFileName(TROOT::GetEtcDir(), gitinfo);

   FILE *fp = fopen(filename, "r");
   if (fp) {
      TString s;
      // read branch name
      s.Gets(fp);
      fGitBranch = s;
      // read commit SHA1
      s.Gets(fp);
      fGitCommit = s;
      // read date/time make was run
      s.Gets(fp);
      fGitDate = s;
      fclose(fp);
   }
   delete [] filename;
}

Bool_t &GetReadingObject() {
   TTHREAD_TLS(Bool_t) fgReadingObject = false;
   return fgReadingObject;
}

////////////////////////////////////////////////////////////////////////////////
/// Deprecated (will be removed in next release).

Bool_t TROOT::ReadingObject() const
{
   return GetReadingObject();
}

void TROOT::SetReadingObject(Bool_t flag)
{
   GetReadingObject() = flag;
}


////////////////////////////////////////////////////////////////////////////////
/// Return date/time make was run.

const char *TROOT::GetGitDate()
{
   if (fGitDate == "") {
      Int_t iday,imonth,iyear, ihour, imin;
      static const char *months[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
      Int_t idate = gROOT->GetBuiltDate();
      Int_t itime = gROOT->GetBuiltTime();
      iday   = idate%100;
      imonth = (idate/100)%100;
      iyear  = idate/10000;
      ihour  = itime/100;
      imin   = itime%100;
      fGitDate.Form("%s %02d %4d, %02d:%02d:00", months[imonth-1], iday, iyear, ihour, imin);
   }
   return fGitDate;
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively remove this object from the list of Cleanups.
/// Typically RecursiveRemove is implemented by classes that can contain
/// mulitple references to a same object or shared ownership of the object
/// with others.

void TROOT::RecursiveRemove(TObject *obj)
{
   R__READ_LOCKGUARD(ROOT::gCoreMutex);

   fCleanups->RecursiveRemove(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Refresh all browsers. Call this method when some command line
/// command or script has changed the browser contents. Not needed
/// for objects that have the kMustCleanup bit set. Most useful to
/// update browsers that show the file system or other objects external
/// to the running ROOT session.

void TROOT::RefreshBrowsers()
{
   TIter next(GetListOfBrowsers());
   TBrowser *b;
   while ((b = (TBrowser*) next()))
      b->SetRefreshFlag(kTRUE);
}
////////////////////////////////////////////////////////////////////////////////
/// Insure that the files, canvases and sockets are closed.

static void CallCloseFiles()
{
   if (TROOT::Initialized() && ROOT::Internal::gROOTLocal) {
      gROOT->CloseFiles();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called by static dictionary initialization to register clang modules
/// for headers. Calls TCling::RegisterModule() unless gCling
/// is NULL, i.e. during startup, where the information is buffered in
/// the static GetModuleHeaderInfoBuffer().

void TROOT::RegisterModule(const char* modulename,
                           const char** headers,
                           const char** includePaths,
                           const char* payloadCode,
                           const char* fwdDeclCode,
                           void (*triggerFunc)(),
                           const TInterpreter::FwdDeclArgsToKeepCollection_t& fwdDeclsArgToSkip,
                           const char** classesHeaders,
                           bool hasCxxModule)
{

   // First a side track to insure proper end of process behavior.

   // Register for each loaded dictionary (and thus for each library),
   // that we need to Close the ROOT files as soon as this library
   // might start being unloaded after main.
   //
   // By calling atexit here (rather than directly from within the
   // library) we make sure that this is not called if the library is
   // 'only' dlclosed.

   // On Ubuntu the linker strips the unused libraries.  Eventhough
   // stressHistogram is explicitly linked against libNet, it is not
   // retained and thus is loaded only as needed in the middle part of
   // the execution.  Concretely this also means that it is loaded
   // *after* the construction of the TApplication object and thus
   // after the registration (atexit) of the EndOfProcessCleanups
   // routine.  Consequently, after the end of main, libNet is
   // unloaded before EndOfProcessCleanups is called.  When
   // EndOfProcessCleanups is executed it indirectly needs the TClass
   // for TSocket and its search will use resources that have already
   // been unloaded (technically the function static in TUnixSystem's
   // DynamicPath and the dictionary from libNet).

   // Similarly, the ordering (before this commit) was broken in the
   // following case:

   //    TApplication creation (EndOfProcessCleanups registration)
   //    load UserLibrary
   //    create TFile
   //    Append UserObject to TFile

   // and after the end of main the order of execution was

   //    unload UserLibrary
   //    call EndOfProcessCleanups
   //       Write the TFile
   //         attempt to write the user object.
   //    ....

   // where what we need is to have the files closen/written before
   // the unloading of the library.

   // To solve the problem we now register an atexit function for
   // every dictionary thus making sure there is at least one executed
   // before the first library tear down after main.

   // If atexit is called directly within a library's code, the
   // function will called *either* when the library is 'dlclose'd or
   // after then end of main (whichever comes first).  We do *not*
   // want the files to be closed whenever a library is unloaded via
   // dlclose.  To avoid this, we add the function (CallCloseFiles)
   // from the dictionary indirectly (via ROOT::RegisterModule).  In
   // this case the function will only only be called either when
   // libCore is 'dlclose'd or right after the end of main.

   atexit(CallCloseFiles);

   // Now register with TCling.
   if (TROOT::Initialized()) {
      gCling->RegisterModule(modulename, headers, includePaths, payloadCode, fwdDeclCode, triggerFunc,
                             fwdDeclsArgToSkip, classesHeaders, false, hasCxxModule);
   } else {
      GetModuleHeaderInfoBuffer().push_back(ModuleHeaderInfo_t(modulename, headers, includePaths, payloadCode,
                                                               fwdDeclCode, triggerFunc, fwdDeclsArgToSkip,
                                                               classesHeaders, hasCxxModule));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove an object from the in-memory list.
///    Since TROOT is global resource, this is lock protected.

TObject *TROOT::Remove(TObject* obj)
{
   R__LOCKGUARD(gROOTMutex);
   return TDirectory::Remove(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a class from the list and map of classes.
/// This routine is deprecated, use TClass::RemoveClass directly.

void TROOT::RemoveClass(TClass *oldcl)
{
   TClass::RemoveClass(oldcl);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete all global interpreter objects created since the last call to Reset
///
/// If option="a" is set reset to startup context (i.e. unload also
/// all loaded files, classes, structs, typedefs, etc.).
///
/// This function is typically used at the beginning (or end) of an unnamed macro
/// to clean the environment.
///
/// IMPORTANT WARNING:
/// Do not use this call from within any function (neither compiled nor
/// interpreted.  This should only be used from a unnamed macro
/// (which starts with a { (curly braces)  ).  For example, using TROOT::Reset
/// from within an interpreted function will lead to the unloading of the
/// dictionary and source file, including the one defining the function being
/// executed.
///

void TROOT::Reset(Option_t *option)
{
   if (IsExecutingMacro()) return;  //True when TMacro::Exec runs
   if (fInterpreter) {
      if (!strncmp(option, "a", 1)) {
         fInterpreter->Reset();
         fInterpreter->SaveContext();
      } else
         gInterpreter->ResetGlobals();

      if (fGlobals) fGlobals->Unload();
      if (fGlobalFunctions) fGlobalFunctions->Unload();

      SaveContext();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save the current interpreter context.

void TROOT::SaveContext()
{
   if (fInterpreter)
      gInterpreter->SaveGlobalsContext();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the default graphical cut class name for the graphics editor
/// By default the graphics editor creates an instance of a class TCutG.
/// This function may be called to specify a different class that MUST
/// derive from TCutG

void TROOT::SetCutClassName(const char *name)
{
   if (!name) {
      Error("SetCutClassName","Invalid class name");
      return;
   }
   TClass *cl = TClass::GetClass(name);
   if (!cl) {
      Error("SetCutClassName","Unknown class:%s",name);
      return;
   }
   if (!cl->InheritsFrom("TCutG")) {
      Error("SetCutClassName","Class:%s does not derive from TCutG",name);
      return;
   }
   fCutClassName = name;
}

////////////////////////////////////////////////////////////////////////////////
/// Set editor mode

void TROOT::SetEditorMode(const char *mode)
{
   fEditorMode = 0;
   if (!mode[0]) return;
   if (!strcmp(mode,"Arc"))      {fEditorMode = kArc;        return;}
   if (!strcmp(mode,"Line"))     {fEditorMode = kLine;       return;}
   if (!strcmp(mode,"Arrow"))    {fEditorMode = kArrow;      return;}
   if (!strcmp(mode,"Button"))   {fEditorMode = kButton;     return;}
   if (!strcmp(mode,"Diamond"))  {fEditorMode = kDiamond;    return;}
   if (!strcmp(mode,"Ellipse"))  {fEditorMode = kEllipse;    return;}
   if (!strcmp(mode,"Pad"))      {fEditorMode = kPad;        return;}
   if (!strcmp(mode,"Pave"))     {fEditorMode = kPave;       return;}
   if (!strcmp(mode,"PaveLabel")){fEditorMode = kPaveLabel;  return;}
   if (!strcmp(mode,"PaveText")) {fEditorMode = kPaveText;   return;}
   if (!strcmp(mode,"PavesText")){fEditorMode = kPavesText;  return;}
   if (!strcmp(mode,"PolyLine")) {fEditorMode = kPolyLine;   return;}
   if (!strcmp(mode,"CurlyLine")){fEditorMode = kCurlyLine;  return;}
   if (!strcmp(mode,"CurlyArc")) {fEditorMode = kCurlyArc;   return;}
   if (!strcmp(mode,"Text"))     {fEditorMode = kText;       return;}
   if (!strcmp(mode,"Marker"))   {fEditorMode = kMarker;     return;}
   if (!strcmp(mode,"CutG"))     {fEditorMode = kCutG;       return;}
}

////////////////////////////////////////////////////////////////////////////////
/// Change current style to style with name stylename

void TROOT::SetStyle(const char *stylename)
{
   TString style_name = stylename;

   TStyle *style = GetStyle(style_name);
   if (style) style->cd();
   else       Error("SetStyle","Unknown style:%s",style_name.Data());
}


//-------- Static Member Functions ---------------------------------------------


////////////////////////////////////////////////////////////////////////////////
/// Decrease the indentation level for ls().

Int_t TROOT::DecreaseDirLevel()
{
   return --fgDirLevel;
}

////////////////////////////////////////////////////////////////////////////////
///return directory level

Int_t TROOT::GetDirLevel()
{
   return fgDirLevel;
}

////////////////////////////////////////////////////////////////////////////////
/// Get macro search path. Static utility function.

const char *TROOT::GetMacroPath()
{
   TString &macroPath = ROOT::GetMacroPath();

   if (macroPath.Length() == 0) {
      macroPath = gEnv->GetValue("Root.MacroPath", (char*)0);
#if defined(R__WIN32)
      macroPath.ReplaceAll("; ", ";");
#else
      macroPath.ReplaceAll(": ", ":");
#endif
      if (macroPath.Length() == 0)
#if !defined(R__WIN32)
         macroPath = ".:" + TROOT::GetMacroDir();
#else
         macroPath = ".;" + TROOT::GetMacroDir();
#endif
   }

   return macroPath;
}

////////////////////////////////////////////////////////////////////////////////
/// Set or extend the macro search path. Static utility function.
/// If newpath=0 or "" reset to value specified in the rootrc file.

void TROOT::SetMacroPath(const char *newpath)
{
   TString &macroPath = ROOT::GetMacroPath();

   if (!newpath || !*newpath)
      macroPath = "";
   else
      macroPath = newpath;
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Specify where web graphics shall be rendered
///
/// The input parameter `webdisplay` defines where web graphics is rendered.
/// `webdisplay` parameter may contain:
///
///  - "off": turns off the web display and comes back to normal graphics in
///    interactive mode.
///  - "batch": turns the web display in batch mode. It can be prepended with
///    another string which is considered as the new current web display.
///  - "nobatch": turns the web display in interactive mode. It can be
///    prepended with another string which is considered as the new current web display.
///
/// If the option "off" is not set, this method turns the normal graphics to
/// "Batch" to avoid the loading of local graphics libraries.

void TROOT::SetWebDisplay(const char *webdisplay)
{
   const char *wd = webdisplay;
   if (!wd)
      wd = "";

   if (!strcmp(wd, "off")) {
      fIsWebDisplay = kFALSE;
      fIsWebDisplayBatch = kFALSE;
      fWebDisplay = "";
   } else {
      fIsWebDisplay = kTRUE;
      if (!strncmp(wd, "batch", 5)) {
         fIsWebDisplayBatch = kTRUE;
         wd += 5;
      } else if (!strncmp(wd, "nobatch", 7)) {
         fIsWebDisplayBatch = kFALSE;
         wd += 7;
      } else {
         fIsWebDisplayBatch = kFALSE;
      }
      fWebDisplay = wd;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Increase the indentation level for ls().

Int_t TROOT::IncreaseDirLevel()
{
   return ++fgDirLevel;
}

////////////////////////////////////////////////////////////////////////////////
/// Functions used by ls() to indent an object hierarchy.

void TROOT::IndentLevel()
{
   for (int i = 0; i < fgDirLevel; i++) std::cout.put(' ');
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize ROOT explicitly.

void TROOT::Initialize() {
   (void) gROOT;
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if the TROOT object has been initialized.

Bool_t TROOT::Initialized()
{
   return fgRootInit;
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if the memory leak checker is on.

Bool_t TROOT::MemCheck()
{
   return fgMemCheck;
}

////////////////////////////////////////////////////////////////////////////////
/// Return Indentation level for ls().

void TROOT::SetDirLevel(Int_t level)
{
   fgDirLevel = level;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert version code to an integer, i.e. 331527 -> 51507.

Int_t TROOT::ConvertVersionCode2Int(Int_t code)
{
   return 10000*(code>>16) + 100*((code&65280)>>8) + (code&255);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert version as an integer to version code as used in RVersion.h.

Int_t TROOT::ConvertVersionInt2Code(Int_t v)
{
   int a = v/10000;
   int b = (v - a*10000)/100;
   int c = v - a*10000 - b*100;
   return (a << 16) + (b << 8) + c;
}

////////////////////////////////////////////////////////////////////////////////
/// Return ROOT version code as defined in RVersion.h.

Int_t TROOT::RootVersionCode()
{
   return ROOT_VERSION_CODE;
}
////////////////////////////////////////////////////////////////////////////////
/// Provide command line arguments to the interpreter construction.
/// These arguments are added to the existing flags (e.g. `-DNDEBUG`).
/// They are evaluated once per process, at the time where TROOT (and thus
/// TInterpreter) is constructed.
/// Returns the new flags.

const std::vector<std::string> &TROOT::AddExtraInterpreterArgs(const std::vector<std::string> &args) {
   static std::vector<std::string> sArgs = {};
   sArgs.insert(sArgs.begin(), args.begin(), args.end());
   return sArgs;
}

////////////////////////////////////////////////////////////////////////////////
/// INTERNAL function!
/// Used by rootcling to inject interpreter arguments through a C-interface layer.

const char**& TROOT::GetExtraInterpreterArgs() {
   static const char** extraInterpArgs = 0;
   return extraInterpArgs;
}

////////////////////////////////////////////////////////////////////////////////

#ifdef ROOTPREFIX
static Bool_t IgnorePrefix() {
   static Bool_t ignorePrefix = gSystem->Getenv("ROOTIGNOREPREFIX");
   return ignorePrefix;
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// Get the rootsys directory in the installation. Static utility function.

const TString& TROOT::GetRootSys() {
   // Avoid returning a reference to a temporary because of the conversion
   // between std::string and TString.
   const static TString rootsys = ROOT::FoundationUtils::GetRootSys();
   return rootsys;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the binary directory in the installation. Static utility function.

const TString& TROOT::GetBinDir() {
#ifdef ROOTBINDIR
   if (IgnorePrefix()) {
#endif
      static TString rootbindir;
      if (rootbindir.IsNull()) {
         rootbindir = "bin";
         gSystem->PrependPathName(GetRootSys(), rootbindir);
      }
      return rootbindir;
#ifdef ROOTBINDIR
   } else {
      const static TString rootbindir = ROOTBINDIR;
      return rootbindir;
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Get the library directory in the installation. Static utility function.

const TString& TROOT::GetLibDir() {
#ifdef ROOTLIBDIR
   if (IgnorePrefix()) {
#endif
      static TString rootlibdir;
      if (rootlibdir.IsNull()) {
         rootlibdir = "lib";
         gSystem->PrependPathName(GetRootSys(), rootlibdir);
      }
      return rootlibdir;
#ifdef ROOTLIBDIR
   } else {
      const static TString rootlibdir = ROOTLIBDIR;
      return rootlibdir;
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Get the include directory in the installation. Static utility function.

const TString& TROOT::GetIncludeDir() {
   // Avoid returning a reference to a temporary because of the conversion
   // between std::string and TString.
   const static TString includedir = ROOT::FoundationUtils::GetIncludeDir();
   return includedir;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the sysconfig directory in the installation. Static utility function.

const TString& TROOT::GetEtcDir() {
   // Avoid returning a reference to a temporary because of the conversion
   // between std::string and TString.
   const static TString etcdir = ROOT::FoundationUtils::GetEtcDir();
   return etcdir;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the data directory in the installation. Static utility function.

const TString& TROOT::GetDataDir() {
#ifdef ROOTDATADIR
   if (IgnorePrefix()) {
#endif
      return GetRootSys();
#ifdef ROOTDATADIR
   } else {
      const static TString rootdatadir = ROOTDATADIR;
      return rootdatadir;
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Get the documentation directory in the installation. Static utility function.

const TString& TROOT::GetDocDir() {
#ifdef ROOTDOCDIR
   if (IgnorePrefix()) {
#endif
      return GetRootSys();
#ifdef ROOTDOCDIR
   } else {
      const static TString rootdocdir = ROOTDOCDIR;
      return rootdocdir;
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Get the macro directory in the installation. Static utility function.

const TString& TROOT::GetMacroDir() {
#ifdef ROOTMACRODIR
   if (IgnorePrefix()) {
#endif
      static TString rootmacrodir;
      if (rootmacrodir.IsNull()) {
         rootmacrodir = "macros";
         gSystem->PrependPathName(GetRootSys(), rootmacrodir);
      }
      return rootmacrodir;
#ifdef ROOTMACRODIR
   } else {
      const static TString rootmacrodir = ROOTMACRODIR;
      return rootmacrodir;
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Get the tutorials directory in the installation. Static utility function.

const TString& TROOT::GetTutorialDir() {
#ifdef ROOTTUTDIR
   if (IgnorePrefix()) {
#endif
      static TString roottutdir;
      if (roottutdir.IsNull()) {
         roottutdir = "tutorials";
         gSystem->PrependPathName(GetRootSys(), roottutdir);
      }
      return roottutdir;
#ifdef ROOTTUTDIR
   } else {
      const static TString roottutdir = ROOTTUTDIR;
      return roottutdir;
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Shut down ROOT.

void TROOT::ShutDown()
{
   if (gROOT)
      gROOT->EndOfProcessCleanups();
   else if (gInterpreter)
      gInterpreter->ShutDown();
}

////////////////////////////////////////////////////////////////////////////////
/// Get the source directory in the installation. Static utility function.

const TString& TROOT::GetSourceDir() {
#ifdef ROOTSRCDIR
   if (IgnorePrefix()) {
#endif
      static TString rootsrcdir;
      if (rootsrcdir.IsNull()) {
         rootsrcdir = "src";
         gSystem->PrependPathName(GetRootSys(), rootsrcdir);
      }
      return rootsrcdir;
#ifdef ROOTSRCDIR
   } else {
      const static TString rootsrcdir = ROOTSRCDIR;
      return rootsrcdir;
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Get the icon path in the installation. Static utility function.

const TString& TROOT::GetIconPath() {
#ifdef ROOTICONPATH
   if (IgnorePrefix()) {
#endif
      static TString rooticonpath;
      if (rooticonpath.IsNull()) {
         rooticonpath = "icons";
         gSystem->PrependPathName(GetRootSys(), rooticonpath);
      }
      return rooticonpath;
#ifdef ROOTICONPATH
   } else {
      const static TString rooticonpath = ROOTICONPATH;
      return rooticonpath;
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Get the fonts directory in the installation. Static utility function.

const TString& TROOT::GetTTFFontDir() {
#ifdef TTFFONTDIR
   if (IgnorePrefix()) {
#endif
      static TString ttffontdir;
      if (ttffontdir.IsNull()) {
         ttffontdir = "fonts";
         gSystem->PrependPathName(GetRootSys(), ttffontdir);
      }
      return ttffontdir;
#ifdef TTFFONTDIR
   } else {
      const static TString ttffontdir = TTFFONTDIR;
      return ttffontdir;
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Get the tutorials directory in the installation. Static utility function.
/// Backward compatibility function - do not use for new code

const char *TROOT::GetTutorialsDir() {
   return GetTutorialDir();
}
