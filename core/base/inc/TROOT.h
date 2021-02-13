// @(#)root/base:$Id$
// Author: Rene Brun   08/12/94

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TROOT
#define ROOT_TROOT


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TROOT                                                                //
//                                                                      //
// The TROOT object is the entry point to the system.                   //
// The single instance of TROOT is accessible via the global gROOT.     //
// Using the gROOT pointer one has access to basically every object     //
// created in a ROOT based program. The TROOT object is essentially a   //
// "dispatcher" with several lists pointing to the ROOT main objects.   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TDirectory.h"
// #include "TList.h"  // included in TDirectory.h
// #include "RConfigure.h"  // included via Rtypes.h

#include <atomic>
#include <string>
#include <vector>
#include <utility>

class TClass;
class TCanvas;
class TColor;
class TDataType;
class TFile;
class TStyle;
class TVirtualPad;
class TApplication;
class TInterpreter;
class TBrowser;
class TGlobal;
class TFunction;
class TFolder;
class TPluginManager;
class TProcessUUID;
class TClassGenerator;
class TVirtualMutex;
class TROOT;
class TListOfDataMembers;
class TListOfEnums;
class TListOfFunctions;
class TListOfFunctionTemplates;
class TSeqCollection;
class TFunctionTemplate;
class TGlobalMappedFunction;

R__EXTERN TVirtualMutex *gROOTMutex;

namespace ROOT {
namespace Internal {
   class TROOTAllocator;

   TROOT *GetROOT2();

   // Manage parallel branch processing
   void EnableParBranchProcessing();
   void DisableParBranchProcessing();
   Bool_t IsParBranchProcessingEnabled();
   class TParBranchProcessingRAII {
   public:
      TParBranchProcessingRAII()  { EnableParBranchProcessing();  }
      ~TParBranchProcessingRAII() { DisableParBranchProcessing(); }
   };
} } // End ROOT::Internal

namespace ROOT {
   // Enable support for multi-threading within the ROOT code,
   // in particular, enables the global mutex to make ROOT thread safe/aware.
   void EnableThreadSafety();
   /// \brief Enable ROOT's implicit multi-threading for all objects and methods that provide an internal
   /// parallelisation mechanism.
   void EnableImplicitMT(UInt_t numthreads = 0);
   void DisableImplicitMT();
   Bool_t IsImplicitMTEnabled();
   UInt_t GetThreadPoolSize();
}

class TROOT : public TDirectory {

friend class TCling;
friend TROOT *ROOT::Internal::GetROOT2();

private:
   Int_t           fLineIsProcessing;     //To synchronize multi-threads

   static Int_t    fgDirLevel;            //Indentation level for ls()
   static Bool_t   fgRootInit;            //Singleton initialization flag
   static Bool_t   fgMemCheck;            //Turn on memory leak checker

   TROOT(const TROOT&) = delete;
   TROOT& operator=(const TROOT&) = delete;

protected:
   typedef std::atomic<TListOfEnums*> AListOfEnums_t;

   TString         fConfigOptions;        //ROOT ./configure set build options
   TString         fConfigFeatures;       //ROOT ./configure detected build features
   TString         fVersion;              //ROOT version (from CMZ VERSQQ) ex 0.05/01
   Int_t           fVersionInt;           //ROOT version in integer format (501)
   Int_t           fVersionCode;          //ROOT version code as used in RVersion.h
   Int_t           fVersionDate;          //Date of ROOT version (ex 951226)
   Int_t           fVersionTime;          //Time of ROOT version (ex 1152)
   Int_t           fBuiltDate;            //Date of ROOT built
   Int_t           fBuiltTime;            //Time of ROOT built
   TString         fGitCommit;            //Git commit SHA1 of built
   TString         fGitBranch;            //Git branch
   TString         fGitDate;              //Date and time when make was run
   Int_t           fTimer;                //Timer flag
   std::atomic<TApplication*> fApplication;  //Pointer to current application
   TInterpreter    *fInterpreter;         //Command interpreter
   Bool_t          fBatch;                //True if session without graphics
   TString         fWebDisplay;           //If not empty it defines where web graphics should be rendered (cef, qt5, browser...)
   Bool_t          fIsWebDisplay;         //True if session with graphics on web
   Bool_t          fIsWebDisplayBatch;    //True if session with graphics on web and batch mode
   Bool_t          fEditHistograms;       //True if histograms can be edited with the mouse
   Bool_t          fFromPopUp;            //True if command executed from a popup menu
   Bool_t          fMustClean;            //True if object destructor scans canvases
   Bool_t          fForceStyle;           //Force setting of current style when reading objects
   Bool_t          fInterrupt;            //True if macro should be interrupted
   Bool_t          fEscape;               //True if ESC has been pressed
   Bool_t          fExecutingMacro;       //True while executing a TMacro
   Int_t           fEditorMode;           //Current Editor mode
   const TObject   *fPrimitive;           //Currently selected primitive
   TVirtualPad     *fSelectPad;           //Currently selected pad
   TCollection     *fClasses;             //List of classes definition
   TCollection     *fTypes;               //List of data types definition
   TListOfFunctionTemplates *fFuncTemplate; //List of global function templates
   TListOfDataMembers*fGlobals;             //List of global variables
   TListOfFunctions*fGlobalFunctions;     //List of global functions
   TSeqCollection  *fClosedObjects;       //List of closed objects from the list of files and sockets, so we can delete them if neededCl.
   TSeqCollection  *fFiles;               //List of files
   TSeqCollection  *fMappedFiles;         //List of memory mapped files
   TSeqCollection  *fSockets;             //List of network sockets
   TSeqCollection  *fCanvases;            //List of canvases
   TSeqCollection  *fStyles;              //List of styles
   TCollection     *fFunctions;           //List of analytic functions
   TSeqCollection  *fTasks;               //List of tasks
   TSeqCollection  *fColors;              //List of colors
   TSeqCollection  *fGeometries;          //List of geometries
   TSeqCollection  *fBrowsers;            //List of browsers
   TSeqCollection  *fSpecials;            //List of special objects
   TSeqCollection  *fCleanups;            //List of recursiveRemove collections
   TSeqCollection  *fMessageHandlers;     //List of message handlers
   TSeqCollection  *fStreamerInfo;        //List of active StreamerInfo classes
   TCollection     *fClassGenerators;     //List of user defined class generators;
   TSeqCollection  *fSecContexts;         //List of security contexts (TSecContext)
   TSeqCollection  *fProofs;              //List of proof sessions
   TSeqCollection  *fClipboard;           //List of clipboard objects
   TSeqCollection  *fDataSets;            //List of data sets (TDSet or TChain)
   AListOfEnums_t   fEnums;               //List of enum types
   TProcessUUID    *fUUIDs;               //Pointer to TProcessID managing TUUIDs
   TFolder         *fRootFolder;          //top level folder //root
   TList           *fBrowsables;          //List of browsables
   TPluginManager  *fPluginManager;       //Keeps track of plugin library handlers
   TString         fCutClassName;         //Name of default CutG class in graphics editor
   TString         fDefCanvasName;        //Name of default canvas

                  TROOT();                //Only used by Dictionary
   void           InitSystem();           //Operating System interface
   void           InitThreads();          //Initialize threads library
   void           InitInterpreter();      //Initialize interpreter (cling)
   void           ReadGitInfo();          //Read Git commit SHA1 and branch name
   void          *operator new(size_t l) { return TObject::operator new(l); }
   void          *operator new(size_t l, void *ptr) { return TObject::operator new(l,ptr); }

   friend class ::ROOT::Internal::TROOTAllocator;

   TListOfFunctions*GetGlobalFunctions();

public:

   typedef std::vector<std::pair<std::string, int> > FwdDeclArgsToKeepCollection_t;

                     TROOT(const char *name, const char *title, VoidFuncPtr_t *initfunc = 0);
   virtual           ~TROOT();
   void              AddClass(TClass *cl);
   void              AddClassGenerator(TClassGenerator *gen);
   virtual void      Append(TObject *obj, Bool_t replace = kFALSE);
   void              Browse(TBrowser *b);
   Bool_t            ClassSaved(TClass *cl);
   void              CloseFiles();
   void              EndOfProcessCleanups();
   virtual TObject  *FindObject(const char *name) const;
   virtual TObject  *FindObject(const TObject *obj) const;
   virtual TObject  *FindObjectAny(const char *name) const;
   virtual TObject  *FindObjectAnyFile(const char *name) const;
   TObject          *FindSpecialObject(const char *name, void *&where);
   const char       *FindObjectClassName(const char *name) const;
   const char       *FindObjectPathName(const TObject *obj) const;
   TClass           *FindSTLClass(const char *name, Bool_t load, Bool_t silent = kFALSE) const;
   void              ForceStyle(Bool_t force = kTRUE) { fForceStyle = force; }
   Bool_t            FromPopUp() const { return fFromPopUp; }
   TPluginManager   *GetPluginManager() const { return fPluginManager; }
   TApplication     *GetApplication() const { return fApplication; }
   TInterpreter     *GetInterpreter() const { return fInterpreter; }
   TClass           *GetClass(const char *name, Bool_t load = kTRUE, Bool_t silent = kFALSE) const;
   TClass           *GetClass(const std::type_info &typeinfo, Bool_t load = kTRUE, Bool_t silent = kFALSE) const;
   TColor           *GetColor(Int_t color) const;
   const char       *GetConfigOptions() const { return fConfigOptions; }
   const char       *GetConfigFeatures() const { return fConfigFeatures; }
   const char       *GetCutClassName() const { return fCutClassName; }
   const char       *GetDefCanvasName() const { return fDefCanvasName; }
   Bool_t            GetEditHistograms() const { return fEditHistograms; }
   Int_t             GetEditorMode() const { return fEditorMode; }
   Bool_t            GetForceStyle() const { return fForceStyle; }
   Int_t             GetBuiltDate() const { return fBuiltDate; }
   Int_t             GetBuiltTime() const { return fBuiltTime; }
   const char       *GetGitCommit() const { return fGitCommit; }
   const char       *GetGitBranch() const { return fGitBranch; }
   const char       *GetGitDate();
   Int_t             GetVersionDate() const { return fVersionDate; }
   Int_t             GetVersionTime() const { return fVersionTime; }
   Int_t             GetVersionInt() const { return fVersionInt; }
   Int_t             GetVersionCode() const { return fVersionCode; }
   const char       *GetVersion() const { return fVersion; }
   TCollection      *GetListOfClasses() const { return fClasses; }
   TSeqCollection   *GetListOfColors() const { return fColors; }
   TCollection      *GetListOfTypes(Bool_t load = kFALSE);
   TCollection      *GetListOfGlobals(Bool_t load = kFALSE);
   TCollection      *GetListOfGlobalFunctions(Bool_t load = kFALSE);
   TSeqCollection   *GetListOfClosedObjects() const { return fClosedObjects; }
   TSeqCollection   *GetListOfFiles() const       { return fFiles; }
   TSeqCollection   *GetListOfMappedFiles() const { return fMappedFiles; }
   TSeqCollection   *GetListOfSockets() const     { return fSockets; }
   TSeqCollection   *GetListOfCanvases() const    { return fCanvases; }
   TSeqCollection   *GetListOfStyles() const      { return fStyles; }
   TCollection      *GetListOfFunctions() const   { return fFunctions; }
   TCollection      *GetListOfFunctionOverloads(const char* name) const;
   TSeqCollection   *GetListOfGeometries() const  { return fGeometries; }
   TSeqCollection   *GetListOfBrowsers() const    { return fBrowsers; }
   TSeqCollection   *GetListOfSpecials() const    { return fSpecials; }
   TSeqCollection   *GetListOfTasks() const       { return fTasks; }
   TSeqCollection   *GetListOfCleanups() const    { return fCleanups; }
   TSeqCollection   *GetListOfStreamerInfo() const { return fStreamerInfo; }
   TSeqCollection   *GetListOfMessageHandlers() const { return fMessageHandlers; }
   TCollection      *GetListOfClassGenerators() const { return fClassGenerators; }
   TSeqCollection   *GetListOfSecContexts() const { return fSecContexts; }
   TSeqCollection   *GetListOfProofs() const { return fProofs; }
   TSeqCollection   *GetClipboard() const { return fClipboard; }
   TSeqCollection   *GetListOfDataSets() const { return fDataSets; }
   TCollection      *GetListOfEnums(Bool_t load = kFALSE);
   TCollection      *GetListOfFunctionTemplates();
   TList            *GetListOfBrowsables() const { return fBrowsables; }
   TDataType        *GetType(const char *name, Bool_t load = kFALSE) const;
   TFile            *GetFile() const { if (gDirectory != this) return gDirectory->GetFile(); else return 0;}
   TFile            *GetFile(const char *name) const;
   TFunctionTemplate*GetFunctionTemplate(const char *name);
   TStyle           *GetStyle(const char *name) const;
   TObject          *GetFunction(const char *name) const;
   TGlobal          *GetGlobal(const char *name, Bool_t load = kFALSE) const;
   TGlobal          *GetGlobal(const TObject *obj, Bool_t load = kFALSE) const;
   TFunction        *GetGlobalFunction(const char *name, const char *params = 0, Bool_t load = kFALSE);
   TFunction        *GetGlobalFunctionWithPrototype(const char *name, const char *proto = 0, Bool_t load = kFALSE);
   TObject          *GetGeometry(const char *name) const;
   const TObject    *GetSelectedPrimitive() const { return fPrimitive; }
   TVirtualPad      *GetSelectedPad() const { return fSelectPad; }
   Int_t             GetNclasses() const;
   Int_t             GetNtypes() const;
   TFolder          *GetRootFolder() const { return fRootFolder; }
   TProcessUUID     *GetUUIDs() const { return fUUIDs; }
   const TString    &GetWebDisplay() const { return fWebDisplay; }
   void              Idle(UInt_t idleTimeInSec, const char *command = 0);
   Int_t             IgnoreInclude(const char *fname, const char *expandedfname);
   Bool_t            IsBatch() const { return fBatch; }
   Bool_t            IsExecutingMacro() const { return fExecutingMacro; }
   Bool_t            IsFolder() const { return kTRUE; }
   Bool_t            IsInterrupted() const { return fInterrupt; }
   Bool_t            IsEscaped() const { return fEscape; }
   Bool_t            IsLineProcessing() const { return fLineIsProcessing ? kTRUE : kFALSE; }
   Bool_t            IsProofServ() const { return fName == "proofserv" ? kTRUE : kFALSE; }
   Bool_t            IsRootFile(const char *filename) const;
   Bool_t            IsWebDisplay() const { return fIsWebDisplay; }
   Bool_t            IsWebDisplayBatch() const { return fIsWebDisplayBatch; }
   void              ls(Option_t *option = "") const;
   Int_t             LoadClass(const char *classname, const char *libname, Bool_t check = kFALSE);
   TClass           *LoadClass(const char *name, Bool_t silent = kFALSE) const;
   Int_t             LoadMacro(const char *filename, Int_t *error = 0, Bool_t check = kFALSE);
   Long_t            Macro(const char *filename, Int_t *error = 0, Bool_t padUpdate = kTRUE);
   TCanvas          *MakeDefCanvas() const;
   void              Message(Int_t id, const TObject *obj);
   Bool_t            MustClean() const { return fMustClean; }
   Long_t            ProcessLine(const char *line, Int_t *error = 0);
   Long_t            ProcessLineSync(const char *line, Int_t *error = 0);
   Long_t            ProcessLineFast(const char *line, Int_t *error = 0);
   Bool_t            ReadingObject() const;
   void              RecursiveRemove(TObject *obj);
   void              RefreshBrowsers();
   static void       RegisterModule(const char* modulename,
                                    const char** headers,
                                    const char** includePaths,
                                    const char* payLoadCode,
                                    const char* fwdDeclCode,
                                    void (*triggerFunc)(),
                                    const FwdDeclArgsToKeepCollection_t& fwdDeclsArgToSkip,
                                    const char** classesHeaders,
                                    bool hasCxxModule = false);
   TObject          *Remove(TObject*);
   void              RemoveClass(TClass *);
   void              Reset(Option_t *option="");
   void              SaveContext();
   void              SetApplication(TApplication *app) { fApplication = app; }
   void              SetBatch(Bool_t batch = kTRUE) { fBatch = batch; }
   void              SetWebDisplay(const char *webdisplay);
   void              SetCutClassName(const char *name = "TCutG");
   void              SetDefCanvasName(const char *name = "c1") { fDefCanvasName = name; }
   void              SetEditHistograms(Bool_t flag = kTRUE) { fEditHistograms = flag; }
   void              SetEditorMode(const char *mode = "");
   void              SetExecutingMacro(Bool_t flag = kTRUE) { fExecutingMacro = flag; }
   void              SetFromPopUp(Bool_t flag = kTRUE) { fFromPopUp = flag; }
   void              SetInterrupt(Bool_t flag = kTRUE) { fInterrupt = flag; }
   void              SetEscape(Bool_t flag = kTRUE) { fEscape = flag; }
   void              SetLineIsProcessing() { fLineIsProcessing++; }
   void              SetLineHasBeenProcessed() { if (fLineIsProcessing) fLineIsProcessing--; }
   void              SetReadingObject(Bool_t flag = kTRUE);
   void              SetMustClean(Bool_t flag = kTRUE) { fMustClean=flag; }
   void              SetSelectedPrimitive(const TObject *obj) { fPrimitive = obj; }
   void              SetSelectedPad(TVirtualPad *pad) { fSelectPad = pad; }
   void              SetStyle(const char *stylename = "Default");
   void              Time(Int_t casetime=1) { fTimer = casetime; }
   Int_t             Timer() const { return fTimer; }

   //---- static functions
   static Int_t       DecreaseDirLevel();
   static Int_t       GetDirLevel();
   static const char *GetMacroPath();
   static void        SetMacroPath(const char *newpath);
   static Int_t       IncreaseDirLevel();
   static void        IndentLevel();
   static void        Initialize();
   static Bool_t      Initialized();
   static Bool_t      MemCheck();
   static void        SetDirLevel(Int_t level = 0);
   static Int_t       ConvertVersionCode2Int(Int_t code);
   static Int_t       ConvertVersionInt2Code(Int_t v);
   static Int_t       RootVersionCode();
   static const std::vector<std::string> &AddExtraInterpreterArgs(const std::vector<std::string> &args);
   static const char**&GetExtraInterpreterArgs();

   static const TString& GetRootSys();
   static const TString& GetBinDir();
   static const TString& GetLibDir();
   static const TString& GetIncludeDir();
   static const TString& GetEtcDir();
   static const TString& GetDataDir();
   static const TString& GetDocDir();
   static const TString& GetMacroDir();
   static const TString& GetTutorialDir();
   static const TString& GetSourceDir();
   static const TString& GetIconPath();
   static const TString& GetTTFFontDir();

   // Backward compatibility function - do not use for new code
   static const char *GetTutorialsDir();
   static void ShutDown();

   ClassDef(TROOT,0)  //Top level (or root) structure for all classes
};


namespace ROOT {
   TROOT *GetROOT();
   namespace Internal {
      R__EXTERN TROOT *gROOTLocal;

      inline void SetRequireCleanup(TObject &obj) {
         obj.SetBit(kIsReferenced);
         obj.SetUniqueID(0);
      }

      inline Bool_t RequiresCleanup(TObject &obj) {
         return obj.TestBit(kIsReferenced) && obj.GetUniqueID() == 0;
      }
   }

   /// \brief call RecursiveRemove for obj if gROOT is valid
   /// and obj.TestBit(kMustCleanup) is true.
   /// Note: this reset the kMustCleanup bit to allow
   /// harmless multiple call to this function.
   inline void CallRecursiveRemoveIfNeeded(TObject &obj)
   {
      if (obj.TestBit(kMustCleanup)) {
         TROOT *root = ROOT::Internal::gROOTLocal;
         if (root && root != &obj && (root->MustClean() || Internal::RequiresCleanup(obj))) {
            root->RecursiveRemove(&obj);
            obj.ResetBit(kMustCleanup);
         }
      }
   }
}
#define gROOT (ROOT::GetROOT())

#endif
