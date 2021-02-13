/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooWorkspace.h,v 1.3 2007/07/16 21:04:28 wouter Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_WORKSPACE
#define ROO_WORKSPACE

#include "RooPrintable.h"
#include "RooArgSet.h"
#include "RooLinkedList.h"
#include "RooCmdArg.h"
#include "RooExpensiveObjectCache.h" 
#include "TUUID.h"
#include <map>
#include <list>
#include <string>
#include "ROOT/RMakeUnique.hxx"

class TClass ;
class RooAbsPdf ;
class RooAbsData ;
class RooRealVar ;
class RooCategory ;
class RooAbsReal ;
class RooAbsCategory ;
class RooFactoryWSTool ;
class RooAbsStudy ;

#include "TNamed.h"
#include "TDirectoryFile.h"

class RooWorkspace : public TNamed {
public:

  RooWorkspace() ;
  RooWorkspace(const char* name, Bool_t doCINTExport) ;
  RooWorkspace(const char* name, const char* title=0) ;
  RooWorkspace(const RooWorkspace& other) ;
  ~RooWorkspace() ;

  void exportToCint(const char* namespaceName=0) ;

  Bool_t importClassCode(const char* pat="*", Bool_t doReplace=kFALSE) ;
  Bool_t importClassCode(TClass* theClass, Bool_t doReplace=kFALSE) ;

  // Import functions for dataset, functions, generic objects
  Bool_t import(const RooAbsArg& arg, 
		const RooCmdArg& arg1=RooCmdArg(),const RooCmdArg& arg2=RooCmdArg(),const RooCmdArg& arg3=RooCmdArg(),
		const RooCmdArg& arg4=RooCmdArg(),const RooCmdArg& arg5=RooCmdArg(),const RooCmdArg& arg6=RooCmdArg(),
		const RooCmdArg& arg7=RooCmdArg(),const RooCmdArg& arg8=RooCmdArg(),const RooCmdArg& arg9=RooCmdArg()) ;
  Bool_t import(const RooArgSet& args, 
		const RooCmdArg& arg1=RooCmdArg(),const RooCmdArg& arg2=RooCmdArg(),const RooCmdArg& arg3=RooCmdArg(),
		const RooCmdArg& arg4=RooCmdArg(),const RooCmdArg& arg5=RooCmdArg(),const RooCmdArg& arg6=RooCmdArg(),
		const RooCmdArg& arg7=RooCmdArg(),const RooCmdArg& arg8=RooCmdArg(),const RooCmdArg& arg9=RooCmdArg()) ;
  Bool_t import(RooAbsData& data, 
		const RooCmdArg& arg1=RooCmdArg(),const RooCmdArg& arg2=RooCmdArg(),const RooCmdArg& arg3=RooCmdArg(),
		const RooCmdArg& arg4=RooCmdArg(),const RooCmdArg& arg5=RooCmdArg(),const RooCmdArg& arg6=RooCmdArg(),
		const RooCmdArg& arg7=RooCmdArg(),const RooCmdArg& arg8=RooCmdArg(),const RooCmdArg& arg9=RooCmdArg()) ;
  Bool_t import(const char *fileSpec, 
		const RooCmdArg& arg1=RooCmdArg(),const RooCmdArg& arg2=RooCmdArg(),const RooCmdArg& arg3=RooCmdArg(),
		const RooCmdArg& arg4=RooCmdArg(),const RooCmdArg& arg5=RooCmdArg(),const RooCmdArg& arg6=RooCmdArg(),
		const RooCmdArg& arg7=RooCmdArg(),const RooCmdArg& arg8=RooCmdArg(),const RooCmdArg& arg9=RooCmdArg()) ;
  Bool_t import(TObject& object, Bool_t replaceExisting=kFALSE) ;
  Bool_t import(TObject& object, const char* aliasName, Bool_t replaceExisting=kFALSE) ;

  // Transaction management interface for multi-step import operations
  Bool_t startTransaction() ;
  Bool_t cancelTransaction() ;
  Bool_t commitTransaction() ;

  // Named set management
  Bool_t defineSet(const char* name, const RooArgSet& aset, Bool_t importMissing=kFALSE) ;
  Bool_t defineSet(const char* name, const char* contentList) ;
  Bool_t extendSet(const char* name, const char* newContents) ;
  Bool_t renameSet(const char* name, const char* newName) ;
  Bool_t removeSet(const char* name) ;
  const RooArgSet* set(const char* name) ;

  // Import, load and save parameter value snapshots
  Bool_t saveSnapshot(const char* name, const char* paramNames) ;
  Bool_t saveSnapshot(const char* name, const RooArgSet& params, Bool_t importValues=kFALSE) ;
  Bool_t loadSnapshot(const char* name) ;  
  const RooArgSet* getSnapshot(const char* name) const ;

  void merge(const RooWorkspace& /*other*/) {} ;

  // Join p.d.f.s and datasets for simultaneous analysis
  //   RooAbsPdf* joinPdf(const char* jointPdfName, const char* indexName, const char* inputMapping) ; 
  //   RooAbsData* joinData(const char* jointDataName, const char* indexName, const char* inputMapping) ; 

  // Accessor functions 
  RooAbsPdf* pdf(const char* name) const ;
  RooAbsReal* function(const char* name) const ;
  RooRealVar* var(const char* name) const ;
  RooCategory* cat(const char* name) const ;
  RooAbsCategory* catfunc(const char* name) const ;
  RooAbsData* data(const char* name) const ;
  RooAbsData* embeddedData(const char* name) const ;
  RooAbsArg* arg(const char* name) const ;
  RooAbsArg* fundArg(const char* name) const ;
  RooArgSet argSet(const char* nameList) const ;
  TIterator* componentIterator() const { return _allOwnedNodes.createIterator() ; }
  const RooArgSet& components() const { return _allOwnedNodes ; }
  TObject* genobj(const char* name) const ;
  TObject* obj(const char* name) const ;

  // Group accessors
  RooArgSet allVars() const;
  RooArgSet allCats() const ;
  RooArgSet allFunctions() const ;
  RooArgSet allCatFunctions() const ;
  RooArgSet allPdfs() const ;
  RooArgSet allResolutionModels() const ;
  std::list<RooAbsData*> allData() const ;
  std::list<RooAbsData*> allEmbeddedData() const ;
  std::list<TObject*> allGenericObjects() const ;

  Bool_t makeDir() ; 
  Bool_t cd(const char* path = 0) ;

  Bool_t writeToFile(const char* fileName, Bool_t recreate=kTRUE) ;

  /// Make internal collection use an unordered_map for
  /// faster searching. Important when large trees are
  /// imported / or modified in the workspace.
  /// Note that RooAbsCollection may eventually switch
  /// this on by itself.
  void useFindsWithHashLookup(bool flag) {
    _allOwnedNodes.useHashMapForFind(flag);
  }

  virtual void RecursiveRemove(TObject *obj);

  // Tools management
  RooFactoryWSTool& factory() ;
  RooAbsArg* factory(const char* expr) ;

  // RooStudyManager modules
  Bool_t addStudy(RooAbsStudy& study) ;  
  TIterator* studyIterator() { return _studyMods.MakeIterator() ; }
  void clearStudies() ;

  // Print function
  void Print(Option_t* opts=0) const ;

  static void autoImportClassCode(Bool_t flag) ;
 
  static void addClassDeclImportDir(const char* dir) ;
  static void addClassImplImportDir(const char* dir) ;
  static void setClassFileExportDir(const char* dir=0) ; 

  const TUUID& uuid() const { return _uuid ; }

  RooExpensiveObjectCache& expensiveObjectCache() { return _eocache ; }

  class CodeRepo : public TObject {
  public:
    CodeRepo(RooWorkspace* wspace=0) : _wspace(wspace), _compiledOK(kTRUE) {} ;

    CodeRepo(const CodeRepo& other, RooWorkspace* wspace=0) : TObject(other) ,
          _wspace(wspace?wspace:other._wspace), 
          _c2fmap(other._c2fmap),
          _fmap(other._fmap),
          _ehmap(other._ehmap),
          _compiledOK(other._compiledOK) {} ;

    virtual ~CodeRepo() {} ;

    Bool_t autoImportClass(TClass* tc, Bool_t doReplace=kFALSE) ;
    Bool_t compileClasses() ;

    Bool_t compiledOK() const { return _compiledOK ; }

    std::string listOfClassNames() const ;

       

    class ClassRelInfo {
    public:
      TString _baseName;      
      TString _fileBase ;
    } ;

    class ClassFiles {
    public:
      ClassFiles() : _extracted(kFALSE) {}
      TString _hext ;
      TString _hfile ;
      TString _cxxfile ;
      Bool_t _extracted ; 
    } ;    
    

    class ExtraHeader {     
    public:
      TString _hname ;
      TString _hfile ;
   } ;

  protected:
    RooWorkspace* _wspace ; // owning workspace
    std::map<TString,ClassRelInfo> _c2fmap ; // List of contained classes
    std::map<TString,ClassFiles> _fmap ; // List of contained files
    std::map<TString,ExtraHeader> _ehmap ; // List of extra header files
    Bool_t _compiledOK ; //! Flag indicating that classes compiled OK

    ClassDef(CodeRepo,2) ; // Code repository for RooWorkspace
  } ;


  class WSDir : public TDirectoryFile {    
  public:
    WSDir(const char* name, const char* title, RooWorkspace* wspace) : 
      TDirectoryFile(name,title,"RooWorkspace::WSDir",0), 
      _wspace(wspace) 
      {
      }

    virtual ~WSDir() { Clear("nodelete") ; } ; 


    virtual void Add(TObject*,Bool_t) ; 
    virtual void Append(TObject*,Bool_t) ; 

  protected:
    friend class RooWorkspace ;
    void InternalAppend(TObject* obj) ;
    RooWorkspace* _wspace ; //! do not persist

    ClassDef(WSDir,1) ; // TDirectory representation of RooWorkspace
  } ;


 private:
    friend class RooAbsArg;
    friend class RooAbsPdf;
    Bool_t defineSetInternal(const char *name, const RooArgSet &aset);

    Bool_t isValidCPPID(const char *name);
    void exportObj(TObject *obj);
    void unExport();

    friend class CodeRepo;
    static std::list<std::string> _classDeclDirList;
    static std::list<std::string> _classImplDirList;
    static std::string _classFileExportDir;

    TUUID _uuid; // Unique workspace ID

    static Bool_t _autoClass; // Automatic import of non-distribution class code

    CodeRepo _classes; // Repository of embedded class code. This data member _must_ be first

    RooArgSet _allOwnedNodes;                    // List of owned pdfs and components
    RooLinkedList _dataList;                     // List of owned datasets
    RooLinkedList _embeddedDataList;             // List of owned datasets that are embedded in pdfs
    RooLinkedList _views;                        // List of model views
    RooLinkedList _snapshots;                    // List of parameter snapshots
    RooLinkedList _genObjects;                   // List of generic objects
    RooLinkedList _studyMods;                    // List if StudyManager modules
    std::map<std::string, RooArgSet> _namedSets; // Map of named RooArgSets

    WSDir *_dir; //! Transient ROOT directory representation of workspace

    RooExpensiveObjectCache _eocache; // Cache for expensive objects

    std::unique_ptr<RooFactoryWSTool> _factory; //! Factory tool associated with workspace

    Bool_t _doExport;          //! Export contents of workspace to CINT?
    std::string _exportNSName; //! Name of CINT namespace to which contents are exported

    Bool_t _openTrans;       //! Is there a transaction open?
    RooArgSet _sandboxNodes; //! Sandbox for incoming objects in a transaction

    ClassDef(RooWorkspace, 8) // Persistable project container for (composite) pdfs, functions, variables and datasets
} ;

#endif
